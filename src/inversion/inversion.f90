module inversion_mod
    use inversion_utils
    use parameters, only : nx, ny, nz, dxi
    use physics, only : f_cor
    use constants, only : zero, two, f12
    use timer, only : start_timer, stop_timer
    implicit none

    integer :: vor2vel_timer,   &
               vtend_timer

    contains

        ! Given the vorticity vector field (svortg) in spectral space, this
        ! returns the associated velocity field (velog) and the velocity
        ! gradient tensor (velgradg).  Note: the
        ! vorticity is modified to be solenoidal and spectrally filtered.
        subroutine vor2vel(svortg,  velog,  velgradg)
            double precision, intent(in)    :: svortg(0:nz, 0:nx-1, 0:ny-1, 3)
            double precision, intent(out)   :: velog(-1:nz+1, 0:ny-1, 0:nx-1, 3)
            double precision, intent(out)   :: velgradg(-1:nz+1, 0:ny-1, 0:nx-1, 5)
            double precision                :: svelog(0:nz, 0:nx-1, 0:ny-1, 3)
            double precision                :: as(0:nz, 0:nx-1, 0:ny-1)
            double precision                :: bs(0:nz, 0:nx-1, 0:ny-1)
            double precision                :: cs(0:nz, 0:nx-1, 0:ny-1)
            double precision                :: ds(0:nz, 0:nx-1, 0:ny-1)
            double precision                :: es(0:nz, 0:nx-1, 0:ny-1)
            double precision                :: ubar(0:nz), vbar(0:nz)
            double precision                :: uavg, vavg
            integer                         :: iz

            call start_timer(vor2vel_timer)

            !-------------------------------------------------------------
            ! Apply 2D Hou and Li filter:
            !$omp parallel shared(as, bs, cs, filt, svortg, nz) private(iz) default(none)
            !$omp do
            do iz = 0, nz
                as(iz, :, :) = filt * svortg(iz, :, :, 1)
                bs(iz, :, :) = filt * svortg(iz, :, :, 2)
                cs(iz, :, :) = filt * svortg(iz, :, :, 3)
            enddo
            !$omp end do
            !$omp end parallel


            !Define horizontally-averaged flow by integrating horizontal vorticity:
            ubar(0) = zero
            vbar(0) = zero
            do iz = 0, nz-1
                ubar(iz+1) = ubar(iz) + dz2 * (bs(iz, 0, 0) + bs(iz+1, 0, 0))
                vbar(iz+1) = vbar(iz) - dz2 * (as(iz, 0, 0) + as(iz+1, 0, 0))
            enddo

            ! remove the mean value to have zero net momentum
            uavg = sum(ubar(1:nz-1) + f12 * ubar(nz)) / dble(nz)
            vavg = sum(vbar(1:nz-1) + f12 * vbar(nz)) / dble(nz)
            do iz = 0, nz
                ubar(iz) = ubar(iz) - uavg
                vbar(iz) = vbar(iz) - vavg
            enddo

            !Form source term for inversion of vertical velocity:
            call diffy(as, ds)
            call diffx(bs, es)

            !$omp parallel
            !$omp workshare
            ds = ds - es
            !$omp end workshare
            !$omp end parallel

            !as & bs are now free to re-use

            !Invert to find vertical velocity \hat{w} (store in ds, spectrally):
            call lapinv0(ds)

            !Find \hat{w}' (store in es, spectrally):
            call diffz(ds, es)

            !Find x velocity component \hat{u}:
            call diffx(es, as)
            call diffy(cs, bs)

            !$omp parallel do
            do iz = 0, nz
                as(iz, :, :) = k2l2i * (as(iz, :, :) + bs(iz, :, :))
            enddo
            !$omp end parallel do

            !Add horizontally-averaged flow:
            as(:, 0, 0) = ubar

            svelog(:, :, :, 1) = as

            !Get u in physical space:
            call fftxys2p(as, velog(0:nz, :, :, 1))

            !Find y velocity component \hat{v}:
            call diffy(es, as)
            call diffx(cs, bs)

            !$omp parallel do
            do iz = 0, nz
                as(iz, :, :) = k2l2i * (as(iz, :, :) - bs(iz, :, :))
            enddo
            !$omp end parallel do

            !Add horizontally-averaged flow:
            as(:, 0, 0) = vbar

            svelog(:, :, :, 2) = as

            !Get v in physical space:
            call fftxys2p(as, velog(0:nz, :, :, 2))

            svelog(:, :, :, 3) = ds

            !Get w in physical space:
            call fftxys2p(ds, velog(0:nz, :, :, 3))

            ! compute the velocity gradient tensor
            call vel2vgrad(svelog, velgradg)

            call stop_timer(vor2vel_timer)

        end subroutine


        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Compute the gridded velocity gradient tensor
        subroutine vel2vgrad(svelog, velgradg)
            double precision, intent(in)  :: svelog(0:nz, 0:nx-1, 0:ny-1, 3)
            double precision, intent(out) :: velgradg(-1:nz+1, 0:ny-1, 0:nx-1, 5)
            double precision              :: ds(0:nz, 0:nx-1, 0:ny-1) ! spectral derivatives

            ! x component:
            call diffx(svelog(:, :, :, 1), ds)         ! u_x = du/dx in spectral space
            call fftxys2p(ds, velgradg(0:nz, :, :, 1)) ! u_x in physical space

            call diffy(svelog(:, :, :, 1), ds)         ! u_y = du/dy in spectral space
            call fftxys2p(ds, velgradg(0:nz, :, :, 2)) ! u_y in physical space

            call diffx(svelog(:, :, :, 3), ds)         ! w_x = dw/dx in spectral space
            call fftxys2p(ds, velgradg(0:nz, :, :, 4)) ! w_x in physical space

            ! use symmetry in du/dx and du/dy to fill z grid points outside domain:
            velgradg(  -1, :, :, 1) =  velgradg(   1, :, :, 1) ! lower boundary du/dx
            velgradg(nz+1, :, :, 1) =  velgradg(nz-1, :, :, 1) ! upper boundary du/dx
            velgradg(  -1, :, :, 2) =  velgradg(   1, :, :, 2) ! lower boundary du/dy
            velgradg(nz+1, :, :, 2) =  velgradg(nz-1, :, :, 2) ! upper boundary du/dy

            ! use anti-symmetry for dw/dx to fill z grid points outside domain:
            velgradg(  -1, :, :, 4) = -velgradg(   1, :, :, 4) ! lower boundary dw/dx
            velgradg(nz+1, :, :, 4) = -velgradg(nz-1, :, :, 4) ! upper boundary dw/dx

            ! y & z components:
            call diffy(svelog(:, :, :, 2), ds)         ! v_y = dv/dy in spectral space
            call fftxys2p(ds, velgradg(0:nz, :, :, 3)) ! v_y in physical space

            call diffy(svelog(:, :, :, 3), ds)         ! w_y = dw/dy in spectral space
            call fftxys2p(ds, velgradg(0:nz, :, :, 5)) ! w_y in physical space

            ! use symmetry in dv/dy to fill z grid points outside domain:
            velgradg(  -1, :, :, 3) = velgradg(   1, :, :, 3) ! lower boundary dv/dy
            velgradg(nz+1, :, :, 3) = velgradg(nz-1, :, :, 3) ! upper boundary dv/dy

            ! use anti-symmetry in dw/dy to fill z grid points outside domain:
            ! w_y(-1) = -w_y(1) and w_y(nz+1) = -w_y(nz-1)
            velgradg(  -1, :, :, 5) = -velgradg(   1, :, :, 5) ! lower boundary dw/dy
            velgradg(nz+1, :, :, 5) = -velgradg(nz-1, :, :, 5) ! upper boundary dw/dy

        end subroutine vel2vgrad

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Compute the gridded vorticity tendency:
        subroutine vorticity_tendency(vortg, velog, buoyg, vtend)
            double precision, intent(in)  :: vortg(0:nz, 0:ny-1, 0:nx-1, 3)
            double precision, intent(in)  :: velog(0:nz, 0:ny-1, 0:nx-1, 3)
            double precision, intent(in)  :: buoyg(0:nz, 0:ny-1, 0:nx-1)
            double precision, intent(out) :: vtend(0:nz, 0:ny-1, 0:nx-1, 3)
            double precision              :: f(0:nz, 0:ny-1, 0:nx-1, 3)

            call start_timer(vtend_timer)

            ! Eqs. 10 and 11 of MPIC paper
            f(:, : , :, 1) = (vortg(:, :, :, 1) + f_cor(1)) * velog(:, :, :, 1)
            f(:, : , :, 2) = (vortg(:, :, :, 2) + f_cor(2)) * velog(:, :, :, 1) + buoyg
            f(:, : , :, 3) = (vortg(:, :, :, 3) + f_cor(3)) * velog(:, :, :, 1)

            call divergence(f, vtend(0:nz, :, :, 1))

            f(:, : , :, 1) = (vortg(:, :, :, 1) + f_cor(1)) * velog(:, :, :, 2) - buoyg
            f(:, : , :, 2) = (vortg(:, :, :, 2) + f_cor(2)) * velog(:, :, :, 2)
            f(:, : , :, 3) = (vortg(:, :, :, 3) + f_cor(3)) * velog(:, :, :, 2)

            call divergence(f, vtend(0:nz, :, :, 2))

            f(:, : , :, 1) = (vortg(:, :, :, 1) + f_cor(1)) * velog(:, :, :, 3)
            f(:, : , :, 2) = (vortg(:, :, :, 2) + f_cor(2)) * velog(:, :, :, 3)
            f(:, : , :, 3) = (vortg(:, :, :, 3) + f_cor(3)) * velog(:, :, :, 3)

            call divergence(f, vtend(0:nz, :, :, 3))


!             ! Fill boundary values:
!             ! \omegax * du/dx + \omegay * dv/dx (where dv/dx = \omegaz + du/dy)
!             vtend(0, :, :, 1) = vortg(0, :, :, 1) * velgradg(0, :, :, 1) &
!                               + vortg(0, :, :, 2) * (vortg(0, :, :, 3 ) + velgradg(0, :, :, 2))
!
!             vtend(nz, :, :, 1) = vortg(nz, :, :, 1) * velgradg(nz, :, :, 1) &
!                                + vortg(nz, :, :, 2) * (vortg(nz, :, :, 3 ) + velgradg(nz, :, :, 2))
!
!             ! \omegax * du/dy + \omegay * dv/dy
!             vtend(0, :, :, 2) = vortg(0, :, :, 1) * velgradg(0, :, :, 2) &
!                               + vortg(0, :, :, 2) * velgradg(0, :, :, 3)
!
!             vtend(nz, :, :, 2) = vortg(nz, :, :, 1) * velgradg(nz, :, :, 2) &
!                                + vortg(nz, :, :, 2) * velgradg(nz, :, :, 3)
!
!             ! - \omegaz * (du/dx + dv/dy)
!             vtend(0,  :, :, 3) = - vortg(0,  :, :, 3) * (velgradg(0,  :, :, 1) + velgradg(0,  :, :, 3))
!             vtend(nz, :, :, 3) = - vortg(nz, :, :, 3) * (velgradg(nz, :, :, 1) + velgradg(nz, :, :, 3))

            call stop_timer(vtend_timer)

        end subroutine vorticity_tendency

        subroutine divergence(f, div)
            double precision, intent(in)  :: f(-1:nz+1, 0:ny-1, 0:nx-1, 3)
            double precision, intent(out) :: div(0:nz, 0:ny-1, 0:nx-1)
            double precision              :: df(0:nz, 0:ny-1, 0:nx-1)
            integer                       :: i

            ! calculate df/dx with central differencing
            do i = 1, nx-2
                div(0:nz, :, i) = f12 * dxi(1) * (f(0:nz, :, i+1, 1) - f(0:nz, :, i-1, 1))
            enddo
            div(0:nz, :, 0)    = f12 * dxi(1) * (f(0:nz, :, 1, 1) - f(0:nz, :, nx-1, 1))
            div(0:nz, :, nx-1) = f12 * dxi(1) * (f(0:nz, :, 0, 1) - f(0:nz, :, nx-2, 1))

            ! calculate df/dy with central differencing
            do i = 1, ny-2
                df(0:nz, i, :) = f12 * dxi(2) * (f(0:nz, i+1, :, 2) - f(0:nz, i-1, :, 2))
            enddo
            df(0:nz, 0,    :) = f12 * dxi(2) * (f(0:nz, 1, :, 2) - f(0:nz, ny-1, :, 2))
            df(0:nz, ny-1, :) = f12 * dxi(2) * (f(0:nz, 0, :, 2) - f(0:nz, ny-2, :, 2))

            div = div + df

            ! calculate df/dz with central differencing
            do i = 0, nz
                df(i, :, :) = f12 * dxi(3) * (f(i+1, :, :, 3) - f(i-1, :, :, 3))
            enddo

            div = div + df

        end subroutine divergence

end module inversion_mod
