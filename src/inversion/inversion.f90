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
        ! returns the associated velocity field (velog).
        subroutine vor2vel(svortg,  velog, svelog)
            double precision, intent(in)    :: svortg(0:nz, 0:nx-1, 0:ny-1, 3)
            double precision, intent(out)   :: velog(-1:nz+1, 0:ny-1, 0:nx-1, 3)
            double precision, intent(out)   :: svelog(0:nz, 0:nx-1, 0:ny-1, 3)
            double precision                :: as(0:nz, 0:nx-1, 0:ny-1)
            double precision                :: bs(0:nz, 0:nx-1, 0:ny-1)
            double precision                :: ds(0:nz, 0:nx-1, 0:ny-1)
            double precision                :: es(0:nz, 0:nx-1, 0:ny-1)
            double precision                :: ubar(0:nz), vbar(0:nz)
            double precision                :: uavg, vavg
            integer                         :: iz

            call start_timer(vor2vel_timer)

            !Define horizontally-averaged flow by integrating horizontal vorticity:
            ubar(0) = zero
            vbar(0) = zero
            do iz = 0, nz-1
                ubar(iz+1) = ubar(iz) + dz2 * (svortg(iz, 0, 0, 2) + svortg(iz+1, 0, 0, 2))
                vbar(iz+1) = vbar(iz) - dz2 * (svortg(iz, 0, 0, 1) + svortg(iz+1, 0, 0, 1))
            enddo

            ! remove the mean value to have zero net momentum
            uavg = sum(ubar(1:nz-1) + f12 * ubar(nz)) / dble(nz)
            vavg = sum(vbar(1:nz-1) + f12 * vbar(nz)) / dble(nz)
            do iz = 0, nz
                ubar(iz) = ubar(iz) - uavg
                vbar(iz) = vbar(iz) - vavg
            enddo

            !Form source term for inversion of vertical velocity:
            call diffy(svortg(:, :, :, 1), ds)
            call diffx(svortg(:, :, :, 2), es)

            !$omp parallel
            !$omp workshare
            ds = ds - es
            !$omp end workshare
            !$omp end parallel

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
            call diffx(svortg(:, :, :, 3), bs)

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

            call stop_timer(vor2vel_timer)

        end subroutine


        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Compute the gridded vorticity tendency: (excluding buoyancy effects)
        subroutine vorticity_tendency(svortg, velog, svtend)
            double precision, intent(in)  :: svortg(0:nz, 0:nx-1, 0:ny-1, 3)
            double precision, intent(in)  :: velog(0:nz, 0:ny-1, 0:nx-1, 3)
            double precision, intent(out) :: svtend(0:nz, 0:nx-1, 0:ny-1, 3)
            double precision              :: xs(0:nz, 0:nx-1, 0:ny-1)
            double precision              :: ys(0:nz, 0:nx-1, 0:ny-1)
            double precision              :: zs(0:nz, 0:nx-1, 0:ny-1)
            double precision              :: xp(0:nz, 0:ny-1, 0:nx-1)
            double precision              :: yp(0:nz, 0:ny-1, 0:nx-1)
            double precision              :: zp(0:nz, 0:ny-1, 0:nx-1)
            integer                       :: nc

            call start_timer(vtend_timer)
            
            do nc = 1, 3
                as = svortg(:, :, :, nc)
                call fftxys2p(as, vortg(:, :, :, nc)
                vortg(:, :, :, nc) = vortg(:, :, :, nc) + f_cor(nc)
            enddo
            
            ! x-component of vorticity tendency:
            yp = vortg(:, :, :, 2) * velog(:, :, :, 1) - velog(:, :, :, 2) * vortg(:, :, :, 1)   ! eta * u - v * xi
            zp = vortg(:, :, :, 3) * velog(:, :, :, 1) - velog(:, :, :, 3) * vortg(:, :, :, 1)   ! zeta * u - w * xi
            
            call fftxys2p(yp, ys)
            call fftxys2p(zp, zs)
            
            call diffy(ys, xs)
            call diffz(zs, svtend(:, :, :, 1))
            svtend(:, :, :, 1) = svtend(:, :, :, 1) + xs
            
            ! y-component of vorticity tendency:
            xp = vortg(:, :, :, 1) * velog(:, :, :, 2) - velog(:, :, :, 1) * vortg(:, :, :, 2)  ! xi * v - u * eta
            zp = vortg(:, :, :, 3) * velog(:, :, :, 2) - velog(:, :, :, 3) * vortg(:, :, :, 2)  ! zeta * v - w * eta
            
            call fftxys2p(xp, xs)
            call fftxys2p(zp, zs)
            
            call diffx(xs, ys)
            call diffz(zs, svtend(:, :, :, 2))
            svtend(:, :, :, 2) = svtend(:, :, :, 2) + ys
            
            ! z-component of vorticity tendency:
            xp = vortg(:, :, :, 1) * velog(:, :, :, 3) - velog(:, :, :, 1) * vortg(:, :, :, 3) ! xi * w - u * zeta
            yp = vortg(:, :, :, 2) * velog(:, :, :, 3) - velog(:, :, :, 2) * vortg(:, :, :, 3) ! eta * w - v * zeta
            
            call fftxys2p(xp, xs)
            call fftxys2p(yp, ys)
            
            call diffx(xs, zs)
            call diffy(ys, svtend(:, :, :, 3))
            svtend(:, :, :, 3) = svtend(:, :, :, 3) + zs

            call stop_timer(vtend_timer)

        end subroutine vorticity_tendency

end module inversion_mod
