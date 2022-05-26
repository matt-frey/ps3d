module inversion_mod
    use inversion_utils
    use parameters, only : nx, ny, nz, dxi
    use physics, only : f_cor
    use constants, only : zero, two, f12
    use sta2dfft, only : dct, dst
    use timer, only : start_timer, stop_timer
    implicit none

    integer :: vor2vel_timer,   &
               vtend_timer

    contains

        ! Given the vorticity vector field (svortg) in spectral space, this
        ! returns the associated velocity field (velog) as well as vorticity
        ! in physical space (vortg)
        subroutine vor2vel(svortg, vortg,  svelog, velog)
            double precision, intent(out)   :: vortg(0:nz, 0:ny-1, 0:nx-1, 3)
            double precision, intent(in)    :: svortg(0:nz, 0:nx-1, 0:ny-1, 3)  ! fully spectral
            double precision, intent(out)   :: velog(0:nz, 0:ny-1, 0:nx-1, 3)
            double precision, intent(out)   :: svelog(0:nz, 0:nx-1, 0:ny-1, 3)  ! semi-spectral
            double precision                :: as(0:nz, 0:nx-1, 0:ny-1)         ! semi-spectral
            double precision                :: bs(0:nz, 0:nx-1, 0:ny-1)         ! semi-spectral
            double precision                :: cs(0:nz, 0:nx-1, 0:ny-1)         ! semi-spectral
            double precision                :: ds(0:nz, 0:nx-1, 0:ny-1)
            double precision                :: es(0:nz, 0:nx-1, 0:ny-1)
            double precision                :: ss(1:nz, 0:nx-1, 0:ny-1)         ! sine transform in z
            double precision                :: ubar(0:nz), vbar(0:nz)
            double precision                :: uavg, vavg
            double precision                :: wtop(0:nx-1, 0:ny-1), wbot(0:nx-1, 0:ny-1)
            integer                         :: iz, nc, kx, ky, kz

            call start_timer(vor2vel_timer)

            !Compute vorticity in physical space:
            do nc = 1, 3
                as = svortg(:, :, :, nc)
                call fftczs2p(as, vortg(:, :, :, nc))
            enddo

            !Form source term for inversion of vertical velocity:
            call diffy(svortg(:, :, :, 1), ds)
            call diffx(svortg(:, :, :, 2), es)

            !$omp parallel
            !$omp workshare
            ds = green * (ds - es)
            !$omp end workshare
            !$omp end parallel

            ! FFT back the particular solution to semi-spectral space:
            do ky = 0, ny-1
                do kx = 0, nx-1
                    call dct(1, nz, ds(:, kx, ky), ztrig, zfactors)
                enddo
            enddo

            wbot = ds(0,  :, :)
            wtop = ds(nz, :, :)

            ! Define the complete vertical velocity in semi-spectral space:
            do iz = 1, nz-1
                ds(iz, :, :) = ds(iz, :, :) - (wbot * decz(nz-iz, :, :) + wtop * decz(iz, :, :))
            enddo
            ds(0,  :, :) = zero
            ds(nz, :, :) = zero

!            ! FFT to fully spectral space (sine transform) as the array ss:
!            do ky = 0, ny-1
!                do kx = 0, nx-1
!                    ss(:, kx, ky) = ds(1:nz, kx, ky)
!                    call dst(1, nz, ss(:, kx, ky), ztrig, zfactors)
!                enddo
!            enddo
!
!            ! Derivative in z (dw/dz in fully spectral space):
!            do kz = 1, nz
!                es(kz, :, :) = rkz(kz) * ss(kz, :, :)
!            enddo
            
!            es(0,  :, :) = zero
!            es(nz, :, :) = zero

            call diffz(ds, es)

            es(0,  :, :) = zero
            es(nz, :, :) = zero

            !! FFT back to semi-spectral space:
            !do ky = 0, ny-1
            !    do kx = 0, nx-1
            !        call dct(1, nz, es(:, kx, ky), ztrig, zfactors)
            !    enddo
            ! enddo

            !print *, 'x-svort', minval(abs(svortg(:, :, :, 1)))
            !print *, 'y-svort', minval(abs(svortg(:, :, :, 2)))
            !print *, 'z-svort', minval(abs(svortg(:, :, :, 3))) 

            !Convert vorticity components to semi-spectral space
            call fftfs2ss(svortg(:, :, :, 1), as)
            call fftfs2ss(svortg(:, :, :, 2), bs)
            call fftfs2ss(svortg(:, :, :, 3), cs)

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

            call stop_timer(vor2vel_timer)

        end subroutine


        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Compute the gridded vorticity tendency: (excluding buoyancy effects)
        subroutine vorticity_tendency(svortg, velog, vortg, svtend)
            double precision, intent(in)  :: svortg(0:nz, 0:nx-1, 0:ny-1, 3)    ! fully sectral
            double precision, intent(in)  :: velog(0:nz, 0:ny-1, 0:nx-1, 3)
            double precision, intent(out) :: vortg(0:nz, 0:ny-1, 0:nx-1, 3)
            double precision, intent(out) :: svtend(0:nz, 0:nx-1, 0:ny-1, 3)    ! fully spectral
            double precision              :: xs(0:nz, 0:nx-1, 0:ny-1)
            double precision              :: ys(0:nz, 0:nx-1, 0:ny-1)
            double precision              :: zs(0:nz, 0:nx-1, 0:ny-1)
            double precision              :: xp(0:nz, 0:ny-1, 0:nx-1)
            double precision              :: yp(0:nz, 0:ny-1, 0:nx-1)
            double precision              :: zp(0:nz, 0:ny-1, 0:nx-1)
            integer                       :: nc

            call start_timer(vtend_timer)

            do nc = 1, 3
                xs = svortg(:, :, :, nc)
                call fftczs2p(xs, vortg(:, :, :, nc))
                vortg(:, :, :, nc) = vortg(:, :, :, nc) + f_cor(nc)
            enddo

            ! x-component of vorticity tendency:
            yp = vortg(:, :, :, 2) * velog(:, :, :, 1) - velog(:, :, :, 2) * vortg(:, :, :, 1)   ! eta * u - v * xi
            zp = vortg(:, :, :, 3) * velog(:, :, :, 1) - velog(:, :, :, 3) * vortg(:, :, :, 1)   ! zeta * u - w * xi

            call fftxyp2s(yp, ys)
            call fftxyp2s(zp, zs)

            call diffy(ys, xs)
            call diffz(zs, svtend(:, :, :, 1))
            svtend(:, :, :, 1) = svtend(:, :, :, 1) + xs

            ! y-component of vorticity tendency:
            xp = vortg(:, :, :, 1) * velog(:, :, :, 2) - velog(:, :, :, 1) * vortg(:, :, :, 2)  ! xi * v - u * eta
            zp = vortg(:, :, :, 3) * velog(:, :, :, 2) - velog(:, :, :, 3) * vortg(:, :, :, 2)  ! zeta * v - w * eta

            call fftxyp2s(xp, xs)
            call fftxyp2s(zp, zs)

            call diffx(xs, ys)
            call diffz(zs, svtend(:, :, :, 2))
            svtend(:, :, :, 2) = svtend(:, :, :, 2) + ys

            ! z-component of vorticity tendency:
            xp = vortg(:, :, :, 1) * velog(:, :, :, 3) - velog(:, :, :, 1) * vortg(:, :, :, 3) ! xi * w - u * zeta
            yp = vortg(:, :, :, 2) * velog(:, :, :, 3) - velog(:, :, :, 2) * vortg(:, :, :, 3) ! eta * w - v * zeta

            call fftxyp2s(xp, xs)
            call fftxyp2s(yp, ys)

            call diffx(xs, zs)
            call diffy(ys, svtend(:, :, :, 3))
            svtend(:, :, :, 3) = svtend(:, :, :, 3) + zs

            !Convert vorticity tendency from semi-spectral to fully-spectral space
            do nc = 1, 3
                xs = svtend(:, :, :, nc)
                call fftss2fs(xs, svtend(:, :, :, nc))
            enddo

            call stop_timer(vtend_timer)

        end subroutine vorticity_tendency

end module inversion_mod
