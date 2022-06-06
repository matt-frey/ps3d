module inversion_mod
    use inversion_utils
    use parameters, only : nx, ny, nz
    use physics, only : f_cor
    use constants, only : zero, two, f12
    use sta2dfft, only : dct, dst
    use timer, only : start_timer, stop_timer
    use fields
    implicit none

    integer :: vor2vel_timer,   &
               vtend_timer

    contains

        ! Given the vorticity vector field (svor) in spectral space, this
        ! returns the associated velocity field (vel) as well as vorticity
        ! in physical space (vor)
        subroutine vor2vel
            double precision :: as(0:nz, 0:nx-1, 0:ny-1)         ! semi-spectral
            double precision :: bs(0:nz, 0:nx-1, 0:ny-1)         ! semi-spectral
            double precision :: ds(0:nz, 0:nx-1, 0:ny-1)         ! semi-spectral
            double precision :: es(0:nz, 0:nx-1, 0:ny-1)         ! semi-spectral
            double precision :: ubar(0:nz), vbar(0:nz)
            double precision :: uavg, vavg
            double precision :: wtop(0:nx-1, 0:ny-1)
            double precision :: wbot(0:nx-1, 0:ny-1)
            integer          :: iz, nc, kx, ky

            call start_timer(vor2vel_timer)

            !Combine vorticity in physical space:
            do nc = 1, 3
                call field_combine(svor(:, :, :, nc), vor(:, :, :, nc))
            enddo

            !----------------------------------------------------------
            !Form source term for inversion of vertical velocity -> ds:
            call diffy(svor(:, :, :, 1), ds)
            call diffx(svor(:, :, :, 2), es)
            !$omp parallel
            !$omp workshare
            ds = ds - es
            !$omp end workshare
            !$omp end parallel

            !Invert Laplacian to find "w_p", the particular solution:
            !$omp parallel
            !$omp workshare
            ds(1:nz-1, :, :) = green * ds(1:nz-1, :, :)
            !$omp end workshare
            !$omp end parallel


            ! Calculate d/dz of sine-series part of vertical velocity
            do iz = 1, nz-1
                as(iz, :, :) = rkz(iz) * ds(iz, :, :)
            enddo

            !Store boundary values of "w_p" to correct "w" below:
            wbot = ds(0,  :, :)
            wtop = ds(nz, :, :)

            !FFT back the particular solution to semi-spectral space:
            do ky = 0, ny-1
                do kx = 0, nx-1
                    call dct(1, nz-1, as(1:nz-1, kx, ky), ztrig, zfactors)
                    call dst(1, nz,   ds(1:nz, kx, ky), ztrig, zfactors)
                enddo
            enddo

            ! Calculate the linear part (bs) and its derivative (es) of vertical velocity in semi-spectral space
            do iz = 1, nz-1
                bs(iz, :, :) =   wbot *  psi(nz-iz, :, :) + wtop *  psi(iz, :, :)
                es(iz, :, :) = - wbot * dpsi(nz-iz, :, :) + wtop * dpsi(iz, :, :)
            enddo

            ! Combine vertical velocity (ds) and its derivative (es) given the sine and linear parts:
            ds(1:nz-1, :, :) = ds(1:nz-1, :, :) + bs(1:nz-1, :, :)
            es(1:nz-1, :, :) = es(1:nz-1, :, :) + as(1:nz-1, :, :)

            !----------------------------------------------------------------------
            !Define horizontally-averaged flow by integrating horizontal vorticity:
            ubar(0) = zero
            vbar(0) = zero
            do iz = 0, nz-1
                ubar(iz+1) = ubar(iz) + dz2 * (svor(iz, 0, 0, 2) + svor(iz+1, 0, 0, 2))
                vbar(iz+1) = vbar(iz) - dz2 * (svor(iz, 0, 0, 1) + svor(iz+1, 0, 0, 1))
            enddo

            !Remove the mean value to have zero net horizontal momentum:
            uavg = sum(ubar(1:nz-1) + f12 * ubar(nz)) / dble(nz)
            vavg = sum(vbar(1:nz-1) + f12 * vbar(nz)) / dble(nz)
            do iz = 0, nz
                ubar(iz) = ubar(iz) - uavg
                vbar(iz) = vbar(iz) - vavg
            enddo

            !-------------------------------------------------------
            !Find x velocity component "u":
            call diffx(es, as)
            call diffy(svor(:, :, :, 3), bs)

            !$omp parallel do
            do iz = 0, nz
                as(iz, :, :) = k2l2i * (as(iz, :, :) + bs(iz, :, :))
            enddo
            !$omp end parallel do

            !Add horizontally-averaged flow:
            as(:, 0, 0) = ubar

            !Store spectral form of "u":
            svel(:, :, :, 1) = as

            !Get "u" in physical space:
            call fftxys2p(as, vel(0:nz, :, :, 1))

            !-------------------------------------------------------
            !Find y velocity component "v":
            call diffy(es, as)
            call diffx(svor(:, :, :, 3), bs)

            !$omp parallel do
            do iz = 0, nz
                as(iz, :, :) = k2l2i * (as(iz, :, :) - bs(iz, :, :))
            enddo
            !$omp end parallel do

            !Add horizontally-averaged flow:
            as(:, 0, 0) = vbar

            !Store spectral form of "v":
            svel(:, :, :, 2) = as

            !Get "v" in physical space:
            call fftxys2p(as, vel(0:nz, :, :, 2))

            !-------------------------------------------------------
            !Store spectral form of "w":
            svel(:, :, :, 3) = ds

            !Get "w" in physical space:
            call fftxys2p(ds, vel(0:nz, :, :, 3))

            call stop_timer(vor2vel_timer)

        end subroutine


        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Compute the gridded vorticity tendency: (excluding buoyancy effects)
        subroutine vorticity_tendency
            double precision :: xs(0:nz, 0:nx-1, 0:ny-1)
            double precision :: ys(0:nz, 0:nx-1, 0:ny-1)
            double precision :: zs(0:nz, 0:nx-1, 0:ny-1)
            double precision :: xp(0:nz, 0:ny-1, 0:nx-1)
            double precision :: yp(0:nz, 0:ny-1, 0:nx-1)
            double precision :: zp(0:nz, 0:ny-1, 0:nx-1)
            integer          :: nc

            call start_timer(vtend_timer)

            !-------------------------------------------------------
            ! First store absolute vorticity in physical space:
            do nc = 1, 3
                call field_combine(svor(:, :, :, nc), vor(:, :, :, nc))
                vor(:, :, :, nc) = vor(:, :, :, nc) + f_cor(nc)
            enddo

            !-------------------------------------------------------
            ! x-component of vorticity tendency:
            yp = vor(:, :, :, 2) * vel(:, :, :, 1) - vel(:, :, :, 2) * vor(:, :, :, 1)   ! eta * u - v * xi
            zp = vor(:, :, :, 3) * vel(:, :, :, 1) - vel(:, :, :, 3) * vor(:, :, :, 1)   ! zeta * u - w * xi

            call fftxyp2s(yp, ys)
            call fftxyp2s(zp, zs)

            call diffy(ys, xs)
            call diffz(zs, svtend(:, :, :, 1))

            ! Extrapolate
            svtend(0,  :, :, 1) = two * svtend(1,    :, :, 1) - svtend(2,    :, :, 1)
            svtend(nz, :, :, 1) = two * svtend(nz-1, :, :, 1) - svtend(nz-2, :, :, 1)

            svtend(:, :, :, 1) = svtend(:, :, :, 1) + xs

            !-------------------------------------------------------
            ! y-component of vorticity tendency:
            xp = vor(:, :, :, 1) * vel(:, :, :, 2) - vel(:, :, :, 1) * vor(:, :, :, 2)  ! xi * v - u * eta
            zp = vor(:, :, :, 3) * vel(:, :, :, 2) - vel(:, :, :, 3) * vor(:, :, :, 2)  ! zeta * v - w * eta

            call fftxyp2s(xp, xs)
            call fftxyp2s(zp, zs)

            call diffx(xs, ys)
            call diffz(zs, svtend(:, :, :, 2))

            ! Extrapolate
            svtend(0,  :, :, 2) = two * svtend(1,    :, :, 2) - svtend(2,    :, :, 2)
            svtend(nz, :, :, 2) = two * svtend(nz-1, :, :, 2) - svtend(nz-2, :, :, 2)

            svtend(:, :, :, 2) = svtend(:, :, :, 2) + ys

            !-------------------------------------------------------
            ! z-component of vorticity tendency:
            xp = vor(:, :, :, 1) * vel(:, :, :, 3) - vel(:, :, :, 1) * vor(:, :, :, 3) ! xi * w - u * zeta
            yp = vor(:, :, :, 2) * vel(:, :, :, 3) - vel(:, :, :, 2) * vor(:, :, :, 3) ! eta * w - v * zeta

            call fftxyp2s(xp, xs)
            call fftxyp2s(yp, ys)

            call diffx(xs, zs)
            call diffy(ys, svtend(:, :, :, 3))
            svtend(:, :, :, 3) = svtend(:, :, :, 3) + zs


            !-------------------------------------------------------
            ! Apply z-filter:
            do nc = 1, 3
                call apply_zfilter(svtend(:, :, :, nc))
            enddo

            call stop_timer(vtend_timer)

        end subroutine vorticity_tendency

end module inversion_mod
