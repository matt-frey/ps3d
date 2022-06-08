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
            double precision :: cs(0:nz, 0:nx-1, 0:ny-1)         ! semi-spectral
            double precision :: ubar(0:nz), vbar(0:nz)
            integer          :: iz, nc, kx, ky, kz

            call start_timer(vor2vel_timer)

            !Filter the vorticity and combine vorticity in physical space:
            do nc = 1, 3
                svor(:, :, :, nc) = filt * svor(:, :, :, nc)
                call field_combine_physical(svor(:, :, :, nc), vor(:, :, :, nc))
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

            !Calculate the boundary contributions of the source to the vertical velocity (bs)
            !and its derivative (es) in semi-spectral space:
            do iz = 1, nz-1
                bs(iz, :, :) = ds(nz, :, :) *  psi(iz, :, :) + ds(0, :, :) *  psi(nz-iz, :, :)
            enddo
            do iz = 0, nz
                es(iz, :, :) = ds(nz, :, :) * dpsi(iz, :, :) - ds(0, :, :) * dpsi(nz-iz, :, :)
            enddo

            !Invert Laplacian to find the part of w expressible as a sine series:
            !$omp parallel
            !$omp workshare
            ds(1:nz-1, :, :) = green * ds(1:nz-1, :, :)
            !$omp end workshare
            !$omp end parallel

            ! Calculate d/dz of this sine series:
            as(0, :, :) = zero
            do kz = 1, nz-1
                as(kz, :, :) = rkz(kz) * ds(kz, :, :)
            enddo
            as(nz, :, :) = zero

            !FFT these quantities back to semi-spectral space:
            do ky = 0, ny-1
                do kx = 0, nx-1
                    call dct(1, nz, as(0:nz, kx, ky), ztrig, zfactors)
                    call dst(1, nz, ds(1:nz, kx, ky), ztrig, zfactors)
                enddo
            enddo

            ! Combine vertical velocity (ds) and its derivative (es) given the sine and linear parts:
            ds(0     , :, :) = zero
            ds(1:nz-1, :, :) = ds(1:nz-1, :, :) + bs(1:nz-1, :, :)
            ds(nz    , :, :) = zero
            es = es + as
            ! OMP the above????

            ! Get complete zeta field in semi-spectral space
            cs = svor(:, :, :, 3)
            call field_combine_semi_spectral(cs)

            !----------------------------------------------------------------------
            !Define horizontally-averaged flow by integrating the horizontal vorticity:

            !First integrate the sine series in svor(1:nz-1, 0, 0, 1 & 2):
            ubar(0) = zero
            vbar(0) = zero
            ubar(1:nz-1) = -rkzi * svor(1:nz-1, 0, 0, 2)
            vbar(1:nz-1) =  rkzi * svor(1:nz-1, 0, 0, 1)
            ubar(nz) = zero
            vbar(nz) = zero

            !Transform to semi-spectral space as a cosine series:
            call dct(1, nz, ubar, ztrig, zfactors)
            call dct(1, nz, vbar, ztrig, zfactors)

            !Add contribution from the linear function connecting the boundary values:
            ubar = ubar + svor(nz, 0, 0, 2) * gamtop - svor(0, 0, 0, 2) * gambot
            vbar = vbar - svor(nz, 0, 0, 1) * gamtop + svor(0, 0, 0, 1) * gambot

            !-------------------------------------------------------
            !Find x velocity component "u":
            call diffx(es, as)
            call diffy(cs, bs)

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
            call fftxys2p(as, vel(:, :, :, 1))

            !-------------------------------------------------------
            !Find y velocity component "v":
            call diffy(es, as)
            call diffx(cs, bs)

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
            call fftxys2p(as, vel(:, :, :, 2))

            !-------------------------------------------------------
            !Store spectral form of "w":
            svel(:, :, :, 3) = ds

            !Get "w" in physical space:
            call fftxys2p(ds, vel(:, :, :, 3))

            call stop_timer(vor2vel_timer)

        end subroutine


        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Compute the gridded vorticity tendency: (excluding buoyancy effects)
        subroutine vorticity_tendency
            double precision :: fp(0:nz, 0:ny-1, 0:nx-1)    ! physical space
            double precision :: gp(0:nz, 0:ny-1, 0:nx-1)    ! physical space
            double precision :: p(0:nz, 0:nx-1, 0:ny-1)     ! mixed spectral space
            double precision :: q(0:nz, 0:nx-1, 0:ny-1)     ! mixed spectral space
            double precision :: r(0:nz, 0:nx-1, 0:ny-1)     ! mixed spectral space
            integer          :: nc

            call start_timer(vtend_timer)

            !-------------------------------------------------------
            ! First store absolute vorticity in physical space:
            do nc = 1, 3
                call field_combine_physical(svor(:, :, :, nc), vor(:, :, :, nc))
                vor(:, :, :, nc) = vor(:, :, :, nc) + f_cor(nc)
            enddo

            call field_combine_physical(sbuoy, buoy)

            !-------------------------------------------------------
            ! Tendency in flux form:
            !   dxi/dt  = dr/dy - dq/dz
            !   deta/dt = dp/dz - dr/dx
            !  dzeta/dt = dq/dx - dp/dy

            ! r = u * eta - v * xi + b
            fp = vel(:, :, :, 1) * vor(:, :, :, 2) - vel(:, :, :, 2) * vor(:, :, :, 1) + buoy
            call field_decompose_physical(fp, r)

            ! q = w * xi - u * zeta
            fp = vel(:, :, :, 3) * vor(:, :, :, 1) - vel(:, :, :, 1) * vor(:, :, :, 3)
            call field_decompose_physical(fp, q)

            ! dxi/dt  = dr/dy - dq/dz
            call diffy(r, svtend(:, :, :, 1))
            call diffz(fp, gp)
            call field_decompose_physical(gp, p)
            svtend(:, :, :, 1) = svtend(:, :, :, 1) - p     ! here: p = dq/dz

            ! p = v * zeta - w * eta
            fp = vel(:, :, :, 2) * vor(:, :, :, 3) - vel(:, :, :, 3) * vor(:, :, :, 2)
            call field_decompose_physical(fp, p)

            ! deta/dt = dp/dz - dr/dx
            call diffx(r, svtend(:, :, :, 2))
            call diffz(fp, gp)
            call field_decompose_physical(gp, r)
            svtend(:, :, :, 2) = r - svtend(:, :, :, 2)     ! here: r = dp/dz

            ! dzeta/dt = dq/dx - dp/dy
            call diffx(q, svtend(:, :, :, 3))
            call diffy(p, r)                                ! here: r = dp/dy
            svtend(:, :, :, 3) = svtend(:, :, :, 3) - r

            call stop_timer(vtend_timer)

        end subroutine vorticity_tendency

end module inversion_mod
