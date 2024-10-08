module inversion_mod
    use inversion_utils
    use parameters, only : nx, ny, nz
    use physics, only : f_cor
#ifdef ENABLE_BUOYANCY
    use physics, only : bfsq
#endif
    use constants, only : zero, two
    use sta2dfft, only : dct, dst
    use sta3dfft, only : rkz, rkzi, ztrig, zfactors, diffx, diffy, fftxyp2s, fftxys2p
    use mpi_timer, only : start_timer, stop_timer
    use fields
    implicit none

    integer :: vor2vel_timer,   &
               vtend_timer,     &
               pres_timer

    contains

        ! Given the vorticity vector field (svor) in spectral space, this
        ! returns the associated velocity field (vel) as well as vorticity
        ! in physical space (vor)
        subroutine vor2vel
            double precision :: as(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))  ! semi-spectral
            double precision :: bs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))  ! semi-spectral
            double precision :: ds(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))  ! semi-spectral
            double precision :: es(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))  ! semi-spectral
            double precision :: cs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))  ! semi-spectral
            double precision :: ubar(0:nz), vbar(0:nz)
            integer          :: iz, nc, kx, ky, kz

            call start_timer(vor2vel_timer)

            !----------------------------------------------------------
            ! Enforce solenoidality
            ! A, B, C are vorticities
            ! D = B_x - A_y; E = C_z
            ! A = k2l2i * (E_x + D_y) and B = k2l2i * (E_y - D_x) --> A_x + B_y + C_z = zero
            call diffx(svor(:, :, :, 2), as) ! as = B_x
            call diffy(svor(:, :, :, 1), bs) ! bs = A_y
            !$omp parallel workshare
            ds = as - bs                     ! ds = D
            cs = svor(:, :, :, 3)
            !$omp end parallel workshare
            call field_combine_semi_spectral(cs)
            call central_diffz(cs, es)                     ! es = E
            call field_decompose_semi_spectral(es)

            ! ubar and vbar are used here to store the mean x and y components of the vorticity
            if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
                ubar = svor(:, 0, 0, 1)
                vbar = svor(:, 0, 0, 2)
            endif

            call diffx(es, svor(:, :, :, 1)) ! E_x
            call diffy(ds, cs)               ! cs = D_y
            !$omp parallel do private(iz)  default(shared)
            do iz = 0, nz
               svor(iz, :, :, 1) = k2l2i * (svor(iz, :, :, 1) + cs(iz, :, :))
            enddo
            !$omp end parallel do

            call diffy(es, svor(:, :, :, 2)) ! E_y
            call diffx(ds, cs)               ! D_x

            !$omp parallel do private(iz)  default(shared)
            do iz = 0, nz
               svor(iz, :, :, 2) = k2l2i * (svor(iz, :, :, 2) - cs(iz, :, :))
            enddo
            !$omp end parallel do

            ! bring back the mean x and y components of the vorticity
            if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
                svor(:, 0, 0, 1) = ubar
                svor(:, 0, 0, 2) = vbar
            endif

            !----------------------------------------------------------
            !Combine vorticity in physical space:
            do nc = 1, 3
                call field_combine_physical(svor(:, :, :, nc), vor(:, :, :, nc))
            enddo

            !----------------------------------------------------------
            !Form source term for inversion of vertical velocity -> ds:
            call diffy(svor(:, :, :, 1), ds)
            call diffx(svor(:, :, :, 2), es)
            !$omp parallel workshare
            ds = ds - es
            !$omp end parallel workshare

            !Calculate the boundary contributions of the source to the vertical velocity (bs)
            !and its derivative (es) in semi-spectral space:
            !$omp parallel do private(iz)  default(shared)
            do iz = 1, nz-1
                bs(iz, :, :) = ds(0, :, :) *  thetam(iz, :, :) + ds(nz, :, :) *  thetap(iz, :, :)
            enddo
            !$omp end parallel do

            !$omp parallel do private(iz)  default(shared)
            do iz = 0, nz
                es(iz, :, :) = ds(0, :, :) * dthetam(iz, :, :) + ds(nz, :, :) * dthetap(iz, :, :)
            enddo
            !$omp end parallel do

            !Invert Laplacian to find the part of w expressible as a sine series:
            !$omp parallel workshare
            ds(1:nz-1, :, :) = green(1:nz-1, :, :) * ds(1:nz-1, :, :)
            !$omp end parallel workshare

            ! Calculate d/dz of this sine series:
            !$omp parallel workshare
            as(0, :, :) = zero
            !$omp end parallel workshare
            !$omp parallel do private(iz)  default(shared)
            do kz = 1, nz-1
                as(kz, :, :) = rkz(kz) * ds(kz, :, :)
            enddo
            !$omp end parallel do
            !$omp parallel workshare
            as(nz, :, :) = zero
            !$omp end parallel workshare

            !FFT these quantities back to semi-spectral space:
            !$omp parallel do collapse(2) private(kx, ky)
            do kx = box%lo(1), box%hi(1)
                do ky = box%lo(2), box%hi(2)
                    call dct(1, nz, as(0:nz, ky, kx), ztrig, zfactors)
                    call dst(1, nz, ds(1:nz, ky, kx), ztrig, zfactors)
                enddo
            enddo
            !$omp end parallel do

            ! Combine vertical velocity (ds) and its derivative (es) given the sine and linear parts:
            !$omp parallel workshare
            ds(0     , :, :) = zero
            ds(1:nz-1, :, :) = ds(1:nz-1, :, :) + bs(1:nz-1, :, :)
            ds(nz    , :, :) = zero
            es = es + as

            ! Get complete zeta field in semi-spectral space
            cs = svor(:, :, :, 3)
            !$omp end parallel workshare
            call field_combine_semi_spectral(cs)

            !----------------------------------------------------------------------
            !Define horizontally-averaged flow by integrating the horizontal vorticity:

            !First integrate the sine series in svor(1:nz-1, 0, 0, 1 & 2):
            if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
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
            endif

            !-------------------------------------------------------
            !Find x velocity component "u":
            call diffx(es, as)
            call diffy(cs, bs)

            !$omp parallel do private(iz) default(shared)
            do iz = 0, nz
                as(iz, :, :) = k2l2i * (as(iz, :, :) + bs(iz, :, :))
            enddo
            !$omp end parallel do

            !Add horizontally-averaged flow:
            if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
                as(:, 0, 0) = ubar
            endif

            !Store spectral form of "u":
            !$omp parallel workshare
            svel(:, :, :, 1) = as
            !$omp end parallel workshare

            !Get "u" in physical space:
            call fftxys2p(as, vel(:, :, :, 1))

            !-------------------------------------------------------
            !Find y velocity component "v":
            call diffy(es, as)
            call diffx(cs, bs)

            !$omp parallel do private(iz) default(shared)
            do iz = 0, nz
                as(iz, :, :) = k2l2i * (as(iz, :, :) - bs(iz, :, :))
            enddo
            !$omp end parallel do

            !Add horizontally-averaged flow:
            if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
                as(:, 0, 0) = vbar
            endif

            !Store spectral form of "v":
            !$omp parallel workshare
            svel(:, :, :, 2) = as
            !$omp end parallel workshare

            !Get "v" in physical space:
            call fftxys2p(as, vel(:, :, :, 2))

            !-------------------------------------------------------
            !Store spectral form of "w":
            !$omp parallel workshare
            svel(:, :, :, 3) = ds
            !$omp end parallel workshare

            !Get "w" in physical space:
            call fftxys2p(ds, vel(:, :, :, 3))

            call stop_timer(vor2vel_timer)

        end subroutine


        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#ifdef ENABLE_BUOYANCY
        ! Compute the gridded buoyancy tendency (using flux form):
        subroutine buoyancy_tendency
            double precision :: fp(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)) ! flux component in phys space
            double precision :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)) ! flux component in spec space
            double precision :: ds(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)) ! buoyancy deriv in spec space
            double precision :: btend(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)) ! buoy source in phys space

            !--------------------------------------------------------------
            ! We calculate the buoyancy source b_t = -(u,v,w)*grad(b)
            ! in flux form b_t = - div(F) where F = (u*b, v*b, w*b);
            ! hence div(F) = u*b_x + v*b_y + w*b_z + b * (u_x + v_y + w_z)
            ! but u_x + v_y + w_z = 0 as we assume incompressibility.

            call field_combine_physical(sbuoy, buoy)

            ! Define the x-component of the flux
            fp = vel(:, :, :, 1) * buoy

            ! Differentiate
            call field_decompose_physical(fp, fs)
            call diffx(fs, ds)
            call field_combine_physical(ds, btend)

            ! Define the y-component of the flux
            fp = vel(:, :, :, 2) * buoy

            call field_decompose_physical(fp, fs)
            call diffy(fs, ds)
            call field_combine_physical(ds, fp)

            btend = - btend - fp

            ! Define the z-component of the flux
            fp = vel(:, :, :, 3) * buoy

            ! Differentiate
            call field_decompose_physical(fp, fs)
            call diffz(fs, ds)
            call field_combine_physical(ds, fp)

            ! b = N^2 * z + b'
            ! db/dt = db/dz * dz/dt + db'/dt
            ! db/dt = N^2 * w + db'/dt
            ! here: buoy = b'
            ! --> we must subtract N^2 * w to get total buoyancy tendency
            ! (note: we calculate -grad(b), therefore - N^2 * w)
            ! Note: We calculate the tendency in flux form:
            !       b_t = - div(F) where F = (u*b, v*b, w*b);
            ! which is in z:
            !
            !   d(w*b)/dz = dw/dz * b + w * db/dz
            !
            ! but dw/dz = 0 as the velocity is solenoidal, hence
            !
            !   d(w*b)/dz = w * db/dz
            !
            btend = btend - bfsq * vel(:, :, :, 3) - fp

            ! Convert to mixed-spectral space:
            call field_decompose_physical(btend, sbuoys)

        end subroutine buoyancy_tendency
#endif

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Compute the gridded vorticity tendency:
        subroutine vorticity_tendency
            double precision :: fp(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))    ! physical space
            double precision :: gp(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))    ! physical space
            double precision :: p(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))     ! mixed spectral space
            double precision :: q(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))     ! mixed spectral space
            double precision :: r(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))     ! mixed spectral space
            integer          :: nc

            call start_timer(vtend_timer)

            !-------------------------------------------------------
            ! First store absolute vorticity in physical space:
            do nc = 1, 3
                !$omp parallel workshare
                vor(:, :, :, nc) = vor(:, :, :, nc) + f_cor(nc)
                !$omp end parallel workshare
            enddo

#ifdef ENABLE_BUOYANCY
            call field_combine_physical(sbuoy, buoy)
#endif
            !-------------------------------------------------------
            ! Tendency in flux form:
            !   dxi/dt  = dr/dy - dq/dz
            !   deta/dt = dp/dz - dr/dx
            !  dzeta/dt = dq/dx - dp/dy

            ! r = u * eta - v * xi + b
            !$omp parallel workshare
            fp = vel(:, :, :, 1) * vor(:, :, :, 2) - vel(:, :, :, 2) * vor(:, :, :, 1)
#ifdef ENABLE_BUOYANCY
            fp = fp + buoy
#endif
            !$omp end parallel workshare
            call field_decompose_physical(fp, r)

            ! q = w * xi - u * zeta
            !$omp parallel workshare
            fp = vel(:, :, :, 3) * vor(:, :, :, 1) - vel(:, :, :, 1) * vor(:, :, :, 3)
            !$omp end parallel workshare
            call field_decompose_physical(fp, q)

            ! dxi/dt  = dr/dy - dq/dz
            call diffy(r, svorts(:, :, :, 1))
            call central_diffz(fp, gp)
            call field_decompose_physical(gp, p)
            !$omp parallel workshare
            svorts(:, :, :, 1) = svorts(:, :, :, 1) - p     ! here: p = dq/dz
            !$omp end parallel workshare

            ! p = v * zeta - w * eta
            !$omp parallel workshare
            fp = vel(:, :, :, 2) * vor(:, :, :, 3) - vel(:, :, :, 3) * vor(:, :, :, 2)
            !$omp end parallel workshare
            call field_decompose_physical(fp, p)

            ! deta/dt = dp/dz - dr/dx
            call diffx(r, svorts(:, :, :, 2))
            call central_diffz(fp, gp)
            call field_decompose_physical(gp, r)
            !$omp parallel workshare
            svorts(:, :, :, 2) = r - svorts(:, :, :, 2)     ! here: r = dp/dz
            !$omp end parallel workshare

            ! dzeta/dt = dq/dx - dp/dy
            call diffx(q, svorts(:, :, :, 3))
            call diffy(p, r)                                ! here: r = dp/dy
            !$omp parallel workshare
            svorts(:, :, :, 3) = svorts(:, :, :, 3) - r
            !$omp end parallel workshare

            call stop_timer(vtend_timer)

        end subroutine vorticity_tendency

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Compute pressure with buoyancy anomaly b' (if --enable-buoyancy)
        ! If buoyancy mode is disabled we solve:
        ! Solves Lap(p) = R + f * zeta + b'_z with dp/dz = b' on each boundary where
        ! R = 2[J_xy(u,v) + J_yz(v,w) + J_zx(w,u)) where J_ab(f,g) = f_a * g_b - f_b * g_a
        ! is the Jacobian.
        ! If buoyancy mode is enabled we solve:
        ! Solves Lap(p) = R with dp/dz = 0 on each boundary
        subroutine pressure(dudx, dudy, dvdy, dwdx, dwdy)
            double precision, intent(in) :: dudx(0:nz, box%lo(2):box%hi(2), &
                                                       box%lo(1):box%hi(1)) ! du/dx in physical space
            double precision, intent(in) :: dudy(0:nz, box%lo(2):box%hi(2), &
                                                       box%lo(1):box%hi(1)) ! du/dy in physical space
            double precision, intent(in) :: dvdy(0:nz, box%lo(2):box%hi(2), &
                                                       box%lo(1):box%hi(1)) ! dv/dy in physical space
            double precision, intent(in) :: dwdx(0:nz, box%lo(2):box%hi(2), &
                                                       box%lo(1):box%hi(1)) ! dw/dx in physical space
            double precision, intent(in) :: dwdy(0:nz, box%lo(2):box%hi(2), &
                                                       box%lo(1):box%hi(1)) ! dw/dy in physical space
            double precision             :: dwdz(0:nz, box%lo(2):box%hi(2), &
                                                       box%lo(1):box%hi(1)) ! dw/dz in physical space
            double precision             :: rs(0:nz, box%lo(2):box%hi(2),   &
                                                     box%lo(1):box%hi(1))   ! rhs in spectral space
#ifdef ENABLE_BUOYANCY
            double precision             :: dbdz(0:nz, box%lo(2):box%hi(2), &
                                                       box%lo(1):box%hi(1))
#endif
            integer                      :: kx, ky

            call start_timer(pres_timer)

            !-------------------------------------------------------
            ! Compute rhs (and store in pres) of Poisson equation to determine the pressure:
            ! J_xy(u, v) = du/dx * dv/dy - du/dy * dv/dx
            ! J_yz(v, w) = dv/dy * dw/dz - dv/dz * dw/dy
            ! J_zx(w, u) = dw/dz * du/dx - dw/dx * du/dz
            ! dv/dx = \zeta + du/dy
            ! du/dz = \eta + dw/dx
            ! dv/dz = dw/dy - \xi
            ! dw/dz = - (du/dx + dv/dy)

            !$omp parallel workshare
            dwdz = - (dudx + dvdy)

            pres = two * (dudx * dvdy - dudy * (vor(:, :, :, 3) + dudy)  + &
                          dvdy * dwdz - dwdy * (dwdy - vor(:, :, :, 1))  + &
                          dwdz * dudx - dwdx * (dwdx + vor(:, :, :, 2)))
            !$omp end parallel workshare

#ifdef ENABLE_BUOYANCY
            call field_combine_physical(sbuoy, buoy)
            call central_diffz(buoy, dbdz)
            pres = pres + dbdz + f_cor(3) * vor(:, :, :, 3)
#endif

            !-------------------------------------------------------
            ! Transform to full spectral space:
            call fftxyp2s(pres, rs)

            !$omp parallel do collapse(2) private(kx, ky)
            do kx = box%lo(1), box%hi(1)
                do ky = box%lo(2), box%hi(2)
                    call dct(1, nz, rs(0:nz, ky, kx), ztrig, zfactors)
                enddo
            enddo
            !$omp end parallel do

            !-------------------------------------------------------
            !Invert Laplacian:
            !$omp parallel workshare
            rs = green * rs
            !$omp end parallel workshare

            !-------------------------------------------------------
            ! Transform to physical space:
            !$omp parallel do collapse(2) private(kx, ky)
            do kx = box%lo(1), box%hi(1)
                do ky = box%lo(2), box%hi(2)
                    call dct(1, nz, rs(0:nz, ky, kx), ztrig, zfactors)
                enddo
            enddo
            !$omp end parallel do

#ifdef ENABLE_BUOYANCY
            ! now in semi-spectral space, note k2l2i(0, 0) = 0
            do kx = box%lo(1), box%hi(1)
                do ky = box%lo(2), box%hi(2)
                    rs(:, ky, kx) = rs(:, ky, kx) &
                                  + sbuoy(nz, ky, kx) * k2l2i(ky, kx) * dphip(:, ky, kx) &
                                  + sbuoy(0,  ky, kx) * k2l2i(ky, kx) * dphim(:, ky, kx)
                enddo
            enddo
#endif

            call fftxys2p(rs, pres)

            call stop_timer(pres_timer)

        end subroutine pressure

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Gets the source terms for vorticity and buoyancy in mixed-spectral space.
        ! Note, vel obtained by vor2vel before calling this
        ! routine is spectrally truncated.
        subroutine source
#ifdef ENABLE_BUOYANCY
            !------------------------------------
            !Buoyancy source:
            call buoyancy_tendency
#endif
            !------------------------------------
            !Vorticity source:
            call vorticity_tendency

        end subroutine source

end module inversion_mod
