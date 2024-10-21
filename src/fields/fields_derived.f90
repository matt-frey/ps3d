module fields_derived
    use constants, only : zero, two
    use parameters, only : nz
    use mpi_layout, only : box, l_mpi_layout_initialised
    use mpi_utils, only : mpi_exit_on_error
    use mpi_timer, only : start_timer, stop_timer
    use fields, only : svel, vor
    use inversion_utils
#ifdef ENABLE_BUOYANCY
    use fields, only : sbuoy, buoy
#endif
    use sta2dfft, only : dct
    use sta3dfft, only : ztrig, zfactors, diffx, diffy, fftxyp2s, fftxys2p
    use physics, only : f_cor
    implicit none

    integer :: pres_timer       &
             , delta_timer

    double precision, allocatable, dimension(:, :, :) :: &
        pres,   &   ! pressure field (physical space)
        delta       ! horizontal divergence (physical space)

    contains

        ! Allocate all derived fields
        subroutine field_derived_alloc
            integer :: lo(3), hi(3)

            if (.not. l_mpi_layout_initialised) then
                call mpi_exit_on_error
            endif

            if (allocated(pres)) then
                return
            endif

            lo = box%lo
            hi = box%hi

            allocate(pres(0:nz, lo(2):hi(2), lo(1):hi(1)))
            allocate(delta(0:nz, lo(2):hi(2), lo(1):hi(1)))

        end subroutine field_derived_alloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Reset derived fields to zero
        subroutine field_derived_default

            call field_derived_alloc

            pres  = zero
            delta = zero

        end subroutine field_derived_default

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

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine horizontal_divergence
            double precision :: uxs(0:nz,                &
                                    box%lo(2):box%hi(2), &
                                    box%lo(1):box%hi(1))
            double precision :: vys(0:nz,                &
                                    box%lo(2):box%hi(2), &
                                    box%lo(1):box%hi(1))

            call start_timer(delta_timer)

            call diffx(svel(:, :, :, 1), uxs)    ! du/dx (in spectral space)
            call diffy(svel(:, :, :, 2), vys)    ! dv/dy (in spectral space)

            !$omp parallel workshare
            uxs = uxs + vys
            !$omp end parallel workshare

            call fftxys2p(uxs, delta)

            call stop_timer(delta_timer)

        end subroutine horizontal_divergence

end module fields_derived
