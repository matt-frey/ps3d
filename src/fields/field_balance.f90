module field_balance
    use constants, only : zero, one, two
    use parameters, only : nz, extent, dx, lower
    use physics, only : bfsq, f_cor
    use inversion_utils, only : k2l2                            &
                              , field_decompose_semi_spectral   &
                              , field_combine_semi_spectral
    use sta3dfft, only : diffx, diffy, fftxys2p
    use fields, only : vel, sbuoy, buoy
    use mpi_layout
    use field_diagnostics, only : get_available_potential_energy    &
                                , get_horizontal_kinetic_energy     &
                                , get_kinetic_energy
    implicit none

    private

    double precision :: depth

    ! buoyancy frequency, N
    double precision :: bf

    ! z' = N * z / f, K = sqrt(k^2 + l^2) (wavenumber magnitude)
    !
    !   pq(z', k, l) = cosh(K * (z' + D)) / (N * K * sinh(K*D))
    !   bq(z', k, l) = sinh(K * (z' + D)) / sinh(K * D)
    !
    !   if K = 0 --> pq(z', 0, 0) = 0 and bq(z', 0, 0) = z' / D + 1
    !
    double precision, allocatable :: pq(:, :, :)    ! in semi-spectral space
    double precision, allocatable :: bq(:, :, :)    ! in semi_spectral space

    double precision :: kebal   ! balanced kinetic energy
    double precision :: keubal  ! imbalanced kinetic energy
    double precision :: apebal  ! balanced available potential energy
    double precision :: apeubal ! imbalanced available potential energy

    public :: initialise_balance        &
            , finalise_balance          &
            , balance_fields            &
            , kebal, keubal             &
            , apebal, apeubal

    contains

        subroutine initialise_balance
            integer          :: kx, ky, iz
            double precision :: zz(0:nz), z

            if (allocated(pq)) then
                return
            endif

            ! buoyancy frequency, N
            bf = dsqrt(bfsq)

            ! D = N * H / f where H is the original depth of the domain
            depth = bf * extent(3) / f_cor(3)

            allocate(pq(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
            allocate(bq(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

            !------------------------------------------------------------------
            ! Define scaled height:
            do iz = 0, nz
                z = lower(3) + dble(iz) * dx(3)
                zz(iz) = bf * z / f_cor(3)
            enddo

            !------------------------------------------------------------------
            ! Fill pq and bq:
            !
            ! bq = sinh(K * (z + D)) / sinh(K * D) = (exp[K * (D+z)] - exp[-K * (D+z)]) / (exp[DK] - exp[-DK])
            !    = (exp[Kz] - exp[-2DK - Kz]) / (1 - exp[-2DK])
            !    = exp[Kz] * (1 - exp[-2DK - 2Kz]) / (1 - exp[-2DK])
            !
            ! pq = cosh(K * (z+D)) / (NK * sinh(KD))
            !    = (exp[-K*(D+z)] + exp[K*(D+z)]) / (KN * (exp[DK] - exp[-DK]))
            !    = (exp[Kz] + exp[-2DK-Kz]) / (KN * (1 - exp[-2DK]))
            !    = exp[Kz] * (1 + exp[-2DK-2Kz]) / (KN * (1 - exp[-2DK]))
            do kx = box%lo(1), box%hi(1)
                do ky = max(1, box%lo(2)), box%hi(2)
                    call set_hyperbolic_functions(kx, ky, zz)
                enddo
            enddo

            ! ky = 0
            if (box%lo(2) == 0) then
                do kx = max(1, box%lo(1)), box%hi(1)
                    call set_hyperbolic_functions(kx, 0, zz)
                enddo
            endif

            if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
                pq(:, 0, 0) = zero
                bq(:, 0, 0) = zz / depth + one
            endif

        end subroutine initialise_balance

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine finalise_balance

            if (allocated(pq)) then
                deallocate(pq)
                deallocate(bq)
            endif

        end subroutine finalise_balance

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! for kx > 0 and ky >= 0 or kx >= 0 and ky > 0
        subroutine set_hyperbolic_functions(kx, ky, zz)
            integer,          intent(in) :: kx, ky
            double precision, intent(in) :: zz(0:nz)
            double precision             :: ep(0:nz), ed(0:nz)
            double precision             :: ef, kl

            kl = dsqrt(k2l2(ky, kx))

            ed = dexp(-two * kl * (zz + depth))  ! exp[-2DK - 2Kz]
            ep = dexp(kl * zz)                   ! exp[Kz]
#ifndef NDEBUG
            ! To avoid "Floating-point exception - erroneous arithmetic operation"
            ! when ep and ed are really small.
            ed = max(ed, dsqrt(tiny(ed)))
            ep = max(ep, dsqrt(tiny(ep)))
#endif
            ! ef = 1 / (1 - exp[-2DK])
            ef = one / (one - dexp(- two * kl * depth))

#ifndef NDEBUG
            ! To avoid "Floating-point exception - erroneous arithmetic operation"
            ! when ef is really small.
            ef = max(ef, dsqrt(tiny(ef)))
#endif

            !pq = exp[Kz] * (1 + exp[-2DK-2Kz]) / (KN * (1 - exp[-2DK]))
            pq(:, ky, kx) = ef * ep * (one + ed) / (kl * bf)

            !bq = exp[Kz] * (1 - exp[-2DK - 2Kz]) / (1 - exp[-2DK])
            bq(:, ky, kx) = ef * ep * (one - ed)

        end subroutine set_hyperbolic_functions

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! psi: streamfunction
        ! Since the horizontal divergence u_x + v_y = 0, we can use a streamfunction to get u and v.
        subroutine balance_fields(l_global)
            logical, intent(in) :: l_global
            double precision    :: psi(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)) ! semi-spectral
            double precision    :: ds(0:nz,  box%lo(2):box%hi(2), box%lo(1):box%hi(1)) ! semi-spectral
            double precision    :: velbal(0:nz,  box%lo(2):box%hi(2), box%lo(1):box%hi(1), 3)
            double precision    :: bbal(0:nz,  box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            integer             :: iz

            !------------------------------------------------------------------
            ! Define streamfunction and obtain balanced buoyancy anomaly:
            call field_combine_semi_spectral(sbuoy)
            do iz = 0, nz
                psi(iz, :, :) = pq(iz, :, :) * sbuoy(nz, :, :)
                ds(iz, :, :)  = bq(iz, :, :) * sbuoy(nz, :, :)
            enddo
            call field_decompose_semi_spectral(sbuoy)
            call fftxys2p(ds, bbal)

            !------------------------------------------------------------------
            ! Obtain balanced velocity (u, v). Note: w = 0.
            call diffy(psi, ds)
            call fftxys2p(ds, velbal(:, :, :, 1))

            ! u = - dpsi/dy
            velbal(:, :, :, 1) = -velbal(:, :, :, 1)

            ! v = dpsi/dx
            call diffx(psi, ds)
            call fftxys2p(ds, velbal(:, :, :, 2))

            !------------------------------------------------------------------
            ! Calculate balanced properties:

            ! domain-averaged balanced kinetic energy: (note: for balanced fields w = 0)
            kebal = get_horizontal_kinetic_energy(velbal, l_global)

            ! domain-averaged balanced ape:
            apebal = get_available_potential_energy(bbal, l_global, l_allreduce=.false.)

            !------------------------------------------------------------------
            ! Calculate imbalanced properties:

            velbal(:, :, :, 3) = zero
            velbal = velbal - vel

            ! domain-averaged imbalanced kinetic energy:
            keubal = get_kinetic_energy(velbal, l_global, l_allreduce=.false.)

            bbal = buoy - bbal

            ! domain-averaged imbalanced ape:
            apeubal = get_available_potential_energy(bbal, l_global, l_allreduce=.false.)

        end subroutine balance_fields

end module field_balance
