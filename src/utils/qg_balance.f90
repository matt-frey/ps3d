module qg_balance
    use constants, only : zero, one
    use parameters, only : nz, extent
    use physics, only : bfsq, f_cor
    use inversion_utils, only : k2l2
    use fields, only : svel, sbuoy
    use mpi_layout
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

    contains

        subroutine initialise_qg_balance
            integer          :: kz, ky, iz
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
                z = lower(3) + dble(iz)
                zz(iz) = bf * z / f_cor(3)
            enddo

            !------------------------------------------------------------------
            ! Fill pq and bq:
            !
            ! bq = sinh(K * (z + D)) / sinh(K * D) = (exp[K * (D+z)] - exp[-K * (D+z)]) / (exp[DK] - exp[-DK])
            !    = (exp[2DK + Kz] - exp[-Kz]) / (exp[2DK] - 1)

            ! pq = cosh(K * (z+D)) / (NK * sinh(KD))
            !    = (exp[-K*(D+z)] + exp[K*(D+z)]) / (KN * (exp[DK] - exp[-DK]))
            !    = (exp[2DK+Kz] + exp[-Kz]) / (KN * (exp¯[2DK] - 1))
            do kx = box%lo(1), box%hi(1)
                do ky = max(1, box%lo(2)), box%hi(2)
                    call set_qg_functions(kx, ky, zz)
                enddo
            enddo

            ! ky = 0
            if (box%lo(2) == 0) then
                do kx = max(1, box%lo(1)), box%hi(1)
                    call set_qg_functions(kx, 0, zz)
                enddo
            endif

            if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
                pq(:, 0, 0) = zero
                bq(:, 0, 0) = zz / depth + one
            endif

        end subroutine initialise_qg_balance

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine finalise_qg_balance

            if (allocated(pq)) then
                deallocate(pq)
                deallocate(bq)
            endif

        end subroutine finalise_qg_balance

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! for kx > 0 and ky >= 0 or kx >= 0 and ky > 0
        subroutine set_qg_functions(kx, ky, zz)
            integer,          intent(in) :: kx, ky
            double precision, intent(in) :: zz(0:nz)
            double precision             :: em(0:nz), ep(0:nz)
            double precision             :: ef, ed

            kl = dsqrt(k2l2(ky, kx))

            ed = dexp(two * depth * kl)  ! exp[2DK]
            ep = dexp(kl * zz)           ! exp[Kz]
#ifndef NDEBUG
            ! To avoid "Floating-point exception - erroneous arithmetic operation"
            ! when ep and em are really small.
            ed = max(ed, dsqrt(tiny(ed)))
            ep = max(ep, dsqrt(tiny(ep)))
#endif
            em = one / ep   ! exp[-Kz]

            ! ef = 1 / (exp[2DK] - 1)
            ef = one / (ed - one)

#ifndef NDEBUG
            ! To avoid "Floating-point exception - erroneous arithmetic operation"
            ! when ep and em are really small.
            ef = max(ef, dsqrt(tiny(ef)))
#endif

            !pq = (exp[2DK+Kz] + exp[-Kz]) / (KN * (exp[2DK] - 1))
            !bq = (exp[2DK+Kz] - exp[-Kz]) / (exp[2DK] - 1)
            pq(:, ky, kx) = ef * (ed * ep + em) / (kl * bf)
            bq(:, ky, kx) = ef * (ed * ep - em)

        end subroutine set_qg_functions

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! psi: streamfunction
        subroutine balance_fields
            double precision :: psi(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)) ! semi-spectral
            double precision :: ds(0:nz,  box%lo(2):box%hi(2), box%lo(1):box%hi(1)) ! semi-spectral
            integer          :: iz

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
            call fftxys2p(ds, ubal)

            ! u = - dpsi/dy
            ubal = -ubal

            ! v = dpsi/dx
            call diffx(psi, ds)
            call fftxys2p(ds, vbal)

        subroutine balance_fields



end module qg_balance
