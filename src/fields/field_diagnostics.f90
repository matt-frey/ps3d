module field_diagnostics
    use parameters, only : nz, ncelli, dx, lower, extent, fnzi, ncell
    use constants, only : zero, f12, f14, one, small
    use merge_sort
    use sta3dfft, only : fftxys2p, diffx, diffy, ztrig, zfactors
    use stafft, only : dst
    use inversion_utils, only : central_diffz                           &
                              , field_combine_semi_spectral             &
                              , field_decompose_semi_spectral           &
                              , field_combine_physical
    use mpi_environment
    use mpi_layout, only : box
    use mpi_collectives, only : mpi_blocking_reduce
    use fields, only : vor, vel, svor, svel, ini_vor_mean
    use physics, only : f_cor
#ifdef ENABLE_BUOYANCY
    use ape_density, only : ape_den
    use fields, only : buoy, sbuoy
    use fields, only : bbarz
    use physics, only : bfsq
#endif
    use mpi_utils, only : mpi_check_for_error
    implicit none

    contains

#ifdef ENABLE_BUOYANCY
        ! domain-averaged available potential energy
        function get_available_potential_energy(bb, l_global, l_allreduce) result(ape)
            double precision, intent(inout) :: bb(box%lo(3):box%hi(3), &
                                                  box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1))
            logical,          intent(in) :: l_global
            logical,          intent(in) :: l_allreduce
            double precision             :: ape
#if !defined(NDEBUG)
            logical                      :: l_dummy
#endif
            integer                      :: i, j, k
            double precision             :: z(0:nz)

            do k = 0, nz
                z(k) = lower(3) + dble(k) * dx(3)
            enddo


            ape = zero
            do i = box%lo(1), box%hi(1)
                do j = box%lo(2), box%hi(2)
                    ape = ape + sum(ape_den(bb(1:nz-1, j, i), z(1:nz-1))) &
                        + f12 *     ape_den(bb(0,      j, i), z(0))       &
                        + f12 *     ape_den(bb(nz,     j, i), z(nz))
                enddo
            enddo

            ape = ape * ncelli

            if (l_global) then
                if (l_allreduce) then
                    call MPI_Allreduce(MPI_IN_PLACE,            &
                                       ape,                     &
                                       1,                       &
                                       MPI_DOUBLE_PRECISION,    &
                                       MPI_SUM,                 &
                                       world%comm,              &
                                       world%err)

                    call mpi_check_for_error(world, &
                        "in MPI_Allreduce of field_diagnostics::get_available_potential_energy.")
                else
                    call mpi_blocking_reduce(ape, MPI_SUM, world)
                endif
            endif

#ifndef NDEBUG
            l_dummy = l_allreduce ! just to avoid unused variable compiler error in debug mode
            l_dummy = l_global
#endif
        end function get_available_potential_energy
#endif

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! domain-averaged kinetic energy
        function get_kinetic_energy(vv, l_global, l_allreduce) result(ke)
            double precision, intent(in) :: vv(box%lo(3):box%hi(3),      &
                                               box%lo(2):box%hi(2),      &
                                               box%lo(1):box%hi(1), 3)
            logical,          intent(in) :: l_global
            logical,          intent(in) :: l_allreduce
            double precision             :: ke

            ke = f12 * sum(vv(1:nz-1, :, :, 1) ** 2      &
                         + vv(1:nz-1, :, :, 2) ** 2      &
                         + vv(1:nz-1, :, :, 3) ** 2)     &
               + f14 * sum(vv(0,  :, :, 1) ** 2          &
                         + vv(0,  :, :, 2) ** 2          &
                         + vv(0,  :, :, 3) ** 2)         &
               + f14 * sum(vv(nz, :, :, 1) ** 2          &
                         + vv(nz, :, :, 2) ** 2          &
                         + vv(nz, :, :, 3) ** 2)

            ke = ke * ncelli

            if (l_global) then
                if (l_allreduce) then
                    call MPI_Allreduce(MPI_IN_PLACE,            &
                                       ke,                      &
                                       1,                       &
                                       MPI_DOUBLE_PRECISION,    &
                                       MPI_SUM,                 &
                                       world%comm,              &
                                       world%err)

                    call mpi_check_for_error(world, &
                        "in MPI_Allreduce of field_diagnostics::get_kinetic_energy.")
                else
                    call mpi_blocking_reduce(ke, MPI_SUM, world)
                endif
            endif

        end function get_kinetic_energy


        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! domain-averaged horizontal kinetic energy, i.e. (u^2 + v^2) / 2
        function get_horizontal_kinetic_energy(vv, l_global) result(ke)
            double precision, intent(in) :: vv(box%lo(3):box%hi(3),      &
                                               box%lo(2):box%hi(2),      &
                                               box%lo(1):box%hi(1), 3)
            logical,          intent(in) :: l_global
            double precision             :: ke

            ke = f12 * sum(vv(1:nz-1, :, :, 1) ** 2    &
                         + vv(1:nz-1, :, :, 2) ** 2)   &
               + f14 * sum(vv(0,      :, :, 1) ** 2    &
                         + vv(0,      :, :, 2) ** 2)   &
               + f14 * sum(vv(nz,     :, :, 1) ** 2    &
                         + vv(nz,     :, :, 2) ** 2)

            ke = ke * ncelli

            if (l_global) then
                call mpi_blocking_reduce(ke, MPI_SUM, world)
            endif

        end function get_horizontal_kinetic_energy

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! domain-averaged vertical kinetic energy, i.e. w^2 / 2
        function get_vertical_kinetic_energy(l_global) result(ke)
            logical, intent(in) :: l_global
            double precision    :: ke

            ke = f12 * sum(vel(1:nz-1, :, :, 3) ** 2)   &
               + f14 * sum(vel(0,      :, :, 3) ** 2)   &
               + f14 * sum(vel(nz,     :, :, 3) ** 2)

            ke = ke * ncelli

            if (l_global) then
                call mpi_blocking_reduce(ke, MPI_SUM, world)
            endif

        end function get_vertical_kinetic_energy

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! domain-averaged enstrophy
        function get_enstrophy(l_global, l_allreduce) result(en)
            logical, intent(in) :: l_global
            logical, intent(in) :: l_allreduce
            double precision    :: en

            en = f12 * sum(vor(1:nz-1, :, :, 1) ** 2    &
                         + vor(1:nz-1, :, :, 2) ** 2    &
                         + vor(1:nz-1, :, :, 3) ** 2)   &
               + f14 * sum(vor(0,      :, :, 1) ** 2    &
                         + vor(0,      :, :, 2) ** 2    &
                         + vor(0,      :, :, 3) ** 2)   &
               + f14 * sum(vor(nz,     :, :, 1) ** 2    &
                         + vor(nz,     :, :, 2) ** 2    &
                         + vor(nz,     :, :, 3) ** 2)

            en = en * ncelli

            if (l_global) then
                if (l_allreduce) then
                    call MPI_Allreduce(MPI_IN_PLACE,            &
                                       en,                      &
                                       1,                       &
                                       MPI_DOUBLE_PRECISION,    &
                                       MPI_SUM,                 &
                                       world%comm,              &
                                       world%err)

                    call mpi_check_for_error(world, &
                        "in MPI_Allreduce of field_diagnostics::get_enstrophy.")
                else
                    call mpi_blocking_reduce(en, MPI_SUM, world)
                endif
            endif

        end function get_enstrophy

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! domain-averaged enstrophy, i.e. (xi^2 + eta^2) / 2
        function get_horizontal_enstrophy(l_global) result(en)
            logical, intent(in) :: l_global
            double precision    :: en

            en = f12 * sum(vor(1:nz-1, :, :, 1) ** 2    &
                         + vor(1:nz-1, :, :, 2) ** 2)   &
               + f14 * sum(vor(0,      :, :, 1) ** 2    &
                         + vor(0,      :, :, 2) ** 2)   &
               + f14 * sum(vor(nz,     :, :, 1) ** 2    &
                         + vor(nz,     :, :, 2) ** 2)

            en = en * ncelli

            if (l_global) then
                call mpi_blocking_reduce(en, MPI_SUM, world)
            endif

        end function get_horizontal_enstrophy

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! maximum horizontal enstrophy, i.e. max(sqrt(xi^2 + eta^2))
        function get_max_horizontal_enstrophy(l_global) result(hemax)
            logical, intent(in) :: l_global
            double precision    :: hemax

            hemax = maxval(dsqrt(vor(:, :, :, 1) ** 2 + vor(:, :, :, 2) ** 2))

            if (l_global) then
                call mpi_blocking_reduce(hemax, MPI_MAX, world)
            endif

        end function get_max_horizontal_enstrophy

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! domain-averaged enstrophy, i.e. zeta^2 / 2
        function get_vertical_enstrophy(l_global) result(en)
            logical, intent(in) :: l_global
            double precision    :: en

            en = f12 * sum(vor(1:nz-1, :, :, 3) ** 2)  &
               + f14 * sum(vor(0,      :, :, 3) ** 2)  &
               + f14 * sum(vor(nz,     :, :, 3) ** 2)

            en = en * ncelli

            if (l_global) then
                call mpi_blocking_reduce(en, MPI_SUM, world)
            endif

        end function get_vertical_enstrophy

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#ifdef ENABLE_BUOYANCY
        ! Ri = b_z/(xi^2+eta^2)
        function get_min_richardson_number(l_global) result(ri)
            logical, intent(in) :: l_global
            double precision    :: ri
            double precision    :: dbdz(box%lo(3):box%hi(3), &
                                        box%lo(2):box%hi(2), &
                                        box%lo(1):box%hi(1))

            call central_diffz(buoy, dbdz)

            ! As we use the pertubation mode, we only have b'_z, i.e. we must
            ! add N^2 because b_z = N^2 + b'_z
            dbdz = (bfsq + dbdz) / (vor(:, :, :, 1) ** 2 + vor(:, :, :, 2) ** 2)

            ri = minval(dbdz)

            if (l_global) then
                call mpi_blocking_reduce(ri, MPI_MIN, world)
            endif

        end function get_min_richardson_number
#endif

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Ro = zeta / f (f is the vertical Coriolis frequency)
        function get_min_rossby_number(l_global) result(ro)
            logical, intent(in) :: l_global
            double precision    :: ro

            ro = minval(vor(:, :, :, 3))
            ro = ro / f_cor(3)

            if (l_global) then
                call mpi_blocking_reduce(ro, MPI_MIN, world)
            endif

        end function get_min_rossby_number

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Ro = zeta / f (f is the vertical Coriolis frequency)
        function get_max_rossby_number(l_global) result(ro)
            logical, intent(in) :: l_global
            double precision    :: ro

            ro = maxval(vor(:, :, :, 3))
            ro = ro / f_cor(3)

            if (l_global) then
                call mpi_blocking_reduce(ro, MPI_MAX, world)
            endif

        end function get_max_rossby_number

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#ifdef ENABLE_BUOYANCY
        ! minimum static stability value, 1 + min(b'_z)/N^2 (if < 0 the flow is overturning)
        ! #pre Assumes we already have the buoyancy anomaly in physical space
        function get_minimum_static_stability(l_global) result(mss)
            logical, intent(in) :: l_global
            double precision    :: mss
            double precision    :: dbdz(box%lo(3):box%hi(3), &
                                        box%lo(2):box%hi(2), &
                                        box%lo(1):box%hi(1))

            call central_diffz(buoy, dbdz)

            mss = minval(dbdz)

            mss = one + mss / bfsq

            if (l_global) then
                call mpi_blocking_reduce(mss, MPI_MIN, world)
            endif

        end function get_minimum_static_stability
#endif

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#ifdef ENABLE_BUOYANCY
        ! domain-averaged grad(b) integral
        function get_gradb_integral(l_global, l_allreduce) result(enb)
            logical, intent(in) :: l_global
            logical, intent(in) :: l_allreduce
            double precision    :: enb
            double precision    :: ds(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision    :: mag(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision    :: dbdx(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision    :: dbdy(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))

            !------------------------------------
            !Obtain magnitude of buoyancy gradient
            call field_combine_semi_spectral(sbuoy)
            call diffx(sbuoy, ds)
            call fftxys2p(ds, dbdx)

            call diffy(sbuoy, ds)
            call fftxys2p(ds, dbdy)

            call central_diffz(sbuoy, mag)
            call fftxys2p(ds, mag)
            call field_decompose_semi_spectral(sbuoy)

            ! mag = |gradb|
            mag = dsqrt(dbdx ** 2 + dbdy ** 2 + mag ** 2)

            !------------------------------------
            ! Calculate domain integral of |gradb|
            enb =       sum(mag(1:nz-1, :, :)) &
                + f12 * sum(mag(0,      :, :)) &
                + f12 * sum(mag(nz,     :, :))

            enb = enb * ncelli

            if (l_global) then
                if (l_allreduce) then
                    call MPI_Allreduce(MPI_IN_PLACE,            &
                                       enb,                     &
                                       1,                       &
                                       MPI_DOUBLE_PRECISION,    &
                                       MPI_SUM,                 &
                                       world%comm,              &
                                       world%err)

                    call mpi_check_for_error(world, &
                        "in MPI_Allreduce of field_diagnostics::get_gradb_integral.")
                else
                    call mpi_blocking_reduce(enb, MPI_SUM, world)
                endif
            endif
        end function get_gradb_integral
#endif

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        function get_mean(ff, l_allreduce) result(mean)
            double precision, intent(in) :: ff(box%lo(3):box%hi(3), &
                                               box%lo(2):box%hi(2), &
                                               box%lo(1):box%hi(1))
            logical,          intent(in) :: l_allreduce
            double precision              :: mean

            ! (divide by ncell since lower and upper edge weights are halved)
            mean = (f12 * sum(ff(0,      box%lo(2):box%hi(2), box%lo(1):box%hi(1))  &
                            + ff(nz,     box%lo(2):box%hi(2), box%lo(1):box%hi(1))) &
                        + sum(ff(1:nz-1, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))) / dble(ncell)

            if (l_allreduce) then
                call MPI_Allreduce(MPI_IN_PLACE,            &
                                   mean,                    &
                                   1,                       &
                                   MPI_DOUBLE_PRECISION,    &
                                   MPI_SUM,                 &
                                   world%comm,              &
                                   world%err)

                call mpi_check_for_error(world, &
                    "in MPI_Allreduce of field_diagnostics::get_mean.")
            else
                call mpi_blocking_reduce(mean, MPI_SUM, world)
            endif

        end function get_mean

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        function get_rms(ff, l_allreduce) result(rms)
            double precision, intent(in) :: ff(box%lo(3):box%hi(3), &
                                               box%lo(2):box%hi(2), &
                                               box%lo(1):box%hi(1))
            logical,          intent(in) :: l_allreduce
            double precision             :: rms

            rms = (f12 * sum(ff(0,      box%lo(2):box%hi(2), box%lo(1):box%hi(1)) ** 2  &
                           + ff(nz,     box%lo(2):box%hi(2), box%lo(1):box%hi(1)) ** 2) &
                       + sum(ff(1:nz-1, box%lo(2):box%hi(2), box%lo(1):box%hi(1)) ** 2)) / dble(ncell)

            if (l_allreduce) then
                call MPI_Allreduce(MPI_IN_PLACE,            &
                                   rms,                     &
                                   1,                       &
                                   MPI_DOUBLE_PRECISION,    &
                                   MPI_SUM,                 &
                                   world%comm,              &
                                   world%err)

                call mpi_check_for_error(world, &
                    "in MPI_Allreduce of field_diagnostics::get_rms.")
            else
                call mpi_blocking_reduce(rms, MPI_SUM, world)
            endif

            rms = dsqrt(rms)

        end function get_rms

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        function get_abs_max(ff, l_allreduce) result(abs_max)
            double precision, intent(in) :: ff(box%lo(3):box%hi(3), &
                                               box%lo(2):box%hi(2), &
                                               box%lo(1):box%hi(1))
            logical,          intent(in) :: l_allreduce
            double precision             :: abs_max

            abs_max = maxval(dabs(ff(box%lo(3):box%hi(3),   &
                                     box%lo(2):box%hi(2),   &
                                     box%lo(1):box%hi(1))))


            if (l_allreduce) then
                call MPI_Allreduce(MPI_IN_PLACE,            &
                                   abs_max,                 &
                                   1,                       &
                                   MPI_DOUBLE_PRECISION,    &
                                   MPI_MAX,                 &
                                   world%comm,              &
                                   world%err)

                call mpi_check_for_error(world, &
                    "in MPI_Allreduce of field_diagnostics::get_abs_max.")
            else
                call mpi_blocking_reduce(abs_max, MPI_MAX, world)
            endif

        end function get_abs_max

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Characteristic vorticity,  <vor^2>/<|vor|> for |vor| > vor_rms:
        !@pre vor array updated
        function get_char_vorticity(vortrms, l_allreduce) result(vorch)
            double precision, intent(in) :: vortrms
            logical,          intent(in) :: l_allreduce
            double precision             :: vorch, vorl1, vorl2
            double precision             :: vortmp1, vortmp2, vortmp3
            double precision             :: buf(2)
            integer                      :: ix, iy, iz

            vorl1 = small
            vorl2 = zero
            do ix = box%lo(1), box%hi(1)
                do iy = box%lo(2), box%hi(2)
                    do iz = 1, nz
                        vortmp1 = f12 * abs(vor(iz-1, iy, ix, 1) + vor(iz, iy, ix, 1))
                        vortmp2 = f12 * abs(vor(iz-1, iy, ix, 2) + vor(iz, iy, ix, 2))
                        vortmp3 = f12 * abs(vor(iz-1, iy, ix, 3) + vor(iz, iy, ix, 3))
                        if (vortmp1 + vortmp2 + vortmp3 .gt. vortrms) then
                            vorl1 = vorl1 + vortmp1 + vortmp2 + vortmp3
                            vorl2 = vorl2 + vortmp1 ** 2 + vortmp2 ** 2 + vortmp3 ** 2
                        endif
                    enddo
                enddo
            enddo

            buf(1) = vorl1
            buf(2) = vorl2

            if (l_allreduce) then
                call MPI_Allreduce(MPI_IN_PLACE,            &
                                   buf(1:2),                &
                                   2,                       &
                                   MPI_DOUBLE_PRECISION,    &
                                   MPI_SUM,                 &
                                   world%comm,              &
                                   world%err)
            else
                call mpi_blocking_reduce(buf, MPI_SUM, world)
            endif

            vorl1 = buf(1)
            vorl2 = buf(2)

            vorch = vorl2 / vorl1

        end function get_char_vorticity

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        function get_mean_vorticity(l_allreduce) result(vormean)
            logical, intent(in) :: l_allreduce
            double precision    :: vormean(3)
            integer             :: nc

            do nc = 1, 3
                !$omp parallel workshare
                vormean(nc) =       sum(vor(1:nz-1, :, :, nc)) &
                            + f12 * sum(vor(0,      :, :, nc)) &
                            + f12 * sum(vor(nz,     :, :, nc))
                !$omp end parallel workshare
            enddo

            vormean = vormean * ncelli

            if (l_allreduce) then
                call MPI_Allreduce(MPI_IN_PLACE,            &
                                   vormean(1:3),            &
                                   3,                       &
                                   MPI_DOUBLE_PRECISION,    &
                                   MPI_SUM,                 &
                                   world%comm,              &
                                   world%err)

                call mpi_check_for_error(world, &
                    "in MPI_Allreduce of field_diagnostics::get_mean_vorticity.")
            else
                call mpi_blocking_reduce(vormean, MPI_SUM, world)
            endif

        end function get_mean_vorticity

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! This is only calculated on the MPI rank having kx = ky = 0
        function calc_vorticity_mean() result(savg)
            double precision :: wk(1:nz)
            integer          :: nc
            double precision :: savg(2)

            if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
                do nc = 1, 2
                    ! Cast svor_S = svor - svor_L onto the z grid as wk for kx = ky = 0:
                    wk(1:nz-1) = svor(1:nz-1, 0, 0, nc)
                    wk(nz) = zero
                    call dst(1, nz, wk(1:nz), ztrig, zfactors)
                    ! Compute average (first part is the part due to svor_L):
                    savg(nc) = f12 * (svor(0, 0, 0, nc) + svor(nz, 0, 0, nc)) + fnzi * sum(wk(1:nz-1))
                enddo
            endif
        end function calc_vorticity_mean

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! This is only calculated on the MPI rank having kx = ky = 0
        subroutine adjust_vorticity_mean
            double precision :: savg(2)
            integer          :: nc

            savg = calc_vorticity_mean()

            if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
                ! Ensure zero global mean horizontal vorticity conservation:
                do nc = 1, 2
                    ! Remove from boundary values (0 & nz):
                    svor(0 , 0, 0, nc) = svor(0 , 0, 0, nc) + ini_vor_mean(nc) - savg(nc)
                    svor(nz, 0, 0, nc) = svor(nz, 0, 0, nc) + ini_vor_mean(nc) - savg(nc)
                enddo
            endif

        end subroutine adjust_vorticity_mean

end module field_diagnostics
