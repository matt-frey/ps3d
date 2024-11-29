module field_diagnostics
    use model, only : layout
    use parameters, only : nz, ncelli, dx, lower, extent   &
                         , fnzi, ncell
    use constants, only : zero, f12, f14, one, small
    use merge_sort
    use sta3dfft, only : fftxys2p, diffx, diffy, ztrig, zfactors
    use stafft, only : dst
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
            integer                      :: i, j
            double precision             :: z(0:nz)

            z = layout%get_z_axis()

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
        ! Note: If called with l_global = .false., the result must be normalised with "ncell"
        !       after a global reduction by the user.
        function get_kinetic_energy(vv, l_global, l_allreduce) result(ke)
            double precision, intent(in) :: vv(box%lo(3):box%hi(3),      &
                                               box%lo(2):box%hi(2),      &
                                               box%lo(1):box%hi(1), 3)
            logical,          intent(in) :: l_global
            logical,          intent(in) :: l_allreduce
            double precision             :: fke(box%lo(3):box%hi(3),      &
                                                box%lo(2):box%hi(2),      &
                                                box%lo(1):box%hi(1))
            double precision             :: ke

            fke = vv(:, :, :, 1) ** 2   &
                + vv(:, :, :, 2) ** 2   &
                + vv(:, :, :, 3) ** 2


            if (l_global) then
                ! already includes normalisation with "ncell"
                ke = f12 * layout%get_mean(fke, l_allreduce)
            else
                ke = f12 * layout%get_local_sum(fke)
            endif

        end function get_kinetic_energy


        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! domain-averaged horizontal kinetic energy, i.e. (u^2 + v^2) / 2
        ! Note: If called with l_global = .false., the result must be normalised with "ncell"
        !       after a global reduction by the user.
        function get_horizontal_kinetic_energy(vv, l_global) result(ke)
            double precision, intent(in) :: vv(box%lo(3):box%hi(3),      &
                                               box%lo(2):box%hi(2),      &
                                               box%lo(1):box%hi(1), 3)
            logical,          intent(in) :: l_global
            double precision             :: fke(box%lo(3):box%hi(3),      &
                                                box%lo(2):box%hi(2),      &
                                                box%lo(1):box%hi(1))
            double precision             :: ke

            fke = vv(:, :, :, 1) ** 2   &
                + vv(:, :, :, 2) ** 2

            if (l_global) then
                ! already includes normalisation with "ncell"
                ke = f12 * layout%get_mean(fke, l_allreduce=.false.)
            else
                ke = f12 * layout%get_local_sum(fke)
            endif

        end function get_horizontal_kinetic_energy

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! domain-averaged vertical kinetic energy, i.e. w^2 / 2
        ! Note: If called with l_global = .false., the result must be normalised with "ncell"
        !       after a global reduction by the user.
        function get_vertical_kinetic_energy(l_global) result(ke)
            logical,         intent(in) :: l_global
            double precision            :: ke
            double precision            :: fke(box%lo(3):box%hi(3),      &
                                               box%lo(2):box%hi(2),      &
                                               box%lo(1):box%hi(1))

            fke = f12 * vel(:, :, :, 3) ** 2

            if (l_global) then
                ! already includes normalisation with "ncell"
                ke = layout%get_mean(fke, l_allreduce=.false.)
            else
                ke = layout%get_local_sum(fke)
            endif

        end function get_vertical_kinetic_energy

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! domain-averaged enstrophy
        ! Note: If called with l_global = .false., the result must be normalised with "ncell"
        !       after a global reduction by the user.
        function get_enstrophy(l_global, l_allreduce) result(en)
            logical,         intent(in) :: l_global
            logical,         intent(in) :: l_allreduce
            double precision            :: fen(box%lo(3):box%hi(3),      &
                                               box%lo(2):box%hi(2),      &
                                               box%lo(1):box%hi(1))
            double precision            :: en

            fen = vor(:, :, :, 1) ** 2  &
                + vor(:, :, :, 2) ** 2  &
                + vor(:, :, :, 3) ** 2

            if (l_global) then
                ! already includes normalisation with "ncell"
                en = f12 * layout%get_mean(fen, l_allreduce)
            else
                en = f12 * layout%get_local_sum(fen)
            endif

        end function get_enstrophy

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! domain-averaged enstrophy, i.e. (xi^2 + eta^2) / 2
        ! Note: If called with l_global = .false., the result must be normalised with "ncell"
        !       after a global reduction by the user.
        function get_horizontal_enstrophy(l_global) result(en)
            logical,         intent(in) :: l_global
            double precision            :: fen(box%lo(3):box%hi(3),      &
                                               box%lo(2):box%hi(2),      &
                                               box%lo(1):box%hi(1))
            double precision            :: en

            fen = vor(:, :, :, 1) ** 2  &
                + vor(:, :, :, 2) ** 2

            if (l_global) then
                ! already includes normalisation with "ncell"
                en = f12 * layout%get_mean(fen, l_allreduce=.false.)
            else
                en = f12 * layout%get_local_sum(fen)
            endif

        end function get_horizontal_enstrophy

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! maximum horizontal enstrophy, i.e. max(sqrt(xi^2 + eta^2))
        function get_max_horizontal_enstrophy(l_global) result(hemax)
            logical, intent(in) :: l_global
            double precision    :: hemax

            hemax = maxval(sqrt(vor(:, :, :, 1) ** 2 + vor(:, :, :, 2) ** 2))

            if (l_global) then
                call mpi_blocking_reduce(hemax, MPI_MAX, world)
            endif

        end function get_max_horizontal_enstrophy

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! domain-averaged enstrophy, i.e. zeta^2 / 2
        ! Note: If called with l_global = .false., the result must be normalised with "ncell"
        !       after a global reduction by the user.
        function get_vertical_enstrophy(l_global) result(en)
            logical,         intent(in) :: l_global
            double precision            :: fen(box%lo(3):box%hi(3),      &
                                               box%lo(2):box%hi(2),      &
                                               box%lo(1):box%hi(1))
            double precision             :: en

            fen = f12 * vor(:, :, :, 3) ** 2

            if (l_global) then
                ! already includes normalisation with "ncell"
                en = layout%get_mean(fen, l_allreduce=.false.)
            else
                en = layout%get_local_sum(fen)
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

            call layout%diffz(buoy, dbdz)

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

            if (f_cor(3) == zero) then
                ro = huge(zero)
                return
            endif

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

            if (f_cor(3) == zero) then
                ro = huge(zero)
                return
            endif

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

            call layout%diffz(buoy, dbdz)

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
            call flayout%combine_semi_spectral(sbuoy)
            call diffx(sbuoy, ds)
            call fftxys2p(ds, dbdx)

            call diffy(sbuoy, ds)
            call fftxys2p(ds, dbdy)

            call layout%diffz(sbuoy, mag)
            call fftxys2p(ds, mag)
            call flayout%decompose_semi_spectral(sbuoy)

            ! mag = |gradb|
            mag = sqrt(dbdx ** 2 + dbdy ** 2 + mag ** 2)

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
                vormean(nc) =  layout%get_local_mean(vor(:, :, :, nc))
            enddo

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

end module field_diagnostics
