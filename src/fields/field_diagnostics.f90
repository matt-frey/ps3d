module field_diagnostics
    use parameters, only : nz, ncelli, dx, lower, extent, fnzi, ncell, acell
    use constants, only : zero, f12, f14, one, small
    use merge_sort
    use sta3dfft, only : fftxys2p, diffx, diffy, ztrig, zfactors, fftxyp2s
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
    use inversion_utils, only : zderiv, zg, zccw
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

            ape = zero
            do i = box%lo(1), box%hi(1)
                do j = box%lo(2), box%hi(2)
                    ape = ape + sum(ape_den(bb(1:nz-1, j, i), zg(1:nz-1))) &
                        + f12 *     ape_den(bb(0,      j, i), zg(0))       &
                        + f12 *     ape_den(bb(nz,     j, i), zg(nz))
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
            integer                      :: iz

            ke = zero

            do iz = 0, nz
                ke = ke + zccw(iz) * sum(vv(iz, :, :, 1) ** 2      &
                                       + vv(iz, :, :, 2) ** 2      &
                                       + vv(iz, :, :, 3) ** 2)
            enddo

            ! The factor dx(1) * dx(2) comes from the trapezoidal rule in x and y
            ! (note that the problem is periodic in x and y, simplifying the
            ! trapezoidal rule, i.e. no 1/2 factor)
            ! The factor f12 * extent(3) comes from the mapping [-1, 1] to [a, b]
            ! where the Chebyshev points are given in [-1, 1]
            ! z = (b-a) / 2 * t + (a+b)/2 for [a, b] --> dz = (b-a) / 2 * dt
            ke = f12 * ke * f12 * extent(3) * dx(1) * dx(2)

            ! divide by domain volume to get domain-average
            ke = ke / product(extent)

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
            integer                      :: iz

            ke = zero

            do iz = 0, nz
                ke = ke + zccw(iz) * sum(vv(iz, :, :, 1) ** 2       &
                                       + vv(iz, :, :, 2) ** 2)
            enddo

            ke = f12 * ke * f12 * extent(3) * dx(1) * dx(2)

            ! divide by domain volume to get domain-average
            ke = ke / product(extent)

            if (l_global) then
                call mpi_blocking_reduce(ke, MPI_SUM, world)
            endif

        end function get_horizontal_kinetic_energy

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! domain-averaged vertical kinetic energy, i.e. w^2 / 2
        function get_vertical_kinetic_energy(l_global) result(ke)
            logical, intent(in) :: l_global
            double precision    :: ke
            integer             :: iz

            ke = zero

            do iz = 0, nz
                ke = ke + zccw(iz) * sum(vel(iz, :, :, 3) ** 2)
            enddo

            ke = f12 * ke * f12 * extent(3) * dx(1) * dx(2)

            ! divide by domain volume to get domain-average
            ke = ke / product(extent)

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
            integer             :: iz

            en = zero

            do iz = 0, nz
                en = en + zccw(iz) * sum(vor(iz, :, :, 1) ** 2      &
                                       + vor(iz, :, :, 2) ** 2      &
                                       + vor(iz, :, :, 3) ** 2)
            enddo

            ! See get_kinetic_energy for explanation of factors
            en = f12 * en * f12 * extent(3) * dx(1) * dx(2)

            ! divide by domain volume to get domain-average
            en = en / product(extent)

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
            integer             :: iz

            en = zero

            do iz = 0, nz
                en = en + zccw(iz) * sum(vor(iz, :, :, 1) ** 2      &
                                       + vor(iz, :, :, 2) ** 2)
            enddo

            ! See get_kinetic_energy for explanation of factors
            en = f12 * en * f12 * extent(3) * dx(1) * dx(2)

            ! divide by domain volume to get domain-average
            en = en / product(extent)

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
            integer             :: iz

            en = zero

            do iz = 0, nz
                en = en + zccw(iz) * sum(vor(iz, :, :, 3) ** 2)
            enddo

            ! See get_kinetic_energy for explanation of factors
            en = f12 * en * f12 * extent(3) * dx(1) * dx(2)

            ! divide by domain volume to get domain-average
            en = en / product(extent)

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
            call diffx(sbuoy, ds)
            call fftxys2p(ds, dbdx)

            call diffy(sbuoy, ds)
            call fftxys2p(ds, dbdy)

            call zderiv(sbuoy, mag)
            call fftxys2p(ds, mag)

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
            double precision             :: mean
            integer                      :: iz


            mean = zero

            do iz = 0, nz
                mean = mean + zccw(iz) * sum(ff(iz, :, :))
            enddo

            mean = mean * f12 * extent(3) * dx(1) * dx(2)

            mean = mean / product(extent)

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
            integer                      :: iz

            rms = zero
            do iz = 0, nz
                rms = rms + zccw(iz) * sum(ff(iz, :, :) ** 2)
            enddo

            rms = rms * f12 * extent(3) * dx(1) * dx(2)

            rms = rms / product(extent)

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
            integer             :: nc, iz


            vormean = zero
            do nc = 1, 3
                do iz = 0, nz
                    vormean(nc) = vormean(nc)  &
                                + zccw(iz) * sum(vor(iz, :, :, nc))
                enddo
                vormean(nc) = vormean(nc) * f12 * extent(3) * dx(1) * dx(2)
                vormean(nc) = vormean(nc) / product(extent)
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

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! This is only calculated on the MPI rank having kx = ky = 0
        function calc_vorticity_mean() result(savg)
            integer          :: nc, iz
            double precision :: savg(3)

            if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
                do nc = 1, 3
                    savg(nc) = zero
                    do iz = 0, nz
                        savg(nc) = savg(nc) + zccw(iz) * svor(iz, 0, 0, nc)
                    enddo
                    ! The factor f12 * extent(3) comes from the mapping [-1, 1] to [a, b]
                    ! where the Chebyshev points are given in [-1, 1]
                    ! z = (b-a) / 2 * t + (a+b)/2 for [a, b] --> dz = (b-a) / 2 * dt
                    ! However, we must divide by extent(3) again in order to get the vertical domain-average.
                    ! Hence, we only scale by f12.
                    savg(nc) = savg(nc) * f12
                enddo
            endif
        end function calc_vorticity_mean

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! This is only calculated on the MPI rank having kx = ky = 0
        subroutine adjust_vorticity_mean
            double precision :: savg(3)
            integer          :: nc

            savg = calc_vorticity_mean()


             if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
                ! Ensure zero global mean horizontal vorticity conservation:
                do nc = 1, 3
                    svor(: , 0, 0, nc) = svor(:, 0, 0, nc) + ini_vor_mean(nc) - savg(nc)
                enddo
             endif

        end subroutine adjust_vorticity_mean

end module field_diagnostics
