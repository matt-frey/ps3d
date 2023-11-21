module field_diagnostics
    use parameters, only : nz, ncelli, dx, lower, extent, fnzi, ncell
    use constants, only : zero, f12, f14, one, small
    use merge_sort
    use sta3dfft, only : fftxys2p, diffx, diffy, ztrig, zfactors
    use stafft, only : dst
    use inversion_utils, only : central_diffz                           &
                              , field_combine_semi_spectral             &
                              , field_decompose_semi_spectral
    use ape_density, only : ape_den
    use mpi_environment
    use mpi_layout, only : box
    use fields, only : vor, vel, svor, ini_vor_mean
#ifdef ENABLE_BUOYANCY
    use fields, only : buoy, sbuoy
#ifdef ENABLE_PERTURBATION_MODE
    use fields, only : bbarz
#endif
#endif
    use mpi_collectives
    use mpi_utils, only : mpi_check_for_error
    implicit none

    contains

        ! domain-averaged available potential energy
        function get_available_potential_energy() result(ape)
            double precision :: ape
#ifdef ENABLE_BUOYANCY
            integer          :: i, j, k
            double precision :: z(0:nz)

            do k = 0, nz
                z(k) = lower(3) + dble(k) * dx(3)
            enddo


            ape = zero
            do i = box%lo(1), box%hi(1)
                do j = box%lo(2), box%hi(2)
#ifdef ENABLE_PERTURBATION_MODE
                    buoy(:, j, i) = buoy(:, j, i) + bbarz
#endif
                    ape = ape + sum(ape_den(buoy(1:nz-1, j, i), z(1:nz-1))) &
                        + f12 *     ape_den(buoy(0,      j, i), z(0))       &
                        + f12 *     ape_den(buoy(nz,     j, i), z(nz))

#ifdef ENABLE_PERTURBATION_MODE
                    buoy(:, j, i) = buoy(:, j, i) - bbarz
#endif
                enddo
            enddo

            ape = ape * ncelli

            call mpi_blocking_reduce(ape, MPI_SUM, world)
#else
            ape = zero
#endif
        end function get_available_potential_energy

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! domain-averaged kinetic energy
        function get_kinetic_energy() result(ke)
            double precision :: ke

            ke = f12 * sum(vel(1:nz-1, :, :, 1) ** 2      &
                         + vel(1:nz-1, :, :, 2) ** 2      &
                         + vel(1:nz-1, :, :, 3) ** 2)     &
               + f14 * sum(vel(0,  :, :, 1) ** 2          &
                         + vel(0,  :, :, 2) ** 2          &
                         + vel(0,  :, :, 3) ** 2)         &
               + f14 * sum(vel(nz, :, :, 1) ** 2          &
                         + vel(nz, :, :, 2) ** 2          &
                         + vel(nz, :, :, 3) ** 2)

            ke = ke * ncelli

            call mpi_blocking_reduce(ke, MPI_SUM, world)

        end function get_kinetic_energy

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! domain-averaged enstrophy
        function get_enstrophy() result(en)
            double precision :: en

            en = f12 * sum(vor(1:nz-1, :, :, 1) ** 2      &
                         + vor(1:nz-1, :, :, 2) ** 2      &
                         + vor(1:nz-1, :, :, 3) ** 2)     &
               + f14 * sum(vor(0,  :, :, 1) ** 2          &
                         + vor(0,  :, :, 2) ** 2          &
                         + vor(0,  :, :, 3) ** 2)         &
               + f14 * sum(vor(nz, :, :, 1) ** 2          &
                         + vor(nz, :, :, 2) ** 2          &
                         + vor(nz, :, :, 3) ** 2)

            en = en * ncelli

            call mpi_blocking_reduce(en, MPI_SUM, world)

        end function get_enstrophy

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#ifdef ENABLE_BUOYANCY
        ! domain-averaged grad(b) integral
        function get_gradb_integral() result(enb)
            double precision :: enb
            double precision :: ds(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision :: mag(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision :: dbdx(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision :: dbdy(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))

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

            call MPI_Allreduce(MPI_IN_PLACE,            &
                               enb,                     &
                               1,                       &
                               MPI_DOUBLE_PRECISION,    &
                               MPI_SUM,                 &
                               world%comm,              &
                               world%err)

            call mpi_check_for_error(world, &
                "in MPI_Allreduce of field_diagnostics::get_gradb_integral.")

        end function get_gradb_integral
#endif

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        function get_rms(ff) result(rms)
            double precision, intent(in) :: ff(box%lo(3):box%hi(3), &
                                               box%lo(2):box%hi(2), &
                                               box%lo(1):box%hi(1))
            double precision :: rms

            rms = (f12 * sum(ff(0,      box%lo(2):box%hi(2), box%lo(1):box%hi(1)) ** 2  &
                           + ff(nz,     box%lo(2):box%hi(2), box%lo(1):box%hi(1)) ** 2) &
                       + sum(ff(1:nz-1, box%lo(2):box%hi(2), box%lo(1):box%hi(1)) ** 2)) / dble(ncell)

            call MPI_Allreduce(MPI_IN_PLACE,            &
                               rms,                     &
                               1,                       &
                               MPI_DOUBLE_PRECISION,    &
                               MPI_SUM,                 &
                               world%comm,              &
                               world%err)

            call mpi_check_for_error(world, &
                "in MPI_Allreduce of field_diagnostics::get_rms.")

            rms = dsqrt(rms)

        end function get_rms

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        function get_abs_max(ff) result(abs_max)
            double precision, intent(in) :: ff(box%lo(3):box%hi(3), &
                                               box%lo(2):box%hi(2), &
                                               box%lo(1):box%hi(1))
            double precision :: abs_max

            abs_max = maxval(dabs(ff(box%lo(3):box%hi(3),   &
                                     box%lo(2):box%hi(2),   &
                                     box%lo(1):box%hi(1))))


            call MPI_Allreduce(MPI_IN_PLACE,            &
                               abs_max,                 &
                               1,                       &
                               MPI_DOUBLE_PRECISION,    &
                               MPI_MAX,                 &
                               world%comm,              &
                               world%err)

            call mpi_check_for_error(world, &
                "in MPI_Allreduce of field_diagnostics::get_abs_max.")

        end function get_abs_max

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Characteristic vorticity,  <vor^2>/<|vor|> for |vor| > vor_rms:
        !@pre vor array updated
        function get_char_vorticity(vortrms) result(vorch)
            double precision, intent(in) :: vortrms
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

            call MPI_Allreduce(MPI_IN_PLACE,            &
                               buf(1:2),                &
                               2,                       &
                               MPI_DOUBLE_PRECISION,    &
                               MPI_SUM,                 &
                               world%comm,              &
                               world%err)

            vorl1 = buf(1)
            vorl2 = buf(2)

            vorch = vorl2 / vorl1

        end function get_char_vorticity

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        function get_mean_vorticity() result(vormean)
            double precision :: vormean(3)
            integer          :: nc

            do nc = 1, 3
                !$omp parallel workshare
                vormean(nc) =       sum(vor(1:nz-1, :, :, nc)) &
                            + f12 * sum(vor(0,      :, :, nc)) &
                            + f12 * sum(vor(nz,     :, :, nc))
                !$omp end parallel workshare
            enddo

            vormean = vormean * ncelli

            call MPI_Allreduce(MPI_IN_PLACE,            &
                               vormean(1:3),            &
                               3,                       &
                               MPI_DOUBLE_PRECISION,    &
                               MPI_SUM,                 &
                               world%comm,              &
                               world%err)

            call mpi_check_for_error(world, &
                "in MPI_Allreduce of field_diagnostics::get_mean_vorticity.")

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
