module field_diagnostics
    use parameters, only : nz, ncelli, dx, lower, extent, fnzi
    use constants, only : zero, f12, f14, one
    use merge_sort
    use sta3dfft, only : fftxys2p, diffx, diffy, ztrig, zfactors
    use stafft, only : dst
    use inversion_utils, only : central_diffz                           &
                              , field_combine_semi_spectral             &
                              , field_decompose_semi_spectral
    use ape_density, only : ape_den
    use mpi_environment
    use mpi_layout, only : box
    use fields, only : buoy, vor, vel, svor, sbuoy, bbarz, ini_vor_mean
    use mpi_collectives
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

            call mpi_blocking_reduce(enb, MPI_SUM, world)

        end function get_gradb_integral
#endif

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

            call mpi_blocking_reduce(vormean, MPI_SUM, world)

        end function get_mean_vorticity

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
