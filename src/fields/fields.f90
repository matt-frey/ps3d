! =============================================================================
!     This module specifies all fields and implements specific subroutines
!     and functions.
! =============================================================================
module fields
    use parameters, only : nx, ny, nz, ncelli, dx, lower, extent
    use constants, only : zero, f12, f14, one
    use merge_sort
    use inversion_utils, only : fftxys2p, diffx, diffy, central_diffz   &
                              , field_combine_semi_spectral             &
                              , field_decompose_semi_spectral
    use ape_density, only : ape_den
    implicit none

    ! x: zonal
    ! y: meridional
    ! z: vertical
    ! Due to periodicity in x and y, the grid points in x go from 0 to nx-1
    ! and from 0 to ny-1 in y
    double precision, allocatable, dimension(:, :, :, :) :: &
        svor,   &   ! full-spectral vorticity for 1:nz-1, semi-spectral for iz = 0 and iz = nz
        vor,    &   ! vorticity vector field (\omegax, \omegay, \omegaz) in physical space
        vel,    &   ! velocity vector field (u, v, w)
        svel,   &   ! velocity vector field (u, v, w) (semi-spectral)
        svorts      ! vorticity source in mixed spectral space

    double precision, allocatable, dimension(:, :, :) :: &
        pres        ! pressure field (physical space)

#ifdef ENABLE_BUOYANCY
    double precision, allocatable, dimension(:, :, :) :: &
        buoy,   &   ! buoyancy (physical)
        sbuoy,  &   ! full-spectral buoyancy for 1:nz-1, semi-spectral for iz = 0 and iz = nz
        sbuoys      ! buoyancy source in mixed spectral space
#endif

    double precision, allocatable, dimension(:, :) :: &
        diss        ! dissipation operator

    ! initial \xi and \eta mean
    double precision :: ini_vor_mean(2)

#if defined(ENABLE_BUOYANCY) && defined(ENABLE_PERTURBATION_MODE)
    ! buoyancy frequency squared
    double precision :: bfsq        ! N**2
    double precision, allocatable :: bbarz(:) ! N**2 * z
#endif

    contains

        ! Allocate all fields
        subroutine field_alloc
            if (allocated(vel)) then
                return
            endif

            allocate(vel(0:nz, 0:ny-1, 0:nx-1, 3))
            allocate(svel(0:nz, 0:nx-1, 0:ny-1, 3))

            allocate(vor(0:nz, 0:ny-1, 0:nx-1, 3))
            allocate(svor(0:nz, 0:nx-1, 0:ny-1, 3))

            allocate(svorts(0:nz, 0:nx-1, 0:ny-1, 3))

#ifdef ENABLE_BUOYANCY
            allocate(buoy(0:nz, 0:ny-1, 0:nx-1))
            allocate(sbuoy(0:nz, 0:nx-1, 0:ny-1))
            allocate(sbuoys(0:nz, 0:nx-1, 0:ny-1))

#ifdef ENABLE_PERTURBATION_MODE
            allocate(bbarz(0:nz))
#endif
#endif

            allocate(pres(0:nz, 0:ny-1, 0:nx-1))
            allocate(diss(0:nx-1, 0:ny-1))



        end subroutine field_alloc

        ! Reset fields to zero
        subroutine field_default
            call field_alloc

            svor   = zero
            vor    = zero
            vel    = zero
            svel   = zero
            svorts = zero
#ifdef ENABLE_BUOYANCY
            buoy   = zero
            sbuoy  = zero
            sbuoys = zero
#ifdef ENABLE_PERTURBATION_MODE
            bfsq = zero
#endif
#endif
            pres   = zero
            diss   = zero

            ini_vor_mean = zero
        end subroutine field_default

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
            do i = 0, nx-1
                do j = 0, ny-1
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

        end function get_enstrophy

#ifdef ENABLE_BUOYANCY
        ! domain-averaged grad(b) integral
        function get_gradb_integral() result(enb)
            double precision :: enb
            double precision :: ds(0:nz, 0:nx-1, 0:ny-1)
            double precision :: mag(0:nz, 0:ny-1, 0:nx-1)
            double precision :: dbdx(0:nz, 0:ny-1, 0:nx-1)
            double precision :: dbdy(0:nz, 0:ny-1, 0:nx-1)

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
        end function get_mean_vorticity

end module fields
