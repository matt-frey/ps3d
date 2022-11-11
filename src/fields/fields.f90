! =============================================================================
!     This module specifies all fields and implements specific subroutines
!     and functions.
! =============================================================================
module fields
    use parameters, only : nx, ny, nz, vcell, ncell, ncelli, dx, lower, extent, ngrid
    use constants, only : zero, f12, f14, one, two
    use merge_sort
    use inversion_utils, only : fftxys2p, diffx, diffy, central_diffz   &
                              , field_combine_semi_spectral             &
                              , field_decompose_semi_spectral
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

    double precision :: peref

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
#endif
            pres   = zero
            diss   = zero

            ini_vor_mean = zero
        end subroutine field_default

        subroutine calculate_peref
            integer          :: ii(ngrid), i, j, k, n, m
            double precision :: b(ngrid)
            double precision :: gam, zmean

            n = 1
            m = nz+1
            do i = 0, nx-1
                do j = 0, ny-1
                    b(n:m) = buoy(:, j, i)
                    n = m + 1
                    m = n + nz
                enddo
            enddo

            call msort(b, ii)

            gam = one / (extent(1) * extent(2))
            zmean = gam * vcell

            peref = - b(1) * vcell * zmean

            do k = 2, ngrid
                zmean = zmean + gam * two * vcell
                peref = peref - b(k) * vcell * zmean
            enddo

        end subroutine calculate_peref

        function get_potential_energy() result(pe)
            double precision :: pe
            integer          :: i, j, k
            double precision :: z(0:nz)

            do k = 0, nz
                z(k) = dble(k) * dx(3)
            enddo

            pe = zero
            do i = 0, nx-1
                do j = 0, ny-1
                    pe = pe - sum(buoy(:, j, i) * z(:))
                enddo
            enddo

            pe = pe * vcell

            pe = pe - peref

        end function get_potential_energy

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

            ! multiply with total volume
            ke = ke * vcell !* ncell

        end function get_kinetic_energy

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

            ! multiply with total volume
            en = en * vcell !* dble(ncell)

        end function get_enstrophy

#ifdef ENABLE_BUOYANCY
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
            enb = f12 * sum(mag(1:nz-1, :, :) ** 2) &
                + f14 * sum(mag(0,      :, :) ** 2) &
                + f14 * sum(mag(nz,     :, :) ** 2)

            enb = enb * vcell

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
