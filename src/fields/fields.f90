! =============================================================================
!     This module specifies all fields and implements specific subroutines
!     and functions.
! =============================================================================
module fields
    use parameters, only : nx, ny, nz, vcell, ncell, ncelli
    use constants, only : zero, f12, f14
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

            allocate(diss(0:nx-1, 0:ny-1))

        end subroutine field_alloc

        ! Reset fields to zero
        subroutine field_default
            call field_alloc

            vel    = zero
            vor    = zero
            svor   = zero
            svorts = zero
#ifdef ENABLE_BUOYANCY
            buoy   = zero
            sbuoy  = zero
#endif
            diss   = zero

            ini_vor_mean = zero
        end subroutine field_default

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
            ke = ke * vcell * dble(ncell)
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
            en = en * vcell * dble(ncell)

        end function get_enstrophy

        function get_mean_vorticity() result(vormean)
            double precision :: vormean(3)
            integer          :: nc

            do nc = 1, 3
                vormean(nc) =       sum(vor(1:nz-1, :, :, nc)) &
                            + f12 * sum(vor(0,      :, :, nc)) &
                            + f12 * sum(vor(nz,     :, :, nc))
            enddo

            vormean = vormean * ncelli
        end function get_mean_vorticity

end module fields
