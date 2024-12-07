! =============================================================================
!     This module specifies all fields and implements specific subroutines
!     and functions.
! =============================================================================
module fields
    use parameters, only : nz
    use constants, only : zero
    use mpi_layout, only : box, l_mpi_layout_initialised
    use mpi_utils, only : mpi_exit_on_error
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
        svorts, &   ! vorticity source in mixed spectral space
        vortsm      ! used for time stepping

#ifdef ENABLE_BUOYANCY
    double precision, allocatable, dimension(:, :, :) :: &
        buoy,   &   ! buoyancy (physical)
        sbuoy,  &   ! full-spectral buoyancy for 1:nz-1, semi-spectral for iz = 0 and iz = nz
        sbuoys, &   ! buoyancy source in mixed spectral space
        bsm         ! used for time stepping
#endif

    double precision, allocatable, dimension(:, :) :: &
        vdiss,  &   ! dissipation operator
        bdiss

    ! initial \xi and \eta mean
    double precision :: ini_vor_mean(2)

    double precision, allocatable :: bbarz(:) ! N**2 * z

contains

    ! Allocate all fields
    subroutine field_alloc
        integer :: lo(3), hi(3)

        if (.not. l_mpi_layout_initialised) then
            call mpi_exit_on_error
        endif

        if (allocated(vel)) then
            return
        endif

        lo = box%lo
        hi = box%hi

        allocate(vel(0:nz,  lo(2):hi(2), lo(1):hi(1), 3))
        allocate(svel(0:nz, lo(2):hi(2), lo(1):hi(1), 3))

        allocate(vor(0:nz,  lo(2):hi(2), lo(1):hi(1), 3))
        allocate(svor(0:nz, lo(2):hi(2), lo(1):hi(1), 3))

        allocate(svorts(0:nz, lo(2):hi(2), lo(1):hi(1), 3))

#ifdef ENABLE_BUOYANCY
        allocate(buoy(0:nz,   lo(2):hi(2), lo(1):hi(1)))
        allocate(sbuoy(0:nz,  lo(2):hi(2), lo(1):hi(1)))
        allocate(sbuoys(0:nz, lo(2):hi(2), lo(1):hi(1)))

        allocate(bbarz(0:nz))
#endif

        allocate(vdiss(lo(2):hi(2), lo(1):hi(1)))
        allocate(bdiss(lo(2):hi(2), lo(1):hi(1)))

        ! Spectral fields needed in time stepping:
        allocate(vortsm(0:nz, lo(2):hi(2), lo(1):hi(1), 3))
#ifdef ENABLE_BUOYANCY
        allocate(bsm(0:nz, lo(2):hi(2), lo(1):hi(1)))
#endif

    end subroutine field_alloc

    ! Reset fields to zero
    subroutine field_default
        call field_alloc

        svor   = zero
        vor    = zero
        vel    = zero
        svel   = zero
        svorts = zero
        vortsm = zero
#ifdef ENABLE_BUOYANCY
        buoy   = zero
        sbuoy  = zero
        sbuoys = zero
        bsm    = zero
#endif
        vdiss  = zero
        bdiss  = zero

        ini_vor_mean = zero
    end subroutine field_default

end module fields
