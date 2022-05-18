! =============================================================================
!     This module specifies all fields and implements specific subroutines
!     and functions.
! =============================================================================
module fields
    use parameters, only : nx, ny, nz
    use constants, only : zero
    implicit none

    ! x: zonal
    ! y: meridional
    ! z: vertical
    ! Due to periodicity in x and y, the grid points in x go from 0 to nx-1
    ! and from 0 to ny-1 in y
    double precision, allocatable, dimension(:, :, :, :) :: &
        svortg,    &   ! vorticity vector field in spectral space
        vortg,     &   ! vorticity vector field (\omegax, \omegay, \omegaz) in physical space
        velog,     &   ! velocity vector field (u, v, w)
        vtend          ! vorticity tendency

    double precision, allocatable, dimension(:, :, :) :: &
        buoyg,     &   ! buoyancy (physical)
        sbuoyg,    &   ! buoyancy (spectral)
        diss           ! dissipation operator (spectral)

    contains

        ! Allocate all fields
        subroutine field_alloc
            if (allocated(velog)) then
                return
            endif

            allocate(velog(0:nz, 0:ny-1, 0:nx-1, 3))

            allocate(vortg(0:nz, 0:ny-1, 0:nx-1, 3))
            allocate(svortg(0:nz, 0:nx-1, 0:ny-1, 3))

            allocate(vtend(0:nz, 0:ny-1, 0:nx-1, 3))

            allocate(buoyg(0:nz, 0:ny-1, 0:nx-1))
            allocate(sbuoyg(0:nz, 0:nx-1, 0:ny-1))

            allocate(diss(0:nz, 0:nx-1, 0:ny-1))

        end subroutine field_alloc

        ! Reset fields to zero
        subroutine field_default
            call field_alloc

            velog    = zero
            vortg    = zero
            svortg   = zero

            vtend    = zero
            buoyg    = zero
            sbuoyg   = zero
            diss     = zero
        end subroutine field_default

end module fields
