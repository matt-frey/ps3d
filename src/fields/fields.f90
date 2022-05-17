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
    ! Halo grid points in vertical direction z are -1 and nz+1,
    ! hence the valid regrion is from 0 to nz
    ! Due to periodicity in x and y, the grid points in x go from 0 to nx-1
    ! and from 0 to ny-1 in y
    double precision, allocatable, dimension(:, :, :, :) :: &
        velog,     &   ! velocity vector field (u, v, w)
        vortg,     &   ! vorticity vector field (\omegax, \omegay, \omegaz)
        vtend,     &   ! vorticity tendency
        velgradg       ! velocity gradient tensor
                       ! ordering: du/dx, du/dy,
                       !                  dv/dy,
                       !           dw/dx, dw/dy
                       ! the derivatives dv/dx, du/dz, dv/dz and dw/dz
                       ! are calculated on the fly with vorticity
                       ! or the assumption of incompressibility (du/dx + dv/dy + dw/dz = 0):
                       !    dv/dx = \omegaz + du/dy
                       !    du/dz = \omegay + dw/dx
                       !    dv/dz = dw/dy - \omegax
                       !    dw/dz = - (du/dx + dv/dy)

    double precision, allocatable, dimension(:, :, :) :: &
        buoyg          ! buoyancy

    contains

        ! Allocate all fields
        subroutine field_alloc
            if (allocated(velog)) then
                return
            endif

            allocate(velog(-1:nz+1, 0:ny-1, 0:nx-1, 3))
            allocate(velgradg(-1:nz+1, 0:ny-1, 0:nx-1, 5))

            allocate(vortg(-1:nz+1, 0:ny-1, 0:nx-1, 3))

            allocate(vtend(-1:nz+1, 0:ny-1, 0:nx-1, 3))

            allocate(buoyg(-1:nz+1, 0:ny-1, 0:nx-1))

        end subroutine field_alloc

        ! Reset fields to zero
        subroutine field_default
            call field_alloc

            velog    = zero
            velgradg = zero
            vortg    = zero
            vtend    = zero
            buoyg    = zero
        end subroutine

end module fields
