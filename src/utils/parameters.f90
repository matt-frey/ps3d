! =============================================================================
! This module contains global parameters that stay constant throughout a
! simulation.
! =============================================================================
module parameters
    use constants
    use mpi_utils, only : mpi_stop
    implicit none

    ! mesh spacing
    double precision :: dx(3)

    ! inverse mesh spacing
    double precision :: dxi(3)

    ! grid type: 'uniform', 'chebyshev'
    character(len=9) :: gridtype

    ! grid cell volume, really area in 2D:
    double precision :: vcell

    ! inverse grid cell volume
    double precision :: vcelli

    ! number of grid cells in each dimension
    integer :: nx, ny, nz

    ! total number of grid cells
    integer :: ncell

    ! inverse of total number of grid cells
    double precision :: ncelli

    ! total number of grid points
    integer :: ngrid

    ! inverse of total number of grid points
    double precision :: ngridi

    ! domain size
    double precision :: extent(3)

    double precision :: vdomaini

    ! domain centre
    double precision :: center(3)

    ! domain half widths values
    double precision :: hl(3)

    double precision :: hli(3)

    ! domain origin
    double precision :: lower(3)

    ! domain upper boundary
    double precision :: upper(3)

    double precision :: fnzi

contains

    ! Update all parameters according to the
    ! user-defined global options.
    subroutine update_parameters

        upper = lower + extent

        dx = extent / dble((/nx, ny, nz/))
        dxi = one / dx;

        vdomaini = one / product(extent)

        vcell = product(dx)
        vcelli = one / vcell

        ncell = nx * ny * nz
        ncelli = one / dble(ncell)

        ! due to x periodicity it is only nx
        ngrid = nx * ny * (nz + 1)
        ngridi = one / dble(ngrid)

        ! domain
        center = f12 * (lower + upper)
        hl = extent / two
        hli = one / hl

        fnzi = one / dble(nz)

    end subroutine update_parameters

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_x_axis()
        double precision :: get_x_axis(0:nx-1)
        integer          :: i

        do i = 0, nx-1
            get_x_axis(i) = lower(1) + dble(i) * dx(1)
        enddo

    end function get_x_axis

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_y_axis()
        double precision :: get_y_axis(0:ny-1)
        integer          :: i

        do i = 0, nz-1
            get_y_axis(i) = lower(2) + dble(i) * dx(2)
        enddo

    end function get_y_axis

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_z_axis()
        double precision :: get_z_axis(0:nz)
        integer          :: i

        select case (gridtype)
            case ('regular')
                do i = 0, nz
                    get_z_axis(i) = lower(3) + dble(i) * dx(3)
                enddo
            case ('chebyshev')
                get_z_axis = 0.0d0 !FIXME
            case default
                call mpi_stop('Error: invalid grid type.')
        end select

    end function get_z_axis

end module parameters
