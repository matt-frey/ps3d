! =============================================================================
! This module contains global parameters that stay constant throughout a
! simulation.
! =============================================================================
module parameters
    use constants
    use netcdf_writer
    implicit none

    ! mesh spacing
    double precision :: dx(3)

    ! inverse mesh spacing
    double precision :: dxi(3)

    ! grid type: uniform or chebyshev
    character(len=16) :: grid_type

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

    subroutine write_netcdf_parameters(ncid)
        integer, intent(in) :: ncid
        integer             :: gid

        ncerr = nf90_inq_ncid(ncid, 'parameters', gid)
        if (ncerr /= 0) then
            ncerr = nf90_def_grp(ncid, 'parameters', gid)
            call check_netcdf_error("Failed to define or group 'parameters'.")
        endif

        ncerr = nf90_put_att(ncid=gid, varid=NF90_GLOBAL, name="ncells", values=(/nx, ny, nz/))
        call check_netcdf_error("Failed to define 'ncells' attribute.")

        ncerr = nf90_put_att(ncid=gid, varid=NF90_GLOBAL, name="extent", values=extent)
        call check_netcdf_error("Failed to define 'extent' attribute.")

        ncerr = nf90_put_att(ncid=gid, varid=NF90_GLOBAL, name="origin", values=lower)
        call check_netcdf_error("Failed to define 'origin' attribute.")

        ncerr = nf90_put_att(ncid=gid, varid=NF90_GLOBAL, name="grid_type", values=grid_type)
        call check_netcdf_error("Failed to define 'grid_type' attribute.")

    end subroutine write_netcdf_parameters

end module parameters
