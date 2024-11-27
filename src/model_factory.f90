module model_factory
    use field_layout, only : layout_t
    use cheby_layout, only : cheby_layout_t
    use mss_layout, only : mss_layout_t
    use mss_ops, only : mss_ops_t
    use cheby_ops, only : cheby_ops_t
    use field_ops, only : ops_t
    use netcdf_reader
    implicit none

    private

    class(layout_t),  allocatable :: layout
    class(ops_t),     allocatable :: ops

    public :: create_model, layout, ops

contains

    subroutine create_model(fname)
        character(len=*), intent(in) :: fname
        integer                      :: ncid
        character(len=16)            :: grid_type

        call open_netcdf_file(fname, NF90_NOWRITE, ncid)

        call read_netcdf_attribute(ncid, 'grid_type', grid_type)

        call close_netcdf_file(ncid)

        select case(grid_type)
            case('Chebyshev')
                allocate(cheby_layout_t :: layout)
                allocate(cheby_ops_t :: ops)
            case('Uniform')
                allocate(mss_layout_t :: layout)
                allocate(mss_ops_t :: ops)
            case default
        end select

    end subroutine create_model

end module model_factory
