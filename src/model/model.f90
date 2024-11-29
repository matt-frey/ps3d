module model
    use field_layout, only : layout_t
    use field_filter, only : filter_t
    use cheby_layout, only : cheby_layout_t
    use mss_layout, only : mss_layout_t
    use mpi_utils, only : mpi_stop
    implicit none

    private

    class(layout_t),  allocatable :: layout
    class(filter_t),  allocatable :: filter

    public :: create_model, layout, filter

contains

    subroutine create_model(grid_type, filter_type)
        character(len=*), intent(in) :: grid_type
        character(len=*), intent(in) :: filter_type

        select case(grid_type)
            case('chebyshev')
                allocate(cheby_layout_t :: layout)
            case('uniform')
                allocate(mss_layout_t :: layout)
            case default
                call mpi_stop(&
                    "Error in model creation. No grid type '" // grid_type // "'.")
        end select

        call layout%initialise

        call filter%initialise(filter_type)

    end subroutine create_model

end module model
