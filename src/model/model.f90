module model
    use options, only : verbose
    use field_layout, only : layout_t
    use field_filter, only : filter_t
    use cheby_layout, only : cheby_layout_t
    use mss_layout, only : mss_layout_t
    use mpi_utils, only : mpi_stop
    use cheby_filter, only : cheby_filter_t
    use mss_filter, only : mss_filter_t
    use mpi_utils, only : mpi_print
    implicit none

    private

    class(layout_t),  allocatable :: layout
    class(filter_t),  allocatable :: filter

    public :: create_model, layout, filter

contains

    subroutine create_model(grid_type, filter_type)
        character(len=*), intent(in) :: grid_type
        character(len=*), intent(in) :: filter_type

        if (allocated(layout)) then
            call layout%finalise
            deallocate(layout)
        endif

        if (allocated(filter)) then
            call filter%finalise
            deallocate(filter)
        endif

        select case(grid_type)
            case('chebyshev')
                allocate(cheby_layout_t :: layout)
                allocate(cheby_filter_t :: filter)
            case('uniform')
                allocate(mss_layout_t :: layout)
                allocate(mss_filter_t :: filter)
            case default
                call mpi_stop(&
                    "Error in model creation. No grid type '" // grid_type // "'.")
        end select

#ifdef ENABLE_VERBOSE
        if (verbose) then
            call mpi_print("Using " // grid_type // " model.")
        endif
#endif

        call layout%initialise

        call filter%initialise(filter_type)

    end subroutine create_model

end module model
