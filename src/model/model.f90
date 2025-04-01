module model
    use options, only : verbose
    use field_layout, only : layout_t
    use cheby_layout, only : cheby_layout_t
    use mss_layout, only : mss_layout_t
    use mpi_utils, only : mpi_stop
    use mpi_utils, only : mpi_print
    implicit none

    private

    class(layout_t),  allocatable :: layout

    public :: create_model, layout

contains

    subroutine create_model(grid_type, filter_type)
        character(len=*), intent(in) :: grid_type
        character(len=*), intent(in) :: filter_type

        if (allocated(layout)) then
            call layout%finalise
            deallocate(layout)
        endif

        select case(grid_type)
            case('chebyshev')
                allocate(cheby_layout_t :: layout)
            case('uniform')
                allocate(mss_layout_t :: layout)
            case default
                call mpi_stop(&
                    "Error in model creation. No grid type '" // grid_type // "'.")
        end select

#ifdef ENABLE_VERBOSE
        if (verbose) then
            call mpi_print("Using " // grid_type // " model.")
        endif
#endif

        call layout%initialise(filter_type)

    end subroutine create_model

end module model
