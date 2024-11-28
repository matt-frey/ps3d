module model_factory
    use field_layout, only : layout_t
    use cheby_layout, only : cheby_layout_t
    use mss_layout, only : mss_layout_t
    use mss_ops, only : mss_ops_t
    use cheby_ops, only : cheby_ops_t
    use field_ops, only : ops_t
    use mpi_utils, only : mpi_stop
    implicit none

    private

    class(layout_t),  allocatable :: layout
    class(ops_t),     allocatable :: ops

    public :: create_model, layout, ops

contains

    subroutine create_model(grid_type)
        character(len=*), intent(in) :: grid_type

        select case(grid_type)
            case('chebyshev')
                allocate(cheby_layout_t :: layout)
                allocate(cheby_ops_t :: ops)
            case('uniform')
                allocate(mss_layout_t :: layout)
                allocate(mss_ops_t :: ops)
            case default
                call mpi_stop(&
                    "Error in model creation. No grid type '" // grid_type // "'.")
        end select

    end subroutine create_model

end module model_factory
