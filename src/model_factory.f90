module model_factory
    use field_layout, only : layout_t
    use cheby_layout, only : cheby_layout_t
    use mss_layout, only : mss_layout_t
    use mss_ops, only : mss_ops_t
    use cheby_ops, only : cheby_ops_t
    use field_ops, only : ops_t
    implicit none

    private


    type :: model_info_t
        character(len=16) :: grid
    end type

    class(layout_t),  allocatable :: layout
    class(ops_t),     allocatable :: ops

    public :: model_info_t, create_model, layout, ops

contains

    subroutine create_model(model_info)
        type(model_info_t), intent(in) :: model_info

        select case(model_info%grid)
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
