module model_factory
    use cn2_mod, only : cn2
    use impl_rk4_mod, only : impl_rk4
    use field_layout, only : flayout_t
    use field_cheby, only : field_cheby_t
    use field_mss, only : field_mss_t
    use advance_mod, only : stepper_t
    use mpi_utils, only : mpi_print, mpi_stop
    implicit none

    private

    type :: model_info_t
        character(len=16) :: grid
        character(len=16) :: stepper
    end type

    type :: model_t
        class(flayout_t),  allocatable :: layout
        class(stepper_t),  allocatable :: stepper
!         class(ops_t),      allocatable :: ops
    end type

    public :: model_info_t, model_t, create_model

contains

    function create_model(model_info) result(model)
        type(model_info_t), intent(in) :: model_info
        type(model_t)                  :: model

        select case(model_info%grid)
            case('Chebyshev')
                allocate(field_cheby_t :: model%layout)
!                 allocate(ops_cheby_t :: model%ops)
            case('Uniform')
                allocate(field_mss_t :: model%layout)
!                 allocate(ops_mss_t :: model%ops)
            case default
        end select

        model%stepper = create_stepper(model_info%stepper)

    end function create_model

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function create_stepper(stepper_type) result(stepper)
        character(len=9),  intent(in) :: stepper_type
        class(stepper_t), allocatable :: stepper

        ! 27 March 2024
        ! https://stackoverflow.com/a/72958237
        select case (stepper_type)
            case ('cn2')
                call mpi_print('Using Crank-Nicholson 2nd order stepper.')
                stepper = cn2()
            case ('impl-diff-rk4')
                call mpi_print('Using implicit diffusion Runge-Kutta 4th order stepper.')
                stepper = impl_rk4()
            case default
                call mpi_stop("No stepper called '" // stepper_type // "' available.")
        end select

        call stepper%setup

    end function create_stepper

end module model_factory
