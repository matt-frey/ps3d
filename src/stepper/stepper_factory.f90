module stepper_factory
    use stepper_mod, only : stepper_t
    use cn2_mod, only : cn2
    use impl_rk4_mod, only : impl_rk4
    use mpi_utils, only : mpi_print, mpi_stop
    implicit none

contains

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function create_stepper(name) result(stepper)
        character(len=9),  intent(in) :: name
        class(stepper_t), allocatable :: stepper

        ! 27 March 2024
        ! https://stackoverflow.com/a/72958237
        select case (name)
            case ('cn2')
                call mpi_print('Using Crank-Nicholson 2nd order stepper.')
                stepper = cn2()
            case ('impl-diff-rk4')
                call mpi_print('Using implicit diffusion Runge-Kutta 4th order stepper.')
                stepper = impl_rk4()
            case default
                call mpi_stop("No stepper called '" // name // "' available.")
        end select

        call stepper%setup

    end function create_stepper

end module stepper_factory
