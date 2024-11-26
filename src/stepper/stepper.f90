module stepper_mod
    implicit none

    type, abstract :: stepper_t
        contains
            procedure(base_diffusion),  deferred :: set_diffusion
            procedure(base_setup),      deferred :: setup
            procedure(base_step),       deferred :: step
    end type

    abstract interface
        subroutine base_diffusion(self, dt, vorch, bf)
            import stepper_t
            class(stepper_t), intent(inout) :: self
            double precision, intent(in)    :: dt
            double precision, intent(in)    :: vorch
            double precision, intent(in)    :: bf
        end subroutine base_diffusion

        subroutine base_setup(self)
            import stepper_t
            class(stepper_t), intent(inout) :: self
        end subroutine base_setup

        subroutine base_step(self, t, dt)
            import stepper_t
            class(stepper_t), intent(inout) :: self
            double precision, intent(inout) :: t
            double precision, intent(in)    :: dt
        end subroutine base_step

    end interface

end module stepper_mod
