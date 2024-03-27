module impl_rk4_mod
    use advance_mod, only : base_stepper
    use constants, only : f12, f13, f16
    use parameters, only : nz
    use fields
    use inversion_utils
    use inversion_mod, only : vor2vel, source
    use field_diagnostics
    implicit none

    double precision :: dt2, dt3, dt6

    type, extends(base_stepper) :: impl_rk4
        ! epq = exp( D * (t-t0))
        ! emq = exp(-D * (t-t0))
        double precision, allocatable :: epq(:, :), emq(:, :)
        double precision, allocatable :: fsvor(:, :, :, :)

        contains
            procedure :: setup  => impl_rk4_setup
            procedure :: step => impl_rk4_step
    end type


    contains

        subroutine impl_rk4_setup(self)
            class(impl_rk4), intent(inout) :: self

            allocate(self%epq(box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
            allocate(self%emq(box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
            allocate(self%fsvor(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1), 3))
        end subroutine impl_rk4_setup

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine impl_rk4_step(self, t, dt)
            class(impl_rk4),  intent(inout) :: self
            double precision, intent(inout) :: t
            double precision, intent(in)    :: dt
            integer                         :: nc
            integer                         :: iz

            dt2 = f12 * dt
            dt3 = f13 * dt
            dt6 = f16 * dt

            ! set integrating factors
            self%emq = dexp(- dt2 * diss)
            self%epq = 1.0d0 / self%emq

            vortsm = svor
            svor = (vortsm + dt2 * svorts)
            do nc = 1, 3
                call field_combine_semi_spectral(svor(:, :, :, nc))
                !$omp parallel do private(iz)  default(shared)
                do iz = 0, nz
                    svor(iz, :, :, nc) = svor(iz, :, :, nc) * self%emq
                enddo
                !$omp end parallel do
                call field_decompose_semi_spectral(svor(:, :, :, nc))
            enddo


            self%fsvor = vortsm + dt6 * svorts

            !------------------------------------------------------------------
            !RK4 corrector step at time t0 + dt/2:
            t = t + dt2

            call vor2vel

            call source

            ! apply integrating factors to source
            do nc = 1, 3
                call field_combine_semi_spectral(svorts(:, :, :, nc))
                !$omp parallel do private(iz)  default(shared)
                do iz = 0, nz
                    svorts(iz, :, :, nc) = svorts(iz, :, :, nc) * self%epq
                enddo
                !$omp end parallel do
                call field_decompose_semi_spectral(svorts(:, :, :, nc))
            enddo

            svor = (vortsm + dt2 * svorts)
            do nc = 1, 3
                call field_combine_semi_spectral(svor(:, :, :, nc))
                !$omp parallel do private(iz)  default(shared)
                do iz = 0, nz
                    svor(iz, :, :, nc) = svor(iz, :, :, nc) * self%emq
                enddo
                !$omp end parallel do
                call field_decompose_semi_spectral(svor(:, :, :, nc))
            enddo

            self%fsvor = self%fsvor + dt3 * svorts

            !------------------------------------------------------------------
            !RK4 predictor step at time t0 + dt:

            !Invert PV and compute velocity:
            call vor2vel

            call source

            ! apply integrating factors to source
            do nc = 1, 3
                call field_combine_semi_spectral(svorts(:, :, :, nc))
                !$omp parallel do private(iz)  default(shared)
                do iz = 0, nz
                    svorts(iz, :, :, nc) = svorts(iz, :, :, nc) * self%epq
                enddo
                !$omp end parallel do
                call field_decompose_semi_spectral(svorts(:, :, :, nc))
            enddo

            self%emq = self%emq**2
            svor = vortsm + dt * svorts
            do nc = 1, 3
                call field_combine_semi_spectral(svor(:, :, :, nc))
                !$omp parallel do private(iz)  default(shared)
                do iz = 0, nz
                    svor(iz, :, :, nc) = svor(iz, :, :, nc) * self%emq
                enddo
                !$omp end parallel do
                call field_decompose_semi_spectral(svor(:, :, :, nc))
            enddo

            self%fsvor = self%fsvor + dt3 * svorts

            !------------------------------------------------------------------
            !RK4 corrector step at time t0 + dt:
            t = t + dt2

            call vor2vel
            call source

            self%epq = self%epq**2
            ! apply integrating factors to source
            do nc = 1, 3
                call field_combine_semi_spectral(svorts(:, :, :, nc))
                !$omp parallel do private(iz)  default(shared)
                do iz = 0, nz
                    svorts(iz, :, :, nc) = svorts(iz, :, :, nc) * self%epq
                enddo
                !$omp end parallel do
                call field_decompose_semi_spectral(svorts(:, :, :, nc))
            enddo

            svor = self%fsvor + dt6 * svorts
            do nc = 1, 3
                call field_combine_semi_spectral(svor(:, :, :, nc))
                !$omp parallel do private(iz)  default(shared)
                do iz = 0, nz
                    svor(iz, :, :, nc) = svor(iz, :, :, nc) * self%emq
                enddo
                !$omp end parallel do
                call field_decompose_semi_spectral(svor(:, :, :, nc))
            enddo

        end subroutine impl_rk4_step

end module impl_rk4_mod
