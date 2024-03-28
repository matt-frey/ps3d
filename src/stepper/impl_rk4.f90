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
        double precision, allocatable :: svorf(:, :, :, :), svori(:, :, :, :)
#ifdef ENABLE_BUOYANCY
        double precision, allocatable :: sbuoyf(:, :, :), sbuoyi(:, :, :)
#endif

        contains
            procedure :: set_diffusion => impl_rk4_set_diffusion
            procedure :: setup  => impl_rk4_setup
            procedure :: step => impl_rk4_step

            procedure, private :: impl_rk4_substep_one
            procedure, private :: impl_rk4_substep_two
            procedure, private :: impl_rk4_substep_three
            procedure, private :: impl_rk4_substep_four
    end type


    contains

        subroutine impl_rk4_set_diffusion(self, df, vorch)
            class(impl_rk4),  intent(inout) :: self
            double precision, intent(in)    :: dt
            double precision, intent(in)    :: vorch
            double precision                :: dfac

            dfac = f12 * vorch * dt

            !$omp parallel workshare
            diss = dfac * hdis
            !$omp end parallel workshare

        end subroutine impl_rk4_set_diffusion

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine impl_rk4_setup(self)
            class(impl_rk4), intent(inout) :: self

            allocate(self%epq(box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
            allocate(self%emq(box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
            allocate(self%svorf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1), 3))
            allocate(self%svori(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1), 3))

#ifdef ENABLE_BUOYANCY
            allocate(self%sbuoyf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
            allocate(self%sbuoyi(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
#endif

        end subroutine impl_rk4_setup

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine impl_rk4_step(self, t, dt)
            class(impl_rk4),  intent(inout) :: self
            double precision, intent(inout) :: t
            double precision, intent(in)    :: dt
            integer                         :: nc

            dt2 = f12 * dt
            dt3 = f13 * dt
            dt6 = f16 * dt

            ! set integrating factors
            self%emq = dexp(- diss)
            self%epq = 1.0d0 / self%emq

            !------------------------------------------------------------------
            ! RK4 predictor step at time t0 + dt/2:
#ifdef ENABLE_BUOYANCY
            call self%impl_rk4_substep_one(q=sbuoy,          &
                                           sqs=sbuoys,       &
                                           qdi=self%sbuoyi,  &
                                           qdf=self%sbuoyf)
#endif

            do nc = 1, 3
                call self%impl_rk4_substep_one(q=svor(:, :, :, nc),          &
                                               sqs=svorts(:, :, :, nc),      &
                                               qdi=self%svori(:, :, :, nc),  &
                                               qdf=self%svorf(:, :, :, nc))
            enddo

            !------------------------------------------------------------------
            ! Invert and get new sources:
            call vor2vel
            call source

            !------------------------------------------------------------------
            !RK4 corrector step at time t0 + dt/2:
            t = t + dt2

#ifdef ENABLE_BUOYANCY
            call self%impl_rk4_substep_two(q=sbuoy,          &
                                           sqs=sbuoys,       &
                                           qdi=self%sbuoyi,  &
                                           qdf=self%sbuoyf)
#endif

            do nc = 1, 3
                call self%impl_rk4_substep_two(q=svor(:, :, :, nc),          &
                                               sqs=svorts(:, :, :, nc),      &
                                               qdi=self%svori(:, :, :, nc),  &
                                               qdf=self%svorf(:, :, :, nc))
            enddo

            !------------------------------------------------------------------
            ! Invert and get new sources:
            call vor2vel
            call source

            !------------------------------------------------------------------
            !RK4 predictor step at time t0 + dt:
            t = t + dt2

            self%emq = self%emq ** 2

#ifdef ENABLE_BUOYANCY
            call self%impl_rk4_substep_three(q=sbuoy,          &
                                             sqs=sbuoys,       &
                                             qdi=self%sbuoyi,  &
                                             qdf=self%sbuoyf,  &
                                             dt=dt)
#endif

            do nc = 1, 3
                call self%impl_rk4_substep_three(q=svor(:, :, :, nc),          &
                                                 sqs=svorts(:, :, :, nc),      &
                                                 qdi=self%svori(:, :, :, nc),  &
                                                 qdf=self%svorf(:, :, :, nc),  &
                                                 dt=dt)
            enddo


            !------------------------------------------------------------------
            ! Invert and get new sources:
            call vor2vel
            call source

            !------------------------------------------------------------------
            !RK4 corrector step at time t0 + dt:

            self%epq = self%epq ** 2

#ifdef ENABLE_BUOYANCY
            call self%impl_rk4_substep_four(q=sbuoy,          &
                                            sqs=sbuoys,       &
                                            qdf=self%sbuoyf)
#endif

            do nc = 1, 3
                call self%impl_rk4_substep_four(q=svor(:, :, :, nc),          &
                                                sqs=svorts(:, :, :, nc),      &
                                                qdf=self%svorf(:, :, :, nc))
            enddo

            call adjust_vorticity_mean

        end subroutine impl_rk4_step

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Initialisation step (t = t0) (predictor):
        subroutine impl_rk4_substep_one(self, q, sqs, qdi, qdf)
            class(impl_rk4),  intent(inout) :: self
            double precision, intent(inout) :: q(0:nz, box%lo(2):box%hi(2), &
                                                       box%lo(1):box%hi(1))
            double precision, intent(inout) :: sqs(0:nz, box%lo(2):box%hi(2), &
                                                         box%lo(1):box%hi(1))
            double precision, intent(inout) :: qdi(0:nz, box%lo(2):box%hi(2), &
                                                         box%lo(1):box%hi(1))
            double precision, intent(inout) :: qdf(0:nz, box%lo(2):box%hi(2), &
                                                         box%lo(1):box%hi(1))
            integer                         :: iz

            ! Filter source:
            sqs = filt * sqs

            qdi = q
            q = (qdi + dt2 * sqs)
            call field_combine_semi_spectral(q)
            !$omp parallel do private(iz)  default(shared)
            do iz = 0, nz
                q(iz, :, :) = q(iz, :, :) * self%emq
            enddo
            !$omp end parallel do
            call field_decompose_semi_spectral(q)

            qdf = qdi + dt6 * sqs

        end subroutine impl_rk4_substep_one

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! First corrector step (t = t0 + dt/2):
        subroutine impl_rk4_substep_two(self, q, sqs, qdi, qdf)
            class(impl_rk4),  intent(inout) :: self
            double precision, intent(inout) :: q(0:nz, box%lo(2):box%hi(2),   &
                                                       box%lo(1):box%hi(1))
            double precision, intent(inout) :: sqs(0:nz, box%lo(2):box%hi(2), &
                                                         box%lo(1):box%hi(1))
            double precision, intent(in)    :: qdi(0:nz, box%lo(2):box%hi(2), &
                                                         box%lo(1):box%hi(1))
            double precision, intent(inout) :: qdf(0:nz, box%lo(2):box%hi(2), &
                                                         box%lo(1):box%hi(1))
            integer                         :: iz

            ! Filter source:
            sqs = filt * sqs

            ! apply integrating factors to source
            call field_combine_semi_spectral(sqs)
            !$omp parallel do private(iz)  default(shared)
            do iz = 0, nz
                sqs(iz, :, :) = sqs(iz, :, :) * self%epq
            enddo
            !$omp end parallel do
            call field_decompose_semi_spectral(sqs)

            q = (qdi + dt2 * sqs)

            call field_combine_semi_spectral(q)
            !$omp parallel do private(iz)  default(shared)
            do iz = 0, nz
                q(iz, :, :) = q(iz, :, :) * self%emq
            enddo
            !$omp end parallel do
            call field_decompose_semi_spectral(q)

            qdf = qdf + dt3 * sqs

        end subroutine impl_rk4_substep_two

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Second predictor step (t = t0 + dt):
        subroutine impl_rk4_substep_three(self, q, sqs, qdi, qdf, dt)
            class(impl_rk4),  intent(inout) :: self
            double precision, intent(inout) :: q(0:nz, box%lo(2):box%hi(2),   &
                                                       box%lo(1):box%hi(1))
            double precision, intent(inout) :: sqs(0:nz, box%lo(2):box%hi(2), &
                                                         box%lo(1):box%hi(1))
            double precision, intent(in)    :: qdi(0:nz, box%lo(2):box%hi(2), &
                                                         box%lo(1):box%hi(1))
            double precision, intent(inout) :: qdf(0:nz, box%lo(2):box%hi(2), &
                                                         box%lo(1):box%hi(1))
            double precision, intent(in)    :: dt
            integer                         :: iz

            ! Filter source:
            sqs = filt * sqs

            ! apply integrating factors to source
            call field_combine_semi_spectral(sqs)
            !$omp parallel do private(iz)  default(shared)
            do iz = 0, nz
                sqs(iz, :, :) = sqs(iz, :, :) * self%epq
            enddo
            !$omp end parallel do
            call field_decompose_semi_spectral(sqs)

            q = qdi + dt * sqs

            call field_combine_semi_spectral(q)
            !$omp parallel do private(iz)  default(shared)
            do iz = 0, nz
                q(iz, :, :) = q(iz, :, :) * self%emq
            enddo
            !$omp end parallel do
            call field_decompose_semi_spectral(q)

            qdf = qdf + dt3 * sqs
        end subroutine impl_rk4_substep_three

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Second corrector step (t = t0 + dt):
        subroutine impl_rk4_substep_four(self, q, sqs, qdf)
            class(impl_rk4),  intent(inout) :: self
            double precision, intent(inout) :: q(0:nz, box%lo(2):box%hi(2),   &
                                                       box%lo(1):box%hi(1))
            double precision, intent(inout) :: sqs(0:nz, box%lo(2):box%hi(2), &
                                                         box%lo(1):box%hi(1))
            double precision, intent(in)    :: qdf(0:nz, box%lo(2):box%hi(2), &
                                                         box%lo(1):box%hi(1))
            integer                         :: iz

            ! Filter source:
            sqs = filt * sqs

            ! apply integrating factors to source
            call field_combine_semi_spectral(sqs)
            !$omp parallel do private(iz)  default(shared)
            do iz = 0, nz
                sqs(iz, :, :) = sqs(iz, :, :) * self%epq
            enddo
            !$omp end parallel do
            call field_decompose_semi_spectral(sqs)

            q = qdf + dt6 * sqs

            call field_combine_semi_spectral(q)
            !$omp parallel do private(iz)  default(shared)
            do iz = 0, nz
                q(iz, :, :) = q(iz, :, :) * self%emq
            enddo
            !$omp end parallel do
            call field_decompose_semi_spectral(q)

        end subroutine impl_rk4_substep_four

end module impl_rk4_mod
