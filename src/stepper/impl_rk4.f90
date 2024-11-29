module impl_rk4_mod
    use model, only : layout, filter
    use stepper_mod, only : stepper_t
    use constants, only : f12, f13, f16
    use parameters, only : nz
    use fields
    use diffusion
    use inversion_mod, only : vor2vel, source
    use field_diagnostics
    implicit none

    double precision :: dt2, dt3, dt6

    type, extends(stepper_t) :: impl_rk4
        ! epq = exp( D * (t-t0))
        ! emq = exp(-D * (t-t0))
        double precision, allocatable :: epq(:, :), emq(:, :)
        double precision, allocatable :: svorf(:, :, :, :), svori(:, :, :, :)
#ifdef ENABLE_BUOYANCY
        double precision, allocatable :: bpq(:, :), bmq(:, :)
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

        subroutine impl_rk4_set_diffusion(self, dt, vorch, bf)
            class(impl_rk4),  intent(inout) :: self
            double precision, intent(in)    :: dt
            double precision, intent(in)    :: vorch, bf
            double precision                :: dfac, dbac

            dfac = f12 * vorch * dt
            dbac = f12 * bf * dt

            !$omp parallel workshare
            vdiss = dfac * vhdis
#ifdef ENABLE_BUOYANCY
            bdiss = dbac * bhdis
#endif
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
            allocate(self%bpq(box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
            allocate(self%bmq(box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
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
            self%epq = dexp(vdiss)
            self%emq = 1.0d0 / self%epq
            call filter%apply2d(self%epq)

#ifdef ENABLE_BUOYANCY
            self%bpq = dexp(bdiss)
            self%bmq = 1.0d0 / self%bpq
            call filter%apply2d(self%bpq)
#endif

            !------------------------------------------------------------------
            ! RK4 predictor step at time t0 + dt/2:
#ifdef ENABLE_BUOYANCY
            call self%impl_rk4_substep_one(q=sbuoy,          &
                                           sqs=sbuoys,       &
                                           qdi=self%sbuoyi,  &
                                           qdf=self%sbuoyf,  &
                                           mq=self%bmq)
#endif

            do nc = 1, 3
                call self%impl_rk4_substep_one(q=svor(:, :, :, nc),          &
                                               sqs=svorts(:, :, :, nc),      &
                                               qdi=self%svori(:, :, :, nc),  &
                                               qdf=self%svorf(:, :, :, nc),  &
                                               mq=self%emq)
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
                                           qdf=self%sbuoyf,  &
                                           mq=self%bmq,      &
                                           pq=self%bpq)
#endif

            do nc = 1, 3
                call self%impl_rk4_substep_two(q=svor(:, :, :, nc),          &
                                               sqs=svorts(:, :, :, nc),      &
                                               qdi=self%svori(:, :, :, nc),  &
                                               qdf=self%svorf(:, :, :, nc),  &
                                               mq=self%emq,                  &
                                               pq=self%epq)
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
            self%bmq = self%bmq ** 2

            call self%impl_rk4_substep_three(q=sbuoy,          &
                                             sqs=sbuoys,       &
                                             qdi=self%sbuoyi,  &
                                             qdf=self%sbuoyf,  &
                                             mq=self%bmq,      &
                                             pq=self%bpq,      &
                                             dt=dt)
#endif

            do nc = 1, 3
                call self%impl_rk4_substep_three(q=svor(:, :, :, nc),          &
                                                 sqs=svorts(:, :, :, nc),      &
                                                 qdi=self%svori(:, :, :, nc),  &
                                                 qdf=self%svorf(:, :, :, nc),  &
                                                 mq=self%emq,                  &
                                                 pq=self%epq,                  &
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
            self%bpq = self%bpq ** 2

            call self%impl_rk4_substep_four(q=sbuoy,          &
                                            sqs=sbuoys,       &
                                            qdf=self%sbuoyf,  &
                                            mq=self%bmq,      &
                                            pq=self%bpq)
#endif

            do nc = 1, 3
                call self%impl_rk4_substep_four(q=svor(:, :, :, nc),          &
                                                sqs=svorts(:, :, :, nc),      &
                                                qdf=self%svorf(:, :, :, nc),  &
                                                mq=self%emq,                  &
                                                pq=self%epq)
            enddo

            ! Ensure zero global mean horizontal vorticity conservation:
            do nc = 1, 2
                call layout%adjust_decomposed_mean(svor(:, :, :, nc), ini_vor_mean(nc))
            enddo

        end subroutine impl_rk4_step

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Initialisation step (t = t0) (predictor):
        subroutine impl_rk4_substep_one(self, q, sqs, qdi, qdf, mq)
            class(impl_rk4),  intent(inout) :: self
            double precision, intent(inout) :: q(0:nz, box%lo(2):box%hi(2), &
                                                       box%lo(1):box%hi(1))
            double precision, intent(inout) :: sqs(0:nz, box%lo(2):box%hi(2), &
                                                         box%lo(1):box%hi(1))
            double precision, intent(inout) :: qdi(0:nz, box%lo(2):box%hi(2), &
                                                         box%lo(1):box%hi(1))
            double precision, intent(inout) :: qdf(0:nz, box%lo(2):box%hi(2), &
                                                         box%lo(1):box%hi(1))
            double precision, intent(in)    :: mq(box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1))
            integer                         :: iz

            ! Filter source
            do iz = 0, nz
                call filter%apply2d(sqs(iz, :, :))
            enddo
            qdi = q
            q = (qdi + dt2 * sqs)
            !$omp parallel do private(iz)  default(shared)
            do iz = 0, nz
                q(iz, :, :) = q(iz, :, :) * mq
            enddo
            !$omp end parallel do

            qdf = qdi + dt6 * sqs

        end subroutine impl_rk4_substep_one

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! First corrector step (t = t0 + dt/2):
        subroutine impl_rk4_substep_two(self, q, sqs, qdi, qdf, mq, pq)
            class(impl_rk4),  intent(inout) :: self
            double precision, intent(inout) :: q(0:nz, box%lo(2):box%hi(2),   &
                                                       box%lo(1):box%hi(1))
            double precision, intent(inout) :: sqs(0:nz, box%lo(2):box%hi(2), &
                                                         box%lo(1):box%hi(1))
            double precision, intent(in)    :: qdi(0:nz, box%lo(2):box%hi(2), &
                                                         box%lo(1):box%hi(1))
            double precision, intent(inout) :: qdf(0:nz, box%lo(2):box%hi(2), &
                                                         box%lo(1):box%hi(1))
            double precision, intent(in)    :: mq(box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1))
            double precision, intent(in)    :: pq(box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1))
            integer                         :: iz

            ! apply integrating factors to source
            !$omp parallel do private(iz)  default(shared)
            do iz = 0, nz
                sqs(iz, :, :) = pq * sqs(iz, :, :)
            enddo
            !$omp end parallel do

            q = qdi + dt2 * sqs

            !$omp parallel do private(iz)  default(shared)
            do iz = 0, nz
                q(iz, :, :) = mq * q(iz, :, :)
            enddo
            !$omp end parallel do

            qdf = qdf + dt3 * sqs

        end subroutine impl_rk4_substep_two

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Second predictor step (t = t0 + dt):
        subroutine impl_rk4_substep_three(self, q, sqs, qdi, qdf, mq, pq, dt)
            class(impl_rk4),  intent(inout) :: self
            double precision, intent(inout) :: q(0:nz, box%lo(2):box%hi(2),   &
                                                       box%lo(1):box%hi(1))
            double precision, intent(inout) :: sqs(0:nz, box%lo(2):box%hi(2), &
                                                         box%lo(1):box%hi(1))
            double precision, intent(in)    :: qdi(0:nz, box%lo(2):box%hi(2), &
                                                         box%lo(1):box%hi(1))
            double precision, intent(inout) :: qdf(0:nz, box%lo(2):box%hi(2), &
                                                         box%lo(1):box%hi(1))
            double precision, intent(in)    :: mq(box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1))
            double precision, intent(in)    :: pq(box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1))
            double precision, intent(in)    :: dt
            integer                         :: iz

            ! apply integrating factors to source
            !$omp parallel do private(iz)  default(shared)
            do iz = 0, nz
                sqs(iz, :, :) = pq * sqs(iz, :, :)
            enddo
            !$omp end parallel do

            q = qdi + dt * sqs

            !$omp parallel do private(iz)  default(shared)
            do iz = 0, nz
                q(iz, :, :) = mq * q(iz, :, :)
            enddo
            !$omp end parallel do

            qdf = qdf + dt3 * sqs
        end subroutine impl_rk4_substep_three

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Second corrector step (t = t0 + dt):
        subroutine impl_rk4_substep_four(self, q, sqs, qdf, mq, pq)
            class(impl_rk4),  intent(inout) :: self
            double precision, intent(inout) :: q(0:nz, box%lo(2):box%hi(2),   &
                                                       box%lo(1):box%hi(1))
            double precision, intent(inout) :: sqs(0:nz, box%lo(2):box%hi(2), &
                                                         box%lo(1):box%hi(1))
            double precision, intent(in)    :: qdf(0:nz, box%lo(2):box%hi(2), &
                                                         box%lo(1):box%hi(1))
            double precision, intent(in)    :: mq(box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1))
            double precision, intent(in)    :: pq(box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1))
            integer                         :: iz

            ! apply integrating factors to source
            !$omp parallel do private(iz)  default(shared)
            do iz = 0, nz
                sqs(iz, :, :) = pq * sqs(iz, :, :)
            enddo
            !$omp end parallel do

            q = qdf + dt6 * sqs

            !$omp parallel do private(iz)  default(shared)
            do iz = 0, nz
                q(iz, :, :) = mq * q(iz, :, :)
            enddo
            !$omp end parallel do

        end subroutine impl_rk4_substep_four

end module impl_rk4_mod
