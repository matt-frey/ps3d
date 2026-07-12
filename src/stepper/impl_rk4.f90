module impl_rk4_mod
    use model, only : layout
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
        ! vep = exp( D * (t-t0))
        ! vem = exp(-D * (t-t0))
        double precision, allocatable :: vep(:, :, :), vem(:, :, :)
        double precision, allocatable :: svorf(:, :, :, :), svori(:, :, :, :)
#ifdef ENABLE_BUOYANCY
        double precision, allocatable :: bep(:, :, :), bem(:, :, :)
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
        double precision                :: dfac

        dfac = f12 * vorch * dt

        !$omp parallel workshare
        vdop = dfac * vdiss
        !$omp end parallel workshare

#ifdef ENABLE_BUOYANCY
        dfac = f12 * bf * dt

        !$omp parallel workshare
        bdop = dfac * bdiss
        !$omp end parallel workshare
#endif

    end subroutine impl_rk4_set_diffusion

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine impl_rk4_setup(self)
        class(impl_rk4), intent(inout) :: self

        allocate(self%vep(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(self%vem(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(self%svorf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1), 3))
        allocate(self%svori(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1), 3))

#ifdef ENABLE_BUOYANCY
        allocate(self%bep(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(self%bem(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
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

        !Define integrating factors
        self%vep = exp(vdop)
        self%vem = 1.0d0 / self%vep

#ifdef ENABLE_BUOYANCY
        self%bep = exp(bdop)
        self%bem = 1.0d0 / self%bep
#endif

        !------------------------------------------------------------------
        ! RK4 predictor step at time t0 + dt/2:
#ifdef ENABLE_BUOYANCY
        call self%impl_rk4_substep_one(q=sbuoy,          &
                                       sqs=sbuoys,       &
                                       qdi=self%sbuoyi,  &
                                       qdf=self%sbuoyf,  &
                                       mq=self%bem)
#endif

        do nc = 1, 3
            call self%impl_rk4_substep_one(q=svor(:, :, :, nc),          &
                                           sqs=svorts(:, :, :, nc),      &
                                           qdi=self%svori(:, :, :, nc),  &
                                           qdf=self%svorf(:, :, :, nc),  &
                                           mq=self%vem)
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
                                       mq=self%bem,      &
                                       pq=self%bep)
#endif

        do nc = 1, 3
            call self%impl_rk4_substep_two(q=svor(:, :, :, nc),          &
                                           sqs=svorts(:, :, :, nc),      &
                                           qdi=self%svori(:, :, :, nc),  &
                                           qdf=self%svorf(:, :, :, nc),  &
                                           mq=self%vem,                  &
                                           pq=self%vep)
        enddo

        !------------------------------------------------------------------
        ! Invert and get new sources:
        call vor2vel
        call source

        !------------------------------------------------------------------
        !RK4 predictor step at time t0 + dt:
        t = t + dt2

#ifdef ENABLE_BUOYANCY
        self%bem = self%bem ** 2

        call self%impl_rk4_substep_three(q=sbuoy,          &
                                         sqs=sbuoys,       &
                                         qdi=self%sbuoyi,  &
                                         qdf=self%sbuoyf,  &
                                         mq=self%bem,      &
                                         pq=self%bep,      &
                                         dt=dt)
#endif

        self%vem = self%vem ** 2

        do nc = 1, 3
            call self%impl_rk4_substep_three(q=svor(:, :, :, nc),          &
                                             sqs=svorts(:, :, :, nc),      &
                                             qdi=self%svori(:, :, :, nc),  &
                                             qdf=self%svorf(:, :, :, nc),  &
                                             mq=self%vem,                  &
                                             pq=self%vep,                  &
                                             dt=dt)
        enddo


        !------------------------------------------------------------------
        ! Invert and get new sources:
        call vor2vel
        call source

        !------------------------------------------------------------------
        !RK4 corrector step at time t0 + dt:

#ifdef ENABLE_BUOYANCY
        self%bep = self%bep ** 2

        call self%impl_rk4_substep_four(q=sbuoy,          &
                                        sqs=sbuoys,       &
                                        qdf=self%sbuoyf,  &
                                        mq=self%bem,      &
                                        pq=self%bep)
#endif

        self%vep = self%vep ** 2

        do nc = 1, 3
            call self%impl_rk4_substep_four(q=svor(:, :, :, nc),          &
                                            sqs=svorts(:, :, :, nc),      &
                                            qdf=self%svorf(:, :, :, nc),  &
                                            mq=self%vem,                  &
                                            pq=self%vep)
        enddo

        ! Ensure zero global mean horizontal vorticity conservation:
        do nc = 1, 2
           call layout%adjust_semi_spectral_mean(svor(:, :, :, nc), &
                                                 ini_vor_mean(nc))
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
        double precision, intent(in)    ::  mq(0:nz, box%lo(2):box%hi(2), &
                                                     box%lo(1):box%hi(1))

        !Use mixed-spectral space to apply 3d diffusion operator:
        call layout%decompose_semi_spectral(q)
        call layout%decompose_semi_spectral(sqs)

        qdi = q

        !Apply integrating factor to source
        q = (qdi + dt2 * sqs) * mq
        !qdi & sqs are in mixed-spectral space, so q is automatically

        !Return field q to semi-spectral space for use in vor2vel & source:
        call layout%combine_semi_spectral(q)

        qdf = qdi + dt6 * sqs
        !qdf is in mixed-spectral space on exit

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
        double precision, intent(in)    ::  mq(0:nz, box%lo(2):box%hi(2), &
                                                     box%lo(1):box%hi(1))
        double precision, intent(in)    ::  pq(0:nz, box%lo(2):box%hi(2), &
                                                     box%lo(1):box%hi(1))

        !Use mixed-spectral space to apply 3d diffusion operator:
        call layout%decompose_semi_spectral(sqs)

        !Apply integrating factor to source
        sqs = pq * sqs

        !qdi & sqs are in mixed-spectral space, so q is automatically
        q = mq * (qdi + dt2 * sqs)

        !Return field q to semi-spectral space for use in vor2vel & source:
        call layout%combine_semi_spectral(q)

        qdf = qdf + dt3 * sqs
        !qdf is in mixed-spectral space on exit

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
        double precision, intent(in)    ::  mq(0:nz, box%lo(2):box%hi(2), &
                                                     box%lo(1):box%hi(1))
        double precision, intent(in)    ::  pq(0:nz, box%lo(2):box%hi(2), &
                                                     box%lo(1):box%hi(1))
        double precision, intent(in)    :: dt

        !Use mixed-spectral space to apply 3d diffusion operator:
        call layout%decompose_semi_spectral(sqs)

        !Apply integrating factor to source
        sqs = pq * sqs

        !qdi & sqs are in mixed-spectral space, so q is automatically
        q = mq * (qdi + dt * sqs)

        !Return field q to semi-spectral space for use in vor2vel & source:
        call layout%combine_semi_spectral(q)

        qdf = qdf + dt3 * sqs
        !qdf is in mixed-spectral space on exit

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
        double precision, intent(in)    ::  mq(0:nz, box%lo(2):box%hi(2), &
                                                     box%lo(1):box%hi(1))
        double precision, intent(in)    ::  pq(0:nz, box%lo(2):box%hi(2), &
                                                     box%lo(1):box%hi(1))

        !Use mixed-spectral space to apply 3d diffusion operator:
        call layout%decompose_semi_spectral(sqs)

        !Apply integrating factor to source
        sqs = pq * sqs

        !qdf & sqs are in mixed-spectral space, so q is automatically
        q = mq * (qdf + dt6 * sqs)

        !Return field q to semi-spectral space for use elsewhere:
        call layout%combine_semi_spectral(q)

    end subroutine impl_rk4_substep_four

end module impl_rk4_mod
