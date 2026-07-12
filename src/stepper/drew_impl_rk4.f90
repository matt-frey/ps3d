module drew_impl_rk4
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

    ! vep = exp( D * (t-t0))
    ! vem = exp(-D * (t-t0))
    double precision, allocatable :: vep(:, :, :), vem(:, :, :)
    double precision, allocatable :: svorf(:, :, :, :), svori(:, :, :, :)
#ifdef ENABLE_BUOYANCY
    double precision, allocatable :: bep(:, :, :), bem(:, :, :)
    double precision, allocatable :: sbuoyf(:, :, :), sbuoyi(:, :, :)
#endif

contains

    subroutine set_diffusion(dt, vorch, bf)
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
        !$omp end parallel workshare

    end subroutine set_diffusion

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine impl_rk4_setup
        allocate(vep(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(vem(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(svorf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1), 3))
        allocate(svori(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1), 3))

#ifdef ENABLE_BUOYANCY
        allocate(bep(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(bem(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(sbuoyf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(sbuoyi(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
#endif

    end subroutine impl_rk4_setup

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine impl_rk4(t, dt)
        double precision, intent(inout) :: t
        double precision, intent(in)    :: dt
        integer                         :: nc

        if (.not. allocated(bep)) then
            call impl_rk4_setup
        endif

        dt2 = f12 * dt
        dt3 = f13 * dt
        dt6 = f16 * dt

        !Define integrating factors
        vep = exp(vdop)
        vem = 1.0d0 / vep

#ifdef ENABLE_BUOYANCY
        bep = exp(bdop)
        bem = 1.0d0 / bep
#endif

        !------------------------------------------------------------------
        ! RK4 predictor step at time t0 + dt/2:
#ifdef ENABLE_BUOYANCY
        call impl_rk4_substep_one(q=sbuoy,     &
                                  sqs=sbuoys,  &
                                  qdi=sbuoyi,  &
                                  qdf=sbuoyf,  &
                                  mq=bem)
#endif

        do nc = 1, 3
            call impl_rk4_substep_one(q=svor(:, :, :, nc),     &
                                      sqs=svorts(:, :, :, nc), &
                                      qdi=svori(:, :, :, nc),  &
                                      qdf=svorf(:, :, :, nc),  &
                                      mq=vem)
        enddo


        !------------------------------------------------------------------
        ! Invert and get new sources:
        call vor2vel
        call source

        !------------------------------------------------------------------
        !RK4 corrector step at time t0 + dt/2:
        t = t + dt2

#ifdef ENABLE_BUOYANCY
        call impl_rk4_substep_two(q=sbuoy,     &
                                  sqs=sbuoys,  &
                                  qdi=sbuoyi,  &
                                  qdf=sbuoyf,  &
                                  mq=bem,      &
                                  pq=bep)
#endif

        do nc = 1, 3
            call impl_rk4_substep_two(q=svor(:, :, :, nc),     &
                                      sqs=svorts(:, :, :, nc), &
                                      qdi=svori(:, :, :, nc),  &
                                      qdf=svorf(:, :, :, nc),  &
                                      mq=vem,                  &
                                      pq=vep)
        enddo

        !------------------------------------------------------------------
        ! Invert and get new sources:
        call vor2vel
        call source

        !------------------------------------------------------------------
        !RK4 predictor step at time t0 + dt:
        t = t + dt2

#ifdef ENABLE_BUOYANCY
        bem = bem ** 2

        call impl_rk4_substep_three(q=sbuoy,     &
                                    sqs=sbuoys,  &
                                    qdi=sbuoyi,  &
                                    qdf=sbuoyf,  &
                                    mq=bem,      &
                                    pq=bep,      &
                                    dt=dt)
#endif

        vem = vem ** 2

        do nc = 1, 3
            call impl_rk4_substep_three(q=svor(:, :, :, nc),     &
                                        sqs=svorts(:, :, :, nc), &
                                        qdi=svori(:, :, :, nc),  &
                                        qdf=svorf(:, :, :, nc),  &
                                        mq=vem,                  &
                                        pq=vep,                  &
                                        dt=dt)
        enddo

        !------------------------------------------------------------------
        ! Invert and get new sources:
        call vor2vel
        call source

        !------------------------------------------------------------------
        !RK4 corrector step at time t0 + dt:

#ifdef ENABLE_BUOYANCY
        bep = bep ** 2

        call impl_rk4_substep_four(q=sbuoy,     &
                                   sqs=sbuoys,  &
                                   qdf=sbuoyf,  &
                                   mq=bem,      &
                                   pq=bep)
#endif

        vep = vep ** 2

        do nc = 1, 3
            call impl_rk4_substep_four(q=svor(:, :, :, nc),     &
                                       sqs=svorts(:, :, :, nc), &
                                       qdf=svorf(:, :, :, nc),  &
                                       mq=vem,                  &
                                       pq=vep)
        enddo
        !!!!!   DONE ADVECTION STEP

        ! Ensure zero global mean horizontal vorticity conservation:
       do nc = 1, 2
          call layout%adjust_semi_spectral_mean(svor(:, :, :, nc), &
                                                ini_vor_mean(nc))
       enddo

    end subroutine impl_rk4

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Initialisation step (t = t0) (predictor):
    subroutine impl_rk4_substep_one(q, sqs, qdi, qdf, mq)
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
    subroutine impl_rk4_substep_two(q, sqs, qdi, qdf, mq, pq)
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
    subroutine impl_rk4_substep_three(q, sqs, qdi, qdf, mq, pq, dt)
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
    subroutine impl_rk4_substep_four(q, sqs, qdf, mq, pq)
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

end module drew_impl_rk4
