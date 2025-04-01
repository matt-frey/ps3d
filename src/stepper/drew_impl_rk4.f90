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

    ! epq = exp( D * (t-t0))
    ! emq = exp(-D * (t-t0))
    double precision, allocatable :: epq(:, :), emq(:, :)
    double precision, allocatable :: svorf(:, :, :, :), svori(:, :, :, :)
#ifdef ENABLE_BUOYANCY
    double precision, allocatable :: bpq(:, :), bmq(:, :)
    double precision, allocatable :: sbuoyf(:, :, :), sbuoyi(:, :, :)
#endif


contains

    subroutine set_diffusion(dt, vorch, bf)
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

    end subroutine set_diffusion

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine impl_rk4_setup
        allocate(epq(box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(emq(box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(svorf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1), 3))
        allocate(svori(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1), 3))

#ifdef ENABLE_BUOYANCY
        allocate(bpq(box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(bmq(box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(sbuoyf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(sbuoyi(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
#endif

    end subroutine impl_rk4_setup

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine impl_rk4(t, dt)
        double precision, intent(inout) :: t
        double precision, intent(in)    :: dt
        double precision                :: nu, kappa
        integer                         :: nc

        if (.not. allocated(epq)) then
            call impl_rk4_setup
        endif

        ! set viscocity/diffusivity
        nu =    vvisc*(1.d0 + 0.0*(0.07/dt - 1.0d0))
        kappa = bvisc*(1.d0 + 0.0*(0.07/dt - 1.0d0))

        dt2 = f12 * dt
        dt3 = f13 * dt
        dt6 = f16 * dt

        ! set integrating factors
        epq = 1.d0  !!exp(vdiss)
        emq = 1.0d0 / epq

#ifdef ENABLE_BUOYANCY
        bpq = 1.0d0 !!exp(bdiss)
        bmq = 1.0d0 / bpq
#endif

        !------------------------------------------------------------------
        ! RK4 predictor step at time t0 + dt/2:
#ifdef ENABLE_BUOYANCY
        call impl_rk4_substep_one(q=sbuoy,     &
                                  sqs=sbuoys,  &
                                  qdi=sbuoyi,  &
                                  qdf=sbuoyf,  &
                                  mq=bmq)
#endif

        do nc = 1, 3
            call impl_rk4_substep_one(q=svor(:, :, :, nc),     &
                                      sqs=svorts(:, :, :, nc), &
                                      qdi=svori(:, :, :, nc),  &
                                      qdf=svorf(:, :, :, nc),  &
                                      mq=emq)
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
                                  mq=bmq,      &
                                  pq=bpq)
#endif

        do nc = 1, 3
            call impl_rk4_substep_two(q=svor(:, :, :, nc),     &
                                      sqs=svorts(:, :, :, nc), &
                                      qdi=svori(:, :, :, nc),  &
                                      qdf=svorf(:, :, :, nc),  &
                                      mq=emq,                  &
                                      pq=epq)
        enddo

        !------------------------------------------------------------------
        ! Invert and get new sources:
        call vor2vel
        call source

        !------------------------------------------------------------------
        !RK4 predictor step at time t0 + dt:
        t = t + dt2

        emq = emq ** 2


#ifdef ENABLE_BUOYANCY
        bmq = bmq ** 2

        call impl_rk4_substep_three(q=sbuoy,     &
                                    sqs=sbuoys,  &
                                    qdi=sbuoyi,  &
                                    qdf=sbuoyf,  &
                                    mq=bmq,      &
                                    pq=bpq,      &
                                    dt=dt)
#endif

        do nc = 1, 3
            call impl_rk4_substep_three(q=svor(:, :, :, nc),     &
                                        sqs=svorts(:, :, :, nc), &
                                        qdi=svori(:, :, :, nc),  &
                                        qdf=svorf(:, :, :, nc),  &
                                        mq=emq,                  &
                                        pq=epq,                  &
                                        dt=dt)
        enddo


        !------------------------------------------------------------------
        ! Invert and get new sources:
        call vor2vel
        call source

        !------------------------------------------------------------------
        !RK4 corrector step at time t0 + dt:

        epq = epq ** 2

#ifdef ENABLE_BUOYANCY
        bpq = bpq ** 2

        call impl_rk4_substep_four(q=sbuoy,     &
                                   sqs=sbuoys,  &
                                   qdf=sbuoyf,  &
                                   mq=bmq,      &
                                   pq=bpq)
#endif

        do nc = 1, 3
            call impl_rk4_substep_four(q=svor(:, :, :, nc),     &
                                       sqs=svorts(:, :, :, nc), &
                                       qdf=svorf(:, :, :, nc),  &
                                       mq=emq,                  &
                                       pq=epq)
        enddo
        !!!!!   DONE ADVECTION STEP
#ifdef ENABLE_BUOYANCY
          !call layout%zdiffNF(sbuoy,dt,kappa,kappa)
          call layout%zdiffuse(sbuoy,dt,kappa,kappa)
#endif
        do nc = 1, 3
          call layout%zdiffuse(svor(:,:,:,nc),dt,nu,1.0*nu)
        enddo

        
        ! Ensure zero global mean horizontal vorticity conservation:
 !       do nc = 1, 2
 !          call layout%adjust_decomposed_mean(svor(:, :, :, nc), ini_vor_mean(nc))
 !      enddo

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
        double precision, intent(in)    :: mq(box%lo(2):box%hi(2), &
                                              box%lo(1):box%hi(1))
        integer                         :: iz

        qdi = q
        q = (qdi + dt2 * sqs)
        !$omp parallel do private(iz)  default(shared)
        !!do iz = 0, nz
        !!    q(iz, :, :) = q(iz, :, :) * mq
        !!enddo
        !$omp end parallel do

        qdf = qdi + dt6 * sqs

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
    subroutine impl_rk4_substep_three(q, sqs, qdi, qdf, mq, pq, dt)
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
    subroutine impl_rk4_substep_four(q, sqs, qdf, mq, pq)
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

end module drew_impl_rk4
