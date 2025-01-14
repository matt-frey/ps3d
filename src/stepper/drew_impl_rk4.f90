module drew_impl_rk4
    use inversion_mod, only : vor2vel, vorticity_tendency
    use model, only : layout
    use mpi_layout, only : box
    use parameters, only : nz
    use fields
    use diffusion
    use constants
    implicit none

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

    subroutine impl_rk4(t, dt)
        double precision, intent(inout) :: t
        double precision, intent(in)    :: dt
        double precision                :: dt2, dt3, dt6
        double precision                :: epq(box%lo(2):box%hi(2), &       ! exp(D * (t-t0))
                                               box%lo(1):box%hi(1))
        double precision                :: emq(box%lo(2):box%hi(2), &       ! exp(- D * (t-t0))
                                               box%lo(1):box%hi(1))
        double precision                :: qdi(0:nz, box%lo(2):box%hi(2), &
                                                     box%lo(1):box%hi(1), 3)
        double precision                :: qdf(0:nz, box%lo(2):box%hi(2), &
                                                     box%lo(1):box%hi(1), 3)
        integer                         :: iz, nc


        dt2 = dt / 2.0d0
        dt3 = dt / 3.0d0
        dt6 = dt / 6.0d0

        ! set integrating factors
        epq = exp(vdiss)
        emq = 1.0d0 / epq

        !------------------------------------------------------------------
        ! RK4 predictor step at time t0 + dt/2:
        qdi = svor
        svor = (qdi + dt2 * svorts)
        do nc = 1, 3
            !$omp parallel do private(iz)  default(shared)
            do iz = 0, nz
                svor(iz, :, :, nc) = svor(iz, :, :, nc) * emq
            enddo
            !$omp end parallel do
        enddo

        qdf = qdi + dt6 * svorts

        !------------------------------------------------------------------
        ! Invert and get new sources:
        call vor2vel
        call source

        !------------------------------------------------------------------
        !RK4 corrector step at time t0 + dt/2:
        t = t + dt2

        ! apply integrating factors to source
        do nc = 1, 3
            !$omp parallel do private(iz)  default(shared)
            do iz = 0, nz
                svorts(iz, :, :, nc) = epq * svorts(iz, :, :, nc)
            enddo
            !$omp end parallel do
        enddo

        svor = (qdi + dt2 * svorts)

        do nc = 1, 3
            !$omp parallel do private(iz)  default(shared)
            do iz = 0, nz
                svor(iz, :, :, nc) = emq * svor(iz, :, :, nc)
            enddo
            !$omp end parallel do
        enddo

        qdf = qdf + dt3 * svorts

        !------------------------------------------------------------------
        ! Invert and get new sources:
        call vor2vel
        call source

        !------------------------------------------------------------------
        !RK4 predictor step at time t0 + dt:
        t = t + dt2

        emq = emq ** 2

        ! apply integrating factors to source
        do nc = 1, 3
            !$omp parallel do private(iz)  default(shared)
            do iz = 0, nz
                svorts(iz, :, :, nc) = epq * svorts(iz, :, :, nc)
            enddo
            !$omp end parallel do
        enddo

        svor = qdi + dt * svorts

        do nc = 1, 3
            !$omp parallel do private(iz)  default(shared)
            do iz = 0, nz
                svor(iz, :, :, nc) = emq * svor(iz, :, :, nc)
            enddo
            !$omp end parallel do
        enddo

        qdf = qdf + dt3 * svorts

        !------------------------------------------------------------------
        ! Invert and get new sources:
        call vor2vel
        call source

        !------------------------------------------------------------------
        !RK4 corrector step at time t0 + dt:
        epq = epq**2

        ! apply integrating factors to source
        do nc = 1, 3
            !$omp parallel do private(iz)  default(shared)
            do iz = 0, nz
                svorts(iz, :, :, nc) = epq * svorts(iz, :, :, nc)
            enddo
            !$omp end parallel do
        enddo

        svor = qdf + dt6 * svorts

        do nc = 1, 3
            !$omp parallel do private(iz)  default(shared)
            do iz = 0, nz
                svor(iz, :, :, nc) = emq * svor(iz, :, :, nc)
            enddo
            !$omp end parallel do
        enddo

        ! Ensure zero global mean horizontal vorticity conservation:
        do nc = 1, 2
            call layout%adjust_decomposed_mean(svor(:, :, :, nc), ini_vor_mean(nc))
        enddo

    end subroutine impl_rk4

end module drew_impl_rk4
