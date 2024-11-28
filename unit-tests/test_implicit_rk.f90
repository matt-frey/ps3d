! =============================================================================
!                               Test implicit RK
!
!   This test assumes S = constant and D = constant. The analytic solution
!   is q(t_0+dt)=q(t_0)E + (1-E)S/D where E = exp(-D dt). The numerical solution
!   (by hand) is q(t_0)E + (dt/6)(1 + 4E^{1/2} + E)S where of course
!   E^{1/2} = exp(-D dt/2).  The leading terms agree exactly (i.e. q(t_0)E),
!   and we find the difference is O(dt^5) as expected.
!   This test uses S = D = 1 and dt = 0.1.
!   For time steps dt = 0.01, 0.02, 0.05 and 0.1, a convergence of dt^5 is
!   achieved.
! =============================================================================
program test_implicit_rk
    use unit_test
    use options, only : vor_visc
    use constants, only : zero, one, f12
    use parameters, only : lower, update_parameters, nx, ny, nz, extent
    use fields
    use inversion_mod, only : vor2vel, vor2vel_timer, vorticity_tendency, vtend_timer
    use field_diagnostics
    use mpi_timer
    use mpi_environment
    use mpi_layout
    use mpi_collectives, only : mpi_blocking_reduce
    use mpi_utils, only : mpi_stop
    use model_factory, only : layout
    implicit none

    double precision :: error
    integer          :: iz, nc, n
    double precision, allocatable :: ref(:, :, :, :), ed(:, :)
    double precision, allocatable :: src(:, :, :, :)
    double precision :: time, time_step

    call mpi_env_initialise

    call register_timer('vorticity', vor2vel_timer)
    call register_timer('vtend', vtend_timer)

    vor_visc%nnu = 3
    vor_visc%prediss = 30

    nx = 32
    ny = 32
    nz = 32

    lower  = (/zero, zero, zero/)
    extent = (/one, one, one/)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    allocate(ref(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1), 3))
    allocate(src(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1), 3))
    allocate(ed(box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

    call update_parameters

    call field_default

!     call init_inversion

    time_step = 0.1d0
    do n = 1, 4
        time = zero
        svor = 2.0d0
        ref = svor
        src = one
        svorts = src

        vdiss = one

        ed = dexp(-time_step * vdiss)
        do nc = 1, 3
            call layout%combine_semi_spectral(ref(:, :, :, nc))
            call layout%combine_semi_spectral(src(:, :, :, nc))
            do iz = 0, nz
                ref(iz, :, :, nc) = ed * ref(iz, :, :, nc) + (one - ed) * src(iz, :, :, nc) / vdiss
            enddo
            call layout%decompose_semi_spectral(ref(:, :, :, nc))
            call layout%decompose_semi_spectral(src(:, :, :, nc))
        enddo

        call impl_rk4(time, time_step)

        error = maxval(dabs(ref - svor))

        time_step = f12 * time_step

    enddo

    call mpi_blocking_reduce(error, MPI_MAX, world)

    if (world%rank == world%root) then
        call print_result_dp('Test implicit RK', error, atol=1.0e-12)
    endif

    deallocate(ref)
    deallocate(src)
    deallocate(ed)

    call mpi_env_finalise


    contains

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


            dt2 = dt / 2.0d0
            dt3 = dt / 3.0d0
            dt6 = dt / 6.0d0

            ! set integrating factors
            emq = dexp(- dt2 * vdiss)
            epq = 1.0d0 / emq

            qdi = svor
            svor = (qdi + dt2 * svorts)
            do nc = 1, 3
                call layout%combine_semi_spectral(svor(:, :, :, nc))
                !$omp parallel do private(iz)  default(shared)
                do iz = 0, nz
                    svor(iz, :, :, nc) = svor(iz, :, :, nc) * emq
                enddo
                !$omp end parallel do
                call layout%decompose_semi_spectral(svor(:, :, :, nc))
            enddo


            qdf = qdi + dt6 * svorts

            !------------------------------------------------------------------
            !RK4 corrector step at time t0 + dt/2:
            t = t + dt2

            call vor2vel

            call vorticity_tendency
            svorts = src
            ! apply integrating factors to source
            do nc = 1, 3
                call layout%combine_semi_spectral(svorts(:, :, :, nc))
                !$omp parallel do private(iz)  default(shared)
                do iz = 0, nz
                    svorts(iz, :, :, nc) = svorts(iz, :, :, nc) * epq
                enddo
                !$omp end parallel do
                call layout%decompose_semi_spectral(svorts(:, :, :, nc))
            enddo

            svor = (qdi + dt2 * svorts)
            do nc = 1, 3
                call layout%combine_semi_spectral(svor(:, :, :, nc))
                !$omp parallel do private(iz)  default(shared)
                do iz = 0, nz
                    svor(iz, :, :, nc) = svor(iz, :, :, nc) * emq
                enddo
                !$omp end parallel do
                call layout%decompose_semi_spectral(svor(:, :, :, nc))
            enddo

            qdf = qdf + dt3 * svorts

            !------------------------------------------------------------------
            !RK4 predictor step at time t0 + dt:

            !Invert PV and compute velocity:
            call vor2vel

            call vorticity_tendency
            svorts = src
            ! apply integrating factors to source
            do nc = 1, 3
                call layout%combine_semi_spectral(svorts(:, :, :, nc))
                !$omp parallel do private(iz)  default(shared)
                do iz = 0, nz
                    svorts(iz, :, :, nc) = svorts(iz, :, :, nc) * epq
                enddo
                !$omp end parallel do
                call layout%decompose_semi_spectral(svorts(:, :, :, nc))
            enddo

            emq = emq**2
            svor = qdi + dt * svorts
            do nc = 1, 3
                call layout%combine_semi_spectral(svor(:, :, :, nc))
                !$omp parallel do private(iz)  default(shared)
                do iz = 0, nz
                    svor(iz, :, :, nc) = svor(iz, :, :, nc) * emq
                enddo
                !$omp end parallel do
                call layout%decompose_semi_spectral(svor(:, :, :, nc))
            enddo

            qdf = qdf + dt3 * svorts

            !------------------------------------------------------------------
            !RK4 corrector step at time t0 + dt:
            t = t + dt2

            call vor2vel
            call vorticity_tendency
            svorts = src

            epq = epq**2
            ! apply integrating factors to source
            do nc = 1, 3
                call layout%combine_semi_spectral(svorts(:, :, :, nc))
                !$omp parallel do private(iz)  default(shared)
                do iz = 0, nz
                    svorts(iz, :, :, nc) = svorts(iz, :, :, nc) * epq
                enddo
                !$omp end parallel do
                call layout%decompose_semi_spectral(svorts(:, :, :, nc))
            enddo

            svor = qdf + dt6 * svorts
            do nc = 1, 3
                call layout%combine_semi_spectral(svor(:, :, :, nc))
                !$omp parallel do private(iz)  default(shared)
                do iz = 0, nz
                    svor(iz, :, :, nc) = svor(iz, :, :, nc) * emq
                enddo
                !$omp end parallel do
                call layout%decompose_semi_spectral(svor(:, :, :, nc))
            enddo

        end subroutine impl_rk4

end program test_implicit_rk
