! Advances fields from time t to t+dt using an iterative implicit
! trapezoidal method of the form
!
!     (F^{n+1}-F^n)/dt = (L^{n+1}+L^n)/2 + (S^{n+1}+S^n)/2
!
! for a field F, where n refers to the time level, L[F] refers to
! the linear dissipation terms (hyperdiffusion), and S[F] refers to
! the remaining source terms.

! We start with the guess S^{n+1} = S^n and iterate  niter  times
! (see parameter statement below).
module cn2
    use constants, only : f12
    use parameters, only : nz
    use fields
    use inversion_utils
    use inversion_mod, only : vor2vel, source
    use field_diagnostics
    implicit none

    ! Number of iterations of above scheme:
    integer, parameter:: niter = 2

    double precision :: dt2

    contains

        subroutine cn2_step(t, dt)
            double precision, intent(inout) :: t
            double precision, intent(in)    :: dt
            integer                         :: iter
            integer                         :: nc

            !Update value of dt/2:
            dt2 = f12 * dt

            !------------------------------------------------------------------
            !Start with a guess for F^{n+1} for all fields:

            !Initialise iteration (dt = dt/2 below):
#ifdef ENABLE_BUOYANCY
            bsm = sbuoy + dt2 * sbuoys
            sbuoy = filt * (bsm + dt2 * sbuoys)
            !call field_combine_semi_spectral(sbuoy)
            !! omp parallel do private(iz)  default(shared)
            !do iz = 0, nz
            !    sbuoy(iz, :, :) = diss * sbuoy(iz, :, :)
            !enddo
            !! omp end parallel do
            !call field_decompose_semi_spectral(sbuoy)
#endif

            ! Advance interior and boundary values of vorticity
            !$omp parallel workshare
            vortsm = svor + dt2 * svorts
            !$omp end parallel workshare

            do nc = 1, 3
                !$omp parallel workshare
                svor(:, :, :, nc) = filt * (vortsm(:, :, :, nc) + dt2 * svorts(:, :, :, nc))
                !$omp end parallel workshare
                !call field_combine_semi_spectral(svor(:, :, :, nc))
                !!omp parallel do private(iz)  default(shared)
                !do iz = 0, nz
                !    svor(iz, :, :, nc) = diss * svor(iz, :, :, nc)
                !enddo
                !!omp end parallel do
                !call field_decompose_semi_spectral(svor(:, :, :, nc))
            enddo

            call adjust_vorticity_mean

            !diss is related to the hyperdiffusive operator (see end of adapt)

            !------------------------------------------------------------------
            !Iterate to improve estimates of F^{n+1}:
            do iter = 1, niter
                !Perform inversion at t^{n+1} from estimated quantities:
                call vor2vel

                !Calculate the source terms (sbuoys,svorts):
                call source

                !Update fields:
#ifdef ENABLE_BUOYANCY
                sbuoy = filt * (bsm + dt2 * sbuoys)
                !call field_combine_semi_spectral(sbuoy)
                !! omp parallel do private(iz)  default(shared)
                !do iz = 0, nz
                !    sbuoy(iz, :, :) = diss * sbuoy(iz, :, :)
                !enddo
                !! omp end parallel do
                !call field_decompose_semi_spectral(sbuoy)
#endif

                do nc = 1, 3
                    !$omp parallel workshare
                    svor(:, :, :, nc) = filt * (vortsm(:, :, :, nc) + dt2 * svorts(:, :, :, nc))
                    !$omp end parallel workshare
                    !call field_combine_semi_spectral(svor(:, :, :, nc))
                    !!omp parallel do private(iz)  default(shared)
                    !do iz = 0, nz
                    !    svor(iz, :, :, nc) = diss * svor(iz, :, :, nc)
                    !enddo
                    !!omp end parallel do
                    !call field_decompose_semi_spectral(svor(:, :, :, nc))
                enddo

                call adjust_vorticity_mean

            enddo

            t = t + dt

        end subroutine cn2_step

end module cn2
