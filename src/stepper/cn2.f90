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
module cn2_mod
#ifndef ENABLE_SMAGORINSKY
    use options, only : viscosity
    use diffusion, only : hdis
#endif
    use advance_mod, only : base_stepper
    use constants, only : f12, one
    use parameters, only : nz
    use fields
    use inversion_utils
    use inversion_mod, only : vor2vel, source
    use field_diagnostics
    implicit none

    type, extends(base_stepper) :: cn2
        contains
#ifndef ENABLE_SMAGORINSKY
            procedure :: set_diffusion => cn2_set_diffusion
#endif
            procedure :: setup  => cn2_setup
            procedure :: step => cn2_step
    end type

    ! Number of iterations of above scheme:
    integer, parameter:: niter = 2

    double precision :: dt2

    contains

#ifndef ENABLE_SMAGORINSKY
        subroutine cn2_set_diffusion(self, dt, vorch)
            class(cn2),       intent(inout) :: self
            double precision, intent(in)    :: dt
            double precision, intent(in)    :: vorch
            double precision                :: dfac

            !---------------------------------------------------------------------
            if (viscosity%nnu .eq. 1) then
                !Update diffusion operator used in time stepping:
                dfac = dt
            else
                !Update hyperdiffusion operator used in time stepping:
                dfac = vorch * dt
             endif

            !$omp parallel workshare
            diss = one / (one + dfac * hdis)
            !$omp end parallel workshare
            !(see inversion_utils.f90)

        end subroutine cn2_set_diffusion
#endif

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine cn2_setup(self)
            class(cn2), intent(inout) :: self

            ! do nothing here

        end subroutine cn2_setup

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine cn2_step(self, t, dt)
            class(cn2),       intent(inout) :: self
            double precision, intent(inout) :: t
            double precision, intent(in)    :: dt
            integer                         :: iter
            integer                         :: nc
#ifndef ENABLE_SMAGORINSKY
            integer                         :: iz
#endif

            !Update value of dt/2:
            dt2 = f12 * dt

            !------------------------------------------------------------------
            !Start with a guess for F^{n+1} for all fields:

            !Initialise iteration (dt = dt/2 below):
#ifdef ENABLE_BUOYANCY
            bsm = sbuoy + dt2 * sbuoys
            !$omp parallel do private(iz)  default(shared)
            for iz = 0, nz
                sbuoy(iz, :, :) = filt * (bsm(iz, :, :) + dt2 * sbuoys(iz, :, :))
#ifndef ENABLE_SMAGORINSKY
                sbuoy(iz, :, :) = diss(iz, :, :) * sbuoy(iz, :, :)
#endif
            enddo
            !$omp end parallel do
#endif

            ! Advance interior and boundary values of vorticity
            !$omp parallel workshare
            vortsm = svor + dt2 * svorts
            !$omp end parallel workshare

            !$omp parallel do collapse(2) private(iz) default(shared)
            do nc = 1, 3
                do iz = 0, nz
                    svor(iz, :, :, nc) = filt * (vortsm(iz, :, :, nc) + dt2 * svorts(iz, :, :, nc))
#ifndef ENABLE_SMAGORINSKY
                    svor(iz, :, :, nc) = diss(iz, :, :) * svor(iz, :, :, nc)
#endif
                enddo
            enddo
            !$omp end parallel do

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
                !$omp parallel do private(iz)  default(shared)
                do iz = 0, nz
                    sbuoy(iz, :, :) = filt * (bsm(iz, :, :) + dt2 * sbuoys(iz, :, :))
#ifndef ENABLE_SMAGORINSKY
                    sbuoy(iz, :, :) = diss(iz, :, :) * sbuoy(iz, :, :)
#endif
                enddo
                !$omp end parallel do
#endif

                !$omp parallel do collapse(2) private(iz) default(shared)
                do nc = 1, 3
                    do iz = 0, nz
                        svor(iz, :, :, nc) = filt * (vortsm(iz, :, :, nc) + dt2 * svorts(iz, :, :, nc))
#ifndef ENABLE_SMAGORINSKY
                        svor(iz, :, :, nc) = diss(iz, :, :) * svor(iz, :, :, nc)
#endif
                    enddo
                enddo
                !$omp end parallel do

                call adjust_vorticity_mean

            enddo

            t = t + dt

        end subroutine cn2_step

end module cn2_mod
