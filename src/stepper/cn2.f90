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
    use model, only : layout
    use options, only : vor_visc
#ifdef ENABLE_BUOYANCY
    use options, only : buoy_visc
#endif
    use stepper_mod, only : stepper_t
    use constants, only : f12, one
    use parameters, only : nz
    use fields
    use diffusion
    use inversion_mod, only : vor2vel, source
    use field_diagnostics
    implicit none

    type, extends(stepper_t) :: cn2
        contains
            procedure :: set_diffusion => cn2_set_diffusion
            procedure :: setup  => cn2_setup
            procedure :: step => cn2_step
    end type

    ! Number of iterations of above scheme:
    integer, parameter:: niter = 2

    double precision :: dt2

contains

    subroutine cn2_set_diffusion(self, dt, vorch, bf)
        class(cn2),       intent(inout) :: self
        double precision, intent(in)    :: dt
        double precision, intent(in)    :: vorch, bf
        double precision                :: dfac

        !---------------------------------------------------------------------
        if (vor_visc%nnu .eq. 1) then
            !Update diffusion operator used in time stepping:
            dfac = dt
        else
            !Update hyperdiffusion operator used in time stepping:
            dfac = vorch * dt
        endif

        !$omp parallel workshare
        vdop = one / (one + dfac * vdiss)
        !$omp end parallel workshare

#ifdef ENABLE_BUOYANCY
        if (buoy_visc%nnu .eq. 1) then
            !Update diffusion operator used in time stepping:
            dfac = dt
        else
            !Update hyperdiffusion operator used in time stepping:
            dfac = bf * dt
        endif

        !$omp parallel workshare
        bdop = one / (one + dfac * bdiss)
        !$omp end parallel workshare
#endif
        !(see inversion_utils.f90)

    end subroutine cn2_set_diffusion

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
        integer                         :: iz

        !Update value of dt/2:
        dt2 = f12 * dt

        !------------------------------------------------------------------
        !Start with a guess for F^{n+1} for all fields:

        !Initialise iteration (dt = dt/2 below):
#ifdef ENABLE_BUOYANCY
        !$omp parallel workshare
        bsm = sbuoy + dt2 * sbuoys
        sbuoy = bdop * (bsm + dt2 * sbuoys)
        !$omp end parallel workshare
#endif

        ! Advance interior and boundary values of vorticity
        !$omp parallel workshare
        vortsm = svor + dt2 * svorts
        !$omp end parallel workshare

        !$omp parallel do private(nc)  default(shared)
        do nc = 1, 3
            svor(:, :, :, nc) = vdop * (vortsm(:, :, :, nc) + &
                                  dt2 * svorts(:, :, :, nc))
        enddo
        !$omp end parallel workshare

        ! Ensure zero global mean horizontal vorticity conservation:
        do nc = 1, 2
            call layout%adjust_semi_spectral_mean(svor(:, :, :, nc), ini_vor_mean(nc))
        enddo

        !bdop and vdop are the full (3d) hyperdiffusive operators
        !(see inversion/diffusion.f90)

        !------------------------------------------------------------------
        !Iterate to improve estimates of F^{n+1}:
        do iter = 1, niter
            !Perform inversion at t^{n+1} from estimated quantities:
            call vor2vel

            !Calculate the source terms (sbuoys,svorts):
            call source

            !Update fields:
#ifdef ENABLE_BUOYANCY
            !$omp parallel workshare
            sbuoy = bdop * (bsm + dt2 * sbuoys)
            !$omp end parallel workshare
#endif

            !$omp parallel do private(nc)  default(shared)
            do nc = 1, 3
                svor(:, :, :, nc) = vdop * (vortsm(:, :, :, nc) + &
                                      dt2 * svorts(:, :, :, nc))
            enddo
            !$omp end parallel workshare

            ! Ensure zero global mean horizontal vorticity conservation:
            do nc = 1, 2
                call layout%adjust_semi_spectral_mean(svor(:, :, :, nc), ini_vor_mean(nc))
            enddo

            !bdop and vdop are the full (3d) hyperdiffusive operators
            !(see inversion/diffusion.f90)

        enddo

        t = t + dt

    end subroutine cn2_step

end module cn2_mod
