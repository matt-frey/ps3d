! =============================================================================
!           Low-storage 3rd and 4th order Runge-Kutta method
!            (see https://doi.org/10.5194/gmd-10-3145-2017)
! =============================================================================
module ls_rk
    use constants, only : zero
    use fields
    use inversion_utils
    use inversion_mod, only : vor2vel, source
    use field_diagnostics
    use mpi_utils, only : mpi_stop
    implicit none

    private

    integer, parameter :: dp=kind(zero)           ! double precision

    ! fourth order RK coefficients:
    double precision, dimension(5), target ::             &
        cas4 = (/- 567301805773.0_dp/1357537059087.0_dp,  &
                 -2404267990393.0_dp/2016746695238.0_dp,  &
                 -3550918686646.0_dp/2091501179385.0_dp,  &
                 -1275806237668.0_dp/842570457699.0_dp,   &
                 0.0/) !dummy value, not actually used

    double precision, dimension(5), target ::             &
        cbs4 =  (/1432997174477.0_dp/9575080441755.0_dp,  &
                  5161836677717.0_dp/13612068292357.0_dp, &
                  1720146321549.0_dp/2090206949498.0_dp,  &
                  3134564353537.0_dp/4481467310338.0_dp,  &
                  2277821191437.0_dp/14882151754819.0_dp/)

    ! thrird order RK coefficients:
    double precision, dimension(3), target :: &
        cas3 = (/  -5.0d0 / 9.0d0,            &
                 -153.0d0 / 128.d0,           &
                    0.0 /) !dummy value, not actually used

    double precision, dimension(3), target :: &
        cbs3 =  (/ 1.0d0 / 3.0d0,             &
                  15.0d0 / 16.0d0,            &
                   8.0d0 / 15.0d0 /)

    double precision, dimension(:), pointer :: captr, cbptr

    integer :: n_stages

    public :: ls_rk_setup, ls_rk_step

    contains

        subroutine ls_rk_setup(order)
            integer, intent(in) :: order

            select case (order)
                case (3)
                    captr => cas3
                    cbptr => cbs3
                    n_stages = 3
                case (4)
                    captr => cas4
                    cbptr => cbs4
                    n_stages = 5
                case default
                    call mpi_stop('Only third and fourth order RK supported.')
            end select

        end subroutine ls_rk_setup

        ! Advances the parcels by a single ls-RK-4 step. It calls a
        ! function to obtain the current time step based on the velocity
        ! strain and the buoyancy gradient.
        ! @param[in] t is the time
        ! Precondition: this routine assumes that the fields are
        ! up-to-date for the first sub-step
        subroutine ls_rk_step(t, dt)
            double precision, intent(inout) :: t
            double precision, intent(in)    :: dt
            integer                         :: n

            do n = 1, n_stages-1
                call ls_rk_substep(dt, n)
            enddo
            call ls_rk_substep(dt, n_stages)

            t = t + dt
        end subroutine ls_rk_step

        ! Do a ls-RK-4 substep.
        ! @param[in] dt is the time step
        ! @param[in] step is the number of the substep (1 to 5 or 1 to 3)
        subroutine ls_rk_substep(dt, step)
            double precision, intent(in) :: dt
            integer,          intent(in) :: step
            double precision             :: ca, cb
            integer                      :: nc

            ca = captr(step)
            cb = cbptr(step)

            if (step == 1) then
                !$omp parallel workshare
#ifdef ENABLE_BUOYANCY
                bsm    = sbuoys
#endif
                vortsm = svorts
                !$omp end parallel workshare
            else
                call vor2vel

                call source

                !$omp parallel workshare
#ifdef ENABLE_BUOYANCY
                bsm = bsm + sbuoys
#endif
                vortsm = vortsm + svorts
                !$omp end parallel workshare
            endif

#ifdef ENABLE_BUOYANCY
            sbuoy = filt * (sbuoy + cb * dt * bsm)
#endif
            do nc = 1, 3
                !$omp parallel workshare
                svor(:, :, :, nc) = filt * (svor(:, :, :, nc) + cb * dt * vortsm(:, :, :, nc))
                !$omp end parallel workshare
            enddo

            call adjust_vorticity_mean

            if (step == n_stages) then
               return
            endif

            !$omp parallel workshare
#ifdef ENABLE_BUOYANCY
            bsm = ca * bsm
#endif
            vortsm = ca * vortsm
            !$omp end parallel workshare

        end subroutine ls_rk_substep

end module ls_rk
