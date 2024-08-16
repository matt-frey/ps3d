module diffusion
    use mpi_environment, only : world
    use options, only : viscosity
    use parameters, only : extent
    use constants, only : zero, f13, one
    use mpi_layout, only : box
    use sta3dfft, only : initialise_fft &
                       , rkx            &
                       , rky            &
                       , k2l2
    implicit none

    private

    ! Spectral dissipation operator
    double precision, allocatable :: hdis(:, :)

    public :: init_diffusion        &
            , hdis

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine init_diffusion(bbdif, te, en)
            double precision, intent(in) :: bbdif ! (bbdif = max(b) - min(b) at t = 0):
            double precision, intent(in) :: te ! total energy
            double precision, intent(in) :: en ! enstrophy
            double precision             :: rkxmax, rkymax, K2max
            double precision             :: visc, wfac

            allocate(hdis(box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

            ! ensure FFT is initialised
            ! (this won't do anything if already initialised)
            call initialise_fft(extent)

            rkxmax = maxval(rkx)
            rkymax = maxval(rky)

            !---------------------------------------------------------------------
            ! Damping, viscous or hyperviscous:
            if (viscosity%nnu .eq. 1) then
                !Define viscosity:
                visc = viscosity%prediss
                if (bbdif > zero) then
                    visc = visc * sqrt(bbdif / rkxmax ** 3)
                endif

                if (world%rank == world%root) then
                    write(*,'(a,1p,e14.7)') ' Viscosity nu = ', visc
                endif

                !Define spectral dissipation operator:
                !$omp parallel workshare
                hdis = visc * k2l2
                !$omp end parallel workshare
             else
                !Define hyperviscosity:
                K2max = max(rkxmax, rkymax) ** 2
                wfac = one / K2max
                visc = viscosity%prediss *  (K2max * te /en) ** f13
                if (world%rank == world%root) then
                    write(*,'(a,1p,e14.7)') ' Hyperviscosity nu = ', visc * wfac ** viscosity%nnu
                endif

                !Define dissipation operator:
                !$omp parallel workshare
                hdis = visc * (wfac * k2l2) ** viscosity%nnu
                !$omp end parallel workshare

                if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
                    !Ensure average is not modified by hyperviscosity:
                    hdis(0, 0) = zero
                endif
             endif
        end subroutine init_diffusion

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine finalise_diffusion
            if (allocated(hdis)) then
                deallocate(hdis)
            endif
        end subroutine finalise_diffusion

end module diffusion
