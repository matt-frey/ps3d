module diffusion
    use mpi_environment, only : world
    use options, only : viscosity
    use parameters, only : nz, extent, lower, upper
    use constants, only : zero, f13, f12, one
    use mpi_layout, only : box
    use sta3dfft, only : initialise_fft &
                       , rkx            &
                       , rky            &
                       , k2l2
    use zops, only : init_zops          &
                   , finalise_zops      &
                   , zg                 &
                   , d1z
    implicit none

    private

    logical :: l_hdis_initialised = .false.

    ! Spectral dissipation operator
    double precision, allocatable :: hdis(:, :, :)

    ! Vertically varying viscosity
    double precision, allocatable :: visc(:)

    public :: init_diffusion    &
            , hdis              &
            , visc

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine init_diffusion(bbdif, te, en)
            double precision, intent(in) :: bbdif ! (bbdif = max(b) - min(b) at t = 0):
            double precision, intent(in) :: te ! total energy
            double precision, intent(in) :: en ! enstrophy
            double precision             :: rkxmax, rkymax, K2max
            double precision             :: wfac, delta
            integer                      :: iz

            if (l_hdis_initialised) then
                return
            endif

            l_hdis_initialised = .true.

            allocate(hdis(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
            allocate(visc(0:nz))

            !------------------------------------------------------------------
            ! Ensure required modules are initialised:

            ! Ensure FFT is initialised
            ! (this won't do anything if already initialised)
            call initialise_fft(extent)

            ! Ensure Chebyshev is initialised:
            call init_zops

            !------------------------------------------------------------------
            ! N(z) = tanh((z-z_min)/delta) * tanh((z_max-z)/delta)
            ! 0 <= N(z) <= 1
            ! N vanishes on each boundary, while reaching 1 at the domain centre (in z)


            ! specified boundary layer width
            ! delta = (upper(3) - lower(3))/dble((2*nz))
            delta = (upper(3) - lower(3)) / dble((nz))
            do iz = 0, nz
                ! N(z)
                visc(iz) = tanh((zg(iz) - lower(3)) / delta)   &
                         * tanh((upper(3) - zg(iz)) / delta)

                ! nu(z) = nu_max * N(z)
                visc(iz) = viscosity%prediss * visc(iz)
            enddo

            rkxmax = maxval(rkx)
            rkymax = maxval(rky)

            !---------------------------------------------------------------------
            ! Damping, viscous or hyperviscous:
            if (viscosity%nnu .eq. 1) then
                !Define viscosity:
                if (bbdif > zero) then
                    visc = visc * sqrt(bbdif / rkxmax ** 3)
                endif

                if (world%rank == world%root) then
                    write(*,'(a,1p,e14.7)') ' Max viscosity nu = ', maxval(visc)
                endif

                !Define spectral dissipation operator:
                !$omp parallel do
                do iz = 0, nz
                    hdis(iz, :, :) = visc(iz) * k2l2
                enddo
                !$omp end parallel do
             else
                !Define hyperviscosity:
                K2max = max(rkxmax, rkymax) ** 2
                wfac = one / K2max
                visc = visc *  (K2max * te /en) ** f13
                if (world%rank == world%root) then
                    write(*,'(a,1p,e14.7)') ' Max hyperviscosity nu = ', maxval(visc) * wfac ** viscosity%nnu
                endif

                !Define dissipation operator:
                !$omp parallel do
                do iz = 0, nz
                    hdis(iz, :, :) = visc(iz) * (wfac * k2l2) ** viscosity%nnu
                enddo
                !$omp end parallel do

                if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
                    !Ensure average is not modified by hyperviscosity:
                    hdis(:, 0, 0) = zero
                endif
             endif
        end subroutine init_diffusion

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine apply_zdiffusion(fs, dt)
            double precision, intent(inout) :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision, intent(in)    :: dt
            double precision                :: u(0:nz), Am(0:nz, 0:nz), Lm(0:nz, 0:nz), eye(0:nz, 0:nz)
            integer                         :: ipiv(0:nz)
            integer                         :: kx, ky, info, iz

            eye = zero
            Am = zero
            do iz = 0, nz
                eye(iz, iz) = one
                Am(iz, iz) = visc(iz)
            enddo

            Am = matmul(matmul(d1z, Am), d1z)
            Lm = eye - f12 * dt * Am

            do kx = box%lo(1), box%hi(1)
                do ky = box%hi(2), box%hi(2)

                    ! right-hand side u^{n}
                    u = fs(:, ky, kx)

                    ! solve the linear system to get \bar{u}
                    Am = Lm
                    ipiv = 0
                    call dgesv(nz+1, 1, Am, nz+1, ipiv, u, nz+1, info)

                    ! update u^{n+1} = 2\bar{u} - u^{n}
                    fs(:, ky, kx) = 2.0d0 * u - fs(:, ky, kx)
                enddo
            enddo

        end subroutine apply_zdiffusion

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine finalise_diffusion

            if (l_hdis_initialised) then
                deallocate(hdis)
                deallocate(visc)
                l_hdis_initialised = .false.
            endif

            call finalise_fft
            call finalise_zops

        end subroutine finalise_diffusion

end module diffusion
