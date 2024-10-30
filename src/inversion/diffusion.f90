module diffusion
    use mpi_environment, only : world
    use parameters, only : nz, extent, lower, upper
    use constants, only : zero, f13, f12, one
    use mpi_layout, only : box
    use sta3dfft, only : initialise_fft &
                       , finalise_fft   &
                       , rkx            &
                       , rky            &
                       , k2l2
    use zops, only : init_zops          &
                   , finalise_zops      &
                   , zg                 &
                   , d1z                &
                   , zcheb
    use options, only : vor_visc
#ifdef ENABLE_BUOYANCY
    use options, only : buoy_visc
#endif
    use mpi_utils, only : mpi_stop
    implicit none

    private

    logical :: l_hdis_initialised = .false.

    ! Spectral dissipation operator for vorticity
    double precision, allocatable :: vhdis(:, :)

#ifdef ENABLE_BUOYANCY
    ! Spectral dissipation operator for buoyancy
    double precision, allocatable :: bhdis(:, :)
#endif

    double precision :: vvisc
#ifdef ENABLE_BUOYANCY
    double precision :: bvisc
#endif

    ! Vertically varying viscosity
    double precision, allocatable :: visc(:)

    public :: init_diffusion        &
            , finalise_diffusion    &
            , apply_zdiffusion      &
            , vhdis                 &
#ifdef ENABLE_BUOYANCY
            , bhdis                 &
            , bvisc                 &
#endif
            , vvisc                 &
            , visc

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine init_diffusion(te, en)
            double precision, intent(in) :: te ! total energy
            double precision, intent(in) :: en ! enstrophy
            double precision             :: rkxmax, rkymax, K2max
            double precision             :: wfac, delta, hvisc
            integer                      :: iz

            if (l_hdis_initialised) then
                return
            endif

            l_hdis_initialised = .true.

            allocate(visc(0:nz))

            !------------------------------------------------------------------
            ! Ensure required modules are initialised:

            ! Ensure FFT is initialised
            ! (this won't do anything if already initialised)
            call initialise_fft(extent)

            ! Ensure Chebyshev is initialised:
            call init_zops

            rkxmax = maxval(rkx)
            rkymax = maxval(rky)

            !------------------------------------------------------------------
            ! Define viscosity:
            K2max = max(rkxmax, rkymax) ** 2
            wfac = one / K2max
            hvisc = vor_visc%prediss *  (K2max * te /en) ** f13 * wfac ** vor_visc%nnu

            !------------------------------------------------------------------
            ! N(z) = tanh((z-z_min)/delta) * tanh((z_max-z)/delta)
            ! 0 <= N(z) <= 1
            ! N vanishes on each boundary, while reaching 1 at the domain centre (in z)


            ! specified boundary layer width
            ! delta = (upper(3) - lower(3))/dble((2*nz))
            delta = (upper(3) - lower(3)) / dble((nz))
            do iz = 0, nz
                ! N(z)
                visc(iz) = (tanh(5.0d0 * (zcheb(iz) + 1.0d0)) &
                          * tanh(5.0d0 * (1.0d0 - zcheb(iz)))) ** 2
!                 visc(iz) = tanh((zg(iz) - lower(3)) / delta)   &
!                          * tanh((upper(3) - zg(iz)) / delta)

                ! nu(z) = nu_max * N(z)
                visc(iz) = hvisc * visc(iz)
            enddo

            vvisc = get_viscosity(vor_visc%length_scale, &
                                  vor_visc%prediss,      &
                                  vor_visc%nnu, te, en)

            allocate(vhdis(box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

            call init_dissipation('Vorticity', vvisc, vor_visc%nnu, vhdis)

#ifdef ENABLE_BUOYANCY
            bvisc = get_viscosity(buoy_visc%length_scale, &
                                  buoy_visc%prediss,      &
                                  buoy_visc%nnu, te, en)

            allocate(bhdis(box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

            call init_dissipation('Buoyancy', bvisc, buoy_visc%nnu, bhdis)
#endif

        end subroutine init_diffusion

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        function get_viscosity(lscale, prediss, p, te, en) result(vis)
            character(len=11), intent(in) :: lscale
            double precision,  intent(in) :: prediss
            integer,           intent(in) :: p
            double precision,  intent(in) :: te ! total energy
            double precision,  intent(in) :: en ! enstrophy
            double precision              :: rkmsi, vis
            double precision              :: rkxmax, rkymax, K2max

            rkxmax = maxval(rkx)
            rkymax = maxval(rky)


            ! Define viscosity:
            K2max = max(rkxmax, rkymax) ** 2
            rkmsi = one / K2max

            select case (lscale)
                case ('Kolmogorov')
                    vis = prediss *  (K2max * te /en) ** f13 * rkmsi ** p
                case ('geophysical')
                    vis = prediss * rkmsi ** p
                case default
                    call mpi_stop(&
                        "We only support 'Kolmogorov' or 'geophysical'")
            end select

        end function get_viscosity

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine init_dissipation(label, visc, p, hdis)
            character(len=*), intent(in)  :: label
            double precision, intent(in)  :: visc
            integer,          intent(in)  :: p
            double precision, intent(out) :: hdis(box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1))

            !---------------------------------------------------------------------
            ! Damping, viscous or hyperviscous:
            if (p .eq. 1) then
                if (world%rank == world%root) then
                    write(*,'(a,1p,e14.7)') label // ' molecular viscosity nu = ', visc
                endif

                !Define spectral dissipation operator:
                !$omp parallel workshare
                hdis = visc * k2l2
                !$omp end parallel workshare
             else
                !Define hyperviscosity:
                if (world%rank == world%root) then
                    write(*,'(a,1p,e14.7)') label // ' hyperviscosity nu = ', visc
                endif

                !Define dissipation operator:
                !$omp parallel workshare
                hdis = visc * k2l2 ** p
                !$omp end parallel workshare
             endif

        end subroutine init_dissipation

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
                deallocate(vhdis)
#ifdef ENABLE_BUOYANCY
                deallocate(bhdis)
#endif
                deallocate(visc)
                l_hdis_initialised = .false.
            endif

            call finalise_fft
            call finalise_zops

        end subroutine finalise_diffusion

end module diffusion
