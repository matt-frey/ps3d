module inversion_utils
    use constants
    use parameters, only : nx, ny, nz, dx, dxi, extent, ncelli, upper, lower
    use mpi_layout
    use mpi_environment
    use sta3dfft, only : initialise_fft &
                       , finalise_fft   &
                       , rkx            &
                       , rky            &
                       , rkz            &
                       , fftxyp2s       &
                       , fftxys2p       &
                       , fftsine        &
                       , fftcosine      &
                       , xfactors       &
                       , xtrig          &
                       , yfactors       &
                       , ytrig          &
                       , zfactors       &
                       , ztrig          &
                       , k2l2
    use stafft, only : dst
    use sta2dfft, only : ptospc
    use deriv1d, only : init_deriv
    use options, only : vor_visc, filtering
#ifdef ENABLE_BUOYANCY
    use options, only : buoy_visc
#endif
    use mpi_utils, only : mpi_print, mpi_stop
    implicit none

    private

    ! Ordering in physical space: z, y, x
    ! Ordering in spectral space: z, y, x

    ! Spectral dissipation operator for vorticity
    double precision, allocatable :: vhdis(:, :)

#ifdef ENABLE_BUOYANCY
    ! Spectral dissipation operator for buoyancy
    double precision, allocatable :: bhdis(:, :)
#endif


    double precision :: dzi, hdzi

    double precision :: vvisc
#ifdef ENABLE_BUOYANCY
    double precision :: bvisc
#endif

    logical :: is_initialised = .false.

    public :: init_inversion        &
            , init_diffusion        &
#ifdef ENABLE_BUOYANCY
            , bvisc                 &
            , bhdis                 &
#endif
            , vvisc                 &
            , hdzi                  &
            , vhdis                 &
            , call_ptospc

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! this routine is call by the genspec2d program
        subroutine call_ptospc(pp, ss)
            double precision, intent(inout) :: pp(box%lo(2):box%hi(2), box%lo(1):box%lo(2))
            double precision, intent(out)   :: ss(box%lo(2):box%hi(2), box%lo(1):box%lo(2))
            call ptospc(nx, ny, pp, ss, xfactors, yfactors, xtrig, ytrig)
        end subroutine call_ptospc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine init_diffusion(te, en)
            double precision, intent(in) :: te ! total energy
            double precision, intent(in) :: en ! enstrophy


            ! check if initialised
            if (.not. is_initialised) then
                call mpi_print("Error: Inversion not initialised!")
            endif


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

        subroutine init_inversion
            integer          :: kx, ky, kz

            if (is_initialised) then
                return
            endif

            is_initialised = .true.

            dzi = dxi(3)
            hdzi = f12 * dxi(3)

            call initialise_fft(extent)

        end subroutine init_inversion

end module inversion_utils
