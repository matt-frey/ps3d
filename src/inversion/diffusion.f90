module diffusion
    use parameters, only : extent, nz
    use constants
    use mpi_layout
    use mpi_environment
    use sta3dfft, only : is_fft_initialised &
                       , rkx                &
                       , rky                &
                       , rkz                &
                       , k2l2
    use options, only : vor_visc
#ifdef ENABLE_BUOYANCY
    use options, only : buoy_visc
#endif
    use mpi_utils, only : mpi_print, mpi_stop
    implicit none

    private

    ! Ordering in physical space: z, y, x
    ! Ordering in spectral space: z, y, x

    ! Spectral dissipation operator for vorticity
    double precision, allocatable :: vhdis(:, :, :)

#ifdef ENABLE_BUOYANCY
    ! Spectral dissipation operator for buoyancy
    double precision, allocatable :: bhdis(:, :, :)
#endif

    double precision :: vvisc
#ifdef ENABLE_BUOYANCY
    double precision :: bvisc
#endif

    logical, protected :: is_diffusion_initialised = .false.

    public :: init_diffusion            &
            , is_diffusion_initialised  &
#ifdef ENABLE_BUOYANCY
            , bvisc                     &
            , bhdis                     &
#endif
            , vvisc                     &
            , vhdis

contains

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine init_diffusion(te, en)
        double precision, intent(in) :: te ! total energy
        double precision, intent(in) :: en ! enstrophy

        ! check if initialised
        if (is_diffusion_initialised) then
            return
        endif

        is_diffusion_initialised = .true.

        if (.not. is_fft_initialised) then
            call mpi_stop("Error: FFT not initialised.")
        endif

        vvisc = get_viscosity(vor_visc%length_scale, &
                              vor_visc%prediss,      &
                              vor_visc%nnu,          &
                              vor_visc%l_use_3d,     &
                              te, en)

        allocate(vhdis(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

        call init_dissipation('Vorticity', vvisc, vor_visc%nnu, vor_visc%l_use_3d, vhdis)

#ifdef ENABLE_BUOYANCY
        bvisc = get_viscosity(buoy_visc%length_scale, &
                              buoy_visc%prediss,      &
                              buoy_visc%nnu,          &
                              buoy_visc%l_use_3d,     &
                              te, en)

        allocate(0:nz, bhdis(box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

        call init_dissipation('Buoyancy', bvisc, buoy_visc%nnu, buoy_visc%l_use_3d, bhdis)
#endif

    end subroutine init_diffusion

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_viscosity(lscale, prediss, p, l_use_3d, te, en) result(vis)
        character(len=11), intent(in) :: lscale
        double precision,  intent(in) :: prediss
        integer,           intent(in) :: p
        logical,           intent(in) :: l_use_3d
        double precision,  intent(in) :: te ! total energy
        double precision,  intent(in) :: en ! enstrophy
        double precision              :: rkmsi, vis
        double precision              :: rkxmax, rkymax, rkzmax, K2max

        rkxmax = maxval(rkx)
        rkymax = maxval(rky)
        rkzmax = maxval(rkz)

        ! Define viscosity:
        K2max = max(rkxmax, rkymax) ** 2

        if (l_use_3d) then
            K2max = max(K2max, rkzmax ** 2)
        endif

        rkmsi = one / K2max

        select case (lscale)
            case ('Kolmogorov')
                vis = prediss *  (K2max * te /en) ** f13 * rkmsi ** p
            case ('geophysical')
                vis = prediss * rkmsi ** p
            case ('constant')
                vis = prediss
            case default
                call mpi_stop(&
                    "We only support 'Kolmogorov', 'geophysical' or 'constant'.")
        end select

    end function get_viscosity

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine init_dissipation(label, visc, p, l_use_3d, hdis)
        character(len=*), intent(in)  :: label
        double precision, intent(in)  :: visc
        integer,          intent(in)  :: p
        logical,          intent(in)  :: l_use_3d
        double precision, intent(out) :: hdis(0:nz,                &
                                              box%lo(2):box%hi(2), &
                                              box%lo(1):box%hi(1))
        integer                       :: kx, ky, kz

        !---------------------------------------------------------------------
        ! Damping, viscous or hyperviscous:
        if (p .eq. 1) then
            if (world%rank == world%root) then
                write(*,'(a,1p,e14.7)') label // ' molecular viscosity nu = ', visc
            endif
        else
            !Define hyperviscosity:
            if (world%rank == world%root) then
                write(*,'(a,1p,e14.7)') label // ' hyperviscosity nu = ', visc
            endif

        endif

        if (l_use_3d) then
            !Define 3D dissipation operator:
            do kx = box%lo(1), box%hi(1)
                do ky = box%lo(2), box%hi(2)
                    hdis(0, ky, kx) = visc * k2l2(ky, kx) ** p
                    hdis(nz, ky, kx) = hdis(0, ky, kx)
                    do kz = 1, nz-1
                        hdis(kz, ky, kx) = visc * (k2l2(ky, kx) + rkz(kz) ** 2) ** p
                    enddo
                enddo
            enddo
        else
            !Define 2D dissipation operator:
            do kx = box%lo(1), box%hi(1)
                do ky = box%lo(2), box%hi(2)
                    hdis(:, ky, kx) = visc * k2l2(ky, kx) ** p
                enddo
            enddo
        endif

    end subroutine init_dissipation

end module diffusion
