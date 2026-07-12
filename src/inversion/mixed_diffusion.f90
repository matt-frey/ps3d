module diffusion
    use parameters, only : extent, vhr2, nz
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

    ! Spectral dissipation operator for vorticity (3d):
    double precision, allocatable :: vdiss(:, :, :)

#ifdef ENABLE_BUOYANCY
    ! Spectral dissipation operator for buoyancy (3d):
    double precision, allocatable :: bdiss(:, :, :)
#endif

    logical, protected :: is_diffusion_initialised = .false.

    public :: init_diffusion            &
            , is_diffusion_initialised  &
#ifdef ENABLE_BUOYANCY
            , bdiss                     &
#endif
            , vdiss

contains

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine init_diffusion(te, en)
        double precision, intent(in) :: te ! total energy
        double precision, intent(in) :: en ! enstrophy
        ! Horizontal and vertical (hyper)viscosity coefficients:
        double precision             :: vis(2)

        ! check if initialised
        if (is_diffusion_initialised) then
            return
        endif

        is_diffusion_initialised = .true.

        if (.not. is_fft_initialised) then
            call mpi_stop("Error: FFT not initialised.")
        endif

        vis = get_viscosity(vor_visc%length_scale, &
                            vor_visc%prediss,      &
                            vor_visc%vweight,      &
                            vor_visc%nnu, te, en)

        allocate(vdiss(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

        call init_dissipation('Vorticity', vis, vor_visc%vweight, &
                              vor_visc%nnu, vdiss)

#ifdef ENABLE_BUOYANCY
        vis = get_viscosity(buoy_visc%length_scale, &
                            buoy_visc%prediss,      &
                            buoy_visc%vweight,      &
                            buoy_visc%nnu, te, en)

        allocate(bdiss(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

        call init_dissipation('Buoyancy', vis, buoy_visc%vweight, &
                              buoy_visc%nnu, bdiss)
#endif

    end subroutine init_diffusion

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_viscosity(lscale, prediss, vwt, p, te, en) result(vis)
        character(len=11), intent(in) :: lscale
        double precision,  intent(in) :: prediss
        double precision,  intent(in) :: vwt
        integer,           intent(in) :: p
        double precision,  intent(in) :: te ! total energy
        double precision,  intent(in) :: en ! enstrophy
        double precision              :: rkmsi, hvis, zvis
        double precision              :: rkxmax, rkymax, K2max
        double precision              :: vis(2)

        rkxmax = maxval(rkx)
        rkymax = maxval(rky)

        ! Define horizontal viscosity (hvis):
        K2max = max(rkxmax, rkymax) ** 2
        rkmsi = one / K2max

        select case (lscale)
            case ('Kolmogorov')
                hvis = prediss *  (K2max * te /en) ** f13 * rkmsi ** p
            case ('geophysical')
                hvis = prediss * rkmsi ** p
            case ('constant')
                hvis = prediss
            case default
                call mpi_stop(&
                 "We only support 'Kolmogorov', 'geophysical' or 'constant'.")
        end select

        ! Define vertical viscosity (zvis):
        zvis = hvis * vwt * vhr2 * rkmsi * K2max ** p
        ! Here vhr2 = L_z^2/(L_x*L_y) & rkmsi * K2max ** p = K2max ** (p-1)

        vis(1) = hvis
        vis(2) = zvis

    end function get_viscosity

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine init_dissipation(label, vis, vwt, p, fdis)
        character(len=*), intent(in)  :: label
        double precision, intent(in)  :: vis(2)
        double precision, intent(in)  :: vwt
        integer,          intent(in)  :: p
        double precision, intent(out) :: fdis(0:nz, box%lo(2):box%hi(2), &
                                                    box%lo(1):box%hi(1))
        double precision              :: hdis(box%lo(2):box%hi(2), &
                                              box%lo(1):box%hi(1))
        double precision              :: skz2(nz-1), hvisc
        integer                       :: m

        !---------------------------------------------------------------------
        ! Extract horizontal (hyper)viscosity coefficient:
        hvisc = vis(1)

        ! Scaled squared vertical wavenumber:
        skz2 = vwt * vhr2 * rkz(1:nz-1)**2
        ! vhr2 = L_z^2/(L_x*L_y) is set in parameters.f90

        ! Damping, viscous or hyperviscous:
        if (p .eq. 1) then
           !Molecular viscosity:

           if (world%rank == world%root) then
              write(*,'(a,1p,e14.7)') label // &
                   ' horizontal molecular viscosity nu_h = ', hvisc
           endif
           !Define 2d spectral dissipation operator (hdis):
           !$omp parallel workshare
           hdis = hvisc * k2l2
           fdis( 0, :, :) = hdis
           fdis(nz, :, :) = hdis
           !$omp end parallel workshare

           !$omp parallel do private(m)  default(shared)
           do m = 1, nz-1
              fdis(m, :, :) = hvisc * (k2l2 + skz2(m))
           enddo
           !$omp end parallel do

        else
           !Hyperviscosity:
           
           if (world%rank == world%root) then
              write(*,'(a,1p,e14.7)') label // &
                   ' horizontal hyperviscosity nu_h = ', hvisc
           endif
           !Define 2d spectral dissipation operator (hdis):
           !$omp parallel workshare
           hdis = hvisc * k2l2 ** p
           fdis( 0, :, :) = hdis
           fdis(nz, :, :) = hdis
           !$omp end parallel workshare

           !$omp parallel do private(m)  default(shared)
           do m = 1, nz-1
              fdis(m, :, :) = hvisc * (k2l2 + skz2(m)) ** p
           enddo
           !$omp end parallel do
        endif

    end subroutine init_dissipation

end module diffusion
