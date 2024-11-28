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

    ! Spectral filter:
    double precision, allocatable :: filt(:, :, :)

    double precision :: dzi, hdzi

    double precision :: vvisc
#ifdef ENABLE_BUOYANCY
    double precision :: bvisc
#endif

    logical :: is_initialised = .false.

    public :: init_inversion        &
            , finalise_inversion    &
            , init_diffusion        &
#ifdef ENABLE_BUOYANCY
            , bvisc                 &
            , bhdis                 &
#endif
            , vvisc                 &
            , filt                  &
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

            !----------------------------------------------------------
            !Define de-aliasing filter:

            allocate(filt(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

            select case (filtering)
                case ("Hou & Li")
                    call init_hou_and_li_filter
                case ("2/3-rule")
                    call init_23rd_rule_filter
                case default
                    call init_hou_and_li_filter
            end select

            !Ensure filter does not change domain mean:
            if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
                filt(:, 0, 0) = one
            endif

        end subroutine init_inversion

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Define Hou and Li filter (2D and 3D):
        subroutine init_hou_and_li_filter
            integer          :: kx, ky, kz
            double precision :: kxmaxi, kymaxi, kzmaxi
            double precision :: skx(box%lo(1):box%hi(1)), &
                                sky(box%lo(2):box%hi(2)), &
                                skz(0:nz)

            call mpi_print("Using Hou & Li de-aliasing filter.")

            kxmaxi = one / maxval(rkx)
            skx = -36.d0 * (kxmaxi * rkx(box%lo(1):box%hi(1))) ** 36
            kymaxi = one/maxval(rky)
            sky = -36.d0 * (kymaxi * rky(box%lo(2):box%hi(2))) ** 36
            kzmaxi = one/maxval(rkz)
            skz = -36.d0 * (kzmaxi * rkz) ** 36

            do kx = box%lo(1), box%hi(1)
                do ky = box%lo(2), box%hi(2)
                    filt(0,  ky, kx) = dexp(skx(kx) + sky(ky))
                    filt(nz, ky, kx) = filt(0, ky, kx)
                    do kz = 1, nz-1
                        filt(kz, ky, kx) = filt(0, ky, kx) * dexp(skz(kz))
                    enddo
                enddo
            enddo

        end subroutine init_hou_and_li_filter

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Define de-aliasing filter (2/3 rule):
        subroutine init_23rd_rule_filter
            integer          :: kx, ky, kz
            double precision :: rkxmax, rkymax, rkzmax
            double precision :: skx(box%lo(1):box%hi(1)), &
                                sky(box%lo(2):box%hi(2)), &
                                skz(0:nz)

            call mpi_print("Using 2/3-rule de-aliasing filter.")

            rkxmax = maxval(rkx)
            rkymax = maxval(rky)
            rkzmax = maxval(rkz)

            do kx = box%lo(1), box%hi(1)
                if (rkx(kx) <= f23 * rkxmax) then
                    skx(kx) = one
                else
                    skx(kx) = zero
                endif
            enddo

            do ky = box%lo(2), box%hi(2)
                if (rky(ky) <= f23 * rkymax) then
                    sky(ky) = one
                else
                    sky(ky) = zero
                endif
            enddo

            do kz = 0, nz
                if (rkz(kz) <= f23 * rkzmax) then
                    skz(kz) = one
                else
                    skz(kz) = zero
                endif
            enddo

            ! Take product of 1d filters:
            do kx = box%lo(1), box%hi(1)
                do ky = box%lo(2), box%hi(2)
                  filt(0,  ky, kx) = skx(kx) * sky(ky)
                  filt(nz, ky, kx) = filt(0, ky, kx)
                  do kz = 1, nz-1
                     filt(kz, ky, kx) = filt(0, ky, kx) * skz(kz)
                  enddo
               enddo
            enddo
        end subroutine init_23rd_rule_filter

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine finalise_inversion
            deallocate(filt)

            call finalise_fft

        end subroutine finalise_inversion

end module inversion_utils
