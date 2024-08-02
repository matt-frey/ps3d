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
                       , ztrig
    use stafft, only : dst
    use sta2dfft, only : ptospc
    use deriv1d, only : init_deriv
    use options, only : viscosity
    use mpi_utils, only : mpi_print
    implicit none

    private

    ! Ordering in physical space: z, y, x
    ! Ordering in spectral space: z, y, x

    ! Note k2l2i = 1/(k^2+l^2) (except k = l = 0, then k2l2i(0, 0) = 0)
    double precision, allocatable :: k2l2i(:, :)

    ! Note k2l2 = k^2+l^2
    double precision, allocatable :: k2l2(:, :)

    double precision, allocatable :: green(:, :, :)

    ! Spectral dissipation operator
    double precision, allocatable :: hdis(:, :)

    ! Spectral filter:
    double precision, allocatable :: filt(:, :, :)

    logical :: is_initialised = .false.

    public :: init_inversion        &
            , finalise_inversion    &
            , init_diffusion        &
            , filt                  &
            , k2l2                  &
            , k2l2i                 &
            , hdis                  &
            , green                 &
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

        subroutine init_diffusion(bbdif, te, en)
            double precision, intent(in) :: bbdif ! (bbdif = max(b) - min(b) at t = 0):
            double precision, intent(in) :: te ! total energy
            double precision, intent(in) :: en ! enstrophy
            double precision             :: rkxmax, rkymax, K2max
            double precision             :: visc, wfac

            allocate(hdis(box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

            ! check if initialised
            if (.not. is_initialised) then
                call mpi_print("Error: Inversion not initialised!")
            endif

            rkxmax = maxval(rkx)
            rkymax = maxval(rky)

            !---------------------------------------------------------------------
            ! Damping, viscous or hyperviscous:
            if (viscosity%nnu .eq. 1) then
                !Define viscosity:
                visc = viscosity%prediss * sqrt(bbdif / rkxmax ** 3)
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

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine init_inversion
            integer          :: kx, ky, kz
            double precision :: kxmaxi, kymaxi, kzmaxi
            double precision :: skx(box%lo(1):box%hi(1)), &
                                sky(box%lo(2):box%hi(2)), &
                                skz(0:nz)

            if (is_initialised) then
                return
            endif

            is_initialised = .true.

            call initialise_fft(extent)

            !----------------------------------------------------------
            !Squared wavenumber array:
            allocate(k2l2i(box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
            allocate(k2l2(box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

            do kx = box%lo(1), box%hi(1)
                do ky = box%lo(2), box%hi(2)
                    k2l2(ky, kx) = rkx(kx) ** 2 + rky(ky) ** 2
                enddo
            enddo

            if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
                k2l2(0, 0) = one
            endif

            k2l2i = one / k2l2

            if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
                k2l2(0, 0) = zero
                k2l2i(0, 0) = zero
            endif

            !----------------------------------------------------------
            !Define Hou and Li filter (2D and 3D):
            allocate(filt(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

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

            !Ensure filter does not change domain mean:
            if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
                filt(:, 0, 0) = one
            endif

            !---------------------------------------------------------------------
            !Define Green function:
            allocate(green(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

            !$omp parallel do
            do kz = 1, nz
                green(kz, :, :) = - one / (k2l2 + rkz(kz) ** 2)
            enddo
            !$omp end parallel do
            green(0, :, :) = - k2l2i

        end subroutine init_inversion

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine finalise_inversion
            deallocate(green)
            deallocate(k2l2i)
            deallocate(k2l2)
            deallocate(filt)

            call finalise_fft

        end subroutine finalise_inversion

end module inversion_utils
