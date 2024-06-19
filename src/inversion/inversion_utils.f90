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

    double precision, allocatable :: gamtop(:), gambot(:)

    ! See for definitions in
    ! Dritschel D, Frey M. The stability of inviscid Beltrami flow between parallel free-slip impermeable
    ! boundaries. Journal of Fluid Mechanics. 2023;954:A31. doi:10.1017/jfm.2022.1007
    double precision, allocatable :: thetam(:, :, :)    ! theta_{-}         (eq. 3.10)
    double precision, allocatable :: thetap(:, :, :)    ! theta_{+}         (eq. 3.11)
    double precision, allocatable :: dthetam(:, :, :)   ! dtheta_{-}/dz
    double precision, allocatable :: dthetap(:, :, :)   ! dtheta_{+}/dz
    double precision, allocatable :: phim(:, :, :)      ! phi_{-}           (eq. 3.4a)
    double precision, allocatable :: phip(:, :, :)      ! phi_{+}           (eq. 3.4b)
#ifdef ENABLE_BUOYANCY
    double precision, allocatable :: dphim(:, :, :)     ! dphi_{-}/dz
    double precision, allocatable :: dphip(:, :, :)     ! dphi_{+}/dz
#endif

    double precision :: dzi, hdzi

    logical :: is_initialised = .false.

    public :: init_inversion        &
            , finalise_inversion    &
            , init_diffusion        &
#ifdef ENABLE_BUOYANCY
            , diffz                 &
            , dphim                 &
            , dphip                 &
#endif
            , central_diffz         &
            , filt                  &
            , hdzi                  &
            , k2l2                  &
            , k2l2i                 &
            , hdis                  &
            , green                 &
            , thetap                &
            , thetam                &
            , dthetap               &
            , dthetam               &
            , gambot                &
            , gamtop                &
            , call_ptospc

    public :: field_combine_semi_spectral   &
            , field_combine_physical        &
            , field_decompose_semi_spectral &
            , field_decompose_physical

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
            integer          :: kx, ky, iz, kz
            double precision :: z, zm(0:nz), zp(0:nz)
            double precision :: phip00(0:nz)
            double precision :: rkxmax, rkymax, rkzmax
            double precision :: skx(box%lo(1):box%hi(1)), &
                                sky(box%lo(2):box%hi(2)), &
                                skz(0:nz)

            if (is_initialised) then
                return
            endif

            is_initialised = .true.

            dzi = dxi(3)
            hdzi = f12 * dxi(3)

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
            !Define de-aliasing filter (2/3 rule):
            allocate(filt(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

            rkxmax = maxval(rkx)
            rkymax = maxval(rky)
            rkzmax = maxval(rkz)

            do kx = box%lo(1), box%hi(1)
                if (rkx(kx) <= f23 * rkxmax)
                    skx(kx) = one
                else
                    skx(kx) = zero
                endif
            enddo

            do ky = box%lo(2), box%hi(2)
                if (rky(ky) <= f23 * rkymax)
                    sky(ky) = one
                else
                    sky(ky) = zero
                endif
            enddo

            do kz = 0, nz
                if (rkz(kz) <= f23 * rkzmax)
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

            !---------------------------------------------------------------------
            !Define zm = zmax - z, zp = z - zmin
            !$omp parallel do private(z)
            do iz = 0, nz
                z = lower(3) + dx(3) * dble(iz)
                zm(iz) = upper(3) - z
                zp(iz) = z - lower(3)
            enddo
            !$omp end parallel do

            !---------------------------------------------------------------------
            !Hyperbolic functions used for solutions of Laplace's equation:
            allocate(phim(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
            allocate(phip(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
#ifdef ENABLE_BUOYANCY
            allocate(dphim(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
            allocate(dphip(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
#endif
            allocate(thetam(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
            allocate(thetap(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
            allocate(dthetam(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
            allocate(dthetap(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

            do kx = box%lo(1), box%hi(1)
                do ky = max(1, box%lo(2)), box%hi(2)
                    call set_hyperbolic_functions(kx, ky, zm, zp)
                enddo
            enddo

            ! ky = 0
            if (box%lo(2) == 0) then
                do kx = max(1, box%lo(1)), box%hi(1)
                    call set_hyperbolic_functions(kx, 0, zm, zp)
                enddo
            endif

            phip00 = zero
            if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
                !$omp parallel workshare
                ! kx = ky = 0
                phim(:, 0, 0) = zm / extent(3)
                phip(:, 0, 0) = zp / extent(3)

#ifdef ENABLE_BUOYANCY
                dphim(:, 0, 0) = - one / extent(3)
                dphip(:, 0, 0) =   one / extent(3)
#endif

                thetam(:, 0, 0) = zero
                thetap(:, 0, 0) = zero

                dthetam(:, 0, 0) = zero
                dthetap(:, 0, 0) = zero

                phip00 = phip(:, 0, 0)
                !$omp end parallel workshare
            endif

            !---------------------------------------------------------------------
            !Define gamtop as the integral of phip(iz, 0, 0) with zero average:
            allocate(gamtop(0:nz))
            allocate(gambot(0:nz))

            call MPI_Allreduce(MPI_IN_PLACE,            &
                               phip00(0:nz),            &
                               nz+1,                    &
                               MPI_DOUBLE_PRECISION,    &
                               MPI_SUM,                 &
                               world%comm,              &
                               world%err)

            !$omp parallel workshare
            gamtop = f12 * extent(3) * (phip00 ** 2 - f13)
            !$omp end parallel workshare

            !$omp parallel do
            do iz = 0, nz
                gambot(iz) = gamtop(nz-iz)
            enddo
            !$omp end parallel do
            !Here gambot is the complement of gamtop.

        end subroutine init_inversion

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine finalise_inversion
            deallocate(green)
            deallocate(gamtop)
            deallocate(gambot)
            deallocate(phim)
            deallocate(phip)
#ifdef ENABLE_BUOYANCY
            deallocate(dphim)
            deallocate(dphip)
#endif
            deallocate(thetam)
            deallocate(thetap)
            deallocate(dthetam)
            deallocate(dthetap)
            deallocate(k2l2i)
            deallocate(k2l2)
            deallocate(filt)

            call finalise_fft

        end subroutine finalise_inversion

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! for kx > 0 and ky >= 0 or kx >= 0 and ky > 0
        subroutine set_hyperbolic_functions(kx, ky, zm, zp)
            integer,          intent(in) :: kx, ky
            double precision, intent(in) :: zm(0:nz), zp(0:nz)
            double precision             :: R(0:nz), Q(0:nz), k2ifac
#ifndef ENABLE_BUOYANCY
            double precision             :: dphim(0:nz), dphip(0:nz)
#endif
            double precision             :: ef, em(0:nz), ep(0:nz), Lm(0:nz), Lp(0:nz)
            double precision             :: fac, div, kl

            kl = dsqrt(k2l2(ky, kx))
            fac = kl * extent(3)
            ef = dexp(- fac)
#ifndef NDEBUG
            ! To avoid "Floating-point exception - erroneous arithmetic operation"
            ! when ef is really small.
            ef = max(ef, dsqrt(tiny(ef)))
#endif
            div = one / (one - ef**2)
            k2ifac = f12 * k2l2i(ky, kx)

            Lm = kl * zm
            Lp = kl * zp

            ep = dexp(- Lp)
            em = dexp(- Lm)

#ifndef NDEBUG
            ! To avoid "Floating-point exception - erroneous arithmetic operation"
            ! when ep and em are really small.
            ep = max(ep, dsqrt(tiny(ep)))
            em = max(em, dsqrt(tiny(em)))
#endif

            phim(:, ky, kx) = div * (ep - ef * em)
            phip(:, ky, kx) = div * (em - ef * ep)

#ifdef ENABLE_BUOYANCY
            dphim(:, ky, kx) = - kl * div * (ep + ef * em)
            dphip(:, ky, kx) =   kl * div * (em + ef * ep)
#else
            dphim = - kl * div * (ep + ef * em)
            dphip =   kl * div * (em + ef * ep)
#endif

            Q = div * (one + ef**2)
            R = div * two * ef

            thetam(:, ky, kx) = k2ifac * (R * Lm * phip(:, ky, kx) - Q * Lp * phim(:, ky, kx))
            thetap(:, ky, kx) = k2ifac * (R * Lp * phim(:, ky, kx) - Q * Lm * phip(:, ky, kx))

#ifdef ENABLE_BUOYANCY
            dthetam(:, ky, kx) = - k2ifac * ((Q * Lp - one) * dphim(:, ky, kx) - R * Lm * dphip(:, ky, kx))
            dthetap(:, ky, kx) = - k2ifac * ((Q * Lm - one) * dphip(:, ky, kx) - R * Lp * dphim(:, ky, kx))
#else
            dthetam(:, ky, kx) = - k2ifac * ((Q * Lp - one) * dphim - R * Lm * dphip)
            dthetap(:, ky, kx) = - k2ifac * ((Q * Lm - one) * dphip - R * Lp * dphim)
#endif
        end subroutine set_hyperbolic_functions

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! fc  - complete field (physical space)
        ! sf  - full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
        ! cfc - copy of complete field (physical space)
        subroutine field_decompose_physical(fc, sf)
            double precision, intent(in)  :: fc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision, intent(out) :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))

            call fftxyp2s(fc, sf)

            call field_decompose_semi_spectral(sf)

        end subroutine field_decompose_physical

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! in : complete field (semi-spectral space)
        ! out: full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
        subroutine field_decompose_semi_spectral(sfc)
            double precision, intent(inout) :: sfc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision                :: sfctop(box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            integer                         :: iz, kx, ky

            ! subtract harmonic part
            !$omp parallel do
            do iz = 1, nz-1
                sfc(iz, :, :) = sfc(iz, :, :) - (sfc(0, :, :) * phim(iz, :, :) + sfc(nz, :, :) * phip(iz, :, :))
            enddo
            !$omp end parallel do

            !$omp parallel workshare
            sfctop = sfc(nz, :, :)
            !$omp end parallel workshare

            ! transform interior to fully spectral
            !$omp parallel do collapse(2)
            do kx = box%lo(1), box%hi(1)
                do ky = box%lo(2), box%hi(2)
                    call dst(1, nz, sfc(1:nz, ky, kx), ztrig, zfactors)
                enddo
            enddo
            !$omp end parallel do

            !$omp parallel workshare
            sfc(nz, :, :) = sfctop
            !$omp end parallel workshare

        end subroutine field_decompose_semi_spectral

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! sf  - full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
        ! fc  - complete field (physical space)
        ! sfc - complete field (semi-spectral space)
        subroutine field_combine_physical(sf, fc)
            double precision, intent(in)  :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision, intent(out) :: fc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision              :: sfc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))

            sfc = sf

            call field_combine_semi_spectral(sfc)

            ! transform to physical space as fc:
            call fftxys2p(sfc, fc)

        end subroutine field_combine_physical

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! in : full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
        ! out: complete field (semi-spectral space)
        subroutine field_combine_semi_spectral(sf)
            double precision, intent(inout) :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision                :: sftop(box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            integer                         :: iz, kx, ky

            ! transform sf(1:nz-1, :, :) to semi-spectral space (sine transform) as the array sf:
            !$omp parallel workshare
            sftop = sf(nz, :, :)
            !$omp end parallel workshare

            !$omp parallel do collapse(2)
            do kx = box%lo(1), box%hi(1)
                do ky = box%lo(2), box%hi(2)
                    sf(nz, ky, kx) = zero
                    call dst(1, nz, sf(1:nz, ky, kx), ztrig, zfactors)
                enddo
            enddo
            !$omp end parallel do

            !$omp parallel workshare
            sf(nz, :, :) = sftop
            !$omp end parallel workshare

            ! add harmonic part to sfc:
            !$omp parallel do
            do iz = 1, nz-1
                sf(iz, :, :) = sf(iz, :, :) + sf(0, :, :) * phim(iz, :, :) + sf(nz, :, :) * phip(iz, :, :)
            enddo
            !$omp end parallel do

        end subroutine field_combine_semi_spectral

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Calculates df/dz for a field f using 2nd-order differencing.
        !Here fs = f, ds = df/dz. In physical or semi-spectral space.
        subroutine central_diffz(fs, ds)
            double precision, intent(in)  :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision, intent(out) :: ds(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            integer                       :: iz

            ! Linear extrapolation at the boundaries:
            ! iz = 0:  (fs(1) - fs(0)) / dz
            ! iz = nz: (fs(nz) - fs(nz-1)) / dz
            !$omp parallel workshare
            ds(0,  :, :) = dzi * (fs(1,    :, :) - fs(0,    :, :))
            ds(nz, :, :) = dzi * (fs(nz,   :, :) - fs(nz-1, :, :))
            !$omp end parallel workshare

            ! central differencing for interior cells
            !$omp parallel do private(iz) default(shared)
            do iz = 1, nz-1
                ds(iz, :, :) = (fs(iz+1, :, :) - fs(iz-1, :, :)) * hdzi
            enddo
            !$omp end parallel do

        end subroutine central_diffz

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#ifdef ENABLE_BUOYANCY
        !Calculates df/dz for a field f in mixed-spectral space
        !Here fs = f, ds = df/dz. Both fields are in mixed-spectral space.
        ! fs - mixed-spectral space
        ! ds - derivative linear part
        ! as - derivative sine part
        subroutine diffz(fs, ds)
            double precision, intent(in)  :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision, intent(out) :: ds(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision              :: as(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            integer                       :: kz, iz

            !Calculate the derivative of the linear part (ds) in semi-spectral space:
            !$omp parallel do private(iz)  default(shared)
            do iz = 0, nz
                ds(iz, :, :) = fs(0, :, :) * dphim(iz, :, :) + fs(nz, :, :) * dphip(iz, :, :)
            enddo
            !$omp end parallel do

            ! Calculate d/dz of this sine series:
            !$omp parallel workshare
            as(0, :, :) = zero
            !$omp end parallel workshare
            !$omp parallel do private(kz)  default(shared)
            do kz = 1, nz-1
                as(kz, :, :) = rkz(kz) * fs(kz, :, :)
            enddo
            !$omp end parallel do
            !$omp parallel workshare
            as(nz, :, :) = zero
            !$omp end parallel workshare

            !FFT these quantities back to semi-spectral space:
            call fftcosine(as)

            ! Combine vertical derivative given the sine (as) and linear (ds) parts:
            !omp parallel workshare
            ds = ds + as
            !omp end parallel workshare

            call field_decompose_semi_spectral(ds)

        end subroutine diffz
#endif

end module inversion_utils
