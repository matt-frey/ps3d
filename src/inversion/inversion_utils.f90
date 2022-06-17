module inversion_utils
    use constants
    use parameters, only : nx, ny, nz, dx, dxi, extent, ncelli, upper, lower
    use stafft
    use sta2dfft
    use deriv1d, only : init_deriv
    use options, only : viscosity
    implicit none

    private

    ! Ordering in physical space: z, y, x
    ! Ordering in spectral space: z, x, y

    ! Wavenumbers:
    double precision, allocatable :: rkx(:), hrkx(:), rky(:), hrky(:), rkz(:), rkzi(:)

    ! Note k2l2i = 1/(k^2+l^2) (except k = l = 0, then k2l2i(0, 0) = 0)
    double precision, allocatable :: k2l2i(:, :)

    ! Note k2l2 = k^2+l^2
    double precision, allocatable :: k2l2(:, :)

    !Quantities needed in FFTs:
    double precision, allocatable :: xtrig(:), ytrig(:), ztrig(:)
    integer :: xfactors(5), yfactors(5), zfactors(5)

    double precision, allocatable :: green(:, :, :)

    ! Spectral dissipation operator
    double precision, allocatable :: hdis(:, :, :)

    ! Spectral filter:
    double precision, allocatable :: filt(:, :, :)

    double precision, allocatable :: gamtop(:), gambot(:)

    double precision, allocatable :: thetam(:, :, :)    ! theta_{-}
    double precision, allocatable :: thetap(:, :, :)    ! theta_{+}
    double precision, allocatable :: dthetam(:, :, :)   ! dtheta_{-}/dz
    double precision, allocatable :: dthetap(:, :, :)   ! dtheta_{+}/dz
    double precision, allocatable :: phim(:, :, :)      ! phi_{-}
    double precision, allocatable :: phip(:, :, :)      ! phi_{+}


    private :: xtrig, ytrig, xfactors, yfactors, & !zfactors, &
               hrkx, hrky!, rkz


    double precision :: dz, dzi, dz2, dz6, dz24, hdzi, dzisq
    integer :: nwx, nwy, nxp2, nyp2

    logical :: is_fft_initialised = .false.

    public :: init_inversion        &
            , init_diffusion        &
            , diffx                 &
            , diffy                 &
            , diffz                 &
            , fftxyp2s              &
            , fftxys2p              &
            , dz2                   &
            , filt                  &
            , hdzi                  &
            , k2l2i                 &
            , hdis                  &
            , fftczp2s              &
            , fftczs2p              &
            , fftss2fs              &
            , fftfs2ss              &
            , green                 &
            , zfactors              &
            , ztrig                 &
            , rkx                   &
            , rky                   &
            , rkz                   &
            , rkzi                  &
            , thetap                &
            , thetam                &
            , dthetap               &
            , dthetam               &
            , gambot                &
            , gamtop

    public :: field_combine_semi_spectral   &
            , field_combine_physical        &
            , field_decompose_semi_spectral &
            , field_decompose_physical

    contains


        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine init_diffusion(bbdif, ke, en)
            double precision, intent(in) :: bbdif ! (bbdif = max(b) - min(b) at t = 0):
            double precision, intent(in) :: ke ! kinetic energy
            double precision, intent(in) :: en ! enstrophy
            double precision             :: rkxmax, rkymax, rkzmax, K2max
            double precision             :: visc, wfac
            integer                      :: kx, ky, kz

            allocate(hdis(0:nz, 0:nx-1, 0:ny-1))

            ! check if FFT is initialised
            if (.not. is_fft_initialised) then
                print *, "Error: FFT not initialised!"
                stop
            endif

            rkxmax = maxval(rkx)
            rkymax = maxval(rky)
            rkzmax = maxval(rkz)

            !---------------------------------------------------------------------
            ! Damping, viscous or hyperviscous:
            if (viscosity%nnu .eq. 1) then
                !Define viscosity:
                visc = viscosity%prediss * sqrt(bbdif / rkxmax ** 3)
                write(*,'(a,1p,e14.7)') ' Viscosity nu = ', visc

                !Define spectral dissipation operator:
                do ky = 0, ny-1
                    do kx = 0, nx-1
                        hdis(0,  kx, ky) = visc * k2l2(kx, ky)
                        hdis(nz, kx, ky) = hdis(0,  kx, ky)
                        do kz = 1, nz-1
                            hdis(kz, kx, ky) = visc * (k2l2(kx, ky) + rkz(kz) ** 2)
                        enddo
                    enddo
                enddo
            else
                !Define hyperviscosity:
                K2max = max(rkxmax, rkymax, rkzmax) ** 2
                wfac = one / K2max
                visc = viscosity%prediss *  (K2max * ke /en) ** f13
                write(*,'(a,1p,e14.7)') ' Hyperviscosity nu = ', visc * wfac ** viscosity%nnu

                !Define dissipation operator:
                do ky = 0, ny-1
                    do kx = 0, nx-1
                        hdis(0,  kx, ky) = visc * (wfac * k2l2(kx, ky)) ** viscosity%nnu
                        hdis(nz, kx, ky) = hdis(0,  kx, ky)
                        do kz = 1, nz-1
                            hdis(kz, kx, ky) = visc * (wfac * (k2l2(kx, ky) + rkz(kz) ** 2)) ** viscosity%nnu
                        enddo
                    enddo
                enddo
                
                !Ensure average is not modified by hyperviscosity:
                hdis(:, 0, 0) = zero
            endif
        end subroutine init_diffusion

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine init_inversion
            integer          :: kx, ky, iz, kz
            double precision :: z, zm(0:nz), zp(0:nz), fac

            call init_fft

            allocate(green(1:nz-1, 0:nx-1, 0:ny-1))
            allocate(gamtop(0:nz))
            allocate(gambot(0:nz))

            allocate(phim(0:nz, 0:nx-1, 0:ny-1))
            allocate(phip(0:nz, 0:nx-1, 0:ny-1))
            allocate(thetam(0:nz, 0:nx-1, 0:ny-1))
            allocate(thetap(0:nz, 0:nx-1, 0:ny-1))
            allocate(dthetam(0:nz, 0:nx-1, 0:ny-1))
            allocate(dthetap(0:nz, 0:nx-1, 0:ny-1))

            !---------------------------------------------------------------------
            !Define Green function
            do kz = 1, nz-1
                green(kz, :, :) = - one / (k2l2 + rkz(kz) ** 2)
            enddo

            !---------------------------------------------------------------------
            !Define zm = zmax - z, zp = z - zmin
            do iz = 0, nz
                z = lower(3) + dx(3) * dble(iz)
                zm(iz) = upper(3) - z
                zp(iz) = z - lower(3)
            enddo

            !Hyperbolic functions used for solutions of Laplace's equation:
            do ky = 1, ny-1
                do kx = 0, nx-1
                    call set_hyperbolic_functions(kx, ky, zm, zp)
                enddo
            enddo

            ! ky = 0
            do kx = 1, nx-1
                call set_hyperbolic_functions(kx, 0, zm, zp)
            enddo

            ! kx = ky = 0
            phim(:, 0, 0) = zm / extent(3)
            phip(:, 0, 0) = zp / extent(3)

            thetam(:, 0, 0) = zero
            thetap(:, 0, 0) = zero

            dthetam(:, 0, 0) = zero
            dthetap(:, 0, 0) = zero

            !---------------------------------------------------------------------
            !Define gamtop as the integral of phip(iz, 0, 0) with zero average:
            gamtop = f12 * extent(3) * (phip(:, 0, 0) ** 2 - f13)
            do iz = 0, nz
                gambot(iz) = gamtop(nz-iz)
            enddo
            !Here gambot is the complement of gamtop.

        end subroutine init_inversion

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! for kx > 0 and ky >= 0 or kx >= 0 and ky > 0
        subroutine set_hyperbolic_functions(kx, ky, zm, zp)
            integer,          intent(in) :: kx, ky
            double precision, intent(in) :: zm(0:nz), zp(0:nz)
            double precision             :: R(0:nz), Q(0:nz), k2ifac, dphim(0:nz), dphip(0:nz)
            double precision             :: ef(0:nz), em(0:nz), ep(0:nz), Lm(0:nz), Lp(0:nz)
            double precision             :: fac, div, kl

            kl = dsqrt(k2l2(kx, ky))
            fac = kl * extent(3)
            ef = dexp(- fac)
            div = one / (one - ef**2)
            k2ifac = f12 * k2l2i(kx, ky)

            Lm = kl * zm
            Lp = kl * zp

            ep = dexp(- Lp)
            em = dexp(- Lm)

            phim(:, kx, ky) = div * (ep - ef * em)
            phip(:, kx, ky) = div * (em - ef * ep)

            dphim = - kl * div * (ep + ef * em)
            dphip =   kl * div * (em + ef * ep)

            Q = div * (one + ef**2)
            R = div * two * ef

            thetam(:, kx, ky) = k2ifac * (R * Lm * phip(:, kx, ky) - Q * Lp * phim(:, kx, ky))
            thetap(:, kx, ky) = k2ifac * (R * Lp * phim(:, kx, ky) - Q * Lm * phip(:, kx, ky))

            dthetam(:, kx, ky) = - k2ifac * ((Q * Lp - one) * dphim - R * Lm * dphip)
            dthetap(:, kx, ky) = - k2ifac * ((Q * Lm - one) * dphip - R * Lp * dphim)

        end subroutine set_hyperbolic_functions

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Initialises this module (FFTs, x & y wavenumbers, tri-diagonal
        !coefficients, etc).
        subroutine init_fft
            double precision              :: rkxmax, rkymax!, rkzmax
            double precision              :: rksqmax
            double precision              :: kxmaxi, kymaxi, kzmaxi
            integer                       :: kx, ky, kz
            double precision              :: skx(0:nx-1), sky(0:ny-1), skz(0:nz)

            if (is_fft_initialised) then
                return
            endif

            is_fft_initialised = .true.

            dz = dx(3)
            dzi = dxi(3)
            dz6  = f16 * dx(3)
            dz2  = f12 * dx(3)
            dz24 = f124 * dx(3)
            dzisq = dxi(3) ** 2
            hdzi = f12 * dxi(3)
            nwx = nx / 2
            nwy = ny / 2
            nyp2 = ny + 2
            nxp2 = nx + 2

            allocate(k2l2i(0:nx-1, 0:ny-1))
            allocate(k2l2(0:nx-1, 0:ny-1))

            allocate(filt(0:nz, 0:nx-1, 0:ny-1))
            allocate(rkx(0:nx-1))
            allocate(hrkx(nx))
            allocate(rky(0:ny-1))
            allocate(hrky(ny))
            allocate(rkz(0:nz))
            allocate(rkzi(1:nz-1))
            allocate(xtrig(2 * nx))
            allocate(ytrig(2 * ny))
            allocate(ztrig(2*nz))

            !----------------------------------------------------------------------
            ! Initialise FFTs and wavenumber arrays:
            call init2dfft(nx, ny, extent(1), extent(2), xfactors, yfactors, xtrig, ytrig, hrkx, hrky)
            call initfft(nz, zfactors, ztrig)

            !Define x wavenumbers:
            rkx(0) = zero
            do kx = 1, nwx-1
                rkx(kx)    = hrkx(2 * kx)
                rkx(nx-kx) = hrkx(2 * kx)
            enddo
            rkx(nwx) = hrkx(nx)
            rkxmax = hrkx(nx)

            !Define y wavenumbers:
            rky(0) = zero
            do ky = 1, nwy-1
                rky(ky)    = hrky(2 * ky)
                rky(ny-ky) = hrky(2 * ky)
            enddo
            rky(nwy) = hrky(ny)
            rkymax = hrky(ny)

            !Define z wavenumbers:
            rkz(0) = zero
            call init_deriv(nz, extent(3), rkz(1:nz))
            rkzi(1:nz-1) = one / rkz(1:nz-1)

            !Squared maximum total wavenumber:
            rksqmax = rkxmax ** 2 + rkymax ** 2

            !Squared wavenumber array (used in tridiagonal solve):
            do ky = 0, ny-1
                do kx = 0, nx-1
                    k2l2(kx, ky) = rkx(kx) ** 2 + rky(ky) ** 2
                enddo
            enddo

            k2l2(0, 0) = one
            k2l2i = one / k2l2
            k2l2(0, 0) = zero

             !----------------------------------------------------------
             !Define Hou and Li filter (2D and 3D):
             kxmaxi = one / maxval(rkx)
             skx = -36.d0 * (kxmaxi * rkx) ** 36
             kymaxi = one/maxval(rky)
             sky = -36.d0 * (kymaxi * rky) ** 36
             kzmaxi = one/maxval(rkz)
             skz = -36.d0 * (kzmaxi * rkz) ** 36
             do ky = 0, ny-1
                 do kx = 0, nx-1
                     filt(0,  kx, ky) = dexp(skx(kx) + sky(ky))
                     filt(nz, kx, ky) = filt(0, kx, ky)
                     do kz = 1, nz-1
                         filt(kz, kx, ky) = filt(0, kx, ky) * dexp(skz(kz))
                     enddo
                 enddo
             enddo
             
             !Ensure filter does not change domain mean:
             filt(:, 0, 0) = one


!            !---------------------------------------------------------------------
!            !Define de-aliasing filter (2/3 rule):
!            skx(0) = one
!            do kx = 1, nx-1
!                if (rkx(kx) .lt. f23 * rkxmax) then
!                    skx(kx) = one
!                else
!                    skx(kx) = zero
!                endif
!            enddo
!
!            sky(0) = one
!            do ky = 1, ny-1
!                if (rky(ky) .lt. f23 * rkymax) then
!                    sky(ky) = one
!                else
!                    sky(ky) = zero
!                endif
!            enddo
!
!            skz(0) = one
!            rkzmax = maxval(rkz)
!            do kz = 1, nz
!                if (rkz(kz) .lt. f23 * rkzmax) then
!                    skz(kz) = one
!                else
!                    skz(kz) = zero
!                endif
!            enddo
!
!            !Take product of 1d filters:
!            do ky = 0, ny-1
!                do kx = 0, nx-1
!                    filt(0,  kx, ky) = skx(kx) * sky(ky)
!                    filt(nz, kx, ky) = filt(0, kx, ky)
!                    do kz = 1, nz-1
!                        filt(kz, kx, ky) = filt(0, kx, ky) * skz(kz)
!                    enddo
!                enddo
!            enddo
             
!            !Ensure filter does not change domain mean:
!            filt(:, 0, 0) = one

        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_decompose_physical(fc, sf)
            double precision, intent(in)  :: fc(0:nz, 0:ny-1, 0:nx-1)    ! complete field (physical space)
            double precision, intent(out) :: sf(0:nz, 0:nx-1, 0:ny-1)    ! full-spectral (1:nz-1),
                                                                         ! semi-spectral at iz = 0 and iz = nz
            double precision              :: cfc(0:nz, 0:ny-1, 0:nx-1)   ! copy of complete field (physical space)

            cfc = fc
            call fftxyp2s(cfc, sf)

            call field_decompose_semi_spectral(sf)

        end subroutine field_decompose_physical

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_decompose_semi_spectral(sfc)
            double precision, intent(inout) :: sfc(0:nz, 0:nx-1, 0:ny-1) ! in : complete field (semi-spectral space)
                                                                         ! out: full-spectral (1:nz-1),
                                                                         !      semi-spectral at iz = 0 and iz = nz
            double precision                :: sfctop(0:nx-1, 0:ny-1)
            integer                         :: iz, kx, ky

            ! subtract harmonic part
            do iz = 1, nz-1
                sfc(iz, :, :) = sfc(iz, :, :) - (sfc(0, :, :) * phim(iz, :, :) + sfc(nz, :, :) * phip(iz, :, :))
            enddo

            sfctop = sfc(nz, :, :)

            ! transform interior to fully spectral
            do ky = 0, ny-1
                do kx = 0, nx-1
                    call dst(1, nz, sfc(1:nz, kx, ky), ztrig, zfactors)
                enddo
            enddo

            sfc(nz, :, :) = sfctop

        end subroutine field_decompose_semi_spectral

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_combine_physical(sf, fc)
            double precision, intent(in)  :: sf(0:nz, 0:nx-1, 0:ny-1)    ! full-spectral (1:nz-1),
                                                                         ! semi-spectral at iz = 0 and iz = nz
            double precision, intent(out) :: fc(0:nz, 0:ny-1, 0:nx-1)    ! complete field (physical space)
            double precision              :: sfc(0:nz, 0:nx-1, 0:ny-1)   ! complete field (semi-spectral space)

            sfc = sf

            call field_combine_semi_spectral(sfc)

            ! transform to physical space as fc:
            call fftxys2p(sfc, fc)

        end subroutine field_combine_physical

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_combine_semi_spectral(sf)
            double precision, intent(inout) :: sf(0:nz, 0:nx-1, 0:ny-1) ! in: full-spectral (1:nz-1),
                                                                        !     semi-spectral at iz = 0 and iz = nz
                                                                        ! out: complete field (semi-spectral space)
            double precision                :: sftop(0:nx-1, 0:ny-1)
            integer                         :: iz, kx, ky

            ! transform sf(1:nz-1, :, :) to semi-spectral space (sine transform) as the array sf:
            sftop = sf(nz, :, :)
            do ky = 0, ny-1
                do kx = 0, nx-1
                    sf(nz, kx, ky) = zero
                    call dst(1, nz, sf(1:nz, kx, ky), ztrig, zfactors)
                enddo
            enddo
            sf(nz, :, :) = sftop

            ! add harmonic part to sfc:
            do iz = 1, nz-1
                sf(iz, :, :) = sf(iz, :, :) + sf(0, :, :) * phim(iz, :, :) + sf(nz, :, :) * phip(iz, :, :)
            enddo

        end subroutine field_combine_semi_spectral

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Given fs in spectral space (at least in x & y), this returns dfs/dx
        ! (partial derivative).  The result is returned in ds, again
        ! spectral.  Uses exact form of the derivative in spectral space.
        subroutine diffx(fs,ds)
            double precision, intent(in)  :: fs(0:nz, nx, ny)
            double precision, intent(out) :: ds(0:nz, nx, ny)
            integer                       :: kx, dkx, kxc

            !Carry out differentiation by wavenumber multiplication:
            ds(:, 1, :) = zero
            do kx = 2, nx - nwx
                dkx = 2 * (kx - 1)
                kxc = nxp2 - kx
                ds(:, kx,  :) = -hrkx(dkx) * fs(:,kxc,:)
                ds(:, kxc, :) =  hrkx(dkx) * fs(:,kx ,:)
            enddo

            if (mod(nx, 2) .eq. 0) then
                kxc = nwx + 1
                ds(:, kxc, :) = zero
            endif
        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Given fs in spectral space (at least in x & y), this returns dfs/dy
        ! (partial derivative).  The result is returned in ds, again
        ! spectral.  Uses exact form of the derivative in spectral space.
        subroutine diffy(fs,ds)
            double precision, intent(in)  :: fs(0:nz, nx, ny)
            double precision, intent(out) :: ds(0:nz, nx, ny)
            double precision              :: fac
            integer                       :: ky, kyc

            !Carry out differentiation by wavenumber multiplication:
            ds(:, :, 1) = zero

            do ky = 2, ny - nwy
                kyc = nyp2 - ky
                fac = hrky(2 * (ky - 1))
                ds(:, :, ky) = -fac * fs(:, :, kyc)
                ds(:, :, kyc) = fac * fs(:, : , ky)
            enddo

            if (mod(ny, 2) .eq. 0) then
                kyc = nwy + 1
                ds(:, :, kyc) = zero
            endif
        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Calculates df/dz for a field f using 2nd-order differencing.
        !Here fs = f, ds = df/dz.
        subroutine diffz(fs, ds)
            double precision, intent(in)  :: fs(0:nz, 0:ny-1, 0:nx-1)
            double precision, intent(out) :: ds(0:nz, 0:ny-1, 0:nx-1)
            integer                       :: iz

            ! Quadratic extrapolation at boundaries:
            ds(0,  :, :) = hdzi * (four * fs(1,  :, :) - three * fs(0,    :, :) - fs(2, :, :))
            ds(nz, :, :) = hdzi * (three * fs(nz, :, :) + fs(nz-2, :, :) - four * fs(nz-1, :, :))

            ! central differencing for interior cells
            !$omp parallel shared(ds, fs, hdzi, nz) private(iz) default(none)
            !$omp do
            do iz = 1, nz-1
                ds(iz, :, :) = (fs(iz+1, :, :) - fs(iz-1, :, :)) * hdzi
            enddo
            !$omp end do
            !$omp end parallel

        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Computes a 2D FFT (in x & y) of a 3D array fp in physical space
        ! and returns the result as fs in spectral space (in x & y).
        ! Only FFTs over the x and y directions are performed.
        ! *** fp is destroyed upon exit ***
        subroutine fftxyp2s(fp, fs)
            double precision, intent(inout) :: fp(:, :, :)       !Physical
            double precision, intent(out)   :: fs(:, :, :)       !Spectral
            integer                         :: kx, iy, nzval, nxval, nyval

            nzval = size(fp, 1)
            nyval = size(fp, 2)
            nxval = size(fp, 3)

            ! Carry out a full x transform first:
            call forfft(nzval * nyval, nxval, fp, xtrig, xfactors)

            ! Transpose array:
            !$omp parallel do shared(fs, fp) private(kx, iy)
            do kx = 1, nxval
                do iy = 1, nyval
                    fs(:, kx, iy) = fp(:, iy, kx)
                enddo
            enddo
            !$omp end parallel do

            ! Carry out a full y transform on transposed array:
            call forfft(nzval * nxval, nyval, fs, ytrig, yfactors)
        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Computes an *inverse* 2D FFT (in x & y) of a 3D array fs in spectral
        ! space and returns the result as fp in physical space (in x & y).
        ! Only inverse FFTs over the x and y directions are performed.
        ! *** fs is destroyed upon exit ***
        subroutine fftxys2p(fs, fp)
            double precision, intent(inout):: fs(:, :, :)  !Spectral
            double precision, intent(out)  :: fp(:, :, :)  !Physical
            integer                        :: kx, iy, nzval, nxval, nyval

            nzval = size(fs, 1)
            nxval = size(fs, 2)
            nyval = size(fs, 3)

            ! Carry out a full inverse y transform first:
            call revfft(nzval * nxval, nyval, fs, ytrig, yfactors)

            ! Transpose array:
            !$omp parallel do shared(fs, fp) private(kx, iy)
            do kx = 1, nxval
                do iy = 1, nyval
                    fp(:, iy, kx) = fs(:, kx, iy)
                enddo
            enddo
            !$omp end parallel do

            ! Carry out a full inverse x transform:
            call revfft(nzval * nyval, nxval, fp, xtrig, xfactors)
        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Computes a 3D FFT of an array fp in physical space and
        !returns the result as fs in spectral space.  It is assumed that
        !fp is generally non-zero at the z boundaries (so a cosine
        !transform is used in z).
        !*** fp is destroyed upon exit ***
        subroutine fftczp2s(fp, fs)
            double precision, intent(inout) :: fp(:, :, :)  !Physical
            double precision, intent(out)   :: fs(:, :, :)  !Spectral
            integer                         :: kx, ky, iy, nzval, nxval, nyval

            nzval = size(fp, 1)
            nyval = size(fp, 2)
            nxval = size(fp, 3)

            !Carry out a full x transform first:
            call forfft(nzval * nyval, nxval, fp, xtrig, xfactors)

            !Transpose array:
            do kx = 1, nxval
                do iy = 1, nyval
                    fs(:, kx, iy) = fp(:, iy, kx)
                enddo
            enddo

            !Carry out a full y transform on transposed array:
            call forfft(nzval * nxval, nyval, fs, ytrig, yfactors)

            !Carry out z FFT for each kx and ky:
            do ky = 1, nyval
                do kx = 1, nxval
                    call dct(1, nzval-1, fs(:, kx, ky), ztrig, zfactors)
                enddo
            enddo
        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Computes an *inverse* 3D FFT of an array fs in spectral space and
        !returns the result as fp in physical space.  It is assumed that
        !fp is generally non-zero at the z boundaries (so a cosine
        !transform is used in z).
        !*** fs is destroyed upon exit ***
        subroutine fftczs2p(fs, fp)
            double precision, intent(inout) :: fs(:, :, :)  !Spectral
            double precision, intent(out)   :: fp(:, :, :)  !Physical
            integer                         :: ix, iy, kx, nzval, nxval, nyval

            nzval = size(fs, 1)
            nxval = size(fs, 2)
            nyval = size(fs, 3)

            !Carry out a full inverse y transform first:
            call revfft(nzval * nxval, nyval,fs,ytrig,yfactors)

            !Transpose array:
            do kx = 1, nxval
                do iy = 1, nyval
                    fp(:, iy, kx) = fs(:, kx, iy)
                enddo
            enddo

            !Carry out a full inverse x transform:
            call revfft(nzval * nyval, nxval, fp, xtrig, xfactors)

            !Carry out z FFT for each ix and iy:
            do ix = 1, nxval
                do iy = 1, nyval
                    call dct(1, nzval-1, fp(:, iy, ix), ztrig, zfactors)
                enddo
            enddo
        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Computes a 3D FFT of an array fp in semi-spectral space and
        !returns the result as fs in fully spectral space.  It is assumed that
        !fp is generally non-zero at the z boundaries (so a cosine
        !transform is used in z).
        subroutine fftss2fs(fp, fs)
            double precision, intent(in)  :: fp(:, :, :)  !semi-spectral
            double precision, intent(out) :: fs(:, :, :)  !fully-spectral
            integer                       :: kx, ky, nzval, nxval, nyval

            nzval = size(fp, 1)
            nxval = size(fp, 2)
            nyval = size(fp, 3)

            !Carry out z FFT for each kx and ky:
            do ky = 1, nyval
                do kx = 1, nxval
                    fs(:, kx, ky) = fp(:, kx, ky)
                    call dct(1, nzval-1, fs(:, kx, ky), ztrig, zfactors)
                enddo
            enddo
        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Computes an *inverse* 3D FFT of an array fs in fully spectral space and
        !returns the result as fp in semi-spectral space.  It is assumed that
        !fp is generally non-zero at the z boundaries (so a cosine
        !transform is used in z).
        subroutine fftfs2ss(fs, fp)
            double precision, intent(in)  :: fs(:, :, :)  !fully spectral
            double precision, intent(out) :: fp(:, :, :)  !semi-spectral
            integer                       :: ix, iy, nzval, nxval, nyval

            nzval = size(fs, 1)
            nxval = size(fs, 2)
            nyval = size(fs, 3)

            !Carry out z FFT for each ix and iy:
            do iy = 1, nyval
                do ix = 1, nxval
                    fp(:, ix, iy) = fs(:, ix, iy)
                    call dct(1, nzval-1, fp(:, ix, iy), ztrig, zfactors)
                enddo
            enddo
        end subroutine

end module inversion_utils
