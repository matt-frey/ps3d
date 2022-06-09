module inversion_utils
    use constants
    use parameters, only : nx, ny, nz, dx, dxi, extent, ncelli
    use stafft
    use sta2dfft
    use deriv1d, only : init_deriv
    use options, only : prefilt
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
    double precision, allocatable :: phitop(:), phibot(:)
    double precision, allocatable :: psi(:, :, :), dpsi(:, :, :)

    private :: xtrig, ytrig, xfactors, yfactors, & !zfactors, &
               hrkx, hrky!, rkz



    double precision :: dz, dzi, dz2, dz6, dz24, hdzi, dzisq
    integer :: nwx, nwy, nxp2, nyp2

    logical :: is_fft_initialised = .false.

    public :: init_inversion        &
            , init_hyperdiffusion   &
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
            , phibot                &
            , phitop                &
            , gambot                &
            , gamtop                &
            , psi                   &
            , dpsi

    public :: field_combine_semi_spectral   &
            , field_combine_physical        &
            , field_decompose_semi_spectral &
            , field_decompose_physical

    contains


        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine init_hyperdiffusion(bbdif, nnu, prediss, ke, en)
            double precision, intent(in) :: bbdif ! (bbdif = max(b) - min(b) at t = 0):
            integer,          intent(in) :: nnu
            double precision, intent(in) :: prediss
            double precision, intent(in) :: ke ! kinetic energy
            double precision, intent(in) :: en ! enstrophy
            double precision             :: visc, rkxmax, rkymax, rkzmax, K2max
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
            if (nnu .eq. 1) then
                !Define viscosity:
                visc = prediss * sqrt(bbdif / rkxmax ** 3)
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
                K2max = rkxmax ** 2 + rkymax ** 2 + rkzmax ** 2
                ! multiply ke with ncelli to make it the mean kinetic energy
                visc = prediss *  (K2max * ke /en) ** f13 / (K2max ** nnu)
                !visc = prediss / max(rkxmax, rkymax, rkzmax) ** (2 * nnu)
                write(*,'(a,1p,e14.7)') ' Hyperviscosity nu = ', visc

                !Define dissipation operator:
                do ky = 0, ny-1
                    do kx = 0, nx-1
                        hdis(0,  kx, ky) = visc * k2l2(kx, ky) ** nnu
                        hdis(nz, kx, ky) = hdis(0,  kx, ky)
                        do kz = 1, nz-1
                            hdis(kz, kx, ky) = visc * (k2l2(kx, ky) + rkz(kz) ** 2) ** nnu
                        enddo
                    enddo
                enddo
            endif
        end subroutine init_hyperdiffusion

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine init_inversion
            integer                      :: kx, ky, iz, kz
            double precision             :: fac, faci, div, kl, kli
            double precision             :: em(1:nz-1), ep(1:nz-1)

            call init_fft

            allocate(green(1:nz-1, 0:nx-1, 0:ny-1))
            allocate(psi(1:nz-1, 0:nx-1, 0:ny-1))
            allocate(dpsi(0:nz, 0:nx-1, 0:ny-1))
            allocate(phitop(1:nz-1))
            allocate(phibot(1:nz-1))
            allocate(gamtop(0:nz))
            allocate(gambot(0:nz))

            !---------------------------------------------------------------------
            !Define Green function
            do kz = 1, nz-1
                green(kz, :, :) = - one / (k2l2 + rkz(kz) ** 2)
            enddo

            !---------------------------------------------------------------------
            !Define phitop = (z-lower(3)) / extent(3) --> phitop goes from 0 to 1
            fac = one / dble(nz)
            do iz = 1, nz-1
                phitop(iz) = fac * dble(iz)
                phibot(iz) = fac * dble(nz-iz)
            enddo
            !Here phibot is the complement of phitop.

            !Define gamtop as the integral of phitop with zero average:
            gamtop(0)  = -f16 * extent(3)
            gamtop(1:nz-1) = f12 * extent(3) * (phitop ** 2 - f13)
            gamtop(nz) =  f13 * extent(3)
            do iz = 0, nz
                gambot(iz)  = gamtop(nz-iz)
            enddo
            !Here gambot is the complement of gamtop.

            !Hyperbolic functions used for solutions of Laplace's equation:
            do ky = 1, ny-1
                do kx = 0, nx-1
                    kl  = dsqrt(k2l2(kx, ky))
                    kli = one / kl
                    fac = kl * extent(3)
                    faci = one / fac
                    div = one / (one - dexp(-two * fac))
                    em = dexp(-fac * (one - phitop))
                    ep = dexp(-fac * (one + phitop))
                    psi(:      , kx, ky) = k2l2i(kx, ky) * (div * (em - ep) - phitop)
                    dpsi(1:nz-1, kx, ky) =           kli * (div * (em + ep) - faci)

                    ! iz = 0  --> phitop = 0
                    dpsi(0 , kx, ky) = kli * (div * two * dexp(-fac) - faci)
                    ! iz = nz --> phitop = 1
                    dpsi(nz, kx, ky) = kli * (div * (one + dexp(-two * fac)) - faci)
                enddo
            enddo

            ! ky = 0
            do kx = 1, nx-1
                kli = one / rkx(kx)
                fac = rkx(kx) * extent(3)
                faci = one / fac
                div = one / (one - dexp(-two * fac))
                em = dexp(-fac * (one - phitop))
                ep = dexp(-fac * (one + phitop))
                psi(:      , kx, 0) = k2l2i(kx, 0) * (div * (em - ep) - phitop)
                dpsi(1:nz-1, kx, 0) =          kli * (div * (em + ep) - faci)

                ! iz = 0  --> phitop = 0
                dpsi(0 ,  kx, 0) = kli * (div * two * dexp(-fac) - faci)
                ! iz = nz --> phitop = 1
                dpsi(nz, kx, 0) = kli * (div * (one + dexp(-two * fac)) - faci)
            enddo

            ! kx = ky = 0
            psi(: , 0, 0) = zero
            dpsi(:, 0, 0) = zero

          end subroutine init_inversion

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Initialises this module (FFTs, x & y wavenumbers, tri-diagonal
        !coefficients, etc).
        subroutine init_fft
            double precision              :: rkxmax, rkymax, rkzmax
            double precision              :: rksqmax
!             double precision              :: kxmaxi, kymaxi, kzmaxi
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

!             !----------------------------------------------------------
!             !Define Hou and Li filter (2D and 3D):
!             kxmaxi = one / maxval(rkx)
!             skx = -36.d0 * (kxmaxi * rkx) ** 36
!             kymaxi = one/maxval(rky)
!             sky = -36.d0 * (kymaxi * rky) ** 36
!             kzmaxi = one/maxval(rkz)
!             skz = -36.d0 * (kzmaxi * rkz) ** 36
!             do ky = 0, ny-1
!                 do kx = 0, nx-1
!                     filt(0,  kx, ky) = dexp(skx(kx) + sky(ky))
!                     filt(nz, kx, ky) = filt(0, kx, ky)
!                     do kz = 1, nz-1
!                         filt(kz, kx, ky) = filt(0, kx, ky) * dexp(skz(kz))
!                     enddo
!                 enddo
!             enddo


            !---------------------------------------------------------------------
            !Define de-aliasing filter (2/3 rule):
            skx(0) = one
            do kx = 1, nx-1
                if (rkx(kx) .lt. f23 * rkxmax) then
                    skx(kx) = one
                else
                    skx(kx) = zero
                endif
            enddo

            sky(0) = one
            do ky = 1, ny-1
                if (rky(ky) .lt. f23 * rkymax) then
                    sky(ky) = one
                else
                    sky(ky) = zero
                endif
            enddo

            skz(0) = one
            rkzmax = maxval(rkz)
            do kz = 1, nz
                if (rkz(kz) .lt. f23 * rkzmax) then
                    skz(kz) = one
                else
                    skz(kz) = zero
                endif
            enddo

            !Take product of 1d filters:
            do ky = 0, ny-1
                do kx = 0, nx-1
                    filt(0,  kx, ky) = skx(kx) * sky(ky)
                    filt(nz, kx, ky) = filt(0, kx, ky)
                    do kz = 1, nz-1
                        filt(kx, kx, ky) = filt(0, kx, ky) * skz(kz)
                    enddo
                enddo
            enddo

        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_decompose_semi_spectral(sfc)
            double precision, intent(inout) :: sfc(0:nz, 0:nx-1, 0:ny-1) ! in : complete field (semi-spectral space)
                                                                         ! out: full-spectral (1:nz-1),
                                                                         !      semi-spectral at iz = 0 and iz = nz
            double precision                :: sfctop(0:nx-1, 0:ny-1)
            integer                         :: iz, kx, ky

            ! get linear part
            do iz = 1, nz-1
                sfc(iz, :, :) = sfc(iz, :, :) - (sfc(0, :, :) * phibot(iz) + sfc(nz, :, :) * phitop(iz))
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

            ! get linear part and add to sfc:
            do iz = 1, nz-1
                sf(iz, :, :) = sf(iz, :, :) + sf(0, :, :) * phibot(iz) + sf(nz, :, :) * phitop(iz)
            enddo

        end subroutine field_combine_semi_spectral

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
