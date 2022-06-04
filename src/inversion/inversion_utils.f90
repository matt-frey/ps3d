module inversion_utils
    use constants
    use parameters, only : nx, ny, nz, dx, dxi, extent, ncelli
    use stafft
    use sta2dfft
    use deriv1d, only : init_deriv
    implicit none

    private

    ! Ordering in physical space: z, y, x
    ! Ordering in spectral space: z, x, y

    ! Tridiagonal arrays for the horizontal vorticity components:
    double precision, allocatable :: etdh(:, :, :), htdh(:, :, :)

    ! Tridiagonal arrays for the vertical vorticity component:
    double precision, allocatable :: etdv(:, :, :), htdv(:, :, :)

    ! Wavenumbers:
    double precision, allocatable :: rkx(:), hrkx(:), rky(:), hrky(:), rkz(:)

    ! Note k2l2i = 1/(k^2+l^2) (except k = l = 0, then k2l2i(0, 0) = 0)
    double precision, allocatable :: k2l2i(:, :)

    !Quantities needed in FFTs:
    double precision, allocatable :: xtrig(:), ytrig(:), ztrig(:)
    integer :: xfactors(5), yfactors(5), zfactors(5)
    integer, parameter :: nsubs_tri = 8 !number of blocks for openmp
    integer :: nxsub

    double precision, allocatable :: green(:, :, :)
    double precision, allocatable :: decz(:, :, :)

    ! Spectral dissipation operator
    double precision, allocatable :: hdis(:, :)

    ! Spectral filter:
    double precision, allocatable :: filt(:, :)

    ! Tridiagonal arrays for the vertical filter:
    double precision, allocatable :: etdf(:), htdf(:), am(:)

    ! Squared magnitude of horizontal wavevector, k^2 + l^2::
    double precision, parameter :: kksq = 8.0d0

    private :: xtrig, ytrig, xfactors, yfactors, & !zfactors, &
               hrkx, hrky!, rkz



    double precision :: dz, dzi, dz2, dz6, dz24, hdzi, dzisq, ap
    integer :: nwx, nwy, nxp2, nyp2

    logical :: is_fft_initialised = .false.

    public :: init_inversion        &
            , init_hyperdiffusion   &
            , diffx                 &
            , diffy                 &
            , diffz                 &
            , lapinv0               &
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
            , decz                  &
            , zfactors              &
            , ztrig                 &
            , rkx                   &
            , rky                   &
            , rkz                   &
            , apply_zfilter

    contains


        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine init_hyperdiffusion(bbdif, nnu, prediss, ke, en)
            double precision, intent(in) :: bbdif ! (bbdif = max(b) - min(b) at t = 0):
            integer,          intent(in) :: nnu
            double precision, intent(in) :: prediss
            double precision, intent(in) :: ke ! kinetic energy
            double precision, intent(in) :: en ! enstrophy
            double precision             :: visc, rkxmax, rkymax, rkzmax, K2max
            integer                      :: kx, ky, iz, kz

            allocate(hdis(0:nx-1, 0:ny-1))

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
                      hdis(kx, ky) = visc * (rkx(kx+1) ** 2 + rky(ky+1) ** 2)
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
                        hdis(kx, ky) = visc * (rkx(kx+1) ** 2 + rky(ky+1) ** 2) ** nnu
                    enddo
                enddo
            endif
        end subroutine init_hyperdiffusion

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine init_inversion
            double precision             :: fac, div
            integer                      :: kx, ky, iz, kz
            double precision             :: zh1(0:nz), zh0(0:nz)

            call init_fft

            allocate(green(0:nz, 0:nx-1, 0:ny-1))
            allocate(decz(0:nz, 0:nx-1, 0:ny-1))


            call init_tridiagonal

            !---------------------------------------------------------------------
            !Define Green function
            do ky = 0, ny-1
                do kx = 0, nx-1
                    do kz = 1, nz
                        green(kz, kx, ky) = - one / (rkx(kx+1) ** 2 + rky(ky+1) ** 2 + rkz(kz) ** 2)
                    enddo
                enddo
            enddo

            do ky = 0, ny-1
                do kx = 1, nx-1
                    green(0, kx, ky) = - one / (rkx(kx+1) ** 2 + rky(ky+1) ** 2)
                enddo
            enddo

            do ky = 1, ny-1
                green(0, 0, ky) = - one / rky(ky+1) ** 2
            enddo

            green(0, 0, 0) = zero

            !---------------------------------------------------------------------
            !Fractional z grid values:
            fac = one / dble(nz)
            do iz = 0, nz
                zh1(iz) = fac * dble(iz)
                zh0(iz) = one - zh1(iz)
            enddo

            !Hyperbolic functions used for solutions of Laplace's equation:
            do ky = 1, ny-1
                do kx = 0, nx-1
                    fac = dsqrt(rky(ky+1) ** 2 + rkx(kx+1) ** 2) * extent(3)
                    div = one / (one - dexp(-two * fac))
                    decz(:, kx, ky) = div * (exp(-fac * (one - zh1)) - &
                                             exp(-fac * (one + zh1)))
                enddo
            enddo

            do kx = 1, nx-1
                fac = rkx(kx+1) * extent(3)
                div = one / (one - dexp(-two * fac))
                decz(:, kx, 0) = div * (exp(-fac * (one - zh1)) - &
                                        exp(-fac * (one + zh1)))
            enddo

            decz(:, 0, 0) = zero

          end subroutine init_inversion

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Initialises this module (FFTs, x & y wavenumbers, tri-diagonal
        !coefficients, etc).
        subroutine init_fft
            double precision, allocatable :: a0(:, :), ksq(:, :)
            double precision              :: rkxmax, rkymax
            double precision              :: rksqmax
            double precision              :: kxmaxi, kymaxi
            integer                       :: kx, ky, iz, isub, ib_sub, ie_sub
            double precision              :: skx(nx), sky(ny)

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

            allocate(a0(nx, ny))
            allocate(ksq(nx, ny))
            allocate(k2l2i(nx, ny))

            allocate(filt(nx, ny))
            allocate(etdh(nz-1, nx, ny))
            allocate(htdh(nz-1, nx, ny))
            allocate(etdv(0:nz, nx, ny))
            allocate(htdv(0:nz, nx, ny))
            allocate(rkx(nx))
            allocate(hrkx(nx))
            allocate(rky(ny))
            allocate(hrky(ny))
            allocate(rkz(0:nz))
            allocate(xtrig(2 * nx))
            allocate(ytrig(2 * ny))
            allocate(ztrig(2*nz))

            nxsub = nx / nsubs_tri

            !----------------------------------------------------------------------
            ! Initialise FFTs and wavenumber arrays:
            call init2dfft(nx, ny, extent(1), extent(2), xfactors, yfactors, xtrig, ytrig, hrkx, hrky)
            call initfft(nz, zfactors, ztrig)

            !Define x wavenumbers:
            rkx(1) = zero
            do kx = 1, nwx-1
                rkx(kx+1)    = hrkx(2 * kx)
                rkx(nx+1-kx) = hrkx(2 * kx)
            enddo
            rkx(nwx+1) = hrkx(nx)
            rkxmax = hrkx(nx)

            !Define y wavenumbers:
            rky(1) = zero
            do ky = 1, nwy-1
                rky(ky+1)    = hrky(2 * ky)
                rky(ny+1-ky) = hrky(2 * ky)
            enddo
            rky(nwy+1) = hrky(ny)
            rkymax = hrky(ny)

            !Define z wavenumbers:
            rkz(0) = zero
            call init_deriv(nz, extent(3), rkz(1:nz))

            !Squared maximum total wavenumber:
            rksqmax = rkxmax ** 2 + rkymax ** 2

            !Squared wavenumber array (used in tridiagonal solve):
            do ky = 1, ny
                do kx = 1, nx
                    ksq(kx, ky) = rkx(kx) ** 2 + rky(ky) ** 2
                enddo
            enddo

            ksq(1, 1) = one
            k2l2i = one / ksq
            ksq(1, 1) = zero

            !----------------------------------------------------------
            !Define Hou and Li filter:
            kxmaxi = one / maxval(rkx)
            skx = -36.d0 * (kxmaxi * rkx) ** 36
            kymaxi = one/maxval(rky)
            sky = -36.d0 * (kymaxi * rky) ** 36
            do ky = 1, ny
                do kx = 1, nx
                    filt(kx, ky) = dexp(skx(kx) + sky(ky))
                enddo
            enddo

            !-----------------------------------------------------------------------
            ! Fixed coefficients used in the tridiagonal problems:
            a0 = -two * dzisq - ksq
            ap = dzisq

            !-----------------------------------------------------------------------
            ! Tridiagonal arrays for the horizontal vorticity components:
            htdh(1, :, :) = one / a0
            etdh(1, :, :) = -ap * htdh(1, :, :)
            !$omp parallel shared(a0, ap, etdh, htdh, nz, nxsub) private(isub, ib_sub, ie_sub, iz) default(none)
            !$omp do
            do isub = 0, nsubs_tri-1
                ib_sub = isub * nxsub + 1
                ie_sub = (isub + 1) * nxsub
                do iz = 2, nz-2
                    htdh(iz, ib_sub:ie_sub, :) = one / (a0(ib_sub:ie_sub, :) &
                                               + ap * etdh(iz-1, ib_sub:ie_sub, :))
                    etdh(iz, ib_sub:ie_sub, :) = -ap * htdh(iz, ib_sub:ie_sub, :)
                enddo
            enddo
            !$omp end do
            !$omp end parallel
            htdh(nz-1, :, :) = one / (a0 + ap * etdh(nz-2, :, :))
            ! Remove horizontally-averaged part (done separately):
            htdh(:, 1, 1) = zero
            etdh(:, 1, 1) = zero

            ! Tridiagonal arrays for the vertical vorticity component:
            htdv(0, :, :) = one / a0
            etdv(0, :, :) = -two * ap * htdv(0, :, :)
            !$omp parallel shared(a0, ap, etdv, htdv, nz, nxsub) private(isub, ib_sub, ie_sub, iz) default(none)
            !$omp do
            do isub = 0, nsubs_tri-1
                ib_sub = isub * nxsub + 1
                ie_sub = (isub + 1) * nxsub
                do iz = 1, nz-1
                    htdv(iz, ib_sub:ie_sub, :) = one / (a0(ib_sub:ie_sub, :) &
                                               + ap * etdv(iz-1, ib_sub:ie_sub, :))
                    etdv(iz, ib_sub:ie_sub, :) = -ap * htdv(iz, ib_sub:ie_sub, :)
                enddo
            enddo
            !$omp end do
            !$omp end parallel

            etdv(nz-1, 1, 1) = zero

            htdv(nz, :, :) = one / (a0 + two * ap * etdv(nz-1, :, :))
            ! Remove horizontally-averaged part (done separately):
            htdv(:, 1, 1) = zero
            etdv(:, 1, 1) = zero

            deallocate(a0)
            deallocate(ksq)
        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Initialises the tridiagonal problem for z-filtering
        subroutine init_tridiagonal
            double precision :: a0(nz-1)
            double precision :: pf, cf, diffmax, z, s
            integer          :: j

            allocate(etdf(nz-2))
            allocate(htdf(nz-1))
            allocate(am(nz))

            !-----------------------------------------------------------------------
            pf = dlog(two) / (one - (one - two / dble(nz)) ** 2)
            write(*,*) ' p = ', pf
            write(*,*)
            write(*,*) ' We take K_max = c * dz^2; enter c: '
            read(*,*) cf
            diffmax = cf * dx(3) ** 2
            write(*,*)
            write(*,*) ' K_min/K_max = ', dexp(-pf)

            ! Set up the tridiagonal system A x = u:

            !   | a0(1) ap(1)   0    ...   ...  ...   0    ||x(1)|   |u(1)|
            !   | am(2) a0(2) ap(2)   0    ...  ...   0    ||x(2)|   |u(2)|
            !   |   0   am(3) a0(3) ap(3)   0   ...   0    ||x(3)| = |u(3)|
            !   |                    ...                   || ...|   | ...|
            !   |   0    ...   ...   ...    0  am(n) a0(n) ||x(n)|   |u(n)|

            ! where n = nz-1 below, and here ap(i) = am(i+1) (symmetric system).

            do j=1,nz
                z = dx(3) * (dble(j) - f12)
                s = (two * z - extent(3)) / extent(3)
                am(j) = -cf * dexp(-pf * (one - s ** 2)) !note diffmax/dz^2 = cf
            enddo

            do j=1,nz-1
                z = dx(3) * dble(j)
                s = (two * z - extent(3)) / extent(3)
                a0(j) = one + kksq * diffmax * dexp(-pf * (one - s ** 2)) - am(j+1) - am(j)
            enddo

            htdf(1) = one / a0(1)
            etdf(1) = -am(2) * htdf(1)
            do j = 2, nz-2
                htdf(j) = one / (a0(j) + am(j) * etdf(j-1))
                etdf(j) = -am(j+1) * htdf(j)
            enddo
            htdf(nz-1) = one / (a0(nz-1) + am(nz-1) * etdf(nz-2))

        end subroutine init_tridiagonal

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Applies the z filter to fs while keeping
        ! the boundaries (iz = 0 and iz = nz).
        subroutine apply_zfilter(fs)
            double precision, intent(inout) :: fs(0:nz, 0:nx-1, 0:ny-1) ! semi-spectral space
            double precision                :: phi(0:nz)
            integer                         :: j, ix, iy

            do iy = 0, ny-1
                do ix = 0, nx-1

                    phi = fs(:, ix, iy)

                    fs(1, ix, iy) = (phi(1) - am(1) * phi(0)) * htdf(1)

                    do j = 2, nz-2
                        fs(j, ix, iy) = (phi(j) - am(j) * fs(j-1, ix, iy)) * htdf(j)
                    enddo

                    fs(nz-1, ix, iy) = (phi(nz-1) - am(nz) * phi(nz) - am(nz-1) * fs(nz-2, ix, iy)) * htdf(nz-1)

                    do j = nz-2, 1, -1
                        fs(j, ix, iy) = etdf(j) * fs(j+1, ix, iy) + fs(j, ix, iy)
                    enddo
                enddo
            enddo

        end subroutine apply_zfilter

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
            double precision, intent(in)  :: fs(0:nz, nx, ny)
            double precision, intent(out) :: ds(0:nz, nx, ny)
            integer                       :: iz

            ! linear extrapolation to fill boundary cells:
            ! iz = 0:  (fs(1) - fs(0)) / dz
            ! iz = nz: (fs(nz) - fs(nz-1)) / dz
            ! could try other method: one-sided differencing
            ds(0,  :, :) = dzi * (fs(1,    :, :) - fs(0,    :, :))
            ds(nz, :, :) = dzi * (fs(nz,   :, :) - fs(nz-1, :, :))

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

        !Inverts Laplace's operator on fs in semi-spectral space.
        !Here fs = 0 on the z boundaries.
        !Uses 2nd-order differencing
        !*** Overwrites fs ***
        subroutine lapinv0(fs)
            double precision, intent(inout) :: fs(0:nz, nx, ny)
            integer                         :: iz, isub, ib_sub, ie_sub

            fs(0, :, :) = zero
            fs(1, :, :) = fs(1, :, :) * htdh(1, :, :)

            !$omp parallel shared(fs, ap, htdh, nz, nxsub) private(isub, ib_sub, ie_sub, iz) default(none)
            !$omp do
            do isub = 0, nsubs_tri-1
                ib_sub = isub * nxsub + 1
                ie_sub = (isub + 1) * nxsub
                do iz = 2, nz-1
                    fs(iz, ib_sub:ie_sub, :) = (fs(iz, ib_sub:ie_sub, :) &
                                             - ap * fs(iz-1, ib_sub:ie_sub, :)) &
                                             * htdh(iz, ib_sub:ie_sub, :)
                enddo
            enddo
            !$omp end do
            !$omp end parallel

            fs(nz, :, :) = zero

            !$omp parallel shared(fs, etdh, nz, nxsub) private(isub, ib_sub, ie_sub, iz)  default(none)
            !$omp do
            do isub = 0, nsubs_tri-1
                ib_sub = isub * nxsub + 1
                ie_sub = (isub + 1) * nxsub
                do iz = nz-2, 1, -1
                    fs(iz, ib_sub:ie_sub, :) = etdh(iz, ib_sub:ie_sub, :) &
                                             * fs(iz+1, ib_sub:ie_sub, :) + fs(iz, ib_sub:ie_sub, :)
                enddo
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
