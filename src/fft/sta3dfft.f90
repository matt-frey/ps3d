module sta3dfft
    use constants, only : zero, f12
    use stafft
    use sta2dfft
    use deriv1d
    implicit none

    ! Spectral filter:
    double precision, allocatable :: filt(:, :, :)

    ! Wavenumbers::
    double precision, allocatable :: rkx(:), hrkx(:), rky(:), hrky(:), rkz(:)

    ! Quantities needed in FFTs:
    double precision, allocatable :: xtrig(:), ytrig(:), ztrig(:)
    integer                       :: xfactors(5), yfactors(5), zfactors(5)


    double precision, allocatable :: skx(:), sky(:), skz(:)

    private :: xtrig, ytrig, xfactors, yfactors,   &
               rkx, hrkx, rky, hrky

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        subroutine init3dfft(nx, ny, nz, extent)
            integer,          intent(in) :: nx, ny, nz
            double precision, intent(in) :: extent(3)
            double precision             :: kxmaxi, kymaxi, kzmaxi
!             double precision             :: kmax, kc !cfilt,
            integer                      :: nwx, nwy, kx, ky, kz

            if (.not. allocated(filt)) then
                allocate(filt(0:nz, nx, ny))
                allocate(rkx(nx))
                allocate(hrkx(nx))
                allocate(rky(ny))
                allocate(hrky(ny))
                allocate(rkz(nz))
                allocate(xtrig(2*nx))
                allocate(ytrig(2*ny))
                allocate(ztrig(2*nz))
                allocate(skx(nx))
                allocate(sky(ny))
                allocate(skz(0:nz))
            endif

            nwx = nx / 2
            nwy = ny / 2

            !----------------------------------------------------------
            ! Set up FFTs:
            ! Initialise FFTs and wavenumber arrays:
            call init2dfft(nx, ny, extent(1), extent(2), xfactors, yfactors, xtrig, ytrig, hrkx, hrky)
            call initfft(nz, zfactors, ztrig)

            !Define x wavenumbers:
            rkx(1) = zero
            do kx = 1, nwx-1
                rkx(kx+1)    = hrkx(2*kx)
                rkx(nx+1-kx) = hrkx(2*kx)
            enddo
            rkx(nwx+1) = hrkx(nx)

            !Define y wavenumbers:
            rky(1) = zero
            do ky = 1, nwy-1
                rky(ky+1)    = hrky(2*ky)
                rky(ny+1-ky) = hrky(2*ky)
            enddo
            rky(nwy+1) = hrky(ny)

            !Define z wavenumbers:
            call init_deriv(nz, extent(3), rkz)


            !----------------------------------------------------------
            !Define Hou and Li filter:
            kxmaxi = one / maxval(rkx)
            skx = -36.d0 * (kxmaxi * rkx) ** 36
            kymaxi = one/maxval(rky)
            sky = -36.d0 * (kymaxi * rky) ** 36
            kzmaxi = one / maxval(rkz)
            skz(0) = zero
            skz(1:nz)=-36.d0 * (kzmaxi * rkz) ** 36
            do ky = 1, ny
                do kx = 1, nx
                    do kz = 0, nz
                        filt(kz, kx, ky) = dexp(skx(kx) + sky(ky) + skz(kz))
                    enddo
                enddo
            enddo

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
                    call dct(1, nzval, fs(:, kx, ky), ztrig, zfactors)
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
                    call dct(1, nzval, fp(:, iy, ix), ztrig, zfactors)
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
            nyval = size(fp, 2)
            nxval = size(fp, 3)

            !Carry out z FFT for each kx and ky:
            do ky = 1, nyval
                do kx = 1, nxval
                    fs(:, kx, ky) = fp(:, kx, ky)
                    call dct(1, nzval, fs(:, kx, ky), ztrig, zfactors)
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
                    call dct(1, nzval, fp(:, ix, iy), ztrig, zfactors)
                enddo
            enddo
        end subroutine

end module sta3dfft
