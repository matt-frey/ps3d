module sta3dfft
    use constants, only : zero, f12
    use stafft
    use sta2dfft
    use deriv1d
    implicit none

    ! Wavenumbers::
    double precision, allocatable :: rkx(:), hrkx(:), rky(:), hrky(:), rkz(:)

    ! Quantities needed in FFTs:
    double precision, allocatable :: xtrig(:), ytrig(:), ztrig(:)
    integer                       :: xfactors(5), yfactors(5), zfactors(5)


    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        subroutine init3dfft(nx, ny, nz, extent)
            integer,          intent(in) :: nx, ny, nz
            double precision, intent(in) :: extent(3)
            integer                      :: nwx, nwy, kx, ky, kz
            integer                      :: kxc, kyc

            if (.not. allocated(rkx)) then
                allocate(rkx(0:nx-1))
                allocate(hrkx(nx))
                allocate(rky(0:ny-1))
                allocate(hrky(ny))
                allocate(rkz(nz))
                allocate(xtrig(2*nx))
                allocate(ytrig(2*ny))
                allocate(ztrig(2*nz))
            endif

            nwx = nx / 2
            nwy = ny / 2

            !----------------------------------------------------------
            ! Set up FFTs:
            ! Initialise FFTs and wavenumber arrays:
            call init2dfft(nx, ny, extent(1), extent(2), xfactors, yfactors, xtrig, ytrig, hrkx, hrky)
            call initfft(nz, zfactors, ztrig)

            !Define x wavenumbers:
            rkx(0) = zero
            do kx = 0, nwx-1
                kxc = nx - kx
                rkx(kx)  = hrkx(2*kx)
                rkx(kxc) = hrkx(2*kx)
            enddo
            rkx(nwx) = hrkx(nx)

            !Define y wavenumbers:
            rky(0) = zero
            do ky = 1, nwy-1
                kyc = ny - ky
                rky(ky)  = hrky(2*ky)
                rky(kyc) = hrky(2*ky)
            enddo
            rky(nwy) = hrky(ny)

            !Define z wavenumbers:
            call init_deriv(nz, extent(3), rkz)

        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Computes a 3D FFT of an array fp in physical space and
        !returns the result as fs in spectral space.  It is assumed that
        !fp is generally non-zero at the z boundaries (so a cosine
        !transform is used in z).
        !*** fp is destroyed upon exit ***
        subroutine fftczp2s(fp, fs, nx, ny, nz)
            integer,          intent(in)    :: nx, ny, nz
            double precision, intent(inout) :: fp(0:nz, ny, nx)  !Physical
            double precision, intent(out)   :: fs(0:nz, nx, ny)  !Spectral
            integer                         :: kx, ky, iy

            !Carry out a full x transform first:
            call forfft((nz+1)*ny, nx, fp, xtrig, xfactors)

            !Transpose array:
            do kx = 1, nx
                do iy = 1, ny
                    fs(:, kx, iy) = fp(:, iy, kx)
                enddo
            enddo

            !Carry out a full y transform on transposed array:
            call forfft((nz+1)*nx, ny, fs, ytrig, yfactors)

            !Carry out z FFT for each kx and ky:
            do ky = 1, ny
                do kx = 1, nx
                    call dct(1, nz, fs(:, kx, ky), ztrig, zfactors)
                enddo
            enddo
        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Computes an *inverse* 3D FFT of an array fs in spectral space and
        !returns the result as fp in physical space.  It is assumed that
        !fp is generally non-zero at the z boundaries (so a cosine
        !transform is used in z).
        !*** fs is destroyed upon exit ***
        subroutine fftczs2p(fs, fp, nx, ny, nz)
            integer,          intent(in)    :: nx, ny, nz
            double precision, intent(inout) :: fs(0:nz, nx, ny)  !Spectral
            double precision, intent(out)   :: fp(0:nz, ny, nx)  !Physical
            integer                         :: ix, iy, kx

            !Carry out a full inverse y transform first:
            call revfft((nz+1)*nx,ny,fs,ytrig,yfactors)

            !Transpose array:
            do kx = 1, nx
                do iy = 1, ny
                    fp(:, iy, kx) = fs(:, kx, iy)
                enddo
            enddo

            !Carry out a full inverse x transform:
            call revfft((nz+1)*ny, nx, fp, xtrig, xfactors)

            !Carry out z FFT for each ix and iy:
            do ix = 1, nx
                do iy = 1, ny
                    call dct(1, nz, fp(:, iy, ix), ztrig, zfactors)
                enddo
            enddo
        end subroutine

end module sta3dfft
