module inversion_utils
    use constants
    use parameters, only : nx, ny, nz, dx, dxi, extent, ncelli
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

        subroutine init_hyperdiffusion(te, en)
            double precision, intent(in) :: te ! kinetic + potential energy
            double precision, intent(in) :: en ! enstrophy
            double precision             :: char_len, char_vel, kolm_len, visc
            integer                      :: kx, ky, kz

            allocate(hdis(0:nz, 0:nx-1, 0:ny-1))

            !Check if FFT is initialised
            if (.not. is_fft_initialised) then
                print *, "Error: FFT not initialised!"
                stop
            endif

            char_len = dsqrt(te / en)
            char_vel = dsqrt(two * te * ncelli) ! ncelli = 1 / (nx * ny * nz)
            kolm_len = viscosity%kolm_fac * maxval(dx)
            visc = char_vel * char_len ** (2 * viscosity%nnu - 1) * (kolm_len / char_len) ** f43
            write(*,'(a,1p,e14.7)') ' (hyper)viscosity nu = ', visc

            !Define spectral dissipation operator:
            do ky = 0, ny-1
                do kx = 0, nx-1
                    hdis(0,  kx, ky) = visc * k2l2(kx, ky) ** viscosity%nnu
                    hdis(nz, kx, ky) = hdis(0,  kx, ky)
                    do kz = 1, nz-1
                        hdis(kz, kx, ky) = visc * (k2l2(kx, ky) + rkz(kz) ** 2) ** viscosity%nnu
                    enddo
                enddo
            enddo

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

end module inversion_utils
