module zops
    use constants, only : zero, f12, one
    use mpi_layout
    use parameters, only : nz       &
                         , extent   &
                         , hl, hli  &
                         , center   &
                         , upper    &
                         , lower
    use cheby ! Import Chebyshev module to set up various matrices needed below
    use sta3dfft, only : initialise_fft
    use inversion_utils, only : k2l2
!     use options, only : zfiltering
    implicit none

    private

    double precision, allocatable :: d1z(:, :), d2z(:, :) &
                                   , zcheb(:)             & ! Chebyshev grid points
                                   , zg(:)                & ! Chebyshev domains points
                                   , zccw(:)              & ! Clenshaw-Curtis weights
                                   , zfilt(:)

    integer :: nxym1 = 0
    logical :: l_zops_initialised = .false.

    public :: init_zops             &
            , finalise_zops         &
            , zderiv                &
            , zzderiv               &
            , zinteg                &
            , vertvel               &
            , zcheb                 &
            , zg                    &
            , zccw                  &
            , apply_zfilter         &
            , d1z, d2z

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine init_zops
            double precision :: fdz1, fdz2, rkmax
            integer          :: iz

            if (l_zops_initialised) then
                return
            endif

            l_zops_initialised = .true.

            !------------------------------------------------------------------
            ! Ensure FFT module is initialised:
            ! (this call does nothing if already initialised)
            call initialise_fft(extent)

            !------------------------------------------------------------------
            ! Allocate arrays:
            allocate(d1z(0:nz, 0:nz))
            allocate(d2z(0:nz, 0:nz))
            allocate(zcheb(0:nz))
            allocate(zg(0:nz))
            allocate(zccw(0:nz))
            allocate(zfilt(0:nz))

            nxym1 = box%size(1) * box%size(2) - 1

            ! Note: hli = two / extent
            fdz1 = - hli(3)
            fdz2 = fdz1 * fdz1

            !-----------------------------------------------------------------
            ! Get Chebyshev points & 1st & 2nd order differentiation matrices:
            call init_cheby(nz, zcheb, d1z, d2z)

            zg = center(3) - hl(3) * zcheb

            ! Scale d1z & d2z for the actual z limits:
            d1z = fdz1 * d1z
            d2z = fdz2 * d2z

            ! Get Clenshaw-Curtis weights:
            call clencurt(nz, zccw)

!             !------------------------------------------------------------------
!             ! Dembenek filter:
!             rkmax = zfiltering%kmax * dble(nz)
!             do iz = 0, nz
!                 zfilt(iz) = dembenek_filter(iz,                 &
!                                             rkmax,              &
!                                             zfiltering%alpha,   &
!                                             zfiltering%beta)
!             enddo

        end subroutine init_zops

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        function dembenek_filter(k, rkmax, alpha, beta) result(res)
            integer,          intent(in) :: k
            double precision, intent(in) :: alpha, beta, rkmax
            double precision             :: res, x

            x = dble(k)/rkmax
            res = dexp(-alpha*x**beta)
        end function dembenek_filter

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine finalise_zops

            if (l_zops_initialised) then
                deallocate(d1z)
                deallocate(d2z)
                deallocate(zcheb)
                deallocate(zg)
                deallocate(zfilt)
                l_zops_initialised = .false.
            endif

        end subroutine finalise_zops

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Calculates g = df/dz
        subroutine zderiv(f, g)
            double precision, intent(in)  :: f(0:nz, 0:nxym1)
            double precision, intent(out) :: g(0:nz, 0:nxym1)

            g = matmul(d1z, f)

        end subroutine zderiv

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Calculates g = d^2f/dz^2
        subroutine zzderiv(f, g)
            double precision, intent(in)  :: f(0:nz, 0:nxym1)
            double precision, intent(out) :: g(0:nz, 0:nxym1)

            g = matmul(d2z, f)

        end subroutine zzderiv

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Solves dg/dz = f for g, either with g(0) at z = zmin when
        ! noavg = .false. or with the average of g over z equal to 0
        ! when noavg = .true.
        ! f & g are 1D arrays over z,
        ! *** Uses dgesv from LAPACK/BLAS ***
        subroutine zinteg(f, g, noavg)
            double precision, intent(in)  :: f(0:nz)
            double precision, intent(out) :: g(0:nz)
            logical,          intent(in)  :: noavg
            double precision              :: dmat(0:nz-1, 0:nz-1), h(0:nz), gavg
            integer                       :: ipiv(0:nz-1), info

            !-----------------------------------------------------
            ! Integrate starting from g = 0 at z = zmin:
            g = f
            dmat = d1z(0:nz-1, 0:nz-1)
            call dgesv(nz, 1, dmat, nz, ipiv, g(0:nz-1), nz, info)
            g(nz) = zero

            if (.not. noavg) then
                return
            endif

            !-----------------------------------------------------
            ! Remove average of g:
            h = g
            dmat = d1z(0:nz-1, 0:nz-1)
            call dgesv(nz, 1, dmat, nz, ipiv, h(0:nz-1), nz, info)
            gavg = h(0) / extent(3)

            g = g + gavg

        end subroutine zinteg

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Solves (d^2/dz^2 - K^2)[ws] = S in semi-spectral space where
        ! K^2 = k^2 + l^2 is the squared horizontal wavenumber and
        ! where ws initially contains the source S (this is overwritten
        ! by the solution).
        ! *** Uses dgesv from LAPACK/BLAS ***
        subroutine vertvel(ws)
            double precision, intent(inout) :: ws(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision                :: dmat(nz-1, nz-1), sol(nz-1)
            integer                         :: ipiv(nz-1), info
            integer                         :: kx, ky, iz

            !-----------------------------------------------------------------
            ! Loop over horizontal wavenumbers and solve linear system:
            do kx = box%lo(1), box%hi(1)
                do ky = box%lo(2), box%hi(2)
                    ! Inner part of d^2/dz^2 matrix:
                    dmat = d2z(1:nz-1, 1:nz-1)

                    ! Remove K^2 down the diagonal:
                    do iz = 1, nz-1
                        dmat(iz, iz) = dmat(iz, iz) - k2l2(ky, kx)
                    enddo

                    ! Linear solve with LAPACK:
                    sol = ws(1:nz-1, ky, kx)
                    call dgesv(nz-1, 1, dmat, nz-1, ipiv, sol, nz-1, info)
                    ws(1:nz-1, ky, kx) = sol

                    ! Add zero boundary values:
                    ws(0,  ky, kx) = zero
                    ws(nz, ky, kx) = zero
                enddo
            enddo

        end subroutine vertvel

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine apply_zfilter(fs)
            double precision, intent(inout) :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision                :: coeffs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision                :: err_e(box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision                :: err_o(box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            integer                         :: iz

            ! get Chebyshev coefficients
            call get_cheb_poly(fs, coeffs)

            ! apply filter on coefficients
            do iz = 0, nz
                coeffs(iz, :, :) = zfilt(iz) * coeffs(iz, :, :)
            enddo

            ! Boundary-Preserving Filter:
            err_e = coeffs(0, :, :)
            err_o = coeffs(1, :, :)

            do iz = 1, nz/2
                err_e  = err_e +  coeffs(2*iz, :, :)
            enddo

            do iz = 1, nz/2-1
                err_o  = err_o +  coeffs(2*iz+1, :, :)
            enddo

            ! Adjust mean value and linear slope to insure 0 BC's
            coeffs(0, :, :) = coeffs(0, :, :) - err_e
            coeffs(1, :, :) = coeffs(1, :, :) - err_o

            ! Return filtered field with 0 bc's
            call cheb_eval(coeffs, fs)

        end subroutine apply_zfilter

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Input:
        ! fs - a vector of length N+1 containing function values at Chebyshev nodes in [-1, 1]
        ! Output:
        ! c - a vector of length N+1 containing the coefficients of the Chebyshev polynomials
        subroutine get_cheb_poly(fs, c)
            double precision, intent(in)  :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision, intent(out) :: c(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            integer                       :: kx, ky

            do kx = box%lo(1), box%hi(1)
                do ky = box%lo(2), box%hi(2)
                    call cheb_poly(nz, fs(:, ky, kx), c(:, ky, kx))
                enddo
            enddo

        end subroutine get_cheb_poly

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Input:
        ! c - a vector of length N+1 containing the coefficients of the Chebyshev polynomials
        ! Output:
        ! fs - a vector of length N+1 containing the values of f(x) at the points in y
        subroutine cheb_eval(c, fs)
            double precision, intent(in)  :: c(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision, intent(out) :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            integer                       :: kx, ky

            do kx = box%lo(1), box%hi(1)
                do ky = box%lo(2), box%hi(2)
                    call cheb_fun(nz, c(:, ky, kx), fs(:, ky, kx))
                enddo
            enddo

        end subroutine cheb_eval

end module zops
