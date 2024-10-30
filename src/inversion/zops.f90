module zops
    use constants, only : f12
    use mpi_layout
    use parameters, only : nz       &
                         , extent   &
                         , hl, hli  &
                         , center   &
                         , upper    &
                         , lower
    use cheby ! Import Chebyshev module to set up various matrices needed below
    use sta3dfft, only : initialise_fft &
                       , k2l2           &
                       , zfactors       &
                       , ztrig
    use stafft, only : dct
    implicit none

    private

    double precision, allocatable :: d1z(:, :), d2z(:, :) &
                                   , zcheb(:)             & ! Chebyshev grid points
                                   , zg(:)                & ! Chebyshev domains points
                                   , zccw(:)              & ! Clenshaw-Curtis weights
                                   , zfilt(:)             &
                                   , tm(:, :)

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
            , d1z, d2z              &
            , cheb_poly

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine init_zops
            double precision :: fdz1, fdz2, c, cutoff
            integer          :: iz, p, i, j

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
            allocate(tm(0:nz, 0:nz))
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

            !------------------------------------------------------------------
            ! Butterworth filter:
            cutoff = 0.9d0
            p = 4 ! filter order
            do iz = 0, nz
                c = dble(iz) / dble(nz)
                zfilt(iz) = one / sqrt(one + (c / cutoff)**(2*p))
            enddo

            ! Construct the T matrix using the cosine formulation
            do i = 0, nz
                do j = 0, nz
                    tm(i, j) = cos(dble(j) * acos(zcheb(i)))
                enddo
            enddo

        end subroutine init_zops

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine finalise_zops

            if (l_zops_initialised) then
                deallocate(d1z)
                deallocate(d2z)
                deallocate(zcheb)
                deallocate(zg)
                deallocate(tm)
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
            integer                         :: iz
            double precision                :: fsbot(box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision                :: fstop(box%lo(2):box%hi(2), box%lo(1):box%hi(1))

            ! get Chebyshev coefficients
            call cheb_poly(fs, coeffs)

            ! apply filter on coefficients
            do iz = 0, nz
                coeffs(iz, :, :) = zfilt(iz) * coeffs(iz, :, :)
            enddo

            fsbot = fs(0,  :, :)
            fstop = fs(nz, :, :)
            ! only works on interior points: fix the endpoints afterwards
            call cheb_eval(fs, coeffs)
            fs(0,  :, :) = fsbot
            fs(nz, :, :) = fstop

        end subroutine apply_zfilter

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Input:
        ! fs - a vector of length N+1 containing function values at Chebyshev nodes in [-1, 1]
        ! Output:
        ! c - a vector of length N+1 containing the coefficients of the Chebyshev polynomials
        subroutine cheb_poly(fs, c)
            double precision, intent(in)  :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision, intent(out) :: c(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            integer                       :: kx, ky
            double precision              :: valsUnitDisc(0:2*nz)

            ! Initialize the coefficients vector
            c = fs

            do kx = box%lo(1), box%hi(1)
                do ky = box%lo(2), box%hi(2)

                    ! String out values into values on equally spaced theta grid:
                    ! valsUnitDisc = [fvals ; flipud(fvals(2:end-1))];
                    valsUnitDisc(0:nz) = fs(0:nz, ky, kx)
                    valsUnitDisc(nz+1:) = fs(2:nz-1, ky, kx)
                    call flipud(valsUnitDisc(nz+1:))



!                     c(:, ky, kx) =

                    ! Forward cosine transform:
                    call dct(1, 2*nz, valsUnitDisc(0:2*nz), ztrig, zfactors)
                    c(0:nz, ky, kx) = valsUnitDisc(0:nz)

                    ! Get Chebyshev coefficients:
                    c(:,  ky, kx) =       c(:,  ky, kx) / dble(nz)
                    c(0,  ky, kx) = f12 * c(0,  ky, kx)
                    c(nz, ky, kx) = f12 * c(nz, ky, kx)
                enddo
            enddo

        end subroutine cheb_poly

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Input:
        ! c - a vector of length N+1 containing the coefficients of the Chebyshev polynomials
        ! Output:
        ! f_values - a vector of length M containing the values of f(x) at the points in y
        subroutine cheb_eval(fs, c)
            double precision, intent(inout) :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision, intent(in)    :: c(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            integer                         :: kx, ky

            ! Compute fs = tm * c using matrix multiplication
            do kx = box%lo(1), box%hi(1)
                do ky = box%lo(2), box%hi(2)
                    fs(:, ky, kx) = matmul(tm, c(:, ky, kx))
                enddo
            enddo

        end subroutine cheb_eval


!         !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!         ! Given N+1 coefficients of the Chebyshev series
!         ! Returns f(x) evaluated at the N+1 chebyshev nodes in [-1,1];
!         subroutine cheb_fun(c, fs)
!             double precision, intent(in)  :: c(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
!             double precision, intent(out) :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
!             integer                       :: kx, ky
!
!
!             do kx = box%lo(1), box%hi(1)
!                 do ky = box%lo(2), box%hi(2)
!                     fs(0,      ky, kx) = two * c(0,      ky, kx)
!                     fs(1:nz-1, ky, kx) =       c(1:nz-1, ky, kx)
!                     fs(nz,     ky, kx) = two * c(nz,     ky, kx)
!                     fs(:,      ky, kx) = nz  * fs(:,     ky, kx)
!
!                     !
!                     ! FourierCoeffs = [ChebC; flipud(ChebC(2:n))];
!                     !
!
!
!                     ! Backward (i.e. inverse) cosine transform:
!                     call dct(1, nz, fs(0:nz, ky, kx), ztrig, zfactors)
!
!                     ! Get Chebyshev coefficients:
!                     ! fvals = fvals(1:n+1);
!                 enddo
!             enddo
!
!         end subroutine cheb_fun

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Flip the vertical direction
        subroutine flipud(fs)
            double precision, intent(inout) :: fs(0:nz)
            double precision                :: gs
            integer                         :: iz, n

            n = int(nz / 2)
            do iz = 0, n
                gs = fs(iz)
                fs(iz) = fs(nz - iz)
                fs(nz - iz) = gs
            enddo

        end subroutine flipud

end module zops
