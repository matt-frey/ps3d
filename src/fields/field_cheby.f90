module field_cheby
    use constants, only : zero, f12, one
    use parameters, only : nz       &
                         , extent   &
                         , hl, hli  &
                         , center   &
                         , upper    &
                         , lower    &
                         , dx
    use cheby ! Import Chebyshev module to set up various matrices needed below
    use field_layout
    use mpi_layout, only : box
    use inversion_utils, only : phim, phip, filt, k2l2
    use sta3dfft, only : fftxyp2s, fftxys2p, initialise_fft
    use mpi_collectives, only : mpi_blocking_reduce
    implicit none

    type, extends (flayout_t) :: field_cheby_t

        private

        integer :: nxym1 = 0
        logical :: l_initialised = .false.

        double precision, allocatable :: d1z(:, :), d2z(:, :) &
                                       , zcheb(:)             & ! Chebyshev grid points
                                       , zccw(:)              & ! Clenshaw-Curtis weights
                                       , zfilt(:)

    contains
        procedure :: initialise
        procedure :: finalise

        procedure :: get_z_axis

        ! Field decompositions:
        procedure :: decompose_physical
        procedure :: combine_physical
        procedure :: decompose_semi_spectral
        procedure :: combine_semi_spectral

        ! Field diagnostics:
        procedure :: get_local_sum

        ! Field operations:
        procedure :: diffz
        procedure :: calc_decomposed_mean
        procedure :: adjust_decomposed_mean
        procedure :: apply_filter

        ! Private procedures:
        procedure, private :: get_cheb_poly
        procedure, private :: cheb_eval

        procedure :: vertvel
        procedure :: zderiv
        procedure :: zinteg

    end type field_cheby_t

contains

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine initialise(this)
        class (field_cheby_t), intent(inout)  :: this
        double precision :: fdz1, fdz2, rkmax
        integer          :: iz

        if (this%l_initialised) then
            return
        endif

        this%l_initialised = .true.

        !------------------------------------------------------------------
        ! Ensure FFT module is initialised:
        ! (this call does nothing if already initialised)
        call initialise_fft(extent)

        !------------------------------------------------------------------
        ! Allocate arrays:
        allocate(this%d1z(0:nz, 0:nz))
        allocate(this%d2z(0:nz, 0:nz))
        allocate(this%zcheb(0:nz))
        allocate(this%zccw(0:nz))
        allocate(this%zfilt(0:nz))

        this%nxym1 = box%size(1) * box%size(2) - 1

        ! Note: hli = two / extent
        fdz1 = - hli(3)
        fdz2 = fdz1 * fdz1

        !-----------------------------------------------------------------
        ! Get Chebyshev points & 1st & 2nd order differentiation matrices:
        call init_cheby(nz, this%zcheb, this%d1z, this%d2z)

        ! Scale d1z & d2z for the actual z limits:
        this%d1z = fdz1 * this%d1z
        this%d2z = fdz2 * this%d2z

        ! Get Clenshaw-Curtis weights:
        call clencurt(nz, this%zccw)

!             !------------------------------------------------------------------
!             ! Dembenek filter:
!             rkmax = zfiltering%kmax * dble(nz)
!             do iz = 0, nz
!                 zfilt(iz) = dembenek_filter(iz,                 &
!                                             rkmax,              &
!                                             zfiltering%alpha,   &
!                                             zfiltering%beta)
!             enddo


    end subroutine

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine finalise(this)
        class (field_cheby_t), intent(inout)  :: this

        if (this%l_initialised) then
            deallocate(this%d1z)
            deallocate(this%d2z)
            deallocate(this%zcheb)
            deallocate(this%zfilt)
            this%l_initialised = .false.
        endif

    end subroutine finalise

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_z_axis(this)
        class (field_cheby_t), intent(in)  :: this
        double precision :: get_z_axis(0:nz)

        get_z_axis = center(3) - hl(3) * this%zcheb

    end function get_z_axis

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! fc  - complete field (physical space)
    ! sf  - full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
    ! cfc - copy of complete field (physical space)
    subroutine decompose_physical(this, fc, sf)
        class (field_cheby_t), intent(in)  :: this
        double precision,      intent(in)  :: fc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision,      intent(out) :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))

        call fftxyp2s(fc, sf)

        call this%decompose_semi_spectral(sf)

    end subroutine decompose_physical

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! in : complete field (semi-spectral space)
    ! out: full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
    subroutine decompose_semi_spectral(this, sfc)
        class (field_cheby_t), intent(in)    :: this
        double precision,      intent(inout) :: sfc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        integer                              :: iz

        ! subtract harmonic part
        !$omp parallel do
        do iz = 1, nz-1
            sfc(iz, :, :) = sfc(iz, :, :) - (sfc(0, :, :) * phim(iz, :, :) + sfc(nz, :, :) * phip(iz, :, :))
        enddo
        !$omp end parallel do

    end subroutine decompose_semi_spectral

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! sf  - full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
    ! fc  - complete field (physical space)
    ! sfc - complete field (semi-spectral space)
    subroutine combine_physical(this, sf, fc)
        class (field_cheby_t), intent(in)  :: this
        double precision,      intent(in)  :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision,      intent(out) :: fc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision                   :: sfc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))

        sfc = sf

        call this%combine_semi_spectral(sfc)

        ! transform to physical space as fc:
        call fftxys2p(sfc, fc)

    end subroutine combine_physical

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! in : full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
    ! out: complete field (semi-spectral space)
    subroutine combine_semi_spectral(this, sf)
        class (field_cheby_t), intent(in)    :: this
        double precision,      intent(inout) :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        integer                              :: iz

        ! add harmonic part to sfc:
        !$omp parallel do
        do iz = 1, nz-1
            sf(iz, :, :) = sf(iz, :, :) + sf(0, :, :) * phim(iz, :, :) + sf(nz, :, :) * phip(iz, :, :)
        enddo
        !$omp end parallel do

    end subroutine combine_semi_spectral

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_local_sum(this, ff) result(res)
        class (field_cheby_t), intent(in) :: this
        double precision,      intent(in) :: ff(box%lo(3):box%hi(3), &
                                                box%lo(2):box%hi(2), &
                                                box%lo(1):box%hi(1))
        double precision                  :: res
        integer                           :: iz

        res = zero
        do iz = box%lo(3), box%hi(3)
            res = res + this%zccw(iz) * sum(ff(iz, :, :))
        enddo

        res = res * f12 * dble(nz)

    end function get_local_sum

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine diffz(this, fs, ds)
        class (field_cheby_t), intent(in)  :: this
        double precision,      intent(in)  :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision,      intent(out) :: ds(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))

        call this%zderiv(fs, ds)

    end subroutine diffz

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! This is only calculated on the MPI rank having kx = ky = 0
    function calc_decomposed_mean(this, fs) result(savg)
        class (field_cheby_t), intent(in) :: this
        double precision,      intent(in) :: fs(box%lo(3):box%hi(3), &
                                                box%lo(2):box%hi(2), &
                                                box%lo(1):box%hi(1))
        double precision                   :: savg
        integer                            :: iz

        if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
            savg = zero
            do iz = box%lo(3), box%hi(3)
                savg = savg + this%zccw(iz) * fs(iz, 0, 0)
            enddo
            ! The factor f12 * extent(3) comes from the mapping [-1, 1] to [a, b]
            ! where the Chebyshev points are given in [-1, 1]
            ! z = (b-a) / 2 * t + (a+b)/2 for [a, b] --> dz = (b-a) / 2 * dt
            ! However, we must divide by extent(3) again in order to get the vertical domain-average.
            ! Hence, we only scale by f12.
            savg = savg * f12
        endif

    end function calc_decomposed_mean

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! This is only calculated on the MPI rank having kx = ky = 0
    subroutine adjust_decomposed_mean(this, fs, avg)
        class (field_cheby_t), intent(in)    :: this
        double precision,      intent(inout) :: fs(box%lo(3):box%hi(3), &
                                                    box%lo(2):box%hi(2), &
                                                    box%lo(1):box%hi(1))
        double precision,      intent(in)    :: avg
        double precision                     :: savg

        savg = this%calc_decomposed_mean(fs)

        if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
            fs(: , 0, 0) = fs(:, 0, 0) + avg - savg
        endif

    end subroutine adjust_decomposed_mean

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine apply_filter(this, fs)
        class (field_cheby_t), intent(in)    :: this
        double precision,      intent(inout) :: fs(box%lo(3):box%hi(3), &
                                                    box%lo(2):box%hi(2), &
                                                    box%lo(1):box%hi(1))
        double precision                     :: coeffs(0:nz,                &
                                                        box%lo(2):box%hi(2), &
                                                        box%lo(1):box%hi(1))
        double precision                     :: err_e(box%lo(2):box%hi(2),  &
                                                        box%lo(1):box%hi(1))
        double precision                     :: err_o(box%lo(2):box%hi(2),  &
                                                        box%lo(1):box%hi(1))
        integer                              :: iz
        double precision                     :: fstop(box%lo(2):box%hi(2), &
                                                        box%lo(1):box%hi(1))
        double precision                     :: fsbot(box%lo(2):box%hi(2), &
                                                        box%lo(1):box%hi(1))

        ! Temporarily store the surfaces
        fsbot = fs(0,  :, :)
        fstop = fs(nz, :, :)

        ! Ensure surface are zero before applying filter in Chebyshev space
        fs(0,  :, :) = zero
        fs(nz, :, :) = zero

        ! Get Chebyshev coefficients
        call this%get_cheb_poly(fs, coeffs)

        ! Apply filter on coefficients
        do iz = 0, nz
            coeffs(iz, :, :) = filt(iz, :, :) * coeffs(iz, :, :)
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
        call this%cheb_eval(coeffs, fs)

        ! Restore filtered surfaces
        fs(0,  :, :) = filt(0,  :, :) * fsbot
        fs(nz, :, :) = filt(nz, :, :) * fstop

    end subroutine apply_filter

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Input:
    ! fs - a vector of length N+1 containing function values at Chebyshev nodes in [-1, 1]
    ! Output:
    ! c - a vector of length N+1 containing the coefficients of the Chebyshev polynomials
    subroutine get_cheb_poly(this, fs, c)
        class (field_cheby_t), intent(in)  :: this
        double precision,      intent(in)  :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision,      intent(out) :: c(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        integer                            :: kx, ky

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
    subroutine cheb_eval(this, c, fs)
        class (field_cheby_t), intent(in)  :: this
        double precision,      intent(in)  :: c(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision,      intent(out) :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        integer                            :: kx, ky

        do kx = box%lo(1), box%hi(1)
            do ky = box%lo(2), box%hi(2)
                call cheb_fun(nz, c(:, ky, kx), fs(:, ky, kx))
            enddo
        enddo

    end subroutine cheb_eval

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Solves dg/dz = f for g, either with g(0) at z = zmin when
    ! noavg = .false. or with the average of g over z equal to 0
    ! when noavg = .true.
    ! f & g are 1D arrays over z,
    ! *** Uses dgesv from LAPACK/BLAS ***
    subroutine zinteg(this, f, g, noavg)
        class (field_cheby_t), intent(in)  :: this
        double precision, intent(in)  :: f(0:nz)
        double precision, intent(out) :: g(0:nz)
        logical,          intent(in)  :: noavg
        double precision              :: dmat(0:nz-1, 0:nz-1), h(0:nz), gavg
        integer                       :: ipiv(0:nz-1), info

        !-----------------------------------------------------
        ! Integrate starting from g = 0 at z = zmin:
        g = f
        dmat = this%d1z(0:nz-1, 0:nz-1)
        call dgesv(nz, 1, dmat, nz, ipiv, g(0:nz-1), nz, info)
        g(nz) = zero

        if (.not. noavg) then
            return
        endif

        !-----------------------------------------------------
        ! Remove average of g:
        h = g
        dmat = this%d1z(0:nz-1, 0:nz-1)
        call dgesv(nz, 1, dmat, nz, ipiv, h(0:nz-1), nz, info)
        gavg = h(0) / extent(3)

        g = g + gavg

    end subroutine zinteg

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Solves (d^2/dz^2 - K^2)[ds] = S in semi-spectral space where
    ! K^2 = k^2 + l^2 is the squared horizontal wavenumber and
    ! where ds initially contains the source S (this is overwritten
    ! by the solution).
    ! *** Uses dgesv from LAPACK/BLAS ***
    subroutine vertvel(this, ds, es)
        class (field_cheby_t), intent(in)  :: this
        double precision, intent(inout) :: ds(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision, intent(out)   :: es(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision                :: dmat(nz-1, nz-1), sol(nz-1)
        integer                         :: ipiv(nz-1), info
        integer                         :: kx, ky, iz

        !-----------------------------------------------------------------
        ! Loop over horizontal wavenumbers and solve linear system:
        do kx = box%lo(1), box%hi(1)
            do ky = box%lo(2), box%hi(2)
                ! Inner part of d^2/dz^2 matrix:
                dmat = this%d2z(1:nz-1, 1:nz-1)

                ! Remove K^2 down the diagonal:
                do iz = 1, nz-1
                    dmat(iz, iz) = dmat(iz, iz) - k2l2(ky, kx)
                enddo

                ! Linear solve with LAPACK:
                sol = ds(1:nz-1, ky, kx)
                call dgesv(nz-1, 1, dmat, nz-1, ipiv, sol, nz-1, info)
                ds(1:nz-1, ky, kx) = sol

                ! Add zero boundary values:
                ds(0,  ky, kx) = zero
                ds(nz, ky, kx) = zero
            enddo
        enddo

        ! Calculate z-derivative of vertical velocity:
        call this%diffz(ds, es)

    end subroutine vertvel

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        function dembenek_filter(k, rkmax, alpha, beta) result(res)
            integer,          intent(in) :: k
            double precision, intent(in) :: alpha, beta, rkmax
            double precision             :: res, x

            x = dble(k)/rkmax
            res = dexp(-alpha*x**beta)
        end function dembenek_filter

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Calculates g = df/dz
        subroutine zderiv(this, f, g)
            class (field_cheby_t), intent(in)  :: this
            double precision, intent(in)  :: f(0:nz, 0:this%nxym1)
            double precision, intent(out) :: g(0:nz, 0:this%nxym1)

            g = matmul(this%d1z, f)

        end subroutine zderiv

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Calculates g = d^2f/dz^2
        subroutine zzderiv(this, f, g)
            class (field_cheby_t), intent(in)  :: this
            double precision, intent(in)  :: f(0:nz, 0:this%nxym1)
            double precision, intent(out) :: g(0:nz, 0:this%nxym1)

            g = matmul(this%d2z, f)

        end subroutine zzderiv

end module field_cheby
