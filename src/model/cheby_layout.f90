module cheby_layout
    use cheby ! Import Chebyshev module to set up various matrices needed below
    use constants, only : zero, f12, f23, one
    use parameters, only : nz       &
                         , hl, hli  &
                         , center   &
                         , upper    &
                         , lower    &
                         , extent   &
                         , dx
    use field_layout
    use mpi_layout, only : box
    use sta3dfft, only : is_fft_initialised, fftxyp2s, fftxys2p
    use mpi_collectives, only : mpi_blocking_reduce
    use mpi_utils, only : mpi_stop
    use sta3dfft, only : k2l2, rkx, rky, rkz
    implicit none

    type, extends (layout_t) :: cheby_layout_t

        integer :: nxym1 = 0
        logical :: l_initialised = .false.

        double precision, allocatable :: d1z(:, :), d2z(:, :) &
                                       , zcheb(:)             & ! Chebyshev grid points
                                       , zccw(:)                ! Clenshaw-Curtis weights
!                                        , zfilt(:)

        double precision, allocatable ::  eye(:, :), D2(:, :)
        double precision, allocatable ::  eyeNF(:, :), D2NF(:, :)

        ! Filter for the Chebyshev cofficients:
        double precision, allocatable :: zfilt(:, :, :)

        ! Filter for the surfaces:
        double precision, allocatable :: filt(:, :)

    contains

        procedure :: initialise
        procedure :: finalise
        procedure :: finalise_filter

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

        procedure :: vertvel
        procedure :: zinteg
        procedure :: zdiffuse
        procedure :: zdiffNF

        ! Filters:
        procedure :: init_filter
        procedure :: apply_filter
        procedure :: apply_hfilter
        procedure, private :: init_hou_and_li_filter
        procedure, private :: init_23rd_rule_filter
        procedure, private :: init_no_filter
        procedure, private :: get_cheb_poly
        procedure, private :: cheb_eval


        ! Routines only available in this class
        procedure :: zderiv
        procedure :: zzderiv

        procedure :: decomposed_diffz
        procedure :: semi_spectral_diffz

    end type cheby_layout_t

contains

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine initialise(this, filter_method)
        class (cheby_layout_t), intent(inout) :: this
        character(*), optional, intent(in)    :: filter_method
        double precision                      :: fdz1, fdz2!, rkmax
        double precision                      :: Am(0:nz, 0:nz)
        integer                               :: iz

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
!         allocate(this%zfilt(0:nz))

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

        ! Call parent class initialise
        call this%init_decomposition

        ! Initialise arrays for zdiffuse:
        allocate(this%D2(1:nz-1, 1:nz-1))
        allocate(this%eye(1:nz-1, 1:nz-1))

        allocate(this%D2NF(0:nz, 0:nz))
        allocate(this%eyeNF(0:nz, 0:nz))

        this%eye = zero
        Am = zero
        do iz = 1, nz-1
            this%eye(iz, iz) = one
        enddo
        Am = matmul(this%d1z, this%d1z)

        this%D2   = Am(1:nz-1,1:nz-1)
        this%D2NF = Am

        this%eyeNF = zero
        do iz = 0, nz
            this%eyeNF(iz, iz) = one
        enddo

        if (present(filter_method)) then
            call this%init_filter(filter_method)
        else
            call this%init_filter("none")
        endif

    end subroutine

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine finalise(this)
        class (cheby_layout_t), intent(inout) :: this

        if (this%l_initialised) then
            deallocate(this%d1z)
            deallocate(this%d2z)
            deallocate(this%zcheb)
!             deallocate(this%zfilt)
            deallocate(this%D2)
            deallocate(this%eye)
            this%l_initialised = .false.
        endif

        call finalise_cheby

        call this%finalise_filter

    end subroutine finalise

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine finalise_filter(this)
        class (cheby_layout_t), intent(inout) :: this

        if (allocated(this%zfilt)) then
            deallocate(this%zfilt)
        endif

        if (allocated(this%filt)) then
            deallocate(this%filt)
        endif

    end subroutine finalise_filter

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_z_axis(this)
        class (cheby_layout_t), intent(in) :: this
        double precision                   :: get_z_axis(0:nz)

        get_z_axis = center(3) - hl(3) * this%zcheb

    end function get_z_axis

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! fc  - complete field (physical space)
    ! sf  - full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
    ! cfc - copy of complete field (physical space)
    subroutine decompose_physical(this, fc, sf)
        class (cheby_layout_t), intent(in)  :: this
        double precision,       intent(in)  :: fc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision,       intent(out) :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))

        call fftxyp2s(fc, sf)

        call this%decompose_semi_spectral(sf)

    end subroutine decompose_physical

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! in : complete field (semi-spectral space)
    ! out: full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
    subroutine decompose_semi_spectral(this, sfc)
        class (cheby_layout_t), intent(in)    :: this
        double precision,       intent(inout) :: sfc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        integer                               :: iz

        ! subtract harmonic part
        !$omp parallel do
        do iz = 1, nz-1
            sfc(iz, :, :) = sfc(iz, :, :) - (sfc(0,  :, :) * this%phim(iz, :, :) + &
                                             sfc(nz, :, :) * this%phip(iz, :, :))
        enddo
        !$omp end parallel do

    end subroutine decompose_semi_spectral

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! sf  - full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
    ! fc  - complete field (physical space)
    ! sfc - complete field (semi-spectral space)
    subroutine combine_physical(this, sf, fc)
        class (cheby_layout_t), intent(in)  :: this
        double precision,       intent(in)  :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision,       intent(out) :: fc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision                    :: sfc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))

        sfc = sf

        call this%combine_semi_spectral(sfc)

        ! transform to physical space as fc:
        call fftxys2p(sfc, fc)

    end subroutine combine_physical

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! in : full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
    ! out: complete field (semi-spectral space)
    subroutine combine_semi_spectral(this, sf)
        class (cheby_layout_t), intent(in)    :: this
        double precision,       intent(inout) :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        integer                               :: iz

        ! add harmonic part to sfc:
        !$omp parallel do
        do iz = 1, nz-1
            sf(iz, :, :) = sf(iz, :, :) + sf(0,  :, :) * this%phim(iz, :, :) &
                                        + sf(nz, :, :) * this%phip(iz, :, :)
        enddo
        !$omp end parallel do

    end subroutine combine_semi_spectral

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_local_sum(this, ff) result(res)
        class (cheby_layout_t), intent(in) :: this
        double precision,       intent(in) :: ff(box%lo(3):box%hi(3), &
                                                 box%lo(2):box%hi(2), &
                                                 box%lo(1):box%hi(1))
        double precision                   :: res
        integer                            :: iz

        res = zero
        do iz = box%lo(3), box%hi(3)
            res = res + this%zccw(iz) * sum(ff(iz, :, :))
        enddo

        res = res * f12 * dble(nz)

    end function get_local_sum

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine diffz(this, fs, ds, l_decomposed)
        class (cheby_layout_t), intent(in)  :: this
        double precision,       intent(in)  :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision,       intent(out) :: ds(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        logical,                intent(in)  :: l_decomposed

        if (l_decomposed) then
            call this%decomposed_diffz(fs, ds)
        else
            call this%semi_spectral_diffz(fs, ds)
        endif

    end subroutine diffz

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine decomposed_diffz(this, fs, ds)
        class (cheby_layout_t), intent(in)  :: this
        double precision,       intent(in)  :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision,       intent(out) :: ds(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision                    :: as(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        integer                             :: iz

        ! Calculate the derivative in fully decomposed form:
        !$omp parallel workshare
        as(0,      :, :) = zero
        as(1:nz-1, :, :) = fs(1:nz-1, :, :)
        as(nz,     :, :) = zero
        !$omp end parallel workshare

        call this%zderiv(as, ds)

        !Calculate the derivative of the linear part in semi-spectral space
        !and add both contributions:
        !$omp parallel do private(iz)  default(shared)
        do iz = 0, nz
            ds(iz, :, :) = ds(iz, :, :) + fs(0,  :, :) * this%dphim(iz, :, :)  &
                                        + fs(nz, :, :) * this%dphip(iz, :, :)
        enddo
        !$omp end parallel do

        call this%decompose_semi_spectral(ds)

    end subroutine decomposed_diffz

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine semi_spectral_diffz(this, fs, ds)
        class (cheby_layout_t), intent(in)  :: this
        double precision,       intent(in)  :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision,       intent(out) :: ds(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision                    :: gs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))

        gs = fs
        call this%decompose_semi_spectral(gs)

        call this%decomposed_diffz(gs, ds)

        call this%combine_semi_spectral(ds)

    end subroutine semi_spectral_diffz

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! This is only calculated on the MPI rank having kx = ky = 0
    function calc_decomposed_mean(this, fs) result(savg)
        class (cheby_layout_t), intent(in) :: this
        double precision,       intent(in) :: fs(0:nz,                &
                                                 box%lo(2):box%hi(2), &
                                                 box%lo(1):box%hi(1))
        double precision                   :: savg, c
        integer                            :: iz

        if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
            savg = this%zccw(0) * fs(0, 0, 0)
            do iz = 1, nz-1
                ! add back the harmonic part
                c = fs(iz, 0, 0) + fs(0,  0, 0) * this%phim(iz, 0, 0) &
                                 + fs(nz, 0, 0) * this%phip(iz, 0, 0)

                savg = savg + this%zccw(iz) * c
            enddo
            savg = savg + this%zccw(nz) * fs(nz, 0, 0)

            ! The factor f12 * extent(3) comes from the mapping [-1, 1] to [a, b]
            ! where the Chebyshev points are given in [-1, 1]
            ! z = (b-a) / 2 * t + (a+b)/2 for [a, b] --> dz = (b-a) / 2 * dt
            ! However, we must divide by extent(3) again in order to get the vertical domain-average.
            ! Hence, we only scale by f12.
            savg = savg * f12
        endif

    end function calc_decomposed_mean

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! This is only calculated on the MPI rank having kx = ky = 0
    subroutine adjust_decomposed_mean(this, fs, avg)
        class (cheby_layout_t), intent(in)    :: this
        double precision,       intent(inout) :: fs(0:nz,                &
                                                    box%lo(2):box%hi(2), &
                                                    box%lo(1):box%hi(1))
        double precision,       intent(in)    :: avg
        double precision                      :: savg, cor
        integer                               :: iz

        savg = this%calc_decomposed_mean(fs)

        cor = avg - savg

        if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
            do iz = 1, nz-1
                ! add back the harmonic part
                fs(iz, 0, 0) = fs(iz, 0, 0) + fs(0,  0, 0) * this%phim(iz, 0, 0) &
                                            + fs(nz, 0, 0) * this%phip(iz, 0, 0)

                ! adjust the mean
                fs(iz, 0, 0) = fs(iz, 0, 0) + cor
            enddo

            ! adjust the mean at the surfaces
            fs(0 , 0, 0) = fs(0,  0, 0) + cor
            fs(nz, 0, 0) = fs(nz, 0, 0) + cor

            do iz = 1, nz-1
                ! remove the harmonic part
                fs(iz, 0, 0) = fs(iz, 0, 0) - (fs(0,  0, 0) * this%phim(iz, 0, 0) + &
                                               fs(nz, 0, 0) * this%phip(iz, 0, 0))
            enddo
        endif

    end subroutine adjust_decomposed_mean

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Solves dg/dz = f for g, either with g(0) at z = zmin when
    ! noavg = .false. or with the average of g over z equal to 0
    ! when noavg = .true.
    ! f & g are 1D arrays over z
    ! Note: Assumes input function for zeroth modes (kx = ky = 0).
    ! *** Uses dgesv from LAPACK/BLAS ***
    subroutine zinteg(this, f, g, noavg)
        class (cheby_layout_t), intent(in)  :: this
        double precision,       intent(in)  :: f(0:nz)
        double precision,       intent(out) :: g(0:nz)
        logical,                intent(in)  :: noavg
        double precision                    :: dmat(0:nz-1, 0:nz-1), h(0:nz), gavg
        integer                             :: ipiv(0:nz-1), info, iz

        g = f

        !-----------------------------------------------------
        ! Go to semi-spectral space, uses kx = ky = 0
        do iz = 1, nz-1
            g(iz) = g(iz) + g(0)  * this%phim(iz, 0, 0) &
                          + g(nz) * this%phip(iz, 0, 0)
        enddo

        !-----------------------------------------------------
        ! Integrate starting from g = 0 at z = zmin:
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

    subroutine zdiffNF(this, fs, dt, alpha_h, alpha_v)
        class (cheby_layout_t), intent(in)    :: this
        double precision,       intent(inout) :: fs(0:nz,                 &
                                                    box%lo(2):box%hi(2),  &
                                                    box%lo(1):box%hi(1))
        double precision,       intent(in)    :: dt
        double precision,       intent(in)    :: alpha_h
        double precision,       intent(in)    :: alpha_v
        double precision                      :: Lm(0:nz, 0:nz)
        double precision                      :: Rm(0:nz, 0:nz)
        double precision                      :: rhs(0:nz)
        integer                               :: ipiv(0:nz)
        integer                               :: kx, ky, info


        call this%combine_semi_spectral(fs)

        Lm = this%eyeNF -  dt *  alpha_v * this%D2NF
        Lm(0,:)  = this%d1z(0,:)
        Lm(nz,:) = this%d1z(nz,:)
        Rm = this%eyeNF + f12 * dt *  alpha_v * this%D2NF

        call dgetrf(nz+1, nz+1, Lm, nz+1, ipiv, info)

        do kx = box%lo(1), box%hi(1)
            do ky = box%lo(2), box%hi(2)

                rhs = fs(:, ky, kx)
                !rhs = matmul(Rm, b)
                rhs(0) = 0.0d0
                rhs(nz) = 0.0d0
                !! Linear Solve in z
                call dgetrs('N', nz+1, 1, Lm, nz+1, ipiv, rhs, nz+1, info)
                !! Diffuse in x-y
                !! fs(:, ky, kx) = rhs / (one + dt * alpha_h * k2l2(ky, kx))
                fs(:, ky, kx) = exp(-alpha_h * k2l2(ky, kx) * dt) * rhs
            enddo
        enddo

        call this%decompose_semi_spectral(fs)

    end subroutine zdiffNF

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine zdiffuse(this, fs, dt, alpha_h, alpha_v)
        class (cheby_layout_t), intent(in)    :: this
        double precision,       intent(inout) :: fs(0:nz,                 &
                                                    box%lo(2):box%hi(2),  &
                                                    box%lo(1):box%hi(1))
        double precision,       intent(in)    :: dt
        double precision,       intent(in)    :: alpha_h
        double precision,       intent(in)    :: alpha_v
        double precision                      :: Lm(1:nz-1, 1:nz-1)
        double precision                      :: Rm(1:nz-1, 1:nz-1)
        double precision                      :: rhs(1:nz-1)
        integer                               :: ipiv(0:nz)
        integer                               :: kx, ky, info

        Lm = this%eye - f12 * dt *  alpha_v * this%D2
        Rm = this%eye + f12 * dt *  alpha_v * this%D2

        call dgetrf(nz-1, nz-1, Lm, nz-1, ipiv, info)

        do kx = box%lo(1), box%hi(1)
            do ky = box%lo(2), box%hi(2)
                !fs(0,  ky, kx) = exp(-alpha_h * k2l2(ky, kx) * dt) * fs(0,  ky, kx)
                !fs(nz, ky, kx) = exp(-alpha_h * k2l2(ky, kx) * dt) * fs(nz, ky, kx)
                fs(:, ky, kx) = exp(-alpha_h * k2l2(ky, kx) * dt) * fs(:, ky, kx)
            enddo
        enddo

        do kx = box%lo(1), box%hi(1)
            do ky = box%lo(2), box%hi(2)
                rhs = matmul(Rm, fs(1:nz-1, ky, kx))
                !! Linear Solve in z
                !call dgesv(nz-1, 1, Lm, nz-1, ipiv, rhs, nz-1, info)
                call dgetrs('N', nz-1, 1, Lm, nz-1, ipiv, rhs, nz-1, info)
                !! Diffuse in x-y
                !! fs(1:nz-1, ky, kx) = rhs / (one + dt * alpha_h * k2l2(ky, kx))
            enddo
        enddo
    end subroutine zdiffuse



    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Solves (d^2/dz^2 - K^2)[ds] = S in semi-spectral space where
    ! K^2 = k^2 + l^2 is the squared horizontal wavenumber and
    ! where ds initially contains the source S (this is overwritten
    ! by the solution).
    ! *** Uses dgesv from LAPACK/BLAS ***
    subroutine vertvel(this, ds, es)
        class (cheby_layout_t), intent(in)    :: this
        double precision,       intent(inout) :: ds(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision,       intent(out)   :: es(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision                      :: dmat(nz-1, nz-1), sol(nz-1)
        integer                               :: ipiv(nz-1), info
        integer                               :: kx, ky, iz

        call this%combine_semi_spectral(ds)

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
        call this%diffz(ds, es, l_decomposed=.false.)

    end subroutine vertvel

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine init_filter(this, method)
        class (cheby_layout_t),  intent(inout) :: this
        character(*),            intent(in)    :: method
        character(len=64)                      :: used_method

        !----------------------------------------------------------
        !Define de-aliasing filter:
        select case (method)
            case ("Hou & Li")
                used_method = method
                call this%init_hou_and_li_filter(l_disable_vertical=.false.)
            case ("2/3-rule")
                call this%init_23rd_rule_filter(l_disable_vertical=.false.)
                used_method = method
            case ("Hou & Li (no vertical)")
                used_method = method
                call this%init_hou_and_li_filter(l_disable_vertical=.true.)
            case ("2/3-rule (no vertical)")
                call this%init_23rd_rule_filter(l_disable_vertical=.true.)
                used_method = method
            case ("none")
                call this%init_no_filter
                used_method = "no"
            case default
                call this%init_hou_and_li_filter(l_disable_vertical=.false.)
                used_method = "Hou & Li"
        end select

#ifdef ENABLE_VERBOSE
        if (verbose) then
            call mpi_print("Using " // trim(used_method) // " de-aliasing filter.")
        endif
#endif

    end subroutine init_filter

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine apply_filter(this, fs)
        class (cheby_layout_t), intent(in)    :: this
        double precision,       intent(inout) :: fs(box%lo(3):box%hi(3), &
                                                    box%lo(2):box%hi(2), &
                                                    box%lo(1):box%hi(1))
        double precision                      :: coeffs(0:nz,                &
                                                        box%lo(2):box%hi(2), &
                                                        box%lo(1):box%hi(1))
        double precision                      :: err_e(box%lo(2):box%hi(2),  &
                                                       box%lo(1):box%hi(1))
        double precision                      :: err_o(box%lo(2):box%hi(2),  &
                                                       box%lo(1):box%hi(1))
        integer                               :: iz
        double precision                      :: fstop(box%lo(2):box%hi(2), &
                                                       box%lo(1):box%hi(1))
        double precision                      :: fsbot(box%lo(2):box%hi(2), &
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
        coeffs = this%zfilt * coeffs

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
        fs(0,  :, :) = this%filt * fsbot
        fs(nz, :, :) = this%filt * fstop

    end subroutine apply_filter

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine apply_hfilter(this, fs)
        class (cheby_layout_t), intent(in)    :: this
        double precision,       intent(inout) :: fs(box%lo(3):box%hi(3), &
                                                    box%lo(2):box%hi(2), &
                                                    box%lo(1):box%hi(1))
        integer                               :: kz

        do kz = 0, nz
         fs(:,:,kz) = this%filt * fs(:,:,kz)
        enddo

    end subroutine apply_hfilter

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    !Define Hou and Li filter (2D and 3D):
    subroutine init_hou_and_li_filter(this, l_disable_vertical)
        class(cheby_layout_t), intent(inout) :: this
        logical,               intent(in)    :: l_disable_vertical
        integer                              :: kx, ky, kz
        double precision                     :: kxmaxi, kymaxi, kzmaxi
        double precision                     :: skx(box%lo(1):box%hi(1)), &
                                                sky(box%lo(2):box%hi(2)), &
                                                skz(0:nz)

        allocate(this%zfilt(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(this%filt(box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

        kxmaxi = one / maxval(rkx)
        skx = -36.d0 * (kxmaxi * rkx(box%lo(1):box%hi(1))) ** 36
        kymaxi = one/maxval(rky)
        sky = -36.d0 * (kymaxi * rky(box%lo(2):box%hi(2))) ** 36
        kzmaxi = one/maxval(rkz)
        if (l_disable_vertical) then
            skz = zero
        else
            skz = -36.d0 * (kzmaxi * rkz) ** 36
        endif

        do kx = box%lo(1), box%hi(1)
            do ky = box%lo(2), box%hi(2)
                this%filt(ky, kx) = exp(skx(kx) + sky(ky))
                do kz = 0, nz
                    this%zfilt(kz, ky, kx) = this%filt(ky, kx) * exp(skz(kz))
                enddo
            enddo
        enddo

        !Ensure filter does not change domain mean:
        if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
            this%filt(0, 0) = one
            this%zfilt(:, 0, 0) = exp(skz)
        endif

    end subroutine init_hou_and_li_filter

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    !Define de-aliasing filter (2/3 rule):
    subroutine init_23rd_rule_filter(this, l_disable_vertical)
        class(cheby_layout_t), intent(inout) :: this
        logical,               intent(in)    :: l_disable_vertical
        integer                              :: kx, ky, kz
        double precision                     :: rkxmax, rkymax, rkzmax
        double precision                     :: skx(box%lo(1):box%hi(1)), &
                                                sky(box%lo(2):box%hi(2)), &
                                                skz(0:nz)

        allocate(this%zfilt(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(this%filt(box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

        rkxmax = maxval(rkx)
        rkymax = maxval(rky)
        rkzmax = maxval(rkz)

        do kx = box%lo(1), box%hi(1)
            if (rkx(kx) <= f23 * rkxmax) then
                skx(kx) = one
            else
                skx(kx) = zero
            endif
        enddo

        do ky = box%lo(2), box%hi(2)
            if (rky(ky) <= f23 * rkymax) then
                sky(ky) = one
            else
                sky(ky) = zero
            endif
        enddo

        if (l_disable_vertical) then
            skz = one
        else
            do kz = 0, nz
                if (rkz(kz) <= f23 * rkzmax) then
                    skz(kz) = one
                else
                    skz(kz) = zero
                endif
            enddo
        endif

        ! Take product of 1d filters:
        do kx = box%lo(1), box%hi(1)
            do ky = box%lo(2), box%hi(2)
                this%filt(ky, kx) = skx(kx) * sky(ky)
                do kz = 0, nz
                    this%zfilt(kz, ky, kx) = this%filt(ky, kx) * skz(kz)
                enddo
            enddo
        enddo

        !Ensure filter does not change domain mean:
        if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
            this%filt(0, 0) = one
            this%zfilt(:, 0, 0) = skz
        endif

    end subroutine init_23rd_rule_filter

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    !Define no filter:
    subroutine init_no_filter(this)
        class(cheby_layout_t), intent(inout) :: this

        allocate(this%zfilt(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(this%filt(box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

        this%filt = one
        this%zfilt = one

    end subroutine init_no_filter

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Input:
    ! fs - a vector of length N+1 containing function values at Chebyshev nodes in [-1, 1]
    ! Output:
    ! c - a vector of length N+1 containing the coefficients of the Chebyshev polynomials
    subroutine get_cheb_poly(this, fs, c)
        class (cheby_layout_t), intent(in)  :: this
        double precision,       intent(in)  :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision,       intent(out) :: c(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        integer                             :: kx, ky

        do kx = box%lo(1), box%hi(1)
            do ky = box%lo(2), box%hi(2)
                call cheb_poly(nz, fs(:, ky, kx), c(:, ky, kx))
            enddo
        enddo

    end subroutine get_cheb_poly

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Input:
    ! c - a vector of length N+1 containing the coefficients of the Chebyshev polynomials
    ! Output:
    ! fs - a vector of length N+1 containing the values of f(x) at the points in y
    subroutine cheb_eval(this, c, fs)
        class (cheby_layout_t), intent(in)  :: this
        double precision,       intent(in)  :: c(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision,       intent(out) :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        integer                             :: kx, ky

        do kx = box%lo(1), box%hi(1)
            do ky = box%lo(2), box%hi(2)
                call cheb_fun(nz, c(:, ky, kx), fs(:, ky, kx))
            enddo
        enddo

    end subroutine cheb_eval


    !             !------------------------------------------------------------------
!             ! Dembenek filter:
!             rkmax = zfiltering%kmax * dble(nz)
!             do iz = 0, nz
!                 zfilt(iz) = dembenek_filter(iz,                 &
!                                             rkmax,              &
!                                             zfiltering%alpha,   &
!                                             zfiltering%beta)
!             enddo
!     !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     function dembenek_filter(k, rkmax, alpha, beta) result(res)
!         integer,          intent(in) :: k
!         double precision, intent(in) :: alpha, beta, rkmax
!         double precision             :: res, x
!
!         x = dble(k)/rkmax
!         res = exp(-alpha*x**beta)
!     end function dembenek_filter

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Calculates g = df/dz
    subroutine zderiv(this, f, g)
        class (cheby_layout_t), intent(in)  :: this
        double precision,       intent(in)  :: f(0:nz, 0:this%nxym1)
        double precision,       intent(out) :: g(0:nz, 0:this%nxym1)

        g = matmul(this%d1z, f)

    end subroutine zderiv

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Calculates g = d^2f/dz^2
    subroutine zzderiv(this, f, g)
        class (cheby_layout_t), intent(in)  :: this
        double precision,       intent(in)  :: f(0:nz, 0:this%nxym1)
        double precision,       intent(out) :: g(0:nz, 0:this%nxym1)

        g = matmul(this%d2z, f)

    end subroutine zzderiv

end module cheby_layout
