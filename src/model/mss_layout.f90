module mss_layout
    use constants, only : zero, f12, f23, one
    use field_layout
    use parameters, only : nz, ncell, dxi, fnzi
    use mpi_layout, only : box
    use sta3dfft, only : ztrig      &
                       , zfactors   &
                       , fftxyp2s   &
                       , fftxys2p   &
                       , rkzi       &
                       , rkz        &
                       , rkx        &
                       , rky        &
                       , green      &
                       , fftcosine
    use stafft, only : dst, dct
    use mpi_utils, only : mpi_check_for_error
    implicit none

    type, extends (layout_t) :: mss_layout_t

    private
        double precision, allocatable :: gamtop(:), gambot(:)

        ! See for definitions in
        ! Dritschel D, Frey M. The stability of inviscid Beltrami flow between parallel free-slip impermeable
        ! boundaries. Journal of Fluid Mechanics. 2023;954:A31. doi:10.1017/jfm.2022.1007
        double precision, allocatable :: thetam(:, :, :)    ! theta_{-}         (eq. 3.10)
        double precision, allocatable :: thetap(:, :, :)    ! theta_{+}         (eq. 3.11)
        double precision, allocatable :: dthetam(:, :, :)   ! dtheta_{-}/dz
        double precision, allocatable :: dthetap(:, :, :)   ! dtheta_{+}/dz
        double precision, allocatable :: phim(:, :, :)      ! phi_{-}           (eq. 3.4a)
        double precision, allocatable :: phip(:, :, :)      ! phi_{+}           (eq. 3.4b)
        double precision, allocatable :: dphim(:, :, :)     ! dphi_{-}/dz
        double precision, allocatable :: dphip(:, :, :)     ! dphi_{+}/dz

        ! Spectral filter:
        double precision, allocatable :: filt(:, :, :)

        logical :: l_initialised = .false.

    contains

        procedure :: initialise
        procedure :: finalise

        procedure :: get_z_axis

        ! Field decompositions:
        procedure :: decompose_semi_spectral
        procedure :: combine_semi_spectral

        ! Field diagnostics:
        procedure :: get_local_sum

        ! Field operations:
        procedure :: diffz
        procedure :: get_semi_spectral_mean
        procedure :: adjust_semi_spectral_mean

        ! Filters:
        procedure :: init_exp_filter
        procedure :: init_cutoff_filter
        procedure :: apply_filter
        procedure :: apply_hfilter

        ! Specific routines:
        procedure :: vertvel
        procedure :: zinteg
        procedure :: zdiffuse
        procedure :: zdiffNF

        procedure :: central_diffz
        procedure :: decomposed_diffz

        procedure, private :: set_hyperbolic_functions

    end type mss_layout_t

contains

    subroutine initialise(this)
        class (mss_layout_t), intent(inout) :: this
        double precision                    :: z(0:nz), zm(0:nz), zp(0:nz)
        double precision                    :: phip00(0:nz)
        integer                             :: kx, ky, iz

        !------------------------------------------------------------------
        ! Ensure FFT module is initialised:
        ! (this call does nothing if already initialised)
        call initialise_fft(extent)

        !---------------------------------------------------------------------
        !Define zm = zmax - z, zp = z - zmin
        z = this%get_z_axis()
        !$omp parallel do private(z)
        do iz = 0, nz
            zm(iz) = upper(3) - z(iz)
            zp(iz) = z(iz) - lower(3)
        enddo
        !$omp end parallel do

        !---------------------------------------------------------------------
        !Hyperbolic functions used for solutions of Laplace's equation:
        allocate(this%phim(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(this%phip(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(this%dphim(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(this%dphip(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(this%thetam(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(this%thetap(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(this%dthetam(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(this%dthetap(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

        do kx = box%lo(1), box%hi(1)
            do ky = max(1, box%lo(2)), box%hi(2)
                call this%set_hyperbolic_functions(kx, ky, zm, zp)
            enddo
        enddo

        ! ky = 0
        if (box%lo(2) == 0) then
            do kx = max(1, box%lo(1)), box%hi(1)
                call this%set_hyperbolic_functions(kx, 0, zm, zp)
            enddo
        endif

        phip00 = zero
        if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
            !$omp parallel workshare
            ! kx = ky = 0
            this%phim(:, 0, 0) = zm / extent(3)
            this%phip(:, 0, 0) = zp / extent(3)

            this%dphim(:, 0, 0) = - one / extent(3)
            this%dphip(:, 0, 0) =   one / extent(3)

            this%thetam(:, 0, 0) = zero
            this%thetap(:, 0, 0) = zero

            this%dthetam(:, 0, 0) = zero
            this%dthetap(:, 0, 0) = zero

            phip00 = this%phip(:, 0, 0)
            !$omp end parallel workshare
        endif

        !---------------------------------------------------------------------
        !Define gamtop as the integral of phip(iz, 0, 0) with zero average:
        allocate(this%gamtop(0:nz))
        allocate(this%gambot(0:nz))

        call MPI_Allreduce(MPI_IN_PLACE,            &
                            phip00(0:nz),           &
                            nz+1,                   &
                            MPI_DOUBLE_PRECISION,   &
                            MPI_SUM,                &
                            world%comm,             &
                            world%err)

        !$omp parallel workshare
        this%gamtop = f12 * extent(3) * (phip00 ** 2 - f13)
        !$omp end parallel workshare

        !$omp parallel do
        do iz = 0, nz
            this%gambot(iz) = this%gamtop(nz-iz)
        enddo
        !$omp end parallel do
        !Here gambot is the complement of gamtop.

        allocate(this%filt(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

        !Default: No filtering
        this%filt = one

    end subroutine initialise

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine finalise(this)
        class (mss_layout_t), intent(inout) :: this

        if (this%l_initialised) then
            this%l_initialised = .false.
            deallocate(this%filt)
        endif

        deallocate(this%gamtop)
        deallocate(this%gambot)
        deallocate(this%phim)
        deallocate(this%phip)
        deallocate(this%dphim)
        deallocate(this%dphip)
        deallocate(this%thetam)
        deallocate(this%thetap)
        deallocate(this%dthetam)

    end subroutine finalise

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_z_axis(this)
        class (mss_layout_t), intent(in) :: this
        double precision                 :: get_z_axis(0:nz)
        integer                          :: i

        do i = 0, nz
            get_z_axis(i) = lower(3) + dble(i) * dx(3)
        enddo

    end function get_z_axis

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! in : complete field (semi-spectral space)
    ! out: full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
    subroutine decompose_semi_spectral(this, sfc)
        class (mss_layout_t), intent(in)    :: this
        double precision,     intent(inout) :: sfc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision                    :: sfctop(box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        integer                             :: iz, kx, ky

        ! subtract harmonic part
        !$omp parallel do
        do iz = 1, nz-1
            sfc(iz, :, :) = sfc(iz, :, :) - (sfc(0,  :, :) * this%phim(iz, :, :) + &
                                             sfc(nz, :, :) * this%phip(iz, :, :))
        enddo
        !$omp end parallel do

        !$omp parallel workshare
        sfctop = sfc(nz, :, :)
        !$omp end parallel workshare

        ! transform interior to fully spectral
        !$omp parallel do collapse(2)
        do kx = box%lo(1), box%hi(1)
            do ky = box%lo(2), box%hi(2)
                call dst(1, nz, sfc(1:nz, ky, kx), ztrig, zfactors)
            enddo
        enddo
        !$omp end parallel do

        !$omp parallel workshare
        sfc(nz, :, :) = sfctop
        !$omp end parallel workshare

    end subroutine decompose_semi_spectral

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! in : full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
    ! out: complete field (semi-spectral space)
    subroutine combine_semi_spectral(this, sf)
        class (mss_layout_t), intent(in)    :: this
        double precision,     intent(inout) :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision                    :: sftop(box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        integer                             :: iz, kx, ky

        ! transform sf(1:nz-1, :, :) to semi-spectral space (sine transform) as the array sf:
        !$omp parallel workshare
        sftop = sf(nz, :, :)
        !$omp end parallel workshare

        !$omp parallel do collapse(2)
        do kx = box%lo(1), box%hi(1)
            do ky = box%lo(2), box%hi(2)
                sf(nz, ky, kx) = zero
                call dst(1, nz, sf(1:nz, ky, kx), ztrig, zfactors)
            enddo
        enddo
        !$omp end parallel do

        !$omp parallel workshare
        sf(nz, :, :) = sftop
        !$omp end parallel workshare

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
        class (mss_layout_t), intent(in) :: this
        double precision,     intent(in) :: ff(box%lo(3):box%hi(3), &
                                               box%lo(2):box%hi(2), &
                                               box%lo(1):box%hi(1))
        double precision                 :: res

        res = f12 * sum(ff(0,      box%lo(2):box%hi(2), box%lo(1):box%hi(1))  &
                      + ff(nz,     box%lo(2):box%hi(2), box%lo(1):box%hi(1))) &
                  + sum(ff(1:nz-1, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

    end function get_local_sum

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    !Calculates df/dz for a field f using 2nd-order differencing.
    !Here fs = f, ds = df/dz. In physical or semi-spectral space.
    subroutine diffz(this, fs, ds, l_decomposed)
        class (mss_layout_t), intent(in)  :: this
        double precision,     intent(in)  :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision,     intent(out) :: ds(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        logical,              intent(in)  :: l_decomposed

        if (l_decomposed) then
            call this%decomposed_diffz(fs, ds)
        else
            call this%central_diffz(fs, ds)
        endif

    end subroutine diffz

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    !Calculates df/dz for a field f using 2nd-order differencing.
    !Here fs = f, ds = df/dz. In physical or semi-spectral space.
    subroutine central_diffz(this, fs, ds)
        class (mss_layout_t), intent(in)  :: this
        double precision,     intent(in)  :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision,     intent(out) :: ds(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        integer                           :: iz
        double precision                  :: hdzi

        hdzi = f12 * dxi(3)

!         ! Linear extrapolation at the boundaries:
!         ! iz = 0:  (fs(1) - fs(0)) / dz
!         ! iz = nz: (fs(nz) - fs(nz-1)) / dz
!         !$omp parallel workshare
!         ds(0,  :, :) = dxi(3) * (fs(1,    :, :) - fs(0,    :, :))
!         ds(nz, :, :) = dxi(3) * (fs(nz,   :, :) - fs(nz-1, :, :))
!         !$omp end parallel workshare

        ! One-sided second order differentiation:
        ds(0,  :, :) = (4.0d0 * fs(1,    :, :) - fs(2,    :, :) - 3.0d0 * fs(0,    :, :)) * hdzi
        ds(nz, :, :) = (3.0d0 * fs(nz,   :, :) + fs(nz-2, :, :) - 4.0d0 * fs(nz-1, :, :)) * hdzi

        ! central differencing for interior cells
        !$omp parallel do private(iz) default(shared)
        do iz = 1, nz-1
            ds(iz, :, :) = (fs(iz+1, :, :) - fs(iz-1, :, :)) * hdzi
        enddo
        !$omp end parallel do

    end subroutine central_diffz

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    !Calculates df/dz for a field f in mixed-spectral space
    !Here fs = f, ds = df/dz. Both fields are in mixed-spectral space.
    ! fs - mixed-spectral space
    ! ds - derivative linear part
    ! as - derivative sine part
    subroutine decomposed_diffz(this, fs, ds)
        class(mss_layout_t), intent(in)  :: this
        double precision,    intent(in)  :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision,    intent(out) :: ds(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision                 :: as(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        integer                          :: kz, iz

        !Calculate the derivative of the linear part (ds) in semi-spectral space:
        !$omp parallel do private(iz)  default(shared)
        do iz = 0, nz
            ds(iz, :, :) = fs(0,  :, :) * this%dphim(iz, :, :)  &
                         + fs(nz, :, :) * this%dphip(iz, :, :)
        enddo
        !$omp end parallel do

        ! Calculate d/dz of this sine series:
        !$omp parallel workshare
        as(0, :, :) = zero
        !$omp end parallel workshare
        !$omp parallel do private(kz)  default(shared)
        do kz = 1, nz-1
            as(kz, :, :) = rkz(kz) * fs(kz, :, :)
        enddo
        !$omp end parallel do
        !$omp parallel workshare
        as(nz, :, :) = zero
        !$omp end parallel workshare

        !FFT these quantities back to semi-spectral space:
        call fftcosine(as)

        ! Combine vertical derivative given the sine (as) and linear (ds) parts:
        !omp parallel workshare
        ds = ds + as
        !omp end parallel workshare

        call this%decompose_semi_spectral(ds)

    end subroutine decomposed_diffz

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! This is only calculated on the MPI rank having kx = ky = 0
    function get_semi_spectral_mean(this, fs) result(savg)
        class (mss_layout_t), intent(in) :: this
        double precision,     intent(in) :: fs(0:nz,                &
                                               box%lo(2):box%hi(2), &
                                               box%lo(1):box%hi(1))
        double precision                 :: savg
        integer                          :: iz

       if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
            savg = (f12 * (fs(0, 0, 0) + fs(nz, 0, 0)) + sum(fs(1:nz-1, 0, 0))) / dble(nz)
       endif
    end function get_semi_spectral_mean

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! This is only calculated on the MPI rank having kx = ky = 0
    subroutine adjust_semi_spectral_mean(this, fs, avg)
        class (mss_layout_t), intent(in)    :: this
        double precision,     intent(inout) :: fs(0:nz,                &
                                                  box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1))
        double precision,     intent(in)    :: avg
        double precision                    :: savg, cor

        savg = this%get_semi_spectral_mean(fs)

        cor = avg - savg

        if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
            ! Ensure zero global mean horizontal vorticity conservation:
            fs(: , 0, 0) = fs(: , 0, 0) + cor
        endif

    end subroutine adjust_semi_spectral_mean


    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! @pre Expects a field in semi-spectral space!
    subroutine apply_filter(this, fs)
        class (mss_layout_t), intent(in)     :: this
        double precision,     intent(inout) :: fs(box%lo(3):box%hi(3), &
                                                  box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1))

        call this%decompose_semi_spectral(fs)

        fs = this%filt * fs

        call this%combine_semi_spectral(fs)

    end subroutine apply_filter

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! @pre Expects a field in semi-spectral space!
    subroutine apply_hfilter(this, fs)
        class (mss_layout_t), intent(in)     :: this
        double precision,     intent(inout) :: fs(box%lo(3):box%hi(3), &
                                                  box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1))

        call this%decompose_semi_spectral(fs)

        fs = this%filt * fs

        call this%combine_semi_spectral(fs)

    end subroutine apply_hfilter

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine init_exp_filter(this, alpha, beta)
        class(mss_layout_t), intent(inout) :: this
        double precision,    intent(in)    :: alpha, beta
        integer                            :: kx, ky, kz
        double precision                   :: kxmaxi, kymaxi, kzmaxi
        double precision                   :: skx(box%lo(1):box%hi(1)), &
                                              sky(box%lo(2):box%hi(2)), &
                                              skz(0:nz)

        kxmaxi = one / maxval(rkx)
        skx = - alpha * (kxmaxi * rkx(box%lo(1):box%hi(1))) ** beta
        kymaxi = one/maxval(rky)
        sky = - alpha * (kymaxi * rky(box%lo(2):box%hi(2))) ** beta
        kzmaxi = one/maxval(rkz)
        skz = - alpha * (kzmaxi * rkz) ** beta

        do kx = box%lo(1), box%hi(1)
            do ky = box%lo(2), box%hi(2)
                this%filt(0,  ky, kx) = exp(skx(kx) + sky(ky))
                this%filt(nz, ky, kx) = this%filt(0, ky, kx)
                do kz = 1, nz-1
                    this%filt(kz, ky, kx) = this%filt(0, ky, kx) * exp(skz(kz))
                enddo
            enddo
        enddo

        !Ensure filter does not change domain mean:
        if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
            this%filt(:, 0, 0) = one
        endif

    end subroutine init_exp_filter

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine init_cutoff_filter(this, cutoff)
        class(mss_layout_t), intent(inout) :: this
        double precision,    intent(in)    :: cutoff
        integer                            :: kx, ky, kz
        double precision                   :: rkxmax, rkymax, rkzmax
        double precision                   :: skx(box%lo(1):box%hi(1)), &
                                              sky(box%lo(2):box%hi(2)), &
                                              skz(0:nz)

        rkxmax = maxval(rkx)
        rkymax = maxval(rky)
        rkzmax = maxval(rkz)

        do kx = box%lo(1), box%hi(1)
            if (rkx(kx) <= cutoff * rkxmax) then
                skx(kx) = one
            else
                skx(kx) = zero
            endif
        enddo

        do ky = box%lo(2), box%hi(2)
            if (rky(ky) <= cutoff * rkymax) then
                sky(ky) = one
            else
                sky(ky) = zero
            endif
        enddo

        do kz = 0, nz
            if (rkz(kz) <= cutoff * rkzmax) then
                skz(kz) = one
            else
                skz(kz) = zero
            endif
        enddo

        ! Take product of 1d filters:
        do kx = box%lo(1), box%hi(1)
            do ky = box%lo(2), box%hi(2)
                this%filt(0,  ky, kx) = skx(kx) * sky(ky)
                this%filt(nz, ky, kx) = this%filt(0, ky, kx)
                do kz = 1, nz-1
                    this%filt(kz, ky, kx) = this%filt(0, ky, kx) * skz(kz)
                enddo
            enddo
        enddo

        !Ensure filter does not change domain mean:
        if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
            this%filt(:, 0, 0) = one
        endif

    end subroutine init_cutoff_filter

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine zinteg(this, f, g, noavg)
        class (mss_layout_t), intent(in)  :: this
        double precision,     intent(in)  :: f(0:nz)
        double precision,     intent(out) :: g(0:nz)
        logical,              intent(in)  :: noavg
        integer                           :: iz

        !--------------------------------------------------
        ! Decompose to mixed-spectral:

        ! subtract harmonic part
        !$omp parallel do
        do iz = 1, nz-1
            g(iz) = f(iz) - (f(0)  * this%phim(iz, 0, 0) + &
                             f(nz) * this%phip(iz, 0, 0))
        enddo
        !$omp end parallel do

        ! transform interior to fully spectral
        call dst(1, nz, g(1:nz), ztrig, zfactors)

        !--------------------------------------------------
        !First integrate the sine series in f(1:nz-1):
        g(0) = zero
        g(1:nz-1) = -rkzi * g(1:nz-1)
        g(nz) = zero

        !Transform to semi-spectral space as a cosine series:
        call dct(1, nz, g, ztrig, zfactors)

        !Add contribution from the linear function connecting the boundary values:
        g = g + f(nz) * this%gamtop - f(0) * this%gambot

    end subroutine zinteg

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine zdiffuse(this, fs, dt, alpha_h, alpha_v)
        class (mss_layout_t), intent(in)    :: this
        double precision,     intent(inout) :: fs(0:nz,                 &
                                                  box%lo(2):box%hi(2),  &
                                                  box%lo(1):box%hi(1))
        double precision,     intent(in)    :: dt
        double precision,     intent(in)    :: alpha_h
        double precision,     intent(in)    :: alpha_v
        ! Do nothing here
    end subroutine zdiffuse

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine zdiffNF(this, fs, dt, alpha_h, alpha_v)
        class (mss_layout_t), intent(in)    :: this
        double precision,     intent(inout) :: fs(0:nz,                 &
                                                  box%lo(2):box%hi(2),  &
                                                  box%lo(1):box%hi(1))
        double precision,     intent(in)    :: dt
        double precision,     intent(in)    :: alpha_h
        double precision,     intent(in)    :: alpha_v
        ! Do nothing here
    end subroutine zdiffNF

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine vertvel(this, ds, es)
        class (mss_layout_t), intent(in)    :: this
        double precision,     intent(inout) :: ds(box%lo(3):box%hi(3),  &
                                                  box%lo(2):box%hi(2),  &
                                                  box%lo(1):box%hi(1))
        double precision,     intent(out)   :: es(0:nz,                 &
                                                  box%lo(2):box%hi(2),  &
                                                  box%lo(1):box%hi(1))  ! semi-spectral
        double precision                    :: as(0:nz,                 &
                                                  box%lo(2):box%hi(2),  &
                                                  box%lo(1):box%hi(1))  ! semi-spectral
        double precision                    :: bs(0:nz,                 &
                                                  box%lo(2):box%hi(2),  &
                                                  box%lo(1):box%hi(1))  ! semi-spectral
        integer                             :: iz, kx, ky, kz


        call this%decompose_semi_spectral(ds)

        !Calculate the boundary contributions of the source to the vertical velocity (bs)
        !and its derivative (es) in semi-spectral space:
        !$omp parallel do private(iz)  default(shared)
        do iz = 1, nz-1
            bs(iz, :, :) = ds(0,  :, :) * this%thetam(iz, :, :) &
                         + ds(nz, :, :) * this%thetap(iz, :, :)
        enddo
        !$omp end parallel do

        !$omp parallel do private(iz)  default(shared)
        do iz = 0, nz
            es(iz, :, :) = ds(0,  :, :) * this%dthetam(iz, :, :) &
                         + ds(nz, :, :) * this%dthetap(iz, :, :)
        enddo
        !$omp end parallel do

        !Invert Laplacian to find the part of w expressible as a sine series:
        !$omp parallel workshare
        ds(1:nz-1, :, :) = green(1:nz-1, :, :) * ds(1:nz-1, :, :)
        !$omp end parallel workshare

        ! Calculate d/dz of this sine series:
        !$omp parallel workshare
        as(0, :, :) = zero
        !$omp end parallel workshare
        !$omp parallel do private(iz)  default(shared)
        do kz = 1, nz-1
            as(kz, :, :) = rkz(kz) * ds(kz, :, :)
        enddo
        !$omp end parallel do
        !$omp parallel workshare
        as(nz, :, :) = zero
        !$omp end parallel workshare

        !FFT these quantities back to semi-spectral space:
        !$omp parallel do collapse(2) private(kx, ky)
        do kx = box%lo(1), box%hi(1)
            do ky = box%lo(2), box%hi(2)
                call dct(1, nz, as(0:nz, ky, kx), ztrig, zfactors)
                call dst(1, nz, ds(1:nz, ky, kx), ztrig, zfactors)
            enddo
        enddo
        !$omp end parallel do

        !Combine vertical velocity (ds) and its derivative (es) given the sine and linear parts:
        !$omp parallel workshare
        ds(0     , :, :) = zero
        ds(1:nz-1, :, :) = ds(1:nz-1, :, :) + bs(1:nz-1, :, :)
        ds(nz    , :, :) = zero
        es = es + as
        !$omp end parallel workshare

    end subroutine vertvel

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! for kx > 0 and ky >= 0 or kx >= 0 and ky > 0
    subroutine set_hyperbolic_functions(this, kx, ky, zm, zp)
        class(mss_layout_t),  intent(inout) :: this
        integer,              intent(in)    :: kx, ky
        double precision,     intent(in)    :: zm(0:nz), zp(0:nz)
        double precision                    :: R(0:nz), Q(0:nz), k2ifac
        double precision                    :: ef, em(0:nz), ep(0:nz), Lm(0:nz), Lp(0:nz)
        double precision                    :: fac, div, kl

        kl = sqrt(k2l2(ky, kx))
        fac = kl * extent(3)
        ef = exp(- fac)
#ifndef NDEBUG
        ! To avoid "Floating-point exception - erroneous arithmetic operation"
        ! when ef is really small.
        ef = max(ef, sqrt(tiny(ef)))
#endif
        div = one / (one - ef**2)
        k2ifac = f12 * k2l2i(ky, kx)

        Lm = kl * zm
        Lp = kl * zp

        ep = exp(- Lp)
        em = exp(- Lm)

#ifndef NDEBUG
        ! To avoid "Floating-point exception - erroneous arithmetic operation"
        ! when ep and em are really small.
        ep = max(ep, sqrt(tiny(ep)))
        em = max(em, sqrt(tiny(em)))
#endif

        this%phim(:, ky, kx) = div * (ep - ef * em)
        this%phip(:, ky, kx) = div * (em - ef * ep)

        this%dphim(:, ky, kx) = - kl * div * (ep + ef * em)
        this%dphip(:, ky, kx) =   kl * div * (em + ef * ep)

        Q = div * (one + ef**2)
        R = div * two * ef

        this%thetam(:, ky, kx) = k2ifac * (R * Lm * this%phip(:, ky, kx) - &
                                           Q * Lp * this%phim(:, ky, kx))
        this%thetap(:, ky, kx) = k2ifac * (R * Lp * this%phim(:, ky, kx) - &
                                           Q * Lm * this%phip(:, ky, kx))

        this%dthetam(:, ky, kx) = - k2ifac * ((Q * Lp - one) * this%dphim(:, ky, kx) - &
                                                      R * Lm * this%dphip(:, ky, kx))
        this%dthetap(:, ky, kx) = - k2ifac * ((Q * Lm - one) * this%dphip(:, ky, kx) - &
                                                      R * Lp * this%dphim(:, ky, kx))
    end subroutine set_hyperbolic_functions

end module mss_layout
