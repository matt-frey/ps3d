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
        ! Spectral filter:
        double precision, allocatable :: filt(:, :, :)

    contains

        procedure :: get_z_axis

        procedure :: finalise_filter

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

        ! Filters:
        procedure :: init_filter
        procedure :: apply_filter
        procedure :: apply_hfilter
        procedure, private :: init_hou_and_li_filter
        procedure, private :: init_23rd_rule_filter
        procedure, private :: init_no_filter

        ! Specific routines:
        procedure :: vertvel
        procedure :: zinteg
        procedure :: zdiffuse
        procedure :: zdiffNF

        procedure :: central_diffz
        procedure :: decomposed_diffz

    end type mss_layout_t

contains

    function get_z_axis(this)
        class (mss_layout_t), intent(in) :: this
        double precision                 :: get_z_axis(0:nz)
        integer                          :: i

        do i = 0, nz
            get_z_axis(i) = lower(3) + dble(i) * dx(3)
        enddo

    end function get_z_axis

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine finalise_filter(this)
        class (mss_layout_t), intent(inout) :: this

        if (allocated(this%filt)) then
            deallocate(this%filt)
        endif

    end subroutine finalise_filter

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! fc  - complete field (physical space)
    ! sf  - full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
    ! cfc - copy of complete field (physical space)
    subroutine decompose_physical(this, fc, sf)
        class (mss_layout_t), intent(in)  :: this
        double precision,     intent(in)  :: fc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision,     intent(out) :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))

        call fftxyp2s(fc, sf)

        call this%decompose_semi_spectral(sf)

    end subroutine decompose_physical

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

    ! sf  - full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
    ! fc  - complete field (physical space)
    ! sfc - complete field (semi-spectral space)
    subroutine combine_physical(this, sf, fc)
        class (mss_layout_t), intent(in)  :: this
        double precision,     intent(in)  :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision,     intent(out) :: fc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision                  :: sfc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))

        sfc = sf

        call this%combine_semi_spectral(sfc)

        ! transform to physical space as fc:
        call fftxys2p(sfc, fc)

    end subroutine combine_physical

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
    function calc_decomposed_mean(this, fs) result(savg)
        class (mss_layout_t), intent(in) :: this
        double precision,     intent(in) :: fs(0:nz,                &
                                               box%lo(2):box%hi(2), &
                                               box%lo(1):box%hi(1))
        double precision                 :: wk(1:nz)
        double precision                 :: savg

        if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
            ! Cast fs_S = fs - fs_L onto the z grid as wk for kx = ky = 0:
            wk(1:nz-1) = fs(1:nz-1, 0, 0)
            wk(nz) = zero
            call dst(1, nz, wk(1:nz), ztrig, zfactors)
            ! Compute average (first part is the part due to svor_L):
            savg = f12 * (fs(0, 0, 0) + fs(nz, 0, 0)) + fnzi * sum(wk(1:nz-1))
        endif
    end function calc_decomposed_mean

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! This is only calculated on the MPI rank having kx = ky = 0
    subroutine adjust_decomposed_mean(this, fs, avg)
        class (mss_layout_t), intent(in)    :: this
        double precision,     intent(inout) :: fs(0:nz,                &
                                                  box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1))
        double precision,     intent(in)    :: avg
        double precision                    :: savg

        savg = this%calc_decomposed_mean(fs)

        if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
            ! Ensure zero global mean horizontal vorticity conservation:
            ! Remove from boundary values (0 & nz):
            fs(0 , 0, 0) = fs(0 , 0, 0) + avg - savg
            fs(nz, 0, 0) = fs(nz, 0, 0) + avg - savg
        endif

    end subroutine adjust_decomposed_mean

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine init_filter(this, method)
        class (mss_layout_t),  intent(inout) :: this
        character(*),          intent(in)    :: method
        character(len=64)                    :: used_method

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
        class (mss_layout_t), intent(in)     :: this
        double precision,     intent(inout) :: fs(box%lo(3):box%hi(3), &
                                                  box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1))

        fs = this%filt * fs

    end subroutine apply_filter

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine apply_hfilter(this, fs)
        class (mss_layout_t), intent(in)     :: this
        double precision,     intent(inout) :: fs(box%lo(3):box%hi(3), &
                                                  box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1))

        fs = this%filt * fs

    end subroutine apply_hfilter

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    !Define Hou and Li filter (2D and 3D):
    subroutine init_hou_and_li_filter(this, l_disable_vertical)
        class(mss_layout_t), intent(inout) :: this
        logical,             intent(in)    :: l_disable_vertical
        integer                            :: kx, ky, kz
        double precision                   :: kxmaxi, kymaxi, kzmaxi
        double precision                   :: skx(box%lo(1):box%hi(1)), &
                                              sky(box%lo(2):box%hi(2)), &
                                              skz(0:nz)

        allocate(this%filt(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

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

    end subroutine init_hou_and_li_filter

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    !Define de-aliasing filter (2/3 rule):
    subroutine init_23rd_rule_filter(this, l_disable_vertical)
        class(mss_layout_t), intent(inout) :: this
        logical,             intent(in)    :: l_disable_vertical
        integer                            :: kx, ky, kz
        double precision                   :: rkxmax, rkymax, rkzmax
        double precision                   :: skx(box%lo(1):box%hi(1)), &
                                              sky(box%lo(2):box%hi(2)), &
                                              skz(0:nz)

        allocate(this%filt(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

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

    end subroutine init_23rd_rule_filter

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    !Define no filter:
    subroutine init_no_filter(this)
        class(mss_layout_t), intent(inout) :: this

        allocate(this%filt(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

        this%filt = one

    end subroutine init_no_filter

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine zinteg(this, f, g, noavg)
        class (mss_layout_t), intent(in)  :: this
        double precision,     intent(in)  :: f(0:nz)
        double precision,     intent(out) :: g(0:nz)
        logical,              intent(in)  :: noavg

        !First integrate the sine series in f(1:nz-1):
        g(0) = zero
        g(1:nz-1) = -rkzi * f(1:nz-1)
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

end module mss_layout
