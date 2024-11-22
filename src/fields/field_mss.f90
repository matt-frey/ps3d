module field_mss
    use constants, only : zero
    use field_layout
    use parameters, only : nz, ncell, dxi, fnzi
    use mpi_layout, only : box
    use inversion_utils, only : phim, phip, filt
    use sta3dfft, only : ztrig, zfactors, fftxyp2s, fftxys2p
    use stafft, only : dst
    use mpi_utils, only : mpi_check_for_error
    implicit none

    type, extends (flayout_t) :: field_mss_t

    contains
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

    end type field_mss_t

contains

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! fc  - complete field (physical space)
    ! sf  - full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
    ! cfc - copy of complete field (physical space)
    subroutine decompose_physical(this, fc, sf)
        class (field_mss_t), intent(in)  :: this
        double precision,    intent(in)  :: fc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision,    intent(out) :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))

        call fftxyp2s(fc, sf)

        call this%decompose_semi_spectral(sf)

    end subroutine decompose_physical

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! in : complete field (semi-spectral space)
    ! out: full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
    subroutine decompose_semi_spectral(this, sfc)
        class (field_mss_t), intent(in)    :: this
        double precision,    intent(inout) :: sfc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision                   :: sfctop(box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        integer                            :: iz, kx, ky

        ! subtract harmonic part
        !$omp parallel do
        do iz = 1, nz-1
            sfc(iz, :, :) = sfc(iz, :, :) - (sfc(0, :, :) * phim(iz, :, :) + sfc(nz, :, :) * phip(iz, :, :))
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
        class (field_mss_t), intent(in)  :: this
        double precision,    intent(in)  :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision,    intent(out) :: fc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision                 :: sfc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))

        sfc = sf

        call this%combine_semi_spectral(sfc)

        ! transform to physical space as fc:
        call fftxys2p(sfc, fc)

    end subroutine combine_physical

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! in : full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
    ! out: complete field (semi-spectral space)
    subroutine combine_semi_spectral(this, sf)
        class (field_mss_t), intent(in)    :: this
        double precision,    intent(inout) :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision                   :: sftop(box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        integer                            :: iz, kx, ky

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
            sf(iz, :, :) = sf(iz, :, :) + sf(0, :, :) * phim(iz, :, :) + sf(nz, :, :) * phip(iz, :, :)
        enddo
        !$omp end parallel do

    end subroutine combine_semi_spectral

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_local_sum(this, ff) result(res)
        class (field_mss_t), intent(in) :: this
        double precision,    intent(in) :: ff(box%lo(3):box%hi(3), &
                                                box%lo(2):box%hi(2), &
                                                box%lo(1):box%hi(1))
        double precision                :: res

        res = f12 * sum(ff(0,      box%lo(2):box%hi(2), box%lo(1):box%hi(1))  &
                        + ff(nz,     box%lo(2):box%hi(2), box%lo(1):box%hi(1))) &
                    + sum(ff(1:nz-1, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

    end function get_local_sum

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    !Calculates df/dz for a field f using 2nd-order differencing.
    !Here fs = f, ds = df/dz. In physical or semi-spectral space.
    subroutine diffz(this, fs, ds)
        class (field_mss_t), intent(in)  :: this
        double precision,    intent(in)  :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision,    intent(out) :: ds(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        integer                          :: iz
        double precision                 :: hdzi

        hdzi = f12 * dxi(3)

        ! Linear extrapolation at the boundaries:
        ! iz = 0:  (fs(1) - fs(0)) / dz
        ! iz = nz: (fs(nz) - fs(nz-1)) / dz
        !$omp parallel workshare
        ds(0,  :, :) = dxi(3) * (fs(1,    :, :) - fs(0,    :, :))
        ds(nz, :, :) = dxi(3) * (fs(nz,   :, :) - fs(nz-1, :, :))
        !$omp end parallel workshare

        ! central differencing for interior cells
        !$omp parallel do private(iz) default(shared)
        do iz = 1, nz-1
            ds(iz, :, :) = (fs(iz+1, :, :) - fs(iz-1, :, :)) * hdzi
        enddo
        !$omp end parallel do

    end subroutine diffz

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! This is only calculated on the MPI rank having kx = ky = 0
    function calc_decomposed_mean(this, fs) result(savg)
        class (field_mss_t), intent(in) :: this
        double precision,    intent(in) :: fs(box%lo(3):box%hi(3), &
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
        class (field_mss_t), intent(in)    :: this
        double precision,    intent(inout) :: fs(box%lo(3):box%hi(3), &
                                                    box%lo(2):box%hi(2), &
                                                    box%lo(1):box%hi(1))
        double precision,    intent(in)    :: avg
        double precision                   :: savg

        savg = this%calc_decomposed_mean(fs)

        if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
            ! Ensure zero global mean horizontal vorticity conservation:
            ! Remove from boundary values (0 & nz):
            fs(0 , 0, 0) = fs(0 , 0, 0) + avg - savg
            fs(nz, 0, 0) = fs(nz, 0, 0) + avg - savg
        endif

    end subroutine adjust_decomposed_mean

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine apply_filter(this, fs)
        class (field_mss_t), intent(in)    :: this
        double precision,    intent(inout) :: fs(box%lo(3):box%hi(3), &
                                                    box%lo(2):box%hi(2), &
                                                    box%lo(1):box%hi(1))

        fs = filt * fs

    end subroutine apply_filter

end module field_mss
