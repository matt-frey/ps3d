module mss_ops
    use field_ops
    use mss_layout, only : mss_layout_t
    use constants, only : zero, f12
    use parameters, only : nz, ncell, dxi, fnzi
    use sta3dfft, only : ztrig, zfactors, rkzi, rkz, green
    use mpi_utils, only : mpi_check_for_error
    use stafft, only : dst, dct
    implicit none

    type, extends (ops_t) :: mss_ops_t

        class(mss_layout_t), pointer :: layout

    contains
        procedure :: get_z_axis

        ! Field diagnostics:
        procedure :: get_local_sum

        ! Field operations:
        procedure :: diffz
        procedure :: calc_decomposed_mean
        procedure :: adjust_decomposed_mean

        ! Specific routines:
        procedure :: vertvel
        procedure :: zinteg
    end type

contains

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_z_axis(this)
        class (mss_ops_t), intent(in)  :: this
        double precision :: get_z_axis(0:nz)
        integer          :: i

        do i = 0, nz
            get_z_axis(i) = lower(3) + dble(i) * dx(3)
        enddo

    end function get_z_axis

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_local_sum(this, ff) result(res)
        class (mss_ops_t), intent(in) :: this
        double precision,  intent(in) :: ff(box%lo(3):box%hi(3), &
                                            box%lo(2):box%hi(2), &
                                            box%lo(1):box%hi(1))
        double precision              :: res

        res = f12 * sum(ff(0,      box%lo(2):box%hi(2), box%lo(1):box%hi(1))  &
                      + ff(nz,     box%lo(2):box%hi(2), box%lo(1):box%hi(1))) &
                  + sum(ff(1:nz-1, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

    end function get_local_sum

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    !Calculates df/dz for a field f using 2nd-order differencing.
    !Here fs = f, ds = df/dz. In physical or semi-spectral space.
    subroutine diffz(this, fs, ds)
        class (mss_ops_t), intent(in)  :: this
        double precision,  intent(in)  :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision,  intent(out) :: ds(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        integer                        :: iz
        double precision               :: hdzi

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
        class (mss_ops_t), intent(in) :: this
        double precision,  intent(in) :: fs(box%lo(3):box%hi(3), &
                                            box%lo(2):box%hi(2), &
                                            box%lo(1):box%hi(1))
        double precision              :: wk(1:nz)
        double precision              :: savg

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
        class (mss_ops_t), intent(in)    :: this
        double precision,  intent(inout) :: fs(box%lo(3):box%hi(3), &
                                               box%lo(2):box%hi(2), &
                                               box%lo(1):box%hi(1))
        double precision,  intent(in)    :: avg
        double precision                 :: savg

        savg = this%calc_decomposed_mean(fs)

        if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
            ! Ensure zero global mean horizontal vorticity conservation:
            ! Remove from boundary values (0 & nz):
            fs(0 , 0, 0) = fs(0 , 0, 0) + avg - savg
            fs(nz, 0, 0) = fs(nz, 0, 0) + avg - savg
        endif

    end subroutine adjust_decomposed_mean

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine zinteg(this, f, g, noavg)
        class (mss_ops_t), intent(in)  :: this
        double precision,  intent(in)  :: f(0:nz)
        double precision,  intent(out) :: g(0:nz)
        logical,           intent(in)  :: noavg

        !First integrate the sine series in f(1:nz-1):
        g(0) = zero
        g(1:nz-1) = -rkzi * f(1:nz-1)
        g(nz) = zero

        !Transform to semi-spectral space as a cosine series:
        call dct(1, nz, g, ztrig, zfactors)

        !Add contribution from the linear function connecting the boundary values:
        g = g + f(nz) * this%layout%gamtop - f(0) * this%layout%gambot

    end subroutine zinteg

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine vertvel(this, ds, es)
        class (mss_ops_t), intent(in)    :: this
        double precision,  intent(inout) :: ds(box%lo(3):box%hi(3), &
                                               box%lo(2):box%hi(2), &
                                               box%lo(1):box%hi(1))
        double precision,  intent(out)   :: es(0:nz,                  &
                                               box%lo(2):box%hi(2),   &
                                               box%lo(1):box%hi(1))  ! semi-spectral
        double precision :: as(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))  ! semi-spectral
        double precision :: bs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))  ! semi-spectral
        integer          :: iz, kx, ky, kz


        !Calculate the boundary contributions of the source to the vertical velocity (bs)
        !and its derivative (es) in semi-spectral space:
        !$omp parallel do private(iz)  default(shared)
        do iz = 1, nz-1
            bs(iz, :, :) = ds(0,  :, :) * this%layout%thetam(iz, :, :) &
                         + ds(nz, :, :) * this%layout%thetap(iz, :, :)
        enddo
        !$omp end parallel do

        !$omp parallel do private(iz)  default(shared)
        do iz = 0, nz
            es(iz, :, :) = ds(0,  :, :) * this%layout%dthetam(iz, :, :) &
                         + ds(nz, :, :) * this%layout%dthetap(iz, :, :)
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

end module mss_ops
