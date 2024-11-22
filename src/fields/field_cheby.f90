module field_cheby
    use constants, only : zero
    use field_layout
    use parameters, only : nz, extent, dx
    use zops, only : zccw
    use mpi_layout, only : box
    use inversion_utils, only : phim, phip
    use sta3dfft, only : fftxyp2s, fftxys2p
    use mpi_collectives, only : mpi_blocking_reduce
    implicit none

    type, extends (flayout_t) :: field_cheby_t
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

    end type field_cheby_t

    contains

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
                res = res + zccw(iz) * sum(ff(iz, :, :))
            enddo

            res = res * f12 * dble(nz)

        end function get_local_sum

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine diffz(this, fs, ds)
            class (field_cheby_t), intent(in)  :: this
            double precision,      intent(in)  :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision,      intent(out) :: ds(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))


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
                    savg = savg + zccw(iz) * fs(iz, 0, 0)
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

end module field_cheby
