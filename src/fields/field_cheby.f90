module field_cheby
    use field_layout
    use parameters, only : nz
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
            procedure :: get_sum
            procedure :: get_local_mean
            procedure :: get_mean

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

            res = 0.0d0

        end function get_local_sum

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        function get_sum(this, ff, l_allreduce) result(res)
            class (field_cheby_t), intent(in) :: this
            double precision,      intent(in) :: ff(box%lo(3):box%hi(3), &
                                                    box%lo(2):box%hi(2), &
                                                    box%lo(1):box%hi(1))
            logical,               intent(in) :: l_allreduce
            double precision                  :: res

            res = this%get_local_sum(ff)

            if (l_allreduce) then
                call MPI_Allreduce(MPI_IN_PLACE,            &
                                   res,                     &
                                   1,                       &
                                   MPI_DOUBLE_PRECISION,    &
                                   MPI_SUM,                 &
                                   world%comm,              &
                                   world%err)

                call mpi_check_for_error(world, &
                    "in MPI_Allreduce of field_diagnostics::get_sum.")
            else
                call mpi_blocking_reduce(res, MPI_SUM, world)
            endif

        end function get_sum

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        function get_local_mean(this, ff) result(res)
            class (field_cheby_t), intent(in) :: this
            double precision,      intent(in) :: ff(box%lo(3):box%hi(3), &
                                                    box%lo(2):box%hi(2), &
                                                    box%lo(1):box%hi(1))
            double precision                  :: res

            res = 0.0d0


        end function get_local_mean

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        function get_mean(this, ff, l_allreduce) result(mean)
            class (field_cheby_t), intent(in) :: this
            double precision,      intent(in) :: ff(box%lo(3):box%hi(3), &
                                                    box%lo(2):box%hi(2), &
                                                    box%lo(1):box%hi(1))
            logical,               intent(in) :: l_allreduce
            double precision                  :: mean

            mean = 0.0d0

        end function get_mean

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

            savg = 0.0d0

        end function calc_decomposed_mean

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! This is only calculated on the MPI rank having kx = ky = 0
        subroutine adjust_decomposed_mean(this, fs, avg)
            class (field_cheby_t), intent(in)    :: this
            double precision,      intent(inout) :: fs(box%lo(3):box%hi(3), &
                                                       box%lo(2):box%hi(2), &
                                                       box%lo(1):box%hi(1))
            double precision,      intent(in)    :: avg

        end subroutine adjust_decomposed_mean

end module field_cheby
