module field_mss
    use constants, only : zero
    use field_layout
    use parameters, only : nz, ncell
    use mpi_layout, only : box
    use inversion_utils, only : phim, phip
    use sta3dfft, only : ztrig, zfactors, fftxyp2s, fftxys2p
    use stafft, only : dst
    use mpi_utils, only : mpi_check_for_error
    implicit none

    type, extends (flayout_t) :: field_mss_t
        contains
            ! Field decompositions:
            procedure :: field_decompose_physical
            procedure :: field_combine_physical
            procedure :: field_decompose_semi_spectral
            procedure :: field_combine_semi_spectral

            ! Field diagnostics:
            procedure :: get_mean

    end type field_mss_t

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! fc  - complete field (physical space)
        ! sf  - full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
        ! cfc - copy of complete field (physical space)
        subroutine field_decompose_physical(this, fc, sf)
            class (field_mss_t), intent(in)  :: this
            double precision,    intent(in)  :: fc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision,    intent(out) :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))

            call fftxyp2s(fc, sf)

            call this%field_decompose_semi_spectral(sf)

        end subroutine field_decompose_physical

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! in : complete field (semi-spectral space)
        ! out: full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
        subroutine field_decompose_semi_spectral(this, sfc)
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

        end subroutine field_decompose_semi_spectral

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! sf  - full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
        ! fc  - complete field (physical space)
        ! sfc - complete field (semi-spectral space)
        subroutine field_combine_physical(this, sf, fc)
            class (field_mss_t), intent(in)  :: this
            double precision,    intent(in)  :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision,    intent(out) :: fc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision                 :: sfc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))

            sfc = sf

            call this%field_combine_semi_spectral(sfc)

            ! transform to physical space as fc:
            call fftxys2p(sfc, fc)

        end subroutine field_combine_physical

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! in : full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
        ! out: complete field (semi-spectral space)
        subroutine field_combine_semi_spectral(this, sf)
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

        end subroutine field_combine_semi_spectral

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        function get_mean(this, ff, l_allreduce) result(mean)
            class (field_mss_t), intent(in) :: this
            double precision,    intent(in) :: ff(box%lo(3):box%hi(3), &
                                                  box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1))
            logical,             intent(in) :: l_allreduce
            double precision                :: mean

            ! (divide by ncell since lower and upper edge weights are halved)
            mean = (f12 * sum(ff(0,      box%lo(2):box%hi(2), box%lo(1):box%hi(1))  &
                            + ff(nz,     box%lo(2):box%hi(2), box%lo(1):box%hi(1))) &
                        + sum(ff(1:nz-1, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))) / dble(ncell)

            if (l_allreduce) then
                call MPI_Allreduce(MPI_IN_PLACE,            &
                                   mean,                    &
                                   1,                       &
                                   MPI_DOUBLE_PRECISION,    &
                                   MPI_SUM,                 &
                                   world%comm,              &
                                   world%err)

                call mpi_check_for_error(world, &
                    "in MPI_Allreduce of field_diagnostics::get_mean.")
            else
                call mpi_blocking_reduce(mean, MPI_SUM, world)
            endif

        end function get_mean

end module field_mss
