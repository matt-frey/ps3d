module field_decomposition
    implicit none

    type, abstract :: field_decomp_t
        contains
            procedure (decomp_phys), deferred :: field_decompose_physical
            procedure (comb_phys),   deferred :: field_combine_physical
            procedure (decomp_semi), deferred :: field_decompose_semi_spectral
            procedure (comb_semi),   deferred :: field_combine_semi_spectral
    end type field_decomp_t

    interface
        subroutine decomp_phys(this, fc, sf)
            use parameters, only : nz
            use mpi_layout, only : box
            import :: field_decomp_t
            class (field_decomp_t), intent(in)  :: this
            double precision,       intent(in)  :: fc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision,       intent(out) :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        end subroutine

        subroutine comb_phys(this, sf, fc)
            use parameters, only : nz
            use mpi_layout, only : box
            import :: field_decomp_t
            class (field_decomp_t), intent(in)  :: this
            double precision,       intent(in)  :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision,       intent(out) :: fc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        end subroutine

        subroutine decomp_semi(this, sfc)
            use parameters, only : nz
            use mpi_layout, only : box
            import :: field_decomp_t
            class (field_decomp_t), intent(in)    :: this
            double precision,       intent(inout) :: sfc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        end subroutine

        subroutine comb_semi(this, sf)
            use parameters, only : nz
            use mpi_layout, only : box
            import :: field_decomp_t
            class (field_decomp_t), intent(in)    :: this
            double precision,       intent(inout) :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        end subroutine
    end interface

end module field_decomposition

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module field_mixed_spectral
    use constants, only : zero
    use field_decomposition
    use parameters, only : nz
    use mpi_layout, only : box
    use inversion_utils, only : phim, phip
    use sta3dfft, only : ztrig, zfactors, fftxyp2s, fftxys2p
    use stafft, only : dst
    implicit none
    private
    public :: field_decomp_mss_t

    type, extends (field_decomp_t) :: field_decomp_mss_t
        contains
            procedure :: field_decompose_physical
            procedure :: field_combine_physical
            procedure :: field_decompose_semi_spectral
            procedure :: field_combine_semi_spectral
    end type field_decomp_mss_t

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! fc  - complete field (physical space)
        ! sf  - full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
        ! cfc - copy of complete field (physical space)
        subroutine field_decompose_physical(this, fc, sf)
            class (field_decomp_mss_t), intent(in)  :: this
            double precision,           intent(in)  :: fc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision,           intent(out) :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))

            call fftxyp2s(fc, sf)

            call this%field_decompose_semi_spectral(sf)

        end subroutine field_decompose_physical

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! in : complete field (semi-spectral space)
        ! out: full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
        subroutine field_decompose_semi_spectral(this, sfc)
            class (field_decomp_mss_t), intent(in)    :: this
            double precision,           intent(inout) :: sfc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision                          :: sfctop(box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            integer                                   :: iz, kx, ky

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
            class (field_decomp_mss_t), intent(in)  :: this
            double precision,           intent(in)  :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision,           intent(out) :: fc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision                        :: sfc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))

            sfc = sf

            call this%field_combine_semi_spectral(sfc)

            ! transform to physical space as fc:
            call fftxys2p(sfc, fc)

        end subroutine field_combine_physical

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! in : full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
        ! out: complete field (semi-spectral space)
        subroutine field_combine_semi_spectral(this, sf)
            class (field_decomp_mss_t), intent(in)    :: this
            double precision,           intent(inout) :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision                          :: sftop(box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            integer                                   :: iz, kx, ky

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

end module field_mixed_spectral

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module field_chebyshev
    use field_decomposition
    use parameters, only : nz
    use mpi_layout, only : box
    use inversion_utils, only : phim, phip
    use sta3dfft, only : fftxyp2s, fftxys2p
    implicit none
    private
    public :: field_decomp_cheb_t

    type, extends (field_decomp_t) :: field_decomp_cheb_t
        contains
            procedure :: field_decompose_physical
            procedure :: field_combine_physical
            procedure :: field_decompose_semi_spectral
            procedure :: field_combine_semi_spectral
    end type field_decomp_cheb_t

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! fc  - complete field (physical space)
        ! sf  - full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
        ! cfc - copy of complete field (physical space)
        subroutine field_decompose_physical(this, fc, sf)
            class (field_decomp_cheb_t), intent(in)  :: this
            double precision,            intent(in)  :: fc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision,            intent(out) :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))

            call fftxyp2s(fc, sf)

            call this%field_decompose_semi_spectral(sf)

        end subroutine field_decompose_physical

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! in : complete field (semi-spectral space)
        ! out: full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
        subroutine field_decompose_semi_spectral(this, sfc)
            class (field_decomp_cheb_t), intent(in)    :: this
            double precision,            intent(inout) :: sfc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            integer                                    :: iz

            ! subtract harmonic part
            !$omp parallel do
            do iz = 1, nz-1
                sfc(iz, :, :) = sfc(iz, :, :) - (sfc(0, :, :) * phim(iz, :, :) + sfc(nz, :, :) * phip(iz, :, :))
            enddo
            !$omp end parallel do

        end subroutine field_decompose_semi_spectral

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! sf  - full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
        ! fc  - complete field (physical space)
        ! sfc - complete field (semi-spectral space)
        subroutine field_combine_physical(this, sf, fc)
            class (field_decomp_cheb_t), intent(in)  :: this
            double precision,            intent(in)  :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision,            intent(out) :: fc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision                         :: sfc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))

            sfc = sf

            call this%field_combine_semi_spectral(sfc)

            ! transform to physical space as fc:
            call fftxys2p(sfc, fc)

        end subroutine field_combine_physical

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! in : full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
        ! out: complete field (semi-spectral space)
        subroutine field_combine_semi_spectral(this, sf)
            class (field_decomp_cheb_t), intent(in)    :: this
            double precision,            intent(inout) :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            integer                                    :: iz

            ! add harmonic part to sfc:
            !$omp parallel do
            do iz = 1, nz-1
                sf(iz, :, :) = sf(iz, :, :) + sf(0, :, :) * phim(iz, :, :) + sf(nz, :, :) * phip(iz, :, :)
            enddo
            !$omp end parallel do

        end subroutine field_combine_semi_spectral

end module field_chebyshev

