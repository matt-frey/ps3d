module field_cheby
    use constants, only : zero
    use field_layout
    use parameters, only : nz, extent, dx
    use zops, only : zccw
    use mpi_layout, only : box
    use inversion_utils, only : phim, phip, filt
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
            procedure :: apply_filter

            ! Private procedures:
            procedure, private :: get_cheb_poly
            procedure, private :: cheb_eval

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

end module field_cheby
