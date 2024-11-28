module cheby_filter
    use field_filter, only : field_filter_t
    implicit none

    type, extends(field_filter_t) :: cheby_filter_t


    contains

        procedure :: apply
        procedure, private :: init_hou_and_li
        procedure, private :: init_23rd_rule

    end type

contains

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine apply(this, fs)
        class (cheby_filter_t), intent(in)    :: this
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

    end subroutine apply

end module
