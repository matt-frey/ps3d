module cheby_filter
    use cheby
    use field_filter, only : filter_t
    use mpi_layout, only : box
    use parameters, only : nz
    use constants, only : zero
    implicit none

    type, extends(filter_t) :: cheby_filter_t

        ! Filter for the Chebyshev cofficients:
        double precision, allocatable :: zfilt(:, :, :)

        ! Filter for the surfaces:
        double precision, allocatable :: filt(:, :)

    contains

        procedure :: apply
        procedure :: apply2d
        procedure, private :: init_hou_and_li
        procedure, private :: init_23rd_rule
        procedure, private :: get_cheb_poly
        procedure, private :: cheb_eval

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
        coeffs = this%zfilt * coeffs

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
        fs(0,  :, :) = this%filt * fsbot
        fs(nz, :, :) = this%filt * fstop

    end subroutine apply

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine apply2d(this, fs)
        class (cheby_filter_t), intent(in) :: this
        double precision,       intent(inout) :: fs(box%lo(2):box%hi(2), &
                                                    box%lo(1):box%hi(1))

    end subroutine apply2d

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    !Define Hou and Li filter (2D and 3D):
    subroutine init_hou_and_li(this)
        class(cheby_filter_t), intent(inout) :: this

        allocate(this%zfilt(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(this%filt(box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

    end subroutine init_hou_and_li

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    !Define de-aliasing filter (2/3 rule):
    subroutine init_23rd_rule(this)
        class(cheby_filter_t), intent(inout) :: this

        allocate(this%zfilt(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(this%filt(box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

    end subroutine init_23rd_rule

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Input:
    ! fs - a vector of length N+1 containing function values at Chebyshev nodes in [-1, 1]
    ! Output:
    ! c - a vector of length N+1 containing the coefficients of the Chebyshev polynomials
    subroutine get_cheb_poly(this, fs, c)
        class (cheby_filter_t), intent(in)  :: this
        double precision,       intent(in)  :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision,       intent(out) :: c(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        integer                             :: kx, ky

        do kx = box%lo(1), box%hi(1)
            do ky = box%lo(2), box%hi(2)
                call cheb_poly(nz, fs(:, ky, kx), c(:, ky, kx))
            enddo
        enddo

    end subroutine get_cheb_poly

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Input:
    ! c - a vector of length N+1 containing the coefficients of the Chebyshev polynomials
    ! Output:
    ! fs - a vector of length N+1 containing the values of f(x) at the points in y
    subroutine cheb_eval(this, c, fs)
        class (cheby_filter_t), intent(in)  :: this
        double precision,       intent(in)  :: c(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision,       intent(out) :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        integer                             :: kx, ky

        do kx = box%lo(1), box%hi(1)
            do ky = box%lo(2), box%hi(2)
                call cheb_fun(nz, c(:, ky, kx), fs(:, ky, kx))
            enddo
        enddo

    end subroutine cheb_eval

end module
