module cheby_filter
    use cheby
    use field_filter, only : filter_t
    use mpi_layout, only : box
    use parameters, only : nz
    use constants, only : zero, f23, one
    use sta3dfft, only : rkx, rky, rkz
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

        fs = filt * fs

    end subroutine apply2d

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    !Define Hou and Li filter (2D and 3D):
    subroutine init_hou_and_li(this)
        class(cheby_filter_t), intent(inout) :: this
        integer                              :: kx, ky, kz
        double precision                     :: kxmaxi, kymaxi, kzmaxi
        double precision                     :: skx(box%lo(1):box%hi(1)), &
                                                sky(box%lo(2):box%hi(2)), &
                                                skz(0:nz)

        allocate(this%zfilt(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(this%filt(box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

        kxmaxi = one / maxval(rkx)
        skx = -36.d0 * (kxmaxi * rkx(box%lo(1):box%hi(1))) ** 36
        kymaxi = one/maxval(rky)
        sky = -36.d0 * (kymaxi * rky(box%lo(2):box%hi(2))) ** 36
        kzmaxi = one/maxval(rkz)
        skz = -36.d0 * (kzmaxi * rkz) ** 36

        do kx = box%lo(1), box%hi(1)
            do ky = box%lo(2), box%hi(2)
                this%filt(ky, kx) = exp(skx(kx) + sky(ky))
                do kz = 0, nz
                    this%zfilt(kz, ky, kx) = this%filt(ky, kx) * exp(skz(kz))
                enddo
            enddo
        enddo

        !Ensure filter does not change domain mean:
        if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
            this%filt(0, 0) = one
            this%zfilt(:, 0, 0) = exp(skz)
        endif

    end subroutine init_hou_and_li

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    !Define de-aliasing filter (2/3 rule):
    subroutine init_23rd_rule(this)
        class(cheby_filter_t), intent(inout) :: this
        integer                              :: kx, ky, kz
        double precision                     :: rkxmax, rkymax, rkzmax
        double precision                     :: skx(box%lo(1):box%hi(1)), &
                                                sky(box%lo(2):box%hi(2)), &
                                                skz(0:nz)

        allocate(this%zfilt(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(this%filt(box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

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

        do kz = 0, nz
            if (rkz(kz) <= f23 * rkzmax) then
                skz(kz) = one
            else
                skz(kz) = zero
            endif
        enddo

        ! Take product of 1d filters:
        do kx = box%lo(1), box%hi(1)
            do ky = box%lo(2), box%hi(2)
                this%filt(ky, kx) = skx(kx) * sky(ky)
                do kz = 0, nz
                    this%zfilt(kz, ky, kx) = this%filt(ky, kx) * skz(kz)
                enddo
            enddo
        enddo

        !Ensure filter does not change domain mean:
        if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
            this%filt(0, 0) = one
            this%zfilt(:, 0, 0) = skz
        endif

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


    !             !------------------------------------------------------------------
!             ! Dembenek filter:
!             rkmax = zfiltering%kmax * dble(nz)
!             do iz = 0, nz
!                 zfilt(iz) = dembenek_filter(iz,                 &
!                                             rkmax,              &
!                                             zfiltering%alpha,   &
!                                             zfiltering%beta)
!             enddo
!     !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     function dembenek_filter(k, rkmax, alpha, beta) result(res)
!         integer,          intent(in) :: k
!         double precision, intent(in) :: alpha, beta, rkmax
!         double precision             :: res, x
!
!         x = dble(k)/rkmax
!         res = dexp(-alpha*x**beta)
!     end function dembenek_filter

end module
