module mss_filter
    use field_filter, only : field_filter_t
    implicit none

    type, extends(field_filter_t) :: mss_filter_t


    contains

        procedure :: apply
        procedure, private :: init_hou_and_li
        procedure, private :: init_23rd_rule

    end type

contains

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine apply(this, fs)
        class (mss_filter_t), intent(in) :: this
        double precision,  intent(inout) :: fs(box%lo(3):box%hi(3), &
                                               box%lo(2):box%hi(2), &
                                               box%lo(1):box%hi(1))

        fs = filt * fs

    end subroutine apply

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    !Define Hou and Li filter (2D and 3D):
    subroutine init_hou_and_li(this)
        class(mss_filter_t), intent(inout) :: this
        integer                            :: kx, ky, kz
        double precision                   :: kxmaxi, kymaxi, kzmaxi
        double precision                   :: skx(box%lo(1):box%hi(1)), &
                                              sky(box%lo(2):box%hi(2)), &
                                              skz(0:nz)

        call mpi_print("Using Hou & Li de-aliasing filter.")

        kxmaxi = one / maxval(rkx)
        skx = -36.d0 * (kxmaxi * rkx(box%lo(1):box%hi(1))) ** 36
        kymaxi = one/maxval(rky)
        sky = -36.d0 * (kymaxi * rky(box%lo(2):box%hi(2))) ** 36
        kzmaxi = one/maxval(rkz)
        skz = -36.d0 * (kzmaxi * rkz) ** 36

        do kx = box%lo(1), box%hi(1)
            do ky = box%lo(2), box%hi(2)
                filt(0,  ky, kx) = dexp(skx(kx) + sky(ky))
                filt(nz, ky, kx) = filt(0, ky, kx)
                do kz = 1, nz-1
                    filt(kz, ky, kx) = filt(0, ky, kx) * dexp(skz(kz))
                enddo
            enddo
        enddo

    end subroutine init_hou_and_li

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    !Define de-aliasing filter (2/3 rule):
    subroutine init_23rd_rule
        class(mss_filter_t), intent(inout) :: this
        integer                            :: kx, ky, kz
        double precision                   :: rkxmax, rkymax, rkzmax
        double precision                   :: skx(box%lo(1):box%hi(1)), &
                                              sky(box%lo(2):box%hi(2)), &
                                              skz(0:nz)

        call mpi_print("Using 2/3-rule de-aliasing filter.")

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
                filt(0,  ky, kx) = skx(kx) * sky(ky)
                filt(nz, ky, kx) = filt(0, ky, kx)
                do kz = 1, nz-1
                    filt(kz, ky, kx) = filt(0, ky, kx) * skz(kz)
                enddo
            enddo
        enddo
    end subroutine init_23rd_rule

end module
