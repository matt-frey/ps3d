! This module contains all subroutines needed for z differentiation,
! integration, and solving Poisson's equation for the vertical velocity.
module inversion_utils
    use constants
    use parameters, only : nx, ny, nz   &
                         , dx, dxi      &
                         , extent       &
                         , ncelli
    use mpi_layout
    use mpi_environment
    use sta3dfft, only : initialise_fft &
                       , finalise_fft   &
                       , rkx            &
                       , rky            &
                       , rkz            &
                       , k2l2           &
                       , k2l2i          &
                       , fftxyp2s       &
                       , fftxys2p       &
                       , xfactors       &
                       , xtrig          &
                       , yfactors       &
                       , ytrig
    use sta2dfft, only : ptospc
    use zops
    use options, only : filtering
#ifdef ENABLE_VERBOSE
    use options, only : verbose
#endif
    use mpi_utils, only : mpi_print
    implicit none

    private

    ! Ordering in physical space: z, y, x
    ! Ordering in spectral space: z, y, x

    ! Spectral filter:
    double precision, allocatable :: filt(:, :)

    logical :: is_initialised = .false.

    public :: init_inversion        &
            , finalise_inversion    &
            , filt                  &
            , call_ptospc


    ! Pass zops public members on
    public :: zderiv                &
            , zzderiv               &
            , zinteg                &
            , vertvel               &
            , zcheb                 &
            , zg                    &
            , zccw                  &
            , apply_zfilter

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! this routine is call by the genspec2d program
        subroutine call_ptospc(pp, ss)
            double precision, intent(inout) :: pp(box%lo(2):box%hi(2), box%lo(1):box%lo(2))
            double precision, intent(out)   :: ss(box%lo(2):box%hi(2), box%lo(1):box%lo(2))
            call ptospc(nx, ny, pp, ss, xfactors, yfactors, xtrig, ytrig)
        end subroutine call_ptospc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine init_inversion

            if (is_initialised) then
                return
            endif

            is_initialised = .true.

            call initialise_fft(extent)

            !----------------------------------------------------------
            !Define de-aliasing filter:

            allocate(filt(box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

            select case (filtering)
                case ("Hou & Li")
                    call init_hou_and_li_filter
                case ("2/3-rule")
                    call init_23rd_rule_filter
                case default
                    call init_hou_and_li_filter
            end select

            !Ensure filter does not change domain mean:
            if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
                filt(0, 0) = one
            endif

            !------------------------------------------------------------------
            !Initialise vertical operations:
            call init_zops

        end subroutine init_inversion

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Define Hou and Li filter (2D and 3D):
        subroutine init_hou_and_li_filter
            integer          :: kx, ky
            double precision :: kxmaxi, kymaxi
            double precision :: skx(box%lo(1):box%hi(1)), &
                                sky(box%lo(2):box%hi(2))

#ifdef ENABLE_VERBOSE
            if (verbose) then
                call mpi_print("Using Hou & Li de-aliasing filter.")
            endif
#endif

            kxmaxi = one / maxval(rkx)
            skx = -36.d0 * (kxmaxi * rkx(box%lo(1):box%hi(1))) ** 36
            kymaxi = one/maxval(rky)
            sky = -36.d0 * (kymaxi * rky(box%lo(2):box%hi(2))) ** 36

            do kx = box%lo(1), box%hi(1)
                do ky = box%lo(2), box%hi(2)
                    filt(ky, kx) = dexp(skx(kx) + sky(ky))
                enddo
            enddo

        end subroutine init_hou_and_li_filter

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Define de-aliasing filter (2/3 rule):
        subroutine init_23rd_rule_filter
            integer          :: kx, ky
            double precision :: rkxmax, rkymax
            double precision :: skx(box%lo(1):box%hi(1)), &
                                sky(box%lo(2):box%hi(2))

#ifdef ENABLE_VERBOSE
            if (verbose) then
                call mpi_print("Using 2/3-rule de-aliasing filter.")
            endif
#endif

            rkxmax = maxval(rkx)
            rkymax = maxval(rky)

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

            ! Take product of 1d filters:
            do kx = box%lo(1), box%hi(1)
                do ky = box%lo(2), box%hi(2)
                  filt(ky, kx) = skx(kx) * sky(ky)
               enddo
            enddo
        end subroutine init_23rd_rule_filter

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine finalise_inversion

            if (is_initialised) then
                deallocate(filt)
                is_initialised = .false.
            endif

            call finalise_zops
            call finalise_fft

        end subroutine finalise_inversion

end module inversion_utils
