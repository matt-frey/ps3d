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
    implicit none

    private

    ! Ordering in physical space: z, y, x
    ! Ordering in spectral space: z, y, x

    double precision, allocatable :: green(:, :, :)

    ! Spectral filter:
    double precision, allocatable :: filt(:, :)

    logical :: is_initialised = .false.

    public :: init_inversion        &
            , finalise_inversion    &
            , filt                  &
            , green                 &
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
            integer          :: kx, ky, kz
            double precision :: kxmaxi, kymaxi
            double precision :: skx(box%lo(1):box%hi(1)), &
                                sky(box%lo(2):box%hi(2))

            if (is_initialised) then
                return
            endif

            is_initialised = .true.

            call initialise_fft(extent)

            !----------------------------------------------------------
            !Define Hou and Li filter (2D and 3D):
            allocate(filt(box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

            kxmaxi = one / maxval(rkx)
            skx = -36.d0 * (kxmaxi * rkx(box%lo(1):box%hi(1))) ** 36
            kymaxi = one/maxval(rky)
            sky = -36.d0 * (kymaxi * rky(box%lo(2):box%hi(2))) ** 36

            do kx = box%lo(1), box%hi(1)
                do ky = box%lo(2), box%hi(2)
                     filt(ky, kx) = dexp(skx(kx) + sky(ky))
               enddo
            enddo

            !Ensure filter does not change domain mean:
            if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
                filt(0, 0) = one
            endif

            !---------------------------------------------------------------------
            !Define Green function:
            allocate(green(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

            !$omp parallel do
            do kz = 1, nz
                green(kz, :, :) = - one / (k2l2 + rkz(kz) ** 2)
            enddo
            !$omp end parallel do
            green(0, :, :) = - k2l2i

            call init_zops

        end subroutine init_inversion

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine finalise_inversion
            deallocate(green)
            deallocate(filt)

            call finalise_zops
            call finalise_fft

        end subroutine finalise_inversion

end module inversion_utils
