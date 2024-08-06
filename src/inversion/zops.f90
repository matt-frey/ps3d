! This module contains all subroutines needed for z differentiation,
! integration, and solving Poisson's equation for the vertical velocity.
module zops
    use mpi_layout, only : box
    use constants, only : zero
    use parameters, only : nz, extent, hli, center, hl
    use cheby ! Import Chebyshev module to set up various matrices needed below
    use inversion_utils, only : k2l2, init_inversion
    implicit none


    private

    double precision, allocatable :: d1z(:, :), d2z(:, :) &
                                    , zcheb(:)            & ! Chebyshev grid points
                                    , zg(:)               & ! Chebyshev domains points
                                    , zccw(:)               ! Clenshaw-Curtis weights

    logical          :: l_initialised = .false.
    integer          :: nxym1 = 0

    public :: init_zops, finalise_zops, zderiv, zzderiv, zinteg, vertvel, zcheb, zg, zccw

    contains

        subroutine init_zops
            double precision :: fdz1, fdz2

            if (l_initialised) then
                return
            endif

            l_initialised = .true.

            allocate(d1z(0:nz, 0:nz))
            allocate(d2z(0:nz, 0:nz))
            allocate(zcheb(0:nz))
            allocate(zg(0:nz))
            allocate(zccw(0:nz))

            nxym1 = box%size(1) * box%size(2) - 1

            ! Note: hli = two / extent
            fdz1 = - hli(3)
            fdz2 = fdz1 * fdz1

            !-----------------------------------------------------------------
            ! Get Chebyshev points & 1st & 2nd order differentiation matrices:
            call init_cheby(zcheb, d1z, d2z)

            zg = center(3) - hl(3) * zcheb

            ! Scale d1z & d2z for the actual z limits:
            d1z = fdz1 * d1z
            d2z = fdz2 * d2z

            ! Get Clenshaw-Curtis weights:
            call clencurt(zccw)

            ! Ensure wave numbers etc are initialised:
            call init_inversion

        end subroutine init_zops

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine finalise_zops

            if (.not. l_initialised) then
                return
            endif

            deallocate(d1z)
            deallocate(d2z)
            deallocate(zcheb)
            deallocate(zg)

        end subroutine finalise_zops

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Calculates g = df/dz
        subroutine zderiv(f, g)
            double precision, intent(in)  :: f(0:nz, 0:nxym1)
            double precision, intent(out) :: g(0:nz, 0:nxym1)

            g = matmul(d1z, f)

        end subroutine zderiv

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Calculates g = d^2f/dz^2
        subroutine zzderiv(f, g)
            double precision, intent(in)  :: f(0:nz, 0:nxym1)
            double precision, intent(out) :: g(0:nz, 0:nxym1)

            g = matmul(d2z, f)

        end subroutine zzderiv

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Solves dg/dz = f for g, either with g(0) at z = zmin when
        ! noavg = .false. or with the average of g over z equal to 0
        ! when noavg = .true.
        ! f & g are 1D arrays over z,
        ! *** Uses dgesv from LAPACK/BLAS ***
        subroutine zinteg(f, g, noavg)
            double precision, intent(in)  :: f(0:nz)
            double precision, intent(out) :: g(0:nz)
            logical,          intent(in)  :: noavg
            double precision              :: dmat(0:nz-1, 0:nz-1), h(0:nz), gavg
            integer                       :: ipiv(0:nz-1), info

            !-----------------------------------------------------
            ! Integrate starting from g = 0 at z = zmin:
            g = f
            dmat = d1z(0:nz-1, 0:nz-1)
            call dgesv(nz, 1, dmat, nz, ipiv, g(0:nz-1), nz, info)
            g(nz) = zero

            if (.not. noavg) then
                return
            endif

            !-----------------------------------------------------
            ! Remove average of g:
            h = g
            dmat = d1z(0:nz-1, 0:nz-1)
            call dgesv(nz, 1, dmat, nz, ipiv, h(0:nz-1), nz, info)
            gavg = h(0) / extent(3)

            g = g + gavg

        end subroutine zinteg

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Solves (d^2/dz^2 - K^2)[ws] = S in semi-spectral space where
        ! K^2 = k^2 + l^2 is the squared horizontal wavenumber and
        ! where ws initially contains the source S (this is overwritten
        ! by the solution).
        ! *** Uses dgesv from LAPACK/BLAS ***
        subroutine vertvel(ws)
            double precision, intent(inout) :: ws(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision                :: dmat(nz-1, nz-1), sol(nz-1)
            integer                         :: ipiv(nz-1), info
            integer                         :: kx, ky, iz

            !-----------------------------------------------------------------
            ! Loop over horizontal wavenumbers and solve linear system:
            do kx = box%lo(1), box%hi(1)
                do ky = box%lo(2), box%hi(2)
                    ! Inner part of d^2/dz^2 matrix:
                    dmat = d2z(1:nz-1, 1:nz-1)

                    ! Remove K^2 down the diagonal:
                    do iz = 1, nz-1
                        dmat(iz, iz) = dmat(iz, iz) - k2l2(ky, kx)
                    enddo

                    ! Linear solve with LAPACK:
                    sol = ws(1:nz-1, ky, kx)
                    call dgesv(nz-1, 1, dmat, nz-1, ipiv, sol, nz-1, info)
                    ws(1:nz-1, ky, kx) = sol

                    ! Add zero boundary values:
                    ws(0,  ky, kx) = zero
                    ws(nz, ky, kx) = zero
                enddo
            enddo

        end subroutine vertvel

end module zops
