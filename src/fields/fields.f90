! =============================================================================
!     This module specifies all fields and implements specific subroutines
!     and functions.
! =============================================================================
module fields
    use parameters, only : nx, ny, nz, vcell, ncell, ncelli
    use constants, only : zero, f12, f14
    implicit none

    ! x: zonal
    ! y: meridional
    ! z: vertical
    ! Due to periodicity in x and y, the grid points in x go from 0 to nx-1
    ! and from 0 to ny-1 in y
    double precision, allocatable, dimension(:, :, :, :) :: &
        svori,     &   ! full-spectral vorticity in the interior of the domain
        vor,     &   ! vorticity vector field (\omegax, \omegay, \omegaz) in physical space
        vel,       &   ! velocity vector field (u, v, w)
        svel,      &   ! velocity vector field (u, v, w) (semi-spectral)
        svtend

    double precision, allocatable, dimension(:, :, :) :: &
        buoyg,     &   ! buoyancy (physical)
        sbuoy,     &   ! buoyancy (semi-spectral)
        diss3d,    &   ! dissipation operator (spectral)
        svortop,   &   ! semi-spectral vorticity at the top z-boundary (iz = nz)
        svorbot        ! semi-spectral vorticity at the bottom z-boundary (iz = 0)

    double precision, allocatable, dimension(:, :) :: &
        diss2d         ! dissipation operator (spectral)


    contains

        ! Allocate all fields
        subroutine field_alloc
            if (allocated(vel)) then
                return
            endif

            allocate(vel(0:nz, 0:ny-1, 0:nx-1, 3))
            allocate(svel(0:nz, 0:nx-1, 0:ny-1, 3))

            allocate(vor(0:nz, 0:ny-1, 0:nx-1, 3))
            allocate(svori(0:nz, 0:nx-1, 0:ny-1, 3))

            allocate(svtend(0:nz, 0:nx-1, 0:ny-1, 3))

            allocate(buoyg(0:nz, 0:ny-1, 0:nx-1))
            allocate(sbuoy(0:nz, 0:nx-1, 0:ny-1))

            allocate(diss2d(0:nx-1, 0:ny-1))
            allocate(diss3d(0:nz, 0:nx-1, 0:ny-1))

            allocate(svortop(0:nx-1, 0:ny-1, 3))
            allocate(svorbot(0:nx-1, 0:ny-1, 3))

        end subroutine field_alloc

        ! Reset fields to zero
        subroutine field_default
            call field_alloc

            vel     = zero
            vor   = zero
            svori   = zero
            svtend  = zero
            buoyg   = zero
            sbuoy   = zero
            diss2d  = zero
            diss3d  = zero
            svortop = zero
            svorbot = zero
        end subroutine field_default

        subroutine field_decompose(sfc, sfi, sftop, sfbot)
            double precision, intent(in)    :: sfc(0:nz, 0:nx-1, 0:ny-1) ! semi-spectral complete field
            double precision, intent(inout) :: sfi(0:nz, 0:nx-1, 0:ny-1) ! semi-spectral interior field
            double precision, intent(inout) :: sftop(0:nx-1, 0:ny-1)     ! semi-spectral top z-boundary layer
            double precision, intent(inout) :: sfbot(0:nx-1, 0:ny-1)     ! semi-spectral bottom z-boundary layer
            double precision                :: sfl(0:nz, 0:nx-1, 0:ny-1) ! linear part in z
            integer                         :: iz

            ! get top and bottom layer
            sftop = sfc(nz, :, :)
            sfbot = sfc(0,  :, :)

            ! get linear part
            do iz = 0, nz
                sfl(iz, :, :) = sfbot * phi00(nz - iz) + sftop * phi00(iz)
            enddo

            sfi = sfc - sfl

        end subroutine field_decompose

        function get_kinetic_energy() result(ke)
            double precision :: ke

            ke = f12 * sum(vel(1:nz-1, :, :, 1) ** 2      &
                         + vel(1:nz-1, :, :, 2) ** 2      &
                         + vel(1:nz-1, :, :, 3) ** 2)     &
               + f14 * sum(vel(0,  :, :, 1) ** 2          &
                         + vel(0,  :, :, 2) ** 2          &
                         + vel(0,  :, :, 3) ** 2)         &
               + f14 * sum(vel(nz, :, :, 1) ** 2          &
                         + vel(nz, :, :, 2) ** 2          &
                         + vel(nz, :, :, 3) ** 2)

            ! multiply with total volume
            ke = ke * vcell * dble(ncell)
        end function get_kinetic_energy

        function get_enstrophy() result(en)
            double precision :: en

            en = f12 * sum(vor(1:nz-1, :, :, 1) ** 2      &
                         + vor(1:nz-1, :, :, 2) ** 2      &
                         + vor(1:nz-1, :, :, 3) ** 2)     &
               + f14 * sum(vor(0,  :, :, 1) ** 2          &
                         + vor(0,  :, :, 2) ** 2          &
                         + vor(0,  :, :, 3) ** 2)         &
               + f14 * sum(vor(nz, :, :, 1) ** 2          &
                         + vor(nz, :, :, 2) ** 2          &
                         + vor(nz, :, :, 3) ** 2)

            ! multiply with total volume
            en = en * vcell * dble(ncell)

        end function get_enstrophy

        function get_mean_vorticity() result(vormean)
            double precision :: vormean(3)
            integer          :: nc

            do nc = 1, 3
                vormean(nc) =       sum(vor(1:nz-1, :, :, nc)) &
                            + f12 * sum(vor(0,      :, :, nc)) &
                            + f12 * sum(vor(nz,     :, :, nc))
            enddo

            vormean = vormean * ncelli
        end function get_mean_vorticity

end module fields
