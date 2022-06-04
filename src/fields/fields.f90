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
        svortg,    &   ! vorticity vector field in semi-spectral space
        vortg,     &   ! vorticity vector field (\omegax, \omegay, \omegaz) in physical space
        velog,     &   ! velocity vector field (u, v, w)
        svelog,    &   ! velocity vector field (u, v, w) (semi-spectral)
        svtend

    double precision, allocatable, dimension(:, :, :) :: &
        buoyg,     &   ! buoyancy (physical)
        sbuoyg         ! buoyancy (semi-spectral)

    double precision, allocatable, dimension(:, :) :: &
        diss           ! dissipation operator (spectral)


    contains

        ! Allocate all fields
        subroutine field_alloc
            if (allocated(velog)) then
                return
            endif

            allocate(velog(0:nz, 0:ny-1, 0:nx-1, 3))
            allocate(svelog(0:nz, 0:nx-1, 0:ny-1, 3))

            allocate(vortg(0:nz, 0:ny-1, 0:nx-1, 3))
            allocate(svortg(0:nz, 0:nx-1, 0:ny-1, 3))

            allocate(svtend(0:nz, 0:nx-1, 0:ny-1, 3))

            allocate(buoyg(0:nz, 0:ny-1, 0:nx-1))
            allocate(sbuoyg(0:nz, 0:nx-1, 0:ny-1))

            allocate(diss(0:nx-1, 0:ny-1))

        end subroutine field_alloc

        ! Reset fields to zero
        subroutine field_default
            call field_alloc

            velog    = zero
            vortg    = zero
            svortg   = zero
            svtend   = zero

            buoyg    = zero
            sbuoyg   = zero
            diss     = zero
        end subroutine field_default


        function get_kinetic_energy() result(ke)
            double precision :: ke

            ke = f12 * sum(velog(1:nz-1, :, :, 1) ** 2      &
                         + velog(1:nz-1, :, :, 2) ** 2      &
                         + velog(1:nz-1, :, :, 3) ** 2)     &
               + f14 * sum(velog(0,  :, :, 1) ** 2          &
                         + velog(0,  :, :, 2) ** 2          &
                         + velog(0,  :, :, 3) ** 2)         &
               + f14 * sum(velog(nz, :, :, 1) ** 2          &
                         + velog(nz, :, :, 2) ** 2          &
                         + velog(nz, :, :, 3) ** 2)

            ! multiply with total volume
            ke = ke * vcell * dble(ncell)
        end function get_kinetic_energy

        function get_enstrophy() result(en)
            double precision :: en

            en = f12 * sum(vortg(1:nz-1, :, :, 1) ** 2      &
                         + vortg(1:nz-1, :, :, 2) ** 2      &
                         + vortg(1:nz-1, :, :, 3) ** 2)     &
               + f14 * sum(vortg(0,  :, :, 1) ** 2          &
                         + vortg(0,  :, :, 2) ** 2          &
                         + vortg(0,  :, :, 3) ** 2)         &
               + f14 * sum(vortg(nz, :, :, 1) ** 2          &
                         + vortg(nz, :, :, 2) ** 2          &
                         + vortg(nz, :, :, 3) ** 2)

            ! multiply with total volume
            en = en * vcell * dble(ncell)

        end function get_enstrophy

        function get_mean_vorticity() result(vormean)
            double precision :: vormean(3)
            integer          :: nc

            do nc = 1, 3
                vormean(nc) =       sum(vortg(1:nz-1, :, :, nc)) &
                            + f12 * sum(vortg(0,      :, :, nc)) &
                            + f12 * sum(vortg(nz,     :, :, nc))
            enddo

            vormean = vormean * ncelli
        end function get_mean_vorticity

end module fields
