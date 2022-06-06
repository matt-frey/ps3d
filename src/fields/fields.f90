! =============================================================================
!     This module specifies all fields and implements specific subroutines
!     and functions.
! =============================================================================
module fields
    use parameters, only : nx, ny, nz, vcell, ncell, ncelli
    use constants, only : zero, f12, f14
    use sta2dfft, only : dst
    use inversion_utils
    implicit none

    ! x: zonal
    ! y: meridional
    ! z: vertical
    ! Due to periodicity in x and y, the grid points in x go from 0 to nx-1
    ! and from 0 to ny-1 in y
    double precision, allocatable, dimension(:, :, :, :) :: &
        svor,   &   ! full-spectral vorticity for 1:nz-1, semi-spectral for iz = 0 and iz = nz
        vor,    &   ! vorticity vector field (\omegax, \omegay, \omegaz) in physical space
        vel,    &   ! velocity vector field (u, v, w)
        svel,   &   ! velocity vector field (u, v, w) (semi-spectral)
        svtend

    double precision, allocatable, dimension(:, :, :) :: &
        buoy,   &   ! buoyancy (physical)
        sbuoy,  &   ! full-spectral buoyancy for 1:nz-1, semi-spectral for iz = 0 and iz = nz
        diss        ! dissipation operator (fully spectral in iz=1, nz-1, semi-spectral at iz = 0 and iz = nz)

    contains

        ! Allocate all fields
        subroutine field_alloc
            if (allocated(vel)) then
                return
            endif

            allocate(vel(0:nz, 0:ny-1, 0:nx-1, 3))
            allocate(svel(0:nz, 0:nx-1, 0:ny-1, 3))

            allocate(vor(0:nz, 0:ny-1, 0:nx-1, 3))
            allocate(svor(0:nz, 0:nx-1, 0:ny-1, 3))

            allocate(svtend(0:nz, 0:nx-1, 0:ny-1, 3))

            allocate(buoy(0:nz, 0:ny-1, 0:nx-1))
            allocate(sbuoy(0:nz, 0:nx-1, 0:ny-1))

            allocate(diss(0:nz, 0:nx-1, 0:ny-1))

        end subroutine field_alloc

        ! Reset fields to zero
        subroutine field_default
            call field_alloc

            vel    = zero
            vor    = zero
            svor   = zero
            svtend = zero
            buoy   = zero
            sbuoy  = zero
            diss   = zero
        end subroutine field_default

        subroutine field_decompose(fc, sf)
            double precision, intent(in)  :: fc(0:nz, 0:ny-1, 0:nx-1)    ! complete field (physical space)
            double precision, intent(out) :: sf(0:nz, 0:nx-1, 0:ny-1)    ! full-spectral (1:nz-1),
                                                                         ! semi-spectral at iz = 0 and iz = nz
            double precision              :: sfc(0:nz, 0:ny-1, 0:nx-1)   ! complete field (semi-spectral space)
            double precision              :: sfl(1:nz-1, 0:nx-1, 0:ny-1) ! linear part in z (semi-spectral)
            double precision              :: cfc(0:nz, 0:ny-1, 0:nx-1)   ! copy of complete field (physical space)
            integer                       :: iz, kx, ky

            cfc = fc
            call fftxyp2s(cfc, sfc)

            ! get linear part
            do iz = 1, nz-1
                sfl(iz, :, :) = sfc(0, :, :) * phi00(nz - iz) + sfc(nz, :, :) * phi00(iz)
            enddo

            ! bottom z-boundary
            sf(0, :, :) = sfc(0, :, :)

            ! interior
            sf(1:nz-1, :, :) = sfc(1:nz-1, :, :) - sfl

            ! transform interior to fully spectral
            do ky = 0, ny-1
                do kx = 0, nx-1
                    call dst(1, nz, sf(1:nz, kx, ky), ztrig, zfactors)
                enddo
            enddo

            ! top z-boundary
            sf(nz, :, :) = sfc(nz, :, :)

        end subroutine field_decompose

        subroutine field_combine(sf, fc)
            double precision, intent(in)  :: sf(0:nz, 0:nx-1, 0:ny-1)    ! full-spectral (1:nz-1),
                                                                         ! semi-spectral at iz = 0 and iz = nz
            double precision, intent(out) :: fc(0:nz, 0:ny-1, 0:nx-1)    ! complete field (physical space)
            double precision              :: sfc(0:nz, 0:ny-1, 0:nx-1)   ! complete field (semi-spectral space)
            double precision              :: sfl(1:nz-1, 0:nx-1, 0:ny-1) ! linear part in z (semi-spectral)
            integer                       :: iz, kx, ky

            ! transform sf(1:nz-1, :, :) to semi-spectral space (sine transform) as the array sfc:
            do ky = 0, ny-1
                do kx = 0, nx-1
                    sfc(1:nz-1, kx, ky) = sf(1:nz-1, kx, ky)
                    sfc(nz    , kx, ky) = zero
                    call dst(1, nz, sfc(1:nz, kx, ky), ztrig, zfactors)
                enddo
            enddo
            sfc(0,  :, :) = sf(0,  :, :)
            sfc(nz, :, :) = sf(nz, :, :)

            ! get linear part and add to sfc:
            do iz = 1, nz-1
                sfl(iz, :, :) = sfc(0, :, :) * phi00(nz - iz) + sfc(nz, :, :) * phi00(iz)
            enddo

            sfc(1:nz-1, :, :) = sfc(1:nz-1, :, :) + sfl

            ! transform to physical space as fc:
            call fftxys2p(sfc, fc)

        end subroutine field_combine

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
