! =============================================================================
!     This module specifies all fields and implements specific subroutines
!     and functions.
! =============================================================================
module fields
    use parameters, only : nx, ny, nz, ncelli, dx, lower, extent, upper
    use constants, only : zero, f12, f14, one
    use merge_sort
    use inversion_utils, only : fftxys2p, diffx, diffy, central_diffz   &
                              , field_combine_semi_spectral             &
                              , field_decompose_semi_spectral           &
                              , fftxyp2s                                &
                              , field_combine_physical                  &
                              , field_decompose_physical                &
                              , integrate_decomposed_field              &
                              , surf_fftxys2p, kh, phip
    use ape_density, only : ape_den
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
        svorts      ! vorticity source in mixed spectral space

    double precision, allocatable, dimension(:, :, :) :: &
        zeta,   &   ! surface zeta in physical space
        szeta,  &   ! surface zeta in semi-spectral space
        szetas, &   ! surface zeta source in semi-spectral space
        pres        ! pressure field (physical space)

#ifdef ENABLE_BUOYANCY
    double precision, allocatable, dimension(:, :, :) :: &
        buoy,   &   ! buoyancy (physical)
        sbuoy,  &   ! full-spectral buoyancy for 1:nz-1, semi-spectral for iz = 0 and iz = nz
        sbuoys      ! buoyancy source in mixed spectral space
#endif

    double precision, allocatable, dimension(:, :) :: &
        diss        ! dissipation operator

    ! initial \xi and \eta mean
    double precision :: ini_vor_mean(2)

#ifdef ENABLE_BUOYANCY
    ! buoyancy frequency squared
    double precision :: bfsq        ! N**2
#ifdef ENABLE_PERTURBATION_MODE
    double precision, allocatable :: bbarz(:) ! N**2 * z
#endif
#endif
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
            allocate(svorts(0:nz, 0:nx-1, 0:ny-1, 2))

            allocate(zeta(0:1, 0:nx-1, 0:ny-1))
            allocate(szeta(0:1, 0:nx-1, 0:ny-1))
            allocate(szetas(0:1, 0:nx-1, 0:ny-1))

#ifdef ENABLE_BUOYANCY
            allocate(buoy(0:nz, 0:ny-1, 0:nx-1))
            allocate(sbuoy(0:nz, 0:nx-1, 0:ny-1))
            allocate(sbuoys(0:nz, 0:nx-1, 0:ny-1))

#ifdef ENABLE_PERTURBATION_MODE
            allocate(bbarz(0:nz))
#endif
#endif

            allocate(pres(0:nz, 0:ny-1, 0:nx-1))
            allocate(diss(0:nx-1, 0:ny-1))



        end subroutine field_alloc

        ! Reset fields to zero
        subroutine field_default
            call field_alloc

            svor   = zero
            vor    = zero
            vel    = zero
            svel   = zero
            svorts = zero
            zeta   = zero
            szeta  = zero
            szetas = zero
#ifdef ENABLE_BUOYANCY
            buoy   = zero
            sbuoy  = zero
            sbuoys = zero
            bfsq = zero
#endif
            pres   = zero
            diss   = zero

            ini_vor_mean = zero
        end subroutine field_default

        ! domain-averaged available potential energy
        function get_available_potential_energy() result(ape)
            double precision :: ape
#ifdef ENABLE_BUOYANCY
            integer          :: i, j, k
            double precision :: z(0:nz)

            do k = 0, nz
                z(k) = lower(3) + dble(k) * dx(3)
            enddo


            ape = zero
            do i = 0, nx-1
                do j = 0, ny-1
#ifdef ENABLE_PERTURBATION_MODE
                    buoy(:, j, i) = buoy(:, j, i) + bbarz
#endif
                    ape = ape + sum(ape_den(buoy(1:nz-1, j, i), z(1:nz-1))) &
                        + f12 *     ape_den(buoy(0,      j, i), z(0))       &
                        + f12 *     ape_den(buoy(nz,     j, i), z(nz))

#ifdef ENABLE_PERTURBATION_MODE
                    buoy(:, j, i) = buoy(:, j, i) - bbarz
#endif
                enddo
            enddo


            ape = ape * ncelli
#else
            ape = zero
#endif
        end function get_available_potential_energy

        ! domain-averaged kinetic energy
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

            ke = ke * ncelli

        end function get_kinetic_energy

        ! domain-averaged enstrophy
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

            en = en * ncelli

        end function get_enstrophy

#ifdef ENABLE_BUOYANCY
        ! domain-averaged grad(b) integral
        function get_gradb_integral() result(enb)
            double precision :: enb
            double precision :: ds(0:nz, 0:nx-1, 0:ny-1)
            double precision :: mag(0:nz, 0:ny-1, 0:nx-1)
            double precision :: dbdx(0:nz, 0:ny-1, 0:nx-1)
            double precision :: dbdy(0:nz, 0:ny-1, 0:nx-1)

            !------------------------------------
            !Obtain magnitude of buoyancy gradient
            call field_combine_semi_spectral(sbuoy)
            call diffx(sbuoy, ds)
            call fftxys2p(ds, dbdx)

            call diffy(sbuoy, ds)
            call fftxys2p(ds, dbdy)

            call central_diffz(sbuoy, mag)
            call fftxys2p(ds, mag)
            call field_decompose_semi_spectral(sbuoy)

            ! mag = |gradb|
            mag = dsqrt(dbdx ** 2 + dbdy ** 2 + mag ** 2)

            !------------------------------------
            ! Calculate domain integral of |gradb|
            enb =       sum(mag(1:nz-1, :, :)) &
                + f12 * sum(mag(0,      :, :)) &
                + f12 * sum(mag(nz,     :, :))

            enb = enb * ncelli

        end function get_gradb_integral
#endif

        function get_mean_vorticity() result(vormean)
            double precision :: vormean(3)
            integer          :: nc

            do nc = 1, 3
                !$omp parallel workshare
                vormean(nc) =       sum(vor(1:nz-1, :, :, nc)) &
                            + f12 * sum(vor(0,      :, :, nc)) &
                            + f12 * sum(vor(nz,     :, :, nc))
                !$omp end parallel workshare
            enddo

            vormean = vormean * ncelli
        end function get_mean_vorticity

        ! Obtain complete zeta in physical (vor(:, :, :, 3))
        ! and semi-spectral space (svor(:, :, :, 3)) by integrating from zeta_min (iz = 0)
        subroutine combine_zeta
            integer          :: iz, kx, ky
            double precision :: ds(0:nz, 0:nx-1, 0:ny-1)    ! mixed-spectral space
            double precision :: fp(0:nz, 0:ny-1, 0:nx-1)    ! physical space
            double precision :: as(0:nz, 0:nx-1, 0:ny-1)    ! mixed-spectral space
            double precision :: psi(0:nz, 0:nx-1, 0:ny-1)
            double precision :: psi_z(0:nz, 0:nx-1, 0:ny-1)
            double precision :: psi_x(0:nz, 0:nx-1, 0:ny-1)
            double precision :: psi_y(0:nz, 0:nx-1, 0:ny-1)
            double precision :: c, d, y, kl, z

            ! d\xi/dx in mixed spectral space
            call diffx(svor(:, :, :, 1), ds)

            ! d\eta/dy in mixed spectral space
            call diffy(svor(:, :, :, 2), as)

            ! d\xi/dx + d\eta/dy
            ds = ds + as

            call integrate_decomposed_field(ds)

            ! correction
            ! ds is semi-spectral here
            do iz = 0, nz
                z = lower(3) + dble(iz) * dx(3)
                do kx = 0, nx-1
                    do ky = 0, ny-1
                        kl = kh(kx, ky)
                        y = szeta(1, kx, ky) - ds(iz, kx, ky)
                        if (kx == 0 .and. ky == 0) then
                            kl = 1.0d0
                            y = 0.0d0
                        endif
!                         psi = (y/K) * exp{K(z - z_max)} * [1 - exp{-2K(z - z_min)}] / [1 - exp{-2KL_z}]
                        psi(iz, kx, ky) = (y/kl) * dexp(kl*(z - upper(3))) &
                                        * (1.0d0 - dexp(-2.0d0* kl*(z - lower(3)))) &
                                        / (1.0d0 - dexp(-2.0d0*kl*extent(3)))
                        psi_z(iz, kx, ky) = y * phip(iz, kx, ky)
                    enddo
                enddo
            enddo
            call diffx(psi, psi_x)
            call diffy(psi, psi_y)


            ds = ds + psi_z

            call field_decompose_semi_spectral(ds)
            call field_decompose_semi_spectral(psi_x)
            call field_decompose_semi_spectral(psi_y)

            svor(:, :, :, 1) = svor(:, :, :, 1) + psi_x
            svor(:, :, :, 2) = svor(:, :, :, 2) + psi_y

            call field_combine_physical(ds, fp)

            ! get surface zeta in physical space
            svor(0, :, :, 3) = szeta(0, :, :)
            call surf_fftxys2p(svor(0, :, :, 3), zeta(0, :, :))
            vor(0, :, :, 3) = zeta(0, :, :)

            ! get complete zeta in physical space
            do iz = 1, nz
                vor(iz, :, :, 3) = zeta(0, :, :) - fp(iz, :, :)
            enddo

            ! get complete zeta in semi-spectral space
            fp = vor(:, :, :, 3)
            call fftxyp2s(fp, svor(:, :, :, 3))

        end subroutine combine_zeta

end module fields
