! Advances fields from time t to t+dt using an iterative implicit
! trapezoidal method of the form
!
!     (F^{n+1}-F^n)/dt = (L^{n+1}+L^n)/2 + (S^{n+1}+S^n)/2
!
! for a field F, where n refers to the time level, L[F] refers to
! the linear dissipation terms (hyperdiffusion), and S[F] refers to
! the remaining source terms.

! We start with the guess S^{n+1} = S^n and iterate  niter  times
! (see parameter statement below).
module advance_mod
    use options, only : time, viscosity
    use constants
    use parameters, only : nx, ny, nz, glmin, cflpf, ncelli, fnzi
    use inversion_mod, only : vor2vel, vorticity_tendency
    use inversion_utils
    use utils, only : write_step
    use sta2dfft, only : dst
    use fields
    use jacobi, only : jacobi_eigenvalues
    implicit none

    integer :: advance_timer

    integer, parameter :: WRITE_VOR = 1234
    integer, parameter :: WRITE_ECOMP = 1235

    ! Number of iterations of above scheme:
    integer, parameter:: niter = 2

    double precision :: dt, dt2

    !Diagnostic quantities:
    double precision :: bfmax, vortmax, vortrms, ggmax, velmax, dfac
    double precision :: vortmp1, vortmp2, vortmp3, vorl1, vorl2, vorch
    integer          :: ix, iy, iz

    contains

        subroutine advance(t)
            double precision, intent(inout) :: t
            integer                         :: iter
            integer                         :: nc
            ! Spectral fields needed in time stepping:
            double precision                :: vortsm(0:nz, 0:nx-1, 0:ny-1, 3)
#ifdef ENABLE_BUOYANCY
            double precision                :: bsm(0:nz, 0:nx-1, 0:ny-1)
#endif

            !-------------------------------------------------------------------
            !Invert vorticity for velocity at current time level, say t=t^n:
            !Also, returns vorticity in physical space for use everywhere
            call vor2vel

            !Adapt the time step
            call adapt(t)

            !Write fields
            call write_step(t)

            !------------------------------------------------------------------
            !Start with a guess for F^{n+1} for all fields:

            !Calculate the source terms (sbuoys, svorts) for buoyancy (sbuoy) and
            !vorticity in spectral space:
            call source

            !Initialise iteration (dt = dt/2 below):
#ifdef ENABLE_BUOYANCY
            bsm = sbuoy + dt2 * sbuoys
            sbuoy = filt * (bsm + dt2 * sbuoys)
            call field_combine_semi_spectral(sbuoy)
            !$omp parallel private(iz)  default(shared)
            !$omp do
            do iz = 0, nz
                sbuoy(iz, :, :) = diss * sbuoy(iz, :, :)
            enddo
            !$omp end do
            !$omp end parallel
            call field_decompose_semi_spectral(sbuoy)
#endif

            ! Advance interior and boundary of vorticity
            !$omp parallel
            !$omp workshare
            vortsm = svor + dt2 * svorts
            !$omp end workshare
            !$omp end parallel

            do nc = 1, 3
                !$omp parallel
                !$omp workshare
                svor(:, :, :, nc) = filt * (vortsm(:, :, :, nc) + dt2 * svorts(:, :, :, nc))
                !$omp end workshare
                !$omp end parallel
                call field_combine_semi_spectral(svor(:, :, :, nc))
                !$omp parallel private(iz)  default(shared)
                !$omp do
                do iz = 0, nz
                    svor(iz, :, :, nc) = diss * svor(iz, :, :, nc)
                enddo
                !$omp end do
                !$omp end parallel
                call field_decompose_semi_spectral(svor(:, :, :, nc))
            enddo

            call adjust_vorticity_mean

            !diss is related to the hyperdiffusive operator (see end of adapt)

            !------------------------------------------------------------------
            !Iterate to improve estimates of F^{n+1}:
            do iter = 1, niter
                !Perform inversion at t^{n+1} from estimated quantities:
                call vor2vel

                !Calculate the source terms (sbuoys,svorts):
                call source

                !Update fields:
#ifdef ENABLE_BUOYANCY
                sbuoy = filt * (bsm + dt2 * sbuoys)
                call field_combine_semi_spectral(sbuoy)
                !$omp parallel private(iz)  default(shared)
                !$omp do
                do iz = 0, nz
                    sbuoy(iz, :, :) = diss * sbuoy(iz, :, :)
                enddo
                !$omp end do
                !$omp end parallel
                call field_decompose_semi_spectral(sbuoy)
#endif

                do nc = 1, 3
                    !$omp parallel
                    !$omp workshare
                    svor(:, :, :, nc) = filt * (vortsm(:, :, :, nc) + dt2 * svorts(:, :, :, nc))
                    !$omp end workshare
                    !$omp end parallel
                    call field_combine_semi_spectral(svor(:, :, :, nc))
                    !$omp parallel private(iz)  default(shared)
                    !$omp do
                    do iz = 0, nz
                        svor(iz, :, :, nc) = diss * svor(iz, :, :, nc)
                    enddo
                    !$omp end do
                    !$omp end parallel
                    call field_decompose_semi_spectral(svor(:, :, :, nc))
                enddo

                call adjust_vorticity_mean

            enddo

            !Advance time:
            print *, "At time", t, "and time step", dt
            t = t + dt
        end subroutine advance

        subroutine adjust_vorticity_mean
            double precision :: wk(1:nz), savg
            integer          :: nc

            ! Ensure zero global mean horizontal vorticity conservation:
            do nc = 1, 2
                ! Cast svor_S = svor - svor_L onto the z grid as wk for kx = ky = 0:
                wk(1:nz-1) = svor(1:nz-1, 0, 0, nc)
                wk(nz) = zero
                call dst(1, nz, wk(1:nz), ztrig, zfactors)
                ! Compute average, savg (first part is the part due to svor_L):
                savg = f12 * (svor(0, 0, 0, nc) + svor(nz, 0, 0, nc)) + fnzi * sum(wk(1:nz-1))
                ! Remove from boundary values (0 & nz):
                svor(0 , 0, 0, nc) = svor(0 , 0, 0, nc) + ini_vor_mean(nc) - savg
                svor(nz, 0, 0, nc) = svor(nz, 0, 0, nc) + ini_vor_mean(nc) - savg
            enddo

        end subroutine adjust_vorticity_mean


        ! Gets the source terms for vorticity and buoyancy in spectral space.
        ! The spectral fields sbuoy and svor are all spectrally truncated.
        ! Note, vel obtained by vor2vel before calling this
        ! routine are spectrally truncated as well.
        subroutine source
#ifdef ENABLE_BUOYANCY
            double precision :: xs(0:nz, 0:nx-1, 0:ny-1)        ! db/dx or x-vtend in spectral space
            double precision :: ys(0:nz, 0:nx-1, 0:ny-1)        ! db/dy or y-vtend in spectral space
            double precision :: zs(0:nz, 0:nx-1, 0:ny-1)        ! db/dz or z-vtend in spectral space
            double precision :: dbdx(0:nz, 0:ny-1, 0:nx-1)      ! db/dx in physical space
            double precision :: dbdy(0:nz, 0:ny-1, 0:nx-1)      ! db/dy in physical space
            double precision :: dbdz(0:nz, 0:ny-1, 0:nx-1)      ! db/dz in physical space

            !--------------------------------------------------------------
            !Buoyancy source bb_t = -(u,v,w)*grad(bb): (might be computed in flux form)

            !Obtain x, y & z derivatives of buoyancy -> xs, ys, zs
            call diffx(sbuoy, xs)
            call diffy(sbuoy, ys)
            call diffz(sbuoy, zs)

            !Obtain gradient of buoyancy in physical space
            call fftxys2p(xs, dbdx)
            call fftxys2p(ys, dbdy)
            call fftxys2p(zs, dbdz)

            !Compute (u,v,w)*grad(bb) -> dbdx in physical space:
            !$omp parallel workshare
            dbdx = vel(:, :, :, 1) * dbdx &   ! u * db/dx
                 + vel(:, :, :, 2) * dbdy &   ! v * db/dy
                 + vel(:, :, :, 3) * dbdz     ! w * db/dz
            !$omp end parallel workshare

            !Convert to semi-spectral space and apply de-aliasing filter:
            call fftxyp2s(dbdx, sbuoys)
#endif

            !--------------------------------------------------------------
            !Vorticity source:

            call vorticity_tendency

        end subroutine source

        !=======================================================================

        ! Adapts the time step and computes various diagnostics
        subroutine adapt(t)
            double precision, intent(in) :: t
            double precision             :: xs(0:nz, 0:nx-1, 0:ny-1)        ! derivatives in x in spectral space
            double precision             :: ys(0:nz, 0:nx-1, 0:ny-1)        ! derivatives in y in spectral space
            double precision             :: xp(0:nz, 0:ny-1, 0:nx-1)        ! derivatives in x in physical space
#ifdef ENABLE_BUOYANCY
            double precision             :: yp(0:nz, 0:ny-1, 0:nx-1)        ! derivatives in y physical space
            double precision             :: zp(0:nz, 0:ny-1, 0:nx-1)        ! derivatives in z physical space
#endif
            double precision             :: strain(3, 3), eigs(3)
            double precision             :: dudx(0:nz, 0:ny-1, 0:nx-1)      ! du/dx in physical space
            double precision             :: dudy(0:nz, 0:ny-1, 0:nx-1)      ! du/dy in physical space
            double precision             :: dvdy(0:nz, 0:ny-1, 0:nx-1)      ! dv/dy in physical space
            double precision             :: dwdx(0:nz, 0:ny-1, 0:nx-1)      ! dw/dx in physical space
            double precision             :: dwdy(0:nz, 0:ny-1, 0:nx-1)      ! dw/dy in physical space
            double precision             :: ke, en, vormean(3)

            bfmax = zero

#ifdef ENABLE_BUOYANCY
            !Obtain x, y & z derivatives of buoyancy -> xs, ys, zs
            call diffx(sbuoy, xs)
            call diffy(sbuoy, ys)

            call field_combine_physical(sbuoy, yp)
            call diffz(yp, zp)

            !Obtain gradient of buoyancy in physical space -> xp, yp, zp
            call fftxys2p(xs, xp)
            call fftxys2p(ys, yp)

            !Compute (db/dx)^2 + (db/dy)^2 + (db/dz)^2 -> xp in physical space:
            !$omp parallel
            !$omp workshare
            xp = xp ** 2 + yp ** 2 + zp ** 2

            !Maximum buoyancy frequency:
            bfmax = dsqrt(dsqrt(maxval(xp)))
            !$omp end workshare
            !$omp end parallel
#endif

            !$omp parallel
            !$omp workshare

            !Compute enstrophy: (reuse xp)
            xp = vor(:, :, :, 1) ** 2 + vor(:, :, :, 2) ** 2 + vor(:, :, :, 3) ** 2

            !Maximum vorticity magnitude:
            vortmax = dsqrt(maxval(xp))

            !R.m.s. vorticity:
            vortrms = dsqrt(ncelli*(f12*sum(xp(0, :, :)+xp(nz, :, :))+sum(xp(1:nz-1, :, :))))
            !$omp end workshare
            !$omp end parallel

            !Characteristic vorticity,  <vor^2>/<|vor|> for |vor| > vor_rms:
            vorl1 = small
            vorl2 = zero
            !$omp parallel private(ix, iy, iz, vortmp1, vortmp2, vortmp3)  default(shared)
            !$omp do reduction(+:vorl1,vorl2) collapse(3)
            do ix = 0, nx-1
                do iy = 0, ny-1
                    do iz = 1, nz
                        vortmp1 = f12 * abs(vor(iz-1, iy, ix, 1) + vor(iz, iy, ix, 1))
                        vortmp2 = f12 * abs(vor(iz-1, iy, ix, 2) + vor(iz, iy, ix, 2))
                        vortmp3 = f12 * abs(vor(iz-1, iy, ix, 3) + vor(iz, iy, ix, 3))
                        if (vortmp1 + vortmp2 + vortmp3 .gt. vortrms) then
                            vorl1 = vorl1 + vortmp1 + vortmp2 + vortmp3
                            vorl2 = vorl2 + vortmp1 ** 2 + vortmp2 ** 2 + vortmp3 ** 2
                        endif
                    enddo
                enddo
            enddo
            !$omp end do
            !$omp end parallel
            vorch = vorl2 / vorl1

            vormean = get_mean_vorticity()

            ! Save vorticity diagnostics to vorticity.asc:
            write(WRITE_VOR, '(1x,f13.6,6(1x,1p,e14.7))') t , vortmax, vortrms, vorch, vormean

            ! Save energy and enstrophy
            ke = get_kinetic_energy()
            en = get_enstrophy()
            write(WRITE_ECOMP, '(1x,f13.6,2(1x,1p,e14.7))') t , ke, en

            !
            ! velocity strain
            !

            ! du/dx
            call diffx(svel(:, :, :, 1), xs)
            call fftxys2p(xs, dudx)

            ! du/dy
            call diffy(svel(:, :, :, 1), ys)
            call fftxys2p(ys, dudy)

            ! dw/dx
            call diffx(svel(:, :, :, 3), xs)
            call fftxys2p(xs, dwdx)

            ! dv/dy
            call diffy(svel(:, :, :, 2), ys)
            call fftxys2p(ys, dvdy)

            ! dw/dy
            call diffy(svel(:, :, :, 3), ys)
            call fftxys2p(ys, dwdy)

            ! find largest stretch -- this corresponds to largest
            ! eigenvalue over all local symmetrised strain matrices.
            ggmax = epsilon(ggmax)
            !$omp parallel private(ix, iy, iz, strain, eigs)  default(shared)
            !$omp do reduction(max:ggmax) collapse(3)
            do ix = 0, nx-1
                do iy = 0, ny-1
                    do iz = 0, nz
                        ! get local symmetrised strain matrix, i.e. 1/ 2 * (S + S^T)
                        ! where
                        !     /u_x u_y u_z\
                        ! S = |v_x v_y v_z|
                        !     \w_x w_y w_z/
                        ! with u_* = du/d* (also derivatives of v and w).
                        ! The derivatives dv/dx, du/dz, dv/dz and dw/dz are calculated
                        ! with vorticity or the assumption of incompressibility
                        ! (du/dx + dv/dy + dw/dz = 0):
                        !    dv/dx = \omegaz + du/dy
                        !    du/dz = \omegay + dw/dx
                        !    dv/dz = dw/dy - \omegax
                        !    dw/dz = - (du/dx + dv/dy)
                        !
                        !                         /  2 * u_x  u_y + v_x u_z + w_x\
                        ! 1/2 * (S + S^T) = 1/2 * |u_y + v_x   2 * v_y  v_z + w_y|
                        !                         \u_z + w_x  v_z + w_y   2 * w_z/
                        !
                        ! S11 = du/dx
                        ! S12 = 1/2 * (du/dy + dv/dx) = 1/2 * (2 * du/dy + \omegaz) = du/dy + 1/2 * \omegaz
                        ! S13 = 1/2 * (du/dz + dw/dx) = 1/2 * (\omegay + 2 * dw/dx) = 1/2 * \omegay + dw/dx
                        ! S22 = dv/dy
                        ! S23 = 1/2 * (dv/dz + dw/dy) = 1/2 * (2 * dw/dy - \omegax) = dw/dy - 1/2 * \omegax
                        ! S33 = dw/dz = - (du/dx + dv/dy)
                        strain(1, 1) = dudx(iz, iy, ix)                            ! S11
                        strain(1, 2) = dudy(iz, iy, ix) + f12 * vor(iz, iy, ix, 3) ! S12
                        strain(1, 3) = dwdx(iz, iy, ix) + f12 * vor(iz, iy, ix, 2) ! S13
                        strain(2, 2) = dvdy(iz, iy, ix)                            ! S22
                        strain(2, 3) = dwdy(iz, iy, ix) - f12 * vor(iz, iy, ix, 1) ! S23
                        strain(3, 3) = -(dudx(iz, iy, ix) + dvdy(iz, iy, ix))      ! S33

                        ! calculate its eigenvalues. The Jacobi solver
                        ! requires the upper triangular matrix only.
                        call jacobi_eigenvalues(strain, eigs)

                        ! we must take the largest eigenvalue in magnitude (absolute value)
                        ggmax = max(ggmax, maxval(abs(eigs)))
                    enddo
                enddo
            enddo
            !$omp end do

            !$omp workshare

            !Maximum speed:
            velmax = dsqrt(maxval(vel(:, :, :, 1) ** 2   &
                                + vel(:, :, :, 2) ** 2   &
                                + vel(:, :, :, 3) ** 2))

            !$omp end workshare
            !$omp end parallel

            !Choose new time step:
            dt = min(time%alpha / (ggmax + small),  &
                     time%alpha / (bfmax + small),  &
                     cflpf / (velmax + small),      &
                     time%limit - t)

            !Update value of dt/2:
            dt2 = f12 * dt

            !---------------------------------------------------------------------
            if (viscosity%nnu .eq. 1) then
                !Update diffusion operator used in time stepping:
                dfac = dt
                !$omp parallel
                !$omp workshare
                diss = one / (one + dfac * hdis)
                !$omp end workshare
                !$omp end parallel
                !(see inversion_utils.f90)
            else
                !Update hyperdiffusion operator used in time stepping:
                dfac = vorch * dt
                !$omp parallel
                !$omp workshare
                diss = one / (one + dfac * hdis)
                !$omp end workshare
                !$omp end parallel
                !(see inversion_utils.f90)
             endif

        end subroutine adapt

end module advance_mod
