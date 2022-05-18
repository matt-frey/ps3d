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
    use options, only : time, nnu
    use constants
    use parameters, only : nx, ny, nz, glmin, cflpf, ncelli
    use spectral, only : vor2vel
    use inversion_utils
    use fields
    implicit none

    ! Number of iterations of above scheme:
    integer, parameter:: niter = 2

    double precision :: dt, dt4

    !Diagnostic quantities:
    double precision :: bfmax, zzmax, zzrms, ggmax, uumax, dfac
    double precision :: zztmp, zzl1, zzl2, zzch
    integer          :: ix, iy, iz

    contains

        subroutine advance(t)
            double precision, intent(inout) :: t
            integer                         :: iter
            ! Spectral fields needed in time stepping:
            double precision                :: bsi(0:nz, 0:nx-1, 0:ny-1)
            double precision                :: bsm(0:nz, 0:nx-1, 0:ny-1)
            double precision                :: sbs(0:nz, 0:nx-1, 0:ny-1)
            double precision                :: zsi(0:nz, 0:nx-1, 0:ny-1)
            double precision                :: zsm(0:nz, 0:nx-1, 0:ny-1)
            double precision                :: szs(0:nz, 0:nx-1, 0:ny-1)

            !-------------------------------------------------------------------
            !Invert vorticity for velocity at current time level, say t=t^n:
            call vor2vel(xi, eta, zeta,  uu, vv, ww,  velgradg)  !call main_invert(zs,uu,vv,zz) !zs is spectral here

            !Adapt the time step and save various diagnostics each time step:
            call adapt(t)

!             !Possibly save field data (gsave is set by adapt):
!             if (gsave) call savegrid

            !------------------------------------------------------------------
            !Start with a guess for F^{n+1} for all fields:

            !Calculate the source terms (sbs, szs) for buoyancy (sbuoyg) and
            !vorticity (svortg) in spectral space:
            call source(sbs, szs)

            !Initialise iteration (dt = dt/4 below):
            bsi = sbuoyg
            bsm = sbuoyg + dt4 * sbs
            sbuoyg = diss * (bsm + dt4 * sbs) - bsi
            zsi = svortg(:, :, :, 1)
            zsm = svortg(:, :, :, 1) + dt4 * szs
            svortg(:, :, :, 1) = diss * (zsm + dt4 * szs) - zsi
            !diss is related to the hyperdiffusive operator (see end of adapt)

            !------------------------------------------------------------------
            !Iterate to improve estimates of F^{n+1}:
            do iter = 1, niter
                !Perform inversion at t^{n+1} from estimated quantities:
                call vor2vel(xi, eta, zeta,  uu, vv, ww,  velgradg) !main_invert(zs,uu,vv,zz) !zs is spectral here

                !Calculate the source terms (sbs,szs):
                call source(sbs, szs)

                !Update fields:
                sbuoyg = diss * (bsm + dt4 * sbs) - bsi
                svortg(:, :, :, 1) = diss * (zsm + dt4 * szs) - zsi
            enddo

            !Advance time:
            t = t + dt
        end subroutine advance


        ! Gets the source terms for vorticity and buoyancy in spectral space.
        ! The spectral fields sbuoyg and svortg are all spectrally truncated.
        ! Note, uu and vv obtained by main_invert before calling this
        ! routine are spectrally truncated as well.
        subroutine source(sbs, szs)
            double precision, intent(inout) :: sbs(0:nz, 0:nx-1, 0:ny-1) ! in spectral space
            double precision, intent(inout) :: szs(0:nz, 0:nx-1, 0:ny-1) ! in spectral space
            double precision                :: px(0:nz, 0:ny-1, 0:nx-1)  ! in physical space
            double precision                :: py(nz, 0:ny-1, 0:nx-1)    ! in physical space
            double precision                :: sx(0:nz, 0:nx-1, 0:ny-1)  ! in spectral space
            double precision                :: sy(nz, 0:nx-1, 0:ny-1)    ! in spectral space

            !--------------------------------------------------------------
            !Buoyancy source bb_t = -(u,v)*grad(bb):

            !Obtain x & y derivatives of buoyancy -> px, py (physical):
            call xderiv_fc(nx, ny, hrkx, sbuoyg, sx)

            !Store spectral db/dx in szs for use in vorticity source below:
            szs = sx

            call spctop_fc(nx, ny, sx, px, xfactors, yfactors, xtrig, ytrig)

            call yderiv_fc(nx, ny, rky, sbuoyg, sy)

            call spctop_fs(nx, ny, sy, py, xfactors, yfactors, xtrig, ytrig)

            !Compute (u, v)*grad(bb) -> px in physical space:
            px(0, :, :) = uu(0, :, :) * px(0, :, :)

            px(1:nz-1, :, :) = uu(1:nz-1, :, :) * px(1:nz-1, :, :)  &
                             + vv(1:nz-1, :, :) * py(1:nz-1, :, :)

            px(nz, :, :) = uu(nz, :, :) * px(nz, :, :)

            !Convert to spectral space as sbs and apply de-aliasing filter:
            call ptospc_fc(nx, ny, px, sbs, xfactors, yfactors, xtrig, ytrig)

            do iz = 0, nz
                sbs(iz, :, :) = -filt * sbs(iz, :, :)
            enddo

            !--------------------------------------------------------------
            !Vorticity source zz_t = bb_x - (u, v)*grad(zz):

            !Obtain x & y derivatives of vorticity:
            call xderiv_fc(nx, ny, hrkx, svortg(:, :, :, 1), sx)
            call spctop_fc(nx, ny, sx, px, xfactors, yfactors, xtrig, ytrig)
            call yderiv_fc(nx, ny, rky, svortg(:, :, :, 1), sy)
            call spctop_fs(nx, ny, sy, py, xfactors, yfactors, xtrig, ytrig)

            !Compute (u, v)*grad(zz) -> px in physical space:
            px(0, :, :) = uu(0, :, :) * px(0, :, :)

            px(1:nz-1, :, :) = uu(1:nz-1, :, :) * px(1:nz-1, :, :) &
                             + vv(1:nz-1, :, :) * py(1:nz-1, :, :)

            px(nz, :, :) = uu(nz, :, :) * px(nz, :, :)

            !Convert to spectral space as sx and apply de-aliasing filter:
            call ptospc_fc(nx, ny, px, sx, xfactors, yfactors, xtrig, ytrig)

            do iz = 0, nz
                szs(iz, :, :) = szs(iz, :, :) - filt * sx(iz, :, :)
            enddo

        end subroutine source

        !=======================================================================

        ! Adapts the time step and computes various diagnostics
        subroutine adapt(t)
            double precision, intent(in) :: t
            !For defining the max strain & buoyancy frequency based time step:
            double precision, parameter:: alpha = 0.1d0
            !Note: EPIC-2D paper recommends alpha = 0.2 for ls-rk4 method

            !Local variables (physical):
            double precision :: px(0:nz, 0:ny-1, 0:nx-1), py(nz, 0:ny-1, 0:nx-1)
            double precision :: fp(0:nz, 0:ny-1, 0:nx-1)

            !Local variables (spectral):
            double precision :: sx(0:nz, 0:nx-1, 0:ny-1), sy(nz, 0:nx-1, 0:ny-1)
            double precision :: fs(0:nz, 0:nx-1, 0:ny-1)

            !----------------------------------------------------------
            !Obtain x & y derivatives of buoyancy -> px, py (physical):
            call xderiv_fc(nx, ny, hrkx, sbuoyg, sx)
            call spctop_fc(nx, ny, sx, px, xfactors, yfactors, xtrig, ytrig)
            call yderiv_fc(nx, ny, rky, sbuoyg, sy)
            call spctop_fs(nx, ny, sy, py, xfactors, yfactors, xtrig, ytrig)

            !Compute px^2 + py^2 -> px in physical space:
            px(0, :, :) = px(0, :, :) ** 2
            px(1:nz-1, :, :) = px(1:nz-1, :, :) ** 2 + py(1:nz-1, :, :) ** 2
            px(nz, :, :) = px(nz, :, :) ** 2

            !Maximum buoyancy frequency:
            bfmax = sqrt(sqrt(maxval(px)))

            !Maximum vorticity:
            px = zeta ** 2
            zzmax = sqrt(maxval(px))

            !R.m.s. vorticity:
            zzrms = sqrt(ncelli*(f12*sum(px(0, :, :)+px(nz, :, :))+sum(px(1:nz-1, :, :))))

            !Characteristic vorticity,  <zz^2>/<|zz|> for |zz| > zz_rms:
            zzl1 = small
            zzl2 = zero
            do ix = 0, nx-1
                do iy = 0, ny-1
                    do iz = 1, nz
                        zztmp = f12 * (zeta(iz-1, iy, ix) + zeta(iz, iy, ix))
                        if (abs(zztmp) .gt. zzrms) then
                            zzl1 = zzl1 + abs(zztmp)
                            zzl2 = zzl2 + zztmp ** 2
                        endif
                    enddo
                enddo
            enddo
            zzch = zzl2 / zzl1

            !Compute x derivative of velocity components:
            fp = uu
            call ptospc_fc(nx, ny, fp, fs, xfactors, yfactors, xtrig, ytrig)
            call xderiv_fc(nx, ny, hrkx, fs, sx)
            call spctop_fc(nx, ny, sx, fp, xfactors, yfactors, xtrig, ytrig)
            !fp = u_x
            px = vv
            call ptospc_fc(nx, ny, px, fs, xfactors, yfactors, xtrig, ytrig)
            call xderiv_fc(nx, ny, hrkx, fs, sx)
            call spctop_fc(nx, ny, sx, px, xfactors, yfactors, xtrig, ytrig)
            !px = v_x

            !Strain rate squared, u_x^2 + (v_x - zz/2)^2:
            fp = fp ** 2 + (px - f12 * xi) ** 2

            !Maximum strain rate:
            ggmax = sqrt(maxval(fp))

            !Maximum speed:
            uumax = sqrt(maxval(uu ** 2 + vv ** 2))

            !Choose new time step:
            dt = min(alpha / (ggmax + small), alpha / (bfmax + small), cflpf / (uumax + small), time%limit - t)

            !Update value of dt/4:
            dt4 = dt / four

            !---------------------------------------------------------------------
            if (nnu .eq. 1) then
                !Update diffusion operator used in time stepping:
                dfac = dt / two
                diss = two / (one + dfac * hdis)
                !hdis = nu*(k_x^2+k_y^2) where nu is the viscosity coefficient
                !(see spectral.f90 and parameters.f90).
            else
                !Update hyperdiffusion operator used in time stepping:
                dfac = zzch * dt / two
                !zzch is the characteristic vorticity defined above.
                diss = two / (one + dfac * hdis)
                !hdis = C*(K/K_max)^{2p} where K^2 = k_x^2+k_y^2, p is the order,
                !K_max is the maximum x or y wavenumber and C is a dimensionless
                !prefactor (see spectral.f90 and parameters.f90 where C = prediss).
            endif

            !---------------------------------------------------------------------
            !Save |u|_max,  N_max and gamma_max to monitor.asc:
            write(22, '(1x,f13.6,3(1x,1p,e14.7))') t, uumax, bfmax, ggmax

            !Save vorticity diagnostics to vorticity.asc:
            write(23, '(1x,f13.6,3(1x,1p,e14.7))') t, zzmax, zzrms, zzch

            !---------------------------------------------------------------------
            !Set flag to save data:
!             gsave=(int(t/tgsave) .eq. igrids)
            !Gridded data will be saved at time t if gsave is true.
        end subroutine adapt

end module advance_mod
