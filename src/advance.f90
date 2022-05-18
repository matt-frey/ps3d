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
    use inversion_mod, only : vor2vel
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
            double precision                :: sbuoys(0:nz, 0:nx-1, 0:ny-1) ! source of buoyancy (spectral)
            double precision                :: vortsi(0:nz, 0:nx-1, 0:ny-1, 3)
            double precision                :: vortsm(0:nz, 0:nx-1, 0:ny-1, 3)
            double precision                :: svorts(0:nz, 0:nx-1, 0:ny-1, 3)  ! source of vorticiy (spectral)
            integer                         :: nc

            !-------------------------------------------------------------------
            !Invert vorticity for velocity at current time level, say t=t^n:
            call vor2vel(svortg, velog)

            !Adapt the time step and save various diagnostics each time step:
            call adapt(t)

!             !Possibly save field data (gsave is set by adapt):
!             if (gsave) call savegrid

            !------------------------------------------------------------------
            !Start with a guess for F^{n+1} for all fields:

            !Calculate the source terms (sbuoys, svorts) for buoyancy (sbuoyg) and
            !vorticity (svortg) in spectral space:
            call source(sbuoys, svorts)

            !Initialise iteration (dt = dt/4 below):
            bsi = sbuoyg
            bsm = sbuoyg + dt4 * sbuoys
            sbuoyg = diss * (bsm + dt4 * sbuoys) - bsi
            vortsi = svortg
            vortsm = svortg + dt4 * svorts

            do nc = 1, 3
                svortg(:, :, :, nc) = diss * (vortsm(:, :, :, nc) + dt4 * svorts(:, :, :, nc)) &
                                    - vortsi(:, :, :, nc)
            enddo
            !diss is related to the hyperdiffusive operator (see end of adapt)

            !------------------------------------------------------------------
            !Iterate to improve estimates of F^{n+1}:
            do iter = 1, niter
                !Perform inversion at t^{n+1} from estimated quantities:
                call vor2vel(svortg, velog)

                !Calculate the source terms (sbuoys,svorts):
                call source(sbuoys, svorts)

                !Update fields:
                sbuoyg = diss * (bsm + dt4 * sbuoys) - bsi

                do nc = 1, 3
                    svortg(:, :, :, nc) = diss * (vortsm(:, :, :, nc) + dt4 * svorts(:, :, :, nc)) &
                                        - vortsi(:, :, :, nc)
                enddo
            enddo

            !Advance time:
            t = t + dt
        end subroutine advance


        ! Gets the source terms for vorticity and buoyancy in spectral space.
        ! The spectral fields sbuoyg and svortg are all spectrally truncated.
        ! Note, velog obtained by vor2vel before calling this
        ! routine are spectrally truncated as well.
        subroutine source(sbuoys, svorts)
            double precision, intent(inout) :: sbuoys(0:nz, 0:nx-1, 0:ny-1)    ! in spectral space
            double precision, intent(inout) :: svorts(0:nz, 0:nx-1, 0:ny-1, 3) ! in spectral space
            double precision                :: dbdxs(0:nz, 0:ny-1, 0:nx-1)     ! db/dx in spectral space
            double precision                :: dbdys(0:nz, 0:ny-1, 0:nx-1)     ! db/dy in spectral space
            double precision                :: dbdzs(0:nz, 0:ny-1, 0:nx-1)     ! db/dz in spectral space
            double precision                :: dbdx(0:nz, 0:ny-1, 0:nx-1)      ! db/dx in physical space
            double precision                :: dbdy(0:nz, 0:ny-1, 0:nx-1)      ! db/dy in physical space
            double precision                :: dbdz(0:nz, 0:ny-1, 0:nx-1)      ! db/dz in physical space

            !--------------------------------------------------------------
            !Buoyancy source bb_t = -(u,v,w)*grad(bb):

            !Obtain x, y & z derivatives of buoyancy -> dbdxs, dbdys, dbdzs
            call diffx(sbuoys, dbdxs)
            call diffy(sbuoys, dbdys)
            call diffz(sbuoys, dbdzs)

            !Store spectral db/dx, db/dy and db/dz in svorts for use in vorticity source below:
            svorts(:, :, :, 1) = dbdxs
            svorts(:, :, :, 2) = dbdys
            svorts(:, :, :, 3) = dbdzs

            !Obtain gradient of buoyancy in physical space
            call fftxys2p(dbdxs, dbdx)
            call fftxys2p(dbdys, dbdy)
            call fftxys2p(dbdzs, dbdz)

            !Compute (u,v,w)*grad(bb) -> dbdx in physical space:
            dbdx = velog(:, :, :, 1) * dbdx &   ! u * db/dx
                 + velog(:, :, :, 2) * dbdy &   ! v * db/dy
                 + velog(:, :, :, 3) * dbdz     ! w * db/dz

            !Convert to spectral space as sbuoys and apply de-aliasing filter:
            call fftxyp2s(dbdx, sbuoys)

            do iz = 0, nz
                sbuoys(iz, :, :) = -filt * sbuoys(iz, :, :)
            enddo

            !--------------------------------------------------------------
            !Vorticity source:

            call vorticity_tendency(vortg, velog, buoyg, vtend)

            !Convert to spectral space and apply de-aliasing filter:
            call fftxyp2s(vtend(:, :, :, 1), dbdxs)
            call fftxyp2s(vtend(:, :, :, 2), dbdys)
            call fftxyp2s(vtend(:, :, :, 3), dbdzs)

            do iz = 0, nz
                svorts(iz, :, :, 1) = svorts(iz, :, :, 1) - filt * dbdxs(iz, :, :)
                svorts(iz, :, :, 2) = svorts(iz, :, :, 2) - filt * dbdys(iz, :, :)
                svorts(iz, :, :, 3) = svorts(iz, :, :, 3) - filt * dbdzs(iz, :, :)
            enddo

        end subroutine source

        !=======================================================================

        ! Adapts the time step and computes various diagnostics
        subroutine adapt(t)
            double precision, intent(in) :: t
            !For defining the max strain & buoyancy frequency based time step:
            double precision, parameter  :: alpha = 0.1d0
            !Note: EPIC-2D paper recommends alpha = 0.2 for ls-rk4 method
            double precision             :: dbdxs(0:nz, 0:ny-1, 0:nx-1)     ! db/dx in spectral space
            double precision             :: dbdys(0:nz, 0:ny-1, 0:nx-1)     ! db/dy in spectral space
            double precision             :: dbdzs(0:nz, 0:ny-1, 0:nx-1)     ! db/dz in spectral space
            double precision             :: dbdx(0:nz, 0:ny-1, 0:nx-1)      ! db/dx in physical space
            double precision             :: dbdy(0:nz, 0:ny-1, 0:nx-1)      ! db/dy in physical space
            double precision             :: dbdz(0:nz, 0:ny-1, 0:nx-1)      ! db/dz in physical space

            !Local variables (physical):
            double precision :: px(0:nz, 0:ny-1, 0:nx-1), py(nz, 0:ny-1, 0:nx-1)
            double precision :: fp(0:nz, 0:ny-1, 0:nx-1)

            !Local variables (spectral):
            double precision :: sx(0:nz, 0:nx-1, 0:ny-1), sy(nz, 0:nx-1, 0:ny-1)
            double precision :: fs(0:nz, 0:nx-1, 0:ny-1)

            !Obtain x, y & z derivatives of buoyancy -> dbdxs, dbdys, dbdzs
            call diffx(sbuoys, dbdxs)
            call diffy(sbuoys, dbdys)
            call diffz(sbuoys, dbdzs)

            !Obtain gradient of buoyancy in physical space
            call fftxys2p(dbdxs, dbdx)
            call fftxys2p(dbdys, dbdy)
            call fftxys2p(dbdzs, dbdz)


            !Compute (db/dx)^2 + (db/dy)^2 + (db/dz)^2 -> dbdx in physical space:
            dbdx = dbdx ** 2 + dbdy ** 2 + dbdz ** 2

            !Maximum buoyancy frequency:
            bfmax = sqrt(sqrt(maxval(dbdx)))

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
                !(see inversion_utils.f90 and parameters.f90).
            else
                !Update hyperdiffusion operator used in time stepping:
                dfac = zzch * dt / two
                !zzch is the characteristic vorticity defined above.
                diss = two / (one + dfac * hdis)
                !hdis = C*(K/K_max)^{2p} where K^2 = k_x^2+k_y^2, p is the order,
                !K_max is the maximum x or y wavenumber and C is a dimensionless
                !prefactor (see inversion_utils.f90 and parameters.f90 where C = prediss).
            endif

        end subroutine adapt

end module advance_mod
