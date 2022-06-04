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
    use inversion_mod, only : vor2vel, vorticity_tendency
    use inversion_utils
    use utils, only : write_step
    use fields
    use jacobi, only : jacobi_eigenvalues
    implicit none

    integer :: advance_timer

    integer, parameter :: WRITE_VOR = 1234
    integer, parameter :: WRITE_ECOMP = 1235

    ! Number of iterations of above scheme:
    integer, parameter:: niter = 2

    double precision :: dt, dt4

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
            double precision                :: bsi(0:nz, 0:nx-1, 0:ny-1)
            double precision                :: bsm(0:nz, 0:nx-1, 0:ny-1)
            double precision                :: sbuoys(0:nz, 0:nx-1, 0:ny-1)     ! source of buoyancy (spectral)
            double precision                :: vortsi(0:nz, 0:nx-1, 0:ny-1, 3)
            double precision                :: vortsm(0:nz, 0:nx-1, 0:ny-1, 3)
            double precision                :: svorts(0:nz, 0:nx-1, 0:ny-1, 3)  ! source of vorticity (spectral)

            !-------------------------------------------------------------------
            !Invert vorticity for velocity at current time level, say t=t^n:
            !Also, returns vorticity in physical space for use everywhere
            call vor2vel(svortg, vortg, svelog, velog)

            !Adapt the time step
            call adapt(t, sbuoyg, velog, svelog)

            !Write fields
            call write_step(t)

            !------------------------------------------------------------------
            !Start with a guess for F^{n+1} for all fields:

            !Calculate the source terms (sbuoys, svorts) for buoyancy (sbuoyg) and
            !vorticity (svortg) in spectral space:
            call source(sbuoys, svorts)

            !Initialise iteration (dt = dt/4 below):
            bsi = sbuoyg
            bsm = sbuoyg + dt4 * sbuoys
            do iz = 0, nz
                sbuoyg(iz, :, :) = diss * (bsm(iz, :, :) + dt4 * sbuoys(iz, :, :)) - bsi(iz, :, :)
            enddo



            vortsi = svortg
            vortsm = svortg + dt4 * svorts

            do nc = 1, 3
                do iz = 0, nz
                    svortg(iz, :, :, nc) = diss * (vortsm(iz, :, :, nc) + dt4 * svorts(iz, :, :, nc)) &
                                         - vortsi(iz, :, :, nc)
                enddo
            enddo
            !diss is related to the hyperdiffusive operator (see end of adapt)

            !------------------------------------------------------------------
            !Iterate to improve estimates of F^{n+1}:
            do iter = 1, niter
                !Perform inversion at t^{n+1} from estimated quantities:
                call vor2vel(svortg, vortg, svelog, velog)

                !Calculate the source terms (sbuoys,svorts):
                call source(sbuoys, svorts)

                !Update fields:
                do iz = 0, nz
                    sbuoyg(iz, :, :) = diss * (bsm(iz, :, :) + dt4 * sbuoys(iz, :, :)) - bsi(iz, :, :)
                enddo

                do nc = 1, 3
                   do iz = 0, nz
                        svortg(iz, :, :, nc) = diss * (vortsm(iz, :, :, nc) + dt4 * svorts(iz, :, :, nc)) &
                                             - vortsi(iz, :, :, nc)
                    enddo
                enddo
            enddo

            !Advance time:
            print *, "At time", t, "and time step", dt
            t = t + dt
        end subroutine advance


        ! Gets the source terms for vorticity and buoyancy in spectral space.
        ! The spectral fields sbuoyg and svortg are all spectrally truncated.
        ! Note, velog obtained by vor2vel before calling this
        ! routine are spectrally truncated as well.
        subroutine source(sbuoys, svorts)
            double precision, intent(inout) :: sbuoys(0:nz, 0:nx-1, 0:ny-1)    ! in spectral space
            double precision, intent(inout) :: svorts(0:nz, 0:nx-1, 0:ny-1, 3) ! in spectral space
            double precision                :: xs(0:nz, 0:nx-1, 0:ny-1)        ! db/dx or x-vtend in spectral space
            double precision                :: ys(0:nz, 0:nx-1, 0:ny-1)        ! db/dy or y-vtend in spectral space
            double precision                :: zs(0:nz, 0:nx-1, 0:ny-1)        ! db/dz or z-vtend in spectral space
            double precision                :: dbdx(0:nz, 0:ny-1, 0:nx-1)      ! db/dx in physical space
            double precision                :: dbdy(0:nz, 0:ny-1, 0:nx-1)      ! db/dy in physical space
            double precision                :: dbdz(0:nz, 0:ny-1, 0:nx-1)      ! db/dz in physical space
            double precision                :: svtend(0:nz, 0:nx-1, 0:ny-1, 3) ! vtend in spectral space

            !--------------------------------------------------------------
            !Buoyancy source bb_t = -(u,v,w)*grad(bb): (might be computed in flux form)

            !Obtain x, y & z derivatives of buoyancy -> xs, ys, zs
            call diffx(sbuoyg, xs)
            call diffy(sbuoyg, ys)
            call diffz(sbuoyg, zs)

            !Store spectral db/dx and db/dy in svorts for use in vorticity source below:
            call fftss2fs(ys, svorts(:, :, :, 1))
            call fftss2fs(xs, svorts(:, :, :, 2))

            !Obtain gradient of buoyancy in physical space
            call fftxys2p(xs, dbdx)
            call fftxys2p(ys, dbdy)
            call fftxys2p(zs, dbdz)

            !Compute (u,v,w)*grad(bb) -> dbdx in physical space:
            dbdx = velog(:, :, :, 1) * dbdx &   ! u * db/dx
                 + velog(:, :, :, 2) * dbdy &   ! v * db/dy
                 + velog(:, :, :, 3) * dbdz     ! w * db/dz

            !Convert to semi-spectral space and apply de-aliasing filter:
            call fftxyp2s(dbdx, sbuoys)

            do iz = 0, nz
                sbuoys(iz, :, :) = -filt * sbuoys(iz, :, :)
            enddo

            !--------------------------------------------------------------
            !Vorticity source (excluding buoyancy effects):

            call vorticity_tendency(svortg, velog, vortg, svtend)

            !Add filtered vorticity tendency to vorticity source:
            do iz = 0, nz
                svorts(iz, :, :, 1) =   svorts(iz, :, :, 1) + filt * svtend(iz, :, :, 1)
                svorts(iz, :, :, 2) = - svorts(iz, :, :, 2) + filt * svtend(iz, :, :, 2)
                svorts(iz, :, :, 3) =                       + filt * svtend(iz, :, :, 3)
            enddo

        end subroutine source

        !=======================================================================

        ! Adapts the time step and computes various diagnostics
        subroutine adapt(t, sbuoyg, velog, svelog)
            double precision, intent(in) :: t
            double precision, intent(in) :: sbuoyg(0:nz, 0:nx-1, 0:ny-1)
            double precision, intent(in) :: velog(0:nz, 0:ny-1, 0:nx-1, 3)
            double precision, intent(in) :: svelog(0:nz, 0:nx-1, 0:ny-1, 3)
            double precision             :: xs(0:nz, 0:nx-1, 0:ny-1)        ! derivatives in x in spectral space
            double precision             :: ys(0:nz, 0:nx-1, 0:ny-1)        ! derivatives in y in spectral space
            double precision             :: zs(0:nz, 0:nx-1, 0:ny-1)        ! derivatives in z in spectral space
            double precision             :: xp(0:nz, 0:ny-1, 0:nx-1)        ! derivatives in x in physical space
            double precision             :: yp(0:nz, 0:ny-1, 0:nx-1)        ! derivatives in y physical space
            double precision             :: zp(0:nz, 0:ny-1, 0:nx-1)        ! derivatives in z physical space
            double precision             :: strain(3, 3), eigs(3)
            double precision             :: dudx(0:nz, 0:ny-1, 0:nx-1)      ! du/dx in physical space
            double precision             :: dudy(0:nz, 0:ny-1, 0:nx-1)      ! du/dy in physical space
            double precision             :: dvdy(0:nz, 0:ny-1, 0:nx-1)      ! dv/dy in physical space
            double precision             :: dwdx(0:nz, 0:ny-1, 0:nx-1)      ! dw/dx in physical space
            double precision             :: dwdy(0:nz, 0:ny-1, 0:nx-1)      ! dw/dy in physical space
            double precision             :: ke, en, vormean(3)


            !Obtain x, y & z derivatives of buoyancy -> xs, ys, zs
            call diffx(sbuoyg, xs)
            call diffy(sbuoyg, ys)
            call diffz(sbuoyg, zs)

            !Obtain gradient of buoyancy in physical space -> xp, yp, zp
            call fftxys2p(xs, xp)
            call fftxys2p(ys, yp)
            call fftxys2p(zs, zp)


            !Compute (db/dx)^2 + (db/dy)^2 + (db/dz)^2 -> xp in physical space:
            xp = xp ** 2 + yp ** 2 + zp ** 2

            !Maximum buoyancy frequency:
            bfmax = sqrt(sqrt(maxval(xp)))

            !Compute enstrophy: (reuse xp)
            xp = vortg(:, :, :, 1) ** 2 + vortg(:, :, :, 2) ** 2 + vortg(:, :, :, 3) ** 2

            !Maximum vorticity magnitude:
            vortmax = sqrt(maxval(xp))

            !R.m.s. vorticity:
            vortrms = sqrt(ncelli*(f12*sum(xp(0, :, :)+xp(nz, :, :))+sum(xp(1:nz-1, :, :))))

            !Characteristic vorticity,  <vor^2>/<|vor|> for |vor| > vor_rms:
            vorl1 = small
            vorl2 = zero
            do ix = 0, nx-1
                do iy = 0, ny-1
                    do iz = 1, nz
                        vortmp1 = f12 * abs(vortg(iz-1, iy, ix, 1) + vortg(iz, iy, ix, 1))
                        vortmp2 = f12 * abs(vortg(iz-1, iy, ix, 2) + vortg(iz, iy, ix, 2))
                        vortmp3 = f12 * abs(vortg(iz-1, iy, ix, 3) + vortg(iz, iy, ix, 3))
                        if (vortmp1 + vortmp2 + vortmp3 .gt. vortrms) then
                            vorl1 = vorl1 + vortmp1 + vortmp2 + vortmp3
                            vorl2 = vorl2 + vortmp1 ** 2 + vortmp2 ** 2 + vortmp3 ** 2
                        endif
                    enddo
                enddo
            enddo
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
            call diffx(svelog(:, :, :, 1), xs)
            call fftxys2p(xs, dudx)

            ! du/dy
            call diffy(svelog(:, :, :, 1), ys)
            call fftxys2p(ys, dudy)

            ! dw/dx
            call diffx(svelog(:, :, :, 3), xs)
            call fftxys2p(xs, dwdx)

            ! dv/dy
            call diffy(svelog(:, :, :, 2), ys)
            call fftxys2p(ys, dvdy)

            ! dw/dy
            call diffy(svelog(:, :, :, 3), ys)
            call fftxys2p(ys, dwdy)

            ! find largest stretch -- this corresponds to largest
            ! eigenvalue over all local symmetrised strain matrices.
            ggmax = epsilon(ggmax)
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
                        strain(1, 1) = dudx(iz, iy, ix)                              ! S11
                        strain(1, 2) = dudy(iz, iy, ix) + f12 * vortg(iz, iy, ix, 3) ! S12
                        strain(1, 3) = dwdx(iz, iy, ix) + f12 * vortg(iz, iy, ix, 2) ! S13
                        strain(2, 2) = dvdy(iz, iy, ix)                              ! S22
                        strain(2, 3) = dwdy(iz, iy, ix) - f12 * vortg(iz, iy, ix, 1) ! S23
                        strain(3, 3) = -(dudx(iz, iy, ix) + dvdy(iz, iy, ix))        ! S33

                        ! calculate its eigenvalues. The Jacobi solver
                        ! requires the upper triangular matrix only.
                        call jacobi_eigenvalues(strain, eigs)

                        ! we must take the largest eigenvalue in magnitude (absolute value)
                        ggmax = max(ggmax, maxval(abs(eigs)))
                    enddo
                enddo
            enddo

            !Maximum speed:
            velmax = sqrt(maxval(velog(:, :, :, 1) ** 2   &
                               + velog(:, :, :, 2) ** 2   &
                               + velog(:, :, :, 3) ** 2))

            !Choose new time step:
            dt = min(time%alpha / (ggmax + small),  &
                     time%alpha / (bfmax + small),  &
                     cflpf / (velmax + small),      &
                     time%limit - t)

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
                dfac = vorch * dt / two
                diss = two / (one + dfac * hdis)
                !hdis = C*(K/K_max)^{2p} where K^2 = k_x^2+k_y^2, p is the order,
                !K_max is the maximum x or y wavenumber and C is a dimensionless
                !prefactor (see inversion_utils.f90 and parameters.f90 where C = prediss).
            endif

        end subroutine adapt

end module advance_mod
