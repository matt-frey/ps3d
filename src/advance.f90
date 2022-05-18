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
    use fields
    use jacobi, only : jacobi_eigenvalues
    implicit none

    integer :: advance_timer

    ! Number of iterations of above scheme:
    integer, parameter:: niter = 2

    double precision :: dt, dt4

    !Diagnostic quantities:
    double precision :: bfmax, vortmax, vortrms, ggmax, velmax, dfac
    double precision :: vortmp, vorl1, vorl2, vorch
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
            call adapt(t, sbuoys, velog)

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
        subroutine adapt(t, sbuoys, velog)
            double precision, intent(in) :: t
            double precision, intent(in) :: sbuoys(0:nz, 0:nx-1, 0:ny-1)
            double precision, intent(in) :: velog(0:nz, 0:ny-1, 0:nx-1, 3)
            double precision             :: svelog(0:nz, 0:nx-1, 0:ny-1, 3)
            !For defining the max strain & buoyancy frequency based time step:
            double precision, parameter  :: alpha = 0.1d0
            !Note: EPIC-2D paper recommends alpha = 0.2 for ls-rk4 method
            double precision             :: dbdxs(0:nz, 0:nx-1, 0:ny-1)     ! db/dx in spectral space
            double precision             :: dbdys(0:nz, 0:nx-1, 0:ny-1)     ! db/dy in spectral space
            double precision             :: dbdzs(0:nz, 0:nx-1, 0:ny-1)     ! db/dz in spectral space
            double precision             :: dbdx(0:nz, 0:ny-1, 0:nx-1)      ! db/dx in physical space
            double precision             :: dbdy(0:nz, 0:ny-1, 0:nx-1)      ! db/dy in physical space
            double precision             :: dbdz(0:nz, 0:ny-1, 0:nx-1)      ! db/dz in physical space
            double precision             :: strain(3, 3), eigs(3)
            double precision             :: dudx(0:nz, 0:ny-1, 0:nx-1)
            double precision             :: dudy(0:nz, 0:ny-1, 0:nx-1)
            double precision             :: dvdy(0:nz, 0:ny-1, 0:nx-1)
            double precision             :: dwdx(0:nz, 0:ny-1, 0:nx-1)
            double precision             :: dwdy(0:nz, 0:ny-1, 0:nx-1)

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

            !Maximum vorticity: (reuse dbdx)
            dbdx = vortg(:, :, :, 1) ** 2 &
                 + vortg(:, :, :, 2) ** 2 &
                 + vortg(:, :, :, 3) ** 2

            vortmax = sqrt(maxval(dbdx))

            !R.m.s. vorticity:
            vortrms = sqrt(ncelli*(f12*sum(dbdx(0, :, :)+dbdx(nz, :, :))+sum(dbdx(1:nz-1, :, :))))

            !Characteristic vorticity,  <vor^2>/<|vor|> for |vor| > vor_rms:
            vorl1 = small
            vorl2 = zero
            do ix = 0, nx-1
                do iy = 0, ny-1
                    do iz = 1, nz
                        vortmp = f12 * (vortg(iz-1, iy, ix, 1) + vortg(iz, iy, ix, 1)) &
                               + f12 * (vortg(iz-1, iy, ix, 2) + vortg(iz, iy, ix, 2)) &
                               + f12 * (vortg(iz-1, iy, ix, 3) + vortg(iz, iy, ix, 3))
                        if (abs(vortmp) .gt. vortrms) then
                            vorl1 = vorl1 + abs(vortmp)
                            vorl2 = vorl2 + vortmp ** 2
                        endif
                    enddo
                enddo
            enddo
            vorch = vorl2 / vorl1

            !Compute x derivative of velocity components:

            !
            ! velocity strain
            !
            dbdx = velog(:, :, :, 1)
            call fftxyp2s(dbdx, svelog(:, :, :, 1))
            dbdx = velog(:, :, :, 2)
            call fftxyp2s(dbdx, svelog(:, :, :, 2))
            dbdx = velog(:, :, :, 3)
            call fftxyp2s(dbdx, svelog(:, :, :, 3))


            ! du/dx
            call diffx(svelog(:, :, :, 1), dbdxs)
            call fftxys2p(dbdxs, dudx)

            ! du/dy
            call diffy(svelog(:, :, :, 1), dbdys)
            call fftxys2p(dbdys, dudy)

            ! dw/dx
            call diffx(svelog(:, :, :, 3), dbdxs)
            call fftxys2p(dbdxs, dwdx)

            ! dv/dy
            call diffy(svelog(:, :, :, 2), dbdys)
            call fftxys2p(dbdys, dvdy)

            ! dw/dy
            call diffy(svelog(:, :, :, 3), dbdys)
            call fftxys2p(dbdys, dwdy)

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
                        strain(3, 3) = -(dudx(iz, iy, ix) + dvdy(iz, iy, ix)) ! S33

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
            dt = min(alpha / (ggmax + small),  &
                     alpha / (bfmax + small),  &
                     cflpf / (velmax + small), &
                     time%limit - t)

            !Update value of dt/4:
            dt4 = dt / four

            !---------------------------------------------------------------------
            if (nnu .eq. 1) then
                !Update diffusion operator used in time stepping:
                dfac = dt / two
                do iz = 0, nz
                    diss(iz, :, :) = two / (one + dfac * hdis)
                enddo
                !hdis = nu*(k_x^2+k_y^2) where nu is the viscosity coefficient
                !(see inversion_utils.f90 and parameters.f90).
            else
                !Update hyperdiffusion operator used in time stepping:
                dfac = vorch * dt / two
                !vorch is the characteristic vorticity defined above.
                do iz = 0, nz
                    diss(iz, :, :) = two / (one + dfac * hdis)
                enddo
                !hdis = C*(K/K_max)^{2p} where K^2 = k_x^2+k_y^2, p is the order,
                !K_max is the maximum x or y wavenumber and C is a dimensionless
                !prefactor (see inversion_utils.f90 and parameters.f90 where C = prediss).
            endif

        end subroutine adapt

end module advance_mod
