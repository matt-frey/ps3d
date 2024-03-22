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
    use options, only : time, viscosity, stepper
    use constants
    use parameters, only : nx, ny, nz, glmin, cflpf, ncelli
    use inversion_mod, only : vor2vel, pressure
#ifdef ENABLE_BUOYANCY
#ifdef ENABLE_PERTURBATION_MODE
    use physics, only : bfsq
#endif
#endif
    use inversion_utils
    use utils, only : write_step
    use sta2dfft, only : dst
    use fields
    use field_diagnostics
    use jacobi, only : jacobi_eigenvalues
    use mpi_environment, only : world
    use field_diagnostics_netcdf, only : set_netcdf_field_diagnostic    &
                                       , NC_OMAX, NC_ORMS, NC_OCHAR     &
                                       , NC_OXMEAN, NC_OYMEAN, NC_OZMEAN
    use cn2, only : cn2_step
    use ls_rk, only : ls_rk_step
    use mpi_utils, only : mpi_stop
    implicit none

    integer :: advance_timer

    double precision :: dt

    !Diagnostic quantities:
    double precision :: bfmax, vortmax, vortrms, ggmax, velmax, dfac
    double precision :: vorch
    integer          :: ix, iy, iz

    contains

        subroutine advance(t)
            double precision, intent(inout) :: t

            !-------------------------------------------------------------------
            !Invert vorticity for velocity at current time level, say t=t^n:
            !Also, returns vorticity in physical space for use everywhere
            call vor2vel

            !Adapt the time step
            call adapt(t)

            !Advance time:
            if (world%rank == world%root) then
                print *, "At time", t, "and time step", dt
            endif

            !Write fields
            call write_step(t)

            select case (stepper)
                case ('RK4')
                    call ls_rk_step(t, dt)
                case ('RK3')
                    call ls_rk_step(t, dt)
                case ('CN2')
                    call cn2_step(t, dt)
                case default
                    call mpi_stop('Only RK or CN time integrators.')
            end select

        end subroutine advance

        !=======================================================================

        ! Adapts the time step and computes various diagnostics
        subroutine adapt(t)
            double precision, intent(in) :: t
            double precision             :: xs(0:nz, box%lo(2):box%hi(2), &
                                                     box%lo(1):box%hi(1)) ! derivatives in x in spectral space
            double precision             :: ys(0:nz, box%lo(2):box%hi(2), &
                                                     box%lo(1):box%hi(1)) ! derivatives in y in spectral space
            double precision             :: xp(0:nz, box%lo(2):box%hi(2), &
                                                     box%lo(1):box%hi(1)) ! derivatives in x in physical space
#ifdef ENABLE_BUOYANCY
            double precision             :: yp(0:nz, box%lo(2):box%hi(2), &
                                                     box%lo(1):box%hi(1)) ! derivatives in y physical space
            double precision             :: zp(0:nz, box%lo(2):box%hi(2), &
                                                     box%lo(1):box%hi(1)) ! derivatives in z physical space
#endif
            double precision             :: strain(3, 3), eigs(3)
            double precision             :: dudx(0:nz, box%lo(2):box%hi(2), &
                                                       box%lo(1):box%hi(1)) ! du/dx in physical space
            double precision             :: dudy(0:nz, box%lo(2):box%hi(2), &
                                                       box%lo(1):box%hi(1)) ! du/dy in physical space
            double precision             :: dvdy(0:nz, box%lo(2):box%hi(2), &
                                                       box%lo(1):box%hi(1)) ! dv/dy in physical space
            double precision             :: dwdx(0:nz, box%lo(2):box%hi(2), &
                                                       box%lo(1):box%hi(1)) ! dw/dx in physical space
            double precision             :: dwdy(0:nz, box%lo(2):box%hi(2), &
                                                       box%lo(1):box%hi(1)) ! dw/dy in physical space
            double precision             :: vormean(3)
            double precision             :: buf(3)

            bfmax = zero

#ifdef ENABLE_BUOYANCY
            !Obtain x, y & z derivatives of buoyancy -> xs, ys, zs
            !Obtain gradient of buoyancy in physical space -> xp, yp, zp
            call field_combine_semi_spectral(sbuoy)
            call diffx(sbuoy, xs)
            call fftxys2p(xs, xp)

            call diffy(sbuoy, ys)
            call fftxys2p(ys, yp)

            call central_diffz(sbuoy, xs)
            call fftxys2p(xs, zp)
            call field_decompose_semi_spectral(sbuoy)


            !Compute (db/dx)^2 + (db/dy)^2 + (db/dz)^2 -> xp in physical space:
            !$omp parallel workshare
            xp = xp ** 2 + yp ** 2 + (zp + bfsq) ** 2
            !$omp end parallel workshare

            !Maximum buoyancy frequency:
            bfmax = dsqrt(dsqrt(maxval(xp)))
#endif


            !Compute enstrophy: (reuse xp)
            !$omp parallel workshare
            xp = vor(:, :, :, 1) ** 2 + vor(:, :, :, 2) ** 2 + vor(:, :, :, 3) ** 2
            !$omp end parallel workshare

            !Maximum vorticity magnitude:
            vortmax = dsqrt(get_abs_max(xp))

            !R.m.s. vorticity: (note that xp is already squared, hence, we only need get_mean)
            vortrms = dsqrt(get_mean(xp))

            !Characteristic vorticity,  <vor^2>/<|vor|> for |vor| > vor_rms:
            vorch = get_char_vorticity(vortrms)

            vormean = get_mean_vorticity()

            ! update diagnostics in netCDF data structure (avoids the re-evaluation)
            call set_netcdf_field_diagnostic(vortmax, NC_OMAX)
            call set_netcdf_field_diagnostic(vortrms, NC_ORMS)
            call set_netcdf_field_diagnostic(vorch, NC_OCHAR)
            call set_netcdf_field_diagnostic(vormean(1), NC_OXMEAN)
            call set_netcdf_field_diagnostic(vormean(2), NC_OYMEAN)
            call set_netcdf_field_diagnostic(vormean(3), NC_OZMEAN)

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
            do ix = box%lo(1), box%hi(1)
                do iy = box%lo(2), box%hi(2)
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

            ! Pressure calculation:
            call pressure(dudx, dudy, dvdy, dwdx, dwdy)

            !Maximum speed:
            !$omp parallel workshare
            velmax = maxval(vel(:, :, :, 1) ** 2   &
                          + vel(:, :, :, 2) ** 2   &
                          + vel(:, :, :, 3) ** 2)
            !$omp end parallel workshare

            buf(1) = bfmax
            buf(2) = ggmax
            buf(3) = velmax

            call MPI_Allreduce(MPI_IN_PLACE,            &
                               buf(1:3),                &
                               3,                       &
                               MPI_DOUBLE_PRECISION,    &
                               MPI_MAX,                 &
                               world%comm,              &
                               world%err)

            bfmax = buf(1)
            ggmax = buf(2)
            velmax = buf(3)

            velmax = dsqrt(velmax)

            !Choose new time step:
            dt = min(time%alpha / (ggmax + small),  &
                     time%alpha / (bfmax + small),  &
                     cflpf / (velmax + small),      &
                     time%limit - t)

            !---------------------------------------------------------------------
            if (viscosity%nnu .eq. 1) then
                !Update diffusion operator used in time stepping:
                dfac = dt
                !$omp parallel workshare
                diss = one / (one + dfac * hdis)
                !$omp end parallel workshare
                !(see inversion_utils.f90)
            else
                !Update hyperdiffusion operator used in time stepping:
                dfac = vorch * dt
                !$omp parallel workshare
                diss = one / (one + dfac * hdis)
                !$omp end parallel workshare
                !(see inversion_utils.f90)
             endif

        end subroutine adapt

end module advance_mod
