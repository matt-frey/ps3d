module advance_mod
    use options, only : time, visc_type, vor_visc, stepper
#ifdef ENABLE_VERBOSE
    use options, only : output
#endif
    use constants
    use parameters, only : nx, ny, nz, ncelli
    use inversion_mod, only : vor2vel, source
    use fields_derived, only : pressure, horizontal_divergence
#ifdef ENABLE_BUOYANCY
    use physics, only : bfsq
    use options, only : buoy_visc
#endif
    use inversion_utils
    use utils, only : write_step
    use sta2dfft, only : dst
    use fields
    use field_diagnostics
    use jacobi, only : jacobi_eigenvalues
    use mpi_environment, only : world
    use mpi_utils, only : mpi_stop
    use field_diagnostics_netcdf, only : set_netcdf_field_diagnostic        &
                                       , NC_OMAX, NC_ORMS, NC_OCHAR         &
                                       , NC_OXMEAN, NC_OYMEAN, NC_OZMEAN    &
                                       , NC_GMAX, NC_RGMAX, NC_RBFMAX       &
                                       , NC_BFMAX, NC_UMAX, NC_VMAX         &
                                       , NC_WMAX, NC_USGMAX, NC_LSGMAX
    use rolling_mean_mod, only : rolling_mean_t
    implicit none

    type, abstract :: base_stepper
        contains
            procedure(base_diffusion),  deferred :: set_diffusion
            procedure(base_setup),      deferred :: setup
            procedure(base_step),       deferred :: step
    end type

    abstract interface
        subroutine base_diffusion(self, dt, vorch, bf)
            import base_stepper
            class(base_stepper), intent(inout) :: self
            double precision,    intent(in)    :: dt
            double precision,    intent(in)    :: vorch
            double precision,    intent(in)    :: bf
        end subroutine base_diffusion

        subroutine base_setup(self)
            import base_stepper
            class(base_stepper), intent(inout) :: self
        end subroutine base_setup

        subroutine base_step(self, t, dt)
            import base_stepper
            class(base_stepper), intent(inout) :: self
            double precision,    intent(inout) :: t
            double precision,    intent(in)    :: dt
        end subroutine base_step

    end interface

    type(rolling_mean_t) :: rollmean
    type(rolling_mean_t) :: buoy_rollmean

    integer :: advance_timer

    !Diagnostic quantities:
    double precision :: bfmax, vortmax, vortrms, ggmax, velmax
    double precision :: usggmax, lsggmax, rmv
#ifdef ENABLE_BUOYANCY
    double precision :: rmb
#endif
    double precision :: vorch
    integer          :: ix, iy, iz

    contains

        subroutine advance(bstep, t)
            class(base_stepper), intent(inout) :: bstep
            double precision,    intent(inout) :: t
            double precision                   :: dt

            !-------------------------------------------------------------------
            !Invert vorticity for velocity at current time level, say t=t^n:
            !Also, returns vorticity in physical space for use everywhere
            call vor2vel

            !Adapt the time step
            call adapt(bstep, t, dt)

            !Write fields
            call write_step(t)

            !Calculate the source terms (sbuoys, svorts) for buoyancy (sbuoy) and
            !vorticity in spectral space:
            call source

            !Advance time:
            if (world%rank == world%root) then
                print *, "At time", t, "and time step", dt
            endif

            call bstep%step(t, dt)

        end subroutine advance

        !=======================================================================

        ! Adapts the time step and computes various diagnostics
        subroutine adapt(bstep, t, dt)
            class(base_stepper), intent(inout) :: bstep
            double precision,    intent(in)    :: t
            double precision,    intent(inout) :: dt
            double precision                :: xs(0:nz, box%lo(2):box%hi(2), &
                                                        box%lo(1):box%hi(1)) ! derivatives in x in spectral space
            double precision                :: ys(0:nz, box%lo(2):box%hi(2), &
                                                        box%lo(1):box%hi(1)) ! derivatives in y in spectral space
            double precision                :: xp(0:nz, box%lo(2):box%hi(2), &
                                                        box%lo(1):box%hi(1)) ! derivatives in x in physical space
#ifdef ENABLE_BUOYANCY
            double precision                :: yp(0:nz, box%lo(2):box%hi(2), &
                                                        box%lo(1):box%hi(1)) ! derivatives in y physical space
            double precision                :: zp(0:nz, box%lo(2):box%hi(2), &
                                                        box%lo(1):box%hi(1)) ! derivatives in z physical space
#endif
            double precision                :: strain(3, 3), eigs(3)
            double precision                :: dudx(0:nz, box%lo(2):box%hi(2), &
                                                          box%lo(1):box%hi(1)) ! du/dx in physical space
            double precision                :: dudy(0:nz, box%lo(2):box%hi(2), &
                                                          box%lo(1):box%hi(1)) ! du/dy in physical space
            double precision                :: dvdy(0:nz, box%lo(2):box%hi(2), &
                                                          box%lo(1):box%hi(1)) ! dv/dy in physical space
            double precision                :: dwdx(0:nz, box%lo(2):box%hi(2), &
                                                          box%lo(1):box%hi(1)) ! dw/dx in physical space
            double precision                :: dwdy(0:nz, box%lo(2):box%hi(2), &
                                                          box%lo(1):box%hi(1)) ! dw/dy in physical space
            double precision                :: vormean(3)
            double precision                :: buf(7)
            double precision                :: umax, vmax, wmax, dtcfl
            double precision                :: vval, bval, lmax
#ifdef ENABLE_VERBOSE
            logical                         :: l_exist = .false.
            character(512)                  :: fname
#endif

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
            vortmax = dsqrt(get_abs_max(xp, l_allreduce=.false.))

            !R.m.s. vorticity: (note that xp is already squared, hence, we only need get_mean)
            vortrms = dsqrt(get_mean(xp, l_allreduce=.true.))

            !Characteristic vorticity,  <vor^2>/<|vor|> for |vor| > vor_rms:
            vorch = get_char_vorticity(vortrms, l_allreduce=.true.)

            vormean = get_mean_vorticity(l_allreduce=.false.)

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
            usggmax = zero
            lsggmax = zero
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
                        lmax = maxval(abs(eigs))
                        ggmax = max(ggmax, lmax)

                        if (iz == nz) then
                            usggmax = max(usggmax, lmax)
                        endif

                        if (iz == 0) then
                            lsggmax = max(lsggmax, lmax)
                        endif
                    enddo
                enddo
            enddo

            ! Pressure calculation:
            call pressure(dudx, dudy, dvdy, dwdx, dwdy)

            ! Horizontal divergence calculation:
            call horizontal_divergence

            !Maximum speed:
            !$omp parallel workshare
            umax = maxval(vel(:, :, :, 1))
            vmax = maxval(vel(:, :, :, 2))
            wmax = maxval(vel(:, :, :, 3))
            !$omp end parallel workshare

            buf(1) = bfmax
            buf(2) = ggmax
            buf(3) = umax
            buf(4) = vmax
            buf(5) = wmax
            buf(6) = usggmax
            buf(7) = lsggmax

            call MPI_Allreduce(MPI_IN_PLACE,            &
                               buf(1:7),                &
                               7,                       &
                               MPI_DOUBLE_PRECISION,    &
                               MPI_MAX,                 &
                               world%comm,              &
                               world%err)

            bfmax = buf(1)
            ggmax = buf(2)
            umax = buf(3)
            vmax = buf(4)
            wmax = buf(5)
            usggmax = buf(6)
            lsggmax = buf(7)

            call set_netcdf_field_diagnostic(bfmax, NC_BFMAX)
            call set_netcdf_field_diagnostic(ggmax, NC_GMAX)
            call set_netcdf_field_diagnostic(umax, NC_UMAX)
            call set_netcdf_field_diagnostic(vmax, NC_VMAX)
            call set_netcdf_field_diagnostic(wmax, NC_WMAX)
            call set_netcdf_field_diagnostic(usggmax, NC_USGMAX)
            call set_netcdf_field_diagnostic(lsggmax, NC_LSGMAX)


            ! CFL time step constraint:
            dtcfl = cflmax * min(dx(1) / (umax + small), &
                                 dx(2) / (vmax + small), &
                                 dx(3) / (wmax + small))

            !Choose new time step:
            dt = min(time%alpha / (ggmax + small),  &
                     time%alpha / (bfmax + small),  &
                     dtcfl,                         &
                     time%limit - t)

#ifdef ENABLE_VERBOSE
            if (world%rank == world%root) then
                fname = trim(output%basename) // '_alpha_time_step.asc'
                inquire(file=trim(fname), exist=l_exist)
                if ((t /= zero) .and. l_exist) then
                    open(unit=1235, file=trim(fname), status='old', position='append')
                else
                    open(unit=1235, file=trim(fname), status='replace')
                    write(1235, *) '  # time (s)                \alpha_s/\gamma_{max}     ' // &
                                   '\alpha_b/N_{max}          CFL'
                endif

                write(1235, *) t, time%alpha / (ggmax + small), time%alpha / (bfmax + small), &
                               dtcfl

                close(1235)
            endif
#endif


            vval = zero
            bval = zero

#ifdef ENABLE_BUOYANCY
            call buoy_rollmean%alloc(buoy_visc%roll_mean_win_size)
            rmb = buoy_rollmean%get_next(bfmax)
            call set_netcdf_field_diagnostic(rmb, NC_RBFMAX)
#endif

            call rollmean%alloc(vor_visc%roll_mean_win_size)
            rmv = rollmean%get_next(ggmax)
            call set_netcdf_field_diagnostic(rmv, NC_RGMAX)

            vval = get_diffusion_pre_factor(vor_visc)

#ifdef ENABLE_BUOYANCY

            bval = get_diffusion_pre_factor(buoy_visc)
#endif

            call bstep%set_diffusion(dt, vval, bval)

        end subroutine adapt

        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        function get_diffusion_pre_factor(visc) result(val)
            type(visc_type),  intent(in) :: visc
            double precision             :: val

            select case (visc%pretype)
                case ('constant')
                    ! diffusion is only controlled by the viscosity: diss = -nu*(k^2+l^2)*dt/2
                    val = one
                case ('vorch')
                    val = vorch
                case ('bfmax')
                    val = bfmax
                case ('roll-mean-max-strain')
                    val = rmv
#ifdef ENABLE_BUOYANCY
                case ('roll-mean-bfmax')
                    val = rmb
#endif
                case ('max-strain')
                    val = ggmax
                case ('us-max-strain')
                    val = usggmax
                case default
                    call mpi_stop(&
                        "We only support 'constant', 'vorch', 'bfmax', " // &
                        "'roll-mean-max-strain', 'roll-mean-bfmax', " // &
                        "'max-strain' and us-max-strain")
            end select

        end function get_diffusion_pre_factor

end module advance_mod
