module model_manager
    use options, only : time                &
                      , visc_type           &
                      , vor_visc            &
                      , time_stepper        &
                      , read_config_file    &
                      , output              &
                      , verbose             &
                      , field_file          &
                      , field_step
    use constants
    use parameters, only : nx, ny, nz, ncelli   &
                         , update_parameters    &
                         , grid_type
    use physics, only : read_physical_quantities    &
                      , print_physical_quantities
#ifdef ENABLE_BUOYANCY
    use physics, only : bfsq, calculate_basic_reference_state
    use options, only : buoy_visc
#endif
    use field_diagnostics
    use jacobi, only : jacobi_eigenvalues
    use mpi_environment, only : world
    use mpi_utils, only : mpi_stop  &
                        , mpi_exit_on_error
    use mpi_layout, only : mpi_layout_init
    use field_diagnostics_netcdf, only : set_netcdf_field_diagnostic        &
                                       , NC_OMAX, NC_ORMS, NC_OCHAR         &
                                       , NC_OXMEAN, NC_OYMEAN, NC_OZMEAN    &
                                       , NC_GMAX, NC_RGMAX, NC_RBFMAX       &
                                       , NC_BFMAX, NC_UMAX, NC_VMAX         &
                                       , NC_WMAX, NC_USGMAX, NC_LSGMAX
    use rolling_mean_mod, only : rolling_mean_t
    use model_factory, only : create_model, layout
    use constants, only : zero
    use mpi_timer
    use fields
    use field_netcdf, only : field_io_timer             &
                           , read_netcdf_fields         &
                           , write_netcdf_fields        &
                           , create_netcdf_field_file
    use field_diagnostics_netcdf, only : field_stats_io_timer           &
                                       , write_netcdf_field_stats       &
                                       , create_netcdf_field_stats_file
    use inversion_mod, only : vor2vel_timer &
                            , vtend_timer   &
                            , vor2vel       &
                            , source
    use fields_derived, only : pres_timer               &
                             , delta_timer              &
                             , field_derived_default    &
                             , pressure                 &
                             , horizontal_divergence
    use diffusion, only : init_diffusion
#ifdef ENABLE_BALANCE
    use field_balance, only : initialise_balance, finalise_balance
#endif
    use stepper_factory, only : stepper_t, create_stepper
    use netcdf_utils
    use netcdf_reader, only : get_file_type         &
                            , get_num_steps         &
                            , get_time_at_step      &
                            , get_time              &
                            , get_netcdf_box        &
                            , read_netcdf_attribute
    implicit none

    private

    integer              :: nfw  = 0    ! number of field writes
    integer              :: nsfw = 0    ! number of field diagnostics writes
    integer              :: ps_timer
    integer              :: advance_timer
    type(rolling_mean_t) :: rollmean
#ifdef ENABLE_BUOYANCY
    type(rolling_mean_t) :: buoy_rollmean
#endif
    class(stepper_t), allocatable :: stepper

    !Diagnostic quantities:
    double precision :: bfmax, vortmax, vortrms, ggmax
    double precision :: usggmax, lsggmax, rmv
#ifdef ENABLE_BUOYANCY
    double precision :: rmb
#endif
    double precision :: vorch
    integer          :: ix, iy, iz

    public :: pre_run, run, post_run

contains

    subroutine pre_run
        call register_timer('ps', ps_timer)
        call register_timer('field I/O', field_io_timer)
        call register_timer('field diagnostics I/O', field_stats_io_timer)
        call register_timer('vor2vel', vor2vel_timer)
        call register_timer('vorticity tendency', vtend_timer)
        call register_timer('advance', advance_timer)
        call register_timer('pressure calculation', pres_timer)
        call register_timer('horizontal divergence calculation', delta_timer)

        call start_timer(ps_timer)

        ! parse the config file
        call read_config_file

        ! read domain dimensions
        call setup_domain_and_parameters(trim(field_file), field_step)

#ifdef ENABLE_BALANCE
        if (output%l_balanced) then
            call initialise_balance
        endif
#endif

        call setup_fields(trim(field_file), field_step)

        call setup_output_files

        stepper = create_stepper(time_stepper)

    end subroutine

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine run
        double precision :: t = zero    ! current time

        t = time%initial

        call start_timer(advance_timer)
        do while (t < time%limit)
            call advance(t)
        enddo
        call stop_timer(advance_timer)

        ! write final step (we only write if we really advanced in time)
        if (t > time%initial) then
            call write_last_step(t)
        endif

    end subroutine run

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine post_run
#ifdef ENABLE_BALANCE
        if (output%l_balanced) then
            call finalise_balance
        endif
#endif

        call stop_timer(ps_timer)
        call write_time_to_csv(output%basename)
        call print_timer

    end subroutine

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine advance(t)
        double precision, intent(inout) :: t
        double precision                :: dt

        !-------------------------------------------------------------------
        !Invert vorticity for velocity at current time level, say t=t^n:
        !Also, returns vorticity in physical space for use everywhere
        call vor2vel

        !Adapt the time step
        call adapt(t, dt)

        !Write fields
        call write_step(t)

        !Calculate the source terms (sbuoys, svorts) for buoyancy (sbuoy) and
        !vorticity in spectral space:
        call source

        !Advance time:
        if (world%rank == world%root) then
            print *, "At time", t, "and time step", dt
        endif

        call stepper%step(t, dt)

    end subroutine advance

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Adapts the time step and computes various diagnostics
    subroutine adapt(t, dt)
        double precision, intent(in)    :: t
        double precision, intent(inout) :: dt
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
        call combine_semi_spectral(sbuoy)
        call diffx(sbuoy, xs)
        call fftxys2p(xs, xp)

        call diffy(sbuoy, ys)
        call fftxys2p(ys, yp)

        call ops%diffz(sbuoy, xs)
        call fftxys2p(xs, zp)
        call decompose_semi_spectral(sbuoy)

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
        vortmax = dsqrt(ops%get_absmax(xp, l_allreduce=.false.))

        !R.m.s. vorticity: (note that xp is already squared, hence, we only need get_mean)
        vortrms = dsqrt(ops%get_mean(xp, l_allreduce=.true.))

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

        call stepper%set_diffusion(dt, vval, bval)

    end subroutine adapt

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

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

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Create NetCDF files and set the step number
    subroutine setup_output_files

        if (output%write_fields) then
            call create_netcdf_field_file(trim(output%basename), &
                                          output%overwrite)
        endif

        if (output%write_field_stats) then
            call create_netcdf_field_stats_file(trim(output%basename),   &
                                                output%overwrite)
        endif

    end subroutine setup_output_files

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Write last step to the NetCDF files. For the time step dt, it
    ! writes zero.
    ! @param[in] t is the time
    subroutine write_last_step(t)
        double precision,  intent(in) :: t

        ! need to be called in order to set initial time step;
        call vor2vel

        call write_step(t, .true.)
    end subroutine write_last_step

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Write step to the NetCDF files.
    ! @param[in] t is the time
    ! @param[in] l_force a logical to force a write (optional)
    subroutine write_step(t, l_force)
        double precision,  intent(in) :: t
        logical, optional, intent(in) :: l_force
        double precision              :: neg = one

        if (present(l_force)) then
            if (l_force) then
                neg = -one
            endif
        endif

        ! make sure we always write initial setup
        if (output%write_fields .and. &
            (t + epsilon(zero) >= neg * dble(nfw) * output%field_freq)) then
            call write_netcdf_fields(t)
            nfw = nfw + 1
        endif

        if (output%write_field_stats .and. &
            (t + epsilon(zero) >= neg * dble(nsfw) * output%field_stats_freq)) then
            call write_netcdf_field_stats(t)
            nsfw = nsfw + 1
        endif

    end subroutine write_step

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine setup_domain_and_parameters(fname, step)
        character(*), intent(in) :: fname
        integer,      intent(in) :: step
        integer                  :: ncid, gid
        integer                  :: ncells(3)
        double precision         :: ini_time = zero

        time%initial = zero ! make sure user cannot start at arbitrary time

        ! set axis and dimension names for the NetCDF output
        call set_netcdf_dimensions((/'x', 'y', 'z', 't'/))
        call set_netcdf_axes((/'X', 'Y', 'Z', 'T'/))

        call open_netcdf_file(fname, NF90_NOWRITE, ncid)

        call read_physical_quantities(ncid)

        if (step < 1) then
            call get_time(ncid, ini_time)
        else
            call get_time_at_step(ncid, step, ini_time)
        endif
        time%initial = ini_time

        if (time%initial > zero) then
            nfw = int(time%initial / output%field_freq)

        endif

        ncerr = nf90_inq_ncid(ncid, 'parameters', gid)

        call check_netcdf_error("No group 'parameters'.")

        call get_netcdf_box(gid, lower, extent, ncells)

        nx = ncells(1)
        ny = ncells(2)
        nz = ncells(3)

        call mpi_layout_init(lower, extent, nx, ny, nz)

        ! update global parameters
        call update_parameters

#ifdef ENABLE_VERBOSE
        if (verbose) then
            call print_physical_quantities
        endif
#endif


        call read_netcdf_attribute(gid, 'grid_type', grid_type)

        call create_model(grid_type)

        call close_netcdf_file(ncid)

    end subroutine setup_domain_and_parameters

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine setup_fields(fname, step)
        character(*), intent(in) :: fname
        integer,      intent(in) :: step
        double precision         :: ke, ape, te, en
        integer                  :: nc
#ifdef ENABLE_BUOYANCY
        integer                  :: iz
        double precision         :: z
#endif

        call field_default
        call field_derived_default

        call read_netcdf_fields(fname, step)

        ! decompose initial fields
#ifdef ENABLE_BUOYANCY
        call calculate_basic_reference_state(nx, ny, nz, extent(3), buoy)

        ! remove basic state from buoyancy
        do iz = 0, nz
            z = lower(3) + dble(iz) * dx(3)
            bbarz(iz) = bfsq * z
            buoy(iz, :, :) = buoy(iz, :, :) - bbarz(iz)
        enddo
        call layout%decompose_physical(buoy, sbuoy)
#endif
        call layout%decompose_physical(vor(:, :, :, 1), svor(:, :, :, 1))
        call layout%decompose_physical(vor(:, :, :, 2), svor(:, :, :, 2))
        call layout%decompose_physical(vor(:, :, :, 3), svor(:, :, :, 3))

        ! calculate the initial \xi and \eta mean and save it in ini_vor_mean:
        do nc = 1, 2
            ini_vor_mean(nc) = ops%calc_decomposed_mean(svor(:, :, :, nc))
        enddo

        call vor2vel
        ke = get_kinetic_energy(vel, l_global=.true., l_allreduce=.true.)
#ifdef ENABLE_BUOYANCY
        ape = get_available_potential_energy(buoy, l_global=.true., l_allreduce=.true.)
#else
        ape = zero
#endif
        te = ke + ape
        en = get_enstrophy(l_global=.true., l_allreduce=.true.)

#ifdef ENABLE_BUOYANCY
        ! add buoyancy term to enstrophy
        en = en + get_gradb_integral(l_global=.true., l_allreduce=.true.)
#endif

        call init_diffusion(te, en)

    end subroutine setup_fields

end module model_manager
