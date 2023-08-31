! =============================================================================
!                       PS3D - Pseudo-Spectral code in 3D
! =============================================================================
program ps3d
    use constants, only : zero
    use timer
    use fields
    use field_netcdf, only : field_io_timer, read_netcdf_fields
    use field_diagnostics_netcdf, only : field_stats_io_timer
    use inversion_mod, only : vor2vel_timer &
                            , vtend_timer   &
                            , vor2vel       &
                            , pres_timer
    use inversion_utils, only : init_inversion          &
                              , init_diffusion          &
                              , field_decompose_physical
    use advance_mod, only : advance             &
                          , calc_vorticity_mean &
                          , advance_timer       &
                          , WRITE_VOR           &
                          , WRITE_ECOMP
    use utils, only : write_last_step, setup_output_files,  &
                      setup_domain_and_parameters
    implicit none

    integer          :: ps_timer

    ! Read command line (verbose, filename, etc.)
    call parse_command_line

    ! Create the model
    call pre_run

    ! Run the model
    call run

    ! Deallocate memory
    call post_run

    contains

        subroutine pre_run
            use options, only : field_file          &
                              , field_step          &
                              , output              &
                              , read_config_file    &
                              , time
            double precision :: bbdif, ke, ape, te, en
#if defined(ENABLE_BUOYANCY) && defined(ENABLE_PERTURBATON_MODE)
            integer          :: iz
            double precision :: z
#endif

            call register_timer('ps', ps_timer)
            call register_timer('field I/O', field_io_timer)
            call register_timer('field diagnostics I/O', field_stats_io_timer)
            call register_timer('vor2vel', vor2vel_timer)
            call register_timer('vorticity tendency', vtend_timer)
            call register_timer('advance', advance_timer)
            call register_timer('pressure calculation', pres_timer)

            call start_timer(ps_timer)

            ! parse the config file
            call read_config_file

            time%initial = zero ! make sure user cannot start at arbitrary time

            ! read domain dimensions
            call setup_domain_and_parameters(trim(field_file), field_step)

            call init_inversion

            call field_default

            call read_netcdf_fields(trim(field_file), field_step)

            ! decompose initial fields
#ifdef ENABLE_BUOYANCY

#ifdef ENABLE_PERTURBATON_MODE
            ! N^2 = (db/dz)^2
            bfsq = sum(buoy(nz, :, :) - buoy(0, :, :)) / (dble(nx * ny) * extent(3))
            print *, "Calculated squared buoyancy frequency:", bfsq

            ! remove basic state from buoyancy
            do iz = 0, nz
                z = lower(3) + dble(iz) * dx(3)
                buoy(iz, :, :) = buoy(iz, :, :) - bfsq * z
            enddo
#endif
            call field_decompose_physical(buoy, sbuoy)
#endif
            call field_decompose_physical(vor(:, :, :, 1), svor(:, :, :, 1))
            call field_decompose_physical(vor(:, :, :, 2), svor(:, :, :, 2))
            call field_decompose_physical(vor(:, :, :, 3), svor(:, :, :, 3))

            ! calculate the initial \xi and \eta mean and save it in ini_vor_mean:
            ini_vor_mean = calc_vorticity_mean()

            call vor2vel
#ifdef ENABLE_BUOYANCY
            bbdif = maxval(buoy) - minval(buoy)
#else
            bbdif = zero
#endif
            ke = get_kinetic_energy()
            ape = get_available_potential_energy()
            te = ke + ape
            en = get_enstrophy()

#ifdef ENABLE_BUOYANCY
            ! add buoyancy term to enstrophy
            en = en + get_gradb_integral()
#endif

            call init_diffusion(bbdif, te, en)


            call setup_output_files

            open(WRITE_VOR, file= trim(output%basename) // '_vorticity.asc', status='replace')
            write(WRITE_VOR, '(a2, a2, a4, a4, a5, a5, a6, a6)') '# ', 't ', 'max ', 'rms ', 'char ', &
                 '<xi> ', '<eta> ', '<zeta>'

            open(WRITE_ECOMP, file= trim(output%basename) // '_ecomp.asc', status='replace')
#ifdef ENABLE_BUOYANCY
            write(WRITE_ECOMP, '(a2, a2, a15, a17, a9)') '# ', 't ', 'kinetic energy ', &
                                                         'potential energy ', 'enstrophy'
#else
            write(WRITE_ECOMP, '(a2, a2, a15, a9)') '# ', 't ', 'kinetic energy ', 'enstrophy'
#endif
        end subroutine


        subroutine run
            use options, only : time
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

        subroutine post_run
            use options, only : output

            close(WRITE_VOR)
            close(WRITE_ECOMP)

            call stop_timer(ps_timer)
            call write_time_to_csv(output%basename)
            call print_timer
        end subroutine


    ! Get the file name provided via the command line
    subroutine parse_command_line
        use options, only : filename
#ifdef ENABLE_VERBOSE
        use options, only : verbose
#endif
        integer                          :: i
        character(len=512)               :: arg

        filename = ''
        i = 0
        do
            call get_command_argument(i, arg)
            if (len_trim(arg) == 0) then
                exit
            endif

            if (arg == '--config') then
                i = i + 1
                call get_command_argument(i, arg)
                filename = trim(arg)
            else if (arg == '--help') then
                print *, 'Run code with "./ps3d --config [config file]"'
                stop
#ifdef ENABLE_VERBOSE
            else if (arg == '--verbose') then
                verbose = .true.
#endif
            endif
            i = i+1
        end do

        if (filename == '') then
            print *, 'No configuration file provided. Run code with "./ps3d --config [config file]"'
            stop
        endif

#ifdef ENABLE_VERBOSE
        ! This is the main application of PS
        if (verbose) then
            print *, 'Running PS3D with "', trim(filename), '"'
        endif
#endif
    end subroutine parse_command_line
end program ps3d
