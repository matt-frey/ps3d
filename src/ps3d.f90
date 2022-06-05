! =============================================================================
!                       PS3D - Pseudo-Spectral code in 3D
! =============================================================================
program ps3d
    use constants, only : zero
    use timer
    use fields
    use field_netcdf, only : field_io_timer, read_netcdf_fields
    use inversion_mod, only : vor2vel_timer, vtend_timer, vor2vel
    use inversion_utils, only : init_inversion, fftxyp2s, filt, init_hyperdiffusion
    use advance_mod, only : advance, advance_timer, WRITE_VOR, WRITE_ECOMP
    use utils, only : write_last_step, setup_output_files,       &
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
                              , output              &
                              , read_config_file    &
                              , time                &
                              , nnu                 &
                              , prediss
            double precision  :: bbdif, ke, en
            integer           :: iz

            call register_timer('ps', ps_timer)
            call register_timer('field I/O', field_io_timer)
            call register_timer('vor2vel', vor2vel_timer)
            call register_timer('vorticity tendency', vtend_timer)
            call register_timer('advance', advance_timer)

            call start_timer(ps_timer)

            ! parse the config file
            call read_config_file

            ! read domain dimensions
            call setup_domain_and_parameters(trim(field_file))

            time%initial = zero ! make sure user cannot start at arbirtrary time

            call field_default

            call init_inversion

            call read_netcdf_fields(trim(field_file))

!             ! apply Hou and Li de-aliasing filter
!             do iz = 0, nz
!                 svori(iz, :, :, 1) = filt * svori(iz, :, :, 1)
!                 svori(iz, :, :, 2) = filt * svori(iz, :, :, 2)
!                 svori(iz, :, :, 3) = filt * svori(iz, :, :, 3)
!                 sbuoy(iz, :, :)    = filt * sbuoy(iz, :, :)
!             enddo

            call vor2vel(svor, vor,  svel, vel)
            bbdif = maxval(buoyg) - minval(buoyg)
            ke = get_kinetic_energy()
            en = get_enstrophy()
            call init_hyperdiffusion(bbdif, nnu, prediss, ke, en)

            call setup_output_files

            open(WRITE_VOR, file= trim(output%basename) // '_vorticity.asc', status='replace')
            write(WRITE_VOR, '(a2, a2, a4, a4, a5, a5, a6, a6)') '# ', 't ', 'max ', 'rms ', 'char ', &
                 '<xi> ', '<eta> ', '<zeta>'

            open(WRITE_ECOMP, file= trim(output%basename) // '_ecomp.asc', status='replace')
            write(WRITE_ECOMP, '(a2, a2, a15, a9)') '# ', 't ', 'kinetic energy ', 'enstrophy'

        end subroutine


        subroutine run
            use options, only : time
            double precision :: t = zero    ! current time

            t = time%initial

            do while (t < time%limit)
                call advance(t)
            enddo

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
