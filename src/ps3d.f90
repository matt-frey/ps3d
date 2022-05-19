! =============================================================================
!                       PS3D - Pseudo-Spectral code in 3D
! =============================================================================
program ps3d
    use constants, only : zero
    use timer
    use fields
    use field_netcdf, only : field_io_timer, read_netcdf_fields
    use inversion_mod, only : vor2vel_timer, vtend_timer
    use inversion_utils, only : init_inversion
    use advance_mod, only : advance, advance_timer
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
            double precision  :: bbdif

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

            call init_inversion(bbdif, nnu, prediss)

            call field_default

            call read_netcdf_fields(trim(field_file))

            call setup_output_files

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
