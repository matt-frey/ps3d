! =============================================================================
!                       PS3D - Pseudo-Spectral code in 3D
! =============================================================================
program ps3d
    use constants, only : zero
    use mpi_timer
    use fields
    use field_netcdf, only : field_io_timer
    use field_diagnostics_netcdf, only : field_stats_io_timer
    use inversion_mod, only : vor2vel_timer &
                            , vtend_timer   &
                            , vor2vel       &
                            , pres_timer
    use inversion_utils, only : init_inversion, finalise_inversion
    use advance_mod, only : advance             &
                          , calc_vorticity_mean &
                          , advance_timer
    use utils, only : write_last_step, setup_output_files,   &
                      setup_domain_and_parameters            &
                    , setup_fields
    use mpi_environment, only : mpi_env_initialise, mpi_env_finalise
    use mpi_utils, only : mpi_print, mpi_stop
    use ls_rk, only : ls_rk_setup
    implicit none

    integer          :: ps_timer

    call mpi_env_initialise

    ! Read command line (verbose, filename, etc.)
    call parse_command_line

    ! Create the model
    call pre_run

    ! Run the model
    call run

    ! Deallocate memory
    call post_run

    call mpi_env_finalise

    contains

        subroutine pre_run
            use options, only : output              &
                              , read_config_file    &
                              , stepper
            integer :: rk_order

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

            ! read domain dimensions
            call setup_domain_and_parameters

            call init_inversion

            call setup_fields

            call setup_output_files

            select case (stepper)
                case ('RK4')
                    rk_order = 4
                case ('RK3')
                    rk_order = 3
                case default
                    rk_order = 4
            end select

            call ls_rk_setup(rk_order)

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

            call finalise_inversion

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
                call mpi_stop('Run code with "./ps3d --config [config file]"')
#ifdef ENABLE_VERBOSE
            else if (arg == '--verbose') then
                verbose = .true.
#endif
            endif
            i = i+1
        end do

        if (filename == '') then
            call mpi_stop(&
                'No configuration file provided. Run code with "./ps3d --config [config file]"')
        endif

#ifdef ENABLE_VERBOSE
        ! This is the main application of PS
        if (verbose) then
            call mpi_print('Running PS3D with "' // trim(filename) // '"')
        endif
#endif
    end subroutine parse_command_line
end program ps3d
