! =============================================================================
!                       PS3D - Pseudo-Spectral code in 3D
! =============================================================================
program ps3d
    use model_factory
    use model_manager, only : pre_run   &
                            , run       &
                            , post_run
    use mpi_environment, only : mpi_env_initialise  &
                              , mpi_env_finalise
    use mpi_utils, only : mpi_print, mpi_stop
    implicit none

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

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

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
        enddo

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
