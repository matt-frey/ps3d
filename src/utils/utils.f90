module utils
    use constants, only : one
    use options, only : output, verbose, time
    use field_netcdf
    use inversion_mod, only : vor2vel
    use netcdf_reader, only : get_file_type, get_num_steps, get_time_at_step, get_time, get_netcdf_box
    use parameters, only : lower, extent, update_parameters
    use fields
    use physics, only : read_physical_quantities, print_physical_quantities
    implicit none

    integer :: nfw  = 0    ! number of field writes

    private :: nfw

    contains

        ! Create NetCDF files and set the step number
        subroutine setup_output_files
            use options, only : output

            if (output%write_fields) then
                call create_netcdf_field_file(trim(output%basename), &
                                              output%overwrite)
            endif

        end subroutine setup_output_files

        ! Write last step to the NetCDF files. For the time step dt, it
        ! writes zero.
        ! @param[in] t is the time
        subroutine write_last_step(t)
            double precision,  intent(in) :: t

            ! need to be called in order to set initial time step;
            ! this is also needed for the first ls-rk4 substep
            call vor2vel

            call write_step(t, .true.)
        end subroutine write_last_step

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


        end subroutine write_step

        subroutine setup_domain_and_parameters(fname, step)
            character(*), intent(in) :: fname
            integer,      intent(in) :: step
            integer                  :: ncid, n_steps
            integer                  :: ncells(3)
            double precision         :: ini_time

            call open_netcdf_file(fname, NF90_NOWRITE, ncid)

            call get_netcdf_box(ncid, lower, extent, ncells)
            call read_physical_quantities(ncid)
            if (step < 1) then
                call get_time(ncid, ini_time)
            else
                call get_time_at_step(ncid, step, ini_time)
            endif
            time%initial = ini_time


            ! we must add +1 to nw if the initial time is not zero (a restart)
            ! to avoid that the same time is written twice
            if (time%initial > zero) then
                nfw = int(time%initial / output%field_freq) + 1
            endif

            call close_netcdf_file(ncid)

            nx = ncells(1)
            ny = ncells(2)
            nz = ncells(3)

            ! update global parameters
            call update_parameters

#ifdef ENABLE_VERBOSE
            if (verbose) then
                call print_physical_quantities
            endif
#endif
        end subroutine setup_domain_and_parameters

end module utils
