module utils
    use constants, only : one
    use options, only : output      &
                      , verbose     &
                      , time        &
                      , field_file  &
                      , field_step
    use netcdf_utils
    use field_netcdf
    use field_diagnostics_netcdf
    use inversion_utils, only : init_diffusion
    use inversion_mod, only : vor2vel
    use netcdf_reader, only : get_file_type, get_num_steps, get_time_at_step, get_time, get_netcdf_box
    use parameters, only : nx, ny, lower, extent, update_parameters, dx
    use fields
    use field_diagnostics
    use field_netcdf, only : read_netcdf_fields
    use physics, only : read_physical_quantities, print_physical_quantities
#ifdef ENABLE_BUOYANCY
    use physics, only : bfsq, calculate_basic_reference_state
#endif
    use mpi_layout, only : mpi_layout_init
    use mpi_utils, only : mpi_exit_on_error
    use sta3dfft, only : fftxyp2s
    implicit none

    integer :: nfw  = 0    ! number of field writes
    integer :: nsfw = 0    ! number of field diagnostics writes

    private :: nfw, nsfw


    contains

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

            if (output%write_field_stats .and. &
                (t + epsilon(zero) >= neg * dble(nsfw) * output%field_stats_freq)) then
                call write_netcdf_field_stats(t)
                nsfw = nsfw + 1
            endif

        end subroutine write_step

        subroutine setup_domain_and_parameters
            integer          :: ncid
            integer          :: ncells(3)
            double precision :: ini_time = zero

            time%initial = zero ! make sure user cannot start at arbitrary time

            ! set axis and dimension names for the NetCDF output
            call set_netcdf_dimensions((/'x', 'y', 'z', 't'/))
            call set_netcdf_axes((/'X', 'Y', 'Z', 'T'/))

            call open_netcdf_file(trim(field_file), NF90_NOWRITE, ncid)

            call get_netcdf_box(ncid, lower, extent, ncells)
            call read_physical_quantities(ncid)
            if (field_step < 1) then
                call get_time(ncid, ini_time)
            else
                call get_time_at_step(ncid, field_step, ini_time)
            endif
            time%initial = ini_time


            if (time%initial > zero) then
                nfw = int(time%initial / output%field_freq)
            endif

            call close_netcdf_file(ncid)

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
        end subroutine setup_domain_and_parameters

        subroutine setup_fields
            double precision :: bbdif, ke, ape, te, en
#ifdef ENABLE_BUOYANCY
            integer          :: iz
            double precision :: z
#endif

            call field_default

            call read_netcdf_fields(trim(field_file), field_step)

            ! decompose initial fields
#ifdef ENABLE_BUOYANCY
            bbdif = maxval(buoy) - minval(buoy)

            call calculate_basic_reference_state(nx, ny, nz, extent(3), buoy)

            ! remove basic state from buoyancy
            do iz = 0, nz
                z = lower(3) + dble(iz) * dx(3)
                bbarz(iz) = bfsq * z
                buoy(iz, :, :) = buoy(iz, :, :) - bbarz(iz)
            enddo
            call fftxyp2s(buoy, sbuoy)
#else
            bbdif = zero
#endif
            call fftxyp2s(vor(:, :, :, 1), svor(:, :, :, 1))
            call fftxyp2s(vor(:, :, :, 2), svor(:, :, :, 2))
            call fftxyp2s(vor(:, :, :, 3), svor(:, :, :, 3))

            ! calculate the initial \xi and \eta mean and save it in ini_vor_mean:
            ini_vor_mean = calc_vorticity_mean()

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

#ifndef ENABLE_SMAGORINSKY
            call init_diffusion(bbdif, te, en)
#endif

        end subroutine setup_fields

end module utils
