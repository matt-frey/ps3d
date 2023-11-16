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
    use inversion_utils, only : init_diffusion, field_decompose_physical
    use inversion_mod, only : vor2vel
    use netcdf_reader, only : get_file_type, get_num_steps, get_time_at_step, get_time, get_netcdf_box
    use parameters, only : lower, extent, update_parameters
    use fields
    use field_netcdf, only : read_netcdf_fields
    use physics, only : read_physical_quantities, print_physical_quantities, bfsq
    use mpi_layout, only : mpi_layout_init
    use mpi_utils, only : mpi_exit_on_error
    implicit none

    integer, parameter :: WRITE_VOR = 1234
    integer, parameter :: WRITE_ECOMP = 1235

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
            double precision :: ini_time

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
#if defined(ENABLE_BUOYANCY) && defined(ENABLE_PERTURBATION_MODE)
            integer          :: iz
            double precision :: z
#endif

            call field_default

            call read_netcdf_fields(trim(field_file), field_step)

            ! decompose initial fields
#ifdef ENABLE_BUOYANCY
            bbdif = maxval(buoy) - minval(buoy)

#ifdef ENABLE_PERTURBATION_MODE
            ! remove basic state from buoyancy
            do iz = 0, nz
                z = lower(3) + dble(iz) * dx(3)
                bbarz(iz) = bfsq * z
                buoy(iz, :, :) = buoy(iz, :, :) - bbarz(iz)
            enddo
#endif
            call field_decompose_physical(buoy, sbuoy)
#else
            bbdif = zero
#endif
            call field_decompose_physical(vor(:, :, :, 1), svor(:, :, :, 1))
            call field_decompose_physical(vor(:, :, :, 2), svor(:, :, :, 2))
            call field_decompose_physical(vor(:, :, :, 3), svor(:, :, :, 3))

            ! calculate the initial \xi and \eta mean and save it in ini_vor_mean:
            ini_vor_mean = calc_vorticity_mean()

            call vor2vel
            ke = get_kinetic_energy()
            ape = get_available_potential_energy()
            te = ke + ape
            en = get_enstrophy()

#ifdef ENABLE_BUOYANCY
            ! add buoyancy term to enstrophy
            en = en + get_gradb_integral()
#endif

            call init_diffusion(bbdif, te, en)

        end subroutine setup_fields

end module utils
