! =============================================================================
!                   Write field diagnostics to NetCDF.
! =============================================================================
module field_diagnostics_netcdf
    use fields
    use inversion_utils
    use netcdf_utils
    use netcdf_writer
    use netcdf_reader
    use constants, only : one
    use parameters, only : lower, extent, nx, ny, nz
    use config, only : package_version, cf_version
    use timer, only : start_timer, stop_timer
    use options, only : write_netcdf_options
    use physics, only : write_physical_quantities
    implicit none

    private

    character(len=512) :: ncfname
    integer            :: ncid
    integer            :: t_axis_id, t_dim_id, n_writes, ke_id, en_id
#ifdef ENABLE_BUOYANCY
    integer            :: pe_id, bmax_id, bmin_id
#endif

    double precision   :: restart_time

    integer :: field_stats_io_timer

    public :: create_netcdf_field_stats_file,   &
              write_netcdf_field_stats,         &
              field_stats_io_timer


    contains

        ! Create the NetCDF field diagnostic file.
        ! @param[in] basename of the file
        ! @param[in] overwrite the file
        subroutine create_netcdf_field_stats_file(basename, overwrite)
            character(*), intent(in)  :: basename
            logical,      intent(in)  :: overwrite
            logical                   :: l_exist

            ncfname =  basename // '_field_stats.nc'

            restart_time = -one
            n_writes = 1

            call exist_netcdf_file(ncfname, l_exist)

            if (l_exist) then
                call open_netcdf_file(ncfname, NF90_NOWRITE, ncid)
                call get_num_steps(ncid, n_writes)
                call get_time(ncid, restart_time)
                call read_netcdf_field_stats_content
                call close_netcdf_file(ncid)
                n_writes = n_writes + 1
                return
            endif

            call create_netcdf_file(ncfname, overwrite, ncid)

            call write_netcdf_info(ncid=ncid,                       &
                                   ps3d_version=package_version,    &
                                   file_type='field_stats',         &
                                   cf_version=cf_version)

            call write_netcdf_box(ncid, lower, extent, (/nx, ny, nz/))

            call write_physical_quantities(ncid)

            call write_netcdf_options(ncid)

            call define_netcdf_temporal_dimension(ncid, t_dim_id, t_axis_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='ke',                                                  &
                long_name='kinetic energy',                                 &
                std_name='',                                                &
                unit='m^5/s^2',                                             &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=ke_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='en',                                                  &
                long_name='enstrophy',                                      &
                std_name='',                                                &
                unit='m^3/s^2',                                             &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=en_id)

#ifdef ENABLE_BUOYANCY
            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='pe',                                                  &
                long_name='potential energy',                               &
                std_name='',                                                &
                unit='m^5/s^2',                                             &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=pe_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='min_buoyancy',                                        &
                long_name='minimum gridded buoyancy',                       &
                std_name='',                                                &
                unit='m/s^2',                                               &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=bmin_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='max_buoyancy',                                        &
                long_name='maximum gridded buoyancy',                       &
                std_name='',                                                &
                unit='m/s^2',                                               &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=bmax_id)
#endif
            call close_definition(ncid)

        end subroutine create_netcdf_field_stats_file

        ! Pre-condition: Assumes an open file
        subroutine read_netcdf_field_stats_content

            call get_dim_id(ncid, 't', t_dim_id)

            call get_var_id(ncid, 't', t_axis_id)

            call get_var_id(ncid, 'ke', ke_id)

            call get_var_id(ncid, 'en', en_id)

#ifdef ENABLE_BUOYANCY
            call get_var_id(ncid, 'pe', pe_id)

            call get_var_id(ncid, 'min_buoyancy', bmin_id)

            call get_var_id(ncid, 'max_buoyancy', bmax_id)
#endif
        end subroutine read_netcdf_field_stats_content

        ! Write a step in the field diagnostic file.
        ! @param[in] t is the time
        ! @param[in] dt is the time step
        subroutine write_netcdf_field_stats(t)
            double precision, intent(in) :: t
            double precision             :: ke, en
#ifdef ENABLE_BUOYANCY
            double precision             :: bmin, bmax, pe
#endif

            call start_timer(field_stats_io_timer)

            if (t <= restart_time) then
                call stop_timer(field_stats_io_timer)
                return
            endif

            ke = get_kinetic_energy()
            en = get_enstrophy()

#ifdef ENABLE_BUOYANCY
            call field_combine_physical(sbuoy, buoy)
            pe = get_potential_energy()
            bmin = minval(buoy)
            bmax = maxval(buoy)
#endif

            call open_netcdf_file(ncfname, NF90_WRITE, ncid)

            ! write time
            call write_netcdf_scalar(ncid, t_axis_id, t, n_writes)

            !
            ! write diagnostics
            !
            call write_netcdf_scalar(ncid, ke_id, ke, n_writes)
            call write_netcdf_scalar(ncid, en_id, en, n_writes)
#ifdef ENABLE_BUOYANCY
            call write_netcdf_scalar(ncid, pe_id, pe, n_writes)
            call write_netcdf_scalar(ncid, bmin_id, bmin, n_writes)
            call write_netcdf_scalar(ncid, bmax_id, bmax, n_writes)
#endif

            ! increment counter
            n_writes = n_writes + 1

            call close_netcdf_file(ncid)

            call stop_timer(field_stats_io_timer)

        end subroutine write_netcdf_field_stats

end module field_diagnostics_netcdf