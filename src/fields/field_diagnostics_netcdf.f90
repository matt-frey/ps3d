! =============================================================================
!                   Write field diagnostics to NetCDF.
! =============================================================================
module field_diagnostics_netcdf
    use fields
    use inversion_utils
    use netcdf_utils
    use netcdf_writer
    use netcdf_reader
    use constants, only : zero, one
    use parameters, only : lower, extent, nx, ny, nz
    use config, only : package_version, cf_version
    use mpi_timer, only : start_timer, stop_timer
    use options, only : write_netcdf_options
    use physics, only : write_physical_quantities
    use mpi_collectives, only : mpi_blocking_reduce
    use field_diagnostics
    implicit none

    private

    character(len=512) :: ncfname
    integer            :: ncid
    integer            :: t_axis_id, t_dim_id, n_writes

    type netcdf_stat_info
        character(32)    :: name      = ''
        character(128)   :: long_name = ''
        character(128)   :: std_name  = ''
        character(16)    :: unit      = ''
        integer          :: dtype     = -1
        integer          :: varid     = -1
        double precision :: val       = zero   ! data to be written
    end type netcdf_stat_info

    integer, parameter :: NC_KE     = 1     &
                        , NC_EN     = 2     &
                        , NC_OMAX   = 3     &
                        , NC_ORMS   = 4     &
                        , NC_OCHAR  = 5     &
                        , NC_OXMEAN = 6     &
                        , NC_OYMEAN = 7     &
                        , NC_OZMEAN = 8     &
                        , NC_KEXY   = 9     &
                        , NC_KEZ    = 10    &
                        , NC_DIVXY2 = 11    &
                        , NC_ENXY   = 12    &
                        , NC_ENZ    = 13
#ifdef ENABLE_BUOYANCY
    integer, parameter :: NC_APE    = 14    &
                        , NC_BMAX   = 15    &
                        , NC_BMIN   = 16

#ifdef ENABLE_PERTURBATION_MODE
    integer, parameter :: NC_BASQ   = 17
    type(netcdf_stat_info) :: nc_dset(NC_BASQ)
#else
    type(netcdf_stat_info) :: nc_dset(NC_BMIN)
#endif
#else
    type(netcdf_stat_info) :: nc_dset(NC_ENZ)
#endif


    double precision   :: restart_time

    integer :: field_stats_io_timer

    public :: create_netcdf_field_stats_file,   &
              write_netcdf_field_stats,         &
              field_stats_io_timer,             &
              set_netcdf_field_diagnostic,      &
              NC_OMAX,                          &
              NC_ORMS,                          &
              NC_OCHAR,                         &
              NC_OXMEAN,                        &
              NC_OYMEAN,                        &
              NC_OZMEAN

    contains

        ! Create the NetCDF field diagnostic file.
        ! @param[in] basename of the file
        ! @param[in] overwrite the file
        subroutine create_netcdf_field_stats_file(basename, overwrite)
            character(*), intent(in)  :: basename
            logical,      intent(in)  :: overwrite
            logical                   :: l_exist
            integer                   :: n

            if (world%rank .ne. world%root) then
                return
            endif

            call set_netcdf_stat_info

            ncfname =  basename // '_field_stats.nc'

            restart_time = -one
            n_writes = 1

            call exist_netcdf_file(ncfname, l_exist)

            if (l_exist) then
                call open_netcdf_file(ncfname, NF90_NOWRITE, ncid, l_serial=.true.)
                call get_num_steps(ncid, n_writes)
                if (n_writes > 0) then
                   call get_time(ncid, restart_time)
                   call read_netcdf_field_stats_content
                   call close_netcdf_file(ncid, l_serial=.true.)
                   n_writes = n_writes + 1
                   return
                else
                   call close_netcdf_file(ncid, l_serial=.true.)
                   call delete_netcdf_file(ncfname)
                endif
            endif

            call create_netcdf_file(ncfname, overwrite, ncid, l_serial=.true.)

            call write_netcdf_info(ncid=ncid,                       &
                                   version_tag=package_version,     &
                                   file_type='field_stats',         &
                                   cf_version=cf_version)

            call write_netcdf_box(ncid, lower, extent, (/nx, ny, nz/))

            call write_physical_quantities(ncid)

            call write_netcdf_options(ncid)

            call define_netcdf_temporal_dimension(ncid, t_dim_id, t_axis_id)

            ! define statitics
            do n = 1, size(nc_dset)
                call define_netcdf_dataset(ncid=ncid,                       &
                                           name=nc_dset(n)%name,            &
                                           long_name=nc_dset(n)%long_name,  &
                                           std_name=nc_dset(n)%std_name,    &
                                           unit=nc_dset(n)%unit,            &
                                           dtype=nc_dset(n)%dtype,          &
                                           dimids=(/t_dim_id/),             &
                                           varid=nc_dset(n)%varid)
            enddo
            call close_definition(ncid)

            call close_netcdf_file(ncid, l_serial=.true.)

        end subroutine create_netcdf_field_stats_file

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Pre-condition: Assumes an open file
        subroutine read_netcdf_field_stats_content
            integer :: n

            call get_dim_id(ncid, 't', t_dim_id)

            call get_var_id(ncid, 't', t_axis_id)

            do n = 1, size(nc_dset)
                call get_var_id(ncid, nc_dset(n)%name, nc_dset(n)%varid)
            enddo

        end subroutine read_netcdf_field_stats_content

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Write a step in the field diagnostic file.
        ! @param[in] t is the time
        ! @param[in] dt is the time step
        subroutine write_netcdf_field_stats(t)
            double precision, intent(in) :: t
            integer                      :: n

            call start_timer(field_stats_io_timer)

            call update_netcdf_field_diagnostics

            if (world%rank /= world%root) then
                call stop_timer(field_stats_io_timer)
                return
            endif

            if (t <= restart_time) then
                call stop_timer(field_stats_io_timer)
                return
            endif

            call open_netcdf_file(ncfname, NF90_WRITE, ncid, l_serial=.true.)

            ! write time
            call write_netcdf_scalar(ncid, t_axis_id, t, n_writes, l_serial=.true.)

            !
            ! write diagnostics
            !
            do n = 1, size(nc_dset)
                call write_netcdf_scalar(ncid,              &
                                         nc_dset(n)%varid,  &
                                         nc_dset(n)%val,    &
                                         n_writes,          &
                                         l_serial=.true.)
            enddo

            ! increment counter
            n_writes = n_writes + 1

            call close_netcdf_file(ncid, l_serial=.true.)

            call stop_timer(field_stats_io_timer)

        end subroutine write_netcdf_field_stats

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine update_netcdf_field_diagnostics
#ifdef ENABLE_BUOYANCY
            double precision             :: bmin, bmax
#ifdef ENABLE_PERTURBATION_MODE
            integer                      :: iz
#endif
#endif
            nc_dset(NC_KE)%val = get_kinetic_energy()
            nc_dset(NC_EN)%val = get_enstrophy()

#ifdef ENABLE_BUOYANCY
            call field_combine_physical(sbuoy, buoy)
            nc_dset(NC_APE)%val = get_available_potential_energy()
#ifdef ENABLE_PERTURBATION_MODE
            do iz = 0, nz
                buoy(iz, :, :) = buoy(iz, :, :) + bbarz(iz)
            enddo
#endif
            bmin = minval(buoy)
            bmax = maxval(buoy)

            call mpi_blocking_reduce(bmin, MPI_MIN, world)
            call mpi_blocking_reduce(bmax, MPI_MAX, world)

            nc_dset(NC_BMAX)%val = bmax
            nc_dset(NC_BMIN)%val = bmin
#endif

        nc_dset(NC_KEXY)%val   = get_horizontal_kinetic_energy()
        nc_dset(NC_KEZ)%val    = get_vertical_kinetic_energy()
        nc_dset(NC_ENXY)%val   = get_horizontal_enstrophy()
        nc_dset(NC_ENZ)%val    = get_vertical_enstrophy()
        nc_dset(NC_DIVXY2)%val = get_squared_horizontal_divergence()

#ifdef ENABLE_BUOYANCY_PERTURBATION_MODE
        nc_dset(NC_BASQ)%val = get_squared_buoyancy_anomaly()
#endif

        end subroutine update_netcdf_field_diagnostics

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine set_netcdf_field_diagnostic(val, n)
            double precision, intent(in) :: val
            integer,          intent(in) :: n

            nc_dset(n)%val = val

        end subroutine set_netcdf_field_diagnostic

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine set_netcdf_stat_info
            nc_dset(NC_KE) = netcdf_stat_info(                          &
                name='ke',                                              &
                long_name='domain-averaged kinetic energy',             &
                std_name='',                                            &
                unit='m^2/s^2',                                         &
                dtype=NF90_DOUBLE)

            nc_dset(NC_EN) = netcdf_stat_info(                          &
                name='en',                                              &
                long_name='domain-averaged enstrophy',                  &
                std_name='',                                            &
                unit='1/s^2',                                           &
                dtype=NF90_DOUBLE)

            nc_dset(NC_OMAX) = netcdf_stat_info(                        &
                name='vortmax',                                         &
                long_name='maximum vorticity magnitude',                &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            nc_dset(NC_ORMS) = netcdf_stat_info(                        &
                name='vortrms',                                         &
                long_name='root-mean square vorticity',                 &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            nc_dset(NC_OCHAR) = netcdf_stat_info(                       &
                name='vorch',                                           &
                long_name='characteristic vorticity',                   &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            nc_dset(NC_OXMEAN) = netcdf_stat_info(                      &
                name='x_vormean',                                       &
                long_name='mean x-vorticity',                           &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            nc_dset(NC_OYMEAN) = netcdf_stat_info(                      &
                name='y_vormean',                                       &
                long_name='mean y-vorticity',                           &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            nc_dset(NC_OZMEAN) = netcdf_stat_info(                      &
                name='z_vormean',                                       &
                long_name='mean z-vorticity',                           &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            nc_dset(NC_KEXY) = netcdf_stat_info(                        &
                name='kexy',                                            &
                long_name='domain-averaged horizontal kinetic energy',  &
                std_name='',                                            &
                unit='m^2/s^2',                                         &
                dtype=NF90_DOUBLE)

            nc_dset(NC_KEZ) = netcdf_stat_info(                         &
                name='kez',                                             &
                long_name='domain-averaged vertical kinetic energy',    &
                std_name='',                                            &
                unit='m^2/s^2',                                         &
                dtype=NF90_DOUBLE)

            nc_dset(NC_ENXY) = netcdf_stat_info(                        &
                name='enxy',                                            &
                long_name='domain-averaged horizontal enstrophy',       &
                std_name='',                                            &
                unit='1/s^2',                                           &
                dtype=NF90_DOUBLE)

            nc_dset(NC_ENZ) = netcdf_stat_info(                         &
                name='enz',                                             &
                long_name='domain-averaged vertical enstrophy',         &
                std_name='',                                            &
                unit='1/s^2',                                           &
                dtype=NF90_DOUBLE)

            nc_dset(NC_DIVXY2) = netcdf_stat_info(                          &
                name='divxy2',                                              &
                long_name='domain-averaged squared horizontal divergence',  &
                std_name='',                                                &
                unit='1/s^2',                                               &
                dtype=NF90_DOUBLE)

#ifdef ENABLE_BUOYANCY
            nc_dset(NC_APE) = netcdf_stat_info(                         &
                name='ape',                                             &
                long_name='domain-averaged available potential energy', &
                std_name='',                                            &
                unit='m^2/s^2',                                         &
                dtype=NF90_DOUBLE)

            nc_dset(NC_BMIN) = netcdf_stat_info(                        &
                name='min_buoyancy',                                    &
                long_name='minimum buoyancy',                           &
                std_name='',                                            &
                unit='m/s^2',                                           &
                dtype=NF90_DOUBLE)

            nc_dset(NC_BMAX) = netcdf_stat_info(                        &
                name='max_buoyancy',                                    &
                long_name='maximum buoyancy',                           &
                std_name='',                                            &
                unit='m/s^2',                                           &
                dtype=NF90_DOUBLE)

#ifdef ENABLE_PERTURBATION_MODE
            nc_dset(NC_BASQ) = netcdf_stat_info(                        &
                name='squared_buoyancy_anomaly',                        &
                long_name='domain-averaged squared buoyany anomaly',    &
                std_name='',                                            &
                unit='m^2/s^4',                                         &
                dtype=NF90_DOUBLE)
#endif
#endif

        end subroutine set_netcdf_stat_info

end module field_diagnostics_netcdf
