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
#ifdef ENABLE_BALANCE
    use field_balance, only : balance_fields      &
                            , kebal, keubal       &
                            , apebal, apeubal
#endif
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
                        , NC_ENZ    = 13    &
                        , NC_OXMIN  = 14    &
                        , NC_OYMIN  = 15    &
                        , NC_OZMIN  = 16    &
                        , NC_OXMAX  = 17    &
                        , NC_OYMAX  = 18    &
                        , NC_OZMAX  = 19    &
                        , NC_HEMAX  = 20    &
                        , NC_GMAX   = 21    &
                        , NC_RGMAX  = 22    &
                        , NC_RIMIN  = 23    &
                        , NC_ROMIN  = 24
#ifdef ENABLE_BUOYANCY
    integer, parameter :: NC_APE     = 25    &
                        , NC_BMAX    = 26    &
                        , NC_BMIN    = 27    &
                        , NC_MSS     = 28    &  ! mss = minimum static stability
                        , NC_KEBAL   = 29    &
                        , NC_KEUBAL  = 30    &
                        , NC_APEBAL  = 31    &
                        , NC_APEUBAL = 32
    type(netcdf_stat_info) :: nc_dset(NC_APEUBAL)
#else
    type(netcdf_stat_info) :: nc_dset(NC_ROMIN)
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
              NC_OZMEAN,                        &
              NC_GMAX,                          &
              NC_RGMAX

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

        !@Pre Assumes all fields are up-to-date
        subroutine update_netcdf_field_diagnostics
#ifdef ENABLE_BALANCE
            use options, only : output
#endif
            integer          :: nc
#ifdef ENABLE_BUOYANCY
            double precision :: tbuoy(0:nz,                 & ! total buoyancy
                                      box%lo(2):box%hi(2),  &
                                      box%lo(1):box%hi(1))
            double precision :: bmin, bmax
            double precision :: buf(12) = zero
            integer          :: iz
#else
            double precision :: buf(7) = zero
#endif

#ifdef ENABLE_BUOYANCY
            call field_combine_physical(sbuoy, buoy)

            ! get total buoyancy:
            do iz = 0, nz
                tbuoy(iz, :, :) = buoy(iz, :, :) + bbarz(iz)
            enddo

            bmin = minval(tbuoy)
            bmax = maxval(tbuoy)
#endif

            !
            ! Summed values
            !

            buf(1) = get_kinetic_energy(vel, l_global=.false., l_allreduce=.false.)
            buf(2) = get_enstrophy(l_global=.false., l_allreduce=.false.)
            buf(3) = get_horizontal_kinetic_energy(vel, l_global=.false.)
            buf(4) = get_vertical_kinetic_energy(l_global=.false.)
            buf(5) = get_horizontal_enstrophy(l_global=.false.)
            buf(6) = get_vertical_enstrophy(l_global=.false.)
            buf(7) = get_squared_horizontal_divergence(l_global=.false.)
            buf(8) = get_available_potential_energy(buoy, l_global=.false., l_allreduce=.false.)

#ifdef ENABLE_BALANCE
            if (output%l_balanced) then
                call balance_fields(l_global=.false.)
                buf(9) = kebal
                buf(10) = keubal
                buf(11) = apebal
                buf(12) = apeubal
            endif
#endif

            call mpi_blocking_reduce(buf, MPI_SUM, world)

            nc_dset(NC_KE)%val     = buf(1)
            nc_dset(NC_EN)%val     = buf(2)
            nc_dset(NC_KEXY)%val   = buf(3)
            nc_dset(NC_KEZ)%val    = buf(4)
            nc_dset(NC_ENXY)%val   = buf(5)
            nc_dset(NC_ENZ)%val    = buf(6)
            nc_dset(NC_DIVXY2)%val = buf(7)
#ifdef ENABLE_BUOYANCY
            nc_dset(NC_APE)%val    = buf(8)

            if (output%l_balanced) then
                nc_dset(NC_KEBAL)%val   = buf(9)
                nc_dset(NC_KEUBAL)%val  = buf(10)
                nc_dset(NC_APEBAL)%val  = buf(11)
                nc_dset(NC_APEUBAL)%val = buf(12)
            endif
#endif


            !
            ! Minimum values
            !

            do nc = 1, 3
                buf(nc) = minval(vor(:, :, :, nc))
            enddo

            buf(4) = get_min_rossby_number(l_global=.false.)

#ifdef ENABLE_BUOYANCY
            buf(5) = get_min_richardson_number(l_global=.false.)

            buf(6) = bmin

            buf(7) = get_minimum_static_stability(l_global=.false.)
#endif

            call mpi_blocking_reduce(buf(1:7), MPI_MIN, world)


            nc_dset(NC_OXMIN)%val = buf(1)
            nc_dset(NC_OYMIN)%val = buf(2)
            nc_dset(NC_OZMIN)%val = buf(3)
            nc_dset(NC_ROMIN)%val = buf(4)
#ifdef ENABLE_BUOYANCY
            nc_dset(NC_RIMIN)%val = buf(5)
            nc_dset(NC_BMIN)%val  = buf(6)
            nc_dset(NC_MSS)%val   = buf(7)
#endif

            !
            ! Maximum values
            !

            do nc = 1, 3
                buf(nc) = maxval(vor(:, :, :, nc))
            enddo

            buf(4) = get_max_horizontal_enstrophy(l_global=.false.)
#ifdef ENABLE_BUOYANCY
            buf(5) = bmax
#endif

            call mpi_blocking_reduce(buf(1:5), MPI_MAX, world)

            nc_dset(NC_OXMAX)%val = buf(1)
            nc_dset(NC_OYMAX)%val = buf(2)
            nc_dset(NC_OZMAX)%val = buf(3)
            nc_dset(NC_HEMAX)%val = buf(4)
#ifdef ENABLE_BUOYANCY
            nc_dset(NC_BMAX)%val  = buf(5)
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
#ifdef ENABLE_BUOYANCY
            use options, only : output
#endif

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

            nc_dset(NC_OXMIN) = netcdf_stat_info(                       &
                name='x_vormin',                                        &
                long_name='min x-vorticity',                            &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            nc_dset(NC_OYMIN) = netcdf_stat_info(                       &
                name='y_vormin',                                        &
                long_name='min y-vorticity',                            &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            nc_dset(NC_OZMIN) = netcdf_stat_info(                       &
                name='z_vormin',                                        &
                long_name='min z-vorticity',                            &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            nc_dset(NC_OXMAX) = netcdf_stat_info(                       &
                name='x_vormax',                                        &
                long_name='max x-vorticity',                            &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            nc_dset(NC_OYMAX) = netcdf_stat_info(                       &
                name='y_vormax',                                        &
                long_name='max y-vorticity',                            &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            nc_dset(NC_OZMAX) = netcdf_stat_info(                       &
                name='z_vormax',                                        &
                long_name='max z-vorticity',                            &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            nc_dset(NC_HEMAX) = netcdf_stat_info(                       &
                name='enxy_max',                                        &
                long_name='max horizontal enstrophy',                   &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            nc_dset(NC_GMAX) = netcdf_stat_info(                        &
                name='gmax',                                            &
                long_name='maximum gamma',                              &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

#ifndef ENABLE_SMAGORINSKY
            nc_dset(NC_RGMAX) = netcdf_stat_info(                       &
                name='rolling_mean_gmax',                               &
                long_name='rolling mean maximum gamma',                 &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)
#endif

            nc_dset(NC_RIMIN) = netcdf_stat_info(                       &
                name='ri_min',                                          &
                long_name='minimum Richardson number',                  &
                std_name='',                                            &
                unit='1',                                               &
                dtype=NF90_DOUBLE)

            nc_dset(NC_ROMIN) = netcdf_stat_info(                       &
                name='ro_min',                                          &
                long_name='minimum Rossby number',                      &
                std_name='',                                            &
                unit='1',                                               &
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

            nc_dset(NC_MSS) = netcdf_stat_info(                         &
                name='minimum_static_stability',                        &
                long_name='minimum static stability',                   &
                std_name='',                                            &
                unit='1',                                               &
                dtype=NF90_DOUBLE)

            if (output%l_balanced) then
                nc_dset(NC_KEBAL) = netcdf_stat_info(                       &
                    name='kebal',                                           &
                    long_name='domain-averaged balanced kinetic energy',    &
                    std_name='',                                            &
                    unit='m^2/s^2',                                         &
                    dtype=NF90_DOUBLE)

                nc_dset(NC_KEUBAL) = netcdf_stat_info(                      &
                    name='keubal',                                          &
                    long_name='domain-averaged imbalanced kinetic energy',  &
                    std_name='',                                            &
                    unit='m^2/s^2',                                         &
                    dtype=NF90_DOUBLE)

                nc_dset(NC_APEBAL) = netcdf_stat_info(                      &
                    name='apebal',                                          &
                    long_name='domain-averaged balanced APE',               &
                    std_name='',                                            &
                    unit='m^2/s^2',                                         &
                    dtype=NF90_DOUBLE)

                nc_dset(NC_APEUBAL) = netcdf_stat_info(                     &
                    name='apeubal',                                         &
                    long_name='domain-averaged imbalanced APE',             &
                    std_name='',                                            &
                    unit='m^2/s^2',                                         &
                    dtype=NF90_DOUBLE)
            endif
#endif

        end subroutine set_netcdf_stat_info

end module field_diagnostics_netcdf
