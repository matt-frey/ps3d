! =============================================================================
!                   Write field diagnostics to NetCDF.
! =============================================================================
module field_diagnostics_netcdf
    use fields
    use netcdf_utils
    use netcdf_writer
    use netcdf_reader
    use constants, only : zero, one
    use parameters, only : nx, ny, nz    &
                         , write_netcdf_parameters
    use config, only : package_version, cf_version
    use mpi_timer, only : start_timer, stop_timer
    use options, only : write_netcdf_options
    use physics, only : write_physical_quantities
    use mpi_collectives, only : mpi_blocking_reduce
    use field_diagnostics
    use fields_derived, only : delta
#if defined(ENABLE_BALANCE) && defined(ENABLE_BUOYANCY)
    use field_balance, only : balance_fields      &
                            , kebal, keubal       &
                            , apebal, apeubal
#endif
    implicit none

    private

    character(len=512) :: ncfname
    integer            :: ncid
    integer            :: t_axis_id, t_dim_id, n_writes

    type, extends(netcdf_info) :: netcdf_stat_info
        double precision :: val       = zero   ! data to be written
    end type netcdf_stat_info

    integer, parameter :: NC_KE       = 1   &
                        , NC_EN       = 2   &
                        , NC_OMAX     = 3   &
                        , NC_ORMS     = 4   &
                        , NC_OCHAR    = 5   &
                        , NC_OXMEAN   = 6   &
                        , NC_OYMEAN   = 7   &
                        , NC_OZMEAN   = 8   &
                        , NC_KEXY     = 9   &
                        , NC_KEZ      = 10  &
                        , NC_ENXY     = 11  &
                        , NC_ENZ      = 12  &
                        , NC_OXMIN    = 13  &
                        , NC_OYMIN    = 14  &
                        , NC_OZMIN    = 15  &
                        , NC_OXMAX    = 16  &
                        , NC_OYMAX    = 17  &
                        , NC_OZMAX    = 18  &
                        , NC_HEMAX    = 19  &
                        , NC_GMAX     = 20  &
                        , NC_BFMAX    = 21  &
                        , NC_UMAX     = 22  &
                        , NC_VMAX     = 23  &
                        , NC_WMAX     = 24  &
                        , NC_USOXMAX  = 25  &
                        , NC_LSOXMAX  = 26  &
                        , NC_USOYMAX  = 27  &
                        , NC_LSOYMAX  = 28  &
                        , NC_USOZMAX  = 29  &
                        , NC_LSOZMAX  = 30  &
                        , NC_USUHMAX  = 31  &
                        , NC_USGMAX   = 32  &
                        , NC_LSGMAX   = 33  &
                        , NC_USZRMS   = 34  &
                        , NC_USDELRMS = 35  &
                        , NC_RGMAX    = 36  &
                        , NC_RBFMAX   = 37  &
                        , NC_RIMIN    = 38  &
                        , NC_ROMIN    = 39  &
                        , NC_ROMAX    = 40
#ifdef ENABLE_BUOYANCY
    integer, parameter :: NC_APE      = 41  &
                        , NC_BMAX     = 42  &
                        , NC_BMIN     = 43  &
                        , NC_BUSMIN   = 44  &
                        , NC_BUSMAX   = 45  &
                        , NC_BLSMIN   = 46  &
                        , NC_BLSMAX   = 47  &
                        , NC_MSS      = 48     ! mss = minimum static stability
#ifdef ENABLE_BALANCE
                        , NC_KEBAL    = 49  &
                        , NC_KEUBAL   = 50  &
                        , NC_APEBAL   = 51  &
                        , NC_APEUBAL  = 52
    type(netcdf_stat_info) :: nc_dset(NC_APEUBAL)
#else
    type(netcdf_stat_info) :: nc_dset(NC_MSS)
#endif
#else
    type(netcdf_stat_info) :: nc_dset(NC_ROMAX)
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
              NC_BFMAX,                         &
              NC_UMAX,                          &
              NC_VMAX,                          &
              NC_WMAX,                          &
              NC_RGMAX,                         &
              NC_RBFMAX,                        &
              NC_USGMAX,                        &
              NC_LSGMAX

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

            call write_netcdf_parameters(ncid)

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
            double precision :: busmin, busmax, blsmin, blsmax
            double precision :: buf(15) = zero
            integer          :: iz
#else
            double precision :: buf(13) = zero
#endif

#ifdef ENABLE_BUOYANCY
            call layout%combine_physical(sbuoy, buoy)

            ! get total buoyancy:
            do iz = 0, nz
                tbuoy(iz, :, :) = buoy(iz, :, :) + bbarz(iz)
            enddo

            bmin = minval(tbuoy)
            bmax = maxval(tbuoy)
            busmin = minval(tbuoy(nz, :, :))
            busmax = maxval(tbuoy(nz, :, :))
            blsmin = minval(tbuoy(0,  :, :))
            blsmax = maxval(tbuoy(0,  :, :))
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


            ! rms of upper surface z-vorticity
            ! (we must take the square-root after the MPI reduction)
            buf(7) = sum(vor(nz, :, :, 3) ** 2) / dble(nx * ny)

            ! rms surface divergence
            ! (we must take the square-root after the MPI reduction)
            buf(8) = sum(delta(nz, :, :) ** 2) / dble(nx * ny)



#ifdef ENABLE_BUOYANCY
            buf(9) = get_available_potential_energy(buoy, l_global=.false., l_allreduce=.false.)

#ifdef ENABLE_BALANCE
            if (output%l_balanced) then
                call balance_fields(l_global=.false.)
                buf(10) = kebal
                buf(11) = keubal
                buf(12) = apebal
                buf(13) = apeubal
            endif
#endif
#endif

            call mpi_blocking_reduce(buf, MPI_SUM, world)

            nc_dset(NC_KE)%val       = buf(1) * ncelli
            nc_dset(NC_EN)%val       = buf(2) * ncelli
            nc_dset(NC_KEXY)%val     = buf(3) * ncelli
            nc_dset(NC_KEZ)%val      = buf(4) * ncelli
            nc_dset(NC_ENXY)%val     = buf(5) * ncelli
            nc_dset(NC_ENZ)%val      = buf(6) * ncelli
            nc_dset(NC_USZRMS)%val   = sqrt(buf(7))
            nc_dset(NC_USDELRMS)%val = sqrt(buf(8))
#ifdef ENABLE_BUOYANCY
            nc_dset(NC_APE)%val    = buf(9)

#ifdef ENABLE_BALANCE
            if (output%l_balanced) then
                nc_dset(NC_KEBAL)%val   = buf(10)
                nc_dset(NC_KEUBAL)%val  = buf(11)
                nc_dset(NC_APEBAL)%val  = buf(12)
                nc_dset(NC_APEUBAL)%val = buf(13)
            endif
#endif
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

            buf(7) = busmin
            buf(8) = blsmin

            buf(9) = get_minimum_static_stability(l_global=.false.)

            call mpi_blocking_reduce(buf(1:9), MPI_MIN, world)
#else
            call mpi_blocking_reduce(buf(1:4), MPI_MIN, world)
#endif



            nc_dset(NC_OXMIN)%val = buf(1)
            nc_dset(NC_OYMIN)%val = buf(2)
            nc_dset(NC_OZMIN)%val = buf(3)
            nc_dset(NC_ROMIN)%val = buf(4)
#ifdef ENABLE_BUOYANCY
            nc_dset(NC_RIMIN)%val  = buf(5)
            nc_dset(NC_BMIN)%val   = buf(6)
            nc_dset(NC_BUSMIN)%val = buf(7)
            nc_dset(NC_BLSMIN)%val = buf(8)
            nc_dset(NC_MSS)%val    = buf(9)
#endif

            !
            ! Maximum values
            !

            do nc = 1, 3
                buf(nc) = maxval(vor(:, :, :, nc))
            enddo

            buf(4) = get_max_horizontal_enstrophy(l_global=.false.)

            buf(5)  = maxval(vor(nz, :, :, 1))
            buf(6)  = maxval(vor(0,  :, :, 1))
            buf(7)  = maxval(vor(nz, :, :, 2))
            buf(8)  = maxval(vor(0,  :, :, 2))
            buf(9)  = maxval(vor(nz, :, :, 3))
            buf(10) = maxval(vor(0,  :, :, 3))
            buf(11) = sqrt(maxval(vel(nz, :, :, 1) ** 2 + &
                                   vel(nz, :, :, 2) ** 2))

            buf(12) = get_max_rossby_number(l_global=.false.)

#ifdef ENABLE_BUOYANCY
            buf(13) = bmax
            buf(14) = busmax
            buf(15) = blsmax

            call mpi_blocking_reduce(buf(1:15), MPI_MAX, world)
#else
            call mpi_blocking_reduce(buf(1:12), MPI_MAX, world)
#endif

            nc_dset(NC_OXMAX)%val    = buf(1)
            nc_dset(NC_OYMAX)%val    = buf(2)
            nc_dset(NC_OZMAX)%val    = buf(3)
            nc_dset(NC_HEMAX)%val    = buf(4)
            nc_dset(NC_USOXMAX)%val  = buf(5)
            nc_dset(NC_LSOXMAX)%val  = buf(6)
            nc_dset(NC_USOYMAX)%val  = buf(7)
            nc_dset(NC_LSOYMAX)%val  = buf(8)
            nc_dset(NC_USOZMAX)%val  = buf(9)
            nc_dset(NC_LSOZMAX)%val  = buf(10)
            nc_dset(NC_USUHMAX)%val  = buf(11)
            nc_dset(NC_ROMAX)%val    = buf(12)

#ifdef ENABLE_BUOYANCY
            nc_dset(NC_BMAX)%val   = buf(13)
            nc_dset(NC_BUSMAX)%val = buf(14)
            nc_dset(NC_BLSMAX)%val = buf(15)
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

            call nc_dset(NC_KE)%set_info(                               &
                name='ke',                                              &
                long_name='domain-averaged kinetic energy',             &
                std_name='',                                            &
                unit='m^2/s^2',                                         &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_EN)%set_info(                               &
                name='en',                                              &
                long_name='domain-averaged enstrophy',                  &
                std_name='',                                            &
                unit='1/s^2',                                           &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_OMAX)%set_info(                             &
                name='vortmax',                                         &
                long_name='maximum vorticity magnitude',                &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_ORMS)%set_info(                             &
                name='vortrms',                                         &
                long_name='root-mean square vorticity',                 &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_OCHAR)%set_info(                            &
                name='vorch',                                           &
                long_name='characteristic vorticity',                   &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_OXMEAN)%set_info(                           &
                name='x_vormean',                                       &
                long_name='mean x-vorticity',                           &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_OYMEAN)%set_info(                           &
                name='y_vormean',                                       &
                long_name='mean y-vorticity',                           &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_OZMEAN)%set_info(                           &
                name='z_vormean',                                       &
                long_name='mean z-vorticity',                           &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_OXMIN)%set_info(                            &
                name='x_vormin',                                        &
                long_name='min x-vorticity',                            &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_OYMIN)%set_info(                            &
                name='y_vormin',                                        &
                long_name='min y-vorticity',                            &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_OZMIN)%set_info(                            &
                name='z_vormin',                                        &
                long_name='min z-vorticity',                            &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_OXMAX)%set_info(                            &
                name='x_vormax',                                        &
                long_name='max x-vorticity',                            &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_OYMAX)%set_info(                            &
                name='y_vormax',                                        &
                long_name='max y-vorticity',                            &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_OZMAX)%set_info(                            &
                name='z_vormax',                                        &
                long_name='max z-vorticity',                            &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_HEMAX)%set_info(                            &
                name='enxy_max',                                        &
                long_name='max horizontal enstrophy',                   &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_GMAX)%set_info(                             &
                name='gmax',                                            &
                long_name='maximum gamma',                              &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

             call nc_dset(NC_BFMAX)%set_info(                           &
                name='bfmax',                                           &
                long_name='maximum buoyancy frequency',                 &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_UMAX)%set_info(                             &
                name='umax',                                            &
                long_name='maximum x-velocity',                         &
                std_name='',                                            &
                unit='m/s',                                             &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_VMAX)%set_info(                             &
                name='vmax',                                            &
                long_name='maximum y-velocity',                         &
                std_name='',                                            &
                unit='m/s',                                             &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_WMAX)%set_info(                             &
                name='wmax',                                            &
                long_name='maximum z-velocity',                         &
                std_name='',                                            &
                unit='m/s',                                             &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_USOXMAX)%set_info(                          &
                name='uoxmax',                                          &
                long_name='upper surface maximum x-vorticity',          &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_LSOXMAX)%set_info(                          &
                name='loxmax',                                          &
                long_name='lower surface maximum x-vorticity',          &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_USOYMAX)%set_info(                          &
                name='uoymax',                                          &
                long_name='upper surface maximum y-vorticity',          &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_LSOYMAX)%set_info(                          &
                name='loymax',                                          &
                long_name='lower surface maximum y-vorticity',          &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_USOZMAX)%set_info(                          &
                name='uozmax',                                          &
                long_name='upper surface maximum z-vorticity',          &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_LSOZMAX)%set_info(                          &
                name='lozmax',                                          &
                long_name='lower surface maximum z-vorticity',          &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_USUHMAX)%set_info(                          &
                name='usuhmax',                                         &
                long_name='upper surface maximum horizontal speed',     &
                std_name='',                                            &
                unit='m/s',                                             &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_USGMAX)%set_info(                           &
                name='usgmax',                                          &
                long_name='upper surface maximum strain',               &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_LSGMAX)%set_info(                           &
                name='lsgmax',                                          &
                long_name='lower surface maximum strain',               &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_USZRMS)%set_info(                           &
                name='uszrms',                                          &
                long_name='rms of upper surface z-vorticity',           &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_USDELRMS)%set_info(                         &
                name='usdeltarms',                                      &
                long_name='rms of upper surface horizontal divergence', &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_RGMAX)%set_info(                            &
                name='rolling_mean_gmax',                               &
                long_name='rolling mean maximum gamma',                 &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_RBFMAX)%set_info(                           &
                name='rolling_mean_bfmax',                              &
                long_name='rolling mean maximum buoyancy frequency',    &
                std_name='',                                            &
                unit='1/s',                                             &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_RIMIN)%set_info(                            &
                name='ri_min',                                          &
                long_name='minimum Richardson number',                  &
                std_name='',                                            &
                unit='1',                                               &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_ROMIN)%set_info(                            &
                name='ro_min',                                          &
                long_name='minimum Rossby number',                      &
                std_name='',                                            &
                unit='1',                                               &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_ROMAX)%set_info(                            &
                name='ro_max',                                          &
                long_name='maximum Rossby number',                      &
                std_name='',                                            &
                unit='1',                                               &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_KEXY)%set_info(                             &
                name='kexy',                                            &
                long_name='domain-averaged horizontal kinetic energy',  &
                std_name='',                                            &
                unit='m^2/s^2',                                         &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_KEZ)%set_info(                              &
                name='kez',                                             &
                long_name='domain-averaged vertical kinetic energy',    &
                std_name='',                                            &
                unit='m^2/s^2',                                         &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_ENXY)%set_info(                             &
                name='enxy',                                            &
                long_name='domain-averaged horizontal enstrophy',       &
                std_name='',                                            &
                unit='1/s^2',                                           &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_ENZ)%set_info(                              &
                name='enz',                                             &
                long_name='domain-averaged vertical enstrophy',         &
                std_name='',                                            &
                unit='1/s^2',                                           &
                dtype=NF90_DOUBLE)

#ifdef ENABLE_BUOYANCY
            call nc_dset(NC_APE)%set_info(                              &
                name='ape',                                             &
                long_name='domain-averaged available potential energy', &
                std_name='',                                            &
                unit='m^2/s^2',                                         &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_BMIN)%set_info(                             &
                name='min_buoyancy',                                    &
                long_name='minimum buoyancy',                           &
                std_name='',                                            &
                unit='m/s^2',                                           &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_BUSMIN)%set_info(                           &
                name='busmin',                                          &
                long_name='minimum upper surface buoyancy',             &
                std_name='',                                            &
                unit='m/s^2',                                           &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_BUSMAX)%set_info(                           &
                name='busmax',                                          &
                long_name='maximum upper surface buoyancy',             &
                std_name='',                                            &
                unit='m/s^2',                                           &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_BLSMIN)%set_info(                           &
                name='blsmin',                                          &
                long_name='minimum lower surface buoyancy',             &
                std_name='',                                            &
                unit='m/s^2',                                           &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_BLSMAX)%set_info(                           &
                name='blsmax',                                          &
                long_name='maximum lower surface buoyancy',             &
                std_name='',                                            &
                unit='m/s^2',                                           &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_BMAX)%set_info(                             &
                name='max_buoyancy',                                    &
                long_name='maximum buoyancy',                           &
                std_name='',                                            &
                unit='m/s^2',                                           &
                dtype=NF90_DOUBLE)

            call nc_dset(NC_MSS)%set_info(                              &
                name='minimum_static_stability',                        &
                long_name='minimum static stability',                   &
                std_name='',                                            &
                unit='1',                                               &
                dtype=NF90_DOUBLE)

#ifdef ENABLE_BALANCE
            if (output%l_balanced) then
                call nc_dset(NC_KEBAL)%set_info(                            &
                    name='kebal',                                           &
                    long_name='domain-averaged balanced kinetic energy',    &
                    std_name='',                                            &
                    unit='m^2/s^2',                                         &
                    dtype=NF90_DOUBLE)

                call nc_dset(NC_KEUBAL)%set_info(                           &
                    name='keubal',                                          &
                    long_name='domain-averaged imbalanced kinetic energy',  &
                    std_name='',                                            &
                    unit='m^2/s^2',                                         &
                    dtype=NF90_DOUBLE)

                call nc_dset(NC_APEBAL)%set_info(                           &
                    name='apebal',                                          &
                    long_name='domain-averaged balanced APE',               &
                    std_name='',                                            &
                    unit='m^2/s^2',                                         &
                    dtype=NF90_DOUBLE)

                call nc_dset(NC_APEUBAL)%set_info(                          &
                    name='apeubal',                                         &
                    long_name='domain-averaged imbalanced APE',             &
                    std_name='',                                            &
                    unit='m^2/s^2',                                         &
                    dtype=NF90_DOUBLE)
            endif
#endif
#endif

        end subroutine set_netcdf_stat_info

end module field_diagnostics_netcdf
