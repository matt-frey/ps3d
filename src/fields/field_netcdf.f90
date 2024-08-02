module field_netcdf
    use options, only : output, verbose
    use constants, only : one
    use netcdf_utils
    use netcdf_writer
    use netcdf_reader
    use fields
    use config, only : package_version, cf_version
    use mpi_timer, only : start_timer, stop_timer
    use options, only : write_netcdf_options
    use physics, only : write_physical_quantities
    use parameters, only : lower, extent, dx, nx, ny, nz
    implicit none

    private

    integer :: field_io_timer

    character(len=512) :: ncfname
    integer            :: ncid
    integer            :: dimids(4)     ! = (x, y, z)
    integer            :: coord_ids(3)  ! = (x, y, z)
    integer            :: t_axis_id
    double precision   :: restart_time
    integer            :: n_writes

    type netcdf_field_info
        character(32)  :: name      = ''
        character(128) :: long_name = ''
        character(128) :: std_name  = ''
        character(16)  :: unit      = ''
        integer        :: dtype     = -1
        integer        :: varid     = -1
        logical        :: l_enabled = .false.
    end type netcdf_field_info

    integer, parameter :: NC_X_VEL   = 1        &
                        , NC_Y_VEL   = 2        &
                        , NC_Z_VEL   = 3        &
                        , NC_X_VOR   = 4        &
                        , NC_Y_VOR   = 5        &
                        , NC_Z_VOR   = 6        &
                        , NC_PRES    = 7

#ifdef ENABLE_BUOYANCY
    integer, parameter :: NC_BUOY    = 8

    integer, parameter :: NC_BUOY_AN = 9

    type(netcdf_field_info) :: nc_dset(NC_BUOY_AN)

#else
    type(netcdf_field_info) :: nc_dset(NC_PRES)
#endif

    public :: create_netcdf_field_file  &
            , write_netcdf_fields       &
            , read_netcdf_fields        &
            , field_io_timer


    contains

        ! Create the field file.
        ! @param[in] basename of the file
        ! @param[in] overwrite the file
        subroutine create_netcdf_field_file(basename, overwrite)
            character(*), intent(in)  :: basename
            logical,      intent(in)  :: overwrite
            logical                   :: l_exist
            integer                   :: n

            call set_netcdf_field_output

            ncfname =  basename // '_fields.nc'

            restart_time = -one
            n_writes = 1

            call exist_netcdf_file(ncfname, l_exist)

            if (l_exist) then
                call open_netcdf_file(ncfname, NF90_NOWRITE, ncid)
                call get_num_steps(ncid, n_writes)
                if (n_writes > 0) then
                   call get_time(ncid, restart_time)
                   call read_netcdf_field_content
                   call close_netcdf_file(ncid)
                   n_writes = n_writes + 1
                   return
                else
                   call close_netcdf_file(ncid)
                   if (world%rank == world%root) then
                        call delete_netcdf_file(ncfname)
                   endif
                endif
            endif

            call create_netcdf_file(ncfname, overwrite, ncid)

            ! define global attributes
            call write_netcdf_info(ncid=ncid,                    &
                                   version_tag=package_version,  &
                                   file_type='fields',           &
                                   cf_version=cf_version)

            call write_netcdf_box(ncid, lower, extent, (/nx, ny, nz/))

            call write_physical_quantities(ncid)

            call write_netcdf_options(ncid)

            ! define dimensions
            call define_netcdf_spatial_dimensions_3d(ncid=ncid,                &
                                                     ngps=(/nx, ny, nz+1/),    &
                                                     dimids=dimids(1:3),       &
                                                     axids=coord_ids)


            call define_netcdf_temporal_dimension(ncid, dimids(4), t_axis_id)

            ! define fields
            do n = 1, size(nc_dset)
                if (nc_dset(n)%l_enabled) then
                    call define_netcdf_dataset(ncid=ncid,                       &
                                               name=nc_dset(n)%name,            &
                                               long_name=nc_dset(n)%long_name,  &
                                               std_name=nc_dset(n)%std_name,    &
                                               unit=nc_dset(n)%unit,            &
                                               dtype=nc_dset(n)%dtype,          &
                                               dimids=dimids,                   &
                                               varid=nc_dset(n)%varid)

                endif
            enddo

            call close_definition(ncid)

            call close_netcdf_file(ncid)

        end subroutine create_netcdf_field_file

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Pre-condition: Assumes an open file
        subroutine read_netcdf_field_content
            integer :: n

            call get_dim_id(ncid, 'x', dimids(1))

            call get_dim_id(ncid, 'y', dimids(2))

            call get_dim_id(ncid, 'z', dimids(3))

            call get_dim_id(ncid, 't', dimids(4))


            call get_var_id(ncid, 'x', coord_ids(1))

            call get_var_id(ncid, 'y', coord_ids(2))

            call get_var_id(ncid, 'z', coord_ids(3))

            call get_var_id(ncid, 't', t_axis_id)

            do n = 1, size(nc_dset)
                if (nc_dset(n)%l_enabled) then
                    call get_var_id(ncid, nc_dset(n)%name, nc_dset(n)%varid)
                endif
            enddo

        end subroutine read_netcdf_field_content

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Write a step in the field file.
        ! @param[in] t is the time
        ! @param[in] dt is the time step
        subroutine write_netcdf_fields(t)
            double precision, intent(in) :: t
            integer                      :: cnt(4), start(4)
#if ENABLE_BUOYANCY
            double precision             :: tbuoy(0:nz,                & ! total buoyancy
                                                  box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1))
            integer                      :: iz
#endif

            call start_timer(field_io_timer)

            if (t <= restart_time) then
                call stop_timer(field_io_timer)
                return
            endif

            call open_netcdf_file(ncfname, NF90_WRITE, ncid)

            ! time step to write [step(4) is the time]
            ! need to add 1 since start must begin with index 1
            start(1:3) = box%lo + 1
            start(4) = n_writes

            cnt(1:3) = box%hi - box%lo + 1
            cnt(4)   = 1

            if (n_writes == 1) then
                call write_netcdf_axis_3d(ncid, dimids(1:3), box%lower, dx, &
                                          box%size, start(1:3), cnt(1:3))
            endif

            ! write time
            call write_netcdf_scalar(ncid, t_axis_id, t, n_writes)


            !
            ! write fields
            !
            call write_field_double(NC_X_VEL, vel(:, :, :, 1), start, cnt)
            call write_field_double(NC_Y_VEL, vel(:, :, :, 2), start, cnt)
            call write_field_double(NC_Z_VEL, vel(:, :, :, 3), start, cnt)

            call write_field_double(NC_X_VOR, vor(:, :, :, 1), start, cnt)
            call write_field_double(NC_Y_VOR, vor(:, :, :, 2), start, cnt)
            call write_field_double(NC_Z_VOR, vor(:, :, :, 3), start, cnt)


            call write_field_double(NC_PRES, pres, start, cnt)

#ifdef ENABLE_BUOYANCY
            call fftxys2p(sbuoy, buoy)

            call write_field_double(NC_BUOY_AN, buoy, start, cnt)

            if (nc_dset(NC_BUOY)%l_enabled) then
                ! get total buoyancy
                do iz = 0, nz
                    tbuoy(iz, :, :) = buoy(iz, :, :) + bbarz(iz)
                enddo

                call write_field_double(NC_BUOY, tbuoy, start, cnt)
            endif
#endif

            ! increment counter
            n_writes = n_writes + 1

            call close_netcdf_file(ncid)

            call stop_timer(field_io_timer)

        end subroutine write_netcdf_fields

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine write_field_double(id, fdata, start, cnt)
            integer,          intent(in) :: id
            double precision, intent(in) :: fdata(box%lo(3):box%hi(3), &
                                                  box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1))
            integer,          intent(in) :: cnt(4), start(4)

            if (nc_dset(id)%l_enabled) then
                call write_netcdf_dataset(ncid, nc_dset(id)%varid,      &
                                          fdata(box%lo(3):box%hi(3),    &
                                                box%lo(2):box%hi(2),    &
                                                box%lo(1):box%hi(1)),   &
                                          start, cnt)
            endif
        end subroutine write_field_double

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine read_netcdf_fields(ncfname, step)
            character(*),      intent(in) :: ncfname
            integer, optional, intent(in) :: step
            integer                       :: n_steps, start(4), cnt(4)

            call start_timer(field_io_timer)

            call set_netcdf_field_info

            call open_netcdf_file(ncfname, NF90_NOWRITE, ncid)

            call get_num_steps(ncid, n_steps)

            if (present(step)) then
                if (step == -1) then
                    ! do nothing
                    call mpi_print("Warning: Read the last time step.")
                else if (step > n_steps) then
                    call mpi_print(&
                        "Warning: NetCDF file has not enough records. The last step is chosen.")
                else
                    n_steps = step
                endif
            endif

            start(1:3) = box%lo + 1      ! need to add 1 since start must begin with index 1
            start(4) = n_steps
            cnt(1:3) = box%size
            cnt(4) = 1


            if (has_dataset(ncid, nc_dset(NC_X_VOR)%name)) then
                ! 19 May 2022
                ! https://stackoverflow.com/questions/45984672/print-values-without-new-line
                if (world%rank == world%root) then
                    write(*, "(a63)", advance="no") &
                        "Found " // nc_dset(NC_X_VOR)%name // " field input, reading ..."
                endif
                call read_netcdf_dataset(ncid,                      &
                                         nc_dset(NC_X_VOR)%name,    &
                                         vor(box%lo(3):box%hi(3),   &
                                             box%lo(2):box%hi(2),   &
                                             box%lo(1):box%hi(1),   &
                                             1),                    &
                                         start=start, cnt=cnt)
                if (world%rank == world%root) then
                    write(*, *) "done"
                endif
            endif

            if (has_dataset(ncid, nc_dset(NC_Y_VOR)%name)) then
                if (world%rank == world%root) then
                    write(*, "(a63)", advance="no") &
                        "Found " // nc_dset(NC_Y_VOR)%name // " field input, reading ..."
                endif
                call read_netcdf_dataset(ncid,                      &
                                         nc_dset(NC_Y_VOR)%name,    &
                                         vor(box%lo(3):box%hi(3),   &
                                             box%lo(2):box%hi(2),   &
                                             box%lo(1):box%hi(1),   &
                                             2),                    &
                                         start=start, cnt=cnt)
                if (world%rank == world%root) then
                    write(*, *) "done"
                endif
            endif

            if (has_dataset(ncid, nc_dset(NC_Z_VOR)%name)) then
                if (world%rank == world%root) then
                    write(*, "(a63)", advance="no") &
                        "Found " // nc_dset(NC_Z_VOR)%name // " field input, reading ..."
                endif
                call read_netcdf_dataset(ncid,                      &
                                         nc_dset(NC_Z_VOR)%name,    &
                                         vor(box%lo(3):box%hi(3),   &
                                             box%lo(2):box%hi(2),   &
                                             box%lo(1):box%hi(1),   &
                                             3),                    &
                                         start=start, cnt=cnt)
                if (world%rank == world%root) then
                    write(*, *) "done"
                endif
            endif

            if (has_dataset(ncid, nc_dset(NC_X_VEL)%name)) then
                if (world%rank == world%root) then
                    write(*, "(a63)", advance="no") &
                        "Found " // nc_dset(NC_X_VEL)%name // " field input, reading ..."
                endif
                call read_netcdf_dataset(ncid,                      &
                                         nc_dset(NC_X_VEL)%name,    &
                                         vel(box%lo(3):box%hi(3),   &
                                             box%lo(2):box%hi(2),   &
                                             box%lo(1):box%hi(1),   &
                                             1),                    &
                                         start=start, cnt=cnt)
                if (world%rank == world%root) then
                    write(*, *) "done"
                endif
            endif

            if (has_dataset(ncid, nc_dset(NC_Y_VEL)%name)) then
                if (world%rank == world%root) then
                    write(*, "(a63)", advance="no") &
                        "Found " // nc_dset(NC_Y_VEL)%name // " field input, reading ..."
                endif
                call read_netcdf_dataset(ncid,                      &
                                         nc_dset(NC_Y_VEL)%name,    &
                                         vel(box%lo(3):box%hi(3),   &
                                             box%lo(2):box%hi(2),   &
                                             box%lo(1):box%hi(1),   &
                                             2),                    &
                                         start=start, cnt=cnt)
                if (world%rank == world%root) then
                    write(*, *) "done"
                endif
            endif

            if (has_dataset(ncid, nc_dset(NC_Z_VEL)%name)) then
                if (world%rank == world%root) then
                    write(*, "(a63)", advance="no") &
                        "Found " // nc_dset(NC_Z_VEL)%name // " field input, reading ..."
                endif
                call read_netcdf_dataset(ncid,                      &
                                         nc_dset(NC_Z_VEL)%name,    &
                                         vel(box%lo(3):box%hi(3),   &
                                             box%lo(2):box%hi(2),   &
                                             box%lo(1):box%hi(1),   &
                                             3),                    &
                                         start=start, cnt=cnt)
                if (world%rank == world%root) then
                    write(*, *) "done"
                endif
            endif

#ifdef ENABLE_BUOYANCY
            if (has_dataset(ncid, nc_dset(NC_BUOY)%name)) then
                if (world%rank == world%root) then
                    write(*, "(a63)", advance="no") &
                        "Found " // nc_dset(NC_BUOY)%name // " field input, reading ..."
                endif
                call read_netcdf_dataset(ncid,                      &
                                         nc_dset(NC_BUOY)%name,     &
                                         buoy(box%lo(3):box%hi(3),  &
                                             box%lo(2):box%hi(2),   &
                                             box%lo(1):box%hi(1)),  &
                                         start=start, cnt=cnt)
                if (world%rank == world%root) then
                    write(*, *) "done"
                endif
            endif
#endif

            call close_netcdf_file(ncid)

            call stop_timer(field_io_timer)

        end subroutine read_netcdf_fields

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine set_netcdf_field_output
            integer :: n

            call set_netcdf_field_info

            ! check custom tags
            if (any('all' == output%field_list(:))) then
                nc_dset(:)%l_enabled = .true.
            else if (any('default' == output%field_list(:))) then
                nc_dset(NC_X_VOR)%l_enabled = .true.
                nc_dset(NC_Y_VOR)%l_enabled = .true.
                nc_dset(NC_Z_VOR)%l_enabled = .true.
                nc_dset(NC_X_VEL)%l_enabled = .true.
                nc_dset(NC_Y_VEL)%l_enabled = .true.
                nc_dset(NC_Z_VEL)%l_enabled = .true.
                nc_dset(NC_PRES)%l_enabled  = .true.
#ifdef ENABLE_BUOYANCY
                nc_dset(NC_BUOY)%l_enabled  = .true.
                nc_dset(NC_BUOY_AN)%l_enabled = .true.
#endif
            else
                ! check individual fields
                do n = 1, size(nc_dset)
                    nc_dset(n)%l_enabled = any(nc_dset(n)%name == output%field_list(:))
                enddo
            endif

            if (count(nc_dset(:)%l_enabled) == 0) then
                if (world%rank == world%root) then
                    print *, "WARNING: No fields are actively selected. PS3D is going to write"
                    print *, "         the default fields. Stop the simulation now if this is"
                    print *, "         not your intention. Fields can be provided to the list"
                    print *, "         'output%field_list' in the configuration file."
                    print *, "         The following fields are available:"
                    do n = 1, size(nc_dset)
                        print *, "         " // nc_dset(n)%name // " : " // trim(nc_dset(n)%long_name)
                    enddo
                    print *, "         " // "all"     // repeat(" ", 29) // " : write all fields"
                    print *, "         " // "default" // repeat(" ", 25) // " : write default fields"
                    print *, ""
                endif
                nc_dset(NC_X_VOR)%l_enabled = .true.
                nc_dset(NC_Y_VOR)%l_enabled = .true.
                nc_dset(NC_Z_VOR)%l_enabled = .true.
                nc_dset(NC_X_VEL)%l_enabled = .true.
                nc_dset(NC_Y_VEL)%l_enabled = .true.
                nc_dset(NC_Z_VEL)%l_enabled = .true.
                nc_dset(NC_PRES)%l_enabled  = .true.
#ifdef ENABLE_BUOYANCY
                nc_dset(NC_BUOY)%l_enabled  = .true.
                nc_dset(NC_BUOY_AN)%l_enabled = .true.
#endif
            endif

#ifdef ENABLE_VERBOSE
            if (verbose .and. (world%rank == world%root)) then
                print *, "PS3D is going to write the following fields:"
                do n = 1, size(nc_dset)
                    if (nc_dset(n)%l_enabled) then
                        print *, repeat(" ", 4) // trim(nc_dset(n)%name)
                    endif
                enddo
                print *, ""
            endif
#endif

        end subroutine set_netcdf_field_output

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine set_netcdf_field_info
            nc_dset(NC_X_VEL) = netcdf_field_info(name='x_velocity',                    &
                                                  long_name='x velocity component',     &
                                                  std_name='',                          &
                                                  unit='m/s',                           &
                                                  dtype=NF90_DOUBLE)

            nc_dset(NC_Y_VEL) = netcdf_field_info(name='y_velocity',                    &
                                                  long_name='y velocity component',     &
                                                  std_name='',                          &
                                                  unit='m/s',                           &
                                                  dtype=NF90_DOUBLE)

            nc_dset(NC_Z_VEL) = netcdf_field_info(name='z_velocity',                    &
                                                  long_name='z velocity component',     &
                                                  std_name='',                          &
                                                  unit='m/s',                           &
                                                  dtype=NF90_DOUBLE)

            nc_dset(NC_X_VOR) = netcdf_field_info(name='x_vorticity',                   &
                                                  long_name='x vorticity component',    &
                                                  std_name='',                          &
                                                  unit='1/s',                           &
                                                  dtype=NF90_DOUBLE)

            nc_dset(NC_Y_VOR) = netcdf_field_info(name='y_vorticity',                   &
                                                  long_name='y vorticity component',    &
                                                  std_name='',                          &
                                                  unit='1/s',                           &
                                                  dtype=NF90_DOUBLE)

            nc_dset(NC_Z_VOR) = netcdf_field_info(name='z_vorticity',                   &
                                                  long_name='z vorticity component',    &
                                                  std_name='',                          &
                                                  unit='1/s',                           &
                                                  dtype=NF90_DOUBLE)

#ifdef ENABLE_BUOYANCY
            nc_dset(NC_PRES) = netcdf_field_info(name='pressure_anomaly',               &
                                                 long_name='pressure anomaly',          &
#else
            nc_dset(NC_PRES) = netcdf_field_info(name='pressure',                       &
                                                 long_name='pressure',                  &
#endif
                                                 std_name='',                           &
                                                 unit='m^2/s^2',                        &
                                                 dtype=NF90_DOUBLE)

#ifdef ENABLE_BUOYANCY
            nc_dset(NC_BUOY) = netcdf_field_info(name='buoyancy',                       &
                                                 long_name='buoyancy',                  &
                                                 std_name='',                           &
                                                 unit='m/s^2',                          &
                                                 dtype=NF90_DOUBLE)

            nc_dset(NC_BUOY_AN) = netcdf_field_info(name='buoyancy_anomaly',            &
                                                    long_name='buoyancy anomaly',       &
                                                    std_name='',                        &
                                                    unit='m/s^2',                       &
                                                    dtype=NF90_DOUBLE)
#endif

        end subroutine set_netcdf_field_info

end module field_netcdf
