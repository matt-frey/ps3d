module field_netcdf
    use constants, only : one
    use netcdf_utils
    use netcdf_writer
    use netcdf_reader
    use fields
    use config, only : package_version, cf_version
    use timer, only : start_timer, stop_timer
    use options, only : write_netcdf_options
    use physics, only : write_physical_quantities
    use parameters, only : lower, extent, dx, nx, ny, nz
    use inversion_utils, only : fftxys2p
    implicit none

    integer :: field_io_timer

    character(len=512) :: ncfname
    integer            :: ncid
    integer            :: dimids(4)     ! = (x, y, z)
    integer            :: coord_ids(3)  ! = (x, y, z)
    integer            :: t_axis_id

    integer            :: x_vel_id, y_vel_id, z_vel_id, &
                          x_vor_id, y_vor_id, z_vor_id, &
                          buoy_id, n_writes, diss_id    &
                          xvtend_id, yvtend_id, zvtend_id

    private :: ncid, ncfname,                   &
               dimids,                          &
               coord_ids, t_axis_id,            &
               x_vel_id, y_vel_id, z_vel_id,    &
               x_vor_id, y_vor_id, z_vor_id,    &
               buoy_id, n_writes, diss_id       &
               xvtend_id, yvtend_id, zvtend_id

    contains

        ! Create the field file.
        ! @param[in] basename of the file
        ! @param[in] overwrite the file
        subroutine create_netcdf_field_file(basename, overwrite)
            character(*), intent(in)  :: basename
            logical,      intent(in)  :: overwrite

            ncfname =  basename // '_fields.nc'

            n_writes = 1

            call create_netcdf_file(ncfname, overwrite, ncid)

            ! define global attributes
            call write_netcdf_info(ncid=ncid,                    &
                                   ps3d_version=package_version, &
                                   file_type='fields',           &
                                   cf_version=cf_version)

            call write_netcdf_box(ncid, lower, extent, (/nx, ny, nz/))

            call write_physical_quantities(ncid)

            call write_netcdf_options(ncid)

            ! define dimensions
            call define_netcdf_spatial_dimensions_3d(ncid=ncid,                &
                                                     ncells=(/nx, ny, nz/),    &
                                                     dimids=dimids(1:3),       &
                                                     axids=coord_ids)


            call define_netcdf_temporal_dimension(ncid, dimids(4), t_axis_id)

            ! define fields
            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='x_velocity',                   &
                                       long_name='x velocity component',    &
                                       std_name='',                         &
                                       unit='m/s',                          &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=x_vel_id)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='y_velocity',                   &
                                       long_name='y velocity component',    &
                                       std_name='',                         &
                                       unit='m/s',                          &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=y_vel_id)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='z_velocity',                   &
                                       long_name='z velocity component',    &
                                       std_name='',                         &
                                       unit='m/s',                          &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=z_vel_id)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='x_vorticity',                  &
                                       long_name='x vorticity component',   &
                                       std_name='',                         &
                                       unit='1/s',                          &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=x_vor_id)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='y_vorticity',                  &
                                       long_name='y vorticity component',   &
                                       std_name='',                         &
                                       unit='1/s',                          &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=y_vor_id)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='z_vorticity',                  &
                                       long_name='z vorticity component',   &
                                       std_name='',                         &
                                       unit='1/s',                          &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=z_vor_id)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='buoyancy',                     &
                                       long_name='total buoyancy',          &
                                       std_name='',                         &
                                       unit='m/s^2',                        &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=buoy_id)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='dissipation',                  &
                                       long_name='dissipation operator',    &
                                       std_name='',                         &
                                       unit='1',                            &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=diss_id)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='x_vorticity_tendency',         &
                                       long_name='x vort. tend. component', &
                                       std_name='',                         &
                                       unit='',                             &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=xvtend_id)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='y_vorticity_tendency',         &
                                       long_name='y vort. tend. component', &
                                       std_name='',                         &
                                       unit='',                             &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=yvtend_id)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='z_vorticity_tendency',         &
                                       long_name='z vort. tend. component', &
                                       std_name='',                         &
                                       unit='',                             &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=zvtend_id)

            call close_definition(ncid)

        end subroutine create_netcdf_field_file

        ! Pre-condition: Assumes an open file
        subroutine read_netcdf_field_content

            call get_dim_id(ncid, 'x', dimids(1))

            call get_dim_id(ncid, 'y', dimids(2))

            call get_dim_id(ncid, 'z', dimids(3))

            call get_dim_id(ncid, 't', dimids(4))


            call get_var_id(ncid, 'x', coord_ids(1))

            call get_var_id(ncid, 'y', coord_ids(2))

            call get_var_id(ncid, 'z', coord_ids(3))

            call get_var_id(ncid, 't', t_axis_id)

            call get_var_id(ncid, 'x_velocity', x_vel_id)

            call get_var_id(ncid, 'y_velocity', y_vel_id)

            call get_var_id(ncid, 'z_velocity', z_vel_id)

            call get_var_id(ncid, 'x_vorticity', x_vor_id)

            call get_var_id(ncid, 'y_vorticity', y_vor_id)

            call get_var_id(ncid, 'z_vorticity', z_vor_id)

            call get_var_id(ncid, 'buoyancy', buoy_id)

            call get_var_id(ncid, 'dissipation', diss_id)

            call get_var_id(ncid, 'x_vorticity_tendency', xvtend_id)

            call get_var_id(ncid, 'y_vorticity_tendency', yvtend_id)

            call get_var_id(ncid, 'z_vorticity_tendency', zvtend_id)

        end subroutine read_netcdf_field_content

        ! Write a step in the field file.
        ! @param[in] t is the time
        ! @param[in] dt is the time step
        subroutine write_netcdf_fields(t)
            double precision, intent(in) :: t
            integer                      :: cnt(4), start(4)
            double precision             :: bs(0:nz, 0:nx-1, 0:ny-1) ! buoyancy in semi-spectral space (temporary)
            double precision             :: vtend(0:nz, 0:ny-1, 0:nx-1)

            call start_timer(field_io_timer)

            call open_netcdf_file(ncfname, NF90_WRITE, ncid)

            if (n_writes == 1) then
                call write_netcdf_axis_3d(ncid, dimids(1:3), lower, dx, (/nx, ny, nz/))
            endif

            ! write time
            call write_netcdf_scalar(ncid, t_axis_id, t, n_writes)

            ! time step to write [step(4) is the time]
            cnt   = (/ nx, ny, nz+1, 1        /)
            start = (/ 1,  1,  1,    n_writes /)

            !
            ! write fields (do not write halo cells)
            !
            call write_netcdf_dataset(ncid, x_vel_id, vel(0:nz, 0:ny-1, 0:nx-1, 1), &
                                      start, cnt)
            call write_netcdf_dataset(ncid, y_vel_id, vel(0:nz, 0:ny-1, 0:nx-1, 2), &
                                      start, cnt)
            call write_netcdf_dataset(ncid, z_vel_id, vel(0:nz, 0:ny-1, 0:nx-1, 3), &
                                      start, cnt)

            call write_netcdf_dataset(ncid, x_vor_id, vor(0:nz, 0:ny-1, 0:nx-1, 1), &
                                      start, cnt)
            call write_netcdf_dataset(ncid, y_vor_id, vor(0:nz, 0:ny-1, 0:nx-1, 2), &
                                      start, cnt)
            call write_netcdf_dataset(ncid, z_vor_id, vor(0:nz, 0:ny-1, 0:nx-1, 3), &
                                      start, cnt)

            bs = sbuoy
            call fftxys2p(bs, buoy)
            call write_netcdf_dataset(ncid, buoy_id, buoy(0:nz, 0:ny-1, 0:nx-1),    &
                 start, cnt)

            call write_netcdf_dataset(ncid, diss_id, diss(0:nz, 0:ny-1, 0:nx-1),    &
                                      start, cnt)

            bs = svtend(:, :, :, 1)
            call fftxys2p(bs, vtend)
            call write_netcdf_dataset(ncid, xvtend_id, vtend(0:nz, 0:ny-1, 0:nx-1), &
                 start, cnt)

            bs = svtend(:, :, :, 2)
            call fftxys2p(bs, vtend)
            call write_netcdf_dataset(ncid, yvtend_id, vtend(0:nz, 0:ny-1, 0:nx-1), &
                 start, cnt)

            bs = svtend(:, :, :, 3)
            call fftxys2p(bs, vtend)
            call write_netcdf_dataset(ncid, zvtend_id, vtend(0:nz, 0:ny-1, 0:nx-1), &
                                      start, cnt)

            ! increment counter
            n_writes = n_writes + 1

            call close_netcdf_file(ncid)

            call stop_timer(field_io_timer)

        end subroutine write_netcdf_fields


        subroutine read_netcdf_fields(ncfname, step)
            character(*),      intent(in) :: ncfname
            integer, optional, intent(in) :: step
            integer                       :: n_steps, start(4), cnt(4)

            call start_timer(field_io_timer)

            call open_netcdf_file(ncfname, NF90_NOWRITE, ncid)

            call get_num_steps(ncid, n_steps)

            if (present(step)) then
                if (step > n_steps) then
                    print *, "Warning: NetCDF file has not enough records. The last step is chosen."
                else
                    n_steps = step
                endif
            endif

            cnt  =  (/ nx, ny, nz+1, 1       /)
            start = (/ 1,  1,  1,    n_steps /)

            if (has_dataset(ncid, 'x_vorticity')) then
                ! 19 May 2022
                ! https://stackoverflow.com/questions/45984672/print-values-without-new-line
                write(*, "(a42)", advance="no") "Found x-vorticity field input, reading ..."
                call read_netcdf_dataset(ncid, 'x_vorticity',            &
                                         vor(0:nz, 0:ny-1, 0:nx-1, 1),   &
                                         start=start, cnt=cnt)
                write(*, *) "done"
            endif

            if (has_dataset(ncid, 'y_vorticity')) then
                write(*, "(a42)", advance="no") "Found y-vorticity field input, reading ..."
                call read_netcdf_dataset(ncid, 'y_vorticity',            &
                                         vor(0:nz, 0:ny-1, 0:nx-1, 2),   &
                                         start=start, cnt=cnt)
                 write(*, *) "done"
            endif

            if (has_dataset(ncid, 'z_vorticity')) then
                write(*, "(a42)", advance="no") "Found z-vorticity field input, reading ..."
                call read_netcdf_dataset(ncid, 'z_vorticity',            &
                                         vor(0:nz, 0:ny-1, 0:nx-1, 3),   &
                                         start=start, cnt=cnt)
                 write(*, *) "done"
            endif

            if (has_dataset(ncid, 'buoyancy')) then
                write(*, "(a39)", advance="no") "Found buoyancy field input, reading ..."
                call read_netcdf_dataset(ncid, 'buoyancy',            &
                                         buoy(0:nz, 0:ny-1, 0:nx-1),  &
                                         start=start, cnt=cnt)
                 write(*, *) "done"
            endif

            call close_netcdf_file(ncid)

            call stop_timer(field_io_timer)

        end subroutine read_netcdf_fields

end module field_netcdf
