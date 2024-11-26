! ====================================================================================
!                           3D BELTRAMI FLOW
!
! Beltrami flow:
!     u(x, y, z) = (k^2 + l^2)^(-1) * [k*m*sin(mz) - l*alpha*cos(m*z) * sin(k*x + l*y)]
!     v(x, y, z) = (k^2 + l^2)^(-1) * [l*m*sin(mz) + k*alpha*cos(m*z) * sin(k*x + l*y)]
!     w(x, y, z) = cos(m*z) * cos(k*x + l*y)
! The vorticity of this flow is
!    xi(x, y, z) = alpha * u(x, y, z)
!   eta(x, y, z) = alpha * v(x, y, z)
!  zeta(x, y, z) = alpha * w(x, y, z)
!
!                   The code uses following variable names:
!                       beltrami_flow%k = k
!                       beltrami_flow%l = l
!                       beltrami_flow%m = m
! ====================================================================================
program beltrami
    use constants, only : pi, f12, zero, one, two
    use parameters, only : nx, ny, nz, dx, lower, extent
    use netcdf_utils
    use netcdf_writer
    use config, only : package_version, cf_version
    use physics, only : read_physical_quantities_from_namelist
    use model_factory, only : ops
    implicit none

    logical                     :: verbose = .false.
    character(len=512)          :: filename = ''
    character(len=512)          :: ncfname = ''
    integer                     :: ncid
    integer                     :: dimids(4), axids(4)
    double precision            :: fk2l2, kk, ll, mm, alpha
    integer                     :: x_vor_id, y_vor_id, z_vor_id
    double precision, parameter :: hpi = f12 * pi

    type box_type
        integer          :: ncells(3)   ! number of cells
        double precision :: extent(3)   ! size of domain
        double precision :: origin(3)   ! origin of domain (lower left corner)
    end type box_type

    type(box_type) :: box

    type beltrami_type
        integer          :: k, l, m
    end type beltrami_type

    type(beltrami_type) :: beltrami_flow


    ! Read command line (verbose, filename, etc.)
    call parse_command_line

    call read_config_file

    call read_physical_quantities_from_namelist(trim(filename))

    call generate_fields

    contains

        subroutine generate_fields

            call create_netcdf_file(ncfname, .false., ncid, l_serial=.true.)

            dx = box%extent / dble(box%ncells)
            nx = box%ncells(1)
            ny = box%ncells(2)
            nz = box%ncells(3)

            ! define global attributes
            call write_netcdf_info(ncid=ncid,                    &
                                   version_tag=package_version,  &
                                   file_type='fields',           &
                                   cf_version=cf_version)

            call write_netcdf_box(ncid, lower, extent, box%ncells)

            call define_netcdf_spatial_dimensions_3d(ncid=ncid,             &
                                                     ngps=(/nx, ny, nz+1/), &
                                                     dimids=dimids(1:3),    &
                                                     axids=axids(1:3))

            call define_netcdf_temporal_dimension(ncid, dimids(4), axids(4))

            ! make origin and extent always a multiple of pi
            box%origin = pi * box%origin
            box%extent = pi * box%extent
            dx = dx * pi

            ! write box
            lower = box%origin
            extent = box%extent
            call write_netcdf_box(ncid, lower, extent, box%ncells)

            call beltrami_init

            call write_netcdf_axis(ncid, dimids(1), ops%get_x_axis())
            call write_netcdf_axis(ncid, dimids(2), ops%get_y_axis())
            call write_netcdf_axis(ncid, dimids(3), ops%get_z_axis())

            ! write time
            call write_netcdf_scalar(ncid, axids(4), zero, 1)

            call close_netcdf_file(ncid, l_serial=.true.)
        end subroutine generate_fields

        subroutine beltrami_init
            double precision                :: pos(3)
            double precision                :: vor(0:nz, 0:ny-1, 0:nx-1, 3)
            integer                         :: i, j, k

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

            call close_definition(ncid)


            kk = dble(beltrami_flow%k)
            ll = dble(beltrami_flow%l)
            mm = dble(beltrami_flow%m)
            alpha = dsqrt(kk ** 2 + ll ** 2 + mm ** 2)
            fk2l2 = alpha / dble(beltrami_flow%k ** 2 + beltrami_flow%l ** 2)

            do i = 0, nx - 1
                do j = 0, ny - 1
                    do k = 0, nz
                        pos = box%origin + dx * dble((/i, j, k/))
                        vor(k, j, i, :) = get_flow_vorticity(pos)
                    enddo
                enddo
            enddo

            call write_netcdf_dataset(ncid, x_vor_id, vor(:, :, :, 1))
            call write_netcdf_dataset(ncid, y_vor_id, vor(:, :, :, 2))
            call write_netcdf_dataset(ncid, z_vor_id, vor(:, :, :, 3))

        end subroutine beltrami_init

        function get_flow_vorticity(pos) result(omega)
            double precision, intent(in) :: pos(3)
            double precision             :: x, y, z
            double precision             :: omega(3)
            double precision             :: cosmz, sinmz, sinkxly, coskxly

            x = pos(1)
            y = pos(2)
            z = pos(3)

            cosmz = dcos(mm * z)
            sinmz = dsin(mm * z)
            sinkxly = dsin(kk * x + ll * y)
            coskxly = dcos(kk * x + ll * y)

            omega(1) = fk2l2 * (kk * mm * sinmz - ll * alpha * cosmz) * sinkxly
            omega(2) = fk2l2 * (ll * mm * sinmz + kk * alpha * cosmz) * sinkxly
            omega(3) = alpha * cosmz * coskxly

        end function get_flow_vorticity


        ! parse configuration file
        ! (see https://cyber.dabamos.de/programming/modernfortran/namelists.html [8 March 2021])
        subroutine read_config_file
            integer :: ios
            integer :: fn = 1
            logical :: exists = .false.

            ! namelist definitions
            namelist /MODELS/ ncfname, box, beltrami_flow

            ! check whether file exists
            inquire(file=filename, exist=exists)

            if (exists .eqv. .false.) then
                print *, 'Error: input file "', trim(filename), '" does not exist.'
                stop
            endif

            ! open and read Namelist file.
            open(action='read', file=filename, iostat=ios, newunit=fn)

            read(nml=MODELS, iostat=ios, unit=fn)

            if (ios /= 0) then
                print *, 'Error: invalid Namelist format.'
                stop
            end if

            close(fn)

            ! check whether NetCDF file already exists
            inquire(file=ncfname, exist=exists)

            if (exists) then
                print *, 'Error: output file "', trim(ncfname), '" already exists.'
                stop
            endif
        end subroutine read_config_file

        ! Get the file name provided via the command line
        subroutine parse_command_line
            integer                          :: i
            character(len=512)               :: arg

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
                else if (arg == '--verbose') then
                    verbose = .true.
                else if (arg == '--help') then
                    print *, 'Run code with "beltrami --config [config file]"'
                    stop
                endif
                i = i+1
            end do

            if (filename == '') then
                print *, 'No configuration file provided. Run code with "beltrami --config [config file]"'
                stop
            endif
        end subroutine parse_command_line
end program beltrami
