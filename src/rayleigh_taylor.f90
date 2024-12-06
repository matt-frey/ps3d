! ====================================================================================
!                           3D Rayleigh-Taylor flow
!
! Rayleigh-Taylor flow with buoyancy
!
!       b(x, y, z) = = -sin(z) + eps * cos^2(z) * h(x,y)
!
! and horizontal perturbation
!
!       h(x, y) = cos(4x) * cos(2y + pi/6) + sin(2x + pi/6) * sin(4y)
!
! in the domain [-pi/2, pi/2]^3.
!
! ====================================================================================
program rayleigh_taylor
    use constants, only : zero, f12, fpi6, one, two, pi
    use parameters, only : nx, ny, nz, lower, extent    &
                         , write_netcdf_parameters      &
                         , update_parameters, grid_type
    use netcdf_utils
    use netcdf_writer
    use mpi_environment
    use mpi_layout
    use config, only : package_version, cf_version
    use physics, only : read_physical_quantities_from_namelist &
                      , write_physical_quantities
    use model, only : layout, create_model
    use mpi_utils, only : mpi_stop
    implicit none

    logical                     :: verbose = .false.
    character(len=512)          :: filename = ''
    character(len=512)          :: ncfname = ''
    integer                     :: ncid
    integer                     :: dimids(4), axids(4)
    integer                     :: buoy_id

    type mesh_type
        integer          :: ncells(3)   ! number of cells
        double precision :: extent(3)   ! size of domain
        double precision :: origin(3)   ! origin of domain (lower left corner)
        character(len=9) :: layout      ! "uniform" or "chebyshev"
    end type mesh_type

    type(mesh_type) :: grid

    double precision :: eps

    call mpi_env_initialise

    ! Read command line (verbose, filename, etc.)
    call parse_command_line

    call read_config_file

    nx = grid%ncells(1)
    ny = grid%ncells(2)
    nz = grid%ncells(3)

    ! make origin and extent always a multiple of pi
    grid%origin = pi * grid%origin
    grid%extent = pi * grid%extent

    ! write box
    lower = grid%origin
    extent = grid%extent
    grid_type = grid%layout

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    call read_physical_quantities_from_namelist(trim(filename))

    ! Filter is being ignored here
    call create_model(grid%layout, "Hou & Li")

    call generate_fields

    call mpi_env_finalise

contains

    subroutine generate_fields

        call create_netcdf_file(ncfname, .false., ncid)


        ! define global attributes
        call write_netcdf_info(ncid=ncid,                     &
                                version_tag=package_version,  &
                                file_type='fields',           &
                                cf_version=cf_version)

        call write_netcdf_parameters(ncid)

        call write_physical_quantities(ncid)

        call define_netcdf_spatial_dimensions_3d(ncid=ncid,             &
                                                 ngps=(/nx, ny, nz+1/), &
                                                 dimids=dimids(1:3),    &
                                                 axids=axids(1:3))

        call define_netcdf_temporal_dimension(ncid, dimids(4), axids(4))

        call define_netcdf_dataset(ncid=ncid,                           &
                                    name='buoyancy',                    &
                                    long_name='buoyancy',               &
                                    std_name='',                        &
                                    unit='m/s^2',                       &
                                    dtype=NF90_DOUBLE,                  &
                                    dimids=dimids,                      &
                                    varid=buoy_id)

        call close_definition(ncid)

        call close_netcdf_file(ncid)

        call rt_init

        ! write time
        call write_netcdf_scalar(ncid, axids(4), zero, 1)

        call close_netcdf_file(ncid)
    end subroutine generate_fields

    subroutine rt_init
        double precision, allocatable :: x(:), y(:), z(:)
        double precision              :: buoy(0:nz,                 &
                                              box%lo(2):box%hi(2),  &
                                              box%lo(1):box%hi(1))
        integer                       :: i, j, k
        integer                       :: cnt(4), start(4)
        double precision              :: h(box%lo(2):box%hi(2), &
                                           box%lo(1):box%hi(1))

        allocate(x(0:nx-1), y(0:ny-1), z(0:nz))
        x = layout%get_x_axis()
        y = layout%get_y_axis()
        z = layout%get_z_axis()
        do i = box%lo(1), box%hi(1)
            do j = box%lo(2), box%hi(2)
                h(j, i) = cos(4.0d0 * x(i)) * cos(2.0d0 * y(j) + fpi6) &
                        + sin(2.0d0 * x(i) + fpi6) * sin(4.0d0 * y(j))
            enddo
        enddo

        do k = 0, nz
            buoy(k, :, :) = -sin(z(k)) + eps * cos(z(k)) ** 2 * h
        enddo

        call open_netcdf_file(ncfname, NF90_WRITE, ncid)

        call write_netcdf_axis(ncid, dimids(1), x)
        call write_netcdf_axis(ncid, dimids(2), y)
        call write_netcdf_axis(ncid, dimids(3), z)

        deallocate(x, y, z)

        ! time step to write [step(4) is the time]
        ! need to add 1 since start must begin with index 1
        start(1:3) = box%lo + 1
        start(4) = 1

        cnt(1:3) = box%hi - box%lo + 1
        cnt(4)   = 1

        call write_netcdf_dataset(ncid, buoy_id, buoy(0:nz,                 &
                                                      box%lo(2):box%hi(2),  &
                                                      box%lo(1):box%hi(1)), &
                                                      start, cnt)
    end subroutine rt_init

    ! parse configuration file
    ! (see https://cyber.dabamos.de/programming/modernfortran/namelists.html [8 March 2021])
    subroutine read_config_file
        integer :: ios
        integer :: fn = 1
        logical :: exists = .false.

        ! namelist definitions
        namelist /MODELS/ ncfname, grid, eps

        ! check whether file exists
        inquire(file=filename, exist=exists)

        if (exists .eqv. .false.) then
            call mpi_stop('Error: input file "' // trim(filename) // '" does not exist.')
            stop
        endif

        ! open and read Namelist file.
        open(action='read', file=filename, iostat=ios, newunit=fn)

        read(nml=MODELS, iostat=ios, unit=fn)

        if (ios /= 0) then
            call mpi_stop('Error: invalid Namelist format.')
        end if

        close(fn)

        ! check whether NetCDF file already exists
        inquire(file=ncfname, exist=exists)

        if (exists) then
            call mpi_stop('Error: output file "' // trim(ncfname) // '" already exists.')
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
                call mpi_stop('Run code with "rayleigh_taylor --config [config file]"')
            endif
            i = i+1
        end do

        if (filename == '') then
            call mpi_stop('No configuration file provided. Run code with "rayleigh_taylor --config [config file]"')
        endif
    end subroutine parse_command_line
end program rayleigh_taylor
