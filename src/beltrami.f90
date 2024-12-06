! ====================================================================================
!                           3D BELTRAMI FLOW
!
! Beltrami flow:
!     u(x, y, z) = (k^2 + l^2)^(-1) * [k*m*sin(mz) - l*alpha*cos(m*z) * sin(k*x + l*y)]
!     v(x, y, z) = (k^2 + l^2)^(-1) * [l*m*sin(mz) + k*alpha*cos(m*z) * sin(k*x + l*y)]
!     w(x, y, z) = cos(m*z) * cos(k*x + l*y)
! The vorticity of this flow is
!    xi(x, y, z) = alpha * u(x, y, z) + xi_pert
!   eta(x, y, z) = alpha * v(x, y, z) + eta_pert
!  zeta(x, y, z) = alpha * w(x, y, z) + zeta_pert
!
! where
!    xi_pert = a * cos(2y) * cos(z)
!   eta_pert = b * cos(2x) * cos(z)
!  zeta_pert = 0.0
!
!                   The code uses following variable names:
!                       beltrami_flow%k = k
!                       beltrami_flow%l = l
!                       beltrami_flow%m = m
!                       beltrami_flow%a = a
!                       beltrami_flow%b = b
! ====================================================================================
program beltrami
    use constants, only : pi, f12, zero, one, two
    use parameters, only : nx, ny, nz, lower, extent    &
                         , write_netcdf_parameters      &
                         , update_parameters, grid_type
    use netcdf_utils
    use netcdf_writer
    use mpi_environment
    use mpi_layout
    use config, only : package_version, cf_version
    use physics, only : read_physical_quantities_from_namelist
    use model, only : layout, create_model
    use mpi_utils, only : mpi_stop
    implicit none

    logical                     :: verbose = .false.
    character(len=512)          :: filename = ''
    character(len=512)          :: ncfname = ''
    integer                     :: ncid
    integer                     :: dimids(4), axids(4)
    double precision            :: fk2l2, kk, ll, mm, alpha, aa, bb
    integer                     :: x_vor_id, y_vor_id, z_vor_id
    double precision, parameter :: hpi = f12 * pi

    type mesh_type
        integer          :: ncells(3)   ! number of cells
        double precision :: extent(3)   ! size of domain
        double precision :: origin(3)   ! origin of domain (lower left corner)
        character(len=9) :: layout      ! "uniform" or "chebyshev"
    end type mesh_type

    type(mesh_type) :: grid

    type beltrami_type
        integer          :: k, l, m
        double precision :: a, b
    end type beltrami_type

    type(beltrami_type) :: beltrami_flow


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

        call define_netcdf_spatial_dimensions_3d(ncid=ncid,             &
                                                 ngps=(/nx, ny, nz+1/), &
                                                 dimids=dimids(1:3),    &
                                                 axids=axids(1:3))

        call define_netcdf_temporal_dimension(ncid, dimids(4), axids(4))

        call define_netcdf_dataset(ncid=ncid,                           &
                                    name='x_vorticity',                 &
                                    long_name='x vorticity component',  &
                                    std_name='',                        &
                                    unit='1/s',                         &
                                    dtype=NF90_DOUBLE,                  &
                                    dimids=dimids,                      &
                                    varid=x_vor_id)

        call define_netcdf_dataset(ncid=ncid,                           &
                                    name='y_vorticity',                 &
                                    long_name='y vorticity component',  &
                                    std_name='',                        &
                                    unit='1/s',                         &
                                    dtype=NF90_DOUBLE,                  &
                                    dimids=dimids,                      &
                                    varid=y_vor_id)

        call define_netcdf_dataset(ncid=ncid,                           &
                                    name='z_vorticity',                 &
                                    long_name='z vorticity component',  &
                                    std_name='',                        &
                                    unit='1/s',                         &
                                    dtype=NF90_DOUBLE,                  &
                                    dimids=dimids,                      &
                                    varid=z_vor_id)

        call close_definition(ncid)

        call close_netcdf_file(ncid)

        call beltrami_init

        ! write time
        call write_netcdf_scalar(ncid, axids(4), zero, 1)

        call close_netcdf_file(ncid)
    end subroutine generate_fields

    subroutine beltrami_init
        double precision, allocatable :: x(:), y(:), z(:)
        double precision              :: vor(0:nz,                      &
                                             box%lo(2):box%hi(2),       &
                                             box%lo(1):box%hi(1), 3)
        integer                       :: i, j, k
        integer                       :: cnt(4), start(4)


        kk = dble(beltrami_flow%k)
        ll = dble(beltrami_flow%l)
        mm = dble(beltrami_flow%m)
        aa = dble(beltrami_flow%a)
        bb = dble(beltrami_flow%b)
        alpha = sqrt(kk ** 2 + ll ** 2 + mm ** 2)
        fk2l2 = alpha / dble(beltrami_flow%k ** 2 + beltrami_flow%l ** 2)


        allocate(x(0:nx-1), y(0:ny-1), z(0:nz))
        x = layout%get_x_axis()
        y = layout%get_y_axis()
        z = layout%get_z_axis()
        do i = box%lo(1), box%hi(1)
            do j = box%lo(2), box%hi(2)
                do k = 0, nz
                    vor(k, j, i, :) = get_flow_vorticity(x(i), y(j), z(k))
                enddo
            enddo
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

        call write_netcdf_dataset(ncid, x_vor_id, vor(0:nz,                    &
                                                      box%lo(2):box%hi(2),     &
                                                      box%lo(1):box%hi(1), 1), &
                                                      start, cnt)

        call write_netcdf_dataset(ncid, y_vor_id, vor(0:nz,                    &
                                                      box%lo(2):box%hi(2),     &
                                                      box%lo(1):box%hi(1), 2), &
                                                      start, cnt)

        call write_netcdf_dataset(ncid, z_vor_id, vor(0:nz,                    &
                                                      box%lo(2):box%hi(2),     &
                                                      box%lo(1):box%hi(1), 3), &
                                                      start, cnt)

    end subroutine beltrami_init

    function get_flow_vorticity(x, y, z) result(omega)
        double precision, intent(in) :: x, y, z
        double precision             :: omega(3)
        double precision             :: cosmz, sinmz, sinkxly, coskxly

        cosmz = cos(mm * z)
        sinmz = sin(mm * z)
        sinkxly = sin(kk * x + ll * y)
        coskxly = cos(kk * x + ll * y)

        omega(1) = fk2l2 * (kk * mm * sinmz - ll * alpha * cosmz) * sinkxly
        omega(2) = fk2l2 * (ll * mm * sinmz + kk * alpha * cosmz) * sinkxly
        omega(3) = alpha * cosmz * coskxly

        ! add horizontal perturbation
        omega(1) = omega(1) + aa * cos(2.0d0 * y) * cos(z)
        omega(2) = omega(2) + bb * cos(2.0d0 * x) * cos(z)

    end function get_flow_vorticity


    ! parse configuration file
    ! (see https://cyber.dabamos.de/programming/modernfortran/namelists.html [8 March 2021])
    subroutine read_config_file
        integer :: ios
        integer :: fn = 1
        logical :: exists = .false.

        ! namelist definitions
        namelist /MODELS/ ncfname, grid, beltrami_flow

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
                call mpi_stop('Run code with "beltrami --config [config file]"')
            endif
            i = i+1
        end do

        if (filename == '') then
            call mpi_stop('No configuration file provided. Run code with "beltrami --config [config file]"')
        endif
    end subroutine parse_command_line
end program beltrami
