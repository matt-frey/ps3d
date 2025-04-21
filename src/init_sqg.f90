! ====================================================================================
!               Initialize SQG Balanced Flow from given surface bouyancy field:
!
! In DGD formulation: Parameters are: N, f, H. K = sqrt(k^2 + ky^2)
!
!       b(kx, ky, z) = sinh( sigma K (z + H))/sinh( sigma K H) b_0(kx,ky); K > 0
!       b(0 , 0 , z) = (z + H)/H b_0(0, 0);
!
!       psi(kx, ky, z) = (N K)^(-1) cosh( sigma K (z + H))/sinh( sigma K H) b_0(kx,ky); K >0
!       psi( 0, 0 , z) = (1/f) z b_0(0, 0);
!
! in the domain (x,y) in (-pi,pi); z in (-H,0);
!
!       From this:
!                  1) ox(kx,ky,z) = (-1/f) i kx b(kx, ky, z)    (2.3)
!                  2) oy(kx,ky,z) = (-1/f) i ky b(kx, ky, z)    (2.3)
!                  3) oz(kx,ky,z) = -K^2 psi(kx, ky, z)         (2.2)
!
! Note:  use the following identities to compute hyperbolic functions:
!        sinh(a+b)/sinh(a) = exp(b) ( 1 - exp(-2(a+b)) )/(1 - exp(-2a))
!        cosh(a+b)/sinh(a) = exp(b) ( 1 + exp(-2(a+b)) )/(1 - exp(-2a))
! ====================================================================================
program init_sqg
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
                      , write_physical_quantities, f_cor, bfsq
    use model, only : layout, create_model
    use mpi_utils, only : mpi_stop
    use sta3dfft, only : initialise_fft &
                       , k2l2           &
                       , xfactors       &
                       , xtrig          &
                       , yfactors       &
                       , ytrig          &
                       , fftxys2p       &
                       , diffx          &
                       , diffy
    use sta2dfft, only : ptospc
    implicit none

    logical                     :: verbose = .false.
    character(len=512)          :: filename = ''
    character(len=512)          :: b0fname = ''
    character(len=512)          :: ncfname = ''
    integer                     :: ncid
    integer                     :: dimids(4), axids(4)
    integer                     :: buoy_id
    integer                     :: x_vor_id, y_vor_id, z_vor_id

    type mesh_type
        integer          :: ncells(3)   ! number of cells
        double precision :: extent(3)   ! size of domain
        double precision :: origin(3)   ! origin of domain (lower left corner)
        character(len=9) :: layout      ! "uniform" or "chebyshev"
    end type mesh_type

    type(mesh_type) :: grid

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

    call initialise_fft(extent)

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
                                    name='buoyancy_anomaly',            &
                                    long_name='buoyancy_anomaly',       &
                                    std_name='',                        &
                                    unit='m/s^2',                       &
                                    dtype=NF90_DOUBLE,                  &
                                    dimids=dimids,                      &
                                    varid=buoy_id)

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

        call sqg_project

        ! write time
        call write_netcdf_scalar(ncid, axids(4), zero, 1)

        call close_netcdf_file(ncid)
    end subroutine generate_fields

    subroutine sqg_project
        double precision, allocatable :: x(:), y(:), z(:)
        double precision              :: sbuoy(0:nz,                 &
                                               box%lo(2):box%hi(2),  &
                                               box%lo(1):box%hi(1))
        double precision              :: psi(0:nz,                   &
                                               box%lo(2):box%hi(2),  &
                                               box%lo(1):box%hi(1))
        double precision              :: buoy(0:nz,                 &
                                              box%lo(2):box%hi(2),  &
                                              box%lo(1):box%hi(1))
        double precision              :: svor(0:nz,                     &
                                             box%lo(2):box%hi(2),       &
                                             box%lo(1):box%hi(1), 3)
        double precision              :: vor(0:nz,                     &
                                             box%lo(2):box%hi(2),      &
                                             box%lo(1):box%hi(1), 3)
        integer                       :: iz, kx, ky
        integer                       :: cnt(4), start(4)
        double precision              :: t, f, bf, sigma, nkinv, kl, a, b, expb, expm2a, expm2ab, ss, cs
        double precision              :: b0(0:ny-1, 0:nx-1)
        double precision              :: wkc(0:ny-1, 0:nx-1)
        double precision              :: sb0(0:ny-1, 0:nx-1)
        integer                       :: nhbytes
        double precision              :: H

        if ((f_cor(1) /= zero) .and. (f_cor(2) /= zero)) then
            call mpi_stop("Expecting only vertical Coriolis frequency.")
        endif


        f = f_cor(3)
        bf = sqrt(bfsq)

        allocate(x(0:nx-1), y(0:ny-1), z(0:nz))
        x = layout%get_x_axis()
        y = layout%get_y_axis()
        z = layout%get_z_axis()

        H = grid%extent(3)

        if (verbose .and. (world%rank == world%root)) then
            print *, "Coriolis frequency:", f
            print *, "buoyancy frequency:", bf
            print *, "Domain depth:", H
            print *, "Domain origin:", grid%origin
            print *, "Domain extent:", grid%extent
        endif



 !!!!!   What this should do:
 !!      either take in the surface field b0, or read it here
 !!      Given b0:  b0 -> fft2(b0)

 !!      Given H, N, f: sig = N/f;
 !!      Make 2D array: AK = sqrt(k2l2);
 !!      Make 2D array: SK = sig*sqrt(k2l2);
 !!!     Loop over z and get b',ox,oy,oz (x,y,z) as per initial comments

        !Read dimensionless surface bouyancy, b0:
        ! (note: each MPI rank reads the whole field)
        nhbytes = 8*(nx*ny+1)
        open(11, file=trim(b0fname), form='unformatted', &
             access='direct', status='old', recl=nhbytes)
        read(11, rec=1) t, b0
        close(11)

        ! 2D FFT of b0 --> sb0
        wkc = b0
        call ptospc(nx, ny, wkc, sb0, xfactors, yfactors, xtrig, ytrig)

        if (verbose .and. (world%rank == world%root)) then
            print *, "Surface buoyancy anomaly: min = ", minval(b0), "max = ", maxval(b0)
        endif

        sigma = bf / f
        ! need to go from serial to parallel storage
        do kx = box%lo(1), box%hi(1)
            do ky = box%lo(2), box%hi(2)
                if ((kx == 0) .and. (ky == 0)) then
                    do iz = 0, nz
                        sbuoy(iz, 0, 0) = (z(iz) + H) / H * sb0(0, 0)

                        psi(iz, 0, 0) = (one / f) * z(iz) * sb0(0, 0)
                    enddo
                else
                    kl = sqrt(k2l2(ky, kx))

                    nkinv = one / (bf * kl)

                    a = sigma * kl * H
                    expm2a = one / (one - exp(-two * a))

                    do iz = 0, nz
                        b = sigma * kl * z(iz)

                        expm2ab = exp(-two * (a + b))
                        expb = exp(b)

                        ! sinh(a+b)/sinh(a) = exp(b) ( 1 - exp(-2(a+b)) )/(1 - exp(-2a))
                        ss = expb * (one - expm2ab) * expm2a

                        ! cosh(a+b)/sinh(a) = exp(b) ( 1 + exp(-2(a+b)) )/(1 - exp(-2a))
                        cs = expb * (one + expm2ab) * expm2a

                        sbuoy(iz, ky, kx) = ss * sb0(ky, kx)
                        psi(iz, ky, kx) = nkinv * cs * sb0(ky, kx)

                        svor(iz, ky, kx, 3) = - k2l2(ky, kx) * psi(iz, ky, kx)

                    enddo
                endif
            enddo
        enddo

        ! ox(kx,ky,z) = (-1/f) i kx b(kx, ky, z)
        call diffx(sbuoy, svor(:, :, :, 1))

        svor(:, :, :, 1) = -one / f * svor(:, :, :, 1)

        ! oy(kx,ky,z) = (-1/f) i ky b(kx, ky, z)
        call diffy(sbuoy, svor(:, :, :, 2))

        svor(:, :, :, 2) = -one / f * svor(:, :, :, 2)

        call fftxys2p(sbuoy, buoy)

        call fftxys2p(svor(:, :, :, 1), vor(:, :, :, 1))
        call fftxys2p(svor(:, :, :, 2), vor(:, :, :, 2))
        call fftxys2p(svor(:, :, :, 3), vor(:, :, :, 3))

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


    end subroutine sqg_project

    ! parse configuration file
    ! (see https://cyber.dabamos.de/programming/modernfortran/namelists.html [8 March 2021])
    subroutine read_config_file
        integer :: ios
        integer :: fn = 1
        logical :: l_exist = .false.

        ! namelist definitions
        namelist /MODELS/ ncfname, grid, b0fname

        ! check whether file l_exist
        inquire(file=filename, exist=l_exist)

        if (l_exist .eqv. .false.) then
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
        inquire(file=ncfname, exist=l_exist)

        if (l_exist) then
            call mpi_stop('Error: output file "' // trim(ncfname) // '" already exists.')
        endif


        inquire(file=b0fname, exist=l_exist)

        if (.not. l_exist) then
            call mpi_stop('Error: Surface buoyancy "' // trim(b0fname) // '" does not exist.')
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
                call mpi_stop('Run code with "init_sqg --config [config file]"')
            endif
            i = i+1
        end do

        if (filename == '') then
            call mpi_stop('No configuration file provided. Run code with "init_sqg --config [config file]"')
        endif
    end subroutine parse_command_line
end program init_sqg
