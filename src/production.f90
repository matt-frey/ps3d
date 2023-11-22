! This program computes the enstrophy production rates
! according to
!   Vincent, A. & Meneguzzi, M. (1994).
!   The dynamics of vorticity tubes in homogeneous turbulence.
!   Journal of Fluid Mechanics, 258, 245-254.
!   doi:10.1017/S0022112094003319
program production
    use constants, only : zero, f12
    use field_netcdf, only : read_netcdf_fields
    use fields
    use netcdf_reader
    use parameters, only : lower, extent, nx, ny, nz, update_parameters, ncelli, dx
    use inversion_utils, only : init_inversion
    use sta3dfft, only : fftxyp2s, fftxys2p, diffx, diffy
    use fields, only : vor, svor, vel, svel
    use jacobi, only : jacobi_diagonalise
    use mpi_environment
    use mpi_layout
    implicit none

    character(512)                :: fname = ""
    character(512)                :: ncfname = ""
    integer                       :: step = 0, ncells(3), n_steps, ncid, i
    logical                       :: l_vertical_profile = .false.
    character(9)                  :: s_step
    double precision              :: t ! time
    double precision              :: z
    ! strain components
    double precision, allocatable :: dudx(:, :, :)  &
                                   , dudy(:, :, :)  &
                                   , dwdx(:, :, :)  &
                                   , dvdy(:, :, :)  &
                                   , dwdy(:, :, :)
    ! enstrophy production rates
    double precision              :: etas(3)

    call mpi_env_initialise

    if (world%size > 1) then
        call mpi_stop("This program does not run in parallel due to call_ptospc!")
    endif


    call parse_command_line

    call read_netcdf_domain(ncfname, lower, extent, ncells)
    nx = ncells(1)
    ny = ncells(2)
    nz = ncells(3)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    call init_inversion

    allocate(dudx(0:nz, 0:ny-1, 0:nx-1))
    allocate(dudy(0:nz, 0:ny-1, 0:nx-1))
    allocate(dwdx(0:nz, 0:ny-1, 0:nx-1))
    allocate(dvdy(0:nz, 0:ny-1, 0:nx-1))
    allocate(dwdy(0:nz, 0:ny-1, 0:nx-1))

    call field_default

    if (l_vertical_profile) then
        !
        ! Height profile of production rates at a specific time
        ! (averaged over x and y)
        !
        call open_netcdf_file(ncfname, NF90_NOWRITE, ncid)
        call get_time_at_step(ncid, step, t)
        call close_netcdf_file(ncid)

        print *, 'Calculate enstrophy production rates at time', t

        ! read field
        call read_netcdf_fields(ncfname, step)

        fname = ncfname(1:len(trim(ncfname))-3) // '_enstrophy_production_rates_step_' // trim(s_step) // '.asc'
        open(unit=1235, file=fname, status='replace')
        write(1235, *) '  # time', t
        write(1235, *) '  # height eta_1 eta_2 eta_3'

        ! calculate the strain components
        call calc_strain

        do i = 0, nz
            ! calculate enstrophy production rates at specific height
            ! (averaged over x and y)
            call calc_production_rate_at_height(i)

            ! write to file
            z = lower(3) + dble(i) * dx(3)
            write(1235, *) z, etas
        enddo
        close(1235)
    else
        !
        !   Time evolution of individual production rates
        !

        call open_netcdf_file(ncfname, NF90_NOWRITE, ncid)
        call get_num_steps(ncid, n_steps)
        call close_netcdf_file(ncid)

        ! 3 is subtracted to remove '.nc'
        fname = ncfname(1:len(trim(ncfname))-3) // '_enstrophy_production_rates.asc'
        open(unit=1235, file=fname, status='replace')
        write(1235, *) '  # time (s) eta_1 eta_2 eta_3'
        do step = 1, n_steps

            ! read field
            call read_netcdf_fields(ncfname, step)

            call open_netcdf_file(ncfname, NF90_NOWRITE, ncid)
            call get_time_at_step(ncid, step, t)
            call close_netcdf_file(ncid)

            ! calculate the strain components
            call calc_strain

            ! calculate enstrophy production rates
            call calc_producition_rates

            ! write to file
            write(1235, *) t, etas
        enddo
        close(1235)
    endif

    deallocate(dudx)
    deallocate(dudy)
    deallocate(dwdx)
    deallocate(dvdy)
    deallocate(dwdy)

    call mpi_env_finalise

    contains

        subroutine calc_production_rate_at_height(iz)
            integer, intent(in) :: iz
            integer             :: ix, iy, nc
            double precision    :: V(3, 3), D(3)

            etas(:) = zero
            D(:) = zero
            V(:, :) = zero

            do ix = 0, nx-1
                do iy = 0, ny-1
                    ! calculate eigenvalues and vectors of symmetrised strain matrix (1/2 * (S + S^T)
                    call calc_eigs(ix, iy, iz, D, V)

                    do nc = 1, 3
                        etas(nc) = etas(nc) + D(nc) * sum(vor(iz, iy, ix, :) * V(:, nc)) ** 2
                    enddo
                enddo
            enddo

            etas(:) = etas(:) / (nx * ny)

        end subroutine calc_production_rate_at_height

        subroutine calc_producition_rates
            integer          :: ix, iy, iz, nc
            double precision :: V(3, 3), D(3)

            etas(:) = zero

            do ix = 0, nx-1
                do iy = 0, ny-1

                    ! interior
                    do iz = 1, nz-1
                        ! calculate eigenvalues and vectors of symmetrised strain matrix (1/2 * (S + S^T)
                        call calc_eigs(ix, iy, iz, D, V)

                        do nc = 1, 3
                            etas(nc) = etas(nc) + D(nc) * sum(vor(iz, iy, ix, :) * V(:, nc)) ** 2
                        enddo
                    enddo

                    ! lower boundary
                    call calc_eigs(ix, iy, 0, D, V)
                    do nc = 1, 3
                        etas(nc) = etas(nc) + f12 * D(nc) * sum(vor(0, iy, ix, :) * V(:, nc)) ** 2
                    enddo

                    ! upper boundary
                    call calc_eigs(ix, iy, nz, D, V)
                    do nc = 1, 3
                        etas(nc) = etas(nc) + f12 * D(nc) * sum(vor(nz, iy, ix, :) * V(:, nc)) ** 2
                    enddo
                enddo
            enddo

            ! ncelli = 1 / (nx * ny * nz)
            etas = etas * ncelli

        end subroutine calc_producition_rates

        subroutine calc_strain
            integer          :: nc
            double precision :: ds(0:nz, 0:nx-1, 0:ny-1)

            ! -----------------------------------------------------------------
            ! Transform velocity components to semi-spectral space:
            ! (note: vel is overwritten and invalid after this operation)
            do nc = 1, 3
                call fftxyp2s(vel(:, :, :, nc), svel(:, :, :, nc))
            enddo

            ! -----------------------------------------------------------------
            ! Calculate velocity derivatives:

            ! du/dx
            call diffx(svel(:, :, :, 1), ds)
            call fftxys2p(ds, dudx)

            ! du/dy
            call diffy(svel(:, :, :, 1), ds)
            call fftxys2p(ds, dudy)

            ! dw/dx
            call diffx(svel(:, :, :, 3), ds)
            call fftxys2p(ds, dwdx)

            ! dv/dy
            call diffy(svel(:, :, :, 2), ds)
            call fftxys2p(ds, dvdy)

            ! dw/dy
            call diffy(svel(:, :, :, 3), ds)
            call fftxys2p(ds, dwdy)

        end subroutine calc_strain


        subroutine calc_eigs(ix, iy, iz, D, V)
            integer, intent(in)           :: ix, iy, iz
            double precision, intent(out) :: V(3, 3), D(3)
            double precision              :: strain(3, 3)


            ! get local symmetrised strain matrix, i.e. 1/ 2 * (S + S^T)
            ! where
            !     /u_x u_y u_z\
            ! S = |v_x v_y v_z|
            !     \w_x w_y w_z/
            ! with u_* = du/d* (also derivatives of v and w).
            ! The derivatives dv/dx, du/dz, dv/dz and dw/dz are calculated
            ! with vorticity or the assumption of incompressibility
            ! (du/dx + dv/dy + dw/dz = 0):
            !    dv/dx = \zeta + du/dy
            !    du/dz = \eta + dw/dx
            !    dv/dz = dw/dy - \xi
            !    dw/dz = - (du/dx + dv/dy)
            !
            !                         /  2 * u_x  u_y + v_x u_z + w_x\
            ! 1/2 * (S + S^T) = 1/2 * |u_y + v_x   2 * v_y  v_z + w_y|
            !                         \u_z + w_x  v_z + w_y   2 * w_z/
            !
            ! S11 = du/dx
            ! S12 = 1/2 * (du/dy + dv/dx) = 1/2 * (2 * du/dy + \zeta) = du/dy + 1/2 * \zeta
            ! S13 = 1/2 * (du/dz + dw/dx) = 1/2 * (\eta + 2 * dw/dx) = 1/2 * \eta + dw/dx
            ! S22 = dv/dy
            ! S23 = 1/2 * (dv/dz + dw/dy) = 1/2 * (2 * dw/dy - \xi) = dw/dy - 1/2 * \xi
            ! S33 = dw/dz = - (du/dx + dv/dy)
            strain(1, 1) = dudx(iz, iy, ix)                               ! S11
            strain(1, 2) = dudy(iz, iy, ix) + f12 * vor(iz, iy, ix, 3)    ! S12
            strain(1, 3) = dwdx(iz, iy, ix) + f12 * vor(iz, iy, ix, 2)    ! S13
            strain(2, 2) = dvdy(iz, iy, ix)                               ! S22
            strain(2, 3) = dwdy(iz, iy, ix) - f12 * vor(iz, iy, ix, 1)    ! S23
            strain(3, 3) = -(dudx(iz, iy, ix) + dvdy(iz, iy, ix))         ! S33

            call jacobi_diagonalise(strain, D, V)

        end subroutine calc_eigs


        subroutine parse_command_line
            integer                          :: i
            character(len=512)               :: arg

            i = 0
            do
                call get_command_argument(i, arg)
                if (len_trim(arg) == 0) then
                    exit
                endif

                if (arg == '--ncfile') then
                    i = i + 1
                    call get_command_argument(i, arg)
                    ncfname = trim(arg)
                else if (arg == '--step') then
                    i = i + 1
                    call get_command_argument(i, arg)
                    ! 21 Oct 2022
                    ! https://stackoverflow.com/questions/9900417/character-to-integer-conversion-in-fortran
                    read(arg,'(i9)') step
                    s_step = trim(arg)
                    l_vertical_profile = .true.
                else if (arg == '--help') then
                    print *, 'Run code with "production --ncfile [NetCDF file]"'
                    print *, '--step [step] (optional) calculate vertical profile of the'
                    print *, '                         enstrophy production rates at a'
                    print *, '                         specific time step'
                    stop
                endif
                i = i+1
            end do

            if (ncfname == '') then
                print *, 'No NetCDF file provided. Run code with "production --ncfile [NetCDF file]"'
                stop
            endif
        end subroutine parse_command_line
end program production
