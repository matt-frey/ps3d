! =============================================================================
!                               Test filter
!
!  This unit test checks the filter using the Beltrami vorticity field
!    xi(x, y, z) = alpha * u(x, y, z)
!   eta(x, y, z) = alpha * v(x, y, z)
!  zeta(x, y, z) = alpha * w(x, y, z)
! and adds some spatial noise to it.
! =============================================================================
program test_filter
    use unit_test
    use constants, only : one, two, pi, f12
    use parameters, only : lower, update_parameters, nx, ny, nz, extent
    use fields
    use mpi_timer
    use mpi_environment
    use mpi_layout
    use mpi_collectives, only : mpi_blocking_reduce
    use model, only : layout, create_model
    implicit none

    character(len=9) :: grid_types(2)
    character(len=8) :: filter_types(2)
    integer          :: i, j

    call mpi_env_initialise

    call parse_command_line

    grid_types = (/"uniform  ", "chebyshev"/)
    filter_types = (/"Hou & Li", "2/3-rule"/)

    nx = 32
    ny = 32
    nz = 128

    lower  = -f12 * pi * (/one, one, one/)
    extent =  pi * (/one, one, one/)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    call field_default

    do i = 1, size(grid_types)
        do j = 1, size(filter_types)
            call run_test(grid_types(i), filter_types(j))
        enddo
    enddo

    call mpi_env_finalise

contains

    subroutine run_test(grid_type, filter_type)
        character(*), intent(in)      :: grid_type
        character(*), intent(in)      :: filter_type
        double precision              :: error
        integer                       :: ix, iy, iz
        double precision              :: alpha, fk2l2, k, l, m, u, v, w
        double precision              :: cosmz, sinmz, sinkxly, coskxly
        double precision, allocatable :: x(:), y(:), z(:), unfiltered(:, :, :, :)
        double precision              :: mean, stddev, noise


        !--------------------------------------------------
        ! Set up initial vorticity field:
        k = two
        l = two
        m = one

        mean = 0.0d0
        stddev = 0.025d0

        allocate(unfiltered(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1), 3))

        call create_model(grid_type, filter_type)

        alpha = dsqrt(k ** 2 + l ** 2 + m ** 2)
        fk2l2 = one / dble(k ** 2 + l ** 2)

        allocate(x(0:nx-1), y(0:ny-1), z(0:nz))

        x = layout%get_x_axis()
        y = layout%get_y_axis()
        z = layout%get_z_axis()
        do ix = box%lo(1), box%hi(1)
            do iy = box%lo(2), box%hi(2)
                do iz = 0, nz
                    cosmz = cos(m * z(iz))
                    sinmz = sin(m * z(iz))
                    sinkxly = sin(k * x(ix) + l * y(iy))
                    coskxly = cos(k * x(ix) + l * y(iy))

                    ! velocity
                    u = fk2l2 * (k * m * sinmz - l * alpha * cosmz) * sinkxly
                    v = fk2l2 * (l * m * sinmz + k * alpha * cosmz) * sinkxly
                    w = cosmz * coskxly

                    call random_gaussian(mean, stddev, noise)

                    ! vorticity
                    vor(iz, iy, ix, 1) = alpha * u + noise
                    vor(iz, iy, ix, 2) = alpha * v + noise
                    vor(iz, iy, ix, 3) = alpha * w + noise

                enddo
            enddo
        enddo


        call layout%decompose_physical(vor(:, :, :, 1), svor(:, :, :, 1))
        call layout%decompose_physical(vor(:, :, :, 2), svor(:, :, :, 2))
        call layout%decompose_physical(vor(:, :, :, 3), svor(:, :, :, 3))

        unfiltered = vor

        call layout%apply_filter(svor(:, :, :, 1))
        call layout%apply_filter(svor(:, :, :, 2))
        call layout%apply_filter(svor(:, :, :, 3))

        call layout%combine_physical(svor(:, :, :, 1), vor(:, :, :, 1))
        call layout%combine_physical(svor(:, :, :, 2), vor(:, :, :, 2))
        call layout%combine_physical(svor(:, :, :, 3), vor(:, :, :, 3))

        if (verbose .and. (world%rank == world%root)) then
            print *, ""
            print *, "FILTERED FIELD"
            print *, "z                  f(x)                  g(x)                       ||g(x)-f(x)|| "
            do iz = 0, nz
                print *, z(iz), unfiltered(iz, 1, 1, 1), vor(iz, 1, 1, 1), &
                                abs(unfiltered(iz, 1, 1, 1) - vor(iz, 1, 1, 1))
            enddo
        endif

        deallocate(x, y, z)
        deallocate(unfiltered)

        call mpi_blocking_reduce(error, MPI_MAX, world)

        if (world%rank == world%root) then
            call print_result_logical('Test filter ' // filter_type // " and " // trim(grid_type), .true.)
        endif
    end subroutine

    subroutine random_gaussian(mean, stddev, noise)
        double precision, intent(in)  :: mean       ! Mean of the Gaussian noise
        double precision, intent(in)  :: stddev     ! Standard deviation of the Gaussian noise
        double precision, intent(out) :: noise      ! Output noise value
        double precision              :: u1, u2, z0

        ! Generate two uniform random numbers in the range (0, 1)
        call random_number(u1)
        call random_number(u2)

        ! Apply the Box-Muller transform to generate a Gaussian random variable
        z0 = sqrt(-2.0d0 * log(u1)) * cos(2.0d0 * 3.141592653589793d0 * u2)

        ! Scale by the standard deviation and add the mean
        noise = z0 * stddev + mean

    end subroutine random_gaussian

end program test_filter
