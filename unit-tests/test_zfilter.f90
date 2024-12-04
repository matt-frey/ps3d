! =============================================================================
!                               Test vertical filter module
! =============================================================================
program test_zfilter
    use unit_test
    use constants, only : zero, f12, one, pi, twopi, f14
    use parameters, only : lower, nx, ny, nz, extent &
                         , update_parameters
    use fields
    use inversion_utils
    use mpi_timer
    use mpi_environment
    use mpi_layout
    use mpi_collectives, only : mpi_blocking_reduce
    use options, only : zfiltering
    implicit none

    double precision, allocatable :: f(:, :, :)
    double precision, allocatable :: g(:, :, :)
    integer                       :: iz

    call mpi_env_initialise

    call parse_command_line

    nx = 16
    ny = 16
    nz = 128

    lower  = (/-pi, -pi, -1.0d0/)
    extent = (/twopi, twopi, 2.0d0/)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    allocate(f(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
    allocate(g(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

    zfiltering%kmax = 1.0d0 / 8.0d0

    call update_parameters

    !-----------------------------------------------------------------
    ! Initialise inversion module:
    call init_inversion


    do iz = 0, nz
        f(iz, :, :) = myf(zcheb(iz), 0.025d0)
    enddo

    g = f

    call apply_zfilter(g)


!     call mpi_blocking_reduce(, MPI_MAX, world)

    if (verbose .and. (world%rank == world%root)) then
        print *, ""
        print *, "FILTERED FIELD"
        print *, "            x                  f(x)                  gf(x)                       ||g(x)-f(x)|| "
        do iz = 0, nz
            print *, zcheb(iz),  f(iz, 0, 0), g(iz, 0, 0), dabs(f(iz, 0, 0) - g(iz, 0, 0))
        enddo
    endif

    if (world%rank == world%root) then
        call print_result_dp('Test zfilter', 0.0d0, atol=6.0e-7)
    endif

    call finalise_inversion

    deallocate(f)
    deallocate(g)

    call mpi_env_finalise

    contains

        function myf(z, stddev) result(res)
            double precision, intent(in) :: z      ! Input variable z
            double precision, intent(in) :: stddev
            double precision             :: res    ! Result of the function
            double precision             :: noise  ! Gaussian noise

            ! Generate Gaussian noise with mean 0 and standard deviation 'stddev'
            call random_gaussian(0.0d0, stddev, noise)

            res = abs(cos(8.0d0 * pi * z)) * (1.0d0 - z **2) * exp(2.0d0 * z)

        end function myf

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

end program test_zfilter
