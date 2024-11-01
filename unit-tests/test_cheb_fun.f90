! =============================================================================
!                               Test cheb_fun routine
! =============================================================================
program test_cheb_poly
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
    use cheby, only : cheb_fun
    implicit none

    double precision, allocatable :: f(:)
    double precision, allocatable :: coeff(:)
    integer                       :: iz
    double precision, allocatable :: sol(:)
    double precision              :: error

    call mpi_env_initialise

    call parse_command_line

    nx = 16
    ny = 16
    nz = 8

    lower  = (/0.0d0, 0.0d0, -1.0d0/)
    extent = (/1.0d0, 1.0d0, 1.0d0/)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    allocate(f(0:nz))
    allocate(coeff(0:nz))
    allocate(sol(0:nz))

    call update_parameters

    !-----------------------------------------------------------------
    ! Initialise inversion module:
    call init_inversion


    coeff = zero
    coeff(0) =  0.0d0
    coeff(1) =  8.0d0
    coeff(2) = -9.0d0
    coeff(3) =  3.0d0
    coeff(4) = -4.0d0
    coeff(5) =  5.0d0
    coeff(6) =  0.0d0
    coeff(7) =  2.5d0

    do iz = 0, nz
        sol(iz) =   8.0d0 * get_cheb_poly(zcheb(iz), 1) &
                  - 9.0d0 * get_cheb_poly(zcheb(iz), 2) &
                  + 3.0d0 * get_cheb_poly(zcheb(iz), 3) &
                  - 4.0d0 * get_cheb_poly(zcheb(iz), 4) &
                  + 5.0d0 * get_cheb_poly(zcheb(iz), 5) &
                  + 2.5d0 * get_cheb_poly(zcheb(iz), 9)
    enddo

    call cheb_fun(nz, coeff, f)

    do iz = 0, nz
        print *, iz, sol(iz), f(iz)
    enddo

    error = maxval(abs(sol - f))

    call mpi_blocking_reduce(error, MPI_MAX, world)

    if (world%rank == world%root) then
        call print_result_dp('Test cheb_fun', error, atol=9.0d-15)
    endif

    call finalise_inversion

    deallocate(f, coeff, sol)

    call mpi_env_finalise

    contains

        function get_cheb_poly(z, n) result(res)
            double precision, intent(in) :: z
            integer,          intent(in) :: n
            double precision             :: res

            res = cos(dble(n) * acos(z))
        end function get_cheb_poly

end program test_cheb_poly
