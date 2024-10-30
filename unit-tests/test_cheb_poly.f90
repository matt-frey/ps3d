! =============================================================================
!                               Test cheb_poly routine
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
    use zops, only : cheb_poly
    implicit none

    double precision, allocatable :: f(:, :, :)
    double precision, allocatable :: coeff(:, :, :)
    integer                       :: iz

    call mpi_env_initialise

    call parse_command_line

    nx = 16
    ny = 16
    nz = 64

    lower  = (/-pi, -pi, -1.0d0/)
    extent = (/twopi, twopi, 2.0d0/)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    allocate(f(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
    allocate(coeff(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

    call update_parameters

    !-----------------------------------------------------------------
    ! Initialise inversion module:
    call init_inversion


    do iz = 0, nz
        f(iz, :, :) =         get_cheb_poly(zcheb(iz), 1) &
                    - 2.0d0 * get_cheb_poly(zcheb(iz), 2) &
                    + 3.0d0 * get_cheb_poly(zcheb(iz), 3) &
                    - 4.0d0 * get_cheb_poly(zcheb(iz), 4) &
                    + 5.0d0 * get_cheb_poly(zcheb(iz), 5)
    enddo

    call cheb_poly(f, coeff)

    do iz = 0, 8
        print *, coeff(iz, 0, 0)
    enddo

!     call mpi_blocking_reduce(, MPI_MAX, world)


    if (world%rank == world%root) then
        call print_result_dp('Test cheb_poly', 0.0d0, atol=6.0e-7)
    endif

    call finalise_inversion

    deallocate(f)
    deallocate(coeff)

    call mpi_env_finalise

    contains

        function get_cheb_poly(z, n) result(res)
            double precision, intent(in) :: z
            integer,          intent(in) :: n
            double precision             :: res

            res = cos(dble(n) * acos(z))
        end function get_cheb_poly

end program test_cheb_poly
