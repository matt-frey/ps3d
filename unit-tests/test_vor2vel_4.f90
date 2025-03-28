! =============================================================================
!                               Test vor2vel
!
!  This unit test checks the calculation of the velocity field using the
!  following flow:
!     u(x, y, z) = z^2 - 1/3
!     v(x, y, z) = 0
!     w(x, y, z) = 0
!  and vorticity
!    xi(x, y, z) = 0
!   eta(x, y, z) = 2 * z
!  zeta(x, y, z) = 0
! with x, y in [-1/2, 1/2] and z in [0, 1].
! =============================================================================
program test_vor2vel_4
    use unit_test
    use constants, only : one, f12, f13, two
    use parameters, only : lower, update_parameters, dx, nx, ny, nz, extent
    use fields
    use inversion_utils
    use inversion_mod, only : vor2vel, vor2vel_timer
    use mpi_timer
    use mpi_environment
    use mpi_layout
    use mpi_collectives, only : mpi_blocking_reduce
    implicit none

    double precision              :: error
    double precision, allocatable :: vel_ref(:, :, :, :)
    integer                       :: iz
    double precision              :: z

    call mpi_env_initialise

    call register_timer('vorticity', vor2vel_timer)

    nx = 32
    ny = nx
    nz = nx

    lower  = (/-f12, -f12, zero/)
    extent = (/one, one, one/)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    allocate(vel_ref(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1), 3))

    call update_parameters

    call field_default

    call init_inversion

    do iz = 0, nz
        z = lower(3) + iz * dx(3)

        ! velocity
        vel_ref(iz, :, :, 1) = z ** 2 - f13
        vel_ref(iz, :, :, 2) = zero
        vel_ref(iz, :, :, 3) = zero

        ! vorticity
        vor(iz, :, :, 1) = zero
        vor(iz, :, :, 2) = two * z
        vor(iz, :, :, 3) = zero
    enddo

    call field_decompose_physical(vor(:, :, :, 1), svor(:, :, :, 1))
    call field_decompose_physical(vor(:, :, :, 2), svor(:, :, :, 2))
    call field_decompose_physical(vor(:, :, :, 3), svor(:, :, :, 3))

    call vor2vel

    error = maxval(dabs(vel_ref - vel))

    call mpi_blocking_reduce(error, MPI_MAX, world)

    if (world%rank == world%root) then
        call print_result_dp('Test vor2vel', error, atol=1.0d-15)
    endif

    deallocate(vel_ref)

    call mpi_env_finalise

end program test_vor2vel_4
