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
    use parameters, only : lower, update_parameters, nx, ny, nz, extent
    use fields
    use inversion_mod, only : vor2vel, vor2vel_timer
    use mpi_timer
    use mpi_environment
    use mpi_layout
    use mpi_collectives, only : mpi_blocking_reduce
    use model, only : layout, create_model
    implicit none

    call mpi_env_initialise

    call register_timer('vorticity', vor2vel_timer)

    nx = 32
    ny = nx
    nz = nx

    lower  = (/-f12, -f12, zero/)
    extent = (/one, one, one/)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    call field_default

    call run_test("uniform")

    call run_test("chebyshev")

    call mpi_env_finalise

contains

    subroutine run_test(grid_type)
        character(*), intent(in)      :: grid_type
        double precision              :: error
        double precision, allocatable :: vel_ref(:, :, :, :)
        double precision, allocatable :: z(:)
        integer                       :: iz

        call create_model(grid_type, "Hou & Li")

        allocate(vel_ref(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1), 3))
        allocate(z(0:nz))

        z = layout%get_z_axis()
        do iz = 0, nz
            ! velocity
            vel_ref(iz, :, :, 1) = z(iz) ** 2 - f13
            vel_ref(iz, :, :, 2) = zero
            vel_ref(iz, :, :, 3) = zero

            ! vorticity
            vor(iz, :, :, 1) = zero
            vor(iz, :, :, 2) = two * z(iz)
            vor(iz, :, :, 3) = zero
        enddo

        call layout%decompose_physical(vor(:, :, :, 1), svor(:, :, :, 1))
        call layout%decompose_physical(vor(:, :, :, 2), svor(:, :, :, 2))
        call layout%decompose_physical(vor(:, :, :, 3), svor(:, :, :, 3))

        call vor2vel

        error = maxval(abs(vel_ref - vel))

        print *, error

        call mpi_blocking_reduce(error, MPI_MAX, world)

        if (world%rank == world%root) then
            call print_result_dp('Test vor2vel', error, atol=1.3d-15)
        endif

        deallocate(z)
        deallocate(vel_ref)

    end subroutine run_test

end program test_vor2vel_4
