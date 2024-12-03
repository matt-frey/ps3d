! =============================================================================
!                               Test diffz
!
!  This unit test checks the calculation of the vertical differentiation of the
!  function
!       f(x, y, z) = z
! =============================================================================
program test_diffz_1
    use unit_test
    use constants, only : zero, one, two, pi, f12
    use parameters, only : lower, update_parameters, dx, nx, ny, nz, extent
    use mpi_timer
    use mpi_environment
    use mpi_layout
    use mpi_collectives
    use model, only : layout, create_model
    implicit none

    call mpi_env_initialise

    nx = 32
    ny = 32
    nz = 32

    lower  = (/zero, zero, zero/)
    extent = (/one, one, one/)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    call run_test("uniform")

    call run_test("chebyshev")

    call mpi_env_finalise

contains

    subroutine run_test(grid_type)
        character(*), intent(in)      :: grid_type
        double precision              :: error
        double precision, allocatable :: dfdz_ref(:, :, :), dfdz(:, :, :)
        double precision, allocatable :: fp(:, :, :)
        double precision, allocatable :: z(:)
        integer                       :: iz

        call create_model(grid_type, "Hou & Li")

        allocate(fp(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(dfdz(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(dfdz_ref(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(z(0:nz))

        dfdz = zero

        z = layout%get_z_axis()
        do iz = 0, nz
            fp(iz, :, :) = z(iz)
            dfdz_ref(iz, :, :) = one
        enddo

        ! calculate z-derivative (dfdz)
        call layout%diffz(fp, dfdz)

        error = maxval(dabs(dfdz_ref - dfdz))

        call mpi_blocking_reduce(error, MPI_MAX, world)

        if (world%rank == world%root) then
            call print_result_dp('Test diffz ' // grid_type, error, atol=3.0e-14)
        endif

        deallocate(fp, z)
        deallocate(dfdz)
        deallocate(dfdz_ref)

    end subroutine run_test

end program test_diffz_1
