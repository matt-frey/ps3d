! =============================================================================
!                               Test diffz
!
!  This unit test checks the calculation of the vertical differentiation of the
!  function
!       f(x, y, z) = z * cos(k*x) * sin(l*y)
! =============================================================================
program test_diffz_2
    use unit_test
    use constants, only : zero, one, two, pi, twopi, f12
    use parameters, only : lower, update_parameters, dx, nx, ny, nz, extent
    use mpi_timer
    use mpi_environment
    use mpi_layout
    use mpi_collectives
    use model, only : layout, create_model
    implicit none


    call mpi_env_initialise

    nx = 16
    ny = 32
    nz = 32

    lower  = (/zero, -f12*pi, zero/)
    extent = (/pi, twopi, twopi/)

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
        double precision, allocatable :: x(:), y(:), z(:)
        integer                       :: ix, iy, iz
        double precision              :: k, l, coskx, sinly

        call create_model(grid_type, "Hou & Li")

        allocate(fp(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(dfdz(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(dfdz_ref(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(x(0:nx-1), y(0:ny-1), z(0:nz))

        dfdz = zero

        k = two
        l = one

        x = layout%get_x_axis()
        y = layout%get_y_axis()
        z = layout%get_z_axis()
        do ix = box%lo(1), box%hi(1)
            coskx = cos(k * x(ix))
            do iy = box%lo(2), box%hi(2)
                sinly = sin(l * y(iy))
                do iz = 0, nz
                    fp(iz, iy, ix) = z(iz) * coskx * sinly
                    dfdz_ref(iz, iy, ix) = coskx * sinly
                enddo
            enddo
        enddo

        call layout%diffz(fp, dfdz)

        error = maxval(dabs(dfdz_ref - dfdz))

        call mpi_blocking_reduce(error, MPI_MAX, world)

        if (world%rank == world%root) then
            call print_result_dp('Test diffz 2 ' // grid_type, error, atol=3.0e-13)
        endif

        deallocate(x, y, z)
        deallocate(fp)
        deallocate(dfdz)
        deallocate(dfdz_ref)

    end subroutine run_test

end program test_diffz_2
