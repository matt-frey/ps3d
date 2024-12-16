! =============================================================================
!                               Test diffz
!
!  This unit test checks the calculation of the vertical differentiation of the
!  function
!       f(x, y, z) = 1 - 3*z^2
! =============================================================================
program test_diffz_4
    use unit_test
    use constants, only : zero, one, three, six
    use parameters, only : lower, update_parameters, nx, ny, nz, extent
    use mpi_timer
    use mpi_environment
    use mpi_layout
    use mpi_collectives
    use model, only : layout, create_model
    use sta3dfft, only : fftxyp2s, fftxys2p
    implicit none

    call mpi_env_initialise

    nx = 32
    ny = nx
    nz = nx

    lower  = (/zero, zero, zero/)
    extent = (/one, one, one/)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    call run_test("uniform", atol=1.0d-1)

    call run_test("chebyshev", atol=2.0d-13)

    call mpi_env_finalise

contains

    subroutine run_test(grid_type, atol)
        character(*),     intent(in)  :: grid_type
        double precision, intent(in)  :: atol
        double precision              :: error
        double precision, allocatable :: dfdz_ref(:, :, :), dfdz(:, :, :)
        double precision, allocatable :: fp(:, :, :), fs(:, :, :)
        double precision, allocatable :: z(:)
        integer                       :: iz

        call create_model(grid_type, "Hou & Li")

        allocate(fp(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(dfdz(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(dfdz_ref(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(z(0:nz))

        dfdz = zero

        z = layout%get_z_axis()
        do iz = 0, nz
            fp(iz, :, :) = one - three * z(iz) ** 2
            dfdz_ref(iz, :, :) = - six * z(iz)
        enddo

        call fftxyp2s(fp, fs)

        call layout%diffz(fs, dfdz, l_decomposed=.false.)

        fs = dfdz

        call fftxys2p(fs, dfdz)

        error = maxval(abs(dfdz_ref - dfdz))

        call mpi_blocking_reduce(error, MPI_MAX, world)

        if (world%rank == world%root) then
            call print_result_dp('Test diffz 4 ' // grid_type, error, atol=atol)
        endif

        deallocate(fp, fs, z)
        deallocate(dfdz)
        deallocate(dfdz_ref)

    end subroutine

end program test_diffz_4
