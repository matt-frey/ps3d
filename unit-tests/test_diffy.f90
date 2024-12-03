! =============================================================================
!                     Test horizontal spectral differentiation
!
!                   This unit test checks the diffy subroutine.
! =============================================================================
program test_diffy
    use unit_test
    use constants, only : pi, twopi, f12, zero, four, two
    use sta3dfft, only : fftxyp2s, fftxys2p, diffy
    use parameters, only : update_parameters, dx, nx, ny, nz, lower, extent
    use mpi_environment
    use mpi_layout
    use model, only : layout, create_model
    implicit none

    logical :: passed = .false.

    call mpi_env_initialise

    passed = (world%err == 0)

    nx = 128
    ny = 64
    nz = 64
    lower = (/-pi, -pi, -pi/)
    extent = (/twopi, twopi, twopi/)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    call run_test("uniform")

    call run_test("chebyshev")

    call mpi_env_finalise


contains

    subroutine run_test(grid_type)
        character(*), intent(in)      :: grid_type
        double precision, allocatable :: fp(:, :, :), &
                                         fs(:, :, :), &
                                         ds(:, :, :)
        double precision, allocatable :: y(:)
        integer                       :: i, j, k

        call create_model(grid_type, "Hou & Li")

        allocate(fp(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(fs(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(ds(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(y(0:ny-1))

        fp(:, :, :) = zero

        y = layout%get_y_axis()

        ! setup test field
        do i = box%lo(1), box%hi(1)
            do j = box%lo(2), box%hi(2)
                do k = box%lo(3), box%hi(3)
                    fp(k, j, i) = cos(four * y(j))
                enddo
            enddo
        enddo

        fs = zero

        ! forward FFT
        call fftxyp2s(fp, fs)

        fp = zero

        call diffy(fs, ds)

        ! inverse FFT
        call fftxys2p(ds, fp)

        ! check result test field
        do i = box%lo(1), box%hi(1)
            do j = box%lo(2), box%hi(2)
                do k = box%lo(3), box%hi(3)
                    passed = (passed .and. (fp(k, j, i) - (-four * sin(four * y(j))) < 1.0e-12))
                enddo
            enddo
        enddo

        deallocate(fp, fs, ds, y)

        if (world%rank == world%root) then
            call MPI_Reduce(MPI_IN_PLACE, passed, 1, MPI_LOGICAL, MPI_LAND, world%root, world%comm, world%err)
        else
            call MPI_Reduce(passed, passed, 1, MPI_LOGICAL, MPI_LAND, world%root, world%comm, world%err)
        endif

        passed = (passed .and. (world%err == 0))

        if (world%rank == world%root) then
            call print_result_logical('Test MPI diffy ' // grid_type, passed)
        endif

    end subroutine run_test

end program test_diffy
