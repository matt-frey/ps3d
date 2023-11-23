! =============================================================================
!                     Test horizontal spectral differentiation
!
!                   This unit test checks the diff2dx soubroutine.
! =============================================================================
program test_diff2dx
    use unit_test
    use constants, only : pi, twopi, f12, zero, four, two
    use sta3dfft, only : initialise_fft, finalise_fft, diff2dx, fft2d, ifft2d
    use parameters, only : update_parameters, dx, nx, ny, nz, lower, extent
    use mpi_environment
    use mpi_layout
    implicit none

    double precision, allocatable :: fp(:, :, :), &
                                     fs(:, :, :), &
                                     ds(:, :, :)
    integer                       :: i, j
    double precision              :: x, y
    logical                       :: passed = .false.

    call mpi_env_initialise

    passed = (world%err == 0)

    nx = 64
    ny = 128
    nz = 64
    lower = (/-pi, -pi, -pi/)
    extent = (/twopi, twopi, twopi/)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    allocate(fp(1, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
    allocate(fs(1, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
    allocate(ds(1, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

    call initialise_fft(extent)

    fp(:, :, :) = zero

    ! setup test field
    do i = box%lo(1), box%hi(1)
        x = lower(1) + dble(i) * dx(1)
        do j = box%lo(2), box%hi(2)
            y = lower(2) + dble(j) * dx(2)
            fp(1, j, i) = dcos(four * x)
        enddo
    enddo

    fs = zero

    ! forward FFT
    call fft2d(fp, fs)

    fp = zero

    call diff2dx(fs, ds)

    ! inverse FFT
    call ifft2d(ds, fp)

    ! check result test field
    do i = box%lo(1), box%hi(1)
        x = lower(1) + dble(i) * dx(1)
        do j = box%lo(2), box%hi(2)
            y = lower(2) + dble(j) * dx(2)
            passed = (passed .and. (fp(1, j, i) - (-four * dsin(four * x)) < 1.0e-12))
        enddo
    enddo

    call finalise_fft

    if (world%rank == world%root) then
        call MPI_Reduce(MPI_IN_PLACE, passed, 1, MPI_LOGICAL, MPI_LAND, world%root, world%comm, world%err)
    else
        call MPI_Reduce(passed, passed, 1, MPI_LOGICAL, MPI_LAND, world%root, world%comm, world%err)
    endif

    call mpi_env_finalise

    passed = (passed .and. (world%err == 0))

    if (world%rank == world%root) then
        call print_result_logical('Test MPI diffx', passed)
    endif

end program test_diff2dx
