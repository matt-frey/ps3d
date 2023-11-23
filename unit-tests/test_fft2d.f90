! =============================================================================
!                       Test FFT module
! =============================================================================
program test_fft2d
    use unit_test
    use constants, only : pi, twopi, f12, zero, four, two
    use sta3dfft, only : initialise_fft, finalise_fft, fft2d, ifft2d, fftxyp2s
    use parameters, only : update_parameters, dx, nx, ny, nz, lower, extent
    use mpi_environment
    use mpi_layout
    implicit none

    double precision              :: error = zero
    double precision, allocatable :: fp1(:, :, :),      &
                                     fp2(:, :, :),      &
                                     fs(:, :, :),       &
                                     fp_ref(:, :, :),   &
                                     fs_ref(:, :, :)
    integer                       :: i, j
    double precision              :: x, y
    logical                       :: passed

    call mpi_env_initialise

    passed = (world%err == 0)

    nx = 32
    ny = 64
    nz = 128
    lower = (/-pi, f12 * pi, zero/)
    extent = (/twopi, two * twopi, four * twopi/)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    allocate(fp1(1, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
    allocate(fp2(1, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
    allocate(fs(1, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
    allocate(fp_ref(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
    allocate(fs_ref(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

    call initialise_fft(extent)

    fp1(:, :, :) = zero
    fp2(:, :, :) = zero

    ! setup test field
    do i = box%lo(1), box%hi(1)
        x = lower(1) + dble(i) * dx(1)
        do j = box%lo(2), box%hi(2)
            y = lower(2) + dble(j) * dx(2)
            fp1(1, j, i) = dcos(four * x) + dsin(y)
        enddo
    enddo

    fp_ref = zero
    fp_ref(0, :, :) = fp1(1, :, :)

    fs = zero

    ! forward FFT
    call fft2d(fp1, fs)

    call fftxyp2s(fp_ref, fs_ref)

    error = maxval(dabs(fs(1,     box%lo(2):box%hi(2), box%lo(1):box%hi(1)) - &
                        fs_ref(0, box%lo(2):box%hi(2), box%lo(1):box%hi(1))))

    ! inverse FFT
    call ifft2d(fs, fp2)

    call finalise_fft

    error = max(error, maxval(dabs(fp1(1, box%lo(2):box%hi(2), box%lo(1):box%hi(1)) - &
                                   fp2(1, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))))

    call mpi_env_finalise

    passed = (passed .and. (world%err == 0))

    if (world%rank == world%root) then
        if (.not. passed) then
            call print_result_logical('Test FFT 2D transform', passed)
        else
            call print_result_dp('Test FFT 2D transform', error, atol=dble(1.0e-14))
        endif
    endif

end program test_fft2d
