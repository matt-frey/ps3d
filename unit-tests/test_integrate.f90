! =============================================================================
!                          Test 3D field integration
!
!   This unit test integrates the function sin(x) * cos(x+y) * z
!   on the interval x = [0, pi/4], y = [0, pi/4] and z = [1/2, 3/2]
!
!   In x and y we use the trapezoidal rule where
!       integrate f(x) from a to b
!       --> dx / 2 * (f(x(0)) + 2f(x(1)) + 2f(x(2)) + ... + 2f(x(N-1)) + f(x(N))
!   Because f is periodic in x and y, we have x(0) = x(N), and thus the
!   rapezoidal rule simplifies to
!       dx * (f(x(0)) + f(x(1)) + f(x(2)) + ... + f(x(N-1)))
!
!   Note: The Chebyshev points are defined on the interval [-1, 1].
!         We therefore need to apply the proper rescaling in z.
! =============================================================================
program test_integrate
    use unit_test
    use constants, only : zero, f12, one, two, pi
    use parameters, only : lower, nx, ny, nz, extent, dx &
                         , update_parameters
    use mpi_environment
    use mpi_layout
    use zops
    implicit none

    double precision :: integral_approx, exact_integral
    double precision :: rsum, error, x, y, z
    integer          :: n, ix, iy, iz
    double precision, allocatable :: f(:, :, :)

    call mpi_env_initialise

    call parse_command_line

    nx = 32
    ny = 32
    nz = 32

    lower  = (/-one, -one, f12/)
    extent = (/two, two, one/)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    allocate(f(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

    call update_parameters

    !-----------------------------------------------------------------
    ! Compute Chebyshev nodes and Clenshaw-Curtis weights
    call init_zops

    ! Define the exact integral of f(x, y, z) = sin(x) * cos(x+y) * z
    exact_integral = 8.0d0 / 3.0d0

    ! Define the function f(x, y, z)
    do ix = box%lo(1), box%hi(1)
        x = lower(1) + dble(ix) * dx(1)
        do iy = box%lo(2), box%hi(2)
            y = lower(2) + dble(iy) * dx(2)
            do iz = 0, nz
                ! Map from Chebyhshev points t in [-1, 1] to [1/2, 3/2]
                ! z = (b-a) / 2 * t + (a+b)/2 for [a, b]
                z = f12 * zcheb(iz) + one

                f(iz, iy, ix) = x ** 2 + y ** 2 * z
            enddo
        enddo
    enddo

    ! Compute the integral using Clenshaw-Curtis quadrature
    rsum = zero

    do iz = 0, nz
        rsum = rsum + zccw(iz) * sum(f(iz, :, :))
    enddo

    ! dz = (b-a)/ 2 * dt = 1/2 * dz
    integral_approx = f12 * rsum * dx(1) * dx(2)

    error = abs(exact_integral-integral_approx)

    ! Output the results
    if (verbose .and. (world%rank == world%root)) then
        print *, 'Clenshaw-Curtis nodes and weights for N =', nz
        print *, '  Node (x)       Weight (w)'
        do iz = 0, nz
            print *, zcheb(iz), zccw(iz)
        enddo
        print *, 'Exact integral of f(x, y, z):', exact_integral
        print *, 'Integral using Clenshaw-Curtis quadrature:', integral_approx
        print *, ''
        print *, 'error: ', error
    endif

    if (world%rank == world%root) then
        call print_result_dp('Test field integration', error, atol=6.0e-3)
    endif

    call finalise_zops

    deallocate(f)

    call mpi_env_finalise

end program test_integrate

