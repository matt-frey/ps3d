! =============================================================================
!                 Test integration using Clenshaw-Curtis weights
!                 for arbitrary integration bounds.
!
!   This unit test integrates the function f(z) = 1
!   on the interval [0, 1].
!
!   Note: The Chebyshev points are defined on the interval [-1, 1].
!         We therefore need to apply the proper rescaling.
! =============================================================================
program test_clenshaw
    use unit_test
    use constants, only : zero, one, f12
    use parameters, only : lower, nx, ny, nz, extent &
                         , update_parameters
    use mpi_environment
    use mpi_layout
    use inversion_utils
    implicit none

    double precision :: integral_approx, exact_integral
    double precision :: rsum, f, error
    integer          :: i

    call mpi_env_initialise

    call parse_command_line

    nx = 16
    ny = 16
    nz = 16

    lower  = (/zero, zero, zero/)
    extent = (/one, one, one/)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    !-----------------------------------------------------------------
    ! Compute Chebyshev nodes and Clenshaw-Curtis weights
    call init_inversion

    ! Define the exact integral of f(z) = 1 over [0, 1]
    exact_integral = one

    ! Compute the integral using Clenshaw-Curtis quadrature
    rsum = zero
    do i = 0, nz
        ! Define the function f(z) = 1
        f = 1.0d0

        rsum = rsum + zccw(i) * f
    enddo

    ! Map from Chebyhshev points t in [-1, 1] to [0, 1]
    ! z = (b-a) / 2 * t + (a+b)/2 for [a, b]
    ! z = 1 / 2 * t + 1/2 --> dz = 1/ 2 * dt
    integral_approx = f12 * rsum

    error = abs(exact_integral-integral_approx)

    ! Output the results
    if (verbose .and. (world%rank == world%root)) then
        print *, 'Clenshaw-Curtis nodes and weights for N =', nz
        print *, '  Node (z)       Weight (w)'
        do i = 0, nz
            print *, zcheb(i), zccw(i)
        enddo
        print *, 'Exact integral of f(z) = 1 over [0, 1]:', exact_integral
        print *, 'Integral using Clenshaw-Curtis quadrature:', integral_approx
        print *, ''
        print *, 'error: ', error
    endif

    if (world%rank == world%root) then
        call print_result_dp('Test Clenshaw-Curtis weights', error, atol=1.0e-13)
    endif

    call finalise_inversion

    call mpi_env_finalise

end program test_clenshaw

