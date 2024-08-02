! =============================================================================
!                               Test vor2vel
!
!  This unit test checks the calculation of the velocity field using the
!  Beltrami flow:
!     u(x, y, z) = (k^2 + l^2)^(-1) * [k*m*sin(mz) - l*alpha*cos(m*z) * sin(k*x + l*y)]
!     v(x, y, z) = (k^2 + l^2)^(-1) * [l*m*sin(mz) + k*alpha*cos(m*z) * sin(k*x + l*y)]
!     w(x, y, z) = cos(m*z) * cos(k*x + l*y)
!  The vorticity of this flow is
!    xi(x, y, z) = alpha * u(x, y, z)
!   eta(x, y, z) = alpha * v(x, y, z)
!  zeta(x, y, z) = alpha * w(x, y, z)
! =============================================================================
program test_vor2vel_1
    use unit_test
    use constants, only : one, two, pi, f12
    use parameters, only : lower, update_parameters, dx, nx, ny, nz, extent
    use fields
    use sta3dfft, only : fftxyp2s
    use inversion_utils, only : init_inversion
    use inversion_mod, only : vor2vel, vor2vel_timer
    use mpi_timer
    use mpi_environment
    use mpi_layout
    use mpi_collectives, only : mpi_blocking_reduce
    implicit none

    double precision              :: error
    double precision, allocatable :: vel_ref(:, :, :, :)
    integer                       :: ix, iy, iz
    double precision              :: x, y, z, alpha, fk2l2, k, l, m
    double precision              :: cosmz, sinmz, sinkxly, coskxly

    call mpi_env_initialise

    call register_timer('vorticity', vor2vel_timer)

    nx = 32
    ny = 32
    nz = 32

    lower  = -f12 * pi * (/one, one, one/)
    extent =  pi * (/one, one, one/)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    allocate(vel_ref(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1), 3))

    call update_parameters

    call field_default

    k = two
    l = two
    m = one

    call init_inversion


    alpha = dsqrt(k ** 2 + l ** 2 + m ** 2)
    fk2l2 = one / dble(k ** 2 + l ** 2)

    do ix = box%lo(1), box%hi(1)
        x = lower(1) + ix * dx(1)
        do iy = box%lo(2), box%hi(2)
            y = lower(2) + iy * dx(2)
            do iz = 0, nz
                z = lower(3) + iz * dx(3)

                cosmz = dcos(m * z)
                sinmz = dsin(m * z)
                sinkxly = dsin(k * x + l * y)
                coskxly = dcos(k * x + l * y)

                ! velocity
                vel_ref(iz, iy, ix, 1) = fk2l2 * (k * m * sinmz - l * alpha * cosmz) * sinkxly
                vel_ref(iz, iy, ix, 2) = fk2l2 * (l * m * sinmz + k * alpha * cosmz) * sinkxly
                vel_ref(iz, iy, ix, 3) = cosmz * coskxly

                ! vorticity
                vor(iz, iy, ix, 1) = alpha * vel_ref(iz, iy, ix, 1)
                vor(iz, iy, ix, 2) = alpha * vel_ref(iz, iy, ix, 2)
                vor(iz, iy, ix, 3) = alpha * vel_ref(iz, iy, ix, 3)

            enddo
        enddo
    enddo

    call fftxyp2s(vor(:, :, :, 1), svor(:, :, :, 1))
    call fftxyp2s(vor(:, :, :, 2), svor(:, :, :, 2))
    call fftxyp2s(vor(:, :, :, 3), svor(:, :, :, 3))

    call vor2vel

    error = maxval(dabs(vel_ref - vel))

    call mpi_blocking_reduce(error, MPI_MAX, world)

    if (world%rank == world%root) then
        call print_result_dp('Test vor2vel', error, atol=1.0e-14)
    endif

    deallocate(vel_ref)

    call mpi_env_finalise

end program test_vor2vel_1
