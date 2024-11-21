! =============================================================================
!                               Test vor2vel
!
!  This unit test checks the calculation of the velocity field using the
!  following flow:
!     u(x, y, z) =  k * df/dz * cos(k*x) * sin(l*y)
!     v(x, y, z) = -l * df/dz * sin(k*x) * cos(l*y)
!     w(x, y, z) = (k^2 - l^2) * f * sin(k*x) * sin(l*y)
!  and vorticity
!    xi(x, y, z) = l * [(k^2 - l^2) * f + d^2f/dz^2] * sin(k*x) * cos(l*y)
!   eta(x, y, z) = k * [d^2f/dz^2 - (k^2 - l^2) * f] * cos(k*x) * sin(l*y)
!  zeta(x, y, z) = - 2 * k * l * df/dz * cos(k*x) * cos(l*y)
!  where
!           f(z) = 2*z - z^2 - z^3
!          df/dz = 2 - 2*z - 3*z^2
!      d^2f/dz^2 = -2 - 6*z
! =============================================================================
program test_vor2vel_2
    use unit_test
    use constants, only : one, two, three, six, pi, twopi, f12
    use parameters, only : lower, update_parameters, dx, nx, ny, nz, extent
    use fields
    use inversion_utils
    use inversion_mod, only : vor2vel, vor2vel_timer
    use mpi_timer
    use mpi_environment
    use mpi_layout
    use mpi_collectives, only : mpi_blocking_reduce
    implicit none

    double precision              :: error
    double precision, allocatable :: vel_ref(:, :, :, :)
    integer                       :: ix, iy, iz
    double precision              :: x, y, z, klsq, k, l, f, dfdz, d2fdz2
    double precision              :: coskx, sinkx, cosly, sinly

    call mpi_env_initialise

    call register_timer('vorticity', vor2vel_timer)

    nx = 32
    ny = 32
    nz = 32

    lower  = (/-f12, -f12, zero/)
    extent = (/one, one, one/)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    allocate(vel_ref(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1), 3))

    call update_parameters

    call field_default

    l = twopi
    k = two * l

    call init_inversion


    klsq = k ** 2 - l ** 2

    do ix = box%lo(1), box%hi(1)
        x = lower(1) + ix * dx(1)
        do iy = box%lo(2), box%hi(2)
            y = lower(2) + iy * dx(2)
            do iz = 0, nz
                z = lower(3) + iz * dx(3)

                f = two * z - z ** 2 - z ** 3
                dfdz = two - two * z - three * z ** 2
                d2fdz2 = -two - six * z
                sinkx = dsin(k * x)
                coskx = dcos(k * x)
                sinly = dsin(l * y)
                cosly = dcos(l * y)

                ! velocity
                vel_ref(iz, iy, ix, 1) =  k * dfdz * coskx * sinly
                vel_ref(iz, iy, ix, 2) = -l * dfdz * sinkx * cosly
                vel_ref(iz, iy, ix, 3) =  klsq * f * sinkx * sinly

                ! vorticity
                vor(iz, iy, ix, 1) = l * (klsq * f + d2fdz2) * sinkx * cosly
                vor(iz, iy, ix, 2) = k * (d2fdz2 - klsq * f) * coskx * sinly
                vor(iz, iy, ix, 3) = - two * k * l * dfdz * coskx * cosly

            enddo
        enddo
    enddo

    call flayout%decompose_physical(vor(:, :, :, 1), svor(:, :, :, 1))
    call flayout%decompose_physical(vor(:, :, :, 2), svor(:, :, :, 2))
    call flayout%decompose_physical(vor(:, :, :, 3), svor(:, :, :, 3))

    call vor2vel

    error = maxval(dabs(vel_ref - vel))

    call mpi_blocking_reduce(error, MPI_MAX, world)

    if (world%rank == world%root) then
        call print_result_dp('Test vor2vel', error, atol=1.2e-2)
    endif

    deallocate(vel_ref)

    call mpi_env_finalise

end program test_vor2vel_2
