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
    use constants, only : one, two, three, six, pi, twopi
    use parameters, only : lower, update_parameters, dx, nx, ny, nz, extent
    use fields
    use inversion_utils
    use inversion_mod, only : vor2vel, vor2vel_timer
    use utils, only : setup_domain_and_parameters, write_step
    use field_netcdf, only : field_io_timer, create_netcdf_field_file
    use timer
    implicit none

    double precision              :: error
    double precision, allocatable :: vel_ref(:, :, :, :)
    integer                       :: ix, iy, iz
    double precision              :: x, y, z, k2l2, k, l, f, dfdz, d2fdz2
    double precision              :: coskx, sinkx, cosly, sinly

    call register_timer('vorticity', vor2vel_timer)
    call register_timer('field I/O', field_io_timer)

    nx = 32
    ny = 32
    nz = 32

    lower  = (/-f12, -f12, zero/)
    extent = (/one, one, one/)

    allocate(vel_ref(0:nz, 0:ny-1, 0:nx-1, 3))

    call update_parameters

    call field_default

    l = twopi
    k = two * l

    call init_inversion


    k2l2 = k ** 2 - l ** 2

    do ix = 0, nx-1
        x = lower(1) + ix * dx(1)
        do iy = 0, ny-1
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
                vel_ref(iz, iy, ix, 3) =  k2l2 * f * sinkx * sinly

                ! vorticity
                vor(iz, iy, ix, 1) = l * (k2l2 * f + d2fdz2) * sinkx * cosly
                vor(iz, iy, ix, 2) = k * (d2fdz2 - k2l2 * f) * coskx * sinly
                vor(iz, iy, ix, 3) = - two * k * l * dfdz * coskx * cosly

            enddo
        enddo
    enddo

    call field_decompose_physical(vor(:, :, :, 1), svor(:, :, :, 1))
    call field_decompose_physical(vor(:, :, :, 2), svor(:, :, :, 2))
    call field_decompose_physical(vor(:, :, :, 3), svor(:, :, :, 3))

    call vor2vel

    error = maxval(dabs(vel_ref - vel))

    call create_netcdf_field_file('test', .true.)

    call write_step(zero)

    call print_result_dp('Test vor2vel', error, atol=6.0e-3)

    deallocate(vel_ref)

end program test_vor2vel_2
