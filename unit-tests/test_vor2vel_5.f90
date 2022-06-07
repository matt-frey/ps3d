! =============================================================================
!                               Test vor2vel
!
!  This unit test checks the calculation of the velocity field using the
!  following flow:
!     u(x, y, z) = z^3 - 1/4
!     v(x, y, z) = 0
!     w(x, y, z) = 0
!  and vorticity
!    xi(x, y, z) = 0
!   eta(x, y, z) = 3 * z^2
!  zeta(x, y, z) = 0
! with x, y in [-1/2, 1/2] and z in [0, 1].
! =============================================================================
program test_vor2vel_5
    use unit_test
    use constants, only : one, three, f12, f14, two
    use parameters, only : lower, update_parameters, dx, nx, ny, nz, extent, upper
    use fields
    use inversion_utils
    use inversion_mod, only : vor2vel, vor2vel_timer
    use utils, only : setup_domain_and_parameters, write_step
    use field_netcdf, only : field_io_timer, create_netcdf_field_file
    use timer
    implicit none

    double precision              :: error
    double precision, allocatable :: vel_ref(:, :, :, :)
    integer                       :: iz
    double precision              :: z

    call register_timer('vorticity', vor2vel_timer)
    call register_timer('field I/O', field_io_timer)

    nx = 32
    ny = nx
    nz = nx

    lower  = (/-f12, -f12, zero/)
    extent = (/one, one, one/)

    allocate(vel_ref(0:nz, 0:ny-1, 0:nx-1, 3))

    call update_parameters

    call field_default

    call init_inversion

    do iz = 0, nz
        z = lower(3) + iz * dx(3)

        ! velocity
        vel_ref(iz, :, :, 1) = z ** 3 - f14
        vel_ref(iz, :, :, 2) = zero
        vel_ref(iz, :, :, 3) = zero

        ! vorticity
        vor(iz, :, :, 1) = zero
        vor(iz, :, :, 2) = three * z ** 2
        vor(iz, :, :, 3) = zero
    enddo

    call field_decompose_physical(vor(:, :, :, 1), svor(:, :, :, 1))
    call field_decompose_physical(vor(:, :, :, 2), svor(:, :, :, 2))
    call field_decompose_physical(vor(:, :, :, 3), svor(:, :, :, 3))

    call vor2vel

    error = maxval(dabs(vel_ref - vel))

    call create_netcdf_field_file('test', .true.)

    call write_step(zero)

    call print_result_dp('Test vor2vel', error, atol=4.0e-6)

    deallocate(vel_ref)

end program test_vor2vel_5
