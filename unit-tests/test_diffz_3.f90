! =============================================================================
!                               Test diffz
!
!  This unit test checks the calculation of the vertical differentiation of the
!  function
!       f(x, y, z) = 6*z + 3*z^2 * cos(k*x) * sin(l*y)
! =============================================================================
program test_diffz_3
    use unit_test
    use constants, only : zero, one, two, pi, twopi, f12, six, three
    use parameters, only : lower, update_parameters, dx, nx, ny, nz, extent
    use inversion_utils
    use timer
    implicit none

    double precision              :: error
    double precision, allocatable :: dfdz_ref(:, :, :), dfdz(:, :, :)
    double precision, allocatable :: fp(:, :, :)
    integer                       :: ix, iy, iz
    double precision              :: x, y, z, k, l, coskx, sinly

    nx = 16
    ny = 2*nx
    nz = 2*nx

    lower  = (/zero, -f12*pi, zero/)
    extent = (/pi, twopi, twopi/)

    allocate(fp(0:nz, 0:ny-1, 0:nx-1))
    allocate(dfdz(0:nz, 0:ny-1, 0:nx-1))
    allocate(dfdz_ref(0:nz, 0:ny-1, 0:nx-1))

    dfdz = zero

    call update_parameters

    call init_inversion

    k = two
    l = one

    do ix = 0, nx-1
        x = lower(1) + ix * dx(1)
        coskx = dcos(k * x)
        do iy = 0, ny-1
            y = lower(2) + iy * dx(2)
            sinly = dsin(l * y)
            do iz = 0, nz
                z = lower(3) + iz * dx(3)

                fp(iz, iy, ix) = six * z + three * coskx * sinly * z
                dfdz_ref(iz, iy, ix) = six + three * coskx * sinly
            enddo
        enddo
    enddo

    call diffz(fp, dfdz)

    error = maxval(dabs(dfdz_ref - dfdz))

    call print_result_dp('Test diffz', error, atol=2.0e-13)

    deallocate(fp)
    deallocate(dfdz)
    deallocate(dfdz_ref)

end program test_diffz_3