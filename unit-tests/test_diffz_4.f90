! =============================================================================
!                               Test diffz
!
!  This unit test checks the calculation of the vertical differentiation of the
!  function
!       f(x, y, z) = 1 - 3*z^2
! =============================================================================
program test_diffz_3
    use unit_test
    use constants, only : zero, one, three, six
    use parameters, only : lower, update_parameters, dx, nx, ny, nz, extent
    use inversion_utils
    use timer
    implicit none

    double precision              :: error
    double precision, allocatable :: dfdz_ref(:, :, :), dfdz(:, :, :), ds(:, :, :)
    double precision, allocatable :: fp(:, :, :)
    integer                       :: iz
    double precision              :: z

    nx = 16
    ny = nx
    nz = nx

    lower  = (/zero, zero, zero/)
    extent = (/one, one, one/)

    allocate(fp(0:nz, 0:ny-1, 0:nx-1))
    allocate(ds(0:nz, 0:nx-1, 0:ny-1))
    allocate(dfdz(0:nz, 0:ny-1, 0:nx-1))
    allocate(dfdz_ref(0:nz, 0:ny-1, 0:nx-1))

    ds = zero
    dfdz = zero

    call update_parameters

    call init_inversion

    do iz = 0, nz
        z = lower(3) + iz * dx(3)
        fp(iz, :, :) = one - three * z ** 2
        dfdz_ref(iz, :, :) = - six * z
    enddo

    call diffz(fp, ds)

    call field_combine_physical(ds, dfdz)

    error = maxval(dabs(dfdz_ref - dfdz))

    print *, error

!     call print_result_dp('Test diffz', error, atol=)

    deallocate(fp)
    deallocate(ds)
    deallocate(dfdz)
    deallocate(dfdz_ref)

end program test_diffz_3
