! =============================================================================
!                               Test diffz
!
!  This unit test checks the calculation of the vertical differentiation of the
!  function
!       f(x, y, z) = z * cos(k*x) * sin(l*y)
! =============================================================================
program test_diffz_2
    use unit_test
    use constants, only : zero, one, two, pi, twopi, f12
    use parameters, only : lower, update_parameters, dx, nx, ny, nz, extent
    use inversion_utils
    use mpi_timer
    use mpi_environment
    use mpi_layout
    use mpi_collectives
    use fields, only : flayout
    implicit none

    double precision              :: error
    double precision, allocatable :: dfdz_ref(:, :, :), dfdz(:, :, :)
    double precision, allocatable :: fp(:, :, :)
    integer                       :: ix, iy, iz
    double precision              :: x, y, z, k, l, coskx, sinly

    call mpi_env_initialise

    nx = 16
    ny = 32
    nz = 32

    lower  = (/zero, -f12*pi, zero/)
    extent = (/pi, twopi, twopi/)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    allocate(fp(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
    allocate(dfdz(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
    allocate(dfdz_ref(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

    dfdz = zero

    call update_parameters

    call init_inversion

    k = two
    l = one

    do ix = box%lo(1), box%hi(1)
        x = lower(1) + ix * dx(1)
        coskx = dcos(k * x)
        do iy = box%lo(2), box%hi(2)
            y = lower(2) + iy * dx(2)
            sinly = dsin(l * y)
            do iz = 0, nz
                z = lower(3) + iz * dx(3)

                fp(iz, iy, ix) = z * coskx * sinly
                dfdz_ref(iz, iy, ix) = coskx * sinly
            enddo
        enddo
    enddo

    call flayout%diffz(fp, dfdz)

    error = maxval(dabs(dfdz_ref - dfdz))

    call mpi_blocking_reduce(error, MPI_MAX, world)

    if (world%rank == world%root) then
        call print_result_dp('Test diffz', error, atol=2.0e-14)
    endif

    deallocate(fp)
    deallocate(dfdz)
    deallocate(dfdz_ref)

    call mpi_env_finalise

end program test_diffz_2
