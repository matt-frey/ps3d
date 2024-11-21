! =============================================================================
!                               Test zeta combination
!
! =============================================================================
program test_zeta
    use unit_test
    use constants, only : one, two, pi, f12, twopi
    use parameters, only : lower, update_parameters, dx, nx, ny, nz, extent
    use fields
    use inversion_utils
    use mpi_timer
    use mpi_environment
    use mpi_layout
    implicit none

    double precision              :: error
    double precision, allocatable :: zeta_ref(:, :, :)
    integer                       :: ix, iy, iz
    double precision              :: x, y, z

    call mpi_env_initialise

    nx = 32
    ny = 32
    nz = 32

    lower  = -f12 * pi * (/one, one, one/)
    extent =  twopi * (/one, one, one/)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    allocate(zeta_ref(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

    call update_parameters

    call field_default

    call init_inversion

    do ix = box%lo(1), box%hi(1)
        x = lower(1) + ix * dx(1)
        do iy = box%lo(2), box%hi(2)
            y = lower(2) + iy * dx(2)
            do iz = 0, nz
                z = lower(3) + iz * dx(3)

                zeta_ref(iz, iy, ix) = dsin(x) * dsin(z)

                ! vorticity
                vor(iz, iy, ix, 1) = dcos(x) * dcos(z)
                vor(iz, iy, ix, 2) = 0.0d0
            enddo
            zeta(0, iy, ix) = dsin(x) * dsin(z)
        enddo
    enddo

    call flayout%field_decompose_physical(vor(:, :, :, 1), svor(:, :, :, 1))
    call flayout%field_decompose_physical(vor(:, :, :, 2), svor(:, :, :, 2))
    vor(0, :, :, 3) = zeta(0, :, :)
    call surf_fftxyp2s(vor(0, :, :, 3), szeta(0, :, :))

    call combine_zeta

    error = maxval(dabs(zeta_ref - vor(:, :, :, 3)))

    call mpi_blocking_reduce(error, MPI_MAX, world)

    if (world%rank == world%root) then
        call print_result_dp('Test zeta', error, atol=1.0e-14)
    endif

    deallocate(zeta_ref)

    call mpi_env_finalise

end program test_zeta
