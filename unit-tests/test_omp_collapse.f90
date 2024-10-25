! =============================================================================
!                               Test OMP collapse
! =============================================================================
program test_omp_collapse
    use unit_test
    use constants, only : one, two, pi, f12
    use parameters, only : lower, update_parameters, nx, ny, nz, extent
    use sta2dfft, only : dst
    use sta3dfft, only : zfactors, ztrig
    use fields
    use omp_lib
    use mpi_environment
    use mpi_layout
    use mpi_collectives, only : mpi_blocking_reduce
    use model, only : layout, create_model
    implicit none

    double precision              :: error
    double precision, allocatable :: x(:), y(:), z(:)
    integer                       :: ix, iy, iz, nc, kx, ky
    double precision              :: alpha, fk2l2, k, l, m
    double precision              :: cosmz, sinmz, sinkxly, coskxly

    call mpi_env_initialise

    nx = 64
    ny = 64
    nz = 64

    lower  = -f12 * pi * (/one, one, one/)
    extent =  pi * (/one, one, one/)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    call field_default

    call create_model('uniform', "Hou & Li")

    allocate(x(0:nx-1), y(0:ny-1), z(0:nz))

    k = two
    l = two
    m = one

    alpha = dsqrt(k ** 2 + l ** 2 + m ** 2)
    fk2l2 = one / dble(k ** 2 + l ** 2)

    x = layout%get_x_axis()
    y = layout%get_y_axis()
    z = layout%get_z_axis()
    do ix = box%lo(1), box%hi(1)
        do iy = box%lo(2), box%hi(2)
            do iz = 0, nz
                cosmz = cos(m * z(iz))
                sinmz = sin(m * z(iz))
                sinkxly = sin(k * x(ix) + l * y(iy))
                coskxly = cos(k * x(ix) + l * y(iy))

                ! velocity
                vel(iz, iy, ix, 1) = fk2l2 * (k * m * sinmz - l * alpha * cosmz) * sinkxly
                vel(iz, iy, ix, 2) = fk2l2 * (l * m * sinmz + k * alpha * cosmz) * sinkxly
                vel(iz, iy, ix, 3) = cosmz * coskxly
            enddo
        enddo
    enddo


    do nc = 1, 3
        call layout%decompose_physical(vel(:, :, :, nc), svel(:, :, :, nc))
    enddo

    svor = svel

    do kx = box%lo(1), box%hi(1)
        do ky = box%lo(2), box%hi(2)
            call dst(1, nz, svor(1:nz, ky, kx, 1), ztrig, zfactors)
            call dst(1, nz, svor(1:nz, ky, kx, 2), ztrig, zfactors)
            call dst(1, nz, svor(1:nz, ky, kx, 3), ztrig, zfactors)
        enddo
    enddo

    !$omp parallel do collapse(2)
    do kx = box%lo(1), box%hi(1)
        do ky = box%lo(2), box%hi(2)
            call dst(1, nz, svel(1:nz, ky, kx, 1), ztrig, zfactors)
            call dst(1, nz, svel(1:nz, ky, kx, 2), ztrig, zfactors)
            call dst(1, nz, svel(1:nz, ky, kx, 3), ztrig, zfactors)
        enddo
    enddo
    !$omp end parallel do


    error = maxval(abs(svel - svor))

    call mpi_blocking_reduce(error, MPI_MAX, world)

    if (world%rank == world%root) then
        call print_result_dp('Test OMP collapse', error, atol=1.0d-15)
    endif

    deallocate(x, y, z)

    call mpi_env_finalise

end program test_omp_collapse
