! =============================================================================
!                               Test diffusion
!
!  This unit test checks the diffusion of the vorticity field using the
!  Beltrami flow:
!     u(x, y, z) = (k^2 + l^2)^(-1) * [k*m*sin(mz) - l*alpha*cos(m*z) * sin(k*x + l*y)]
!     v(x, y, z) = (k^2 + l^2)^(-1) * [l*m*sin(mz) + k*alpha*cos(m*z) * sin(k*x + l*y)]
!     w(x, y, z) = cos(m*z) * cos(k*x + l*y)
!  The vorticity of this flow is
!    xi(x, y, z) = alpha * u(x, y, z)
!   eta(x, y, z) = alpha * v(x, y, z)
!  zeta(x, y, z) = alpha * w(x, y, z)
! =============================================================================
program test_diffusion
    use unit_test
    use constants, only : one, two, pi, f12
    use parameters, only : lower, update_parameters, nx, ny, nz, extent
    use fields
    use mpi_timer
    use mpi_environment
    use diffusion
    use mpi_layout
    use mpi_collectives, only : mpi_blocking_reduce
    use options, only : vor_visc
    use model, only : layout, create_model
    implicit none

    call mpi_env_initialise

    nx = 32
    ny = 32
    nz = 32

    lower  = -f12 * pi * (/one, one, one/)
    extent =  pi * (/one, one, one/)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    vor_visc%length_scale = 'Kolmogorov'
    vor_visc%prediss = 30.0d0
    vor_visc%nnu = 3

    call update_parameters

    call field_default

    call run_test('uniform')

    call run_test('chebyshev')


    call mpi_env_finalise

contains

    subroutine run_test(grid_type)
        character(*), intent(in)      :: grid_type
        double precision              :: error
        double precision, allocatable :: vd1(:, :, :, :), vd2(:, :, :, :)
        double precision, allocatable :: x(:), y(:), z(:)
        integer                       :: ix, iy, iz, nc
        double precision              :: alpha, fk2l2, k, l, m
        double precision              :: cosmz, sinmz, sinkxly, coskxly, te, en

        call create_model(grid_type, "Hou & Li")

        allocate(vd1(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1), 3))
        allocate(vd2(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1), 3))
        allocate(x(0:nx-1), y(0:ny-1), z(0:nz))

        k = two
        l = two
        m = one

        if (.not. is_diffusion_initialised) then
            te = 0.28125d0
            en = 2.53113708762278d0
            call init_diffusion(te, en)
        endif

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

                    ! vorticity
                    vor(iz, iy, ix, 1) = alpha * fk2l2 * (k * m * sinmz - l * alpha * cosmz) * sinkxly
                    vor(iz, iy, ix, 2) = alpha * fk2l2 * (l * m * sinmz + k * alpha * cosmz) * sinkxly
                    vor(iz, iy, ix, 3) = alpha * cosmz * coskxly
                enddo
            enddo
        enddo

        call layout%decompose_physical(vor(:, :, :, 1), svor(:, :, :, 1))
        call layout%decompose_physical(vor(:, :, :, 2), svor(:, :, :, 2))
        call layout%decompose_physical(vor(:, :, :, 3), svor(:, :, :, 3))

        vd1 = svor
        vd2 = svor

        do nc = 1, 3
            call layout%combine_semi_spectral(vd1(:, :, :, nc))
            vd1(:, :, :, nc) = vhdis * vd1(:, :, :, nc)
            vd2(:, :, :, nc) = vhdis * vd2(:, :, :, nc)
            call layout%decompose_semi_spectral(vd1(:, :, :, nc))
        enddo

        error = maxval(abs(vd1 - vd2))

        call mpi_blocking_reduce(error, MPI_MAX, world)

        if (world%rank == world%root) then
            call print_result_dp('Test diffusion ' // grid_type, error, atol=1.0e-14)
        endif

        deallocate(vd1, vd2, x, y, z)

    end subroutine run_test

end program test_diffusion
