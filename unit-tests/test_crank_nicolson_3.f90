! =============================================================================
!                          Test Crank-Nicolson method
!
! This test solves
!       u_t = (alpha_v(z) u_z)_z - alpha_z K^2 u
! where alpha_v = alpha_v(z) is the vertical diffusivity and
! \alpha_h is the horizontal diffusivity, u = u(z, t) and
! subscripts t or z denote
! derivatives with respect to that variable. We assume K = 1.
!
! We use u(z, t=0) = exp(z) and
!
!           alpha_v(z) = [tanh(5 * (1 + z)) * tanh(5 * (1 - z))]^2
! =============================================================================
program test_crank_nicolson
    use unit_test
    use constants, only : zero, f12, one, two, pi, ten
    use parameters, only : lower, nx, ny, nz, extent, hl &
                         , update_parameters, center, upper
    use fields
    use inversion_utils
    use sta3dfft, only : k2l2
    use mpi_timer
    use mpi_environment
    use mpi_layout
    use mpi_collectives, only : mpi_blocking_reduce
    use zops, only : d1z
    implicit none

    double precision, allocatable :: u(:), alpha_v(:), eye(:, :), Lm(:, :), Rm(:, :), Am(:, :), v(:)
    double precision, allocatable :: sol(:, :)
    double precision              :: dt, k2
    double precision              :: gavg, alpha_h, rk, fac, l1norm_initial, l1norm_final
    integer                       :: iz, i, j, nt, info, nsave
    integer, allocatable          :: ipiv(:)

    call mpi_env_initialise

    call parse_command_line

    nx = 16
    ny = 16
    nz = 64
    nsave = 50

    lower  = (/-1.0d0, -1.0d0, -1.0d0/)
    extent = (/2.0d0, 2.0d0, 2.0d0/)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    allocate(u(0:nz))
    allocate(v(0:nz))
    allocate(alpha_v(0:nz))
    allocate(eye(0:nz, 0:nz))
    allocate(Am(0:nz, 0:nz))
    allocate(Lm(0:nz, 0:nz))
    allocate(Rm(0:nz, 0:nz))
    allocate(ipiv(0:nz))

    call update_parameters

    dt = 0.0025d0
    nt = 1.0d0 / dt

    allocate(sol(0:int(nt/nsave), 0:nz))

    call init_inversion

    u = exp(zcheb)

    ! Diffusivity (hyperbolic function)
    alpha_v = (tanh(5.0d0 * (zcheb + 1.0d0)) * tanh(5.0d0 * (1.0d0 - zcheb))) ** 2

    alpha_h = maxval(alpha_v)

    ! Create the system matrix for the Crank-Nicolson method
    eye = zero
    Am = zero
    do iz = 0, nz
        eye(iz, iz) = one
        Am(iz, iz) = alpha_v(iz)
    enddo

    k2 = 1.0d0

    Lm = eye + f12 * dt * k2 * alpha_h * eye
    Rm = eye - f12 * dt * k2 * alpha_h * eye

    Am = matmul(matmul(d1z, Am), d1z)

    Lm = Lm - f12 * dt * Am
    Rm = Rm + f12 * dt * Am

    l1norm_initial = dot_product(zccw, u)

    sol(0, :) = u
    i = 1
    do j = 1, nt
        ! right-hand side
        v = matmul(Rm, u)

        ! solve the linear system
        Am = Lm
        ipiv = 0
        call dgesv(nz+1, 1, Am, nz+1, ipiv, v, nz+1, info)

        ! compute L1 norm
        u = v

        if (mod(j, nsave) == 0) then
            sol(i, :) = u
            i = i + 1
        endif
    enddo

    l1norm_final = dot_product(zccw, u)

    if (world%rank == world%root) then
        call print_result_logical('Test Crank-Nicolson 3', (l1norm_final < l1norm_initial))
    endif

    if (verbose .and. (world%rank == world%root)) then
        do iz = 0, nz
            print *, zcheb(iz), sol(:, iz)
        enddo
    endif

    call finalise_inversion

    deallocate(u, v, sol)
    deallocate(alpha_v, ipiv)
    deallocate(eye, Am, Lm, Rm)

    call mpi_env_finalise

end program test_crank_nicolson
