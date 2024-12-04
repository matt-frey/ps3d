! =============================================================================
!                          Test Crank-Nicolson method
!
! This test solves
!       u_t = (\alpha_v(z) u_z)_z - \alpha_h K^2 u
! where D = D(z) is the vertical diffusivity and
! \alpha_h is the horizontal diffusivity, u = u(z, t) and subscripts t or z denote
! derivatives with respect to that variable. We assume K = 1. Here, we use
! the time-average approach according to
!
!   Dritschel D. G., Jalali M. R.
!   The validity of two-dimensional models of a rotating shallow fluid layer.
!   Journal of Fluid Mechanics. 2020;900:A33. doi:10.1017/jfm.2020.487
!
! We use u(z, t=0) = exp(z) and
!
!           \alpha_v(z) = [tanh(5 * (1 + z)) * tanh(5 * (1 - z))]^2
! =============================================================================
program test_crank_nicolson
    use unit_test
    use constants, only : zero, f12, one, two, pi, ten
    use parameters, only : lower, nx, ny, nz, extent &
                         , update_parameters
    use fields
    use mpi_timer
    use mpi_environment
    use mpi_layout
    use mpi_collectives, only : mpi_blocking_reduce
    use cheby_layout, only : cheby_layout_t
    implicit none

    double precision, allocatable :: u(:), alpha_v(:), eye(:, :), Lm(:, :), Am(:, :), v(:)
    double precision, allocatable :: sol(:, :)
    double precision              :: dt, k2, alpha_h
    double precision              :: l1norm_initial, l1norm_final
    integer                       :: iz, i, j, nt, info, nsave
    integer, allocatable          :: ipiv(:)
    type(cheby_layout_t)          :: layout

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
    allocate(ipiv(0:nz))

    call update_parameters

    dt = 0.0025d0
    nt = int(1.0d0 / dt)

    allocate(sol(0:int(nt/nsave), 0:nz))

    call layout%initialise

    u = exp(layout%zcheb)

    ! Diffusivity (hyperbolic function)
    alpha_v = (tanh(5.0d0 * (layout%zcheb + 1.0d0)) * tanh(5.0d0 * (1.0d0 - layout%zcheb))) ** 2

    alpha_h = maxval(alpha_v)

    ! Create the system matrix for the Crank-Nicolson method
    eye = zero
    Am = zero
    do iz = 0, nz
        eye(iz, iz) = one
        Am(iz, iz) = alpha_v(iz)
    enddo

    k2 = 1.0d0

    Lm = eye + f12 * k2 * dt * alpha_h * eye

    Am = matmul(matmul(layout%d1z, Am), layout%d1z)

    Lm = Lm - f12 * dt * Am

    l1norm_initial = dot_product(layout%zccw, u)

    sol(0, :) = u
    i = 1
    do j = 1, nt
        ! right-hand side u^{n}
        v = u

        ! solve the linear system to get \bar{u}
        Am = Lm
        ipiv = 0
        call dgesv(nz+1, 1, Am, nz+1, ipiv, v, nz+1, info)

        ! get u^{n+1} = 2\bar{u} - u^{n}
        u = 2.0d0 * v - u

        ! compute L1 norm
        if (mod(j, nsave) == 0) then
            sol(i, :) = u
            i = i + 1
        endif
    enddo

    l1norm_final = dot_product(layout%zccw, u)

    if (world%rank == world%root) then
        call print_result_logical('Test Crank-Nicolson 2', (l1norm_final < l1norm_initial))
    endif

    if (verbose .and. (world%rank == world%root)) then
        do iz = 0, nz
            print *, layout%zcheb(iz), sol(:, iz)
        enddo
    endif

    deallocate(u, v, sol)
    deallocate(alpha_v, ipiv)
    deallocate(eye, Am, Lm)

    call mpi_env_finalise

end program test_crank_nicolson
