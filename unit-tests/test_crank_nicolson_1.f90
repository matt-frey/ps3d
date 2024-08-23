! =============================================================================
!                          Test Crank-Nicolson method
!
! This test solves the diffusion equation
!       u_t = (D u_z)_z
!           = D_z u_z + D u_zz
! where D = D(z) is the diffusivity and subscripts t or z denote derivatives
! with respect to that variable. This method uses the time-average approach
! according to
!
!   Dritschel D. G., Jalali M. R.
!   The validity of two-dimensional models of a rotating shallow fluid layer.
!   Journal of Fluid Mechanics. 2020;900:A33. doi:10.1017/jfm.2020.487
!
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
    use zops, only : d1z, d2z
    implicit none

    double precision, allocatable :: u(:), d_x(:), eye(:, :), Lm(:, :), Am(:, :), v(:)
    double precision, allocatable :: sol(:, :)
    double precision              :: a, b, dt
    double precision              :: gavg, alpha, rk, fac, l1norm_initial, l1norm_final
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
    allocate(d_x(0:nz))
    allocate(eye(0:nz, 0:nz))
    allocate(Am(0:nz, 0:nz))
    allocate(Lm(0:nz, 0:nz))
    allocate(ipiv(0:nz))

    call update_parameters

    dt = 0.0025d0
    nt = 1.0d0 / dt

    allocate(sol(0:int(nt/nsave), 0:nz))

    call init_inversion

    a = 0.75d0
    b = 0.75d0
    u = one - two * ( (zcheb - b) / a) * dexp(-( (zcheb - b) / a) ** 2) + 0.1d0 * dsin(3.0d0 * pi * zcheb)

    ! Diffusivity (hyperbolic function)
    d_x = (tanh(5.0d0 * (zcheb + 1.0d0)) * tanh(5.0d0 * (1.0d0 - zcheb))) ** 2

!     do iz = 0, nz
!         print *, zcheb(iz), d_x(iz)
!     enddo

    ! Create the system matrix for the Crank-Nicolson method
    eye = zero
    Am = zero
    do iz = 0, nz
        eye(iz, iz) = one
        Am(iz, iz) = d_x(iz)
    enddo

    Am = matmul(matmul(d1z, Am), d1z)

    Lm = eye - f12 * dt * Am

    l1norm_initial = dot_product(zccw, u)

    sol(0, :) = u
    i = 1
    do j = 1, nt
        ! right-hand side u^{n}
        v = u

        ! solve the linear system to get \bar{u}
        Am = Lm
        ipiv = 0
        call dgesv(nz+1, 1, Am, nz+1, ipiv, v, nz+1, info)

        ! update u^{n+1} = 2\bar{u} - u^{n}
        u = 2.0d0 * v - u

        ! compute L1 norm
!         l1norm = dot_product(zccw, u)
!         print *, j * dt, l1norm

        if (mod(j, nsave) == 0) then
            sol(i, :) = u
            i = i + 1
        endif
    enddo

    l1norm_final = dot_product(zccw, u)

    if (world%rank == world%root) then
        ! L1 norm should remain constant.
        call print_result_dp('Test Crank-Nicolson 1', dabs(l1norm_final - l1norm_initial), atol=1.0e-13)
    endif

    if (verbose .and. (world%rank == world%root)) then
        do iz = 0, nz
            print *, zcheb(iz), sol(:, iz)
        enddo
    endif

    call finalise_inversion

    deallocate(u, v, sol)
    deallocate(d_x, ipiv)
    deallocate(eye, Am, Lm)

    call mpi_env_finalise

end program test_crank_nicolson
