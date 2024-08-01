! =============================================================================
!                               Test zops module
! =============================================================================
program test_zops
    use unit_test
    use constants, only : zero, f12, one, pi, twopi, f14
    use parameters, only : lower, dx, nx, ny, nz, extent, hl &
                         , update_parameters, center, upper
    use fields
    use inversion_utils, only : k2l2
    use mpi_timer
    use mpi_environment
    use mpi_layout
    use mpi_collectives, only : mpi_blocking_reduce
    use zops
    implicit none

    double precision, allocatable :: f(:, :, :), f1(:)
    double precision, allocatable :: g(:, :, :), g1(:)
    double precision, allocatable :: z(:)
    double precision              :: gavg, alpha, rk, fac, maxerr
    integer                       :: kx, ky, iz

    call mpi_env_initialise

    call parse_command_line

    nx = 16
    ny = 16
    nz = 8

    lower  = (/-pi, -pi, -f14 * pi/)
    extent = (/twopi, twopi, f14 * pi/)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    allocate(f(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
    allocate(g(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
    allocate(f1(0:nz))
    allocate(g1(0:nz))
    allocate(z(0:nz))

    call update_parameters

    !-----------------------------------------------------------------
    ! Initialise zops module:
    call init_zops

    ! Set up z grid:
    z = center(3) - hl(3) * zcheb

    if (verbose .and. (world%rank == world%root)) then
        write(*,'(a,i2,a)') '  @@@   Using nz = ',nz,' grid intervals in z   @@@'
        write(*,*)
        write(*,*) ' ==>  Testing f = exp(z)  <=='
    endif

    !-----------------------------------------------------------------
    ! Test first-order differentiation: f = exp(z)
    do kx = box%lo(1), box%hi(1)
        do ky = box%lo(2), box%hi(2)
            f(:, ky, kx) = exp(z)
        enddo
    enddo

    call zderiv(f, g)

    if (verbose) then
        maxerr = maxval(abs(g - f))
        call mpi_blocking_reduce(maxerr, MPI_MAX, world)
        if (world%rank == world%root) then
            write(*,*)
            write(*,*) ' For df/dz, the maximum abs error = ', maxerr
        endif
    endif

    !-----------------------------------------------------------------
    ! Test second-order differentiation:
    do kx = box%lo(1), box%hi(1)
        do ky = box%lo(2), box%hi(2)
            f(:, ky, kx) = exp(z)
        enddo
    enddo

    call zzderiv(f, g)

    if (verbose) then
        maxerr = maxval(abs(g - f))
        call mpi_blocking_reduce(maxerr, MPI_MAX, world)
        if (world%rank == world%root) then
            write(*,*)
            write(*,*) ' For d^2f/dz^2, the maximum abs error = ', maxerr
        endif
    endif

    !-----------------------------------------------------------------
    ! Test integration:

    f1 = exp(z)

    call zinteg(f1, g1, .true.)

    gavg = (one + lower(3) - exp(lower(3))) / extent(3)
    f1 = f1 - one - gavg

    if (verbose .and. (world%size == 1)) then
        write(*,*)
        write(*,*) '    z        exact integal    approx integral       error'
        do iz = 0, nz
            write(*,'(f10.7,2(2x,f15.10),2x,1p,e14.7)') &
                      z(iz), f1(iz), g1(iz), g1(iz)-f1(iz)
        enddo
    endif

    !-----------------------------------------------------------------
    ! Test vertical velocity solve:

    alpha = f12
    do kx = box%lo(1), box%hi(1)
        do ky = box%lo(2), box%hi(2)
            g(:, ky, kx) = exp(alpha*z)
        enddo
    enddo

    call vertvel(g)

    gavg = zero
    do kx = box%lo(1), box%hi(1)
        do ky = box%lo(2), box%hi(2)
            rk = sqrt(k2l2(ky, kx))
            fac = one / (alpha ** 2 - k2l2(ky, kx))
            f(:, ky, kx) = fac * exp(alpha * z) !Particular solution
            fac = one / sinh(rk * extent(3))
            !Remove harmonic functions:
            f(:, ky, kx) = f(:, ky, kx) - fac*(f(0,  ky, kx)*sinh(rk*(upper(3) - z)) + &
                                               f(nz, ky, kx)*sinh(rk*(z - lower(3))))
            maxerr = maxval(abs(g(:, ky, kx) - f(:, ky, kx)))
            gavg = max(gavg, maxerr)

            if (verbose .and. (world%size == 1)) then
                if (maxerr > 1.d-6) write(*,'(2(1x,i3),1x,1p,e14.7)') kx, ky, maxerr
            endif
        enddo
    enddo

    maxerr = gavg

    call mpi_blocking_reduce(maxerr, MPI_MAX, world)

    if (verbose .and. (world%rank == world%root)) then
        write(*,*)
        write(*,*) ' Overall, the maximum abs error = ', maxerr
    endif

    if (world%rank == world%root) then
        call print_result_dp('Test zops', maxerr, atol=6.0e-7)
    endif

    call finalise_zops

    deallocate(f)
    deallocate(g)
    deallocate(f1)
    deallocate(g1)
    deallocate(z)

    call mpi_env_finalise

end program test_zops
