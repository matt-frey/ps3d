! =============================================================================
!                               Test diffz
!
!  This unit test checks the calculation of the vertical differentiation of the
!  function
!       f(x, y, z) = z
! =============================================================================
program test_diffz_1
    use unit_test
    use constants, only : zero, one, two, pi, f12
    use parameters, only : lower, update_parameters, nx, ny, nz, extent, hl, center
    use inversion_utils
    use mpi_timer
    use mpi_environment
    use mpi_layout
    use mpi_collectives
    implicit none

    double precision              :: error
    double precision, allocatable :: dfdz_ref(:, :, :), dfdz(:, :, :)
    double precision, allocatable :: fp(:, :, :)
    integer                       :: iz
    double precision              :: z

    call mpi_env_initialise

    nx = 32
    ny = 32
    nz = 32

    lower  = (/zero, zero, zero/)
    extent = (/one, one, one/)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    allocate(fp(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
    allocate(dfdz(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
    allocate(dfdz_ref(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

    dfdz = zero

    call update_parameters

    call init_inversion

    do iz = 0, nz
        z = center(3) - hl(3) * zcheb(iz)
        fp(iz, :, :) = z
        dfdz_ref(iz, :, :) = one
    enddo

    ! calculate z-derivative (dfdz)
    call zderiv(fp, dfdz)

    error = maxval(dabs(dfdz_ref - dfdz))

    call mpi_blocking_reduce(error, MPI_MAX, world)

    if (world%rank == world%root) then
        call print_result_dp('Test diffz', error, atol=3.0e-14)
    endif

    deallocate(fp)
    deallocate(dfdz)
    deallocate(dfdz_ref)

    call mpi_env_finalise

end program test_diffz_1
