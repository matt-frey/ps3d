! =============================================================================
!                               Test vor2vel
!
!  This unit test checks the calculation of the velocity field using the
!  Nonseperarable flow:
!     u(x, y, z) = sin(kx + ly + mz) = S()
!     v(x, y, z) = sin(kx + ly + mz) = S()
!     w(x, y, z) = -(k+l)/m sin(kx + ly + mz) = alpha S()
!  The vorticity of this flow is
!    xi(x, y, z) = -[ l(k+l)/m + m ] C()
!   eta(x, y, z) =  [ k(k+l)/m + m ] C()
!  zeta(x, y, z) =  [ (k-l)] C()
! =============================================================================
program test_vor2vel_6
    use unit_test
    use constants, only : one, two, pi, f12, f14, f32, f23, f34, two, three, four, five, six, ten
    use parameters, only : lower, update_parameters, nx, ny, nz, extent
    use fields
    use inversion_mod, only : vor2vel, vor2vel_timer, source, vorticity_tendency, vtend_timer
    use mpi_timer
    use mpi_environment
    use mpi_layout
    use mpi_collectives, only : mpi_blocking_reduce
    use sta3dfft, only : fftxyp2s, fftxys2p
    use model, only : layout, create_model!, filter

    implicit none

    call mpi_env_initialise

    call register_timer('vorticity', vor2vel_timer)
    call register_timer('tendency', vtend_timer)

    PRINT *, 'Enter the number of grid points in each direction:'
    READ(*,*) nx, ny, nz

    lower  = -f12 * pi * (/one, one, one/)
    extent =  pi * (/one, one, one/)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    call field_default

    write(*,*) "Uniform Grid"
    call run_test("uniform",nx)

    write(*,*) ""
    write(*,*) ""

    write(*,*) "Chebyshev Grid"
    call run_test("chebyshev",nx)

    call mpi_env_finalise

contains

    subroutine run_test(grid_type,nx)
        character(*), intent(in)      :: grid_type
        integer, intent(in)           :: nx
        double precision              :: error,ex,ey,ez
        integer                       :: ix, iy, iz, npx,npy
        double precision              :: k, l, m
        double precision              :: uu,vv, xx, yy, zz
        double precision              :: u, v, w, xi, eta, zeta
        double precision              :: u_x, u_y, u_z
        double precision              :: v_x, v_y, v_z
        double precision              :: w_x, w_y, w_z
        double precision              :: xi_x, xi_y, xi_z
        double precision              :: eta_x, eta_y, eta_z
        double precision              :: zeta_x, zeta_y, zeta_z
        double precision              :: sinkxly, coskxly, sincosm, sincosp, cossinm, cossinp
        double precision              :: cosx, cosy, sinx, siny, sinz, cosz
        double precision, allocatable :: vel_ref(:, :, :, :)
        double precision, allocatable :: tend_ref(:, :, :, :)
        double precision, allocatable :: x(:), y(:), z(:)
        character(len=100)       :: fname


        fname = TRIM(ADJUSTL(grid_type))
        fname = trim(fname) // '.dat'
        print *, fname

!        i = INDEX(fname, '(')
!IF (i > 0) THEN
!fname = fname(i+1:LEN(fname)-1)
!        ELSE
!PRINT *, 'Error: Invalid input string format'
!RETURN
!        END IF

        open(unit=10,file=fname,status='replace',action='write')
!         open(unit=10,file=fname,status='old',position='append',action='write')

        k = two
        l = two
        m = two

        npx = nx/2
        npy = nx/2

        call create_model(grid_type)

        uu =  1.0/ten
        vv =  -1.0/five

        allocate(vel_ref(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1), 3))
        allocate(tend_ref(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1), 3))
        allocate(x(0:nx-1), y(0:ny-1), z(0:nz))


        !!!!!!!!   Newly Perturbed Beltrami !!!!!! (2/13/25)
        x = layout%get_x_axis()
        y = layout%get_y_axis()
        z = layout%get_z_axis()
        do ix = box%lo(1), box%hi(1)
            do iy = box%lo(2), box%hi(2)
                do iz = 0, nz
                    xx = x(ix)
                    yy = y(iy)
                    zz = z(iz)

                    sinz = sin(zz)
                    cosz = cos(zz)
                    sinkxly = sin(k * x(ix) + l * y(iy) )
                    coskxly = cos(k * x(ix) + l * y(iy) )
                    sincosm = sinz - three * cosz
                    sincosp = sinz + three * cosz
                    cossinm = cosz - three * sinz
                    cossinp = cosz + three * sinz
                    cosx = cos(two*xx)
                    cosy = cos(two*yy)
                    sinx = sin(two*xx)
                    siny = sin(two*yy)



                    ! velocity
                    u = f14 * sincosm * sinkxly + uu * cosx * sinz
                    v = f14 * sincosp * sinkxly + vv * cosy * sinz
                    w = cosz * coskxly - two * cosz * (uu * sinx + vv * siny)

                    vel_ref(iz, iy, ix, 1) = u
                    vel_ref(iz, iy, ix, 2) = v
                    vel_ref(iz, iy, ix, 3) = w

                    ! vorticity
                    xi   = f34 * sincosm * sinkxly - five * vv * cosy * cosz
                    eta  = f34 * sincosp * sinkxly + five *  uu * cosx * cosz
                    zeta = three * cosz * coskxly

                    vor(iz,iy,ix,1) = xi
                    vor(iz,iy,ix,2) = eta
                    vor(iz,iy,ix,3) = zeta

                    u_x = f12 * sincosm * coskxly - two * uu * sinx * sinz
                    u_y = f12 * sincosm * coskxly
                    u_z = f14 * (cosz + three * sinz) * sinkxly + uu * cosx * cosz

                    v_x = f12 * sincosp * coskxly
                    v_y = f12 * sincosp * coskxly - two * vv * siny * sinz
                    v_z = f14 * cossinm * sinkxly + vv * cosy * cosz

                    w_x = - two * cosz * sinkxly - four * uu * cosz * cosx
                    w_y = - two * cosz * sinkxly - four * vv * cosz * cosy
                    w_z = - sinz * coskxly + two * sinz * (uu * sinx + vv * siny)


                    xi_x = f32 * sincosm * coskxly
                    xi_y = f32 * sincosm * coskxly + ten * vv * siny * cosz
                    xi_z = f34 * cossinp * sinkxly + five * vv * cosy * sinz

                    eta_x = f32 * sincosp * coskxly - ten * uu * sinx * cosz
                    eta_y = f32 * sincosp * coskxly
                    eta_z = f34 * cossinm * sinkxly - five * uu * cosx * sinz

                    zeta_x = - six * cosz * sinkxly
                    zeta_y = zeta_x
                    zeta_z = - three * sinz * coskxly

                    ! tendency reference
                    tend_ref(iz,iy,ix,1) = xi * u_x + eta * u_y + zeta * u_z &
                                         - u * xi_x - v * xi_y - w * xi_z

                    tend_ref(iz,iy,ix,2) = xi * v_x + eta * v_y + zeta * v_z &
                                         - u * eta_x - v * eta_y - w * eta_z

                    tend_ref(iz,iy,ix,3) = xi * w_x + eta * w_y + zeta * w_z &
                                         - u * zeta_x - v * zeta_y - w * zeta_z

                enddo
            enddo
        enddo

        call fftxyp2s(vor(:, :, :, 1), svor(:, :, :, 1))
        call fftxyp2s(vor(:, :, :, 2), svor(:, :, :, 2))
        call fftxyp2s(vor(:, :, :, 3), svor(:, :, :, 3))

        call vor2vel

        error = maxval(abs(vel_ref - vel))
        write(*, '(A, I0, A, I0)') "Grid: nx = ", nx, ", nz = ", nz
                print *, error

        write(*, '(A, F8.3, A, F8.3)') "Sample: x0 = ", x(npx), ", y0 = ", y(npy)


        do iz = 0,nz
        ex = vel(iz,npx,npy,1)-vel_ref(iz,npx,npy,1)
        ey = vel(iz,npx,npy,2)-vel_ref(iz,npx,npy,2)
        ez = vel(iz,npx,npy,3)-vel_ref(iz,npx,npy,3)
!         write(10,'(f8.4,4x,e12.5,4x,e12.5,4x,e12.5)')z(iz), ex,ey,ez
        enddo

!         write(10,*)''
!         do iz = 0,nz
!         write(10,'(f8.4,4x,f15.12,4x,f15.12)')z(iz),  vel(iz,npx,npy,1), vel(iz,npx,npy,3)
!         enddo

        call source

        call fftxys2p(svorts(:, :, :, 1), vor(:, :, :, 1))
        call fftxys2p(svorts(:, :, :, 2), vor(:, :, :, 2))
        call fftxys2p(svorts(:, :, :, 3), vor(:, :, :, 3))

!         write(10,*)''
!         do iz = 0,nz
!         write(10,'(f8.4,4x,f15.8,4x,f15.8,2x,f15.8)')z(iz),vor(iz,npx,npy,1),vor(iz,npx,npy,2),vor(iz,npx,npy,3)
!         enddo

        ex = maxval(abs(tend_ref(0, :, :, 1) - vor(0, :, :, 1)))
        ey = maxval(abs(tend_ref(0, :, :, 2) - vor(0, :, :, 2)))
        ez = maxval(abs(tend_ref(0, :, :, 3) - vor(0, :, :, 3)))
        u = maxval(abs(tend_ref(nz/2, :, :, 1) - vor(nz/2, :, :, 1)))
        v = maxval(abs(tend_ref(nz/2, :, :, 2) - vor(nz/2, :, :, 2)))
        w = maxval(abs(tend_ref(nz/2, :, :, 3) - vor(nz/2, :, :, 3)))
        xi = maxval(abs(tend_ref(nz, :, :, 1) - vor(nz, :, :, 1)))
        eta = maxval(abs(tend_ref(nz, :, :, 2) - vor(nz, :, :, 2)))
        zeta = maxval(abs(tend_ref(nz, :, :, 3) - vor(nz, :, :, 3)))
        xi_x = maxval(abs(tend_ref(:, :, :, 1) - vor(:, :, :, 1)))
        xi_y = maxval(abs(tend_ref(:, :, :, 2) - vor(:, :, :, 2)))
        xi_z = maxval(abs(tend_ref(:, :, :, 3) - vor(:, :, :, 3)))
        write(10, *) nz, ex, ey, ez, u, v, w, xi, eta, zeta, xi_x, xi_y, xi_z

!        call fftxyp2s(vor(:, :, :, 1), svor(:, :, :, 1))
!        call fftxyp2s(vor(:, :, :, 2), svor(:, :, :, 2))
!        call fftxyp2s(vor(:, :, :, 3), svor(:, :, :, 3))

!        call filter%apply(svor(:, :, :, 1))
!        call filter%apply(svor(:, :, :, 2))
!        call filter%apply(svor(:, :, :, 3))

!        call fftxys2p(svor(:, :, :, 1), vor(:, :, :, 1))
!        call fftxys2p(svor(:, :, :, 2), vor(:, :, :, 2))
!        call fftxys2p(svor(:, :, :, 3), vor(:, :, :, 3))

!        write(10,*)''
!        do iz = 0,nz
!        write(10,'(f8.4,4x,f15.8,4x,f15.8,2x,f15.8)')z(iz),vor(iz,npy,npx,1),vor(iz,npy,npx,2),vor(iz,npy,npx,3)
!        enddo


        close(10)

        deallocate(x, y, z)
        deallocate(vel_ref)
        deallocate(tend_ref)

        call mpi_blocking_reduce(error, MPI_MAX, world)

        if (world%rank == world%root) then
            call print_result_dp('Test vor2vel ' // grid_type, error, atol=1.0e-14)
        endif
    end subroutine

end program test_vor2vel_6
