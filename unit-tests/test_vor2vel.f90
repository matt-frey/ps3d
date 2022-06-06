! =============================================================================
!                               Test vor2vel
!
!  This unit test checks the calculation of the velocity field using the
!  Beltrami flow:
!     u(x, y, z) = (k^2 + l^2)^(-1) * [k*m*sin(mz) - l*alpha*cos(m*z) * sin(k*x + l*y)]
!     v(x, y, z) = (k^2 + l^2)^(-1) * [l*m*sin(mz) + k*alpha*cos(m*z) * sin(k*x + l*y)]
!     w(x, y, z) = cos(m*z) * cos(k*x + l*y)
!  The vorticity of this flow is
!    xi(x, y, z) = alpha * u(x, y, z)
!   eta(x, y, z) = alpha * v(x, y, z)
!  zeta(x, y, z) = alpha * w(x, y, z)
! =============================================================================
program test_vor2vel
    use unit_test
    use constants, only : one, two, pi, f12, f34, three
    use parameters, only : lower, update_parameters, dx, nx, ny, nz, extent
    use fields
    use inversion_utils
    use inversion_mod, only : vor2vel, vor2vel_timer
    use timer
    implicit none

    double precision              :: error
    double precision, allocatable :: vel_ref(:, :, :, :)
    integer                       :: ix, iy, iz, ik, il, im
    double precision              :: x, y, z, alpha, fk2l2, k, l, m
    double precision              :: cosmz, sinmz, sinkxly, coskxly

    call register_timer('vorticity', vor2vel_timer)

    nx = 128
    ny = 128
    nz = 128

    lower  = -f12 * pi * (/one, one, one/)
    extent =  pi * (/one, one, one/)


    allocate(vel_ref(0:nz, 0:ny-1, 0:nx-1, 3))

    call update_parameters

    call field_default

    call init_inversion

    do ik = 1, 3
        k = dble(ik)
        do il = 1, 3
            l = dble(il)
            do im = 1, 3, 2
                m = dble(im)

                alpha = dsqrt(k ** 2 + l ** 2 + m ** 2)
                fk2l2 = one / dble(k ** 2 + l ** 2)

                do ix = 0, nx-1
                    x = lower(1) + ix * dx(1)
                    do iy = 0, ny-1
                        y = lower(2) + iy * dx(2)
                        do iz = -1, nz+1
                            z = lower(3) + iz * dx(3)

                            cosmz = dcos(m * z)
                            sinmz = dsin(m * z)
                            sinkxly = dsin(k * x + l * y)
                            coskxly = dcos(k * x + l * y)

                            ! velocity
                            vel_ref(iz, iy, ix, 1) = fk2l2 * (k * m * sinmz - l * alpha * cosmz) * sinkxly
                            vel_ref(iz, iy, ix, 2) = fk2l2 * (l * m * sinmz + k * alpha * cosmz) * sinkxly
                            vel_ref(iz, iy, ix, 3) = cosmz * coskxly

                            ! vorticity
                            vor(iz, iy, ix, 1) = alpha * vel_ref(iz, iy, ix, 1)
                            vor(iz, iy, ix, 2) = alpha * vel_ref(iz, iy, ix, 2)
                            vor(iz, iy, ix, 3) = alpha * vel_ref(iz, iy, ix, 3)

                        enddo
                    enddo
                enddo

                call vor2vel

                error = max(error, maxval(dabs(vel_ref - vel)))
            enddo
        enddo
    enddo

    call print_result_dp('Test vor2vel', error, atol=5.0e-2)

    deallocate(vel_ref)

end program test_vor2vel
