! =============================================================================
!                       Test convergence of vor2vel
! =============================================================================
program test_vor2vel
    use constants, only : f12, f13, one, two, three, six, pi, twopi, f14
    use parameters, only : lower, update_parameters, dx, nx, ny, nz, extent, upper
    use fields
    use inversion_utils
    use inversion_mod, only : vor2vel, vor2vel_timer
    use mpi_timer
    use mpi_environment
    use mpi_layout
    use mpi_utils, only : mpi_stop
    use sta3dfft, only : fftxyp2s
    implicit none

    double precision              :: emax, erms
    double precision, allocatable :: vel_ref(:, :, :, :)
    double precision, allocatable :: mag(:, :, :)
    integer                       :: ix, iy, iz, casenum
    double precision              :: x, y, z, k, l, m, alpha, klsq, coskx, cosly, sinkx, sinly
    double precision              :: f, dfdz, d2fdz2, coskxly, sinkxly, cosmz, sinmz, fk2l2

    call mpi_env_initialise

    if (world%size > 1) then
        call mpi_stop("This program only works with 1 MPI rank!")
    endif

    call register_timer('vorticity', vor2vel_timer)

    call get_input_arguments

    call mpi_layout_init(lower, extent, nx, ny, nz)

    allocate(vel_ref(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1), 3))
    allocate(mag(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
    call field_default

    call init_inversion

    if (casenum == 1) then
        k = two
        l = two
        m = one

        alpha = dsqrt(k ** 2 + l ** 2 + m ** 2)
        fk2l2 = one / dble(k ** 2 + l ** 2)

        do ix = box%lo(1), box%hi(1)
            x = lower(1) + ix * dx(1)
            do iy = box%lo(2), box%hi(2)
                y = lower(2) + iy * dx(2)
                do iz = 0, nz
                    cosmz = dcos(m * zg(iz))
                    sinmz = dsin(m * zg(iz))
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

        call invert

        call write_result('test_vor2vel_1.asc')

    else if (casenum == 2) then

        ! test_vor2vel_2
        l = twopi
        k = two * l

        klsq = k ** 2 - l ** 2

        do ix = box%lo(1), box%hi(1)
            x = lower(1) + ix * dx(1)
            do iy = box%lo(2), box%hi(2)
                y = lower(2) + iy * dx(2)
                do iz = 0, nz
                    z = zg(iz)

                    f = two * z - z ** 2 - z ** 3
                    dfdz = two - two * z - three * z ** 2
                    d2fdz2 = -two - six * z
                    sinkx = dsin(k * x)
                    coskx = dcos(k * x)
                    sinly = dsin(l * y)
                    cosly = dcos(l * y)

                    ! velocity
                    vel_ref(iz, iy, ix, 1) =  k * dfdz * coskx * sinly
                    vel_ref(iz, iy, ix, 2) = -l * dfdz * sinkx * cosly
                    vel_ref(iz, iy, ix, 3) =  klsq * f * sinkx * sinly

                    ! vorticity
                    vor(iz, iy, ix, 1) = l * (klsq * f + d2fdz2) * sinkx * cosly
                    vor(iz, iy, ix, 2) = k * (d2fdz2 - klsq * f) * coskx * sinly
                    vor(iz, iy, ix, 3) = - two * k * l * dfdz * coskx * cosly
                enddo
            enddo
        enddo

        call invert

        call write_result('test_vor2vel_2.asc')

        !----------------------------------------------------------------------

        ! linear -- test_vor2vel_3
        do iz = 0, nz
            z = zg(iz)

            ! velocity
            vel_ref(iz, :, :, 1) = z - f12 * (lower(3) + upper(3))
            vel_ref(iz, :, :, 2) = zero
            vel_ref(iz, :, :, 3) = zero

            ! vorticity
            vor(iz, :, :, 1) = zero
            vor(iz, :, :, 2) = one
            vor(iz, :, :, 3) = zero
        enddo

        call invert

        call write_result('test_vor2vel_3.asc')

        !----------------------------------------------------------------------

        ! quadratic -- test_vor2vel_4
        do iz = 0, nz
            z = zg(iz)

            ! velocity
            vel_ref(iz, :, :, 1) = z ** 2 - f13
            vel_ref(iz, :, :, 2) = zero
            vel_ref(iz, :, :, 3) = zero

            ! vorticity
            vor(iz, :, :, 1) = zero
            vor(iz, :, :, 2) = two * z
            vor(iz, :, :, 3) = zero
        enddo

        call invert

        call write_result('test_vor2vel_4.asc')

        !----------------------------------------------------------------------

        ! cubic -- test_vor2vel_5
        do iz = 0, nz
            z = zg(iz)

            ! velocity
            vel_ref(iz, :, :, 1) = z ** 3 - f14
            vel_ref(iz, :, :, 2) = zero
            vel_ref(iz, :, :, 3) = zero

            ! vorticity
            vor(iz, :, :, 1) = zero
            vor(iz, :, :, 2) = three * z ** 2
            vor(iz, :, :, 3) = zero
        enddo

        call invert

        call write_result('test_vor2vel_5.asc')

    endif

    deallocate(vel_ref)
    deallocate(mag)

    call mpi_env_finalise

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine get_input_arguments
            write(*, *) ' Enter the number of cells:'
            read(*, *) nx
            ny = nx
            nz = nx

            write(*, *) ' Enter the case to study:'
            write(*, *) ' 1 - Beltrami'
            write(*, *) ' 2 - linear/quadratic/cubic'
            read(*, *) casenum

            if (casenum == 1) then
                lower  = -f12 * pi * (/one, one, one/)
                extent =  pi * (/one, one, one/)
            else if (casenum == 2) then
                lower  = (/-f12, -f12, zero/)
                extent = (/one, one, one/)
            else
                write(*, *) ' No such case. Exiting.'
                stop
            endif

            call update_parameters

        end subroutine get_input_arguments

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine invert
            call fftxyp2s(vor(:, :, :, 1), svor(:, :, :, 1))
            call fftxyp2s(vor(:, :, :, 2), svor(:, :, :, 2))
            call fftxyp2s(vor(:, :, :, 3), svor(:, :, :, 3))

            call vor2vel
        end subroutine invert

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine magnitude
            vel_ref = vel_ref - vel

            mag = dsqrt(vel_ref(:, :, :, 1) ** 2 &
                      + vel_ref(:, :, :, 2) ** 2 &
                      + vel_ref(:, :, :, 3) ** 2)
        end subroutine magnitude

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine get_emax
            emax = maxval(mag)
        end subroutine get_emax

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine get_erms
            erms = f12 * sum(mag(0, :, :) ** 2 + mag(nz, :, :) ** 2) + sum(mag(1:nz-1, :, :) ** 2)
            erms = dsqrt(erms / dble(nz*nx*ny))
        end subroutine get_erms

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine write_result(fname)
            character(*), intent(in) :: fname
            logical                  :: l_exist

            call magnitude
            call get_emax
            call get_erms

            inquire(file=fname, exist=l_exist)

            if (l_exist) then
                open(88,file=fname, status='old', position='append')
            else
                open(88,file=fname, status='new')
            endif

            write(88,'(i4,2(2x,1p,e14.7))') nz, emax, erms
            close(88)
        end subroutine write_result

end program test_vor2vel
