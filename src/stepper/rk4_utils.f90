module rk4_utils
    use fields, only : velgradg, tbuoyg, vortg
    use constants, only : zero, one, two, f12
    use parameters, only : nx, ny, nz, dxi, vcell
    use jacobi, only : jacobi_eigenvalues
#ifdef ENABLE_VERBOSE
    use options, only : output
#endif

    implicit none

    contains

        ! Estimate a suitable time step based on the velocity strain
        ! and buoyancy gradient.
        ! @param[in] t is the time
        ! @returns the time step
        function get_time_step(t) result(dt)
            use options, only : time
            double precision, intent(in) :: t
            double precision             :: dt
            double precision             :: gmax, bmax, strain(3, 3), D(3)
            double precision             :: gradb(0:nz, 0:ny-1, 0:nx-1)
            double precision             :: db2(0:nz, 0:ny-1, 0:nx-1)
            integer                      :: ix, iy, iz
#if ENABLE_VERBOSE
            logical                      :: exists = .false.
            character(:), allocatable    :: fname
#endif

            !
            ! velocity strain
            !
            ! find largest stretch -- this corresponds to largest
            ! eigenvalue over all local symmetrised strain matrices.
            gmax = epsilon(gmax)
            do ix = 0, nx-1
                do iy = 0, ny-1
                    do iz = 0, nz
                        ! get local symmetrised strain matrix, i.e. 1/ 2 * (S + S^T)
                        ! where
                        !     /u_x u_y u_z\
                        ! S = |v_x v_y v_z|
                        !     \w_x w_y w_z/
                        ! with u_* = du/d* (also derivatives of v and w).
                        ! The derivatives dv/dx, du/dz, dv/dz and dw/dz are calculated
                        ! with vorticity or the assumption of incompressibility
                        ! (du/dx + dv/dy + dw/dz = 0):
                        !    dv/dx = \omegaz + du/dy
                        !    du/dz = \omegay + dw/dx
                        !    dv/dz = dw/dy - \omegax
                        !    dw/dz = - (du/dx + dv/dy)
                        !
                        !                         /  2 * u_x  u_y + v_x u_z + w_x\
                        ! 1/2 * (S + S^T) = 1/2 * |u_y + v_x   2 * v_y  v_z + w_y|
                        !                         \u_z + w_x  v_z + w_y   2 * w_z/
                        !
                        ! S11 = du/dx
                        ! S12 = 1/2 * (du/dy + dv/dx) = 1/2 * (2 * du/dy + \omegaz) = du/dy + 1/2 * \omegaz
                        ! S13 = 1/2 * (du/dz + dw/dx) = 1/2 * (\omegay + 2 * dw/dx) = 1/2 * \omegay + dw/dx
                        ! S22 = dv/dy
                        ! S23 = 1/2 * (dv/dz + dw/dy) = 1/2 * (2 * dw/dy - \omegax) = dw/dy - 1/2 * \omegax
                        ! S33 = dw/dz = - (du/dx + dv/dy)
                        strain(1, 1) = velgradg(iz, iy, ix, 1)                              ! S11
                        strain(1, 2) = velgradg(iz, iy, ix, 2) + f12 * vortg(iz, iy, ix, 3) ! S12
                        strain(1, 3) = velgradg(iz, iy, ix, 4) + f12 * vortg(iz, iy, ix, 2) ! S13
                        strain(2, 2) = velgradg(iz, iy, ix, 3)                              ! S22
                        strain(2, 3) = velgradg(iz, iy, ix, 5) - f12 * vortg(iz, iy, ix, 1) ! S23
                        strain(3, 3) = -(velgradg(iz, iy, ix, 1) + velgradg(iz, iy, ix, 3)) ! S33

                        ! calculate its eigenvalues. The Jacobi solver
                        ! requires the upper triangular matrix only.
                        call jacobi_eigenvalues(strain, D)

                        ! we must take the largest eigenvalue in magnitude (absolute value)
                        gmax = max(gmax, maxval(abs(D)))
                    enddo
                enddo
            enddo

            !
            ! buoyancy gradient (central difference)
            !

            ! db/dx inner part
            gradb(:, :, 1:nx-2) = f12 * dxi(1) * (tbuoyg(0:nz, :, 2:nx-1) - tbuoyg(0:nz, :, 0:nx-3))

            ! db/dx boundary grid points (make use of periodicity)
            gradb(:, :, 0)    = f12 * dxi(1) * (tbuoyg(0:nz, :, 1) - tbuoyg(0:nz, :, nx-1))
            gradb(:, :, nx-1) = f12 * dxi(1) * (tbuoyg(0:nz, :, 0) - tbuoyg(0:nz, :, nx-2))

            db2 = gradb ** 2

            ! db/dy inner part
            gradb(:, 1:ny-2, :) = f12 * dxi(2) * (tbuoyg(0:nz, 2:ny-1, :) - tbuoyg(0:nz, 0:ny-3, :))

            ! db/dy boundary grid points (make use of periodicity)
            gradb(:, 0, :)    = f12 * dxi(2) * (tbuoyg(0:nz, 1, :) - tbuoyg(0:nz, ny-1, :))
            gradb(:, ny-1, :) = f12 * dxi(2) * (tbuoyg(0:nz, 0, :) - tbuoyg(0:nz, ny-2, :))

            db2 = db2 + gradb ** 2

            ! db/dz
            gradb = f12 * dxi(3) * (tbuoyg(1:nz+1, :, :) - tbuoyg(-1:nz-1, :, :))

            bmax = dsqrt(dsqrt(maxval(db2 + gradb ** 2)))
            bmax = max(epsilon(bmax), bmax)

            dt = min(time%alpha / gmax, time%alpha / bmax)
#ifdef ENABLE_VERBOSE
            fname = trim(output%basename) // '_alpha_time_step.asc'
            inquire(file=fname, exist=exists)
            ! 23 August
            ! https://stackoverflow.com/questions/15526203/single-command-to-open-a-file-or-create-it-and-the-append-data
            if ((t /= zero) .and. exists) then
                open(unit=1235, file=fname, status='old', position='append')
            else
                open(unit=1235, file=fname, status='replace')
                write(1235, *) '  # time (s)                \alpha_s/\gamma_{max}     \alpha_b/N_{max}'
            endif

            write(1235, *) t, time%alpha / gmax, time%alpha / bmax

            close(1235)
#endif
            if (time%precise_stop .and. (t + dt > time%limit)) then
                dt = time%limit - t
            endif

            if (dt <= zero) then
                print "(a10, f0.2, a6)", "Time step ", dt, " <= 0!"
                stop
            endif
        end function get_time_step

end module rk4_utils
