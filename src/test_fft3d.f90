program test_fft3d
    use constants
    use inversion_utils, only : init_inversion, fftxys2p, diffz
    use sta3dfft, only : fftczp2s, fftfs2ss
    use parameters, only : update_parameters, dx, nx, ny, nz, lower, extent
    implicit none

    double precision, allocatable :: fp(:, :, :), fs(:, :, :), ref(:, :, :), ss(:, :, :), ds(:, :, :)
    integer            :: ix, iy, iz
    double precision   :: x, y, z, k, l, m

    nx = 32
    ny = 32
    nz = 32
    lower = (/-pi, -pi, -pi/)
    extent = (/twopi, twopi, twopi/)

    allocate(fp(0:nz, 0:ny-1, 0:nx-1))
    allocate(fs(0:nz, 0:nx-1, 0:ny-1))
    allocate(ss(0:nz, 0:nx-1, 0:ny-1))
    allocate(ds(0:nz, 0:nx-1, 0:ny-1))
    allocate(ref(0:nz, 0:nx-1, 0:ny-1))


    call update_parameters

    call init_inversion(bbdif=zero, nnu=3, prediss=10.0d0)


    fp = zero

    k = two * pi / extent(1)
    l = two * pi / extent(2)
    m =       pi / extent(3)

    ! setup test field
    do ix = 0, nx-1
        x = lower(1) + dble(ix) * dx(1)
        do iy = 0, ny-1
            y = lower(2) + dble(iy) * dx(2)
            do iz = 0, nz
                z = lower(3) + dble(iz) * dx(3)
                fp(iz, iy, ix) = dcos(k * x) * dcos(l * y) * dcos(m * z)
                ref(iz, ix, iy) = -m * dcos(k * x) * dcos(l * y) * dsin(m * z)
            enddo
        enddo
    enddo

    call fftczp2s(fp, fs)

    call fftfs2ss(fs, ss)

    call diffz(ss, ds)

    call fftxys2p(ds, fp)

    print *, "max. abs. error = ", maxval(abs(ref - fp))
    print *, "min. abs. error = ", minval(abs(ref - fp))

    deallocate(fp)
    deallocate(fs)
    deallocate(ss)
    deallocate(ds)
    deallocate(ref)

end program test_fft3d
