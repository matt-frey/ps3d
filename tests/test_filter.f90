!==============================================================================
  ! Tests effect of either the 2/3 de-aliasing rule or the Hou & Li filter
  ! on an smooth field q(z) specified below.

  ! We also examine the effect of the filter on dq/dz.
!==============================================================================
program test_filter
    use stafft, only : initfft, dst, dct
    use deriv1d, only : init_deriv
    use constants, only : zero, one, two, pi
    use cheby, only : cheb_fun, cheb_poly, init_cheby
    implicit none

    ! Resolution in z:
    integer :: nz

    logical        :: l_exist = .true.
    integer        :: file_num = 0
    character(512) :: fname

    ! Generic double precision numerical constants:
    double precision, parameter:: hpi = pi / two

    ! Domain height and grid length:
    double precision, parameter :: Lz = pi
    double precision, parameter :: zmin = -Lz/two, zmax = zmin + Lz
    double precision            :: dz, alpha, beta, rkmax, kmax, x

    ! Arrays:
    double precision, allocatable :: z(:), q(:), qs(:), filt(:), coeff(:)
    double precision, allocatable :: dq(:), dqs(:)
    double precision, allocatable :: d1z(:, :), d2z(:, :), zcheb(:) ! Chebyshev grid points

    ! FFT arrays:
    integer :: zfactors(5)
    double precision, allocatable :: ztrig(:), rkz(:)

    ! Others:
    double precision:: eps, fac, err_e, err_o
    integer:: iz, kz, iopt

    write(*,*) ' Enter number of cells:'
    read(*,*) nz

    allocate(z(0:nz))
    allocate(q(0:nz))
    allocate(qs(0:nz))
    allocate(filt(0:nz))
    allocate(dq(0:nz))
    allocate(dqs(0:nz))
    allocate(ztrig(2*nz))
    allocate(rkz(nz))
    allocate(coeff(0:nz))
    allocate(d1z(0:nz, 0:nz))
    allocate(d2z(0:nz, 0:nz))
    allocate(zcheb(0:nz))

    call init_cheby(nz, zcheb, d1z, d2z)

    dz = Lz / dble(nz)

    !-------------------------------------------------------------
    ! Initialise Fourier transform and wavenumbers:
    call initfft(nz, zfactors, ztrig)
    call init_deriv(nz, Lz, rkz)

    !-------------------------------------------------------------
    ! Set up functions:
    do iz = 0, nz
        z(iz) = zmin + dz*dble(iz)
        call random_number(fac)
        q(iz) = two * fac - one
    enddo
    write(*,*) ' Enter the amplitude of the disturbance:'
    read(*,*) eps
    q = (z - zmin) * (Lz**2 - (z - zmin)**2) * (one + eps*q(iz))

    ! Get derivative dq/dz spectrally:
    dq = q
    call dst(1, nz, dq(1:nz), ztrig, zfactors)
    dq(0) = zero
    dq(1:nz-1) = rkz(1:nz-1)*dq(1:nz-1)
    dq(nz) = zero
    call dct(1, nz, dq, ztrig, zfactors)

    write(*,*) ' Choose (1) 2/3 rule, (2) Hou & Li filter or (3) Dembenek filter:'
    read(*,*) iopt

    if (iopt == 1) then
        fac = (2.d0/3.d0)*rkz(nz)
        do kz = 0, nz
            if (rkz(kz) < fac) then
                filt(kz) = one
            else
                filt(kz) = zero
            endif
        enddo
    else if (iopt == 2) then
        fac = one/rkz(nz)
        filt(0) = one
        do kz = 1, nz
            filt(kz) = exp(-36.d0*(fac*rkz(kz))**36.d0)
        enddo
    else
        alpha = 27.0d0
        beta = 5.3d0
        kmax = 2.0d0 / 3.0d0
        rkmax = kmax * dble(nz)
        do kz = 0, nz
            x = dble(kz) / rkmax
            filt(kz) = dexp(-alpha * x**beta)
        enddo
    endif

    qs = q

    if (iopt == 3) then
        call cheb_poly(nz, qs, coeff)
        coeff = filt * coeff
        err_e = coeff(0)
        err_o = coeff(1)
        do iz = 1, nz/2
            err_e  = err_e +  coeff(2*iz)
        enddo

        do iz = 1, nz/2-1
            err_o  = err_o +  coeff(2*iz+1)
        enddo
        coeff(0) = coeff(0) - err_e
        coeff(1) = coeff(1) - err_o
        call cheb_fun(nz, coeff, qs)

        dqs = matmul(d1z, qs)
    else
        ! Transform to spectral space:
        call dst(1, nz, qs(1:nz), ztrig, zfactors)

        ! Apply filter:
        qs = filt * qs

        ! Differentiate:
        dqs(0) = zero
        dqs(1:nz-1) = rkz(1:nz-1)*qs(1:nz-1)
        dqs(nz) = zero

        ! Transform back to physical space:
        call dst(1, nz, qs(1:nz), ztrig, zfactors)
        call dct(1, nz, dqs, ztrig, zfactors)
    endif

    write (fname, "(a17,i1,a4)") 'q_profile_filter_', file_num, '.asc'
    do while (l_exist)
        file_num = file_num + 1
        write (fname, "(a17,i1,a4)") 'q_profile_filter_', file_num, '.asc'
        inquire(file=fname, exist=l_exist)
    enddo
    open(80, file = fname, status = 'replace')
    if (iopt == 1) then
        write(80,*) '#', nz, eps, '2/3 rule filter'
    else if (iopt == 2) then
        write(80,*) '#', nz, eps, 'Hou and Li filter'
    else
        write(80,*) '#', nz, eps, 'Dembenek filter'
    endif
    do iz = 0, nz
        write(80,*) z(iz), qs(iz), q(iz), qs(iz) - q(iz)
    enddo
    close(80)

    write (fname, "(a18,i1,a4)") 'dq_profile_filter_', file_num, '.asc'
    open(80, file = fname, status = 'replace')
    if (iopt == 1) then
        write(80,*) '#', nz, eps, '2/3 rule filter'
    else if (iopt == 2) then
        write(80,*) '#', nz, eps, 'Hou and Li filter'
    else
        write(80,*) '#', nz, eps, 'Dembenek filter'
    endif
    do iz = 0, nz
        write(80,*) z(iz), dqs(iz), dq(iz), dqs(iz) - dq(iz)
    enddo
    close(80)

    deallocate(z)
    deallocate(q)
    deallocate(qs)
    deallocate(filt)
    deallocate(dq)
    deallocate(dqs)
    deallocate(ztrig)
    deallocate(rkz)
    deallocate(coeff)
    deallocate(d1z, d2z, zcheb)
end program test_filter
