!==============================================================================
  ! Tests effect of either the 2/3 de-aliasing rule or the Hou & Li filter
  ! on an smooth field q(z) specified below.

  ! We also examine the effect of the filter on dq/dz.
!==============================================================================
program test_filter
    use stafft, only : initfft, dst, dct
    use deriv1d, only : init_deriv
    use constants, only : zero, one, two, pi
    implicit none

    ! Resolution in z:
    integer,parameter:: nz = 32

    logical        :: l_exist = .true.
    integer        :: file_num = 0
    character(512) :: fname

    ! Generic double precision numerical constants:
    double precision,parameter:: hpi = pi/two

    ! Domain height and grid length:
    double precision,parameter:: Lz = pi, dz = Lz/dble(nz)
    double precision,parameter:: dzi = one/dz, hdzi = one/(two*dz)
    double precision,parameter:: zmin = -Lz/two, zmax = zmin + Lz

    ! Arrays:
    double precision:: z(0:nz), q(0:nz), qs(0:nz), filt(0:nz)
    double precision:: dq(0:nz), dqs(0:nz)

    ! FFT arrays:
    integer:: zfactors(5)
    double precision:: ztrig(2*nz), rkz(nz)

    ! Others:
    double precision:: eps, t, fac
    integer:: iz, kz, iopt

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

    write(*,*) ' Choose (1) 2/3 rule or (2) Hou & Li filter:'
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
    else
        fac = one/rkz(nz)
        filt(0) = one
        do kz = 1, nz
            filt(kz) = exp(-36.d0*(fac*rkz(kz))**36.d0)
        enddo
    endif

    ! Transform to spectral space:
    qs = q
    call dst(1, nz, qs(1:nz), ztrig, zfactors)

    ! Apply filter:
    qs = filt*qs

    ! Differentiate:
    dqs(0) = zero
    dqs(1:nz-1) = rkz(1:nz-1)*qs(1:nz-1)
    dqs(nz) = zero

    ! Transform back to physical space:
    call dst(1, nz, qs(1:nz), ztrig, zfactors)
    call dct(1, nz, dqs, ztrig, zfactors)

    write (fname, "(a17,i1,a4)") 'q_profile_filter_', file_num, '.asc'
    do while (l_exist)
        file_num = file_num + 1
        write (fname, "(a17,i1,a4)") 'q_profile_filter_', file_num, '.asc'
        inquire(file=fname, exist=l_exist)
    enddo
    open(80, file = fname, status = 'replace')
    if (iopt == 1) then
        write(80,*) '#', eps, '2/3 rule filter'
    else
        write(80,*) '#', eps, 'Hou and Li filter'
    endif
    do iz = 0, nz
        write(80,*) z(iz), qs(iz), q(iz), qs(iz) - q(iz)
    enddo
    close(80)

    write (fname, "(a18,i1,a4)") 'dq_profile_filter_', file_num, '.asc'
    open(80, file = fname, status = 'replace')
    do iz = 0, nz
        write(80,*) z(iz), dqs(iz), dq(iz), dqs(iz) - dq(iz)
    enddo
    close(80)
end program test_filter
