!==============================================================================
! Tests effect of hyperdiffusion on an initial field q(z,0) = q_0(z)
! in [-pi/2,pi/2],  Here, q_t = -nu (-d^2/dz^2)^p(q). We use the exact
! solution in spectral space q_m(0) * exp(-nu * m^{2p} * t).
!
! We also record dq/dz.
!==============================================================================
program test_hyper
    use stafft, only : initfft, dst, dct
    use deriv1d, only : init_deriv
    use constants, only : zero, one, two, pi
    implicit none

    ! Resolution in z:
    integer :: nz

    logical        :: l_exist = .true.
    integer        :: file_num = 0
    character(512) :: fname

    ! Generic double precision numerical constants:
    double precision,parameter:: hpi = pi/two

    ! Domain height and grid length:
    double precision, parameter :: Lz = pi
    double precision, parameter :: zmin = -Lz/two, zmax = zmin + Lz
    double precision            :: dz

    ! Arrays:
    double precision, allocatable :: z(:), q(:), qs(:), diss(:)
    double precision, allocatable :: dq(:), dqs(:), q0(:), dq0(:)


    ! FFT arrays:
    integer :: zfactors(5)
    double precision, allocatable :: ztrig(:), rkz(:)

    ! Others:
    double precision:: eps, t, fac
    integer:: iz, nnu

    write(*,*) ' Enter number of cells:'
    read(*,*) nz

    allocate(z(0:nz))
    allocate(q(0:nz))
    allocate(qs(0:nz))
    allocate(diss(0:nz))
    allocate(dq(0:nz))
    allocate(dqs(0:nz))
    allocate(q0(0:nz))
    allocate(dq0(0:nz))
    allocate(ztrig(2*nz))
    allocate(rkz(nz))

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
    q0 = q

    write(*,*) ' Hyperviscosity order (use 1 for molecular viscosity)?'
    read(*,*) nnu

    ! Get derivative dq0/dz spectrally:
    qs = q
    call dst(1, nz, qs(1:nz), ztrig, zfactors)
    dqs(0) = zero
    dqs(1:nz-1) = rkz(1:nz-1)*qs(1:nz-1)
    dqs(nz) = zero
    call dct(1, nz, dqs, ztrig, zfactors)
    dq0 = dqs

    write(*,*) ' We take nu = 1/kz_max^{2p}.  Enter t:'
    read(*,*) t
    fac = one/rkz(nz)
    diss(0) = zero
    diss(1:nz-1) = exp(-t*(fac*rkz(1:nz-1))**(2*nnu))
    diss(nz) = zero

    ! Transform to spectral space:
    qs = q
    call dst(1, nz, qs(1:nz), ztrig, zfactors)

    ! Apply damping:
    qs = diss*qs

    ! Differentiate:
    dqs(0) = zero
    dqs(1:nz-1) = rkz(1:nz-1)*qs(1:nz-1)
    dqs(nz) = zero

    ! Transform back to physical space:
    call dst(1, nz, qs(1:nz), ztrig, zfactors)
    call dct(1, nz, dqs, ztrig, zfactors)
    q = qs
    dq = dqs

    write (fname, "(a16,i1,a4)") 'q_profile_hyper_', file_num, '.asc'
    do while (l_exist)
        file_num = file_num + 1
        write (fname, "(a16,i1,a4)") 'q_profile_hyper_', file_num, '.asc'
        inquire(file=fname, exist=l_exist)
    enddo
    open(80, file = fname, status = 'replace')
    write(80,*) '#', nz, eps, nnu
    do iz = 0, nz
        write(80,*) z(iz), q(iz), q0(iz), q(iz) - q0(iz)
    enddo
    close(80)

    write (fname, "(a17,i1,a4)") 'dq_profile_hyper_', file_num, '.asc'
    open(80, file = fname, status = 'replace')
    write(80,*) '#', nz, eps, nnu
    do iz = 0, nz
        write(80,*) z(iz), dq(iz), dq0(iz), dq(iz) - dq0(iz)
    enddo
    close(80)

    deallocate(z)
    deallocate(q)
    deallocate(qs)
    deallocate(diss)
    deallocate(dq)
    deallocate(dqs)
    deallocate(q0)
    deallocate(dq0)
    deallocate(ztrig)
    deallocate(rkz)
end program test_hyper
