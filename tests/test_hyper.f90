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
    integer,parameter:: nz = 32

    ! Generic double precision numerical constants:
    double precision,parameter:: hpi = pi/two

    ! Domain height and grid length:
    double precision,parameter:: Lz = pi, dz = Lz/dble(nz)
    double precision,parameter:: dzi = one/dz, hdzi = one/(two*dz)
    double precision,parameter:: zmin = -Lz/two, zmax = zmin + Lz

    ! Arrays:
    double precision:: z(0:nz), q(0:nz), qs(0:nz), diss(0:nz)
    double precision:: dq(0:nz), dqs(0:nz)
    double precision:: q0(0:nz), dq0(0:nz)

    ! FFT arrays:
    integer:: zfactors(5)
    double precision:: ztrig(2*nz), rkz(nz)

    ! Others:
    double precision:: eps, t, fac
    integer:: iz, nnu

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

    open(80, file = 'q_profile_hyper.asc', status = 'replace')
    do iz = 0, nz
        write(80,*) z(iz), q(iz), q0(iz), q(iz) - q0(iz)
    enddo
    close(80)

    open(80, file = 'dq_profile_hyper.asc', status = 'replace')
    do iz = 0, nz
        write(80,*) z(iz), dq(iz), dq0(iz), dq(iz) - dq0(iz)
    enddo
    close(80)
end program test_hyper
