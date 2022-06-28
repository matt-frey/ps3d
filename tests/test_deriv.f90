program test_deriv

!==============================================================================
  ! Tests differentiation of phi(z) = 1 - 2z^3 by first decomposing phi
  ! into phi_L = 1 - 2z and phi_S = phi - phi_L, using a sine series for
  ! phi_S, differentiating this spectrally, returning the result as a
  ! cosine series, and adding back phi_L.

!   gfortran -O3 -o test_deriv stafft.f90 deriv1d.f90 test_deriv.f90
!==============================================================================

! Import the 1D FFT module:
use stafft, only : initfft, dst, dct
use deriv1d, only : init_deriv

implicit none

! Resolution in z:
integer,parameter:: nz = 1024

!Generic double precision numerical constants:
double precision,parameter:: zero = 0.d0, one = 1.d0, two = 2.d0
double precision,parameter:: three = 3.d0, four = 4.d0, six = 6.d0

! Array to differentiate, and the linear part:
double precision:: phi(0:nz), phil(0:nz)

! Corresponding z derivative (exact and approximate):
double precision:: dphie(0:nz), dphia(0:nz)
double precision:: dphif(0:nz), dphiq(0:nz)

! FFT arrays:
integer:: zfactors(5)
double precision:: ztrig(2*nz), rkz(nz)

! Grid length:
double precision,parameter:: dz = one/dble(nz)
double precision,parameter:: dzi = one/dz, hdzi = one/(two*dz)

! Other constants:
double precision:: z, dphil
double precision:: earms, eamax
double precision:: efrms, efmax
double precision:: eqrms, eqmax
integer:: iz

!-------------------------------------------------------------
! Initialise Fourier transform and wavenumbers:
call initfft(nz, zfactors, ztrig)
call init_deriv(nz, one, rkz)

!-------------------------------------------------------------
! Set up function, its linear part and its derivative:
do iz = 0, nz
   z = dz*dble(iz)
   phi(iz) = one - two*z**3
   phil(iz) = one - two*z
   dphie(iz) = -six*z**2
enddo
dphil = -two

! Finite differences:
dphif(1:nz-1) = hdzi*(phi(2:nz) - phi(0:nz-2))
dphiq(1:nz-1) = dphif(1:nz-1)
! Linear extrapolation at boundaries:
dphif(0) = dzi*(phi(1) - phi(0))
dphif(nz) = dzi*(phi(nz) - phi(nz-1))
! Quadratic extrapolation at boundaries:
dphiq(0) = hdzi*(four*phi(1) - three*phi(0) - phi(2))
dphiq(nz) = hdzi*(three*phi(nz) + phi(nz-2) - four*phi(nz-1))

!write(*,*) dphie(0), dphie(nz)
!write(*,*) dphif(0), dphif(nz)
!write(*,*) dphiq(0), dphiq(nz)

!-------------------------------------------------------------
! Do a sine DFT on phi - phil:
dphia = phi - phil
call dst(1, nz, dphia(1:nz), ztrig, zfactors)

! Take derivative spectrally:
dphia(0) = zero
dphia(1:nz) = rkz*dphia(1:nz)

! Invert to physical space (now a cosine DFT is needed):
call dct(1, nz, dphia, ztrig, zfactors)

! Add back derivative of linear part:
dphia = dphia + dphil

!-------------------------------------------------------------
! Compute max and rms error:
eamax = maxval(abs(dphia - dphie))
efmax = maxval(abs(dphif - dphie))
eqmax = maxval(abs(dphiq - dphie))
dphia = (dphia - dphie)**2
dphif = (dphif - dphie)**2
dphiq = (dphiq - dphie)**2
earms = sqrt((0.5d0*(dphia(0)+dphia(nz)) + sum(dphia(1:nz-1)))/dble(nz))
efrms = sqrt((0.5d0*(dphif(0)+dphif(nz)) + sum(dphif(1:nz-1)))/dble(nz))
eqrms = sqrt((0.5d0*(dphiq(0)+dphiq(nz)) + sum(dphiq(1:nz-1)))/dble(nz))

write(*,*)
write(*,*) ' nz, max/rms errors for decompose, linear and quadratic:'
write(*,'(i4,6(2x,1p,e14.7))') nz, eamax, earms, efmax, efrms, eqmax, eqrms
write(*,*)

open(88,file='errors',status='replace')
write(88,'(i4,6(2x,1p,e14.7))') nz, eamax, earms, efmax, efrms, eqmax, eqrms
close(88)

end program test_deriv
