program turb
!-----------------------------------------------------------------
!    Generates a random phased buoyancy distribution with a 
!    spectrum Q(k) = c k^{2p-1} * exp[-(p-1)*(k/k_0)^2], p > 1.
!-----------------------------------------------------------------

use spectral

implicit double precision(a-h,o-z)
double precision:: qs(ng,ng),qa(ng,ng)
integer, dimension(:), allocatable :: seed
integer,parameter:: nw=ng/2
integer:: k

! Initialise inversion constants and arrays:
call init_spectral

! Set exponent p in power spectrum:
pow=three

write(*,*) ' We assume initial potential enstrophy spectrum of the form'
write(*,*)
write(*,*) '   Q(k) = c*k^{2p-1}*exp[-(p-1)*(k/k_0)^2]'
write(*,*)
write(*,*) ' We take p = 3.  Enter k_0:'
read(*,*) ak0

write(*,*) ' Enter the maximum |b|/(f*N):'
read(*,*) qeddy

write(*,*) ' Enter an integer seed for the random number generator:'
read(*,*) ngen
call random_seed(size=k)
allocate(seed(1:k))
seed(:)=ngen
do i=1,ngen
  call random_seed(put=seed)
enddo

! Generate potential enstrophy spectrum / k (actually, its square root):
ak0sqi=one/ak0**2
p1=pow-one
do ky=1,nw+1
  do kx=1,nw+1
    aksq=rk(kx)**2+rk(ky)**2
    s=ak0sqi*aksq
    qs(kx,ky)=sqrt(s**p1*exp(-p1*s))
  enddo
enddo

! Apply to generate full spectrum:
do ky=2,nw
  kyc=ng+2-ky
  do kx=2,nw
    kxc=ng+2-kx
    call random_number(uni)
    phix=twopi*uni-pi
    call random_number(uni)
    phiy=twopi*uni-pi
    cx=cos(phix)
    sx=sin(phix)
    cy=cos(phiy)
    sy=sin(phiy)
    amp=qs(kx,ky)
    qs(kx ,ky )=amp*cx*cy
    qs(kxc,ky )=amp*sx*cy
    qs(kx, kyc)=amp*cx*sy
    qs(kxc,kyc)=amp*sx*sy
  enddo
enddo

ky=1
do kx=2,nw
  kxc=ng+2-kx
  call random_number(uni)
  phix=twopi*uni-pi
  cx=cos(phix)
  sx=sin(phix)
  amp=qs(kx,ky)
  qs(kx ,ky )=amp*cx
  qs(kxc,ky )=amp*sx
enddo

kx=1
do ky=2,nw
  kyc=ng+2-ky
  call random_number(uni)
  phiy=twopi*uni-pi
  cy=cos(phiy)
  sy=sin(phiy)
  amp=qs(kx,ky)
  qs(kx ,ky )=amp*cy
  qs(kx, kyc)=amp*sy
enddo

ky=nw+1
do kx=2,nw
  kxc=ng+2-kx
  call random_number(uni)
  phix=twopi*uni-pi
  cx=cos(phix)
  sx=sin(phix)
  amp=qs(kx,ky)
  qs(kx ,ky )=amp*cx
  qs(kxc,ky )=amp*sx
enddo

kx=nw+1
do ky=2,nw
  kyc=ng+2-ky
  call random_number(uni)
  phiy=twopi*uni-pi
  cy=cos(phiy)
  sy=sin(phiy)
  amp=qs(kx,ky)
  qs(kx ,ky )=amp*cy
  qs(kx, kyc)=amp*sy
enddo

qs(1,1)=zero
qs(nw+1,nw+1)=zero

! Transform to physical space:
call spctop(ng,ng,qs,qa,xfactors,yfactors,xtrig,ytrig)

! Work out max/min values and total pot. enstrophy:
ens=zero
qamin=qa(1,1)
qamax=qamin
do ix=1,ng
  do iy=1,ng
    ens=ens+qa(iy,ix)**2
    qamin=min(qamin,qa(iy,ix))
    qamax=max(qamax,qa(iy,ix))
  enddo
enddo
ens=ens/(two*dble(ng*ng))

! Renormalise PV:
fmult=qeddy/max(abs(qamax),abs(qamin))
do ix=1,ng
  do iy=1,ng
    qa(iy,ix)=fmult*qa(iy,ix)
  enddo
enddo

! Work out max/min values and total pot. enstrophy:
ens=ens*fmult
qamin=qamin*fmult
qamax=qamax*fmult

 !Save average for plotting purposes:
open(44,file='average_qq.asc',status='replace')
write(44,*) zero
close(44)

! Write data:
open(11,file='b0_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nhbytes)
write(11,rec=1) zero,qa
close(11)

write(*,*)
write(*,'(a,f12.5)') ' rms b = ',sqrt(2.d0*ens)
write(*,'(a,f12.7,a,f11.7)') ' min b = ',qamin,'  &  max b = ',qamax

end program turb
