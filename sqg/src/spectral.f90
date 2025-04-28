module spectral

! Module containing subroutines for spectral operations, inversion, etc.

use constants
use sta2dfft

 !Common arrays, constants:

 !Squared wavenumber array:
double precision:: rksq(ng,ng)

 !Tridiagonal arrays for inverting the scaled QG operator with Neumann BCs:
double precision:: etdv(ng,ng,0:nz),htdv(ng,ng,0:nz)
double precision:: ap(ng,ng),apb(ng,ng)
double precision:: zlin(0:nz)

 !Tridiagonal arrays for z differentiation (Neumann BCs):
double precision:: etd1(1:nzm1),htd1(1:nzm1)

 !For 2D FFTs:
double precision:: hrkx(ng),hrky(ng),rk(ng)
double precision:: xtrig(2*ng),ytrig(2*ng)
integer:: xfactors(5),yfactors(5)

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Internal subroutine definitions (inherit global variables):
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

contains

!=============================================================

subroutine init_spectral
! Initialises this module

implicit none

!Local variables:
double precision:: a0(ng,ng),a0b(ng,ng)
integer:: kx,ky,k,iz

!----------------------------------------------------------------------
 !Set up 2D FFTs:
call init2dfft(ng,ng,twopi,twopi,xfactors,yfactors,xtrig,ytrig,hrkx,hrky)

 !Define wavenumbers:
rk(1)=zero
do k=1,ng/2-1
  rk(k+1)   =hrkx(2*k)
  rk(ng+1-k)=hrkx(2*k)
enddo
rk(ng/2+1)=hrkx(ng)

!-----------------------------------------------------------------------
 !Define squared wavenumber array:
do ky=1,ng
   do kx=1,ng
      rksq(kx,ky)=rk(kx)**2+rk(ky)**2
   enddo
enddo

!-----------------------------------------------------------------------
 !Tridiagonal coefficients depending only on kx and ky:
a0=-two*dzisq-f56*rksq
a0b=-dzisq-f13*rksq
ap=dzisq-f112*rksq
apb=dzisq-f16*rksq

 !Tridiagonal arrays for inverting the QG stretched Laplace operator
 !with Neumann boundary conditions (4th-order compact scheme):
htdv(:,:,0)=one/a0b
etdv(:,:,0)=-apb*htdv(:,:,0)
do iz=1,nzm1
   htdv(:,:,iz)=one/(a0+ap*etdv(:,:,iz-1))
   etdv(:,:,iz)=-ap*htdv(:,:,iz)
enddo
htdv(:,:,nz)=one/(a0b+apb*etdv(:,:,nzm1))

 !Horizontal mean part must be solved for separately:
htdv(1,1,:)=zero
etdv(1,1,:)=zero

 !Tridiagonal arrays for computing phi_z with given boundary values
 !(4th-order compact scheme):
htd1(1)=one/f23
etd1(1)=-f16*htd1(1)
do iz=2,nzm1
   htd1(iz)=one/(f23+f16*etd1(iz-1))
   etd1(iz)=-f16*htd1(iz)
enddo

 !Define zlin = z + D (vanishes at bottom, iz = 0):
do iz=0,nz
   zlin(iz)=dz*dble(iz)
enddo

return 
end subroutine init_spectral

!=================================================================

subroutine qgbal(b0,bb,xxi,eta,zz,pp,uu,vv,gg)
! Assuming QG balance (hydrostatic and geostrophic), and given
! only the surface distribution of buoyancy (b0), this routine
! finds the scaled buoyancy anomaly b'/(f*N), and the scaled vorticity
! components xi/N, eta/N and zeta/f. Returns these quantities in bb, xxi,
! eta and zz, respectively. Also returns the perturbation pressure / f^2
! in pp, as well as u/f & v/f in uu & vv, and the static stability (Gamma)
! in gg.
  
implicit none

 !Passed arrays:
double precision:: b0(ng,ng)
double precision:: bb(ng,ng,0:nz),xxi(ng,ng,0:nz),eta(ng,ng,0:nz)
double precision:: zz(ng,ng,0:nz),pp(ng,ng,0:nz)
double precision:: uu(ng,ng,0:nz),vv(ng,ng,0:nz)
double precision:: gg(ng,ng,0:nz)

 !Local variables:
double precision:: pz(ng,ng,0:nz)
double precision:: pz0(ng,ng),wkb(ng,ng),wkc(ng,ng)
double precision:: b0bar,q0,sfac
integer:: iz

!-------------------------------------------------------------------
! Scale b0 by f*N:
b0=b0/(cof*bvf)
! Horizontal mean value of b_0(x,y):
b0bar=dsumi*sum(b0)
! Initial guess for the PV anomaly, q_0:
q0=b0bar*depthi
! depthi = 1/D where D is the scaled depth (NH/f).

write(*,*)
write(*,*) ' QG solution for PV anomaly:'
write(*,'(2(a,1p,e14.7))') ' q_0 = ',q0

! Convert b_0 to spectral space as pz0:
wkc=b0
call ptospc(ng,ng,wkc,pz0,xfactors,yfactors,xtrig,ytrig)
! Factor used to convert a physical mean value to a spectral mean value:
sfac=pz0(1,1)/b0bar
! *** Do not re-use pz0 below!

!-------------------------------------------------------------------
! QG solution: source is zero only at the top boundary:
pp(:,:,nz)=-dzi*pz0*htdv(:,:,nz)
! This is a partial tridiagonal solve (dzi = 1/dz)
do iz=nzm1,0,-1
   pp(:,:,iz)=etdv(:,:,iz)*pp(:,:,iz+1)
enddo

! Compute phi_z -> pz:
call zderiv(pp,pz,pz0)
! Here, phi_z = 0 at iz = 0, and phi_z = pz0 (defined above) at iz = nz.

! Overwrite horizontal-mean part with q_0*(z+D):
do iz=0,nz
   pz(1,1,iz)=sfac*q0*zlin(iz)
enddo

! Save physical copy of pz in bb, the 3d buoyancy anomaly field:
call spctop3d(pz,bb,0,nz)

! Compute xi/N = -phi_zx, eta/N = -phi_zy and zeta/f = Lap_h(phi)
! and store in xxi, eta and zz:
do iz=0,nz
   call xderiv(ng,ng,hrkx,pz(:,:,iz),wkb)
   call spctop(ng,ng,wkb,wkc,xfactors,yfactors,xtrig,ytrig)
   xxi(:,:,iz)=-wkc
   call yderiv(ng,ng,hrky,pz(:,:,iz),wkb)
   call spctop(ng,ng,wkb,wkc,xfactors,yfactors,xtrig,ytrig)
   eta(:,:,iz)=-wkc
   wkb=-rksq*pp(:,:,iz)
   call spctop(ng,ng,wkb,wkc,xfactors,yfactors,xtrig,ytrig)
   zz(:,:,iz)=wkc
enddo

! Compute static stability, Gamma (QG approximation of it):
gg=one-zz

! Add horizontal mean part of phi by vertically integrating \bar\phi_z:
pp(1,1,0)=zero
do iz=1,nz
   pp(1,1,iz)=pp(1,1,iz-1)+dz2*(pz(1,1,iz)+pz(1,1,iz-1))
enddo

! Compute u/f, v/f & p/f^2 and store in uu, vv & pp:
do iz=0,nz
   call yderiv(ng,ng,hrky,pp(:,:,iz),wkb)
   call spctop(ng,ng,wkb,wkc,xfactors,yfactors,xtrig,ytrig)
   uu(:,:,iz)=-wkc

   call xderiv(ng,ng,hrky,pp(:,:,iz),wkb)
   call spctop(ng,ng,wkb,wkc,xfactors,yfactors,xtrig,ytrig)
   vv(:,:,iz)=wkc

   wkb=pp(:,:,iz)
   call spctop(ng,ng,wkb,wkc,xfactors,yfactors,xtrig,ytrig)
   pp(:,:,iz)=wkc
enddo

! Undo scaling of b0:
b0=b0*cof*bvf

return 
end subroutine qgbal

!=================================================================

subroutine zderiv(pp,pz,pztop)
! Computes the z derivative of a function pp given that its
! derivative is zero when iz = 0 and equal to pztop when iz = nz.
! Uses a 4th-order compact difference method.

implicit none

 !Passed arrays:
double precision:: pp(ng,ng,0:nz),pz(ng,ng,0:nz),pztop(ng,ng)

 !Local variable:
integer:: iz

!---------------------------------------------------------
pz(:,:,0)=zero
do iz=1,nzm1
   pz(:,:,iz)=hdzi*(pp(:,:,iz+1)-pp(:,:,iz-1))
enddo
pz(:,:,nz)=pztop

! Get interior derivatives by 4th-order compact difference solution:
pz(:,:,nzm1)=pz(:,:,nzm1)-f16*pztop
pz(:,:,1)=pz(:,:,1)*htd1(1)
do iz=2,nzm1
   pz(:,:,iz)=(pz(:,:,iz)-f16*pz(:,:,iz-1))*htd1(iz)
enddo
do iz=nzm2,1,-1
   pz(:,:,iz)=etd1(iz)*pz(:,:,iz+1)+pz(:,:,iz)
enddo

return
end subroutine zderiv

!=================================================================

subroutine ptospc3d(fp,fs,izbeg,izend)
! Transforms a physical 3d field fp to spectral space (horizontally)
! as the array fs.

implicit none

 !Passed variables:
double precision:: fp(ng,ng,0:nz)  !Physical
double precision:: fs(ng,ng,0:nz)  !Spectral
integer:: izbeg,izend

 !Work arrays:
double precision:: wkp(ng,ng)  !Physical
double precision:: wks(ng,ng)  !Spectral
integer:: iz

!---------------------------------------------------------
do iz=izbeg,izend
   wkp=fp(:,:,iz)
   call ptospc(ng,ng,wkp,wks,xfactors,yfactors,xtrig,ytrig)
   fs(:,:,iz)=wks
enddo

return
end subroutine ptospc3d

!=================================================================

subroutine spctop3d(fs,fp,izbeg,izend)
! Transforms a spectral 3d field fs to physical space (horizontally)
! as the array fp.

implicit none

 !Passed variables:
double precision:: fp(ng,ng,0:nz)  !Physical
double precision:: fs(ng,ng,0:nz)  !Spectral
integer:: izbeg,izend

 !Work arrays:
double precision:: wkp(ng,ng)  !Physical
double precision:: wks(ng,ng)  !Spectral
integer:: iz

!---------------------------------------------------------
do iz=izbeg,izend
   wks=fs(:,:,iz)
   call spctop(ng,ng,wks,wkp,xfactors,yfactors,xtrig,ytrig)
   fp(:,:,iz)=wkp
enddo

return
end subroutine spctop3d

!=================================================================

end module spectral
