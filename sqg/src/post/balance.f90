!#########################################################################
!  Computes the quasi-geostrophic balanced flow for a specified surface
!  distribution of buoyancy b_0(x,y) in b0_init.r8.

!  The results are output to various files in the qgbal subdirectory.

!  Updated 3 February 2024 by D G Dritschel @ St Andrews
!#########################################################################

program balance

 !Import spectral module:
use spectral

implicit none

 !Various quantities needed below (global variables):
double precision:: b0(ng,ng)
double precision:: bb(ng,ng,0:nz),pp(ng,ng,0:nz)
double precision:: ox(ng,ng,0:nz),oy(ng,ng,0:nz),oz(ng,ng,0:nz)
double precision:: ux(ng,ng,0:nz),uy(ng,ng,0:nz)
double precision:: gg(ng,ng,0:nz)
double precision:: t

!---------------------------------------------------------
 !Define fixed arrays and constants and read initial data:
call initialise

!---------------------------------------------------------------
 !Find QG balanced 3D flow:
write(*,*)
write(*,*) ' -------------------------------------------------'
write(*,*) ' Finding fields using Quasi-Geostrophic balance...'
call qgbal(b0,bb,ox,oy,oz,pp,ux,uy,gg)

 !Write data to qgbal subdirectory:
call savedata('qgbal')

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Internal subroutine definitions (inherit global variables):
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

contains

!=============================================================

subroutine initialise
! Initialises spectral module and reads data

implicit none

!---------------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

!---------------------------------------------------------------
 !Read dimensionless surface bouyancy, b0:
open(11,file='b0_init.r8',form='unformatted', &
      access='direct',status='old',recl=2*nhbytes)
read(11,rec=1) t,b0
close(11)

return
end subroutine initialise

!=============================================================

subroutine savedata(directory)
! Writes balanced data to directory (passed)

character(len=5):: directory

double precision:: sfac

 !Undo scaling of variables for direct use in PS3D or EPIC:
sfac=cof*bvf
bb=sfac*bb
ox=bvf*ox
oy=bvf*oy
oz=cof*oz
sfac=cof**2
pp=sfac*pp
ux=cof*ux
uy=cof*uy

!------------------------------------------------------------------
 !Write data:
write(*,*)
write(*,*) ' Writing x vorticity component to '//directory//'/ox.r4'
open(11,file=directory//'/ox.r4',form='unformatted',access='direct', &
     status='replace',recl=ntbytes)
write(11,rec=1) real(zero),real(ox)
close(11)
write(*,'(2(a,f10.6))') ' Min value = ',minval(ox), &
                      '   Max value = ',maxval(ox)

write(*,*)
write(*,*) ' Writing y vorticity component to '//directory//'/oy.r4'
open(11,file=directory//'/oy.r4',form='unformatted',access='direct', &
     status='replace',recl=ntbytes)
write(11,rec=1) real(zero),real(oy)
close(11)
write(*,'(2(a,f10.6))') ' Min value = ',minval(oy), &
                      '   Max value = ',maxval(oy)

write(*,*)
write(*,*) ' Writing z vorticity component to '//directory//'/oz.r4'
open(11,file=directory//'/oz.r4',form='unformatted',access='direct', &
     status='replace',recl=ntbytes)
write(11,rec=1) real(zero),real(oz)
close(11)
write(*,'(2(a,f10.6))') ' Min value = ',minval(oz), &
                      '   Max value = ',maxval(oz)

write(*,*)
write(*,*) ' Writing buoyancy anomaly to '//directory//'/ba.r4'
open(11,file=directory//'/ba.r4',form='unformatted',access='direct', &
     status='replace',recl=ntbytes)
write(11,rec=1) real(zero),real(bb)
close(11)
write(*,'(2(a,f10.6))') ' Min value = ',minval(bb), &
                      '   Max value = ',maxval(bb)

write(*,*)
write(*,*) ' Writing perturbation pressure to '//directory//'/pp.r4'
open(11,file=directory//'/pp.r4',form='unformatted',access='direct', &
     status='replace',recl=ntbytes)
write(11,rec=1) real(zero),real(pp)
close(11)
write(*,'(2(a,f10.6))') ' Min value = ',minval(pp), &
                      '   Max value = ',maxval(pp)

write(*,*)
write(*,*) ' Writing x velocity component to '//directory//'/ux.r4'
open(11,file=directory//'/ux.r4',form='unformatted',access='direct', &
     status='replace',recl=ntbytes)
write(11,rec=1) real(zero),real(ux)
close(11)
write(*,'(2(a,f10.6))') ' Min value = ',minval(ux), &
                      '   Max value = ',maxval(ux)

write(*,*)
write(*,*) ' Writing y velocity component to '//directory//'/uy.r4'
open(11,file=directory//'/uy.r4',form='unformatted',access='direct', &
     status='replace',recl=ntbytes)
write(11,rec=1) real(zero),real(uy)
close(11)
write(*,'(2(a,f10.6))') ' Min value = ',minval(uy), &
                      '   Max value = ',maxval(uy)

write(*,*)
write(*,*) ' Writing static stability Gamma to '//directory//'/gg.r4'
open(11,file=directory//'/gg.r4',form='unformatted',access='direct', &
     status='replace',recl=ntbytes)
write(11,rec=1) real(zero),real(gg)
close(11)
write(*,'(2(a,f10.6))') ' Min value = ',minval(gg), &
                      '   Max value = ',maxval(gg)

 !Compute inverse Richardson number:
sfac=one/bvf**2
gg=sfac*(ox**2+oy**2)/gg

write(*,*)
write(*,*) ' Writing inverse Richardson number to '//directory//'/ri.r4'
open(11,file=directory//'/ri.r4',form='unformatted',access='direct', &
     status='replace',recl=ntbytes)
write(11,rec=1) real(zero),real(gg)
close(11)
write(*,'(2(a,f10.6))') ' Min value = ',minval(gg), &
                      '   Max value = ',maxval(gg)

end subroutine savedata

 !End main program
end program balance
!=======================================================================
