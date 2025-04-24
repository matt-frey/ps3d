program vstrip
! Initialises a scaled surface buoyancy field b_0(x,y)/(f*N) where f is
! the Coriolis frequency and N is the buoyancy frequency of the background
! hydrostatic state.

! Here, we take b_0 = b_max * e^[-(x - c*sin(y))^2/a^2]

use constants

implicit none

double precision:: b0(ng,ng)
double precision:: bm,x0,c,xfac,x,y
integer:: ix,iy

write(*,*)
write(*,*) ' We take b_0 = b_m * e^{-s^2} where s = (x - c*sin(y))/x_0.'
write(*,*) ' Enter b_m/(f*N), x_0 and c:'
read(*,*) bm,x0,c

xfac=one/x0

do ix=1,ng
  x=gl*dble(ix-1)-pi
  do iy=1,ng
    y=gl*dble(iy-1)-pi
    b0(iy,ix)=bm*exp(-((x-c*sin(y))*xfac)**2)
  enddo
enddo

! Write data:
open(11,file='b0_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nhbytes)
write(11,rec=1) zero,b0
close(11)

end program vstrip
