module constants

 !Module containing all non-modifiable parameters.

use parameters

 !Sizes of records used in unformatted writes of real*4 data:
integer,parameter:: nhgp=ng*ng,nzm1=nz-1,nzm2=nz-2
integer,parameter:: nhbytes=4*(nhgp+1),ntbytes=4*(nhgp*(nz+1)+1)

 !Generic double precision numerical constants:
double precision,parameter:: zero=0.d0,one=1.d0,two=2.d0
double precision,parameter:: three=3.d0,four=4.d0,six=6.d0
double precision,parameter:: f12=one/two,f13=one/three,f23=two/three
double precision,parameter:: f16=one/six,f56=5.d0/six,f112=one/12.d0
double precision,parameter:: twopi=two*pi

 !Inverse of domain depth:
double precision,parameter:: depthi=one/depth

 !Grid constants:
double precision,parameter:: gl=twopi/dble(ng),dsumi=one/dble(ng*ng)
double precision,parameter:: dz=depth/dble(nz),dz2=dz/two
double precision,parameter:: dzi=one/dz,hdzi=f12*dzi,dzisq=dzi**2

end module constants
