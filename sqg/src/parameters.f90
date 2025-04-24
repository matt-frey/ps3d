module parameters

! Module containing all the modifiable parameters (except pi below).

! The number pi (*** non-modifiable ***):
double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0

! ==> Numerical parameters <==
integer,parameter:: ng=512,nz=64
! ng     : horizontal grid resolution in both x and y
!          (Note: the domain is a 2*pi periodic box horizontally)
! nz     : vertical grid resolution (nz+1 grid points including edges)

! ==> Physical parameters <==
double precision,parameter:: cof=1.d0,bvf=8.d0
double precision,parameter:: depth=pi/4.d0
! cof    : Coriolis frequency f
! bvf    : Buoyancy frequency N of the background linear stratification
! depth  : NH/f, total fluid depth H stretched by N/f; note the domain
!          width is L = 2*pi. To have an isotropic grid in unscaled
!          coordinates, ensure L/ng = H/nz, i.e.
!              2*pi/ng = cof*depth/(bvf*nz)
!----------------------------------------------------------------

end module parameters
