module spectral
    use constants
    use parameters, only : nx, ny, nz
    use sta3dfft
    implicit none

    ! Spectral filter:
    double precision, allocatable :: filt(:, :, :)
    double precision, allocatable :: skx(:), sky(:), skz(:)

!     !Common arrays, constants:
!     double precision:: yh0(0:ny),yh1(0:ny),pbar(0:ny)
!     double precision:: rkx(0:nxm1),hrkx(nx),rky(ny)
!     double precision:: green(0:nxm1,0:ny),hdis(0:nxm1,0:ny)
!     double precision::  filt(0:nxm1,0:ny),decy(nym1,0:nxm1)
!     double precision:: spmf(0:max(nx,ny)),alk(max(nx,ny))
!     double precision:: xtrig(2*nx),ytrig(2*ny)
!     integer:: xfactors(5),yfactors(5)
!     integer:: kmag(0:nxm1,0:ny),kmax

    !====================================================================!
    ! From main code: call init_spectral(bbdif)   to initialise          !
    ! then            call main_invert(zz,uu,vv)  to perform inversion   !
    !====================================================================!

    contains

        !=====================================================================

        subroutine init_spectral(bbdif)
            double precision, intent(in) :: bbdif ! (bbdif = max(b) - min(b) at t = 0):
            double precision             :: kxmaxi, kymaxi, kzmaxi
!             double precision             :: fac,yg,scx,scy,rkxmax,rkymax
!             double precision             :: delk,delki,snorm,div,visc
!             integer                      :: iy, kx, kxc, ky, k


            allocate(filt(0:nz, 0:nx-1, 0:ny-1))
            allocate(skx(nx))
            allocate(sky(ny))
            allocate(skz(0:nz))
!
!
!             !---------------------------------------------------------------------
!             !Fractional y grid values:
!             fac=one/dble(ny)
!             do iy=0,ny
!                 yh1(iy)=fac*dble(iy)
!                 yh0(iy)=one-yh1(iy)
!             enddo
!
!             !Define part of streamfunction proportional to the mean vorticity:
!             do iy=0,ny
!                 yg=gly*dble(iy)
!                 pbar(iy)=f12*yg*(yg-elly)
!             enddo
!
            !---------------------------------------------------------------------
            !Set up FFTs:
            call init3dfft(nx, ny, nz, extent)


            !----------------------------------------------------------
            !Define Hou and Li filter:
            kxmaxi = one / maxval(rkx)
            skx = -36.d0 * (kxmaxi * rkx) ** 36
            kymaxi = one/maxval(rky)
            sky = -36.d0 * (kymaxi * rky) ** 36
            kzmaxi = one / maxval(rkz)
            skz(0) = zero
            skz(1:nz) = -36.d0 * (kzmaxi * rkz) ** 36
            do kz = 0, nz
                do kx = 0, nx-1
                    do ky = 0, ny-1
                        filt(kx, ky, kz) = dexp(skx(kx) + sky(ky) + skz(kz))
                    enddo
                enddo
            enddo

!             !---------------------------------------------------------------------
!             !Define Green function:
!             green(0, 0) = zero
!             do kx = 1, nx-1
!                 green(kx, 0) = -one / rkx(kx) ** 2
!             enddo
!
!             do ky=1,ny
!                 green(0, ky) = -one / rky(ky) ** 2
!             enddo
!
!             do kz = 1, nz
!                 green(kx, ky, kz) =
!
!             do ky=1,ny
!                 do kx=1,nxm1
!                     green(kx,ky)=-one/(rkx(kx)**2+rky(ky)**2)
!                 enddo
!             enddo
!
!             !---------------------------------------------------------------------
!             !Hyperbolic functions used for solutions of Laplace's equation:
!             decy(1:nym1,0)=yh1(1:nym1)
!             do kx=1,nxm1
!                 fac=rkx(kx)*elly
!                 div=one/(one-exp(-two*fac))
!                 decy(1:nym1,kx)=(exp(-fac*(one-yh1(1:nym1)))- &
!                                 exp(-fac*(one+yh1(1:nym1))))*div
!             enddo
!
!             !---------------------------------------------------------------------
!             ! Damping, viscous or hyperviscous:
!             if (nnu .eq. 1) then
!                 !Define viscosity:
!                 visc=prediss*sqrt(bbdif/rkxmax**3)
!                 write(*,'(a,1p,e14.7)') ' Viscosity nu = ',visc
!
!                 !Define spectral dissipation operator:
!                 hdis(0,0)=zero
!                 do kx=1,nxm1
!                     hdis(kx,0)=visc*rkx(kx)**2
!                 enddo
!                 do ky=1,ny
!                     hdis(0,ky)=visc*rky(ky)**2
!                 enddo
!                 do ky=1,ny
!                     do kx=1,nxm1
!                         hdis(kx,ky)=visc*(rkx(kx)**2+rky(ky)**2)
!                     enddo
!                 enddo
!
!             else
!                 !Define hyperviscosity:
!                 visc=prediss/max(rkxmax,rkymax)**(2*nnu)
!                 write(*,'(a,1p,e14.7)') ' Hyperviscosity nu = ',visc
!
!                 !Define dissipation operator:
!                 hdis(0,0)=zero
!                 do kx=1,nxm1
!                     hdis(kx,0)=visc*rkx(kx)**(2*nnu)
!                 enddo
!                 do ky=1,ny
!                     hdis(0,ky)=visc*rky(ky)**(2*nnu)
!                 enddo
!                 do ky=1,ny
!                     do kx=1,nxm1
!                     hdis(kx,ky)=visc*(rkx(kx)**2+rky(ky)**2)**nnu
!                     enddo
!                 enddo
!             endif
!
        end subroutine init_spectral

        !=====================================================================

        ! Given the vorticity (xs, es, zs) in spectral space, this routine computes
        ! the velocity (uu, vv, ww) and (xi, eta, zeta) in physical space.
        subroutine vor2vel(xs, es, zs,  uu, vv, ww)
            double precision, intent(in)  :: xs(0:nx-1, 0:ny-1, 0:nz)
            double precision, intent(in)  :: es(0:nx-1, 0:ny-1, 0:nz)
            double precision, intent(in)  :: zs(0:nx-1, 0:ny-1, 0:nz)
            double precision, intent(out) :: uu(0:nz, 0:ny-1, 0:nx-1)
            double precision, intent(out) :: vv(0:nz, 0:ny-1, 0:nx-1)
            double precision, intent(out) :: ww(0:nz, 0:ny-1, 0:nx-1)


        end subroutine vor2vel

!         ! Given the vorticity zs in spectral space, this routine computes
!         ! the velocity (uu,vv) and zz in physical space.
!         subroutine main_invert(zs,uu,vv,zz)
!
!             !Input array (spectral):
!             double precision:: zs(0:nxm1,0:ny)
!             !Output arrays (physical):
!             double precision:: zz(0:ny,0:nxm1),uu(0:ny,0:nxm1),vv(0:ny,0:nxm1)
!
!             !Local arrays:
!             double precision:: ss(0:nxm1,0:ny),pp(0:ny,0:nxm1)
!             double precision:: pbot(0:nxm1),ptop(0:nxm1),cppy(nym1,0:nxm1)
!
!             !Other quantities:
!             double precision:: zbar
!             integer:: kx,ky,ix,iy
!
!             !--------------------------------
!             !Solve for psi (pp):
!
!             !(1) compute mean zz (zbar):
!             ss=zs
!             call spctop_fc(nx,ny,ss,zz,xfactors,yfactors,xtrig,ytrig)
!             !Now zz contains the vorticity in physical space
!
!             zbar=(f12*sum(zz(0,:)+zz(ny,:))+sum(zz(1:nym1,:)))*dsumi
!
!             !(2) Remove mean vorticity from zz:
!             pp=zz-zbar
!
!             !(3) Invert vorticity to get uncorrected streamfunction pp:
!             ss=green*zs
!             call spctop_fc(nx,ny,ss,pp,xfactors,yfactors,xtrig,ytrig)
!
!             !(4) Add part of pp due to mean vorticity:
!             do ix=0,nxm1
!                 pp(:,ix)=pp(:,ix)+zbar*pbar
!             enddo
!
!             !(5) Do a sine transform of pp at y = ymin and ymax and obtain the
!             !    interior field (cppy) that must be subtracted to give pp = 0
!             !    at y = ymin and ymax:
!             pbot=pp(0,:)
!             ptop=pp(ny,:)
!             call forfft(1,nx,pbot,xtrig,xfactors)
!             call forfft(1,nx,ptop,xtrig,xfactors)
!
!             !Define the interior semi-spectral field:
!             do kx=0,nxm1
!                 do iy=1,nym1
!                     cppy(iy,kx)=pbot(kx)*decy(ny-iy,kx)+ptop(kx)*decy(iy,kx)
!                 enddo
!             enddo
!             !Invert using a full transform in x:
!             call revfft(nym1,nx,cppy,xtrig,xfactors)
!
!             !(6) Remove cppy to obtain the final streamfunction pp:
!             pp(0,:)=zero
!             pp(ny,:)=zero
!             pp(1:nym1,:)=pp(1:nym1,:)-cppy(1:nym1,:)
!
!             !(7) Compute velocity field from pp:
!             call getvel(pp,uu,vv)
!
!         end subroutine main_invert
!
!         !=====================================================================
!
!         subroutine getvel(pp,uu,vv)
!             ! Computes the velocity components uu & vv from the streamfunction
!             ! pp via uu = -d(pp)/dy and vv = d(pp)/dx.
!             ! *** pp, uu & vv are all in physical space
!             ! *** and include the domain edges.
!
!             !Passed arrays:
!             double precision:: pp(0:ny,0:nxm1),uu(0:ny,0:nxm1),vv(0:ny,0:nxm1)
!
!             !Local arrays:
!             double precision:: ppi(ny  ,0:nxm1),pps(0:nxm1  ,ny)
!             double precision:: ppx(0:nxm1,  ny),vvi(  ny,0:nxm1)
!             double precision:: ppy(0:nxm1,0:ny)
!
!             !-------------------------------------------------------------------
!             !Copy non-zero interior values of pp to ppi:
!             ppi(1:nym1,:)=pp(1:nym1,:)
!
!             !Transform ppi to spectral space:
!             call ptospc_fs(nx,ny,ppi,pps,xfactors,yfactors,xtrig,ytrig)
!
!             !Apply de-aliasing filter:
!             pps(:,1:ny)=pps(:,1:ny)*filt(:,1:ny)
!
!             !Compute d(ppi)/dx = ppx spectrally:
!             call xderiv_fs(nx,ny,hrkx,pps,ppx)
!
!             !Transform ppx back to physical space as vvi:
!             call spctop_fs(nx,ny,ppx,vvi,xfactors,yfactors,xtrig,ytrig)
!
!             !Copy vvi into vv and add on zero edge values at iy = 0 & ny:
!             vv(0,:)=zero
!             vv(1:nym1,:)=vvi(1:nym1,:)
!             vv(ny,:)=zero
!
!             !-------------------------------------------------------------------
!             !Compute d(ppi)/dy = ppy spectrally:
!             call yderiv_fs(nx,ny,rky,pps,ppy)
!
!             !Transform ppy back to physical space as uu:
!             call spctop_fc(nx,ny,ppy,uu,xfactors,yfactors,xtrig,ytrig)
!
!             !Correct sign:
!             uu=-uu
!
!         end subroutine getvel
!
!         !=====================================================================
!
!         ! Computes the 1d spectrum of a spectral field var which is
!         ! periodic in x and represented by a cosine series in y.
!         ! Returns the result in spec.
!         subroutine spec1d_fc(var,spec)
!
!             !Passed arrays:
!             double precision:: var(0:nxm1,0:ny),spec(0:max(nx,ny))
!             !Local indices:
!             integer:: k,kx,ky
!
!             !-------------------------------------------------------------------
!             !Initialise spectrum:
!             spec(0:kmax)=zero
!
!             !x and y-independent mode:
!             k=kmag(0,0)
!             spec(k)=spec(k)+f14*var(0,0)**2
!
!             !y-independent mode:
!             do kx=1,nxm1
!                 k=kmag(kx,0)
!                 spec(k)=spec(k)+f12*var(kx,0)**2
!             enddo
!
!             !x-independent mode:
!             do ky=1,ny
!                 k=kmag(0,ky)
!                 spec(k)=spec(k)+f12*var(0,ky)**2
!             enddo
!
!             !All other modes:
!             do ky=1,ny
!                 do kx=1,nxm1
!                     k=kmag(kx,ky)
!                     spec(k)=spec(k)+var(kx,ky)**2
!                 enddo
!             enddo
!         end subroutine spec1d_fc

end module spectral
