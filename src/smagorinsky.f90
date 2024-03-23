module smagorinsky_mod
    use constants, only : f12, f13, two
    use parameters, only : vcell
    use dimensions, only : I_X, I_Y, I_Z
    use mpi_layout, only : box
    use fields, only : svorts, svel, svor, vor
#ifdef ENABLE_BUOYANCY
    use fields, only : sbuoy, sbuoys
#endif
    use sta3dfft, only : diffx      &
                       , diffy      &
                       , fftxys2p   &
                       , fftxyp2s
    use inversion_utils, only : field_decompose_physical    &
                              , field_combine_physical      &
                              , central_diffz
    implicit none

    private

    ! Lilly's Smagorinsky coefficient for homogeneous isotropic turbulence (HIT):
    double precision, parameter :: c_s = 0.25d0 !0.173d0

    ! velocity strain indices
    integer, parameter :: I_DUDX = 1 & ! index for du/dx strain component
                        , I_DUDY = 2 & ! index for du/dy strain component
                        , I_DVDY = 3 & ! index for dv/dy strain component
                        , I_DWDX = 4 & ! index for dw/dx strain component
                        , I_DWDY = 5   ! index for dw/dy strain component

    public :: apply_smagorinsky
#ifdef ENABLE_BUOYANCY
    public :: apply_smagorinsky_buoyancy
#endif

    contains
        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Get Smagorinsky viscosity *smag* in physical space
        subroutine smagorinsky(smag)
            double precision, intent(out) :: smag(box%lo(3):box%hi(3), &
                                                  box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1))
            double precision              :: velgrad(box%lo(3):box%hi(3), & ! 0:nz
                                                     box%lo(2):box%hi(2), &
                                                     box%lo(1):box%hi(1), &
                                                     5)
            double precision              :: lscale


            call vel2velgrad(velgrad)

            !------------------------------------------------------------------
            ! Obtain Smagorinsky viscosity:
            ! Calculate straing magnitude and store in *smag*
            call strain_magnitude(velgrad, smag)

            ! Multiply with length scale (Smagorinsky coefficient: c_s = 0.173 according to Lilly)
            ! Deardorff: (dx * dy * dz) ** f13
            lscale = (c_s * vcell ** f13) ** 2

            !$omp parallel workshare
            smag = lscale * smag
            !$omp end parallel workshare

        end subroutine smagorinsky

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Compute the gridded velocity gradient tensor
        subroutine vel2velgrad(velgrad)
            double precision, intent(out) :: velgrad(box%lo(3):box%hi(3),  & ! 0:nz
                                                     box%lo(2):box%hi(2),  &
                                                     box%lo(1):box%hi(1),  &
                                                     5)
            double precision              :: ds(box%lo(3):box%hi(3), & ! 0:nz
                                                box%lo(2):box%hi(2), &
                                                box%lo(1):box%hi(1)) ! semi-spectral derivatives

            ! x component:
            call diffx(svel(:, :, :, I_X), ds)           ! u_x = du/dx in semi-spectral space
            call fftxys2p(ds, velgrad(:, :, :, I_DUDX)) ! u_x in physical space

            call diffy(svel(:, :, :, I_X), ds)           ! u_y = du/dy in semi-spectral space
            call fftxys2p(ds, velgrad(:, :, :, I_DUDY)) ! u_y in physical space

            call diffx(svel(:, :, :, I_Z), ds)           ! w_x = dw/dx in semi-spectral space
            call fftxys2p(ds, velgrad(:, :, :, I_DWDX)) ! w_x in physical space

            ! y & z components:
            call diffy(svel(:, :, :, I_Y), ds)           ! v_y = dv/dy in semi-spectral space
            call fftxys2p(ds, velgrad(:, :, :, I_DVDY)) ! v_y in physical space

            call diffy(svel(:, :, :, I_Z), ds)           ! w_y = dw/dy in semi-spectral space
            call fftxys2p(ds, velgrad(:, :, :, I_DWDY)) ! w_y in physical space

        end subroutine vel2velgrad

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Calculate the velocity strain magnitude (and store in smag) for each grid point
        ! @param[in] velocity gradient tensor at grid point
        ! @param[in] vorticity at grid point
        ! @returns 3x3 strain matrix
        subroutine strain_magnitude(velgrad, smag)
            double precision, intent(in) :: velgrad(box%lo(3):box%hi(3),  & ! 0:nz
                                                    box%lo(2):box%hi(2),  &
                                                    box%lo(1):box%hi(1),  &
                                                     5)
            double precision, intent(out) :: smag(box%lo(3):box%hi(3), &
                                                  box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1))
            double precision              :: strain(3,3)
            integer                       :: ix, iy, iz

            do ix = box%lo(1), box%hi(1)
                do iy = box%lo(2), box%hi(2)
                    do iz = box%lo(3), box%hi(3) ! = 0, nz
                        ! get local symmetrised strain matrix, i.e. 1/ 2 * (S + S^T)
                        ! where
                        !     /u_x u_y u_z\
                        ! S = |v_x v_y v_z|
                        !     \w_x w_y w_z/
                        ! with u_* = du/d* (also derivatives of v and w).
                        ! The derivatives dv/dx, du/dz, dv/dz and dw/dz are calculated
                        ! with vorticity or the assumption of incompressibility
                        ! (du/dx + dv/dy + dw/dz = 0):
                        !    dv/dx = \zeta + du/dy
                        !    du/dz = \eta + dw/dx
                        !    dv/dz = dw/dy - \xi
                        !    dw/dz = - (du/dx + dv/dy)
                        !
                        !                         /  2 * u_x  u_y + v_x u_z + w_x\
                        ! 1/2 * (S + S^T) = 1/2 * |u_y + v_x   2 * v_y  v_z + w_y|
                        !                         \u_z + w_x  v_z + w_y   2 * w_z/
                        !
                        ! S11 = du/dx
                        ! S12 = 1/2 * (du/dy + dv/dx) = 1/2 * (2 * du/dy + \zeta) = du/dy + 1/2 * \zeta
                        ! S13 = 1/2 * (du/dz + dw/dx) = 1/2 * (\eta + 2 * dw/dx) = 1/2 * \eta + dw/dx
                        ! S22 = dv/dy
                        ! S23 = 1/2 * (dv/dz + dw/dy) = 1/2 * (2 * dw/dy - \xi) = dw/dy - 1/2 * \xi
                        ! S33 = dw/dz = - (du/dx + dv/dy)
                        strain(1, 1) = velgrad(iz, iy, ix, I_DUDX)                                 ! S11
                        strain(1, 2) = velgrad(iz, iy, ix, I_DUDY) + f12 * vor(iz, iy, ix, I_Z)    ! S12
                        strain(1, 3) = velgrad(iz, iy, ix, I_DWDX) + f12 * vor(iz, iy, ix, I_Y)    ! S13
                        strain(2, 1) = strain(1, 2)
                        strain(2, 2) = velgrad(iz, iy, ix, I_DVDY)                                ! S22
                        strain(2, 3) = velgrad(iz, iy, ix, I_DWDY) - f12 * vor(iz, iy, ix, I_X)   ! S23
                        strain(3, 1) = strain(1, 3)
                        strain(3, 2) = strain(2, 3)
                        strain(3, 3) = -(velgrad(iz, iy, ix, I_DUDX) + velgrad(iz, iy, ix, I_DVDY)) ! S33


                        smag(iz, iy, ix) = dsqrt(two * (strain(1, 1) * strain(1, 1) + &
                                                        strain(1, 2) * strain(1, 2) + &
                                                        strain(1, 3) * strain(1, 3) + &
                                                        strain(2, 1) * strain(2, 1) + &
                                                        strain(2, 2) * strain(2, 2) + &
                                                        strain(2, 3) * strain(2, 3) + &
                                                        strain(3, 1) * strain(3, 1) + &
                                                        strain(3, 2) * strain(3, 2) + &
                                                        strain(3, 3) * strain(3, 3) ))
                   enddo
                enddo
            enddo

        end subroutine strain_magnitude

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Calculate grad(smag * grad(vor)) where "smag" is the (Smagorinsky) eddy viscosity
        subroutine apply_smagorinsky
            double precision :: ds(box%lo(3):box%hi(3),  & ! 0:nz
                                   box%lo(2):box%hi(2),  &
                                   box%lo(1):box%hi(1))
            double precision :: dp(box%lo(3):box%hi(3),  & ! 0:nz
                                   box%lo(2):box%hi(2),  &
                                   box%lo(1):box%hi(1))
            double precision :: smag(box%lo(3):box%hi(3), &
                                     box%lo(2):box%hi(2), &
                                     box%lo(1):box%hi(1))
            double precision :: omg(box%lo(3):box%hi(3),  & ! 0:nz
                                    box%lo(2):box%hi(2),  &
                                    box%lo(1):box%hi(1))
            integer          :: nc


            call smagorinsky(smag)

            ! grad(smag * grad(vor)):
            ! d(smag * dxi/dx) / dx   + d(smag * dxi/dy) / dy   + d(smag * dxi/dz) / dz
            ! d(smag * deta/dx) / dx  + d(smag * deta/dy) / dy  + d(smag * deta/dz) / dz
            ! d(smag * dzeta/dx) / dx + d(smag * dzeta/dy) / dy + d(smag * dzeta/dz) / dz
            do nc = 1, 3

                !--------------------------------------------------------------
                ! x-derivatives:

                ! dxi/dx, deta/dx or dzeta/dx
                call diffx(svor(:, :, :, nc), ds)

                call field_combine_physical(ds, dp)

                ! apply eddy viscosity:
                dp = smag * dp

                call field_decompose_physical(dp, ds)

                ! d^2*/dx^2
                call diffx(ds, omg)

                svorts(:, :, :, nc) = svorts(:, :, :, nc) + omg

                !--------------------------------------------------------------
                ! y-derivatives:

                ! dxi/dy, deta/dy or dzeta/dy
                call diffy(svor(:, :, :, nc), ds)

                call field_combine_physical(ds, dp)

                ! apply eddy viscosity:
                dp = smag * dp

                call field_decompose_physical(dp, ds)

                ! d^2*/dy^2
                call diffy(ds, omg)

                svorts(:, :, :, nc) = svorts(:, :, :, nc) + omg

                !--------------------------------------------------------------
                ! z-derivatives:

                ! dxi/dz, deta/dz or dzeta/dz
                call field_combine_physical(svor(:, :, :, nc), omg)

                call central_diffz(omg, dp)

                ! apply eddy viscosity:
                omg = smag * dp

                ! d^2*/dy^2
                call central_diffz(omg, dp)

                call field_decompose_physical(dp, omg)

                svorts(:, :, :, nc) = svorts(:, :, :, nc) + omg

            enddo

        end subroutine apply_smagorinsky

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#ifdef ENABLE_BUOYANCY
        ! Calculate grad(smag * grad(buoy)) where "smag" is the (Smagorinsky) eddy viscosity
        subroutine apply_smagorinsky_buoyancy
            double precision :: ds(box%lo(3):box%hi(3),  & ! 0:nz
                                   box%lo(2):box%hi(2),  &
                                   box%lo(1):box%hi(1))
            double precision :: dp(box%lo(3):box%hi(3),  & ! 0:nz
                                   box%lo(2):box%hi(2),  &
                                   box%lo(1):box%hi(1))
            double precision :: smag(box%lo(3):box%hi(3), &
                                     box%lo(2):box%hi(2), &
                                     box%lo(1):box%hi(1))
            double precision :: omg(box%lo(3):box%hi(3),  & ! 0:nz
                                    box%lo(2):box%hi(2),  &
                                    box%lo(1):box%hi(1))

            call smagorinsky(smag)

            ! grad(smag * grad(buoy)):
            ! d(smag * dbuoy/dx) / dx   + d(smag * dbuoy/dy) / dy   + d(smag * dbuoy/dz) / dz
            !--------------------------------------------------------------
            ! x-derivatives:
            call diffx(sbuoy, ds)

            call field_combine_physical(ds, dp)

            ! apply eddy viscosity:
            dp = smag * dp

            call field_decompose_physical(dp, ds)

            ! d^2*/dx^2
            call diffx(ds, omg)

            sbuoys(:, :, :) = sbuoys(:, :, :) + omg

            !--------------------------------------------------------------
            ! y-derivatives:

            call diffy(sbuoy, ds)

            call field_combine_physical(ds, dp)

            ! apply eddy viscosity:
            dp = smag * dp

            call field_decompose_physical(dp, ds)

            ! d^2*/dy^2
            call diffy(ds, omg)

            sbuoys = sbuoys + omg

            !--------------------------------------------------------------
            ! z-derivatives:

            call field_combine_physical(sbuoy, omg)

            call central_diffz(omg, dp)

            ! apply eddy viscosity:
            omg = smag * dp

            ! d^2*/dy^2
            call central_diffz(omg, dp)

            call field_decompose_physical(dp, omg)

            sbuoys = sbuoys + omg

       end subroutine apply_smagorinsky_buoyancy
#endif

end module smagorinsky_mod
