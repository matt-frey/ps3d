module field_zeta
    use constants, only : zero
    use parameters, only : nz, lower, dx, extent, upper
    use mpi_layout, only : box
    use fields, only : svor, vor, szeta, zeta
    use inversion_utils, only : psim, psip, kh, phip            &
                              , field_decompose_semi_spectral   &
                              , field_combine_physical
    use stafft, only : dct
    use sta3dfft, only : fftxyp2s   &
                       , fftxys2p   &
                       , xfactors   &
                       , yfactors   &
                       , zfactors   &
                       , xtrig      &
                       , ytrig      &
                       , ztrig      &
                       , rkzi       &
                       , diffx      &
                       , diffy      &
                       , ifft2d
    implicit none

    contains

        ! Obtain complete zeta in physical (vor(:, :, :, 3))
        ! and semi-spectral space (svor(:, :, :, 3)) by integrating from zeta_min (iz = 0)
        subroutine combine_zeta
            integer          :: iz, kx, ky
            double precision :: ds(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))    ! mixed-spectral space
            double precision :: fp(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))    ! physical space
            double precision :: as(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))    ! mixed-spectral space
            double precision :: psi(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision :: psi_z(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision :: psi_x(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision :: psi_y(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision :: y, kl, z

            ! d\xi/dx in mixed spectral space
            call diffx(svor(:, :, :, 1), ds)

            ! d\eta/dy in mixed spectral space
            call diffy(svor(:, :, :, 2), as)

            ! d\xi/dx + d\eta/dy
            ds = ds + as

            call integrate_decomposed_field(ds)


            call field_decompose_semi_spectral(ds)
            call field_combine_physical(ds, fp)

            ! get surface zeta in physical space
            svor(0, :, :, 3) = szeta(0, :, :)
            call ifft2d(svor(0, :, :, 3), zeta(0, :, :))
            vor(0, :, :, 3) = zeta(0, :, :)

            ! get complete zeta in physical space
            do iz = 1, nz
                vor(iz, :, :, 3) = zeta(0, :, :) - fp(iz, :, :)
            enddo

            call fftxyp2s(vor(:, :, :, 3), svor(:, :, :, 3))


            ! correction
            ! ds is semi-spectral here
            do iz = 0, nz
                z = lower(3) + dble(iz) * dx(3)
                do kx = box%lo(1), box%hi(1)
                    do ky = box%lo(2), box%hi(2)
                        y = szeta(1, ky, kx) - svor(nz, ky, kx, 3)
                        kl = kh(ky, kx)
                        if (kx == 0 .and. ky == 0) then
                            kl = 1.0d0
                            y = 0.0d0
                        endif
!                         psi = (y/K) * exp{K(z - z_max)} * [1 + exp{-2K(z - z_min)}] / [1 - exp{-2KL_z}]
                        psi(iz, ky, kx) = (y/kl) * dexp(kl * (z - upper(3))) &
                                        * (1.0d0 + dexp(-2.0d0* kl * (z - lower(3)))) &
                                        / (1.0d0 - dexp(-2.0d0 *kl * extent(3)))
                        psi_z(iz, ky, kx) = y * phip(iz, kx, ky)
                    enddo
                enddo
            enddo

            call diffx(psi, psi_x)
            call diffy(psi, psi_y)

            call field_decompose_semi_spectral(psi_x)
            call field_decompose_semi_spectral(psi_y)
            svor(:, :, :, 1) = svor(:, :, :, 1) + psi_x
            svor(:, :, :, 2) = svor(:, :, :, 2) + psi_y

            svor(:, :, :, 3) = svor(:, :, :, 3) + psi_z

            ds = svor(:, :, :, 3)
            call fftxys2p(ds, vor(:, :, :, 3))

        end subroutine combine_zeta

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine integrate_decomposed_field(fs)
            double precision, intent(inout) :: fs(0:nz, box%lo(2):box%hi(2), &
                                                        box%lo(1):box%hi(1))    ! mixed-spectral space
            double precision                :: es(0:nz, box%lo(2):box%hi(2), &
                                                        box%lo(1):box%hi(1))
            integer                         :: kx, ky, kz, iz

            ! harmonic part
            !$omp parallel do private(iz)  default(shared)
            do iz = 0, nz
                es(iz, :, :) = fs(0, :, :) * psim(iz, :, :) + fs(nz, :, :) * psip(iz, :, :)
            enddo
            !$omp end parallel do

            fs(0,  :, :) = zero
            fs(nz, :, :) = zero

            !$omp parallel do
            do kz = 1, nz-1
                fs(kz, :, :) = rkzi(kz) * fs(kz, :, :)
            enddo
            !$omp end parallel do

            !$omp parallel do collapse(2) private(kx, ky)
            do kx = box%lo(1), box%hi(1)
                do ky = box%lo(2), box%hi(2)
                    call dct(1, nz, fs(0:nz, ky, kx), ztrig, zfactors)
                enddo
            enddo
            !$omp end parallel do

            !$omp parallel do
            do iz = 1, nz
                fs(iz, :, :) = fs(0, :, :) - fs(iz, :, :)
            enddo
            !$omp end parallel do

            fs(0, :, :) = zero

            fs = fs + es

        end subroutine integrate_decomposed_field

end module field_zeta
