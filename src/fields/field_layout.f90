module field_layout
    use mpi_layout
    use mpi_environment
    use constants, only : zero, f12, f13, one, two
    use parameters, only : nz, extent, dx, lower, upper
    use sta3dfft, only : k2l2, k2l2i, initialise_fft
    implicit none

    type, abstract :: layout_t

        double precision, allocatable :: gamtop(:), gambot(:)

        ! See for definitions in
        ! Dritschel D, Frey M. The stability of inviscid Beltrami flow between parallel free-slip impermeable
        ! boundaries. Journal of Fluid Mechanics. 2023;954:A31. doi:10.1017/jfm.2022.1007
        double precision, allocatable :: thetam(:, :, :)    ! theta_{-}         (eq. 3.10)
        double precision, allocatable :: thetap(:, :, :)    ! theta_{+}         (eq. 3.11)
        double precision, allocatable :: dthetam(:, :, :)   ! dtheta_{-}/dz
        double precision, allocatable :: dthetap(:, :, :)   ! dtheta_{+}/dz
        double precision, allocatable :: phim(:, :, :)      ! phi_{-}           (eq. 3.4a)
        double precision, allocatable :: phip(:, :, :)      ! phi_{+}           (eq. 3.4b)
#ifdef ENABLE_BUOYANCY
        double precision, allocatable :: dphim(:, :, :)     ! dphi_{-}/dz
        double precision, allocatable :: dphip(:, :, :)     ! dphi_{+}/dz
#endif

    contains
        procedure :: initialise => m_initialise
        procedure :: finalise => m_finalise

        ! Field decompositions:
        procedure (m_decompose_physical),      deferred :: decompose_physical
        procedure (m_combine_physical),        deferred :: combine_physical
        procedure (m_decompose_semi_spectral), deferred :: decompose_semi_spectral
        procedure (m_combine_semi_spectral),   deferred :: combine_semi_spectral

        procedure, private :: set_hyperbolic_functions

    end type layout_t

    interface
        subroutine m_decompose_physical(this, fc, sf)
            use parameters, only : nz
            use mpi_layout, only : box
            import :: layout_t
            class (layout_t), intent(in)  :: this
            double precision, intent(in)  :: fc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision, intent(out) :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        end subroutine

        subroutine m_combine_physical(this, sf, fc)
            use parameters, only : nz
            use mpi_layout, only : box
            import :: layout_t
            class (layout_t), intent(in)  :: this
            double precision, intent(in)  :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision, intent(out) :: fc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        end subroutine

        subroutine m_decompose_semi_spectral(this, sfc)
            use parameters, only : nz
            use mpi_layout, only : box
            import :: layout_t
            class (layout_t), intent(in)    :: this
            double precision, intent(inout) :: sfc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        end subroutine

        subroutine m_combine_semi_spectral(this, sf)
            use parameters, only : nz
            use mpi_layout, only : box
            import :: layout_t
            class (layout_t), intent(in)    :: this
            double precision, intent(inout) :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        end subroutine
    end interface

contains

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine m_initialise(this)
        class(layout_t), intent(inout) :: this
        double precision               :: z, zm(0:nz), zp(0:nz)
        double precision               :: phip00(0:nz)
        integer                        :: kx, ky, iz

        !------------------------------------------------------------------
        ! Ensure FFT module is initialised:
        ! (this call does nothing if already initialised)
        call initialise_fft(extent)

        !---------------------------------------------------------------------
        !Define zm = zmax - z, zp = z - zmin
        !$omp parallel do private(z)
        do iz = 0, nz
!             z = lower(3) + dx(3) * dble(iz) !FIXME
            zm(iz) = upper(3) - z
            zp(iz) = z - lower(3)
        enddo
        !$omp end parallel do

        !---------------------------------------------------------------------
        !Hyperbolic functions used for solutions of Laplace's equation:
        allocate(this%phim(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(this%phip(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
#ifdef ENABLE_BUOYANCY
        allocate(this%dphim(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(this%dphip(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
#endif
        allocate(this%thetam(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(this%thetap(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(this%dthetam(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(this%dthetap(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

        do kx = box%lo(1), box%hi(1)
            do ky = max(1, box%lo(2)), box%hi(2)
                call this%set_hyperbolic_functions(kx, ky, zm, zp)
            enddo
        enddo

        ! ky = 0
        if (box%lo(2) == 0) then
            do kx = max(1, box%lo(1)), box%hi(1)
                call this%set_hyperbolic_functions(kx, 0, zm, zp)
            enddo
        endif

        phip00 = zero
        if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
            !$omp parallel workshare
            ! kx = ky = 0
            this%phim(:, 0, 0) = zm / extent(3)
            this%phip(:, 0, 0) = zp / extent(3)

#ifdef ENABLE_BUOYANCY
            this%dphim(:, 0, 0) = - one / extent(3)
            this%dphip(:, 0, 0) =   one / extent(3)
#endif

            this%thetam(:, 0, 0) = zero
            this%thetap(:, 0, 0) = zero

            this%dthetam(:, 0, 0) = zero
            this%dthetap(:, 0, 0) = zero

            phip00 = this%phip(:, 0, 0)
            !$omp end parallel workshare
        endif

        !---------------------------------------------------------------------
        !Define gamtop as the integral of phip(iz, 0, 0) with zero average:
        allocate(this%gamtop(0:nz))
        allocate(this%gambot(0:nz))

        call MPI_Allreduce(MPI_IN_PLACE,            &
                            phip00(0:nz),            &
                            nz+1,                    &
                            MPI_DOUBLE_PRECISION,    &
                            MPI_SUM,                 &
                            world%comm,              &
                            world%err)

        !$omp parallel workshare
        this%gamtop = f12 * extent(3) * (phip00 ** 2 - f13)
        !$omp end parallel workshare

        !$omp parallel do
        do iz = 0, nz
            this%gambot(iz) = this%gamtop(nz-iz)
        enddo
        !$omp end parallel do
        !Here gambot is the complement of gamtop.

    end subroutine m_initialise

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine m_finalise(this)
        class(layout_t), intent(inout) :: this

        deallocate(this%gamtop)
        deallocate(this%gambot)
        deallocate(this%phim)
        deallocate(this%phip)
#ifdef ENABLE_BUOYANCY
        deallocate(this%dphim)
        deallocate(this%dphip)
#endif
        deallocate(this%thetam)
        deallocate(this%thetap)
        deallocate(this%dthetam)

    end subroutine m_finalise

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! for kx > 0 and ky >= 0 or kx >= 0 and ky > 0
    subroutine set_hyperbolic_functions(this, kx, ky, zm, zp)
        class(layout_t),  intent(inout) :: this
        integer,          intent(in)    :: kx, ky
        double precision, intent(in)    :: zm(0:nz), zp(0:nz)
        double precision                :: R(0:nz), Q(0:nz), k2ifac
#ifndef ENABLE_BUOYANCY
        double precision                :: dphim(0:nz), dphip(0:nz)
#endif
        double precision                :: ef, em(0:nz), ep(0:nz), Lm(0:nz), Lp(0:nz)
        double precision                :: fac, div, kl

        kl = sqrt(k2l2(ky, kx))
        fac = kl * extent(3)
        ef = dexp(- fac)
#ifndef NDEBUG
        ! To avoid "Floating-point exception - erroneous arithmetic operation"
        ! when ef is really small.
        ef = max(ef, sqrt(tiny(ef)))
#endif
        div = one / (one - ef**2)
        k2ifac = f12 * k2l2i(ky, kx)

        Lm = kl * zm
        Lp = kl * zp

        ep = dexp(- Lp)
        em = dexp(- Lm)

#ifndef NDEBUG
        ! To avoid "Floating-point exception - erroneous arithmetic operation"
        ! when ep and em are really small.
        ep = max(ep, sqrt(tiny(ep)))
        em = max(em, sqrt(tiny(em)))
#endif

        this%phim(:, ky, kx) = div * (ep - ef * em)
        this%phip(:, ky, kx) = div * (em - ef * ep)

#ifdef ENABLE_BUOYANCY
        this%dphim(:, ky, kx) = - kl * div * (ep + ef * em)
        this%dphip(:, ky, kx) =   kl * div * (em + ef * ep)
#else
        dphim = - kl * div * (ep + ef * em)
        dphip =   kl * div * (em + ef * ep)
#endif

        Q = div * (one + ef**2)
        R = div * two * ef

        this%thetam(:, ky, kx) = k2ifac * (R * Lm * this%phip(:, ky, kx) - &
                                           Q * Lp * this%phim(:, ky, kx))
        this%thetap(:, ky, kx) = k2ifac * (R * Lp * this%phim(:, ky, kx) - &
                                           Q * Lm * this%phip(:, ky, kx))

#ifdef ENABLE_BUOYANCY
        this%dthetam(:, ky, kx) = - k2ifac * ((Q * Lp - one) * this%dphim(:, ky, kx) - &
                                                      R * Lm * this%dphip(:, ky, kx))
        this%dthetap(:, ky, kx) = - k2ifac * ((Q * Lm - one) * this%dphip(:, ky, kx) - &
                                                      R * Lp * this%dphim(:, ky, kx))
#else
        this%dthetam(:, ky, kx) = - k2ifac * ((Q * Lp - one) * dphim - R * Lm * dphip)
        this%dthetap(:, ky, kx) = - k2ifac * ((Q * Lm - one) * dphip - R * Lp * dphim)
#endif
    end subroutine set_hyperbolic_functions

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#ifdef ENABLE_BUOYANCY
    !Calculates df/dz for a field f in mixed-spectral space
    !Here fs = f, ds = df/dz. Both fields are in mixed-spectral space.
    ! fs - mixed-spectral space
    ! ds - derivative linear part
    ! as - derivative sine part
    subroutine diffz(fs, ds)
        double precision, intent(in)  :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision, intent(out) :: ds(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        double precision              :: as(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        integer                       :: kz, iz

        !Calculate the derivative of the linear part (ds) in semi-spectral space:
        !$omp parallel do private(iz)  default(shared)
        do iz = 0, nz
            ds(iz, :, :) = fs(0, :, :) * dphim(iz, :, :) + fs(nz, :, :) * dphip(iz, :, :)
        enddo
        !$omp end parallel do

        ! Calculate d/dz of this sine series:
        !$omp parallel workshare
        as(0, :, :) = zero
        !$omp end parallel workshare
        !$omp parallel do private(kz)  default(shared)
        do kz = 1, nz-1
            as(kz, :, :) = rkz(kz) * fs(kz, :, :)
        enddo
        !$omp end parallel do
        !$omp parallel workshare
        as(nz, :, :) = zero
        !$omp end parallel workshare

        !FFT these quantities back to semi-spectral space:
        call fftcosine(as)

        ! Combine vertical derivative given the sine (as) and linear (ds) parts:
        !omp parallel workshare
        ds = ds + as
        !omp end parallel workshare

        call decompose_semi_spectral(ds)

    end subroutine diffz
#endif

end module field_layout
