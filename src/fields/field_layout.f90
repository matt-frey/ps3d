module field_layout
    use mpi_layout
    use mpi_environment
    use constants, only : zero, f12, f13, one, two
    use parameters, only : nx, ny, nz, extent, dx, lower, upper, ncell
    use sta3dfft, only : k2l2, k2l2i, initialise_fft, rkz
    use mpi_collectives, only : mpi_blocking_reduce
    use mpi_utils, only : mpi_check_for_error
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
        double precision, allocatable :: dphim(:, :, :)     ! dphi_{-}/dz
        double precision, allocatable :: dphip(:, :, :)     ! dphi_{+}/dz

    contains
        procedure :: initialise => m_initialise
        procedure :: finalise => m_finalise

        ! Internal routine allowed to be called in child classes
        procedure :: init_decomposition

        ! Axes
        procedure :: get_x_axis => m_get_x_axis
        procedure :: get_y_axis => m_get_y_axis
        procedure (m_get_z_axis), deferred :: get_z_axis

        ! Field decompositions:
        procedure (m_decompose_physical),      deferred :: decompose_physical
        procedure (m_combine_physical),        deferred :: combine_physical
        procedure (m_decompose_semi_spectral), deferred :: decompose_semi_spectral
        procedure (m_combine_semi_spectral),   deferred :: combine_semi_spectral

        procedure, private :: set_hyperbolic_functions

        ! Field diagnostics:
        procedure (get_field_local_sum),  deferred :: get_local_sum
        procedure :: get_sum => get_field_sum
        procedure :: get_local_mean => get_field_local_mean
        procedure :: get_mean => get_field_mean
        procedure :: get_rms => get_field_rms
        procedure :: get_absmax => get_field_absmax

        ! Field operations:
        procedure (m_diffz), deferred :: diffz
        procedure (m_calc_decomposed_mean), deferred :: calc_decomposed_mean
        procedure (m_adjust_decomposed_mean), deferred :: adjust_decomposed_mean

        ! Specific routines:
        procedure (m_vertvel), deferred :: vertvel
        procedure (m_zinteg),  deferred :: zinteg

        procedure (m_zdiffuse), deferred :: zdiffuse

    end type layout_t

    interface
        function m_get_z_axis(this) result(get_z_axis)
            use parameters, only : nz
            import :: layout_t
            class (layout_t), intent(in) :: this
            double precision             :: get_z_axis(0:nz)
        end function

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

        function get_field_local_sum(this, ff) result(res)
            use mpi_layout, only : box
            import :: layout_t
            class (layout_t), intent(in) :: this
            double precision, intent(in) :: ff(box%lo(3):box%hi(3), &
                                               box%lo(2):box%hi(2), &
                                               box%lo(1):box%hi(1))
            double precision             :: res
        end function

        subroutine m_diffz(this, fs, ds, l_decomposed)
            use parameters, only : nz
            use mpi_layout, only : box
            import :: layout_t
            class (layout_t), intent(in)  :: this
            double precision, intent(in)  :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision, intent(out) :: ds(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            logical,          intent(in)  :: l_decomposed
        end subroutine

        function m_calc_decomposed_mean(this, fs) result(savg)
            use parameters, only : nz
            use mpi_layout, only : box
            import :: layout_t
            class (layout_t), intent(in) :: this
            double precision, intent(in) :: fs(0:nz,                &
                                               box%lo(2):box%hi(2), &
                                               box%lo(1):box%hi(1))
            double precision             :: savg
        end function

        subroutine m_adjust_decomposed_mean(this, fs, avg)
            use parameters, only : nz
            use mpi_layout, only : box
            import :: layout_t
            class (layout_t), intent(in) :: this
            double precision, intent(inout) :: fs(0:nz,                &
                                                  box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1))
            double precision, intent(in)    :: avg
        end subroutine

        subroutine m_vertvel(this, ds, es)
            use mpi_layout, only : box
            import :: layout_t
            class (layout_t), intent(in) :: this
            double precision, intent(inout) :: ds(box%lo(3):box%hi(3), &
                                                  box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1))
            double precision, intent(out) :: es(box%lo(3):box%hi(3), &
                                                box%lo(2):box%hi(2), &
                                                box%lo(1):box%hi(1))
        end subroutine

        subroutine m_zinteg(this, f, g, noavg)
            use parameters, only : nz
            import :: layout_t
            class (layout_t), intent(in)  :: this
            double precision, intent(in)  :: f(0:nz)
            double precision, intent(out) :: g(0:nz)
            logical,          intent(in)  :: noavg
        end subroutine

        subroutine m_zdiffuse(this, fs, dt, alpha_h, alpha_v)
            use mpi_layout, only : box
            use parameters, only : nz
            import :: layout_t
            class (layout_t), intent(in)    :: this
            double precision, intent(inout) :: fs(0:nz,                &
                                                  box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1))
            double precision, intent(in)    :: dt
            double precision, intent(in)    :: alpha_h
            double precision, intent(in)    :: alpha_v
        end subroutine
    end interface

contains

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine m_initialise(this)
        class(layout_t), intent(inout) :: this

        call this%init_decomposition

    end subroutine m_initialise

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine init_decomposition(this)
        class(layout_t), intent(inout) :: this
        double precision               :: z(0:nz), zm(0:nz), zp(0:nz)
        double precision               :: phip00(0:nz)
        integer                        :: kx, ky, iz

        !------------------------------------------------------------------
        ! Ensure FFT module is initialised:
        ! (this call does nothing if already initialised)
        call initialise_fft(extent)

        !---------------------------------------------------------------------
        !Define zm = zmax - z, zp = z - zmin
        z = this%get_z_axis()
        !$omp parallel do private(z)
        do iz = 0, nz
            zm(iz) = upper(3) - z(iz)
            zp(iz) = z(iz) - lower(3)
        enddo
        !$omp end parallel do

        !---------------------------------------------------------------------
        !Hyperbolic functions used for solutions of Laplace's equation:
        allocate(this%phim(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(this%phip(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(this%dphim(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
        allocate(this%dphip(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
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

            this%dphim(:, 0, 0) = - one / extent(3)
            this%dphip(:, 0, 0) =   one / extent(3)

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
                            phip00(0:nz),           &
                            nz+1,                   &
                            MPI_DOUBLE_PRECISION,   &
                            MPI_SUM,                &
                            world%comm,             &
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

    end subroutine init_decomposition

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine m_finalise(this)
        class(layout_t), intent(inout) :: this

        deallocate(this%gamtop)
        deallocate(this%gambot)
        deallocate(this%phim)
        deallocate(this%phip)
        deallocate(this%dphim)
        deallocate(this%dphip)
        deallocate(this%thetam)
        deallocate(this%thetap)
        deallocate(this%dthetam)

    end subroutine m_finalise

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function m_get_x_axis(this)
        class (layout_t), intent(in) :: this
        double precision             :: m_get_x_axis(0:nx-1)
        integer                      :: i

        do i = 0, nx-1
            m_get_x_axis(i) = lower(1) + dble(i) * dx(1)
        enddo

    end function m_get_x_axis

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function m_get_y_axis(this)
        class (layout_t), intent(in) :: this
        double precision             :: m_get_y_axis(0:ny-1)
        integer                      :: i

        do i = 0, ny-1
            m_get_y_axis(i) = lower(2) + dble(i) * dx(2)
        enddo

    end function m_get_y_axis

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! for kx > 0 and ky >= 0 or kx >= 0 and ky > 0
    subroutine set_hyperbolic_functions(this, kx, ky, zm, zp)
        class(layout_t),  intent(inout) :: this
        integer,          intent(in)    :: kx, ky
        double precision, intent(in)    :: zm(0:nz), zp(0:nz)
        double precision                :: R(0:nz), Q(0:nz), k2ifac
        double precision                :: ef, em(0:nz), ep(0:nz), Lm(0:nz), Lp(0:nz)
        double precision                :: fac, div, kl

        kl = sqrt(k2l2(ky, kx))
        fac = kl * extent(3)
        ef = exp(- fac)
#ifndef NDEBUG
        ! To avoid "Floating-point exception - erroneous arithmetic operation"
        ! when ef is really small.
        ef = max(ef, sqrt(tiny(ef)))
#endif
        div = one / (one - ef**2)
        k2ifac = f12 * k2l2i(ky, kx)

        Lm = kl * zm
        Lp = kl * zp

        ep = exp(- Lp)
        em = exp(- Lm)

#ifndef NDEBUG
        ! To avoid "Floating-point exception - erroneous arithmetic operation"
        ! when ep and em are really small.
        ep = max(ep, sqrt(tiny(ep)))
        em = max(em, sqrt(tiny(em)))
#endif

        this%phim(:, ky, kx) = div * (ep - ef * em)
        this%phip(:, ky, kx) = div * (em - ef * ep)

        this%dphim(:, ky, kx) = - kl * div * (ep + ef * em)
        this%dphip(:, ky, kx) =   kl * div * (em + ef * ep)

        Q = div * (one + ef**2)
        R = div * two * ef

        this%thetam(:, ky, kx) = k2ifac * (R * Lm * this%phip(:, ky, kx) - &
                                           Q * Lp * this%phim(:, ky, kx))
        this%thetap(:, ky, kx) = k2ifac * (R * Lp * this%phim(:, ky, kx) - &
                                           Q * Lm * this%phip(:, ky, kx))

        this%dthetam(:, ky, kx) = - k2ifac * ((Q * Lp - one) * this%dphim(:, ky, kx) - &
                                                      R * Lm * this%dphip(:, ky, kx))
        this%dthetap(:, ky, kx) = - k2ifac * ((Q * Lm - one) * this%dphip(:, ky, kx) - &
                                                      R * Lp * this%dphim(:, ky, kx))
    end subroutine set_hyperbolic_functions

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_field_sum(this, ff, l_allreduce) result(res)
        class (layout_t), intent(in) :: this
        double precision, intent(in) :: ff(box%lo(3):box%hi(3), &
                                           box%lo(2):box%hi(2), &
                                           box%lo(1):box%hi(1))
        logical,          intent(in) :: l_allreduce
        double precision             :: res

        res = this%get_local_sum(ff)

        if (l_allreduce) then
            call MPI_Allreduce(MPI_IN_PLACE,            &
                                res,                     &
                                1,                       &
                                MPI_DOUBLE_PRECISION,    &
                                MPI_SUM,                 &
                                world%comm,              &
                                world%err)

            call mpi_check_for_error(world, &
                "in MPI_Allreduce of field_layout::get_field_sum.")
        else
            call mpi_blocking_reduce(res, MPI_SUM, world)
        endif

    end function get_field_sum

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_field_local_mean(this, ff) result(res)
        class (layout_t), intent(in) :: this
        double precision, intent(in) :: ff(box%lo(3):box%hi(3), &
                                           box%lo(2):box%hi(2), &
                                           box%lo(1):box%hi(1))
        double precision             :: res

        res = this%get_local_sum(ff) / dble(ncell)

    end function get_field_local_mean

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_field_mean(this, ff, l_allreduce) result(mean)
        class (layout_t), intent(in) :: this
        double precision, intent(in) :: ff(box%lo(3):box%hi(3), &
                                           box%lo(2):box%hi(2), &
                                           box%lo(1):box%hi(1))
        logical,          intent(in) :: l_allreduce
        double precision             :: mean

        ! (divide by ncell since lower and upper edge weights are halved)
        mean = this%get_sum(ff, l_allreduce) / dble(ncell)

        end function get_field_mean

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_field_rms(this, ff, l_allreduce) result(rms)
        class (layout_t), intent(in) :: this
        double precision, intent(in) :: ff(box%lo(3):box%hi(3), &
                                           box%lo(2):box%hi(2), &
                                           box%lo(1):box%hi(1))
        logical,          intent(in) :: l_allreduce
        double precision             :: fsq(box%lo(3):box%hi(3), &
                                            box%lo(2):box%hi(2), &
                                            box%lo(1):box%hi(1))
        double precision             :: rms

        fsq = ff ** 2

        rms = this%get_mean(fsq, l_allreduce)

        rms = sqrt(rms)

    end function get_field_rms

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_field_absmax(this, ff, l_allreduce) result(absmax)
        class (layout_t), intent(in) :: this
        double precision, intent(in) :: ff(box%lo(3):box%hi(3), &
                                           box%lo(2):box%hi(2), &
                                           box%lo(1):box%hi(1))
        logical,          intent(in) :: l_allreduce
        double precision             :: absmax

        absmax = maxval(abs(ff(box%lo(3):box%hi(3),   &
                               box%lo(2):box%hi(2),   &
                               box%lo(1):box%hi(1))))


        if (l_allreduce) then
            call MPI_Allreduce(MPI_IN_PLACE,            &
                               absmax,                  &
                               1,                       &
                               MPI_DOUBLE_PRECISION,    &
                               MPI_MAX,                 &
                               world%comm,              &
                               world%err)

            call mpi_check_for_error(world, &
                "in MPI_Allreduce of field_diagnostics::get_field_absmax.")
        else
            call mpi_blocking_reduce(absmax, MPI_MAX, world)
        endif

    end function get_field_absmax

end module field_layout
