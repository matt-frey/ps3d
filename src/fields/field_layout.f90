module field_layout
    use mpi_layout
    use mpi_environment
    use constants, only : zero, f12, f13, one, two
    use parameters, only : nx, ny, nz, extent, dx, lower, upper, ncell
    use sta3dfft, only : k2l2, k2l2i, initialise_fft, rkz
    use mpi_collectives, only : mpi_blocking_reduce
    use mpi_utils, only : mpi_check_for_error
    use options, only : verbose, filter_type
    use mpi_utils, only : mpi_print
    implicit none

    type, abstract :: layout_t

    contains
        procedure (m_initialise), deferred :: initialise
        procedure (m_finalise),   deferred :: finalise

        ! Axes
        procedure :: get_x_axis => m_get_x_axis
        procedure :: get_y_axis => m_get_y_axis
        procedure (m_get_z_axis), deferred :: get_z_axis

        ! Field decompositions:
        procedure :: decompose_semi_spectral => m_decompose_semi_spectral
        procedure :: combine_semi_spectral => m_combine_semi_spectral

        ! Field diagnostics:
        procedure (get_field_local_sum),  deferred :: get_local_sum
        procedure :: get_sum => get_field_sum
        procedure :: get_local_mean => get_field_local_mean
        procedure :: get_mean => get_field_mean
        procedure :: get_rms => get_field_rms
        procedure :: get_absmax => get_field_absmax

        ! Field operations:
        procedure (m_diffz), deferred :: diffz
        procedure (m_get_semi_spectral_mean), deferred :: get_semi_spectral_mean
        procedure (m_adjust_semi_spectral_mean), deferred :: adjust_semi_spectral_mean

        ! Filters:
        procedure :: set_filter
        procedure (m_init_exp_filter),    deferred :: init_exp_filter
        procedure (m_init_cutoff_filter), deferred :: init_cutoff_filter
        procedure (m_apply_filter),       deferred :: apply_filter
        procedure (m_apply_hfilter),      deferred :: apply_hfilter

        ! Specific routines:
        procedure (m_vertvel), deferred :: vertvel
        procedure (m_zinteg),  deferred :: zinteg

        procedure (m_zdiffuse),   deferred :: zdiffuse
        procedure (m_zdiffNF),    deferred :: zdiffNF

    end type layout_t

    interface
        subroutine m_initialise(this)
            import :: layout_t
            class(layout_t),  intent(inout) :: this
        end subroutine m_initialise

        subroutine m_finalise(this)
            import :: layout_t
            class(layout_t),  intent(inout) :: this
        end subroutine m_finalise

        function m_get_z_axis(this) result(get_z_axis)
            use parameters, only : nz
            import :: layout_t
            class (layout_t), intent(in) :: this
            double precision             :: get_z_axis(0:nz)
        end function

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

        function m_get_semi_spectral_mean(this, fs) result(savg)
            use parameters, only : nz
            use mpi_layout, only : box
            import :: layout_t
            class (layout_t), intent(in) :: this
            double precision, intent(in) :: fs(0:nz,                &
                                               box%lo(2):box%hi(2), &
                                               box%lo(1):box%hi(1))
            double precision             :: savg
        end function

        subroutine m_adjust_semi_spectral_mean(this, fs, avg)
            use parameters, only : nz
            use mpi_layout, only : box
            import :: layout_t
            class (layout_t), intent(in) :: this
            double precision, intent(inout) :: fs(0:nz,                &
                                                  box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1))
            double precision, intent(in)    :: avg
        end subroutine

        subroutine m_init_exp_filter(this, alpha, beta)
            import :: layout_t
            class(layout_t),  intent(inout) :: this
            double precision, intent(in)    :: alpha, beta
        end subroutine m_init_exp_filter

        subroutine m_init_cutoff_filter(this, cutoff)
            import :: layout_t
            class(layout_t),  intent(inout) :: this
            double precision, intent(in)    :: cutoff
        end subroutine m_init_cutoff_filter

        subroutine m_apply_filter(this, fs)
            use mpi_layout, only : box
            import :: layout_t
            class(layout_t),  intent(in)    :: this
            double precision, intent(inout) :: fs(box%lo(3):box%hi(3), &
                                                  box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1))
        end subroutine m_apply_filter

        subroutine m_apply_hfilter(this, fs)
            use mpi_layout, only : box
            import :: layout_t
            class(layout_t),  intent(in)    :: this
            double precision, intent(inout) :: fs(box%lo(3):box%hi(3), &
                                                  box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1))
        end subroutine m_apply_hfilter

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

        subroutine m_zdiffNF(this, fs, dt, alpha_h, alpha_v)
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

    subroutine m_decompose_semi_spectral(this, sfc)
        class (layout_t), intent(in)    :: this
        double precision, intent(inout) :: sfc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))

        ! Do nothing!

    end subroutine m_decompose_semi_spectral

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine m_combine_semi_spectral(this, sf)
        class (layout_t), intent(in)    :: this
        double precision, intent(inout) :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))

        ! Do nothing!

    end subroutine m_combine_semi_spectral

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

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    !Define de-aliasing filter:
    subroutine set_filter(this, filter)
        class (layout_t),            intent(inout) :: this
        type(filter_type), optional, intent(in)    :: filter
        character(len=6)                           :: family

        if (.not. present(filter)) then
            family = "no"
        else
            family = filter%family
        endif

        select case (family)
            case ("exp")
                call this%init_exp_filter(filter%alpha, filter%beta)
            case ("cutoff")
                call this%init_cutoff_filter(filter%cutoff)
            case ("no")
                ! do nothing
            case default
                ! do nothing
        end select

        if (verbose) then
            call mpi_print("Using " // trim(family) // " de-aliasing filter.")
        endif

    end subroutine set_filter

end module field_layout
