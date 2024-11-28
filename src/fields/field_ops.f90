module field_ops
    use parameters, only : nx, ny, nz, dx, lower
    use mpi_environment
    use mpi_layout, only : box
    use parameters, only : ncell
    use mpi_collectives, only : mpi_blocking_reduce
    use mpi_utils, only : mpi_check_for_error
    implicit none

    type, abstract :: ops_t

    contains
        procedure :: initialise => m_initialise
        procedure :: finalise => m_finalise

        ! Axes
        procedure :: get_x_axis => m_get_x_axis
        procedure :: get_y_axis => m_get_y_axis
        procedure (m_get_z_axis), deferred :: get_z_axis

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
        procedure (m_apply_filter), deferred :: apply_filter

        ! Specific routines:
        procedure (m_vertvel), deferred :: vertvel
        procedure (m_zinteg),  deferred :: zinteg
    end type

    interface
        function m_get_z_axis(this) result(get_z_axis)
            use parameters, only : nz
            import :: ops_t
            class (ops_t), intent(in)  :: this
            double precision :: get_z_axis(0:nz)
        end function

        function get_field_local_sum(this, ff) result(res)
            use mpi_layout, only : box
            import :: ops_t
            class (ops_t),    intent(in) :: this
            double precision, intent(in) :: ff(box%lo(3):box%hi(3), &
                                               box%lo(2):box%hi(2), &
                                               box%lo(1):box%hi(1))
            double precision             :: res
        end function

        subroutine m_diffz(this, fs, ds)
            use parameters, only : nz
            use mpi_layout, only : box
            import :: ops_t
            class (ops_t),    intent(in)  :: this
            double precision, intent(in)  :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision, intent(out) :: ds(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        end subroutine

        function m_calc_decomposed_mean(this, fs) result(savg)
            use mpi_layout, only : box
            import :: ops_t
            class (ops_t),    intent(in) :: this
            double precision, intent(in) :: fs(box%lo(3):box%hi(3), &
                                               box%lo(2):box%hi(2), &
                                               box%lo(1):box%hi(1))
            double precision             :: savg
        end function

        subroutine m_adjust_decomposed_mean(this, fs, avg)
            use mpi_layout, only : box
            import :: ops_t
            class (ops_t),    intent(in)    :: this
            double precision, intent(inout) :: fs(box%lo(3):box%hi(3), &
                                                  box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1))
            double precision, intent(in)    :: avg
        end subroutine

        subroutine m_apply_filter(this, fs)
            use mpi_layout, only : box
            import :: ops_t
            class (ops_t),    intent(in)    :: this
            double precision, intent(inout) :: fs(box%lo(3):box%hi(3), &
                                                  box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1))
        end subroutine

        subroutine m_vertvel(this, ds, es)
            use mpi_layout, only : box
            import :: ops_t
            class (ops_t),    intent(in)    :: this
            double precision, intent(inout) :: ds(box%lo(3):box%hi(3), &
                                                  box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1))
            double precision, intent(out) :: es(box%lo(3):box%hi(3), &
                                                box%lo(2):box%hi(2), &
                                                box%lo(1):box%hi(1))
        end subroutine

        subroutine m_zinteg(this, f, g, noavg)
            use parameters, only : nz
            import :: ops_t
            class (ops_t),    intent(in)  :: this
            double precision, intent(in)  :: f(0:nz)
            double precision, intent(out) :: g(0:nz)
            logical,          intent(in)  :: noavg
        end subroutine
    end interface

contains

    subroutine m_initialise(this)
        class (ops_t), intent(inout)  :: this

        ! The base class does nothing here

    end subroutine m_initialise

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine m_finalise(this)
        class (ops_t), intent(inout)  :: this

        ! The base class does nothing here

    end subroutine m_finalise

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function m_get_x_axis(this)
        class (ops_t), intent(in) :: this
        double precision             :: m_get_x_axis(0:nx-1)
        integer                      :: i

        do i = 0, nx-1
            m_get_x_axis(i) = lower(1) + dble(i) * dx(1)
        enddo

    end function m_get_x_axis

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function m_get_y_axis(this)
        class (ops_t), intent(in) :: this
        double precision             :: m_get_y_axis(0:ny-1)
        integer                      :: i

        do i = 0, nz-1
            m_get_y_axis(i) = lower(2) + dble(i) * dx(2)
        enddo

    end function m_get_y_axis

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_field_sum(this, ff, l_allreduce) result(res)
        class (ops_t),    intent(in) :: this
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
        class (ops_t),    intent(in) :: this
        double precision, intent(in) :: ff(box%lo(3):box%hi(3), &
                                           box%lo(2):box%hi(2), &
                                           box%lo(1):box%hi(1))
        double precision             :: res

        res = this%get_local_sum(ff) / dble(ncell)

    end function get_field_local_mean

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_field_mean(this, ff, l_allreduce) result(mean)
        class (ops_t),    intent(in) :: this
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
        class (ops_t),    intent(in) :: this
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
        class (ops_t),    intent(in) :: this
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

end module field_ops
