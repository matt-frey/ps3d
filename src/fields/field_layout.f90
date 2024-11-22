module field_layout
    use constants, only : f12
    use parameters, only : ncell
    use mpi_layout, only : box
    use mpi_environment
    use mpi_collectives, only : mpi_blocking_reduce
    use mpi_utils, only : mpi_check_for_error
    implicit none

    type, abstract :: flayout_t
        contains
            ! Field decompositions:
            procedure (field_decompose_physical),      deferred :: decompose_physical
            procedure (field_combine_physical),        deferred :: combine_physical
            procedure (field_decompose_semi_spectral), deferred :: decompose_semi_spectral
            procedure (field_combine_semi_spectral),   deferred :: combine_semi_spectral

            ! Field diagnostics:
            procedure (get_field_local_sum),  deferred :: get_local_sum
            procedure :: get_sum => get_field_sum
            procedure :: get_local_mean => get_field_local_mean
            procedure :: get_mean => get_field_mean
            procedure :: get_rms => get_field_rms
            procedure :: get_absmax => get_field_absmax

            ! Field operations:
            procedure (field_diffz), deferred :: diffz
            procedure (calc_field_decomposed_mean), deferred :: calc_decomposed_mean
            procedure (adjust_field_decomposed_mean), deferred :: adjust_decomposed_mean
            procedure (apply_field_filter), deferred :: apply_filter

    end type flayout_t

    interface
        subroutine field_decompose_physical(this, fc, sf)
            use parameters, only : nz
            use mpi_layout, only : box
            import :: flayout_t
            class (flayout_t), intent(in)  :: this
            double precision,  intent(in)  :: fc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision,  intent(out) :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        end subroutine

        subroutine field_combine_physical(this, sf, fc)
            use parameters, only : nz
            use mpi_layout, only : box
            import :: flayout_t
            class (flayout_t), intent(in)  :: this
            double precision,  intent(in)  :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision,  intent(out) :: fc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        end subroutine

        subroutine field_decompose_semi_spectral(this, sfc)
            use parameters, only : nz
            use mpi_layout, only : box
            import :: flayout_t
            class (flayout_t), intent(in)    :: this
            double precision,  intent(inout) :: sfc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        end subroutine

        subroutine field_combine_semi_spectral(this, sf)
            use parameters, only : nz
            use mpi_layout, only : box
            import :: flayout_t
            class (flayout_t), intent(in)    :: this
            double precision,  intent(inout) :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        end subroutine

        function get_field_local_sum(this, ff) result(res)
            use mpi_layout, only : box
            import :: flayout_t
            class (flayout_t), intent(in) :: this
            double precision,  intent(in) :: ff(box%lo(3):box%hi(3), &
                                                box%lo(2):box%hi(2), &
                                                box%lo(1):box%hi(1))
            double precision              :: res
        end function

        subroutine field_diffz(this, fs, ds)
            use parameters, only : nz
            use mpi_layout, only : box
            import :: flayout_t
            class (flayout_t), intent(in)  :: this
            double precision,  intent(in)  :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision,  intent(out) :: ds(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        end subroutine

        function calc_field_decomposed_mean(this, fs) result(savg)
            use mpi_layout, only : box
            import :: flayout_t
            class (flayout_t), intent(in) :: this
            double precision,  intent(in) :: fs(box%lo(3):box%hi(3), &
                                                box%lo(2):box%hi(2), &
                                                box%lo(1):box%hi(1))
            double precision              :: savg
        end function

        subroutine adjust_field_decomposed_mean(this, fs, avg)
            use mpi_layout, only : box
            import :: flayout_t
            class (flayout_t), intent(in)    :: this
            double precision,  intent(inout) :: fs(box%lo(3):box%hi(3), &
                                                   box%lo(2):box%hi(2), &
                                                   box%lo(1):box%hi(1))
            double precision,  intent(in)    :: avg
        end subroutine

        subroutine apply_field_filter(this, fs)
            use mpi_layout, only : box
            import :: flayout_t
            class (flayout_t), intent(in)    :: this
            double precision,  intent(inout) :: fs(box%lo(3):box%hi(3), &
                                                   box%lo(2):box%hi(2), &
                                                   box%lo(1):box%hi(1))
        end subroutine

    end interface

contains

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_field_sum(this, ff, l_allreduce) result(res)
        class (flayout_t), intent(in) :: this
        double precision,  intent(in) :: ff(box%lo(3):box%hi(3), &
                                            box%lo(2):box%hi(2), &
                                            box%lo(1):box%hi(1))
        logical,           intent(in) :: l_allreduce
        double precision              :: res

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
        class (flayout_t), intent(in) :: this
        double precision,  intent(in) :: ff(box%lo(3):box%hi(3), &
                                            box%lo(2):box%hi(2), &
                                            box%lo(1):box%hi(1))
        double precision              :: res

        res = this%get_local_sum(ff) / dble(ncell)

    end function get_field_local_mean

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_field_mean(this, ff, l_allreduce) result(mean)
        class (flayout_t), intent(in) :: this
        double precision,  intent(in) :: ff(box%lo(3):box%hi(3), &
                                            box%lo(2):box%hi(2), &
                                            box%lo(1):box%hi(1))
        logical,           intent(in) :: l_allreduce
        double precision              :: mean

        ! (divide by ncell since lower and upper edge weights are halved)
        mean = this%get_sum(ff, l_allreduce) / dble(ncell)

        end function get_field_mean

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_field_rms(this, ff, l_allreduce) result(rms)
        class (flayout_t), intent(in) :: this
        double precision,  intent(in) :: ff(box%lo(3):box%hi(3), &
                                            box%lo(2):box%hi(2), &
                                            box%lo(1):box%hi(1))
        logical,           intent(in) :: l_allreduce
        double precision              :: fsq(box%lo(3):box%hi(3), &
                                             box%lo(2):box%hi(2), &
                                             box%lo(1):box%hi(1))
        double precision              :: rms

        fsq = ff ** 2

        rms = this%get_mean(fsq, l_allreduce)

        rms = sqrt(rms)

    end function get_field_rms

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_field_absmax(this, ff, l_allreduce) result(absmax)
        class (flayout_t), intent(in) :: this
        double precision,  intent(in) :: ff(box%lo(3):box%hi(3), &
                                            box%lo(2):box%hi(2), &
                                            box%lo(1):box%hi(1))
        logical,           intent(in) :: l_allreduce
        double precision              :: absmax

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
