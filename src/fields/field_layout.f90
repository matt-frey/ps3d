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
            procedure (get_field_sum),  deferred :: get_sum
            procedure (get_field_mean), deferred :: get_mean
            procedure :: get_rms => get_field_rms
            procedure :: get_absmax => get_field_absmax

            ! Field operations:
            procedure (field_diffz), deferred :: diffz

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

        function get_field_sum(this, ff, l_allreduce) result(res)
            use mpi_layout, only : box
            import :: flayout_t
            class (flayout_t), intent(in) :: this
            double precision,  intent(in) :: ff(box%lo(3):box%hi(3), &
                                                box%lo(2):box%hi(2), &
                                                box%lo(1):box%hi(1))
            logical,           intent(in) :: l_allreduce
            double precision              :: res
        end function

        function get_field_mean(this, ff, l_allreduce) result(mean)
            use mpi_layout, only : box
            import :: flayout_t
            class (flayout_t), intent(in) :: this
            double precision,  intent(in) :: ff(box%lo(3):box%hi(3), &
                                                box%lo(2):box%hi(2), &
                                                box%lo(1):box%hi(1))
            logical,           intent(in) :: l_allreduce
            double precision              :: mean
        end function

        subroutine field_diffz(this, fs, ds)
            use parameters, only : nz
            use mpi_layout, only : box
            import :: flayout_t
            class (flayout_t), intent(in)  :: this
            double precision,  intent(in)  :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision,  intent(out) :: ds(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
        end subroutine

    end interface

contains

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
