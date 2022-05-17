! =============================================================================
!               Module for array manipulation routines
! =============================================================================
module armanip

    contains

        subroutine array_1d_resize(src, newsize)
            double precision, intent(inout) :: src
            integer,          intent(in)    :: newsize
            double precision, allocatable   :: buf(:)
            integer                         :: oldsize

            oldsize = size(src)

            if (oldsize == newsize) then
                return
            endif

            allocate(buf(newsize))

            ! copy data
            if (newsize > oldsize) then
                buf(1:oldsize) = src
            else
                buf = src(1:newsize)
            endif

            call move_alloc(from=buf, to=src)

        end subroutine array_1d_resize

        subroutine array_2d_resize(src, newsize)
            double precision, intent(inout) :: src
            integer,          intent(in)    :: newsize
            double precision, allocatable   :: buf(:, :)
            integer                         :: oldsize, ncomp

            (/ncomp, oldsize/) = shape(src)

            if (oldsize == newsize) then
                return
            endif

            allocate(buf(ncomp, newsize))

            ! copy data
            if (newsize > oldsize) then
                buf(:, 1:oldsize) = src
            else
                buf = src(:, 1:newsize)
            endif

            call move_alloc(from=buf, to=src)

        end subroutine array_2d_resize

module armanip
