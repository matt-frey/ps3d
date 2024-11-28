module field_filter
    implicit none

    type, abstract :: field_filter_t

        ! Spectral filter:
        double precision, allocatable :: filt(:, :, :)

    contains
        procedure :: initialise => m_initialise
        procedure :: finalise => m_finalise

        procedure (m_init_hou_and_li), private, deferred :: init_hou_and_li
        procedure (m_init_23rd_rule),  private, deferred :: init_23rd_rule
    end type

    interface
        subroutine m_init_hou_and_li(this)
            class(field_filter_t), intent(in) :: this
        end subroutine m_init_hou_and_li

        subroutine m_init_23rd_rule(this)
            class(field_filter_t), intent(in) :: this
        end subroutine m_init_23rd_rule
    end interface

contains

    subroutine m_initialise(this)
        class(field_filter_t), intent(inout) :: this

        !----------------------------------------------------------
        !Define de-aliasing filter:

        allocate(this%filt(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

        select case (filtering)
            case ("Hou & Li")
                call this%init_hou_and_li
            case ("2/3-rule")
                call this%init_23rd_rule
            case default
                call this%init_hou_and_li
        end select

        !Ensure filter does not change domain mean:
        if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
            filt(:, 0, 0) = one
        endif

    end subroutine m_initialise

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine m_finalise(this)
        class(field_filter_t), intent(inout) :: this

        deallocate(this%filt)

    end subroutine m_finalise

end module
