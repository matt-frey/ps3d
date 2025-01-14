module field_filter
    use options, only : verbose
    use mpi_utils, only : mpi_print
    implicit none

    type, abstract :: filter_t

    contains
        procedure :: initialise => m_initialise
        procedure :: finalise => m_finalise


        procedure (m_apply),   deferred :: apply
        procedure (m_apply2d), deferred :: apply2d
        procedure (m_init_hou_and_li), private, deferred :: init_hou_and_li
        procedure (m_init_23rd_rule),  private, deferred :: init_23rd_rule
        procedure (m_init_none),       private, deferred :: init_none
    end type

    interface
        subroutine m_apply(this, fs)
            use mpi_layout, only : box
            import filter_t
            class(filter_t),  intent(in)    :: this
            double precision, intent(inout) :: fs(box%lo(3):box%hi(3), &
                                                  box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1))
        end subroutine m_apply

        subroutine m_apply2d(this, fs)
            use mpi_layout, only : box
            import filter_t
            class(filter_t),  intent(in)    :: this
            double precision, intent(inout) :: fs(box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1))
        end subroutine m_apply2d

        subroutine m_init_hou_and_li(this)
            import filter_t
            class(filter_t), intent(inout) :: this
        end subroutine m_init_hou_and_li

        subroutine m_init_23rd_rule(this)
            import filter_t
            class(filter_t), intent(inout) :: this
        end subroutine m_init_23rd_rule

        subroutine m_init_none(this)
            import filter_t
            class(filter_t), intent(inout) :: this
        end subroutine m_init_none
    end interface

contains

    subroutine m_initialise(this, method)
        character(8),    intent(in)    :: method
        class(filter_t), intent(inout) :: this
        character(8)                   :: used_method

        !----------------------------------------------------------
        !Define de-aliasing filter:
        select case (method)
            case ("Hou & Li")
                used_method = method
                call this%init_hou_and_li
            case ("2/3-rule")
                call this%init_23rd_rule
                used_method = method
            case ("none")
                call this%init_none
                used_method = "no"
            case default
                call this%init_hou_and_li
                used_method = "Hou & Li"
        end select

#ifdef ENABLE_VERBOSE
        if (verbose) then
            call mpi_print("Using " // used_method // " de-aliasing filter.")
        endif
#endif

    end subroutine m_initialise

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine m_finalise(this)
        class(filter_t), intent(inout) :: this

    end subroutine m_finalise

end module
