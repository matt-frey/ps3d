module field_filter
    use options, only : verbose
    use mpi_utils, only : mpi_print
    implicit none

    type, abstract :: filter_t

    contains
        procedure :: initialise => m_initialise
        procedure :: finalise => m_finalise


        procedure (m_apply),   deferred :: apply
        procedure (m_init_hou_and_li), private, deferred :: init_hou_and_li
        procedure (m_init_23rd_rule),  private, deferred :: init_23rd_rule
        procedure (m_init_23rd_rule_circular),  private, deferred :: init_23rd_rule_circular
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

        subroutine m_init_hou_and_li(this, l_disable_vertical)
            import filter_t
            class(filter_t), intent(inout) :: this
            logical,         intent(in)    :: l_disable_vertical
        end subroutine m_init_hou_and_li

        subroutine m_init_23rd_rule(this, l_disable_vertical)
            import filter_t
            class(filter_t), intent(inout) :: this
            logical,         intent(in)    :: l_disable_vertical
        end subroutine m_init_23rd_rule

        subroutine m_init_23rd_rule_circular(this)
            import filter_t
            class(filter_t), intent(inout) :: this
        end subroutine m_init_23rd_rule_circular

        subroutine m_init_none(this)
            import filter_t
            class(filter_t), intent(inout) :: this
        end subroutine m_init_none
    end interface

contains

    subroutine m_initialise(this, method)
        character(8),    intent(in)    :: method
        class(filter_t), intent(inout) :: this
        character(24)                  :: used_method

        !----------------------------------------------------------
        !Define de-aliasing filter:
        select case (method)
            case ("Hou & Li")
                used_method = method
                call this%init_hou_and_li(l_disable_vertical=.false.)
            case ("2/3-rule")
                call this%init_23rd_rule(l_disable_vertical=.false.)
                used_method = method
            case ("2/3-rule (circular)")
                call this%init_23rd_rule_circular
                used_method = method
            case ("Hou & Li (no vertical)")
                used_method = method
                call this%init_hou_and_li(l_disable_vertical=.true.)
            case ("2/3-rule (no vertical)")
                call this%init_23rd_rule(l_disable_vertical=.true.)
                used_method = method
            case ("none")
                call this%init_none
                used_method = "no"
            case default
                call this%init_hou_and_li(l_disable_vertical=.false.)
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
