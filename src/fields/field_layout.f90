module field_layout
    implicit none

    type, abstract :: layout_t

    contains
        ! Field decompositions:
        procedure (m_decompose_physical),      deferred :: decompose_physical
        procedure (m_combine_physical),        deferred :: combine_physical
        procedure (m_decompose_semi_spectral), deferred :: decompose_semi_spectral
        procedure (m_combine_semi_spectral),   deferred :: combine_semi_spectral

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
end module field_layout
