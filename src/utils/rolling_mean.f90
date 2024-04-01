module rolling_mean_mod
    implicit none

    ! rolling mean type for dissipation prefactor calculation:
    ! See e.g. https://en.wikipedia.org/wiki/Moving_average#Simple_moving_average
    type rolling_mean_t
        double precision, allocatable, private :: history(:)
        integer,                       private :: inew = 1
        integer,                       private :: iold = 1
        double precision,              private :: sma = 0.0d0 ! simple moving average
        logical,                       private :: l_filled = .false.
        integer,                       private :: length = 0

        contains
            procedure :: alloc
            procedure :: get_next

    end type rolling_mean_t

    contains

        subroutine alloc(self, n)
            class(rolling_mean_t), intent(inout) :: self
            integer,               intent(in)    :: n

            if (.not. allocated(self%history)) then
                allocate(self%history(n))
            endif

            self%length = n

        end subroutine alloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        function get_next(self, vnew) result(rm)
            class(rolling_mean_t), intent(inout) :: self
            double precision,      intent(in)    :: vnew
            double precision                     :: vold
            double precision                     :: rm

            if (self%l_filled) then

                ! get oldest value
                vold = self%history(self%iold)
                self%iold = mod(self%iold + 1, self%length+1)

                self%sma = self%sma + (vnew - vold) / dble(self%length)

                ! add new value
                self%history(self%inew) = vnew
                self%inew = mod(self%inew + 1, self%length+1)

            else
                ! add newest value
                self%history(self%inew) = vnew

                ! if we have not yet fully filled, we just
                ! return the mean value up to this point
                self%sma = sum(self%history(1:self%inew)) / dble(self%inew)

                self%l_filled = (self%length == self%inew)

                self%inew = mod(self%inew + 1, self%length+1)
            endif

            rm = self%sma

        end function get_next

end module rolling_mean_mod
