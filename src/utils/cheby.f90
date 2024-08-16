module cheby
    use constants, only : zero, one, two, pi, four
    implicit none

    contains

        subroutine init_cheby(n, x, d1, d2)
            integer,          intent(in)  :: n
            double precision, intent(out) :: x(0:n)
            double precision, intent(out) :: d1(0:n, 0:n), d2(0:n, 0:n)
            integer                       :: i, j
            double precision              :: c(0:n)
            double precision              :: rc(0:n, 0:n)
            double precision              :: fac, ss

            ! Compute Chebyshev points:
            fac = pi / dble(n)   ! acos(-1) = pi
            do i = 0, n
                x(i) = cos(fac * dble(i))
            enddo

            ! Compute weights:
            c(0) = two
            do i = 1, n-1
                c(i) = one
            enddo
            c(n) = two

            ss = one
            do i = 0, n
                c(i) = c(i) * ss
                ss = -ss
            enddo

            ! Compute differentiation matrix:
            do i = 0, n
                do j = 0, n
                    rc(i, j) = c(i) / c(j)
                    d1(i, j) = x(i) - x(j)
                enddo
            enddo
            do i = 0, n
                d1(i, i) = d1(i, i) + one
            enddo

            d1 = rc / d1
            do i = 0, n
                d1(i, i) = one - sum(d1(i, :))
            enddo

            ! Set up second-order differentiation matrix:
            d2 = matmul(d1, d1)

        end subroutine init_cheby

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! w:  Clenshaw-Curtis weights
        subroutine clencurt(n, w)
            integer,          intent(in)  :: n
            double precision, intent(out) :: w(0:n)
            double precision              :: v(0:n-2)
            double precision              :: theta
            integer                       :: i, k

            ! Initialize weights
            w = zero
            w(0) = one / (dble(n) ** 2 - one)
            w(n) = w(0)

            ! Initialize vector v
            v = one

            ! Update vector v using the given formula
            do k = 1, n/2 - 1
                do i = 1, n-1
                    theta = pi * dble(i) / dble(n)
                    v(i-1) = v(i-1) - two * cos(two * dble(k) * theta) / (four * dble(k) ** 2 - one)
                end do
            end do

            ! Finalize vector v
            do i = 1, n-1
                theta = pi * dble(i) / dble(n)
                v(i-1) = v(i-1) - cos(dble(n) * theta) / (dble(n) ** 2 - one)
            end do

            ! Compute the weights
            w(1:n-1) = two * v / dble(n)

        end subroutine clencurt

end module cheby
