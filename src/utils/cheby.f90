module cheby
    use constants, only : zero, one, two, pi, four
    use parameters, only : nz
    implicit none

    contains

        subroutine init_cheby(x, d1, d2)
            double precision, intent(out) :: x(0:nz)
            double precision, intent(out) :: d1(0:nz, 0:nz), d2(0:nz, 0:nz)
            integer                       :: i, j
            double precision              :: c(0:nz)
            double precision              :: rc(0:nz, 0:nz)
            double precision              :: fac, ss

            ! Compute Chebyshev points:
            fac = pi / dble(nz)   ! acos(-1) = pi
            do i = 0, nz
                x(i) = cos(fac * dble(i))
            enddo

            ! Compute weights:
            c(0) = two
            do i = 1, nz-1
                c(i) = one
            enddo
            c(nz) = two

            ss = one
            do i = 0, nz
                c(i) = c(i) * ss
                ss = -ss
            enddo

            ! Compute differentiation matrix:
            do i = 0, nz
                do j = 0, nz
                    rc(i, j) = c(i) / c(j)
                    d1(i, j) = x(i) - x(j)
                enddo
            enddo
            do i = 0, nz
                d1(i, i) = d1(i, i) + one
            enddo

            d1 = rc / d1
            do i = 0, nz
                d1(i, i) = one - sum(d1(i, :))
            enddo

            ! Set up second-order differentiation matrix:
            d2 = matmul(d1, d1)

        end subroutine init_cheby

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! w:  Clenshaw-Curtis weights
        subroutine clencurt(w)
            double precision, intent(out) :: w(0:nz)
            double precision              :: v(0:nz-2)
            double precision              :: theta
            integer                       :: i, k

            ! Initialize weights
            w = zero
            w(0) = one / (dble(nz) ** 2 - one)
            w(nz) = w(0)

            ! Initialize vector v
            v = one

            ! Update vector v using the given formula
            do k = 1, nz/2 - 1
                do i = 1, nz-1
                    theta = pi * dble(i) / dble(nz)
                    v(i-1) = v(i-1) - two * cos(two * dble(k) * theta) / (four * dble(k) ** 2 - one)
                end do
            end do

            ! Finalize vector v
            do i = 1, nz-1
                theta = pi * dble(i) / dble(nz)
                v(i-1) = v(i-1) - cos(dble(nz) * theta) / (dble(nz) ** 2 - one)
            end do

            ! Compute the weights
            w(1:nz-1) = two * v / dble(nz)

        end subroutine clencurt

end module cheby
