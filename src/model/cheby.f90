module cheby
    use constants, only : zero, one, two, pi, four, f12
    use stafft, only : initfft, forfft, revfft
    implicit none

    private

        double precision, allocatable :: trig(:)
        integer                       :: factors(5)

    public :: init_cheby        &
            , finalise_cheby    &
            , clencurt          &
            , cheb_poly         &
            , cheb_fun

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

        allocate(trig(4*n))
        call initfft(2*n, factors, trig)

    end subroutine init_cheby

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine finalise_cheby
        if (allocated(trig)) then
            deallocate(trig)
        endif
    end subroutine finalise_cheby

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

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

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Input:
    ! f - a vector of length n+1 containing function values at Chebyshev nodes in [-1, 1]
    ! Output:
    ! c - a vector of length n+1 containing the coefficients of the Chebyshev polynomials
    subroutine cheb_poly(n, f, c)
        integer,          intent(in)  :: n
        double precision, intent(in)  :: f(0:n)
        double precision, intent(out) :: c(0:n)
        double precision              :: valsUnitDisc(0:2*n-1)

        ! Fill valsUnitDisc with f and mirrored values:
        valsUnitDisc(0:n) = f(0:n)
        valsUnitDisc(n+1:) = f(n-1:1:-1)

        ! Perform the FFT:
        call forfft(1, 2*n, valsUnitDisc, trig, factors)
        c(0:n) = valsUnitDisc(0:n)

        ! Get Chebyshev coefficients:
        ! Note: forfft normalises with 1 / sqrt(N) where N = 2*n;
        !       we must therefore multiply the output with sqrt(2 / n)
        !       in order to get a normalisation of 1 / n
        c = c * sqrt(two / dble(n))
        c(0) = f12 * c(0)
        c(n) = f12 * c(n)

    end subroutine cheb_poly

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Given n+1 coefficients of the Chebyshev series
    ! Returns f(x) evaluated at the n+1 chebyshev nodes in [-1,1];
    subroutine cheb_fun(n, c, f)
        integer,          intent(in) :: n
        double precision, intent(in)  :: c(0:n)
        double precision, intent(out) :: f(0:n)
        double precision              :: coeffs(0:2*n-1)

        f = c
        f(0) = two * f(0)
        f(n) = two * f(n)
        ! We need to undo the  1 / n scaling, but then get the 1 / sqrt(2 * n)
        ! scaling of forfft --> multiply by sqrt(n / 2)
        f = f * sqrt(dble(n) / two)

        ! Want 2N coefficients: First N+1 are the cosines we have:
        coeffs = zero
        coeffs(0:n) = f(0:n)

        ! Perform the inverse FFT:
        call revfft(1, 2*n, coeffs, trig, factors)

        ! Get Chebyshev coefficients:
        f(0:n) = coeffs(0:n)

    end subroutine cheb_fun

end module cheby
