module linsolve
    use constants, only : one
    implicit none

    contains

        ! Solves a 3x3 linear system Ax = b using Cramer's rule
        ! (see https://en.wikipedia.org/wiki/Cramer%27s_rule#)
        subroutine solve3x3(A, b, x)
            double precision, intent(in)  :: A(3, 3)
            double precision, intent(in)  :: b(3)
            double precision, intent(out) :: x(3)
            double precision              :: C(3, 3)
            double precision              :: detA, det1, det2, det3
            double precision              :: a1, a2, a3


            a1 = A(1, 2) * A(2, 3) - A(1, 3) * A(2, 2)
            a2 = A(1, 3) * A(2, 1) - A(1, 1) * A(2, 3)
            a3 = A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)

            detA = A(3, 1) * a1 + A(3, 2) * a2 + A(3, 3) * a3

            det1 = b(3)    * (a1) &
                 + A(3, 2) * (A(1, 3) * b(2)    - b(1)    * A(2, 3)) &
                 + A(3, 3) * (b(1)    * A(2, 2) - A(1, 2) * b(2)   )

            det2 = A(3, 1) * (b(1)    * A(2, 3) - A(1, 3) * b(2))    &
                 + b(3)    * a2                                      &
                 + A(3, 3) * (A(1, 1) * b(2)    - b(1)    * A(2, 1))

            det3 = A(3, 1) * (A(1, 2) * b(2)    - b(1)    * A(2, 2)) &
                 + A(3, 2) * (b(1)    * A(2, 1) - A(1, 1) * b(2)   ) &
                 + b(3)    * a3

            detA = one / detA
            x(1) = det1 * detA
            x(2) = det2 * detA
            x(3) = det3 * detA

        end subroutine solve3x3

end module linsolve
