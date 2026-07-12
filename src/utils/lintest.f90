program lintest
    use linsolve
    implicit none

    double precision :: A(3, 3)
    double precision :: b(3), x(3)

    A(1, 1) = 1.0d0
    A(1, 2) = 2.0d0
    A(1, 3) = 0.0d0

    A(2, 1) = 4.0d0
    A(2, 2) = 5.0d0
    A(2, 3) = 6.0d0

    A(3, 1) = 2.0d0
    A(3, 2) = 2.0d0
    A(3, 3) = 9.0d0

    b(1) = -1.0d0
    b(2) = 2.0d0
    b(3) = -3.0d0

    call solve3x3(A, b, x)

    print *, x

end program
