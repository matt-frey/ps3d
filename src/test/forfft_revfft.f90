program fft_test
        !!!! To compile: gfortran constants.f90 stafft.f90 forfft_revfft.f90
    use constants, only : pi
    use stafft, only : initfft, forfft, revfft, dct
    implicit none

    integer :: N                                        ! Length of the array
    integer :: i                                        ! Loop variable
    double precision, allocatable :: fcos(:), fsin(:)
    double precision, allocatable :: g(:)
    double precision, allocatable :: trig(:)
    integer                       :: factors(5)
    double precision :: x

    ! Define the size of the array
    N = 16

    ! Allocate input and output arrays
    allocate(fcos(0:N-1))
    allocate(fsin(0:N-1))
    allocate(g(0:N-1))
    allocate(trig(2*N))

    call initfft(N, factors, trig)

    ! Initialize the input array with some known coefficients
    do i = 0, N-1
        x = 2.0d0 * pi / dble(N) * dble(i)
        fcos(i) = 4.0d0 + cos(3.0d0 * x)
        g(i) = 4.0d0 + cos(3.0d0 * x)
        fsin(i) = 4.0d0 + sin(3.0d0 * x)
    end do

    call forfft(1, N, fcos, trig, factors)
    call dct(1, N, g, trig, factors)
    call forfft(1, N, fsin, trig, factors)

    print *, "index, FFT(4 + cos(3*x)), DCT(4 + cos(3*x)), FFT(4 + sin(3*x))"
    do i = 0, N-1
        print "(I2, F18.8, F18.8, F18.8)", i, fcos(i), g(i), fsin(i)
    enddo

    ! Deallocate arrays
    deallocate(fcos, fsin, trig)

end program fft_test

