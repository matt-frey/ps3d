program fftw_example
        !!!! To compile with fftw3 installed: gfortran -o fftw_example fftw_example.f90 -lfftw3 -lm
    use, intrinsic :: iso_c_binding
    implicit none
    include 'fftw3.f03'        ! Include the FFTW Fortran 90 interface

    !integer(C_INT), parameter :: FFTW_ESTIMATE = 64  ! FFTW_ESTIMATE flag
    integer(C_INT) :: N                             ! Length of the array
    real(C_DOUBLE), allocatable :: in(:)            ! Input array
    complex(C_DOUBLE_COMPLEX), allocatable :: out(:) ! Output array
    type(C_PTR) :: plan                             ! FFTW plan pointer
    integer :: i                                    ! Loop variable
    double precision, allocatable :: f(:)
    double precision, allocatable :: z(:)
    double precision, allocatable :: coeff(:)
    integer                       :: iz

    ! Define the size of the array
    N = 32                      ! Example array length, can be any integer

    ! Allocate input and output arrays
    allocate(z(N+1))
    allocate(f(N+1))
    allocate(in(2*N))
    allocate(out(N + 1))     ! For real-to-complex FFT, output is N/2+1 complex numbers

    call get_cheb_nodes(z,N)

    ! Initialize the input array with some known coefficients
    do i = 1, N+1
        f(i)  =       8.0d0 * get_cheb_poly(z(i), 1) &
                    - 9.0d0 * get_cheb_poly(z(i), 2) &
                    + 3.0d0 * get_cheb_poly(z(i), 3) &
                    - 4.0d0 * get_cheb_poly(z(i), 4) &
                    + 5.0d0 * get_cheb_poly(z(i), 5) &
                    + 2.5d0 * get_cheb_poly(z(i), 9)
    end do
   
    !  String these out into theta space:
                    in(1:N+1) =f
                    in(N+2:) = f(N:2:-1)


    ! Create the FFTW plan for real-to-complex transform
    plan = fftw_plan_dft_r2c_1d(2*N, in, out, FFTW_ESTIMATE)

    ! Execute the FFT
    call fftw_execute_dft_r2c(plan, in, out)

    out = out/n
    out(1) = 0.5d0*out(1)
    out(N+1) = 0.5d0*out(N+1)


    ! Print the results
    print *, "1st 10 Chebyshev Coefficients:"
    do i = 1, 10
         write(*, '(I5, F18.8)') i, REAL(out(i))
    end do

    ! Clean up
    call fftw_destroy_plan(plan)
    call fftw_cleanup()

    ! Deallocate arrays
    deallocate(in, out)

    contains

        function get_cheb_poly(z, n) result(res)
            double precision, intent(in) :: z
            integer,          intent(in) :: n
            double precision             :: res

            res = cos(dble(n) * acos(z))
        end function get_cheb_poly

        subroutine get_cheb_nodes(x,n) 
            integer,          intent(in) :: n
            double precision, allocatable,  intent(out) :: x(:)
            integer          :: i 
            double precision  :: pi
            allocate(x(n+1))
            pi = dacos(-1.0d0)

            do i =1,n+1
             x(i) = cos(pi*(i-1)/n);
            enddo

        end subroutine



end program fftw_example

