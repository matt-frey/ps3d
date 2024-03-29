module ape_density
    use constants
    implicit none

    public :: ape_den

    contains

        elemental function ape_den(b, z) result(a)
            double precision, intent(in) :: b       ! buoyancy value
            double precision, intent(in) :: z       ! height
            double precision             :: a       ! APE density
#if defined(ENABLE_IW_TEST_CASE) || defined(ENABLE_RT_TEST_CASE)
            double precision             :: br
#endif

#ifdef ENABLE_IW_TEST_CASE
            br = max(b, -twopi)
            br = min(br, twopi)
            a = f18 * (br - four * z) ** 2
#elif ENABLE_RT_TEST_CASE
            br = max(b, -one)
            br = min(br, one)
            a = br * dasin(br) + dsqrt(one - br ** 2) - z * br - dcos(z)
#else
            ! dummy line to avoid compiler warning of 'unused variables'
            a = b + z

            a = zero
#endif
        end function ape_den

end module ape_density
