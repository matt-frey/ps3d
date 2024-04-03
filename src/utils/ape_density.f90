module ape_density
    use constants
    use physics, only : bfsq
    implicit none

    public :: ape_den

    contains

        elemental function ape_den(b, z) result(a)
            double precision, intent(in) :: b       ! perturbation buoyancy value
            double precision, intent(in) :: z       ! height
            double precision             :: a       ! APE density
            double precision             :: br

            br = b

#ifdef ENABLE_IW_TEST_CASE
            br = br + bfsq * z
            br = max(br, -twopi)
            br = min(br,  twopi)
            br = br - bfsq * z
#endif

#ifdef ENABLE_RT_TEST_CASE
            br = b + bfsq * z  ! convert buoyancy perturbation to total buoyancy
            br = max(br, -one)
            br = min(br,  one)
            a = br * dasin(br) + dsqrt(one - br ** 2) - z * br - dcos(z)
#else
            ! linear stratification:
            a = br ** 2 / (two * bfsq)
#endif
        end function ape_den

end module ape_density
