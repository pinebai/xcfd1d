!-----------------------------------------------------------------------
! Module for declaring single and double precision values.
!-----------------------------------------------------------------------
module realSizes
  !Single precision definition.
  integer, parameter :: sp = kind(0.e0)!selectED_REAL_KIND(p=6,r=37)
  !Double precision definition.
  integer, parameter :: dp = kind(0.d0)!selectED_REAL_KIND(p=13,r=200)
end module realSizes


!-----------------------------------------------------------------------
! Module for declaring useful numbers.
!-----------------------------------------------------------------------
module numbers
  !Include real sizes module.
  use realSizes, only: dp
  !Set prefixes:
  real(dp), parameter :: MEGA  = 1.d6
  real(dp), parameter :: KILO  = 1.d3
  real(dp), parameter :: HECTO = 1.d2
  real(dp), parameter :: DEKA  = 1.d1
  real(dp), parameter :: DECI  = 1.d-1
  real(dp), parameter :: CENTI = 1.d-2
  real(dp), parameter :: MILLI = 1.d-3
  real(dp), parameter :: MICRO = 1.d-6
  real(dp), parameter :: TOLER = 1.d-8 ! Not an SI prefix.
  real(dp), parameter :: NANO  = 1.d-9
  real(dp), parameter :: PICO  = 1.d-12
  real(dp), parameter :: FEMTO = 1.d-15
  !Set pi:
  real(dp), parameter :: PI = 3.14159265358979323846_dp
end module numbers


!-----------------------------------------------------------------------
! Module for useful mathematical functions.
!-----------------------------------------------------------------------
module mathFunctions
  use realSizes, only: dp

contains

  real(dp) function sqr(x)
    implicit none
    real(dp), intent(in) :: x
    sqr = x*x
    return
  end function sqr

  real(dp) function cube(x)
    implicit none
    real(dp), intent(in) :: x
    cube = x*x*x
    return
  end function cube

end module mathFunctions
