!-----------------------------------------------------------------------
! Module defining initial condition types.
!-----------------------------------------------------------------------
module cfdParams
  use realSizes, only: dp

  !Implicit declaration.
  implicit none

  !Define initial conditions:
  integer, parameter :: SOD_PROBLEM = 1
  integer, parameter :: MODIFIED_SOD = 2
  integer, parameter :: STRONG_SOD = 3
  integer, parameter :: PROBLEM_123 = 4
  integer, parameter :: THREE_RIGHT_WAVES = 5
  integer, parameter :: STATIONARY_CONTACT = 6
  integer, parameter :: SUBSONIC_NOZZLE = 7
  integer, parameter :: TRANSONIC_NOZZLE = 8
  integer, parameter :: SQUARE_WAVE = 9
  integer, parameter :: SINE_SQUARED_WAVE = 10

  !Define grid types:
  integer, parameter :: GRID_SHOCK_TUBE = 0
  integer, parameter :: GRID_NOZZLE = 1
  integer, parameter :: GRID_WAVE_PROBLEM = 2

  !Define grid stretching function types:
  integer, parameter :: STR_LINEAR = 0
  integer, parameter :: STR_MIN = 1
  integer, parameter :: STR_MAX = 2
  integer, parameter :: STR_MINMAX = 3
  integer, parameter :: STR_MIDPT = 4
  integer, parameter :: STR_SINE = 5
  integer, parameter :: STR_COSINE = 6

  !Define boundary conditions:
  integer, parameter :: BC_NONE = 0
  integer, parameter :: BC_FIXED = 1
  integer, parameter :: BC_CONSTANT_EXTRAPOLATION = 2
  integer, parameter :: BC_CHARACTERISTIC = 3
  integer, parameter :: BC_REFLECTION = 4
  integer, parameter :: BC_PERIODIC = 5

  !Define time-marching types:
  integer, parameter :: EXPLICIT_EULER = 0
  integer, parameter :: PREDICTOR_CORRECTOR = 1
  integer, parameter :: RUNGE_KUTTA = 2
  integer, parameter :: MULTISTAGE_OPTIMAL_SMOOTHING = 3

  !Define time-stepping types:
  integer, parameter :: GLOBAL_TIME_STEP = 0
  integer, parameter :: LOCAL_TIME_STEP = 1

  !Define flux functions:
  integer, parameter :: FLUX_GODUNOV = 0
  integer, parameter :: FLUX_ISENTROPIC = 1
  integer, parameter :: FLUX_RUSANOV = 2
  integer, parameter :: FLUX_HLLE = 3
  integer, parameter :: FLUX_HLLL = 4
  integer, parameter :: FLUX_HLLC = 5
  integer, parameter :: FLUX_OSHER = 6
  integer, parameter :: FLUX_ROE = 7
  integer, parameter :: FLUX_VANLEER = 8
  integer, parameter :: FLUX_AUSMplus = 9
  integer, parameter :: FLUX_AUSMplusup = 10
  integer, parameter :: FLUX_SLAU = 11

  !Define reconstruction types:
  integer, parameter :: RECONSTRUCTION_GG = 0
  integer, parameter :: RECONSTRUCTION_LSQ = 1
  integer, parameter :: RECONSTRUCTION_WENO = 2
  integer, parameter :: RECONSTRUCTION_PPM = 3
  integer, parameter :: RECONSTRUCTION_GAMMA = 4

  !Define limiter types:
  integer, parameter :: LIMITER_ZERO = 0
  integer, parameter :: LIMITER_ONE = 1
  integer, parameter :: LIMITER_MINMOD = 2
  integer, parameter :: LIMITER_UMIST = 3
  integer, parameter :: LIMITER_DOUBLE_MINMOD = 4
  integer, parameter :: LIMITER_SUPERBEE = 5
  integer, parameter :: LIMITER_PHI = 6
  integer, parameter :: LIMITER_VANLEER = 7
  integer, parameter :: LIMITER_VANALBADA = 8
  integer, parameter :: LIMITER_SINE = 9
  integer, parameter :: LIMITER_BARTH_JESPERSEN = 10
  integer, parameter :: LIMITER_VENKATAKRISHNAN = 11

contains

  !---------------------------------------------------------------------
  ! Min-mod function with two input variables.
  !---------------------------------------------------------------------
  real(dp) function minmod2(a, b) result(v)
    real(dp), intent(in) :: a, b
    v = sign(1.0_dp,a)*max(0.0_dp,min(abs(a),sign(1.0_dp,a)*b))
    return
  end function minmod2

  !---------------------------------------------------------------------
  ! Min-mod function with three input variables.
  !---------------------------------------------------------------------
  real(dp) function minmod3(a, b, c) result(v)
    real(dp), intent(in) :: a, b, c
    v = sign(1.0_dp,a)*max(0.0_dp,min(abs(a),min(sign(1.0_dp,a)*b,sign(1.0_dp,a)*c)))
    return
  end function minmod3

  !---------------------------------------------------------------------
  ! Double min-mod function.
  !---------------------------------------------------------------------
  real(dp) function dminmod(a, b, c, d) result(v)
    real(dp), intent(in) :: a, b, c, d
    v = sign(1.0_dp,a)*max(0.0_dp,min(abs(a),min(sign(1.0_dp,a)*b,min(sign(1.0_dp,a)*c,sign(1.0_dp,a)*d))))
    return
  end function dminmod

  !---------------------------------------------------------------------
  ! Superbee limiter function.
  !---------------------------------------------------------------------
  real(dp) function superbee(a, b) result(v)
    real(dp), intent(in) :: a, b
    v = sign(1.0_dp,a)*max(0.0_dp,max(min(2.0_dp*abs(a),sign(1.0_dp,a)*b),min(abs(a),2.0_dp*sign(1.0_dp,a)*b)))
    return
  end function superbee

  !---------------------------------------------------------------------
  ! Phi-limiter function.
  !---------------------------------------------------------------------
  real(dp) function philimiter(a, b, c) result(v)
    real(dp), intent(in) :: a, b, c
    v = sign(1.0_dp,a)*max(abs(minmod2(c*a,b)),abs(minmod2(a,c*b)))
    return
  end function philimiter

  !---------------------------------------------------------------------
  ! Van Leer limiter routine.
  !---------------------------------------------------------------------
  real(dp) function vanleer1(a, b) result(v)
    use numbers, only: TOLER
    real(dp), intent(in) :: a, b
    v = (abs(a*b)+a*b)/(a+b+sign(1.0_dp,a+b)*(TOLER**2))
    return
  end function vanleer1

  real(dp) function vanleer(ul, um, ur, umin, umax) result(v)
    use numbers, only: TOLER
    real(dp), intent(in) :: ul, um, ur, umin, umax
    real(dp)             :: y
    v = 1.0_dp
    if(ul-um.gt.TOLER) then
      y = (umax-um)/(ul-um)
      v = min(v,0.5_dp*vanleer1(1.0_dp,y))
    else if(ul-um.lt.-TOLER) then
      y = (umin-um)/(ul-um)
      v = min(v,0.5_dp*vanleer1(1.0_dp,y))
    else
      v = min(v,1.0_dp)
    end if
    if(ur-um.gt.TOLER) then
      y = (umax-um)/(ur-um)
      v = min(v,0.5_dp*vanleer1(1.0_dp,y))
    else if(ur-um.lt.-TOLER) then
      y = (umin-um)/(ur-um)
      v = min(v,0.5_dp*vanleer1(1.0_dp,y))
    else
      v = min(v,1.0_dp)
    end if
    return
  end function vanleer

  !---------------------------------------------------------------------
  ! van Albada limiter functions.
  !---------------------------------------------------------------------
  real(dp) function vanalbada1(a, b, eps) result(v)
    real(dp), intent(in) :: a, b, eps
    v = (a*b+(eps*eps))*(a+b)/((a*a)+(b*b)+2.*(eps*eps))
    return
  end function vanalbada1

  real(dp) function vanalbada(ul, um, ur, umin, umax) result(v)
    use numbers, only: TOLER, DECI
    real(dp), intent(in) :: ul, um, ur, umin, umax
    real(dp)             :: y
    v = 1.0_dp
    if(ul-um.gt.TOLER) then
      y = (umax-um)/(ul-um)
      v = min(v,0.5_dp*vanalbada1(1.0_dp,y,DECI))
    else if(ul-um.lt.-TOLER) then
       y = (umin-um)/(ul-um)
       v = min(v,0.5_dp*vanalbada1(1.0_dp,y,DECI))
    else
       v = min(v,1.0_dp)
    end if
    if(ur-um.gt.TOLER) then
       y = (umax-um)/(ur-um)
       v = min(v,0.5_dp*vanalbada1(1.0_dp,y,DECI))
    else if(ur-um.lt.-TOLER) then
       y = (umin-um)/(ur-um)
       v = min(v,0.5_dp*vanalbada1(1.0_dp,y,DECI))
    else
       v = min(v,1.0_dp)
    end if
    return
  end function vanalbada

  !---------------------------------------------------------------------
  ! Barth-Jespersen limiter function.
  !---------------------------------------------------------------------
  real(dp) function BarthJespersen(ul, um, ur, umin, umax) result(v)
    use numbers, only: TOLER
    real(dp), intent(in) :: ul, um, ur, umin, umax
    v = 1.0_dp
    if(ul-um.gt.TOLER) then
      v = min(v,min(1.0_dp,(umax-um)/(ul-um)))
    else if(ul-um.lt.-TOLER) then
      v = min(v,min(1.0_dp,(umin-um)/(ul-um)))
    else
      v = min(v,1.0_dp)
    end if
    if(ur-um.gt.TOLER) then
      v = min(v,min(1.0_dp,(umax-um)/(ur-um)))
    else if(ur-um.lt.-TOLER) then
      v = min(v,min(1.0_dp,(umin-um)/(ur-um)))
    else
      v = min(v,1.0_dp)
    end if
    return
  end function BarthJespersen

  !---------------------------------------------------------------------
  ! Venkatakrishnan limiter function.
  !---------------------------------------------------------------------
  real(dp) function Venkatakrishnan(ul, um, ur, umin, umax) result(v)
    use numbers, only: TOLER
    real(dp), intent(in) :: ul, um, ur, umin, umax
    real(dp)             :: y
    v = 1.0_dp
    if(ul-um.gt.TOLER) then
      y = (umax-um)/(ul-um)
      v = min(v,(y*y+2.0_dp*y)/(y*y+y+2.0_dp))
    else if(ul-um.lt.-TOLER) then
      y = (umin-um)/(ul-um)
      v = min(v,(y*y+2.0_dp*y)/(y*y+y+2.0_dp))
    else
      v = min(v,1.0_dp)
    end if
    if(ur-um.gt.TOLER) then
      y = (umax-um)/(ur-um)
      v = min(v,(y*y+2.0_dp*y)/(y*y+y+2.0_dp))
    else if(ur-um.lt.-TOLER) then
      y = (umin-um)/(ur-um)
      v = min(v,(y*y+2.0_dp*y)/(y*y+y+2.0_dp))
    else
      v = min(v,1.0_dp)
    end if
    return
  end function Venkatakrishnan

end module cfdParams
