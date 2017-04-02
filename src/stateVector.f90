!-----------------------------------------------------------------------
! Module for a (rho, du, E) conserved-variable solution state
!-----------------------------------------------------------------------
module Euler1D_UState
  use realSizes, only: dp

  implicit none

  !State variables
  type Euler1D_U_State
     real(dp) :: rho   !Gas density
     real(dp) :: du    !Gas momentum
     real(dp) :: E     !Gas total energy
  end type Euler1D_U_State

  !Assignment operator
  interface assignment(=)
    module procedure copy_U
  end interface

  !Addition operator
  interface operator(+)
    module procedure add_U
  end interface

  !Subtraction operator
  interface operator(-)
    module procedure subtract_U
  end interface

  !Multiplication operator
  interface operator(*)
    module procedure scalar_multiply1_U
    module procedure scalar_multiply2_U
    module procedure outer_product_U
  end interface

  !Multiplication operator
  interface operator(.dot.)
    module procedure inner_product_U
  end interface

  !Division operator
  interface operator(/)
    module procedure scalar_division_U
  end interface

contains

  type(Euler1D_U_State) function add_U(U1, U2)
    type(Euler1D_U_State), intent(in) :: U1, U2
    add_U%rho = U1%rho + U2%rho
    add_U%du = U1%du + U2%du
    add_U%E = U1%E + U2%E
    return
  end function add_U

  type(Euler1D_U_State) function subtract_U(U1, U2)
    type(Euler1D_U_State), intent(in) :: U1, U2
    subtract_U%rho = U1%rho - U2%rho
    subtract_U%du = U1%du - U2%du
    subtract_U%E = U1%E - U2%E
    return
  end function subtract_U

  type(Euler1D_U_State) function scalar_multiply1_U(a, U)
    real(dp),              intent(in) :: a
    type(Euler1D_U_State), intent(in) :: U
    scalar_multiply1_U%rho = a*U%rho
    scalar_multiply1_U%du = a*U%du
    scalar_multiply1_U%E = a*U%E
    return
  end function scalar_multiply1_U

  type(Euler1D_U_State) function scalar_multiply2_U(U, a)
    type(Euler1D_U_State), intent(in) :: U
    real(dp),              intent(in) :: a
    scalar_multiply2_U%rho = a*U%rho
    scalar_multiply2_U%du = a*U%du
    scalar_multiply2_U%E = a*U%E
    return
  end function scalar_multiply2_U

  type(Euler1D_U_State) function outer_product_U(U1, U2)
    type(Euler1D_U_State), intent(in) :: U1, U2
    outer_product_U%rho = U1%rho*U2%rho
    outer_product_U%du = U1%du*U2%du
    outer_product_U%E = U1%E*U2%E
    return
  end function outer_product_U

  real(dp) function inner_product_U(U1, U2)
    type(Euler1D_U_State), intent(in) :: U1, U2
    inner_product_U = U1%rho*U2%rho + U1%du*U2%du + U1%E*U2%E
    return
  end function inner_product_U

  type(Euler1D_U_State) function scalar_division_U(U, a)
    type(Euler1D_U_State), intent(in) :: U
    real(dp),              intent(in) :: a
    scalar_division_U%rho = U%rho/a
    scalar_division_U%du = U%du/a
    scalar_division_U%E = U%E/a
    return
  end function scalar_division_U

  !Set to vacuum state
  subroutine vacuum_U(U)
    type(Euler1D_U_State), intent(out) :: U
    U%rho = 0.0_dp
    U%du = 0.0_dp
    U%E = 0.0_dp
    return
  end subroutine vacuum_U

  !Set state to standard atmosphere
  subroutine standard_atmosphere_U(U)
    use gasConstants, only: gm1i
    type(Euler1D_U_State), intent(out) :: U
    U%rho = 1.225_dp
    U%du = 0.0_dp
    U%E = gm1i*101325.0_dp
    return
  end subroutine standard_atmosphere_U

  !Set state to specified constant
  subroutine constant_U(U, val)
    type(Euler1D_U_State), intent(out) :: U
    real(dp),              intent(in)  :: val
    U%rho = val
    U%du = val
    U%E = val
    return
  end subroutine constant_U

  !Copy the specified state
  subroutine copy_U(U1, U2)
    type(Euler1D_U_State), intent(out) :: U1
    type(Euler1D_U_State), intent(in)  :: U2
    U1%rho = U2%rho
    U1%du = U2%du
    U1%E = U2%E
    return
  end subroutine copy_U

  !Set the state to the specified values
  type(Euler1D_U_State) function UState(rho, du, E)
    real(dp), intent(in) :: rho, du, E
    UState%rho = rho
    UState%du = du
    UState%E = E
    return
  end function UState

  !Determine the minimum state values
  type(Euler1D_U_State) function min_U(U1, U2)
    type(Euler1D_U_State), intent(in) :: U1
    type(Euler1D_U_State), intent(in) :: U2
    min_U%rho = min(U1%rho,U2%rho)
    min_U%du = min(U1%du,U2%du)
    min_U%E = min(U1%E,U2%E)
    return
  end function min_U

  !Determine the maximum state values
  type(Euler1D_U_State) function max_U(U1, U2)
    type(Euler1D_U_State), intent(in) :: U1
    type(Euler1D_U_State), intent(in) :: U2
    max_U%rho = max(U1%rho,U2%rho)
    max_U%du = max(U1%du,U2%du)
    max_U%E = max(U1%E,U2%E)
    return
  end function max_U

  !Determine the absolute values of the state variables
  type(Euler1D_U_State) function abs_U(U)
    type(Euler1D_U_State), intent(in) :: U
    abs_U%rho = abs(U%rho)
    abs_U%du = abs(U%du)
    abs_U%E = abs(U%E)
    return
  end function abs_U

  !Determine the square of the state variables
  type(Euler1D_U_State) function sqr_U(U)
    type(Euler1D_U_State), intent(in) :: U
    sqr_U%rho = U%rho*U%rho
    sqr_U%du = U%du*U%du
    sqr_U%E = U%E*U%E
    return
  end function sqr_U

  !Determine the square-root of the state variables
  type(Euler1D_U_State) function sqrt_U(U)
    type(Euler1D_U_State), intent(in) :: U
    sqrt_U%rho = sqrt(U%rho)
    sqrt_U%du = sqrt(U%du)
    sqrt_U%E = sqrt(U%E)
    return
  end function sqrt_U

  !Calculate the flow speed
  real(dp) function u_U(U)
    type(Euler1D_U_State), intent(in) :: U
    u_U = U%du/U%rho
    return
  end function u_U

  !Calculate the pressure
  real(dp) function p_U(U)
    use gasConstants, only: gm1
    type(Euler1D_U_State), intent(in) :: U
    p_U = gm1*(U%E - 0.5_dp*U%du*U%du/U%rho)
    return
  end function p_U

  !Calculate the temperature of the gas
  real(dp) function T_U(U)
    use gasConstants, only: R
    type(Euler1D_U_State), intent(in) :: U
    T_U = p_U(U)/(U%rho*R)
    return
  end function T_U

  !Calculate the specific internal energy of the gas
  real(dp) function eps_U(U)
    use gasConstants, only: gm1i
    type(Euler1D_U_State), intent(in) :: U
    eps_U = gm1i*p_U(U)/U%rho
    return
  end function eps_U

  !Calculate the specific internal energy of the gas
  real(dp) function Esp_U(U)
    type(Euler1D_U_State), intent(in) :: U
    Esp_U = U%E/U%rho
    return
  end function Esp_U

  !Calculate the specific internal enthalpy of the gas
  real(dp) function hint_U(U)
    use gasConstants, only: g, gm1i
    type(Euler1D_U_State), intent(in) :: U
    hint_U = g*gm1i*p_U(U)/U%rho
    return
  end function hint_U

  !Calculate the total enthalpy of the gas
  real(dp) function H_U(U)
    type(Euler1D_U_State), intent(in) :: U
    H_U = p_U(U) + U%E
    return
  end function H_U

  !Calculate the specific enthalpy of the gas
  real(dp) function Hsp_U(U)
    type(Euler1D_U_State), intent(in) :: U
    Hsp_U = (U%E + p_U(U))/U%rho
    return
  end function Hsp_U

  !Calculate the speed of sound squared
  real(dp) function a2_U(U)
    use gasConstants, only: g
    type(Euler1D_U_State), intent(in) :: U
    a2_U = g*p_U(U)/U%rho
    return
  end function a2_U

  !Calculate the speed of sound
  real(dp) function a_U(U)
    type(Euler1D_U_State), intent(in) :: U
    a_U = sqrt(a2_U(U))
    return
  end function a_U

  !Calculate the flow Mach number
  real(dp) function M_U(U)
    type(Euler1D_U_State), intent(in) :: U
    M_U = u_U(U)/a_U(U)
    return
  end function M_U

  !Calculate the entropy of the gas
  real(dp) function s_U(U)
    use gasConstants, only: g, gm1i, R
    type(Euler1D_U_State), intent(in) :: U
    s_U = R*gm1i*log(p_U(U)/(U%rho**g))
    return
  end function s_U

  !Check for unphysical properties
  integer function unphysical_properties_U(U) result(ierr)
    type(Euler1D_U_State), intent(in) :: U
    if((U%rho.le.0.0_dp).or.(U%rho.ne.U%rho)) then
     ierr = 1
    else if(p_U(U).lt.0.0_dp) then
     ierr = 1
    else if((U%E.lt.0.0_dp).or.(U%E.ne.U%E)) then
     ierr = 1
    else
     ierr = 0
    end if
  end function unphysical_properties_U

  !Conserved right eigenvector
  subroutine rc_U(U, r, index)
    type(Euler1D_U_State), intent(in)  :: U
    type(Euler1D_U_State), intent(out) :: r
    integer,               intent(in)  :: index
    real(dp)  :: v, a, h
    v = u_U(U)
    a = a_U(U)
    h = Hsp_U(U)
    select case(index)
    case(1)
     r%rho = 1.0_dp
     r%du = v - a
     r%E = h - a*v
    case(2)
     r%rho = 1.0_dp
     r%du = v
     r%E = 0.5_dp*v*v
    case(3)
     r%rho = 1.0_dp
     r%du = v + a
     r%E = h + a*v
    end select
  end subroutine rc_U

  !Area-change source term
  type(Euler1D_U_State) function Sa(U, xA, dA)
    type(Euler1D_U_State), intent(in) :: U
    real(dp),              intent(in) :: xA, dA
    real(dp) :: H
    H = H_U(U)
    Sa%rho = -(U%du)*dA/xA
    Sa%du = -(U%du*U%du/U%rho)*dA/xA
    Sa%E = -(U%du*H/U%rho)*dA/xA
  end function Sa

  !Flux
  subroutine Flux_U(U, F)
    type(Euler1D_U_State), intent(in)  :: U
    type(Euler1D_U_State), intent(out) :: F
    real(dp) :: p
    p = p_U(U)
    F%rho = U%du
    F%du = U%du*U%du/U%rho + p
    F%E = U%du*(U%E+p)/U%rho
  end subroutine Flux_U

  !Flux Jacobian
  subroutine dFdU_U(U, dFdU)
    use gasConstants, only: g, gm1
    type(Euler1D_U_State),    intent(in)  :: U
    real(dp), dimension(3,3), intent(out) :: dFdU
    real(dp) :: v
    v = u_U(U)
    dFdU(1,1) = 0.0_dp
    dFdU(1,2) = 1.0_dp
    dFdU(1,3) = 0.0_dp
    dFdU(2,1) = 0.5_dp*(g-3.0_dp)*v*v
    dFdU(2,2) = (3.0_dp-g)*v
    dFdU(2,3) = gm1
    dFdU(3,1) = gm1*v*v*v - g*v*U%E/U%rho
    dFdU(3,2) = g*U%E/U%rho - 0.5_dp*3.0_dp*gm1*v*v
    dFdU(3,3) = g*v
  end subroutine dFdU_U

end module Euler1D_UState


!-----------------------------------------------------------------------
! Module for a (rho, u, p) primitive-variable solution state
!-----------------------------------------------------------------------
module Euler1D_WState
  use realSizes, only: dp
  use Euler1D_UState

  implicit none

  !State variables
  type Euler1D_W_State
     real(dp) :: rho   !Gas density
     real(dp) :: u     !Gas speed
     real(dp) :: p     !Gas pressure
  end type Euler1D_W_State

  !Assignment operator
  interface assignment(=)
    module procedure copy_W
  end interface

  !Addition operator
  interface operator(+)
    module procedure add_W
  end interface

  !Subtraction operator
  interface operator(-)
    module procedure subtract_W
  end interface

  !Multiplication operator
  interface operator(*)
    module procedure scalar_multiply1_W
    module procedure scalar_multiply2_W
    module procedure outer_product_W
  end interface

  !Multiplication operator
  interface operator(.dot.)
    module procedure inner_product_W
  end interface

  !Division operator
  interface operator(/)
    module procedure scalar_division_W
  end interface

contains

  type(Euler1D_W_State) function add_W(W1, W2)
    type(Euler1D_W_State), intent(in) :: W1, W2
    add_W%rho = W1%rho + W2%rho
    add_W%u = W1%u + W2%u
    add_W%p = W1%p + W2%p
    return
  end function add_W

  type(Euler1D_W_State) function subtract_W(W1, W2)
    type(Euler1D_W_State), intent(in) :: W1, W2
    subtract_W%rho = W1%rho - W2%rho
    subtract_W%u = W1%u - W2%u
    subtract_W%p = W1%p - W2%p
    return
  end function subtract_W

  type(Euler1D_W_State) function scalar_multiply1_W(a, W)
    real(dp),              intent(in) :: a
    type(Euler1D_W_State), intent(in) :: W
    scalar_multiply1_W%rho = a*W%rho
    scalar_multiply1_W%u = a*W%u
    scalar_multiply1_W%p = a*W%p
    return
  end function scalar_multiply1_W

  type(Euler1D_W_State) function scalar_multiply2_W(W, a)
    type(Euler1D_W_State), intent(in) :: W
    real(dp),              intent(in) :: a
    scalar_multiply2_W%rho = a*W%rho
    scalar_multiply2_W%u = a*W%u
    scalar_multiply2_W%p = a*W%p
    return
  end function scalar_multiply2_W

  type(Euler1D_W_State) function outer_product_W(W1, W2)
    type(Euler1D_W_State), intent(in) :: W1, W2
    outer_product_W%rho = W1%rho*W2%rho
    outer_product_W%u = W1%u*W2%u
    outer_product_W%p = W1%p*W2%p
    return
  end function outer_product_W

  real(dp) function inner_product_W(W1, W2)
    type(Euler1D_W_State), intent(in) :: W1, W2
    inner_product_W = W1%rho*W2%rho + W1%u*W2%u + W1%p*W2%p
    return
  end function inner_product_W

  type(Euler1D_W_State) function scalar_division_W(W, a)
    type(Euler1D_W_State), intent(in) :: W
    real(dp),              intent(in) :: a
    scalar_division_W%rho = W%rho/a
    scalar_division_W%u = W%u/a
    scalar_division_W%p = W%p/a
    return
  end function scalar_division_W

  !Set to vacuum state
  subroutine vacuum_W(W)
    type(Euler1D_W_State), intent(out) :: W
    W%rho = 0.0_dp
    W%u = 0.0_dp
    W%p = 0.0_dp
    return
  end subroutine vacuum_W

  !Set state to standard atmosphere
  subroutine standard_atmosphere_W(W)
    type(Euler1D_W_State), intent(out) :: W
    W%rho = 1.225_dp
    W%u = 0.0_dp
    W%p = 101325.0_dp
    return
  end subroutine standard_atmosphere_W

  !Set state to specified constant
  subroutine constant_W(W, val)
    type(Euler1D_W_State), intent(out) :: W
    real(dp),              intent(in)  :: val
    W%rho = val
    W%u = val
    W%p = val
    return
  end subroutine constant_W

  !Copy the specified state
  subroutine copy_W(W1, W2)
    type(Euler1D_W_State), intent(out) :: W1
    type(Euler1D_W_State), intent(in)  :: W2
    W1%rho = W2%rho
    W1%u = W2%u
    W1%p = W2%p
    return
  end subroutine copy_W

  !Set the state to the specified values
  type(Euler1D_W_State) function WState(rho, u, p)
    real(dp), intent(in) :: rho, u, p
    WState%rho = rho
    WState%u = u
    WState%p = p
    return
  end function WState

  !Determine the minimum state values
  type(Euler1D_W_State) function min_W(W1, W2)
    type(Euler1D_W_State), intent(in) :: W1
    type(Euler1D_W_State), intent(in) :: W2
    min_W%rho = min(W1%rho,W2%rho)
    min_W%u = min(W1%u,W2%u)
    min_W%p = min(W1%p,W2%p)
    return
  end function min_W

  !Determine the maximum state values
  type(Euler1D_W_State) function max_W(W1, W2)
    type(Euler1D_W_State), intent(in) :: W1
    type(Euler1D_W_State), intent(in) :: W2
    max_W%rho = max(W1%rho,W2%rho)
    max_W%u = max(W1%u,W2%u)
    max_W%p = max(W1%p,W2%p)
    return
  end function max_W

  !Calculate the momentum of the gas
  real(dp) function du_W(W)
    type(Euler1D_W_State), intent(in) :: W
    du_W = W%rho*W%u
    return
  end function du_W

  !Calculate the temperature of the gas
  real(dp) function T_W(W)
    use gasConstants, only: R
    type(Euler1D_W_State), intent(in) :: W
    T_W = W%p/(W%rho*R)
    return
  end function T_W

  !Calculate the specific internal energy of the gas
  real(dp) function eps_W(W)
    use gasConstants, only: gm1i
    type(Euler1D_W_State), intent(in) :: W
    eps_W = gm1i*W%p/W%rho
    return
  end function eps_W

  !Calculate the total energy of the gas
  real(dp) function E_W(W)
    use gasConstants, only: gm1i
    type(Euler1D_W_State), intent(in) :: W
    E_W = gm1i*W%p + 0.5_dp*W%rho*W%u*W%u
    return
  end function E_W

  !Calculate the specific total energy of the gas
  real(dp) function Esp_W(W)
    use gasConstants, only: gm1i
    type(Euler1D_W_State), intent(in) :: W
    Esp_W = gm1i*W%p/W%rho + 0.5_dp*W%u*W%u
    return
  end function Esp_W

  !Calculate the specific internal enthalpy of the gas
  real(dp) function hint_W(W)
    use gasConstants, only: g, gm1i
    type(Euler1D_W_State), intent(in) :: W
    hint_W = g*gm1i*W%p/W%rho
    return
  end function hint_W

  !Calculate the total enthalpy of the gas
  real(dp) function H_W(W)
    use gasConstants, only: g, gm1i
    type(Euler1D_W_State), intent(in) :: W
    H_W = g*gm1i*W%p + 0.5_dp*W%rho*W%u*W%u
    return
  end function H_W

  !Calculate the specific enthalpy of the gas
  real(dp) function Hsp_W(W)
    use gasConstants, only: g, gm1i
    type(Euler1D_W_State), intent(in) :: W
    Hsp_W = g*gm1i*W%p/W%rho + 0.5_dp*W%u*W%u
    return
  end function Hsp_W

  !Calculate the speed of sound squared
  real(dp) function a2_W(W)
    use gasConstants, only: g
    type(Euler1D_W_State), intent(in) :: W
    a2_W = g*W%p/W%rho
    return
  end function a2_W

  !Calculate the speed of sound
  real(dp) function a_W(W)
    type(Euler1D_W_State), intent(in) :: W
    a_W = sqrt(a2_W(W))
    return
  end function a_W

  !Calculate the flow Mach number
  real(dp) function M_W(W)
    type(Euler1D_W_State), intent(in) :: W
    M_W = abs(W%u)/a_W(W)
    return
  end function M_W

  !Calculate the entropy of the gas
  real(dp) function s_W(W)
    use gasConstants, only: g, gm1i, R
    type(Euler1D_W_State), intent(in) :: W
    s_W = R*gm1i*log(W%p/(W%rho**g))
    return
  end function s_W

  !Check for unphysical properties
  integer function unphysical_properties_W(W) result(ierr)
    type(Euler1D_W_State), intent(in) :: W
    if((W%rho.lt.0.0_dp).or.(W%rho.ne.W%rho)) then
     ierr = 1
    else if((W%p.lt.0.0_dp).or.(W%p.ne.W%p)) then
     ierr = 1
    else if(E_W(W).lt.0.0_dp) then
     ierr = 1
    else
     ierr = 0
    end if
  end function unphysical_properties_W

  !Flux Vector
  subroutine Flux_W(W, F)
    type(Euler1D_W_State), intent(in)  :: W
    type(Euler1D_U_State), intent(out) :: F
    F%rho = W%rho*W%u
    F%du = W%rho*W%u*W%u + W%p
    F%E = W%u*H_W(W)
  end subroutine Flux_W

  !Flux Jacobian
  subroutine dFdU_W(W, dFdU)
    use gasConstants, only: g, gm1
    type(Euler1D_W_State),    intent(in)  :: W
    real(dp), dimension(3,3), intent(out) :: dFdU
    real(dp) :: E
    E = E_W(W)
    dFdU(1,1) = 0.0_dp
    dFdU(1,2) = 1.0_dp
    dFdU(1,3) = 0.0_dp
    dFdU(2,1) = 0.5_dp*(g-3.0_dp)*W%u*W%u
    dFdU(2,2) = (3.0_dp-g)*W%u
    dFdU(2,3) = gm1
    dFdU(3,1) = gm1*W%u*W%u*W%u - g*W%u*E/W%rho
    dFdU(3,2) = g*E/W%rho - 0.5_dp*3.0_dp*gm1*W%u*W%u
    dFdU(3,3) = g*W%u
  end subroutine dFdU_W

  !Eigenvalue(s)
  subroutine lambda_W(W, lambda)
    type(Euler1D_W_State), intent(in)  :: W
    type(Euler1D_W_State), intent(out) :: lambda
    real(dp) :: a
    a = a_W(W)
    lambda%rho = W%u-a
    lambda%u = W%u
    lambda%p = W%u+a
  end subroutine lambda_W

  !Primitive right eigenvector
  subroutine rp_W(W, rp, index)
    type(Euler1D_W_State), intent(in)  :: W
    type(Euler1D_W_State), intent(out) :: rp
    integer,               intent(in)  :: index
    real(dp) :: a, a2
    a2 = a2_W(W)
    a = sqrt(a2)
    select case(index)
    case(1)
     rp%rho = 1.0_dp
     rp%u = -a/W%rho
     rp%p = a2
    case(2)
     rp%rho = 1.0_dp
     rp%u = 0.0_dp
     rp%p = 0.0_dp
    case(3)
     rp%rho = 1.0_dp
     rp%u = a/W%rho
     rp%p = a2
    end select
  end subroutine rp_W

  !Primitive left eigenvector
  subroutine lp_W(W, lp, index)
    type(Euler1D_W_State), intent(in)  :: W
    type(Euler1D_W_State), intent(out) :: lp
    integer,               intent(in)  :: index
    real(dp) :: a, a2
    a2 = a2_W(W)
    a = sqrt(a2)
    select case(index)
    case(1)
     lp%rho = 0.0_dp
     lp%u = -0.5_dp*W%rho/a
     lp%p = 0.5_dp/a2
    case(2)
     lp%rho = 1.0_dp
     lp%u = 0.0_dp
     lp%p = -1.0_dp/a2
    case(3)
     lp%rho = 0.0_dp
     lp%u = 0.5_dp*W%rho/a
     lp%p = 0.5_dp/a2
    end select
  end subroutine lp_W

  !Conserved right eigenvector
  subroutine rc_W(W, rc, index)
    type(Euler1D_W_State), intent(in)  :: W
    type(Euler1D_U_State), intent(out) :: rc
    integer,               intent(in)  :: index
    real(dp) :: a, h
    a = a_W(W)
    h = Hsp_W(W)
    select case(index)
    case(1)
     rc%rho = 1.0_dp
     rc%du = W%u - a
     rc%E = h - a*W%u
    case(2)
     rc%rho = 1.0_dp
     rc%du = W%u
     rc%E = 0.5_dp*W%u*W%u
    case(3)
     rc%rho = 1.0_dp
     rc%du = W%u + a
     rc%E = h + a*W%u
    end select
  end subroutine rc_W

  !Reflection boundary condition
  subroutine Reflect_W(W1, W2)
    type(Euler1D_W_State), intent(out) :: W1
    type(Euler1D_W_State), intent(in)  :: W2
    W1%rho = W2%rho
    W1%u   = - W2%u
    W1%p   = W2%p
  end subroutine Reflect_W

end module Euler1D_WState


!-----------------------------------------------------------------------
! Module for a (p, u, T) primitive-variable solution state
!-----------------------------------------------------------------------
module Euler1D_QState
  use realSizes, only: dp
  use Euler1D_UState

  implicit none

  !State variables
  type Euler1D_Q_State
     real(dp) :: p     !Gas pressure
     real(dp) :: u     !Gas speed
     real(dp) :: T     !Gas temperature
  end type Euler1D_Q_State

  !Assignment operator
  interface assignment(=)
    module procedure copy_Q
  end interface

  !Addition operator
  interface operator(+)
    module procedure add_Q
  end interface

  !Subtraction operator
  interface operator(-)
    module procedure subtract_Q
  end interface

  !Multiplication operator
  interface operator(*)
    module procedure scalar_multiply1_Q
    module procedure scalar_multiply2_Q
    module procedure outer_product_Q
  end interface

  !Multiplication operator
  interface operator(.dot.)
    module procedure inner_product_Q
  end interface

  !Division operator
  interface operator(/)
    module procedure scalar_division_Q
  end interface

contains

  type(Euler1D_Q_State) function add_Q(Q1, Q2)
    type(Euler1D_Q_State), intent(in) :: Q1, Q2
    add_Q%p = Q1%p + Q2%p
    add_Q%u = Q1%u + Q2%u
    add_Q%T = Q1%T + Q2%T
    return
  end function add_Q

  type(Euler1D_Q_State) function subtract_Q(Q1, Q2)
    type(Euler1D_Q_State), intent(in) :: Q1, Q2
    subtract_Q%p = Q1%p - Q2%p
    subtract_Q%u = Q1%u - Q2%u
    subtract_Q%T = Q1%T - Q2%T
    return
  end function subtract_Q

  type(Euler1D_Q_State) function scalar_multiply1_Q(a, Q)
    real(dp),              intent(in) :: a
    type(Euler1D_Q_State), intent(in) :: Q
    scalar_multiply1_Q%p = a*Q%p
    scalar_multiply1_Q%u = a*Q%u
    scalar_multiply1_Q%T = a*Q%T
    return
  end function scalar_multiply1_Q

  type(Euler1D_Q_State) function scalar_multiply2_Q(Q, a)
    type(Euler1D_Q_State), intent(in) :: Q
    real(dp),              intent(in) :: a
    scalar_multiply2_Q%p = a*Q%p
    scalar_multiply2_Q%u = a*Q%u
    scalar_multiply2_Q%T = a*Q%T
    return
  end function scalar_multiply2_Q

  type(Euler1D_Q_State) function outer_product_Q(Q1, Q2)
    type(Euler1D_Q_State), intent(in) :: Q1, Q2
    outer_product_Q%p = Q1%p*Q2%p
    outer_product_Q%u = Q1%u*Q2%u
    outer_product_Q%T = Q1%T*Q2%T
    return
  end function outer_product_Q

  real(dp) function inner_product_Q(Q1, Q2)
    type(Euler1D_Q_State), intent(in) :: Q1, Q2
    inner_product_Q = Q1%p*Q2%p + Q1%u*Q2%u + Q1%T*Q2%T
    return
  end function inner_product_Q

  type(Euler1D_Q_State) function scalar_division_Q(Q, a)
    type(Euler1D_Q_State), intent(in) :: Q
    real(dp),              intent(in) :: a
    scalar_division_Q%p = Q%p/a
    scalar_division_Q%u = Q%u/a
    scalar_division_Q%T = Q%T/a
    return
  end function scalar_division_Q

  !Set to vacuum state
  subroutine vacuum_Q(Q)
    type(Euler1D_Q_State), intent(out) :: Q
    Q%p = 0.0_dp
    Q%u = 0.0_dp
    Q%T = 0.0_dp
  end subroutine vacuum_Q

  !Set state to standard atmosphere
  subroutine standard_atmosphere_Q(Q)
    type(Euler1D_Q_State), intent(out) :: Q
    Q%p = 101325.0_dp
    Q%u = 0.0_dp
    Q%T = 288.15_dp
  end subroutine standard_atmosphere_Q

  !Set state to specified constant
  subroutine constant_Q(Q, val)
    type(Euler1D_Q_State), intent(out) :: Q
    real(dp),              intent(in)  :: val
    Q%p = val
    Q%u = val
    Q%T = val
  end subroutine constant_Q

  !Copy the specified state
  subroutine copy_Q(Q1, Q2)
    type(Euler1D_Q_State), intent(out) :: Q1
    type(Euler1D_Q_State), intent(in)  :: Q2
    Q1%p = Q2%p
    Q1%u = Q2%u
    Q1%T = Q2%T
  end subroutine copy_Q

  !Set the state to the specified values
  type(Euler1D_Q_State) function QState(p, u, T)
    real(dp), intent(in) :: p, u, T
    QState%p = p
    QState%u = u
    QState%T = T
  end function QState

  !Determine the minimum state values
  type(Euler1D_Q_State) function min_Q(Q1, Q2)
    type(Euler1D_Q_State), intent(in) :: Q1
    type(Euler1D_Q_State), intent(in) :: Q2
    min_Q%p = min(Q1%p,Q2%p)
    min_Q%u = min(Q1%u,Q2%u)
    min_Q%T = min(Q1%T,Q2%T)
  end function min_Q

  !Determine the maximum state values
  type(Euler1D_Q_State) function max_Q(Q1, Q2)
    type(Euler1D_Q_State), intent(in) :: Q1
    type(Euler1D_Q_State), intent(in) :: Q2
    max_Q%p = max(Q1%p,Q2%p)
    max_Q%u = max(Q1%u,Q2%u)
    max_Q%T = max(Q1%T,Q2%T)
  end function max_Q

  !Calculate the density of the gas
  real(dp) function rho_Q(Q)
    use gasConstants, only: R
    type(Euler1D_Q_State), intent(in) :: Q
    rho_Q = Q%p/(R*Q%T)
  end function rho_Q

  !Calculate the momentum of the gas
  real(dp) function du_Q(Q)
    type(Euler1D_Q_State), intent(in) :: Q
    real(dp) :: rho
    rho = rho_Q(Q)
    du_Q = rho*Q%u
  end function du_Q

  !Calculate the specific internal energy of the gas
  real(dp) function eps_Q(Q)
    use gasConstants, only: gm1i, R
    type(Euler1D_Q_State), intent(in) :: Q
    eps_Q = R*gm1i*Q%T
    return
  end function eps_Q

  !Calculate the total energy of the gas
  real(dp) function E_Q(Q)
    use gasConstants, only: gm1i
    type(Euler1D_Q_State), intent(in) :: Q
    real(dp) :: rho
    rho = rho_Q(Q)
    E_Q = gm1i*Q%p + 0.5_dp*rho*Q%u*Q%u
  end function E_Q

  !Calculate the specific internal energy of the gas
  real(dp) function Esp_Q(Q)
    use gasConstants, only: gm1i, R
    type(Euler1D_Q_State), intent(in) :: Q
    Esp_Q = R*gm1i*Q%T + 0.5_dp*Q%u*Q%u
    return
  end function Esp_Q

  !Calculate the specific internal enthalpy of the gas
  real(dp) function hint_Q(Q)
    use gasConstants, only: g, gm1i, R
    type(Euler1D_Q_State), intent(in) :: Q
    hint_Q = R*g*gm1i*Q%T
    return
  end function hint_Q

  !Calculate the total enthalpy of the gas
  real(dp) function H_Q(Q)
    use gasConstants, only: g, gm1i
    type(Euler1D_Q_State), intent(in) :: Q
    real(dp) :: rho
    rho = rho_Q(Q)
    H_Q = g*gm1i*Q%p + 0.5_dp*rho*Q%u*Q%u
    return
  end function H_Q

  !Calculate the specific enthalpy of the gas
  real(dp) function Hsp_Q(Q)
    use gasConstants, only: g, gm1i, R
    type(Euler1D_Q_State), intent(in) :: Q
    Hsp_Q = R*g*gm1i*Q%T + 0.5_dp*Q%u*Q%u
    return
  end function Hsp_Q

  !Calculate the speed of sound squared
  real(dp) function a2_Q(Q)
    use gasConstants, only: g, R
    type(Euler1D_Q_State), intent(in) :: Q
    a2_Q = g*R*Q%T
    return
  end function a2_Q

  !Calculate the speed of sound
  real(dp) function a_Q(Q)
    type(Euler1D_Q_State), intent(in) :: Q
    a_Q = sqrt(a2_Q(Q))
    return
  end function a_Q

  !Calculate the flow Mach number
  real(dp) function M_Q(Q)
    type(Euler1D_Q_State), intent(in) :: Q
    M_Q = abs(Q%u)/a_Q(Q)
    return
  end function M_Q

  !Calculate the entropy of the gas
  real(dp) function s_Q(Q)
    use gasConstants, only: g, gm1i, R
    type(Euler1D_Q_State), intent(in) :: Q
    real(dp) :: rho
    rho = rho_Q(Q)
    s_Q = R*gm1i*log(Q%p/(rho**g))
    return
  end function s_Q

  !Check for unphysical properties
  integer function unphysical_properties_Q(Q) result(ierr)
    type(Euler1D_Q_State), intent(in) :: Q
    if((Q%p.lt.0.0_dp).or.(Q%p.ne.Q%p)) then
     ierr = 1
    else if((Q%T.lt.0.0_dp).or.(Q%T.ne.Q%T)) then
     ierr = 1
    else if(E_Q(Q).lt.0.0_dp) then
     ierr = 1
    else
     ierr = 0
    end if
    return
  end function unphysical_properties_Q

  !Flux Vector
  subroutine Flux_Q(Q, F)
    type(Euler1D_Q_State), intent(in)  :: Q
    type(Euler1D_U_State), intent(out) :: F
    real(dp) :: rho, H
    rho = rho_Q(Q)
    H = H_Q(Q)
    F%rho = rho*Q%u
    F%du = rho*Q%u*Q%u + Q%p
    F%E = Q%u*H
    return
  end subroutine Flux_Q

  !Flux Jacobian
  subroutine dFdU_Q(Q, dFdU)
    use gasConstants, only: g, gm1
    type(Euler1D_Q_State),    intent(in)  :: Q
    real(dp), dimension(3,3), intent(out) :: dFdU
    real(dp) :: rho, E
    rho = rho_Q(Q)
    E = E_Q(Q)
    dFdU(1,1) = 0.0_dp
    dFdU(1,2) = 1.0_dp
    dFdU(1,3) = 0.0_dp
    dFdU(2,1) = 0.5_dp*(g-3.0_dp)*Q%u*Q%u
    dFdU(2,2) = (3.0_dp-g)*Q%u
    dFdU(2,3) = gm1
    dFdU(3,1) = gm1*Q%u*Q%u*Q%u - g*Q%u*E/rho
    dFdU(3,2) = g*E/rho - 0.5_dp*3.*gm1*Q%u*Q%u
    dFdU(3,3) = g*Q%u
  end subroutine dFdU_Q

  !Eigenvalue(s)
  subroutine lambda_Q(Q, lambda)
    type(Euler1D_Q_State), intent(in)  :: Q
    type(Euler1D_Q_State), intent(out) :: lambda
    real(dp) :: a
    a = a_Q(Q)
    lambda%p = Q%u-a
    lambda%u = Q%u
    lambda%T = Q%u+a
  end subroutine lambda_Q

  !Primitive right eigenvector
  subroutine rp_Q(Q, rp, index)
    type(Euler1D_Q_State), intent(in)  :: Q
    type(Euler1D_Q_State), intent(out) :: rp
    integer,               intent(in)  :: index
    real(dp) :: rho
    real(dp) :: a, a2
    rho = rho_Q(Q)
    a2 = a2_Q(Q)
    a = sqrt(a2)
    select case(index)
    case(1)
     rp%p = 1.0_dp
     rp%u = -a/rho
     rp%T = a2
    case(2)
     rp%p = 1.0_dp
     rp%u = 0.0_dp
     rp%T = 0.0_dp
    case(3)
     rp%p = 1.0_dp
     rp%u = a/rho
     rp%T = a2
    end select
  end subroutine rp_Q

  !Primitive left eigenvector
  subroutine lp_Q(Q, lp, index)
    type(Euler1D_Q_State), intent(in)  :: Q
    type(Euler1D_Q_State), intent(out) :: lp
    integer,               intent(in)  :: index
    real(dp) :: rho
    real(dp) :: a, a2
    rho = rho_Q(Q)
    a2 = a2_Q(Q)
    a = sqrt(a2)
    select case(index)
    case(1)
     lp%p = 0.0_dp
     lp%u = -0.5_dp*rho/a
     lp%T = 0.5_dp/a2
    case(2)
     lp%p = 1.0_dp
     lp%u = 0.0_dp
     lp%T = -1.0_dp/a2
    case(3)
     lp%p = 0.0_dp
     lp%u = 0.5_dp*rho/a
     lp%T = 0.5_dp/a2
    end select
  end subroutine lp_Q

  !Conserved right eigenvector
  subroutine rc_Q(Q, rc, index)
    type(Euler1D_Q_State), intent(in)  :: Q
    type(Euler1D_U_State), intent(out) :: rc
    integer, intent(in)                :: index
    real(dp)                     :: a, h
    a = a_Q(Q)
    h = Hsp_Q(Q)
    select case(index)
    case(1)
     rc%rho = 1.0_dp
     rc%du = Q%u - a
     rc%E = h - a*Q%u
    case(2)
     rc%rho = 1.0_dp
     rc%du = Q%u
     rc%E = 0.5_dp*Q%u*Q%u
    case(3)
     rc%rho = 1.0_dp
     rc%du = Q%u + a
     rc%E = h + a*Q%u
    end select
  end subroutine rc_Q

  !Reflection boundary condition
  subroutine Reflect_Q(Q1, Q2)
    type(Euler1D_Q_State), intent(out) :: Q1
    type(Euler1D_Q_State), intent(in)  :: Q2
    Q1%p = Q2%p
    Q1%u = - Q2%u
    Q1%T = Q2%T
  end subroutine Reflect_Q

end module Euler1D_QState
