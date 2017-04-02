!-----------------------------------------------------------------------
! Module defining gas-types.
!-----------------------------------------------------------------------
module gasConstants
  use realSizes, only: dp

  !Implicit declaration
  implicit none

  !Define gas types
  integer, parameter :: GAS_AIR = 0
  integer, parameter :: GAS_HE  = 1
  integer, parameter :: GAS_H2  = 2
  integer, parameter :: GAS_N2  = 3
  integer, parameter :: GAS_O2  = 4

  real(dp) :: xmw    !Molecular weight
  real(dp) :: g      !Ratio of specific heats
  real(dp) :: R      !Gas constant
  real(dp) :: gm1    !g-1
  real(dp) :: gm1i   !1/(g-1)
  real(dp) :: cv     !Specific heat at constant volume
  real(dp) :: cp     !Specific heat at constant pressure
  real(dp) :: gp1    !g+1
  real(dp) :: gp1i   !1/(g+1)
  real(dp) :: alpha  !
  real(dp) :: alphai !
  real(dp) :: beta   !
  real(dp) :: betai  !
  real(dp) :: Pr     !Prandtl number
  
  !Constants for viscosity and thermal conductivity correlations
  real(dp) :: c1
  real(dp) :: c2
  real(dp) :: c3
  real(dp) :: c4
  real(dp) :: c5

contains

  !Set the gas constants
  subroutine setGas(gas_type)
    use realSizes, only: dp
    implicit none
    integer, intent(inout) :: gas_type
    real(dp), parameter    :: Runiv = 8314.4598_dp
    if(gas_type.eq.GAS_AIR) then
      g   =    1.4_dp
      xmw =   28.9640_dp
      c1  =    0.000000521920_dp
      c2  =   -3.311320_dp
      c3  =    0.865351_dp
      c4  = 2365.270000_dp
      c5  =    0.000000_dp
    else if(gas_type.eq.GAS_HE) then
      g   =    1.6666666667_dp
      xmw =    4.0026_dp
      c1  =    0.000000395509_dp
      c2  =    0.000000_dp
      c3  =    0.813299_dp
      c4  =    0.000000_dp
      c5  =    0.000000_dp
    else if(gas_type.eq.GAS_H2) then
      g   =    1.405_dp
      xmw =    2.0160_dp
      c1  =    0.000000172256_dp
      c2  =    0.000000_dp
      c3  =    0.808749_dp
      c4  =    0.000000_dp
      c5  =    0.000000_dp
    else if(gas_type.eq.GAS_N2) then
      g   =    1.4_dp
      xmw =   28.0130_dp
      c1  =    0.000000484270_dp
      c2  =   -0.829421_dp
      c3  =    0.851275_dp
      c4  = 1219.260000_dp
      c5  =    0.000000_dp
    else if(gas_type.eq.GAS_O2) then
      g   =    1.395_dp
      xmw =   32.0000_dp
      c1  =    0.000000540279_dp
      c2  =   -0.164235_dp
      c3  =    0.851096_dp
      c4  = 2049.970000_dp
      c5  =    0.000000_dp
    end if
    R  =  Runiv/xmw
    gm1 = g-1.0_dp
    gm1i = 1.0_dp/gm1
    cv = R*gm1i
    cp = g*R*gm1i
    gp1 = g+1.0_dp
    gp1i = 1.0_dp/gp1
    alpha = gp1/gm1
    alphai = gm1/gp1
    beta = 0.5_dp*gm1/g
    betai = 2.0_dp*g/gm1
    Pr = 2.0_dp*g/(3.9_dp*g - 1.5_dp) !Euken's formula
    return
  end subroutine setGas

  !Dynamic viscosity correlation from James Gottlieb [kg/(m-s)]
  real(dp) function viscosity(T)
    implicit none
    real(dp), intent(in) :: T
    viscosity = c1*(T**1.5_dp)/(c2 + T**c3 + c4/T) + c5
    return
  end function viscosity

  !Thermal conductivity correlation from James Gottlieb [W/(m-K)]
  real(dp) function kappa(T)
    implicit none
    real(dp), intent(in) :: T
    kappa = viscosity(T)*cp/Pr
  end function kappa

end module gasConstants
