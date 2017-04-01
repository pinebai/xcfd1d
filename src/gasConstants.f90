!-----------------------------------------------------------------------
! Module defining gas-types.
!-----------------------------------------------------------------------
module gasConstants
  use realSizes, only: dp

  !Implicit declaration.
  implicit none

  !Define boundary conditions:
  integer, parameter :: GAS_AIR = 0

  real(dp) :: R      !Gas constant.
  real(dp) :: g      !Ratio of specific heats.
  real(dp) :: gm1    !g-1
  real(dp) :: gm1i   !1/(g-1)
  real(dp) :: cv     !Specific heat at constant volume.
  real(dp) :: cp     !Specific heat at constant pressure.
  real(dp) :: gp1    !g+1
  real(dp) :: gp1i   !1/(g+1)
  real(dp) :: alpha  !
  real(dp) :: alphai !
  real(dp) :: beta   !
  real(dp) :: betai  !

contains

  subroutine setGas(gas_type)
    use realSizes, only: dp
    implicit none
    integer, intent(in) :: gas_type
    if(gas_type.eq.GAS_AIR) then
      R = 287.06_dp
      g = 1.4_dp
    else
      R = 287.06_dp
      g = 1.4_dp
    end if
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
    return
  end subroutine setGas

end module gasConstants
