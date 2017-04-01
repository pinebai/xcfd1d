module exactSoln_module
  use realsizes, only: dp
  use Euler1D_WState
  implicit none

  integer                                          :: ne
  real(dp),              allocatable, dimension(:) :: xe
  type(Euler1D_W_State), allocatable, dimension(:) :: We

contains

  !Deallocate routine -- allocation done in specific solution routines
  subroutine exactSoln_deallocate
    implicit none
    if(allocated(xe)) deallocate(xe)
    if(allocated(We)) deallocate(We)
    return
  end subroutine exactSoln_deallocate

  !----------------------------------------------------------------------
  !Sod Problem: (x,t)     = (0.50,0.25), x in [0,1]
  !             (d,u,p)_l = (1.000,0,1.0)
  !             (d,u,p)_r = (0.125,0,0.1)
  !---------------------------------------------------------------------
  subroutine exactSod(t)
    use realsizes, only: dp
    use mathFunctions, only: sqr
    use numbers, only: TOLER
    use gasConstants, only: g, gm1i, gp1i, alpha, beta, alphai, betai
    use Euler1D_WState
    implicit none
    !Argument variables:
    real(dp), intent(in)  :: t
    !Local variables:
    integer               :: i
    real(dp)              :: xo, dx
    real(dp)              :: al, ar, als, ars, a
    type(Euler1D_W_State) :: Wl, Wr, Wls, Wrs, Wint
    !Allocate exact solution arrays
    ne = 27
    allocate(xe(ne))
    allocate(We(ne))
    !Location of jump conditions
    xo = 0.5_dp
    !Initialize the left state
    Wl = WState(1.0_dp,0.0_dp,1.0_dp)
    al = a_W(Wl)
    !Initialize the right state
    Wr = WState(0.125_dp,0.0_dp,0.1_dp)
    ar = a_W(Wr)
    !If initial conditions
    if(t.lt.TOLER) then
      xe(1) = 0.0_dp
      xe(2:26) = xo
      xe(27) = 1.0_dp
      We(1:25) = Wl
      We(26:27) = Wr
      return
    end if
    !Solve shocktube/Riemann problem
    call Riemann(Wl,Wr,Wls,Wrs,Wint,0)
    als = a_W(Wls)
    ars = a_W(Wrs)
    !1    - Left-side of domain
    xe(1) = 0.0_dp
    We(1) = Wl
    !2-22 - Head-to-tail of rarefaction wave
    dx = ((Wls%u - als) - (Wl%u - al))*t/20.0_dp
    do i = 2, 22
      xe(i) = xo + (Wl%u - al)*t + dx*(i-2)
      We(i)%u = alphai*(Wl%u + 2.0_dp*gm1i*al) + 2.0_dp*gp1i*(xe(i)-xo)/t
      a = We(i)%u - (xe(i)-xo)/t
      We(i)%p = Wl%p*((a/al)**betai)
      We(i)%rho = g*We(i)%p/sqr(a)
    end do
    !23   - Contact surface (rarefaction side)
    xe(23) = xo + Wls%u*t
    We(23) = Wls
    !24   - Contact surface (shock side)
    xe(24) = xe(23)
    We(24) = Wrs
    !25   - Shock wave (contact surface side)
    xe(25) = xo + (Wrs%u + ars*sqrt(beta*alpha*Wr%p/Wrs%p + beta))*t
    We(25) = Wrs
    !26   - Shock wave (silent side)
    xe(26) = xe(25)
    We(26) = Wr
    !27   - Right-side of domain
    xe(27) = 1.0_dp
    We(27) = Wr
    return
  end subroutine exactSod

  !-------------------------------------------------------------------
  !Modified Sod Problem
  !(x,t) = (0.30,0.20), x in [0,1]
  !Wl = Euler1D_pState(1.000,0.75,1.0)
  !Wr = Euler1D_pState(0.125,0.00,0.1)
  !-------------------------------------------------------------------
  subroutine exactModifiedSod(t)
    use realsizes, only: dp
    use mathFunctions, only: sqr
    use numbers, only: TOLER
    use gasConstants, only: g, gm1i, gp1i, alpha, beta, alphai, betai
    use Euler1D_WState
    implicit none
    !Argument variables:
    real(dp), intent(in)  :: t
    !Local variables:
    integer               :: i
    real(dp)              :: xo, dx
    real(dp)              :: al, ar, als, ars, a
    type(Euler1D_W_State) :: Wl, Wr, Wls, Wrs, Wint
    !Allocate exact solution arrays
    ne = 27
    allocate(xe(ne))
    allocate(We(ne))
    !Initial location of jump conditions
    xo = 0.3_dp
    !Initialize the left state
    Wl = WState(1.0_dp,0.75_dp,1.0_dp)
    !Wl = WState(1.0_dp,280.624_dp,1.0d4)
    al = a_W(Wl)
    !Initialize the right state
    Wr = WState(0.125_dp,0.0_dp,0.1_dp)
    !Wr = WState(0.125_dp,0.0_dp,1.0d4)
    ar = a_W(Wr)
    !If initial conditions
    if(t.lt.TOLER) then
      xe(1) = 0.0_dp
      xe(2:26) = xo
      xe(27) = 1.0_dp
      We(1:25) = Wl
      We(25:27) = Wr
      return
    end if
    !Solve shocktube/Riemann problem
    call Riemann(Wl,Wr,Wls,Wrs,Wint,0)
    als = a_W(Wls)
    ars = a_W(Wrs)
    !1    - Left-side of domain
    xe(1) = 0.0_dp
    We(1) = Wl
    !2-22 - Head-to-tail of rarefaction wave
    dx = ((Wls%u - als) - (Wl%u - al))*t/20.0_dp
    do i = 2, 22
      xe(i) = xo + (Wl%u - al)*t + dx*(i-2)
      We(i)%u = alphai*(Wl%u + 2.0_dp*gm1i*al) + 2.0_dp*gp1i*(xe(i)-xo)/t
      a = We(i)%u - (xe(i)-xo)/t
      We(i)%p = Wl%p*((a/al)**betai)
      We(i)%rho = g*We(i)%p/sqr(a)
    end do
    !23   - Contact surface (rarefaction side)
    xe(23) = xo + Wls%u*t
    We(23) = Wls
    !24   - Contact surface (shock side)
    xe(24) = xe(23)
    We(24) = Wrs
    !25   - Shock wave (contact surface side)
    xe(25) = xo + (Wrs%u + ars*sqrt(beta*alpha*Wr%p/Wrs%p + beta))*t
    We(25) = Wrs
    !26   - Shock wave (silent side)
    xe(26) = xe(25)
    We(26) = Wr
    !27   - Right-side of domain
    xe(27) = 1.0_dp
    We(27) = Wr
    return
  end subroutine exactModifiedSod

  !-------------------------------------------------------------------
  !123 Problem
  !(x,t) = (0.5,0.15), x in [0,1]
  !Wl = Euler1D_pState(1.0,-2.0,0.4)
  !Wr = Euler1D_pState(1.0, 2.0,0.4)
  !-------------------------------------------------------------------
  subroutine Exact123Problem(t)
    use realsizes, only: dp
    use mathFunctions, only: sqr
    use numbers, only: TOLER
    use gasConstants, only: g, gm1i, gp1i, alphai, betai
    use Euler1D_WState
    implicit none
    !Argument variables:
    real(dp), intent(in)  :: t
    !Local variables:
    integer               :: i
    real(dp)              :: xo, dx
    real(dp)              :: al, ar, als, ars, a
    type(Euler1D_W_State) :: Wl, Wr, Wls, Wrs, Wint
    !Allocate exact solution arrays
    ne = 44
    allocate(xe(ne))
    allocate(We(ne))
    !Initial location of jump conditions
    xo = 0.5_dp
    !Initialize the left state
    Wl = WState(1.0_dp,-2.0_dp,0.4_dp)
    al = a_W(Wl)
    !Initialize the right state
    Wr = WState(1.0_dp,2.0_dp,0.4_dp)
    ar = a_W(Wr)
    !If initial conditions
    if(t.lt.TOLER) then
      xe(1) = 0.0_dp
      xe(2:43) = xo
      xe(44) = 1.0_dp
      We(1:22) = Wl
      We(23:44) = Wr
      return
    end if
    !Solve shocktube/Riemann problem
    call Riemann(Wl,Wr,Wls,Wrs,Wint,0)
    als = a_W(Wls)
    ars = a_W(Wrs)
    !1    - Left-side of domain
    xe(1) = 0.0_dp
    We(1) = Wl
    !2-22 - Head-to-tail of rarefation wave
    dx = ((Wls%u - als) - (Wl%u - al))*t/20.0_dp
    do i = 2, 22
      xe(i) = xo + (Wl%u - al)*t + dx*(i-2)
      We(i)%u = alphai*(Wl%u + 2.0_dp*gm1i*al) + 2.0_dp*gp1i*(xe(i)-xo)/t
      a = We(i)%u - (xe(i)-xo)/t
      We(i)%p = Wl%p*((a/al)**betai)
      We(i)%rho = g*We(i)%p/sqr(a)
    end do
    !23-43 - Tail-to-head of rarefation wave
    dx = ((Wr%u + ar) - (Wrs%u + ars))*t/20.0_dp
    do i = 43, 23, -1
      xe(i) = xo + (Wr%u + ar)*t + dx*(i-43)
      We(i)%u = 2.0_dp*gp1i*(xe(i)-xo)/t - alphai*(Wrs%u + 2.0_dp*gm1i*ars)
      a = (xe(i)-xo)/t - We(i)%u
      We(i)%p = Wr%p*((a/ar)**betai)
      We(i)%rho = g*We(i)%p/sqr(a)
    end do
    !44   - Right-side of domain
    xe(44) = 1.0_dp
    We(44) = Wr
    return
  end subroutine exact123Problem

  !-------------------------------------------------------------------
  !Left Rarefaction, Right Shock Problem
  !(x,t) = (0.5,0.012), x in [0,1]
  !Wl = Euler1D_pState(1.0,0.0,0.01)
  !Wr = Euler1D_pState(1.0,0.0,1000.0)
  !-------------------------------------------------------------------
  subroutine exactStrongSod(t)
    use realsizes, only: dp
    use mathFunctions, only: sqr
    use numbers, only: TOLER
    use gasConstants, only: g, gm1i, gp1i, alpha, beta, alphai, betai
    use Euler1D_WState
    implicit none
    !Argument variables:
    real(dp), intent(in)  :: t
    !Local variables:
    integer               :: i
    real(dp)              :: xo, dx
    real(dp)              :: al, ar, als, ars, a
    type(Euler1D_W_State) :: Wl, Wr, Wls, Wrs, Wint
    !Allocate exact solution arrays
    ne = 27
    allocate(xe(ne))
    allocate(We(ne))
    !Initial location of jump conditions
    xo = 0.5_dp
    !Initialize the left state
    Wl = WState(1.0_dp,0.0_dp,1000.0_dp)
    al = a_W(Wl)
    !Initialize the right state
    Wr = WState(0.125_dp,0.0_dp,0.01_dp)
    ar = a_W(Wr)
    !If initial conditions
    if(t.lt.TOLER) then
      xe(1) = 0.0_dp
      xe(2:26) = xo
      xe(27) = 1.0_dp
      We(1:25) = Wl
      We(26:27) = Wr
      return
    end if
    !Solve shocktube/Riemann problem
    call Riemann(Wl,Wr,Wls,Wrs,Wint,0)
    als = a_W(Wls)
    ars = a_W(Wrs)
    !1    - Left-side of domain
    xe(1) = 0.0_dp
    We(1) = Wl
    !2-22 - Head-to-tail of rarefaction wave
    dx = ((Wls%u - als) - (Wl%u - al))*t/20.0_dp
    do i = 2, 22
      xe(i) = xo + (Wl%u - al)*t + dx*(i-2)
      We(i)%u = alphai*(Wl%u + 2.0_dp*gm1i*al) + 2.0_dp*gp1i*(xe(i)-xo)/t
      a = We(i)%u - (xe(i)-xo)/t
      We(i)%p = Wl%p*((a/al)**betai)
      We(i)%rho = g*We(i)%p/sqr(a)
    end do
    !23   - Contact surface (rarefaction side)
    xe(23) = xo + Wls%u*t
    We(23) = Wls
    !24   - Contact surface (shock side)
    xe(24) = xe(23)
    We(24) = Wrs
    !25   - Shock wave (contact surface side)
    xe(25) = xo + (Wrs%u + ars*sqrt(beta*alpha*Wr%p/Wrs%p + beta))*t
    We(25) = Wrs
    !26   - Shock wave (silent side)
    xe(26) = xe(25)
    We(26) = Wr
    !27   - Right-side of domain
    xe(27) = 1.0_dp
    We(27) = Wr
    return
  end subroutine exactStrongSod

  !---------------------------------------------------------------------
  !Left Rarefaction, Stationary Contact, Right Shock Problem
  !(x,t) = (0.8,0.012), x in [0,1]
  !Wl = Euler1D_pState(1.0,-19.59745,1000.0)
  !Wr = Euler1D_pState(1.0,-19.59745,0.01)
  !---------------------------------------------------------------------
  subroutine exactStationaryContact(t)
    use realsizes, only: dp
    use mathFunctions, only: sqr
    use numbers, only: TOLER
    use gasConstants, only: g, gm1i, gp1i, alpha, beta, alphai, betai
    use Euler1D_WState
    implicit none
    !Argument variables:
    real(dp), intent(in)  :: t
    !Local variables:
    integer               :: i
    real(dp)              :: xo, dx
    real(dp)              :: al, ar, als, ars, a
    type(Euler1D_W_State) :: Wl, Wr, Wls, Wrs, Wint
    !Allocate exact solution arrays
    ne = 27
    allocate(xe(ne))
    allocate(We(ne))
    !Initial location of jump conditions
    xo = 0.8_dp
    !Initialize the left state
    Wl = WState(1.0_dp,-19.59745_dp,1000.0_dp)
    al = a_W(Wl)
    !Initialize the right state
    Wr = WState(0.125_dp,-19.59745_dp,0.01_dp)
    ar = a_W(Wr)
    !If initial conditions
    if(t.lt.TOLER) then
      xe(1) = 0.0_dp
      xe(2:26) = xo
      xe(27) = 1.0_dp
      We(1:25) = Wl
      We(26:27) = Wr
      return
    end if
    !Set and solve shocktube/Riemann problem
    call Riemann(Wl,Wr,Wls,Wrs,Wint,0)
    als = a_W(Wls)
    ars = a_W(Wrs)
    !1    - Left-side of domain
    xe(1) = 0.0_dp
    We(1) = Wl
    !2-22 - Head-to-tail of rarefaction wave
    dx = ((Wls%u - als) - (Wl%u - al))*t/20.0_dp
    do i = 2, 22
      xe(i) = xo + (Wl%u - al)*t + dx*(i-2)
      We(i)%u = alphai*(Wl%u + 2.0_dp*gm1i*al) + 2.0_dp*gp1i*(xe(i)-xo)/t
      a = We(i)%u - (xe(i)-xo)/t
      We(i)%p = Wl%p*((a/al)**betai)
      We(i)%rho = g*We(i)%p/sqr(a)
    end do
    !23   - Contact surface (rarefaction side)
    xe(23) = xo + Wls%u*t
    We(23) = Wls
    !24   - Contact surface (shock side)
    xe(24) = xe(23)
    We(24) = Wrs
    !25   - Shock wave (contact surface side)
    xe(25) = xo + (Wrs%u + ars*sqrt(beta*alpha*Wr%p/Wrs%p + beta))*t
    We(25) = Wrs
    !26   - Shock wave (silent side)
    xe(26) = xe(25)
    We(26) = Wr
    !27   - Right-side of domain
    xe(27) = 1.0_dp
    We(27) = Wr
    return
  end subroutine exactStationaryContact

  !---------------------------------------------------------------------
  !Three Right, Discontinuities Problem
  !(x,t) = (0.4,0.035), x in [0,1]
  !Wl = Euler1D_pState(5.99924,19.59750,460.894)
  !Wr = Euler1D_pState(5.99242,-6.19633, 46.095)
  !---------------------------------------------------------------------
  subroutine exact3RightWaves(t)
    use realsizes, only: dp
    use numbers, only: TOLER
    use gasConstants, only: alpha, beta
    use Euler1D_WState
    implicit none
    !Argument variables:
    real(dp), intent(in)  :: t
    !Local variables:
    integer               :: i
    real(dp)              :: xo, dx
    real(dp)              :: al, ar, als, ars, a
    type(Euler1D_W_State) :: Wl, Wr, Wls, Wrs, Wint
    !Allocate exact solution arrays
    ne = 8
    allocate(xe(ne))
    allocate(We(ne))
    !Initial location of jump conditions
    xo = 0.4_dp
    !Initialize the left state
    Wl = WState(5.99924_dp,19.59750_dp,460.894_dp)
    al = a_W(Wl)
    !Initialize the right state
    Wr = WState(5.99242_dp,-6.19633_dp,46.095_dp)
    ar = a_W(Wr)
    !If initial conditions
    if(t.lt.TOLER) then
      xe(1) = 0.0_dp
      xe(2:7) = xo
      xe(8) = 1.0_dp
      We(1:6) = Wl
      We(7:8) = Wr
      return
    end if
    !Solve shocktube/Riemann problem
    call Riemann(Wl,Wr,Wls,Wrs,Wint,0)
    als = a_W(Wls)
    ars = a_W(Wrs)
    !1 - Left-side of domain
    xe(1) = 0.0_dp
    We(1) = Wl
    !2 - Before first shock-wave
    xe(2) = xo + t*(Wl%u - al*sqrt(beta*alpha*Wls%p/Wl%p + beta))
    We(2) = Wl
    !3 - After first shock-wave
    xe(3) = xo + t*(Wl%u - al*sqrt(beta*alpha*Wls%p/Wl%p + beta))
    We(3) = Wls
    !4 - Before contact-wave
    xe(4) = xo + t*Wls%u
    We(4) = Wls
    !5 - After contact-wave
    xe(5) = xo + t*Wls%u
    We(5) = Wrs
    !6 - Before last shock-wave
    xe(6) = xo + t*(Wrs%u + ars*sqrt(beta*alpha*Wr%p/Wrs%p + beta))
    We(6) = Wrs
    !7 - After last shock-wave
    xe(7) = xo + t*(Wrs%u + ars*sqrt(beta*alpha*Wr%p/Wrs%p + beta))
    We(7) = Wr
    !8 - Right-side of domain
    xe(8) = 1.0_dp
    We(8) = Wr
    return
  end subroutine exact3RightWaves

  !---------------------------------------------------------------------
  !Subsonic or Transonic Nozzle
  !Subsonic: (Sthroat,xShock) = (0.8,11.0)
  !Transonic: (Sthroat,xShock) = (1.0,7.0)
  !---------------------------------------------------------------------
  recursive subroutine Nozzle(ic_type,t)
    use realsizes, only: dp
    use numbers, only: TOLER
    use cfdParams, only: IC_SUBSONIC_NOZZLE, IC_TRANSONIC_NOZZLE
    use gasConstants, only: g, gm1, gp1, gm1i, gp1i, R, alpha, alphai, betai
    use Euler1D_WState
    implicit none
    !Argument variables:
    integer,  intent(in)  :: ic_type
    real(dp), intent(in)  :: t
    !Local variables:
    real(dp)              :: xs, Ms !Location and Mach number at sonic point
    type(Euler1D_W_State) :: Ws     !Sonic point state variables
    real(dp)              :: To, po, Mo
    real(dp)              :: f, fp, Sstar
    real(dp)              :: Mr, pr
!!$    !Allocate exact solution arrays:
!!$    if(IC.eq.IC_SUBSONIC_NOZZLE) then
!!$      ne = 101
!!$      Mo = 0.1_dp
!!$      To = 300.0_dp
!!$      po = 100000.0_dp
!!$      Sstar = 0.8_dp
!!$      xs = 11.0_dp
!!$    else
!!$      ne = 102
!!$      Sstar = 1.0_dp
!!$      xs = 7.0_dp
!!$    end if
!!$    allocate(xe(ne))
!!$    allocate(We(ne))
!!$    !Set constants:
!!$    Mo = 0.1_dp
!!$    To = 300.0_dp
!!$    po = 100000.0_dp
!!$    do n = 1, ne
!!$      xe(n) = 0.0_dp + 10.0_dp*(n-1.0_dp)/(ne-1.0_dp)
!!$      S = 1.0_dp + 0.5_dp*(1.0_dp - xe(n)/5.0_dp)**2
!!$      if(IC.eq.IC_TRANSONIC_NOZZLE) then
!!$        Mo = 0.1_dp
!!$        To = 300.0_dp
!!$        po = 100000.0_dp
!!$        if((x.gt.5.0_dp).and.(x.le.xs)) then
!!$          M = 1.2_dp
!!$        else if(x.gt.xs) then
!!$          M = 0.1_dp
!!$          call Nozzle(xs,1.0_dp+0.5_dp*((1.0_dp-xs/5.0_dp)**2),IC,Wl)
!!$          Ml = M_W(Wl)
!!$          Mr = sqrt((2.0_dp + gm1*sqr(Ml))/(2.0_dp*g*sqr(Ml)-gm1))
!!$          pr = Wl%p*(2.0_dp*g*sqr(Ml)-gm1)/gp1
!!$          po = po*(((gp1*sqr(Ml))/(2.0_dp+gm1*sqr(Ml)))**(g*gm1i))/ &
!!$                  ((alphai*(betai*sqr(Ml)-1.0_dp))**gm1i)
!!$          po = pr + 0.5_dp*sqr(Mr)*g*pr
!!$        end if
!!$      end if
!!$      !Iterate to find actual M.
!!$      do
!!$        f  = (S/Sstar)*M - ((2.0_dp*gp1i*(1.0_dp + 0.5_dp*gm1*sqr(M)))**(alpha/2.0_dp))
!!$        fp = (S/Sstar) - M*((2.0_dp*gp1i*(1.0_dp + 0.5_dp*gm1*sqr(M)))**((alpha/2.0_dp) - 1.0_dp))
!!$        M = M - f/fp
!!$        if(abs(f).le.TOLER) exit
!!$      end do
!!$      !Compute the state variables.
!!$      T         = To/(1.0_dp + 0.5_dp*gm1*sqr(M))
!!$      a         = sqrt(g*R*T)
!!$      We(n)%u   = M*a
!!$      We(n)%p   = po*((T/To)**(g*gm1i))
!!$      We(n)%rho = W%p/(R*T)
!!$    end do
    return
  end subroutine Nozzle

  !---------------------------------------------------------------------
  !Square Density Wave
  !Wl = Euler1D_pState(1.225,100.0,101325.0)
  !Wm = Euler1D_pState(2.450,100.0,101325.0)
  !Wr = Euler1D_pState(1.225,100.0,101325.0)
  !---------------------------------------------------------------------
  subroutine Square_Wave(t)
    use realsizes, only: dp
    use Euler1D_WState
    implicit none
    !Argument variables:
    real(dp), intent(in)  :: t
    !Local variables:
    integer               :: n
    real(dp)              :: xl, xr
    type(Euler1D_W_State) :: W, Wm
    !Allocate exact solution arrays
    ne = 101
    allocate(xe(ne))
    allocate(We(ne))
    !Initialize the wave states
    W%rho = 1.0_dp ; W%u = 100.0_dp ; W%p = 100000.0_dp
    Wm%rho = 2.0_dp ; Wm%u = 100.0_dp ; Wm%p = 100000.0_dp
    !Set left/right location of wave
    xl = 0.4_dp + W%u*t
    do while(xl.gt.1.0_dp)
      xl = xl - 1.0_dp
    end do
    xr = 0.6_dp + W%u*t
    do while(xr.gt.1.0_dp)
      xr = xr - 1.0_dp
    end do
    !Calculate the primitive solution state based on the location and time
    do n = 1, ne
      xe(n) = 0.0_dp + (n-1.0_dp)/(ne-1.0_dp)
      We(n) = W
      if((xe(n).gt.xl).or.(xe(n).lt.xr)) We(n) = Wm
    end do
    return
  end subroutine Square_Wave

  !---------------------------------------------------------------------
  !Sine-Squared Density Wave
  !---------------------------------------------------------------------
  subroutine Sine_Squared_Wave(x,t,W)
    use realsizes, only: dp
    use numbers, only: PI
    use Euler1D_WState
    implicit none
    real(dp), intent(in)               :: x, t
    type(Euler1D_W_State), intent(out) :: W
    !Set the primitive solution state
    W%rho = 1.225_dp ; W%u = 100.0_dp ; W%p = 101325.0_dp
    !Update the density based on the location and time
    if(6.0_dp*x.le.1.0_dp + 6.0_dp*W%u*t) then
    else if(x.le.0.5_dp + W%u*t) then
      W%rho = W%rho*(1.0_dp + (sin(0.5*PI*(6.0_dp*(x-W%u*t)-1.0_dp))**2))
    else
    end if
    return
  end subroutine Sine_Squared_Wave

  !---------------------------------------------------------------------
  !Semi-Ellipse Density Wave
  !---------------------------------------------------------------------
  subroutine Semi_Ellipse_Wave(x,t,W)
    use realsizes, only: dp
    use Euler1D_WState
    implicit none
    real(dp), intent(in)               :: x, t
    type(Euler1D_W_State), intent(out) :: W
    !Set the primitive solution state.
    W%rho = 1.225_dp ; W%u = 100.0_dp ; W%p = 101325.0_dp
    !Update the density based on the location and time
    if(6.0_dp*x.le.1.0_dp + 6.0_dp*W%u*t) then
    else if(x.le.0.5_dp + W%u*t) then
      W%rho = W%rho*(1.0_dp + sqrt(1.0_dp-(6.0_dp*(x-W%u*t)-2.0_dp)**2))
    else
    end if
    return
  end subroutine Semi_Ellipse_Wave


end module exactSoln_module
