module exactSoln_module
  use realsizes, only: dp
  use gridBlock_module
  use Euler1D_WState
  implicit none

  integer                                          :: ne
  type(cellGeom)       , allocatable, dimension(:) :: Xe
  type(Euler1D_W_State), allocatable, dimension(:) :: We

contains

  !Deallocate routine -- allocation done in specific solution routines
  subroutine exactSoln_deallocate
    implicit none
    if(allocated(Xe)) deallocate(Xe)
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
    integer               :: n
    real(dp)              :: xo, dx
    real(dp)              :: al, ar, als, ars, a
    type(Euler1D_W_State) :: Wl, Wr, Wls, Wrs, Wint
    !Allocate exact solution arrays
    ne = 27
    allocate(Xe(ne))
    allocate(We(ne))
    do n = 1, ne
      Xe(n)%xc = 0.0_dp ; Xe(n)%dx = 0.0_dp
      Xe(n)%xa = 0.0_dp ; Xe(n)%da = 0.0_dp
      call Vacuum_W(We(n))
    end do
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
      Xe(1)%xc = 0.0_dp
      Xe(2:26)%xc = xo
      Xe(27)%xc = 1.0_dp
      We(1:25) = Wl
      We(26:27) = Wr
      return
    end if
    !Solve shocktube/Riemann problem
    call Riemann(Wl,Wr,Wls,Wrs,Wint,0)
    als = a_W(Wls)
    ars = a_W(Wrs)
    !1    - Left-side of domain
    Xe(1)%xc = 0.0_dp
    We(1) = Wl
    !2-22 - Head-to-tail of rarefaction wave
    dx = ((Wls%u - als) - (Wl%u - al))*t/20.0_dp
    do n = 2, 22
      Xe(n)%xc = xo + (Wl%u - al)*t + dx*(n-2)
      We(n)%u = alphai*(Wl%u + 2.0_dp*gm1i*al) + 2.0_dp*gp1i*(Xe(n)%xc-xo)/t
      a = We(n)%u - (Xe(n)%xc-xo)/t
      We(n)%p = Wl%p*((a/al)**betai)
      We(n)%rho = g*We(n)%p/sqr(a)
    end do
    !23   - Contact surface (rarefaction side)
    xe(23)%xc = xo + Wls%u*t
    We(23) = Wls
    !24   - Contact surface (shock side)
    xe(24)%xc = Xe(23)%xc
    We(24) = Wrs
    !25   - Shock wave (contact surface side)
    xe(25)%xc = xo + (Wrs%u + ars*sqrt(beta*alpha*Wr%p/Wrs%p + beta))*t
    We(25) = Wrs
    !26   - Shock wave (silent side)
    xe(26)%xc = Xe(25)%xc
    We(26) = Wr
    !27   - Right-side of domain
    xe(27)%xc = 1.0_dp
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
    integer               :: n
    real(dp)              :: xo, dx
    real(dp)              :: al, ar, als, ars, a
    type(Euler1D_W_State) :: Wl, Wr, Wls, Wrs, Wint
    !Allocate exact solution arrays
    ne = 27
    allocate(Xe(ne))
    allocate(We(ne))
    do n = 1, ne
      Xe(n)%xc = 0.0_dp ; Xe(n)%dx = 0.0_dp
      Xe(n)%xa = 0.0_dp ; Xe(n)%da = 0.0_dp
      call Vacuum_W(We(n))
    end do
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
      Xe(1)%xc = 0.0_dp
      Xe(2:26)%xc = xo
      Xe(27)%xc = 1.0_dp
      We(1:25) = Wl
      We(25:27) = Wr
      return
    end if
    !Solve shocktube/Riemann problem
    call Riemann(Wl,Wr,Wls,Wrs,Wint,0)
    als = a_W(Wls)
    ars = a_W(Wrs)
    !1    - Left-side of domain
    Xe(1)%xc = 0.0_dp
    We(1) = Wl
    !2-22 - Head-to-tail of rarefaction wave
    dx = ((Wls%u - als) - (Wl%u - al))*t/20.0_dp
    do n = 2, 22
      Xe(n)%xc = xo + (Wl%u - al)*t + dx*(n-2)
      We(n)%u = alphai*(Wl%u + 2.0_dp*gm1i*al) + 2.0_dp*gp1i*(Xe(n)%xc-xo)/t
      a = We(n)%u - (Xe(n)%xc-xo)/t
      We(n)%p = Wl%p*((a/al)**betai)
      We(n)%rho = g*We(n)%p/sqr(a)
    end do
    !23   - Contact surface (rarefaction side)
    Xe(23)%xc = xo + Wls%u*t
    We(23) = Wls
    !24   - Contact surface (shock side)
    Xe(24)%xc = Xe(23)%xc
    We(24) = Wrs
    !25   - Shock wave (contact surface side)
    Xe(25)%xc = xo + (Wrs%u + ars*sqrt(beta*alpha*Wr%p/Wrs%p + beta))*t
    We(25) = Wrs
    !26   - Shock wave (silent side)
    Xe(26)%xc = Xe(25)%xc
    We(26) = Wr
    !27   - Right-side of domain
    Xe(27)%xc = 1.0_dp
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
    integer               :: n
    real(dp)              :: xo, dx
    real(dp)              :: al, ar, als, ars, a
    type(Euler1D_W_State) :: Wl, Wr, Wls, Wrs, Wint
    !Allocate exact solution arrays
    ne = 44
    allocate(Xe(ne))
    allocate(We(ne))
    do n = 1, ne
      Xe(n)%xc = 0.0_dp ; Xe(n)%dx = 0.0_dp
      Xe(n)%xa = 0.0_dp ; Xe(n)%da = 0.0_dp
      call Vacuum_W(We(n))
    end do
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
      Xe(1)%xc = 0.0_dp
      Xe(2:43)%xc = xo
      Xe(44)%xc = 1.0_dp
      We(1:22) = Wl
      We(23:44) = Wr
      return
    end if
    !Solve shocktube/Riemann problem
    call Riemann(Wl,Wr,Wls,Wrs,Wint,0)
    als = a_W(Wls)
    ars = a_W(Wrs)
    !1    - Left-side of domain
    Xe(1)%xc = 0.0_dp
    We(1) = Wl
    !2-22 - Head-to-tail of rarefation wave
    dx = ((Wls%u - als) - (Wl%u - al))*t/20.0_dp
    do n = 2, 22
      Xe(n)%xc = xo + (Wl%u - al)*t + dx*(n-2)
      We(n)%u = alphai*(Wl%u + 2.0_dp*gm1i*al) + 2.0_dp*gp1i*(Xe(n)%xc-xo)/t
      a = We(n)%u - (Xe(n)%xc-xo)/t
      We(n)%p = Wl%p*((a/al)**betai)
      We(n)%rho = g*We(n)%p/sqr(a)
    end do
    !23-43 - Tail-to-head of rarefation wave
    dx = ((Wr%u + ar) - (Wrs%u + ars))*t/20.0_dp
    do n = 43, 23, -1
      Xe(n)%xc = xo + (Wr%u + ar)*t + dx*(n-43)
      We(n)%u = 2.0_dp*gp1i*(Xe(n)%xc-xo)/t - alphai*(Wrs%u + 2.0_dp*gm1i*ars)
      a = (Xe(n)%xc-xo)/t - We(n)%u
      We(n)%p = Wr%p*((a/ar)**betai)
      We(n)%rho = g*We(n)%p/sqr(a)
    end do
    !44   - Right-side of domain
    Xe(44)%xc = 1.0_dp
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
    integer               :: n
    real(dp)              :: xo, dx
    real(dp)              :: al, ar, als, ars, a
    type(Euler1D_W_State) :: Wl, Wr, Wls, Wrs, Wint
    !Allocate exact solution arrays
    ne = 27
    allocate(Xe(ne))
    allocate(We(ne))
    do n = 1, ne
      Xe(n)%xc = 0.0_dp ; Xe(n)%dx = 0.0_dp
      Xe(n)%xa = 0.0_dp ; Xe(n)%da = 0.0_dp
      call Vacuum_W(We(n))
    end do
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
      Xe(1)%xc = 0.0_dp
      Xe(2:26)%xc = xo
      Xe(27)%xc = 1.0_dp
      We(1:25) = Wl
      We(26:27) = Wr
      return
    end if
    !Solve shocktube/Riemann problem
    call Riemann(Wl,Wr,Wls,Wrs,Wint,0)
    als = a_W(Wls)
    ars = a_W(Wrs)
    !1    - Left-side of domain
    Xe(1)%xc = 0.0_dp
    We(1) = Wl
    !2-22 - Head-to-tail of rarefaction wave
    dx = ((Wls%u - als) - (Wl%u - al))*t/20.0_dp
    do n = 2, 22
      Xe(n)%xc = xo + (Wl%u - al)*t + dx*(n-2)
      We(n)%u = alphai*(Wl%u + 2.0_dp*gm1i*al) + 2.0_dp*gp1i*(Xe(n)%xc-xo)/t
      a = We(n)%u - (Xe(n)%xc-xo)/t
      We(n)%p = Wl%p*((a/al)**betai)
      We(n)%rho = g*We(n)%p/sqr(a)
    end do
    !23   - Contact surface (rarefaction side)
    Xe(23)%xc = xo + Wls%u*t
    We(23) = Wls
    !24   - Contact surface (shock side)
    Xe(24)%xc = Xe(23)%xc
    We(24) = Wrs
    !25   - Shock wave (contact surface side)
    Xe(25)%xc = xo + (Wrs%u + ars*sqrt(beta*alpha*Wr%p/Wrs%p + beta))*t
    We(25) = Wrs
    !26   - Shock wave (silent side)
    Xe(26)%xc = Xe(25)%xc
    We(26) = Wr
    !27   - Right-side of domain
    Xe(27)%xc = 1.0_dp
    We(27) = Wr
    return
  end subroutine exactStrongSod

  !---------------------------------------------------------------------
  !Stationary Contact (w/ left rarefaction and right shock)
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
    integer               :: n
    real(dp)              :: xo, dx
    real(dp)              :: al, ar, als, ars, a
    type(Euler1D_W_State) :: Wl, Wr, Wls, Wrs, Wint
    !Allocate exact solution arrays
    ne = 27
    allocate(Xe(ne))
    allocate(We(ne))
    do n = 1, ne
      Xe(n)%xc = 0.0_dp ; Xe(n)%dx = 0.0_dp
      Xe(n)%xa = 0.0_dp ; Xe(n)%da = 0.0_dp
      call Vacuum_W(We(n))
    end do
    !Initial location of jump conditions
    xo = 0.8_dp
    !Initialize the left state
    Wl = WState(1.0_dp,-19.59745_dp,1000.0_dp)
    al = a_W(Wl)
    !Initialize the right state
    Wr = WState(1.0_dp,-19.59745_dp,0.01_dp)
    ar = a_W(Wr)
    !If initial conditions
    if(t.lt.TOLER) then
      Xe(1)%xc = 0.0_dp
      Xe(2:26)%xc = xo
      Xe(27)%xc = 1.0_dp
      We(1:25) = Wl
      We(26:27) = Wr
      return
    end if
    !Set and solve shocktube/Riemann problem
    call Riemann(Wl,Wr,Wls,Wrs,Wint,0)
    als = a_W(Wls)
    ars = a_W(Wrs)
    !1    - Left-side of domain
    Xe(1)%xc = 0.0_dp
    We(1) = Wl
    !2-22 - Head-to-tail of rarefaction wave
    dx = ((Wls%u - als) - (Wl%u - al))*t/20.0_dp
    do n = 2, 22
      Xe(n)%xc = xo + (Wl%u - al)*t + dx*(n-2)
      We(n)%u = alphai*(Wl%u + 2.0_dp*gm1i*al) + 2.0_dp*gp1i*(Xe(n)%xc-xo)/t
      a = We(n)%u - (Xe(n)%xc-xo)/t
      We(n)%p = Wl%p*((a/al)**betai)
      We(n)%rho = g*We(n)%p/sqr(a)
    end do
    !23   - Contact surface (rarefaction side)
    Xe(23)%xc = xo + Wls%u*t
    We(23) = Wls
    !24   - Contact surface (shock side)
    Xe(24)%xc = xe(23)%xc
    We(24) = Wrs
    !25   - Shock wave (contact surface side)
    Xe(25)%xc = xo + (Wrs%u + ars*sqrt(beta*alpha*Wr%p/Wrs%p + beta))*t
    We(25) = Wrs
    !26   - Shock wave (silent side)
    Xe(26)%xc = Xe(25)%xc
    We(26) = Wr
    !27   - Right-side of domain
    Xe(27)%xc = 1.0_dp
    We(27) = Wr
    return
  end subroutine exactStationaryContact

  !---------------------------------------------------------------------
  !Three Right Discontinuities Problem
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
    integer               :: n
    real(dp)              :: xo, dx
    real(dp)              :: al, ar, als, ars, a
    type(Euler1D_W_State) :: Wl, Wr, Wls, Wrs, Wint
    !Allocate exact solution arrays
    ne = 8
    allocate(Xe(ne))
    allocate(We(ne))
    do n = 1, ne
      Xe(n)%xc = 0.0_dp ; Xe(n)%dx = 0.0_dp
      Xe(n)%xa = 0.0_dp ; Xe(n)%da = 0.0_dp
      call Vacuum_W(We(n))
    end do
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
      Xe(1)%xc = 0.0_dp
      Xe(2:7)%xc = xo
      Xe(8)%xc = 1.0_dp
      We(1:6) = Wl
      We(7:8) = Wr
      return
    end if
    !Solve shocktube/Riemann problem
    call Riemann(Wl,Wr,Wls,Wrs,Wint,0)
    als = a_W(Wls)
    ars = a_W(Wrs)
    !1 - Left-side of domain
    Xe(1)%xc = 0.0_dp
    We(1) = Wl
    !2 - Before first shock-wave
    Xe(2)%xc = xo + t*(Wl%u - al*sqrt(beta*alpha*Wls%p/Wl%p + beta))
    We(2) = Wl
    !3 - After first shock-wave
    Xe(3)%xc = xo + t*(Wl%u - al*sqrt(beta*alpha*Wls%p/Wl%p + beta))
    We(3) = Wls
    !4 - Before contact-wave
    Xe(4)%xc = xo + t*Wls%u
    We(4) = Wls
    !5 - After contact-wave
    Xe(5)%xc = xo + t*Wls%u
    We(5) = Wrs
    !6 - Before last shock-wave
    Xe(6)%xc = xo + t*(Wrs%u + ars*sqrt(beta*alpha*Wr%p/Wrs%p + beta))
    We(6) = Wrs
    !7 - After last shock-wave
    Xe(7)%xc = xo + t*(Wrs%u + ars*sqrt(beta*alpha*Wr%p/Wrs%p + beta))
    We(7) = Wr
    !8 - Right-side of domain
    Xe(8)%xc = 1.0_dp
    We(8) = Wr
    return
  end subroutine exact3RightWaves

  !---------------------------------------------------------------------
  !Subsonic or Transonic Nozzle
  !Subsonic: (Sthroat,xShock) = (0.8,11.0)
  !Transonic: (Sthroat,xShock) = (1.0,7.0)
  !---------------------------------------------------------------------
  recursive subroutine Nozzle(t)
    use realsizes, only: dp
    use numbers, only: TOLER
    use inputParams, only: i_problem
    use cfdParams, only: SUBSONIC_NOZZLE, TRANSONIC_NOZZLE
    use gasConstants, only: g, gm1, gp1, gm1i, gp1i, R, alpha, alphai, betai
    use Euler1D_WState
    implicit none
    !Argument variables:
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
!!$    allocate(Xe(ne))
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
  !---------------------------------------------------------------------
  subroutine exactSquareWave(t)
    use realsizes, only: dp
    use inputParams, only: rhom, um, pm, wm
    use Euler1D_WState
    implicit none
    !Argument variables:
    real(dp), intent(in)  :: t
    !Local variables:
    integer               :: n
    real(dp)              :: xl, xr
    type(Euler1D_W_State) :: Wo, Ww
    !Initialize the wave states
    Wo%rho = 1.0_dp*rhom ; Wo%u = um ; Wo%p = pm
    Ww%rho = 2.0_dp*rhom ; Ww%u = um ; Ww%p = pm
    !Get wave limit locations:
    xl = 0.0_dp - 0.5_dp*wm + t*um
    xr = 0.0_dp + 0.5_dp*wm + t*um
    if(um.gt.0.0_dp) then
      do while(xl.gt.0.5_dp)
        xl = xl - 1.0_dp
      end do
      do while(xr.gt.0.5_dp)
        xr = xr - 1.0_dp
      end do
    else
      do while(xl.lt.-0.5_dp)
        xl = xl + 1.0_dp
      end do
      do while(xr.lt.-0.5_dp)
        xr = xr + 1.0_dp
      end do
    end if
    !Allocate exact solution arrays
    ne = 6
    allocate(Xe(ne))
    allocate(We(ne))
    do n = 1, ne
      Xe(n)%xc = 0.0_dp ; Xe(n)%dx = 0.0_dp
      Xe(n)%xa = 0.0_dp ; Xe(n)%da = 0.0_dp
      call Vacuum_W(We(n))
    end do
    !Calculate the primitive solution state based on the location and time
    if(xl.lt.xr) then
      Xe(1)%xc = -0.5_dp ; We(1) = Wo
      Xe(2)%xc = xl      ; We(2) = Wo
      Xe(3)%xc = xl      ; We(3) = Ww
      Xe(4)%xc = xr      ; We(4) = Ww
      Xe(5)%xc = xr      ; We(5) = Wo
      Xe(6)%xc =  0.5_dp ; We(6) = Wo
    else
      Xe(1)%xc = -0.5_dp ; We(1) = Ww
      Xe(2)%xc = xr      ; We(2) = Ww
      Xe(3)%xc = xr      ; We(3) = Wo
      Xe(4)%xc = xl      ; We(4) = Wo
      Xe(5)%xc = xl      ; We(5) = Ww
      Xe(6)%xc =  0.5_dp ; We(6) = Ww
    end if
    return
  end subroutine exactSquareWave

  !---------------------------------------------------------------------
  !Sine-Squared Density Wave
  !---------------------------------------------------------------------
  subroutine exactSineSquaredWave(t)
    use realsizes, only: dp
    use inputParams, only: rhom, um, pm, wm
    use numbers, only: PI
    use Euler1D_WState
    implicit none
    !Argument variables:
    real(dp), intent(in)  :: t
    !Local variables:
    integer               :: n
    real(dp)              :: xl, xr
    type(Euler1D_W_State) :: Wo
    !Initialize the wave states
    Wo%rho = rhom ; Wo%u = um ; Wo%p = pm
    !Get wave limit locations:
    xl = 0.0_dp - 0.5_dp*wm + t*um
    xr = 0.0_dp + 0.5_dp*wm + t*um
    if(um.gt.0.0_dp) then
      do while(xl.gt.0.5_dp)
        xl = xl - 1.0_dp
      end do
      do while(xr.gt.0.5_dp)
        xr = xr - 1.0_dp
      end do
    else
      do while(xl.lt.-0.5_dp)
        xl = xl + 1.0_dp
      end do
      do while(xr.lt.-0.5_dp)
        xr = xr + 1.0_dp
      end do
    end if
    ne = 53
    allocate(Xe(ne))
    allocate(We(ne))
    do n = 1, ne
      Xe(n)%xc = 0.0_dp ; Xe(n)%dx = 0.0_dp
      Xe(n)%xa = 0.0_dp ; Xe(n)%da = 0.0_dp
      call Vacuum_W(We(n))
    end do
    !Calculate the primitive solution state based on the location and time
    if(xl.lt.xr) then
      Xe( 1)%xc = -0.5_dp ; We( 1) = Wo
      do n = 2, 52
        We(n) = Wo
        Xe(n)%xc = xl + (n-2)*(xr-xl)/50
        We(n)%rho = Wo%rho*(1.0_dp + 0.5_dp*(cos(PI*(-0.5_dp + (n-2.0_dp)/50.0_dp))**2))
      end do
      Xe(53)%xc =  0.5_dp ; We(53) = Wo
    else
    end if
    return
  end subroutine exactSineSquaredWave

end module exactSoln_module
