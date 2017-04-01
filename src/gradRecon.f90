!---------------------------------------------------------------------
! Determine the second-order central difference approximation to the
! gradient.
!---------------------------------------------------------------------
subroutine GGreconstruction
  use realSizes, only: dp
  use solnBlock_module, only: W, dWdx, NCl, NCu
  use gridBlock_module, only: Cell, xfaceL, xfaceR
  use Euler1D_WState
  implicit none
  integer               :: nc
  real(dp)              :: dx
  type(Euler1D_W_State) :: dWrl
  do nc = NCl-1, NCu+1
   dx   = Cell(nc+1)%Xc - Cell(nc-1)%Xc
   dWrl = W(nc+1) - W(nc-1)
   dWdx(nc) = dWrl/dx
  end do
  return
end subroutine GGreconstruction


!---------------------------------------------------------------------
! Determine the least squares approximation to the gradient.
!---------------------------------------------------------------------
subroutine LSQreconstruction
  use realSizes, only: dp
  use mathFunctions, only: sqr
  use solnBlock_module, only: W, dWdx, phi, NCl, NCu, Ng
  use gridBlock_module, only: Cell, xfaceL, xfaceR
  use Euler1D_WState
  implicit none
  integer  :: nc
  real(dp) :: dxm, dxp, dxl, dxr
  do nc = NCl-1, NCu+1
   dxm = Cell(nc-1)%Xc - Cell(nc)%Xc
   dxp = Cell(nc+1)%Xc - Cell(nc)%Xc
   dxl = dxm/(sqr(dxp) + sqr(dxm))
   dxr = dxp/(sqr(dxp) + sqr(dxm))
   dWdx(nc) = dxr*(W(nc+1)-W(nc)) + dxl*(W(nc-1)-W(nc))
  end do
  return
end subroutine LSQreconstruction


!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine WENOstencil(v1,v2,v3,v4,v5,v)
  use realSizes, only: dp
  use numbers, only: MICRO
  real(dp), intent(in)  :: v1, v2, v3, v4, v5
  real(dp), intent(out) :: v
  real(dp) :: p1, p2, p3        !ENO stencil approximations
  real(dp) :: s1, s2, s3        !Smoothness estimates
  real(dp) :: a1, a2, a3, atot  !Alpha constants
  real(dp) :: w1, w2, w3        !Stencil weightings
  !ENO approximations:
  p1 = (2.0_dp*v1 - 7.0_dp*v2 + 11.0_dp*v3)/6.0_dp
  p2 = (      -v2 + 5.0_dp*v3 +  2.0_dp*v4)/6.0_dp
  p3 = (2.0_dp*v3 + 5.0_dp*v4 -         v5)/6.0_dp
  !Smoothness Estimates:
  s1 = (13.0_dp/12.0_dp)*(v1-2.0_dp*v2+v3)**2 + 0.25_dp*(v1-4.0_dp*v2+3.0_dp*v3)**2
  s2 = (13.0_dp/12.0_dp)*(v2-2.0_dp*v3+v4)**2 + 0.25_dp*(v2-v4)**2
  s3 = (13.0_dp/12.0_dp)*(v3-2.0_dp*v4+v5)**2 + 0.25_dp*(3.0_dp*v3-4.0_dp*v4+v5)**2
  !Alpha constants:
  a1 = 0.1_dp/((s1 + MICRO)**2)
  a2 = 0.6_dp/((s2 + MICRO)**2)
  a3 = 0.3_dp/((s3 + MICRO)**2)
  atot = at + a2 + a3
  !Stencil weightings:
  w1 = a1/atot
  w2 = a2/atot
  w3 = a3/atot
  !Determine the WENO gradient as a weighted combination of the ENO
  !approximations.
  v = w1*p1 + w2*p2 + w3*p3
  return
end subroutine WENOstencil


!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine WENOreconstruction
  use realSizes, only: dp
  use gridBlock_module, only: Cell
  use solnBlock_module, only: W, dWdx, NCl, NCu
  use Euler1D_WState
  integer  :: nc
  integer  :: k
  real(dp) :: v1, v2, v3, v4, v5
  real(dp) :: dxl, dxr
  real(dp) :: wmin, wmax, wl, wm, wr, dw, p
  real(dp), external :: xfracel, xfacer

  do nc = NCl-1, NCu+1

    if(W(nc)%u.lt.0.0_dp) then
      !Compute right-biased density stencil.
      v1 = (W(nc+3)%rho-W(nc+2)%rho)/(Cell(nc+3)%Xc - Cell(nc+2)%Xc)
      v2 = (W(nc+2)%rho-W(nc+1)%rho)/(Cell(nc+2)%Xc - Cell(nc+1)%Xc)
      v3 = (W(nc+1)%rho-W(nc  )%rho)/(Cell(nc+1)%Xc - Cell(nc  )%Xc)
      v4 = (W(nc  )%rho-W(nc-1)%rho)/(Cell(nc  )%Xc - Cell(nc-1)%Xc)
      v5 = (W(nc-1)%rho-W(nc-2)%rho)/(Cell(nc-1)%Xc - Cell(nc-2)%Xc)
      call WENOstencil(v1,v2,v3,v4,v5,dWdx(nc)%rho)
      !Compute right-biased velocity stencil.
      v1 = (W(nc+3)%u-W(nc+2)%u)/(Cell(nc+3)%Xc - Cell(nc+2)%Xc)
      v2 = (W(nc+2)%u-W(nc+1)%u)/(Cell(nc+2)%Xc - Cell(nc+1)%Xc)
      v3 = (W(nc+1)%u-W(nc  )%u)/(Cell(nc+1)%Xc - Cell(nc  )%Xc)
      v4 = (W(nc  )%u-W(nc-1)%u)/(Cell(nc  )%Xc - Cell(nc-1)%Xc)
      v5 = (W(nc-1)%u-W(nc-2)%u)/(Cell(nc-1)%Xc - Cell(nc-2)%Xc)
      call WENOstencil(v1,v2,v3,v4,v5,dWdx(nc)%u)
      !Compute right-biased pressure stencil.
      v1 = (W(nc+3)%p-W(nc+2)%p)/(Cell(nc+3)%Xc - Cell(nc+2)%Xc)
      v2 = (W(nc+2)%p-W(nc+1)%p)/(Cell(nc+2)%Xc - Cell(nc+1)%Xc)
      v3 = (W(nc+1)%p-W(nc  )%p)/(Cell(nc+1)%Xc - Cell(nc  )%Xc)
      v4 = (W(nc  )%p-W(nc-1)%p)/(Cell(nc  )%Xc - Cell(nc-1)%Xc)
      v5 = (W(nc-1)%p-W(nc-2)%p)/(Cell(nc-1)%Xc - Cell(nc-2)%Xc)
      call WENOstencil(v1,v2,v3,v4,v5,dWdx(nc)%p)

    else !if(W(nc)%u.ge.0.0_dp) then
      !Compute left-biased density stencil.
      v1 = (W(nc-2)%rho-W(nc-3)%rho)/(Cell(nc-2)%Xc - Cell(nc-3)%Xc)
      v2 = (W(nc-1)%rho-W(nc-2)%rho)/(Cell(nc-1)%Xc - Cell(nc-2)%Xc)
      v3 = (W(nc  )%rho-W(nc-1)%rho)/(Cell(nc  )%Xc - Cell(nc-1)%Xc)
      v4 = (W(nc+1)%rho-W(nc  )%rho)/(Cell(nc+1)%Xc - Cell(nc  )%Xc)
      v5 = (W(nc+2)%rho-W(nc+1)%rho)/(Cell(nc+2)%Xc - Cell(nc+1)%Xc)
      call WENOstencil(v1,v2,v3,v4,v5,dWdx(nc)%rho)
      ! Compute left-biased velocity stencil.
      v1 = (W(nc-2)%u-W(nc-3)%u)/(Cell(nc-2)%Xc - Cell(nc-3)%Xc)
      v2 = (W(nc-1)%u-W(nc-2)%u)/(Cell(nc-1)%Xc - Cell(nc-2)%Xc)
      v3 = (W(nc  )%u-W(nc-1)%u)/(Cell(nc  )%Xc - Cell(nc-1)%Xc)
      v4 = (W(nc+1)%u-W(nc  )%u)/(Cell(nc+1)%Xc - Cell(nc  )%Xc)
      v5 = (W(nc+2)%u-W(nc+1)%u)/(Cell(nc+2)%Xc - Cell(nc+1)%Xc)
      call WENOstencil(v1,v2,v3,v4,v5,dWdx(nc)%u)
      ! Compute left-biased pressure stencil.
      v1 = (W(nc-2)%p-W(nc-3)%p)/(Cell(nc-2)%Xc - Cell(nc-3)%Xc)
      v2 = (W(nc-1)%p-W(nc-2)%p)/(Cell(nc-1)%Xc - Cell(nc-2)%Xc)
      v3 = (W(nc  )%p-W(nc-1)%p)/(Cell(nc  )%Xc - Cell(nc-1)%Xc)
      v4 = (W(nc+1)%p-W(nc  )%p)/(Cell(nc+1)%Xc - Cell(nc  )%Xc)
      v5 = (W(nc+2)%p-W(nc+1)%p)/(Cell(nc+2)%Xc - Cell(nc+1)%Xc)
      call WENOstencil(v1,v2,v3,v4,v5,dWdx(nc)%p)
    end if

  end do

  return
end subroutine WENOreconstruction


!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine WENOreconstruction2
  use realSizes, only: dp
  use solnBlock_module, only: W, Ql, Qr
  use Euler1D_WState
  integer  :: nc
  integer  :: k
  real(dp) :: v1, v2, v3, v4, v5

  do nc = NCl-1, NCu+1

    !Compute the right state for the current cell.
    if(nc.ge.NCl-1) then
      !Compute right-biased density stencil.
      v1 = W(nc-2)%rho
      v2 = W(nc-1)%rho
      v3 = W(nc  )%rho
      v4 = W(nc+1)%rho
      v5 = W(nc+2)%rho
      call WENOstencil(v1,v2,v3,v4,v5,Qr(nc)%rho)
      !Compute right-biased speed stencil.
      v1 = W(nc-2)%u
      v2 = W(nc-1)%u
      v3 = W(nc  )%u
      v4 = W(nc+1)%u
      v5 = W(nc+2)%u
      call WENOstencil(v1,v2,v3,v4,v5,Qr(nc)%u)
      !Compute right-biased pressure stencil.
      v1 = W(nc-2)%p
      v2 = W(nc-1)%p
      v3 = W(nc  )%p
      v4 = W(nc+1)%p
      v5 = W(nc+2)%p
      call WENOstencil(v1,v2,v3,v4,v5,Qr(nc)%p)
    end if

    !Compute the left state for the current cell.
    if(nc.le.NCu+1) then
      !Compute left-biased density stencil.
      v1 = W(nc+2)%rho
      v2 = W(nc+1)%rho
      v3 = W(nc  )%rho
      v4 = W(nc-1)%rho
      v5 = W(nc-2)%rho
      call WENOstencil(v1,v2,v3,v4,v5,Ql(nc)%rho)
      !Compute left-biased speed stencil.
      v1 = W(nc+2)%u
      v2 = W(nc+1)%u
      v3 = W(nc  )%u
      v4 = W(nc-1)%u
      v5 = W(nc-2)%u
      call WENOstencil(v1,v2,v3,v4,v5,Ql(nc)%u)
      ! Compute left-biased pressure stencil.
      v1 = W(nc+2)%p
      v2 = W(nc+1)%p
      v3 = W(nc  )%p
      v4 = W(nc-1)%p
      v5 = W(nc-2)%p
      call WENOstencil(v1,v2,v3,v4,v5,Ql(nc)%p)
    end if

  end do
  
  return
end subroutine WENOreconstruction2


!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine PPMstencil(vm,v,vp,vpp,vl,vr)
  use realSizes, only: dp
  implicit none
  real(dp), intent(in)  :: vm, v, vp, vpp
  real(dp), intent(out) :: vl, vr
  real(dp) :: dv, dvm, dvp, S, Sp
  dv = vp - v
  dvm = v - vm
  if(dv*dvm.gt.0.0_dp) then
    S = min(0.5_dp*abs(dv+dvm),2.0_dp*min(abs(dv),abs(dvm)))
    if(0.5_dp*(dv+dvm).lt.0.0_dp) then
      S = -S
    end if
  else
    S = 0.0_dp
  end if
  dvp = vpp - vp
  if(dvp*dv.gt.0.0_dp) then
    Sp = min(0.5_dp*abs(dvp+dv),2.0_dp*min(abs(dvp),abs(dv)))
    if(0.5_dp*(dvp+dv) .lt. 0.0_dp) then
      Sp = -Sp
    end if
  else
    Sp = 0.
  end if
  vl = 0.5_dp*(v+vp) + (S+Sp)/6.0_dp
  vr = vl
  return
end subroutine PPMstencil

subroutine PPMlimiting(v,vr,vl)
  use realSizes, only: dp
  use numbers, only: NANO
  implicit none
  real(dp), intent(in)    :: v
  real(dp), intent(inout) :: vr,vl
  real(dp) :: C, D
  if((vl-v)/(v-vr+NANO) .lt. 0.0_dp) then
    vl = v ; vr = v
  else
    C = vl - vr
    D = 6.0_dp*(v - 0.5_dp*(vl+vr))
    if(D*C .gt. C*C) then
      vr = 3.0_dp*v - 2.0_dp*vl
    else if(-C*C .gt. D*C) then
      vl = 3.0_dp*v - 2.0_dp*vr
    end if
  end if
  return
end subroutine PPMlimiting

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine PPMreconstruction
  use realSizes, only: dp
  use solnBlock_module, only: W, Ql, Qr, NCl, NCu
  use Euler1D_WState
  implicit none
  integer :: nc

  !Compute the left and right states at cell interface i+1/2.
  do nc = NCl-2, NCu+1
    !Compute density stencil.
    call PPMstencil(W(nc-1)%rho,W(nc)%rho,W(nc+1)%rho,W(nc+2)%rho,Qr(nc)%rho,Ql(nc+1)%rho)
    !Compute speed stencil.
    call PPMstencil(W(nc-1)%u,W(nc)%u,W(nc+1)%u,W(nc+2)%u,Qr(nc)%u,Ql(nc+1)%u)
    !Compute pressure stencil.
    call PPMstencil(W(nc-1)%p,W(nc)%p,W(nc+1)%p,W(nc+2)%p,Qr(nc)%p,Ql(nc+1)%p)
  end do

  !Limit the left and right states based on the right state at i-1/2
  !and the left state at i+1/2 (all current cell).
  do nc = NCl-1, NCu+1
    call PPMlimiting(W(nc)%rho,Ql(nc)%rho,Qr(nc)%rho)
    call PPMlimiting(W(nc)%u,Ql(nc)%u,Qr(nc)%u)
    call PPMlimiting(W(nc)%p,Ql(nc)%p,Qr(nc)%p)
  end do
  
  return
end subroutine PPMreconstruction


!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine Gamma_Differencing_Scalar(phiL,phiC,phiR,dphidx,xl,xc,xr,v,phi)
  use realSizes, only: dp
  use inputParams, only: betam
  implicit none
  real(dp), intent(in) :: phiL, phiC, phiR, dphidx, xl, xc, xr, v
  real(dp), intent(inout) :: phi
  real(dp) :: phiU, phiD, xu, xd
  real(dp) :: phitilde, d, fx, g
  if(v.ge.0.) then
    phiU = phiL ; xu = xl ; phiD = phiR ; xd = xr
  else
    phiU = phiR ; xu = xr ; phiD = phiL ; xd = xl
  end if
  d = xd-xc
  fx = 0.5_dp*d/d
  phitilde = 1.0_dp - 0.5_dp*(phiD-phiC)/(dphidx*d)
  if(phitilde.le.0.0_dp) then
    phi = phiC
  else if(phitilde.lt.betam) then
    g = phitilde/betam
    phi = (1.0_dp - g*(1.0_dp-fx))*phiC + g*(1.0_dp-fx)*phiD
  else if(phitilde.lt.1.0_dp) then
    phi = fx*phiC + (1.0_dp-fx)*phiD
  else
    phi = phiC
  end if
  return
end subroutine Gamma_Differencing_Scalar

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine Gamma_Differencing_Scheme
  use realSizes, only: dp
  use gridBlock_module, only: Cell
  use solnBlock_module, only: W, dWdx, Ql, Qr, NCl, NCu
  use Euler1D_WState
  implicit none
  integer :: nc
  real(dp) :: v
  type(Euler1D_W_State) :: phi, dphi, dW
  do nc = NCl-1, NCu+1
    if(nc.ge.NCl) then
      v = 0.5_dp*(W(nc)%u+W(nc-1)%u)
      call Gamma_Differencing_Scalar(W(nc+1)%rho,W(nc)%rho,W(nc-1)%rho,dWdx(nc)%rho, &
           Cell(nc+1)%Xc,Cell(nc)%Xc,Cell(nc-1)%Xc,v,Ql(nc)%rho)
      call Gamma_Differencing_Scalar(W(nc+1)%u,W(nc)%u,W(nc-1)%u,dWdx(nc)%u, &
           Cell(nc+1)%Xc,Cell(nc)%Xc,Cell(nc-1)%Xc,v,Ql(nc)%u)
      call Gamma_Differencing_Scalar(W(nc+1)%p,W(nc)%p,W(nc-1)%p,dWdx(nc)%p, &
           Cell(nc+1)%Xc,Cell(nc)%Xc,Cell(nc-1)%Xc,v,Ql(nc)%p)
    end if
    if(nc.le.NCu) then
      v = 0.5_dp*(W(nc)%u+W(nc+1)%u)
      call Gamma_Differencing_Scalar(W(nc-1)%rho,W(nc)%rho,W(nc+1)%rho,dWdx(nc)%rho, &
           Cell(nc-1)%Xc,Cell(nc)%Xc,Cell(nc+1)%Xc,v,Qr(nc)%rho)
      call Gamma_Differencing_Scalar(W(nc-1)%u,W(nc)%u,W(nc+1)%u,dWdx(nc)%u, &
           Cell(nc-1)%Xc,Cell(nc)%Xc,Cell(nc+1)%Xc,v,Qr(nc)%u)
      call Gamma_Differencing_Scalar(W(nc-1)%p,W(nc)%p,W(nc+1)%p,dWdx(nc)%p, &
           Cell(nc-1)%Xc,Cell(nc)%Xc,Cell(nc+1)%Xc,v,Qr(nc)%p)
    end if
  end do
  return
end subroutine Gamma_Differencing_Scheme


!---------------------------------------------------------------------
! Determine the gradient limiters for each variable.
!---------------------------------------------------------------------
subroutine calcLimiters
  use realSizes, only: dp
  use cfdParams
  use inputParams, only: i_limiter
  use solnBlock_module, only: W, dWdx, phi, NCl, NCu
  use gridBlock_module, only: Cell, xfaceL, xfaceR
  use Euler1D_WState
  implicit none
  integer               :: nc
  real(dp)              :: dxl, dxr
  type(Euler1D_W_State) :: Wmin, Wmax, Wl, Wr

  do nc = NCl-1, NCu+1

   if(i_limiter.eq.LIMITER_ZERO) then
     call constant_W(phi(nc),0.0_dp)

   else if(i_limiter.eq.LIMITER_ONE) then
     call constant_W(phi(nc),1.0_dp)

   else
     dxl = xfaceL(nc) - Cell(nc)%Xc
     dxr = xfaceR(nc) - Cell(nc)%Xc
     Wmin = min_W(W(nc),min_W(W(nc-1),W(nc+1)))
     Wmax = max_W(W(nc),max_W(W(nc-1),W(nc+1)))
     Wl = W(nc) + dxl*dWdx(nc)
     Wr = W(nc) + dxr*dWdx(nc)
     select case(i_limiter)
     case(LIMITER_VANLEER)
       phi(nc)%rho = VanLeer(Wl%rho,W(nc)%rho,Wr%rho,Wmin%rho,Wmax%rho)
       phi(nc)%u   = VanLeer(Wl%u,W(nc)%u,Wr%u,Wmin%u,Wmax%u)
       phi(nc)%p   = VanLeer(Wl%p,W(nc)%p,Wr%p,Wmin%p,Wmax%p)
     case(LIMITER_VANALBADA)
       phi(nc)%rho = vanalbada(Wl%rho,W(nc)%rho,Wr%rho,Wmin%rho,Wmax%rho)
       phi(nc)%u   = vanalbada(Wl%u,W(nc)%u,Wr%u,Wmin%u,Wmax%u)
       phi(nc)%p   = vanalbada(Wl%p,W(nc)%p,Wr%p,Wmin%p,Wmax%p)
     case(LIMITER_BARTH_JESPERSEN)
       phi(nc)%rho = BarthJespersen(Wl%rho,W(nc)%rho,Wr%rho,Wmin%rho,Wmax%rho)
       phi(nc)%u   = BarthJespersen(Wl%u,W(nc)%u,Wr%u,Wmin%u,Wmax%u)
       phi(nc)%p   = BarthJespersen(Wl%p,W(nc)%p,Wr%p,Wmin%p,Wmax%p)
     case(LIMITER_VENKATAKRISHNAN)
       phi(nc)%rho = Venkatakrishnan(Wl%rho,W(nc)%rho,Wr%rho,Wmin%rho,Wmax%rho)
       phi(nc)%u   = Venkatakrishnan(Wl%u,W(nc)%u,Wr%u,Wmin%u,Wmax%u)
       phi(nc)%p   = Venkatakrishnan(Wl%p,W(nc)%p,Wr%p,Wmin%p,Wmax%p)
     end select
   end if

  end do

  return
end subroutine calcLimiters


