!-----------------------------------------------------------------------
!Subroutines for computing fluxes.
!-----------------------------------------------------------------------

!---------------------------------------------------------------------
!Transform a primitive (rho,u,p) state to a conserved variable state.
!---------------------------------------------------------------------
subroutine transform_W_to_U(W, U)
  use Euler1D_UState
  use Euler1D_WState
  implicit none
  type(Euler1D_W_State), intent(in)  :: W
  type(Euler1D_U_State), intent(out) :: U
  U%rho = W%rho
  U%du = du_W(W)
  U%E = E_W(W)
  return
end subroutine transform_W_to_U

!---------------------------------------------------------------------
!Transform a conserved variable state to a primitive (rho,u,p) state.
!---------------------------------------------------------------------
subroutine transform_U_to_W(U, W)
  use Euler1D_UState
  use Euler1D_WState
  implicit none
  type(Euler1D_U_State), intent(in)  :: U
  type(Euler1D_W_State), intent(out) :: W
  W%rho = U%rho
  W%u = u_U(U)
  W%p = p_U(U)
  return
end subroutine transform_U_to_W

!---------------------------------------------------------------------
!Determine the Roe-Average state from the two input primitive states.
!---------------------------------------------------------------------
subroutine roe_average(Wl, Wr, Wa)
  use realSizes, only: dp
  use gasConstants, only: g, gm1
  use Euler1D_WState
  implicit none
  type(Euler1D_W_State), intent(in)  :: Wl, Wr
  type(Euler1D_W_State), intent(out) :: Wa
  real(dp)                     :: hl, hr, srhol, srhor, ha
  !Determine the left and right state specific enthalpies and square
  !roots of the density.
  hl = Hsp_W(Wl)
  hr = Hsp_W(Wr)
  srhol = sqrt(Wl%rho)
  srhor = sqrt(Wr%rho)
  !Determine the appropriate Roe averages.
  Wa%rho = srhol*srhor
  Wa%u = (srhol*Wl%u + srhor*Wr%u)/(srhol+srhor)
  ha = (srhol*hl + srhor*hr)/(srhol+srhor)
  Wa%p = Wa%rho*gm1*(ha - 0.5*Wa%u*Wa%u)/g
  return
end subroutine roe_average

!---------------------------------------------------------------------
!Compute the positive parts of the corrected elemental wave speeds
!(eigenvalues) according to Harten's entropy fix (1983).
!---------------------------------------------------------------------
subroutine HartenFixPos(lambda_a, lambda_l, lambda_r, lambda)
  use realSizes, only: dp
  implicit none
  real(dp), intent(in)  :: lambda_a, lambda_l, lambda_r
  real(dp), intent(out) :: lambda
  real(dp)              :: delta
  delta = max(0.,2.*(lambda_r-lambda_l))
  if(abs(lambda_a).gt.delta.or.delta.lt.0.0000001) then
    delta = 0.
  else
    delta = 0.5*(delta + lambda_a*lambda_a/delta) - abs(lambda_a)
  end if
  lambda = 0.5*(lambda_a + abs(lambda_a) + delta)
  return
end subroutine HartenFixPos

!---------------------------------------------------------------------
!Compute the Harten positive-wave entropy fix.
!---------------------------------------------------------------------
subroutine harten_fix_pos(lambda_a, lambda_l, lambda_r, lambda)
  use Euler1D_WState
  implicit none
  type(Euler1D_W_State), intent(in)  :: lambda_a, lambda_l, lambda_r
  type(Euler1D_W_State), intent(out) :: lambda
  call HartenFixPos(lambda_a%rho,lambda_l%rho,lambda_r%rho,lambda%rho)
  lambda%u = 0.5*(lambda_a%u+abs(lambda_a%u))
  call HartenFixPos(lambda_a%p,lambda_l%p,lambda_r%p,lambda%p)
  return
end subroutine harten_fix_pos

!---------------------------------------------------------------------
!Compute the negative parts of the corrected elemental wave speeds
!(eigenvalues) according to Harten's entropy fix (1983).
!---------------------------------------------------------------------
subroutine HartenFixNeg(lambda_a, lambda_l, lambda_r, lambda)
  use realSizes, only: dp
  implicit none
  real(dp), intent(in)  :: lambda_a, lambda_l, lambda_r
  real(dp), intent(out) :: lambda
  real(dp)              :: delta
  delta = max(0.,2.*(lambda_r-lambda_l))
  if(abs(lambda_a).gt.delta.or.delta.lt.0.0000001) then
    delta = 0.
  else
    delta = 0.5*(delta + lambda_a*lambda_a/delta) - abs(lambda_a)
  end if
  lambda = 0.5*(lambda_a - abs(lambda_a) - delta)
  return
end subroutine HartenFixNeg

!---------------------------------------------------------------------
!Compute the Harten negative-wave entropy fix.
!---------------------------------------------------------------------
subroutine harten_fix_neg(lambda_a, lambda_l, lambda_r, lambda)
  use Euler1D_WState
  implicit none
  type(Euler1D_W_State), intent(in)  :: lambda_a, lambda_l, lambda_r
  type(Euler1D_W_State), intent(out) :: lambda
  call HartenFixNeg(lambda_a%rho,lambda_l%rho,lambda_r%rho,lambda%rho)
  lambda%u = 0.5*(lambda_a%u-abs(lambda_a%u))
  call HartenFixNeg(lambda_a%p,lambda_l%p,lambda_r%p,lambda%p)
  return
end subroutine harten_fix_neg

!---------------------------------------------------------------------
!Compute the absolute values of the corrected elemental wave speeds
!(eigenvalues) according to Harten's entropy fix (1983).
!---------------------------------------------------------------------
subroutine HartenFixAbs(lambda_a, lambda_l, lambda_r, lambda)
  use realSizes, only: dp
  implicit none
  real(dp), intent(in)  :: lambda_a, lambda_l, lambda_r
  real(dp), intent(out) :: lambda
  real(dp)              :: delta
  delta = max(0.,2.*(lambda_r-lambda_l))
  if(abs(lambda_a).gt.delta.or.delta.lt.0.0000001) then
    delta = 0.
  else
    delta = 0.5*(delta + lambda_a*lambda_a/delta) - abs(lambda_a)
  end if
  lambda = abs(lambda_a) + delta
  return
end subroutine HartenFixAbs

!---------------------------------------------------------------------
!Compute the Harten absolute entropy fix.
!---------------------------------------------------------------------
subroutine harten_fix_abs(lambda_a, lambda_l, lambda_r, lambda)
  use Euler1D_WState
  implicit none
  type(Euler1D_W_State), intent(in)  :: lambda_a, lambda_l, lambda_r
  type(Euler1D_W_State), intent(out) :: lambda
  call HartenFixAbs(lambda_a%rho,lambda_l%rho,lambda_r%rho,lambda%rho)
  lambda%u = abs(lambda_a%u)
  call HartenFixAbs(lambda_a%p,lambda_l%p,lambda_r%p,lambda%p)
  return
end subroutine harten_fix_abs

!---------------------------------------------------------------------
!This subroutine uses a Newton-Raphson interative procedure to obtain
!the exact solution to the Riemann problem for the 1D Euler equations
!in the x-direction, returning the intermediate state variables along
!the ray x/t=0.  See Gottlieb and Groth (1987).
!---------------------------------------------------------------------
subroutine Riemann(Wl,Wr,Wls,Wrs,W,istate)
  use realSizes, only: dp
  use numbers, only: NANO
  use mathFunctions, only: sqr, cube
  use GasConstants, only: g, gm1, gm1i, gp1, gp1i, beta, betai
  use Euler1D_WState
  implicit none
  type(Euler1D_W_State), intent(in)  :: Wl, Wr
  type(Euler1D_W_State), intent(out) :: Wls, Wrs, W
  integer, intent(in)                :: istate
  integer        :: n_iterations
  real(dp) :: al, ar, CL, CR, Z
  real(dp) :: dml, dmr, vm, pm, aml, amr
  real(dp) :: msl, pml, dpmldum, msr, pmr, dpmrdum
  real(dp) :: vsl, vhl, vtl, vsr, vhr, vtr
  logical        :: soln_found

  !Determine the left and right state sound speeds.
  al = a_W(Wl)
  ar = a_W(Wr)

  !Compute the left and right state Riemann invariants.
  CL = Wl%u + 2.*gm1i*al
  CR = Wr%u - 2.*gm1i*ar

  !Check for vacuum state.
  if(CL-CR.le.0.) then
    call vacuum_W(Wls)
    call vacuum_W(Wrs)
    call vacuum_W(W)
    return
  end if

  !Make an initial estimate of the intermediate state flow velocity 
  !to begin the Newton-Raphson iterative solution procedure.  The 
  !initial guess state velocity is made based on isentropic flow
  !theory.
  Z = (ar/al)*((Wl%p/Wr%p)**beta)
  vm = (CL*Z + CR)/(1.0_dp + Z)

  !In the case that two rarefaction waves are present, then an exact
  !solution has been found and the iterative procedure is not
  !required
  soln_found = .false.
  if((vm.ge.Wl%u).and.(vm.le.Wr%u)) then
    if(vm.ge.0.) then
      aml = al - 0.5*gm1*(vm - Wl%u)
      pm = Wl%p*((aml/al)**betai)
      vhl = Wl%u - al
      vtl = vm - aml
      if(vhl.ge.0.) then
        Wls = Wl
      else if(vtl.le.0.) then
        dml = g*pm/sqr(aml)
        Wls = WState(dml,vm,pm)
      else
        vm = (gm1*Wl%u + 2.*al)/gp1
        pm = Wl%p*((vm/al)**beta)
        dml = g*pm/sqr(vm)
        Wls = WState(dml,vm,pm)
      end if
      Wrs = Wls
      soln_found = .true.
    else
      amr = ar + 0.5*gm1*(vm - Wr%u)
      pm = Wr%p*((amr/ar)**betai)
      vhr = Wr%u + ar
      vtr = vm + amr
      if(vhr.le.0.) then
        Wrs = Wr
      else if(vtr.ge.0.) then
        dmr = g*pm/sqr(amr)
        Wrs = WState(dmr,vm,pm)
      else
        vm = (gm1*Wr%u - 2.*ar)/gp1
        pm = Wr%p*((-vm/ar)**betai)
        dmr = g*pm/sqr(vm)
        Wrs = WState(dmr,vm,pm)
      end if
      Wls = Wrs
      soln_found = .true.
    end if
  end if

  !Perform the Newton-Raphson iterative procedure and solve for the
  !velocity in the intermediate state.  During this iterative process
  !the pressure in the intermediate state is also found.
  if(.not.soln_found) then

    do n_iterations = 1, 1000

      !Determine solution changes for left wave.
      if(vm.lt.Wl%u) then
        !Shock wave.
        msl = 0.25*gp1*(vm - Wl%u)/al
        msl = msl - sqrt(1.0_dp + sqr(msl))
        pml = Wl%p*(1.0_dp + g*(vm - Wl%u)*msl/al)
        dpmldum = 2.0_dp*g*Wl%p*cube(msl)/(al*(1.0_dp + sqr(msl)))
        aml = al*sqrt((gp1+gm1*pml/Wl%p)/(gp1+gm1*Wl%p/pml))
      else
        !Rarefaction wave.
        aml = al - 0.5*gm1*(vm - Wl%u)
        pml = Wl%p*((aml/al)**betai)
        dpmldum = -g*pml/aml
      end if

      !Determine solution changes for right wave.
      if(vm.gt.Wr%u) then
        !Shock wave.
        msr = 0.25_dp*gp1*(vm - Wr%u)/ar
        msr = msr + sqrt(1.0_dp + sqr(msr))
        pmr = Wr%p*(1.+g*(vm - Wr%u)*msr/ar)
        dpmrdum = 2.0_dp*g*Wr%p*cube(msr)/(ar*(1.0_dp+sqr(msr)))
        amr = ar*sqrt((gp1+gm1*pmr/Wr%p)/(gp1+gm1*Wr%p/pmr))
      else
        !Rarefaction wave.
        amr = ar + 0.5_dp*gm1*(vm - Wr%u)
        pmr = Wr%p*((amr/ar)**betai)
        dpmrdum = g*pmr/amr
      end if

      !Check for convergence (i.e., pml=pmr).
      if(abs(1.-pml/pmr).le.NANO) exit

      !Compute next estimate for the intermediate state velocity, vm.
      vm = vm-(pml-pmr)/(dpmldum-dpmrdum)

      !If at 1000 iterations and not converged, print error and stop.
      if(n_iterations.eq.1000) then
        write(6,*) "ERROR: Newton-Raphson iterations did not converge in 100 steps in the exact Riemann solver!!!"
        stop
      end if

    end do

    pm = 0.5*(pml+pmr)

    Wls = WState(g*pm/sqr(aml),vm,pm)
    Wrs = WState(g*pm/sqr(amr),vm,pm)

  end if

  !Determine the intermediate solution state if required.
  if(istate.eq.1) then
    if(vm.ge.0.0_dp) then
      if(vm.lt.Wl%u) then
        aml = al*sqrt((gp1+gm1*pm/Wl%p)/(gp1+gm1*Wl%p/pm))
        vsl = Wl%u + msl*al
        if(vsl.ge.0.0_dp) then
          W = Wl
        else
          W = Wls
        end if
      else
        vhl = Wl%u - al
        vtl = vm - aml
        if(vhl.ge.0.0_dp) then
          W = Wl
        else if(vtl.le.0.0_dp) then
          W = Wls
        else
          vm = gp1i*(gm1*Wl%u + 2.0_dp*al)
          pm = Wl%p*((vm/al)**betai)
          dml = g*pm/sqr(vm)
          W  = WState(dml,vm,pm)
        end if
      end if
    else
      if(vm.gt.Wr%u) then
        amr = ar*sqrt((gp1+gm1*pm/Wr%p)/(gp1+gm1*Wr%p/pm))
        vsr = Wr%u + msr*ar
        if(vsr.le.0.0_dp) then
          W = Wr
        else
          W = Wrs
        end if
      else
        vhr = Wr%u + ar
        vtr = vm + amr
        if(vhr.le.0.0_dp) then
          W = Wr
        else if(vtr.ge.0.0_dp) then
          W = Wrs
        else
          vm = gp1i*(gm1*Wr%u - 2.0_dp*ar)
          pm = Wr%p*((-vm/ar)**betai)
          dmr = g*pm/sqr(vm)
          W = WState(dmr,vm,pm)
        end if
      end if
    end if

  else
    W = WState(-1.0_dp,-1.0_dp,-1.0_dp)
  end if

  return
end subroutine Riemann

!---------------------------------------------------------------------
!Determine the intermediate state solution flux by using the exact
!Riemann solution as originally proposed by Godunov [1959].
!---------------------------------------------------------------------
subroutine Godunov(Wl, Wr, F)
  use Euler1D_WState
  use Euler1D_UState
  implicit none
  type(Euler1D_W_State), intent(in)  :: Wl, Wr
  type(Euler1D_U_State), intent(out) :: F
  type(Euler1D_W_State)              :: Wls, Wrs, W
  !Determine the exact Riemann solution.
  call Riemann(Wl,Wr,Wls,Wrs,W,1)
  !Determine the intermediate flux.
  call Flux_W(W,F)
  return
end subroutine Godunov

!---------------------------------------------------------------------
!Determine the intermediate state solution flux by using an
!isentropic approximation.
!---------------------------------------------------------------------
subroutine IsentropicFlux(Wl, Wr, F)
  use realSizes, only: dp
  use GasConstants, only: g, gm1, gm1i, beta, betai
  use Euler1D_WState
  use Euler1D_UState
  implicit none
  type(Euler1D_W_State), intent(in)  :: Wl, Wr
  type(Euler1D_U_State), intent(out) :: F
  type(Euler1D_W_State)              :: Wls, Wrs
  real(dp)                     :: al, ar, H
  !Compute left and right state sound speeds quantities:
  al = a_W(Wl)
  ar = a_W(Wr)
  !Determine the intermediate state flux:
  if(Wl%u - al.ge.0.0_dp) then
    call Flux_W(Wl,F)
  else if(Wr%u + ar.le.0.0_dp) then
    call Flux_W(Wr,F)
  else
    !Evaluate the left state:
    H = (Wl%p/Wr%p)**beta
    Wls%p = ((al + ar - 0.5_dp*(Wr%u - Wl%u)*gm1)/(al/((Wl%p**beta)) + ar/((Wr%p**beta))))**betai
    Wls%u = (H*Wl%u/al + Wr%u/ar + 2.0_dp*(H-1.0_dp)*gm1i)/(H/al + 1.0_dp/ar)
    Wls%rho = Wl%rho*((Wls%p/Wl%p)**(1.0_dp/g))
    !Evaluate the right state:
    Wrs%p = Wls%p
    Wrs%u = Wls%u
    Wrs%rho = Wr%rho*((Wrs%p/Wr%p)**(1.0_dp/g))
    !Determine the flux:
    if(Wls%u.ge.0.) then
      call Flux_W(Wls,F)
    else
      call Flux_W(Wrs,F)
    end if
  end if
  return
end subroutine IsentropicFlux

!---------------------------------------------------------------------
!Determine the intermediate state solution flux by using Rusanov's
!single wave approximation (1964).
!---------------------------------------------------------------------
subroutine Rusanov(Wl, Wr, F)
  use realSizes, only: dp
  use Euler1D_WState
  use Euler1D_UState
  implicit none
  type(Euler1D_W_State), intent(in)  :: Wl, Wr
  type(Euler1D_U_State), intent(out) :: F
  type(Euler1D_W_State)              :: Wa
  type(Euler1D_W_State)              :: lambda_l, lambda_r, lambda_a, wavespeeds
  type(Euler1D_U_State)              :: dUrl, Fl, Fr
  real(dp)                     :: lambda_max
  !Evaluate the Roe-average primitive solution state:
  call roe_average(Wl,Wr,Wa)
  !Evaluate the jumps in the conserved solution states:
  dUrl%rho = Wr%rho - Wl%rho
  dUrl%du  = du_W(Wr) - du_W(Wl)
  dUrl%E   = E_W(Wr) - E_W(Wl)
  !Evaluate the left, right, and average state eigenvalues:
  call lambda_W(Wl,lambda_l)
  call lambda_W(Wr,lambda_r)
  call lambda_W(Wa,lambda_a)
  !Determine the average fluxes from the left and right states:
  call Flux_W(Wl,Fl)
  call Flux_W(Wr,Fr)
  F = 0.5_dp*(Fl+Fr)
  !Determine the dissipation from the single-wave approxiimation:
  call harten_fix_abs(lambda_a,lambda_l,lambda_r,wavespeeds)
  lambda_max = 0.5*max(wavespeeds%rho,max(wavespeeds%u,wavespeeds%p))
  F = F - lambda_max*dUrl
  return
end subroutine Rusanov

!---------------------------------------------------------------------
!Determine the intermediate state solution flux by using the
!HLLE-approximation [1983].
!---------------------------------------------------------------------
subroutine HLLE(Wl, Wr, F)
  use realSizes, only: dp
  use Euler1D_WState
  use Euler1D_UState
  implicit none
  type(Euler1D_W_State), intent(in)  :: Wl, Wr
  type(Euler1D_U_State), intent(out) :: F
  type(Euler1D_W_State)              :: Wa, lambda_l, lambda_r, lambda_a
  type(Euler1D_U_State)              :: dUrl, Fl, Fr
  real(dp)                     :: wavespeed_l, wavespeed_r
  !Evaluate the Roe-average primitive solution state:
  call roe_average(Wl,Wr,Wa)
  !Evaluate the jumps in the conserved solution states:
  dUrl%rho = Wr%rho - Wl%rho
  dUrl%du  = du_W(Wr) - du_W(Wl)
  dUrl%E   = E_W(Wr) - E_W(Wl)
  !Evaluate the left, right, and average state eigenvalues:
  call lambda_W(Wl,lambda_l)
  call lambda_W(Wr,lambda_r)
  call lambda_W(Wa,lambda_a)
  !Determine the left and right wave-speeds:
  wavespeed_l = min(0.,min(lambda_l%rho,lambda_a%p))
  wavespeed_r = max(0.,max(lambda_r%rho,lambda_a%p))
  !Determine the intermediate state flux:
  if(wavespeed_l.ge.0.) then
    call Flux_W(Wl,F)
  else if(wavespeed_r.le.0.) then
    call Flux_W(Wr,F)
  else
    call Flux_W(Wl,Fl)
    call Flux_W(Wr,Fr)
    F = ((wavespeed_r*Fl-wavespeed_l*Fr) + (wavespeed_l*wavespeed_r)*dUrl)/(wavespeed_r-wavespeed_l)
  end if
  return
end subroutine HLLE

!---------------------------------------------------------------------
!Determine the intermediate state solution flux by using Linde's
!HLLL-approximation [2006].
!---------------------------------------------------------------------
subroutine HLLL(Wl, Wr, F)
  use realSizes, only: dp
  use Numbers, only: TOLER
  use Euler1D_WState
  use Euler1D_UState
  implicit none
  type(Euler1D_W_State), intent(in)  :: Wl, Wr
  type(Euler1D_U_State), intent(out) :: F
  type(Euler1D_W_State)              :: Wa, lambda_l, lambda_r, lambda_a
  type(Euler1D_U_State)              :: dUrl, Fl, Fr, dFrl, dFwave
  real(dp)                     :: wavespeed_l, wavespeed_r, wavespeed_m
  real(dp)                     :: da, ca, dU, alph, alph2
  !Evaluate the Roe-average primitive solution state:
  call roe_average(Wl,Wr,Wa)
  !Evaluate the left, right, and average state eigenvalues:
  call lambda_W(Wl,lambda_l)
  call lambda_W(Wr,lambda_r)
  call lambda_W(Wa,lambda_a)
  !Determine the left and right wave-speeds:
  wavespeed_l = min(0.,min(lambda_l%rho,lambda_a%p))
  wavespeed_r = max(0.,max(lambda_r%rho,lambda_a%p))
  !Determine the intermediate state flux:
  if(wavespeed_l.ge.0.) then
    call Flux_W(Wl,F)
  else if(wavespeed_r.le.0.) then
    call Flux_W(Wr,F)
  else
    !Evaluate the jumps in the conserved solution states:
    dUrl%rho = Wr%rho - Wl%rho
    dUrl%du  = du_W(Wr) - du_W(Wl)
    dUrl%E   = E_W(Wr) - E_W(Wl)
    !Evaluate the jumps in the fluxes:
    call Flux_W(Wl,Fl)
    call Flux_W(Wr,Fr)
    dFrl = Fr - Fl
    !Determine the middle wavespeed:
    wavespeed_m = Wa%u
    !Determine the switch value:
    da = Wa%rho
    ca = a_W(Wa)
    dU = abs(dUrl%rho)/da + abs(dUrl%du)/(da*ca) + abs(dUrl%E)/(da*ca*ca)
    if(dU.le.TOLER) then
      alph = 0.0_dp
    else 
      dU = 1.0_dp/dU
      dFwave = dFrl - wavespeed_m*dUrl
      alph = 1.0_dp - (abs(dFwave%rho)/(da*ca) + abs(dFwave%du)/(da*ca*ca) + abs(dFwave%E)/(da*ca*ca*ca))*dU
      alph = max(0.0_dp,alph)
    end if
    alph2 = (wavespeed_l*wavespeed_r)*(1.0_dp-(1.0_dp-max(wavespeed_m/wavespeed_r,wavespeed_m/wavespeed_l))*alph)
    !Determine the flux:
    F = ((wavespeed_r*Fl-wavespeed_l*Fr) + alph2*dUrl)/(wavespeed_r-wavespeed_l)
  end if
  return
end subroutine HLLL

!---------------------------------------------------------------------
!Determine the intermediate state solution flux by using Toro's
!HLLC-approximation.
!---------------------------------------------------------------------
subroutine HLLC(Wl, Wr, F)
  use realSizes, only: dp
  use GasConstants, only: g, gm1, gm1i, beta, betai
  use Euler1D_WState
  use Euler1D_UState
  implicit none
  type(Euler1D_W_State), intent(in)  :: Wl, Wr
  type(Euler1D_U_State), intent(out) :: F
  type(Euler1D_U_State)              :: Uls, Urs, dU
  real(dp)                     :: al, ar, CL, CR, Z, ql, qr
  real(dp)                     :: um, pm, als, ars
  real(dp)                     :: wavespeed_l, wavespeed_m, wavespeed_r
  real(dp)                     :: wavespeed_ratio
  !Determine an approximation of the intermediate state properties
  !using an isentropic approximation:
  al  = a_W(Wl)
  ar  = a_W(Wr)
  CL  = Wl%u + 2.0_dp*gm1i*al
  CR  = Wr%u - 2.0_dp*gm1i*ar
  Z   = (ar/al)*((Wl%p/Wr%p)**beta)
  um  = (CL*Z+CR)/(1.0_dp+Z)
  als = al - 0.5_dp*gm1*(um-Wl%u)
  pm  = Wl%p*((als/al)**betai)
  ars = ar + 0.5_dp*gm1*(um-Wr%u)
  pm  = 0.5_dp*(pm + Wr%p*((ars/ar)**betai));
  !Determine the left-wave speed:
  if(pm/Wl%p.le.1.0_dp) then
    ql = 1.0_dp
  else
    ql = sqrt(1.0_dp+beta*(pm/Wl%p))
  end if
  wavespeed_l = Wl%u - ql*al
  !Determine the right-wave speed:
  if(pm/Wr%p.le.1.0_dp) then
    qr = 1.0_dp
  else
    qr = sqrt(1.0_dp+beta*(pm/Wr%p))
  end if
  wavespeed_r = Wr%u + qr*ar
  !Determine the middle-wave speed:
  wavespeed_m = um
  !Determine the intermediate state flux.
  if(wavespeed_l.ge.0.0_dp) then
    call Flux_W(Wl,F)
  else if(wavespeed_r.le.0.0_dp) then
    call Flux_W(Wr,F)
  else if(wavespeed_m.ge.0.0_dp) then
    wavespeed_ratio = (wavespeed_l-Wl%u)/(wavespeed_l-wavespeed_m)
    Uls%rho = Wl%rho*wavespeed_ratio
    Uls%du  = Wl%rho*wavespeed_m*wavespeed_ratio
    Uls%E   = (E_W(Wl) + Wl%rho*(wavespeed_m-Wl%u)*(wavespeed_m+Wl%p/(Wl%rho*(wavespeed_l-Wl%u))))*wavespeed_ratio
    call Flux_W(Wl,F)
    dU%rho = Uls%rho - Wl%rho
    dU%du  = Uls%du  - du_W(Wl)
    dU%E   = Uls%E   - E_W(Wl)
    F = F + wavespeed_l*dU
  else
    wavespeed_ratio = (wavespeed_r-Wr%u)/(wavespeed_r-wavespeed_m)
    Urs%rho = Wr%rho*wavespeed_ratio
    Urs%du  = Wr%rho*wavespeed_m*wavespeed_ratio
    Urs%E   = (E_W(Wr) + Wr%rho*(wavespeed_m-Wr%u)*(wavespeed_m+Wr%p/(Wr%rho*(wavespeed_r-Wr%u))))*wavespeed_ratio
    call Flux_W(Wr,F)
    dU%rho = Urs%rho - Wr%rho
    dU%du  = Urs%du  - du_W(Wr)
    dU%E   = Urs%E   - E_W(Wr)
    F = F + wavespeed_r*dU
  end if
  return
end subroutine HLLC

!---------------------------------------------------------------------
!Determine the intermediate state solution flux using Osher's
!integration method with the physical ordering.
!---------------------------------------------------------------------
subroutine Osher(Wl, Wr, F)
  use realSizes, only: dp
  use gasConstants, only: g, gm1, gm1i, gp1, gp1i, beta, betai
  use Euler1D_WState
  use Euler1D_UState
  implicit none
  type(Euler1D_W_State), intent(in)  :: Wl, Wr
  type(Euler1D_U_State), intent(out) :: F
  type(Euler1D_W_State)              :: Wls, Wrs, Wlsp, Wrsp
  type(Euler1D_U_State)              :: Fl, Fr, Fls, Frs, Flsp, Frsp
  real(dp)                     :: al, ar, H, als, ars, alsp, arsp
  !Determine the sound speed of the left and right states:
  al = a_W(Wl)
  ar = a_W(Wr)
  !Determine the average fluxes from the left and right states:
  call Flux_W(Wl,Fl)
  call Flux_W(Wr,Fr)
  !Determine the states and flux of the left intermediate state
  !solution using an isentropic approximation:
  H = ((Wl%p/Wr%p)**beta)
  Wls%p = ((al+ar - 0.5_dp*(Wr%u-Wl%u)*gm1)/(al/(Wl%p**beta) + ar/(Wr%p**beta)))**betai
  Wls%u = (H*Wl%u/al + Wr%u/ar + 2.0_dp*(H-1.0_dp)*gm1i)/(H/al + 1.0_dp/ar)
  Wls%rho = Wl%rho*((Wls%p/Wl%p)**(1.0_dp/g))
  als = a_W(Wls)
  call Flux_W(Wls,Fls)
  !Determine the states and flux of the right intermediate state
  !solution using an isentropic approximation:
  Wrs%p = Wls%p
  Wrs%u = Wls%u
  Wrs%rho = Wr%rho*((Wrs%p/Wr%p)**(1.0_dp/g))
  ars = a_W(Wrs)
  call Flux_W(Wrs,Frs)
  !Determine the left sonic point state and flux:
  Wlsp%u = gp1i*(Wl%u*gm1 + 2.0_dp*al)
  alsp = Wlsp%u
  Wlsp%rho = Wl%rho*((alsp/al)**(2.0_dp*gm1i))
  Wlsp%p = Wl%p*((Wlsp%rho/Wl%rho)**g)
  call Flux_W(Wlsp,Flsp)
  !Determine the right sonic point state and flux:
  Wrsp%u = gp1i*(Wr%u*gm1 - 2.0_dp*ar)
  arsp = -Wrsp%u
  Wrsp%rho = Wr%rho*((arsp/ar)**(2.0_dp*gm1i))
  Wrsp%p = Wr%p*((Wrsp%rho/Wr%rho)**g)
  call Flux_W(Wrsp,Frsp)
  !Determine the intermediate state solution flux:
  F = Fl
  !Integrate over path u-a.
  if((Wl%u-al.ge.0.0_dp).and.(Wls%u-als.ge.0.0_dp)) then
  else if((Wl%u-al.le.0.0_dp).and.(Wls%u-als.le.0.0_dp)) then
     F = F + Fls - Fl
  else if((Wl%u-al.ge.0.0_dp).and.(Wls%u-als.le.0.0_dp)) then
     F = F + Fls - Flsp
  else if((Wl%u-al.le.0.0_dp).and.(Wls%u-als.ge.0.0_dp)) then
     F = F + Flsp - Fl
  end if
  !Integrate over path u.
  if(Wls%u.lt.0.0_dp) then
     F = F + Frs - Fls
  end if
  !Integrate over path u+a.
  if((Wr%u+ar.ge.0.0_dp).and.(Wrs%u+ars.ge.0.0_dp)) then
  else if((Wr%u+ar.le.0.0_dp).and.(Wrs%u+ars.le.0.0_dp)) then
     F = F + Fr - Frs
  else if((Wr%u+ar.le.0.0_dp).and.(Wrs%u+ars.ge.0.0_dp)) then
     F = F + Fr - Frsp
  else if((Wr%u+ar.ge.0.0_dp).and.(Wrs%u+ars.le.0.0_dp)) then
     F = F + Frsp - Frs
  end if
  return
end subroutine Osher

!---------------------------------------------------------------------
!Determine the intermediate state solution flux by using the
!"linearized" approximate Riemann solver of Roe [1981].
!---------------------------------------------------------------------
subroutine Roe(Wl, Wr, F)
  use realSizes, only: dp
  use Euler1D_WState
  use Euler1D_UState
  implicit none
  type(Euler1D_W_State), intent(in)  :: Wl, Wr
  type(Euler1D_U_State), intent(out) :: F
  type(Euler1D_W_State)              :: Wa, lp, dWrl
  type(Euler1D_W_State)              :: lambda_l, lambda_r, lambda_a, wavespeeds
  type(Euler1D_U_State)              :: Fl, Fr, rc
  real(dp)                     :: dU
  !Determine the average fluxes from the left and right states:
  call Flux_W(Wl,Fl)
  call Flux_W(Wr,Fr)
  F = 0.5_dp*(Fl+Fr)
  !Evaluate the jumps in the primitive solution states:
  dWrl = Wr - Wl
  !Evaluate the Roe-average primitive solution state:
  call roe_average(Wl,Wr,Wa)
  !Evaluate the left, right, and average state eigenvalues:
  call lambda_W(Wl,lambda_l)
  call lambda_W(Wr,lambda_r)
  call lambda_W(Wa,lambda_a)
  !
  call harten_fix_abs(lambda_a,lambda_l,lambda_r,wavespeeds)
  !
  call lp_W(Wa,lp,1)
  call rc_W(Wa,rc,1)
  dU = 0.5_dp*wavespeeds%rho*(lp.dot.dWrl)
  F = F - dU*rc
  call lp_W(Wa,lp,2)
  call rc_W(Wa,rc,2)
  dU = wavespeeds%u*(lp.dot.dWrl)
  F = F - dU*rc
  call lp_W(Wa,lp,3)
  call rc_W(Wa,rc,3)
  dU = wavespeeds%p*(lp.dot.dWrl)
  F = F - dU*rc
  return
end subroutine Roe

!---------------------------------------------------------------------
!Determine the intermediate state solution flux by using 
!van Leer's flux-splitting scheme [19??].
!---------------------------------------------------------------------
subroutine VanLeerFlux(Wl, Wr, F)
  use realSizes, only: dp
  use gasConstants, only: g, gm1
  use Euler1D_WState
  use Euler1D_UState
  implicit none
  type(Euler1D_W_State), intent(in)  :: Wl, Wr
  type(Euler1D_U_State), intent(out) :: F
  type(Euler1D_U_State)              :: Fl, Fr
  real(dp)                     :: M, Mp, Mn, a, a2
  !Left flux:
  a2 = a2_W(Wl)
  a = sqrt(a2)
  M = Wl%u/a
  if(M.lt.-1.0_dp) then
    call vacuum_U(Fl)
  else if(abs(M).le.1.0_dp) then
    Mp =  0.25_dp*Wl%rho*a*(1.0_dp+ M)*(1.0_dp + M)
    Fl%rho = Mp
    Fl%du = Mp*((2.*a/g)*(0.5_dp*gm1*M + 1.0_dp))
    Fl%E = Mp*((2.*a2/(g*g-1.0_dp))*(0.5_dp*gm1*M+1.0_dp)*(0.5_dp*gm1*M+1.0_dp))
  else if(M.gt.1.0_dp) then
    call Flux_W(Wl,Fl)
  end if
  !Right flux:
  a2 = a2_W(Wr)
  a = sqrt(a2)
  M = Wr%u/a
  if(-M.lt.-1.0_dp) then
    call vacuum_U(Fr)
  else if(abs(M).le.1.0_dp) then
    Mn = -0.25_dp*Wr%rho*a*(1.0_dp - M)*(1.0_dp - M)
    Fr%rho = Mn
    Fr%du = Mn*((2.0_dp*a/g)*(0.5_dp*gm1*M - 1.0_dp))
    Fr%E = Mn*((2.0_dp*a2/(g*g-1.0_dp))*(0.5_dp*gm1*M-1.0_dp)*(0.5_dp*gm1*M-1.0_dp))
  else if(-M.gt.1.0_dp) then
    call Flux_W(Wr,Fr)
  end if
  !Compute the final flux:
  F = Fl + Fr
  return
end subroutine VanLeerFlux

!---------------------------------------------------------------------
!Determine the intermediate state solution flux by using 
!Liou's AUSM+ scheme [1994].
!---------------------------------------------------------------------
subroutine AUSMplus(Wl, Wr, F)
  use realSizes, only: dp
  use numbers, only: NANO
  use mathFunctions, only: sqr
  use Euler1D_WState
  use Euler1D_UState
  implicit none
  type(Euler1D_W_State), intent(in)    :: Wl, Wr
  type(Euler1D_U_State), intent(inout) :: F
  real(dp)                       :: beta = 0.125
  real(dp)                       :: alpha = 0.1875
  real(dp)                       :: Ml, Mr, Mplus, Mminus, Mhalf
  real(dp)                       :: al, ar, ahalf, Hl, Hr, uhalf
  real(dp)                       :: pplus, pminus, phalf
  !Compute required left and right state quantities:
  al = a_W(Wl)
  ar = a_W(Wr)
  Hl = H_W(Wl)
  Hr = H_W(Wr)
  !Determine the intermediate state sound speed:
  ahalf = 0.5*(al+ar)
  !Determine the left and right state Mach numbers based on the
  !intermediate state sound speed:
  Ml = Wl%u/ahalf
  Mr = Wr%u/ahalf
  !Determine the left state split Mach number:
  if(abs(Ml).le.1.) then
    Mplus = 0.25*sqr(Ml+1.) + beta*sqr(Ml*Ml-1.)
    pplus = 0.25*sqr(Ml+1.)*(2.-Ml) + alpha*Ml*sqr(Ml*Ml-1.)
  else
    Mplus = 0.5*(Ml + abs(Ml))
    pplus = 0.5*(1. + abs(Ml)/max(Ml,NANO))
  end if
  !Determine the right state split Mach number:
  if(abs(Mr).lt.1.) then
    Mminus = -0.25*sqr(Mr-1.) - beta*sqr(Mr*Mr-1.)
    pminus = 0.25*sqr(Mr-1.)*(2.+Mr) - alpha*Mr*sqr(Mr*Mr-1.)
  else
    Mminus = 0.5*(Mr - abs(Mr))
    pminus = 0.5*(1. - abs(Mr)/max(Mr,NANO))
  end if
  !Determine the intermediate state Mach number and pressure:
  Mhalf = Mplus + Mminus
  phalf = pplus*Wl%p + pminus*Wr%p
  uhalf = ahalf*Mhalf
  !Determine the intermediate state solution convective flux:
  F%rho = 0.5*uhalf*(Wl%rho + Wr%rho)
  F%du = 0.5*uhalf*(Wl%rho*Wl%u + Wr%rho*Wr%u)
  F%E = 0.5*uhalf*(Hl + Hr)
  !Add the numerical dissipation to the intermediate state solution flux:
  F%rho = F%rho - 0.5*abs(uhalf)*(Wr%rho - Wl%rho)
  F%du = F%du - 0.5*abs(uhalf)*(Wr%rho*Wr%u - Wl%rho*Wl%u)
  F%E = F%E - 0.5*abs(uhalf)*(Hr - Hl)
  !Add the pressure contribution to the intermediate state solution flux:
  F%du = F%du + phalf
  return
end subroutine AUSMplus

!---------------------------------------------------------------------
!Determine the intermediate state solution flux by using 
!Liou's AUSM+ scheme [1994].
!---------------------------------------------------------------------
subroutine AUSMplusup(Wl, Wr, F)
  use realSizes, only: dp
  use numbers, only: NANO
  use mathFunctions, only: sqr
  use Euler1D_WState
  use Euler1D_UState
  implicit none
  type(Euler1D_W_State), intent(in)    :: Wl, Wr
  type(Euler1D_U_State), intent(inout) :: F
  real(dp)                       :: beta = 0.125
  real(dp)                       :: alpha = 0.1875
  real(dp)                       :: Ml, Mr, Mplus, Mminus, Mhalf
  real(dp)                       :: al, ar, ahalf, Hl, Hr, uhalf
  real(dp)                       :: pplus, pminus, phalf
  !Compute required left and right state quantities:
  al = a_W(Wl)
  ar = a_W(Wr)
  Hl = H_W(Wl)
  Hr = H_W(Wr)
  !Determine the intermediate state sound speed:
  ahalf = 0.5*(al+ar)
  !Determine the left and right state Mach numbers based on the
  !intermediate state sound speed:
  Ml = Wl%u/ahalf
  Mr = Wr%u/ahalf
  !Determine the left state split Mach number:
  if(abs(Ml).le.1.) then
    Mplus = 0.25*sqr(Ml+1.) + beta*sqr(Ml*Ml-1.)
    pplus = 0.25*sqr(Ml+1.)*(2.-Ml) + alpha*Ml*sqr(Ml*Ml-1.)
  else
    Mplus = 0.5*(Ml + abs(Ml))
    pplus = 0.5*(1. + abs(Ml)/max(Ml,NANO))
  end if
  !Determine the right state split Mach number:
  if(abs(Mr).lt.1.) then
    Mminus = -0.25*sqr(Mr-1.) - beta*sqr(Mr*Mr-1.)
    pminus = 0.25*sqr(Mr-1.)*(2.+Mr) - alpha*Mr*sqr(Mr*Mr-1.)
  else
    Mminus = 0.5*(Mr - abs(Mr))
    pminus = 0.5*(1. - abs(Mr)/max(Mr,NANO))
  end if
  !Determine the intermediate state Mach number and pressure:
  Mhalf = Mplus + Mminus
  phalf = pplus*Wl%p + pminus*Wr%p
  uhalf = ahalf*Mhalf
  !Determine the intermediate state solution convective flux:
  F%rho = 0.5*uhalf*(Wl%rho + Wr%rho)
  F%du = 0.5*uhalf*(Wl%rho*Wl%u + Wr%rho*Wr%u)
  F%E = 0.5*uhalf*(Hl + Hr)
  !Add the numerical dissipation to the intermediate state solution flux:
  F%rho = F%rho - 0.5*abs(uhalf)*(Wr%rho - Wl%rho)
  F%du = F%du - 0.5*abs(uhalf)*(Wr%rho*Wr%u - Wl%rho*Wl%u)
  F%E = F%E - 0.5*abs(uhalf)*(Hr - Hl)
  !Add the pressure contribution to the intermediate state solution flux:
  F%du = F%du + phalf
  return
end subroutine AUSMplusup
