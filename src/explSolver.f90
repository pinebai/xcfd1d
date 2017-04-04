!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
module explSolver
  use realSizes, only: dp
  use Euler1D_UState

  implicit none

  !Solution block variable declaration.
  type(Euler1D_U_State), dimension(:),   allocatable :: Uo   !Array of conserved variables
  type(Euler1D_U_State), dimension(:,:), allocatable :: dUdt !Array of solution residuals
  real(dp),              dimension(:),   allocatable :: dt   !Array of local time-steps

contains

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  real(dp) function RKcoef(i_stage, n_stage) result(beta)
    use realSizes, only: dp
    implicit none
    integer, intent(in) :: i_stage
    integer, intent(in) :: n_stage
    select case(n_stage)
    case(1)
      beta = 1.0_dp
    case(2)
      beta = 1.0_dp
      if(i_stage.eq.2) beta = 0.5_dp*beta
    case(4)
      beta = 1.0_dp
      if(i_stage.eq.4) beta = beta/6.0_dp
    end select
    return
  end function RKcoef

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  real(dp) function MSOScoef(i_stage, n_stage) result(beta)
    use realSizes, only: dp
    use cfdParams, only: LIMITER_ZERO
    use inputParams, only: i_limiter
    implicit none
    integer, intent(in) :: i_stage
    integer, intent(in) :: n_stage
    select case(n_stage)
    case(1)
      beta = 1.0_dp
    case(2)
      if(i_limiter.eq.LIMITER_ZERO) then
        beta = 1.0_dp
        if(i_stage.eq.1) beta = beta + 0.3333_dp
      else
        beta = 0.4693_dp
        if(i_stage.eq.1) beta = beta + 0.4242_dp
      end if
    case(3)
      if(i_limiter.eq.LIMITER_ZERO) then
        beta = 1.50_dp
        if(i_stage.eq.1) beta = beta + 0.1481_dp
        if(i_stage.eq.2) beta = beta + 0.4000_dp
      else
        beta = 0.6936_dp
        if(i_stage.eq.1) beta = beta + 0.1918_dp
        if(i_stage.eq.2) beta = beta + 0.4929_dp
      end if
    case(4)
      if(i_limiter.eq.LIMITER_ZERO) then
        beta = 2.0_dp
        if(i_stage.eq.1) beta = beta + 0.0833_dp
        if(i_stage.eq.2) beta = beta + 0.2069_dp
        if(i_stage.eq.3) beta = beta + 0.4265_dp
      else
        beta = 0.9214_dp
        if(i_stage.eq.1) beta = beta + 0.1084_dp
        if(i_stage.eq.2) beta = beta + 0.2602_dp
        if(i_stage.eq.3) beta = beta + 0.5052_dp
      end if
    case(5)
      if(i_limiter.eq.LIMITER_ZERO) then
        beta = 2.5_dp
        if(i_stage.eq.1) beta = beta + 0.0533_dp
        if(i_stage.eq.2) beta = beta + 0.1263_dp
        if(i_stage.eq.3) beta = beta + 0.2375_dp
        if(i_stage.eq.4) beta = beta + 0.4414_dp
      else
        beta = 1.1508_dp
        if(i_stage.eq.1) beta = beta + 0.0695_dp
        if(i_stage.eq.2) beta = beta + 0.1602_dp
        if(i_stage.eq.3) beta = beta + 0.2898_dp
        if(i_stage.eq.4) beta = beta + 0.5060_dp
      end if
    case(6)
      if(i_limiter.eq.LIMITER_ZERO) then
        beta = 3.00_dp
        if(i_stage.eq.1) beta = beta + 0.0370_dp
        if(i_stage.eq.2) beta = beta + 0.0851_dp
        if(i_stage.eq.3) beta = beta + 0.1521_dp
        if(i_stage.eq.4) beta = beta + 0.2562_dp
        if(i_stage.eq.5) beta = beta + 0.4512_dp
      else
        beta = 1.3805_dp
        if(i_stage.eq.1) beta = beta + 0.0482_dp
        if(i_stage.eq.2) beta = beta + 0.1085_dp
        if(i_stage.eq.3) beta = beta + 0.1885_dp
        if(i_stage.eq.4) beta = beta + 0.3050_dp
        if(i_stage.eq.5) beta = beta + 0.5063_dp
      end if
    end select
    return
  end function MSOScoef
  
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  subroutine solnExplicit_allocate
    use realSizes, only: dp
    use inputParams, only: max_time_steps
    use solnBlock_module, only: NCl, NCu, Ng, Nc
    implicit none
    integer :: n, m
    allocate(Uo(Nc))
    allocate(dUdt(Nc,3))
    allocate(dt(Nc))
    do n = NCl-Ng, NCu+Ng
      call standard_atmosphere_U(Uo(n))
      do m = 1, 3
        call vacuum_U(dUdt(n,m))
      end do
      dt(n) = 0.0_dp
    end do
    return
  end subroutine solnExplicit_allocate

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  subroutine solnExplicit_deallocate
    implicit none
    if(allocated(Uo)) deallocate(Uo)
    if(allocated(dUdt)) deallocate(dUdt)
    if(allocated(dt)) deallocate(dt)
    return
  end subroutine solnExplicit_deallocate

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  subroutine setTimeStep(dt_global,cfl_number)
    use realSizes, only: dp
    use numbers, only: MEGA
    use gridBlock_module, only: Cell
    use solnBlock_module, only: NCl, NCu, Ng, W
    use Euler1D_WState
    implicit none
    real(dp), intent(out) :: dt_global
    real(dp), intent(in)  :: cfl_number
    integer               :: n
    real(dp)              :: lambda
    dt_global = MEGA
    do n = NCl-Ng, NCu+Ng
      !Determine the local time-step for the current cell:
      lambda = abs(W(n)%u) + a_W(W(n))
      dt(n) = cfl_number*Cell(n)%dx/lambda
      !Determine the minimum time-step for the block:
      dt_global = min(dt_global,dt(n))
    end do
    return
  end subroutine setTimeStep

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  subroutine setGlobalTimeStep(dt_global)
    use realSizes, only: dp
    use solnBlock_module, only: NCl, NCu, Ng
    implicit none
    real(dp), intent(in) :: dt_global
    integer :: n
    do n = NCl-Ng, NCu+Ng
      dt(n) = dt_global
    end do
    return
  end subroutine setGlobalTimeStep

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  subroutine residual_l1_norm(l1norm)
    use realSizes, only: dp
    use solnBlock_module, only: NCl, NCu
    use Euler1D_UState
    implicit none
    type(Euler1D_U_State), intent(out) :: l1norm
    integer :: n
    call vacuum_U(l1norm)
    do n = NCl, NCu
      l1norm%rho = l1norm%rho + abs(dUdt(n,1)%rho)
      l1norm%du  = l1norm%du  + abs(dUdt(n,1)%du)
      l1norm%E   = l1norm%E   + abs(dUdt(n,1)%E)
    end do
    return
  end subroutine residual_l1_norm

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  subroutine residual_l2_norm(l2norm)
    use realSizes, only: dp
    use solnBlock_module, only: NCl, NCu
    use Euler1D_UState
    implicit none
    type(Euler1D_U_State), intent(out) :: l2norm
    integer :: n
    call vacuum_U(l2norm)
    do n = NCl, NCu
      l2norm%rho = l2norm%rho + (dUdt(n,1)%rho)**2
      l2norm%du  = l2norm%du  + (dUdt(n,1)%du)**2
      l2norm%E   = l2norm%E   + (dUdt(n,1)%E)**2
    end do
    l2norm%rho = sqrt(l2norm%rho)
    l2norm%du  = sqrt(l2norm%du)
    l2norm%E   = sqrt(l2norm%E)
    return
  end subroutine residual_l2_norm

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  subroutine residual_max_norm(maxnorm)
    use realSizes, only: dp
    use solnBlock_module, only: NCl, NCu
    use Euler1D_UState
    implicit none
    type(Euler1D_U_State), intent(out) :: maxnorm
    integer :: n
    call vacuum_U(maxnorm)
    do n = NCl, NCu
      maxnorm%rho = max(maxnorm%rho,abs(dUdt(n,1)%rho))
      maxnorm%du  = max(maxnorm%du,abs(dUdt(n,1)%du))
      maxnorm%E   = max(maxnorm%E,abs(dUdt(n,1)%E))
    end do
    return
  end subroutine residual_max_norm

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  subroutine open_residual_file
    implicit none
    open(unit=4,file='residual.dat',status='REPLACE')
    write(4,'(a5,10(a11))') 'n', 'time', 'l1-rho', 'l1-du', 'l1-E', &
         'l2-rho', 'l2-du', 'l2-E', 'mx-rho', 'mx-du', 'mx-E'
    return
  end subroutine open_residual_file

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  subroutine write_residual_file(n,t,l1,l2,mx)
    use Euler1D_UState
    implicit none
    integer,               intent(in) :: n
    real(dp),              intent(in) :: t
    type(Euler1D_U_State), intent(in) :: l1, l2, mx
    write(4,'(i5,10(es11.4))') n, t, l1%rho, l1%du, l1%E, &
         l2%rho, l2%du, l2%E, mx%rho, mx%du, mx%E
    return
  end subroutine write_residual_file

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  subroutine close_residual_file
    implicit none
    close(4)
    return
  end subroutine close_residual_file

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  subroutine multistageExplicit(i_stage, ierr)
    use cfdParams
    use inputParams, only: i_explicit, n_stage, i_recon, i_flux, cfl_number
    use Euler1D_WState
    use Euler1D_UState
    use gridBlock_module, only: Cell, xfaceL, xfaceR
    use solnBlock_module, only: NCl, NCu, Ng, W, U, phi, dWdx, BCl, BCr, WoL, WoR
    implicit none
    !Argument variables:
    integer, intent(in)   :: i_stage
    integer, intent(out)  :: ierr
    !Local variables:
    integer               :: n
    integer               :: k_residual
    real(dp)              :: dx, sigma
    type(Euler1D_W_State) :: Wl, Wr
    type(Euler1D_U_State) :: Flux

    !Initialize error flag.
    ierr = 0

    !Evaluate the time step fraction and residual storage location for
    !the stage.
    select case(i_explicit)
    case(EXPLICIT_EULER)
      k_residual = 1
    case(PREDICTOR_CORRECTOR)
      k_residual = 1
    case(RUNGE_KUTTA)
      k_residual = 1
      if(n_stage.eq.4) then
        if(i_stage.eq.4) then
          k_residual = 1
        else
          k_residual = i_stage
        end if
      end if
    case(MULTISTAGE_OPTIMAL_SMOOTHING)
      k_residual = 1
    end select

    !Perform the reconstruction within each cell of the computational
    !grid for this stage.
    if(i_recon.eq.RECONSTRUCTION_GG) then
      call GGreconstruction
    else if(i_recon.eq.RECONSTRUCTION_LSQ) then
      call LSQreconstruction
    else if(i_recon.eq.RECONSTRUCTION_WENO) then
      call WENOreconstruction
    else if(i_recon.eq.RECONSTRUCTION_PPM) then
      call PPMreconstruction
    else if(i_recon.eq.RECONSTRUCTION_GAMMA) then
      call GGreconstruction
      call Gamma_Differencing_Scheme
    end if

    !Determine limiters for the gradient reconstruction.
    call calcLimiters

    !Evaluate the time rate of change of the solution (i.e., the
    !solution residuals) using a second-order limited upwind scheme
    !with a variety of flux functions.

    do n = NCl-1, NCu

      if(i_stage.eq.1) then
        Uo(n+1) = U(n+1)
        call vacuum_U(dUdt(n+1,k_residual))
      else
        select case(i_explicit)
        case(RUNGE_KUTTA)
          if(n_stage.eq.2) then
          else if((n_stage.eq.4).and.(i_stage.eq.4)) then
            dUdt(n+1,k_residual) = dUdt(n+1,1) + 2.0_dp*dUdt(n+1,2) + 2.0_dp*dUdt(n+1,3)
          else
            call vacuum_U(dUdt(n+1,k_residual))
          end if
        case(MULTISTAGE_OPTIMAL_SMOOTHING)
          call vacuum_U(dUdt(n+1,k_residual))
        end select
      end if

      if((n.eq.NCl-1).and.(BCl.ne.BC_PERIODIC)) then
        dx = xfaceL(n+1) - Cell(n+1)%Xc
        Wr = W(n+1) + dx*(phi(n+1)*dWdx(n+1))
        select case(BCl)
        case(BC_FIXED)
          Wl = WoL
        case(BC_CONSTANT_EXTRAPOLATION)
          Wl = Wr
        case(BC_CHARACTERISTIC)
          !call characteristic_W(Wl,Wr)
        case(BC_REFLECTION)
          call Reflect_W(Wl,Wr)
        end select
        
      else if((n.eq.NCu).and.(BCr.ne.BC_PERIODIC)) then
        dx = xfaceR(n) - Cell(n)%Xc
        Wl = W(n) + dx*(phi(n)*dWdx(n))
        select case(BCr)
        case(BC_FIXED)
          Wr = WoR
        case(BC_CONSTANT_EXTRAPOLATION)
          Wr = Wl
        case(BC_CHARACTERISTIC)
          !call characteristic_W(Wr,Wl)
        case(BC_REFLECTION)
          call Reflect_W(Wr,Wl)
        end select

      else
        dx = xfaceR(n) - Cell(n)%Xc
        Wl = W(n) + dx*(phi(n)*dWdx(n))
        dx = xfaceL(n+1) - Cell(n+1)%Xc
        Wr = W(n+1) + dx*(phi(n+1)*dWdx(n+1))
      end if

      !Determine RIGHT face HYPERBOLIC flux.
      select case(i_flux)
      case(FLUX_GODUNOV)
        call Godunov(Wl,Wr,Flux)
      case(FLUX_ISENTROPIC)
        call IsentropicFlux(Wl,Wr,Flux)
      case(FLUX_RUSANOV)
        call Rusanov(Wl,Wr,Flux)
      case(FLUX_HLLE)
        call HLLE(Wl,Wr,Flux)
      case(FLUX_HLLL)
        call HLLL(Wl,Wr,Flux)
      case(FLUX_HLLC)
        call HLLC(Wl,Wr,Flux)
      case(FLUX_OSHER)
        call Osher(Wl,Wr,Flux)
      case(FLUX_ROE)
        call Roe(Wl,Wr,Flux)
      case(FLUX_VANLEER)
        call VanLeerFlux(Wl,Wr,Flux)
      case(FLUX_AUSMplus)
        call AUSMplus(Wl,Wr,Flux)
      case(FLUX_AUSMplusup)
        call AUSMplusup(Wl,Wr,Flux)
      end select

      !Evaluate cell-averaged solution changes for the current cell.
      sigma = dt(n)/Cell(n)%dx
      dUdt(n,k_residual) = dUdt(n,k_residual) - sigma*Flux

      !Evaluate cell-averaged solution changes for the neighbour cell.
      sigma = dt(n+1)/Cell(n+1)%dx
      dUdt(n+1,k_residual) = dUdt(n+1,k_residual) + sigma*Flux

      !Include area-change source term.
      dUdt(n,k_residual) = dUdt(n,k_residual) + dt(n)*Sa(U(n),Cell(n)%xA,Cell(n)%dA)

    end do

    call vacuum_U(dUdt(NCl-1,k_residual))
    call vacuum_U(dUdt(NCu+1,k_residual))

    return
  end subroutine multistageExplicit


  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  subroutine updateExplicit(i_stage, ierr)
    use cfdParams
    use inputParams, only: i_limiter, i_explicit, n_stage, i_time_step
    use Euler1D_UState
    use gridBlock_module, only: Cell
    use solnBlock_module, only: NCl, NCu, U, W
    implicit none
    !Argument variables:
    integer, intent(in)  :: i_stage
    integer, intent(out) :: ierr
    integer              :: k_residual, nc, n_residual_reduction
    !Local variables:
    real(dp) :: omega
    real(dp) :: residual_reduction_factor

    ierr = 0

    !Evaluate the time step fraction and residual storage location for
    !the stage.
    select case(i_explicit)
    case(EXPLICIT_EULER)
      omega = RKcoef(i_stage,n_stage)
      k_residual = 1
    case(PREDICTOR_CORRECTOR)
      omega = RKcoef(i_stage,n_stage)
      k_residual = 1
    case(RUNGE_KUTTA)
      omega = RKcoef(i_stage,n_stage)
      k_residual = 1
      if(n_Stage.eq.4) then
        if(i_stage.eq.4) then
          k_residual = 1
        else
          k_residual = i_stage
        end if
      end if
    case(MULTISTAGE_OPTIMAL_SMOOTHING)
      omega = MSOScoef(i_stage,n_stage)
      k_residual = 1
    end select

    !Update solution variables for this stage.
    do nc = NCl, NCu

      !Perform the explicit solution update.
      U(nc) = Uo(nc) + omega*dUdt(nc,k_residual)

      !Perform residual reductions (reduce the time-step) if
      !necessary for scalar or semi-implicit local time-stepping.
      if(i_time_step.eq.LOCAL_TIME_STEP) then
        if(unphysical_properties_U(U(nc)).ne.0) then
          n_residual_reduction = 1
          residual_reduction_factor = 1.0_dp
          do while(n_residual_reduction.lt.11)
            !Reduce the residual by half.
            residual_reduction_factor = 0.5_dp*residual_reduction_factor
            dt(nc) = dt(nc)*residual_reduction_factor
            dUdt(nc,k_residual) = residual_reduction_factor*dUdt(nc,k_residual)
            !Re-perform the explicit solution update.
            U(nc) = Uo(nc) + omega*dUdt(nc,k_residual)
            !Check for unphysical properties.
            if(unphysical_properties_U(U(nc)).ne.0) exit
            n_residual_reduction = n_residual_reduction + 1
          end do
        end if
      end if
      
      !Check for unphysical properties.
      ierr = unphysical_properties_U(U(nc))
      if(ierr.ne.0) then
        write(6,*) ' Unphysical property detected: '
        write(6,*) ' cell = (', nc, ') '
        write(6,*) ' X    = ', Cell(nc)%Xc
        write(6,*) ' U    = ', U(nc)%rho, U(nc)%du, U(nc)%E
        write(6,*) ' W    = ', W(nc)%rho, W(nc)%u, W(nc)%p
        write(6,*) ' Uo   = ', Uo(nc)%rho, Uo(nc)%du, Uo(nc)%E
        write(6,*) ' dUdt = ', dUdt(nc,k_residual)%rho, dUdt(nc,k_residual)%du, dUdt(nc,k_residual)%E
        return
      end if

      !Update the primitive variabels.
      call transform_U_to_W(U(nc),W(nc))

    end do

  end subroutine updateExplicit


end module explSolver
