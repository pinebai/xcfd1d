!-----------------------------------------------------------------------
! Module for reading input parameters
!-----------------------------------------------------------------------
module inputParams
  !Include real sizes module
  use realSizes, only: dp
  !Include the (rho,u,p)-primitive-variable solution state
  use Euler1D_WState

  !Implicit declaration
  implicit none

  !Domain related input parameters:
  integer                 :: num_cells           !Number of cells in the domain
  integer                 :: num_ghost           !Number of ghost cells
  integer                 :: i_grid              !Grid type
  real(dp)                :: xl, xr              !Left and right domain locations
  integer                 :: BCleft, BCright     !Left and right boundary conditions
  integer                 :: Str_fcn             !Stretching function
  real(dp)                :: beta_str, tau_str   !Stretching parameters
  namelist / domain / num_cells, num_ghost, i_grid, xl, xr, BCleft, BCright, &
       Str_fcn, beta_str, tau_str

  !Gas-type input parameters:
  integer                 :: i_gas_type          !Gas type
  namelist / gas / i_gas_type

  !Initial condition input parameters:
  integer                 :: init_type           !Initial condition type
  real(dp)                :: yl, ym, yr          !Wave-locations
  type(Euler1D_W_State)   :: Wl, Wm, Wr          !Left, middle, and right initial states
  namelist / initial / init_type, yl, ym, yr, Wl, Wm, Wr

  !Temporal discretization parameters:
  integer                 :: i_explicit          !Time integration type
  integer                 :: n_stage             !Number of stages in time-marching scheme
  real(dp)                :: cfl_number          !CFL number
  integer                 :: i_time_step         !Global or local time-stepping
  real(dp)                :: time_max            !Maximum solution time
  integer                 :: max_time_steps      !Maximum number of time steps
  integer                 :: i_output_freq       !Frequency of updates on screen
  namelist / temporal / i_explicit, n_stage, cfl_number, time_max, &
       i_time_step, max_time_steps, i_output_freq

  !Spatial discretization parameters:
  integer                 :: i_flux              !Flux function type
  integer                 :: i_limiter           !Limiter type
  integer                 :: i_recon             !Reconstruction type
  real(dp)                :: betam               !Parameter for gamma-differencing scheme
  namelist / spatial / i_flux, i_limiter, i_recon, betam

  !Output formats:
  logical                 :: plot_tecplot        !Write output into Tecplot format
  logical                 :: plot_gnuplot        !Write output into gnuplot format
  logical                 :: plot_eps            !Make encapsulated post-script plots
  namelist / output / plot_tecplot, plot_gnuplot, plot_eps

  !EPS Figures:
  integer,  parameter     :: max_plots = 20
  integer                 :: nplots
  real(dp), dimension(20) :: Xf, Yf, Xo, Yo, Xa, Ya
  real(dp), dimension(20) :: xmin, xmax, dx, xdim, ymin, ymax, dy, ydim
  integer,  dimension(20) :: xvar, yvar, xsig, ysig
  namelist / eps_plots / nplots, Xf, Yf, Xo, Yo, Xa, Ya, &
       xvar, xmin, xmax, dx, yvar, ymin, ymax, dy, xdim, ydim, xsig, ysig

contains

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  subroutine inputDefaults
    use cfdParams
    use gasConstants, only: GAS_AIR, setGas
    implicit none
    !Domain related input parameters:
    num_cells = 10
    num_ghost = 2
    BCleft = BC_REFLECTION
    BCright = BC_REFLECTION
    Xl = 0.0_dp
    Xr = 1.0_dp
    Str_fcn = STR_LINEAR
    beta_str = 1.0_dp
    tau_str = 1.0_dp
    !Gas-type:
    i_gas_type = GAS_AIR
    call setGas(i_gas_type)
    call standard_atmosphere_W(Wl)
    call standard_atmosphere_W(Wm)
    call standard_atmosphere_W(Wr)
    !Initial conditions:
    init_type = IC_UNIFORM
    yl = Xl
    yr = Xr
    ym = 0.5_dp*(yr-yl)
    !Time-integration parameters:
    i_explicit = EXPLICIT_EULER
    n_stage = 1
    cfl_Number = 0.5_dp
    i_time_step = GLOBAL_TIME_STEP
    time_max = 0.0_dp
    max_time_steps = 0
    i_output_freq = 10
    !Spatial discretization parameters:
    i_flux = FLUX_ROE
    i_limiter = LIMITER_ZERO
    i_recon = RECONSTRUCTION_LSQ
    betam = 0.5_dp
    !Output formats:
    plot_tecplot = .false.
    plot_gnuplot = .false.
    plot_eps     = .false.
    !EPS Figures:
    nplots = 0
    Xf(:) = 0.0_dp; Yf(:) = 0.0_dp
    Xo(:) = 0.0_dp; Yo(:) = 0.0_dp
    Xa(:) = 0.0_dp; Ya(:) = 0.0_dp
    xmin(:) = 0.0_dp; xmax(:) = 0.0_dp; dx(:) = 0.0_dp; xdim(:) = 1.0_dp
    ymin(:) = 0.0_dp; ymax(:) = 0.0_dp; dy(:) = 0.0_dp; ydim(:) = 1.0_dp
    xvar(:) = 0; yvar(:) = 0
    xsig(:) = 0; ysig(:) = 0
    return
  end subroutine inputDefaults


  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  subroutine inputNamelists
    implicit none
    !Open the input file.
    open(unit=7,file='input.in')
    !Read in the grid namelist:
    read(unit=7,nml=domain)
    !Read in the initial conditions namelist:
    read(unit=7,nml=initial)
    !Read in the temporal discretization namelist:
    read(unit=7,nml=temporal)
    !Read in the spatial discretization namelist:
    read(unit=7,nml=spatial)
    !Read in the output format namelist:
    read(unit=7,nml=output)
    !Read in the EPS plots namelist:
    read(unit=7,nml=eps_plots)
    !Close the input file.
    close(7)
    return
  end subroutine inputNamelists


  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  subroutine inputConsistency(ierr)
    use Numbers, only: TOLER
    use cfdParams
    use gasConstants, only: GAS_AIR
    implicit none
    integer, intent(out) :: ierr
    ierr = 0

    if(init_type.lt.IC_UNIFORM) then
      write(6,"(a)") "    - Initial condition type specified incorrectly."
      ierr = ierr + 1
    end if

    if((i_gas_type.lt.GAS_AIR).or.(i_gas_type.gt.GAS_AIR)) then
      write(6,"(a)") "    ~- Gas type specified incorrectly."
      ierr = ierr + 1
    end if

    if(num_cells.lt.5) then
      write(6,"(a)") "    - To few cells have been requested."
      ierr = ierr + 1
    end if

    if(num_ghost.lt.2.or.num_ghost.gt.3) then
      write(6,"(a)") "    - An incorrect number of ghost cells have been specified."
      ierr = ierr + 1
    end if

    if((BCleft.lt.BC_FIXED).or.(BCleft.gt.BC_PERIODIC)) then
      write(6,"(a)") "    - An invalid left boundary condition has been specified."
      ierr = ierr + 1
    end if

    if((BCright.lt.BC_FIXED).or.(BCright.gt.BC_PERIODIC)) then
      write(6,"(a)") "    - An invalid right boundary condition has been specified."
      ierr = ierr + 1
    end if

    if((i_grid.lt.GRID_SHOCK_TUBE).or.(i_grid.gt.GRID_NOZZLE)) then
      write(6,"(a)") "    - Grid type specified incorrectly."
      ierr = ierr + 1
    end if

    if((Str_fcn.lt.STR_LINEAR).or.(Str_fcn.gt.STR_COSINE)) then
      write(6,"(a)") "    - An invalid stretching function has been specified."
      ierr = ierr + 1
    end if

    if((Str_fcn.ne.STR_LINEAR).and.(beta_str.le.1.)) then
      write(6,"(a)") "    - An invalid stretching beta-parameter has been specified."
      ierr = ierr + 1
    end if

    if((Str_fcn.ne.STR_LINEAR).and.(tau_str.le.1.)) then
      write(6,"(a)") "    - An invalid stretching tau-parameter has been specified."
      ierr = ierr + 1
    end if

    if((i_limiter.ne.LIMITER_ZERO) .and. &
       (i_limiter.ne.LIMITER_ONE) .and. &
       (i_limiter.ne.LIMITER_VANLEER) .and. &
       (i_limiter.ne.LIMITER_VANALBADA) .and. &
       (i_limiter.ne.LIMITER_BARTH_JESPERSEN) .and. &
       (i_limiter.ne.LIMITER_VENKATAKRISHNAN)) then
      write(6,"(a)") "    - An invalid limiter has been specified."
      ierr = ierr + 1
    end if

    if(i_limiter.ne.LIMITER_ZERO) then
      if((i_recon.ne.RECONSTRUCTION_GG) .and. &
         (i_recon.ne.RECONSTRUCTION_LSQ) .and. &
         (i_recon.ne.RECONSTRUCTION_WENO) .and. &
         (i_recon.ne.RECONSTRUCTION_PPM) .and. &
         (i_recon.ne.RECONSTRUCTION_GAMMA)) then
        write(6,"(a)") "    - An invalid reconstruction type has been specified."
        ierr = ierr + 1
      end if
    end if

    if(i_recon.eq.RECONSTRUCTION_GAMMA) then
      if(i_limiter.ne.LIMITER_ONE) then
        i_limiter = LIMITER_ONE
        write(6,"(a)") "    - The limiter was reset so that the gradient is unlimited for"
        write(6,"(a)") "      the gamma differencing scheme."
      end if
    end if

    if(i_recon.eq.RECONSTRUCTION_GAMMA) then
      if((betam.lt.0.1).or.(betam.gt.0.5)) then
        write(6,"(a)") "    - The betam parameter for the gamma differencing scheme must"
        write(6,"(a)") "      be set within 0.1 and 0.5 (default is 0.5)."
        ierr = ierr + 1
      end if
    end if

    if(i_output_freq.lt.1) then
      write(6,"(a)") "    - The output frequency has been incorrectly specified.  It should"
      write(6,"(a)") "      be greater than or equal to 1."
      ierr = ierr + 1
    end if

    if(ierr.ne.0) then
       call flush(6)
    end if

    return
  end subroutine inputConsistency


  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  subroutine inputWrite
    use cfdParams
    implicit none
    if(init_type.eq.IC_UNIFORM) then
      write(6,"(a)") "    - Uniform"
    else if(init_type.eq.IC_SOD) then
      write(6,"(a)") "    - Sod problem"
    else if(init_type.eq.IC_MODIFIED_SOD) then
      write(6,"(a)") "    - Modified Sod problem"
    else if(init_type.eq.IC_STRONG_SOD) then
      write(6,"(a)") "    - Strong Sod"
    else if(init_type.eq.IC_123_PROBLEM) then
      write(6,"(a)") "    - 123 problem"
    else if(init_type.eq.IC_THREE_RIGHT_WAVES) then
      write(6,"(a)") "    - Three right travelling waves"
    else if(init_type.eq.IC_STATIONARY_CONTACT) then
      write(6,"(a)") "    - Stationary contact surface"
    else if(init_type.eq.IC_SUBSONIC_NOZZLE) then
      write(6,"(a)") "    - Subsonic nozzle flow"
    else if(init_type.eq.IC_TRANSONIC_NOZZLE) then
      write(6,"(a)") "    - Transonic nozzle flow"
    else if(init_type.eq.IC_SQUARE_WAVE) then
      write(6,"(a)") "    - Square Wave"
    else if(init_type.eq.IC_SINE_SQUARED_WAVE) then
      write(6,"(a)") "    - Sine-Squared Wave"
    else if(init_type.eq.IC_SEMI_ELLIPSE_WAVE) then
      write(6,"(a)") "    - Semi-Ellipse Wave"
    end if
    if(i_explicit.eq.EXPLICIT_EULER) then
      write(6,"(a)") "    - Explicit Euler"
    else if(i_explicit.eq.PREDICTOR_CORRECTOR) then
      write(6,"(a)") "    - Predictor-Corrector"
    else if(i_explicit.eq.RUNGE_KUTTA) then
      write(6,"(a)") "    - Runge-Kutta"
    else if(i_explicit.eq.MULTISTAGE_OPTIMAL_SMOOTHING) then
      write(6,"(a)") "    - Multistage Optimal Smoothing"
    end if
    if(i_flux.eq.FLUX_GODUNOV) then
      write(6,"(a)") "    - Fluxes calculated using the exact Riemann solution"
    else if(i_flux.eq.FLUX_ISENTROPIC) then
      write(6,"(a)") "    - Fluxes calculated using an isentropic approximation"
    else if(i_flux.eq.FLUX_RUSANOV) then
      write(6,"(a)") "    - Fluxes calculated using Rusanov's approximate solver"
    else if(i_flux.eq.FLUX_HLLE) then
      write(6,"(a)") "    - Fluxes calculated using the HLLE approximate solver"
    else if(i_flux.eq.FLUX_HLLL) then
      write(6,"(a)") "    - Fluxes calculated using the HLLL approximate solver"
    else if(i_flux.eq.FLUX_HLLC) then
      write(6,"(a)") "    - Fluxes calculated using the HLLC approximate solver"
    else if(i_flux.eq.FLUX_OSHER) then
      write(6,"(a)") "    - Fluxes calculated using Osher's approximate solver with physical ordering"
    else if(i_flux.eq.FLUX_ROE) then
      write(6,"(a)") "    - Fluxes calculated using Roe's approximate solver"
    else if(i_flux.eq.FLUX_VANLEER) then
      write(6,"(a)") "    - Fluxes calculated using van Leer's flux vector splitting scheme"
    else if(i_flux.eq.FLUX_AUSMplus) then
      write(6,"(a)") "    - Fluxes calculated using Liou's AUSM+"
    else if(i_flux.eq.FLUX_AUSMplusup) then
      write(6,"(a)") "    - Fluxes calculated using Liou's AUSM+up"
!   else if(i_flux.eq.FLUX_SLAU) then
!     write(6,"(a)") "    - Fluxes calculated using Shima's version of AUSM"
    end if
    if(i_limiter.ne.LIMITER_ZERO) then
      if(i_recon.eq.RECONSTRUCTION_GG) then
        write(6,"(a)") "    - Green-Gauss gradient reconstruction"
      else if(i_recon.eq.RECONSTRUCTION_LSQ) then
        write(6,"(a)") "    - Least squares gradient reconstruction"
      else if(i_recon.eq.RECONSTRUCTION_WENO) then
        write(6,"(a)") "    - WENO gradient reconstruction"
      else if(i_recon.eq.RECONSTRUCTION_PPM) then
        write(6,"(a)") "    - Piecewise parabolic reconstruction"
      else if(i_recon.eq.RECONSTRUCTION_GAMMA) then
        write(6,"(a,F4.2)") "    - Gamma differencing scheme with betam = ", betam
      end if
    end if
    if(i_limiter.eq.LIMITER_ZERO) then
      write(6,"(a)") "    - Fully limited (first order)"
    else if(i_limiter.eq.LIMITER_ONE) then
      write(6,"(a)") "    - Fully unlimited"
    else if(i_limiter.eq.LIMITER_VANLEER) then
      write(6,"(a)") "    - Van Leer's limiter"
    else if(i_limiter.eq.LIMITER_VANALBADA) then
      write(6,"(a)") "    - Van Albada's limiter"
    else if(i_limiter.eq.LIMITER_BARTH_JESPERSEN) then
      write(6,"(a)") "    - Barth-Jespersen limiter"
    else if(i_limiter.eq.LIMITER_VENKATAKRISHNAN) then
      write(6,"(a)") "    - Venkatakrishnan's limiter"
    end if

    return
  end subroutine inputWrite


end module inputParams
