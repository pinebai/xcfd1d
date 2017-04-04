!-----------------------------------------------------------------------
! Module for reading input parameters
!-----------------------------------------------------------------------
module inputParams
  use realSizes, only: dp
  use Euler1D_WState
  
  implicit none

  !Problem related input parameters:
  integer                 :: i_problem           !Problem specification
  integer                 :: num_cells           !Number of cells in the domain
  integer                 :: num_ghost           !Number of ghost cells
  integer                 :: i_gas_type          !Gas type
  real(dp)                :: rhom, um, pm        !Reference state for wave problems
  real(dp)                :: wm                  !Wave width
  namelist / problem / i_problem, num_cells, num_ghost, i_gas_type, &
       rhom, um, pm, wm

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
  character(len=64)       :: plot_name           !Name of output file
  logical                 :: plot_tecplot        !Write output to Tecplot format
  logical                 :: plot_gnuplot        !Write output to gnuplot format
  logical                 :: plot_python         !Write output to format ready for python (matplotlib)
  logical                 :: plot_eps            !Make encapsulated post-script plots
  namelist / output / plot_name, plot_tecplot, plot_gnuplot, plot_python, plot_eps

  !EPS Figures:
  integer,  parameter     :: max_plots = 20
  integer                 :: nplots
  real(dp), dimension(20) :: xf, yf, xo, yo, xa, ya
  real(dp), dimension(20) :: xmin, xmax, dx, xdim, ymin, ymax, dy, ydim
  integer,  dimension(20) :: xvar, yvar, xsig, ysig
  namelist / eps_plots / nplots, xf, yf, xo, yo, xa, ya, &
       xvar, xmin, xmax, dx, yvar, ymin, ymax, dy, xdim, ydim, xsig, ysig

contains

  !---------------------------------------------------------------------
  ! Set input parameter default values
  !---------------------------------------------------------------------
  subroutine inputDefaults
    use cfdParams
    use gasConstants, only: GAS_AIR, g, xmw, setGas
    implicit none
    !Problem related input parameters:
    i_problem = SOD_PROBLEM
    num_cells = 20
    num_ghost = 2
    i_gas_type = GAS_AIR
    call setGas(i_gas_type)
    rhom = 1.225_dp
    um = 10.0_dp
    pm = 101325.0_dp
    wm = 0.1_dp
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
    plot_name    = 'output'
    plot_tecplot = .false.
    plot_gnuplot = .false.
    plot_eps     = .false.
    plot_python  = .false.
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
  ! Read in input parameters from namelists
  !---------------------------------------------------------------------
  subroutine inputNamelists
    use gasConstants, only: setGas
    implicit none
    !Open the input file
    open(unit=7,file='input.in')
    !Read in the problem namelist:
    read(unit=7,nml=problem)
    !Read in the temporal discretization namelist:
    read(unit=7,nml=temporal)
    !Read in the spatial discretization namelist:
    read(unit=7,nml=spatial)
    !Read in the output format namelist:
    read(unit=7,nml=output)
    !Read in the EPS plots namelist:
    read(unit=7,nml=eps_plots)
    !Set the gas types:
    call setGas(i_gas_type)
    !Close the input file
    close(7)
    return
  end subroutine inputNamelists

  !---------------------------------------------------------------------
  ! Check that input parameters are consistent and error free
  !---------------------------------------------------------------------
  subroutine inputConsistency(ierr)
    use Numbers, only: TOLER
    use cfdParams
    use gasConstants, only: GAS_AIR, GAS_O2, g, xmw, setGas
    implicit none
    integer, intent(out) :: ierr
    ierr = 0

    if(i_problem.lt.SOD_PROBLEM.or.i_problem.gt.SINE_SQUARED_WAVE) then
      write(6,"(a)") "    - Initial condition type specified incorrectly."
      ierr = ierr + 1
    end if

    if(num_cells.lt.5) then
      write(6,"(a)") "    - Too few cells have been requested."
      ierr = ierr + 1
    end if

    if(num_ghost.lt.2.or.num_ghost.gt.3) then
      write(6,"(a)") "    - An incorrect number of ghost cells have been specified."
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
  ! Write the input parameters to standard output
  !---------------------------------------------------------------------
  subroutine inputWrite
    use cfdParams
    implicit none
    if(i_problem.eq.SOD_PROBLEM) then
      write(6,"(a)") "    - Sod problem"
    else if(i_problem.eq.MODIFIED_SOD) then
      write(6,"(a)") "    - Modified Sod problem"
    else if(i_problem.eq.STRONG_SOD) then
      write(6,"(a)") "    - Strong Sod"
    else if(i_problem.eq.PROBLEM_123) then
      write(6,"(a)") "    - 123 problem"
    else if(i_problem.eq.THREE_RIGHT_WAVES) then
      write(6,"(a)") "    - Three right travelling waves"
    else if(i_problem.eq.STATIONARY_CONTACT) then
      write(6,"(a)") "    - Stationary contact surface"
    else if(i_problem.eq.SUBSONIC_NOZZLE) then
      write(6,"(a)") "    - Subsonic nozzle flow"
    else if(i_problem.eq.TRANSONIC_NOZZLE) then
      write(6,"(a)") "    - Transonic nozzle flow"
    else if(i_problem.eq.SQUARE_WAVE) then
      write(6,"(a)") "    - Square density wave"
    else if(i_problem.eq.SINE_SQUARED_WAVE) then
      write(6,"(a)") "    - Sine squared density wave"
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
