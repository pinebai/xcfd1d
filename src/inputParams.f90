!-----------------------------------------------------------------------
! Module for reading input parameters
!-----------------------------------------------------------------------
module inputParams
  use realSizes, only: dp
  use Euler1D_WState
  
  implicit none

  !Define input directives:
  integer, parameter      :: TERMINATE = -1      !Terminates the program
  integer, parameter      :: NOTHING   =  0      !Find the next input
  integer, parameter      :: EXECUTE   =  1      !Start running the program
  integer, parameter      :: CONTINUE  =  2      !Continue running the program after input changes
  integer, parameter      :: PYTHON    =  3      !Write out python plot files
  integer, parameter      :: GNUPLOT   =  4      !Write out gnuplot files
  integer, parameter      :: TECPLOT   =  5      !Write out tecplot files

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
  real(dp)                :: dt_fixed            !Fixed time-step
  namelist / temporal / i_explicit, n_stage, cfl_number, time_max, &
       i_time_step, max_time_steps, dt_fixed

  !Spatial discretization parameters:
  integer                 :: i_flux              !Flux function type
  integer                 :: i_limiter           !Limiter type
  integer                 :: i_recon             !Reconstruction type
  real(dp)                :: betam               !Parameter for gamma-differencing scheme
  namelist / spatial / i_flux, i_limiter, i_recon, betam

  !Output formats:
  integer                 :: i_output_freq       !Frequency of updates on screen
  character(len=64)       :: plot_file           !Name of output file
  logical                 :: plot_tecplot        !Write output to Tecplot format
  logical                 :: plot_gnuplot        !Write output to gnuplot format
  logical                 :: plot_python         !Write output to format ready for python (matplotlib)
  logical                 :: plot_eps            !Make encapsulated post-script plots
  namelist / output / i_output_freq, plot_file, plot_tecplot, plot_gnuplot, plot_python, plot_eps

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
    dt_fixed = 0.0_dp
    i_output_freq = 10
    !Spatial discretization parameters:
    i_flux = FLUX_ROE
    i_limiter = LIMITER_ZERO
    i_recon = RECONSTRUCTION_LSQ
    betam = 0.5_dp
    !Output formats:
    plot_file    = 'output'
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
  ! Open the input file
  !---------------------------------------------------------------------
  subroutine inputOpen(fname)
    implicit none
    character(len=128) :: fname
    open(unit=7,file=trim(fname))
    return
  end subroutine inputOpen

  !---------------------------------------------------------------------
  ! Close the input file
  !---------------------------------------------------------------------
  subroutine inputClose
    close(7)
    return
  end subroutine inputClose

  !---------------------------------------------------------------------
  ! Read input parameters from namelists
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
      write(6,"(a)") "    - Multistage optimal smoothing"
    end if
    if(i_time_step.eq.GLOBAL_TIME_STEP) then
      write(6,"(a)") "    - Global time-stepping"
    else if(i_time_step.eq.LOCAL_TIME_STEP) then
      write(6,"(a)") "    - Local time-stepping"
    else if(i_time_step.eq.GLOBAL_STEADY_STATE) then
      write(6,"(a)") "    - Global steady-state"
    else if(i_time_step.eq.FIXED_TIME_STEP) then
      write(6,"(a)") "    - Fixed time-stepping"
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
    else if(i_flux.eq.FLUX_SLAU) then
      write(6,"(a)") "    - Fluxes calculated using Shima's version of AUSM"
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

  !---------------------------------------------------------------------
  ! Remove equals sign from string
  !---------------------------------------------------------------------
  subroutine remove_char(str,buf)
    implicit none
    character(len=128), intent(inout) :: str
    character(len=1),   intent(in)    :: buf
    integer :: i, ic
    character(len=128) :: strc
    strc = str
    do i = 1, len_trim(strc)
      ic = index(buf,str(i:i))
      if(ic.gt.0) str(i:i) = ' '
    end do
    return
  end subroutine remove_char

  !---------------------------------------------------------------------
  ! Changes string entries to lower case
  !---------------------------------------------------------------------
  subroutine to_lower(str)
    implicit none
    character(len=128), intent(inout):: str
    integer :: i, ic
    character(len=28), parameter :: upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(len=28), parameter :: lower = 'abcdefghijklmnopqrstuvwxyz'
    character(len=128) :: strc
    strc = str
    do i = 1, len_trim(strc)
      ic = index(upper,str(i:i))
      if(ic.gt.0) str(i:i) = lower(ic:ic)
    end do
    return
  end subroutine to_lower

  !---------------------------------------------------------------------
  ! Parse the input file, return the next directive tag
  !---------------------------------------------------------------------
  subroutine inputParse(iotag)
    use cfdParams
    use gasConstants
    implicit none
    !Argument variable:
    integer, intent(out) :: iotag
    !Local variables:
    character(len=128) :: buffer
    iotag = NOTHING
    read(7,'(a)') buffer
    call remove_char(buffer,'=')
    call remove_char(buffer,'-')
    call remove_char(buffer,'_')
    call to_lower(buffer)
    if(trim(buffer).eq.'terminate') then
      iotag = TERMINATE
    else if(trim(buffer).eq.'execute') then
      iotag = EXECUTE
    else if(trim(buffer).eq.'continue') then
      iotag = CONTINUE
    else if(trim(buffer).eq.'plot python') then
      iotag = PYTHON
    else if(trim(buffer).eq.'plot gnuplot') then
      iotag = GNUPLOT
    else if(trim(buffer).eq.'plot tecplot') then
      iotag = TECPLOT
    else if(buffer(1:1).eq.'#') then
      iotag = NOTHING
    else
      iotag = NOTHING
      !--------------------!
      ! Problem parameters !
      !--------------------!
      if(trim(buffer).eq.'problem') then
        read(7,'(a)') buffer
        call to_lower(buffer)
        select case(trim(buffer))
        case('sod problem')
          i_problem = SOD_PROBLEM
        case('modified sod')
          i_problem = MODIFIED_SOD
        case('strong sod')
          i_problem = STRONG_SOD
        case('123 problem')
          i_problem = PROBLEM_123
        case('three right waves')
          i_problem = THREE_RIGHT_WAVES
        case('stationary contact')
          i_problem = STATIONARY_CONTACT
        case('subsonic nozzle')
          i_problem = SUBSONIC_NOZZLE
        case('transonic nozzle')
          i_problem = TRANSONIC_NOZZLE
        case('square wave')
          i_problem = SQUARE_WAVE
        case('sine squared wave')
          i_problem = SINE_SQUARED_WAVE
        end select
      else if(trim(buffer).eq.'wave width') then
        read(7,*) wm
      else if(trim(buffer).eq.'wave speed') then
        read(7,*) um
      else if(trim(buffer).eq.'wave density') then
        read(7,*) rhom
      else if(trim(buffer).eq.'wave pressure') then
        read(7,*) pm
      else if(trim(buffer).eq.'number of cells') then
        read(7,*) num_cells
      else if(trim(buffer).eq.'gas type') then
        read(7,'(a)') buffer
        call to_lower(buffer)
        select case(trim(buffer))
        case('air')
          i_gas_type = GAS_AIR
        case('he')
          i_gas_type = GAS_HE
        case('h2')
          i_gas_type = GAS_H2
        case('n2')
          i_gas_type = GAS_N2
        case('o2')
          i_gas_type = GAS_O2
        case('helium')
          i_gas_type = GAS_HE
        case('hydrogen')
          i_gas_type = GAS_H2
        case('nitrogen')
          i_gas_type = GAS_N2
        case('oxygen')
          i_gas_type = GAS_O2
        end select
      end if
      !-------------------------!
      ! Temporal discretization !
      !-------------------------!
      if(trim(buffer).eq.'time stepping scheme') then
        read(7,'(a)') buffer
        call remove_char(buffer,'-')
        call remove_char(buffer,'_')
        call to_lower(buffer)
        select case(trim(buffer))
        case('explicit euler')
          i_explicit = EXPLICIT_EULER
        case('predictor corrector')
          i_explicit = PREDICTOR_CORRECTOR
        case('runge kutta')
          i_explicit = RUNGE_KUTTA
        case('multistage optimal smoothing')
          i_explicit = MULTISTAGE_OPTIMAL_SMOOTHING
        end select
      else if(trim(buffer).eq.'number of stages') then
        read(7,*) n_stage
      else if(trim(buffer).eq.'cfl number') then
        read(7,*) cfl_number
      else if(trim(buffer).eq.'time stepping') then
        read(7,'(a)') buffer
        call remove_char(buffer,'-')
        call remove_char(buffer,'_')
        call to_lower(buffer)
        select case(trim(buffer))
        case('global')
          i_time_step = GLOBAL_TIME_STEP
        case('local')
          i_time_step = LOCAL_TIME_STEP
        case('global steady state')
          i_time_step = GLOBAL_STEADY_STATE
        case('fixed')
          i_time_step = FIXED_TIME_STEP
        end select
      else if(trim(buffer).eq.'maximum time') then
        read(7,*) time_max
      else if((trim(buffer).eq.'max time steps').or. &
              (trim(buffer).eq.'maximum time steps')) then
        read(7,*) max_time_steps
      else if(trim(buffer).eq.'fixed time step') then
        read(7,*) dt_fixed
      end if
      !-------------------------!
      ! Spatial discretization  !
      !-------------------------!
      if(trim(buffer).eq.'flux scheme') then
        read(7,'(a)') buffer
        call remove_char(buffer,'-')
        call remove_char(buffer,'_')
        call to_lower(buffer)
        select case(trim(buffer))
        case('godunov')
          i_flux = FLUX_GODUNOV
        case('isentropic')
          i_flux = FLUX_ISENTROPIC
        case('rusanov')
          i_flux = FLUX_RUSANOV
        case('hlle')
          i_flux = FLUX_HLLE
        case('hlll')
          i_flux = FLUX_HLLL
        case('hllc')
          i_flux = FLUX_HLLC
        case('osher')
          i_flux = FLUX_OSHER
        case('roe')
          i_flux = FLUX_ROE
        case('van leer')
          i_flux = FLUX_VANLEER
        case('ausm plus')
          i_flux = FLUX_AUSMplus
        case('ausm plus up')
          i_flux = FLUX_AUSMplusup
        case('slau')
          i_flux = FLUX_SLAU
        end select
      else if(trim(buffer).eq.'gradient reconstruction') then
        read(7,'(a)') buffer
        call remove_char(buffer,'-')
        call remove_char(buffer,'_')
        call to_lower(buffer)
        select case(trim(buffer))
        case('green gauss')
          i_recon = RECONSTRUCTION_GG
        case('least squares')
          i_recon = RECONSTRUCTION_LSQ
        case('weno')
          i_recon = RECONSTRUCTION_WENO
        case('ppm')
          i_recon = RECONSTRUCTION_PPM
        case('gamma')
          i_recon = RECONSTRUCTION_GAMMA
        end select
      else if(trim(buffer).eq.'limiter') then
        read(7,'(a)') buffer
        call remove_char(buffer,'-')
        call remove_char(buffer,'_')
        call to_lower(buffer)
        select case(trim(buffer))
        case('zero')
          i_limiter = LIMITER_ZERO
        case('one')
          i_limiter = LIMITER_ONE
        case('minmod')
          i_limiter = LIMITER_MINMOD
        case('umist')
          i_limiter = LIMITER_UMIST
        case('double minmod')
          i_limiter = LIMITER_DOUBLE_MINMOD
        case('superbee')
          i_limiter = LIMITER_SUPERBEE
        case('phi')
          i_limiter = LIMITER_PHI
        case('van leer')
          i_limiter = LIMITER_VANLEER
        case('van albada')
          i_limiter = LIMITER_VANALBADA
        case('sine')
          i_limiter = LIMITER_SINE
        case('barth jespersen')
          i_limiter = LIMITER_BARTH_JESPERSEN
        case('venkatakrisnan')
          i_limiter = LIMITER_VENKATAKRISHNAN
        end select
      end if
      !--------------!
      ! Input/Output !
      !--------------!
      if(trim(buffer).eq.'output frequency') then
        read(7,*) i_output_freq
      else if(trim(buffer).eq.'file name') then
        read(7,'(a)') plot_file
      else if(trim(buffer).eq.'plot format') then
        read(7,'(a)') buffer
        call remove_char(buffer,'-')
        call remove_char(buffer,'_')
        call to_lower(buffer)
        select case(trim(buffer))
        case('python')
          plot_python = .true.
        case('gnuplot')
          plot_gnuplot = .true.
        case('tecplot')
          plot_tecplot = .true.
        case('eps')
          plot_eps = .true.
        end select
      end if

    end if
    return
  end subroutine inputParse

end module inputParams
