!-----------------------------------------------------------------------
! Performs the solution of the one-dimensional Euler equations of 
! compressible gas dynamics for a thermally and calorically perfect gas.
! Spatial disctretization is accomplished using a Godunov-type finite-
! volume scheme and temporal discretization can be achieved using 
! various explicit time-marching schemes.
!-----------------------------------------------------------------------
program xcfd1d
  use inputParams
  use explSolver
  use solnBlock_module
  use gridBlock_module
  use exactSoln_module
  implicit none
  integer            :: ierr
  integer            :: iotag
  character(len=128) :: fname

  write(6,"(a)") '----------------------------------------------------------------------'
  write(6,"(a)") ' xCFD1D SOLVER FOR THE 1D EULER EQUATIONS'
  write(6,"(a)") '----------------------------------------------------------------------'

  call inputDefaults

  fname = 'input.in'
  call inputOpen(fname)

  ierr = 0
  iotag = NOTHING

  do while(iotag.ne.TERMINATE)
    call inputParse(iotag)
    select case(iotag)
    case(NOTHING)
    case(EXECUTE)
      call compute(1)
    case(CONTINUE)
      call compute(0)
    case(PYTHON)
      plot_python = .true.
      plot_gnuplot = .false.
      plot_tecplot = .false.
      plot_eps = .false.
      call plot_data
    case(GNUPLOT)
      plot_python = .false.
      plot_gnuplot = .true.
      plot_tecplot = .false.
      plot_eps = .false.
      call plot_data
    case(TECPLOT)
      plot_python = .false.
      plot_gnuplot = .false.
      plot_tecplot = .true.
      plot_eps = .false.
      call plot_data
    end select
  end do
  
  call inputClose

  write(6,"(a)") '----------------------------------------------------------------------'
  write(6,"(a)") ' xCFD1D SOLVER COMPLETED'
  write(6,"(a)") '----------------------------------------------------------------------'

  !Memory deallocation:  
  write(6,"(a)") ' -> Memory deallocation'
  call gridBlock_deallocate
  call solnBlock_deallocate
  call solnExplicit_deallocate
  call exactSoln_deallocate

  call exit(0)
end program xcfd1d

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
subroutine compute(iopt)
  use realSizes, only: dp
  use inputParams
  use cfdParams
  use Euler1D_UState
  use explSolver
  use solnBlock_module
  use gridBlock_module
  use exactSoln_module
  implicit none
  !Argument variable:
  integer, intent(in) :: iopt
  !Local variables:
  integer               :: ierr      !Error flag
  integer               :: i_stage   !Explicit time-stepping stage number
  integer               :: n_steps   !Number of time-steps
  real(dp)              :: dtime     !Global time-step
  real(dp)              :: time      !Total time
  type(Euler1D_U_State) :: l1norm    !L1-residual norm
  type(Euler1D_U_State) :: l2norm    !L2-residual norm
  type(Euler1D_U_State) :: maxnorm   !Max-residual norm
  character(len=128)    :: fname

  ierr = 0
  
  !Memory deallocation:  
  call gridBlock_deallocate
  call solnBlock_deallocate
  call solnExplicit_deallocate
  call exactSoln_deallocate

  write(6,"(a)") ' -> Checking consistency of the input parameters'
  call inputConsistency(ierr)
  if(ierr.ne.0) call exit(1)

  write(6,"(a)") ' -> Input parameters:'
  call inputWrite

  if(iopt.eq.1) then
    write(6,"(a)") ' -> Memory allocation'
    call gridBlock_allocate(ierr,i_problem,num_cells,num_ghost)
    if(ierr.eq.0) call solnBlock_allocate(ierr,num_cells,num_ghost)
    if(ierr.eq.0) call solnExplicit_allocate
    if(ierr.ne.0) then
      write(6,"(a)") ' ERROR during allocation, error number ', ierr
      call exit(2)
    end if
    write(6,"(a)") ' -> Grid generation'
    call createBlock
    write(6,"(a)") ' -> Setting initial conditions'
    call applyICs
    write(6,"(a)") ' -> Applying boundary conditions'
    call applyBCs
  end if

  call open_residual_file(iopt)

  if(iopt.eq.1) then
    time = 0.
    n_steps = 0
  end if

  !Perform required number of iterations (time steps).
  time_marching: do

    !Consider exit criteria.
    if((i_time_step.eq.GLOBAL_TIME_STEP).and.(time.ge.time_max)) exit
    if((i_time_step.eq.LOCAL_TIME_STEP).and.(n_steps.ge.max_time_steps)) exit
    if((i_time_step.eq.GLOBAL_STEADY_STATE).and.(n_steps.ge.max_time_steps)) exit

    !Determine local and global time steps.
    call setTimeStep(dtime, cfl_number)
    if(i_time_step.eq.GLOBAL_TIME_STEP)then
      if(time + dtime.gt.time_max) then
        dtime = time_max-time
      end if
      call setGlobalTimeStep(dtime)
    end if
    if(i_time_step.eq.GLOBAL_STEADY_STATE) call setGlobalTimeStep(dtime)

    !Update solution for next time step using a multistage time-stepping scheme.
    do i_stage = 1, n_stage
     
      !Apply boundary conditions for stage.
      call applyBCs

      !Compute residual for this stage.
      call multistageExplicit(i_stage, ierr)
      if(ierr.ne.0) call exit(3)

      !Update solution for this stage.
      call updateExplicit(i_stage, ierr)
      if(ierr.ne.0) call exit(4)

    end do

    !Update time and time step counter.
    n_steps = n_steps + 1
    time = time + dtime

    !Determine the L1, L2, and max norms of the solution residual.
    call residual_l1_norm(l1norm)
    call residual_l2_norm(l2norm)
    call residual_max_norm(maxnorm)
    call write_residual_file(n_steps,time,l1norm,l2norm,maxnorm)

    !Output progress information for the calculation.
    if(n_steps-i_output_freq*(n_steps/i_output_freq).eq.0) then
      if(n_steps.lt.1000) then
        write(6,100) n_steps, time, l2norm%rho, l2norm%du, l2norm%E
      else if(n_steps.lt.10000) then
        write(6,101) n_steps, time, l2norm%rho, l2norm%du, l2norm%E
      else
        write(6,102) n_steps, time, l2norm%rho, l2norm%du, l2norm%E
      end if
    end if
100 format(' n = ',I3, ' t = ',F8.6,' l2_norm = ',E12.6,' ',E12.6,' ',E12.6)
101 format(' n = ',I4, ' t = ',F8.6,' l2_norm = ',E12.6,' ',E12.6,' ',E12.6)
102 format(' n = ',I5, ' t = ',F8.6,' l2_norm = ',E12.6,' ',E12.6,' ',E12.6)

  end do time_marching

  !Apply boundary conditions:
  call applyBCs

  !Close residual file:
  call close_residual_file

  return
end subroutine compute
