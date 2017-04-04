!-----------------------------------------------------------------------
! Performs the solution of the one-dimensional Euler equations of 
! compressible gas dynamics for a thermally and calorically perfect gas.
! Spatial disctretization is accomplished using a Godunov-type finite-
! volume scheme and temporal discretization can be achieved using 
! various explicit time-marching schemes.
!-----------------------------------------------------------------------
program xcfd1d
  use realSizes, only: dp
  use inputParams
  use cfdParams
  use Euler1D_UState
  use explSolver
  use solnBlock_module
  use gridBlock_module
  use exactSoln_module

  !Implicit declaration.
  implicit none

  !Local variables:
  integer               :: ierr = 0  !Error flag
  integer               :: i_stage   !Explicit time-stepping stage number
  integer               :: n_steps   !Number of time-steps
  real(dp)              :: dtime     !Global time-step
  real(dp)              :: time      !Total time
  type(Euler1D_U_State) :: l1norm
  type(Euler1D_U_State) :: l2norm
  type(Euler1D_U_State) :: maxnorm

  write(6,"(a)") '----------------------------------------------------------------------'
  write(6,"(a)") ' xCFD1D SOLVER FOR THE 1D EULER EQUATIONS'
  write(6,"(a)") '----------------------------------------------------------------------'
  write(6,"(a)") ' Reading input parameters:'

  write(6,"(a)") ' -> Setting defaults'
  call inputDefaults

  write(6,"(a)") ' -> Reading namelists'
  call inputNamelists

  write(6,"(a)") ' -> Checking consistency of the input parameters'
  call inputConsistency(ierr)
  if(ierr.ne.0) call exit(1)

  write(6,"(a)") ' -> Input parameters:'
  call inputWrite

  write(6,"(a)") '----------------------------------------------------------------------'
  write(6,"(a)") ' Problem setup and initialization:'

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

  call open_residual_file

  write(6,"(a)") '----------------------------------------------------------------------'
  write(6,"(a)") ' Problem solution:'

  time = 0.
  n_steps = 0

  !Perform required number of iterations (time steps).
  time_marching: do

   !Consider exit criteria.
   if(i_time_step.eq.GLOBAL_TIME_STEP.and.time.ge.time_max) exit
   if(i_time_step.eq.LOCAL_TIME_STEP.and.n_steps.ge.max_time_steps) exit

   !Determine local and global time steps.
   call setTimeStep(dtime, cfl_number)
   if(i_time_step.eq.GLOBAL_TIME_STEP) then
     if(time + dtime.gt.time_max) then
       dtime = time_max-time
     end if
     call setGlobalTimeStep(dtime)
   end if

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

  write(6,"(a)") '----------------------------------------------------------------------'
  write(6,"(a)") ' Performing post-processing and memory deallocation:'

  !Close residual file:
  call close_residual_file

  !Output cell-centered data in Tecplot format:
  if(plot_tecplot) then
    write(6,"(a)") ' -> Output cell-centered data for plotting with Tecplot'
    call tecplot_output
    call tecplot_residual
  end if

  !Output cell-centered data in gnuplot format:
  if(plot_gnuplot) then
    write(6,"(a)") ' -> Output cell-centered data for plotting with gnuplot'
    call gnuplot_output
  end if

  !Output plots of cell-centered data directly as eps figures:
  if(plot_eps) then
    write(6,"(a)") ' -> Creating EPS figures with cell-centered data'
    call eps_output
  end if

  !Output cell-centered data in format for python plotting:
  if(plot_python) then
    write(6,"(a)") ' -> Output cell-centered data for plotting with python/matplotlib'
    call python_output
  end if

  !Memory deallocation:  
  write(6,"(a)") ' -> Memory deallocation'
  call gridBlock_deallocate
  call solnBlock_deallocate
  call solnExplicit_deallocate
  call exactSoln_deallocate

  write(6,"(a)") '----------------------------------------------------------------------'
  write(6,"(a)") ' xCFD1D SOLVER COMPLETED'
  write(6,"(a)") '----------------------------------------------------------------------'

  call exit(0)
end program xcfd1d
