!-----------------------------------------------------------------------
! Module for a one-dimensional solution block.
!-----------------------------------------------------------------------
module solnBlock_module
  use Euler1D_UState  !Conserved-variable solution state
  use Euler1D_WState  !Primitive-variable solution state (rho,u,p)

  implicit none

  !Solution block variable declaration.
  type(Euler1D_U_State), allocatable, dimension(:) :: U    !Array of conserved variables
  type(Euler1D_W_State), allocatable, dimension(:) :: W    !Array of (rho,u,p)-primitive variables
  type(Euler1D_W_State), allocatable, dimension(:) :: dWdx !Array of unlimited solution gradient
  type(Euler1D_W_State), allocatable, dimension(:) :: phi  !Array of solution limiters
  type(Euler1D_W_State)                            :: WoL  !Left boundary reference states
  type(Euler1D_W_State)                            :: WoR  !Right boundary reference states
  integer                                          :: NCl  !Lower cell domain index
  integer                                          :: NCu  !Upper cell domain index
  integer                                          :: NC   !Number of cells
  integer                                          :: Ng   !Number of ghost cells (per side)
  integer                                          :: BCl  !Left boundary condition
  integer                                          :: BCr  !Right boundary condition
  type(Euler1D_W_State), allocatable, dimension(:) :: Ql   !Left interface state
  type(Euler1D_W_State), allocatable, dimension(:) :: Qr   !Right interface state

contains

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  subroutine solnBlock_allocate(status,num_cells,num_ghost)
    use cfdParams, only: BC_CONSTANT_EXTRAPOLATION
    implicit none
    integer, intent(out) :: status    !Allocation status flag
    integer, intent(in)  :: num_cells !Number of cells
    integer, intent(in)  :: num_ghost !Number of ghost cells
    integer              :: n
    !Set the cell counters:
    Ng  = num_ghost
    NC  = num_cells + 2*Ng
    NCl = Ng + 1
    NCu = Nc - Ng
    !Set the boundary conditions:
    BCl = BC_CONSTANT_EXTRAPOLATION
    BCr = BC_CONSTANT_EXTRAPOLATION
    !Allocate memory for the array of conserved variables:
    allocate(U(NC),STAT=status)
    if(status.ne.0) return
    !Allocate memory for the array of (rho,u,p)-primitive variables:
    allocate(W(NC),STAT=status)
    if(status.ne.0) return
    !Allocate memory for the array of unlimited solution gradient:
    allocate(dWdx(NC),STAT=status)
    if(status.ne.0) return
    !Allocate memory for the array of solution limiters:
    allocate(phi(NC),STAT=status)
    if(status.ne.0) return
    !Allocate memory for the left/right reconstructed states:
    allocate(Ql(NC),STAT=status)
    if(status.ne.0) return
    allocate(Qr(NC),STAT=status)
    if(status.ne.0) return
    do n = NCl-Ng, NCu+Ng
     call standard_atmosphere_U(U(n))
     call standard_atmosphere_W(W(n))
     call vacuum_W(dWdx(n))
     call vacuum_W(phi(n))
     call vacuum_W(Ql(n))
     call vacuum_W(Qr(n))
    end do
  end subroutine solnBlock_allocate

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  subroutine solnBlock_deallocate
    implicit none
    if(allocated(U)) deallocate(U)
    if(allocated(W)) deallocate(W)
    if(allocated(dWdx)) deallocate(dWdx)
    if(allocated(phi)) deallocate(phi)
    if(allocated(Ql)) deallocate(Ql)
    if(allocated(Qr)) deallocate(Qr)
  end subroutine solnBlock_deallocate

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  subroutine applyICs
    use realSizes, only: dp
    use numbers, only: PI
    use cfdParams
    use inputParams, only: i_problem, time_max, rhom, um, pm, wm
    use gridBlock_module
    use exactSoln_module
    implicit none
    integer               :: n
    type(Euler1D_W_State) :: Wl, Wr
    select case(i_problem)
    case(SOD_PROBLEM)
      Wl = WState(1.0_dp,0.0_dp,1.0_dp)
      Wr = WState(0.125_dp,0.0_dp,0.1_dp)
      do n = NCl-Ng, NCu+Ng
        if(Cell(n)%Xc.le.0.5_dp) then
          W(n) = Wl
        else
          W(n) = Wr
        end if
      end do
      call exactSod(time_max)
    case(MODIFIED_SOD)
      Wl = WState(1.0_dp,0.75_dp,1.0_dp)
      Wr = WState(0.125_dp,0.0_dp,0.1_dp)
      do n = NCl-Ng, NCu+Ng
        if(Cell(n)%Xc.le.0.3_dp) then
          W(n) = Wl
        else
          W(n) = Wr
        end if
      end do
      call exactModifiedSod(time_max)
    case(STRONG_SOD)
      Wl = WState(1.0_dp,0.0_dp,1000.0_dp)
      Wr = WState(0.125_dp,0.0_dp,0.01_dp)
      do n = NCl-Ng, NCu+Ng
        if(Cell(n)%Xc.le.0.5_dp) then
          W(n) = Wl
        else
          W(n) = Wr
        end if
      end do
      call exactStrongSod(time_max)
    case(PROBLEM_123)
      Wl = WState(1.0_dp,-2.0_dp,0.4_dp)
      Wr = WState(1.0_dp, 2.0_dp,0.4_dp)
      do n = NCl-Ng, NCu+Ng
        if(Cell(n)%Xc.le.0.5_dp) then
          W(n) = Wl
        else
          W(n) = Wr
        end if
      end do
      call exact123Problem(time_max)
    case(THREE_RIGHT_WAVES)
      Wl = WState(5.99924_dp,19.5975_dp,460.894_dp)
      Wr = WState(5.99242_dp,-6.19633_dp,46.095_dp)
      do n = NCl-Ng, NCu+Ng
        if(Cell(n)%Xc.le.0.4_dp) then
          W(n) = Wl
        else
          W(n) = Wr
        end if
      end do
      call exact3RightWaves(time_max)
    case(STATIONARY_CONTACT)
      Wl = WState(1.0_dp,-19.59745_dp,1000.0_dp)
      Wr = WState(1.0_dp,-19.59745_dp,0.01_dp)
      do n = NCl-Ng, NCu+Ng
        if(Cell(n)%Xc.le.0.8_dp) then
          W(n) = Wl
        else
          W(n) = Wr
        end if
      end do
      call exactStationaryContact(time_max)
    case(SUBSONIC_NOZZLE)
      Wl = WState(1.1409_dp,65.451_dp,97534.0_dp)
      Wr = WState(1.1008_dp,113.06_dp,92772.0_dp)
      do nc = NCl-Ng, NCu+Ng
        if(Cell(nc)%Xc.le.5.0_dp) then
          W(nc) = Wl
        else
          W(nc) = Wr
        end if
        !call Nozzle(Cell(nc)%Xc,Cell(nc)%Xa,IC_SUBSONIC_NOZZLE,We(nc))
      end do
      BCl = BC_FIXED
      BCr = BC_CONSTANT_EXTRAPOLATION
    case(TRANSONIC_NOZZLE)
      Wl = WState(1.1288_dp,82.693_dp,96085.0_dp)
      Wr = WState(1.0261_dp,151.62_dp,84974.0_dp)
      do nc = NCl-Ng, NCu+Ng
        if(Cell(nc)%Xc.le.5.0_dp) then
          W(nc) = Wl
        else
          W(nc) = Wr
        end if
        !call Nozzle(Cell(nc)%Xc,Cell(nc)%Xa,IC_TRANSONIC_NOZZLE,We(nc))
      end do
      BCl = BC_FIXED
      BCr = BC_CONSTANT_EXTRAPOLATION
    case(SQUARE_WAVE)
      Wl%rho = 1.0_dp*rhom ; Wl%u = um ; Wl%p = pm
      Wr%rho = 2.0_dp*rhom ; Wr%u = um ; Wr%p = pm
      do nc = NCl-Ng, NCu+Ng
        if((Cell(nc)%Xc.le.-0.5_dp*wm).or.(Cell(nc)%Xc.ge.0.5_dp)) then
          W(nc) = Wl
        else
          W(nc) = Wr
        end if
        !call exactSquareWave(Cell(nc)%Xc,time_max,We(nc))
      end do
      BCl = BC_PERIODIC
      BCr = BC_PERIODIC
    case(SINE_SQUARED_WAVE)
      Wl%rho = rhom ; Wl%u = um ; Wl%p = pm
      do nc = NCl-Ng, NCu+Ng
        if(6.0_dp*Cell(nc)%Xc.lt.1.0_dp) then
          W(nc) = Wl
        else if(Cell(nc)%Xc.le.0.5_dp) then
          W(nc) = Wl
          W(nc)%rho = W(nc)%rho*(1.0_dp + (sin(0.5_dp*PI*(6.0_dp*Cell(nc)%Xc-1.0_dp))**2))
        else
          W(nc) = Wl
        end if
        !call Sine_Squared_Wave(Cell(nc)%Xc,time_max,We(nc))
      end do
      BCl = BC_PERIODIC
      BCr = BC_PERIODIC
    case(SEMI_ELLIPSE_WAVE)
      Wl%rho = 1.225_dp ; Wl%u = 100.0_dp ; Wl%p = 101325.0_dp
      do nc = NCl-Ng, NCu+Ng
        if(6.0_dp*Cell(nc)%Xc.lt.1.0_dp) then
          W(nc) = Wl
        else if(Cell(nc)%Xc.le.0.50_dp) then
          W(nc) = Wl
          W(nc)%rho = W(nc)%rho*(1.0_dp + sqrt(1.0_dp-(6.0_dp*Cell(nc)%Xc-2.0_dp)**2))
        else
          W(nc) = Wl
        end if
        !call Semi_Ellipse_Wave(Cell(nc)%Xc,time_max,We(nc))
      end do
      BCl = BC_PERIODIC
      BCr = BC_PERIODIC
    end select
    do n = NCl-Ng, NCu+Ng
      call transform_W_to_U(W(n),U(n))
    end do
    return
  end subroutine applyICs

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  subroutine applyBCs
    use cfdParams
    implicit none
    integer :: n
    !Apply left boundary condition.
    select case(BCl)
    case(BC_NONE)
    case(BC_FIXED)
    case(BC_CONSTANT_EXTRAPOLATION)
      do n = NCl-1, NCl-Ng
        W(n) = W(NCl)
        U(n) = U(NCl)
      end do
    case(BC_REFLECTION)
      do n = 1, Ng
        call reflect_W(W(NCl-n),W(NCl+n-1))
        call transform_W_to_U(W(NCl-n),U(NCl-n))
      end do
    case(BC_PERIODIC)
      do n = 1, Ng
        W(NCl-n) = W(NCu-n+1)
        U(NCl-n) = U(NCu-n+1)
      end do
    end select
    !Apply right boundary condition.
    select case(BCr)
    case(BC_NONE)
    case(BC_FIXED)
    case(BC_CONSTANT_EXTRAPOLATION)
      do n = NCu+1, NCu+Ng
        W(n) = W(NCu)
        U(n) = U(NCu)
      end do
    case(BC_REFLECTION)
      do n = 1, Ng
        call reflect_W(W(NCu+n),W(NCu-n+1))
        call transform_W_to_U(W(NCu+n),U(NCu+n))
      end do
    case(BC_PERIODIC)
      do n = 1, Ng
        W(NCu+n) = W(NCl+n-1)
        U(NCu+n) = U(NCl+n-1)
      end do
    end select
    return
  end subroutine applyBCs

end module solnBlock_module
