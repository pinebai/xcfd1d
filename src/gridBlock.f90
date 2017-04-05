!-----------------------------------------------------------------------
! Module for a one-dimensional grid block.
!-----------------------------------------------------------------------
module gridBlock_module
  !Include real sizes module
  use realSizes, only: dp

  !Implicit declaration
  implicit none

  !Solution block variable declaration
  type nodeGeom
    real(dp) :: x !Node location
  end type nodeGeom

  !Solution block variable declaration
  type cellGeom
    real(dp) :: xc !Cell center
    real(dp) :: dx !Cell length
    real(dp) :: xa !Cross-sectional area at cell center
    real(dp) :: da !Change in the cross-sectional area
  end type cellGeom

  !Solution block variable declaration
  type gridData
    integer :: NN  !Total number of nodes
    integer :: NC  !Total number of cells
    integer :: Ng  !Number of ghost cells (per side)
    integer :: NNl !Lower node domain index
    integer :: NNu !Upper node domain index
    integer :: NCl !Lower cell domain index
    integer :: NCu !Upper cell domain index
  end type gridData

  type(nodeGeom), allocatable, dimension(:) :: Node !Array of node locations
  type(cellGeom), allocatable, dimension(:) :: Cell !Array of cell-centre locations
  type(gridData)                            :: Grid

  integer :: i_grid
  
contains

  !---------------------------------------------------------------------
  ! Computes the value of the various algebraic stretching functions
  ! used to control grid point distributions [Anderson, Tannehill, and
  ! Pletcher, 1984].
  !---------------------------------------------------------------------
  real(dp) function StretchingFcn(x, beta, tau, str_fcn)
    use numbers, only: PI
    use cfdParams
    implicit none
    real(dp), intent(in) :: x         ! Normalized distance
    real(dp), intent(in) :: beta, tau ! Stretching parameters
    integer,  intent(in) :: str_fcn   ! Stretching function
    real(dp)             :: a
    select case(str_fcn)
    case(STR_LINEAR)
      StretchingFcn = x
    case(STR_MIN)
      a = ((beta+1.0_dp)/(beta-1.0_dp))**(1.0_dp-x)
      StretchingFcn = ((beta+1.0_dp)-(beta-1.0_dp)*a)/(a+1.0_dp)
    case(STR_MAX)
      a = ((beta+1.0_dp)/(beta-1.0_dp))**x
      StretchingFcn = beta*(a - 1.0_dp)/(a+1.0_dp)
    case(STR_MINMAX)
      a = ((beta+1.0_dp)/(beta-1.0_dp))**((x-0.5_dp)/0.5_dp)
      StretchingFcn = (a*(beta+1.0_dp) - beta + 1.0_dp)/(2.0_dp*(a+1.0_dp))
    case(STR_MIDPT)
      a = 0.5_dp*log((1.0_dp+(exp(tau)-1.0_dp)*beta)/(1.0_dp+(exp(-tau)-1.0_dp)*beta))/tau
      StretchingFcn = beta*(1.0_dp + sinh(tau*(x-a))/sinh(tau*a))
    case(STR_SINE)
      StretchingFcn = sin(0.5_dp*PI*x)
    case(STR_COSINE)
      StretchingFcn = 1.0_dp - cos(0.5_dp*PI*x)
    end select
    return
  end function StretchingFcn

  !---------------------------------------------------------------------
  ! Allocate grid node and cell arrays.
  !---------------------------------------------------------------------
  subroutine gridBlock_allocate(status,i_problem,num_cells,num_ghost)
    use cfdParams
    implicit none
    integer, intent(out) :: status    !Allocation status flag
    integer, intent(in)  :: i_problem !Problem specification
    integer, intent(in)  :: num_cells !Number of cells
    integer, intent(in)  :: num_ghost !Number of ghost cells
    integer :: n
    i_grid = GRID_SHOCK_TUBE
    if((i_problem.eq.SUBSONIC_NOZZLE).or.(i_problem.eq.TRANSONIC_NOZZLE)) i_grid = GRID_NOZZLE
    if((i_problem.ge.SQUARE_WAVE).and.(i_problem.le.SINE_SQUARED_WAVE)) i_grid = GRID_WAVE_PROBLEM
    Grid%Ng  = num_ghost
    Grid%NC  = num_cells + 2*Grid%Ng
    Grid%NN  = Grid%NC + 1
    Grid%NNl = Grid%Ng + 1
    Grid%NNu = Grid%NN - Grid%Ng
    Grid%NCl = Grid%Ng + 1
    Grid%NCu = Grid%NC - Grid%Ng
    allocate(Node(Grid%NN),STAT=status)
    if(status.ne.0) return
    allocate(Cell(Grid%NC),STAT=status)
    if(status.ne.0) return
    do n = 1, Grid%NN
      Node(n)%x = 0.0_dp
    end do
    do n = 1, Grid%NC
      Cell(n)%xc = 0.0_dp
      Cell(n)%dx = 0.0_dp
      Cell(n)%xa = 1.0_dp
      Cell(n)%da = 0.0_dp
    end do
    return
  end subroutine gridBlock_allocate

  !---------------------------------------------------------------------
  ! Allocate grid node and cell arrays
  !---------------------------------------------------------------------
  subroutine gridBlock_deallocate
    implicit none
    if(allocated(Node)) deallocate(Node)
    if(allocated(Cell)) deallocate(Cell)
    return
  end subroutine gridBlock_deallocate

  !---------------------------------------------------------------------
  ! Return the cell centroid
  !---------------------------------------------------------------------
  real(dp) function centroid(nc)
    implicit none
    integer, intent(in) :: nc
    centroid = 0.5_dp*(Node(nc)%X + Node(nc+1)%X)
    return
  end function centroid

  !---------------------------------------------------------------------
  ! Return the cell length
  !---------------------------------------------------------------------
  real(dp) function length(nc)
    implicit none
    integer, intent(in) :: nc
    length = Node(nc+1)%X - Node(nc)%X
    return
  end function length

  !---------------------------------------------------------------------
  ! Return the position of the left face of the cell
  !---------------------------------------------------------------------
  real(dp) function xfaceL(nc)
    implicit none
    integer, intent(in) :: nc
    xfaceL = Node(nc)%X
    return
  end function xfaceL

  !---------------------------------------------------------------------
  ! Return the position of the right face of the cell
  !---------------------------------------------------------------------
  real(dp) function xfaceR(nc)
    implicit none
    integer, intent(in) :: nc
    xfaceR = Node(nc+1)%X
    return
  end function xfaceR

  !---------------------------------------------------------------------
  ! Update the exterior nodes of the grid block
  !---------------------------------------------------------------------
  subroutine updateNodes(BCl,BCr)
    use cfdParams
    implicit none
    integer, intent(in) :: BCl
    integer, intent(in) :: BCr
    integer :: nn
    !Update left-side exterior nodes:
    select case(BCl)
    case(BC_NONE, BC_FIXED, BC_CONSTANT_EXTRAPOLATION, BC_CHARACTERISTIC)
      do nn = 1, Grid%Ng
        Node(Grid%NNl-nn)%X = Node(Grid%NNl)%X + nn*(Node(Grid%NNl)%X - Node(Grid%NNl+1)%X)
      end do
    case(BC_REFLECTION)
      do nn = 1, Grid%Ng
        Node(Grid%NNl-nn)%X = 2.*Node(Grid%NNl)%X - Node(Grid%NNl+nn)%X
      end do
    case(BC_PERIODIC)
      do nn = 1, Grid%Ng
        Node(Grid%NNl-nn)%X = Node(Grid%NNl)%X - (Node(Grid%NNu)%X-Node(Grid%NNu-nn)%X)
      end do
    end select
    !Update right-side exterior nodes:
    select case(BCr)
    case(BC_NONE, BC_FIXED, BC_CONSTANT_EXTRAPOLATION, BC_CHARACTERISTIC)
      do nn = 1, Grid%Ng
        Node(Grid%NNu+nn)%X = Node(Grid%NNu)%X + nn*(Node(Grid%NNu)%X - Node(Grid%NNu-1)%X)
      end do
    case(BC_REFLECTION)
      do nn = 1, Grid%Ng
        Node(Grid%NNu+nn)%X = 2.*Node(Grid%NNu)%X - Node(Grid%NNu-nn)%X
      end do
    case(BC_PERIODIC)
      do nn = 1, Grid%Ng
        Node(Grid%NNu+nn)%X = Node(Grid%NNu)%X + (Node(Grid%NNl)%X-Node(Grid%NNl-nn)%X)
      end do
    end select
    return
  end subroutine updateNodes

  !---------------------------------------------------------------------
  ! Update the cells of the grid block
  !---------------------------------------------------------------------
  subroutine updateCells
    implicit none
    integer :: nc
    do nc = 1, Grid%NC
      Cell(nc)%xc = centroid(nc)
      Cell(nc)%dx = length(nc)
    end do
    return
  end subroutine updateCells

  !---------------------------------------------------------------------
  ! Create the grid block
  !---------------------------------------------------------------------
  subroutine createBlock
    use cfdParams
    implicit none
    integer  :: nn
    real(dp) :: x, xl, xr
    integer  :: bcl, bcr
    integer  :: str_fcn
    real(dp) :: beta, tau
    !Set parameters for problem:
    select case(i_grid)
    case(GRID_SHOCK_TUBE)
      xl = 0.0_dp ; bcl = BC_REFLECTION
      xr = 1.0_dp ; bcr = BC_REFLECTION
      beta = 1.0_dp ; tau = 1.0_dp
      str_fcn = STR_LINEAR
    case(GRID_NOZZLE)
      xl =  0.0_dp ; bcl = BC_FIXED
      xr = 10.0_dp ; bcr = BC_CONSTANT_EXTRAPOLATION
      beta = 1.0_dp ; tau = 1.0_dp
      str_fcn = STR_MIDPT
    case(GRID_WAVE_PROBLEM)
      xl = -0.5_dp ; bcl = BC_PERIODIC
      xr =  0.5_dp ; bcr = BC_PERIODIC
      beta = 1.0_dp ; tau = 1.0_dp
      str_fcn = STR_LINEAR
    end select
    !Create the internal nodes of the block:
    do nn = Grid%NNl, Grid%NNu
      x = real(nn-Grid%NNl)/real(Grid%NNu-Grid%NNl)
      Node(nn)%x = xl + (xr-xl)*StretchingFcn(x,beta,tau,str_fcn)
    end do
    !Update the exterior nodes:
    call updateNodes(bcl,bcr)
    !Update the cell centroids and widths:
    call updateCells
    !Set cross-sectional area:
    call set_xarea
    return
  end subroutine createBlock

  !---------------------------------------------------------------------
  ! Set the cross sectional area
  !---------------------------------------------------------------------
  subroutine set_xarea
    use cfdParams
    implicit none
    integer :: nc
    select case(i_grid)
    case(GRID_SHOCK_TUBE,GRID_WAVE_PROBLEM)
      do nc = 1, Grid%NC
        Cell(nc)%xA = 1.0_dp
        Cell(nc)%dA = 0.0_dp
      end do
    case(GRID_NOZZLE)
      do nc = 1, Grid%NC
        if(Cell(nc)%Xc.lt.5.0_dp) then
          Cell(nc)%xA = 1.0_dp + 1.5_dp*(1.0_dp - Cell(nc)%Xc/5.0_dp)**2
          Cell(nc)%dA =        - 0.6_dp*(1.0_dp - Cell(nc)%Xc/5.0_dp)
        else
          Cell(nc)%xA = 1.0_dp + 0.5_dp*(1.0_dp - Cell(nc)%Xc/5.0_dp)**2
          Cell(nc)%dA =        - 0.2_dp*(1.0_dp - Cell(nc)%Xc/5.0_dp)
        end if
      end do
    end select
    return
  end subroutine set_xarea


end module gridBlock_module
