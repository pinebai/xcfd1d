!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine plot_data
  use inputParams, only: plot_python, plot_gnuplot, plot_tecplot, plot_eps
  implicit none
  
  !Output cell-centered data in format for python plotting:
  if(plot_python) then
    write(6,"(a)") ' -> Output cell-centered data for plotting with python/matplotlib'
    call python_output
  end if

  !Output cell-centered data in gnuplot format:
  if(plot_gnuplot) then
    write(6,"(a)") ' -> Output cell-centered data for plotting with gnuplot'
    call gnuplot_output
  end if

  !Output cell-centered data in Tecplot format:
  if(plot_tecplot) then
    write(6,"(a)") ' -> Output cell-centered data for plotting with Tecplot'
    call tecplot_output
    call tecplot_residual
  end if

  !Output plots of cell-centered data directly as eps figures:
  if(plot_eps) then
    write(6,"(a)") ' -> Output cell-centered data in EPS format'
    call eps_output
  end if

  return
end subroutine plot_data

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine python_output
  use cfdParams
  use inputParams, only: plot_file
  use solnBlock_module
  use gridBlock_module
  use exactSoln_module
  implicit none
  integer :: n, nvar
  character(len=128) :: fname
  fname = trim(plot_file)//'.dat'
  !Open the input file.
  open(unit=8,file=trim(fname),status='REPLACE')
  !Create the output file header:
  nvar = 11
  if(i_grid.eq.GRID_NOZZLE) nvar = 12
  write(8,*) nvar
  write(8,*) 'Position, x, [m]'
  if(i_grid.eq.GRID_NOZZLE) write(8,*) 'Area, xA, [m2]'
  write(8,*) 'Density, rho, [kg/m3]'
  write(8,*) 'Speed, u, [m/s]'
  write(8,*) 'Pressure, p, [Pa]'
  write(8,*) 'Temperature, T, [K]'
  write(8,*) 'Sound speed, a, [m/s]'
  write(8,*) 'Mach number, M, [-]'
  write(8,*) 'Internal Energy, eps, [J/kg]'
  write(8,*) 'Total Energy, E, [J]'
  write(8,*) 'Total Enthalpy, H, [J]'
  write(8,*) 'Entropy, s, [J/kg]'
  !Output the exact solution if available:
  write(8,*) ne
  do n = 1, ne
    write(8,'(es20.12,$)') Xe(n)%xc
    if(i_grid.eq.GRID_NOZZLE) write(8,'(es20.12,$)') Xe(n)%xa
    write(8,'(es20.12,$)') We(n)%rho
    write(8,'(es20.12,$)') We(n)%u
    write(8,'(es20.12,$)') We(n)%p
    write(8,'(es20.12,$)') T_W(We(n))
    write(8,'(es20.12,$)') a_W(We(n))
    write(8,'(es20.12,$)') M_W(We(n))
    write(8,'(es20.12,$)') Esp_W(We(n))
    write(8,'(es20.12,$)') E_W(We(n))
    write(8,'(es20.12,$)') H_W(We(n))
    write(8,'(es20.12)') s_W(We(n))
  end do
  !Output the solution data:
  write(8,*) NCu - NCl + 1
  do n = NCl, NCu
    write(8,'(es20.12,$)') Cell(n)%xc
    if(i_grid.eq.GRID_NOZZLE) write(8,'(es20.12,$)') Cell(n)%xa
    write(8,'(es20.12,$)') W(n)%rho
    write(8,'(es20.12,$)') W(n)%u
    write(8,'(es20.12,$)') W(n)%p
    write(8,'(es20.12,$)') T_W(W(n))
    write(8,'(es20.12,$)') a_W(W(n))
    write(8,'(es20.12,$)') M_W(W(n))
    write(8,'(es20.12,$)') Esp_W(W(n))
    write(8,'(es20.12,$)') E_W(W(n))
    write(8,'(es20.12,$)') H_W(W(n))
    write(8,'(es20.12)') s_W(W(n))
  end do
  !Close the output file.
  close(8)
  return
end subroutine python_output


!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine eps_output
  use realSizes, only: dp
  use solnBlock_module
  use gridBlock_module
  use exactSoln_module
  use inputParams, only: nplots, Xf, Yf, Xo, Yo, Xa, Ya, &
       xvar, xmin, xmax, dx, yvar, ymin, ymax, dy, xdim, ydim, xsig, ysig
  use epsPlots_module
  implicit none
  integer :: m, n, np
  real(dp), dimension(Nc-2*Ng) :: var
  real(dp), dimension(ne) :: loce
  real(dp), dimension(ne) :: vare
  np = Nc-2*Ng
  do m = 1, nplots
    loce(1:ne) = Xe(1:ne)%xc
    if(yvar(m).eq.PLOT_DENSITY) then
      var(1:np) = W(NCl:NCu)%rho
      vare(1:ne) = We(1:ne)%rho
    else if(yvar(m).eq.PLOT_SPEED) then
      var(1:np) = W(NCl:NCu)%u
      vare(1:ne) = We(1:ne)%u
    else if(yvar(m).eq.PLOT_PRESSURE) then
      var(1:np) = W(NCl:NCu)%p
      vare(1:ne) = We(1:ne)%p
    else if(yvar(m).eq.PLOT_TEMPERATURE) then
      do n = NCl, NCu 
        var(n-Ng) = T_W(W(n))
      end do
      do n = 1, ne
        vare(n) = T_W(We(n))
      end do
    else if(yvar(m).eq.PLOT_INTERNAL_ENERGY) then
      do n = NCl, NCu 
        var(n-Ng) = eps_W(W(n))
      end do
      do n = 1, ne
        vare(n) = eps_W(We(n))
      end do
    else if(yvar(m).eq.PLOT_MACH_NUMBER) then
      do n = NCl, NCu 
        var(n-Ng) = M_W(W(n))
      end do
      do n = 1, ne
        var(n) = M_W(We(n))
      end do
    else
      write(6,*) "  ~~> ERROR: Profile ", m, " is an invalid selection for EPS plotting."
    end if
    call profile_plot(Xf(m), Yf(m), Xo(m), Yo(m), Xa(m), Ya(m), &
         xmin(m), xmax(m), dx(m), xdim(m), &
         ymin(m), ymax(m), dy(m), ydim(m), xsig(m), ysig(m), &
         yvar(m), np, Cell(NCl:NCu)%xc, var, ne, loce, vare)
  end do
  return
end subroutine eps_output

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine tecplot_output
  use cfdParams
  use inputParams, only: plot_file
  use solnBlock_module
  use gridBlock_module
  use exactSoln_module
  implicit none
  integer :: n
  character(len=128) :: fname
  fname = trim(plot_file)//'.dat'
  !Open the input file.
  open(unit=8,file=trim(fname),status='REPLACE')
  !Create the output file header:
  write(8,*) 'TITLE = "xCFD1D Solution"'
  write(8,*) 'VARIABLES = "x" \\'
  if(i_grid.eq.GRID_NOZZLE) write(8,*) '"xA" \\'
  write(8,*) '"rho" \\'
  write(8,*) '"u" \\'
  write(8,*) '"p" \\'
  write(8,*) '"T" \\'
  write(8,*) '"a" \\'
  write(8,*) '"M" \\'
  write(8,*) '"eps" \\'
  write(8,*) '"E" \\'
  write(8,*) '"H" \\'
  write(8,*) '"s" \\'
  !Output the exact solution if available:
  if(ne.gt.0) then
    write(8,*) 'ZONE T = "Exact Solution" \\'
    write(8,*) 'I = ', ne, ' \\'
    write(8,*) 'F = POINT'
    do n = 1, ne
      write(8,'(es20.12,$)') Xe(n)%xc
      if(i_grid.eq.GRID_NOZZLE) write(8,'(es20.12,$)') Xe(n)%xa
      write(8,'(es20.12,$)') We(n)%rho
      write(8,'(es20.12,$)') We(n)%u
      write(8,'(es20.12,$)') We(n)%p
      write(8,'(es20.12,$)') T_W(We(n))
      write(8,'(es20.12,$)') a_W(We(n))
      write(8,'(es20.12,$)') M_W(We(n))
      write(8,'(es20.12,$)') Esp_W(We(n))
      write(8,'(es20.12,$)') E_W(We(n))
      write(8,'(es20.12,$)') H_W(We(n))
      write(8,'(es20.12)') s_W(We(n))
    end do
  end if
  !Output the solution data:
  write(8,*) 'ZONE T = "xCFD1D Solution" \\'
  write(8,*) 'I = ', NCu - NCl + 1, ' \\'
  write(8,*) 'F = POINT'
  do n = NCl, NCu
    write(8,'(es20.12,$)') Cell(n)%xc
    if(i_grid.eq.GRID_NOZZLE) write(8,'(es20.12,$)') Cell(n)%xa
    write(8,'(es20.12,$)') W(n)%rho
    write(8,'(es20.12,$)') W(n)%u
    write(8,'(es20.12,$)') W(n)%p
    write(8,'(es20.12,$)') T_W(W(n))
    write(8,'(es20.12,$)') a_W(W(n))
    write(8,'(es20.12,$)') M_W(W(n))
    write(8,'(es20.12,$)') Esp_W(W(n))
    write(8,'(es20.12,$)') E_W(W(n))
    write(8,'(es20.12,$)') H_W(W(n))
    write(8,'(es20.12)') s_W(W(n))
  end do
  !Close the input file.
  close(8)
  return
end subroutine tecplot_output

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine tecplot_residual
  use inputParams, only: plot_file
  use Euler1D_UState
  use realsizes, only: dp
  implicit none
  integer :: ierr
  integer :: n, npts, ns
  real(dp) :: t
  type(Euler1D_U_State) :: l1, l2, mx
  character(len=255) :: buffer
  character(len=128) :: fname
  fname = trim(plot_file)//'_residual.dat'
  !Open the input file.
  open(unit=8,file=trim(fname),status='REPLACE')
  !Create the output file header:
  write(8,*) 'TITLE = "xCFD1D Residual"'
  write(8,*) 'VARIABLES = "n" \\'
  write(8,*) '"t" \\'
  write(8,*) '"l1-rho" \\'
  write(8,*) '"l1-du" \\'
  write(8,*) '"l1-E" \\'
  write(8,*) '"l2-rho" \\'
  write(8,*) '"l2-du" \\'
  write(8,*) '"l2-E" \\'
  write(8,*) '"mx-rho" \\'
  write(8,*) '"mx-du" \\'
  write(8,*) '"mx-E" \\'
  fname = trim(plot_file)//'.rsd'
  open(unit=4,file=trim(fname),status='old',action='read',position='rewind')
  read(4,'(a)') buffer !Header line
  ierr = 0
  npts = 0
  do while(ierr.eq.0)
    read(4,'(a)',iostat=ierr) buffer
    if(ierr.eq.0) npts = npts + 1
  end do
  write(8,*) 'ZONE T = "xCFD1D Residual" \\'
  write(8,*) 'I = ', npts, ' \\'
  write(8,*) 'F = POINT'
  rewind(4)
  read(4,'(a)') buffer !Header line
  do n = 1, npts
    read(4,'(a5,10(es11.4))') ns, t, l1%rho, l1%du, l1%E, l2%rho, l2%du, l2%E, mx%rho, mx%du, mx%E
    write(8,'(a5,10(es11.4))') ns, t, l1%rho, l1%du, l1%E, l2%rho, l2%du, l2%E, mx%rho, mx%du, mx%E
  end do
  !Close the input file.
  close(4)
  close(8)
  return
end subroutine tecplot_residual

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine gnuplot_output
  use cfdParams
  use inputParams, only: plot_file
  use solnBlock_module
  use gridBlock_module
  use exactSoln_module
  implicit none
  integer :: n
  character(len=128) :: fname
  !Write the exact solution file if available:
  if(ne.gt.0) then
    fname = trim(plot_file)//'_exact.dat'
    open(unit=8,file=trim(fname),status='REPLACE')
    do n = 1, ne
      write(8,'(es20.12,$)') Xe(n)%xc
      if(i_grid.eq.GRID_NOZZLE) write(8,'(es20.12,$)') Xe(n)%xa
      write(8,'(es20.12,$)') We(n)%rho
      write(8,'(es20.12,$)') We(n)%u
      write(8,'(es20.12,$)') We(n)%p
      write(8,'(es20.12,$)') T_W(We(n))
      write(8,'(es20.12,$)') a_W(We(n))
      write(8,'(es20.12,$)') M_W(We(n))
      write(8,'(es20.12,$)') Esp_W(We(n))
      write(8,'(es20.12,$)') E_W(We(n))
      write(8,'(es20.12,$)') H_W(We(n))
      write(8,'(es20.12)') s_W(We(n))
    end do
    close(8)
  end if
  !Write the current solution to the data file.
  fname = trim(plot_file)//'.dat'
  open(unit=8,file='output_gnuplot.dat',status='REPLACE')
  do n = NCl, NCu
    write(8,'(es20.12,$)') Cell(n)%xc
    if(i_grid.eq.GRID_NOZZLE) write(8,'(es20.12,$)') Cell(n)%xa
    write(8,'(es20.12,$)') W(n)%rho
    write(8,'(es20.12,$)') W(n)%u
    write(8,'(es20.12,$)') W(n)%p
    write(8,'(es20.12,$)') T_W(W(n))
    write(8,'(es20.12,$)') a_W(W(n))
    write(8,'(es20.12,$)') M_W(W(n))
    write(8,'(es20.12,$)') Esp_W(W(n))
    write(8,'(es20.12,$)') E_W(W(n))
    write(8,'(es20.12,$)') H_W(W(n))
    write(8,'(es20.12)') s_W(W(n))
  end do
  close(8)
  return
end subroutine gnuplot_output
