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
  real(dp), dimension(ne) :: vare
  np = Nc-2*Ng
  do m = 1, nplots
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
         yvar(m), np, Cell(NCl:NCu)%Xc, var, ne, xe, vare)
  end do
end subroutine eps_output


!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine tecplot_output
  use solnBlock_module
  use gridBlock_module
  use exactSoln_module
  integer :: n
  !Open the input file.
  open(unit=7,file='output_tecplot.dat',status='REPLACE')
  !Create the output file header:
  write(7,*) 'TITLE = "xCFD1D Solution"'
  write(7,*) 'VARIABLES = "x" \\'
  write(7,*) '"rho" \\'
  write(7,*) '"u" \\'
  write(7,*) '"p" \\'
  write(7,*) '"T" \\'
  write(7,*) '"a" \\'
  write(7,*) '"M" \\'
  write(7,*) '"eps" \\'
  write(7,*) '"E" \\'
  write(7,*) '"H" \\'
  write(7,*) '"s" \\'
  !Output the exact solution if available:
  if(ne.gt.0) then
    write(7,*) 'ZONE T = "Exact Solution" \\'
    write(7,*) 'I = ', ne, ' \\'
    write(7,*) 'F = POINT'
    do n = 1, ne
      write(7,('(11(ES14.6))')) &
           xe(n),  &
           We(n)%rho,    &
           We(n)%u,      &
           We(n)%p,      &
           T_W(We(n)),   &
           a_W(We(n)),   &
           M_W(We(n)),   &
           Esp_W(We(n)), &
           E_W(We(n)),   &
           H_W(We(n)),   &
           s_W(We(n))
    end do
  end if
  !Output the solution data:
  write(7,*) 'ZONE T = "xCFD1D Solution" \\'
  write(7,*) 'I = ', NCu - NCl + 1, ' \\'
  write(7,*) 'F = POINT'
  do n = NCl, NCu
    write(7,('(11(ES14.6))')) &
         Cell(n)%Xc,  &
         W(n)%rho,    &
         W(n)%u,      &
         W(n)%p,      &
         T_W(W(n)),   &
         a_W(W(n)),   &
         M_W(W(n)),   &
         Esp_W(W(n)), &
         E_W(W(n)),   &
         H_W(W(n)),   &
         s_W(W(n))
  end do
  !Close the input file.
  close(7)
  return
end subroutine tecplot_output


!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine tecplot_residual
  use Euler1D_UState
  use realsizes, only: dp
  implicit none
  integer :: ierr
  integer :: n, npts, ns
  real(dp) :: t
  type(Euler1D_U_State) :: l1, l2, mx
  character(len=255) :: buffer
  !Open the input file.
  open(unit=7,file='residual_tecplot.dat',status='REPLACE')
  !Create the output file header:
  write(7,*) 'TITLE = "xCFD1D Residual"'
  write(7,*) 'VARIABLES = "n" \\'
  write(7,*) '"t" \\'
  write(7,*) '"l1-rho" \\'
  write(7,*) '"l1-du" \\'
  write(7,*) '"l1-E" \\'
  write(7,*) '"l2-rho" \\'
  write(7,*) '"l2-du" \\'
  write(7,*) '"l2-E" \\'
  write(7,*) '"mx-rho" \\'
  write(7,*) '"mx-du" \\'
  write(7,*) '"mx-E" \\'
  open(unit=4,file='residual.dat',status='old',action='read',position='rewind')
  read(4,'(a)') buffer !Header line
  ierr = 0
  npts = 0
  do while(ierr.eq.0)
    read(4,'(a)',iostat=ierr) buffer
    if(ierr.eq.0) npts = npts + 1
  end do
  write(7,*) 'ZONE T = "xCFD1D Residual" \\'
  write(7,*) 'I = ', npts, ' \\'
  write(7,*) 'F = POINT'
  rewind(4)
  read(4,'(a)') buffer !Header line
  do n = 1, npts
    read(4,('(a5,10(a11))')) ns, t, l1%rho, l1%du, l1%E, l2%rho, l2%du, l2%E, mx%rho, mx%du, mx%E
    write(7,('(a5,10(a11))')) ns, t, l1%rho, l1%du, l1%E, l2%rho, l2%du, l2%E, mx%rho, mx%du, mx%E
  end do
  !Close the input file.
  close(4)
  close(7)
  return
end subroutine tecplot_residual


!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine gnuplot_output
  use solnBlock_module
  use gridBlock_module
  use exactSoln_module
  integer :: n
  !Write the exact solution file if available:
  if(ne.gt.0) then
    open(unit=7,file='exact_gnuplot.dat',status='REPLACE')
    do n = 1, ne
      write(7,('(11(ES14.6))')) &
           xe(n),  &
           We(n)%rho,    &
           We(n)%u,      &
           We(n)%p,      &
           T_W(We(n)),   &
           a_W(We(n)),   &
           M_W(We(n)),   &
           Esp_W(We(n)), &
           E_W(We(n)),   &
           H_W(We(n)),   &
           s_W(We(n))
    end do
    close(7)
  end if
  !Write the current solution to the data file.
  open(unit=7,file='output_gnuplot.dat',status='REPLACE')
  do n = NCl, NCu
    write(7,('(11(ES14.6))')) &
         Cell(n)%Xc,  &
         W(n)%rho,    &
         W(n)%u,      &
         W(n)%p,      &
         T_W(W(n)),   &
         a_W(W(n)),   &
         M_W(W(n)),   &
         Esp_W(W(n)), &
         E_W(W(n)),   &
         H_W(W(n)),   &
         s_W(W(n))
  end do
  close(7)
  !Create the gnuplot file for the density.
  open(UNIT=7,FILE='density.gplt',STATUS='REPLACE')
  write(7,'(a)') 'set xlabel "x"'
  write(7,'(a)') 'set ylabel "density"'
  write(7,'(a)') 'set grid'
  if(ne.gt.0) then
    write(7,'(a)') 'plot "exact_gnuplot.dat" using 1:2 title "exact" with lines, \'
    write(7,'(a)') '     "output_gnuplot.dat" using 1:2 title "cfd" with lines'
  else
    write(7,'(a)') 'plot "output_gnuplot.dat" using 1:2 title "" with lines'
  end if
  write(7,'(a)') 'pause -1  "Hit return to continue"'
  close(7)
  !Create the gnuplot file for the speed.
  open(UNIT=7,FILE='speed.gplt',STATUS='REPLACE')
  write(7,'(a)') 'set xlabel "x"'
  write(7,'(a)') 'set ylabel "speed"'
  write(7,'(a)') 'set grid'
  if(ne.gt.0) then
    write(7,'(a)') 'plot "exact_gnuplot.dat" using 1:3 title "exact" with lines, \'
    write(7,'(a)') '     "output_gnuplot.dat" using 1:3 title "cfd" with lines'
  else
    write(7,'(a)') 'plot "output_gnuplot.dat" using 1:3 title "" with lines'
  end if
  write(7,'(a)') 'pause -1  "Hit return to continue"'
  close(7)
  !Create the gnuplot file for the pressure.
  open(UNIT=7,FILE='pressure.gplt',STATUS='REPLACE')
  write(7,'(a)') 'set xlabel "x"'
  write(7,'(a)') 'set ylabel "pressure"'
  write(7,'(a)') 'set grid'
  if(ne.gt.0) then
    write(7,'(a)') 'plot "exact_gnuplot.dat" using 1:4 title "exact" with lines, \'
    write(7,'(a)') '     "output_gnuplot.dat" using 1:4 title "cfd" with lines'
  else
    write(7,'(a)') 'plot "output_gnuplot.dat" using 1:4 title "" with lines'
  end if
  write(7,'(a)') 'pause -1  "Hit return to continue"'
  close(7)
  !Create the gnuplot file for the temperature.
  open(UNIT=7,FILE='temperature.gplt',STATUS='REPLACE')
  write(7,'(a)') 'set xlabel "x"'
  write(7,'(a)') 'set ylabel "temperature"'
  write(7,'(a)') 'set grid'
  if(ne.gt.0) then
    write(7,'(a)') 'plot "exact_gnuplot.dat" using 1:5 title "exact" with lines, \'
    write(7,'(a)') '     "output_gnuplot.dat" using 1:5 title "cfd" with lines'
  else
    write(7,'(a)') 'plot "output_gnuplot.dat" using 1:5 title "" with lines'
  end if
  write(7,'(a)') 'pause -1  "Hit return to continue"'
  close(7)
  !Create the gnuplot file for the Mach number.
  open(UNIT=7,FILE='mach.gplt',STATUS='REPLACE')
  write(7,'(a)') 'set xlabel "x"'
  write(7,'(a)') 'set ylabel "Mach Number"'
  write(7,'(a)') 'set grid'
  if(ne.gt.0) then
    write(7,'(a)') 'plot "exact_gnuplot.dat" using 1:7 title "exact" with lines, \'
    write(7,'(a)') '     "output_gnuplot.dat" using 1:7 title "cfd" with lines'
  else
    write(7,'(a)') 'plot "output_gnuplot.dat" using 1:7 title "" with lines'
  end if
  write(7,'(a)') 'pause -1  "Hit return to continue"'
  close(7)
  return
end subroutine gnuplot_output


!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine gnuplot_residual
  integer :: n
  !Create the gnuplot file for the l1-norm of the residual.
  open(UNIT=7,FILE='l1resid.gplt',STATUS='REPLACE')
  write(7,'(a)') 'set xlabel "n"'
  write(7,'(a)') 'set ylabel "l1-norm"'
  write(7,'(a)') 'set logscale y'
  write(7,'(a)') 'set grid'
  write(7,'(a)') 'plot "residual.dat" using 1:3 title "rho" with lines, \'
  write(7,'(a)') '     "residual.dat" using 1:4 title "du" with lines, \'
  write(7,'(a)') '     "residual.dat" using 1:5 title "E" with lines'
  write(7,'(a)') 'pause -1  "Hit return to continue"'
  close(7)
  !Create the gnuplot file for the l2-norm of the residual.
  open(UNIT=7,FILE='l1resid.gplt',STATUS='REPLACE')
  write(7,'(a)') 'set xlabel "n"'
  write(7,'(a)') 'set ylabel "l2-norm"'
  write(7,'(a)') 'set logscale y'
  write(7,'(a)') 'set grid'
  write(7,'(a)') 'plot "residual.dat" using 1:6 title "rho" with lines, \'
  write(7,'(a)') '     "residual.dat" using 1:7 title "du" with lines, \'
  write(7,'(a)') '     "residual.dat" using 1:8 title "E" with lines'
  write(7,'(a)') 'pause -1  "Hit return to continue"'
  close(7)
  !Create the gnuplot file for the maxnorm of the residual.
  open(UNIT=7,FILE='l1resid.gplt',STATUS='REPLACE')
  write(7,'(a)') 'set xlabel "n"'
  write(7,'(a)') 'set ylabel "maxnorm"'
  write(7,'(a)') 'set logscale y'
  write(7,'(a)') 'set grid'
  write(7,'(a)') 'plot "residual.dat" using 1:9 title "rho" with lines, \'
  write(7,'(a)') '     "residual.dat" using 1:10 title "du" with lines, \'
  write(7,'(a)') '     "residual.dat" using 1:11 title "E" with lines'
  write(7,'(a)') 'pause -1  "Hit return to continue"'
  close(7)
  return
end subroutine gnuplot_residual
