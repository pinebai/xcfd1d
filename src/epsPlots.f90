!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
module epsPlots_module
  !Include real sizes module.
  use realSizes, only: dp

  !Implicit declaration.
  implicit none

  integer, parameter :: PLOT_DENSITY = 1
  integer, parameter :: PLOT_SPEED = 2
  integer, parameter :: PLOT_PRESSURE = 3
  integer, parameter :: PLOT_TEMPERATURE = 4
  integer, parameter :: PLOT_INTERNAL_ENERGY = 5
  integer, parameter :: PLOT_MACH_NUMBER = 6
  integer, parameter :: PLOT_CROSS_SECTIONAL_AREA = 7
  integer, parameter :: PLOT_CHANGE_IN_CROSS_SECTIONAL_AREA = 8
  integer, parameter :: PLOT_L1_NORM_RESIDUAL = 11
  integer, parameter :: PLOT_L2_NORM_RESIDUAL = 12
  integer, parameter :: PLOT_MAX_NORM_RESIDUAL = 13

contains

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  real(dp) function x_intersection(xa1, ya1, xa3, ya3, xb1, yb1, xb3, yb3)
    implicit none
    real(dp), intent(in) :: xa1, ya1, xa3, ya3
    real(dp), intent(in) :: xb1, yb1, xb3, yb3
    real(dp)             :: xa2, ya2, xb2, yb2
    real(dp)             :: det, s, t
    xa2 = (xa3 - xa1)
    ya2 = (ya3 - ya1)
    xb2 = (xb3 - xb1)
    yb2 = (yb3 - yb1)
    det = ya2*xb2 - xa2*yb2
    s = ((yb1 - ya1)*xb2 - (xb1 - xa1)*yb2)/det
    t = ((yb1 - ya1)*xa2 - (xb1 - xa1)*ya2)/det
    x_intersection = xa1 + s*xa2
    return
  end function x_intersection


  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  subroutine profile_plot(Xf, Yf, Xo, Yo, Xa, Ya, xmin, xmax, dx, xdim, &
       ymin, ymax, dy, ydim, xsig, ysig, m, n, x, y, ne, xe, ye)
    use inputParams, only: plot_file
    implicit none
    !Argument variables:
    real(dp), intent(in) :: Xf, Yf, Xo, Yo, Xa, Ya
    real(dp), intent(in) :: xmin, xmax, dx, xdim
    real(dp), intent(in) :: ymin, ymax, dy, ydim
    integer,  intent(in) :: xsig, ysig
    integer,  intent(in) :: m, n, ne
    real(dp), dimension(n), intent(in) :: x, y
    real(dp), dimension(n), intent(in) :: xe, ye
    !Local variables:
    real(dp) :: xx, yy
    integer  :: nn, nx, ny
    integer  :: i, j
    character(len=28) :: fname
    if(m.eq.PLOT_DENSITY) then
      fname = trim(plot_file)//'_density.eps'
    else if(m.eq.PLOT_SPEED) then
      fname = trim(plot_file)//'_speed.eps'
    else if(m.eq.PLOT_PRESSURE) then
      fname = trim(plot_file)//'_pressure.eps'
    else if(m.eq.PLOT_TEMPERATURE) then
      fname = trim(plot_file)//'_temperature.eps'
    else if(m.eq.PLOT_INTERNAL_ENERGY) then
      fname = trim(plot_file)//'_internal_energy.eps'
    else if(m.eq.PLOT_MACH_NUMBER) then
      fname = trim(plot_file)//'_mach_number.eps'
    else if(m.eq.PLOT_CROSS_SECTIONAL_AREA) then
      fname = trim(plot_file)//'_area.eps'
    else
      fname = trim(plot_file)//'_darea.eps'
    end if
    open(unit=7,file=trim(fname),status='replace')
    write(7,'(a)') "%!EPS-Adobe-3.0"
    write(7,'(2(a,i6),2(a,f6.1))') "%%BoundingBox: ", -5, " ", -5, " ", Xf+5, " ", Yf+5
    write(7,'(a)') "%%Magnification: 1.0000"
    write(7,'(a)') "%%EndComments"
    write(7,'(a)') ""
    write(7,'(a)') "%----------------------------------------------------------------------%";
    write(7,'(a)') "%---------------------------- BEGIN FIGURE ----------------------------%";
    write(7,'(a)') "%----------------------------------------------------------------------%";
    write(7,'(a)') ""
    write(7,'(a)') ""
    !Draw box:
    !write(7,'(a)') "newpath"
    !write(7,*)  0, " ",  0, " moveto"
    !write(7,*) Xf, " ",  0, " lineto"
    !write(7,*) Xf, " ", Yf, " lineto"
    !write(7,*)  0, " ", Yf, " lineto"
    !write(7,*)  0, " ",  0, " lineto 0.1 setlinewidth stroke"
    !X-axis:
    !write(7,'(a)') "newpath"
    !write(7,*) Xo, " ", Yo, " moveto"
    !write(7,*) Xf, " ", Yo, " lineto 0.75 setlinewidth stroke"
    !Y-axis:
    !write(7,'(a)') "newpath"
    !write(7,*) Xo, " ", Yo, " moveto "
    !write(7,*) Xo, " ", Yf, " lineto 0.75 setlinewidth stroke"
    write(7,'(a)') "newpath"
    write(7,*) Xo, " ", Yo, " moveto"
    write(7,*) Xf, " ", Yo, " lineto"
    write(7,*) Xf, " ", Yf, " lineto"
    write(7,*) Xo, " ", Yf, " lineto"
    write(7,*) Xo, " ", Yo, " lineto 0.1 setlinewidth stroke"
    !Vertical dashed lines.
    xx = xmin + dx
    do
      if(Xo+xx*(Xf-Xo)/(xmax-xmin).lt.Xf) then
        write(7,'(a)') "newpath "
        write(7,*) Xo+xx*(Xf-Xo)/(xmax-xmin), " ", Yo, " moveto"
        write(7,*) Xo+xx*(Xf-Xo)/(xmax-xmin), " ", Yf, " lineto 0.25 setlinewidth [3 3] 0 setdash stroke"
      end if
      xx = xx + dx
      if(xx.gt.xmax) exit
    end do
    !Horizontal dashed lines.
    yy = ymin + dy
    do
      if(Yo+(yy-ymin)*(Yf-Yo)/(ymax-ymin).lt.Yf) then
        write(7,'(a)') "newpath "
        write(7,*) Xo, " ", Yo+(yy-ymin)*(Yf-Yo)/(ymax-ymin), " moveto"
        write(7,*) Xf, " ", Yo+(yy-ymin)*(Yf-Yo)/(ymax-ymin), " lineto 0.25 setlinewidth [3 3] 0 setdash stroke"
      end if
      yy = yy + dy
      if(yy.gt.ymax) exit
    end do
    !Labels:
    write(7,'(a)') "/Helvetica findfont 16 scalefont setfont"
    xx = xmin + dx
    do 
      nx = -18
      if(xx.ge.xmax-0.0001) nx = -28
      if(xsig.eq.0) then
        write(7,"(2(F12.4,a),I4,a)") Xo+(xx-xmin)*(Xf-Xo)/(xmax-xmin)+nx, " ", Yo-Ya, " moveto (", int(xx), ") show"
      else if(xsig.eq.1) then
        write(7,"(2(F12.4,a),F4.1,a)") Xo+(xx-xmin)*(Xf-Xo)/(xmax-xmin)+nx, " ", Yo-Ya, " moveto (", xx, ") show"
      else if(xsig.eq.2) then
        write(7,"(2(F12.4,a),F4.2,a)") Xo+(xx-xmin)*(Xf-Xo)/(xmax-xmin)+nx, " ", Yo-Ya, " moveto (", xx, ") show"
      else 
        write(7,"(2(F12.4,a),F4.3,a)") Xo+(xx-xmin)*(Xf-Xo)/(xmax-xmin)+nx, " ", Yo-Ya, " moveto (", xx, ") show"
      end if
      xx = xx + dx
      if(xx.gt.xmax) exit
    end do
    write(7,*) (Xo+Xf/2)-4*16, " ", 4, " moveto (X-Coordinate (m)) show"
    yy = ymin
    do
      !if(yy.lt.ymin+0.5*dy) then
      !  ny = 2
      !else 
      if(yy.gt.ymax-0.0001) then
        ny = -14
      else
        ny = -6
      end if
      if(ysig.eq.0) then
        write(7,"(2(F12.4,a),I5,a)") Xo-Xa, " ", Yo+(yy-ymin)*(Yf-Yo)/(ymax-ymin)+ny, " moveto (", int(yy), ") show"
      else if(ysig.eq.1) then
        write(7,"(2(F12.4,a),F5.1,a)") Xo-Xa, " ", Yo+(yy-ymin)*(Yf-Yo)/(ymax-ymin)+ny, " moveto (", yy, ") show"
      else if(ysig.eq.2) then
        write(7,"(2(F12.4,a),F5.2,a)") Xo-Xa, " ", Yo+(yy-ymin)*(Yf-Yo)/(ymax-ymin)+ny, " moveto (", yy, ") show"
      else if(ysig.eq.3) then
        write(7,"(2(F12.4,a),F5.3,a)") Xo-Xa, " ", Yo+(yy-ymin)*(Yf-Yo)/(ymax-ymin)+ny, " moveto (", yy, ") show"
      else
        write(7,"(2(F12.4,a),F6.4,a)") Xo-Xa, " ", Yo+(yy-ymin)*(Yf-Yo)/(ymax-ymin)+ny, " moveto (", yy, ") show"
      end if
      yy = yy + dy
      if(yy.gt.ymax) exit
    end do
    write(7,'(a)') "90 rotate"
    if(m.eq.PLOT_DENSITY) then
      write(7,'(a)') "160 -15 moveto (Density (kg/m^3)) show"
    else if(m.eq.PLOT_SPEED) then
      write(7,'(a)') "176 -15 moveto (Speed (m/s)) show"
    else if(m.eq.PLOT_PRESSURE) then
      if(ydim.eq.1.) then
        write(7,'(a)') "168 -15 moveto (Pressure (Pa)) show"
      else if(ydim.eq.1.0d3) then
        write(7,'(a)') "168 -15 moveto (Pressure (kPa)) show"
      else
        write(7,'(a)') "160 -15 moveto (Pressure (MPa)) show"
      end if
    else if(m.eq.PLOT_TEMPERATURE) then
      write(7,'(a)') "160 -15 moveto (Temperature (K)) show"
    else if(m.eq.PLOT_INTERNAL_ENERGY) then
      if(ydim.eq.1.) then
        write(7,'(a)') "160 -15 moveto (Internal Energy (J)) show"
      else if(ydim.eq.1.0d3) then
        write(7,'(a)') "160 -15 moveto (Internal Energy (kJ)) show"
      else if(ydim.eq.1.0d6) then
        write(7,'(a)') "160 -15 moveto (Internal Energy (MJ)) show"
      else
        write(7,'(a)') "160 -15 moveto (Internal Energy (J)) show"
      end if
    else if(m.eq.PLOT_MACH_NUMBER) then
      write(7,'(a)') "176 -15 moveto (Mach Number) show"
    else if(m.eq.PLOT_CROSS_SECTIONAL_AREA) then
      write(7,'(a)') "120 -15 moveto (Cross-Section Area (m^2)) show"
    else
      write(7,'(a)') "108 -15 moveto (Change in Cross-Section Area (m)) show"
    end if
    write(7,'(a)') "-90 rotate"
    write(7,'(a)') " "
    write(7,'(a)') "newpath"
    do nn = 1, ne
      if(nn.eq.1) then
        write(7,*) Xo+((xe(nn)-xmin)/xdim)*(Xf-Xo)/(xmax-xmin), " ", Yo+((ye(nn)-ymin)/ydim)*(Yf-Yo)/(ymax-ymin), " moveto"
      else
        write(7,*) Xo+((xe(nn)-xmin)/xdim)*(Xf-Xo)/(xmax-xmin), " ", Yo+((ye(nn)-ymin)/ydim)*(Yf-Yo)/(ymax-ymin), " lineto"
      end if
    end do
    write(7,'(a)') " 1 setlinewidth [12 0] 1 setdash 0 0 0 setrgbcolor stroke"
    write(7,'(a)') " "
    write(7,'(a)') "newpath"
    do nn = 1, N
      if(nn.eq.1) then
        write(7,*) Xo+((x(nn)-xmin)/xdim)*(Xf-Xo)/(xmax-xmin), " ", Yo+((y(nn)-ymin)/ydim)*(Yf-Yo)/(ymax-ymin), " moveto"
      else
        write(7,*) Xo+((x(nn)-xmin)/xdim)*(Xf-Xo)/(xmax-xmin), " ", Yo+((y(nn)-ymin)/ydim)*(Yf-Yo)/(ymax-ymin), " lineto"
      end if
    end do
    write(7,'(a)') " 1 setlinewidth [12 0] 1 setdash 0 0 1 setrgbcolor stroke"
    write(7,'(a)') " "
    write(7,'(a)') " "
    write(7,'(a)') "%% DISPLAY THE FIGURE:"
    write(7,'(a)') "showpage"
    write(7,'(a)') "%----------------------------------------------------------------------%"
    write(7,'(a)') "%----------------------------- end FIGURE -----------------------------%"
    write(7,'(a)') "%----------------------------------------------------------------------%"
    !Close the eps file:
    close(7)
    return
  end subroutine profile_plot


  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  subroutine residual_plot(Xf, Yf, Xo, Yo, Xa, Ya, xmin, xmax, dx, &
       ymin, ymax, dy, m, n, y)
    use inputParams, only: plot_file
    implicit none
    real(dp), intent(in) :: Xf, Yf, Xo, Yo, Xa, Ya
    real(dp), intent(in) :: xmin, xmax, dx, ymin, ymax, dy
    integer, intent(in) :: m, n
    real(dp), dimension(n,3), intent(in) :: y
    real(dp) :: xx, yy
    integer  :: nn, nx, ny, nv
    character(len=28) :: fname
    if(m.eq.PLOT_L1_NORM_RESIDUAL) then
      fname = trim(plot_file)//'_l1norm.eps'
    else if(m.eq.PLOT_L2_NORM_RESIDUAL) then
      fname = trim(plot_file)//'_l2norm.eps'
    else
      fname = trim(plot_file)//'_mxnorm.eps'
    end if
    open(unit=7,file=trim(fname),status='replace')
    write(7,'(a)') "%!EPS-Adobe-3.0"
    write(7,'(2(a,i6),2(a,f6.1))') "%%BoundingBox: ", -5, " ", -5, " ", Xf+5, " ", Yf+5
    write(7,'(a)') "%%Magnification: 1.0000"
    write(7,'(a)') "%%EndComments"
    write(7,'(a)') ""
    write(7,'(a)') "%----------------------------------------------------------------------%";
    write(7,'(a)') "%---------------------------- BEGIN FIGURE ----------------------------%";
    write(7,'(a)') "%----------------------------------------------------------------------%";
    write(7,'(a)') ""
    write(7,'(a)') ""
    !Draw box:
    write(7,'(a)') "newpath"
    write(7,*)  0, " ",  0, " moveto"
    write(7,*) Xf, " ",  0, " lineto"
    write(7,*) Xf, " ", Yf, " lineto"
    write(7,*)  0, " ", Yf, " lineto"
    write(7,*)  0, " ",  0, " lineto 0.1 setlinewidth stroke"
    !X-axis:
    write(7,'(a)') "newpath"
    write(7,*)  0, " ", Yo, " moveto"
    write(7,*) Xf, " ", Yo, " lineto 0.75 setlinewidth stroke"
    !Y-axis:
    write(7,'(a)') "newpath"
    write(7,*) Xo, " ",  0, " moveto "
    write(7,*) Xo, " ", Yf, " lineto 0.75 setlinewidth stroke"
    !Vertical dashed lines.
    xx = xmin + dx
    do
      if(Xo+xx*(Xf-Xo)/(xmax-xmin).lt.Xf) then
        write(7,'(a)') "newpath "
        write(7,*) Xo+xx*(Xf-Xo)/(xmax-xmin), " ", Yo, " moveto"
        write(7,*) Xo+xx*(Xf-Xo)/(xmax-xmin), " ", Yf, " lineto 0.25 setlinewidth [3 3] 0 setdash stroke"
      end if
      xx = xx + dx
      if(xx.gt.xmax) exit
    end do
    !Horizontal dashed lines.
    yy = ymin + dy
    do
      if(Yo+(yy-ymin)*(Yf-Yo)/(ymax-ymin).lt.Yf) then
        write(7,'(a)') "newpath "
        write(7,*) Xo, " ", Yo+(yy-ymin)*(Yf-Yo)/(ymax-ymin), " moveto"
        write(7,*) Xf, " ", Yo+(yy-ymin)*(Yf-Yo)/(ymax-ymin), " lineto 0.25 setlinewidth [3 3] 0 setdash stroke"
      end if
      yy = yy + dy
      if(yy.gt.ymax) exit
    end do
    !Labels:
    write(7,'(a)') "/Helvetica findfont 16 scalefont setfont"
    !do xx = xmin+dx, xmax+dx, dx
    xx = xmin + dx
    do
      nx = -20
      if(xx .ge. xmax-0.0001) nx = -36
      write(7,"(2(F12.4,a),I4,a)") Xo+(xx-xmin)*(Xf-Xo)/(xmax-xmin)+nx, " ", Yo-Ya, " moveto (", int(xx), ") show"
      xx = xx + dx
      if(xx.gt.xmax+dx) return
    end do
    write(7,'(a)') (Xo+Xf/2)-4*16, " ", 4, " moveto (X-Coordinate (m)) show"
    !do yy = ymin, ymax+dy, dy
    yy =  ymin
    do
      if(yy.lt.ymin+0.5*dy) then
        ny = 2
      else if(yy.lt.ymax) then
        ny = -6
      else
        ny = -14
      end if
      write(7,'(a)') "/Helvetica findfont 16 scalefont setfont"
      !write(7,"(2(F12.4,a),I3,a)") Xo-Xa, " ", Yo+(yy-ymin)*(Yf-Yo)/(ymax-ymin)+ny, " moveto (10", int(yy), ") show"
      write(7,"(2(F12.4,a))") Xo-Xa, " ", Yo+(yy-ymin)*(Yf-Yo)/(ymax-ymin)+ny, " moveto (10) show"
      write(7,'(a)') "/Helvetica findfont 10 scalefont setfont"
      if(int(yy) .ge. 0) then
        write(7,"(2(F12.4,a),I3,a)") Xo-Xa+14, " ", Yo+(yy-ymin)*(Yf-Yo)/(ymax-ymin)+ny+6, " moveto (", int(yy), ") show"
      else if(int(yy) .ge. -9) then
        write(7,"(2(F12.4,a),I3,a)") Xo-Xa+16, " ", Yo+(yy-ymin)*(Yf-Yo)/(ymax-ymin)+ny+6, " moveto (", int(yy), ") show"
      else
        write(7,"(2(F12.4,a),I3,a)") Xo-Xa+18, " ", Yo+(yy-ymin)*(Yf-Yo)/(ymax-ymin)+ny+6, " moveto (", int(yy), ") show"
      end if
      yy = yy + dy
      if(yy.gt.ymax+dy) exit
    end do
    write(7,'(a)') "/Helvetica findfont 16 scalefont setfont"
    write(7,'(a)') "90 rotate"
    if(m.eq.PLOT_L1_NORM_RESIDUAL) then
      write(7,'(a)') "160 -15 moveto (L1-Norm Residual) show"
    else if(m.eq.PLOT_L2_NORM_RESIDUAL) then
      write(7,'(a)') "160 -15 moveto (L2-Norm Residual) show"
    else
      write(7,'(a)') "160 -15 moveto (Max-Norm Residual) show"
    end if
    write(7,'(a)') "-90 rotate"
    write(7,'(a)') " "
    do nv = 1, 3
      write(7,'(a)') "newpath"
      do nn = 1, N
        if(nn.eq.1) then
          write(7,*) Xo+(nn-xmin)*(Xf-Xo)/(xmax-xmin), " ", Yo+(log10(y(nn,nv))-ymin)*(Yf-Yo)/(ymax-ymin), " moveto"
        else
          write(7,*) Xo+(nn-xmin)*(Xf-Xo)/(xmax-xmin), " ", Yo+(log10(y(nn,nv))-ymin)*(Yf-Yo)/(ymax-ymin), " lineto"
        end if
      end do
      if(nv.eq.1) then
        write(7,'(a)') " 1 setlinewidth [12 0] 1 setdash 1 0 0 setrgbcolor stroke"
      else if(nv.eq.2) then
        write(7,'(a)') " 1 setlinewidth [12 0] 1 setdash 0 1 0 setrgbcolor stroke"
      else
        write(7,'(a)') " 1 setlinewidth [12 0] 1 setdash 0 0 1 setrgbcolor stroke"
      end if
      write(7,'(a)') " "
    end do
    write(7,'(a)') " "
    write(7,'(a)') "%% DISPLAY THE FIGURE:"
    write(7,'(a)') "showpage"
    write(7,'(a)') "%----------------------------------------------------------------------%"
    write(7,'(a)') "%----------------------------- end FIGURE -----------------------------%"
    write(7,'(a)') "%----------------------------------------------------------------------%"
    !Close the eps file:
    close(7)
  end subroutine residual_plot

end module epsPlots_module
