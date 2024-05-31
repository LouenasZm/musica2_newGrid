! ==============================================================================
module mod_interp1
! ==============================================================================
  !> Module for 1D interpolation
! ==============================================================================
  use precision
  use warnstop
  implicit none

contains

  !===============================================================
  subroutine interp1(xi,fi,nxi,x,f,nx,method)
  !===============================================================
    !> author: Julien Berland
    !> date: July 2007
    !> 1D interpolation
  !===============================================================
    implicit none
    ! ------------------------------------------------------------
    integer :: nxi,nx
    real(wp), dimension(1:nx) :: x, f  ! interpolant
    real(wp), dimension(1:nxi):: xi,fi ! interpolated function
    character(*) :: method
    ! ------------------------------------------------------------
    integer :: i,j,ir(nxi)
    real(wp) :: rdum,c(4,nx)
    ! ------------------------------------------------------------

    select case (trim(method))

    case ('nearest')

       ! Nearest Neighbor interpolation
       !-------------------------------

       ! Seek indice ir(i) such as xi(i) is the closest to x(ir(i))
       j=1

       do i=1,nxi
          rdum = 1.
          do while ((j.lt.nx).and.(rdum.ge.0.))
             j=j+1
             rdum= xi(i)-x(j)
          enddo

          if ((xi(i)-x(j-1)).lt.(x(j)-xi(i))) then
             ir(i)= j-1
          else
             ir(i)= j
          endif
          j=j-1
       enddo

       ! Interpolation itself
       do i=1,nxi
          fi(i)= f(ir(i))
       enddo

    case ('linear')

       ! Two-point linear interpolation
       !-------------------------------

       ! Seek smallest ir(i) such as xi(i)<x(ir(i))
       j=1

       do i=1,nxi
          rdum =1.
          do while ((j.lt.nx).and.(rdum.ge.0.))
             j=j+1
             rdum= xi(i)-x(j)
          enddo

          ir(i)= j-1
          j=j-1
       enddo

       ! Two-point linear interpolation
       if (ir(nxi).lt.nx) then
          do i=1,nxi
             rdum = (xi(i)-x(ir(i)))/(x(ir(i)+1)-x(ir(i)))
             fi(i)= (1.-rdum)*f(ir(i)) + rdum*f(ir(i)+1)
          enddo
       else
          do i=1,nxi-1
             rdum  = (xi(i)-x(ir(i)))/(x(ir(i)+1)-x(ir(i)))
             fi(i) = (1.-rdum)*f(ir(i)) + rdum*f(ir(i)+1)
          enddo

          i = nxi
          rdum = (xi(i)-x(ir(i)))/(x(ir(i)-1)-x(ir(i)))
          fi(i)= (1.-rdum)*f(ir(i)) + rdum*f(ir(i)-1)
       endif

    case ('cubic_B_spline')

       do i=1,nxi
          call spline_B_val(nx,x,f,xi(i),fi(i))
       enddo

    case ('quadratic_spline')

       do i=1,nxi
          call spline_quadratic_val(nx,x,f,xi(i),fi(i))
       enddo

    case ('Hermite_spline')

       ! set up Hermite interpolant
       do i=1,nx-1
          c(1,i)= f(i)
          c(2,i)= (f(i+1)-f(i))/(x(i+1)-x(i))
       enddo
       c(1,nx)= f(nx)
       c(2,nx)= c(1,nx-1)

       call spline_Hermite_set(nx,x,c)

       ! compute values
       do i=1,nxi
          call spline_Hermite_val(nx,x,c,xi(i),fi(i))
       enddo

    case default
       call mpistop('1D Interpolation: Unknown interpolation method in "interp1"',0)

    end select

  end subroutine interp1

  !=========================================================================
  subroutine spline_B_val(ndata,tdata,ydata,tval,yval)
  !=========================================================================
    !> author: John Burkardt
    !> date: April 1999
    !> Evaluates a cubic B spline approximant.
    !    The cubic B spline will approximate the data, but is not
    !    designed to interpolate it.
    !    In effect, two "phantom" data values are appended to the data,
    !    so that the spline will interpolate the first and last data values.
  !=========================================================================
    implicit none
    ! ----------------------------------------------------------------------
    integer ndata ! number of data values
    real(wp), dimension(ndata) :: tdata ! abscissas of the interpolant data
    real(wp), dimension(ndata) :: ydata ! values of interpolant data
    real(wp) :: tval  ! point at which the spline is to be evaluated
    real(wp) :: yval  ! the value of the function at tval
    ! ----------------------------------------------------------------------
    integer :: i,left,right
    real(wp) :: bval,u
    ! ----------------------------------------------------------------------

    ! Find the nearest interval [TDATA(LEFT),TDATA(RIGHT)] to TVAL.
    do i=2,ndata-1
       if (tval<tdata(i)) then
          left=i-1
          right=i
          goto 10
       endif
    enddo
    left =ndata-1
    right=ndata
10  continue

    ! Evaluate the 5 nonzero B spline basis functions in the interval,
    ! weighted by their corresponding data values.
    u= (tval-tdata(left))/(tdata(right)-tdata(left))
    yval=0.0_wp

    ! B function associated with node LEFT-1, (or "phantom node"),
    ! evaluated in its 4th interval.
    bval=(1.0_wp - 3.0_wp*u + 3.0_wp*u**2 - u**3)/6.0_wp
    if (left-1>0) then
       yval=yval+ydata(left-1)*bval
    else
       yval=yval+(2.0_wp*ydata(1)-ydata(2))*bval
    endif

    ! B function associated with node LEFT,
    ! evaluated in its third interval.
    bval=(4.0_wp - 6.0_wp*u**2 + 3.0_wp*u**3)/6.0_wp
    yval=yval+ydata(left)*bval

    ! B function associated with node RIGHT,
    ! evaluated in its second interval.
    bval=(1.0_wp + 3.0_wp*u + 3.0_wp*u**2 - 3.0_wp*u**3)/6.0_wp
    yval=yval+ydata(right)*bval

    ! B function associated with node RIGHT+1, (or "phantom node"),
    ! evaluated in its first interval.
     bval=u**3/6.0_wp
    if (right+1<=ndata) then
       yval=yval+ydata(right+1)*bval
    else
       yval=yval+(2.0_wp*ydata(ndata)-ydata(ndata-1))*bval
    endif

  end subroutine spline_B_val

  !=========================================================================
  subroutine spline_quadratic_val(ndata,tdata,ydata,tval,yval)
  !=========================================================================
    !> author: John Burkardt
    !> date: April 1999
    !> Evaluates a quadratic spline at a specific point.
    !    Because of the simple form of a piecewise quadratic spline,
    !    the raw data points (TDATA(I),YDATA(I)) can be used directly to
    !    evaluate the spline at any point. No processing of the data
    !    is required.
    !
    ! Rq: The values of TDATA should be distinct and increasing
    !=========================================================================
    implicit none
    ! ----------------------------------------------------------------------
    integer ndata ! number of data values
    real(wp), dimension(ndata) :: tdata ! abscissas of the interpolant data
    real(wp), dimension(ndata) :: ydata ! values of interpolant data
    real(wp) :: tval  ! point at which the spline is to be evaluated
    real(wp) :: yval  ! the value of the function at tval
    ! ----------------------------------------------------------------------
    integer :: i,left,right
    real(wp) :: dif1,dif2
    real(wp) :: t1,t2,t3
    real(wp) :: y1,y2,y3
    real(wp) :: ypval
    ! ----------------------------------------------------------------------

    if (mod(ndata,3)==0) then
       call mpistop('1D Interpolation: Number of data must be odd in "interp1" with the "quadratic_spline"',0)
    endif

    ! Find the interval [TDATA(LEFT),TDATA(RIGHT)] that contains, or is nearest to TVAL
    ! ============================================
    do i=2,ndata-1
       if (tval<tdata(i)) then
          left=i-1
          right=i
          goto 10
       endif
    enddo
    left=ndata-1
    right=ndata
10 continue

    !  Force LEFT to be odd.
    if (mod(left,2)== 0) then
       left=left-1
    endif

    ! Copy out the three abscissas
    ! ============================
    ! force to remain within array bounds
    if (left+2.gt.ndata) left=ndata-2

    t1 = tdata(left)
    t2 = tdata(left+1)
    t3 = tdata(left+2)

    ! Construct and evaluate a parabolic interpolant for the data
    ! in each dimension.
    y1 = ydata(left)
    y2 = ydata(left+1)
    y3 = ydata(left+2)

    dif1=(y2-y1)/(t2-t1)
    dif2=((y3-y1)/(t3-t1)-(y2-y1)/(t2-t1))/(t3-t2)

    yval=y1+(tval-t1)*(dif1+(tval-t2)*dif2)
    ypval=dif1+dif2*(2.0E+00*tval-t1-t2)

  end subroutine spline_quadratic_val

  !=========================================================================
  subroutine spline_Hermite_set(ndata,tdata,c)
  !=========================================================================
    !> author: John Burkardt
    !> date: April 1999
    !> Sets up a piecewise cubic Hermite interpolant.
    !    In effect, two "phantom" data values are appended to the data,
    !    so that the spline will interpolate the first and last data values.
    !
    !    On input, C(1,I) and C(2,I) should contain the value of the
    !    function and its derivative at TDATA(I), for I = 1 to NDATA.
    !    These values will not be changed by this routine.
    !
    !    On output, C(3,I) and C(4,I) contain the quadratic
    !    and cubic coefficients of the Hermite polynomial
    !    in the interval (TDATA(I), TDATA(I+1)), for I=1 to NDATA-1.
    !    C(3,NDATA) and C(4,NDATA) are set to 0.
    !
    !    In the interval (TDATA(I), TDATA(I+1)), the interpolating Hermite
    !    polynomial is given by
    !    SVAL(TVAL) = C(1,I)+(TVAL-TDATA(I))*(C(2,I)+(TVAL-TDATA(I))*(C(3,I)+(TVAL-TDATA(I))*C(4,I)))
    !=========================================================================
    implicit none
    ! ----------------------------------------------------------------------
    integer ndata ! number of data points (must be at least 2)
    real(wp), dimension(ndata) :: tdata ! abscissas of the interpolant data (strictly increasing)
    real(wp), dimension(4,ndata) :: c ! input/output
    ! ----------------------------------------------------------------------
    integer :: i
    real(wp) :: divdif1,divdif3
    real(wp) :: dt
    ! ----------------------------------------------------------------------

    do i=1,ndata-1
       dt=tdata(i+1)-tdata(i)
       divdif1= (c(1,i+1)-c(1,i))/dt
       divdif3= c(2,i)+c(2,i+1)-2.0_wp*divdif1
       c(3,i)= (divdif1-c(2,i)-divdif3)/dt
       c(4,i)= divdif3/dt**2
    enddo

    c(3,ndata)= 0.0_wp
    c(4,ndata)= 0.0_wp

  end subroutine spline_Hermite_set

  !=========================================================================
  subroutine spline_Hermite_val(ndata,tdata,c,tval,sval)
  !=========================================================================
    !> author: John Burkardt
    !> date: April 1999
    !> Evaluates a piecewise cubic Hermite interpolant.
    !    SPLINE_HERMITE_SET must be called first, to set up the
    !    spline data from the raw function and derivative data.
    !=========================================================================
    implicit none
    ! ----------------------------------------------------------------------
    integer ndata ! number of data points (must be at least 2)
    real(wp), dimension(ndata) :: tdata ! abscissas of the interpolant data (strictly increasing)
    real(wp), dimension(4,ndata) :: c ! data computed by SPLINE_HERMITE_SET
    real(wp) :: tval ! the point where the interpolant is to be evaluated
    real(wp) :: sval ! the value of the interpolant at TVAL
    ! ----------------------------------------------------------------------
    integer :: i,left,right
    real(wp) :: dt
    ! ----------------------------------------------------------------------

    ! Find the interval [TDATA(LEFT),TDATA(RIGHT)] that contains or is nearest to TVAL
    ! ============================================
    do i= 2,ndata-1
       if (tval<tdata(i)) then
          left=i-1
          right=i
          goto 10
       endif
    enddo
    left=ndata-1
    right=ndata
10  continue

    ! Evaluate the cubic polynomial
    ! =============================
    dt=tval-tdata(left)

    sval=c(1,left)+dt*(c(2,left)+dt*(c(3,left)+dt*c(4,left)))

  end subroutine spline_Hermite_val

end module mod_interp1
