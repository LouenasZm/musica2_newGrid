!===============================================================================
subroutine e1xa(x,e1)
!===============================================================================
  !> author: Shanjie Zhang, Jianming Jin
  !> date: 06 July 2012
  !> Compute exponential integral E1(x)
  !  Licensing:
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin. However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !  Reference: Shanjie Zhang, Jianming Jin,
  !             Computation of Special Functions, Wiley, 1996,
  !             ISBN: 0-471-11963-6, LC: QA351.C45.
!===============================================================================
  use precision
  implicit none
  ! ------------------------------------------------------------
  real(wp) :: x  ! Input , real(kind=8) X, the argument.
  real(wp) :: e1 ! Output, real(kind=8) E1, the function value
  ! ------------------------------------------------------------
  real(wp) :: es1,es2
  ! ------------------------------------------------------------

  if (x==0.0_wp) then

     e1 =1.0e+300_wp

  else if (x<=1.0_wp) then

     e1 = - log(x) + ((((      &
            1.07857e-03_wp*x   &
          - 9.76004e-03_wp)*x  &
          + 5.519968e-02_wp)*x &
          - 0.24991055_wp)*x   &
          + 0.99999193_wp)*x   &
          - 0.57721566_wp

  else

     es1 = ((( x &
          + 8.5733287401_wp)*x &
          +18.059016973_wp )*x &
          + 8.6347608925_wp)*x &
          + 0.2677737343_wp

     es2 = ((( x &
          +  9.5733223454_wp)*x &
          + 25.6329561486_wp)*x &
          + 21.0996530827_wp)*x &
          +  3.9584969228_wp

     e1 =exp(-x)/x*es1/es2

  end if

end subroutine e1xa

!===============================================================================
subroutine e1xb(x,e1)
!===============================================================================
  !> author: Shanjie Zhang, Jianming Jin
  !> date: 06 July 2012
  !> Compute exponential integral E1(x)
  !  Licensing:
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin. However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !  Reference: Shanjie Zhang, Jianming Jin,
  !             Computation of Special Functions, Wiley, 1996,
  !             ISBN: 0-471-11963-6, LC: QA351.C45.
!===============================================================================
  use precision
  implicit none
  ! ------------------------------------------------------------
  real(wp) :: x  ! Input , real(kind=8) X, the argument.
  real(wp) :: e1 ! Output, real(kind=8) E1, the function value
  ! ------------------------------------------------------------
  integer :: k,m
  real(wp) :: ga,r,t,t0
  ! ------------------------------------------------------------

  if (x==0.0_wp) then
     
     e1 =1.0e+300_wp
     
  else if (x<=1.0_wp) then
     
     e1=1.0_wp
     r =1.0_wp

     do k=1,25
        r =-r*k*x/(k+1.0_wp)**2
        e1= e1+r
        if (abs(r)<=abs(e1)*1.0e-15_wp) then
           exit
        end if
     end do

     ga =0.5772156649015328_wp
     e1 =-ga-log(x)+x*e1

  else

     m=20+int(80.0_wp/x)
     t0=0.0_wp
     do k=m,1,-1
        t0=k/(1.0_wp+k/(x+t0))
     end do
     t =1.0_wp/(x+t0)
     e1=exp(-x)*t

  end if

end subroutine e1xb
