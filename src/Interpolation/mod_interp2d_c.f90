! ==============================================================================
module mod_interp2d_c
! ==============================================================================
  !> Module for 2D curvilinear interpolation
! ==============================================================================
  use mod_ngh
  implicit none
  ! ----------------------------------------------------------------------------
  integer :: dw ! stencil size (number of interp points per direction)
  integer :: dl ! number of degrees of freedom (dw*dw) of interpolation
  integer :: nint ! total number of points to be interpolated
  ! (i1i,i2i) are indices of point to be interpolated
  integer, dimension(:), allocatable  :: i1i,i2i
  ! (i1s,i2s) are indices of lower-left stencil point for each interp point
  integer, dimension(:), allocatable  :: i1s,i2s
  ! Lagrangian interpolation coefficients
  real, dimension(:,:), allocatable :: coeff_lag
  ! ----------------------------------------------------------------------------

contains

  !=============================================================================
  subroutine interp2d_c(nx1,ny1,x1,y1,u1,x2,y2,u2)
  !=============================================================================
    !> author: XG
    !> date: April 2021
    !> 2D curvilinear interpolation
  !=============================================================================
    implicit none
    ! --------------------------------------------------------------------------
    ! Interpolant field (Donor) - Nota: extended to ghost cells
    integer, intent(in) :: nx1,ny1
    real(wp), dimension(1-ngh:nx1+ngh,1-ngh:ny1+ngh), intent(in) :: x1,y1
    real(wp), dimension(1-ngh:nx1+ngh,1-ngh:ny1+ngh), intent(in) :: u1
    ! Interpolated field (Target) - Nota: without ghost cells
    real(wp), dimension(:,:), intent(in) :: x2,y2
    real(wp), dimension(:,:), intent(inout) :: u2
    ! --------------------------------------------------------------------------
    ! local variables
    integer :: nx2,ny2
    ! --------------------------------------------------------------------------
    
    ! Get dimensions of donor field from input
    ! ========================================
    nx2=size(u2,1)
    ny2=size(u2,2)

    ! Stencil size (fixed in calling procedure)
    ! ============
    ! -> number of degrees of freedom
    dl=dw*dw
    
    ! Preprocessing for interpolations
    ! ================================
    call calc_stencil(nx1,ny1,x1,y1,nx2,ny2,x2,y2)
    
    call calc_coeff_interp2d(nx1,ny1,x1,y1,nx2,ny2,x2,y2)

    ! Perform curvilinear interpolation
    ! =================================
    call apply_interp2dc(nx1,ny1,u1,nx2,ny2,u2)

    ! Free memory
    ! ===========
    deallocate(i1i,i2i,i1s,i2s,coeff_lag)
    
  end subroutine interp2d_c

  !=============================================================================
  subroutine calc_stencil(nx1,ny1,x1,y1,nx2,ny2,x2,y2)
  !=============================================================================
    !> Determination of lower-left point of the interpolation stencil (reference)
    !> OUTPUTS: indices (i1i,i2i) of points to be interpolated
    !>          indices (i1s,i2s) of lower-left stencil points
  !=============================================================================
    use mod_bc
    implicit none
    ! --------------------------------------------------------------------------
    ! Donor grid (known values)
    integer, intent(in) :: nx1,ny1
    real(wp), dimension(1-ngh:nx1+ngh,1-ngh:ny1+ngh), intent(in) :: x1,y1
    ! New grid (where values are to be interpolated)
    integer, intent(in) :: nx2,ny2
    real(wp), dimension(nx2,ny2), intent(in) :: x2,y2
    ! --------------------------------------------------------------------------
    integer :: i,j,ii,jj,l
    integer :: dws2,i1p,i2p
    real(wp) :: min_d,dist
    real(wp) :: x2i,y2i ! position of interpolation point
    logical :: is_inside
    ! --------------------------------------------------------------------------
    ! Mapping function of original grid
    real(wp), parameter :: one=1.0_wp,zer=0.0_wp
    real(wp), dimension(4) :: qx,qy ! quad coordinates
    real(wp), dimension(4) :: a,b ! coeff of bilinear function
    real(wp) :: xm,xl ! logical coordinates
    real(wp), dimension(4,4) :: AI ! inverse coeff matrix
    real(wp) :: aa,bb,cc,det ! work var
    ! --------------------------------------------------------------------------

    ! Mapping functions for original grid
    ! ===================================
    ! used to determine if an interpolation point is inside or outside a cell
    !
    ! The bilinear mapping function [Hughes, T.J.R., The Finite Element Method,
    ! Dover Publications, 20002] is given by:
    !   x=a1+a2l+a3m+a4lm  & y=b1+b2l+b3m+b4lm
    
    !                                [1 0 0 0 ]
    ! Define matrix AI, inverse of A=[1 1 0 0 ]
    !                                [1 1 1 1 ]
    !                                [1 0 1 0 ]
    AI(1,:)=[ one, zer, zer, zer]
    AI(2,:)=[-one, one, zer, zer]
    AI(3,:)=[-one, zer, zer, one]
    AI(4,:)=[ one,-one, one,-one]
    
    ! Total number of points to be interpolated
    ! =========================================
    nint=nx2*ny2

    ! Determination of interpolation stencil
    ! ======================================
    ! allocations
    allocate(i1i(nint),i2i(nint),i1s(nint),i2s(nint))
    ! half stencil size
    dws2=dw/2
    
    l=0
    do ii=1,nx2
       do jj=1,ny2
    !do ii=48,48
    !   do jj=3,3

          is_inside=.false.
          l=l+1 ! counter

          ! point to be interpolated in the new grid (#2)
          ! ----------------------------------------
          i1i(l)=ii ! ~> stored for later use
          i2i(l)=jj ! ~> stored for later use
          x2i=x2(ii,jj)
          y2i=y2(ii,jj)
          
          ! find cell containing the point
          ! ------------------------------
          ! we use the mapping function for each quadrilateral
          do i=1,nx1-1
             do j=1,ny1-1
         
                ! create quadrilateral from original grid
                qx=[x1(i,j),x1(i+1,j),x1(i+1,j+1),x1(i,j+1)]
                qy=[y1(i,j),y1(i+1,j),y1(i+1,j+1),y1(i,j+1)]

                ! compute coefficients
                a=matmul(AI,qx)
                b=matmul(AI,qy)

                !print *,i,j,a,b

                ! calculate logical coordinates		
                aa= a(4)*b(3)-a(3)*b(4)
                bb= a(4)*b(1)-a(1)*b(4)+a(2)*b(3)-a(3)*b(2)+x2i*b(4)-y2i*a(4)
                cc= a(2)*b(1)-a(1)*b(2)+x2i*b(2)-y2i*a(2)

                ! compute m=(-b+sqrt(b^2-4ac))/(2a)
                det=bb*bb-4.0_wp*aa*cc
                if (det>=0) then
                   det=sqrt(det)   
                   xm =(-bb+det)/(2.0_wp*aa)

                   ! compute l
                   xl =(x2i-a(1)-a(3)*xm)/(a(2)+a(4)*xm)

                   ! is the interpolation point inside the quad?
                   if (.not.(xm<0.0_wp.or.xm>1.0_wp.or.xl<0.0_wp.or.xl>1.0_wp)) then
                      i1p=i
                      i2p=j
                      !print *,i1p,i2p
                      is_inside=.true.
                      go to 18
                   endif
                endif
                
             enddo
          enddo

18        continue
          !print *,i1p,i2p,is_inside

          if (.not.is_inside) then
             ! closest point in the donor grid (#1)
             ! -------------------------------
             min_d=1.e9
             do i=1,nx1
                do j=1,ny1
                   dist=sqrt((x1(i,j)-x2i)**2+(y1(i,j)-y2i)**2)
                   if (dist<min_d) then
                      i1p=i
                      i2p=j
                      min_d=dist
                   endif
                enddo
             enddo
          endif
          
          ! location of reference point of the interpolation stencil (lower-left)
          ! --------------------------------------------------------
          if ((is_boundary(1,1)).and.& ! correction for left wall
               ((ii==1).or.(ii==2).or.(ii==3))) then
             i1s(l)=1
          elseif ((is_boundary(1,2)).and.& ! correction for right wall
               ((ii==nx2-2).or.(ii==nx2-1).or.(ii==nx2))) then
             i1s(l)=nx1-dw+1
          else
          if (x1(i1p,i2p)<=x2i) then
             i1s(l)=i1p-dws2+2
             !print *,'a',i1s(l)
          else
             i1s(l)=i1p-dws2+1
             !print *,'b',i1s(l)
          endif
          endif
          if ((is_boundary(1,1)).and.(i1s(l)==0)) i1s(l)=1
          if ((is_boundary(1,2)).and.(i1s(l)==nx1-dw+2)) i1s(l)=nx1-dw+1

          if ((is_boundary(2,1)).and.& ! correction for bottom wall
               ((jj==1).or.(jj==2).or.(jj==3))) then
             i2s(l)=1
          elseif ((is_boundary(2,2)).and.& ! correction for top wall
               ((jj==ny2-2).or.(jj==ny2-1).or.(jj==ny2))) then
             i2s(l)=ny1-dw+1
          else
             if (y1(i1p,i2p)<=y2i) then
                i2s(l)=i2p-dws2+1
             else
                i2s(l)=i2p-dws2
             endif
          endif
          !if ((is_boundary(2,1)).and.(i2s(l)==0)) print *,'in calc_stencil i2s(l)=0',ii,jj
          if ((is_boundary(2,1)).and.(i2s(l)==0)) i2s(l)=1
          if ((is_boundary(2,2)).and.(i2s(l)==ny1-dw+2)) print *,'in calc_stencil i2s(l)=nx-2',ii,jj
          if ((is_boundary(2,2)).and.(i2s(l)==ny1-dw+2)) i2s(l)=ny1-dw+1

          !if (ii==242.and.jj==38) print *,l,i1i(l),i2i(l),i1s(l),i2s(l)
          !print *,l,i1i(l),i2i(l),i1s(l),i2s(l)
          !stop
       enddo
    enddo

  end subroutine calc_stencil

  !=============================================================================
  subroutine calc_coeff_interp2d(nx1,ny1,x1,y1,nx2,ny2,x2,y2)
  !=============================================================================
    !> Compute interpolation coefficients (Lagrangian polynomials)
  !=============================================================================
    implicit none 
    ! --------------------------------------------------------------------------
    ! Donor grid (known values)
    integer, intent(in) :: nx1,ny1
    real(wp), dimension(1-ngh:nx1+ngh,1-ngh:ny1+ngh), intent(in) :: x1,y1
    ! New grid (where values are to be interpolated)
    integer, intent(in) :: nx2,ny2
    real(wp), dimension(nx2,ny2), intent(in) :: x2,y2
    ! --------------------------------------------------------------------------
    !integer :: ngh
    ! --------------------------------------------------------------------------
    integer :: i,j,l,ptint
    real(wp) :: coef_sum
    ! interpolation point
    integer :: io,jo ! lower-left corner of stencil
    real(wp) :: xint,yint,x0,y0 ! interpolation point
    ! position of stencil points in physical space
    real(wp), dimension(:), allocatable :: xsten,ysten
    ! position of stencil points in unitary stencil
    real(wp), dimension(:), allocatable :: xsten_c,ysten_c
    real(wp), dimension(:), allocatable :: coeff_lag_x,coeff_lag_y
    ! computational grid
    real(wp) :: pasx,pasy
    real(wp), dimension(:), allocatable :: x_c,y_c
    ! position of interpolation point in computational stencil measured from lower-left corner
    real(wp) :: deltax_c,deltay_c
    ! --------------------------------------------------------------------------

    ! Computational grid (unitary regular mesh)
    ! ==================
    allocate(x_c(1-ngh:nx1+ngh))
    pasx=1.0_wp/real(nx1-1)
    x_c(1)=0.0_wp
    do i=2,nx1+ngh
       x_c(i)=x_c(i-1)+pasx
    enddo
    ! potentially interface ??????? TO BE CHANGED
    do i=0,1-ngh,-1
       x_c(i)=x_c(i+1)-pasx
    enddo

    !! allocate(y_c(-1:ny1+5)) TO BE CHANGED
    allocate(y_c(1-ngh:ny1+ngh))
    pasy=1.0_wp/real(ny1-1)
    y_c(1)=0.0_wp
    do i=2,ny1+ngh
       y_c(i)=y_c(i-1)+pasy
    enddo
    do i=0,1-ngh,-1
       y_c(i)=y_c(i+1)-pasy
    enddo

    ! Allocations
    ! ===========
    ! stencil in physical space (curvilinear)
    allocate(xsten(dl),ysten(dl))
    ! stencil in computational space (Cartesian)
    allocate(xsten_c(dw),ysten_c(dw))
    ! interpolation coefficients for all interpolation points
    allocate(coeff_lag_x(dw),coeff_lag_y(dw))
    allocate(coeff_lag(dl,nint))

    ! Start loop over the nint interpolation points
    ! =============================================
    do ptint=1,nint

       ! position (x0,y0) of interpolation point in donor grid
       ! -----------------------------------------------------
       x0=x2(i1i(ptint),i2i(ptint))
       y0=y2(i1i(ptint),i2i(ptint))

       ! indices of lower-left point of stencil in donor grid
       ! ----------------------------------------------------
       io=i1s(ptint)
       jo=i2s(ptint)
      
       ! compute stencils
       ! ----------------
       ! coordinates of stencil points in physical space
       xsten=0.0_wp
       ysten=0.0_wp
       l=0 ! counter
       do i=1,dw
          do j=1,dw
             l=l+1 ! increment counter
             xsten(l)=x1(io+(i-1),jo+(j-1))
             ysten(l)=y1(io+(i-1),jo+(j-1))
          end do
       end do

       ! coordinates of stencil points in computational space
       do i=1,dw
          xsten_c(i)=x_c(io+(i-1))
       end do
       do j=1,dw
          ysten_c(j)=y_c(jo+(j-1))
       end do

       ! compute position of interpolation point in computational stencil
       ! measured from lower-left corner (deltax_c,deltay_c)
       ! ---------------------------------------------------
       ! initial guess (middle of stencil)
       deltax_c=dw/2.0_wp-0.5_wp
       deltay_c=dw/2.0_wp-0.5_wp

       ! check position of interpolation point in computational stencil
       !if (ptint.eq.5000) &
       !     print *,'       initial: deltax_c= ',deltax_c,'deltay_c= ',deltay_c

       ! apply Newton-Raphson method to compute (deltax_c,deltay_c)
       call newt_2D(deltax_c,deltay_c,xsten,ysten,x0,y0)

       ! check
       !if (ptint.eq.5000) &
       !     print *,'  after Newton: deltax_c= ',deltax_c,'deltay_c= ',deltay_c

       ! position in computational stencil
       xint=x_c(io)+deltax_c*pasx
       yint=y_c(jo)+deltay_c*pasy

       ! compute Lagrangian interpolation coefficients
       ! ---------------------------------------------
       ! x-direction
       coeff_lag_x=1.0_wp
       do j=1,dw
          do l=1,dw
             if (l.ne.j) then
                coeff_lag_x(j)=coeff_lag_x(j)*(xint-xsten_c(l))/(xsten_c(j)-xsten_c(l))
             end if
          end do
       end do
       ! y-direction
       coeff_lag_y=1.0_wp
       do j=1,dw
          do l=1,dw
             if (l.ne.j) then
                coeff_lag_y(j)=coeff_lag_y(j)*(yint-ysten_c(l))/(ysten_c(j)-ysten_c(l))
             end if
          end do
       end do
       ! gather coefficients
       do i=1,dw
          do j=1,dw
             coeff_lag((i-1)*dw+j,ptint)=coeff_lag_x(i)*coeff_lag_y(j)
          end do
       end do
       ! check sum
       coef_sum=0.
       do i=1,dl
          coef_sum=coef_sum+coeff_lag(i,ptint)
       enddo

    end do
    ! End loop over the nint interpolation points
    ! ===========================================

    deallocate(xsten,ysten,xsten_c,ysten_c)
    deallocate(coeff_lag_x,coeff_lag_y)
    
  end subroutine calc_coeff_interp2d

  !=============================================================================
  subroutine newt_2D(deltax_c,deltay_c,xsten,ysten,x0,y0)
  !=============================================================================
    !> Newton-Raphson method to compute (deltax_c,deltay_c)
  !=============================================================================
    implicit none 
    ! --------------------------------------------------------------------------
    real(wp), intent(inout) :: deltax_c,deltay_c
    real(wp), intent(in) :: x0,y0
    real(wp), dimension(1:dl), intent(in) :: xsten,ysten 
    ! --------------------------------------------------------------------------
    integer :: it
    real(wp) :: fx,fy,resid
    ! Jacobian matrix
    real(wp) :: J11,J22,J12,J21,det_J
    ! --------------------------------------------------------------------------
    
    ! Initialization of residual
    resid=999.0_wp

    ! 2D Newton's algorithm
    ! =====================
    it=0
    
    do while ((resid>1.0e-2_wp).and.(it<3000))

       it=it+1 ! counter iterations

       ! Jacobian matrix
       ! ---------------
       J11=d1_ifunc(deltax_c,deltay_c,xsten)
       J22=d2_ifunc(deltax_c,deltay_c,ysten)
       J12=d2_ifunc(deltax_c,deltay_c,xsten)
       J21=d1_ifunc(deltax_c,deltay_c,ysten)
       ! its determinant
       det_J=J11*J22-J12*J21

       ! increments
       ! ----------
       fx=ifunc(deltax_c,deltay_c,xsten,x0)
       fy=ifunc(deltax_c,deltay_c,ysten,y0)

       ! update deltax_c & deltay_c
       ! --------------------------
       deltax_c=deltax_c - (fx*J22-fy*J12)/det_J
       deltay_c=deltay_c - (fy*J11-fx*J21)/det_J

       ! residual of increments (L2-norm)
       ! ----------------------
       resid=sqrt(fx**2+fy**2)

    enddo
    
    !print *,'number of iterations',it

  end subroutine newt_2D

  !=============================================================================
  function ifunc(deltax,deltay,xsten,x0)
  !=============================================================================
    !> Compute function F = sum(x')-x0 = 0
  !=============================================================================
    implicit none 
    ! --------------------------------------------------------------------------
    real(wp), intent(in) :: deltax,deltay
    real(wp), intent(in) :: x0
    real(wp), dimension(1:dl), intent(in) :: xsten
    real(wp) :: ifunc ! function output
    ! --------------------------------------------------------------------------
    integer :: i,j
    real(wp) :: productx,producty,alfa,sum_x
    real(wp), dimension(:), allocatable :: coeff_lag_x,coeff_lag_y,coeff_lagf
    ! --------------------------------------------------------------------------
    
    ! Allocations of Lagrangian interpolation coefficients
    ! ====================================================
    allocate(coeff_lag_x(1:dw),coeff_lag_y(1:dw))
    allocate(coeff_lagf(1:dl))

    
    ! Compute Lagrangian interpolation coefficients (unitary grid)
    ! =============================================
    do j=0,dw-1
       productx=1.0_wp
       producty=1.0_wp
       do i=0,dw-1
          if (i.ne.j) then
             productx=productx*(deltax-real(i))
             producty=producty*(deltay-real(i))
          end if
       end do

       alfa=((-1.0_wp)**(dw+j-1))/(fact(dw-(j+1))*fact(j))

       coeff_lag_x(j+1)=alfa*productx
       coeff_lag_y(j+1)=alfa*producty
    end do
    ! gather coefficients
    do i=1,dw
       do j=1,dw
          coeff_lagf((i-1)*dw+j)=coeff_lag_x(i)*coeff_lag_y(j)
       end do
    end do

    ! Compute sum(x')-x0
    ! ==================
    sum_x=0.
    do i=1,dl
       sum_x=sum_x+coeff_lagf(i)*xsten(i)
    end do

    ifunc=sum_x-x0

    ! free memory
    deallocate(coeff_lag_x,coeff_lag_y,coeff_lagf)

  end function ifunc

  !=============================================================================
  function d1_ifunc(deltax,deltay,xsten)
  !=============================================================================
    !> Compute derivative along 1 of function F = sum(x')-x0 = 0
  !=============================================================================
    implicit none 
    ! --------------------------------------------------------------------------
    real(wp), intent(in) :: deltax,deltay
    real(wp), dimension(1:dl), intent(in) :: xsten
    real(wp) :: d1_ifunc ! function output
    ! --------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: productx,dproductx,producty,alfa,sum_x
    real(wp), dimension(:), allocatable :: coeff_lag_x,coeff_lag_y,coeff_lagf
    ! --------------------------------------------------------------------------

    ! Allocations of Lagrangian interpolation coefficients
    ! ====================================================
    allocate(coeff_lag_x(1:dw),coeff_lag_y(1:dw))
    allocate(coeff_lagf(1:dl))

    ! Compute Lagrangian interpolation coefficients (unitary grid)
    ! =============================================
    do j=0,dw-1
       dproductx=0.0_wp
       producty=1.0_wp
       do i=0,dw-1
          if (i.ne.j) then
             productx=1.0_wp
             do k=0,dw-1
                if ((k.ne.i).and.(k.ne.j)) then
                   productx=productx*(deltax-real(k))
                end if
             end do
             dproductx=dproductx+productx
          end if

          if(i.ne.j) then
             producty=producty*(deltay-real(i))
          end if

       end do

       alfa=((-1.0_wp)**(dw+j-1))/(fact(dw-(j+1))*fact(j))
 
       coeff_lag_x(j+1)=alfa*dproductx
       coeff_lag_y(j+1)=alfa*producty
    end do
    ! gather coefficients
    do i=1,dw
       do j=1,dw
          coeff_lagf((i-1)*dw+j)=coeff_lag_x(i)*coeff_lag_y(j)
       end do
    end do

    ! Compute derivative along 1 of function F = sum(x')-x0 = 0
    ! =========================================================
    sum_x=0.
    do i=1,dl
       sum_x=sum_x+coeff_lagf(i)*xsten(i)
    end do

    d1_ifunc=sum_x

    ! free memory
    deallocate(coeff_lag_x,coeff_lag_y,coeff_lagf)
    
  end function d1_ifunc

  !=============================================================================
  function d2_ifunc(deltax,deltay,xsten)
  !=============================================================================
    !> Compute derivative along 2 of function F = sum(x')-x0 = 0
  !=============================================================================
    implicit none 
    ! --------------------------------------------------------------------------
    real(wp), intent(in) :: deltax,deltay
    real(wp), dimension(1:dl), intent(in) :: xsten
    real(wp) :: d2_ifunc ! function output
    ! --------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: productx,producty,dproducty,alfa,sum_x
    real(wp), dimension(:), allocatable :: coeff_lag_x,coeff_lag_y,coeff_lagf
    ! --------------------------------------------------------------------------

    ! Allocations of Lagrangian interpolation coefficients
    ! ====================================================
    allocate(coeff_lag_x(1:dw),coeff_lag_y(1:dw))
    allocate(coeff_lagf(1:dl))

    ! Compute Lagrangian interpolation coefficients (unitary grid)
    ! =============================================
    do j=0,dw-1
       productx=1.0_wp
       dproducty=0.0_wp
       do i=0,dw-1
          if (i.ne.j) then
             productx=productx*(deltax-real(i))
          end if

          if(i.ne.j) then
             producty=1.0_wp
             do k=0,dw-1
                if((k.ne.i).and.(k.ne.j)) then
                   producty=producty*(deltay-real(k))
                end if
             end do
             dproducty=dproducty+producty
          end if
       end do

       alfa=((-1.0_wp)**(dw+j-1))/(fact(dw-(j+1))*fact(j))

       coeff_lag_x(j+1)=alfa*productx
       coeff_lag_y(j+1)=alfa*dproducty
    end do
    ! gather coefficients
    do i=1,dw
       do j=1,dw
          coeff_lagf((i-1)*dw+j)=coeff_lag_x(i)*coeff_lag_y(j)
       end do
    end do

    ! Compute derivative along 2 of function F = sum(x')-x0 = 0
    ! =========================================================
    sum_x=0.
    do i=1,dl
       sum_x=sum_x+coeff_lagf(i)*xsten(i)
    end do

    d2_ifunc=sum_x

    ! free memory
    deallocate(coeff_lag_x,coeff_lag_y,coeff_lagf)

  end function d2_ifunc

  !=============================================================================
  function fact(int_num)
  !=============================================================================
    !> Factorial function
  !=============================================================================
    implicit none
    ! --------------------------------------------------------------------------
    integer, intent(in) :: int_num
    real(wp) :: fact ! function output
    ! --------------------------------------------------------------------------
    integer :: i
    ! --------------------------------------------------------------------------

    fact=1.0_wp
    
    do i=1,int_num
       fact=fact*real(i)
    end do

  end function fact

  !=============================================================================
  subroutine apply_interp2dc(nx1,ny1,u1,nx2,ny2,u2)
  !=============================================================================
    !> apply 2D curvilinear interpolations (based on Lagrangian polynomials)
  !=============================================================================
    implicit none
    ! --------------------------------------------------------------------------
    integer, intent(in) :: nx1,ny1
    integer, intent(in) :: nx2,ny2
    real(wp), dimension(1-ngh:nx1+ngh,1-ngh:ny1+ngh) :: u1
    real(wp), dimension(nx2,ny2) :: u2
    ! --------------------------------------------------------------------------
    integer :: i,j,l,ptint,ii,ji,iv,jv
    ! --------------------------------------------------------------------------

    do ptint=1,nint

       ii=i1i(ptint)
       ji=i2i(ptint)
       iv=i1s(ptint)
       jv=i2s(ptint)

       u2(ii,ji)=0.0_wp
       do i=1,dw
          do j=1,dw
             l=(i-1)*dw+j
             u2(ii,ji)= u2(ii,ji) + coeff_lag(l,ptint)*u1((i-1)+iv,(j-1)+jv)
          end do
       end do

    enddo

  end subroutine apply_interp2dc

!!$  !=============================================================================
!!$  subroutine interp_z(u1,u1i)
!!$  !=============================================================================
!!$    !> Mid-point interpolations along z
!!$    !> if nz2=2*nz1 (doubling regular grid along z)
!!$  !=============================================================================
!!$    implicit none
!!$    ! --------------------------------------------------------------------------
!!$    integer, parameter :: ngh=5
!!$    real(wp), dimension(1-ngh:nx1+ngh,ny1,1-ngh:nz1+ngh) :: u1
!!$    real(wp), dimension(1-ngh:nx1+ngh,ny1,nz2) :: u1i
!!$    ! --------------------------------------------------------------------------
!!$    integer :: i,j,k
!!$    ! --------------------------------------------------------------------------
!!$
!!$    do i=1-ngh,nx1+ngh
!!$       do j=1,ny1
!!$
!!$          do k=1,nz2
!!$             if (mod(k+1,2)==0) then 
!!$                ! coincident planes
!!$                !print *,k,(k+1)/2
!!$                u1i(i,j,k)=u1(i,j,(k+1)/2)
!!$             else
!!$                !print *,k,k/2,k/2+1
!!$                u1i(i,j,k)=0.5_wp*(u1(i,j,k/2)+u1(i,j,k/2+1))
!!$             endif
!!$          enddo
!!$
!!$       enddo
!!$    enddo
!!$
!!$  end subroutine interp_z

!!$  !=============================================================================
!!$  subroutine interp_field(nomvar)
!!$  !=============================================================================
!!$    !> Read / Interpolate / Write field
!!$  !=============================================================================
!!$    implicit none
!!$    ! --------------------------------------------------------------------------
!!$    character(len=4) :: nomvar
!!$    ! --------------------------------------------------------------------------
!!$    integer :: i,j,k,l,i1,j1,k1
!!$    real(wp), dimension(:,:,:), allocatable :: u1,u1i,u2
!!$    ! --------------------------------------------------------------------------
!!$
!!$    allocate( u1(1-ngh:nx1+ngh,ny1,1-ngh:nz1+ngh)) ! <- donor field
!!$    allocate(u1i(1-ngh:nx1+ngh,ny1,nz2)) ! <- intermediate field (z-interpolation)
!!$
!!$    call read_var(u1(1:nx1,1:ny1,1:nz1),nomvar)
!!$    ! enforce periodicity along x
!!$    u1(nx1+1:nx1+ngh,:,:)=u1(1:ngh,:,:)
!!$    u1(1-ngh:0,:,:)=u1(nx1-ngh+1:nx1,:,:)
!!$    ! enforce periodicity along z
!!$    u1(:,:,nz1+1:nz1+ngh)=u1(:,:,1:ngh)
!!$    u1(:,:,1-ngh:0)=u1(:,:,nz1-ngh+1:nz1)
!!$
!!$    call interp_z(u1,u1i)
!!$
!!$    deallocate(u1)
!!$    allocate(u2(nx2,ny2,nz2))           ! <- a interpoler
!!$
!!$    do k=1,nz2
!!$       call interp_xy(u1i(:,:,k),u2(:,:,k))
!!$    enddo
!!$
!!$    ! correct for no-slip-wall condition
!!$    if ((nomvar=='rhou').or.(nomvar=='rhov').or.(nomvar=='rhow') &
!!$         .or.(nomvar=='u1  ').or.(nomvar=='v1  ').or.(nomvar=='w1  ')) then
!!$       u2(:,1,:)=0.0_wp
!!$       u2(:,ny2,:)=0.0_wp
!!$    endif
!!$
!!$    call ecrit_var(u2,nomvar)
!!$
!!$    if (nomvar=='rhou') then
!!$
!!$       i1=nx2/2
!!$       j1=ny2/2
!!$       !j1=1
!!$       k1=nz2/2-1
!!$       open(51,file='rho.bin',form='unformatted',status='unknown')
!!$       rewind(51)
!!$       write(51) nx2
!!$       write(51) ny2
!!$       write(51) nz2
!!$       write(51) ((x2(i,j),i=1,nx2),j=1,ny2)
!!$       write(51) ((y2(i,j),i=1,nx2),j=1,ny2)
!!$       write(51) (z2(k),k=1,nz2)
!!$       write(51) ((u2(i,j,k1),i=1,nx2),j=1,ny2)
!!$       write(51) ((u2(i,j1,k),i=1,nx2),k=1,nz2)
!!$       write(51) ((u2(i1,j,k),j=1,ny2),k=1,nz2)
!!$
!!$       write(51) nint
!!$       write(51) dw
!!$       write(51) (i1i(l),l=1,nint)
!!$       write(51) (i2i(l),l=1,nint)
!!$       write(51) (i1s(l),l=1,nint)
!!$       write(51) (i2s(l),l=1,nint)
!!$       close(51)
!!$    endif
!!$
!!$    deallocate(u1i,u2)
!!$
!!$  end subroutine interp_field
  
end module mod_interp2d_c
