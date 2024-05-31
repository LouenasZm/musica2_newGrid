!===============================================================================
subroutine grid_metrics_cart(dx_,x_,nx_,dir,ngh_)
!===============================================================================
  !> Compute Cartesian metrics
!===============================================================================
  use mod_mpi_part
  use mod_grid
  use mod_bc
  use mod_coeff_deriv
  use warnstop
  implicit none
  ! ---------------------------------------------------------------------------
  integer, intent(in) :: dir  ! direction
  integer, intent(in) :: nx_  ! size
  integer, intent(in) :: ngh_ ! number of ghost points
  real(wp), dimension(1-ngh:nx_+ngh), intent(in)  :: x_  ! grid
  real(wp), dimension(1:nx_), intent(out) :: dx_ ! metrics
  ! ---------------------------------------------------------------------------
  ! local variables
  integer :: i,l
  integer :: nx_min,nx_max
  real(wp), dimension(-ngh_:ngh_) :: a_ ! DF scheme
  ! ---------------------------------------------------------------------------

  ! 1-D metrics along direction dir
  ! ===============================
  
  ! Initializations
  ! ---------------
  !allocate(dx_(-4:nx_+5))
  dx_=0.0_wp

  ! BC i_min
  ! --------
  nx_min=1
  if (BC_face(dir,1)%sort<=0) then
     nx_min=1+ngh_
     select case (ngh_)
     case (5)
        i=1
        do l=0,10
           dx_(i)=dx_(i)+a010(1+l)*x_(i+l)
        enddo
        i=2
        do l=-1,9
           dx_(i)=dx_(i)+a19(2+l)*x_(i+l)
        enddo
        i=3
        do l=-2,8
           dx_(i)=dx_(i)+a28(3+l)*x_(i+l)
        enddo
        i=4
        do l=-3,7
           dx_(i)=dx_(i)+a37(4+l)*x_(i+l)
        enddo
        i=5
        do l=-4,6
           dx_(i)=dx_(i)+a46(5+l)*x_(i+l)
        enddo
     case (4)
        i=1
        do l=0,8
           dx_(i)=dx_(i)+a08(1+l)*x_(i+l)
        enddo
        i=2
        do l=-1,7
           dx_(i)=dx_(i)+a17(2+l)*x_(i+l)
        enddo
        i=3
        do l=-2,6
           dx_(i)=dx_(i)+a26(3+l)*x_(i+l)
        enddo
        i=4
        do l=-3,5
           dx_(i)=dx_(i)+a35(4+l)*x_(i+l)
        enddo
     case (3)
        i=1
        do l=0,6
           dx_(i)=dx_(i)+a06(1+l)*x_(i+l)
        enddo
        i=2
        do l=-1,5
           dx_(i)=dx_(i)+a15(2+l)*x_(i+l)
        enddo
        i=3
        do l=-2,4
           dx_(i)=dx_(i)+a24(3+l)*x_(i+l)
        enddo
     case (2)
        i=1
        do l=0,4
           dx_(i)=dx_(i)+a04(1+l)*x_(i+l)
        enddo
        i=2
        do l=-1,3
           dx_(i)=dx_(i)+a13(2+l)*x_(i+l)
        enddo
     case (1)
        i=1
        do l=0,2
           dx_(i)=dx_(i)+a02(1+l)*x_(i+l)
        enddo
        !dx_(1)=x_(2)-x_(1)
     case default
        call mpistop('max ngh_ is 5 (11pts-stencil)!', 0)
     end select
  end if

  ! BC i_max
  ! --------
  nx_max=nx_
  if (BC_face(dir,2)%sort<=0) then
     nx_max=nx_-ngh_
     select case (ngh_)
     case (5)
        i=nx_-4
        do l=-6,4
           dx_(i)=dx_(i)-a46(5-l)*x_(i+l)
        enddo
        i=nx_-3
        do l=-7,3
           dx_(i)=dx_(i)-a37(4-l)*x_(i+l)
        enddo
        i=nx_-2
        do l=-8,2
           dx_(i)=dx_(i)-a28(3-l)*x_(i+l)
        enddo
        i=nx_-1
        do l=-9,1
           dx_(i)=dx_(i)-a19(2-l)*x_(i+l)
        enddo
        i=nx_
        do l=-10,0
           dx_(i)=dx_(i)-a010(1-l)*x_(i+l)
        enddo
     case (4)
        i=nx_-3
        do l=-5,3
           dx_(i)=dx_(i)-a35(4-l)*x_(i+l)
        enddo
        i=nx_-2
        do l=-6,2
           dx_(i)=dx_(i)-a26(3-l)*x_(i+l)
        enddo
        i=nx_-1
        do l=-7,1
           dx_(i)=dx_(i)-a17(2-l)*x_(i+l)
        enddo
        i=nx_
        do l=-8,0
           dx_(i)=dx_(i)-a08(1-l)*x_(i+l)
        enddo
     case (3)
        i=nx_-2
        do l=-4,2
           dx_(i)=dx_(i)-a24(3-l)*x_(i+l)
        enddo
        i=nx_-1
        do l=-5,1
           dx_(i)=dx_(i)-a15(2-l)*x_(i+l)
        enddo
        i=nx_
        do l=-6,0
           dx_(i)=dx_(i)-a06(1-l)*x_(i+l)
        enddo
     case (2)
        i=nx_-1
        do l=-3,1
           dx_(i)=dx_(i)-a13(2-l)*x_(i+l)
        enddo
        i=nx_
        do l=-4,0
           dx_(i)=dx_(i)-a04(1-l)*x_(i+l)
        enddo
     case (1)
        i=nx_
        do l=-2,0
           dx_(i)=dx_(i)-a02(1-l)*x_(i+l)
        enddo
        !dx_(nx_)=x_(nx_)-x_(nx_-1)
     case default
        call mpistop('max ngh_ is 5 (11pts-stencil)!', 0)
     end select
  end if

  ! Interior points
  ! ---------------
  select case (ngh_)
  case (5)
     a_=a11
  case (4)
     a_=a9
  case (3)
     a_=a7
  case (2)
     a_=a5
  case (1)
     a_=a3        
  case default
     call mpistop('max ngh_ is 5 (11pts-stencil)!', 0)
  end select
  
  do i=nx_min,nx_max
     do l=-ngh_,ngh_
        dx_(i) = dx_(i) + a_(l)*x_(i+l)
     enddo
  enddo
  
  ! Inverse of metrics (*idx instead of /dx)
  ! ==================
  do i=1,nx_
     dx_(i)=1.0_wp/dx_(i)
  enddo
  
end subroutine grid_metrics_cart

!===============================================================================
subroutine grid_metrics_wall_cart
!===============================================================================
  !> Compute near-wall Cartesian metrics
!===============================================================================
  use mod_grid
  use mod_bc
  use mod_coeff_deriv
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j,k
  ! ---------------------------------------------------------------------------

  ! 3-point stencil
  ! ===============
  if ((BC_face(1,1)%sort<=0).and.(.not.is_curv)) then
     i=1
     !idx1_imin= x(i+1)-x(i)
     idx1_imin=a02(1)*x(i)+a02(2)*x(i+1)+a02(3)*x(i+2) 
     idx1_imin=1.0_wp/idx1_imin
  endif
  if ((BC_face(1,2)%sort<=0).and.(.not.is_curv)) then
     i=nx
     !idx1_imax= x(i)-x(i-1)     
     idx1_imax=a20(1)*x(i)+a20(2)*x(i-1)+a20(3)*x(i-2)
     idx1_imax=1.0_wp/idx1_imax
  endif
  if ((BC_face(2,1)%sort<=0).and.(.not.is_curv)) then
     j=1
     !idy1_jmin= y(j+1)-y(j)    
     idy1_jmin=a02(1)*y(j)+a02(2)*y(j+1)+a02(3)*y(j+2) 
     idy1_jmin=1.0_wp/idy1_jmin
  endif
  if ((BC_face(2,2)%sort<=0).and.(.not.is_curv)) then
     j=ny
     !idy1_jmax= y(j)-y(j-1)   
     idy1_jmax=a20(1)*y(j)+a20(2)*y(j-1)+a20(3)*y(j-2)
     idy1_jmax=1.0_wp/idy1_jmax
  endif
  if (BC_face(3,1)%sort<=0) then
     k=1
     !idz1_kmin= z(k+1)-z(k)    
     idz1_kmin=a02(1)*z(k)+a02(2)*z(k+1)+a02(3)*z(k+2) 
     idz1_kmin=1.0_wp/idz1_kmin
  endif
  if (BC_face(3,2)%sort<=0) then
     k=nz
     !idz1_kmax= z(k)-z(k-1)
     idz1_kmax=a20(1)*z(k)+a20(2)*z(k-1)+a20(3)*z(k-2)
     idz1_kmax=1.0_wp/idz1_kmax
  endif

  if (ngh==1) return
  
  ! 3-point stencil
  ! ===============
  if ((BC_face(1,1)%sort<=0).and.(.not.is_curv)) then
     i=2
     idx2_imin= a3(1)*(x(i+1)-x(i-1))
     idx2_imin=1.0_wp/idx2_imin
  endif
  if ((BC_face(1,2)%sort<=0).and.(.not.is_curv)) then
     i=nx-1
     idx2_imax= a3(1)*(x(i+1)-x(i-1))     
     idx2_imax=1.0_wp/idx2_imax
  endif
  if ((BC_face(2,1)%sort<=0).and.(.not.is_curv)) then
     j=2
     idy2_jmin= a3(1)*(y(j+1)-y(j-1))    
     idy2_jmin=1.0_wp/idy2_jmin
  endif
  if ((BC_face(2,2)%sort<=0).and.(.not.is_curv)) then
     j=ny-1
     idy2_jmax= a3(1)*(y(j+1)-y(j-1))    
     idy2_jmax=1.0_wp/idy2_jmax
  endif
  if (BC_face(3,1)%sort<=0) then
     k=2
     idz2_kmin= a3(1)*(z(k+1)-z(k-1))     
     idz2_kmin=1.0_wp/idz2_kmin
  endif
  if (BC_face(3,2)%sort<=0) then
     k=nz-1
     idz2_kmax= a3(1)*(z(k+1)-z(k-1))   
     idz2_kmax=1.0_wp/idz2_kmax
  endif
   
  if (ngh==2) return
     
  ! 5-point stencil
  ! ===============
  if ((BC_face(1,1)%sort<=0).and.(.not.is_curv)) then
     i=3
     idx4_imin= a5(1)*(x(i+1)-x(i-1)) &
              + a5(2)*(x(i+2)-x(i-2))     
     idx4_imin=1.0_wp/idx4_imin
  endif
  if ((BC_face(1,2)%sort<=0).and.(.not.is_curv)) then
     i=nx-2
     idx4_imax= a5(1)*(x(i+1)-x(i-1)) &
              + a5(2)*(x(i+2)-x(i-2))     
     idx4_imax=1.0_wp/idx4_imax
  endif
  if ((BC_face(2,1)%sort<=0).and.(.not.is_curv)) then
     j=3
     idy4_jmin= a5(1)*(y(j+1)-y(j-1)) &
              + a5(2)*(y(j+2)-y(j-2))     
     idy4_jmin=1.0_wp/idy4_jmin
  endif
  if ((BC_face(2,2)%sort<=0).and.(.not.is_curv)) then
     j=ny-2
     idy4_jmax= a5(1)*(y(j+1)-y(j-1)) &
              + a5(2)*(y(j+2)-y(j-2))    
     idy4_jmax=1.0_wp/idy4_jmax
  endif
  if (BC_face(3,1)%sort<=0) then
     k=3
     idz4_kmin= a5(1)*(z(k+1)-z(k-1)) &
              + a5(2)*(z(k+2)-z(k-2))     
     idz4_kmin=1.0_wp/idz4_kmin
  endif
  if (BC_face(3,2)%sort<=0) then
     k=nz-2
     idz4_kmax= a5(1)*(z(k+1)-z(k-1)) &
              + a5(2)*(z(k+2)-z(k-2))    
     idz4_kmax=1.0_wp/idz4_kmax
  endif
   
  if (ngh==3) return
     
  ! 7-point stencil
  ! ===============
  if ((BC_face(1,1)%sort<=0).and.(.not.is_curv)) then
     i=4
     idx6_imin= a7(1)*(x(i+1)-x(i-1)) &
              + a7(2)*(x(i+2)-x(i-2)) &
              + a7(3)*(x(i+3)-x(i-3))     
     idx6_imin=1.0_wp/idx6_imin
  endif
  if ((BC_face(1,2)%sort<=0).and.(.not.is_curv)) then
     i=nx-3
     idx6_imax= a7(1)*(x(i+1)-x(i-1)) &
              + a7(2)*(x(i+2)-x(i-2)) &
              + a7(3)*(x(i+3)-x(i-3))     
     idx6_imax=1.0_wp/idx6_imax
  endif
  if ((BC_face(2,1)%sort<=0).and.(.not.is_curv)) then
     j=4
     idy6_jmin= a7(1)*(y(j+1)-y(j-1)) &
              + a7(2)*(y(j+2)-y(j-2)) &
              + a7(3)*(y(j+3)-y(j-3))     
     idy6_jmin=1.0_wp/idy6_jmin
  endif
  if ((BC_face(2,2)%sort<=0).and.(.not.is_curv)) then
     j=ny-3
     idy6_jmax= a7(1)*(y(j+1)-y(j-1)) &
              + a7(2)*(y(j+2)-y(j-2)) &
              + a7(3)*(y(j+3)-y(j-3))    
     idy6_jmax=1.0_wp/idy6_jmax
  endif
  if (BC_face(3,1)%sort<=0) then
     k=4
     idz6_kmin= a7(1)*(z(k+1)-z(k-1)) &
              + a7(2)*(z(k+2)-z(k-2)) &
              + a7(3)*(z(k+3)-z(k-3))     
     idz6_kmin=1.0_wp/idz6_kmin
  endif
  if (BC_face(3,2)%sort<=0) then
     k=nz-3
     idz6_kmax= a7(1)*(z(k+1)-z(k-1)) &
              + a7(2)*(z(k+2)-z(k-2)) &
              + a7(3)*(z(k+3)-z(k-3))    
     idz6_kmax=1.0_wp/idz6_kmax
  endif
   
  if (ngh==4) return
     
  ! 9-point stencil
  ! ===============
  if ((BC_face(1,1)%sort<=0).and.(.not.is_curv)) then
     i=5
     idx8_imin= a9(1)*(x(i+1)-x(i-1)) &
              + a9(2)*(x(i+2)-x(i-2)) &
              + a9(3)*(x(i+3)-x(i-3)) &
              + a9(4)*(x(i+4)-x(i-4))     
     idx8_imin=1.0_wp/idx8_imin
  endif
  if ((BC_face(1,2)%sort<=0).and.(.not.is_curv)) then
     i=nx-4
     idx8_imax= a9(1)*(x(i+1)-x(i-1)) &
              + a9(2)*(x(i+2)-x(i-2)) &
              + a9(3)*(x(i+3)-x(i-3)) &
              + a9(4)*(x(i+4)-x(i-4))     
     idx8_imax=1.0_wp/idx8_imax
  endif
  if ((BC_face(2,1)%sort<=0).and.(.not.is_curv)) then
     j=5
     idy8_jmin= a9(1)*(y(j+1)-y(j-1)) &
              + a9(2)*(y(j+2)-y(j-2)) &
              + a9(3)*(y(j+3)-y(j-3)) &
              + a9(4)*(y(j+4)-y(j-4))     
     idy8_jmin=1.0_wp/idy8_jmin
  endif
  if ((BC_face(2,2)%sort<=0).and.(.not.is_curv)) then
     j=ny-4
     idy8_jmax= a9(1)*(y(j+1)-y(j-1)) &
              + a9(2)*(y(j+2)-y(j-2)) &
              + a9(3)*(y(j+3)-y(j-3)) &
              + a9(4)*(y(j+4)-y(j-4))     
     idy8_jmax=1.0_wp/idy8_jmax
  endif
  if (BC_face(3,1)%sort<=0) then
     k=5
     idz8_kmin= a9(1)*(z(k+1)-z(k-1)) &
              + a9(2)*(z(k+2)-z(k-2)) &
              + a9(3)*(z(k+3)-z(k-3)) &
              + a9(4)*(z(k+4)-z(k-4))     
     idz8_kmin=1.0_wp/idz8_kmin
  endif
  if (BC_face(3,2)%sort<=0) then
     k=nz-4
     idz8_kmax= a9(1)*(z(k+1)-z(k-1)) &
              + a9(2)*(z(k+2)-z(k-2)) &
              + a9(3)*(z(k+3)-z(k-3)) &
              + a9(4)*(z(k+4)-z(k-4))     
     idz8_kmax=1.0_wp/idz8_kmax
  endif
   
  if (ngh==5) return

end subroutine grid_metrics_wall_cart

!===============================================================================
subroutine grid_metrics_wall_cart_SBP4
!===============================================================================
  !> Compute near-wall Cartesian metrics
  !> version SBP4: Summation by Parts order 4
!===============================================================================
  use mod_grid
  use mod_bc
  use mod_coeff_deriv
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j,k
  ! ---------------------------------------------------------------------------

  ! Point #1: stencil [o x x x]
  ! ===========================
  if ((BC_face(1,1)%sort<=0).and.(.not.is_curv)) then
     !i=1
     idx1_imin= as4p0(1)*x(1)+as4p0(2)*x(2) &
              + as4p0(3)*x(3)+as4p0(4)*x(4)
     idx1_imin=1.0_wp/idx1_imin
  endif
  if ((BC_face(1,2)%sort<=0).and.(.not.is_curv)) then
     !i=nx
     idx1_imax= as4m0(1)*x(nx  )+as4m0(2)*x(nx-1) &
              + as4m0(3)*x(nx-2)+as4m0(4)*x(nx-3)
     idx1_imax=1.0_wp/idx1_imax
  endif
  if ((BC_face(2,1)%sort<=0).and.(.not.is_curv)) then
     !j=1
     idy1_jmin= as4p0(1)*y(1)+as4p0(2)*y(2) &
              + as4p0(3)*y(3)+as4p0(4)*y(4)
     idy1_jmin=1.0_wp/idy1_jmin
  endif
  if ((BC_face(2,2)%sort<=0).and.(.not.is_curv)) then
     !j=ny
     idy1_jmax= as4m0(1)*y(ny  )+as4m0(2)*y(ny-1) &
              + as4m0(3)*y(ny-2)+as4m0(4)*y(ny-3)
     idy1_jmax=1.0_wp/idy1_jmax
  endif
  if (BC_face(3,1)%sort<=0) then
     !k=1
     idz1_kmin= as4p0(1)*z(1)+as4p0(2)*z(2) &
              + as4p0(3)*z(3)+as4p0(4)*z(4)
     idz1_kmin=1.0_wp/idz1_kmin
  endif
  if (BC_face(3,2)%sort<=0) then
     !k=nz
     idz1_kmax= as4m0(1)*z(nz  )+as4m0(2)*z(nz-1) &
              + as4m0(3)*z(nz-2)+as4m0(4)*z(nz-3)
     idz1_kmax=1.0_wp/idz1_kmax
  endif

  ! Point #2: stencil [x o x x]
  ! ===========================
  if ((BC_face(1,1)%sort<=0).and.(.not.is_curv)) then
     !i=2
     idx2_imin= as4p1(1)*x(1)+as4p1(2)*x(2) &
              + as4p1(3)*x(3)+as4p1(4)*x(4)
     idx2_imin=1.0_wp/idx2_imin
  endif
  if ((BC_face(1,2)%sort<=0).and.(.not.is_curv)) then
     !i=nx-1
     idx2_imax= as4m1(1)*x(nx  )+as4m1(2)*x(nx-1) &
              + as4m1(3)*x(nx-2)+as4m1(4)*x(nx-3)
     idx2_imax=1.0_wp/idx2_imax
  endif
  if ((BC_face(2,1)%sort<=0).and.(.not.is_curv)) then
     !j=2
     idy2_jmin= as4p1(1)*y(1)+as4p1(2)*y(2) &
              + as4p1(3)*y(3)+as4p1(4)*y(4)
     idy2_jmin=1.0_wp/idy2_jmin
  endif
  if ((BC_face(2,2)%sort<=0).and.(.not.is_curv)) then
     !j=ny-1
     idy2_jmax= as4m1(1)*y(ny  )+as4m1(2)*y(ny-1) &
              + as4m1(3)*y(ny-2)+as4m1(4)*y(ny-3)
     idy2_jmax=1.0_wp/idy2_jmax
  endif
  if (BC_face(3,1)%sort<=0) then
     !k=2
     idz2_kmin= as4p1(1)*z(1)+as4p1(2)*z(2) &
              + as4p1(3)*z(3)+as4p1(4)*z(4)
     idz2_kmin=1.0_wp/idz2_kmin
  endif
  if (BC_face(3,2)%sort<=0) then
     !k=nz-1
     idz2_kmax= as4m1(1)*z(nz  )+as4m1(2)*z(nz-1) &
              + as4m1(3)*z(nz-2)+as4m1(4)*z(nz-3)
     idz2_kmax=1.0_wp/idz2_kmax
  endif
   
  ! Point #3: stencil [x x o x x]
  ! =============================
  if ((BC_face(1,1)%sort<=0).and.(.not.is_curv)) then
     !i=3
     idx4_imin= as4p2(1)*x(1)+as4p2(2)*x(2) &
              + as4p2(3)*x(3)+as4p2(4)*x(4) &
              + as4p2(5)*x(5)
     idx4_imin=1.0_wp/idx4_imin
  endif
  if ((BC_face(1,2)%sort<=0).and.(.not.is_curv)) then
     !i=nx-2
     idx4_imax= as4m2(1)*x(nx  )+as4m2(2)*x(nx-1) &
              + as4m2(3)*x(nx-2)+as4m2(4)*x(nx-3) &
              + as4m2(5)*x(nx-4)
     idx4_imax=1.0_wp/idx4_imax
  endif
  if ((BC_face(2,1)%sort<=0).and.(.not.is_curv)) then
     !j=3
     idy4_jmin= as4p2(1)*y(1)+as4p2(2)*y(2) &
              + as4p2(3)*y(3)+as4p2(4)*y(4) &
              + as4p2(5)*y(5)
     idy4_jmin=1.0_wp/idy4_jmin
  endif
  if ((BC_face(2,2)%sort<=0).and.(.not.is_curv)) then
     !j=ny-2
     idy4_jmax= as4m2(1)*y(ny  )+as4m2(2)*y(ny-1) &
              + as4m2(3)*y(ny-2)+as4m2(4)*y(ny-3) &
              + as4m2(5)*y(ny-4)
     idy4_jmax=1.0_wp/idy4_jmax
  endif
  if (BC_face(3,1)%sort<=0) then
     !k=3
     idz4_kmin= as4p2(1)*z(1)+as4p2(2)*z(2) &
              + as4p2(3)*z(3)+as4p2(4)*z(4) &
              + as4p2(5)*z(5)
     idz4_kmin=1.0_wp/idz4_kmin
  endif
  if (BC_face(3,2)%sort<=0) then
     !k=nz-2
     idz4_kmax= as4m2(1)*z(nz  )+as4m2(2)*z(nz-1) &
              + as4m2(3)*z(nz-2)+as4m2(4)*z(nz-3) &
              + as4m2(5)*z(nz-4)
     idz4_kmax=1.0_wp/idz4_kmax
  endif
   
  ! Point #4: stencil [x x x o x x]
  ! ===============================
  if ((BC_face(1,1)%sort<=0).and.(.not.is_curv)) then
     !i=4
     idx6_imin= as4p3(1)*x(1)+as4p3(2)*x(2) &
              + as4p3(3)*x(3)+as4p3(4)*x(4) &
              + as4p3(5)*x(5)+as4p3(6)*x(6)
     idx6_imin=1.0_wp/idx6_imin
  endif
  if ((BC_face(1,2)%sort<=0).and.(.not.is_curv)) then
     !i=nx-3
     idx6_imax= as4m3(1)*x(nx  )+as4m3(2)*x(nx-1) &
              + as4m3(3)*x(nx-2)+as4m3(4)*x(nx-3) &
              + as4m3(5)*x(nx-4)+as4m3(6)*x(nx-5)
     idx6_imax=1.0_wp/idx6_imax
  endif
  if ((BC_face(2,1)%sort<=0).and.(.not.is_curv)) then
     !j=4
     idy6_jmin= as4p3(1)*y(1)+as4p3(2)*y(2) &
              + as4p3(3)*y(3)+as4p3(4)*y(4) &
              + as4p3(5)*y(5)+as4p3(6)*y(6)
     idy6_jmin=1.0_wp/idy6_jmin
  endif
  if ((BC_face(2,2)%sort<=0).and.(.not.is_curv)) then
     !j=ny-3
     idy6_jmax= as4m3(1)*y(ny  )+as4m3(2)*y(ny-1) &
              + as4m3(3)*y(ny-2)+as4m3(4)*y(ny-3) &
              + as4m3(5)*y(ny-4)+as4m3(6)*y(ny-5)
     idy6_jmax=1.0_wp/idy6_jmax
  endif
  if (BC_face(3,1)%sort<=0) then
     !k=4
     idz6_kmin= as4p3(1)*z(1)+as4p3(2)*z(2) &
              + as4p3(3)*z(3)+as4p3(4)*z(4) &
              + as4p3(5)*z(5)+as4p3(6)*z(6)
     idz6_kmin=1.0_wp/idz6_kmin
  endif
  if (BC_face(3,2)%sort<=0) then
     !k=nz-3
     idz6_kmax= as4m3(1)*z(nz  )+as4m3(2)*z(nz-1) &
              + as4m3(3)*z(nz-2)+as4m3(4)*z(nz-3) &
              + as4m3(5)*z(nz-4)+as4m3(6)*z(nz-5)
     idz6_kmax=1.0_wp/idz6_kmax
  endif
   
  if (ngh==4) return
     
  ! 9-point stencil
  ! ===============
  if ((BC_face(1,1)%sort<=0).and.(.not.is_curv)) then
     i=5
     idx8_imin= a9(1)*(x(i+1)-x(i-1)) &
              + a9(2)*(x(i+2)-x(i-2)) &
              + a9(3)*(x(i+3)-x(i-3)) &
              + a9(4)*(x(i+4)-x(i-4))     
     idx8_imin=1.0_wp/idx8_imin
  endif
  if ((BC_face(1,2)%sort<=0).and.(.not.is_curv)) then
     i=nx-4
     idx8_imax= a9(1)*(x(i+1)-x(i-1)) &
              + a9(2)*(x(i+2)-x(i-2)) &
              + a9(3)*(x(i+3)-x(i-3)) &
              + a9(4)*(x(i+4)-x(i-4))     
     idx8_imax=1.0_wp/idx8_imax
  endif
  if ((BC_face(2,1)%sort<=0).and.(.not.is_curv)) then
     j=5
     idy8_jmin= a9(1)*(y(j+1)-y(j-1)) &
              + a9(2)*(y(j+2)-y(j-2)) &
              + a9(3)*(y(j+3)-y(j-3)) &
              + a9(4)*(y(j+4)-y(j-4))     
     idy8_jmin=1.0_wp/idy8_jmin
  endif
  if ((BC_face(2,2)%sort<=0).and.(.not.is_curv)) then
     j=ny-4
     idy8_jmax= a9(1)*(y(j+1)-y(j-1)) &
              + a9(2)*(y(j+2)-y(j-2)) &
              + a9(3)*(y(j+3)-y(j-3)) &
              + a9(4)*(y(j+4)-y(j-4))     
     idy8_jmax=1.0_wp/idy8_jmax
  endif
  if (BC_face(3,1)%sort<=0) then
     k=5
     idz8_kmin= a9(1)*(z(k+1)-z(k-1)) &
              + a9(2)*(z(k+2)-z(k-2)) &
              + a9(3)*(z(k+3)-z(k-3)) &
              + a9(4)*(z(k+4)-z(k-4))     
     idz8_kmin=1.0_wp/idz8_kmin
  endif
  if (BC_face(3,2)%sort<=0) then
     k=nz-4
     idz8_kmax= a9(1)*(z(k+1)-z(k-1)) &
              + a9(2)*(z(k+2)-z(k-2)) &
              + a9(3)*(z(k+3)-z(k-3)) &
              + a9(4)*(z(k+4)-z(k-4))     
     idz8_kmax=1.0_wp/idz8_kmax
  endif
   
  if (ngh==5) return
  
end subroutine grid_metrics_wall_cart_SBP4
