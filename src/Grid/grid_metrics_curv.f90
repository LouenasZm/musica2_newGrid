!===============================================================================
subroutine grid_metrics_ksi
!===============================================================================
  !> Compute curvilinear metrics along ksi 
  !> - stencil -ngh:+ngh for inviscid fluxes -
!===============================================================================
  use mod_coeff_deriv
  use mod_grid
  use mod_bc
  use warnstop
  use mod_mpi !! just for check
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j,l
  real(wp), dimension(-ngh:ngh) :: a_ ! interior DF scheme
  ! ---------------------------------------------------------------------------

  ! Initializations
  ! ---------------
  allocate(x_ksi(nx1:nx2,ny1:ny2),y_ksi(nx1:nx2,ny1:ny2))
  x_ksi=0.0_wp
  y_ksi=0.0_wp

  ! BC imin
  ! -------
  if (BC_face(1,1)%sort<=0) then
     select case (ngh)
     case (5)
        i=1
        do j=ndyt,nfyt
           do l=0,10
              x_ksi(i,j)=x_ksi(i,j)+a010(1+l)*xc(i+l,j)
              y_ksi(i,j)=y_ksi(i,j)+a010(1+l)*yc(i+l,j)
           enddo
        enddo
        i=2
        do j=ndyt,nfyt
           do l=-1,9
              x_ksi(i,j)=x_ksi(i,j)+a19(2+l)*xc(i+l,j)
              y_ksi(i,j)=y_ksi(i,j)+a19(2+l)*yc(i+l,j)
           enddo
        enddo
        i=3
        do j=ndyt,nfyt
           do l=-2,8
              x_ksi(i,j)=x_ksi(i,j)+a28(3+l)*xc(i+l,j)
              y_ksi(i,j)=y_ksi(i,j)+a28(3+l)*yc(i+l,j)
           enddo
        enddo
        i=4
        do j=ndyt,nfyt
           do l=-3,7
              x_ksi(i,j)=x_ksi(i,j)+a37(4+l)*xc(i+l,j)
              y_ksi(i,j)=y_ksi(i,j)+a37(4+l)*yc(i+l,j)
           enddo
        enddo
        i=5
        do j=ndyt,nfyt
           do l=-4,6
              x_ksi(i,j)=x_ksi(i,j)+a46(5+l)*xc(i+l,j)
              y_ksi(i,j)=y_ksi(i,j)+a46(5+l)*yc(i+l,j)
           enddo
        enddo
     case (4)
        i=1
        do j=ndyt,nfyt
           do l=0,8
              x_ksi(i,j)=x_ksi(i,j)+a08(1+l)*xc(i+l,j)
              y_ksi(i,j)=y_ksi(i,j)+a08(1+l)*yc(i+l,j)
           enddo
        enddo
        i=2
        do j=ndyt,nfyt
           do l=-1,7
              x_ksi(i,j)=x_ksi(i,j)+a17(2+l)*xc(i+l,j)
              y_ksi(i,j)=y_ksi(i,j)+a17(2+l)*yc(i+l,j)
           enddo
        enddo
        i=3
        do j=ndyt,nfyt
           do l=-2,6
              x_ksi(i,j)=x_ksi(i,j)+a26(3+l)*xc(i+l,j)
              y_ksi(i,j)=y_ksi(i,j)+a26(3+l)*yc(i+l,j)
           enddo
        enddo
        i=4
        do j=ndyt,nfyt
           do l=-3,5
              x_ksi(i,j)=x_ksi(i,j)+a35(4+l)*xc(i+l,j)
              y_ksi(i,j)=y_ksi(i,j)+a35(4+l)*yc(i+l,j)
           enddo
        enddo
     case (3)
        i=1
        do j=ndyt,nfyt
           do l=0,6
              x_ksi(i,j)=x_ksi(i,j)+a06(1+l)*xc(i+l,j)
              y_ksi(i,j)=y_ksi(i,j)+a06(1+l)*yc(i+l,j)
           enddo
        enddo
        i=2
        do j=ndyt,nfyt
           do l=-1,5
              x_ksi(i,j)=x_ksi(i,j)+a15(2+l)*xc(i+l,j)
              y_ksi(i,j)=y_ksi(i,j)+a15(2+l)*yc(i+l,j)
           enddo
        enddo
        i=3
        do j=ndyt,nfyt
           do l=-2,4
              x_ksi(i,j)=x_ksi(i,j)+a24(3+l)*xc(i+l,j)
              y_ksi(i,j)=y_ksi(i,j)+a24(3+l)*yc(i+l,j)
           enddo
        enddo
     case (2)
        i=1
        do j=ndyt,nfyt
           do l=0,4
              x_ksi(i,j)=x_ksi(i,j)+a04(1+l)*xc(i+l,j)
              y_ksi(i,j)=y_ksi(i,j)+a04(1+l)*yc(i+l,j)
           enddo
        enddo
        i=2
        do j=ndyt,nfyt
           do l=-1,3
              x_ksi(i,j)=x_ksi(i,j)+a13(2+l)*xc(i+l,j)
              y_ksi(i,j)=y_ksi(i,j)+a13(2+l)*yc(i+l,j)
           enddo
        enddo
     case (1)
        i=1
        do j=ndyt,nfyt
           do l=0,2
              x_ksi(i,j)=x_ksi(i,j)+a02(1+l)*xc(i+l,j)
              y_ksi(i,j)=y_ksi(i,j)+a02(1+l)*yc(i+l,j)
           enddo
        enddo
        !do j=ndyt,nfyt
        !   x_ksi(1,j)=xc(2,j)-xc(1,j)
        !   y_ksi(1,j)=yc(2,j)-yc(1,j)
        !enddo
     case default
        call mpistop('max ngh is 5 (11pts-stencil)!', 0)
     end select
  endif
  
  ! BC imax
  ! -------
  if (BC_face(1,2)%sort<=0) then
     select case (ngh)
     case (5)
        i=nx-4
        do j=ndyt,nfyt
           do l=-6,4
              x_ksi(i,j)=x_ksi(i,j)-a46(5-l)*xc(i+l,j)
              y_ksi(i,j)=y_ksi(i,j)-a46(5-l)*yc(i+l,j)
           enddo
        enddo
        i=nx-3
        do j=ndyt,nfyt
           do l=-7,3
              x_ksi(i,j)=x_ksi(i,j)-a37(4-l)*xc(i+l,j)
              y_ksi(i,j)=y_ksi(i,j)-a37(4-l)*yc(i+l,j)
           enddo
        enddo
        i=nx-2
        do j=ndyt,nfyt
           do l=-8,2
              x_ksi(i,j)=x_ksi(i,j)-a28(3-l)*xc(i+l,j)
              y_ksi(i,j)=y_ksi(i,j)-a28(3-l)*yc(i+l,j)
           enddo
        enddo
        i=nx-1
        do j=ndyt,nfyt
           do l=-9,1
              x_ksi(i,j)=x_ksi(i,j)-a19(2-l)*xc(i+l,j)
              y_ksi(i,j)=y_ksi(i,j)-a19(2-l)*yc(i+l,j)
           enddo
        enddo
        i=nx
        do j=ndyt,nfyt
           do l=-10,0
              x_ksi(i,j)=x_ksi(i,j)-a010(1-l)*xc(i+l,j)
              y_ksi(i,j)=y_ksi(i,j)-a010(1-l)*yc(i+l,j)
           enddo
        enddo
     case (4)
        i=nx-3
        do j=ndyt,nfyt
           do l=-5,3
              x_ksi(i,j)=x_ksi(i,j)-a35(4-l)*xc(i+l,j)
              y_ksi(i,j)=y_ksi(i,j)-a35(4-l)*yc(i+l,j)
           enddo
        enddo
        i=nx-2
        do j=ndyt,nfyt
           do l=-6,2
              x_ksi(i,j)=x_ksi(i,j)-a26(3-l)*xc(i+l,j)
              y_ksi(i,j)=y_ksi(i,j)-a26(3-l)*yc(i+l,j)
           enddo
        enddo
        i=nx-1
        do j=ndyt,nfyt
           do l=-7,1
              x_ksi(i,j)=x_ksi(i,j)-a17(2-l)*xc(i+l,j)
              y_ksi(i,j)=y_ksi(i,j)-a17(2-l)*yc(i+l,j)
           enddo
        enddo
        i=nx
        do j=ndyt,nfyt
           do l=-8,0
              x_ksi(i,j)=x_ksi(i,j)-a08(1-l)*xc(i+l,j)
              y_ksi(i,j)=y_ksi(i,j)-a08(1-l)*yc(i+l,j)
           enddo
        enddo
     case (3)
        i=nx-2
        do j=ndyt,nfyt
           do l=-4,2
              x_ksi(i,j)=x_ksi(i,j)-a24(3-l)*xc(i+l,j)
              y_ksi(i,j)=y_ksi(i,j)-a24(3-l)*yc(i+l,j)
           enddo
        enddo
        i=nx-1
        do j=ndyt,nfyt
           do l=-5,1
              x_ksi(i,j)=x_ksi(i,j)-a15(2-l)*xc(i+l,j)
              y_ksi(i,j)=y_ksi(i,j)-a15(2-l)*yc(i+l,j)
           enddo
        enddo
        i=nx
        do j=ndyt,nfyt
           do l=-6,0
              x_ksi(i,j)=x_ksi(i,j)-a06(1-l)*xc(i+l,j)
              y_ksi(i,j)=y_ksi(i,j)-a06(1-l)*yc(i+l,j)
           enddo
        enddo
     case (2)
        i=nx-1
        do j=ndyt,nfyt
           do l=-3,1
              x_ksi(i,j)=x_ksi(i,j)-a13(2-l)*xc(i+l,j)
              y_ksi(i,j)=y_ksi(i,j)-a13(2-l)*yc(i+l,j)
           enddo
        enddo
        i=nx
        do j=ndyt,nfyt
           do l=-4,0
              x_ksi(i,j)=x_ksi(i,j)-a04(1-l)*xc(i+l,j)
              y_ksi(i,j)=y_ksi(i,j)-a04(1-l)*yc(i+l,j)
           enddo
        enddo
     case (1)
        i=nx
        do j=ndyt,nfyt
           do l=-2,0
              x_ksi(i,j)=x_ksi(i,j)-a02(1-l)*xc(i+l,j)
              y_ksi(i,j)=y_ksi(i,j)-a02(1-l)*yc(i+l,j)
           enddo
        enddo
        !do j=ndyt,nfyt
        !   x_ksi(nx,j)=xc(nx,j)-xc(nx-1,j)
        !   y_ksi(nx,j)=yc(nx,j)-yc(nx-1,j)
        !enddo
     case default
        call mpistop('max ngh is 5 (11pts-stencil)!', 0)
     end select
  endif

  ! Interior points
  ! ---------------
  select case (ngh)
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
     call mpistop('max ngh is 5 (11pts-stencil)!', 0)
  end select
 
  do i=ndx,nfx
     do j=ndyt,nfyt
        do l=-ngh,ngh
           x_ksi(i,j)= x_ksi(i,j) + a_(l)*xc(i+l,j)
           y_ksi(i,j)= y_ksi(i,j) + a_(l)*yc(i+l,j)
        enddo
     enddo
  enddo
        
end subroutine grid_metrics_ksi

!===============================================================================
subroutine grid_metrics_eta
!===============================================================================
  !> Compute curvilinear metrics along eta 
  !> - stencil -ngh:+ngh for inviscid fluxes -
!===============================================================================
  use mod_coeff_deriv
  use mod_grid
  use mod_bc
  use warnstop
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j,l
  real(wp), dimension(-ngh:ngh) :: a_ ! interior DF scheme
  ! ---------------------------------------------------------------------------

  ! Initializations
  ! ---------------
  allocate(x_eta(nx1:nx2,ny1:ny2),y_eta(nx1:nx2,ny1:ny2))
  x_eta=0.0_wp
  y_eta=0.0_wp

  ! BC jmin
  ! -------
  if (BC_face(2,1)%sort<=0) then
     select case (ngh)
     case (5)
        j=1
        do i=ndxt,nfxt
           do l=0,10
              x_eta(i,j)=x_eta(i,j)+a010(1+l)*xc(i,j+l)
              y_eta(i,j)=y_eta(i,j)+a010(1+l)*yc(i,j+l)
           enddo
        enddo
        j=2
        do i=ndxt,nfxt
           do l=-1,9
              x_eta(i,j)=x_eta(i,j)+a19(2+l)*xc(i,j+l)
              y_eta(i,j)=y_eta(i,j)+a19(2+l)*yc(i,j+l)
           enddo
        enddo
        j=3
        do i=ndxt,nfxt
           do l=-2,8
              x_eta(i,j)=x_eta(i,j)+a28(3+l)*xc(i,j+l)
              y_eta(i,j)=y_eta(i,j)+a28(3+l)*yc(i,j+l)
           enddo
        enddo
        j=4
        do i=ndxt,nfxt
           do l=-3,7
              x_eta(i,j)=x_eta(i,j)+a37(4+l)*xc(i,j+l)
              y_eta(i,j)=y_eta(i,j)+a37(4+l)*yc(i,j+l)
           enddo
        enddo
        j=5
        do i=ndxt,nfxt
           do l=-4,6
              x_eta(i,j)=x_eta(i,j)+a46(5+l)*xc(i,j+l)
              y_eta(i,j)=y_eta(i,j)+a46(5+l)*yc(i,j+l)
           enddo
        enddo
     case (4)
        j=1
        do i=ndxt,nfxt
           do l=0,8
              x_eta(i,j)=x_eta(i,j)+a08(1+l)*xc(i,j+l)
              y_eta(i,j)=y_eta(i,j)+a08(1+l)*yc(i,j+l)
           enddo
        enddo
        j=2
        do i=ndxt,nfxt
           do l=-1,7
              x_eta(i,j)=x_eta(i,j)+a17(2+l)*xc(i,j+l)
              y_eta(i,j)=y_eta(i,j)+a17(2+l)*yc(i,j+l)
           enddo
        enddo
        j=3
        do i=ndxt,nfxt
           do l=-2,6
              x_eta(i,j)=x_eta(i,j)+a26(3+l)*xc(i,j+l)
              y_eta(i,j)=y_eta(i,j)+a26(3+l)*yc(i,j+l)
           enddo
        enddo
        j=4
        do i=ndxt,nfxt
           do l=-3,5
              x_eta(i,j)=x_eta(i,j)+a35(4+l)*xc(i,j+l)
              y_eta(i,j)=y_eta(i,j)+a35(4+l)*yc(i,j+l)
           enddo
        enddo
     case (3)
        j=1
        do i=ndxt,nfxt
           do l=0,6
              x_eta(i,j)=x_eta(i,j)+a06(1+l)*xc(i,j+l)
              y_eta(i,j)=y_eta(i,j)+a06(1+l)*yc(i,j+l)
           enddo
        enddo
        j=2
        do i=ndxt,nfxt
           do l=-1,5
              x_eta(i,j)=x_eta(i,j)+a15(2+l)*xc(i,j+l)
              y_eta(i,j)=y_eta(i,j)+a15(2+l)*yc(i,j+l)
           enddo
        enddo
        j=3
        do i=ndxt,nfxt
           do l=-2,4
              x_eta(i,j)=x_eta(i,j)+a24(3+l)*xc(i,j+l)
              y_eta(i,j)=y_eta(i,j)+a24(3+l)*yc(i,j+l)
           enddo
        enddo
     case (2)
        j=1
        do i=ndxt,nfxt
           do l=0,4
              x_eta(i,j)=x_eta(i,j)+a04(1+l)*xc(i,j+l)
              y_eta(i,j)=y_eta(i,j)+a04(1+l)*yc(i,j+l)
           enddo
        enddo
        j=2
        do i=ndxt,nfxt
           do l=-1,3
              x_eta(i,j)=x_eta(i,j)+a13(2+l)*xc(i,j+l)
              y_eta(i,j)=y_eta(i,j)+a13(2+l)*yc(i,j+l)
           enddo
        enddo
     case (1)
        j=1
        do i=ndxt,nfxt
           do l=0,2
              x_eta(i,j)=x_eta(i,j)+a02(1+l)*xc(i,j+l)
              y_eta(i,j)=y_eta(i,j)+a02(1+l)*yc(i,j+l)
           enddo
        enddo
        !do i=ndxt,nfxt
        !   x_eta(i,1)=xc(i,2)-xc(i,1)
        !   y_eta(i,1)=yc(i,2)-yc(i,1)
        !enddo
     case default
        call mpistop('max ngh is 5 (11pts-stencil)!', 0)
     end select
  endif

  ! BC jmax
  ! -------
  if (BC_face(2,2)%sort<=0) then
     select case (ngh)
     case (5)
        j=ny-4
        do i=ndxt,nfxt
           do l=-6,4
              x_eta(i,j)=x_eta(i,j)-a46(5-l)*xc(i,j+l)
              y_eta(i,j)=y_eta(i,j)-a46(5-l)*yc(i,j+l)
           enddo
        enddo
        j=ny-3
        do i=ndxt,nfxt
           do l=-7,3
              x_eta(i,j)=x_eta(i,j)-a37(4-l)*xc(i,j+l)
              y_eta(i,j)=y_eta(i,j)-a37(4-l)*yc(i,j+l)
           enddo
        enddo
        j=ny-2
        do i=ndxt,nfxt
           do l=-8,2
              x_eta(i,j)=x_eta(i,j)-a28(3-l)*xc(i,j+l)
              y_eta(i,j)=y_eta(i,j)-a28(3-l)*yc(i,j+l)
           enddo
        enddo
        j=ny-1
        do i=ndxt,nfxt
           do l=-9,1
              x_eta(i,j)=x_eta(i,j)-a19(2-l)*xc(i,j+l)
              y_eta(i,j)=y_eta(i,j)-a19(2-l)*yc(i,j+l)
           enddo
        enddo
        j=ny
        do i=ndxt,nfxt
           do l=-10,0
              x_eta(i,j)=x_eta(i,j)-a010(1-l)*xc(i,j+l)
              y_eta(i,j)=y_eta(i,j)-a010(1-l)*yc(i,j+l)
           enddo
        enddo
     case (4)
        j=ny-3
        do i=ndxt,nfxt
           do l=-5,3
              x_eta(i,j)=x_eta(i,j)-a35(4-l)*xc(i,j+l)
              y_eta(i,j)=y_eta(i,j)-a35(4-l)*yc(i,j+l)
           enddo
        enddo
        j=ny-2
        do i=ndxt,nfxt
           do l=-6,2
              x_eta(i,j)=x_eta(i,j)-a26(3-l)*xc(i,j+l)
              y_eta(i,j)=y_eta(i,j)-a26(3-l)*yc(i,j+l)
           enddo
        enddo
        j=ny-1
        do i=ndxt,nfxt
           do l=-7,1
              x_eta(i,j)=x_eta(i,j)-a17(2-l)*xc(i,j+l)
              y_eta(i,j)=y_eta(i,j)-a17(2-l)*yc(i,j+l)
           enddo
        enddo
        j=ny
        do i=ndxt,nfxt
           do l=-8,0
              x_eta(i,j)=x_eta(i,j)-a08(1-l)*xc(i,j+l)
              y_eta(i,j)=y_eta(i,j)-a08(1-l)*yc(i,j+l)
           enddo
        enddo
     case (3)
        j=ny-2
        do i=ndxt,nfxt
           do l=-4,2
              x_eta(i,j)=x_eta(i,j)-a24(3-l)*xc(i,j+l)
              y_eta(i,j)=y_eta(i,j)-a24(3-l)*yc(i,j+l)
           enddo
        enddo
        j=ny-1
        do i=ndxt,nfxt
           do l=-5,1
              x_eta(i,j)=x_eta(i,j)-a15(2-l)*xc(i,j+l)
              y_eta(i,j)=y_eta(i,j)-a15(2-l)*yc(i,j+l)
           enddo
        enddo
        j=ny
        do i=ndxt,nfxt
           do l=-6,0
              x_eta(i,j)=x_eta(i,j)-a06(1-l)*xc(i,j+l)
              y_eta(i,j)=y_eta(i,j)-a06(1-l)*yc(i,j+l)
           enddo
        enddo
     case (2)
        j=ny-1
        do i=ndxt,nfxt
           do l=-3,1
              x_eta(i,j)=x_eta(i,j)-a13(2-l)*xc(i,j+l)
              y_eta(i,j)=y_eta(i,j)-a13(2-l)*yc(i,j+l)
           enddo
        enddo
        j=ny
        do i=ndxt,nfxt
           do l=-4,0
              x_eta(i,j)=x_eta(i,j)-a04(1-l)*xc(i,j+l)
              y_eta(i,j)=y_eta(i,j)-a04(1-l)*yc(i,j+l)
           enddo
        enddo
     case (1)
        j=ny
        do i=ndxt,nfxt
           do l=-2,0
              x_eta(i,j)=x_eta(i,j)-a02(1-l)*xc(i,j+l)
              y_eta(i,j)=y_eta(i,j)-a02(1-l)*yc(i,j+l)
           enddo
        enddo
        !do i=ndxt,nfxt
        !   x_eta(i,ny)=xc(i,ny)-xc(i,ny-1)
        !   y_eta(i,ny)=yc(i,ny)-yc(i,ny-1)
        !enddo
     case default
        call mpistop('max ngh is 5 (11pts-stencil)!', 0)
     end select
  endif
 
  ! Interior points
  ! ---------------
  select case (ngh)
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
     call mpistop('max ngh is 5 (11pts-stencil)!', 0)
  end select
  
  do j=ndy,nfy
     do i=ndxt,nfxt
        do l=-ngh,ngh
           x_eta(i,j)= x_eta(i,j) + a_(l)*xc(i,j+l)
           y_eta(i,j)= y_eta(i,j) + a_(l)*yc(i,j+l)
        enddo
     enddo
  enddo

end subroutine grid_metrics_eta

!===============================================================================
subroutine grid_metrics_ijacob
!===============================================================================
  !> Compute Jacobians of metrics
  !> - stencil -ngh:+ngh for inviscid fluxes -
!===============================================================================
  use mod_grid
  use mod_mpi ! for temporary check
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j
  ! ---------------------------------------------------------------------------

!!$  allocate(ijacob(ndx1:nfx1,ndy1:nfy1))
!!$
!!$  ! Jacobian of metrics transformation
!!$  ! ----------------------------------
!!$  ijacob=1.0_wp
!!$  
!!$  do i=ndx1,nfx1
!!$     do j=1,ny
!!$        ijacob(i,j)=x_ksi(i,j)*y_eta(i,j)-x_eta(i,j)*y_ksi(i,j)
!!$        if (ijacob(i,j)==0) print *,'jac1',iproc,i,j
!!$     enddo
!!$  enddo
!!$  
!!$  do i=1,nx
!!$     do j=ndy1,nfy1
!!$        ijacob(i,j)=x_ksi(i,j)*y_eta(i,j)-x_eta(i,j)*y_ksi(i,j)
!!$        if (ijacob(i,j)==0) print *,'jac2',iproc,i,j
!!$     enddo
!!$  enddo
!!$
!!$  ! Inverse of Jacobian
!!$  ! -------------------
!!$  ijacob=1.0_wp/ijacob

  allocate(ijacob(0:nfx1,0:nfy1))

  ! Jacobian of metrics transformation
  ! ----------------------------------
  ijacob=1.0_wp
  
  do i=ndx1,nfx1
     do j=1,ny
        ijacob(i,j)=x_ksi(i,j)*y_eta(i,j)-x_eta(i,j)*y_ksi(i,j)
        if (ijacob(i,j)==0) print *,'jac1',iproc,i,j
     enddo
  enddo
  
  if (ndx1==1) then
     ijacob(0,:)=ijacob(1,:)
  endif
  
  do i=1,nx
     do j=ndy1,nfy1
        ijacob(i,j)=x_ksi(i,j)*y_eta(i,j)-x_eta(i,j)*y_ksi(i,j)
        if (ijacob(i,j)==0) print *,'jac2',iproc,i,j
     enddo
  enddo

  if (ndy1==1) then
     ijacob(:,0)=ijacob(:,1)
  endif

  ! Inverse of Jacobian
  ! -------------------
  ijacob=1.0_wp/ijacob

end subroutine grid_metrics_ijacob

!===============================================================================
subroutine grid_metrics_gradients
!===============================================================================
  !> Compute norms of metrics gradients
  !> - stencil -ngh:+ngh for inviscid fluxes -
!===============================================================================
  use mod_grid
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j
  ! ---------------------------------------------------------------------------

!!$  !allocate(g_ksi(ndx1:nfx1,1:ny),g_eta(1:nx,ndy1:nfy1))
!!$
!!$  allocate(g_ksi(ndx1:nfx1,ndy1:nfy1),g_eta(ndx1:nfx1,ndy1:nfy1))
!!$
!!$  ! Norm of metrics gradients
!!$  ! -------------------------
!!$  ! sqrt(ksi_x^2+ksi_y^2)=sqrt[grad(ksi).grad(ksi)]
!!$  !                      =sqrt(x_eta^2+y_eta^2)/|J|
!!$  do i=ndx1,nfx1
!!$     do j=1,ny
!!$        g_ksi(i,j)=sqrt(x_eta(i,j)**2+y_eta(i,j)**2)*abs(ijacob(i,j))
!!$        g_eta(i,j)=sqrt(x_ksi(i,j)**2+y_ksi(i,j)**2)*abs(ijacob(i,j))
!!$     enddo
!!$  enddo
!!$
!!$  ! sqrt(eta_x^2+eta_y^2)=sqrt[grad(eta).grad(eta)]
!!$  !                      =sqrt(x_ksi^2+y_ksi^2)/|J|
!!$  do i=1,nx
!!$     do j=ndy1,nfy1
!!$        g_ksi(i,j)=sqrt(x_eta(i,j)**2+y_eta(i,j)**2)*abs(ijacob(i,j))
!!$        g_eta(i,j)=sqrt(x_ksi(i,j)**2+y_ksi(i,j)**2)*abs(ijacob(i,j))
!!$     enddo
!!$  enddo
  
  allocate(g_ksi(0:nfx1,0:nfy1),g_eta(0:nfx1,0:nfy1))

  ! Norm of metrics gradients
  ! -------------------------
  ! sqrt(ksi_x^2+ksi_y^2)=sqrt[grad(ksi).grad(ksi)]
  !                      =sqrt(x_eta^2+y_eta^2)/|J|
  do i=ndx1,nfx1
     do j=1,ny
        g_ksi(i,j)=sqrt(x_eta(i,j)**2+y_eta(i,j)**2)*abs(ijacob(i,j))
        g_eta(i,j)=sqrt(x_ksi(i,j)**2+y_ksi(i,j)**2)*abs(ijacob(i,j))
     enddo
  enddo

  if (ndx1==1) then
     g_ksi(0,:)=g_ksi(1,:)
     g_eta(0,:)=g_eta(1,:)
  endif

  ! sqrt(eta_x^2+eta_y^2)=sqrt[grad(eta).grad(eta)]
  !                      =sqrt(x_ksi^2+y_ksi^2)/|J|
  do i=1,nx
     do j=ndy1,nfy1
        g_ksi(i,j)=sqrt(x_eta(i,j)**2+y_eta(i,j)**2)*abs(ijacob(i,j))
        g_eta(i,j)=sqrt(x_ksi(i,j)**2+y_ksi(i,j)**2)*abs(ijacob(i,j))
     enddo
  enddo
  
  if (ndy1==1) then
     g_ksi(:,0)=g_ksi(:,1)
     g_eta(:,0)=g_eta(:,1)
  endif
  
end subroutine grid_metrics_gradients

!===============================================================================
subroutine grid_metrics_ksi_v
!===============================================================================
  !> Compute curvilinear metrics along ksi 
  !> - stencil -ngh_v:+ngh_v for viscous fluxes -
!===============================================================================
  use mod_coeff_deriv
  use mod_grid
  use mod_bc
  use warnstop
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j,l
  integer :: ndxv,nfxv
  real(wp), dimension(-ngh_v:ngh_v) :: a_ ! interior DF scheme
  ! ---------------------------------------------------------------------------

  ! Initializations
  ! ---------------
  allocate(x_ksi_v(nx1_v:nx2_v,ny1_v:ny2_v),y_ksi_v(nx1_v:nx2_v,ny1_v:ny2_v))
  x_ksi_v=0.0_wp
  y_ksi_v=0.0_wp
  
  ! BC imin
  ! -------
  ndxv=1
  if (BC_face(1,1)%sort<=0) then
     select case (ngh_v)
     case (2)
        ndxv=3
        i=1
        !do j=1,ny
        do j=ndy_v1,nfy_v1
           do l=0,4
              x_ksi_v(i,j)=x_ksi_v(i,j)+a04(1+l)*xc(i+l,j)
              y_ksi_v(i,j)=y_ksi_v(i,j)+a04(1+l)*yc(i+l,j)
           enddo
        enddo
        i=2
        !do j=1,ny
        do j=ndy_v1,nfy_v1
           do l=-1,3
              x_ksi_v(i,j)=x_ksi_v(i,j)+a13(2+l)*xc(i+l,j)
              y_ksi_v(i,j)=y_ksi_v(i,j)+a13(2+l)*yc(i+l,j)
           enddo
        enddo
     case (1)
        ndxv=2
        i=1
        !do j=1,ny
        do j=ndy_v1,nfy_v1
           do l=0,2
              x_ksi_v(i,j)=x_ksi_v(i,j)+a02(1+l)*xc(i+l,j)
              y_ksi_v(i,j)=y_ksi_v(i,j)+a02(1+l)*yc(i+l,j)
           enddo
        enddo
        !!do j=1,ny
        !do j=ndy_v1,nfy_v1
        !   x_ksi_v(1,j)=xc(2,j)-xc(1,j)
        !   y_ksi_v(1,j)=yc(2,j)-yc(1,j)
        !enddo
     case default
        call mpistop('max ngh_v is 2 (5pts-stencil)!', 0)
     end select
  endif
  
  ! BC imax
  ! -------
  nfxv=nx
  if (BC_face(1,2)%sort<=0) then
     select case (ngh_v)
     case (2)
        nfxv=nx-2
        i=nx-1
        !do j=1,ny
        do j=ndy_v1,nfy_v1
           do l=-3,1
              x_ksi_v(i,j)=x_ksi_v(i,j)-a13(2-l)*xc(i+l,j)
              y_ksi_v(i,j)=y_ksi_v(i,j)-a13(2-l)*yc(i+l,j)
           enddo
        enddo
        i=nx
        !do j=1,ny
        do j=ndy_v1,nfy_v1
           do l=-4,0
              x_ksi_v(i,j)=x_ksi_v(i,j)-a04(1-l)*xc(i+l,j)
              y_ksi_v(i,j)=y_ksi_v(i,j)-a04(1-l)*yc(i+l,j)
           enddo
        enddo
     case (1)
        nfxv=nx-1
        i=nx
        !do j=1,ny
        do j=ndy_v1,nfy_v1
           do l=-2,0
              x_ksi_v(i,j)=x_ksi_v(i,j)-a02(1-l)*xc(i+l,j)
              y_ksi_v(i,j)=y_ksi_v(i,j)-a02(1-l)*yc(i+l,j)
           enddo
        enddo
        !!do j=1,ny
        !do j=ndy_v1,nfy_v1
        !   x_ksi_v(nx,j)=xc(nx,j)-xc(nx-1,j)
        !   y_ksi_v(nx,j)=yc(nx,j)-yc(nx-1,j)
        !enddo
     case default
        call mpistop('max ngh_v is 2 (5pts-stencil)!', 0)
     end select
  endif

  ! Interior points
  ! ---------------
  select case (ngh_v)
  case (2)
     a_=a5
  case (1)
     a_=a3        
  case default
     call mpistop('max ngh_v is 2 (5pts-stencil)!', 0)
  end select
  
  do i=ndxv,nfxv
     !do j=1,ny
     do j=ndy_v1,nfy_v1
        do l=-ngh_v,ngh_v
           x_ksi_v(i,j)= x_ksi_v(i,j) + a_(l)*xc(i+l,j)
           y_ksi_v(i,j)= y_ksi_v(i,j) + a_(l)*yc(i+l,j)
        enddo
     enddo
  enddo

end subroutine grid_metrics_ksi_v

!===============================================================================
subroutine grid_metrics_eta_v
!===============================================================================
  !> Compute curvilinear metrics along eta 
  !> - stencil -ngh_v:+ngh_v for viscous fluxes -
!===============================================================================
  use mod_coeff_deriv
  use mod_grid
  use mod_bc
  use warnstop
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j,l
  integer :: ndyv,nfyv
  real(wp), dimension(-ngh_v:ngh_v) :: a_ ! interior DF scheme
  ! ---------------------------------------------------------------------------

  ! Initializations
  ! ---------------
  allocate(x_eta_v(nx1_v:nx2_v,ny1_v:ny2_v),y_eta_v(nx1_v:nx2_v,ny1_v:ny2_v))
  x_eta_v=0.0_wp
  y_eta_v=0.0_wp
  
  ! BC jmin
  ! -------
  ndyv=1
  if (BC_face(2,1)%sort<=0) then
     select case (ngh_v)
     case (2)
        ndyv=3
        j=1
        !do i=1,nx
        do i=ndx_v1,nfx_v1
           do l=0,4
              x_eta_v(i,j)=x_eta_v(i,j)+a04(1+l)*xc(i,j+l)
              y_eta_v(i,j)=y_eta_v(i,j)+a04(1+l)*yc(i,j+l)
           enddo
        enddo
        j=2
        !do i=1,nx
        do i=ndx_v1,nfx_v1
           do l=-1,3
              x_eta_v(i,j)=x_eta_v(i,j)+a13(2+l)*xc(i,j+l)
              y_eta_v(i,j)=y_eta_v(i,j)+a13(2+l)*yc(i,j+l)
           enddo
        enddo
     case (1)
        ndyv=2
        j=1
        !do i=1,nx
        do i=ndx_v1,nfx_v1
           do l=0,2
              x_eta_v(i,j)=x_eta_v(i,j)+a02(1+l)*xc(i,j+l)
              y_eta_v(i,j)=y_eta_v(i,j)+a02(1+l)*yc(i,j+l)
           enddo
        enddo
        !!do i=1,nx
        !do i=ndx_v1,nfx_v1
        !   x_eta_v(i,1)=xc(i,2)-xc(i,1)
        !   y_eta_v(i,1)=yc(i,2)-yc(i,1)
        !enddo
     case default
        call mpistop('max ngh_v is 2 (5pts-stencil)!', 0)
     end select
  endif

  ! BC jmax
  ! -------
  nfyv=ny
  if (BC_face(2,2)%sort<=0) then
     select case (ngh_v)
     case (2)
        nfyv=ny-2
        j=ny-1
        !do i=1,nx
        do i=ndx_v1,nfx_v1
           do l=-3,1
              x_eta_v(i,j)=x_eta_v(i,j)-a13(2-l)*xc(i,j+l)
              y_eta_v(i,j)=y_eta_v(i,j)-a13(2-l)*yc(i,j+l)
           enddo
        enddo
        j=ny
        !do i=1,nx
        do i=ndx_v1,nfx_v1
           do l=-4,0
              x_eta_v(i,j)=x_eta_v(i,j)-a04(1-l)*xc(i,j+l)
              y_eta_v(i,j)=y_eta_v(i,j)-a04(1-l)*yc(i,j+l)
           enddo
        enddo
     case (1)
        nfyv=ny-1
        j=ny
        !do i=1,nx
        do i=ndx_v1,nfx_v1
           do l=-2,0
              x_eta_v(i,j)=x_eta_v(i,j)-a02(1-l)*xc(i,j+l)
              y_eta_v(i,j)=y_eta_v(i,j)-a02(1-l)*yc(i,j+l)
           enddo
        enddo
        !!do i=1,nx
        !do i=ndx_v1,nfx_v1
        !   x_eta_v(i,ny)=xc(i,ny)-xc(i,ny-1)
        !   y_eta_v(i,ny)=yc(i,ny)-yc(i,ny-1)
        !enddo
     case default
        call mpistop('max ngh_v is 2 (5pts-stencil)!', 0)
     end select
  endif

  ! Interior points
  ! ---------------
  select case (ngh_v)
  case (2)
     a_=a5
  case (1)
     a_=a3        
  case default
     call mpistop('max ngh_v is 2 (5pts-stencil)!', 0)
  end select
  
  do j=ndyv,nfyv
     !do i=1,nx
     do i=ndx_v1,nfx_v1
        do l=-ngh_v,ngh_v
           x_eta_v(i,j)= x_eta_v(i,j) + a_(l)*xc(i,j+l)
           y_eta_v(i,j)= y_eta_v(i,j) + a_(l)*yc(i,j+l)
        enddo
     enddo
  enddo

end subroutine grid_metrics_eta_v

!===============================================================================
subroutine grid_metrics_ijacob_v
!===============================================================================
  !> Compute Jacobians of metrics for viscous fluxes
  !> - stencil -ngh_v:+ngh_v for viscous fluxes -
!===============================================================================
  use mod_grid
  use mod_mpi ! for temporary check
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j
  ! ---------------------------------------------------------------------------

  allocate(ijacob_v(nx1_v:nx2_v,ny1_v:ny2_v))

  ! Jacobian of metrics transformation
  ! ----------------------------------
  ijacob_v=1.0_wp

  !do i=1,nx
  !   do j=1,ny
  do j=ndy_v1,nfy_v1
     do i=1,nx
        ijacob_v(i,j)=x_ksi_v(i,j)*y_eta_v(i,j)-x_eta_v(i,j)*y_ksi_v(i,j)
        if (ijacob_v(i,j)==0) print *,'jac1_v',iproc,i,j
     enddo
  enddo

  do j=1,ny
     do i=ndx_v1,nfx_v1
        ijacob_v(i,j)=x_ksi_v(i,j)*y_eta_v(i,j)-x_eta_v(i,j)*y_ksi_v(i,j)
        if (ijacob_v(i,j)==0) print *,'jac2_v',iproc,i,j
     enddo
  enddo

  ! Inverse of Jacobian
  ! -------------------
  ijacob_v=1.0_wp/ijacob_v

end subroutine grid_metrics_ijacob_v

!===============================================================================
subroutine grid_metrics_wall_imin
!===============================================================================
  !> Compute curvilinear metrics for wall BC at first rows of i-indices (imin)
  !> (wall BC with reduced stencil approaching the wall)
!===============================================================================
  use mod_grid
  use mod_coeff_deriv
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j
  ! ---------------------------------------------------------------------------

  ! Second-order (3-pt stencil) for point i=1
  ! =========================================
  i=1
  do j=ny1,ny2
     x_ksi(i,j)=a02(1)*xc(i,j)+a02(2)*xc(i+1,j)+a02(3)*xc(i+2,j) 
     y_ksi(i,j)=a02(1)*yc(i,j)+a02(2)*yc(i+1,j)+a02(3)*yc(i+2,j) 
  enddo

  if (ngh<2) return

  ! Second-order (3-pt stencil) for point i=2
  ! =========================================
  i=2
  do j=ny1,ny2
     x_ksi(i,j)=0.5_wp*(xc(i+1,j)-xc(i-1,j))
     y_ksi(i,j)=0.5_wp*(yc(i+1,j)-yc(i-1,j))
  enddo

  if (ngh<3) return
  
  ! Fourth-order (5-pt stencil) for point i=3
  ! =========================================
  i=3
  do j=ny1,ny2
     x_ksi(i,j)= a5(1)*(xc(i+1,j)-xc(i-1,j)) &
               + a5(2)*(xc(i+2,j)-xc(i-2,j))
     y_ksi(i,j)= a5(1)*(yc(i+1,j)-yc(i-1,j)) &
               + a5(2)*(yc(i+2,j)-yc(i-2,j))
  enddo

  if (ngh<4) return

  ! Sixth-order (7-pt stencil) for point i=4
  ! ========================================
  i=4
  do j=ny1,ny2
     x_ksi(i,j)= a7(1)*(xc(i+1,j)-xc(i-1,j)) &
               + a7(2)*(xc(i+2,j)-xc(i-2,j)) &
               + a7(3)*(xc(i+3,j)-xc(i-3,j))
     y_ksi(i,j)= a7(1)*(yc(i+1,j)-yc(i-1,j)) &
               + a7(2)*(yc(i+2,j)-yc(i-2,j)) &
               + a7(3)*(yc(i+3,j)-yc(i-3,j))
  enddo

  if (ngh<5) return
  
  ! Eighth-order (9-pt stencil) for point i=5
  ! =========================================
  i=5
  do j=ny1,ny2
     x_ksi(i,j)= a9(1)*(xc(i+1,j)-xc(i-1,j)) &
               + a9(2)*(xc(i+2,j)-xc(i-2,j)) &
               + a9(3)*(xc(i+3,j)-xc(i-3,j)) &
               + a9(4)*(xc(i+4,j)-xc(i-4,j))
     y_ksi(i,j)= a9(1)*(yc(i+1,j)-yc(i-1,j)) &
               + a9(2)*(yc(i+2,j)-yc(i-2,j)) &
               + a9(3)*(yc(i+3,j)-yc(i-3,j)) &
               + a9(4)*(yc(i+4,j)-yc(i-4,j))
  enddo

end subroutine grid_metrics_wall_imin

!===============================================================================
subroutine grid_metrics_wall_imax
!===============================================================================
  !> Compute curvilinear metrics for wall BC at last rows of i-indices (imax)
  !> (wall BC with reduced stencil approaching the wall)
!===============================================================================
  use mod_grid
  use mod_coeff_deriv
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j
  ! ---------------------------------------------------------------------------

  ! Second-order (3-pt stencil) for point i=nx
  ! ==========================================
  i=nx
  do j=ny1,ny2
     x_ksi(i,j)=a20(1)*xc(i,j)+a20(2)*xc(i-1,j)+a20(3)*xc(i-2,j)
     y_ksi(i,j)=a20(1)*yc(i,j)+a20(2)*yc(i-1,j)+a20(3)*yc(i-2,j)
  enddo

  if (ngh<2) return

  ! Second-order (3-pt stencil) for point i=nx-1
  ! ============================================
  i=nx-1
  do j=ny1,ny2
     x_ksi(i,j)=0.5_wp*(xc(i+1,j)-xc(i-1,j))
     y_ksi(i,j)=0.5_wp*(yc(i+1,j)-yc(i-1,j))
  enddo

  if (ngh<3) return

  ! Fourth-order (5-pt stencil) for point i=nx-2
  ! ============================================
  i=nx-2
  do j=ny1,ny2
     x_ksi(i,j)= a5(1)*(xc(i+1,j)-xc(i-1,j)) &
               + a5(2)*(xc(i+2,j)-xc(i-2,j))
     y_ksi(i,j)= a5(1)*(yc(i+1,j)-yc(i-1,j)) &
               + a5(2)*(yc(i+2,j)-yc(i-2,j))
  enddo

  if (ngh<4) return

  ! Sixth-order (7-pt stencil) for point i=nx-3
  ! ===========================================
  i=nx-3
  do j=ny1,ny2
     x_ksi(i,j)= a7(1)*(xc(i+1,j)-xc(i-1,j)) &
               + a7(2)*(xc(i+2,j)-xc(i-2,j)) &
               + a7(3)*(xc(i+3,j)-xc(i-3,j))
     y_ksi(i,j)= a7(1)*(yc(i+1,j)-yc(i-1,j)) &
               + a7(2)*(yc(i+2,j)-yc(i-2,j)) &
               + a7(3)*(yc(i+3,j)-yc(i-3,j))
  enddo

  if (ngh<5) return

  ! Eighth-order (9-pt stencil) for point i=nx-4
  ! ============================================
  i=nx-4
  do j=ny1,ny2
     x_ksi(i,j)= a9(1)*(xc(i+1,j)-xc(i-1,j)) &
               + a9(2)*(xc(i+2,j)-xc(i-2,j)) &
               + a9(3)*(xc(i+3,j)-xc(i-3,j)) &
               + a9(4)*(xc(i+4,j)-xc(i-4,j))
     y_ksi(i,j)= a9(1)*(yc(i+1,j)-yc(i-1,j)) &
               + a9(2)*(yc(i+2,j)-yc(i-2,j)) &
               + a9(3)*(yc(i+3,j)-yc(i-3,j)) &
               + a9(4)*(yc(i+4,j)-yc(i-4,j))
  enddo

end subroutine grid_metrics_wall_imax

!===============================================================================
subroutine grid_metrics_wall_jmin
!===============================================================================
  !> Compute curvilinear metrics for wall BC at first rows of j-indices (jmin)
  !> (wall BC with reduced stencil approaching the wall)
!===============================================================================
  use mod_grid
  use mod_coeff_deriv
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j
  ! ---------------------------------------------------------------------------

  ! Second-order (3-pt stencil) for point j=1
  ! =========================================
  j=1
  do i=nx1,nx2
     x_eta(i,j)=a02(1)*xc(i,j)+a02(2)*xc(i,j+1)+a02(3)*xc(i,j+2) 
     y_eta(i,j)=a02(1)*yc(i,j)+a02(2)*yc(i,j+1)+a02(3)*yc(i,j+2) 
  enddo

  if (ngh<2) return

  ! Second-order (3-pt stencil) for point j=2
  ! =========================================
  j=2
  do i=nx1,nx2
     x_eta(i,j)=0.5_wp*(xc(i,j+1)-xc(i,j-1))
     y_eta(i,j)=0.5_wp*(yc(i,j+1)-yc(i,j-1))
  enddo

  if (ngh<3) return

  ! Fourth-order (5-pt stencil) for point j=3
  ! =========================================
  j=3
  do i=nx1,nx2
     x_eta(i,j)= a5(1)*(xc(i,j+1)-xc(i,j-1)) &
               + a5(2)*(xc(i,j+2)-xc(i,j-2))
     y_eta(i,j)= a5(1)*(yc(i,j+1)-yc(i,j-1)) &
               + a5(2)*(yc(i,j+2)-yc(i,j-2))
  enddo
     
  if (ngh<4) return

  ! Sixth-order (7-pt stencil) for point j=4
  ! ========================================
  j=4
  do i=nx1,nx2
     x_eta(i,j)= a7(1)*(xc(i,j+1)-xc(i,j-1)) &
               + a7(2)*(xc(i,j+2)-xc(i,j-2)) &
               + a7(3)*(xc(i,j+3)-xc(i,j-3))
     y_eta(i,j)= a7(1)*(yc(i,j+1)-yc(i,j-1)) &
               + a7(2)*(yc(i,j+2)-yc(i,j-2)) &
               + a7(3)*(yc(i,j+3)-yc(i,j-3))
  enddo

  if (ngh<5) return

  ! Eighth-order (9-pt stencil) for point j=5
  ! =========================================
  j=5
  do i=nx1,nx2
     x_eta(i,j)= a9(1)*(xc(i,j+1)-xc(i,j-1)) &
               + a9(2)*(xc(i,j+2)-xc(i,j-2)) &
               + a9(3)*(xc(i,j+3)-xc(i,j-3)) &
               + a9(4)*(xc(i,j+4)-xc(i,j-4))
     y_eta(i,j)= a9(1)*(yc(i,j+1)-yc(i,j-1)) &
               + a9(2)*(yc(i,j+2)-yc(i,j-2)) &
               + a9(3)*(yc(i,j+3)-yc(i,j-3)) &
               + a9(4)*(yc(i,j+4)-yc(i,j-4))
  enddo

end subroutine grid_metrics_wall_jmin

!===============================================================================
subroutine grid_metrics_wall_jmax
!===============================================================================
  !> Compute curvilinear metrics for wall BC at last rows of j-indices (jmax)
  !> (wall BC with reduced stencil approaching the wall)
!===============================================================================
  use mod_grid
  use mod_coeff_deriv
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j
  ! ---------------------------------------------------------------------------

  ! Second-order (3-pt stencil) for point j=ny
  ! ==========================================
  j=ny
  do i=nx1,nx2
     x_eta(i,j)=a20(1)*xc(i,j)+a20(2)*xc(i,j-1)+a20(3)*xc(i,j-2)
     y_eta(i,j)=a20(1)*yc(i,j)+a20(2)*yc(i,j-1)+a20(3)*yc(i,j-2)
  enddo

  if (ngh<2) return

  ! Second-order (3-pt stencil) for point j=ny-1
  ! ============================================
  j=ny-1
  do i=nx1,nx2
     x_eta(i,j)=0.5_wp*(xc(i,j+1)-xc(i,j-1))
     y_eta(i,j)=0.5_wp*(yc(i,j+1)-yc(i,j-1))
  enddo

  if (ngh<3) return

  ! Fourth-order (5-pt stencil) for point j=ny-2
  ! ============================================
  j=ny-2
  do i=nx1,nx2
     x_eta(i,j)= a5(1)*(xc(i,j+1)-xc(i,j-1)) &
               + a5(2)*(xc(i,j+2)-xc(i,j-2))
     y_eta(i,j)= a5(1)*(yc(i,j+1)-yc(i,j-1)) &
               + a5(2)*(yc(i,j+2)-yc(i,j-2))
  enddo

  if (ngh<4) return

  ! Sixth-order (7-pt stencil) for point j=ny-3
  ! ===========================================
  j=ny-3
  do i=nx1,nx2
     x_eta(i,j)= a7(1)*(xc(i,j+1)-xc(i,j-1)) &
               + a7(2)*(xc(i,j+2)-xc(i,j-2)) &
               + a7(3)*(xc(i,j+3)-xc(i,j-3))
     y_eta(i,j)= a7(1)*(yc(i,j+1)-yc(i,j-1)) &
               + a7(2)*(yc(i,j+2)-yc(i,j-2)) &
               + a7(3)*(yc(i,j+3)-yc(i,j-3))
  enddo

  if (ngh<5) return

  ! Eighth-order (9-pt stencil) for point j=ny-4
  ! ============================================
  j=ny-4
  do i=nx1,nx2
     x_eta(i,j)= a9(1)*(xc(i,j+1)-xc(i,j-1)) &
               + a9(2)*(xc(i,j+2)-xc(i,j-2)) &
               + a9(3)*(xc(i,j+3)-xc(i,j-3)) &
               + a9(4)*(xc(i,j+4)-xc(i,j-4))
     y_eta(i,j)= a9(1)*(yc(i,j+1)-yc(i,j-1)) &
               + a9(2)*(yc(i,j+2)-yc(i,j-2)) &
               + a9(3)*(yc(i,j+3)-yc(i,j-3)) &
               + a9(4)*(yc(i,j+4)-yc(i,j-4))
  enddo

end subroutine grid_metrics_wall_jmax

!===============================================================================
subroutine grid_metrics_wall_imin_SBP4
!===============================================================================
  !> Compute curvilinear metrics for wall BC at first rows of i-indices (imin)
  !> (wall BC with SBP4: Summation by Parts order 4)
!===============================================================================
  use mod_grid
  use mod_coeff_deriv
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j
  ! ---------------------------------------------------------------------------

  ! Point #1: stencil [o x x x]
  ! ===========================
  !i=1
  do j=ny1,ny2     
     x_ksi(1,j)= as4p0(1)*xc(1,j)+as4p0(2)*xc(2,j) &
               + as4p0(3)*xc(3,j)+as4p0(4)*xc(4,j)
     y_ksi(1,j)= as4p0(1)*yc(1,j)+as4p0(2)*yc(2,j) &
               + as4p0(3)*yc(3,j)+as4p0(4)*yc(4,j)
  enddo

  ! Point #2: stencil [x o x x]
  ! ===========================
  !i=2
  do j=ny1,ny2
     x_ksi(2,j)= as4p1(1)*xc(1,j)+as4p1(2)*xc(2,j) &
               + as4p1(3)*xc(3,j)+as4p1(4)*xc(4,j)
     y_ksi(2,j)= as4p1(1)*yc(1,j)+as4p1(2)*yc(2,j) &
               + as4p1(3)*yc(3,j)+as4p1(4)*yc(4,j)
  enddo

  ! Point #3: stencil [x x o x x]
  ! =============================
  !i=3
  do j=ny1,ny2
     x_ksi(3,j)= as4p2(1)*xc(1,j)+as4p2(2)*xc(2,j) &
               + as4p2(3)*xc(3,j)+as4p2(4)*xc(4,j) &
               + as4p2(5)*xc(5,j)
     y_ksi(3,j)= as4p2(1)*yc(1,j)+as4p2(2)*yc(2,j) &
               + as4p2(3)*yc(3,j)+as4p2(4)*yc(4,j) &
               + as4p2(5)*yc(5,j)
  enddo

  ! Point #4: stencil [x x x o x x]
  ! ===============================
  !i=4
  do j=ny1,ny2
     x_ksi(4,j)= as4p3(1)*xc(1,j)+as4p3(2)*xc(2,j) &
               + as4p3(3)*xc(3,j)+as4p3(4)*xc(4,j) &
               + as4p3(5)*xc(5,j)+as4p3(6)*xc(6,j)
     y_ksi(4,j)= as4p3(1)*yc(1,j)+as4p3(2)*yc(2,j) &
               + as4p3(3)*yc(3,j)+as4p3(4)*yc(4,j) &
               + as4p3(5)*yc(5,j)+as4p3(6)*yc(6,j)
  enddo
  
  ! Eighth-order (9-pt stencil) for point i=5
  ! =========================================
  i=5
  do j=ny1,ny2
     x_ksi(i,j)= a9(1)*(xc(i+1,j)-xc(i-1,j)) &
               + a9(2)*(xc(i+2,j)-xc(i-2,j)) &
               + a9(3)*(xc(i+3,j)-xc(i-3,j)) &
               + a9(4)*(xc(i+4,j)-xc(i-4,j))
     y_ksi(i,j)= a9(1)*(yc(i+1,j)-yc(i-1,j)) &
               + a9(2)*(yc(i+2,j)-yc(i-2,j)) &
               + a9(3)*(yc(i+3,j)-yc(i-3,j)) &
               + a9(4)*(yc(i+4,j)-yc(i-4,j))
  enddo

end subroutine grid_metrics_wall_imin_SBP4

!===============================================================================
subroutine grid_metrics_wall_imax_SBP4
!===============================================================================
  !> Compute curvilinear metrics for wall BC at last rows of i-indices (imax)
  !> (wall BC with SBP4: Summation by Parts order 4)
!===============================================================================
  use mod_grid
  use mod_coeff_deriv
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j
  ! ---------------------------------------------------------------------------

  ! Point #1: stencil [o x x x]
  ! ===========================
  i=nx
  do j=ny1,ny2
     x_ksi(i,j)= as4m0(1)*xc(nx  ,j)+as4m0(2)*xc(nx-1,j) &
               + as4m0(3)*xc(nx-2,j)+as4m0(4)*xc(nx-3,j)
     y_ksi(i,j)= as4m0(1)*yc(nx  ,j)+as4m0(2)*yc(nx-1,j) &
               + as4m0(3)*yc(nx-2,j)+as4m0(4)*yc(nx-3,j)
  enddo

  ! Point #2: stencil [x o x x]
  ! ===========================
  i=nx-1
  do j=ny1,ny2
     x_ksi(i,j)= as4m1(1)*xc(nx  ,j)+as4m1(2)*xc(nx-1,j) &
               + as4m1(3)*xc(nx-2,j)+as4m1(4)*xc(nx-3,j)
     y_ksi(i,j)= as4m1(1)*yc(nx  ,j)+as4m1(2)*yc(nx-1,j) &
               + as4m1(3)*yc(nx-2,j)+as4m1(4)*yc(nx-3,j)
  enddo

  ! Point #3: stencil [x x o x x]
  ! =============================
  i=nx-2
  do j=ny1,ny2
     x_ksi(i,j)= as4m2(1)*xc(nx  ,j)+as4m2(2)*xc(nx-1,j) &
               + as4m2(3)*xc(nx-2,j)+as4m2(4)*xc(nx-3,j) &
               + as4m2(5)*xc(nx-4,j)
     y_ksi(i,j)= as4m2(1)*yc(nx  ,j)+as4m2(2)*yc(nx-1,j) &
               + as4m2(3)*yc(nx-2,j)+as4m2(4)*yc(nx-3,j) &
               + as4m2(5)*yc(nx-4,j)
  enddo

  ! Point #4: stencil [x x x o x x]
  ! ===============================
  i=nx-3
  do j=ny1,ny2
     x_ksi(i,j)= as4m3(1)*xc(nx  ,j)+as4m3(2)*xc(nx-1,j) &
               + as4m3(3)*xc(nx-2,j)+as4m3(4)*xc(nx-3,j) &
               + as4m3(5)*xc(nx-4,j)+as4m3(6)*xc(nx-5,j)
     y_ksi(i,j)= as4m3(1)*yc(nx  ,j)+as4m3(2)*yc(nx-1,j) &
               + as4m3(3)*yc(nx-2,j)+as4m3(4)*yc(nx-3,j) &
               + as4m3(5)*yc(nx-4,j)+as4m3(6)*yc(nx-5,j)
  enddo

  ! Eighth-order (9-pt stencil) for point i=nx-4
  ! ============================================
  i=nx-4
  do j=ny1,ny2
     x_ksi(i,j)= a9(1)*(xc(i+1,j)-xc(i-1,j)) &
               + a9(2)*(xc(i+2,j)-xc(i-2,j)) &
               + a9(3)*(xc(i+3,j)-xc(i-3,j)) &
               + a9(4)*(xc(i+4,j)-xc(i-4,j))
     y_ksi(i,j)= a9(1)*(yc(i+1,j)-yc(i-1,j)) &
               + a9(2)*(yc(i+2,j)-yc(i-2,j)) &
               + a9(3)*(yc(i+3,j)-yc(i-3,j)) &
               + a9(4)*(yc(i+4,j)-yc(i-4,j))
  enddo

end subroutine grid_metrics_wall_imax_SBP4

!===============================================================================
subroutine grid_metrics_wall_jmin_SBP4
!===============================================================================
  !> Compute curvilinear metrics for wall BC at first rows of j-indices (jmin)
  !> (wall BC with SBP4: Summation by Parts order 4)
!===============================================================================
  use mod_grid
  use mod_coeff_deriv
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j
  ! ---------------------------------------------------------------------------

  ! Point #1: stencil [o x x x]
  ! ===========================
  !j=1
  do i=nx1,nx2
     x_eta(i,1)= as4p0(1)*xc(i,1)+as4p0(2)*xc(i,2) &
               + as4p0(3)*xc(i,3)+as4p0(4)*xc(i,4)
     y_eta(i,1)= as4p0(1)*yc(i,1)+as4p0(2)*yc(i,2) &
               + as4p0(3)*yc(i,3)+as4p0(4)*yc(i,4)
  enddo

  ! Point #2: stencil [x o x x]
  ! ===========================
  !j=2
  do i=nx1,nx2
     x_eta(i,2)= as4p1(1)*xc(i,1)+as4p1(2)*xc(i,2) &
               + as4p1(3)*xc(i,3)+as4p1(4)*xc(i,4)
     y_eta(i,2)= as4p1(1)*yc(i,1)+as4p1(2)*yc(i,2) &
               + as4p1(3)*yc(i,3)+as4p1(4)*yc(i,4)
  enddo

  ! Point #3: stencil [x x o x x]
  ! =============================
  !j=3
  do i=nx1,nx2
     x_eta(i,3)= as4p2(1)*xc(i,1)+as4p2(2)*xc(i,2) &
               + as4p2(3)*xc(i,3)+as4p2(4)*xc(i,4) &
               + as4p2(5)*xc(i,5)
     y_eta(i,3)= as4p2(1)*yc(i,1)+as4p2(2)*yc(i,2) &
               + as4p2(3)*yc(i,3)+as4p2(4)*yc(i,4) &
               + as4p2(5)*yc(i,5)
  enddo
     
  ! Point #4: stencil [x x x o x x]
  ! ===============================
  !j=4
  do i=nx1,nx2
     x_eta(i,4)= as4p3(1)*xc(i,1)+as4p3(2)*xc(i,2) &
               + as4p3(3)*xc(i,3)+as4p3(4)*xc(i,4) &
               + as4p3(5)*xc(i,5)+as4p3(6)*xc(i,6)
     y_eta(i,4)= as4p3(1)*yc(i,1)+as4p3(2)*yc(i,2) &
               + as4p3(3)*yc(i,3)+as4p3(4)*yc(i,4) &
               + as4p3(5)*yc(i,5)+as4p3(6)*yc(i,6)
  enddo

  ! Eighth-order (9-pt stencil) for point j=5
  ! =========================================
  j=5
  do i=nx1,nx2
     x_eta(i,j)= a9(1)*(xc(i,j+1)-xc(i,j-1)) &
               + a9(2)*(xc(i,j+2)-xc(i,j-2)) &
               + a9(3)*(xc(i,j+3)-xc(i,j-3)) &
               + a9(4)*(xc(i,j+4)-xc(i,j-4))
     y_eta(i,j)= a9(1)*(yc(i,j+1)-yc(i,j-1)) &
               + a9(2)*(yc(i,j+2)-yc(i,j-2)) &
               + a9(3)*(yc(i,j+3)-yc(i,j-3)) &
               + a9(4)*(yc(i,j+4)-yc(i,j-4))
  enddo

end subroutine grid_metrics_wall_jmin_SBP4

!===============================================================================
subroutine grid_metrics_wall_jmax_SBP4
!===============================================================================
  !> Compute curvilinear metrics for wall BC at last rows of j-indices (jmax)
  !> (wall BC with SBP4: Summation by Parts order 4)
!===============================================================================
  use mod_grid
  use mod_coeff_deriv
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j
  ! ---------------------------------------------------------------------------

  ! Point #1: stencil [o x x x]
  ! ===========================
  j=ny
  do i=nx1,nx2
     x_eta(i,j)= as4m0(1)*xc(i,ny  )+as4m0(2)*xc(i,ny-1) &
               + as4m0(3)*xc(i,ny-2)+as4m0(4)*xc(i,ny-3)
     y_eta(i,j)= as4m0(1)*yc(i,ny  )+as4m0(2)*yc(i,ny-1) &
               + as4m0(3)*yc(i,ny-2)+as4m0(4)*yc(i,ny-3)
  enddo

  ! Point #2: stencil [x o x x]
  ! ===========================
  j=ny-1
  do i=nx1,nx2
     x_eta(i,j)= as4m1(1)*xc(i,ny  )+as4m1(2)*xc(i,ny-1) &
               + as4m1(3)*xc(i,ny-2)+as4m1(4)*xc(i,ny-3)
     y_eta(i,j)= as4m1(1)*yc(i,ny  )+as4m1(2)*yc(i,ny-1) &
               + as4m1(3)*yc(i,ny-2)+as4m1(4)*yc(i,ny-3)
  enddo

  ! Point #3: stencil [x x o x x]
  ! =============================
  j=ny-2
  do i=nx1,nx2
     x_eta(i,j)= as4m2(1)*xc(i,ny  )+as4m2(2)*xc(i,ny-1) &
               + as4m2(3)*xc(i,ny-2)+as4m2(4)*xc(i,ny-3) &
               + as4m2(5)*xc(i,ny-4)
     y_eta(i,j)= as4m2(1)*yc(i,ny  )+as4m2(2)*yc(i,ny-1) &
               + as4m2(3)*yc(i,ny-2)+as4m2(4)*yc(i,ny-3) &
               + as4m2(5)*yc(i,ny-4)
  enddo

  ! Point #4: stencil [x x x o x x]
  ! ===============================
  j=ny-3
  do i=nx1,nx2
     x_eta(i,j)= as4m3(1)*xc(i,ny  )+as4m3(2)*xc(i,ny-1) &
               + as4m3(3)*xc(i,ny-2)+as4m3(4)*xc(i,ny-3) &
               + as4m3(5)*xc(i,ny-4)+as4m3(6)*xc(i,ny-5)
     y_eta(i,j)= as4m3(1)*yc(i,ny  )+as4m3(2)*yc(i,ny-1) &
               + as4m3(3)*yc(i,ny-2)+as4m3(4)*yc(i,ny-3) &
               + as4m3(5)*yc(i,ny-4)+as4m3(6)*yc(i,ny-5)
  enddo

  ! Eighth-order (9-pt stencil) for point j=ny-4
  ! ============================================
  j=ny-4
  do i=nx1,nx2
     x_eta(i,j)= a9(1)*(xc(i,j+1)-xc(i,j-1)) &
               + a9(2)*(xc(i,j+2)-xc(i,j-2)) &
               + a9(3)*(xc(i,j+3)-xc(i,j-3)) &
               + a9(4)*(xc(i,j+4)-xc(i,j-4))
     y_eta(i,j)= a9(1)*(yc(i,j+1)-yc(i,j-1)) &
               + a9(2)*(yc(i,j+2)-yc(i,j-2)) &
               + a9(3)*(yc(i,j+3)-yc(i,j-3)) &
               + a9(4)*(yc(i,j+4)-yc(i,j-4))
  enddo

end subroutine grid_metrics_wall_jmax_SBP4

!===============================================================================
subroutine grid_metrics_wall_imin_jmin
!===============================================================================
  !> Compute curvilinear metrics for wall edge imin-jmin
  !> (wall BC with reduced stencil approaching the wall)
!===============================================================================
  use mod_grid
  use mod_coeff_deriv
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j
  ! ---------------------------------------------------------------------------

  allocate(x_eta_imin_jmin(nx1:2*ngh,1:ngh),y_eta_imin_jmin(nx1:2*ngh,1:ngh))
  allocate(x_ksi_imin_jmin(1:ngh,ny1:2*ngh),y_ksi_imin_jmin(1:ngh,ny1:2*ngh))

  ! First-order (3-pt stencil) for point j=1
  ! ========================================
  j=1
  do i=nx1,2*ngh
     x_eta_imin_jmin(i,j)=a02(1)*xc(i,j)+a02(2)*xc(i,j+1)+a02(3)*xc(i,j+2) 
     y_eta_imin_jmin(i,j)=a02(1)*yc(i,j)+a02(2)*yc(i,j+1)+a02(3)*yc(i,j+2) 
  enddo
  !do i=nx1,2*ngh
  !   x_eta_imin_jmin(i,j)=xc(i,j+1)-xc(i,j)
  !   y_eta_imin_jmin(i,j)=yc(i,j+1)-yc(i,j)
  !enddo

  ! First-order (3-pt stencil) for point i=1
  ! ========================================
  i=1
  do j=ny1,2*ngh
     x_ksi_imin_jmin(i,j)=a02(1)*xc(i,j)+a02(2)*xc(i+1,j)+a02(3)*xc(i+2,j) 
     y_ksi_imin_jmin(i,j)=a02(1)*yc(i,j)+a02(2)*yc(i+1,j)+a02(3)*yc(i+2,j) 
  enddo
  !do j=ny1,2*ngh
  !   x_ksi_imin_jmin(i,j)=xc(i+1,j)-xc(i,j)
  !   y_ksi_imin_jmin(i,j)=yc(i+1,j)-yc(i,j)
  !enddo

  if (ngh<2) return

  ! Second-order (3-pt stencil) for point j=2
  ! =========================================
  j=2
  do i=nx1,2*ngh
     x_eta_imin_jmin(i,j)=0.5_wp*(xc(i,j+1)-xc(i,j-1))
     y_eta_imin_jmin(i,j)=0.5_wp*(yc(i,j+1)-yc(i,j-1))
  enddo
  
  ! Second-order (3-pt stencil) for point i=2
  ! =========================================
  i=2
  do j=ny1,2*ngh
     x_ksi_imin_jmin(i,j)=0.5_wp*(xc(i+1,j)-xc(i-1,j))
     y_ksi_imin_jmin(i,j)=0.5_wp*(yc(i+1,j)-yc(i-1,j))
  enddo

  if (ngh<3) return

  ! Fourth-order (5-pt stencil) for point j=3
  ! =========================================
  j=3
  do i=nx1,2*ngh
     x_eta_imin_jmin(i,j)= a5(1)*(xc(i,j+1)-xc(i,j-1)) &
                         + a5(2)*(xc(i,j+2)-xc(i,j-2))
     y_eta_imin_jmin(i,j)= a5(1)*(yc(i,j+1)-yc(i,j-1)) &
                         + a5(2)*(yc(i,j+2)-yc(i,j-2))
  enddo
     
  ! Fourth-order (5-pt stencil) for point i=3
  ! =========================================
  i=3
  do j=ny1,2*ngh
     x_ksi_imin_jmin(i,j)= a5(1)*(xc(i+1,j)-xc(i-1,j)) &
                         + a5(2)*(xc(i+2,j)-xc(i-2,j))
     y_ksi_imin_jmin(i,j)= a5(1)*(yc(i+1,j)-yc(i-1,j)) &
                         + a5(2)*(yc(i+2,j)-yc(i-2,j))
  enddo
     
  if (ngh<4) return

  ! Sixth-order (7-pt stencil) for point j=4
  ! ========================================
  j=4
  do i=nx1,2*ngh
     x_eta_imin_jmin(i,j)= a7(1)*(xc(i,j+1)-xc(i,j-1)) &
                         + a7(2)*(xc(i,j+2)-xc(i,j-2)) &
                         + a7(3)*(xc(i,j+3)-xc(i,j-3))
     y_eta_imin_jmin(i,j)= a7(1)*(yc(i,j+1)-yc(i,j-1)) &
                         + a7(2)*(yc(i,j+2)-yc(i,j-2)) &
                         + a7(3)*(yc(i,j+3)-yc(i,j-3))
  enddo

  ! Sixth-order (7-pt stencil) for point i=4
  ! ========================================
  i=4
  do j=ny1,2*ngh
     x_ksi_imin_jmin(i,j)= a7(1)*(xc(i+1,j)-xc(i-1,j)) &
                         + a7(2)*(xc(i+2,j)-xc(i-2,j)) &
                         + a7(3)*(xc(i+3,j)-xc(i-3,j))
     y_ksi_imin_jmin(i,j)= a7(1)*(yc(i+1,j)-yc(i-1,j)) &
                         + a7(2)*(yc(i+2,j)-yc(i-2,j)) &
                         + a7(3)*(yc(i+3,j)-yc(i-3,j))
  enddo

  if (ngh<5) return

  ! Eighth-order (9-pt stencil) for point j=5
  ! =========================================
  j=5
  do i=nx1,2*ngh
     x_eta_imin_jmin(i,j)= a9(1)*(xc(i,j+1)-xc(i,j-1)) &
                         + a9(2)*(xc(i,j+2)-xc(i,j-2)) &
                         + a9(3)*(xc(i,j+3)-xc(i,j-3)) &
                         + a9(4)*(xc(i,j+4)-xc(i,j-4))
     y_eta_imin_jmin(i,j)= a9(1)*(yc(i,j+1)-yc(i,j-1)) &
                         + a9(2)*(yc(i,j+2)-yc(i,j-2)) &
                         + a9(3)*(yc(i,j+3)-yc(i,j-3)) &
                         + a9(4)*(yc(i,j+4)-yc(i,j-4))
  enddo

  ! Eighth-order (9-pt stencil) for point i=5
  ! =========================================
  i=5
  do j=ny1,2*ngh
     x_ksi_imin_jmin(i,j)= a9(1)*(xc(i+1,j)-xc(i-1,j)) &
                         + a9(2)*(xc(i+2,j)-xc(i-2,j)) &
                         + a9(3)*(xc(i+3,j)-xc(i-3,j)) &
                         + a9(4)*(xc(i+4,j)-xc(i-4,j))
     y_ksi_imin_jmin(i,j)= a9(1)*(yc(i+1,j)-yc(i-1,j)) &
                         + a9(2)*(yc(i+2,j)-yc(i-2,j)) &
                         + a9(3)*(yc(i+3,j)-yc(i-3,j)) &
                         + a9(4)*(yc(i+4,j)-yc(i-4,j))
  enddo

end subroutine grid_metrics_wall_imin_jmin

!===============================================================================
subroutine grid_metrics_wall_imin_jmax
!===============================================================================
  !> Compute curvilinear metrics for wall edge imin-jmax
  !> (wall BC with reduced stencil approaching the wall)
!===============================================================================
  use mod_grid
  use mod_coeff_deriv
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j
  ! ---------------------------------------------------------------------------

  allocate(x_eta_imin_jmax(nx1:2*ngh,ny-ngh+1:ny),y_eta_imin_jmax(nx1:2*ngh,ny-ngh+1:ny))
  allocate(x_ksi_imin_jmax(1:ngh,ny-2*ngh+1:ny2),y_ksi_imin_jmax(1:ngh,ny-2*ngh+1:ny2))
  
  ! First-order (3-pt stencil) for point j=ny
  ! =========================================
  j=ny
  do i=nx1,2*ngh
     x_eta_imin_jmax(i,j)=a20(1)*xc(i,j)+a20(2)*xc(i,j-1)+a20(3)*xc(i,j-2)
     y_eta_imin_jmax(i,j)=a20(1)*yc(i,j)+a20(2)*yc(i,j-1)+a20(3)*yc(i,j-2)
  enddo
  !do i=nx1,2*ngh
  !   x_eta_imin_jmax(i,j)=xc(i,j)-xc(i,j-1)
  !   y_eta_imin_jmax(i,j)=yc(i,j)-yc(i,j-1)
  !enddo

  ! First-order (3-pt stencil) for point i=1
  ! ========================================
  i=1
  do j=ny-2*ngh+1,ny2
     x_ksi_imin_jmax(i,j)=a02(1)*xc(i,j)+a02(2)*xc(i+1,j)+a02(3)*xc(i+2,j) 
     y_ksi_imin_jmax(i,j)=a02(1)*yc(i,j)+a02(2)*yc(i+1,j)+a02(3)*yc(i+2,j) 
  enddo
  !do j=ny-2*ngh+1,ny2
  !   x_ksi_imin_jmax(i,j)=xc(i+1,j)-xc(i,j)
  !   y_ksi_imin_jmax(i,j)=yc(i+1,j)-yc(i,j)
  !enddo

  if (ngh<2) return

  ! Second-order (3-pt stencil) for point j=ny-1
  ! ============================================
  j=ny-1
  do i=nx1,2*ngh
     x_eta_imin_jmax(i,j)=0.5_wp*(xc(i,j+1)-xc(i,j-1))
     y_eta_imin_jmax(i,j)=0.5_wp*(yc(i,j+1)-yc(i,j-1))
  enddo

  ! Second-order (3-pt stencil) for point i=2
  ! =========================================
  i=2
  do j=ny-2*ngh+1,ny2
     x_ksi_imin_jmax(i,j)=0.5_wp*(xc(i+1,j)-xc(i-1,j))
     y_ksi_imin_jmax(i,j)=0.5_wp*(yc(i+1,j)-yc(i-1,j))
  enddo

  if (ngh<3) return

  ! Fourth-order (5-pt stencil) for point j=ny-2
  ! ============================================
  j=ny-2
  do i=nx1,2*ngh
     x_eta_imin_jmax(i,j)= a5(1)*(xc(i,j+1)-xc(i,j-1)) &
                         + a5(2)*(xc(i,j+2)-xc(i,j-2))
     y_eta_imin_jmax(i,j)= a5(1)*(yc(i,j+1)-yc(i,j-1)) &
                         + a5(2)*(yc(i,j+2)-yc(i,j-2))
  enddo
     
  ! Fourth-order (5-pt stencil) for point i=3
  ! =========================================
  i=3
  do j=ny-2*ngh+1,ny2
     x_ksi_imin_jmax(i,j)= a5(1)*(xc(i+1,j)-xc(i-1,j)) &
                         + a5(2)*(xc(i+2,j)-xc(i-2,j))
     y_ksi_imin_jmax(i,j)= a5(1)*(yc(i+1,j)-yc(i-1,j)) &
                         + a5(2)*(yc(i+2,j)-yc(i-2,j))
  enddo
     
  if (ngh<4) return

  ! Sixth-order (7-pt stencil) for point j=ny-3
  ! ===========================================
  j=ny-3
  do i=nx1,2*ngh
     x_eta_imin_jmax(i,j)= a7(1)*(xc(i,j+1)-xc(i,j-1)) &
                         + a7(2)*(xc(i,j+2)-xc(i,j-2)) &
                         + a7(3)*(xc(i,j+3)-xc(i,j-3))
     y_eta_imin_jmax(i,j)= a7(1)*(yc(i,j+1)-yc(i,j-1)) &
                         + a7(2)*(yc(i,j+2)-yc(i,j-2)) &
                         + a7(3)*(yc(i,j+3)-yc(i,j-3))
  enddo

  ! Sixth-order (7-pt stencil) for point i=4
  ! ========================================
  i=4
  do j=ny-2*ngh+1,ny2
     x_ksi_imin_jmax(i,j)= a7(1)*(xc(i+1,j)-xc(i-1,j)) &
                         + a7(2)*(xc(i+2,j)-xc(i-2,j)) &
                         + a7(3)*(xc(i+3,j)-xc(i-3,j))
     y_ksi_imin_jmax(i,j)= a7(1)*(yc(i+1,j)-yc(i-1,j)) &
                         + a7(2)*(yc(i+2,j)-yc(i-2,j)) &
                         + a7(3)*(yc(i+3,j)-yc(i-3,j))
  enddo

  if (ngh<5) return

  ! Eighth-order (9-pt stencil) for point j=ny-4
  ! ============================================
  j=ny-4
  do i=nx1,2*ngh
     x_eta_imin_jmax(i,j)= a9(1)*(xc(i,j+1)-xc(i,j-1)) &
                         + a9(2)*(xc(i,j+2)-xc(i,j-2)) &
                         + a9(3)*(xc(i,j+3)-xc(i,j-3)) &
                         + a9(4)*(xc(i,j+4)-xc(i,j-4))
     y_eta_imin_jmax(i,j)= a9(1)*(yc(i,j+1)-yc(i,j-1)) &
                         + a9(2)*(yc(i,j+2)-yc(i,j-2)) &
                         + a9(3)*(yc(i,j+3)-yc(i,j-3)) &
                         + a9(4)*(yc(i,j+4)-yc(i,j-4))
  enddo
  
  ! Eighth-order (9-pt stencil) for point i=5
  ! =========================================
  i=5
  do j=ny-2*ngh+1,ny2
     x_ksi_imin_jmax(i,j)= a9(1)*(xc(i+1,j)-xc(i-1,j)) &
                         + a9(2)*(xc(i+2,j)-xc(i-2,j)) &
                         + a9(3)*(xc(i+3,j)-xc(i-3,j)) &
                         + a9(4)*(xc(i+4,j)-xc(i-4,j))
     y_ksi_imin_jmax(i,j)= a9(1)*(yc(i+1,j)-yc(i-1,j)) &
                         + a9(2)*(yc(i+2,j)-yc(i-2,j)) &
                         + a9(3)*(yc(i+3,j)-yc(i-3,j)) &
                         + a9(4)*(yc(i+4,j)-yc(i-4,j))
  enddo

end subroutine grid_metrics_wall_imin_jmax

!===============================================================================
subroutine grid_metrics_wall_imax_jmin
!===============================================================================
  !> Compute curvilinear metrics for wall  edge imax-jmin
  !> (wall BC with reduced stencil approaching the wall)
!===============================================================================
  use mod_grid
  use mod_coeff_deriv
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j
  ! ---------------------------------------------------------------------------

  allocate(x_eta_imax_jmin(nx-2*ngh+1:nx2,1:ngh),y_eta_imax_jmin(nx-2*ngh+1:nx2,1:ngh))
  allocate(x_ksi_imax_jmin(nx-ngh+1:nx,nx1:2*ngh),y_ksi_imax_jmin(nx-ngh+1:nx,nx1:2*ngh))
 
  ! First-order (3-pt stencil) for point j=1
  ! ========================================
  j=1
  do i=nx-2*ngh+1,nx2
     x_eta_imax_jmin(i,j)=a02(1)*xc(i,j)+a02(2)*xc(i,j+1)+a02(3)*xc(i,j+2) 
     y_eta_imax_jmin(i,j)=a02(1)*yc(i,j)+a02(2)*yc(i,j+1)+a02(3)*yc(i,j+2) 
  enddo
  !do i=nx-2*ngh+1,nx2
  !   x_eta_imax_jmin(i,j)=xc(i,j+1)-xc(i,j)
  !   y_eta_imax_jmin(i,j)=yc(i,j+1)-yc(i,j)
  !enddo

  ! First-order (3-pt stencil) for point i=nx
  ! =========================================
  i=nx
  do j=ny1,2*ngh
     x_ksi_imax_jmin(i,j)=a20(1)*xc(i,j)+a20(2)*xc(i-1,j)+a20(3)*xc(i-2,j)
     y_ksi_imax_jmin(i,j)=a20(1)*yc(i,j)+a20(2)*yc(i-1,j)+a20(3)*yc(i-2,j)
  enddo
  !do j=ny1,2*ngh
  !   x_ksi_imax_jmin(i,j)=xc(i,j)-xc(i-1,j)
  !   y_ksi_imax_jmin(i,j)=yc(i,j)-yc(i-1,j)
  !enddo

  if (ngh<2) return

  ! Second-order (3-pt stencil) for point j=2
  ! =========================================
  j=2
  do i=nx-2*ngh+1,nx2
     x_eta_imax_jmin(i,j)=0.5_wp*(xc(i,j+1)-xc(i,j-1))
     y_eta_imax_jmin(i,j)=0.5_wp*(yc(i,j+1)-yc(i,j-1))
  enddo

  ! Second-order (3-pt stencil) for point i=nx-1
  ! ============================================
  i=nx-1
  do j=ny1,2*ngh
     x_ksi_imax_jmin(i,j)=0.5_wp*(xc(i+1,j)-xc(i-1,j))
     y_ksi_imax_jmin(i,j)=0.5_wp*(yc(i+1,j)-yc(i-1,j))
  enddo

  if (ngh<3) return

  ! Fourth-order (5-pt stencil) for point j=3
  ! =========================================
  j=3
  do i=nx-2*ngh+1,nx2
     x_eta_imax_jmin(i,j)= a5(1)*(xc(i,j+1)-xc(i,j-1)) &
                         + a5(2)*(xc(i,j+2)-xc(i,j-2))
     y_eta_imax_jmin(i,j)= a5(1)*(yc(i,j+1)-yc(i,j-1)) &
                         + a5(2)*(yc(i,j+2)-yc(i,j-2))
  enddo
     
  ! Fourth-order (5-pt stencil) for point i=nx-2
  ! ============================================
  i=nx-2
  do j=ny1,2*ngh
     x_ksi_imax_jmin(i,j)= a5(1)*(xc(i+1,j)-xc(i-1,j)) &
                         + a5(2)*(xc(i+2,j)-xc(i-2,j))
     y_ksi_imax_jmin(i,j)= a5(1)*(yc(i+1,j)-yc(i-1,j)) &
                         + a5(2)*(yc(i+2,j)-yc(i-2,j))
  enddo
     
  if (ngh<4) return

  ! Sixth-order (7-pt stencil) for point j=4
  ! ========================================
  j=4
  do i=nx-2*ngh+1,nx2
     x_eta_imax_jmin(i,j)= a7(1)*(xc(i,j+1)-xc(i,j-1)) &
                         + a7(2)*(xc(i,j+2)-xc(i,j-2)) &
                         + a7(3)*(xc(i,j+3)-xc(i,j-3))
     y_eta_imax_jmin(i,j)= a7(1)*(yc(i,j+1)-yc(i,j-1)) &
                         + a7(2)*(yc(i,j+2)-yc(i,j-2)) &
                         + a7(3)*(yc(i,j+3)-yc(i,j-3))
  enddo

  ! Sixth-order (7-pt stencil) for point i=nx-3
  ! ===========================================
  i=nx-3
  do j=ny1,2*ngh
     x_ksi_imax_jmin(i,j)= a7(1)*(xc(i+1,j)-xc(i-1,j)) &
                         + a7(2)*(xc(i+2,j)-xc(i-2,j)) &
                         + a7(3)*(xc(i+3,j)-xc(i-3,j))
     y_ksi_imax_jmin(i,j)= a7(1)*(yc(i+1,j)-yc(i-1,j)) &
                         + a7(2)*(yc(i+2,j)-yc(i-2,j)) &
                         + a7(3)*(yc(i+3,j)-yc(i-3,j))
  enddo

  if (ngh<5) return

  ! Eighth-order (9-pt stencil) for point j=5
  ! =========================================
  j=5
  do i=nx-2*ngh+1,nx2
     x_eta_imax_jmin(i,j)= a9(1)*(xc(i,j+1)-xc(i,j-1)) &
                         + a9(2)*(xc(i,j+2)-xc(i,j-2)) &
                         + a9(3)*(xc(i,j+3)-xc(i,j-3)) &
                         + a9(4)*(xc(i,j+4)-xc(i,j-4))
     y_eta_imax_jmin(i,j)= a9(1)*(yc(i,j+1)-yc(i,j-1)) &
                         + a9(2)*(yc(i,j+2)-yc(i,j-2)) &
                         + a9(3)*(yc(i,j+3)-yc(i,j-3)) &
                         + a9(4)*(yc(i,j+4)-yc(i,j-4))
  enddo

  ! Eighth-order (9-pt stencil) for point i=nx-4
  ! ============================================
  i=nx-4
  do j=ny1,2*ngh
     x_ksi_imax_jmin(i,j)= a9(1)*(xc(i+1,j)-xc(i-1,j)) &
                         + a9(2)*(xc(i+2,j)-xc(i-2,j)) &
                         + a9(3)*(xc(i+3,j)-xc(i-3,j)) &
                         + a9(4)*(xc(i+4,j)-xc(i-4,j))
     y_ksi_imax_jmin(i,j)= a9(1)*(yc(i+1,j)-yc(i-1,j)) &
                         + a9(2)*(yc(i+2,j)-yc(i-2,j)) &
                         + a9(3)*(yc(i+3,j)-yc(i-3,j)) &
                         + a9(4)*(yc(i+4,j)-yc(i-4,j))
  enddo

end subroutine grid_metrics_wall_imax_jmin

!===============================================================================
subroutine grid_metrics_wall_imax_jmax
!===============================================================================
  !> Compute curvilinear metrics for wall edge imax-jmax
  !> (wall BC with reduced stencil approaching the wall)
!===============================================================================
  use mod_grid
  use mod_coeff_deriv
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j
  ! ---------------------------------------------------------------------------

  allocate(x_eta_imax_jmax(nx-2*ngh+1:nx2,ny-ngh+1:ny),y_eta_imax_jmax(nx-2*ngh+1:nx2,ny-ngh+1:ny))
  allocate(x_ksi_imax_jmax(nx-ngh+1:nx,ny-2*ngh+1:ny2),y_ksi_imax_jmax(nx-ngh+1:nx,ny-2*ngh+1:ny2))
 
  ! First-order (3-pt stencil) for point j=ny
  ! =========================================
  j=ny
  do i=nx-2*ngh+1,nx2
     x_eta_imax_jmax(i,j)=a20(1)*xc(i,j)+a20(2)*xc(i,j-1)+a20(3)*xc(i,j-2)
     y_eta_imax_jmax(i,j)=a20(1)*yc(i,j)+a20(2)*yc(i,j-1)+a20(3)*yc(i,j-2)
  enddo
  !do i=nx-2*ngh+1,nx2
  !   x_eta_imax_jmax(i,j)=xc(i,j)-xc(i,j-1)
  !   y_eta_imax_jmax(i,j)=yc(i,j)-yc(i,j-1)
  !enddo

  ! First-order (3-pt stencil) for point i=nx
  ! =========================================
  i=nx
  do j=ny-2*ngh+1,ny2
     x_ksi_imax_jmax(i,j)=a20(1)*xc(i,j)+a20(2)*xc(i-1,j)+a20(3)*xc(i-2,j)
     y_ksi_imax_jmax(i,j)=a20(1)*yc(i,j)+a20(2)*yc(i-1,j)+a20(3)*yc(i-2,j)
  enddo
  !do j=ny-2*ngh+1,ny2
  !   x_ksi_imax_jmax(i,j)=xc(i,j)-xc(i-1,j)
  !   y_ksi_imax_jmax(i,j)=yc(i,j)-yc(i-1,j)
  !enddo

  if (ngh<2) return

  ! Second-order (3-pt stencil) for point j=ny-1
  ! ============================================
  j=ny-1
  do i=nx-2*ngh+1,nx2
     x_eta_imax_jmax(i,j)=0.5_wp*(xc(i,j+1)-xc(i,j-1))
     y_eta_imax_jmax(i,j)=0.5_wp*(yc(i,j+1)-yc(i,j-1))
  enddo

  ! Second-order (3-pt stencil) for point i=nx-1
  ! ============================================
  i=nx-1
  do j=ny-2*ngh+1,ny2
     x_ksi_imax_jmax(i,j)=0.5_wp*(xc(i+1,j)-xc(i-1,j))
     y_ksi_imax_jmax(i,j)=0.5_wp*(yc(i+1,j)-yc(i-1,j))
  enddo

  if (ngh<3) return

  ! Fourth-order (5-pt stencil) for point j=ny-2
  ! ============================================
  j=ny-2
  do i=nx-2*ngh+1,nx2
     x_eta_imax_jmax(i,j)= a5(1)*(xc(i,j+1)-xc(i,j-1)) &
                         + a5(2)*(xc(i,j+2)-xc(i,j-2))
     y_eta_imax_jmax(i,j)= a5(1)*(yc(i,j+1)-yc(i,j-1)) &
                         + a5(2)*(yc(i,j+2)-yc(i,j-2))
  enddo
     
  ! Fourth-order (5-pt stencil) for point i=nx-2
  ! ============================================
  i=nx-2
  do j=ny-2*ngh+1,ny2
     x_ksi_imax_jmax(i,j)= a5(1)*(xc(i+1,j)-xc(i-1,j)) &
                         + a5(2)*(xc(i+2,j)-xc(i-2,j))
     y_ksi_imax_jmax(i,j)= a5(1)*(yc(i+1,j)-yc(i-1,j)) &
                         + a5(2)*(yc(i+2,j)-yc(i-2,j))
  enddo
     
  if (ngh<4) return

  ! Sixth-order (7-pt stencil) for point j=ny-3
  ! ===========================================
  j=ny-3
  do i=nx-2*ngh+1,nx2
     x_eta_imax_jmax(i,j)= a7(1)*(xc(i,j+1)-xc(i,j-1)) &
                         + a7(2)*(xc(i,j+2)-xc(i,j-2)) &
                         + a7(3)*(xc(i,j+3)-xc(i,j-3))
     y_eta_imax_jmax(i,j)= a7(1)*(yc(i,j+1)-yc(i,j-1)) &
                         + a7(2)*(yc(i,j+2)-yc(i,j-2)) &
                         + a7(3)*(yc(i,j+3)-yc(i,j-3))
  enddo

  ! Sixth-order (7-pt stencil) for point i=nx-3
  ! ===========================================
  i=nx-3
  do j=ny-2*ngh+1,ny2
     x_ksi_imax_jmax(i,j)= a7(1)*(xc(i+1,j)-xc(i-1,j)) &
                         + a7(2)*(xc(i+2,j)-xc(i-2,j)) &
                         + a7(3)*(xc(i+3,j)-xc(i-3,j))
     y_ksi_imax_jmax(i,j)= a7(1)*(yc(i+1,j)-yc(i-1,j)) &
                         + a7(2)*(yc(i+2,j)-yc(i-2,j)) &
                         + a7(3)*(yc(i+3,j)-yc(i-3,j))
  enddo

  if (ngh<5) return

  ! Eighth-order (9-pt stencil) for point j=ny-4
  ! ============================================
  j=ny-4
  do i=nx-2*ngh+1,nx2
     x_eta_imax_jmax(i,j)= a9(1)*(xc(i,j+1)-xc(i,j-1)) &
                         + a9(2)*(xc(i,j+2)-xc(i,j-2)) &
                         + a9(3)*(xc(i,j+3)-xc(i,j-3)) &
                         + a9(4)*(xc(i,j+4)-xc(i,j-4))
     y_eta_imax_jmax(i,j)= a9(1)*(yc(i,j+1)-yc(i,j-1)) &
                         + a9(2)*(yc(i,j+2)-yc(i,j-2)) &
                         + a9(3)*(yc(i,j+3)-yc(i,j-3)) &
                         + a9(4)*(yc(i,j+4)-yc(i,j-4))
  enddo

  ! Eighth-order (9-pt stencil) for point i=nx-4
  ! ============================================
  i=nx-4
  do j=ny-2*ngh+1,ny2
     x_ksi_imax_jmax(i,j)= a9(1)*(xc(i+1,j)-xc(i-1,j)) &
                         + a9(2)*(xc(i+2,j)-xc(i-2,j)) &
                         + a9(3)*(xc(i+3,j)-xc(i-3,j)) &
                         + a9(4)*(xc(i+4,j)-xc(i-4,j))
     y_ksi_imax_jmax(i,j)= a9(1)*(yc(i+1,j)-yc(i-1,j)) &
                         + a9(2)*(yc(i+2,j)-yc(i-2,j)) &
                         + a9(3)*(yc(i+3,j)-yc(i-3,j)) &
                         + a9(4)*(yc(i+4,j)-yc(i-4,j))
  enddo

end subroutine grid_metrics_wall_imax_jmax

!===============================================================================
subroutine grid_metrics_wall_imin_jmin_SBP4
!===============================================================================
  !> Compute curvilinear metrics for wall edge imin-jmin
  !> (wall BC with SBP4: Summation by Parts order 4)
!===============================================================================
  use mod_grid
  use mod_coeff_deriv
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j
  ! ---------------------------------------------------------------------------

  allocate(x_eta_imin_jmin(nx1:2*ngh,1:ngh),y_eta_imin_jmin(nx1:2*ngh,1:ngh))
  allocate(x_ksi_imin_jmin(1:ngh,ny1:2*ngh),y_ksi_imin_jmin(1:ngh,ny1:2*ngh))

  ! Point #1: stencil [o x x x] for point j=1
  ! ===========================
  !j=1
  do i=nx1,2*ngh
     x_eta_imin_jmin(i,1)= as4p0(1)*xc(i,1)+as4p0(2)*xc(i,2) &
                         + as4p0(3)*xc(i,3)+as4p0(4)*xc(i,4)
     y_eta_imin_jmin(i,1)= as4p0(1)*yc(i,1)+as4p0(2)*yc(i,2) &
                         + as4p0(3)*yc(i,3)+as4p0(4)*yc(i,4)
  enddo

  ! Point #1: stencil [o x x x] for point i=1
  ! ===========================
  !i=1
  do j=ny1,2*ngh
     x_ksi_imin_jmin(1,j)= as4p0(1)*xc(1,j)+as4p0(2)*xc(2,j) &
                         + as4p0(3)*xc(3,j)+as4p0(4)*xc(4,j)
     y_ksi_imin_jmin(1,j)= as4p0(1)*yc(1,j)+as4p0(2)*yc(2,j) &
                         + as4p0(3)*yc(3,j)+as4p0(4)*yc(4,j)
  enddo

  ! Point #2: stencil [x o x x] for point j=2
  ! ===========================
  !j=2
  do i=nx1,2*ngh
     x_eta_imin_jmin(i,2)= as4p1(1)*xc(i,1)+as4p1(2)*xc(i,2) &
                         + as4p1(3)*xc(i,3)+as4p1(4)*xc(i,4)
     y_eta_imin_jmin(i,2)= as4p1(1)*yc(i,1)+as4p1(2)*yc(i,2) &
                         + as4p1(3)*yc(i,3)+as4p1(4)*yc(i,4)
  enddo

  ! Point #2: stencil [x o x x] for point i=2
  ! ===========================
  !i=2
  do j=ny1,2*ngh
     x_ksi_imin_jmin(2,j)= as4p1(1)*xc(1,j)+as4p1(2)*xc(2,j) &
                         + as4p1(3)*xc(3,j)+as4p1(4)*xc(4,j)
     y_ksi_imin_jmin(2,j)= as4p1(1)*yc(1,j)+as4p1(2)*yc(2,j) &
                         + as4p1(3)*yc(3,j)+as4p1(4)*yc(4,j)
  enddo

  ! Point #3: stencil [x x o x x] for point j=3
  ! =============================
  !j=3
  do i=nx1,2*ngh
     x_eta_imin_jmin(i,3)= as4p2(1)*xc(i,1)+as4p2(2)*xc(i,2) &
                         + as4p2(3)*xc(i,3)+as4p2(4)*xc(i,4) &
                         + as4p2(5)*xc(i,5)
     y_eta_imin_jmin(i,3)= as4p2(1)*yc(i,1)+as4p2(2)*yc(i,2) &
                         + as4p2(3)*yc(i,3)+as4p2(4)*yc(i,4) &
                         + as4p2(5)*yc(i,5)
  enddo
     
  ! Point #3: stencil [x x o x x] for point i=3
  ! =============================
  !i=3
  do j=ny1,2*ngh
     x_ksi_imin_jmin(3,j)= as4p2(1)*xc(1,j)+as4p2(2)*xc(2,j) &
                         + as4p2(3)*xc(3,j)+as4p2(4)*xc(4,j) &
                         + as4p2(5)*xc(5,j)
     y_ksi_imin_jmin(3,j)= as4p2(1)*yc(1,j)+as4p2(2)*yc(2,j) &
                         + as4p2(3)*yc(3,j)+as4p2(4)*yc(4,j) &
                         + as4p2(5)*yc(5,j)
  enddo
     
  ! Point #4: stencil [x x x o x x] for point j=4
  ! ===============================
  !j=4
  do i=nx1,2*ngh
     x_eta_imin_jmin(i,4)= as4p3(1)*xc(i,1)+as4p3(2)*xc(i,2) &
                         + as4p3(3)*xc(i,3)+as4p3(4)*xc(i,4) &
                         + as4p3(5)*xc(i,5)+as4p3(6)*xc(i,6)
     y_eta_imin_jmin(i,4)= as4p3(1)*yc(i,1)+as4p3(2)*yc(i,2) &
                         + as4p3(3)*yc(i,3)+as4p3(4)*yc(i,4) &
                         + as4p3(5)*yc(i,5)+as4p3(6)*yc(i,6)
  enddo

  ! Point #4: stencil [x x x o x x] for point i=4
  ! ===============================
  !i=4
  do j=ny1,2*ngh
     x_ksi_imin_jmin(4,j)= as4p3(1)*xc(1,j)+as4p3(2)*xc(2,j) &
                         + as4p3(3)*xc(3,j)+as4p3(4)*xc(4,j) &
                         + as4p3(5)*xc(5,j)+as4p3(6)*xc(6,j)
     y_ksi_imin_jmin(4,j)= as4p3(1)*yc(1,j)+as4p3(2)*yc(2,j) &
                         + as4p3(3)*yc(3,j)+as4p3(4)*yc(4,j) &
                         + as4p3(5)*yc(5,j)+as4p3(6)*yc(6,j)
  enddo

  ! Eighth-order (9-pt stencil) for point j=5
  ! =========================================
  j=5
  do i=nx1,2*ngh
     x_eta_imin_jmin(i,j)= a9(1)*(xc(i,j+1)-xc(i,j-1)) &
                         + a9(2)*(xc(i,j+2)-xc(i,j-2)) &
                         + a9(3)*(xc(i,j+3)-xc(i,j-3)) &
                         + a9(4)*(xc(i,j+4)-xc(i,j-4))
     y_eta_imin_jmin(i,j)= a9(1)*(yc(i,j+1)-yc(i,j-1)) &
                         + a9(2)*(yc(i,j+2)-yc(i,j-2)) &
                         + a9(3)*(yc(i,j+3)-yc(i,j-3)) &
                         + a9(4)*(yc(i,j+4)-yc(i,j-4))
  enddo

  ! Eighth-order (9-pt stencil) for point i=5
  ! =========================================
  i=5
  do j=ny1,2*ngh
     x_ksi_imin_jmin(i,j)= a9(1)*(xc(i+1,j)-xc(i-1,j)) &
                         + a9(2)*(xc(i+2,j)-xc(i-2,j)) &
                         + a9(3)*(xc(i+3,j)-xc(i-3,j)) &
                         + a9(4)*(xc(i+4,j)-xc(i-4,j))
     y_ksi_imin_jmin(i,j)= a9(1)*(yc(i+1,j)-yc(i-1,j)) &
                         + a9(2)*(yc(i+2,j)-yc(i-2,j)) &
                         + a9(3)*(yc(i+3,j)-yc(i-3,j)) &
                         + a9(4)*(yc(i+4,j)-yc(i-4,j))
  enddo

end subroutine grid_metrics_wall_imin_jmin_SBP4

!===============================================================================
subroutine grid_metrics_wall_imin_jmax_SBP4
!===============================================================================
  !> Compute curvilinear metrics for wall edge imin-jmax
  !> (wall BC with SBP4: Summation by Parts order 4)
!===============================================================================
  use mod_grid
  use mod_coeff_deriv
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j
  ! ---------------------------------------------------------------------------

  allocate(x_eta_imin_jmax(nx1:2*ngh,ny-ngh+1:ny),y_eta_imin_jmax(nx1:2*ngh,ny-ngh+1:ny))
  allocate(x_ksi_imin_jmax(1:ngh,ny-2*ngh+1:ny2),y_ksi_imin_jmax(1:ngh,ny-2*ngh+1:ny2))
  
  ! Point #1: stencil [o x x x] for point j=ny
  ! ===========================
  j=ny
  do i=nx1,2*ngh
     x_eta_imin_jmax(i,j)= as4m0(1)*xc(i,ny  )+as4m0(2)*xc(i,ny-1) &
                         + as4m0(3)*xc(i,ny-2)+as4m0(4)*xc(i,ny-3)
     y_eta_imin_jmax(i,j)= as4m0(1)*yc(i,ny  )+as4m0(2)*yc(i,ny-1) &
                         + as4m0(3)*yc(i,ny-2)+as4m0(4)*yc(i,ny-3)
  enddo

  ! Point #1: stencil [o x x x] for point i=1
  ! ===========================
  !i=1
  do j=ny-2*ngh+1,ny2
     x_ksi_imin_jmax(1,j)= as4p0(1)*xc(1,j)+as4p0(2)*xc(2,j) &
                         + as4p0(3)*xc(3,j)+as4p0(4)*xc(4,j)
     y_ksi_imin_jmax(1,j)= as4p0(1)*yc(1,j)+as4p0(2)*yc(2,j) &
                         + as4p0(3)*yc(3,j)+as4p0(4)*yc(4,j)
  enddo

  ! Point #2: stencil [x o x x] for point j=ny-1
  ! ===========================
  j=ny-1
  do i=nx1,2*ngh
     x_eta_imin_jmax(i,j)= as4m1(1)*xc(i,ny  )+as4m1(2)*xc(i,ny-1) &
                         + as4m1(3)*xc(i,ny-2)+as4m1(4)*xc(i,ny-3)
     y_eta_imin_jmax(i,j)= as4m1(1)*yc(i,ny  )+as4m1(2)*yc(i,ny-1) &
                         + as4m1(3)*yc(i,ny-2)+as4m1(4)*yc(i,ny-3)
  enddo

  ! Point #2: stencil [x o x x] for point i=2
  ! ===========================
  !i=2
  do j=ny-2*ngh+1,ny2
     x_ksi_imin_jmax(2,j)= as4p1(1)*xc(1,j)+as4p1(2)*xc(2,j) &
                         + as4p1(3)*xc(3,j)+as4p1(4)*xc(4,j)
     y_ksi_imin_jmax(2,j)= as4p1(1)*yc(1,j)+as4p1(2)*yc(2,j) &
                         + as4p1(3)*yc(3,j)+as4p1(4)*yc(4,j)
  enddo

  ! Point #3: stencil [x x o x x] for point j=ny-2
  ! =============================
  j=ny-2
  do i=nx1,2*ngh
     x_eta_imin_jmax(i,j)= as4m2(1)*xc(i,ny  )+as4m2(2)*xc(i,ny-1) &
                         + as4m2(3)*xc(i,ny-2)+as4m2(4)*xc(i,ny-3) &
                         + as4m2(5)*xc(i,ny-4)
     y_eta_imin_jmax(i,j)= as4m2(1)*yc(i,ny  )+as4m2(2)*yc(i,ny-1) &
                         + as4m2(3)*yc(i,ny-2)+as4m2(4)*yc(i,ny-3) &
                         + as4m2(5)*yc(i,ny-4)
  enddo
     
  ! Point #3: stencil [x x o x x] for point i=3
  ! =============================
  !i=3
  do j=ny-2*ngh+1,ny2
     x_ksi_imin_jmax(3,j)= as4p2(1)*xc(1,j)+as4p2(2)*xc(2,j) &
                         + as4p2(3)*xc(3,j)+as4p2(4)*xc(4,j) &
                         + as4p2(5)*xc(5,j)
     y_ksi_imin_jmax(3,j)= as4p2(1)*yc(1,j)+as4p2(2)*yc(2,j) &
                         + as4p2(3)*yc(3,j)+as4p2(4)*yc(4,j) &
                         + as4p2(5)*yc(5,j)
  enddo

  ! Point #4: stencil [x x x o x x] for point j=ny-3
  ! ===============================
  j=ny-3
  do i=nx1,2*ngh
     x_eta_imin_jmax(i,j)= as4m3(1)*xc(i,ny  )+as4m3(2)*xc(i,ny-1) &
                         + as4m3(3)*xc(i,ny-2)+as4m3(4)*xc(i,ny-3) &
                         + as4m3(5)*xc(i,ny-4)+as4m3(6)*xc(i,ny-5)
     y_eta_imin_jmax(i,j)= as4m3(1)*yc(i,ny  )+as4m3(2)*yc(i,ny-1) &
                         + as4m3(3)*yc(i,ny-2)+as4m3(4)*yc(i,ny-3) &
                         + as4m3(5)*yc(i,ny-4)+as4m3(6)*yc(i,ny-5)
  enddo

  ! Point #4: stencil [x x x o x x] for point i=4
  ! ===============================
  !i=4
  do j=ny-2*ngh+1,ny2
     x_ksi_imin_jmax(4,j)= as4p3(1)*xc(1,j)+as4p3(2)*xc(2,j) &
                         + as4p3(3)*xc(3,j)+as4p3(4)*xc(4,j) &
                         + as4p3(5)*xc(5,j)+as4p3(6)*xc(6,j)
     y_ksi_imin_jmax(4,j)= as4p3(1)*yc(1,j)+as4p3(2)*yc(2,j) &
                         + as4p3(3)*yc(3,j)+as4p3(4)*yc(4,j) &
                         + as4p3(5)*yc(5,j)+as4p3(6)*yc(6,j)
  enddo

  ! Eighth-order (9-pt stencil) for point j=ny-4
  ! ============================================
  j=ny-4
  do i=nx1,2*ngh
     x_eta_imin_jmax(i,j)= a9(1)*(xc(i,j+1)-xc(i,j-1)) &
                         + a9(2)*(xc(i,j+2)-xc(i,j-2)) &
                         + a9(3)*(xc(i,j+3)-xc(i,j-3)) &
                         + a9(4)*(xc(i,j+4)-xc(i,j-4))
     y_eta_imin_jmax(i,j)= a9(1)*(yc(i,j+1)-yc(i,j-1)) &
                         + a9(2)*(yc(i,j+2)-yc(i,j-2)) &
                         + a9(3)*(yc(i,j+3)-yc(i,j-3)) &
                         + a9(4)*(yc(i,j+4)-yc(i,j-4))
  enddo
  
  ! Eighth-order (9-pt stencil) for point i=5
  ! =========================================
  i=5
  do j=ny-2*ngh+1,ny2
     x_ksi_imin_jmax(i,j)= a9(1)*(xc(i+1,j)-xc(i-1,j)) &
                         + a9(2)*(xc(i+2,j)-xc(i-2,j)) &
                         + a9(3)*(xc(i+3,j)-xc(i-3,j)) &
                         + a9(4)*(xc(i+4,j)-xc(i-4,j))
     y_ksi_imin_jmax(i,j)= a9(1)*(yc(i+1,j)-yc(i-1,j)) &
                         + a9(2)*(yc(i+2,j)-yc(i-2,j)) &
                         + a9(3)*(yc(i+3,j)-yc(i-3,j)) &
                         + a9(4)*(yc(i+4,j)-yc(i-4,j))
  enddo

end subroutine grid_metrics_wall_imin_jmax_SBP4

!===============================================================================
subroutine grid_metrics_wall_imax_jmin_SBP4
!===============================================================================
  !> Compute curvilinear metrics for wall  edge imax-jmin
  !> (wall BC with SBP4: Summation by Parts order 4)
!===============================================================================
  use mod_grid
  use mod_coeff_deriv
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j
  ! ---------------------------------------------------------------------------

  allocate(x_eta_imax_jmin(nx-2*ngh+1:nx2,1:ngh),y_eta_imax_jmin(nx-2*ngh+1:nx2,1:ngh))
  allocate(x_ksi_imax_jmin(nx-ngh+1:nx,nx1:2*ngh),y_ksi_imax_jmin(nx-ngh+1:nx,nx1:2*ngh))
 
  ! Point #1: stencil [o x x x] for point j=1
  ! ===========================
  j=1
  do i=nx-2*ngh+1,nx2
     x_eta_imax_jmin(i,1)= as4p0(1)*xc(i,1)+as4p0(2)*xc(i,2) &
                         + as4p0(3)*xc(i,3)+as4p0(4)*xc(i,4)
     y_eta_imax_jmin(i,1)= as4p0(1)*yc(i,1)+as4p0(2)*yc(i,2) &
                         + as4p0(3)*yc(i,3)+as4p0(4)*yc(i,4)
  enddo

  ! Point #1: stencil [o x x x] for point i=nx
  ! ===========================
  i=nx
  do j=ny1,2*ngh
     x_ksi_imax_jmin(i,j)= as4m0(1)*xc(nx  ,j)+as4m0(2)*xc(nx-1,j) &
                         + as4m0(3)*xc(nx-2,j)+as4m0(4)*xc(nx-3,j)
     y_ksi_imax_jmin(i,j)= as4m0(1)*yc(nx  ,j)+as4m0(2)*yc(nx-1,j) &
                         + as4m0(3)*yc(nx-2,j)+as4m0(4)*yc(nx-3,j)
  enddo

  ! Point #2: stencil [x o x x] for point j=2
  ! ===========================
  j=2
  do i=nx-2*ngh+1,nx2
     x_eta_imax_jmin(i,2)= as4p1(1)*xc(i,1)+as4p1(2)*xc(i,2) &
                         + as4p1(3)*xc(i,3)+as4p1(4)*xc(i,4)
     y_eta_imax_jmin(i,2)= as4p1(1)*yc(i,1)+as4p1(2)*yc(i,2) &
                         + as4p1(3)*yc(i,3)+as4p1(4)*yc(i,4)
  enddo

  ! Point #2: stencil [x o x x] for point i=nx-1
  ! ===========================
  i=nx-1
  do j=ny1,2*ngh
     x_ksi_imax_jmin(i,j)= as4m1(1)*xc(nx  ,j)+as4m1(2)*xc(nx-1,j) &
                         + as4m1(3)*xc(nx-2,j)+as4m1(4)*xc(nx-3,j)
     y_ksi_imax_jmin(i,j)= as4m1(1)*yc(nx  ,j)+as4m1(2)*yc(nx-1,j) &
                         + as4m1(3)*yc(nx-2,j)+as4m1(4)*yc(nx-3,j)
  enddo

  ! Point #3: stencil [x x o x x] for point j=3
  ! =============================
  j=3
  do i=nx-2*ngh+1,nx2
     x_eta_imax_jmin(i,3)= as4p2(1)*xc(i,1)+as4p2(2)*xc(i,2) &
                         + as4p2(3)*xc(i,3)+as4p2(4)*xc(i,4) &
                         + as4p2(5)*xc(i,5)
     y_eta_imax_jmin(i,3)= as4p2(1)*yc(i,1)+as4p2(2)*yc(i,2) &
                         + as4p2(3)*yc(i,3)+as4p2(4)*yc(i,4) &
                         + as4p2(5)*yc(i,5)
  enddo
     
  ! Point #3: stencil [x x o x x] for point i=nx-2
  ! =============================
  i=nx-2
  do j=ny1,2*ngh
     x_ksi_imax_jmin(i,j)= as4m2(1)*xc(nx  ,j)+as4m2(2)*xc(nx-1,j) &
                         + as4m2(3)*xc(nx-2,j)+as4m2(4)*xc(nx-3,j) &
                         + as4m2(5)*xc(nx-4,j)
     y_ksi_imax_jmin(i,j)= as4m2(1)*yc(nx  ,j)+as4m2(2)*yc(nx-1,j) &
                         + as4m2(3)*yc(nx-2,j)+as4m2(4)*yc(nx-3,j) &
                         + as4m2(5)*yc(nx-4,j)
  enddo
     
  ! Point #4: stencil [x x x o x x] for point j=4
  ! ===============================
  j=4
  do i=nx-2*ngh+1,nx2
     x_eta_imax_jmin(i,4)= as4p3(1)*xc(i,1)+as4p3(2)*xc(i,2) &
                         + as4p3(3)*xc(i,3)+as4p3(4)*xc(i,4) &
                         + as4p3(5)*xc(i,5)+as4p3(6)*xc(i,6)
     y_eta_imax_jmin(i,4)= as4p3(1)*yc(i,1)+as4p3(2)*yc(i,2) &
                         + as4p3(3)*yc(i,3)+as4p3(4)*yc(i,4) &
                         + as4p3(5)*yc(i,5)+as4p3(6)*yc(i,6)
  enddo

  ! Point #4: stencil [x x x o x x] for point i=nx-3
  ! ===============================
  i=nx-3
  do j=ny1,2*ngh
     x_ksi_imax_jmin(i,j)= as4m3(1)*xc(nx  ,j)+as4m3(2)*xc(nx-1,j) &
                         + as4m3(3)*xc(nx-2,j)+as4m3(4)*xc(nx-3,j) &
                         + as4m3(5)*xc(nx-4,j)+as4m3(6)*xc(nx-5,j)
     y_ksi_imax_jmin(i,j)= as4m3(1)*yc(nx  ,j)+as4m3(2)*yc(nx-1,j) &
                         + as4m3(3)*yc(nx-2,j)+as4m3(4)*yc(nx-3,j) &
                         + as4m3(5)*yc(nx-4,j)+as4m3(6)*yc(nx-5,j)
  enddo

  ! Eighth-order (9-pt stencil) for point j=5
  ! =========================================
  j=5
  do i=nx-2*ngh+1,nx2
     x_eta_imax_jmin(i,j)= a9(1)*(xc(i,j+1)-xc(i,j-1)) &
                         + a9(2)*(xc(i,j+2)-xc(i,j-2)) &
                         + a9(3)*(xc(i,j+3)-xc(i,j-3)) &
                         + a9(4)*(xc(i,j+4)-xc(i,j-4))
     y_eta_imax_jmin(i,j)= a9(1)*(yc(i,j+1)-yc(i,j-1)) &
                         + a9(2)*(yc(i,j+2)-yc(i,j-2)) &
                         + a9(3)*(yc(i,j+3)-yc(i,j-3)) &
                         + a9(4)*(yc(i,j+4)-yc(i,j-4))
  enddo

  ! Eighth-order (9-pt stencil) for point i=nx-4
  ! ============================================
  i=nx-4
  do j=ny1,2*ngh
     x_ksi_imax_jmin(i,j)= a9(1)*(xc(i+1,j)-xc(i-1,j)) &
                         + a9(2)*(xc(i+2,j)-xc(i-2,j)) &
                         + a9(3)*(xc(i+3,j)-xc(i-3,j)) &
                         + a9(4)*(xc(i+4,j)-xc(i-4,j))
     y_ksi_imax_jmin(i,j)= a9(1)*(yc(i+1,j)-yc(i-1,j)) &
                         + a9(2)*(yc(i+2,j)-yc(i-2,j)) &
                         + a9(3)*(yc(i+3,j)-yc(i-3,j)) &
                         + a9(4)*(yc(i+4,j)-yc(i-4,j))
  enddo

end subroutine grid_metrics_wall_imax_jmin_SBP4

!===============================================================================
subroutine grid_metrics_wall_imax_jmax_SBP4
!===============================================================================
  !> Compute curvilinear metrics for wall edge imax-jmax
  !> (wall BC with SBP4: Summation by Parts order 4)
!===============================================================================
  use mod_grid
  use mod_coeff_deriv
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j
  ! ---------------------------------------------------------------------------

  allocate(x_eta_imax_jmax(nx-2*ngh+1:nx2,ny-ngh+1:ny),y_eta_imax_jmax(nx-2*ngh+1:nx2,ny-ngh+1:ny))
  allocate(x_ksi_imax_jmax(nx-ngh+1:nx,ny-2*ngh+1:ny2),y_ksi_imax_jmax(nx-ngh+1:nx,ny-2*ngh+1:ny2))
 
  ! Point #1: stencil [o x x x] for point j=ny
  ! ===========================
  j=ny
  do i=nx-2*ngh+1,nx2
     x_eta_imax_jmax(i,j)= as4m0(1)*xc(i,ny  )+as4m0(2)*xc(i,ny-1) &
                         + as4m0(3)*xc(i,ny-2)+as4m0(4)*xc(i,ny-3)
     y_eta_imax_jmax(i,j)= as4m0(1)*yc(i,ny  )+as4m0(2)*yc(i,ny-1) &
                         + as4m0(3)*yc(i,ny-2)+as4m0(4)*yc(i,ny-3)
  enddo

  ! Point #1: stencil [o x x x] for point i=nx
  ! ===========================
  i=nx
  do j=ny-2*ngh+1,ny2
     x_ksi_imax_jmax(i,j)= as4m0(1)*xc(nx  ,j)+as4m0(2)*xc(nx-1,j) &
                         + as4m0(3)*xc(nx-2,j)+as4m0(4)*xc(nx-3,j)
     y_ksi_imax_jmax(i,j)= as4m0(1)*yc(nx  ,j)+as4m0(2)*yc(nx-1,j) &
                         + as4m0(3)*yc(nx-2,j)+as4m0(4)*yc(nx-3,j)
  enddo

  ! Point #2: stencil [x o x x] for point j=ny-1
  ! ===========================
  j=ny-1
  do i=nx-2*ngh+1,nx2
     x_eta_imax_jmax(i,j)= as4m1(1)*xc(i,ny  )+as4m1(2)*xc(i,ny-1) &
                         + as4m1(3)*xc(i,ny-2)+as4m1(4)*xc(i,ny-3)
     y_eta_imax_jmax(i,j)= as4m1(1)*yc(i,ny  )+as4m1(2)*yc(i,ny-1) &
                         + as4m1(3)*yc(i,ny-2)+as4m1(4)*yc(i,ny-3)
  enddo

  ! Point #2: stencil [x o x x] for point i=nx-1
  ! ===========================
  i=nx-1
  do j=ny-2*ngh+1,ny2
     x_ksi_imax_jmax(i,j)= as4m1(1)*xc(nx  ,j)+as4m1(2)*xc(nx-1,j) &
                         + as4m1(3)*xc(nx-2,j)+as4m1(4)*xc(nx-3,j)
     y_ksi_imax_jmax(i,j)= as4m1(1)*yc(nx  ,j)+as4m1(2)*yc(nx-1,j) &
                         + as4m1(3)*yc(nx-2,j)+as4m1(4)*yc(nx-3,j)
  enddo

  ! Point #3: stencil [x x o x x] for point j=ny-2
  ! =============================
  j=ny-2
  do i=nx-2*ngh+1,nx2
     x_eta_imax_jmax(i,j)= as4m2(1)*xc(i,ny  )+as4m2(2)*xc(i,ny-1) &
                         + as4m2(3)*xc(i,ny-2)+as4m2(4)*xc(i,ny-3) &
                         + as4m2(5)*xc(i,ny-4)
     y_eta_imax_jmax(i,j)= as4m2(1)*yc(i,ny  )+as4m2(2)*yc(i,ny-1) &
                         + as4m2(3)*yc(i,ny-2)+as4m2(4)*yc(i,ny-3) &
                         + as4m2(5)*yc(i,ny-4)
  enddo
     
  ! Point #3: stencil [x x o x x] for point i=nx-2
  ! =============================
  i=nx-2
  do j=ny-2*ngh+1,ny2
     x_ksi_imax_jmax(i,j)= as4m2(1)*xc(nx  ,j)+as4m2(2)*xc(nx-1,j) &
                         + as4m2(3)*xc(nx-2,j)+as4m2(4)*xc(nx-3,j) &
                         + as4m2(5)*xc(nx-4,j)
     y_ksi_imax_jmax(i,j)= as4m2(1)*yc(nx  ,j)+as4m2(2)*yc(nx-1,j) &
                         + as4m2(3)*yc(nx-2,j)+as4m2(4)*yc(nx-3,j) &
                         + as4m2(5)*yc(nx-4,j)
  enddo
     
  ! Point #4: stencil [x x x o x x] for point j=ny-3
  ! ===============================
  j=ny-3
  do i=nx-2*ngh+1,nx2
    x_eta_imax_jmax(i,j)= as4m3(1)*xc(i,ny  )+as4m3(2)*xc(i,ny-1) &
                         + as4m3(3)*xc(i,ny-2)+as4m3(4)*xc(i,ny-3) &
                         + as4m3(5)*xc(i,ny-4)+as4m3(6)*xc(i,ny-5)
     y_eta_imax_jmax(i,j)= as4m3(1)*yc(i,ny  )+as4m3(2)*yc(i,ny-1) &
                         + as4m3(3)*yc(i,ny-2)+as4m3(4)*yc(i,ny-3) &
                         + as4m3(5)*yc(i,ny-4)+as4m3(6)*yc(i,ny-5)
  enddo

  ! Point #4: stencil [x x x o x x] for point i=nx-3
  ! ===============================
  i=nx-3
  do j=ny-2*ngh+1,ny2
     x_ksi_imax_jmax(i,j)= as4m3(1)*xc(nx  ,j)+as4m3(2)*xc(nx-1,j) &
                         + as4m3(3)*xc(nx-2,j)+as4m3(4)*xc(nx-3,j) &
                         + as4m3(5)*xc(nx-4,j)+as4m3(6)*xc(nx-5,j)
     y_ksi_imax_jmax(i,j)= as4m3(1)*yc(nx  ,j)+as4m3(2)*yc(nx-1,j) &
                         + as4m3(3)*yc(nx-2,j)+as4m3(4)*yc(nx-3,j) &
                         + as4m3(5)*yc(nx-4,j)+as4m3(6)*yc(nx-5,j)
  enddo

  ! Eighth-order (9-pt stencil) for point j=ny-4
  ! ============================================
  j=ny-4
  do i=nx-2*ngh+1,nx2
     x_eta_imax_jmax(i,j)= a9(1)*(xc(i,j+1)-xc(i,j-1)) &
                         + a9(2)*(xc(i,j+2)-xc(i,j-2)) &
                         + a9(3)*(xc(i,j+3)-xc(i,j-3)) &
                         + a9(4)*(xc(i,j+4)-xc(i,j-4))
     y_eta_imax_jmax(i,j)= a9(1)*(yc(i,j+1)-yc(i,j-1)) &
                         + a9(2)*(yc(i,j+2)-yc(i,j-2)) &
                         + a9(3)*(yc(i,j+3)-yc(i,j-3)) &
                         + a9(4)*(yc(i,j+4)-yc(i,j-4))
  enddo

  ! Eighth-order (9-pt stencil) for point i=nx-4
  ! ============================================
  i=nx-4
  do j=ny-2*ngh+1,ny2
     x_ksi_imax_jmax(i,j)= a9(1)*(xc(i+1,j)-xc(i-1,j)) &
                         + a9(2)*(xc(i+2,j)-xc(i-2,j)) &
                         + a9(3)*(xc(i+3,j)-xc(i-3,j)) &
                         + a9(4)*(xc(i+4,j)-xc(i-4,j))
     y_ksi_imax_jmax(i,j)= a9(1)*(yc(i+1,j)-yc(i-1,j)) &
                         + a9(2)*(yc(i+2,j)-yc(i-2,j)) &
                         + a9(3)*(yc(i+3,j)-yc(i-3,j)) &
                         + a9(4)*(yc(i+4,j)-yc(i-4,j))
  enddo

end subroutine grid_metrics_wall_imax_jmax_SBP4
