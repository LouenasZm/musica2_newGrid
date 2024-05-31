! ==============================================================================
module mod_deriv2d
! ==============================================================================
  !> Module for 2-D derivatives
! ==============================================================================
  use mod_grid
  use mod_bc
  use mod_coeff_deriv
  implicit none
  ! -----------------------------------------------------------------------------
  integer, private :: i,j,l
  ! -----------------------------------------------------------------------------

contains

  !==============================================================================
  subroutine deriv2_x_11pts(var,dvar)
  !==============================================================================
    !> derivatives of a 2-D variable along x on 11-point stencil
  !==============================================================================  
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(1:nx,1:ny), intent(in) :: var
    real(wp), dimension(1:nx,1:ny), intent(out) :: dvar
    ! ---------------------------------------------------------------------------

    dvar = 0.0_wp

    if (is_boundary(1,1)) then
       i=1
       do j=1,ny
          do l=0,10
             dvar(i,j)=dvar(i,j)+a010(1+l)*var(i+l,j)*idx(i)
          enddo
       enddo

       i=2
       do j=1,ny
          do l=-1,9
             dvar(i,j)=dvar(i,j)+a19(2+l)*var(i+l,j)*idx(i)
          enddo
       enddo

       i=3
       do j=1,ny
          do l=-2,8
             dvar(i,j)=dvar(i,j)+a28(3+l)*var(i+l,j)*idx(i)
          enddo
       enddo

       i=4
       do j=1,ny
          do l=-3,7
             dvar(i,j)=dvar(i,j)+a37(4+l)*var(i+l,j)*idx(i)
          enddo
       enddo

       i=5
       do j=1,ny
          do l=-4,6
             dvar(i,j)=dvar(i,j)+a46(5+l)*var(i+l,j)*idx(i)
          enddo
       enddo
    endif

    do j=1,ny
       do i=ndx,nfx
          dvar(i,j)=( a11(1)*(var(i+1,j)-var(i-1,j)) &
               + a11(2)*(var(i+2,j)-var(i-2,j)) &
               + a11(3)*(var(i+3,j)-var(i-3,j)) &
               + a11(4)*(var(i+4,j)-var(i-4,j)) &
               + a11(5)*(var(i+5,j)-var(i-5,j)) )*idx(i)
       enddo
    enddo

    if (is_boundary(1,2)) then
       i=nx-4
       do j=1,ny
          do l=-6,4
             dvar(i,j)=dvar(i,j)-a46(5-l)*var(i+l,j)*idx(i)
          enddo
       enddo

       i=nx-3
       do j=1,ny
          do l=-7,3
             dvar(i,j)=dvar(i,j)-a37(4-l)*var(i+l,j)*idx(i)
          enddo
       enddo

       i=nx-2
       do j=1,ny
          do l=-8,2
             dvar(i,j)=dvar(i,j)-a28(3-l)*var(i+l,j)*idx(i)
          enddo
       enddo

       i=nx-1
       do j=1,ny
          do l=-9,1
             dvar(i,j)=dvar(i,j)-a19(2-l)*var(i+l,j)*idx(i)
          enddo
       enddo

       i=nx
       do j=1,ny
          do l=-10,0
             dvar(i,j)=dvar(i,j)-a010(1-l)*var(i+l,j)*idx(i)
          enddo
       enddo
    endif

  end subroutine deriv2_x_11pts

  !==============================================================================
  subroutine deriv2_y_11pts(var,dvar)
  !==============================================================================
    !> derivatives of a 2-D variable along y on 11-point stencil
  !==============================================================================  
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(1:nx,1:ny), intent(in) :: var
    real(wp), dimension(1:nx,1:ny), intent(out) :: dvar
    ! ---------------------------------------------------------------------------

    dvar = 0.0_wp

    if (is_boundary(2,1)) then
       j=1
       do i=1,nx
          do l=0,10
             dvar(i,j)=dvar(i,j)+a010(1+l)*var(i,j+l)*idy(j)
          enddo
       enddo

       j=2
       do i=1,nx
          do l=-1,9
             dvar(i,j)=dvar(i,j)+a19(2+l)*var(i,j+l)*idy(j)
          enddo
       enddo

       j=3
       do i=1,nx
          do l=-2,8
             dvar(i,j)=dvar(i,j)+a28(3+l)*var(i,j+l)*idy(j)
          enddo
       enddo

       j=4
       do i=1,nx
          do l=-3,7
             dvar(i,j)=dvar(i,j)+a37(4+l)*var(i,j+l)*idy(j)
          enddo
       enddo

       j=5
       do i=1,nx
          do l=-4,6
             dvar(i,j)=dvar(i,j)+a46(5+l)*var(i,j+l)*idy(j)
          enddo
       enddo
    endif

    do j=ndy,nfy
       do i=1,nx
          dvar(i,j)=( a11(1)*(var(i,j+1)-var(i,j-1)) &
               + a11(2)*(var(i,j+2)-var(i,j-2)) &
               + a11(3)*(var(i,j+3)-var(i,j-3)) &
               + a11(4)*(var(i,j+4)-var(i,j-4)) &
               + a11(5)*(var(i,j+5)-var(i,j-5)) )*idy(j)
       enddo
    enddo

    if (is_boundary(2,2)) then
       j=ny-4
       do i=1,nx
          do l=-6,4
             dvar(i,j)=dvar(i,j)-a46(5-l)*var(i,j+l)*idy(j)
          enddo
       enddo

       j=ny-3
       do i=1,nx
          do l=-7,3
             dvar(i,j)=dvar(i,j)-a37(4-l)*var(i,j+l)*idy(j)
          enddo
       enddo

       j=ny-2
       do i=1,nx
          do l=-8,2
             dvar(i,j)=dvar(i,j)-a28(3-l)*var(i,j+l)*idy(j)
          enddo
       enddo

       j=ny-1
       do i=1,nx
          do l=-9,1
             dvar(i,j)=dvar(i,j)-a19(2-l)*var(i,j+l)*idy(j)
          enddo
       enddo

       j=ny
       do i=1,nx
          do l=-10,0
             dvar(i,j)=dvar(i,j)-a010(1-l)*var(i,j+l)*idy(j)
          enddo
       enddo
    endif

  end subroutine deriv2_y_11pts

  !==============================================================================
  subroutine deriv2_x_9pts(var,dvar)
  !==============================================================================
    !> derivatives of a 2-D variable along x on 9-point stencil
  !==============================================================================  
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(1:nx,1:ny), intent(in) :: var
    real(wp), dimension(1:nx,1:ny), intent(out) :: dvar
    ! ---------------------------------------------------------------------------

    dvar = 0.0_wp

    if (is_boundary(1,1)) then
       i=1
       do j=1,ny
          do l=0,8
             dvar(i,j)=dvar(i,j)+a08(1+l)*var(i+l,j)*idx(i)
          enddo
       enddo

       i=2
       do j=1,ny
          do l=-1,7
             dvar(i,j)=dvar(i,j)+a17(2+l)*var(i+l,j)*idx(i)
          enddo
       enddo

       i=3
       do j=1,ny
          do l=-2,6
             dvar(i,j)=dvar(i,j)+a26(3+l)*var(i+l,j)*idx(i)
          enddo
       enddo

       i=4
       do j=1,ny
          do l=-3,5
             dvar(i,j)=dvar(i,j)+a35(4+l)*var(i+l,j)*idx(i)
          enddo
       enddo
    endif

    do j=1,ny
       do i=ndx,nfx
          dvar(i,j)=( a9(1)*(var(i+1,j)-var(i-1,j)) &
               + a9(2)*(var(i+2,j)-var(i-2,j)) &
               + a9(3)*(var(i+3,j)-var(i-3,j)) &
               + a9(4)*(var(i+4,j)-var(i-4,j)) )*idx(i)
       enddo
    enddo

    if (is_boundary(1,2)) then
       i=nx-3
       do j=1,ny
          do l=-5,3
             dvar(i,j)=dvar(i,j)-a35(4-l)*var(i+l,j)*idx(i)
          enddo
       enddo

       i=nx-2
       do j=1,ny
          do l=-6,2
             dvar(i,j)=dvar(i,j)-a26(3-l)*var(i+l,j)*idx(i)
          enddo
       enddo

       i=nx-1
       do j=1,ny
          do l=-7,1
             dvar(i,j)=dvar(i,j)-a17(2-l)*var(i+l,j)*idx(i)
          enddo
       enddo

       i=nx
       do j=1,ny
          do l=-8,0
             dvar(i,j)=dvar(i,j)-a08(1-l)*var(i+l,j)*idx(i)
          enddo
       enddo
    endif

  end subroutine deriv2_x_9pts

  !==============================================================================
  subroutine deriv2_y_9pts(var,dvar)
  !==============================================================================
    !> derivatives of a 2-D variable along y on 9-point stencil
  !==============================================================================  
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(1:nx,1:ny), intent(in) :: var
    real(wp), dimension(1:nx,1:ny), intent(out) :: dvar
    ! ---------------------------------------------------------------------------

    dvar = 0.0_wp

    if (is_boundary(2,1)) then
       j=1
       do i=1,nx
          do l=0,8
             dvar(i,j)=dvar(i,j)+a08(1+l)*var(i,j+l)*idy(j)
          enddo
       enddo

       j=2
       do i=1,nx
          do l=-1,7
             dvar(i,j)=dvar(i,j)+a17(2+l)*var(i,j+l)*idy(j)
          enddo
       enddo

       j=3
       do i=1,nx
          do l=-2,6
             dvar(i,j)=dvar(i,j)+a26(3+l)*var(i,j+l)*idy(j)
          enddo
       enddo

       j=4
       do i=1,nx
          do l=-3,5
             dvar(i,j)=dvar(i,j)+a35(4+l)*var(i,j+l)*idy(j)
          enddo
       enddo
    endif

    do j=ndy,nfy
       do i=1,nx
          dvar(i,j)=( a9(1)*(var(i,j+1)-var(i,j-1)) &
               + a9(2)*(var(i,j+2)-var(i,j-2)) &
               + a9(3)*(var(i,j+3)-var(i,j-3)) &
               + a9(4)*(var(i,j+4)-var(i,j-4)) )*idy(j)
       enddo
    enddo

    if (is_boundary(2,2)) then
       j=ny-3
       do i=1,nx
          do l=-5,3
             dvar(i,j)=dvar(i,j)-a35(4-l)*var(i,j+l)*idy(j)
          enddo
       enddo

       j=ny-2
       do i=1,nx
          do l=-6,2
             dvar(i,j)=dvar(i,j)-a26(3-l)*var(i,j+l)*idy(j)
          enddo
       enddo

       j=ny-1
       do i=1,nx
          do l=-7,1
             dvar(i,j)=dvar(i,j)-a17(2-l)*var(i,j+l)*idy(j)
          enddo
       enddo

       j=ny
       do i=1,nx
          do l=-8,0
             dvar(i,j)=dvar(i,j)-a08(1-l)*var(i,j+l)*idy(j)
          enddo
       enddo
    endif

  end subroutine deriv2_y_9pts

  !==============================================================================
  subroutine deriv2_x_7pts(var,dvar)
  !==============================================================================
    !> derivatives of a 2-D variable along x on 7-point stencil
  !==============================================================================  
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(1:nx,1:ny), intent(in) :: var
    real(wp), dimension(1:nx,1:ny), intent(out) :: dvar
    ! ---------------------------------------------------------------------------

    dvar = 0.0_wp

    if (is_boundary(1,1)) then
       i=1
       do j=1,ny
          do l=0,6
             dvar(i,j)=dvar(i,j)+a06(1+l)*var(i+l,j)*idx(i)
          enddo
       enddo

       i=2
       do j=1,ny
          do l=-1,5
             dvar(i,j)=dvar(i,j)+a15(2+l)*var(i+l,j)*idx(i)
          enddo
       enddo

       i=3
       do j=1,ny
          do l=-2,4
             dvar(i,j)=dvar(i,j)+a24(3+l)*var(i+l,j)*idx(i)
          enddo
       enddo
    endif

    do j=1,ny
       do i=ndx,nfx
          dvar(i,j)=( a7(1)*(var(i+1,j)-var(i-1,j)) &
               + a7(2)*(var(i+2,j)-var(i-2,j)) &
               + a7(3)*(var(i+3,j)-var(i-3,j)) )*idx(i)
       enddo
    enddo

    if (is_boundary(1,2)) then
       i=nx-2
       do j=1,ny
          do l=-4,2
             dvar(i,j)=dvar(i,j)-a24(3-l)*var(i+l,j)*idx(i)
          enddo
       enddo

       i=nx-1
       do j=1,ny
          do l=-5,1
             dvar(i,j)=dvar(i,j)-a15(2-l)*var(i+l,j)*idx(i)
          enddo
       enddo

       i=nx
       do j=1,ny
          do l=-6,0
             dvar(i,j)=dvar(i,j)-a06(1-l)*var(i+l,j)*idx(i)
          enddo
       enddo
    endif

  end subroutine deriv2_x_7pts

  !==============================================================================
  subroutine deriv2_y_7pts(var,dvar)
  !==============================================================================
    !> derivatives of a 2-D variable along y on 7-point stencil
  !==============================================================================  
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(1:nx,1:ny), intent(in) :: var
    real(wp), dimension(1:nx,1:ny), intent(out) :: dvar
    ! ---------------------------------------------------------------------------

    dvar = 0.0_wp

    if (is_boundary(2,1)) then
       j=1
       do i=1,nx
          do l=0,6
             dvar(i,j)=dvar(i,j)+a06(1+l)*var(i,j+l)*idy(j)
          enddo
       enddo

       j=2
       do i=1,nx
          do l=-1,5
             dvar(i,j)=dvar(i,j)+a15(2+l)*var(i,j+l)*idy(j)
          enddo
       enddo

       j=3
       do i=1,nx
          do l=-2,4
             dvar(i,j)=dvar(i,j)+a24(3+l)*var(i,j+l)*idy(j)
          enddo
       enddo
    endif

    do j=ndy,nfy
       do i=1,nx
          dvar(i,j)=( a7(1)*(var(i,j+1)-var(i,j-1)) &
               + a7(2)*(var(i,j+2)-var(i,j-2)) &
               + a7(3)*(var(i,j+3)-var(i,j-3)) )*idy(j)
       enddo
    enddo

    if (is_boundary(2,2)) then
       j=ny-2
       do i=1,nx
          do l=-4,2
             dvar(i,j)=dvar(i,j)-a24(3-l)*var(i,j+l)*idy(j)
          enddo
       enddo

       j=ny-1
       do i=1,nx
          do l=-5,1
             dvar(i,j)=dvar(i,j)-a15(2-l)*var(i,j+l)*idy(j)
          enddo
       enddo

       j=ny
       do i=1,nx
          do l=-6,0
             dvar(i,j)=dvar(i,j)-a06(1-l)*var(i,j+l)*idy(j)
          enddo
       enddo
    endif

  end subroutine deriv2_y_7pts

  !==============================================================================
  subroutine deriv2_x_5pts(var,dvar)
  !==============================================================================
    !> derivatives of a 2-D variable along x on 5-point stencil
  !==============================================================================  
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(1:nx,1:ny), intent(in) :: var
    real(wp), dimension(1:nx,1:ny), intent(out) :: dvar
    ! ---------------------------------------------------------------------------

    dvar = 0.0_wp

    if (is_boundary(1,1)) then
       i=1
       do j=1,ny
          do l=0,4
             dvar(i,j)=dvar(i,j)+a04(1+l)*var(i+l,j)*idx(i)
          enddo
       enddo

       i=2
       do j=1,ny
          do l=-1,3
             dvar(i,j)=dvar(i,j)+a13(2+l)*var(i+l,j)*idx(i)
          enddo
       enddo
    endif

    do j=1,ny
       do i=ndx,nfx
          dvar(i,j)=( a5(1)*(var(i+1,j)-var(i-1,j)) &
               + a5(2)*(var(i+2,j)-var(i-2,j)) )*idx(i)
       enddo
    enddo

    if (is_boundary(1,2)) then
       i=nx-1
       do j=1,ny
          do l=-3,1
             dvar(i,j)=dvar(i,j)-a13(2-l)*var(i+l,j)*idx(i)
          enddo
       enddo

       i=nx
       do j=1,ny
          do l=-4,0
             dvar(i,j)=dvar(i,j)-a04(1-l)*var(i+l,j)*idx(i)
          enddo
       enddo
    endif

  end subroutine deriv2_x_5pts

  !==============================================================================
  subroutine deriv2_y_5pts(var,dvar)
  !==============================================================================
    !> derivatives of a 2-D variable along y on 5-point stencil
  !==============================================================================  
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(1:nx,1:ny), intent(in) :: var
    real(wp), dimension(1:nx,1:ny), intent(out) :: dvar
    ! ---------------------------------------------------------------------------

    dvar = 0.0_wp

    if (is_boundary(2,1)) then
       j=1
       do i=1,nx
          do l=0,4
             dvar(i,j)=dvar(i,j)+a04(1+l)*var(i,j+l)*idy(j)
          enddo
       enddo

       j=2
       do i=1,nx
          do l=-1,3
             dvar(i,j)=dvar(i,j)+a13(2+l)*var(i,j+l)*idy(j)
          enddo
       enddo
    endif

    do j=ndy,nfy
       do i=1,nx
          dvar(i,j)=( a5(1)*(var(i,j+1)-var(i,j-1)) &
               + a5(2)*(var(i,j+2)-var(i,j-2)) )*idy(j)
       enddo
    enddo

    if (is_boundary(2,2)) then
       j=ny-1
       do i=1,nx
          do l=-3,1
             dvar(i,j)=dvar(i,j)-a13(2-l)*var(i,j+l)*idy(j)
          enddo
       enddo

       j=ny
       do i=1,nx
          do l=-4,0
             dvar(i,j)=dvar(i,j)-a04(1-l)*var(i,j+l)*idy(j)
          enddo
       enddo
    endif

  end subroutine deriv2_y_5pts

  !==============================================================================
  subroutine deriv2_x_3pts(var,dvar)
  !==============================================================================
    !> derivatives of a 2-D variable along x on 3-point stencil
  !==============================================================================  
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(1:nx,1:ny), intent(in) :: var
    real(wp), dimension(1:nx,1:ny), intent(out) :: dvar
    ! ---------------------------------------------------------------------------

    dvar = 0.0_wp

    if (is_boundary(1,1)) then
       i=1
       do j=1,ny
          dvar(i,j)=(var(i+1,j)-var(i,j))*idx(i)
       enddo
    endif

    do j=1,ny
       do i=ndx,nfx
          dvar(i,j)=0.5_wp*(var(i+1,j)-var(i-1,j))*idx(i)
       enddo
    enddo

    if (is_boundary(1,2)) then
       i=nx
       do j=1,ny
          dvar(i,j)=(var(i,j)-var(i-1,j))*idx(i)
       enddo
    endif

  end subroutine deriv2_x_3pts

  !==============================================================================
  subroutine deriv2_y_3pts(var,dvar)
  !==============================================================================
    !> derivatives of a 2-D variable along y on 3-point stencil
  !==============================================================================  
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(1:nx,1:ny), intent(in) :: var
    real(wp), dimension(1:nx,1:ny), intent(out) :: dvar
    ! ---------------------------------------------------------------------------

    dvar = 0.0_wp

    if (is_boundary(2,1)) then
       j=1
       do i=1,nx
          dvar(i,j)=(var(i,j+1)-var(i,j))*idy(j)
       enddo
    endif

    do j=ndy,nfy
       do i=1,nx
          dvar(i,j)=0.5_wp*(var(i,j+1)-var(i,j-1))*idy(j)
       enddo
    enddo

    if (is_boundary(2,2)) then
       j=ny
       do i=1,nx
          dvar(i,j)=(var(i,j)-var(i,j-1))*idy(j)
       enddo
    endif

  end subroutine deriv2_y_3pts

end module mod_deriv2d
