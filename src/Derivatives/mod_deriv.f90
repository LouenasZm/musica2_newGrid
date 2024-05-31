! ==============================================================================
module mod_deriv
! ==============================================================================
  !> Module for derivatives
! ==============================================================================
  use mod_grid
  use mod_bc
  use mod_coeff_deriv
  implicit none
  ! -----------------------------------------------------------------------------
  integer, private :: i,j,k,l
  ! -----------------------------------------------------------------------------

contains

  !==============================================================================
  subroutine deriv_x_11pts(var,dvar)
  !==============================================================================
    !> derivatives of a variable along x on 11-point stencil
  !==============================================================================  
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(in) :: var
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(out) :: dvar
    ! ---------------------------------------------------------------------------
    integer :: ndx_,nfx_

    dvar = 0.0_wp
    ndx_=1;nfx_=nx

    if (is_boundary(1,1)) then
       ndx_=6

       i=1
       do k=1,nz
          do j=1,ny
             do l=0,10
                dvar(i,j,k)=dvar(i,j,k)+a010(1+l)*var(i+l,j,k)*idx(i)
             enddo
          enddo
       enddo

       i=2
       do k=1,nz
          do j=1,ny
             do l=-1,9
                dvar(i,j,k)=dvar(i,j,k)+a19(2+l)*var(i+l,j,k)*idx(i)
             enddo
          enddo
       enddo

       i=3
       do k=1,nz
          do j=1,ny
             do l=-2,8
                dvar(i,j,k)=dvar(i,j,k)+a28(3+l)*var(i+l,j,k)*idx(i)
             enddo
          enddo
       enddo

       i=4
       do k=1,nz
          do j=1,ny
             do l=-3,7
                dvar(i,j,k)=dvar(i,j,k)+a37(4+l)*var(i+l,j,k)*idx(i)
             enddo
          enddo
       enddo

       i=5
       do k=1,nz
          do j=1,ny
             do l=-4,6
                dvar(i,j,k)=dvar(i,j,k)+a46(5+l)*var(i+l,j,k)*idx(i)
             enddo
          enddo
       enddo
    endif

    if (is_boundary(1,2)) then
       nfx_=nx-5

       i=nx-4
       do k=1,nz
          do j=1,ny
             do l=-6,4
                dvar(i,j,k)=dvar(i,j,k)-a46(5-l)*var(i+l,j,k)*idx(i)
             enddo
          enddo
       enddo

       i=nx-3
       do k=1,nz
          do j=1,ny
             do l=-7,3
                dvar(i,j,k)=dvar(i,j,k)-a37(4-l)*var(i+l,j,k)*idx(i)
             enddo
          enddo
       enddo

       i=nx-2
       do k=1,nz
          do j=1,ny
             do l=-8,2
                dvar(i,j,k)=dvar(i,j,k)-a28(3-l)*var(i+l,j,k)*idx(i)
             enddo
          enddo
       enddo

       i=nx-1
       do k=1,nz
          do j=1,ny
             do l=-9,1
                dvar(i,j,k)=dvar(i,j,k)-a19(2-l)*var(i+l,j,k)*idx(i)
             enddo
          enddo
       enddo

       i=nx
       do k=1,nz
          do j=1,ny
             do l=-10,0
                dvar(i,j,k)=dvar(i,j,k)-a010(1-l)*var(i+l,j,k)*idx(i)
             enddo
          enddo
       enddo
    endif

    do k=1,nz
       do j=1,ny
          do i=ndx_,nfx_
             dvar(i,j,k)=( a11(1)*(var(i+1,j,k)-var(i-1,j,k)) &
                         + a11(2)*(var(i+2,j,k)-var(i-2,j,k)) &
                         + a11(3)*(var(i+3,j,k)-var(i-3,j,k)) &
                         + a11(4)*(var(i+4,j,k)-var(i-4,j,k)) &
                         + a11(5)*(var(i+5,j,k)-var(i-5,j,k)) )*idx(i)
          enddo
       enddo
    enddo

  end subroutine deriv_x_11pts

  !==============================================================================
  subroutine deriv_y_11pts(var,dvar)
  !==============================================================================
    !> derivatives of a variable along y on 11-point stencil
  !==============================================================================  
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(in) :: var
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(out) :: dvar
    ! ---------------------------------------------------------------------------
    integer :: ndy_,nfy_

    dvar = 0.0_wp
    ndy_=1;nfy_=ny

    if (is_boundary(2,1)) then
       ndy_=6

       j=1
       do k=1,nz
          do i=1,nx
             do l=0,10
                dvar(i,j,k)=dvar(i,j,k)+a010(1+l)*var(i,j+l,k)*idy(j)
             enddo
          enddo
       enddo

       j=2
       do k=1,nz
          do i=1,nx
             do l=-1,9
                dvar(i,j,k)=dvar(i,j,k)+a19(2+l)*var(i,j+l,k)*idy(j)
             enddo
          enddo
       enddo

       j=3
       do k=1,nz
          do i=1,nx
             do l=-2,8
                dvar(i,j,k)=dvar(i,j,k)+a28(3+l)*var(i,j+l,k)*idy(j)
             enddo
          enddo
       enddo

       j=4
       do k=1,nz
          do i=1,nx
             do l=-3,7
                dvar(i,j,k)=dvar(i,j,k)+a37(4+l)*var(i,j+l,k)*idy(j)
             enddo
          enddo
       enddo

       j=5
       do k=1,nz
          do i=1,nx
             do l=-4,6
                dvar(i,j,k)=dvar(i,j,k)+a46(5+l)*var(i,j+l,k)*idy(j)
             enddo
          enddo
       enddo
    endif

    if (is_boundary(2,2)) then
       nfy_=ny-5

       j=ny-4
       do k=1,nz
          do i=1,nx
             do l=-6,4
                dvar(i,j,k)=dvar(i,j,k)-a46(5-l)*var(i,j+l,k)*idy(j)
             enddo
          enddo
       enddo

       j=ny-3
       do k=1,nz
          do i=1,nx
             do l=-7,3
                dvar(i,j,k)=dvar(i,j,k)-a37(4-l)*var(i,j+l,k)*idy(j)
             enddo
          enddo
       enddo

       j=ny-2
       do k=1,nz
          do i=1,nx
             do l=-8,2
                dvar(i,j,k)=dvar(i,j,k)-a28(3-l)*var(i,j+l,k)*idy(j)
             enddo
          enddo
       enddo

       j=ny-1
       do k=1,nz
          do i=1,nx
             do l=-9,1
                dvar(i,j,k)=dvar(i,j,k)-a19(2-l)*var(i,j+l,k)*idy(j)
             enddo
          enddo
       enddo

       j=ny
       do k=1,nz
          do i=1,nx
             do l=-10,0
                dvar(i,j,k)=dvar(i,j,k)-a010(1-l)*var(i,j+l,k)*idy(j)
             enddo
          enddo
       enddo
    endif

    do k=1,nz
       do j=ndy_,nfy_
          do i=1,nx
             dvar(i,j,k)=( a11(1)*(var(i,j+1,k)-var(i,j-1,k)) &
                         + a11(2)*(var(i,j+2,k)-var(i,j-2,k)) &
                         + a11(3)*(var(i,j+3,k)-var(i,j-3,k)) &
                         + a11(4)*(var(i,j+4,k)-var(i,j-4,k)) &
                         + a11(5)*(var(i,j+5,k)-var(i,j-5,k)) )*idy(j)
          enddo
       enddo
    enddo

  end subroutine deriv_y_11pts

  !==============================================================================
  subroutine deriv_z_11pts(var,dvar)
  !==============================================================================
    !> derivatives of a variable along z on 11-point stencil
  !==============================================================================  
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(in) :: var
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(out) :: dvar
    ! ---------------------------------------------------------------------------
    integer :: ndz_,nfz_

    dvar = 0.0_wp
    ndz_=1;nfz_=nz

    if (is_boundary(3,1)) then
       ndz_=5

       k=1
       do j=1,ny
          do i=1,nx
             do l=0,10
                dvar(i,j,k)=dvar(i,j,k)+a010(1+l)*var(i,j,k+l)*idz(k)
             enddo
          enddo
       enddo

       k=2
       do j=1,ny
          do i=1,nx
             do l=-1,9
                dvar(i,j,k)=dvar(i,j,k)+a19(2+l)*var(i,j,k+l)*idz(k)
             enddo
          enddo
       enddo

       k=3
       do j=1,ny
          do i=1,nx
             do l=-2,8
                dvar(i,j,k)=dvar(i,j,k)+a28(3+l)*var(i,j,k+l)*idz(k)
             enddo
          enddo
       enddo

       k=4
       do j=1,ny
          do i=1,nx
             do l=-3,7
                dvar(i,j,k)=dvar(i,j,k)+a37(4+l)*var(i,j,k+l)*idz(k)
             enddo
          enddo
       enddo

       k=5
       do j=1,ny
          do i=1,nx
             do l=-4,6
                dvar(i,j,k)=dvar(i,j,k)+a46(5+l)*var(i,j,k+l)*idz(k)
             enddo
          enddo
       enddo
    endif

    if (is_boundary(3,2)) then
       nfz_=nz-5

       k=nz-4
       do j=1,ny
          do i=1,nx
             do l=-6,4
                dvar(i,j,k)=dvar(i,j,k)-a46(5-l)*var(i,j,k+l)*idz(k)
             enddo
          enddo
       enddo

       k=nz-3
       do j=1,ny
          do i=1,nx
             do l=-7,3
                dvar(i,j,k)=dvar(i,j,k)-a37(4-l)*var(i,j,k+l)*idz(k)
             enddo
          enddo
       enddo

       k=nz-2
       do j=1,ny
          do i=1,nx
             do l=-8,2
                dvar(i,j,k)=dvar(i,j,k)-a28(3-l)*var(i,j,k+l)*idz(k)
             enddo
          enddo
       enddo

       k=nz-1
       do j=1,ny
          do i=1,nx
             do l=-9,1
                dvar(i,j,k)=dvar(i,j,k)-a19(2-l)*var(i,j,k+l)*idz(k)
             enddo
          enddo
       enddo

       k=nz
       do j=1,ny
          do i=1,nx
             do l=-10,0
                dvar(i,j,k)=dvar(i,j,k)-a010(1-l)*var(i,j,k+l)*idz(k)
             enddo
          enddo
       enddo
    endif

    do k=ndz_,nfz_
       do j=1,ny
          do i=1,nx
             dvar(i,j,k)=( a11(1)*(var(i,j,k+1)-var(i,j,k-1)) &
                         + a11(2)*(var(i,j,k+2)-var(i,j,k-2)) &
                         + a11(3)*(var(i,j,k+3)-var(i,j,k-3)) &
                         + a11(4)*(var(i,j,k+4)-var(i,j,k-4)) &
                         + a11(5)*(var(i,j,k+5)-var(i,j,k-5)) )*idz(k)
          enddo
       enddo
    enddo

  end subroutine deriv_z_11pts

  !==============================================================================
  subroutine deriv_x_9pts(var,dvar)
  !==============================================================================
    !> derivatives of a variable along x on 9-point stencil
  !==============================================================================  
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(in) :: var
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(out) :: dvar
    ! ---------------------------------------------------------------------------
    integer :: ndx_,nfx_

    dvar = 0.0_wp
    ndx_=1;nfx_=nx

    if (is_boundary(1,1)) then
       ndx_=5

       i=1
       do k=1,nz
          do j=1,ny
             do l=0,8
                dvar(i,j,k)=dvar(i,j,k)+a08(1+l)*var(i+l,j,k)*idx(i)
             enddo
          enddo
       enddo

       i=2
       do k=1,nz
          do j=1,ny
             do l=-1,7
                dvar(i,j,k)=dvar(i,j,k)+a17(2+l)*var(i+l,j,k)*idx(i)
             enddo
          enddo
       enddo

       i=3
       do k=1,nz
          do j=1,ny
             do l=-2,6
                dvar(i,j,k)=dvar(i,j,k)+a26(3+l)*var(i+l,j,k)*idx(i)
             enddo
          enddo
       enddo

       i=4
       do k=1,nz
          do j=1,ny
             do l=-3,5
                dvar(i,j,k)=dvar(i,j,k)+a35(4+l)*var(i+l,j,k)*idx(i)
             enddo
          enddo
       enddo
    endif

    if (is_boundary(1,2)) then
       nfx_=nx-4

       i=nx-3
       do k=1,nz
          do j=1,ny
             do l=-5,3
                dvar(i,j,k)=dvar(i,j,k)-a35(4-l)*var(i+l,j,k)*idx(i)
             enddo
          enddo
       enddo

       i=nx-2
       do k=1,nz
          do j=1,ny
             do l=-6,2
                dvar(i,j,k)=dvar(i,j,k)-a26(3-l)*var(i+l,j,k)*idx(i)
             enddo
          enddo
       enddo

       i=nx-1
       do k=1,nz
          do j=1,ny
             do l=-7,1
                dvar(i,j,k)=dvar(i,j,k)-a17(2-l)*var(i+l,j,k)*idx(i)
             enddo
          enddo
       enddo

       i=nx
       do k=1,nz
          do j=1,ny
             do l=-8,0
                dvar(i,j,k)=dvar(i,j,k)-a08(1-l)*var(i+l,j,k)*idx(i)
             enddo
          enddo
       enddo
    endif

    do k=1,nz
       do j=1,ny
          do i=ndx_,nfx_
             dvar(i,j,k)=( a9(1)*(var(i+1,j,k)-var(i-1,j,k)) &
                         + a9(2)*(var(i+2,j,k)-var(i-2,j,k)) &
                         + a9(3)*(var(i+3,j,k)-var(i-3,j,k)) &
                         + a9(4)*(var(i+4,j,k)-var(i-4,j,k)) )*idx(i)
          enddo
       enddo
    enddo

  end subroutine deriv_x_9pts

  !==============================================================================
  subroutine deriv_y_9pts(var,dvar)
  !==============================================================================
    !> derivatives of a variable along y on 9-point stencil
  !==============================================================================  
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(in) :: var
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(out) :: dvar
    ! ---------------------------------------------------------------------------
    integer :: ndy_,nfy_

    dvar = 0.0_wp
    ndy_=1;nfy_=ny

    if (is_boundary(2,1)) then
       ndy_=5

       j=1
       do k=1,nz
          do i=1,nx
             do l=0,8
                dvar(i,j,k)=dvar(i,j,k)+a08(1+l)*var(i,j+l,k)*idy(j)
             enddo
          enddo
       enddo

       j=2
       do k=1,nz
          do i=1,nx
             do l=-1,7
                dvar(i,j,k)=dvar(i,j,k)+a17(2+l)*var(i,j+l,k)*idy(j)
             enddo
          enddo
       enddo

       j=3
       do k=1,nz
          do i=1,nx
             do l=-2,6
                dvar(i,j,k)=dvar(i,j,k)+a26(3+l)*var(i,j+l,k)*idy(j)
             enddo
          enddo
       enddo

       j=4
       do k=1,nz
          do i=1,nx
             do l=-3,5
                dvar(i,j,k)=dvar(i,j,k)+a35(4+l)*var(i,j+l,k)*idy(j)
             enddo
          enddo
       enddo
    endif

    if (is_boundary(2,2)) then
       nfy_=ny-4

       j=ny-3
       do k=1,nz
          do i=1,nx
             do l=-5,3
                dvar(i,j,k)=dvar(i,j,k)-a35(4-l)*var(i,j+l,k)*idy(j)
             enddo
          enddo
       enddo

       j=ny-2
       do k=1,nz
          do i=1,nx
             do l=-6,2
                dvar(i,j,k)=dvar(i,j,k)-a26(3-l)*var(i,j+l,k)*idy(j)
             enddo
          enddo
       enddo

       j=ny-1
       do k=1,nz
          do i=1,nx
             do l=-7,1
                dvar(i,j,k)=dvar(i,j,k)-a17(2-l)*var(i,j+l,k)*idy(j)
             enddo
          enddo
       enddo

       j=ny
       do k=1,nz
          do i=1,nx
             do l=-8,0
                dvar(i,j,k)=dvar(i,j,k)-a08(1-l)*var(i,j+l,k)*idy(j)
             enddo
          enddo
       enddo
    endif

    do k=1,nz
       do j=ndy_,nfy_
          do i=1,nx
             dvar(i,j,k)=( a9(1)*(var(i,j+1,k)-var(i,j-1,k)) &
                         + a9(2)*(var(i,j+2,k)-var(i,j-2,k)) &
                         + a9(3)*(var(i,j+3,k)-var(i,j-3,k)) &
                         + a9(4)*(var(i,j+4,k)-var(i,j-4,k)) )*idy(j)
          enddo
       enddo
    enddo

  end subroutine deriv_y_9pts

  !==============================================================================
  subroutine deriv_z_9pts(var,dvar)
  !==============================================================================
    !> derivatives of a variable along z on 9-point stencil
  !==============================================================================  
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(in) :: var
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(out) :: dvar
    ! ---------------------------------------------------------------------------
    integer :: ndz_,nfz_

    dvar = 0.0_wp
    ndz_=1;nfz_=nz

    if (is_boundary(3,1)) then
       ndz_=5

       k=1
       do j=1,ny
          do i=1,nx
             do l=0,8
                dvar(i,j,k)=dvar(i,j,k)+a08(1+l)*var(i,j,k+l)*idz(k)
             enddo
          enddo
       enddo

       k=2
       do j=1,ny
          do i=1,nx
             do l=-1,7
                dvar(i,j,k)=dvar(i,j,k)+a17(2+l)*var(i,j,k+l)*idz(k)
             enddo
          enddo
       enddo

       k=3
       do j=1,ny
          do i=1,nx
             do l=-2,6
                dvar(i,j,k)=dvar(i,j,k)+a26(3+l)*var(i,j,k+l)*idz(k)
             enddo
          enddo
       enddo

       k=4
       do j=1,ny
          do i=1,nx
             do l=-3,5
                dvar(i,j,k)=dvar(i,j,k)+a35(4+l)*var(i,j,k+l)*idz(k)
             enddo
          enddo
       enddo
    endif

    if (is_boundary(3,2)) then
       nfz_=nz-4

       k=nz-3
       do j=1,ny
          do i=1,nx
             do l=-5,3
                dvar(i,j,k)=dvar(i,j,k)-a35(4-l)*var(i,j,k+l)*idz(k)
             enddo
          enddo
       enddo

       k=nz-2
       do j=1,ny
          do i=1,nx
             do l=-6,2
                dvar(i,j,k)=dvar(i,j,k)-a26(3-l)*var(i,j,k+l)*idz(k)
             enddo
          enddo
       enddo

       k=nz-1
       do j=1,ny
          do i=1,nx
             do l=-7,1
                dvar(i,j,k)=dvar(i,j,k)-a17(2-l)*var(i,j,k+l)*idz(k)
             enddo
          enddo
       enddo

       k=nz
       do j=1,ny
          do i=1,nx
             do l=-8,0
                dvar(i,j,k)=dvar(i,j,k)-a08(1-l)*var(i,j,k+l)*idz(k)
             enddo
          enddo
       enddo
    endif

    do k=ndz_,nfz_
       do j=1,ny
          do i=1,nx
             dvar(i,j,k)=( a9(1)*(var(i,j,k+1)-var(i,j,k-1)) &
                         + a9(2)*(var(i,j,k+2)-var(i,j,k-2)) &
                         + a9(3)*(var(i,j,k+3)-var(i,j,k-3)) &
                         + a9(4)*(var(i,j,k+4)-var(i,j,k-4)) )*idz(k)
          enddo
       enddo
    enddo

  end subroutine deriv_z_9pts

  !==============================================================================
  subroutine deriv_x_7pts(var,dvar)
  !==============================================================================
    !> derivatives of a variable along x on 7-point stencil
  !==============================================================================  
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(in) :: var
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(out) :: dvar
    ! ---------------------------------------------------------------------------
    integer :: ndx_,nfx_

    dvar = 0.0_wp
    ndx_=1;nfx_=nx

    if (is_boundary(1,1)) then
       ndx_=4

       i=1
       do k=1,nz
          do j=1,ny
             do l=0,6
                dvar(i,j,k)=dvar(i,j,k)+a06(1+l)*var(i+l,j,k)*idx(i)
             enddo
          enddo
       enddo

       i=2
       do k=1,nz
          do j=1,ny
             do l=-1,5
                dvar(i,j,k)=dvar(i,j,k)+a15(2+l)*var(i+l,j,k)*idx(i)
             enddo
          enddo
       enddo

       i=3
       do k=1,nz
          do j=1,ny
             do l=-2,4
                dvar(i,j,k)=dvar(i,j,k)+a24(3+l)*var(i+l,j,k)*idx(i)
             enddo
          enddo
       enddo
    endif

    if (is_boundary(1,2)) then
       nfx_=nz-3

       i=nx-2
       do k=1,nz
          do j=1,ny
             do l=-4,2
                dvar(i,j,k)=dvar(i,j,k)-a24(3-l)*var(i+l,j,k)*idx(i)
             enddo
          enddo
       enddo

       i=nx-1
       do k=1,nz
          do j=1,ny
             do l=-5,1
                dvar(i,j,k)=dvar(i,j,k)-a15(2-l)*var(i+l,j,k)*idx(i)
             enddo
          enddo
       enddo

       i=nx
       do k=1,nz
          do j=1,ny
             do l=-6,0
                dvar(i,j,k)=dvar(i,j,k)-a06(1-l)*var(i+l,j,k)*idx(i)
             enddo
          enddo
       enddo
    endif

    do k=1,nz
       do j=1,ny
          do i=ndx_,nfx_
             dvar(i,j,k)=( a7(1)*(var(i+1,j,k)-var(i-1,j,k)) &
                         + a7(2)*(var(i+2,j,k)-var(i-2,j,k)) &
                         + a7(3)*(var(i+3,j,k)-var(i-3,j,k)) )*idx(i)
          enddo
       enddo
    enddo

  end subroutine deriv_x_7pts

  !==============================================================================
  subroutine deriv_y_7pts(var,dvar)
  !==============================================================================
    !> derivatives of a variable along y on 7-point stencil
  !==============================================================================  
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(in) :: var
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(out) :: dvar
    ! ---------------------------------------------------------------------------
    integer :: ndy_,nfy_

    dvar = 0.0_wp
    ndy_=1;nfy_=ny

    if (is_boundary(2,1)) then
       ndy_=4

       j=1
       do k=1,nz
          do i=1,nx
             do l=0,6
                dvar(i,j,k)=dvar(i,j,k)+a06(1+l)*var(i,j+l,k)*idy(j)
             enddo
          enddo
       enddo

       j=2
       do k=1,nz
          do i=1,nx
             do l=-1,5
                dvar(i,j,k)=dvar(i,j,k)+a15(2+l)*var(i,j+l,k)*idy(j)
             enddo
          enddo
       enddo

       j=3
       do k=1,nz
          do i=1,nx
             do l=-2,4
                dvar(i,j,k)=dvar(i,j,k)+a24(3+l)*var(i,j+l,k)*idy(j)
             enddo
          enddo
       enddo
    endif

    if (is_boundary(2,2)) then
       nfy_=ny-3

       j=ny-2
       do k=1,nz
          do i=1,nx
             do l=-4,2
                dvar(i,j,k)=dvar(i,j,k)-a24(3-l)*var(i,j+l,k)*idy(j)
             enddo
          enddo
       enddo

       j=ny-1
       do k=1,nz
          do i=1,nx
             do l=-5,1
                dvar(i,j,k)=dvar(i,j,k)-a15(2-l)*var(i,j+l,k)*idy(j)
             enddo
          enddo
       enddo

       j=ny
       do k=1,nz
          do i=1,nx
             do l=-6,0
                dvar(i,j,k)=dvar(i,j,k)-a06(1-l)*var(i,j+l,k)*idy(j)
             enddo
          enddo
       enddo
    endif

    do k=1,nz
       do j=ndy_,nfy_
          do i=1,nx
             dvar(i,j,k)=( a7(1)*(var(i,j+1,k)-var(i,j-1,k)) &
                         + a7(2)*(var(i,j+2,k)-var(i,j-2,k)) &
                         + a7(3)*(var(i,j+3,k)-var(i,j-3,k)) )*idy(j)
          enddo
       enddo
    enddo

  end subroutine deriv_y_7pts

  !==============================================================================
  subroutine deriv_z_7pts(var,dvar)
  !==============================================================================
    !> derivatives of a variable along z on 7-point stencil
  !==============================================================================  
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(in) :: var
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(out) :: dvar
    ! ---------------------------------------------------------------------------
    integer :: ndz_,nfz_

    dvar = 0.0_wp
    ndz_=1;nfz_=nz

    if (is_boundary(3,1)) then
       ndz_=4

       k=1
       do j=1,ny
          do i=1,nx
             do l=0,6
                dvar(i,j,k)=dvar(i,j,k)+a06(1+l)*var(i,j,k+l)*idz(k)
             enddo
          enddo
       enddo

       k=2
       do j=1,ny
          do i=1,nx
             do l=-1,5
                dvar(i,j,k)=dvar(i,j,k)+a15(2+l)*var(i,j,k+l)*idz(k)
             enddo
          enddo
       enddo

       k=3
       do j=1,ny
          do i=1,nx
             do l=-2,4
                dvar(i,j,k)=dvar(i,j,k)+a24(3+l)*var(i,j,k+l)*idz(k)
             enddo
          enddo
       enddo
    endif

    if (is_boundary(3,2)) then
       nfz_=nz-3

       k=nz-2
       do j=1,ny
          do i=1,nx
             do l=-4,2
                dvar(i,j,k)=dvar(i,j,k)-a24(3-l)*var(i,j,k+l)*idz(k)
             enddo
          enddo
       enddo

       k=nz-1
       do j=1,ny
          do i=1,nx
             do l=-5,1
                dvar(i,j,k)=dvar(i,j,k)-a15(2-l)*var(i,j,k+l)*idz(k)
             enddo
          enddo
       enddo

       k=nz
       do j=1,ny
          do i=1,nx
             do l=-6,0
                dvar(i,j,k)=dvar(i,j,k)-a06(1-l)*var(i,j,k+l)*idz(k)
             enddo
          enddo
       enddo
    endif

    do k=ndz_,nfz_
       do j=1,ny
          do i=1,nx
             dvar(i,j,k)=( a7(1)*(var(i,j,k+1)-var(i,j,k-1)) &
                         + a7(2)*(var(i,j,k+2)-var(i,j,k-2)) &
                         + a7(3)*(var(i,j,k+3)-var(i,j,k-3)) )*idz(k)
          enddo
       enddo
    enddo

  end subroutine deriv_z_7pts

  !==============================================================================
  subroutine deriv_x_5pts(var,dvar)
  !==============================================================================
    !> derivatives of a variable along x on 5-point stencil
  !==============================================================================  
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(in) :: var
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(out) :: dvar
    ! ---------------------------------------------------------------------------
    integer :: ndx_,nfx_

    dvar = 0.0_wp
    ndx_=1;nfx_=nx

    if (is_boundary(1,1)) then
       ndx_=3

       i=1
       do k=1,nz
          do j=1,ny
             do l=0,4
                dvar(i,j,k)=dvar(i,j,k)+a04(1+l)*var(i+l,j,k)*idx(i)
             enddo
          enddo
       enddo

       i=2
       do k=1,nz
          do j=1,ny
             do l=-1,3
                dvar(i,j,k)=dvar(i,j,k)+a13(2+l)*var(i+l,j,k)*idx(i)
             enddo
          enddo
       enddo
    endif

    if (is_boundary(1,2)) then
       nfx_=nx-2

       i=nx-1
       do k=1,nz
          do j=1,ny
             do l=-3,1
                dvar(i,j,k)=dvar(i,j,k)-a13(2-l)*var(i+l,j,k)*idx(i)
             enddo
          enddo
       enddo

       i=nx
       do k=1,nz
          do j=1,ny
             do l=-4,0
                dvar(i,j,k)=dvar(i,j,k)-a04(1-l)*var(i+l,j,k)*idx(i)
             enddo
          enddo
       enddo
    endif

    do k=1,nz
       do j=1,ny
          do i=ndx_,nfx_
             dvar(i,j,k)=( a5(1)*(var(i+1,j,k)-var(i-1,j,k)) &
                         + a5(2)*(var(i+2,j,k)-var(i-2,j,k)) )*idx(i)
          enddo
       enddo
    enddo

  end subroutine deriv_x_5pts

  !==============================================================================
  subroutine deriv_y_5pts(var,dvar)
  !==============================================================================
    !> derivatives of a variable along y on 5-point stencil
  !==============================================================================  
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(in) :: var
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(out) :: dvar
    ! ---------------------------------------------------------------------------
    integer :: ndy_,nfy_

    dvar = 0.0_wp
    ndy_=1;nfy_=ny

    if (is_boundary(2,1)) then
       ndy_=3

       j=1
       do k=1,nz
          do i=1,nx
             do l=0,4
                dvar(i,j,k)=dvar(i,j,k)+a04(1+l)*var(i,j+l,k)*idy(j)
             enddo
          enddo
       enddo

       j=2
       do k=1,nz
          do i=1,nx
             do l=-1,3
                dvar(i,j,k)=dvar(i,j,k)+a13(2+l)*var(i,j+l,k)*idy(j)
             enddo
          enddo
       enddo
    endif

    if (is_boundary(2,2)) then
       nfy_=ny-2

       j=ny-1
       do k=1,nz
          do i=1,nx
             do l=-3,1
                dvar(i,j,k)=dvar(i,j,k)-a13(2-l)*var(i,j+l,k)*idy(j)
             enddo
          enddo
       enddo

       j=ny
       do k=1,nz
          do i=1,nx
             do l=-4,0
                dvar(i,j,k)=dvar(i,j,k)-a04(1-l)*var(i,j+l,k)*idy(j)
             enddo
          enddo
       enddo
    endif

    do k=1,nz
       do j=ndy_,nfy_
          do i=1,nx
             dvar(i,j,k)=( a5(1)*(var(i,j+1,k)-var(i,j-1,k)) &
                         + a5(2)*(var(i,j+2,k)-var(i,j-2,k)) )*idy(j)
          enddo
       enddo
    enddo

  end subroutine deriv_y_5pts

  !==============================================================================
  subroutine deriv_z_5pts(var,dvar)
  !==============================================================================
    !> derivatives of a variable along z on 5-point stencil
  !==============================================================================  
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(in) :: var
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(out) :: dvar
    ! ---------------------------------------------------------------------------
    integer :: ndz_,nfz_

    dvar = 0.0_wp
    ndz_=1;nfz_=nz

    if (is_boundary(3,1)) then
       ndz_=3

       k=1
       do j=1,ny
          do i=1,nx
             do l=0,4
                dvar(i,j,k)=dvar(i,j,k)+a04(1+l)*var(i,j,k+l)*idz(k)
             enddo
          enddo
       enddo

       k=2
       do j=1,ny
          do i=1,nx
             do l=-1,3
                dvar(i,j,k)=dvar(i,j,k)+a13(2+l)*var(i,j,k+l)*idz(k)
             enddo
          enddo
       enddo
    endif

    if (is_boundary(3,2)) then
       nfz_=nz-2

       k=nz-1
       do j=1,ny
          do i=1,nx
             do l=-3,1
                dvar(i,j,k)=dvar(i,j,k)-a13(2-l)*var(i,j,k+l)*idz(k)
             enddo
          enddo
       enddo

       k=nz
       do j=1,ny
          do i=1,nx
             do l=-4,0
                dvar(i,j,k)=dvar(i,j,k)-a04(1-l)*var(i,j,k+l)*idz(k)
             enddo
          enddo
       enddo
    endif

    do k=ndz_,nfz_
       do j=1,ny
          do i=1,nx
             dvar(i,j,k)=( a5(1)*(var(i,j,k+1)-var(i,j,k-1)) &
                         + a5(2)*(var(i,j,k+2)-var(i,j,k-2)) )*idz(k)
          enddo
       enddo
    enddo

  end subroutine deriv_z_5pts

  !==============================================================================
  subroutine deriv_x_3pts(var,dvar)
  !==============================================================================
    !> derivatives of a variable along x on 3-point stencil
  !==============================================================================  
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(in) :: var
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(out) :: dvar
    ! ---------------------------------------------------------------------------
    integer :: ndx_,nfx_

    dvar = 0.0_wp
    ndx_=1;nfx_=nx

    if (is_boundary(1,1)) then
       ndx_=2

       i=1
       do k=1,nz
          do j=1,ny
             dvar(i,j,k)=(var(i+1,j,k)-var(i,j,k))*idx(i)
          enddo
       enddo
    endif

    if (is_boundary(1,2)) then
       nfx_=nx-1

       i=nx
       do k=1,nz
          do j=1,ny
             dvar(i,j,k)=(var(i,j,k)-var(i-1,j,k))*idx(i)
          enddo
       enddo
    endif

    do k=1,nz
       do j=1,ny
          do i=ndx_,nfx_
             dvar(i,j,k)=0.5_wp*(var(i+1,j,k)-var(i-1,j,k))*idx(i)
          enddo
       enddo
    enddo

  end subroutine deriv_x_3pts

  !==============================================================================
  subroutine deriv_y_3pts(var,dvar)
  !==============================================================================
    !> derivatives of a variable along y on 3-point stencil
  !==============================================================================  
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(in) :: var
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(out) :: dvar
    ! ---------------------------------------------------------------------------
    integer :: ndy_,nfy_

    dvar = 0.0_wp
    ndy_=1;nfy_=ny

    if (is_boundary(2,1)) then
       ndy_=2

       j=1
       do k=1,nz
          do i=1,nx
             dvar(i,j,k)=(var(i,j+1,k)-var(i,j,k))*idy(j)
          enddo
       enddo
    endif

    if (is_boundary(2,2)) then
       nfy_=ny-1

       j=ny
       do k=1,nz
          do i=1,nx
             dvar(i,j,k)=(var(i,j,k)-var(i,j-1,k))*idy(j)
          enddo
       enddo
    endif

    do k=1,nz
       do j=ndy_,nfy_
          do i=1,nx
             dvar(i,j,k)=0.5_wp*(var(i,j+1,k)-var(i,j-1,k))*idy(j)
          enddo
       enddo
    enddo

  end subroutine deriv_y_3pts

  !==============================================================================
  subroutine deriv_z_3pts(var,dvar)
  !==============================================================================
    !> derivatives of a variable along z on 3-point stencil
  !==============================================================================  
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(in) :: var
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(out) :: dvar
    ! ---------------------------------------------------------------------------
    integer :: ndz_,nfz_

    dvar = 0.0_wp
    ndz_=1;nfz_=nz

    if (is_boundary(3,1)) then
       ndz_=2

       k=1
       do j=1,ny
          do i=1,nx
             dvar(i,j,k)=(var(i,j,k+1)-var(i,j,k))*idz(k)
          enddo
       enddo
    endif

    if (is_boundary(3,2)) then
       nfz_=nz-1

       k=nz
       do j=1,ny
          do i=1,nx
             dvar(i,j,k)=(var(i,j,k)-var(i,j,k-1))*idz(k)
          enddo
       enddo
    endif

    do k=ndz_,nfz_
       do j=1,ny
          do i=1,nx
             dvar(i,j,k)=0.5_wp*(var(i,j,k+1)-var(i,j,k-1))*idz(k)
          enddo
       enddo
    enddo

  end subroutine deriv_z_3pts

end module mod_deriv
