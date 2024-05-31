! ==============================================================================
module mod_deriv_c
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
  subroutine deriv_ksi_5pts_c(var,dvar,nx_1,nx_2,ny_1,ny_2,nz_1,nz_2)
  !==============================================================================
    !> derivatives of a variable along ksi on 5-point stencil
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer, intent(in) :: nx_1,nx_2,ny_1,ny_2,nz_1,nz_2
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(in) :: var
    real(wp), dimension(nx_1:nx_2,ny_1:ny_2,nz_1:nz_2), intent(out) :: dvar
    ! ---------------------------------------------------------------------------
    integer :: ndx_,nfx_
    ! ---------------------------------------------------------------------------

    dvar = 0.0_wp
    ndx_=nx_1;nfx_=nx_2

    if ((is_boundary(1,1)).and.(nx_1.ne.nx_2)) then
       ndx_=3

       i=1
       do k=nz_1,nz_2
          do j=ny_1,ny_2
             do l=0,4
                dvar(i,j,k)=dvar(i,j,k)+a04(1+l)*var(i+l,j,k)
             enddo
          enddo
       enddo

       i=2
       do k=nz_1,nz_2
          do j=ny_1,ny_2
             do l=-1,3
                dvar(i,j,k)=dvar(i,j,k)+a13(2+l)*var(i+l,j,k)
             enddo
          enddo
       enddo
    endif

    if ((is_boundary(1,2)).and.(nx_1.ne.nx_2))then
       nfx_=nx-2

       i=nx-1
       do k=nz_1,nz_2
          do j=ny_1,ny_2
             do l=-3,1
                dvar(i,j,k)=dvar(i,j,k)-a13(2-l)*var(i+l,j,k)
             enddo
          enddo
       enddo

       i=nx
       do k=nz_1,nz_2
          do j=ny_1,ny_2
             do l=-4,0
                dvar(i,j,k)=dvar(i,j,k)-a04(1-l)*var(i+l,j,k)
             enddo
          enddo
       enddo
    endif

    do k=nz_1,nz_2
       do j=ny_1,ny_2
          do i=ndx_,nfx_
             dvar(i,j,k)=( a5(1)*(var(i+1,j,k)-var(i-1,j,k)) &
                         + a5(2)*(var(i+2,j,k)-var(i-2,j,k)) )
          enddo
       enddo
    enddo

  end subroutine deriv_ksi_5pts_c

  !==============================================================================
  subroutine deriv_eta_5pts_c(var,dvar,nx_1,nx_2,ny_1,ny_2,nz_1,nz_2)
  !==============================================================================
    !> derivatives of a variable along eta on 5-point stencil
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer, intent(in) :: nx_1,nx_2,ny_1,ny_2,nz_1,nz_2
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(in) :: var
    real(wp), dimension(nx_1:nx_2,ny_1:ny_2,nz_1:nz_2), intent(out) :: dvar
    ! ---------------------------------------------------------------------------
    integer :: ndy_,nfy_
    ! ---------------------------------------------------------------------------

    dvar = 0.0_wp
    ndy_=ny_1;nfy_=ny_2

    if ((is_boundary(2,1)).and.(ny_1.ne.ny_2)) then
       ndy_=3

       j=1
       do k=nz_1,nz_2
          do i=nx_1,nx_2
             do l=0,4
                dvar(i,j,k)=dvar(i,j,k)+a04(1+l)*var(i,j+l,k)
             enddo
          enddo
       enddo

       j=2
       do k=nz_1,nz_2
          do i=nx_1,nx_2
             do l=-1,3
                dvar(i,j,k)=dvar(i,j,k)+a13(2+l)*var(i,j+l,k)
             enddo
          enddo
       enddo
    endif

    if ((is_boundary(2,2)).and.(ny_1.ne.ny_2)) then
       nfy_=ny-2

       j=ny-1
       do k=nz_1,nz_2
          do i=nx_1,nx_2
             do l=-3,1
                dvar(i,j,k)=dvar(i,j,k)-a13(2-l)*var(i,j+l,k)
             enddo
          enddo
       enddo

       j=ny
       do k=nz_1,nz_2
          do i=nx_1,nx_2
             do l=-4,0
                dvar(i,j,k)=dvar(i,j,k)-a04(1-l)*var(i,j+l,k)
             enddo
          enddo
       enddo
    endif

    do k=nz_1,nz_2
       do j=ndy_,nfy_
          do i=nx_1,nx_2
             dvar(i,j,k)=( a5(1)*(var(i,j+1,k)-var(i,j-1,k)) &
                         + a5(2)*(var(i,j+2,k)-var(i,j-2,k)) )
          enddo
       enddo
    enddo

  end subroutine deriv_eta_5pts_c

  !==============================================================================
  subroutine deriv_x_5pts_c(var,dvar,nx_1,nx_2,ny_1,ny_2,nz_1,nz_2)
  !==============================================================================
    !> derivatives of a variable along x on 5-point stencil
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer, intent(in) :: nx_1,nx_2,ny_1,ny_2,nz_1,nz_2
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(in) :: var
    real(wp), dimension(nx_1:nx_2,ny_1:ny_2,nz_1:nz_2), intent(out) :: dvar
    ! ---------------------------------------------------------------------------
    real(wp) :: dvarksi_y_eta,dvareta_y_ksi
    real(wp), dimension(nx_1:nx_2,ny_1:ny_2,nz_1:nz_2) :: dvarksi,dvareta
    ! ---------------------------------------------------------------------------

    call deriv_ksi_5pts_c(var,dvarksi,nx_1,nx_2,ny_1,ny_2,nz_1,nz_2)
    call deriv_eta_5pts_c(var,dvareta,nx_1,nx_2,ny_1,ny_2,nz_1,nz_2)

    ! Compute derivatives in Cartesian coordinates
    ! ============================================
    do k=nz_1,nz_2
       do j=ny_1,ny_2
          do i=nx_1,nx_2
             dvarksi_y_eta=dvarksi(i,j,k)*y_eta_v(i,j)
             dvareta_y_ksi=dvareta(i,j,k)*y_ksi_v(i,j)
             ! d(var)/dx
             dvar(i,j,k)=(dvarksi_y_eta-dvareta_y_ksi)*ijacob_v(i,j)
          enddo
       enddo
    enddo

  end subroutine deriv_x_5pts_c

  !==============================================================================
  subroutine deriv_y_5pts_c(var,dvar,nx_1,nx_2,ny_1,ny_2,nz_1,nz_2)
  !==============================================================================
    !> derivatives of a variable along y on 5-point stencil
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer, intent(in) :: nx_1,nx_2,ny_1,ny_2,nz_1,nz_2
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(in) :: var
    real(wp), dimension(nx_1:nx_2,ny_1:ny_2,nz_1:nz_2), intent(out) :: dvar
    ! ---------------------------------------------------------------------------
    real(wp) :: dvareta_x_ksi,dvarksi_x_eta
    real(wp), dimension(nx_1:nx_2,ny_1:ny_2,nz_1:nz_2) :: dvarksi,dvareta
    ! ---------------------------------------------------------------------------

    call deriv_ksi_5pts_c(var,dvarksi,nx_1,nx_2,ny_1,ny_2,nz_1,nz_2)
    call deriv_eta_5pts_c(var,dvareta,nx_1,nx_2,ny_1,ny_2,nz_1,nz_2)

    ! Compute derivatives in Cartesian coordinates
    ! ============================================
    do k=nz_1,nz_2
       do j=ny_1,ny_2
          do i=nx_1,nx_2
             dvareta_x_ksi=dvareta(i,j,k)*x_ksi_v(i,j)
             dvarksi_x_eta=dvarksi(i,j,k)*x_eta_v(i,j)
             ! d(var)/dy
             dvar(i,j,k)=(dvareta_x_ksi-dvarksi_x_eta)*ijacob_v(i,j)
          enddo
       enddo
    enddo

  end subroutine deriv_y_5pts_c

  !==============================================================================
  subroutine deriv_x_mz_5pts_c(var,dvar,nx_1,nx_2,ny_1,ny_2,nz_1,nz_2)
  !==============================================================================
    !> derivatives of a variable along x on 5-point stencil
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer, intent(in) :: nx_1,nx_2,ny_1,ny_2,nz_1,nz_2
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(in) :: var
    real(wp), dimension(nx_1:nx_2,ny_1:ny_2), intent(out) :: dvar
    ! ---------------------------------------------------------------------------
    real(wp) :: dvarksi_y_eta,dvareta_y_ksi
    real(wp), dimension(nx_1:nx_2,ny_1:ny_2,nz_1:nz_2) :: dvarksi,dvareta
    ! ---------------------------------------------------------------------------

    call deriv_ksi_5pts_c(var,dvarksi,nx_1,nx_2,ny_1,ny_2,nz_1,nz_2)
    call deriv_eta_5pts_c(var,dvareta,nx_1,nx_2,ny_1,ny_2,nz_1,nz_2)

    ! Compute derivatives in Cartesian coordinates
    ! ============================================
    dvar=0.0_wp
    do k=nz_1,nz_2
       do j=ny_1,ny_2
          do i=nx_1,nx_2
             dvarksi_y_eta=dvarksi(i,j,k)*y_eta_v(i,j)
             dvareta_y_ksi=dvareta(i,j,k)*y_ksi_v(i,j)
             ! d(var)/dx
             dvar(i,j)=dvar(i,j)+(dvarksi_y_eta-dvareta_y_ksi)*ijacob_v(i,j)
          enddo
       enddo
    enddo
    dvar=dvar/(nz_2-nz_1+1)

  end subroutine deriv_x_mz_5pts_c

  !==============================================================================
  subroutine deriv_y_mz_5pts_c(var,dvar,nx_1,nx_2,ny_1,ny_2,nz_1,nz_2)
  !==============================================================================
    !> derivatives of a variable along y on 5-point stencil
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer, intent(in) :: nx_1,nx_2,ny_1,ny_2,nz_1,nz_2
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(in) :: var
    real(wp), dimension(nx_1:nx_2,ny_1:ny_2), intent(out) :: dvar
    ! ---------------------------------------------------------------------------
    real(wp) :: dvareta_x_ksi,dvarksi_x_eta
    real(wp), dimension(nx_1:nx_2,ny_1:ny_2,nz_1:nz_2) :: dvarksi,dvareta
    ! ---------------------------------------------------------------------------

    call deriv_ksi_5pts_c(var,dvarksi,nx_1,nx_2,ny_1,ny_2,nz_1,nz_2)
    call deriv_eta_5pts_c(var,dvareta,nx_1,nx_2,ny_1,ny_2,nz_1,nz_2)

    ! Compute derivatives in Cartesian coordinates
    ! ============================================
    dvar=0.0_wp
    do k=nz_1,nz_2
       do j=ny_1,ny_2
          do i=nx_1,nx_2
             dvareta_x_ksi=dvareta(i,j,k)*x_ksi_v(i,j)
             dvarksi_x_eta=dvarksi(i,j,k)*x_eta_v(i,j)
             ! d(var)/dy
             dvar(i,j)=dvar(i,j)+(dvareta_x_ksi-dvarksi_x_eta)*ijacob_v(i,j)
          enddo
       enddo
    enddo
    dvar=dvar/(nz_2-nz_1+1)

  end subroutine deriv_y_mz_5pts_c

end module mod_deriv_c
