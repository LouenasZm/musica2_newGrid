! ===============================================================================
subroutine grad_scalar_3pts_c3(var,dvarx,dvary,dvarz,dumy)
! ===============================================================================
  !> Compute scalar derivatives
  !> 3-point stencil - curvilinear version -
! ===============================================================================
  use mod_constant
  use mod_coeff_deriv
  use mod_flow
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j,k
  real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(in)  :: var,dumy
  real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v), intent(out) :: dvarx,dvary,dvarz
  real(wp) :: dvardx,dvardy,dvardz
  ! ---------------------------------------------------------------------------

  dvarx = 0.0_wp
  dvary = 0.0_wp
  dvarz = 0.0_wp

  ! Derivatives along ksi (Caution: dvarx is dvar/dksi, etc ...)
  ! =====================
  if (is_boundary(1,1)) then
     i=1
     do k=1,nz
        do j=1,ny
           dvarx(i,j,k) = var(i+1,j,k)-var(i,j,k)
        enddo
     enddo
  endif

  do k=1,nz
     do j=1,ny
        do i=ndx_v3_r,nfx_v3_r
           dvarx(i,j,k) = a3(1)*(var(i+1,j,k)-var(i-1,j,k))
        enddo
     enddo
  enddo

  if (is_boundary(1,2)) then
     i=nx
     do k=1,nz
        do j=1,ny
           dvarx(i,j,k) = var(i,j,k)-var(i-1,j,k)
        enddo
     enddo
  endif

  ! Derivatives along eta (Caution: dvary is dvar/deta, etc ...)
  ! =====================
  if (is_boundary(2,1)) then
     j=1
     do k=1,nz
        do i=1,nx
           dvary(i,j,k) = var(i,j+1,k)-var(i,j,k)
        enddo
     enddo
  endif

  do k=1,nz
     do j=ndy_v3_r,nfy_v3_r
        do i=1,nx
           dvary(i,j,k) = a3(1)*(var(i,j+1,k)-var(i,j-1,k))
        enddo
     enddo
  enddo

  if (is_boundary(2,2)) then
     j=ny
     do k=1,nz
        do i=1,nx
           dvary(i,j,k) = var(i,j,k)-var(i,j-1,k)
        enddo
     enddo
  endif

  ! Derivatives along phi (Caution: dvarz is dvar/dphi, etc ...)
  ! =====================
  if (is_boundary(3,1)) then
     k=1
     do j=1,ny
        do i=1,nx
           dvarz(i,j,k) = var(i,j,k+1)-var(i,j,k)
        enddo
     enddo
  endif

  do k=ndz_v3_r,nfz_v3_r
     do j=1,ny
        do i=1,nx
           dvarz(i,j,k) = a3(1)*(var(i,j,k+1)-var(i,j,k-1))
        enddo
     enddo
  enddo

  if (is_boundary(3,2)) then
     k=nz
     do j=1,ny
        do i=1,nx
           dvarz(i,j,k) = var(i,j,k)-var(i,j,k-1)
        enddo
     enddo
  endif

  ! Compute derivatives in Cartesian coordinates
  ! ============================================
  do k=1,nz
     do j=1,ny
        do i=1,nx
           dvardx=dvarx(i,j,k)*ksi_x_v(i,j,k)+dvary(i,j,k)*eta_x_v(i,j,k)+dvarz(i,j,k)*phi_x_v(i,j,k)
           dvardy=dvarx(i,j,k)*ksi_y_v(i,j,k)+dvary(i,j,k)*eta_y_v(i,j,k)+dvarz(i,j,k)*phi_y_v(i,j,k)
           dvardz=dvarx(i,j,k)*ksi_z_v(i,j,k)+dvary(i,j,k)*eta_z_v(i,j,k)+dvarz(i,j,k)*phi_z_v(i,j,k)
           ! dvar/dx
           dvarx(i,j,k)=dvardx*ijacob3_v(i,j,k)
           ! dvar/dy
           dvary(i,j,k)=dvardy*ijacob3_v(i,j,k)
           ! dvar/dz
           dvarz(i,j,k)=dvardz*ijacob3_v(i,j,k)
        enddo
     enddo
  enddo

end subroutine grad_scalar_3pts_c3
