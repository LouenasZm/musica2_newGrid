! ===============================================================================
subroutine grad_scalar_3pts_c(var,dvarx,dvary,dvarz,dumy)
! ===============================================================================
  !> Compute scalar derivatives
  !> 3-point stencil - curvilinear version -
! ===============================================================================
  use mod_mpi
  use mod_constant
  use mod_coeff_deriv
  use mod_flow
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j,k
  real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(in)  :: var,dumy
  real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(out) :: dvarx,dvary,dvarz
  real(wp) :: dvarksi_y_eta,dvareta_y_ksi,dvareta_x_ksi,dvarksi_x_eta
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
           dvarx(i,j,k)= a02(1)*var(i  ,j,k)+a02(2)*var(i+1,j,k) &
                        +a02(3)*var(i+2,j,k)
        enddo
     enddo
  endif

  do k=1,nz
     do j=1,ny
        do i=ndx_v3,nfx_v3
           dvarx(i,j,k) = a3(1)*(var(i+1,j,k)-var(i-1,j,k))
        enddo
     enddo
  enddo

  if (is_boundary(1,2)) then
     i=nx
     do k=1,nz
        do j=1,ny
           dvarx(i,j,k) = a20(3)*var(i-2,j,k)+a20(2)*var(i-1,j,k) &
                         +a20(1)*var(i  ,j,k)
        enddo
     enddo
  endif

  ! Derivatives along eta (Caution: dvary is dvar/deta, etc ...)
  ! =====================
  if (is_boundary(2,1)) then
     j=1
     do k=1,nz
        do i=1,nx
           dvary(i,j,k) = a02(1)*var(i,j  ,k)+a02(2)*var(i,j+1,k) &
                         +a02(3)*var(i,j+2,k)
        enddo
     enddo
  endif

  do k=1,nz
     do j=ndy_v3,nfy_v3
        do i=1,nx
           dvary(i,j,k) = a3(1)*(var(i,j+1,k)-var(i,j-1,k))
        enddo
     enddo
  enddo

  if (is_boundary(2,2)) then
     j=ny
     do k=1,nz
        do i=1,nx
           dvary(i,j,k) = a20(3)*var(i,j-2,k)+a20(2)*var(i,j-1,k) &
                         +a20(1)*var(i,j  ,k)
        enddo
     enddo
  endif

  ! Compute derivatives in Cartesian coordinates
  ! ============================================
  do k=1,nz
     do j=1,ny
        do i=1,nx
           dvarksi_y_eta=dvarx(i,j,k)*y_eta_v(i,j)
           dvareta_y_ksi=dvary(i,j,k)*y_ksi_v(i,j)
           dvareta_x_ksi=dvary(i,j,k)*x_ksi_v(i,j)
           dvarksi_x_eta=dvarx(i,j,k)*x_eta_v(i,j)
           ! dvar/dx
           dvarx(i,j,k)=(dvarksi_y_eta-dvareta_y_ksi)*ijacob_v(i,j)
           ! dvar/dy
           dvary(i,j,k)=(dvareta_x_ksi-dvarksi_x_eta)*ijacob_v(i,j)
        enddo
     enddo
  enddo

  !****************
  if (is_2D) return
  !****************

  ! Derivatives along z
  ! ===================
  if (is_boundary(3,1)) then
     k=1
     do j=1,ny
        do i=1,nx
           dvarz(i,j,k) = ( a02(1)*var(i,j,k  )+a02(2)*var(i,j,k+1) &
                           +a02(3)*var(i,j,k+2))*idz_v(k)
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
           dvarz(i,j,k) = ( a20(3)*var(i,j,k-2)+a20(2)*var(i,j,k-1) &
                           +a20(1)*var(i,j,k  ) )*idz_v(k)
        enddo
     enddo
  endif

end subroutine grad_scalar_3pts_c
