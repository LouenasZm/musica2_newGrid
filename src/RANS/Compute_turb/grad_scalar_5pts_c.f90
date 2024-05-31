! ===============================================================================
subroutine grad_scalar_5pts_c(var,dvarx,dvary,dvarz)
! ===============================================================================
  !> Compute varde derivatives
  !> 5-point stencil - curvilinear version -
! ===============================================================================
  use mod_constant
  use mod_coeff_deriv
  use mod_flow
  use mod_mpi
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j,k
  real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(in)  :: var
  real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v), intent(out) :: dvarx,dvary,dvarz
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
           dvarx(i,j,k) = (a04(1)*var(i  ,j,k)+a04(2)*var(i+1,j,k) &
                          +a04(3)*var(i+2,j,k)+a04(4)*var(i+3,j,k) &
                          +a04(5)*var(i+4,j,k))
        enddo
     enddo

     i=2
     do k=1,nz
        do j=1,ny
           dvarx(i,j,k) = (a13(1)*var(i-1,j,k)+a13(2)*var(i,  j,k) &
                          +a13(3)*var(i+1,j,k)+a13(4)*var(i+2,j,k) &
                          +a13(5)*var(i+3,j,k))
        enddo
     enddo
  endif

  do k=1,nz
     do j=1,ny
        do i=ndx_v3_r,nfx_v3_r
           dvarx(i,j,k) = a5(1)*(var(i+1,j,k)-var(i-1,j,k)) &
                         +a5(2)*(var(i+2,j,k)-var(i-2,j,k))
        enddo
     enddo
  enddo

  if (is_boundary(1,2)) then
     i=nx-1
     do k=1,nz
        do j=1,ny
           dvarx(i,j,k) = (a31(5)*var(i-3,j,k)+a31(4)*var(i-2,j,k) &
                          +a31(3)*var(i-1,j,k)+a31(2)*var(i  ,j,k) &
                          +a31(1)*var(i+1,j,k))
        enddo
     enddo

     i=nx
     do k=1,nz
        do j=1,ny
           dvarx(i,j,k) = (a40(5)*var(i-4,j,k)+a40(4)*var(i-3,j,k) &
                          +a40(3)*var(i-2,j,k)+a40(2)*var(i-1,j,k) &
                          +a40(1)*var(i  ,j,k) )
        enddo
     enddo
  endif

  ! Derivatives along eta (Caution: dvary is dvar/deta, etc ...)
  ! =====================
  if (is_boundary(2,1)) then
     j=1
     do k=1,nz
        do i=1,nx
 
           dvary(i,j,k) = (a04(1)*var(i,j  ,k)+a04(2)*var(i,j+1,k) &
                          +a04(3)*var(i,j+2,k)+a04(4)*var(i,j+3,k) &
                          +a04(5)*var(i,j+4,k))
        enddo
     enddo

     j=2
     do k=1,nz
        do i=1,nx
           dvary(i,j,k) = (a13(1)*var(i,j-1,k)+a13(2)*var(i,j  ,k) &
                          +a13(3)*var(i,j+1,k)+a13(4)*var(i,j+2,k) &
                          +a13(5)*var(i,j+3,k))
        enddo
     enddo
  endif

  do k=1,nz
     do j=ndy_v3_r,nfy_v3_r
        do i=1,nx
           dvary(i,j,k) = a5(1)*(var(i,j+1,k)-var(i,j-1,k)) &
                         +a5(2)*(var(i,j+2,k)-var(i,j-2,k))
        enddo
     enddo
  enddo

  if (is_boundary(2,2)) then
     j=ny-1
     do k=1,nz
        do i=1,nx
           dvary(i,j,k) = (a31(5)*var(i,j-3,k)+a31(4)*var(i,j-2,k) &
                          +a31(3)*var(i,j-1,k)+a31(2)*var(i,j  ,k) &
                          +a31(1)*var(i,j+1,k))
        enddo
     enddo

     j=ny
     do k=1,nz
        do i=1,nx
           dvary(i,j,k) = (a40(5)*var(i,j-4,k)+a40(4)*var(i,j-3,k) &
                          +a40(3)*var(i,j-2,k)+a40(2)*var(i,j-1,k) &
                          +a40(1)*var(i,j  ,k) )
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
           dvarz(i,j,k) = (a04(1)*var(i,j,k  )+a04(2)*var(i,j,k+1) &
                          +a04(3)*var(i,j,k+2)+a04(4)*var(i,j,k+3) &
                          +a04(5)*var(i,j,k+4))*idz_v(k)
        enddo
     enddo
    
     k=2
     do j=1,ny
        do i=1,nx
           dvarz(i,j,k) = (a13(1)*var(i,j,k-1)+a13(2)*var(i,j,k  ) &
                          +a13(3)*var(i,j,k+1)+a13(4)*var(i,j,k+2) &
                          +a13(5)*var(i,j,k+3))*idz_v(k)
        enddo
     enddo
  endif

  do k=ndz_v3_r,nfz_v3_r
     do j=1,ny
        do i=1,nx
           dvarz(i,j,k) = (a5(1)*( var(i,j,k+1)-var(i,j,k-1)) &
                         + a5(2)*( var(i,j,k+2)-var(i,j,k-2)))*idz_v(k)
        enddo
     enddo
  enddo

  if (is_boundary(3,2)) then
     k=nz-1
     do j=1,ny
        do i=1,nx
           dvarz(i,j,k) = (a31(5)*var(i,j,k-3)+a31(4)*var(i,j,k-2) &
                          +a31(3)*var(i,j,k-1)+a31(2)*var(i,j,k  ) &
                          +a31(1)*var(i,j,k+1))*idz_v(k)
        enddo
     enddo

     k=nz
     do j=1,ny
        do i=1,nx
           dvarz(i,j,k) = (a40(5)*var(i,j,k-4)+a40(4)*var(i,j,k-3) &
                            +a40(3)*var(i,j,k-2)+a40(2)*var(i,j,k-1) &
                            +a40(1)*var(i,j,k))*idz_v(k)
        enddo
     enddo
  endif

end subroutine grad_scalar_5pts_c
