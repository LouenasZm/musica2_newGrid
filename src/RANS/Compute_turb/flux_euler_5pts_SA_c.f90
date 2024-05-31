!==============================================================================
subroutine flux_euler_5pts_SA_c
!==============================================================================
  !> Derivatives of Eulerian fluxes (inviscid part) - 5-point stencil -
  !> - curvilinear version -
!==============================================================================
  use mod_mpi
  use mod_coeff_deriv
  use mod_flow
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j,k
  real(wp) :: vc
  ! ---------------------------------------------------------------------------

  ! ==================================== !
  ! /!\ stencil Euler = stencil visc /!\ !
  ! ==================================== !

  ! Initialize flux derivative arrays at wall points
  ! ================================================
  ! Wall BC at imin
  ! ---------------
  if (is_bc_wall(1,1)) then
      Knutil(1,ndy_e_r:nfy_e_r,ndz_e_r:nfz_e_r)=0.0_wp
  endif
  ! Wall BC at imax
  ! ---------------
  if (is_bc_wall(1,2)) then
      Knutil(nx,ndy_e_r:nfy_e_r,ndz_e_r:nfz_e_r)=0.0_wp
  endif
  ! Wall BC at jmin
  ! ---------------
  if (is_bc_wall(2,1)) then
      Knutil(ndx_e_r:nfx_e_r,1,ndz_e_r:nfz_e_r)=0.0_wp
  endif
  ! Wall BC at jmax
  ! ---------------
  if (is_bc_wall(2,2)) then
      Knutil(ndx_e_r:nfx_e_r,ny,ndz_e_r:nfz_e_r)=0.0_wp
  endif
  ! Wall BC at kmin
  ! ---------------
  if (is_bc_wall(3,1)) then
      Knutil(ndx_e_r:nfx_e_r,ndy_e_r:nfy_e_r,1)=0.0_wp
  endif
  ! Wall BC at kmax
  ! ---------------
  if (is_bc_wall(3,2)) then
      Knutil(ndx_e_r:nfx_e_r,ndy_e_r:nfy_e_r,nz)=0.0_wp
  endif

  ! Along ksi
  ! ==========================
  ! Wall BC at imin
  ! ---------------
  if (is_bc_1pt(1,1)) then
     i=2
     do k=ndz_e_r,nfz_e_r
        do j=ndy_e_r,nfy_e_r
           ! contravariant velocity
           vc=uu(i,j,k)*y_eta_v(i,j)-vv(i,j,k)*x_eta_v(i,j)
           Knutil(i,j,k) = 0.5_wp*(nutil_n(i+1,j,k)-nutil_n(i-1,j,k))*vc
        enddo
     enddo
  endif

  ! Interior points
  ! ---------------
  do k=ndz_e_r,nfz_e_r
     do j=ndy_e_r,nfy_e_r
        do i=ndx_r,nfx_r
           ! contravariant velocity
           vc=uu(i,j,k)*y_eta_v(i,j)-vv(i,j,k)*x_eta_v(i,j)
           Knutil(i,j,k) = (a5(1)*(nutil_n(i+1,j,k)-nutil_n(i-1,j,k)) &
                           +a5(2)*(nutil_n(i+2,j,k)-nutil_n(i-2,j,k)))*vc
        enddo
     enddo
  enddo
 
  ! Wall BC at imax
  ! ---------------
  if (is_bc_1pt(1,2)) then
     i=nx-1
     do k=ndz_e_r,nfz_e_r
        do j=ndy_e_r,nfy_e_r
           ! contravariant velocity
           vc=uu(i,j,k)*y_eta_v(i,j)-vv(i,j,k)*x_eta_v(i,j)
           Knutil(i,j,k) = 0.5_wp*(nutil_n(i+1,j,k)-nutil_n(i-1,j,k))*vc
        enddo
     enddo
  endif

  ! Along eta
  ! ==========================
  ! Wall BC at jmin
  ! ---------------
  if (is_bc_1pt(2,1)) then
     j=2
     do k=ndz_e_r,nfz_e_r  
        do i=ndx_e_r,nfx_e_r
           ! contravariant velocity
           vc=vv(i,j,k)*x_ksi_v(i,j)-uu(i,j,k)*y_ksi_v(i,j)
           Knutil(i,j,k) = 0.5_wp*(nutil_n(i,j+1,k)-nutil_n(i,j-1,k))*vc + Knutil(i,j,k)
        enddo
     enddo
  endif

  ! Interior points
  ! ---------------
  do k=ndz_e_r,nfz_e_r
     do j=ndy_r,nfy_r
        do i=ndx_e_r,nfx_e_r
           ! contravariant velocity
           vc=vv(i,j,k)*x_ksi_v(i,j)-uu(i,j,k)*y_ksi_v(i,j)
           Knutil(i,j,k) = (a5(1)*(nutil_n(i,j+1,k)-nutil_n(i,j-1,k)) &
                           +a5(2)*(nutil_n(i,j+2,k)-nutil_n(i,j-2,k)))*vc + Knutil(i,j,k)
        enddo
     enddo
  enddo

  ! Wall BC at jmax
  ! ---------------
  if (is_bc_1pt(2,2)) then
     j=ny-1
     do k=ndz_e_r,nfz_e_r
        do i=ndx_e_r,nfx_e_r
           ! contravariant velocity
           vc=vv(i,j,k)*x_ksi_v(i,j)-uu(i,j,k)*y_ksi_v(i,j)
           Knutil(i,j,k) = 0.5_wp*(nutil_n(i,j+1,k)-nutil_n(i,j-1,k))*vc + Knutil(i,j,k)
        enddo
     enddo
  endif

  ! Multiply by inverse Jacobian
  ! ============================
  do k=ndz_e_r,nfz_e_r
     do j=ndy_e_r,nfy_e_r
        do i=ndx_e_r,nfx_e_r
            Knutil(i,j,k)=  Knutil(i,j,k)*ijacob_v(i,j)
        enddo
     enddo
  enddo

  !****************
  if (is_2D) return
  !****************

  ! Flux derivatives along z
  ! ========================
  
  ! BC at kmin
  ! ----------
  if (is_bc_1pt(3,1)) then
     k=2
     do j=ndy_e_r,nfy_e_r
        do i=ndx_e_r,nfx_e_r
           Knutil(i,j,k) = Knutil(i,j,k) + (nutil_n(i,j,k+1)-nutil_n(i,j,k-1))*idz2_kmin*ww(i,j,k)
        enddo
     enddo
  endif

  ! Interior points
  ! ---------------
  do k=ndz_r,nfz_r
     do j=ndy_e_r,nfy_e_r
        do i=ndx_e_r,nfx_e_r
           Knutil(i,j,k) = (a5(1)*(nutil_n(i,j,k+1)-nutil_n(i,j,k-1)) &
                           +a5(2)*(nutil_n(i,j,k+2)-nutil_n(i,j,k-2)))*idz4_kmin*ww(i,j,k) + Knutil(i,j,k)
        enddo
     enddo
  enddo
 
  ! BC at kmax
  ! ----------
  if (is_bc_1pt(3,2)) then
     k=nz-1
     do j=ndy_e_r,nfy_e_r
        do i=ndx_e_r,nfx_e_r
           Knutil(i,j,k) = Knutil(i,j,k) + (nutil_n(i,j,k+1)-nutil_n(i,j,k-1))*idz2_kmax*ww(i,j,k)
        enddo
     enddo
  endif

end subroutine flux_euler_5pts_SA_c
