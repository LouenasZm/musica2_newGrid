!==============================================================================
subroutine flux_euler_3pts_SA_transition_c
!==============================================================================
  !> Derivatives of Eulerian fluxes (inviscid part) - 3-point stencil -
  !> - curvilinear version -
!==============================================================================
  use mod_mpi
  use mod_coeff_deriv
  use mod_flow
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j,k
  real(wp) :: vc,rovc
  ! ---------------------------------------------------------------------------

  ! Initialize flux derivative arrays at wall points
  ! ================================================
  ! Wall BC at imin
  ! ---------------
  if (is_bc_wall(1,1)) then
      Knutil(1,ndy_e_r:nfy_e_r,ndz_e_r:nfz_e_r)=0.0_wp
      Kgamma(1,ndy_e_r:nfy_e_r,ndz_e_r:nfz_e_r)=0.0_wp
      Kre_theta(1,ndy_e_r:nfy_e_r,ndz_e_r:nfz_e_r)=0.0_wp
  endif
  ! Wall BC at imax
  ! ---------------
  if (is_bc_wall(1,2)) then
      Knutil(nx,ndy_e_r:nfy_e_r,ndz_e_r:nfz_e_r)=0.0_wp
      Kgamma(1,ndy_e_r:nfy_e_r,ndz_e_r:nfz_e_r)=0.0_wp
      Kre_theta(1,ndy_e_r:nfy_e_r,ndz_e_r:nfz_e_r)=0.0_wp
  endif
  ! Wall BC at jmin
  ! ---------------
  if (is_bc_wall(2,1)) then
      Knutil(ndx_e_r:nfx_e_r,1,ndz_e_r:nfz_e_r)=0.0_wp
      Kgamma(ndx_e_r:nfx_e_r,1,ndz_e_r:nfz_e_r)=0.0_wp
      Kre_theta(ndx_e_r:nfx_e_r,1,ndz_e_r:nfz_e_r)=0.0_wp
  endif
  ! Wall BC at jmax
  ! ---------------
  if (is_bc_wall(2,2)) then
      Knutil(ndx_e_r:nfx_e_r,ny,ndz_e_r:nfz_e_r)=0.0_wp
      Kgamma(ndx_e_r:nfx_e_r,1,ndz_e_r:nfz_e_r)=0.0_wp
      Kre_theta(ndx_e_r:nfx_e_r,1,ndz_e_r:nfz_e_r)=0.0_wp
  endif
  ! Wall BC at kmin
  ! ---------------
  if (is_bc_wall(3,1)) then
      Knutil(ndx_e_r:nfx_e_r,ndy_e_r:nfy_e_r,1)=0.0_wp
      Kgamma(ndx_e_r:nfx_e_r,1,ndz_e_r:nfz_e_r)=0.0_wp
      Kre_theta(ndx_e_r:nfx_e_r,1,ndz_e_r:nfz_e_r)=0.0_wp
  endif
  ! Wall BC at kmax
  ! ---------------
  if (is_bc_wall(3,2)) then
      Knutil(ndx_e_r:nfx_e_r,ndy_e_r:nfy_e_r,nz)=0.0_wp
      Kgamma(ndx_e_r:nfx_e_r,1,ndz_e_r:nfz_e_r)=0.0_wp
      Kre_theta(ndx_e_r:nfx_e_r,1,ndz_e_r:nfz_e_r)=0.0_wp
  endif

  ! Along ksi
  ! ==========================
  ! Interior points
  ! ---------------
  do k=ndz_e_r,nfz_e_r
     do j=ndy_e_r,nfy_e_r
        do i=ndx_r,nfx_r
           ! contravariant velocity
           vc=uu(i,j,k)*y_eta_v(i,j)-vv(i,j,k)*x_eta_v(i,j)
           Knutil(i,j,k) = 0.5_wp*(nutil_n(i+1,j,k)-nutil_n(i-1,j,k))*vc
           Kgamma(i,j,k) = 0.5_wp*(rhogamma_n(i+1,j,k)-rhogamma_n(i-1,j,k))*vc
           Knutil(i,j,k) = 0.5_wp*(nutil_n(i+1,j,k)-nutil_n(i-1,j,k))*vc
        enddo
     enddo
  enddo
 
  ! Along eta
  ! ==========================
  ! Interior points
  ! ---------------
  do k=ndz_e_r,nfz_e_r
     do i=ndx_e_r,nfx_e_r
        do j=ndy_r,nfy_r
           ! contravariant velocity
           vc=vv(i,j,k)*x_ksi_v(i,j)-uu(i,j,k)*y_ksi_v(i,j)
           Knutil(i,j,k) = 0.5_wp*(nutil_n(i,j+1,k)-nutil_n(i,j-1,k))*vc + Knutil(i,j,k)
           Kgamma(i,j,k) = 0.5_wp*(rhogamma_n(i,j+1,k)-rhogamma_n(i,j-1,k))*vc + Kgamma(i,j,k)
           Kre_theta(i,j,k) = 0.5_wp*(rhore_theta_n(i,j+1,k)-rhore_theta_n(i,j-1,k))*vc + Kre_theta(i,j,k)
        enddo
     enddo
  enddo

  ! Multiply by inverse Jacobian
  ! ============================
  do k=ndz_e_r,nfz_e_r
     do j=ndy_e_r,nfy_e_r
        do i=ndx_e_r,nfx_e_r
            Knutil(i,j,k)=  Knutil(i,j,k)*ijacob_v(i,j)
            Kgamma(i,j,k)=  Kgamma(i,j,k)*ijacob_v(i,j)
            Kre_theta(i,j,k)=  Kre_theta(i,j,k)*ijacob_v(i,j)
        enddo
     enddo
  enddo

  !****************
  if (is_2D) return
  !****************

  ! Flux derivatives along z
  ! ========================
  ! Interior points
  ! ---------------
  do k=ndz_r,nfz_r
     do j=ndy_e_r,nfy_e_r
        do i=ndx_e_r,nfx_e_r
           Knutil(i,j,k) = 0.5_wp*(nutil_n(i,j,k+1)-nutil_n(i,j,k-1))*ww(i,j,k) + Knutil(i,j,k)
           Kgamma(i,j,k) = 0.5_wp*(rhogamma_n(i,j,k+1)-rhogamma_n(i,j,k-1))*ww(i,j,k) + Kgamma(i,j,k)
           Kre_theta(i,j,k) = 0.5_wp*(rhore_theta_n(i,j,k+1)-rhore_theta_n(i,j,k-1))*ww(i,j,k) + Kre_theta(i,j,k)
        enddo
     enddo
  enddo

end subroutine flux_euler_3pts_SA_transition_c
