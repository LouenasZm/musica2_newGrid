!===============================================================================
subroutine flux_euler_w_imin_jmin_5pts_SA_c
!===============================================================================
  !> Derivatives of Euler fluxes for edge I_MIN/J_MIN close to two walls
  !> - curvilinear coordinate - 5-point stencil -
!===============================================================================
  use mod_flow
  use mod_coeff_deriv
  implicit none
  ! ----------------------------------------------------------------------------
  integer :: i,j,k
  real(wp) :: vc,rovc
  ! ----------------------------------------------------------------------------

  ! ~> First pass: 5-pt metrics & order reduction for derivatives
  ! =============================================================
 
  ! fluxes along ksi
  ! ----------------
  do k=ndz_e,nfz_e
     do j=1,2
        do i=-1,4
           ! contravariant velocity
           vc=uu(i,j,k)*y_eta(i,j)-vv(i,j,k)*x_eta(i,j)
           rovc=rho_n(i,j,k)*vc

           Fnutil(i,j,k) = rovc*nutil(i,j,k)
        enddo
     enddo
  enddo
  
  ! derivatives along ksi
  ! ---------------------
  i=1
  do j=1,2
     do k=ndz_e,nfz_e
        Knutil(i,j,k) = (Fnutil(i+1,j,k) - Fnutil(i,j,k))
     enddo
  enddo

  i=2
  do j=1,2
     do k=ndz_e,nfz_e
        Knutil(i,j,k) = 0.5_wp*(Fnutil(i+1,j,k) - Fnutil(i-1,j,k))
     enddo
  enddo

  ! modified fluxes along eta
  ! -------------------------
  do k=ndz_e,nfz_e
     do i=1,2
        do j=-1,4
           ! contravariant velocity
           vc=vv(i,j,k)*x_ksi(i,j)-uu(i,j,k)*y_ksi(i,j)
           rovc=rho_n(i,j,k)*vc

           Fnutil(i,j,k) = rovc*nutil(i,j,k)
        enddo
     enddo
  enddo

  ! derivatives along eta
  ! ---------------------
  j=1
  do i=1,2
     do k=ndz_e,nfz_e  
        Knutil(i,j,k) = (Fnutil(i,j+1,k) - Fnutil(i,j,k)) + Knutil(i,j,k)
     enddo
  enddo

  j=2
  do i=1,2
     do k=ndz_e,nfz_e  
        Knutil(i,j,k) = 0.5_wp*(Fnutil(i,j+1,k) - Fnutil(i,j-1,k)) + Knutil(i,j,k)
     enddo
  enddo

  ! ~> Second pass: order reduction for metrics & 5-pt derivatives
  ! ==============================================================

  ! modified fluxes along ksi
  ! -------------------------
  do k=ndz_e,nfz_e
     do j=1,2
        do i=-1,4
           ! contravariant velocity
           vc=uu(i,j,k)*y_eta_imin_jmin(i,j)-vv(i,j,k)*x_eta_imin_jmin(i,j)
           rovc=rho_n(i,j,k)*vc

           Fnutil(i,j,k) = rovc*nutil(i,j,k)
        enddo
     enddo
  enddo
  
  ! derivatives along ksi
  ! ---------------------
  do k=ndz_e,nfz_e
     do j=1,2
        do i=1,2
           Knutil(i,j,k) = a5(1) * ( Fnutil(i+1,j,k)-Fnutil(i-1,j,k) ) &
                            + a5(2) * ( Fnutil(i+2,j,k)-Fnutil(i-2,j,k) ) + Knutil(i,j,k)
        enddo
     enddo
  enddo

  ! modified fluxes along eta
  ! -------------------------
  do k=ndz_e,nfz_e
     do i=1,2
        do j=-1,4
           ! contravariant velocity
           vc=vv(i,j,k)*x_ksi_imin_jmin(i,j)-uu(i,j,k)*y_ksi_imin_jmin(i,j)
           rovc=rho_n(i,j,k)*vc

           Fnutil(i,j,k) = rovc*nutil(i,j,k)
        enddo
     enddo
  enddo

  ! derivatives along eta
  ! ---------------------
  do k=ndz_e,nfz_e
     do j=1,2
        do i=1,2
           Knutil(i,j,k) = a5(1) * ( Fnutil(i,j+1,k)-Fnutil(i,j-1,k) ) &
                            + a5(2) * ( Fnutil(i,j+2,k)-Fnutil(i,j-2,k) ) + Knutil(i,j,k)
        enddo
     enddo
  enddo

  ! Average of the two derivative evaluations
  ! =========================================
  do k=ndz_e,nfz_e
     do j=1,2
        do i=1,2
           Knutil(i,j,k) = Knutil(i,j,k)*0.5
        enddo
     enddo
  enddo

end subroutine flux_euler_w_imin_jmin_5pts_SA_c

!===============================================================================
subroutine flux_euler_w_imin_jmax_5pts_SA_c
!===============================================================================
  !> Derivatives of Euler fluxes for edge I_MIN/J_MAX close to two walls
  !> - curvilinear coordinate - 5-point stencil -
!===============================================================================
  use mod_flow
  use mod_coeff_deriv
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j,k
  real(wp) :: vc,rovc
  ! ---------------------------------------------------------------------------

  ! ~> First pass: 5-pt metrics & order reduction for derivatives
  ! =============================================================

  ! fluxes along ksi
  ! ----------------
  do k=ndz_e,nfz_e
     do j=ny-1,ny
        do i=-1,4
           ! contravariant velocity
           vc=uu(i,j,k)*y_eta(i,j)-vv(i,j,k)*x_eta(i,j)
           rovc=rho_n(i,j,k)*vc

           Fnutil(i,j,k) = rovc*nutil(i,j,k)
        enddo
     enddo
  enddo
  
  ! derivatives along ksi
  ! ---------------------
  i=1
  do j=ny-2,ny
     do k=ndz_e,nfz_e
        Knutil(i,j,k) = ( Fnutil(i+1,j,k)-Fnutil(i,j,k))
     enddo
  enddo

  i=2
  do j=ny-2,ny
     do k=ndz_e,nfz_e
        Knutil(i,j,k) = 0.5_wp*( Fnutil(i+1,j,k)-Fnutil(i-1,j,k))
     enddo
  enddo

  ! modified fluxes along eta
  ! -------------------------
  do k=ndz_e,nfz_e
     do i=1,2
        do j=ny-3,ny+2
           ! contravariant velocity
           vc=vv(i,j,k)*x_ksi(i,j)-uu(i,j,k)*y_ksi(i,j)
           rovc=rho_n(i,j,k)*vc

           Fnutil(i,j,k) = rovc*nutil(i,j,k)
        enddo
     enddo
  enddo

  ! derivatives along eta
  ! ---------------------
  j=ny
  do i=1,2
     do k=ndz_e,nfz_e  
        Knutil(i,j,k) = (Fnutil(i,j,k)-Fnutil(i,j-1,k)) + Knutil(i,j,k)
     enddo
  enddo

  j=ny-1
  do i=1,2
     do k=ndz_e,nfz_e  
        Knutil(i,j,k) = 0.5_wp*(Fnutil(i,j+1,k)-Fnutil(i,j-1,k)) + Knutil(i,j,k)
     enddo
  enddo

  ! ~> Second pass: order reduction for metrics & 9-pt derivatives
  ! ==============================================================

  ! modified fluxes along ksi
  ! -------------------------
  do k=ndz_e,nfz_e
     do j=ny-1,ny
        do i=-1,2
           ! contravariant velocity
           vc=uu(i,j,k)*y_eta_imin_jmax(i,j)-vv(i,j,k)*x_eta_imin_jmax(i,j)
           rovc=rho_n(i,j,k)*vc

           Fnutil(i,j,k) = rovc*nutil(i,j,k)
        enddo
     enddo
  enddo
  
  ! derivatives along ksi
  ! ---------------------
  do k=ndz_e,nfz_e
     do j=ny-1,ny
        do i=1,2
           Knutil(i,j,k) = a5(1) * ( Fnutil(i+1,j,k)-Fnutil(i-1,j,k) ) &
                            + a5(2) * ( Fnutil(i+2,j,k)-Fnutil(i-2,j,k) ) + Knutil(i,j,k)
        enddo
     enddo
  enddo

  ! modified fluxes along eta
  ! -------------------------
  do k=ndz_e,nfz_e
     do i=1,2
        do j=ny-3,ny+2
           ! contravariant velocity
           vc=vv(i,j,k)*x_ksi_imin_jmax(i,j)-uu(i,j,k)*y_ksi_imin_jmax(i,j)
           rovc=rho_n(i,j,k)*vc

           Fnutil(i,j,k) = rovc*nutil(i,j,k)
        enddo
     enddo
  enddo

  ! derivatives along eta
  ! ---------------------
  do k=ndz_e,nfz_e
     do j=ny-1,ny
        do i=1,2
           Knutil(i,j,k) = a5(1) * ( Fnutil(i,j+1,k)-Fnutil(i,j-1,k) ) &
                            + a5(2) * ( Fnutil(i,j+2,k)-Fnutil(i,j-2,k) ) + Knutil(i,j,k)
        enddo
     enddo
  enddo

  ! Average of the two derivative evaluations
  ! =========================================
  do k=ndz_e,nfz_e
     do j=ny-1,ny
        do i=1,2
           Knutil(i,j,k) = Knutil(i,j,k)*0.5
        enddo
     enddo
  enddo

end subroutine flux_euler_w_imin_jmax_5pts_SA_c

!===============================================================================
subroutine flux_euler_w_imax_jmin_5pts_SA_c
!===============================================================================
  !> Derivatives of Euler fluxes for edge I_MAX/J_MIN close to two walls
  !> - curvilinear coordinate - 5-point stencil -
!===============================================================================
  use mod_flow
  use mod_coeff_deriv
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j,k
  real(wp) :: vc,rovc
  ! ---------------------------------------------------------------------------

  ! ~> First pass: 5-pt metrics & order reduction for derivatives
  ! =============================================================

  ! fluxes along ksi
  ! ----------------
  do k=ndz_e,nfz_e
     do j=1,2
        do i=ny-3,nx+2
           ! contravariant velocity
           vc=uu(i,j,k)*y_eta(i,j)-vv(i,j,k)*x_eta(i,j)
           rovc=rho_n(i,j,k)*vc

           Fnutil(i,j,k) = rovc*nutil(i,j,k)
        enddo
     enddo
  enddo
  
  ! derivatives along ksi
  ! ---------------------
  i=nx-1
  do j=1,2
     do k=ndz_e,nfz_e
        Knutil(i,j,k) = ( Fnutil(i+1,j,k)-Fnutil(i,j,k))
     enddo
  enddo

  i=nx
  do j=1,2
     do k=ndz_e,nfz_e
        Knutil(i,j,k) = 0.5_wp*( Fnutil(i+1,j,k)-Fnutil(i-1,j,k))
     enddo
  enddo

  ! modified fluxes along eta
  ! -------------------------
  do k=ndz_e,nfz_e
     do i=nx-1,nx
        do j=-1,4
           ! contravariant velocity
           vc=vv(i,j,k)*x_ksi(i,j)-uu(i,j,k)*y_ksi(i,j)
           rovc=rho_n(i,j,k)*vc

           Fnutil(i,j,k) = rovc*nutil(i,j,k)
        enddo
     enddo
  enddo

  ! derivatives along eta
  ! ---------------------
  j=1
  do i=nx-1,nx
     do k=ndz_e,nfz_e  
        Knutil(i,j,k) = (Fnutil(i,j,k)-Fnutil(i,j-1,k)) + Knutil(i,j,k)
     enddo
  enddo

  j=2
  do i=1,2
     do k=ndz_e,nfz_e  
        Knutil(i,j,k) = 0.5_wp*(Fnutil(i,j+1,k)-Fnutil(i,j-1,k)) + Knutil(i,j,k)
     enddo
  enddo

  ! ~> Second pass: order reduction for metrics & 9-pt derivatives
  ! ==============================================================

  ! modified fluxes along ksi
  ! -------------------------
  do k=ndz_e,nfz_e
     do j=1,2
        do i=nx-3,nx+2
           ! contravariant velocity
           vc=uu(i,j,k)*y_eta_imin_jmax(i,j)-vv(i,j,k)*x_eta_imin_jmax(i,j)
           rovc=rho_n(i,j,k)*vc

           Fnutil(i,j,k) = rovc*nutil(i,j,k)
        enddo
     enddo
  enddo
  
  ! derivatives along ksi
  ! ---------------------
  do k=ndz_e,nfz_e
     do j=1,2
        do i=nx-1,nx
           Knutil(i,j,k) = a5(1) * ( Fnutil(i+1,j,k)-Fnutil(i-1,j,k) ) &
                            + a5(2) * ( Fnutil(i+2,j,k)-Fnutil(i-2,j,k) ) + Knutil(i,j,k)
        enddo
     enddo
  enddo

  ! modified fluxes along eta
  ! -------------------------
  do k=ndz_e,nfz_e
     do i=nx-1,nx
        do j=-1,4
           ! contravariant velocity
           vc=vv(i,j,k)*x_ksi_imin_jmax(i,j)-uu(i,j,k)*y_ksi_imin_jmax(i,j)
           rovc=rho_n(i,j,k)*vc

           Fnutil(i,j,k) = rovc*nutil(i,j,k)
        enddo
     enddo
  enddo

  ! derivatives along eta
  ! ---------------------
  do k=ndz_e,nfz_e
     do j=1,2
        do i=nx-1,nx
           Knutil(i,j,k) = a5(1) * ( Fnutil(i,j+1,k)-Fnutil(i,j-1,k) ) &
                            + a5(2) * ( Fnutil(i,j+2,k)-Fnutil(i,j-2,k) ) + Knutil(i,j,k)
        enddo
     enddo
  enddo

  ! Average of the two derivative evaluations
  ! =========================================
  do k=ndz_e,nfz_e
     do j=1,2
        do i=nx-1,nx
           Knutil(i,j,k) = Knutil(i,j,k)*0.5
        enddo
     enddo
  enddo

end subroutine flux_euler_w_imax_jmin_5pts_SA_c

!===============================================================================
subroutine flux_euler_w_imax_jmax_5pts_SA_c
!===============================================================================
  !> Derivatives of Euler fluxes for edge I_MAX/J_MAX close to two walls
  !> - curvilinear coordinate - 5-point stencil -
!===============================================================================
  use mod_flow
  use mod_coeff_deriv
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j,k
  real(wp) :: vc,rovc
  ! ---------------------------------------------------------------------------

  ! ~> First pass: 5-pt metrics & order reduction for derivatives
  ! =============================================================

  ! modified fluxes along ksi
  ! -------------------------
  do k=ndz_e,nfz_e
     do j=ny-1,ny
        do i=nx-3,nx+2
           ! contravariant velocity
           vc=uu(i,j,k)*y_eta(i,j)-vv(i,j,k)*x_eta(i,j)
           rovc=rho_n(i,j,k)*vc

           Fnutil(i,j,k) = rovc*nutil(i,j,k)
        enddo
     enddo
  enddo
  
  ! derivatives along ksi
  ! ---------------------
  i=nx-1
  do j=ny-1,ny
     do k=ndz_e,nfz_e
        Knutil(i,j,k) = 0.5_wp*( Fnutil(i+1,j,k)-Fnutil(i-1,j,k))
     enddo
  enddo

  i=nx
  do j=ny-1,ny
     do k=ndz_e,nfz_e
        Knutil(i,j,k) = ( Fnutil(i,j,k)-Fnutil(i-1,j,k))
     enddo
  enddo

  ! modified fluxes along eta
  ! -------------------------
  do k=ndz_e,nfz_e
     do i=nx-1,nx
        do j=ny-3,ny+2
           ! contravariant velocity
           vc=vv(i,j,k)*x_ksi(i,j)-uu(i,j,k)*y_ksi(i,j)
           rovc=rho_n(i,j,k)*vc

           Fnutil(i,j,k)  =  rovc*nutil(i,j,k)
        enddo
     enddo
  enddo

  ! derivatives along eta
  ! ---------------------
  j=ny
  do i=nx-1,nx
     do k=ndz_e,nfz_e  
        Knutil(i,j,k)  = (Fnutil(i,j,k)-Fnutil(i,j-1,k)) + Knutil(i,j,k)
     enddo
  enddo

  j=ny-1
  do i=nx-1,nx
     do k=ndz_e,nfz_e  
        Knutil(i,j,k)  = 0.5_wp*(Fnutil(i,j+1,k)-Fnutil(i,j-1,k)) + Knutil(i,j,k)
     enddo
  enddo

  ! ~> Second pass: order reduction for metrics & 9-pt derivatives
  ! ==============================================================

  ! modified fluxes along ksi
  ! -------------------------
  do k=ndz_e,nfz_e
     do j=ny-1,ny
        do i=nx-3,nx+2
           ! contravariant velocity
           vc=uu(i,j,k)*y_eta_imax_jmax(i,j)-vv(i,j,k)*x_eta_imax_jmax(i,j)
           rovc=rho_n(i,j,k)*vc

           Fnutil(i,j,k) = rovc*nutil(i,j,k)
        enddo
     enddo
  enddo
  
  ! derivatives along ksi
  ! ---------------------
  do k=ndz_e,nfz_e
     do j=ny-1,ny
        do i=nx-1,nx
           Knutil(i,j,k) = a5(1) * ( Fnutil(i+1,j,k)-Fnutil(i-1,j,k) ) &
                            + a5(2) * ( Fnutil(i+2,j,k)-Fnutil(i-2,j,k) ) + Knutil(i,j,k)
        enddo
     enddo
  enddo

  ! modified fluxes along eta
  ! -------------------------
  do k=ndz_e,nfz_e
     do i=nx-1,nx
        do j=ny-3,ny+2
           ! contravariant velocity
           vc=vv(i,j,k)*x_ksi_imax_jmax(i,j)-uu(i,j,k)*y_ksi_imax_jmax(i,j)
           rovc=rho_n(i,j,k)*vc

           Fnutil(i,j,k) = rovc*nutil(i,j,k)
        enddo
     enddo
  enddo

  ! derivatives along eta
  ! ---------------------
  do k=ndz_e,nfz_e
     do j=ny-1,ny
        do i=nx-1,nx
           Knutil(i,j,k) = a5(1) * ( Fnutil(i,j+1,k)-Fnutil(i,j-1,k) ) &
                            + a5(2) * ( Fnutil(i,j+2,k)-Fnutil(i,j-2,k) ) + Knutil(i,j,k)
        enddo
     enddo
  enddo

  ! Average of the two derivative evaluations
  ! =========================================
  do k=ndz_e,nfz_e
     do j=ny-1,ny
        do i=nx-1,nx
           Knutil(i,j,k) = Knutil(i,j,k)*0.5
        enddo
     enddo
  enddo

end subroutine flux_euler_w_imax_jmax_5pts_SA_c
