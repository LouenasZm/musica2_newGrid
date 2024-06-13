!==============================================================================
subroutine flux_visc_3pts_SA_transition_c
!==============================================================================
  !> Compute derivatives of viscous fluxes (3-point stencil - order 4)
  !> - curvilinear version -
!==============================================================================
  use mod_coeff_deriv
  use mod_flow
  use mod_rans
  implicit none
  ! ---------------------------------------------------------------------------
  integer  :: i,j,k
  real(wp) :: nu,dnutilc_ksi,dnutilc_eta, mu
  real(wp) :: dgammac_eta, dgammac_ksi, dre_thetac_eta, dre_thetac_ksi 
  real(wp), parameter :: TWO_THIRDS=2.0_wp/3.0_wp
  ! ---------------------------------------------------------------------------

  ! Model constants
  ! ===============
  sig = TWO_THIRDS

  ! Compute viscous fluxes
  ! ======================
  do k=ndzt_v_r,nfzt_v_r
     do j=ndyt_v_r,nfyt_v_r
        do i=ndxt_v_r,nfxt_v_r          
           ! kinematic viscosity
           nu   = visc(i,j,k)/rho_n(i,j,k)
           mu   = visc(i,j,k)
           ! contravariant grad(Nutilde)
           dnutilc_ksi =  dnutilx(i,j,k)*y_eta_v(i,j)-dnutily(i,j,k)*x_eta_v(i,j)
           dnutilc_eta = -dnutilx(i,j,k)*y_ksi_v(i,j)+dnutily(i,j,k)*x_ksi_v(i,j)
           
           ! contravariant grad(gamma) 
           dgammac_ksi =  dgammax(i,j,k)*y_eta_v(i,j)-dgammay(i,j,k)*x_eta_v(i,j)
           dgammac_eta = -dgammax(i,j,k)*y_ksi_v(i,j)+dgammay(i,j,k)*x_ksi_v(i,j)
           
           ! contravariant grad(re_theta)
           dre_thetac_ksi =  dre_thetax(i,j,k)*y_eta_v(i,j)-dre_thetay(i,j,k)*x_eta_v(i,j)
           dre_thetac_eta = -dre_thetax(i,j,k)*y_ksi_v(i,j)+dre_thetay(i,j,k)*x_ksi_v(i,j)
           
           ! viscous fluxes along ksi
            Fnutil(i,j,k)       = -(nu+nutil_n(i,j,k))/sig*dnutilc_ksi
            Fgamma(i,j,k)       = -(mu+mut(i,j,k))*dgammac_ksi
            Fre_theta(i,j,k)    = -sigma_theta_t*(mu+mut(i,j,k))*dre_thetac_ksi
           ! viscous fluxes along eta
           Gnutil(i,j,k)    = -(nu+nutil_n(i,j,k))/sig*dnutilc_eta
           Ggamma(i,j,k)    = -(mu+mut(i,j,k))*dgammac_eta
           Gre_theta(i,j,k) = -sigma_theta_t*(mu+mut(i,j,k))*dre_thetac_eta
           ! viscous fluxes along z
           Hnutil(i,j,k)    = -(nu+nutil_n(i,j,k))/sig*dnutilz(i,j,k)
           Hgamma(i,j,k)    = -(mu+mut(i,j,k))*dgammaz(i,j,k)
           Hre_theta(i,j,k) = -sigma_theta_t*(mu+mut(i,j,k))*dre_thetaz(i,j,k)
        enddo
     enddo
  enddo
  
  ! Compute derivatives of viscous fluxes along ksi
  ! ===============================================
  if (is_bc_wall(1,1)) then
     i=1
     do k=ndz_v_r,nfz_v_r
        do j=ndy_v_r,nfy_v_r
           Knutil(i,j,k) = (Fnutil(i+1,j,k)-Fnutil(i,j,k))*ijacob_v(i,j) + Knutil(i,j,k)
           Kgamma(i,j,k) = (Fgamma(i+1,j,k)-Fgamma(i,j,k))*ijacob_v(i,j) + Kgamma(i,j,k)
           Kre_theta(i,j,k) = (Fre_theta(i+1,j,k)-Fre_theta(i,j,k))*ijacob_v(i,j) + Kre_theta(i,j,k)
        enddo
     enddo
  endif

  do k=ndz_v_r,nfz_v_r
     do j=ndy_v_r,nfy_v_r
        do i=ndx_vi_r,nfx_vi_r
!           Knutil(i,j,k) = (Fnutil(i+1,j,k)-2.0_wp*Fnutil(i,j,k)+Fnutil(i-1,j,k))*ijacob_v(i,j) &
            Knutil(i,j,k) = (Fnutil(i+1,j,k)-Fnutil(i-1,j,k))*ijacob_v(i,j) &
            + Knutil(i,j,k)
            Kgamma(i,j,k) = (Fgamma(i+1,j,k)-Fgamma(i-1,j,k))*ijacob_v(i,j) &
                             + Kgamma(i,j,k)
            Kre_theta(i,j,k) = (Fre_theta(i+1,j,k)-Fre_theta(i-1,j,k))*ijacob_v(i,j) &
                             + Kre_theta(i,j,k)
        enddo
     enddo
  enddo

  if (is_bc_wall(1,2)) then
     i=nx
     do k=ndz_v_r,nfz_v_r
        do j=ndy_v_r,nfy_v_r
           Knutil(i,j,k) = (Fnutil(i,j,k)-Fnutil(i-1,j,k))*ijacob_v(i,j) + Knutil(i,j,k)
           Kgamma(i,j,k) = (Fgamma(i,j,k)-Fgamma(i-1,j,k))*ijacob_v(i,j) + Kgamma(i,j,k)
           Kre_theta(i,j,k) = (Fre_theta(i,j,k)-Fre_theta(i-1,j,k))*ijacob_v(i,j) + Kre_theta(i,j,k)
        enddo
     enddo
  endif

  ! Compute derivatives of viscous fluxes along eta
  ! ===============================================
  if (is_bc_wall(2,1)) then
     j=1
     do k=ndz_v_r,nfz_v_r
        do i=ndx_v_r,nfx_v_r
           Knutil(i,j,k) = (Gnutil(i,j+1,k)-Gnutil(i,j,k))*ijacob_v(i,j) + Knutil(i,j,k)
           Kgamma(i,j,k) = (Ggamma(i,j+1,k)-Ggamma(i,j,k))*ijacob_v(i,j) + Kgamma(i,j,k)
           Kre_theta(i,j,k) = (Gre_theta(i,j+1,k)-Gre_theta(i,j,k))*ijacob_v(i,j) + Kre_theta(i,j,k)
        enddo
     enddo
  endif

  do k=ndz_v_r,nfz_v_r
     do j=ndy_vi_r,nfy_vi_r
        do i=ndx_v_r,nfx_v_r
!           Knutil(i,j,k) = (Gnutil(i,j+1,k)-2.0_wp*Gnutil(i,j,k)+Gnutil(i,j-1,k))*ijacob_v(i,j) & 
            Knutil(i,j,k) = (Gnutil(i,j+1,k)-Gnutil(i,j-1,k))*ijacob_v(i,j) & 
                             + Knutil(i,j,k)
            Kgamma(i,j,k) = (Ggamma(i,j+1,k)-Ggamma(i,j-1,k))*ijacob_v(i,j) & 
                             + Kgamma(i,j,k)
            Kre_theta(i,j,k) = (Gre_theta(i,j+1,k)-Gre_theta(i,j-1,k))*ijacob_v(i,j) & 
                             + Kre_theta(i,j,k)
        enddo
     enddo
  enddo

  if (is_bc_wall(2,2)) then
     j=ny
     do k=ndz_v_r,nfz_v_r
        do i=ndx_v_r,nfx_v_r
           Knutil(i,j,k) = (Gnutil(i,j,k)-Gnutil(i,j-1,k))*ijacob_v(i,j) + Knutil(i,j,k)
           Kgamma(i,j,k) = (Ggamma(i,j,k)-Ggamma(i,j-1,k))*ijacob_v(i,j) + Kgamma(i,j,k)
           Kre_theta(i,j,k) = (Gre_theta(i,j,k)-Gre_theta(i,j-1,k))*ijacob_v(i,j) + Kre_theta(i,j,k)
        enddo
     enddo
  endif

  if (is_2D) return

  ! Compute derivatives of viscous fluxes along z
  ! =============================================
  if (is_bc_wall(3,1)) then
     k=1
     do j=ndy_v_r,nfy_v_r
        do i=ndx_v_r,nfx_v_r
           Knutil(i,j,k) = (Hnutil(i,j,k+1)-Hnutil(i,j,k))*ijacob_v(i,j) + Knutil(i,j,k)
           Kgamma(i,j,k) = (Hgamma(i,j,k+1)-Hgamma(i,j,k))*ijacob_v(i,j) + Kgamma(i,j,k)
           Kre_theta(i,j,k) = (Hre_theta(i,j,k+1)-Hre_theta(i,j,k))*ijacob_v(i,j) + Kre_theta(i,j,k)
        enddo
     enddo
  endif

  do k=ndz_vi_r,nfz_vi_r
     do j=ndy_v_r,nfy_v_r
        do i=ndx_v_r,nfx_v_r
!           Knutil(i,j,k) = (Hnutil(i,j,k+1)-2.0_wp*Hnutil(i,j,k)+Hnutil(i,j,k-1))*ijacob_v(i,j) &
            Knutil(i,j,k) = (Hnutil(i,j,k+1)-Hnutil(i,j,k-1))*ijacob_v(i,j) &
                             + Knutil(i,j,k)
            Kgamma(i,j,k) = (Hgamma(i,j,k+1)-Hgamma(i,j,k-1))*ijacob_v(i,j) &
                             + Kgamma(i,j,k)
            Kre_theta(i,j,k) = (Hre_theta(i,j,k+1)-Hre_theta(i,j,k-1))*ijacob_v(i,j) &
                             + Kre_theta(i,j,k)                                                          
        enddo
     enddo
  enddo

  if (is_bc_wall(3,2)) then
     k=nz
     do j=ndy_v_r,nfy_v_r
        do i=ndx_v_r,nfx_v_r
           Knutil(i,j,k) = (Hnutil(i,j,k)-Hnutil(i,j,k-1))*ijacob_v(i,j) + Knutil(i,j,k)
           Kgamma(i,j,k) = (Hgamma(i,j,k)-Hgamma(i,j,k-1))*ijacob_v(i,j) + Kgamma(i,j,k)
           Kre_theta(i,j,k) = (Hre_theta(i,j,k)-Hre_theta(i,j,k-1))*ijacob_v(i,j) + Kre_theta(i,j,k)
        enddo
     enddo
  endif

end subroutine flux_visc_3pts_SA_transition_c

