!==============================================================================
subroutine flux_visc_3pts_SA_transition_c3
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
  real(wp) :: nu,dnutilc_ksi,dnutilc_eta,dnutilc_phi, mu
  real(wp) :: dgammac_eta, dgammac_ksi, dgammac_phi 
  real(wp) :: dre_thetac_eta, dre_thetac_ksi, dre_thetac_phi
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
           nu    = visc(i,j,k)/rho_n(i,j,k)
           mu    = visc(i,j,k)
           ! compute derivatives 
           dnutilc_ksi = dnutilx(i,j,k)*ksi_x_v(i,j,k)+dnutily(i,j,k)*ksi_y_v(i,j,k)+dnutilz(i,j,k)*ksi_z_v(i,j,k)
           dnutilc_eta = dnutilx(i,j,k)*eta_x_v(i,j,k)+dnutily(i,j,k)*eta_y_v(i,j,k)+dnutilz(i,j,k)*eta_z_v(i,j,k)
           dnutilc_phi = dnutilx(i,j,k)*phi_x_v(i,j,k)+dnutily(i,j,k)*phi_y_v(i,j,k)+dnutilz(i,j,k)*phi_z_v(i,j,k)

           dgammac_ksi = dgammax(i,j,k)*ksi_x_v(i,j,k)+dgammay(i,j,k)*ksi_y_v(i,j,k)+dgammaz(i,j,k)*ksi_z_v(i,j,k)
           dgammac_eta = dgammax(i,j,k)*eta_x_v(i,j,k)+dgammay(i,j,k)*eta_y_v(i,j,k)+dgammaz(i,j,k)*eta_z_v(i,j,k)
           dgammac_phi = dgammax(i,j,k)*phi_x_v(i,j,k)+dgammay(i,j,k)*phi_y_v(i,j,k)+dgammaz(i,j,k)*phi_z_v(i,j,k)

           dre_thetac_ksi = dre_thetax(i,j,k)*ksi_x_v(i,j,k)+dre_thetay(i,j,k)*ksi_y_v(i,j,k)+dre_thetaz(i,j,k)*ksi_z_v(i,j,k)
           dre_thetac_eta = dre_thetax(i,j,k)*eta_x_v(i,j,k)+dre_thetay(i,j,k)*eta_y_v(i,j,k)+dre_thetaz(i,j,k)*eta_z_v(i,j,k)
           dre_thetac_phi = dre_thetax(i,j,k)*phi_x_v(i,j,k)+dre_thetay(i,j,k)*phi_y_v(i,j,k)+dre_thetaz(i,j,k)*phi_z_v(i,j,k)
           ! viscous fluxes along ksi
           Fnutil(i,j,k)    = -(nu+nutil_n(i,j,k))/sig*dnutilc_ksi
           Fgamma(i,j,k)    = -(mu+mut(i,j,k))*dgammac_ksi
           Fre_theta(i,j,k) = -sigma_theta_t*(mu+mut(i,j,k))*dre_thetac_ksi

           ! viscous fluxes along eta
           Gnutil(i,j,k)    = -(nu+nutil_n(i,j,k))/sig*dnutilc_eta
           Ggamma(i,j,k)    = -(mu+mut(i,j,k))*dgammac_eta
           Gre_theta(i,j,k) = -sigma_theta_t*(mu+mut(i,j,k))*dre_thetac_eta
           
           ! viscous fluxes along phi
           Hnutil(i,j,k)    = -(nu+nutil_n(i,j,k))/sig*dnutilc_phi
           Hgamma(i,j,k)    = -(mu+mut(i,j,k))*dgammac_phi
           Hre_theta(i,j,k) = -sigma_theta_t*(mu+mut(i,j,k))*dre_thetac_phi

        enddo
     enddo
  enddo
  
  ! Compute derivatives of viscous fluxes along ksi
  ! ===============================================
  if (is_bc_wall(1,1)) then
     i=1
     do k=ndz_v_r,nfz_v_r
        do j=ndy_v_r,nfy_v_r
            Knutil(i,j,k) = (a02(1)*Fnutil(i,j,k) + a02(2)*Fnutil(i+1,j,k) + &
                            a02(3)*Fnutil(i+2,j,k))*ijacob3_v(i,j,k) + Knutil(i,j,k)
            Kgamma(i,j,k) = (a02(1)*Fgamma(i,j,k) + a02(2)*Fgamma(i+1,j,k) + &
                            a02(3)*Fgamma(i+2,j,k))*ijacob3_v(i,j,k) + Kgamma(i,j,k)
            Kre_theta(i,j,k) = (a02(1)*Fre_theta(i,j,k) + a02(2)*Fre_theta(i+1,j,k) + &
                            a02(3)*Fre_theta(i+2,j,k))*ijacob3_v(i,j,k) + Kre_theta(i,j,k)
        enddo
     enddo
  endif

  do k=ndz_v_r,nfz_v_r
     do j=ndy_v_r,nfy_v_r
        do i=ndx_vi_r,nfx_vi_r
            Knutil(i,j,k) = a3(1)*(Fnutil(i+1,j,k)-Fnutil(i-1,j,k))*ijacob3_v(i,j,k) &
                             + Knutil(i,j,k)
            Kgamma(i,j,k) = a3(1)*(Fgamma(i+1,j,k)-Fgamma(i-1,j,k))*ijacob3_v(i,j,k) &
                             + Kgamma(i,j,k)
            Kre_theta(i,j,k) = a3(1)*(Fre_theta(i+1,j,k)-Fre_theta(i-1,j,k))*ijacob3_v(i,j,k) &
                             + Kre_theta(i,j,k)
        enddo
     enddo
  enddo

  if (is_bc_wall(1,2)) then
     i=nx
     do k=ndz_v_r,nfz_v_r
        do j=ndy_v_r,nfy_v_r
            Knutil(i,j,k) = (a20(1)*Fnutil(i,j,k) + a20(2)*Fnutil(i-1,j,k) + &
                            a20(3)*Fnutil(i-2,j,k))*ijacob3_v(i,j,k) + Knutil(i,j,k)
            Kgamma(i,j,k) = (a20(1)*Fgamma(i,j,k) + a20(2)*Fgamma(i-1,j,k) + &
                            a20(3)*Fgamma(i-2,j,k))*ijacob3_v(i,j,k) + Kgamma(i,j,k)
            Kre_theta(i,j,k) = (a20(1)*Fre_theta(i,j,k) + a20(2)*Fre_theta(i-1,j,k) + &
                            a20(3)*Fre_theta(i-2,j,k))*ijacob3_v(i,j,k) + Kre_theta(i,j,k)
        enddo
     enddo
  endif

  ! Compute derivatives of viscous fluxes along eta
  ! ===============================================
  if (is_bc_wall(2,1)) then
     j=1
     do k=ndz_v_r,nfz_v_r
        do i=ndx_v_r,nfx_v_r
            Knutil(i,j,k) = (a02(1)*Gnutil(i,j,k) + a02(2)*Gnutil(i,j+1,k) + &
                            a02(3)*Gnutil(i,j+2,k))*ijacob3_v(i,j,k) + Knutil(i,j,k)
            Kgamma(i,j,k) = (a02(1)*Ggamma(i,j,k) + a02(2)*Ggamma(i,j+1,k) + &
                            a02(3)*Ggamma(i,j+2,k))*ijacob3_v(i,j,k) + Kgamma(i,j,k)
            Kre_theta(i,j,k) = (a02(1)*Gre_theta(i,j,k) + a02(2)*Gre_theta(i,j+1,k) + &
                            a02(3)*Gre_theta(i,j+2,k))*ijacob3_v(i,j,k) + Kre_theta(i,j,k)
        enddo
     enddo
  endif

  do k=ndz_v_r,nfz_v_r
     do j=ndy_vi_r,nfy_vi_r
        do i=ndx_v_r,nfx_v_r
            Knutil(i,j,k) = a3(1)*(Gnutil(i,j+1,k)-Gnutil(i,j-1,k))*ijacob3_v(i,j,k) &
                             + Knutil(i,j,k)
            Kgamma(i,j,k) = a3(1)*(Ggamma(i,j+1,k)-Ggamma(i,j-1,k))*ijacob3_v(i,j,k) &
                             + Kgamma(i,j,k)
            Kre_theta(i,j,k) = a3(1)*(Gre_theta(i,j+1,k)-Gre_theta(i,j-1,k))*ijacob3_v(i,j,k) &
                             + Kre_theta(i,j,k)
        enddo
     enddo
  enddo

  if (is_bc_wall(2,2)) then
     j=ny
     do k=ndz_v_r,nfz_v_r
        do i=ndx_v_r,nfx_v_r
            Knutil(i,j,k) = (a20(1)*Gnutil(i,j,k) + a20(2)*Gnutil(i,j-1,k) + &
                            a20(3)*Gnutil(i,j-2,k))*ijacob3_v(i,j,k) + Knutil(i,j,k)
            Kgamma(i,j,k) = (a20(1)*Ggamma(i,j,k) + a20(2)*Ggamma(i,j-1,k) + &
                            a20(3)*Ggamma(i,j-2,k))*ijacob3_v(i,j,k) + Kgamma(i,j,k)
            Kre_theta(i,j,k) = (a20(1)*Gre_theta(i,j,k) + a20(2)*Gre_theta(i,j-1,k) + &
                            a20(3)*Gre_theta(i,j-2,k))*ijacob3_v(i,j,k) + Kre_theta(i,j,k)
        enddo
     enddo
  endif

  ! Compute derivatives of viscous fluxes along phi
  ! ===============================================
  if (is_bc_wall(3,1)) then
     k=1
     do j=ndy_v_r,nfy_v_r
        do i=ndx_v_r,nfx_v_r
            Knutil(i,j,k) = (a02(1)*Hnutil(i,j,k) + a02(2)*Hnutil(i,j,k+1) + &
                            a02(3)*Hnutil(i,j,k+2))*ijacob3_v(i,j,k) + Knutil(i,j,k)
            Kgamma(i,j,k) = (a02(1)*Hgamma(i,j,k) + a02(2)*Hgamma(i,j,k+1) + &
                            a02(3)*Hgamma(i,j,k+2))*ijacob3_v(i,j,k) + Kgamma(i,j,k)
            Kre_theta(i,j,k) = (a02(1)*Hre_theta(i,j,k) + a02(2)*Hre_theta(i,j,k+1) + &
                            a02(3)*Hre_theta(i,j,k+2))*ijacob3_v(i,j,k) + Kre_theta(i,j,k)
        enddo
     enddo
  endif

  do k=ndz_vi_r,nfz_vi_r
     do j=ndy_v_r,nfy_v_r
        do i=ndx_v_r,nfx_v_r
            Knutil(i,j,k) = a3(1)*(Hnutil(i,j,k+1)-Hnutil(i,j,k-1))*ijacob3_v(i,j,k) &
                             + Knutil(i,j,k)
            Kgamma(i,j,k) = a3(1)*(Hgamma(i,j,k+1)-Hgamma(i,j,k-1))*ijacob3_v(i,j,k) &
                             + Kgamma(i,j,k)
            Kre_theta(i,j,k) = a3(1)*(Hre_theta(i,j,k+1)-Hre_theta(i,j,k-1))*ijacob3_v(i,j,k) &
                             + Kre_theta(i,j,k)
        enddo
     enddo
  enddo

  if (is_bc_wall(3,2)) then
     k=nz
     do j=ndy_v_r,nfy_v_r
        do i=ndx_v_r,nfx_v_r
            Knutil(i,j,k) = (a20(1)*Hnutil(i,j,k) + a20(2)*Hnutil(i,j,k-1) + &
                            a20(3)*Hnutil(i,j,k-2))*ijacob3_v(i,j,k) + Knutil(i,j,k)
            Kgamma(i,j,k) = (a20(1)*Hgamma(i,j,k) + a20(2)*Hgamma(i,j,k-1) + &
                            a20(3)*Hgamma(i,j,k-2))*ijacob3_v(i,j,k) + Kgamma(i,j,k)
            Kre_theta(i,j,k) = (a20(1)*Hre_theta(i,j,k) + a20(2)*Hre_theta(i,j,k-1) + &
                            a20(3)*Hre_theta(i,j,k-2))*ijacob3_v(i,j,k) + Kre_theta(i,j,k)
        enddo
     enddo
  endif

end subroutine flux_visc_3pts_SA_transition_c3