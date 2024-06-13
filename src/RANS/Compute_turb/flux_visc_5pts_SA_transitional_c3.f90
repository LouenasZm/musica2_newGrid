!==============================================================================
subroutine flux_visc_5pts_SA_transition_c3
!==============================================================================
  !> Compute derivatives of viscous fluxes (5-point stencil - order 4)
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
           nu = visc(i,j,k)/rho_n(i,j,k)
           mu    = visc(i,j,k)
           ! contravariant grad(Nutilde)
           dnutilc_ksi = dnutilx(i,j,k)*ksi_x_v(i,j,k)+dnutily(i,j,k)*ksi_y_v(i,j,k)+dnutilz(i,j,k)*ksi_z_v(i,j,k)
           dnutilc_eta = dnutilx(i,j,k)*eta_x_v(i,j,k)+dnutily(i,j,k)*eta_y_v(i,j,k)+dnutilz(i,j,k)*eta_z_v(i,j,k)
           dnutilc_eta = dnutilx(i,j,k)*phi_x_v(i,j,k)+dnutily(i,j,k)*phi_y_v(i,j,k)+dnutilz(i,j,k)*phi_z_v(i,j,k)

           dgammac_ksi = dgammax(i,j,k)*ksi_x_v(i,j,k)+dgammay(i,j,k)*ksi_y_v(i,j,k)+dgammaz(i,j,k)*ksi_z_v(i,j,k)
           dgammac_eta = dgammax(i,j,k)*eta_x_v(i,j,k)+dgammay(i,j,k)*eta_y_v(i,j,k)+dgammaz(i,j,k)*eta_z_v(i,j,k)
           dgammac_phi = dgammax(i,j,k)*phi_x_v(i,j,k)+dgammay(i,j,k)*phi_y_v(i,j,k)+dgammaz(i,j,k)*phi_z_v(i,j,k)

           dre_thetac_ksi = dre_thetax(i,j,k)*ksi_x_v(i,j,k)+dre_thetay(i,j,k)*ksi_y_v(i,j,k)+dre_thetaz(i,j,k)*ksi_z_v(i,j,k)
           dre_thetac_eta = dre_thetax(i,j,k)*eta_x_v(i,j,k)+dre_thetay(i,j,k)*eta_y_v(i,j,k)+dre_thetaz(i,j,k)*eta_z_v(i,j,k)
           dre_thetac_phi = dre_thetax(i,j,k)*phi_x_v(i,j,k)+dre_thetay(i,j,k)*phi_y_v(i,j,k)+dre_thetaz(i,j,k)*phi_z_v(i,j,k)

           ! viscous fluxes along ksi
           Fnutil(i,j,k)    = -(nu+nutil_n(i,j,k))/sig*dnutilc_ksi
           Fgamma(i,j,k)    = -(nu+nut(i,j,k))*dgammac_ksi
           Fre_theta(i,j,k) = -sigma_theta_t*(nu+nut(i,j,k))*dre_thetac_ksi

           ! viscous fluxes along eta
           Gnutil(i,j,k)    = -(nu+nutil_n(i,j,k))/sig*dnutilc_eta
           Ggamma(i,j,k)    = -(mu+mut(i,j,k))*dgammac_eta
           Gre_theta(i,j,k) = -sigma_theta_t*(mu+mut(i,j,k))*dre_thetac_eta
           
           ! viscous fluxes along phi
           Hnutil(i,j,k) = -(nu+nutil_n(i,j,k))/sig*dnutilc_phi
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
            Knutil(i,j,k) = (a04(1)*Fnutil(i  ,j,k)+a04(2)*Fnutil(i+1,j,k) &
                              +a04(3)*Fnutil(i+2,j,k)+a04(4)*Fnutil(i+3,j,k) &
                              +a04(5)*Fnutil(i+4,j,k) )*ijacob3_v(i,j,k) + Knutil(i,j,k)
            Kgamma(i,j,k) = (a04(1)*Fgamma(i  ,j,k)+a04(2)*Fgamma(i+1,j,k) &
                              +a04(3)*Fgamma(i+2,j,k)+a04(4)*Fgamma(i+3,j,k) &
                              +a04(5)*Fgamma(i+4,j,k) )*ijacob3_v(i,j,k) + Kgamma(i,j,k)
            Kre_theta(i,j,k) = (a04(1)*Fre_theta(i  ,j,k)+a04(2)*Fre_theta(i+1,j,k) &
                              +a04(3)*Fre_theta(i+2,j,k)+a04(4)*Fre_theta(i+3,j,k) &
                              +a04(5)*Fre_theta(i+4,j,k) )*ijacob3_v(i,j,k) + Kre_theta(i,j,k)
        enddo
     enddo
  endif
  if (is_bc_1pt(1,1)) then
     i=2
     do k=ndz_v_r,nfz_v_r
        do j=ndy_v_r,nfy_v_r
            Knutil(i,j,k) = (a13(1)*Fnutil(i-1,j,k)+a13(2)*Fnutil(i  ,j,k) &
                              +a13(3)*Fnutil(i+1,j,k)+a13(4)*Fnutil(i+2,j,k) &
                              +a13(5)*Fnutil(i+3,j,k) )*ijacob3_v(i,j,k) + Knutil(i,j,k)
            Kgamma(i,j,k) = (a13(1)*Fgamma(i-1,j,k)+a13(2)*Fgamma(i  ,j,k) &
                              +a13(3)*Fgamma(i+1,j,k)+a13(4)*Fgamma(i+2,j,k) &
                              +a13(5)*Fgamma(i+3,j,k) )*ijacob3_v(i,j,k) + Kgamma(i,j,k)
            Kre_theta(i,j,k) = (a13(1)*Fre_theta(i-1,j,k)+a13(2)*Fre_theta(i  ,j,k) &
                              +a13(3)*Fre_theta(i+1,j,k)+a13(4)*Fre_theta(i+2,j,k) &
                              +a13(5)*Fre_theta(i+3,j,k) )*ijacob3_v(i,j,k) + Kre_theta(i,j,k)
        enddo
     enddo
  endif

  do k=ndz_v_r,nfz_v_r
     do j=ndy_v_r,nfy_v_r
        do i=ndx_vi_r,nfx_vi_r
            Knutil(i,j,k) = (a5(1)*(Fnutil(i+1,j,k)-Fnutil(i-1,j,k) ) &
                              +a5(2)*(Fnutil(i+2,j,k)-Fnutil(i-2,j,k)))*ijacob3_v(i,j,k) + Knutil(i,j,k)
            Kgamma(i,j,k) = (a5(1)*(Fgamma(i+1,j,k)-Fgamma(i-1,j,k) ) &
                              +a5(2)*(Fgamma(i+2,j,k)-Fgamma(i-2,j,k)))*ijacob3_v(i,j,k) + Kgamma(i,j,k)
            Kre_theta(i,j,k) = (a5(1)*(Fre_theta(i+1,j,k)-Fre_theta(i-1,j,k) ) &
                              +a5(2)*(Fre_theta(i+2,j,k)-Fre_theta(i-2,j,k)))*ijacob3_v(i,j,k) + Kre_theta(i,j,k)
        enddo
     enddo
  enddo

  if (is_bc_1pt(1,2)) then
     i=nx-1
     do k=ndz_v_r,nfz_v_r
        do j=ndy_v_r,nfy_v_r
            Knutil(i,j,k) = (a31(1)*Fnutil(i+1,j,k)+a31(2)*Fnutil(i,j,k)   &  
                              +a31(3)*Fnutil(i-1,j,k)+a31(4)*Fnutil(i-2,j,k) &
                              +a31(5)*Fnutil(i-3,j,k))*ijacob3_v(i,j,k) + Knutil(i,j,k)
            Kgamma(i,j,k) = (a31(1)*Fgamma(i+1,j,k)+a31(2)*Fgamma(i,j,k)   &  
                              +a31(3)*Fgamma(i-1,j,k)+a31(4)*Fgamma(i-2,j,k) &
                              +a31(5)*Fgamma(i-3,j,k))*ijacob3_v(i,j,k) + Kgamma(i,j,k)
            Kre_theta(i,j,k) = (a31(1)*Fre_theta(i+1,j,k)+a31(2)*Fre_theta(i,j,k)   &  
                              +a31(3)*Fre_theta(i-1,j,k)+a31(4)*Fre_theta(i-2,j,k) &
                              +a31(5)*Fre_theta(i-3,j,k))*ijacob3_v(i,j,k) + Kre_theta(i,j,k)
        enddo
     enddo
  endif
  if (is_bc_wall(1,2)) then
     i=nx
     do k=ndz_v_r,nfz_v_r
        do j=ndy_v_r,nfy_v_r
            Knutil(i,j,k) = (a40(1)*Fnutil(i  ,j,k)+a40(2)*Fnutil(i-1,j,k) &
                              +a40(3)*Fnutil(i-2,j,k)+a40(4)*Fnutil(i-3,j,k) &
                              +a40(5)*Fnutil(i-4,j,k))*ijacob3_v(i,j,k) + Knutil(i,j,k)
            Kgamma(i,j,k) = (a40(1)*Fgamma(i  ,j,k)+a40(2)*Fgamma(i-1,j,k) &
                              +a40(3)*Fgamma(i-2,j,k)+a40(4)*Fgamma(i-3,j,k) &
                              +a40(5)*Fgamma(i-4,j,k))*ijacob3_v(i,j,k) + Kgamma(i,j,k)
            Kre_theta(i,j,k) = (a40(1)*Fre_theta(i  ,j,k)+a40(2)*Fre_theta(i-1,j,k) &
                              +a40(3)*Fre_theta(i-2,j,k)+a40(4)*Fre_theta(i-3,j,k) &
                              +a40(5)*Fre_theta(i-4,j,k))*ijacob3_v(i,j,k) + Kre_theta(i,j,k)
        enddo
     enddo
  endif

  ! Compute derivatives of viscous fluxes along eta
  ! ===============================================
  if (is_bc_wall(2,1)) then
     j=1
     do k=ndz_v_r,nfz_v_r
        do i=ndx_v_r,nfx_v_r
            Knutil(i,j,k) = (a04(1)*Gnutil(i,j  ,k)+a04(2)*Gnutil(i,j+1,k) &
                              +a04(3)*Gnutil(i,j+2,k)+a04(4)*Gnutil(i,j+3,k) &
                              +a04(5)*Gnutil(i,j+4,k))*ijacob3_v(i,j,k) + Knutil(i,j,k)
            Kgamma(i,j,k) = (a04(1)*Ggamma(i,j  ,k)+a04(2)*Ggamma(i,j+1,k) &
                              +a04(3)*Ggamma(i,j+2,k)+a04(4)*Ggamma(i,j+3,k) &
                              +a04(5)*Ggamma(i,j+4,k))*ijacob3_v(i,j,k) + Kgamma(i,j,k)
            Kre_theta(i,j,k) = (a04(1)*Gre_theta(i,j  ,k)+a04(2)*Gre_theta(i,j+1,k) &
                              +a04(3)*Gre_theta(i,j+2,k)+a04(4)*Gre_theta(i,j+3,k) &
                              +a04(5)*Gre_theta(i,j+4,k))*ijacob3_v(i,j,k) + Kre_theta(i,j,k)
        enddo
     enddo
  endif
  if (is_bc_1pt(2,1)) then
     j=2
     do k=ndz_v_r,nfz_v_r
        do i=ndx_v_r,nfx_v_r
            Knutil(i,j,k) = (a13(1)*Gnutil(i,j-1,k)+a13(2)*Gnutil(i,j,k)   &
                              +a13(3)*Gnutil(i,j+1,k)+a13(4)*Gnutil(i,j+2,k) &
                              +a13(5)*Gnutil(i,j+3,k))*ijacob3_v(i,j,k) + Knutil(i,j,k)
            Kgamma(i,j,k) = (a13(1)*Ggamma(i,j-1,k)+a13(2)*Ggamma(i,j,k)   &
                              +a13(3)*Ggamma(i,j+1,k)+a13(4)*Ggamma(i,j+2,k) &
                              +a13(5)*Ggamma(i,j+3,k))*ijacob3_v(i,j,k) + Kgamma(i,j,k)
            Kre_theta(i,j,k) = (a13(1)*Gre_theta(i,j-1,k)+a13(2)*Gre_theta(i,j,k)   &
                              +a13(3)*Gre_theta(i,j+1,k)+a13(4)*Gre_theta(i,j+2,k) &
                              +a13(5)*Gre_theta(i,j+3,k))*ijacob3_v(i,j,k) + Kre_theta(i,j,k)
        enddo
     enddo
  endif

  do k=ndz_v_r,nfz_v_r
     do j=ndy_vi_r,nfy_vi_r
        do i=ndx_v_r,nfx_v_r
            Knutil(i,j,k) = (a5(1)*(Gnutil(i,j+1,k)-Gnutil(i,j-1,k)) &
                              +a5(2)*(Gnutil(i,j+2,k)-Gnutil(i,j-2,k)) )*ijacob3_v(i,j,k) + Knutil(i,j,k)
            Kgamma(i,j,k) = (a5(1)*(Ggamma(i,j+1,k)-Ggamma(i,j-1,k)) &
                              +a5(2)*(Ggamma(i,j+2,k)-Ggamma(i,j-2,k)) )*ijacob3_v(i,j,k) + Kgamma(i,j,k)
            Kre_theta(i,j,k) = (a5(1)*(Gre_theta(i,j+1,k)-Gre_theta(i,j-1,k)) &
                              +a5(2)*(Gre_theta(i,j+2,k)-Gre_theta(i,j-2,k)) )*ijacob3_v(i,j,k) + Kre_theta(i,j,k)
        enddo
     enddo
  enddo
  
  if (is_bc_1pt(2,2)) then
     j=ny-1
     do k=ndz_v_r,nfz_v_r
        do i=ndx_v_r,nfx_v_r
            Knutil(i,j,k) = (a31(1)*Gnutil(i,j+1,k)+a31(2)*Gnutil(i,j,k)   &
                              +a31(3)*Gnutil(i,j-1,k)+a31(4)*Gnutil(i,j-2,k) &
                              +a31(5)*Gnutil(i,j-3,k))*ijacob3_v(i,j,k) + Knutil(i,j,k)
            Kgamma(i,j,k) = (a31(1)*Ggamma(i,j+1,k)+a31(2)*Ggamma(i,j,k)   &
                              +a31(3)*Ggamma(i,j-1,k)+a31(4)*Ggamma(i,j-2,k) &
                              +a31(5)*Ggamma(i,j-3,k))*ijacob3_v(i,j,k) + Kgamma(i,j,k)
            Kre_theta(i,j,k) = (a31(1)*Gre_theta(i,j+1,k)+a31(2)*Gre_theta(i,j,k)   &
                              +a31(3)*Gre_theta(i,j-1,k)+a31(4)*Gre_theta(i,j-2,k) &
                              +a31(5)*Gre_theta(i,j-3,k))*ijacob3_v(i,j,k) + Kre_theta(i,j,k)
        enddo
     enddo
  endif
  if (is_bc_wall(2,2)) then
     j=ny
     do k=ndz_v_r,nfz_v_r
        do i=ndx_v_r,nfx_v_r
            Knutil(i,j,k) = (a40(1)*Gnutil(i,j  ,k)+a40(2)*Gnutil(i,j-1,k) &
                              +a40(3)*Gnutil(i,j-2,k)+a40(4)*Gnutil(i,j-3,k) &
                              +a40(5)*Gnutil(i,j-4,k))*ijacob3_v(i,j,k) + Knutil(i,j,k)
            Kgamma(i,j,k) = (a40(1)*Ggamma(i,j  ,k)+a40(2)*Ggamma(i,j-1,k) &
                              +a40(3)*Ggamma(i,j-2,k)+a40(4)*Ggamma(i,j-3,k) &
                              +a40(5)*Ggamma(i,j-4,k))*ijacob3_v(i,j,k) + Kgamma(i,j,k)
            Kre_theta(i,j,k) = (a40(1)*Gre_theta(i,j  ,k)+a40(2)*Gre_theta(i,j-1,k) &
                              +a40(3)*Gre_theta(i,j-2,k)+a40(4)*Gre_theta(i,j-3,k) &
                              +a40(5)*Gre_theta(i,j-4,k))*ijacob3_v(i,j,k) + Kre_theta(i,j,k)
        enddo
     enddo
  endif

  ! Compute derivatives of viscous fluxes along phi
  ! ===============================================
  if (is_bc_wall(3,1)) then
     k=1
     do j=ndy_v_r,nfy_v_r
        do i=ndx_v_r,nfx_v_r
            Knutil(i,j,k) = (a04(1)*Hnutil(i,j,k  )+a04(2)*Hnutil(i,j,k+1) &
                              +a04(3)*Hnutil(i,j,k+2)+a04(4)*Hnutil(i,j,k+3) &
                              +a04(5)*Hnutil(i,j,k+4))*ijacob3_v(i,j,k) + Knutil(i,j,k)
            Kgamma(i,j,k) = (a04(1)*Hgamma(i,j,k  )+a04(2)*Hgamma(i,j,k+1) &
                              +a04(3)*Hgamma(i,j,k+2)+a04(4)*Hgamma(i,j,k+3) &
                              +a04(5)*Hgamma(i,j,k+4))*ijacob3_v(i,j,k) + Kgamma(i,j,k)
            Kre_theta(i,j,k) = (a04(1)*Hre_theta(i,j,k  )+a04(2)*Hre_theta(i,j,k+1) &
                              +a04(3)*Hre_theta(i,j,k+2)+a04(4)*Hre_theta(i,j,k+3) &
                              +a04(5)*Hre_theta(i,j,k+4))*ijacob3_v(i,j,k) + Kre_theta(i,j,k)
        enddo
     enddo
  endif
  if (is_bc_1pt(3,1)) then
     k=2
     do j=ndy_v_r,nfy_v_r
        do i=ndx_v_r,nfx_v_r
            Knutil(i,j,k) = (a13(1)*Hnutil(i,j,k-1)+a13(2)*Hnutil(i,j,k)   &
                              +a13(3)*Hnutil(i,j,k+1)+a13(4)*Hnutil(i,j,k+2) &
                              +a13(5)*Hnutil(i,j,k+3))*ijacob3_v(i,j,k) + Knutil(i,j,k)
            Kgamma(i,j,k) = (a13(1)*Hgamma(i,j,k-1)+a13(2)*Hgamma(i,j,k)   &
                              +a13(3)*Hgamma(i,j,k+1)+a13(4)*Hgamma(i,j,k+2) &
                              +a13(5)*Hgamma(i,j,k+3))*ijacob3_v(i,j,k) + Kgamma(i,j,k)
            Kre_theta(i,j,k) = (a13(1)*Hre_theta(i,j,k-1)+a13(2)*Hre_theta(i,j,k)   &
                              +a13(3)*Hre_theta(i,j,k+1)+a13(4)*Hre_theta(i,j,k+2) &
                              +a13(5)*Hre_theta(i,j,k+3))*ijacob3_v(i,j,k) + Kre_theta(i,j,k)
        enddo
     enddo
  endif

  do k=ndz_vi_r,nfz_vi_r
     do j=ndy_v_r,nfy_v_r
        do i=ndx_v_r,nfx_v_r
            Knutil(i,j,k) = (a5(1)*(Hnutil(i,j ,k+1)-Hnutil(i,j,k-1)) &
                              +a5(2)*( Hnutil(i,j,k+2)-Hnutil(i,j,k-2)))*ijacob3_v(i,j,k) + Knutil(i,j,k)
            Kgamma(i,j,k) = (a5(1)*(Hgamma(i,j ,k+1)-Hgamma(i,j,k-1)) &
                              +a5(2)*( Hgamma(i,j,k+2)-Hgamma(i,j,k-2)))*ijacob3_v(i,j,k) + Kgamma(i,j,k)
            Kre_theta(i,j,k) = (a5(1)*(Hre_theta(i,j ,k+1)-Hre_theta(i,j,k-1)) &
                              +a5(2)*( Hre_theta(i,j,k+2)-Hre_theta(i,j,k-2)))*ijacob3_v(i,j,k) + Kre_theta(i,j,k)
        enddo
     enddo
  enddo

  if (is_bc_1pt(3,2)) then
     k=nz-1
     do j=ndy_v_r,nfy_v_r
        do i=ndx_v_r,nfx_v_r
            Knutil(i,j,k) = (a31(1)*Hnutil(i,j,k+1)+a31(2)*Hnutil(i,j,k)   &
                              +a31(3)*Hnutil(i,j,k-1)+a31(4)*Hnutil(i,j,k-2) &
                              +a31(5)*Hnutil(i,j,k-3))*ijacob3_v(i,j,k) + Knutil(i,j,k)
            Kgamma(i,j,k) = (a31(1)*Hgamma(i,j,k+1)+a31(2)*Hgamma(i,j,k)   &
                              +a31(3)*Hgamma(i,j,k-1)+a31(4)*Hgamma(i,j,k-2) &
                              +a31(5)*Hgamma(i,j,k-3))*ijacob3_v(i,j,k) + Kgamma(i,j,k)
            Kre_theta(i,j,k) = (a31(1)*Hre_theta(i,j,k+1)+a31(2)*Hre_theta(i,j,k)   &
                              +a31(3)*Hre_theta(i,j,k-1)+a31(4)*Hre_theta(i,j,k-2) &
                              +a31(5)*Hre_theta(i,j,k-3))*ijacob3_v(i,j,k) + Kre_theta(i,j,k)
        enddo
     enddo
  endif
  if (is_bc_wall(3,2)) then
     k=nz
     do j=ndy_v_r,nfy_v_r
        do i=ndx_v_r,nfx_v_r
            Knutil(i,j,k) = (a40(1)*Hnutil(i,j,k)  +a40(2)*Hnutil(i,j,k-1) &
                              +a40(3)*Hnutil(i,j,k-2)+a40(4)*Hnutil(i,j,k-3) &
                              +a40(5)*Hnutil(i,j,k-4))*ijacob3_v(i,j,k) + Knutil(i,j,k)
            Kgamma(i,j,k) = (a40(1)*Hgamma(i,j,k)  +a40(2)*Hgamma(i,j,k-1) &
                              +a40(3)*Hgamma(i,j,k-2)+a40(4)*Hgamma(i,j,k-3) &
                              +a40(5)*Hgamma(i,j,k-4))*ijacob3_v(i,j,k) + Kgamma(i,j,k)
            Kre_theta(i,j,k) = (a40(1)*Hre_theta(i,j,k)  +a40(2)*Hre_theta(i,j,k-1) &
                              +a40(3)*Hre_theta(i,j,k-2)+a40(4)*Hre_theta(i,j,k-3) &
                              +a40(5)*Hre_theta(i,j,k-4))*ijacob3_v(i,j,k) + Kre_theta(i,j,k)
        enddo
     enddo
  endif

end subroutine flux_visc_5pts_SA_transition_c3

