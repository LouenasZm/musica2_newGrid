!==============================================================================
subroutine source_SA_criv
!==============================================================================
  !> Compute Spalart-Allmaras source terms -> nutilde transport equation RHS 
!==============================================================================
  use mod_flow
  use mod_rans
  use mod_mpi
  use mod_wall_dist
  use mod_turb_model_length_scale
  implicit none
  ! ---------------------------------------------------------------------------
  integer  :: i,j,k
  real(wp) :: nu,D11,D12,D13,D22,D23,D33,D_om,ksi_nu,g_n
  real(wp) :: cw36,inv_6,cnu13,inv_kap2
  real(wp), parameter :: TWO_THIRDS=2.0_wp/3.0_wp,large_St1=5.0E3_wp
  ! ---------------------------------------------------------------------------

  ! Model constants & parameters
  ! ============================
  sig  = TWO_THIRDS
  kap  = 0.41_wp
  cb1  = 0.1355_wp
  cb2  = 0.622_wp
  cnu1 = 7.1_wp
  cnu13=cnu1**3
  cnu2 = 0.7_wp
  cnu3 = 0.9_wp
  inv_6 = 1.0_wp/6.0_wp
  inv_kap2 = 1.0_wp/kap**2
  cw1  = cb1*inv_kap2+(1.0_wp+cb2)/sig
  cw2  = 0.3_wp
  cw3  = 2.0_wp
  
  do k=ndz_s_r,nfz_s_r 
     do j=ndy_s_r,nfy_s_r
        do i=ndx_s_r,nfx_s_r
           ! Kinematic & turbulent viscosities
           nu = visc(i,j,k)/rho_n(i,j,k)
         
           ! Vorticity tensor component
           om12 = 0.5_wp*(duy(i,j,k)-dvx(i,j,k))
           om13 = 0.5_wp*(duz(i,j,k)-dwx(i,j,k))
           om23 = 0.5_wp*(dvz(i,j,k)-dwy(i,j,k))
           magn_om = sqrt(4.0_wp*(om12**2+om13**2+om23**2))

           ! SA-R correction
           D11 = dux(i,j,k)
           D22 = dvy(i,j,k)
           D33 = dwz(i,j,k)
           D12 = 0.5_wp*(duy(i,j,k)+dvx(i,j,k))
           D13 = 0.5_wp*(duz(i,j,k)+dwx(i,j,k))
           D23 = 0.5_wp*(dvx(i,j,k)+dwy(i,j,k))
           D_om = sqrt(2.0_wp*(D11**2+D22**2+D33**2 + &
                       2.0_wp*(D12**2+D13**2+D23**2))) - magn_om


           S = magn_om + 2.0_wp*min(0.0_wp,D_om)

           ! Modified vorticity
           Stil = nutil(i,j,k)*inv_kap2/lengthscale(i,j,k)**2*fnu2
           Stil = S+Stil

           ! Model parameters
           khi  = nutil(i,j,k)/nu
           fnu1 = khi**3/(khi**3+cnu13)
           fnu2 = 1.0_wp-khi/(1.0_wp+khi*fnu1)

           ! Crivellini et. al. modif pour nutil<0
           if (khi.ge.0.0_wp) then
              ! Model parameters evaluated only in this case
              r = min(nutil(i,j,k)/Stil*inv_kap2/lengthscale(i,j,k)**2,10.0_wp)
              if (r.le.0.0_wp) r = 10.0_wp
              g_sa = r + cw2*(r**6-r)
              cw36=cw3**6
              fw = g_sa*((1.0_wp+cw36)/(g_sa**6+cw36))**inv_6
              ksi_nu = nu*(1.0_wp+khi)

              ! Source term 1 -> Production (same as before, just prevents negative values)
              St1 = cb1*nutil(i,j,k)**2/r*inv_kap2/lengthscale(i,j,k)**2

              ! Source term 2 -> Dissipation
              St2 = -cw1*fw*(nutil(i,j,k)/lengthscale(i,j,k))**2
            else
              ksi_nu = nu*(1.0_wp+khi+0.5_wp*khi**2)
              g_n = 1.0_wp-10.0_wp**3*khi**2/(1.0_wp+khi**2)

              ! Source term 1 -> Production
              St1 = cb1*S*nutil(i,j,k)*g_n

              ! Source term 2 -> Dissipation
              St2 = cw1*(nutil(i,j,k)/lengthscale(i,j,k))**2
            endif

           if (St1.gt.large_St1) St1 = large_St1
          
           ! Source term 3 -> Diffusion
           dnutili2 = dnutilx(i,j,k)**2+dnutily(i,j,k)**2+dnutilz(i,j,k)**2
           St3      = cb2/sig*dnutili2

           ! Total source term
           Sterm(i,j,k) = St1+St2+St3

           ! Plot source terms
!            uvar(i,j,k,3) = St1
!            uvar(i,j,k,4) = St2
!            uvar(i,j,k,5) = St3
!            uvar(i,j,k,6) = Sterm(i,j,k)
        enddo
     enddo
  enddo
 
end subroutine source_SA_criv
