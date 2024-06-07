!==============================================================================
subroutine source_SA_transition_algebraic
!==============================================================================
  !> Compute Spalart-Allmaras source terms -> nutilde transport equation RHS
  !  for transitional flows using Crivellini's algebraic intermittency model  
!==============================================================================
  use mod_flow
  use mod_rans
  use mod_mpi
  use mod_wall_dist
  use mod_turb_model_length_scale
  implicit none
  ! ---------------------------------------------------------------------------
  
  ! ---------------------------------------------------------------------------
  integer  :: i,j,k
  real(wp) :: nu,D11,D12,D13,D21,D22,D23,D31,D32,D33,D_om,ksi_nu,g_n
  real(wp) :: cw36,inv_6,cnu13,inv_kap2
  ! ---------------------------------------------------------------------------

  ! Model constants & parameters
  ! ============================
  sig       = TWO_THIRDS
  kap       = 0.41_wp
  ksi1      = 0.002_wp 
  ksi2      = 50.0_wp 
  cb1       = 0.1355_wp
  cb2       = 0.622_wp
  cnu1      = 7.1_wp
  cnu13     = cnu1**3
  inv_6     = 1.0_wp/6.0_wp
  inv_kap2  = 1.0_wp/kap**2
  cw1       = cb1*inv_kap2+(1.0_wp+cb2)/sig
  cw3       = 2.0_wp
  cw4       = 0.21_wp 
  cw5       = 1.5_wp 
  
  do k=ndz_s_r,nfz_s_r 
     do j=ndy_s_r,nfy_s_r
        do i=ndx_s_r,nfx_s_r
           ! Kinematic & turbulent viscosities
           nu           = visc(i,j,k)/rho_n(i,j,k)
           nut(i,j,k)   = nutil(i,j,k)*fnu1
         
           ! Vorticity tensor component
           om12     = 0.5_wp*(duy(i,j,k)-dvx(i,j,k))
           om13     = 0.5_wp*(duz(i,j,k)-dwx(i,j,k))
           om23     = 0.5_wp*(dvz(i,j,k)-dwy(i,j,k))
           magn_om  = sqrt(4.0_wp*(om12**2+om13**2+om23**2))
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
        
           ! Transitional onset: 
           re_vort  = lengthscale(i,j,k)**2*magn_om/nu 
           re_theta = re_vort/2.193_wp 
           re_crit  = 803.73_wp *( Tu_inlet + 0.6067_wp )**(-1.027_wp)
           term1    = max( re_theta - re_crit, 0.0_wp )/(ksi1*re_crit)
           term2    = max( nut(i,j,k)/(ksi2*nu) , 0.0_wp)
           gamma_cb = 1.0_wp - exp( -sqrt(term1) -sqrt(term2) )
           if(gamma_cb .lt. 0.0_wp) gamma_cb = 0.0_wp 
        
           ! Model parameters
           khi  = nutil_n(i,j,k)/nu
           fnu1 = khi**3/(khi**3+cnu13)
           fnu2 = 1.0_wp-khi/(1.0_wp+khi*fnu1)

           ! Modified vorticity
           Stil = nutil_n(i,j,k)*inv_kap2/lengthscale(i,j,k)**2*fnu2
           Stil = S+Stil
           ! Crivellini et. al. modif pour nutil<0
           if (khi.ge.0.0_wp) then
              ! Model parameters evaluated only in this case
              r = min(nutil_n(i,j,k)/Stil*inv_kap2/lengthscale(i,j,k)**2,10.0_wp)
              if (r.le.0.0_wp) r = 10.0_wp
              cw2_lre   = cw4 + cw5/(khi)
              g_sa      = r + cw2_lre*(r**6-r)
              cw36      = cw3**6
              fw        = g_sa*((1.0_wp+cw36)/(g_sa**6+cw36))**inv_6
              ksi_nu    = nu*(1.0_wp+khi)

              ! Source term 1 -> Production (same as before, just prevents negative values)
              St1 = cb1*nutil_n(i,j,k)**2/r*inv_kap2/lengthscale(i,j,k)**2

              ! Source term 2 -> Dissipation
              St2 = -cw1*fw*(nutil_n(i,j,k)/lengthscale(i,j,k))**2
            else
              ksi_nu = nu*(1.0_wp+khi+0.5_wp*khi**2)
              g_n = 1.0_wp-10.0_wp**3*khi**2/(1.0_wp+khi**2)

              ! Source term 1 -> Production
              St1 = gamma_cb*cb1*S*nutil_n(i,j,k)*g_n

              ! Source term 2 -> Dissipation
              St2 = cw1*(nutil_n(i,j,k)/lengthscale(i,j,k))**2
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
 
end subroutine source_SA_transition_algebraic