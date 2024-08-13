!==============================================================================
subroutine source_SA_transition
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
  real(wp) :: nu,D11,D12,D13,D21,D22,D23,D31,D32,D33,D_om,ksi_nu,g_n, D_mag
  real(wp) :: cw36,inv_6,inv_kap2, cnu13
  real(wp), parameter :: TWO_THIRDS=2.0_wp/3.0_wp, large_St1=5.0E3_wp
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

            ! 
            ! Kinematic & turbulent viscosities
            nu  = visc(i,j,k)/rho_n(i,j,k)
            if (intermittency(i,j,k) .lt. 0.0_wp)then
                intermittency(i,j,k) = 0.0_wp 
            elseif(intermittency(i,j,k) .ge. 1.0_wp)then
                intermittency(i,j,k) = 1.0_wp 
            endif
            ! ========================================
            ! Source terms in the equation of Re_theta
            !
            ! Vorticity tensor component
            om12 = 0.5_wp*(duy(i,j,k)-dvx(i,j,k))
            om13 = 0.5_wp*(duz(i,j,k)-dwx(i,j,k))
            om23 = 0.5_wp*(dvz(i,j,k)-dwy(i,j,k))
            magn_om = sqrt(4.0_wp*(om12**2+om13**2+om23**2))

            ! U magnitude of velocity:
            magn_u  = sqrt( uu(i,j,k)**2 + vv(i,j,k)**2 + ww(i,j,k)**2 )
            magn_u2 = uu(i,j,k)**2 + vv(i,j,k)**2 + ww(i,j,k)**2

            ! T:
            T_stheta = 500.0_wp * nu / magn_u2
            ! delta:
            delta_stheta = 375.0_wp *magn_om*reynolds_theta(i,j,k)*lengthscale(i,j,k)/magn_u2
            if (delta_stheta .eq. 0.0_wp) then 
                delta_stheta = 1e-10
            endif
            ! F_theta:
            F_theta = min( max( exp(-(lengthscale(i,j,k)/delta_stheta)**4 ), 1.0_wp-((intermittency(i,j,k) - 1.0_wp/ce2)/(1.0_wp - 1.0_wp/ce2))**2 ), 1.0_wp)

            ! Re_theta_t in production term of reynolds_theta:
            du_ds           = (uu(i,j,k)*uu(i,j,k)*dux(i,j,k) + vv(i,j,k)*uu(i,j,k)*dvx(i,j,k) + uu(i,j,k)*vv(i,j,k)*duy(i,j,k) + vv(i,j,k)*vv(i,j,k)*dvy(i,j,k))/magn_u2
            if(.not. is_2D)then 
                du_ds = du_ds + (ww(i,j,k)*uu(i,j,k)*dwx(i,j,k) + ww(i,j,k)*vv(i,j,k)*dwy(i,j,k) + &
                                uu(i,j,k)*ww(i,j,k)*duz(i,j,k) + vv(i,j,k)*ww(i,j,k)*dvz(i,j,k) + ww(i,j,k)*ww(i,j,k)*dwz(i,j,k))/magn_u2
            endif
            
            theta           = reynolds_theta(i,j,k)*nu/magn_u
            lambda_theta    = theta**2*du_ds/nu 
            ! Limiter for Lambda_theta:
            if (lambda_theta .lt. -0.1_wp)then
                lambda_theta = -0.1_wp 
            elseif(lambda_theta .gt. 0.1_wp)then
                lambda_theta = 0.1_wp 
            endif
            ! F_lambda_t
            if (lambda_theta .le. 0.0_wp)then
                f_lambda_t  = 1.0_wp + ( 12.986_wp*lambda_theta + 123.66_wp*lambda_theta**2 + 405.689_wp*lambda_theta**3 )*exp(-(tu_inlet/1.5_wp)**1.5) 
            else
                f_lambda_t  = 1.0_wp + 0.275_wp*(1.0_wp - exp(-35.0_wp*lambda_theta))*exp(-tu_inlet/0.5_wp)
            endif
            ! Re_theta eq:
            if(tu_inlet .le. 1.3)then
                Re_theta_t  = (1173.51_wp - 589.428_wp*tu_inlet + 0.2196_wp/tu_inlet**2)*f_lambda_t
            else
                Re_theta_t  = 331.5_wp*(tu_inlet-0.5668_wp)**(-0.671_wp)*f_lambda_t
            endif
            ! Limiter for Re_theta_t:
            if (Re_theta_t .lt. 20.0_wp) Re_theta_t = 20.0_wp 
            ! Production of re theta:
            Stheta(i,j,k) = ctheta_t/T_stheta*(Re_theta_t - reynolds_theta(i,j,k))*(1 - F_theta)

            ! =====================================
            ! Source terms in the equation of Gamma
            !
            re_theta_crit = min(0.615_wp*reynolds_theta(i,j,k) + 61.5_wp, reynolds_theta(i,j,k))
            if (re_theta_crit == 0.0_wp)then
                re_theta_crit = 1e-10
            endif
            ! Mean strain rate tensor:
            D11 = dux(i,j,k)
            D22 = dvy(i,j,k)
            D33 = dwz(i,j,k)
            D12 = 0.5_wp*(duy(i,j,k)+dvx(i,j,k))
            D13 = 0.5_wp*(duz(i,j,k)+dwx(i,j,k))
            D23 = 0.5_wp*(dvx(i,j,k)+dwy(i,j,k))
            D_mag   = sqrt(2.0_wp*(D11**2+D22**2+D33**2 + &
                        2.0_wp*(D12**2+D13**2+D23**2)))
            D_om = D_mag - magn_om

            ! Reynolds_nu and reynolds T:
            re_nu   = D_mag*lengthscale(i,j,k)**2/nu 
            re_turb = nut(i,j,k)/nu
            ! F_turb:
            F_turb      = exp( - (0.25_wp*re_turb)**4 )
            ! F_onsets:
            F_onset1    = re_nu/(2.193_wp*re_theta_crit)
            F_onset2    = min( max( F_onset1, F_onset1**4 ), 4.0_wp )
            F_onset3    = max( (2.0_wp - (re_turb/2.5_wp)**3 ), 0.0_wp)
            F_onset     = max( ( F_onset2 - F_onset3 ), 0.0_wp )

            ! Production and destruction:
            F_length        = min( exp( 7.168_wp - 0.01173_wp*reynolds_theta(i,j,k) ) + 0.5_wp, 300.0_wp )
            p_gamma         = ca1*D_mag*(intermittency(i,j,k)*F_onset)**(0.5_wp)*(1.0_wp - intermittency(i,j,k))*F_length
            d_gamma         = ca2*magn_om*intermittency(i,j,k)*F_turb*(ce2*intermittency(i,j,k) - 1.0_wp)
            Sgamma(i,j,k)   = p_gamma - d_gamma

             ! =========================
            ! Source terms in the equation of Nutil
            !
            ! Modified vorticity
            S    = magn_om + min( 0.0_wp, D_om )
            Stil = nutil(i,j,k)*inv_kap2/lengthscale(i,j,k)**2*fnu2
            Stil = S+Stil

            ! Model parameters
            
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
            
            
            ! Take into consideration transition mechanisms:
            F_reattach          = exp( - (re_turb/20)**4 )
            intermittency_sep   = min( 2.0_wp*max( 0.0_wp, (re_nu/(3.235_wp*re_theta_crit) - 1.0_wp))*F_reattach, 2.0_wp)*F_theta
            
            intermittency_eff(i,j,k)   = max(intermittency(i,j,k), intermittency_sep)

            if (intermittency_eff(i,j,k) .lt. 0.0_wp)then
                intermittency_eff(i,j,k) = 0.0_wp 
            elseif(intermittency_eff(i,j,k) .ge. 1.0_wp)then
                intermittency_eff(i,j,k) = 1.0_wp 
            endif

            ! Production and destruction terms of nutil:
            St1 = intermittency_eff(i,j,k)*St1

            dnutili2 = dnutilx(i,j,k)**2+dnutily(i,j,k)**2+dnutilz(i,j,k)**2
            St3      = cb2/sig*dnutili2
            
            ! Total source term
            Sterm(i,j,k) = St1+St2+St3

        enddo
     enddo
  enddo
 
end subroutine source_SA_transition