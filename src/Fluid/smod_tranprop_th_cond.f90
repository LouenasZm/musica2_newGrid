!=================================================================================
submodule (mod_tranprop) smod_tranprop_th_cond
!=================================================================================
  !> author: XG
  !> date: January 2024
  !> Submodule for routines computing thermal conductivity
!=================================================================================

contains

  !===============================================================================
  module function thconductivity_Prcost(mu,T,ro)
  !===============================================================================
    !> Constant Prandtl number law for thermal conductivity
    !> (only for pfg gases)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: mu,T,ro
    real(wp) :: thconductivity_Prcost ! OUTPUT
    ! ----------------------------------------------------------------------------
    
    thconductivity_Prcost= var2*mu
    
  end function thconductivity_Prcost

  !===============================================================================
  module function thconductivity_refprop(mu,T,ro)
  !===============================================================================
    !> Thermal conductivity law of RefProp
  !===============================================================================
    use mod_ineos_ref
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: mu,T,ro
    real(wp) :: thconductivity_refprop ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: d,eta,tcx
    ! ----------------------------------------------------------------------------
    
    ! Calculate viscosity (eta) and thermal conductivity (tcx)
    ! ========================================================
    d=ro/pmol
    call TRNPRP(t,d,1.0,eta,tcx,ierr,herr)
    thconductivity_refprop=tcx
    
  end function thconductivity_refprop

  !===============================================================================
  module function thconductivity_chung(mu,T,ro)
  !===============================================================================
    !> Chung-Lee-Starling law for thermal conductivity
    !> (for real gases)
    !===============================================================================
    use mod_eos  ! for: cvcalc_tro   
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: mu,T,ro
    real(wp) :: thconductivity_chung ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: lambda
    real(wp) :: t_st,omegav,Yv
    real(wp) :: alfav,betav,zetav,psiv
    real(wp) :: G1,H2,eta_0,lambda_0,lambda_k,lambda_p
    ! ----------------------------------------------------------------------------

    ! Compute viscosity first
    ! =======================
    t_st= 1.2593_wp*T/Tc
    omegav = avisc/t_st**(bvisc) + cvisc/exp(dvisc*t_st)   &
           + evisc/exp(fvisc*t_st) &
           + gvisc*t_st**(bvisc)*sin(svisc*t_st**(wvisc) - hvisc)

    eta_0 = k_chung*sqrt(T)/omegav

    ! Yv = ro*vmolc/6.0_wp
    Yv= ro/(6.0_wp*roc0)

    ! Compute thermal conductivity
    ! ============================

    ! rotational coefficient
    alfav= cvcalc_tro(T,ro)/rg - 1.5_wp
    ! diffusion term empirically linked to the acentric factor
    betav= 0.7862_wp - 0.7109_wp*om + 1.3168_wp*om**2
    ! number of collisions required to interchange a quantum of rotational energy with translational energy
    zetav= 2.0_wp + 10.5_wp*(T/Tc)**2
    ! modified Eucken-type correlation based on kinetic theory extended to polyatomic gases
    psiv = 1.0_wp + alfav*(0.215_wp + 0.28288_wp*alfav - 1.061_wp*betav       &
         + 0.26665_wp*zetav)/(0.6366_wp+betav*zetav+1.061_wp*alfav*betav)

    lambda_0= 7.452_wp*eta_0/pmol*psiv

    G1= (1.0_wp-0.5_wp*Yv)/(1.0_wp-Yv)**3
    H2= (BB1(1)/Yv*(1.0_wp-exp(-BB1(4)*Yv)) + BB1(2)*G1*exp(BB1(5)*Yv) + BB1(3)*G1) / &
         (BB1(1)*BB1(4) + BB1(2) + BB1(3))

    ! Dilute-Gas component
    lambda_k = lambda_0 * (1.0_wp/H2 + BB1(6)*Yv)

    ! Dense-Gas component
    lambda_p=3.039e-4_wp*sqrt(T/pmol)/vmolc**(2.0_wp/3.0_wp)*BB1(7)*Yv**2*H2

    ! Convert from cal/(cm s K) to W/m/K
    lambda= (lambda_k+lambda_p)*418.68_wp

    thconductivity_chung= lambda

  end function thconductivity_chung

end submodule smod_tranprop_th_cond
