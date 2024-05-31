!=================================================================================
submodule (mod_tranprop) smod_tranprop_visc
!=================================================================================
  !> author: XG
  !> date: January 2024
  !> Submodule for routines computing dynamic viscosity
!=================================================================================

contains

  !===============================================================================
  module function viscosity_law_sutherland(T,ro)
  !===============================================================================
    !> Sutherland's law for viscosity (for perfect gas only)
    !> (rho is not used)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: viscosity_law_sutherland ! OUTPUT
    ! ----------------------------------------------------------------------------
    
    !viscosity_law_sutherland= mu0*sqrt(T/T0)*var1/(1.0_wp+S0/T)
    viscosity_law_sutherland= c0*sqrt(T**3)/(T+S0)
    
  end function viscosity_law_sutherland

  !===============================================================================
  module function viscosity_law_power(T,ro)
  !===============================================================================
    !> Power law for viscosity (exponent defined in feos_FLUID.ini)
    !> (rho is not used)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: viscosity_law_power ! OUTPUT
    ! ----------------------------------------------------------------------------
    
    viscosity_law_power= var3*T**nexp
    
  end function viscosity_law_power

  !===============================================================================
  module function viscosity_law_refprop(T,ro)
  !===============================================================================
    !> Viscosity law of RefProp
    !===============================================================================
    use mod_ineos_ref
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: viscosity_law_refprop ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: d,eta,tcx
    ! ----------------------------------------------------------------------------

    ! Calculate viscosity (eta) and thermal conductivity (tcx)
    ! ========================================================
    d=ro/pmol
    call TRNPRP(t,d,1.0,eta,tcx,ierr,herr)
    viscosity_law_refprop=eta*1e-6_wp
    
  end function viscosity_law_refprop

  !===============================================================================
  module function viscosity_law_chung(T,ro)
  !===============================================================================
    !> Chung-Lee-Starling law for viscosity (real gas)
    !> from Chung, Ajlan, Lee & Starling (Ind. Eng. Chem. Res., 1988, 27, 671-679)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: viscosity_law_chung ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: eta,T_st,omegav,Yv
    real(wp) :: G1,G2,eta_0,eta_k,eta_p
    ! ----------------------------------------------------------------------------

    ! Dilute-gas limit viscosity (Chapman-Enskog)
    ! ==========================
    ! dimensionless temperature T_st=k*T/eps
    !                           with k: Boltzmann cste & eps: potential energy param.    
    T_st=1.2593_wp*T/Tc
    
    ! reduced collision integral
    ! (empirical eq. from Neufeld et al., 1972)
    omegav = avisc/T_st**(bvisc) + cvisc/exp(dvisc*T_st)   &
           + evisc/exp(fvisc*T_st) &
           + gvisc*T_st**(bvisc)*sin(svisc*T_st**(wvisc)-hvisc)

    eta_0= k_chung*sqrt(T)/omegav

    ! Yv = ro*vmolc/6.0_wp
    Yv= ro/(6.0_wp*roc0)

    G1= (1.0_wp-0.5_wp*Yv)/(1.0_wp-Yv)**3
    G2= (AA1(1)/Yv*(1.0_wp-exp(-AA1(4)*Yv)) + AA1(2)*G1*exp(AA1(5)*Yv) + AA1(3)*G1)/ &
         (AA1(1)*AA1(4) + AA1(2) + AA1(3))
    
    ! Dilute-Gas component with density dependence
    ! ============================================
    eta_k=eta_0*(1.0_wp/G2 + AA1(6)*Yv)

    ! Dense-Gas component (correction)
    ! ===================
    eta_p=(36.344e-6_wp*sqrt(pmol*Tc)/(vmolc)**(2.0_wp/3.0_wp)) * &
         AA1(7)*Yv**2*G2*exp(AA1(8)+AA1(9)/T_st+AA1(10)/T_st**2)

    eta=0.1_wp*(eta_k+eta_p) ! 0.1 to convert P in Pa.s

!!$    eta=0.1_wp*eta_k
    
    viscosity_law_chung= eta

  end function viscosity_law_chung

  !===============================================================================
  module function viscosity_law_wen(T,ro)
  !===============================================================================
    !> Wen et al.'s law for viscosity (reference law of RefProp for Novec649)
    !> from Wen, Meng, Huber & Wu (J. Chem. Eng. Data, 2017, 62, 3603-3609)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: viscosity_law_wen ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: eta,T_st,omegav
    real(wp) :: eta_0,B_1,deta,ror,Tr,den
    ! ----------------------------------------------------------------------------

    ! Dilute-gas limit viscosity (given in P=Poise=0.1 Pa.s in Chung et al. 1988)
    ! ==========================
    ! dimensionless temperature T_st=k*T/eps
    !                           with k: Boltzmann cste & eps: potential energy param.    
    T_st=1.2593_wp*T/Tc

    ! reduced collision integral
    ! (empirical eq. from Neufeld et al.,1972)
    omegav=avisc/T_st**(bvisc) + cvisc/exp(dvisc*T_st) &
          +evisc/exp(fvisc*T_st) &
          +gvisc*T_st**(bvisc)*sin(svisc*T_st**(wvisc)-hvisc)

    eta_0=k_chung*sqrt(T)/omegav*0.1_wp ! convert in Pa.s

    ! First-density coefficient (initial density dependence of viscosity)
    ! =========================
    ! Vogel et al. (Journal of Physical and Chemical Reference Data 27, 947, 1998)
    B_1=bw(1)+bw(2)*T_st**(-0.25_wp)+bw(3)*T_st**(-0.5_wp)+bw(4)*T_st**(-0.75_wp) &
       +bw(5)*T_st**(-1.0_wp)+bw(6)*T_st**(-1.25_wp)+bw(7)*T_st**(-1.5_wp) &
       +bw(8)*T_st**(-2.5_wp)+bw(9)*T_st**(-5.5_wp)

    ! Residual term (empirical fit, given in micro Pa.s in Wen et al.)
    ! =============
    Tr=T/Tc ! reduced temperature
    ror=ro/roc0 ! reduced density
    den=cw(3)+cw(4)*ror+cw(5)*Tr+cw(6)*ror*Tr+cw(7)*ror**2*Tr
    deta=ror**(2.0_wp/3.0_wp)*sqrt(Tr)*(cw(1)+cw(2)/den)

    ! Composite solution
    ! ==================    
    eta=eta_0*(1.0_wp+B_1*ro)+deta*1.e-6_wp

    viscosity_law_wen= eta

  end function viscosity_law_wen

end submodule smod_tranprop_visc
