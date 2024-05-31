!=================================================================================
module mod_ineos_prs
!=================================================================================
  !> Module to share constants of PENG-ROBINSON EoS [for inlining purposes]
!=================================================================================
  use precision
  implicit none
  ! ------------------------------------------------------------------------------
  real(wp), parameter :: tol=1.e-6_wp,tol1=1.e-4_wp
  real(wp), parameter :: sqr8=sqrt(8.0_wp)
  real(wp) :: Kpr
  real(wp) :: apr,bpr
  ! ------------------------------------------------------------------------------
end module mod_ineos_prs

!=================================================================================
module mod_eos_prs
!=================================================================================
  !> Module to define subroutines for PENG-ROBINSON Equation of State (EoS)
!=================================================================================
  use mod_fluid
  use mod_ineos_prs
  use warnstop
  implicit none
  ! ------------------------------------------------------------------------------
  ! ------------------------------------------------------------------------------
  
contains
  
  !===============================================================================
  subroutine init_eos_prs
  !===============================================================================
    !> Initializations of gas coefficients
    !> - PRS EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp) :: L0,L1,L2,L3
    ! ----------------------------------------------------------------------------

    call read_eos
    
    ! The general parameters of the feos.ini file are corrected in order to be
    ! compatible with the PRSV equation:
    ! ----------------------------------
    zc = 0.3112_wp  ! Novec649
    !zc = 0.27795_wp ! MDM
    roc0=roc
    roc= pc/(zc*rg*Tc)
    apr= 0.457235_wp*(Tc*rg)**2/pc
    bpr= 0.077796_wp*rg*Tc/pc
    L0 = 0.378893_wp
    L1 = 1.4897153_wp
    L2 = 0.17131848_wp
    L3 = 0.0196554_wp
    
    ! gamma - 1
    ! ---------
    gam1= gam-1.0_wp
    !cvfg = 1.0_wp/gam1*rg
    !cpfg = gam/gam1*rg
    !s_crit = scalc_tro_prs(Tc,roc)
    
    Kpr=L0+L1*om-L2*om**2+L3*om**3
    
  end subroutine init_eos_prs

  !===============================================================================
  function alp_prs(T)
  !===============================================================================
    !> Compute Peng-Robinson function alpha
    !> - PRS EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T
    real(wp) :: alp_prs ! OUTPUT
    ! ----------------------------------------------------------------------------

    alp_prs=(1.0_wp+Kpr*(1.0_wp-sqrt(abs(T/Tc))))**2

  end function alp_prs

  !===============================================================================
  function avcalc_tro_prs(T,ro)
  !===============================================================================
    !> Compute isobaric expansion coefficient from T and rho
    !> - PRS EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: avcalc_tro_prs ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: v,dpdv,dpdT,dalpdT
    ! ----------------------------------------------------------------------------

    ! specific volume
    v= 1.0_wp/ro

    dpdv= -rg*T/(v-bpr)**2 &
        + 2.0_wp*apr*alp_prs(T)*(v+bpr)/(v**2+2.0_wp*bpr*v-bpr**2)**2

    dalpdT=2.0_wp*(1.0_wp+Kpr*(1.0_wp-sqrt(T/Tc)))*(-Kpr/(2.0_wp*(sqrt(T*Tc))))

    dpdT= rg/(v-bpr)-(dalpdT)*apr/(v**2+2.0_wp*bpr*v-bpr**2)

    avcalc_tro_prs= -ro*dpdT/dpdv

  end function avcalc_tro_prs

  !===============================================================================
  function c2calc_tro_prs(T,ro)
  !===============================================================================
    !> Compute speed of sound from T and rho
    !> - PRS EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: c2calc_tro_prs ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: v,dalpdT
    ! ----------------------------------------------------------------------------

    ! specific volume
    v= 1.0_wp/ro

    dalpdT=2.0_wp*(1.0_wp+Kpr*(1.0_wp-sqrt(T/Tc)))*(-Kpr/(2.0_wp*(sqrt(T*Tc))))

    c2calc_tro_prs= (-v**2)*(-rg*T/(v-bpr)**2 &
                  + 2.0_wp*apr*alp_prs(T)*(v+bpr)/(v**2+2.0_wp*bpr*v-bpr**2)**2 &
                  - T/cvcalc_id_tro_prs(T,ro)*(rg/(v-bpr)-dalpdT*apr/(v*v+2*v*bpr-bpr**2))**2)

  end function c2calc_tro_prs

  !===============================================================================
  function cpcalc_tro_prs(T,ro)
  !===============================================================================
    !> Compute heat capacity at constant pressure
    !> - PRS EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: cpcalc_tro_prs ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: v,dalpdT,d2alpdT2,cv,dpdv,dpdT,den
    ! ----------------------------------------------------------------------------

    ! specific volume
    v= 1.0_wp/ro
    
    dalpdT= 2.0_wp*(1.0_wp+Kpr*(1.0_wp-sqrt(T/Tc)))*(-Kpr/(2.0_wp*(sqrt(T*Tc))))

    d2alpdT2=(Kpr**2)/(2.0_wp*Tc*T)+Kpr/(4.0_wp*sqrt(Tc*T**3))*2.0_wp*(1+Kpr*(1.0_wp-sqrt(T/Tc)))

    cv= cvinf*abs(T/Tc)**nexp + apr*d2alpdT2*T/(bpr*sqr8) &
         *log(abs(2.0_wp*v+2.0_wp*bpr+bpr*sqr8)/abs(2.0_wp*v+2.0_wp*bpr-bpr*sqr8))

    den=1.0_wp/(v**2+2.0_wp*bpr*v-bpr**2)
    
    dpdv=-rg*T/(v-bpr)**2 + 2.0_wp*apr*alp_prs(T)*(v+bpr)*den**2

    dpdT= rg/(v-bpr)-(dalpdT)*apr*den

    cpcalc_tro_prs= cv - T*dpdT**2/dpdv

  end function cpcalc_tro_prs

  !===============================================================================
   function cvcalc_id_tro_prs(T,ro)
  !===============================================================================
    !> Compute heat capacity at constant volume in the ideal (dilute) limit
    !> - PRS EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: cvcalc_id_tro_prs ! OUTPUT
    ! ----------------------------------------------------------------------------

    cvcalc_id_tro_prs= cvinf*abs(T/Tc)**nexp

  end function cvcalc_id_tro_prs

  !===============================================================================
  function cvcalc_tro_prs(T,ro)
  !===============================================================================
    !> Compute heat capacity at constant volume
    !> - PRS EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: cvcalc_tro_prs ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: v,dalpdT,d2alpdT2
    ! ----------------------------------------------------------------------------

    ! cv = de/dT
    
    ! specific volume
    v = 1.0_wp/ro

    d2alpdT2= (Kpr**2)/(2.0_wp*Tc*T) + Kpr/(4.0_wp*sqrt(Tc*T**3))*2.0_wp*(1.0_wp+Kpr*(1.0_wp-sqrt(T/Tc)))

    cvcalc_tro_prs= cvinf*abs(T/Tc)**nexp + apr*d2alpdT2*T/(bpr*sqr8) &
         *log(abs(2.0_wp*v+2.0_wp*bpr+bpr*sqr8)/abs(2.0_wp*v+2.0_wp*bpr-bpr*sqr8))

  end function cvcalc_tro_prs
  
  !===============================================================================
  function dpdicalc_tro_prs(T,ro)
  !===============================================================================
    !> Compute pressure derivative w.r.t temperature
    !> - PRS EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: dpdicalc_tro_prs ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: v,dalpdT,d2alpdT2,cv,dpdT
    ! ----------------------------------------------------------------------------

    ! specific volume
    v = 1.0_wp/ro
    
    dalpdT = 2.0_wp*(1.0_wp+Kpr*(1.0_wp-sqrt(T/Tc)))*(-Kpr/(2.0_wp*(sqrt(T*Tc))))

    d2alpdT2 = (Kpr**2)/(2.0_wp*Tc*T) + Kpr/(4.0_wp*sqrt(Tc*T**3))*2.0_wp*(1.0_wp+Kpr*(1.0_wp-sqrt(T/Tc)))

    cv = cvinf*abs(T/Tc)**nexp + apr*d2alpdT2*T/(bpr*sqr8) &
         * log(abs(2.0_wp*v+2.0_wp*bpr+bpr*sqr8)/abs(2.0_wp*v+2.0_wp*bpr-bpr*sqr8))

    dpdT = rg/(v-bpr)-dalpdT*apr/(v**2+2.0_wp*bpr*v-bpr**2)

    dpdicalc_tro_prs= dpdT/cv

  end function dpdicalc_tro_prs

  !===============================================================================
  function dpdTcalc_tro_prs(T,ro)
  !===============================================================================
    !> Compute pressure derivative w.r.t temperature
    !> - PRS EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: dpdTcalc_tro_prs ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: v,dalpdT
    ! ----------------------------------------------------------------------------

    ! specific volume
    v = 1.0_wp/ro

    dalpdT = 2.0_wp*(1.0_wp+Kpr*(1.0_wp-sqrt(T/Tc)))*(-Kpr/(2.0_wp*(sqrt(T*Tc))))

    dpdTcalc_tro_prs = rg/(v-bpr)-dalpdT*apr/(v**2+2.0_wp*bpr*v-bpr**2)

  end function dpdTcalc_tro_prs

 !===============================================================================
  function dpdvcalc_tro_prs(T,ro)
  !===============================================================================
    !> Compute first-order pressure derivative w.r.t volume
    !> - PRS EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: dpdvcalc_tro_prs ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: v
    ! ----------------------------------------------------------------------------

    ! specific volume
    v= 1.0_wp/ro

    dpdvcalc_tro_prs=-rg*T/(v-bpr)**2 &
                    + 2.0_wp*apr*alp_prs(T)*(v+bpr)/(v**2+2.0_wp*bpr*v-bpr**2)**2

  end function dpdvcalc_tro_prs

  !===============================================================================
  function d2pdv2calc_tro_prs(T,ro)
  !===============================================================================
    !> Compute second-order pressure derivative w.r.t volume
    !> - PRS EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: d2pdv2calc_tro_prs ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: v
    ! ----------------------------------------------------------------------------

    ! specific volume
    v= 1.0_wp/ro

    d2pdv2calc_tro_prs=2.0_wp*rg*T/(v-bpr)**3 &
                      -2.0_wp*apr*alp_prs(T)*(3.0_wp*v**2+6.0_wp*bpr*v+5*bpr**2) &
                                             /(v**2+2.0_wp*bpr*v-bpr**2)**3

  end function d2pdv2calc_tro_prs

  !===============================================================================
  function d3pdv3calc_tro_prs(T,ro)
  !===============================================================================
    !> Compute third-order pressure derivative w.r.t volume
    !> - PRS EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: d3pdv3calc_tro_prs ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: v
    ! ----------------------------------------------------------------------------

    v= 1.0_wp/ro

    d3pdv3calc_tro_prs=-6.0_wp*rg*T/(v-bpr)**4 &
                      + 24.0_wp*apr*alp_prs(T)*(v+bpr)*(v**2+2.0_wp*bpr*v+3.0_wp*bpr**2) &
                                              /(v**2+2.0_wp*bpr*v-bpr**2)**4

  end function d3pdv3calc_tro_prs

  !===============================================================================
  function ecalc_pro_prs(p,ro,Ttent)
  !===============================================================================
    !> Compute internal energy from p and rho
    !> - PRS EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: p,ro,Ttent
    real(wp) :: ecalc_pro_prs ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: T,T1,fn,der_fn,v,err,dalpdT,den
    ! ----------------------------------------------------------------------------

    ! specific volume
    v= 1.0_wp / ro
    
    ! initial guess
    T= Ttent

    ! Newton's algorithm
    ! ------------------
    do i=1,100

       dalpdT=2.0_wp*(1.0_wp+Kpr*(1.0_wp-sqrt(T/Tc)))*(-Kpr/(2.0_wp*(sqrt(T*Tc))))

       ! function & derivative
       den=1.0_wp/(v**2+2.0_wp*bpr*v-bpr**2)
       fn= rg*T/(v-bpr) - apr*alp_prs(T)*den-p
       der_fn= rg/(v-bpr)-dalpdT*apr*den

       ! update solution
       T1= T - fn/der_fn

       err= abs(T1-T)/T
       if (err.le.tol) then
          ecalc_pro_prs= ecalc_tro_prs(T1,ro)
          return
       endif
       T= T1
    enddo

    print *,'in ecalc_pro: at iteration',i,'error=',err
    print *,'p :',p,'ro:',ro
    call mpistop('function ecalc_pro not converged',0)

  end function ecalc_pro_prs

  !===============================================================================
  function dedrocalc_tro_prs(T,ro)
  !===============================================================================
    !> Compute internal energy derivative w.r.t rho
    !> - PRS EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: dedrocalc_tro_prs ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: dalpdT
    ! ----------------------------------------------------------------------------

    dalpdT=2.0_wp*(1.0_wp+Kpr*(1.0_wp-sqrt(T/Tc)))*(-Kpr/(2.0_wp*(sqrt(T*Tc))))

    dedrocalc_tro_prs=-apr/bpr/sqr8*(alp_prs(T)-T*dalpdt)* &
                      (2.0_wp*bpr+bpr*sqr8)/((2.0_wp*bpr+bpr*sqr8)*ro+2)- &
                      (2.0_wp*bpr-bpr*sqr8)/((2.0_wp*bpr-bpr*sqr8)*ro+2)

  end function dedrocalc_tro_prs

  !===============================================================================
  function dedTcalc_tro_prs(T,ro)
  !===============================================================================
    !> Compute internal energy derivative w.r.t temperature
    !> - PRS EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: dedTcalc_tro_prs ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: d2alpdT2,v
    ! ----------------------------------------------------------------------------

    v=1.0_wp/ro

    ! second derivative of alpha
    d2alpdT2= (-Kpr/sqrt(T*Tc))*((Kpr/(-2.0_wp*sqrt(T*Tc)))-(1+Kpr*(1.0_wp-sqrt(T/Tc)))/(2.0_wp*T))

    dedTcalc_tro_prs=cvinf*(T/Tc)**nexp-apr/bpr/sqr8*T*d2alpdT2* &
                     log(abs(2.0_wp*v+2.0_wp*bpr+bpr*sqr8)/abs(2.0_wp*v+2.0_wp*bpr-bpr*sqr8))

  end function dedTcalc_tro_prs

  !===============================================================================
  function ecalc_tro_prs(T,ro)
  !===============================================================================
    !> Compute internal energy from T and rho
    !> - PRS EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: ecalc_tro_prs ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: v,dalpdT
    ! ----------------------------------------------------------------------------

    ! de = cvinf dT + (T*dp/dt-p) dv

    ! specific volume
    v= 1.0_wp/ro

    dalpdT=2.0_wp*(1.0_wp+Kpr*(1.0_wp-sqrt(T/Tc)))*(-Kpr/(2.0_wp*(sqrt(T*Tc))))

    ecalc_tro_prs= cvinf/(nexp+1.0_wp)*(abs(T/Tc))**(nexp+1.0_wp)*Tc    &
                 - apr*(alp_prs(T)-T*dalpdT)/(bpr*sqr8) &
                 *log(abs(2.0_wp*v+2.0_wp*bpr+bpr*sqr8)/abs(2.0_wp*v+2.0_wp*bpr-bpr*sqr8))

  end function ecalc_tro_prs

  !===============================================================================
  function gcalc_tro_prs(T,ro)
  !===============================================================================
    !> Compute fundamental derivative of gas dynamics from T and ro
    !> - PRS EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: gcalc_tro_prs ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: v,dpdv,d2pdv2,dpdT,d2pdT2,d2pdvdT
    real(wp) :: alpr,dalpdT,d2alpdT2,d3alpdT3,cv,dcvdT,c2calc,den,cc
    ! ----------------------------------------------------------------------------

    ! specific volume
    v= 1.0_wp/ro

    ! alpha
    alpr=alp_prs(T)
    
    ! first derivative of alpha
    dalpdT= -Kpr*(1+Kpr*(1.0_wp-sqrt(T/Tc)))/sqrt(T*Tc)

    ! seonc derivative of alpha
    d2alpdT2= (-Kpr/sqrt(T*Tc))*((Kpr/(-2.0_wp*sqrt(T*Tc)))-(1+Kpr*(1.0_wp-sqrt(T/Tc)))/(2.0_wp*T))

    ! third derivative of alpha
    d3alpdT3= -Kpr/sqrt(Tc*T)*(Kpr/(4.0_wp*sqrt(Tc*T**3)) - (-Kpr/(2.0_wp*sqrt(T*Tc)))/T &
            + 3.0_wp*(1+Kpr*(1.0_wp-sqrt(T/Tc)))/(4.0_wp*T*T))

    den=1.0_wp/(v**2+2.0_wp*bpr*v-bpr**2)

    dpdT= rg/(v-bpr)-dalpdT*apr*den

    d2pdT2= -apr*d2alpdT2*den

    d2pdvdT=-rg/(v-bpr)**2+2.0_wp*apr*dalpdT*(v+bpr)*den**2

    d2pdv2= 2.0_wp*rg*T/(v-bpr)**3 &
          - 2.0_wp*apr*alpr*(3.0_wp*v**2+6.0_wp*bpr*v+5.0_wp*bpr**2)*den**3

    cc=log(abs(2.0_wp*v+2.0_wp*bpr+bpr*sqr8)/abs(2.0_wp*v+2.0_wp*bpr-bpr*sqr8))
    
    cv = cvinf*abs(T/Tc)**nexp + apr*d2alpdT2*T/(bpr*sqr8)*cc

    dcvdT= nexp*cvinf*abs(T/Tc)**(nexp-1)/Tc + apr*(d2alpdT2+d3alpdT3*T)*1.0_wp/(bpr*sqr8)*cc

    dpdv=-rg*T/(v-bpr)**2 + 2.0_wp*apr*alpr*(v+bpr)*den**2

    c2calc = v**2*(T/cv*dpdT**2 - dpdv)

    gcalc_tro_prs= 0.5_wp*v**3/c2calc*( d2pdv2 - 3.0_wp*T/cv*dpdT*d2pdvdT &
                 + (T/cv*dpdT)**2*(3.0_wp*d2pdT2+dpdT/T*(1.0_wp-T/cv*dcvdT)) )

  end function gcalc_tro_prs

  !===============================================================================
  function pcalc_roero_prs(roe,ro,Ttent)
  !===============================================================================
    !> Compute pressure from rhoe and rho
    !> - PRS EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: roe,ro,Ttent
    real(wp) :: pcalc_roero_prs ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: v,T,T1,alpr,dalpdT,d2alpdT2
    real(wp) ::err,fn,der_fn,cc
    ! ----------------------------------------------------------------------------

    ! specific volume
    v= 1.0_wp/ro

    ! initial guess
    T= Ttent

    ! Newton's algorithm
    ! ------------------
    do i=1,100
       alpr= alp_prs(T)
       dalpdT= 2.0_wp*(1.0_wp+Kpr*(1.0_wp-sqrt(T/Tc)))*(-Kpr/(2.0_wp*(sqrt(T*Tc))))
       d2alpdT2= (Kpr**2)/(2.0_wp*Tc*T) + Kpr/(4.0_wp*sqrt(Tc*T**3))*2.0_wp*(1+Kpr*(1.0_wp-sqrt(T/Tc)))

       ! function & derivative
       cc=log(abs(2.0_wp*v+2.0_wp*bpr+bpr*sqr8)/abs(2.0_wp*v+2.0_wp*bpr-bpr*sqr8))
       
       fn =  cvinf/(nexp+1.0_wp)*(abs(T/Tc))**(nexp+1.0_wp)*Tc &
            - apr*(alpr - T*dalpdT)* 1.0_wp/(bpr*sqr8)*cc - roe/ro

       der_fn = cvinf*abs(T/Tc)**nexp + apr*d2alpdT2*T/(bpr*sqr8)*cc

       ! update solution
       T1= T - fn/der_fn

       err= abs(T1-T)/T
       if (err.le.tol) then
          pcalc_roero_prs= rg*T/(v-bpr)-apr*alpr/(v**2+2.0_wp*bpr*v-bpr**2)
          return
       endif
       T= T1

    enddo

    print *,'in pcalc_roero: at iteration',i,'error=',err
    print *,'roe:',roe,'ro:',ro
    pcalc_roero_prs= 0.0_wp
    call mpistop('function pcalc_roero not converged',0)
  
  end function pcalc_roero_prs

  !===============================================================================
  function pcalc_tro_prs(T,ro)
  !===============================================================================
    !> Compute pressure from T and rho
    !> - PRS EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: pcalc_tro_prs ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: v
    ! ----------------------------------------------------------------------------
    
    ! specific volume
    v= 1.0_wp/ro

    pcalc_tro_prs= rg*T/(v-bpr) - apr*alp_prs(T)/(v**2+2.0_wp*v*bpr-bpr**2)

  end function pcalc_tro_prs

  !===============================================================================
  function rocalc_ep_prs(e,p,Ttent)
  !===============================================================================
    !> Compute density from internal energy e and p
    !> - PRS EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: e,p,Ttent
    real(wp) :: rocalc_ep_prs ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: T,T1,ro,err
    ! ----------------------------------------------------------------------------

    ! initial guess
    T1= Ttent
    ro= 0.5_wp*roc
    
    ! iterations
    do i=1,100
       T = T1
       ro= rocalc_pt_prs(p,T,ro)
       T1= tcalc_roero_prs(ro*e,ro,T1)
       err= abs(T1-T)/T

       if (err.lt.tol) then
          rocalc_ep_prs= ro
          return
       endif
    enddo

    print *,'in rocalc_ep: at iteration',i,'error=',err
    call mpistop('function rocalc_ep not converged',0)

  end function rocalc_ep_prs

  !===============================================================================
  function rocalc_ps_prs(p,s,Ttent)
  !===============================================================================
    !> Compute density from p and entropy s
    !> - PRS EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: p,s,Ttent
    real(wp) :: rocalc_ps_prs ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: T,T1,ro,err
    ! ----------------------------------------------------------------------------

    ! initial guess
    T1= Ttent
    ro= 0.1_wp*roc
    
    ! iterations
    do i=1,100
       T = T1
       ro= rocalc_pt_prs(p,T,ro)
       T1= tcalc_sro_prs(s,ro,T1)
       err= abs(T1-T)/T
       if (err.lt.tol) then
          rocalc_ps_prs= ro
          return
       endif
    enddo

    print *,'in rocalc_ps: at iteration',i,'error=',err
    print *,'p:',p,'s:',s
    rocalc_ps_prs= 0.0_wp
    call mpistop('function rocalc_ps not converged',0)

  end function rocalc_ps_prs

  !===============================================================================
  function rocalc_pt_prs(p,T,rotent)
  !===============================================================================
    !> Compute density from p and T
    !> - PRS EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: p,T,rotent
    real(wp) :: rocalc_pt_prs ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: v,v1,fn,der_fn,err
    ! ----------------------------------------------------------------------------

    ! initial guess (specific volume)
    v= 1.0_wp/rotent
    
    ! Newton's algorithm
    ! ------------------
    do i=1,1000

       ! function & derivative
       fn = rg*T/(v-bpr)- apr*alp_prs(T)/(v**2+2.0_wp*bpr*v-bpr**2) - p
       der_fn = -rg*T/(v-bpr)**2 + 2.0_wp*apr*alp_prs(T)*(v+bpr)/(v**2+2.0_wp*v*bpr-bpr**2)**2

       ! update solution
       v1 = v - fn/der_fn

       err = abs(v1-v)/v

       if (err.le.tol) then
          rocalc_pt_prs= 1.0_wp/v1
          return
       endif
       v= v1
    enddo

    print *,'in rocalc_pt: at iteration',i,'error=',err
    print *,'p:',p/pc,'T:',T/tc
    rocalc_pt_prs= 0.0_wp
    call mpistop('function rocalc_pt not converged',0)

  end function rocalc_pt_prs

  !===============================================================================
  function rocalc_st_prs(s,T,rotent)
  !===============================================================================
    !> Compute density from entropy and T
    !> - PRS EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: s,T,rotent
    real(wp) :: rocalc_st_prs ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: v,v1,dalpdT
    real(wp) :: fn,der_fn,err
    ! ----------------------------------------------------------------------------

    dalpdT = 2.0_wp*(1.0_wp+Kpr*(1.0_wp-sqrt(T/Tc)))*(-Kpr/(2.0_wp*(sqrt(T*Tc))))

    ! initial guess (specific volume)
    v = 1.0_wp/rotent

    ! Newton's algorithm
    ! ------------------
    do i=1,1000

       ! function & derivative
       fn = cvinf/nexp*abs(T/Tc)**nexp + rg*log(v-bpr)+ apr*dalpdT* 1.0_wp/(bpr*sqr8)&
            *log(abs(2.0_wp*v+2.0_wp*bpr+bpr*sqr8)/abs(2.0_wp*v+2.0_wp*bpr-bpr*sqr8)) - s
       der_fn = rg/(v-bpr)-apr*dalpdT/(v**2+2.0_wp*bpr*v-bpr**2)

       ! update solution
       v1= v - fn/der_fn

       err = abs(v1-v)/v
       if (err.le.tol) then
          rocalc_st_prs= 1.0_wp/v1
          return
       endif
       v= v1

    enddo

    print *,'in rocalc_st: at iteration',i,'error=',err
    print *,'s:',s,'T:',T
    rocalc_st_prs= 0.0_wp
    call mpistop('function rocalc_st not converged',0)

  end function rocalc_st_prs

  !===============================================================================
  function scalc_tro_prs(T,ro)
  !===============================================================================
    !> Compute entropy s from T and rho
    !> - PRS EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: scalc_tro_prs ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: dalpdT,v
    ! ----------------------------------------------------------------------------

    ! ds = de/T + pdv/T

    ! specific volume
    v= 1.0_wp/ro

    dalpdT = 2.0_wp*(1.0_wp+Kpr*(1.0_wp-sqrt(T/Tc)))*(-Kpr/(2.0_wp*(sqrt(T*Tc))))

    scalc_tro_prs= cvinf/nexp*abs(T/Tc)**nexp + rg*log(v-bpr)+ apr*dalpdT* 1.0_wp/(bpr*sqr8)&
         * log(abs(2.0_wp*v+2.0_wp*bpr+bpr*sqr8)/abs(2.0_wp*v+2.0_wp*bpr-bpr*sqr8))

  end function scalc_tro_prs

  !===============================================================================
  function tcalc_ph_prs(p,h,rotent,Ttent)
  !===============================================================================
    !> Compute temperature from p and enthalpy h
    !> - PRS EOS - 
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: p,h,Ttent,rotent
    real(wp) :: tcalc_ph_prs ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: j
    real(wp) :: err,T,T1,ro,e
    ! ----------------------------------------------------------------------------

    ! initial guess
    T1= Ttent
    
    ! iterations
    do j=1,1000
       T  = T1
       ro = rocalc_pt_prs(p,T,rotent)

       e = h-p/ro
       T1= tcalc_roero_prs(ro*e,ro,T)
       ! e=  ecalc_pro(p,ro,Ttent)
       ! h1= e+p/ro
       ! err2= abs(h-h1)
       err= abs(T1-T)/T

       if (err.le.tol) then
          tcalc_ph_prs= T
          return
       endif
    enddo

    print *,'in tcalc_ph: at iteration',j,'error=',err
    print *,'T',T,'p:',p,'ro:',ro
    tcalc_ph_prs= 0.0_wp
    call mpistop('function tcalc_ph not converged',1)
    
  end function tcalc_ph_prs

  !===============================================================================
  function tcalc_pro_prs(p,ro,Ttent)
  !===============================================================================
    !> Compute temperature from p and rho
    !> - PRS EOS - 
    ! Tn+1 = Tn -(P(Tn)-Pref)/(dPdT(Tn))
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: p,ro,Ttent
    real(wp) :: tcalc_pro_prs ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: dalpdT,T,T1,v,err
    ! ----------------------------------------------------------------------------

    ! specific volume
    v= 1.0_wp / ro

    ! initial guess
    T = Ttent

    ! iterations
    do i=1,100

       dalpdT= -Kpr*(1.0_wp+Kpr*(1.0_wp-sqrt(T/Tc)))/sqrt(T*Tc)

       T1= T - (rg*T/(v-bpr)-apr*alp_prs(T)/(v**2+2.0_wp*bpr*v-bpr**2) - p) &
              /(rg/(v-bpr)-dalpdT*apr/(v**2 +2.0_wp*bpr*v-bpr**2))

       err= abs(T1-T)/T
       if (err.le.tol) then
          tcalc_pro_prs= T1
          return
       endif
       T= T1
    enddo

    print *,'in tcalc_pro: at iteration',i,'error=',err
    print *,'p:',p,'ro:',ro
    tcalc_pro_prs= 0.0_wp
    call mpistop('function tcalc_pro not converged',0)

  end function tcalc_pro_prs

  !===============================================================================
  function tcalc_roero_prs(roe,ro,Ttent)
  !===============================================================================
    !> Compute temperature from rho*e and rho
    !> - PRS EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: roe,ro,Ttent
    real(wp) :: tcalc_roero_prs ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: T,T1,err,v,fn,der_fn,dalpdT,d2alpdT2,cc
    ! ----------------------------------------------------------------------------

    ! specific volume
    v= 1.0_wp/ro
    cc=log(abs(2.0_wp*v+2.0_wp*bpr+bpr*sqr8)/abs(2.0_wp*v+2.0_wp*bpr-bpr*sqr8))

    ! initial guess
    T= Ttent

    ! Newton's algorithm
    ! ------------------
    do i=1,100

       dalpdT = 2.0_wp*(1.0_wp+Kpr*(1.0_wp-sqrt(T/Tc)))*(-Kpr/(2.0_wp*(sqrt(T*Tc))))
       d2alpdT2 = (Kpr**2)/(2.0_wp*Tc*T) + Kpr/(4.0_wp*sqrt(Tc*T**3))*2.0_wp*(1+Kpr*(1.0_wp-sqrt(T/Tc)))

       ! function & derivative       
       fn= cvinf/(nexp+1.0_wp)*(abs(T/Tc))**(nexp+1.0_wp)*Tc &
         - apr*(alp_prs(T)-T*dalpdT)/(bpr*sqr8)*cc  - roe/ro

       der_fn= cvinf*abs(T/Tc)**nexp + apr*d2alpdT2*T/(bpr*sqr8)*cc

       ! update solution
       T1= T - fn/der_fn

       err= abs(T1-T)/T
       if (err.le.tol) then
          tcalc_roero_prs= T1
          return
       endif
       T= T1
    enddo

    print *,'in tcalc_roero: at iteration',i,'error=',err
    print *,'roe:',roe,'ro:',ro,'T_t :',Ttent
    tcalc_roero_prs= 0.0_wp
    call mpistop('function tcalc_roero not converged',0)

  end function tcalc_roero_prs

  !===============================================================================
  function tcalc_sro_prs(s,ro,Ttent)
  !===============================================================================
    !> Compute temperature from entropy s and rho
    !> - PRS EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: s,ro,Ttent
    real(wp) :: tcalc_sro_prs ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: v,dalpdT,d2alpdT2,T,T1,fn,der_fn,err,cc
    ! ----------------------------------------------------------------------------

    ! specific volume
    v= 1.0_wp/ro

    T= Ttent
    
    ! Newton's algorithm
    ! ------------------
    do i=1,100

       dalpdT= 2.0_wp*(1.0_wp+Kpr*(1.0_wp-sqrt(T/Tc)))*(-Kpr/(2.0_wp*(sqrt(T*Tc))))

       d2alpdT2= (Kpr**2)/(2.0_wp*Tc*T) + Kpr/(4.0_wp*sqrt(Tc*T**3))*2.0_wp*(1+Kpr*(1.0_wp-sqrt(T/Tc)))

       ! function & derivative
       cc=log(abs(2.0_wp*v+2.0_wp*bpr+bpr*sqr8)/abs(2.0_wp*v+2.0_wp*bpr-bpr*sqr8))

       fn= cvinf/nexp*abs(T/Tc)**nexp+rg*log(v-bpr)+apr*dalpdT/(bpr*sqr8)*cc - s
       der_fn= cvinf*abs(T/Tc)**(nexp-1)/Tc + apr*d2alpdT2 * 1.0_wp/(bpr*sqr8)*cc

       ! update solution
       T1= T - fn/der_fn

       err= abs(T1-T)/T
       if (err.le.tol) then
          tcalc_sro_prs= T1
          return
       endif
       T = T1

    enddo

    print *,'in tcalc_sro: at iteration',i,'error=',err
    print *,'s :',s,'ro:',ro
    tcalc_sro_prs= 0.0_wp
    call mpistop('function tcalc_sro not converged',0)

  end function tcalc_sro_prs

  !===============================================================================
  function tcalc_Hs_prs(H,s,Ttent,rotent,ptent)
  !===============================================================================
    !> Compute temperature from total enthalpy and entropy s
    !> - PRS EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: H,s,Ttent,rotent,ptent
    real(wp) :: tcalc_Hs_prs ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: p,T,T1,ro,fn,der_fn,dpdT,err
    ! ----------------------------------------------------------------------------

    p = pTent
    T = Ttent
    ro= rotent

    ! Newton's algorithm
    ! ------------------
    do i=1,1000

       ! update density from total entropy
       ro=rocalc_st_prs(s,T1,ro)

       ! update pressure
       p=pcalc_tro_prs(T1,ro)

       ! update pressure derivative wrt T
       dpdT=dpdicalc_tro_prs(T1,ro)*cvcalc_tro_prs(T1,ro)

       ! update temperature
       fn = ecalc_tro_prs(T1,ro) + p/ro - H
       der_fn = dedTcalc_tro_prs(T1,ro) + dpdT/ro
       T= T1-fn/der_fn

       err= abs(T1-T)/T

       if (err.le.tol) then
          tcalc_Hs_prs = T
          return
       endif

       ! store previous solution
       T1=T

    enddo

    print *,'in tcalc_Hs_prs: at iteration',i,'error=',err
    print *,'s :',s,'H :',H,'T:',T
    tcalc_Hs_prs= 0.0_wp
    call mpistop('function tcalc_Hs_prs not converged',0)

  end function tcalc_Hs_prs

  !===============================================================================
  subroutine tcalc_Hstot_prs(T,ro,p,H_tot,s_tot,am0,coeff)
  !===============================================================================
    !> Compute temperature from total stagnation enthalpy (Newton iteration)
    !> used for imposing Riemann invariant in inflow BC
    !> - PRS EOS -
  !===============================================================================
    use warnstop
    implicit none
    ! ----------------------------------------------------------------------------
    ! input/output: initial and updated thermo. variables
    real(wp), intent(inout) :: T,ro,p
    ! input: total enthalpy and entropy to be conserved
     real(wp), intent(in) :: H_tot,s_tot
   ! input: interior part of Riemann invariant am0
    real(wp), intent(in) :: am0
    ! input: coeff taking into account flow direction
    real(wp), intent(in) :: coeff
    ! ----------------------------------------------------------------------------
    integer :: i,its
    real(wp) :: fn,der_fn,err,errs
    real(wp) :: T1,ro1,dpdT,e,dedT
    real(wp) :: v,v1,sqrt_Tr,log_v
    real(wp) :: alp,dalpdT,d2alpdT2
    ! ----------------------------------------------------------------------------

    T1=T
    ro1=ro

    ! **********************
    ! v1 : not inlined
    ! **********************

!!$    ! Newton's algorithm
!!$    ! ------------------
!!$    do i=1,1000
!!$
!!$       ! update density from total entropy
!!$       ro=rocalc_st_prs(s_tot,T1,ro1)
!!$
!!$       ! update pressure
!!$       p=pcalc_tro_prs(T1,ro)
!!$
!!$       ! update pressure derivative wrt T
!!$       dpdT=dpdicalc_tro_prs(T1,ro)*cvcalc_tro_prs(T1,ro)
!!$
!!$       ! update temperature
!!$
!!$       fn=ecalc_tro_prs(T1,ro)+p/ro+0.5_wp*coeff*(am0-p)**2-H_tot
!!$
!!$       der_fn=dedTcalc_tro_prs(T1,ro)+dpdT/ro-coeff*(am0-p)*dpdT
!!$
!!$       T= T1-fn/der_fn
!!$
!!$       err= abs(T1-T)/T
!!$
!!$       if (err.le.tol*10.) then
!!$          ro=rocalc_st_prs(s_tot,T,ro1)
!!$          p=pcalc_tro_prs(T,ro)
!!$          return
!!$       endif
!!$
!!$       ! store previous solution
!!$       T1=T
!!$       ro1=ro
!!$
!!$    enddo

    ! **********************
    ! v2 : partially inlined
    ! **********************

!!$    ! Newton's algorithm
!!$    ! ------------------
!!$    do i=1,50
!!$
!!$       ! compute intermediate quantities
!!$       sqrt_Tr=sqrt(T1/Tc)
!!$       alp=(1.0_wp+Kpr*(1.0_wp-sqrt_Tr))**2
!!$       dalpdT= 2.0_wp*(1.0_wp+Kpr*(1.0_wp-sqrt_Tr))*(-Kpr/(2.0_wp*(sqrt(T1*Tc))))
!!$       d2alpdT2= (Kpr**2)/(2.0_wp*Tc*T1) + Kpr/(4.0_wp*sqrt(Tc*T1**3))*2.0_wp*(1.0_wp+Kpr*(1.0_wp-sqrt_Tr))
!!$
!!$       ! update density from total entropy
!!$       ro=rocalc_st_prs(s_tot,T1,ro1)
!!$
!!$       ! specific volume
!!$       v= 1.0_wp/ro
!!$
!!$       ! update pressure [p=pcalc_tro_prs(T1,ro) inlined]
!!$       p=rg*T1/(v-bpr)-apr*alp/(v**2+2.0_wp*v*bpr-bpr**2)
!!$
!!$       ! update pressure derivative wrt T
!!$       ! [dpdT= dpdicalc_tro_prs(T1,ro)*cvcalc_tro_prs(T1,ro) inlined]
!!$       dpdT=rg/(v-bpr)-dalpdT*apr/(v**2+2.0_wp*bpr*v-bpr**2)
!!$
!!$       ! update internal energy (de=cvinf dT + (T*dp/dt-p) dv)
!!$       ! [ecalc_tro_prs(T1,ro) inlined]
!!$       log_v=log(abs(2.0_wp*v+2.0_wp*bpr+bpr*sqr8)/abs(2.0_wp*v+2.0_wp*bpr-bpr*sqr8))
!!$
!!$       e= cvinf/(nexp+1.0_wp)*(abs(T1/Tc))**(nexp+1.0_wp)*Tc &
!!$          - apr*(alp-T1*dalpdT)/(bpr*sqr8)*log_v
!!$
!!$       ! update internal energy derivative wrt T
!!$       ! [dedTcalc_tro_prs(T1,ro) inlined]
!!$       dedT=cvinf*(T1/Tc)**nexp-apr/bpr/sqr8*T1*d2alpdT2*log_v
!!$
!!$       ! update temperature
!!$
!!$       fn= e+p/ro+0.5_wp*coeff*(am0-p)**2-H_tot
!!$
!!$       der_fn= dedT+dpdT/ro-coeff*(am0-p)*dpdT
!!$
!!$       T= T1-fn/der_fn
!!$
!!$       err= abs(T1-T)/T
!!$
!!$       if (err.le.tol) then
!!$          ro=rocalc_st_prs(s_tot,T,ro1)
!!$          p=pcalc_tro_prs(T,ro)
!!$          return
!!$       endif
!!$
!!$       ! store previous solution
!!$       T1=T
!!$       ro1=ro
!!$
!!$    enddo

    ! **********************
    ! v3 : fully inlined
    ! **********************

    ! Newton's algorithm
    ! ------------------
    do i=1,50

       ! compute intermediate quantities
       sqrt_Tr=sqrt(T1/Tc)
       alp=(1.0_wp+Kpr*(1.0_wp-sqrt_Tr))**2
       dalpdT= 2.0_wp*(1.0_wp+Kpr*(1.0_wp-sqrt_Tr))*(-Kpr/(2.0_wp*(sqrt(T1*Tc))))
       d2alpdT2= (Kpr**2)/(2.0_wp*Tc*T1) + Kpr/(4.0_wp*sqrt(Tc*T1**3))*2.0_wp*(1.0_wp+Kpr*(1.0_wp-sqrt_Tr))

       ! update density from total entropy
       ! * initial guess (specific volume)
       v= 1.0_wp/ro1
       errs=1.0_wp
       its=0
       ! * Newton's algorithm
       do while ((errs>tol).and.(its<50))
          ! * function & derivative
          fn = cvinf/nexp*abs(T/Tc)**nexp + rg*log(v-bpr)+ apr*dalpdT* 1.0_wp/(bpr*sqr8)&
               *log(abs(2.0_wp*v+2.0_wp*bpr+bpr*sqr8)/abs(2.0_wp*v+2.0_wp*bpr-bpr*sqr8)) - s_tot
          der_fn = rg/(v-bpr)-apr*dalpdT/(v**2+2.0_wp*bpr*v-bpr**2)
          ! * update solution
          v1= v-fn/der_fn
          errs=abs(v1-v)/v
          its=its+1
          v= v1
       enddo

       ! update pressure
       !p=pcalc_tro_prs(T1,ro)
       p=rg*T1/(v-bpr)-apr*alp/(v**2+2.0_wp*v*bpr-bpr**2)

       ! update pressure derivative wrt T
       !dpdT= dpdicalc_tro_prs(T1,ro)*cvcalc_tro(T1,ro)
       dpdT=rg/(v-bpr)-dalpdT*apr/(v**2+2.0_wp*bpr*v-bpr**2)

       ! update internal energy (de=cvinf dT + (T*dp/dt-p) dv)
       log_v=log(abs(2.0_wp*v+2.0_wp*bpr+bpr*sqr8)/abs(2.0_wp*v+2.0_wp*bpr-bpr*sqr8))

       e= cvinf/(nexp+1.0_wp)*(abs(T1/Tc))**(nexp+1.0_wp)*Tc &
          - apr*(alp-T1*dalpdT)/(bpr*sqr8)*log_v

       ! update internal energy derivative wrt T
       dedT=cvinf*(T1/Tc)**nexp-apr/bpr/sqr8*T1*d2alpdT2*log_v

       ! update temperature

       fn= e+p*v+0.5_wp*coeff*(am0-p)**2-H_tot

       der_fn= dedT+dpdT*v-coeff*(am0-p)*dpdT

       !fn= ecalc_tro_prs(T1,ro)+p/ro+0.5_wp*coeff*(am0-p)**2-H_tot

       !der_fn= dedTcalc_tro_prs(T1,ro)+dpdT/ro-coeff*(am0-p)*dpdT

       T= T1-fn/der_fn

       err= abs(T1-T)/T

       if (err.le.tol) then
          ro=rocalc_st_prs(s_tot,T,ro1)
          p=pcalc_tro_prs(T,ro)
          return
       endif

       ! store previous solution
       T1=T
       ro1=ro

    enddo

    print *,'in tcalc_Hstot: at iteration',i,'error=',err
    print *,'T :',T,'ro:',ro
    call mpistop('function tcalc_Hstot not converged',0)

  end subroutine tcalc_Hstot_prs

  !===============================================================================
  function vvol_prs(pi,Ti,vtent)
  !===============================================================================
    !> Compute vvol such that p(vvol)=p
    !> - PRS EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: pi,Ti,vtent
    real(wp) :: vvol_prs ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: v0,v1,err,p,T
    ! ----------------------------------------------------------------------------

    T= Ti*Tc
    p= pi*pc
    
    ! initial guess (specific volume)
    v0 = vtent/roc
   
    ! iterations
    do i=1,15000
       v1 = v0 - (p-pcalc_tro_prs(T,1.0_wp/v0))/(-dpdvcalc_tro_prs(T,1.0_wp/v0))

       err = abs(v1-v0)/v0
       if (err.le.tol1) then
          vvol_prs= v1*roc
          return
       endif

       v0=v1
    enddo

    print *,'in vvol: at iteration',i,'error=',err
    call mpistop('function vvol not converged',0)
    
  end function vvol_prs

  !===============================================================================
  function vvol_d1_prs(Ti,vtent)
    implicit none
  !===============================================================================
    !> Compute vvol_d1 from Ti and tentative v
    !> - PRS EOS -
  !===============================================================================
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: Ti,vtent
    real(wp) :: vvol_d1_prs ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: v0,v1,err,T
    ! ----------------------------------------------------------------------------

    T= Ti*Tc
    
    ! initial guess (specific volume)
    v0 = vtent/roc

    ! iterations
    do i=1,15000
       v1=v0-dpdvcalc_tro_prs(T,1.0_wp/v0)/(d2pdv2calc_tro_prs(T,1.0_wp/v0)+1.e-16_wp)

       err = abs(v1-v0)/v0
       if (err.le.tol1) then
          vvol_d1_prs= v1*roc
          return
       endif

       v0=v1
    enddo

    print *,'in vvol_d1: at iteration',i,'error=',err
    call mpistop('function vvol_d1 not converged',0)

  end function vvol_d1_prs

  !===============================================================================
  function vvol_d2_prs(Ti,vtent)
  !===============================================================================
    !> Compute vvol_d2 from Ti and tentative v
    !> - PRS EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: Ti,vtent
    real(wp) :: vvol_d2_prs ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: v0,v1,err,T
    ! ----------------------------------------------------------------------------

    T = Ti*Tc
    ! initial guess (specific volume)
    v0 = vtent/roc

    ! iterations
    do i=1,15000
       v1= v0-d2pdv2calc_tro_prs(T,1.0_wp/v0)/(d3pdv3calc_tro_prs(T,1.0_wp/v0)+1.e-16_wp)

       err= abs(v1-v0)/v0
       if (err.le.tol1) then
          vvol_d2_prs= v1*roc
          return
       endif

       v0=v1
    enddo

    print *,'in vvol_d1: at iteration',i,'error=',err
    call mpistop('function vvol_d2 not converged',0)

  end function vvol_d2_prs

  !===============================================================================
  function intpcalc_tro_prs(Ti,roi)
  !===============================================================================
    !> Compute pressure integral from T and rho
    !> - PRS EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: Ti,roi
    real(wp) :: intpcalc_tro_prs ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: v,ro,T
    ! ----------------------------------------------------------------------------

    T = Ti*Tc
    ro = roi*roc

    ! specific volume
    v= 1.0_wp/ro

    intpcalc_tro_prs= rg*T*log(v-bpr) + apr*alp_prs(T)/(bpr*sqr8)&
         *log(abs(2.0_wp*v+2.0_wp*bpr+bpr*sqr8)/abs(2.0_wp*v+2.0_wp*bpr-bpr*sqr8))

    intpcalc_tro_prs= intpcalc_tro_prs/pc*roc

  end function intpcalc_tro_prs


  ! !===============================================================================
  ! subroutine stagnation_prs(h0,s0,e0,T0,p0,ro0)
  ! !===============================================================================
  !   !> Compute stagnation quantities from static quantities & h0 & s0
  !   !> - PRS EOS -
  ! !===============================================================================
  !   use mod_mpi
  !   implicit none
  !   ! ----------------------------------------------------------------------------
  !   real(wp), intent(in) :: h0,s0
  !   real(wp), intent(inout) :: e0,T0,p0,ro0
  !   ! ----------------------------------------------------------------------------
  !   integer :: j
  !   real(wp) :: err,p1
  !   ! ----------------------------------------------------------------------------

  !   ! ! Version 1
  !   ! p1 = p0
  !   ! err = 1.0_wp
  !   ! do j=1,1000
  !   !    p0 = p1
  !   !    T0 = tcalc_ph_prs(p0,h0,ro0,T0)
  !   !    ro0 = rocalc_st_prs(s0,T0,ro0)
  !   !    e0 = h0 - p0/ro0
  !   !    p1 = pcalc_roero_prs(ro0*e0,ro0,T0)
  !   !    err = abs(p1-p0)/p0

  !   !    if (err.le.tol) return
  !   ! enddo


  !   ! Version 2
  !   p1 = p0
  !   err = 1.0_wp
  !   do j=1,1000
  !      p0 = p1
  !      T0 = tcalc_sro_prs(s0,ro0,T0)
  !      ro0 = rocalc_st_prs(s0,T0,ro0)
  !      e0 = h0 - p0/ro0
  !      p1 = pcalc_roero_prs(ro0*e0,ro0,T0)
  !      err = abs(p1-p0)/p0

  !      ! if (err.le.tol) return
  !      if (err.le.10*tol) return
  !   enddo

  !   print *,'in stagnation_calc: at iteration',j,'error=',err
  !   print *,"h0,s0,e0,T0,p0,ro0",h0,s0,e0,T0,p0,ro0
  !   call mpistop('function stagnation_calc not converged',1)

  ! end subroutine stagnation_prs

  !===============================================================================
  subroutine stagnation_prs(h0,s,e0,T0,p0,ro0)
  !===============================================================================
    !> Compute stagnation quantities from h0 & s
    !> - PRS EOS -
  !===============================================================================
    use warnstop
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: h0,s
    real(wp), intent(inout) :: T0,ro0
    real(wp), intent(out) :: e0,p0
    real(wp) :: tcalc_Hs_prs ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: T1,fn,der_fn,dpdT,err
    ! ----------------------------------------------------------------------------

    T1 = T0

    ! Newton's algorithm
    ! ------------------
    do i=1,1000
       ! update density from total entropy
       ro0=rocalc_st_prs(s,T1,ro0)

       ! update pressure
       p0=pcalc_tro_prs(T1,ro0)

       ! update pressure derivative wrt T0
       dpdT=dpdicalc_tro_prs(T1,ro0)*cvcalc_tro_prs(T1,ro0)

       ! update temperature
       fn = ecalc_tro_prs(T1,ro0) + p0/ro0 - h0
       der_fn = dedTcalc_tro_prs(T1,ro0) + dpdT/ro0
       T0= T1-fn/der_fn

       err= abs(T1-T0)/T0

       if (err.le.tol) then
          ro0=rocalc_st_prs(s,T0,ro0)
          p0=pcalc_tro_prs(T0,ro0)
          e0 = h0 - p0/ro0
          return
       endif

       ! store previous solution
       T1=T0

    enddo

    print *,'in stagnation_prs: at iteration',i,'error=',err
    print *,'s :',s,'h0 :',h0,'T0:',T0,'p0:',ro0
    call mpistop('subroutine stagnation_prs not converged',0)

  end subroutine stagnation_prs

end module mod_eos_prs
