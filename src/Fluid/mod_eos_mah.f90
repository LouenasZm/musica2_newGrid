!=================================================================================
module mod_ineos_mah
!=================================================================================
  !> Module to share constants of MARTIN-HOU EoS [for inlining purposes]
!=================================================================================
  use precision
  implicit none
  ! ------------------------------------------------------------------------------
  real(wp), parameter :: kmh=5.475_wp,tol=1.e-6_wp
  real(wp) :: a2mh,a3mh,a4mh
  real(wp) :: b2mh,b3mh,b4mh,b5mh
  real(wp) :: c2mh,c3mh,c4mh,c5mh
  real(wp) :: bmah,s_crit
  ! ------------------------------------------------------------------------------
end module mod_ineos_mah

!=================================================================================
module mod_eos_mah
!=================================================================================
  !> Module to define subroutines for MARTIN-HOU Equation of State (EoS)
!=================================================================================
  use mod_fluid
  use mod_ineos_mah
  use warnstop
  implicit none
  ! ------------------------------------------------------------------------------
  ! ------------------------------------------------------------------------------

contains

  !===============================================================================
  subroutine init_eos_mah
  !===============================================================================
    !> Initializations of gas coefficients
    !> - MAH EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp) :: beta,vc,Tp,Tb,tebr,psib,ac
    real(wp) :: f2,f3,f4,f5
    ! ----------------------------------------------------------------------------

    call read_eos

    ! for compatibility with cubic EoS that change roc value
    roc0=roc

    ! Compute of cpfg for compatibility with mod_eos_pfg in mod_tranprop
    ! ------------------------------------------------------------------
    cpfg= gam/(gam-1.0_wp)*rg

    vc= 1.0_wp/roc
    ! Values of coefficients can be found in Martin-Hou (1955)
    ! --------------------------------------------------------
    beta= zc*(20.533_wp - 31.883_wp*zc)
    bmah= vc*(1.0_wp  - beta/(15.0_wp*zc))
    Tp  = Tc*(0.9869_wp - 0.6751_wp*zc)

    ! Tb: Boyle Temperature
    ! ---------------------
    Tb= 9.654291_wp + 2.505308_wp*Tc - 6.492324e-4_wp*Tc**2

    ! tebr: Reduced boiling point temperature
    tebr= teb/Tc

    ! Riedel relation [ac is the Riedel coefficient (e.g. Velasco 2008)]
    ! ---------------
    psib= -35.0_wp + 36.0_wp/tebr + 42.0_wp*log(tebr) - tebr**6
    ac  = (0.315_wp *psib - log(101325.0_wp/pc)) / &
         (0.0838_wp*psib - log(tebr))

    ! Computation of Martin-Hou coefficients
    ! --------------------------------------
    f2=  9.0_wp*pc   *(vc-bmah)**2 -  3.8_wp*rg*Tc*(vc-bmah)
    f3=  5.4_wp*rg*Tc*(vc-bmah)**2 - 17.0_wp*pc   *(vc-bmah)**3
    f4= 12.0_wp*pc   *(vc-bmah)**4 -  3.4_wp*rg*Tc*(vc-bmah)**3
    f5=  0.8_wp*rg*Tc*(vc-bmah)**4 -  3.0_wp*pc   *(vc-bmah)**5

    c2mh= ((f2 + bmah*rg*Tp + (rg*Tp)**2/pc*(1.0_wp-zc))*(Tb-Tc) + &
           (f2 + bmah*rg*Tb                            )*(Tc-Tp))/ &
             ((Tb-Tc)*(exp(-kmh) - exp(-kmh*Tp/Tc)) -              &
              (Tc-Tp)*(exp(-kmh*Tb/Tc) - exp(-kmh)))
    b2mh= (-f2 - bmah*rg*Tb - c2mh*(exp(-kmh*Tb/Tc) - exp(-kmh))) / (Tb-Tc)
    a2mh= f2 - b2mh*Tc - c2mh*exp(-kmh)

    ! ------------------------------------------------------------------------
    ! Modification were proposed in Martin-Hou to terms c3mh,b3mh and c5mh.
    ! Although,they have to be validated yet. Invert the commented lines to test.
    ! rn = 1.7
    ! c3mh = (c2mh*((vc-b)**3 - ((vc/rn-b)**3))) / ((vc/rn-b)**2 - (vc-b)**2)
    c3mh= -(vc-bmah)*c2mh
    ! c5mh = -c2mh*(vc-b)**3 - c3mh*(vc-b)**2
    c5mh= 0.0_wp
    ! b5mh = (f5 - c5mh*exp(-kmh))/Tc
    b5mh= f5/Tc
    ! -------------------------------------------------------------------------
    
    ! Transform ac to m = -ac*Pc/Tc
    ! -----------------------------
    b3mh= ac*pc/Tc*(vc-bmah)**3 - rg*(vc-bmah)**2 - b2mh*(vc-bmah) - b5mh/(vc-bmah)**2
    a3mh= f3 - b3mh*Tc - c3mh*exp(-kmh)
    a4mh= f4

    ! critical entropy 
    ! ----------------
    s_crit= scalc_tro_mah(Tc,1.0_wp/vc)

  end subroutine init_eos_mah

  !===============================================================================
  function avcalc_tro_mah(T,ro)
  !===============================================================================
    !> Compute isobaric expansion coefficient from T and rho
    !> - MAH EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: avcalc_tro_mah ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: dpdv,dpdT,v,vbm,ekTr
    ! ----------------------------------------------------------------------------

    ! specific volume
    v= 1.0_wp/ro
    
    vbm = v - bmah
    ekTr = exp(-kmh*T/Tc)

    dpdv=-                              rg*T/vbm**2 &
        - 2.0_wp*(a2mh + b2mh*T + c2mh*ekTr)/vbm**3 &
        - 3.0_wp*(a3mh + b3mh*T + c3mh*ekTr)/vbm**4 &
        - 4.0_wp*(a4mh                     )/vbm**5 &
        - 5.0_wp*(       b5mh*T + c5mh*ekTr)/vbm**6

    dpdT=                        rg/vbm    &
        + (b2mh - c2mh*kmh/Tc*ekTr)/vbm**2 &
        + (b3mh - c3mh*kmh/Tc*ekTr)/vbm**3 &
        + (b5mh - c5mh*kmh/Tc*ekTr)/vbm**5

    avcalc_tro_mah= -ro*dpdT/dpdv

  end function avcalc_tro_mah

  !===============================================================================
  function c2calc_tro_mah(T,ro)
  !===============================================================================
    !> Compute speed of sound from T and rho
    !> - MAH EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: c2calc_tro_mah ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: dpdv,dpdT,v,vbm,cv,cvres,ekTr
    ! ----------------------------------------------------------------------------

    ! specific volume
    v= 1.0_wp/ro
    
    vbm = v - bmah
    ekTr = exp(-kmh*T/Tc)
    cvres= c2mh/vbm + c3mh/(2.0_wp*vbm**2) + c5mh/(4.0_wp*vbm**4)

    cv= cvinf*abs(T/Tc)**nexp - T*kmh**2/Tc**2*ekTr*cvres

    dpdv=-                              rg*T/vbm**2 &
        - 2.0_wp*(a2mh + b2mh*T + c2mh*ekTr)/vbm**3 &
        - 3.0_wp*(a3mh + b3mh*T + c3mh*ekTr)/vbm**4 &
        - 4.0_wp*(a4mh                     )/vbm**5 &
        - 5.0_wp*(       b5mh*T + c5mh*ekTr)/vbm**6

    dpdT=                        rg/vbm    &
        + (b2mh - c2mh*kmh/Tc*ekTr)/vbm**2 &
        + (b3mh - c3mh*kmh/Tc*ekTr)/vbm**3 &
        + (b5mh - c5mh*kmh/Tc*ekTr)/vbm**5

    c2calc_tro_mah= v**2*(T/cv*dpdT**2 - dpdv)
    
  end function c2calc_tro_mah

  !===============================================================================
  function cpcalc_tro_mah(T,ro)
  !===============================================================================
    !> Compute heat capacity at constant pressure
    !> - MAH EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: cpcalc_tro_mah ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: v,vbm,ekTr,cvres,cv,dpdv,dpdT
    ! ----------------------------------------------------------------------------

    ! specific volume
    v= 1.0_wp/ro
    
    vbm = v - bmah
    ekTr = exp(-kmh*T/Tc)
    cvres= c2mh/vbm + c3mh/(2.0_wp*vbm**2) + c5mh/(4.0_wp*vbm**4)
    cv = cvinf*abs(T/Tc)**nexp - T*kmh**2/Tc**2*ekTr*cvres

    dpdv= -                             rg*T/vbm**2 &
        - 2.0_wp*(a2mh + b2mh*T + c2mh*ekTr)/vbm**3 &
        - 3.0_wp*(a3mh + b3mh*T + c3mh*ekTr)/vbm**4 &
        - 4.0_wp*(a4mh                     )/vbm**5 &
        - 5.0_wp*(       b5mh*T + c5mh*ekTr)/vbm**6

    dpdT=                        rg/vbm    &
        + (b2mh - c2mh*kmh/Tc*ekTr)/vbm**2 &
        + (b3mh - c3mh*kmh/Tc*ekTr)/vbm**3 &
        + (b5mh - c5mh*kmh/Tc*ekTr)/vbm**5

    cpcalc_tro_mah= cv - T*dpdT**2/dpdv
    
  end function cpcalc_tro_mah

  !===============================================================================
  function cvcalc_id_tro_mah(T,ro)
  !===============================================================================
    !> Compute heat capacity at constant volume in the ideal (dilute) limit
    !> - MAH EOS -
  !===============================================================================
     implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: cvcalc_id_tro_mah ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: v
    ! ----------------------------------------------------------------------------

    ! specific volume
    v= 1.0_wp/ro
    
    cvcalc_id_tro_mah= cvinf*abs(T/Tc)**nexp
    
  end function cvcalc_id_tro_mah

  !===============================================================================
  function cvcalc_tro_mah(T,ro)
  !===============================================================================
    !> Compute heat capacity at constant volume
    !> - MAH EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: cvcalc_tro_mah ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: v,ekTr,cvres
    ! ----------------------------------------------------------------------------

    ! specific volume
    v= 1.0_wp/ro
    
    ekTr = exp(-kmh*T/Tc)
    cvres= c2mh/(v-bmah) + c3mh/(2.0_wp*(v-bmah)**2) + c5mh/(4.0_wp*(v-bmah)**4)
    
    cvcalc_tro_mah= cvinf*abs(T/Tc)**nexp - T*kmh**2/Tc**2*ekTr*cvres

  end function cvcalc_tro_mah

  !===============================================================================
  function dpdicalc_tro_mah(T,ro)
  !===============================================================================
    !> Compute pressure derivative w.r.t temperature
    !> - MAH EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: dpdicalc_tro_mah ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: v,vbm,ekTr,cvres,cv,dpdT
    ! ----------------------------------------------------------------------------

    ! specific volume
    v= 1.0_wp/ro
    
    vbm = v - bmah
    ekTr = exp(-kmh*T/Tc)
    cvres= c2mh/vbm + c3mh/(2.0_wp*vbm**2) + c5mh/(4.0_wp*vbm**4)

    cv = cvinf*abs(T/Tc)**nexp - T*kmh**2/Tc**2*ekTr*cvres

    dpdT=                        rg/vbm    &
        + (b2mh - c2mh*kmh/Tc*ekTr)/vbm**2 &
        + (b3mh - c3mh*kmh/Tc*ekTr)/vbm**3 &
        + (b5mh - c5mh*kmh/Tc*ekTr)/vbm**5

    dpdicalc_tro_mah= dpdT/cv
    
  end function dpdicalc_tro_mah

  !===============================================================================
  function dpdTcalc_tro_mah(T,ro)
  !===============================================================================
    !> Compute pressure derivative w.r.t temperature
    !> - MAH EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: dpdTcalc_tro_mah ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: v,vbm,ekTr
    ! ----------------------------------------------------------------------------

    ! specific volume
    v= 1.0_wp/ro

    vbm = v - bmah
    ekTr = exp(-kmh*T/Tc)

    dpdTcalc_tro_mah=            rg/vbm    &
        + (b2mh - c2mh*kmh/Tc*ekTr)/vbm**2 &
        + (b3mh - c3mh*kmh/Tc*ekTr)/vbm**3 &
        + (b5mh - c5mh*kmh/Tc*ekTr)/vbm**5

   end function dpdTcalc_tro_mah

  !===============================================================================
  function dpdvcalc_tro_mah(T,ro)
  !===============================================================================
    !> Compute first-order pressure derivative w.r.t volume
    !> - MAH EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: dpdvcalc_tro_mah ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: v,vbm,ekTr,f1,f2,f3,f4,f5
    ! ----------------------------------------------------------------------------

    ! specific volume
    v= 1.0_wp/ro
    
    vbm = 1.0_wp/(v - bmah)
    ekTr= exp(-kmh*T/Tc)

    f1= rg*T
    f2= a2mh + b2mh*T + c2mh*ekTr
    f3= a3mh + b3mh*T + c3mh*ekTr
    f4= a4mh
    f5=        b5mh*T + c5mh*ekTr

    dpdvcalc_tro_mah= -      f1*vbm**2 &
                    - 2.0_wp*f2*vbm**3 &
                    - 3.0_wp*f3*vbm**4 &
                    - 4.0_wp*f4*vbm**5 &
                    - 5.0_wp*f5*vbm**6
    
  end function dpdvcalc_tro_mah

  !===============================================================================
  function d2pdv2calc_tro_mah(T,ro)
  !===============================================================================
    !> Compute second-order pressure derivative w.r.t volume
    !> - MAH EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: d2pdv2calc_tro_mah ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: v,vbm,ekTr,f1,f2,f3,f4,f5
    ! ----------------------------------------------------------------------------

    ! specific volume
    v= 1.0_wp/ro

    vbm = 1.0_wp/(v - bmah)
    ekTr= exp(-kmh*T/Tc)

    f1= rg*T
    f2= a2mh + b2mh*T + c2mh*ekTr
    f3= a3mh + b3mh*T + c3mh*ekTr
    f4= a4mh
    f5=       b5mh*T + c5mh*ekTr

    d2pdv2calc_tro_mah=2.0_wp*f1*vbm**3 &
                    +  6.0_wp*f2*vbm**4 &
                    + 12.0_wp*f3*vbm**5 &
                    + 20.0_wp*f4*vbm**6 &
                    + 30.0_wp*f5*vbm**7
    
  end function d2pdv2calc_tro_mah

  !===============================================================================
  function d3pdv3calc_tro_mah(T,ro)
  !===============================================================================
    !> Compute third-order pressure derivative w.r.t volume
    !> - MAH EOS -
  !===============================================================================
     implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: d3pdv3calc_tro_mah ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: v,vbm,ekTr,f1,f2,f3,f4,f5
    ! ----------------------------------------------------------------------------

    ! specific volume
    v= 1.0_wp/ro
    
    vbm = 1.0_wp/(v - bmah)
    ekTr= exp(-kmh*T/Tc)

    f1= rg*T
    f2= a2mh + b2mh*T + c2mh*ekTr
    f3= a3mh + b3mh*T + c3mh*ekTr
    f4= a4mh
    f5=        b5mh*T + c5mh*ekTr

    d3pdv3calc_tro_mah= - 6.0_wp*f1*vbm**4 &
                      -  24.0_wp*f2*vbm**5 &
                      -  60.0_wp*f3*vbm**6 &
                      - 120.0_wp*f4*vbm**7 &
                      - 210.0_wp*f5*vbm**8
    
  end function d3pdv3calc_tro_mah

  !===============================================================================
  function ecalc_pro_mah(p,ro,Ttent)
  !===============================================================================
    !> Compute internal energy from p and rho
    !> - MAH EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: p,ro,Ttent
    real(wp) :: ecalc_pro_mah ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: a,b,c,T,T1,ekTr,v,vbm,err
    ! ----------------------------------------------------------------------------

    ! specific volume
    v= 1.0_wp / ro
    vbm= v - bmah

    a= - p    + a2mh/vbm**2 + a3mh/vbm**3 + a4mh/vbm**4
    b= rg/vbm + b2mh/vbm**2 + b3mh/vbm**3 + b5mh/vbm**5
    c=          c2mh/vbm**2 + c3mh/vbm**3 + c5mh/vbm**5

    ! Newton's algorithm
    T= Ttent
    do i=1,100
       ekTr= exp(-kmh*T/Tc)
       T1= T - (a + b*T + c*ekTr)/(b - kmh/Tc*c*ekTr)

       err= abs(T1-T)/T
       if (err.le.tol) then
          ecalc_pro_mah= ecalc_tro_mah(T1,ro)
          return
       endif
       T= T1
    enddo

    print *,'in ecalc_pro: at iteration',i,'error=',err
    print *,'p:',p,'ro:',ro
    ecalc_pro_mah= 0.0_wp
    call mpistop('function ecalc_pro not converged',0)

  end function ecalc_pro_mah

  !===============================================================================
  function ecalc_tro_mah(T,ro)
  !===============================================================================
    !> Compute internal energy from T and rho
    !> - MAH EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: ecalc_tro_mah ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: v,vbm,Tr,ekTr! ,cvres
    ! ----------------------------------------------------------------------------

    ! specific volume
    v= 1.0_wp/ro
    vbm= v - bmah

    Tr=T/Tc
    ekTr = exp(-kmh*Tr)
    ! cvres= c2mh/vbm + c3mh/(2.0_wp*vbm**2) + c5mh/(4.0_wp*vbm**4)

    ! /!\ check why cvres lines are commented

    ecalc_tro_mah= cvinf/(nexp+1.0_wp)*(abs(Tr))**(nexp+1.0_wp)*Tc              &
                ! - cvres* (1.0_wp - ekTr*(1.0_wp+kmh*Tr))                       &
                ! + cvres*(exp(-kmh*Tr)*(1.0_wp+kmh*Tr)-exp(-kmh)*(1.0_wp+kmh))  &
                 + (c2mh*Tr*kmh*ekTr + a2mh + c2mh*ekTr)/      vbm              &
                 + (c3mh*Tr*kmh*ekTr + a3mh + c3mh*ekTr)/(2.0_wp*vbm**2)        &
                 + (                   a4mh            )/(3.0_wp*vbm**3)        &
                 + (c5mh*Tr*kmh*ekTr +        c5mh*ekTr)/(4.0_wp*vbm**4)
    
  end function ecalc_tro_mah

  !===============================================================================
  function dedTcalc_tro_mah(T,ro)
  !===============================================================================
    !> Compute internal energy from T and rho
    !> - MAH EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: dedTcalc_tro_mah ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: v,vbm,Tr,ekTr
    ! ----------------------------------------------------------------------------

    ! specific volume
    v= 1.0_wp/ro
    vbm= v - bmah

    Tr=T/Tc
    ekTr = exp(-kmh*Tr)

    ! [/!\ without cvres lines]
    dedTcalc_tro_mah= cvinf*(abs(Tr))**nexp                   &
                    - (c2mh*Tr*(kmh**2)*ekTr)/        vbm     &
                    - (c3mh*Tr*(kmh**2)*ekTr)/(2.0_wp*vbm**2) &
                    - (c5mh*Tr*(kmh**2)*ekTr)/(4.0_wp*vbm**4)

  end function dedTcalc_tro_mah

  !===============================================================================
  function gcalc_tro_mah(T,ro)
  !===============================================================================
    !> Compute fundamental derivative of gas dynamics from T and ro
    !> - MAH EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: gcalc_tro_mah ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: v,vbm,c2calc,cv,cvres,ekTr
    real(wp) :: dpdv,dpdT,d2pdT2,d2pdv2,d2pdvdT,dcvdT
    ! ----------------------------------------------------------------------------

    ! specific volume
    v= 1.0_wp/ro    
    vbm= v - bmah
    
    ekTr = exp(-kmh*T/Tc)
    cvres= c2mh/vbm + c3mh/(2.0_wp*vbm**2) + c5mh/(4.0_wp*vbm**4)

    cv= cvinf*abs(T/Tc)**nexp - T*kmh**2/Tc**2*ekTr*cvres

    dpdv= -                             rg*T/vbm**2 &
        - 2.0_wp*(a2mh + b2mh*T + c2mh*ekTr)/vbm**3 &
        - 3.0_wp*(a3mh + b3mh*T + c3mh*ekTr)/vbm**4 &
        - 4.0_wp*(a4mh                     )/vbm**5 &
        - 5.0_wp*(       b5mh*T + c5mh*ekTr)/vbm**6

    !if (dpdv.gt.0_wp) dpdv = 0.0_wp

    dpdT=                        rg/vbm    &
        + (b2mh - c2mh*kmh/Tc*ekTr)/vbm**2 &
        + (b3mh - c3mh*kmh/Tc*ekTr)/vbm**3 &
        + (b5mh - c5mh*kmh/Tc*ekTr)/vbm**5

    d2pdT2= kmh**2/Tc**2*ekTr * (c2mh/vbm**2 + c3mh/vbm**3 + c5mh/vbm**5)

    d2pdv2= 2.0_wp*                       rg*T/vbm**3 &
          + 6.0_wp*(a2mh + b2mh*T + c2mh*ekTr)/vbm**4 &
          +12.0_wp*(a3mh + b3mh*T + c3mh*ekTr)/vbm**5 &
          +20.0_wp*(a4mh                     )/vbm**6 &
          +30.0_wp*(       b5mh*T + c5mh*ekTr)/vbm**7

    d2pdvdT= -                             rg/vbm**2 &
           - 2.0_wp*(b2mh - c2mh*kmh/Tc*ekTr)/vbm**3 &
           - 3.0_wp*(b3mh - c3mh*kmh/Tc*ekTr)/vbm**4 &
           - 5.0_wp*(b5mh - c5mh*kmh/Tc*ekTr)/vbm**6

    ! dcvdT = 0.0_wp   !! This is in case of MH EoS for a calorically perfect gas

    dcvdT= nexp*cvinf*(abs(T/Tc))**(nexp-1.0_wp)/Tc + &
         kmh**2/Tc**2*ekTr*cvres*(kmh*T/Tc-1.0_wp)

    c2calc = v**2*(T/cv*dpdT**2 - dpdv)

    gcalc_tro_mah= 0.5_wp*v**3/c2calc*(d2pdv2 - 3.0_wp*T/cv*dpdT*d2pdvdT + &
                  (T/cv*dpdT)**2 * (3.0_wp*d2pdT2 + dpdT/T*(1.0_wp-T/cv*dcvdT)))

  end function gcalc_tro_mah

  !===============================================================================
  function pcalc_roero_mah(roe,ro,Ttent)
  !===============================================================================
    !> Compute pressure from rhoe and rho
    !> - MAH EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: roe,ro,Ttent
    real(wp) :: pcalc_roero_mah ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: a,cc1,cc2,T,T1,err,v,vbm,ekTr,fn,der_fn
    ! ----------------------------------------------------------------------------

    ! specific volume
    v= 1.0_wp/ro
    vbm= v - bmah

    cc1= (c2mh/vbm + c3mh/(2.0_wp*vbm**2) + c5mh/(4.0_wp*vbm**4))*kmh/Tc
    cc2= (c2mh/vbm + c3mh/(2.0_wp*vbm**2) + c5mh/(4.0_wp*vbm**4))
    a  = (a2mh/vbm + a3mh/(2.0_wp*vbm**2) + a4mh/(3.0_wp*vbm**3)) - roe/ro

    ! Newton's algorithm
    ! ------------------
    ! initial guess
    T= Ttent

    do i=1,100
       ekTr = exp(-kmh*T/Tc)
       ! function & derivative
       fn= a + cc1*T*ekTr + cc2*ekTr + cvinf*Tc/(nexp+1.0_wp)*(abs(T/Tc))**(nexp+1.0_wp)
       der_fn= (cc1 - kmh/Tc*cc2)*ekTr - kmh/Tc*cc1*T*ekTr + cvinf*(abs(T/Tc))**nexp

       !update solution
       T1= T - fn/der_fn

       err = abs(T1-T)/T
       if (err.le.tol) then
          ekTr= exp(-kmh*T1/Tc)
          pcalc_roero_mah=                        rg*T/vbm    &
                         + (a2mh + b2mh*T + c2mh*ekTr)/vbm**2 &
                         + (a3mh + b3mh*T + c3mh*ekTr)/vbm**3 &
                         + (a4mh                     )/vbm**4 &
                         + (       b5mh*T + c5mh*ekTr)/vbm**5
          return
       endif
       T= T1
    enddo

    print *,'in pcalc_roero: at iteration',i,'error=',err
    print *,'roe:',roe,'ro:',ro
    pcalc_roero_mah= 0.0_wp
    call mpistop('function pcalc_roero not converged',0)

  end function pcalc_roero_mah

  !===============================================================================
  function pcalc_tro_mah(T,ro)
  !===============================================================================
    !> Compute pressure from T and rho
    !> - MAH EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: pcalc_tro_mah ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: v,vbm,ekTr
    ! ----------------------------------------------------------------------------

    ! specific volume
    v= 1.0_wp/ro
    vbm = v - bmah
    ekTr= exp(-kmh*T/Tc)

    pcalc_tro_mah=                        rg*T/vbm    &
                 + (a2mh + b2mh*T + c2mh*ekTr)/vbm**2 &
                 + (a3mh + b3mh*T + c3mh*ekTr)/vbm**3 &
                 + (a4mh                     )/vbm**4 &
                 + (       b5mh*T + c5mh*ekTr)/vbm**5

  end function pcalc_tro_mah

  !===============================================================================
  function rocalc_ep_mah(e,p,Ttent)
  !===============================================================================
    !> Compute density from internal energy e and p
    !> - MAH EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: e,p,Ttent
    real(wp) :: rocalc_ep_mah ! OUTPUT
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
       ro= rocalc_pt_mah(p,T,ro)
       T1= tcalc_roero_mah(ro*e,ro,T1)
       err= abs(T1-T)/T

       if (err.lt.tol) then
          rocalc_ep_mah= ro
          return
       endif
    enddo

    print *,'in rocalc_ep: at iteration',i,'error=',err
    call mpistop('function rocalc_ep not converged',0)
   
  end function rocalc_ep_mah

  !===============================================================================
  function rocalc_ps_mah(p,s,Ttent)
  !===============================================================================
    !> Compute density from p and entropy s
    !> - MAH EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: p,s,Ttent
    real(wp) :: rocalc_ps_mah ! OUTPUT
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
       ro= rocalc_pt_mah(p,T,ro)
       T1= tcalc_sro_mah(s,ro,T1)
       err= abs(T1-T)/T
       if (err.lt.tol) then
          rocalc_ps_mah= ro
          return
       endif
    enddo
   
    print *,'in rocalc_ps: at iteration',i,'error=',err
    print *,'p:',p,'s:',s
    rocalc_ps_mah= 0.0_wp
    call mpistop('function rocalc_ps not converged',0)

  end function rocalc_ps_mah

  !===============================================================================
  function rocalc_pt_mah(p,T,rotent)
  !===============================================================================
    !> Compute density from p and T
    !> - MAH EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: p,T,rotent
    real(wp) :: rocalc_pt_mah ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: fn,der_fn,err
    real(wp) :: ekTr,f2,f3,f4,f5,v,v1,inv1,inv2,inv3,inv4,inv5,inv6
    ! ----------------------------------------------------------------------------

    ekTr= exp(-kmh*T/Tc)

    f2= a2mh + b2mh*T + c2mh*ekTr
    f3= a3mh + b3mh*T + c3mh*ekTr
    f4= a4mh
    f5=        b5mh*T + c5mh*ekTr

    ! initial guess (specific volume)
    v= 1.0_wp/rotent

    ! Newton's algorithm
    ! ------------------
    do i=1,1000

       inv1= 1.0_wp/(v-bmah)
       inv2= inv1*inv1
       inv3= inv2*inv1
       inv4= inv3*inv1
       inv5= inv4*inv1
       inv6= inv5*inv1

       ! function & derivative
       fn= rg*T*inv1 + f2*inv2 + f3*inv3 + f4*inv4 + f5*inv5 - p
       der_fn=    - rg*T*inv2 &
             - 2.0_wp*f2*inv3 &
             - 3.0_wp*f3*inv4 &
             - 4.0_wp*f4*inv5 &
             - 5.0_wp*f5*inv6

       ! update solution
       v1 = v - fn/der_fn

       err = abs(v1-v)/v

       if (err.le.tol) then
          rocalc_pt_mah= 1.0_wp/v1
          return
       endif
       v = v1
    enddo

    print *,'in rocalc_pt: at iteration',i,'error=',err
    print *,'p:',p/pc,'T:',T/tc
    rocalc_pt_mah= 0.0_wp
    call mpistop('function rocalc_pt not converged',0)

  end function rocalc_pt_mah

  !===============================================================================
  function rocalc_st_mah(s,T,rotent)
  !===============================================================================
    !> Compute density from entropy s, T
    !> - MAH EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: s,T,rotent
    real(wp) :: rocalc_st_mah ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: fn,der_fn,err
    real(wp) :: v,v1,ekTr,f2,f3,f5,fc
    ! ----------------------------------------------------------------------------

    ekTr = exp(-kmh*T/Tc)
    f2 = -b2mh + c2mh*kmh/Tc*ekTr
    f3 = -b3mh + c3mh*kmh/Tc*ekTr
    f5 = -b5mh + c5mh*kmh/Tc*ekTr
    fc = cvinf/nexp*(abs(T/Tc)**nexp - 1.0_wp)

    ! specific volume (initial guess)
    v = 1.0_wp/rotent

    ! Newton's algorithm
    ! ------------------    
    do i=1,100
       ! function & derivative
       fn= fc + rg*log((v-bmah)/(1.0_wp/roc-bmah)) + f2/(v-bmah) &
         + f3/(2.0_wp*(v-bmah)**2) + f5/(4.0_wp*(v-bmah)**4) - s_crit - s

       der_fn= fc + rg*1.0_wp/(v-bmah) &
             - f2/(v-bmah)**2  &
             - f3/(v-bmah)**3  &
             - f5/(v-bmah)**5

       ! update solution
       v1 = v - fn/der_fn

       err = abs(v1-v)/v
       if (err.le.tol) then
          rocalc_st_mah= 1.0_wp/v1
          return
       endif
       v= v1
    enddo

    print *,'in rocalc_st: at iteration',i,'error=',err
    print *,'s:',s,'T:',T
    rocalc_st_mah= 0.0_wp
    call mpistop('function rocalc_st not converged',0)
   
  end function rocalc_st_mah

  !===============================================================================
  function scalc_tro_mah(T,ro)
  !===============================================================================
    !> Compute entropy s from T and rho
    !> - MAH EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: scalc_tro_mah ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: cv,ekTr,cvres,v
    ! ----------------------------------------------------------------------------

    ! specific volume
    v= 1.0_wp/ro
    
    ekTr = exp(-kmh*T/Tc)
    cvres= c2mh/(v-bmah) + c3mh/(2.0_wp*(v-bmah)**2) + c5mh/(4.0_wp*(v-bmah)**4)
    cv = cvinf*abs(T/Tc)**nexp - T*kmh**2/Tc**2*ekTr*cvres

    scalc_tro_mah= cvinf/nexp*(abs(T/Tc)**nexp - 1.0_wp)          &
                !- k/Tc*(exp(-kmh)-exp(-k*T/Tc))       & 
                + rg*log((v-bmah)/(1.0_wp/roc-bmah))              &
                + (-b2mh + c2mh*kmh/Tc*ekTr)/(v-bmah)             &
                + (-b3mh + c3mh*kmh/Tc*ekTr)/(2.0_wp*(v-bmah)**2) &
                + (-b5mh + c5mh*kmh/Tc*ekTr)/(4.0_wp*(v-bmah)**4) - s_crit

  end function scalc_tro_mah

  !===============================================================================
  function tcalc_ph_mah(p,h,rotent,Ttent)
  !===============================================================================
    !> Compute temperature from p and enthalpy h
    !> - MAH EOS - 
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: p,h,Ttent,rotent
    real(wp) :: tcalc_ph_mah ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: j
    real(wp) :: err,T,T1,ro,e
    ! ----------------------------------------------------------------------------

    ! initial guess
    T1= Ttent
    
    ! Newton's algorithm
    ! ------------------
    do j=1,100
       T = T1
       ro= rocalc_pt_mah(p,T,rotent)

       e= h-p/ro
       T1= tcalc_roero_mah(ro*e,ro,T)
       ! e= ecalc_pro_mah(p,ro,Ttent)
       ! h1= e+p/ro
       ! err2= abs(h-h1)
       err= abs(T1-T)/T

       if (err.le.tol) then
          tcalc_ph_mah= T
          return
       endif
    enddo
    
  end function tcalc_ph_mah

  !===============================================================================
  function tcalc_pro_mah(p,ro,Ttent)
  !===============================================================================
    !> Compute temperature from p and rho
    !> - MAH EOS - 
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: p,ro,Ttent
    real(wp) :: tcalc_pro_mah ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: a,b,c,T,T1,ekTr,v,err
    ! ----------------------------------------------------------------------------

    ! specific volume
    v= 1.0_wp / ro

    a= - p         + a2mh/(v-bmah)**2 + a3mh/(v-bmah)**3 + a4mh/(v-bmah)**4
    b= rg/(v-bmah) + b2mh/(v-bmah)**2 + b3mh/(v-bmah)**3 + b5mh/(v-bmah)**5
    c=               c2mh/(v-bmah)**2 + c3mh/(v-bmah)**3 + c5mh/(v-bmah)**5

    ! initial guess
    T = Ttent

    ! Newton's algorithm
    ! ------------------
    do i=1,100
       ekTr= exp(-kmh*T/Tc)
       T1 = T-(a+b*T+c*ekTr)/(b-kmh/Tc*c*ekTr)

       err = abs(T1-T)/T
       if (err.le.tol) then
          tcalc_pro_mah= T1
          return
       endif
       T= T1
    enddo

    print *,'in tcalc_pro: at iteration',i,'error=',err
    print *,'p:',p,'ro:',ro
    tcalc_pro_mah= 0.0_wp
    call mpistop('function tcalc_pro not converged',0)

  end function tcalc_pro_mah

  !===============================================================================
  function tcalc_roero_mah(roe,ro,Ttent)
  !===============================================================================
    !> Compute temperature from rho*e and rho
    !> - MAH EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: roe,ro,Ttent
    real(wp) :: tcalc_roero_mah ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: a,cc1,cc2,T,T1,err,v,ekTr,fn,der_fn
    ! ----------------------------------------------------------------------------

    ! specific volume
    v= 1.0_wp/ro

    cc1= (c2mh/(v-bmah) + c3mh/(2.0_wp*(v-bmah)**2) + c5mh/(4.0_wp*(v-bmah)**4))*kmh/Tc
    cc2= (c2mh/(v-bmah) + c3mh/(2.0_wp*(v-bmah)**2) + c5mh/(4.0_wp*(v-bmah)**4))
    a  = (a2mh/(v-bmah) + a3mh/(2.0_wp*(v-bmah)**2) + a4mh/(3.0_wp*(v-bmah)**3)) - roe/ro

    ! initial guess
    T= Ttent

    ! Newton's algorithm
    ! ------------------
    do i=1,100
       ekTr= exp(-kmh*T/Tc)
       ! function & derivative
       fn = a + cc1*T*ekTr + cc2*ekTr + cvinf*Tc/(nexp+1.0_wp)*(abs(T/Tc))**(nexp+1.0_wp)
       der_fn = (cc1 - kmh/Tc*cc2)*ekTr - kmh/Tc*cc1*T*ekTr + cvinf*(abs(T/Tc))**nexp

       ! update solution
       T1= T - fn/der_fn

       err = abs(T1-T)/T
       if (err.le.tol) then
          tcalc_roero_mah= T1
          return
       endif
       T= T1
    enddo

    print *,'in tcalc_roero: at iteration',i,'error=',err
    print *,'roe:',roe,'ro:',ro,'T_t :',Ttent
    tcalc_roero_mah= 0.0_wp
    call mpistop('function tcalc_roero not converged',0)

  end function tcalc_roero_mah

  !===============================================================================
  function tcalc_sro_mah(s,ro,Ttent)
  !===============================================================================
    !> Compute temperature from entropy s and rho
    !> - MAH EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: s,ro,Ttent
    real(wp) :: tcalc_sro_mah ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: v,ekTr,T,T1,b,c,fn,der_fn,err
    ! ----------------------------------------------------------------------------

    ! specific volume
    v= 1.0_wp/ro

    b= rg*log((v-bmah)/(1.0_wp/roc-bmah))          &
     -  (b2mh/(v-bmah) + b3mh/(2.0_wp*(v-bmah)**2) + b5mh/(4.0_wp*(v-bmah)**4))
    c= kmh/Tc*(c2mh/(v-bmah) + c3mh/(2.0_wp*(v-bmah)**2) + c5mh/(4.0_wp*(v-bmah)**4))

    ! initial guess
    T= Ttent
    
    ! Newton's algorithm
    ! ------------------
    do i=1,100
       ekTr= exp(-kmh*T/Tc)

       ! function & derivative
       fn= cvinf/nexp*(abs(T/Tc)**nexp - 1.0_wp) + b + c*ekTr - s_crit - s
       der_fn= cvinf/Tc*abs(T/Tc)**(nexp-1.0_wp) - kmh/Tc*ekTr*c

       ! update solution
       T1= T - fn/der_fn

       err = abs(T1-T)/T
       if (err.le.tol) then
          tcalc_sro_mah= T1
          return
       endif
       T= T1

    enddo

    print *,'in tcalc_sro: at iteration',i,'error=',err
    print *,'s :',s,'ro:',ro
    tcalc_sro_mah= 0.0_wp
    call mpistop('function tcalc_sro not converged',0)

  end function tcalc_sro_mah

  !===============================================================================
  subroutine tcalc_Hstot_mah(T,ro,p,H_tot,s_tot,am0,coeff)
  !===============================================================================
    !> Compute temperature from total stagnation enthalpy (Newton iteration)
    !> used for imposing Riemann invariant in inflow BC
    !> - MAH EOS -
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
    integer :: i
    real(wp) :: fn,der_fn,err
    real(wp) :: T1,ro1,dpdT
    ! ----------------------------------------------------------------------------

    T1=T
    ro1=ro

    ! **********************
    ! v1 : not inlined
    ! **********************
    ! inlined version not written [cf mod_eos_prs.f90 for ex inlined version]

    ! Newton's algorithm
    ! ------------------
    do i=1,1000

       ! update density from total entropy
       ro=rocalc_st_mah(s_tot,T1,ro1)

       ! update pressure
       p=pcalc_tro_mah(T1,ro)

       ! update pressure derivative wrt T
       dpdT=dpdicalc_tro_mah(T1,ro)*cvcalc_tro_mah(T1,ro)

       ! update temperature

       fn=ecalc_tro_mah(T1,ro)+p/ro+0.5_wp*coeff*(am0-p)**2-H_tot

       der_fn=dedTcalc_tro_mah(T1,ro)+dpdT/ro-coeff*(am0-p)*dpdT

       T= T1-fn/der_fn

       err= abs(T1-T)/T

       if (err.le.tol*10.) then
          ro=rocalc_st_mah(s_tot,T,ro1)
          p=pcalc_tro_mah(T,ro)
          return
       endif

       ! store previous solution
       T1=T
       ro1=ro

    enddo

    print *,'in tcalc_Hstot: at iteration',i,'error=',err
    print *,'T :',T,'ro:',ro
    call mpistop('function tcalc_Hstot not converged',0)

  end subroutine tcalc_Hstot_mah

  !===============================================================================
  function vvol_mah(pi,Ti,vtent)
  !===============================================================================
    !> Compute vvol such that p(vvol)=0
    !> - MAH EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: pi,Ti,vtent
    real(wp) :: vvol_mah ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: v0,v1,err,p,T
    ! ----------------------------------------------------------------------------

    T = Ti*Tc
    p = pi*pc
    
    ! initial guess (specific volume)
    v0 = vtent/roc

    ! iterations
    do i=1,15000
       v1= v0 - (p-pcalc_tro_mah(T,1.0_wp/v0))/(-dpdvcalc_tro_mah(T,1.0_wp/v0))

       err = abs(v1-v0)/v0
       if (err.le.tol) then
          vvol_mah= v1*roc
          return
       endif

       v0=v1
    enddo

    print *,'in vvol: at iteration',i,'error=',err
    call mpistop('function vvol not converged',0)
    
  end function vvol_mah

  !===============================================================================
  function vvol_d1_mah(Ti,vtent)
  !===============================================================================
    !> Compute vvol_d1 from Ti and tentative v
    !> - MAH EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: Ti,vtent
    real(wp) :: vvol_d1_mah ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: v0,v1,err,T
    ! ----------------------------------------------------------------------------

    T = Ti*Tc
     
    ! initial guess (specific volume)
   v0 = vtent/roc

    ! iterations
    do i=1,15000
       v1= v0 - dpdvcalc_tro_mah(T,1.0_wp/v0)/d2pdv2calc_tro_mah(T,1.0_wp/v0)

       err= abs(v1-v0)/v0
       if (err.le.tol) then
          vvol_d1_mah= v1*roc
          return
       endif

       v0=v1
    enddo

    print *,'in vvol_d1: at iteration',i,'error=',err
    call mpistop('function vvol_d1 not converged',0)

  end function vvol_d1_mah

  !===============================================================================
  function vvol_d2_mah(Ti,vtent)
  !===============================================================================
    !> Compute vvol_d2 from Ti and tentative v
    !> - MAH EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: Ti,vtent
    real(wp) :: vvol_d2_mah ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: v0,v1,err,T
    ! ----------------------------------------------------------------------------

    T= Ti*Tc
    
    ! initial guess (specific volume)
    v0= vtent/roc

    ! iterations
    do i=1,15000
       v1 = v0 - d2pdv2calc_tro_mah(T,1.0_wp/v0)/d3pdv3calc_tro_mah(T,1.0_wp/v0)

       err= abs(v1-v0)/v0
       if (err.le.tol) then
          vvol_d2_mah= v1*roc
          return
       endif

       v0=v1
    enddo

    print *,'in vvol_d2: at iteration',i,'error=',err
    call mpistop('function vvol_d2 not converged',0)

  end function vvol_d2_mah

  !===============================================================================
  function intpcalc_tro_mah(Ti,roi)
  !===============================================================================
    !> Compute pressure integral from T and rho
    !> - MAH EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: Ti,roi
    real(wp) :: intpcalc_tro_mah ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: v,vbm,ekTr,ro,T
    ! ----------------------------------------------------------------------------

    T = Ti*Tc
    ro= roi*roc
    
    ! specific volume
    v = 1.0_wp/ro
    vbm = v - bmah
    ekTr= exp(-kmh*T/Tc)

    intpcalc_tro_mah=                     rg*T*log(vbm)           &
                    - (a2mh + b2mh*T + c2mh*ekTr)/(vbm)           &
                    - (a3mh + b3mh*T + c3mh*ekTr)/(2.0_wp*vbm**2) &
                    - (a4mh                     )/(3.0_wp*vbm**3) &
                    - (       b5mh*T + c5mh*ekTr)/(4.0_wp*vbm**4)

    intpcalc_tro_mah= intpcalc_tro_mah/pc*roc

  end function intpcalc_tro_mah

  !===============================================================================
  subroutine dpcalc_tro_mah(T,ro,dpdro,dpdT)
  !===============================================================================
    !> Compute pressure derivatives w.r.t. density and temperature using T and rho
    !> - MAH EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in)  :: T,ro
    real(wp), intent(out) :: dpdro,dpdT ! OUTPUTS
    ! ----------------------------------------------------------------------------
    real(wp) :: v,vbm,ekTr,dpdv
    ! ----------------------------------------------------------------------------

    ! specific volume
    v= 1.0_wp/ro
    
    vbm = v - bmah
    ekTr= exp(-kmh*T/Tc)

    dpdv= -                             rg*T/vbm**2 &
        - 2.0_wp*(a2mh + b2mh*T + c2mh*ekTr)/vbm**3 &
        - 3.0_wp*(a3mh + b3mh*T + c3mh*ekTr)/vbm**4 &
        - 4.0_wp*(a4mh                     )/vbm**5 &
        - 5.0_wp*(     b5mh*T + c5mh*ekTr)/vbm**6

    dpdro= -dpdv*v**2

    dpdT=                        rg/vbm    &
        + (b2mh - c2mh*kmh/Tc*ekTr)/vbm**2 &
        + (b3mh - c3mh*kmh/Tc*ekTr)/vbm**3 &
        + (b5mh - c5mh*kmh/Tc*ekTr)/vbm**5

  end subroutine dpcalc_tro_mah

  !===============================================================================
  subroutine stagnation_mah(h0,s0,e0,T0,p0,ro0)
  !===============================================================================
    !> Compute stagnation quantities from static quantities & h0 & s0
    !> - MAH EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: h0,s0
    real(wp), intent(inout) :: e0,T0,p0,ro0
    ! ----------------------------------------------------------------------------
    integer :: j
    real(wp) :: err,p1
    ! ----------------------------------------------------------------------------

    p1 = p0
    err = 1.0_wp
    do j=1,100
       p0 = p1
       T0 = tcalc_ph_mah(p0,h0,ro0,T0)
       ro0 = rocalc_st_mah(s0,T0,ro0)
       e0 = h0 - p0/ro0
       p1 = pcalc_roero_mah(ro0*e0,ro0,T0)
       err = abs(p1-p0)/p0

       if (err.le.tol) return
    enddo

  end subroutine stagnation_mah


!!$  !===============================================================================
!!$  function p0calc_mah(T,p,ro,M)
!!$  !===============================================================================
!!$    !> Compute total pressure from static quantities
!!$    !> - MAH EOS -
!!$  !===============================================================================
!!$    implicit none
!!$    ! ----------------------------------------------------------------------------
!!$    real(wp), intent(in) :: T,M,ro,p
!!$    real(wp) :: p0calc_mah ! OUTPUT
!!$    ! ----------------------------------------------------------------------------
!!$    integer :: j
!!$    real(wp) :: s0,ro0,e0,p0,h0
!!$    real(wp) :: err,c2mh,e,T0,p1,h
!!$    ! ----------------------------------------------------------------------------
!!$
!!$    c2mh = c2calc_tro_mah(T,ro)
!!$    e = ecalc_tro_mah(T,ro)
!!$    h = e + p/ro
!!$
!!$    s0 = scalc_tro_mah(T,ro)
!!$    h0 = h + (M**2*c2mh)/2
!!$    T0 = T
!!$    ro0=ro
!!$    p1 = p
!!$    err = 1.0_wp
!!$    do j=1,100
!!$       p0 = p1
!!$       T0 = tcalc_ph_mah(p0,h0,ro0,T0)
!!$       ro0 = rocalc_st_mah(s0,T0,ro0)
!!$       e0 = h0 - p0/ro0
!!$       p1 = pcalc_roero_mah(ro0*e0,ro0,T)
!!$       err = abs(p1-p0)/p0
!!$
!!$       if (err.le.tol) then
!!$          p0calc_mah= p0
!!$          return
!!$       endif
!!$    enddo
!!$    
!!$  end function p0calc_mah
  
end module mod_eos_mah
