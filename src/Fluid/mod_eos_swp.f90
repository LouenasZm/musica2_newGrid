!=================================================================================
module mod_ineos_swp
!=================================================================================
  !> Module to share constants of SPAN-WAGNER-polar EoS [for inlining purposes]
!=================================================================================
  use precision
  implicit none
  ! ------------------------------------------------------------------------------
  real(wp), parameter :: tol=1.e-5_wp
  real(wp) :: zci,Tci,roci,pcvc
  ! ------------------------------------------------------------------------------
end module mod_ineos_swp

!=================================================================================
module mod_eos_swp
!=================================================================================
  !> Module to define subroutines for polar SPAN-WAGNER Equation of State (EoS)
!=================================================================================
  use mod_fluid
  use mod_ineos_swp
  use warnstop
  implicit none
  ! ------------------------------------------------------------------------------
  real(wp), private :: pci
  real(wp), private :: s_crit
  ! ------------------------------------------------------------------------------

contains

  !===============================================================================
  pure function psir(tau,del)
  !===============================================================================
    !>  function Ψr
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: tau,del
    real(wp) :: psir ! OUTPUT
    ! ----------------------------------------------------------------------------

    psir = n01*del   *tau**(0.25_wp)                 &
         + n02*del   *tau**(1.25_wp)                 &
         + n03*del   *tau**(1.5_wp)                  &
         + n04*del**3*tau**(0.25_wp)                 &
         + n05*del**7*tau**(0.875_wp)                &
         + n06*del   *tau**(2.375_wp)*exp(-del)      &
         + n07*del**2*tau**(2.0_wp)  *exp(-del)      &
         + n08*del**5*tau**(2.125_wp)*exp(-del)      &
         + n09*del   *tau**(3.5_wp)  *exp(-del**2)   &
         + n10*del   *tau**(6.5_wp)  *exp(-del**2)   &
         + n11*del**4*tau**(4.75_wp) *exp(-del**2)   &
         + n12*del**2*tau**(12.5_wp) *exp(-del**3)

  end function psir

  !===============================================================================
  pure function psi0(tau,del)
  !===============================================================================
    !>  function Ψ₀
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: tau,del
    real(wp) :: psi0 ! OUTPUT
    ! ----------------------------------------------------------------------------

    psi0= (eta1 - 1.0_wp)*(-tau-log(1.0_wp/tau))              &
         + eta2*Tc   *(-1.0_wp/(2.0_wp *tau)    - tau/2.0_wp) &
         + eta3*Tc**2*(-1.0_wp/(6.0_wp *tau**2) - tau/3.0_wp) &
         + eta4*Tc**3*(-1.0_wp/(12.0_wp*tau**3) - tau/4.0_wp) + log(del)

  end function psi0

  !===============================================================================
  pure function d1psir_del1(tau,del)
  !===============================================================================
    !>  function  ∂¹ Ψr/∂δ¹
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: tau,del
    real(wp) :: d1psir_del1 ! OUTPUT
    ! ----------------------------------------------------------------------------

    d1psir_del1=                                                                 &
         + n01              *tau**(0.25_wp )                                     &
         + n02              *tau**(1.25_wp )                                     &
         + n03              *tau**(1.5_wp  )                                     &
         + n04*3.0_wp*del**2*tau**(0.25_wp )                                     &
         + n05*7.0_wp*del**6*tau**(0.875_wp)                                     &
         + n06              *tau**(2.375_wp)*exp(-del   )*(1.0_wp-     del   )   &
         + n07       *del   *tau**(2       )*exp(-del   )*(2.0_wp-     del   )   &
         + n08       *del**4*tau**(2.125_wp)*exp(-del   )*(5.0_wp-     del   )   &
         + n09              *tau**(3.5_wp  )*exp(-del**2)*(1.0_wp-2.0_wp*del**2) &
         + n10              *tau**(6.5_wp  )*exp(-del**2)*(1.0_wp-2.0_wp*del**2) &
         + n11       *del**3*tau**(4.75_wp )*exp(-del**2)*(4.0_wp-2.0_wp*del**2) &
         + n12       *del   *tau**(12.5_wp )*exp(-del**3)*(2.0_wp-3.0_wp*del**3)

  end function d1psir_del1

  !===============================================================================
  pure function d2psir_del2(tau,del)
  !===============================================================================
    !>  function  ∂² Ψr/∂δ²
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: tau,del
    real(wp) :: d2psir_del2 ! OUTPUT
    ! ----------------------------------------------------------------------------
 
    d2psir_del2=                                                                                 &
        + n04*6.0_wp *del   *tau**(0.25_wp )                                                     &
        + n05*42.0_wp*del**5*tau**(0.875_wp)                                                     &
        + n06               *tau**(2.375_wp)*exp(-del   )*(-2.0_wp +        del)                 &
        + n07               *tau**(2      ) *exp(-del   )*( 2.0_wp -4.0_wp *del+del**2)          &
        + n08        *del**3*tau**(2.125_wp)*exp(-del   )*( 20.0_wp-10.0_wp*del+del**2)          &
        + n09        *del   *tau**(3.5_wp  )*exp(-del**2)*(-6.0_wp +4.0_wp *del**2)              &
        + n10        *del   *tau**(6.5_wp  )*exp(-del**2)*(-6.0_wp +4.0_wp *del**2)              &
        + n11        *del**2*tau**(4.75_wp )*exp(-del**2)*( 12.0_wp-18.0_wp*del**2+4.0_wp*del**4)&
        + n12               *tau**(12.5_wp )*exp(-del**3)*( 2.0_wp -18.0_wp*del**3+9.0_wp*del**6)

  end function d2psir_del2

  !===============================================================================
  pure function d3psir_del3(tau,del)
  !===============================================================================
    !>  function ∂³ Ψr/∂δ³
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: tau,del
    real(wp) :: d3psir_del3 ! OUTPUT
    ! ----------------------------------------------------------------------------

    d3psir_del3=                                                                                       &
      + n04       *tau**(0.25_wp ) *6.0_wp                                                             &
      + n05*del**4*tau**(0.875_wp)*210.0_wp                                                            &
      + n06       *tau**(2.375_wp)*exp(-del   )*( 3.0_wp -       del)                                  &
      + n07       *tau**(2.0_wp  )*exp(-del   )*(-6.0_wp +6.0_wp  *del-del**2)                         &
      + n08*del**2*tau**(2.125_wp)*exp(-del   )*( 60.0_wp-60.0_wp *del+15.0_wp*del**2-1.0_wp*del**3)   &
      + n09       *tau**(3.5_wp  )*exp(-del**2)*(-6.0_wp +24.0_wp *del**2-8.0_wp*del**4)               &
      + n10       *tau**(6.5_wp  )*exp(-del**2)*(-6.0_wp +24.0_wp *del**2-8.0_wp*del**4)               &
      + n11*del   *tau**(4.75_wp )*exp(-del**2)*( 24.0_wp-96.0_wp *del**2+60.0_wp*del**4-8.0_wp*del**6)&
      + n12*del**2*tau**(12.5_wp )*exp(-del**3)*(-60.0_wp+108.0_wp*del**3-27.0_wp*del**6)
    
  end function d3psir_del3


  !===============================================================================
  pure function d4psir_del4(tau,del)
  !===============================================================================
    !>  function ∂⁴ Ψr/∂δ⁴
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: tau,del
    real(wp) :: d4psir_del4 ! OUTPUT
    ! ----------------------------------------------------------------------------

    d4psir_del4=                                                                           &
         + n05*del**3*tau**(0.875_wp)*840.0_wp                                             &
         + n06       *tau**(2.375_wp)*exp(-del   )*(-4.0_wp + del)                         &
         + n07       *tau**(2.0_wp)  *exp(-del   )*(12.0_wp-8.0_wp*del+del**2)             &
         + n08*del   *tau**(2.125_wp)*exp(-del   )                                         &
         *(120.0_wp-240.0_wp*del+120.0_wp*del**2-20.0_wp*del**3+del**4)                    &
         + n09*del   *tau**(3.5_wp)  *exp(-del**2)*(60.0_wp-80.0_wp*del**2+16.0_wp*del**4) &
         + n10*del   *tau**(6.5_wp)  *exp(-del**2)*(60.0_wp-80.0_wp*del**2+16.0_wp*del**4) &
         + n11       *tau**(4.75_wp) *exp(-del**2)                                         &
         *(24.0_wp-336.0_wp*del**2+492.0_wp*del**4-176.0_wp*del**6+16.0_wp*del**8)         &
         + n12*del   *tau**(12.5_wp) *exp(-del**3)                                         &
         *(-120.0_wp+720.0_wp*del**3-540.0_wp*del**6+81.0_wp*del**9)

  end function d4psir_del4

  !===============================================================================
  pure function d1psir_tau1(tau,del)
  !===============================================================================
    !>  function ∂¹ Ψr/ ∂τ¹
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: tau,del
    real(wp) :: d1psir_tau1 ! OUTPUT
    ! ----------------------------------------------------------------------------

    d1psir_tau1=                                             &
         + n01*0.25_wp *del   *tau**(-0.75_wp )              &
         + n02*1.25_wp *del   *tau**( 0.25_wp )              &
         + n03*1.5_wp  *del   *tau**( 0.5_wp  )              &
         + n04*0.25_wp *del**3*tau**(-0.75_wp )              &
         + n05*0.875_wp*del**7*tau**(-0.125_wp)              &
         + n06*2.375_wp*del   *tau**( 1.375_wp)*exp(-del)    &
         + n07*2.0_wp  *del**2*tau             *exp(-del)    &
         + n08*2.125_wp*del**5*tau**( 1.125_wp)*exp(-del)    &
         + n09*3.5_wp  *del   *tau**( 2.5_wp)  *exp(-del**2) &
         + n10*6.5_wp  *del   *tau**( 5.5_wp)  *exp(-del**2) &
         + n11*4.75_wp *del**4*tau**( 3.75_wp) *exp(-del**2) &
         + n12*12.5_wp *del**2*tau**( 11.5_wp) *exp(-del**3)

  end function d1psir_tau1

  !===============================================================================
  pure  function d2psir_tau2(tau,del)
  !===============================================================================
    !>  function ∂² Ψr/∂τ²
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: tau,del
    real(wp) :: d2psir_tau2 ! OUTPUT
    ! ----------------------------------------------------------------------------

    d2psir_tau2=                                                         &
         + n01*0.25_wp *(-0.75_wp) *del   *tau**(-1.75_wp)               &
         + n02*1.25_wp *( 0.25_wp) *del   *tau**(-0.75_wp)               &
         + n03*1.5_wp  *( 0.5_wp)  *del   *tau**(-0.5_wp)                &
         + n04*0.25_wp *(-0.75_wp) *del**3*tau**(-1.75_wp)               &
         + n05*0.875_wp*(-0.125_wp)*del**7*tau**(-1.125_wp)              &
         + n06*2.375_wp*( 1.375_wp)*del   *tau**( 0.375_wp)*exp(-del)    &
         + n07*2.0_wp              *del**2                 *exp(-del)    &
         + n08*2.125_wp*( 1.125_wp)*del**5*tau**( 0.125_wp)*exp(-del)    &
         + n09*3.5_wp  *( 2.5_wp)  *del   *tau**( 1.5_wp)  *exp(-del**2) &
         + n10*6.5_wp  *( 5.5_wp)  *del   *tau**( 4.5_wp)  *exp(-del**2) &
         + n11*4.75_wp *( 3.75_wp) *del**4*tau**( 2.75_wp) *exp(-del**2) &
         + n12*12.5_wp *( 11.5_wp) *del**2*tau**( 10.5_wp) *exp(-del**3)

  end function d2psir_tau2

  !===============================================================================
  pure function d3psir_tau3(tau,del)
  !===============================================================================
    !>  function ∂³ Ψr/∂τ³
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: tau,del
    real(wp) :: d3psir_tau3 ! OUTPUT
    ! ----------------------------------------------------------------------------

    d3psir_tau3=                                                                     &
         + n01*0.25_wp *(-0.75_wp )*(-1.75_wp )*del   *tau**(-2.75_wp )              &
         + n02*1.25_wp *( 0.25_wp )*(-0.75_wp )*del   *tau**(-1.75_wp )              &
         + n03*1.5_wp  *( 0.5_wp  )*(-0.5_wp  )*del   *tau**(-1.5_wp  )              &
         + n04*0.25_wp *(-0.75_wp )*(-1.75_wp )*del**3*tau**(-2.75_wp )              &
         + n05*0.875_wp*(-0.125_wp)*(-1.125_wp)*del**7*tau**(-2.125_wp)              &
         + n06*2.375_wp*( 1.375_wp)*( 0.375_wp)*del   *tau**(-0.625_wp)*exp(-del   ) &
         + n08*2.125_wp*( 1.125_wp)*( 0.125_wp)*del**5*tau**(-0.875_wp)*exp(-del   ) &
         + n09*3.5_wp  *( 2.5_wp  )*( 1.5_wp  )*del   *tau**( 0.5_wp  )*exp(-del**2) &
         + n10*6.5_wp  *( 5.5_wp  )*( 4.5_wp  )*del   *tau**( 3.5_wp  )*exp(-del**2) &
         + n11*4.75_wp *( 3.75_wp )*( 2.75_wp )*del**4*tau**( 1.75_wp )*exp(-del**2) &
         + n12*12.5_wp *( 11.5_wp )*( 10.5_wp )*del**2*tau**( 9.5_wp  )*exp(-del**3)

  end function d3psir_tau3

  !===============================================================================
  pure function d2psir_tau1del1(tau,del)
  !===============================================================================
    !>  function  ∂² Ψr/∂τ ∂δ
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: tau,del
    real(wp) :: d2psir_tau1del1 ! OUTPUT
    ! ----------------------------------------------------------------------------

    d2psir_tau1del1=                                                                       &
         + n01*0.25_wp               *tau**(-0.75_wp )                                     &
         + n02*1.25_wp               *tau**( 0.25_wp )                                     &
         + n03*1.5_wp                *tau**( 0.5_wp  )                                     &
         + n04*3.0_wp*0.25_wp *del**2*tau**(-0.75_wp )                                     &
         + n05*7.0_wp*0.875_wp*del**6*tau**(-0.125_wp)                                     &
         + n06*2.375_wp              *tau**( 1.375_wp)*exp(-del   )*(1.0_wp-del        )   &
         + n07*2.0_wp         *del   *tau             *exp(-del   )*(2.0_wp-del        )   &
         + n08*2.125_wp       *del**4*tau**( 1.125_wp)*exp(-del   )*(5.0_wp-del        )   &
         + n09*3.5_wp                *tau**( 2.5_wp  )*exp(-del**2)*(1.0_wp-2.0_wp*del**2) &
         + n10*6.5_wp                *tau**( 5.5_wp  )*exp(-del**2)*(1.0_wp-2.0_wp*del**2) &
         + n11*4.75_wp        *del**3*tau**( 3.75_wp )*exp(-del**2)*(4.0_wp-2.0_wp*del**2) &
         + n12*12.5_wp        *del   *tau**( 11.5_wp )*exp(-del**3)*(2.0_wp-3.0_wp*del**3)

  end function d2psir_tau1del1

  !===============================================================================
  pure function d3psir_tau2del1(tau,del)
  !===============================================================================
    !>  function ∂³ Ψr/∂τ²∂δ
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: tau,del
    real(wp) :: d3psir_tau2del1 ! OUTPUT
    ! ----------------------------------------------------------------------------

    d3psir_tau2del1=                                                                                  &
         + n01*0.25_wp        *(-0.75_wp )       *tau**(-1.75_wp )                                    &
         + n02*1.25_wp        *( 0.25_wp )       *tau**(-0.75_wp )                                    &
         + n03*1.5_wp         *( 0.5_wp  )       *tau**(-0.5_wp  )                                    &
         + n04*3.0_wp*0.25_wp *(-0.75_wp )*del**2*tau**(-1.75_wp )                                    &
         + n05*7.0_wp*0.875_wp*(-0.125_wp)*del**6*tau**(-1.125_wp)                                    &
         + n06*2.375_wp       *( 1.375_wp)       *tau**( 0.375_wp)*exp(-del   )*(1.0_wp-     del   )  &
         + n07*2.0_wp                     *del                    *exp(-del   )*(2.0_wp-     del   )  &
         + n08*2.125_wp       *( 1.125_wp)*del**4*tau**( 0.125_wp)*exp(-del   )*(5.0_wp-     del   )  &
         + n09*3.5_wp         *( 2.5_wp  )       *tau**( 1.5_wp  )*exp(-del**2)*(1.0_wp-2.0_wp*del**2)&
         + n10*6.5_wp         *( 5.5_wp  )       *tau**( 4.5_wp  )*exp(-del**2)*(1.0_wp-2.0_wp*del**2)&
         + n11*4.75_wp        *( 3.75_wp )*del**3*tau**( 2.75_wp )*exp(-del**2)*(4.0_wp-2.0_wp*del**2)&
         + n12*12.5_wp        *( 11.5_wp )*del   *tau**( 10.5_wp )*exp(-del**3)*(2.0_wp-3.0_wp*del**3)

  end function d3psir_tau2del1

  !===============================================================================
  pure function d3psir_tau1del2(tau,del)
  !===============================================================================
    !>  function ∂³ Ψr/∂τ∂δ²
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: del,tau
    real(wp) :: d3psir_tau1del2 ! OUTPUT
    ! ----------------------------------------------------------------------------

    d3psir_tau1del2=                                                                                     &
      + n04*6.0_wp *0.25_wp *del   *tau**(-0.75_wp )                                                     &
      + n05*42.0_wp*0.875_wp*del**5*tau**(-0.125_wp)                                                     &
      + n06*2.375_wp               *tau**( 1.375_wp)*exp(-del   )*(-2.0_wp +del)                         &
      + n07*2.0_wp                 *tau             *exp(-del   )*( 2.0_wp -4.0_wp *del+del**2        )  &
      + n08*2.125_wp        *del**3*tau**( 1.125_wp)*exp(-del   )*( 20.0_wp-10.0_wp*del+del**2        )  &
      + n09*3.5_wp          *del   *tau**( 2.5_wp  )*exp(-del**2)*(-6.0_wp +4.0_wp *del**2            )  &
      + n10*2.0_wp *6.5_wp  *del   *tau**( 5.5_wp  )*exp(-del**2)*(-3.0_wp +2.0_wp *del**2            )  &
      + n11*4.75_wp         *del**2*tau**( 3.75_wp )*exp(-del**2)*( 12.0_wp-18.0_wp*del**2+4.0_wp*del**4)&
      + n12*12.5_wp                *tau**( 11.5_wp )*exp(-del**3)*( 2.0_wp -18.0_wp*del**3+9.0_wp*del**6)

  end function d3psir_tau1del2

  !===============================================================================
  pure function d1psi0_del1(del)
  !===============================================================================
    !>  function ∂¹ Ψ₀/∂δ¹
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: del
    real(wp) :: d1psi0_del1 ! OUTPUT
    ! ----------------------------------------------------------------------------

    d1psi0_del1 = 1.0_wp/del

  end function d1psi0_del1

  !===============================================================================
  pure function d1psi0_tau1(tau)
  !===============================================================================
    !>  function ∂¹ Ψ₀/∂τ¹
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: tau
    real(wp) :: d1psi0_tau1 ! OUTPUT
    ! ----------------------------------------------------------------------------

    d1psi0_tau1= (eta1-1.0_wp)*(1.0_wp/tau - 1.0_wp)         &
                + eta2/2.0_wp*Tc   *(1.0_wp/tau**2 - 1.0_wp) &
                + eta3/3.0_wp*Tc**2*(1.0_wp/tau**3 - 1.0_wp) &
                + eta4/4.0_wp*Tc**3*(1.0_wp/tau**4 - 1.0_wp)

  end function d1psi0_tau1

  !===============================================================================
  pure function d2psi0_tau2(tau)
  !===============================================================================
    !>  function ∂² Ψ₀/∂τ²
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: tau
    real(wp) :: d2psi0_tau2 ! OUTPUT
    ! ----------------------------------------------------------------------------

    d2psi0_tau2=-(eta1-1.0_wp)/tau**2         &
                + eta2*Tc   *(-1.0_wp/tau**3) &
                + eta3*Tc**2*(-1.0_wp/tau**4) &
                + eta4*Tc**3*(-1.0_wp/tau**5)

  end function d2psi0_tau2

  !===============================================================================
  pure function d3psi0_tau3(tau)
  !===============================================================================
    !>  function ∂³ Ψ₀/∂τ³
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: tau
    real(wp) :: d3psi0_tau3 ! OUTPUT
    ! ----------------------------------------------------------------------------

    d3psi0_tau3= 2.0_wp*(eta1-1.0_wp     )/tau**3  &
               + 3.0_wp*eta2*Tc   *(1.0_wp/tau**4) &
               + 4.0_wp*eta3*Tc**2*(1.0_wp/tau**5) &
               + 5.0_wp*eta4*Tc**3*(1.0_wp/tau**6)
    
  end function d3psi0_tau3

  !===============================================================================
  subroutine init_eos_swp
  !===============================================================================
    !> Initializations of gas coefficients
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    call read_eos

    rg = 8.31451d3/pmol
    zc = pc/(roc*rg*Tc)
    zci = 1.0_wp/zc
    Tci = 1.0_wp/Tc
    pci = 1.0_wp/pc
    roci = 1.0_wp/roc
    pcvc = pc/roc

    ! for compatibility with cubic EoS that change roc value
    roc0=roc

    ! critical entropy 
    ! ----------------
    s_crit = 0.0_wp ! Initialization
    s_crit = scalc_tro_swp(Tc,roc)

  end subroutine init_eos_swp

  !===============================================================================
  function avcalc_tro_swp(T,ro)
  !===============================================================================
    !> Compute isobaric expansion coefficient from T and rho
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: avcalc_tro_swp ! OUTPUT
    ! ----------------------------------------------------------------------------
    !real(wp) :: dpdv,dpdT,v,vbm,ekTr
    ! ----------------------------------------------------------------------------

    call mpistop('procedure avcalc_tro_swp not written yet !!',0)

    avcalc_tro_swp= 1.0_wp

  end function avcalc_tro_swp

  !===============================================================================
  function c2calc_tro_swp(Ti,roi)
  !===============================================================================
    !> Compute speed of sound from T and rho
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: Ti,roi
    real(wp) :: c2calc_tro_swp ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: T,ro,tau,del,dpdv,dpdT,cv
    ! ----------------------------------------------------------------------------

    T = Ti*Tci
    ro= roi*roci

    tau= 1.0_wp/T
    del= ro

    cv=-tau**2*zci*( d2psi0_tau2(tau) &
                   + d2psir_tau2(tau,del) )

    dpdv=-del**2/tau*zci*( 1.0_wp+2.0_wp*del*d1psir_del1(tau,del) &
                                    + del**2*d2psir_del2(tau,del) )

    dpdT= del*zci*(1.0_wp+del    *d1psir_del1(tau,del)    &
                         -del*tau*d2psir_tau1del1(tau,del))

    c2calc_tro_swp= 1.0_wp/ro**2*(T/cv*dpdT**2 - dpdv)

    c2calc_tro_swp= c2calc_tro_swp*pcvc

  end function c2calc_tro_swp

  !===============================================================================
  function cpcalc_tro_swp(Ti,roi)
  !===============================================================================
    !> Compute heat capacity at constant pressure
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: Ti,roi
    real(wp) :: cpcalc_tro_swp ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: T,ro,tau,del,cv,dpdv,dpdT
    ! ----------------------------------------------------------------------------

    T = Ti*Tci
    ro= roi*roci

    tau= 1.0_wp/T
    del= ro

    cv=-tau**2*zci*(d2psi0_tau2(tau)+d2psir_tau2(tau,del))

    dpdv=-del**2/tau*zci*( 1.0_wp+2.0_wp*del*d1psir_del1(tau,del) &
                                    + del**2*d2psir_del2(tau,del) )

    dpdT= del*zci*( 1.0_wp+del*d1psir_del1(tau,del) &
                 - del*tau*d2psir_tau1del1(tau,del) )

    cpcalc_tro_swp= cv - T*dpdT**2/dpdv

    cpcalc_tro_swp= cpcalc_tro_swp*zc*rg

  end function cpcalc_tro_swp

  !===============================================================================
  function cvcalc_tro_swp(Ti,roi)
  !===============================================================================
    !> Compute heat capacity at constant volume
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: Ti,roi
    real(wp) :: cvcalc_tro_swp ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: T,ro,tau,del
    ! ----------------------------------------------------------------------------

    T = Ti*Tci
    ro= roi*roci

    tau= 1.0_wp/T
    del= ro

    cvcalc_tro_swp= -tau**2*zci*(d2psi0_tau2(tau)+d2psir_tau2(tau,del))

    cvcalc_tro_swp= cvcalc_tro_swp*zc*rg

  end function cvcalc_tro_swp

  !===============================================================================
  function dpdicalc_tro_swp(Ti,roi)
  !===============================================================================
    !> Compute pressure derivative w.r.t temperature
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: Ti,roi
    real(wp) :: dpdicalc_tro_swp ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: T,ro,dpdT,cv,tau,del
    ! ----------------------------------------------------------------------------

    T = Ti*Tci
    ro= roi*roci

    tau= 1.0_wp/T
    del= ro

    cv= -tau**2*zci*(d2psi0_tau2(tau)+d2psir_tau2(tau,del))
    
    dpdT= del*zci*(1.0_wp + del*d1psir_del1(tau,del) &
                  - del*tau*d2psir_tau1del1(tau,del) )

    dpdicalc_tro_swp= dpdT/cv  ! dimensionless

  end function dpdicalc_tro_swp

  !===============================================================================
  function dpdTcalc_tro_swp(Ti,roi)
  !===============================================================================
    !> Compute pressure derivative w.r.t temperature
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: Ti,roi
    real(wp) :: dpdTcalc_tro_swp ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: T,ro,tau,del
    real(wp) :: dpdT
    ! ----------------------------------------------------------------------------

    T = Ti*Tci
    ro= roi*roci

    tau= 1.0_wp/T
    del= ro

    dpdT= del*zci*( 1.0_wp+del*d1psir_del1(tau,del) &
                 - del*tau*d2psir_tau1del1(tau,del) )  ! dimensionless

  end function dpdTcalc_tro_swp

!!$  !===============================================================================
!!$  function dpdvcalc_tro_swp(Ti,roi)
!!$  !===============================================================================
!!$    !> Compute first-order pressure derivative w.r.t volume
!!$    !> - SWP EOS -
!!$  !===============================================================================
!!$    implicit none
!!$    ! ----------------------------------------------------------------------------
!!$    real(wp), intent(in) :: Ti,roi
!!$    real(wp) :: dpdvcalc_tro_swp ! OUTPUT
!!$    ! ----------------------------------------------------------------------------
!!$    real(wp) :: T,ro,tau,del,cv,dpdv,dpdT
!!$    ! ----------------------------------------------------------------------------
!!$
!!$    T = Ti*Tci
!!$    ro= roi*roci
!!$
!!$    tau= 1.0_wp/T
!!$    del= ro
!!$
!!$    dpdvcalc_tro_swp=-del**2/tau*zci*( 1.0_wp+2.0_wp*del*d1psir_del1(tau,del) &
!!$                                                + del**2*d2psir_del2(tau,del) )
!!$    !Dimensionalization
!!$    dpdvcalc_tro_swp= dpdvcalc_tro_swp/(roc*pc)
!!$
!!$  end function dpdvcalc_tro_swp
!!$
!!$  !===========================================================================
!!$  !                  PRESSURE 2ND DERIVATIVE W.R.T. VOLUME(T,RO)
!!$  !===========================================================================
!!$  !===============================================================================
!!$  function d2pdv2calc_tro_swp(Ti,roi)
!!$  !===============================================================================
!!$    !> Compute second-order pressure derivative w.r.t volume
!!$    !> - SWP EOS -
!!$  !===============================================================================
!!$    implicit none
!!$    ! ----------------------------------------------------------------------------
!!$    real(wp), intent(in) :: Ti,roi
!!$    real(wp) :: d2pdv2calc_tro_swp ! OUTPUT
!!$    ! ----------------------------------------------------------------------------
!!$    real(wp) :: T,ro,tau,del,cv,dpdv,dpdT
!!$    ! ----------------------------------------------------------------------------
!!$
!!$    T = Ti*Tci
!!$    ro= roi*roci
!!$
!!$    tau= 1.0_wp/T
!!$    del= ro
!!$
!!$    d2pdv2calc_tro_swp=-del**2/tau*zci*( 1.0_wp+2.0_wp*del*d1psir_del1(tau,del) &
!!$                                                  + del**2*d2psir_del2(tau,del) )
!!$    ! Dimensionalization
!!$    d2pdv2calc_tro _swp= d2pdv2calc_tro_swp/(roc*pc)
!!$
!!$  end function d2pdv2calc_tro_swp
!!$
!!$  !===============================================================================
!!$  function d3pdv3calc_tro_swp(Ti,roi)
!!$  !===============================================================================
!!$    !> Compute third-order pressure derivative w.r.t volume
!!$    !> - SWP EOS -
!!$  !===============================================================================
!!$    implicit none
!!$    ! ----------------------------------------------------------------------------
!!$    real(wp), intent(in) :: Ti,roi
!!$    real(wp) :: d3pdv3calc_tro_swp ! OUTPUT
!!$    ! ----------------------------------------------------------------------------
!!$    real(wp) :: T,ro,tau,del,cv,dpdv,dpdT
!!$    ! ----------------------------------------------------------------------------
!!$    
!!$    T = Ti*Tci
!!$    ro= roi*roci
!!$
!!$    tau= 1.0_wp/T
!!$    del= ro
!!$
!!$    d3pdv3calc_tro_swp=-del**2/tau*zci*( 1.0_wp+2.0_wp*del*d1psir_del1(tau,del) &
!!$                                                  + del**2*d2psir_del2(tau,del) )
!!$    ! Dimensionalization
!!$    d3pdv3calc_tro_swp= d3pdv3calc_tro_swp/(roc*pc)
!!$
!!$  end function d3pdv3calc_tro_swp

  !===============================================================================
  function ecalc_pro_swp(p,ro,Ttent)
  !===============================================================================
    !> Compute internal energy from p and rho
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: p,ro,Ttent
    real(wp) :: ecalc_pro_swp ! OUTPUT
    ! ----------------------------------------------------------------------------
    ! integer :: i
    !real(wp) :: T,T1,v,err
    ! ----------------------------------------------------------------------------

    call mpistop('procedure ecalc_pro_swp not written yet !!',0)

    ecalc_pro_swp= 1.0_wp

    !print *,'in ecalc_pro: at iteration',i,'error=',err
    !print *,'p:',p,'ro:',ro
    !ecalc_pro_swp= 0.0_wp
    !call mpistop('function ecalc_pro not converged',0)

  end function ecalc_pro_swp

  !===============================================================================
  function ecalc_tro_swp(Ti,roi)
  !===============================================================================
    !> Compute internal energy from T and rho
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: Ti,roi
    real(wp) :: ecalc_tro_swp ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: T,ro,tau,del
    ! ----------------------------------------------------------------------------

    T = Ti*Tci
    ro= roi*roci

    tau= 1.0_wp/T
    del= ro

    ecalc_tro_swp= 100.0_wp + zci*(d1psir_tau1(tau,del) + d1psi0_tau1(tau))

    ecalc_tro_swp= ecalc_tro_swp*zc*pcvc
    
  end function ecalc_tro_swp

  !===============================================================================
  function dedTcalc_tro_swp(Ti,roi)
  !===============================================================================
    !> Compute internal energy from T and rho
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: Ti,roi
    real(wp) :: dedTcalc_tro_swp ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: T,ro,tau,del
    ! ----------------------------------------------------------------------------

    T = Ti*Tci
    ro= roi*roci

    tau= 1.0_wp/T
    del= ro

    dedTcalc_tro_swp= -zci*(d2psir_tau2(tau,del)+d2psi0_tau2(tau))/T**2 ! dimensionless

    !! check dimensionality
    !!dedTcalc_tro_swp= dedTcalc_tro_swp*zc*pcvc/Tc
    !!dedTcalc_tro_swp= dedTcalc_tro_swp*zc*zc

  end function dedTcalc_tro_swp

  !===============================================================================
  function gcalc_tro_swp(Ti,roi)
  !===============================================================================
    !> Compute fundamental derivative of gas dynamics from T and ro
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: Ti,roi
    real(wp) :: gcalc_tro_swp ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: T,ro,tau,del
    real(wp) :: cv,dpdv,dpdT,d2pdT2,d2pdv2,d2pdvdT,dcvdT,c2calc
    ! ----------------------------------------------------------------------------

    T = Ti*Tci
    ro= roi*roci

    tau= 1.0_wp/T
    del= ro

    cv= -tau**2*zci*(d2psi0_tau2(tau)+d2psir_tau2(tau,del))

    dpdv= -del**2/tau*zci*( 1.0_wp+2.0_wp*del*d1psir_del1(tau,del) &
                                     + del**2*d2psir_del2(tau,del) )

    dpdT= del*zci*( 1.0_wp+del*d1psir_del1(tau,del) &
                 - del*tau*d2psir_tau1del1(tau,del) )

    d2pdT2= del**2*tau**3*zci*d3psir_tau2del1(tau,del)

    d2pdv2=        del**3/tau*zci *( 2.0_wp    &
          + 6.0_wp*del   *d1psir_del1(tau,del) &
          + 6.0_wp*del**2*d2psir_del2(tau,del) &
          +        del**3*d3psir_del3(tau,del) )

    d2pdvdT=        del**2*zci*( -1.0_wp                &
           + 2.0_wp*del   *tau*d2psir_tau1del1(tau,del) &
           +        del**2*tau*d3psir_tau1del2(tau,del) &
           - 2.0_wp*del       *d1psir_del1(tau,del)     &
           -        del**2    *d2psir_del2(tau,del)     )

    dcvdT= tau**3*zci * (2.0_wp*(d2psi0_tau2(tau) + d2psir_tau2(tau,del)) + &
           tau*(d3psi0_tau3(tau) + d3psir_tau3(tau,del))  )

    c2calc= 1.0_wp/ro**2*(T/cv*dpdT**2 - dpdv)

    gcalc_tro_swp= 1.0_wp/(2.0_wp*ro**3*c2calc)*( d2pdv2-3.0_wp*T/cv*dpdT*d2pdvdT+ &
                         (T/cv*dpdT)**2*(3.0_wp*d2pdT2+dpdT/T*(1.0_wp-T/cv*dcvdT)) )

  end function gcalc_tro_swp

  !===============================================================================
  function pcalc_roero_swp(roei,roi,Ttent)
  !===============================================================================
    !> Compute pressure from rhoe and rho
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: roei,roi,Ttent
    real(wp) :: pcalc_roero_swp ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: e,ro,tau,del,err,T,T1,dedT,func,ecalc
    ! ----------------------------------------------------------------------------
    
    ! initial guess
    T = Ttent*Tci
    ro= roi*roci
    e = roei/roi/(zc*pcvc)

    del=ro

    ! Newton's algorithm
    ! ------------------
    do i=1,20

       tau= 1.0_wp/T

       ! function & derivative
       ecalc= 100.0_wp + zci*(d1psir_tau1(tau,del) + d1psi0_tau1(tau))
       dedT =          - zci*(d2psir_tau2(tau,del) + d2psi0_tau2(tau))/T**2
       func= ecalc - e
       
       ! update solution
       T1 = T - func/dedT

       err = abs(T1-T)/T
       if (err.le.tol) then
          pcalc_roero_swp= del*zci/tau*(1.0_wp + del*d1psir_del1(tau,del))
          pcalc_roero_swp= pcalc_roero_swp*pc
          return
       endif
       T=T1
    enddo

    print *,'in pcalc_roero: at iteration',i,'error=',err
    print *,'roe:',roei,'ro:',roi
    pcalc_roero_swp= 0.0_wp
    call mpistop('function pcalc_roero not converged',0)

  end function pcalc_roero_swp

  !===============================================================================
  function pcalc_tro_swp(Ti,roi)
  !===============================================================================
    !> Compute pressure from T and rho
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: Ti,roi
    real(wp) :: pcalc_tro_swp ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: T,ro,tau,del
    ! ----------------------------------------------------------------------------

    T = Ti*Tci
    ro= roi*roci

    tau= 1.0_wp/T
    del= ro

    pcalc_tro_swp= del*zci/tau*(1.0_wp+del*d1psir_del1(tau,del))
    pcalc_tro_swp= pcalc_tro_swp*pc

  end function pcalc_tro_swp

  !===============================================================================
  function rocalc_ep_swp(e,p,Ttent)
  !===============================================================================
    !> Compute density from internal energy e and p
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: e,p,Ttent
    real(wp) :: rocalc_ep_swp ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: T,T1,ro,err
    ! ----------------------------------------------------------------------------

    ! initial guess
    T1= Ttent
    ro= roc

    do i=1,50
       T = T1
       ro= rocalc_pt_swp(p,T,ro)
       T1= Tcalc_roero_swp(ro*e,ro,T)
       
       err= abs(T1-T)/T
       if (err.lt.1.e-10_wp) then
          rocalc_ep_swp= ro
          return
       endif
    enddo

    print *,'in rocalc_ep: at iteration',i,'error=',err
    call mpistop('function rocalc_ep not converged',0)

  end function rocalc_ep_swp

  !===============================================================================
  function rocalc_ps_swp(p,s,Ttent)
  !===============================================================================
    !> Compute density from p and entropy s
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: p,s,Ttent
    real(wp) :: rocalc_ps_swp ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: T,T1,ro,err
    ! ----------------------------------------------------------------------------

    ! initial guess
    T1= Ttent
    ro= 0.5_wp*roc
    
    ! iterations
    do i=1,50
       T = T1
       ro= rocalc_pt_swp(p,T,ro)
       T1= Tcalc_sro_swp(s,ro,T1)
       
       err= abs(T1-T)/T
       if (err.lt.tol) then
          rocalc_ps_swp= ro
          return
       endif
    enddo

    print *,'in rocalc_ps: at iteration',i,'error=',err
    print *,'p:',p,'s:',s
    rocalc_ps_swp= 0.0_wp
    call mpistop('function rocalc_ps not converged',0)

  end function rocalc_ps_swp

  !===============================================================================
  function rocalc_pt_swp(pi,Ti,rotent)
  !===============================================================================
    !> Compute density from p and T
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: pi,Ti,rotent
    real(wp) :: rocalc_pt_swp ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: p,T,ro,ro1,tau,del,func,pcalc,err,dpdro
    ! ----------------------------------------------------------------------------

    ! initial guess
    p = pi*pci
    T = Ti*Tci
    ro= rotent*roci

    tau= 1.0_wp/T

    ! Newton's algorithm
    ! ------------------
    do i=1,20

       del= ro
       pcalc= del*zci/tau*(1.0_wp + del*d1psir_del1(tau,del))
       dpdro=     zci/tau*(1.0_wp + del*( 2.0_wp*d1psir_del1(tau,del)  &
                                           + del*d2psir_del2(tau,del)) )
       func= pcalc - p
       
       ! update solution
       ro1= ro - func/dpdro

       err= abs(ro1 - ro)/ro
       if (err.le.tol) then
          rocalc_pt_swp= ro1*roc
          return
       endif
       ro=ro1
    enddo

    print *,'in rocalc_pt: at iteration',i,'error=',err
    print *,'p:',p/pc,'T:',T/tc
    rocalc_pt_swp= 0.0_wp
    call mpistop('function rocalc_pt not converged',0)

  end function rocalc_pt_swp

  !===============================================================================
  function rocalc_st_swp(si,Ti,rotent)
  !===============================================================================
    !> Compute density from entropy s, T
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: si,Ti,rotent
    real(wp) :: rocalc_st_swp ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: s,T,ro,ro1,tau,del,scalc,err,dsdro
    ! ----------------------------------------------------------------------------

    ! initial guess
    T = Ti*Tci
    s = si/(zc*rg)
    ro= rotent*roci

    tau = 1.0_wp/T

    ! Newton's algorithm
    ! ------------------
    do i=1,20
       del= ro

       scalc= zci*(tau*(d1psir_tau1(tau,del)+d1psi0_tau1(tau)) &
            - psir(tau,del) - psi0(tau,del))

       dsdro= zci*(tau*d2psir_tau1del1(tau,del) &
            - d1psir_del1(tau,del) - d1psi0_del1(del))

       ! update solution
       ro1= ro - (scalc-s_crit-s)/dsdro

       err = abs(ro1 - ro)/ro
       if (err.le.tol) then
          rocalc_st_swp= ro1*roc
          return
       endif
       ro= ro1
    enddo

    print *,'in rocalc_st: at iteration',i,'error=',err
    print *,'s:',s,'T:',T
    rocalc_st_swp= 0.0_wp
    call mpistop('function rocalc_st not converged',0)

  end function rocalc_st_swp

  !===============================================================================
  function scalc_tro_swp(Ti,roi)
  !===============================================================================
    !> Compute entropy s from T and rho
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: Ti,roi
    real(wp) :: scalc_tro_swp ! OUTPUT
     ! ----------------------------------------------------------------------------
    real(wp) :: T,ro,tau,del
    ! ----------------------------------------------------------------------------

    T = Ti*Tci
    ro = roi*roci

    tau = 1.0_wp/T
    del = ro

    scalc_tro_swp= zci*(tau*(d1psir_tau1(tau,del) + d1psi0_tau1(tau)) &
         - psir(tau,del) - psi0(tau,del)) - s_crit

    scalc_tro_swp= scalc_tro_swp*zc*rg

  end function scalc_tro_swp

  !===============================================================================
  function tcalc_ph_swp(p,h,rotent,Ttent)
  !===============================================================================
    !> Compute temperature from p and enthalpy h
    !> - SWP EOS - 
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: p,h,Ttent,rotent
    real(wp) :: tcalc_ph_swp ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: j
    real(wp) :: ro,e,T1,err,T
    ! ----------------------------------------------------------------------------

    ! initial guess
    T1 = Ttent
    
    ! Newton's algorithm
    ! ------------------
    do j=1,100
       T = T1
       ro = rocalc_pt_swp(p,T,rotent)

       e= h-p/ro
       T1= tcalc_roero_swp(ro*e,ro,T)
       err= abs(T1-T)/T

       if (err.le.tol) then
          tcalc_ph_swp= T
          return
       endif
    enddo

  end function tcalc_ph_swp

  !===============================================================================
  function tcalc_pro_swp(pi,roi,Ttent)
  !===============================================================================
    !> Compute temperature from p and rho
    !> - SWP EOS - 
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: pi,roi,Ttent
    real(wp) :: tcalc_pro_swp ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: p,ro,tau,del,pcalc,err,T,T1,func,dpdT
    ! ----------------------------------------------------------------------------

    ! initial guess
    T = Ttent*Tci
    ro= roi*roci
    p = pi*pci

    del= ro

    ! Newton's algorithm
    ! ------------------
    do i=1,20
       tau= 1.0_wp/T

       ! function & derivative
       pcalc= del*zci/tau*(1.0_wp+del*d1psir_del1(tau,del))

       dpdT= del*zci*(1.0_wp+del*d1psir_del1(tau,del)     &
                        - del*tau*d2psir_tau1del1(tau,del) )
       func= pcalc - p
       
       ! update solution
       T1= T - func/dpdT

       err = abs(T1-T)/T
       if (err.le.tol) then
          tcalc_pro_swp= T1*Tc
          return
       endif
       T= T1
    enddo

    print *,'in tcalc_pro: at iteration',i,'error=',err
    print *,'p:',p,'ro:',ro
    tcalc_pro_swp=-1.0_wp
    call mpistop('function tcalc_pro not converged',0)

  end function tcalc_pro_swp

  !===============================================================================
  function tcalc_roero_swp(roei,roi,Ttent)
  !===============================================================================
    !> Compute temperature from rho*e and rho
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: roei,roi,Ttent
    real(wp) :: tcalc_roero_swp ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: e,ro,tau,del,err,T,T1,dedT,func,ecalc
    ! ----------------------------------------------------------------------------

    ! initial guess
    T = Ttent*Tci
    ro= roi*roci
    e = roei/roi/(zc*pcvc)

    del= ro

    ! Newton's algorithm
    ! ------------------
    do i=1,20

       tau= 1.0_wp/T

       ecalc= 100.0_wp + zci*(d1psir_tau1(tau,del)+d1psi0_tau1(tau))
       dedT = - zci*(d2psir_tau2(tau,del)+d2psi0_tau2(tau))/T**2
       func = ecalc - e
       
       ! update solution
       T1= T - func/dedT

       err= abs(T1-T)/T
       if (err.le.tol) then
          tcalc_roero_swp= T1*Tc
          return
       endif
       T = T1

    enddo

    print *,'in tcalc_roero: at iteration',i,'error=',err
    print *,'roe:',roei,'ro:',roi,'T_t :',Ttent
    tcalc_roero_swp= 0.0_wp
    call mpistop('function tcalc_roero not converged',0)
  
  end function tcalc_roero_swp

  !===============================================================================
  function tcalc_sro_swp(si,roi,Ttent)
  !===============================================================================
    !> Compute temperature from entropy s and rho
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: si,roi,Ttent
    real(wp) :: tcalc_sro_swp ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: s,ro,T,tau,tau1,del,func,der_func,scalc,err
    ! ----------------------------------------------------------------------------

    ! initial guess
    T = Ttent*Tci
    ro= roi*roci
    s = si/(zc*rg)

    tau= 1.0_wp/T
    del= ro

    ! Newton's algorithm
    ! ------------------
    do i=1,20

       ! function & derivative
       scalc=zci*(tau*(d1psir_tau1(tau,del)+d1psi0_tau1(tau)) &
            - psir(tau,del) - psi0(tau,del))
       func= scalc - s_crit - s

       der_func= zci*((d1psir_tau1(tau,del)+d1psi0_tau1(tau)) &
               + tau*(d2psir_tau2(tau,del)+d2psi0_tau2(tau)))

       ! update solution
       tau1= tau - func/der_func

       err= abs(tau1 - tau)/tau
       if (err.le.tol) then
          tcalc_sro_swp= 1.0_wp/tau
          tcalc_sro_swp= tcalc_sro_swp*Tc
          return
       endif
       tau= tau1
    enddo

    print *,'in tcalc_sro: at iteration',i,'error=',err
    print *,'s :',s,'ro:',ro
    tcalc_sro_swp= 0.0_wp
    call mpistop('function tcalc_sro not converged',0)

  end function tcalc_sro_swp

  !===============================================================================
  subroutine tcalc_Hstot_swp(T,ro,p,H_tot,s_tot,am0,coeff)
  !===============================================================================
    !> Compute temperature from total stagnation enthalpy (Newton iteration)
    !> used for imposing Riemann invariant in inflow BC
    !> - SWP EOS -
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

    ! ************************
    ! /!\ check dimensionality : work in dim or adim in loop
    ! ************************

    ! **********************
    ! v1 : not inlined
    ! **********************
    ! inlined version not written [cf mod_eos_prs.f90 for ex inlined version]

    ! Newton's algorithm
    ! ------------------
    do i=1,1000

       ! update density from total entropy
       ro=rocalc_st_swp(s_tot,T1,ro1)

       ! update pressure
       p=pcalc_tro_swp(T1,ro)

       ! update pressure derivative wrt T
       dpdT=dpdicalc_tro_swp(T1,ro)*cvcalc_tro_swp(T1,ro)

       ! update temperature

       fn=ecalc_tro_swp(T1,ro)+p/ro+0.5_wp*coeff*(am0-p)**2-H_tot

       der_fn=dedTcalc_tro_swp(T1,ro)+dpdT/ro-coeff*(am0-p)*dpdT

       T= T1-fn/der_fn

       err= abs(T1-T)/T

       if (err.le.tol*10.) then
          ro=rocalc_st_swp(s_tot,T,ro1)
          p=pcalc_tro_swp(T,ro)
          return
       endif

       ! store previous solution
       T1=T
       ro1=ro

    enddo

    print *,'in tcalc_Hstot: at iteration',i,'error=',err
    print *,'T :',T,'ro:',ro
    call mpistop('function tcalc_Hstot not converged',0)

  end subroutine tcalc_Hstot_swp

  !===============================================================================
  function vvol_swp(pi,Ti,vtent)
  !===============================================================================
    !> Compute vvol such that p(vvol)=0
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: pi,Ti,vtent
    real(wp) :: vvol_swp ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: p,T,del,tau,v0,v1,dpdv,err
    ! ----------------------------------------------------------------------------

    T= Ti
    p= pi
    
    ! initial guess (specific volume)
    v0=vtent
    tau=1.0_wp/T

    ! iterations
    do i=1,15000

       del=1.0_wp/v0

       dpdv=-del**2/tau/zc*( 1.0_wp+2.0_wp*del*d1psir_del1(tau,del) &
                                      + del**2*d2psir_del2(tau,del) )

       v1= v0 - (p-pcalc_tro_swp(Ti*Tc,roc/v0)/pc)/(-dpdv)

       err= abs(v1-v0)/v0
       if (err.le.tol) then
          vvol_swp= v1
          return
       endif

       v0=v1
    enddo

    print *,'in vvol: at iteration',i,'error=',err
    call mpistop('function vvol not converged',0)

  end function vvol_swp

  !===============================================================================
  function vvol_d1_swp(Ti,vtent)
  !===============================================================================
    !> Compute vvol_d1 from Ti and tentative v
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: Ti,vtent
    real(wp) :: vvol_d1_swp ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: T,del,tau,err
    real(wp) :: v0,v1,dpdv,d2pdv2,dpdro,d2pdro2
    ! ----------------------------------------------------------------------------

    T= Ti
    tau=1.0_wp/T
    
    ! initial guess (specific volume)
    v0=vtent

    ! iterations
    aloop: do i=1,15000

       del=1.0_wp/v0

       dpdro= 1.0_wp/tau/zc*( 1.0_wp+2.0_wp*del*d1psir_del1(tau,del) &
                                       + del**2*d2psir_del2(tau,del) )
       d2pdro2= 1.0_wp/tau/zc*( 2.0_wp*d1psir_del1(tau,del) &
                          + 4.0_wp*del*d2psir_del2(tau,del) &
                              + del**2*d3psir_del3(tau,del) )

       dpdv=-del**2*dpdro
       d2pdv2=2.0_wp*del**3*dpdro + del**4*d2pdro2

       v1= v0 - dpdv/d2pdv2

       err = abs(v1-v0)/v0
       if (err.le.tol) then
          vvol_d1_swp= v1
          return
       endif

       !----- Returns negative value if not converged -----
       if (v1.lt.0_wp) then
          v0= -1.0_wp
          exit aloop
       elseif (v1.gt.1.e7_wp) then
          v0= -1.0_wp
          exit aloop
       endif

       v0=v1
    enddo aloop

    print *,'in vvol_d1: at iteration',i,'error=',err
    call mpistop('function vvol_d1 not converged',0)

  end function vvol_d1_swp

  !===============================================================================
  function vvol_d2_swp(Ti,vtent)
  !===============================================================================
    !> Compute vvol_d2 from Ti and tentative v
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: Ti,vtent
    real(wp) :: vvol_d2_swp ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: T,del,tau,err
    real(wp) :: v0,v1,d2pdv2,d3pdv3,dpdro,d2pdro2,d3pdro3
    ! ----------------------------------------------------------------------------

    T= Ti
    tau=1.0_wp/T
   
    ! initial guess (specific volume)
    v0=vtent

    ! iterations
    do i=1,15000

       del=1.0_wp/v0

       dpdro= 1.0_wp/tau/zc*( 1.0_wp+2.0_wp*del*d1psir_del1(tau,del) &
                                       + del**2*d2psir_del2(tau,del) )

       d2pdro2= 1.0_wp/tau/zc*( 2.0_wp*d1psir_del1(tau,del) &
                          + 4.0_wp*del*d2psir_del2(tau,del) &
                              + del**2*d3psir_del3(tau,del) )

       d3pdro3= 1.0_wp/tau/zc*( 6.0_wp*d2psir_del2(tau,del) &
                          + 6.0_wp*del*d3psir_del3(tau,del) &
                              + del**2*d4psir_del4(tau,del) )

       d2pdv2= 2.0_wp*del**3*dpdro+del**4*d2pdro2
       d3pdv3=-6.0_wp*del**4*dpdro-6.0_wp*del**5*d2pdro2-del**6*d3pdro3

       v1= v0 - d2pdv2/d3pdv3

       err = abs(v1-v0)/v0
       if (err.le.tol) then
          vvol_d2_swp= v1
          return
       endif

       v0=v1
    enddo

    print *,'in vvol_d2: at iteration',i,'error=',err
    call mpistop('function vvol_d2 not converged',0)
    
  end function vvol_d2_swp

  !===============================================================================
  function intpcalc_tro_swp(Ti,roi)
  !===============================================================================
    !> Compute pressure integral from T and rho
    !> - SWP EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: Ti,roi
    real(wp) :: intpcalc_tro_swp ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: T,ro,tau,del
    ! ----------------------------------------------------------------------------

    T = Ti
    ro= roi

    tau= 1.0_wp/T
    del= ro
    
    intpcalc_tro_swp=-1.0_wp/(zc*tau)*(log(del)+psir(tau,del))

  end function intpcalc_tro_swp

end module mod_eos_swp
