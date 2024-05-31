!=================================================================================
module mod_tranprop
!=================================================================================
  !> author: XG
  !> date: January 2024
  !> Module to define transport properties
!=================================================================================
  use mod_fluid ! for: eos_type
  use warnstop  ! for: mpistop
  implicit none
  ! ------------------------------------------------------------------------------
  ! Interfaces for procedure pointer
  ! ------------------------------------------------------------------------------
  abstract interface
     real(wp) function thermo_type1(T,ro)
       use precision
       implicit none
       real(wp), intent(in) :: T,ro
     end function thermo_type1
  end interface
  ! ------------------------------------------------------------------------------
  procedure(thermo_type1), pointer :: viscosity_law => null()
  ! ------------------------------------------------------------------------------
  ! Interfaces for procedure pointer
  ! ------------------------------------------------------------------------------
  abstract interface
     real(wp) function thermo_type2(mu,T,ro)
       use precision
       implicit none
       real(wp), intent(in) :: mu,T,ro
     end function thermo_type2
  end interface
  ! ------------------------------------------------------------------------------
  procedure(thermo_type2), pointer :: thconductivity => null()
  ! ------------------------------------------------------------------------------
  ! Perfect gas variables
  ! =====================
  real(wp) :: T0,S0,mu0,var1,var2,var3,c0
  ! ------------------------------------------------------------------------------
  ! Dense gas variables
  ! ===================
  ! -> parameters for viscosity computations
  ! ========================================
  real(wp) :: vmolc
  real(wp) :: k_chung
  real(wp), parameter :: avisc = 1.16145e+0_wp
  real(wp), parameter :: bvisc = 1.48740e-1_wp
  real(wp), parameter :: cvisc = 5.24870e-1_wp
  real(wp), parameter :: dvisc = 7.73200e-1_wp
  real(wp), parameter :: evisc = 2.16178e+0_wp
  real(wp), parameter :: fvisc = 2.43787e+0_wp
  real(wp), parameter :: gvisc =-6.43500e-4_wp
  real(wp), parameter :: hvisc = 7.27371e+0_wp
  real(wp), parameter :: svisc = 1.80323e+1_wp
  real(wp), parameter :: wvisc =-7.68300e-1_wp

  real(wp), dimension(10), private :: av = (/ 6.32402e+0_wp,  1.21020e-3_wp, &
                                              5.28346e+0_wp,  6.62263e+0_wp, &
                                              1.97454e+1_wp, -1.89992e+0_wp, &
                                              2.42745e+1_wp,  7.97160e-1_wp, &
                                             -2.38160e-1_wp,  6.86290e-2_wp /)
  real(wp), dimension(10), private :: bv = (/ 5.04119e+1_wp, -1.15360e-3_wp, &
                                              2.54209e+2_wp,  3.80957e+1_wp, &
                                              7.63034e+0_wp, -1.25367e+1_wp, &
                                              3.44945e+0_wp,  1.11764e+0_wp, &
                                              6.76950e-2_wp,  3.47930e+0_wp /)
  real(wp), dimension(10), private :: cv = (/-5.16801e+1_wp, -6.25710e-3_wp, &
                                             -1.68481e+2_wp, -8.46414e+0_wp, &
                                             -1.43544e+1_wp,  4.98529e+0_wp, &
                                             -1.12913e+1_wp,  1.23480e-2_wp, &
                                             -8.16300e-1_wp,  5.92560e-1_wp /)
  real(wp), dimension(10), private :: dv = (/ 1.18902e+3_wp,  3.72830e-2_wp, &
                                              3.89827e+3_wp,  3.14178e+1_wp, &
                                              3.15267e+1_wp, -1.81507e+1_wp, &
                                              6.93466e+1_wp, -4.11661e+0_wp, &
                                              4.02528e+0_wp, -7.26630e-1_wp  /)
  real(wp), dimension(10) :: AA1
 
  ! -> law of Wen, Meng, Huber & Wu (J. Chem. Eng. Data, 2017)
  ! ==========================================================
  real(wp), dimension(9), private :: bw = (/ -19.572881_wp, 219.73999_wp, &
                                             -1015.3226_wp, 2471.0125_wp, &
                                             -3375.1717_wp, 2491.6597_wp, &
                                             -787.26086_wp, 14.085455_wp, &
                                             -0.34664158_wp /)
  real(wp), dimension(7), private :: cw = (/ 22.0057_wp,  231.063_wp, &
                                            0.423359_wp,-0.122057_wp, &
                                             18.4610_wp, -11.1393_wp, &
                                             1.67777_wp /)
  
  ! ------------------------------------------------------------------------------
  ! -> parameters for thermal conductivity computations
  ! ===================================================
  real(wp), dimension(7), private :: atc = (/ 2.41657e+0_wp, -5.09240e-1_wp, &
                                              6.61069e+0_wp,  1.45425e+1_wp, &
                                              7.92740e-1_wp, -5.86340e+0_wp, &
                                              8.11710e+1_wp /)
  real(wp), dimension(7), private :: btc = (/ 7.48240e-1_wp, -1.50936e+0_wp, &
                                              5.62073e+0_wp, -8.91387e+0_wp, &
                                              8.20190e-1_wp,  1.28005e+1_wp, &
                                              1.14158e+2_wp /)
  real(wp), dimension(7), private :: ctc = (/-9.18580e-1_wp, -4.99912e+1_wp, &
                                              6.47599e+1_wp, -5.63794e+0_wp, &
                                             -6.93690e-1_wp,  9.58926e+0_wp, &
                                             -6.08410e+1_wp /)
  real(wp), dimension(7), private :: dtc = (/ 1.21721e+2_wp,  6.99834e+1_wp, &
                                              2.70389e+1_wp,  7.43435e+1_wp, &
                                              6.31734e+0_wp,  6.55292e+1_wp, &
                                              4.66775e+2_wp /)
  real(wp), dimension(7) :: BB1
  ! ------------------------------------------------------------------------------

  interface
     !===============================================================================
     ! Initialization
     !===============================================================================
     module subroutine init_viscosity
     end subroutine init_viscosity
     !===============================================================================
     ! Dynamic viscosity laws
     !===============================================================================
     module function viscosity_law_sutherland(T,ro)
       ! ----------------------------------------------------------------------------
       real(wp), intent(in) :: T,ro
       real(wp) :: viscosity_law_sutherland
       ! ----------------------------------------------------------------------------
     end function viscosity_law_sutherland
     !===============================================================================
     module function viscosity_law_power(T,ro)
       ! ----------------------------------------------------------------------------
       real(wp), intent(in) :: T,ro
       real(wp) :: viscosity_law_power
       ! ----------------------------------------------------------------------------
     end function viscosity_law_power
     !===============================================================================
     module function viscosity_law_refprop(T,ro)
       ! ----------------------------------------------------------------------------
       real(wp), intent(in) :: T,ro
       real(wp) :: viscosity_law_refprop
       ! ----------------------------------------------------------------------------
     end function viscosity_law_refprop
     !===============================================================================
     module function viscosity_law_chung(T,ro)
       ! ----------------------------------------------------------------------------
       real(wp), intent(in) :: T,ro
       real(wp) :: viscosity_law_chung
       ! ----------------------------------------------------------------------------
     end function viscosity_law_chung
     !===============================================================================
     module function viscosity_law_wen(T,ro)
       ! ----------------------------------------------------------------------------
       real(wp), intent(in) :: T,ro
       real(wp) :: viscosity_law_wen
       ! ----------------------------------------------------------------------------
     end function viscosity_law_wen
     !===============================================================================
     ! Thermal conductivity laws
     !===============================================================================
     module function thconductivity_Prcost(mu,T,ro)
       ! ----------------------------------------------------------------------------
       real(wp), intent(in) :: mu,T,ro
       real(wp) :: thconductivity_Prcost
       ! ----------------------------------------------------------------------------
     end function thconductivity_Prcost
     !===============================================================================
     module function thconductivity_refprop(mu,T,ro)
       ! ----------------------------------------------------------------------------
       real(wp), intent(in) :: mu,T,ro
       real(wp) :: thconductivity_refprop
       ! ----------------------------------------------------------------------------
     end function thconductivity_refprop
     !===============================================================================
     module function thconductivity_chung(mu,T,ro)
       ! ----------------------------------------------------------------------------
       real(wp), intent(in) :: mu,T,ro
       real(wp) :: thconductivity_chung
       ! ----------------------------------------------------------------------------
     end function thconductivity_chung
     !===============================================================================
  end interface

end module mod_tranprop
