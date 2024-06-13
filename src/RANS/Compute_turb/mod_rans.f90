!===============================================================================
module mod_rans
!===============================================================================
  !> Module for RANS modeling
!===============================================================================
  use precision
  use mod_constant
  implicit none
  ! -----------------------------------------------------------------------------
  ! Flow stuff
  ! ==========
  ! Arrays
  real(wp), dimension(:,:,:), allocatable :: Om
  
  ! Components 
  ! real (wp) :: om12,om13,om23,magn_om,tauR11,tauR12,tauR13,tauR22,tauR23,tauR33 ! Test perfo
  real (wp) :: om12,om13,om23,magn_om ! Test perfo
  ! ------------------------------------------------------------------------------
  ! Common to every model
  ! =====================
  real(wp) :: prrt ! turbulent Prandtl and heat transfer coefficient

  ! Bradshaw constant
  real(wp) :: a1

  ! El famoso mut
  real(wp), dimension(:,:,:), allocatable :: mut,nut
  ! ------------------------------------------------------------------------------
  ! Spalart-Allmaras
  ! ================

  ! Reference nutilde
  real(wp) :: Nutref
  
  ! Source terms
  real(wp), dimension(:,:,:), allocatable :: Sterm
  real(wp) :: St1,St2,St3,St4

  ! Model constants & auxilary relations
  real(wp) :: cb1,cb2,ft2,cnu1,cnu2,cnu3,sig,cw1,cw2,cw3,kap,fw,fnu1,fnu2,khi,g_sa,r,Stil,S
  ! real(wp) :: prrt,cokt ! turbulent Prandtl and heat transfer coefficient ! Test perfo
  
  ! Additional shortcut expressions
  real(wp) :: dnutili2 ! (dnutil/dxi)*(dnutil/dxi)

  !
  real(wp), dimension(:,:,:), allocatable :: lengthscale,delta_max

  ! Transitional models:
  ! ================
  ! Inlet turbulent itnensity:
    real(wp)  :: tu_inlet
  ! --------------------------------------------------------------------------------------------------------------------
  ! Crivellini algebraic transitinoal model (Computers and Fluids 2023: https://doi.org/10.1016/j.compfluid.2023.105791)
  ! =======================================
    real(wp) :: gamma_cb, term1, term2, cw4,cw5, cw2_lre, ksi1, ksi2, re_vort, re_theta, re_crit

  ! ------------------------------------------------------------------------------
  ! =====================
  !     SA-Gamma-Re_theta
  ! =====================
  
  ! Source terms of gamma and Re_theta:
  real(wp)  :: magn_u, du_ds, magn_u2, T_stheta, delta_stheta, F_theta
  real(wp)  :: theta, lambda_theta, f_lambda_t, Re_theta_t, re_theta_crit
  real(wp)  :: re_nu, re_turb, F_turb, F_onset1, F_onset2, F_onset3, F_onset, F_length, F_reattach
  real(wp)  :: p_gamma, d_gamma
  ! Intermittency in equation of nutil:
  real(wp)  :: intermittency_sep
  real(wp), dimension(:,:,:), allocatable :: Sgamma, Stheta
  ! Model parameters:
  real(wp), parameter  :: sigma_theta_t= 2.0_wp, ce2=50.0_wp, ctheta_t=0.03_wp, ca1=2.0_wp, ca2=0.06_wp


end module mod_rans
