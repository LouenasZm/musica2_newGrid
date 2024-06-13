!==============================================================================================
module mod_flow
!==============================================================================================
  !> Module for flow arrays
!==============================================================================================
  use mod_grid
  use mod_bc
  implicit none
  ! -------------------------------------------------------------------------------------------
  ! Flow variables
  ! -------------------------------------------------------------------------------------------
  ! conservative variables
  real(wp), dimension(:,:,:), allocatable, target :: rho,rhou,rhov,rhow,rhoe
  ! updated conservative variables
  real(wp), dimension(:,:,:), allocatable :: rho_n,rhou_n,rhov_n,rhow_n,rhoe_n
  real(wp), dimension(:,:,:), allocatable :: rhon,rhoun,rhovn,rhown,rhoen
  ! primitive variables
  real(wp), dimension(:,:,:), allocatable, target :: uu,vv,ww,prs,Tmp
  ! velocity derivatives
  real(wp), dimension(:,:,:), allocatable :: dux,dvx,dwx,duy,dvy,dwy,duz,dvz,dwz
  ! temperature derivatives
  real(wp), dimension(:,:,:), allocatable :: dTx,dTy,dTz
  ! sound speed
  real(wp), dimension(:,:,:), allocatable, target :: c_
  ! viscosity and thermal conductivity
  real(wp), dimension(:,:,:), allocatable :: visc,cok
  ! -------------------------------------------------------------------------------------------
  ! Numerical Fluxes
  ! -------------------------------------------------------------------------------------------
  ! direction 1 (x or ksi)
  real(wp), dimension(:,:,:), allocatable :: Frho,Frhou,Frhov,Frhow,Frhoe
  ! direction 2 (y or eta)
  real(wp), dimension(:,:,:), allocatable :: Grho,Grhou,Grhov,Grhow,Grhoe
  ! direction 3 (z or zet)
  real(wp), dimension(:,:,:), allocatable :: Hrho,Hrhou,Hrhov,Hrhow,Hrhoe
  ! derivatives of fluxes: solution increments
  real(wp), dimension(:,:,:), allocatable :: Krho,Krhou,Krhov,Krhow,Krhoe
  real(wp), dimension(:,:,:), allocatable :: Krhobis,Krhoubis,Krhovbis,Krhowbis,Krhoebis
  ! strain tensor
  real(wp), dimension(:,:,:), allocatable :: S11,S12,S13,S22,S23,S33
  real(wp), dimension(:,:,:), allocatable :: S11f,S12f,S13f,S22f,S23f,S33f
  ! ===========================================================================================
  ! RANS arrays (Spalart-Allmaras)
  ! -------------------------------------------------------------------------------------------
  ! updated variable
  real(wp), dimension(:,:,:), allocatable, target :: nutil_n,nutil
  ! nutilde derivatives
  real(wp), dimension(:,:,:), allocatable :: dnutilx,dnutily,dnutilz
  ! -------------------------------------------------------------------------------------------
  ! Numerical Fluxes
  ! -------------------------------------------------------------------------------------------
  ! directions 1,2,3 (x/y/z or ksi/eta/z), and derivatives of fluxes: solution increments
  real(wp), dimension(:,:,:), allocatable :: Fnutil,Gnutil,Hnutil,Knutil
  ! -------------------------------------------------------------------------------------------
  ! Statistics and mean quantities
  ! -------------------------------------------------------------------------------------------
  ! space- and time-averaged quantities [plane output]
  real(wp), dimension(:,:,:), allocatable, target :: avg_s,avg_t
  ! global time-averaged quantities [plane output]
  real(wp), dimension(:,:,:), allocatable, target :: avg_tg
  ! time-averaged quantities [volume output]
  real(wp), dimension(:,:,:,:), allocatable, target :: avg_v
  ! wall time-averaged quantities [plane output]
  real(wp), dimension(:,:,:), allocatable, target :: avg_w
  ! arrays for particular outputs: uvar (user variable)
  real(wp), dimension(:,:,:,:), allocatable :: uvar
  ! ===================
  ! SA-Gamma-Re_theta
  ! -------------------------------------------------------------------------------------------
   ! updated variable
  real(wp), dimension(:,:,:), allocatable, target :: rhore_theta, rhore_theta_n, rho_gamma, rhogamma_n ! rho_gama is conservative form of intermittency
  real(wp), dimension(:,:,:), allocatable, target :: reynolds_theta, intermittency, intermittency_eff
  ! derivatives
  real(wp), dimension(:,:,:), allocatable :: dre_thetax,dre_thetay,dre_thetaz, dgammax, dgammay, dgammaz
  real(wp), dimension(:,:,:), allocatable :: drhore_thetax, drhore_thetay, drhore_thetaz
  real(wp), dimension(:,:,:), allocatable :: drho_gammax, drho_gammay, drho_gammaz
  ! -------------------------------------------------------------------------------------------
  ! Numerical Fluxes
  ! -------------------------------------------------------------------------------------------
  ! directions 1,2,3 (x/y/z or ksi/eta/z), and derivatives of fluxes: solution increments
  real(wp), dimension(:,:,:), allocatable :: Fre_theta,Gre_theta, Hre_theta, Kre_theta
  real(wp), dimension(:,:,:), allocatable :: Fgamma, Ggamma, Hgamma,Kgamma
  !dir$ attributes align:32 :: rho
  !dir$ attributes align:32 :: rhou
  !dir$ attributes align:32 :: rhov
  !dir$ attributes align:32 :: rhow
  !dir$ attributes align:32 :: rhoe
  !dir$ attributes align:32 :: rhore_theta
  !dir$ attributes align:32 :: rho_gamma
  
  
  !dir$ attributes align:32 :: uu
  !dir$ attributes align:32 :: vv
  !dir$ attributes align:32 :: ww
  !dir$ attributes align:32 :: prs
  !dir$ attributes align:32 :: Tmp
  !dir$ attributes align:32 :: c_
  !dir$ attributes align:32 :: visc
  !dir$ attributes align:32 :: cok
  !dir$ attributes align:32 :: reynolds_theta
  !dir$ attributes align:32 :: intermittency
  !dir$ attributes align:32 :: intermittency_eff

  
  !dir$ attributes align:32 :: dux
  !dir$ attributes align:32 :: dvx
  !dir$ attributes align:32 :: dwx
  !dir$ attributes align:32 :: dTx
  !dir$ attributes align:32 :: dre_thetax
  !dir$ attributes align:32 :: dgammax
  !dir$ attributes align:32 :: drhore_thetax
  !dir$ attributes align:32 :: drho_gammax

  !dir$ attributes align:32 :: duy
  !dir$ attributes align:32 :: dvy
  !dir$ attributes align:32 :: dwy
  !dir$ attributes align:32 :: dTy
  !dir$ attributes align:32 :: dre_thetay
  !dir$ attributes align:32 :: dgammay
  !dir$ attributes align:32 :: drhore_thetay
  !dir$ attributes align:32 :: drho_gammay

  !dir$ attributes align:32 :: duz
  !dir$ attributes align:32 :: dvz
  !dir$ attributes align:32 :: dwz
  !dir$ attributes align:32 :: dTz
  !dir$ attributes align:32 :: dre_thetaz
  !dir$ attributes align:32 :: dgammaz
  !dir$ attributes align:32 :: drhore_thetaz
  !dir$ attributes align:32 :: drho_gammaz

  !dir$ attributes align:32 :: rho_n
  !dir$ attributes align:32 :: rhou_n
  !dir$ attributes align:32 :: rhov_n
  !dir$ attributes align:32 :: rhow_n
  !dir$ attributes align:32 :: rhoe_n
  !dir$ attributes align:32 :: rhore_theta_n
  !dir$ attributes align:32 :: rhogamma_n

  !dir$ attributes align:32 :: Krho
  !dir$ attributes align:32 :: Krhou
  !dir$ attributes align:32 :: Krhov
  !dir$ attributes align:32 :: Krhow
  !dir$ attributes align:32 :: Krhoe
  !dir$ attributes align:32 :: Kre_theta
  !dir$ attributes align:32 :: Kgamma


  !dir$ attributes align:32 :: Frho
  !dir$ attributes align:32 :: Frhou
  !dir$ attributes align:32 :: Frhov
  !dir$ attributes align:32 :: Frhow
  !dir$ attributes align:32 :: Frhoe
  !dir$ attributes align:32 :: Fre_theta
  !dir$ attributes align:32 :: Fgamma

  !dir$ attributes align:32 :: Grho
  !dir$ attributes align:32 :: Grhou
  !dir$ attributes align:32 :: Grhov
  !dir$ attributes align:32 :: Grhow
  !dir$ attributes align:32 :: Grhoe
  !dir$ attributes align:32 :: Gre_theta
  !dir$ attributes align:32 :: Ggamma

  !dir$ attributes align:32 :: Hrho
  !dir$ attributes align:32 :: Hrhou
  !dir$ attributes align:32 :: Hrhov
  !dir$ attributes align:32 :: Hrhow
  !dir$ attributes align:32 :: Hrhoe
  !dir$ attributes align:32 :: Hre_theta
  !dir$ attributes align:32 :: Hgamma

end module mod_flow
