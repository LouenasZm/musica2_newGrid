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
  ! Spalart-Allmaras
  ! ================
  ! El famoso mut
  real(wp), dimension(:,:,:), allocatable :: mut,nut

  ! Reference nutilde
  real(wp) :: Nutref
  
  ! Source terms
  real(wp), dimension(:,:,:), allocatable :: Sterm
  real(wp) :: St1,St2,St3,St4

  ! Model constants & auxilary relations
  real(wp) :: cb1,cb2,ft2,cnu1,cnu2,cnu3,sig,cw1,cw2,cw3,kap,fw,fnu1,fnu2,khi,g_sa,r,Stil,S
  ! real(wp) :: prrt,cokt ! turbulent Prandtl and heat transfer coefficient ! Test perfo
  real(wp) :: prrt ! turbulent Prandtl and heat transfer coefficient

  ! Bradshaw constant
  real(wp) :: a1
  
  ! Additional shortcut expressions
  real(wp) :: dnutili2 ! (dnutil/dxi)*(dnutil/dxi)

  !
  real(wp), dimension(:,:,:), allocatable :: lengthscale,delta_max


end module mod_rans
