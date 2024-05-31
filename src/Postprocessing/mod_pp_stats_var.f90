!==============================================================================
module mod_pp_stats_var
!==============================================================================
  !> Module for post-processing variables
!==============================================================================
  use precision
  implicit none
  ! ---------------------------------------------------------------------------
  !integer :: type_data
  !character(len=60) :: dirDATA,dirRESU
  !character(len=30) :: name_output
  !real(wp) :: time_ini ! initial time for samples [useful if is_timestamp=T]
  
  ! wall quantities
  ! ---------------
  real(wp), dimension(:), allocatable :: rhowall,muwall,nuwall,pwall,cwall,ewall,hwall
  real(wp), dimension(:), allocatable :: Twall,dTwall,lawall,Qwall,tauwall,dudywall
  
  ! bulk quantities (channel flow)
  ! ---------------
  real(wp) :: rhobulk,ubulk,Rebulkb

  ! friction quantities
  ! -------------------
  real(wp), dimension(:), allocatable :: c_f,u_tau
  
  ! thicknesses
  ! -----------
  integer, dimension(:), allocatable :: j99,j_inn,j_log
  real(wp), dimension(:), allocatable :: delt99,deltas,deltatheta,Hfac
  real(wp), dimension(:), allocatable :: delt99_,deltas_,deltatheta_,Hfac_
  real(wp), dimension(:), allocatable :: deltas_c,deltatheta_c,Hfac_c ! compressible version
  
  ! Reynolds numbers
  ! ----------------
  real(wp), dimension(:), allocatable :: Re_d99,Re_t,Re_ds,Re_theta,Re_d2,Re_ds_c,Re_theta_c
  
  ! Variables for mean quantitities
  ! ===============================
  ! Main variables
  ! --------------
  ! primitives variables 
  real(wp), dimension(:,:), allocatable :: rhom,um,vm,wm,pm,Tm
  ! TKE (turbulent kinetic energy))
  real(wp), dimension(:,:), allocatable :: ktm
  ! velocity gradients
  real(wp), dimension(:,:), allocatable :: duxm,duym,duzm ! <- not in stats
  real(wp), dimension(:,:), allocatable :: dvxm,dvym,dvzm ! <- not in stats
  real(wp), dimension(:,:), allocatable :: dwxm,dwym,dwzm ! <- not in stats
  ! vorticity
  real(wp), dimension(:,:), allocatable :: vrtm,vrtxm,vrtym,vrtzm ! <- vrtm not in stats
  
  ! Second-order quantities (to compute rms)
  ! -----------------------
  ! second-order turbulent intensities
  real(wp), dimension(:,:), allocatable :: u2m,v2m,w2m,uvm,uwm,vwm
  ! second-order thermodynamic variables
  real(wp), dimension(:,:), allocatable :: rho2m,p2m,T2m
  ! second-order velocity gradients
  real(wp), dimension(:,:), allocatable :: dux2m,duy2m,duz2m
  real(wp), dimension(:,:), allocatable :: dvx2m,dvy2m,dvz2m
  real(wp), dimension(:,:), allocatable :: dwx2m,dwy2m,dwz2m
  ! second-order vorticity
  real(wp), dimension(:,:), allocatable :: vrt2m,vrtx2m,vrty2m,vrtz2m ! <- vrt2m not in stats

  ! Skewness and flatness
  ! ---------------------
  real(wp), dimension(:,:), allocatable :: u3m,u4m ! u-velocity
  real(wp), dimension(:,:), allocatable :: p3m,p4m ! pressure

  ! Viscous terms
  ! -------------
  ! strain rate components
  real(wp), dimension(:,:), allocatable :: tau11m,tau12m,tau13m,tau22m,tau23m,tau33m
  ! tau_ij times velocity components
  real(wp), dimension(:,:), allocatable :: tau11um,tau12um,tau12vm,tau22um,tau22vm,tau13wm,tau23wm
  ! tau_ij times velocity gradients
  real(wp), dimension(:,:), allocatable :: tau11duxm,tau11dvxm,tau12duxm
  real(wp), dimension(:,:), allocatable :: tau12duym,tau12dvxm,tau12dvym,tau13duzm
  real(wp), dimension(:,:), allocatable :: tau13dvzm,tau13dwxm,tau22duym,tau22dvym
  real(wp), dimension(:,:), allocatable :: tau23duzm,tau23dvzm,tau23dwym,tau33dwzm
  
  ! Velocity divergence
  ! -------------------
  real(wp), dimension(:,:), allocatable :: divm,div2m ! first/second-order
  real(wp), dimension(:,:), allocatable :: rhodivm,rhodiv2m ! Favre 

  ! Heat exchange (Nusselt)
  ! -------------
  real(wp), dimension(:,:), allocatable :: ladTxm,ladTym,ladTzm 
  
  ! Terms specific to TKE budgets
  ! -----------------------------
  real(wp), dimension(:,:), allocatable :: rhoduydvxm,rhoduzdwxm,rhodvzdwym
  real(wp), dimension(:,:), allocatable :: rhouum,rhovvm,rhowwm,rhouvm,rhouwm,rhovwm
  real(wp), dimension(:,:), allocatable :: rhoTTm,rhovTm
  ! triple
  real(wp), dimension(:,:), allocatable :: rhouuum,rhovvvm,rhowwwm 
  real(wp), dimension(:,:), allocatable :: rhouuvm,rhouvvm,rhowwvm,rhowwum

  ! Favre averaging (compressibility effects)
  ! ---------------
  ! velocity
  real(wp), dimension(:,:), allocatable :: rhoum,rhovm,rhowm
  ! energy eq.
  real(wp), dimension(:,:), allocatable :: rhoTm,rhoem,rhohm
  ! velocity gradients
  real(wp), dimension(:,:), allocatable :: rhoduxm,rhoduym,rhoduzm
  real(wp), dimension(:,:), allocatable :: rhodvxm,rhodvym,rhodvzm
  real(wp), dimension(:,:), allocatable :: rhodwxm,rhodwym,rhodwzm
  ! second-order velocity gradients
  real(wp), dimension(:,:), allocatable :: rhodux2m,rhoduy2m,rhoduz2m
  real(wp), dimension(:,:), allocatable :: rhodvx2m,rhodvy2m,rhodvz2m
  real(wp), dimension(:,:), allocatable :: rhodwx2m,rhodwy2m,rhodwz2m  
  ! vorticity
  real(wp), dimension(:,:), allocatable :: rhovrtm,rhovrtxm,rhovrtym,rhovrtzm
  ! second-order vorticity
  real(wp), dimension(:,:), allocatable :: rhovrt2m,rhovrtx2m,rhovrty2m,rhovrtz2m  

  ! Additional thermodynamic variables (mean and second-order)
  ! ----------------------------------
  ! dynamic viscosity & thermal conductivity
  real(wp), dimension(:,:), allocatable :: mum,lam,mu2m,la2m
  ! sound speed and Mach number
  real(wp), dimension(:,:), allocatable :: cm,Mm,ccm,M2m
  ! internal energy, enthalpy & entropy
  real(wp), dimension(:,:), allocatable :: em,hm,sm,e2m,h2m,s2m

  ! Non-ideal thermodynamic variables
  ! ---------------------------------
  ! fundamental derivative of gas dynamics Gamma
  real(wp), dimension(:,:), allocatable :: Gm,G2m
  ! heat capacities
  real(wp), dimension(:,:), allocatable :: cpm,cvm,cv2m,cp2m
  ! Prandtl and Eckert numbers
  real(wp), dimension(:,:), allocatable :: prm,eckm,pr2m,eck2m
  
  ! Correlations between two variables
  ! ----------------------------------
  ! pressure-divergence
  real(wp), dimension(:,:), allocatable :: pdivm
  ! 
  real(wp), dimension(:,:), allocatable :: upm,vpm,uTm,vTm,usm,vsm,rhopm,pTm,spm,sTm,rhosm
  ! with pressure ------> budgets
  real(wp), dimension(:,:), allocatable :: pduxm,pdvym,pdwzm,pduym,pdvxm
  ! with enthalpy ------> budgets
  real(wp), dimension(:,:), allocatable :: hum,hvm,hwm,rhohum,rhohvm,rhohwm
  ! with derivative of gas dynamics Gamma (non-ideal)
  real(wp), dimension(:,:), allocatable :: Grhom,Gpm,Gsm,GTm,Gum,Gvm
  
  ! Numerical dissipation
  ! ---------------------
  real(wp), dimension(:,:), allocatable :: Drhom,Drhoum,Drhovm,Drhowm,Drhoem
  real(wp), dimension(:,:), allocatable :: uDrom,vDrom,wDrom,uDroum,vDrovm,wDrowm

  ! Channel flow
  ! ============
  ! total temperature
  real(wp), dimension(:,:), allocatable :: Ttotm,Ttot2m
  ! Ducros sensor
  real(wp), dimension(:,:), allocatable :: ducrm,ducr2m
  ! bulk forcing term
  real(wp), dimension(:,:), allocatable :: rhofm
  
  ! Others
  ! ---------------------
  real(wp), dimension(:,:), allocatable :: rhoukuk,kolmog,l_tayl,L_f_scale,dissip,Re_fst

  ! Favre-averaged quantitities
  ! ===========================
  ! mean variables
  real(wp), dimension(:,:), allocatable :: Tmfavre,umfavre,vmfavre,wmfavre  
  ! fluctuating variables
  real(wp), dimension(:,:), allocatable :: Tffavre,uffavre,vffavre,wffavre
  
  ! non-dimensional coordinates
  ! ===========================
  real(wp), dimension(:,:), allocatable :: ygw,ygp,ygs,xgp,zgp,xgs,zgs
  ! ---------------------------------------------------------------------------

end module mod_pp_stats_var
