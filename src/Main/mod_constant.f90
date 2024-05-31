!=============================================================================
module mod_constant
!=============================================================================
  !> Module for Constants (TO BE CHANGED)
!=============================================================================
  use precision
  implicit none
  !---------------------------------------------------------------------------
  real(wp), parameter :: pi=3.141592653589793_wp
  real(wp), parameter :: twopi=2.0_wp*pi
  real(wp), parameter :: nep=2.718281828459045_wp
  real(wp), parameter :: cutoff=1.e-16_wp
  real(wp), parameter :: tiny=epsilon(1.0_wp)
  real(wp), parameter :: ONE_THIRD=1.0_wp/3.0_wp
  !---------------------------------------------------------------------------

  integer, parameter  :: uni=69

  !---------------------------------------------------------------------------
  ! run mode : type of simulation
  integer, parameter :: PRE_PROCESSING = 0
  integer, parameter :: FROM_SCRATCH   = 1
  integer, parameter :: FROM_FILE      = 2
  integer, parameter :: FROM_INTERP    = 3
  integer, parameter :: POST_PROCESSING= 4  
  !---------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  ! Simulation parameters
  logical  ::  TGV  = .false.   &
             , CHIT = .false.   &
             , CHAN = .false.   &
             , CAV  = .false.   &
             , ACT  = .false.   &
             , PHILL= .false.   &
             , SRC  = .false.   &
             , STBL = .false.   &
             , CYL  = .false.   &
             , SHIT = .false.   &
             , LE = .false.   &
             , TE = .false.   &
             , TURB = .false.
  !---------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  ! Directories
  character(len=200) :: dirDATA,dirRESU,dirGRID
  character(len=200) :: nameGRID

  ! Reference thermodynamic quantities
  real(wp) :: rho_ref,T_ref,p_ref,c_ref
  real(wp) :: e_ref ! <- USED only in setupref, residuals, init_chit
  real(wp) :: g_ref ! <- USED only in setupref
  real(wp) :: s_ref,h_ref ! <- USED only in mod_blasius
  real(wp) :: mu_ref,muw_ref
  real(wp) :: la_ref ! <- NOT USED
  ! Reference flow quantities
  real(wp) :: Re_ref,Mach ! <- given in param
  real(wp) :: u_ref,L_ref ! <- a priori deduced from Re_ref,Mach
  ! Reference inlet profile
  real(wp), dimension(:), allocatable :: Tin_ref ! kept for temporary compatibility with mean0 TO BE CHANGED

  ! components of constant velocity vector (normalized by reference velocity)
  real(wp) :: Uc1,Uc2,Uc3 ! <- given in param

  ! imposed directions of velocity at inlet
  real(wp) :: theta_ref,phi_ref
  ! imposed exit pressure
  real(wp) :: p_exit
  !------------------------------------------------------------------------------
  ! Type of display at screen: if F, strict necessary / if T, display everything
  logical :: verbose

  ! TO BE CHANGED
  real(wp) :: diffscale
  
  !---------------------------------------------------------------------------

  ! number of stats TO BE CHANGED -> future module mod_stats.f90
  integer  :: nstat,nstat_v,nstat_w

  ! indicator TO BE CHANGED
  integer  :: idepart
  ! indicator TO BE CHANGED
  integer :: flowtype
  ! indicator TO BE CHANGED
  logical  :: is_residue
  ! indicator TO BE CHANGED
  logical  :: is_mean0
  ! indicator TO BE CHANGED
  logical  :: is_adiab
  ! indicator TO BE CHANGED
  logical  :: is_init_2D3D
  ! indicator TO BE CHANGED
  logical  :: is_SBP
  ! indicator TO BE CHANGED
  logical  :: is_wall2
  ! indicator TO BE CHANGED
  logical  :: is_comm_onesided
  ! indicator TO BE CHANGED
  logical  :: is_dissip_in_increments
  ! indicator for RFM
  logical :: is_RFM,is_RFM_block,is_init1_RFM
  ! indicator for suction & blowing
  logical :: is_suction_blowing
  ! indicator for user-defined grid
  logical :: is_def_grid
  ! indicator for old way of defining grid <~ To be suppressed
  logical :: is_grid_old
  ! indicator TO BE CHANGED
  logical :: is_stagnation
  ! indicator TO BE CHANGED
  logical :: is_rea
  character(len=5) :: loc_rea
  ! Pre-processing grid indicator TO BE CHANGED
  logical  :: is_read_ex,is_half_cell,is_coarse_grid,is_add_sponge
  ! Pre-processing modes TO BE CHANGED
  logical  :: is_satur_curve,is_LST
  !---------------------------------------------------------------------------

  ! Definition of numerical schemes
  ! ===============================
  ! Scheme parameters
  integer  :: stencil
  logical  :: is_DRP
  ! for viscous terms
  integer  :: iorder_visc
  ! Numerical dissipation
  integer  :: iorder_SF ! obsolete TO BE CHANGED
  logical  :: is_SF
  real(wp) :: dissip_coeff,dissip_shock
  ! Shock capturing  *** NOT IMPLEMENTED YET ***
  logical  :: is_shock
  logical :: is_ducros,is_sw_edoh
  ! -> mix of Jameson's sensor and TVD version of Swanson et al. (1998)
  ! csens=0.01_wp --> activate Jameson's sensor
  ! csens=0.5_wp  --> mix of Jameson/TVD sensors [recommended]
  ! csens=1.0_wp  --> full TVD sensor
  real(wp) :: csens
  ! obsolete TO BE CHANGED
  character(1) :: shock_type        ! obsolete TO BE CHANGED
  logical  :: is_supersonic         ! obsolete TO BE CHANGED

  ! Tam & Dong radiation center (temp) TO BE CHANGED
  real(wp) :: xcr_,ycr_,zcr_
  !
  ! Boundary layer paramaters STBL TO BE CHANGED
  real(wp) :: jdel,deltas_in,deltas_forc
  real(wp) :: Re_inlet
  ! wall forcing (suction & blowing) TO BE CHANGED
  real(wp) :: ampl_sb,omeg_sb,beta_sb,sigma2_sb,x_forc,kr_sb,ki_sb
  real(wp), dimension(:), allocatable :: gx_sb,gz_sb
  
  ! Initial constants for CHAN TO BE CHANGED
  real(wp) :: Rec,Reb,Re_tau,utau,hc,rho_wall
  real(wp) :: rho_bulk,U_bulk
  real(wp) :: forc_rhou_ref
  
  ! constants for TGV
  real(wp) :: mgtot0,ektot0,wrtot0 ! TO BE CHANGED

  !---------------------------------------------------------------------------
  ! Turbulence modelling
  ! ====================
  ! General parameter
  character(len=10) :: turb_model
  ! Indicators for RANS
  logical :: is_RANS ! Activation of RANS
  logical :: is_RANS_adv ! advanced RANS settings
  character(len=10) :: model_RANS,simulation_RANS
  integer :: ndeb_RANS
  logical :: is_slip_in
  ! Indicators for DES
  logical :: is_DES ! Activation of DES
  logical :: is_SLA ! Shear-layer adapted
  character(len=10) :: model_DES
  ! Indicators for LES
  logical :: is_sgs_model ! Subgrid scale modelling
  logical :: is_scale_sim ! Scale-similarity
  integer :: is_filt_Sij, is_filt_nu
  real(wp) :: Cs_SM,Ci_SM
  character(len=6) :: model_LES
  ! Indicators for wall-modeled LES
  logical :: is_wall_model
  character(len=3) :: wm_model_type
  real(wp), dimension(:,:), allocatable :: utau_jmin,utau_jmax

  ! Definition of numerical schemes
  ! ===============================
  ! Scheme parameters
  integer :: stencil_RANS

  !---------------------------------------------------------------------------
  ! Post-processing: parameters in param_pp.ini TO BE CHANGED and included in param.ini
  ! ================
  integer :: type_pp
  integer :: idem_pp,iend_pp,jdem_pp,jend_pp ! if post-processing not over the whole block

end module mod_constant
