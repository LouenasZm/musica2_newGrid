!================================================================================
module mod_grid
!================================================================================
  !> Module for grids and metrics
!================================================================================
  use mod_ngh
  implicit none
  !------------------------------------------------------------------------------
  ! Grid definition
  ! ===============
  ! Global grid sizes for a block
  ! -----------------------------
  integer :: ngx,ngy,ngz
  
  ! Local grid sizes for proc within a block
  ! ----------------------------------------
  integer :: nx,ny,nz
  
  ! Local grid sizes with ghost points ngh (inviscid part)
  ! ------------------------------------------------------
  integer :: nx1,nx2,ny1,ny2,nz1,nz2
  
  ! Local grid sizes with ghost points ngh_v (viscous part)
  ! -------------------------------------------------------
  integer :: nx1_v,nx2_v,ny1_v,ny2_v,nz1_v,nz2_v
  
  ! Local grid sizes with ghost points ngh_irs (IRS)
  ! ------------------------------------------------
  integer :: nx1_irs,nx2_irs,ny1_irs,ny2_irs,nz1_irs,nz2_irs
  
  ! Cartesian grid (global to a block)
  ! ---------------------------------------------
  real(wp), dimension(:), allocatable :: xg,yg,zg

  ! Cartesian grid (local to proc within a block)
  ! ---------------------------------------------
  real(wp), dimension(:), allocatable :: x,y,z

  ! 2-D Curvilinear grid (global to a block)
  ! ----------------------------------------
  real(wp), dimension(:,:), allocatable :: xgc,ygc
  ! extended global grid (redundant with xgc, to be simplified)
  real(wp), dimension(:,:), allocatable :: xgce,ygce

  ! 2-D Curvilinear grid (local to proc within a block)
  ! ---------------------------------------------------
  real(wp), dimension(:,:), allocatable :: xc,yc

  ! 3-D Curvilinear grid (global to a block)
  ! ----------------------------------------
  real(wp), dimension(:,:,:), allocatable :: xgc3,ygc3,zgc3
  ! extended global grid (redundant with xgc3, to be simplified)
  real(wp), dimension(:,:,:), allocatable :: xgc3e,ygc3e,zgc3e

  ! 3-D Curvilinear grid (local to proc within a block)
  ! ---------------------------------------------------
  real(wp), dimension(:,:,:), allocatable :: xc3,yc3,zc3
  !------------------------------------------------------------------------------
  ! Cartesian grid metrics
  ! =======================
  ! Local metrics for proc within a block (inviscid part)
  real(wp), dimension(:), allocatable :: idx,idy,idz
  ! Local metrics for proc within a block (viscous part)
  real(wp), dimension(:), allocatable :: idx_v,idy_v,idz_v
  ! boundary metrics at imin or imax    
  real(wp) :: idx1_imin,idx2_imin,idx4_imin,idx6_imin,idx8_imin
  real(wp) :: idx1_imax,idx2_imax,idx4_imax,idx6_imax,idx8_imax
  ! boundary metrics at kmin or kmax    
  real(wp) :: idy1_jmin,idy2_jmin,idy4_jmin,idy6_jmin,idy8_jmin
  real(wp) :: idy1_jmax,idy2_jmax,idy4_jmax,idy6_jmax,idy8_jmax
  ! boundary metrics at kmin or kmax    
  real(wp) :: idz1_kmin,idz2_kmin,idz4_kmin,idz6_kmin,idz8_kmin
  real(wp) :: idz1_kmax,idz2_kmax,idz4_kmax,idz6_kmax,idz8_kmax
  !------------------------------------------------------------------------------
  ! Boundary normals
  ! ================
  ! BC imin
  real(wp) :: sgn_imin
  real(wp), dimension(:,:), allocatable :: nxn_imin,nyn_imin,nzn_imin,inksi_imin
  ! BC imax
  real(wp) :: sgn_imax
  real(wp), dimension(:,:), allocatable :: nxn_imax,nyn_imax,nzn_imax,inksi_imax
  ! BC jmin
  real(wp) :: sgn_jmin
  real(wp), dimension(:,:), allocatable :: nxn_jmin,nyn_jmin,nzn_jmin,ineta_jmin
  ! BC jmax
  real(wp) :: sgn_jmax
  real(wp), dimension(:,:), allocatable :: nxn_jmax,nyn_jmax,nzn_jmax,ineta_jmax
  ! BC kmin
  real(wp) :: sgn_kmin
  real(wp), dimension(:,:), allocatable :: nxn_kmin,nyn_kmin,nzn_kmin,inphi_kmin
  ! BC kmax
  real(wp) :: sgn_kmax
  real(wp), dimension(:,:), allocatable :: nxn_kmax,nyn_kmax,nzn_kmax,inphi_kmax
  !------------------------------------------------------------------------------
  ! 2-D curvilinear grid metrics
  ! ============================
  ! 2-D curvilinear local variables
  real(wp), dimension(:,:), allocatable :: x_ksi,y_ksi,x_eta,y_eta,ijacob
  real(wp), dimension(:,:), allocatable :: x_ksi_v,y_ksi_v,x_eta_v,y_eta_v,ijacob_v
  ! gradients of metrics : sqrt[||grad(ksi)||²] & sqrt[||grad(eta)||²]
  real(wp), dimension(:,:), allocatable :: g_ksi,g_eta
  ! wall parameters BC imin
  real(wp), dimension(:), allocatable :: nxndl_imin,nyndl_imin,dl_imin
  real(wp), dimension(:), allocatable :: gksigeta_imin,txeta_imin,tyeta_imin
  ! wall parameters BC imax
  real(wp), dimension(:), allocatable :: nxndl_imax,nyndl_imax,dl_imax
  real(wp), dimension(:), allocatable :: gksigeta_imax,txeta_imax,tyeta_imax
  ! wall parameters BC jmin
  real(wp), dimension(:), allocatable :: nxndl_jmin,nyndl_jmin,dl_jmin
  real(wp), dimension(:), allocatable :: gksigeta_jmin,txksi_jmin,tyksi_jmin
  ! wall parameters BC jmax
  real(wp), dimension(:), allocatable :: nxndl_jmax,nyndl_jmax,dl_jmax
  real(wp), dimension(:), allocatable :: gksigeta_jmax,txksi_jmax,tyksi_jmax
  ! particular metrics for edges
  real(wp), dimension(:,:), allocatable :: x_ksi_imin_jmin,y_ksi_imin_jmin
  real(wp), dimension(:,:), allocatable :: x_eta_imin_jmin,y_eta_imin_jmin
  real(wp), dimension(:,:), allocatable :: x_ksi_imax_jmin,y_ksi_imax_jmin
  real(wp), dimension(:,:), allocatable :: x_eta_imax_jmin,y_eta_imax_jmin
  real(wp), dimension(:,:), allocatable :: x_ksi_imin_jmax,y_ksi_imin_jmax
  real(wp), dimension(:,:), allocatable :: x_eta_imin_jmax,y_eta_imin_jmax
  real(wp), dimension(:,:), allocatable :: x_ksi_imax_jmax,y_ksi_imax_jmax
  real(wp), dimension(:,:), allocatable :: x_eta_imax_jmax,y_eta_imax_jmax
  ! check curvilinear metrics commutations
  real(wp), dimension(:,:,:), allocatable :: x_ksi_eta,y_ksi_eta
  real(wp), dimension(:,:,:), allocatable :: x_eta_ksi,y_eta_ksi
  !------------------------------------------------------------------------------
  ! 3-D curvilinear grid metrics
  ! ============================
  ! 3-D curvilinear local metrics for inviscid fluxes
  real(wp), dimension(:,:,:), allocatable :: ksi_x,ksi_y,ksi_z
  real(wp), dimension(:,:,:), allocatable :: eta_x,eta_y,eta_z
  real(wp), dimension(:,:,:), allocatable :: phi_x,phi_y,phi_z
  real(wp), dimension(:,:,:), allocatable :: ijacob3 ! inverse Jacobian
  ! 3-D curvilinear local metrics for viscous fluxes
  real(wp), dimension(:,:,:), allocatable :: ksi_x_v,ksi_y_v,ksi_z_v
  real(wp), dimension(:,:,:), allocatable :: eta_x_v,eta_y_v,eta_z_v
  real(wp), dimension(:,:,:), allocatable :: phi_x_v,phi_y_v,phi_z_v
  real(wp), dimension(:,:,:), allocatable :: ijacob3_v ! inverse Jacobian
  ! gradients of metrics : sqrt[||grad(ksi)||²],sqrt[||grad(eta)||²],sqrt[||grad(phi)||²]
  real(wp), dimension(:,:,:), allocatable :: g3_ksi,g3_eta,g3_phi
  ! wall parameters BC imin
  real(wp), dimension(:,:), allocatable :: nxnds_imin,nynds_imin,nznds_imin
  real(wp), dimension(:,:), allocatable :: getagksi3_imin,gphigksi3_imin
  ! wall parameters BC imax
  real(wp), dimension(:,:), allocatable :: nxnds_imax,nynds_imax,nznds_imax
  real(wp), dimension(:,:), allocatable :: getagksi3_imax,gphigksi3_imax
  ! wall parameters BC jmin
  real(wp), dimension(:,:), allocatable :: nxnds_jmin,nynds_jmin,nznds_jmin
  real(wp), dimension(:,:), allocatable :: gksigeta3_jmin,gphigeta3_jmin
  ! wall parameters BC jmax
  real(wp), dimension(:,:), allocatable :: nxnds_jmax,nynds_jmax,nznds_jmax
  real(wp), dimension(:,:), allocatable :: gksigeta3_jmax,gphigeta3_jmax
  ! wall parameters BC kmin
  real(wp), dimension(:,:), allocatable :: nxnds_kmin,nynds_kmin,nznds_kmin
  real(wp), dimension(:,:), allocatable :: gksigphi3_kmin,getagphi3_kmin
  ! wall parameters BC kmax
  real(wp), dimension(:,:), allocatable :: nxnds_kmax,nynds_kmax,nznds_kmax
  real(wp), dimension(:,:), allocatable :: gksigphi3_kmax,getagphi3_kmax
  ! tangent BC imin
  real(wp), dimension(:,:), allocatable :: txn_eta_imin,tyn_eta_imin,tzn_eta_imin
  real(wp), dimension(:,:), allocatable :: txn_phi_imin,tyn_phi_imin,tzn_phi_imin
  real(wp), dimension(:,:), allocatable :: tx_eta_imin,ty_eta_imin,tz_eta_imin
  real(wp), dimension(:,:), allocatable :: tx_phi_imin,ty_phi_imin,tz_phi_imin
  ! tangent BC imax
  real(wp), dimension(:,:), allocatable :: txn_eta_imax,tyn_eta_imax,tzn_eta_imax
  real(wp), dimension(:,:), allocatable :: txn_phi_imax,tyn_phi_imax,tzn_phi_imax
  real(wp), dimension(:,:), allocatable :: tx_eta_imax,ty_eta_imax,tz_eta_imax
  real(wp), dimension(:,:), allocatable :: tx_phi_imax,ty_phi_imax,tz_phi_imax
  ! tangent BC jmin
  real(wp), dimension(:,:), allocatable :: txn_ksi_jmin,tyn_ksi_jmin,tzn_ksi_jmin
  real(wp), dimension(:,:), allocatable :: txn_phi_jmin,tyn_phi_jmin,tzn_phi_jmin
  real(wp), dimension(:,:), allocatable :: tx_ksi_jmin,ty_ksi_jmin,tz_ksi_jmin
  real(wp), dimension(:,:), allocatable :: tx_phi_jmin,ty_phi_jmin,tz_phi_jmin
  ! tangent BC jmax
  real(wp), dimension(:,:), allocatable :: txn_ksi_jmax,tyn_ksi_jmax,tzn_ksi_jmax
  real(wp), dimension(:,:), allocatable :: txn_phi_jmax,tyn_phi_jmax,tzn_phi_jmax
  real(wp), dimension(:,:), allocatable :: tx_ksi_jmax,ty_ksi_jmax,tz_ksi_jmax
  real(wp), dimension(:,:), allocatable :: tx_phi_jmax,ty_phi_jmax,tz_phi_jmax
  ! tangent BC kmin
  real(wp), dimension(:,:), allocatable :: txn_ksi_kmin,tyn_ksi_kmin,tzn_ksi_kmin
  real(wp), dimension(:,:), allocatable :: txn_eta_kmin,tyn_eta_kmin,tzn_eta_kmin
  real(wp), dimension(:,:), allocatable :: tx_ksi_kmin,ty_ksi_kmin,tz_ksi_kmin
  real(wp), dimension(:,:), allocatable :: tx_eta_kmin,ty_eta_kmin,tz_eta_kmin
  ! tangent BC kmin
  real(wp), dimension(:,:), allocatable :: txn_ksi_kmax,tyn_ksi_kmax,tzn_ksi_kmax
  real(wp), dimension(:,:), allocatable :: txn_eta_kmax,tyn_eta_kmax,tzn_eta_kmax
  real(wp), dimension(:,:), allocatable :: tx_ksi_kmax,ty_ksi_kmax,tz_ksi_kmax
  real(wp), dimension(:,:), allocatable :: tx_eta_kmax,ty_eta_kmax,tz_eta_kmax

  !------------------------------------------------------------------------------
  ! Index bounds for MPI
  ! ==================== 
  ! index bounds for inviscid fluxes
  ! --------------------------------
  ! index for Eulerian fluxes
  integer :: ndx_e,nfx_e,ndy_e,nfy_e,ndz_e,nfz_e  
  ! index for Eulerian fluxes limited to interior points
  integer :: ndx,nfx,ndy,nfy,ndz,nfz
  ! index for Eulerian fluxes extended to ghost points
  integer :: ndxt,nfxt,ndyt,nfyt,ndzt,nfzt
  
  ! index bounds for viscous fluxes
  ! -------------------------------
  ! index for viscous fluxes
  integer :: ndx_v,nfx_v,ndy_v,nfy_v,ndz_v,nfz_v
  ! index for viscous fluxes limited to interior points
  integer :: ndx_vi,nfx_vi,ndy_vi,nfy_vi,ndz_vi,nfz_vi
  ! index for viscous fluxes extended to ghost points
  integer :: ndxt_v,nfxt_v,ndyt_v,nfyt_v,ndzt_v,nfzt_v
  ! index bounds for velocity gradients (double derivative)
  integer :: ndx_v1,nfx_v1,ndy_v1,nfy_v1,ndz_v1,nfz_v1
  ! index bounds for velocity gradients limited to interior points
  integer :: ndx_v2,nfx_v2,ndy_v2,nfy_v2,ndz_v2,nfz_v2
  ! index bounds for velocity gradients limited to interior points (without ghost cells for double derivatives)
  integer :: ndx_v3,nfx_v3,ndy_v3,nfy_v3,ndz_v3,nfz_v3

  ! index bounds for numerical dissipation
  ! --------------------------------------
  ! all points where dissipation is applied
  integer :: ndx_d,nfx_d,ndy_d,nfy_d,ndz_d,nfz_d
  ! same +/- 1pt for sensor calculation
  integer :: ndx_d1,nfx_d1,ndy_d1,nfy_d1,ndz_d1,nfz_d1
  ! interior points where dissipation is applied
  integer :: ndx_di,nfx_di,ndy_di,nfy_di,ndz_di,nfz_di

  ! index bounds for wall BCs
  ! -------------------------
  ! -> wall at i=cste (parallel dir. are j in 2D curv + z in 3D curv)
  integer :: ndy_imin,nfy_imin,ndy_imax,nfy_imax
  integer :: ndz_imin,nfz_imin,ndz_imax,nfz_imax
  ! -> wall at j=cste (parallel dir. are i in 2D curv + z in 3D curv)
  integer :: ndx_jmin,nfx_jmin,ndx_jmax,nfx_jmax
  integer :: ndz_jmin,nfz_jmin,ndz_jmax,nfz_jmax
  ! -> wall at k=cste (parallel dir. are i,j only in 3D curv)
  integer :: ndx_kmin,nfx_kmin,ndx_kmax,nfx_kmax
  integer :: ndy_kmin,nfy_kmin,ndy_kmax,nfy_kmax

  ! index bounds for boundaries
  ! ---------------------------
  ! all boundary points where TamDong 1 point conditions are applied (except corners)
  integer :: ndx_td1,nfx_td1,ndy_td1,nfy_td1
  integer :: ndx_td1m1,nfx_td1p1,ndy_td1m1,nfy_td1p1
  ! all boundary points where characteristic conditions are applied
  integer :: ndx_c,nfx_c,ndy_c,nfy_c,ndz_c,nfz_c

  !!!!!!!!!!!!!!!!!!!!!!!!
  integer :: ndx1,nfx1,ndy1,nfy1,ndz1,nfz1 ! <- in MOD_ART_VISC *** TO BE CHANGED
  integer :: ndx3,nfx3,ndy3,nfy3,ndz3,nfz3 ! <- in MOD_CHARAC   *** TO BE CHANGED  
  ! degradation de l'ordre du filtrage *** non regression *** TO BE CHANGED
  integer :: ndx_f,nfx_f,ndy_f,nfy_f,ndz_f,nfz_f

  !------------------------------------------------------------------------------
  ! RANS indices
  ! Index bounds for MPI
  ! ====================
  ! index bounds for inviscid fluxes
  ! --------------------------------
  ! index for Eulerian fluxes
  integer :: ndx_e_r,nfx_e_r,ndy_e_r,nfy_e_r,ndz_e_r,nfz_e_r
  ! index for Eulerian fluxes limited to interior points
  integer :: ndx_r,nfx_r,ndy_r,nfy_r,ndz_r,nfz_r
  ! index for Eulerian fluxes extended to ghost points
  integer :: ndxt_r,nfxt_r,ndyt_r,nfyt_r,ndzt_r,nfzt_r

  ! index bounds for viscous fluxes
  ! -------------------------------
  ! index for viscous fluxes
  integer :: ndx_v_r,nfx_v_r,ndy_v_r,nfy_v_r,ndz_v_r,nfz_v_r
  ! index for viscous fluxes limited to interior points
  integer :: ndx_vi_r,nfx_vi_r,ndy_vi_r,nfy_vi_r,ndz_vi_r,nfz_vi_r
  ! index for viscous fluxes extended to ghost points
  integer :: ndxt_v_r,nfxt_v_r,ndyt_v_r,nfyt_v_r,ndzt_v_r,nfzt_v_r
  ! index bounds for velocity gradients (double derivative)
  integer :: ndx_v1_r,nfx_v1_r,ndy_v1_r,nfy_v1_r,ndz_v1_r,nfz_v1_r
  ! index bounds for velocity gradients limited to interior points
  integer :: ndx_v2_r,nfx_v2_r,ndy_v2_r,nfy_v2_r,ndz_v2_r,nfz_v2_r
  ! index bounds for velocity gradients limited to interior points (without ghost cells for double derivatives)
  integer :: ndx_v3_r,nfx_v3_r,ndy_v3_r,nfy_v3_r,ndz_v3_r,nfz_v3_r

  ! index bounds for numerical dissipation
  ! --------------------------------------
  ! all points where dissipation is applied
  integer :: ndx_d_r,nfx_d_r,ndy_d_r,nfy_d_r,ndz_d_r,nfz_d_r
  ! same +/- 1pt for sensor calculation
  integer :: ndx_d1_r,nfx_d1_r,ndy_d1_r,nfy_d1_r,ndz_d1_r,nfz_d1_r
  ! interior points where dissipation is applied
  integer :: ndx_di_r,nfx_di_r,ndy_di_r,nfy_di_r,ndz_di_r,nfz_di_r

  ! index bounds for source terms
  ! -----------------------------
  integer :: ndx_s_r,nfx_s_r,ndy_s_r,nfy_s_r,ndz_s_r,nfz_s_r

  
  !------------------------------------------------------------------------------
  ! Grid parameters
  ! ===============

  real(wp) :: dy0p
  
  real(wp) :: xmin,xmax,deltax,longx
  real(wp) :: ymin,ymax,deltay,longy
  real(wp) :: zmin,zmax,deltaz,longz
  real(wp) :: dxmax

  ! Grid spacing & location at origin of block
  real(wp), dimension(:), allocatable :: deltax0,deltay0,x0,y0
  
  ! block center used in .....
  integer  :: i0,j0,k0 ! used to determine normal in init_bc_inlet_outlet

  ! Number for grid generated in solver
  integer  :: num_grid

  ! Streching parameters in param.ini & param_grid.ini
  integer :: nstretchz
  integer, dimension(:), allocatable :: nrz1,nrz2
  real(wp), dimension(:), allocatable :: rsz
  integer, dimension(:), allocatable :: nstretchx,nstretchy
  
  ! grid type for SRC, CAV, ... TO BE CHANGED
  integer :: igrd

  ! Scaling value for grid
  real(wp) :: Lgrid

  ! For post-processing: To be changed
  ! necessary for reading grid for fstt pp
  integer :: iblc_pp ! Local block to be post-processed
  integer :: nplr,nsr ! vol/pl/line/pt number & snapshot number read for post-processing


end module mod_grid
