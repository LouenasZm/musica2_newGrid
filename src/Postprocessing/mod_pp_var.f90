!==============================================================================
module mod_pp_var
!==============================================================================
  !> Module for post-processing variables
!==============================================================================
  use precision
  use mod_constant
  use mod_block         ! <- for bloc type
  use mod_grid          ! <- for nplr, ipr
  implicit none
  ! ---------------------------------------------------------------------------
  ! parameters in param_pp.ini
  ! integer :: type_pp,type_data
  integer :: type_data
  ! integer :: nplr,ipr ! plane number read
  integer :: nl1,nl2,nl3 ! local plane dimensions
  integer :: ng1,ng2,ng3 ! global plane dimensions
  integer :: ndom1,ndom2 ! number of procs in spatial direction
  integer, dimension(:), allocatable :: coordx1,coordx2
  character(len=30) :: name_output
  real(wp) :: time_ini ! initial time for samples [useful if is_timestamp=T]
  real(wp) :: dt_spec
  ! ---------------------------------------------------------------------------
  
  ! ---------------------------------------------------------------------------
  ! Type for spectra paramaters in one direction
  ! ---------------------------------------------------------------------------
  type param_dir
     character :: name           ! name of direction: t,x,y,z
     integer :: i                ! index range for direction
     integer :: lbloc            ! block length  
     character :: type_win       ! name of window
     real(wp) :: param_win       ! parameter for windowing
     real(wp) :: Cw              ! windowing factors     
  end type param_dir

  ! ---------------------------------------------------------------------------
  ! Type for spectra paramaters in one direction
  ! ---------------------------------------------------------------------------
  type spectra
     integer :: dim                ! 1D, 2D or 3D spectra
     integer :: nvar               ! number of variables
     logical :: is_onesided=.true. ! one-sided or two-sided spectra
     logical :: is_overlap         ! time segment overlapping for Welch's method
     integer :: loverlap           ! overlap length for Welsh's method
     integer :: nbloc              ! number of blocks
     logical :: is_capon=.false.   ! use Capon's spectral estimator
     integer :: ncapon             ! Capon's filter order length
     real(wp) :: dk                ! spectral resolution: product dk1*dk2*..
     character(6), dimension(20) :: varname ! name of variables (max:20)
     type(param_dir), dimension(:), allocatable :: d ! spectra paramaters per direction
  end type spectra

  ! define spectra type
  type(spectra) :: sp
  ! define spectra averaging directions
  integer :: ndir_av ! number of directions for averaging spectra
  character, dimension(:), allocatable :: dir_av !  name dir averaging
  integer, dimension(:), allocatable :: i_av     ! index dir averaging
  integer :: n_av     ! nb samples for averaging
  ! non-inhomogeneous direction
  integer :: i_in,n_in,ng_in ! direction, local size & global size
  integer :: COMM_in ! communicator
  character :: name_in ! name
  ! new MPI dim for transposition (parallel multi-D spectra)
  integer :: ngt,nt_txt,nt_tzt
  ! number of procs/index in spectral direction (1D spectrum)
  integer :: ndom_sp,i_sp
  ! number of samples for averaging
  integer :: ng_av

  ! indicator for modal analysis (write selected modes in mod_pp_sp_modal)
  logical :: is_modal

  ! slices for multiD-spectra
  logical :: is_slice
  integer :: n_sl(3) ! number for each directions
  integer, dimension(:), allocatable :: i1_sl,i2_sl,i3_sl ! slice indices

  
  ! ---------------------------------------------------------------------------
  ! real variable arrays
  real(wp), dimension(:,:,:), allocatable :: var_r,varm_r,var0_r ! real 3D var
  real(wp), dimension(:,:,:), allocatable :: var_rt ! transpose real 3D var
  ! real variable arrays
  complex(wp), dimension(:,:,:), allocatable :: var_c  ! cplx 3D var
  complex(wp), dimension(:,:,:), allocatable :: var_ct ! transpose cplx 3D var
  ! global variance for multiD normalization
  real(wp) :: rms2,rms2m ! rms squared
  real(wp), dimension(:), allocatable :: rms2_in,rms2m_in ! rms with inhomogeneous dir.
  real(wp), dimension(:), allocatable :: skew_in,skewm_in ! skewness
  real(wp), dimension(:), allocatable :: kurt_in,kurtm_in ! kurtosis (or flatness)
  real(wp), dimension(:,:), allocatable :: rms2_1var ! rms for substract 1 var
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  ! time param
  integer :: nclass
  ! ---------------------------------------------------------------------------
  ! Variables for PP FST
  logical :: is_check_fstt_vl,is_snap_pp_fstt,is_check_stats,is_pp_discr,&
             is_pp_streaks,is_restart_pp,is_lim_fst
  integer :: nout ! for loop
  integer :: nsnap_pp_fstt ! number of snapshots for outputs
  integer :: npoints_pp,nlines_pp,nplanes_pp,nvolumes_pp     ! number of points, lines, planes & volumes
  integer :: nit_beg,nit_end,nit_restart_pp,nbeg_pp,nend_pp,noutputs ! number of outputs to be post-processed
  integer :: lvl_th ! levels for thresold
  integer :: ngh_pp,nkernel_k ! Ghost points + kernel extension in k direction
  integer :: ndx_pp,nfx_pp,ndy_pp,nfy_pp,ndz_pp,nfz_pp ! interior index bounds
  integer :: ndxt_pp,nfxt_pp,ndyt_pp,nfyt_pp,ndzt_pp,nfzt_pp ! index bounds extended
  integer :: ndxtpngh_pp,nfxtmngh_pp ! index bounds extended+ngh_pp
  integer :: nit_interp,nil_interp,nj_interp,nkt_interp,nk_interp,ndivx,ndivy
  integer :: nstreaks_loc,nstreaks_merge,nstreaks_final
  integer :: nx_full,nx_f_deb,idem_f_pp,iend_f_pp
  integer(8) :: ivol_2read
  real(wp) :: u_lim, h_kern, max_D, dist_max_x, dist_max_yz, dist_max_extr, l_min_str
  real(wp) :: prct_kernel ! Pourcentage of kernel limit
  type(bloc), dimension(:), allocatable :: bl_glob ! save block dimensions to read stats
  integer, dimension(:), allocatable :: ipt_2_isn_pp,ili_2_isn_pp,ipl_2_isn_pp,ivl_2_isn_pp
  integer, dimension(:), allocatable :: j_1p5d99 ! for limit of 2*d99
  integer, dimension(:), allocatable :: neighbor_E_fstt ! for EAST mpi neighbors
  integer, dimension(:), allocatable :: nit_interp_g,nil_interp_g,nj_interp_g,nkt_interp_g,nk_interp_g
  integer, dimension(:,:), allocatable :: j_edge,j_edge2 ! for limit of BL
  integer, dimension(:,:), allocatable :: born_kern_j,born_kern_j2 ! kernel extension in j direction
  integer, dimension(:,:,:), allocatable :: bk_ic,bk_jc,bk_jc2 ! kernel extension in i & j direction in curvilinear
  real(wp), dimension(:), allocatable :: U0e, d99 ! for 99% BL thickness
  real(wp), dimension(:), allocatable :: x_interp,y_interp,z_interp ! grid for interpolation
  real(wp), dimension(:,:), allocatable :: xc_interp,yc_interp ! grid for interpolation in curvilinear
  real(wp), dimension(:,:), allocatable :: xgc_f,ygc_f
  real(wp), dimension(:,:), allocatable :: coeff_kernel2 ! coeff kernel on computationnal grid for thresold smoothing
  real(wp), dimension(:,:,:), allocatable :: uut,uun ! tangential and normal velocity fluctuations
  real(wp), dimension(:,:,:), allocatable :: coeff_kc2 ! coeff kernel on computationnal grid for thresold smoothing in curvilinear
  real(wp), dimension(:,:,:), allocatable :: stats_proc,stats_interp,stats_full ! for stats from general stats
  real(wp), dimension(:,:,:), allocatable :: stats_lam,stats_turb,avg_s_lam,avg_s_turb,stats_cpt ! lam/turb stats
  real(wp), dimension(:,:,:), allocatable :: stats_tot,avg_s_tot ! lam/turb stats check_up
  real(wp), dimension(:,:,:), allocatable :: uu_interp, vv_interp, ww_interp ! for interpolated var
  real(wp), dimension(:,:,:), allocatable :: uut_interp, uun_interp ! for interpolated var in curv
  real(wp), dimension(:,:,:), allocatable :: uu_fluct ! for fluctuating streamwise velocity
  real(wp), dimension(:,:,:), allocatable :: extr_density,extr_density2,ltbl ! for lamin-turb discrimination
  real(wp), dimension(:,:,:,:), allocatable :: coeff_kernel
  real(wp), dimension(:,:,:,:,:), allocatable :: coeff_kc
  logical, dimension(:), allocatable :: perio_z_loc

end module mod_pp_var
