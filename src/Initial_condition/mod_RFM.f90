!=================================================================================
module mod_RFM
!=================================================================================
  !> author: AB & XG
  !> date: February 2024
  !> Generation of Random Fourier Modes (RFM) for injection of:
  !>   I - Synthetic freestream turbulence
  !>   II - Turbulent boundary layers
  !> Restricted to BC imin
!=================================================================================
  use mod_constant
  use mod_mpi_part
  use mod_flow
  use mod_eos
  use mod_block
  use mod_utils
  implicit none
  !-------------------------------------------------------------------------------
  ! RFM application
  logical :: is_RFM_FST
  logical, dimension(2) :: is_RFM_TBL
  ! parameters to define RFM modes
  integer :: Nmode             ! number of RFM modes
  integer :: ndy_RFM,nfy_rfm ! global indices to specify the zone where to inject RFM
  real(wp) :: xkmin,xkmax,xkmx ! min, max and peak wavenumbers for turbulence spectrum xkkol
  real(wp) :: xkkol            ! wavenumbers associated to Kolmogorov's length for turbulence spectrum
  real(wp) :: Lf_int           ! Integral length scale [for von Karman spectrum]
  character(3) :: kdist        ! wavenumber distribution
  real(wp) :: Tu_RFM, ampl_RFM ! amplitude of RFM modes in FST (ampl_RFM is dimensional)
  real(wp) :: Tu_TBL, ampl_TBL ! amplitude of RFM modes in TBL (ampl_RFM is dimensional)
  character(1) :: time_turb    ! turbulence time evolution
  ! logical :: is_convec         ! mean flow convection
  real(wp) :: u_m,v_m,w_m      ! mean flow convection
  real(wp) :: u_m_TBL,v_m_TBL,w_m_TBL ! mean flow convection for RFM injected in TBL
  logical :: is_xkm_given
  ! parameters to define turbulence anisotropy
  character(1) :: anisotropy
  character(3) :: base_FST, base_TBL
  real(wp) :: Re_prof_FST
  real(wp), dimension(2) :: Re_prof_TBL
  real(wp), dimension(:), allocatable :: u_rms,v_rms,w_rms
  ! check convergence of rms profiles
  logical :: is_check_conv
  ! velocity field initialization
  logical :: is_field
  ! velocity field initialization
  logical :: is_init_modes
  ! TBL application
  real(wp) :: Lf_TBL
  !-------------------------------------------------------------------------------
  ! RFM modes
  real(wp), dimension(:), allocatable :: omn ! turbulence angular frequency
  real(wp), dimension(:), allocatable :: psi ! mode phase
  real(wp), dimension(:), allocatable :: xk1,xk2,xk3 ! unitary wavenumber components
  real(wp), dimension(:), allocatable :: sigma1,sigma2,sigma3 ! mode directions
  ! anisotropic field
  real(wp), dimension(:,:), allocatable :: xk1s,xk2s,xk3s ! scaled wavenumber components
  real(wp), dimension(:,:), allocatable :: sigma1s,sigma2s,sigma3s ! scaled mode directions
  real(wp), dimension(:,:), allocatable :: wr   ! eigenvalues
  real(wp), dimension(:,:,:), allocatable :: vr ! right eigenvectors
  ! TBL RFM modes
  real(wp), dimension(:,:), allocatable :: psi_TBL ! mode phase
  real(wp), dimension(:,:), allocatable :: omn_TBL ! turb. evolution frequency
  real(wp), dimension(:,:,:), allocatable :: xk1s_TBL,xk2s_TBL,xk3s_TBL ! scaled wavenumber components
  real(wp), dimension(:,:,:), allocatable :: sigma1s_TBL,sigma2s_TBL,sigma3s_TBL ! scaled mode directions
  real(wp), dimension(:,:,:,:), allocatable :: vr_TBL ! right eigenvectors
  ! writting of inlet planes
  integer, dimension(:), allocatable :: inlet_ip
  integer :: nob_inp
  ! Damping if inlet injection zone restricted
  integer :: nrfm_proc
  logical :: is_damping
  real(wp), dimension(2) :: h_damp
  real(wp), dimension(:), allocatable :: damping_coeff
  real(wp), dimension(:,:), allocatable :: damping_coeff3
  ! Reference origin for inlet turbulence
  real(wp) :: xg_in
  ! Wall coordinates
  real(wp), dimension(2) :: yw_BCj
  ! TBL
  real(wp) :: delta ! ?????????????? TO BE CHANGED ?????????????????? local variable ???????
  real(wp), dimension(2) :: d99_TBL
  !-------------------------------------------------------------------------------
  integer :: iproc_leader_rfm
  logical, dimension(:), allocatable :: is_RFM_blocks

  interface
     !===============================================================================
     ! initialization routines
     !===============================================================================
     module subroutine init_RFM
     end subroutine init_RFM
     !===============================================================================
     module subroutine display_RFM_summary
     end subroutine display_RFM_summary
     !===============================================================================
     module subroutine read_param_RFM
     end subroutine read_param_RFM
     !===============================================================================
     module subroutine init_RFM_planes
     end subroutine init_RFM_planes
     !===============================================================================
     ! routines to enter RFM modes in inlet plane(s)
     !===============================================================================
     ! -> for Tam & Dong's boundary conditions
     !===============================================================================
     module subroutine disturb_inlet_RFM_TamDong_imin(vg,ut_in,vt_in,wt_in)
       real(wp), dimension(ngh,ny,nz), intent(in) :: vg
       real(wp), dimension(ngh,ny,nz), intent(out) :: ut_in,vt_in,wt_in
     end subroutine disturb_inlet_RFM_TamDong_imin
     !===============================================================================
     module subroutine disturb_inlet_RFM_TamDong_imin_jmin(vg,ut_in,vt_in,wt_in)
       real(wp), dimension(ngh,ngh,nz), intent(in) :: vg
       real(wp), dimension(ngh,ngh,nz), intent(out) :: ut_in,vt_in,wt_in
     end subroutine disturb_inlet_RFM_TamDong_imin_jmin
     !===============================================================================
     module subroutine disturb_inlet_RFM_TamDong_imin_jmax(vg,ut_in,vt_in,wt_in)
       real(wp), dimension(ngh,ngh,nz), intent(in) :: vg
       real(wp), dimension(ngh,ngh,nz), intent(out) :: ut_in,vt_in,wt_in
     end subroutine disturb_inlet_RFM_TamDong_imin_jmax
     !===============================================================================
     module subroutine disturb_inlet_RFM_TamDong1pt_imin(vg,ut_in,vt_in,wt_in)
       real(wp), dimension(ny,nz), intent(in) :: vg
       real(wp), dimension(ny,nz), intent(out) :: ut_in,vt_in,wt_in
     end subroutine disturb_inlet_RFM_TamDong1pt_imin
     !===============================================================================
     ! -> for characteristics & Riemann (turbomachinery)
     !===============================================================================
     module subroutine disturb_inlet_RFM_charac(ut_in,vt_in,wt_in)
       real(wp), dimension(ny,nz), intent(out) :: ut_in,vt_in,wt_in
     end subroutine disturb_inlet_RFM_charac
     !===============================================================================
     module subroutine disturb_inlet_RFM_turb(u_in,v_in,w_in)
       real(wp), dimension(ny,nz), intent(out) :: u_in,v_in,w_in
     end subroutine disturb_inlet_RFM_turb
     !===============================================================================
     ! initialization of field with RFM modes
     !===============================================================================
     module subroutine init_vel_RFM
     end subroutine init_vel_RFM
     !===============================================================================
     module subroutine compute_RFM_check
     end subroutine compute_RFM_check
     !===============================================================================
     module subroutine write_RFM_rms(varms,ind)
       integer, intent(in) :: ind
       real(wp), dimension(ny) :: varms
     end subroutine write_RFM_rms
     !===============================================================================
  end interface

end module mod_RFM
