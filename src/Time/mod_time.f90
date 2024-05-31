! ==============================================================================
module mod_time
! ==============================================================================
  !> Module for Time variables
! ==============================================================================
  use precision
  use PaScaL_TDMA
  implicit none
  ! ============================================================================
  ! Time loop
  ! =========
  ! number of iterations OR max nondimensional time
  integer  :: nmax
  real(wp) :: timemax
  ! max allowed CPU time (for time-limited classes)
  real(wp) :: cpumax
  ! current time
  real(wp) :: time   ! dimensional
  real(wp) :: tstar  ! nondimensional
  real(wp) :: tscale ! time scale
  ! counters
  integer  :: ntime,ntotal
  integer  :: ntotal_old ! <- solver & info ? useful in module ? TO BE CHANGED
  !!!! real(wp) :: time_threshold ! <- commented in setupref (STBL case)
  ! ============================================================================
  ! Time step
  ! =========
  real(wp) :: CFL    ! Courant-Friedrichs-Lewy (CFL) number
  real(wp) :: deltat ! dimensional
  real(wp) :: dtstar ! nondimensional
  ! indicator for variable time step or local timestep
  logical  :: is_dtvar,is_dtlocal
  integer  :: neval_dt ! if is_dtvar, evaluation of dt each neval_dt
  real(wp) :: c_dtloc ! coeff. between 0 and 1, 0~>global, 1~>fully local
  ! local timestep
  real(wp), dimension(:,:,:), allocatable :: dt_local
  real(wp), dimension(:,:,:), allocatable :: cfl_i,cfl_j,cfl_k
  ! ============================================================================
  ! Time for I/O
  ! ============
  ! iteration after which stats are computed
  integer  :: ndeb,ndeb2
  ! frequency of print (screen/stats/plane/restart)
  integer  :: nprint,freq_stats,freq_volume,freq_plane,freq_line,freq_point,freq_field,freq_plane_check
  logical  :: printstep(3)
  real(wp) :: dtprint(3)
  ! ============================================================================
  ! measure of the elapsed CPU time
  ! ===============================
  real(wp) :: t_start,t_end,cputot,cpurun,cputime
  real(wp) :: time_after,time_before ! <- in solver ? useful in module ? TO BE CHANGED
  ! ============================================================================
  ! Runge-Kutta (RK) schemes
  ! ========================
  ! number of RK substeps
  integer  :: nrk
  ! counters for RK steps
  integer  :: irk 
  ! version of RK
  integer  :: vrk
  ! coefficients of RK schemes
  real(wp), dimension(:), allocatable :: crk,ck
  real(wp), dimension(:), allocatable :: Ark,Brk
  ! theoretical CFL limit
  real(wp) :: cfl_limit
  ! ============================================================================
  ! Implicit Residual Smoothing (IRS)
  ! ===========================
  ! IRS order
  integer  :: iirs
  ! IRS indicators
  logical  :: is_irs ! .true. if IRS in at least one direction
  logical  :: is_irs_i,is_irs_j,is_irs_k ! indicators per direction
  logical  :: is_irs_i_master,is_irs_j_master,is_irs_k_master ! per direction
              ! for "master" block (read in param.ini)       
  ! IRS coefficients
  real(wp) :: theta_irs1,theta_irs2,theta_irs4,theta_irs6,theta_irs8
  ! IRS parallelisation
  integer :: type_para
  ! local CFL
  real(wp), dimension(:,:,:), allocatable :: cfl_l
  ! ============================================================================
  
  real(wp), dimension(:), allocatable :: idx_irs,idy_irs,idz_irs
  real(wp), dimension(:,:), allocatable :: ijacob_irs
  ! SCALAPACK variables
  real(wp), dimension(:), allocatable       :: WORK
  real(wp), dimension(:), allocatable       :: AF_x,D_x,DL_x,DU_x
  real(wp), dimension(:), allocatable       :: AF_y,D_y,DL_y,DU_y
  real(wp), dimension(:), allocatable       :: AF_z,D_z,DL_z,DU_z
  real(wp), dimension(:,:), allocatable     :: A_x,A_y,A_z
  integer, dimension(3) :: ictxt
  integer, dimension(9) :: DESCA_x, DESCA_y, DESCA_z, DESCB_x, DESCB_y, DESCB_z
  integer :: LAF_x,LAF_y,LAF_z
  integer :: IB,JA,LWORK,bw,NRHS
  ! PaScaL tridiagonal solver
  real(wp), dimension(:,:,:), allocatable   :: Dp_x,DLp_x,DUp_x,Dp_y,DLp_y,DUp_y,Dp_z,DLp_z,DUp_z
  real(wp), dimension(:,:,:), allocatable   :: Dp_ksi,DLp_ksi,DUp_ksi,Dp_eta,DLp_eta,DUp_eta
  type(ptdma_plan_many) :: p_many_x, p_many_y, p_many_z

contains

  !===============================================================
  subroutine init_time
  !===============================================================
    !> Initialize variables to measure the elapsed CPU time
  !===============================================================
    use mpi
    implicit none
    ! ------------------------------------------------------------
    ! ------------------------------------------------------------

    t_start = MPI_WTIME()
    cpurun  = 0.0_wp
    cputime = 0.0_wp
    cputot  = 0.0_wp

  end subroutine init_time

  !===============================================================
  subroutine init_RK
  !===============================================================
    !> Initialize Runge-Kutta scheme coefficients
  !===============================================================
    use warnstop
    implicit none
    ! ------------------------------------------------------------
    real(wp) :: g1,g2,g3,g4,g5,g6
    ! ------------------------------------------------------------

    ! nrk is read in param.ini
    ! ------------------------

    ! allocations
    ! -----------
    allocate(crk(nrk),ck(nrk))
    allocate(Ark(nrk),Brk(nrk))
    
    ! Coefficients of Runge-Kutta schemes
    ! -----------------------------------
    select case (nrk)
    case (3)
       ! RK 4 sub-steps (o4 linear, o2 nonlinear)
       crk(1) = 1.0_wp/2.0_wp
       crk(2) = 1.0_wp/2.0_wp
       crk(3) = 1.0_wp

       ck(1) = 0.0_wp
       ck(2) = crk(1)
       ck(3) = crk(2)

       ! theoretical CFL limit for FDo10/SFo10
       ! [value=1.598457286432161 in stabmax.m]
       cfl_limit=1.59
       
       ! SSP-RK3(3) in Williamson 2-N form (close to Shu & Osher RK3)
       Ark(1) =  0.
       Ark(2) = -2.915493957701923
       Ark(3) =  0.

       Brk(1) = 0.924574112262461
       Brk(2) = 0.287712943868770
       Brk(3) = 0.626538293270800

    case (4)
       ! RK 4 sub-steps (o4 linear, o2 nonlinear)
       crk(1) = 1.0_wp/4.0_wp
       crk(2) = 1.0_wp/3.0_wp
       crk(3) = 1.0_wp/2.0_wp
       crk(4) = 1.0_wp

       ck(1) = 0.0_wp
       ck(2) = crk(1)
       ck(3) = crk(2)
       ck(4) = crk(3)

       ! theoretical CFL limit for FDo10/SFo10
       ! [value=1.598457286432161 in stabmax.m]
       cfl_limit=1.59
       
       ! RK4 in Williamson 2-N form: does not work ???
       Ark(1) =  0.
       Ark(2) = -3./4.
       Ark(3) = -2./3.
       Ark(4) = -1./2.

       Brk(1) = 1./4.
       Brk(2) = 1./3.
       Brk(3) = 1./2.
       Brk(4) = 1.
       
    case (5)
       ! RK 5 sub-steps (o4 lin. optimized, o2 nonlinear)
       ! [Bogey & Bailly JCP 2004 ~> DRP optimization]
       g1 =   1.000000000000_wp
       g2 =   0.500000000000_wp
       g3 =   0.165435011584000_wp
       g4 =   3.941077299200001e-2_wp
       g5 =   7.268879360000001e-3_wp
       
!!$       g1 =   1.000000000000_wp
!!$       g2 =   0.500000000000_wp
!!$       g3 =   0.165250353664000_wp
!!$       g4 =   3.937258598400000e-2_wp
!!$       g5 =   7.149096448000001e-3_wp
!!$
!!$       ! RK5 JST
!!$       g1 =   1.0_wp
!!$       g2 =   1.0_wp/2.0_wp
!!$       g3 =   3.0_wp/16.0_wp
!!$       g4 =   1.0_wp/32.0_wp
!!$       g5 =   1.0_wp/128.0_wp
!!$g3 = 1.0_wp/6.0_wp       
!!$g4 = 3.999457612800000E-002;
!!$g5 = 7.829006335999998E-003;
              
       crk(5) = g1
       crk(4) = g2
       crk(3) = g3/ crk(4)
       crk(2) = g4/(crk(3)*crk(4))
       crk(1) = g5/(crk(2)*crk(3)*crk(4))

!!$       ! RK5 HALE JST
!!$       crk(1) = 1.0_wp/4.0_wp
!!$       crk(2) = 1.0_wp/6.0_wp
!!$       crk(3) = 3.0_wp/8.0_wp
!!$       crk(4) = 1.0_wp/2.0_wp
!!$       crk(5) = 1.0_wp

       ck(1) = 0.0_wp
       ck(2) = crk(1)
       ck(3) = crk(2)
       ck(4) = crk(3)
       ck(5) = crk(4)

       ! theoretical CFL limit for FDo10/SFo10
       ! [value=1.719020100502513 in stabmax.m]
       cfl_limit=1.72
       
       ! SSP-RK5(4) in Williamson 2-N form
       Ark(1) =  0.
       Ark(2) = -4.344339134485095
       Ark(3) =  0.
       Ark(4) = -3.770024161386381
       Ark(5) = -0.046347284573284

       Brk(1) = 0.713497331193829
       Brk(2) = 0.133505249805329
       Brk(3) = 0.713497331193829
       Brk(4) = 0.149579395628565
       Brk(5) = 0.384471116121269

    case (6)
       ! RK 6 sub-steps (o4 lin. optimized, o2 nonlinear)
       ! [Bogey & Bailly JCP 2004 ~> DRP optimization]
       g1 = 1.000000000000_wp
       g2 = 0.500000000000_wp
       g3 = 0.165919771368_wp
       g4 = 0.040919732041_wp
       g5 = 0.007555704391_wp
       g6 = 0.000891421261_wp

!!$       ! [Berland, Bogey & Bailly ~> gamma_i optimized 4th-order non-linear]
!!$       g1 = 1.0_wp
!!$       g2 = 1.0_wp/2.0_wp
!!$       g3 = 1.0_wp/6.0_wp
!!$       g4 = 1.0_wp/24.0_wp
!!$       g5 = 0.00781005_wp
!!$       g6 = 0.00132141_wp
!!$
!!$       ! RK6 HALE Hixon ~> gamma_i
!!$       g1 = 1.0_wp
!!$       g2 = 1.0_wp/2.0_wp
!!$       g3 = 1.0_wp/6.0_wp
!!$       g4 = 1.0_wp/24.0_wp
!!$       g5 = 0.00556351_wp
!!$       g6 = 0.00092671_wp
!!$
!!$       ! RK6 Calvo ~> gamma_i
!!$       g1 = 1.0_wp
!!$       g2 = 1.0_wp/2.0_wp
!!$       g3 = 1.0_wp/6.0_wp
!!$       g4 = 1.0_wp/24.0_wp
!!$       g5 = 0.00785333;
!!$       g6 = 0.00094889;

       crk(6) = g1
       crk(5) = g2
       crk(4) = g3/ crk(5)
       crk(3) = g4/(crk(4)*crk(5))
       crk(2) = g5/(crk(4)*crk(3)*crk(5))
       crk(1) = g6/(crk(4)*crk(3)*crk(5)*crk(2))

       ck(1) = 0.0_wp
       ck(2) = crk(1)
       ck(3) = crk(2)
       ck(4) = crk(3)
       ck(5) = crk(4)
       ck(6) = crk(5)

       ! theoretical CFL limit for FDo10/SFo10
       ! [value=1.899864321608040 in stabmax.m]
       cfl_limit=1.89
       
       ! rk46-Berland
       Ark(1) =  0.
       Ark(2) = -0.737101392796
       Ark(3) = -1.634740794341
       Ark(4) = -0.744739003780
       Ark(5) = -1.469897351522
       Ark(6) = -2.813971388035

       Brk(1) = 0.032918605146
       Brk(2) = 0.823256998200
       Brk(3) = 0.381530948900
       Brk(4) = 0.200092213184
       Brk(5) = 1.718581042715
       Brk(6) = 0.27

       !ck(1) = 0.
       !ck(2) = 0.032918605146
       !ck(3) = 0.249351723343
       !ck(4) = 0.466911705055
       !ck(5) = 0.582030414044
       !ck(6) = 0.847252983783
       
    case default
       call mpistop('Wrong value for nrk!', 0)
    end select

  end subroutine init_RK

end module mod_time
