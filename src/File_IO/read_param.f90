!===============================================================================
subroutine read_param(paramfile)
!===============================================================================
  !> author: XG & AB
  !> date: February 2024
  !> Reading solver & case parameters in param.ini
!===============================================================================
  use mod_mpi
  use mod_constant
  use mod_coeff_deriv
  use mod_flow
  use mod_rans
  use mod_time
  use mod_io
  use mod_eos
  use mod_tranprop
  use mod_init_flow
  use mod_routines
  use mod_forcing_bulk  ! <~ for forcing parameters
  use mod_eigenmode     ! <~ for entering eigenmodes
  use mod_RFM           ! <~ for using Random Fourier Modes
  use mod_init_irs      ! <~ for irs_version !!!!!!!!!! TEMPORARY !!!!!!!!!!!!!!!!
  use mod_bc_periodicity! <~ to fill Lxp,Lyp,Lzp or theta_period
  use mod_sponge        ! for alpha_sz
  use mod_init_TamDong  ! for xcr, ycr, zcr
  use warnstop
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,isolver,ibc_scheme,type_wall_bc
  character(len=*), intent(in) :: paramfile
  logical :: iexist
  logical :: is_param_case ! for canonical case settings
  ! ---------------------------------------------------------------------------
  ! Temp
  real(wp), dimension(:,:), allocatable :: temp_array
  character(100) :: test_str

  ! Initializations
  ! ===============
  is_2d           = .false.
  is_curv         = .false.
  is_curv3        = .false.
  is_shock        = .false.
  is_IOtec_read   = .false.
  is_IOtec_write  = .false.
  is_timestamp    = .false.
  is_src          = .false.
  is_pulse        = .false.
  is_mean0        = .false.
  is_similarity   = .false.
  is_RFM          = .false.
  is_forcing_bulk = .false.
  is_supersonic   = .false.
  is_eigenmode    = .false.
  is_residue      = .false.
  is_init_2D3D    = .false.
  is_RANS         = .false.
  is_RANS_adv     = .false.
  is_SLA          = .false.
  is_scale_sim    = .false.
  is_wall_model   = .false.
  is_slip_in      = .false.
  is_param_case   = .false.
  is_def_grid     = .false.
  is_BC_ref       = .false.
  is_stagnation   = .false.
  is_rea          = .false.
  is_curv         = .false.
  is_adiab        = .false.
  is_curv3         = .false.
  is_init1_RFM    = .true. ! RFM initialization cut in two part, needs to be .true. initially
  is_LBL          = .false.
  is_adjoint_block = .false.
  is_SBP = .false.
  is_wall2 = .false.
  is_comm_onesided = .false.
  is_dissip_in_increments=.false.
  is_dtlocal=.false.
  is_irs_i=.false.; is_irs_j=.false.; is_irs_k=.false.; is_irs=.false.

  ! Read param.ini file
  ! ===================
  inquire(file=trim(paramfile), exist=iexist)
  if (.not.iexist) then
     call mpistop('Paramfile does not exist!', 0)
  endif

  open(30,file=trim(paramfile))
  rewind(30)
  read(30,*) !===============================================================================================
  read(30,*) !===============================================================================================
  read(30,*) ! MUSICA2 : fill Solver & Case parameters
  read(30,*) !===============================================================================================
  read(30,*) !===============================================================================================
  read(30,*) ! 1/ Specify mode:
  read(30,*) !    --> start, restart, pre-processing, post-processing
  read(30,*) ! 2/ Specify numerical discretization & turbulence modelling:
  read(30,*) !    --> Time discretization & parameters
  read(30,*) !    --> Space discretization
  read(30,*) !    --> Turbulence modelling
  read(30,*) ! 3/ Select proper equation of state & transport equations
  read(30,*) ! 4/ Adjust input/output management parameters
  read(30,*) ! 5/ Specify case parameters:
  read(30,*) !    --> Reference quantities
  read(30,*) !    --> General boundary condition (BC) parameters
  read(30,*) !        ~> BC of block(s) specified in param_blocks.ini
  read(30,*) !    --> Flow initialization
  read(30,*) !    --> Grid:
  read(30,*) !        ~> a) Read from namefile prescribed
  read(30,*) !        ~> b) Generated in the solver (for basic grid)
  read(30,*) !    --> Select flow type indicator
  read(30,*) ! 6/ If pre/post-processing, additional parameters
  read(30,*) !===============================================================================================
  read(30,*) !===============================================================================================
  read(30,*) !                                       Solver parameters
  read(30,*) !===============================================================================================
  read(30,*) ! Restart indicator [0:pre-processing; 1:from_scratch;
  read(30,*) !                    2:from_field; 3:from_interp; 4:post-processing]
  read(30,*) idepart
  read(30,*) ! Solver version [1:2D/3D Cartesian; 2:2D/3D curvilinear 2D; 3:full 3D curvilinear]
  read(30,*) isolver
  if (isolver.eq.3) then
     is_curv3=.true.; is_curv=.false.
  else if (isolver.eq.2) then
     is_curv3=.false.; is_curv=.true.
  else
     is_curv3=.false.; is_curv=.false.
  endif
  read(30,*) ! Multiblock: Adjoint blocks [T: neighboring block interfaces coincide, F: half cell approach]
  read(30,*) is_adjoint_block
  if (is_adjoint_block) call mpistop('Adjoint-block mode needs to be checked',0)
  read(30,*) ! Extended verbose mode: [T: Display of more information, F: summary of vital information]
  read(30,*) verbose
  read(30,*) !-----------------------------------------------------------------------------------------------
  read(30,*) ! Time discretization
  read(30,*) !-----------------------------------------------------------------------------------------------
  read(30,*) ! Runge-Kutta scheme: sub-steps (ss) [3-6]  version [see comments below]
  read(30,*) ! 3 ss: **** NOT IMPLEMENTED yet ****
  read(30,*) ! 4 ss: classic low-storage (4 1), order 4 (4 2), DRP-optimized (4 3)
  read(30,*) ! 5 ss: **** NOT IMPLEMENTED yet ****
  read(30,*) ! 6 ss: DRP-optimized (6 1)
  read(30,*) nrk, vrk
  if ((vrk.ne.1).or.((nrk.ne.4).and.(nrk.ne.6))) call mpistop("Version "//trim(numchar(vrk))//" of Runge-Kutta with "//&
                                                               trim(numchar(nrk))//" sub-steps not implemented yet",0)
  read(30,*) ! Implicit Residual Smoothing [0:explicit; 2:IRS2; 4:IRS4, 6:IRS6, 8:IRS8]
  read(30,*) iirs
  read(30,*) ! Residual smoothing parameter: IRS1 / IRS2 / IRS4 / IRS6 / IRS8
  read(30,*) theta_irs1, theta_irs2, theta_irs4, theta_irs6, theta_irs8
  read(30,*) ! Direction of implicitation [T:yes;F:no]: i / j / k
  read(30,*) is_irs_i, is_irs_j, is_irs_k
  if (iirs==0) then
     ! default initialization if IRS is not activated
     is_irs_i=.false.; is_irs_j=.false.; is_irs_k=.false.; is_irs=.false.
  else
     ! global indicator for IRS
     is_irs=is_irs_i.or.is_irs_j.or.is_irs_k
  endif
  read(30,*) !-----------------------------------------------------------------------------------------------
  read(30,*) ! Time parameters
  read(30,*) !-----------------------------------------------------------------------------------------------
  read(30,*) ! Max number of temporal iterations or Final time
  read(30,*) nmax, timemax
  read(30,*) ! Deltat: is_dtvar, neval_dt (evaluation of dt each neval_dt iteration)
  read(30,*) is_dtvar, neval_dt
  if (is_dtvar) call mpistop('is_dtvar not implemented yet',0)
  read(30,*) ! Deltat: is_dtlocal (indep. eval. for each cells), c_dtloc (0 to 1, 0~>global, 1~>fully local)
  read(30,*) is_dtlocal, c_dtloc
  if ((is_dtlocal).and.(c_dtloc.ne.1.0)) call mpistop('c_dtloc from CM solver version, needs to be implemented back',0)
  ! c_dtloc equivalent of sigma in CM solver version in timestep_new.f90
  ! call MPI_ALLREDUCE(minval(dt_local),dt_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,COMM_global,info)
  ! dt_local(i,j,k) = dt_min*(1.0_wp-sigma) + sigma*alpha*dt_min
  read(30,*) ! Maximal CFL number targeted
  read(30,*) CFL
  read(30,*) ! Hours of simulations (max. elapsed time for limited computer classes)
  read(30,*) cpumax
  cpumax = cpumax*3600.0_wp
  read(30,*) !-----------------------------------------------------------------------------------------------
  read(30,*) ! Space discretization
  read(30,*) !-----------------------------------------------------------------------------------------------
  read(30,*) ! Finite Difference & Dissipation: stencil [3;5;7;9;11 pts]; is_DRP [T: DRP; F:standard FD]
  read(30,*) stencil, is_DRP
  read(30,*) ! Boundary schemes: 0: reduced order; 1: SBP (Summation by Parts) [only 4th-order -> stencil>=9]
  read(30,*) ibc_scheme
  if (ibc_scheme.eq.0) then
     is_SBP=.false.
  else if (ibc_scheme.eq.1) then
     is_SBP=.true.
  else
     call mpistop("Wrong selection of boundary scheme in param.ini",0)
  endif
  read(30,*) ! Order of viscous fluxes [2:second-order on 3pts; 4:fourth-order on 5pts; 0:Euler]
  read(30,*) iorder_visc
  read(30,*) ! Selective Filtering: is_SF [T:SF; F:artifical viscosity]
  read(30,*) is_SF
  read(30,*) ! Switch of Edoh for selective filtering (if is_SF) [T:yes;F:no]
  read(30,*) is_sw_edoh
  read(30,*) ! Filtering or Artificial viscosity amplitude
  read(30,*) ! (between 0 and 1 for SF, recommended value 0.1 / around 1.0 for art.visc)
  read(30,*) dissip_coeff
  read(30,*) ! Indicator is_shock and Coefficient of low-order term
  read(30,*) is_shock, dissip_shock
  read(30,*) ! Shock sensor: Ducros sensor [T:yes;F:no], pressure sensor [0:Jameson; 1:TVD-like; 0.5:mix]
  read(30,*) is_ducros, csens
  if ((csens.lt.0.0_wp).or.(csens.gt.1.0_wp)) call mpistop("Wrong choice for pressure sensor coefficient, must be between 0 and 1.",0)
  csens=max(csens,0.01_wp) ! Protection against divide by zero
  read(30,*) !-----------------------------------------------------------------------------------------------
  read(30,*) ! Turbulence modelling
  read(30,*) !-----------------------------------------------------------------------------------------------
  read(30,*) ! Turbulence modelling: ['N':none; 'RANS'; 'LES'; 'DES'; 'DES-sgs': DES and LES modeling]
  read(30,*) ! --> DES without subgrid scale modelling ('DES') or with sgs modelling ('DES-sgs')
  read(30,*) turb_model
  if (turb_model.eq.'DES-sgs') then
     is_RANS=.true.; is_DES=.true.; is_sgs_model=.true.
  else if (turb_model.eq.'DES') then
     is_RANS=.true.; is_DES=.true.; is_sgs_model=.false.
  else if (turb_model.eq.'RANS') then
     is_RANS=.true.; is_DES=.false.; is_sgs_model=.false.
  else if (turb_model.eq.'LES') then
     is_RANS=.false.; is_DES=.false.; is_sgs_model=.true.
  else if (turb_model.eq.'N') then
     is_RANS=.false.; is_DES=.false.; is_sgs_model=.false.
  endif
  read(30,*) ! RANS model (if RANS or DES modelling): ['SA':Spalart-Allmaras; 'KO':k-omega; ...]
  read(30,*) model_RANS
  read(30,*) ! Iteration where RANS modelling is plugged: ndeb_RANS
  read(30,*) ndeb_RANS
  if (.not.is_RANS) ndeb_RANS=99999999
  read(30,*) ! Convective, diffusive and dissipative RANS fluxes: [3(Rusanov);5 pts]
  read(30,*) stencil_RANS
  read(30,*) ! Advanced settings of RANS [if T, needs param_RANS.ini to define advanced parameters]
  read(30,*) is_RANS_adv
  ! else, do nothing for the moment but at term, initialize all RANS parameters with default values
  read(30,*) ! DES model: ['DES': DES97; 'DDES': detached DES; 'IDDES': improved DDES; ...]
  read(30,*) model_DES
  read(30,*) ! Option: SLA (Shear-Layer-Adaptive) [T:yes;F:no]
  read(30,*) is_SLA
  read(30,*) ! LES model: ['SM':Smagorinsky;'DSM':dynamic SM;'WALE': WALE;'MSM':multiscale SM;
  read(30,*) !             'MSM-ls':MSM large-small; 'MSM-ss':MSM small-small]
  read(30,*) model_LES
  ! To be corrected & completed
  if ((model_LES.eq."SM").or.(model_LES.eq."DSM")) then
     is_filt_Sij=.false.; is_filt_nu=.false.
  else if ((model_LES.eq."MSM").or.(model_LES.eq."MDSM")) then
     is_filt_Sij=.true.; is_filt_nu=.true.
  endif
  read(30,*) ! Smagorinsky constants: Cs and Ci (used for 'SM' or 'MSM')
  read(30,*) Cs_SM, Ci_SM
  read(30,*) ! Options: 1/Scale-similarity [T:yes;F:no]; 2/
  read(30,*) is_scale_sim
  read(30,*) ! Test filter for LES: **** NOT IMPLEMENTED yet ****
  read(30,*) !
  read(30,*) ! Wall-model for LES ['N':none; 'ALG':algebraic; 'ODE':ordinary differential equation; ...]
  read(30,*) wm_model_type
  if (wm_model_type.ne."N") is_wall_model=.true.
  !if (is_wall_model) call mpistop('Wall-model not implemented yet',0)
  read(30,*) !-----------------------------------------------------------------------------------------------
  read(30,*) ! Fluid thermo-physical properties
  read(30,*) !-----------------------------------------------------------------------------------------------
  read(30,*) ! Fluid (-> fluid parameters in feos_(name_of_fluid).ini)
  read(30,*) fluidname
  read(30,*) ! Equation of State (EOS)
  read(30,*) ! Perfect gas          : pfg
  read(30,*) ! van der Waals        : vdw
  read(30,*) ! Martin-Hou           : mah
  read(30,*) ! Span-Wagner polar    : swp
  read(30,*) ! Span-Wagner non-polar: swn
  read(30,*) ! Peng-Robinson        : prs
  read(30,*) ! NIST REFPROP library : ref
  read(30,*) eos_type
  ! Check
  if (fluidname.eq.'air' .and. eos_type.ne.'pfg') then
     call mpistop('Only PFG EoS for air!',0)
  elseif (fluidname.eq.'pp11' .and. eos_type.eq.'pfg') then
     call mpistop('Do not use PFG EoS with PP11!',0) ! TO BE CHANGED
  endif
  read(30,*) ! Viscosity law [S:Sutherland; P:powerlaw; C:Chung-Lee-Starling] (only chung for dense gas)
  read(30,*) visc_type
  ! Check
  if ( visc_type.eq.'C'.and.eos_type.eq.'pfg' .or. &
       visc_type.eq.'S'.and.eos_type.ne.'pfg') then
     call mpistop('Wrong combination of visc law and EoS!', 0)! TO BE CHANGED
  endif
  read(30,*) ! Thermal conductivity law  [C:constant Prandtl; ...] **** NOT IMPLEMENTED yet ****
  read(30,*) ! **** NOT IMPLEMENTED yet **** conductivity_type
  read(30,*) !===============================================================================================
  read(30,*) !                                    Input/Output management
  read(30,*) !===============================================================================================
  read(30,*) ! Filestamp of initial restart file (if not specified or default 0000_0000, restart.bin is used)
  read(30,*) filestamp
  if (idepart.eq.1) filestamp='0000_0000'
  read(30,*) ! Format of inputs: is_IOtec_read [T:Tecplot; F:fortran binary]  **** Att. binary ENDIANESS ****
  read(30,*) is_IOtec_read
  filext_read ='.bin'
  if (is_IOtec_read) filext_read='.plt'
  read(30,*) ! Format of outputs: is_IOtec_write [T:Tecplot; F:fortran binary]  **** Att. binary ENDIANESS ****
  read(30,*) is_IOtec_write
  filext_write='.bin'
  if (is_IOtec_write) filext_write='.plt'
  read(30,*) ! is_timestamp for planes [T:multiple files with timestamp; F:append a single file]
  read(30,*) is_timestamp
  read(30,*) ! Output frequencies: screen / stats / fields
  read(30,*) nprint, freq_stats, freq_field
  read(30,*) ! Snapshot frequencies [if not imposed in param_blocks.ini]: points / lines / planes / volumes
  read(30,*) freq_point, freq_line, freq_plane, freq_volume
  read(30,*) ! Compute residuals (for steady computations) [T:yes;F:no]
  read(30,*) is_residue
  read(30,*) ! Iteration number to start statistics: ndeb
  read(30,*) ndeb
  read(30,*) !===============================================================================================
  read(30,*) !                                        Case parameters
  read(30,*) !===============================================================================================
  read(30,*) ! Flowtype indicator [0:not predefined; 1:TGV; 2:CHIT; 3:CHAN; 4:Periodic hill; 5:STBL;
  read(30,*) !                     6:Cavity flow; 7:Actuator; 8:Cylinder; 9:SHIT; 10:Turbine ...]
  read(30,*) flowtype
  if (flowtype.eq.1) then !!!   TO BE CHANGED
     TGV  = .true.
  elseif (flowtype.eq.2) then
     CHIT = .true.
  elseif (flowtype.eq.3) then
     CHAN = .true.
  elseif (flowtype.eq.4) then
     PHILL = .true.
  elseif (flowtype.eq.5) then
     STBL = .true.
  elseif (flowtype.eq.6) then
     CAV = .true.
  elseif (flowtype.eq.7) then
     ACT = .true.
  elseif (flowtype.eq.8) then
     CYL = .true.
  elseif (flowtype.eq.9) then
     SHIT = .true.
  elseif (flowtype.eq.10) then
     TURB = .true.
  elseif (flowtype.eq.11) then
     LE = .true.
  elseif (flowtype.eq.12) then
     TE = .true.
  elseif (flowtype.eq.13) then 
     T3C = .true.
  elseif (flowtype.eq.0) then
     SRC = .true.
  endif
  read(30,*) ! Parameters for canonical flows (source, pulse, ...) [T:yes;F:no] [needs param_case.ini]
  read(30,*) is_param_case
  read(30,*) !-----------------------------------------------------------------------------------------------
  read(30,*) ! Reference quantities
  read(30,*) !-----------------------------------------------------------------------------------------------
  read(30,*) ! stagnation (T) or static (F) freestream conditions **** NOT IMPLEMENTED yet **** only F
  read(30,*) is_stagnation
  read(30,*) ! Reference temperature (could be freestream, wall, ...)
  read(30,*) T_ref
  read(30,*) ! Reference density (freestream, bulk, ...) (used only for static conditions, easier than p for dense gases)
  read(30,*) rho_ref
  read(30,*) ! Reference pressure (used only for stagnation conditions)
  read(30,*) p_ref
  read(30,*) ! Reference Mach number (freestream, bulk, turbulent ...)
  read(30,*) Mach
  read(30,*) ! Reference Reynolds number
  read(30,*) Re_ref
  read(30,*) ! Reference length [optional; deduced from reference Re unless stated otherwise]
  read(30,*) L_ref
  read(30,*) ! Reference velocity [optional; deduced from reference Mach unless stated otherwise]
  read(30,*) u_ref
  read(30,*) !-----------------------------------------------------------------------------------------------
  read(30,*) ! Boundary conditions [type of BCs are set in param_blocks.ini]
  read(30,*) !-----------------------------------------------------------------------------------------------
  read(30,*) ! Damping coefficient for sponge zone [defined in param_blocks.ini]
  read(30,*) alpha_sz
  read(30,*) ! Wall BC type [1:all variables imposed (dpdn=0);
  read(30,*) !               2:rho advanced on walls (compatible with slip wall)]
  read(30,*) type_wall_bc
  if (type_wall_bc.eq.2) is_wall2=.true.
  read(30,*) ! Adiabatic wall [T:adiabatic; F:isothermal], T_wall (isothermal, if not specified set to T_ref)
  read(30,*) is_adiab, T_wall
  read(30,*) ! Non-reflective characteristic BC: Relaxation coefficients **** NOT IMPLEMENTED yet ****
  read(30,*) !
  read(30,*) ! Tam & Dong BC: coordinates of radiation center [xcr ycr zcr]
  read(30,*) xcr, ycr, zcr
  read(30,*) ! Back-pressure outflow BC: Radial Equilibrium [T:yes;F:no]; ref. location ('min';'max';'mean')
  read(30,*) is_rea, loc_rea
  read(30,*) !-----------------------------------------------------------------------------------------------
  read(30,*) ! Flow initialization & forcing
  read(30,*) !-----------------------------------------------------------------------------------------------
  read(30,*) ! Inlet velocity vector (normalized by ref. velocity, if 0. 0. 0. then flow angles is used)
  read(30,*) Uc1, Uc2, Uc3
  read(30,*) ! Flow angles: theta_ref & phi_ref [degrees]
  read(30,*) theta_ref, phi_ref
  if (((Uc1**2+Uc2**2+Uc3**2).ne.1.0_wp).and.((Uc1**2+Uc2**2+Uc3**2).ne.0.0_wp)) call mpistop('Problem with inlet velocity vector specified',0)
  if ((Uc1**2+Uc2**2+Uc3**2).eq.0.0_wp) then
     theta_ref=theta_ref*pi/180.0_wp; phi_ref=phi_ref*pi/180.0_wp
     Uc1 = 1.0_wp/sqrt(cos(theta_ref)**2-sin(phi_ref)**2)
     Uc2 = sin(theta_ref)
     Uc3 = sin(phi_ref)
  else
     theta_ref=asin(Uc2)
       phi_ref=asin(Uc3)
  endif
  read(30,*) ! For back-pressure outflow BC: exit static pressure [Pa] (p_exit) [if 0. then p_ref is used]
  read(30,*) p_exit
  read(30,*) ! laminar boundary layer (LBL): is_LBL, is_similarity, Re_inlet, jdel
  read(30,*) is_LBL, is_similarity, Re_inlet, jdel
  read(30,*) ! is_forcing_bulk (to enforce mass flow rate for channel flow/periodic hill)
  read(30,*) is_forcing_bulk
  read(30,*) ! is_eigenmode (to enter eigenmodes at inlet)
  read(30,*) ! [needs param_stab.ini to define eigenmode parameters]
  read(30,*) is_eigenmode
  read(30,*) ! is_RFM (Random Fourier Modes) [needs param_RFM.ini to define RFM parameters]
  read(30,*) !  --> enter RFM at inlet OR initialize a field with RFM
  read(30,*) !  --> for freestream turbulence or turbulent boundary layers
  read(30,*) is_RFM
  read(30,*) ! is_suction_blowing [T:yes;F:no] [needs param_??.ini to define suction & blowing parameters]
  read(30,*) is_suction_blowing
  if (is_suction_blowing) call mpistop("Suction & blowing needs to be implemented",0)
  read(30,*) !-----------------------------------------------------------------------------------------------
  read(30,*) ! Geometry parameters
  read(30,*) !-----------------------------------------------------------------------------------------------
  read(30,*) ! Scaling value for the grid Lgrid (coordinates multiplied by Lgrid, if 0.0 L_ref is used)
  read(30,*) Lgrid
  read(30,*) ! Grid size for the spanwise extrusion if 3D (cartesian & 2D curvilinear solver): deltaz
  read(30,*) deltaz
  if ((nz.gt.1).and.(.not.is_curv3).and.(deltaz.le.0.0_wp)) &
      call mpistop("A grid size deltaz must be prescribed in param.ini for the spanwise extrusion",0)
  read(30,*) ! Stretching for the spanwise extrusion: nstrechz,  nrz1,nrz2,rsz (X nstrechi)
  read(30,*) ! ~> 1 stretching: begin index (nrz1) and stop index (nrz2), along with stretching coeff (rsz)
  read(30,*) ! ~> [nrz1, nrz2, rsz] replicated as the number of stretching (nstrechz)
  ! read(30, '(A)') test_str
  ! print *,test_str
  ! call mpistop('',0)
  read(30,*) nstretchz
  ! If stretching specified in the z direction
  if (nstretchz.gt..0) then
     rewind(30)
     do i=1,225 !number of lines before this one
        read(30,*) !
     enddo
     ! Allocation of temporary array
     allocate(temp_array(3,nstretchz))
     ! Read line
     read(30,*) nstretchz, temp_array
     ! Rearrange array with new ones
     allocate(nrz1(nstretchz),nrz2(nstretchz),rsz(nstretchz))
     do i=1,nstretchz
        nrz1(i)=temp_array(1,i); nrz2(i)=temp_array(2,i); rsz(i)=temp_array(3,i)
     enddo
     deallocate(temp_array)
     ! Check
     do i=1,nstretchz
        if (nrz2(i).le.nrz1(i)) call mpistop("Streching number "//trim(numchar(i))//" in z: nrz2 must be greater than nrz1",0)
        if (rsz(i).le.0.0_wp) call mpistop("Streching number "//trim(numchar(i))//" in z: stretching ratio rsz must be > 0",0)
        if (rsz(i).gt.2.0_wp) call mpistop("Streching number "//trim(numchar(i))//" in z: stretching ratio rsz seems relatively elevated...",0)
        if (i.gt.1) then
           if (nrz1(i).lt.nrz2(i-1)) call mpistop("Beginning index of streching num."//trim(numchar(i))//" ("//trim(numchar(nrz1(i)))//&
                                                  ")  must not be lower than end index of streching num."//trim(numchar(i-1))//" ("//trim(numchar(nrz2(i-1)))//")",0)
        endif
     enddo
  endif
  read(30,*) ! Translation vector for planar periodicity: Lxp, Lyp, Lzp (if Lzp=0.0, directly determined)
  read(30,*) Lxp, Lyp, Lzp
  read(30,*) ! Angle for angular periodicity: theta_period [degrees]
  read(30,*) theta_period
  read(30,*) ! Parameters for grid: creation/reading of grid (.x) & writting in .bin
  read(30,*) ! user-defined grid [T:yes;F:no] [needs param_grid.ini]
  read(30,*) is_def_grid
  read(30,*) ! If not user-defined and directly prescribed, needs .x grid files if idepart=1:
  read(30,*) ! Directory for grid files: dirGRID
  read(30,*) dirGRID
  read(30,*) ! Name for grid files: nameGRID
  read(30,*) nameGRID
  read(30,*) !===============================================================================================
  read(30,*) !                                  Pre/post-processing options
  read(30,*) !===============================================================================================
  read(30,*) ! [needs param_pp.ini to define post-processing parameters]
  read(30,*) ! Directory (for post-processing mode)
  read(30,*) dirDATA
  read(30,*) ! Pre-processing for grid:
  read(30,*) ! Half-cell suppression (if not adjoint blocks): is_half_cell
  read(30,*) is_half_cell
  read(30,*) ! Coarse grid on half the points: is_coarse_grid
  read(30,*) is_coarse_grid
  read(30,*) ! create extended grid (boolean, only for full 3D) (obsolete now ?)
  read(30,*) is_read_ex
  read(30,*) ! Add stretching zone for exit blocks: is_add_sponge
  read(30,*) is_add_sponge
  read(30,*) ! Other pre-processing modes:
  read(30,*) ! Saturation curve: is_satur_curve
  read(30,*) is_satur_curve
  read(30,*) ! Linear Stability solver: is_LST
  read(30,*) is_LST
  close(30)

  ! Set number of ghost points
  ! --------------------------
  ! Euler fluxes
  ngh  =(stencil-1)/2   ! (11pts-stencil)
  ! Viscous fluxes
  ngh_v=iorder_visc/2
  if (iorder_visc.eq.0) ngh_v=2 ! for Euler calculations keep 5-pt stencil USEFUL?? TO BE CHANGED
  ! RANS fluxes
  ngh_r=(stencil_RANS-1)/2   ! (5/3pts-stencil)
  ngh_v_r=ngh_r

  ! Read additional advanced setting files
  ! --------------------------------------
  if (is_param_case) call read_param_case
  if (is_RANS_adv) call read_param_RANS

  ! Extra setting (obsolete or very advanced)
  ! -------------
  ! Type of parallelisation for IRS (1: ngh, 2: scalapack band/tri, 3 : pascal)
  ! For the moment, only mode 1 is working
  type_para=1

  ! for what ????????????? TO BE CHANGED
  dtprint = 10000.0_wp
  
  ! Use of one-sided or two sided communications
  ! **** BOOLEAN NOT USED yet ****
  is_two_sided_comm = .false.; is_comm_onesided = .false.

  ! Numerical dissipation added to increment before IRS smoothing
  ! is_dissip_in_increments=.true.

  ! Old way of defining grid: to be let .false. by default !
  is_grid_old=.true.

  if (is_residue) then
     open(31,file='residuals.bin',form='unformatted',status='unknown')
     rewind(31)
  endif
  
end subroutine read_param


!===============================================================================
subroutine read_param_RANS
!===============================================================================
  !> author: (To be written)
  !> date: (To be written)
  !> Reading advanced setting for RANS
!===============================================================================
    use mod_constant
    use mod_rans
  use warnstop
  implicit none
  ! ---------------------------------------------------------------------------
  logical :: iexist
  ! ---------------------------------------------------------------------------

  ! Initializaiton of the options:
  is_transition=.false.
   ! Read param.ini file
  ! ===================
  inquire(file=trim("param_rans.ini"), exist=iexist)
  if (.not.iexist) then
     call mpistop('"param_rans.ini" does not exist!', 0)
  endif

  open(30,file=trim("param_rans.ini"))
  rewind(30)
  read(30,*) !===============================================================================================
  read(30,*) !===============================================================================================
  read(30,*) ! MUSICA2 : Read advanced RANS settings in param_RANS.ini
  read(30,*) !===============================================================================================
  read(30,*) !===============================================================================================
  read(30,*) ! 1)   Specify transition parameters
  read(30,*) !      - Algebraic or transport model
  read(30,*) !      - Inlet conditions
  read(30,*) !===============================================================================================
  read(30,*) !===============================================================================================
  read(30,*) !              1) Transitional Modelling: 
  read(30,*) !===============================================================================================
  read(30,*) ! Is it transitional simulation or not, T:yes, F:no
  read(30,*) is_transition
  read(30,*) ! What model ( "ALG": algebraic model, "GRE": Gamma-Re_theta transport model )
  read(30,*) model_transition
  model_transition = trim(model_transition)
  read(30,*) ! Inlet turbulent intensity:
  read(30,*) tu_inlet
  tu_inlet=tu_inlet/100.0_wp
end subroutine read_param_RANS


!===============================================================================
subroutine read_param_case
!===============================================================================
  !> author: (To be written)
  !> date: (To be written)
  !> Reading case settings, for canonical flows and/or basic geometries
!===============================================================================
  use warnstop
  use mod_constant
  use mod_init_flow
  use mod_routines
  implicit none
  ! ---------------------------------------------------------------------------
  logical :: iexist
  ! ---------------------------------------------------------------------------

  ! Read param_case.ini file
  ! ========================
  inquire(file="param_case.ini", exist=iexist)
  if (.not.iexist) then
     call mpistop('param_case file for canonical case does not exist!', 0)
  endif
  ! ---------------
  ! TO BE COMPLETED
  ! ---------------
  open(30,file="param_case.ini")
  rewind(30)
  read(30,*) !===============================================================================================
  read(30,*) !===============================================================================================
  read(30,*) ! MUSICA2 : fill advanced case parameters for canonical flow
  read(30,*) !===============================================================================================
  read(30,*) !===============================================================================================
  read(30,*) ! ***** WORK IN PROGRESS *****  TO BE DEVELOPED
  read(30,*) ! ----------------------------------------------------------------------------------------------
  read(30,*) ! Flow initialization: generate canonical flow initialization
  read(30,*) ! ----------------------------------------------------------------------------------------------
  read(30,*) !   -> Accoustic source
  read(30,*) !   -> Accoustic pulse
  read(30,*) !   -> Vortex
  read(30,*) !   -> And others to be implemented...
  read(30,*) !===============================================================================================
  read(30,*) !===============================================================================================
  read(30,*) !                                      Flow initialization
  read(30,*) !===============================================================================================
  read(30,*) ! source definition: ampl(Pa), lambda/deltax, half-width b/deltax, x_src, y_src
  read(30,*) is_src, ampl_src, omeg_src, b_src, x_src, y_src
  read(30,*) ! pulse definition: ampl(Pa), half-width b/deltax, x_pulse, y_pulse, z_pulse
  read(30,*) is_pulse, ampl_pulse, b_pulse, x_pulse, y_pulse, z_pulse
  read(30,*) ! vortex definition: type_vortex, x_vortex, y_vortex, z_vortex
  read(30,*) is_vortex, type_vortex, x_vortex, y_vortex, z_vortex
  read(30,*) ! is_supersonic **** NOT IMPLEMENTED yet ****
  read(30,*) is_supersonic
  close(30)

end subroutine read_param_case


