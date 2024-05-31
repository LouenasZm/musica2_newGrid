!===============================================================================
submodule (mod_RFM) smod_RFM_init
!===============================================================================
  !> author: AB & XG
  !> date: February 2024
  !> Generation of Random Fourier Modes (RFM)
  !> Initialization routines
!=============================================================================== 

contains

  !===============================================================================
  module subroutine init_RFM
  !===============================================================================
    !> Initialization of RFM modes
  !===============================================================================
    use mod_database  ! <- module for TBL/CHAN databases [Initial_condition]
    use mod_interp1   ! <- module for 1D interpolation [Mathematics]
    use mod_utils
    use mod_tranprop
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: k,j,jbc
    ! ----------------------------------------------------------------------------
    ! random number generator
    !real(wp) :: randomnb
    !integer nseed
    !integer, dimension(:), allocatable :: initseed
    ! ----------------------------------------------------------------------------
    ! turbulent spectrum
    real(wp) :: uf
    real(wp) :: xke,xeps
    real(wp) :: dkl,xkk,norm
    real(wp), dimension(:), allocatable :: xkn,dxkn,vkm,utilde
    ! ----------------------------------------------------------------------------
    ! divergence-free stochastic field
    real(wp) :: var
    real(wp) :: alpha,calpha,salpha,phi,cphi,sphi,ctheta,stheta
    ! real(wp), dimension(:), allocatable :: um_m
    ! ----------------------------------------------------------------------------
    ! anisotropy transformation
    real(wp) :: Ldim_a
    real(wp), dimension(:), allocatable :: uu_rms,vv_rms,ww_rms,uv_rms,uw_rms,vw_rms
    real(wp), dimension(3,3) :: tau
    ! for LAPACK's diagonalization
    integer :: info
    integer, parameter :: ldvl=1,ldvr=3,lwork=4*3
    real(wp), dimension(lwork) :: work
    real(wp), dimension(3) :: wi
    real(wp), dimension(ldvl,3) :: vl
    ! ----------------------------------------------------------------------------
    logical :: iexist
    ! Wall coordinates
    real(wp) :: yw_BCj_proc
    ! Damping if inlet injection zone restricted
    integer :: cinj_jmin,cinj_jmax
    real(wp) :: xinj_jmin,xinj_jmax,xinj_jmin_proc,xinj_jmax_proc
    real(wp) :: yinj_jmin,yinj_jmax,yinj_jmin_proc,yinj_jmax_proc
    real(wp), dimension(1:nz) :: xinj3_jmin,yinj3_jmin,zinj3_jmin
    real(wp), dimension(1:nz) :: xinj3_jmin_proc,yinj3_jmin_proc,zinj3_jmin_proc
    real(wp), dimension(1:nz) :: xinj3_jmax,yinj3_jmax,zinj3_jmax
    real(wp), dimension(1:nz) :: xinj3_jmax_proc,yinj3_jmax_proc,zinj3_jmax_proc
    real(wp), dimension(:,:), allocatable :: damping_coeff_all
    ! TBL
    real(wp), dimension(2) :: d99_TBL_proc
    ! To determine mean convection velocity
    logical :: from_BCref
    integer :: ncount_,jav1,jav2
    real(wp) :: rho_m_proc,rhou_m_proc,rhov_m_proc,rhow_m_proc
    real(wp) :: rho_m_block,rhou_m_block,rhov_m_block,rhow_m_block
    real(wp) :: rho_m,rhou_m,rhov_m,rhow_m
    ! To determine reference wall quantities
    real(wp) :: mu_wm_proc,rho_wm_proc,prs_wm_proc,rho_wm,prs_wm,T_wm,mu_wm
    real(wp) :: utau_m
    ! Temporary, to be changed
    real(wp) :: norm_,norm_proc,xg_in_

    if (is_init1_RFM) then
       
       ! Read RFM mode parameters
       ! ------------------------
       call read_param_RFM

       ! Determine origin of inlet plane <~ for init_vel_RFM
       ! -------------------------------
       xg_in_=0.0_wp
       if (iproc.eq.0) then ! XG VALID ONLY IF proc 0 is in block inlet ??????????????
          if (is_curv3) then
             xg_in_= xc3(1,1,1)
          else if (is_curv) then
             xg_in_= xc(1,1)
          else
             xg_in_= x(1)
          endif
       endif
       ! to be replace by BCAST
       call MPI_BCAST(xg_in_,1,MPI_DOUBLE_PRECISION,0,COMM_global,info)
       !call MPI_ALLREDUCE(xg_in_,xg_in,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info) ! XG WHY MPI_SUM other are zero ??????????????

       ! Definition of a reference velocity of convection ~> stored in u_m, v_m, w_m
       ! ------------------------------------------------
       ncount_=0
       rho_m_proc=0.0_wp
       rhou_m_proc=0.0_wp
       rhov_m_proc=0.0_wp
       rhow_m_proc=0.0_wp
       rho_m_block=0.0_wp
       rhou_m_block=0.0_wp
       rhov_m_block=0.0_wp
       rhow_m_block=0.0_wp
       ! other initialization
       yw_BCj=0.0_wp
       d99_TBL=0.0_wp
       d99_TBL_proc=0.0_wp

       ! ------------------------
       ! Mean velocity convection
       ! ------------------------
       ! REX: mean field needs to be constant along the injection plane

       ! Determination of if convection mean field based on BC_ref
       ! ---------------------------------------------------------
       ! ~> if yes, u_m, v_m & w_m calculated based on reference values
       ! ~> if no, u_m, v_m & w_m calculated based on U_ref, theta_ref and phi_ref

       ! XG WHY is_RFM ?????? should be true at this stage
       call MPI_Allreduce(is_RFM.and.BC_face(1,1)%is_mean_ref,from_BCref,1,MPI_LOGICAL,MPI_LOR,COMM_global,info)

       ! Initialization of indicator is_RFM_block
       is_RFM_block=is_RFM
       ! is_RFM put to false if not at imin
       if (BC_face(1,1)%sort.ge.0) is_RFM=.false.

       if ((from_BCref).and.(idepart.ne.FROM_SCRATCH)) then
          ! Mean convection field based on BC_ref
          ! -------------------------------------
          ! Density-averaged
          ! Reference velocity vector at i_min where RFM injected

          ! Average on injection zone
          ! -------------------------
          jav1=0; jav2=0
          if (is_RFM) then
             jav1=1; jav2=ny
             if (ndy_RFM.gt.1)    jav1 = min(max(ndy_RFM-coord(2)*ny,1),ny)
             if (nfy_rfm.lt.ngy) jav2 = max(min(nfy_rfm-coord(2)*ny,ny),1)
          endif

          if (jav1.ne.jav2) then
             ! Averaged per proc
             ! -----------------
             ! rho, rhou, rhov & rhow averaged
              rho_m_proc = sum(BC_face(1,1)%Uref(jav1:jav2,1:nz,1))/((jav2+1-jav1)*nz)
             rhou_m_proc = sum(BC_face(1,1)%Uref(jav1:jav2,1:nz,1)*BC_face(1,1)%Uref(jav1:jav2,1:nz,2))/((jav2+1-jav1)*nz)
             rhov_m_proc = sum(BC_face(1,1)%Uref(jav1:jav2,1:nz,1)*BC_face(1,1)%Uref(jav1:jav2,1:nz,3))/((jav2+1-jav1)*nz)
             rhow_m_proc = sum(BC_face(1,1)%Uref(jav1:jav2,1:nz,1)*BC_face(1,1)%Uref(jav1:jav2,1:nz,4))/((jav2+1-jav1)*nz)

             ! Averaged per block
             ! ------------------
             ! number of processor involved in the average
             ncount_ = 0
             call MPI_ALLREDUCE(1,ncount_,1,MPI_INTEGER,MPI_SUM,COMM_intrablock,info)
             ! rho, rhou, rhov & rhow averaged
             call MPI_ALLREDUCE(rho_m_proc,rho_m_block,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_intrablock,info)
             call MPI_ALLREDUCE(rhou_m_proc,rhou_m_block,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_intrablock,info)
             call MPI_ALLREDUCE(rhov_m_proc,rhov_m_block,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_intrablock,info)
             call MPI_ALLREDUCE(rhow_m_proc,rhow_m_block,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_intrablock,info)
             rho_m_block=rho_m_block/ncount_; rhou_m_block=rhou_m_block/ncount_
             rhov_m_block=rhov_m_block/ncount_; rhow_m_block=rhow_m_block/ncount_

          else if (is_RFM_block) then
             ! Averaged per block
             ! ------------------
             ! All processors involved, necessary if iproc_leader
             ! not at BC with RFM (~> here, imin so it is not
             ! the case but it is just to have some genericity)
             ! ------------------
             ! number of processor involved in the average
             call MPI_ALLREDUCE(0,ncount_,1,MPI_INTEGER,MPI_SUM,COMM_intrablock,info)
             ! rho, rhou, rhov & rhow averaged
             call MPI_ALLREDUCE(rho_m_proc,rho_m_block,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_intrablock,info)
             call MPI_ALLREDUCE(rhou_m_proc,rhou_m_block,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_intrablock,info)
             call MPI_ALLREDUCE(rhov_m_proc,rhov_m_block,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_intrablock,info)
             call MPI_ALLREDUCE(rhow_m_proc,rhow_m_block,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_intrablock,info)
             rho_m_block=rho_m_block/ncount_; rhou_m_block=rhou_m_block/ncount_
             rhov_m_block=rhov_m_block/ncount_; rhow_m_block=rhow_m_block/ncount_
          endif

          ! Averaged accross all blocks
          ! ---------------------------
          ncount_=0; rho_m=0.0_wp; rhou_m=0.0_wp; rhov_m=0.0_wp; rhow_m=0.0_wp
          if (iproc.eq.iproc_leader(nob(iproc))) then
             ! number of blocks involved in the average
             if (is_RFM_block) then
                call MPI_ALLREDUCE(1,ncount_,1,MPI_INTEGER,MPI_SUM,COMM_interblock,info)
             else
                call MPI_ALLREDUCE(0,ncount_,1,MPI_INTEGER,MPI_SUM,COMM_interblock,info)
             endif

             ! rho, rhou, rhov & rhow averaged
             call MPI_ALLREDUCE(rho_m_block,rho_m,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_interblock,info)
             call MPI_ALLREDUCE(rhou_m_block,rhou_m,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_interblock,info)
             call MPI_ALLREDUCE(rhov_m_block,rhov_m,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_interblock,info)
             call MPI_ALLREDUCE(rhow_m_block,rhow_m,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_interblock,info)
             rho_m=rho_m/ncount_; rhou_m=rhou_m/ncount_; rhov_m=rhov_m/ncount_; rhow_m=rhow_m/ncount_
          endif

          ! Distribution accross all processors
          ! -----------------------------------
          if (iproc.eq.0) then
             ! Temporary stored in rho_m_proc, ...
             rho_m_proc=rho_m; rhou_m_proc=rhou_m; rhov_m_proc=rhov_m; rhow_m_proc=rhow_m
             ! Distribution
             call MPI_ALLREDUCE(rho_m_proc,rho_m,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
             call MPI_ALLREDUCE(rhou_m_proc,rhou_m,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
             call MPI_ALLREDUCE(rhov_m_proc,rhov_m,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
             call MPI_ALLREDUCE(rhow_m_proc,rhow_m,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
          else
             call MPI_ALLREDUCE(0.0_wp,rho_m,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
             call MPI_ALLREDUCE(0.0_wp,rhou_m,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
             call MPI_ALLREDUCE(0.0_wp,rhov_m,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
             call MPI_ALLREDUCE(0.0_wp,rhow_m,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
          endif

          ! Computation of u_m, v_m & w_m
          u_m = rhou_m/rho_m
          v_m = rhov_m/rho_m
          w_m = rhow_m/rho_m
       else
          ! Mean convection field based on U_ref, theta_ref & phi_ref (specified in param.ini)
          ! ---------------------------------------------------------
          u_m = U_ref*sqrt(cos(theta_ref*pi/180.0_wp)**2-sin(phi_ref*pi/180.0_wp)**2)
          v_m = U_ref*sin(theta_ref*pi/180.0_wp)
          w_m = U_ref*sin(phi_ref*pi/180.0_wp)
          rho_m = rho_ref
       endif

       ! Mean velocity convection in the TBL
       ! -----------------------------------
       ! REX: needs to be constant accross the boundary layer
       ! Arbitrary choice ~> 0.6 of reference convection velocity
       if ((is_RFM_TBL(1)).or.(is_RFM_TBL(2))) then
          u_m_TBL = u_m*0.6; v_m_TBL = v_m*0.6; w_m_TBL = w_m*0.6
       endif

       ! Amplitude of fluctuations
       ! -------------------------
       ampl_RFM = Tu_RFM*u_m ! ampl_RFM = Tu*U_ref [m/s]
       ampl_TBL = Tu_TBL*u_m

       ! Determination of the list of blocks using RFM & the new proc leader for RFM display
       ! -----------------------------------------------------------------------------------
       allocate(is_RFM_blocks(nbloc))
       iproc_leader_rfm=-1; ncount_=0
       if (iproc.eq.iproc_leader(nob(iproc))) then
          call MPI_ALLGATHER(is_RFM_block,1,MPI_LOGICAL,is_RFM_blocks,1,MPI_LOGICAL,COMM_interblock,info)

          do while ((iproc_leader_rfm.eq.-1).and.(ncount_.lt.nbloc))
             ncount_=ncount_+1
             if (is_RFM_blocks(ncount_)) iproc_leader_rfm=iproc_leader(ncount_)
          enddo
          ! Temporary stored in ncount_
          ncount_ = iproc_leader_rfm
       endif
       ! Communication on all procs of each block
       call MPI_ALLREDUCE(ncount_,iproc_leader_rfm,1,MPI_INTEGER,MPI_SUM,COMM_intrablock,info)

       ! ----------------------------------------------
       ! End of subroutine if is_RFM_block put to false
       ! ----------------------------------------------
       if (.not.is_RFM_block) then
          is_init1_RFM=.false.
          return
       endif

       ! ----------------------------------------
       !      II/ Turbulent Boundary Layer
       ! ----------------------------------------
       ! ~> Reading database TBL profile target
       ! ~> if is_damping (injection zone):
       !    _ Determination of indices for FST
       !      injection based on TBL thickness
       ! ~> Discretization of wavenumbers
       !  & Random number generation
       !  & Computation of spectrum
       !  & Anisotropic transformation of modes
       !    _ Spectra computed based on TBL thickness

       ! if processor at imin, all steps
       if (is_boundary(1,1)) then

          ! Loop on the number of TBL to inject (jmin and/or jmax)
          ! -----------------------------------
          do jbc=1,2
             if (is_RFM_TBL(jbc)) then
                
                ! Read Reynolds stresses
                ! ----------------------
                call read_database(base_TBL,Re_prof_TBL(jbc))

                ! Determination of reference value at the wall
                ! --------------------------------------------
                if (is_bc_wall(2,jbc)) then
                   ! If BC_ref, averaged at wall in each block
                   if ((BC_face(1,1)%is_mean_ref).and.(idepart.ne.FROM_SCRATCH)) then
                      ! If BC_ref, communication of prs and rho at the wall
                      if (jbc.eq.1) then
                         rho_wm_proc = sum(BC_face(1,1)%Uref(1,1:nz,1))/nz
                         prs_wm_proc = sum(BC_face(1,1)%Uref(1,1:nz,5))/nz
                      else
                         rho_wm_proc = sum(BC_face(1,1)%Uref(ny,1:nz,1))/nz
                         prs_wm_proc = sum(BC_face(1,1)%Uref(ny,1:nz,5))/nz
                      endif
                      call MPI_ALLREDUCE(rho_wm_proc,rho_wm,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXY,info)
                      rho_wm=rho_wm/ndomz
                      call MPI_ALLREDUCE(prs_wm_proc,prs_wm,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXY,info)
                      prs_wm=prs_wm/ndomz
                      ! Determination of T and then mu based on averaged wall values
                      T_wm = tcalc_pro(prs_wm,rho_wm,T_ref)
                      mu_wm = viscosity_law(T_wm,rho_wm)
                   ! Else, based on general reference value
                   else
                      mu_wm=muw_ref
                      rho_wm=rho_ref
                   endif

                   ! Communication to other procs of the BC
                   mu_wm_proc=mu_wm; rho_wm_proc=rho_wm
                   call MPI_ALLREDUCE(mu_wm_proc,mu_wm,1,MPI_DOUBLE_PRECISION,MPI_MAX,COMMX,info)
                   call MPI_ALLREDUCE(rho_wm_proc,rho_wm,1,MPI_DOUBLE_PRECISION,MPI_MAX,COMMX,info)
                else
                   ! Reception from processes of the BC at wall
                   mu_wm_proc=0.0_wp; rho_wm_proc=0.0_wp
                   call MPI_ALLREDUCE(mu_wm_proc,mu_wm,1,MPI_DOUBLE_PRECISION,MPI_MAX,COMMX,info)
                   call MPI_ALLREDUCE(rho_wm_proc,rho_wm,1,MPI_DOUBLE_PRECISION,MPI_MAX,COMMX,info)
                endif

                ! Determination of BL thickness
                ! -----------------------------
                utau_m=(u_m**2+v_m**2)**0.5/uinf_utau_db
                d99_TBL(jbc)=Re_tau_db*mu_wm/(rho_wm*utau_m)

                ! ! TEMPORARY /!\
                ! d99_TBL(1) = d99_TBL(1)/6

                ! ! If FST injected, needs to determine injection zone
                ! ! --------------------------------------------------
                ! if (is_RFM_FST) then
                   ! Damping injection height
                   ! ------------------------
                   h_damp(1)=0.1*d99_TBL(jbc)

                   ! Scatter wall coordinates (from coord(2)==0/ndomj-1 process to all processes)
                   ! ------------------------
                   if (is_curv3) then
                      call mpistop('Not tested in mod_RFM.f90',1)
                      ! xinj3_jmin=0.0_wp; xinj3_jmax=0.0_wp; xinj3_jmin_proc = 0.0_wp; xinj3_jmax_proc = 0.0_wp
                      ! yinj3_jmin=0.0_wp; yinj3_jmax=0.0_wp; yinj3_jmin_proc = 0.0_wp; yinj3_jmax_proc = 0.0_wp
                      ! zinj3_jmin=0.0_wp; zinj3_jmax=0.0_wp; zinj3_jmin_proc = 0.0_wp; zinj3_jmax_proc = 0.0_wp
                      ! ! Location at jmin
                      ! if (coord(2)==cinj_jmin) then
                      !   xinj3_jmin_proc(1:nz)=xc3(1,ndy_RFM-coord(2)*ny,1:nz)
                      !   yinj3_jmin_proc(1:nz)=yc3(1,ndy_RFM-coord(2)*ny,1:nz)
                      !   zinj3_jmin_proc(1:nz)=zc3(1,ndy_RFM-coord(2)*ny,1:nz)
                      ! endif
                      ! call MPI_ALLREDUCE(xinj3_jmin_proc(1:nz),xinj3_jmin(1:nz),nz,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXZ,info)
                      ! call MPI_ALLREDUCE(yinj3_jmin_proc(1:nz),yinj3_jmin(1:nz),nz,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXZ,info)
                      ! call MPI_ALLREDUCE(zinj3_jmin_proc(1:nz),zinj3_jmin(1:nz),nz,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXZ,info)
                      ! ! Location at jmax
                      ! if (coord(2)==cinj_jmax) then
                      !   xinj3_jmax_proc(1:nz)=xc3(1,nfy_rfm-coord(2)*ny,1:nz)
                      !   yinj3_jmax_proc(1:nz)=yc3(1,nfy_rfm-coord(2)*ny,1:nz)
                      !   zinj3_jmax_proc(1:nz)=zc3(1,nfy_rfm-coord(2)*ny,1:nz)
                      ! endif
                      ! call MPI_ALLREDUCE(xinj3_jmax_proc(1:nz),xinj3_jmax(1:nz),nz,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXZ,info)
                      ! call MPI_ALLREDUCE(yinj3_jmax_proc(1:nz),yinj3_jmax(1:nz),nz,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXZ,info)
                      ! call MPI_ALLREDUCE(zinj3_jmax_proc(1:nz),zinj3_jmax(1:nz),nz,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXZ,info)
                   else if (is_curv) then
                      call mpistop('Not tested in mod_RFM.f90',1)
                      ! xinj_jmin=0.0_wp; xinj_jmax=0.0_wp; xinj_jmin_proc = 0.0_wp; xinj_jmax_proc = 0.0_wp
                      ! yinj_jmin=0.0_wp; yinj_jmax=0.0_wp; yinj_jmin_proc = 0.0_wp; yinj_jmax_proc = 0.0_wp
                      ! ! Location at jmin
                      ! if (coord(2)==cinj_jmin) then
                      !   xinj_jmin_proc=xc(1,ndy_RFM-coord(2)*ny); yinj_jmin_proc=yc(1,ndy_RFM-coord(2)*ny)
                      ! endif
                      ! call MPI_ALLREDUCE(xinj_jmin_proc,xinj_jmin,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXZ,info)
                      ! call MPI_ALLREDUCE(yinj_jmin_proc,yinj_jmin,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXZ,info)
                      ! ! Location at jmax
                      ! if (coord(2)==cinj_jmax) then
                      !   xinj_jmax_proc=xc(1,nfy_rfm-coord(2)*ny); yinj_jmax_proc=yc(1,nfy_rfm-coord(2)*ny)
                      ! endif
                      ! call MPI_ALLREDUCE(xinj_jmax_proc,xinj_jmax,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXZ,info)
                      ! call MPI_ALLREDUCE(yinj_jmax_proc,yinj_jmax,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXZ,info)
                   else
                      yw_BCj(jbc)=0.0_wp; yw_BCj_proc = 0.0_wp
                      ! Location at jmin
                      if ((coord(2)==0).and.(jbc.eq.1)) then
                         yw_BCj_proc=y(1)
                      else if ((coord(2)==ndomy-1).and.(jbc.eq.2)) then
                         yw_BCj_proc=y(ny)
                      endif
                      call MPI_ALLREDUCE(yw_BCj_proc,yw_BCj(jbc),1,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXZ,info)
                   endif

                   ! Determination of injection position (ndy_RFM,nfy_rm)
                   ! -----------------------------------
                   ! Taken above TBL at 1.2*d99
                   if (is_curv3) then
                      call mpistop('Not implemented in mod_RFM.f90',1)
                   else if (is_curv) then
                      call mpistop('Not implemented in mod_RFM.f90',1)
                   else
                      ! Evaluation for jmin
                      if (jbc.eq.1) then
                         ! Variation of y along j increasing
                         if (y(2)-y(1).gt.0.0_wp) then
                            ! 1.2*d99 TBL ends after the processor
                            if ((y(ny)-yw_BCj(jbc)).lt.1.2*d99_TBL(jbc)) then
                               ndy_RFM=0
                            ! 1.2*d99 TBL ends before the processor
                            elseif ((y(1)-yw_BCj(jbc)).gt.1.2*d99_TBL(jbc)) then
                               ndy_RFM=0
                            ! 1.2*d99 TBL in the proc domain
                            else
                               j=1
                               do while (y(j)-yw_BCj(jbc).lt.1.2*d99_TBL(jbc))
                                  ndy_RFM=j
                                  j=j+1
                               enddo
                               ndy_RFM=ndy_RFM+coord(2)*ny
                            endif
                         ! Variation of y along j decreasing
                         else
                            ! 1.2*d99 TBL ends after the processor
                            if ((y(1)-yw_BCj(jbc)).lt.1.2*d99_TBL(jbc)) then
                               ndy_RFM=0
                            ! 1.2*d99 TBL ends before the processor
                            else if ((y(ny)-yw_BCj(jbc)).gt.1.2*d99_TBL(jbc)) then
                               ndy_RFM=0
                            ! 1.2*d99 TBL in the proc domain
                            else
                               j=ny
                               do while (y(j)-yw_BCj(jbc).lt.1.2*d99_TBL(jbc))
                                  ndy_RFM=j
                                  j=j-1
                               enddo
                               ndy_RFM=ndy_RFM+coord(2)*ny
                            endif
                         endif
                      ! Evaluation for jmax
                      else
                         ! Variation of y along j increasing
                         if (y(2)-y(1).lt.0.0_wp) then
                            ! 1.2*d99 TBL ends after the processor
                            if ((y(ny)-yw_BCj(jbc)).lt.1.2*d99_TBL(jbc)) then
                               nfy_rfm=0
                            ! 1.2*d99 TBL ends before the processor
                            else if ((y(1)-yw_BCj(jbc)).gt.1.2*d99_TBL(jbc)) then
                               nfy_rfm=0
                            ! 1.2*d99 TBL in the proc domain
                            else
                               j=1
                               do while (y(j)-yw_BCj(jbc).lt.1.2*d99_TBL(jbc))
                                  nfy_rfm=j
                                  j=j+1
                               enddo
                               nfy_rfm=nfy_rfm+coord(2)*ny
                            endif
                         ! Variation of y along j decreasing
                         else
                            ! 1.2*d99 TBL ends after the processor
                            if ((y(1)-yw_BCj(jbc)).lt.1.2*d99_TBL(jbc)) then
                               nfy_rfm=0
                            ! 1.2*d99 TBL ends before the processor
                            else if ((y(ny)-yw_BCj(jbc)).gt.1.2*d99_TBL(jbc)) then
                               nfy_rfm=0
                            ! 1.2*d99 TBL in the proc domain
                            else
                               j=ny
                               do while (y(j)-yw_BCj(jbc).lt.1.2*d99_TBL(jbc))
                                  nfy_rfm=j
                                  j=j-1
                               enddo
                               nfy_rfm=nfy_rfm+coord(2)*ny
                            endif
                         endif
                      endif
                   endif

                   ! Communication of injection position
                   if (jbc.eq.1) then
                      nrfm_proc = ndy_RFM
                      call MPI_ALLREDUCE(nrfm_proc,ndy_RFM,1,MPI_INTEGER,MPI_MAX,COMMXZ,info)
                   else
                      nrfm_proc = nfy_rfm
                      call MPI_ALLREDUCE(nrfm_proc,nfy_rfm,1,MPI_INTEGER,MPI_MAX,COMMXZ,info)
                   endif

                ! endif

                ! Free memory allocated for database
                deallocate(y_db,um_db,uu_db,vv_db,ww_db,uv_db,uw_db,vw_db)
             endif
          enddo
       endif

       ! Communication of ndy_RFM & nfy_rfm
       ! ----------------------------------
       nrfm_proc = ndy_RFM
       call MPI_ALLREDUCE(nrfm_proc,ndy_RFM,1,MPI_INTEGER,MPI_MAX,COMM_intrablock,info)
       nrfm_proc = nfy_RFM
       call MPI_ALLREDUCE(nrfm_proc,nfy_RFM,1,MPI_INTEGER,MPI_MAX,COMM_intrablock,info)

       ! Communication of d99_TBL to all other processors
       !   (necessary for BL initialization at the block)
       d99_TBL_proc = d99_TBL
       call MPI_ALLREDUCE(d99_TBL_proc,d99_TBL,2,MPI_DOUBLE_PRECISION,MPI_MAX,COMM_intrablock,info)

       ! End of first initialization of RFM
       ! ----------------------------------
       is_init1_RFM = .false.
       if (idepart.eq.FROM_SCRATCH) return
    endif

    ! Determination of processors really involved in generation of TBL
    ! ----------------------------------------------------------------
    ! Depends on injection indices: ndy_RFM & nfy_RFM
    if ((is_RFM_TBL(1)).and.(1+coord(2)*ny.ge.ndy_RFM))   is_RFM_TBL(1)=.false.
    if ((is_RFM_TBL(2)).and.((coord(2)+1)*ny.le.nfy_RFM)) is_RFM_TBL(2)=.false.

    ! Generation of random Fourier modes for TBL for these processors
    ! ---------------------------------------------------------------
    ! Only processors with TBL to inject generate the random Fourier modes
    ! Unfrozen turbulence evolution not applied XG: WHY ??? needed for time derivative
    
    if ((is_RFM_TBL(1)).or.(is_RFM_TBL(2))) then
       allocate(sigma1s_TBL(Nmode,ny,2),sigma2s_TBL(Nmode,ny,2),sigma3s_TBL(Nmode,ny,2))
       allocate(xk1s_TBL(Nmode,ny,2),xk2s_TBL(Nmode,ny,2),xk3s_TBL(Nmode,ny,2))
       allocate(psi_TBL(Nmode,2),omn_TBL(Nmode,2))
       allocate(vr_TBL(ldvr,3,ny,2))

       ! Loop on the number of TBL to inject (jmin and/or jmax)
       ! -----------------------------------
       do jbc=1,2
          ! if ((Re_prof_TBL(jbc).gt.0).and.(((jbc.eq.1).and.((coord(2)+1)*ny.le.ndy_RFM)) .or. &
          !                                  ((jbc.eq.2).and.(1+coord(2)*ny.ge.nfy_rfm)))) then
          if (is_RFM_TBL(jbc)) then
             ! Memory allocation
             ! -----------------
             allocate(xkn(Nmode),dxkn(Nmode),vkm(Nmode),utilde(Nmode))

             ! Read Reynolds stresses
             ! ----------------------
             call read_database(base_TBL,Re_prof_TBL(jbc))

             ! Calculation of kmin, kmax and kmn
             ! ---------------------------------
             Lf_tbl = 0.3*d99_TBL(jbc)
             delta = Lf_tbl/(sqrt(12.0_wp/5.0_wp)*0.74684_wp)
             xkmx = 1/delta
             xkmax = 2.0_wp*pi/Lf_tbl*7
             xkmin = 2.0_wp*pi/(Lf_tbl*9)

             ! linear or logarithmic wavenumber discretization
             ! -----------------------------------------------
             if (kdist=='log') then
                dkl=(log(xkmax)-log(xkmin))/dble(Nmode-1)
                do k=1,Nmode
                   xkn(k)=exp(log(xkmin)+dble(k-1)*dkl)
                enddo
             else if (kdist=='lin') then
                dkl=(xkmax-xkmin)/dble(Nmode-1)
                do k=1,Nmode
                   xkn(k) = xkmin+dble(k-1)*dkl
                enddo
             else
                call mpistop('Wrong choice for the wavenumber discretization in param_RFM !',0)
             endif

             ! Wavenumber increments
             ! ---------------------
             do k=2,Nmode-1
                dxkn(k)=0.5_wp*(xkn(k+1)-xkn(k-1))
             enddo
             dxkn(1) = (xkn(2)-xkn(1))
             dxkn(Nmode) = (xkn(Nmode)-xkn(Nmode-1))

             ! Spectra parameter
             ! -----------------
             ! rms amplitude of fluctuations
             uf=ampl_TBL
             ! max and peak wavenumbers
             xke=xkmx*sqrt(5.0_wp/12.0_wp)
             ! dissipation (von Karman cst kappa=0.41)
             if (CHAN) then ! TO BE CHANGED
                xeps=2.0_wp*rho_ref*uf**3/(0.41_wp*delta)
             else
                xeps=0.9_wp*uf**3/Lf_tbl
             endif
             ! Kolmogorov scale
             if (xkkol.eq.0.0) then
                xkkol=(xeps/(mu_ref/rho_ref)**3)**0.25_wp
             else
                xkkol=1.0_wp/xkkol
             endif

             ! Energy spectrum
             ! ---------------
             do k=1,Nmode
                xkk=xkn(k)/xke
                vkm(k)=1.453_wp*uf**2/xke*xkk**4/exp(17.0_wp/6.0_wp*log(1.0_wp+xkk**2)) &
                      * exp(-1.5_wp*1.613_wp*((xkn(k)/xkkol)**2)) &
                      * (1.0_wp+0.522_wp*(0.5_wp+atan(10.0_wp*log10(xkn(k)/xkkol)+12.58_wp)/pi))
             enddo

             ! Normalization factor
             ! --------------------
             norm=0.0_wp
             do k=1,Nmode
                norm = norm + vkm(k)*dxkn(k)
             enddo

             ! Deduce mode amplitude from energy spectrum
             ! ------------------------------------------
             utilde = sqrt(vkm*dxkn)
             norm=0.0_wp
             do k=1,Nmode
                norm=norm+utilde(k)**2
             enddo

             ! Compute divergence-free stochastic field
             ! ----------------------------------------
             ! memory allocation
             allocate(sigma1(Nmode),sigma2(Nmode),sigma3(Nmode))
             allocate(xk1(Nmode),xk2(Nmode),xk3(Nmode))

             ! Random number generation for angles and directions
             ! --------------------------------------------------
             var = 0.0_wp
             do k=1,Nmode
                ! angle alpha (and sine & cosine)
                call random_number(var)
                alpha = 2.0_wp*pi*var
                calpha = cos(alpha)
                salpha = sin(alpha)

                ! angle phi (and sine & cosine)
                call random_number(var)
                phi = 2.0_wp*pi*var
                cphi = cos(phi)
                sphi = sin(phi)

                ! sine & cosine of theta
                call random_number(var)
                ctheta = 1.0_wp-2.0_wp*var
                stheta = sqrt(1.0_wp-ctheta*ctheta)

                ! phase psi
                call random_number(var)
                psi_TBL(k,jbc) = 2.0_wp*pi*var

                ! orientation vector (sigma1,sigma2,sigma3)
                sigma1(k) = calpha*ctheta*cphi - salpha*sphi
                sigma2(k) = calpha*ctheta*sphi + salpha*cphi
                sigma3(k) =-calpha*stheta

                ! wavenumber components on unitary sphere
                xk1(k) = stheta*cphi
                xk2(k) = stheta*sphi
                xk3(k) = ctheta
             enddo

             ! multiply directions by mode amplitude (from spectrum)
             ! -------------------------------------
             do k=1,Nmode
                ! sigma = 2 * û * an
                sigma1(k)=2.0_wp*utilde(k)*sigma1(k)
                sigma2(k)=2.0_wp*utilde(k)*sigma2(k)
                sigma3(k)=2.0_wp*utilde(k)*sigma3(k)
             enddo

             ! multiply wavenumber components by wavenumber amplitude (from spectrum)
             ! ------------------------------------------------------
             do k=1,Nmode
                xk1(k)=xkn(k)*xk1(k)
                xk2(k)=xkn(k)*xk2(k)
                xk3(k)=xkn(k)*xk3(k)
             enddo

             ! Interpolate on current global grid
             ! ----------------------------------
             allocate(uu_rms(ny),vv_rms(ny),ww_rms(ny))
             allocate(uv_rms(ny),uw_rms(ny),vw_rms(ny))

             if (is_curv3) then
                call mpistop('Not implemented in mod_RFM.f90',1)
             else if (is_curv) then
                call mpistop('Not implemented in mod_RFM.f90',1)
             else
                call interp1(y(1:ny)/d99_TBL(jbc),uu_rms,ny,y_db,uu_db,ny_db,'linear')
                call interp1(y(1:ny)/d99_TBL(jbc),vv_rms,ny,y_db,vv_db,ny_db,'linear')
                call interp1(y(1:ny)/d99_TBL(jbc),ww_rms,ny,y_db,ww_db,ny_db,'linear')
                call interp1(y(1:ny)/d99_TBL(jbc),uv_rms,ny,y_db,uv_db,ny_db,'linear')
                call interp1(y(1:ny)/d99_TBL(jbc),uw_rms,ny,y_db,uw_db,ny_db,'linear')
                call interp1(y(1:ny)/d99_TBL(jbc),vw_rms,ny,y_db,vw_db,ny_db,'linear')
             endif
             
             ! Find transformation for which the off-diagonal terms in the tensor are equal to zero
             ! ------------------------------------------------------------------------------------
             ! allocate eigenvalues and right eigenvector array for each j of local grid
             allocate(wr(3,ny))

             ! Find normalization factor
             norm_=0.0_wp
             do j=1,ny
                norm_=MAX(norm_,(uu_rms(j)+vv_rms(j)+ww_rms(j))/3.0_wp)
             enddo
             norm_proc=norm_
             call MPI_ALLREDUCE(norm_proc,norm_,1,MPI_DOUBLE_PRECISION,MPI_MAX,COMMX,info)

             do j=1,ny
                if (uu_rms(j)<0.0_wp) uu_rms(j)=0.0_wp
                if (vv_rms(j)<0.0_wp) vv_rms(j)=0.0_wp
                if (ww_rms(j)<0.0_wp) ww_rms(j)=0.0_wp
             enddo
             
             ! form Reynolds-stress tensor tau for each j on local grid (symmetric 3x3 matrix)
             do j=1,ny
                tau(1,1) = uu_rms(j)
                tau(1,2) = uv_rms(j)
                tau(1,3) = uw_rms(j)
                tau(2,1) = uv_rms(j)
                tau(2,2) = vv_rms(j)
                tau(2,3) = vw_rms(j)
                tau(1,3) = uw_rms(j)
                tau(2,3) = vw_rms(j)
                tau(3,3) = ww_rms(j)
                ! Normalization
                tau(:,:) = tau(:,:)/norm_
                ! Diagonalize Reynolds-stress tensor
                call dgeev('N','V',3,tau,3,wr(:,j),wi,vl,ldvl,vr_TBL(:,:,j,jbc),ldvr,work,lwork,info)
             enddo

             ! Take square root of eigenvalues
             ! -------------------------------
             ! (the velocity field is scaled by wr^1/2)
             wr=sqrt(wr)

             ! Transformed (scaled) wavenumbers and directions
             ! -----------------------------------------------
             ! -> sigma_scaled=R^T sigma
             ! -> k_scaled=R^T k
             do j=1,ny
                do k=1,Nmode
                   xk1s_TBL(k,j,jbc)= xk1(k)*vr_TBL(1,1,j,jbc)+xk2(k)*vr_TBL(2,1,j,jbc)+xk3(k)*vr_TBL(3,1,j,jbc)
                   xk2s_TBL(k,j,jbc)= xk1(k)*vr_TBL(1,2,j,jbc)+xk2(k)*vr_TBL(2,2,j,jbc)+xk3(k)*vr_TBL(3,2,j,jbc)
                   xk3s_TBL(k,j,jbc)= xk1(k)*vr_TBL(1,3,j,jbc)+xk2(k)*vr_TBL(2,3,j,jbc)+xk3(k)*vr_TBL(3,3,j,jbc)
                   sigma1s_TBL(k,j,jbc)= sigma1(k)*vr_TBL(1,1,j,jbc)+sigma2(k)*vr_TBL(2,1,j,jbc)+sigma3(k)*vr_TBL(3,1,j,jbc)
                   sigma2s_TBL(k,j,jbc)= sigma1(k)*vr_TBL(1,2,j,jbc)+sigma2(k)*vr_TBL(2,2,j,jbc)+sigma3(k)*vr_TBL(3,2,j,jbc)
                   sigma3s_TBL(k,j,jbc)= sigma1(k)*vr_TBL(1,3,j,jbc)+sigma2(k)*vr_TBL(2,3,j,jbc)+sigma3(k)*vr_TBL(3,3,j,jbc)
                enddo
             enddo

             ! Scale directions by normalized stress (wr^1/2)
             ! -------------------------------------
             do j=1,ny
                do k=1,Nmode
                   sigma1s_TBL(k,j,jbc) = sigma1s_TBL(k,j,jbc)*wr(1,j)
                   sigma2s_TBL(k,j,jbc) = sigma2s_TBL(k,j,jbc)*wr(2,j)
                   sigma3s_TBL(k,j,jbc) = sigma3s_TBL(k,j,jbc)*wr(3,j)
                enddo
             enddo

             ! Determine angular frequency of turbulence evolution
             ! ---------------------------------------------------
             time_turb='K'
             if (time_turb=='K') then
                ! angular frequency based on Kolmogorov's analysis
                do k=1,Nmode
                   omn_TBL(k,jbc)=xeps**(1.0_wp/3.0_wp)*xkn(k)**(2.0_wp/3.0_wp)
                enddo
             elseif (time_turb=='H') then
                ! angular frequency based on Heisenberg's time
                do k=1,Nmode
                   omn_TBL(k,jbc)=2.0_wp*pi*uf*xkn(k)
                enddo
             elseif (time_turb=='N') then
                do k=1,Nmode
                   omn_TBL(k,jbc)=0.0_wp
                enddo
             else
                call mpistop('problem with definition of turbulence time, please check!',0)
             endif

             ! Free unnecessary memory
             ! -----------------------
             deallocate(xkn,dxkn,vkm,utilde)
             deallocate(xk1,xk2,xk3,sigma1,sigma2,sigma3)
             deallocate(uu_rms,vv_rms,ww_rms,uv_rms,uw_rms,vw_rms)
             deallocate(wr)
             deallocate(y_db,um_db,uu_db,vv_db,ww_db,uv_db,uw_db,vw_db)
          endif
       enddo

    else
       ! To be able to use COMMX
       norm_proc=0.0_wp; norm_=0.0_wp
       call MPI_ALLREDUCE(norm_proc,norm_,1,MPI_DOUBLE_PRECISION,MPI_MAX,COMMX,info)
    endif

    ! ----------------------------------------
    !         I/ Freestream turbulence
    ! ----------------------------------------
    ! ~> if is_damping or TBL (injection zone):
    !    _ Definition of damping coefficient
    ! ~> if is_init_modes:
    !    _ Discretization of wavenumbers
    !    _ Random number generation
    ! ~> Computation of spectrum
    ! ~> if anisotropy:
    !    _ generation of anistropic field
    !    _ only working for CHAN for the moment
    if (is_RFM_FST) then
       ! Reading wavenumbers discretization (if necessary)
       ! ----------------------------------
       ! Memory allocation
       allocate(xkn(Nmode),dxkn(Nmode),vkm(Nmode),utilde(Nmode))
       ! Reading
       if (.not.is_init_modes) then
          inquire(file=trim('vkm.bin'), exist=iexist)
          if (.not.iexist) then
             call mpistop('vkm.bin does not exist!', 0)
          endif
          open(30,file='vkm.bin',form='unformatted',status='old',action='read')
          rewind(30)
          read(30) Nmode
          read(30) (xkn(k),k=1,Nmode)
          read(30) (dxkn(k),k=1,Nmode)
          close(30)
       endif

       ! Damping height
       ! --------------
       is_damping = .false.
       if ((ndy_RFM.ne.1).or.(nfy_rfm.ne.ngy)) is_damping = .true.
       if (is_damping) then
          if (ndy_RFM.ne.1)   h_damp(1) = h_damp(1)/(0.137_wp*(-log(1.0_wp-(0.01_wp)**0.1)))
          if (nfy_rfm.ne.ngy) h_damp(2) = h_damp(2)/(0.137_wp*(-log(1.0_wp-(0.01_wp)**0.1)))
       endif

       ! Scatter injection limit coordinates
       ! -----------------------------------
       cinj_jmin = int(ndy_RFM/ny); cinj_jmax = int(nfy_rfm/ny)
       if (is_curv3) then
          xinj3_jmin=0.0_wp; xinj3_jmax=0.0_wp; xinj3_jmin_proc = 0.0_wp; xinj3_jmax_proc = 0.0_wp
          yinj3_jmin=0.0_wp; yinj3_jmax=0.0_wp; yinj3_jmin_proc = 0.0_wp; yinj3_jmax_proc = 0.0_wp
          zinj3_jmin=0.0_wp; zinj3_jmax=0.0_wp; zinj3_jmin_proc = 0.0_wp; zinj3_jmax_proc = 0.0_wp
          ! Location at jmin
          if (coord(2)==cinj_jmin) then
            xinj3_jmin_proc(1:nz)=xc3(1,ndy_RFM-coord(2)*ny,1:nz)
            yinj3_jmin_proc(1:nz)=yc3(1,ndy_RFM-coord(2)*ny,1:nz)
            zinj3_jmin_proc(1:nz)=zc3(1,ndy_RFM-coord(2)*ny,1:nz)
          endif
          call MPI_ALLREDUCE(xinj3_jmin_proc(1:nz),xinj3_jmin(1:nz),nz,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXZ,info)
          call MPI_ALLREDUCE(yinj3_jmin_proc(1:nz),yinj3_jmin(1:nz),nz,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXZ,info)
          call MPI_ALLREDUCE(zinj3_jmin_proc(1:nz),zinj3_jmin(1:nz),nz,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXZ,info)
          ! Location at jmax
          if (coord(2)==cinj_jmax) then
            xinj3_jmax_proc(1:nz)=xc3(1,nfy_rfm-coord(2)*ny,1:nz)
            yinj3_jmax_proc(1:nz)=yc3(1,nfy_rfm-coord(2)*ny,1:nz)
            zinj3_jmax_proc(1:nz)=zc3(1,nfy_rfm-coord(2)*ny,1:nz)
          endif
          call MPI_ALLREDUCE(xinj3_jmax_proc(1:nz),xinj3_jmax(1:nz),nz,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXZ,info)
          call MPI_ALLREDUCE(yinj3_jmax_proc(1:nz),yinj3_jmax(1:nz),nz,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXZ,info)
          call MPI_ALLREDUCE(zinj3_jmax_proc(1:nz),zinj3_jmax(1:nz),nz,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXZ,info)
       else if (is_curv) then
          xinj_jmin=0.0_wp; xinj_jmax=0.0_wp; xinj_jmin_proc = 0.0_wp; xinj_jmax_proc = 0.0_wp
          yinj_jmin=0.0_wp; yinj_jmax=0.0_wp; yinj_jmin_proc = 0.0_wp; yinj_jmax_proc = 0.0_wp
          ! Location at jmin
          if (coord(2)==cinj_jmin) then
            xinj_jmin_proc=xc(1,ndy_RFM-coord(2)*ny); yinj_jmin_proc=yc(1,ndy_RFM-coord(2)*ny)
          endif
          call MPI_ALLREDUCE(xinj_jmin_proc,xinj_jmin,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXZ,info)
          call MPI_ALLREDUCE(yinj_jmin_proc,yinj_jmin,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXZ,info)
          ! Location at jmax
          if (coord(2)==cinj_jmax) then
            xinj_jmax_proc=xc(1,nfy_rfm-coord(2)*ny); yinj_jmax_proc=yc(1,nfy_rfm-coord(2)*ny)
          endif
          call MPI_ALLREDUCE(xinj_jmax_proc,xinj_jmax,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXZ,info)
          call MPI_ALLREDUCE(yinj_jmax_proc,yinj_jmax,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXZ,info)
       else
          yinj_jmin=0.0_wp; yinj_jmax=0.0_wp; yinj_jmin_proc = 0.0_wp; yinj_jmax_proc = 0.0_wp
          ! Location at jmin
          if (coord(2)==cinj_jmin)  yinj_jmin_proc=y(ndy_RFM-coord(2)*ny)
          call MPI_ALLREDUCE(yinj_jmin_proc,yinj_jmin,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXZ,info)
          ! Location at jmax
          if (coord(2)==cinj_jmax)  yinj_jmax_proc=y(nfy_rfm-coord(2)*ny)
          call MPI_ALLREDUCE(yinj_jmax_proc,yinj_jmax,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXZ,info)
       endif

       ! Calculation of damping coefficient
       ! ----------------------------------
       if (is_damping) then
          ! Allocation
          if (is_curv3) then
             allocate(damping_coeff3(ndyt:nfyt,1:nz))
             damping_coeff = 1.0_wp
          else
             allocate(damping_coeff(ndyt:nfyt))
             damping_coeff = 1.0_wp
          endif
          ! Computation
          if ((h_damp(1).gt.0.0_wp).and.(h_damp(2).gt.0.0_wp)) then
             if (is_curv3) then
                ! Distance to the injection location: ||X-Xinj||_2  with  X=(x y z)^T,  Xinj(k)=(xinj yinj zinj)^T
                do j=ndyt,nfyt
                   if ((j+coord(2)*ny.gt.ndy_RFM).and.(j+coord(2)*ny.lt.nfy_rfm)) then
                      do k=1,nz
                         damping_coeff3(j,k) = (1 - exp(-((xc3(1,j,k)-xinj3_jmin(k))**2 + (yc3(1,j,k)-yinj3_jmin(k))**2 + (zc3(1,j,k)-zinj3_jmin(k))**2)**0.5/(0.137_wp*h_damp(1))))**10 * &
                                               (1 - exp(-((xc3(1,j,k)-xinj3_jmax(k))**2 + (yc3(1,j,k)-yinj3_jmax(k))**2 + (zc3(1,j,k)-zinj3_jmax(k))**2)**0.5/(0.137_wp*h_damp(1))))**10
                      enddo
                   else
                      damping_coeff3(j,:) = 0.0_wp
                   endif
                enddo
             else if (is_curv) then
                ! Distance to the injection location: ||X-Xinj||_2  with  X=(x y)^T,  Xinj(k)=(xinj yinj)^T
                do j=ndyt,nfyt
                   if ((j+coord(2)*ny.gt.ndy_RFM).and.(j+coord(2)*ny.lt.nfy_rfm)) then
                      damping_coeff(j) = (1 - exp(-((xc(1,j)-xinj_jmin)**2 + (yc(1,j)-yinj_jmin)**2)**0.5/(0.137_wp*h_damp(1))))**10 * &
                                         (1 - exp(-((xinj_jmax-xc(1,j))**2 + (yinj_jmax-yc(1,j))**2)**0.5/(0.137_wp*h_damp(2))))**10
                   else
                      damping_coeff(j) = 0.0_wp
                   endif
                enddo
             else
                ! Distance to the injection location: |x-xinj|
                do j=ndyt,nfyt
                   if ((j+coord(2)*ny.gt.ndy_RFM).and.(j+coord(2)*ny.lt.nfy_rfm)) then
                      damping_coeff(j) = (1 - exp(-abs(y(j)-yinj_jmin)/(0.137_wp*h_damp(1))))**10 * &
                                         (1 - exp(-abs(yinj_jmax-y(j))/(0.137_wp*h_damp(2))))**10
                   else
                      damping_coeff(j) = 0.0_wp
                   endif
                enddo
             endif
          elseif (h_damp(1).gt.0.0_wp) then
             if (is_curv3) then
                do j=ndyt,nfyt
                   if (j+coord(2)*ny.gt.ndy_RFM) then
                      do k=1,nz
                         damping_coeff3(j,k) = (1 - exp(-((xc3(1,j,k)-xinj3_jmin(k))**2 + (yc3(1,j,k)-yinj3_jmin(k))**2 + (zc3(1,j,k)-zinj3_jmin(k))**2)**0.5/(0.137_wp*h_damp(1))))**10
                      enddo
                   else
                      damping_coeff3(j,:) = 0.0_wp
                   endif
                enddo
             else if (is_curv) then
                do j=ndyt,nfyt
                   if (j+coord(2)*ny.gt.ndy_RFM) then
                      damping_coeff(j) = (1 - exp(-((xc(1,j)-xinj_jmin)**2 + (yc(1,j)-yinj_jmin)**2)**0.5/(0.137_wp*h_damp(1))))**10
                   else
                      damping_coeff(j) = 0.0_wp
                   endif
                enddo
             else
                do j=ndyt,nfyt
                   if (j+coord(2)*ny.gt.ndy_RFM) then
                      damping_coeff(j) = (1 - exp(-abs(y(j)-yinj_jmin)/(0.137_wp*h_damp(1))))**10
                   else
                      damping_coeff(j) = 0.0_wp
                   endif
                enddo
             endif
          elseif (h_damp(2).gt.0.0_wp) then
             if (is_curv3) then
                do j=ndyt,nfyt
                   if (j+coord(2)*ny.lt.nfy_rfm) then
                      do k=1,nz
                         damping_coeff3(j,k) = (1 - exp(-((xc3(1,j,k)-xinj3_jmax(k))**2 + (yc3(1,j,k)-yinj3_jmax(k))**2 + (zc3(1,j,k)-zinj3_jmax(k))**2)**0.5/(0.137_wp*h_damp(1))))**10
                      enddo
                   else
                      damping_coeff3(j,:) = 0.0_wp
                   endif
                enddo
             else if (is_curv) then
                do j=ndyt,nfyt
                   if (j+coord(2)*ny.lt.nfy_rfm) then
                      damping_coeff(j) = (1 - exp(-((xinj_jmax-xc(1,j))**2 + (yinj_jmax-yc(1,j))**2)**0.5/(0.137_wp*h_damp(2))))**10
                   else
                      damping_coeff(j) = 0.0_wp
                   endif
                enddo
             else
                do j=ndyt,nfyt
                   if (j+coord(2)*ny.lt.nfy_rfm) then
                      damping_coeff(j) = (1 - exp(-abs(yinj_jmax-y(j))/(0.137_wp*h_damp(2))))**10
                   else
                      damping_coeff(j) = 0.0_wp
                   endif
                enddo
             endif
          endif
       endif

       ! Output file to check damping coeff
       ! ----------------------------------
       if (is_damping) then
          if (is_curv3) then
             ! To be generalized for ngz ~> needs to communicate data between proc in direction k
             allocate(damping_coeff_all(1:ngy,1:nz))
             do k=1,nz
                call MPI_Gather(damping_coeff3(1:ny,k), ny, MPI_DOUBLE_PRECISION, damping_coeff_all(1:ngy,k), ny, MPI_DOUBLE_PRECISION, iproc_leader(nob(iproc)), COMMXZ, info)
             enddo
          else
             allocate(damping_coeff_all(1:ngy,1))
             call MPI_Gather(damping_coeff(1:ny), ny, MPI_DOUBLE_PRECISION, damping_coeff_all(1:ngy,1), ny, MPI_DOUBLE_PRECISION, iproc_leader(nob(iproc)), COMMXZ, info)
          endif
          if (iproc.eq.iproc_leader(nob(iproc))) then
             open(30,file='damping_coeff_bl'//trim(numchar(nob(iproc)))//'.bin',form='unformatted',status='unknown')
             rewind(30)
             if (is_curv3) then
                write(30) ngy, nz
                write(30) ((damping_coeff_all(j,k),j=1,ngy),k=1,nz)
             else
                write(30) ngy
                write(30) (damping_coeff_all(j,1),j=1,ngy)
             endif
             close(30)
          endif
          deallocate(damping_coeff_all)
       endif

       ! Calculation of kmin, kmax and kmn
       ! ---------------------------------
       delta = Lf_int/(sqrt(12.0_wp/5.0_wp)*0.74684_wp)
       xkmx = 1/delta
       if (is_xkm_given) then
          if (iproc.eq.0) print *,"xkmin and xkmax specified in the param file"
       else
          if (iproc.eq.0) print *,"xkmin and xkmax calculated based on grid dimension"
          if ((nbloc.eq.1).and.(.not.is_curv)) then
             if (is_2D) then
                xkmax = 2.0_wp*pi*0.21739_wp / abs((yg(ngy/2) - yg(ngy/2-1))) ! TO BE CHANGED ~> depends of the scheme, here 11-pts DRP Bogey & Bailly
                xkmin = max(min(2.0_wp*pi/(abs(yg(ngy)-yg(1))),0.8_wp*xkmx),0.237_wp*xkmx)
             else
                xkmax = 2.0_wp*pi*0.21739_wp / max(abs((zg(ngz/2) - zg(ngz/2-1))),abs((xg(ngx/2) - xg(ngx/2-1))))! TO BE CHANGED ~> depends of the scheme, here 11-pts DRP Bogey & Bailly
                xkmin = max(min(2.0_wp*pi/(abs(zg(ngz)-zg(1))),0.8_wp*xkmx),0.237_wp*xkmx)
             endif
          else
                if (iproc==0) then
                   print *,"============================================"
                   print *,"                 ATTENTION"
                   print *,"              KMIN MIS EN DUR"
                   print *,"                 ATTENTION"
                   print *,"============================================"
                endif
                ! xkmin = 2.0_wp*pi/(Lf_int*9)
                ! ! dxmax defined in grid.f90
                ! xkmax = 2.0_wp*pi*0.21739_wp / dxmax ! TO BE CHANGED ~> depends of the scheme, here 11-pts DRP Bogey & Bailly
                ! xkmax = xkmin*64*0.21739

                ! ! Temporary
                ! ! xkmin = 0.8_wp*xkmx
                ! xkmax = 2.0_wp*pi/(Lf_int*9)*64*0.21739

                ! Temporary ~> inlet LS89 turb
                xkmax = 2.0_wp*pi/(Lf_int*9)*64*0.21739
                xkmin = 2.0_wp*pi/(Lf_int*9)
          endif
       endif

       ! Print RFM parameters
       ! --------------------
       if (iproc.eq.iproc_leader_rfm) then
          if (kdist=='log') then
             write(6,20) Nmode
          elseif (kdist=='lin') then
             write(6,21) Nmode
          endif
          write(6,22) xkmin,xkmax,xkmx
          if (xkmin>xkmx) print *,'  ~> xkmin too big compared to the maximum position of the energy spectrum'
          if (xkmax<xkmx) print *,'  ~> xkmax too small compared to the maximum position of the energy spectrum'
          if (.not.is_curv) then
             if (is_2D) then
                if (1.0_wp/(abs(ymax-ymin)).ge.xkmx) &
                     print *,'  ~> direction y too small compared to the maximum position of the energy spectrum'
             else
                if (1.0_wp/(abs(zmax-zmin)).ge.xkmx) &
                     print *,'  ~> direction z too small compared to the maximum position of the energy spectrum'
             endif
          endif
          if (time_turb=='H') then
             write(6,*) '    unfrozen field with Heisenberg''s time'
          elseif (time_turb=='K') then
             write(6,*) '    unfrozen field with Kolmogorov''s time'
          endif
       endif
   20  format(5x,'logarithmic distribution of ',i4,' RFM modes')
   21  format(5x,'linear distribution of ',i4,' RFM modes')
   22  format(5x,'between k_min=',f11.3,' and k_max=',f11.3,' with peak at k=',f11.3)


       ! ! Initialization of random number generator
       ! !------------------------------------------
       ! call random_seed(size=nseed)
       ! allocate(initseed(nseed))
       ! do k=1,nmax
       !    call random_number(randomnb)
       ! enddo
       ! do k=1,nseed
       !    call random_number(randomnb)
       !    initseed(k) = int(randomnb*10**8) * (-1)**k
       ! enddo
       ! call random_seed(put=initseed)
       ! deallocate(initseed)

       ! Wavenumber discretization
       ! =========================
       if (is_init_modes) then
          ! linear or logarithmic wavenumbers discretization
          ! ------------------------------------------------
          if (kdist=='log') then
             dkl=(log(xkmax)-log(xkmin))/dble(Nmode-1)
             do k=1,Nmode
                xkn(k)=exp(log(xkmin)+dble(k-1)*dkl)
             enddo
          else if (kdist=='lin') then
             dkl=(xkmax-xkmin)/dble(Nmode-1)
             do k=1,Nmode
                xkn(k) = xkmin+dble(k-1)*dkl
             enddo
          else
             call mpistop('Wrong choice for the wavenumber discretization in param_RFM !',0)
          endif

          ! wavenumber increments
          ! ---------------------
          do k=2,Nmode-1
             dxkn(k)=0.5_wp*(xkn(k+1)-xkn(k-1))
          enddo
          dxkn(1) = (xkn(2)-xkn(1))
          dxkn(Nmode) = (xkn(Nmode)-xkn(Nmode-1))
       endif

       ! Compute von Karman-Pao energy spectrum
       ! ======================================

       ! Parameters
       ! ----------
       ! rms amplitude of fluctuations
       uf=ampl_RFM
       ! max and peak wavenumbers
       xke=xkmx*sqrt(5.0_wp/12.0_wp)
       ! dissipation (von Karman cst kappa=0.41)
       if (CHAN) then ! TO BE CHANGED
          xeps=2.0_wp*rho_ref*uf**3/(0.41_wp*delta)
       else
          xeps=0.9_wp*uf**3/Lf_int
       endif
       ! Kolmogorov scale
       if (xkkol.eq.0.0) then
          xkkol=(xeps/(mu_ref/rho_ref)**3)**0.25_wp
       else
          xkkol=1.0_wp/xkkol
       endif

       ! Energy spectrum
       ! ---------------
       do k=1,Nmode
          xkk=xkn(k)/xke
          ! vkm(k)=1.453_wp*uf**2/xke*xkk**4/exp(17.0_wp/6.0_wp*log(1.0_wp+xkk**2)) &
          !       * exp(-2.0_wp*((xkn(k)/xkkol)**2))
          ! vkm(k)=1.453_wp*uf**2/xke*xkk**4/exp(17.0_wp/6.0_wp*log(1.0_wp+xkk**2)) &
          !         * exp(-1.5_wp*1.613_wp*((xkn(k)/xkkol)**2))
          vkm(k)=1.453_wp*uf**2/xke*xkk**4/exp(17.0_wp/6.0_wp*log(1.0_wp+xkk**2)) &
                * exp(-1.5_wp*1.613_wp*((xkn(k)/xkkol)**2)) &
                * (1.0_wp+0.522_wp*(0.5_wp+atan(10.0_wp*log10(xkn(k)/xkkol)+12.58_wp)/pi))
       enddo

       ! Check spectrum
       ! --------------
       if (iproc.eq.iproc_leader_rfm) then
          open(30,file='vkm.bin',form='unformatted',status='unknown')
          rewind(30)
          write(30) Nmode
          write(30) (xkn(k),k=1,Nmode)
          write(30) (dxkn(k),k=1,Nmode)
          write(30) (vkm(k),k=1,Nmode)
          close(30)
       endif

       ! normalization factor
       ! --------------------
       norm=0.0_wp
       do k=1,Nmode
          norm = norm + vkm(k)*dxkn(k)
       enddo
       ! vkm=1.5_wp*vkm/norm
       ! vkm = vkm/norm
       ! if (iproc==0) print *,'    normalization factor:',norm

       ! Deduce mode amplitude from energy spectrum
       ! ------------------------------------------
       utilde = sqrt(vkm*dxkn)
       norm=0.0_wp
       do k=1,Nmode
          norm=norm+utilde(k)**2
       enddo
       ! if (iproc==0) print *,'    kinetic energy:',norm

       ! Compute divergence-free stochastic field
       ! ========================================

       ! memory allocation
       allocate(sigma1(Nmode),sigma2(Nmode),sigma3(Nmode))
       allocate(xk1(Nmode),xk2(Nmode),xk3(Nmode))
       allocate(psi(Nmode),omn(Nmode))

       if (is_init_modes) then
          ! Random number generation for angles and directions
          ! --------------------------------------------------
          var = 0.0_wp
          do k=1,Nmode

             ! angle alpha (and sine & cosine)
             call random_number(var)
             alpha = 2.0_wp*pi*var
             calpha = cos(alpha)
             salpha = sin(alpha)

             ! angle phi (and sine & cosine)
             call random_number(var)
             phi = 2.0_wp*pi*var
             cphi = cos(phi)
             sphi = sin(phi)

             ! sine & cosine of theta
             call random_number(var)
             ctheta = 1.0_wp-2.0_wp*var
             stheta = sqrt(1.0_wp-ctheta*ctheta)

             ! phase psi
             call random_number(var)
             psi(k) = 2.0_wp*pi*var

             ! orientation vector (sigma1,sigma2,sigma3)
             sigma1(k) = calpha*ctheta*cphi - salpha*sphi
             sigma2(k) = calpha*ctheta*sphi + salpha*cphi
             sigma3(k) =-calpha*stheta

             ! wavenumber components on unitary sphere
             xk1(k) = stheta*cphi
             xk2(k) = stheta*sphi
             xk3(k) = ctheta
          enddo

          ! Check wavenumbers and orientations (should be on unitary sphere)
          ! ----------------------------------
          if (iproc.eq.iproc_leader_rfm) then
             open(30,file='k_random.bin',form='unformatted',status='unknown')
             rewind(30)
             write(30) Nmode
             write(30) (xk1(k),k=1,Nmode)
             write(30) (xk2(k),k=1,Nmode)
             write(30) (xk3(k),k=1,Nmode)
             write(30) (sigma1(k),k=1,Nmode)
             write(30) (sigma2(k),k=1,Nmode)
             write(30) (sigma3(k),k=1,Nmode)
             write(30) (psi(k),k=1,Nmode)
             close(30)
          endif

       else
          ! Reading k_random.bin file
          ! -------------------------
          inquire(file=trim('k_random.bin'), exist=iexist)
          if (.not.iexist) then
             call mpistop('k_random.bin does not exist!', 0)
          endif
          open(30,file='k_random.bin',form='unformatted',status='old',action='read')
          rewind(30)
          read(30) Nmode
          read(30) (xk1(k),k=1,Nmode)
          read(30) (xk2(k),k=1,Nmode)
          read(30) (xk3(k),k=1,Nmode)
          read(30) (sigma1(k),k=1,Nmode)
          read(30) (sigma2(k),k=1,Nmode)
          read(30) (sigma3(k),k=1,Nmode)
          read(30) (psi(k),k=1,Nmode)
          close(30)
       endif

       ! multiply directions by mode amplitude (from spectrum)
       ! -------------------------------------
       do k=1,Nmode
          ! sigma = 2 * û * an
          sigma1(k)=2.0_wp*utilde(k)*sigma1(k)
          sigma2(k)=2.0_wp*utilde(k)*sigma2(k)
          sigma3(k)=2.0_wp*utilde(k)*sigma3(k)
       enddo

       ! multiply wavenumber components by wavenumber amplitude (from spectrum)
       ! ------------------------------------------------------
       do k=1,Nmode
          xk1(k)=xkn(k)*xk1(k)
          xk2(k)=xkn(k)*xk2(k)
          xk3(k)=xkn(k)*xk3(k)
       enddo

       ! Transformation to take anisotropy into account
       ! ==============================================
       ! ~> based on 'Random Flow Generation', Smirnov et al., J.Fluid Eng. 123 (2001)]
       ! ~> method to introduce anisotropy described in:
       !    Billson, Eriksson & Davidson,  AIAA paper 2004-2857
       if (anisotropy.ne.'N') then
          call mpistop('This part needs to be cleaned',1)

          ! Read Reynolds stresses
          ! ----------------------
          call read_database(base_FST,Re_prof_FST)

          ! Interpolate on current global grid
          ! ----------------------------------
          allocate(uu_rms(ngy),vv_rms(ngy),ww_rms(ngy))
          allocate(uv_rms(ngy),uw_rms(ngy),vw_rms(ngy))

          if (CHAN) then
             Ldim_a = hc
          else if (STBL) then
             Ldim_a = L_ref
             call mpistop('Erreur ici mod_RFM.f90 l.642',1)
          else if (TE) then
             Ldim_a = L_ref
             call mpistop('Why is it L_ref ??? To be changed',1)
          else
             call mpistop('Anistropy in mod_RFM.f90 not defined for cases other than CHAN, STBL and TE',1)
          endif

          call interp1(yg(1:ngy)/Ldim_a,uu_rms,ngy,y_db,uu_db,ny_db,'linear')
          call interp1(yg(1:ngy)/Ldim_a,vv_rms,ngy,y_db,vv_db,ny_db,'linear')
          call interp1(yg(1:ngy)/Ldim_a,ww_rms,ngy,y_db,ww_db,ny_db,'linear')
          call interp1(yg(1:ngy)/Ldim_a,uv_rms,ngy,y_db,uv_db,ny_db,'linear')
          call interp1(yg(1:ngy)/Ldim_a,uw_rms,ngy,y_db,uw_db,ny_db,'linear')
          call interp1(yg(1:ngy)/Ldim_a,vw_rms,ngy,y_db,vw_db,ny_db,'linear')

          deallocate(uu_db,vv_db,ww_db,uv_db,uw_db,vw_db)

          ! ! Define mean field (for mean flow convection or field initialization)
          ! ! -----------------
          ! allocate(u_m(ny,nz))
          ! if ((is_convec).or.(is_field)) then
          !    allocate(um_m(ngy))
          !    call interp1(yg(1:ngy)/hc,um_m,ngy,y_db,um_db,ny_db,'linear')
          !    do k=1,nz
          !       do j=1,ny
          !          u_m=u_ref*um_m(j+coord(2)*ny)
          !       enddo
          !    enddo
          !    deallocate(um_m,um_db)
          ! else
          !    u_m=0.0_wp
          !    deallocate(um_db)
          ! endif
          ! -> done in mod_init_flow now

          ! ! Mean velocity put to 0.6*U_ref
          ! if (.not.allocated(u_m)) allocate(u_m(ny,nz),v_m(ny,nz),w_m(ny,nz))
          ! u_m(:,:) = 0.6*U_ref
          ! v_m(:,:) = 0.0_wp
          ! w_m(:,:) = 0.0_wp

          ! Store local rms profiles to be used as weighting functions
          ! ----------------------------------------------------------
          if (anisotropy=='W') then

             allocate(u_rms(ny),v_rms(ny),w_rms(ny))
             do j=1,ny
                u_rms(j)=sqrt(uu_rms(j+coord(2)*ny))
                v_rms(j)=sqrt(vv_rms(j+coord(2)*ny))
                w_rms(j)=sqrt(ww_rms(j+coord(2)*ny))
             enddo

          elseif (anisotropy=='S') then

             ! Find transformation for which the off-diagonal terms in the tensor are equal to zero
             ! ------------------------------------------------------------------------------------
             ! tau^*=R^T tau R
             ! -> tau^* contains eigenvalues of tau
             ! -> R is a rotation matrix formed by right eigenvectors

             ! allocate eigenvalues and right eigenvector array for each j of local grid
             allocate(wr(3,ny))
             allocate(vr(ldvr,3,ny))
             ! Nota: The right eigenvector v(j) of A satisfies
             !              A * v(j) = lambda(j) * v(j)
             !       where lambda(j) is its eigenvalue.

             ! Find normalization factor
             norm_ = 0.0_wp
             do j=1,ngy
                norm_ = MAX(norm_,(uu_rms(j)+vv_rms(j)+ww_rms(j))/3.0_wp)
             enddo

             ! form Reynolds-stress tensor tau for each j on local grid (symmetric 3x3 matrix)
             do j=1,ny
                tau(1,1)=uu_rms(j+coord(2)*ny)
                tau(1,2)=uv_rms(j+coord(2)*ny)
                tau(1,3)=uw_rms(j+coord(2)*ny)
                tau(2,1)=uv_rms(j+coord(2)*ny)
                tau(2,2)=vv_rms(j+coord(2)*ny)
                tau(2,3)=vw_rms(j+coord(2)*ny)
                tau(1,3)=uw_rms(j+coord(2)*ny)
                tau(2,3)=vw_rms(j+coord(2)*ny)
                tau(3,3)=ww_rms(j+coord(2)*ny)

                tau(:,:) = tau(:,:)/norm_

                ! diagonalize Reynolds-stress tensor
                call dgeev('N','V',3,tau,3,wr(:,j),wi,vl,ldvl,vr(:,:,j),ldvr,work,lwork,info)

                ! ===========================================================================
                ! Recall syntax of LAPACK's dgeev
                ! ===========================================================================
                ! DGEEV computes for an N-by-N real nonsymmetric matrix A, the
                ! eigenvalues and, optionally, the left and/or right eigenvectors.
                ! ===========================================================================
                ! call dgeev(JOBVL,JOBVR,N,A,LDA,WR,WI,VL,LDVL,VR,LDVR,WORK,LWORK,INFO)
                ! ===========================================================================
                ! JOBVL='N': left eigenvectors of A are not computed;
                ! JOBVR='V': right eigenvectors of A are computed;
                ! N: order of the matrix A
                ! A: N-by-N matrix A (overwritten on exit)
                ! LDA: leading dimension of A
                ! WR(N): real part of the computed eigenvalues
                ! WI(N): imaginary part of the computed eigenvalues
                !        [not referenced since A symmetric -> real eigenvalues]
                ! VL(LDVL,N): not referenced since JOBVL='N'
                ! LDVL: leading dim VL
                ! VR(LDVR,N): right eigenvectors v(j) stored one after another in
                !             the columns of VR, in the same order as their eigenvalues.
                ! LDVR: leading dim VR
                ! WORK: work array of size LWORK
                !       [if INFO = 0, WORK(1) returns the optimal LWORK]
                ! LWORK: size of WORK [LWORK >= 4*N]
                ! INFO= 0:  successful exit
                ! ===========================================================================
             enddo

             ! Take square root of eigenvalues
             ! -------------------------------
             ! (the velocity field is scaled by wr^1/2)
             wr=sqrt(wr)

             !The computed eigenvectors are normalized to have Euclidean norm
             !equal to 1 and largest component real.

             ! Transformed (scaled) wavenumbers and directions
             ! -----------------------------------------------
             ! -> sigma_scaled=R^T sigma
             ! -> k_scaled=R^T k
             allocate(sigma1s(Nmode,ny),sigma2s(Nmode,ny),sigma3s(Nmode,ny))
             allocate(xk1s(Nmode,ny),xk2s(Nmode,ny),xk3s(Nmode,ny))

             do j=1,ny
                do k=1,Nmode
                   xk1s(k,j)=xk1(k)*vr(1,1,j)+xk2(k)*vr(2,1,j)+xk3(k)*vr(3,1,j)
                   xk2s(k,j)=xk1(k)*vr(1,2,j)+xk2(k)*vr(2,2,j)+xk3(k)*vr(3,2,j)
                   xk3s(k,j)=xk1(k)*vr(1,3,j)+xk2(k)*vr(2,3,j)+xk3(k)*vr(3,3,j)
                   sigma1s(k,j)=sigma1(k)*vr(1,1,j)+sigma2(k)*vr(2,1,j)+sigma3(k)*vr(3,1,j)
                   sigma2s(k,j)=sigma1(k)*vr(1,2,j)+sigma2(k)*vr(2,2,j)+sigma3(k)*vr(3,2,j)
                   sigma3s(k,j)=sigma1(k)*vr(1,3,j)+sigma2(k)*vr(2,3,j)+sigma3(k)*vr(3,3,j)
                enddo
             enddo

             deallocate(xk1,xk2,xk3,sigma1,sigma2,sigma3)

             ! scale directions by normalized stress (wr^1/2)
             ! -------------------------------------
             do j=1,ny
                do k=1,Nmode
                   sigma1s(k,j)=sigma1s(k,j)*wr(1,j)
                   sigma2s(k,j)=sigma2s(k,j)*wr(2,j)
                   sigma3s(k,j)=sigma3s(k,j)*wr(3,j)
                enddo
             enddo

             ! ! scale wavenumbers by normalized stress (wr^-1/2)
             ! ! --------------------------------------
             ! ! (to ensure divergence-free field such that xks.sigmas=0
             ! do j=1,ny
             !    do k=1,Nmode
             !       if (wr(1,j)==0.0_wp) wr(1,j)=1.e-16_wp
             !       if (wr(2,j)==0.0_wp) wr(2,j)=1.e-16_wp
             !       if (wr(3,j)==0.0_wp) wr(3,j)=1.e-16_wp
             !       xk1s(k,j)=xk1s(k,j)/wr(1,j)
             !       xk2s(k,j)=xk2s(k,j)/wr(2,j)
             !       xk3s(k,j)=xk3s(k,j)/wr(3,j)
             !    enddo
             ! enddo

          endif

          deallocate(uu_rms,vv_rms,ww_rms,uv_rms,uw_rms,vw_rms)
       endif

       ! Determine angular frequency of turbulence evolution
       ! ---------------------------------------------------
       if (time_turb=='K') then
          ! angular frequency based on Kolmogorov's analysis
          do k=1,Nmode
             omn(k) = xeps**(1.0_wp/3.0_wp)*xkn(k)**(2.0_wp/3.0_wp)
          enddo
       elseif (time_turb=='H') then
          ! angular frequency based on Heisenberg's time
          do k=1,Nmode
             omn(k) = 2.0_wp*pi*uf*xkn(k)
          enddo
       elseif (time_turb=='N') then
          do k=1,Nmode
             omn(k)=0.0_wp
          enddo
       else
          call mpistop('problem with definition of turbulence time, please check!',0)
       endif

       ! Free unnecessary memory
       ! -----------------------
       deallocate(xkn,dxkn,vkm,utilde)

       ! Convergence check
       ! -----------------
       if (is_check_conv) then
          if (is_bc_1pt(1,1)) then
             nob_inp = 1
          else
             nob_inp = 5
          endif
         allocate(inlet_ip(nob_inp))
          if (iproc.eq.iproc_leader_rfm) then
             call display_RFM_summary
             print *,''
             print *,'Check convergence of rms profiles ...'
          endif

          call compute_RFM_check
          call mpistop('', 0)
       endif

    endif

    ! Add synthetic turbulence in initial field
    ! -----------------------------------------
    if ((is_field).and.((idepart.eq.1).or.(is_init_2D3D))) then
       if (CHAN) call init_vel_RFM
       if ((SHIT).or.(STBL).or.(TURB).or.(LE).or.(TE)) call init_vel_RFM
    endif

    if (iproc.eq.iproc_leader_rfm) call display_RFM_summary

    ! Check on what is the processor really doing
    ! -------------------------------------------
    if ((coord(2)+1)*ny.lt.ndy_RFM) is_RFM_FST=.false.
    if (1+coord(2)*ny.ge.ndy_RFM)   is_RFM_TBL(1)=.false.
    if ((coord(2)+1)*ny.le.nfy_RFM) is_RFM_TBL(2)=.false.
    if (1+coord(2)*ny.gt.nfy_RFM) is_RFM_FST=.false.

  end subroutine init_RFM

  !===============================================================================
  module subroutine display_RFM_summary
  !===============================================================================
    !> Display on screen information summary on RFM
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    character(len=100) :: out_str
    integer :: ib
    ! ----------------------------------------------------------------------------

    print *,""
    print *,repeat('=',70)
    print *,repeat(' ',22),'Random Fourier Modes (RFM)'
    print *,repeat('=',70)

    if (nbloc.gt.1) then
       out_str = " RFM applied in block(s): "
       if (is_RFM_blocks(1))  write(out_str(len_trim(out_str)+1:), '(A,I0)') '',1
       do ib=2,nbloc
           if (is_RFM_blocks(ib))  write(out_str(len_trim(out_str)+1:), '(A,I0)') ', ',ib
       end do
       print *,out_str
    endif

    out_str = " Application of RFM to generate "
    out_str = trim(out_str) ! Remove trailing spaces
    if ((is_RFM_FST).and.((is_RFM_TBL(1)).or.(is_RFM_TBL(2)))) then
       write(out_str(len_trim(out_str)+1:), '(A)') ' FST & TBL'
    else if (is_RFM_FST) then
       write(out_str(len_trim(out_str)+1:), '(A)') " FST"
    else if (is_RFM_TBL(1).or.is_RFM_TBL(2)) then
       write(out_str(len_trim(out_str)+1:), '(A)') " TBL without FST"
    endif
    write(*, '(A)') out_str

    if ((is_RFM_FST).or.(is_RFM_TBL(1)).or.(is_RFM_TBL(2))) then
       write(*, '(A)') ' ------------------------'
       write(*, '(A)') ' General parameterization'
       write(*, '(A)') ' ------------------------'
       write(*, '(A,I0)') ' _ Number of modes: ',Nmode
       write(*, '(A,I0,A,I0)') ' _ Injection :   ndy_RFM=',ndy_RFM,', nfy_rfm=',nfy_rfm
    endif

    if (is_RFM_FST) then
       write(*, '(A)') ' -------------------'
       write(*, '(A)') ' FST parameterization'
       write(*, '(A)') ' --------------------'
       if (is_init_modes) then
          write(*, '(A)') ' _ Generation of modes discretization:'
          if (is_xkm_given) then
             write(*, '(A)') '   ~> kmin and kmax specified'
          else
             write(*, '(A)') '   ~> kmin and kmax automatically calculated'
          endif
          if (kdist.eq."log") then
             write(*, '(A)') '   ~> logarithmic distribution'
          else
             write(*, '(A)') '   ~> linear distribution'
          endif
       else
          write(*, '(A)') ' _ Modes discretization prescribed with vkm.bin file'
       endif
       ! Inlet spectra
       write(*, '(A)') ' _ Characteristics of inlet target spectra:'
       write(*, '(A,F8.5)') '   ~> Integral length scale [mm]: Lf=',Lf_int*1000
       write(*, '(A,I0)')   '   ~> Reynolds based on Lf and reference values: ',int(Lf_int*rho_ref*(u_m**2+v_m**2+w_m**2)**0.5/mu_ref)
       ! FST intensity
       write(*, '(A,F4.1,A)') ' _ Freestream turbulence intensity: ',Tu_RFM*100,' %'
    endif

    if ((is_RFM_TBL(1)).or.(is_RFM_TBL(2))) then
       write(*, '(A)') ' --------------------'
       write(*, '(A)') ' TBL parameterization'
       write(*, '(A)') ' --------------------'
       ! Turbulence intensity
       write(*, '(A,F4.1,A)') ' _ Turbulence intensity: ',Tu_TBL*100,' %'
       if (is_RFM_TBL(1)) then
          write(*, '(A)') ' _ Characteristics of inlet TBL at imin:'
          write(*, '(A,F8.5)') '   ~> delta_99 [mm]: d99=',d99_TBL(1)*1000
          write(*, '(A,I0)')   '   ~> Reynolds based on d99 and reference values: ',int(d99_TBL(1)*rho_ref*(u_m**2+v_m**2+w_m**2)**0.5/mu_ref)
       if (is_RFM_TBL(2)) then
          write(*, '(A)') ' _ Characteristics of inlet TBL at imax:'
          write(*, '(A,F8.5)') '   ~> delta_99 [mm]: d99=',d99_TBL(2)*1000
          write(*, '(A,I0)')   '   ~> Reynolds based on d99 and reference values: ',int(d99_TBL(2)*rho_ref*(u_m**2+v_m**2+w_m**2)**0.5/mu_ref)
       endif
       endif
    endif

    ! call mpistop('TEST',1)

  end subroutine display_RFM_summary

  ! !===============================================================================
  ! subroutine read_param_RFM
  ! !===============================================================================
  !   !> Read parameters in param_RFM.ini
  ! !===============================================================================
  !   implicit none
  !   ! ----------------------------------------------------------------------------
  !   integer :: iblc
  !   logical :: iexist
  !   ! ----------------------------------------------------------------------------

  !   ! Read param_RFM.ini
  !   ! ===================
  !   inquire(file=trim(dirDATA)//'param_RFM.ini', exist=iexist)
  !   if (.not.iexist) then
  !      call mpistop('Paramfile param_RFM.ini does not exist!', 0)
  !   endif

  !   !=============================================================================
  !   open(30,file=trim('param_RFM.ini'))
  !   !=============================================================================
  !   read(30,*)! ==================================================================
  !   read(30,*)! RANDOM FOURIER MODES (RFM) PARAMETERS
  !   read(30,*)! ==================================================================
  !   read(30,*)! Main parameters for stochastic field
  !   read(30,*)! ==================================================================
  !   read(30,*) ! Integral length (for von-Karman spectrum)
  !   read(30,*) Lf_int
  !   delta = Lf_int/(sqrt(12.0_wp/5.0_wp)*0.74684_wp)
  !   read(30,*) ! Kolmogorov microscale [for dissip. range and bottleneck, estimated if put to 0.]
  !   read(30,*) xkkol ! Modified after as xkkol = 1.0_wp/xkkol
  !   read(30,*)! Turbulence intensity Tu
  !   read(30,*) Tu_RFM
  !   ampl_RFM = Tu_RFM*U_ref ! ampl_RFM = Tu*U_ref [m/s]
  !   read(30,*)! Injection zone imin: ndy_RFM nfy_rfm h_damp(1) h_damp(2)     (1 line per block)
  !   read(30,*)! ndy_RFM nfy_rfm: if 0 0 -> 1 ngy, if <0 -> is_RFM .false.
  !   read(30,*)! If wall: ndy_RFM such as yg(ndy_RFM)/d_99~1.0  and  h_damp~0.15*d_99
  !   read(30,*)! If non reflective BC: nfy_rfm<ngy-30 & outside sponge/streching zone  and  h_damp~2*dy
  !   h_damp = 0.0_wp
  !   do iblc=1,nbloc
  !      if (iblc==nob(iproc)) then
  !         read(30,*) ndy_RFM, nfy_rfm, h_damp(1), h_damp(2)
  !      else
  !         read(30,*)
  !      endif
  !   enddo
  !   if (ndy_RFM.eq.0) ndy_RFM=1
  !   if (nfy_rfm.eq.0) nfy_rfm=ngy
  !   if ((nfy_rfm.lt.0).or.(nfy_rfm.lt.0).or.(is_boundary(1,1).ge.0)) is_RFM=.false.
  !   if (nfy_rfm.gt.ngy) call mpistop("nfy_rfm can't be superior to ngy ~> param_RFM.ini to be modified...",1)
  !   if ((h_damp(1).lt.0.0_wp).or.(h_damp(2).lt.0.0_wp)) h_damp = 0.0_wp
  !   if (ndy_RFM.eq.1) h_damp(1)=0.0_wp
  !   if (nfy_rfm.eq.ngy) h_damp(2)=0.0_wp
  !   read(30,*)! ==================================================================
  !   read(30,*)! Choices for initialization of RFM orientation and discretization (if is_init_modes is T)
  !   read(30,*)! ==================================================================
  !   read(30,*)! Init. of fourier modes orientation/k-components and discretization
  !   read(30,*) is_init_modes
  !   read(30,*) ! is_xkm_given xkmin xkmax [only if is_init_modes, and calculated if is_xkm_given put to F]
  !   read(30,*) is_xkm_given, xkmin, xkmax
  !   read(30,*)! Number of RFM modes
  !   read(30,*) Nmode
  !   read(30,*)! wavenumber discretization: linear ('lin') or logarithmic ('log') distributions
  !   read(30,*) kdist
  !   read(30,*)! ==================================================================
  !   read(30,*)! Turbulence time evolution [unfrozen turbulence]
  !   read(30,*)! ==================================================================
  !   read(30,*)! Choice of method for turbulence time evolution
  !   read(30,*)! ['K': Kolmogorov time; 'H': Heisenberg time; 'N': none]
  !   read(30,*) time_turb
  !   ! read(30,*)! include mean flow convection
  !   ! read(30,*) is_convec
  !   read(30,*)! ==================================================================
  !   read(30,*)! Turbulence anisotropy
  !   read(30,*)! ==================================================================
  !   read(30,*)! Method to impose anisotropy
  !   read(30,*)! ['S': Smirnov et al. transformation; 'W': weighting function; 'N': none]
  !   read(30,*) anisotropy
  !   read(30,*)! database for TBL (Turbulent Boundary Layer) or CHAN (channel flow): base
  !   read(30,*)! [TBL: 'KTH'; 'LES'; 'Jim']
  !   read(30,*)! [CHAN: 'KMM'; 'VK']
  !   read(30,*) base
  !   read(30,*)! local Reynolds number
  !   read(30,*)! [based on momentum boundary layer thickness for TBL]
  !   read(30,*)! [based on half-width and u_tau for channel flow]
  !   read(30,*) Re_prof
  !   read(30,*)! ==================================================================
  !   read(30,*)! Additional choices
  !   read(30,*)! ==================================================================
  !   read(30,*)! check convergence of rms profiles
  !   read(30,*) is_check_conv
  !   read(30,*)! velocity field initialization
  !   read(30,*) is_field

  ! end subroutine read_param_RFM

  !===============================================================================
  module subroutine read_param_RFM
  !===============================================================================
    !> Read parameters in param_RFM.ini
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: iblc
    logical :: iexist
    ! ----------------------------------------------------------------------------
    ! Temp
    character(100) :: test_str

    ! Read param_RFM.ini
    ! ===================
    inquire(file=trim(dirDATA)//'param_RFM.ini', exist=iexist)
    if (.not.iexist) then
       call mpistop('Paramfile param_RFM.ini does not exist!', 0)
    endif

    !=============================================================================
    open(30,file=trim('param_RFM.ini'))
    !=============================================================================
    read(30,*) !=============================================================
    read(30,*) !=============================================================
    read(30,*) ! RANDOM FOURIER MODES (RFM): fill parameters
    read(30,*) !=============================================================
    read(30,*) !=============================================================
    read(30,*) ! RFM imposed only at boundary condition imin (temporary)
    read(30,*) ! 2 applications for RFM (both can be done at the same time):
    read(30,*) !   I/  Generate synthetic Freestream Turbulence (FST)
    read(30,*) !   II/ Generate Turbulent Boundary Layers (TBL)
    read(30,*) !
    read(30,*) !
    read(30,*) ! --> Number of RFM modes:
    read(30,*) !    ~> between 80 and 150, important computational overhead
    read(30,*) !       if Nmode is too large
    read(30,*) !    ~> parameter shared by the 2 applications
    read(30,*) !
    read(30,*) ! ------------------------------------------------------------
    read(30,*) ! I/ Freestream Turbulence parameters:
    read(30,*) ! ------------------------------------------------------------
    read(30,*) ! --> Integral length scale & FST intensity Tu
    read(30,*) ! --> Kolmogorov microscale
    read(30,*) !    ~> for dissipation range and bottleneck
    read(30,*) !    ~> directly estimated if put to 0
    read(30,*) ! ---------------------------
    read(30,*) ! --> Injection zone for FST:
    read(30,*) !    ~> Restriction possible only in direction j
    read(30,*) !    ~> "ndy_RFM nfy_RFM h_damp(1) h_damp(2)"
    read(30,*) !    ~> 1 line per block
    read(30,*) !    ~> If FST no activated in a block, specify "-1 -1 -1 -1"
    read(30,*) !    ~> If no injection restriction, specify "0 0 0 0"
    read(30,*) !    ~> If wall & laminar boundary layer, specify:
    read(30,*) !       _ ndy_RFM such as yg(ndy_RFM)/d_99 > 1.0
    read(30,*) !       _ h_damp ~ 0.15*d_99
    read(30,*) !    ~> If wall & turbulent boundary layer, specify:
    read(30,*) !       _ ndy_RFM=1,h_damp(1)=0 or nfy_RFM=ngy,h_damp(2)=0
    read(30,*) !    ~> If non reflective BC, specify:
    read(30,*) !       _ ndy_RFM > 20 & outside streching zone
    read(30,*) !       _ nfy_RFM < ngy-20 & outside streching zone
    read(30,*) !       _ h_damp ~ 2*delta_y
    read(30,*) ! ---------------------------
    read(30,*) ! --> Modes parametrization:
    read(30,*) ! Init. of modes, orientation/k-components, discretization
    read(30,*) !    a) Directly prescribed in vkm.bin file
    read(30,*) !      ~> put is_init_modes to F
    read(30,*) !    b) Calculated in solver with specified borns kmin/kmax
    read(30,*) !      ~> put is_init_modes to T & is_km_given to T
    read(30,*) !      ~> Give kmin & kmax in dimensional
    read(30,*) !    c) Calculated in solver without specifying kmin/kmax
    read(30,*) !      ~> put is_init_modes to T & is_km_given to F
    read(30,*) !      ~> user must verify dimensions of the domain and
    read(30,*) !         resolution at inlet are sufficient:
    read(30,*) !        _ Ly/Lz > =9*Lf & max(dx,dy,dz) ~ 0.25*kmax
    read(30,*) !      ~> kmin calculated as 2*pi/(9*Lf)
    read(30,*) !      ~> kmax calculated as 7*(2*pi/Lf)
    read(30,*) ! ---------------------------
    read(30,*) ! --> Turbulence time evolution:
    read(30,*) !    ~> hypothesis of unfrozen turbulence
    read(30,*) ! ---------------------------
    read(30,*) ! --> Turbulence anisotropy
    read(30,*) !
    read(30,*) ! ------------------------------------------------------------
    read(30,*) ! II/ Turbulent Boundary Layer parameters:
    read(30,*) ! ------------------------------------------------------------
    read(30,*) ! --> TBL database to be used:
    read(30,*) !     ~> KTH (only one implemented for the moment)
    read(30,*) !     ~> Personal DNS/LES
    read(30,*) ! --> Reynolds number:
    read(30,*) !    ~> based on momentum boundary layer thickness
    read(30,*) !    ~> 1 line per block, with jmin and jmax (ex: 1410 1000)
    read(30,*) !    ~> only applied if is_wall at jmin and/or jmax
    read(30,*) !    ~> if no TBL for block at jmin/jmax, put 0
    read(30,*) !
    read(30,*) ! ------------------------------------------------------------
    read(30,*) ! Additional choices
    read(30,*) ! ------------------------------------------------------------
    read(30,*) ! --> Check convergence of rms profiles
    read(30,*) !    ~> Preprocessing step to ensure RFM is properly
    read(30,*) !       generating the FST
    read(30,*) ! --> Velocity field initialization
    read(30,*) !    ~> For CHAN, all the field is initialized with turbulence
    read(30,*) !    ~> For other cases, only the first processor at imin
    read(30,*) !    ~> Help reduce the simulation transient
    read(30,*) !=============================================================
    read(30,*) !=============================================================
    read(30,*) ! Number of RFM modes: Nmode
    read(30,*) Nmode
    read(30,*) !=============================================================
    read(30,*) !                   I/ FREESTREAM TURBULENCE
    read(30,*) !=============================================================
    read(30,*) ! Activation of freestream turbulence: is_RFM_FST
    read(30,*) is_RFM_FST
    read(30,*) ! ---------------------------
    read(30,*) ! Main parameters
    read(30,*) ! ---------------------------
    read(30,*) ! Turbulence intensity Tu
    read(30,*) Tu_RFM
    read(30,*) ! Integral length scale [m]
    read(30,*) Lf_int
    delta = Lf_int/(sqrt(12.0_wp/5.0_wp)*0.74684_wp)
    read(30,*) ! Kolmogorov microscale [m] (estimated if put to 0)
    read(30,*) xkkol
    read(30,*) ! ---------------------------
    read(30,*) ! Injection zone
    read(30,*) ! ---------------------------
    read(30,*) ! ndy_RFM nfy_RFM h_damp(1) h_damp(2)    # 1 block per line
    ndy_RFM=0; nfy_RFM=0; h_damp = 0.0_wp
    do iblc=1,nbloc
       if (iblc==nob(iproc)) then
          read(30,*) ndy_RFM, nfy_RFM, h_damp(1), h_damp(2)
       else
          read(30,*)
       endif
    enddo
    read(30,*) ! ---------------------------
    read(30,*) ! Modes parametrization
    read(30,*) ! ---------------------------
    read(30,*) ! is_init_modes
    read(30,*) is_init_modes
    read(30,*) ! is_km_given kmin [/m] kmax [/m] (only if is_init_modes)
    read(30,*) is_xkm_given, xkmin, xkmax
    read(30,*) ! Discretization: linear ('lin') or logarithmic ('log')
    read(30,*) kdist
    read(30,*) ! ---------------------------
    read(30,*) ! Turbulence time evolution
    read(30,*) ! ---------------------------
    read(30,*) ! Choice of method for turbulence time evolution: time_turb
    read(30,*) ! ['K': Kolmogorov time; 'H': Heisenberg time; 'N': none]
    read(30,*) time_turb
    read(30,*) ! ---------------------------
    read(30,*) ! Turbulence anisotropy
    read(30,*) ! ---------------------------
    read(30,*) ! Method to impose anisotropy
    read(30,*) ! ['S': Smirnov et al. transformation; 'W': simple weighting function; 'N': none]
    read(30,*) anisotropy
    read(30,*) ! Database for CHAN (channel flow): base
    read(30,*) ! [CHAN: 'MKM'; 'V&K']        **** only MKM implemented ****
    read(30,*) base_FST
    read(30,*) ! local Reynolds number: Re_prof_FST
    read(30,*) ! [based on half-width and u_tau]
    read(30,*) Re_prof_FST
    read(30,*) !=============================================================
    read(30,*) !                II/ Turbulent Boundary Layers
    read(30,*) !=============================================================
    read(30,*) ! Activation of turbulent boundary layers: is_RFM_TBL
    read(30,*) is_RFM_TBL(1)
    is_RFM_TBL(2)=is_RFM_TBL(1)
    read(30,*) ! Turbulence intensity Tu in TBL
    read(30,*) Tu_TBL
    read(30,*) ! Database for TBL: base
    read(30,*) ! ['KTH'; 'LES'; 'Jim']  **** only KTH implemented ****
    read(30,*) base_TBL
    read(30,*) ! Local Reynolds number
    read(30,*) ! Re_prof_TBL at jmin and jmax            # 1 block per line
    Re_prof_TBL=0.0_wp
    do iblc=1,nbloc
       if (iblc==nob(iproc)) then
          read(30,*) Re_prof_TBL(1), Re_prof_TBL(2)
       else
          read(30,*)
       endif
    enddo
    read(30,*) !=============================================================
    read(30,*) !                     ADDITIONNAL CHOICES
    read(30,*) !=============================================================
    read(30,*) ! check convergence of rms profiles: is_check_conv
    read(30,*) is_check_conv
    read(30,*) ! velocity field initialization: is_field
    read(30,*) is_field

    ! Verifications
    ! -------------
    ! No injection zone
    if (ndy_RFM.eq.0) ndy_RFM=1
    if (nfy_RFM.eq.0) nfy_RFM=ngy
    ! No FST injected in the block
    if (bl(nob(iproc))%BC(1).ge.0) is_RFM=.false.
    if ((nfy_RFM.lt.0).or.(nfy_RFM.lt.0)) is_RFM=.false.
    ! Bad specification of injection zone
    if (nfy_RFM.gt.ngy) call mpistop("ndy_RFM can't be superior to ngy ~> param_RFM.ini to be modified...",1)
    ! Check on damping height
    if ((h_damp(1).lt.0.0_wp).or.(h_damp(2).lt.0.0_wp)) h_damp = 0.0_wp
    if (ndy_RFM.eq.1) h_damp(1)=0.0_wp
    if (nfy_RFM.eq.ngy) h_damp(2)=0.0_wp
    ! For the moment, anisotropy in FST mode only available for CHAN
    ! If is_RFM_FST and not CHAN, anisotropy must be put to F
    ! Anisotropy other than based on BL needs to be implemented
    if ((.not.CHAN).and.(anisotropy.ne."N")) call mpistop("Anistropy can only be imposed on CHAN for the moment",0)
    ! No damping if is_RFM_TBL ~> the Reynolds stresses naturally go to 0 in the FST
    if (.not.is_RFM_FST) is_damping=.false.
    ! Check BC for which is_RFM_TBL is .true.
    if ((is_RFM_TBL(1)).and.(Re_prof_TBL(1).le.0)) is_RFM_TBL(1)=.false.
    if ((is_RFM_TBL(2)).and.(Re_prof_TBL(2).le.0)) is_RFM_TBL(2)=.false.
    ! Check on TBL specified
    if ((coord(2).eq.0).and.(.not.is_bc_wall(2,1)).and.(is_RFM_TBL(1))) &
         call mpistop("TBL cannot be applied at jmin in block "//trim(numchar(nob(iproc)))//" because it is not a wall",1)
    if ((coord(2).eq.ndomy-1).and.(.not.is_bc_wall(2,2)).and.(is_RFM_TBL(2))) &
         call mpistop("TBL cannot be applied at jmax in block "//trim(numchar(nob(iproc)))//" because it is not a wall",1)

  end subroutine read_param_RFM

  !===============================================================================
  module subroutine init_RFM_planes
  !===============================================================================
    !> routine to detect the number plane of the inlet planes
  !===============================================================================
    use mod_io_snapshots
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: ip
    ! ----------------------------------------------------------------------------

    inlet_ip(:) = 0

    do ip=1,nsnapshots
       if ((snapshots(ip)%normal.eq.1).and.(snapshots(ip)%nvar==3).and.(snapshots(ip)%index<6) &
           .and.(.not.snapshots(ip)%stamp)) then
          if ((snapshots(ip)%var(1)=='udf')&
             .and.(snapshots(ip)%var(2)=='udf')&
             .and.(snapshots(ip)%var(3)=='udf')) inlet_ip(snapshots(ip)%index) = ip
       endif
    enddo

  end subroutine init_RFM_planes

end submodule smod_RFM_init
