!=================================================================================
module mod_init_flow
!=================================================================================
  !> Module to initialize flow
!=================================================================================
  use mod_mpi
  use mod_flow
  use mod_constant
  implicit none
  !-------------------------------------------------------------------------------
  ! Laminar boundary layer & Blasius solution
  logical  :: is_LBL,is_similarity
  !-------------------------------------------------------------------------------
  ! Pulse parameters
  logical :: is_pulse
  real(wp) :: ampl_pulse,b_pulse
  real(wp) :: x_pulse,y_pulse,z_pulse
  !-------------------------------------------------------------------------------
  ! Vortex parameters
  logical :: is_vortex
  integer :: type_vortex
  real(wp) :: x_vortex,y_vortex,z_vortex
  !-------------------------------------------------------------------------------

contains
  
  !===============================================================================
  subroutine init_flow
  !===============================================================================
    !> Initialization of flow variables
  !===============================================================================
    use mod_eos
    use mod_RFM
    use mod_flow0
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k
    ! ----------------------------------------------------------------------------

    ! Initialization of thermodynamic variables to reference values
    ! =============================================================
    rho= rho_ref
    prs= p_ref
    Tmp= T_ref
    if (ACT) then
       prs= 100000.0_wp!p_ref
       Tmp= 288.0_wp!T_ref
       rho= prs/(rg*Tmp)
    endif

    ! First initialization of velocity components
    ! ===========================================
    uu = 0.0_wp
    vv = 0.0_wp
    ww = 0.0_wp

    ! General initialization
    ! ----------------------
    if (is_RFM) call init_RFM

    ! Particular initializations
    ! ==========================
    if (TGV) call add_tgv   ! Taylor-Green-Vortex
    
    ! if not is_RFM, random velocities in channel flow
    ! if is_RFM, mean field added here, based on database
    if (CHAN) call init_vel_chan
    
    if (PHILL) call init_vel_phill ! random velocities in periodic hill flow
    
    if (CYL) call init_vel_cyl ! potential flow past cylinder

    if (TURB) call init_turbine ! turbine flow

    if ((SRC).or.(CAV).or.(SHIT)) then
       ! constant u
       uu= Uc1*u_ref
       vv= Uc2*u_ref
       ww= Uc3*u_ref
      !call add_shock
    endif

    ! cavity flow (4 block case)
    if ((CAV).and.(iorder_visc.ne.0)) call add_BL_over_cavity
    ! if (STBL) call add_lam_BL
    if (STBL) then
       ! if (is_similarity) then
       !     call add_lam_BL
       !  else
       !     call add_turb_BL
       !     ! call add_lam_BL
       !  endif
       !call add_BL
        uu= Uc1*u_ref
        vv= Uc2*u_ref
        ww= Uc3*u_ref
     endif
    if(T3C)then 
        uu= Uc1*u_ref
        vv= Uc2*u_ref
        ww= Uc3*u_ref

        ! Damp velocity near the wall
          ! ===========================
          ! (to avoid violation of non-penetrability condition)
        !if (is_bc_wall(2,1)) then
        !    do k=1,nz
        !       do j=1,10
        !          do i=1,nx
        !             uu(i,j,k)=uu(i,j,k)/10.0_wp*dble(j-1)
        !             vv(i,j,k)=vv(i,j,k)/10.0_wp*dble(j-1)
        !          enddo
        !       enddo
        !    enddo
        ! endif
        ! if (is_bc_wall(2,2)) then
        !    do k=1,nz
        !       do j=ny-9,ny
        !          do i=1,nx
        !             uu(i,j,k)=uu(i,j,k)/10.0_wp*dble(ny-j)
        !             vv(i,j,k)=vv(i,j,k)/10.0_wp*dble(ny-j)
        !          enddo
        !       enddo
        !    enddo
        ! endif

    endif
    ! actuator case
    if (ACT) call init_vel_act

    ! LE configuration  /!\ Completely obsolete, LE case needs to be reimplemented
    if (LE) then
       ! Initialization of the field based on RANS simu
       call init_vel_vane
       ! Addition of RFM close to the inlet
       !!!if ((BC_face(1,1)%sort.lt.0).and.(is_RFM)) call init_vel_RFM ! Not here anymore -> in init_RFM
    endif

    ! Trailing-edge configuration
    if (TE) then
       ! Initialization of the field based on full blade simu or from scratch
       call init_vel_TE
       ! return before calculating rho, rhou, ... again
       return
    endif

    ! Initialize mean primitive variables [to avoid using NaN]
    ! ===================================
    if (is_mean0) then
       rho0=rho
       u0=uu
       v0=vv
       w0=ww
       p0=prs
       T0=Tmp
    endif

    ! Add pulse
    ! =========
    if (is_pulse) call add_pulse

    ! Add vortex
    ! ==========
    if (is_vortex) call add_vortex

    ! Initialization of conservative variables from primitive variables
    ! =================================================================
    ! momenta
    rhou = rho*uu
    rhov = rho*vv
    rhow = rho*ww
    ! total energy
    do k=ndzt,nfzt
       do j=ndyt,nfyt
          do i=ndxt,nfxt
             rhoe(i,j,k)= rho(i,j,k)*(ecalc_tro(Tmp(i,j,k),rho(i,j,k)) &
                        + 0.5_wp*(uu(i,j,k)**2+vv(i,j,k)**2+ww(i,j,k)**2))
          enddo
       enddo
    enddo
    
  end subroutine init_flow


  !===============================================================================
  subroutine init_thermo
  !===============================================================================
    !> Generic initialization of primitive variables
    !  [needed to be sure that unused ghost cells are filled once]
  !===============================================================================
    use mod_eos
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k
    ! ----------------------------------------------------------------------------

    ! Initialization of arrays to reference values
    ! ============================================
    rho= rho_ref
    prs= p_ref
    Tmp= T_ref
    do k=ndzt,nfzt
       do j=ndyt,nfyt
          do i=ndxt,nfxt
             rhoe(i,j,k)= rho(i,j,k)*ecalc_tro(Tmp(i,j,k),rho(i,j,k))
          enddo
       enddo
    enddo

  end subroutine init_thermo

  !=============================================================================
  subroutine init_BC_ref
  !=============================================================================
    !> initialization of reference field for BCs
  !=============================================================================
    use mod_utils
    use mod_eos
    use mod_block      !temporary
    use mod_mpi_part   !temporary
    use mod_io_restart_BCref
    implicit none
    !---------------------------------------------------------------------------
    integer :: l,l1,l2,ng1_,ng2_,n1_,n2_,i1,i2,i3,ig1,ig2,d1,d2,nproc_l1
    integer :: fh
    logical :: iexist
    integer, dimension(MPI_STATUS_SIZE) :: statut
    character(len=31) :: ref_name='BC_profiles/bc_prof_imin_bl.bin'
    character(:), allocatable :: ref_file
    real(wp) :: Temp ! temperature
    real(wp), dimension(:,:), allocatable :: ro_ref,rou_ref,rov_ref,row_ref,roe_ref
    real(wp), dimension(:,:), allocatable :: c_ref,prs_ref ! for T&D
    real(wp), dimension(:,:), allocatable :: H0_ref,u_ref,v_ref,w_ref,s_ref
    !---------------------------------------------------------------------------

    is_BCref_init = .true.

    ! Check if reference values has been saved previously
    ! ===================================================
    ! Not read if is_read_ref is true
    ref_file='restart_BCref.bin'
    inquire(file=ref_file, exist=iexist)
    if ((is_read_ref).or.(iexist)) is_BCref_init=.false.
    if ((idepart.eq.1).and.(.not.is_read_ref)) is_BCref_init=.true.
    if ((is_init_2D3D).and.(.not.is_read_ref)) then
       is_BCref_init=.true.
       if ((iexist).and.(iproc.eq.0)) call system('rm '//trim(ref_file))
    endif

    ! Display type of application for reference values
    ! ================================================
    if (iproc.eq.0) then
       print *,repeat('=',70)
       print *,'Init reference values in BC (init_BC_ref)'
       print *,repeat('=',70)
       if (is_read_ref) then
          print *,'Reference values imposed from a BC file:'
       else if (is_BCref_init) then
          print *,'Reference values initialized from field'
       else
          print *,'Reference values read from previous run: restart_BCref.bin'
       endif
    endif

    ! Loop on each faces
    ! ==================
    do l=1,6
       
       ! Initialization of BC reference values on profiles
       ! -------------------------------------------------
       ! face 1 (imin): (l1,l2)==(1,1)
       ! face 2 (imax): (l1,l2)==(1,2)
       ! face 3 (jmin): (l1,l2)==(2,1)
       ! face 4 (jmax): (l1,l2)==(2,2)
       ! face 5 (kmin): (l1,l2)==(3,1)
       ! face 6 (kmax): (l1,l2)==(3,2)
       l1 = (l+1)/2
       l2 = l-(l1-1)*2

       if (l1.eq.1) then
          ng1_=ngy; ng2_=ngz; n1_=ny; n2_=nz; d1=2; d2=3; nproc_l1=ndomx
       elseif (l1.eq.2) then
          ng1_=ngx; ng2_=ngz; n1_=nx; n2_=nz; d1=1; d2=3; nproc_l1=ndomy
       elseif (l1.eq.3) then
          ng1_=ngx; ng2_=ngy; n1_=nx; n2_=ny; d1=1; d2=2; nproc_l1=ndomz
       endif

       ! Read from reference BC profile files
       ! ------------------------------------
       if (is_read_ref) then
          ! Reference file name
          ref_file = trim(ref_name)
          if (l.eq.1) then
             ref_file(21:24) = 'imin'
          else if (l.eq.2) then
             ref_file(21:24) = 'imax'
          else if (l.eq.3) then
             ref_file(21:24) = 'jmin'
          else if (l.eq.4) then
             ref_file(21:24) = 'jmax'
          else if (l.eq.5) then
             ref_file(21:24) = 'kmin'
          else if (l.eq.6) then
             ref_file(21:24) = 'kmax'
          endif
          ref_file = ref_file(1:27)//trim(numchar(nob(iproc)))//'.bin'

          ! Verification of existence (needed if BC_face(l1,l2)%is_mean_ref)
          inquire(file=ref_file, exist=iexist)
          if ((.not.iexist).and.(BC_face(l1,l2)%is_mean_ref)) then
             call mpistop('Reference profile '//ref_file//' does not exist!',1)
          endif

          ! Read from a reference file
          ! --------------------------
          if (iexist) then
             ! Each proc of the block open the file
             call MPI_FILE_OPEN(COMM_intrablock,ref_file, &
               MPI_MODE_RDONLY, MPI_INFO_NULL,fh,info)

             ! Only the proc concerned allocate ref variables
             if (BC_face(l1,l2)%is_mean_ref) then
                ! Indication of what is read
                if ((coord(l1).eq.0).or.(coord(l1).eq.nproc_l1-1)) &
                     print *,"     Block ",trim(numchar(nob(iproc)))," ~> reading ref. for BC ",ref_file(21:24)

                ! If turbomachine condition, only 5 reference variables needed
                if (BC_face(l1,l2)%sort==-4) then
                   ! Temporary arrays for the ref bc block
                   allocate( H0_ref(ng1_,ng2_),u_ref(ng1_,ng2_), &
                             v_ref(ng1_,ng2_),w_ref(ng1_,ng2_), &
                             s_ref(ng1_,ng2_))
                   allocate(BC_face(l1,l2)%Uref(1:n1_,1:n2_,5))

                   ! Reading variables
                   call MPI_FILE_READ(fh,H0_ref,size(H0_ref),&
                        MPI_DOUBLE_PRECISION,statut,info)
                   call MPI_FILE_READ(fh,u_ref,size(u_ref),&
                        MPI_DOUBLE_PRECISION,statut,info)
                   call MPI_FILE_READ(fh,v_ref,size(v_ref),&
                        MPI_DOUBLE_PRECISION,statut,info)
                   call MPI_FILE_READ(fh,w_ref,size(w_ref),&
                        MPI_DOUBLE_PRECISION,statut,info)
                   call MPI_FILE_READ(fh,s_ref,size(s_ref),&
                        MPI_DOUBLE_PRECISION,statut,info)

                else
                   ! Temporary arrays for the ref bc block
                   allocate( ro_ref(ng1_,ng2_),rou_ref(ng1_,ng2_), &
                             rov_ref(ng1_,ng2_),row_ref(ng1_,ng2_), &
                             prs_ref(ng1_,ng2_),roe_ref(ng1_,ng2_),c_ref(ng1_,ng2_))
                   allocate(BC_face(l1,l2)%Uref(1:n1_,1:n2_,10))

                   ! Reading variables
                   call MPI_FILE_READ(fh,ro_ref,size(ro_ref),&
                        MPI_DOUBLE_PRECISION,statut,info)
                   call MPI_FILE_READ(fh,rou_ref,size(rou_ref),&
                        MPI_DOUBLE_PRECISION,statut,info)
                   call MPI_FILE_READ(fh,rov_ref,size(rov_ref),&
                        MPI_DOUBLE_PRECISION,statut,info)
                   call MPI_FILE_READ(fh,row_ref,size(row_ref),&
                        MPI_DOUBLE_PRECISION,statut,info)
                   call MPI_FILE_READ(fh,roe_ref,size(roe_ref),&
                        MPI_DOUBLE_PRECISION,statut,info)
                   call MPI_FILE_READ(fh,prs_ref,size(prs_ref),&
                        MPI_DOUBLE_PRECISION,statut,info)
                   call MPI_FILE_READ(fh,c_ref,size(c_ref),&
                        MPI_DOUBLE_PRECISION,statut,info)
                endif

             endif

             ! Each proc of the block close the file
             call MPI_FILE_CLOSE(fh,info)

          endif

          if (BC_face(l1,l2)%is_mean_ref) then
             
             ! Filling reference values
             if ((BC_face(l1,l2)%sort==-4).and.(iexist)) then
                do i1=1,n1_
                   do i2=1,n2_
                      ig1 = i1 + coord(d1)*n1_
                      ig2 = i2 + coord(d2)*n2_
                      BC_face(l1,l2)%Uref(i1,i2,1) =H0_ref(ig1,ig2)
                      BC_face(l1,l2)%Uref(i1,i2,2) = u_ref(ig1,ig2)
                      BC_face(l1,l2)%Uref(i1,i2,3) = v_ref(ig1,ig2)
                      BC_face(l1,l2)%Uref(i1,i2,4) = w_ref(ig1,ig2)
                      BC_face(l1,l2)%Uref(i1,i2,5) = s_ref(ig1,ig2)
                   enddo
                enddo

                deallocate(H0_ref,u_ref,v_ref,w_ref,s_ref)
             elseif (iexist) then
                do i1=1,n1_
                   do i2=1,n2_
                      ig1 = i1 + coord(d1)*n1_
                      ig2 = i2 + coord(d2)*n2_
                      BC_face(l1,l2)%Uref(i1,i2,1) =  ro_ref(ig1,ig2)
                      BC_face(l1,l2)%Uref(i1,i2,2) = rou_ref(ig1,ig2)/ro_ref(ig1,ig2)
                      BC_face(l1,l2)%Uref(i1,i2,3) = rov_ref(ig1,ig2)/ro_ref(ig1,ig2)
                      BC_face(l1,l2)%Uref(i1,i2,4) = row_ref(ig1,ig2)/ro_ref(ig1,ig2)
                      BC_face(l1,l2)%Uref(i1,i2,5) = prs_ref(ig1,ig2)
                      BC_face(l1,l2)%Uref(i1,i2,6) = c_ref(ig1,ig2)**2
                      BC_face(l1,l2)%Uref(i1,i2,7) = rou_ref(ig1,ig2)
                      BC_face(l1,l2)%Uref(i1,i2,8) = rov_ref(ig1,ig2)
                      BC_face(l1,l2)%Uref(i1,i2,9) = row_ref(ig1,ig2)
                      BC_face(l1,l2)%Uref(i1,i2,10)= roe_ref(ig1,ig2)
                   enddo
                enddo
                deallocate(ro_ref,rou_ref,rov_ref,row_ref,prs_ref,roe_ref,c_ref)
             else
                call mpistop('Reference profiles for BCs do not exist !',1)
             endif
          endif

       ! Based on field values
       ! ---------------------
       elseif (is_BCref_init) then
          if (BC_face(l1,l2)%is_mean_ref) then
             ! BC_face%Uref is allocated here <- TO BE CHANGED ?
             allocate(BC_face(l1,l2)%Uref(1:n1_,1:n2_,10))

             if (l1.eq.1) then
                i3 = 1
                if (l2.eq.2) i3=nx
                ! Filling reference values
                if ((idepart.eq.FROM_SCRATCH).or.(.not.is_bc_TD(l1,l2))) then
                   do i1=1,n1_
                      do i2=1,n2_
                         BC_face(l1,l2)%Uref(i1,i2,1) = rho(i3,i1,i2)
                         BC_face(l1,l2)%Uref(i1,i2,2) = rhou(i3,i1,i2)/rho(i3,i1,i2)
                         BC_face(l1,l2)%Uref(i1,i2,3) = rhov(i3,i1,i2)/rho(i3,i1,i2)
                         BC_face(l1,l2)%Uref(i1,i2,4) = rhow(i3,i1,i2)/rho(i3,i1,i2)
                         BC_face(l1,l2)%Uref(i1,i2,7) = rhou(i3,i1,i2)
                         BC_face(l1,l2)%Uref(i1,i2,8) = rhov(i3,i1,i2)
                         BC_face(l1,l2)%Uref(i1,i2,9) = rhow(i3,i1,i2)
                         BC_face(l1,l2)%Uref(i1,i2,10) =rhoe(i3,i1,i2)
                      enddo
                   enddo
                   if ((is_similarity).and.(.not.is_curv).and.(.not.is_curv3)) then
                      ! c2 and P_ref based on reference values
                      do i1=1,n1_
                         do i2=1,n2_
                            BC_face(l1,l2)%Uref(i1,i2,5) = p_ref
                            BC_face(l1,l2)%Uref(i1,i2,6) = c2calc_tro(T_ref,rho_ref)
                            !!XG BC_face(l1,l2)%Uref(i1,i2,5) = pcalc_roero(rhoe(i3,i1,i2),rho(i3,i1,i2),Tmp(i3,i1,i2))
                            !!XG Temp=tcalc_roero(rhoe(i3,i1,i2),rho(i3,i1,i2),Tmp(i3,i1,i2))
                            !!XG BC_face(l1,l2)%Uref(i1,i2,6) = Temp
                         enddo
                      enddo
                   else
                      ! c2 and P_ref based on initial values
                      do i1=1,n1_
                         do i2=1,n2_
                            BC_face(l1,l2)%Uref(i1,i2,5) = pcalc_roero(rhoe(i3,i1,i2),rho(i3,i1,i2),Tmp(i3,i1,i2))
                            Temp=tcalc_roero(rhoe(i3,i1,i2),rho(i3,i1,i2),Tmp(i3,i1,i2))
                            !BC_face(l1,l2)%Uref(i1,i2,6)= c2calc_tro(Temp,rho(i3,i1,i2))
                            BC_face(l1,l2)%Uref(i1,i2,6)= Temp
                         enddo
                      enddo
                   endif
                else
                   do i1=1,n1_
                      do i2=1,n2_
                         BC_face(l1,l2)%Uref(i1,i2,1) = BC_face(l1,l2)%U0(1,i1,i2,1)
                         BC_face(l1,l2)%Uref(i1,i2,2) = BC_face(l1,l2)%U0(1,i1,i2,2)
                         BC_face(l1,l2)%Uref(i1,i2,3) = BC_face(l1,l2)%U0(1,i1,i2,3)
                         BC_face(l1,l2)%Uref(i1,i2,4) = BC_face(l1,l2)%U0(1,i1,i2,4)
                         BC_face(l1,l2)%Uref(i1,i2,5) = BC_face(l1,l2)%U0(1,i1,i2,5)
                         BC_face(l1,l2)%Uref(i1,i2,6) = BC_face(l1,l2)%U0(1,i1,i2,6)
                         BC_face(l1,l2)%Uref(i1,i2,7) = BC_face(l1,l2)%U0(1,i1,i2,1)*BC_face(l1,l2)%U0(1,i1,i2,2)
                         BC_face(l1,l2)%Uref(i1,i2,8) = BC_face(l1,l2)%U0(1,i1,i2,1)*BC_face(l1,l2)%U0(1,i1,i2,3)
                         BC_face(l1,l2)%Uref(i1,i2,9) = BC_face(l1,l2)%U0(1,i1,i2,1)*BC_face(l1,l2)%U0(1,i1,i2,4)
                         BC_face(l1,l2)%Uref(i1,i2,10)= BC_face(l1,l2)%U0(1,i1,i2,1) &
                              *ecalc_pro(BC_face(l1,l2)%U0(1,i1,i2,6),BC_face(l1,l2)%U0(1,i1,i2,1),T_ref)
                      enddo
                   enddo
                endif

             elseif (l1.eq.2) then
                i3 = 1
                if (l2.eq.2) i3=ny
                ! Filling reference values
                if ((idepart.eq.FROM_SCRATCH).or.(.not.is_bc_TD(l1,l2))) then
                   do i1=1,n1_
                      do i2=1,n2_
                         BC_face(l1,l2)%Uref(i1,i2,1) = rho(i1,i3,i2)
                         BC_face(l1,l2)%Uref(i1,i2,2) = rhou(i1,i3,i2)/rho(i1,i3,i2)
                         BC_face(l1,l2)%Uref(i1,i2,3) = rhov(i1,i3,i2)/rho(i1,i3,i2)
                         BC_face(l1,l2)%Uref(i1,i2,4) = rhow(i1,i3,i2)/rho(i1,i3,i2)
                         BC_face(l1,l2)%Uref(i1,i2,5) = pcalc_roero(rhoe(i1,i3,i2),rho(i1,i3,i2),Tmp(i3,i1,i2))
                         BC_face(l1,l2)%Uref(i1,i2,6) = c2calc_tro( tcalc_roero(rhoe(i1,i3,i2),rho(i1,i3,i2),Tmp(i1,i3,i2)) , rho(i1,i3,i2))
                         ! BC_face(l1,l2)%Uref(i1,i2,6) = c2calc_tro(Tmp(i1,i3,i2),rho(i1,i3,i2)) XG WHY ????
                         BC_face(l1,l2)%Uref(i1,i2,7) = rhou(i1,i3,i2)
                         BC_face(l1,l2)%Uref(i1,i2,8) = rhov(i1,i3,i2)
                         BC_face(l1,l2)%Uref(i1,i2,9) = rhow(i1,i3,i2)
                         BC_face(l1,l2)%Uref(i1,i2,10) =rhoe(i1,i3,i2)
                      enddo
                   enddo
                else
                   do i1=1,n1_
                      do i2=1,n2_
                         BC_face(l1,l2)%Uref(i1,i2,1) = BC_face(l1,l2)%U0(1,i1,i2,1)
                         BC_face(l1,l2)%Uref(i1,i2,2) = BC_face(l1,l2)%U0(1,i1,i2,2)
                         BC_face(l1,l2)%Uref(i1,i2,3) = BC_face(l1,l2)%U0(1,i1,i2,3)
                         BC_face(l1,l2)%Uref(i1,i2,4) = BC_face(l1,l2)%U0(1,i1,i2,4)
                         BC_face(l1,l2)%Uref(i1,i2,5) = BC_face(l1,l2)%U0(1,i1,i2,5)
                         BC_face(l1,l2)%Uref(i1,i2,6) = BC_face(l1,l2)%U0(1,i1,i2,6)
                         BC_face(l1,l2)%Uref(i1,i2,7) = BC_face(l1,l2)%U0(1,i1,i2,1)*BC_face(l1,l2)%U0(1,i1,i2,2)
                         BC_face(l1,l2)%Uref(i1,i2,8) = BC_face(l1,l2)%U0(1,i1,i2,1)*BC_face(l1,l2)%U0(1,i1,i2,3)
                         BC_face(l1,l2)%Uref(i1,i2,9) = BC_face(l1,l2)%U0(1,i1,i2,1)*BC_face(l1,l2)%U0(1,i1,i2,4)
                         BC_face(l1,l2)%Uref(i1,i2,10)= BC_face(l1,l2)%U0(1,i1,i2,1)*ecalc_pro(BC_face(l1,l2)%U0(1,i1,i2,6),BC_face(l1,l2)%U0(1,i1,i2,1),T_ref)
                      enddo
                   enddo
                endif
             else if (l1.eq.3) then
                i3 = 1
                if (l2.eq.2) i3=nz
                ! Filling reference values
                if ((idepart.eq.FROM_SCRATCH).or.(.not.is_bc_TD(l1,l2))) then
                   do i1=1,n1_
                      do i2=1,n2_
                         BC_face(l1,l2)%Uref(i1,i2,1) = rho(i1,i2,i3)
                         BC_face(l1,l2)%Uref(i1,i2,2) = rhou(i1,i2,i3)/rho(i1,i2,i3)
                         BC_face(l1,l2)%Uref(i1,i2,3) = rhov(i1,i2,i3)/rho(i1,i2,i3)
                         BC_face(l1,l2)%Uref(i1,i2,4) = rhow(i1,i2,i3)/rho(i1,i2,i3)
                         BC_face(l1,l2)%Uref(i1,i2,5) = pcalc_roero(rhoe(i1,i2,i3),rho(i1,i2,i3),Tmp(i3,i1,i2))
                         BC_face(l1,l2)%Uref(i1,i2,6) = c2calc_tro( tcalc_roero(rhoe(i1,i2,i3),rho(i1,i2,i3),Tmp(i1,i2,i3)) , rho(i1,i2,i3))
                         ! BC_face(l1,l2)%Uref(i1,i2,6) = c2calc_tro(Tmp(i1,i2,i3),rho(i1,i2,i3))
                         BC_face(l1,l2)%Uref(i1,i2,7) = rhou(i1,i2,i3)
                         BC_face(l1,l2)%Uref(i1,i2,8) = rhov(i1,i2,i3)
                         BC_face(l1,l2)%Uref(i1,i2,9) = rhow(i1,i2,i3)
                         BC_face(l1,l2)%Uref(i1,i2,10) =rhoe(i1,i2,i3)
                      enddo
                   enddo
                else
                   do i1=1,n1_
                      do i2=1,n2_
                         BC_face(l1,l2)%Uref(i1,i2,1) = BC_face(l1,l2)%U0(1,i1,i2,1)
                         BC_face(l1,l2)%Uref(i1,i2,2) = BC_face(l1,l2)%U0(1,i1,i2,2)
                         BC_face(l1,l2)%Uref(i1,i2,3) = BC_face(l1,l2)%U0(1,i1,i2,3)
                         BC_face(l1,l2)%Uref(i1,i2,4) = BC_face(l1,l2)%U0(1,i1,i2,4)
                         BC_face(l1,l2)%Uref(i1,i2,5) = BC_face(l1,l2)%U0(1,i1,i2,5)
                         BC_face(l1,l2)%Uref(i1,i2,6) = BC_face(l1,l2)%U0(1,i1,i2,6)
                         BC_face(l1,l2)%Uref(i1,i2,7) = BC_face(l1,l2)%U0(1,i1,i2,1)*BC_face(l1,l2)%U0(1,i1,i2,2)
                         BC_face(l1,l2)%Uref(i1,i2,8) = BC_face(l1,l2)%U0(1,i1,i2,1)*BC_face(l1,l2)%U0(1,i1,i2,3)
                         BC_face(l1,l2)%Uref(i1,i2,9) = BC_face(l1,l2)%U0(1,i1,i2,1)*BC_face(l1,l2)%U0(1,i1,i2,4)
                         BC_face(l1,l2)%Uref(i1,i2,10)= BC_face(l1,l2)%U0(1,i1,i2,1)*ecalc_pro(BC_face(l1,l2)%U0(1,i1,i2,6),BC_face(l1,l2)%U0(1,i1,i2,1),T_ref)
                      enddo
                   enddo
                endif
             endif

          endif

       ! Read from restart_BCref.bin file
       ! --------------------------------
       else
          if (BC_face(l1,l2)%is_mean_ref) then
             ! BC_face%Uref is allocated here <- TO BE CHANGED ?
             allocate(BC_face(l1,l2)%Uref(1:n1_,1:n2_,10))
          endif
       endif
    enddo

    ! Everybody read at the same time
    ! -------------------------------
    if ((.not.is_read_ref).and.(.not.is_BCref_init)) call read_restart_BCref('restart_BCref.bin')

    if (iproc.eq.0) print *,"======================================================================"

  end subroutine init_BC_ref

  !===============================================================================
  subroutine add_tgv
  !===============================================================================
    !> Initialization for Taylor-Green Vortex (TGV)
  !===============================================================================
    use mod_eos
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    ! ---------------------------------------------------------------------------
    
    ! Incompressible TGV field
    if (is_curv) then
       do k=1,nz
          do j=1,ny
             do i=1,nx
                uu(i,j,k) =  u_ref*sin(xc(i,j))*cos(yc(i,j))*cos(z(k))
                vv(i,j,k) = -u_ref*cos(xc(i,j))*sin(yc(i,j))*cos(z(k))
                ww(i,j,k) = 0.0_wp
                prs(i,j,k)= p_ref + rho_ref*u_ref**2/16.0_wp*(cos(2.0_wp*z(k))+ 2.0_wp) &
                                                            *(cos(2.0_wp*xc(i,j))+cos(2.0_wp*yc(i,j)))

                Tmp(i,j,k)= prs(i,j,k)/p_ref*T_ref
                rho(i,j,k)= rocalc_pt(prs(i,j,k),Tmp(i,j,k),Tmp(i,j,k))
             enddo
          enddo
       enddo
    else
       do k=1,nz
          do j=1,ny
             do i=1,nx
                uu(i,j,k) =  u_ref*sin(x(i))*cos(y(j))*cos(z(k))
                vv(i,j,k) = -u_ref*cos(x(i))*sin(y(j))*cos(z(k))
                ww(i,j,k) = 0.0_wp
                prs(i,j,k)= p_ref + rho_ref*u_ref**2/16.0_wp*(cos(2.0_wp*z(k))+ 2.0_wp) &
                                                            *(cos(2.0_wp*x(i))+cos(2.0_wp*y(j)))

                Tmp(i,j,k)= prs(i,j,k)/p_ref*T_ref
                rho(i,j,k)= rocalc_pt(prs(i,j,k),Tmp(i,j,k),Tmp(i,j,k))
             enddo
          enddo
       enddo
    endif

  end subroutine add_tgv

  ! modifs with Camille for STBL RANS
!!$  !===============================================================================
!!$  subroutine add_lam_BL
!!$  !===============================================================================
!!$    !> Initialization of laminar boundary layer
!!$    !> only for (bottom wall) (Cartesian)
!!$    !> NEEDS Re_inlet OR jdel
!!$  !===============================================================================
!!$    use mod_blasius ! <- module Blasius solution [Initial_condition]
!!$    use mod_interp1 ! <- module for 1D interpolation [Mathematics]
!!$    use mod_comm    ! <- communications SEND/RECV [Parallel]
!!$    use mod_deriv   ! <- module for 1D derivative [Derivatives]
!!$    implicit none
!!$    ! ---------------------------------------------------------------------------
!!$    integer :: i,j,k
!!$    ! ---------------------------------------------------------------------------    
!!$    ! Definition of initial boundary layer
!!$    real(wp) :: delta
!!$    real(wp) :: x_0,Re_x0,xi,x_start,Re_outlet
!!$    real(wp), dimension(:),     allocatable :: u_b,rho_b,T_b,Psi_b
!!$    real(wp), dimension(:,:,:), allocatable :: psib,dpsib
!!$    ! ---------------------------------------------------------------------------
!!$
!!$    ! Definition of initial boundary layer
!!$    ! ====================================
!!$    ! if is_similarity=.false. : compressible similarity solution
!!$    ! if is_similarity= .true. : polynomial approximation
!!$
!!$    if (is_similarity) then
!!$       ! Compressible similarity solution
!!$       ! ================================
!!$       ! allocate var
!!$       allocate(u_bl(neta),v_bl(neta),rho_bl(neta),T_bl(neta),mu_bl(neta),G_bl(neta),M_bl(neta))
!!$       allocate(cp_bl(neta),lambda_bl(neta),Pr_bl(neta),feta(neta),ffeta(neta))
!!$       allocate(eta(neta),etai(neta))
!!$
!!$       ! compressible similarity solution
!!$       call compute_blasius
!!$
!!$       ! var. in global grid for interpolation
!!$       allocate(u_b(ngy),rho_b(ngy),T_b(ngy),Psi_b(ngy))
!!$       ! var. in local grid
!!$       allocate(psib(nx1:nx2,ny1:ny2,nz1:nz2),dpsib(nx1:nx2,ny1:ny2,nz1:nz2))
!!$       psib = 0.0_wp
!!$       dpsib= 0.0_wp
!!$
!!$       ! x_start = Re_inlet**2*mu_ref/(rho_ref*u_ref)/matto**2
!!$       ! based on Re_L
!!$       ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!!$       x_start = Re_inlet**2*mu_ref/(rho_ref*u_ref)
!!$       deltas_in = Re_inlet*mu_ref/(rho_ref*u_ref)
!!$
!!$       x_start=0.0_wp
!!$       uu=u_ref
!!$       do k=1,nz
!!$          do i=50,nx
!!$             x_0=xc(i,1)+x_start
!!$             Re_x0=rho_ref*u_ref*x_0/mu_ref
!!$             xi=rho_ref*u_ref*mu_ref*x_0
!!$
!!$             ! interpolate similarity solution on global grid
!!$             call interp1(ygc(i,1:ngy),  u_b,ngy,etai*x_0/sqrt(Re_x0),  u_bl,neta,'linear')
!!$             call interp1(ygc(i,1:ngy),rho_b,ngy,etai*x_0/sqrt(Re_x0),rho_bl,neta,'linear')
!!$             call interp1(ygc(i,1:ngy),  T_b,ngy,etai*x_0/sqrt(Re_x0),  T_bl,neta,'linear')
!!$             call interp1(ygc(i,1:ngy),Psi_b,ngy,etai*x_0/sqrt(Re_x0), feta ,neta,'linear')
!!$
!!$             ! fill local var. and dimensionalize
!!$             do j=ndyt,nfyt
!!$                uu (i,j,k) =  u_b(j+coord(2)*ny)*u_ref
!!$                rho(i,j,k) =rho_b(j+coord(2)*ny)*rho_ref
!!$                tmp(i,j,k) =  T_b(j+coord(2)*ny)*T_ref
!!$                psib(i,j,k)=Psi_b(j+coord(2)*ny)*sqrt(2.0_wp*xi)
!!$             enddo
!!$          enddo
!!$       enddo
!!$
!!$       Re_outlet = matto*sqrt(u_ref*rho_ref*(x_start+xg(ngx))/mu_ref)
!!$       if (iproc.eq.0) then
!!$          write(*,*) '==========================================='
!!$          write(*,*) 'Domain:'
!!$          write(*,*) 'deltas_in      = ', deltas_in
!!$          write(*,*) 'x_start        = ', x_start
!!$          write(*,*) 'xmax/deltas_ref= ', xg(ngx)/deltas_in
!!$          write(*,*) 'ymax/deltas_ref= ', yg(ngy)/deltas_in
!!$          write(*,*) 'zmax/deltas_ref= ', (zg(ngz)-zg(1))/deltas_in
!!$          write(*,*) 'Red*_in        = ', Re_inlet
!!$          write(*,*) 'Red*_out       = ', Re_outlet
!!$          write(*,*) '==========================================='
!!$       endif
!!$
!!$       ! compute vertical velocity
!!$       if (is_2d) then
!!$          call communic2d(psib)
!!$       else
!!$          call communic3d(psib)
!!$       endif
!!$       call deriv_x_11pts(psib,dpsib)
!!$       vv=-dpsib/rho
!!$
!!$       ! deallocate var
!!$       deallocate(u_bl,v_bl,rho_bl,T_bl,mu_bl,G_bl,M_bl,cp_bl,lambda_bl,Pr_bl,eta,etai,feta)
!!$       deallocate(psib,dpsib)
!!$       deallocate(u_b,rho_b,T_b,Psi_b)
!!$
!!$    else
!!$       ! Polynomial approx of Blasius solution (Att. jdel to be defined) in param.ini ?????
!!$       ! =====================================
!!$       !delta = yg(jdel)
!!$       delta = ygc(1,jdel)
!!$       ! var. in global grid
!!$       allocate(u_b(ngy))     
!!$       ! polynomial approx
!!$       do j=1,jdel
!!$          !u_b(j)=((yg(j)/delta)*(2.0_wp-2.0_wp*(yg(j)/delta)**2+(yg(j)/delta)**3))
!!$          u_b(j)=((ygc(1,j)/delta)*(2.0_wp-2.0_wp*(ygc(1,j)/delta)**2+(ygc(1,j)/delta)**3))
!!$       enddo
!!$       do j=jdel+1,ngy
!!$          u_b(j)=1.0_wp
!!$       enddo
!!$       ! fill local var. and dimensionalize
!!$       do k=1,nz
!!$          do i=1,nx
!!$             do j=ndyt,nfyt
!!$                uu(i,j,k)=u_b(j+coord(2)*ny)*u_ref
!!$             enddo
!!$          enddo
!!$       enddo
!!$       ! deallocate var
!!$       deallocate(u_b)
!!$       
!!$    end if
!!$    
!!$  end subroutine add_lam_BL
  
  !===============================================================================
  subroutine add_BL
  !===============================================================================
    !> Initialization of boundary layers
    !> only for jmin & jmax (cartesian)
  !===============================================================================
    use mod_RFM
    use mod_database
    use mod_interp1  ! <- module for 1D interpolation [Mathematics]
    use mod_blasius  ! <- module Blasius solution [Initial_condition]
    use mod_comm    ! <- communications SEND/RECV [Parallel]
    use mod_deriv   ! <- module for 1D derivative [Derivatives]
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k,l,jbc
    real(wp) :: x_0,Re_x0,xi,x_start,Re_outlet
    integer, dimension(2) :: ndeb_j,nend_j
    real(wp), dimension(2) :: yw_BC,yw_BC_proc
    real(wp), dimension(:), allocatable :: um1,y_int1
    real(wp), dimension(:),     allocatable :: u_b,rho_b,T_b,Psi_b
    real(wp), dimension(:,:,:), allocatable :: psib,dpsib
    ! ---------------------------------------------------------------------------

    if ((is_bc_wall(2,1)).and.(is_bc_wall(2,2))) then
       yw_BC(1)=y(1); yw_BC(2)=y(ny)
       ! Indices where to begin & end
       ndeb_j(1)=1; nend_j(1)=ny/2
       ndeb_j(2)=ny/2+1; nend_j(2)=ny
    else
       ! Indices where to begin & end
       ndeb_j=1; nend_j=ny
       ! if 2 walls, which one to consider to init the flow for the processor
       if ((bl(nob(iproc))%BC(3).eq.0).and.(bl(nob(iproc))%BC(4).eq.0)) then
          if ((coord(2)+1)*ny.le.ngy/2) then
             ndeb_j(2)=0; nend_j(2)=0
          else
             ndeb_j(1)=0; nend_j(1)=0
          endif
       else if (bl(nob(iproc))%BC(3).eq.0) then
          ndeb_j(2)=0; nend_j(2)=0
       else if (bl(nob(iproc))%BC(3).eq.0) then
          ndeb_j(1)=0; nend_j(1)=0
       else
          ndeb_j=0; nend_j=0
       endif
       ! Scatter wall coordinates (from coord(2)==0/ndomj-1 process to all processes)
       ! ------------------------
       if (is_curv3) then
          call mpistop('Not implemented in mod_RFM.f90',1)
       else if (is_curv) then
          call mpistop('Not implemented in mod_RFM.f90',1)
       else
          yw_BC=0.0_wp; yw_BC_proc=0.0_wp
          ! Location at jmin & jmax (if wall)
          if (is_bc_wall(2,1)) yw_BC_proc(1)=y(1)
          if (is_bc_wall(2,2)) yw_BC_proc(2)=y(ny)
          call MPI_ALLREDUCE(yw_BC_proc,yw_BC,2,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXZ,info)
       endif
    endif

    ! ---------------
    ! Define BL field
    ! ---------------
    loopjbc: do jbc=1,2
       ! Test if necessary to add BL
       ! ---------------------------
       if ((bl(nob(iproc))%BC(2+jbc).eq.0).or.&  ! Wall for this block BC
           (ndeb_j(jbc).ne.0)) then              ! This wall is "closer" than the opposite one (in term of indices)

          ! Determination of interpolation profile
          ! --------------------------------------
          ! Profil put in ascendant order
          allocate(y_int1(1:ny))
          if (jbc.eq.1) then
             if (y(2)-y(1).gt.0.0_wp) then
                do j=1,ny
                   y_int1(j) = y(j)-yw_BC(jbc)
                enddo
             else
                do j=1,ny
                   ! l=nend_j(jbc)+1-j
                   y_int1(j) = yw_BC(jbc)-y(j)
                enddo
             endif
          else if (jbc.eq.2) then
             if (y(2)-y(1).lt.0.0_wp) then
                do j=1,ny
                   l=ny+1-j
                   y_int1(j) = y(l)-yw_BC(jbc)
                enddo
             else
                do j=1,ny
                   l=ny+1-j
                   y_int1(j) = yw_BC(jbc)-y(l)
                enddo
             endif
          endif

          ! Computation of laminar boundary layer
          ! -------------------------------------
          if ((is_LBL).and.(.not.is_RFM_TBL(jbc))) then
             
             ! Compressible similarity solution
             ! --------------------------------
             if (is_similarity) then
                ! allocate var
                allocate(u_bl(neta),v_bl(neta),rho_bl(neta),T_bl(neta),mu_bl(neta),G_bl(neta),M_bl(neta))
                allocate(cp_bl(neta),lambda_bl(neta),Pr_bl(neta),feta(neta),ffeta(neta))
                allocate(eta(neta),etai(neta))

                ! compressible similarity solution
                call compute_blasius

                ! var. in global grid for interpolation
                allocate(u_b(ny),rho_b(ny),T_b(ny),Psi_b(ny))
                ! var. in local grid
                allocate(psib(nx1:nx2,ny1:ny2,nz1:nz2),dpsib(nx1:nx2,ny1:ny2,nz1:nz2))
                psib = 0.0_wp
                dpsib= 0.0_wp

                ! x_start = Re_inlet**2*mu_ref/(rho_ref*u_ref)/matto**2
                ! based on Re_L
                x_start = Re_inlet**2*mu_ref/(rho_ref*u_ref)
                deltas_in = Re_inlet*mu_ref/(rho_ref*u_ref)

                do k=1,nz
                   do i=1,nx
                      x_0=x(i)+x_start
                      Re_x0=rho_ref*u_ref*x_0/mu_ref
                      xi=rho_ref*u_ref*mu_ref*x_0

                      ! interpolate similarity solution on global grid
                      call interp1(y_int1(1:ny),  u_b,ny,etai*x_0/sqrt(Re_x0),  u_bl,neta,'linear')
                      call interp1(y_int1(1:ny),rho_b,ny,etai*x_0/sqrt(Re_x0),rho_bl,neta,'linear')
                      call interp1(y_int1(1:ny),  T_b,ny,etai*x_0/sqrt(Re_x0),  T_bl,neta,'linear')
                      call interp1(y_int1(1:ny),Psi_b,ny,etai*x_0/sqrt(Re_x0), feta ,neta,'linear')

                      ! fill local var. and dimensionalize
                      if (jbc.eq.1) then
                         do j=ndeb_j(jbc),nend_j(jbc)
                            uu (i,j,k) =  u_b(j)*u_ref
                            rho(i,j,k) =rho_b(j)*rho_ref
                            tmp(i,j,k) =  T_b(j)*T_ref
                            psib(i,j,k)=Psi_b(j)*sqrt(2.0_wp*xi)
                         enddo
                      else
                         do j=1,nend_j(jbc)+1-ndeb_j(jbc)
                            l=ny+1-j
                            uu (i,l,k) =  u_b(j)*u_ref
                            rho(i,l,k) =rho_b(j)*rho_ref
                            tmp(i,l,k) =  T_b(j)*T_ref
                            psib(i,l,k)=Psi_b(j)*sqrt(2.0_wp*xi)
                         enddo
                      endif
                   enddo
                enddo

                ! Re_outlet = matto*sqrt(u_ref*rho_ref*(x_start+xg(ngx))/mu_ref)
                ! if (iproc.eq.0) then
                !    write(*,*) '==========================================='
                !    write(*,*) 'Domain:'
                !    write(*,*) 'deltas_in      = ', deltas_in
                !    write(*,*) 'x_start        = ', x_start
                !    write(*,*) 'xmax/deltas_ref= ', xg(ngx)/deltas_in
                !    write(*,*) 'ymax/deltas_ref= ', yg(ngy)/deltas_in
                !    write(*,*) 'zmax/deltas_ref= ', (zg(ngz)-zg(1))/deltas_in
                !    write(*,*) 'Red*_in        = ', Re_inlet
                !    write(*,*) 'Red*_out       = ', Re_outlet
                !    write(*,*) '==========================================='
                ! endif

                ! compute vertical velocity
                if (is_2d) then
                   call communic2d(psib)
                else
                   call communic3d(psib)
                endif
                call deriv_x_11pts(psib,dpsib)
                vv(:,ndeb_j(jbc):nend_j(jbc),:)=-dpsib(:,ndeb_j(jbc):nend_j(jbc),:)/rho(:,ndeb_j(jbc):nend_j(jbc),:)

                ! deallocate var
                deallocate(u_bl,v_bl,rho_bl,T_bl,mu_bl,G_bl,M_bl,cp_bl,lambda_bl,Pr_bl,eta,etai,feta,ffeta)
                deallocate(psib,dpsib)
                deallocate(u_b,rho_b,T_b,Psi_b)
                deallocate(y_int1)
             else
                ! Polynomial approx of Blasius solution
                ! =====================================
                ! (Att. jdel to be defined) in param.ini ?????
                ! Maybe defined delta differently ????
                ! delta = yg(jdel)
                call mpistop('Polynomial approx. to be implemented in a general manner',1)
             end if

          ! Computation of turbulent boundary layer
          ! ---------------------------------------
          elseif (is_RFM_TBL(jbc)) then
             ! Read Reynolds stresses
             call read_database(base_TBL,Re_prof_TBL(jbc))
             deallocate(uu_db,vv_db,ww_db,uv_db,uw_db,vw_db)
             ! Define mean field
             allocate(um1(ny))
             ! Interp
             call interp1(y_int1/d99_TBL(jbc),um1,ny,y_db,um_db,ny_db,'linear')
             if (jbc.eq.1) then
                do j=ndeb_j(jbc),nend_j(jbc)
                      uu(:,j,:) = u_ref*um1(j)
                enddo
             else
                do j=1,nend_j(jbc)+1-ndeb_j(jbc)
                   l=ny+1-j
                      uu(:,l,:) = u_ref*um1(j)
                enddo
             endif
             deallocate(um1,y_int1,um_db,y_db)
            
          ! Damp velocity near the wall
          ! ---------------------------
          else
             if (is_bc_wall(2,jbc)) then
                if (jbc.eq.1) then
                   do j=1,10
                      uu(:,j,:)=uu(:,j,:)/10.0_wp*dble(j-1)
                   enddo
                else
                   do j=1,10
                      l=ny+1-j
                      uu(:,l,:)=uu(:,l,:)/10.0_wp*dble(j-1)
                   enddo
                endif
             endif
          endif

       endif

    enddo loopjbc

  end subroutine add_BL

  !===============================================================================
  subroutine add_lam_BL
  !===============================================================================
    !> Initialization of laminar boundary layer
    !> only for (bottom wall) (Cartesian)
    !> NEEDS Re_inlet OR jdel
  !===============================================================================
    use mod_blasius ! <- module Blasius solution [Initial_condition]
    use mod_interp1 ! <- module for 1D interpolation [Mathematics]
    use mod_comm    ! <- communications SEND/RECV [Parallel]
    use mod_deriv   ! <- module for 1D derivative [Derivatives]
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    ! ---------------------------------------------------------------------------    
    ! Definition of initial boundary layer
    real(wp) :: delta
    real(wp) :: x_0,Re_x0,xi,x_start,Re_outlet
    real(wp), dimension(:),     allocatable :: u_b,rho_b,T_b,Psi_b
    real(wp), dimension(:,:,:), allocatable :: psib,dpsib
    ! ---------------------------------------------------------------------------

    ! Definition of initial boundary layer
    ! ====================================
    ! if is_similarity=.false. : compressible similarity solution
    ! if is_similarity= .true. : polynomial approximation

    if (is_similarity) then
       ! Compressible similarity solution
       ! ================================
       ! allocate var
       allocate(u_bl(neta),v_bl(neta),rho_bl(neta),T_bl(neta),mu_bl(neta),G_bl(neta),M_bl(neta))
       allocate(cp_bl(neta),lambda_bl(neta),Pr_bl(neta),feta(neta),ffeta(neta))
       allocate(eta(neta),etai(neta))

       ! compressible similarity solution
       call compute_blasius

       ! var. in global grid for interpolation
       allocate(u_b(ngy),rho_b(ngy),T_b(ngy),Psi_b(ngy))
       ! var. in local grid
       allocate(psib(nx1:nx2,ny1:ny2,nz1:nz2),dpsib(nx1:nx2,ny1:ny2,nz1:nz2))
       psib = 0.0_wp
       dpsib= 0.0_wp

       ! x_start = Re_inlet**2*mu_ref/(rho_ref*u_ref)/matto**2
       ! based on Re_L
       ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       x_start = Re_inlet**2*mu_ref/(rho_ref*u_ref)
       deltas_in = Re_inlet*mu_ref/(rho_ref*u_ref)

       do k=1,nz
          do i=1,nx
             x_0=x(i)+x_start
             Re_x0=rho_ref*u_ref*x_0/mu_ref
             xi=rho_ref*u_ref*mu_ref*x_0

             ! interpolate similarity solution on global grid
             call interp1(yg(1:ngy),  u_b,ngy,etai*x_0/sqrt(Re_x0),  u_bl,neta,'linear')
             call interp1(yg(1:ngy),rho_b,ngy,etai*x_0/sqrt(Re_x0),rho_bl,neta,'linear')
             call interp1(yg(1:ngy),  T_b,ngy,etai*x_0/sqrt(Re_x0),  T_bl,neta,'linear')
             call interp1(yg(1:ngy),Psi_b,ngy,etai*x_0/sqrt(Re_x0), feta ,neta,'linear')

             ! fill local var. and dimensionalize
             do j=ndyt,nfyt
                uu (i,j,k) =  u_b(j+coord(2)*ny)*u_ref
                rho(i,j,k) =rho_b(j+coord(2)*ny)*rho_ref
                tmp(i,j,k) =  T_b(j+coord(2)*ny)*T_ref
                psib(i,j,k)=Psi_b(j+coord(2)*ny)*sqrt(2.0_wp*xi)
             enddo
          enddo
       enddo

       Re_outlet = matto*sqrt(u_ref*rho_ref*(x_start+xg(ngx))/mu_ref)
       if (iproc.eq.0) then
          write(*,*) '==========================================='
          write(*,*) 'Domain:'
          write(*,*) 'deltas_in      = ', deltas_in
          write(*,*) 'x_start        = ', x_start
          write(*,*) 'xmax/deltas_ref= ', xg(ngx)/deltas_in
          write(*,*) 'ymax/deltas_ref= ', yg(ngy)/deltas_in
          write(*,*) 'zmax/deltas_ref= ', (zg(ngz)-zg(1))/deltas_in
          write(*,*) 'Red*_in        = ', Re_inlet
          write(*,*) 'Red*_out       = ', Re_outlet
          write(*,*) '==========================================='
       endif

       ! compute vertical velocity
       if (is_2d) then
          call communic2d(psib)
       else
          call communic3d(psib)
       endif
       call deriv_x_11pts(psib,dpsib)
       vv=-dpsib/rho

       ! deallocate var
       deallocate(u_bl,v_bl,rho_bl,T_bl,mu_bl,G_bl,M_bl,cp_bl,lambda_bl,Pr_bl,eta,etai,feta)
       deallocate(psib,dpsib)
       deallocate(u_b,rho_b,T_b,Psi_b)

    else
       ! Polynomial approx of Blasius solution (Att. jdel to be defined) in param.ini ?????
       ! =====================================
       delta = yg(jdel)
       ! var. in global grid
       allocate(u_b(ngy))     
       ! polynomial approx
       do j=1,jdel
          u_b(j)=((yg(j)/delta)*(2.0_wp-2.0_wp*(yg(j)/delta)**2+(yg(j)/delta)**3))
       enddo
       do j=jdel+1,ngy
          u_b(j)=1.0_wp
       enddo
       ! fill local var. and dimensionalize
       do k=1,nz
          do i=1,nx
             do j=ndyt,nfyt
                uu(i,j,k)=u_b(j+coord(2)*ny)*u_ref
             enddo
          enddo
       enddo
       ! deallocate var
       deallocate(u_b)
       
    end if
    
  end subroutine add_lam_BL

  !===============================================================================
  subroutine add_turb_BL
  !===============================================================================
    !> Initialization of turbulent boundary layer
    !> only for (bottom wall) (Cartesian)
    !> NEEDS turbulent BL database
  !===============================================================================
    use mod_RFM
    use mod_database
    use mod_interp1   ! <- module for 1D interpolation [Mathematics]
    implicit none
    ! ---------------------------------------------------------------------------
    ! integer :: i,j,k
    ! real(wp), dimension(:), allocatable :: um1
    ! ---------------------------------------------------------------------------

    call mpistop('add_turb_BL needs to be updated',1)

    ! ! Read RFM mode parameters
    ! ! ------------------------
    ! call read_param_RFM

    ! ! Read Reynolds stresses
    ! ! ----------------------
    ! call read_database(base,Re_prof)
    ! deallocate(uu_db,vv_db,ww_db,uv_db,uw_db,vw_db)

    ! ! Define mean field (for mean flow convection or field initialization)
    ! ! -----------------
    ! if (is_field) then
    !    allocate(um1(ngy))
    !    ! L_ref corresponds to d99 = u_inf^+*Re_tau of Schlatter database
    !    call interp1(yg(1:ngy)/L_ref,um1,ngy,y_db,um_db,ny_db,'linear')
    !    do k=1,nz
    !       do j=1,ny
    !          do i=1,nx
    !             uu(i,j,k) = u_ref*um1(j+coord(2)*ny)
    !          enddo
    !       enddo
    !    enddo
    !    deallocate(um1,um_db,y_db)
    ! else
    !    uu = 0.0_wp
    !    deallocate(um_db,y_db)
    ! endif

  end subroutine add_turb_BL

  !===============================================================================
  subroutine add_BL_over_cavity
  !===============================================================================
    !> Initialization of laminar boundary layer spanning a cavity
    !> NEEDS jdel
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    ! ---------------------------------------------------------------------------    
    ! Definition of boundary layer (polynomial approx)
    real(wp) :: delta
    real(wp), dimension(:), allocatable :: u_b
    ! ---------------------------------------------------------------------------

    if (is_curv) then
       if (nob(iproc)<4) then
          delta = yg(jdel)
          ! var. in global grid
          allocate(u_b(ngy))     
          ! polynomial approx
          do j=1,jdel
             u_b(j)=((yg(j)/delta)*(2.0_wp-2.0_wp*(yg(j)/delta)**2+(yg(j)/delta)**3))
          enddo
          do j=jdel+1,ngy
             u_b(j)=1.0_wp
          enddo
          ! fill local var. and dimensionalize
          do k=1,nz
             do i=1,nx
                do j=1,ny
                   uu(i,j,k)=u_b(j+coord(2)*ny)*u_ref
                enddo
             enddo
          enddo
          ! free memory
          deallocate(u_b)
       else
          uu=0.0_wp
       endif
    else
       if (nob(iproc)<4) then
          delta = yg(jdel)
          ! var. in global grid
          allocate(u_b(ngy))     
          ! polynomial approx
          do j=1,jdel
             u_b(j)=((yg(j)/delta)*(2.0_wp-2.0_wp*(yg(j)/delta)**2+(yg(j)/delta)**3))
          enddo
          do j=jdel+1,ngy
             u_b(j)=1.0_wp
          enddo
          ! fill local var. and dimensionalize
          do k=1,nz
             do i=1,nx
                do j=1,ny
                   uu(i,j,k)=u_b(j+coord(2)*ny)*u_ref
                enddo
             enddo
          enddo
          ! free memory
          deallocate(u_b)
       else
          uu=0.0_wp
       endif
    endif

  end subroutine add_BL_over_cavity

  !===============================================================================
  subroutine init_vel_chan
  !===============================================================================
    !> Initialization of velocity field in plane channel flow
  !===============================================================================
    use mod_RFM
    use mod_database
    use mod_interp1   ! <- module for 1D interpolation [Mathematics]
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: yy,eps,epsilon1,pis2
    real(wp), dimension(:), allocatable :: um1
    ! ----------------------------------------------------------------------------

    ! ~> ref velocity is u_ref=U_bulk
    ! ~> ref length is channel half-width: L_ref=hc

    if (is_RFM) then
       ! Read Reynolds stresses
       ! ----------------------
       call read_database(base_FST,Re_prof_FST)
       deallocate(uu_db,vv_db,ww_db,uv_db,uw_db,vw_db)

       ! Define mean field (for mean flow convection or field initialization)
       ! -----------------
       if (is_field) then
          allocate(um1(ngy))
          call interp1(yg(1:ngy)/hc,um1,ngy,y_db,um_db,ny_db,'linear')
          do k=1,nz
             do j=1,ny
                do i=1,nx
                   uu(i,j,k) = u_ref*um1(j+coord(2)*ny)
                enddo
             enddo
          enddo
          deallocate(um1,um_db,y_db)
       else
          uu = 0.0_wp
          deallocate(um_db,y_db)
       endif
    else
       ! Poiseuille profile (between -hc and hc)
       ! ==================
       allocate(um1(ngy))
       if (is_2D) then
          pis2=0.5*acos(-1.)
          do j=1,ngy
             if (is_curv) then
                yy=ygc(10,j)/L_ref
             else
                yy=yg(j)/L_ref
             endif
             um1(j)=(cos(pis2*yy))**2/1.5
          enddo
       else
          do j=1,ngy
             yy=yg(j)/L_ref
             um1(j)=1.0_wp-yy**2
          enddo
       endif

       ! Initialization of velocity components
       ! =====================================
       epsilon1=0.2_wp
       !epsilon1=0.0_wp
       if (is_2D) epsilon1=0.0_wp
       !call random_seed()
   !    do k=1,nz
   !       do j=1,ny
   !          do i=1,nx
   !             call random_number(eps)
   !             eps = epsilon1*2.0_wp*(eps-0.5_wp)
   !             uu(i,j,k)= um1(j+coord(2)*ny)*(1.0_wp+eps)*1.5_wp*u_ref
   !             vv(i,j,k)= uu(i,j,k)*eps
   !             ww(i,j,k)= uu(i,j,k)*eps
   !          enddo
   !       enddo
   !    enddo

       ! Power law for turbulent velocity profile: n=7
       ! =============================================
       do k=1,nz
          do j=1,ny
             do i=1,nx

                if (yc(i,j).lt.0.0_wp) then
                   uu(i,j,k)= 8.0_wp/7.0_wp*u_ref*(1.0_wp+yc(i,j)/L_ref)**(1.0_wp/7.0_wp)
                else
                   uu(i,j,k)= 8.0_wp/7.0_wp*u_ref*(1.0_wp-yc(i,j)/L_ref)**(1.0_wp/7.0_wp)
                endif

                call random_number(eps)
                eps = epsilon1*2.0_wp*(eps-0.5_wp)
                uu(i,j,k)= um1(j+coord(2)*ny)*(1.0_wp+eps)*1.5_wp*u_ref
                vv(i,j,k)= uu(i,j,k)*eps
                ww(i,j,k)= uu(i,j,k)*eps
             enddo
          enddo
       enddo

    endif

  end subroutine init_vel_chan

  !===============================================================================
  subroutine init_vel_phill
  !===============================================================================
    !> Initialization periodic hill flow
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer  :: i,j,k
    real(wp) :: yy,eps,epsilon1,pis2,Ly
    real(wp), dimension(:), allocatable :: um1
    ! ----------------------------------------------------------------------------

    ! ~> ref velocity is u_ref=U_bulk
    ! ~> ref length is hill height: L_ref=hc
    ! ~> we need the height of the inlet section 2.035*L_ref
    if (iproc==0) Ly=ygc(1,ngy)-ygc(1,1)
    call MPI_BCAST(Ly,1,MPI_DOUBLE_PRECISION,0,COMM_global,info)
    if (iproc==0) print *,'in init vel_phill: Ly=',Ly,2.035*L_ref

    ! Poiseuille profile
    ! ==================
    allocate(um1(ngy))
    if (is_2D) then
       pis2=0.5*acos(-1.)
       do j=1,ngy
          !yy=yg(j)/L_ref
          !yy=ygc(1,j)/L_ref
          yy=2.0_wp*(ygc(1,j)-L_ref)/Ly-1.0_wp
          !um1(j)=(cos(pis2*yy))**2
          !um1(j)=1.0_wp-yy**2
          um1(j)=1.0_wp-yy**2
       enddo
    else
       do j=1,ngy
          yy=yg(j)/L_ref
          um1(j)=1.0_wp-yy**2
       enddo
    endif

!!$    if (iproc==0) then
!!$       print *,'L_ref',L_ref
!!$       print *,yg/L_ref
!!$       print *,yg(2)-yg(1),yg(3)-yg(2),yg(4)-yg(3)
!!$       print *,''
!!$       print *,um1
!!$    endif
!!$    stop
!!$    Ubp=0.0_wp
!!$    if (coord(1)==0) then
!!$       do j=1,ny
!!$          !Ub=Ub+uu(1,j,1)*dyi(j)
!!$          Ubp=Ubp+um1(j+coord(2)*ny)*dyi(j)
!!$       enddo
!!$    endif
!!$    call MPI_ALLREDUCE(Ubp,Ub,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
!!$    print *,L_ref,Ly
!!$    Ub=Ub/Ly
!!$    !Ub=Ub*deltaz
!!$    !print *,'Qi=',Qi,Qi/Ly/Lz/rho_infty!/c_infty

    ! Initialization of velocity components
    ! =====================================
    !epsilon1 = 0.2_wp
    epsilon1 = 0.0_wp
    !call random_seed()
    
    ! ~> ref velocity is U_bulk
    ! ~> ref length is hill height L_ref=hc
   
    
    if (is_2D) then

!!$       do k=1,nz
!!$          do j=1,ny
!!$             do i=1,nx
!!$                !call random_number(eps)
!!$                eps = epsilon1*2.0_wp*(eps-0.5_wp)
!!$                uu(i,j,k) = um1(j+coord(2)*ny)*(1.0_wp+eps)*2.0_wp*u_ref
!!$                vv(i,j,k) = uu(i,j,k)*eps
!!$             enddo
!!$          enddo
!!$       enddo

       do k=1,nz
          do i=1,nx
             do j=1,ny
                yy=2.0_wp*(ygc(i+coord(1)*nx,j+coord(2)*ny)-L_ref)/Ly-1.0_wp
                ! yy=2.0_wp*(ygc(i,j)-L_ref)/Ly-1.0_wp
                uu(i,j,k)=(1.0_wp-yy**2)*1.5_wp*u_ref
                !              if (ygc(i+coord(1)*nx,j+coord(2)*ny)>L_ref) then
                !                 uu(i,j,k)=(cos(pis2*yy))**2*2.0_wp*u_ref
                !              else
                !                 uu(i,j,k)=0.0_wp
                !              endif
             enddo
          enddo
       enddo

!!$       call read_init2D_phill
!!$       Q_ip=0.
!!$       if (coord(1)==0) then
!!$          do j=1,ny
!!$             Q_ip=Q_ip+rho_ref*uu(1,j,1)*dyi(j)
!!$          enddo
!!$       endif
!!$       call MPI_ALLREDUCE(Q_ip,Q_i,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
!!$       uu(:,:,1)=uu(:,:,1)*u_ref*Ly*Lz*rho_ref/Q_i

    else

       do k=1,nz
          do j=1,ny
             do i=1,nx
                call random_number(eps)
                eps = epsilon1*2.0_wp*(eps-0.5_wp)
                uu(i,j,k) = um1(j+coord(2)*ny)*(1.0_wp+eps)*1.5_wp*u_ref
                vv(i,j,k) = uu(i,j,k)*eps
                ww(i,j,k) = uu(i,j,k)*eps
             enddo
          enddo
       enddo

    endif

  end subroutine init_vel_phill

  !===============================================================================
  subroutine read_init2D_phill
  !===============================================================================
    !> read 'init.bin' for initialization of periodic hill from low-Re 2D sol.
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,nnn
    real(wp) :: t
    real(wp), dimension(ngx,ngy) :: varg2d
    ! ----------------------------------------------------------------------------

    open(47,file='init.bin',form='unformatted',status='unknown')
    rewind(47)
    read(47) nnn
    read(47) t

    ! Read rho for each proc
    ! ----------------------
    read(47) ((varg2d(i,j),i=1,ngx),j=1,ngy)
    call MPI_BARRIER(COMM_global,info)
    do j=1,ny
       do i=1,nx
          rho(i,j,1)=varg2d(i+coord(1)*nx,j+coord(2)*ny)
       enddo
    enddo
    call MPI_BARRIER(COMM_global,info)

    ! Read u for each proc
    ! --------------------
    read(47) ((varg2d(i,j),i=1,ngx),j=1,ngy)
    call MPI_BARRIER(COMM_global,info)
    do j=1,ny
       do i=1,nx
          uu(i,j,1)=varg2d(i+coord(1)*nx,j+coord(2)*ny)
       enddo
    enddo
    call MPI_BARRIER(COMM_global,info)

    ! Read v for each proc
    ! --------------------
    read(47) ((varg2d(i,j),i=1,ngx),j=1,ngy)
    call MPI_BARRIER(COMM_global,info)
    do j=1,ny
       do i=1,nx
          vv(i,j,1)=varg2d(i+coord(1)*nx,j+coord(2)*ny)
       enddo
    enddo
    call MPI_BARRIER(COMM_global,info)

    ! Read p for each proc
    ! --------------------
    read(47) ((varg2d(i,j),i=1,ngx),j=1,ngy)
    call MPI_BARRIER(COMM_global,info)
    do j=1,ny
       do i=1,nx
          prs(i,j,1)=varg2d(i+coord(1)*nx,j+coord(2)*ny)
       enddo
    enddo
    call MPI_BARRIER(COMM_global,info)

    close(47)

  end subroutine read_init2D_phill

  !===============================================================================
  subroutine init_vel_cyl
  !===============================================================================
    !> Initialization of velocity field in cylinder flow
  !===============================================================================
    use mod_block
    use mod_bc_periodicity ! <~ needed for Lzp  TO BE CHANGED
    use warnstop
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: radius,aa,rr
    complex(wp), parameter :: ii=(0.0_wp,1.0_wp)
    complex(wp) :: zg_c
    complex(wp), dimension(ngx,ngy) :: wg_c
    ! ----------------------------------------------------------------------------
    real(wp) :: r,theta,eps
    ! ----------------------------------------------------------------------------

    if (is_curv3) then
    !if ((is_curv3).or.(is_curv)) then

       ! Sphere radius =diameter(Lref)/2
       ! -------------------------------
       radius=L_ref/2.0_wp

       ! Amplitude
       ! ---------
       aa=-u_ref*radius**3/2.0_wp

       ! Compute velocity components
       ! ---------------------------
       do k=1,nz
          do j=1,ny
             do i=1,nx
                rr=sqrt(xc3(i,j,k)**2+yc3(i,j,k)**2+zc3(i,j,k)**2)
                !rr=sqrt(xc(i,j)**2+yc(i,j)**2+z(k)**2)

!!$                uu(i,j,k)=u_ref
!!$                vv(i,j,k)=0.0_wp
!!$                ww(i,j,k)=0.0_wp
!!$                uu(i,j,k)=0.0_wp
!!$                vv(i,j,k)=0.0_wp
!!$                ww(i,j,k)=u_ref
                !uu(i,j,k)=aa*(3.0_wp*xc3(i,j,k)*xc3(i,j,k)/rr**3-1.0_wp/rr**5)
                uu(i,j,k)=aa*(3.0_wp*xc3(i,j,k)*xc3(i,j,k)/rr**5-1.0_wp/rr**3)+u_ref
                vv(i,j,k)=aa* 3.0_wp*xc3(i,j,k)*yc3(i,j,k)/rr**5
                ww(i,j,k)=aa* 3.0_wp*xc3(i,j,k)*zc3(i,j,k)/rr**5
             enddo
          enddo
       enddo

       ! Damp velocity near the cylinder wall
       ! ------------------------------------
       ! (to avoid violation of no-slip condition)
!!$       if (is_bc_wall(2,1)) then
!!$          do k=1,nz
!!$             do j=1,10
!!$                do i=1,nx
!!$                   uu(i,j,k)=uu(i,j,k)/10.0_wp*dble(j-1)
!!$                   vv(i,j,k)=vv(i,j,k)/10.0_wp*dble(j-1)
!!$                   ww(i,j,k)=ww(i,j,k)/10.0_wp*dble(j-1)
!!$                enddo
!!$             enddo
!!$          enddo
!!$       endif
       if (is_bc_wall(3,1)) then
          do k=1,10
             do j=1,ny
                do i=1,nx
                   uu(i,j,k)=uu(i,j,k)/10.0_wp*dble(k-1)
                   vv(i,j,k)=vv(i,j,k)/10.0_wp*dble(k-1)
                   ww(i,j,k)=ww(i,j,k)/10.0_wp*dble(k-1)
                enddo
             enddo
          enddo
       endif

       return
    endif
    
    if (iorder_visc==0) then

       ! Cylinder radius =diameter(Lref)/2 
       ! =================================
       radius=L_ref/2.0_wp

       ! Compute complex potential velocity field
       ! ========================================
       do j=1,ngy
          do i=1,ngx
             zg_c=xgc(i,j)+ii*ygc(i,j)
             wg_c(i,j)=u_ref*(1.0_wp-(radius)**2/zg_c**2)
          enddo
       enddo

       ! Compute velocity components
       ! ===========================
       do k=1,nz
          do j=1,ny
             do i=1,nx
                uu(i,j,k)=  real(wg_c(i+coord(1)*nx,j+coord(2)*ny))
                vv(i,j,k)=-aimag(wg_c(i+coord(1)*nx,j+coord(2)*ny))
             enddo
          enddo
       enddo

       return
    endif

    if ((nbloc==12).or.(nbloc==6)) then

       ! Lzp=zg(ngz)-zg(1)+deltaz
       Lzp=zmax-zmin+deltaz
       !print *,'Lzp',Lzp/L_ref
       
       ! Damp velocity near the cylinder wall using the AS1 smoothing function
       ! =====================================================================
       ! (to avoid violation of no-slip condition)
       uu=u_ref
       vv=0.0_wp
       ww=0.0_wp
       eps=0.1_wp

       do k=1,nz
          do j=1,ny
             do i=1,nx
                r=sqrt(xc(i,j)**2+yc(i,j)**2)/L_ref

                if (is_bc_wall(2,1).and.j==1) then
                   uu(i,1,k)=0.0_wp
                elseif (r<1.5_wp.and.r>0.5_wp) then
                   uu(i,j,k)=u_ref*AS1(r,1.0_wp,0.5_wp)

                   if (r>0.5_wp.and.r<=1.0_wp) then
                      ww(i,j,k)=eps*u_ref*AS1(r,0.75_wp,0.25_wp)*sin(pi*z(k)/Lzp)**2
                   else
                      ww(i,j,k)=eps*u_ref*(1.0_wp-AS1(r,1.25_wp,0.25_wp))*sin(pi*z(k)/Lzp)**2
                   endif

                   if (yc(i,j)>=0.0_wp) then
                      !if ((xc(i,j)/(r*L_ref)>1.0_wp).or.(xc(i,j)/(r*L_ref)<-1.0_wp)) print *,i,j,k,r,xc(i,j),xc(i,j)/(r*L_ref)
                      aa=xc(i,j)/(r*L_ref)
                      if (aa> 1.0_wp) aa=1.0_wp
                      if (aa<-1.0_wp) aa=-1.0_wp
                      theta=acos(aa)
                      if ((theta.ne.pi).and.(theta.ne.0.0_wp)) then
                         if (r>0.5_wp.and.r<=1.0_wp) then
                            vv(i,j,k)=eps*u_ref*AS1(r,0.75_wp,0.25_wp) &
                                 *AS1(theta,pi/2.0_wp,pi/2.0_wp)
                         else
                            vv(i,j,k)=eps*u_ref*(1.0_wp-AS1(r,1.25_wp,0.25_wp)) &
                                 *AS1(theta,pi/2.0_wp,pi/2.0_wp)
                         endif
                      endif
                   else
                      !theta=2.0_wp*pi-acos(xc(i,j)/(r*L_ref))
                      theta=pi+acos(xc(i,j)/(r*L_ref))
                      if ((theta.ne.pi).and.(theta.ne.2.0_wp*pi)) then
                         if (r>0.5_wp.and.r<=1.0_wp) then
                            vv(i,j,k)=eps*u_ref*AS1(r,0.75_wp,0.25_wp) &
                                 *(1.0_wp-AS1(theta,3.0_wp*pi/2.0_wp,pi/2.0_wp))
                         else
                            vv(i,j,k)=eps*u_ref*(1-AS1(r,1.25_wp,0.25_wp)) &
                                 *(1.0_wp-AS1(theta,3.0_wp*pi/2.0_wp,pi/2.0_wp))
                         endif
                      endif
                   endif
                endif
             enddo
          enddo
       enddo

    else

       ! Cylinder radius =diameter(Lref)/2 
       ! =================================
       radius=L_ref/2.0_wp

       ! Compute complex potential velocity field
       ! ========================================
       ! do j=1,ngy
       !    do i=1,ngx
       !       zg_c=xgc(i,j)+ii*ygc(i,j)
       !       wg_c(i,j)=u_ref*(1.0_wp-(radius)**2/zg_c**2)
       !    enddo
       ! enddo
       ! Using local grid
       do j=1,ny
          do i=1,nx
             zg_c=xc(i,j)+ii*yc(i,j)
             wg_c(i,j)=u_ref*(1.0_wp-(radius)**2/zg_c**2)
          enddo
       enddo

       ! Compute velocity components
       ! ===========================
       do k=1,nz
          do j=1,ny
             do i=1,nx
                ! uu(i,j,k)=  real(wg_c(i+coord(1)*nx,j+coord(2)*ny))
                ! vv(i,j,k)=-aimag(wg_c(i+coord(1)*nx,j+coord(2)*ny))
                uu(i,j,k)=  real(wg_c(i,j))
                vv(i,j,k)=-aimag(wg_c(i,j))
             enddo
          enddo
       enddo

       ! Damp velocity near the cylinder wall
       ! ====================================
       ! (to avoid violation of no-slip condition)
       !if (coord(2)==0) then
       if (is_bc_wall(2,1)) then
          do k=1,nz
             do j=1,10
                do i=1,nx
                   uu(i,j,k)=uu(i,j,k)/10.0_wp*dble(j-1)
                   vv(i,j,k)=vv(i,j,k)/10.0_wp*dble(j-1)
                enddo
             enddo
          enddo
       endif

    endif
    
  end subroutine init_vel_cyl

  !==============================================================================
  function AS1(r,rm,ra)
  !==============================================================================
    !> Wall function for potential flow solution around a circular cylinder (AS1)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp) :: r,rm,ra,AS1
    ! ---------------------------------------------------------------------------

    AS1 = 0.5_wp*(1.0_wp+tanh((2.0_wp*ra*(r-rm)/(ra**2-(r-rm)**2))))
    
  end function AS1

  !===============================================================================
  subroutine init_turbine
  !===============================================================================
    !> Initialization of turbine flow
  !===============================================================================
    use mod_eos  ! for rg (gas cst)
    use mod_block
    use warnstop
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: rr,ks,uax,utheta
    complex(wp), parameter :: ii=(0.0_wp,1.0_wp)
    real(wp) :: angl
    ! real(wp), dimension(:), allocatable :: xmin,xmax
    real(wp) :: M_in,M_out
    real(wp), dimension(:,:,:), allocatable :: M_loc
    ! ----------------------------------------------------------------------------
    ! real(wp) :: r,theta
    ! ----------------------------------------------------------------------------
    ! real(wp) :: Lx
    real(wp) :: fac,den1,den2,Ei1,Ei2,R1,R2,alpha_1

    ! ~> ref velocity is u_ref=incoming freestram
    ! ~> ref length is blade chord: L_ref

    !if (iorder_visc==0) then

       ! Annular duct (test of radial equilibrium)
       ! =========================================
       if ((nbloc==1.or.nbloc==2).and.(is_curv3)) then

          prs= p_ref
          Tmp= T_ref

          !dimensions
          R1=0.2_wp
          R2=0.28_wp

          ! free vortex flow
          ! ================
          ! strength
          ks=20.11_wp
          ! axial velocity: 5 m/s
          uax=5.0_wp*4.

          fac=ks**2/(2.0_wp*rg*T_ref)
          ! compute exponential integrals
          call e1xb(fac/R1**2,Ei1)
          call e1xb(fac/R2**2,Ei2)
          den1=-fac*Ei1+R1**2*exp(-fac/R1**2)
          den2=-fac*Ei2+R2**2*exp(-fac/R2**2)
          alpha_1=p_ref*(R2**2-R1**2)/(den2-den1)

          do k=1,nz
             do j=1,ny
                do i=1,nx
                   rr=sqrt(yc3(i,j,k)**2+zc3(i,j,k)**2)
                   !rr=sqrt(xc3(i,j,k)**2+yc3(i,j,k)**2+zc3(i,j,k)**2)
                   utheta=ks/rr
                   uu(i,j,k)= uax
                   vv(i,j,k)=-utheta*zc3(i,j,k)/rr
                   ww(i,j,k)= utheta*yc3(i,j,k)/rr
                   prs(i,j,k)=alpha_1*exp(-fac/rr**2)
                enddo
             enddo
          enddo

          rho= prs/(rg*Tmp)

       ! Ailette
       ! =======
       elseif (nbloc==6) then

          ! Initialization of thermodynamic variables to reference values
          ! =============================================================
          prs= 0.92*p_ref
          !prs= p_ref
          Tmp= T_ref
          rho= prs/(rg*Tmp)

          angl=10.0_wp
          angl=36.0_wp

          ! Compute velocity components
          ! ===========================
          do k=1,nz
             do j=1,ny
                do i=1,nx
                   uu(i,j,k)=u_ref*cos(angl/180.0_wp*pi)
                   vv(i,j,k)=u_ref*sin(angl/180.0_wp*pi)
                enddo
             enddo
          enddo

          ! Damp velocity near the wall
          ! ===========================
          ! (to avoid violation of non-penetrability condition)
          if (is_bc_wall(2,1)) then
             do k=1,nz
                do j=1,10
                   do i=1,nx
                      uu(i,j,k)=uu(i,j,k)/10.0_wp*dble(j-1)
                      vv(i,j,k)=vv(i,j,k)/10.0_wp*dble(j-1)
                   enddo
                enddo
             enddo
          endif
          if (is_bc_wall(2,2)) then
             do k=1,nz
                do j=ny-9,ny
                   do i=1,nx
                      uu(i,j,k)=uu(i,j,k)/10.0_wp*dble(ny-j)
                      vv(i,j,k)=vv(i,j,k)/10.0_wp*dble(ny-j)
                   enddo
                enddo
             enddo
          endif

       ! test of inlet
       ! =============
       elseif (nbloc==4) then

          ! Compute velocity components
          ! ===========================
          do k=1,nz
             do j=1,ny
                do i=1,nx
                   uu(i,j,k)=u_ref!*0.0_wp
                   vv(i,j,k)=0.0_wp
                enddo
             enddo
          enddo

!!$          ! Damp velocity near the  wall
!!$          ! ============================
!!$          ! (to avoid violation of non-penetrability condition)
!!$          if (is_bc_wall(2,1)) then
!!$             do k=1,nz
!!$                do j=1,10
!!$                   do i=1,nx
!!$                      uu(i,j,k)=uu(i,j,k)/10.0_wp*dble(j-1)
!!$                      vv(i,j,k)=vv(i,j,k)/10.0_wp*dble(j-1)
!!$                   enddo
!!$                enddo
!!$             enddo
!!$          endif
!!$          if (is_bc_wall(2,2)) then
!!$             do k=1,nz
!!$                do j=ny-9,ny
!!$                   do i=1,nx
!!$                      uu(i,j,k)=uu(i,j,k)/10.0_wp*dble(ny-j)
!!$                      vv(i,j,k)=vv(i,j,k)/10.0_wp*dble(ny-j)
!!$                   enddo
!!$                enddo
!!$             enddo
!!$          endif

       ! LS59 rotor blade cascade / Baumgartner
       ! ========================
       elseif ((nbloc==9).or.(nbloc==12).or.(nbloc==13)) then
       ! TEST elseif (nbloc==2) then

          ! Initialization of thermodynamic variables to reference values
          ! =============================================================
          ! LS59
          prs= 0.534460302143549_wp*p_ref
          !prs= 0.423153711961457_wp*p_ref
          !prs= 0.655370258341356_wp*p_ref
          prs= 0.85_wp*p_ref
          prs= 0.80_wp*p_ref
          ! LS59 HSU
          !prs= 0.53_wp*p_ref
          ! Baumgartner
          !prs= p_ref/4.35_wp
          !prs= p_ref/2.0_wp
          !prs= p_ref
          Tmp= T_ref
          rho= prs/(rg*Tmp)

          ! Compute velocity components
          ! ===========================
          ! LS59
          angl=30.0_wp
          ! Baumgartner
          !angl=0.0_wp

!!$          do k=1,nz
!!$             do j=1,ny
!!$                do i=1,nx
!!$                   uu(i,j,k)=u_ref*cos(angl/180.0_wp*pi)
!!$                   vv(i,j,k)=u_ref*sin(angl/180.0_wp*pi)
!!$                enddo
!!$             enddo
!!$          enddo

          ! Impose linear increase of Mach number through vane
          ! ==================================================
          allocate(M_loc(nx,ny,nz))
          M_in=Mach
          !M_out=1.2-M_in ! Baumgartner
          M_out=0.98-M_in ! LS-59
          !M_out=0.6-M_in ! LS-59 Stefan (Sanz et.al -> subsonic)

          if (is_curv3) then
             do k=1,nz
                do j=1,ny
                   do i=1,nx
                      if (xc3(i,j,k).le.-0.5_wp*L_ref) then
                         M_loc(i,j,k)=M_in
                      elseif (xc3(i,j,k).gt.0.5_wp*L_ref) then
                         M_loc(i,j,k)=M_out+M_in
                      else
                         M_loc(i,j,k)=M_out*((xc3(i,j,k)+0.5_wp*L_ref)/L_ref + M_in/M_out)
                      endif
                   enddo
                enddo
             enddo
          else
             do k=1,nz
                do j=1,ny
                   do i=1,nx
                      if (xc(i,j).le.0._wp) then
                         M_loc(i,j,k)=M_in
                      elseif (xc(i,j).gt.L_ref) then
                         M_loc(i,j,k)=M_out+M_in
                      else
                         M_loc(i,j,k)=M_out*(xc(i,j)/L_ref + M_in/M_out)
                      endif
                   enddo
                enddo
             enddo
          endif

          ! Compute thermos and velocity
          !angl=-65.0_wp
          angl=0.0_wp
          do k=1,nz
             do j=1,ny
                do i=1,nx
                   ! /!\ p_ref and T_ref are stagnation quantities /!\
                   prs(i,j,k)= p_ref/(1.0_wp+(gam-1.0_wp)/2.0_wp*M_loc(i,j,k)**2)**(gam/(gam-1.0_wp))
                   Tmp(i,j,k)= T_ref/(1.0_wp+(gam-1.0_wp)/2.0_wp*M_loc(i,j,k)**2)
                   rho(i,j,k)= prs(i,j,k)/(rg*Tmp(i,j,k))
                   uu(i,j,k)=M_loc(i,j,k)*sqrt(c2calc_tro(Tmp(i,j,k),rho(i,j,k)))*cos(angl*pi/180.0_wp)
                   vv(i,j,k)=M_loc(i,j,k)*sqrt(c2calc_tro(Tmp(i,j,k),rho(i,j,k)))*sin(angl*pi/180.0_wp)
                enddo
             enddo
          enddo
          deallocate(M_loc)

          ! Damp velocity near the wall
          ! ===========================
          ! (to avoid violation of non-penetrability condition)
          if (is_bc_wall(1,2)) then
             do k=1,nz
                do j=1,ny
                   do i=nx-9,nx
                      uu(i,j,k)=uu(i,j,k)*dble(nx-i)/10.0_wp
                      vv(i,j,k)=vv(i,j,k)*dble(nx-i)/10.0_wp
                   enddo
                enddo
             enddo
          endif
          if (is_bc_wall(2,1)) then
             do k=1,nz
                do j=1,10
                   do i=1,nx
                      uu(i,j,k)=uu(i,j,k)*dble(j-1)/10.0_wp
                      vv(i,j,k)=vv(i,j,k)*dble(j-1)/10.0_wp
                   enddo
                enddo
             enddo
          endif

       ! LS89 vane passage
       ! =================
       elseif (nbloc==14) then

          ! Initialization of thermodynamic variables to reference values
          ! =============================================================
          prs= p_ref/5.12_wp
          Tmp= T_ref
          rho= prs/(rg*Tmp)

          ! Compute velocity components
          ! ===========================
          angl=0.0_wp

          do k=1,nz
             do j=1,ny
                do i=1,nx
                   uu(i,j,k)=u_ref*cos(angl/180.0_wp*pi)
                   vv(i,j,k)=u_ref*sin(angl/180.0_wp*pi)
                enddo
             enddo
          enddo

          ! Damp velocity near the wall
          ! ===========================
          ! (to avoid violation of non-penetrability condition)
          if (is_bc_wall(1,2)) then
             do k=1,nz
                do j=1,ny
                   do i=nx-9,nx
                      uu(i,j,k)=uu(i,j,k)/10.0_wp*dble(nx-i)
                      vv(i,j,k)=vv(i,j,k)/10.0_wp*dble(nx-i)
                   enddo
                enddo
             enddo
          endif
          if (is_bc_wall(2,1)) then
             do k=1,nz
                do j=1,10
                   do i=1,nx
                      uu(i,j,k)=uu(i,j,k)/10.0_wp*dble(j-1)
                      vv(i,j,k)=vv(i,j,k)/10.0_wp*dble(j-1)
                   enddo
                enddo
             enddo
          endif

       ! Idealized blade vane configuration
       ! ==================================
       elseif (nbloc==32) then

          ! Initialization of thermodynamic variables to reference values
          ! =============================================================
          ! prs= 0.2283*p_ref
          ! prs= 0.5*p_ref
          Tmp= T_ref
          ! rho= prs/(rg*Tmp)

          prs = p_ref
          rho= prs/(rg*Tmp)

          ! Lx = 3.55*L_ref
          ! do k=1,nz
          !    do j=1,ny
          !       do i=1,nx
          !          prs(i,j,k) = p_ref - 0.75*p_ref*(xc(i,j) + 1.0_wp*L_ref)/Lx
          !          rho(i,j,k) = prs(i,j,k)/(rg*Tmp(i,j,k))
          !          if ((iproc.eq.10).and.(i.eq.45).and.(j.eq.25)) print *,'ICI proc 10',prs(i,j,k)
          !          if ((iproc.eq.17).and.(i.eq.20).and.(j.eq.25)) print *,'ICI proc 17',prs(i,j,k)
          !       enddo
          !    enddo
          ! enddo

          ! Compute velocity components
          ! ===========================
          angl=0.0_wp
          vv = 0.0_wp
          M_out = 1.6_wp

          uu = u_ref
          ! do k=1,nz
          !    do j=1,ny
          !       do i=1,nx
          !          uu(i,j,k) = (M_out-Mach)*c_ref*(xc(i,j) + 1.0_wp*L_ref)/Lx + c_ref*Mach
          !       enddo
          !    enddo
          ! enddo

          ! Damp velocity near the wall
          ! ===========================
          ! (to avoid violation of non-penetrability condition)
          if (is_bc_wall(1,1)) then
             do k=1,nz
                do j=1,ny
                   do i=1,5
                      uu(i,j,k)=uu(i,j,k)/30.0_wp*dble(i-1)
                   enddo
                enddo
             enddo
          endif
          if (is_bc_wall(2,1)) then
             do k=1,nz
                do j=1,5
                   do i=1,nx
                      uu(i,j,k)=uu(i,j,k)/30.0_wp*dble(j-1)
                   enddo
                enddo
             enddo
          endif

       endif

    !endif

  end subroutine init_turbine

  !===============================================================================
  subroutine init_src
  !===============================================================================
    !> Initialization source/pulse/vortex
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    !integer :: i,j,k,ip
    ! ---------------------------------------------------------------------------

!!$    ! determination des proc qui contiennent les coupes ! TO BE CHANGED needed for what ?????????????????????
!!$    ! =================================================
!!$    allocate(ixy(0:nproc-1))
!!$    ixy = .false.
!!$
!!$    if (coord(3).eq.0) ixy(iproc) = .true.
!!$
!!$    do ip=0,nproc-1
!!$       call MPI_BCAST(ixy(ip), 1, MPI_LOGICAL, ip, COMM_global, info)
!!$    enddo
!!$
!!$    call MPI_BARRIER(COMM_global,info)

    if (iproc.eq.0) write (*,*) 'Initialisation OK'

  end subroutine init_src

  !===============================================================================
  subroutine add_shock
  !===============================================================================
    !> Add shock
  !===============================================================================
    use mod_eos
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: xc ! shock location 
    real(wp) :: p1,rho1,u1 ! left state
    real(wp) :: p2,rho2,u2 ! right state
    ! ---------------------------------------------------------------------------
    ! real(wp) :: M1_2,g

    ! Rankine Hugoniot 
    ! ================

    ! State 1 (left)
    ! --------------
    u1=u_ref
    rho1=rho_ref
    p1=p_ref

    ! State 2 (right)
    ! --------------
    p2=p1*(1.0_wp+2.0_wp*gam*(Mach**2-1.0_wp)/(gam+1.0_wp))
    rho2=rho1*(gam+1.0_wp)*Mach**2/(2+(gam-1.0_wp)*Mach**2)   
    u2=u1*rho1/rho2

    ! shock location
    ! --------------
    xc=1.0_wp
    
    ! initialize field
    ! ----------------
    do i=1,nx
       do j=1,ny
          do k=1,nz

             if (x(i)<xc) then  ! left state
                rho(i,j,k)=rho1
                uu(i,j,k)=u1
                vv(i,j,k)=0.0_wp
                prs(i,j,k)=p1
             else               ! right state                          
                rho(i,j,k)=rho2
                uu(i,j,k)=u2
                vv(i,j,k)=0.0_wp
                prs(i,j,k)=p2
             endif

          enddo
       enddo
    enddo

  end subroutine add_shock

  !===============================================================================
  subroutine add_pulse
  !===============================================================================
    !> Add pressure pulse
  !===============================================================================
    use mod_eos
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: ar
    ! ----------------------------------------------------------------------------

    ! Add pressure pulse
    ! ==================
    
    ! Gaussian width
    ar=-log(2.0_wp)/(b_pulse*deltax)**2

    if (is_2D) then
       if (is_curv) then
          do j=1,ny
             do i=1,nx
                prs(i,j,1)=p_ref+ampl_pulse*exp(ar*((xc(i,j)-x_pulse)**2+(yc(i,j)-y_pulse)**2))
             enddo
          enddo
       else
          do j=1,ny
             do i=1,nx
                prs(i,j,1)=p_ref+ampl_pulse*exp(ar*((x(i)-x_pulse)**2+(y(j)-y_pulse)**2))
             enddo
          enddo
       endif
    else
       if (is_curv3) then
          do k=1,nz
             do j=1,ny
                do i=1,nx
                   prs(i,j,k)=p_ref+ampl_pulse &
                             *exp(ar*((xc3(i,j,k)-x_pulse)**2+(yc3(i,j,k)-y_pulse)**2+(zc3(i,j,k)-z_pulse)**2))
                enddo
             enddo
          enddo
       else
          if (is_curv) then
             do k=1,nz
                do j=1,ny
                   do i=1,nx
                      prs(i,j,k)=p_ref+ampl_pulse*exp(ar*((xc(i,j)-x_pulse)**2+(yc(i,j)-y_pulse)**2+(z(k)-z_pulse)**2))
                   enddo
                enddo
             enddo
          else
             do k=1,nz
                do j=1,ny
                   do i=1,nx
                      prs(i,j,k)=p_ref+ampl_pulse*exp(ar*((x(i)-x_pulse)**2+(y(j)-y_pulse)**2+(z(k)-z_pulse)**2))
                   enddo
                enddo
             enddo
          endif
       endif
    endif

    ! we need Tmp and rho instead of prs
    do k=1,nz
       do j=1,ny
          do i=1,nx
!!$             Tmp(i,j,k)= prs(i,j,k)/p_ref*T_ref
!!$             rho(i,j,k)= rocalc_pt(prs(i,j,k),Tmp(i,j,k),Tmp(i,j,k))
             !Tmp(i,j,k)= T_ref
             !rho(i,j,k)= rocalc_pt(prs(i,j,k),Tmp(i,j,k),rho_ref)
             !rho(i,j,k)=prs(i,j,k)/c_ref**2
             !rho(i,j,k)=rho_ref**2
             Tmp(i,j,k)=prs(i,j,k)/(rg*rho(i,j,k))
          enddo
       enddo
    enddo
          
    !if (iproc.eq.0) print *,'Add pulse OK'
    
  end subroutine add_pulse

  !==============================================================================
  subroutine add_vortex
  !==============================================================================
    !> Add vortex
    !==============================================================================
    use mod_vortex_model ! <- vortex models
    use mod_eos          ! <- rocalc_pt
    use warnstop
    use mod_utils ! TEMPORAIRE test comm1d
    use mod_interface ! TEMPORAIRE test comm1d
    use mod_comm1 ! TEMPORAIRE test comm1d
    use mod_grid_directions ! TEMPORAIRE test comm1d
    use mod_mpi_types_one_sided
    use mod_comm ! TEMPORAIRE test comm1d
    use mod_mpi_types_two_sided
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    real(wp), dimension(nx,ny) :: u_vortex,v_vortex,p_vortex
    real(wp), dimension(nx,ny,nz) :: u3_vortex,v3_vortex,p3_vortex
    ! ---------------------------------------------------------------------------
!!$    integer :: nghex,dim

    if (is_curv3) then
       ! Define vortex
       ! =============
       call vortex_3d_model(type_vortex,x_vortex,y_vortex,u3_vortex,v3_vortex,p3_vortex)

       ! Add vortex
       ! ==========
       do k=1,nz
          do j=1,ny
             do i=1,nx
                uu(i,j,k)= uu(i,j,k) + u3_vortex(i,j,k)
                vv(i,j,k)= vv(i,j,k) + v3_vortex(i,j,k)
                prs(i,j,k) =  p_ref  + p3_vortex(i,j,k)
             enddo
          enddo
       enddo
    else
       ! Define vortex
       ! =============
       call vortex_2d_model(type_vortex,x_vortex,y_vortex,u_vortex,v_vortex,p_vortex)

       ! Add vortex
       ! ==========
       do k=1,nz
          do j=1,ny
             do i=1,nx
                uu(i,j,k)= uu(i,j,k) + u_vortex(i,j)
                vv(i,j,k)= vv(i,j,k) + v_vortex(i,j)
                prs(i,j,k) =  p_ref  + p_vortex(i,j)
             enddo
          enddo
       enddo
    endif

!!$    ! TEST ONE-SIDED COMMS
!!$    rho_n=rho
!!$    rhou_n=uu
!!$    rhov_n=vv
!!$    rhoe_n=prs
!!$    call communication
!!$    
!!$    open(194,file='init1_p'//trim(numchar(iproc))//'_ex.bin',form='unformatted',status='unknown')
!!$    rewind(194)
!!$    write(194) nx+2*ngh
!!$    write(194) ny+2*ngh
!!$    if (is_curv) then
!!$       write(194) ((xc(i,j),i=1-ngh,nx+ngh),j=1-ngh,ny+ngh)
!!$       write(194) ((yc(i,j),i=1-ngh,nx+ngh),j=1-ngh,ny+ngh)
!!$    else
!!$       write(194) (x(i),i=1-ngh,nx+ngh)
!!$       write(194) (y(j),j=1-ngh,ny+ngh)
!!$    endif
!!$    write(194) ((rho_n(i,j,1),i=1-ngh,nx+ngh),j=1-ngh,ny+ngh)
!!$    write(194) ((rhou_n(i,j,1),i=1-ngh,nx+ngh),j=1-ngh,ny+ngh)
!!$    write(194) ((rhov_n(i,j,1),i=1-ngh,nx+ngh),j=1-ngh,ny+ngh)
!!$    write(194) ((rhoe_n(i,j,1),i=1-ngh,nx+ngh),j=1-ngh,ny+ngh)
!!$    close(194)
!!$
!!$    call mpistop('init vortex+comm',0)

!!$    ! TEST ONE-SIDED COMMS
!!$    ! Add vortex
!!$    ! ==========
!!$    do k=1,nz
!!$       do j=1,ny
!!$          do i=1,nx
!!$             uu(i,j,k)= uu(i,j,k) + u_vortex(i,j)*(2.0_wp*dble(k-1)/dble(nz-1)-1)
!!$             vv(i,j,k)= vv(i,j,k) + v_vortex(i,j)*(2.0_wp*dble(k-1)/dble(nz-1)-1)
!!$             prs(i,j,k) =  p_ref  + p_vortex(i,j)*(2.0_wp*dble(k-1)/dble(nz-1)-1)
!!$          enddo
!!$       enddo
!!$    enddo
!!$
!!$    j=ny/2
!!$    !if (iproc==0) print *,vv(10,j,:)
!!$    rho_n=rho
!!$    rhou_n=uu
!!$    rhov_n=vv
!!$    rhoe_n=prs
!!$    call communication
!!$
!!$    !do i=1,4
!!$    !   if (neighbor(i)>=0) print *,'moi proc',iproc,' echange',i,' avec',neighbor(i),' dans sa direction',ndir(i)
!!$    !enddo
!!$    !z(-4:0)=0.
!!$    !z(nx+1:nx+5)=0.
!!$    k=49
!!$    print *,nz,k
!!$    
!!$    open(194,file='init1_p'//trim(numchar(iproc))//'_ex.bin',form='unformatted',status='unknown')
!!$    rewind(194)
!!$    write(194) nx+2*ngh
!!$    write(194) ny+2*ngh
!!$    write(194) nz+2*ngh
!!$    if (is_curv) then
!!$       write(194) ((xc(i,j),i=1-ngh,nx+ngh),j=1-ngh,ny+ngh)
!!$       write(194) ((yc(i,j),i=1-ngh,nx+ngh),j=1-ngh,ny+ngh)
!!$    else
!!$       write(194) (x(i),i=1-ngh,nx+ngh)
!!$       write(194) (y(j),j=1-ngh,ny+ngh)
!!$    endif
!!$    write(194) (z(k),k=1-ngh,nz+ngh)
!!$    !write(194) (z(k),k=1,nz)
!!$    ! check YZ
!!$    i=nx/2
!!$    write(194) ((rho_n(i,j,k),j=1-ngh,ny+ngh),k=1-ngh,nz+ngh)
!!$    write(194) ((rhou_n(i,j,k),j=1-ngh,ny+ngh),k=1-ngh,nz+ngh)
!!$    write(194) ((rhov_n(i,j,k),j=1-ngh,ny+ngh),k=1-ngh,nz+ngh)
!!$    write(194) ((rhoe_n(i,j,k),j=1-ngh,ny+ngh),k=1-ngh,nz+ngh)
!!$!!    write(194) ((rho_n(i,j,k),j=1-ngh,ny+ngh),k=1,nz)
!!$!!    write(194) ((rhou_n(i,j,k),j=1-ngh,ny+ngh),k=1,nz)
!!$!!    write(194) ((rhov_n(i,j,k),j=1-ngh,ny+ngh),k=1,nz)
!!$!!    write(194) ((rhoe_n(i,j,k),j=1-ngh,ny+ngh),k=1,nz)
!!$!!    write(194) ((rho_n(i,j,k),i=1-ngh,nx+ngh),k=1,nz)
!!$!!    write(194) ((rhou_n(i,j,k),i=1-ngh,nx+ngh),k=1,nz)
!!$!!    write(194) ((rhov_n(i,j,k),i=1-ngh,nx+ngh),k=1,nz)
!!$!!    write(194) ((rhoe_n(i,j,k),i=1-ngh,nx+ngh),k=1,nz)
!!$!!    write(194) ((rho_n(i,j,k),i=1-ngh,nx+ngh),k=1-ngh,nz+ngh)
!!$!!    write(194) ((rhou_n(i,j,k),i=1-ngh,nx+ngh),k=1-ngh,nz+ngh)
!!$!!    write(194) ((rhov_n(i,j,k),i=1-ngh,nx+ngh),k=1-ngh,nz+ngh)
!!$!!    write(194) ((rhoe_n(i,j,k),i=1-ngh,nx+ngh),k=1-ngh,nz+ngh)
!!$!!    write(194) ((rho_n(i,j,k),i=1-ngh,nx+ngh),j=1-ngh,ny+ngh)
!!$!!    write(194) ((rhou_n(i,j,k),i=1-ngh,nx+ngh),j=1-ngh,ny+ngh)
!!$!!    write(194) ((rhov_n(i,j,k),i=1-ngh,nx+ngh),j=1-ngh,ny+ngh)
!!$!!    write(194) ((rhoe_n(i,j,k),i=1-ngh,nx+ngh),j=1-ngh,ny+ngh)
!!$    close(194)
!!$
!!$    call mpistop('init vortex+comm',0)

!!$    open(194,file='init1_p'//trim(numchar(iproc))//'_ex.bin',form='unformatted',status='unknown')
!!$    rewind(194)
!!$    write(194) nx+2*ngh
!!$    write(194) ny+2*ngh
!!$    if (is_curv) then
!!$       write(194) ((xc(i,j),i=1-ngh,nx+ngh),j=1-ngh,ny+ngh)
!!$       write(194) ((yc(i,j),i=1-ngh,nx+ngh),j=1-ngh,ny+ngh)
!!$    else
!!$       write(194) (x(i),i=1-ngh,nx+ngh)
!!$       write(194) (y(j),j=1-ngh,ny+ngh)
!!$    endif
!!$    write(194) ((rho_n(i,j,k),i=1-ngh,nx+ngh),j=1-ngh,ny+ngh)
!!$    write(194) ((rhou_n(i,j,k),i=1-ngh,nx+ngh),j=1-ngh,ny+ngh)
!!$    write(194) ((rhov_n(i,j,k),i=1-ngh,nx+ngh),j=1-ngh,ny+ngh)
!!$    write(194) ((rhoe_n(i,j,k),i=1-ngh,nx+ngh),j=1-ngh,ny+ngh)
!!$    close(194)
!!$
!!$    call mpistop('init vortex+comm',0)


!!$    if (iproc==0) print *,'TEST ONE-SIDED COMMS for increments'
!!$    ! TEST ONE-SIDED COMMS for increments
!!$    ngh_irs=5
!!$    nghex=8
!!$
!!$    if (BC_face(1,1)%sort>0) ngh_irs(1)=nghex
!!$    if (BC_face(1,2)%sort>0) ngh_irs(2)=nghex
!!$    if (BC_face(2,1)%sort>0) ngh_irs(3)=nghex
!!$    if (BC_face(2,2)%sort>0) ngh_irs(4)=nghex
!!$
!!$    !tv_bl in grid.f90
!!$!!    if (iproc==0) then
!!$!!       ! DEFAULT ngh_irs=[2,7,8,5,2,2]
!!$!!       ngh_irs=[2,7,8,5,2,2]
!!$!!    endif
!!$!!    if (iproc==1) then
!!$!!       ! DEFAULT ngh_irs=[2,7,12,8,2,2]
!!$!!       ngh_irs=[8,12,7,2,2,2]
!!$!!    endif
!!$    
!!$    !tr_bl in grid.f90
!!$!!    if (iproc==0) then
!!$!!       ! DEFAULT ngh_irs=[5,8,2,7,2,2]
!!$!!       ngh_irs=[7,2,8,5,2,2]
!!$!!    endif
!!$!!    if (iproc==1) then
!!$!!       ! DEFAULT ngh_irs=[8,5,2,7,2,2]
!!$!!       ngh_irs=[8,5,2,7,2,2]
!!$!!    endif
!!$
!!$!!    i=1; j=1
!!$!!    print *,'(1,1)',2*(i-1)+j
!!$!!    i=1; j=2
!!$!!    print *,'(1,2)',2*(i-1)+j
!!$!!    i=2; j=1
!!$!!    print *,'(2,1)',2*(i-1)+j
!!$!!    i=2; j=2
!!$!!    print *,'(2,2)',2*(i-1)+j
!!$!!    i=3; j=1
!!$!!    print *,'(3,1)',2*(i-1)+j
!!$!!    i=3; j=2
!!$!!    print *,'(3,2)',2*(i-1)+j
!!$!!    call mpistop('essai',0)
!!$
!!$    ! initialization of ngh_irs in neighboring domain
!!$    do i=1,3
!!$       do j=1,2
!!$          BC_face(i,j)%ngh_irs=[2,2,2,2,2,2]
!!$       enddo
!!$    enddo
!!$
!!$    if (is_2D) then
!!$       dim=2
!!$    else
!!$       dim=3
!!$    endif
!!$    
!!$    do i=1,dim
!!$       do j=1,2
!!$          k=2*(i-1)+j
!!$          if (neighbor(k)>=0) then
!!$             call MPI_SENDRECV(ngh_irs,6,MPI_INTEGER,neighbor(k),tag, &
!!$                 BC_face(i,j)%ngh_irs,6,MPI_INTEGER,neighbor(k),tag,COMM_global,status,info)  
!!$          endif
!!$       enddo
!!$    enddo
!!$
!!$    !!call mpistop('essai',0)
!!$    !!print *,ngh_irs(1:4),iproc
!!$    
!!$    nx1_irs= 1-ngh_irs(1)
!!$    nx2_irs=nx+ngh_irs(2)
!!$    ny1_irs= 1-ngh_irs(3)
!!$    ny2_irs=ny+ngh_irs(4)
!!$    nz1_irs= 1-ngh_irs(5)
!!$    nz2_irs=nz+ngh_irs(6)
!!$    if (is_2d) then
!!$       nz1_irs=1
!!$       nz2_irs=1
!!$    endif
!!$    ! Reallocate increments arrays with extended ghost points
!!$    !!deallocate(Krho,Krhou,Krhov,Krhow,Krhoe)
!!$    allocate( Krho(nx1_irs:nx2_irs,ny1_irs:ny2_irs,nz1_irs:nz2_irs))
!!$    allocate(Krhou(nx1_irs:nx2_irs,ny1_irs:ny2_irs,nz1_irs:nz2_irs))
!!$    allocate(Krhov(nx1_irs:nx2_irs,ny1_irs:ny2_irs,nz1_irs:nz2_irs))
!!$    allocate(Krhow(nx1_irs:nx2_irs,ny1_irs:ny2_irs,nz1_irs:nz2_irs))
!!$    allocate(Krhoe(nx1_irs:nx2_irs,ny1_irs:ny2_irs,nz1_irs:nz2_irs))
!!$    ! Definition of MPI types for communications
!!$    !call mpi_types_comm_inc
!!$    call mpi_types_comm1_inc
!!$   
!!$    ! Define windows on allocated arrays (for one-sided comms)
!!$    call mpi_win_comm1_inc       
!!$     Krho=rho_ref
!!$    Krhou=0.0_wp
!!$    Krhov=-u_ref
!!$    Krhoe=p_ref
!!$     Krho(1:nx,1:ny,1:nz)=rho(1:nx,1:ny,1:nz)
!!$    Krhou(1:nx,1:ny,1:nz)= uu(1:nx,1:ny,1:nz)
!!$    Krhov(1:nx,1:ny,1:nz)= vv(1:nx,1:ny,1:nz)
!!$    Krhoe(1:nx,1:ny,1:nz)=prs(1:nx,1:ny,1:nz)
!!$    !call communication_inc(Krho,Krhou,Krhov,Krhow,Krhoe)
!!$    call communication1_inc          
!!$    open(194,file='init1_p'//trim(numchar(iproc))//'_ex.bin',form='unformatted',status='unknown')
!!$    rewind(194)
!!$    write(194) ngh_irs
!!$    write(194) is_swapij2
!!$    write(194) is_rev2(:,1)
!!$    write(194) is_rev2(:,2)
!!$    write(194) nx+ngh_irs(1)+ngh_irs(2)
!!$    write(194) ny+ngh_irs(3)+ngh_irs(4)
!!$    write(194) ((Krho(i,j,1),i=nx1_irs,nx2_irs),j=ny1_irs,ny2_irs)
!!$    write(194) ((Krhou(i,j,1),i=nx1_irs,nx2_irs),j=ny1_irs,ny2_irs)
!!$    write(194) ((Krhov(i,j,1),i=nx1_irs,nx2_irs),j=ny1_irs,ny2_irs)
!!$    write(194) ((Krhoe(i,j,1),i=nx1_irs,nx2_irs),j=ny1_irs,ny2_irs)
!!$    close(194)
!!$
!!$    call mpistop('init vortex+comm',0)
!!$    if (iproc==0) print *,'TEST ONE-SIDED COMMS for increments'
    
!!$    ! TEST ONE-SIDED COMMS for 3D increments
!!$    ngh_irs=5
!!$    nghex=8
!!$
!!$    if (BC_face(1,1)%sort>0) ngh_irs(1)=nghex
!!$    if (BC_face(1,2)%sort>0) ngh_irs(2)=nghex
!!$    if (BC_face(2,1)%sort>0) ngh_irs(3)=nghex
!!$    if (BC_face(2,2)%sort>0) ngh_irs(4)=nghex
!!$
!!$    ! initialization of ngh_irs in neighboring domain
!!$    do i=1,3
!!$       do j=1,2
!!$          BC_face(i,j)%ngh_irs=[2,2,2,2,2,2]
!!$       enddo
!!$    enddo
!!$
!!$    if (is_2D) then
!!$       dim=2
!!$    else
!!$       dim=3
!!$    endif
!!$    
!!$    ! send ngh_irs to existing neighbors
!!$    do i=1,dim
!!$       do j=1,2
!!$          k=2*(i-1)+j
!!$          if (neighbor(k)>=0) call MPI_SEND(ngh_irs,6,MPI_INTEGER,neighbor(k),tag,COMM_global,info)
!!$       enddo
!!$    enddo
!!$
!!$    ! receive ngh_irs from existing neighbors
!!$    do i=1,dim
!!$       do j=1,2
!!$          k=2*(i-1)+j
!!$          if (neighbor(k)>=0) call MPI_RECV(BC_face(i,j)%ngh_irs,6,MPI_INTEGER,neighbor(k),tag,COMM_global,status,info)
!!$       enddo
!!$    enddo
!!$
!!$    !if (iproc==1) then
!!$    !   do i=1,dim
!!$    !      do j=1,2
!!$    !         k=2*(i-1)+j
!!$    !         print*,iproc,'iproc',BC_face(i,j)%ngh_irs,k
!!$    !      enddo
!!$    !   enddo
!!$    !endif
!!$         
!!$    nx1_irs= 1-ngh_irs(1)
!!$    nx2_irs=nx+ngh_irs(2)
!!$    ny1_irs= 1-ngh_irs(3)
!!$    ny2_irs=ny+ngh_irs(4)
!!$    nz1_irs= 1-ngh_irs(5)
!!$    nz2_irs=nz+ngh_irs(6)
!!$    if (is_2d) then
!!$       nz1_irs=1
!!$       nz2_irs=1
!!$    endif
!!$    ! Reallocate increments arrays with extended ghost points
!!$    !!deallocate(Krho,Krhou,Krhov,Krhow,Krhoe)
!!$    allocate( Krho(nx1_irs:nx2_irs,ny1_irs:ny2_irs,nz1_irs:nz2_irs))
!!$    allocate(Krhou(nx1_irs:nx2_irs,ny1_irs:ny2_irs,nz1_irs:nz2_irs))
!!$    allocate(Krhov(nx1_irs:nx2_irs,ny1_irs:ny2_irs,nz1_irs:nz2_irs))
!!$    allocate(Krhow(nx1_irs:nx2_irs,ny1_irs:ny2_irs,nz1_irs:nz2_irs))
!!$    allocate(Krhoe(nx1_irs:nx2_irs,ny1_irs:ny2_irs,nz1_irs:nz2_irs))
!!$    ! Definition of MPI types for communications
!!$    !call mpi_types_comm_inc
!!$    call mpi_types_comm1_inc
!!$   
!!$    ! Define windows on allocated arrays (for one-sided comms)
!!$    call mpi_win_comm1_inc       
!!$     Krho=rho_ref
!!$    Krhou=0.0_wp
!!$    Krhov=-u_ref
!!$    Krhoe=p_ref
!!$     Krho(1:nx,1:ny,1:nz)=rho(1:nx,1:ny,1:nz)
!!$    Krhou(1:nx,1:ny,1:nz)= uu(1:nx,1:ny,1:nz)
!!$    Krhov(1:nx,1:ny,1:nz)= vv(1:nx,1:ny,1:nz)
!!$    Krhoe(1:nx,1:ny,1:nz)=prs(1:nx,1:ny,1:nz)
!!$    !call communication_inc(Krho,Krhou,Krhov,Krhow,Krhoe)
!!$    call communication1_inc          
!!$    open(194,file='init1_p'//trim(numchar(iproc))//'_ex.bin',form='unformatted',status='unknown')
!!$    rewind(194)
!!$    write(194) ngh_irs
!!$    write(194) is_swapij2
!!$    write(194) is_rev2(:,1)
!!$    write(194) is_rev2(:,2)
!!$    write(194) nx+ngh_irs(1)+ngh_irs(2)
!!$    write(194) ny+ngh_irs(3)+ngh_irs(4)
!!$    write(194) nz+ngh_irs(5)+ngh_irs(6)
!!$    k=nz
!!$    write(194) ((Krho(i,j,k),i=nx1_irs,nx2_irs),j=ny1_irs,ny2_irs)
!!$    write(194) ((Krhou(i,j,k),i=nx1_irs,nx2_irs),j=ny1_irs,ny2_irs)
!!$    write(194) ((Krhov(i,j,k),i=nx1_irs,nx2_irs),j=ny1_irs,ny2_irs)
!!$    write(194) ((Krhoe(i,j,k),i=nx1_irs,nx2_irs),j=ny1_irs,ny2_irs)
!!$    if (iproc==1) then
!!$       j=ny/2
!!$       print *,'iproc',iproc,j,nx,nz
!!$       write(194) ((Krho(i,j,k),i=nx1_irs,nx2_irs),k=nz1_irs,nz2_irs)
!!$       write(194) ((Krhou(i,j,k),i=nx1_irs,nx2_irs),k=nz1_irs,nz2_irs)
!!$       write(194) ((Krhov(i,j,k),i=nx1_irs,nx2_irs),k=nz1_irs,nz2_irs)
!!$       write(194) ((Krhoe(i,j,k),i=nx1_irs,nx2_irs),k=nz1_irs,nz2_irs)
!!$    endif
!!$    if (iproc==0) then
!!$       i=nx/2
!!$       print *,'iproc',iproc,i,ny,nz
!!$       write(194) ((Krho(i,j,k),j=ny1_irs,ny2_irs),k=nz1_irs,nz2_irs)
!!$       write(194) ((Krhou(i,j,k),j=ny1_irs,ny2_irs),k=nz1_irs,nz2_irs)
!!$       write(194) ((Krhov(i,j,k),j=ny1_irs,ny2_irs),k=nz1_irs,nz2_irs)
!!$       write(194) ((Krhoe(i,j,k),j=ny1_irs,ny2_irs),k=nz1_irs,nz2_irs)
!!$    endif
!!$    close(194)
!!$
!!$    call mpistop('init vortex+comm',0)
    
    ! Thermodynamic variables
    ! =======================
    ! we need Tmp and rho instead of prs
    do k=1,nz
       do j=1,ny
          do i=1,nx
             Tmp(i,j,k)= prs(i,j,k)/p_ref*T_ref
             rho(i,j,k)= rocalc_pt(prs(i,j,k),Tmp(i,j,k),Tmp(i,j,k))
          enddo
       enddo
    enddo

    !if (iproc.eq.0) write (*,*) 'Add vortex OK'

  end subroutine add_vortex
  
  !===============================================================================
  subroutine init_vel_act
  !===============================================================================
    !> Initialization of velocity field inside the actuator
  !===============================================================================
    use mod_block
    use warnstop 
    use mod_coeff_deriv
    implicit none
    ! ----------------------------------------------------------------------------
    ! integer :: i,j,k
    ! real(wp) :: uksi,ueta
    ! ----------------------------------------------------------------------------

    ! set velocity as zero everywhere
    uu = 0.0_wp
    vv = 0.0_wp
    ww = 0.0_wp
    
!     ! FOR NOW (TO BE CHANGED)
!     vv=u_ref
!     uksi=0.0_wp
!     ueta=u_ref
!     if(nob(iproc)==9) then
!        uksi=u_ref
!        ueta=0.0_wp
!     elseif(nob(iproc)==10) then
!        uksi=-u_ref
!        ueta=0.0_wp
!     endif
!     
!     ! orient the flow direction parallel to each wall inside the actuator
!     do k=1,nz
!        do j=1,ny
!           do i=1,nx
!              uu(i,j,k) = uksi*x_ksi(i,j) + ueta*x_eta(i,j)
!              vv(i,j,k) = uksi*y_ksi(i,j) + ueta*y_eta(i,j)
!           enddo
!        enddo
!     enddo
!   
!     ! Damp velocity near the actuator walls
!     ! ====================================
!     ! (to avoid violation of no-slip condition)
!     if (is_bc_wall(1,1)) then ! imin face
!        do k=1,nz
!           do j=1,ny
!              do i=1,10
!                 uu(i,j,k)=uu(i,j,k)/10.0_wp*dble(i-1)
!                 vv(i,j,k)=vv(i,j,k)/10.0_wp*dble(i-1)
!              enddo
!           enddo
!        enddo
!     endif
!     if (is_bc_wall(1,2)) then ! imax face
!        do k=1,nz
!           do j=1,ny
!              do i=nx-9,nx
!                 uu(i,j,k)=uu(i,j,k)/10.0_wp*dble(nx-i)
!                 vv(i,j,k)=vv(i,j,k)/10.0_wp*dble(nx-i)
!              enddo
!           enddo
!        enddo
!     endif
!     if (is_bc_wall(2,1)) then ! jmin face
!        do k=1,nz
!           do j=1,10
!              do i=1,nx
!                 uu(i,j,k)=uu(i,j,k)/10.0_wp*dble(j-1)
!                 vv(i,j,k)=vv(i,j,k)/10.0_wp*dble(j-1)
!              enddo
!           enddo
!        enddo
!     endif
!     if (is_bc_wall(2,2)) then ! jmax face
!        do k=1,nz
!           do j=ny-9,ny
!              do i=1,nx
!                 uu(i,j,k)=uu(i,j,k)/10.0_wp*dble(ny-j)
!                 vv(i,j,k)=vv(i,j,k)/10.0_wp*dble(ny-j)
!              enddo
!           enddo
!        enddo
!     endif
    
  end subroutine init_vel_act

  !===============================================================================
  subroutine init_vel_vane
  !===============================================================================
    !> Initialization of velocity field inside the vane configuration
  !===============================================================================
    use mod_block
    use warnstop
    use mod_coeff_deriv
    use mod_utils
    use mod_eos
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,ig,jg
    integer :: fh
    logical :: iexist
    integer, dimension(MPI_STATUS_SIZE) :: statut
    real(wp), dimension(:,:), allocatable :: ro_rans,rou_rans,rov_rans,prs_rans
    ! ----------------------------------------------------------------------------

    ! ! set velocity as zero everywhere
    ! uu = 0.0_wp
    ! vv = 0.0_wp
    ww = 0.0_wp

    ! Temporary arrays for the RANS block solution
    allocate( ro_rans(ngx,ngy))
    allocate(rou_rans(ngx,ngy))
    allocate(rov_rans(ngx,ngy))
    allocate(prs_rans(ngx,ngy))

    ! Read from RANS initialisation
    inquire(file='RANS_sol/RANS_start_bl'//trim(numchar(nob(iproc)))//'.bin', exist=iexist)
    if (.not.iexist) then
       call mpistop('RANS initialisation files does not exist!',0)
    endif

    ! Reading of RANS file
    ! --------------------
    call MPI_FILE_OPEN(COMM_global,'RANS_sol/RANS_start_bl'//trim(numchar(nob(iproc)))//'.bin', &
      MPI_MODE_RDONLY, MPI_INFO_NULL,fh,info)

    ! Lecture de var
    call MPI_FILE_READ(fh,ro_rans,size(ro_rans),&
         MPI_DOUBLE_PRECISION,statut,info)
    call MPI_FILE_READ(fh,rou_rans,size(rou_rans),&
         MPI_DOUBLE_PRECISION,statut,info)
    call MPI_FILE_READ(fh,rov_rans,size(rov_rans),&
         MPI_DOUBLE_PRECISION,statut,info)
    call MPI_FILE_READ(fh,prs_rans,size(prs_rans),&
         MPI_DOUBLE_PRECISION,statut,info)

    ! Fermeture du fichier
    call MPI_FILE_CLOSE(fh,info)

    ! Initialisation of the conservative variables
    do j=1,ny
       jg = j + coord(2)*ny
       do i=1,nx
          ig = i + coord(1)*nx
          rho(i,j,:)  =  ro_rans(ig,jg)
          ! rou(i,j,:) = rou_rans(ig,jg)
          ! rov(i,j,:) = rov_rans(ig,jg)
          prs(i,j,:) = prs_rans(ig,jg)
       enddo
    enddo

    ! Computation of primitives variables
    do j=1,ny
       jg = j + coord(2)*ny
       do i=1,nx
          ig = i + coord(1)*nx
          uu(i,j,:) = rou_rans(ig,jg)/ro_rans(ig,jg)
          vv(i,j,:) = rov_rans(ig,jg)/ro_rans(ig,jg)
       enddo
    enddo

    ! ! compute pressure from rhoe and rho
    ! do j=1,ny
    !    jg = j + coord(2)*ny
    !    do i=1,nx
    !       ig = i + coord(1)*nx
    !       prs(i,j,:) = pcalc_roero(prs_rans(ig,jg),ro_rans(ig,jg),T_ref)
    !    enddo
    ! enddo

    ! compute temperature from pressure and rho
    do j=1,ny
       jg = j + coord(2)*ny
       do i=1,nx
          ig = i + coord(1)*nx
          ! Tmp(i,j,:) = tcalc_roero(prs_rans(ig,jg),ro_rans(ig,jg),T_ref)
          Tmp(i,j,:) = tcalc_pro(prs_rans(ig,jg),ro_rans(ig,jg),T_ref)
       enddo
    enddo

    ! Deallocation of temporary array
    deallocate(ro_rans,rou_rans,rov_rans,prs_rans)

  end subroutine init_vel_vane

  !===============================================================================
  subroutine init_vel_TE
  !===============================================================================
    !> Initialization of velocity field inside the TE configuration
  !===============================================================================
    ! use mod_block
    use warnstop
    use mod_utils
    use mod_eos
    use mod_io
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k
    logical :: iexist
    character(len=30) :: interpfile='interp_restart/restart_bl1.bin'
    ! ----------------------------------------------------------------------------

    interpfile = 'interp_restart/restart_bl'//trim(numchar(nob(iproc)))//'.bin'
    inquire(file=interpfile, exist=iexist)

    call init_thermo

    ! Read from restart interpolation
    ! -------------------------------
    if (iexist) then

       call read_write_volume(interpfile,"None",READ)

       ! Damp velocity near the wall
       ! ===========================
       ! (to avoid violation of non-penetrability condition)
       if (is_bc_wall(1,1)) then
          rhou(1,:,:)  = 0.0_wp
          rhov(1,:,:)  = 0.0_wp
          rhow(1,:,:)  = 0.0_wp
       endif
       if (is_bc_wall(1,2)) then
          rhou(nx,:,:) = 0.0_wp
          rhov(nx,:,:) = 0.0_wp
          rhow(nx,:,:) = 0.0_wp
       endif
       if (is_bc_wall(2,1)) then
          rhou(:,1,:)  = 0.0_wp
          rhov(:,1,:)  = 0.0_wp
          rhow(:,1,:)  = 0.0_wp
       endif
       if (is_bc_wall(2,2)) then
          rhou(:,ny,:) = 0.0_wp
          rhov(:,ny,:) = 0.0_wp
          rhow(:,ny,:) = 0.0_wp
       endif
       if (is_bc_wall(3,1)) then
          rhou(:,:,1)  = 0.0_wp
          rhov(:,:,1)  = 0.0_wp
          rhow(:,:,1)  = 0.0_wp
       endif
       if (is_bc_wall(3,2)) then
          rhou(:,:,nz) = 0.0_wp
          rhov(:,:,nz) = 0.0_wp
          rhow(:,:,nz) = 0.0_wp
       endif

       ! if presence of stagnation point, /!\ negative value with interp
       do k=1,nz
          do j=1,ny
             do i=1,nx
                if (rho(i,j,k).le.0.0_wp) then
                   rho(i,j,k) = rho(i,j+1,k)
                   rhoe(i,j,k)= rhoe(i,j+1,k)
                endif
             enddo
          enddo
       enddo

    ! Initialize from scratch
    ! -----------------------
    else
       ! Initialization of thermodynamic variables to reference values
       ! =============================================================
       Tmp= T_ref
       prs= p_ref
       rho= rho_ref

       ! Compute velocity components
       ! ===========================
       vv = 0.0_wp
       uu = u_ref

       ! Damp velocity near the wall
       ! ===========================
       ! (to avoid violation of non-penetrability condition)
       if (is_bc_wall(1,1)) then
          do k=1,nz
             do j=1,ny
                do i=1,5
                   uu(i,j,k)=uu(i,j,k)/30.0_wp*dble(i-1)
                enddo
             enddo
          enddo
       endif
       if (is_bc_wall(2,1)) then
          do k=1,nz
             do j=1,5
                do i=1,nx
                   uu(i,j,k)=uu(i,j,k)/30.0_wp*dble(j-1)
                enddo
             enddo
          enddo
       endif
    endif


  end subroutine init_vel_TE

end module mod_init_flow
