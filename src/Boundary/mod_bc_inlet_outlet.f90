!================================================================================
module mod_bc_inlet_outlet
!================================================================================
  !> author: XG,OY,CM
  !> date: jan-mar 2023
  !> 1/ Inlet boundary condition based on Riemann invariants
  !> 2/ Outlet boundary condition imposing back pressure
!================================================================================
  use mod_flow      ! for: flow variables
  use mod_eos       ! for: comp. thermodynamic quantities
  use mod_bc_rad_eq ! for: radial equilibirum approximation
  implicit none
  !------------------------------------------------------------------------------
  ! total (stagnation) quantities
  real(wp) :: p_tot,T_tot,rho_tot
  real(wp), dimension(:,:), allocatable :: H_tot,s_tot
  real(wp) :: H_to,s_to ! for compatibility with bc_inlet_test_imin
  ! direction relative to (x,y,z) coordinates 
  real(wp), dimension(:,:,:), allocatable :: db
  ! sign correction to have inwards inlet velocity
  real(wp) :: sign_imin,sign_imax,sign_jmin,sign_jmax,sign_kmin,sign_kmax
  ! flow direction
  real(wp), dimension(:,:), allocatable :: cdir_imin,cdir_imax
  real(wp), dimension(:,:), allocatable :: cdir_jmin,cdir_jmax
  real(wp), dimension(:,:), allocatable :: cdir_kmin,cdir_kmax
  ! normal component of velocity for Riemann invariants
  real(wp) :: Uin,Ubn
  !------------------------------------------------------------------------------

contains

  !==============================================================================
  subroutine init_bc_inlet_outlet
  !==============================================================================
    !> Initialization of Inlet boundary conditions
    !> [only 2D BC, extruded along z]
  !==============================================================================
    use mod_mpi ! for iproc (print screen)
    use mod_fluid ! for gam,gam1 ! TEMP for PFG
    use warnstop ! TEMP for check
    implicit none
    !----------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: dbnb,gsgm1
    ! local velocity components & Mach number
    real(wp) :: ui,vi,wi,Unorm,Mi
    ! local direction angles of velocity vector
    real(wp) :: theta_vel,phi_vel
    !----------------------------------------------------------------------------

    ! Constant direction of the velocity at the boundary
    ! ==================================================
    
    ! ===================================================
    !                     /!\ CAREFUL /!\
    ! theta_ref, phi_ref & p_exit now given in param.ini
    ! Directly transformed in radians in read_param
    ! ===================================================

    ! 2D flow angle TO BE CHANGED (should be from param.ini & setupref)
    !theta_ref = 29.92_wp
    !theta_ref= 0.0_wp
    ! 3D flow angle TO BE CHANGED (should be from param.ini & setupref)
    ! theta_ref = 30.0_wp
    ! phi_ref= 0.0_wp

    ! Set exit pressure for back-pressure BC TO BE CHANGED (should be from param.ini & setupref)
    ! ======================================
    !p_exit = 10943.65_wp
    !p_exit = p_ref*0.528281787717174_wp
    !p_exit = p_ref*0.534460302143549_wp
    !p_exit = p_ref*0.5373_wp
    !p_exit = p_ref*0.423153711961457_wp
    !p_exit = p_ref*0.655370258341356_wp
    ! case M2is=0.7
    !p_exit = p_ref*0.720927860991916_wp
    ! case M2is=0.624
    !p_exit = p_ref*0.769148719011164_wp
    !p_exit = p_ref*0.523_wp
    !p_exit= p_ref*0.51_wp
    !p_exit= p_ref/1.82_wp
    ! calculs Novec up to 28/08
    !!p_exit = p_ref*0.78_wp
    !p_exit = p_ref*0.80_wp
    ! p_exit = p_ref*0.79_wp
    !!p_exit = p_ref*0.85_wp
    !!p_exit = p_ref*0.93_wp
    !!p_exit = p_ref*0.96_wp
    !!p_exit = p_ref*0.98_wp
    ! cases Graz University
    ! ---------------------
    ! case M2is=0.7
    !p_exit = p_ref*0.720927860991916_wp
    !p_exit = p_ref*0.73_wp
    ! case M2is=0.624
    !p_exit = p_ref*0.769148719011164_wp
    !p_exit = p_ref*0.77_wp

    ! test-case Conical Nozzle (conv3)
    ! ------------------------
    !p_exit= p_ref*0.8_wp
    
    ! test-case Converging-Diverging Nozzle NASA
    ! ------------------------------------------
    !!!p_exit= p_ref ! subsonic, isentropic flow (pexit/pt=0.89)
    !p_exit= p_ref*0.89_wp ! subsonic, isentropic flow (pexit/pt=0.89)
    !p_exit= p_ref*0.75_wp ! supersonic flow with a normal shock in the diffusing section (pexit/pt=0.75)
    !p_exit= p_ref*0.16_wp ! supersonic, isentropic flow (pexit/pt=0.16)

    ! annular duct
    ! ------------
    !p_exit= p_ref

    ! LS89 vane
    ! ---------
    !p_exit= p_ref/5.12_wp

    ! Baumgartner
    ! -----------
    !p_exit= p_ref/4.35_wp
    !p_exit= p_ref/2.0_wp

    ! idealized blade vane configuration
    ! ------------------------------------------
    ! p_exit= p_ref*0.25_wp ! config passman (pexit/pt=0.25)
    !p_exit= p_ref*0.2_wp ! config passman (pexit/pt=0.25)
      
    ! test-case TROVA Nozzle
    ! ----------------------
    !p_exit= 200000.0_wp
    
    ! Init procedure for radial equilibrium (axis is Ox)
    ! =====================================
    if (is_rea) call init_rea

    ! Define unitary normal vector to the BC face
    ! ===========================================
    ! normals are computed in Grid/grid_normals.f90
    ! /!\ normals from generalized coordinates
    !     grad(ksi)/||grad(ksi)|| for ksi=cst BC
    !     grad(eta)/||grad(eta)|| for eta=cst BC
    ! /!\ directon depends on (i,j) orientation
    !     NEEDS to set INWARDS normals for inlet
    !     not necessary for outlet (only direction)

    ! Nota: for Inlet BCs: INWARDS normals
    ! =====
    ! -> by taking scalar product with inlet direction at ny/2 (or nx/2)
    !    correct_sign=+1 if (db.nb)>0
    !    correct_sign=-1 if (db.nb)<0
    
    ! Inlet BC at imin
    ! ================
    if (BC_face(1,1)%sort==-4) then

       ! Allocations
       ! -----------
       ! total enthalpy and entropy
       allocate(H_tot(ny,nz),s_tot(ny,nz))
       ! inlet velocity direction vector
       allocate(db(ny,nz,3))

       ! Set total (stagnation) quantities & imposed velocity direction vector
       ! ---------------------------------------------------------------------
       if ((BC_face(1,1)%is_mean_ref).and.(.not.is_stagnation)) then

          ! compute total quantities and direction from local variables
          ! -----------------------------------------------------------
          gsgm1=gam/gam1 ! exponent gam/(gam-1)

          do k=1,nz
             do j=1,ny
                ! local velocity components
                ui=BC_face(1,1)%Uref(j,k,2)
                vi=BC_face(1,1)%Uref(j,k,3)
                wi=BC_face(1,1)%Uref(j,k,4)
                ! local velocity norm
                Unorm=sqrt(ui**2+vi**2+wi**2)
                ! local Mach number from inlet velocities
                Mi=Unorm/c_ref

                ! total temperature /!\ PFG (TEMP -> generalization dense gas)
                T_tot=T_ref*(1.0_wp+0.5_wp*gam1*Mi**2)
       
                ! total pressure /!\ PFG (TEMP -> generalization dense gas)
                p_tot=p_ref*(1.0_wp+0.5_wp*gam1*Mi**2)**gsgm1

                rho_tot=rocalc_pt(p_tot,T_tot,rho_ref)
                H_tot(j,k)= ecalc_tro(T_tot,rho_tot) + p_tot/rho_tot
                s_tot(j,k)= scalc_tro(T_tot,rho_tot)
                
                ! local flow angles
                theta_vel=asin(vi/Unorm)
                  phi_vel=asin(wi/Unorm)
                
!!$                if (j==ny) then
!!$                   print *,'mach',Mi
!!$                   print *,'tot',T_tot,p_tot,rho_tot,H_tot(j,k),s_tot(j,k)
!!$                   print *,'angles',theta_vel*180./pi,phi_vel*180./pi
!!$                endif
                  
                ! direction relative to (x,y,z) coordinates
                dbnb=sqrt(cos(theta_vel)**2-sin(phi_vel)**2)
                db(j,k,:)= (/dbnb,sin(theta_vel),sin(phi_vel)/)
                
!!$                if (j==ny) then
!!$                   print *,'dir',db(j,k,:)
!!$                endif
             enddo
          enddo
!!$       else if ((BC_face(1,1)%is_mean_ref).and.(is_stagnation)) then
!!$          do k=1,nz
!!$             do j=1,ny
!!$                ! local velocity components
!!$                ui=BC_face(1,1)%Uref(j,k,2)
!!$                vi=BC_face(1,1)%Uref(j,k,3)
!!$                wi=BC_face(1,1)%Uref(j,k,4)
!!$                ! local velocity norm
!!$                Unorm=sqrt(ui**2+vi**2+wi**2)
!!$                if ((Unorm.eq.0).and.(k.eq.1).and.(coord(3).eq.0)) print *,"/!\ Unorm=0 bc inlet_outlet",nob(iproc),j+coord(2)*ny,k+coord(3)*nz
!!$                if (Unorm.eq.0) Unorm=1.0_wp
!!$
!!$                ! Reference H_tot & s
!!$                H_tot(j,k)=BC_face(1,1)%Uref(j,k,1)
!!$                s_tot(j,k)=BC_face(1,1)%Uref(j,k,5)
!!$
!!$                ! local flow angles
!!$                theta_vel=asin(vi/Unorm)
!!$                  phi_vel=asin(wi/Unorm)
!!$
!!$                T_tot=T_ref
!!$                p_tot=p_ref
!!$                rho_tot=rocalc_pt(p_tot,T_tot,rho_ref)
!!$
!!$                ! direction relative to (x,y,z) coordinates
!!$                dbnb=sqrt(cos(theta_vel)**2-sin(phi_vel)**2)
!!$                db(j,k,:)= (/dbnb,sin(theta_vel),sin(phi_vel)/)
!!$             enddo
!!$          enddo

       else
          ! imposed total quantities and direction from param.ini
          ! -----------------------------------------------------

          ! total quantities
          T_tot=T_ref
          p_tot=p_ref
          rho_tot= rocalc_pt(p_tot,T_tot,rho_ref)
          H_tot= ecalc_tro(T_tot,rho_tot) + p_tot/rho_tot
          s_tot= scalc_tro(T_tot,rho_tot)
          
          ! constant values & convert angles in radians
          ! ==============================================
          !               /!\ CAREFUL /!\
          ! theta_ref & phi_ref now given in param.ini
          ! Directly transformed in radians in read_param
          ! ==============================================
          ! theta_ref=theta_ref*pi/180.0_wp
          !   phi_ref=  phi_ref*pi/180.0_wp
          ! direction relative to (x,y,z) coordinates
          dbnb=sqrt(cos(theta_ref)**2-sin(phi_ref)**2)
          do k=1,nz
             do j=1,ny
                db(j,k,:)= (/dbnb,sin(theta_ref),sin(phi_ref)/)                
             enddo
          enddo
       endif

       ! Projection on normal direction
       ! ------------------------------       
       allocate(cdir_imin(ny,nz))

       ! determine correct_sign at (j0=ny/2,k0=nz/2)
       dbnb=db(j0,k0,1)*nxn_imin(j0,k0)+db(j0,k0,2)*nyn_imin(j0,k0)+db(j0,k0,3)*nzn_imin(j0,k0)
       sign_imin=dbnb/abs(dbnb)

       ! flow direction relative to the boundary
       do k=1,nz
          do j=1,ny
             cdir_imin(j,k)=sign_imin/(db(j,k,1)*nxn_imin(j,k)+db(j,k,2)*nyn_imin(j,k)+db(j,k,3)*nzn_imin(j,k))
          enddo
       enddo
    endif

    ! Inlet BC at imax
    ! ================
    if (BC_face(1,2)%sort==-4) then
       
       ! Allocations
       ! -----------
       ! total enthalpy and entropy
       allocate(H_tot(ny,nz),s_tot(ny,nz))
       ! inlet velocity direction vector
       allocate(db(ny,nz,3))

       ! Set total (stagnation) quantities & imposed velocity direction vector
       ! ---------------------------------------------------------------------
       if ((BC_face(1,2)%is_mean_ref).and.(.not.is_stagnation)) then

          ! compute total quantities and direction from local variables
          ! -----------------------------------------------------------
          gsgm1=gam/gam1 ! exponent gam/(gam-1)

          do k=1,nz
             do j=1,ny
                ! local velocity components
                ui=BC_face(1,2)%Uref(j,k,2)
                vi=BC_face(1,2)%Uref(j,k,3)
                wi=BC_face(1,2)%Uref(j,k,4)
                ! local velocity norm
                Unorm=sqrt(ui**2+vi**2+wi**2)
                ! local Mach number from inlet velocities
                Mi=Unorm/c_ref

                ! total temperature /!\ PFG (TEMP -> generalization dense gas)
                T_tot=T_ref*(1.0_wp+0.5_wp*gam1*Mi**2)
       
                ! total pressure /!\ PFG (TEMP -> generalization dense gas)
                p_tot=p_ref*(1.0_wp+0.5_wp*gam1*Mi**2)**gsgm1

                rho_tot=rocalc_pt(p_tot,T_tot,rho_ref)
                H_tot(j,k)= ecalc_tro(T_tot,rho_tot) + p_tot/rho_tot
                s_tot(j,k)= scalc_tro(T_tot,rho_tot)
                
                ! local flow angles
                theta_vel=asin(vi/Unorm)
                  phi_vel=asin(wi/Unorm)

                ! direction relative to (x,y,z) coordinates
                dbnb=sqrt(cos(theta_vel)**2-sin(phi_vel)**2)
                db(j,k,:)= (/dbnb,sin(theta_vel),sin(phi_vel)/)                
             enddo
          enddo
       else
          ! imposed total quantities and direction from param.ini
          ! -----------------------------------------------------

          ! total quantities
          T_tot=T_ref
          p_tot=p_ref
          rho_tot= rocalc_pt(p_tot,T_tot,rho_ref)
          H_tot= ecalc_tro(T_tot,rho_tot) + p_tot/rho_tot
          s_tot= scalc_tro(T_tot,rho_tot)
          
          ! constant values & convert angles in radians
          ! ==============================================
          !               /!\ CAREFUL /!\
          ! theta_ref & phi_ref now given in param.ini
          ! Directly transformed in radians in read_param
          ! ==============================================
          ! theta_ref=theta_ref*pi/180.0_wp
          !   phi_ref=  phi_ref*pi/180.0_wp
          ! direction relative to (x,y,z) coordinates
          dbnb=sqrt(cos(theta_ref)**2-sin(phi_ref)**2)
          do k=1,nz
             do j=1,ny
                db(j,k,:)= (/dbnb,sin(theta_ref),sin(phi_ref)/)                
             enddo
          enddo
       endif

       ! Projection on normal direction
       ! ------------------------------       
       allocate(cdir_imax(ny,nz))
       
       ! determine correct_sign (at j0=ny/2)
       dbnb=db(j0,k0,1)*nxn_imax(j0,k0)+db(j0,k0,2)*nyn_imax(j0,k0)+db(j0,k0,3)*nzn_imax(j0,k0)
       sign_imax=dbnb/abs(dbnb)
       
       ! flow direction relative to the boundary
       do k=1,nz
          do j=1,ny
             cdir_imax(j,k)=sign_imax/(db(j,k,1)*nxn_imax(j,k)+db(j,k,2)*nyn_imax(j,k)+db(j,k,3)*nzn_imax(j,k))
          enddo
       enddo       
    endif

    ! Inlet BC at jmin
    ! ================
    if (BC_face(2,1)%sort==-4) then
       
       ! Allocations
       ! -----------
       ! total enthalpy and entropy
       allocate(H_tot(nx,nz),s_tot(nx,nz))
       ! inlet velocity direction vector
       allocate(db(nx,nz,3))

       ! Set total (stagnation) quantities & imposed velocity direction vector
       ! ---------------------------------------------------------------------
       if ((BC_face(2,1)%is_mean_ref).and.(.not.is_stagnation)) then

          ! compute total quantities and direction from local variables
          ! -----------------------------------------------------------
          gsgm1=gam/gam1 ! exponent gam/(gam-1)

          do k=1,nz
             do i=1,nx
                ! local velocity components
                ui=BC_face(2,1)%Uref(i,k,2)
                vi=BC_face(2,1)%Uref(i,k,3)
                wi=BC_face(2,1)%Uref(i,k,4)
                ! local velocity norm
                Unorm=sqrt(ui**2+vi**2+wi**2)
                ! local Mach number from inlet velocities
                Mi=Unorm/c_ref

                ! total temperature /!\ PFG (TEMP -> generalization dense gas)
                T_tot=T_ref*(1.0_wp+0.5_wp*gam1*Mi**2)
       
                ! total pressure /!\ PFG (TEMP -> generalization dense gas)
                p_tot=p_ref*(1.0_wp+0.5_wp*gam1*Mi**2)**gsgm1

                rho_tot=rocalc_pt(p_tot,T_tot,rho_ref)
                H_tot(i,k)= ecalc_tro(T_tot,rho_tot) + p_tot/rho_tot
                s_tot(i,k)= scalc_tro(T_tot,rho_tot)
                
                ! local flow angles
                theta_vel=asin(vi/Unorm)
                  phi_vel=asin(wi/Unorm)
                
                ! direction relative to (x,y,z) coordinates
                dbnb=sqrt(cos(theta_vel)**2-sin(phi_vel)**2)
                db(i,k,:)= (/dbnb,sin(theta_vel),sin(phi_vel)/)                
             enddo
          enddo
       else
          ! imposed total quantities and direction from param.ini
          ! -----------------------------------------------------

          ! total quantities
          T_tot=T_ref
          p_tot=p_ref
          rho_tot= rocalc_pt(p_tot,T_tot,rho_ref)
          H_tot= ecalc_tro(T_tot,rho_tot) + p_tot/rho_tot
          s_tot= scalc_tro(T_tot,rho_tot)
          
          ! constant values & convert angles in radians
          ! ==============================================
          !               /!\ CAREFUL /!\
          ! theta_ref & phi_ref now given in param.ini
          ! Directly transformed in radians in read_param
          ! ==============================================
          ! theta_ref=theta_ref*pi/180.0_wp
          !   phi_ref=  phi_ref*pi/180.0_wp
          ! direction relative to (x,y,z) coordinates
          dbnb=sqrt(cos(theta_ref)**2-sin(phi_ref)**2)
          do k=1,nz
             do i=1,nx
                db(i,k,:)= (/dbnb,sin(theta_ref),sin(phi_ref)/)                
             enddo
          enddo
       endif

       ! Projection on normal direction
       ! ------------------------------       
       allocate(cdir_jmin(nx,nz))
       
       ! determine correct_sign at (i0=nx/2,k0=nz/2)
       dbnb=db(i0,k0,1)*nxn_jmin(i0,k0)+db(i0,k0,2)*nyn_jmin(i0,k0)+db(i0,k0,3)*nzn_jmin(i0,k0)
       sign_jmin=dbnb/abs(dbnb)
       
       ! flow direction relative to the boundary
       do k=1,nz
          do i=1,nx
             cdir_jmin(i,k)=sign_jmin/(db(i,k,1)*nxn_jmin(i,k)+db(i,k,2)*nyn_jmin(i,k)+db(i,k,3)*nzn_jmin(i,k))
          enddo
       enddo
    endif

    ! Inlet BC at jmax
    ! ================
    if (BC_face(2,2)%sort==-4) then
       
       ! Allocations
       ! -----------
       ! total enthalpy and entropy
       allocate(H_tot(nx,nz),s_tot(nx,nz))
       ! inlet velocity direction vector
       allocate(db(nx,nz,3))

       ! Set total (stagnation) quantities & imposed velocity direction vector
       ! ---------------------------------------------------------------------
       if ((BC_face(2,2)%is_mean_ref).and.(.not.is_stagnation)) then

          ! compute total quantities and direction from local variables
          ! -----------------------------------------------------------
          gsgm1=gam/gam1 ! exponent gam/(gam-1)

          do k=1,nz
             do i=1,nx
                ! local velocity components
                ui=BC_face(2,2)%Uref(i,k,2)
                vi=BC_face(2,2)%Uref(i,k,3)
                wi=BC_face(2,2)%Uref(i,k,4)
                ! local velocity norm
                Unorm=sqrt(ui**2+vi**2+wi**2)
                ! local Mach number from inlet velocities
                Mi=Unorm/c_ref

                ! total temperature /!\ PFG (TEMP -> generalization dense gas)
                T_tot=T_ref*(1.0_wp+0.5_wp*gam1*Mi**2)
       
                ! total pressure /!\ PFG (TEMP -> generalization dense gas)
                p_tot=p_ref*(1.0_wp+0.5_wp*gam1*Mi**2)**gsgm1

                rho_tot=rocalc_pt(p_tot,T_tot,rho_ref)
                H_tot(i,k)= ecalc_tro(T_tot,rho_tot) + p_tot/rho_tot
                s_tot(i,k)= scalc_tro(T_tot,rho_tot)
                
                ! local flow angles
                theta_vel=asin(vi/Unorm)
                  phi_vel=asin(wi/Unorm)
                
                ! direction relative to (x,y,z) coordinates
                dbnb=sqrt(cos(theta_vel)**2-sin(phi_vel)**2)
                db(i,k,:)= (/dbnb,sin(theta_vel),sin(phi_vel)/)                
             enddo
          enddo
       else
          ! imposed total quantities and direction from param.ini
          ! -----------------------------------------------------

          ! total quantities
          T_tot=T_ref
          p_tot=p_ref
          rho_tot= rocalc_pt(p_tot,T_tot,rho_ref)
          H_tot= ecalc_tro(T_tot,rho_tot) + p_tot/rho_tot
          s_tot= scalc_tro(T_tot,rho_tot)
          
          ! constant values & convert angles in radians
          ! ==============================================
          !               /!\ CAREFUL /!\
          ! theta_ref & phi_ref now given in param.ini
          ! Directly transformed in radians in read_param
          ! ==============================================
          ! theta_ref=theta_ref*pi/180.0_wp
          !   phi_ref=  phi_ref*pi/180.0_wp
          ! direction relative to (x,y,z) coordinates
          dbnb=sqrt(cos(theta_ref)**2-sin(phi_ref)**2)
          do k=1,nz
             do i=1,nx
                db(i,k,:)= (/dbnb,sin(theta_ref),sin(phi_ref)/)                
             enddo
          enddo
       endif

       ! Projection on normal direction
       ! ------------------------------       
       allocate(cdir_jmax(nx,nz))
       
       ! determine correct_sign at (i0=nx/2,k0=nz/2)
       dbnb=db(i0,k0,1)*nxn_jmax(i0,k0)+db(i0,k0,2)*nyn_jmax(i0,k0)+db(i0,k0,3)*nzn_jmax(i0,k0)
       sign_jmax=dbnb/abs(dbnb)
       
       ! flow direction relative to the boundary
       do k=1,nz
          do i=1,nx
             cdir_jmax(i,k)=sign_jmax/(db(i,k,1)*nxn_jmax(i,k)+db(i,k,2)*nyn_jmax(i,k)+db(i,k,3)*nzn_jmax(i,k))
          enddo
       enddo       
    endif

    ! Inlet BC at kmin
    ! ================
    if (BC_face(3,1)%sort==-4) then
             
       ! Allocations
       ! -----------
       ! total enthalpy and entropy
       allocate(H_tot(nx,ny),s_tot(nx,ny))
       ! inlet velocity direction vector
       allocate(db(nx,ny,3))

       ! Set total (stagnation) quantities & imposed velocity direction vector
       ! ---------------------------------------------------------------------
       if ((BC_face(3,1)%is_mean_ref).and.(.not.is_stagnation)) then

          ! compute total quantities and direction from local variables
          ! -----------------------------------------------------------
          gsgm1=gam/gam1 ! exponent gam/(gam-1)

          do j=1,ny
             do i=1,nx
                ! local velocity components
                ui=BC_face(3,1)%Uref(i,j,2)
                vi=BC_face(3,1)%Uref(i,j,3)
                wi=BC_face(3,1)%Uref(i,j,4)
                ! local velocity norm
                Unorm=sqrt(ui**2+vi**2+wi**2)
                ! local Mach number from inlet velocities
                Mi=Unorm/c_ref

                ! total temperature /!\ PFG (TEMP -> generalization dense gas)
                T_tot=T_ref*(1.0_wp+0.5_wp*gam1*Mi**2)
       
                ! total pressure /!\ PFG (TEMP -> generalization dense gas)
                p_tot=p_ref*(1.0_wp+0.5_wp*gam1*Mi**2)**gsgm1

                rho_tot=rocalc_pt(p_tot,T_tot,rho_ref)
                H_tot(i,j)= ecalc_tro(T_tot,rho_tot) + p_tot/rho_tot
                s_tot(i,j)= scalc_tro(T_tot,rho_tot)
                
                ! local flow angles
                theta_vel=asin(vi/Unorm)
                  phi_vel=asin(wi/Unorm)
                
                ! direction relative to (x,y,z) coordinates
                dbnb=sqrt(cos(theta_vel)**2-sin(phi_vel)**2)
                db(i,j,:)= (/dbnb,sin(theta_vel),sin(phi_vel)/)                
             enddo
          enddo
       else
          ! imposed total quantities and direction from param.ini
          ! -----------------------------------------------------

          ! total quantities
          T_tot=T_ref
          p_tot=p_ref
          rho_tot= rocalc_pt(p_tot,T_tot,rho_ref)
          H_tot= ecalc_tro(T_tot,rho_tot) + p_tot/rho_tot
          s_tot= scalc_tro(T_tot,rho_tot)
          
          ! constant values & convert angles in radians
          ! ==============================================
          !               /!\ CAREFUL /!\
          ! theta_ref & phi_ref now given in param.ini
          ! Directly transformed in radians in read_param
          ! ==============================================
          ! theta_ref=theta_ref*pi/180.0_wp
          !   phi_ref=  phi_ref*pi/180.0_wp
          ! direction relative to (x,y,z) coordinates
          dbnb=sqrt(cos(theta_ref)**2-sin(phi_ref)**2)
          do j=1,ny
             do i=1,nx
                db(i,j,:)= (/dbnb,sin(theta_ref),sin(phi_ref)/)                
             enddo
          enddo
       endif

       ! Projection on normal direction
       ! ------------------------------       
       allocate(cdir_kmin(nx,ny))
       
       ! determine correct_sign at (i0=nx/2,j0=ny/2)
       dbnb=db(i0,j0,1)*nxn_kmin(i0,j0)+db(i0,j0,2)*nyn_kmin(i0,j0)+db(i0,j0,3)*nzn_kmin(i0,j0)
       sign_kmin=dbnb/abs(dbnb)
              
       ! flow direction relative to the boundary
       do j=1,ny
          do i=1,nx
             cdir_kmin(i,j)=sign_kmin/(db(i,j,1)*nxn_kmin(i,j)+db(i,j,2)*nyn_kmin(i,j)+db(i,j,3)*nzn_kmin(i,j))
          enddo
       enddo
    endif

    ! Inlet BC at kmax
    ! ================
    if (BC_face(3,2)%sort==-4) then
       
       ! Allocations
       ! -----------
       ! total enthalpy and entropy
       allocate(H_tot(nx,ny),s_tot(nx,ny))
       ! inlet velocity direction vector
       allocate(db(nx,ny,3))

       ! Set total (stagnation) quantities & imposed velocity direction vector
       ! ---------------------------------------------------------------------
       if ((BC_face(3,1)%is_mean_ref).and.(.not.is_stagnation)) then

          ! compute total quantities and direction from local variables
          ! -----------------------------------------------------------
          gsgm1=gam/gam1 ! exponent gam/(gam-1)

          do j=1,ny
             do i=1,nx
                ! local velocity components
                ui=BC_face(3,1)%Uref(i,j,2)
                vi=BC_face(3,1)%Uref(i,j,3)
                wi=BC_face(3,1)%Uref(i,j,4)
                ! local velocity norm
                Unorm=sqrt(ui**2+vi**2+wi**2)
                ! local Mach number from inlet velocities
                Mi=Unorm/c_ref

                ! total temperature /!\ PFG (TEMP -> generalization dense gas)
                T_tot=T_ref*(1.0_wp+0.5_wp*gam1*Mi**2)
       
                ! total pressure /!\ PFG (TEMP -> generalization dense gas)
                p_tot=p_ref*(1.0_wp+0.5_wp*gam1*Mi**2)**gsgm1

                rho_tot=rocalc_pt(p_tot,T_tot,rho_ref)
                H_tot(i,j)= ecalc_tro(T_tot,rho_tot) + p_tot/rho_tot
                s_tot(i,j)= scalc_tro(T_tot,rho_tot)
                
                ! local flow angles
                theta_vel=asin(vi/Unorm)
                  phi_vel=asin(wi/Unorm)
                
                ! direction relative to (x,y,z) coordinates
                dbnb=sqrt(cos(theta_vel)**2-sin(phi_vel)**2)
                db(i,j,:)= (/dbnb,sin(theta_vel),sin(phi_vel)/)                
             enddo
          enddo
       else
          ! imposed total quantities and direction from param.ini
          ! -----------------------------------------------------

          ! total quantities
          T_tot=T_ref
          p_tot=p_ref
          rho_tot= rocalc_pt(p_tot,T_tot,rho_ref)
          H_tot= ecalc_tro(T_tot,rho_tot) + p_tot/rho_tot
          s_tot= scalc_tro(T_tot,rho_tot)
          
          ! constant values & convert angles in radians
          ! ==============================================
          !               /!\ CAREFUL /!\
          ! theta_ref & phi_ref now given in param.ini
          ! Directly transformed in radians in read_param
          ! ==============================================
          ! theta_ref=theta_ref*pi/180.0_wp
          !   phi_ref=  phi_ref*pi/180.0_wp
          ! direction relative to (x,y,z) coordinates
          dbnb=sqrt(cos(theta_ref)**2-sin(phi_ref)**2)
          do j=1,ny
             do i=1,nx
                db(i,j,:)= (/dbnb,sin(theta_ref),sin(phi_ref)/)                
             enddo
          enddo
       endif

       ! Projection on normal direction
       ! ------------------------------       
       allocate(cdir_kmax(nx,ny))
       
       ! determine correct_sign at (i0=nx/2,j0=ny/2)
       dbnb=db(i0,j0,1)*nxn_kmax(i0,j0)+db(i0,j0,2)*nyn_kmax(i0,j0)+db(i0,j0,3)*nzn_kmax(i0,j0)
       sign_kmax=dbnb/abs(dbnb)
       
       ! flow direction relative to the boundary
       do j=1,ny
          do i=1,nx
             cdir_kmax(i,j)=sign_kmax/(db(i,j,1)*nxn_kmax(i,j)+db(i,j,2)*nyn_kmax(i,j)+db(i,j,3)*nzn_kmax(i,j))
          enddo
       enddo       
    endif
    
    if (iproc.eq.0) print *,"=================="
    if (iproc.eq.0) print *,"p_tot",p_tot
    if (iproc.eq.0) print *,"T_tot",T_tot
    if (iproc.eq.0) print *,"rho_tot",rho_tot
    if (iproc.eq.0) print *,"=================="

!!$    if (iproc==0) then
!!$       open(200,file='inlet1.dat',form='formatted',status='replace')
!!$       rewind(200)
!!$    elseif (iproc==1) then
!!$       open(201,file='inlet2.dat',form='formatted',status='replace')
!!$       rewind(201)
!!$    endif

!!$    call mpistop('check in init_bc_inlet_outlet',0)

  end subroutine init_bc_inlet_outlet
    
  !==============================================================================
  subroutine bc_inflow_imin
  !==============================================================================
    !> Subsonic Inflow Boundary Condition at imin: extrapolate Riemann invariant
    !> to impose direction, total pressure and total temperature
  !==============================================================================
    use mod_time !for check
    use mod_RFM
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: j,k
    ! Riemann invariant dp-(roc)_0 dU
    real(wp) :: am0,roc0,coeff
    ! quantities at boundary points
    real(wp) :: rob,Tb,pb
    ! ---------------------------------------------------------------------------
    ! param for fluctuations
    real(wp), dimension(ny,nz) :: u_in,v_in,w_in
    ! ---------------------------------------------------------------------------

    ! Update inlet points
    ! ===================
    if (is_RFM) call disturb_inlet_RFM_turb(u_in,v_in,w_in)

    ! Update inlet points
    ! ===================
    do k=1,nz
       do j=1,ny

          ! Quantities from interior nodes
          ! ==============================
          
          ! normal velocity at interior nodes Uin=Ui.nb
          Uin= uu(2,j,k)*nxn_imin(j,k)+vv(2,j,k)*nyn_imin(j,k)+ww(2,j,k)*nzn_imin(j,k)

          !if ((j==1.or.j==ny).and.(k==1).and.(irk==nrk)) print *,ntime,j,Uin
          
          ! correct sign to have inwards velocity
          Uin=sign_imin*Uin
          
          ! product (roc)_0 (extrapolated from interior)
          roc0=rho_n(2,j,k)*c_(2,j,k)
          !roc0=rho_n(1,j,k)*c_(1,j,k)

          ! compute outgoing characteristic from interior nodes
          ! dp-(roc)_0 dU=cste <=> pb-(roc)_0 Ubn=pi-(roc)_0 Uin=am0
          am0= prs(2,j,k)-roc0*Uin

          ! coefficient taken flow direction into account
          coeff=(cdir_imin(j,k)/roc0)**2

          ! Solve p_b-roc0*U_b-am0=0 with Newton's method
          ! =============================================

          ! initial values (taken as interior values)
          ! --------------
          !Tb =T_tot
          !pb =p_tot
          !rob=rho_tot
          Tb =Tmp(2,j,k)
          pb =prs(2,j,k)
          rob=rho_n(2,j,k)


          ! Solve using approximate Newton's method
          ! ---------------------------------------
          ! input/output: initial and updated thermo. variables
          ! input: total enthalpy and entropy to be conserved
          ! input: interior part of Riemann invariant am0
          ! input: coeff taking into account flow direction
          !print *, 'call tcalc_Hstot',H_tot(j,k),s_tot(j,k)
          !print *, 'init',Tb,rob,pb
          !print *, 'var',Uin,roc0,am0,coeff
          call tcalc_Hstot(Tb,rob,pb,H_tot(j,k),s_tot(j,k),am0,coeff)

          !print *, 'call tcalc_Hstot ok',j,k
          ! update thermo variables
          ! -----------------------
          Tmp(1,j,k)=Tb
          prs(1,j,k)=pb
          rho_n(1,j,k)=rob     

          ! update normal velocity
          ! ----------------------
          Ubn=(pb-am0)/roc0

          ! compute velocity components by projecting modulus
          ! -------------------------------------------------
          if (is_RFM) then
             uu(1,j,k)= Ubn*db(j,k,1)*cdir_imin(j,k) + u_in(j,k)
             vv(1,j,k)= Ubn*db(j,k,2)*cdir_imin(j,k) + v_in(j,k)
             ww(1,j,k)= Ubn*db(j,k,3)*cdir_imin(j,k) + w_in(j,k)
          else
             uu(1,j,k)= Ubn*db(j,k,1)*cdir_imin(j,k)
             vv(1,j,k)= Ubn*db(j,k,2)*cdir_imin(j,k)
             ww(1,j,k)= Ubn*db(j,k,3)*cdir_imin(j,k)
          endif

          ! for annulus case, impose all variables
          !rob=BC_face(1,1)%Uref(j,k,1)
          !rho_n(1,j,k)=BC_face(1,1)%Uref(j,k,1)
          !uu(1,j,k)=BC_face(1,1)%Uref(j,k,2)
          !vv(1,j,k)=BC_face(1,1)%Uref(j,k,3)
          !ww(1,j,k)=BC_face(1,1)%Uref(j,k,4)
          !prs(1,j,k)=BC_face(1,1)%Uref(j,k,5)
          !Tmp(1,j,k)=prs(1,j,k)/rg/rob

!!$          if (k==nz/2) then
!!$             print *,j,uu(1,j,k)
!!$          endif
          !if ((j==1.or.j==ny).and.(k==1).and.(irk==nrk)) print *,ntime,j,uu(1,j,k),vv(1,j,k),ww(1,j,k),pb,rob,Tb
                 
!!$          if ((j==ny/2).and.(k==nz/2).and.(irk==nrk)) then
!!$             if (iproc==0) then
!!$                write(200,'(i8,1x,f15.10,1x,f17.5,1x,f15.10,1x)') ntime,uu(1,j,k),pb,Tb
!!$             elseif (iproc==1) then
!!$                write(201,'(i8,1x,f15.10,1x,f17.5,1x,f15.10,1x)') ntime,uu(1,j,k),pb,Tb
!!$             endif
!!$          endif
                
          ! update conservative variables
          ! -----------------------------
          rhou_n(1,j,k)= rob*uu(1,j,k)
          rhov_n(1,j,k)= rob*vv(1,j,k)
          rhow_n(1,j,k)= rob*ww(1,j,k)
          rhoe_n(1,j,k)= rob*( ecalc_tro(Tmp(1,j,k),rob) &
                       + 0.5_wp*(uu(1,j,k)**2+vv(1,j,k)**2+ww(1,j,k)**2) )

          ! nullify increments
          ! ------------------
          Krho(1,j,k)=0.0_wp
          Krhou(1,j,k)=0.0_wp
          Krhov(1,j,k)=0.0_wp
          Krhow(1,j,k)=0.0_wp
          Krhoe(1,j,k)=0.0_wp
       enddo
    enddo

    !print *, 'call bc_inflow_imin ok',ntime,irk
 
  end subroutine bc_inflow_imin

  !==============================================================================
  subroutine bc_inflow_imax
  !==============================================================================
    !> Subsonic Inflow Boundary Condition at imax: extrapolate Riemann invariant
    !> to impose direction, total pressure and total temperature
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: j,k
    ! Riemann invariant dp-(roc)_0 dU
    real(wp) :: am0,roc0,coeff
    ! quantities at boundary points
    real(wp) :: rob,Tb,pb
    ! ---------------------------------------------------------------------------

    ! Update inlet points
    ! ===================
    do k=1,nz
       do j=1,ny

          ! Quantities from interior nodes
          ! ==============================
          
          ! normal velocity at interior nodes Uin=Ui.nb
          Uin= uu(nx-1,j,k)*nxn_imax(j,k)+vv(nx-1,j,k)*nyn_imax(j,k)+ww(nx-1,j,k)*nzn_imax(j,k)

          ! correct sign to have inwards velocity
          Uin=sign_imax*Uin
          
          ! product (roc)_0 (extrapolated from interior)
          roc0=rho_n(nx-1,j,k)*c_(nx-1,j,k)

          ! compute outgoing characteristic from interior nodes
          ! dp-(roc)_0 dU=cste <=> pb-(roc)_0 Ubn=pi-(roc)_0 Uin=am0
          am0= prs(nx-1,j,k)-roc0*Uin

          ! coefficient taken flow direction into account
          coeff=(cdir_imax(j,k)/roc0)**2

          ! Solve p_b-roc0*U_b-am0=0 with Newton's method
          ! =============================================

          ! initial values (taken as interior values)
          ! --------------
          Tb =T_tot
          pb =p_tot
          rob=rho_tot

          ! Solve using approximate Newton's method
          ! ---------------------------------------
          ! input/output: initial and updated thermo. variables
          ! input: total enthalpy and entropy to be conserved
          ! input: interior part of Riemann invariant am0
          ! input: coeff taking into account flow direction
          call tcalc_Hstot(Tb,rob,pb,H_tot(j,k),s_tot(j,k),am0,coeff)

          ! update thermo variables
          ! -----------------------
          Tmp(nx,j,k)=Tb
          prs(nx,j,k)=pb
          rho_n(nx,j,k)=rob     

          ! update normal velocity
          ! ----------------------
          Ubn=(pb-am0)/roc0

          ! compute velocity components by projecting modulus
          ! -------------------------------------------------
          uu(nx,j,k)= Ubn*db(j,k,1)*cdir_imax(j,k)
          vv(nx,j,k)= Ubn*db(j,k,2)*cdir_imax(j,k)
          ww(nx,j,k)= Ubn*db(j,k,3)*cdir_imax(j,k)

          ! update conservative variables
          ! -----------------------------
          rhou_n(nx,j,k)= rob*uu(nx,j,k)
          rhov_n(nx,j,k)= rob*vv(nx,j,k)
          rhow_n(nx,j,k)= rob*ww(nx,j,k)
          rhoe_n(nx,j,k)= rob*( ecalc_tro(Tmp(nx,j,k),rob) &
                        + 0.5_wp*(uu(nx,j,k)**2+vv(nx,j,k)**2+ww(nx,j,k)**2) )

          ! nullify increments
          ! ------------------
          Krho(nx,j,k)=0.0_wp
          Krhou(nx,j,k)=0.0_wp
          Krhov(nx,j,k)=0.0_wp
          Krhow(nx,j,k)=0.0_wp
          Krhoe(nx,j,k)=0.0_wp
       enddo
    enddo
 
  end subroutine bc_inflow_imax

  !==============================================================================
  subroutine bc_inflow_jmin
  !==============================================================================
    !> Subsonic Inflow Boundary Condition at jmin: extrapolate Riemann invariant
    !> to impose direction, total pressure and total temperature
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,k
    ! Riemann invariant dp-(roc)_0 dU
    real(wp) :: am0,roc0,coeff
    ! quantities at boundary points
    real(wp) :: rob,Tb,pb
    ! ---------------------------------------------------------------------------

    ! Update inlet points
    ! ===================
    do k=1,nz
       do i=1,nx

          ! Quantities from interior nodes
          ! ==============================
          
          ! normal velocity at interior nodes Uin=Ui.nb
          Uin= uu(i,2,k)*nxn_jmin(i,k)+vv(i,2,k)*nyn_jmin(i,k)+ww(i,2,k)*nzn_jmin(i,k)

          ! correct sign to have inwards velocity
          Uin=sign_jmin*Uin
          
          ! product (roc)_0 (extrapolated from interior)
          roc0=rho_n(i,2,k)*c_(i,2,k)

          ! compute outgoing characteristic from interior nodes
          ! dp-(roc)_0 dU=cste <=> pb-(roc)_0 Ubn=pi-(roc)_0 Uin=am0
          am0= prs(i,2,k)-roc0*Uin

          ! coefficient taken flow direction into account
          coeff=(cdir_jmin(i,k)/roc0)**2

          ! Solve p_b-roc0*U_b-am0=0 with Newton's method
          ! =============================================

          ! initial values (taken as interior values)
          ! --------------
          Tb =T_tot
          pb =p_tot
          rob=rho_tot

          ! Solve using approximate Newton's method
          ! ---------------------------------------
          ! input/output: initial and updated thermo. variables
          ! input: total enthalpy and entropy to be conserved
          ! input: interior part of Riemann invariant am0
          ! input: coeff taking into account flow direction
          call tcalc_Hstot(Tb,rob,pb,H_tot(i,k),s_tot(i,k),am0,coeff)

          ! update thermo variables
          ! -----------------------
          Tmp(i,1,k)=Tb
          prs(i,1,k)=pb
          rho_n(i,1,k)=rob     

          ! update normal velocity
          ! ----------------------
          Ubn=(pb-am0)/roc0

          ! compute velocity components by projecting modulus
          ! -------------------------------------------------
          uu(i,1,k)= Ubn*db(i,k,1)*cdir_jmin(i,k)
          vv(i,1,k)= Ubn*db(i,k,2)*cdir_jmin(i,k)
          ww(i,1,k)= Ubn*db(i,k,3)*cdir_jmin(i,k)

          ! update conservative variables
          ! -----------------------------
          rhou_n(i,1,k)= rob*uu(i,1,k)
          rhov_n(i,1,k)= rob*vv(i,1,k)
          rhow_n(i,1,k)= rob*ww(i,1,k)
          rhoe_n(i,1,k)= rob*( ecalc_tro(Tmp(i,1,k),rob) &
                       + 0.5_wp*(uu(i,1,k)**2+vv(i,1,k)**2+ww(i,1,k)**2) )

          ! nullify increments
          ! ------------------
          Krho(i,1,k)=0.0_wp
          Krhou(i,1,k)=0.0_wp
          Krhov(i,1,k)=0.0_wp
          Krhow(i,1,k)=0.0_wp
          Krhoe(i,1,k)=0.0_wp
       enddo
    enddo
 
  end subroutine bc_inflow_jmin

  !==============================================================================
  subroutine bc_inflow_jmax
  !==============================================================================
    !> Subsonic Inflow Boundary Condition at jmax: extrapolate Riemann invariant
    !> to impose direction, total pressure and total temperature
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,k
    ! Riemann invariant dp-(roc)_0 dU
    real(wp) :: am0,roc0,coeff
    ! quantities at boundary points
    real(wp) :: rob,Tb,pb
    ! ---------------------------------------------------------------------------

     ! Update inlet points
    ! ===================
    do k=1,nz
       do i=1,nx

          ! Quantities from interior nodes
          ! ==============================
          
          ! normal velocity at interior nodes Uin=Ui.nb
          Uin= uu(i,ny-1,k)*nxn_jmax(i,k)+vv(i,ny-1,k)*nyn_jmax(i,k)+ww(i,ny-1,k)*nzn_jmax(i,k)

          ! correct sign to have inwards velocity
          Uin=sign_jmax*Uin
          
          ! product (roc)_0 (extrapolated from interior)
          roc0=rho_n(i,ny-1,k)*c_(i,ny-1,k)

          ! compute outgoing characteristic from interior nodes
          ! dp-(roc)_0 dU=cste <=> pb-(roc)_0 Ubn=pi-(roc)_0 Uin=am0
          am0= prs(i,ny-1,k)-roc0*Uin

          ! coefficient taken flow direction into account
          coeff=(cdir_jmax(i,k)/roc0)**2

          ! Solve p_b-roc0*U_b-am0=0 with Newton's method
          ! =============================================

          ! initial values (taken as interior values)
          ! --------------
          Tb =T_tot
          pb =p_tot
          rob=rho_tot

          ! Solve using approximate Newton's method
          ! ---------------------------------------
          ! input/output: initial and updated thermo. variables
          ! input: total enthalpy and entropy to be conserved
          ! input: interior part of Riemann invariant am0
          ! input: coeff taking into account flow direction
          call tcalc_Hstot(Tb,rob,pb,H_tot(i,k),s_tot(i,k),am0,coeff)

          ! update thermo variables
          ! -----------------------
          Tmp(i,ny,k)=Tb
          prs(i,ny,k)=pb
          rho_n(i,ny,k)=rob     

          ! update normal velocity
          ! ----------------------
          Ubn=(pb-am0)/roc0

          ! compute velocity components by projecting modulus
          ! -------------------------------------------------
          uu(i,ny,k)= Ubn*db(i,k,1)*cdir_jmax(i,k)
          vv(i,ny,k)= Ubn*db(i,k,2)*cdir_jmax(i,k)
          ww(i,ny,k)= Ubn*db(i,k,3)*cdir_jmax(i,k)

          ! update conservative variables
          ! -----------------------------
          rhou_n(i,ny,k)= rob*uu(i,ny,k)
          rhov_n(i,ny,k)= rob*vv(i,ny,k)
          rhow_n(i,ny,k)= rob*ww(i,ny,k)
          rhoe_n(i,ny,k)= rob*( ecalc_tro(Tmp(i,ny,k),rob) &
                        + 0.5_wp*(uu(i,ny,k)**2+vv(i,ny,k)**2+ww(i,ny,k)**2) )

          ! nullify increments
          ! ------------------
          Krho(i,ny,k)=0.0_wp
          Krhou(i,ny,k)=0.0_wp
          Krhov(i,ny,k)=0.0_wp
          Krhow(i,ny,k)=0.0_wp
          Krhoe(i,ny,k)=0.0_wp
       enddo
    enddo

  end subroutine bc_inflow_jmax

  !==============================================================================
  subroutine bc_inflow_kmin
  !==============================================================================
    !> Subsonic Inflow Boundary Condition at kmin: extrapolate Riemann invariant
    !> to impose direction, total pressure and total temperature
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j
    ! Riemann invariant dp-(roc)_0 dU
    real(wp) :: am0,roc0,coeff
    ! quantities at boundary points
    real(wp) :: rob,Tb,pb
    ! ---------------------------------------------------------------------------

    ! Update inlet points
    ! ===================
    do j=1,ny
       do i=1,nx

          ! Quantities from interior nodes
          ! ==============================
          
          ! normal velocity at interior nodes Uin=Ui.nb
          Uin= uu(i,j,2)*nxn_kmin(i,j)+vv(i,j,2)*nyn_kmin(i,j)+ww(i,j,2)*nzn_kmin(i,j)

          ! correct sign to have inwards velocity
          Uin=sign_kmin*Uin
          
          ! product (roc)_0 (extrapolated from interior)
          roc0=rho_n(i,j,2)*c_(i,j,2)

          ! compute outgoing characteristic from interior nodes
          ! dp-(roc)_0 dU=cste <=> pb-(roc)_0 Ubn=pi-(roc)_0 Uin=am0
          am0= prs(i,j,2)-roc0*Uin

          ! coefficient taken flow direction into account
          coeff=(cdir_kmin(i,j)/roc0)**2

          ! Solve p_b-roc0*U_b-am0=0 with Newton's method
          ! =============================================

          ! initial values (taken as interior values)
          ! --------------
          Tb =T_tot
          pb =p_tot
          rob=rho_tot

          ! Solve using approximate Newton's method
          ! ---------------------------------------
          ! input/output: initial and updated thermo. variables
          ! input: total enthalpy and entropy to be conserved
          ! input: interior part of Riemann invariant am0
          ! input: coeff taking into account flow direction
          call tcalc_Hstot(Tb,rob,pb,H_tot(i,j),s_tot(i,j),am0,coeff)

          ! update thermo variables
          ! -----------------------
          Tmp(i,j,1)=Tb
          prs(i,j,1)=pb
          rho_n(i,j,1)=rob     

          ! update normal velocity
          ! ----------------------
          Ubn=(pb-am0)/roc0

          ! compute velocity components by projecting modulus
          ! -------------------------------------------------
          uu(i,j,1)= Ubn*db(i,j,1)*cdir_kmin(i,j)
          vv(i,j,1)= Ubn*db(i,j,2)*cdir_kmin(i,j)
          ww(i,j,1)= Ubn*db(i,j,3)*cdir_kmin(i,j)

          ! update conservative variables
          ! -----------------------------
          rhou_n(i,j,1)= rob*uu(i,j,1)
          rhov_n(i,j,1)= rob*vv(i,j,1)
          rhow_n(i,j,1)= rob*ww(i,j,1)
          rhoe_n(i,j,1)= rob*( ecalc_tro(Tmp(i,j,1),rob) &
                       + 0.5_wp*(uu(i,j,1)**2+vv(i,j,1)**2+ww(i,j,1)**2) )

          ! nullify increments
          ! ------------------
          Krho(i,j,1)=0.0_wp
          Krhou(i,j,1)=0.0_wp
          Krhov(i,j,1)=0.0_wp
          Krhow(i,j,1)=0.0_wp
          Krhoe(i,j,1)=0.0_wp
       enddo
    enddo
 
  end subroutine bc_inflow_kmin

  !==============================================================================
  subroutine bc_inflow_kmax
  !==============================================================================
    !> Subsonic Inflow Boundary Condition at kmax: extrapolate Riemann invariant
    !> to impose direction, total pressure and total temperature
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j
    ! Riemann invariant dp-(roc)_0 dU
    real(wp) :: am0,roc0,coeff
    ! quantities at boundary points
    real(wp) :: rob,Tb,pb
    ! ---------------------------------------------------------------------------

     ! Update inlet points
    ! ===================
    do j=1,ny
       do i=1,nx

          ! Quantities from interior nodes
          ! ==============================
          
          ! normal velocity at interior nodes Uin=Ui.nb
          Uin= uu(i,j,nz-1)*nxn_kmax(i,j)+vv(i,j,nz-1)*nyn_kmax(i,j)+ww(i,j,nz-1)*nzn_kmax(i,j)

          ! correct sign to have inwards velocity
          Uin=sign_kmax*Uin
          
          ! product (roc)_0 (extrapolated from interior)
          roc0=rho_n(i,j,nz-1)*c_(i,j,nz-1)

          ! compute outgoing characteristic from interior nodes
          ! dp-(roc)_0 dU=cste <=> pb-(roc)_0 Ubn=pi-(roc)_0 Uin=am0
          am0= prs(i,j,nz-1)-roc0*Uin

          ! coefficient taken flow direction into account
          coeff=(cdir_kmax(i,j)/roc0)**2

          ! Solve p_b-roc0*U_b-am0=0 with Newton's method
          ! =============================================

          ! initial values (taken as interior values)
          ! --------------
          Tb =T_tot
          pb =p_tot
          rob=rho_tot

          ! Solve using approximate Newton's method
          ! ---------------------------------------
          ! input/output: initial and updated thermo. variables
          ! input: total enthalpy and entropy to be conserved
          ! input: interior part of Riemann invariant am0
          ! input: coeff taking into account flow direction
          call tcalc_Hstot(Tb,rob,pb,H_tot(i,j),s_tot(i,j),am0,coeff)

          ! update thermo variables
          ! -----------------------
          Tmp(i,j,nz)=Tb
          prs(i,j,nz)=pb
          rho_n(i,j,nz)=rob     

          ! update normal velocity
          ! ----------------------
          Ubn=(pb-am0)/roc0

          ! compute velocity components by projecting modulus
          ! -------------------------------------------------
          uu(i,j,nz)= Ubn*db(i,j,1)*cdir_kmax(i,j)
          vv(i,j,nz)= Ubn*db(i,j,2)*cdir_kmax(i,j)
          ww(i,j,nz)= Ubn*db(i,j,3)*cdir_kmax(i,j)

          ! update conservative variables
          ! -----------------------------
          rhou_n(i,j,nz)= rob*uu(i,j,nz)
          rhov_n(i,j,nz)= rob*vv(i,j,nz)
          rhow_n(i,j,nz)= rob*ww(i,j,nz)
          rhoe_n(i,j,nz)= rob*( ecalc_tro(Tmp(i,j,nz),rob) &
                        + 0.5_wp*(uu(i,j,nz)**2+vv(i,j,nz)**2+ww(i,j,nz)**2) )

          ! nullify increments
          ! ------------------
          Krho(i,j,nz)=0.0_wp
          Krhou(i,j,nz)=0.0_wp
          Krhov(i,j,nz)=0.0_wp
          Krhow(i,j,nz)=0.0_wp
          Krhoe(i,j,nz)=0.0_wp
       enddo
    enddo

  end subroutine bc_inflow_kmax

  !==============================================================================
  subroutine bc_backpressure_imin
  !==============================================================================
    !> Outflow Boundary Condition at imin: set back pressure
    !> and extrapolate entropy from interior
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: j,k
    ! normal Mach number
    real(wp) :: Min
    ! quantities at boundary points
    real(wp) :: rob,sb
    ! ---------------------------------------------------------------------------
    ! radial equilibrium correction
    real(wp), dimension(ny,nz) :: dpr
    ! ---------------------------------------------------------------------------
   
    ! Compute Radial Equilibrium Approximation (REA): dp/dr=ro*(vth^2)/r
    ! ===============================================  
    if (is_rea) then
       call bc_rea_ox(dpr)
    else
       dpr=0.0_wp
    endif
   
    ! Update outlet points by specifying a back-pressure
    ! ==================================================
    
    do k=1,nz
       do j=1,ny

          ! velocity normal to the boundary
          ! -------------------------------
          ! normal velocity at interior nodes Uin=Ui.nb
          Uin= uu(2,j,k)*nxn_imin(j,k)+vv(2,j,k)*nyn_imin(j,k)+ww(2,j,k)*nzn_imin(j,k)
          ! normal Mach number at interior nodes
          ! (absolute value to suppress particular sign of the normal direction)
          Min= abs(Uin)/c_(2,j,k)

          ! set back pressure when Min<1
          ! ----------------------------
          ! if Min>=1, pressure is extrapolated from interior nodes
          if (Min.lt.1.0_wp) then
             prs(1,j,k)= p_exit+dpr(j,k)
          else
             prs(1,j,k)= prs(2,j,k)
          endif
       
          ! extrapolate entropy from interior nodes
          ! ---------------------------------------
          sb= scalc_tro(Tmp(2,j,k),rho_n(2,j,k))

          ! update thermo variables
          ! -----------------------
          rob= rocalc_ps(prs(1,j,k),sb,Tmp(2,j,k))
          rho_n(1,j,k)= rob
          Tmp(1,j,k)= tcalc_pro(prs(1,j,k),rob,Tmp(2,j,k))

          ! extrapolate velocity components from interior nodes
          ! ---------------------------------------------------
          uu(1,j,k) = uu(2,j,k)
          vv(1,j,k) = vv(2,j,k)
          ww(1,j,k) = ww(2,j,k)

          ! update conservative variables
          ! -----------------------------
          rhou_n(1,j,k)= rob*uu(1,j,k)
          rhov_n(1,j,k)= rob*vv(1,j,k)
          rhow_n(1,j,k)= rob*ww(1,j,k)
          rhoe_n(1,j,k)= rob*( ecalc_tro(Tmp(1,j,k),rob) &
                        + 0.5_wp*(uu(1,j,k)**2+vv(1,j,k)**2+ww(1,j,k)**2) )

          ! nullify increments
          ! ------------------
          Krho(1,j,k)=0.0_wp
          Krhou(1,j,k)=0.0_wp
          Krhov(1,j,k)=0.0_wp
          Krhow(1,j,k)=0.0_wp
          Krhoe(1,j,k)=0.0_wp
         
       enddo
    enddo
  
  end subroutine bc_backpressure_imin
  
  !==============================================================================
  subroutine bc_backpressure_imax
  !==============================================================================
    !> Outflow Boundary Condition at imax: set back pressure
    !> and extrapolate entropy from interior
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: j,k
    ! normal Mach number
    real(wp) :: Min
    ! quantities at boundary points
    real(wp) :: rob,sb
    ! ---------------------------------------------------------------------------
    ! radial equilibrium correction
    real(wp), dimension(ny,nz) :: dpr
    ! ---------------------------------------------------------------------------

    ! Compute Radial Equilibrium Approximation (REA): dp/dr=ro*(vth^2)/r
    ! ===============================================  
    if (is_rea) then
       call bc_rea_ox(dpr)
       !call bc_rea_analytical(dpr,1)
    else
       dpr=0.0_wp
    endif
   
    ! Update outlet points by specifying a back-pressure
    ! ==================================================   
    do k=1,nz
       do j=1,ny

          ! velocity normal to the boundary
          ! -------------------------------
          ! normal velocity at interior nodes Uin=Ui.nb
          Uin= uu(nx-1,j,k)*nxn_imax(j,k)+vv(nx-1,j,k)*nyn_imax(j,k)+ww(nx-1,j,k)*nzn_imax(j,k)
          ! normal Mach number at interior nodes
          ! (absolute value to suppress particular sign of the normal direction)
          Min= abs(Uin)/c_(nx-1,j,k)

          ! set back pressure when Min<1
          ! ----------------------------
          ! if Min>=1, pressure is extrapolated from interior nodes
          if (Min.lt.1.0_wp) then
             prs(nx,j,k)= p_exit+dpr(j,k)
          else
             prs(nx,j,k)= prs(nx-1,j,k)
          endif
       
          ! extrapolate entropy from interior nodes
          ! ---------------------------------------
          sb= scalc_tro(Tmp(nx-1,j,k),rho_n(nx-1,j,k))

          ! update thermo variables
          ! -----------------------
          rob= rocalc_ps(prs(nx,j,k),sb,Tmp(nx-1,j,k))
          rho_n(nx,j,k)= rob
          Tmp(nx,j,k)= tcalc_pro(prs(nx,j,k),rob,Tmp(nx-1,j,k))

          ! extrapolate velocity components from interior nodes
          ! ---------------------------------------------------
          uu(nx,j,k) = uu(nx-1,j,k)
          vv(nx,j,k) = vv(nx-1,j,k)
          ww(nx,j,k) = ww(nx-1,j,k)

          ! update conservative variables
          ! -----------------------------
          rhou_n(nx,j,k)= rob*uu(nx,j,k)
          rhov_n(nx,j,k)= rob*vv(nx,j,k)
          rhow_n(nx,j,k)= rob*ww(nx,j,k)
          rhoe_n(nx,j,k)= rob*( ecalc_tro(Tmp(nx,j,k),rob) &
                        + 0.5_wp*(uu(nx,j,k)**2+vv(nx,j,k)**2+ww(nx,j,k)**2) )

          ! nullify increments
          ! ------------------
          Krho(nx,j,k)=0.0_wp
          Krhou(nx,j,k)=0.0_wp
          Krhov(nx,j,k)=0.0_wp
          Krhow(nx,j,k)=0.0_wp
          Krhoe(nx,j,k)=0.0_wp
         
       enddo
    enddo
  
  end subroutine bc_backpressure_imax
  
  !==============================================================================
  subroutine bc_backpressure_jmin
  !==============================================================================
    !> Outflow Boundary Condition at jmin: set back pressure
    !> and extrapolate entropy from interior
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,k
    ! normal Mach number
    real(wp) :: Min
    ! quantities at boundary points
    real(wp) :: rob,sb
    ! ---------------------------------------------------------------------------
    ! radial equilibrium correction
    real(wp), dimension(nx,nz) :: dpr
    ! ---------------------------------------------------------------------------
   
    ! Compute Radial Equilibrium Approximation (REA): dp/dr=ro*(vth^2)/r
    ! ===============================================  
    if (is_rea) then
       call bc_rea_ox(dpr)
    else
       dpr=0.0_wp
    endif
   
    ! Update outlet points by specifying a back-pressure
    ! ==================================================
    
    do k=1,nz
       do i=1,nx

          ! velocity normal to the boundary
          ! -------------------------------
          ! normal velocity at interior nodes Uin=Ui.nb
          Uin= uu(i,2,k)*nxn_jmin(i,k)+vv(i,2,k)*nyn_jmin(i,k)+ww(i,2,k)*nzn_jmin(i,k)
          ! normal Mach number at interior nodes
          ! (absolute value to suppress particular sign of the normal direction)
          Min= abs(Uin)/c_(i,2,k)

          ! set back pressure when Min<1
          ! ----------------------------
          ! if Min>=1, pressure is extrapolated from interior nodes
          if (Min.lt.1.0_wp) then
             prs(i,1,k)= p_exit+dpr(i,k)
          else
             prs(i,1,k)= prs(i,2,k)
          endif
       
          ! extrapolate entropy from interior nodes
          ! ---------------------------------------
          sb= scalc_tro(Tmp(i,2,k),rho_n(i,2,k))

          ! update thermo variables
          ! -----------------------
          rob= rocalc_ps(prs(i,1,k),sb,Tmp(i,2,k))
          rho_n(i,1,k)= rob
          Tmp(i,1,k)= tcalc_pro(prs(i,1,k),rob,Tmp(i,2,k))

          ! extrapolate velocity components from interior nodes
          ! ---------------------------------------------------
          uu(i,1,k) = uu(i,2,k)
          vv(i,1,k) = vv(i,2,k)
          ww(i,1,k) = ww(i,2,k)

          ! update conservative variables
          ! -----------------------------
          rhou_n(i,1,k)= rob*uu(i,1,k)
          rhov_n(i,1,k)= rob*vv(i,1,k)
          rhow_n(i,1,k)= rob*ww(i,1,k)
          rhoe_n(i,1,k)= rob*( ecalc_tro(Tmp(i,1,k),rob) &
                        + 0.5_wp*(uu(i,1,k)**2+vv(i,1,k)**2+ww(i,1,k)**2) )

          ! nullify increments
          ! ------------------
          Krho(i,1,k)=0.0_wp
          Krhou(i,1,k)=0.0_wp
          Krhov(i,1,k)=0.0_wp
          Krhow(i,1,k)=0.0_wp
          Krhoe(i,1,k)=0.0_wp
         
       enddo
    enddo
  
  end subroutine bc_backpressure_jmin
  
  !==============================================================================
  subroutine bc_backpressure_jmax
  !==============================================================================
    !> Outflow Boundary Condition at jmax: set back pressure
    !> and extrapolate entropy from interior
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,k
    ! normal Mach number
    real(wp) :: Min
    ! quantities at boundary points
    real(wp) :: rob,sb
    ! ---------------------------------------------------------------------------
    ! radial equilibrium correction
    real(wp), dimension(nx,nz) :: dpr
    ! ---------------------------------------------------------------------------
   
    ! Compute Radial Equilibrium Approximation (REA): dp/dr=ro*(vth^2)/r
    ! ===============================================  
    if (is_rea) then
       call bc_rea_ox(dpr)
    else
       dpr=0.0_wp
    endif
   
    ! Update outlet points by specifying a back-pressure
    ! ==================================================
      
    do k=1,nz
       do i=1,nx

          ! velocity normal to the boundary
          ! -------------------------------
          ! normal velocity at interior nodes Uin=Ui.nb
          Uin= uu(i,ny-1,k)*nxn_jmax(i,k)+vv(i,ny-1,k)*nyn_jmax(i,k)+ww(i,ny-1,k)*nzn_jmax(i,k)
          ! normal Mach number at interior nodes
          ! (absolute value to suppress particular sign of the normal direction)
          Min= abs(Uin)/c_(i,ny-1,k)

          ! set back pressure when Min<1
          ! ----------------------------
          ! if Min>=1, pressure is extrapolated from interior nodes
          if (Min.lt.1.0_wp) then
             prs(i,ny,k)= p_exit+dpr(i,k)
          else
             prs(i,ny,k)= prs(i,ny-1,k)
          endif

!!$          ! ************************************************
!!$          ! used temporary for MAH (blow-up in tcalc_sro_mah)
!!$          ! ************************************************
!!$          
!!$          ! extrapolate temperature from interior nodes
!!$          ! -------------------------------------------
!!$          Tmp(i,ny,k) = Tmp(i,ny-1,k)!2.0_wp*Tmp(i,ny-1,k)-Tmp(i,ny-2,k)
!!$
!!$          ! update thermo variables
!!$          ! -----------------------
!!$          rob= rocalc_pt(prs(i,ny,k),Tmp(i,ny,k),rho_n(i,ny-1,k))
!!$          rho_n(i,ny,k)= rob
!!$          Tmp(i,ny,k)= tcalc_pro(prs(i,ny,k),rob,Tmp(i,ny-1,k))

          ! extrapolate entropy from interior nodes
          ! ---------------------------------------
          sb= scalc_tro(Tmp(i,ny-1,k),rho_n(i,ny-1,k))

          ! update thermo variables
          ! -----------------------
          rob= rocalc_ps(prs(i,ny,k),sb,Tmp(i,ny-1,k))
          rho_n(i,ny,k)= rob
          Tmp(i,ny,k)= tcalc_pro(prs(i,ny,k),rob,Tmp(i,ny-1,k))

          ! extrapolate velocity components from interior nodes
          ! ---------------------------------------------------
          uu(i,ny,k) = uu(i,ny-1,k)
          vv(i,ny,k) = vv(i,ny-1,k)
          ww(i,ny,k) = ww(i,ny-1,k)

          ! update conservative variables
          ! -----------------------------
          rhou_n(i,ny,k)= rob*uu(i,ny,k)
          rhov_n(i,ny,k)= rob*vv(i,ny,k)
          rhow_n(i,ny,k)= rob*ww(i,ny,k)
          rhoe_n(i,ny,k)= rob*( ecalc_tro(Tmp(i,ny,k),rob) &
                        + 0.5_wp*(uu(i,ny,k)**2+vv(i,ny,k)**2+ww(i,ny,k)**2) )

          ! nullify increments
          ! ------------------
          Krho(i,ny,k)=0.0_wp
          Krhou(i,ny,k)=0.0_wp
          Krhov(i,ny,k)=0.0_wp
          Krhow(i,ny,k)=0.0_wp
          Krhoe(i,ny,k)=0.0_wp
         
       enddo
    enddo

  end subroutine bc_backpressure_jmax
  
  !==============================================================================
  subroutine bc_backpressure_kmin
  !==============================================================================
    !> Outflow Boundary Condition at kmin: set back pressure
    !> and extrapolate entropy from interior
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j
    ! normal Mach number
    real(wp) :: Min
    ! quantities at boundary points
    real(wp) :: rob,sb
    ! ---------------------------------------------------------------------------
    ! radial equilibrium correction
    real(wp), dimension(nx,ny) :: dpr
    ! ---------------------------------------------------------------------------
   
    ! Compute Radial Equilibrium Approximation (REA): dp/dr=ro*(vth^2)/r
    ! ===============================================  
    if (is_rea) then
       call bc_rea_ox(dpr)
    else
       dpr=0.0_wp
    endif
   
    ! Update outlet points by specifying a back-pressure
    ! ==================================================
    
    do j=1,ny
       do i=1,nx

          ! velocity normal to the boundary
          ! -------------------------------
          ! normal velocity at interior nodes Uin=Ui.nb
          Uin= uu(i,j,2)*nxn_kmin(i,j)+vv(i,j,2)*nyn_kmin(i,j)+ww(i,j,2)*nzn_kmin(i,j)
          ! normal Mach number at interior nodes
          ! (absolute value to suppress particular sign of the normal direction)
          Min= abs(Uin)/c_(i,j,2)

          ! set back pressure when Min<1
          ! ----------------------------
          ! if Min>=1, pressure is extrapolated from interior nodes
          if (Min.lt.1.0_wp) then
             prs(i,j,1)= p_exit+dpr(i,j)
          else
             prs(i,j,1)= prs(i,j,2)
          endif
       
          ! extrapolate entropy from interior nodes
          ! ---------------------------------------
          sb= scalc_tro(Tmp(i,j,2),rho_n(i,j,2))

          ! update thermo variables
          ! -----------------------
          rob= rocalc_ps(prs(i,j,1),sb,Tmp(i,j,2))
          rho_n(i,j,1)= rob
          Tmp(i,j,1)= tcalc_pro(prs(i,j,1),rob,Tmp(i,j,2))

          ! extrapolate velocity components from interior nodes
          ! ---------------------------------------------------
          uu(i,j,1) = uu(i,j,2)
          vv(i,j,1) = vv(i,j,2)
          ww(i,j,1) = ww(i,j,2)

          ! update conservative variables
          ! -----------------------------
          rhou_n(i,j,1)= rob*uu(i,j,1)
          rhov_n(i,j,1)= rob*vv(i,j,1)
          rhow_n(i,j,1)= rob*ww(i,j,1)
          rhoe_n(i,j,1)= rob*( ecalc_tro(Tmp(i,j,1),rob) &
                        + 0.5_wp*(uu(i,j,1)**2+vv(i,j,1)**2+ww(i,j,1)**2) )

          ! nullify increments
          ! ------------------
          Krho(i,j,1)=0.0_wp
          Krhou(i,j,1)=0.0_wp
          Krhov(i,j,1)=0.0_wp
          Krhow(i,j,1)=0.0_wp
          Krhoe(i,j,1)=0.0_wp
         
       enddo
    enddo
  
  end subroutine bc_backpressure_kmin
  
  !==============================================================================
  subroutine bc_backpressure_kmax
  !==============================================================================
    !> Outflow Boundary Condition at kmax: set back pressure
    !> and extrapolate entropy from interior
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j
    ! normal Mach number
    real(wp) :: Min
    ! quantities at boundary points
    real(wp) :: rob,sb
    ! ---------------------------------------------------------------------------
    ! radial equilibrium correction
    real(wp), dimension(nx,ny) :: dpr
    ! ---------------------------------------------------------------------------
   
    ! Compute Radial Equilibrium Approximation (REA): dp/dr=ro*(vth^2)/r
    ! ===============================================  
    if (is_rea) then
       call bc_rea_ox(dpr)
    else
       dpr=0.0_wp
    endif
   
    ! Update outlet points by specifying a back-pressure
    ! ==================================================
      
    do j=1,ny
       do i=1,nx

          ! velocity normal to the boundary
          ! -------------------------------
          ! normal velocity at interior nodes Uin=Ui.nb
          Uin= uu(i,j,nz-1)*nxn_kmax(i,j)+vv(i,j,nz-1)*nyn_kmax(i,j)+ww(i,j,nz-1)*nzn_kmax(i,j)
          ! normal Mach number at interior nodes
          ! (absolute value to suppress particular sign of the normal direction)
          Min= abs(Uin)/c_(i,j,nz-1)

          ! set back pressure when Min<1
          ! ----------------------------
          ! if Min>=1, pressure is extrapolated from interior nodes
          if (Min.lt.1.0_wp) then
             prs(i,j,nz)= p_exit+dpr(i,j)
          else
             prs(i,j,nz)= prs(i,j,nz-1)
          endif

          ! extrapolate entropy from interior nodes
          ! ---------------------------------------
          sb= scalc_tro(Tmp(i,j,nz-1),rho_n(i,j,nz-1))

          ! update thermo variables
          ! -----------------------
          rob= rocalc_ps(prs(i,j,nz),sb,Tmp(i,j,nz-1))
          rho_n(i,j,nz)= rob
          Tmp(i,j,nz)= tcalc_pro(prs(i,j,nz),rob,Tmp(i,j,nz-1))

          ! extrapolate velocity components from interior nodes
          ! ---------------------------------------------------
          uu(i,j,nz) = uu(i,j,nz-1)
          vv(i,j,nz) = vv(i,j,nz-1)
          ww(i,j,nz) = ww(i,j,nz-1)

          ! update conservative variables
          ! -----------------------------
          rhou_n(i,j,nz)= rob*uu(i,j,nz)
          rhov_n(i,j,nz)= rob*vv(i,j,nz)
          rhow_n(i,j,nz)= rob*ww(i,j,nz)
          rhoe_n(i,j,nz)= rob*( ecalc_tro(Tmp(i,j,nz),rob) &
                        + 0.5_wp*(uu(i,j,nz)**2+vv(i,j,nz)**2+ww(i,j,nz)**2) )

          ! nullify increments
          ! ------------------
          Krho(i,j,nz)=0.0_wp
          Krhou(i,j,nz)=0.0_wp
          Krhov(i,j,nz)=0.0_wp
          Krhow(i,j,nz)=0.0_wp
          Krhoe(i,j,nz)=0.0_wp
         
       enddo
    enddo

  end subroutine bc_backpressure_kmax
  
  !==============================================================================
  subroutine bc_inlet_test_imin
  !==============================================================================
    !> Subsonic Inflow Boundary Condition at imin: extrapolate Riemann invariant
    !> to impose direction, total pressure and total temperature
  !==============================================================================
    use mod_time ! TEMP for printing purpose
    use mod_mpi  ! TEMP for printing purpose
    use warnstop ! TEMP for not defined iversion  
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: j,k

    ! Integrated Riemann invariant Rp=U-2c/(gam-1)
    real(wp) :: Rp
    ! Riemann invariant dp-(roc)_0 du
    real(wp) :: am0,roc0,iroc0,coeff1
    ! quantities at interior points
    real(wp) :: V2,ci,ci2,Uit
    ! quantities at boundary points
    real(wp) :: b,Tb,pb,rob
    ! quadratic sol. Carlson BC
    real(wp) :: atmp,btmp,ctmp,dtmp
   
    real(wp) :: c2t ! cos^2(theta)
    real(wp) :: d_b(3)
    
    real(wp) :: Htmp,dHtmpdT,dpdT
        
    real(wp) :: ro1

    real(wp) :: gamt,gam2t,gsgm1
    real(wp) :: gam2 ! =2/(gam-1)
    real(wp) :: gam3 ! =1+2/(gam-1)=(gam+1)/(gam-1)
    real(wp) :: igam3 !=1/gam3=(gam-1)/(gam+1)
    
    ! Newton's algorithm
    integer :: it
    real(wp) :: f,df,dU
    real(wp) :: err,tol
    ! ---------------------------------------------------------------------------
    integer :: iversion
    ! ---------------------------------------------------------------------------

    ! Set total pressure and temperature TO BE CHANGED (should be from param.ini & setupref)
    ! ==================================
!!    T_tot = 287.00_wp
!!    !p_ref = pcalc_tro(T_ref,rho_ref)
!!    !p_tot = p_ref*(1.0_wp+gam1/2.0_wp*Mach**2)**(gam/gam1)
!!    p_tot=129000.0_wp

    T_tot=T_ref
    p_tot=p_ref
    
    !p_tot = p_ref*(1.0_wp+gam1/2.0_wp*Mach**2)**(gam/gam1)
    !T_tot = T_ref*(1.0_wp+gam1/2.0_wp*Mach**2)
    !print *,p_tot,T_tot,p_ref,T_ref
    !stop
  
    ! Imposed direction of the velocity at the boundary TO BE CHANGED (should be from param.ini & setupref)
    ! =================================================
    !theta_ref = 30.0_wp
    theta_ref = 0.0_wp
    theta_ref = theta_ref*pi/180.0_wp
    d_b = (/cos(theta_ref),sin(theta_ref),0.0_wp/)
    ! cos²(theta)
    c2t=cos(theta_ref)**2

    ! Definitions TO BE CHANGED -> to be factorized
    ! ===========
    ! 2/(gam-1) for integrated Riemann invariant
    gam2=2.0_wp/gam1
    ! 1+2/(gam-1)=(gam+1)/(gam-1)
    gam3=1.0_wp+2.0_wp/gam1
    ! 1/gam3=(gam-1)/(gam+1)
    igam3=1.0_wp/gam3
    ! ratio of gam/(gam-1) ! already defined in mod_fluid
    gsgm1=gam/gam1

    ! Version
    ! =======
    ! iversion=1:  original BC of J.-R. Carlson (NASA TM-2011-217181)
    ! iversion=2: corrected BC of J.-R. Carlson (NASA TM-2011-217181)
    ! iversion=3:  original BC of Blazek (Book)
    ! iversion=4: corrected BC of Blazek (Book)
    ! iversion=5: Vuillot BC (adapted from DEGAS/CANARI)
    ! iversion=6: Congedo BC (adapted from DEGAS/CANARI) for real gas (air version)
    ! iversion=7: Camille BC for real gas (air version)
    
    iversion=2
    
    select case (iversion)
       
    case (1) ! original BC of J.-R. Carlson (NASA TM-2011-217181)
             ! **************************************************
       
       ! Update inlet points
       ! ===================
       do k=1,nz
          do j=1,ny

             ! square of velocity vector norm at interior nodes
             V2=uu(2,j,k)*uu(2,j,k)+vv(2,j,k)*vv(2,j,k)+ww(2,j,k)*ww(2,j,k)

             ! compute total enthalpy from interior nodes
             ! total enthalpy is conserved at boundary due to isentropic assumption
             ! Eq.(36) of Carlson11
             H_to=prs(2,j,k)/rho_n(2,j,k)*(gam/gam1)+0.5_wp*V2

             ! normal velocity at interior nodes Uin=Ui.nb
             ! /!\ outwards normal: change sign first component of normal
             Uin=-uu(2,j,k)*nxn_imin(j,k)+vv(2,j,k)*nyn_imin(j,k)+ww(2,j,k)*nzn_imin(j,k)

             ! compute outgoing Riemann invariant R+ at interior nodes to be extrapolated
             ! Eq.(37) of Carlson11
             Rp=-Uin-gam2*c_(2,j,k)

             ! obtain speed of sound at the boundary by solving quadratic equation Eq(41)
             ! cf Eq.(43) of Carlson11
             atmp = 1.0_wp+2.0_wp/gam1
             btmp = 2.0_wp*Rp
             ctmp = 0.5_wp*gam1*(Rp**2-2.0_wp*H_to)
             dtmp = sqrt(btmp**2-4.0_wp*atmp*ctmp+1.0e-20_wp)
             ! use of max Eq.(44) of Carlson11
             c_(1,j,k)=max(-btmp+dtmp,-btmp-dtmp)/(2.0_wp*atmp)

             ! Update step
             ! ===========

             ! update thermodynamic quantities at the boundary
             ! -----------------------------------------------
             ! temperature from new sound speed c_b and stagnation temperature
             Tmp(1,j,k)=T_tot*(c_(1,j,k)**2/gam1/H_to)
             ! pressure from stagnation pressure & temperature
             prs(1,j,k)=p_tot*(Tmp(1,j,k)/T_tot)**(gam/gam1)
             ! density from EoS
             ro1=prs(1,j,k)/rg/Tmp(1,j,k)

             ! update velocity modulus at the boundary
             ! ---------------------------------------
             ! update velocity at the boundary by extrapolating outgoing invariant
             ! Eq(45) of Carlson11
             !Ubn=gam2*c_(1,j,k)-Rp ! first sign corrected
             Ubn=-gam2*c_(1,j,k)-Rp
             ! change sign because of outward normal definition
             Ubn=-Ubn

             if ((mod(ntotal,1000)==0).and.(irk==nrk).and.(j==25).and.(iproc<=3)) print *,'iproc',iproc,'M1',Ubn/c_(1,j,k)

             ! compute velocity components by projecting modulus
             ! -------------------------------------------------
             uu(1,j,k)= Ubn*d_b(1)
             vv(1,j,k)= Ubn*d_b(2)
             ww(1,j,k)= Ubn*d_b(3)

             ! update conservative variables
             ! -----------------------------
              rho_n(1,j,k)= ro1     
             rhou_n(1,j,k)= ro1*uu(1,j,k)
             rhov_n(1,j,k)= ro1*vv(1,j,k)
             rhow_n(1,j,k)= ro1*ww(1,j,k)
             rhoe_n(1,j,k)= ro1*( ecalc_tro(Tmp(1,j,k),ro1) &
                          + 0.5_wp*(uu(1,j,k)**2+vv(1,j,k)**2+ww(1,j,k)**2) )

             ! nullify increments
             ! ------------------
              Krho(1,j,k)=0.0_wp
             Krhou(1,j,k)=0.0_wp
             Krhov(1,j,k)=0.0_wp
             Krhow(1,j,k)=0.0_wp
             Krhoe(1,j,k)=0.0_wp

          enddo
       enddo

    case (2) ! corrected BC of J.-R. Carlson (NASA TM-2011-217181) [tangential comp. added from interior]
             ! ***************************************************

       ! inwards normal: modified version
       ! ===============
       
       ! Update inlet points
       ! ===================
       do k=1,nz
          do j=1,ny

             ! square of velocity vector norm at interior nodes
             V2=uu(2,j,k)*uu(2,j,k)+vv(2,j,k)*vv(2,j,k)+ww(2,j,k)*ww(2,j,k)

             ! compute total enthalpy from interior nodes
             ! total enthalpy is conserved at boundary due to isentropic assumption
             ! Eq.(36) of Carlson11
             H_to=prs(2,j,k)/rho_n(2,j,k)*(gam/gam1)+0.5_wp*V2

             ! normal velocity at interior nodes Uin=Ui.nb
             Uin=uu(2,j,k)*nxn_imin(j,k)+vv(2,j,k)*nyn_imin(j,k)+ww(2,j,k)*nzn_imin(j,k)

             ! correct sign to have inwards velocity
             Uin=sign_imin*Uin
             
             ! tangent velocity at interior nodes Uit
             Uit=0.!sqrt(V2-Uin**2)

             ! compute outgoing Riemann invariant R+ at interior nodes to be extrapolated
             ! /!\ changed sign wrt Eq.(37) of Carlson11
             Rp=Uin-gam2*c_(2,j,k)

             ! obtain speed of sound at the boundary by solving quadratic equation Eq(41)
             ! version with tangential velocity component
             ! simplified b'=2b and gam3=1.0_wp+2.0_wp/gam1
             ! select solution (-b'+sqrt(delta))/a
             ! different from Carlson11 (use of max in Eq.(44))
             c_(1,j,k)=(-Rp+sqrt(Rp**2-0.5_wp*gam3*gam1*(Rp**2+Uit**2-2.0_wp*H_to)))*igam3

             ! Update step
             ! ===========
             
             ! update thermodynamic quantities at the boundary
             ! -----------------------------------------------
             ! temperature from new sound speed c_b and stagnation temperature
             Tmp(1,j,k)=T_tot*(c_(1,j,k)**2/gam1/H_to)
             ! pressure from stagnation pressure & temperature
             prs(1,j,k)=p_tot*(Tmp(1,j,k)/T_tot)**(gam/gam1)
             ! density from EoS
             ro1=prs(1,j,k)/rg/Tmp(1,j,k)

             ! update velocity modulus at the boundary
             ! ---------------------------------------
             ! update velocity at the boundary by extrapolating outgoing invariant
             ! Eq(45) of Carlson11 /!\ sign of gam2*c changed
             Ubn=-gam2*c_(1,j,k)-Rp
             ! add tangential component
             Ubn=sqrt(Ubn**2+Uit**2)

             if ((mod(ntotal,1000)==0).and.(irk==nrk).and.(j==25).and.(iproc<=3)) print *,'iproc',iproc,'M1',Ubn/c_(1,j,k)

             ! compute velocity components by projecting modulus
             ! -------------------------------------------------
             uu(1,j,k)= Ubn*d_b(1)
             vv(1,j,k)= Ubn*d_b(2)
             ww(1,j,k)= Ubn*d_b(3)

!!$             if ((j==ny/2).and.(k==nz/2).and.(irk==nrk)) then
!!$                if (iproc==0) then
!!$                   write(200,'(i8,1x,f15.10,1x,f17.5,1x,f15.10,1x)') ntime,uu(1,j,k),prs(1,j,k),Tmp(1,j,k)
!!$                elseif (iproc==1) then
!!$                   write(201,'(i8,1x,f15.10,1x,f17.5,1x,f15.10,1x)') ntime,uu(1,j,k),prs(1,j,k),Tmp(1,j,k)
!!$                endif
!!$             endif
             
             ! update conservative variables
             ! -----------------------------
              rho_n(1,j,k)= ro1     
             rhou_n(1,j,k)= ro1*uu(1,j,k)
             rhov_n(1,j,k)= ro1*vv(1,j,k)
             rhow_n(1,j,k)= ro1*ww(1,j,k)
             rhoe_n(1,j,k)= ro1*( ecalc_tro(Tmp(1,j,k),ro1) &
                          + 0.5_wp*(uu(1,j,k)**2+vv(1,j,k)**2+ww(1,j,k)**2) )

             ! nullify increments
             ! ------------------
              Krho(1,j,k)=0.0_wp
             Krhou(1,j,k)=0.0_wp
             Krhov(1,j,k)=0.0_wp
             Krhow(1,j,k)=0.0_wp
             Krhoe(1,j,k)=0.0_wp

          enddo
       enddo
       
    case (3) ! original BC of Blazek (Book)
             ! *****************************

       ! Update inlet points
       ! ===================
       do k=1,nz
          !do j=ndy_c,nfy_c
          do j=1,ny

             ! square of velocity vector norm at interior nodes
             V2=uu(2,j,k)*uu(2,j,k)+vv(2,j,k)*vv(2,j,k)+ww(2,j,k)*ww(2,j,k)

             ! compute total enthalpy from the interior nodes
             ! total enthalpy is conserved at boundary due to adiabatic assumption
             H_to = prs(2,j,k)/rho_n(2,j,k)*(gam/gam1)+0.5_wp*V2

             ! normal velocity at interior nodes Uin=Ui.nb
             ! /!\ outwards normal: change sign first component of normal
             Uin = -uu(2,j,k)*nxn_imin(j,k)+vv(2,j,k)*nyn_imin(j,k)+ww(2,j,k)*nzn_imin(j,k)

             ! compute outgoing Riemann invariant R+ at interior nodes to be extrapolated
             ! /!\ changed sign wrt Eq.() of Blazek
             Rp = Uin-gam2*c_(2,j,k)

             ! obtain speed of sound at the boundary
             ! Eq() of Blazek
             atmp=c2t+2.0_wp/gam1
             c_(1,j,k)=-Rp/atmp*(1+cos(theta_ref)*sqrt(gam1*(atmp*H_to/Rp**2-0.5_wp)))

             ! Update step
             ! ===========

             ! update thermodynamic quantities at the boundary
             ! -----------------------------------------------
             ! temperature from new sound speed c_b and stagnation temperature
             Tmp(1,j,k)=T_tot*(c_(1,j,k)**2/gam1/H_to)
             ! pressure from stagnation pressure & temperature
             prs(1,j,k)=p_tot*(Tmp(1,j,k)/T_tot)**(gam/gam1)
             ! density from EoS
             ro1=prs(1,j,k)/rg/Tmp(1,j,k)

             ! update velocity modulus at the boundary
             ! ---------------------------------------
             ! update velocity at the boundary by extrapolating outgoing invariant
             ! /!\ sign of gam2*c changed wrt Blazek
             Ubn=gam2*c_(1,j,k)+Rp
             ! correct modulus with flow angle
             Ubn=Ubn/cos(theta_ref)
             ! change sign because of outward normal definition
             Ubn=-Ubn

             if ((mod(ntotal,1000)==0).and.(irk==nrk).and.(j==25).and.(iproc<=3)) print *,'iproc',iproc,'M1',Ubn/c_(1,j,k)

             ! compute velocity components by projecting modulus
             ! -------------------------------------------------
             uu(1,j,k)= Ubn*d_b(1)
             vv(1,j,k)= Ubn*d_b(2)
             ww(1,j,k)= Ubn*d_b(3)

             ! update conservative variables
             ! -----------------------------
              rho_n(1,j,k)= ro1     
             rhou_n(1,j,k)= ro1*uu(1,j,k)
             rhov_n(1,j,k)= ro1*vv(1,j,k)
             rhow_n(1,j,k)= ro1*ww(1,j,k)
             rhoe_n(1,j,k)= ro1*( ecalc_tro(Tmp(1,j,k),ro1) &
                          + 0.5_wp*(uu(1,j,k)**2+vv(1,j,k)**2+ww(1,j,k)**2) )

             ! nullify increments
             ! ------------------
              Krho(1,j,k)=0.0_wp
             Krhou(1,j,k)=0.0_wp
             Krhov(1,j,k)=0.0_wp
             Krhow(1,j,k)=0.0_wp
             Krhoe(1,j,k)=0.0_wp

          enddo
       enddo
           
    case (4) ! corrected BC of Blazek (Book)
             ! ******************************

       ! Update inlet points
       ! ===================
       do k=1,nz
          !do j=ndy_c,nfy_c
          do j=1,ny

             ! square of velocity vector norm at interior nodes
             V2=uu(2,j,k)*uu(2,j,k)+vv(2,j,k)*vv(2,j,k)+ww(2,j,k)*ww(2,j,k)

             ! compute total enthalpy from the interior nodes
             ! total enthalpy is conserved at boundary due to adiabatic assumption
             H_to = prs(2,j,k)/rho_n(2,j,k)*(gam/gam1)+0.5_wp*V2

             ! normal velocity at interior nodes Uin=Ui.nb
             Uin = uu(2,j,k)*nxn_imin(j,k)+vv(2,j,k)*nyn_imin(j,k)+ww(2,j,k)*nzn_imin(j,k)

             ! correct sign to have inwards velocity
             Uin=sign_imin*Uin
             
             ! compute outgoing Riemann invariant R+ at interior nodes to be extrapolated
             ! /!\ changed sign wrt Eq.() of Blazek
             Rp = Uin-gam2*c_(2,j,k)

             ! obtain speed of sound at the boundary
             ! Eq() of Blazek
             atmp=c2t+2.0_wp/gam1
             c_(1,j,k)=-Rp/atmp*(1+cos(theta_ref)*sqrt(gam1*(atmp*H_to/Rp**2-0.5_wp)))

             ! Update step
             ! ===========

             ! update thermodynamic quantities at the boundary
             ! -----------------------------------------------
             ! temperature from new sound speed c_b and stagnation temperature
             Tmp(1,j,k)=T_tot*(c_(1,j,k)**2/gam1/H_to)
             ! pressure from stagnation pressure & temperature
             prs(1,j,k)=p_tot*(Tmp(1,j,k)/T_tot)**(gam/gam1)
             ! density from EoS
             ro1=prs(1,j,k)/rg/Tmp(1,j,k)

             ! update velocity modulus at the boundary
             ! ---------------------------------------
             ! update velocity at the boundary by extrapolating outgoing invariant
             ! /!\ sign of gam2*c changed wrt Blazek
             Ubn=gam2*c_(1,j,k)+Rp
             ! correct modulus with flow angle
             Ubn=Ubn/cos(theta_ref)

             if ((mod(ntotal,1000)==0).and.(irk==nrk).and.(j==25).and.(iproc<=3)) print *,'iproc',iproc,'M1',Ubn/c_(1,j,k)

             ! compute velocity components by projecting modulus
             ! -------------------------------------------------
             uu(1,j,k)= Ubn*d_b(1)
             vv(1,j,k)= Ubn*d_b(2)
             ww(1,j,k)= Ubn*d_b(3)

             ! update conservative variables
             ! -----------------------------
              rho_n(1,j,k)= ro1     
             rhou_n(1,j,k)= ro1*uu(1,j,k)
             rhov_n(1,j,k)= ro1*vv(1,j,k)
             rhow_n(1,j,k)= ro1*ww(1,j,k)
             rhoe_n(1,j,k)= ro1*( ecalc_tro(Tmp(1,j,k),ro1) &
                          + 0.5_wp*(uu(1,j,k)**2+vv(1,j,k)**2+ww(1,j,k)**2) )

             ! nullify increments
             ! ------------------
              Krho(1,j,k)=0.0_wp
             Krhou(1,j,k)=0.0_wp
             Krhov(1,j,k)=0.0_wp
             Krhow(1,j,k)=0.0_wp
             Krhoe(1,j,k)=0.0_wp

          enddo
       enddo
    
    case (5) ! Vuillot BC (adapted from DEGAS/CANARI)
       ! **************************************
       
       ! Update inlet points
       ! ===================
       do k=1,nz
          !do j=ndy_c,nfy_c
          do j=1,ny

             ! square of velocity vector norm at interior nodes
             V2=uu(2,j,k)*uu(2,j,k)+vv(2,j,k)*vv(2,j,k)+ww(2,j,k)*ww(2,j,k)

             ! compute total enthalpy from the interior nodes
             ! total enthalpy is conserved at boundary due to adiabatic assumption
             H_to = prs(2,j,k)/rho_n(2,j,k)*(gam/gam1)+0.5_wp*V2

             ! normal velocity at interior nodes Uin=Ui.nb
             Uin = uu(2,j,k)*nxn_imin(j,k)+vv(2,j,k)*nyn_imin(j,k)+ww(2,j,k)*nzn_imin(j,k)

             ! correct sign to have inwards velocity
             Uin=sign_imin*Uin
             
             ! interior sound speed & roc0 (extrapolated from interior)
             ci2=gam*rg*Tmp(2,j,k)
             roc0=rho_n(2,j,k)*sqrt(ci2)

             ! compute outgoing characteristic R+ at interior nodes to be extrapolated
             ! dp-(roc)_0 du=cste
             am0= prs(2,j,k)-roc0*Uin

             ! Solve p_b-roc0*U_b-am0=0
             ! ========================
             Ubn=Uin
             err=1.e5_wp
             tol=1.e-6_wp
             it=0
             gamt=gam/ci2
             gam2t=0.5_wp*gam1/ci2

             ! Newton
             do while (err>tol.and.it<10)
                b=1.0_wp-gam2t*Ubn**2
                Tb=T_ref*b
                pb=p_ref*b**gsgm1
                rob=gamt*pb/b

                f=pb-roc0*Ubn-am0
                df=-roc0-rob*Ubn
                dU=-f/df
                it=it+1
                !print *,'it',it,'dU',dU,j,k
                err=abs(dU)
                Ubn=Ubn+dU
             enddo

             ! Update step
             ! ===========

             ! update thermodynamic quantities at the boundary
             ! -----------------------------------------------
             b=1.0_wp-gam2t*Ubn**2
             ! temperature from new sound speed c_b and stagnation temperature
             Tmp(1,j,k)=T_tot*b
             ! pressure from stagnation pressure & temperature
             prs(1,j,k)=p_tot*b**gsgm1
             ! density from EoS
             ro1=prs(1,j,k)/rg/Tmp(1,j,k)

!!$             ! update velocity modulus at the boundary
!!$             ! ---------------------------------------
!!$             ! correct modulus with flow angle
!!$             Ubn=Ubn/cos(theta_ref)

             if ((mod(ntotal,1000)==0).and.(irk==nrk).and.(j==10).and.(iproc<=3)) print *,'iproc',iproc,'M1',Ubn/c_(1,j,k),it

             ! compute velocity components by projecting modulus
             ! -------------------------------------------------
             uu(1,j,k)= Ubn*d_b(1)/d_b(1)
             vv(1,j,k)= Ubn*d_b(2)/d_b(1)
             ww(1,j,k)= Ubn*d_b(3)/d_b(1)

             ! update conservative variables
             ! -----------------------------
              rho_n(1,j,k)= ro1     
             rhou_n(1,j,k)= ro1*uu(1,j,k)
             rhov_n(1,j,k)= ro1*vv(1,j,k)
             rhow_n(1,j,k)= ro1*ww(1,j,k)
             rhoe_n(1,j,k)= ro1*( ecalc_tro(Tmp(1,j,k),ro1) &
                          + 0.5_wp*(uu(1,j,k)**2+vv(1,j,k)**2+ww(1,j,k)**2) )

             ! nullify increments
             ! ------------------
              Krho(1,j,k)=0.0_wp
             Krhou(1,j,k)=0.0_wp
             Krhov(1,j,k)=0.0_wp
             Krhow(1,j,k)=0.0_wp
             Krhoe(1,j,k)=0.0_wp

          enddo
       enddo
    
    case (6) ! Congedo BC (adapted from DEGAS/CANARI) for real gas (air version)
             ! ***************************************************

       rho_tot= rocalc_pt(p_tot,T_tot,rho_n(1,1,1))
       H_to= ecalc_tro(T_tot,rho_tot) + p_tot/rho_tot
       s_to= scalc_tro(T_tot,rho_tot)

       ! Update inlet points
       ! ===================
       do k=1,nz
          !do j=ndy_c,nfy_c
          do j=1,ny

             ! square of velocity vector norm at interior nodes
             V2=uu(2,j,k)*uu(2,j,k)+vv(2,j,k)*vv(2,j,k)+ww(2,j,k)*ww(2,j,k)

             ! compute total enthalpy from the interior nodes
             ! total enthalpy is conserved at boundary due to adiabatic assumption
             !H_to= prs(2,j,k)/rho_n(2,j,k)*(gam/gam1)+0.5_wp*V2

             ! normal velocity at interior nodes Uin=Ui.nb
             Uin= uu(2,j,k)*nxn_imin(j,k)+vv(2,j,k)*nyn_imin(j,k)+ww(2,j,k)*nzn_imin(j,k)

             ! correct sign to have inwards velocity
             Uin=sign_imin*Uin
             
             ! interior sound speed & roc0 (extrapolated from interior)
             ci2=gam*rg*Tmp(2,j,k)
             roc0=rho_n(2,j,k)*sqrt(ci2)

             ! compute outgoing characteristic R+ at interior nodes to be extrapolated
             ! dp-(roc)_0 du=cste
             am0= prs(2,j,k)-roc0*Uin

             ! direction
             !usdn(m) =1./(d0x(m)*nxn(mn)+d0y(m)*nyn(mn)+d0z(m)*nzn(mn))
             !usdn2(m)=usdn(m)**2
                          
             !coeff1=(usdn(m)/roc0(m))**2
             coeff1=(1.0_wp/roc0)**2
             
             ! Solve p_b-roc0*U_b-am0=0 with Newton's method
             ! ========================
             Tb=tcalc_Hstot_(H_to,s_to,Tmp(2,j,k),rho_n(2,j,k),coeff1,am0)

!!$             Ubn=Uin
!!$             err=1.e5_wp
!!$             tol=1.e-6_wp
!!$             it=0
!!$
!!$             Tb=Tmp(2,j,k)
!!$             rob=rho_n(2,j,k)
!!$             pb=prs(2,j,k)
!!$             
!!$             
!!$             do while ((err>tol).and.(it<50))
!!$
!!$                it=it+1
!!$
!!$                Tb1=Tb
!!$                rob=rocalc_st(s_to,Tb1,rob)
!!$                !cb=c2calc_tro(Tb1,rob)
!!$                
!!$                !roc0=rob*sqrt(cb)
!!$
!!$                ! Riemann invariant
!!$                !am0=pb-roc0*Ubn
!!$
!!$                !coeff1=(usdn(m)/roc0(m))**2
!!$                coeff1=(1.0_wp/roc0)**2
!!$
!!$                Tb1=tcalc_Htot(H_to,Tb,rob,coeff1,am0)
!!$                err=abs((Tb1-Tb)/Tb1)
!!$
!!$                Tb=Tb1
!!$             enddo
             
!!$             rob=rocalc_st(s_to,Tb,rob)
!!$             pb=pcalc_tro(Tb,rob)
!!$
!!$             if ((mod(ntotal,100)==0).and.(irk==nrk).and.(j==10).and.(iproc<=3)) print *,'iproc',iproc,'M1',Ubn/c_(1,j,k)

             ! update thermo variables
             ! -----------------------
             !!Tb=tcalc_Hstot(H_to,s_to,Tmp(2,j,k),rho_n(2,j,k),coeff1,am0)

             ro1=rocalc_st(s_to,Tb,rho_n(2,j,k))
             prs(1,j,k)=pcalc_tro(Tb,ro1)
             Tmp(1,j,k)=Tb
             
             ! update normal velocity
             ! ----------------------
             Ubn=(prs(1,j,k)-am0)/roc0

             if ((mod(ntotal,1000)==0).and.(irk==nrk).and.(j==10).and.(iproc<=3)) print *,'iproc',iproc,'M1',Ubn/c_(1,j,k)

             ! compute velocity components by projecting modulus
             ! -------------------------------------------------
             uu(1,j,k)= Ubn*d_b(1)
             vv(1,j,k)= Ubn*d_b(2)
             ww(1,j,k)= Ubn*d_b(3)

             ! update conservative variables
             ! -----------------------------
              rho_n(1,j,k)= ro1     
             rhou_n(1,j,k)= ro1*uu(1,j,k)
             rhov_n(1,j,k)= ro1*vv(1,j,k)
             rhow_n(1,j,k)= ro1*ww(1,j,k)
             rhoe_n(1,j,k)= ro1*( ecalc_tro(Tmp(1,j,k),ro1) &
                          + 0.5_wp*(uu(1,j,k)**2+vv(1,j,k)**2+ww(1,j,k)**2) )

             ! nullify increments
             ! ------------------
              Krho(1,j,k)=0.0_wp
             Krhou(1,j,k)=0.0_wp
             Krhov(1,j,k)=0.0_wp
             Krhow(1,j,k)=0.0_wp
             Krhoe(1,j,k)=0.0_wp
           enddo
       enddo
        
    case (7) ! Camille BC for real gas (air version)
             ! *************************************
       
       rho_tot = rocalc_pt(p_tot,T_tot,rho_n(1,1,1))
       H_to = ecalc_tro(T_tot,rho_tot) + p_tot/rho_tot
       s_to = scalc_tro(T_tot,rho_tot) ! s0 = sinf
    
       ! Initialize
       tol = 1.e-6_wp
       
       do k=1,nz
          do j=ndy_c,nfy_c
          !do j=1,ny
             
             err=1.0_wp
             it=0
             
             ! use neighbour point to initialize boundary state
             rob = rho_n(2,j,k)
             pb   = prs(2,j,k)
             Tb   = Tmp(2,j,k)

             ! interior nodes
             ci = sqrt(c2calc_tro(Tmp(2,j,k),rho_n(2,j,k)))
             iroc0=1.0_wp/rho_n(2,j,k)/ci

             ! normal velocity at interior nodes Uin=Ui.nb
             Uin = uu(2,j,k)*nxn_imin(j,k)+vv(2,j,k)*nyn_imin(j,k)+ww(2,j,k)*nzn_imin(j,k)

             ! correct sign to have inwards velocity
             Uin=sign_imin*Uin
             
             ! Newton iteration: rho(n+1) = rho(n) - (Htmp(rho(n)-H_to))/dHtmp/drho
             ! Htmp(T,rho) = e(T,rho) + 0.5u**2(rho) + P(T,rho)/rho

             ! Newton iteration: T(n+1) = T(n) - (Htmp(T(n)-H_to))/dHtmp/dT
             ! Htmp(T,rho) = e(T,rho) + 0.5u**2(rho) + P(T,rho)/rho
             
             !do while (err1.gt.tol.or.err.gt.tol)
             do while (err.gt.tol)

                ! d(P/rho)/dT = dP/dT/rho
                dpdT = dpdicalc_tro(Tb,rob)*cvcalc_tro(Tb,rob)

                ! extrapolate R+ from the interior nodes, so dR=0
                ! dR = u_i-u_bc + 1/(rho_i*c_i)(P_i-P_bc)
                Ubn = Uin + iroc0*(pb-prs(2,j,k)) ! /!\ changed sign wrt to Camille (direction normal)

                ! total enthalpy derivative w.r.t temperature
                !dHtmpdT = dedTcalc_tro(Tb,rob) + 0.5_wp*dveldT + dpdT/rob
                dHtmpdT =  cvfg + dpdT/rob + Ubn*iroc0*dpdT 

                ! total enthalpy
                Htmp = ecalc_tro(Tb,rob)+pb/rob + 0.5_wp*Ubn**2 

                ! new thermos value
                df=-(Htmp-H_to)/dHtmpdT
                Tb = Tb+df
                rob= rocalc_st(s_to,Tb,rob)
                pb = pcalc_tro(Tb,rob)

                ! check convergence 
                err=abs(df)/Tb
                it=it+1
             enddo

             if ((mod(ntotal,1000)==0).and.(irk==nrk).and.(j==10).and.(iproc<=3)) print *,'iproc',iproc,'M1',Ubn/c_(1,j,k),it

             ! update thermo variables
             ! -----------------------
             ro1=rob
             prs(1,j,k)=pb 
             Tmp(1,j,k)=Tb
             
             ! compute velocity components by projecting modulus
             ! -------------------------------------------------
             uu(1,j,k)= Ubn*d_b(1)
             vv(1,j,k)= Ubn*d_b(2)
             ww(1,j,k)= Ubn*d_b(3)

             ! update conservative variables
             ! -----------------------------
              rho_n(1,j,k)= ro1     
             rhou_n(1,j,k)= ro1*uu(1,j,k)
             rhov_n(1,j,k)= ro1*vv(1,j,k)
             rhow_n(1,j,k)= ro1*ww(1,j,k)
             rhoe_n(1,j,k)= ro1*( ecalc_tro(Tmp(1,j,k),ro1) &
                          + 0.5_wp*(uu(1,j,k)**2+vv(1,j,k)**2+ww(1,j,k)**2) )

             ! nullify increments
             ! ------------------
              Krho(1,j,k)=0.0_wp
             Krhou(1,j,k)=0.0_wp
             Krhov(1,j,k)=0.0_wp
             Krhow(1,j,k)=0.0_wp
             Krhoe(1,j,k)=0.0_wp
          enddo
       enddo
       
    end select
 
  end subroutine bc_inlet_test_imin

end module mod_bc_inlet_outlet
