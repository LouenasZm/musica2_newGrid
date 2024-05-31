!===============================================================================
subroutine setupref
!===============================================================================
  !> Compute reference quantities depending on the parameters and the case
  !> TO BE CHANGED
!===============================================================================
   use mod_mpi
   use mod_constant
   use mod_eos
   use mod_tranprop

   use mod_forcing_bulk

   use mod_flow
   use mod_io
   use mod_time

   use mod_routines ! <- for src param (source) TO BE CHANGED

   use mod_init_hit ! <- for ak0 TO BE CHANGED
   implicit none
   ! ---------------------------------------------------------------------------
   ! ! for CHIT
   ! real(wp) :: Re_lambda,L_taylor,Atau,tau_eddy
   ! real(wp) :: vn2_init,dudx2_init,vort2_init
   ! ---------------------------------------------------------------------------

!!$   print *,'T_ref',T_ref,'rho_ref',rho_ref,'rho_c',roc
!!$   p_ref= 904388.0_wp
!!$   rho_ref=rocalc_pt(p_ref,T_ref,roc)
!!$   print *,'p_ref',p_ref,'rho_ref',rho_ref
!!$   stop

   ! Compute reference thermodynamic quantities (default)
   ! ==========================================
   
   ! pressure from temperature and density (in param.ini)
   ! -------------------------------------
   p_ref= pcalc_tro(T_ref,rho_ref)

   !p_ref=2.3e5_wp
   !rho_ref=rocalc_pt(p_ref,T_ref,rho_ref)
   !if (iproc==0) print *,p_ref,T_ref,rho_ref
   !call mpistop('prop',0)

   ! internal energy
   ! ---------------
   e_ref= ecalc_tro(T_ref,rho_ref)
   ! sound speed
   ! -----------
   c_ref= sqrt(c2calc_tro(T_ref,rho_ref))
   ! dynamic viscosity from thermodynamic conditions
   ! -----------------------------------------------
   mu_ref= viscosity_law(T_ref,rho_ref)
   ! ~> stored in muw_ref before eventual rescaling
   muw_ref= mu_ref
   ! thermal conductivity
   ! --------------------
   la_ref= thconductivity(mu_ref,T_ref,rho_ref) ! <- NOT USED hereafter
   ! fundamental derivative of gas dynamics
   ! --------------------------------------
   g_ref= gcalc_tro(T_ref,rho_ref)  
   ! default reference velocity based on Mach number (in param.ini)
   ! -----------------------------------------------
   u_ref= Mach*c_ref
   ! default reference length (if not specified) is deduced from ref Reynolds number
   ! -------------------------------------------------------------------------------
   if (L_ref==0.0_wp) then
      L_ref=Re_ref*mu_ref/(rho_ref*u_ref)
   else
      Re_ref=L_ref*rho_ref*u_ref/mu_ref
   endif

   if (CHAN) then
      ! Predefined case #3: turbulent channel flow
      ! ==========================================
      
      ! ~> isothermal wall: T_ref is T_wall [in mod_bc]
      !    --------------------------------
      T_ref=T_wall
      
      ! ~> rho_ref given in param is rho_bulk
      !    ----------------------------------
      ! rho_ref not equal to rho_wall [in mod_constant]
      ! but for low Mach (quasi-incompressible)
      rho_wall=rho_ref
      
      ! ~> c_ref is sound speed at the wall
      !    --------------------------------
      c_ref=sqrt(c2calc_tro(T_wall,rho_wall))

      ! ref viscosity is wall viscosity
      muw_ref=viscosity_law(T_wall,rho_wall)
      mu_ref = muw_ref

      ! ref Mach number is based on U_bulk and c_wall
      ! ~> ref velocity is U_bulk
      u_ref=Mach*c_ref
      
      if (.not.is_2D) then
         ! Att: for incompressible channel flow, the Mach
         ! number is based on centerline velocity
         ! values of Reynolds centerline
         Rec=3300.0_wp ! Re180
         ! Rec=7178.0_wp ! Re360
         ! Rec=15494.4_wp ! Re720
         ! correct u_ref
         u_ref=u_ref*Re_ref/Rec
      endif

      ! ref length is deduced from Re_ref
      L_ref=Re_ref*mu_ref/(rho_ref*u_ref)
      ! ~> channel half-width for plane channel flow
      hc=L_ref
      ! bulk quantities
      U_bulk=u_ref
      rho_bulk=rho_ref

      if (.not.is_2D) then
         ! for incompressible channel flow the Mach number
         ! is based on centerline velocity
         ! values of Reynolds centerline
         Re_tau=180.0_wp
         ! Re_tau=360.0_wp
         ! Re_tau=720.0_wp
         ! used to have a first guess of u_tau
         utau = Re_tau*mu_ref/(rho_wall*L_ref)
         ! non-dimensional time scale         
         tscale = hc/utau
      else
         ! for incompressible channel flow the Mach number
         ! is based on centerline velocity
         ! values of Reynolds centerline
         Re_tau=180.0_wp
         ! used to have a first guess of u_tau
         utau = Re_tau*mu_ref/(rho_wall*L_ref)
         ! non-dimensional time scale
         tscale = hc/utau

!         utau=0.0_wp
!         tscale = hc/u_ref
      endif

      if (iproc==0) then
         write(*,*) repeat('=',40)
         write(*,*) 'Ubulk    [m/s]:', U_bulk
         write(*,*) 'utau     [m/s]:', Utau
         write(*,*) 'Rbulk [kg/m^3]:', rho_bulk
         write(*,*) 'Rwall [kg/m^3]:', rho_wall
         write(*,*) 'cwall    [m/s]:', c_ref
         write(*,*) 'Hcanal     [m]:', hc
         write(*,*) 'Tscale     [s]:', tscale
         write(*,*) repeat('=',40)
      endif
      
      ! Streamwise dimension: Lx (times pi*hc)
      longx=4.0_wp ! Re180 & MFU & Re360
      ! longx=8.0_wp ! Re720
      !longx=1.0_wp
      !longx=2.0_wp
      ! Wall scaled deltay: dy0
      dy0p=0.8_wp ! Re180
      ! dy0p=0.95_wp ! Re360
      dy0p=deltay ! Wall-model
      ! dy0p=0.49_wp ! MFU
      !dy0p=1.8_wp
      !dy0p=3.0_wp
      ! Spanwise dimension: Lz (times pi*hc)
      longz=2.0_wp ! Re180 & Re360
      ! longz=4.0_wp ! Re720
      ! longz=0.3333333333333_wp ! MFU
      !longz=1.0_wp

      ! ! To impose forc_rhou
      ! forc_rhou_ref = -0.300153324098E+05   ! Re360, from ref Xavier
      
   elseif (PHILL) then
      ! Predefined case #4: periodic hill flow
      ! =======================================
      
      ! ~> isothermal wall: T_ref is T_wall [in mod_bc]
      T_ref=T_wall
      ! ~> rho_ref given in param is rho_bulk
      ! rho_ref not equal to rho_wall [in mod_constant]
      ! but for low Mach (quasi-incompressible)
      rho_wall=rho_ref
      ! ~> c_ref is sound speed at the wall
      c_ref=sqrt(c2calc_tro(T_wall,rho_wall))

      ! ref viscosity is wall viscosity
      muw_ref = viscosity_law(T_wall,rho_wall)
      mu_ref = muw_ref

      ! ref Mach number is based on U_bulk and c_wall
      ! ~> ref velocity is U_bulk
      u_ref=Mach*c_ref

      ! ref length is deduced from Re_ref
      L_ref=Re_ref*mu_ref/(rho_ref*u_ref)
      ! ~> hill height for periodic hill flow
      hc=L_ref
      U_bulk=u_ref
      rho_bulk=rho_ref
      utau=0.0_wp

      tscale=9.*hc/u_ref

      if (iproc==0) then
         write(*,*) repeat('=',40)
         write(*,*) 'Ubulk    [m/s]:', U_bulk
         write(*,*) 'Rbulk [kg/m^3]:', rho_bulk
         write(*,*) 'Rwall [kg/m^3]:', rho_wall
         write(*,*) 'cwall    [m/s]:', c_ref
         write(*,*) 'H hill     [m]:', hc
         write(*,*) 'Tscale     [s]:', tscale
         write(*,*) repeat('=',40)
      endif

   elseif (CHIT) then
      ! Predefined case #2: Homogeneous Isotropic Turbulence (HIT)
      ! ==========================================================

      ! ref Mach number is turbulent Mach number
      ! ----------------------------------------
      ! ~> Mach = M_t = q/c_ref
      ! with
      ! q^2=<u'^2+v'^2+w'^2>=3 urms^2=2 K (K:kinetic energy)
      ! q=sqrt(3)*urms=sqrt(2*K)
      ! K=3/2*urms^2
      
      ! ref velocity is based on turbulent Mach number
      ! ----------------------------------------------
      ! u_ref=Mach*c_ref
      ! ~> u_ref=q=sqrt(3)*urms=sqrt(2*K)
      
      ! ref length if not specified is deduced from ref Reynolds number
      ! ---------------------------------------------------------------
      if (L_ref==0.0_wp) L_ref=Re_ref*mu_ref/(rho_ref*u_ref)
      

!!$      if (idepart.eq.FROM_SCRATCH) return
!!$
!!$      open(60, file='vn2dudx2.dat', action='read')
!!$      read(60,*) vn2_init
!!$      read(60,*) dudx2_init
!!$      read(60,*) u_ref
!!$      read(60,*) vort2_init
!!$      close(60)
!!$
!!$      Re_lambda = Re_ref
!!$      L_taylor  = sqrt(vn2_init/dudx2_init)
!!$      mu_ref = sqrt(vn2_init/3.0_wp)*L_taylor*rho_ref/Re_lambda
!!$      Atau = 64.0_wp*(0.5_wp*vn2_init)/(3.0_wp*sqrt(2.0_wp*pi)*ak0**5)
!!$      tau_eddy = sqrt(32.0_wp/Atau)*(2.0_wp*pi)**0.25*ak0**(-3.5_wp)
!!$
!!$      tscale = tau_eddy

   elseif (TGV) then
      ! Predefined case #1: Taylor-Green vortex
      ! =======================================

!!$    if (L_ref-two_pi.lt.1e-3) then
!!$       ! if L_ref is 2*pi, we need to rescale viscosity

      tscale = 1.0_wp/u_ref
      mu_ref = rho_ref*u_ref*1.0_wp/Re_ref

   elseif (STBL) then
      ! Predefined case #5: boundary layers
      ! ===================================
      deltas_in = Re_inlet*mu_ref/(rho_ref*c_ref) !!! TO BE CHANGED (Re_inlet is Re_ref redundant)
      
      ! tscale = xg(ngx)/u_ref
      ! tscale = deltax*ngx/u_ref
      tscale = abs(xmax-xmin)/u_ref
      if (is_RANS) tscale=L_ref/u_ref

      utau = 8.0_wp ! Calculated based on database after is is_rfm

      
   elseif ((CAV).and.(iorder_visc.ne.0)) then
      
      tscale =2.0_wp*2.561992378660378e-4_wp/c_ref
      
   elseif (CYL) then
      ! Predefined case #8: flow past a cylinder
      ! ========================================

      if (u_ref==0.0_wp) then
         tscale= L_ref/c_ref
         !u_ref=c_ref
      else
         ! Shedding frequency of roughly St=fD/U=0.2 for cylinders
         tscale= L_ref/u_ref
      endif

      ! omeg_src = (2.0_wp*pi*c_ref)/(omeg_src*deltax)
      ! if (iproc==0) print *,'omeg_src',omeg_src!,u_ref,deltas_forc

   elseif ((TURB).or.(TE)) then
      ! Predefined case #10: flow past a turbine
      ! ========================================

      if (u_ref==0.0_wp) then
         tscale= L_ref/c_ref
         !u_ref=c_ref
      else
         ! Shedding frequency of roughly St=fD/U=0.2 for cylinders
         tscale= L_ref/u_ref
      endif

   elseif (SHIT) then
      ! Predefined case #10: flow past a turbine
      ! ========================================

      if (u_ref==0.0_wp) then
         tscale= L_ref/c_ref
      else
         tscale= L_ref/u_ref
      endif
   
   elseif ((SRC)) then
      !tscale = xg(ngx)/u_ref
      !tscale = xg(ngx)/c_ref
      tscale = 1.0_wp/c_ref
      !tscale =2.0_wp*2.561992378660378e-4_wp/c_ref
      
      ! Reference length
      ! ----------------
      deltas_forc = Re_ref*mu_ref/rho_ref/c_ref
      !!deltas_forc = Re_ref*mu_ref/rho_ref/u_ref
      if (iproc.eq.0) write(*,*) 'deltas_forc:', deltas_forc

      omeg_sb = omeg_sb/deltas_forc*u_ref
      if (is_2d) then
         beta_sb = 0.0_wp
      else
         beta_sb = beta_sb/deltas_forc
      endif
      
      omeg_src = (2.0_wp*pi*c_ref)/(omeg_src*deltax)
      if (iproc==0) print *,'omeg_src',omeg_src!,u_ref,deltas_forc
      
   elseif (ACT) then   
   
      tscale = 1.0_wp/226.0_wp !freq=226Hz in exp !L_ref/u_ref

   elseif (LE) then
      ! Put L_ref = chord = 61.2mm in param.ini
      Re_ref = L_ref*u_ref*rho_ref/mu_ref
      tscale = L_ref/u_ref
      
   endif

    ! Rescaling factor for viscosity
    ! ==============================
   diffscale = mu_ref/muw_ref

   !!!!time_threshold = 0.0_wp*tscale

   ! Grid scaling
   ! ------------
   if (Lgrid.eq.0.0_wp) Lgrid=L_ref
   
   if (iorder_visc==0) then ! Euler equations ~> no Re_ref / L_ref instead
      if (L_ref==0.0_wp) L_ref=1.0_wp
   end if

   ! Print reference quantities at screeen
   ! =====================================
   if (iproc.eq.0) then
      print *,repeat('=',70)
      print *,'Reference quantities'
      print *,repeat('=',70)      
      write(*,*) '    Mach        [-]:', Mach
      write(*,*) '    Re_ref      [-]:', Re_ref
      write(*,*) '    L_ref       [m]:', L_ref
      write(*,*) '    rho_ref [kg/m3]:', rho_ref
      write(*,*) '    T_ref       [K]:', T_ref
      write(*,*) '    p_ref      [Pa]:', p_ref
      write(*,*) '    c_ref     [m/s]:', c_ref
      write(*,*) '    u_ref     [m/s]:', u_ref
      write(*,*) '    mu_ref  [kg/ms]:', mu_ref
      write(*,*) '    muw_ref [kg/ms]:', muw_ref
      write(*,*) '    Gamma_ref   [-]:', g_ref
      write(*,*) '    tscale      [s]:', tscale
      write(*,*) '    diffscale   [-]:', diffscale
      print *,repeat('=',70)
   endif
   
end subroutine setupref
