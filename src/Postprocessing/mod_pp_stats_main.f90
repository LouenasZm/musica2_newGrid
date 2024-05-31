!==============================================================================
module mod_pp_stats_main
!==============================================================================
  !> Module for stats post-processing
  !> Nota: not parallel -> to be run on a single proc
  !>       only for Cartesian coordinate
!==============================================================================
  use mod_pp_stats_var
  use mod_pp_stats_read
  use mod_pp_stats_TKE_budgets
  use mod_pp_stats_skin_friction
  implicit none
  ! ---------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------

contains

  !============================================================================
  subroutine pp_stats_main
  !============================================================================
    !> main routine for stats post-processing
    !> - only for cases CHAN, STBL -
  !============================================================================
    use mod_mpi
    use mod_constant
    implicit none
    ! -------------------------------------------------------------------------
    ! -------------------------------------------------------------------------

    if (iproc==0) then
       print *,repeat('=',70)
       print *,'Post-processing of stats'
       print *,repeat('=',70)
    endif

    ! Post-processing for channel flow
    ! ================================
    if (CHAN) then
       call pp_stats_read_chan

       print *,'  ~> Statistics ..'
       call pp_stats_chan

       print *,'  ~> Reynolds stress budgets ..'
       call pp_TKE_budgets_chan

       print *,'  ~> Mean skin friction decomposition ..'
       !call pp_FIK_chan
    endif

    ! Post-processing for boundary layer
    ! ==================================
    if (STBL) then
       call pp_stats_read_xy

       print *,'  ~> Statistics ..'
       call pp_stats_stbl

       print *,'  ~> Reynolds stress budgets ..'
       call pp_TKE_budgets_xy

       print *,'  ~> Mean skin friction decomposition ..'
       !call pp_FIK_stbl
       call pp_RD_stbl
       call pp_FIK2_stbl
    endif

    if (TURB) then
       ! ! Calculate stagnation & isentropic values for stats
       ! ! --------------------------------------------------
       ! call pp_stats_turb

       ! Calculate mixed state for a line at the outlet
       ! ----------------------------------------------
       call pp_mixed_out_line
    endif

    call mpistop('stop in pp_stats_main', 0)

  end subroutine pp_stats_main

  !============================================================================
  subroutine pp_stats_chan
  !============================================================================
    !> main routine for stats post-processing
    !> - only for cases CHAN, STBL -
  !============================================================================
    use mod_grid
    implicit none
    ! -------------------------------------------------------------------------
    ! -------------------------------------------------------------------------
    
    ! not included yet
    
  end subroutine pp_stats_chan
  
  !============================================================================
  subroutine pp_stats_stbl
  !============================================================================
    !> main routine for stats post-processing
    !> - only for cases CHAN, STBL -
  !============================================================================
    use mod_grid
    use mod_deriv2d
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i,j
    real(wp) :: arg1,arg2,Rem
    ! -------------------------------------------------------------------------
    
    ! Computation of velocity derivatives
    ! ===================================
  
    allocate(duxm(ngx,ngy),dvxm(ngx,ngy),dwxm(ngx,ngy))
    allocate(duym(ngx,ngy),dvym(ngx,ngy),dwym(ngx,ngy))
    allocate(duzm(ngx,ngy),dvzm(ngx,ngy),dwzm(ngx,ngy))

    call deriv2_x_11pts(um,duxm)
    call deriv2_x_11pts(vm,dvxm)
    call deriv2_x_11pts(wm,dwxm)

    call deriv2_y_11pts(um,duym)
    call deriv2_y_11pts(vm,dvym)
    call deriv2_y_11pts(wm,dwym)

    ! homogeneous in z
    duzm=0.0_wp
    dvzm=0.0_wp
    dwzm=0.0_wp

    ! Computation of Favre-averaged quantities
    ! ========================================

    ! mean quantities: [f]=<rhof>/<rho>
    ! ---------------------------------
    ! velocity components [ui] ! <- used in TKE budgets
    allocate(umfavre(ngx,ngy),vmfavre(ngx,ngy),wmfavre(ngx,ngy))
    umfavre= rhoum/rhom
    vmfavre= rhovm/rhom
    wmfavre= rhowm/rhom
    ! temperature [T] ! <- used for what ???????????????????????????
    allocate(Tmfavre(ngx,ngy))
    Tmfavre = rhoTm/rhom

    ! fluctuating quantities: <f''> = <f> - <rhof>/<rho> + <rho'f'>/<rho>
    ! -------------------------------------------------------------------
    ! velocity components <ui''> ! <- used in TKE budgets
    allocate(uffavre(ngx,ngy),vffavre(ngx,ngy),wffavre(ngx,ngy))
    uffavre= um-umfavre !+rhour/rhom
    vffavre= vm-vmfavre !+rhovr/rhom
    wffavre= wm-wmfavre !+rhowr/rhom
    ! temperature <T''> ! <- used for what ???????????????????????????
    allocate(Tffavre(ngx,ngy))
    Tffavre= Tm-Tmfavre !+rhoT1r/rho1m

    ! Determination of wall quantities & friction <- used for some normalization (budgets,..)
    ! ===========================================

    !!!! TO BE CLEANED !!!! What do I need ??
    !real(wp), dimension(:), allocatable :: rhowall,muwall,nuwall,pwall,cwall,ewall,hwall
    !real(wp), dimension(:), allocatable :: Twall,dTwall,lawall,Qwall,tauwall,dudywall

    !! TKE budgets: normb=rhowall*u_tau**4/nuwall
    
    allocate(rhowall(ngx),muwall(ngx),Twall(ngx),dudywall(ngx))
    allocate(pwall(ngx),cwall(ngx),Qwall(ngx))

    rhowall =rhom(:,1)
    muwall  = mum(:,1)
    Twall   =  Tm(:,1)
    !dTwall  =dTmdy(:,1)
    dudywall=duym(:,1)
    pwall   =  pm(:,1)
    cwall   =  cm(:,1)
    !ewall   =  em(:,1)
    !hwall   =  hm(:,1)
    !lawall  = lam(:,1)
    Qwall   =ladTym(:,1)
    !print *,'wall q.',rhowall,muwall,Twall,dudywall,pwall,cwall,tauwall

    ! Friction coefficient & Reynolds number
    ! ======================================
    allocate(tauwall(ngx),nuwall(ngx),u_tau(ngx),c_f(ngx))

    ! wall shear stress
    tauwall=muwall*dudywall
    ! kinematic viscosity
    nuwall =muwall/rhowall
    ! friction velocity
    u_tau   =sqrt(tauwall/rhowall)
    ! friction coefficient
    c_f=tauwall/(0.5_wp*rhowall*u_ref**2)

    ! Determination of thicknesses
    ! ============================
    allocate(j99(ngx),j_inn(ngx),j_log(ngx))
    allocate(delt99(ngx),deltas(ngx),deltatheta(ngx),Hfac(ngx))
    allocate(delt99_(ngx),deltas_(ngx),deltatheta_(ngx),Hfac_(ngx))
    allocate(deltas_c(ngx),deltatheta_c(ngx),Hfac_c(ngx))

    ! 99% thickness
    ! -------------
    do i=1,ngx
       j=1
       do while ((um(i,j)<0.99_wp*u_ref).and.(j<ny))
          j=j+1
       enddo
       delt99(i)=yg(j)
       j99(i)=j
       delt99_(i)=yg(j99(i)-1)+(0.99_wp*u_ref-um(i,j99(i)-1))/(um(i,j99(i))-um(i,j99(i)-1))*(yg(j99(i))-yg(j99(i)-1))
    enddo

    ! index for y^+=30
    ! ----------------
    do i=1,ngx
       j=1
       do while ((yg(j)*u_tau(i)/nuwall(i)<30).and.(j<ny))
          j=j+1
       enddo
       j_inn(i)=j
    enddo

    ! index for 0.3 delta
    ! -------------------
    do i=1,ngx
       j=1
       do while ((yg(j)<0.3_wp*delt99_(i)).and.(j<ny))
          j=j+1
       enddo
       j_log(i)=j
    enddo

    ! Displacement thickness
    ! ----------------------
    ! incompressible version
    deltas=0.0_wp
    do i=1,ngx
       do j=1,j99(i)
          arg1=1.0_wp-um(i,j)/u_ref
          arg2=1.0_wp-um(i,j+1)/u_ref
          deltas(i)= deltas(i)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo
    enddo

    do i=1,ngx
       do j=1,j99(i)-1
          arg1=1.0_wp-um(i,j)/u_ref
          arg2=1.0_wp-um(i,j+1)/u_ref
          deltas_(i)= deltas_(i)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo
       arg1=1.0_wp-um(i,j99(i)-1)/u_ref
       arg2=1.0_wp-0.99_wp
       deltas_(i)= deltas_(i)+0.5_wp*(arg1+arg2)*(delt99_(i)-yg(j99(i)-1))
    enddo
    
    ! compressible version
    deltas_c=0.0_wp
    do i=1,ngx
       do j=1,j99(i)
          arg1=1.0_wp-rhom(i,j)*um(i,j)/rho_ref/u_ref
          arg2=1.0_wp-rhom(i,j+1)*um(i,j+1)/rho_ref/u_ref
          deltas_c(i)= deltas_c(i)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo
    enddo

    ! Momentum thickness
    ! ------------------
    ! incompressible version
    deltatheta=0.0_wp
    do i=1,ngx
       do j=1,j99(i)+10
          arg1=um(i,j)/u_ref*(1.0_wp-um(i,j)/u_ref)
          arg2=um(i,j+1)/u_ref*(1.0_wp-um(i,j+1)/u_ref)
          deltatheta(i)=deltatheta(i)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo
    enddo
    
    do i=1,ngx
       do j=1,j99(i)-1
          arg1=um(i,j)/u_ref*(1.0_wp-um(i,j)/u_ref)
          arg2=um(i,j+1)/u_ref*(1.0_wp-um(i,j+1)/u_ref)
          deltatheta_(i)= deltatheta_(i)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo
       arg1=um(i,j99(i)-1)/u_ref*(1.0_wp-um(i,j99(i)-1)/u_ref)
       arg2=0.99_wp*(1.0_wp-0.99_wp)
       deltatheta_(i)= deltatheta_(i)+0.5_wp*(arg1+arg2)*(delt99_(i)-yg(j99(i)-1))
    enddo
    
    ! compressible version
    deltatheta_c=0.0_wp
    do i=1,ngx
       do j=1,j99(i)+10
          arg1=rhom(i,j)*um(i,j)/rho_ref/u_ref*(1.0_wp-um(i,j)/u_ref)
          arg2=rhom(i,j+1)*um(i,j+1)/rho_ref/u_ref*(1.0_wp-um(i,j+1)/u_ref)
          deltatheta_c(i)=deltatheta_c(i)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo
    enddo

    ! Shape factor
    ! ------------
    ! incompressible version
    Hfac=deltas/deltatheta
    ! incompressible version
    Hfac_=deltas_/deltatheta_
    ! compressible version
    Hfac_c=deltas_c/deltatheta_c

    ! Definition of Reynolds numbers
    ! ==============================
    allocate(Re_d99(ngx),Re_t(ngx),Re_ds(ngx),Re_theta(ngx),Re_d2(ngx),Re_ds_c(ngx),Re_theta_c(ngx))

    Rem=rho_ref*u_ref/mu_ref
    Re_d99=Rem*delt99_
    Re_ds=Rem*deltas_
    Re_theta=Rem*deltatheta_
    Re_t=u_tau*delt99_/nuwall
    ! with wall viscosity
    Re_d2=rho_ref*u_ref*deltatheta/muwall    
    ! compressible version
    Re_ds_c=Rem*deltas_c
    Re_theta_c=Rem*deltatheta_c
  
    open(69,file=trim(dirRESU)//'delt'//trim(name_output)//'.bin',form='unformatted',status='unknown')
    rewind(69)
    write(69) ngx
    write(69) (delt99(i),i=1,ngx)
    write(69) (delt99_(i),i=1,ngx)
    write(69) (deltas(i),i=1,ngx)
    write(69) (deltas_(i),i=1,ngx)
    write(69) (deltatheta(i),i=1,ngx)
    write(69) (deltatheta_(i),i=1,ngx)
    write(69) (Hfac(i),i=1,ngx)
    write(69) (Hfac_(i),i=1,ngx)
    write(69) (Re_theta(i),i=1,ngx)
    Re_theta=Rem*deltatheta_
    write(69) (Re_theta(i),i=1,ngx)
    close(69)
!!$    stop
    
!!$    i=6943
!!$    print *,'Re_theta',Re_theta(i)
    i=7009
    i=7185
    !print *,'Re_theta',Re_theta(i)
    !print *,'Re_theta',Re_theta(5000),Re_theta(6000),Re_theta(7000),Re_theta(8000)
    !stop

   !print *,'Re_theta',Re_theta
!!$    i=5720
!!$    print *,'Re_theta',Re_theta(i)
   
!!$    i=2945
!!$    print *,'Re_theta',Re_theta(i)

!!$  ! Vorticity
!!$  allocate(vrtm(ngy),vrt2m(ngy))
!!$  vrtm=sqrt(vrtxm**2+vrtym**2+vrtzm**2)
!!$  vrt2m=     vrtx2m   +vrty2m   +vrtz2m

  !  delta ??
  !! Re_tau=utau*delta*rhowall/muwall

!!$  ! semi-local
!!$  allocate(lstar(ngy),utaus(ngy),Re_taus(ngy))
!!$  utaus(:) =sqrt(tauwall/rho1m(:))
!!$  Re_taus=rho1m(:)*utaus(:)*hc/mu1m(:)

!!$  c0=sqrt(c2calc_tro(Twall,rhowall))
!!$  cp=cpcalc_tro(Twall,rhowall)

!!$  ! Wall and semi-local scalings
!!$  ! ----------------------------
!!$  tplus=muwall/(rhowall*utau**2)
!!$  lplus=muwall/(rhowall*utau)
!!$  lstar=mu1m/(rho1m*utaus)
!!$  ! wall scaling
!!$  xgp=xg*hc/lplus
!!$  ygp=yg*hc/lplus
!!$  zgp=zg*hc/lplus
!!$  ! semi-local scaling
!!$  ! N.B. -> xgs and zgs are only y-dependent!!
!!$  ygs=yg*hc/lstar
!!$  xgs=(xg(2)-xg(1))*hc/lstar ! /!\ dy non constant
!!$  zgs=(zg(2)-zg(1))*hc/lstar

  end subroutine pp_stats_stbl
  
  !============================================================================
  subroutine pp_stats_turb
  !============================================================================
    !> main routine for stats post-processing for TURB
  !============================================================================
    use mod_grid
    use mod_eos
    use mod_utils
    implicit none
    ! -------------------------------------------------------------------------
    ! total (stagnation) quantities
    integer :: i,j
    real(wp) :: p_tot,T_tot,rho_tot,H_tot,s_tot,g_eq
    real(wp), dimension(:,:), allocatable :: prs_pp,T_is,rho_is,h_is,e_is,M_is,c_is
    real(wp), dimension(:,:), allocatable :: h_0,e_0,prs_0,T_0,rho_0
    ! -------------------------------------------------------------------------

    print *,'  ~> Isentropic evaluation ...'

    ! Determination of isentropic quantities, based on inlet isentropy
    ! ================================================================
    if (is_stagnation) then
       ! total quantities of reference (inlet)
       T_tot=T_ref
       p_tot=p_ref
       rho_tot= rocalc_pt(p_tot,T_tot,rho_ref)
       H_tot= ecalc_tro(T_tot,rho_tot) + p_tot/rho_tot
       s_tot= scalc_tro(T_tot,rho_tot)
       print *,"T_tot",T_tot
       print *,"p_tot",p_tot
       print *,"rho_tot",rho_tot
       print *,"H_tot",H_tot
       print *,"s_tot",s_tot
    else
       call mpistop('is_stagnation must be True, with inlet total quantities in param.ini',0)
    endif

    ! Read stats I/O file
    ! ===================
    call read_write_stats_xy(READ,iblc_pp)

    ! Association for stats_turb
    ! rho <-> avg_t(:,:,1)
    !  p  <-> avg_t(:,:,5)
    !  T  <-> avg_t(:,:,6)
    !  u2 <-> avg_t(:,:,12)
    !  v2 <-> avg_t(:,:,13)
    !  w2 <-> avg_t(:,:,14)
    !  e  <-> avg_t(:,:,24)
    !  h  <-> avg_t(:,:,25)
    !  s  <-> avg_t(:,:,27)
    !  M  <-> avg_t(:,:,28)
    !  Cp <-> avg_t(:,:,32)
    !  Cv <-> avg_t(:,:,33)

    ! Allocation
    ! ==========
    allocate(prs_pp(ngx,ngy))
    allocate(rho_is(ngx,ngy))
    allocate(T_is(ngx,ngy))
    allocate(e_is(ngx,ngy))
    allocate(h_is(ngx,ngy))
    allocate(M_is(ngx,ngy))
    allocate(c_is(ngx,ngy))

    allocate(h_0(ngx,ngy))
    allocate(e_0(ngx,ngy))
    allocate(prs_0(ngx,ngy))
    allocate(T_0(ngx,ngy))
    allocate(rho_0(ngx,ngy))

    ! Calculation of isentropic values
    ! --------------------------------
    do j=1,ngy
       do i=1,ngx
          ! Mean pressure
          prs_pp(i,j) = avg_t(i,j,5)
          ! Isentropic density computed on s_tot and p
          rho_is(i,j) = rocalc_ps(prs_pp(i,j),s_tot,T_tot)
          ! Isentropic temperature computed on s_tot and rho_is
          T_is(i,j) = tcalc_sro(s_tot,rho_is(i,j),T_tot)
          ! Isentropic internal energy based on rho_is and T_is <~ store in rhoe
          e_is(i,j) = ecalc_tro(T_is(i,j),rho_is(i,j))
          ! Isentropic enthalpy based on e_is, p and rho_is <~ store in rhoe_n
          h_is(i,j) = e_is(i,j) + prs_pp(i,j)/rho_is(i,j)
          if (h_is(i,j)>H_tot) then
             print *,"/!\ h_is > H_tot",i,j,(h_is(i,j)-H_tot)/H_tot
             h_is(i,j) = H_tot
          endif
          ! Isentropic sound speed based on T_is and rho_is
          c_is(i,j) = (c2calc_tro(T_is(i,j),rho_is(i,j)))**0.5
          ! Isentropic mach based on c_is and U_is = sqrt(2*(H_tot - h_is))
          M_is(i,j) = (2*(H_tot - h_is(i,j)))**0.5/c_is(i,j)
       enddo
   enddo

    ! Calculation of stagnation values
    ! --------------------------------
    do j=1,ngy
       do i=1,ngx
          ! h0 = h + V**2/2
          h_0(i,j) = avg_t(i,j,25) + 0.5_wp*(avg_t(i,j,12) + avg_t(i,j,13) + avg_t(i,j,14))

          ! Calculation of stagnation quantities
          ! try based on static or pfg relation, sometimes one works and the other don't ~> not robust !
          ! First guess for T_0 & rho_0 based on perfect gas relations
          g_eq = avg_t(i,j,32)/avg_t(i,j,33)
          T_0(i,j) = avg_t(i,j,6)*(1 + (g_eq-1)/2 * avg_t(i,j,28)**2)
          rho_0(i,j) = avg_t(i,j,1)*(1 + (g_eq-1)/2 * avg_t(i,j,28)**2)**(1/(g_eq-1))
          ! ! First guess for T_0 & rho_0 based on static values
          ! T_0(i,j) = avg_t(i,j,6)
          ! rho_0(i,j) = avg_t(i,j,1)
          call stagnation_calc(h_0(i,j),avg_t(i,j,27),e_0(i,j),T_0(i,j),prs_0(i,j),rho_0(i,j))
       enddo
    enddo


    ! Writting isentropic values in file
    ! ==================================
    ! order: rho_is, T_is, e_is and h_is
    open(194,file='isentropic_bl'//trim(numchar(iblc_pp))//'.bin',form='unformatted',status='unknown')
    rewind(194)
    write(194) ((rho_is(i,j),i=1,ngx),j=1,ngy)
    write(194) ((T_is(i,j),i=1,ngx),j=1,ngy)
    write(194) ((e_is(i,j),i=1,ngx),j=1,ngy)
    write(194) ((h_is(i,j),i=1,ngx),j=1,ngy)
    write(194) ((c_is(i,j),i=1,ngx),j=1,ngy)
    write(194) ((M_is(i,j),i=1,ngx),j=1,ngy)
    close(194)

    ! Writting total values in file
    ! =============================
    ! order: rho_is, T_is, e_is and h_is
    open(194,file='stagnation_bl'//trim(numchar(iblc_pp))//'.bin',form='unformatted',status='unknown')
    rewind(194)
    write(194) ((rho_0(i,j),i=1,ngx),j=1,ngy)
    write(194) ((T_0(i,j),i=1,ngx),j=1,ngy)
    write(194) ((e_0(i,j),i=1,ngx),j=1,ngy)
    write(194) ((h_0(i,j),i=1,ngx),j=1,ngy)
    write(194) ((prs_0(i,j),i=1,ngx),j=1,ngy)
    close(194)

  end subroutine pp_stats_turb

  !============================================================================
  subroutine pp_mixed_out_line
  !============================================================================
    !> routine to calculate mixed state on a line at the outlet
  !============================================================================
    use mod_eos
    use mod_utils
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i,n_pts
    real(wp) :: sign_n,normal_x,normal_y,norm_n ,l_line  ! normal to the line
    real(wp) :: m_pt,theta_n,theta_x,theta_y,theta_z,h0_mean   ! quantities conserved
    real(wp) :: Un_bar,u_bar,v_bar,w_bar,rho_bar,p_bar,T_bar,e_bar,h_bar  ! mixed-out state
    real(wp) :: fn,dpdrho,dedrho_p,dhdrho_p,dfdrho,err1,err2  ! Newton loop
    real(wp) :: p_test,rho_test  ! Newton loop
    real(wp) :: M_is,rho_is,T_is,c_is,h_is,e_is ! isentropic mixed-out quantities
    real(wp) :: H_tot,s_tot,p_tot,T_tot,rho_tot ! reference inlet quantities
    real(wp) :: p0_bar,rho0_bar,T0_bar,e0_bar ! stagnation mixed-out quantities
    real(wp) :: c_bar,M_bar,s_bar,g_pv,Gamma,g_eq ! other mixed-out quantities
    real(wp), dimension(:), allocatable :: x_line,y_line,u_line,v_line,w_line,Un_line
    real(wp), dimension(:), allocatable :: rho_line,p_line,T_line,h0_line
    real(wp), dimension(:), allocatable :: dl_line
    ! -------------------------------------------------------------------------
    ! real(wp) :: Q,L,C ! PFG test /!\ Not working...

    ! Reading line of data
    ! --------------------
    ! Should contain:
    ! -> x, y, rho, u, v, w, h0, p
    open(50,file='meas_line.dat',form='formatted')
    rewind(50)
    read(50,*) n_pts
    read(50,*)      ! x y u v w rho p T h0
    allocate(x_line(n_pts)); allocate(y_line(n_pts)); allocate(rho_line(n_pts)); allocate(u_line(n_pts))
    allocate(v_line(n_pts)); allocate(w_line(n_pts)); allocate(h0_line(n_pts)); allocate(p_line(n_pts))
    allocate(T_line(n_pts))
    do i=1,n_pts
        read(50,*)  x_line(i), y_line(i), u_line(i), v_line(i), w_line(i), rho_line(i), p_line(i), T_line(i), h0_line(i)
    enddo
    close(50)

    ! Compute abscisse along the line
    ! -------------------------------
    allocate(dl_line(n_pts-1))
    do i=1,n_pts-1
       dl_line(i) = ((x_line(i+1)-x_line(i))**2 +  (y_line(i+1)-y_line(i))**2)**0.5_wp
    enddo

    l_line = dl_line(1)
    do i=2,n_pts-1
       l_line = l_line + dl_line(i)
    enddo

    ! Compute normal to line
    ! ----------------------
    if (y_line(2)-y_line(1).lt.0.0_wp) then
       sign_n = 1.0_wp
    else
       sign_n = -1.0_wp
    endif
    normal_x = - (y_line(n_pts)-y_line(1))*sign_n
    normal_y = (x_line(n_pts)-x_line(1))*sign_n
    ! Unitary vector
    norm_n = (normal_x**2 + normal_y**2)**0.5
    normal_x = normal_x/norm_n; normal_y = normal_y/norm_n

    ! Compute velocity normal to line
    ! -------------------------------
    allocate(Un_line(n_pts))
    do i=1,n_pts
       Un_line(i) = u_line(i)*normal_x + v_line(i)*normal_y
    enddo

    ! Compute mass, momentum and enthalpy integrals
    ! ---------------------------------------------
    m_pt=0.0_wp; theta_n=0.0_wp;  theta_x=0.0_wp; theta_y=0.0_wp; theta_z=0.0_wp; h0_mean=0.0_wp
    do i=1,n_pts-1
       m_pt = m_pt + (rho_line(i+1)*Un_line(i+1) + rho_line(i)*Un_line(i))*0.5_wp*dl_line(i)
       theta_n = theta_n + (rho_line(i+1)*Un_line(i+1)*Un_line(i+1) + rho_line(i)*Un_line(i)*Un_line(i))*0.5_wp*dl_line(i) + &
                           (p_line(i+1) + p_line(i))*0.5_wp*dl_line(i)
       theta_x = theta_x + (rho_line(i+1)*Un_line(i+1)*u_line(i+1) + rho_line(i)*Un_line(i)*u_line(i))*0.5_wp*dl_line(i) + &
                           (p_line(i+1) + p_line(i))*0.5_wp*normal_x*dl_line(i)
       theta_y = theta_y + (rho_line(i+1)*Un_line(i+1)*v_line(i+1) + rho_line(i)*Un_line(i)*v_line(i))*0.5_wp*dl_line(i) + &
                           (p_line(i+1) + p_line(i))*0.5_wp*normal_y*dl_line(i)
       theta_z = theta_z + (rho_line(i+1)*Un_line(i+1)*w_line(i+1) + rho_line(i)*Un_line(i)*w_line(i))*0.5_wp*dl_line(i)
       h0_mean = h0_mean + (rho_line(i+1)*Un_line(i+1)*h0_line(i+1) + rho_line(i)*Un_line(i)*h0_line(i))*0.5_wp*dl_line(i)
    enddo
    ! integrals normalized by line length
    m_pt=m_pt/l_line; theta_n=theta_n/l_line; theta_x=theta_x/l_line; theta_y=theta_y/l_line;  theta_z=theta_z/l_line; h0_mean=h0_mean/l_line
    ! momentum averaged for h
    h0_mean = h0_mean/m_pt

    ! ---------------------------
    ! Newton iteration: p and rho
    ! ---------------------------
    ! Intialization with momentum averaged quantities for u,v,w,Un and averaged quantities for rho,p,T
    Un_bar=0.0_wp; u_bar=0.0_wp; v_bar=0.0_wp; rho_bar=0.0_wp; p_bar=0.0_wp; T_bar=0.0_wp
    do i=1,n_pts-1
       Un_bar = Un_bar + (rho_line(i+1)*Un_line(i+1)*Un_line(i+1) + rho_line(i)*Un_line(i)*Un_line(i))*0.5_wp*dl_line(i)
       u_bar = u_bar + (rho_line(i+1)*Un_line(i+1)*u_line(i+1) + rho_line(i)*Un_line(i)*u_line(i))*0.5_wp*dl_line(i)
       v_bar = v_bar + (rho_line(i+1)*Un_line(i+1)*v_line(i+1) + rho_line(i)*Un_line(i)*v_line(i))*0.5_wp*dl_line(i)
       rho_bar = rho_bar + (rho_line(i+1) + rho_line(i))*0.5_wp*dl_line(i)
       p_bar = p_bar + (p_line(i+1) + p_line(i))*0.5_wp*dl_line(i)
       T_bar = T_bar + (T_line(i+1) + T_line(i))*0.5_wp*dl_line(i)
    enddo
    Un_bar=Un_bar/m_pt/l_line; u_bar=u_bar/m_pt/l_line; v_bar=v_bar/m_pt/l_line
    rho_bar=rho_bar/l_line; p_bar=p_bar/l_line; T_bar=T_bar/l_line

    ! Initialization of e and h
    e_bar = ecalc_pro(p_bar,rho_bar,T_bar)
    h_bar = e_bar + p_bar/rho_bar
    w_bar = theta_z/m_pt

    ! Newton loop
    ! -----------
    loop_Newton: do i=1,1000
       ! Initialize with previous iteration
       rho_test=rho_bar; p_test=p_bar

       ! print *,"rho_test",rho_test,"p_test",p_test

       ! function to cancel
       ! ------------------
       ! f = h + 0.5*(u**2 + v**2 + w**2) - h0
       fn = h_bar + 0.5_wp*(u_bar**2 + v_bar**2 + w_bar**2) - h0_mean

       ! print *,"fn",fn
       ! print *,"h0_mean",h0_mean
       ! print *,"h_bar + 0.5_wp*(u_bar**2 + v_bar**2 + w_bar**2)",h_bar + 0.5_wp*(u_bar**2 + v_bar**2 + w_bar**2)

       ! derivative against density
       ! --------------------------
       ! derivative of pressure w.r.t density
       dpdrho = -dpdvcalc_tro(T_bar,rho_test)/rho_test**2
       ! derivative of energy w.r.t density at constant pressure
       dedrho_p = dedrocalc_tro(T_bar,rho_test) - dedTcalc_tro(T_bar,rho_test)*dpdrho/dpdTcalc_tro(T_bar,rho_test)
       ! derivative of enthalpy w.r.t density at constant pressure
       dhdrho_p = dedrho_p - p_test/rho_test**2
       ! derivative w.r.t density
       dfdrho = dhdrho_p - m_pt**2/rho_test**3

       ! print *,"fn/dfdrho",fn/dfdrho

       ! update
       ! ------
       ! newton
       rho_bar = rho_test - fn/dfdrho


       ! print *,"rho_bar",rho_bar

       ! normal velocity from mass flow rate
       Un_bar = m_pt/rho_bar
       ! pressure from momentum
       p_bar = theta_n-Un_bar*m_pt
       ! velocity components from momentum components
       u_bar = (theta_x - p_bar*normal_x)/m_pt
       v_bar = (theta_y - p_bar*normal_y)/m_pt
       ! rest from pressure and density
       T_bar = tcalc_pro(p_bar,rho_bar,T_bar)
       e_bar = ecalc_pro(p_bar,rho_bar,T_bar)
       h_bar = e_bar + p_bar/rho_bar

       err1 = abs(p_bar-p_test)/p_test
       err2 = abs(rho_bar-rho_test)/rho_test

       if ((err1.le.1e-6_wp).and.(err2.le.1e-6_wp)) exit loop_Newton

    enddo loop_Newton



    ! ! For perfect gas
    ! ! ---------------
    ! ! Prasad (2015)
    ! e_bar=0
    ! do i=1,n_pts-1
    !    e_bar = e_bar + cpfg*(rho_line(i+1)*Un_line(i+1)*T_line(i+1) + rho_line(i)*Un_line(i)*T_line(i))*0.5_wp*dl_line(i) + &
    !            0.5_wp*(rho_line(i+1)*Un_line(i+1)*(u_line(i+1)**2 + v_line(i+1)**2 + w_line(i+1)**2) + &
    !                    rho_line(i  )*Un_line(i  )*(u_line(i  )**2 + v_line(i  )**2 + w_line(i  )**2))*0.5_wp*dl_line(i)
    ! enddo
    ! e_bar=e_bar/l_line

    ! ! Coefficients Q, L, C
    ! Q = 1/m_pt**2*(1.0_wp - 2*gam/(gam1))
    ! L = 2/m_pt**2*(gam/(gam1)*theta_n - theta_x*normal_x - theta_y*normal_y)
    ! C = 1/m_pt**2*(theta_x**2 + theta_y**2 + theta_z**2) - 2*e_bar/m_pt

    ! ! /!\ sub/supersonic solution in x-dir only
    ! p_bar = (-L - (L**2 - 4*Q*C)**0.5_wp)/2/Q ! subsonic
    ! ! p_bar = (-L + (L**2-4*Q*C)**0.5_wp)/2/Q ! supersonic
    ! Un_bar = (theta_n - p_bar)/m_pt
    ! u_bar = (theta_x - p_bar*normal_x)/m_pt
    ! v_bar = (theta_y - p_bar*normal_y)/m_pt
    ! w_bar = theta_z/m_pt
    ! rho_bar = m_pt/Un_bar
    ! T_bar = p_bar/rho_bar/rg
    ! M_bar = ((u_bar**2+v_bar**2+w_bar**2)/gam/rg/T_bar)**0.5_wp
    ! h_bar = e_bar + p_bar/rho_bar

    ! Total quantities of reference (inlet)
    ! -------------------------------------
    if (is_stagnation) then
       T_tot=T_ref; p_tot=p_ref
       rho_tot= rocalc_pt(p_tot,T_tot,rho_ref)
       H_tot= ecalc_tro(T_tot,rho_tot) + p_tot/rho_tot
       s_tot= scalc_tro(T_tot,rho_tot)
    else
       call mpistop('is_stagnation must be True, with inlet total quantities in param.ini',0)
    endif

    ! Compute different mixed-out values
    ! ----------------------------------
    c_bar = (c2calc_tro(T_bar,rho_bar))**0.5
    M_bar = sqrt(u_bar**2 + v_bar**2 + w_bar**2)/c_bar
    Gamma = gcalc_tro(T_bar,rho_bar)
    g_eq = cpcalc_tro(T_bar,rho_bar)/cvcalc_tro(T_bar,rho_bar)  ! gamma = Cp/Cv
    ! isentropic exponent: gamma_pv = - gamma * v/p * (dp/dv)_Tcst
    g_pv = -g_eq*dpdvcalc_tro(T_bar,rho_bar)/(rho_bar*p_bar)
    s_bar = scalc_tro(T_bar,rho_bar)

    ! Compute isentropic values
    ! -------------------------
    rho_is = rocalc_ps(p_bar,s_tot,T_tot)
    T_is = tcalc_sro(s_tot,rho_is,T_tot)
    e_is = ecalc_tro(T_is,rho_is)
    h_is = e_is + p_bar/rho_is
    if (h_is>H_tot) then
       print *,"/!\ h_is > H_tot",(h_is-H_tot)/H_tot
       call mpistop("/!\ h_is > H_tot",0)
    endif
    c_is = (c2calc_tro(T_is,rho_is))**0.5
    M_is = (2*(H_tot - h_is))**0.5/c_is

    ! Compute stagnation values
    ! -------------------------
    g_eq = cpcalc_tro(T_tot,rho_tot)/cvcalc_tro(T_tot,rho_tot)
    ! Initialized with gamma equivalent
    ! T0_bar = T_bar*(1 + (g_eq-1)/2 * M_bar**2)
    ! rho0_bar = rho_bar*(1 + (g_eq-1)/2 * M_bar**2)**(1/(g_eq-1))
    ! Initialized with static quantities
    T0_bar = T_bar; rho0_bar = rho_bar
    e0_bar=e_bar; p0_bar=p_bar; rho0_bar=rho_bar
    call stagnation_calc(h0_mean,s_bar,e0_bar,T0_bar,p0_bar,rho0_bar)

    ! Write mixed state
    ! -----------------
    ! Will contain:
    ! -> u, v, w, c, M, ....
    ! -> s: entropy
    ! -> g_pv: isentropic exponent
    ! -> Gamma: fundamental derivative of gas dynamics
    open(50,file='mixed_out_state.dat',form='formatted',position="append")
    rewind(50)
    write(50,"(A, 10F20.12)") "u   v   w   c   c_is   M   M_is   rho  rho0   rho_is   p   p0   T   T0   T_is   e   e0   e_is   h   h0  h_is   s   g_pv   Gamma"
    write(50,"(24F20.8)") u_bar, v_bar, w_bar, c_bar, c_is, M_bar, M_is, rho_bar, rho0_bar, rho_is, p_bar, p0_bar, T_bar, T0_bar, T_is, e_bar, e0_bar, e_is, h_bar, h0_mean, h_is, s_bar, g_pv, Gamma
    close(50)



  end subroutine pp_mixed_out_line

end module mod_pp_stats_main
