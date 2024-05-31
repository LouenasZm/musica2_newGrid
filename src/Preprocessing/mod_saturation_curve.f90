!=================================================================================
module mod_saturation_curve
!=================================================================================
  !> Module to initialize HIT (Homogeneous Isotropic Turbulence)
!=================================================================================
  use mod_eos
  use mod_tranprop
  implicit none
  !-------------------------------------------------------------------------------
  integer :: nb_pts,i_stop,change_num,vnum,pnum
  integer :: icont,icont2,icont3,i1up,i2up,i1up2,i2up2,i1up3,i2up3
  integer :: ni,n_long
  real(wp) :: vmin,vmax,pmin,pmax
  real(wp) :: dpdi,p_tro,e_tro,g_tro,s_tro,e_pro,p_roero
  real(wp) :: ro_ep,ro_ps,ro_pt,ro_st,T_roero,T_sro,T_pro
  real(wp) :: T_tent,ro_tent
  real(wp), dimension(:), allocatable :: x_c,y_c,v_inf,v_sup,p_cur,p_satcurve
  real(wp), dimension(:,:), allocatable :: sat_curve
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------  

contains
  
  !===============================================================================
  subroutine pre_thermo
  !===============================================================================
    !> Subroutine to compute the saturation curve for different fluids
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer  :: choose, i
    real(wp) :: ro_in,T_in,p_in
    !!real(wp) :: zc
    ! ----------------------------------------------------------------------------

!!$    ! fluidname='r245fa'
!!$    ! fluidname='pp11'
!!$    ! fluidname='mdm'
!!$    ! visc_type = 'C'
!!$    ! fluidname='air'
!!$    ! visc_type = 'P'
!!$    ! Initialization of EoS and transport properties
!!$    if (trim(fluidname).eq.'air' .and. eos_type.ne.'pfg') then
!!$       print *,'For air only PFG model is possible.'
!!$       stop
!!$    endif
!!$    call init_eos
!!$    call init_viscosity

    print *,'Fluid and eos:', fluidname, eos_type

    ! Number of points
    ! ================
    nb_pts=1000
    
    ! Maximum and minimum reduced volume
    ! ==================================
    vnum=nb_pts
    vmin=0.55_wp
    vmin=0.4_wp
    vmin=0.1_wp
    vmax=15.0_wp
    !vmin=0.55_wp
    vmax=6.0_wp
    
    ! Maximum and minimum reduced pressure
    ! ====================================
    pnum=nb_pts
    pmin=0.4_wp
    !pmin=0.1_wp
    pmax=3.0_wp

    ! Simple test of values
    ! ---------------------
    p_in = 1.0_wp; T_in = 1.0_wp; ro_in = 1.0_wp

    if (trim(fluidname).eq.'air'.or.trim(fluidname).eq.'nitrogen') then
       T_in  = 298.15_wp
       ro_in = 1.0_wp
    elseif (fluidname(1:2).eq.'pp') then
       if (eos_type.eq.'pfg') then              ! VISCOUS PFG PP11 - THI
          T_in  = 1.01_wp*tc
          ro_in = 0.618_wp*roc
       elseif (eos_type.eq.'vdw') then          ! VISCOUS VDW PP11 - THI
          T_in  = 1.01_wp*tc
          ro_in = 0.618_wp*roc
       elseif (eos_type.eq.'mah') then          ! VISCOUS MAH PP11 - THI
          ! Viscous MAH - THI
          T_in  = 656.65_wp           !1.01_wp*tc
          ro_in = 387.56_wp           !0.618_wp*roc
          ! Viscous MAH2 - THI
          !T_in = 650.80015_wp
          !ro_in = 387.567715.0_wp
       endif
    elseif (trim(fluidname).eq.'novec649') then  ! VISCOUS PRS NOV - THI
      !T_in  = 1.01_wp*tc
      !ro_in = 0.618_wp*roc
      ! operating conditions A (nominal)
      T_in  = 373.15
      ro_in = 48.51

      ! operating conditions B
      T_in  = 408.15
      ro_in = 131.58

!!$      ! operating conditions C
!!$      T_in  = 433.15
!!$      ro_in = 263.16

      T_tent  = 1.01*tc
    elseif (trim(fluidname).eq.'mdm') then  ! VISCOUS PRS NOV - THI
      ! operating conditions TROVA
      T_in  = 542.13
      !ro_in = 251.709648550126
      !p_in = 904388.0

      ro_in =72.7678427548977

      T_tent  = 1.01*tc
   elseif (trim(fluidname).eq.'d6') then        ! VISCOUS SWN D6 - THI
      T_in  = 652.24_wp
      ro_in = 172.48_wp
   elseif (trim(fluidname).eq.'r245fa') then    ! VISCOUS SWP R245fa - THI
      T_in  = 431.43_wp
      ro_in = 318.94_wp
   elseif (trim(fluidname).eq.'r134a') then    ! VISCOUS SWP R245fa - THI
      T_in  = 431.43_wp
      ro_in = 318.94_wp
   elseif (trim(fluidname).eq.'r1233zde') then    ! VISCOUS SWP R245fa - THI
      T_in  = 439.6_wp
      ro_in = 480.2_wp
   endif

!!$   print *,tc,roc,roc0
   print *,'cv_c',cvcalc_tro(tc,1e-10)

   !call test_eos_Novec
   !call mpistop('',0)

    choose = 1 ! compute  p from T and ro
    ! choose = 2 ! compute  T from p and ro
    ! choose = 3 ! compute ro from T and p
    if (choose.eq.1) p_in = pcalc_tro(T_in,ro_in)
    if (choose.eq.2) T_in = tcalc_pro(p_in,ro_in,T_tent)
    if (choose.eq.3) ro_in= rocalc_pt(p_in,T_in,ro_tent)

    call test_values(T_in, ro_in, p_in)

    !call mpistop('',0)
    
    !call test_visc_Novec

    !call mpistop('',0)

    ! ---------------------------------------------------------------------------
    ! Computation of Saturation curve
    ! -------------------------------
    allocate(p_satcurve(nb_pts))
    p_satcurve = pmin
    if (eos_type.ne.'pfg') then
       if (eos_type.ne.'ref') then
          call compute_sat_curve
       else
          call compute_sat_curve_refprop
       endif
       
       ! Writing of saturation curve on file
       ! -----------------------------------
       print *,'Writing of saturation curve..'
       open(unit=73,file='satcurve_'//trim(fluidname)//'_'//trim(eos_type)//'.dat', status='replace')
       write(73,*) 'VARIABLES= "v/vc" "p/pc" "T/Tc"'
       write(73,'(a8,i5,a31)') ' ZONE I=',nb_pts,', T="Saturation curve"'

       do i=1,nb_pts
          write(73,'(3(1x,e21.14))') sat_curve(i,1),sat_curve(i,2),sat_curve(i,3)
       enddo
       close(73)
    endif
    
    ! ---------------------------------------------------------------------------
    ! Computation of Clapeyron diagram
    ! --------------------------------
    print *,'compute Clapeyron'
    
    vmin = max(vmin,minval(sat_curve(:,1)))
    vmax = min(vmax,maxval(sat_curve(:,1)))
    p_satcurve = sat_curve(:,2)
    call compute_clapeyron
    ! call compute_clapeyron_blasius
   
  end subroutine pre_thermo
 
  !===============================================================================
  subroutine test_eos_Novec
  !===============================================================================
    !> Check viscosity law for Novec649
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: n
    real(wp) :: err
    real(wp), dimension(5) :: T_ref,ro_ref,p_ref,cv_ref,cp_ref
    real(wp), dimension(5) :: p_,cv_,cp_
    ! ----------------------------------------------------------------------------

    ! Table 9 of McLinden et al. (J. Chem. Eng. Data, 2015, 60, 3646-3659)
    ! =========================
    ! temperature (K)
    T_ref=[250.0_wp,400.0_wp,442.0_wp,250.0_wp,450.0_wp]
    ! density (mol/l -> kg/m^3)
    ro_ref=[5.6_wp,4.0_wp,1.92_wp,0.001_wp,4.6_wp]*pmol
    ! pressure (MPa)
    p_ref=[11.459869_wp,2.9272365_wp,1.8757729_wp,0.0020718017_wp,42.305705_wp]*1.e6_wp
    ! isochoric heat capacity (J/mol/K -> J/kg/K)
    cv_ref=[277.136_wp,308.183_wp,351.688_wp,254.272_wp,319.295_wp]*3.16411!*pmol*1.0e-2_wp
    ! isobaric heat capacity (J/mol/K -> J/kg/K)
    cp_ref=[341.656_wp,386.271_wp,45430.0_wp,262.743_wp,365.304_wp]*3.16411!*pmol*1.0e-2_wp

    ! Check EoS law for Novec649
    ! ==========================
    do n=1,5
       p_(n) = pcalc_tro(T_ref(n),ro_ref(n))
       cp_(n)=cpcalc_tro(T_ref(n),ro_ref(n))
       cv_(n)=cvcalc_tro(T_ref(n),ro_ref(n))
       write(6,'(A)') '=========================================='
       write(6,'(A,i3,A,f8.3,A,f10.3)') 'point',n,' T=',T_ref(n),' ro=',ro_ref(n)
       err=(p_ref(n)-p_(n))/p_ref(n)*100.0_wp
       write(6,'(A,2(1x,f10.3),A,f8.3,A)') 'p :',1.e-6_wp*p_(n),1.e-6_wp*p_ref(n),' error:',err,'%'
       err=(cv_ref(n)-cv_(n))/cv_ref(n)*100.0_wp
       write(6,'(A,2(1x,f10.3),A,f8.3,A)') 'cv:',cv_(n),cv_ref(n),' error:',err,'%'
       err=(cp_ref(n)-cp_(n))/cp_ref(n)*100.0_wp
       write(6,'(A,2(1x,f10.3),A,f8.3,A)') 'cp:',cp_(n),cp_ref(n),' error:',err,'%'
    enddo

  end subroutine test_eos_Novec

  !===============================================================================
  subroutine test_visc_Novec
  !===============================================================================
    !> Check viscosity law for Novec649
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: n
    real(wp), dimension(9) :: ro,T,mu,mu_ref
    real(wp) :: ro_in,T_in,p_in,mu_in,mu_in_ref,dp
    ! ----------------------------------------------------------------------------

    ! Table 4 of Wen et al. (J. Chem. Eng. Data, 2017, 62, 3603-3609)
    ! =====================
    ! temperature (K)
    T=[250.0_wp,250.0_wp,250.0_wp,300.0_wp,300.0_wp,300.0_wp,350.0_wp,350.0_wp,350.0_wp]
    ! density (kg/m^3)
    ro=[0.0_wp,0.41_wp,1809.77_wp,0.0_wp,3.89_wp,1701.48_wp,0.0_wp,4.42_wp,1595.99_wp]
    ! dynamic viscosity (micro Pa.s)
    mu_ref=[8.09_wp,8.33_wp,2377.5_wp,9.77_wp,10.85_wp,1059.7_wp,11.43_wp,12.65_wp,587.87_wp]

!!$    ! temperature (K)
!!$    T=[250.0_wp,250.0_wp,250.0_wp,300.0_wp,300.0_wp,300.0_wp,350.0_wp,350.0_wp,350.0_wp]
!!$    ! density (kg/m^3)
!!$    ro=[0.0_wp,0.41_wp,1809.77_wp,0.0_wp,3.89_wp,1701.48_wp,0.0_wp,4.42_wp,1595.99_wp]
!!$    ! dynamic viscosity (micro Pa.s)
!!$    mu_ref=[8.09_wp,8.33_wp,2377.5_wp,9.77_wp,10.85_wp,1059.7_wp,11.43_wp,12.65_wp,587.87_wp]

    ! Check viscosity law for Novec649
    ! ================================
    do n=1,9
       if (ro(n)==0.0_wp) ro(n)=1.0e-6_wp
       mu(n)=viscosity_law(T(n),ro(n))
    enddo

    ! Display results at screen
    ! =========================
    do n=1,9
       print *,mu(n)*1.e6_wp,mu_ref(n)
    enddo

    mu_in_ref=18.3006940754905
    T_in =373.15_wp
    ro_in=48.51_wp
    mu_in=viscosity_law(T_in,ro_in)
    print *,mu_in*1.e6_wp,mu_in_ref,(mu_in*1.e6_wp-mu_in_ref)/mu_in_ref*100.0_wp

!!$    T_in =373.15_wp
!!$    ro_in=roc
!!$    open(73,file='visc_'//trim(fluidname)//'_'//trim(eos_type)//'.dat', status='replace')
!!$    dp=0.1e6_wp
!!$    p_in=0.0_wp
!!$    do n=1,200
!!$       p_in=p_in+dp
!!$       ro_in=rocalc_pt(p_in,T_in,ro_in)
!!$       mu_in=viscosity_law(T_in,ro_in)
!!$       print *,p_in*1e-6,ro_in,mu_in*1000.0_wp
!!$       write(73,'(3(1x,e21.14))') p_in,ro_in,mu_in
!!$    enddo
!!$    close(73)
!!$    stop

  end subroutine test_visc_Novec

  !===============================================================================
  subroutine test_values(T,ro,p)
  !===============================================================================
    !> Test simple values
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(inout) :: ro,T,p
    ! ----------------------------------------------------------------------------
    real(wp) :: s,g,c,mu,th,cvv,cpp,c2
    ! ----------------------------------------------------------------------------

    T_tent = tc
    ro_tent = roc
!!$    T=tc
!!$    ro=roc
!!$    T  = 542.13*1.05
!!$    ro = 251.709648550126
!!$    !ro = 230.0

!!$    print *,'cv',cvcalc_tro(378.9727445147505,1.e-10_wp)
!!$    print *,'cv',cvcalc_tro(398.4030926226526,1.e-10_wp)
!!$    print *,'cv',cvcalc_tro(418.8296559810146,1.e-10_wp)
!!$    print *,'cv',cvcalc_tro(440.3035116379542,1.e-10_wp)
!!$    print *,'cv',cvcalc_tro(462.8783554178455,1.e-10_wp)
!!$    print *,'cv',cvcalc_tro(486.6106361888495,1.e-10_wp)
!!$    print *,'cv',cvcalc_tro(511.5596970144863,1.e-10_wp)
!!$    print *,'cv',cvcalc_tro(537.7879235422057,1.e-10_wp)
!!$    print *,'cv',cvcalc_tro(565.3609000000000,1.e-10_wp)

    ! Critical quantities computation
    ! ===============================
    p  = pcalc_tro(T,ro)
    s = scalc_tro(t,ro)
    print *,t,ro,p,s
    g = gcalc_tro(t,ro)
    print *,g,c2calc_tro(t,ro)
    c = sqrt(c2calc_tro(t,ro))
    mu = viscosity_law(t,ro)
    th = thconductivity(mu,t,ro)

    ro_tent = 1.5_wp*roc
    T_tent  = 1.1_wp*tc

    if (eos_type.ne.'air') then
       print *,'Reduced T, p, ro'
       print *,T/tc,p/pc,ro/roc
    endif
    print *,'Dimensional T, p, ro'
    print *,T,p,ro

    ! Thermodynamic properties
    ! ========================
    c2     = sqrt(c2calc_tro(T,ro));          print *,'c2_tro  :', c2
    cpp    = cpcalc_tro(T,ro);                print *,'cp_tro  :', cpp
    cvv    = cvcalc_tro(T,ro);                print *,'cv_tro  :', cvv
    if (eos_type.ne.'ref') &
    dpdi   = dpdicalc_tro(T,ro);              print *,'dpdi_tro:', dpdi
    p_tro  = pcalc_tro(T,ro);                 print *,'p_Tro   :', p_tro
    e_tro  = ecalc_tro(T,ro);                 print *,'e_tro   :', e_tro
    g_tro  = gcalc_tro(T,ro);                 print *,'G_tro   :', g_tro
    s_tro  = scalc_tro(T,ro);                 print *,'s_tro   :', s_tro/(rg*tc)
    p_roero= pcalc_roero(ro*e_tro,ro,T_tent); print *,'p_roero :', p_roero
    !e_pro  = ecalc_pro(p,ro,T_tent);         print *,'e_pro   :', e_pro
    !ro_ep  = rocalc_ep(e_tro,p,T_tent);      print *,'ro_ep   :', ro_ep
    !ro_ps  = rocalc_ps(p_tro,s_tro,T_tent);  print *,'ro_ps   :', ro_ps
    ro_pt  = rocalc_pt(p,T,ro_tent);          print *,'ro_pt   :', ro_pt
    ro_st  = rocalc_st(s_tro,T,ro_tent);      print *,'ro_st   :', ro_st
    T_roero= tcalc_roero(ro*e_tro,ro,T_tent); print *,'T_roero :', T_roero
    T_sro  = tcalc_sro(s_tro,ro,T_tent);      print *,'T_sro   :', T_sro
    T_pro  = tcalc_pro(p,ro,T_tent);          print *,'T_pro   :', T_pro

    ! Transport properties
    ! ====================
    mu = viscosity_law(T,ro);                 print *,'mu_tro  :', mu
    th = thconductivity(mu,T,ro);             print *,'la_tro  :', th
    cpp = cpcalc_tro(t,ro);                   print *,'prandtl :', mu*cpp/th

    print *,'OK test values'

    call compute_clapeyron_single_point(T,ro)

  end subroutine test_values
    
  !===============================================================================
  subroutine compute_sat_curve
  !===============================================================================
    !> Computation of the saturation curve
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp) :: t
    ! real(wp) :: kap0,gam0,kap1,gam1,a_coef,b_coef,c_coef,d_coef
    ! real(wp) :: a_coef1,b_coef1,c_coef1,a_coef2,b_coef2,c_coef2
    ! integer  :: nb_sat_curve
    integer :: i,j,ii
    real(wp) :: value0,value1,value2,oldval2,val2_change,v1,v2
    real(wp) :: top_coeff,p_val,p_value
    real(wp) :: limit_inf,limit_sup,delta_pts
    ! saturation quantities for transport equation
    real(wp), dimension(:,:), allocatable :: sat_curve2
    ! External functions
    real(wp) :: deltav,x,y,tp
    integer :: i1up_sup,i2up_sup,i1up_inf,i2up_inf,isup,iinf,nbez
    real(wp), dimension(4) :: xbez,ybez
    real(wp), dimension(:), allocatable :: xinf,xsup,yinf,ysup
    ! ----------------------------------------------------------------------------

    ! Computation of left and right parts of the saturation curve
    ! ===========================================================
    t = 1.0_wp
    icont=0
    i_stop = 1
    top_coeff = 1.5_wp
    ! For VDW, use top_coeff ~ 1.1!
    ! For R245fa, use top_coeff = 1.55.0_wp
    if (eos_type.eq.'vdw') top_coeff = 1.1_wp
    if (eos_type.eq.'mah') top_coeff = 1.5_wp
    if (eos_type.eq.'swp') top_coeff = 1.55_wp
    if (eos_type.eq.'swn') top_coeff = 1.50_wp
    if (eos_type.eq.'prs') top_coeff = 1.1_wp
    ni=2000
    n_long=2000

    value0 = vvol_d2(t,vmin)
    !print *,'T, value0',t,value0
    !print *,top_coeff*value0
    value2 = vvol_d1(t,top_coeff*value0)

    allocate(v_inf(ni),v_sup(ni),p_cur(ni))
    
    ! Temperature cycle: under the saturation curve, the isotherms are horizontal lines
    ! =================================================================================
    ! N.B. vvol functions take dimensionless input parameters!

    iloop: do i=1,ni
       t = t - 0.0005_wp

       ! Computation of initialization value for the convergence loop
       value0 = vvol_d2(t,vmin)
       !print *,'T, value0',t,value0

       ! Check on the minimum pressure blocking the isotherms
       ! If the pressure corresponding to v=value0 is negative, the minimum is reached -> exit
       if (pcalc_tro(t*tc,1.0_wp/(value0*roc)).lt.0_wp) then !! Att BUG
          write(*,*), 'Minimum reached. T,v:', t,value0
          i_stop=i
          exit iloop
       endif

       oldval2 = value2
       value1=vvol_d1(t,vmin)
       !print *,'T,top_coeff',i,t,top_coeff
       value2=vvol_d1(t,top_coeff*value0)

       ! Check of top_coeff
       change_num = 0
       top_loop: do j=1,20

          if (i.ne.1) then
             val2_change=abs(value2-oldval2)/oldval2
          else
             val2_change=0.0_wp   ! Ignore on first iter
          endif

          if (value2.gt.0.0_wp.and.val2_change.lt.0.2_wp) then
             exit top_loop
          else
             top_coeff=top_coeff + 0.2_wp

             value2=vvol_d1(t,top_coeff*value0)

             change_num=change_num+1

             if (change_num.gt.20) then
                stop 'ERROR. Check top_coeff'
             endif
          endif
       enddo top_loop

       !print *,'value',value1,value2
       !print *,'p',pcalc_tro(t*tc,roc/value1),pcalc_tro(t*tc,roc/value2)

       ! Computation of pressure initialization value (pressure average)  ! modif XG 11/11/22
       if (i>100) then
          p_val=p_cur(i-1)
       else
          p_val=(pcalc_tro(t*tc,roc/value1) + &
                 pcalc_tro(t*tc,roc/value2))/2.0_wp/pc
       endif
       
       ! Check if average pressure value is under the zone under investigation
       !if (p_val.lt.0.001_wp) then ! modif XG 11/11/22
       if (p_val.lt.0.1_wp) then
          !print *,'Average pressure under 0.001, stop' ! modif XG 11/11/22
          print *,'Average pressure under 0.1, stop'
          i_stop=i
          exit iloop
       endif

       ! Loop for pressure and volume convergence with fixed temperature
       jloop: do j=1,1000

          ! Initialization
          p_value=p_val

          ! Computation of 2 volume values with fixed p and T
          v1=vvol(p_value,t,value1*0.8_wp)
          v2=vvol(p_value,t,value2*1.1_wp)
          !print *,i,j,v1,v2

          ! Computation of pressure values with new volume values
          p_val=(intpcalc_tro(T,1.0_wp/v2)-intpcalc_tro(T,1.0_wp/v1))/(v2-v1)

          !if (abs(p_val-p_value)/p_value.lt.1.e-10_wp) then
          if (abs(p_val-p_value)/p_value.lt.1.e-5_wp) then
             p_value=p_val
             icont=1

             p_cur(i)=p_value
             v_inf(i)=v1
             v_sup(i)=v2

             exit jloop
          endif
       enddo jloop
       !print *,i,'p_cur(i),v_sup(i)',p_cur(i),v_sup(i)

       if (icont.ne.1) stop 'Not converged in computation of saturation curve'
    enddo iloop

    ! ---------------------------------------------------------------------------

    if (i_stop==1) i_stop=ni+1
    print *,'i_stop',i_stop

    print *,'Upper and lower limits'
    print *,'Minimum Pres:',p_cur(i_stop-1)
    print *,'Minimum Temp:',t
    print *,'Minimum Volm:',v_inf(i_stop-1)
    print *,'Maximum Volm:',v_sup(i_stop-1)

    ! Computation of the top of the saturation curve
    ! ==============================================
    !
    ! Use of 3rd-order polynomial function
    ! ------------------------------------
    ! kap0=v_inf(1)**3.-3.*v_inf(1)+2.
    ! gam0=v_inf(1)**2.-2.*v_inf(1)+1.
    ! kap1=v_sup(1)**3.-3.*v_sup(1)+2.
    ! gam1=v_sup(1)**2.-2.*v_sup(1)+1.
    ! a_coef=(1./kap0*gam1-gam0*kap1)*(gam1-gam0)*(p_cur(i)-1.)
    ! b_coef=(1./kap0*gam1-gam0*kap1)(-kap1+kap0)*(p_cur(i)-1.)
    ! c_coef=-3.*a_coef-2.*b_coef
    ! d_coef=1.-a_coef-b_coef-c_coef
    !
    ! Use of 2nd-order polynomial function
    ! ------------------------------------
    ! kap0 = v_inf(1)**2-1.0_wp
    ! gam0 = v_inf(1)   -1.0_wp
    ! kap1 = v_sup(1)**2-1.0_wp
    ! gam1 = v_sup(1)   -1.0_wp
    ! a_coef = 1.0_wp/(kap0*gam1-gam0*kap1)*( gam1-gam0)*(p_cur(1)-1.0_wp)
    ! b_coef = 1.0_wp/(kap0*gam1-gam0*kap1)*(-kap1+kap0)*(p_cur(1)-1.0_wp)
    ! c_coef = 1.0_wp-a_coef-b_coef
    ! The coefficients below allow to ensure that the maximum is in (1,1)
    !a_coef1 = (p_cur(1)-1.0_wp)/(v_inf(1)-1.0_wp)**2
    !b_coef1 = -2.0_wp*a_coef1
    !c_coef1 = a_coef1 + 1.0_wp
    !a_coef2 = (p_cur(1)-1.0_wp)/(v_sup(1)-1.0_wp)**2
    !b_coef2 = -2.0_wp*a_coef2
    !c_coef2 = a_coef2 + 1.0_wp

    ! Fill array for saturation curve
    ! ===============================
    ii=0
    iinf=0

    allocate(x_c(n_long),y_c(n_long))
    allocate(xinf(n_long),xsup(n_long))
    allocate(yinf(n_long),ysup(n_long))
    do i=1,i_stop-1
       ii=ii+1
       x_c(ii)=v_inf(i_stop-i)
       y_c(ii)=p_cur(i_stop-i)
       
       iinf=iinf+1
       xinf(iinf)=v_inf(i_stop-i)
       yinf(iinf)=p_cur(i_stop-i)
    enddo

    !do i=1,5
    !  ! 2-nd order representation
    !  ii = ii+1
    !  x_c(ii) = v_inf(1) + i*(1.0_wp-v_inf(1))/6.0_wp
    !  !y_c(ii) = a_coef*x_c(ii)**2 + b_coef*x_c(ii) + c_coef
    !  y_c(ii) = a_coef1*x_c(ii)**2 + b_coef1*x_c(ii) + c_coef1
    !  !
    !  iinf = iinf+1
    !  xinf(iinf) = v_inf(1) + i*(1.0_wp-v_inf(1))/6.0_wp
    !  yinf(iinf) = a_coef1*xinf(iinf)**2 + b_coef1*x_c(iinf) + c_coef1
    !  print*, ii,iinf
    !  print*, x_c(ii),y_c(ii),yinf(iinf)
    !enddo
    nbez=20
    deltav=v_inf(1)-v_sup(1)

    xbez(1)=v_inf(2)
    xbez(2)=v_inf(1)
    xbez(3)=1.0_wp - (1.0_wp-v_inf(1))/2.0_wp
    xbez(4)=1.0_wp
    ybez(1)=p_cur(2)
    ybez(2)=p_cur(1)
    ybez(3)=1.0_wp
    ybez(4)=1.0_wp
    
    ! The last point is used as Bezier control point
    ii=ii-1
    do i=1,nbez
       ii=ii+1
       tp=dble(i)/dble(nbez+1)
       call bez3(tp,xbez,ybez,x,y)
       x_c(ii)=x
       y_c(ii)=y
    enddo

    ii=ii+1
    x_c(ii)=1.0_wp
    y_c(ii)=1.0_wp
    
    iinf=iinf+1
    xinf(iinf)=1.0_wp
    yinf(iinf)=1.0_wp
    
    isup=1
    xsup(1)=1.0_wp
    ysup(1)=1.0_wp

    xbez(1)=1.0_wp
    xbez(2)=1.0_wp + (v_sup(1)-1.0_wp)/2.0_wp
    xbez(3)=v_sup(1)
    xbez(4)=v_sup(2)
    ybez(1)=1.0_wp
    ybez(2)=1.0_wp
    ybez(3)=p_cur(1)
    ybez(4)=p_cur(2)

    do i=1,nbez
       ii=ii+1
       tp=dble(i)/dble(nbez+1)
       call bez3(tp,xbez,ybez,x,y)
       x_c(ii)=x
       y_c(ii)=y
    enddo
    !do i=1,5
    !  ! 2-nd order representation
    !  ii = ii+1
    !  x_c(ii) = 1.0_wp + i*(v_sup(1)-1.0_wp)/6.0_wp
    !  y_c(ii) = a_coef2*x_c(ii)**2 + b_coef2*x_c(ii) + c_coef2
    !  !
    !  isup = isup+1
    !  xsup(isup) = 1.0_wp + i*(v_sup(1)-1.0_wp)/6.0_wp
    !  ysup(isup) = a_coef2*xsup(isup)**2 + b_coef2*xsup(isup) + c_coef2
    !enddo

    ! The first point is used as Bezier control point
    do i=2,i_stop-1
       ii=ii+1
       x_c(ii)=v_sup(i)
       y_c(ii)=p_cur(i)

       isup=isup+1
       xsup(isup)=v_sup(i)
       ysup(isup)=p_cur(i)
    enddo

    ! Spline limits
    ! -------------
    i1up=1
    i2up=ii
    i1up_inf=1
    i2up_inf=iinf
    i1up_sup=1
    i2up_sup=isup

    ! Allocation of sat_curve
    ! -----------------------
    if (allocated(sat_curve)) deallocate(sat_curve)
    allocate(sat_curve(nb_pts,3))
    sat_curve=0.0_wp
    ! sat_curve(:,1)=v_sat
    ! sat_curve(:,2)=p_sat
    ! sat_curve(:,3)=T_sat
    allocate(sat_curve2(nb_pts,2))
    sat_curve2=0.0_wp

!!$    ! Fill with raw data
!!$    ! ==================
!!$    nb_sat_curve = icont3
!!$    allocate(sat_curve(nb_sat_curve,2))
!!$    sat_curve = 0.0_wp
!!$    do i=1,nb_sat_curve
!!$       sat_curve(i,1) = x_c(i)
!!$       sat_curve(i,2) = y_c(i)
!!$    enddo

    ! Fill with spline
    ! ----------------
    !print *,ii,x_c(1),x_c(ii)

    limit_inf=max(vmin,x_c(1))
    limit_sup=min(vmax,x_c(ii))

    delta_pts=(limit_sup-limit_inf)/dble(nb_pts-1)

    ii=0
    do i=1,nb_pts
       sat_curve(i,1)=limit_inf + (i-1)*delta_pts
       sat_curve(i,2)=spline(sat_curve(i,1),i1up,i2up,x_c,y_c,n_long)
       sat_curve(i,3)=tcalc_pro(sat_curve(i,2)*pc, 1.0_wp/sat_curve(i,1)*roc, Tc)/tc

       ! viscosity
       sat_curve2(i,1)=viscosity_law(sat_curve(i,3)*tc,1.0_wp/sat_curve(i,1)*roc)
       ! thermal conductivity
       sat_curve2(i,2)=thconductivity(sat_curve2(i,1),sat_curve(i,3)*tc,1.0_wp/sat_curve(i,1)*roc)

       if ((ii.eq.0).and.(sat_curve(i,1).ge.1.0_wp)) ii = i
    enddo

    ! Write dimensional saturation temperature T_sat
    ! ----------------------------------------------
    open(unit=73,file='Tsat_'//trim(fluidname)//'_'//trim(eos_type)//'.dat', status='replace')
    write(73,*) 'VARIABLES= "v" "T"'
    write(73,'(a8,i5,a31)') ' ZONE I=',nb_pts-ii+1,', T="Saturation curve"'
    do i=ii,nb_pts
       write(73,'(2(1x,e21.14))') sat_curve(i,1)/roc,sat_curve(i,3)*tc
    enddo
    close(73)
    
    open(unit=73,file='Tsat_mu_'//trim(fluidname)//'_'//trim(eos_type)//'.dat', status='replace')
    write(73,*) 'VARIABLES= "mu" "T"'
    write(73,'(a8,i5,a31)') ' ZONE I=',nb_pts-ii+1,', T="Saturation curve"'
    do i=ii,nb_pts
       write(73,'(2(1x,e21.14))') sat_curve2(i,1),sat_curve(i,3)*tc
    enddo
    close(73)

    open(unit=73,file='Tsat_k_'//trim(fluidname)//'_'//trim(eos_type)//'.dat', status='replace')
    write(73,*) 'VARIABLES= "k" "T"'
    write(73,'(a8,i5,a31)') ' ZONE I=',nb_pts-ii+1,', T="Saturation curve"'
    do i=ii,nb_pts
       write(73,'(2(1x,e21.14))') sat_curve2(i,2),sat_curve(i,3)*tc
    enddo
    close(73)

  end subroutine compute_sat_curve

  !===============================================================================
  subroutine compute_sat_curve_refprop
  !===============================================================================
    !> Computation of the saturation curve using REFPROP
  !===============================================================================
    use mod_ineos_ref
    use mod_ineos_ref
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,kr
    real(wp) :: D,v,T,P,Dl,Dv,x,y,deltav
    ! saturation quantities for transport equation
    real(wp), dimension(:,:), allocatable :: sat_curve2
    ! ----------------------------------------------------------------------------
    external SATD
    ! ----------------------------------------------------------------------------

    allocate(sat_curve(nb_pts,3),sat_curve2(nb_pts,2))
    sat_curve =0.0_wp
    sat_curve2=0.0_wp

    deltav=(vmax-vmin)/dble(vnum-1)

    do i=1,vnum
       v=(vmin+(i-1)*deltav)/roc

       D=1.0_wp/v/pmol
       call SATD(D,1.0_wp,0,kr,T,P,Dl,Dv,x,y,ierr,herr)
       !p=p*1000.0_wp
       sat_curve(i,1)=v*roc
       sat_curve(i,2)=p/pc
       sat_curve(i,3)=T/Tc

       ! viscosity
       sat_curve2(i,1)=viscosity_law(sat_curve(i,3)*Tc,1.0_wp/sat_curve(i,1)*roc)
       ! thermal conductivity
       sat_curve2(i,2)=thconductivity(sat_curve2(i,1),sat_curve(i,3)*tc,1.0_wp/sat_curve(i,1)*roc)
    enddo

    ! Write dimensional saturation temperature T_sat
    ! ----------------------------------------------
    open(unit=73,file='Tsat_'//trim(fluidname)//'_'//trim(eos_type)//'.dat', status='replace')
    write(73,*) 'VARIABLES= "v" "T"'
    write(73,'(a8,i5,a31)') ' ZONE I=',nb_pts,', T="Saturation curve"'
    do i=1,nb_pts
       write(73,'(2(1x,e21.14))') sat_curve(i,1)/roc,sat_curve(i,3)*tc
    enddo
    close(73)

    open(unit=73,file='Tsat_mu_'//trim(fluidname)//'_'//trim(eos_type)//'.dat', status='replace')
    write(73,*) 'VARIABLES= "mu" "T"'
    write(73,'(a8,i5,a31)') ' ZONE I=',nb_pts,', T="Saturation curve"'
    do i=1,nb_pts
       write(73,'(2(1x,e21.14))') sat_curve2(i,1),sat_curve(i,3)*tc
    enddo
    close(73)

    open(unit=73,file='Tsat_k_'//trim(fluidname)//'_'//trim(eos_type)//'.dat', status='replace')
    write(73,*) 'VARIABLES= "k" "T"'
    write(73,'(a8,i5,a31)') ' ZONE I=',nb_pts,', T="Saturation curve"'
    do i=1,nb_pts
       write(73,'(2(1x,e21.14))') sat_curve2(i,2),sat_curve(i,3)*tc
    enddo
    close(73)

  end subroutine compute_sat_curve_refprop

  !===============================================================================
  subroutine compute_clapeyron
  !===============================================================================
    !> Computation of the Clapeyron diagram
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer  :: i,j
    real(wp) :: p,v,T,ro,e,g,c,s,mu,th,cvv,cpp,adim,deltap,deltav,p_c
    ! ----------------------------------------------------------------------------

    ! /!\ if RefProp, pc is in [kPa]
    if (eos_type=='ref') then
       p_c=1000.0_wp*pc
    else
       p_c=pc
    endif

    print *,'Clapeyron diagram computation'
    open(unit=14,file='pv_'//trim(fluidname)//'_'//trim(eos_type)//'.dat', status='replace')
    write(14,*) 'TITLE="Clapeyron Diagram"'
    write(14,'(A)')'VARIABLES = "v/v_c" "p/p_c" "T/T_c" "(s-s_c)/RT_c" "e/RT_c" "`G" "c" "c_v" "c_p" "`m" "`l" "Pr"'
    write(14,'(a16,i5,a4,i5)') 'ZONE F=POINT, I=',vnum,' ,J=',pnum

    adim=pc/roc
    deltav=(vmax-vmin)/dble(vnum-1)

    print *,'tc,rg',tc,rg,pc

    do j=1,pnum
       do i=1,vnum
          v=(vmin+(i-1)*deltav)/roc
          ro=1.0_wp/v
          pmin=p_satcurve(i)
          deltap=(pmax-pmin)/dble(pnum-1)

          p=(pmin+(j-1)*deltap)*p_c
          T=tcalc_pro(p,ro,tc)
          e=ecalc_tro(t,ro)
          g=gcalc_tro(t,ro)
          c=sqrt(c2calc_tro(t,ro))
          s=scalc_tro(t,ro)
          e=ecalc_tro(t,ro)
          cvv=cvcalc_tro(t,ro)
          cpp=cpcalc_tro(t,ro)
          mu=viscosity_law(t,ro)
          th=thconductivity(mu,t,ro)
          write(14,'(1p11(E19.12,1X))') v*roc,p/p_c,t/tc,s/(rg*tc),e/(rg*Tc),g, &
                                        c/sqrt(adim),cvv/rg,cpp/rg,mu,th,mu*cpp/th
       enddo
    enddo

    close(14)
    
  end subroutine compute_clapeyron

!!$  !===============================================================================
!!$  subroutine compute_clapeyron_blasius
!!$  !===============================================================================
!!$    !> Computation of the Clapeyron diagram
!!$  !===============================================================================
!!$    implicit none
!!$    use mod_blasius
!!$    use mod_mpi
!!$    ! ----------------------------------------------------------------------------
!!$    integer :: istart
!!$    real(wp) :: v,th_ref,cvv,cpp,adim,deltap,deltav
!!$    ! ----------------------------------------------------------------------------
!!$
!!$    ! Allocation for Blasius program
!!$    allocate(u_bl(neta),v_bl(neta),rho_bl(neta),T_bl(neta),mu_bl(neta),G_bl(neta),M_bl(neta))
!!$    allocate(cp_bl(neta),lambda_bl(neta),Pr_bl(neta),feta(neta) )
!!$    allocate(eta(neta),etai(neta))
!!$
!!$    adim=pc/roc
!!$
!!$    deltav=(vmax-vmin)/dble(vnum-1)
!!$
!!$    iloop: do i=1, vnum
!!$       v = (vmin + (i-1)*deltav)
!!$       if (v.gt.1.05_wp) then
!!$          istart = i
!!$          exit iloop
!!$       endif
!!$    enddo iloop
!!$
!!$    print *,'Clapeyron diagram computation'
!!$    open(unit=14,file='pv_'//trim(fluidname)//'_'//trim(eos_type)//'_blasius.dat', status='replace')
!!$    write(14,*) 'TITLE="Clapeyron Diagram"'
!!$    write(14,'(A)')'VARIABLES = "v/v_c" "p/p_c" "T/T_c" "(s-s_c)/RT_c" "`G" "c"&
!!$         & "c_v" "c_p" "`m" "`l" "Pr" "`G_w" "`G_M_1"'
!!$    write(14,'(a16,i5,a4,i5)') 'ZONE F=POINT, I=',vnum-istart+1,' ,J=',pnum
!!$
!!$    Mach = 6.0_wp
!!$
!!$    do j=1,pnum
!!$       do i=istart,vnum
!!$          v = (vmin + (i-1)*deltav)/roc
!!$          rho_ref = 1.0_wp/v
!!$          pmin = p_satcurve(i)
!!$          deltap = (pmax-pmin)/dble(pnum-1)
!!$
!!$          p_ref = (pmin + (j-1)*deltap)*pc
!!$          T_ref = tcalc_pro(p_ref,rho_ref,tc)
!!$          e_ref = ecalc_tro(t_ref,rho_ref)
!!$          g_ref = gcalc_tro(t_ref,rho_ref)
!!$          c_ref = sqrt(c2calc_tro(t_ref,rho_ref))
!!$          s_ref = scalc_tro(t_ref,rho_ref)
!!$          cvv= cvcalc_tro(t_ref,rho_ref)
!!$          cpp= cpcalc_tro(t_ref,rho_ref)
!!$          mu_ref = viscosity_law(t_ref,rho_ref)
!!$          th_ref = thconductivity(mu_ref,t_ref,rho_ref)
!!$
!!$          write(*,'(A,i0,2X,i0,2X,f12.5,2X,f12.5)') 'I,J, V/Vc, P/Pc: ' &
!!$               , i, j, v*roc, p_ref/pc
!!$          call compute_blasius
!!$
!!$          if (is_nan) then
!!$             G_w  = -1.0_wp
!!$             G_M1 = -1.0_wp
!!$          endif
!!$          write(14,'(13(1pE19.12,1X))') v*roc, p_ref/pc, t_ref/tc, s_ref/(rg*tc), g_ref, &
!!$               c_ref/sqrt(adim), cvv/rg, cpp/rg,    &
!!$               mu_ref, th_ref, mu_ref*cpp/th_ref, G_w, G_M1
!!$       enddo
!!$    enddo
!!$    close(14)
!!$    
!!$  end subroutine compute_clapeyron_blasius

  !===============================================================================
  subroutine compute_clapeyron_single_point(T,ro)
  !===============================================================================
    !> Computation of the Clapeyron diagram
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    ! ----------------------------------------------------------------------------
    real(wp) :: p,v,e,g,c,s,mu,th,cvv,cpp,adim
    ! ----------------------------------------------------------------------------

    adim=pc/roc

    print *,'Initial conditions'
    open(unit=14,file='pv_'//trim(fluidname)//'_'//trim(eos_type)//'_ic.dat', status='replace')
    write(14,'(A)')'VARIABLES = "v/v_c" "p/p_c" "T/T_c" "(s-s_c)/RT_c" "`G" "c" "c_v" "c_p" "`m" "`l" "Pr"'

    v=1.0_wp/ro
    p=pcalc_tro(T,ro)
    e=ecalc_tro(t,ro)
    g=gcalc_tro(t,ro)
    c=sqrt(c2calc_tro(t,ro))
    s=scalc_tro(t,ro)
    cvv=cvcalc_tro(t,ro)
    cpp=cpcalc_tro(t,ro)
    mu=viscosity_law(t,ro)
    th=thconductivity(mu,t,ro)
    
    write(14,'(1p11(E19.12,1X))') v*roc,p/pc,t/tc,s/(rg*tc),g,c/sqrt(adim), &
                                  cvv/rg,cpp/rg,mu,th,mu*cpp/th
    close(14)

  end subroutine compute_clapeyron_single_point

  !===============================================================================
  function spline(xcur,i1up,i2up,xpro,ypro,n_long)
  !===============================================================================
    !> Spline computation
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer, intent(in) :: i1up,i2up,n_long
    real(wp), intent(in) :: xcur
    real(wp), dimension(n_long), intent(in) :: xpro,ypro
    real(wp) :: spline ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer  :: i,k,k1,k2
    real(wp) :: y
    real(wp), dimension(n_long) :: a,b,c,d,r,t
    ! ----------------------------------------------------------------------------

    a(1)=1.
    b(1)=1.
    c(1)=1.
    d(1)=1.
    do i=i1up+1,i2up-1
       a(i)= xpro(i)-xpro(i-1)
       b(i)= 2.0_wp*(xpro(i+1)-xpro(i-1))
       c(i)= xpro(i+1)-xpro(i)
       d(i)= 6.0_wp*((ypro(i+1)-ypro(i))/(xpro(i+1)-xpro(i))+ &
                     (ypro(i-1)-ypro(i))/(xpro(i)-xpro(i-1)))
    enddo

    k1=i1up
    k2=i2up

    r(k2)=0.0_wp
    t(k2)=0.0_wp

    do k=k2-1,k1,-1
       r(k)=-a(k)/(c(k)*r(k+1)+b(k))
       t(k)=(d(k)-c(k)*t(k+1))/(c(k)*r(k+1)+b(k))
    enddo

    d(k1)=0.0_wp

    do k=k1,k2-1
       d(k+1)=r(k+1)*d(k)+t(k+1)
    enddo

    do i=i1up+1,i2up-1
       if ((xcur>=xpro(i)).and.(xcur<=xpro(i+1))) then
          y = (d(i)*((xpro(i+1)-xcur)**3)+d(i+1)*((xcur-xpro(i))**3))/  &
              (6.0_wp*(xpro(i+1)-xpro(i)))+(ypro(i)/(xpro(i+1)-xpro(i)) &
              -d(i)*(xpro(i+1)-xpro(i))/6.0_wp)*  &
               (xpro(i+1)-xcur)+(ypro(i+1)/(xpro(i+1)-xpro(i))-d(i+1)*  &
               (xpro(i+1)-xpro(i))/6.0_wp)*(xcur-xpro(i))
       endif
    enddo

    spline=y

  end function spline

  !===============================================================================
  subroutine bez3(t,px,py,x,y)
  !===============================================================================
    !> Bezier curve computation
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! Input:
    ! ------
    ! t variable in the interval [0,1]
    real(wp), intent(in) :: t
    ! px,py coordinates of the four control points
    real(wp), dimension(4), intent(in) :: px,py
    ! Output:
    ! -------
    ! x,y returned point on the bezier curve
    real(wp), intent(out) :: x,y
    ! ----------------------------------------------------------------------------
    real(wp) :: b0,b
    ! ----------------------------------------------------------------------------

    b0= 1-t
    x= 0.0_wp
    y= 0.0_wp
    
    b= b0*b0*b0
    x= x + b*px(1)
    y= y + b*py(1)
    
    b= 3.0_wp*b0*b0*t
    x= x + b*px(2)
    y= y + b*py(2)
    
    b= 3.0_wp*b0*t*t
    x= x + b*px(3)
    y= y + b*py(3)
    
    b= t*t*t
    x= x + b*px(4)
    y= y + b*py(4)
    
  end subroutine bez3

end module mod_saturation_curve
