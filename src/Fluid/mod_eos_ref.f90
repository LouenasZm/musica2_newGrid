!=================================================================================
module mod_ineos_ref
!=================================================================================
  !> Module to share constants of PENG-ROBINSON EoS [for inlining purposes]
!=================================================================================
  use precision
  implicit none
  ! ------------------------------------------------------------------------------
  real(wp), parameter :: tol=1.e-6_wp
  ! variables for RefProp
  integer :: ierr
  character(len=255) :: herr
  ! ------------------------------------------------------------------------------
end module mod_ineos_ref

!=================================================================================
module mod_eos_ref
!=================================================================================
  !> Module to define subroutines for PENG-ROBINSON Equation of State (EoS)
!=================================================================================
  use mod_fluid
  use mod_ineos_ref
  use warnstop
  implicit none
  ! ------------------------------------------------------------------------------
  ! ------------------------------------------------------------------------------
  
contains
  
  !===============================================================================
  subroutine init_eos_ref
  !===============================================================================
    !> Initializations of gas coefficients
    !> - RefProp EOS -
    !===============================================================================
    use mod_mpi
    implicit none
    ! ----------------------------------------------------------------------------
    character(len=255) :: fileout,hf(1),hmix
    character(len=30)  :: fluid
    real(wp) :: ttp,acf
    ! TEST
    real(wp) :: p,p1,p2,cp,ro,e1,e2,h,T,c
    ! ----------------------------------------------------------------------------

    ! RefProp setup
    ! =============
    call SETPATH('../../src/Fluid/refprop/fluids')
    
    ! Call SETUP to initialize the program and set the pure fluid component name
    ! ==========================================================================
    ! call SETUP(1,trim(fluidname)//'.FLD','hmx.bnc','DEF',ierr,herr)
    
    hf(1)=trim(fluidname)//'.fld'

    hmix ='hmix.bnc'
    call SETUP(1,hf,hmix,'DEF',ierr,herr)
    if (ierr.ne.0) print *,herr
    
    ! Display RefProp info
    ! ====================
    
    call INFOr(1,pmol,ttp,teb,tc,pc,roc,zc,acf,mdm,rg)

    roc=roc*pmol

    if (iproc==0) then
       write (*,'(a,a)') 'RefProp properties for fluid: ', trim(fluidname)
       write (*,'(a,f11.4)') 'pmol [g/mol]  =', pmol
       write (*,'(a,f11.4)') 'acf  [-]      =', acf
       write (*,'(a,f11.4)') 'mdm  [-]      =', mdm
       write (*,'(a,f11.4)') 'Tc   [K]      =', Tc
       write (*,'(a,f11.4)') 'pc   [kPa]    =', pc
       write (*,'(a,f11.4)') 'rhoc [mol/L]  =', roc
       write (*,'(a,f11.4)') 'Rgas [J/mol-K]=', Rg
    endif
 
    ! Test RefProp routines
    ! =====================
    
    T =1.1_wp*Tc
    ro=1.1_wp*roc
    c =sqrt(c2calc_tro_ref(T,ro))
    p1=pcalc_tro_ref(T,ro)
    cp=cpcalc_tro_ref(T,ro)
    ro=rocalc_pt_ref(p1,T,1.0)
    e1=ecalc_tro_ref(T,ro)
    e2=ecalc_pro_ref(p1,ro,300.0)
    T =tcalc_roero_ref(ro*e1,ro,300.0)
    p2=pcalc_roero_ref(ro*e1,ro,300.0)
    h =e1 + p1/ro
    T =tcalc_ph_ref(p1,h,1.0,300.0)
    
    if (iproc==0) then
       print *,'T :', T
       print *,'ro:', ro
       print *,'cp:', cp
       print *,'c :', c
       print *,'p1:', p1
       print *,'p2:', p2
       print *,'e1:', e1
       print *,'e2:', e2
       print *,'h :', h
    endif

  end subroutine init_eos_ref

  !===============================================================================
  function avcalc_tro_ref(T,ro)
  !===============================================================================
    !> Compute isobaric expansion coefficient from T and rho
    !> - RefProp EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: avcalc_tro_ref ! OUTPUT
    ! ----------------------------------------------------------------------------

    call mpistop('RefProp function avcalc_tro not implemented',0)
    
  end function avcalc_tro_ref

  !===============================================================================
  function c2calc_tro_ref(T,ro)
  !===============================================================================
    !> Compute speed of sound from T and rho
    !> - RefProp EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: c2calc_tro_ref ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: d,p,e,h,s,cv,cp,w,hjt
    ! ----------------------------------------------------------------------------
    
    d=ro/pmol
    call THERM(T,d,1.0_wp,p,e,h,s,cv,cp,w,hjt)
    c2calc_tro_ref=w**2
    
  end function c2calc_tro_ref

  !===============================================================================
  function cpcalc_tro_ref(T,ro)
  !===============================================================================
    !> Compute heat capacity at constant pressure
    !> - RefProp EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: cpcalc_tro_ref ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: d,p,e,h,s,cv,cp,w,hjt
    ! ----------------------------------------------------------------------------
    
    d=ro/pmol
    call THERM(T,d,1.0_wp,p,e,h,s,cv,cp,w,hjt)
    cpcalc_tro_ref=cp*1000.0_wp/pmol
    
  end function cpcalc_tro_ref

  !===============================================================================
  function cvcalc_tro_ref(T,ro)
  !===============================================================================
    !> Compute heat capacity at constant volume
    !> - RefProp EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: cvcalc_tro_ref ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: d,p,e,h,s,cv,cp,w,hjt
    ! ----------------------------------------------------------------------------
    
    d=ro/pmol
    call THERM(T,d,1.0_wp,p,e,h,s,cv,cp,w,hjt)
    cvcalc_tro_ref=cv*1000.0_wp/pmol
    
  end function cvcalc_tro_ref

  !===============================================================================
  function dpdicalc_tro_ref(T,ro)
  !===============================================================================
    !> Compute pressure derivative w.r.t temperature
    !> /!\ not used for computation of saturation curve
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: dpdicalc_tro_ref ! OUTPUT
    ! ----------------------------------------------------------------------------

    call mpistop('RefProp function dpdicalc_tro not implemented',0)

  end function dpdicalc_tro_ref

  !===============================================================================
  subroutine dpcalc_tro_ref(T,ro,dpdro,dpdT)
  !===============================================================================
    !> Compute pressure derivatives w.r.t. density and temperature
    !> /!\ not used for computation of saturation curve
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in)  :: T,ro
    real(wp), intent(out) :: dpdro,dpdT
    ! ----------------------------------------------------------------------------

    call mpistop('RefProp function dpcalc_tro not implemented',0)

  end subroutine dpcalc_tro_ref

  !===============================================================================
  function ecalc_pro_ref(p,ro,Ttent)
  !===============================================================================
    !> Compute internal energy from p and rho
    !> - RefProp EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: p,ro,Ttent
    real(wp) :: ecalc_pro_ref ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: pp,D,T,Dl,Dv,x,y,q,e,h,s,Cv,Cp,w
    ! ----------------------------------------------------------------------------
    
    pp=p/1000.0_wp
    d =ro/pmol
    call PDFLSH(pp,D,1.0,T,Dl,Dv,x,y,q,e,h,s,Cv,Cp,w,ierr,herr)
    ecalc_pro_ref=e*1000.0_wp/pmol
    
  end function ecalc_pro_ref

  !===============================================================================
  function ecalc_tro_ref(T,ro)
  !===============================================================================
    !> Compute internal energy from T and rho
    !> - RefProp EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: ecalc_tro_ref ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: d,p,e,h,s,cv,cp,w,hjt
    ! ----------------------------------------------------------------------------
    
    d=ro/pmol
    call THERM(T,d,1.0_wp,p,e,h,s,cv,cp,w,hjt)
    ecalc_tro_ref=e*1000.0_wp/pmol
    
  end function ecalc_tro_ref

  !===============================================================================
  function gcalc_tro_ref(T,ro)
  !===============================================================================
    !> Compute fundamental derivative of gas dynamics from T and ro
    !> - RefProp EOS -
  !===============================================================================
    implicit none
    ! -----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: gcalc_tro_ref ! OUTPUT
    ! -----------------------------------------------------------------------------

    !call mpistop('RefProp function gcalc_tro not implemented',0)
    gcalc_tro_ref=9999.0_wp
    
  end function gcalc_tro_ref

  !===============================================================================
  function pcalc_roero_ref(roe,ro,Ttent)
  !===============================================================================
    !> Compute pressure from rhoe and rho
    !> - RefProp EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: roe,ro,Ttent
    real(wp) :: pcalc_roero_ref ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: D,e,T,P,Dl,Dv,x,y,q,h,s,Cv,Cp,w
    ! ----------------------------------------------------------------------------
    
    D=ro/pmol
    e=roe/ro/1000.0_wp*pmol
    call DEFLSH(D,e,1.0,T,P,Dl,Dv,x,y,q,h,s,Cv,Cp,w,ierr,herr)
    pcalc_roero_ref=p*1000.0_wp
    
  end function pcalc_roero_ref

  !===============================================================================
  function pcalc_tro_ref(T,ro)
  !===============================================================================
    !> Compute pressure from T and rho
    !> - RefProp EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: pcalc_tro_ref ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: d,p,e,h,s,cv,cp,w,hjt
    ! ----------------------------------------------------------------------------
    
    d=ro/pmol
    call THERM(T,d,1.0_wp,p,e,h,s,cv,cp,w,hjt)
    pcalc_tro_ref=p*1000.0_wp
    
  end function pcalc_tro_ref

  !===============================================================================
  function rocalc_ep_ref(e,p,Ttent)
  !===============================================================================
    !> Compute density from internal energy e and p
    !> - RefProp EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: e,p,Ttent
    real(wp) :: rocalc_ep_ref ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: ee,pp,T,D,Dl,Dv,x,y,q,h,s,Cv,Cp,w
    ! ----------------------------------------------------------------------------
    
    ee=e/1000.0_wp*pmol
    pp=p/1000.0_wp
    call PEFLSH(pp,ee,1.0,T,D,Dl,Dv,x,y,q,h,s,Cv,Cp,w,ierr,herr)
    rocalc_ep_ref=D*pmol
    
  end function rocalc_ep_ref

  !===============================================================================
  function rocalc_ps_ref(p,s,Ttent)
  !===============================================================================
    !> Compute density from p and entropy s
    !> - RefProp EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: p,s,Ttent
    real(wp) :: rocalc_ps_ref ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: pp,ss,T,D,Dl,Dv,x,y,q,e,h,Cv,Cp,w
    ! ----------------------------------------------------------------------------
    
    pp=p/1000.0_wp
    ss=s/1000.0_wp*pmol
    call PSFLSH(pp,ss,1.0,T,D,Dl,Dv,x,y,q,e,h,Cv,Cp,w,ierr,herr)
    rocalc_ps_ref=D*pmol
    
  end function rocalc_ps_ref

  !===============================================================================
  function rocalc_pt_ref(p,T,rotent)
  !===============================================================================
    !> Compute density from p and T
    !> - RefProp EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: p,T,rotent
    real(wp) :: rocalc_pt_ref ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: ro_guess,pp
    ! ----------------------------------------------------------------------------
    
    pp=p/1000.0_wp
    ro_guess=rotent/pmol
    call TPRHO(T,pp,1.0,2,0,ro_guess,ierr,herr)
    rocalc_pt_ref=ro_guess*pmol
    
  end function rocalc_pt_ref

  !===============================================================================
  function rocalc_st_ref(s,T,ro)
  !===============================================================================
    !> Compute density from entropy and T
    !> - RefProp EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: s,T,ro
    real(wp) :: rocalc_st_ref ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: ss,P,D,Dl,Dv,x,y,q,e,h,Cv,Cp,w
    ! ----------------------------------------------------------------------------
    
    ss=s/1000.0_wp*pmol
    call TSFLSH(T,ss,1.0,1,P,D,Dl,Dv,x,y,q,e,h,Cv,Cp,w,ierr,herr)
    rocalc_st_ref=D*pmol
    
  end function rocalc_st_ref

  !===============================================================================
  function scalc_tro_ref(T,ro)
  !===============================================================================
    !> Compute entropy s from T and rho
    !> - RefProp EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: scalc_tro_ref ! OUTPUT 
    ! ----------------------------------------------------------------------------
    real(wp) :: d,p,e,h,s,cv,cp,w,hjt
    ! ----------------------------------------------------------------------------
    
    d=ro/pmol
    call THERM(T,d,1.0_wp,p,e,h,s,cv,cp,w,hjt)
    scalc_tro_ref=s*1000.0_wp/pmol
    
  end function scalc_tro_ref

  !===============================================================================
  function tcalc_ph_ref(p,h,rotent,Ttent)
  !===============================================================================
    !> Compute temperature from p and enthalpy h
    !> - RefProp EOS - 
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: p,h,rotent,Ttent
    real(wp) :: tcalc_ph_ref ! OUTPUT
    !-----------------------------------------------------------------------------
    real(wp) :: pp,hh,T,D,Dl,Dv,x,y,q,e,s,Cv,Cp,w
    !-----------------------------------------------------------------------------
    
    pp=p/1000.0_wp
    hh=h/1000.0_wp*pmol
    call PHFLSH(pp,hh,1.0,T,D,Dl,Dv,x,y,q,e,s,Cv,Cp,w,ierr,herr)
    tcalc_ph_ref=T
    
  end function tcalc_ph_ref

  !===============================================================================
  function tcalc_pro_ref(p,ro,Ttent)
  !===============================================================================
    !> Compute temperature from p and rho
    !> - RefProp EOS - 
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: p,ro,Ttent
    real(wp) :: tcalc_pro_ref ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: pp,D,T
    ! ----------------------------------------------------------------------------
    
    pp=p/1000.0_wp
    D =ro/pmol
    call PDFL1(pp,D,1.0,T,ierr,herr)
    tcalc_pro_ref=T
    
  end function tcalc_pro_ref

  !===============================================================================
  function tcalc_roero_ref(roe,ro,Ttent)
  !===============================================================================
    !> Compute temperature from rho*e and rho
    !> - RefProp EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: roe,ro,Ttent
    real(wp) :: tcalc_roero_ref ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: D,e,T
    ! ----------------------------------------------------------------------------
    
    D=ro/pmol
    e=roe/ro/1000.0_wp*pmol
    call DEFL1(D,e,1.0,T,ierr,herr)
    tcalc_roero_ref=T
    
  end function tcalc_roero_ref

  !===============================================================================
  function tcalc_sro_ref(s,ro,Ttent)
  !===============================================================================
    !> Compute temperature from entropy s and rho
    !> - RefProp EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: s,ro,Ttent
    real(wp) :: tcalc_sro_ref ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: T,D,ss
    ! ----------------------------------------------------------------------------
    
    ss=s/1000.0_wp*pmol
    D =ro/pmol
    call DSFL1(D,ss,1.0,T,ierr,herr)
    tcalc_sro_ref=T
    
  end function tcalc_sro_ref

  !===============================================================================
  function vvol_ref(a,b,c)
  !===============================================================================
    !> Compute vvol such that p(vvol)=p
    !> /!\ not used for computation of saturation curve
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: a,b,c
    real(wp) :: vvol_ref ! OUTPUT
    ! ----------------------------------------------------------------------------

    call mpistop('RefProp function vvol not implemented',0)

  end function vvol_ref

  !===============================================================================
  function vvol_d1_ref(a,b)
  !===============================================================================
    !> Compute vvol_d1 from Ti and tentative v
    !> /!\ not used for computation of saturation curve
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: a,b
    real(wp) :: vvol_d1_ref ! OUTPUT
    ! ----------------------------------------------------------------------------

    call mpistop('RefProp function vvol_d1 not implemented',0)

  end function vvol_d1_ref

  !===============================================================================
  function vvol_d2_ref(a,b)
  !===============================================================================
    !> Compute vvol_d2 from Ti and tentative v
    !> /!\ not used for computation of saturation curve
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: a,b
    real(wp) :: vvol_d2_ref ! OUTPUT
    ! ----------------------------------------------------------------------------

    call mpistop('RefProp function vvol_d2 not implemented',0)

  end function vvol_d2_ref

  !===============================================================================
  function intpcalc_tro_ref(a,b)
  !===============================================================================
    !> Compute pressure integral from T and rho
    !> /!\ not used for computation of saturation curve
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: a,b
    real(wp) :: intpcalc_tro_ref ! OUTPUT
    ! ----------------------------------------------------------------------------

    call mpistop('RefProp function intpcalc_tro not implemented',0)

  end function intpcalc_tro_ref

end module mod_eos_ref
