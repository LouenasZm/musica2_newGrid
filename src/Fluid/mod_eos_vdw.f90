!=================================================================================
module mod_ineos_vdw
!=================================================================================
  !> Module to share constants of VAN DER WAALS EoS [for inlining purposes]
!=================================================================================
  use precision
  implicit none
  ! ------------------------------------------------------------------------------
  real(wp) :: avw,bvw
  ! ------------------------------------------------------------------------------
end module mod_ineos_vdw
  
!=================================================================================
module mod_eos_vdw
!=================================================================================
  !> Module to define subroutines for VAN DER WAALS Equation of State (EoS)
!=================================================================================
  use mod_fluid
  use mod_ineos_vdw
  use warnstop
  implicit none
  ! ------------------------------------------------------------------------------
  real(wp), private, parameter :: tol=1.e-6_wp
  ! ------------------------------------------------------------------------------
  
contains
  
  !===============================================================================
  subroutine init_eos_vdw
  !===============================================================================
    !> Initializations of gas coefficients
    !> - VDW EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp) :: s_crit
    ! ----------------------------------------------------------------------------

    call read_eos
    
    ! The critical parameters of the feos.ini file are corrected in order to be
    ! compatible with the Van der Waals equation:
    ! -------------------------------------------
    zc  = 0.375_wp
    roc0=roc
    roc = pc/(zc*rg*Tc)

    ! constants van der Waals
    ! -----------------------
    ! molecular interaction
    avw = 27.0_wp/64.0_wp*(rg*tc)**2/pc
    ! covolume
    bvw = rg*tc/(8.0_wp*pc)
    
    ! gamma - 1
    ! ---------
    gam1= gam-1.0_wp
    
    ! specific heats
    ! --------------
    cvfg= 1.0_wp/gam1*rg
    cpfg= gam/gam1*rg

    ! critical entropy 
    ! ----------------
    s_crit = scalc_tro_vdw(Tc,roc)

  end subroutine init_eos_vdw

  !===============================================================================
  function avcalc_tro_vdw(T,ro)
  !===============================================================================
    !> Compute isobaric expansion coefficient from T and rho
    !> - VDW EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: avcalc_tro_vdw ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: dpdv,dpdT,v
    ! ----------------------------------------------------------------------------

    ! specific volume
    v = 1.0_wp/ro

    dpdv=-rg*T/(v-bvw)**2 + 2.0_wp*avw/v**3
    dpdT= rg*ro/(1.0_wp-ro*bvw)

    avcalc_tro_vdw= -ro*dpdT/dpdv

  end function avcalc_tro_vdw

  !===============================================================================
  function c2calc_tro_vdw(T,ro)
  !===============================================================================
    !> Compute speed of sound from T and rho
    !> - VDW EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: c2calc_tro_vdw ! OUTPUT
    ! ----------------------------------------------------------------------------

    c2calc_tro_vdw= gam*rg*T/(1.0_wp-ro*bvw)**2-2.0_wp*avw*ro
    
  end function c2calc_tro_vdw

  !===============================================================================
  function cpcalc_tro_vdw(T,ro)
  !===============================================================================
    !> Compute heat capacity at constant pressure
    !> - VDW EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: cpcalc_tro_vdw ! OUTPUT
    ! ----------------------------------------------------------------------------

    cpcalc_tro_vdw= cpfg
    !cpcalc = Rvw/(gam - 1.0_wp) + Rvw/(1.0_wp-2.0_wp*ro*avw/(Rvw*tr)*(1.0_wp-ro*bvw)**2)
    
  end function cpcalc_tro_vdw

  !===============================================================================
  function cvcalc_tro_vdw(T,ro)
  !===============================================================================
    !> Compute heat capacity at constant volume
    !> - VDW EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: cvcalc_tro_vdw ! OUTPUT
    ! ----------------------------------------------------------------------------
    
    cvcalc_tro_vdw= cvfg
    
  end function cvcalc_tro_vdw

  !===============================================================================
  function dpdicalc_tro_vdw(T,ro)
  !===============================================================================
    !> Compute pressure derivative w.r.t temperature
    !> - VDW EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: dpdicalc_tro_vdw ! OUTPUT
    ! ----------------------------------------------------------------------------

    dpdicalc_tro_vdw= gam1*ro/(1.0_wp-ro*bvw)
    
  end function dpdicalc_tro_vdw

  !===============================================================================
  function dpdTcalc_tro_vdw(T,ro)
  !===============================================================================
    !> Compute pressure derivative w.r.t temperature
    !> - VDW EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: dpdTcalc_tro_vdw ! OUTPUT
    ! ----------------------------------------------------------------------------

    dpdTcalc_tro_vdw= rg*ro/(1.0_wp-ro*bvw)

  end function dpdTcalc_tro_vdw

  !===============================================================================
  function dpdvcalc_tro_vdw(T,ro)
  !===============================================================================
    !> Compute first-order pressure derivative w.r.t volume
    !> - VDW EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: dpdvcalc_tro_vdw ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: v
    ! ----------------------------------------------------------------------------

    v = 1.0_wp/ro
    dpdvcalc_tro_vdw= -rg*T/(v-bvw)**2 + 2.0_wp*avw/v**3
    
  end function dpdvcalc_tro_vdw

  !===============================================================================
  function d2pdv2calc_tro_vdw(T,ro)
  !===============================================================================
    !> Compute second-order pressure derivative w.r.t volume
    !> - VDW EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: d2pdv2calc_tro_vdw ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: v
    ! ----------------------------------------------------------------------------

    v = 1.0_wp/ro
    d2pdv2calc_tro_vdw= 2.0_wp*rg*T/(v-bvw)**3 - 6.0_wp*avw/v**4
    
  end function d2pdv2calc_tro_vdw

  !===============================================================================
  function d3pdv3calc_tro_vdw(T,ro)
  !===============================================================================
    !> Compute third-order pressure derivative w.r.t volume
    !> - VDW EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: d3pdv3calc_tro_vdw ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: v
    ! ----------------------------------------------------------------------------

    v = 1.0_wp/ro
    d3pdv3calc_tro_vdw= -6.0_wp*rg*T/(v-bvw)**4 + 24.0_wp*avw/v**5
    
  end function d3pdv3calc_tro_vdw

  !===============================================================================
  function ecalc_pro_vdw(p,ro,Ttent)
  !===============================================================================
    !> Compute internal energy from p and rho
    !> (the tentative temperature is not used for vdw)
    !> - VDW EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: p,ro,Ttent
    real(wp) :: ecalc_pro_vdw ! OUTPUT
    ! ----------------------------------------------------------------------------

    ecalc_pro_vdw= (p+ro**2*avw)*(1.0_wp-ro*bvw)/(ro*gam1) - avw*ro
    
  end function ecalc_pro_vdw

  !===============================================================================
  function ecalc_tro_vdw(T,ro)
  !===============================================================================
    !> Compute internal energy from T and rho
    !> - VDW EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: ecalc_tro_vdw ! OUTPUT
    ! ----------------------------------------------------------------------------

    ecalc_tro_vdw= rg*T/gam1 - avw*ro
    
  end function ecalc_tro_vdw

  !===============================================================================
  function dedTcalc_tro_vdw(T,ro)
  !===============================================================================
    !> Compute internal energy from T and rho
    !> - VDW EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: dedTcalc_tro_vdw ! OUTPUT
    ! ----------------------------------------------------------------------------

    dedTcalc_tro_vdw= rg/gam1

  end function dedTcalc_tro_vdw

  !===============================================================================
  function gcalc_tro_vdw(T,ro)
  !===============================================================================
    !> Compute fundamental derivative of gas dynamics from T and ro
    !> - VDW EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: gcalc_tro_vdw ! OUTPUT
    ! ----------------------------------------------------------------------------

    gcalc_tro_vdw= 0.5_wp*(gam*(gam+1.0_wp)*rg*T/(1.0_wp-bvw*ro)**3-6.0_wp*avw*ro) &
         / (gam*rg*T/(1.0_wp-bvw*ro)**2-2.0_wp*avw*ro)
    
  end function gcalc_tro_vdw

  !===============================================================================
  function pcalc_roero_vdw(roe,ro,Ttent)
  !===============================================================================
    !> Compute pressure from rhoe and rho
    !> (the tentative temperature is not used for vdw)
    !> - VDW EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: roe,ro,Ttent
    real(wp) :: pcalc_roero_vdw ! OUTPUT
    ! ----------------------------------------------------------------------------

    pcalc_roero_vdw= gam1*(roe+ro**2*avw)/(1.0_wp-ro*bvw)-ro**2*avw
    
  end function pcalc_roero_vdw

  !===============================================================================
  function pcalc_tro_vdw(T,ro)
  !===============================================================================
    !> Compute pressure from T and rho
    !> - VDW EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: pcalc_tro_vdw ! OUTPUT
    ! ----------------------------------------------------------------------------

    pcalc_tro_vdw= ro*rg*T/(1.0_wp-ro*bvw)-ro**2*avw
    
  end function pcalc_tro_vdw

  !===============================================================================
  function rocalc_ep_vdw(e,p,Ttent)
  !===============================================================================
    !> Compute density from internal energy e and p
    !> (the tentative temperature is not used for vdw)
    !> - VDW EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: e,p,Ttent
    real(wp) :: rocalc_ep_vdw ! OUTPUT
    ! ----------------------------------------------------------------------------

    rocalc_ep_vdw= 1.0_wp
    
  end function rocalc_ep_vdw

  !===============================================================================
  function rocalc_ps_vdw(p,s,Ttent)
  !===============================================================================
    !> Compute density from p and entropy s
    !> (the tentative temperature is not used for vdw)
    !> - VDW EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: s,p,Ttent
    real(wp) :: rocalc_ps_vdw ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: ro,ro1,func,der_func,err
    ! ----------------------------------------------------------------------------

    ro = 0.5_wp*roc

    do i=1,100
       func    = (p + ro**2*avw)*(1.0_wp/ro-bvw)**gam - s
       der_func= (1.0_wp/ro-bvw)**gam1 * (-gam*p/ro**2 + (2.0_wp-gam)*avw - 2.0_wp*ro*avw*bvw)

       ro1 = ro - func/der_func
       err = abs(ro1-ro)/ro1
       if (abs(err).le.tol) then
          rocalc_ps_vdw= ro1
          return
       endif
       ro = ro1
    enddo

    print *,'in rocalc_ps: at iteration',i,'error=',err
    call mpistop('function rocalc_ps not converged',0)
    
  end function rocalc_ps_vdw

  !===============================================================================
  function rocalc_pt_vdw(p,T,rotent)
  !===============================================================================
    !> Compute density from p and T
    !> - VDW EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: p,T,rotent
    real(wp) :: rocalc_pt_vdw ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: ro,ro1,func,der_func,err
    ! ----------------------------------------------------------------------------

    ro = rotent

    do i=1,200
       func    = ro*rg*T/(1.0_wp-ro*bvw) - ro**2 * avw - p
       der_func= rg*T/(1.0_wp-ro*bvw)**2 - 2.0_wp*ro*avw

       ro1 = ro - func/der_func
       err = (ro1-ro)/ro1
       if (abs(err).le.tol) then
          rocalc_pt_vdw= ro1
          return
       endif
       ro = ro1
    enddo

    print *,'in rocalc_pt: at iteration',i,'error=',err
    print *,'p:',p,'T:',T
    rocalc_pt_vdw= 0.0_wp
    call mpistop('function rocalc_pt not converged',0)
    
  end function rocalc_pt_vdw

  !===============================================================================
  function rocalc_st_vdw(s,T,rotent)
  !===============================================================================
    !> Compute density from entropy s, T
    !> (the tentative density is not used for vdw)
    !> - VDW EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: s,T,rotent
    real(wp) :: rocalc_st_vdw ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: ro,ro1,func,der_func,err
    ! ----------------------------------------------------------------------------

    ro= rotent

    do i=1,100
       func    = s - rg*T*(1.0_wp/ro-bvw)**gam1
       der_func= rg*T/ro**2*gam1*(1.0_wp/ro-bvw)**(gam-2.0_wp)

       ro1 = ro - func/der_func
       err = (ro1-ro)/ro1

       if (abs(err).le.tol) then
          rocalc_st_vdw= ro1
          return
       endif
       ro = ro1
    enddo

    print *,'in rocalc_st: at iteration',i,'error=',err
    call mpistop('function rocalc_st not converged',0)
    
  end function rocalc_st_vdw

  !===============================================================================
  function scalc_tro_vdw(T,ro)
  !===============================================================================
    !> Compute entropy s from T and rho
    !> - VDW EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: scalc_tro_vdw ! OUTPUT
    ! ----------------------------------------------------------------------------

    scalc_tro_vdw= rg*T*(1.0_wp/ro-bvw)**gam1
    
  end function scalc_tro_vdw

  !===============================================================================
  function tcalc_ph_vdw(p,h,rotent,Ttent)
  !===============================================================================
    !> Compute temperature from p and enthalpy h
    !> (the tentative density are not used for vdw)
    !> - VDW EOS - 
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: p,h,rotent,Ttent
    real(wp) :: tcalc_ph_vdw ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: j
    real(wp) :: T,T1,ro,e,err
    ! ----------------------------------------------------------------------------
 
    T1= Ttent
    do j=1,10000
       T =T1
       ro= rocalc_pt_vdw(p,T,rotent)

       e = h - p/ro
       T1 = tcalc_roero_vdw(ro*e,ro,T)
       err = abs(T1-T)/T

       if (err.le.tol) then
          tcalc_ph_vdw= T
          return
       endif
    enddo
    
  end function tcalc_ph_vdw

  !===============================================================================
  function tcalc_pro_vdw(p,ro,Ttent)
  !===============================================================================
    !> Compute temperature from p and rho
    !> (the tentative temperature are not used for vdw)
    !> - VDW EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: p,ro,Ttent
    real(wp) :: tcalc_pro_vdw ! OUTPUT
    ! ----------------------------------------------------------------------------

    tcalc_pro_vdw= (p+ro**2*avw)*(1.0_wp/ro-bvw)/rg
    
  end function tcalc_pro_vdw

  !===============================================================================
  function tcalc_roero_vdw(roe,ro,Ttent)
  !===============================================================================
    !> Compute temperature from rho*e and rho
    !> (the tentative temperature are not used for vdw)
    !> - VDW EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: roe,ro,Ttent
    real(wp) :: tcalc_roero_vdw ! OUTPUT
    ! ----------------------------------------------------------------------------

    tcalc_roero_vdw= gam1/rg*(roe/ro+ro*avw)
    
  end function tcalc_roero_vdw

  !===============================================================================
  function tcalc_sro_vdw(s,ro,Ttent)
  !===============================================================================
    !> Compute temperature from entropy s and rho
    !> (the tentative temperature are not used for vdw)
    !> - VDW EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: s,ro,Ttent
    real(wp) :: tcalc_sro_vdw ! OUTPUT
    ! ----------------------------------------------------------------------------

    tcalc_sro_vdw= s/(rg*(1.0_wp/ro-bvw)**gam1)
    
  end function tcalc_sro_vdw

  !===============================================================================
  subroutine tcalc_Hstot_vdw(T,ro,p,H_tot,s_tot,am0,coeff)
  !===============================================================================
    !> Compute temperature from total stagnation enthalpy (Newton iteration)
    !> used for imposing Riemann invariant in inflow BC
    !> - VDW EOS -
  !===============================================================================
    use warnstop
    implicit none
    ! ----------------------------------------------------------------------------
    ! input/output: initial and updated thermo. variables
    real(wp), intent(inout) :: T,ro,p
    ! input: total enthalpy and entropy to be conserved
     real(wp), intent(in) :: H_tot,s_tot
   ! input: interior part of Riemann invariant am0
    real(wp), intent(in) :: am0
    ! input: coeff taking into account flow direction
    real(wp), intent(in) :: coeff
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: fn,der_fn,err
    real(wp) :: T1,ro1,dpdT
    ! ----------------------------------------------------------------------------

    T1=T
    ro1=ro

    ! **********************
    ! v1 : not inlined
    ! **********************
    ! inlined version not written [cf mod_eos_prs.f90 for ex inlined version]

    ! Newton's algorithm
    ! ------------------
    do i=1,1000

       ! update density from total entropy
       ro=rocalc_st_vdw(s_tot,T1,ro1)

       ! update pressure
       p=pcalc_tro_vdw(T1,ro)

       ! update pressure derivative wrt T
       dpdT=dpdicalc_tro_vdw(T1,ro)*cvcalc_tro_vdw(T1,ro)

       ! update temperature

       fn=ecalc_tro_vdw(T1,ro)+p/ro+0.5_wp*coeff*(am0-p)**2-H_tot

       der_fn=dedTcalc_tro_vdw(T1,ro)+dpdT/ro-coeff*(am0-p)*dpdT

       T= T1-fn/der_fn

       err= abs(T1-T)/T

       if (err.le.tol*10.) then
          ro=rocalc_st_vdw(s_tot,T,ro1)
          p=pcalc_tro_vdw(T,ro)
          return
       endif

       ! store previous solution
       T1=T
       ro1=ro

    enddo

    print *,'in tcalc_Hstot: at iteration',i,'error=',err
    print *,'T :',T,'ro:',ro
    call mpistop('function tcalc_Hstot not converged',0)

  end subroutine tcalc_Hstot_vdw

  !===============================================================================
  function vvol_vdw(pi,Ti,vtent)
  !===============================================================================
    !> Compute vvol such that p(vvol)=0
    !> - VDW EOS -
  !===============================================================================
     implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: pi,Ti,vtent
    real(wp) :: vvol_vdw ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer  :: i
    real(wp) :: v0,v1,err,p,T
    ! ----------------------------------------------------------------------------

    T = Ti*tc
    p = pi*pc
    v0= vtent/roc

    do i=1,15000
       v1=v0-(pcalc_tro_vdw(T,1.0_wp/v0)-p)/dpdvcalc_tro_vdw(T,1.0_wp/v0)

       err = abs(v1-v0)/v0
       if (err.le.tol) then
          vvol_vdw= v1*roc
          return
       endif

       v0=v1
    enddo

    print *,'in vvol: at iteration',i,'error=',err
    call mpistop('function vvol not converged',0)
    
  end function vvol_vdw

  !===============================================================================
  function vvol_d1_vdw(Ti,vtent)
  !===============================================================================
    !> Compute vvol_d1 from Ti and tentative v
    !> - VDW EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: Ti,vtent
    real(wp) :: vvol_d1_vdw ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: v0,v1,err,T
    ! ----------------------------------------------------------------------------

    T = Ti*tc
    v0 = vtent/roc

    do i=1,15000
       v1 = v0 - dpdvcalc_tro_vdw(T,1.0_wp/v0)/d2pdv2calc_tro_vdw(T,1.0_wp/v0)

       err = abs(v1-v0)/v0
       if (err.le.tol) then
          vvol_d1_vdw= v1*roc
          return
       endif

       v0=v1
    enddo

    print *,'in vvol_d1: at iteration',i,'error=',err
    call mpistop('function vvol_d1 not converged',0)

  end function vvol_d1_vdw

  !===============================================================================
  function vvol_d2_vdw(Ti,vtent)
  !===============================================================================
    !> Compute vvol_d2 from Ti and tentative v
    !> - VDW EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: Ti,vtent
    real(wp) :: vvol_d2_vdw ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: v0,v1,err,T
    ! ----------------------------------------------------------------------------

    T= Ti*tc
    v0= vtent/roc

    do i=1,15000
       v1 = v0 - d2pdv2calc_tro_vdw(T,1.0_wp/v0)/d3pdv3calc_tro_vdw(T,1.0_wp/v0)

       err = abs(v1-v0)/v0
       if (err.le.tol) then
          vvol_d2_vdw= v1*roc
          return
       endif

       v0=v1
    enddo

    print *,'in vvol_d2: at iteration',i,'error=',err
    call mpistop('function vvol_d2 not converged',0)

  end function vvol_d2_vdw

  !===============================================================================
  function intpcalc_tro_vdw(Ti,roi)
  !===============================================================================
    !> Compute pressure integral from T and rho
    !> - VDW EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: Ti,roi
    real(wp) :: intpcalc_tro_vdw ! OUTPUT
    ! ----------------------------------------------------------------------------
    real(wp) :: v,ro,T
    ! ----------------------------------------------------------------------------

    T = Ti*tc
    ro= roi*roc
    v = 1.0_wp/ro

    intpcalc_tro_vdw= rg*T*log(v-bvw) + avw/v

    intpcalc_tro_vdw= intpcalc_tro_vdw/pc*roc

  end function intpcalc_tro_vdw

  !===============================================================================
  subroutine dpcalc_tro_vdw(T,ro,dpdro,dpdT)
  !===============================================================================
    !> Compute pressure derivatives w.r.t. density and temperature using T and rho
    !> - VDW EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in)  :: T,ro
    real(wp), intent(out) :: dpdro,dpdT ! OUTPUTS
    ! ----------------------------------------------------------------------------

    ! Derivative w.r.t. density
    dpdro = rg*T /(1.0_wp-ro*bvw)**2 - 2.0_wp*ro*avw
    
    ! Derivative w.r.t. temperature
    dpdT  = rg*ro/(1.0_wp-ro*bvw)
    
  end subroutine dpcalc_tro_vdw

  !===============================================================================
  function p0calc_vdw(T,p,ro,M)
  !===============================================================================
    !> Compute total pressure from static quantities
    !> - VDW EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,M,ro,p
    real(wp) :: p0calc_vdw ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: j
    real(wp) :: s0,ro0,e0,p0,h0 ! total quantities
    real(wp) :: err,c2,e,T0,p1,h
    ! ----------------------------------------------------------------------------

    c2 = c2calc_tro_vdw(T,ro)
    e = ecalc_tro_vdw(T,ro)
    h = e + p/ro

    h0 = h + (M**2*c2)/2.0_wp
    s0 = scalc_tro_vdw(T,ro)

    p1 = p
    T0 = T
    ro0=ro
    
    err = 1.0_wp
    do j=1,10000
       p0 = p1
       T0 = tcalc_ph_vdw(p0,h0,ro0,T0)
       ro0 = rocalc_st_vdw(s0,T0,ro0)
       e0 = h0 - p0/ro0
       p1 = pcalc_roero_vdw(ro0*e0,ro0,T0)
       err = abs(p1-p0)/p0

       if (err.le.tol) then
          p0calc_vdw = p0
          return
       endif
    enddo

  end function p0calc_vdw

!!$  !===============================================================================
!!$  function tscalc_vdw(T0,p0,ro0,M,Ttent)
!!$  !===============================================================================
!!$    !> Compute static temperature from total quantities
!!$    !> - VDW EOS -
!!$  !===============================================================================
!!$    implicit none
!!$    ! ----------------------------------------------------------------------------
!!$    real(wp), intent(in) :: T0,M,ro0,p0,Ttent
!!$    real(wp) :: tscalc_vdw ! OUTPUT
!!$    ! ----------------------------------------------------------------------------
!!$    integer :: j
!!$    real(wp) :: c2,s0,ro,e,p,h,h0,err,T1,T
!!$    ! ----------------------------------------------------------------------------
!!$
!!$    h0 = ecalc_tro(T0,ro0) + p0/ro0
!!$    s0 = scalc_tro(T0,ro0)
!!$
!!$    T1 = Ttent*0.9_wp
!!$    do j=1,100
!!$       T = T1
!!$       ro = rocalc_st_vdw(s0,T,ro0)
!!$       !e = ecalc_tro_vdw(T,ro)
!!$       p = pcalc_tro_vdw(T,ro)
!!$       c2 = c2calc_tro_vdw(T,ro)
!!$
!!$       !T1 = tcalc_sro_vdw(s0,ro,T)
!!$       !T1 = tcalc_pro_vdw(p,ro,T)
!!$       !T1 = tcalc_roero_vdw(ro*e,ro,T)
!!$       h = h0 - (M**2*c2)/2
!!$       T1 = tcalc_ph_vdw(p,h,ro,T)
!!$
!!$       err = abs(T1-T)/T
!!$       write(*,*) j
!!$       if (err.le.tol) then
!!$          !h = h0 - (M**2*c2)/2
!!$          tscalc_vdw= tcalc_ph_vdw(p,h,ro,T)
!!$          return
!!$       endif
!!$    enddo
!!$
!!$  end function tscalc_vdw
!!$
!!$  !===============================================================================
!!$  function pscalc_vdw(T0,p0,Ts)
!!$  !===============================================================================
!!$    !> Compute static pressure from total quantities
!!$    !> - VDW EOS -
!!$  !===============================================================================
!!$    implicit none
!!$    ! ----------------------------------------------------------------------------
!!$    real(wp), intent(in) :: p0,T0,Ts
!!$    real(wp) :: pscalc_vdw ! OUTPUT
!!$    ! ----------------------------------------------------------------------------
!!$
!!$    pscalc_vdw = p0*(Ts/T0)**(gam/gam1)
!!$
!!$  end function pscalc_vdw
!!$
end module mod_eos_vdw
