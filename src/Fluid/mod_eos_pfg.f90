!=================================================================================
module mod_eos_pfg
!=================================================================================
  !> Module to define subroutines for PERFECT GAS Equation of State (EoS)
!=================================================================================
  use mod_fluid
  implicit none
  ! ------------------------------------------------------------------------------
  real(wp), private, parameter :: tol=1.e-6_wp
  ! ------------------------------------------------------------------------------

contains
  
  !===============================================================================
  subroutine init_eos_pfg ! [called by musica_main]
  !===============================================================================
    !> Initializations of gas coefficients
    !> - PFG EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    ! read gam (ratio of specific heats)
    ! [read_eos in mod_fluid.f90]
    call read_eos
    
    gam1= gam-1.0_wp
    igm1= 1.0_wp/gam1

    ! specific heats (at constant pressure and volume)
    cpfg= gam/gam1*rg
    cvfg= 1.0_wp/gam1*rg
    
  end subroutine init_eos_pfg

  !===============================================================================
  function avcalc_tro_pfg(T,ro)
  !===============================================================================
    !> Compute isobaric expansion coefficient from T (and rho, not used)
    !> - PFG EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: avcalc_tro_pfg ! OUTPUT
    ! ----------------------------------------------------------------------------

    avcalc_tro_pfg= 1.0_wp/T
    
  end function avcalc_tro_pfg

  !===============================================================================
  function c2calc_tro_pfg(T,ro)
  !===============================================================================
    !> Compute speed of sound from T (and rho, not used)
    !> - PFG EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: c2calc_tro_pfg ! OUTPUT 
    ! ----------------------------------------------------------------------------

    c2calc_tro_pfg= gam*rg*T
    
  end function c2calc_tro_pfg

  !===============================================================================
  function cpcalc_tro_pfg(T,ro)
  !===============================================================================
    !> Compute heat capacity at constant pressure
    !> (a constant for a calorically perfect gas)
    !> - PFG EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: cpcalc_tro_pfg ! OUTPUT 
    ! ----------------------------------------------------------------------------
    
    cpcalc_tro_pfg= cpfg
    
  end function cpcalc_tro_pfg

  !===============================================================================
  function cvcalc_tro_pfg(T,ro)
  !===============================================================================
    !> Compute heat capacity at constant volume
    !> (a constant for a calorically perfect gas)
    !> - PFG EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: cvcalc_tro_pfg ! OUTPUT 
    ! ----------------------------------------------------------------------------

    cvcalc_tro_pfg= cvfg
    
  end function cvcalc_tro_pfg

  !===============================================================================
  function dpdicalc_tro_pfg(T,ro)
  !===============================================================================
    !> Compute pressure derivative w.r.t temperature
    !> - PFG EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: dpdicalc_tro_pfg ! OUTPUT 
    ! ----------------------------------------------------------------------------
    
    dpdicalc_tro_pfg= gam1*ro
    
  end function dpdicalc_tro_pfg

  !===============================================================================
  function dpdTcalc_tro_pfg(T,ro)
  !===============================================================================
    !> Compute pressure derivative w.r.t temperature
    !> - PFG EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: dpdTcalc_tro_pfg ! OUTPUT
    ! ----------------------------------------------------------------------------

    dpdTcalc_tro_pfg = ro*rg

  end function dpdTcalc_tro_pfg

  !===============================================================================
  function dpdvcalc_tro_pfg(T,ro)
  !===============================================================================
    !> Compute pressure derivative w.r.t volume
    !> - PFG EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: dpdvcalc_tro_pfg ! OUTPUT
    ! ----------------------------------------------------------------------------

    dpdvcalc_tro_pfg = -ro**2*rg*T

  end function dpdvcalc_tro_pfg

  !===============================================================================
  function ecalc_pro_pfg(p,ro,Ttent)
  !===============================================================================
    !> Compute internal energy from p and rho
    !> (the tentative temperature is not used for pfg)
    !> - PFG EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: p,ro ! INPUT pressure and density
    real(wp), intent(in), optional :: Ttent ! first guess for more complex EOS
    real(wp) :: ecalc_pro_pfg ! OUTPUT computed internal energy
    ! ----------------------------------------------------------------------------

    ecalc_pro_pfg= igm1*p/ro

  end function ecalc_pro_pfg

  !===============================================================================
  function ecalc_tro_pfg(T,ro)
  !===============================================================================
    !> Compute internal energy from T (and rho, not used)
    !> 1st Joule's law for pfg
    !> - PFG EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: ecalc_tro_pfg ! OUTPUT 
    ! ----------------------------------------------------------------------------
    
    ecalc_tro_pfg= cvfg*T
    
  end function ecalc_tro_pfg

  !===============================================================================
  function dedTcalc_tro_pfg(T,ro)
  !===============================================================================
    !> Compute internal energy derivative w.r.t rho
    !> - PFG EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: dedTcalc_tro_pfg ! OUTPUT
    ! ----------------------------------------------------------------------------

    dedTcalc_tro_pfg= cvfg

  end function dedTcalc_tro_pfg

  !===============================================================================
  function dedrocalc_tro_pfg(T,ro)
  !===============================================================================
    !> Compute internal energy derivative w.r.t rho
    !> - PFG EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: dedrocalc_tro_pfg ! OUTPUT
    ! ----------------------------------------------------------------------------

    dedrocalc_tro_pfg= -cvfg*T/ro

  end function dedrocalc_tro_pfg

  !===============================================================================
  function gcalc_tro_pfg(T,ro)
  !===============================================================================
    !> Compute fundamental derivative of gas dynamics (from T and ro, not used)
    !> constant (gam+1)/2 for pfg
    !> - PFG EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: gcalc_tro_pfg ! OUTPUT
    ! ----------------------------------------------------------------------------
    
    gcalc_tro_pfg= 0.5_wp*(gam + 1.0_wp)
    
  end function gcalc_tro_pfg

  !===============================================================================
  function pcalc_roero_pfg(roe,ro,Ttent)
  !===============================================================================
    !> Compute pressure from rhoe (and rho, not used)
    !> (the tentative temperature is not used for pfg)
    !> - PFG EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: roe,ro,Ttent
    real(wp) :: pcalc_roero_pfg ! OUTPUT
    ! ----------------------------------------------------------------------------
    
    pcalc_roero_pfg= gam1*roe
    
  end function pcalc_roero_pfg

  !===============================================================================
  function pcalc_tro_pfg(T,ro)
  !===============================================================================
    !> Compute pressure from T and rho
    !> - PFG EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,ro
    real(wp) :: pcalc_tro_pfg ! OUTPUT
    ! ----------------------------------------------------------------------------
    
    pcalc_tro_pfg= ro*rg*T
    
  end function pcalc_tro_pfg

  !===============================================================================
  function rocalc_ep_pfg(e,p,Ttent)
  !===============================================================================
    !> Compute density from internal energy e and p
    !> (the tentative temperature is not used for pfg)
    !> - PFG EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: e,p,Ttent
    real(wp) :: rocalc_ep_pfg ! OUTPUT
    ! ----------------------------------------------------------------------------
    
    rocalc_ep_pfg= p/(e*gam1)
    
  end function rocalc_ep_pfg

  !===============================================================================
  function rocalc_ps_pfg(p,s,Ttent)
  !===============================================================================
    !> Compute density from p and entropy s
    !> (the tentative temperature is not used for pfg)
    !> - PFG EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: p,s,Ttent
    real(wp) :: rocalc_ps_pfg ! OUTPUT
    ! ----------------------------------------------------------------------------
    
    rocalc_ps_pfg= (p/s)**(1.0_wp/gam)
    
  end function rocalc_ps_pfg

  !===============================================================================
  function rocalc_pt_pfg(p,T,rotent)
  !===============================================================================
    !> Compute density from p and T
    !> (the tentative density is not used for pfg)
    !> - PFG EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: p,T,rotent
    real(wp) :: rocalc_pt_pfg ! OUTPUT
    ! ----------------------------------------------------------------------------
    
    rocalc_pt_pfg= p/(rg*T)
    
  end function rocalc_pt_pfg

  !===============================================================================
  function rocalc_st_pfg(s,T,rotent)
  !===============================================================================
    !> Compute density from entropy s, T
    !> (the tentative density is not used for pfg)
    !> - PFG EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: s,T,rotent
    real(wp) :: rocalc_st_pfg ! OUTPUT
    ! ----------------------------------------------------------------------------
    
    rocalc_st_pfg= (rg*T/s)**(1.0_wp/gam1)
    
  end function rocalc_st_pfg

  !===============================================================================
  function scalc_tro_pfg(tr,ro) ! TO BE CHANGED Why tr ???
  !===============================================================================
    !> Compute entropy s from T and rho
    !> - PFG EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: tr,ro
    real(wp) :: scalc_tro_pfg ! OUTPUT
    ! ----------------------------------------------------------------------------

    scalc_tro_pfg= rg*tr/ro**gam1
    
  end function scalc_tro_pfg

  !===============================================================================
  function tcalc_ph_pfg(p,h,rotent,Ttent)
  !===============================================================================
    !> Compute temperature from p and enthalpy h
    !> (the tentative density & temperature are not used for pfg)
    !> only h is used (2nd Joule's law)
    !> - PFG EOS - 
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: p,h,rotent,Ttent
    real(wp) :: tcalc_ph_pfg ! OUTPUT
    ! ----------------------------------------------------------------------------
    
    tcalc_ph_pfg= h/cpfg
    
  end function tcalc_ph_pfg

  !===============================================================================
  function tcalc_pro_pfg(p,ro,Ttent)
  !===============================================================================
    !> Compute temperature from p and rho
    !> (the tentative temperature are not used for pfg)
    !> - PFG EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: p,ro,Ttent
    real(wp) :: tcalc_pro_pfg ! OUTPUT
    ! ----------------------------------------------------------------------------
    
    tcalc_pro_pfg= p/(rg*ro)
    
  end function tcalc_pro_pfg

  !===============================================================================
  function tcalc_roero_pfg(roe,ro,Ttent)
  !===============================================================================
    !> Compute temperature from rho*e and rho
    !> (the tentative temperature are not used for pfg)
    !> - PFG EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: roe,ro,Ttent
    real(wp) :: tcalc_roero_pfg ! OUTPUT
    ! ----------------------------------------------------------------------------
    
    tcalc_roero_pfg= roe/(ro*cvfg)
    
  end function tcalc_roero_pfg

  !===============================================================================
  function tcalc_sro_pfg(s,ro,Ttent)
  !===============================================================================
    !> Compute temperature from entropy s and rho
    !> (the tentative temperature are not used for pfg)
    !> - PFG EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: s,ro,Ttent
    real(wp) :: tcalc_sro_pfg ! OUTPUT
    ! ----------------------------------------------------------------------------
    
    tcalc_sro_pfg= s*ro**gam1/rg
    
  end function tcalc_sro_pfg

  !===============================================================================
  function tcalc_Hstot_pfg_(H_tot,s_tot,Ttent,rotent,coeff1,coeff2)
  !===============================================================================
    !> Compute temperature from total stagnation enthalpy (Newton iteration)
    !> used for imposing Riemann invariant in inflow BC
    !> - PFG EOS -
  !===============================================================================
    use warnstop
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: H_tot,s_tot,Ttent,rotent,coeff1,coeff2
    real(wp) :: tcalc_Hstot_pfg_ ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: T,ro,T1,p,dpdT,fn,der_fn,err!,ro1
    ! ----------------------------------------------------------------------------

    T1=Ttent
    !ro1=rotent not used for PFG

    ! Newton's algorithm
    ! ------------------
    do i=1,100

       ! update density from total entropy
       ro=(rg*T1/s_tot)**(1.0_wp/gam1)

       ! update pressure and its derivative wrt T
       p=ro*rg*T1
       dpdT=ro*rg

       ! update temperature
       fn= cvfg*T1+p/ro+0.5_wp*coeff1*(coeff2-p)**2-H_tot
       der_fn= cvfg+dpdT/ro-coeff1*(coeff2-p)*dpdT
       T= T1 - fn/der_fn

       err= abs(T1-T)/T

       if (err.le.tol) then
          tcalc_Hstot_pfg_= T
          return
       endif

       ! store previous solution
       T1=T
       !ro1=ro not used for PFG

    enddo

    print *,'in tcalc_Hstot: at iteration',i,'error=',err
    print *,'T :',T,'ro:',ro
    tcalc_Hstot_pfg_= 0.0_wp
    call mpistop('function tcalc_Hstot not converged',0)

  end function tcalc_Hstot_pfg_

  !===============================================================================
  subroutine tcalc_Hstot_pfg(T,ro,p,H_tot,s_tot,am0,coeff)
  !===============================================================================
    !> Compute temperature from total stagnation enthalpy (Newton iteration)
    !> used for imposing Riemann invariant in inflow BC
    !> - PFG EOS -
  !===============================================================================
    use warnstop
    use mod_mpi
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
    real(wp) :: T1,dpdT,fn,der_fn,err!,ro1
    ! ----------------------------------------------------------------------------

    T1=T*0.999
    !ro1=ro not used for PFG

    ! Newton's algorithm
    ! ------------------
    do i=1,1000

       ! update density from total entropy
       ro=(rg*T1/s_tot)**(1.0_wp/gam1)

       ! update pressure and its derivative wrt T
       p=ro*rg*T1
       dpdT=ro*rg

       ! update temperature
       fn= cvfg*T1+p/ro+0.5_wp*coeff*(am0-p)**2-H_tot
       der_fn= cvfg+dpdT/ro-coeff*(am0-p)*dpdT
       T= T1-fn/der_fn
       !print *,'dfn terms',cvfg,dpdT/ro,-coeff*(am0-p)*dpdT
       !print *,'fn dfn',fn,der_fn,T1,fn/der_fn
       !print *,'update T',T,ro,i

       err= abs(T1-T)/T

       if ((err.le.tol).or.(i>=1000)) then
          if ((i==1000).and.(iproc.eq.0)) print *,'i1000 -> found T',T,err
          !if (err<0) T=300.0_wp
          if (err<0) go to 10
          !print *,'found T',T
          ro=(rg*T/s_tot)**(1.0_wp/gam1)
          p=ro*rg*T
          return
       endif

       ! store previous solution
       T1=T
       !ro1=ro not used for PFG

    enddo

10  print *,'in tcalc_Hstot: at iteration',i,'error=',err
    print *,'T :',T,'ro:',ro
    !T=300.0_wp
    !ro=(rg*T/s_tot)**(1.0_wp/gam1)
    !p=ro*rg*T
    call mpistop('function tcalc_Hstot not converged',0)

  end subroutine tcalc_Hstot_pfg

  !===============================================================================
  function vvol_pfg(a,b,c)
  !===============================================================================
    !> Compute vvol such that p(vvol)=0
    !> vvol_=1 (a,b,c not used)
    !> - PFG EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp) :: a,b,c
    real(wp) :: vvol_pfg ! OUTPUT
    ! ----------------------------------------------------------------------------
    
    vvol_pfg= 1.0_wp
    
  end function vvol_pfg

  !===============================================================================
  function vvol_d1_pfg(a,b)
  !===============================================================================
    !> vvol_d1=1 (a,b not used)
    !> - PFG EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp) :: a,b
    real(wp) :: vvol_d1_pfg ! OUTPUT
    ! ----------------------------------------------------------------------------
    
    vvol_d1_pfg= 1.0_wp
    
  end function vvol_d1_pfg

  !===============================================================================
  function vvol_d2_pfg(a,b)
  !===============================================================================
    !> vvol_d2=1 (a,b not used)
    !> - PFG EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp) :: a,b
    real(wp) :: vvol_d2_pfg ! OUTPUT
    ! ----------------------------------------------------------------------------
    
    vvol_d2_pfg= 1.0_wp
    
  end function vvol_d2_pfg

  !===============================================================================
  function intpcalc_tro_pfg(a,b)
  !===============================================================================
    !> Compute pressure integral
    !> intpcalc_tro=1 (a,b not used)
    !> - PFG EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp) :: a,b
    real(wp) :: intpcalc_tro_pfg ! OUTPUT
    ! ----------------------------------------------------------------------------
    
    intpcalc_tro_pfg=1.0_wp
    
  end function intpcalc_tro_pfg

  !===============================================================================
  subroutine dpcalc_tro_pfg(T,ro,dpdro,dpdT)
  !===============================================================================
    !> Compute pressure derivatives w.r.t. density and temperature using T and rho
    !> - PFG EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in)  :: T,ro
    real(wp), intent(out) :: dpdro,dpdT ! OUTPUTS
    ! ----------------------------------------------------------------------------
    
    ! Derivative w.r.t. density
    dpdro= rg*T
    
    ! Derivative w.r.t. temperature
    dpdT = rg*ro
    
  end subroutine dpcalc_tro_pfg

  !===============================================================================
  function p0calc_pfg(T,p,ro,M)
  !===============================================================================
    !> Compute total pressure from static quantities
    !> - PFG EOS -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: T,M,ro,p
    real(wp) :: p0calc_pfg ! OUTPUT
    ! ----------------------------------------------------------------------------
    integer :: j
    real(wp) :: s0,ro0,e0,p0,h0 ! total quantities
    real(wp) :: err,c2,e,T0,p1,h
    ! ----------------------------------------------------------------------------

    c2= c2calc_tro_pfg(T,ro)
    e = ecalc_tro_pfg(T,ro)
    h = e+p/ro

    h0 = h+(M**2*c2)*0.5_wp
    s0 = scalc_tro_pfg(T,ro)

    p1 = p
    T0 = T
    ro0=ro
    
    err = 1.0_wp
    do j=1,10000
       p0 = p1
       T0 = tcalc_ph_pfg(p0,h0,ro0,T0)
       ro0= rocalc_st_pfg(s0,T0,ro0)
       e0 = h0 - p0/ro0
       p1 = pcalc_roero_pfg(ro0*e0,ro0,T0)
       err= abs(p1-p0)/p0

       if (err.le.tol) then
          p0calc_pfg = p0
          return
       endif
    enddo

  end function p0calc_pfg

  ! !===============================================================================
  ! subroutine stagnation_pfg(h0,s0,e0,T0,p0,ro0)
  ! !===============================================================================
  !   !> Compute stagnation quantities from static quantities & h0 & s0
  !   !> - PFG EOS -
  ! !===============================================================================
  !   implicit none
  !   ! ----------------------------------------------------------------------------
  !   real(wp), intent(in) :: h0,s0
  !   real(wp), intent(inout) :: e0,T0,p0,ro0
  !   ! ----------------------------------------------------------------------------
  !   integer :: j
  !   real(wp) :: err,p1
  !   ! ----------------------------------------------------------------------------

  !   p1 = p0
  !   err = 1.0_wp
  !   do j=1,100
  !      p0 = p1
  !      T0 = tcalc_ph_pfg(p0,h0,ro0,T0)
  !      ro0 = rocalc_st_pfg(s0,T0,ro0)
  !      e0 = h0 - p0/ro0
  !      p1 = pcalc_roero_pfg(ro0*e0,ro0,T0)
  !      err = abs(p1-p0)/p0

  !      if (err.le.tol) return
  !   enddo

  ! end subroutine stagnation_pfg

  !===============================================================================
  subroutine stagnation_pfg(h0,s,e0,T0,p0,ro0)
  !===============================================================================
    !> Compute stagnation quantities from h0 & s
    !> - PRS EOS -
  !===============================================================================
    use warnstop
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: h0,s
    real(wp), intent(inout) :: T0,ro0
    real(wp), intent(out) :: e0,p0
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: T1,fn,der_fn,dpdT,err
    ! ----------------------------------------------------------------------------

    T1 = T0

    ! Newton's algorithm
    ! ------------------
    do i=1,1000
       ! update density from total entropy
       ro0=rocalc_st_pfg(s,T1,ro0)

       ! update pressure
       p0=pcalc_tro_pfg(T1,ro0)

       ! update pressure derivative wrt T0
       dpdT=dpdicalc_tro_pfg(T1,ro0)*cvcalc_tro_pfg(T1,ro0)

       ! update temperature
       fn = ecalc_tro_pfg(T1,ro0) + p0/ro0 - h0
       der_fn = cvfg + dpdT/ro0
       T0= T1-fn/der_fn

       err= abs(T1-T0)/T0

       if (err.le.tol) then
          ro0=rocalc_st_pfg(s,T0,ro0)
          p0=pcalc_tro_pfg(T0,ro0)
          e0 = h0 - p0/ro0
          return
       endif

       ! store previous solution
       T1=T0

    enddo

    print *,'in stagnation_pfg: at iteration',i,'error=',err
    print *,'s :',s,'h0 :',h0,'T0:',T0,'p0:',ro0
    call mpistop('subroutine stagnation_pfg not converged',0)

  end subroutine stagnation_pfg


!!$  !===============================================================================
!!$  function tscalc_pfg(T0,p0,ro0,M,Ttent)
!!$  !===============================================================================
!!$    !> Compute static temperature from total quantities
!!$    !> - PFG EOS -
!!$  !===============================================================================
!!$    implicit none
!!$    ! ----------------------------------------------------------------------------
!!$    real(wp), intent(in) :: T0,M,p0,ro0,Ttent
!!$    real (wp) :: tscalc_pfg ! OUTPUT
!!$    ! ----------------------------------------------------------------------------
!!$    
!!$    tscalc = T0/(1.0_wp+0.5_wp*gam1*M**2)
!!$    
!!$  end  function tscalc
!!$
!!$  !===============================================================================
!!$  function pscalc_pfg(T0,p0,Ts)
!!$  !===============================================================================
!!$    !> Compute static pressure from total quantities
!!$    !> - PFG EOS -
!!$  !===============================================================================
!!$    implicit none
!!$    ! ----------------------------------------------------------------------------
!!$    real(wp), intent(in) :: p0,T0,Ts
!!$    real(wp) :: pscalc_pfg ! OUTPUT
!!$    ! ----------------------------------------------------------------------------
!!$    
!!$    pscalc = p0*(Ts/T0)**(gam/gam1)
!!$    !pscalc = p/(1.0_wp + 0.5_wp*gam1*M**2)**(gam/gam1)
!!$    
!!$  end function pscalc
    
end module mod_eos_pfg
