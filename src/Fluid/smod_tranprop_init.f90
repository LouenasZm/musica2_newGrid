!=================================================================================
submodule (mod_tranprop) smod_tranprop_init
!=================================================================================
  !> author: XG
  !> date: January 2024
  !> Submodule to init/assign transport properties
!=================================================================================

contains

  !===============================================================================
  module subroutine init_viscosity
  !===============================================================================
    !> Initializations for viscosity computation
    !===============================================================================
    use mod_mpi
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i
    real(wp) :: n_OH,kv,mdmr,fct
    !!real(wp) :: sigma_c
    !!real(wp), parameter :: NA=6.022140857e23_wp ! Avogadro's constant
    ! ----------------------------------------------------------------------------
!!$    character(len=255) :: hf(1),hmix,herr
!!$    integer :: ierr
!!$    character(len=255) :: herr
!!$    real(wp) :: T_in,ro_in

    if (eos_type.eq.'pfg'.or.eos_type.eq.'vdw') then
       
       ! Power-Law viscosity parameters for PFG equation
       ! -----------------------------------------------
       if( fluidname.eq.'air') then
          ! T0 = 65.0_wp
          ! mu0= 4.3559e-6_wp
          T0  = 273.0_wp
          mu0 = 1.716e-5_wp
          S0 = 110.4_wp
          var1 = 1.0_wp + S0/T0
          c0 = 1.4579326545176254e-06_wp
       elseif (fluidname.eq.'pp11') then
          T0 = Tc
          mu0= 2.5e-5_wp
       else
          call mpistop('no ref viscosity values for this fluid!',0)
       endif

       var2 = cpfg/Pr           ! cp/Pr = lambda/mu
       var3 = mu0/T0**nexp

       ! Define viscosity law
       ! --------------------
       if (visc_type.eq.'S') then
          viscosity_law => viscosity_law_sutherland
       elseif (visc_type.eq.'P') then
          viscosity_law => viscosity_law_power
       else
          call mpistop('Problem in init_viscosity',0)
       endif
       
       ! Define thermal conductivity law
       ! -------------------------------
       thconductivity => thconductivity_Prcost
       
    elseif (eos_type.eq.'ref') then
       
       viscosity_law => viscosity_law_refprop
       thconductivity => thconductivity_refprop
       
    else
       
       ! Choose original roc or roc modified by cubic law
       ! ------------------------------------------------
       ! MANUAL SWITCH
       !if (eos_type.eq.'prs') then
       !   roc0=roc
       !endif
       !print *,roc0,roc
       
       ! This loop computes the array needed for Chung-Lee-Starling law init
       ! -------------------------------------------------------------------
       ! To be modified, n_OH is the number of -OH groups in the fluid formula
       n_OH=0.0_wp
       
       ! computation of critical volume in cm^3/mol
       ! ------------------------------------------
       vmolc=1.e3_wp*pmol/roc0
       
       ! association factor
       ! ------------------
       kv=0.0_wp ! 0.0682.0_wp + 4.704.0_wp * n_OH/pmol
       
       ! reduced dipole moment
       ! ---------------------
       mdmr=131.3_wp*mdm/sqrt(Tc*vmolc)
       
       ! Fc coefficient
       ! --------------
       fct=1.0_wp-2.756e-1_wp*om+5.9035e-2_wp*mdmr**4+kv

       k_chung=4.0785e-5_wp*sqrt(pmol)/vmolc**(2.0_wp/3.0_wp)*fct

       do i=1,10
          AA1(i)=av(i)+bv(i)*om+cv(i)*mdmr**4+dv(i)*kv
       enddo
       do i=1,7
          BB1(i)=atc(i)+btc(i)*om+ctc(i)*mdmr**4+dtc(i)*kv
       enddo

       ! Coefficient for first-density contribution in Chen et al. law
       ! -------------------------------------------------------------
       ! Lennard-Jones collision diameter (converted nm -> m)
       !!sigma_c=0.809*vmolc**(1.0_wp/3.0_wp)*1.e-9_wp
       ! reducing bw coefficients
       !!bw=bw*NA*sigma_c**3/pmol
       bw=bw*0.6022140857_wp*0.809**3*vmolc*1.e-3_wp/pmol

       if (iproc==0) print *,'reducing bw coefficient:',0.6022140857_wp*0.809**3*vmolc*1.e-3_wp
              
       if (visc_type.eq.'C') then
          viscosity_law  => viscosity_law_chung
          thconductivity => thconductivity_chung
       elseif (visc_type.eq.'W') then
          viscosity_law  => viscosity_law_wen
          !thconductivity => thconductivity_wen
          thconductivity => thconductivity_chung
       else
          call mpistop('Problem in init_viscosity',0)
       endif
    endif

!!$    !test T_in =373.15_wp; ro_in=48.51_wp !!REFPROP:1.830069407549051E-005
!!$    T_in =373.15_wp
!!$    ro_in=48.51_wp
!!$    T_in =350.0_wp
!!$    ro_in=4.42_wp
!!$    print *,'Chung',viscosity_law_chung(T_in,ro_in)
!!$    print *,'Wen',viscosity_law_wen(T_in,ro_in)
!!$
!!$    ! RefProp setup
!!$    ! =============
!!$    call SETPATH('../../src/Fluid/refprop/fluids')
!!$    hf(1)=trim(fluidname)//'.fld'
!!$    hmix ='hmix.bnc'
!!$    call SETUP(1,hf,hmix,'DEF',ierr,herr)
!!$    !call INFOr(1,pmol,ttp,teb,tc,pc,roc,zc,acf,mdm,rg)
!!$
!!$    roc=roc*pmol
!!$    
!!$    print *,'Refprop',viscosity_law_refprop(T_in,ro_in)
!!$    
!!$    call mpistop('dev',0)
    
  end subroutine init_viscosity

end submodule smod_tranprop_init
