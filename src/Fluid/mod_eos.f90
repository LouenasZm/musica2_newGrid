!=================================================================================
module mod_eos
!=================================================================================
  !> Module to define pointer procedures for the chosen Equation of State (EoS)
!=================================================================================
  use mod_fluid
  use mod_constant
  use mod_primitives ! for: primitives subroutines
  use warnstop
  implicit none
  ! ------------------------------------------------------------------------------
  ! Interfaces for procedure pointer
  ! ------------------------------------------------------------------------------
  abstract interface
     real(wp) function thermo_typ1(var1,var2)
       use precision
       implicit none
       real(wp), intent(in) :: var1,var2
     end function thermo_typ1
  end interface
  ! ------------------------------------------------------------------------------
  abstract interface
     real(wp) function thermo_typ2(var1,var2,var3)
       use precision
       implicit none
       real(wp), intent(in) :: var1,var2,var3
     end function thermo_typ2
  end interface
  ! ------------------------------------------------------------------------------
  abstract interface
     real(wp) function thermo_typ3(var1,var2,var3,var4)
       use precision
       implicit none
       real(wp), intent(in) :: var1,var2,var3,var4
     end function thermo_typ3
  end interface
  ! ------------------------------------------------------------------------------
  abstract interface
     real(wp) function thermo_typ4(var1,var2,var3,var4,var5,var6)
       use precision
       implicit none
       real(wp), intent(in) :: var1,var2,var3,var4,var5,var6
     end function thermo_typ4
  end interface
  ! ------------------------------------------------------------------------------
  abstract interface
     subroutine thermo_typ5(var1,var2,var3,var4,var5,var6,var7)
       use precision
       implicit none
       real(wp), intent(inout) :: var1,var2,var3
       real(wp), intent(in) :: var4,var5,var6,var7
     end subroutine thermo_typ5
  end interface
  ! ------------------------------------------------------------------------------
  abstract interface
     subroutine thermo_typ6(var1,var2,var3,var4,var5,var6)
       use precision
       implicit none
       real(wp), intent(in) :: var1,var2
       real(wp), intent(inout) :: var4,var6
       real(wp), intent(out) :: var3,var5
     end subroutine thermo_typ6
  end interface
  ! ------------------------------------------------------------------------------
  abstract interface
     subroutine thermo_noarg
       implicit none
     end subroutine thermo_noarg
  end interface
  ! ------------------------------------------------------------------------------
  procedure(thermo_typ1), pointer ::  avcalc_tro => null() &
                                   ,  c2calc_tro => null() &
                                   ,  cpcalc_tro => null() &
                                   ,  cvcalc_tro => null() &
                                   ,   ecalc_tro => null() &
                                   ,   gcalc_tro => null() &
                                   ,   pcalc_tro => null() &
                                   ,   scalc_tro => null() &
                                   ,dpdicalc_tro => null() &
                                   ,   vvol_d1   => null() &
                                   ,   vvol_d2   => null() &
                                   ,intpcalc_tro => null() &
                                   ,dedrocalc_tro=> null() &
                                   ,dedTcalc_tro => null() &
                                   ,dpdvcalc_tro => null() &
                                   ,dpdTcalc_tro => null()
  procedure(thermo_typ2), pointer ::   ecalc_pro => null() &
                                   , pcalc_roero => null() &
                                   ,   rocalc_ep => null() &
                                   ,   rocalc_ps => null() &
                                   ,   rocalc_pt => null() &
                                   ,   rocalc_st => null() &
                                   ,   tcalc_sro => null() &
                                   ,   tcalc_pro => null() &
                                   , tcalc_roero => null() &
                                   ,     vvol    => null()

  procedure(thermo_typ3), pointer ::   tcalc_ph  => null()
  
  procedure(thermo_typ4), pointer ::   tcalc_Hstot_=> null()
  
  procedure(thermo_typ5), pointer ::   tcalc_Hstot=> null()

  procedure(thermo_typ6), pointer ::   stagnation_calc => null()

  procedure(thermo_noarg), pointer :: primitives     => null(), &
                                      primitives_visc=> null()
  
  procedure(thermo_noarg), pointer :: primitives_rans=> null()

  ! ------------------------------------------------------------------------------
  external :: primitives_SA
  ! ------------------------------------------------------------------------------

contains

  !===============================================================================
  subroutine init_eos ! [called by musica_main]
  !===============================================================================
    !> Initializations of EoS & assign procedure pointers for EoS functions
  !===============================================================================
    use mod_eos_pfg
    use mod_eos_vdw
    use mod_eos_mah
    use mod_eos_prs
    use mod_eos_swn
    use mod_eos_swp
    use mod_eos_ref
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------
    
    if (eos_type.eq.'pfg') then
       call init_eos_pfg
    elseif (eos_type.eq.'vdw') then
       call init_eos_vdw
    elseif (eos_type.eq.'mah') then
       call init_eos_mah
    elseif (eos_type.eq.'prs') then
       call init_eos_prs
    elseif (eos_type.eq.'swn') then
       call init_eos_swn
    elseif (eos_type.eq.'swp') then
       call init_eos_swp
    elseif (eos_type.eq.'ref') then
       call init_eos_ref
    endif
    
    call assign_procedures_eos
    
  end subroutine init_eos
  
  !===============================================================================
  subroutine assign_procedures_eos
  !===============================================================================
    !> Initializations of EoS & assign procedure pointers for EoS functions
  !===============================================================================
    use mod_eos_pfg
    use mod_eos_vdw
    use mod_eos_mah
    use mod_eos_prs
    use mod_eos_swn
    use mod_eos_swp
    use mod_eos_ref
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------
    
    if (eos_type.eq.'pfg') then
       
       ! primitive variables from conservative ones
       primitives => primitives_pfg
       primitives_visc => primitives_visc_pfg

       avcalc_tro => avcalc_tro_pfg
       c2calc_tro => c2calc_tro_pfg
       cpcalc_tro => cpcalc_tro_pfg
       cvcalc_tro => cvcalc_tro_pfg
       ecalc_tro => ecalc_tro_pfg
       gcalc_tro => gcalc_tro_pfg
       pcalc_tro => pcalc_tro_pfg
       scalc_tro => scalc_tro_pfg
       ecalc_pro => ecalc_pro_pfg
       pcalc_roero => pcalc_roero_pfg
       rocalc_pt => rocalc_pt_pfg
       rocalc_st => rocalc_st_pfg
       tcalc_ph  => tcalc_ph_pfg

       ! used for inflow BC (H_tot, s_tot, flow angle imposed)
       tcalc_Hstot_ => tcalc_Hstot_pfg_
       tcalc_Hstot => tcalc_Hstot_pfg

       ! NOT USED yet
       dpdicalc_tro => dpdicalc_tro_pfg
       dpdvcalc_tro => dpdvcalc_tro_pfg
       dpdTcalc_tro => dpdTcalc_tro_pfg
       dedrocalc_tro => dedrocalc_tro_pfg
       dedTcalc_tro => dedTcalc_tro_pfg
       rocalc_ep => rocalc_ep_pfg
       rocalc_ps => rocalc_ps_pfg
       tcalc_sro => tcalc_sro_pfg
       tcalc_pro => tcalc_pro_pfg
       tcalc_roero => tcalc_roero_pfg   

       ! functions used to compute saturation curves
       !vvol    => vvol_pfg
       !vvol_d1 => vvol_d1_pfg
       !vvol_d2 => vvol_d2_pfg
       !intpcalc_tro => intpcalc_tro_pfg

       ! function used to compute total quantities from static quantities
       stagnation_calc => stagnation_pfg

    elseif (eos_type.eq.'vdw') then
       
       ! primitive variables from conservative ones
       primitives => primitives_vdw
       primitives_visc => primitives_visc_vdw
       
       avcalc_tro => avcalc_tro_vdw
       c2calc_tro => c2calc_tro_vdw
       cpcalc_tro => cpcalc_tro_vdw
       cvcalc_tro => cvcalc_tro_vdw
       ecalc_tro => ecalc_tro_vdw
       gcalc_tro => gcalc_tro_vdw
       pcalc_tro => pcalc_tro_vdw
       scalc_tro => scalc_tro_vdw
       ecalc_pro => ecalc_pro_vdw
       pcalc_roero => pcalc_roero_vdw
       rocalc_pt => rocalc_pt_vdw
       rocalc_st => rocalc_st_vdw
       tcalc_ph  => tcalc_ph_vdw
       dedTcalc_tro => dedTcalc_tro_vdw
       dpdTcalc_tro => dpdTcalc_tro_vdw

       ! used for inflow BC (H_tot, s_tot, flow angle imposed)
       tcalc_Hstot => tcalc_Hstot_vdw

       ! NOT USED yet
       dpdicalc_tro => dpdicalc_tro_vdw
       rocalc_ep => rocalc_ep_vdw
       rocalc_ps => rocalc_ps_vdw
       tcalc_sro => tcalc_sro_vdw
       tcalc_pro => tcalc_pro_vdw
       tcalc_roero => tcalc_roero_vdw   

       ! functions used to compute saturation curves
       vvol    => vvol_vdw
       vvol_d1 => vvol_d1_vdw
       vvol_d2 => vvol_d2_vdw
       intpcalc_tro => intpcalc_tro_vdw

       ! function used to compute total quantities from static quantities
       ! stagnation_calc => stagnation_vdw
       !call mpistop('Subroutine stagnation_vdw not implemented yet !',0)

    elseif (eos_type.eq.'mah') then
       
       ! primitive variables from conservative ones
       primitives => primitives_mah
       primitives_visc => primitives_visc_mah
       
       avcalc_tro => avcalc_tro_mah
       c2calc_tro => c2calc_tro_mah
       cpcalc_tro => cpcalc_tro_mah
       cvcalc_tro => cvcalc_tro_mah
       ecalc_tro => ecalc_tro_mah
       gcalc_tro => gcalc_tro_mah
       pcalc_tro => pcalc_tro_mah
       scalc_tro => scalc_tro_mah
       ecalc_pro => ecalc_pro_mah
       pcalc_roero => pcalc_roero_mah
       rocalc_pt => rocalc_pt_mah
       rocalc_st => rocalc_st_mah
       tcalc_ph  => tcalc_ph_mah
       dedTcalc_tro => dedTcalc_tro_mah
       dpdTcalc_tro => dpdTcalc_tro_mah

       ! used for inflow BC (H_tot, s_tot, flow angle imposed)
       tcalc_Hstot => tcalc_Hstot_mah

       ! NOT USED yet
       dpdicalc_tro => dpdicalc_tro_mah
       rocalc_ep => rocalc_ep_mah
       rocalc_ps => rocalc_ps_mah
       tcalc_sro => tcalc_sro_mah
       tcalc_pro => tcalc_pro_mah
       tcalc_roero => tcalc_roero_mah   

       ! functions used to compute saturation curves
       vvol    => vvol_mah
       vvol_d1 => vvol_d1_mah
       vvol_d2 => vvol_d2_mah
       intpcalc_tro => intpcalc_tro_mah

       ! function used to compute total quantities from static quantities
       !stagnation_calc => stagnation_mah
       !call mpistop('Subroutine stagnation_mah not tested yet !',0)

    elseif (eos_type.eq.'prs') then
       
       ! primitive variables from conservative ones
       primitives => primitives_prs
       primitives_visc => primitives_visc_prs
       
       avcalc_tro => avcalc_tro_prs
       c2calc_tro => c2calc_tro_prs
       cpcalc_tro => cpcalc_tro_prs
       cvcalc_tro => cvcalc_tro_prs
       ecalc_tro => ecalc_tro_prs
       gcalc_tro => gcalc_tro_prs
       pcalc_tro => pcalc_tro_prs
       scalc_tro => scalc_tro_prs
       ecalc_pro => ecalc_pro_prs
       pcalc_roero => pcalc_roero_prs
       rocalc_pt => rocalc_pt_prs
       rocalc_st => rocalc_st_prs
       tcalc_ph  => tcalc_ph_prs
       dedrocalc_tro => dedrocalc_tro_prs
       dedTcalc_tro => dedTcalc_tro_prs
       dpdvcalc_tro => dpdvcalc_tro_prs
       dpdTcalc_tro => dpdTcalc_tro_prs

       ! used for inflow BC (H_tot, s_tot, flow angle imposed)
       tcalc_Hstot => tcalc_Hstot_prs

       ! NOT USED yet
       dpdicalc_tro => dpdicalc_tro_prs
       rocalc_ep => rocalc_ep_prs
       rocalc_ps => rocalc_ps_prs
       tcalc_sro => tcalc_sro_prs
       tcalc_pro => tcalc_pro_prs
       tcalc_roero => tcalc_roero_prs   

       ! functions used to compute saturation curves
       vvol    => vvol_prs
       vvol_d1 => vvol_d1_prs
       vvol_d2 => vvol_d2_prs
       intpcalc_tro => intpcalc_tro_prs

       ! function used to compute total quantities from static quantities
       stagnation_calc => stagnation_prs

    elseif (eos_type.eq.'swn') then
       
       ! primitive variables from conservative ones
       primitives => primitives_swn
       primitives_visc => primitives_visc_swn
       
       avcalc_tro => avcalc_tro_swn
       c2calc_tro => c2calc_tro_swn
       cpcalc_tro => cpcalc_tro_swn
       cvcalc_tro => cvcalc_tro_swn
       ecalc_tro => ecalc_tro_swn
       gcalc_tro => gcalc_tro_swn
       pcalc_tro => pcalc_tro_swn
       scalc_tro => scalc_tro_swn
       ecalc_pro => ecalc_pro_swn
       pcalc_roero => pcalc_roero_swn
       rocalc_pt => rocalc_pt_swn
       rocalc_st => rocalc_st_swn
       tcalc_ph  => tcalc_ph_swn
       dedTcalc_tro => dedTcalc_tro_swn
       dpdTcalc_tro => dpdTcalc_tro_swn

       ! used for inflow BC (H_tot, s_tot, flow angle imposed)
       tcalc_Hstot => tcalc_Hstot_swn

       ! NOT USED yet
       dpdicalc_tro => dpdicalc_tro_swn
       rocalc_ep => rocalc_ep_swn
       rocalc_ps => rocalc_ps_swn
       tcalc_sro => tcalc_sro_swn
       tcalc_pro => tcalc_pro_swn
       tcalc_roero => tcalc_roero_swn   

       ! functions used to compute saturation curves
       vvol    => vvol_swn
       vvol_d1 => vvol_d1_swn
       vvol_d2 => vvol_d2_swn
       intpcalc_tro => intpcalc_tro_swn

       ! function used to compute total quantities from static quantities
       ! stagnation_calc => stagnation_swn
       call mpistop('Subroutine stagnation_swn not implemented yet !',0)

    elseif (eos_type.eq.'swp') then
       
       ! primitive variables from conservative ones
       primitives => primitives_swp
       primitives_visc => primitives_visc_swp
       
       avcalc_tro => avcalc_tro_swp
       c2calc_tro => c2calc_tro_swp
       cpcalc_tro => cpcalc_tro_swp
       cvcalc_tro => cvcalc_tro_swp
       ecalc_tro => ecalc_tro_swp
       gcalc_tro => gcalc_tro_swp
       pcalc_tro => pcalc_tro_swp
       scalc_tro => scalc_tro_swp
       ecalc_pro => ecalc_pro_swp
       pcalc_roero => pcalc_roero_swp
       rocalc_pt => rocalc_pt_swp
       rocalc_st => rocalc_st_swp
       tcalc_ph  => tcalc_ph_swp
       dedTcalc_tro => dedTcalc_tro_swp
       dpdTcalc_tro => dpdTcalc_tro_swp

       ! used for inflow BC (H_tot, s_tot, flow angle imposed)
       tcalc_Hstot => tcalc_Hstot_swp

       ! NOT USED yet
       dpdicalc_tro => dpdicalc_tro_swp
       rocalc_ep => rocalc_ep_swp
       rocalc_ps => rocalc_ps_swp
       tcalc_sro => tcalc_sro_swp
       tcalc_pro => tcalc_pro_swp
       tcalc_roero => tcalc_roero_swp   

       ! functions used to compute saturation curves
       vvol    => vvol_swp
       vvol_d1 => vvol_d1_swp
       vvol_d2 => vvol_d2_swp
       intpcalc_tro => intpcalc_tro_swp

       ! function used to compute total quantities from static quantities
       ! stagnation_calc => stagnation_swp
       call mpistop('Subroutine stagnation_swp not implemented yet !',0)

    elseif (eos_type.eq.'ref') then

       ! primitive variables from conservative ones
       ! NOT WRITTEN YET
       !! primitives => primitives_ref
       !! primitives_visc => primitives_visc_ref

       avcalc_tro => avcalc_tro_ref
       c2calc_tro => c2calc_tro_ref
       cpcalc_tro => cpcalc_tro_ref
       cvcalc_tro => cvcalc_tro_ref
       ecalc_tro => ecalc_tro_ref
       gcalc_tro => gcalc_tro_ref
       pcalc_tro => pcalc_tro_ref
       scalc_tro => scalc_tro_ref
       ecalc_pro => ecalc_pro_ref
       pcalc_roero => pcalc_roero_ref
       rocalc_pt => rocalc_pt_ref
       rocalc_st => rocalc_st_ref
       tcalc_ph  => tcalc_ph_ref

       ! NOT USED yet
       dpdicalc_tro => dpdicalc_tro_ref
       rocalc_ep => rocalc_ep_ref
       rocalc_ps => rocalc_ps_ref
       tcalc_sro => tcalc_sro_ref
       tcalc_pro => tcalc_pro_ref
       tcalc_roero => tcalc_roero_ref

       ! functions used to compute saturation curves
       ! /!\ not used for computation of saturation curve
       vvol    => vvol_ref
       vvol_d1 => vvol_d1_ref
       vvol_d2 => vvol_d2_ref
       intpcalc_tro => intpcalc_tro_ref

       ! function used to compute total quantities from static quantities
       ! stagnation_calc => stagnation_ref
       !call mpistop('Subroutine stagnation_ref not implemented yet !',0)

    endif
  
    if (is_RANS) then
       if (model_RANS.eq.'SA') then
          primitives_rans => primitives_SA
       endif
    endif

  end subroutine assign_procedures_eos
  
end module mod_eos
