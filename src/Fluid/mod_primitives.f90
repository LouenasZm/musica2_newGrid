!===============================================================================
module mod_primitives
!===============================================================================
  !> author: XG
  !> date: January 2024
  !> Module with routines to compute primitives variables from conservative ones
  !> [Nota: vectorization directive "!dir$ simd" commented]
!=============================================================================== 
  use mod_flow  ! for: primitives & conservative variables
  use mod_fluid ! for: fluid constants
  implicit none
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

  interface
     !===============================================================================
     ! PFG EoS (perfect gas eq.)
     !===============================================================================
     module subroutine primitives_pfg
     end subroutine primitives_pfg
     !===============================================================================
     module subroutine primitives_visc_pfg
     end subroutine primitives_visc_pfg
     !===============================================================================
     ! VDW EoS (van der Waals eq.)
     !===============================================================================
     module subroutine primitives_vdw
     end subroutine primitives_vdw
     !===============================================================================
     module subroutine primitives_visc_vdw
     end subroutine primitives_visc_vdw
     !===============================================================================
     ! PRSV EoS (Peng-Robinson-Stryjek-Vera eq.)
     !===============================================================================
     module subroutine primitives_prs
     end subroutine primitives_prs
     !===============================================================================
     module subroutine primitives_visc_prs
     end subroutine primitives_visc_prs
     !===============================================================================
     ! MAH EoS (Martin & Hou eq.)
     !===============================================================================
     module subroutine primitives_mah
     end subroutine primitives_mah
     !===============================================================================
     module subroutine primitives_visc_mah
     end subroutine primitives_visc_mah
     !===============================================================================
     ! SW non-polar EoS (Span-Wagner eq. for non-polar compounds)
     !===============================================================================
     module subroutine primitives_swn
     end subroutine primitives_swn
     !===============================================================================
     module subroutine primitives_visc_swn
     end subroutine primitives_visc_swn
     !===============================================================================
     ! SW polar EoS (Span-Wagner eq. for polar compounds)
     !===============================================================================
     module subroutine primitives_swp
     end subroutine primitives_swp
     !===============================================================================
     module subroutine primitives_visc_swp
     end subroutine primitives_visc_swp
     !===============================================================================
  end interface

end module mod_primitives
