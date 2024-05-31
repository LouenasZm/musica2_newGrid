!================================================================================
module mod_bc_wall
!================================================================================
  !> Module to apply wall Boundary Conditions on conservative variables
  !> |_ version #1: all variables are imposed
  !> |_ version #2: advancement of rho at the wall => allow slip walls
!================================================================================
  use mod_flow  ! for: flow variables
  use mod_eos   ! for: rocalc_pt, ecalc_tro
  use mod_coeff_deriv ! for cextrp2, cextrp3
  implicit none
  !------------------------------------------------------------------------------
  real(wp), parameter :: onethird=1.0_wp/3.0_wp
  !------------------------------------------------------------------------------

  interface
     !==============================================================================
     ! Cartesian version #1 (all variables are imposed)
     !==============================================================================
     module subroutine bc_wall_imin_adiabatic_dpdn
     end subroutine bc_wall_imin_adiabatic_dpdn
     !==============================================================================
     module subroutine bc_wall_imin_isotherm_dpdn
     end subroutine bc_wall_imin_isotherm_dpdn
     !==============================================================================
     module subroutine bc_wall_imax_adiabatic_dpdn
     end subroutine bc_wall_imax_adiabatic_dpdn
     !==============================================================================
     module subroutine bc_wall_imax_isotherm_dpdn
     end subroutine bc_wall_imax_isotherm_dpdn
     !==============================================================================
     module subroutine bc_wall_jmin_adiabatic_dpdn
     end subroutine bc_wall_jmin_adiabatic_dpdn
     !==============================================================================
     module subroutine bc_wall_jmin_isotherm_dpdn
     end subroutine bc_wall_jmin_isotherm_dpdn
     !==============================================================================
     module subroutine bc_wall_jmax_adiabatic_dpdn
     end subroutine bc_wall_jmax_adiabatic_dpdn
     !==============================================================================
     module subroutine bc_wall_jmax_isotherm_dpdn
     end subroutine bc_wall_jmax_isotherm_dpdn
     !==============================================================================
     module subroutine bc_wall_kmin_adiabatic_dpdn
     end subroutine bc_wall_kmin_adiabatic_dpdn
     !==============================================================================
     module subroutine bc_wall_kmin_isotherm_dpdn
     end subroutine bc_wall_kmin_isotherm_dpdn
     !==============================================================================
     module subroutine bc_wall_kmax_adiabatic_dpdn
     end subroutine bc_wall_kmax_adiabatic_dpdn
     !==============================================================================
     module subroutine bc_wall_kmax_isotherm_dpdn
     end subroutine bc_wall_kmax_isotherm_dpdn
     !==============================================================================
     ! 2.5D curvilinear version #1 (all variables are imposed)
     !==============================================================================
     module subroutine bc_wall_imin_adiabatic_dpdn_c
     end subroutine bc_wall_imin_adiabatic_dpdn_c
     !==============================================================================
     module subroutine bc_wall_imin_isotherm_dpdn_c
     end subroutine bc_wall_imin_isotherm_dpdn_c
     !==============================================================================
     module subroutine bc_wall_imax_adiabatic_dpdn_c
     end subroutine bc_wall_imax_adiabatic_dpdn_c
     !==============================================================================
     module subroutine bc_wall_imax_isotherm_dpdn_c
     end subroutine bc_wall_imax_isotherm_dpdn_c
     !==============================================================================
     module subroutine bc_wall_jmin_adiabatic_dpdn_c
     end subroutine bc_wall_jmin_adiabatic_dpdn_c
     !==============================================================================
     module subroutine bc_wall_jmin_isotherm_dpdn_c
     end subroutine bc_wall_jmin_isotherm_dpdn_c
     !==============================================================================
     module subroutine bc_wall_jmax_adiabatic_dpdn_c
     end subroutine bc_wall_jmax_adiabatic_dpdn_c
     !==============================================================================
     module subroutine bc_wall_jmax_isotherm_dpdn_c
     end subroutine bc_wall_jmax_isotherm_dpdn_c
     !==============================================================================
     ! 3D curvilinear version #1 (all variables are imposed)
     !==============================================================================
     module subroutine bc_wall_imin_adiabatic_dpdn_c3
     end subroutine bc_wall_imin_adiabatic_dpdn_c3
     !==============================================================================
     module subroutine bc_wall_imin_isotherm_dpdn_c3
     end subroutine bc_wall_imin_isotherm_dpdn_c3
     !==============================================================================
     module subroutine bc_wall_imax_adiabatic_dpdn_c3
     end subroutine bc_wall_imax_adiabatic_dpdn_c3
     !==============================================================================
     module subroutine bc_wall_imax_isotherm_dpdn_c3
     end subroutine bc_wall_imax_isotherm_dpdn_c3
     !==============================================================================
     module subroutine bc_wall_jmin_adiabatic_dpdn_c3
     end subroutine bc_wall_jmin_adiabatic_dpdn_c3
     !==============================================================================
     module subroutine bc_wall_jmin_isotherm_dpdn_c3
     end subroutine bc_wall_jmin_isotherm_dpdn_c3
     !==============================================================================
     module subroutine bc_wall_jmax_adiabatic_dpdn_c3
     end subroutine bc_wall_jmax_adiabatic_dpdn_c3
     !==============================================================================
     module subroutine bc_wall_jmax_isotherm_dpdn_c3
     end subroutine bc_wall_jmax_isotherm_dpdn_c3
     !==============================================================================
     module subroutine bc_wall_kmin_adiabatic_dpdn_c3
     end subroutine bc_wall_kmin_adiabatic_dpdn_c3
     !==============================================================================
     module subroutine bc_wall_kmin_isotherm_dpdn_c3
     end subroutine bc_wall_kmin_isotherm_dpdn_c3
     !==============================================================================
     module subroutine bc_wall_kmax_adiabatic_dpdn_c3
     end subroutine bc_wall_kmax_adiabatic_dpdn_c3
     !==============================================================================
     module subroutine bc_wall_kmax_isotherm_dpdn_c3
     end subroutine bc_wall_kmax_isotherm_dpdn_c3
     !==============================================================================
     ! Cartesian version #2 (advancement of rho) => allow slip walls
     !==============================================================================
     module subroutine bc_wall_imin_slip
     end subroutine bc_wall_imin_slip
     !==============================================================================
     module subroutine bc_wall_imin_adiabatic
     end subroutine bc_wall_imin_adiabatic
     !==============================================================================
     module subroutine bc_wall_imin_isotherm
     end subroutine bc_wall_imin_isotherm
     !==============================================================================
     module subroutine bc_wall_imax_slip
     end subroutine bc_wall_imax_slip
     !==============================================================================
     module subroutine bc_wall_imax_adiabatic
     end subroutine bc_wall_imax_adiabatic
     !==============================================================================
     module subroutine bc_wall_imax_isotherm
     end subroutine bc_wall_imax_isotherm
     !==============================================================================
     module subroutine bc_wall_jmin_slip
     end subroutine bc_wall_jmin_slip
     !==============================================================================
     module subroutine bc_wall_jmin_adiabatic
     end subroutine bc_wall_jmin_adiabatic
     !==============================================================================
     module subroutine bc_wall_jmin_isotherm
     end subroutine bc_wall_jmin_isotherm
     !==============================================================================
     module subroutine bc_wall_jmax_slip
     end subroutine bc_wall_jmax_slip
     !==============================================================================
     module subroutine bc_wall_jmax_adiabatic
     end subroutine bc_wall_jmax_adiabatic
     !==============================================================================
     module subroutine bc_wall_jmax_isotherm
     end subroutine bc_wall_jmax_isotherm
     !==============================================================================
     module subroutine bc_wall_kmin_slip
     end subroutine bc_wall_kmin_slip
     !==============================================================================
     module subroutine bc_wall_kmin_adiabatic
     end subroutine bc_wall_kmin_adiabatic
     !==============================================================================
     module subroutine bc_wall_kmin_isotherm
     end subroutine bc_wall_kmin_isotherm
     !==============================================================================
     module subroutine bc_wall_kmax_slip
     end subroutine bc_wall_kmax_slip
     !==============================================================================
     module subroutine bc_wall_kmax_adiabatic
     end subroutine bc_wall_kmax_adiabatic
     !==============================================================================
     module subroutine bc_wall_kmax_isotherm
     end subroutine bc_wall_kmax_isotherm
     !==============================================================================
     ! 2.5D curvilinear version #2 (advancement of rho) => allow slip walls
     !==============================================================================
     module subroutine bc_wall_imin_slip_c
     end subroutine bc_wall_imin_slip_c
     !==============================================================================
     module subroutine bc_wall_imin_adiabatic_c
     end subroutine bc_wall_imin_adiabatic_c
     !==============================================================================
     module subroutine bc_wall_imax_slip_c
     end subroutine bc_wall_imax_slip_c
     !==============================================================================
     module subroutine bc_wall_imax_adiabatic_c
     end subroutine bc_wall_imax_adiabatic_c
     !==============================================================================
     module subroutine bc_wall_jmin_slip_c
     end subroutine bc_wall_jmin_slip_c
     !==============================================================================
     module subroutine bc_wall_jmin_adiabatic_c
     end subroutine bc_wall_jmin_adiabatic_c
     !==============================================================================
     module subroutine bc_wall_jmax_slip_c
     end subroutine bc_wall_jmax_slip_c
     !==============================================================================
     module subroutine bc_wall_jmax_adiabatic_c
     end subroutine bc_wall_jmax_adiabatic_c
     !==============================================================================
     ! 3D curvilinear version #2 (advancement of rho) => allow slip walls
     !==============================================================================
     module subroutine bc_wall_imin_slip_c3
     end subroutine bc_wall_imin_slip_c3
     !==============================================================================
     module subroutine bc_wall_imin_adiabatic_c3
     end subroutine bc_wall_imin_adiabatic_c3
     !==============================================================================
     module subroutine bc_wall_imax_slip_c3
     end subroutine bc_wall_imax_slip_c3
     !==============================================================================
     module subroutine bc_wall_imax_adiabatic_c3
     end subroutine bc_wall_imax_adiabatic_c3
     !==============================================================================
     module subroutine bc_wall_jmin_slip_c3
     end subroutine bc_wall_jmin_slip_c3
     !==============================================================================
     module subroutine bc_wall_jmin_adiabatic_c3
     end subroutine bc_wall_jmin_adiabatic_c3
     !==============================================================================
     module subroutine bc_wall_jmax_slip_c3
     end subroutine bc_wall_jmax_slip_c3
     !==============================================================================
     module subroutine bc_wall_jmax_adiabatic_c3
     end subroutine bc_wall_jmax_adiabatic_c3
     !==============================================================================
     module subroutine bc_wall_kmin_slip_c3
     end subroutine bc_wall_kmin_slip_c3
     !==============================================================================
     module subroutine bc_wall_kmin_adiabatic_c3
     end subroutine bc_wall_kmin_adiabatic_c3
     !==============================================================================
     module subroutine bc_wall_kmax_slip_c3
     end subroutine bc_wall_kmax_slip_c3
     !==============================================================================
     module subroutine bc_wall_kmax_adiabatic_c3
     end subroutine bc_wall_kmax_adiabatic_c3
     !==============================================================================
  end interface

end module mod_bc_wall
