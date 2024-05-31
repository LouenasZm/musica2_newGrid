!===============================================================================
module mod_flux_euler
!===============================================================================
  !> author: XG
  !> date: January 2024
  !> Module with routines to compute Eulerian fluxes (inviscid part)
!=============================================================================== 
  use mod_coeff_deriv
  use mod_flow
  implicit none
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

  interface
     !===============================================================================
     ! Cartesian version
     !===============================================================================
     module subroutine flux_euler_11pts
     end subroutine flux_euler_11pts
     !===============================================================================
     module subroutine flux_euler_11pts_SBP4
     end subroutine flux_euler_11pts_SBP4
     !===============================================================================
     module subroutine flux_euler_9pts
     end subroutine flux_euler_9pts
     !===============================================================================
     module subroutine flux_euler_7pts
     end subroutine flux_euler_7pts
     !===============================================================================
     module subroutine flux_euler_5pts
     end subroutine flux_euler_5pts
     !===============================================================================
     ! Cartesian periodic [for tests only]
     !===============================================================================
     module subroutine flux_euler_11pts_perio
     end subroutine flux_euler_11pts_perio
     !===============================================================================
     ! 2.5D curvilinear version
     !===============================================================================
     module subroutine flux_euler_11pts_c
     end subroutine flux_euler_11pts_c
     !===============================================================================
     module subroutine flux_euler_11pts_SBP4_c
     end subroutine flux_euler_11pts_SBP4_c
     !===============================================================================
     module subroutine flux_euler_9pts_c
     end subroutine flux_euler_9pts_c
     !===============================================================================
     module subroutine flux_euler_7pts_c
     end subroutine flux_euler_7pts_c
     !===============================================================================
     module subroutine flux_euler_5pts_c
     end subroutine flux_euler_5pts_c
     !===============================================================================
     module subroutine flux_euler_5pts_SBP4_c
     end subroutine flux_euler_5pts_SBP4_c
     !===============================================================================
     ! wedge correction for 2.5D curvilinear version
     !===============================================================================
     module subroutine flux_euler_w_imin_jmin_11pts_c
     end subroutine flux_euler_w_imin_jmin_11pts_c
     !===============================================================================
     module subroutine flux_euler_w_imin_jmax_11pts_c
     end subroutine flux_euler_w_imin_jmax_11pts_c
     !===============================================================================
     module subroutine flux_euler_w_imax_jmin_11pts_c
     end subroutine flux_euler_w_imax_jmin_11pts_c
     !===============================================================================
     module subroutine flux_euler_w_imax_jmax_11pts_c
     end subroutine flux_euler_w_imax_jmax_11pts_c
     !===============================================================================
     module subroutine flux_euler_w_imin_jmin_11pts_SBP4_c
     end subroutine flux_euler_w_imin_jmin_11pts_SBP4_c
     !===============================================================================
     module subroutine flux_euler_w_imin_jmax_11pts_SBP4_c
     end subroutine flux_euler_w_imin_jmax_11pts_SBP4_c
     !===============================================================================
     module subroutine flux_euler_w_imax_jmin_11pts_SBP4_c
     end subroutine flux_euler_w_imax_jmin_11pts_SBP4_c
     !===============================================================================
     module subroutine flux_euler_w_imax_jmax_11pts_SBP4_c
     end subroutine flux_euler_w_imax_jmax_11pts_SBP4_c
     !===============================================================================
     module subroutine flux_euler_w_imin_jmin_9pts_c
     end subroutine flux_euler_w_imin_jmin_9pts_c
     !===============================================================================
     module subroutine flux_euler_w_imin_jmax_9pts_c
     end subroutine flux_euler_w_imin_jmax_9pts_c
     !===============================================================================
     module subroutine flux_euler_w_imax_jmin_9pts_c
     end subroutine flux_euler_w_imax_jmin_9pts_c
     !===============================================================================
     module subroutine flux_euler_w_imax_jmax_9pts_c
     end subroutine flux_euler_w_imax_jmax_9pts_c
     !===============================================================================
     module subroutine flux_euler_w_imin_jmin_5pts_c
     end subroutine flux_euler_w_imin_jmin_5pts_c
     !===============================================================================
     module subroutine flux_euler_w_imin_jmax_5pts_c
     end subroutine flux_euler_w_imin_jmax_5pts_c
     !===============================================================================
     module subroutine flux_euler_w_imax_jmin_5pts_c
     end subroutine flux_euler_w_imax_jmin_5pts_c
     !===============================================================================
     module subroutine flux_euler_w_imax_jmax_5pts_c
     end subroutine flux_euler_w_imax_jmax_5pts_c
     !===============================================================================
     ! Full 3D curvilinear version
     !===============================================================================
     module subroutine flux_euler_11pts_c3
     end subroutine flux_euler_11pts_c3
     !===============================================================================
     module subroutine flux_euler_11pts_SBP4_c3
     end subroutine flux_euler_11pts_SBP4_c3
     !===============================================================================
     module subroutine flux_euler_9pts_c3
     end subroutine flux_euler_9pts_c3
     !===============================================================================
     module subroutine flux_euler_7pts_c3
     end subroutine flux_euler_7pts_c3
     !===============================================================================
     module subroutine flux_euler_5pts_c3
     end subroutine flux_euler_5pts_c3
     !===============================================================================
  end interface

end module mod_flux_euler
