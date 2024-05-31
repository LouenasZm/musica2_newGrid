!===============================================================================
module mod_gradient
!===============================================================================
  !> author: XG
  !> date: January 2024
  !> Module with routines to compute gradients of T & velocity components
!=============================================================================== 
  use mod_mpi
  use mod_constant
  use mod_coeff_deriv
  use mod_flow
  implicit none
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

  interface
     !===============================================================================
     ! Cartesian version
     !===============================================================================
     module subroutine grad_T_5pts
     end subroutine grad_T_5pts
     !===============================================================================
     module subroutine grad_T_3pts
     end subroutine grad_T_3pts
     !===============================================================================
     module subroutine grad_vel_5pts
     end subroutine grad_vel_5pts
     !===============================================================================
     module subroutine grad_vel_3pts
     end subroutine grad_vel_3pts
     !===============================================================================
     ! 2.5D curvilinear version
     !===============================================================================
     module subroutine grad_T_5pts_c
     end subroutine grad_T_5pts_c
     !===============================================================================
     module subroutine grad_T_3pts_c
     end subroutine grad_T_3pts_c
     !===============================================================================
     module subroutine grad_vel_5pts_c
     end subroutine grad_vel_5pts_c
     !===============================================================================
     module subroutine grad_vel_3pts_c
     end subroutine grad_vel_3pts_c
     !===============================================================================
     ! Full 3D curvilinear version
     !===============================================================================
     module subroutine grad_T_5pts_c3
     end subroutine grad_T_5pts_c3
     !===============================================================================
     module subroutine grad_T_3pts_c3
     end subroutine grad_T_3pts_c3
     !===============================================================================
     module subroutine grad_vel_5pts_c3
     end subroutine grad_vel_5pts_c3
     !===============================================================================
     module subroutine grad_vel_3pts_c3
     end subroutine grad_vel_3pts_c3
     !===============================================================================
  end interface

end module mod_gradient
