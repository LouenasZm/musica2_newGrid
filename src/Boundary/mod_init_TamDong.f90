!===============================================================================
module mod_init_TamDong
!===============================================================================
  !> author: XG
  !> date: January 2022
  !> Initialization of Tam and Dong's boundary conditions
!===============================================================================
  use mod_flow        ! for: grid & primitive variables
  !!use mod_init_flow   ! for: is_vortex ** COMMENTED **
  use mod_constant    ! for: is_mean0
  use mod_mpi         ! for: iproc
  use mod_eos         ! for: calc c2
  implicit none
  !-----------------------------------------------------------------------------
  integer :: nghp3,nxmngh,nymngh,nzmngh,nxmnghp1,nymnghp1,nzmnghp1
  real(wp) :: xcr,ycr,zcr ! radiation center
  !-----------------------------------------------------------------------------

  interface
     !=============================================================================
     module subroutine init_bc_TD
     end subroutine init_bc_TD
     !=============================================================================
     module subroutine init_U0_TD
     end subroutine init_U0_TD
     !=============================================================================
     module subroutine comm_U0_TD
     end subroutine comm_U0_TD
     !=============================================================================
     module subroutine mean0_U0_TD
     end subroutine mean0_U0_TD
     !=============================================================================
     module subroutine init_param_TD2d
     end subroutine init_param_TD2d
     !=============================================================================
     module subroutine init_param_TD2d_1pt
     end subroutine init_param_TD2d_1pt
     !=============================================================================
     module subroutine init_param_TD3d
     end subroutine init_param_TD3d
     !=============================================================================
     module subroutine free_bc_TD2d
     end subroutine free_bc_TD2d
     !=============================================================================
     module subroutine free_bc_TD3d
     end subroutine free_bc_TD3d
     !=============================================================================
  end interface

end module mod_init_TamDong
