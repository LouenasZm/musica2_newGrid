!================================================================================
module mod_ngh
!================================================================================
  !> Module to share ngh, ngh_v, ngh_irs, is_2D, is_curv
!================================================================================
  use precision
  implicit none
  ! -----------------------------------------------------------------------------
  ! Ghost points
  ! -----------------------------------------------------------------------------
  ! number of ghost points corresponding to the stencil for Euler terms (inviscid part)
  integer :: ngh,ngh_r ! (also the total number of ghost points for the code) RANS ghost cells too
  ! ----------------------------------------------------------------------------- 
  ! number of ghost points corresponding to the stencil for viscous terms (viscous part)
  integer :: ngh_v,ngh_v_r ! (should be lower than ngh) RANS ghost cells too
  ! -----------------------------------------------------------------------------
  ! number of ghost points for IRS parallelization (per face; corresponding to 2*CFL+1)
  integer :: ngh_irs(6) ! (should be greater than 2)
  integer :: nghirs ! (should be greater than 2)
  ! -----------------------------------------------------------------------------
  ! Solver options
  ! -----------------------------------------------------------------------------
  ! dimension: 2-D or 3D
  integer :: dim
  ! indicator to call 2-D curvilinear solver
  logical  :: is_curv
  ! indicator to call 3-D curvilinear solver
  logical  :: is_curv3
  ! indicator to restrict to 2D grid
  logical  :: is_2d
  ! -----------------------------------------------------------------------------
end module mod_ngh
