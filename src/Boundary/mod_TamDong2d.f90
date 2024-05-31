!===============================================================================
module mod_TamDong2d
!===============================================================================
  !> author: XG
  !> date: January 2024
  !> Non-reflecting boundary conditions of Tam and Dong
  !> - 2D / 2.5D (periodic) Cartesian version -
!=============================================================================== 
  use mod_coeff_deriv  ! for: coeff Tam & Webb scheme
  use mod_time         ! for: irk
  use mod_init_TamDong ! already include:
                       !   * mod_flow (for grid & flow variables)
                       !   * mod_constant (for is_eigenmode,is_RFM,is_mean_ref)
                       !   * mod_mpi (for iproc)
                       !   * mod_eos (for calc c2,av,cp)
                       ! for: useful indices
                       !      nghp3=ngh+3
                       !      nxmngh=nx-ngh
                       !      nymngh=ny-ngh
                       !      nxmnghp1=nx-ngh+1
                       !      nymnghp1=ny-ngh+1
  implicit none
  !-----------------------------------------------------------------------------
  real(wp), dimension(:,:), pointer :: ir,cosphi,sinphi
  !-----------------------------------------------------------------------------

  interface
     !===============================================================================
     ! Faces - version with Tam & Webb's DRP boundary schemes
     !===============================================================================
     module subroutine bc_TD2d_imin
     end subroutine bc_TD2d_imin
     !===============================================================================
     module subroutine bc_TD2d_imax
     end subroutine bc_TD2d_imax
     !===============================================================================
     module subroutine bc_TD2d_jmin
     end subroutine bc_TD2d_jmin
     !===============================================================================
     module subroutine bc_TD2d_jmax
     end subroutine bc_TD2d_jmax
     !===============================================================================
     ! Edges - version with Tam & Webb's DRP boundary schemes
     !===============================================================================
     module subroutine bc_TD2d_imin_jmin
     end subroutine bc_TD2d_imin_jmin
     !===============================================================================
     module subroutine bc_TD2d_imin_jmax
     end subroutine bc_TD2d_imin_jmax
     !===============================================================================
     module subroutine bc_TD2d_imax_jmin
     end subroutine bc_TD2d_imax_jmin
     !===============================================================================
     module subroutine bc_TD2d_imax_jmax
     end subroutine bc_TD2d_imax_jmax
     !===============================================================================
     ! Faces - version with SBP4 boundary schemes [dev, not for regular use]
     !===============================================================================
      module subroutine bc_TD2d_imin_SBP4
     end subroutine bc_TD2d_imin_SBP4
     !===============================================================================
     module subroutine bc_TD2d_imax_SBP4
     end subroutine bc_TD2d_imax_SBP4
     !===============================================================================
     module subroutine bc_TD2d_jmin_SBP4
     end subroutine bc_TD2d_jmin_SBP4
     !===============================================================================
     module subroutine bc_TD2d_jmax_SBP4
     end subroutine bc_TD2d_jmax_SBP4
     !===============================================================================
     ! Edges - version with SBP4 boundary schemes [dev, not for regular use]
     !===============================================================================
     module subroutine bc_TD2d_imin_jmin_SBP4
     end subroutine bc_TD2d_imin_jmin_SBP4
     !===============================================================================
     module subroutine bc_TD2d_imin_jmax_SBP4
     end subroutine bc_TD2d_imin_jmax_SBP4
     !===============================================================================
     module subroutine bc_TD2d_imax_jmin_SBP4
     end subroutine bc_TD2d_imax_jmin_SBP4
     !===============================================================================
     module subroutine bc_TD2d_imax_jmax_SBP4
     end subroutine bc_TD2d_imax_jmax_SBP4
     !===============================================================================
     ! Outflow T&D conditions ~~~> restricted to imax (exit = right condition)
     !===============================================================================
     module subroutine bc_TD2d_outflow_imax
     end subroutine bc_TD2d_outflow_imax
     !===============================================================================
  end interface

end module mod_TamDong2d
