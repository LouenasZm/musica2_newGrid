!===============================================================================
module mod_TamDong2d_1pt
!===============================================================================
  !> author: AB
  !> date: April 2022 (created from 5 points version)
  !> Non-reflecting boundary conditions of Tam and Dong
  !> - 2D / 2.5D (periodic) Cartesian version - version on 1 point
  !> /!\ NOT REGULAR - kept for dev (tests inflow conditions)
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
     !==========================================================================================
     module subroutine bc_TD2d_1pt_imin
     end subroutine bc_TD2d_1pt_imin
     !===========================================================================================
     module subroutine bc_TD2d_1pt_imax
     end subroutine bc_TD2d_1pt_imax
     !============================================================================================
     module subroutine bc_TD2d_1pt_jmin
     end subroutine bc_TD2d_1pt_jmin
     !=========================================================================================
     module subroutine bc_TD2d_1pt_jmax
     end subroutine bc_TD2d_1pt_jmax
     !=========================================================================================
  end interface

end module mod_TamDong2d_1pt
