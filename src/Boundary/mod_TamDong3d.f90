!===============================================================================
module mod_TamDong3d
!===============================================================================
  !> author: XG
  !> date: January 2024
  !> Non-reflecting boundary conditions of Tam and Dong
  !> - 3D Cartesian version -
!=============================================================================== 
  use mod_coeff_deriv  ! for: coeff Tam & Webb scheme
  use mod_time         ! for: irk
  use mod_init_TamDong ! already include:
                       !   * mod_flow (for grid & flow variables)
                       !   * mod_constant (for is_eigenmode,is_RFM)
                       !   * mod_mpi (for iproc)
                       !   * mod_eos (for calc c2,av,cp)
                       ! for: useful indices
                       !      nghp3=ngh+3
                       !      nxmngh=nx-ngh
                       !      nymngh=ny-ngh
                       !      nzmngh=nz-ngh
                       !      nxmnghp1=nx-ngh+1
                       !      nymnghp1=ny-ngh+1
                       !      nzmnghp1=nz-ngh+1
  implicit none
  !-----------------------------------------------------------------------------
  real(wp), dimension(:,:,:), pointer :: ir,cosphi,sinphi,costeta,sinteta
  real(wp), dimension(:,:,:), pointer :: costcosp,costsinp,sintcosp,sintsinp
  !-----------------------------------------------------------------------------

  interface
     !===============================================================================
     ! Faces
     !===============================================================================
     module subroutine bc_TD3d_imin
     end subroutine bc_TD3d_imin
     !===============================================================================
     module subroutine bc_TD3d_imax
     end subroutine bc_TD3d_imax
     !===============================================================================
     module subroutine bc_TD3d_jmin
     end subroutine bc_TD3d_jmin
     !===============================================================================
     module subroutine bc_TD3d_jmax
     end subroutine bc_TD3d_jmax
     !===============================================================================
     module subroutine bc_TD3d_kmin
     end subroutine bc_TD3d_kmin
     !===============================================================================
     module subroutine bc_TD3d_kmax
     end subroutine bc_TD3d_kmax
     !===============================================================================
     ! Edges along x
     !===============================================================================
     module subroutine bc_TD3d_jmin_kmin
     end subroutine bc_TD3d_jmin_kmin
     !===============================================================================
     module subroutine bc_TD3d_jmin_kmax
     end subroutine bc_TD3d_jmin_kmax
     !===============================================================================
     module subroutine bc_TD3d_jmax_kmin
     end subroutine bc_TD3d_jmax_kmin 
     !===============================================================================
     module subroutine bc_TD3d_jmax_kmax
     end subroutine bc_TD3d_jmax_kmax
     !===============================================================================
     ! Edges along y
     !===============================================================================
     module subroutine bc_TD3d_kmin_imin
     end subroutine bc_TD3d_kmin_imin
     !===============================================================================
     module subroutine bc_TD3d_kmin_imax
     end subroutine bc_TD3d_kmin_imax
     !===============================================================================
     module subroutine bc_TD3d_kmax_imin
     end subroutine bc_TD3d_kmax_imin
     !===============================================================================
     module subroutine bc_TD3d_kmax_imax
     end subroutine bc_TD3d_kmax_imax
     !===============================================================================
     ! Edges along z
     !===============================================================================
     module subroutine bc_TD3d_imin_jmin
     end subroutine bc_TD3d_imin_jmin
     !===============================================================================
     module subroutine bc_TD3d_imin_jmax
     end subroutine bc_TD3d_imin_jmax
     !===============================================================================
     module subroutine bc_TD3d_imax_jmin
     end subroutine bc_TD3d_imax_jmin
     !===============================================================================
     module subroutine bc_TD3d_imax_jmax
     end subroutine bc_TD3d_imax_jmax
     !===============================================================================
     ! Corners
     !===============================================================================
     module subroutine bc_TD3d_imin_jmin_kmin
     end subroutine bc_TD3d_imin_jmin_kmin
     !===============================================================================
     module subroutine bc_TD3d_imin_jmin_kmax
     end subroutine bc_TD3d_imin_jmin_kmax
     !===============================================================================
     module subroutine bc_TD3d_imin_jmax_kmin
     end subroutine bc_TD3d_imin_jmax_kmin
     !===============================================================================
     module subroutine bc_TD3d_imin_jmax_kmax
     end subroutine bc_TD3d_imin_jmax_kmax
     !===============================================================================
     module subroutine bc_TD3d_imax_jmin_kmin
     end subroutine bc_TD3d_imax_jmin_kmin
     !===============================================================================
     module subroutine bc_TD3d_imax_jmin_kmax
     end subroutine bc_TD3d_imax_jmin_kmax
     !===============================================================================
     module subroutine bc_TD3d_imax_jmax_kmin
     end subroutine bc_TD3d_imax_jmax_kmin
     !===============================================================================
     module subroutine bc_TD3d_imax_jmax_kmax
     end subroutine bc_TD3d_imax_jmax_kmax
     !===============================================================================
  end interface

end module mod_TamDong3d
