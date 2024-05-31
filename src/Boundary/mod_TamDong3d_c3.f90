!===============================================================================
module mod_TamDong3d_c3
!===============================================================================
  !> author: XG
  !> date: January 2024
  !> Non-reflecting boundary conditions of Tam and Dong
  !> - 3D curvilinear version -
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
     module subroutine bc_TD3d_imin_c3
     end subroutine bc_TD3d_imin_c3
     !===============================================================================
     module subroutine bc_TD3d_imax_c3
     end subroutine bc_TD3d_imax_c3
     !===============================================================================
     module subroutine bc_TD3d_jmin_c3
     end subroutine bc_TD3d_jmin_c3
     !===============================================================================
     module subroutine bc_TD3d_jmax_c3
     end subroutine bc_TD3d_jmax_c3
     !===============================================================================
     module subroutine bc_TD3d_kmin_c3
     end subroutine bc_TD3d_kmin_c3
     !===============================================================================
     module subroutine bc_TD3d_kmax_c3
     end subroutine bc_TD3d_kmax_c3
     !===============================================================================
     ! Edges along x 
     !===============================================================================
     module subroutine bc_TD3d_jmin_kmin_c3
     end subroutine bc_TD3d_jmin_kmin_c3
     !===============================================================================
     module subroutine bc_TD3d_jmin_kmax_c3
     end subroutine bc_TD3d_jmin_kmax_c3
     !===============================================================================
     module subroutine bc_TD3d_jmax_kmin_c3
     end subroutine bc_TD3d_jmax_kmin_c3
     !===============================================================================
     module subroutine bc_TD3d_jmax_kmax_c3
     end subroutine bc_TD3d_jmax_kmax_c3
     !===============================================================================
     ! Edges along y
     !===============================================================================
     module subroutine bc_TD3d_kmin_imin_c3
     end subroutine bc_TD3d_kmin_imin_c3
     !===============================================================================
     module subroutine bc_TD3d_kmin_imax_c3
     end subroutine bc_TD3d_kmin_imax_c3
     !===============================================================================
     module subroutine bc_TD3d_kmax_imin_c3
     end subroutine bc_TD3d_kmax_imin_c3
     !===============================================================================
     module subroutine bc_TD3d_kmax_imax_c3
     end subroutine bc_TD3d_kmax_imax_c3
     !===============================================================================
     ! Edges along z 
     !===============================================================================
     module subroutine bc_TD3d_imin_jmin_c3
     end subroutine bc_TD3d_imin_jmin_c3
     !===============================================================================
     module subroutine bc_TD3d_imin_jmax_c3
     end subroutine bc_TD3d_imin_jmax_c3
     !===============================================================================
     module subroutine bc_TD3d_imax_jmin_c3
     end subroutine bc_TD3d_imax_jmin_c3
     !===============================================================================
     module subroutine bc_TD3d_imax_jmax_c3
     end subroutine bc_TD3d_imax_jmax_c3
     !===============================================================================
     ! Corners
     !===============================================================================
     module subroutine bc_TD3d_imin_jmin_kmin_c3
     end subroutine bc_TD3d_imin_jmin_kmin_c3
     !===============================================================================
     module subroutine bc_TD3d_imin_jmin_kmax_c3
     end subroutine bc_TD3d_imin_jmin_kmax_c3
     !===============================================================================
     module subroutine bc_TD3d_imin_jmax_kmin_c3
     end subroutine bc_TD3d_imin_jmax_kmin_c3
     !===============================================================================
     module subroutine bc_TD3d_imin_jmax_kmax_c3
     end subroutine bc_TD3d_imin_jmax_kmax_c3
     !===============================================================================
     module subroutine bc_TD3d_imax_jmin_kmin_c3
     end subroutine bc_TD3d_imax_jmin_kmin_c3
     !===============================================================================
     module subroutine bc_TD3d_imax_jmin_kmax_c3
     end subroutine bc_TD3d_imax_jmin_kmax_c3
     !===============================================================================
     module subroutine bc_TD3d_imax_jmax_kmin_c3
     end subroutine bc_TD3d_imax_jmax_kmin_c3
     !===============================================================================
     module subroutine bc_TD3d_imax_jmax_kmax_c3
     end subroutine bc_TD3d_imax_jmax_kmax_c3
     !===============================================================================
  end interface

end module mod_TamDong3d_c3
