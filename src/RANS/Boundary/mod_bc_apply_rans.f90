!=================================================================================
module mod_bc_apply_rans
!=================================================================================
  !> Module to apply Boundary Conditions on turbulent quantities for RANS
!=================================================================================
  implicit none
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

contains

  !===============================================================================
  subroutine bc_rans
  !===============================================================================
    !> author: CM
    !> date: February 2022
    !> Apply Boundary Conditions 
  !===============================================================================
    use mod_interface
    use mod_charac_rans
    use mod_neumann_rans
    use mod_inlet_rans
    implicit none
    ! ---------------------------------------------------------------------------
    ! Apply Boundary Conditions on turbulent quantities for RANS
    ! ==========================================================
    ! Impose value at inlet
    if (BC_face(1,1)%sort==-1.or.BC_face(1,1)%sort==-4) call inlet_imin_rans !!!! should be -4

    if (BC_face(1,2)%sort==-1.or.BC_face(1,2)%sort==-4) call inlet_imax_rans

    if (BC_face(2,1)%sort==-1.or.BC_face(2,1)%sort==-4) call inlet_jmin_rans

    if (BC_face(2,2)%sort==-1.or.BC_face(2,2)%sort==-4) call inlet_jmax_rans

    if (BC_face(3,1)%sort==-1.or.BC_face(3,1)%sort==-4) call inlet_kmin_rans

    if (BC_face(3,2)%sort==-1.or.BC_face(3,2)%sort==-4) call inlet_kmax_rans

    ! Neumann
    if (BC_face(1,1)%sort==-3.or.BC_face(1,1)%sort==-2 & !!!! should be -1
                             .or.BC_face(1,1)%sort==-5) call neu_imin_rans

    if (BC_face(1,2)%sort==-3.or.BC_face(1,2)%sort==-2 & 
                             .or.BC_face(1,2)%sort==-5) call neu_imax_rans

    if (BC_face(2,1)%sort==-3.or.BC_face(2,1)%sort==-2 & 
                             .or.BC_face(2,1)%sort==-5) call neu_jmin_rans

    if (BC_face(2,2)%sort==-3.or.BC_face(2,2)%sort==-2 & 
                             .or.BC_face(2,2)%sort==-5) call neu_jmax_rans

    if (BC_face(3,1)%sort==-3.or.BC_face(3,1)%sort==-2 & 
                             .or.BC_face(3,1)%sort==-5) call neu_kmin_rans

    if (BC_face(3,2)%sort==-3.or.BC_face(3,2)%sort==-2 & 
                             .or.BC_face(3,2)%sort==-5) call neu_kmax_rans

  end subroutine bc_rans

  !===============================================================================
  subroutine bc_wall_rans
  !===============================================================================
    !> author: CM
    !> date: February 2022
    !> Apply Wall Boundary Conditions on turbulent quantities for RANS
  !===============================================================================
    use mod_bc
    use mod_interface
    implicit none
    ! ---------------------------------------------------------------------------
    ! ---------------------------------------------------------------------------

    if (BC_face(1,1)%sort==0) call bc_wall_imin_rans

    if (BC_face(1,2)%sort==0) call bc_wall_imax_rans

    if (BC_face(2,1)%sort==0) call bc_wall_jmin_rans

    if (BC_face(2,2)%sort==0) call bc_wall_jmax_rans

    if (BC_face(3,1)%sort==0) call bc_wall_kmin_rans

    if (BC_face(3,2)%sort==0) call bc_wall_kmax_rans

  end subroutine bc_wall_rans

end module mod_bc_apply_rans
