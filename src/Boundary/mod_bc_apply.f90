
!=================================================================================
module mod_bc_apply
!=================================================================================
  !> Module to apply Boundary Conditions
!=================================================================================
  implicit none
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

contains

  !===============================================================================
  subroutine bc_non_reflecting ! TO BE CHANGED -> bc_apply (except wall)
    !!    --> call bc_imin, .... with mod_interface ??
    !!    --> pb with corner and edges -> voidp ?? not used for bc_1pt ??
  !===============================================================================
    !> author: XG
    !> date: March 2021
    !> Apply Non Reflecting Boundary Conditions
  !===============================================================================
    use mod_interface
    use mod_charac
    use mod_bc_inlet_outlet
    use mod_turb_inlet_outlet
    !use mod_inlet_outlet
    implicit none
    ! ---------------------------------------------------------------------------
    ! ---------------------------------------------------------------------------

    if (is_TamDong3D) then
       
       ! 3-D Tam & Dong's Non Reflecting Boundary Conditions
       ! ===================================================
       
       ! Boundary Conditions for faces
       ! ----------------------------- 
       if (BC_face(1,1)%sort==-1) call bc_TD_imin
       if (BC_face(1,1)%sort==-2) call bc_TD_outflow_imin

       if (BC_face(1,2)%sort==-1) call bc_TD_imax
       if (BC_face(1,2)%sort==-2) call bc_TD_outflow_imax

       if (BC_face(2,1)%sort==-1) call bc_TD_jmin
       if (BC_face(2,1)%sort==-2) call bc_TD_outflow_jmin

       if (BC_face(2,2)%sort==-1) call bc_TD_jmax
       if (BC_face(2,2)%sort==-2) call bc_TD_outflow_jmax

       if (BC_face(3,1)%sort==-1) call bc_TD_kmin
       if (BC_face(3,1)%sort==-2) call bc_TD_outflow_kmin

       if (BC_face(3,2)%sort==-1) call bc_TD_kmax
       if (BC_face(3,2)%sort==-2) call bc_TD_outflow_kmax

       ! Boundary Conditions for edges along z
       ! -------------------------------------
       ! 1,1,1 : imin-jmin
       if (BC_edge(1,1,1)%sort==-1) call bc_TD_imin_jmin
       ! 1,1,2 : imin-jmax
       if (BC_edge(1,1,2)%sort==-1) call bc_TD_imin_jmax
       ! 1,2,1 : imax-jmin
       if (BC_edge(1,2,1)%sort==-1) call bc_TD_imax_jmin
       ! 1,2,2 : imax-jmax
       if (BC_edge(1,2,2)%sort==-1) call bc_TD_imax_jmax

       ! Boundary Conditions for edges along x
       ! -------------------------------------
       ! 2,1,1 : jmin-kmin
       if (BC_edge(2,1,1)%sort==-1) call bc_TD_jmin_kmin
       ! 2,1,2 : jmin-kmax
       if (BC_edge(2,1,2)%sort==-1) call bc_TD_jmin_kmax
       ! 2,2,1 : jmax-kmin
       if (BC_edge(2,2,1)%sort==-1) call bc_TD_jmax_kmin
       ! 2,2,2 : jmax-kmax
       if (BC_edge(2,2,2)%sort==-1) call bc_TD_jmax_kmax

       ! Boundary Conditions for edges along y
       ! -------------------------------------
       ! 3,1,1 : kmin-imin
       if (BC_edge(3,1,1)%sort==-1) call bc_TD_kmin_imin
       ! 3,1,2 : kmin-imax
       if (BC_edge(3,1,2)%sort==-1) call bc_TD_kmin_imax
       ! 3,2,1 : kmax-imin
       if (BC_edge(3,2,1)%sort==-1) call bc_TD_kmax_imin
       ! 3,2,2 : kmax-imax
       if (BC_edge(3,2,2)%sort==-1) call bc_TD_kmax_imax
       
       ! Boundary Conditions for corners
       ! -------------------------------
       ! 1,1,1 : imin-jmin-kmin
       if (BC_corner(1,1,1)%sort==-1) call bc_TD_imin_jmin_kmin
       ! 1,1,2 : imin-jmin-kmax
       if (BC_corner(1,1,2)%sort==-1) call bc_TD_imin_jmin_kmax
       ! 1,2,1 : imin-jmax-kmin
       if (BC_corner(1,2,1)%sort==-1) call bc_TD_imin_jmax_kmin
       ! 1,2,2 : imin-jmax-kmax
       if (BC_corner(1,2,2)%sort==-1) call bc_TD_imin_jmax_kmax
       ! 2,1,1 : imax-jmin-kmin
       if (BC_corner(2,1,1)%sort==-1) call bc_TD_imax_jmin_kmin
       ! 2,1,2 : imax-jmin-kmax
       if (BC_corner(2,1,2)%sort==-1) call bc_TD_imax_jmin_kmax
       ! 2,2,1 : imax-jmax-kmin
       if (BC_corner(2,2,1)%sort==-1) call bc_TD_imax_jmax_kmin
       ! 2,2,2 : imax-jmax-kmax
       if (BC_corner(2,2,2)%sort==-1) call bc_TD_imax_jmax_kmax
       
    else
       
       ! 2-D Tam & Dong's Non Reflecting Boundary Conditions
       ! ===================================================
       
       ! Boundary Conditions for faces (& edges for TamDong 1 point (-11))
       ! ----------------------------- 
       if ((BC_face(1,1)%sort==-1).or.(BC_face(1,1)%sort==-11)) call bc_TD_imin
       if (BC_face(1,1)%sort==-2) call bc_TD_outflow_imin

       if ((BC_face(1,2)%sort==-1).or.(BC_face(1,2)%sort==-11)) call bc_TD_imax
       if (BC_face(1,2)%sort==-2) call bc_TD_outflow_imax

       if ((BC_face(2,1)%sort==-1).or.(BC_face(2,1)%sort==-11)) call bc_TD_jmin
       if (BC_face(2,1)%sort==-2) call bc_TD_outflow_jmin

       if ((BC_face(2,2)%sort==-1).or.(BC_face(2,2)%sort==-11)) call bc_TD_jmax
       if (BC_face(2,2)%sort==-2) call bc_TD_outflow_jmax

       ! Boundary Conditions for edges
       ! ----------------------------- 
       ! 1,1,1 : imin-jmin
       if ((BC_edge(1,1,1)%sort==-1).or.(BC_edge(1,1,1)%sort==-2)) &
            call bc_TD_imin_jmin
       ! 1,1,2 : imin-jmax
       if ((BC_edge(1,1,2)%sort==-1).or.(BC_edge(1,1,2)%sort==-2)) &
            call bc_TD_imin_jmax
       ! 1,2,1 : imax-jmin
       if ((BC_edge(1,2,1)%sort==-1).or.(BC_edge(1,2,1)%sort==-2)) &
            call bc_TD_imax_jmin
       ! 1,2,2 : imax-jmax
       if ((BC_edge(1,2,2)%sort==-1).or.(BC_edge(1,2,2)%sort==-2)) &
            call bc_TD_imax_jmax
       
    endif
    
    ! NRCBC: Non Reflecting Characteristic Boundary Conditions (Cartesian only)
    ! ========================================================
    !if (.not.is_SBP) then
       if (BC_face(1,1)%sort==-3) then
          if (Mach.ge.1) then
             call bc_supersonic_inlet
          else
             call bc_charac_imin
          endif
             
!!$          if (Mach.ge.1) then
!!$             call bc_supersonic_inlet
!!$          else
!!$             call bc_charac_imin
!!$          endif
       endif

       if (BC_face(1,2)%sort==-3) call bc_charac_imax

       if (BC_face(2,1)%sort==-3) call bc_charac_jmin
       !if ((BC_face(2,1)%sort==-3).and.(ACT)) call bc_inlet1_jmin

       if (BC_face(2,2)%sort==-3) call bc_charac_jmax
!!$    else
!!$       if (BC_face(1,1)%sort==-3) call bc_charac_imin_SBP4
!!$
!!$       if (BC_face(1,2)%sort==-3) call bc_charac_imax_SBP4
!!$
!!$       if (BC_face(2,1)%sort==-3) call bc_charac_jmin_SBP4
!!$
!!$       if (BC_face(2,2)%sort==-3) call bc_charac_jmax_SBP4
!!$    endif

       if (BC_face(3,2)%sort==-3) call bc_charac_kmax

    ! Turbomachine BCs
    ! ================
       !if (TURB) then

    if (BC_face(1,1)%sort==-4) call bc_inlet_imin
    !if (BC_face(1,1)%sort==-4) call bc_inlet_prs_imin
    !if (BC_face(1,1)%sort==-4) call bc_inlet_test_imin
    !if (BC_face(1,1)%sort==-4) call bc_turb_inlet_imin
    if (BC_face(1,1)%sort==-5) call bc_outlet_imin

    if (BC_face(1,2)%sort==-4) call bc_inlet_imax
    if (BC_face(1,2)%sort==-5) call bc_outlet_imax
    !if (BC_face(1,2)%sort==-5) call bc_backpressure_prs_imax

    if (BC_face(2,1)%sort==-4) call bc_inlet_jmin
    if (BC_face(2,1)%sort==-5) call bc_outlet_jmin

    if (BC_face(2,2)%sort==-4) call bc_inlet_jmax
    if (BC_face(2,2)%sort==-5) call bc_outlet_jmax

    if (BC_face(3,1)%sort==-4) call bc_inlet_kmin
    if (BC_face(3,1)%sort==-5) call bc_outlet_kmin

    if (BC_face(3,2)%sort==-4) call bc_inlet_kmax
    if (BC_face(3,2)%sort==-5) call bc_outlet_kmax

!!$    else
!!$
!!$!!! TEMP FOR RANS !!!
!!$!!! REFLECTIVE RIEMANN BCs !!!
!!$
!!$       ! Riemann based BCs for external flows
!!$       ! ====================================
!!$       if (BC_face(1,1)%sort==-4) call bc_inlet_imin
!!$       if (BC_face(1,1)%sort==-5) call bc_outlet_imin
!!$
!!$       if (BC_face(1,2)%sort==-4) call bc_inlet_imax
!!$       if (BC_face(1,2)%sort==-5) call bc_outlet_imax
!!$
!!$       if (BC_face(2,1)%sort==-4) call bc_inlet_jmin
!!$       if (BC_face(2,1)%sort==-5) call bc_outlet_jmin
!!$
!!$       if (BC_face(2,2)%sort==-4) call bc_inlet_jmax
!!$       if (BC_face(2,2)%sort==-5) call bc_outlet_jmax
!!$    endif

    ! Axis BCs
    ! ========
    !if (BC_face(1,1)%sort==-8) call bc_axis_imin
    !if (BC_face(1,2)%sort==-8) call bc_axis_imax
    !if (BC_face(2,1)%sort==-8) call bc_axis_jmin
    !if (BC_face(2,2)%sort==-8) call bc_axis_jmax
!!$    if (BC_face(3,1)%sort==-8) call bc_axis_kmin
!!$    if (BC_face(3,2)%sort==-8) call bc_axis_kmax

  end subroutine bc_non_reflecting

  !===============================================================================
  subroutine bc_wall
  !===============================================================================
    !> author: XG
    !> date: March 2021
    !> Apply Wall Boundary Conditions
  !===============================================================================
    use mod_bc
    use mod_interface
    implicit none
    ! ---------------------------------------------------------------------------
    ! ---------------------------------------------------------------------------

    if (BC_face(1,1)%sort==0) call bc_wall_imin

    if (BC_face(1,2)%sort==0) call bc_wall_imax

    if (BC_face(2,1)%sort==0) call bc_wall_jmin

    if (BC_face(2,2)%sort==0) call bc_wall_jmax

    if (BC_face(3,1)%sort==0) call bc_wall_kmin

    if (BC_face(3,2)%sort==0) call bc_wall_kmax

  end subroutine bc_wall

end module mod_bc_apply
