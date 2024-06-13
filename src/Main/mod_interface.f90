!=================================================================================
module mod_interface
!=================================================================================
  !> Module to define interfaces and specification for procedure pointers
!=================================================================================
  implicit none
  ! ------------------------------------------------------------------------------
  ! Interfaces for procedure pointer
  ! ------------------------------------------------------------------------------
  abstract interface
     subroutine noarg_typ
       implicit none
     end subroutine noarg_typ
  end interface
  ! ------------------------------------------------------------------------------
  abstract interface
     subroutine comm_typ(var1,var2,var3,var4,var5)
       use mod_grid
       implicit none
       real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(inout) :: var1,var2,var3,var4,var5
     end subroutine comm_typ
  end interface
  ! ------------------------------------------------------------------------------
  abstract interface
     subroutine comm_typ_irs(var1,var2,var3,var4,var5,var6)
       use mod_grid
       implicit none
       real(wp), dimension(nx1_irs:nx2_irs,ny1_irs:ny2_irs,nz1_irs:nz2_irs)&
               , intent(inout) :: var1,var2,var3,var4,var5,var6
     end subroutine comm_typ_irs
  end interface
  ! ------------------------------------------------------------------------------


!!!!!!!!!!!!!!!!!!! Modif for RANS !!!!!!!!!!!!!!!!!!!!!!
  ! ------------------------------------------------------------------------------
  abstract interface ! -> made for 3D scalar derivatives
     subroutine comm_typ_rans(var1,var2,var3)
       use mod_grid
       implicit none
       real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(in) :: var1,var2,var3
     end subroutine comm_typ_rans
  end interface
  ! ------------------------------------------------------------------------------
  abstract interface ! -> made for single variable communication
     subroutine comm_single_rans(var1)
       use mod_grid
       implicit none
       real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(inout) :: var1
     end subroutine comm_single_rans
  end interface
  ! ------------------------------------------------------------------------------
  abstract interface
     subroutine grad_typ(var1,var2,var3,var4)
       use mod_grid
       implicit none
       real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(in)  :: var1
       real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v), intent(out) :: var2,var3,var4
     end subroutine grad_typ
  end interface
  ! ------------------------------------------------------------------------------

  procedure(noarg_typ), pointer :: stats           => null() &
                                 , num_dissip      => null() &
                                 , flux_euler      => null() &
                                 , flux_euler_rans => null() &
                                 , flux_visc       => null() &
                                 , flux_visc_rans  => null() &
                                 , grad_vel        => null() &
                                 , grad_T          => null() &
                                 , irs             => null() &
                                 , irs_solver      => null() &
                                 ! wall BC 
                                 , bc_wall_imin => null() &
                                 , bc_wall_imax => null() &
                                 , bc_wall_jmin => null() &
                                 , bc_wall_jmax => null() &
                                 , bc_wall_kmin => null() &
                                 , bc_wall_kmax => null() &
                                 ! T&D BC: faces 
                                 , bc_TD_imin => null() &
                                 , bc_TD_imax => null() &
                                 , bc_TD_jmin => null() &
                                 , bc_TD_jmax => null() &
                                 , bc_TD_kmin => null() &
                                 , bc_TD_kmax => null() &
                                 ! T&D BC: outflow faces 
                                 , bc_TD_outflow_imin => null() &
                                 , bc_TD_outflow_imax => null() &
                                 , bc_TD_outflow_jmin => null() &
                                 , bc_TD_outflow_jmax => null() &
                                 , bc_TD_outflow_kmin => null() &
                                 , bc_TD_outflow_kmax => null() &
                                 ! T&D BC: edges 
                                 , bc_TD_imin_jmin => null() &
                                 , bc_TD_imin_jmax => null() &
                                 , bc_TD_imax_jmin => null() &
                                 , bc_TD_imax_jmax => null() &
                                 , bc_TD_jmin_kmin => null() &
                                 , bc_TD_jmin_kmax => null() &
                                 , bc_TD_jmax_kmin => null() &
                                 , bc_TD_jmax_kmax => null() &
                                 , bc_TD_kmin_imin => null() &
                                 , bc_TD_kmin_imax => null() &
                                 , bc_TD_kmax_imin => null() &
                                 , bc_TD_kmax_imax => null() &
                                 ! T&D BC: corners 
                                 , bc_TD_imin_jmin_kmin => null() &
                                 , bc_TD_imin_jmin_kmax => null() &
                                 , bc_TD_imin_jmax_kmin => null() &
                                 , bc_TD_imin_jmax_kmax => null() &
                                 , bc_TD_imax_jmin_kmin => null() &
                                 , bc_TD_imax_jmin_kmax => null() &
                                 , bc_TD_imax_jmax_kmin => null() &
                                 , bc_TD_imax_jmax_kmax => null() &
                                 ! Inlet BC
                                 , bc_inlet_imin => null() &
                                 , bc_inlet_imax => null() &
                                 , bc_inlet_jmin => null() &
                                 , bc_inlet_jmax => null() &
                                 , bc_inlet_kmin => null() &
                                 , bc_inlet_kmax => null() &
                                 ! Outlet BC
                                 , bc_outlet_imin => null() &
                                 , bc_outlet_imax => null() &
                                 , bc_outlet_jmin => null() &
                                 , bc_outlet_jmax => null() &
                                 , bc_outlet_kmin => null() &
                                 , bc_outlet_kmax => null() &
                                 ! angular_periodicity BC 
                                 , bc_angular_periodicity => null() &
                                 , bc_angular_periodicity_n => null() &
                                 , bc_angular_periodicity_grad => null() &
                                 ! communications 
                                 , communication => null() &
                                 , communication_v => null() &
                                 ! IRS 
                                 , apply_irs_1 => null() &
                                 , apply_irs_2 => null() &
                                 , apply_irs_3 => null() &
                                 ! update numerical dissipation 
                                 , update_varn => null() &
                                 ! RANS
                                 , init_rans       => null() &
                                 , source_rans     => null() &
                                 , num_dissip_rans => null() &
                                 , update_var_rans => null() & 
                                 , length_scale_rans => null() &
                                 , start_runge_kutta_rans => null() &
                                 , runge_kutta_rans => null() &
                                 ! wall BC RANS
                                 , bc_wall_imin_rans => null() &
                                 , bc_wall_imax_rans => null() &
                                 , bc_wall_jmin_rans => null() &
                                 , bc_wall_jmax_rans => null() &
                                 , bc_wall_kmin_rans => null() &
                                 , bc_wall_kmax_rans => null() &
                                 ! Inlet RANS BC:
                                 , inlet_imin_rans => null() &
                                 , inlet_imax_rans => null() &
                                 , inlet_jmin_rans => null() &
                                 , inlet_jmax_rans => null() &
                                 , inlet_kmin_rans => null() &
                                 , inlet_kmax_rans => null() &
                                 ! Neumann RANS BC:
                                 , neu_imin_rans => null() &
                                 , neu_imax_rans => null() &
                                 , neu_jmin_rans => null() &
                                 , neu_jmax_rans => null() &
                                 , neu_kmin_rans => null() &
                                 , neu_kmax_rans => null() &
                                 ! IRS RANS
                                 , irs_rans         => null() &
                                 , irs_solver_rans  => null() &
                                 , apply_irs_1_rans => null() &
                                 , apply_irs_2_rans => null() &
                                 , apply_irs_3_rans => null() &
                                 , irs2_ngh_i_rans  => null() &
                                 , irs2_ngh_j_rans  => null() &
                                 , irs2_ngh_k_rans  => null() &
                                 , irs4_ngh_i_rans  => null() &
                                 , irs4_ngh_j_rans  => null() &
                                 , irs4_ngh_k_rans  => null() &
                                 ! wall-model
                                 , bc_wm_jmin => null() &
                                 , bc_wm_jmax => null() &
                                 ! update numerical dissipation
                                 , update_varn_rans => null()
  ! ------------------------------------------------------------------------------
  procedure(comm_typ), pointer :: communication_,communication_edges
  procedure(comm_typ), pointer :: communication_NS,communication_EW,communication_FB
  procedure(comm_typ_irs), pointer :: communication_inc
  procedure(comm_single_rans), pointer :: communication_inc_rans
  procedure(grad_typ), pointer :: grad_rans ! -> RANS gradient computation
  procedure(comm_typ_rans), pointer :: communication_grad_rans ! -> RANS gradient comm
  procedure(comm_single_rans), pointer :: communication_rans ! -> RANS variable comm
  ! ------------------------------------------------------------------------------

end module mod_interface
