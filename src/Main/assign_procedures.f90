!=================================================================================
subroutine assign_procedures
!=================================================================================
  !> Routine for procedure pointer assignment
!=================================================================================
  use mod_interface ! <- procedure pointer interfaces
  use mod_utils     ! <- procedure voidp
  use mod_io
  use mod_constant
  use mod_tranprop
  ! Numerical dissipation (at the last substep)
  use mod_filtering
  use mod_filtering_shock
  use mod_artvisc
  use mod_artvisc_shock
  ! Numerical dissipation (added to increments)
  use mod_filtering_inc
  use mod_artvisc_inc
  use mod_artvisc_shock_inc
  ! Communications
  use mod_comm
  use mod_comm1
  ! Gradients and Fluxes
  ! ====================
  ! inviscid fluxes
  use mod_flux_euler
  ! viscous fluxes
  use mod_flux_visc
  ! temperature & velocity gradients
  use mod_gradient
  ! Boundary conditions
  ! ===================
  ! set periodicities
  use mod_bc_periodicity
  ! wall BCs
  use mod_bc_wall
  ! 2D/2.5D versions of TD
  use mod_TamDong2d
  use mod_TamDong2d_1pt ! [still plugged but for tests only]
  use mod_TamDong2d_c
  ! 3D versions of TD [Nota: outflow cond. not implemented]
  use mod_TamDong3d
  use mod_TamDong3d_c3
  ! 1D degenerescence of TD [test XG saved as _dev]
  !use mod_TamDong1d_faces_c3
  !use mod_TamDong1d_edgez_c3
  ! Riemann turbomachinery BBC
  use mod_bc_inlet_outlet
  !use mod_turb_inlet_outlet
  !use mod_riemann
  ! RANS modules
  use mod_irs_d_rans
  use mod_runge_kutta_rans
  use mod_bc_wall_rans
  use mod_inlet_rans
  use mod_neumann_rans
  use mod_artvisc_rans
  use mod_artvisc_shock_rans
  use mod_filtering_rans
  use mod_filtering_shock_rans
  use warnstop
  use mod_turb_model_length_scale
  use mod_wall_model
  implicit none
  ! ------------------------------------------------------------------------------
  external :: flux_visc_5pts_SM &
       , grad_scalar_3pts_c, grad_scalar_3pts_c3 &
       , grad_scalar_5pts_c, grad_scalar_5pts_c3 &
       , stats_tgv, stats_chit, stats_chan, stats_stbl, stats_cyl, stats_sphere, stats_shit, stats_act, stats_turb, stats_turb3 &
       , filtering_fluct_11pts, filtering_11pts_overlap &
       , update_varn_filter, update_varn_artvisc &
       ! Spalart Allmaras procedures
       , flux_euler_3pts_SA_c, flux_visc_3pts_SA_c, flux_euler_5pts_SA_c, flux_visc_5pts_SA_c &
       , flux_euler_3pts_SA_c3, flux_visc_3pts_SA_c3, flux_euler_5pts_SA_c3, flux_visc_5pts_SA_c3 &
       , source_SA, source_SA_criv, source_SA_transition_algebraic, init_SA, init_SA_transition_algebraic ,update_var_SA, update_varn_artvisc_SA, update_varn_filter_SA &
       ! SA-Gamma-Re_theta procedures
       , flux_euler_3pts_SA_transition_c, flux_visc_3pts_SA_transition_c, flux_euler_5pts_SA_transition_c, flux_visc_5pts_SA_transition_c &
       , flux_euler_3pts_SA_transition_c3, flux_visc_3pts_SA_transition_c3, flux_euler_5pts_SA_transition_c3, flux_visc_5pts_SA_transition_c3 &
       , source_SA_transition, init_SA_transition, update_var_SA_transition, update_varn_artvisc_SA_transition

  ! ------------------------------------------------------------------------------

  ! Define cases
  ! ============
  if (TGV) then
     stats => stats_tgv
  elseif (CHIT) then
     stats => stats_chit
  elseif (CHAN) then
     stats => stats_chan
  elseif (STBL) then
     stats => stats_stbl
  elseif (TURB) then
     if (is_curv3) then
        stats => stats_turb3
     else
        stats => stats_turb
     endif
  elseif (CYL) then
     if (is_curv3) then
        stats => stats_sphere
     else
        stats => stats_cyl
     endif
  elseif (SHIT) then
     stats => stats_shit
  elseif (ACT) then
     stats => stats_act
  elseif (SRC.or.CAV.or.PHILL) then
     stats => voidp
  elseif (LE.or.T3C) then
     stats => stats_stbl
  elseif (TE) then
     stats => stats_turb
  endif
  
  ! Define numerical dissipation
  ! ============================
  if (is_dissip_in_increments) then
     if (is_shock) then
        ! Artificial viscosity (DNC-Jameson)
        ! ----------------------------------
        select case (stencil)
        case (11)
           num_dissip => artvisc_o9_shock_inc
        case (9)
           call mpistop('artvisc_o7_shock to be written!', 0)
           !num_dissip => artvisc_o7_shock_inc
        case (7)
           call mpistop('artvisc_o5_shock to be written!', 0)
           !num_dissip => artvisc_o5_shock_inc
        case (5)
           call mpistop('artvisc_o3_shock to be written!', 0)
           !num_dissip => artvisc_o3_shock_inc
        case default
           call mpistop('Wrong choice for artificial dissipation !', 0)
        end select
     else  
        if (is_SF) then ! Selective filtering
           ! -------------------
           select case (stencil)
           case (11)
              num_dissip => filtering_11pts_inc
           case (9)
              num_dissip => filtering_9pts_inc
           case (7)
              num_dissip => filtering_7pts_inc
           case (5)
              num_dissip => filtering_5pts_inc
           case default
              call mpistop('Wrong choice for filtering !', 0)
           end select
           !if (CHIT.or.TGV) num_dissip => filtering_11pts_periodic
        else ! Artificial viscosity (DNC)
           ! --------------------------
           select case (stencil)
           case (11)
              num_dissip => artvisc_o9_inc
           case (9)
              num_dissip => artvisc_o7_inc
           case (7)
              num_dissip => artvisc_o5_inc
           case (5)
              num_dissip => artvisc_o3_inc
           case default
              call mpistop('Wrong choice for artificial dissipation !', 0)
           end select
        endif
     endif
  else
     if (is_shock) then
        if (is_SF) then ! Selective filtering
                        ! -------------------
           select case (stencil)
           case (11)
              num_dissip => filter_o10_shock
           case (9)
              num_dissip => filter_o8_shock
           case (7)
              num_dissip => filter_o6_shock
           case (5)
              num_dissip => filter_o4_shock
           case default
              call mpistop('Wrong choice for artificial dissipation !', 0)
           end select
           ! updating procedure
           ! ------------------
           update_varn => update_varn_filter
        else
           ! Artificial viscosity (DNC-Jameson)
           ! ----------------------------------
           select case (stencil)
           case (11)
              num_dissip => artvisc_o9_shock
           case (9)
              num_dissip => artvisc_o7_shock
           case (7)
              num_dissip => artvisc_o5_shock
           case (5)
              num_dissip => artvisc_o3_shock
           case default
              call mpistop('Wrong choice for artificial dissipation !', 0)
           end select
           ! updating procedure
           ! ------------------
           update_varn => update_varn_artvisc
        endif
     else  
        if (is_SF) then ! Selective filtering
           ! -------------------
           select case (stencil)
           case (11)
              num_dissip => filtering_11pts
              !num_dissip => filtering_11pts_overlap
              !num_dissip => filtering_fluct_11pts
              !num_dissip => filtering_11pts_periodic
           case (9)
              num_dissip => filtering_9pts
           case (7)
              num_dissip => filtering_7pts
           case (5)
              num_dissip => filtering_5pts
           case default
              call mpistop('Wrong choice for filtering !', 0)
           end select
           !if (CHIT.or.TGV) num_dissip => filtering_11pts_periodic
           ! updating procedure
           ! ------------------
           update_varn => update_varn_filter
        else ! Artificial viscosity (DNC)
           ! --------------------------
           select case (stencil)
           case (11)
              num_dissip => artvisc_o9
           case (9)
              num_dissip => artvisc_o7
           case (7)
              num_dissip => artvisc_o5
           case (5)
              num_dissip => artvisc_o3
           case default
              call mpistop('Wrong choice for artificial dissipation !', 0)
           end select
           ! updating procedure
           ! ------------------
           update_varn => update_varn_artvisc
        endif
     endif
  endif

  ! Define order of Euler fluxes
  ! ============================
  select case (stencil)
  case (11)
     if (is_curv3) then
        if (is_SBP) then
           flux_euler => flux_euler_11pts_SBP4_c3
        else
           flux_euler => flux_euler_11pts_c3
        endif
     elseif (is_curv) then
        if (is_SBP) then
           flux_euler => flux_euler_11pts_SBP4_c
        else
           flux_euler => flux_euler_11pts_c
        endif
     else
        if (is_SBP) then
           flux_euler => flux_euler_11pts_SBP4
        else
           flux_euler => flux_euler_11pts
        endif
     endif
  case (9)
     if (is_curv3) then
        flux_euler => flux_euler_9pts_c3
     elseif (is_curv) then
        flux_euler => flux_euler_9pts_c
     else
        flux_euler => flux_euler_9pts
     endif
  case (7)
     if (is_curv3) then
        flux_euler => flux_euler_7pts_c3
     elseif (is_curv) then
        flux_euler => flux_euler_7pts_c
     else
        flux_euler => flux_euler_7pts
     endif
  case (5)
     if (is_curv3) then
        flux_euler => flux_euler_5pts_c3
     elseif (is_curv) then
        if (is_SBP) then
           flux_euler => flux_euler_5pts_SBP4_c
        else
           flux_euler => flux_euler_5pts_c
        endif
     else
        flux_euler => flux_euler_5pts
     endif
  case default
     call mpistop('Assign flux_euler not implemented !', 0)
  end select

  ! Define order of viscous fluxes
  ! ==============================
  select case (iorder_visc)
  case (0) ! <-- Euler
     grad_vel  => voidp
     grad_T    => voidp
     flux_visc => voidp
     if (is_shock) then
        if (is_curv) then
           grad_vel  => grad_vel_5pts_c
        else
           grad_vel  => grad_vel_5pts  
        endif
     endif
     if (idepart==POST_PROCESSING) then
        if (is_curv3) then
           grad_vel  => grad_vel_5pts_c3
        elseif (is_curv) then
           grad_vel  => grad_vel_5pts_c
        else
           grad_vel  => grad_vel_5pts
        endif
     endif
  case (2) ! <-- order 2
     if (is_curv3) then
        grad_vel  => grad_vel_3pts_c3
        grad_T    => grad_T_3pts_c3
        flux_visc => flux_visc_3pts_c3
     else
        if (is_curv) then
           grad_vel  => grad_vel_3pts_c
           grad_T    => grad_T_3pts_c
           if (is_sgs_model) then
              call mpistop('LES models not implemented in curvilinear', 0)
           else
              flux_visc => flux_visc_3pts_c
           endif
        else
           grad_vel  => grad_vel_3pts
           grad_T    => grad_T_3pts
           if (is_sgs_model) then
              call mpistop('LES models not implemented for 3pts stencil', 0)
           else
              flux_visc => flux_visc_3pts
           endif
        endif
     endif
  case (4) ! <-- order 4
     if (is_curv3) then
        grad_vel  => grad_vel_5pts_c3
        grad_T    => grad_T_5pts_c3
        if (is_sgs_model) then
           call mpistop('LES models not implemented in 3D curvilinear', 0)
        else
           flux_visc => flux_visc_5pts_c3
        endif
     else
        if (is_curv) then
           grad_vel  => grad_vel_5pts_c
           grad_T    => grad_T_5pts_c
           if (is_sgs_model) then
              call mpistop('LES models not implemented in curvilinear', 0)
           else
              flux_visc => flux_visc_5pts_c
           end if
        else
           grad_vel  => grad_vel_5pts
           grad_T    => grad_T_5pts
           if (is_sgs_model) then
              if (model_LES.eq."SM") then
                 flux_visc => flux_visc_5pts_SM
              else
                  call mpistop('No LES model other than "SM" implemented for 5pts stencil', 0)
              endif
           else
              flux_visc => flux_visc_5pts
           endif
        endif
     endif
  case default
     call mpistop('Assign flux_visc: bad choice, only values [0,2,4]!', 0)
  end select

  ! Define one-sided/two-sided or 2D/3D communications
  ! ==================================================
  ! for inviscid fluxes [conservative variables rho_n,rhou_n,...]
  ! -------------------
  if (is_comm_onesided) then
     communication => communication1
  else
     if (is_2D) then
        communication => communication2d
     else
        communication => communication3d
     endif
  endif

  if (is_2D) then
     communication_ => communication_2d
  else
     communication_ => communication_3d
  endif
  
  !!!!!!! ONE-SIDED COMM VISC NOT WRITTEN YET 
  
  ! for viscous fluxes [velocity derivatives]
  ! ------------------
  if (iorder_visc==0) then
     communication_v => voidp
  else
     if (is_2D) then
        communication_v => communication_2dv
     else
        communication_v => communication_3dv
     endif
  endif

  ! Define wall boundary conditions
  ! ===============================
  if (is_wall2) then
     if (iorder_visc==0) then
        if (is_curv3) then
           bc_wall_imin => bc_wall_imin_slip_c3
           bc_wall_imax => bc_wall_imax_slip_c3
           bc_wall_jmin => bc_wall_jmin_slip_c3
           bc_wall_jmax => bc_wall_jmax_slip_c3
           bc_wall_kmin => bc_wall_kmin_slip_c3
           bc_wall_kmax => bc_wall_kmax_slip_c3
        else
           if (is_curv) then
              bc_wall_imin => bc_wall_imin_slip_c
              bc_wall_imax => bc_wall_imax_slip_c
              bc_wall_jmin => bc_wall_jmin_slip_c
              bc_wall_jmax => bc_wall_jmax_slip_c
           else
              bc_wall_imin => bc_wall_imin_slip
              bc_wall_imax => bc_wall_imax_slip
              bc_wall_jmin => bc_wall_jmin_slip
              bc_wall_jmax => bc_wall_jmax_slip
           endif
           bc_wall_kmin => bc_wall_kmin_slip
           bc_wall_kmax => bc_wall_kmax_slip
        endif
     else
        if (is_adiab) then
           if (is_curv3) then
              if (is_slip(1,1)) then
                 bc_wall_imin => bc_wall_imin_slip_c3
              else
                 bc_wall_imin => bc_wall_imin_adiabatic_c3
              endif
              if (is_slip(1,2)) then
                 bc_wall_imax => bc_wall_imax_slip_c3
              else
                 bc_wall_imax => bc_wall_imax_adiabatic_c3
              endif
              if (is_slip(2,1)) then
                 bc_wall_jmin => bc_wall_jmin_slip_c3
              else
                 bc_wall_jmin => bc_wall_jmin_adiabatic_c3
              endif
              if (is_slip(2,2)) then
                 bc_wall_jmax => bc_wall_jmax_slip_c3
              else
                 bc_wall_jmax => bc_wall_jmax_adiabatic_c3
              endif
              if (is_slip(3,1)) then
                 bc_wall_kmin => bc_wall_kmin_slip_c3
              else
                 bc_wall_kmin => bc_wall_kmin_adiabatic_c3
              endif
              if (is_slip(3,2)) then
                 bc_wall_kmax => bc_wall_kmax_slip_c3
              else
                 bc_wall_kmax => bc_wall_kmax_adiabatic_c3
              endif
           else
              if (is_curv) then
                if (is_slip(1,1)) then
                    bc_wall_imin => bc_wall_imin_slip_c
                else
                    bc_wall_imin => bc_wall_imin_adiabatic_c
                endif
                if (is_slip(1,2)) then
                    bc_wall_imax => bc_wall_imax_slip_c
                else
                    bc_wall_imax => bc_wall_imax_adiabatic_c
                endif
                if (is_slip(2,1)) then
                    bc_wall_jmin => bc_wall_jmin_slip_c
                else
                    bc_wall_jmin => bc_wall_jmin_adiabatic_c
                endif
                if (is_slip(2,2)) then
                    bc_wall_jmax => bc_wall_jmax_slip_c
                else
                    bc_wall_jmax => bc_wall_jmax_adiabatic_c
                endif
              else
                if (is_slip(1,1)) then
                    bc_wall_imin => bc_wall_imin_slip
                else
                    bc_wall_imin => bc_wall_imin_adiabatic
                endif
                if (is_slip(1,2)) then
                    bc_wall_imax => bc_wall_imax_slip
                else
                    bc_wall_imax => bc_wall_imax_adiabatic
                endif
                if (is_slip(2,1)) then
                    bc_wall_jmin => bc_wall_jmin_slip
                else
                    bc_wall_jmin => bc_wall_jmin_adiabatic
                endif
                if (is_slip(2,2)) then
                    bc_wall_jmax => bc_wall_jmax_slip
                else
                    bc_wall_jmax => bc_wall_jmax_adiabatic
                endif
              endif
              bc_wall_kmin => bc_wall_kmin_adiabatic
              bc_wall_kmax => bc_wall_kmax_adiabatic
           endif
        else
           bc_wall_imin => bc_wall_imin_isotherm
           bc_wall_imax => bc_wall_imax_isotherm
           bc_wall_jmin => bc_wall_jmin_isotherm
           bc_wall_jmax => bc_wall_jmax_isotherm
           bc_wall_kmin => bc_wall_kmin_isotherm
           bc_wall_kmax => bc_wall_kmax_isotherm
        endif
     endif
  else
     if (is_curv3) then
        if (is_adiab) then
           bc_wall_imin => bc_wall_imin_adiabatic_dpdn_c3
           bc_wall_imax => bc_wall_imax_adiabatic_dpdn_c3
           bc_wall_jmin => bc_wall_jmin_adiabatic_dpdn_c3
           bc_wall_jmax => bc_wall_jmax_adiabatic_dpdn_c3
           bc_wall_kmin => bc_wall_kmin_adiabatic_dpdn_c3
           bc_wall_kmax => bc_wall_kmax_adiabatic_dpdn_c3
        else
           bc_wall_imin => bc_wall_imin_isotherm_dpdn_c3
           bc_wall_imax => bc_wall_imax_isotherm_dpdn_c3
           bc_wall_jmin => bc_wall_jmin_isotherm_dpdn_c3
           bc_wall_jmax => bc_wall_jmax_isotherm_dpdn_c3
           bc_wall_kmin => bc_wall_kmin_isotherm_dpdn_c3
           bc_wall_kmax => bc_wall_kmax_isotherm_dpdn_c3
        endif
     else
        if (is_curv) then
           if (is_adiab) then
              bc_wall_imin => bc_wall_imin_adiabatic_dpdn_c
              bc_wall_imax => bc_wall_imax_adiabatic_dpdn_c
              bc_wall_jmin => bc_wall_jmin_adiabatic_dpdn_c
              bc_wall_jmax => bc_wall_jmax_adiabatic_dpdn_c
           else
              bc_wall_imin => bc_wall_imin_isotherm_dpdn_c
              bc_wall_imax => bc_wall_imax_isotherm_dpdn_c
              bc_wall_jmin => bc_wall_jmin_isotherm_dpdn_c
              bc_wall_jmax => bc_wall_jmax_isotherm_dpdn_c
           endif
        else
           if (is_adiab) then
              bc_wall_imin => bc_wall_imin_adiabatic_dpdn
              bc_wall_imax => bc_wall_imax_adiabatic_dpdn
              bc_wall_jmin => bc_wall_jmin_adiabatic_dpdn
              bc_wall_jmax => bc_wall_jmax_adiabatic_dpdn
           else
              bc_wall_imin => bc_wall_imin_isotherm_dpdn
              bc_wall_imax => bc_wall_imax_isotherm_dpdn
              bc_wall_jmin => bc_wall_jmin_isotherm_dpdn
              bc_wall_jmax => bc_wall_jmax_isotherm_dpdn
           endif
        endif
        if (is_adiab) then
           bc_wall_kmin => bc_wall_kmin_adiabatic_dpdn
           bc_wall_kmax => bc_wall_kmax_adiabatic_dpdn
        else
           bc_wall_kmin => bc_wall_kmin_isotherm_dpdn
           bc_wall_kmax => bc_wall_kmax_isotherm_dpdn
        endif
     endif
  endif

  ! Angular periodicity boundary conditions
  ! =======================================
  ! Nota: rotation axis should be ox
  bc_angular_periodicity => voidp
  bc_angular_periodicity_n => voidp
  bc_angular_periodicity_grad=> voidp

  if (theta_period.ne.0.0_wp) then
     bc_angular_periodicity => bc_angular_period
     bc_angular_periodicity_n => bc_angular_period_n
     if ((iorder_visc.ne.0).or.(idepart==POST_PROCESSING)) &
          bc_angular_periodicity_grad => bc_angular_period_grad
  endif

  ! Define Tam & Dong's boundary conditions
  ! =======================================
  if (is_TamDong3D) then
     if (is_curv3) then
        bc_TD_imin => bc_TD3d_imin_c3
        !bc_TD_imin => bc_TD1d_imin_c3 ! only test XG (save as _dev)
        bc_TD_imax => bc_TD3d_imax_c3
        bc_TD_jmin => bc_TD3d_jmin_c3
        bc_TD_jmax => bc_TD3d_jmax_c3
        bc_TD_kmin => bc_TD3d_kmin_c3
        bc_TD_kmax => bc_TD3d_kmax_c3
        ! edges along z
        bc_TD_imin_jmin => bc_TD3d_imin_jmin_c3
        bc_TD_imin_jmax => bc_TD3d_imin_jmax_c3
        !bc_TD_imin_jmin => bc_TD1d_imin_jmin_c3 ! only test XG (save as _dev)
        !bc_TD_imin_jmax => bc_TD1d_imin_jmax_c3 ! only test XG (save as _dev)
        bc_TD_imax_jmin => bc_TD3d_imax_jmin_c3
        bc_TD_imax_jmax => bc_TD3d_imax_jmax_c3
        ! edges along x
        bc_TD_jmin_kmin => bc_TD3d_jmin_kmin_c3
        bc_TD_jmin_kmax => bc_TD3d_jmin_kmax_c3
        bc_TD_jmax_kmin => bc_TD3d_jmax_kmin_c3
        bc_TD_jmax_kmax => bc_TD3d_jmax_kmax_c3
        ! edges along y
        bc_TD_kmin_imin => bc_TD3d_kmin_imin_c3
        bc_TD_kmin_imax => bc_TD3d_kmin_imax_c3
        bc_TD_kmax_imin => bc_TD3d_kmax_imin_c3
        bc_TD_kmax_imax => bc_TD3d_kmax_imax_c3
        ! corners
        bc_TD_imin_jmin_kmin => bc_TD3d_imin_jmin_kmin_c3
        bc_TD_imin_jmin_kmax => bc_TD3d_imin_jmin_kmax_c3
        bc_TD_imin_jmax_kmin => bc_TD3d_imin_jmax_kmin_c3
        bc_TD_imin_jmax_kmax => bc_TD3d_imin_jmax_kmax_c3
        bc_TD_imax_jmin_kmin => bc_TD3d_imax_jmin_kmin_c3
        bc_TD_imax_jmin_kmax => bc_TD3d_imax_jmin_kmax_c3
        bc_TD_imax_jmax_kmin => bc_TD3d_imax_jmax_kmin_c3
        bc_TD_imax_jmax_kmax => bc_TD3d_imax_jmax_kmax_c3
     else
        if (is_curv) then
           call mpistop('curvilinear TamDong3D not written yet. Shutting down...',0)
!!$           ! faces
!!$           bc_TD_imin => bc_TD3d_imin_c
!!$           bc_TD_imax => bc_TD3d_imax_c
!!$           bc_TD_jmin => bc_TD3d_jmin_c
!!$           bc_TD_jmax => bc_TD3d_jmax_c
!!$           bc_TD_kmin => bc_TD3d_kmin_c
!!$           bc_TD_kmax => bc_TD3d_kmax_c
!!$           ! edges along z
!!$           bc_TD_imin_jmin => bc_TD3d_imin_jmin_c
!!$           bc_TD_imin_jmax => bc_TD3d_imin_jmax_c
!!$           bc_TD_imax_jmin => bc_TD3d_imax_jmin_c
!!$           bc_TD_imax_jmax => bc_TD3d_imax_jmax_c
!!$           ! edges along x
!!$           bc_TD_imin_jmin => bc_TD3d_jmin_kmin_c
!!$           bc_TD_imin_jmax => bc_TD3d_jmin_kmax_c
!!$           bc_TD_imax_jmin => bc_TD3d_jmax_kmin_c
!!$           bc_TD_imax_jmax => bc_TD3d_jmax_kmax_c
!!$           ! edges along y
!!$           bc_TD_imin_jmin => bc_TD3d_kmin_imin_c
!!$           bc_TD_imin_jmax => bc_TD3d_kmin_imax_c
!!$           bc_TD_imax_jmin => bc_TD3d_kmax_imin_c
!!$           bc_TD_imax_jmax => bc_TD3d_kmax_imax_c
!!$           ! corners
!!$           bc_TD_imin_jmin_kmin => bc_TD3d_imin_jmin_kmin_c
!!$           bc_TD_imin_jmin_kmax => bc_TD3d_imin_jmin_kmax_c
!!$           bc_TD_imin_jmax_kmin => bc_TD3d_imin_jmax_kmin_c
!!$           bc_TD_imin_jmax_kmax => bc_TD3d_imin_jmax_kmax_c
!!$           bc_TD_imax_jmin_kmin => bc_TD3d_imax_jmin_kmin_c
!!$           bc_TD_imax_jmin_kmax => bc_TD3d_imax_jmin_kmax_c
!!$           bc_TD_imax_jmax_kmin => bc_TD3d_imax_jmax_kmin_c
!!$           bc_TD_imax_jmax_kmax => bc_TD3d_imax_jmax_kmax_c
!!$           ! outflow NOT WRITTEN
!!$           bc_TD_outflow_imin => bc_TD3d_imin_c
!!$           bc_TD_outflow_imax => bc_TD3d_imax_c
!!$           !bc_TD_outflow_imax => bc_TD3d_outflow_imax_c
!!$           bc_TD_outflow_jmin => bc_TD3d_jmin_c
!!$           bc_TD_outflow_jmax => bc_TD3d_jmax_c
!!$           bc_TD_outflow_kmin => bc_TD3d_kmin_c
!!$           bc_TD_outflow_kmax => bc_TD3d_kmax_c
        else
           ! faces
           bc_TD_imin => bc_TD3d_imin
           bc_TD_imax => bc_TD3d_imax
           bc_TD_jmin => bc_TD3d_jmin
           bc_TD_jmax => bc_TD3d_jmax
           bc_TD_kmin => bc_TD3d_kmin
           bc_TD_kmax => bc_TD3d_kmax
           ! edges along z
           bc_TD_imin_jmin => bc_TD3d_imin_jmin
           bc_TD_imin_jmax => bc_TD3d_imin_jmax
           bc_TD_imax_jmin => bc_TD3d_imax_jmin
           bc_TD_imax_jmax => bc_TD3d_imax_jmax
           ! edges along x
           bc_TD_jmin_kmin => bc_TD3d_jmin_kmin
           bc_TD_jmin_kmax => bc_TD3d_jmin_kmax
           bc_TD_jmax_kmin => bc_TD3d_jmax_kmin
           bc_TD_jmax_kmax => bc_TD3d_jmax_kmax
           ! edges along y
           bc_TD_kmin_imin => bc_TD3d_kmin_imin
           bc_TD_kmin_imax => bc_TD3d_kmin_imax
           bc_TD_kmax_imin => bc_TD3d_kmax_imin
           bc_TD_kmax_imax => bc_TD3d_kmax_imax
           ! corners
           bc_TD_imin_jmin_kmin => bc_TD3d_imin_jmin_kmin
           bc_TD_imin_jmin_kmax => bc_TD3d_imin_jmin_kmax
           bc_TD_imin_jmax_kmin => bc_TD3d_imin_jmax_kmin
           bc_TD_imin_jmax_kmax => bc_TD3d_imin_jmax_kmax
           bc_TD_imax_jmin_kmin => bc_TD3d_imax_jmin_kmin
           bc_TD_imax_jmin_kmax => bc_TD3d_imax_jmin_kmax
           bc_TD_imax_jmax_kmin => bc_TD3d_imax_jmax_kmin
           bc_TD_imax_jmax_kmax => bc_TD3d_imax_jmax_kmax
           ! outflow NOT WRITTEN
           bc_TD_outflow_imin => bc_TD3d_imin
           bc_TD_outflow_imax => bc_TD3d_imax
           !bc_TD_outflow_imax => bc_TD3d_outflow_imax
           bc_TD_outflow_jmin => bc_TD3d_jmin
           bc_TD_outflow_jmax => bc_TD3d_jmax
           bc_TD_outflow_kmin => bc_TD3d_kmin
           bc_TD_outflow_kmax => bc_TD3d_kmax
        endif
     endif
  else
     if (is_TamDong_1pt) then
        if (is_curv) then
           ! bc_TD_imin => bc_TD2d_1pt_imin_c
           ! bc_TD_imax => bc_TD2d_1pt_imax_c
           ! bc_TD_jmin => bc_TD2d_1pt_jmin_c
           ! bc_TD_jmax => bc_TD2d_1pt_jmax_c
           ! bc_TD_imin_jmin => bc_TD2d_1pt_imin_jmin_c
           ! bc_TD_imin_jmax => bc_TD2d_1pt_imin_jmax_c
           ! bc_TD_imax_jmin => bc_TD2d_1pt_imax_jmin_c
           ! bc_TD_imax_jmax => bc_TD2d_1pt_imax_jmax_c
           ! bc_TD_outflow_imin => bc_TD2d_1pt_outflow_imin_c
           ! bc_TD_outflow_imax => bc_TD2d_1pt_outflow_imax_c
           ! bc_TD_outflow_jmin => bc_TD2d_1pt_outflow_jmin_c
           ! bc_TD_outflow_jmax => bc_TD2d_1pt_outflow_jmax_c
           call mpistop('Curvilinear TamDong BCs on 1 point not implemented yet. Shutting down...',0)
        else
           bc_TD_imin => bc_TD2d_1pt_imin
           bc_TD_imax => bc_TD2d_1pt_imax
           bc_TD_jmin => bc_TD2d_1pt_jmin
           bc_TD_jmax => bc_TD2d_1pt_jmax
           if (iproc==0) print *,'~> bc_TD2d_outflow_imax on 1 point not implemented yet'
           bc_TD_outflow_imax => bc_TD2d_outflow_imax
           bc_TD_imin_jmin => bc_TD2d_imin_jmin
           bc_TD_imin_jmax => bc_TD2d_imin_jmax
           bc_TD_imax_jmin => bc_TD2d_imax_jmin
           bc_TD_imax_jmax => bc_TD2d_imax_jmax
           bc_TD_outflow_imin => bc_TD2d_1pt_imin
           bc_TD_outflow_jmin => bc_TD2d_1pt_jmin
           bc_TD_outflow_jmax => bc_TD2d_1pt_jmax
        endif
     else
        if (is_curv) then
           bc_TD_imin => bc_TD2d_imin_c!_SBP4
           bc_TD_imax => bc_TD2d_imax_c!_SBP4
           bc_TD_jmin => bc_TD2d_jmin_c!_SBP4
           bc_TD_jmax => bc_TD2d_jmax_c!_SBP4
           bc_TD_imin_jmin => bc_TD2d_imin_jmin_c!_SBP4
           bc_TD_imin_jmax => bc_TD2d_imin_jmax_c!_SBP4
           bc_TD_imax_jmin => bc_TD2d_imax_jmin_c!_SBP4
           bc_TD_imax_jmax => bc_TD2d_imax_jmax_c!_SBP4
           bc_TD_outflow_imin => bc_TD2d_outflow_imin_c
           bc_TD_outflow_imax => bc_TD2d_outflow_imax_c
           bc_TD_outflow_jmin => bc_TD2d_outflow_jmin_c
           bc_TD_outflow_jmax => bc_TD2d_outflow_jmax_c
        else
           bc_TD_imin => bc_TD2d_imin!_SBP4
           bc_TD_imax => bc_TD2d_imax!_SBP4
           bc_TD_jmin => bc_TD2d_jmin!_SBP4
           bc_TD_jmax => bc_TD2d_jmax!_SBP4
           bc_TD_imin_jmin => bc_TD2d_imin_jmin!_SBP4
           bc_TD_imin_jmax => bc_TD2d_imin_jmax!_SBP4
           bc_TD_imax_jmin => bc_TD2d_imax_jmin!_SBP4
           bc_TD_imax_jmax => bc_TD2d_imax_jmax!_SBP4
           bc_TD_outflow_imin => bc_TD2d_imin_SBP4
           bc_TD_outflow_imax => bc_TD2d_outflow_imax
           bc_TD_outflow_jmin => bc_TD2d_jmin_SBP4
           bc_TD_outflow_jmax => bc_TD2d_jmax_SBP4
        endif
     endif
  endif

  ! Define turbomachine boundary conditions
  ! =======================================
!!$  if (TURB) then
!!$     if (eos_type.eq.'pfg') then
        ! Inlet
        bc_inlet_imin  => bc_inflow_imin
        bc_inlet_imax  => bc_inflow_imax
        bc_inlet_jmin  => bc_inflow_jmin
        bc_inlet_jmax  => bc_inflow_jmax
        bc_inlet_kmin  => bc_inflow_kmin
        bc_inlet_kmax  => bc_inflow_kmax
        ! Outlet
        bc_outlet_imin => bc_backpressure_imin
        bc_outlet_imax => bc_backpressure_imax
        bc_outlet_jmin => bc_backpressure_jmin
        bc_outlet_jmax => bc_backpressure_jmax
        bc_outlet_kmin => bc_backpressure_kmin
        bc_outlet_kmax => bc_backpressure_kmax
!!$     elseif (eos_type.eq.'prs') then
!!$        ! Inlet
!!$        bc_inlet_imin  => bc_inlet_prs_imin
!!$        bc_inlet_imax  => bc_inlet_prs_imax
!!$        bc_inlet_jmin  => bc_inlet_prs_jmin
!!$        bc_inlet_jmax  => bc_inlet_prs_jmax
!!$        ! Outlet
!!$        bc_outlet_imin => bc_backpressure_prs_imin
!!$        bc_outlet_imax => bc_backpressure_prs_imax
!!$        bc_outlet_jmin => bc_backpressure_prs_jmin
!!$        bc_outlet_jmax => bc_backpressure_prs_jmax
!!$     endif
!!$  else
!!$  ! Define Riemann based boundary conditions (Euler flow)
!!$  ! ========================================
!!$     if (eos_type.eq.'pfg') then
!!$        ! Inlet
!!$        bc_inlet_imin  => bc_inlet_riemann_pfg_imin
!!$        bc_inlet_imax  => bc_inlet_riemann_pfg_imax
!!$        bc_inlet_jmin  => bc_inlet_riemann_pfg_jmin
!!$        bc_inlet_jmax  => bc_inlet_riemann_pfg_jmax
!!$        ! Outlet
!!$        bc_outlet_imin => bc_outlet_riemann_pfg_imin
!!$        bc_outlet_imax => bc_outlet_riemann_pfg_imax
!!$        bc_outlet_jmin => bc_outlet_riemann_pfg_jmin
!!$        bc_outlet_jmax => bc_outlet_riemann_pfg_jmax
!!$     else
!!$        call mpistop('In/outflow Riemann based conditions not implemented for PRS',0)
!!$     endif
!!$  endif

  ! Same as earlier, only for RANS
  ! ==============================
  if (is_RANS) then

  !******************!
  ! Spalart-Allmaras !
  !******************!
     if (model_RANS.eq.'SA') then
        if(is_transition .and. model_transition == "GRE")then ! GRE: Gamma-Re_theta
                        ! ======= Transition with the Gamma-Re_theta transport model
            start_runge_kutta_rans => start_runge_kutta_SA_transition
            runge_kutta_rans       => runge_kutta_SA_transition
            ! Initialise turbulent variable field
             init_rans => init_SA_transition
    
            ! Define order of Euler and viscous fluxes and derivatives
            if (is_curv3) then
                if (stencil_RANS.eq.5) then
                   flux_euler_rans => flux_euler_5pts_SA_transition_c3
                   grad_rans       => grad_scalar_5pts_c3
                   flux_visc_rans  => flux_visc_5pts_SA_transition_c3
                elseif (stencil_RANS.eq.3) then
                   flux_euler_rans => flux_euler_3pts_SA_transition_c3
                   grad_rans      => grad_scalar_3pts_c3
                   flux_visc_rans => flux_visc_3pts_SA_transition_c3
                endif
            else
               if (stencil_RANS.eq.5) then
                  flux_euler_rans => flux_euler_5pts_SA_transition_c
                  grad_rans       => grad_scalar_5pts_c
                  flux_visc_rans  => flux_visc_5pts_SA_transition_c
               elseif (stencil_RANS.eq.3) then
                  flux_euler_rans => flux_euler_3pts_SA_transition_c
                  grad_rans      => grad_scalar_3pts_c
                  flux_visc_rans => flux_visc_3pts_SA_transition_c
               endif
            endif
            ! Define source term expression
            source_rans => source_SA_transition
    
            ! Define communication
            if (is_2D) then
               communication_grad_rans => communication_2d_grad_rans
               communication_rans      => communic2d
            else
               communication_grad_rans => communication_3d_grad_rans
               communication_rans      => communic3d
            endif
    
            ! Define wall boundary conditions
            bc_wall_imin_rans => bc_wall_imin_SA_transition
            bc_wall_imax_rans => bc_wall_imax_SA_transition
            bc_wall_jmin_rans => bc_wall_jmin_SA_transition
            bc_wall_jmax_rans => bc_wall_jmax_SA_transition
            bc_wall_kmin_rans => bc_wall_kmin_SA_transition
            bc_wall_kmax_rans => bc_wall_kmax_SA_transition
            
            ! Define inlet boundary conditions
            inlet_imin_rans => inlet_imin_SA_transition
            inlet_imax_rans => inlet_imax_SA_transition
            inlet_jmin_rans => inlet_jmin_SA_transition
            inlet_jmax_rans => inlet_jmax_SA_transition
            inlet_kmin_rans => inlet_kmin_SA_transition
            inlet_kmax_rans => inlet_kmax_SA_transition
    
            ! Neumann boundary condition:
            neu_imin_rans => neu_imin_SA_transition
            neu_imax_rans => neu_imax_SA_transition
            neu_jmin_rans => neu_jmin_SA_transition
            neu_jmax_rans => neu_jmax_SA_transition
            neu_kmin_rans => neu_kmin_SA_transition
            neu_kmax_rans => neu_kmax_SA_transition
            ! Define numerical dissipation: DNC, for now only applied to nutild, needs to be changed
            if (is_shock) then
               num_dissip_rans => artvisc_o3_shock_SA
               !num_dissip_rans => artvisc_rus_SA_transition
            else
               num_dissip_rans => artvisc_o3_SA_transition
               !if (stencil_RANS.eq.3) num_dissip_rans => artvisc_rus_SA_transition
            endif
    
            ! Update conservaive variables
            update_var_rans => update_var_SA_transition
        else
            !========= Fully turbulent SA or transitional with algebraic model
            start_runge_kutta_rans => start_runge_kutta_SA
            runge_kutta_rans       => runge_kutta_SA
            ! Initialise turbulent variable field
            if(.not. is_transition)then ! Fully turbulent
                init_rans => init_SA
            elseif(is_transition .and. model_transition == "ALG")then ! Transitional with algebraic model, maybe redundant
                init_rans => init_SA_transition_algebraic
            endif

            ! Define order of Euler and viscous fluxes and derivatives
            if (is_curv3) then
               if (stencil_RANS.eq.5) then
                  flux_euler_rans => flux_euler_5pts_SA_c3
                  grad_rans       => grad_scalar_5pts_c3
                  flux_visc_rans  => flux_visc_5pts_SA_c3
               elseif (stencil_RANS.eq.3) then
                  flux_euler_rans => flux_euler_3pts_SA_c3
                  grad_rans      => grad_scalar_3pts_c3
                  flux_visc_rans => flux_visc_3pts_SA_c3
               endif
            else
               if (stencil_RANS.eq.5) then
                  flux_euler_rans => flux_euler_5pts_SA_c
                  grad_rans       => grad_scalar_5pts_c
                  flux_visc_rans  => flux_visc_5pts_SA_c
               elseif (stencil_RANS.eq.3) then
                  flux_euler_rans => flux_euler_3pts_SA_c
                  grad_rans      => grad_scalar_3pts_c
                  flux_visc_rans => flux_visc_3pts_SA_c
               endif
            endif
            ! Define IRS procedures:
            irs2_ngh_i_rans => irs2_ngh_i_SA
            irs2_ngh_j_rans => irs2_ngh_j_SA
            irs2_ngh_k_rans => irs2_ngh_k_SA
            irs4_ngh_i_rans => irs4_ngh_i_SA
            irs4_ngh_j_rans => irs4_ngh_j_SA
            irs4_ngh_k_rans => irs4_ngh_k_SA
            
            ! Define source term expression
            if(.not. is_transition)then ! Fully turbulent
                source_rans => source_SA_criv
            elseif(is_transition .and. model_transition == "ALG")then ! Transitional with algebraic model, maybe redundant
                source_rans => source_SA_transition_algebraic
            endif

            ! Define communication
            if (is_2D) then
               communication_grad_rans => communication_2d_grad_rans
               communication_rans      => communic2d
            else
               communication_grad_rans => communication_3d_grad_rans
               communication_rans      => communic3d
            endif

            ! Define wall boundary conditions
            bc_wall_imin_rans => bc_wall_imin_SA
            bc_wall_imax_rans => bc_wall_imax_SA
            bc_wall_jmin_rans => bc_wall_jmin_SA
            bc_wall_jmax_rans => bc_wall_jmax_SA
            bc_wall_kmin_rans => bc_wall_kmin_SA
            bc_wall_kmax_rans => bc_wall_kmax_SA

            ! Define inlet boundary conditions
            inlet_imin_rans => inlet_imin_sa
            inlet_imax_rans => inlet_imax_sa
            inlet_jmin_rans => inlet_jmin_sa
            inlet_jmax_rans => inlet_jmax_sa
            inlet_kmin_rans => inlet_kmin_sa
            inlet_kmax_rans => inlet_kmax_sa
    
            ! Neumann boundary condition:
            neu_imin_rans => neu_imin_sa
            neu_imax_rans => neu_imax_sa
            neu_jmin_rans => neu_jmin_sa
            neu_jmax_rans => neu_jmax_sa
            neu_kmin_rans => neu_kmin_sa
            neu_kmax_rans => neu_kmax_sa

           ! Define numerical dissipation
        ! ----------------------------
        ! Selective filtering
        if (is_SF) then
           if (is_shock) then
              num_dissip_rans => filter_o4_shock_SA
              if (stencil_RANS.eq.3) num_dissip_rans => filter_o2_shock_SA
           else
              num_dissip_rans => filtering_5pts_SA
              if (stencil_RANS.eq.3) num_dissip_rans => filtering_3pts_SA
           endif
        else
        ! Artificial viscosity DNC
           if (is_shock) then
              num_dissip_rans => artvisc_o3_shock_SA
              if (stencil_RANS.eq.3) num_dissip_rans => artvisc_rus_SA
           else
              num_dissip_rans => artvisc_o3_SA
              if (stencil_RANS.eq.3) num_dissip_rans => artvisc_rus_SA
           endif
        endif

        ! Updating procedure
        ! ------------------
        ! Temporal update
        update_var_rans => update_var_SA
        ! Dissipation update
        if (.not.is_dissip_in_increments) then
           if (is_SF) then
              update_varn_rans => update_varn_filter_SA
           else
              update_varn_rans => update_varn_artvisc_SA
           endif
        endif

        endif
       ! Define turb model length scale
       ! ------------------------------
       if (.not.is_DES) then
          length_scale_rans => RANS_length_scale
       elseif (model_DES.eq.'DES97') then
          if (is_SLA) then
             call mpistop('SLA option cannot be activated with DES97', 0)
          else
             ! Not to be used, just for validation test case
             length_scale_rans => DES97_length_scale
          endif
       elseif (model_DES.eq.'DES') then
          ! AB: /!\ DES97 is not DES. DES97 kept because original paper with test case
          !     I believe what we should call DES is in:
          !     Shur et al., 1999, "Detached-eddy simulation of an airfoil at high angle of attack."
          !     Maybe it is just the correct calibration of coefficient c_DES in the paper ?
          !     Needs to be verified at least !
          call mpistop('DES approach not implemented yet', 0)
       elseif (model_DES.eq.'DDES') then
          if (is_SLA) then
             length_scale_rans => DDES_SLA_length_scale
          else
             length_scale_rans => DDES_length_scale
          endif
       elseif (model_DES.eq.'IDDES') then
          if (is_SLA) then
             call mpistop('SLA option not implemented yet with IDDES', 0)
          else
             length_scale_rans => IDDES_length_scale
          endif
       else
          call mpistop('invalid simulation option', 0)
       endif

     endif
  endif

  ! Assign for wall-model
  ! =====================
  if (is_wall_model) then
     if (wm_model_type.eq."ODE") then
       bc_wm_jmin => bc_wm_ODE_jmin
       bc_wm_jmax => bc_wm_ODE_jmax
     else if (wm_model_type.eq."ALG") then
       bc_wm_jmin => bc_wm_alg_jmin
       bc_wm_jmax => bc_wm_alg_jmax
     else if (wm_model_type.eq."ALG_WALE") then 
        bc_wm_jmin => bc_wm_wale_jmin
        bc_wm_jmax => bc_wm_wale_jmax 
     else
       call mpistop('Wrong choice for wm_model_type...',0)
     endif
  endif

end subroutine assign_procedures
