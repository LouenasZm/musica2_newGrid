!===============================================================================
subroutine runge_kutta
!===============================================================================
  !> Runge-Kutta integration (compute viscous terms at the last substep)
!===============================================================================
  use mod_flow          ! for: conservative variables
  use mod_time          ! for: irk,nrk
  use mod_constant      ! for: is_mean0, is_residue
  use mod_eos           ! for: pointer procedures primitives
  use mod_interface     ! for: pointer procedures flux_euler,flux_visc,filter,communication ...
  use mod_bc_apply      ! for: bc_wall, bc_non_reflecting
  use mod_routines      ! for: source, residuals
  use mod_sponge        ! for: sponge zone
  use mod_flow0         ! for: mean0
  use mod_bc_periodicity! for: bc_angular_periodicity
  implicit none
  ! ---------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------
  
  ! First  RK steps
  ! ===============
  do irk=1,nrk-1
     if (irk.gt.1) then ! Nota: these routines were already
                        ! called at the end of the previous RK
        ! apply wall boundary conditions
        call bc_wall
        ! compute primitive variables (ui,p,T)
        call primitives
        ! apply non-reflecting boundary conditions
        call bc_non_reflecting
     endif

     ! add source
     if (is_src) call source

     ! compute Eulerian fluxes
     call flux_euler

     ! if selectionned, application of implicit residual smoothing
     if (is_irs) then
        call irs
     else
        ! update conservative variables
        call update_var
     endif

     ! communicate conservative variables (MPI)
     call communication
     ! impose angular periodicity
     call bc_angular_periodicity_n
  enddo
  
  ! Last RK step
  ! ============
  irk=nrk

  ! apply wall boundary conditions
  call bc_wall
  ! compute primitive variables (ui,p,T)
  call primitives
  ! advance temporal mean of primitive variables
  if (is_mean0) call mean0
  ! apply non-reflecting boundary conditions
  call bc_non_reflecting

  ! add source
  if (is_src) call source

  ! compute Eulerian fluxes
  call flux_euler

  ! compute viscous fluxes + viscosity
  call primitives_visc!(rho,rhou,rhov,rhow,rhoe)
  call grad_vel
  call grad_T
  call communication_v
  ! impose angular periodicity
  call bc_angular_periodicity_grad

  ! compute RANS field -> Boussinesq
  irk=1
  if (is_RANS) call runge_kutta_rans
  irk=nrk

  ! compute viscous fluxes
  call flux_visc
  
  ! add numerical dissipation to increment
  if (is_dissip_in_increments) then
     call num_dissip
     ! apply sponge zone
     if (isponge(iproc)) then
        if (is_curv3) then
           call apply_sponge_3d
        else
           call apply_sponge
        endif
     endif
  endif

  ! if selectionned, application of implicit residual smoothing
  if (is_irs) then
     call irs
  else
     ! update conservative variables
     call update_var
  end if

  call communication
  ! impose angular periodicity
  call bc_angular_periodicity_n

  ! compute numerical dissipation terms
  if (.not.is_dissip_in_increments) then
     call num_dissip
     ! apply sponge zone
     if (isponge(iproc)) then
        if (is_curv3) then
           call apply_sponge_3d
        else
           call apply_sponge
        endif
     endif
     call update_varn
     ! communicate conservative variables
     call communication
     ! impose angular periodicity
     call bc_angular_periodicity_n
  endif
  
  !if (CHAN) call forcing_rho

  irk=1

  ! apply wall boundary conditions
  call bc_wall
  ! compute primitive variables (ui,p,T)
  call primitives
  ! apply non-reflecting boundary conditions
  call bc_non_reflecting

  ! compute residuals (for steady solutions)
  if ((is_residue).and.(mod(ntime,nprint).eq.0)) call residuals

  ! update solution for next iteration
  rho=rho_n
  rhou=rhou_n
  rhov=rhov_n
  rhow=rhow_n
  rhoe=rhoe_n

end subroutine runge_kutta

!===============================================================================
subroutine runge_kutta_
!===============================================================================
  !> Runge-Kutta integration (compute viscous terms at the last substep)
!===============================================================================
  use mod_flow          ! for: conservative variables
  use mod_time          ! for: irk,nrk
  use mod_constant      ! for: is_mean0, is_residue
  use mod_eos           ! for: pointer procedures primitives
  use mod_interface     ! for: pointer procedures flux_euler,flux_visc,filter,communication ...
  use mod_bc_apply      ! for: bc_wall, bc_non_reflecting
  use mod_routines      ! for: source, residuals
  use mod_sponge        ! for: sponge zone
  use mod_flow0         ! for: mean0
  use mod_bc_periodicity! for: bc_angular_periodicity
  implicit none
  ! ---------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------

  ! First  RK steps
  ! ===============
  do irk=1,nrk-1
     if (irk.gt.1) then ! Nota: these routines were already
                        ! called at the end of the previous RK
        ! apply wall boundary conditions
        call bc_wall
        ! compute primitive variables (ui,p,T)
        call primitives
        ! apply non-reflecting boundary conditions
        call bc_non_reflecting
     endif

     ! add source
     if (is_src) call source

     ! compute Eulerian fluxes
     call flux_euler

     ! compute viscous fluxes + viscosity
     call grad_vel
     call grad_T
     call communication_v
     ! impose angular periodicity
     call bc_angular_periodicity_grad
     
     call flux_visc

     ! if selectionned, application of implicit residual smoothing
     if (is_irs) then
        call irs
     else
        ! update conservative variables
        call update_var
        !call update_var2N
     end if

     ! communicate conservative variables (MPI)
     call communication
     ! impose angular periodicity
     call bc_angular_periodicity_n

     ! compute numerical dissipation terms
     if (.not.is_dissip_in_increments) then
        call num_dissip
        ! apply sponge zone
        if (isponge(iproc)) call apply_sponge
        call update_varn
        ! communicate conservative variables
        call communication
        ! impose angular periodicity
        call bc_angular_periodicity_n
     endif
  enddo

  ! Last RK step
  ! ============
  irk=nrk

  ! apply wall boundary conditions
  call bc_wall
  ! compute primitive variables (ui,p,T)
  call primitives
  ! advance temporal mean of primitive variables
  if (is_mean0) call mean0
  ! apply non-reflecting boundary conditions
  call bc_non_reflecting

  ! add source
  if (is_src) call source

  ! compute Eulerian fluxes
  call flux_euler

  ! compute viscous fluxes + viscosity
  call grad_vel
  call grad_T
  call communication_v
  ! impose angular periodicity
  call bc_angular_periodicity_grad
       
  ! compute RANS field -> Boussinesq
  irk=1
  if (is_RANS) call runge_kutta_rans
  irk=nrk
  
  call flux_visc

  ! add numerical dissipation to increment
  if (is_dissip_in_increments) then
     call num_dissip
     ! apply sponge zone
     if (isponge(iproc)) call apply_sponge
  endif
  
  ! if selectionned, application of implicit residual smoothing
  if (is_irs) then
     call irs
  else
     ! update conservative variables
     call update_var
     !call update_var2N
  end if

  call communication
  ! impose angular periodicity
  call bc_angular_periodicity_n
  
  ! compute numerical dissipation terms
  if (.not.is_dissip_in_increments) then
     call num_dissip
     ! apply sponge zone
     if (isponge(iproc)) call apply_sponge
     call update_varn
     ! communicate conservative variables
     call communication
     ! impose angular periodicity
     call bc_angular_periodicity_n
  endif

  !if (CHAN) call forcing_rho

  irk=1
  ! apply wall boundary conditions
  call bc_wall
  ! compute primitive variables (ui,p,T)
  call primitives
  ! apply non-reflecting boundary conditions
  call bc_non_reflecting

  ! compute residuals (for steady solutions)
  if ((is_residue).and.(mod(ntime,nprint).eq.0)) call residuals

  ! update solution for next iteration
  rho=rho_n
  rhou=rhou_n
  rhov=rhov_n
  rhow=rhow_n
  rhoe=rhoe_n

end subroutine runge_kutta_

!===============================================================================
subroutine start_runge_kutta
!===============================================================================
  !> Start first step of Runge-Kutta
!===============================================================================
  use mod_flow      ! for: conservative variables
  use mod_time      ! for: irk
  use mod_eos       ! for: pointer procedures primitives
  use mod_interface ! for: pointer procedures grad_vel grad_T, ..
  use mod_bc_apply  ! for: bc_wall, bc_non_reflecting
  implicit none
  ! ---------------------------------------------------------------------------
  
  ! Start first step of Runge-Kutta
  ! ===============================
  irk=1
  ! apply wall boundary conditions
  call bc_wall
  ! compute primitive variables (ui,p,T,c,...)
  !call primitives_visc(rho_n,rhou_n,rhov_n,rhow_n,rhoe_n)
  call primitives
  ! compute velocity & temperature gradients
  call grad_vel
  call grad_T
  ! apply non-reflecting boundary conditions
  call bc_non_reflecting
  
end subroutine start_runge_kutta
