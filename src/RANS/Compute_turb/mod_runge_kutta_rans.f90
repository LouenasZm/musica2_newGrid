
module mod_runge_kutta_rans

  !-------------------------------------------------------------
  use mod_flow            ! for: conservative variables
  use mod_time            ! for: irk,nrk
  use mod_constant        ! for: is_mean0, is_residue
  use mod_eos             ! for: pointer procedures primitives
  use mod_interface       ! for: pointer procedures flux_euler,flux_visc,filter...
  use mod_bc_apply_rans   ! for: bc_wall, bc_non_reflecting
  use mod_neumann_rans    ! for: Neumann bc
  use mod_routines_rans   ! for: source, residuals
  use mod_sponge          ! for: sponge zone
  use mod_flow0           ! for: mean0
  !use mod_source_rans     ! for: source terms
  implicit none
  !-------------------------------------------------------------
  
  contains

   !******************!
   ! Spalart Allmaras !
   !******************!

!===============================================================================
subroutine runge_kutta_SA
!===============================================================================
  !> runge-kutta integration for rans equations 
!===============================================================================
  use mod_flow           ! for: conservative variables
  use mod_time           ! for: irk,nrk
  use mod_constant       ! for: is_mean0, is_residue
  use mod_eos            ! for: pointer procedures primitives
  use mod_interface      ! for: pointer procedures flux_euler,flux_visc,filter...
  use mod_bc_apply_rans  ! for: bc_wall, bc_non_reflecting
  use mod_neumann_rans   ! for: Neumann bc
  use mod_routines_rans  ! for: source, residuals
  use mod_sponge         ! for: sponge zone
  use mod_flow0          ! for: mean0
  implicit none
  ! ---------------------------------------------------------------------------
 
  ! first  rk steps
  ! ===============
  do irk=1,1
  !do irk=1,nrk
     !if (irk.gt.1) then
        ! compute primitive variables (nutilde,k,omega,...)
        call primitives_SA

        ! apply wall and non-reflecting boundary conditions
        call bc_wall_rans
        call bc_rans
     !endif

     ! compute gradients
     call grad_rans(nutil,dnutilx,dnutily,dnutilz) ! --> grad_scalar_Xpts_c(nutil,dnutilx,...)
     call communication_grad_rans(dnutilx,dnutily,dnutilz)

     ! compute the turbulent model length scale (to determine RANS/DES options)
     call length_scale_rans

     ! compute turbulent source terms (transport eq rhs w/o flux_visc)
     call source_SA_criv

     ! compute eulerian fluxes
     call flux_euler_rans
  
     ! compute viscous fluxes
     call flux_visc_rans
    
     ! compute numerical dissipation
     if (is_dissip_in_increments) then
        call num_dissip_rans
     endif
   
     ! if selectionned, application of implicit residual smoothing
     if (is_irs) then
        call irs_rans
     else
        ! update conservative variables
        call update_var_SA
     end if 
     
     ! correct farfield turbulent bc
     call bc_rans
   
     ! communicate conservative variables
     call communication_rans(nutil_n) ! points towards communic2d(nutil_n)
   
  enddo

  ! Apply dissipation outside the temporal increment
  if (.not.is_dissip_in_increments) then
     call num_dissip_rans
     ! add dissipation to variable
     call update_varn_rans
     ! communicate dissipated variable
     call communication_rans(nutil_n)
  endif

  ! Print residuals to screen & file
  if (is_residue.and.mod(ntime,nprint).eq.0) call residuals_SA

  irk=1
  !call primitives_SA
  call bc_wall_rans
  call bc_rans
  ! update solution for next PHYSICAL iteration
  nutil=nutil_n

end subroutine runge_kutta_SA

!===============================================================================
subroutine start_runge_kutta_SA
!===============================================================================
  !> Start first step of Runge-Kutta
!===============================================================================
  use mod_flow           ! for: conservative variables
  use mod_eos            ! for: pointer procedures primitives
  use mod_interface      ! for: pointer procedures grad_velT
  use mod_bc_apply_rans  ! for: bc_wall, bc_non_reflecting
  use mod_mpi
  implicit none
  ! ---------------------------------------------------------------------------
  ! Start first step of Runge-Kutta
  ! ===============================
  call bc_wall_rans
  call bc_rans
  call grad_rans(nutil,dnutilx,dnutily,dnutilz) ! --> grad_scalar_5pts_c(nutil,dnutilx,...)
  call communication_grad_rans(dnutilx,dnutily,dnutilz)
  if (iproc.eq.0) print*, 'RANS init RK OK'

  end subroutine start_runge_kutta_SA

    !*******************!
    ! SA-Gamma-Re_theta !
    !*******************!

    !===============================================================================
    subroutine runge_kutta_SA_transition
    !===============================================================================
      !> runge-kutta integration for rans equations 
    !===============================================================================
      use mod_flow           ! for: conservative variables
      use mod_time           ! for: irk,nrk
      use mod_constant       ! for: is_mean0, is_residue
      use mod_eos            ! for: pointer procedures primitives
      use mod_interface      ! for: pointer procedures flux_euler,flux_visc,filter...
      use mod_bc_apply_rans  ! for: bc_wall, bc_non_reflecting
      use mod_neumann_rans   ! for: Neumann bc
      use mod_routines_rans  ! for: source, residuals
      use mod_sponge         ! for: sponge zone
      use mod_flow0          ! for: mean0
      implicit none
      ! ---------------------------------------------------------------------------
      ! first  rk steps
      ! ===============
        do irk=1,1
        ! if (irk.gt.1) then
            ! compute primitive variables (nutilde,k,omega,...)
            call primitives_SA_transition
            ! apply wall and non-reflecting boundary conditions
            call bc_wall_rans
            call bc_rans
        ! endif
    
         ! compute gradients
         call grad_rans(nutil,dnutilx,dnutily,dnutilz) ! --> grad_scalar_Xpts_c(nutil,dnutilx,...)
         call communication_grad_rans(dnutilx,dnutily,dnutilz)

         call grad_rans(intermittency, dgammax, dgammay, dgammaz)
         call communication_grad_rans(dgammax, dgammay, dgammaz)

         call grad_rans(reynolds_theta, dre_thetax, dre_thetay, dre_thetaz)
         call communication_grad_rans(dre_thetax, dre_thetay, dre_thetaz)
    
         ! compute the turbulent model length scale (to determine RANS/DES options)
         call length_scale_rans
    
         ! compute turbulent source terms (transport eq rhs w/o flux_visc)
         call source_SA_transition
    
         
         ! compute eulerian fluxes
         call flux_euler_rans
         ! compute viscous fluxes
         call flux_visc_rans
        
         ! compute numerical dissipation
         call num_dissip_rans
       
         ! if selectionned, application of implicit residual smoothing
         if (is_irs) then
            call irs_rans
         else
            ! update conservative variables
            call update_var_SA_transition ! -> points towards update_var_SA_transition
         end if 
         
         ! correct farfield turbulent bc
         call bc_rans
       
         ! communicate conservative variables
         call communication_rans(nutil_n) ! points towards communic2d(nutil_n)
       
      enddo
    
      ! Print residuals to screen & file
      if (is_residue.and.mod(ntime,nprint).eq.0) call residuals_SA
    
      irk=1
      call primitives_SA_transition
      call bc_wall_rans
      call bc_rans
      ! update solution for next PHYSICAL iteration
      nutil=nutil_n
    
    end subroutine runge_kutta_SA_transition
    
    !===============================================================================
    subroutine start_runge_kutta_SA_transition
    !===============================================================================
      !> Start first step of Runge-Kutta
    !===============================================================================
      use mod_flow           ! for: conservative variables
      use mod_eos            ! for: pointer procedures primitives
      use mod_interface      ! for: pointer procedures grad_velT
      use mod_bc_apply_rans  ! for: bc_wall, bc_non_reflecting
      use mod_mpi
      implicit none
      ! ---------------------------------------------------------------------------
      ! Start first step of Runge-Kutta
      ! ===============================
      call bc_wall_rans
      call bc_rans
      call grad_rans(nutil,dnutilx,dnutily,dnutilz) ! --> grad_scalar_5pts_c(nutil,dnutilx,...)
      call grad_rans(intermittency, dgammax, dgammay, dgammaz) 
      call grad_rans(reynolds_theta, dre_thetax, dre_thetay, dre_thetaz)

      call communication_grad_rans(dnutilx,dnutily,dnutilz)
      call communication_grad_rans(dgammax,dgammay,dgammaz)
      call communication_grad_rans(dre_thetax,dre_thetay,dre_thetaz)

      if (iproc.eq.0) print*, 'RANS init RK OK'
    
    end subroutine start_runge_kutta_SA_transition
    
end module mod_runge_kutta_rans
