!================================================================================
subroutine primitives_SA_transition
!================================================================================
  !> Computation of primitive variables for perfect gas EoS for RANS equations
!================================================================================
  use mod_flow
  use mod_fluid
  use mod_mpi
  implicit none
  ! ----------------------------------------------------------------------------
  integer :: i,j,k
  ! ----------------------------------------------------------------------------

  ! Compute primitive variable from conservative one
  ! ================================================
  ! for all points + ghost cells
  do k=ndzt,nfzt
     do j=ndyt,nfyt
        do i=ndxt,nfxt
           ! Turbulent kinematic viscosity
           ! -----------------------------
           nutil(i,j,k)=nutil_n(i,j,k)
           ! Intermittency factor and Re_theta
           ! ---------------------------------
           Intermittency(i,j,k)     = rhogamma_n(i,j,k)/rho_n(i,j,k)
           reynolds_theta(i,j,k)    = rhore_theta_n(i,j,k)/rho_n(i,j,k)

            if (intermittency(i,j,k) .lt. 0.0_wp)then
                intermittency(i,j,k) = 0.0_wp 
            elseif(intermittency(i,j,k) .ge. 1.0_wp)then
                intermittency(i,j,k) = 1.0_wp
            endif
        enddo
     enddo
  enddo

end subroutine primitives_SA_transition
