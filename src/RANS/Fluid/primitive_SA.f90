!================================================================================
subroutine primitives_SA
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
            
        enddo
     enddo
  enddo

end subroutine primitives_SA

