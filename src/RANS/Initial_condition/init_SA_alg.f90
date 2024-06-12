!===============================================================================
subroutine init_SA_transition_algebraic
!===============================================================================
  !> Initialise SA model hen coupled with algebraic transitional model-> mut
!===============================================================================
  use mod_flow
  use mod_rans
  use mod_mpi
  use mod_init_flow
  use mod_block
  use mod_time
  use mod_wall_dist
  use mod_interface
  implicit none
  ! ------------------------------------------------------------------------------
  integer :: i,j,k
  ! real(wp) :: radius
  real(wp) :: low_value
  ! ------------------------------------------------------------------------------

!  Nutref = mu_ref/rho_ref*3.0_wp ! fully turbulent
   Nutref = 0.79_wp*mu_ref/rho_ref ! laminar
  if ((idepart.eq.FROM_SCRATCH).and.(ndeb_RANS.gt.0)) then
     nutil  = Nutref
     low_value = 1.0E-20_wp

     ! Initialize with tanh
     do k=1,nz
        do j=1,ny
           do i=1,nx
              nutil(i,j,k) = Nutref*tanh(d(i,j,k)/L_ref*2.0_wp)
           enddo
        enddo
     enddo
     call communication_rans(nutil)

  
  elseif ((ntotal.le.ndeb_RANS).and.(ndeb_RANS.gt.0)) then
     nutil = Nutref
  endif

  nut = nutil

  if (iproc.eq.0) print*, 'Spalart-Allmaras init OK'

end subroutine init_SA_transition_algebraic
