!===============================================================================
subroutine init_SA_transition
!===============================================================================
  !> Initialise SA model -> mut
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
  real(wp) :: radius,low_value, rhogamma_ref
  ! ------------------------------------------------------------------------------
  
  Nutref = 3.0_wp*mu_ref/rho_ref ! fully turbulent
  !Nutref = mu_ref/rho_ref*1e-3_wp ! laminar
  if (idepart.eq.FROM_SCRATCH) then
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
     ! Intialize gamma at 1: to be changed later
     rhogamma_ref = rho_ref*1.0_wp
     do k=1,nz
        do j=1,ny
           do i=1,nx
                rho_gamma(i,j,k)        = rhogamma_ref
                rhogamma_n(i,j,k)       = rhogamma_ref
                intermittency(i,j,k)    = 1.0_wp
           enddo
        enddo
     enddo
     ! init Re_theta with correlation considering du/ds=0
     if(tu_inlet .lt. 1.3_wp)then
        reynolds_theta  = 1173.51_wp - 589.428*tu_inlet + 0.2196_wp/tu_inlet**2
     else
        reynolds_theta  = 331.5_wp * (tu_inlet - 0.5668_wp)**(-0.671_wp)
     endif
     
     ! Initialize Re_theta with analytical solution (check papers first)
     call communication_rans(nutil)
     call communication_rans(intermittency)
     call communication_rans(rho_gamma)
     call communication_rans(rhogamma_n)

  
  elseif (ntotal.le.ndeb_RANS) then
     nutil = Nutref
  endif

  nut = nutil

  if (iproc.eq.0) print*, 'SA-Gamma-Re_theta init OK'

end subroutine init_SA_transition