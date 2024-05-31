!===============================================================================
subroutine init_SA
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
  ! real(wp) :: radius
  real(wp) :: low_value
  ! ------------------------------------------------------------------------------

  Nutref = mu_ref/rho_ref*3.0_wp ! fully turbulent
  !Nutref = mu_ref/rho_ref*1e-3_wp ! laminar
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

     !if (.not.CYL) then
   
     !   do k=1,nz
     !      do j=1,ny
     !         do i=1,nx
     !            radius=sqrt(xc(i,j)**2+yc(i,j)**2)/L_ref
     ! 
     !            if (is_bc_wall(2,1).and.j==1) then
     !               nutil(i,1,k)=0.0_wp
     !            elseif (radius<1.5_wp.and.radius>0.5_wp) then
     !               nutil(i,j,k)=nutil(i,j,k)*AS1(radius,1.0_wp,0.5_wp)
     !            endif
     !            
     !            ! set min value above 0
     !            if (nutil(i,j,k).eq.0.0_wp.and.j.gt.1) nutil(i,j,k)=low_value
     ! 
     !         enddo
     !      enddo
     !   enddo
   
     !endif
  
  elseif ((ntotal.le.ndeb_RANS).and.(ndeb_RANS.gt.0)) then
     nutil = Nutref
  endif

  nut = nutil

  if (iproc.eq.0) print*, 'Spalart-Allmaras init OK'

end subroutine init_SA
