!===============================================================
subroutine update_var_SA
!===============================================================
! Update conservative nutilde at next RK step
!===============================================================
  use mod_flow
  use mod_time
  use mod_rans
  implicit none
  ! ---------------------------------------------------------------------------
  integer  :: i,j,k
  real(wp) :: alpha,large_value,low_value
  ! ---------------------------------------------------------------------------

  ! Update variable
  ! =============== 
  large_value=30000.0_wp
  low_value  =1.0E-20

  if (is_dtlocal) then
     do k=1,nz
        do j=1,ny
           do i=1,nx
!               uvar(i,j,k,2) = nutil_n(i,j,k)
              nutil_n(i,j,k) = nutil(i,j,k)-dt_local(i,j,k)*crk(irk)*(Knutil(i,j,k)-Sterm(i,j,k))
              !if (irk.eq.nrk.and.nutil_n(i,j,k)/mu_ref.gt.large_value) then
              !   !print*,nutil_n(i,j,k),'mut too high, i,j,k=',i,j,k,'iproc, nbl=',iproc,nbl
              !   nutil_n(i,j,k)=mu_ref*large_value
              !endif
              if (nutil_n(i,j,k).lt.0.0_wp) nutil_n(i,j,k)=low_value
           enddo
        enddo
     enddo
  else

     alpha=deltat*crk(irk)
     do k=1,nz
        do j=1,ny
           do i=1,nx
!               uvar(i,j,k,2) = nutil_n(i,j,k)
              nutil_n(i,j,k) = nutil(i,j,k)-alpha*(Knutil(i,j,k)-Sterm(i,j,k))
              !if (irk.eq.nrk.and.nutil_n(i,j,k)/mu_ref.gt.large_value) then
              !   print*,nutil_n(i,j,k),'mut too high, i,j,k=',i,j,k,'iproc, nbl=',iproc,nbl
              !   nutil_n(i,j,k)=mu_ref*large_value
              !endif
              if (nutil_n(i,j,k).lt.0.0_wp) nutil_n(i,j,k)=low_value
           enddo
        enddo
     enddo
  endif

end subroutine update_var_SA

!========================================================================
subroutine update_varn_artvisc_SA
!========================================================================
! Update conservative nutilde after application of numerical dissipation
!========================================================================
  use mod_flow
  use mod_time
  use mod_rans
  implicit none
  ! ---------------------------------------------------------------------------
  integer  :: i,j,k
  real(wp) :: large_value,low_value
  ! ---------------------------------------------------------------------------

  ! Update variable with added dissipation
  ! =============== 

  if (is_dtlocal) then
     do k=1,nz
        do j=1,ny
           do i=1,nx
              nutil_n(i,j,k) = nutil(i,j,k)-dt_local(i,j,k)*Knutil(i,j,k)
              if (nutil_n(i,j,k).lt.0.0_wp) nutil_n(i,j,k)=low_value
           enddo
        enddo
     enddo
  else
     do k=1,nz
        do j=1,ny
           do i=1,nx
              nutil_n(i,j,k) = nutil(i,j,k)-deltat*Knutil(i,j,k)
              if (nutil_n(i,j,k).lt.0.0_wp) nutil_n(i,j,k)=low_value
           enddo
        enddo
     enddo
  endif

end subroutine update_varn_artvisc_SA

!===============================================================================
subroutine update_varn_filter_SA
!===============================================================================
  !> Update conservative variables after application of numerical dissipation
  !> * version to be used with filtering *
!===============================================================================
  use mod_flow         ! for: conservative variables
  use mod_rans         ! for: conservative variables
  implicit none
  ! ----------------------------------------------------------------------------
  integer  :: i,j,k
  ! ----------------------------------------------------------------------------

  ! Update variables
  ! ================

  do k=1,nz
     do j=1,ny
        do i=1,nx
           nutil_n(i,j,k)  = nutil_n(i,j,k)  + Knutil(i,j,k)
        enddo
     enddo
  enddo

end subroutine update_varn_filter_SA

!===============================================================
! Update conservative variables to have dwi before application
!  of IRS operator
!===============================================================
subroutine update_var_in_dw_SA
  use mod_flow
  use mod_time
  use mod_forcing_bulk
  implicit none
  ! ---------------------------------------------------------------------------
  integer  :: i,j,k
  real(wp) :: alpha
  ! ---------------------------------------------------------------------------

  ! Update variables
  ! ================
  alpha=deltat*crk(irk)

  !do k=ndz0,nfz0
  do k=1,nz
     do j=1,ny
        do i=1,nx
           Krho(i,j,k)  = - alpha*Krho(i,j,k)
           Krhou(i,j,k) = - alpha*Krhou(i,j,k)
           Krhov(i,j,k) = - alpha*Krhov(i,j,k)
           Krhow(i,j,k) = - alpha*Krhow(i,j,k)
           Krhoe(i,j,k) = - alpha*Krhoe(i,j,k)
        enddo
     enddo
  enddo

  if (is_2D) Krhow = 0.0_wp

  ! Enforce mass-flow rate
  ! ======================
  if (is_forcing_bulk) then
     ! compute forcing term
     if (irk==nrk) then
        if (PHILL) then
           call forcing_rhou_c
        else
           call forcing_rhou
        endif
     endif
     ! update increment
     alpha=-alpha*forc_rhou
     do k=1,nz
        do j=1,ny
           do i=1,nx
              Krhou(i,j,k) = Krhou(i,j,k) + alpha
              Krhoe(i,j,k) = Krhoe(i,j,k) + alpha*uu(i,j,k)
           enddo
        enddo
     enddo
  endif

end subroutine update_var_in_dw_SA

!===============================================================================
subroutine update_var_of_dw_SA
!===============================================================================
  !> Update conservative turbulent variables after application of IRS operator
!===============================================================================
  use mod_flow         ! for: conservative variables
  use mod_time         ! for: irk,crk,deltat
  use mod_rans
  implicit none
  ! ----------------------------------------------------------------------------
  integer  :: i,j,k
  real(wp) :: alpha
  ! ----------------------------------------------------------------------------

  ! Update variables
  ! ================
  if (is_dtlocal) then
     do k=1,nz
        do j=1,ny
           do i=1,nx
              nutil_n(i,j,k) = nutil(i,j,k)-dt_local(i,j,k)*crk(irk)*(Knutil(i,j,k)-Sterm(i,j,k))
           enddo
        enddo
     enddo
  else
     alpha=-deltat*crk(irk)

     do k=1,nz
        do j=1,ny
           do i=1,nx
              nutil_n(i,j,k) = nutil(i,j,k)-alpha*(Knutil(i,j,k)-Sterm(i,j,k))
           enddo
        enddo
     enddo
  endif

  if (is_2D) rhow_n=0.0_wp

end subroutine update_var_of_dw_SA
