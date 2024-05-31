!================================================================================
module mod_bc_wall_rans
!================================================================================
  !> Module to apply wall Boundary Conditions
!================================================================================
  use mod_flow      ! for: flow variables
  implicit none
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------

contains

  !******************!
  ! Spalart Allmaras !
  !******************!

  !==============================================================================
  subroutine bc_wall_imin_SA
  !==============================================================================
    !> Apply wall Boundary Conditions at imin:
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: j,k
    ! ---------------------------------------------------------------------------

    ! Update wall points
    ! ==================
    do k=1,nz
       do j=1,ny
          ! No turb content at wall
          nutil(1,j,k) = 0.0_wp
          nutil_n(1,j,k) = 0.0_wp
       enddo
    enddo

  end subroutine bc_wall_imin_SA

  !==============================================================================
  subroutine bc_wall_imax_SA
  !==============================================================================
    !> Apply wall Boundary Conditions at imax:
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: j,k
    ! ---------------------------------------------------------------------------

    ! Update wall points
    ! ==================
    do k=1,nz
       do j=1,ny
          ! No turb content at wall
          nutil(nx,j,k) = 0.0_wp
          nutil_n(nx,j,k) = 0.0_wp
       enddo
    enddo

  end subroutine bc_wall_imax_SA

  !==============================================================================
  subroutine bc_wall_jmin_SA
  !==============================================================================
    !> Apply wall Boundary Conditions at jmin:
  !==============================================================================
    use mod_constant ! for is_slip_in
    use mod_mpi      ! for iproc
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,k
    ! ---------------------------------------------------------------------------

    ! Update wall points
    ! ==================
    if ((nob(iproc).eq.1).and.is_slip_in) then ! apply slip wall BC
       do k=1,nz
          do i=1,nx
             nutil(i,1,k) = nutil(i,2,k)
             nutil_n(i,1,k) = nutil_n(i,2,k)
          enddo
       enddo
    else
       do k=1,nz
          do i=1,nx
             ! No turb content at wall
             nutil(i,1,k) = 0.0_wp
             nutil_n(i,1,k) = 0.0_wp
          enddo
       enddo
    endif

  end subroutine bc_wall_jmin_SA
 
  !==============================================================================
  subroutine bc_wall_jmax_SA
  !==============================================================================
    !> Apply wall Boundary Conditions at jmax:
  !==============================================================================
    use mod_constant ! for is_slip_in
    use mod_mpi      ! for iproc
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,k
    ! ---------------------------------------------------------------------------

    ! Update wall points
    ! ==================
    if ((nob(iproc).eq.1).and.is_slip_in) then ! apply slip wall BC
       do k=1,nz
          do i=1,nx
             nutil(i,ny,k) = nutil(i,ny-1,k)
             nutil_n(i,ny,k) = nutil_n(i,ny-1,k)
          enddo
       enddo
    else
       do k=1,nz
          do i=1,nx
             ! No turb content at wall
             nutil(i,ny,k) = 0.0_wp
             nutil_n(i,ny,k) = 0.0_wp
          enddo
       enddo
    endif
  
  end subroutine bc_wall_jmax_SA 

  !==============================================================================
  subroutine bc_wall_kmin_SA
  !==============================================================================
    !> Apply wall Boundary Conditions at kmin:
  !==============================================================================
    use mod_constant ! for is_slip_in
    use mod_mpi      ! for iproc
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j
    ! ---------------------------------------------------------------------------

    ! Update wall points
    ! ==================
    if ((nob(iproc).eq.1).and.is_slip_in) then ! apply slip wall BC
       do j=1,ny
          do i=1,nx
             nutil(i,j,1) = nutil(i,j,2)
             nutil_n(i,j,1) = nutil_n(i,j,2)
          enddo
       enddo
    else
       do j=1,ny
          do i=1,nx
             ! No turb content at wall
             nutil(i,j,1) = 0.0_wp
             nutil_n(i,j,1) = 0.0_wp
          enddo
       enddo
    endif
  
  end subroutine bc_wall_kmin_SA
 
  !==============================================================================
  subroutine bc_wall_kmax_SA
  !==============================================================================
    !> Apply wall Boundary Conditions at kmax:
  !==============================================================================
    use mod_constant ! for is_slip_in
    use mod_mpi      ! for iproc
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j
    ! ---------------------------------------------------------------------------

    ! Update wall points
    ! ==================
    if ((nob(iproc).eq.1).and.is_slip_in) then ! apply slip wall BC
       do j=1,ny
          do i=1,nx
             nutil(i,j,nz) = nutil(i,j,nz-1)
             nutil_n(i,j,nz) = nutil_n(i,j,nz-1)
          enddo
       enddo
    else
       do j=1,ny
          do i=1,nx
             ! No turb content at wall
             nutil(i,j,nz) = 0.0_wp
             nutil_n(i,j,nz) = 0.0_wp
          enddo
       enddo
    endif
  
  end subroutine bc_wall_kmax_SA

end module mod_bc_wall_rans
