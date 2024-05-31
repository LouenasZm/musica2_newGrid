!============================================================================
module mod_neumann_rans
!============================================================================
  !> author: CM
  !> date: March 2022
  !> Homogeneous Neumann boundary condition
!=============================================================================
  use mod_flow
  use mod_rans
  use mod_constant
  use mod_flow
  implicit none
  !---------------------------------------------------------------------------

contains

  !===============================================================================
  subroutine neu_imin_rans
  !===============================================================================
    !> BC: boundary condition at imax (right)
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k,a
    ! ---------------------------------------------------------------------------
    ! Index of left boundary
    ! ======================
    i=1

    do k=1,nz
       do j=1,ny
          nutil_n(i,j,k) = nutil_n(i+1,j,k)
          nutil(i,j,k) = nutil(i+1,j,k)
          ! fill ghost cells
          do a=-ngh_r,-1
             nutil_n(i+a,j,k) = nutil_n(i,j,k)
          enddo
       enddo
    enddo

  end subroutine neu_imin_rans

  !===============================================================================
  subroutine neu_imax_rans
  !===============================================================================
    !> BC: boundary condition at imax (right)
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k,a
    ! ---------------------------------------------------------------------------
    ! Index of right boundary
    ! =======================
    i=nx

    do k=1,nz
       do j=1,ny
          nutil_n(i,j,k) = nutil_n(i-1,j,k)
          nutil(i,j,k) = nutil(i-1,j,k)
          ! fill ghost cells
          do a=1,ngh_r
             nutil_n(i+a,j,k) = nutil_n(i,j,k)
          enddo
       enddo
    enddo

  end subroutine neu_imax_rans

  !===============================================================================
  subroutine neu_jmin_rans
  !===============================================================================
    !> BC: boundary condition at jmin (bottom)
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k,a
    ! ---------------------------------------------------------------------------
    ! Index of bottom boundary
    ! ========================
    j=1

    do k=1,nz
       do i=1,nx
          nutil_n(i,j,k) = nutil_n(i,j+1,k)
          nutil(i,j,k) = nutil(i,j+1,k)
          ! fill ghost cells
          do a=-1,-ngh_r
             nutil_n(i,j+a,k) = nutil_n(i,j,k)
          enddo
       enddo
    enddo

  end subroutine neu_jmin_rans

  !===============================================================================
  subroutine neu_jmax_rans
  !===============================================================================
    !> BC: boundary condition at jmax (top) 
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k,a
    ! ---------------------------------------------------------------------------
    ! Index of top boundary
    ! =====================
    j=ny

    do k=1,nz
       do i=1,nx
          nutil_n(i,j,k) = nutil_n(i,j-1,k)
          nutil(i,j,k) = nutil(i,j-1,k)
          ! fill ghost cells
          do a=1,ngh_r
             nutil_n(i,j+a,k) = nutil_n(i,j,k)
          enddo
       enddo
    enddo

  end subroutine neu_jmax_rans

  !===============================================================================
  subroutine neu_kmin_rans
  !===============================================================================
    !> BC: boundary condition at kmin (top) 
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k,a
    ! ---------------------------------------------------------------------------
    ! Index of top boundary
    ! =====================
    k=1

    do j=1,ny
       do i=1,nx
          nutil_n(i,j,k) = nutil_n(i,j,k+1)
          nutil(i,j,k) = nutil(i,j,k+1)
          ! fill ghost cells
          do a=1,ngh_r
             nutil_n(i,j,k-a) = nutil_n(i,j,k)
          enddo
       enddo
    enddo

  end subroutine neu_kmin_rans

  !===============================================================================
  subroutine neu_kmax_rans
  !===============================================================================
    !> BC: boundary condition at kmax (top) 
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k,a
    ! ---------------------------------------------------------------------------
    ! Index of top boundary
    ! =====================
    k=nz

    do j=1,ny
       do i=1,nx
          nutil_n(i,j,k) = nutil_n(i,j,k-1)
          nutil(i,j,k) = nutil(i,j,k-1)
          ! fill ghost cells
          do a=1,ngh_r
             nutil_n(i,j,k+a) = nutil_n(i,j,k)
          enddo
       enddo
    enddo

  end subroutine neu_kmax_rans

end module mod_neumann_rans

