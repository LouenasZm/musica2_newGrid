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

  !******************!
  ! Spalart Allmaras !
  !******************!

  !===============================================================================
  subroutine neu_imin_sa
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

  end subroutine neu_imin_sa

  !===============================================================================
  subroutine neu_imax_sa
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

  end subroutine neu_imax_sa

  !===============================================================================
  subroutine neu_jmin_sa
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

  end subroutine neu_jmin_sa

  !===============================================================================
  subroutine neu_jmax_sa
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

  end subroutine neu_jmax_sa

  !===============================================================================
  subroutine neu_kmin_sa
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

  end subroutine neu_kmin_sa

  !===============================================================================
  subroutine neu_kmax_sa
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

  end subroutine neu_kmax_sa

 !*******************!
   ! SA-Gamma-Re_theta !
   !*******************!

  !===============================================================================
  subroutine neu_imin_SA_transition
    !===============================================================================
      !> BC: boundary condition at imax (right)
    !===============================================================================
      implicit none
      ! ---------------------------------------------------------------------------
      integer :: i,j,k,a
      real(wp) :: C_neu_ksi,C_neu_eta,C_neu,C_rho
      ! ---------------------------------------------------------------------------
      ! Index of left boundary
      ! ======================
      i=1
  
      do k=1,nz
         do j=1,ny
            ! nu_til
            nutil_n(i,j,k) = nutil_n(i+1,j,k)
            nutil(i,j,k) = nutil(i+1,j,k)
            ! Intermittency factor
            rhogamma_n(i,j,k) = rhogamma_n(i+1,j,k)
            rho_gamma(i,j,k) = rho_gamma(i+1,j,k)
            ! Re_theta
            rhore_theta_n(i,j,k) = rhore_theta_n(i+1,j,k)
            rhore_theta(i,j,k) = rhore_theta(i+1,j,k)
            ! fill ghost cells
            do a=-ngh_r,-1
               nutil_n(i+a,j,k) = nutil_n(i,j,k)
               rhogamma_n(i+a,j,k) = rhogamma_n(i,j,k)
               rhore_theta_n(i+a,j,k) = rhore_theta_n(i,j,k)
            enddo
         enddo
      enddo
  
    end subroutine neu_imin_SA_transition
  
    !===============================================================================
    subroutine neu_imax_SA_transition
    !===============================================================================
      !> BC: boundary condition at imax (right)
    !===============================================================================
      implicit none
      ! ---------------------------------------------------------------------------
      integer :: i,j,k,a
      real(wp) :: C_neu_ksi,C_neu_eta,C_neu,C_rho
      ! ---------------------------------------------------------------------------
      ! Index of right boundary
      ! =======================
      i=nx
  
      do k=1,nz
         do j=1,ny
            ! Nu_tilde
            nutil_n(i,j,k) = nutil_n(i-1,j,k)
            nutil(i,j,k) = nutil(i-1,j,k)
            ! Intermittency factor:
            rhogamma_n(i,j,k) = rhogamma_n(i-1,j,k)
            rho_gamma(i,j,k) = rho_gamma(i-1,j,k)
            ! Re_theta
            rhore_theta_n(i,j,k) = rhore_theta_n(i-1,j,k)
            rhore_theta(i,j,k) = rhore_theta(i-1,j,k)
            ! fill ghost cells
            do a=1,ngh_r
               nutil_n(i+a,j,k) = nutil_n(i,j,k)
               rhogamma_n(i+a,j,k) = rhogamma_n(i,j,k)
               rhore_theta_n(i+a,j,k) = rhore_theta_n(i,j,k)
            enddo
         enddo
      enddo
  
    end subroutine neu_imax_SA_transition
  
    !===============================================================================
    subroutine neu_jmin_SA_transition
    !===============================================================================
      !> BC: boundary condition at jmin (bottom)
    !===============================================================================
      implicit none
      ! ---------------------------------------------------------------------------
      integer :: i,j,k,a
      real(wp) :: C_neu_ksi,C_neu_eta,C_neu,C_rho
      ! ---------------------------------------------------------------------------
      ! Index of bottom boundary
      ! ========================
      j=1
  
      do k=1,nz
         do i=1,nx         
            ! Nu_tilde
            nutil_n(i,j,k) = nutil_n(i,j+1,k)
            nutil(i,j,k) = nutil(i,j+1,k)
            ! Intermittency factor:
            rhogamma_n(i,j,k) = rhogamma_n(i,j+1,k)
            rho_gamma(i,j,k) = rho_gamma(i,j+1,k)
            ! Re_theta
            rhore_theta_n(i,j,k) = rhore_theta_n(i-1,j,k)
            rhore_theta(i,j,k) = rhore_theta(i-1,j,k)
            ! fill ghost cells
            do a=-1,-ngh_r
                nutil_n(i,j+a,k) = nutil_n(i,j,k)
                rhogamma_n(i,j+a,k) = rhogamma_n(i,j,k)
                rhore_theta_n(i,j+a,k) = rhore_theta_n(i,j,k)
            enddo
         enddo
      enddo
  
    end subroutine neu_jmin_SA_transition
  
    !===============================================================================
    subroutine neu_jmax_SA_transition
    !===============================================================================
      !> BC: boundary condition at jmax (top) 
    !===============================================================================
      implicit none
      ! ---------------------------------------------------------------------------
      integer :: i,j,k,a
      real(wp) :: C_neu_ksi,C_neu_eta,C_neu,C_rho
      ! ---------------------------------------------------------------------------
      ! Index of top boundary
      ! =====================
      j=ny
  
      do k=1,nz
         do i=1,nx
            ! Nu_tilde
            nutil_n(i,j,k) = nutil_n(i,j-1,k)
            nutil(i,j,k) = nutil(i,j-1,k)
            ! Intermittency factor:
            rhogamma_n(i,j,k) = rhogamma_n(i,j-1,k)
            rho_gamma(i,j,k) = rho_gamma(i,j-1,k)
            ! Re_theta:
            rhore_theta_n(i,j,k) = rhore_theta_n(i,j-1,k)
            rhore_theta(i,j,k) = rhore_theta(i,j-1,k)
            ! fill ghost cells
            do a=1,ngh_r
               nutil_n(i,j+a,k) = nutil_n(i,j,k)
               rhogamma_n(i,j+a,k) = rhogamma_n(i,j,k)
               rhore_theta_n(i,j+a,k) = rhore_theta_n(i,j,k)
            enddo
         enddo
      enddo
  
    end subroutine neu_jmax_SA_transition
  
    !===============================================================================
    subroutine neu_kmin_SA_transition
    !===============================================================================
      !> BC: boundary condition at kmin (top) 
    !===============================================================================
      implicit none
      ! ---------------------------------------------------------------------------
      integer :: i,j,k,a
      real(wp) :: C_neu_ksi,C_neu_eta,C_neu,C_rho
      ! ---------------------------------------------------------------------------
      ! Index of top boundary
      ! =====================
      k=1
  
      do j=1,ny
         do i=1,nx
            ! Nu_tilde
            nutil_n(i,j,k) = nutil_n(i,j,k+1)
            nutil(i,j,k) = nutil(i,j,k+1)
            ! Intermittency factor
            rhogamma_n(i,j,k) = rhogamma_n(i,j,k+1)
            rho_gamma(i,j,k) = rho_gamma(i,j,k+1)
            ! Re_theta
            rhore_theta_n(i,j,k) = rhore_theta_n(i,j,k+1)
            rhore_theta(i,j,k) = rhore_theta(i,j,k+1)
            ! fill ghost cells
            do a=1,ngh_r
               nutil_n(i,j,k-a) = nutil_n(i,j,k)
               rhogamma_n(i,j,k-a) = rhogamma_n(i,j,k)
               rhore_theta_n(i,j,k-a) = rhore_theta_n(i,j,k)
            enddo
         enddo
      enddo
  
    end subroutine neu_kmin_SA_transition
  
    !===============================================================================
    subroutine neu_kmax_SA_transition
    !===============================================================================
      !> BC: boundary condition at kmax (top) 
    !===============================================================================
      implicit none
      ! ---------------------------------------------------------------------------
      integer :: i,j,k,a
      real(wp) :: C_neu_ksi,C_neu_eta,C_neu,C_rho
      ! ---------------------------------------------------------------------------
      ! Index of top boundary
      ! =====================
      k=nz
  
      do j=1,ny
         do i=1,nx
            ! Nu_tilde
            nutil_n(i,j,k) = nutil_n(i,j,k-1)
            nutil(i,j,k) = nutil(i,j,k-1)
            ! Intermittency factor:
            rhogamma_n(i,j,k) = rhogamma_n(i,j,k-1)
            rho_gamma(i,j,k) = rho_gamma(i,j,k-1)
            ! Re_theta:
            rhore_theta_n(i,j,k) = rhore_theta_n(i,j,k-1)
            rhore_theta(i,j,k) = rhore_theta(i,j,k-1)
            ! fill ghost cells
            do a=1,ngh_r
               nutil_n(i,j,k+a) = nutil_n(i,j,k)
               rhogamma_n(i,j,k+a) = rhogamma_n(i,j,k)
               rhore_theta_n(i,j,k+a) = rhore_theta_n(i,j,k)
            enddo
         enddo
      enddo
  
    end subroutine neu_kmax_SA_transition

end module mod_neumann_rans

