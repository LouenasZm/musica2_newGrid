!============================================================================
module mod_inlet_rans
!============================================================================
  !> author: CM
  !> date: March 2022
  !> Inlet boundary condition
!=============================================================================
  use mod_flow
  use mod_rans
  use mod_constant
  use mod_flow
  use mod_mpi
  implicit none
  !---------------------------------------------------------------------------

contains

  !******************!
  ! Spalart Allmaras !
  !******************!

  !===============================================================================
  subroutine inlet_imin_sa
  !===============================================================================
    !> BC: boundary condition at imax (right)
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    ! ---------------------------------------------------------------------------
    ! Index of left boundary
    ! ======================
    i=1

    do k=1,nz
       do j=1,ny
          nutil_n(i,j,k)=Nutref
          Knutil(i,j,k) = 0.0_wp
       enddo
    enddo

  end subroutine inlet_imin_sa

  !===============================================================================
  subroutine inlet_imax_sa
  !===============================================================================
    !> BC: boundary condition at imax (right)
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    ! ---------------------------------------------------------------------------
    ! Index of right boundary
    ! =======================
    i=nx

    do k=1,nz
       do j=1,ny
          nutil_n(i,j,k)=Nutref
          Knutil(i,j,k) = 0.0_wp
       enddo
    enddo

    if (STBL) then
       do k=1,nz
          do j=1,ny
             nutil_n(i,j,k)=nutil_n(i-1,j,k)
             Knutil(i,j,k) = 0.0_wp
          enddo
       enddo
    endif

  end subroutine inlet_imax_sa

  !===============================================================================
  subroutine inlet_jmin_sa
  !===============================================================================
    !> BC: boundary condition at jmin (bottom)
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    ! ---------------------------------------------------------------------------
    ! Index of bottom boundary
    ! ========================
    j=1

    do k=1,nz
       do i=1,nx
          nutil_n(i,j,k)=Nutref
          Knutil(i,j,k) = 0.0_wp
       enddo
    enddo

  end subroutine inlet_jmin_sa

  !===============================================================================
  subroutine inlet_jmax_sa
  !===============================================================================
    !> BC: boundary condition at jmax (top) 
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    ! ---------------------------------------------------------------------------
    ! Index of top boundary
    ! =====================
    j=ny

    do k=1,nz
       do i=1,nx
          nutil_n(i,j,k)=Nutref
          Knutil(i,j,k) = 0.0_wp
       enddo
    enddo

    if (STBL) then
       do k=1,nz
          do i=1,nx
             nutil_n(i,j,k)=nutil_n(i,j-1,k)
             Knutil(i,j,k) = 0.0_wp
          enddo
       enddo
    endif

  end subroutine inlet_jmax_sa

  !===============================================================================
  subroutine inlet_kmin_sa
  !===============================================================================
    !> BC: boundary condition at kmin (bottom)
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    ! ---------------------------------------------------------------------------
    ! Index of bottom boundary
    ! ========================
    k=1

    do j=1,ny
       do i=1,nx
          nutil_n(i,j,k)=Nutref
          Knutil(i,j,k) = 0.0_wp
       enddo
    enddo

  end subroutine inlet_kmin_sa

  !===============================================================================
  subroutine inlet_kmax_sa
  !===============================================================================
    !> BC: boundary condition at kmax (bottom)
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    ! ---------------------------------------------------------------------------
    ! Index of bottom boundary
    ! ========================
    k=nz

    do j=1,ny
       do i=1,nx
          nutil_n(i,j,k)=Nutref
          Knutil(i,j,k) = 0.0_wp
       enddo
    enddo

  end subroutine inlet_kmax_sa

    !****************************!
    !       SA-Gamma-Re_theta    !
    !****************************!

  !===============================================================================
  subroutine inlet_imin_SA_transition
    !===============================================================================
      !> BC: boundary condition at imax (right)
    !===============================================================================
      implicit none
      ! ---------------------------------------------------------------------------
      integer :: i,j,k,a,nbl
      ! ---------------------------------------------------------------------------
      ! Index of left boundary
      ! ======================
      i=1
  
      do k=1,nz
         do j=1,ny
            nutil_n(i,j,k)      = Nutref
            Knutil(i,j,k)       = 0.0_wp
            rhogamma_n(i,j,k)  = rho_n(i,j,k)*1.0_wp
            ! Re_theta
            if(tu_inlet .lt. 1.3_wp)then
                reynolds_theta(i,j,k)  = 1173.51_wp - 589.428*tu_inlet + 0.2196_wp/tu_inlet**2
             else
                reynolds_theta(i,j,k)  = 331.5_wp * (tu_inlet - 0.5668_wp)**(-0.671_wp)
             endif
         enddo
      enddo
  
      
      ! Inlet of 
    end subroutine inlet_imin_SA_transition
  
    !===============================================================================
    subroutine inlet_imax_SA_transition
    !===============================================================================
      !> BC: boundary condition at imax (right)
    !===============================================================================
      implicit none
      ! ---------------------------------------------------------------------------
      integer :: i,j,k,a,nbl
      ! ---------------------------------------------------------------------------
      ! Index of right boundary
      ! =======================
      i=nx
  
      do k=1,nz
         do j=1,ny
            nutil_n(i,j,k)=Nutref
            Knutil(i,j,k) = 0.0_wp
            rhogamma_n(i,j,k)  = rho_n(i,j,k)*1.0_wp
            ! Re_theta
            if(tu_inlet .lt. 1.3_wp)then
                reynolds_theta(i,j,k)  = 1173.51_wp - 589.428*tu_inlet + 0.2196_wp/tu_inlet**2
             else
                reynolds_theta(i,j,k)  = 331.5_wp * (tu_inlet - 0.5668_wp)**(-0.671_wp)
             endif
         enddo
      enddo
  
      if (STBL) then
         do k=1,nz
            do j=1,ny
               nutil_n(i,j,k)=nutil_n(i-1,j,k)
               Knutil(i,j,k) = 0.0_wp
            enddo
         enddo
      endif
  
    end subroutine inlet_imax_SA_transition
  
    !===============================================================================
    subroutine inlet_jmin_SA_transition
    !===============================================================================
      !> BC: boundary condition at jmin (bottom)
    !===============================================================================
      implicit none
      ! ---------------------------------------------------------------------------
      integer :: i,j,k
      ! ---------------------------------------------------------------------------
      ! Index of bottom boundary
      ! ========================
      j=1
  
      do k=1,nz
         do i=1,nx
            nutil_n(i,j,k)=Nutref
            Knutil(i,j,k) = 0.0_wp
         enddo
      enddo
  
    end subroutine inlet_jmin_SA_transition
  
    !===============================================================================
    subroutine inlet_jmax_SA_transition
    !===============================================================================
      !> BC: boundary condition at jmax (top) 
    !===============================================================================
      implicit none
      ! ---------------------------------------------------------------------------
      integer :: i,j,k,a,nbl
      ! ---------------------------------------------------------------------------
      ! Index of top boundary
      ! =====================
      j=ny
  
      do k=1,nz
         do i=1,nx
            nutil_n(i,j,k)=Nutref
            Knutil(i,j,k) = 0.0_wp
            rhogamma_n(i,j,k)  = rho_n(i,j,k)*1.0_wp
            ! Re_theta
            if(tu_inlet .lt. 1.3_wp)then
                reynolds_theta(i,j,k)  = 1173.51_wp - 589.428*tu_inlet + 0.2196_wp/tu_inlet**2
             else
                reynolds_theta(i,j,k)  = 331.5_wp * (tu_inlet - 0.5668_wp)**(-0.671_wp)
             endif
         enddo
      enddo
  
      if (STBL) then
         do k=1,nz
            do i=1,nx
               nutil_n(i,j,k)=nutil_n(i,j-1,k)
               Knutil(i,j,k) = 0.0_wp
            enddo
         enddo
      endif
  
    end subroutine inlet_jmax_SA_transition
  
    !===============================================================================
    subroutine inlet_kmin_SA_transition
    !===============================================================================
      !> BC: boundary condition at kmin (bottom)
    !===============================================================================
      implicit none
      ! ---------------------------------------------------------------------------
      integer :: i,j,k
      ! ---------------------------------------------------------------------------
      ! Index of bottom boundary
      ! ========================
      k=1
  
      do j=1,ny
         do i=1,nx
            nutil_n(i,j,k)=Nutref
            Knutil(i,j,k) = 0.0_wp
            rhogamma_n(i,j,k)  = rho_n(i,j,k)*1.0_wp
            ! Re_theta
            if(tu_inlet .lt. 1.3_wp)then
                reynolds_theta(i,j,k)  = 1173.51_wp - 589.428*tu_inlet + 0.2196_wp/tu_inlet**2
             else
                reynolds_theta(i,j,k)  = 331.5_wp * (tu_inlet - 0.5668_wp)**(-0.671_wp)
             endif
         enddo
      enddo
  
    end subroutine inlet_kmin_SA_transition
  
    !===============================================================================
    subroutine inlet_kmax_SA_transition
    !===============================================================================
      !> BC: boundary condition at kmax (bottom)
    !===============================================================================
      implicit none
      ! ---------------------------------------------------------------------------
      integer :: i,j,k
      ! ---------------------------------------------------------------------------
      ! Index of bottom boundary
      ! ========================
      k=nz
  
      do j=1,ny
         do i=1,nx
            nutil_n(i,j,k)=Nutref
            Knutil(i,j,k) = 0.0_wp
            rhogamma_n(i,j,k)  = rho_n(i,j,k)*1.0_wp
            ! Re_theta
            if(tu_inlet .lt. 1.3_wp)then
                reynolds_theta(i,j,k)  = 1173.51_wp - 589.428*tu_inlet + 0.2196_wp/tu_inlet**2
             else
                reynolds_theta(i,j,k)  = 331.5_wp * (tu_inlet - 0.5668_wp)**(-0.671_wp)
             endif
         enddo
      enddo
  
    end subroutine inlet_kmax_SA_transition      



end module mod_inlet_rans

