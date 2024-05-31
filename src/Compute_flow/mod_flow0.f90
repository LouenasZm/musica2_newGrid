!==============================================================================================
module mod_flow0
!==============================================================================================
  !> Module for time-averaged flow arrays [obsolete]
!==============================================================================================
  use precision
  implicit none
  ! -------------------------------------------------------------------------------------------
  ! averaged primitive variables [used in Tam & Dong's BC]
  real(wp), dimension(:,:,:), allocatable :: rho0,u0,v0,w0,p0,T0
  ! -------------------------------------------------------------------------------------------

contains

    !===============================================================================
    subroutine mean0
    !===============================================================================
      !> Compute online time-averaged primitive variables
    !===============================================================================
      use mod_flow
      use mod_time
      use mod_constant
      implicit none
      ! ---------------------------------------------------------------------------
      real(wp) :: inn,ntm1
      real(wp) :: eps
      ! ---------------------------------------------------------------------------
      ! Temporary : TO BE CHANGED
      integer :: n_moy
      integer :: i,j,k

      ! inn = 1.0_wp/dble(ntotal)
      ! ntm1=dble(ntotal-1)

      ! rho0= (ntm1*rho0 +rho_n) *inn
      ! u0  = (ntm1*  u0 + uu  ) *inn
      ! v0  = (ntm1*  v0 + vv  ) *inn
      ! w0  = (ntm1*  w0 + ww  ) *inn
      ! p0  = (ntm1*  p0 + prs ) *inn
      ! T0  = (ntm1*  T0 + Tmp ) *inn

      ! Test 2: TO BE CHANGED
      if (BC_face(1,1)%is_mean_ref) then
         if (BC_face(1,1)%sort==-1) then
            eps = 0.1_wp
            n_moy = 100
            inn = 1.0_wp/dble(n_moy)
            ntm1=dble(n_moy-1)
            do k=nz1,nz2
               do j=ny1,ny2
                  do i=1,8
                     rho0(i,j,k)= (ntm1*rho0(i,j,k) +rho_n(i,j,k))*eps*inn + (1.0_wp-eps)*BC_face(1,1)%Uref(j,k,1)
                     u0(i,j,k)  = (ntm1*  u0(i,j,k) +   uu(i,j,k))*eps*inn + (1.0_wp-eps)*BC_face(1,1)%Uref(j,k,2)
                     v0(i,j,k)  = (ntm1*  v0(i,j,k) +   vv(i,j,k))*eps*inn + (1.0_wp-eps)*BC_face(1,1)%Uref(j,k,3)
                     w0(i,j,k)  = (ntm1*  w0(i,j,k) +   ww(i,j,k))*eps*inn + (1.0_wp-eps)*BC_face(1,1)%Uref(j,k,4)
                     p0(i,j,k)  = (ntm1*  p0(i,j,k) +  prs(i,j,k))*eps*inn + (1.0_wp-eps)*BC_face(1,1)%Uref(j,k,5)
                     T0(i,j,k)  = (ntm1*  T0(i,j,k) +  Tmp(i,j,k))*eps*inn + (1.0_wp-eps)*Tin_ref(j)
                  enddo
               enddo
            enddo

            inn = 1.0_wp/dble(ntotal)
            ntm1=dble(ntotal-1)
            do k=nz1,nz2
               do j=ny1,ny2
                  do i=9,nx2
                     rho0(i,j,k) = (ntm1*rho0(i,j,k) + rho_n(i,j,k))*inn
                     u0(i,j,k)   = (ntm1*  u0(i,j,k) +    uu(i,j,k))*inn
                     v0(i,j,k)   = (ntm1*  v0(i,j,k) +    vv(i,j,k))*inn
                     w0(i,j,k)   = (ntm1*  w0(i,j,k) +    ww(i,j,k))*inn
                     p0(i,j,k)   = (ntm1*  p0(i,j,k) +   prs(i,j,k))*inn
                     T0(i,j,k)   = (ntm1*  T0(i,j,k) +   Tmp(i,j,k))*inn
                  enddo
               enddo
            enddo
         else
            inn = 1.0_wp/dble(ntotal)
            ntm1=dble(ntotal-1)
            rho0= (ntm1*rho0 +rho_n) *inn
            u0  = (ntm1*  u0 + uu  ) *inn
            v0  = (ntm1*  v0 + vv  ) *inn
            w0  = (ntm1*  w0 + ww  ) *inn
            p0  = (ntm1*  p0 + prs ) *inn
            T0  = (ntm1*  T0 + Tmp ) *inn
         endif
      else
         inn = 1.0_wp/dble(ntotal)
         ntm1=dble(ntotal-1)
         rho0= (ntm1*rho0 +rho_n) *inn
         u0  = (ntm1*  u0 + uu  ) *inn
         v0  = (ntm1*  v0 + vv  ) *inn
         w0  = (ntm1*  w0 + ww  ) *inn
         p0  = (ntm1*  p0 + prs ) *inn
         T0  = (ntm1*  T0 + Tmp ) *inn
      endif

    end subroutine mean0

end module mod_flow0
