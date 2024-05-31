!==============================================================================
module mod_artvisc_shock_rans
!==============================================================================
  !> Module for artificial viscosity with Jameson's type shock capturing
!==============================================================================
  use mod_flow
  use mod_artvisc
  use mod_rans
  implicit none
  ! ---------------------------------------------------------------------------
  ! declaration of coefficients in mod_artvisc.f90
  ! amplitude dissipation
  real(wp), private :: xnu0 ! low-order smoothing term
  real(wp), private :: xnu10,xnu8,xnu6,xnu4,xnu2
  ! ---------------------------------------------------------------------------

contains

  !============================================================================
  subroutine init_artvisc_shock_rans(dissip_coeff,dissip_shock)
  !============================================================================
    !> Initialization of artificial viscosity
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    real(wp), intent(inout) :: dissip_coeff,dissip_shock
    ! -------------------------------------------------------------------------


    !!!! /!\ Temporary ??????
    dissip_coeff=0.15_wp
!     dissip_shock=1.0_wp


    ! Artificial viscosity coefficients (DNC-Jameson)
    ! =================================
    ! /!\ Change sign because added to increment
    ! order 1
    d2(1)=-1.0_wp

    ! order 3
    d4(1)=-3.0_wp
    d4(2)= 1.0_wp

    ! Dissipation coefficients [from param.ini]
    ! ========================
    xnu4= dissip_coeff/12.0_wp
    xnu2= 0.0_wp*dissip_coeff

    ! low-order dissipation term [from param.ini]
    ! ==========================
    ! /!\ Change sign because added to increment
    xnu0=-dissip_shock
    
  end subroutine init_artvisc_shock_rans

  !==============================================================================
  subroutine artvisc_o3_shock_SA
  !==============================================================================
    !> Apply artificial viscosity - Conservative version DNC order 3
    !> Modified version CM 08/2022
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: vc1,vc2,vc,srm12,srp12,sr
    real(wp), dimension(0:nx,0:ny,0:nz) :: Dnutil
    ! ---------------------------------------------------------------------------
    real(wp) :: xk0,xk2,xk4
    ! real(wp), dimension(ndx_d1:nfx_d1,ndy_d1:nfy_d1,ndz_d1:nfz_d1) :: sens
    real(wp), dimension(ndx_d1:nfx_d1,ndy_d1:nfy_d1,ndz_d1:nfz_d1) :: psens
    ! ---------------------------------------------------------------------------
    real(wp) :: psens_trad,psens_tvd,csens

    csens=1.0_wp
    ! Compute pressure sensor along i-direction
    ! =========================================
    do k=1,nz
       do j=1,ny
          do i=ndx_d1_r,nfx_d1_r
!              psens(i,j,k)=sens(i,j,k)*abs(prs(i-1,j,k)-2.0_wp*prs(i,j,k)+prs(i+1,j,k)) / &
!                                       abs(prs(i-1,j,k)+2.0_wp*prs(i,j,k)+prs(i+1,j,k))

             ! an alternative sensor (implemented by OY)
             ! a mix of a TVD concept and the traditional one (suggested by Swanson et al. [1998])
             psens_trad = prs(i-1,j,k)+2.0_wp*prs(i,j,k)+prs(i+1,j,k)
             psens_tvd = abs(prs(i+1,j,k)-prs(i,j,k))+abs(prs(i,j,k)-prs(i-1,j,k))
             psens(i,j,k)=abs(prs(i-1,j,k)-2.0_wp*prs(i,j,k)+prs(i+1,j,k)) / &
                                      (csens*psens_trad + (1.0_wp-csens)*psens_tvd)

             !!uvar(i,j,k,2)=psens(i,j,k)
          enddo
       enddo
    enddo
    
    ! Artificial viscosity in i-direction [compute Drho at i+1/2]
    ! ===================================
    
    if (is_boundary(1,1)) then
       ! no dissipation at i=1 -> initialize increment array
       do j=1,ny
           Knutil(1,j,:)=0.0_wp
       enddo
       
       ! 1st-order centered dissipation
       ! ------------------------------
       i=1
       do k=1,nz
          do j=1,ny             
             xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))
             xk2=max(0.0_wp,xnu2+xk0)
             !uvar(i,j,k,1)=xk2
             
             Dnutil(i,j,k) = (nutil_n(i+1,j,k)-nutil_n(i,j,k))*xk0 &
                         + d2(1)*(nutil_n(i+1,j,k)- nutil_n(i,j,k))*xk2       
          enddo
       enddo
    endif
       
    ! 3rd-order centered dissipation
    ! ------------------------------
    do k=1,nz
       do j=1,ny 
          do i=ndx_di_r,nfx_di_r
             xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))
             xk4=max(0.0_wp,xnu4+xk0)
             !uvar(i,j,k,1)=xk4/xnu4
             
             Dnutil(i,j,k) = (nutil_n(i+1,j,k)-nutil_n(i,j,k))*xk0      &
                         + (d4(1)*(nutil_n(i+1,j,k)-nutil_n(i  ,j,k)) &
                           +d4(2)*(nutil_n(i+2,j,k)-nutil_n(i-1,j,k)))*xk4
          enddo
       enddo
    enddo

    if (is_boundary(1,2)) then           
       ! 1st-order centered dissipation
       ! ------------------------------
       i=nx-1
       do k=1,nz
          do j=1,ny             
             xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))
             xk2=max(0.0_wp,xnu2+xk0)
             !uvar(i,j,k,1)=xk2
             
             Dnutil(i,j,k) = (nutil_n(i+1,j,k)-nutil_n(i,j,k))*xk0 &
                         + d2(1)*(nutil_n(i+1,j,k)-nutil_n(i,j,k))*xk2
          enddo
       enddo
       
       ! no dissipation at i=nx -> initialize increment array
       i=nx
       do j=1,ny
           Knutil(nx,j,:)=0.0_wp
       enddo
    endif
    
    ! Compute artificial dissipation terms
    ! ------------------------------------
    if (is_curv) then
       do k=1,nz
          do j=1,ny
             do i=ndx_d_r,nfx_d_r
                ! contravariant velocity at i+1/2 (vc2) and i-1/2 (vc1)
                vc=(uu(i,j,k)*y_eta_v(i,j)-vv(i,j,k)*x_eta_v(i,j))*ijacob_v(i,j)
                vc2=(uu(i+1,j,k)*y_eta_v(i+1,j)-vv(i+1,j,k)*x_eta_v(i+1,j))*ijacob_v(i+1,j)
                vc1=(uu(i-1,j,k)*y_eta_v(i-1,j)-vv(i-1,j,k)*x_eta_v(i-1,j))*ijacob_v(i-1,j)
                ! compute spectral radius at i+1/2 (srp12) and i-1/2 (srm12)
                sr=abs(vc)+c_(i,j,k)*g_ksi(i,j)
                srp12=0.5_wp*(abs(vc2)+c_(i+1,j,k)*g_ksi(i+1,j)+sr)
                srm12=0.5_wp*(abs(vc1)+c_(i-1,j,k)*g_ksi(i-1,j)+sr)
                ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dksi
                Knutil(i,j,k) = Knutil(i,j,k) -(srm12*Dnutil(i-1,j,k) -srp12*Dnutil(i,j,k))
             enddo
          enddo
       enddo
    elseif (is_curv3) then
       do k=1,nz
          do j=1,ny
             do i=ndx_d_r,nfx_d_r
                ! contravariant velocity at i+1/2 (vc2) and i-1/2 (vc1)
                vc=(uu(i,j,k)*ksi_x_v(i,j,k)+vv(i,j,k)*ksi_y_v(i,j,k)+ww(i,j,k)*ksi_z_v(i,j,k))*ijacob3_v(i,j,k)
                vc2=(uu(i+1,j,k)*ksi_x_v(i+1,j,k)+vv(i+1,j,k)*ksi_y_v(i+1,j,k)+ww(i+1,j,k)*ksi_z_v(i+1,j,k))*ijacob3_v(i+1,j,k)
                vc1=(uu(i-1,j,k)*ksi_x_v(i-1,j,k)+vv(i-1,j,k)*ksi_y_v(i-1,j,k)+ww(i-1,j,k)*ksi_z_v(i-1,j,k))*ijacob3_v(i-1,j,k)
                ! compute spectral radius at i+1/2 (srp12) and i-1/2 (srm12)
                sr=abs(vc)+c_(i,j,k)*g3_ksi(i,j,k)
                srp12=0.5_wp*(abs(vc2)+c_(i+1,j,k)*g3_ksi(i+1,j,k)+sr)
                srm12=0.5_wp*(abs(vc1)+c_(i-1,j,k)*g3_ksi(i-1,j,k)+sr)
                ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dksi
                Knutil(i,j,k) = Knutil(i,j,k) -(srm12*Dnutil(i-1,j,k) -srp12*Dnutil(i,j,k))
             enddo
          enddo
       enddo
    else
       do k=1,nz
          do j=1,ny
             do i=ndx_d_r,nfx_d_r
                ! compute spectral radius at i+1/2 (srp12) and i-1/2 (srm12)
                sr=abs(uu(i,j,k))+c_(i,j,k)
                srp12=0.5_wp*(abs(uu(i+1,j,k))+c_(i+1,j,k)+sr)
                srm12=0.5_wp*(abs(uu(i-1,j,k))+c_(i-1,j,k)+sr)
                ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dx
                Knutil(i,j,k) = -(srm12*Dnutil(i-1,j,k) -srp12*Dnutil(i,j,k))*idx(i)
                ! store dissipation coefficient
                !!uvar(i,j,k,1)=sr
             enddo
          enddo
       enddo
    endif

    ! Compute pressure sensor along j-direction
    ! =========================================
    do k=1,nz
       do j=ndy_d1_r,nfy_d1_r
          do i=1,nx
!              psens(i,j,k)=sens(i,j,k)*abs(prs(i,j-1,k)-2.0_wp*prs(i,j,k)+prs(i,j+1,k)) / &
!                                       abs(prs(i,j-1,k)+2.0_wp*prs(i,j,k)+prs(i,j+1,k))

             ! an alternative sensor (implemented by OY)
             ! a mix of a TVD concept and the traditional one (suggested by Swanson et al. [1998])
             psens_trad = prs(i,j-1,k)+2.0_wp*prs(i,j,k)+prs(i,j+1,k)
             psens_tvd = abs(prs(i,j+1,k)-prs(i,j,k))+abs(prs(i,j,k)-prs(i,j-1,k))
             psens(i,j,k)=abs(prs(i,j-1,k)-2.0_wp*prs(i,j,k)+prs(i,j+1,k)) / &
                                      (csens*psens_trad + (1.0_wp-csens)*psens_tvd)
          enddo
       enddo
    enddo
 
    ! Artificial viscosity in y-direction [compute Drho at j+1/2]
    ! ===================================

    if (is_boundary(2,1)) then
       ! 1st-order centered dissipation
       ! ------------------------------
       j=1
       do k=1,nz
          do i=1,nx
             xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))
             xk2=max(0.0_wp,xnu2+xk0)
             !uvar(i,j,k,2)=xk2
             
             Dnutil(i,j,k) = (nutil_n(i,j+1,k)-nutil_n(i,j,k))*xk0 &
                         + d2(1)*(nutil_n(i,j+1,k)-nutil_n(i,j,k))*xk2
          enddo
       enddo
    endif
       
    ! 3rd-order centered dissipation
    ! ------------------------------
    do k=1,nz
       do j=ndy_di_r,nfy_di_r
          do i=1,nx
             xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))
             xk4=max(0.0_wp,xnu4+xk0)
             !uvar(i,j,k,2)=xk4/xnu4
             
             Dnutil(i,j,k) = (nutil_n(i,j+1,k)-nutil_n(i,j,k))*xk0      &
                         + (d4(1)*(nutil_n(i,j+1,k)-nutil_n(i,j  ,k)) &
                           +d4(2)*(nutil_n(i,j+2,k)-nutil_n(i,j-1,k)))*xk4
          enddo
       enddo
    enddo

    if (is_boundary(2,2)) then
       ! 1st-order centered dissipation
       ! ------------------------------
       j=ny-1
       do k=1,nz
          do i=1,nx
             xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))
             xk2=max(0.0_wp,xnu2+xk0)
             !uvar(i,j,k,2)=xk2
             
             Dnutil(i,j,k) = (nutil_n(i,j+1,k)-nutil_n(i,j,k))*xk0 &
                         + d2(1)*(nutil_n(i,j+1,k)-nutil_n(i,j,k))*xk2
          enddo
       enddo
    endif

    ! Compute artificial dissipation terms
    ! ------------------------------------
    if (is_curv) then
       do k=1,nz
          do j=ndy_d_r,nfy_d_r
             do i=1,nx
                ! contravariant velocity at j+1/2 (vc2) and j-1/2 (vc1)
                vc=(vv(i,j,k)*x_ksi_v(i,j)-uu(i,j,k)*y_ksi_v(i,j))*ijacob_v(i,j)
                vc2=(vv(i,j+1,k)*x_ksi_v(i,j+1)-uu(i,j+1,k)*y_ksi_v(i,j+1))*ijacob_v(i,j+1)
                vc1=(vv(i,j-1,k)*x_ksi_v(i,j-1)-uu(i,j-1,k)*y_ksi_v(i,j-1))*ijacob_v(i,j-1)
                ! compute spectral radius at j+1/2 (srp12) and j-1/2 (srm12)
                sr=abs(vc)+c_(i,j,k)*g_eta(i,j)
                srp12=0.5_wp*(abs(vc2)+c_(i,j+1,k)*g_eta(i,j+1)+sr)
                srm12=0.5_wp*(abs(vc1)+c_(i,j-1,k)*g_eta(i,j-1)+sr)
                ! compute dissipation term: delta_j/delta_y=(D_{j+1/2}-D_{j-1/2})/dy
                Knutil(i,j,k) = Knutil(i,j,k) -(srm12*Dnutil(i,j-1,k) -srp12*Dnutil(i,j,k))
             enddo
          enddo
       enddo   
    elseif (is_curv3) then
       do k=1,nz
          do j=1,ny
             do i=ndx_d_r,nfx_d_r
                ! contravariant velocity at i+1/2 (vc2) and i-1/2 (vc1)
                vc=(uu(i,j,k)*eta_x_v(i,j,k)+vv(i,j,k)*eta_y_v(i,j,k)+ww(i,j,k)*eta_z_v(i,j,k))*ijacob3_v(i,j,k)
                vc2=(uu(i,j+1,k)*eta_x_v(i,j+1,k)+vv(i,j+1,k)*eta_y_v(i,j+1,k)+ww(i,j+1,k)*eta_z_v(i,j+1,k))*ijacob3_v(i,j+1,k)
                vc1=(uu(i,j-1,k)*eta_x_v(i,j-1,k)+vv(i,j-1,k)*eta_y_v(i,j-1,k)+ww(i,j-1,k)*eta_z_v(i,j-1,k))*ijacob3_v(i,j-1,k)
                ! compute spectral radius at i+1/2 (srp12) and i-1/2 (srm12)
                sr=abs(vc)+c_(i,j,k)*g3_eta(i,j,k)
                srp12=0.5_wp*(abs(vc2)+c_(i,j+1,k)*g3_eta(i,j+1,k)+sr)
                srm12=0.5_wp*(abs(vc1)+c_(i,j-1,k)*g3_eta(i,j-1,k)+sr)
                ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dksi
                Knutil(i,j,k) = Knutil(i,j,k) -(srm12*Dnutil(i,j-1,k) -srp12*Dnutil(i,j,k))
             enddo
          enddo
       enddo 
    else
       do k=1,nz
          do j=ndy_d_r,nfy_d_r
             do i=1,nx
                ! compute spectral radius at j+1/2 (srp12) and j-1/2 (srm12)
                sr=abs(vv(i,j,k))+c_(i,j,k)
                srp12=0.5_wp*(abs(vv(i,j+1,k))+c_(i,j+1,k)+sr)
                srm12=0.5_wp*(abs(vv(i,j-1,k))+c_(i,j-1,k)+sr)
                ! compute dissipation term: delta_j/delta_y=(D_{j+1/2}-D_{j-1/2})/dy
                Knutil(i,j,k) = Knutil(i,j,k) -(srm12*Dnutil(i,j-1,k) -srp12*Dnutil(i,j,k) )*idy(j)
             enddo
          enddo
       enddo
    endif
 
    !*******************
    if (.not.is_2D) then
    !*******************

       ! Compute pressure sensor along k-direction
       ! =========================================
       do i=1,nx
          do j=1,ny
             do k=ndz_d1_r,nfz_d1_r
!                 psens(i,j,k)=sens(i,j,k)*abs(prs(i,j,k-1)-2.0_wp*prs(i,j,k)+prs(i,j,k+1)) / &
!                                          abs(prs(i,j,k-1)+2.0_wp*prs(i,j,k)+prs(i,j,k+1))

                ! an alternative sensor (implemented by OY)
                ! a mix of a TVD concept and the traditional one (suggested by Swanson et al. [1998])
                psens_trad = prs(i,j,k-1)+2.0_wp*prs(i,j,k)+prs(i,j,k+1)
                psens_tvd = abs(prs(i,j,k+1)-prs(i,j,k))+abs(prs(i,j,k)-prs(i,j,k-1))
                psens(i,j,k)=abs(prs(i,j,k-1)-2.0_wp*prs(i,j,k)+prs(i,j,k+1)) / &
                                      (csens*psens_trad + (1.0_wp-csens)*psens_tvd)
             enddo
          enddo
       enddo
 
       ! Artificial viscosity in k-direction [compute Drho at k+1/2]
       ! ===================================
       if (is_boundary(3,1)) then
          ! 1st-order centered dissipation
          ! ------------------------------
          k=1
          do j=1,ny
             do i=1,nx
                xk0=xnu0*max(psens(i,j,k+1),psens(i,j,k))
                xk2=max(0.0_wp,xnu2+xk0)
             
                Dnutil(i,j,k) = (nutil_n(i,j,k+1)-nutil_n(i,j,k))*xk0 &
                            + d2(1)*(nutil_n(i,j,k+1)-nutil_n(i,j,k))*xk2
             enddo
          enddo
       endif

       ! 3rd-order centered dissipation
       ! ------------------------------
       do k=ndz_di_r,nfz_di_r
          do j=1,ny
             do i=1,nx
                xk0=xnu0*max(psens(i,j,k+1),psens(i,j,k))
                xk4=max(0.0_wp,xnu4+xk0)
             
                Dnutil(i,j,k) = (nutil_n(i,j,k+1)-nutil_n(i,j,k))*xk0      &
                            + (d4(1)*(nutil_n(i,j,k+1)-nutil_n(i,j,k  )) &
                              +d4(2)*(nutil_n(i,j,k+2)-nutil_n(i,j,k-1)))*xk4
             enddo
          enddo
       enddo
       
       if (is_boundary(3,2)) then
          ! 1st-order centered dissipation
          ! ------------------------------
          k=nz-1
          do j=1,ny
             do i=1,nx
                xk0=xnu0*max(psens(i,j,k+1),psens(i,j,k))
                xk2=max(0.0_wp,xnu2+xk0)
             
                Dnutil(i,j,k) = (nutil_n(i,j,k+1)-nutil_n(i,j,k))*xk0 &
                            + d2(1)*(nutil_n(i,j,k+1)-nutil_n(i,j,k))*xk2
             enddo
          enddo
       endif

       ! Compute artificial dissipation terms
       ! ------------------------------------
       if (is_curv3) then
          do k=1,nz
             do j=1,ny
                do i=ndx_d_r,nfx_d_r
                   ! contravariant velocity at i+1/2 (vc2) and i-1/2 (vc1)
                   vc=(uu(i,j,k)*phi_x_v(i,j,k)+vv(i,j,k)*phi_y_v(i,j,k)+ww(i,j,k)*phi_z_v(i,j,k))*ijacob3_v(i,j,k)
                   vc2=(uu(i,j,k+1)*phi_x_v(i,j,k+1)+vv(i,j,k+1)*phi_y_v(i,j,k+1)+ww(i,j,k+1)*phi_z_v(i,j,k+1))*ijacob3_v(i,j,k+1)
                   vc1=(uu(i,j,k-1)*phi_x_v(i,j,k-1)+vv(i,j,k-1)*phi_y_v(i,j,k-1)+ww(i,j,k-1)*phi_z_v(i,j,k-1))*ijacob3_v(i,j,k-1)
                   ! compute spectral radius at i+1/2 (srp12) and i-1/2 (srm12)
                   sr=abs(vc)+c_(i,j,k)*g3_phi(i,j,k)
                   srp12=0.5_wp*(abs(vc2)+c_(i,j,k+1)*g3_phi(i,j,k+1)+sr)
                   srm12=0.5_wp*(abs(vc1)+c_(i,j,k-1)*g3_phi(i,j,k+1)+sr)
                   ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dksi
                Knutil(i,j,k) = Knutil(i,j,k) -(srm12*Dnutil(i,j,k-1) -srp12*Dnutil(i,j,k))
                enddo
             enddo
          enddo 
       else
          do k=ndz_d_r,nfz_d_r
             do j=1,ny
                do i=1,nx
                   ! compute spectral radius at k+1/2 (srp12) and k-1/2 (srm12)
                   sr=abs(ww(i,j,k))+c_(i,j,k)
                   srp12=0.5_wp*(abs(ww(i,j,k+1))+c_(i,j,k+1)+sr)
                   srm12=0.5_wp*(abs(ww(i,j,k-1))+c_(i,j,k-1)+sr)
                   ! compute dissipation term: delta_k/delta_z=(D_{k+1/2}-D_{k-1/2})/dz
                   Knutil(i,j,k) = Knutil(i,j,k) -(srm12*Dnutil(i,j,k-1) -srp12*Dnutil(i,j,k) )*idz(k)
                enddo
             enddo
          enddo
       endif

    endif
    
  end subroutine artvisc_o3_shock_SA

end module mod_artvisc_shock_rans
