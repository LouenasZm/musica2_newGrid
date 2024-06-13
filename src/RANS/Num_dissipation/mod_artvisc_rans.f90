!==============================================================================
module mod_artvisc_rans
!==============================================================================
  !> Module for artificial viscosity
!==============================================================================
  use mod_flow
  use mod_artvisc
  implicit none

contains
  
  !==============================================================================
  subroutine init_artvisc_rans(dissip_coeff)
  !==============================================================================
    !> Initialization of artificial viscosity
  !==============================================================================
    use mod_time
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), intent(inout) :: dissip_coeff
    real(wp) :: xnu4,xnu2
    ! ---------------------------------------------------------------------------

    !!!! /!\ Temporary ??????
    dissip_coeff=0.6_wp

    ! Artificial viscosity coefficients (DNC-Jameson)
    ! =================================
    ! order 1
    d2(1)= 1.0_wp

    ! order 3
    d4(1)= 3.0_wp
    d4(2)=-1.0_wp

    ! dissipation coefficient [from param.ini]
    xnu4= dissip_coeff/12.0_wp
    xnu2= 0.*dissip_coeff/2.0_wp

    ! Multiply coeff by dissipation coefficient
    ! =========================================
     d4= d4*xnu4
     d2= d2*xnu2

  end subroutine init_artvisc_rans


  !**************************
  !    Spallart-Allmaras    *
  !**************************

  !===============================================================================
  subroutine artvisc_o3_SA
  !===============================================================================
    !> Apply artificial viscosity - Conservative version DNC order 3
    !> Modified version XG 06/2021
  !===============================================================================
    use mod_time ! for: deltat
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: vc1,vc2,vc,srm12,srp12,sr
    real(wp), dimension(0:nx,0:ny,0:nz) :: Dnutil
    ! ---------------------------------------------------------------------------

    ! Artificial viscosity in x-direction [compute Dnutil at i+1/2]
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
             Dnutil(i,j,k) = d2(1)*(nutil_n(i+1,j,k)-nutil_n(i,j,k))            
          enddo
       enddo
    endif

    ! Interior points: 3rd-order centered dissipation
    ! -----------------------------------------------
    do k=1,nz
       do j=1,ny
          do i=ndx_di_r,nfx_di_r
             Dnutil(i,j,k) = d4(1)*(nutil_n(i+1,j,k)-nutil_n(i  ,j,k)) &
                              + d4(2)*(nutil_n(i+2,j,k)-nutil_n(i-1,j,k))
          enddo
       enddo
    enddo

    if (is_boundary(1,2)) then
       ! 1st-order centered dissipation
       ! ------------------------------
       i=nx-1
       do k=1,nz
          do j=1,ny
             Dnutil(i,j,k) = d2(1)*(nutil_n(i+1,j,k)-nutil_n(i,j,k))            
          enddo
       enddo

       ! no dissipation at i=nx -> initialize increment array
       do j=1,ny
          Knutil(nx,j,:)=0.0_wp
       enddo
    endif
    
    ! Add dissipation to numerical fluxes
    ! -----------------------------------
    if (is_curv) then
       do k=1,nz
          do j=1,ny
             do i=ndx_d_r,nfx_d_r
                ! contravariant velocity at i+1/2 (vc2) and i-1/2 (vc1)
                vc=(uu(i,j,k)*y_eta(i,j)-vv(i,j,k)*x_eta(i,j))*ijacob(i,j)
                vc2=(uu(i+1,j,k)*y_eta(i+1,j)-vv(i+1,j,k)*x_eta(i+1,j))*ijacob(i+1,j)
                vc1=(uu(i-1,j,k)*y_eta(i-1,j)-vv(i-1,j,k)*x_eta(i-1,j))*ijacob(i-1,j)
                ! compute spectral radius at i+1/2 (srp12) and i-1/2 (srm12)
                sr=abs(vc)+c_(i,j,k)*g_ksi(i,j)
                srp12=0.5_wp*(abs(vc2)+c_(i+1,j,k)*g_ksi(i+1,j)+sr)
                srm12=0.5_wp*(abs(vc1)+c_(i-1,j,k)*g_ksi(i-1,j)+sr)
                ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dksi
                Knutil(i,j,k) = Knutil(i,j,k) -(srp12*Dnutil(i,j,k) -srm12*Dnutil(i-1,j,k))
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
                Knutil(i,j,k) = Knutil(i,j,k) -(srp12*Dnutil(i,j,k) -srm12*Dnutil(i-1,j,k))
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
                Knutil(i,j,k) = Knutil(i,j,k) -(srp12*Dnutil(i,j,k) -srm12*Dnutil(i-1,j,k) )*idx(i)
             enddo
          enddo
       enddo
    endif

    ! Artificial viscosity in y-direction [compute Dnutil at j+1/2]
    ! ===================================
    if (is_boundary(2,1)) then
       ! 1st-order centered dissipation
       ! ------------------------------
       j=1
       do k=1,nz
          do i=1,nx
             Dnutil(i,j,k) = d2(1)*(nutil_n(i,j+1,k)-nutil_n(i,j,k))
          enddo
       enddo
    endif
    
    ! Interior points: 3rd-order centered dissipation
    ! -----------------------------------------------
    do k=1,nz
       do j=ndy_di_r,nfy_di_r
          do i=1,nx
             Dnutil(i,j,k) = d4(1)*(nutil_n(i,j+1,k)-nutil_n(i,j  ,k)) &
                         + d4(2)*(nutil_n(i,j+2,k)-nutil_n(i,j-1,k))
          enddo
       enddo
    enddo

    if (is_boundary(2,2)) then
       ! 1st-order centered dissipation
       ! ------------------------------
       j=ny-1
       do k=1,nz
          do i=1,nx
             Dnutil(i,j,k) = d2(1)*(nutil_n(i,j+1,k)-nutil_n(i,j,k))
          enddo
       enddo
    endif

    ! Add dissipation to numerical fluxes
    ! -----------------------------------
    if (is_curv) then
       do k=1,nz
          do j=ndy_d_r,nfy_d_r
             do i=1,nx
                ! contravariant velocity at j+1/2 (vc2) and j-1/2 (vc1)
                vc=(vv(i,j,k)*x_ksi(i,j)-uu(i,j,k)*y_ksi(i,j))*ijacob(i,j)
                vc2=(vv(i,j+1,k)*x_ksi(i,j+1)-uu(i,j+1,k)*y_ksi(i,j+1))*ijacob(i,j+1)
                vc1=(vv(i,j-1,k)*x_ksi(i,j-1)-uu(i,j-1,k)*y_ksi(i,j-1))*ijacob(i,j-1)
                ! compute spectral radius at j+1/2 (srp12) and j-1/2 (srm12)
                sr=abs(vc)+c_(i,j,k)*g_eta(i,j)
                srp12=0.5_wp*(abs(vc2)+c_(i,j+1,k)*g_eta(i,j+1)+sr)
                srm12=0.5_wp*(abs(vc1)+c_(i,j-1,k)*g_eta(i,j-1)+sr)
                ! compute dissipation term: delta_j/delta_y=(D_{j+1/2}-D_{j-1/2})/dy
                Knutil(i,j,k) = Knutil(i,j,k) -(srp12*Dnutil(i,j,k) -srm12*Dnutil(i,j-1,k) )
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
                Knutil(i,j,k) = Knutil(i,j,k) -(srp12*Dnutil(i,j,k) -srm12*Dnutil(i,j-1,k))
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
                Knutil(i,j,k) = Knutil(i,j,k) -(srp12*Dnutil(i,j,k) -srm12*Dnutil(i,j-1,k) )*idy(j)
             enddo
          enddo
       enddo
    endif
 
    !*******************
    if (.not.is_2D) then
    !*******************

       ! Artificial viscosity in z-direction [compute Dnutil at k+1/2]
       ! ===================================
       if (is_boundary(3,1)) then
          ! 1st-order centered dissipation
          ! ------------------------------
          k=1
          do j=1,ny
             do i=1,nx
                Dnutil(i,j,k) = d2(1)*(nutil_n(i,j,k+1)-nutil_n(i,j,k))
             enddo
          enddo
       endif
       
       ! Interior points: 3rd-order centered dissipation
       ! -----------------------------------------------
       do k=ndz_di_r,nfz_di_r
          do j=1,ny
             do i=1,nx
                Dnutil(i,j,k) = d4(1)*(nutil_n(i,j,k+1)-nutil_n(i,j,k  )) &
                            + d4(2)*(nutil_n(i,j,k+2)-nutil_n(i,j,k-1))
             enddo
          enddo
       enddo

       if (is_boundary(3,2)) then
          ! 1st-order centered dissipation
          ! ------------------------------
          k=nz-1
          do j=1,ny
             do i=1,nx
                Dnutil(i,j,k) = d2(1)*(nutil_n(i,j,k+1)-nutil_n(i,j,k))
             enddo
          enddo
       endif

       ! Add dissipation to numerical fluxes
       ! -----------------------------------
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
                   Knutil(i,j,k) = Knutil(i,j,k) -(srp12*Dnutil(i,j,k) -srm12*Dnutil(i,j,k-1))
                enddo
             enddo
          enddo 
       else
          do k=ndz_d,nfz_d
             do j=1,ny
                do i=1,nx
                   ! compute spectral radius at k+1/2 (srp12) and k-1/2 (srm12)
                   sr=abs(ww(i,j,k))+c_(i,j,k)
                   srp12=0.5_wp*(abs(ww(i,j,k+1))+c_(i,j,k+1)+sr)
                   srm12=0.5_wp*(abs(ww(i,j,k-1))+c_(i,j,k-1)+sr)
                   ! compute dissipation term: delta_k/delta_z=(D_{k+1/2}-D_{k-1/2})/dz
                   Knutil(i,j,k) = Knutil(i,j,k) -(srp12*Dnutil(i,j,k) -srm12*Dnutil(i,j,k-1) )*idz(k)
                enddo
             enddo
          enddo
       endif

    endif
    
  end subroutine artvisc_o3_SA

  !===============================================================================
  subroutine artvisc_rus_SA
  !===============================================================================
    !> Apply artificial viscosity - Conservative version Rusanov scheme
    !> CM 04/2022 
  !===============================================================================
    use mod_time ! for: deltat
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: vc1,vc2,vc,srm12,srp12,sr
    real(wp), dimension(0:nx,0:ny,0:nz) :: Dnutil
    ! ---------------------------------------------------------------------------

    ! Artificial viscosity in x-direction [compute Dnutil at i+1/2]
    ! ===================================
    ! 1st-order centered dissipation
    ! ------------------------------
    if (is_boundary(1,1)) then
       i=1
       do k=1,nz
          do j=1,ny
             Dnutil(i,j,k) = d2(1)*(nutil_n(i+1,j,k)-nutil_n(i,j,k))
          enddo
       enddo
    endif

    do k=1,nz
       do j=1,ny
          do i=ndx_di_r,nfx_di_r
             Dnutil(i,j,k) = 0.5_wp*(nutil_n(i+1,j,k)-nutil_n(i,j,k))
          enddo
       enddo
    enddo

    if (is_boundary(1,2)) then
       i=nx-1
       do k=1,nz
          do j=1,ny
             Dnutil(i,j,k) = d2(1)*(nutil_n(i+1,j,k)-nutil_n(i,j,k))            
          enddo
       enddo
    endif

    ! Add dissipation to numerical fluxes
    ! -----------------------------------
    if (is_curv) then
       do k=1,nz
          do j=1,ny
             do i=ndx_d_r,nfx_d_r
                ! contravariant velocity at i+1 (vc2) and i-1 (vc1)
                vc =(uu(i  ,j,k)*y_eta_v(i  ,j)-vv(i  ,j,k)*x_eta_v(i  ,j))*ijacob_v(i  ,j)
                vc2=(uu(i+1,j,k)*y_eta_v(i+1,j)-vv(i+1,j,k)*x_eta_v(i+1,j))*ijacob_v(i+1,j)
                vc1=(uu(i-1,j,k)*y_eta_v(i-1,j)-vv(i-1,j,k)*x_eta_v(i-1,j))*ijacob_v(i-1,j)
                ! compute spectral radius at i+1 (srp12) and i-1 (srm12)
                sr   =abs(vc )+c_(i  ,j,k)*g_ksi(i  ,j)
                srp12=abs(vc2)+c_(i+1,j,k)*g_ksi(i+1,j)
                srm12=abs(vc1)+c_(i-1,j,k)*g_ksi(i-1,j)
                ! set sr +/- 1/2 to max value between i and i+1, i and i-1 (Rusanov)
                srp12=max(sr,srp12)
                srm12=max(sr,srm12)
                ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dksi
                Knutil(i,j,k) = Knutil(i,j,k) -(srp12*Dnutil(i,j,k) -srm12*Dnutil(i-1,j,k))
             enddo
          enddo
       enddo
    else
       do k=1,nz
          do j=1,ny
             do i=ndx_d_r,nfx_d_r
                ! compute spectral radius at i+1 (srp1) and i-1 (srm1)
                sr=abs(uu(i,j,k))+c_(i,j,k)
                srp12=abs(uu(i+1,j,k))+c_(i+1,j,k)
                srm12=abs(uu(i-1,j,k))+c_(i-1,j,k)
                ! set sr +/- 1/2 to max value between i and i+1, i and i-1 (Rusanov)
                srp12=max(sr,srp12)
                srm12=max(sr,srm12)
                ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dx
                Knutil(i,j,k) = Knutil(i,j,k) -(srp12*Dnutil(i,j,k) -srm12*Dnutil(i-1,j,k))*idx(i)
             enddo
          enddo
       enddo
    endif

    ! Artificial viscosity in y-direction [compute Dnutil at j+1/2]
    ! ===================================
    ! 1st-order centered dissipation
    ! ------------------------------
    if (is_boundary(2,1)) then
       j=1
       do k=1,nz
          do i=1,nx
             Dnutil(i,j,k) = d2(1)*(nutil_n(i,j+1,k)-nutil_n(i,j,k))
          enddo
       enddo
    endif

    do k=1,nz
       do j=ndy_di_r,nfy_di_r
          do i=1,nx
             Dnutil(i,j,k) = 0.5_wp*(nutil_n(i,j+1,k)-nutil_n(i,j,k))
          enddo
       enddo
    enddo

    if (is_boundary(2,2)) then
       j=ny-1
       do k=1,nz
          do i=1,nx
             Dnutil(i,j,k) = d2(1)*(nutil_n(i,j+1,k)-nutil_n(i,j,k))
          enddo
       enddo
    endif

    ! Add dissipation to numerical fluxes
    ! -----------------------------------
    if (is_curv) then
       do k=1,nz
          do j=ndy_d_r,nfy_d_r
             do i=1,nx
                ! contravariant velocity at j+1 (vc2) and j-1 (vc1)
                vc =(vv(i,j  ,k)*x_ksi_v(i,j  )-uu(i,j  ,k)*y_ksi_v(i,j  ))*ijacob_v(i,j  )
                vc2=(vv(i,j+1,k)*x_ksi_v(i,j+1)-uu(i,j+1,k)*y_ksi_v(i,j+1))*ijacob_v(i,j+1)
                vc1=(vv(i,j-1,k)*x_ksi_v(i,j-1)-uu(i,j-1,k)*y_ksi_v(i,j-1))*ijacob_v(i,j-1)
                ! compute spectral radius at j+1 (srp12) and j-1 (srm12)
                sr   =abs(vc )+c_(i,j  ,k)*g_eta(i,j  )
                srp12=abs(vc2)+c_(i,j+1,k)*g_eta(i,j+1)
                srm12=abs(vc1)+c_(i,j-1,k)*g_eta(i,j-1)
                ! set sr +/- 1/2 to max value between j and j+1, i and j-1 (Rusanov)
                srp12=max(sr,srp12)
                srm12=max(sr,srm12)
                ! compute dissipation term: delta_j/delta_x=(D_{j+1/2}-D_{j-1/2})/deta
                Knutil(i,j,k) = Knutil(i,j,k) -(srp12*Dnutil(i,j,k) -srm12*Dnutil(i,j-1,k))
             enddo
          enddo
       enddo
    else
       do k=1,nz
          do j=ndy_d_r,nfy_d_r
             do i=1,nx
                ! compute spectral radius at j+1 (srp1) and j-1 (srm1)
                sr=abs(vv(i,j,k))+c_(i,j,k)
                srp12=abs(vv(i,j+1,k))+c_(i,j+1,k)
                srm12=abs(vv(i,j-1,k))+c_(i,j-1,k)
                ! set sr +/- 1/2 to max value between j and j+1, i and j-1 (Rusanov)
                srp12=max(sr,srp12)
                srm12=max(sr,srm12)
                ! compute dissipation term: delta_j/delta_x=(D_{j+1/2}-D_{j-1/2})/dx
                Knutil(i,j,k) = Knutil(i,j,k) -(srp12*Dnutil(i,j,k) -srm12*Dnutil(i,j-1,k) )*idx(i)
             enddo
          enddo
       enddo
    endif
 
    !*******************
    if (.not.is_2D) then
    !*******************

       ! Artificial viscosity in z-direction [compute Dnutil at k+1/2]
       ! ===================================
       ! 1st-order centered dissipation
       ! ------------------------------
       do k=ndz_di_r,nfz_di_r
          do j=1,ny
             do i=1,nx
                Dnutil(i,j,k) = 0.5_wp*(nutil_n(i,j,k)-nutil_n(i,j,k+1))            
             enddo
          enddo
       enddo

       ! Add dissipation to numerical fluxes
       ! -----------------------------------
       do k=ndy_d_r,nfy_d_r
          do j=1,ny
             do i=1,nx
                ! compute spectral radius at j+1 (srp1) and j-1 (srm1)
                sr=abs(uu(i,j,k))+c_(i,j,k)
                srp12=abs(uu(i,j,k+1))+c_(i,j,k+1)
                srm12=abs(uu(i,j,k-1))+c_(i,j,k-1)
                ! set sr +/- 1/2 to max value between j and j+1, i and j-1 (Rusanov)
                srp12=max(sr,srp12)
                srm12=max(sr,srm12)
                ! compute dissipation term: delta_j/delta_x=(D_{j+1/2}-D_{j-1/2})/dx
                Knutil(i,j,k) = Knutil(i,j,k) -(srp12*Dnutil(i,j,k) -srm12*Dnutil(i,j,k-1) )*idz(i)
             enddo
          enddo
       enddo
    endif
    
  end subroutine artvisc_rus_SA

    !**************************
    ! SA-noft2-Gamma-Re_theta *
    !**************************
  
    !===============================================================================
  subroutine artvisc_o3_SA_transition
    !===============================================================================
      !> Apply artificial viscosity - Conservative version DNC order 3
      !> Modified version XG 06/2021
    !===============================================================================
      use mod_time ! for: deltat
      implicit none
      ! ---------------------------------------------------------------------------
      integer :: i,j,k
      real(wp) :: vc1,vc2,vc,srm12,srp12,sr
      real(wp), dimension(0:nx,0:ny,0:nz) :: Dnutil, Dgamma, Dre_theta
      ! ---------------------------------------------------------------------------
  
      ! Artificial viscosity in x-direction [compute Dnutil at i+1/2]
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
               Dnutil(i,j,k)   = d2(1)*(nutil_n(i+1,j,k)-nutil_n(i,j,k))   
               Dgamma(i,j,k) = d2(1)*(rhogamma_n(i+1,j,k)-rhogamma_n(i,j,k))            
               Dre_theta(i,j,k) = d2(1)*(rhore_theta_n(i+1,j,k)-rhore_theta_n(i,j,k))       
            enddo
         enddo
      endif
  
      ! Interior points: 3rd-order centered dissipation
      ! -----------------------------------------------
      do k=1,nz
         do j=1,ny
            do i=ndx_di_r,nfx_di_r
               Dnutil(i,j,k) = d4(1)*(nutil_n(i+1,j,k)-nutil_n(i  ,j,k)) &
                                + d4(2)*(nutil_n(i+2,j,k)-nutil_n(i-1,j,k))
               Dgamma(i,j,k) =    d4(1)*(rhogamma_n(i+1,j,k)-rhogamma_n(i  ,j,k)) &
                                + d4(2)*(rhogamma_n(i+2,j,k)-rhogamma_n(i-1,j,k))
               Dre_theta(i,j,k) = d4(1)*(rhore_theta_n(i+1,j,k)-rhore_theta_n(i  ,j,k)) &
                                + d4(2)*(rhore_theta_n(i+2,j,k)-rhore_theta_n(i-1,j,k))
            enddo
         enddo
      enddo
  
      if (is_boundary(1,2)) then
         ! 1st-order centered dissipation
         ! ------------------------------
         i=nx-1
         do k=1,nz
            do j=1,ny
               Dnutil(i,j,k) = d2(1)*(nutil_n(i+1,j,k)-nutil_n(i,j,k))            
               Dgamma(i,j,k) = d2(1)*(rhogamma_n(i+1,j,k)-rhogamma_n(i,j,k))            
               Dre_theta(i,j,k) = d2(1)*(rhore_theta_n(i+1,j,k)-rhore_theta_n(i,j,k))            
            enddo
         enddo
  
         ! no dissipation at i=nx -> initialize increment array
         do j=1,ny
            Knutil(nx,j,:)=0.0_wp
            Kgamma(nx,j,:)=0.0_wp
            Kre_theta(nx,j,:)=0.0_wp
         enddo
      endif
      
      ! Add dissipation to numerical fluxes
      ! -----------------------------------
      if (is_curv) then
         do k=1,nz
            do j=1,ny
               do i=ndx_d_r,nfx_d_r
                  ! contravariant velocity at i+1/2 (vc2) and i-1/2 (vc1)
                  vc=(uu(i,j,k)*y_eta(i,j)-vv(i,j,k)*x_eta(i,j))*ijacob(i,j)
                  vc2=(uu(i+1,j,k)*y_eta(i+1,j)-vv(i+1,j,k)*x_eta(i+1,j))*ijacob(i+1,j)
                  vc1=(uu(i-1,j,k)*y_eta(i-1,j)-vv(i-1,j,k)*x_eta(i-1,j))*ijacob(i-1,j)
                  ! compute spectral radius at i+1/2 (srp12) and i-1/2 (srm12)
                  sr=abs(vc)+c_(i,j,k)*g_ksi(i,j)
                  srp12=0.5_wp*(abs(vc2)+c_(i+1,j,k)*g_ksi(i+1,j)+sr)
                  srm12=0.5_wp*(abs(vc1)+c_(i-1,j,k)*g_ksi(i-1,j)+sr)
                  ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dksi
                  Knutil(i,j,k) = Knutil(i,j,k) -(srp12*Dnutil(i,j,k) -srm12*Dnutil(i-1,j,k))
                  Kgamma(i,j,k) = Kgamma(i,j,k) -(srp12*Dgamma(i,j,k) -srm12*Dgamma(i-1,j,k))
                  Kre_theta(i,j,k) = Kre_theta(i,j,k) -(srp12*Dre_theta(i,j,k) -srm12*Dre_theta(i-1,j,k))
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
                  Knutil(i,j,k) = Knutil(i,j,k) -(srp12*Dnutil(i,j,k) -srm12*Dnutil(i-1,j,k))
                  Kgamma(i,j,k) = Kgamma(i,j,k) -(srp12*Dgamma(i,j,k) -srm12*Dgamma(i-1,j,k))
                  Kre_theta(i,j,k) = Kre_theta(i,j,k) -(srp12*Dre_theta(i,j,k) -srm12*Dre_theta(i-1,j,k))
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
                  Knutil(i,j,k) = Knutil(i,j,k) -(srp12*Dnutil(i,j,k) -srm12*Dnutil(i-1,j,k) )*idx(i)
                  Kgamma(i,j,k) = Kgamma(i,j,k) -(srp12*Dgamma(i,j,k) -srm12*Dgamma(i-1,j,k) )*idx(i)
                  Kre_theta(i,j,k) = Kre_theta(i,j,k) -(srp12*Dre_theta(i,j,k) -srm12*Dre_theta(i-1,j,k) )*idx(i)
               enddo
            enddo
         enddo
      endif
  
      ! Artificial viscosity in y-direction [compute Dnutil at j+1/2]
      ! ===================================
      if (is_boundary(2,1)) then
         ! 1st-order centered dissipation
         ! ------------------------------
         j=1
         do k=1,nz
            do i=1,nx
               Dnutil(i,j,k) = d2(1)*(nutil_n(i,j+1,k)-nutil_n(i,j,k))
               Dgamma(i,j,k) = d2(1)*(rhogamma_n(i,j+1,k)-rhogamma_n(i,j,k))
               Dre_theta(i,j,k) = d2(1)*(rhore_theta_n(i,j+1,k)-rhore_theta_n(i,j,k))
            enddo
         enddo
      endif
      
      ! Interior points: 3rd-order centered dissipation
      ! -----------------------------------------------
      do k=1,nz
         do j=ndy_di_r,nfy_di_r
            do i=1,nx
               Dnutil(i,j,k) = d4(1)*(nutil_n(i,j+1,k)-nutil_n(i,j  ,k)) &
                           + d4(2)*(nutil_n(i,j+2,k)-nutil_n(i,j-1,k))
               Dgamma(i,j,k) = d4(1)*(rhogamma_n(i,j+1,k)-rhogamma_n(i,j  ,k)) &
                           + d4(2)*(rhogamma_n(i,j+2,k)-rhogamma_n(i,j-1,k))
               Dre_theta(i,j,k) = d4(1)*(rhore_theta_n(i,j+1,k)-rhore_theta_n(i,j  ,k)) &
                           + d4(2)*(rhore_theta_n(i,j+2,k)-rhore_theta_n(i,j-1,k))
            enddo
         enddo
      enddo
  
      if (is_boundary(2,2)) then
         ! 1st-order centered dissipation
         ! ------------------------------
         j=ny-1
         do k=1,nz
            do i=1,nx
               Dnutil(i,j,k) = d2(1)*(nutil_n(i,j+1,k)-nutil_n(i,j,k))
               Dgamma(i,j,k) = d2(1)*(rhogamma_n(i,j+1,k)-rhogamma_n(i,j,k))
               Dre_theta(i,j,k) = d2(1)*(rhore_theta_n(i,j+1,k)-rhore_theta_n(i,j,k))
            enddo
         enddo
      endif
  
      ! Add dissipation to numerical fluxes
      ! -----------------------------------
      if (is_curv) then
         do k=1,nz
            do j=ndy_d_r,nfy_d_r
               do i=1,nx
                  ! contravariant velocity at j+1/2 (vc2) and j-1/2 (vc1)
                  vc=(vv(i,j,k)*x_ksi(i,j)-uu(i,j,k)*y_ksi(i,j))*ijacob(i,j)
                  vc2=(vv(i,j+1,k)*x_ksi(i,j+1)-uu(i,j+1,k)*y_ksi(i,j+1))*ijacob(i,j+1)
                  vc1=(vv(i,j-1,k)*x_ksi(i,j-1)-uu(i,j-1,k)*y_ksi(i,j-1))*ijacob(i,j-1)
                  ! compute spectral radius at j+1/2 (srp12) and j-1/2 (srm12)
                  sr=abs(vc)+c_(i,j,k)*g_eta(i,j)
                  srp12=0.5_wp*(abs(vc2)+c_(i,j+1,k)*g_eta(i,j+1)+sr)
                  srm12=0.5_wp*(abs(vc1)+c_(i,j-1,k)*g_eta(i,j-1)+sr)
                  ! compute dissipation term: delta_j/delta_y=(D_{j+1/2}-D_{j-1/2})/dy
                  Knutil(i,j,k) = Knutil(i,j,k) -(srp12*Dnutil(i,j,k) -srm12*Dnutil(i,j-1,k) )
                  Kgamma(i,j,k) = Kgamma(i,j,k) -(srp12*Dgamma(i,j,k) -srm12*Dgamma(i,j-1,k) )
                  Kre_theta(i,j,k) = Kre_theta(i,j,k) -(srp12*Dre_theta(i,j,k) -srm12*Dre_theta(i,j-1,k) )
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
                  Knutil(i,j,k) = Knutil(i,j,k) -(srp12*Dnutil(i,j,k) -srm12*Dnutil(i,j-1,k))
                  Kgamma(i,j,k) = Kgamma(i,j,k) -(srp12*Dgamma(i,j,k) -srm12*Dgamma(i,j-1,k))
                  Kre_theta(i,j,k) = Kre_theta(i,j,k) -(srp12*Dre_theta(i,j,k) -srm12*Dre_theta(i,j-1,k))
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
                  Knutil(i,j,k) = Knutil(i,j,k) -(srp12*Dnutil(i,j,k) -srm12*Dnutil(i,j-1,k) )*idy(j)
                  Kgamma(i,j,k) = Kgamma(i,j,k) -(srp12*Dgamma(i,j,k) -srm12*Dgamma(i,j-1,k) )*idy(j)
                  Kre_theta(i,j,k) = Kre_theta(i,j,k) -(srp12*Dre_theta(i,j,k) -srm12*Dre_theta(i,j-1,k) )*idy(j)
               enddo
            enddo
         enddo
      endif
   
      !*******************
      if (.not.is_2D) then
      !*******************
  
         ! Artificial viscosity in z-direction [compute Dnutil at k+1/2]
         ! ===================================
         if (is_boundary(3,1)) then
            ! 1st-order centered dissipation
            ! ------------------------------
            k=1
            do j=1,ny
               do i=1,nx
                  Dnutil(i,j,k) = d2(1)*(nutil_n(i,j,k+1)-nutil_n(i,j,k))
                  Dgamma(i,j,k) = d2(1)*(rhogamma_n(i,j,k+1)-rhogamma_n(i,j,k))
                  Dre_theta(i,j,k) = d2(1)*(rhore_theta_n(i,j,k+1)-rhore_theta_n(i,j,k))
               enddo
            enddo
         endif
         
         ! Interior points: 3rd-order centered dissipation
         ! -----------------------------------------------
         do k=ndz_di_r,nfz_di_r
            do j=1,ny
               do i=1,nx
                  Dnutil(i,j,k) = d4(1)*(nutil_n(i,j,k+1)-nutil_n(i,j,k  )) &
                              + d4(2)*(nutil_n(i,j,k+2)-nutil_n(i,j,k-1))
                  Dgamma(i,j,k) = d4(1)*(rhogamma_n(i,j,k+1)-rhogamma_n(i,j,k  )) &
                              + d4(2)*(rhogamma_n(i,j,k+2)-rhogamma_n(i,j,k-1))
                  Dre_theta(i,j,k) = d4(1)*(rhore_theta_n(i,j,k+1)-rhore_theta_n(i,j,k  )) &
                              + d4(2)*(rhore_theta_n(i,j,k+2)-rhore_theta_n(i,j,k-1))
               enddo
            enddo
         enddo
  
         if (is_boundary(3,2)) then
            ! 1st-order centered dissipation
            ! ------------------------------
            k=nz-1
            do j=1,ny
               do i=1,nx
                  Dnutil(i,j,k) = d2(1)*(nutil_n(i,j,k+1)-nutil_n(i,j,k))
                  Dgamma(i,j,k) = d2(1)*(rhogamma_n(i,j,k+1)-rhogamma_n(i,j,k))
                  Dre_theta(i,j,k) = d2(1)*(rhore_theta_n(i,j,k+1)-rhore_theta_n(i,j,k))
               enddo
            enddo
         endif
  
         ! Add dissipation to numerical fluxes
         ! -----------------------------------
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
                     Knutil(i,j,k) = Knutil(i,j,k) -(srp12*Dnutil(i,j,k) -srm12*Dnutil(i,j,k-1))
                     Kgamma(i,j,k) = Kgamma(i,j,k) -(srp12*Dgamma(i,j,k) -srm12*Dgamma(i,j,k-1))
                     Kre_theta(i,j,k) = Kre_theta(i,j,k) -(srp12*Dre_theta(i,j,k) -srm12*Dre_theta(i,j,k-1))
                  enddo
               enddo
            enddo 
         else
            do k=ndz_d,nfz_d
               do j=1,ny
                  do i=1,nx
                     ! compute spectral radius at k+1/2 (srp12) and k-1/2 (srm12)
                     sr=abs(ww(i,j,k))+c_(i,j,k)
                     srp12=0.5_wp*(abs(ww(i,j,k+1))+c_(i,j,k+1)+sr)
                     srm12=0.5_wp*(abs(ww(i,j,k-1))+c_(i,j,k-1)+sr)
                     ! compute dissipation term: delta_k/delta_z=(D_{k+1/2}-D_{k-1/2})/dz
                     Knutil(i,j,k) = Knutil(i,j,k) -(srp12*Dnutil(i,j,k) -srm12*Dnutil(i,j,k-1) )*idz(k)
                     Kgamma(i,j,k) = Kgamma(i,j,k) -(srp12*Dgamma(i,j,k) -srm12*Dgamma(i,j,k-1) )*idz(k)
                     Kre_theta(i,j,k) = Kre_theta(i,j,k) -(srp12*Dre_theta(i,j,k) -srm12*Dre_theta(i,j,k-1) )*idz(k)
                  enddo
               enddo
            enddo
         endif
  
      endif
      
    end subroutine artvisc_o3_SA_transition
  
    !===============================================================================
    subroutine artvisc_rus_SA_transition
    !===============================================================================
      !> Apply artificial viscosity - Conservative version Rusanov scheme
      !> CM 04/2022 
    !===============================================================================
      use mod_time ! for: deltat
      implicit none
      ! ---------------------------------------------------------------------------
      integer :: i,j,k
      real(wp) :: vc1,vc2,vc,srm12,srp12,sr
      real(wp), dimension(0:nx,0:ny,0:nz) :: Dnutil, Dgamma, Dre_theta
      ! ---------------------------------------------------------------------------
  
      ! Artificial viscosity in x-direction [compute Dnutil at i+1/2]
      ! ===================================
      ! 1st-order centered dissipation
      ! ------------------------------
      if (is_boundary(1,1)) then
         i=1
         do k=1,nz
            do j=1,ny
               Dnutil(i,j,k) = d2(1)*(nutil_n(i+1,j,k)-nutil_n(i,j,k))
               Dgamma(i,j,k) = d2(1)*(rhogamma_n(i+1,j,k)-rhogamma_n(i,j,k))
               Dre_theta(i,j,k) = d2(1)*(rhore_theta_n(i+1,j,k)-rhore_theta_n(i,j,k))
            enddo
         enddo
      endif
  
      do k=1,nz
         do j=1,ny
            do i=ndx_di_r,nfx_di_r
               Dnutil(i,j,k) = 0.5_wp*(nutil_n(i+1,j,k)-nutil_n(i,j,k))
               Dgamma(i,j,k) = 0.5_wp*(rhogamma_n(i+1,j,k)-rhogamma_n(i,j,k))
               Dre_theta(i,j,k) = 0.5_wp*(rhore_theta_n(i+1,j,k)-rhore_theta_n(i,j,k))
            enddo
         enddo
      enddo
  
      if (is_boundary(1,2)) then
         i=nx-1
         do k=1,nz
            do j=1,ny
               Dnutil(i,j,k) = d2(1)*(nutil_n(i+1,j,k)-nutil_n(i,j,k))            
               Dgamma(i,j,k) = d2(1)*(rhogamma_n(i+1,j,k)-rhogamma_n(i,j,k))            
               Dre_theta(i,j,k) = d2(1)*(rhore_theta_n(i+1,j,k)-rhore_theta_n(i,j,k))            
            enddo
         enddo
      endif
  
      ! Add dissipation to numerical fluxes
      ! -----------------------------------
      if (is_curv) then
         do k=1,nz
            do j=1,ny
               do i=ndx_d_r,nfx_d_r
                  ! contravariant velocity at i+1 (vc2) and i-1 (vc1)
                  vc =(uu(i  ,j,k)*y_eta_v(i  ,j)-vv(i  ,j,k)*x_eta_v(i  ,j))*ijacob_v(i  ,j)
                  vc2=(uu(i+1,j,k)*y_eta_v(i+1,j)-vv(i+1,j,k)*x_eta_v(i+1,j))*ijacob_v(i+1,j)
                  vc1=(uu(i-1,j,k)*y_eta_v(i-1,j)-vv(i-1,j,k)*x_eta_v(i-1,j))*ijacob_v(i-1,j)
                  ! compute spectral radius at i+1 (srp12) and i-1 (srm12)
                  sr   =abs(vc )+c_(i  ,j,k)*g_ksi(i  ,j)
                  srp12=abs(vc2)+c_(i+1,j,k)*g_ksi(i+1,j)
                  srm12=abs(vc1)+c_(i-1,j,k)*g_ksi(i-1,j)
                  ! set sr +/- 1/2 to max value between i and i+1, i and i-1 (Rusanov)
                  srp12=max(sr,srp12)
                  srm12=max(sr,srm12)
                  ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dksi
                  Knutil(i,j,k) = Knutil(i,j,k) -(srp12*Dnutil(i,j,k) -srm12*Dnutil(i-1,j,k))
                  Kgamma(i,j,k) = Kgamma(i,j,k) -(srp12*Dgamma(i,j,k) -srm12*Dgamma(i-1,j,k))
                  Kre_theta(i,j,k) = Kre_theta(i,j,k) -(srp12*Dre_theta(i,j,k) -srm12*Dre_theta(i-1,j,k))
               enddo
            enddo
         enddo
      else
         do k=1,nz
            do j=1,ny
               do i=ndx_d_r,nfx_d_r
                  ! compute spectral radius at i+1 (srp1) and i-1 (srm1)
                  sr=abs(uu(i,j,k))+c_(i,j,k)
                  srp12=abs(uu(i+1,j,k))+c_(i+1,j,k)
                  srm12=abs(uu(i-1,j,k))+c_(i-1,j,k)
                  ! set sr +/- 1/2 to max value between i and i+1, i and i-1 (Rusanov)
                  srp12=max(sr,srp12)
                  srm12=max(sr,srm12)
                  ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dx
                  Knutil(i,j,k) = Knutil(i,j,k) -(srp12*Dnutil(i,j,k) -srm12*Dnutil(i-1,j,k))*idx(i)
                  Kgamma(i,j,k) = Kgamma(i,j,k) -(srp12*Dgamma(i,j,k) -srm12*Dgamma(i-1,j,k))*idx(i)
                  Kre_theta(i,j,k) = Kre_theta(i,j,k) -(srp12*Dre_theta(i,j,k) -srm12*Dre_theta(i-1,j,k))*idx(i)
               enddo
            enddo
         enddo
      endif
  
      ! Artificial viscosity in y-direction [compute Dnutil at j+1/2]
      ! ===================================
      ! 1st-order centered dissipation
      ! ------------------------------
      if (is_boundary(2,1)) then
         j=1
         do k=1,nz
            do i=1,nx
               Dnutil(i,j,k) = d2(1)*(nutil_n(i,j+1,k)-nutil_n(i,j,k))
               Dgamma(i,j,k) = d2(1)*(rhogamma_n(i,j+1,k)-rhogamma_n(i,j,k))
               Dre_theta(i,j,k) = d2(1)*(rhore_theta_n(i,j+1,k)-rhore_theta_n(i,j,k))
            enddo
         enddo
      endif
  
      do k=1,nz
         do j=ndy_di_r,nfy_di_r
            do i=1,nx
               Dnutil(i,j,k) = 0.5_wp*(nutil_n(i,j+1,k)-nutil_n(i,j,k))
               Dgamma(i,j,k) = 0.5_wp*(rhogamma_n(i,j+1,k)-rhogamma_n(i,j,k))
               Dre_theta(i,j,k) = 0.5_wp*(rhore_theta_n(i,j+1,k)-rhore_theta_n(i,j,k))
            enddo
         enddo
      enddo
  
      if (is_boundary(2,2)) then
         j=ny-1
         do k=1,nz
            do i=1,nx
               Dnutil(i,j,k) = d2(1)*(nutil_n(i,j+1,k)-nutil_n(i,j,k))
               Dgamma(i,j,k) = d2(1)*(rhogamma_n(i,j+1,k)-rhogamma_n(i,j,k))
               Dre_theta(i,j,k) = d2(1)*(rhore_theta_n(i,j+1,k)-rhore_theta_n(i,j,k))
            enddo
         enddo
      endif
  
      ! Add dissipation to numerical fluxes
      ! -----------------------------------
      if (is_curv) then
         do k=1,nz
            do j=ndy_d_r,nfy_d_r
               do i=1,nx
                  ! contravariant velocity at j+1 (vc2) and j-1 (vc1)
                  vc =(vv(i,j  ,k)*x_ksi_v(i,j  )-uu(i,j  ,k)*y_ksi_v(i,j  ))*ijacob_v(i,j  )
                  vc2=(vv(i,j+1,k)*x_ksi_v(i,j+1)-uu(i,j+1,k)*y_ksi_v(i,j+1))*ijacob_v(i,j+1)
                  vc1=(vv(i,j-1,k)*x_ksi_v(i,j-1)-uu(i,j-1,k)*y_ksi_v(i,j-1))*ijacob_v(i,j-1)
                  ! compute spectral radius at j+1 (srp12) and j-1 (srm12)
                  sr   =abs(vc )+c_(i,j  ,k)*g_eta(i,j  )
                  srp12=abs(vc2)+c_(i,j+1,k)*g_eta(i,j+1)
                  srm12=abs(vc1)+c_(i,j-1,k)*g_eta(i,j-1)
                  ! set sr +/- 1/2 to max value between j and j+1, i and j-1 (Rusanov)
                  srp12=max(sr,srp12)
                  srm12=max(sr,srm12)
                  ! compute dissipation term: delta_j/delta_x=(D_{j+1/2}-D_{j-1/2})/deta
                  Knutil(i,j,k) = Knutil(i,j,k) -(srp12*Dnutil(i,j,k) -srm12*Dnutil(i,j-1,k))
                  Kgamma(i,j,k) = Kgamma(i,j,k) -(srp12*Dgamma(i,j,k) -srm12*Dgamma(i,j-1,k))
                  Kre_theta(i,j,k) = Kre_theta(i,j,k) -(srp12*Dre_theta(i,j,k) -srm12*Dre_theta(i,j-1,k))
               enddo
            enddo
         enddo
      else
         do k=1,nz
            do j=ndy_d_r,nfy_d_r
               do i=1,nx
                  ! compute spectral radius at j+1 (srp1) and j-1 (srm1)
                  sr=abs(vv(i,j,k))+c_(i,j,k)
                  srp12=abs(vv(i,j+1,k))+c_(i,j+1,k)
                  srm12=abs(vv(i,j-1,k))+c_(i,j-1,k)
                  ! set sr +/- 1/2 to max value between j and j+1, i and j-1 (Rusanov)
                  srp12=max(sr,srp12)
                  srm12=max(sr,srm12)
                  ! compute dissipation term: delta_j/delta_x=(D_{j+1/2}-D_{j-1/2})/dx
                  Knutil(i,j,k) = Knutil(i,j,k) -(srp12*Dnutil(i,j,k) -srm12*Dnutil(i,j-1,k) )*idx(i)
                  Kgamma(i,j,k) = Kgamma(i,j,k) -(srp12*Dgamma(i,j,k) -srm12*Dgamma(i,j-1,k) )*idx(i)
                  Kre_theta(i,j,k) = Kre_theta(i,j,k) -(srp12*Dre_theta(i,j,k) -srm12*Dre_theta(i,j-1,k) )*idx(i)
               enddo
            enddo
         enddo
      endif
   
      !*******************
      if (.not.is_2D) then
      !*******************
  
         ! Artificial viscosity in z-direction [compute Dnutil at k+1/2]
         ! ===================================
         ! 1st-order centered dissipation
         ! ------------------------------
         do k=ndz_di_r,nfz_di_r
            do j=1,ny
               do i=1,nx
                  Dnutil(i,j,k) = 0.5_wp*(nutil_n(i,j,k)-nutil_n(i,j,k+1))    
                  Dgamma(i,j,k) = 0.5_wp*(rhogamma_n(i,j,k)-rhogamma_n(i,j,k+1))    
                  Dre_theta(i,j,k) = 0.5_wp*(rhore_theta_n(i,j,k)-rhore_theta_n(i,j,k+1))            
               enddo
            enddo
         enddo
  
         ! Add dissipation to numerical fluxes
         ! -----------------------------------
         do k=ndy_d_r,nfy_d_r
            do j=1,ny
               do i=1,nx
                  ! compute spectral radius at j+1 (srp1) and j-1 (srm1)
                  sr=abs(uu(i,j,k))+c_(i,j,k)
                  srp12=abs(uu(i,j,k+1))+c_(i,j,k+1)
                  srm12=abs(uu(i,j,k-1))+c_(i,j,k-1)
                  ! set sr +/- 1/2 to max value between j and j+1, i and j-1 (Rusanov)
                  srp12=max(sr,srp12)
                  srm12=max(sr,srm12)
                  ! compute dissipation term: delta_j/delta_x=(D_{j+1/2}-D_{j-1/2})/dx
                  Knutil(i,j,k) = Knutil(i,j,k) -(srp12*Dnutil(i,j,k) -srm12*Dnutil(i,j,k-1) )*idz(i)
                  Kgamma(i,j,k) = Kgamma(i,j,k) -(srp12*Dgamma(i,j,k) -srm12*Dgamma(i,j,k-1) )*idz(i)
                  Kre_theta(i,j,k) = Kre_theta(i,j,k) -(srp12*Dre_theta(i,j,k) -srm12*Dre_theta(i,j,k-1) )*idz(i)
               enddo
            enddo
         enddo
      endif
      
    end subroutine artvisc_rus_SA_transition  

end module mod_artvisc_rans
