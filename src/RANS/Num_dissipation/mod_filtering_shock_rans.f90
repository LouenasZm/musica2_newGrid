!==============================================================================
module mod_filtering_shock_rans
!==============================================================================
  !> Module for filtering with Jameson-type shock capturing
!==============================================================================
  use mod_flow
  use mod_artvisc
  use mod_filtering_shock
  implicit none
  ! ---------------------------------------------------------------------------
  ! amplitude dissipation
  real(wp), private :: xnu0 ! low-order smoothing term
  real(wp), private :: xnu4,xnu2
  ! ---------------------------------------------------------------------------

contains

  !============================================================================
  subroutine init_filter_shock_rans(dissip_coeff,dissip_shock)
  !============================================================================
    !> Initialization of filtering with Jameson-type shock capturing
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    real(wp), intent(in) :: dissip_coeff,dissip_shock
    real(wp) :: dissip_coeff_rans,dissip_shock_rans
    ! -------------------------------------------------------------------------

!!$    ! Filtering coefficients (DNC-Jameson)
!!$    ! =================================
!!$    ! order 1
!!$    d2(1)= 1.0_wp
!!$
!!$    ! order 3
!!$    d4(1)= 3.0_wp
!!$    d4(2)=-1.0_wp
!!$
!!$    ! order 5
!!$    d6(1)= 10.0_wp
!!$    d6(2)=-5.0_wp
!!$    d6(3)= 1.0_wp
!!$
!!$    ! order 7
!!$    d8(1)= 35.0_wp
!!$    d8(2)=-21.0_wp
!!$    d8(3)=  7.0_wp
!!$    d8(4)= -1.0_wp
!!$
!!$    ! order 9
!!$    d10(1)= 126.0_wp
!!$    d10(2)=-84.0_wp
!!$    d10(3)= 36.0_wp
!!$    d10(4)=-9.0_wp
!!$    d10(5)= 1.0_wp
!!$    
!!$    ! Dissipation coefficients [from param.ini]
!!$    ! ========================
!!$    xnu10=dissip_coeff/1260.0_wp
!!$    xnu8= dissip_coeff/280.0_wp
!!$    xnu6= dissip_coeff/60.0_wp
!!$    xnu4= dissip_coeff/12.0_wp
!!$    xnu2= 0.0_wp*dissip_coeff/2.0_wp
!!$
!!$    ! low-order dissipation term [from param.ini]
!!$    ! ==========================
!!$    xnu0=dissip_shock



    ! TEMPORARY only for RANS
    dissip_coeff_rans=0.1_wp
    dissip_shock_rans=0.1_wp
   
    chi=dissip_coeff_rans

    ! Filtering coefficients (DNC-Jameson)
    ! =================================
    ! /!\ Change sign because added to increment
    ! order 1
    d2(1)=-1.0_wp

    ! order 3
    d4(1)=-3.0_wp
    d4(2)= 1.0_wp

    ! Dissipation coefficients [from param.ini]
    ! ========================
    xnu4= dissip_coeff_rans/12.0_wp
    xnu2= 0.0_wp*dissip_coeff_rans/2.0_wp

    ! low-order dissipation term [from param.ini]
    ! ==========================
    ! /!\ Change sign because added to increment
    xnu0=-dissip_shock_rans

    ! Options TO BE CHANGED ---> in param.ini
    ! =======
    
    ! Ducros sensor
    ! -------------
    is_ducros=.true.

    ! Pressure sensor
    ! ---------------
    ! -> mix of Jameson's sensor and TVD version of Swanson et al. (1998)
    ! csens=0.01_wp --> activate Jameson's sensor
    ! csens=0.5_wp  --> mix of Jameson/TVD sensors [recommended]
    ! csens=1.0_wp  --> full TVD sensor
    csens=0.5_wp
    
    ! Switch toward artificial dissipation-like coefficient (Edoh et al., JCP, 2018)
    ! -----------------------------------------------------
    is_sw_edoh=.false.
    !is_sw_edoh=.true.
    
  end subroutine init_filter_shock_rans

  !==============================================================================
  subroutine filter_o4_shock_SA
  !==============================================================================
    !> Apply 4th-order filtering with Jameson-type shock capturing
    !> (Conservative version inspired from DNC order 9)
    !> modified version XG 2023
  !==============================================================================
    use mod_time
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: vc,cfll,mul
    real(wp), dimension(0:nx,0:ny,0:nz) :: Dnutil
    ! ---------------------------------------------------------------------------
    real(wp) :: div2,om1,om2,om3
    real(wp) :: xk0,xk2,xk4
    real(wp) :: psens_trad,psens_tvd
    real(wp), dimension(ndx_d1:nfx_d1,ndy_d1:nfy_d1,ndz_d1:nfz_d1) :: sens,psens
    ! ---------------------------------------------------------------------------

    ! Compute Ducros sensor
    ! =====================
    if (is_ducros) then
       do k=ndz_d1,nfz_d1
          do j=ndy_d1,nfy_d1
             do i=ndx_d1,nfx_d1
                div2=(dux(i,j,k)+dvy(i,j,k)+dwz(i,j,k))**2
                om1 =(dwy(i,j,k)-dvz(i,j,k))**2
                om2 =(duz(i,j,k)-dwx(i,j,k))**2
                om3 =(dvx(i,j,k)-duy(i,j,k))**2
                ! Ducros' sensor
                sens(i,j,k)=div2/(div2+om1+om2+om3+cutoff)
                !uvar(i,j,k,1)=sens(i,j,k)
             enddo
          enddo
       enddo
    else
       sens=1.0_wp
    endif
    
    ! Compute pressure sensor along i-direction
    ! =========================================
    do k=1,nz
       do j=1,ny
          do i=ndx_d1,nfx_d1
!              psens(i,j,k)=sens(i,j,k)*abs(prs(i-1,j,k)-2.0_wp*prs(i,j,k)+prs(i+1,j,k)) / &
!                                       abs(prs(i-1,j,k)+2.0_wp*prs(i,j,k)+prs(i+1,j,k))

             ! an alternative sensor (implemented by OY)
             ! a mix of a TVD concept and the traditional one (suggested by Swanson et al. [1998])
             psens_trad = prs(i-1,j,k)+2.0_wp*prs(i,j,k)+prs(i+1,j,k)
             psens_tvd = abs(prs(i+1,j,k)-prs(i,j,k))+abs(prs(i,j,k)-prs(i-1,j,k))
             psens(i,j,k)=sens(i,j,k)*abs(prs(i-1,j,k)-2.0_wp*prs(i,j,k)+prs(i+1,j,k)) / &
                                      (csens*psens_trad + (1.0_wp-csens)*psens_tvd)

          enddo
       enddo
    enddo

    ! Filtering in i-direction [compute Dnutil at i+1/2]
    ! ========================

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
       do j=1,ny
           Knutil(nx,j,:)=0.0_wp
       enddo
    endif

    ! Compute filtering terms
    ! ------------------------
    if (is_sw_edoh) then
       if (is_curv3) then
          do k=1,nz
             do j=1,ny
                do i=ndx_d,nfx_d
                   ! contravariant velocity
                   vc=(uu(i,j,k)*ksi_x(i,j,k)+vv(i,j,k)*ksi_y(i,j,k)+ww(i,j,k)*ksi_z(i,j,k))*ijacob3(i,j,k)
                   ! compute local CFL
                   cfll=(abs(vc)+c_(i,j,k)*g3_ksi(i,j,k))*deltat
                   ! switch Edoh et al.
                   mul=min(1.0_wp,cfll/chi)
                   ! If call to local_cfl
                   !mul=min(1.0_wp,cfl_j(i,j,k)/chi)

                   ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dx
                   Knutil(i,j,k) = mul*(Dnutil(i-1,j,k) -Dnutil(i,j,k))
                enddo
             enddo
          enddo
       elseif (is_curv) then
          do k=1,nz
             do j=1,ny
                do i=ndx_d,nfx_d
                   ! contravariant velocity
                   vc=(uu(i,j,k)*y_eta(i,j)-vv(i,j,k)*x_eta(i,j))*ijacob(i,j)
                   ! compute local CFL
                   cfll=(abs(vc)+c_(i,j,k)*g_ksi(i,j))*deltat
                   ! switch Edoh et al.
                   mul=min(1.0_wp,cfll/chi)
                   ! If call to local_cfl
                   !mul=min(1.0_wp,cfl_j(i,j,k)/chi)

                   ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dx
                   Knutil(i,j,k) = mul*(Dnutil(i-1,j,k) -Dnutil(i,j,k))
                enddo
             enddo
          enddo
       else
          do k=1,nz
             do j=1,ny
                do i=ndx_d,nfx_d
                   ! compute local CFL
                   cfll=(abs(uu(i,j,k))+c_(i,j,k))*deltat*idx(i)
                   ! switch Edoh et al.
                   mul=min(1.0_wp,cfll/chi)!(chi,cfll)
                   
                   ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dx
                   Knutil(i,j,k) = mul*(Dnutil(i-1,j,k) -Dnutil(i,j,k))
                enddo
             enddo
          enddo
       endif
    else
       do k=1,nz
          do j=1,ny
             do i=ndx_d,nfx_d
                ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dx
                Knutil(i,j,k) = Dnutil(i-1,j,k) -Dnutil(i,j,k)
             enddo
          enddo
       enddo
    endif

    ! Compute pressure sensor along j-direction
    ! =========================================
    do k=1,nz
       do j=ndy_d1,nfy_d1
          do i=1,nx
!              psens(i,j,k)=sens(i,j,k)*abs(prs(i,j-1,k)-2.0_wp*prs(i,j,k)+prs(i,j+1,k)) / &
!                                       abs(prs(i,j-1,k)+2.0_wp*prs(i,j,k)+prs(i,j+1,k))

             ! an alternative sensor (implemented by OY)
             ! a mix of a TVD concept and the traditional one (suggested by Swanson et al. [1998])
             psens_trad = prs(i,j-1,k)+2.0_wp*prs(i,j,k)+prs(i,j+1,k)
             psens_tvd = abs(prs(i,j+1,k)-prs(i,j,k))+abs(prs(i,j,k)-prs(i,j-1,k))
             psens(i,j,k)=sens(i,j,k)*abs(prs(i,j-1,k)-2.0_wp*prs(i,j,k)+prs(i,j+1,k)) / &
                                      (csens*psens_trad + (1.0_wp-csens)*psens_tvd)
          enddo
       enddo
    enddo

    ! Filtering in y-direction [compute Dnutil at j+1/2]
    ! ========================

    if (is_boundary(2,1)) then

       ! no dissipation at j=1 -> initialize increment array
       do i=1,nx
           Knutil(i,1,:)=0.0_wp
       enddo

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

       ! no dissipation at j=ny -> initialize increment array
       do i=1,nx
           Knutil(i,ny,:)=0.0_wp
       enddo

    endif

    ! Compute filtering terms
    ! ------------------------
    if (is_sw_edoh) then
       if (is_curv3) then
          do k=1,nz
             do j=ndy_d,nfy_d
                do i=1,nx
                   ! contravariant velocity
                   vc=(uu(i,j,k)*eta_x(i,j,k)+vv(i,j,k)*eta_y(i,j,k)+ww(i,j,k)*eta_z(i,j,k))*ijacob3(i,j,k)
                   ! compute local CFL
                   cfll=(abs(vc)+c_(i,j,k)*g3_eta(i,j,k))*deltat
                   ! switch Edoh et al.
                   mul=min(1.0_wp,cfll/chi)!(chi,cfll)
                   ! If call to local_cfl
                   !mul=min(1.0_wp,cfl_j(i,j,k)/chi)

                   ! compute dissipation term: delta_j/delta_y=(D_{j+1/2}-D_{j-1/2})/dy
                   Knutil(i,j,k) = Knutil(i,j,k) +mul*(Dnutil(i,j-1,k) -Dnutil(i,j,k))
                enddo
             enddo
          enddo
       elseif (is_curv) then
          do k=1,nz
             do j=ndy_d,nfy_d
                do i=1,nx
                   ! contravariant velocity
                   vc=(vv(i,j,k)*x_ksi(i,j)-uu(i,j,k)*y_ksi(i,j))*ijacob(i,j)
                   ! compute local CFL
                   cfll=(abs(vc)+c_(i,j,k)*g_eta(i,j))*deltat
                   ! switch Edoh et al.
                   mul=min(1.0_wp,cfll/chi)!(chi,cfll)
                   ! If call to local_cfl
                   !mul=min(1.0_wp,cfl_j(i,j,k)/chi)

                   ! compute dissipation term: delta_j/delta_y=(D_{j+1/2}-D_{j-1/2})/dy
                   Knutil(i,j,k) = Knutil(i,j,k) +mul*(Dnutil(i,j-1,k) -Dnutil(i,j,k))
                enddo
             enddo
          enddo
       else
          do k=1,nz
             do j=ndy_d,nfy_d
                do i=1,nx
                   ! compute local CFL
                   cfll=(abs(vv(i,j,k))+c_(i,j,k))*deltat*idy(j)
                   ! switch Edoh et al.
                   mul=min(1.0_wp,cfll/chi)!(chi,cfll)
                   ! If call to local_cfl
                   !mul=min(1.0_wp,cfl_j(i,j,k)/chi)

                   ! compute dissipation term: delta_j/delta_y=(D_{j+1/2}-D_{j-1/2})/dy
                   Knutil(i,j,k) = Knutil(i,j,k) +mul*(Dnutil(i,j-1,k) -Dnutil(i,j,k))
                enddo
             enddo
          enddo
       endif
    else
       do k=1,nz
          do j=ndy_d,nfy_d
             do i=1,nx
                ! compute dissipation term: delta_j/delta_x=(D_{j+1/2}-D_{j-1/2})/dx
                Knutil(i,j,k) = Knutil(i,j,k) +Dnutil(i,j-1,k) -Dnutil(i,j,k)
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
             do k=ndz_d1,nfz_d1
!                 psens(i,j,k)=sens(i,j,k)*abs(prs(i,j,k-1)-2.0_wp*prs(i,j,k)+prs(i,j,k+1)) / &
!                                          abs(prs(i,j,k-1)+2.0_wp*prs(i,j,k)+prs(i,j,k+1))

                ! an alternative sensor (implemented by OY)
                ! a mix of a TVD concept and the traditional one (suggested by Swanson et al. [1998])
                psens_trad = prs(i,j,k-1)+2.0_wp*prs(i,j,k)+prs(i,j,k+1)
                psens_tvd = abs(prs(i,j,k+1)-prs(i,j,k))+abs(prs(i,j,k)-prs(i,j,k-1))
                psens(i,j,k)=sens(i,j,k)*abs(prs(i,j,k-1)-2.0_wp*prs(i,j,k)+prs(i,j,k+1)) / &
                                      (csens*psens_trad + (1.0_wp-csens)*psens_tvd)
             enddo
          enddo
       enddo

       ! Filtering in k-direction [compute Dnutil at k+1/2]
       ! ========================
       if (is_boundary(3,1)) then

          ! no dissipation at k=1 -> initialize increment array
          Knutil(:,:,1)=0.0_wp

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

       ! Interior points: 3rd-order centered dissipation
       ! -----------------------------------------------
       do k=ndz_di_r,nfz_di_r
          do j=1,ny
             do i=1,nx
                xk0=xnu0*max(psens(i,j,k+1),psens(i,j,k))
                xk4=max(0.0_wp,xnu4+xk0)

                Dnutil(i,j,k) = (nutil_n(i,j,k+1)-nutil_n(i,j,k))*xk0      &
                            + (d4(1)*(nutil_n(i,j,k+1)-nutil_n(i,j,k  )) &
                            +  d4(2)*(nutil_n(i,j,k+2)-nutil_n(i,j,k-1)))*xk4
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

          ! no dissipation at k=nk -> initialize increment array
          Knutil(:,:,nz)=0.0_wp
       
       endif
       
       ! Compute filtering terms
       ! ------------------------
       if (is_sw_edoh) then
          if (is_curv3) then
             do k=ndz_d,nfz_d
                do j=1,ny
                   do i=1,nx
                      ! contravariant velocity
                      vc=(uu(i,j,k)*phi_x(i,j,k)+vv(i,j,k)*phi_y(i,j,k)+ww(i,j,k)*phi_z(i,j,k))*ijacob3(i,j,k)
                      ! compute local CFL
                      cfll=(abs(vc)+c_(i,j,k)*g3_phi(i,j,k))*deltat
                      ! switch Edoh et al.
                      mul=min(1.0_wp,cfll/chi)!(chi,cfll)
                      ! If call to local_cfl
                      !mul=min(1.0_wp,cfl_k(i,j,k)/chi)
                      
                      ! compute dissipation term: delta_k/delta_z=(D_{k+1/2}-D_{k-1/2})/dz
                      Knutil(i,j,k) = Knutil(i,j,k) +mul*(Dnutil(i,j,k-1) -Dnutil(i,j,k))
                   enddo
                enddo
             enddo
          else
             do k=ndz_d,nfz_d
                do j=1,ny
                   do i=1,nx
                      ! compute local CFL
                      cfll=(abs(ww(i,j,k))+c_(i,j,k))*deltat*idz(k)
                      ! switch Edoh et al.
                      mul=min(1.0_wp,cfll/chi)!(chi,cfll)
                      ! If call to local_cfl
                      !mul=min(1.0_wp,cfl_k(i,j,k)/chi)

                      ! compute dissipation term: delta_k/delta_z=(D_{k+1/2}-D_{k-1/2})/dz
                      Knutil(i,j,k) = Knutil(i,j,k) +mul*(Dnutil(i,j,k-1) -Dnutil(i,j,k))
                   enddo
                enddo
             enddo
          endif
       else
          do k=ndz_d,nfz_d
             do j=1,ny
                do i=1,nx
                   ! compute dissipation term: delta_k/delta_z=(D_{k+1/2}-D_{k-1/2})/dz
                   Knutil(i,j,k) = Knutil(i,j,k) +Dnutil(i,j,k-1) -Dnutil(i,j,k)
                enddo
             enddo
          enddo
       endif

    endif

  end subroutine filter_o4_shock_SA

  !==============================================================================
  subroutine filter_o2_shock_SA
  !==============================================================================
    !> Apply 2nd-order filtering with Jameson-type shock capturing
    !> (Conservative version inspired from DNC order 9)
    !> modified version XG 2023
  !==============================================================================
    use mod_time
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: vc,cfll,mul
    real(wp), dimension(0:nx,0:ny,0:nz) :: Dnutil
    ! ---------------------------------------------------------------------------
    real(wp) :: div2,om1,om2,om3
    real(wp) :: xk0,xk2,xk4
    real(wp) :: psens_trad,psens_tvd
    real(wp), dimension(ndx_d1:nfx_d1,ndy_d1:nfy_d1,ndz_d1:nfz_d1) :: sens,psens
    ! ---------------------------------------------------------------------------

    ! Compute Ducros sensor
    ! =====================
    if (is_ducros) then
       do k=ndz_d1,nfz_d1
          do j=ndy_d1,nfy_d1
             do i=ndx_d1,nfx_d1
                div2=(dux(i,j,k)+dvy(i,j,k)+dwz(i,j,k))**2
                om1 =(dwy(i,j,k)-dvz(i,j,k))**2
                om2 =(duz(i,j,k)-dwx(i,j,k))**2
                om3 =(dvx(i,j,k)-duy(i,j,k))**2
                ! Ducros' sensor
                sens(i,j,k)=div2/(div2+om1+om2+om3+cutoff)
                !uvar(i,j,k,1)=sens(i,j,k)
             enddo
          enddo
       enddo
    else
       sens=1.0_wp
    endif
    
    ! Compute pressure sensor along i-direction
    ! =========================================
    do k=1,nz
       do j=1,ny
          do i=ndx_d1,nfx_d1
!              psens(i,j,k)=sens(i,j,k)*abs(prs(i-1,j,k)-2.0_wp*prs(i,j,k)+prs(i+1,j,k)) / &
!                                       abs(prs(i-1,j,k)+2.0_wp*prs(i,j,k)+prs(i+1,j,k))

             ! an alternative sensor (implemented by OY)
             ! a mix of a TVD concept and the traditional one (suggested by Swanson et al. [1998])
             psens_trad = prs(i-1,j,k)+2.0_wp*prs(i,j,k)+prs(i+1,j,k)
             psens_tvd = abs(prs(i+1,j,k)-prs(i,j,k))+abs(prs(i,j,k)-prs(i-1,j,k))
             psens(i,j,k)=sens(i,j,k)*abs(prs(i-1,j,k)-2.0_wp*prs(i,j,k)+prs(i+1,j,k)) / &
                                      (csens*psens_trad + (1.0_wp-csens)*psens_tvd)

          enddo
       enddo
    enddo

    ! Filtering in i-direction [compute Dnutil at i+1/2]
    ! ========================

    if (is_boundary(1,1)) then
       ! no dissipation at i=1 -> initialize increment array
       do j=1,ny
           Knutil(1,j,:)=0.0_wp
       enddo
    endif

    ! 1st-order centered dissipation
    ! ------------------------------
    do k=1,nz
       do j=1,ny
          do i=ndx_di_r,nfx_di_r
             xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))
             xk2=max(0.0_wp,xnu2+xk0)
             !uvar(i,j,k,1)=xk2

             Dnutil(i,j,k) = (nutil_n(i+1,j,k)-nutil_n(i,j,k))*xk0 &
                         + d2(1)*(nutil_n(i+1,j,k)- nutil_n(i,j,k))*xk2
          enddo
       enddo
    enddo

    ! no dissipation at i=nx -> initialize increment array
    if (is_boundary(1,2)) then
       do j=1,ny
           Knutil(nx,j,:)=0.0_wp
       enddo
    endif

    ! Compute filtering terms
    ! ------------------------
    if (is_sw_edoh) then
       if (is_curv3) then
          do k=1,nz
             do j=1,ny
                do i=ndx_d_r,nfx_d_r
                   ! contravariant velocity
                   vc=(uu(i,j,k)*ksi_x(i,j,k)+vv(i,j,k)*ksi_y(i,j,k)+ww(i,j,k)*ksi_z(i,j,k))*ijacob3(i,j,k)
                   ! compute local CFL
                   cfll=(abs(vc)+c_(i,j,k)*g3_ksi(i,j,k))*deltat
                   ! switch Edoh et al.
                   mul=min(1.0_wp,cfll/chi)
                   ! If call to local_cfl
                   !mul=min(1.0_wp,cfl_j(i,j,k)/chi)

                   ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dx
                   Knutil(i,j,k) = mul*(Dnutil(i-1,j,k) -Dnutil(i,j,k))
                enddo
             enddo
          enddo
       elseif (is_curv) then
          do k=1,nz
             do j=1,ny
                do i=ndx_d_r,nfx_d_r
                   ! contravariant velocity
                   vc=(uu(i,j,k)*y_eta(i,j)-vv(i,j,k)*x_eta(i,j))*ijacob(i,j)
                   ! compute local CFL
                   cfll=(abs(vc)+c_(i,j,k)*g_ksi(i,j))*deltat
                   ! switch Edoh et al.
                   mul=min(1.0_wp,cfll/chi)
                   ! If call to local_cfl
                   !mul=min(1.0_wp,cfl_j(i,j,k)/chi)

                   ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dx
                   Knutil(i,j,k) = mul*(Dnutil(i-1,j,k) -Dnutil(i,j,k))
                enddo
             enddo
          enddo
       else
          do k=1,nz
             do j=1,ny
                do i=ndx_d_r,nfx_d_r
                   ! compute local CFL
                   cfll=(abs(uu(i,j,k))+c_(i,j,k))*deltat*idx(i)
                   ! switch Edoh et al.
                   mul=min(1.0_wp,cfll/chi)!(chi,cfll)
                   
                   ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dx
                   Knutil(i,j,k) = mul*(Dnutil(i-1,j,k) -Dnutil(i,j,k))
                enddo
             enddo
          enddo
       endif
    else
       do k=1,nz
          do j=1,ny
             do i=ndx_d_r,nfx_d_r
                ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dx
                Knutil(i,j,k) = Dnutil(i-1,j,k) -Dnutil(i,j,k)
             enddo
          enddo
       enddo
    endif

    ! Compute pressure sensor along j-direction
    ! =========================================
    do k=1,nz
       do j=ndy_d1,nfy_d1
          do i=1,nx
!              psens(i,j,k)=sens(i,j,k)*abs(prs(i,j-1,k)-2.0_wp*prs(i,j,k)+prs(i,j+1,k)) / &
!                                       abs(prs(i,j-1,k)+2.0_wp*prs(i,j,k)+prs(i,j+1,k))

             ! an alternative sensor (implemented by OY)
             ! a mix of a TVD concept and the traditional one (suggested by Swanson et al. [1998])
             psens_trad = prs(i,j-1,k)+2.0_wp*prs(i,j,k)+prs(i,j+1,k)
             psens_tvd = abs(prs(i,j+1,k)-prs(i,j,k))+abs(prs(i,j,k)-prs(i,j-1,k))
             psens(i,j,k)=sens(i,j,k)*abs(prs(i,j-1,k)-2.0_wp*prs(i,j,k)+prs(i,j+1,k)) / &
                                      (csens*psens_trad + (1.0_wp-csens)*psens_tvd)
          enddo
       enddo
    enddo

    ! Filtering in y-direction [compute Dnutil at j+1/2]
    ! ========================

    if (is_boundary(2,1)) then
       ! no dissipation at j=1 -> initialize increment array
       do i=1,nx
           Knutil(i,1,:)=0.0_wp
       enddo
    endif

    ! 1st-order centered dissipation
    ! ------------------------------
    do k=1,nz
       do j=ndy_di_r,nfy_di_r
          do i=1,nx
             xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))
             xk2=max(0.0_wp,xnu2+xk0)
             !uvar(i,j,k,2)=xk2

             Dnutil(i,j,k) = (nutil_n(i,j+1,k)-nutil_n(i,j,k))*xk0 &
                         + d2(1)*(nutil_n(i,j+1,k)-nutil_n(i,j,k))*xk2
          enddo
       enddo
    enddo

    if (is_boundary(2,2)) then
       ! no dissipation at j=ny -> initialize increment array
       do i=1,nx
           Knutil(i,ny,:)=0.0_wp
       enddo
    endif

    ! Compute filtering terms
    ! ------------------------
    if (is_sw_edoh) then
       if (is_curv3) then
          do k=1,nz
             do j=ndy_d_r,nfy_d_r
                do i=1,nx
                   ! contravariant velocity
                   vc=(uu(i,j,k)*eta_x(i,j,k)+vv(i,j,k)*eta_y(i,j,k)+ww(i,j,k)*eta_z(i,j,k))*ijacob3(i,j,k)
                   ! compute local CFL
                   cfll=(abs(vc)+c_(i,j,k)*g3_eta(i,j,k))*deltat
                   ! switch Edoh et al.
                   mul=min(1.0_wp,cfll/chi)!(chi,cfll)
                   ! If call to local_cfl
                   !mul=min(1.0_wp,cfl_j(i,j,k)/chi)

                   ! compute dissipation term: delta_j/delta_y=(D_{j+1/2}-D_{j-1/2})/dy
                   Knutil(i,j,k) = Knutil(i,j,k) +mul*(Dnutil(i,j-1,k) -Dnutil(i,j,k))
                enddo
             enddo
          enddo
       elseif (is_curv) then
          do k=1,nz
             do j=ndy_d_r,nfy_d_r
                do i=1,nx
                   ! contravariant velocity
                   vc=(vv(i,j,k)*x_ksi(i,j)-uu(i,j,k)*y_ksi(i,j))*ijacob(i,j)
                   ! compute local CFL
                   cfll=(abs(vc)+c_(i,j,k)*g_eta(i,j))*deltat
                   ! switch Edoh et al.
                   mul=min(1.0_wp,cfll/chi)!(chi,cfll)
                   ! If call to local_cfl
                   !mul=min(1.0_wp,cfl_j(i,j,k)/chi)

                   ! compute dissipation term: delta_j/delta_y=(D_{j+1/2}-D_{j-1/2})/dy
                   Knutil(i,j,k) = Knutil(i,j,k) +mul*(Dnutil(i,j-1,k) -Dnutil(i,j,k))
                enddo
             enddo
          enddo
       else
          do k=1,nz
             do j=ndy_d_r,nfy_d_r
                do i=1,nx
                   ! compute local CFL
                   cfll=(abs(vv(i,j,k))+c_(i,j,k))*deltat*idy(j)
                   ! switch Edoh et al.
                   mul=min(1.0_wp,cfll/chi)!(chi,cfll)
                   ! If call to local_cfl
                   !mul=min(1.0_wp,cfl_j(i,j,k)/chi)

                   ! compute dissipation term: delta_j/delta_y=(D_{j+1/2}-D_{j-1/2})/dy
                   Knutil(i,j,k) = Knutil(i,j,k) +mul*(Dnutil(i,j-1,k) -Dnutil(i,j,k))
                enddo
             enddo
          enddo
       endif
    else
       do k=1,nz
          do j=ndy_d_r,nfy_d_r
             do i=1,nx
                ! compute dissipation term: delta_j/delta_x=(D_{j+1/2}-D_{j-1/2})/dx
                Knutil(i,j,k) = Knutil(i,j,k) +Dnutil(i,j-1,k) -Dnutil(i,j,k)
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
             do k=ndz_d1,nfz_d1
!                 psens(i,j,k)=sens(i,j,k)*abs(prs(i,j,k-1)-2.0_wp*prs(i,j,k)+prs(i,j,k+1)) / &
!                                          abs(prs(i,j,k-1)+2.0_wp*prs(i,j,k)+prs(i,j,k+1))

                ! an alternative sensor (implemented by OY)
                ! a mix of a TVD concept and the traditional one (suggested by Swanson et al. [1998])
                psens_trad = prs(i,j,k-1)+2.0_wp*prs(i,j,k)+prs(i,j,k+1)
                psens_tvd = abs(prs(i,j,k+1)-prs(i,j,k))+abs(prs(i,j,k)-prs(i,j,k-1))
                psens(i,j,k)=sens(i,j,k)*abs(prs(i,j,k-1)-2.0_wp*prs(i,j,k)+prs(i,j,k+1)) / &
                                      (csens*psens_trad + (1.0_wp-csens)*psens_tvd)
             enddo
          enddo
       enddo

       ! Filtering in k-direction [compute Dnutil at k+1/2]
       ! ========================
       if (is_boundary(3,1)) then
          ! no dissipation at k=1 -> initialize increment array
          Knutil(:,:,1)=0.0_wp
       endif

       ! Interior points: 1rd-order centered dissipation
       ! -----------------------------------------------
       do k=ndz_di_r,nfz_di_r
          do j=1,ny
             do i=1,nx
                xk0=xnu0*max(psens(i,j,k+1),psens(i,j,k))
                xk4=max(0.0_wp,xnu4+xk0)

                Dnutil(i,j,k) = (nutil_n(i,j,k+1)-nutil_n(i,j,k))*xk0      &
                            + d2(1)*(nutil_n(i,j,k+1)-nutil_n(i,j,k  ))
             enddo
          enddo
       enddo

       if (is_boundary(3,2)) then
          ! no dissipation at k=nz -> initialize increment array
          Knutil(:,:,nz)=0.0_wp
       endif
       
       ! Compute filtering terms
       ! ------------------------
       if (is_sw_edoh) then
          if (is_curv3) then
             do k=ndz_d_r,nfz_d_r
                do j=1,ny
                   do i=1,nx
                      ! contravariant velocity
                      vc=(uu(i,j,k)*phi_x(i,j,k)+vv(i,j,k)*phi_y(i,j,k)+ww(i,j,k)*phi_z(i,j,k))*ijacob3(i,j,k)
                      ! compute local CFL
                      cfll=(abs(vc)+c_(i,j,k)*g3_phi(i,j,k))*deltat
                      ! switch Edoh et al.
                      mul=min(1.0_wp,cfll/chi)!(chi,cfll)
                      ! If call to local_cfl
                      !mul=min(1.0_wp,cfl_k(i,j,k)/chi)
                      
                      ! compute dissipation term: delta_k/delta_z=(D_{k+1/2}-D_{k-1/2})/dz
                      Knutil(i,j,k) = Knutil(i,j,k) +mul*(Dnutil(i,j,k-1) -Dnutil(i,j,k))
                   enddo
                enddo
             enddo
          else
             do k=ndz_d_r,nfz_d_r
                do j=1,ny
                   do i=1,nx
                      ! compute local CFL
                      cfll=(abs(ww(i,j,k))+c_(i,j,k))*deltat*idz(k)
                      ! switch Edoh et al.
                      mul=min(1.0_wp,cfll/chi)!(chi,cfll)
                      ! If call to local_cfl
                      !mul=min(1.0_wp,cfl_k(i,j,k)/chi)

                      ! compute dissipation term: delta_k/delta_z=(D_{k+1/2}-D_{k-1/2})/dz
                      Knutil(i,j,k) = Knutil(i,j,k) +mul*(Dnutil(i,j,k-1) -Dnutil(i,j,k))
                   enddo
                enddo
             enddo
          endif
       else
          do k=ndz_d_r,nfz_d_r
             do j=1,ny
                do i=1,nx
                   ! compute dissipation term: delta_k/delta_z=(D_{k+1/2}-D_{k-1/2})/dz
                   Knutil(i,j,k) = Knutil(i,j,k) +Dnutil(i,j,k-1) -Dnutil(i,j,k)
                enddo
             enddo
          enddo
       endif

    endif

  end subroutine filter_o2_shock_SA


end module mod_filtering_shock_rans
