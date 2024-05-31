!==============================================================================
module mod_filtering_shock
!==============================================================================
  !> Module for filtering with Jameson-type shock capturing
!==============================================================================
  use mod_flow
  use mod_artvisc
  use mod_constant
  implicit none
  ! ---------------------------------------------------------------------------
  real(wp) :: chi
  ! declaration of coefficients in mod_artvisc.f90
  ! real(wp), parameter, private :: cutoff = 1.e-16_wp
  ! amplitude dissipation
  real(wp), private :: xnu0 ! low-order smoothing term
  real(wp), private :: xnu10,xnu8,xnu6,xnu4,xnu2
  ! ---------------------------------------------------------------------------

contains

  !============================================================================
  subroutine init_filter_shock(dissip_coeff,dissip_shock)
  !============================================================================
    !> Initialization of filtering with Jameson-type shock capturing
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    real(wp), intent(in) :: dissip_coeff,dissip_shock
    ! -------------------------------------------------------------------------
    integer, parameter :: iSF=0

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

    chi=dissip_coeff
    
    ! Filtering coefficients (DNC-Jameson)
    ! =================================
    ! /!\ Change sign because added to increment
    ! order 1
    d2(1)=-1.0_wp

    ! order 3
    d4(1)=-3.0_wp
    d4(2)= 1.0_wp

    ! order 5
    d6(1)=-10.0_wp
    d6(2)= 5.0_wp
    d6(3)=-1.0_wp

    ! order 7
    d8(1)=-35.0_wp
    d8(2)= 21.0_wp
    d8(3)= -7.0_wp
    d8(4)=  1.0_wp

    ! order 9
    d10(1)=-126.0_wp
    d10(2)= 84.0_wp
    d10(3)=-36.0_wp
    d10(4)= 9.0_wp
    d10(5)=-1.0_wp
    
    ! Dissipation coefficients [from param.ini]
    ! ========================
    xnu10=dissip_coeff/1260.0_wp
    xnu8= dissip_coeff/280.0_wp
    xnu6= dissip_coeff/60.0_wp
    xnu4= dissip_coeff/12.0_wp
    xnu2= 0.0_wp*dissip_coeff/2.0_wp

    if (iSF==1) then
       d4(1)=-1.0_wp
       d4(2)= 0.0_wp

       ! order 5
       d6(1)=-1.0_wp
       d6(2)= 0.0_wp
       d6(3)= 0.0_wp

       ! order 7
       d8(1)=-1.0_wp
       d8(2)= 0.0_wp
       d8(3)= 0.0_wp
       d8(4)= 0.0_wp

       ! order 9
       d10(1)=-1.0_wp
       d10(2)= 0.0_wp
       d10(3)=-0.0_wp
       d10(4)= 0.0_wp
       d10(5)=-0.0_wp

       ! Dissipation coefficients [from param.ini]
       ! ========================
       xnu10=dissip_coeff/2.0_wp
       xnu8= dissip_coeff/2.0_wp
       xnu6= dissip_coeff/2.0_wp
       xnu4= dissip_coeff/2.0_wp

    else if (iSF==2) then
       ! 7 points => 5 points
       d6(1)=-3.0_wp
       d6(2)= 1.0_wp
       d6(3)= 0.0_wp

       ! 9 points => 5 points
       d8(1)=-3.0_wp
       d8(2)= 1.0_wp
       d8(3)= 0.0_wp
       d8(4)= 0.0_wp

       ! 11 points => 5 points
       d10(1)=-3.0_wp
       d10(2)= 1.0_wp
       d10(3)=-0.0_wp
       d10(4)= 0.0_wp
       d10(5)=-0.0_wp

       ! Dissipation coefficients [from param.ini]
       ! ========================
       xnu10=dissip_coeff/12.0_wp
       xnu8= dissip_coeff/12.0_wp
       xnu6= dissip_coeff/12.0_wp
       xnu4= dissip_coeff/12.0_wp

    else if (iSF==3) then
       ! 9 points => 7 points
       d8(1)=-10.0_wp
       d8(2)= 5.0_wp
       d8(3)=-1.0_wp
       d8(4)= 0.0_wp

       ! 11 points => 7 points
       d10(1)=-10.0_wp
       d10(2)= 5.0_wp
       d10(3)=-1.0_wp
       d10(4)= 0.0_wp
       d10(5)=-0.0_wp

       ! Dissipation coefficients [from param.ini]
       ! ========================
       xnu10=dissip_coeff/60.0_wp
       xnu8= dissip_coeff/60.0_wp

    else if (iSF==4) then
       ! 11 points => 9 points
       d10(1)=-35.0_wp
       d10(2)= 21.0_wp
       d10(3)=-7.0_wp
       d10(4)= 1.0_wp
       d10(5)=-0.0_wp

       ! Dissipation coefficients [from param.ini]
       ! ========================
       xnu10=dissip_coeff/280.0_wp
    endif

    ! low-order dissipation term [from param.ini]
    ! ==========================
    ! /!\ Change sign because added to increment
    xnu0=-dissip_shock
    
  end subroutine init_filter_shock

  !==============================================================================
  subroutine filter_o10_shock
  !==============================================================================
    !> Apply 10th-order filtering with Jameson-type shock capturing
    !> (Conservative version inspired from DNC order 9)
    !> modified version XG 2023
  !==============================================================================
    use mod_time
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: vc,cfll,mul
    real(wp), dimension(0:nx,0:ny,0:nz) :: Drho,Drhou,Drhov,Drhow,Drhoe
    ! ---------------------------------------------------------------------------
    real(wp) :: div2,om1,om2,om3
    real(wp) :: xk0,xk2,xk4,xk6,xk8,xk10
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

             !uvar(i,j,k,2)=psens(i,j,k)
          enddo
       enddo
    enddo
    
    ! Filtering in i-direction [compute Drho at i+1/2]
    ! ========================
    
    if (is_boundary(1,1)) then
       ! no dissipation at i=1 -> initialize increment array
       do j=1,ny
           Krho(1,j,:)=0.0_wp
          Krhou(1,j,:)=0.0_wp
          Krhov(1,j,:)=0.0_wp
          Krhow(1,j,:)=0.0_wp
          Krhoe(1,j,:)=0.0_wp
       enddo
       
       ! 1st-order centered dissipation
       ! ------------------------------
       if (is_bc_wall(1,1)) then
          i=1
          do k=1,nz
             do j=1,ny
                xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)
                !uvar(i,j,k,1)=xk2

                Drho(i,j,k) = (rho_n(i+1,j,k)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i+1,j,k)- rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
       else
          i=1
          do k=1,nz
             do j=1,ny
                xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)
                !uvar(i,j,k,1)=xk2

                Drho(i,j,k) = (rho_n(i+1,j,k)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i+1,j,k)- rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
       endif
       
       ! 3rd-order centered dissipation
       ! ------------------------------
       i=2
       do k=1,nz
          do j=1,ny            
             xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))
             xk4=max(0.0_wp,xnu4+xk0)
             !uvar(i,j,k,1)=xk4/xnu4
             
             Drho(i,j,k) = (rho_n(i+1,j,k)-rho_n(i,j,k))*xk0      &
                         + (d4(1)*(rho_n(i+1,j,k)-rho_n(i  ,j,k)) &
                           +d4(2)*(rho_n(i+2,j,k)-rho_n(i-1,j,k)))*xk4
             Drhou(i,j,k)= (rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk0      &
                         + (d4(1)*(rhou_n(i+1,j,k)-rhou_n(i  ,j,k)) &
                           +d4(2)*(rhou_n(i+2,j,k)-rhou_n(i-1,j,k)))*xk4
             Drhov(i,j,k)= (rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk0      &
                         + (d4(1)*(rhov_n(i+1,j,k)-rhov_n(i  ,j,k)) &
                           +d4(2)*(rhov_n(i+2,j,k)-rhov_n(i-1,j,k)))*xk4
             Drhow(i,j,k)= (rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk0      &
                         + (d4(1)*(rhow_n(i+1,j,k)-rhow_n(i  ,j,k)) &
                           +d4(2)*(rhow_n(i+2,j,k)-rhow_n(i-1,j,k)))*xk4
             Drhoe(i,j,k)= (rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk0      &
                         + (d4(1)*(rhoe_n(i+1,j,k)-rhoe_n(i  ,j,k)) &
                           +d4(2)*(rhoe_n(i+2,j,k)-rhoe_n(i-1,j,k)))*xk4
          enddo
       enddo
       
       ! 5th-order centered dissipation
       ! ------------------------------
       i=3
       do k=1,nz
          do j=1,ny             
             xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))
             !xk6=max(0.0_wp,xnu6+xk0/5.0_wp)
             xk6=max(0.0_wp,xnu6+xk0)
             !uvar(i,j,k,1)=xk6/xnu6
             
             Drho(i,j,k) = (rho_n(i+1,j,k)-rho_n(i,j,k))*xk0      &
                         + (d6(1)*(rho_n(i+1,j,k)-rho_n(i  ,j,k)) &
                           +d6(2)*(rho_n(i+2,j,k)-rho_n(i-1,j,k)) &
                           +d6(3)*(rho_n(i+3,j,k)-rho_n(i-2,j,k)))*xk6
             Drhou(i,j,k)= (rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk0      &
                         + (d6(1)*(rhou_n(i+1,j,k)-rhou_n(i  ,j,k)) &
                           +d6(2)*(rhou_n(i+2,j,k)-rhou_n(i-1,j,k)) &
                           +d6(3)*(rhou_n(i+3,j,k)-rhou_n(i-2,j,k)))*xk6
             Drhov(i,j,k)= (rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk0      &
                         + (d6(1)*(rhov_n(i+1,j,k)-rhov_n(i  ,j,k)) &
                           +d6(2)*(rhov_n(i+2,j,k)-rhov_n(i-1,j,k)) &
                           +d6(3)*(rhov_n(i+3,j,k)-rhov_n(i-2,j,k)))*xk6
             Drhow(i,j,k)= (rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk0      &
                         + (d6(1)*(rhow_n(i+1,j,k)-rhow_n(i  ,j,k)) &
                           +d6(2)*(rhow_n(i+2,j,k)-rhow_n(i-1,j,k)) &
                           +d6(3)*(rhow_n(i+3,j,k)-rhow_n(i-2,j,k)))*xk6
             Drhoe(i,j,k)= (rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk0      &
                         + (d6(1)*(rhoe_n(i+1,j,k)-rhoe_n(i  ,j,k)) &
                           +d6(2)*(rhoe_n(i+2,j,k)-rhoe_n(i-1,j,k)) &
                           +d6(3)*(rhoe_n(i+3,j,k)-rhoe_n(i-2,j,k)))*xk6
          enddo
       enddo
       
       ! 7th-order centered dissipation
       ! ------------------------------
       i=4
       do k=1,nz
          do j=1,ny             
             xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))
             !xk8=max(0.0_wp,xnu8+3.0_wp/70.0_wp*xk0)
             xk8=max(0.0_wp,xnu8+xk0)
             !uvar(i,j,k,1)=xk8/xnu8
             
             Drho(i,j,k) = (rho_n(i+1,j,k)-rho_n(i,j,k))*xk0      &
                         + (d8(1)*(rho_n(i+1,j,k)-rho_n(i  ,j,k)) &
                           +d8(2)*(rho_n(i+2,j,k)-rho_n(i-1,j,k)) &
                           +d8(3)*(rho_n(i+3,j,k)-rho_n(i-2,j,k)) &
                           +d8(4)*(rho_n(i+4,j,k)-rho_n(i-3,j,k)))*xk8
             Drhou(i,j,k)= (rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk0      &
                         + (d8(1)*(rhou_n(i+1,j,k)-rhou_n(i  ,j,k)) &
                           +d8(2)*(rhou_n(i+2,j,k)-rhou_n(i-1,j,k)) &
                           +d8(3)*(rhou_n(i+3,j,k)-rhou_n(i-2,j,k)) &
                           +d8(4)*(rhou_n(i+4,j,k)-rhou_n(i-3,j,k)))*xk8
             Drhov(i,j,k)= (rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk0      &
                         + (d8(1)*(rhov_n(i+1,j,k)-rhov_n(i  ,j,k)) &
                           +d8(2)*(rhov_n(i+2,j,k)-rhov_n(i-1,j,k)) &
                           +d8(3)*(rhov_n(i+3,j,k)-rhov_n(i-2,j,k)) &
                           +d8(4)*(rhov_n(i+4,j,k)-rhov_n(i-3,j,k)))*xk8
             Drhow(i,j,k)= (rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk0      &
                         + (d8(1)*(rhow_n(i+1,j,k)-rhow_n(i  ,j,k)) &
                           +d8(2)*(rhow_n(i+2,j,k)-rhow_n(i-1,j,k)) &
                           +d8(3)*(rhow_n(i+3,j,k)-rhow_n(i-2,j,k)) &
                           +d8(4)*(rhow_n(i+4,j,k)-rhow_n(i-3,j,k)))*xk8
             Drhoe(i,j,k)= (rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk0      &
                         + (d8(1)*(rhoe_n(i+1,j,k)-rhoe_n(i  ,j,k)) &
                           +d8(2)*(rhoe_n(i+2,j,k)-rhoe_n(i-1,j,k)) &
                           +d8(3)*(rhoe_n(i+3,j,k)-rhoe_n(i-2,j,k)) &
                           +d8(4)*(rhoe_n(i+4,j,k)-rhoe_n(i-3,j,k)))*xk8
          enddo
       enddo
    endif

    ! Interior points: 9th-order centered dissipation
    ! -----------------------------------------------
    do k=1,nz
       do j=1,ny
          do i=ndx_di,nfx_di             
             xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))
             xk10=max(0.0_wp,xnu10+xk0/105.0_wp)
             !if (i>=1) uvar(i,j,k,1)=xk10/xnu10
             
             Drho(i,j,k) = (rho_n(i+1,j,k)-rho_n(i,j,k))*xk0       &
                         + (d10(1)*(rho_n(i+1,j,k)-rho_n(i  ,j,k)) &
                           +d10(2)*(rho_n(i+2,j,k)-rho_n(i-1,j,k)) &
                           +d10(3)*(rho_n(i+3,j,k)-rho_n(i-2,j,k)) &
                           +d10(4)*(rho_n(i+4,j,k)-rho_n(i-3,j,k)) &
                           +d10(5)*(rho_n(i+5,j,k)-rho_n(i-4,j,k)))*xk10
             
             Drhou(i,j,k)= (rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk0       &
                         + (d10(1)*(rhou_n(i+1,j,k)-rhou_n(i  ,j,k)) &
                           +d10(2)*(rhou_n(i+2,j,k)-rhou_n(i-1,j,k)) &
                           +d10(3)*(rhou_n(i+3,j,k)-rhou_n(i-2,j,k)) &
                           +d10(4)*(rhou_n(i+4,j,k)-rhou_n(i-3,j,k)) &
                           +d10(5)*(rhou_n(i+5,j,k)-rhou_n(i-4,j,k)))*xk10

             Drhov(i,j,k)= (rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk0       &
                         + (d10(1)*(rhov_n(i+1,j,k)-rhov_n(i  ,j,k)) &
                           +d10(2)*(rhov_n(i+2,j,k)-rhov_n(i-1,j,k)) &
                           +d10(3)*(rhov_n(i+3,j,k)-rhov_n(i-2,j,k)) &
                           +d10(4)*(rhov_n(i+4,j,k)-rhov_n(i-3,j,k)) &
                           +d10(5)*(rhov_n(i+5,j,k)-rhov_n(i-4,j,k)))*xk10

             Drhow(i,j,k)= (rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk0       &
                         + (d10(1)*(rhow_n(i+1,j,k)-rhow_n(i  ,j,k)) &
                           +d10(2)*(rhow_n(i+2,j,k)-rhow_n(i-1,j,k)) &
                           +d10(3)*(rhow_n(i+3,j,k)-rhow_n(i-2,j,k)) &
                           +d10(4)*(rhow_n(i+4,j,k)-rhow_n(i-3,j,k)) &
                           +d10(5)*(rhow_n(i+5,j,k)-rhow_n(i-4,j,k)))*xk10

             Drhoe(i,j,k)= (rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk0       &
                         + (d10(1)*(rhoe_n(i+1,j,k)-rhoe_n(i  ,j,k)) &
                           +d10(2)*(rhoe_n(i+2,j,k)-rhoe_n(i-1,j,k)) &
                           +d10(3)*(rhoe_n(i+3,j,k)-rhoe_n(i-2,j,k)) &
                           +d10(4)*(rhoe_n(i+4,j,k)-rhoe_n(i-3,j,k)) &
                           +d10(5)*(rhoe_n(i+5,j,k)-rhoe_n(i-4,j,k)))*xk10
          enddo
       enddo
    enddo

    if (is_boundary(1,2)) then
       ! 7th-order centered dissipation
       ! ------------------------------
       i=nx-4
       do k=1,nz
          do j=1,ny             
             xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))
             !xk8=max(0.0_wp,xnu8+3.0_wp/70.0_wp*xk0)
             xk8=max(0.0_wp,xnu8+xk0)
             !uvar(i,j,k,1)=xk8/xnu8
             
             Drho(i,j,k) = (rho_n(i+1,j,k)-rho_n(i,j,k))*xk0      &
                         + (d8(1)*(rho_n(i+1,j,k)-rho_n(i  ,j,k)) &
                           +d8(2)*(rho_n(i+2,j,k)-rho_n(i-1,j,k)) &
                           +d8(3)*(rho_n(i+3,j,k)-rho_n(i-2,j,k)) &
                           +d8(4)*(rho_n(i+4,j,k)-rho_n(i-3,j,k)))*xk8
             Drhou(i,j,k)= (rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk0      &
                         + (d8(1)*(rhou_n(i+1,j,k)-rhou_n(i  ,j,k)) &
                           +d8(2)*(rhou_n(i+2,j,k)-rhou_n(i-1,j,k)) &
                           +d8(3)*(rhou_n(i+3,j,k)-rhou_n(i-2,j,k)) &
                           +d8(4)*(rhou_n(i+4,j,k)-rhou_n(i-3,j,k)))*xk8
             Drhov(i,j,k)= (rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk0      &
                         + (d8(1)*(rhov_n(i+1,j,k)-rhov_n(i  ,j,k)) &
                           +d8(2)*(rhov_n(i+2,j,k)-rhov_n(i-1,j,k)) &
                           +d8(3)*(rhov_n(i+3,j,k)-rhov_n(i-2,j,k)) &
                           +d8(4)*(rhov_n(i+4,j,k)-rhov_n(i-3,j,k)))*xk8
             Drhow(i,j,k)= (rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk0      &
                         + (d8(1)*(rhow_n(i+1,j,k)-rhow_n(i  ,j,k)) &
                           +d8(2)*(rhow_n(i+2,j,k)-rhow_n(i-1,j,k)) &
                           +d8(3)*(rhow_n(i+3,j,k)-rhow_n(i-2,j,k)) &
                           +d8(4)*(rhow_n(i+4,j,k)-rhow_n(i-3,j,k)))*xk8
             Drhoe(i,j,k)= (rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk0      &
                         + (d8(1)*(rhoe_n(i+1,j,k)-rhoe_n(i  ,j,k)) &
                           +d8(2)*(rhoe_n(i+2,j,k)-rhoe_n(i-1,j,k)) &
                           +d8(3)*(rhoe_n(i+3,j,k)-rhoe_n(i-2,j,k)) &
                           +d8(4)*(rhoe_n(i+4,j,k)-rhoe_n(i-3,j,k)))*xk8
          enddo 
       enddo
       
       ! 5th-order centered dissipation
       ! ------------------------------
       i=nx-3
       do k=1,nz
          do j=1,ny            
             xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))
             !xk6=max(0.0_wp,xnu6+xk0/5.0_wp)
             xk6=max(0.0_wp,xnu6+xk0)
             !uvar(i,j,k,1)=xk6/xnu6
             
             Drho(i,j,k) = (rho_n(i+1,j,k)-rho_n(i,j,k))*xk0      &
                         + (d6(1)*(rho_n(i+1,j,k)-rho_n(i  ,j,k)) &
                           +d6(2)*(rho_n(i+2,j,k)-rho_n(i-1,j,k)) &
                           +d6(3)*(rho_n(i+3,j,k)-rho_n(i-2,j,k)))*xk6
             Drhou(i,j,k)= (rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk0      &
                         + (d6(1)*(rhou_n(i+1,j,k)-rhou_n(i  ,j,k)) &
                           +d6(2)*(rhou_n(i+2,j,k)-rhou_n(i-1,j,k)) &
                           +d6(3)*(rhou_n(i+3,j,k)-rhou_n(i-2,j,k)))*xk6
             Drhov(i,j,k)= (rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk0      &
                         + (d6(1)*(rhov_n(i+1,j,k)-rhov_n(i  ,j,k)) &
                           +d6(2)*(rhov_n(i+2,j,k)-rhov_n(i-1,j,k)) &
                           +d6(3)*(rhov_n(i+3,j,k)-rhov_n(i-2,j,k)))*xk6
             Drhow(i,j,k)= (rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk0      &
                         + (d6(1)*(rhow_n(i+1,j,k)-rhow_n(i  ,j,k)) &
                           +d6(2)*(rhow_n(i+2,j,k)-rhow_n(i-1,j,k)) &
                           +d6(3)*(rhow_n(i+3,j,k)-rhow_n(i-2,j,k)))*xk6
             Drhoe(i,j,k)= (rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk0      &
                         + (d6(1)*(rhoe_n(i+1,j,k)-rhoe_n(i  ,j,k)) &
                           +d6(2)*(rhoe_n(i+2,j,k)-rhoe_n(i-1,j,k)) &
                           +d6(3)*(rhoe_n(i+3,j,k)-rhoe_n(i-2,j,k)))*xk6
          enddo
       enddo
       
       ! 3rd-order centered dissipation
       ! ------------------------------
       i=nx-2
       do k=1,nz
          do j=1,ny
             xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))
             xk4=max(0.0_wp,xnu4+xk0)
             !uvar(i,j,k,1)=xk4/xnu4
            
             Drho(i,j,k) = (rho_n(i+1,j,k)-rho_n(i,j,k))*xk0      &
                         + (d4(1)*(rho_n(i+1,j,k)-rho_n(i  ,j,k)) &
                           +d4(2)*(rho_n(i+2,j,k)-rho_n(i-1,j,k)))*xk4
             Drhou(i,j,k)= (rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk0      &
                         + (d4(1)*(rhou_n(i+1,j,k)-rhou_n(i  ,j,k)) &
                           +d4(2)*(rhou_n(i+2,j,k)-rhou_n(i-1,j,k)))*xk4
             Drhov(i,j,k)= (rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk0      &
                         + (d4(1)*(rhov_n(i+1,j,k)-rhov_n(i  ,j,k)) &
                           +d4(2)*(rhov_n(i+2,j,k)-rhov_n(i-1,j,k)))*xk4
             Drhow(i,j,k)= (rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk0      &
                         + (d4(1)*(rhow_n(i+1,j,k)-rhow_n(i  ,j,k)) &
                           +d4(2)*(rhow_n(i+2,j,k)-rhow_n(i-1,j,k)))*xk4
             Drhoe(i,j,k)= (rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk0      &
                         + (d4(1)*(rhoe_n(i+1,j,k)-rhoe_n(i  ,j,k)) &
                           +d4(2)*(rhoe_n(i+2,j,k)-rhoe_n(i-1,j,k)))*xk4

          enddo
       enddo
       
       ! 1st-order centered dissipation
       ! ------------------------------
       if (is_bc_wall(1,2)) then
          i=nx-1
          do k=1,nz
             do j=1,ny
                xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)
                !uvar(i,j,k,1)=xk2

                Drho(i,j,k) = (rho_n(i+1,j,k)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i+1,j,k)-rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
       else
          i=nx-1
          do k=1,nz
             do j=1,ny
                xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)
                !uvar(i,j,k,1)=xk2

                Drho(i,j,k) = (rho_n(i+1,j,k)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i+1,j,k)-rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
       endif
       
       ! no dissipation at i=nx -> initialize increment array
       do j=1,ny
           Krho(nx,j,:)=0.0_wp
          Krhou(nx,j,:)=0.0_wp
          Krhov(nx,j,:)=0.0_wp
          Krhow(nx,j,:)=0.0_wp
          Krhoe(nx,j,:)=0.0_wp
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
                   !mul=min(1.0_wp,cfl_i(i,j,k)/chi)

                   ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dx
                   Krho(i,j,k) = mul*(Drho(i-1,j,k) -Drho(i,j,k))
                   Krhou(i,j,k)= mul*(Drhou(i-1,j,k)-Drhou(i,j,k))
                   Krhov(i,j,k)= mul*(Drhov(i-1,j,k)-Drhov(i,j,k))
                   Krhow(i,j,k)= mul*(Drhow(i-1,j,k)-Drhow(i,j,k))
                   Krhoe(i,j,k)= mul*(Drhoe(i-1,j,k)-Drhoe(i,j,k))
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
                   !mul=min(1.0_wp,cfl_i(i,j,k)/chi)

                   ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dx
                   Krho(i,j,k) = mul*(Drho(i-1,j,k) -Drho(i,j,k))
                   Krhou(i,j,k)= mul*(Drhou(i-1,j,k)-Drhou(i,j,k))
                   Krhov(i,j,k)= mul*(Drhov(i-1,j,k)-Drhov(i,j,k))
                   Krhow(i,j,k)= mul*(Drhow(i-1,j,k)-Drhow(i,j,k))
                   Krhoe(i,j,k)= mul*(Drhoe(i-1,j,k)-Drhoe(i,j,k))
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
                   ! If call to local_cfl
                   !mul=min(1.0_wp,cfl_i(i,j,k)/chi)

                   ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dx
                   Krho(i,j,k) = mul*(Drho(i-1,j,k) -Drho(i,j,k))
                   Krhou(i,j,k)= mul*(Drhou(i-1,j,k)-Drhou(i,j,k))
                   Krhov(i,j,k)= mul*(Drhov(i-1,j,k)-Drhov(i,j,k))
                   Krhow(i,j,k)= mul*(Drhow(i-1,j,k)-Drhow(i,j,k))
                   Krhoe(i,j,k)= mul*(Drhoe(i-1,j,k)-Drhoe(i,j,k))
                enddo
             enddo
          enddo
       endif
    else
       do k=1,nz
          do j=1,ny
             do i=ndx_d,nfx_d
                ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dx
                Krho(i,j,k) = Drho(i-1,j,k) -Drho(i,j,k)
                Krhou(i,j,k)= Drhou(i-1,j,k)-Drhou(i,j,k)
                Krhov(i,j,k)= Drhov(i-1,j,k)-Drhov(i,j,k)
                Krhow(i,j,k)= Drhow(i-1,j,k)-Drhow(i,j,k)
                Krhoe(i,j,k)= Drhoe(i-1,j,k)-Drhoe(i,j,k)
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
 
    ! Filtering in y-direction [compute Drho at j+1/2]
    ! ========================

    if (is_boundary(2,1)) then
       ! 1st-order centered dissipation
       ! ------------------------------
       if (is_bc_wall(2,1)) then
          j=1
          do k=1,nz
             do i=1,nx
                xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)
                !uvar(i,j,k,2)=xk2

                Drho(i,j,k) = (rho_n(i,j+1,k)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i,j+1,k)-rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
       else
          j=1
          do k=1,nz
             do i=1,nx
                xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)
                !uvar(i,j,k,2)=xk2

                Drho(i,j,k) = (rho_n(i,j+1,k)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i,j+1,k)-rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
       endif
       
       ! 3rd-order centered dissipation
       ! ------------------------------
       j=2
       do k=1,nz
          do i=1,nx
             xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))
             xk4=max(0.0_wp,xnu4+xk0)
             !uvar(i,j,k,2)=xk4/xnu4
             
             Drho(i,j,k) = (rho_n(i,j+1,k)-rho_n(i,j,k))*xk0      &
                         + (d4(1)*(rho_n(i,j+1,k)-rho_n(i,j  ,k)) &
                           +d4(2)*(rho_n(i,j+2,k)-rho_n(i,j-1,k)))*xk4
             Drhou(i,j,k)= (rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk0      &
                         + (d4(1)*(rhou_n(i,j+1,k)-rhou_n(i,j  ,k)) &
                           +d4(2)*(rhou_n(i,j+2,k)-rhou_n(i,j-1,k)))*xk4
             Drhov(i,j,k)= (rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk0      &
                         + (d4(1)*(rhov_n(i,j+1,k)-rhov_n(i,j  ,k)) &
                           +d4(2)*(rhov_n(i,j+2,k)-rhov_n(i,j-1,k)))*xk4
             Drhow(i,j,k)= (rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk0      &
                         + (d4(1)*(rhow_n(i,j+1,k)-rhow_n(i,j  ,k)) &
                           +d4(2)*(rhow_n(i,j+2,k)-rhow_n(i,j-1,k)))*xk4
             Drhoe(i,j,k)= (rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk0      &
                         + (d4(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j  ,k)) &
                           +d4(2)*(rhoe_n(i,j+2,k)-rhoe_n(i,j-1,k)))*xk4
          enddo
       enddo
       
       ! 5th-order centered dissipation
       ! ------------------------------
       j=3
       do k=1,nz
          do i=1,nx
             xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))
             !xk6=max(0.0_wp,xnu6+xk0/5.0_wp)
             xk6=max(0.0_wp,xnu6+xk0)
             !uvar(i,j,k,2)=xk6/xnu6
            
             Drho(i,j,k) = (rho_n(i,j+1,k)-rho_n(i,j,k))*xk0      &
                         + (d6(1)*(rho_n(i,j+1,k)-rho_n(i,j  ,k)) &
                           +d6(2)*(rho_n(i,j+2,k)-rho_n(i,j-1,k)) &
                           +d6(3)*(rho_n(i,j+3,k)-rho_n(i,j-2,k)))*xk6
             Drhou(i,j,k)= (rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk0      &
                         + (d6(1)*(rhou_n(i,j+1,k)-rhou_n(i,j  ,k)) &
                           +d6(2)*(rhou_n(i,j+2,k)-rhou_n(i,j-1,k)) &
                           +d6(3)*(rhou_n(i,j+3,k)-rhou_n(i,j-2,k)))*xk6
             Drhov(i,j,k)= (rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk0      &
                         + (d6(1)*(rhov_n(i,j+1,k)-rhov_n(i,j  ,k)) &
                           +d6(2)*(rhov_n(i,j+2,k)-rhov_n(i,j-1,k)) &
                           +d6(3)*(rhov_n(i,j+3,k)-rhov_n(i,j-2,k)))*xk6
             Drhow(i,j,k)= (rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk0      &
                         + (d6(1)*(rhow_n(i,j+1,k)-rhow_n(i,j  ,k)) &
                           +d6(2)*(rhow_n(i,j+2,k)-rhow_n(i,j-1,k)) &
                           +d6(3)*(rhow_n(i,j+3,k)-rhow_n(i,j-2,k)))*xk6
             Drhoe(i,j,k)= (rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk0      &
                         + (d6(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j  ,k)) &
                           +d6(2)*(rhoe_n(i,j+2,k)-rhoe_n(i,j-1,k)) &
                           +d6(3)*(rhoe_n(i,j+3,k)-rhoe_n(i,j-2,k)))*xk6
          enddo
       enddo
       
       ! 7th-order centered dissipation
       ! ------------------------------
       j=4
       do k=1,nz
          do i=1,nx
             xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))
             !xk8=max(0.0_wp,xnu8+3.0_wp/70.0_wp*xk0)
             xk8=max(0.0_wp,xnu8+xk0)
             !uvar(i,j,k,2)=xk8/xnu8
             
             Drho(i,j,k) = (rho_n(i,j+1,k)-rho_n(i,j,k))*xk0      &
                         + (d8(1)*(rho_n(i,j+1,k)-rho_n(i,j  ,k)) &
                           +d8(2)*(rho_n(i,j+2,k)-rho_n(i,j-1,k)) &
                           +d8(3)*(rho_n(i,j+3,k)-rho_n(i,j-2,k)) &
                           +d8(4)*(rho_n(i,j+4,k)-rho_n(i,j-3,k)))*xk8
             Drhou(i,j,k)= (rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk0      &
                         + (d8(1)*(rhou_n(i,j+1,k)-rhou_n(i,j  ,k)) &
                           +d8(2)*(rhou_n(i,j+2,k)-rhou_n(i,j-1,k)) &
                           +d8(3)*(rhou_n(i,j+3,k)-rhou_n(i,j-2,k)) &
                           +d8(4)*(rhou_n(i,j+4,k)-rhou_n(i,j-3,k)))*xk8
             Drhov(i,j,k)= (rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk0      &
                         + (d8(1)*(rhov_n(i,j+1,k)-rhov_n(i,j  ,k)) &
                           +d8(2)*(rhov_n(i,j+2,k)-rhov_n(i,j-1,k)) &
                           +d8(3)*(rhov_n(i,j+3,k)-rhov_n(i,j-2,k)) &
                           +d8(4)*(rhov_n(i,j+4,k)-rhov_n(i,j-3,k)))*xk8
             Drhow(i,j,k)= (rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk0      &
                         + (d8(1)*(rhow_n(i,j+1,k)-rhow_n(i,j  ,k)) &
                           +d8(2)*(rhow_n(i,j+2,k)-rhow_n(i,j-1,k)) &
                           +d8(3)*(rhow_n(i,j+3,k)-rhow_n(i,j-2,k)) &
                           +d8(4)*(rhow_n(i,j+4,k)-rhow_n(i,j-3,k)))*xk8
             Drhoe(i,j,k)= (rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk0      &
                         + (d8(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j  ,k)) &
                           +d8(2)*(rhoe_n(i,j+2,k)-rhoe_n(i,j-1,k)) &
                           +d8(3)*(rhoe_n(i,j+3,k)-rhoe_n(i,j-2,k)) &
                           +d8(4)*(rhoe_n(i,j+4,k)-rhoe_n(i,j-3,k)))*xk8
          enddo
       enddo
    endif
    
    ! Interior points: 9th-order centered dissipation
    ! -----------------------------------------------
    do k=1,nz
       do j=ndy_di,nfy_di
          do i=1,nx
             xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))
             xk10=max(0.0_wp,xnu10+xk0/105.0_wp)
             !if (j>=1) uvar(i,j,k,2)=xk10/xnu10
             
             Drho(i,j,k) = (rho_n(i,j+1,k)-rho_n(i,j,k))*xk0       &
                         + (d10(1)*(rho_n(i,j+1,k)-rho_n(i,j  ,k)) &
                           +d10(2)*(rho_n(i,j+2,k)-rho_n(i,j-1,k)) &
                           +d10(3)*(rho_n(i,j+3,k)-rho_n(i,j-2,k)) &
                           +d10(4)*(rho_n(i,j+4,k)-rho_n(i,j-3,k)) &
                           +d10(5)*(rho_n(i,j+5,k)-rho_n(i,j-4,k)))*xk10

             Drhou(i,j,k)= (rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk0       &
                         + (d10(1)*(rhou_n(i,j+1,k)-rhou_n(i,j  ,k)) &
                           +d10(2)*(rhou_n(i,j+2,k)-rhou_n(i,j-1,k)) &
                           +d10(3)*(rhou_n(i,j+3,k)-rhou_n(i,j-2,k)) &
                           +d10(4)*(rhou_n(i,j+4,k)-rhou_n(i,j-3,k)) &
                           +d10(5)*(rhou_n(i,j+5,k)-rhou_n(i,j-4,k)))*xk10

             Drhov(i,j,k)= (rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk0       &
                         + (d10(1)*(rhov_n(i,j+1,k)-rhov_n(i,j  ,k)) &
                           +d10(2)*(rhov_n(i,j+2,k)-rhov_n(i,j-1,k)) &
                           +d10(3)*(rhov_n(i,j+3,k)-rhov_n(i,j-2,k)) &
                           +d10(4)*(rhov_n(i,j+4,k)-rhov_n(i,j-3,k)) &
                           +d10(5)*(rhov_n(i,j+5,k)-rhov_n(i,j-4,k)))*xk10

             Drhow(i,j,k)= (rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk0       &
                         + (d10(1)*(rhow_n(i,j+1,k)-rhow_n(i,j  ,k)) &
                           +d10(2)*(rhow_n(i,j+2,k)-rhow_n(i,j-1,k)) &
                           +d10(3)*(rhow_n(i,j+3,k)-rhow_n(i,j-2,k)) &
                           +d10(4)*(rhow_n(i,j+4,k)-rhow_n(i,j-3,k)) &
                           +d10(5)*(rhow_n(i,j+5,k)-rhow_n(i,j-4,k)))*xk10

             Drhoe(i,j,k)= (rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk0       &
                         + (d10(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j  ,k)) &
                           +d10(2)*(rhoe_n(i,j+2,k)-rhoe_n(i,j-1,k)) &
                           +d10(3)*(rhoe_n(i,j+3,k)-rhoe_n(i,j-2,k)) &
                           +d10(4)*(rhoe_n(i,j+4,k)-rhoe_n(i,j-3,k)) &
                           +d10(5)*(rhoe_n(i,j+5,k)-rhoe_n(i,j-4,k)))*xk10
          enddo
       enddo
    enddo

    if (is_boundary(2,2)) then
       ! 7th-order centered dissipation
       ! ------------------------------
       j=ny-4
       do k=1,nz
          do i=1,nx
             xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))
             !xk8=max(0.0_wp,xnu8+3.0_wp/70.0_wp*xk0)
             xk8=max(0.0_wp,xnu8+xk0)
             !uvar(i,j,k,2)=xk8/xnu8
             
             Drho(i,j,k) = (rho_n(i,j+1,k)-rho_n(i,j,k))*xk0      &
                         + (d8(1)*(rho_n(i,j+1,k)-rho_n(i,j  ,k)) &
                           +d8(2)*(rho_n(i,j+2,k)-rho_n(i,j-1,k)) &
                           +d8(3)*(rho_n(i,j+3,k)-rho_n(i,j-2,k)) &
                           +d8(4)*(rho_n(i,j+4,k)-rho_n(i,j-3,k)))*xk8
             Drhou(i,j,k)= (rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk0      &
                         + (d8(1)*(rhou_n(i,j+1,k)-rhou_n(i,j  ,k)) &
                           +d8(2)*(rhou_n(i,j+2,k)-rhou_n(i,j-1,k)) &
                           +d8(3)*(rhou_n(i,j+3,k)-rhou_n(i,j-2,k)) &
                           +d8(4)*(rhou_n(i,j+4,k)-rhou_n(i,j-3,k)))*xk8
             Drhov(i,j,k)= (rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk0      &
                         + (d8(1)*(rhov_n(i,j+1,k)-rhov_n(i,j  ,k)) &
                           +d8(2)*(rhov_n(i,j+2,k)-rhov_n(i,j-1,k)) &
                           +d8(3)*(rhov_n(i,j+3,k)-rhov_n(i,j-2,k)) &
                           +d8(4)*(rhov_n(i,j+4,k)-rhov_n(i,j-3,k)))*xk8
             Drhow(i,j,k)= (rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk0      &
                         + (d8(1)*(rhow_n(i,j+1,k)-rhow_n(i,j  ,k)) &
                           +d8(2)*(rhow_n(i,j+2,k)-rhow_n(i,j-1,k)) &
                           +d8(3)*(rhow_n(i,j+3,k)-rhow_n(i,j-2,k)) &
                           +d8(4)*(rhow_n(i,j+4,k)-rhow_n(i,j-3,k)))*xk8
             Drhoe(i,j,k)= (rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk0      &
                         + (d8(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j  ,k)) &
                           +d8(2)*(rhoe_n(i,j+2,k)-rhoe_n(i,j-1,k)) &
                           +d8(3)*(rhoe_n(i,j+3,k)-rhoe_n(i,j-2,k)) &
                           +d8(4)*(rhoe_n(i,j+4,k)-rhoe_n(i,j-3,k)))*xk8
          enddo
       enddo
       
       ! 5th-order centered dissipation
       ! ------------------------------
       j=ny-3
       do k=1,nz
          do i=1,nx
             xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))
             !xk6=max(0.0_wp,xnu6+xk0/5.0_wp)
             xk6=max(0.0_wp,xnu6+xk0)
             !uvar(i,j,k,2)=xk6/xnu6
             
             Drho(i,j,k) = (rho_n(i,j+1,k)-rho_n(i,j,k))*xk0      &
                         + (d6(1)*(rho_n(i,j+1,k)-rho_n(i,j  ,k)) &
                           +d6(2)*(rho_n(i,j+2,k)-rho_n(i,j-1,k)) &
                           +d6(3)*(rho_n(i,j+3,k)-rho_n(i,j-2,k)))*xk6
             Drhou(i,j,k)= (rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk0      &
                         + (d6(1)*(rhou_n(i,j+1,k)-rhou_n(i,j  ,k)) &
                           +d6(2)*(rhou_n(i,j+2,k)-rhou_n(i,j-1,k)) &
                           +d6(3)*(rhou_n(i,j+3,k)-rhou_n(i,j-2,k)))*xk6
             Drhov(i,j,k)= (rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk0      &
                         + (d6(1)*(rhov_n(i,j+1,k)-rhov_n(i,j  ,k)) &
                           +d6(2)*(rhov_n(i,j+2,k)-rhov_n(i,j-1,k)) &
                           +d6(3)*(rhov_n(i,j+3,k)-rhov_n(i,j-2,k)))*xk6
             Drhow(i,j,k)= (rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk0      &
                         + (d6(1)*(rhow_n(i,j+1,k)-rhow_n(i,j  ,k)) &
                           +d6(2)*(rhow_n(i,j+2,k)-rhow_n(i,j-1,k)) &
                           +d6(3)*(rhow_n(i,j+3,k)-rhow_n(i,j-2,k)))*xk6
             Drhoe(i,j,k)= (rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk0      &
                         + (d6(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j  ,k)) &
                           +d6(2)*(rhoe_n(i,j+2,k)-rhoe_n(i,j-1,k)) &
                           +d6(3)*(rhoe_n(i,j+3,k)-rhoe_n(i,j-2,k)))*xk6
          enddo
       enddo
       
       ! 3rd-order centered dissipation
       ! ------------------------------
       j=ny-2
       do k=1,nz
          do i=1,nx
             xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))
             xk4=max(0.0_wp,xnu4+xk0)
             !uvar(i,j,k,2)=xk4/xnu4
             
             Drho(i,j,k) = (rho_n(i,j+1,k)-rho_n(i,j,k))*xk0      &
                         + (d4(1)*(rho_n(i,j+1,k)-rho_n(i,j  ,k)) &
                           +d4(2)*(rho_n(i,j+2,k)-rho_n(i,j-1,k)))*xk4
             Drhou(i,j,k)= (rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk0      &
                         + (d4(1)*(rhou_n(i,j+1,k)-rhou_n(i,j  ,k)) &
                           +d4(2)*(rhou_n(i,j+2,k)-rhou_n(i,j-1,k)))*xk4
             Drhov(i,j,k)= (rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk0      &
                         + (d4(1)*(rhov_n(i,j+1,k)-rhov_n(i,j  ,k)) &
                           +d4(2)*(rhov_n(i,j+2,k)-rhov_n(i,j-1,k)))*xk4
             Drhow(i,j,k)= (rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk0      &
                         + (d4(1)*(rhow_n(i,j+1,k)-rhow_n(i,j  ,k)) &
                           +d4(2)*(rhow_n(i,j+2,k)-rhow_n(i,j-1,k)))*xk4
             Drhoe(i,j,k)= (rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk0      &
                         + (d4(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j  ,k)) &
                           +d4(2)*(rhoe_n(i,j+2,k)-rhoe_n(i,j-1,k)))*xk4
          enddo
       enddo
       
       ! 1st-order centered dissipation
       ! ------------------------------
       if (is_bc_wall(2,2)) then
          j=ny-1
          do k=1,nz
             do i=1,nx
                xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)
                !uvar(i,j,k,2)=xk2

                Drho(i,j,k) = (rho_n(i,j+1,k)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i,j+1,k)-rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
       else
          j=ny-1
          do k=1,nz
             do i=1,nx
                xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)
                !uvar(i,j,k,2)=xk2

                Drho(i,j,k) = (rho_n(i,j+1,k)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i,j+1,k)-rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
       endif
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
                   Krho(i,j,k) = Krho(i,j,k) +mul*(Drho(i,j-1,k) -Drho(i,j,k))
                   Krhou(i,j,k)= Krhou(i,j,k)+mul*(Drhou(i,j-1,k)-Drhou(i,j,k))
                   Krhov(i,j,k)= Krhov(i,j,k)+mul*(Drhov(i,j-1,k)-Drhov(i,j,k))
                   Krhow(i,j,k)= Krhow(i,j,k)+mul*(Drhow(i,j-1,k)-Drhow(i,j,k))
                   Krhoe(i,j,k)= Krhoe(i,j,k)+mul*(Drhoe(i,j-1,k)-Drhoe(i,j,k))
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
                   Krho(i,j,k) = Krho(i,j,k) +mul*(Drho(i,j-1,k) -Drho(i,j,k))
                   Krhou(i,j,k)= Krhou(i,j,k)+mul*(Drhou(i,j-1,k)-Drhou(i,j,k))
                   Krhov(i,j,k)= Krhov(i,j,k)+mul*(Drhov(i,j-1,k)-Drhov(i,j,k))
                   Krhow(i,j,k)= Krhow(i,j,k)+mul*(Drhow(i,j-1,k)-Drhow(i,j,k))
                   Krhoe(i,j,k)= Krhoe(i,j,k)+mul*(Drhoe(i,j-1,k)-Drhoe(i,j,k))
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
                   Krho(i,j,k) = Krho(i,j,k) +mul*(Drho(i,j-1,k) -Drho(i,j,k))
                   Krhou(i,j,k)= Krhou(i,j,k)+mul*(Drhou(i,j-1,k)-Drhou(i,j,k))
                   Krhov(i,j,k)= Krhov(i,j,k)+mul*(Drhov(i,j-1,k)-Drhov(i,j,k))
                   Krhow(i,j,k)= Krhow(i,j,k)+mul*(Drhow(i,j-1,k)-Drhow(i,j,k))
                   Krhoe(i,j,k)= Krhoe(i,j,k)+mul*(Drhoe(i,j-1,k)-Drhoe(i,j,k))
                enddo
             enddo
          enddo
       endif
    else
       do k=1,nz
          do j=ndy_d,nfy_d
             do i=1,nx
                ! compute dissipation term: delta_j/delta_y=(D_{j+1/2}-D_{j-1/2})/dy
                Krho(i,j,k) = Krho(i,j,k) +Drho(i,j-1,k) -Drho(i,j,k)
                Krhou(i,j,k)= Krhou(i,j,k)+Drhou(i,j-1,k)-Drhou(i,j,k)
                Krhov(i,j,k)= Krhov(i,j,k)+Drhov(i,j-1,k)-Drhov(i,j,k)
                Krhow(i,j,k)= Krhow(i,j,k)+Drhow(i,j-1,k)-Drhow(i,j,k)
                Krhoe(i,j,k)= Krhoe(i,j,k)+Drhoe(i,j-1,k)-Drhoe(i,j,k)
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
 
       ! Filtering in k-direction [compute Drho at k+1/2]
       ! ========================
       if (is_boundary(3,1)) then
          ! 1st-order centered dissipation
          ! ------------------------------
          k=1
          do j=1,ny
             do i=1,nx
                xk0=xnu0*max(psens(i,j,k+1),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)
             
                Drho(i,j,k) = (rho_n(i,j,k+1)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i,j,k+1)-rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i,j,k+1)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i,j,k+1)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i,j,k+1)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i,j,k+1)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
          
          ! 3rd-order centered dissipation
          ! ------------------------------
          k=2
          do j=1,ny
             do i=1,nx
                xk0=xnu0*max(psens(i,j,k+1),psens(i,j,k))
                xk4=max(0.0_wp,xnu4+xk0)
             
                Drho(i,j,k) = (rho_n(i,j,k+1)-rho_n(i,j,k))*xk0      &
                            + (d4(1)*(rho_n(i,j,k+1)-rho_n(i,j,k  )) &
                              +d4(2)*(rho_n(i,j,k+2)-rho_n(i,j,k-1)))*xk4
                Drhou(i,j,k)= (rhou_n(i,j,k+1)-rhou_n(i,j,k))*xk0      &
                            + (d4(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k  )) &
                              +d4(2)*(rhou_n(i,j,k+2)-rhou_n(i,j,k-1)))*xk4
                Drhov(i,j,k)= (rhov_n(i,j,k+1)-rhov_n(i,j,k))*xk0      &
                            + (d4(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k  )) &
                              +d4(2)*(rhov_n(i,j,k+2)-rhov_n(i,j,k-1)))*xk4
                Drhow(i,j,k)= (rhow_n(i,j,k+1)-rhow_n(i,j,k))*xk0      &
                            + (d4(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k  )) &
                              +d4(2)*(rhow_n(i,j,k+2)-rhow_n(i,j,k-1)))*xk4
                Drhoe(i,j,k)= (rhoe_n(i,j,k+1)-rhoe_n(i,j,k))*xk0      &
                            + (d4(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k  )) &
                              +d4(2)*(rhoe_n(i,j,k+2)-rhoe_n(i,j,k-1)))*xk4
             enddo
          enddo
          
          ! 5th-order centered dissipation
          ! ------------------------------
          k=3
          do j=1,ny
             do i=1,nx
                xk0=xnu0*max(psens(i,j,k+1),psens(i,j,k))
                !xk6=max(0.0_wp,xnu6+xk0/5.0_wp)
                xk6=max(0.0_wp,xnu6+xk0)
             
                Drho(i,j,k) = (rho_n(i,j,k+1)-rho_n(i,j,k))*xk0      &
                            + (d6(1)*(rho_n(i,j,k+1)-rho_n(i,j,k  )) &
                              +d6(2)*(rho_n(i,j,k+2)-rho_n(i,j,k-1)) &
                              +d6(3)*(rho_n(i,j,k+3)-rho_n(i,j,k-2)))*xk6
                Drhou(i,j,k)= (rhou_n(i,j,k+1)-rhou_n(i,j,k))*xk0      &
                            + (d6(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k  )) &
                              +d6(2)*(rhou_n(i,j,k+2)-rhou_n(i,j,k-1)) &
                              +d6(3)*(rhou_n(i,j,k+3)-rhou_n(i,j,k-2)))*xk6
                Drhov(i,j,k)= (rhov_n(i,j,k+1)-rhov_n(i,j,k))*xk0      &
                            + (d6(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k  )) &
                              +d6(2)*(rhov_n(i,j,k+2)-rhov_n(i,j,k-1)) &
                              +d6(3)*(rhov_n(i,j,k+3)-rhov_n(i,j,k-2)))*xk6
                Drhow(i,j,k)= (rhow_n(i,j,k+1)-rhow_n(i,j,k))*xk0      &
                            + (d6(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k  )) &
                              +d6(2)*(rhow_n(i,j,k+2)-rhow_n(i,j,k-1)) &
                              +d6(3)*(rhow_n(i,j,k+3)-rhow_n(i,j,k-2)))*xk6
                Drhoe(i,j,k)= (rhoe_n(i,j,k+1)-rhoe_n(i,j,k))*xk0      &
                            + (d6(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k  )) &
                              +d6(2)*(rhoe_n(i,j,k+2)-rhoe_n(i,j,k-1)) &
                              +d6(3)*(rhoe_n(i,j,k+3)-rhoe_n(i,j,k-2)))*xk6
             enddo
          enddo
          
          ! 7th-order centered dissipation
          ! ------------------------------
          k=4
          do j=1,ny
             do i=1,nx
                xk0=xnu0*max(psens(i,j,k+1),psens(i,j,k))
                !xk8=max(0.0_wp,xnu8+3.0_wp/70.0_wp*xk0)
                xk8=max(0.0_wp,xnu8+xk0)
             
                Drho(i,j,k) = (rho_n(i,j,k+1)-rho_n(i,j,k))*xk0      &
                            + (d8(1)*(rho_n(i,j,k+1)-rho_n(i,j,k  )) &
                              +d8(2)*(rho_n(i,j,k+2)-rho_n(i,j,k-1)) &
                              +d8(3)*(rho_n(i,j,k+3)-rho_n(i,j,k-2)) &
                              +d8(4)*(rho_n(i,j,k+4)-rho_n(i,j,k-3)))*xk8
                Drhou(i,j,k)= (rhou_n(i,j,k+1)-rhou_n(i,j,k))*xk0      &
                            + (d8(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k  )) &
                              +d8(2)*(rhou_n(i,j,k+2)-rhou_n(i,j,k-1)) &
                              +d8(3)*(rhou_n(i,j,k+3)-rhou_n(i,j,k-2)) &
                              +d8(4)*(rhou_n(i,j,k+4)-rhou_n(i,j,k-3)))*xk8
                Drhov(i,j,k)= (rhov_n(i,j,k+1)-rhov_n(i,j,k))*xk0      &
                            + (d8(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k  )) &
                              +d8(2)*(rhov_n(i,j,k+2)-rhov_n(i,j,k-1)) &
                              +d8(3)*(rhov_n(i,j,k+3)-rhov_n(i,j,k-2)) &
                              +d8(4)*(rhov_n(i,j,k+4)-rhov_n(i,j,k-3)))*xk8
                Drhow(i,j,k)= (rhow_n(i,j,k+1)-rhow_n(i,j,k))*xk0      &
                            + (d8(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k  )) &
                              +d8(2)*(rhow_n(i,j,k+2)-rhow_n(i,j,k-1)) &
                              +d8(3)*(rhow_n(i,j,k+3)-rhow_n(i,j,k-2)) &
                              +d8(4)*(rhow_n(i,j,k+4)-rhow_n(i,j,k-3)))*xk8
                Drhoe(i,j,k)= (rhoe_n(i,j,k+1)-rhoe_n(i,j,k))*xk0      &
                            + (d8(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k  )) &
                              +d8(2)*(rhoe_n(i,j,k+2)-rhoe_n(i,j,k-1)) &
                              +d8(3)*(rhoe_n(i,j,k+3)-rhoe_n(i,j,k-2)) &
                              +d8(4)*(rhoe_n(i,j,k+4)-rhoe_n(i,j,k-3)))*xk8
             enddo
          enddo
       endif
       
       ! Interior points: 9th-order centered dissipation
       ! -----------------------------------------------
       do k=ndz_di,nfz_di
          do j=1,ny
             do i=1,nx
                xk0=xnu0*max(psens(i,j,k+1),psens(i,j,k))
                xk10=max(0.0_wp,xnu10+xk0/105.0_wp)
             
                Drho(i,j,k) = (rho_n(i,j,k+1)-rho_n(i,j,k))*xk0       &
                            + (d10(1)*(rho_n(i,j,k+1)-rho_n(i,j,k  )) &
                              +d10(2)*(rho_n(i,j,k+2)-rho_n(i,j,k-1)) &
                              +d10(3)*(rho_n(i,j,k+3)-rho_n(i,j,k-2)) &
                              +d10(4)*(rho_n(i,j,k+4)-rho_n(i,j,k-3)) &
                              +d10(5)*(rho_n(i,j,k+5)-rho_n(i,j,k-4)))*xk10

                Drhou(i,j,k)= (rhou_n(i,j,k+1)-rhou_n(i,j,k))*xk0       &
                            + (d10(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k  )) &
                              +d10(2)*(rhou_n(i,j,k+2)-rhou_n(i,j,k-1)) &
                              +d10(3)*(rhou_n(i,j,k+3)-rhou_n(i,j,k-2)) &
                              +d10(4)*(rhou_n(i,j,k+4)-rhou_n(i,j,k-3)) &
                              +d10(5)*(rhou_n(i,j,k+5)-rhou_n(i,j,k-4)))*xk10

                Drhov(i,j,k)= (rhov_n(i,j,k+1)-rhov_n(i,j,k))*xk0       &
                            + (d10(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k  )) &
                              +d10(2)*(rhov_n(i,j,k+2)-rhov_n(i,j,k-1)) &
                              +d10(3)*(rhov_n(i,j,k+3)-rhov_n(i,j,k-2)) &
                              +d10(4)*(rhov_n(i,j,k+4)-rhov_n(i,j,k-3)) &
                              +d10(5)*(rhov_n(i,j,k+5)-rhov_n(i,j,k-4)))*xk10

                Drhow(i,j,k)= (rhow_n(i,j,k+1)-rhow_n(i,j,k))*xk0       &
                            + (d10(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k  )) &
                              +d10(2)*(rhow_n(i,j,k+2)-rhow_n(i,j,k-1)) &
                              +d10(3)*(rhow_n(i,j,k+3)-rhow_n(i,j,k-2)) &
                              +d10(4)*(rhow_n(i,j,k+4)-rhow_n(i,j,k-3)) &
                              +d10(5)*(rhow_n(i,j,k+5)-rhow_n(i,j,k-4)))*xk10

                Drhoe(i,j,k)= (rhoe_n(i,j,k+1)-rhoe_n(i,j,k))*xk0       &
                            + (d10(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k  )) &
                              +d10(2)*(rhoe_n(i,j,k+2)-rhoe_n(i,j,k-1)) &
                              +d10(3)*(rhoe_n(i,j,k+3)-rhoe_n(i,j,k-2)) &
                              +d10(4)*(rhoe_n(i,j,k+4)-rhoe_n(i,j,k-3)) &
                              +d10(5)*(rhoe_n(i,j,k+5)-rhoe_n(i,j,k-4)))*xk10
             enddo
          enddo
       enddo

       if (is_boundary(3,2)) then
          ! 7th-order centered dissipation
          ! ------------------------------
          k=nz-4
          do j=1,ny
             do i=1,nx
             
                xk0 = xnu0*max(psens(i,j,k+1),psens(i,j,k))
                !xk8=max(0.0_wp,xnu8+3.0_wp/70.0_wp*xk0)
                xk8=max(0.0_wp,xnu8+xk0)
             
                Drho(i,j,k) = (rho_n(i,j,k+1)-rho_n(i,j,k))*xk0      &
                            + (d8(1)*(rho_n(i,j,k+1)-rho_n(i,j,k  )) &
                              +d8(2)*(rho_n(i,j,k+2)-rho_n(i,j,k-1)) &
                              +d8(3)*(rho_n(i,j,k+3)-rho_n(i,j,k-2)) &
                              +d8(4)*(rho_n(i,j,k+4)-rho_n(i,j,k-3)))*xk8
                Drhou(i,j,k)= (rhou_n(i,j,k+1)-rhou_n(i,j,k))*xk0      &
                            + (d8(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k  )) &
                              +d8(2)*(rhou_n(i,j,k+2)-rhou_n(i,j,k-1)) &
                              +d8(3)*(rhou_n(i,j,k+3)-rhou_n(i,j,k-2)) &
                              +d8(4)*(rhou_n(i,j,k+4)-rhou_n(i,j,k-3)))*xk8
                Drhov(i,j,k)= (rhov_n(i,j,k+1)-rhov_n(i,j,k))*xk0      &
                            + (d8(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k  )) &
                              +d8(2)*(rhov_n(i,j,k+2)-rhov_n(i,j,k-1)) &
                              +d8(3)*(rhov_n(i,j,k+3)-rhov_n(i,j,k-2)) &
                              +d8(4)*(rhov_n(i,j,k+4)-rhov_n(i,j,k-3)))*xk8
                Drhow(i,j,k)= (rhow_n(i,j,k+1)-rhow_n(i,j,k))*xk0      &
                            + (d8(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k  )) &
                              +d8(2)*(rhow_n(i,j,k+2)-rhow_n(i,j,k-1)) &
                              +d8(3)*(rhow_n(i,j,k+3)-rhow_n(i,j,k-2)) &
                              +d8(4)*(rhow_n(i,j,k+4)-rhow_n(i,j,k-3)))*xk8
                Drhoe(i,j,k)= (rhoe_n(i,j,k+1)-rhoe_n(i,j,k))*xk0      &
                            + (d8(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k  )) &
                              +d8(2)*(rhoe_n(i,j,k+2)-rhoe_n(i,j,k-1)) &
                              +d8(3)*(rhoe_n(i,j,k+3)-rhoe_n(i,j,k-2)) &
                              +d8(4)*(rhoe_n(i,j,k+4)-rhoe_n(i,j,k-3)))*xk8
             enddo
          enddo
          
          ! 5th-order centered dissipation
          ! ------------------------------
          k=nz-3
          do j=1,ny
             do i=1,nx
                xk0=xnu0*max(psens(i,j,k+1),psens(i,j,k))
                !xk6=max(0.0_wp,xnu6+xk0/5.0_wp)
                xk6=max(0.0_wp,xnu6+xk0)
             
                Drho(i,j,k) = (rho_n(i,j,k+1)-rho_n(i,j,k))*xk0      &
                            + (d6(1)*(rho_n(i,j,k+1)-rho_n(i,j,k  )) &
                              +d6(2)*(rho_n(i,j,k+2)-rho_n(i,j,k-1)) &
                              +d6(3)*(rho_n(i,j,k+3)-rho_n(i,j,k-2)))*xk6
                Drhou(i,j,k)= (rhou_n(i,j,k+1)-rhou_n(i,j,k))*xk0      &
                            + (d6(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k  )) &
                              +d6(2)*(rhou_n(i,j,k+2)-rhou_n(i,j,k-1)) &
                              +d6(3)*(rhou_n(i,j,k+3)-rhou_n(i,j,k-2)))*xk6
                Drhov(i,j,k)= (rhov_n(i,j,k+1)-rhov_n(i,j,k))*xk0      &
                            + (d6(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k  )) &
                              +d6(2)*(rhov_n(i,j,k+2)-rhov_n(i,j,k-1)) &
                              +d6(3)*(rhov_n(i,j,k+3)-rhov_n(i,j,k-2)))*xk6
                Drhow(i,j,k)= (rhow_n(i,j,k+1)-rhow_n(i,j,k))*xk0      &
                            + (d6(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k  )) &
                              +d6(2)*(rhow_n(i,j,k+2)-rhow_n(i,j,k-1)) &
                              +d6(3)*(rhow_n(i,j,k+3)-rhow_n(i,j,k-2)))*xk6
                Drhoe(i,j,k)= (rhoe_n(i,j,k+1)-rhoe_n(i,j,k))*xk0      &
                            + (d6(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k  )) &
                              +d6(2)*(rhoe_n(i,j,k+2)-rhoe_n(i,j,k-1)) &
                              +d6(3)*(rhoe_n(i,j,k+3)-rhoe_n(i,j,k-2)))*xk6
             enddo
          enddo
          
          ! 3rd-order centered dissipation
          ! ------------------------------
          k=nz-2
          do j=1,ny
             do i=1,nx
                xk0=xnu0*max(psens(i,j,k+1),psens(i,j,k))
                xk4=max(0.0_wp,xnu4+xk0)
             
                Drho(i,j,k) = (rho_n(i,j,k+1)-rho_n(i,j,k))*xk0      &
                            + (d4(1)*(rho_n(i,j,k+1)-rho_n(i,j,k  )) &
                              +d4(2)*(rho_n(i,j,k+2)-rho_n(i,j,k-1)))*xk4
                Drhou(i,j,k)= (rhou_n(i,j,k+1)-rhou_n(i,j,k))*xk0 &
                            + (d4(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k  )) &
                              +d4(2)*(rhou_n(i,j,k+2)-rhou_n(i,j,k-1)))*xk4
                Drhov(i,j,k)= (rhov_n(i,j,k+1)-rhov_n(i,j,k))*xk0      &
                            + (d4(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k  )) &
                              +d4(2)*(rhov_n(i,j,k+2)-rhov_n(i,j,k-1)))*xk4
                Drhow(i,j,k)= (rhow_n(i,j,k+1)-rhow_n(i,j,k))*xk0      &
                            + (d4(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k  )) &
                              +d4(2)*(rhow_n(i,j,k+2)-rhow_n(i,j,k-1)))*xk4
                Drhoe(i,j,k)= (rhoe_n(i,j,k+1)-rhoe_n(i,j,k))*xk0      &
                            + (d4(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k  )) &
                              +d4(2)*(rhoe_n(i,j,k+2)-rhoe_n(i,j,k-1)))*xk4
             enddo
          enddo
          
          ! 1st-order centered dissipation
          ! ------------------------------
          k=nz-1
          do j=1,ny
             do i=1,nx
                xk0=xnu0*max(psens(i,j,k+1),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)
             
                Drho(i,j,k) = (rho_n(i,j,k+1)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i,j,k+1)-rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i,j,k+1)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i,j,k+1)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i,j,k+1)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i,j,k+1)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
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
                      Krho(i,j,k) = Krho(i,j,k) +mul*(Drho(i,j,k-1) -Drho(i,j,k))
                      Krhou(i,j,k)= Krhou(i,j,k)+mul*(Drhou(i,j,k-1)-Drhou(i,j,k))
                      Krhov(i,j,k)= Krhov(i,j,k)+mul*(Drhov(i,j,k-1)-Drhov(i,j,k))
                      Krhow(i,j,k)= Krhow(i,j,k)+mul*(Drhow(i,j,k-1)-Drhow(i,j,k))
                      Krhoe(i,j,k)= Krhoe(i,j,k)+mul*(Drhoe(i,j,k-1)-Drhoe(i,j,k))
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
                      Krho(i,j,k) = Krho(i,j,k) +mul*(Drho(i,j,k-1) -Drho(i,j,k))
                      Krhou(i,j,k)= Krhou(i,j,k)+mul*(Drhou(i,j,k-1)-Drhou(i,j,k))
                      Krhov(i,j,k)= Krhov(i,j,k)+mul*(Drhov(i,j,k-1)-Drhov(i,j,k))
                      Krhow(i,j,k)= Krhow(i,j,k)+mul*(Drhow(i,j,k-1)-Drhow(i,j,k))
                      Krhoe(i,j,k)= Krhoe(i,j,k)+mul*(Drhoe(i,j,k-1)-Drhoe(i,j,k))
                   enddo
                enddo
             enddo
          endif
       else
          do k=ndz_d,nfz_d
             do j=1,ny
                do i=1,nx
                   ! compute dissipation term: delta_k/delta_z=(D_{k+1/2}-D_{k-1/2})/dz
                   Krho(i,j,k) = Krho(i,j,k) +Drho(i,j,k-1) -Drho(i,j,k)
                   Krhou(i,j,k)= Krhou(i,j,k)+Drhou(i,j,k-1)-Drhou(i,j,k)
                   Krhov(i,j,k)= Krhov(i,j,k)+Drhov(i,j,k-1)-Drhov(i,j,k)
                   Krhow(i,j,k)= Krhow(i,j,k)+Drhow(i,j,k-1)-Drhow(i,j,k)
                   Krhoe(i,j,k)= Krhoe(i,j,k)+Drhoe(i,j,k-1)-Drhoe(i,j,k)
                enddo
             enddo
          enddo
       endif

    endif
       
  end subroutine filter_o10_shock

  !==============================================================================
  subroutine filter_o8_shock
  !==============================================================================
    !> Apply 8th-order filtering with Jameson-type shock capturing
    !> (Conservative version inspired from DNC order 9)
    !> modified version XG 2023
  !==============================================================================
    use mod_time
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: vc,cfll,mul
    real(wp), dimension(0:nx,0:ny,0:nz) :: Drho,Drhou,Drhov,Drhow,Drhoe
    ! ---------------------------------------------------------------------------
    real(wp) :: div2,om1,om2,om3
    real(wp) :: xk0,xk2,xk4,xk6,xk8
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

             !uvar(i,j,k,2)=psens(i,j,k)
          enddo
       enddo
    enddo

    ! Filtering in i-direction [compute Drho at i+1/2]
    ! ========================

    if (is_boundary(1,1)) then
       ! no dissipation at i=1 -> initialize increment array
       do j=1,ny
           Krho(1,j,:)=0.0_wp
          Krhou(1,j,:)=0.0_wp
          Krhov(1,j,:)=0.0_wp
          Krhow(1,j,:)=0.0_wp
          Krhoe(1,j,:)=0.0_wp
       enddo

       ! 1st-order centered dissipation
       ! ------------------------------
       if (is_bc_wall(1,1)) then
          i=1
          do k=1,nz
             do j=1,ny
                xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)
                !uvar(i,j,k,1)=xk2

                Drho(i,j,k) = (rho_n(i+1,j,k)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i+1,j,k)- rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
       else
          i=1
          do k=1,nz
             do j=1,ny
                xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)
                !uvar(i,j,k,1)=xk2

                Drho(i,j,k) = (rho_n(i+1,j,k)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i+1,j,k)- rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
       endif

       ! 3rd-order centered dissipation
       ! ------------------------------
       i=2
       do k=1,nz
          do j=1,ny
             xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))
             xk4=max(0.0_wp,xnu4+xk0)
             !uvar(i,j,k,1)=xk4/xnu4

             Drho(i,j,k) = (rho_n(i+1,j,k)-rho_n(i,j,k))*xk0      &
                         + (d4(1)*(rho_n(i+1,j,k)-rho_n(i  ,j,k)) &
                           +d4(2)*(rho_n(i+2,j,k)-rho_n(i-1,j,k)))*xk4
             Drhou(i,j,k)= (rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk0      &
                         + (d4(1)*(rhou_n(i+1,j,k)-rhou_n(i  ,j,k)) &
                           +d4(2)*(rhou_n(i+2,j,k)-rhou_n(i-1,j,k)))*xk4
             Drhov(i,j,k)= (rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk0      &
                         + (d4(1)*(rhov_n(i+1,j,k)-rhov_n(i  ,j,k)) &
                           +d4(2)*(rhov_n(i+2,j,k)-rhov_n(i-1,j,k)))*xk4
             Drhow(i,j,k)= (rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk0      &
                         + (d4(1)*(rhow_n(i+1,j,k)-rhow_n(i  ,j,k)) &
                           +d4(2)*(rhow_n(i+2,j,k)-rhow_n(i-1,j,k)))*xk4
             Drhoe(i,j,k)= (rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk0      &
                         + (d4(1)*(rhoe_n(i+1,j,k)-rhoe_n(i  ,j,k)) &
                           +d4(2)*(rhoe_n(i+2,j,k)-rhoe_n(i-1,j,k)))*xk4
          enddo
       enddo

       ! 5th-order centered dissipation
       ! ------------------------------
       i=3
       do k=1,nz
          do j=1,ny
             xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))
             !xk6=max(0.0_wp,xnu6+xk0/5.0_wp)
             xk6=max(0.0_wp,xnu6+xk0)
             !uvar(i,j,k,1)=xk6/xnu6

             Drho(i,j,k) = (rho_n(i+1,j,k)-rho_n(i,j,k))*xk0      &
                         + (d6(1)*(rho_n(i+1,j,k)-rho_n(i  ,j,k)) &
                           +d6(2)*(rho_n(i+2,j,k)-rho_n(i-1,j,k)) &
                           +d6(3)*(rho_n(i+3,j,k)-rho_n(i-2,j,k)))*xk6
             Drhou(i,j,k)= (rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk0      &
                         + (d6(1)*(rhou_n(i+1,j,k)-rhou_n(i  ,j,k)) &
                           +d6(2)*(rhou_n(i+2,j,k)-rhou_n(i-1,j,k)) &
                           +d6(3)*(rhou_n(i+3,j,k)-rhou_n(i-2,j,k)))*xk6
             Drhov(i,j,k)= (rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk0      &
                         + (d6(1)*(rhov_n(i+1,j,k)-rhov_n(i  ,j,k)) &
                           +d6(2)*(rhov_n(i+2,j,k)-rhov_n(i-1,j,k)) &
                           +d6(3)*(rhov_n(i+3,j,k)-rhov_n(i-2,j,k)))*xk6
             Drhow(i,j,k)= (rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk0      &
                         + (d6(1)*(rhow_n(i+1,j,k)-rhow_n(i  ,j,k)) &
                           +d6(2)*(rhow_n(i+2,j,k)-rhow_n(i-1,j,k)) &
                           +d6(3)*(rhow_n(i+3,j,k)-rhow_n(i-2,j,k)))*xk6
             Drhoe(i,j,k)= (rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk0      &
                         + (d6(1)*(rhoe_n(i+1,j,k)-rhoe_n(i  ,j,k)) &
                           +d6(2)*(rhoe_n(i+2,j,k)-rhoe_n(i-1,j,k)) &
                           +d6(3)*(rhoe_n(i+3,j,k)-rhoe_n(i-2,j,k)))*xk6
          enddo
       enddo
    endif

    ! Interior points: 7th-order centered dissipation
    ! -----------------------------------------------
    do k=1,nz
       do j=1,ny
          do i=ndx_di,nfx_di
             xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))
             xk8=max(0.0_wp,xnu8+3.0_wp/70.0_wp*xk0)
             !uvar(i,j,k,1)=xk8/xnu8
             Drho(i,j,k) = (rho_n(i+1,j,k)-rho_n(i,j,k))*xk0      &
                         + (d8(1)*(rho_n(i+1,j,k)-rho_n(i  ,j,k)) &
                           +d8(2)*(rho_n(i+2,j,k)-rho_n(i-1,j,k)) &
                           +d8(3)*(rho_n(i+3,j,k)-rho_n(i-2,j,k)) &
                           +d8(4)*(rho_n(i+4,j,k)-rho_n(i-3,j,k)))*xk8
             Drhou(i,j,k)= (rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk0      &
                         + (d8(1)*(rhou_n(i+1,j,k)-rhou_n(i  ,j,k)) &
                           +d8(2)*(rhou_n(i+2,j,k)-rhou_n(i-1,j,k)) &
                           +d8(3)*(rhou_n(i+3,j,k)-rhou_n(i-2,j,k)) &
                           +d8(4)*(rhou_n(i+4,j,k)-rhou_n(i-3,j,k)))*xk8
             Drhov(i,j,k)= (rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk0      &
                         + (d8(1)*(rhov_n(i+1,j,k)-rhov_n(i  ,j,k)) &
                           +d8(2)*(rhov_n(i+2,j,k)-rhov_n(i-1,j,k)) &
                           +d8(3)*(rhov_n(i+3,j,k)-rhov_n(i-2,j,k)) &
                           +d8(4)*(rhov_n(i+4,j,k)-rhov_n(i-3,j,k)))*xk8
             Drhow(i,j,k)= (rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk0      &
                         + (d8(1)*(rhow_n(i+1,j,k)-rhow_n(i  ,j,k)) &
                           +d8(2)*(rhow_n(i+2,j,k)-rhow_n(i-1,j,k)) &
                           +d8(3)*(rhow_n(i+3,j,k)-rhow_n(i-2,j,k)) &
                           +d8(4)*(rhow_n(i+4,j,k)-rhow_n(i-3,j,k)))*xk8
             Drhoe(i,j,k)= (rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk0      &
                         + (d8(1)*(rhoe_n(i+1,j,k)-rhoe_n(i  ,j,k)) &
                           +d8(2)*(rhoe_n(i+2,j,k)-rhoe_n(i-1,j,k)) &
                           +d8(3)*(rhoe_n(i+3,j,k)-rhoe_n(i-2,j,k)) &
                           +d8(4)*(rhoe_n(i+4,j,k)-rhoe_n(i-3,j,k)))*xk8
          enddo
       enddo
    enddo

    if (is_boundary(1,2)) then
       ! 5th-order centered dissipation
       ! ------------------------------
       i=nx-3
       do k=1,nz
          do j=1,ny
             xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))
             !xk6=max(0.0_wp,xnu6+xk0/5.0_wp)
             xk6=max(0.0_wp,xnu6+xk0)
             !uvar(i,j,k,1)=xk6/xnu6

             Drho(i,j,k) = (rho_n(i+1,j,k)-rho_n(i,j,k))*xk0      &
                         + (d6(1)*(rho_n(i+1,j,k)-rho_n(i  ,j,k)) &
                           +d6(2)*(rho_n(i+2,j,k)-rho_n(i-1,j,k)) &
                           +d6(3)*(rho_n(i+3,j,k)-rho_n(i-2,j,k)))*xk6
             Drhou(i,j,k)= (rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk0      &
                         + (d6(1)*(rhou_n(i+1,j,k)-rhou_n(i  ,j,k)) &
                           +d6(2)*(rhou_n(i+2,j,k)-rhou_n(i-1,j,k)) &
                           +d6(3)*(rhou_n(i+3,j,k)-rhou_n(i-2,j,k)))*xk6
             Drhov(i,j,k)= (rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk0      &
                         + (d6(1)*(rhov_n(i+1,j,k)-rhov_n(i  ,j,k)) &
                           +d6(2)*(rhov_n(i+2,j,k)-rhov_n(i-1,j,k)) &
                           +d6(3)*(rhov_n(i+3,j,k)-rhov_n(i-2,j,k)))*xk6
             Drhow(i,j,k)= (rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk0      &
                         + (d6(1)*(rhow_n(i+1,j,k)-rhow_n(i  ,j,k)) &
                           +d6(2)*(rhow_n(i+2,j,k)-rhow_n(i-1,j,k)) &
                           +d6(3)*(rhow_n(i+3,j,k)-rhow_n(i-2,j,k)))*xk6
             Drhoe(i,j,k)= (rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk0      &
                         + (d6(1)*(rhoe_n(i+1,j,k)-rhoe_n(i  ,j,k)) &
                           +d6(2)*(rhoe_n(i+2,j,k)-rhoe_n(i-1,j,k)) &
                           +d6(3)*(rhoe_n(i+3,j,k)-rhoe_n(i-2,j,k)))*xk6
          enddo
       enddo

       ! 3rd-order centered dissipation
       ! ------------------------------
       i=nx-2
       do k=1,nz
          do j=1,ny
             xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))
             xk4=max(0.0_wp,xnu4+xk0)
             !uvar(i,j,k,1)=xk4/xnu4

             Drho(i,j,k) = (rho_n(i+1,j,k)-rho_n(i,j,k))*xk0      &
                         + (d4(1)*(rho_n(i+1,j,k)-rho_n(i  ,j,k)) &
                           +d4(2)*(rho_n(i+2,j,k)-rho_n(i-1,j,k)))*xk4
             Drhou(i,j,k)= (rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk0      &
                         + (d4(1)*(rhou_n(i+1,j,k)-rhou_n(i  ,j,k)) &
                           +d4(2)*(rhou_n(i+2,j,k)-rhou_n(i-1,j,k)))*xk4
             Drhov(i,j,k)= (rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk0      &
                         + (d4(1)*(rhov_n(i+1,j,k)-rhov_n(i  ,j,k)) &
                           +d4(2)*(rhov_n(i+2,j,k)-rhov_n(i-1,j,k)))*xk4
             Drhow(i,j,k)= (rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk0      &
                         + (d4(1)*(rhow_n(i+1,j,k)-rhow_n(i  ,j,k)) &
                           +d4(2)*(rhow_n(i+2,j,k)-rhow_n(i-1,j,k)))*xk4
             Drhoe(i,j,k)= (rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk0      &
                         + (d4(1)*(rhoe_n(i+1,j,k)-rhoe_n(i  ,j,k)) &
                           +d4(2)*(rhoe_n(i+2,j,k)-rhoe_n(i-1,j,k)))*xk4

          enddo
       enddo

       ! 1st-order centered dissipation
       ! ------------------------------
       if (is_bc_wall(1,2)) then
          i=nx-1
          do k=1,nz
             do j=1,ny
                xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)
                !uvar(i,j,k,1)=xk2

                Drho(i,j,k) = (rho_n(i+1,j,k)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i+1,j,k)-rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
       else
          i=nx-1
          do k=1,nz
             do j=1,ny
                xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)
                !uvar(i,j,k,1)=xk2

                Drho(i,j,k) = (rho_n(i+1,j,k)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i+1,j,k)-rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
       endif

       ! no dissipation at i=nx -> initialize increment array
       do j=1,ny
           Krho(nx,j,:)=0.0_wp
          Krhou(nx,j,:)=0.0_wp
          Krhov(nx,j,:)=0.0_wp
          Krhow(nx,j,:)=0.0_wp
          Krhoe(nx,j,:)=0.0_wp
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
                   !mul=min(1.0_wp,cfl_i(i,j,k)/chi)

                   ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dx
                   Krho(i,j,k) = mul*(Drho(i-1,j,k) -Drho(i,j,k))
                   Krhou(i,j,k)= mul*(Drhou(i-1,j,k)-Drhou(i,j,k))
                   Krhov(i,j,k)= mul*(Drhov(i-1,j,k)-Drhov(i,j,k))
                   Krhow(i,j,k)= mul*(Drhow(i-1,j,k)-Drhow(i,j,k))
                   Krhoe(i,j,k)= mul*(Drhoe(i-1,j,k)-Drhoe(i,j,k))
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
                   !mul=min(1.0_wp,cfl_i(i,j,k)/chi)

                   ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dx
                   Krho(i,j,k) = mul*(Drho(i-1,j,k) -Drho(i,j,k))
                   Krhou(i,j,k)= mul*(Drhou(i-1,j,k)-Drhou(i,j,k))
                   Krhov(i,j,k)= mul*(Drhov(i-1,j,k)-Drhov(i,j,k))
                   Krhow(i,j,k)= mul*(Drhow(i-1,j,k)-Drhow(i,j,k))
                   Krhoe(i,j,k)= mul*(Drhoe(i-1,j,k)-Drhoe(i,j,k))
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
                   ! If call to local_cfl
                   !mul=min(1.0_wp,cfl_i(i,j,k)/chi)

                   ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dx
                   Krho(i,j,k) = mul*(Drho(i-1,j,k) -Drho(i,j,k))
                   Krhou(i,j,k)= mul*(Drhou(i-1,j,k)-Drhou(i,j,k))
                   Krhov(i,j,k)= mul*(Drhov(i-1,j,k)-Drhov(i,j,k))
                   Krhow(i,j,k)= mul*(Drhow(i-1,j,k)-Drhow(i,j,k))
                   Krhoe(i,j,k)= mul*(Drhoe(i-1,j,k)-Drhoe(i,j,k))
                enddo
             enddo
          enddo
       endif
    else
       do k=1,nz
          do j=1,ny
             do i=ndx_d,nfx_d
                ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dx
                Krho(i,j,k) = Drho(i-1,j,k) -Drho(i,j,k)
                Krhou(i,j,k)= Drhou(i-1,j,k)-Drhou(i,j,k)
                Krhov(i,j,k)= Drhov(i-1,j,k)-Drhov(i,j,k)
                Krhow(i,j,k)= Drhow(i-1,j,k)-Drhow(i,j,k)
                Krhoe(i,j,k)= Drhoe(i-1,j,k)-Drhoe(i,j,k)
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

    ! Filtering in y-direction [compute Drho at j+1/2]
    ! ========================

    if (is_boundary(2,1)) then
       ! 1st-order centered dissipation
       ! ------------------------------
       if (is_bc_wall(2,1)) then
          j=1
          do k=1,nz
             do i=1,nx
                xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)
                !uvar(i,j,k,2)=xk2

                Drho(i,j,k) = (rho_n(i,j+1,k)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i,j+1,k)-rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
       else
          j=1
          do k=1,nz
             do i=1,nx
                xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)
                !uvar(i,j,k,2)=xk2

                Drho(i,j,k) = (rho_n(i,j+1,k)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i,j+1,k)-rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
       endif

       ! 3rd-order centered dissipation
       ! ------------------------------
       j=2
       do k=1,nz
          do i=1,nx
             xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))
             xk4=max(0.0_wp,xnu4+xk0)
             !uvar(i,j,k,2)=xk4/xnu4

             Drho(i,j,k) = (rho_n(i,j+1,k)-rho_n(i,j,k))*xk0      &
                         + (d4(1)*(rho_n(i,j+1,k)-rho_n(i,j  ,k)) &
                           +d4(2)*(rho_n(i,j+2,k)-rho_n(i,j-1,k)))*xk4
             Drhou(i,j,k)= (rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk0      &
                         + (d4(1)*(rhou_n(i,j+1,k)-rhou_n(i,j  ,k)) &
                           +d4(2)*(rhou_n(i,j+2,k)-rhou_n(i,j-1,k)))*xk4
             Drhov(i,j,k)= (rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk0      &
                         + (d4(1)*(rhov_n(i,j+1,k)-rhov_n(i,j  ,k)) &
                           +d4(2)*(rhov_n(i,j+2,k)-rhov_n(i,j-1,k)))*xk4
             Drhow(i,j,k)= (rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk0      &
                         + (d4(1)*(rhow_n(i,j+1,k)-rhow_n(i,j  ,k)) &
                           +d4(2)*(rhow_n(i,j+2,k)-rhow_n(i,j-1,k)))*xk4
             Drhoe(i,j,k)= (rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk0      &
                         + (d4(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j  ,k)) &
                           +d4(2)*(rhoe_n(i,j+2,k)-rhoe_n(i,j-1,k)))*xk4
          enddo
       enddo

       ! 5th-order centered dissipation
       ! ------------------------------
       j=3
       do k=1,nz
          do i=1,nx
             xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))
             !xk6=max(0.0_wp,xnu6+xk0/5.0_wp)
             xk6=max(0.0_wp,xnu6+xk0)
             !uvar(i,j,k,2)=xk6/xnu6

             Drho(i,j,k) = (rho_n(i,j+1,k)-rho_n(i,j,k))*xk0      &
                         + (d6(1)*(rho_n(i,j+1,k)-rho_n(i,j  ,k)) &
                           +d6(2)*(rho_n(i,j+2,k)-rho_n(i,j-1,k)) &
                           +d6(3)*(rho_n(i,j+3,k)-rho_n(i,j-2,k)))*xk6
             Drhou(i,j,k)= (rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk0      &
                         + (d6(1)*(rhou_n(i,j+1,k)-rhou_n(i,j  ,k)) &
                           +d6(2)*(rhou_n(i,j+2,k)-rhou_n(i,j-1,k)) &
                           +d6(3)*(rhou_n(i,j+3,k)-rhou_n(i,j-2,k)))*xk6
             Drhov(i,j,k)= (rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk0      &
                         + (d6(1)*(rhov_n(i,j+1,k)-rhov_n(i,j  ,k)) &
                           +d6(2)*(rhov_n(i,j+2,k)-rhov_n(i,j-1,k)) &
                           +d6(3)*(rhov_n(i,j+3,k)-rhov_n(i,j-2,k)))*xk6
             Drhow(i,j,k)= (rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk0      &
                         + (d6(1)*(rhow_n(i,j+1,k)-rhow_n(i,j  ,k)) &
                           +d6(2)*(rhow_n(i,j+2,k)-rhow_n(i,j-1,k)) &
                           +d6(3)*(rhow_n(i,j+3,k)-rhow_n(i,j-2,k)))*xk6
             Drhoe(i,j,k)= (rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk0      &
                         + (d6(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j  ,k)) &
                           +d6(2)*(rhoe_n(i,j+2,k)-rhoe_n(i,j-1,k)) &
                           +d6(3)*(rhoe_n(i,j+3,k)-rhoe_n(i,j-2,k)))*xk6
          enddo
       enddo
    endif

    ! Interior points: 7th-order centered dissipation
    ! -----------------------------------------------
    do k=1,nz
       do j=ndy_di,nfy_di
          do i=1,nx
             xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))
             xk8=max(0.0_wp,xnu8+3.0_wp/70.0_wp*xk0)
             !uvar(i,j,k,2)=xk8/xnu8

             Drho(i,j,k) = (rho_n(i,j+1,k)-rho_n(i,j,k))*xk0      &
                         + (d8(1)*(rho_n(i,j+1,k)-rho_n(i,j  ,k)) &
                           +d8(2)*(rho_n(i,j+2,k)-rho_n(i,j-1,k)) &
                           +d8(3)*(rho_n(i,j+3,k)-rho_n(i,j-2,k)) &
                           +d8(4)*(rho_n(i,j+4,k)-rho_n(i,j-3,k)))*xk8
             Drhou(i,j,k)= (rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk0      &
                         + (d8(1)*(rhou_n(i,j+1,k)-rhou_n(i,j  ,k)) &
                           +d8(2)*(rhou_n(i,j+2,k)-rhou_n(i,j-1,k)) &
                           +d8(3)*(rhou_n(i,j+3,k)-rhou_n(i,j-2,k)) &
                           +d8(4)*(rhou_n(i,j+4,k)-rhou_n(i,j-3,k)))*xk8
             Drhov(i,j,k)= (rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk0      &
                         + (d8(1)*(rhov_n(i,j+1,k)-rhov_n(i,j  ,k)) &
                           +d8(2)*(rhov_n(i,j+2,k)-rhov_n(i,j-1,k)) &
                           +d8(3)*(rhov_n(i,j+3,k)-rhov_n(i,j-2,k)) &
                           +d8(4)*(rhov_n(i,j+4,k)-rhov_n(i,j-3,k)))*xk8
             Drhow(i,j,k)= (rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk0      &
                         + (d8(1)*(rhow_n(i,j+1,k)-rhow_n(i,j  ,k)) &
                           +d8(2)*(rhow_n(i,j+2,k)-rhow_n(i,j-1,k)) &
                           +d8(3)*(rhow_n(i,j+3,k)-rhow_n(i,j-2,k)) &
                           +d8(4)*(rhow_n(i,j+4,k)-rhow_n(i,j-3,k)))*xk8
             Drhoe(i,j,k)= (rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk0      &
                         + (d8(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j  ,k)) &
                           +d8(2)*(rhoe_n(i,j+2,k)-rhoe_n(i,j-1,k)) &
                           +d8(3)*(rhoe_n(i,j+3,k)-rhoe_n(i,j-2,k)) &
                           +d8(4)*(rhoe_n(i,j+4,k)-rhoe_n(i,j-3,k)))*xk8
          enddo
       enddo
    enddo

    if (is_boundary(2,2)) then
       ! 5th-order centered dissipation
       ! ------------------------------
       j=ny-3
       do k=1,nz
          do i=1,nx
             xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))
             !xk6=max(0.0_wp,xnu6+xk0/5.0_wp)
             xk6=max(0.0_wp,xnu6+xk0)
             !uvar(i,j,k,2)=xk6/xnu6

             Drho(i,j,k) = (rho_n(i,j+1,k)-rho_n(i,j,k))*xk0      &
                         + (d6(1)*(rho_n(i,j+1,k)-rho_n(i,j  ,k)) &
                           +d6(2)*(rho_n(i,j+2,k)-rho_n(i,j-1,k)) &
                           +d6(3)*(rho_n(i,j+3,k)-rho_n(i,j-2,k)))*xk6
             Drhou(i,j,k)= (rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk0      &
                         + (d6(1)*(rhou_n(i,j+1,k)-rhou_n(i,j  ,k)) &
                           +d6(2)*(rhou_n(i,j+2,k)-rhou_n(i,j-1,k)) &
                           +d6(3)*(rhou_n(i,j+3,k)-rhou_n(i,j-2,k)))*xk6
             Drhov(i,j,k)= (rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk0      &
                         + (d6(1)*(rhov_n(i,j+1,k)-rhov_n(i,j  ,k)) &
                           +d6(2)*(rhov_n(i,j+2,k)-rhov_n(i,j-1,k)) &
                           +d6(3)*(rhov_n(i,j+3,k)-rhov_n(i,j-2,k)))*xk6
             Drhow(i,j,k)= (rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk0      &
                         + (d6(1)*(rhow_n(i,j+1,k)-rhow_n(i,j  ,k)) &
                           +d6(2)*(rhow_n(i,j+2,k)-rhow_n(i,j-1,k)) &
                           +d6(3)*(rhow_n(i,j+3,k)-rhow_n(i,j-2,k)))*xk6
             Drhoe(i,j,k)= (rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk0      &
                         + (d6(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j  ,k)) &
                           +d6(2)*(rhoe_n(i,j+2,k)-rhoe_n(i,j-1,k)) &
                           +d6(3)*(rhoe_n(i,j+3,k)-rhoe_n(i,j-2,k)))*xk6
          enddo
       enddo

       ! 3rd-order centered dissipation
       ! ------------------------------
       j=ny-2
       do k=1,nz
          do i=1,nx
             xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))
             xk4=max(0.0_wp,xnu4+xk0)
             !uvar(i,j,k,2)=xk4/xnu4

             Drho(i,j,k) = (rho_n(i,j+1,k)-rho_n(i,j,k))*xk0      &
                         + (d4(1)*(rho_n(i,j+1,k)-rho_n(i,j  ,k)) &
                           +d4(2)*(rho_n(i,j+2,k)-rho_n(i,j-1,k)))*xk4
             Drhou(i,j,k)= (rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk0      &
                         + (d4(1)*(rhou_n(i,j+1,k)-rhou_n(i,j  ,k)) &
                           +d4(2)*(rhou_n(i,j+2,k)-rhou_n(i,j-1,k)))*xk4
             Drhov(i,j,k)= (rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk0      &
                         + (d4(1)*(rhov_n(i,j+1,k)-rhov_n(i,j  ,k)) &
                           +d4(2)*(rhov_n(i,j+2,k)-rhov_n(i,j-1,k)))*xk4
             Drhow(i,j,k)= (rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk0      &
                         + (d4(1)*(rhow_n(i,j+1,k)-rhow_n(i,j  ,k)) &
                           +d4(2)*(rhow_n(i,j+2,k)-rhow_n(i,j-1,k)))*xk4
             Drhoe(i,j,k)= (rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk0      &
                         + (d4(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j  ,k)) &
                           +d4(2)*(rhoe_n(i,j+2,k)-rhoe_n(i,j-1,k)))*xk4
          enddo
       enddo

       ! 1st-order centered dissipation
       ! ------------------------------
       if (is_bc_wall(2,2)) then
          j=ny-1
          do k=1,nz
             do i=1,nx
                xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)
                !uvar(i,j,k,2)=xk2

                Drho(i,j,k) = (rho_n(i,j+1,k)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i,j+1,k)-rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
       else
          j=ny-1
          do k=1,nz
             do i=1,nx
                xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)
                !uvar(i,j,k,2)=xk2

                Drho(i,j,k) = (rho_n(i,j+1,k)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i,j+1,k)-rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
       endif
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

                   ! compute dissipation term: delta_i/delta_y=(D_{j+1/2}-D_{j-1/2})/dy
                   Krho(i,j,k) = Krho(i,j,k) +mul*(Drho(i,j-1,k) -Drho(i,j,k))
                   Krhou(i,j,k)= Krhou(i,j,k)+mul*(Drhou(i,j-1,k)-Drhou(i,j,k))
                   Krhov(i,j,k)= Krhov(i,j,k)+mul*(Drhov(i,j-1,k)-Drhov(i,j,k))
                   Krhow(i,j,k)= Krhow(i,j,k)+mul*(Drhow(i,j-1,k)-Drhow(i,j,k))
                   Krhoe(i,j,k)= Krhoe(i,j,k)+mul*(Drhoe(i,j-1,k)-Drhoe(i,j,k))
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
                   ! ! If call to local_cfl
                   !mul=min(1.0_wp,cfl_j(i,j,k)/chi)

                   ! compute dissipation term: delta_i/delta_y=(D_{j+1/2}-D_{j-1/2})/dy
                   Krho(i,j,k) = Krho(i,j,k) +mul*(Drho(i,j-1,k) -Drho(i,j,k))
                   Krhou(i,j,k)= Krhou(i,j,k)+mul*(Drhou(i,j-1,k)-Drhou(i,j,k))
                   Krhov(i,j,k)= Krhov(i,j,k)+mul*(Drhov(i,j-1,k)-Drhov(i,j,k))
                   Krhow(i,j,k)= Krhow(i,j,k)+mul*(Drhow(i,j-1,k)-Drhow(i,j,k))
                   Krhoe(i,j,k)= Krhoe(i,j,k)+mul*(Drhoe(i,j-1,k)-Drhoe(i,j,k))
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
                   !mul=min(1.0_wp,cfl_j(i,j,k)/chi)

                   ! compute dissipation term: delta_i/delta_y=(D_{j+1/2}-D_{j-1/2})/dy
                   Krho(i,j,k) = Krho(i,j,k) +mul*(Drho(i,j-1,k) -Drho(i,j,k))
                   Krhou(i,j,k)= Krhou(i,j,k)+mul*(Drhou(i,j-1,k)-Drhou(i,j,k))
                   Krhov(i,j,k)= Krhov(i,j,k)+mul*(Drhov(i,j-1,k)-Drhov(i,j,k))
                   Krhow(i,j,k)= Krhow(i,j,k)+mul*(Drhow(i,j-1,k)-Drhow(i,j,k))
                   Krhoe(i,j,k)= Krhoe(i,j,k)+mul*(Drhoe(i,j-1,k)-Drhoe(i,j,k))
                enddo
             enddo
          enddo
       endif
    else
       do k=1,nz
          do j=ndy_d,nfy_d
             do i=1,nx
                ! compute dissipation term: delta_j/delta_y=(D_{j+1/2}-D_{j-1/2})/dy
                Krho(i,j,k) = Krho(i,j,k) +Drho(i,j-1,k) -Drho(i,j,k)
                Krhou(i,j,k)= Krhou(i,j,k)+Drhou(i,j-1,k)-Drhou(i,j,k)
                Krhov(i,j,k)= Krhov(i,j,k)+Drhov(i,j-1,k)-Drhov(i,j,k)
                Krhow(i,j,k)= Krhow(i,j,k)+Drhow(i,j-1,k)-Drhow(i,j,k)
                Krhoe(i,j,k)= Krhoe(i,j,k)+Drhoe(i,j-1,k)-Drhoe(i,j,k)
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

       ! Filtering in k-direction [compute Drho at k+1/2]
       ! ========================
       if (is_boundary(3,1)) then
          ! 1st-order centered dissipation
          ! ------------------------------
          k=1
          do j=1,ny
             do i=1,nx
                xk0=xnu0*max(psens(i,j,k+1),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)

                Drho(i,j,k) = (rho_n(i,j,k+1)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i,j,k+1)-rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i,j,k+1)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i,j,k+1)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i,j,k+1)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i,j,k+1)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k))*xk2
             enddo
          enddo

          ! 3rd-order centered dissipation
          ! ------------------------------
          k=2
          do j=1,ny
             do i=1,nx
                xk0=xnu0*max(psens(i,j,k+1),psens(i,j,k))
                xk4=max(0.0_wp,xnu4+xk0)

                Drho(i,j,k) = (rho_n(i,j,k+1)-rho_n(i,j,k))*xk0      &
                            + (d4(1)*(rho_n(i,j,k+1)-rho_n(i,j,k  )) &
                              +d4(2)*(rho_n(i,j,k+2)-rho_n(i,j,k-1)))*xk4
                Drhou(i,j,k)= (rhou_n(i,j,k+1)-rhou_n(i,j,k))*xk0      &
                            + (d4(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k  )) &
                              +d4(2)*(rhou_n(i,j,k+2)-rhou_n(i,j,k-1)))*xk4
                Drhov(i,j,k)= (rhov_n(i,j,k+1)-rhov_n(i,j,k))*xk0      &
                            + (d4(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k  )) &
                              +d4(2)*(rhov_n(i,j,k+2)-rhov_n(i,j,k-1)))*xk4
                Drhow(i,j,k)= (rhow_n(i,j,k+1)-rhow_n(i,j,k))*xk0      &
                            + (d4(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k  )) &
                              +d4(2)*(rhow_n(i,j,k+2)-rhow_n(i,j,k-1)))*xk4
                Drhoe(i,j,k)= (rhoe_n(i,j,k+1)-rhoe_n(i,j,k))*xk0      &
                            + (d4(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k  )) &
                              +d4(2)*(rhoe_n(i,j,k+2)-rhoe_n(i,j,k-1)))*xk4
             enddo
          enddo

          ! 5th-order centered dissipation
          ! ------------------------------
          k=3
          do j=1,ny
             do i=1,nx
                xk0=xnu0*max(psens(i,j,k+1),psens(i,j,k))
                !xk6=max(0.0_wp,xnu6+xk0/5.0_wp)
                xk6=max(0.0_wp,xnu6+xk0)

                Drho(i,j,k) = (rho_n(i,j,k+1)-rho_n(i,j,k))*xk0      &
                            + (d6(1)*(rho_n(i,j,k+1)-rho_n(i,j,k  )) &
                              +d6(2)*(rho_n(i,j,k+2)-rho_n(i,j,k-1)) &
                              +d6(3)*(rho_n(i,j,k+3)-rho_n(i,j,k-2)))*xk6
                Drhou(i,j,k)= (rhou_n(i,j,k+1)-rhou_n(i,j,k))*xk0      &
                            + (d6(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k  )) &
                              +d6(2)*(rhou_n(i,j,k+2)-rhou_n(i,j,k-1)) &
                              +d6(3)*(rhou_n(i,j,k+3)-rhou_n(i,j,k-2)))*xk6
                Drhov(i,j,k)= (rhov_n(i,j,k+1)-rhov_n(i,j,k))*xk0      &
                            + (d6(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k  )) &
                              +d6(2)*(rhov_n(i,j,k+2)-rhov_n(i,j,k-1)) &
                              +d6(3)*(rhov_n(i,j,k+3)-rhov_n(i,j,k-2)))*xk6
                Drhow(i,j,k)= (rhow_n(i,j,k+1)-rhow_n(i,j,k))*xk0      &
                            + (d6(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k  )) &
                              +d6(2)*(rhow_n(i,j,k+2)-rhow_n(i,j,k-1)) &
                              +d6(3)*(rhow_n(i,j,k+3)-rhow_n(i,j,k-2)))*xk6
                Drhoe(i,j,k)= (rhoe_n(i,j,k+1)-rhoe_n(i,j,k))*xk0      &
                            + (d6(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k  )) &
                              +d6(2)*(rhoe_n(i,j,k+2)-rhoe_n(i,j,k-1)) &
                              +d6(3)*(rhoe_n(i,j,k+3)-rhoe_n(i,j,k-2)))*xk6
             enddo
          enddo
       endif

       ! Interior points: 7th-order centered dissipation
       ! -----------------------------------------------
       do k=ndz_di,nfz_di
          do j=1,ny
             do i=1,nx
                xk0=xnu0*max(psens(i,j,k+1),psens(i,j,k))
                xk8=max(0.0_wp,xnu8+3.0_wp/70.0_wp*xk0)

                Drho(i,j,k) = (rho_n(i,j,k+1)-rho_n(i,j,k))*xk0      &
                            + (d8(1)*(rho_n(i,j,k+1)-rho_n(i,j,k  )) &
                              +d8(2)*(rho_n(i,j,k+2)-rho_n(i,j,k-1)) &
                              +d8(3)*(rho_n(i,j,k+3)-rho_n(i,j,k-2)) &
                              +d8(4)*(rho_n(i,j,k+4)-rho_n(i,j,k-3)))*xk8
                Drhou(i,j,k)= (rhou_n(i,j,k+1)-rhou_n(i,j,k))*xk0      &
                            + (d8(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k  )) &
                              +d8(2)*(rhou_n(i,j,k+2)-rhou_n(i,j,k-1)) &
                              +d8(3)*(rhou_n(i,j,k+3)-rhou_n(i,j,k-2)) &
                              +d8(4)*(rhou_n(i,j,k+4)-rhou_n(i,j,k-3)))*xk8
                Drhov(i,j,k)= (rhov_n(i,j,k+1)-rhov_n(i,j,k))*xk0      &
                            + (d8(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k  )) &
                              +d8(2)*(rhov_n(i,j,k+2)-rhov_n(i,j,k-1)) &
                              +d8(3)*(rhov_n(i,j,k+3)-rhov_n(i,j,k-2)) &
                              +d8(4)*(rhov_n(i,j,k+4)-rhov_n(i,j,k-3)))*xk8
                Drhow(i,j,k)= (rhow_n(i,j,k+1)-rhow_n(i,j,k))*xk0      &
                            + (d8(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k  )) &
                              +d8(2)*(rhow_n(i,j,k+2)-rhow_n(i,j,k-1)) &
                              +d8(3)*(rhow_n(i,j,k+3)-rhow_n(i,j,k-2)) &
                              +d8(4)*(rhow_n(i,j,k+4)-rhow_n(i,j,k-3)))*xk8
                Drhoe(i,j,k)= (rhoe_n(i,j,k+1)-rhoe_n(i,j,k))*xk0      &
                            + (d8(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k  )) &
                              +d8(2)*(rhoe_n(i,j,k+2)-rhoe_n(i,j,k-1)) &
                              +d8(3)*(rhoe_n(i,j,k+3)-rhoe_n(i,j,k-2)) &
                              +d8(4)*(rhoe_n(i,j,k+4)-rhoe_n(i,j,k-3)))*xk8
             enddo
          enddo
       enddo

       if (is_boundary(3,2)) then
          ! 5th-order centered dissipation
          ! ------------------------------
          k=nz-3
          do j=1,ny
             do i=1,nx
                xk0=xnu0*max(psens(i,j,k+1),psens(i,j,k))
                !xk6=max(0.0_wp,xnu6+xk0/5.0_wp)
                xk6=max(0.0_wp,xnu6+xk0)

                Drho(i,j,k) = (rho_n(i,j,k+1)-rho_n(i,j,k))*xk0      &
                            + (d6(1)*(rho_n(i,j,k+1)-rho_n(i,j,k  )) &
                              +d6(2)*(rho_n(i,j,k+2)-rho_n(i,j,k-1)) &
                              +d6(3)*(rho_n(i,j,k+3)-rho_n(i,j,k-2)))*xk6
                Drhou(i,j,k)= (rhou_n(i,j,k+1)-rhou_n(i,j,k))*xk0      &
                            + (d6(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k  )) &
                              +d6(2)*(rhou_n(i,j,k+2)-rhou_n(i,j,k-1)) &
                              +d6(3)*(rhou_n(i,j,k+3)-rhou_n(i,j,k-2)))*xk6
                Drhov(i,j,k)= (rhov_n(i,j,k+1)-rhov_n(i,j,k))*xk0      &
                            + (d6(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k  )) &
                              +d6(2)*(rhov_n(i,j,k+2)-rhov_n(i,j,k-1)) &
                              +d6(3)*(rhov_n(i,j,k+3)-rhov_n(i,j,k-2)))*xk6
                Drhow(i,j,k)= (rhow_n(i,j,k+1)-rhow_n(i,j,k))*xk0      &
                            + (d6(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k  )) &
                              +d6(2)*(rhow_n(i,j,k+2)-rhow_n(i,j,k-1)) &
                              +d6(3)*(rhow_n(i,j,k+3)-rhow_n(i,j,k-2)))*xk6
                Drhoe(i,j,k)= (rhoe_n(i,j,k+1)-rhoe_n(i,j,k))*xk0      &
                            + (d6(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k  )) &
                              +d6(2)*(rhoe_n(i,j,k+2)-rhoe_n(i,j,k-1)) &
                              +d6(3)*(rhoe_n(i,j,k+3)-rhoe_n(i,j,k-2)))*xk6
             enddo
          enddo

          ! 3rd-order centered dissipation
          ! ------------------------------
          k=nz-2
          do j=1,ny
             do i=1,nx
                xk0=xnu0*max(psens(i,j,k+1),psens(i,j,k))
                xk4=max(0.0_wp,xnu4+xk0)

                Drho(i,j,k) = (rho_n(i,j,k+1)-rho_n(i,j,k))*xk0      &
                            + (d4(1)*(rho_n(i,j,k+1)-rho_n(i,j,k  )) &
                              +d4(2)*(rho_n(i,j,k+2)-rho_n(i,j,k-1)))*xk4
                Drhou(i,j,k)= (rhou_n(i,j,k+1)-rhou_n(i,j,k))*xk0 &
                            + (d4(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k  )) &
                              +d4(2)*(rhou_n(i,j,k+2)-rhou_n(i,j,k-1)))*xk4
                Drhov(i,j,k)= (rhov_n(i,j,k+1)-rhov_n(i,j,k))*xk0      &
                            + (d4(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k  )) &
                              +d4(2)*(rhov_n(i,j,k+2)-rhov_n(i,j,k-1)))*xk4
                Drhow(i,j,k)= (rhow_n(i,j,k+1)-rhow_n(i,j,k))*xk0      &
                            + (d4(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k  )) &
                              +d4(2)*(rhow_n(i,j,k+2)-rhow_n(i,j,k-1)))*xk4
                Drhoe(i,j,k)= (rhoe_n(i,j,k+1)-rhoe_n(i,j,k))*xk0      &
                            + (d4(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k  )) &
                              +d4(2)*(rhoe_n(i,j,k+2)-rhoe_n(i,j,k-1)))*xk4
             enddo
          enddo

          ! 1st-order centered dissipation
          ! ------------------------------
          k=nz-1
          do j=1,ny
             do i=1,nx
                xk0=xnu0*max(psens(i,j,k+1),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)

                Drho(i,j,k) = (rho_n(i,j,k+1)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i,j,k+1)-rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i,j,k+1)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i,j,k+1)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i,j,k+1)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i,j,k+1)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
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
                      Krho(i,j,k) = Krho(i,j,k) +mul*(Drho(i,j,k-1) -Drho(i,j,k))
                      Krhou(i,j,k)= Krhou(i,j,k)+mul*(Drhou(i,j,k-1)-Drhou(i,j,k))
                      Krhov(i,j,k)= Krhov(i,j,k)+mul*(Drhov(i,j,k-1)-Drhov(i,j,k))
                      Krhow(i,j,k)= Krhow(i,j,k)+mul*(Drhow(i,j,k-1)-Drhow(i,j,k))
                      Krhoe(i,j,k)= Krhoe(i,j,k)+mul*(Drhoe(i,j,k-1)-Drhoe(i,j,k))
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
                      Krho(i,j,k) = Krho(i,j,k) +mul*(Drho(i,j,k-1) -Drho(i,j,k))
                      Krhou(i,j,k)= Krhou(i,j,k)+mul*(Drhou(i,j,k-1)-Drhou(i,j,k))
                      Krhov(i,j,k)= Krhov(i,j,k)+mul*(Drhov(i,j,k-1)-Drhov(i,j,k))
                      Krhow(i,j,k)= Krhow(i,j,k)+mul*(Drhow(i,j,k-1)-Drhow(i,j,k))
                      Krhoe(i,j,k)= Krhoe(i,j,k)+mul*(Drhoe(i,j,k-1)-Drhoe(i,j,k))
                   enddo
                enddo
             enddo
          endif
       else
          do k=ndz_d,nfz_d
             do j=1,ny
                do i=1,nx
                   ! compute dissipation term: delta_k/delta_z=(D_{k+1/2}-D_{k-1/2})/dz
                   Krho(i,j,k) = Krho(i,j,k) +Drho(i,j,k-1) -Drho(i,j,k)
                   Krhou(i,j,k)= Krhou(i,j,k)+Drhou(i,j,k-1)-Drhou(i,j,k)
                   Krhov(i,j,k)= Krhov(i,j,k)+Drhov(i,j,k-1)-Drhov(i,j,k)
                   Krhow(i,j,k)= Krhow(i,j,k)+Drhow(i,j,k-1)-Drhow(i,j,k)
                   Krhoe(i,j,k)= Krhoe(i,j,k)+Drhoe(i,j,k-1)-Drhoe(i,j,k)
                enddo
             enddo
          enddo
       endif

    endif

  end subroutine filter_o8_shock

  !==============================================================================
  subroutine filter_o6_shock
  !==============================================================================
    !> Apply 6th-order filtering with Jameson-type shock capturing
    !> (Conservative version inspired from DNC order 9)
    !> modified version XG 2023
  !==============================================================================
    use mod_time
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: vc,cfll,mul
    real(wp), dimension(0:nx,0:ny,0:nz) :: Drho,Drhou,Drhov,Drhow,Drhoe
    ! ---------------------------------------------------------------------------
    real(wp) :: div2,om1,om2,om3
    real(wp) :: xk0,xk2,xk4,xk6
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

             !uvar(i,j,k,2)=psens(i,j,k)
          enddo
       enddo
    enddo

    ! Filtering in i-direction [compute Drho at i+1/2]
    ! ========================

    if (is_boundary(1,1)) then
       ! no dissipation at i=1 -> initialize increment array
       do j=1,ny
           Krho(1,j,:)=0.0_wp
          Krhou(1,j,:)=0.0_wp
          Krhov(1,j,:)=0.0_wp
          Krhow(1,j,:)=0.0_wp
          Krhoe(1,j,:)=0.0_wp
       enddo

       ! 1st-order centered dissipation
       ! ------------------------------
       if (is_bc_wall(1,1)) then
          i=1
          do k=1,nz
             do j=1,ny
                xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)
                !uvar(i,j,k,1)=xk2

                Drho(i,j,k) = (rho_n(i+1,j,k)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i+1,j,k)- rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
       else
          i=1
          do k=1,nz
             do j=1,ny
                xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)
                !uvar(i,j,k,1)=xk2

                Drho(i,j,k) = (rho_n(i+1,j,k)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i+1,j,k)- rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
       endif

       ! 3rd-order centered dissipation
       ! ------------------------------
       i=2
       do k=1,nz
          do j=1,ny
             xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))
             xk4=max(0.0_wp,xnu4+xk0)
             !uvar(i,j,k,1)=xk4/xnu4

             Drho(i,j,k) = (rho_n(i+1,j,k)-rho_n(i,j,k))*xk0      &
                         + (d4(1)*(rho_n(i+1,j,k)-rho_n(i  ,j,k)) &
                           +d4(2)*(rho_n(i+2,j,k)-rho_n(i-1,j,k)))*xk4
             Drhou(i,j,k)= (rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk0      &
                         + (d4(1)*(rhou_n(i+1,j,k)-rhou_n(i  ,j,k)) &
                           +d4(2)*(rhou_n(i+2,j,k)-rhou_n(i-1,j,k)))*xk4
             Drhov(i,j,k)= (rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk0      &
                         + (d4(1)*(rhov_n(i+1,j,k)-rhov_n(i  ,j,k)) &
                           +d4(2)*(rhov_n(i+2,j,k)-rhov_n(i-1,j,k)))*xk4
             Drhow(i,j,k)= (rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk0      &
                         + (d4(1)*(rhow_n(i+1,j,k)-rhow_n(i  ,j,k)) &
                           +d4(2)*(rhow_n(i+2,j,k)-rhow_n(i-1,j,k)))*xk4
             Drhoe(i,j,k)= (rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk0      &
                         + (d4(1)*(rhoe_n(i+1,j,k)-rhoe_n(i  ,j,k)) &
                           +d4(2)*(rhoe_n(i+2,j,k)-rhoe_n(i-1,j,k)))*xk4
          enddo
       enddo
    endif

    ! Interior points: 5th-order centered dissipation
    ! -----------------------------------------------
    do k=1,nz
       do j=1,ny
          do i=ndx_di,nfx_di
             xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))
             xk6=max(0.0_wp,xnu6+xk0/5.0_wp)
             !uvar(i,j,k,1)=xk6/xnu6

             Drho(i,j,k) = (rho_n(i+1,j,k)-rho_n(i,j,k))*xk0      &
                         + (d6(1)*(rho_n(i+1,j,k)-rho_n(i  ,j,k)) &
                           +d6(2)*(rho_n(i+2,j,k)-rho_n(i-1,j,k)) &
                           +d6(3)*(rho_n(i+3,j,k)-rho_n(i-2,j,k)))*xk6
             Drhou(i,j,k)= (rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk0      &
                         + (d6(1)*(rhou_n(i+1,j,k)-rhou_n(i  ,j,k)) &
                           +d6(2)*(rhou_n(i+2,j,k)-rhou_n(i-1,j,k)) &
                           +d6(3)*(rhou_n(i+3,j,k)-rhou_n(i-2,j,k)))*xk6
             Drhov(i,j,k)= (rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk0      &
                         + (d6(1)*(rhov_n(i+1,j,k)-rhov_n(i  ,j,k)) &
                           +d6(2)*(rhov_n(i+2,j,k)-rhov_n(i-1,j,k)) &
                           +d6(3)*(rhov_n(i+3,j,k)-rhov_n(i-2,j,k)))*xk6
             Drhow(i,j,k)= (rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk0      &
                         + (d6(1)*(rhow_n(i+1,j,k)-rhow_n(i  ,j,k)) &
                           +d6(2)*(rhow_n(i+2,j,k)-rhow_n(i-1,j,k)) &
                           +d6(3)*(rhow_n(i+3,j,k)-rhow_n(i-2,j,k)))*xk6
             Drhoe(i,j,k)= (rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk0      &
                         + (d6(1)*(rhoe_n(i+1,j,k)-rhoe_n(i  ,j,k)) &
                           +d6(2)*(rhoe_n(i+2,j,k)-rhoe_n(i-1,j,k)) &
                           +d6(3)*(rhoe_n(i+3,j,k)-rhoe_n(i-2,j,k)))*xk6
          enddo
       enddo
    enddo

    if (is_boundary(1,2)) then
       ! 3rd-order centered dissipation
       ! ------------------------------
       i=nx-2
       do k=1,nz
          do j=1,ny
             xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))
             xk4=max(0.0_wp,xnu4+xk0)
             !uvar(i,j,k,1)=xk4/xnu4

             Drho(i,j,k) = (rho_n(i+1,j,k)-rho_n(i,j,k))*xk0      &
                         + (d4(1)*(rho_n(i+1,j,k)-rho_n(i  ,j,k)) &
                           +d4(2)*(rho_n(i+2,j,k)-rho_n(i-1,j,k)))*xk4
             Drhou(i,j,k)= (rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk0      &
                         + (d4(1)*(rhou_n(i+1,j,k)-rhou_n(i  ,j,k)) &
                           +d4(2)*(rhou_n(i+2,j,k)-rhou_n(i-1,j,k)))*xk4
             Drhov(i,j,k)= (rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk0      &
                         + (d4(1)*(rhov_n(i+1,j,k)-rhov_n(i  ,j,k)) &
                           +d4(2)*(rhov_n(i+2,j,k)-rhov_n(i-1,j,k)))*xk4
             Drhow(i,j,k)= (rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk0      &
                         + (d4(1)*(rhow_n(i+1,j,k)-rhow_n(i  ,j,k)) &
                           +d4(2)*(rhow_n(i+2,j,k)-rhow_n(i-1,j,k)))*xk4
             Drhoe(i,j,k)= (rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk0      &
                         + (d4(1)*(rhoe_n(i+1,j,k)-rhoe_n(i  ,j,k)) &
                           +d4(2)*(rhoe_n(i+2,j,k)-rhoe_n(i-1,j,k)))*xk4

          enddo
       enddo

       ! 1st-order centered dissipation
       ! ------------------------------
       if (is_bc_wall(1,2)) then
          i=nx-1
          do k=1,nz
             do j=1,ny
                xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)
                !uvar(i,j,k,1)=xk2

                Drho(i,j,k) = (rho_n(i+1,j,k)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i+1,j,k)-rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
       else
          i=nx-1
          do k=1,nz
             do j=1,ny
                xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)
                !uvar(i,j,k,1)=xk2

                Drho(i,j,k) = (rho_n(i+1,j,k)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i+1,j,k)-rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
       endif

       ! no dissipation at i=nx -> initialize increment array
       do j=1,ny
           Krho(nx,j,:)=0.0_wp
          Krhou(nx,j,:)=0.0_wp
          Krhov(nx,j,:)=0.0_wp
          Krhow(nx,j,:)=0.0_wp
          Krhoe(nx,j,:)=0.0_wp
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
                   !mul=min(1.0_wp,cfl_i(i,j,k)/chi)

                   ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dx
                   Krho(i,j,k) = mul*(Drho(i-1,j,k) -Drho(i,j,k))
                   Krhou(i,j,k)= mul*(Drhou(i-1,j,k)-Drhou(i,j,k))
                   Krhov(i,j,k)= mul*(Drhov(i-1,j,k)-Drhov(i,j,k))
                   Krhow(i,j,k)= mul*(Drhow(i-1,j,k)-Drhow(i,j,k))
                   Krhoe(i,j,k)= mul*(Drhoe(i-1,j,k)-Drhoe(i,j,k))
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
                   !mul=min(1.0_wp,cfl_i(i,j,k)/chi)

                   ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dx
                   Krho(i,j,k) = mul*(Drho(i-1,j,k) -Drho(i,j,k))
                   Krhou(i,j,k)= mul*(Drhou(i-1,j,k)-Drhou(i,j,k))
                   Krhov(i,j,k)= mul*(Drhov(i-1,j,k)-Drhov(i,j,k))
                   Krhow(i,j,k)= mul*(Drhow(i-1,j,k)-Drhow(i,j,k))
                   Krhoe(i,j,k)= mul*(Drhoe(i-1,j,k)-Drhoe(i,j,k))
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
                   ! If call to local_cfl
                   !mul=min(1.0_wp,cfl_i(i,j,k)/chi)

                   ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dx
                   Krho(i,j,k) = mul*(Drho(i-1,j,k) -Drho(i,j,k))
                   Krhou(i,j,k)= mul*(Drhou(i-1,j,k)-Drhou(i,j,k))
                   Krhov(i,j,k)= mul*(Drhov(i-1,j,k)-Drhov(i,j,k))
                   Krhow(i,j,k)= mul*(Drhow(i-1,j,k)-Drhow(i,j,k))
                   Krhoe(i,j,k)= mul*(Drhoe(i-1,j,k)-Drhoe(i,j,k))
                enddo
             enddo
          enddo
       endif
    else
       do k=1,nz
          do j=1,ny
             do i=ndx_d,nfx_d
                ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dx
                Krho(i,j,k) = Drho(i-1,j,k) -Drho(i,j,k)
                Krhou(i,j,k)= Drhou(i-1,j,k)-Drhou(i,j,k)
                Krhov(i,j,k)= Drhov(i-1,j,k)-Drhov(i,j,k)
                Krhow(i,j,k)= Drhow(i-1,j,k)-Drhow(i,j,k)
                Krhoe(i,j,k)= Drhoe(i-1,j,k)-Drhoe(i,j,k)
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

    ! Filtering in y-direction [compute Drho at j+1/2]
    ! ========================

    if (is_boundary(2,1)) then
       ! 1st-order centered dissipation
       ! ------------------------------
       if (is_bc_wall(2,1)) then
          j=1
          do k=1,nz
             do i=1,nx
                xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)
                !uvar(i,j,k,2)=xk2

                Drho(i,j,k) = (rho_n(i,j+1,k)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i,j+1,k)-rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
       else
          j=1
          do k=1,nz
             do i=1,nx
                xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)
                !uvar(i,j,k,2)=xk2

                Drho(i,j,k) = (rho_n(i,j+1,k)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i,j+1,k)-rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
       endif

       ! 3rd-order centered dissipation
       ! ------------------------------
       j=2
       do k=1,nz
          do i=1,nx
             xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))
             xk4=max(0.0_wp,xnu4+xk0)
             !uvar(i,j,k,2)=xk4/xnu4

             Drho(i,j,k) = (rho_n(i,j+1,k)-rho_n(i,j,k))*xk0      &
                         + (d4(1)*(rho_n(i,j+1,k)-rho_n(i,j  ,k)) &
                           +d4(2)*(rho_n(i,j+2,k)-rho_n(i,j-1,k)))*xk4
             Drhou(i,j,k)= (rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk0      &
                         + (d4(1)*(rhou_n(i,j+1,k)-rhou_n(i,j  ,k)) &
                           +d4(2)*(rhou_n(i,j+2,k)-rhou_n(i,j-1,k)))*xk4
             Drhov(i,j,k)= (rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk0      &
                         + (d4(1)*(rhov_n(i,j+1,k)-rhov_n(i,j  ,k)) &
                           +d4(2)*(rhov_n(i,j+2,k)-rhov_n(i,j-1,k)))*xk4
             Drhow(i,j,k)= (rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk0      &
                         + (d4(1)*(rhow_n(i,j+1,k)-rhow_n(i,j  ,k)) &
                           +d4(2)*(rhow_n(i,j+2,k)-rhow_n(i,j-1,k)))*xk4
             Drhoe(i,j,k)= (rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk0      &
                         + (d4(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j  ,k)) &
                           +d4(2)*(rhoe_n(i,j+2,k)-rhoe_n(i,j-1,k)))*xk4
          enddo
       enddo
    endif

    ! Interior points: 5th-order centered dissipation
    ! -----------------------------------------------
    do k=1,nz
       do j=ndy_di,nfy_di
          do i=1,nx
             xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))
             xk6=max(0.0_wp,xnu6+xk0/5.0_wp)
             !uvar(i,j,k,2)=xk6/xnu6

             Drho(i,j,k) = (rho_n(i,j+1,k)-rho_n(i,j,k))*xk0      &
                         + (d6(1)*(rho_n(i,j+1,k)-rho_n(i,j  ,k)) &
                           +d6(2)*(rho_n(i,j+2,k)-rho_n(i,j-1,k)) &
                           +d6(3)*(rho_n(i,j+3,k)-rho_n(i,j-2,k)))*xk6
             Drhou(i,j,k)= (rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk0      &
                         + (d6(1)*(rhou_n(i,j+1,k)-rhou_n(i,j  ,k)) &
                           +d6(2)*(rhou_n(i,j+2,k)-rhou_n(i,j-1,k)) &
                           +d6(3)*(rhou_n(i,j+3,k)-rhou_n(i,j-2,k)))*xk6
             Drhov(i,j,k)= (rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk0      &
                         + (d6(1)*(rhov_n(i,j+1,k)-rhov_n(i,j  ,k)) &
                           +d6(2)*(rhov_n(i,j+2,k)-rhov_n(i,j-1,k)) &
                           +d6(3)*(rhov_n(i,j+3,k)-rhov_n(i,j-2,k)))*xk6
             Drhow(i,j,k)= (rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk0      &
                         + (d6(1)*(rhow_n(i,j+1,k)-rhow_n(i,j  ,k)) &
                           +d6(2)*(rhow_n(i,j+2,k)-rhow_n(i,j-1,k)) &
                           +d6(3)*(rhow_n(i,j+3,k)-rhow_n(i,j-2,k)))*xk6
             Drhoe(i,j,k)= (rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk0      &
                         + (d6(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j  ,k)) &
                           +d6(2)*(rhoe_n(i,j+2,k)-rhoe_n(i,j-1,k)) &
                           +d6(3)*(rhoe_n(i,j+3,k)-rhoe_n(i,j-2,k)))*xk6
          enddo
       enddo
    enddo

    if (is_boundary(2,2)) then
       ! 3rd-order centered dissipation
       ! ------------------------------
       j=ny-2
       do k=1,nz
          do i=1,nx
             xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))
             xk4=max(0.0_wp,xnu4+xk0)
             !uvar(i,j,k,2)=xk4/xnu4

             Drho(i,j,k) = (rho_n(i,j+1,k)-rho_n(i,j,k))*xk0      &
                         + (d4(1)*(rho_n(i,j+1,k)-rho_n(i,j  ,k)) &
                           +d4(2)*(rho_n(i,j+2,k)-rho_n(i,j-1,k)))*xk4
             Drhou(i,j,k)= (rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk0      &
                         + (d4(1)*(rhou_n(i,j+1,k)-rhou_n(i,j  ,k)) &
                           +d4(2)*(rhou_n(i,j+2,k)-rhou_n(i,j-1,k)))*xk4
             Drhov(i,j,k)= (rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk0      &
                         + (d4(1)*(rhov_n(i,j+1,k)-rhov_n(i,j  ,k)) &
                           +d4(2)*(rhov_n(i,j+2,k)-rhov_n(i,j-1,k)))*xk4
             Drhow(i,j,k)= (rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk0      &
                         + (d4(1)*(rhow_n(i,j+1,k)-rhow_n(i,j  ,k)) &
                           +d4(2)*(rhow_n(i,j+2,k)-rhow_n(i,j-1,k)))*xk4
             Drhoe(i,j,k)= (rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk0      &
                         + (d4(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j  ,k)) &
                           +d4(2)*(rhoe_n(i,j+2,k)-rhoe_n(i,j-1,k)))*xk4
          enddo
       enddo

       ! 1st-order centered dissipation
       ! ------------------------------
       if (is_bc_wall(2,2)) then
          j=ny-1
          do k=1,nz
             do i=1,nx
                xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)
                !uvar(i,j,k,2)=xk2

                Drho(i,j,k) = (rho_n(i,j+1,k)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i,j+1,k)-rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
       else
          j=ny-1
          do k=1,nz
             do i=1,nx
                xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)
                !uvar(i,j,k,2)=xk2

                Drho(i,j,k) = (rho_n(i,j+1,k)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i,j+1,k)-rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
       endif
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
                   Krho(i,j,k) = Krho(i,j,k) +mul*(Drho(i,j-1,k) -Drho(i,j,k))
                   Krhou(i,j,k)= Krhou(i,j,k)+mul*(Drhou(i,j-1,k)-Drhou(i,j,k))
                   Krhov(i,j,k)= Krhov(i,j,k)+mul*(Drhov(i,j-1,k)-Drhov(i,j,k))
                   Krhow(i,j,k)= Krhow(i,j,k)+mul*(Drhow(i,j-1,k)-Drhow(i,j,k))
                   Krhoe(i,j,k)= Krhoe(i,j,k)+mul*(Drhoe(i,j-1,k)-Drhoe(i,j,k))
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
                   Krho(i,j,k) = Krho(i,j,k) +mul*(Drho(i,j-1,k) -Drho(i,j,k))
                   Krhou(i,j,k)= Krhou(i,j,k)+mul*(Drhou(i,j-1,k)-Drhou(i,j,k))
                   Krhov(i,j,k)= Krhov(i,j,k)+mul*(Drhov(i,j-1,k)-Drhov(i,j,k))
                   Krhow(i,j,k)= Krhow(i,j,k)+mul*(Drhow(i,j-1,k)-Drhow(i,j,k))
                   Krhoe(i,j,k)= Krhoe(i,j,k)+mul*(Drhoe(i,j-1,k)-Drhoe(i,j,k))
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
                   Krho(i,j,k) = Krho(i,j,k) +mul*(Drho(i,j-1,k) -Drho(i,j,k))
                   Krhou(i,j,k)= Krhou(i,j,k)+mul*(Drhou(i,j-1,k)-Drhou(i,j,k))
                   Krhov(i,j,k)= Krhov(i,j,k)+mul*(Drhov(i,j-1,k)-Drhov(i,j,k))
                   Krhow(i,j,k)= Krhow(i,j,k)+mul*(Drhow(i,j-1,k)-Drhow(i,j,k))
                   Krhoe(i,j,k)= Krhoe(i,j,k)+mul*(Drhoe(i,j-1,k)-Drhoe(i,j,k))
                enddo
             enddo
          enddo
       endif
    else
       do k=1,nz
          do j=ndy_d,nfy_d
             do i=1,nx
                ! compute dissipation term: delta_j/delta_y=(D_{j+1/2}-D_{j-1/2})/dy
                Krho(i,j,k) = Krho(i,j,k) +Drho(i,j-1,k) -Drho(i,j,k)
                Krhou(i,j,k)= Krhou(i,j,k)+Drhou(i,j-1,k)-Drhou(i,j,k)
                Krhov(i,j,k)= Krhov(i,j,k)+Drhov(i,j-1,k)-Drhov(i,j,k)
                Krhow(i,j,k)= Krhow(i,j,k)+Drhow(i,j-1,k)-Drhow(i,j,k)
                Krhoe(i,j,k)= Krhoe(i,j,k)+Drhoe(i,j-1,k)-Drhoe(i,j,k)
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

       ! Filtering in k-direction [compute Drho at k+1/2]
       ! ========================
       if (is_boundary(3,1)) then
          ! 1st-order centered dissipation
          ! ------------------------------
          k=1
          do j=1,ny
             do i=1,nx
                xk0=xnu0*max(psens(i,j,k+1),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)

                Drho(i,j,k) = (rho_n(i,j,k+1)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i,j,k+1)-rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i,j,k+1)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i,j,k+1)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i,j,k+1)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i,j,k+1)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k))*xk2
             enddo
          enddo

          ! 3rd-order centered dissipation
          ! ------------------------------
          k=2
          do j=1,ny
             do i=1,nx
                xk0=xnu0*max(psens(i,j,k+1),psens(i,j,k))
                xk4=max(0.0_wp,xnu4+xk0)

                Drho(i,j,k) = (rho_n(i,j,k+1)-rho_n(i,j,k))*xk0      &
                            + (d4(1)*(rho_n(i,j,k+1)-rho_n(i,j,k  )) &
                              +d4(2)*(rho_n(i,j,k+2)-rho_n(i,j,k-1)))*xk4
                Drhou(i,j,k)= (rhou_n(i,j,k+1)-rhou_n(i,j,k))*xk0      &
                            + (d4(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k  )) &
                              +d4(2)*(rhou_n(i,j,k+2)-rhou_n(i,j,k-1)))*xk4
                Drhov(i,j,k)= (rhov_n(i,j,k+1)-rhov_n(i,j,k))*xk0      &
                            + (d4(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k  )) &
                              +d4(2)*(rhov_n(i,j,k+2)-rhov_n(i,j,k-1)))*xk4
                Drhow(i,j,k)= (rhow_n(i,j,k+1)-rhow_n(i,j,k))*xk0      &
                            + (d4(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k  )) &
                              +d4(2)*(rhow_n(i,j,k+2)-rhow_n(i,j,k-1)))*xk4
                Drhoe(i,j,k)= (rhoe_n(i,j,k+1)-rhoe_n(i,j,k))*xk0      &
                            + (d4(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k  )) &
                              +d4(2)*(rhoe_n(i,j,k+2)-rhoe_n(i,j,k-1)))*xk4
             enddo
          enddo
       endif

       ! Interior points: 9th-order centered dissipation
       ! -----------------------------------------------
       do k=ndz_di,nfz_di
          do j=1,ny
             do i=1,nx
                xk0=xnu0*max(psens(i,j,k+1),psens(i,j,k))
                xk6=max(0.0_wp,xnu6+xk0/5.0_wp)

                Drho(i,j,k) = (rho_n(i,j,k+1)-rho_n(i,j,k))*xk0      &
                            + (d6(1)*(rho_n(i,j,k+1)-rho_n(i,j,k  )) &
                              +d6(2)*(rho_n(i,j,k+2)-rho_n(i,j,k-1)) &
                              +d6(3)*(rho_n(i,j,k+3)-rho_n(i,j,k-2)))*xk6
                Drhou(i,j,k)= (rhou_n(i,j,k+1)-rhou_n(i,j,k))*xk0      &
                            + (d6(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k  )) &
                              +d6(2)*(rhou_n(i,j,k+2)-rhou_n(i,j,k-1)) &
                              +d6(3)*(rhou_n(i,j,k+3)-rhou_n(i,j,k-2)))*xk6
                Drhov(i,j,k)= (rhov_n(i,j,k+1)-rhov_n(i,j,k))*xk0      &
                            + (d6(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k  )) &
                              +d6(2)*(rhov_n(i,j,k+2)-rhov_n(i,j,k-1)) &
                              +d6(3)*(rhov_n(i,j,k+3)-rhov_n(i,j,k-2)))*xk6
                Drhow(i,j,k)= (rhow_n(i,j,k+1)-rhow_n(i,j,k))*xk0      &
                            + (d6(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k  )) &
                              +d6(2)*(rhow_n(i,j,k+2)-rhow_n(i,j,k-1)) &
                              +d6(3)*(rhow_n(i,j,k+3)-rhow_n(i,j,k-2)))*xk6
                Drhoe(i,j,k)= (rhoe_n(i,j,k+1)-rhoe_n(i,j,k))*xk0      &
                            + (d6(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k  )) &
                              +d6(2)*(rhoe_n(i,j,k+2)-rhoe_n(i,j,k-1)) &
                              +d6(3)*(rhoe_n(i,j,k+3)-rhoe_n(i,j,k-2)))*xk6
             enddo
          enddo
       enddo

       if (is_boundary(3,2)) then
          ! 3rd-order centered dissipation
          ! ------------------------------
          k=nz-2
          do j=1,ny
             do i=1,nx
                xk0=xnu0*max(psens(i,j,k+1),psens(i,j,k))
                xk4=max(0.0_wp,xnu4+xk0)

                Drho(i,j,k) = (rho_n(i,j,k+1)-rho_n(i,j,k))*xk0      &
                            + (d4(1)*(rho_n(i,j,k+1)-rho_n(i,j,k  )) &
                              +d4(2)*(rho_n(i,j,k+2)-rho_n(i,j,k-1)))*xk4
                Drhou(i,j,k)= (rhou_n(i,j,k+1)-rhou_n(i,j,k))*xk0 &
                            + (d4(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k  )) &
                              +d4(2)*(rhou_n(i,j,k+2)-rhou_n(i,j,k-1)))*xk4
                Drhov(i,j,k)= (rhov_n(i,j,k+1)-rhov_n(i,j,k))*xk0      &
                            + (d4(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k  )) &
                              +d4(2)*(rhov_n(i,j,k+2)-rhov_n(i,j,k-1)))*xk4
                Drhow(i,j,k)= (rhow_n(i,j,k+1)-rhow_n(i,j,k))*xk0      &
                            + (d4(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k  )) &
                              +d4(2)*(rhow_n(i,j,k+2)-rhow_n(i,j,k-1)))*xk4
                Drhoe(i,j,k)= (rhoe_n(i,j,k+1)-rhoe_n(i,j,k))*xk0      &
                            + (d4(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k  )) &
                              +d4(2)*(rhoe_n(i,j,k+2)-rhoe_n(i,j,k-1)))*xk4
             enddo
          enddo

          ! 1st-order centered dissipation
          ! ------------------------------
          k=nz-1
          do j=1,ny
             do i=1,nx
                xk0=xnu0*max(psens(i,j,k+1),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)

                Drho(i,j,k) = (rho_n(i,j,k+1)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i,j,k+1)-rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i,j,k+1)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i,j,k+1)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i,j,k+1)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i,j,k+1)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
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
                      Krho(i,j,k) = Krho(i,j,k) +mul*(Drho(i,j,k-1) -Drho(i,j,k))
                      Krhou(i,j,k)= Krhou(i,j,k)+mul*(Drhou(i,j,k-1)-Drhou(i,j,k))
                      Krhov(i,j,k)= Krhov(i,j,k)+mul*(Drhov(i,j,k-1)-Drhov(i,j,k))
                      Krhow(i,j,k)= Krhow(i,j,k)+mul*(Drhow(i,j,k-1)-Drhow(i,j,k))
                      Krhoe(i,j,k)= Krhoe(i,j,k)+mul*(Drhoe(i,j,k-1)-Drhoe(i,j,k))
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
                      Krho(i,j,k) = Krho(i,j,k) +mul*(Drho(i,j,k-1) -Drho(i,j,k))
                      Krhou(i,j,k)= Krhou(i,j,k)+mul*(Drhou(i,j,k-1)-Drhou(i,j,k))
                      Krhov(i,j,k)= Krhov(i,j,k)+mul*(Drhov(i,j,k-1)-Drhov(i,j,k))
                      Krhow(i,j,k)= Krhow(i,j,k)+mul*(Drhow(i,j,k-1)-Drhow(i,j,k))
                      Krhoe(i,j,k)= Krhoe(i,j,k)+mul*(Drhoe(i,j,k-1)-Drhoe(i,j,k))
                   enddo
                enddo
             enddo
          endif
       else
          do k=ndz_d,nfz_d
             do j=1,ny
                do i=1,nx
                   ! compute dissipation term: delta_k/delta_z=(D_{k+1/2}-D_{k-1/2})/dz
                   Krho(i,j,k) = Krho(i,j,k) +Drho(i,j,k-1) -Drho(i,j,k)
                   Krhou(i,j,k)= Krhou(i,j,k)+Drhou(i,j,k-1)-Drhou(i,j,k)
                   Krhov(i,j,k)= Krhov(i,j,k)+Drhov(i,j,k-1)-Drhov(i,j,k)
                   Krhow(i,j,k)= Krhow(i,j,k)+Drhow(i,j,k-1)-Drhow(i,j,k)
                   Krhoe(i,j,k)= Krhoe(i,j,k)+Drhoe(i,j,k-1)-Drhoe(i,j,k)
                enddo
             enddo
          enddo
       endif

    endif

  end subroutine filter_o6_shock

  !==============================================================================
  subroutine filter_o4_shock
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
    real(wp), dimension(0:nx,0:ny,0:nz) :: Drho,Drhou,Drhov,Drhow,Drhoe
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

    ! Filtering in i-direction [compute Drho at i+1/2]
    ! ========================

    if (is_boundary(1,1)) then
       ! no dissipation at i=1 -> initialize increment array
       do j=1,ny
           Krho(1,j,:)=0.0_wp
          Krhou(1,j,:)=0.0_wp
          Krhov(1,j,:)=0.0_wp
          Krhow(1,j,:)=0.0_wp
          Krhoe(1,j,:)=0.0_wp
       enddo

       ! 1st-order centered dissipation
       ! ------------------------------
       if (is_bc_wall(1,1)) then
          i=1
          do k=1,nz
             do j=1,ny
                xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)
                !uvar(i,j,k,1)=xk2

                Drho(i,j,k) = (rho_n(i+1,j,k)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i+1,j,k)- rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
       else
          i=1
          do k=1,nz
             do j=1,ny
                xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)
                !uvar(i,j,k,1)=xk2

                Drho(i,j,k) = (rho_n(i+1,j,k)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i+1,j,k)- rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
       endif
    endif

    ! 3rd-order centered dissipation
    ! ------------------------------
    do k=1,nz
       do j=1,ny
          do i=ndx_di,nfx_di
             xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))
             xk4=max(0.0_wp,xnu4+xk0)
             !uvar(i,j,k,1)=xk4/xnu4

             Drho(i,j,k) = (rho_n(i+1,j,k)-rho_n(i,j,k))*xk0      &
                         + (d4(1)*(rho_n(i+1,j,k)-rho_n(i  ,j,k)) &
                           +d4(2)*(rho_n(i+2,j,k)-rho_n(i-1,j,k)))*xk4
             Drhou(i,j,k)= (rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk0      &
                         + (d4(1)*(rhou_n(i+1,j,k)-rhou_n(i  ,j,k)) &
                           +d4(2)*(rhou_n(i+2,j,k)-rhou_n(i-1,j,k)))*xk4
             Drhov(i,j,k)= (rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk0      &
                         + (d4(1)*(rhov_n(i+1,j,k)-rhov_n(i  ,j,k)) &
                           +d4(2)*(rhov_n(i+2,j,k)-rhov_n(i-1,j,k)))*xk4
             Drhow(i,j,k)= (rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk0      &
                         + (d4(1)*(rhow_n(i+1,j,k)-rhow_n(i  ,j,k)) &
                           +d4(2)*(rhow_n(i+2,j,k)-rhow_n(i-1,j,k)))*xk4
             Drhoe(i,j,k)= (rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk0      &
                         + (d4(1)*(rhoe_n(i+1,j,k)-rhoe_n(i  ,j,k)) &
                           +d4(2)*(rhoe_n(i+2,j,k)-rhoe_n(i-1,j,k)))*xk4
          enddo
       enddo
    enddo

    if (is_boundary(1,2)) then
       ! 1st-order centered dissipation
       ! ------------------------------
       if (is_bc_wall(1,2)) then
          i=nx-1
          do k=1,nz
             do j=1,ny
                xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)
                !uvar(i,j,k,1)=xk2

                Drho(i,j,k) = (rho_n(i+1,j,k)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i+1,j,k)-rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
       else
          i=nx-1
          do k=1,nz
             do j=1,ny
                xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)
                !uvar(i,j,k,1)=xk2

                Drho(i,j,k) = (rho_n(i+1,j,k)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i+1,j,k)-rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i+1,j,k)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i+1,j,k)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i+1,j,k)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i+1,j,k)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
       endif

       ! no dissipation at i=nx -> initialize increment array
       do j=1,ny
           Krho(nx,j,:)=0.0_wp
          Krhou(nx,j,:)=0.0_wp
          Krhov(nx,j,:)=0.0_wp
          Krhow(nx,j,:)=0.0_wp
          Krhoe(nx,j,:)=0.0_wp
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
                   Krho(i,j,k) = mul*(Drho(i-1,j,k) -Drho(i,j,k))
                   Krhou(i,j,k)= mul*(Drhou(i-1,j,k)-Drhou(i,j,k))
                   Krhov(i,j,k)= mul*(Drhov(i-1,j,k)-Drhov(i,j,k))
                   Krhow(i,j,k)= mul*(Drhow(i-1,j,k)-Drhow(i,j,k))
                   Krhoe(i,j,k)= mul*(Drhoe(i-1,j,k)-Drhoe(i,j,k))
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
                   Krho(i,j,k) = mul*(Drho(i-1,j,k) -Drho(i,j,k))
                   Krhou(i,j,k)= mul*(Drhou(i-1,j,k)-Drhou(i,j,k))
                   Krhov(i,j,k)= mul*(Drhov(i-1,j,k)-Drhov(i,j,k))
                   Krhow(i,j,k)= mul*(Drhow(i-1,j,k)-Drhow(i,j,k))
                   Krhoe(i,j,k)= mul*(Drhoe(i-1,j,k)-Drhoe(i,j,k))
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
                   Krho(i,j,k) = mul*(Drho(i-1,j,k) -Drho(i,j,k))
                   Krhou(i,j,k)= mul*(Drhou(i-1,j,k)-Drhou(i,j,k))
                   Krhov(i,j,k)= mul*(Drhov(i-1,j,k)-Drhov(i,j,k))
                   Krhow(i,j,k)= mul*(Drhow(i-1,j,k)-Drhow(i,j,k))
                   Krhoe(i,j,k)= mul*(Drhoe(i-1,j,k)-Drhoe(i,j,k))
                enddo
             enddo
          enddo
       endif
    else
       do k=1,nz
          do j=1,ny
             do i=ndx_d,nfx_d
                ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dx
                Krho(i,j,k) = Drho(i-1,j,k) -Drho(i,j,k)
                Krhou(i,j,k)= Drhou(i-1,j,k)-Drhou(i,j,k)
                Krhov(i,j,k)= Drhov(i-1,j,k)-Drhov(i,j,k)
                Krhow(i,j,k)= Drhow(i-1,j,k)-Drhow(i,j,k)
                Krhoe(i,j,k)= Drhoe(i-1,j,k)-Drhoe(i,j,k)
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

    ! Filtering in y-direction [compute Drho at j+1/2]
    ! ========================

    if (is_boundary(2,1)) then
       ! 1st-order centered dissipation
       ! ------------------------------
       if (is_bc_wall(2,1)) then
          j=1
          do k=1,nz
             do i=1,nx
                xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)
                !uvar(i,j,k,2)=xk2

                Drho(i,j,k) = (rho_n(i,j+1,k)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i,j+1,k)-rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
       else
          j=1
          do k=1,nz
             do i=1,nx
                xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)
                !uvar(i,j,k,2)=xk2

                Drho(i,j,k) = (rho_n(i,j+1,k)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i,j+1,k)-rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
       endif
    endif

       ! 3rd-order centered dissipation
       ! ------------------------------
    do k=1,nz
       do j=ndy_di,nfy_di
          do i=1,nx
             xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))
             xk4=max(0.0_wp,xnu4+xk0)
             !uvar(i,j,k,2)=xk4/xnu4

             Drho(i,j,k) = (rho_n(i,j+1,k)-rho_n(i,j,k))*xk0      &
                         + (d4(1)*(rho_n(i,j+1,k)-rho_n(i,j  ,k)) &
                           +d4(2)*(rho_n(i,j+2,k)-rho_n(i,j-1,k)))*xk4
             Drhou(i,j,k)= (rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk0      &
                         + (d4(1)*(rhou_n(i,j+1,k)-rhou_n(i,j  ,k)) &
                           +d4(2)*(rhou_n(i,j+2,k)-rhou_n(i,j-1,k)))*xk4
             Drhov(i,j,k)= (rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk0      &
                         + (d4(1)*(rhov_n(i,j+1,k)-rhov_n(i,j  ,k)) &
                           +d4(2)*(rhov_n(i,j+2,k)-rhov_n(i,j-1,k)))*xk4
             Drhow(i,j,k)= (rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk0      &
                         + (d4(1)*(rhow_n(i,j+1,k)-rhow_n(i,j  ,k)) &
                           +d4(2)*(rhow_n(i,j+2,k)-rhow_n(i,j-1,k)))*xk4
             Drhoe(i,j,k)= (rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk0      &
                         + (d4(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j  ,k)) &
                           +d4(2)*(rhoe_n(i,j+2,k)-rhoe_n(i,j-1,k)))*xk4
          enddo
       enddo
    enddo

    if (is_boundary(2,2)) then
       ! 1st-order centered dissipation
       ! ------------------------------
       if (is_bc_wall(2,2)) then
          j=ny-1
          do k=1,nz
             do i=1,nx
                xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))
                xk2=max(0.0_wp,xnu2+xk0)
                !uvar(i,j,k,2)=xk2

                Drho(i,j,k) = (rho_n(i,j+1,k)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i,j+1,k)-rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
       else
          j=ny-1
          do k=1,nz
             do i=1,nx
                xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)
                !uvar(i,j,k,2)=xk2

                Drho(i,j,k) = (rho_n(i,j+1,k)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i,j+1,k)-rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i,j+1,k)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i,j+1,k)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i,j+1,k)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
       endif
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
                   Krho(i,j,k) = Krho(i,j,k) +mul*(Drho(i,j-1,k) -Drho(i,j,k))
                   Krhou(i,j,k)= Krhou(i,j,k)+mul*(Drhou(i,j-1,k)-Drhou(i,j,k))
                   Krhov(i,j,k)= Krhov(i,j,k)+mul*(Drhov(i,j-1,k)-Drhov(i,j,k))
                   Krhow(i,j,k)= Krhow(i,j,k)+mul*(Drhow(i,j-1,k)-Drhow(i,j,k))
                   Krhoe(i,j,k)= Krhoe(i,j,k)+mul*(Drhoe(i,j-1,k)-Drhoe(i,j,k))
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
                   Krho(i,j,k) = Krho(i,j,k) +mul*(Drho(i,j-1,k) -Drho(i,j,k))
                   Krhou(i,j,k)= Krhou(i,j,k)+mul*(Drhou(i,j-1,k)-Drhou(i,j,k))
                   Krhov(i,j,k)= Krhov(i,j,k)+mul*(Drhov(i,j-1,k)-Drhov(i,j,k))
                   Krhow(i,j,k)= Krhow(i,j,k)+mul*(Drhow(i,j-1,k)-Drhow(i,j,k))
                   Krhoe(i,j,k)= Krhoe(i,j,k)+mul*(Drhoe(i,j-1,k)-Drhoe(i,j,k))
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
                   Krho(i,j,k) = Krho(i,j,k) +mul*(Drho(i,j-1,k) -Drho(i,j,k))
                   Krhou(i,j,k)= Krhou(i,j,k)+mul*(Drhou(i,j-1,k)-Drhou(i,j,k))
                   Krhov(i,j,k)= Krhov(i,j,k)+mul*(Drhov(i,j-1,k)-Drhov(i,j,k))
                   Krhow(i,j,k)= Krhow(i,j,k)+mul*(Drhow(i,j-1,k)-Drhow(i,j,k))
                   Krhoe(i,j,k)= Krhoe(i,j,k)+mul*(Drhoe(i,j-1,k)-Drhoe(i,j,k))
                enddo
             enddo
          enddo
       endif
    else
       do k=1,nz
          do j=ndy_d,nfy_d
             do i=1,nx
                ! compute dissipation term: delta_j/delta_x=(D_{j+1/2}-D_{j-1/2})/dx
                Krho(i,j,k) = Krho(i,j,k) +Drho(i,j-1,k) -Drho(i,j,k)
                Krhou(i,j,k)= Krhou(i,j,k)+Drhou(i,j-1,k)-Drhou(i,j,k)
                Krhov(i,j,k)= Krhov(i,j,k)+Drhov(i,j-1,k)-Drhov(i,j,k)
                Krhow(i,j,k)= Krhow(i,j,k)+Drhow(i,j-1,k)-Drhow(i,j,k)
                Krhoe(i,j,k)= Krhoe(i,j,k)+Drhoe(i,j-1,k)-Drhoe(i,j,k)
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

       ! Filtering in k-direction [compute Drho at k+1/2]
       ! ========================
       if (is_boundary(3,1)) then
          ! 1st-order centered dissipation
          ! ------------------------------
          k=1
          do j=1,ny
             do i=1,nx
                xk0=xnu0*max(psens(i,j,k+1),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)

                Drho(i,j,k) = (rho_n(i,j,k+1)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i,j,k+1)-rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i,j,k+1)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i,j,k+1)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i,j,k+1)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i,j,k+1)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
       endif

       ! Interior points: 3rd-order centered dissipation
       ! -----------------------------------------------
       do k=ndz_di,nfz_di
          do j=1,ny
             do i=1,nx
                xk0=xnu0*max(psens(i,j,k+1),psens(i,j,k))
                xk4=max(0.0_wp,xnu4+xk0)

                Drho(i,j,k) = (rho_n(i,j,k+1)-rho_n(i,j,k))*xk0      &
                            + (d4(1)*(rho_n(i,j,k+1)-rho_n(i,j,k  )) &
                              +d4(2)*(rho_n(i,j,k+2)-rho_n(i,j,k-1)))*xk4
                Drhou(i,j,k)= (rhou_n(i,j,k+1)-rhou_n(i,j,k))*xk0      &
                            + (d4(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k  )) &
                              +d4(2)*(rhou_n(i,j,k+2)-rhou_n(i,j,k-1)))*xk4
                Drhov(i,j,k)= (rhov_n(i,j,k+1)-rhov_n(i,j,k))*xk0      &
                            + (d4(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k  )) &
                              +d4(2)*(rhov_n(i,j,k+2)-rhov_n(i,j,k-1)))*xk4
                Drhow(i,j,k)= (rhow_n(i,j,k+1)-rhow_n(i,j,k))*xk0      &
                            + (d4(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k  )) &
                              +d4(2)*(rhow_n(i,j,k+2)-rhow_n(i,j,k-1)))*xk4
                Drhoe(i,j,k)= (rhoe_n(i,j,k+1)-rhoe_n(i,j,k))*xk0      &
                            + (d4(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k  )) &
                              +d4(2)*(rhoe_n(i,j,k+2)-rhoe_n(i,j,k-1)))*xk4
             enddo
          enddo
       enddo

       if (is_boundary(3,2)) then
          ! 1st-order centered dissipation
          ! ------------------------------
          k=nz-1
          do j=1,ny
             do i=1,nx
                xk0=xnu0*max(psens(i,j,k+1),psens(i,j,k))/10.0_wp
                xk2=max(0.0_wp,xnu2+xk0)

                Drho(i,j,k) = (rho_n(i,j,k+1)-rho_n(i,j,k))*xk0 &
                            + d2(1)*(rho_n(i,j,k+1)-rho_n(i,j,k))*xk2
                Drhou(i,j,k)= (rhou_n(i,j,k+1)-rhou_n(i,j,k))*xk0 &
                            + d2(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k))*xk2
                Drhov(i,j,k)= (rhov_n(i,j,k+1)-rhov_n(i,j,k))*xk0 &
                            + d2(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k))*xk2
                Drhow(i,j,k)= (rhow_n(i,j,k+1)-rhow_n(i,j,k))*xk0 &
                            + d2(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k))*xk2
                Drhoe(i,j,k)= (rhoe_n(i,j,k+1)-rhoe_n(i,j,k))*xk0 &
                            + d2(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k))*xk2
             enddo
          enddo
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
                      Krho(i,j,k) = Krho(i,j,k) +mul*(Drho(i,j,k-1) -Drho(i,j,k))
                      Krhou(i,j,k)= Krhou(i,j,k)+mul*(Drhou(i,j,k-1)-Drhou(i,j,k))
                      Krhov(i,j,k)= Krhov(i,j,k)+mul*(Drhov(i,j,k-1)-Drhov(i,j,k))
                      Krhow(i,j,k)= Krhow(i,j,k)+mul*(Drhow(i,j,k-1)-Drhow(i,j,k))
                      Krhoe(i,j,k)= Krhoe(i,j,k)+mul*(Drhoe(i,j,k-1)-Drhoe(i,j,k))
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
                      Krho(i,j,k) = Krho(i,j,k) +mul*(Drho(i,j,k-1) -Drho(i,j,k))
                      Krhou(i,j,k)= Krhou(i,j,k)+mul*(Drhou(i,j,k-1)-Drhou(i,j,k))
                      Krhov(i,j,k)= Krhov(i,j,k)+mul*(Drhov(i,j,k-1)-Drhov(i,j,k))
                      Krhow(i,j,k)= Krhow(i,j,k)+mul*(Drhow(i,j,k-1)-Drhow(i,j,k))
                      Krhoe(i,j,k)= Krhoe(i,j,k)+mul*(Drhoe(i,j,k-1)-Drhoe(i,j,k))
                   enddo
                enddo
             enddo
          endif
       else
          do k=ndz_d,nfz_d
             do j=1,ny
                do i=1,nx
                   ! compute dissipation term: delta_k/delta_z=(D_{k+1/2}-D_{k-1/2})/dz
                   Krho(i,j,k) = Krho(i,j,k) +Drho(i,j,k-1) -Drho(i,j,k)
                   Krhou(i,j,k)= Krhou(i,j,k)+Drhou(i,j,k-1)-Drhou(i,j,k)
                   Krhov(i,j,k)= Krhov(i,j,k)+Drhov(i,j,k-1)-Drhov(i,j,k)
                   Krhow(i,j,k)= Krhow(i,j,k)+Drhow(i,j,k-1)-Drhow(i,j,k)
                   Krhoe(i,j,k)= Krhoe(i,j,k)+Drhoe(i,j,k-1)-Drhoe(i,j,k)
                enddo
             enddo
          enddo
       endif

    endif

  end subroutine filter_o4_shock

end module mod_filtering_shock
