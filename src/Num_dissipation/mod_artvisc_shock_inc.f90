!==============================================================================
module mod_artvisc_shock_inc
!==============================================================================
  !> Module for artificial viscosity with Jameson's type shock capturing
  !> Added to increments (for IRS smoothing)
!==============================================================================
  use mod_flow
  use mod_artvisc
  implicit none
  ! ---------------------------------------------------------------------------
  ! declaration of coefficients in mod_artvisc.f90
  real(wp), parameter, private :: cutoff = 1.e-16_wp
  ! amplitude dissipation
  real(wp), private :: xnu0 ! low-order smoothing term
  real(wp), private :: xnu10,xnu8,xnu6,xnu4,xnu2
  ! ---------------------------------------------------------------------------

contains

  !============================================================================
  subroutine init_artvisc_shock_inc(dissip_coeff,dissip_shock)
  !============================================================================
    !> Initialization of artificial viscosity
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    real(wp), intent(in) :: dissip_coeff,dissip_shock
    ! -------------------------------------------------------------------------

    ! Artificial viscosity coefficients (DNC-Jameson)
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

    ! low-order dissipation term [from param.ini]
    ! ==========================
    ! /!\ Change sign because added to increment
    xnu0=-dissip_shock

  end subroutine init_artvisc_shock_inc

  !==============================================================================
  subroutine artvisc_o9_shock_inc
  !==============================================================================
    !> Apply artificial viscosity - Conservative version DNC order 9
    !> Modified version XG 06/2021
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: uc1,uc2,uc
    real(wp) :: vc1,vc2,vc,srm12,srp12,sr
    real(wp), dimension(0:nx,0:ny,0:nz) :: Drho,Drhou,Drhov,Drhow,Drhoe
    ! ---------------------------------------------------------------------------
    real(wp) :: div2,om1,om2,om3
    real(wp) :: xk0,xk2,xk4,xk6,xk8,xk10
    real(wp), dimension(ndx_d1:nfx_d1,ndy_d1:nfy_d1,ndz_d1:nfz_d1) :: sens,psens
    ! ---------------------------------------------------------------------------
    real(wp) :: psens_trad,psens_tvd,csens

    ! Compute Ducros sensor
    ! =====================
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
    !sens=1.0_wp
   
    csens = 0.01_wp
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
 
    ! Artificial viscosity in i-direction [compute Drho at i+1/2]
    ! ===================================
    
    if (is_boundary(1,1)) then
       ! 1st-order centered dissipation
       ! ------------------------------
       i=1
       do k=1,nz
          do j=1,ny             
             xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))
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
       
       ! 3rd-order centered dissipation
       ! ------------------------------
       i=2
       do k=1,nz
          do j=1,ny            
             xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))
             xk4=max(0.0_wp,xnu4+xk0)
             !uvar(i,j,k,1)=xk4
             
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
             xk6=max(0.0_wp,xnu6+xk0)
             !uvar(i,j,k,1)=xk6
            
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
             xk8=max(0.0_wp,xnu8+xk0)
             !uvar(i,j,k,1)=xk8
             
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
             xk10=max(0.0_wp,xnu10+xk0)
             !uvar(i,j,k,1)=xk10
             
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
             xk8=max(0.0_wp,xnu8+xk0)
             !uvar(i,j,k,1)=xk8
            
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
             xk6=max(0.0_wp,xnu6+xk0)
             !uvar(i,j,k,1)=xk6
             
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
             !uvar(i,j,k,1)=xk4
             
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
       i=nx-1
       do k=1,nz
          do j=1,ny             
             xk0=xnu0*max(psens(i+1,j,k),psens(i,j,k))
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
    
    ! Add dissipation to numerical fluxes
    ! -----------------------------------
    if (is_curv) then
       do k=1,nz
          do j=1,ny
             do i=ndx_d,nfx_d
                ! contravariant velocity at i+1/2 (vc2) and i-1/2 (vc1)
                uc=(uu(i,j,k)*y_eta(i,j)-vv(i,j,k)*x_eta(i,j))*ijacob(i,j)
                uc2=(uu(i+1,j,k)*y_eta(i+1,j)-vv(i+1,j,k)*x_eta(i+1,j))*ijacob(i+1,j)
                uc1=(uu(i-1,j,k)*y_eta(i-1,j)-vv(i-1,j,k)*x_eta(i-1,j))*ijacob(i-1,j)
                vc=(vv(i,j,k)*x_ksi(i,j)-uu(i,j,k)*y_ksi(i,j))*ijacob(i,j)
                vc2=(vv(i+1,j,k)*x_ksi(i+1,j)-uu(i+1,j,k)*y_ksi(i+1,j))*ijacob(i+1,j)
                vc1=(vv(i-1,j,k)*x_ksi(i-1,j)-uu(i-1,j,k)*y_ksi(i-1,j))*ijacob(i-1,j)
                ! compute spectral radius at i+1/2 (srp12) and i-1/2 (srm12)
                sr=abs(uc)+c_(i,j,k)*g_ksi(i,j)
                srp12=0.5_wp*(abs(uc2)+c_(i+1,j,k)*g_ksi(i+1,j)+sr)
                srm12=0.5_wp*(abs(uc1)+c_(i-1,j,k)*g_ksi(i-1,j)+sr)
!!$                sr=abs(uc)+c_(i,j,k)*g_ksi(i,j) &
!!$                  +abs(vc)+c_(i,j,k)*g_eta(i,j)
!!$                srp12=0.5_wp*abs(uc2)+c_(i+1,j,k)*g_ksi(i+1,j)+0.5_wp*sr &
!!$                   +0.5_wp*abs(vc2)+c_(i+1,j,k)*g_eta(i+1,j)
!!$                srm12=0.5_wp*abs(uc1)+c_(i-1,j,k)*g_ksi(i-1,j)+0.5_wp*sr &
!!$                     +0.5_wp*abs(vc1)+c_(i-1,j,k)*g_eta(i-1,j)
                ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dksi
                Krho(i,j,k) = Krho(i,j,k) +(srp12*Drho(i,j,k) -srm12*Drho(i-1,j,k))
                Krhou(i,j,k)= Krhou(i,j,k)+(srp12*Drhou(i,j,k)-srm12*Drhou(i-1,j,k))
                Krhov(i,j,k)= Krhov(i,j,k)+(srp12*Drhov(i,j,k)-srm12*Drhov(i-1,j,k))
                Krhow(i,j,k)= Krhow(i,j,k)+(srp12*Drhow(i,j,k)-srm12*Drhow(i-1,j,k))
                Krhoe(i,j,k)= Krhoe(i,j,k)+(srp12*Drhoe(i,j,k)-srm12*Drhoe(i-1,j,k))
                ! store dissipation coefficient
                !uvar(i,j,k,1)=sr
             enddo
          enddo
       enddo
    else
       do k=1,nz
          do j=1,ny
             do i=ndx_d,nfx_d
                ! compute spectral radius at i+1/2 (srp12) and i-1/2 (srm12)
                sr=abs(uu(i,j,k))+c_(i,j,k)
                srp12=0.5_wp*(abs(uu(i+1,j,k))+c_(i+1,j,k)+sr)
                srm12=0.5_wp*(abs(uu(i-1,j,k))+c_(i-1,j,k)+sr)
                ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dx
                Krho(i,j,k) = Krho(i,j,k) +(srp12*Drho(i,j,k) -srm12*Drho(i-1,j,k) )*idx(i)
                Krhou(i,j,k)= Krhou(i,j,k)+(srp12*Drhou(i,j,k)-srm12*Drhou(i-1,j,k))*idx(i)
                Krhov(i,j,k)= Krhov(i,j,k)+(srp12*Drhov(i,j,k)-srm12*Drhov(i-1,j,k))*idx(i)
                Krhow(i,j,k)= Krhow(i,j,k)+(srp12*Drhow(i,j,k)-srm12*Drhow(i-1,j,k))*idx(i)
                Krhoe(i,j,k)= Krhoe(i,j,k)+(srp12*Drhoe(i,j,k)-srm12*Drhoe(i-1,j,k))*idx(i)
                ! store dissipation coefficient
                !uvar(i,j,k,1)=sr
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
       
       ! 3rd-order centered dissipation
       ! ------------------------------
       j=2
       do k=1,nz
          do i=1,nx
             xk0=xnu0*max(psens(i,j+1,k),psens(i,j,k))
             xk4=max(0.0_wp,xnu4+xk0)
             !uvar(i,j,k,2)=xk4
             
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
             xk6=max(0.0_wp,xnu6+xk0)
             !uvar(i,j,k,2)=xk6
             
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
             xk8=max(0.0_wp,xnu8+xk0)
             !uvar(i,j,k,2)=xk8
             
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
             xk10=max(0.0_wp,xnu10+xk0)
             !uvar(i,j,k,2)=xk10
             
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
             xk8=max(0.0_wp,xnu8+xk0)
             !uvar(i,j,k,2)=xk8
             
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
             xk6=max(0.0_wp,xnu6+xk0)
             !uvar(i,j,k,2)=xk6
             
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
             !uvar(i,j,k,2)=xk4
             
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
    endif

    ! Add dissipation to numerical fluxes
    ! -----------------------------------
    if (is_curv) then
       do k=1,nz
          do j=ndy_d,nfy_d
             do i=1,nx
                ! contravariant velocity at j+1/2 (vc2) and j-1/2 (vc1)
                uc=(uu(i,j,k)*y_eta(i,j)-vv(i,j,k)*x_eta(i,j))*ijacob(i,j)
                uc2=(uu(i,j+1,k)*y_eta(i,j+1)-vv(i,j+1,k)*x_eta(i,j+1))*ijacob(i,j+1)
                uc1=(uu(i,j-1,k)*y_eta(i,j-1)-vv(i,j-1,k)*x_eta(i,j-1))*ijacob(i,j-1)
                vc=(vv(i,j,k)*x_ksi(i,j)-uu(i,j,k)*y_ksi(i,j))*ijacob(i,j)
                vc2=(vv(i,j+1,k)*x_ksi(i,j+1)-uu(i,j+1,k)*y_ksi(i,j+1))*ijacob(i,j+1)
                vc1=(vv(i,j-1,k)*x_ksi(i,j-1)-uu(i,j-1,k)*y_ksi(i,j-1))*ijacob(i,j-1)
                ! compute spectral radius at j+1/2 (srp12) and j-1/2 (srm12)
                sr=abs(vc)+c_(i,j,k)*g_eta(i,j)
                srp12=0.5_wp*(abs(vc2)+c_(i,j+1,k)*g_eta(i,j+1)+sr)
                srm12=0.5_wp*(abs(vc1)+c_(i,j-1,k)*g_eta(i,j-1)+sr)
!!$                sr=abs(uc)+c_(i,j,k)*g_ksi(i,j) &
!!$                  +abs(vc)+c_(i,j,k)*g_eta(i,j)
!!$                srp12=0.5_wp*abs(uc2)+c_(i,j+1,k)*g_ksi(i,j+1)+0.5_wp*sr &
!!$                   +0.5_wp*abs(vc2)+c_(i,j+1,k)*g_eta(i,j+1)
!!$                srm12=0.5_wp*abs(uc1)+c_(i,j-1,k)*g_ksi(i,j-1)+0.5_wp*sr &
!!$                   +0.5_wp*abs(vc1)+c_(i,j-1,k)*g_eta(i,j-1)
                ! compute dissipation term: delta_j/delta_y=(D_{j+1/2}-D_{j-1/2})/dy
                Krho(i,j,k) = Krho(i,j,k) +(srp12*Drho(i,j,k) -srm12*Drho(i,j-1,k) )
                Krhou(i,j,k)= Krhou(i,j,k)+(srp12*Drhou(i,j,k)-srm12*Drhou(i,j-1,k))
                Krhov(i,j,k)= Krhov(i,j,k)+(srp12*Drhov(i,j,k)-srm12*Drhov(i,j-1,k))
                Krhow(i,j,k)= Krhow(i,j,k)+(srp12*Drhow(i,j,k)-srm12*Drhow(i,j-1,k))
                Krhoe(i,j,k)= Krhoe(i,j,k)+(srp12*Drhoe(i,j,k)-srm12*Drhoe(i,j-1,k))
                ! store dissipation coefficient
                !uvar(i,j,k,2)=sr
             enddo
          enddo
       enddo
    else
       do k=1,nz
          do j=ndy_d,nfy_d
             do i=1,nx
                ! compute spectral radius at j+1/2 (srp12) and j-1/2 (srm12)
                sr=abs(vv(i,j,k))+c_(i,j,k)
                srp12=0.5_wp*(abs(vv(i,j+1,k))+c_(i,j+1,k)+sr)
                srm12=0.5_wp*(abs(vv(i,j-1,k))+c_(i,j-1,k)+sr)
                ! compute dissipation term: delta_j/delta_y=(D_{j+1/2}-D_{j-1/2})/dy
                Krho(i,j,k) = Krho(i,j,k) +(srp12*Drho(i,j,k) -srm12*Drho(i,j-1,k) )*idy(j)
                Krhou(i,j,k)= Krhou(i,j,k)+(srp12*Drhou(i,j,k)-srm12*Drhou(i,j-1,k))*idy(j)
                Krhov(i,j,k)= Krhov(i,j,k)+(srp12*Drhov(i,j,k)-srm12*Drhov(i,j-1,k))*idy(j)
                Krhow(i,j,k)= Krhow(i,j,k)+(srp12*Drhow(i,j,k)-srm12*Drhow(i,j-1,k))*idy(j)
                Krhoe(i,j,k)= Krhoe(i,j,k)+(srp12*Drhoe(i,j,k)-srm12*Drhoe(i,j-1,k))*idy(j)
                ! store dissipation coefficient
                !uvar(i,j,k,2)=sr
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
                xk10=max(0.0_wp,xnu10+xk0)
             
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
                xk0=xnu0*max(psens(i,j,k+1),psens(i,j,k))
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

       ! Add dissipation to numerical fluxes
       ! -----------------------------------
       do k=ndz_d,nfz_d
          do j=1,ny
             do i=1,nx
                ! compute spectral radius at k+1/2 (srp12) and k-1/2 (srm12)
                sr=abs(ww(i,j,k))+c_(i,j,k)
                srp12=0.5_wp*(abs(ww(i,j,k+1))+c_(i,j,k+1)+sr)
                srm12=0.5_wp*(abs(ww(i,j,k-1))+c_(i,j,k-1)+sr)
                ! compute dissipation term: delta_k/delta_z=(D_{k+1/2}-D_{k-1/2})/dz
                Krho(i,j,k) = Krho(i,j,k) +(srp12*Drho(i,j,k) -srm12*Drho(i,j,k-1) )*idz(k)
                Krhou(i,j,k)= Krhou(i,j,k)+(srp12*Drhou(i,j,k)-srm12*Drhou(i,j,k-1))*idz(k)
                Krhov(i,j,k)= Krhov(i,j,k)+(srp12*Drhov(i,j,k)-srm12*Drhov(i,j,k-1))*idz(k)
                Krhow(i,j,k)= Krhow(i,j,k)+(srp12*Drhow(i,j,k)-srm12*Drhow(i,j,k-1))*idz(k)
                Krhoe(i,j,k)= Krhoe(i,j,k)+(srp12*Drhoe(i,j,k)-srm12*Drhoe(i,j,k-1))*idz(k)
             enddo
          enddo
       enddo

    endif
    
  end subroutine artvisc_o9_shock_inc

end module mod_artvisc_shock_inc
