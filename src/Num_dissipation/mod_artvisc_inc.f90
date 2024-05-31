!================================================================================
module mod_artvisc_inc
!================================================================================
  !> Module for artificial viscosity
  !> Added to increments (for IRS smoothing)
!================================================================================
  use mod_flow
  use mod_artvisc
  implicit none
  ! -----------------------------------------------------------------------------
  ! declaration of coefficients in mod_artvisc.f90
  ! -----------------------------------------------------------------------------

contains

  !==============================================================================
  subroutine init_artvisc_inc(dissip_coeff)
  !==============================================================================
    !> Initialization of artificial viscosity
  !==============================================================================
    use mod_time
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), intent(in) :: dissip_coeff
    real(wp) :: xnu10,xnu8,xnu6,xnu4,xnu2
    ! ---------------------------------------------------------------------------

    ! Artificial viscosity coefficients (DNC-Jameson)
    ! =================================
    ! order 1
    d2(1)= 1.0_wp

    ! order 3
    d4(1)= 3.0_wp
    d4(2)=-1.0_wp

    ! order 5
    d6(1)=10.0_wp
    d6(2)=-5.0_wp
    d6(3)= 1.0_wp

    ! order 7
    d8(1)= 35.0_wp
    d8(2)=-21.0_wp
    d8(3)=  7.0_wp
    d8(4)= -1.0_wp

    ! order 9
    d10(1)= 126.0_wp
    d10(2)=-84.0_wp
    d10(3)= 36.0_wp
    d10(4)=-9.0_wp
    d10(5)= 1.0_wp
    
    ! dissipation coefficient [from param.ini]
    xnu10=dissip_coeff/1260.0_wp
    xnu8= dissip_coeff/280.0_wp
    xnu6= dissip_coeff/60.0_wp
    xnu4= dissip_coeff/12.0_wp
    xnu2= 0.0_wp*dissip_coeff/2.0_wp

!!$    xnu8= xnu10
!!$    xnu6= xnu10
!!$    xnu4= xnu10
!!$    xnu2= xnu10
   
!!$   xnu8= dissip_coeff/280.0_wp/4.
!!$   xnu6= dissip_coeff/60.0_wp/6.
!!$   xnu4= dissip_coeff/12.0_wp/8.
!!$   xnu2= dissip_coeff/2.0_wp/10.

    ! Multiply coeff by dissipation coefficient and change sign
    ! =========================================================
    ! /!\ Change sign because added to increment
    d10=-d10*xnu10
     d8= -d8*xnu8
     d6= -d6*xnu6
     d4= -d4*xnu4
     d2= -d2*xnu2
   
  end subroutine init_artvisc_inc

  !==============================================================================
  subroutine artvisc_o9_inc
  !==============================================================================
    !> Apply artificial viscosity - Conservative version DNC order 9
    !> Modified version XG 06/2022
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: uc1,uc2,uc
    real(wp) :: vc1,vc2,vc,srm12,srp12,sr
    real(wp), dimension(0:nx,0:ny,0:nz) :: Drho,Drhou,Drhov,Drhow,Drhoe
    ! ---------------------------------------------------------------------------

    ! Artificial viscosity in x-direction [compute Drho at i+1/2]
    ! ===================================
    
    if (is_boundary(1,1)) then
       ! 1st-order centered dissipation
       ! ------------------------------
       i=1
       do k=1,nz
          do j=1,ny
             Drho(i,j,k) = d2(1)*(rho_n(i+1,j,k) -rho_n(i,j,k))            
             Drhou(i,j,k)= d2(1)*(rhou_n(i+1,j,k)-rhou_n(i,j,k))
             Drhov(i,j,k)= d2(1)*(rhov_n(i+1,j,k)-rhov_n(i,j,k))
             Drhow(i,j,k)= d2(1)*(rhow_n(i+1,j,k)-rhow_n(i,j,k))
             Drhoe(i,j,k)= d2(1)*(rhoe_n(i+1,j,k)-rhoe_n(i,j,k))
          enddo
       enddo
       
       ! 3rd-order centered dissipation
       ! ------------------------------
       i=2
       do k=1,nz
          do j=1,ny
             Drho(i,j,k) = d4(1)*(rho_n(i+1,j,k)-rho_n(i  ,j,k)) &
                         + d4(2)*(rho_n(i+2,j,k)-rho_n(i-1,j,k))
             Drhou(i,j,k)= d4(1)*(rhou_n(i+1,j,k)-rhou_n(i  ,j,k)) &
                         + d4(2)*(rhou_n(i+2,j,k)-rhou_n(i-1,j,k))
             Drhov(i,j,k)= d4(1)*(rhov_n(i+1,j,k)-rhov_n(i  ,j,k)) &
                         + d4(2)*(rhov_n(i+2,j,k)-rhov_n(i-1,j,k))
             Drhow(i,j,k)= d4(1)*(rhow_n(i+1,j,k)-rhow_n(i  ,j,k)) &
                         + d4(2)*(rhow_n(i+2,j,k)-rhow_n(i-1,j,k))
             Drhoe(i,j,k)= d4(1)*(rhoe_n(i+1,j,k)-rhoe_n(i  ,j,k)) &
                         + d4(2)*(rhoe_n(i+2,j,k)-rhoe_n(i-1,j,k))
          enddo
       enddo
       
       ! 5th-order centered dissipation
       ! ------------------------------
       i=3
       do k=1,nz
          do j=1,ny
             Drho(i,j,k) = d6(1)*(rho_n(i+1,j,k)-rho_n(i  ,j,k)) &
                         + d6(2)*(rho_n(i+2,j,k)-rho_n(i-1,j,k)) &
                         + d6(3)*(rho_n(i+3,j,k)-rho_n(i-2,j,k))
             Drhou(i,j,k)= d6(1)*(rhou_n(i+1,j,k)-rhou_n(i  ,j,k)) &
                         + d6(2)*(rhou_n(i+2,j,k)-rhou_n(i-1,j,k)) &
                         + d6(3)*(rhou_n(i+3,j,k)-rhou_n(i-2,j,k))
             Drhov(i,j,k)= d6(1)*(rhov_n(i+1,j,k)-rhov_n(i  ,j,k)) &
                         + d6(2)*(rhov_n(i+2,j,k)-rhov_n(i-1,j,k)) &
                         + d6(3)*(rhov_n(i+3,j,k)-rhov_n(i-2,j,k))
             Drhow(i,j,k)= d6(1)*(rhow_n(i+1,j,k)-rhow_n(i  ,j,k)) &
                         + d6(2)*(rhow_n(i+2,j,k)-rhow_n(i-1,j,k)) &
                         + d6(3)*(rhow_n(i+3,j,k)-rhow_n(i-2,j,k))
             Drhoe(i,j,k)= d6(1)*(rhoe_n(i+1,j,k)-rhoe_n(i  ,j,k)) &
                         + d6(2)*(rhoe_n(i+2,j,k)-rhoe_n(i-1,j,k)) &
                         + d6(3)*(rhoe_n(i+3,j,k)-rhoe_n(i-2,j,k))
          enddo
       enddo
       
       ! 7th-order centered dissipation
       ! ------------------------------
       i=4
       do k=1,nz
          do j=1,ny
             Drho(i,j,k) = d8(1)*(rho_n(i+1,j,k)-rho_n(i  ,j,k)) &
                         + d8(2)*(rho_n(i+2,j,k)-rho_n(i-1,j,k)) &
                         + d8(3)*(rho_n(i+3,j,k)-rho_n(i-2,j,k)) &
                         + d8(4)*(rho_n(i+4,j,k)-rho_n(i-3,j,k))
             Drhou(i,j,k)= d8(1)*(rhou_n(i+1,j,k)-rhou_n(i  ,j,k)) &
                         + d8(2)*(rhou_n(i+2,j,k)-rhou_n(i-1,j,k)) &
                         + d8(3)*(rhou_n(i+3,j,k)-rhou_n(i-2,j,k)) &
                         + d8(4)*(rhou_n(i+4,j,k)-rhou_n(i-3,j,k))
             Drhov(i,j,k)= d8(1)*(rhov_n(i+1,j,k)-rhov_n(i  ,j,k)) &
                         + d8(2)*(rhov_n(i+2,j,k)-rhov_n(i-1,j,k)) &
                         + d8(3)*(rhov_n(i+3,j,k)-rhov_n(i-2,j,k)) &
                         + d8(4)*(rhov_n(i+4,j,k)-rhov_n(i-3,j,k))
             Drhow(i,j,k)= d8(1)*(rhow_n(i+1,j,k)-rhow_n(i  ,j,k)) &
                         + d8(2)*(rhow_n(i+2,j,k)-rhow_n(i-1,j,k)) &
                         + d8(3)*(rhow_n(i+3,j,k)-rhow_n(i-2,j,k)) &
                         + d8(4)*(rhow_n(i+4,j,k)-rhow_n(i-3,j,k))
             Drhoe(i,j,k)= d8(1)*(rhoe_n(i+1,j,k)-rhoe_n(i  ,j,k)) &
                         + d8(2)*(rhoe_n(i+2,j,k)-rhoe_n(i-1,j,k)) &
                         + d8(3)*(rhoe_n(i+3,j,k)-rhoe_n(i-2,j,k)) &
                         + d8(4)*(rhoe_n(i+4,j,k)-rhoe_n(i-3,j,k))
          enddo
       enddo
    endif

    ! Interior points: 9th-order centered dissipation
    ! -----------------------------------------------
    do k=1,nz
       do j=1,ny
          do i=ndx_di,nfx_di
             Drho(i,j,k) = d10(1)*(rho_n(i+1,j,k)-rho_n(i  ,j,k)) &
                         + d10(2)*(rho_n(i+2,j,k)-rho_n(i-1,j,k)) &
                         + d10(3)*(rho_n(i+3,j,k)-rho_n(i-2,j,k)) &
                         + d10(4)*(rho_n(i+4,j,k)-rho_n(i-3,j,k)) &
                         + d10(5)*(rho_n(i+5,j,k)-rho_n(i-4,j,k))
             
             Drhou(i,j,k)= d10(1)*(rhou_n(i+1,j,k)-rhou_n(i  ,j,k)) &
                         + d10(2)*(rhou_n(i+2,j,k)-rhou_n(i-1,j,k)) &
                         + d10(3)*(rhou_n(i+3,j,k)-rhou_n(i-2,j,k)) &
                         + d10(4)*(rhou_n(i+4,j,k)-rhou_n(i-3,j,k)) &
                         + d10(5)*(rhou_n(i+5,j,k)-rhou_n(i-4,j,k))

             Drhov(i,j,k)= d10(1)*(rhov_n(i+1,j,k)-rhov_n(i  ,j,k)) &
                         + d10(2)*(rhov_n(i+2,j,k)-rhov_n(i-1,j,k)) &
                         + d10(3)*(rhov_n(i+3,j,k)-rhov_n(i-2,j,k)) &
                         + d10(4)*(rhov_n(i+4,j,k)-rhov_n(i-3,j,k)) &
                         + d10(5)*(rhov_n(i+5,j,k)-rhov_n(i-4,j,k))

             Drhow(i,j,k)= d10(1)*(rhow_n(i+1,j,k)-rhow_n(i  ,j,k)) &
                         + d10(2)*(rhow_n(i+2,j,k)-rhow_n(i-1,j,k)) &
                         + d10(3)*(rhow_n(i+3,j,k)-rhow_n(i-2,j,k)) &
                         + d10(4)*(rhow_n(i+4,j,k)-rhow_n(i-3,j,k)) &
                         + d10(5)*(rhow_n(i+5,j,k)-rhow_n(i-4,j,k))

             Drhoe(i,j,k)= d10(1)*(rhoe_n(i+1,j,k)-rhoe_n(i  ,j,k)) &
                         + d10(2)*(rhoe_n(i+2,j,k)-rhoe_n(i-1,j,k)) &
                         + d10(3)*(rhoe_n(i+3,j,k)-rhoe_n(i-2,j,k)) &
                         + d10(4)*(rhoe_n(i+4,j,k)-rhoe_n(i-3,j,k)) &
                         + d10(5)*(rhoe_n(i+5,j,k)-rhoe_n(i-4,j,k))
          enddo
       enddo
    enddo

    if (is_boundary(1,2)) then
       ! 7th-order centered dissipation
       ! ------------------------------
       i=nx-4
       do k=1,nz
          do j=1,ny
             Drho(i,j,k) = d8(1)*(rho_n(i+1,j,k)-rho_n(i  ,j,k)) &
                         + d8(2)*(rho_n(i+2,j,k)-rho_n(i-1,j,k)) &
                         + d8(3)*(rho_n(i+3,j,k)-rho_n(i-2,j,k)) &
                         + d8(4)*(rho_n(i+4,j,k)-rho_n(i-3,j,k))
             Drhou(i,j,k)= d8(1)*(rhou_n(i+1,j,k)-rhou_n(i  ,j,k)) &
                         + d8(2)*(rhou_n(i+2,j,k)-rhou_n(i-1,j,k)) &
                         + d8(3)*(rhou_n(i+3,j,k)-rhou_n(i-2,j,k)) &
                         + d8(4)*(rhou_n(i+4,j,k)-rhou_n(i-3,j,k))
             Drhov(i,j,k)= d8(1)*(rhov_n(i+1,j,k)-rhov_n(i  ,j,k)) &
                         + d8(2)*(rhov_n(i+2,j,k)-rhov_n(i-1,j,k)) &
                         + d8(3)*(rhov_n(i+3,j,k)-rhov_n(i-2,j,k)) &
                         + d8(4)*(rhov_n(i+4,j,k)-rhov_n(i-3,j,k))
             Drhow(i,j,k)= d8(1)*(rhow_n(i+1,j,k)-rhow_n(i  ,j,k)) &
                         + d8(2)*(rhow_n(i+2,j,k)-rhow_n(i-1,j,k)) &
                         + d8(3)*(rhow_n(i+3,j,k)-rhow_n(i-2,j,k)) &
                         + d8(4)*(rhow_n(i+4,j,k)-rhow_n(i-3,j,k))
             Drhoe(i,j,k)= d8(1)*(rhoe_n(i+1,j,k)-rhoe_n(i  ,j,k)) &
                         + d8(2)*(rhoe_n(i+2,j,k)-rhoe_n(i-1,j,k)) &
                         + d8(3)*(rhoe_n(i+3,j,k)-rhoe_n(i-2,j,k)) &
                         + d8(4)*(rhoe_n(i+4,j,k)-rhoe_n(i-3,j,k))
          enddo
       enddo
       
       ! 5th-order centered dissipation
       ! ------------------------------
       i=nx-3
       do k=1,nz
          do j=1,ny
             Drho(i,j,k) = d6(1)*(rho_n(i+1,j,k)-rho_n(i  ,j,k)) &
                         + d6(2)*(rho_n(i+2,j,k)-rho_n(i-1,j,k)) &
                         + d6(3)*(rho_n(i+3,j,k)-rho_n(i-2,j,k))    
             Drhou(i,j,k)= d6(1)*(rhou_n(i+1,j,k)-rhou_n(i  ,j,k)) &
                         + d6(2)*(rhou_n(i+2,j,k)-rhou_n(i-1,j,k)) &
                         + d6(3)*(rhou_n(i+3,j,k)-rhou_n(i-2,j,k))
             Drhov(i,j,k)= d6(1)*(rhov_n(i+1,j,k)-rhov_n(i  ,j,k)) &
                         + d6(2)*(rhov_n(i+2,j,k)-rhov_n(i-1,j,k)) &
                         + d6(3)*(rhov_n(i+3,j,k)-rhov_n(i-2,j,k))
             Drhow(i,j,k)= d6(1)*(rhow_n(i+1,j,k)-rhow_n(i  ,j,k)) &
                         + d6(2)*(rhow_n(i+2,j,k)-rhow_n(i-1,j,k)) &
                         + d6(3)*(rhow_n(i+3,j,k)-rhow_n(i-2,j,k))
             Drhoe(i,j,k)= d6(1)*(rhoe_n(i+1,j,k)-rhoe_n(i  ,j,k)) &
                         + d6(2)*(rhoe_n(i+2,j,k)-rhoe_n(i-1,j,k)) &
                         + d6(3)*(rhoe_n(i+3,j,k)-rhoe_n(i-2,j,k))
          enddo
       enddo
       
       ! 3rd-order centered dissipation
       ! ------------------------------
       i=nx-2
       do k=1,nz
          do j=1,ny
             Drho(i,j,k) = d4(1)*(rho_n(i+1,j,k)-rho_n(i  ,j,k)) &
                         + d4(2)*(rho_n(i+2,j,k)-rho_n(i-1,j,k))            
             Drhou(i,j,k)= d4(1)*(rhou_n(i+1,j,k)-rhou_n(i  ,j,k)) &
                         + d4(2)*(rhou_n(i+2,j,k)-rhou_n(i-1,j,k))
             Drhov(i,j,k)= d4(1)*(rhov_n(i+1,j,k)-rhov_n(i  ,j,k)) &
                         + d4(2)*(rhov_n(i+2,j,k)-rhov_n(i-1,j,k))
             Drhow(i,j,k)= d4(1)*(rhow_n(i+1,j,k)-rhow_n(i  ,j,k)) &
                         + d4(2)*(rhow_n(i+2,j,k)-rhow_n(i-1,j,k))
             Drhoe(i,j,k)= d4(1)*(rhoe_n(i+1,j,k)-rhoe_n(i  ,j,k)) &
                         + d4(2)*(rhoe_n(i+2,j,k)-rhoe_n(i-1,j,k))
          enddo
       enddo
       
       ! 1st-order centered dissipation
       ! ------------------------------
       i=nx-1
       do k=1,nz
          do j=1,ny
             Drho(i,j,k) = d2(1)*(rho_n(i+1,j,k)-rho_n(i,j,k))            
             Drhou(i,j,k)= d2(1)*(rhou_n(i+1,j,k)-rhou_n(i,j,k))
             Drhov(i,j,k)= d2(1)*(rhov_n(i+1,j,k)-rhov_n(i,j,k))
             Drhow(i,j,k)= d2(1)*(rhow_n(i+1,j,k)-rhow_n(i,j,k))
             Drhoe(i,j,k)= d2(1)*(rhoe_n(i+1,j,k)-rhoe_n(i,j,k))
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

    ! Artificial viscosity in y-direction [compute Drho at j+1/2]
    ! ===================================

    if (is_boundary(2,1)) then
       ! 1st-order centered dissipation
       ! ------------------------------
       j=1
       do k=1,nz
          do i=1,nx
             Drho(i,j,k) = d2(1)*(rho_n(i,j+1,k)-rho_n(i,j,k))
             Drhou(i,j,k)= d2(1)*(rhou_n(i,j+1,k)-rhou_n(i,j,k))
             Drhov(i,j,k)= d2(1)*(rhov_n(i,j+1,k)-rhov_n(i,j,k))
             Drhow(i,j,k)= d2(1)*(rhow_n(i,j+1,k)-rhow_n(i,j,k))
             Drhoe(i,j,k)= d2(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j,k))
          enddo
       enddo
       
       ! 3rd-order centered dissipation
       ! ------------------------------
       j=2
       do k=1,nz
          do i=1,nx
             Drho(i,j,k) = d4(1)*(rho_n(i,j+1,k)-rho_n(i,j  ,k)) &
                         + d4(2)*(rho_n(i,j+2,k)-rho_n(i,j-1,k))
             Drhou(i,j,k)= d4(1)*(rhou_n(i,j+1,k)-rhou_n(i,j  ,k)) &
                         + d4(2)*(rhou_n(i,j+2,k)-rhou_n(i,j-1,k))
             Drhov(i,j,k)= d4(1)*(rhov_n(i,j+1,k)-rhov_n(i,j  ,k)) &
                         + d4(2)*(rhov_n(i,j+2,k)-rhov_n(i,j-1,k))
             Drhow(i,j,k)= d4(1)*(rhow_n(i,j+1,k)-rhow_n(i,j  ,k)) &
                         + d4(2)*(rhow_n(i,j+2,k)-rhow_n(i,j-1,k))
             Drhoe(i,j,k)= d4(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j  ,k)) &
                         + d4(2)*(rhoe_n(i,j+2,k)-rhoe_n(i,j-1,k))
          enddo
       enddo
       
       ! 5th-order centered dissipation
       ! ------------------------------
       j=3
       do k=1,nz
          do i=1,nx
             Drho(i,j,k) = d6(1)*(rho_n(i,j+1,k)-rho_n(i,j  ,k)) &
                         + d6(2)*(rho_n(i,j+2,k)-rho_n(i,j-1,k)) &
                         + d6(3)*(rho_n(i,j+3,k)-rho_n(i,j-2,k))
             Drhou(i,j,k)= d6(1)*(rhou_n(i,j+1,k)-rhou_n(i,j  ,k)) &
                         + d6(2)*(rhou_n(i,j+2,k)-rhou_n(i,j-1,k)) &
                         + d6(3)*(rhou_n(i,j+3,k)-rhou_n(i,j-2,k))
             Drhov(i,j,k)= d6(1)*(rhov_n(i,j+1,k)-rhov_n(i,j  ,k)) &
                         + d6(2)*(rhov_n(i,j+2,k)-rhov_n(i,j-1,k)) &
                         + d6(3)*(rhov_n(i,j+3,k)-rhov_n(i,j-2,k))
             Drhow(i,j,k)= d6(1)*(rhow_n(i,j+1,k)-rhow_n(i,j  ,k)) &
                         + d6(2)*(rhow_n(i,j+2,k)-rhow_n(i,j-1,k)) &
                         + d6(3)*(rhow_n(i,j+3,k)-rhow_n(i,j-2,k))
             Drhoe(i,j,k)= d6(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j  ,k)) &
                         + d6(2)*(rhoe_n(i,j+2,k)-rhoe_n(i,j-1,k)) &
                         + d6(3)*(rhoe_n(i,j+3,k)-rhoe_n(i,j-2,k))
          enddo
       enddo
       
       ! 7th-order centered dissipation
       ! ------------------------------
       j=4
       do k=1,nz
          do i=1,nx
             Drho(i,j,k) = d8(1)*(rho_n(i,j+1,k)-rho_n(i,j  ,k)) &
                         + d8(2)*(rho_n(i,j+2,k)-rho_n(i,j-1,k)) &
                         + d8(3)*(rho_n(i,j+3,k)-rho_n(i,j-2,k)) &
                         + d8(4)*(rho_n(i,j+4,k)-rho_n(i,j-3,k))
             Drhou(i,j,k)= d8(1)*(rhou_n(i,j+1,k)-rhou_n(i,j  ,k)) &
                         + d8(2)*(rhou_n(i,j+2,k)-rhou_n(i,j-1,k)) &
                         + d8(3)*(rhou_n(i,j+3,k)-rhou_n(i,j-2,k)) &
                         + d8(4)*(rhou_n(i,j+4,k)-rhou_n(i,j-3,k))
             Drhov(i,j,k)= d8(1)*(rhov_n(i,j+1,k)-rhov_n(i,j  ,k)) &
                         + d8(2)*(rhov_n(i,j+2,k)-rhov_n(i,j-1,k)) &
                         + d8(3)*(rhov_n(i,j+3,k)-rhov_n(i,j-2,k)) &
                         + d8(4)*(rhov_n(i,j+4,k)-rhov_n(i,j-3,k))
             Drhow(i,j,k)= d8(1)*(rhow_n(i,j+1,k)-rhow_n(i,j  ,k)) &
                         + d8(2)*(rhow_n(i,j+2,k)-rhow_n(i,j-1,k)) &
                         + d8(3)*(rhow_n(i,j+3,k)-rhow_n(i,j-2,k)) &
                         + d8(4)*(rhow_n(i,j+4,k)-rhow_n(i,j-3,k))
             Drhoe(i,j,k)= d8(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j  ,k)) &
                         + d8(2)*(rhoe_n(i,j+2,k)-rhoe_n(i,j-1,k)) &
                         + d8(3)*(rhoe_n(i,j+3,k)-rhoe_n(i,j-2,k)) &
                         + d8(4)*(rhoe_n(i,j+4,k)-rhoe_n(i,j-3,k))
          enddo
       enddo
    endif
    
    ! Interior points: 9th-order centered dissipation
    ! -----------------------------------------------
    do k=1,nz
       do j=ndy_di,nfy_di
          do i=1,nx
             Drho(i,j,k) = d10(1)*(rho_n(i,j+1,k)-rho_n(i,j  ,k)) &
                         + d10(2)*(rho_n(i,j+2,k)-rho_n(i,j-1,k)) &
                         + d10(3)*(rho_n(i,j+3,k)-rho_n(i,j-2,k)) &
                         + d10(4)*(rho_n(i,j+4,k)-rho_n(i,j-3,k)) &
                         + d10(5)*(rho_n(i,j+5,k)-rho_n(i,j-4,k))

             Drhou(i,j,k)= d10(1)*(rhou_n(i,j+1,k)-rhou_n(i,j  ,k)) &
                         + d10(2)*(rhou_n(i,j+2,k)-rhou_n(i,j-1,k)) &
                         + d10(3)*(rhou_n(i,j+3,k)-rhou_n(i,j-2,k)) &
                         + d10(4)*(rhou_n(i,j+4,k)-rhou_n(i,j-3,k)) &
                         + d10(5)*(rhou_n(i,j+5,k)-rhou_n(i,j-4,k))

             Drhov(i,j,k)= d10(1)*(rhov_n(i,j+1,k)-rhov_n(i,j  ,k)) &
                         + d10(2)*(rhov_n(i,j+2,k)-rhov_n(i,j-1,k)) &
                         + d10(3)*(rhov_n(i,j+3,k)-rhov_n(i,j-2,k)) &
                         + d10(4)*(rhov_n(i,j+4,k)-rhov_n(i,j-3,k)) &
                         + d10(5)*(rhov_n(i,j+5,k)-rhov_n(i,j-4,k))

             Drhow(i,j,k)= d10(1)*(rhow_n(i,j+1,k)-rhow_n(i,j  ,k)) &
                         + d10(2)*(rhow_n(i,j+2,k)-rhow_n(i,j-1,k)) &
                         + d10(3)*(rhow_n(i,j+3,k)-rhow_n(i,j-2,k)) &
                         + d10(4)*(rhow_n(i,j+4,k)-rhow_n(i,j-3,k)) &
                         + d10(5)*(rhow_n(i,j+5,k)-rhow_n(i,j-4,k))

             Drhoe(i,j,k)= d10(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j  ,k)) &
                         + d10(2)*(rhoe_n(i,j+2,k)-rhoe_n(i,j-1,k)) &
                         + d10(3)*(rhoe_n(i,j+3,k)-rhoe_n(i,j-2,k)) &
                         + d10(4)*(rhoe_n(i,j+4,k)-rhoe_n(i,j-3,k)) &
                         + d10(5)*(rhoe_n(i,j+5,k)-rhoe_n(i,j-4,k))
          enddo
       enddo
    enddo

    if (is_boundary(2,2)) then
       ! 7th-order centered dissipation
       ! ------------------------------
       j=ny-4
       do k=1,nz
          do i=1,nx
             Drho(i,j,k) = d8(1)*(rho_n(i,j+1,k)-rho_n(i,j  ,k)) &
                         + d8(2)*(rho_n(i,j+2,k)-rho_n(i,j-1,k)) &
                         + d8(3)*(rho_n(i,j+3,k)-rho_n(i,j-2,k)) &
                         + d8(4)*(rho_n(i,j+4,k)-rho_n(i,j-3,k))
             Drhou(i,j,k)= d8(1)*(rhou_n(i,j+1,k)-rhou_n(i,j  ,k)) &
                         + d8(2)*(rhou_n(i,j+2,k)-rhou_n(i,j-1,k)) &
                         + d8(3)*(rhou_n(i,j+3,k)-rhou_n(i,j-2,k)) &
                         + d8(4)*(rhou_n(i,j+4,k)-rhou_n(i,j-3,k))
             Drhov(i,j,k)= d8(1)*(rhov_n(i,j+1,k)-rhov_n(i,j  ,k)) &
                         + d8(2)*(rhov_n(i,j+2,k)-rhov_n(i,j-1,k)) &
                         + d8(3)*(rhov_n(i,j+3,k)-rhov_n(i,j-2,k)) &
                         + d8(4)*(rhov_n(i,j+4,k)-rhov_n(i,j-3,k))
             Drhow(i,j,k)= d8(1)*(rhow_n(i,j+1,k)-rhow_n(i,j  ,k)) &
                         + d8(2)*(rhow_n(i,j+2,k)-rhow_n(i,j-1,k)) &
                         + d8(3)*(rhow_n(i,j+3,k)-rhow_n(i,j-2,k)) &
                         + d8(4)*(rhow_n(i,j+4,k)-rhow_n(i,j-3,k))
             Drhoe(i,j,k)= d8(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j  ,k)) &
                         + d8(2)*(rhoe_n(i,j+2,k)-rhoe_n(i,j-1,k)) &
                         + d8(3)*(rhoe_n(i,j+3,k)-rhoe_n(i,j-2,k)) &
                         + d8(4)*(rhoe_n(i,j+4,k)-rhoe_n(i,j-3,k))
          enddo
       enddo
       
       ! 5th-order centered dissipation
       ! ------------------------------
       j=ny-3
       do k=1,nz
          do i=1,nx
             Drho(i,j,k) = d6(1)*(rho_n(i,j+1,k)-rho_n(i,j  ,k)) &
                         + d6(2)*(rho_n(i,j+2,k)-rho_n(i,j-1,k)) &
                         + d6(3)*(rho_n(i,j+3,k)-rho_n(i,j-2,k))
             Drhou(i,j,k)= d6(1)*(rhou_n(i,j+1,k)-rhou_n(i,j  ,k)) &
                         + d6(2)*(rhou_n(i,j+2,k)-rhou_n(i,j-1,k)) &
                         + d6(3)*(rhou_n(i,j+3,k)-rhou_n(i,j-2,k))
             Drhov(i,j,k)= d6(1)*(rhov_n(i,j+1,k)-rhov_n(i,j  ,k)) &
                         + d6(2)*(rhov_n(i,j+2,k)-rhov_n(i,j-1,k)) &
                         + d6(3)*(rhov_n(i,j+3,k)-rhov_n(i,j-2,k))
             Drhow(i,j,k)= d6(1)*(rhow_n(i,j+1,k)-rhow_n(i,j  ,k)) &
                         + d6(2)*(rhow_n(i,j+2,k)-rhow_n(i,j-1,k)) &
                         + d6(3)*(rhow_n(i,j+3,k)-rhow_n(i,j-2,k))
             Drhoe(i,j,k)= d6(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j  ,k)) &
                         + d6(2)*(rhoe_n(i,j+2,k)-rhoe_n(i,j-1,k)) &
                         + d6(3)*(rhoe_n(i,j+3,k)-rhoe_n(i,j-2,k))
          enddo
       enddo
       
       ! 3rd-order centered dissipation
       ! ------------------------------
       j=ny-2
       do k=1,nz
          do i=1,nx
             Drho(i,j,k) = d4(1)*(rho_n(i,j+1,k)-rho_n(i,j  ,k)) &
                         + d4(2)*(rho_n(i,j+2,k)-rho_n(i,j-1,k))
             Drhou(i,j,k)= d4(1)*(rhou_n(i,j+1,k)-rhou_n(i,j  ,k)) &
                         + d4(2)*(rhou_n(i,j+2,k)-rhou_n(i,j-1,k))
             Drhov(i,j,k)= d4(1)*(rhov_n(i,j+1,k)-rhov_n(i,j  ,k)) &
                         + d4(2)*(rhov_n(i,j+2,k)-rhov_n(i,j-1,k))
             Drhow(i,j,k)= d4(1)*(rhow_n(i,j+1,k)-rhow_n(i,j  ,k)) &
                         + d4(2)*(rhow_n(i,j+2,k)-rhow_n(i,j-1,k))
             Drhoe(i,j,k)= d4(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j  ,k)) &
                         + d4(2)*(rhoe_n(i,j+2,k)-rhoe_n(i,j-1,k))
          enddo
       enddo
       
       ! 1st-order centered dissipation
       ! ------------------------------
       j=ny-1
       do k=1,nz
          do i=1,nx
             Drho(i,j,k) = d2(1)*(rho_n(i,j+1,k)-rho_n(i,j,k))
             Drhou(i,j,k)= d2(1)*(rhou_n(i,j+1,k)-rhou_n(i,j,k))
             Drhov(i,j,k)= d2(1)*(rhov_n(i,j+1,k)-rhov_n(i,j,k))
             Drhow(i,j,k)= d2(1)*(rhow_n(i,j+1,k)-rhow_n(i,j,k))
             Drhoe(i,j,k)= d2(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j,k))
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

       ! Artificial viscosity in z-direction [compute Drho at k+1/2]
       ! ===================================
       if (is_boundary(3,1)) then
          ! 1st-order centered dissipation
          ! ------------------------------
          k=1
          do j=1,ny
             do i=1,nx
                Drho(i,j,k) = d2(1)*(rho_n(i,j,k+1)-rho_n(i,j,k))
                Drhou(i,j,k)= d2(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k))
                Drhov(i,j,k)= d2(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k))
                Drhow(i,j,k)= d2(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k))
                Drhoe(i,j,k)= d2(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k))
             enddo
          enddo
          
          ! 3rd-order centered dissipation
          ! ------------------------------
          k=2
          do j=1,ny
             do i=1,nx
                Drho(i,j,k) = d4(1)*(rho_n(i,j,k+1)-rho_n(i,j,k  )) &
                            + d4(2)*(rho_n(i,j,k+2)-rho_n(i,j,k-1))
                Drhou(i,j,k)= d4(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k  )) &
                            + d4(2)*(rhou_n(i,j,k+2)-rhou_n(i,j,k-1))
                Drhov(i,j,k)= d4(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k  )) &
                            + d4(2)*(rhov_n(i,j,k+2)-rhov_n(i,j,k-1))
                Drhow(i,j,k)= d4(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k  )) &
                            + d4(2)*(rhow_n(i,j,k+2)-rhow_n(i,j,k-1))
                Drhoe(i,j,k)= d4(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k  )) &
                            + d4(2)*(rhoe_n(i,j,k+2)-rhoe_n(i,j,k-1))
             enddo
          enddo
          
          ! 5th-order centered dissipation
          ! ------------------------------
          k=3
          do j=1,ny
             do i=1,nx
                Drho(i,j,k) = d6(1)*(rho_n(i,j,k+1)-rho_n(i,j,k  )) &
                            + d6(2)*(rho_n(i,j,k+2)-rho_n(i,j,k-1)) &
                            + d6(3)*(rho_n(i,j,k+3)-rho_n(i,j,k-2))
                Drhou(i,j,k)= d6(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k  )) &
                            + d6(2)*(rhou_n(i,j,k+2)-rhou_n(i,j,k-1)) &
                            + d6(3)*(rhou_n(i,j,k+3)-rhou_n(i,j,k-2))
                Drhov(i,j,k)= d6(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k  )) &
                            + d6(2)*(rhov_n(i,j,k+2)-rhov_n(i,j,k-1)) &
                            + d6(3)*(rhov_n(i,j,k+3)-rhov_n(i,j,k-2))
                Drhow(i,j,k)= d6(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k  )) &
                            + d6(2)*(rhow_n(i,j,k+2)-rhow_n(i,j,k-1)) &
                            + d6(3)*(rhow_n(i,j,k+3)-rhow_n(i,j,k-2))
                Drhoe(i,j,k)= d6(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k  )) &
                            + d6(2)*(rhoe_n(i,j,k+2)-rhoe_n(i,j,k-1)) &
                            + d6(3)*(rhoe_n(i,j,k+3)-rhoe_n(i,j,k-2))
             enddo
          enddo
          
          ! 7th-order centered dissipation
          ! ------------------------------
          k=4
          do j=1,ny
             do i=1,nx
                Drho(i,j,k) = d8(1)*(rho_n(i,j,k+1)-rho_n(i,j,k  )) &
                            + d8(2)*(rho_n(i,j,k+2)-rho_n(i,j,k-1)) &
                            + d8(3)*(rho_n(i,j,k+3)-rho_n(i,j,k-2)) &
                            + d8(4)*(rho_n(i,j,k+4)-rho_n(i,j,k-3))
                Drhou(i,j,k)= d8(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k  )) &
                            + d8(2)*(rhou_n(i,j,k+2)-rhou_n(i,j,k-1)) &
                            + d8(3)*(rhou_n(i,j,k+3)-rhou_n(i,j,k-2)) &
                            + d8(4)*(rhou_n(i,j,k+4)-rhou_n(i,j,k-3))
                Drhov(i,j,k)= d8(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k  )) &
                            + d8(2)*(rhov_n(i,j,k+2)-rhov_n(i,j,k-1)) &
                            + d8(3)*(rhov_n(i,j,k+3)-rhov_n(i,j,k-2)) &
                            + d8(4)*(rhov_n(i,j,k+4)-rhov_n(i,j,k-3))
                Drhow(i,j,k)= d8(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k  )) &
                            + d8(2)*(rhow_n(i,j,k+2)-rhow_n(i,j,k-1)) &
                            + d8(3)*(rhow_n(i,j,k+3)-rhow_n(i,j,k-2)) &
                            + d8(4)*(rhow_n(i,j,k+4)-rhow_n(i,j,k-3))
                Drhoe(i,j,k)= d8(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k  )) &
                            + d8(2)*(rhoe_n(i,j,k+2)-rhoe_n(i,j,k-1)) &
                            + d8(3)*(rhoe_n(i,j,k+3)-rhoe_n(i,j,k-2)) &
                            + d8(4)*(rhoe_n(i,j,k+4)-rhoe_n(i,j,k-3))
             enddo
          enddo
       endif
       
       ! Interior points: 9th-order centered dissipation
       ! -----------------------------------------------
       do k=ndz_di,nfz_di
          do j=1,ny
             do i=1,nx
                Drho(i,j,k) = d10(1)*(rho_n(i,j,k+1)-rho_n(i,j,k  )) &
                            + d10(2)*(rho_n(i,j,k+2)-rho_n(i,j,k-1)) &
                            + d10(3)*(rho_n(i,j,k+3)-rho_n(i,j,k-2)) &
                            + d10(4)*(rho_n(i,j,k+4)-rho_n(i,j,k-3)) &
                            + d10(5)*(rho_n(i,j,k+5)-rho_n(i,j,k-4))

                Drhou(i,j,k)= d10(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k  )) &
                            + d10(2)*(rhou_n(i,j,k+2)-rhou_n(i,j,k-1)) &
                            + d10(3)*(rhou_n(i,j,k+3)-rhou_n(i,j,k-2)) &
                            + d10(4)*(rhou_n(i,j,k+4)-rhou_n(i,j,k-3)) &
                            + d10(5)*(rhou_n(i,j,k+5)-rhou_n(i,j,k-4))

                Drhov(i,j,k)= d10(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k  )) &
                            + d10(2)*(rhov_n(i,j,k+2)-rhov_n(i,j,k-1)) &
                            + d10(3)*(rhov_n(i,j,k+3)-rhov_n(i,j,k-2)) &
                            + d10(4)*(rhov_n(i,j,k+4)-rhov_n(i,j,k-3)) &
                            + d10(5)*(rhov_n(i,j,k+5)-rhov_n(i,j,k-4))

                Drhow(i,j,k)= d10(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k  )) &
                            + d10(2)*(rhow_n(i,j,k+2)-rhow_n(i,j,k-1)) &
                            + d10(3)*(rhow_n(i,j,k+3)-rhow_n(i,j,k-2)) &
                            + d10(4)*(rhow_n(i,j,k+4)-rhow_n(i,j,k-3)) &
                            + d10(5)*(rhow_n(i,j,k+5)-rhow_n(i,j,k-4))

                Drhoe(i,j,k)= d10(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k  )) &
                            + d10(2)*(rhoe_n(i,j,k+2)-rhoe_n(i,j,k-1)) &
                            + d10(3)*(rhoe_n(i,j,k+3)-rhoe_n(i,j,k-2)) &
                            + d10(4)*(rhoe_n(i,j,k+4)-rhoe_n(i,j,k-3)) &
                            + d10(5)*(rhoe_n(i,j,k+5)-rhoe_n(i,j,k-4))
             enddo
          enddo
       enddo

       if (is_boundary(3,2)) then
          ! 7th-order centered dissipation
          ! ------------------------------
          k=nz-4
          do j=1,ny
             do i=1,nx
                Drho(i,j,k) = d8(1)*(rho_n(i,j,k+1)-rho_n(i,j,k  )) &
                            + d8(2)*(rho_n(i,j,k+2)-rho_n(i,j,k-1)) &
                            + d8(3)*(rho_n(i,j,k+3)-rho_n(i,j,k-2)) &
                            + d8(4)*(rho_n(i,j,k+4)-rho_n(i,j,k-3))
                Drhou(i,j,k)= d8(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k  )) &
                            + d8(2)*(rhou_n(i,j,k+2)-rhou_n(i,j,k-1)) &
                            + d8(3)*(rhou_n(i,j,k+3)-rhou_n(i,j,k-2)) &
                            + d8(4)*(rhou_n(i,j,k+4)-rhou_n(i,j,k-3))
                Drhov(i,j,k)= d8(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k  )) &
                            + d8(2)*(rhov_n(i,j,k+2)-rhov_n(i,j,k-1)) &
                            + d8(3)*(rhov_n(i,j,k+3)-rhov_n(i,j,k-2)) &
                            + d8(4)*(rhov_n(i,j,k+4)-rhov_n(i,j,k-3))
                Drhow(i,j,k)= d8(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k  )) &
                            + d8(2)*(rhow_n(i,j,k+2)-rhow_n(i,j,k-1)) &
                            + d8(3)*(rhow_n(i,j,k+3)-rhow_n(i,j,k-2)) &
                            + d8(4)*(rhow_n(i,j,k+4)-rhow_n(i,j,k-3))
                Drhoe(i,j,k)= d8(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k  )) &
                            + d8(2)*(rhoe_n(i,j,k+2)-rhoe_n(i,j,k-1)) &
                            + d8(3)*(rhoe_n(i,j,k+3)-rhoe_n(i,j,k-2)) &
                            + d8(4)*(rhoe_n(i,j,k+4)-rhoe_n(i,j,k-3))
             enddo
          enddo
          
          ! 5th-order centered dissipation
          ! ------------------------------
          k=nz-3
          do j=1,ny
             do i=1,nx
                Drho(i,j,k) = d6(1)*(rho_n(i,j,k+1)-rho_n(i,j,k  )) &
                            + d6(2)*(rho_n(i,j,k+2)-rho_n(i,j,k-1)) &
                            + d6(3)*(rho_n(i,j,k+3)-rho_n(i,j,k-2))
                Drhou(i,j,k)= d6(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k  )) &
                            + d6(2)*(rhou_n(i,j,k+2)-rhou_n(i,j,k-1)) &
                            + d6(3)*(rhou_n(i,j,k+3)-rhou_n(i,j,k-2))
                Drhov(i,j,k)= d6(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k  )) &
                            + d6(2)*(rhov_n(i,j,k+2)-rhov_n(i,j,k-1)) &
                            + d6(3)*(rhov_n(i,j,k+3)-rhov_n(i,j,k-2))
                Drhow(i,j,k)= d6(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k  )) &
                            + d6(2)*(rhow_n(i,j,k+2)-rhow_n(i,j,k-1)) &
                            + d6(3)*(rhow_n(i,j,k+3)-rhow_n(i,j,k-2))
                Drhoe(i,j,k)= d6(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k  )) &
                            + d6(2)*(rhoe_n(i,j,k+2)-rhoe_n(i,j,k-1)) &
                            + d6(3)*(rhoe_n(i,j,k+3)-rhoe_n(i,j,k-2))
             enddo
          enddo
          
          ! 3rd-order centered dissipation
          ! ------------------------------
          k=nz-2
          do j=1,ny
             do i=1,nx
                Drho(i,j,k) = d4(1)*(rho_n(i,j,k+1)-rho_n(i,j,k  )) &
                            + d4(2)*(rho_n(i,j,k+2)-rho_n(i,j,k-1))
                Drhou(i,j,k)= d4(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k  )) &
                            + d4(2)*(rhou_n(i,j,k+2)-rhou_n(i,j,k-1))
                Drhov(i,j,k)= d4(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k  )) &
                            + d4(2)*(rhov_n(i,j,k+2)-rhov_n(i,j,k-1))
                Drhow(i,j,k)= d4(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k  )) &
                            + d4(2)*(rhow_n(i,j,k+2)-rhow_n(i,j,k-1))
                Drhoe(i,j,k)= d4(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k  )) &
                            + d4(2)*(rhoe_n(i,j,k+2)-rhoe_n(i,j,k-1))
             enddo
          enddo
          
          ! 1st-order centered dissipation
          ! ------------------------------
          k=nz-1
          do j=1,ny
             do i=1,nx
                Drho(i,j,k) = d2(1)*(rho_n(i,j,k+1)-rho_n(i,j,k))
                Drhou(i,j,k)= d2(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k))
                Drhov(i,j,k)= d2(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k))
                Drhow(i,j,k)= d2(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k))
                Drhoe(i,j,k)= d2(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k))
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
    
  end subroutine artvisc_o9_inc

  !==============================================================================
  subroutine artvisc_o7_inc
  !==============================================================================
    !> Apply artificial viscosity - Conservative version DNC order 7
    !> Modified version XG 06/2022
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: vc1,vc2,vc,srm12,srp12,sr
    real(wp), dimension(0:nx,0:ny,0:nz) :: Drho,Drhou,Drhov,Drhow,Drhoe
    ! ---------------------------------------------------------------------------

    ! Artificial viscosity in x-direction [compute Drho at i+1/2]
    ! ===================================
    
    if (is_boundary(1,1)) then
       ! 1st-order centered dissipation
       ! ------------------------------
       i=1
       do k=1,nz
          do j=1,ny
             Drho(i,j,k) = d2(1)*(rho_n(i+1,j,k) -rho_n(i,j,k))            
             Drhou(i,j,k)= d2(1)*(rhou_n(i+1,j,k)-rhou_n(i,j,k))
             Drhov(i,j,k)= d2(1)*(rhov_n(i+1,j,k)-rhov_n(i,j,k))
             Drhow(i,j,k)= d2(1)*(rhow_n(i+1,j,k)-rhow_n(i,j,k))
             Drhoe(i,j,k)= d2(1)*(rhoe_n(i+1,j,k)-rhoe_n(i,j,k))
          enddo
       enddo
       
       ! 3rd-order centered dissipation
       ! ------------------------------
       i=2
       do k=1,nz
          do j=1,ny
             Drho(i,j,k) = d4(1)*(rho_n(i+1,j,k)-rho_n(i  ,j,k)) &
                         + d4(2)*(rho_n(i+2,j,k)-rho_n(i-1,j,k))
             Drhou(i,j,k)= d4(1)*(rhou_n(i+1,j,k)-rhou_n(i  ,j,k)) &
                         + d4(2)*(rhou_n(i+2,j,k)-rhou_n(i-1,j,k))
             Drhov(i,j,k)= d4(1)*(rhov_n(i+1,j,k)-rhov_n(i  ,j,k)) &
                         + d4(2)*(rhov_n(i+2,j,k)-rhov_n(i-1,j,k))
             Drhow(i,j,k)= d4(1)*(rhow_n(i+1,j,k)-rhow_n(i  ,j,k)) &
                         + d4(2)*(rhow_n(i+2,j,k)-rhow_n(i-1,j,k))
             Drhoe(i,j,k)= d4(1)*(rhoe_n(i+1,j,k)-rhoe_n(i  ,j,k)) &
                         + d4(2)*(rhoe_n(i+2,j,k)-rhoe_n(i-1,j,k))
          enddo
       enddo
       
       ! 5th-order centered dissipation
       ! ------------------------------
       i=3
       do k=1,nz
          do j=1,ny
             Drho(i,j,k) = d6(1)*(rho_n(i+1,j,k)-rho_n(i  ,j,k)) &
                         + d6(2)*(rho_n(i+2,j,k)-rho_n(i-1,j,k)) &
                         + d6(3)*(rho_n(i+3,j,k)-rho_n(i-2,j,k))
             Drhou(i,j,k)= d6(1)*(rhou_n(i+1,j,k)-rhou_n(i  ,j,k)) &
                         + d6(2)*(rhou_n(i+2,j,k)-rhou_n(i-1,j,k)) &
                         + d6(3)*(rhou_n(i+3,j,k)-rhou_n(i-2,j,k))
             Drhov(i,j,k)= d6(1)*(rhov_n(i+1,j,k)-rhov_n(i  ,j,k)) &
                         + d6(2)*(rhov_n(i+2,j,k)-rhov_n(i-1,j,k)) &
                         + d6(3)*(rhov_n(i+3,j,k)-rhov_n(i-2,j,k))
             Drhow(i,j,k)= d6(1)*(rhow_n(i+1,j,k)-rhow_n(i  ,j,k)) &
                         + d6(2)*(rhow_n(i+2,j,k)-rhow_n(i-1,j,k)) &
                         + d6(3)*(rhow_n(i+3,j,k)-rhow_n(i-2,j,k))
             Drhoe(i,j,k)= d6(1)*(rhoe_n(i+1,j,k)-rhoe_n(i  ,j,k)) &
                         + d6(2)*(rhoe_n(i+2,j,k)-rhoe_n(i-1,j,k)) &
                         + d6(3)*(rhoe_n(i+3,j,k)-rhoe_n(i-2,j,k))
          enddo
       enddo
    endif

    ! Interior points: 7th-order centered dissipation
    ! -----------------------------------------------
    do k=1,nz
       do j=1,ny
          do i=ndx_di,nfx_di
             Drho(i,j,k) = d8(1)*(rho_n(i+1,j,k)-rho_n(i  ,j,k)) &
                         + d8(2)*(rho_n(i+2,j,k)-rho_n(i-1,j,k)) &
                         + d8(3)*(rho_n(i+3,j,k)-rho_n(i-2,j,k)) &
                         + d8(4)*(rho_n(i+4,j,k)-rho_n(i-3,j,k))
             
             Drhou(i,j,k)= d8(1)*(rhou_n(i+1,j,k)-rhou_n(i  ,j,k)) &
                         + d8(2)*(rhou_n(i+2,j,k)-rhou_n(i-1,j,k)) &
                         + d8(3)*(rhou_n(i+3,j,k)-rhou_n(i-2,j,k)) &
                         + d8(4)*(rhou_n(i+4,j,k)-rhou_n(i-3,j,k))

             Drhov(i,j,k)= d8(1)*(rhov_n(i+1,j,k)-rhov_n(i  ,j,k)) &
                         + d8(2)*(rhov_n(i+2,j,k)-rhov_n(i-1,j,k)) &
                         + d8(3)*(rhov_n(i+3,j,k)-rhov_n(i-2,j,k)) &
                         + d8(4)*(rhov_n(i+4,j,k)-rhov_n(i-3,j,k))

             Drhow(i,j,k)= d8(1)*(rhow_n(i+1,j,k)-rhow_n(i  ,j,k)) &
                         + d8(2)*(rhow_n(i+2,j,k)-rhow_n(i-1,j,k)) &
                         + d8(3)*(rhow_n(i+3,j,k)-rhow_n(i-2,j,k)) &
                         + d8(4)*(rhow_n(i+4,j,k)-rhow_n(i-3,j,k))

             Drhoe(i,j,k)= d8(1)*(rhoe_n(i+1,j,k)-rhoe_n(i  ,j,k)) &
                         + d8(2)*(rhoe_n(i+2,j,k)-rhoe_n(i-1,j,k)) &
                         + d8(3)*(rhoe_n(i+3,j,k)-rhoe_n(i-2,j,k)) &
                         + d8(4)*(rhoe_n(i+4,j,k)-rhoe_n(i-3,j,k))
          enddo
       enddo
    enddo

    if (is_boundary(1,2)) then
       ! 5th-order centered dissipation
       ! ------------------------------
       i=nx-3
       do k=1,nz
          do j=1,ny
             Drho(i,j,k) = d6(1)*(rho_n(i+1,j,k)-rho_n(i  ,j,k)) &
                         + d6(2)*(rho_n(i+2,j,k)-rho_n(i-1,j,k)) &
                         + d6(3)*(rho_n(i+3,j,k)-rho_n(i-2,j,k))    
             Drhou(i,j,k)= d6(1)*(rhou_n(i+1,j,k)-rhou_n(i  ,j,k)) &
                         + d6(2)*(rhou_n(i+2,j,k)-rhou_n(i-1,j,k)) &
                         + d6(3)*(rhou_n(i+3,j,k)-rhou_n(i-2,j,k))
             Drhov(i,j,k)= d6(1)*(rhov_n(i+1,j,k)-rhov_n(i  ,j,k)) &
                         + d6(2)*(rhov_n(i+2,j,k)-rhov_n(i-1,j,k)) &
                         + d6(3)*(rhov_n(i+3,j,k)-rhov_n(i-2,j,k))
             Drhow(i,j,k)= d6(1)*(rhow_n(i+1,j,k)-rhow_n(i  ,j,k)) &
                         + d6(2)*(rhow_n(i+2,j,k)-rhow_n(i-1,j,k)) &
                         + d6(3)*(rhow_n(i+3,j,k)-rhow_n(i-2,j,k))
             Drhoe(i,j,k)= d6(1)*(rhoe_n(i+1,j,k)-rhoe_n(i  ,j,k)) &
                         + d6(2)*(rhoe_n(i+2,j,k)-rhoe_n(i-1,j,k)) &
                         + d6(3)*(rhoe_n(i+3,j,k)-rhoe_n(i-2,j,k))
          enddo
       enddo
       
       ! 3rd-order centered dissipation
       ! ------------------------------
       i=nx-2
       do k=1,nz
          do j=1,ny
             Drho(i,j,k) = d4(1)*(rho_n(i+1,j,k)-rho_n(i  ,j,k)) &
                         + d4(2)*(rho_n(i+2,j,k)-rho_n(i-1,j,k))            
             Drhou(i,j,k)= d4(1)*(rhou_n(i+1,j,k)-rhou_n(i  ,j,k)) &
                         + d4(2)*(rhou_n(i+2,j,k)-rhou_n(i-1,j,k))
             Drhov(i,j,k)= d4(1)*(rhov_n(i+1,j,k)-rhov_n(i  ,j,k)) &
                         + d4(2)*(rhov_n(i+2,j,k)-rhov_n(i-1,j,k))
             Drhow(i,j,k)= d4(1)*(rhow_n(i+1,j,k)-rhow_n(i  ,j,k)) &
                         + d4(2)*(rhow_n(i+2,j,k)-rhow_n(i-1,j,k))
             Drhoe(i,j,k)= d4(1)*(rhoe_n(i+1,j,k)-rhoe_n(i  ,j,k)) &
                         + d4(2)*(rhoe_n(i+2,j,k)-rhoe_n(i-1,j,k))
          enddo
       enddo
       
       ! 1st-order centered dissipation
       ! ------------------------------
       i=nx-1
       do k=1,nz
          do j=1,ny
             Drho(i,j,k) = d2(1)*(rho_n(i+1,j,k)-rho_n(i,j,k))            
             Drhou(i,j,k)= d2(1)*(rhou_n(i+1,j,k)-rhou_n(i,j,k))
             Drhov(i,j,k)= d2(1)*(rhov_n(i+1,j,k)-rhov_n(i,j,k))
             Drhow(i,j,k)= d2(1)*(rhow_n(i+1,j,k)-rhow_n(i,j,k))
             Drhoe(i,j,k)= d2(1)*(rhoe_n(i+1,j,k)-rhoe_n(i,j,k))
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
                vc=(uu(i,j,k)*y_eta(i,j)-vv(i,j,k)*x_eta(i,j))*ijacob(i,j)
                vc2=(uu(i+1,j,k)*y_eta(i+1,j)-vv(i+1,j,k)*x_eta(i+1,j))*ijacob(i+1,j)
                vc1=(uu(i-1,j,k)*y_eta(i-1,j)-vv(i-1,j,k)*x_eta(i-1,j))*ijacob(i-1,j)
                ! compute spectral radius at i+1/2 (srp12) and i-1/2 (srm12)
                sr=abs(vc)+c_(i,j,k)*g_ksi(i,j)
                srp12=0.5_wp*(abs(vc2)+c_(i+1,j,k)*g_ksi(i+1,j)+sr)
                srm12=0.5_wp*(abs(vc1)+c_(i-1,j,k)*g_ksi(i-1,j)+sr)
                ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dksi
                Krho(i,j,k) = Krho(i,j,k) +(srp12*Drho(i,j,k) -srm12*Drho(i-1,j,k))
                Krhou(i,j,k)= Krhou(i,j,k)+(srp12*Drhou(i,j,k)-srm12*Drhou(i-1,j,k))
                Krhov(i,j,k)= Krhov(i,j,k)+(srp12*Drhov(i,j,k)-srm12*Drhov(i-1,j,k))
                Krhow(i,j,k)= Krhow(i,j,k)+(srp12*Drhow(i,j,k)-srm12*Drhow(i-1,j,k))
                Krhoe(i,j,k)= Krhoe(i,j,k)+(srp12*Drhoe(i,j,k)-srm12*Drhoe(i-1,j,k))
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
             enddo
          enddo
       enddo
    endif

    ! Artificial viscosity in y-direction [compute Drho at j+1/2]
    ! ===================================
    if (is_boundary(2,1)) then
       ! 1st-order centered dissipation
       ! ------------------------------
       j=1
       do k=1,nz
          do i=1,nx
             Drho(i,j,k) = d2(1)*(rho_n(i,j+1,k)-rho_n(i,j,k))
             Drhou(i,j,k)= d2(1)*(rhou_n(i,j+1,k)-rhou_n(i,j,k))
             Drhov(i,j,k)= d2(1)*(rhov_n(i,j+1,k)-rhov_n(i,j,k))
             Drhow(i,j,k)= d2(1)*(rhow_n(i,j+1,k)-rhow_n(i,j,k))
             Drhoe(i,j,k)= d2(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j,k))
          enddo
       enddo
       
       ! 3rd-order centered dissipation
       ! ------------------------------
       j=2
       do k=1,nz
          do i=1,nx
             Drho(i,j,k) = d4(1)*(rho_n(i,j+1,k)-rho_n(i,j  ,k)) &
                         + d4(2)*(rho_n(i,j+2,k)-rho_n(i,j-1,k))
             Drhou(i,j,k)= d4(1)*(rhou_n(i,j+1,k)-rhou_n(i,j  ,k)) &
                         + d4(2)*(rhou_n(i,j+2,k)-rhou_n(i,j-1,k))
             Drhov(i,j,k)= d4(1)*(rhov_n(i,j+1,k)-rhov_n(i,j  ,k)) &
                         + d4(2)*(rhov_n(i,j+2,k)-rhov_n(i,j-1,k))
             Drhow(i,j,k)= d4(1)*(rhow_n(i,j+1,k)-rhow_n(i,j  ,k)) &
                         + d4(2)*(rhow_n(i,j+2,k)-rhow_n(i,j-1,k))
             Drhoe(i,j,k)= d4(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j  ,k)) &
                         + d4(2)*(rhoe_n(i,j+2,k)-rhoe_n(i,j-1,k))
          enddo
       enddo
       
       ! 5th-order centered dissipation
       ! ------------------------------
       j=3
       do k=1,nz
          do i=1,nx
             Drho(i,j,k) = d6(1)*(rho_n(i,j+1,k)-rho_n(i,j  ,k)) &
                         + d6(2)*(rho_n(i,j+2,k)-rho_n(i,j-1,k)) &
                         + d6(3)*(rho_n(i,j+3,k)-rho_n(i,j-2,k))
             Drhou(i,j,k)= d6(1)*(rhou_n(i,j+1,k)-rhou_n(i,j  ,k)) &
                         + d6(2)*(rhou_n(i,j+2,k)-rhou_n(i,j-1,k)) &
                         + d6(3)*(rhou_n(i,j+3,k)-rhou_n(i,j-2,k))
             Drhov(i,j,k)= d6(1)*(rhov_n(i,j+1,k)-rhov_n(i,j  ,k)) &
                         + d6(2)*(rhov_n(i,j+2,k)-rhov_n(i,j-1,k)) &
                         + d6(3)*(rhov_n(i,j+3,k)-rhov_n(i,j-2,k))
             Drhow(i,j,k)= d6(1)*(rhow_n(i,j+1,k)-rhow_n(i,j  ,k)) &
                         + d6(2)*(rhow_n(i,j+2,k)-rhow_n(i,j-1,k)) &
                         + d6(3)*(rhow_n(i,j+3,k)-rhow_n(i,j-2,k))
             Drhoe(i,j,k)= d6(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j  ,k)) &
                         + d6(2)*(rhoe_n(i,j+2,k)-rhoe_n(i,j-1,k)) &
                         + d6(3)*(rhoe_n(i,j+3,k)-rhoe_n(i,j-2,k))
          enddo
       enddo
    endif
    
    ! Interior points: 7th-order centered dissipation
    ! -----------------------------------------------
    do k=1,nz
       do j=ndy_di,nfy_di
          do i=1,nx
             Drho(i,j,k) = d8(1)*(rho_n(i,j+1,k)-rho_n(i,j  ,k)) &
                         + d8(2)*(rho_n(i,j+2,k)-rho_n(i,j-1,k)) &
                         + d8(3)*(rho_n(i,j+3,k)-rho_n(i,j-2,k)) &
                         + d8(4)*(rho_n(i,j+4,k)-rho_n(i,j-3,k))

             Drhou(i,j,k)= d8(1)*(rhou_n(i,j+1,k)-rhou_n(i,j  ,k)) &
                         + d8(2)*(rhou_n(i,j+2,k)-rhou_n(i,j-1,k)) &
                         + d8(3)*(rhou_n(i,j+3,k)-rhou_n(i,j-2,k)) &
                         + d8(4)*(rhou_n(i,j+4,k)-rhou_n(i,j-3,k))

             Drhov(i,j,k)= d8(1)*(rhov_n(i,j+1,k)-rhov_n(i,j  ,k)) &
                         + d8(2)*(rhov_n(i,j+2,k)-rhov_n(i,j-1,k)) &
                         + d8(3)*(rhov_n(i,j+3,k)-rhov_n(i,j-2,k)) &
                         + d8(4)*(rhov_n(i,j+4,k)-rhov_n(i,j-3,k))

             Drhow(i,j,k)= d8(1)*(rhow_n(i,j+1,k)-rhow_n(i,j  ,k)) &
                         + d8(2)*(rhow_n(i,j+2,k)-rhow_n(i,j-1,k)) &
                         + d8(3)*(rhow_n(i,j+3,k)-rhow_n(i,j-2,k)) &
                         + d8(4)*(rhow_n(i,j+4,k)-rhow_n(i,j-3,k))

             Drhoe(i,j,k)= d8(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j  ,k)) &
                         + d8(2)*(rhoe_n(i,j+2,k)-rhoe_n(i,j-1,k)) &
                         + d8(3)*(rhoe_n(i,j+3,k)-rhoe_n(i,j-2,k)) &
                         + d8(4)*(rhoe_n(i,j+4,k)-rhoe_n(i,j-3,k))
          enddo
       enddo
    enddo

    if (is_boundary(2,2)) then
       ! 5th-order centered dissipation
       ! ------------------------------
       j=ny-3
       do k=1,nz
          do i=1,nx
             Drho(i,j,k) = d6(1)*(rho_n(i,j+1,k)-rho_n(i,j  ,k)) &
                         + d6(2)*(rho_n(i,j+2,k)-rho_n(i,j-1,k)) &
                         + d6(3)*(rho_n(i,j+3,k)-rho_n(i,j-2,k))
             Drhou(i,j,k)= d6(1)*(rhou_n(i,j+1,k)-rhou_n(i,j  ,k)) &
                         + d6(2)*(rhou_n(i,j+2,k)-rhou_n(i,j-1,k)) &
                         + d6(3)*(rhou_n(i,j+3,k)-rhou_n(i,j-2,k))
             Drhov(i,j,k)= d6(1)*(rhov_n(i,j+1,k)-rhov_n(i,j  ,k)) &
                         + d6(2)*(rhov_n(i,j+2,k)-rhov_n(i,j-1,k)) &
                         + d6(3)*(rhov_n(i,j+3,k)-rhov_n(i,j-2,k))
             Drhow(i,j,k)= d6(1)*(rhow_n(i,j+1,k)-rhow_n(i,j  ,k)) &
                         + d6(2)*(rhow_n(i,j+2,k)-rhow_n(i,j-1,k)) &
                         + d6(3)*(rhow_n(i,j+3,k)-rhow_n(i,j-2,k))
             Drhoe(i,j,k)= d6(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j  ,k)) &
                         + d6(2)*(rhoe_n(i,j+2,k)-rhoe_n(i,j-1,k)) &
                         + d6(3)*(rhoe_n(i,j+3,k)-rhoe_n(i,j-2,k))
          enddo
       enddo
       
       ! 3rd-order centered dissipation
       ! ------------------------------
       j=ny-2
       do k=1,nz
          do i=1,nx
             Drho(i,j,k) = d4(1)*(rho_n(i,j+1,k)-rho_n(i,j  ,k)) &
                         + d4(2)*(rho_n(i,j+2,k)-rho_n(i,j-1,k))
             Drhou(i,j,k)= d4(1)*(rhou_n(i,j+1,k)-rhou_n(i,j  ,k)) &
                         + d4(2)*(rhou_n(i,j+2,k)-rhou_n(i,j-1,k))
             Drhov(i,j,k)= d4(1)*(rhov_n(i,j+1,k)-rhov_n(i,j  ,k)) &
                         + d4(2)*(rhov_n(i,j+2,k)-rhov_n(i,j-1,k))
             Drhow(i,j,k)= d4(1)*(rhow_n(i,j+1,k)-rhow_n(i,j  ,k)) &
                         + d4(2)*(rhow_n(i,j+2,k)-rhow_n(i,j-1,k))
             Drhoe(i,j,k)= d4(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j  ,k)) &
                         + d4(2)*(rhoe_n(i,j+2,k)-rhoe_n(i,j-1,k))
          enddo
       enddo
       
       ! 1st-order centered dissipation
       ! ------------------------------
       j=ny-1
       do k=1,nz
          do i=1,nx
             Drho(i,j,k) = d2(1)*(rho_n(i,j+1,k)-rho_n(i,j,k))
             Drhou(i,j,k)= d2(1)*(rhou_n(i,j+1,k)-rhou_n(i,j,k))
             Drhov(i,j,k)= d2(1)*(rhov_n(i,j+1,k)-rhov_n(i,j,k))
             Drhow(i,j,k)= d2(1)*(rhow_n(i,j+1,k)-rhow_n(i,j,k))
             Drhoe(i,j,k)= d2(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j,k))
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
                vc=(vv(i,j,k)*x_ksi(i,j)-uu(i,j,k)*y_ksi(i,j))*ijacob(i,j)
                vc2=(vv(i,j+1,k)*x_ksi(i,j+1)-uu(i,j+1,k)*y_ksi(i,j+1))*ijacob(i,j+1)
                vc1=(vv(i,j-1,k)*x_ksi(i,j-1)-uu(i,j-1,k)*y_ksi(i,j-1))*ijacob(i,j-1)
                ! compute spectral radius at j+1/2 (srp12) and j-1/2 (srm12)
                sr=abs(vc)+c_(i,j,k)*g_eta(i,j)
                srp12=0.5_wp*(abs(vc2)+c_(i,j+1,k)*g_eta(i,j+1)+sr)
                srm12=0.5_wp*(abs(vc1)+c_(i,j-1,k)*g_eta(i,j-1)+sr)
                ! compute dissipation term: delta_j/delta_y=(D_{j+1/2}-D_{j-1/2})/dy
                Krho(i,j,k) = Krho(i,j,k) +(srp12*Drho(i,j,k) -srm12*Drho(i,j-1,k) )
                Krhou(i,j,k)= Krhou(i,j,k)+(srp12*Drhou(i,j,k)-srm12*Drhou(i,j-1,k))
                Krhov(i,j,k)= Krhov(i,j,k)+(srp12*Drhov(i,j,k)-srm12*Drhov(i,j-1,k))
                Krhow(i,j,k)= Krhow(i,j,k)+(srp12*Drhow(i,j,k)-srm12*Drhow(i,j-1,k))
                Krhoe(i,j,k)= Krhoe(i,j,k)+(srp12*Drhoe(i,j,k)-srm12*Drhoe(i,j-1,k))
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
             enddo
          enddo
       enddo
    endif
 
    !*******************
    if (.not.is_2D) then
    !*******************

       ! Artificial viscosity in z-direction [compute Drho at k+1/2]
       ! ===================================
       if (is_boundary(3,1)) then
          ! 1st-order centered dissipation
          ! ------------------------------
          k=1
          do j=1,ny
             do i=1,nx
                Drho(i,j,k) = d2(1)*(rho_n(i,j,k+1)-rho_n(i,j,k))
                Drhou(i,j,k)= d2(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k))
                Drhov(i,j,k)= d2(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k))
                Drhow(i,j,k)= d2(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k))
                Drhoe(i,j,k)= d2(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k))
             enddo
          enddo
          
          ! 3rd-order centered dissipation
          ! ------------------------------
          k=2
          do j=1,ny
             do i=1,nx
                Drho(i,j,k) = d4(1)*(rho_n(i,j,k+1)-rho_n(i,j,k  )) &
                            + d4(2)*(rho_n(i,j,k+2)-rho_n(i,j,k-1))
                Drhou(i,j,k)= d4(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k  )) &
                            + d4(2)*(rhou_n(i,j,k+2)-rhou_n(i,j,k-1))
                Drhov(i,j,k)= d4(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k  )) &
                            + d4(2)*(rhov_n(i,j,k+2)-rhov_n(i,j,k-1))
                Drhow(i,j,k)= d4(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k  )) &
                            + d4(2)*(rhow_n(i,j,k+2)-rhow_n(i,j,k-1))
                Drhoe(i,j,k)= d4(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k  )) &
                            + d4(2)*(rhoe_n(i,j,k+2)-rhoe_n(i,j,k-1))
             enddo
          enddo
          
          ! 5th-order centered dissipation
          ! ------------------------------
          k=3
          do j=1,ny
             do i=1,nx
                Drho(i,j,k) = d6(1)*(rho_n(i,j,k+1)-rho_n(i,j,k  )) &
                            + d6(2)*(rho_n(i,j,k+2)-rho_n(i,j,k-1)) &
                            + d6(3)*(rho_n(i,j,k+3)-rho_n(i,j,k-2))
                Drhou(i,j,k)= d6(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k  )) &
                            + d6(2)*(rhou_n(i,j,k+2)-rhou_n(i,j,k-1)) &
                            + d6(3)*(rhou_n(i,j,k+3)-rhou_n(i,j,k-2))
                Drhov(i,j,k)= d6(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k  )) &
                            + d6(2)*(rhov_n(i,j,k+2)-rhov_n(i,j,k-1)) &
                            + d6(3)*(rhov_n(i,j,k+3)-rhov_n(i,j,k-2))
                Drhow(i,j,k)= d6(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k  )) &
                            + d6(2)*(rhow_n(i,j,k+2)-rhow_n(i,j,k-1)) &
                            + d6(3)*(rhow_n(i,j,k+3)-rhow_n(i,j,k-2))
                Drhoe(i,j,k)= d6(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k  )) &
                            + d6(2)*(rhoe_n(i,j,k+2)-rhoe_n(i,j,k-1)) &
                            + d6(3)*(rhoe_n(i,j,k+3)-rhoe_n(i,j,k-2))
             enddo
          enddo
       endif
       
       ! Interior points: 7th-order centered dissipation
       ! -----------------------------------------------
       do k=ndz_di,nfz_di
          do j=1,ny
             do i=1,nx
                Drho(i,j,k) = d8(1)*(rho_n(i,j,k+1)-rho_n(i,j,k  )) &
                            + d8(2)*(rho_n(i,j,k+2)-rho_n(i,j,k-1)) &
                            + d8(3)*(rho_n(i,j,k+3)-rho_n(i,j,k-2)) &
                            + d8(4)*(rho_n(i,j,k+4)-rho_n(i,j,k-3))

                Drhou(i,j,k)= d8(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k  )) &
                            + d8(2)*(rhou_n(i,j,k+2)-rhou_n(i,j,k-1)) &
                            + d8(3)*(rhou_n(i,j,k+3)-rhou_n(i,j,k-2)) &
                            + d8(4)*(rhou_n(i,j,k+4)-rhou_n(i,j,k-3))

                Drhov(i,j,k)= d8(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k  )) &
                            + d8(2)*(rhov_n(i,j,k+2)-rhov_n(i,j,k-1)) &
                            + d8(3)*(rhov_n(i,j,k+3)-rhov_n(i,j,k-2)) &
                            + d8(4)*(rhov_n(i,j,k+4)-rhov_n(i,j,k-3))

                Drhow(i,j,k)= d8(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k  )) &
                            + d8(2)*(rhow_n(i,j,k+2)-rhow_n(i,j,k-1)) &
                            + d8(3)*(rhow_n(i,j,k+3)-rhow_n(i,j,k-2)) &
                            + d8(4)*(rhow_n(i,j,k+4)-rhow_n(i,j,k-3))

                Drhoe(i,j,k)= d8(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k  )) &
                            + d8(2)*(rhoe_n(i,j,k+2)-rhoe_n(i,j,k-1)) &
                            + d8(3)*(rhoe_n(i,j,k+3)-rhoe_n(i,j,k-2)) &
                            + d8(4)*(rhoe_n(i,j,k+4)-rhoe_n(i,j,k-3))
             enddo
          enddo
       enddo

       if (is_boundary(3,2)) then
          ! 5th-order centered dissipation
          ! ------------------------------
          k=nz-3
          do j=1,ny
             do i=1,nx
                Drho(i,j,k) = d6(1)*(rho_n(i,j,k+1)-rho_n(i,j,k  )) &
                            + d6(2)*(rho_n(i,j,k+2)-rho_n(i,j,k-1)) &
                            + d6(3)*(rho_n(i,j,k+3)-rho_n(i,j,k-2))
                Drhou(i,j,k)= d6(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k  )) &
                            + d6(2)*(rhou_n(i,j,k+2)-rhou_n(i,j,k-1)) &
                            + d6(3)*(rhou_n(i,j,k+3)-rhou_n(i,j,k-2))
                Drhov(i,j,k)= d6(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k  )) &
                            + d6(2)*(rhov_n(i,j,k+2)-rhov_n(i,j,k-1)) &
                            + d6(3)*(rhov_n(i,j,k+3)-rhov_n(i,j,k-2))
                Drhow(i,j,k)= d6(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k  )) &
                            + d6(2)*(rhow_n(i,j,k+2)-rhow_n(i,j,k-1)) &
                            + d6(3)*(rhow_n(i,j,k+3)-rhow_n(i,j,k-2))
                Drhoe(i,j,k)= d6(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k  )) &
                            + d6(2)*(rhoe_n(i,j,k+2)-rhoe_n(i,j,k-1)) &
                            + d6(3)*(rhoe_n(i,j,k+3)-rhoe_n(i,j,k-2))
             enddo
          enddo
          
          ! 3rd-order centered dissipation
          ! ------------------------------
          k=nz-2
          do j=1,ny
             do i=1,nx
                Drho(i,j,k) = d4(1)*(rho_n(i,j,k+1)-rho_n(i,j,k  )) &
                            + d4(2)*(rho_n(i,j,k+2)-rho_n(i,j,k-1))
                Drhou(i,j,k)= d4(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k  )) &
                            + d4(2)*(rhou_n(i,j,k+2)-rhou_n(i,j,k-1))
                Drhov(i,j,k)= d4(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k  )) &
                            + d4(2)*(rhov_n(i,j,k+2)-rhov_n(i,j,k-1))
                Drhow(i,j,k)= d4(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k  )) &
                            + d4(2)*(rhow_n(i,j,k+2)-rhow_n(i,j,k-1))
                Drhoe(i,j,k)= d4(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k  )) &
                            + d4(2)*(rhoe_n(i,j,k+2)-rhoe_n(i,j,k-1))
             enddo
          enddo
          
          ! 1st-order centered dissipation
          ! ------------------------------
          k=nz-1
          do j=1,ny
             do i=1,nx
                Drho(i,j,k) = d2(1)*(rho_n(i,j,k+1)-rho_n(i,j,k))
                Drhou(i,j,k)= d2(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k))
                Drhov(i,j,k)= d2(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k))
                Drhow(i,j,k)= d2(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k))
                Drhoe(i,j,k)= d2(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k))
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
    
  end subroutine artvisc_o7_inc
  
  !==============================================================================
  subroutine artvisc_o5_inc
  !==============================================================================
    !> Apply artificial viscosity - Conservative version DNC order 5
    !> Modified version XG 06/2022
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: vc1,vc2,vc,srm12,srp12,sr
    real(wp), dimension(0:nx,0:ny,0:nz) :: Drho,Drhou,Drhov,Drhow,Drhoe
    ! ---------------------------------------------------------------------------

    ! Artificial viscosity in x-direction [compute Drho at i+1/2]
    ! ===================================
    
    if (is_boundary(1,1)) then
       ! 1st-order centered dissipation
       ! ------------------------------
       i=1
       do k=1,nz
          do j=1,ny
             Drho(i,j,k) = d2(1)*(rho_n(i+1,j,k) -rho_n(i,j,k))            
             Drhou(i,j,k)= d2(1)*(rhou_n(i+1,j,k)-rhou_n(i,j,k))
             Drhov(i,j,k)= d2(1)*(rhov_n(i+1,j,k)-rhov_n(i,j,k))
             Drhow(i,j,k)= d2(1)*(rhow_n(i+1,j,k)-rhow_n(i,j,k))
             Drhoe(i,j,k)= d2(1)*(rhoe_n(i+1,j,k)-rhoe_n(i,j,k))
          enddo
       enddo
       
       ! 3rd-order centered dissipation
       ! ------------------------------
       i=2
       do k=1,nz
          do j=1,ny
             Drho(i,j,k) = d4(1)*(rho_n(i+1,j,k)-rho_n(i  ,j,k)) &
                         + d4(2)*(rho_n(i+2,j,k)-rho_n(i-1,j,k))
             Drhou(i,j,k)= d4(1)*(rhou_n(i+1,j,k)-rhou_n(i  ,j,k)) &
                         + d4(2)*(rhou_n(i+2,j,k)-rhou_n(i-1,j,k))
             Drhov(i,j,k)= d4(1)*(rhov_n(i+1,j,k)-rhov_n(i  ,j,k)) &
                         + d4(2)*(rhov_n(i+2,j,k)-rhov_n(i-1,j,k))
             Drhow(i,j,k)= d4(1)*(rhow_n(i+1,j,k)-rhow_n(i  ,j,k)) &
                         + d4(2)*(rhow_n(i+2,j,k)-rhow_n(i-1,j,k))
             Drhoe(i,j,k)= d4(1)*(rhoe_n(i+1,j,k)-rhoe_n(i  ,j,k)) &
                         + d4(2)*(rhoe_n(i+2,j,k)-rhoe_n(i-1,j,k))
          enddo
       enddo
    endif

    ! Interior points: 5th-order centered dissipation
    ! -----------------------------------------------
    do k=1,nz
       do j=1,ny
          do i=ndx_di,nfx_di
             Drho(i,j,k) = d6(1)*(rho_n(i+1,j,k)-rho_n(i  ,j,k)) &
                         + d6(2)*(rho_n(i+2,j,k)-rho_n(i-1,j,k)) &
                         + d6(3)*(rho_n(i+3,j,k)-rho_n(i-2,j,k))
             
             Drhou(i,j,k)= d6(1)*(rhou_n(i+1,j,k)-rhou_n(i  ,j,k)) &
                         + d6(2)*(rhou_n(i+2,j,k)-rhou_n(i-1,j,k)) &
                         + d6(3)*(rhou_n(i+3,j,k)-rhou_n(i-2,j,k))

             Drhov(i,j,k)= d6(1)*(rhov_n(i+1,j,k)-rhov_n(i  ,j,k)) &
                         + d6(2)*(rhov_n(i+2,j,k)-rhov_n(i-1,j,k)) &
                         + d6(3)*(rhov_n(i+3,j,k)-rhov_n(i-2,j,k))

             Drhow(i,j,k)= d6(1)*(rhow_n(i+1,j,k)-rhow_n(i  ,j,k)) &
                         + d6(2)*(rhow_n(i+2,j,k)-rhow_n(i-1,j,k)) &
                         + d6(3)*(rhow_n(i+3,j,k)-rhow_n(i-2,j,k))

             Drhoe(i,j,k)= d6(1)*(rhoe_n(i+1,j,k)-rhoe_n(i  ,j,k)) &
                         + d6(2)*(rhoe_n(i+2,j,k)-rhoe_n(i-1,j,k)) &
                         + d6(3)*(rhoe_n(i+3,j,k)-rhoe_n(i-2,j,k))
          enddo
       enddo
    enddo

    if (is_boundary(1,2)) then
       ! 3rd-order centered dissipation
       ! ------------------------------
       i=nx-2
       do k=1,nz
          do j=1,ny
             Drho(i,j,k) = d4(1)*(rho_n(i+1,j,k)-rho_n(i  ,j,k)) &
                         + d4(2)*(rho_n(i+2,j,k)-rho_n(i-1,j,k))            
             Drhou(i,j,k)= d4(1)*(rhou_n(i+1,j,k)-rhou_n(i  ,j,k)) &
                         + d4(2)*(rhou_n(i+2,j,k)-rhou_n(i-1,j,k))
             Drhov(i,j,k)= d4(1)*(rhov_n(i+1,j,k)-rhov_n(i  ,j,k)) &
                         + d4(2)*(rhov_n(i+2,j,k)-rhov_n(i-1,j,k))
             Drhow(i,j,k)= d4(1)*(rhow_n(i+1,j,k)-rhow_n(i  ,j,k)) &
                         + d4(2)*(rhow_n(i+2,j,k)-rhow_n(i-1,j,k))
             Drhoe(i,j,k)= d4(1)*(rhoe_n(i+1,j,k)-rhoe_n(i  ,j,k)) &
                         + d4(2)*(rhoe_n(i+2,j,k)-rhoe_n(i-1,j,k))
          enddo
       enddo
       
       ! 1st-order centered dissipation
       ! ------------------------------
       i=nx-1
       do k=1,nz
          do j=1,ny
             Drho(i,j,k) = d2(1)*(rho_n(i+1,j,k)-rho_n(i,j,k))            
             Drhou(i,j,k)= d2(1)*(rhou_n(i+1,j,k)-rhou_n(i,j,k))
             Drhov(i,j,k)= d2(1)*(rhov_n(i+1,j,k)-rhov_n(i,j,k))
             Drhow(i,j,k)= d2(1)*(rhow_n(i+1,j,k)-rhow_n(i,j,k))
             Drhoe(i,j,k)= d2(1)*(rhoe_n(i+1,j,k)-rhoe_n(i,j,k))
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
                vc=(uu(i,j,k)*y_eta(i,j)-vv(i,j,k)*x_eta(i,j))*ijacob(i,j)
                vc2=(uu(i+1,j,k)*y_eta(i+1,j)-vv(i+1,j,k)*x_eta(i+1,j))*ijacob(i+1,j)
                vc1=(uu(i-1,j,k)*y_eta(i-1,j)-vv(i-1,j,k)*x_eta(i-1,j))*ijacob(i-1,j)
                ! compute spectral radius at i+1/2 (srp12) and i-1/2 (srm12)
                sr=abs(vc)+c_(i,j,k)*g_ksi(i,j)
                srp12=0.5_wp*(abs(vc2)+c_(i+1,j,k)*g_ksi(i+1,j)+sr)
                srm12=0.5_wp*(abs(vc1)+c_(i-1,j,k)*g_ksi(i-1,j)+sr)
                ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dksi
                Krho(i,j,k) = Krho(i,j,k) +(srp12*Drho(i,j,k) -srm12*Drho(i-1,j,k))
                Krhou(i,j,k)= Krhou(i,j,k)+(srp12*Drhou(i,j,k)-srm12*Drhou(i-1,j,k))
                Krhov(i,j,k)= Krhov(i,j,k)+(srp12*Drhov(i,j,k)-srm12*Drhov(i-1,j,k))
                Krhow(i,j,k)= Krhow(i,j,k)+(srp12*Drhow(i,j,k)-srm12*Drhow(i-1,j,k))
                Krhoe(i,j,k)= Krhoe(i,j,k)+(srp12*Drhoe(i,j,k)-srm12*Drhoe(i-1,j,k))
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
             enddo
          enddo
       enddo
    endif

    ! Artificial viscosity in y-direction [compute Drho at j+1/2]
    ! ===================================
    if (is_boundary(2,1)) then
       ! 1st-order centered dissipation
       ! ------------------------------
       j=1
       do k=1,nz
          do i=1,nx
             Drho(i,j,k) = d2(1)*(rho_n(i,j+1,k)-rho_n(i,j,k))
             Drhou(i,j,k)= d2(1)*(rhou_n(i,j+1,k)-rhou_n(i,j,k))
             Drhov(i,j,k)= d2(1)*(rhov_n(i,j+1,k)-rhov_n(i,j,k))
             Drhow(i,j,k)= d2(1)*(rhow_n(i,j+1,k)-rhow_n(i,j,k))
             Drhoe(i,j,k)= d2(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j,k))
          enddo
       enddo
       
       ! 3rd-order centered dissipation
       ! ------------------------------
       j=2
       do k=1,nz
          do i=1,nx
             Drho(i,j,k) = d4(1)*(rho_n(i,j+1,k)-rho_n(i,j  ,k)) &
                         + d4(2)*(rho_n(i,j+2,k)-rho_n(i,j-1,k))
             Drhou(i,j,k)= d4(1)*(rhou_n(i,j+1,k)-rhou_n(i,j  ,k)) &
                         + d4(2)*(rhou_n(i,j+2,k)-rhou_n(i,j-1,k))
             Drhov(i,j,k)= d4(1)*(rhov_n(i,j+1,k)-rhov_n(i,j  ,k)) &
                         + d4(2)*(rhov_n(i,j+2,k)-rhov_n(i,j-1,k))
             Drhow(i,j,k)= d4(1)*(rhow_n(i,j+1,k)-rhow_n(i,j  ,k)) &
                         + d4(2)*(rhow_n(i,j+2,k)-rhow_n(i,j-1,k))
             Drhoe(i,j,k)= d4(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j  ,k)) &
                         + d4(2)*(rhoe_n(i,j+2,k)-rhoe_n(i,j-1,k))
          enddo
       enddo
    endif
    
    ! Interior points: 5th-order centered dissipation
    ! -----------------------------------------------
    do k=1,nz
       do j=ndy_di,nfy_di
          do i=1,nx
             Drho(i,j,k) = d6(1)*(rho_n(i,j+1,k)-rho_n(i,j  ,k)) &
                         + d6(2)*(rho_n(i,j+2,k)-rho_n(i,j-1,k)) &
                         + d6(3)*(rho_n(i,j+3,k)-rho_n(i,j-2,k))

             Drhou(i,j,k)= d6(1)*(rhou_n(i,j+1,k)-rhou_n(i,j  ,k)) &
                         + d6(2)*(rhou_n(i,j+2,k)-rhou_n(i,j-1,k)) &
                         + d6(3)*(rhou_n(i,j+3,k)-rhou_n(i,j-2,k))

             Drhov(i,j,k)= d6(1)*(rhov_n(i,j+1,k)-rhov_n(i,j  ,k)) &
                         + d6(2)*(rhov_n(i,j+2,k)-rhov_n(i,j-1,k)) &
                         + d6(3)*(rhov_n(i,j+3,k)-rhov_n(i,j-2,k))

             Drhow(i,j,k)= d6(1)*(rhow_n(i,j+1,k)-rhow_n(i,j  ,k)) &
                         + d6(2)*(rhow_n(i,j+2,k)-rhow_n(i,j-1,k)) &
                         + d6(3)*(rhow_n(i,j+3,k)-rhow_n(i,j-2,k))
             
             Drhoe(i,j,k)= d6(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j  ,k)) &
                         + d6(2)*(rhoe_n(i,j+2,k)-rhoe_n(i,j-1,k)) &
                         + d6(3)*(rhoe_n(i,j+3,k)-rhoe_n(i,j-2,k))
          enddo
       enddo
    enddo

    if (is_boundary(2,2)) then
       ! 3rd-order centered dissipation
       ! ------------------------------
       j=ny-2
       do k=1,nz
          do i=1,nx
             Drho(i,j,k) = d4(1)*(rho_n(i,j+1,k)-rho_n(i,j  ,k)) &
                         + d4(2)*(rho_n(i,j+2,k)-rho_n(i,j-1,k))
             Drhou(i,j,k)= d4(1)*(rhou_n(i,j+1,k)-rhou_n(i,j  ,k)) &
                         + d4(2)*(rhou_n(i,j+2,k)-rhou_n(i,j-1,k))
             Drhov(i,j,k)= d4(1)*(rhov_n(i,j+1,k)-rhov_n(i,j  ,k)) &
                         + d4(2)*(rhov_n(i,j+2,k)-rhov_n(i,j-1,k))
             Drhow(i,j,k)= d4(1)*(rhow_n(i,j+1,k)-rhow_n(i,j  ,k)) &
                         + d4(2)*(rhow_n(i,j+2,k)-rhow_n(i,j-1,k))
             Drhoe(i,j,k)= d4(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j  ,k)) &
                         + d4(2)*(rhoe_n(i,j+2,k)-rhoe_n(i,j-1,k))
          enddo
       enddo
       
       ! 1st-order centered dissipation
       ! ------------------------------
       j=ny-1
       do k=1,nz
          do i=1,nx
             Drho(i,j,k) = d2(1)*(rho_n(i,j+1,k)-rho_n(i,j,k))
             Drhou(i,j,k)= d2(1)*(rhou_n(i,j+1,k)-rhou_n(i,j,k))
             Drhov(i,j,k)= d2(1)*(rhov_n(i,j+1,k)-rhov_n(i,j,k))
             Drhow(i,j,k)= d2(1)*(rhow_n(i,j+1,k)-rhow_n(i,j,k))
             Drhoe(i,j,k)= d2(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j,k))
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
                vc=(vv(i,j,k)*x_ksi(i,j)-uu(i,j,k)*y_ksi(i,j))*ijacob(i,j)
                vc2=(vv(i,j+1,k)*x_ksi(i,j+1)-uu(i,j+1,k)*y_ksi(i,j+1))*ijacob(i,j+1)
                vc1=(vv(i,j-1,k)*x_ksi(i,j-1)-uu(i,j-1,k)*y_ksi(i,j-1))*ijacob(i,j-1)
                ! compute spectral radius at j+1/2 (srp12) and j-1/2 (srm12)
                sr=abs(vc)+c_(i,j,k)*g_eta(i,j)
                srp12=0.5_wp*(abs(vc2)+c_(i,j+1,k)*g_eta(i,j+1)+sr)
                srm12=0.5_wp*(abs(vc1)+c_(i,j-1,k)*g_eta(i,j-1)+sr)
                ! compute dissipation term: delta_j/delta_y=(D_{j+1/2}-D_{j-1/2})/dy
                Krho(i,j,k) = Krho(i,j,k) +(srp12*Drho(i,j,k) -srm12*Drho(i,j-1,k) )
                Krhou(i,j,k)= Krhou(i,j,k)+(srp12*Drhou(i,j,k)-srm12*Drhou(i,j-1,k))
                Krhov(i,j,k)= Krhov(i,j,k)+(srp12*Drhov(i,j,k)-srm12*Drhov(i,j-1,k))
                Krhow(i,j,k)= Krhow(i,j,k)+(srp12*Drhow(i,j,k)-srm12*Drhow(i,j-1,k))
                Krhoe(i,j,k)= Krhoe(i,j,k)+(srp12*Drhoe(i,j,k)-srm12*Drhoe(i,j-1,k))
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
             enddo
          enddo
       enddo
    endif
 
    !*******************
    if (.not.is_2D) then
    !*******************

       ! Artificial viscosity in z-direction [compute Drho at k+1/2]
       ! ===================================
       if (is_boundary(3,1)) then
          ! 1st-order centered dissipation
          ! ------------------------------
          k=1
          do j=1,ny
             do i=1,nx
                Drho(i,j,k) = d2(1)*(rho_n(i,j,k+1)-rho_n(i,j,k))
                Drhou(i,j,k)= d2(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k))
                Drhov(i,j,k)= d2(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k))
                Drhow(i,j,k)= d2(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k))
                Drhoe(i,j,k)= d2(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k))
             enddo
          enddo
          
          ! 3rd-order centered dissipation
          ! ------------------------------
          k=2
          do j=1,ny
             do i=1,nx
                Drho(i,j,k) = d4(1)*(rho_n(i,j,k+1)-rho_n(i,j,k  )) &
                            + d4(2)*(rho_n(i,j,k+2)-rho_n(i,j,k-1))
                Drhou(i,j,k)= d4(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k  )) &
                            + d4(2)*(rhou_n(i,j,k+2)-rhou_n(i,j,k-1))
                Drhov(i,j,k)= d4(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k  )) &
                            + d4(2)*(rhov_n(i,j,k+2)-rhov_n(i,j,k-1))
                Drhow(i,j,k)= d4(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k  )) &
                            + d4(2)*(rhow_n(i,j,k+2)-rhow_n(i,j,k-1))
                Drhoe(i,j,k)= d4(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k  )) &
                            + d4(2)*(rhoe_n(i,j,k+2)-rhoe_n(i,j,k-1))
             enddo
          enddo
       endif
       
       ! Interior points: 5th-order centered dissipation
       ! -----------------------------------------------
       do k=ndz_di,nfz_di
          do j=1,ny
             do i=1,nx
                Drho(i,j,k) = d6(1)*(rho_n(i,j,k+1)-rho_n(i,j,k  )) &
                            + d6(2)*(rho_n(i,j,k+2)-rho_n(i,j,k-1)) &
                            + d6(3)*(rho_n(i,j,k+3)-rho_n(i,j,k-2))

                Drhou(i,j,k)= d6(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k  )) &
                            + d6(2)*(rhou_n(i,j,k+2)-rhou_n(i,j,k-1)) &
                            + d6(3)*(rhou_n(i,j,k+3)-rhou_n(i,j,k-2))

                Drhov(i,j,k)= d6(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k  )) &
                            + d6(2)*(rhov_n(i,j,k+2)-rhov_n(i,j,k-1)) &
                            + d6(3)*(rhov_n(i,j,k+3)-rhov_n(i,j,k-2))

                Drhow(i,j,k)= d6(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k  )) &
                            + d6(2)*(rhow_n(i,j,k+2)-rhow_n(i,j,k-1)) &
                            + d6(3)*(rhow_n(i,j,k+3)-rhow_n(i,j,k-2))

                Drhoe(i,j,k)= d6(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k  )) &
                            + d6(2)*(rhoe_n(i,j,k+2)-rhoe_n(i,j,k-1)) &
                            + d6(3)*(rhoe_n(i,j,k+3)-rhoe_n(i,j,k-2))
             enddo
          enddo
       enddo

       if (is_boundary(3,2)) then
          ! 3rd-order centered dissipation
          ! ------------------------------
          k=nz-2
          do j=1,ny
             do i=1,nx
                Drho(i,j,k) = d4(1)*(rho_n(i,j,k+1)-rho_n(i,j,k  )) &
                            + d4(2)*(rho_n(i,j,k+2)-rho_n(i,j,k-1))
                Drhou(i,j,k)= d4(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k  )) &
                            + d4(2)*(rhou_n(i,j,k+2)-rhou_n(i,j,k-1))
                Drhov(i,j,k)= d4(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k  )) &
                            + d4(2)*(rhov_n(i,j,k+2)-rhov_n(i,j,k-1))
                Drhow(i,j,k)= d4(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k  )) &
                            + d4(2)*(rhow_n(i,j,k+2)-rhow_n(i,j,k-1))
                Drhoe(i,j,k)= d4(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k  )) &
                            + d4(2)*(rhoe_n(i,j,k+2)-rhoe_n(i,j,k-1))
             enddo
          enddo
          
          ! 1st-order centered dissipation
          ! ------------------------------
          k=nz-1
          do j=1,ny
             do i=1,nx
                Drho(i,j,k) = d2(1)*(rho_n(i,j,k+1)-rho_n(i,j,k))
                Drhou(i,j,k)= d2(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k))
                Drhov(i,j,k)= d2(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k))
                Drhow(i,j,k)= d2(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k))
                Drhoe(i,j,k)= d2(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k))
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
    
  end subroutine artvisc_o5_inc
  
  !==============================================================================
  subroutine artvisc_o3_inc
  !==============================================================================
    !> Apply artificial viscosity - Conservative version DNC order 3
    !> Modified version XG 06/2022
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: vc1,vc2,vc,srm12,srp12,sr
    real(wp), dimension(0:nx,0:ny,0:nz) :: Drho,Drhou,Drhov,Drhow,Drhoe
    ! ---------------------------------------------------------------------------

    ! Artificial viscosity in x-direction [compute Drho at i+1/2]
    ! ===================================
    
    if (is_boundary(1,1)) then
       ! 1st-order centered dissipation
       ! ------------------------------
       i=1
       do k=1,nz
          do j=1,ny
             Drho(i,j,k) = d2(1)*(rho_n(i+1,j,k) -rho_n(i,j,k))            
             Drhou(i,j,k)= d2(1)*(rhou_n(i+1,j,k)-rhou_n(i,j,k))
             Drhov(i,j,k)= d2(1)*(rhov_n(i+1,j,k)-rhov_n(i,j,k))
             Drhow(i,j,k)= d2(1)*(rhow_n(i+1,j,k)-rhow_n(i,j,k))
             Drhoe(i,j,k)= d2(1)*(rhoe_n(i+1,j,k)-rhoe_n(i,j,k))
          enddo
       enddo
    endif

    ! Interior points: 3rd-order centered dissipation
    ! -----------------------------------------------
    do k=1,nz
       do j=1,ny
          do i=ndx_di,nfx_di
             Drho(i,j,k) = d4(1)*(rho_n(i+1,j,k)-rho_n(i  ,j,k)) &
                         + d4(2)*(rho_n(i+2,j,k)-rho_n(i-1,j,k))
             
             Drhou(i,j,k)= d4(1)*(rhou_n(i+1,j,k)-rhou_n(i  ,j,k)) &
                         + d4(2)*(rhou_n(i+2,j,k)-rhou_n(i-1,j,k))

             Drhov(i,j,k)= d4(1)*(rhov_n(i+1,j,k)-rhov_n(i  ,j,k)) &
                         + d4(2)*(rhov_n(i+2,j,k)-rhov_n(i-1,j,k))

             Drhow(i,j,k)= d4(1)*(rhow_n(i+1,j,k)-rhow_n(i  ,j,k)) &
                         + d4(2)*(rhow_n(i+2,j,k)-rhow_n(i-1,j,k))

             Drhoe(i,j,k)= d4(1)*(rhoe_n(i+1,j,k)-rhoe_n(i  ,j,k)) &
                         + d4(2)*(rhoe_n(i+2,j,k)-rhoe_n(i-1,j,k))
          enddo
       enddo
    enddo

    if (is_boundary(1,2)) then
       ! 1st-order centered dissipation
       ! ------------------------------
       i=nx-1
       do k=1,nz
          do j=1,ny
             Drho(i,j,k) = d2(1)*(rho_n(i+1,j,k)-rho_n(i,j,k))            
             Drhou(i,j,k)= d2(1)*(rhou_n(i+1,j,k)-rhou_n(i,j,k))
             Drhov(i,j,k)= d2(1)*(rhov_n(i+1,j,k)-rhov_n(i,j,k))
             Drhow(i,j,k)= d2(1)*(rhow_n(i+1,j,k)-rhow_n(i,j,k))
             Drhoe(i,j,k)= d2(1)*(rhoe_n(i+1,j,k)-rhoe_n(i,j,k))
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
                vc=(uu(i,j,k)*y_eta(i,j)-vv(i,j,k)*x_eta(i,j))*ijacob(i,j)
                vc2=(uu(i+1,j,k)*y_eta(i+1,j)-vv(i+1,j,k)*x_eta(i+1,j))*ijacob(i+1,j)
                vc1=(uu(i-1,j,k)*y_eta(i-1,j)-vv(i-1,j,k)*x_eta(i-1,j))*ijacob(i-1,j)
                ! compute spectral radius at i+1/2 (srp12) and i-1/2 (srm12)
                sr=abs(vc)+c_(i,j,k)*g_ksi(i,j)
                srp12=0.5_wp*(abs(vc2)+c_(i+1,j,k)*g_ksi(i+1,j)+sr)
                srm12=0.5_wp*(abs(vc1)+c_(i-1,j,k)*g_ksi(i-1,j)+sr)
                ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dksi
                Krho(i,j,k) = Krho(i,j,k) +(srp12*Drho(i,j,k) -srm12*Drho(i-1,j,k))
                Krhou(i,j,k)= Krhou(i,j,k)+(srp12*Drhou(i,j,k)-srm12*Drhou(i-1,j,k))
                Krhov(i,j,k)= Krhov(i,j,k)+(srp12*Drhov(i,j,k)-srm12*Drhov(i-1,j,k))
                Krhow(i,j,k)= Krhow(i,j,k)+(srp12*Drhow(i,j,k)-srm12*Drhow(i-1,j,k))
                Krhoe(i,j,k)= Krhoe(i,j,k)+(srp12*Drhoe(i,j,k)-srm12*Drhoe(i-1,j,k))
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
             enddo
          enddo
       enddo
    endif

    ! Artificial viscosity in y-direction [compute Drho at j+1/2]
    ! ===================================
    if (is_boundary(2,1)) then
       ! 1st-order centered dissipation
       ! ------------------------------
       j=1
       do k=1,nz
          do i=1,nx
             Drho(i,j,k) = d2(1)*(rho_n(i,j+1,k)-rho_n(i,j,k))
             Drhou(i,j,k)= d2(1)*(rhou_n(i,j+1,k)-rhou_n(i,j,k))
             Drhov(i,j,k)= d2(1)*(rhov_n(i,j+1,k)-rhov_n(i,j,k))
             Drhow(i,j,k)= d2(1)*(rhow_n(i,j+1,k)-rhow_n(i,j,k))
             Drhoe(i,j,k)= d2(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j,k))
          enddo
       enddo
    endif
    
    ! Interior points: 3rd-order centered dissipation
    ! -----------------------------------------------
    do k=1,nz
       do j=ndy_di,nfy_di
          do i=1,nx
             Drho(i,j,k) = d4(1)*(rho_n(i,j+1,k)-rho_n(i,j  ,k)) &
                         + d4(2)*(rho_n(i,j+2,k)-rho_n(i,j-1,k))

             Drhou(i,j,k)= d4(1)*(rhou_n(i,j+1,k)-rhou_n(i,j  ,k)) &
                         + d4(2)*(rhou_n(i,j+2,k)-rhou_n(i,j-1,k))

             Drhov(i,j,k)= d4(1)*(rhov_n(i,j+1,k)-rhov_n(i,j  ,k)) &
                         + d4(2)*(rhov_n(i,j+2,k)-rhov_n(i,j-1,k))

             Drhow(i,j,k)= d4(1)*(rhow_n(i,j+1,k)-rhow_n(i,j  ,k)) &
                         + d4(2)*(rhow_n(i,j+2,k)-rhow_n(i,j-1,k))
             
             Drhoe(i,j,k)= d4(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j  ,k)) &
                         + d4(2)*(rhoe_n(i,j+2,k)-rhoe_n(i,j-1,k))
          enddo
       enddo
    enddo

    if (is_boundary(2,2)) then
       ! 1st-order centered dissipation
       ! ------------------------------
       j=ny-1
       do k=1,nz
          do i=1,nx
             Drho(i,j,k) = d2(1)*(rho_n(i,j+1,k)-rho_n(i,j,k))
             Drhou(i,j,k)= d2(1)*(rhou_n(i,j+1,k)-rhou_n(i,j,k))
             Drhov(i,j,k)= d2(1)*(rhov_n(i,j+1,k)-rhov_n(i,j,k))
             Drhow(i,j,k)= d2(1)*(rhow_n(i,j+1,k)-rhow_n(i,j,k))
             Drhoe(i,j,k)= d2(1)*(rhoe_n(i,j+1,k)-rhoe_n(i,j,k))
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
                vc=(vv(i,j,k)*x_ksi(i,j)-uu(i,j,k)*y_ksi(i,j))*ijacob(i,j)
                vc2=(vv(i,j+1,k)*x_ksi(i,j+1)-uu(i,j+1,k)*y_ksi(i,j+1))*ijacob(i,j+1)
                vc1=(vv(i,j-1,k)*x_ksi(i,j-1)-uu(i,j-1,k)*y_ksi(i,j-1))*ijacob(i,j-1)
                ! compute spectral radius at j+1/2 (srp12) and j-1/2 (srm12)
                sr=abs(vc)+c_(i,j,k)*g_eta(i,j)
                srp12=0.5_wp*(abs(vc2)+c_(i,j+1,k)*g_eta(i,j+1)+sr)
                srm12=0.5_wp*(abs(vc1)+c_(i,j-1,k)*g_eta(i,j-1)+sr)
                ! compute dissipation term: delta_j/delta_y=(D_{j+1/2}-D_{j-1/2})/dy
                Krho(i,j,k) = Krho(i,j,k) +(srp12*Drho(i,j,k) -srm12*Drho(i,j-1,k) )
                Krhou(i,j,k)= Krhou(i,j,k)+(srp12*Drhou(i,j,k)-srm12*Drhou(i,j-1,k))
                Krhov(i,j,k)= Krhov(i,j,k)+(srp12*Drhov(i,j,k)-srm12*Drhov(i,j-1,k))
                Krhow(i,j,k)= Krhow(i,j,k)+(srp12*Drhow(i,j,k)-srm12*Drhow(i,j-1,k))
                Krhoe(i,j,k)= Krhoe(i,j,k)+(srp12*Drhoe(i,j,k)-srm12*Drhoe(i,j-1,k))
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
             enddo
          enddo
       enddo
    endif
 
    !*******************
    if (.not.is_2D) then
    !*******************

       ! Artificial viscosity in z-direction [compute Drho at k+1/2]
       ! ===================================
       if (is_boundary(3,1)) then
          ! 1st-order centered dissipation
          ! ------------------------------
          k=1
          do j=1,ny
             do i=1,nx
                Drho(i,j,k) = d2(1)*(rho_n(i,j,k+1)-rho_n(i,j,k))
                Drhou(i,j,k)= d2(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k))
                Drhov(i,j,k)= d2(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k))
                Drhow(i,j,k)= d2(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k))
                Drhoe(i,j,k)= d2(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k))
             enddo
          enddo
       endif
       
       ! Interior points: 3rd-order centered dissipation
       ! -----------------------------------------------
       do k=ndz_di,nfz_di
          do j=1,ny
             do i=1,nx
                Drho(i,j,k) = d4(1)*(rho_n(i,j,k+1)-rho_n(i,j,k  )) &
                            + d4(2)*(rho_n(i,j,k+2)-rho_n(i,j,k-1))

                Drhou(i,j,k)= d4(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k  )) &
                            + d4(2)*(rhou_n(i,j,k+2)-rhou_n(i,j,k-1))

                Drhov(i,j,k)= d4(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k  )) &
                            + d4(2)*(rhov_n(i,j,k+2)-rhov_n(i,j,k-1))

                Drhow(i,j,k)= d4(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k  )) &
                            + d4(2)*(rhow_n(i,j,k+2)-rhow_n(i,j,k-1))

                Drhoe(i,j,k)= d4(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k  )) &
                            + d4(2)*(rhoe_n(i,j,k+2)-rhoe_n(i,j,k-1))
             enddo
          enddo
       enddo

       if (is_boundary(3,2)) then
          ! 1st-order centered dissipation
          ! ------------------------------
          k=nz-1
          do j=1,ny
             do i=1,nx
                Drho(i,j,k) = d2(1)*(rho_n(i,j,k+1)-rho_n(i,j,k))
                Drhou(i,j,k)= d2(1)*(rhou_n(i,j,k+1)-rhou_n(i,j,k))
                Drhov(i,j,k)= d2(1)*(rhov_n(i,j,k+1)-rhov_n(i,j,k))
                Drhow(i,j,k)= d2(1)*(rhow_n(i,j,k+1)-rhow_n(i,j,k))
                Drhoe(i,j,k)= d2(1)*(rhoe_n(i,j,k+1)-rhoe_n(i,j,k))
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
    
  end subroutine artvisc_o3_inc

end module mod_artvisc_inc
