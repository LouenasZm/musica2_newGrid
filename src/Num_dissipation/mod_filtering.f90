! ===============================================================================
module mod_filtering
! ===============================================================================
  !> Module for filtering variables
! ===============================================================================
  use precision
  implicit none
  ! -----------------------------------------------------------------------------
  ! coefficients for filters
  real(wp) :: chi
  real(wp) :: d3(0:1),d5(0:2),d7(0:3),d9(0:4),d11(0:5),d13(0:6)
  real(wp) :: d15(7),d24(7),d28(11),d37(11),d46(11)
  ! -----------------------------------------------------------------------------

contains

  !==============================================================================
  subroutine init_coeff_filters(stencil,is_DRP)
  !==============================================================================
    !> Initialization of selective filter coefficients
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    ! Scheme parameters
    integer, intent(in)  :: stencil
    logical, intent(in)  :: is_DRP
    ! ---------------------------------------------------------------------------

    ! Standard filter order 2 (3pts)
    ! -----------------------
    d3(0)=  0.50_wp
    d3(1)= -0.25_wp

    ! Standard filter order 4 (5pts)
    ! -----------------------
    d5(0)=  6.0_wp/16.0_wp
    d5(1)= -4.0_wp/16.0_wp
    d5(2)=  1.0_wp/16.0_wp

    if (stencil<7) return

    ! Standard filter order 6 (7pts)
    ! -----------------------
    d7(0)=  20.0_wp/64.0_wp
    d7(1)= -15.0_wp/64.0_wp
    d7(2)=   6.0_wp/64.0_wp
    d7(3)= - 1.0_wp/64.0_wp

    if (stencil<9) return

    ! Standard filter order 8 (9pts)
    ! -----------------------
    d9(0)=  70.0_wp/256.0_wp
    d9(1)= -56.0_wp/256.0_wp
    d9(2)=  28.0_wp/256.0_wp
    d9(3)= - 8.0_wp/256.0_wp
    d9(4)=   1.0_wp/256.0_wp

    if (stencil<11) return

    ! Standard filter order 10 (11pts)
    ! ------------------------
    d11(0)= 252.0_wp/1024.0_wp
    d11(1)=-210.0_wp/1024.0_wp
    d11(2)= 120.0_wp/1024.0_wp
    d11(3)=- 45.0_wp/1024.0_wp
    d11(4)=  10.0_wp/1024.0_wp
    d11(5)=-  1.0_wp/1024.0_wp

    ! if DRP (Dispersion Relation Preserving) and 11-pt stencil,
    ! replace by optimized scheme
    if ((stencil==11).and.(is_DRP)) then
       ! DRP filter 11 points (DRP order 2)
       ! ----------------------------------
       ! [Bogey & Bailly JCP 2004]
       !d11(0)=  0.215044884112_wp
       !d11(1)= -0.187772883589_wp
       !d11(2)=  0.123755948787_wp
       !d11(3)= -0.059227575576_wp
       !d11(4)=  0.018721609157_wp
       !d11(5)= -0.002999540835_wp

       ! DRP filter 11 points (DRP order 6)
       ! ----------------------------------
       ! [Bogey, de Cacqueray, Bailly JCP 2009]
       d11(0)=  0.234810479761700_wp
       d11(1)= -0.199250131285813_wp
       d11(2)=  0.120198310245186_wp
       d11(3)= -0.049303775636020_wp
       d11(4)=  0.012396449873964_wp
       d11(5)= -0.001446093078167_wp
    endif

    ! Optimized boundary filters
    ! --------------------------
    ! [Berland, Marsden, Bogey & Bailly JCP 2007]

    ! Boundary filter - Non-centered stencil 1-5
    d15(1)= -0.085777408970_wp+ 0.000000000001_wp ! modified
    d15(2)=  0.277628171524_wp
    d15(3)= -0.356848072173_wp
    d15(4)=  0.223119093072_wp
    d15(5)= -0.057347064865_wp
    d15(6)= -0.000747264596_wp
    d15(7)= -0.000027453993_wp

    ! Boundary filter - Non-centered stencil 2-4
    d24(1)=  0.032649010764_wp
    d24(2)= -0.143339502575_wp
    d24(3)=  0.273321177980_wp
    d24(4)= -0.294622121167_wp
    d24(5)=  0.186711738069_wp
    d24(6)= -0.062038376258_wp
    d24(7)=  0.007318073189_wp

    ! Boundary filter - Non-centered stencil 2-8
    d28(1) =  0.0307159855992469_wp
    d28(2) = -0.148395705486028_wp
    d28(3) =  0.312055385963757_wp
    d28(4) = -0.363202245195514_wp
    d28(5) =  0.230145457063431_wp
    d28(6) = -0.0412316564605079_wp
    d28(7) = -0.0531024700805787_wp
    d28(8) =  0.0494343261171287_wp
    d28(9) = -0.0198143585458560_wp
    d28(10)=  0.00339528102492129_wp
    d28(11)=  0.0_wp

    ! Boundary filter - Non-centered stencil 3-7
    d37(1) = -0.000054596010_wp
    d37(2) =  0.042124772446_wp
    d37(3) = -0.173103107841_wp
    d37(4) =  0.299615871352_wp
    d37(5) = -0.276543612935_wp
    d37(6) =  0.131223506571_wp
    d37(7) = -0.023424966418_wp
    d37(8) =  0.013937561779_wp
    d37(9) = -0.024565095706_wp
    d37(10)=  0.013098287852_wp
    d37(11)= -0.002308621090_wp

    ! Boundary filter - Non-centered stencil 4-6
    d46(1) =  0.008391235145_wp
    d46(2) = -0.047402506444_wp
    d46(3) =  0.121438547725_wp
    d46(4) = -0.200063042812_wp
    d46(5) =  0.240069047836_wp
    d46(6) = -0.207269200140_wp- 0.000000000001_wp ! modified
    d46(7) =  0.122263107844_wp- 0.000000000001_wp ! modified
    d46(8) = -0.047121062819_wp
    d46(9) =  0.009014891495_wp
    d46(10)=  0.001855812216_wp
    d46(11)= -0.001176830044_wp

    if (stencil<13) return

    ! DRP filter 13 points (DRP order 4)
    ! ----------------------------------
    ! [Bogey & Bailly JCP 2004]
    d13(0)= 0.1935444903379030_wp
    d13(1)=-0.1733088204161889_wp
    d13(2)= 0.1236749560113578_wp
    d13(3)=-6.87982177922706e-2_wp
    d13(4)= 2.84391878109580e-2_wp
    d13(5)=-7.8929617915404e-3_wp
    d13(6)= 1.1136110087326e-3_wp

  end subroutine init_coeff_filters

  !==============================================================================
  subroutine init_filtering(stencil,is_DRP,X_sf)
  !==============================================================================
    !> Initialization of selective filters
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    ! Scheme parameters
    integer, intent(in)  :: stencil
    logical, intent(in)  :: is_DRP
    real(wp), intent(in) :: X_sf
    ! ---------------------------------------------------------------------------

    !!chi=1.0_wp!/X_sf
    chi=X_sf
    
    ! Fill selective filter coefficients
    ! ==================================
    call init_coeff_filters(stencil,is_DRP)

    ! Multiply coeff by filtering amplitude (chi_SF) and change sign
    ! ==============================================================
    d3 =-d3*X_sf
    d5 =-d5*X_sf
    d7 =-d7*X_sf
    d9 =-d9*X_sf
    d11=-d11*X_sf
    d13=-d13*X_sf
    
    d15=-d15*X_sf
    d24=-d24*X_sf
    d28=-d28*X_sf
    d37=-d37*X_sf
    d46=-d46*X_sf

  end subroutine init_filtering

  !==============================================================================
  subroutine filtering_11pts_
  !==============================================================================
    !> Apply filtering on 11-point stencils (+ boundaries)
    !> * old version without comput./comm. overlap *
  !==============================================================================
    use mod_flow
    use mod_time
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    ! ---------------------------------------------------------------------------
    real(wp) :: uc,vc,cfll,mul
    real(wp), dimension(nx,ny,nz) :: Drho,Drhou,Drhov,Drhow,Drhoe
    
    ! Filtering along i-direction
    ! ===========================
    
    ! BC at imin
    ! ----------
    if (is_boundary(1,1)) then
       ! no filter for i=1 & i=2 -> initialize array
       do j=1,ny
           Drho(1,j,:)=0.0_wp
          Drhou(1,j,:)=0.0_wp
          Drhov(1,j,:)=0.0_wp
          Drhow(1,j,:)=0.0_wp
          Drhoe(1,j,:)=0.0_wp
          
           Drho(2,j,:)=0.0_wp
          Drhou(2,j,:)=0.0_wp
          Drhov(2,j,:)=0.0_wp
          Drhow(2,j,:)=0.0_wp
          Drhoe(2,j,:)=0.0_wp
       enddo
       ! 4-th order centered filter
       i=3
       do k=1,nz
          do j=1,ny
             Drho(i,j,k) = d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i+1,j,k)+d5(2)*rho_n(i+2,j,k) &
                  + d5(1)*rho_n(i-1,j,k)+d5(2)*rho_n(i-2,j,k)
             Drhou(i,j,k)= d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i+1,j,k)+d5(2)*rhou_n(i+2,j,k) &
                  + d5(1)*rhou_n(i-1,j,k)+d5(2)*rhou_n(i-2,j,k)
             Drhov(i,j,k)= d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i+1,j,k)+d5(2)*rhov_n(i+2,j,k) &
                  + d5(1)*rhov_n(i-1,j,k)+d5(2)*rhov_n(i-2,j,k)
             Drhow(i,j,k)= d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i+1,j,k)+d5(2)*rhow_n(i+2,j,k) &
                  + d5(1)*rhow_n(i-1,j,k)+d5(2)*rhow_n(i-2,j,k)
             Drhoe(i,j,k)= d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i+1,j,k)+d5(2)*rhoe_n(i+2,j,k) &
                  + d5(1)*rhoe_n(i-1,j,k)+d5(2)*rhoe_n(i-2,j,k)
          enddo
       enddo
       ! 6-th order centered filter
       i=4
       do k=1,nz
          do j=1,ny
             Drho(i,j,k) = d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i+1,j,k)+d7(2)*rho_n(i+2,j,k)+d7(3)*rho_n(i+3,j,k) &
                  + d7(1)*rho_n(i-1,j,k)+d7(2)*rho_n(i-2,j,k)+d7(3)*rho_n(i-3,j,k)
             Drhou(i,j,k)= d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i+1,j,k)+d7(2)*rhou_n(i+2,j,k)+d7(3)*rhou_n(i+3,j,k) &
                  + d7(1)*rhou_n(i-1,j,k)+d7(2)*rhou_n(i-2,j,k)+d7(3)*rhou_n(i-3,j,k)
             Drhov(i,j,k)= d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i+1,j,k)+d7(2)*rhov_n(i+2,j,k)+d7(3)*rhov_n(i+3,j,k) &
                  + d7(1)*rhov_n(i-1,j,k)+d7(2)*rhov_n(i-2,j,k)+d7(3)*rhov_n(i-3,j,k)
             Drhow(i,j,k)= d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i+1,j,k)+d7(2)*rhow_n(i+2,j,k)+d7(3)*rhow_n(i+3,j,k) &
                  + d7(1)*rhow_n(i-1,j,k)+d7(2)*rhow_n(i-2,j,k)+d7(3)*rhow_n(i-3,j,k)
             Drhoe(i,j,k)= d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i+1,j,k)+d7(2)*rhoe_n(i+2,j,k)+d7(3)*rhoe_n(i+3,j,k) &
                  + d7(1)*rhoe_n(i-1,j,k)+d7(2)*rhoe_n(i-2,j,k)+d7(3)*rhoe_n(i-3,j,k)
          enddo
       enddo
       ! 8-th order centered filter
       i=5
       do k=1,nz
          do j=1,ny
             Drho(i,j,k) = d9(0)*rho_n(i,j,k)   &
                  + d9(1)*rho_n(i+1,j,k)+d9(2)*rho_n(i+2,j,k)+d9(3)*rho_n(i+3,j,k)+d9(4)*rho_n(i+4,j,k) &
                  + d9(1)*rho_n(i-1,j,k)+d9(2)*rho_n(i-2,j,k)+d9(3)*rho_n(i-3,j,k)+d9(4)*rho_n(i-4,j,k)
             Drhou(i,j,k)= d9(0)*rhou_n(i,j,k)   &
                  + d9(1)*rhou_n(i+1,j,k)+d9(2)*rhou_n(i+2,j,k)+d9(3)*rhou_n(i+3,j,k)+d9(4)*rhou_n(i+4,j,k) &
                  + d9(1)*rhou_n(i-1,j,k)+d9(2)*rhou_n(i-2,j,k)+d9(3)*rhou_n(i-3,j,k)+d9(4)*rhou_n(i-4,j,k)
             Drhov(i,j,k)= d9(0)*rhov_n(i,j,k)   &
                  + d9(1)*rhov_n(i+1,j,k)+d9(2)*rhov_n(i+2,j,k)+d9(3)*rhov_n(i+3,j,k)+d9(4)*rhov_n(i+4,j,k) &
                  + d9(1)*rhov_n(i-1,j,k)+d9(2)*rhov_n(i-2,j,k)+d9(3)*rhov_n(i-3,j,k)+d9(4)*rhov_n(i-4,j,k)
             Drhow(i,j,k)= d9(0)*rhow_n(i,j,k)   &
                  + d9(1)*rhow_n(i+1,j,k)+d9(2)*rhow_n(i+2,j,k)+d9(3)*rhow_n(i+3,j,k)+d9(4)*rhow_n(i+4,j,k) &
                  + d9(1)*rhow_n(i-1,j,k)+d9(2)*rhow_n(i-2,j,k)+d9(3)*rhow_n(i-3,j,k)+d9(4)*rhow_n(i-4,j,k)
             Drhoe(i,j,k)= d9(0)*rhoe_n(i,j,k)   &
                  + d9(1)*rhoe_n(i+1,j,k)+d9(2)*rhoe_n(i+2,j,k)+d9(3)*rhoe_n(i+3,j,k)+d9(4)*rhoe_n(i+4,j,k) &
                  + d9(1)*rhoe_n(i-1,j,k)+d9(2)*rhoe_n(i-2,j,k)+d9(3)*rhoe_n(i-3,j,k)+d9(4)*rhoe_n(i-4,j,k)
          enddo
       enddo
    endif

    ! Interior points
    ! ---------------
    ! 10-th order centered filter
    do k=1,nz
       do j=1,ny
          do i=ndx,nfx
             Drho(i,j,k) = d11(0)*rho_n(i,j,k) &
                  + d11(1)*( rho_n(i+1,j,k)+rho_n(i-1,j,k) )&
                  + d11(2)*( rho_n(i+2,j,k)+rho_n(i-2,j,k) )&
                  + d11(3)*( rho_n(i+3,j,k)+rho_n(i-3,j,k) )&
                  + d11(4)*( rho_n(i+4,j,k)+rho_n(i-4,j,k) )&
                  + d11(5)*( rho_n(i+5,j,k)+rho_n(i-5,j,k) )

             Drhou(i,j,k)= d11(0)*rhou_n(i,j,k) &
                  + d11(1)*( rhou_n(i+1,j,k)+rhou_n(i-1,j,k) )&
                  + d11(2)*( rhou_n(i+2,j,k)+rhou_n(i-2,j,k) )&
                  + d11(3)*( rhou_n(i+3,j,k)+rhou_n(i-3,j,k) )&
                  + d11(4)*( rhou_n(i+4,j,k)+rhou_n(i-4,j,k) )&
                  + d11(5)*( rhou_n(i+5,j,k)+rhou_n(i-5,j,k) )

             Drhov(i,j,k)= d11(0)*rhov_n(i,j,k) &
                  + d11(1)*( rhov_n(i+1,j,k)+rhov_n(i-1,j,k) )&
                  + d11(2)*( rhov_n(i+2,j,k)+rhov_n(i-2,j,k) )&
                  + d11(3)*( rhov_n(i+3,j,k)+rhov_n(i-3,j,k) )&
                  + d11(4)*( rhov_n(i+4,j,k)+rhov_n(i-4,j,k) )&
                  + d11(5)*( rhov_n(i+5,j,k)+rhov_n(i-5,j,k) )

             Drhow(i,j,k)= d11(0)*rhow_n(i,j,k) &
                  + d11(1)*( rhow_n(i+1,j,k)+rhow_n(i-1,j,k) )&
                  + d11(2)*( rhow_n(i+2,j,k)+rhow_n(i-2,j,k) )&
                  + d11(3)*( rhow_n(i+3,j,k)+rhow_n(i-3,j,k) )&
                  + d11(4)*( rhow_n(i+4,j,k)+rhow_n(i-4,j,k) )&
                  + d11(5)*( rhow_n(i+5,j,k)+rhow_n(i-5,j,k) )

             Drhoe(i,j,k)= d11(0)*rhoe_n(i,j,k) &
                  + d11(1)*( rhoe_n(i+1,j,k)+rhoe_n(i-1,j,k) )&
                  + d11(2)*( rhoe_n(i+2,j,k)+rhoe_n(i-2,j,k) )&
                  + d11(3)*( rhoe_n(i+3,j,k)+rhoe_n(i-3,j,k) )&
                  + d11(4)*( rhoe_n(i+4,j,k)+rhoe_n(i-4,j,k) )&
                  + d11(5)*( rhoe_n(i+5,j,k)+rhoe_n(i-5,j,k) )
          enddo
       enddo
    enddo

    ! BC at imax
    ! ----------
    if (is_boundary(1,2)) then
       ! 8-th order centered filter
       i=nx-4
       do k=1,nz
          do j=1,ny
             Drho(i,j,k) = d9(0)*rho_n(i,j,k)   &
                  + d9(1)*rho_n(i+1,j,k)+d9(2)*rho_n(i+2,j,k)+d9(3)*rho_n(i+3,j,k)+d9(4)*rho_n(i+4,j,k) &
                  + d9(1)*rho_n(i-1,j,k)+d9(2)*rho_n(i-2,j,k)+d9(3)*rho_n(i-3,j,k)+d9(4)*rho_n(i-4,j,k)
             Drhou(i,j,k)= d9(0)*rhou_n(i,j,k)   &
                  + d9(1)*rhou_n(i+1,j,k)+d9(2)*rhou_n(i+2,j,k)+d9(3)*rhou_n(i+3,j,k)+d9(4)*rhou_n(i+4,j,k) &
                  + d9(1)*rhou_n(i-1,j,k)+d9(2)*rhou_n(i-2,j,k)+d9(3)*rhou_n(i-3,j,k)+d9(4)*rhou_n(i-4,j,k)
             Drhov(i,j,k)= d9(0)*rhov_n(i,j,k)   &
                  + d9(1)*rhov_n(i+1,j,k)+d9(2)*rhov_n(i+2,j,k)+d9(3)*rhov_n(i+3,j,k)+d9(4)*rhov_n(i+4,j,k) &
                  + d9(1)*rhov_n(i-1,j,k)+d9(2)*rhov_n(i-2,j,k)+d9(3)*rhov_n(i-3,j,k)+d9(4)*rhov_n(i-4,j,k)
             Drhow(i,j,k)= d9(0)*rhow_n(i,j,k)   &
                  + d9(1)*rhow_n(i+1,j,k)+d9(2)*rhow_n(i+2,j,k)+d9(3)*rhow_n(i+3,j,k)+d9(4)*rhow_n(i+4,j,k) &
                  + d9(1)*rhow_n(i-1,j,k)+d9(2)*rhow_n(i-2,j,k)+d9(3)*rhow_n(i-3,j,k)+d9(4)*rhow_n(i-4,j,k)
             Drhoe(i,j,k)= d9(0)*rhoe_n(i,j,k)   &
                  + d9(1)*rhoe_n(i+1,j,k)+d9(2)*rhoe_n(i+2,j,k)+d9(3)*rhoe_n(i+3,j,k)+d9(4)*rhoe_n(i+4,j,k) &
                  + d9(1)*rhoe_n(i-1,j,k)+d9(2)*rhoe_n(i-2,j,k)+d9(3)*rhoe_n(i-3,j,k)+d9(4)*rhoe_n(i-4,j,k)
          enddo
       enddo
       ! 6-th order centered filter
       i=nx-3
       do k=1,nz
          do j=1,ny
             Drho(i,j,k) = d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i+1,j,k)+d7(2)*rho_n(i+2,j,k)+d7(3)*rho_n(i+3,j,k) &
                  + d7(1)*rho_n(i-1,j,k)+d7(2)*rho_n(i-2,j,k)+d7(3)*rho_n(i-3,j,k)
             Drhou(i,j,k)= d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i+1,j,k)+d7(2)*rhou_n(i+2,j,k)+d7(3)*rhou_n(i+3,j,k) &
                  + d7(1)*rhou_n(i-1,j,k)+d7(2)*rhou_n(i-2,j,k)+d7(3)*rhou_n(i-3,j,k)
             Drhov(i,j,k)= d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i+1,j,k)+d7(2)*rhov_n(i+2,j,k)+d7(3)*rhov_n(i+3,j,k) &
                  + d7(1)*rhov_n(i-1,j,k)+d7(2)*rhov_n(i-2,j,k)+d7(3)*rhov_n(i-3,j,k)
             Drhow(i,j,k)= d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i+1,j,k)+d7(2)*rhow_n(i+2,j,k)+d7(3)*rhow_n(i+3,j,k) &
                  + d7(1)*rhow_n(i-1,j,k)+d7(2)*rhow_n(i-2,j,k)+d7(3)*rhow_n(i-3,j,k)
             Drhoe(i,j,k)= d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i+1,j,k)+d7(2)*rhoe_n(i+2,j,k)+d7(3)*rhoe_n(i+3,j,k) &
                  + d7(1)*rhoe_n(i-1,j,k)+d7(2)*rhoe_n(i-2,j,k)+d7(3)*rhoe_n(i-3,j,k)
          enddo
       enddo
       ! 4-th order centered filter
       i=nx-2
       do k=1,nz
          do j=1,ny
             Drho(i,j,k) = d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i+1,j,k)+d5(2)*rho_n(i+2,j,k) &
                  + d5(1)*rho_n(i-1,j,k)+d5(2)*rho_n(i-2,j,k)
             Drhou(i,j,k)= d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i+1,j,k)+d5(2)*rhou_n(i+2,j,k) &
                  + d5(1)*rhou_n(i-1,j,k)+d5(2)*rhou_n(i-2,j,k)
             Drhov(i,j,k)= d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i+1,j,k)+d5(2)*rhov_n(i+2,j,k) &
                  + d5(1)*rhov_n(i-1,j,k)+d5(2)*rhov_n(i-2,j,k)
             Drhow(i,j,k)= d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i+1,j,k)+d5(2)*rhow_n(i+2,j,k) &
                  + d5(1)*rhow_n(i-1,j,k)+d5(2)*rhow_n(i-2,j,k)
             Drhoe(i,j,k)= d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i+1,j,k)+d5(2)*rhoe_n(i+2,j,k) &
                  + d5(1)*rhoe_n(i-1,j,k)+d5(2)*rhoe_n(i-2,j,k)
          enddo
       enddo
       ! no filter for i=nx-1 & i=nx -> initialize array
       do j=1,ny
           Drho(nx-1,j,:)=0.0_wp
          Drhou(nx-1,j,:)=0.0_wp
          Drhov(nx-1,j,:)=0.0_wp
          Drhow(nx-1,j,:)=0.0_wp
          Drhoe(nx-1,j,:)=0.0_wp
          
           Drho(nx,j,:)=0.0_wp
          Drhou(nx,j,:)=0.0_wp
          Drhov(nx,j,:)=0.0_wp
          Drhow(nx,j,:)=0.0_wp
          Drhoe(nx,j,:)=0.0_wp
       enddo
    endif
    
    ! Compute artificial dissipation terms
    ! ------------------------------------
    if (is_curv) then
       do k=1,nz
          do j=1,ny
             do i=1,nx
!!$                ! contravariant velocity
!!$                uc=(uu(i,j,k)*y_eta(i,j)-vv(i,j,k)*x_eta(i,j))*ijacob(i,j)
!!$                ! compute local CFL
!!$                cfll=(abs(uc)+c_(i,j,k)*g_ksi(i,j))*deltat
!!$                ! switch
!!$                mul=min(1.0_wp,cfll/chi)                            
                mul=min(1.0_wp,cfl_i(i,j,k)/chi)                            
                ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dx
                Krho(i,j,k) = mul*Drho(i,j,k)
                Krhou(i,j,k)= mul*Drhou(i,j,k)
                Krhov(i,j,k)= mul*Drhov(i,j,k)
                Krhow(i,j,k)= mul*Drhow(i,j,k)
                Krhoe(i,j,k)= mul*Drhoe(i,j,k)
             enddo
          enddo
       enddo
    else
       do k=1,nz
          do j=1,ny
             do i=ndx_d,nfx_d
!!$                ! compute local CFL
!!$                cfll=(abs(uu(i,j,k))+c_(i,j,k))*deltat*idx(i)
!!$                ! switch
!!$                mul=min(1.0_wp,cfll/chi)                            
                mul=min(1.0_wp,cfl_i(i,j,k)/chi)                            
                ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dx
                Krho(i,j,k) = mul*Drho(i,j,k)
                Krhou(i,j,k)= mul*Drhou(i,j,k)
                Krhov(i,j,k)= mul*Drhov(i,j,k)
                Krhow(i,j,k)= mul*Drhow(i,j,k)
                Krhoe(i,j,k)= mul*Drhoe(i,j,k)
             enddo
          enddo
       enddo
    endif
    
    ! Filtering along j-direction
    ! ===========================
    
    ! BC at jmin
    ! ----------
    if (is_boundary(2,1)) then
       ! no filter for j=1 & j=2
       
       ! 4-th order centered filter
       j=3
       do k=1,nz
          do i=1,nx
             Drho(i,j,k) = Drho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i,j+1,k)+d5(2)*rho_n(i,j+2,k) &
                  + d5(1)*rho_n(i,j-1,k)+d5(2)*rho_n(i,j-2,k)
             Drhou(i,j,k)= Drhou(i,j,k)  +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j+1,k)+d5(2)*rhou_n(i,j+2,k) &
                  + d5(1)*rhou_n(i,j-1,k)+d5(2)*rhou_n(i,j-2,k)
             Drhov(i,j,k)= Drhov(i,j,k)  +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j+1,k)+d5(2)*rhov_n(i,j+2,k) &
                  + d5(1)*rhov_n(i,j-1,k)+d5(2)*rhov_n(i,j-2,k)
             Drhow(i,j,k)= Drhow(i,j,k)  +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j+1,k)+d5(2)*rhow_n(i,j+2,k) &
                  + d5(1)*rhow_n(i,j-1,k)+d5(2)*rhow_n(i,j-2,k)
             Drhoe(i,j,k)= Drhoe(i,j,k)  +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i,j+1,k)+d5(2)*rhoe_n(i,j+2,k) &
                  + d5(1)*rhoe_n(i,j-1,k)+d5(2)*rhoe_n(i,j-2,k)
          enddo
       enddo
       ! 6-th order centered filter
       j=4
       do k=1,nz
          do i=1,nx
             Drho(i,j,k) = Drho(i,j,k)  +d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i,j+1,k)+d7(2)*rho_n(i,j+2,k)+d7(3)*rho_n(i,j+3,k) &
                  + d7(1)*rho_n(i,j-1,k)+d7(2)*rho_n(i,j-2,k)+d7(3)*rho_n(i,j-3,k)
             Drhou(i,j,k)= Drhou(i,j,k)  +d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i,j+1,k)+d7(2)*rhou_n(i,j+2,k)+d7(3)*rhou_n(i,j+3,k) &
                  + d7(1)*rhou_n(i,j-1,k)+d7(2)*rhou_n(i,j-2,k)+d7(3)*rhou_n(i,j-3,k)
             Drhov(i,j,k)= Drhov(i,j,k)  +d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i,j+1,k)+d7(2)*rhov_n(i,j+2,k)+d7(3)*rhov_n(i,j+3,k) &
                  + d7(1)*rhov_n(i,j-1,k)+d7(2)*rhov_n(i,j-2,k)+d7(3)*rhov_n(i,j-3,k)
             Drhow(i,j,k)= Drhow(i,j,k)  +d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i,j+1,k)+d7(2)*rhow_n(i,j+2,k)+d7(3)*rhow_n(i,j+3,k) &
                  + d7(1)*rhow_n(i,j-1,k)+d7(2)*rhow_n(i,j-2,k)+d7(3)*rhow_n(i,j-3,k)
             Drhoe(i,j,k)= Drhoe(i,j,k)  +d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i,j+1,k)+d7(2)*rhoe_n(i,j+2,k)+d7(3)*rhoe_n(i,j+3,k) &
                  + d7(1)*rhoe_n(i,j-1,k)+d7(2)*rhoe_n(i,j-2,k)+d7(3)*rhoe_n(i,j-3,k)
          enddo
       enddo
       ! 8-th order centered filter
       j=5
       do k=1,nz
          do i=1,nx
             Drho(i,j,k) = Drho(i,j,k)  +d9(0)*rho_n(i,j,k)   &
                  + d9(1)*rho_n(i,j+1,k)+d9(2)*rho_n(i,j+2,k)+d9(3)*rho_n(i,j+3,k)+d9(4)*rho_n(i,j+4,k) &
                  + d9(1)*rho_n(i,j-1,k)+d9(2)*rho_n(i,j-2,k)+d9(3)*rho_n(i,j-3,k)+d9(4)*rho_n(i,j-4,k)
             Drhou(i,j,k)= Drhou(i,j,k)  +d9(0)*rhou_n(i,j,k)   &
                  + d9(1)*rhou_n(i,j+1,k)+d9(2)*rhou_n(i,j+2,k)+d9(3)*rhou_n(i,j+3,k)+d9(4)*rhou_n(i,j+4,k) &
                  + d9(1)*rhou_n(i,j-1,k)+d9(2)*rhou_n(i,j-2,k)+d9(3)*rhou_n(i,j-3,k)+d9(4)*rhou_n(i,j-4,k)
             Drhov(i,j,k)= Drhov(i,j,k)  +d9(0)*rhov_n(i,j,k)   &
                  + d9(1)*rhov_n(i,j+1,k)+d9(2)*rhov_n(i,j+2,k)+d9(3)*rhov_n(i,j+3,k)+d9(4)*rhov_n(i,j+4,k) &
                  + d9(1)*rhov_n(i,j-1,k)+d9(2)*rhov_n(i,j-2,k)+d9(3)*rhov_n(i,j-3,k)+d9(4)*rhov_n(i,j-4,k)
             Drhow(i,j,k)= Drhow(i,j,k)  +d9(0)*rhow_n(i,j,k)   &
                  + d9(1)*rhow_n(i,j+1,k)+d9(2)*rhow_n(i,j+2,k)+d9(3)*rhow_n(i,j+3,k)+d9(4)*rhow_n(i,j+4,k) &
                  + d9(1)*rhow_n(i,j-1,k)+d9(2)*rhow_n(i,j-2,k)+d9(3)*rhow_n(i,j-3,k)+d9(4)*rhow_n(i,j-4,k)
             Drhoe(i,j,k)= Drhoe(i,j,k)  +d9(0)*rhoe_n(i,j,k)   &
                  + d9(1)*rhoe_n(i,j+1,k)+d9(2)*rhoe_n(i,j+2,k)+d9(3)*rhoe_n(i,j+3,k)+d9(4)*rhoe_n(i,j+4,k) &
                  + d9(1)*rhoe_n(i,j-1,k)+d9(2)*rhoe_n(i,j-2,k)+d9(3)*rhoe_n(i,j-3,k)+d9(4)*rhoe_n(i,j-4,k)
          enddo
       enddo
    endif

    ! Interior points
    ! ---------------
    ! 10-th order centered filter
    do k=1,nz
       do j=ndy,nfy
          do i=1,nx
             Drho(i,j,k) = Drho(i,j,k)  + d11(0)*rho_n(i,j,k) &
                  + d11(1)*( rho_n(i,j+1,k)+rho_n(i,j-1,k) )&
                  + d11(2)*( rho_n(i,j+2,k)+rho_n(i,j-2,k) )&
                  + d11(3)*( rho_n(i,j+3,k)+rho_n(i,j-3,k) )&
                  + d11(4)*( rho_n(i,j+4,k)+rho_n(i,j-4,k) )&
                  + d11(5)*( rho_n(i,j+5,k)+rho_n(i,j-5,k) )

             Drhou(i,j,k)= Drhou(i,j,k)  + d11(0)*rhou_n(i,j,k) &
                  + d11(1)*( rhou_n(i,j+1,k)+rhou_n(i,j-1,k) )&
                  + d11(2)*( rhou_n(i,j+2,k)+rhou_n(i,j-2,k) )&
                  + d11(3)*( rhou_n(i,j+3,k)+rhou_n(i,j-3,k) )&
                  + d11(4)*( rhou_n(i,j+4,k)+rhou_n(i,j-4,k) )&
                  + d11(5)*( rhou_n(i,j+5,k)+rhou_n(i,j-5,k) )

             Drhov(i,j,k)= Drhov(i,j,k)  + d11(0)*rhov_n(i,j,k) &
                  + d11(1)*( rhov_n(i,j+1,k)+rhov_n(i,j-1,k) )&
                  + d11(2)*( rhov_n(i,j+2,k)+rhov_n(i,j-2,k) )&
                  + d11(3)*( rhov_n(i,j+3,k)+rhov_n(i,j-3,k) )&
                  + d11(4)*( rhov_n(i,j+4,k)+rhov_n(i,j-4,k) )&
                  + d11(5)*( rhov_n(i,j+5,k)+rhov_n(i,j-5,k) )

             Drhow(i,j,k)= Drhow(i,j,k)  + d11(0)*rhow_n(i,j,k) &
                  + d11(1)*( rhow_n(i,j+1,k)+rhow_n(i,j-1,k) )&
                  + d11(2)*( rhow_n(i,j+2,k)+rhow_n(i,j-2,k) )&
                  + d11(3)*( rhow_n(i,j+3,k)+rhow_n(i,j-3,k) )&
                  + d11(4)*( rhow_n(i,j+4,k)+rhow_n(i,j-4,k) )&
                  + d11(5)*( rhow_n(i,j+5,k)+rhow_n(i,j-5,k) )

             Drhoe(i,j,k)= Drhoe(i,j,k)  + d11(0)*rhoe_n(i,j,k) &
                  + d11(1)*( rhoe_n(i,j+1,k)+rhoe_n(i,j-1,k) )&
                  + d11(2)*( rhoe_n(i,j+2,k)+rhoe_n(i,j-2,k) )&
                  + d11(3)*( rhoe_n(i,j+3,k)+rhoe_n(i,j-3,k) )&
                  + d11(4)*( rhoe_n(i,j+4,k)+rhoe_n(i,j-4,k) )&
                  + d11(5)*( rhoe_n(i,j+5,k)+rhoe_n(i,j-5,k) )
          enddo
       enddo
    enddo

    ! BC at jmax
    ! ----------
    if (is_boundary(2,2)) then
       ! 8-th order centered filter
       j=ny-4
       do k=1,nz
          do i=1,nx
             Drho(i,j,k) = Drho(i,j,k)  +d9(0)*rho_n(i,j,k)   &
                  + d9(1)*rho_n(i,j+1,k)+d9(2)*rho_n(i,j+2,k)+d9(3)*rho_n(i,j+3,k)+d9(4)*rho_n(i,j+4,k) &
                  + d9(1)*rho_n(i,j-1,k)+d9(2)*rho_n(i,j-2,k)+d9(3)*rho_n(i,j-3,k)+d9(4)*rho_n(i,j-4,k)
             Drhou(i,j,k)= Drhou(i,j,k)  +d9(0)*rhou_n(i,j,k)   &
                  + d9(1)*rhou_n(i,j+1,k)+d9(2)*rhou_n(i,j+2,k)+d9(3)*rhou_n(i,j+3,k)+d9(4)*rhou_n(i,j+4,k) &
                  + d9(1)*rhou_n(i,j-1,k)+d9(2)*rhou_n(i,j-2,k)+d9(3)*rhou_n(i,j-3,k)+d9(4)*rhou_n(i,j-4,k)
             Drhov(i,j,k)= Drhov(i,j,k)  +d9(0)*rhov_n(i,j,k)   &
                  + d9(1)*rhov_n(i,j+1,k)+d9(2)*rhov_n(i,j+2,k)+d9(3)*rhov_n(i,j+3,k)+d9(4)*rhov_n(i,j+4,k) &
                  + d9(1)*rhov_n(i,j-1,k)+d9(2)*rhov_n(i,j-2,k)+d9(3)*rhov_n(i,j-3,k)+d9(4)*rhov_n(i,j-4,k)
             Drhow(i,j,k)= Drhow(i,j,k)  +d9(0)*rhow_n(i,j,k)   &
                  + d9(1)*rhow_n(i,j+1,k)+d9(2)*rhow_n(i,j+2,k)+d9(3)*rhow_n(i,j+3,k)+d9(4)*rhow_n(i,j+4,k) &
                  + d9(1)*rhow_n(i,j-1,k)+d9(2)*rhow_n(i,j-2,k)+d9(3)*rhow_n(i,j-3,k)+d9(4)*rhow_n(i,j-4,k)
             Drhoe(i,j,k)= Drhoe(i,j,k)  +d9(0)*rhoe_n(i,j,k)   &
                  + d9(1)*rhoe_n(i,j+1,k)+d9(2)*rhoe_n(i,j+2,k)+d9(3)*rhoe_n(i,j+3,k)+d9(4)*rhoe_n(i,j+4,k) &
                  + d9(1)*rhoe_n(i,j-1,k)+d9(2)*rhoe_n(i,j-2,k)+d9(3)*rhoe_n(i,j-3,k)+d9(4)*rhoe_n(i,j-4,k)
          enddo
       enddo
       ! 6-th order centered filter
       j=ny-3
       do k=1,nz
          do i=1,nx
             Drho(i,j,k) = Drho(i,j,k)  +d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i,j+1,k)+d7(2)*rho_n(i,j+2,k)+d7(3)*rho_n(i,j+3,k) &
                  + d7(1)*rho_n(i,j-1,k)+d7(2)*rho_n(i,j-2,k)+d7(3)*rho_n(i,j-3,k)
             Drhou(i,j,k)= Drhou(i,j,k)  +d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i,j+1,k)+d7(2)*rhou_n(i,j+2,k)+d7(3)*rhou_n(i,j+3,k) &
                  + d7(1)*rhou_n(i,j-1,k)+d7(2)*rhou_n(i,j-2,k)+d7(3)*rhou_n(i,j-3,k)
             Drhov(i,j,k)= Drhov(i,j,k)  +d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i,j+1,k)+d7(2)*rhov_n(i,j+2,k)+d7(3)*rhov_n(i,j+3,k) &
                  + d7(1)*rhov_n(i,j-1,k)+d7(2)*rhov_n(i,j-2,k)+d7(3)*rhov_n(i,j-3,k)
             Drhow(i,j,k)= Drhow(i,j,k)  +d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i,j+1,k)+d7(2)*rhow_n(i,j+2,k)+d7(3)*rhow_n(i,j+3,k) &
                  + d7(1)*rhow_n(i,j-1,k)+d7(2)*rhow_n(i,j-2,k)+d7(3)*rhow_n(i,j-3,k)
             Drhoe(i,j,k)= Drhoe(i,j,k)  +d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i,j+1,k)+d7(2)*rhoe_n(i,j+2,k)+d7(3)*rhoe_n(i,j+3,k) &
                  + d7(1)*rhoe_n(i,j-1,k)+d7(2)*rhoe_n(i,j-2,k)+d7(3)*rhoe_n(i,j-3,k)
          enddo
       enddo
       ! 4-th order centered filter
       j=ny-2
       do k=1,nz
          do i=1,nx
             Drho(i,j,k)= Drho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i,j+1,k)+d5(2)*rho_n(i,j+2,k) &
                  + d5(1)*rho_n(i,j-1,k)+d5(2)*rho_n(i,j-2,k)
             Drhou(i,j,k)= Drhou(i,j,k)  +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j+1,k)+d5(2)*rhou_n(i,j+2,k) &
                  + d5(1)*rhou_n(i,j-1,k)+d5(2)*rhou_n(i,j-2,k)
             Drhov(i,j,k)= Drhov(i,j,k)  +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j+1,k)+d5(2)*rhov_n(i,j+2,k) &
                  + d5(1)*rhov_n(i,j-1,k)+d5(2)*rhov_n(i,j-2,k)
             Drhow(i,j,k)= Drhow(i,j,k)  +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j+1,k)+d5(2)*rhow_n(i,j+2,k) &
                  + d5(1)*rhow_n(i,j-1,k)+d5(2)*rhow_n(i,j-2,k)
             Drhoe(i,j,k)= Drhoe(i,j,k)  +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i,j+1,k)+d5(2)*rhoe_n(i,j+2,k) &
                  + d5(1)*rhoe_n(i,j-1,k)+d5(2)*rhoe_n(i,j-2,k)
          enddo
       enddo
       ! no filter for j=ny-1 & j=ny
    endif

    if (is_curv) then
       do k=1,nz
          do j=1,ny
             do i=1,nx
!!$                ! contravariant velocity
!!$                vc=(vv(i,j,k)*x_ksi(i,j)-uu(i,j,k)*y_ksi(i,j))*ijacob(i,j)
!!$                ! compute local CFL
!!$                cfll=(abs(vc)+c_(i,j,k)*g_eta(i,j))*deltat
!!$                ! switch
!!$                mul=min(1.0_wp,cfll/chi)                            
                mul=min(1.0_wp,cfl_j(i,j,k)/chi)                            
                ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dx
                Krho(i,j,k) = Krho(i,j,k) +mul*Drho(i,j,k)
                Krhou(i,j,k)= Krhou(i,j,k)+mul*Drhou(i,j,k)
                Krhov(i,j,k)= Krhov(i,j,k)+mul*Drhov(i,j,k)
                Krhow(i,j,k)= Krhow(i,j,k)+mul*Drhow(i,j,k)
                Krhoe(i,j,k)= Krhoe(i,j,k)+mul*Drhoe(i,j,k)
             enddo
          enddo
       enddo
    else
       do k=1,nz
          do j=1,ny
             do i=1,nx
!!$                ! compute local CFL
!!$                cfll=(abs(vv(i,j,k))+c_(i,j,k))*deltat*idy(j)
!!$                ! switch
!!$                mul=min(1.0_wp,cfll/chi)                            
                mul=min(1.0_wp,cfl_j(i,j,k)/chi)                            
                ! compute dissipation term: delta_i/delta_x=(D_{i+1/2}-D_{i-1/2})/dx
                Krho(i,j,k) = Krho(i,j,k) +mul*Drho(i,j,k)
                Krhou(i,j,k)= Krhou(i,j,k)+mul*Drhou(i,j,k)
                Krhov(i,j,k)= Krhov(i,j,k)+mul*Drhov(i,j,k)
                Krhow(i,j,k)= Krhow(i,j,k)+mul*Drhow(i,j,k)
                Krhoe(i,j,k)= Krhoe(i,j,k)+mul*Drhoe(i,j,k)
             enddo
          enddo
       enddo
    endif
 
    !****************
    if (is_2D) return
    !****************
    
    ! Filtering along k-direction
    ! ===========================
    
    ! BC at kmin
    ! ----------
    if (is_boundary(3,1)) then
       ! no filter for k=1 & k=2

       ! 4-th order centered filter
       k=3
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i,j,k+1)+d5(2)*rho_n(i,j,k+2) &
                  + d5(1)*rho_n(i,j,k-1)+d5(2)*rho_n(i,j,k-2)
             Krhou(i,j,k)= Krhou(i,j,k)  +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j,k+1)+d5(2)*rhou_n(i,j,k+2) &
                  + d5(1)*rhou_n(i,j,k-1)+d5(2)*rhou_n(i,j,k-2)
             Krhov(i,j,k)= Krhov(i,j,k)  +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j,k+1)+d5(2)*rhov_n(i,j,k+2) &
                  + d5(1)*rhov_n(i,j,k-1)+d5(2)*rhov_n(i,j,k-2)
             Krhow(i,j,k)= Krhow(i,j,k)  +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j,k+1)+d5(2)*rhow_n(i,j,k+2) &
                  + d5(1)*rhow_n(i,j,k-1)+d5(2)*rhow_n(i,j,k-2)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i,j,k+1)+d5(2)*rhoe_n(i,j,k+2) &
                  + d5(1)*rhoe_n(i,j,k-1)+d5(2)*rhoe_n(i,j,k-2)
          enddo
       enddo
       ! 6-th order centered filter
       k=4
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i,j,k+1)+d7(2)*rho_n(i,j,k+2)+d7(3)*rho_n(i,j,k+3) &
                  + d7(1)*rho_n(i,j,k-1)+d7(2)*rho_n(i,j,k-2)+d7(3)*rho_n(i,j,k-3)
             Krhou(i,j,k)= Krhou(i,j,k)  +d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i,j,k+1)+d7(2)*rhou_n(i,j,k+2)+d7(3)*rhou_n(i,j,k+3) &
                  + d7(1)*rhou_n(i,j,k-1)+d7(2)*rhou_n(i,j,k-2)+d7(3)*rhou_n(i,j,k-3)
             Krhov(i,j,k)= Krhov(i,j,k)  +d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i,j,k+1)+d7(2)*rhov_n(i,j,k+2)+d7(3)*rhov_n(i,j,k+3) &
                  + d7(1)*rhov_n(i,j,k-1)+d7(2)*rhov_n(i,j,k-2)+d7(3)*rhov_n(i,j,k-3)
             Krhow(i,j,k)= Krhow(i,j,k)  +d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i,j,k+1)+d7(2)*rhow_n(i,j,k+2)+d7(3)*rhow_n(i,j,k+3) &
                  + d7(1)*rhow_n(i,j,k-1)+d7(2)*rhow_n(i,j,k-2)+d7(3)*rhow_n(i,j,k-3)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i,j,k+1)+d7(2)*rhoe_n(i,j,k+2)+d7(3)*rhoe_n(i,j,k+3) &
                  + d7(1)*rhoe_n(i,j,k-1)+d7(2)*rhoe_n(i,j,k-2)+d7(3)*rhoe_n(i,j,k-3)
          enddo
       enddo
       ! 8-th order centered filter
       k=5
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d9(0)*rho_n(i,j,k)   &
                  + d9(1)*rho_n(i,j,k+1)+d9(2)*rho_n(i,j,k+2)+d9(3)*rho_n(i,j,k+3)+d9(4)*rho_n(i,j,k+4) &
                  + d9(1)*rho_n(i,j,k-1)+d9(2)*rho_n(i,j,k-2)+d9(3)*rho_n(i,j,k-3)+d9(4)*rho_n(i,j,k-4)
             Krhou(i,j,k)= Krhou(i,j,k)  +d9(0)*rhou_n(i,j,k)   &
                  + d9(1)*rhou_n(i,j,k+1)+d9(2)*rhou_n(i,j,k+2)+d9(3)*rhou_n(i,j,k+3)+d9(4)*rhou_n(i,j,k+4) &
                  + d9(1)*rhou_n(i,j,k-1)+d9(2)*rhou_n(i,j,k-2)+d9(3)*rhou_n(i,j,k-3)+d9(4)*rhou_n(i,j,k-4)
             Krhov(i,j,k)= Krhov(i,j,k)  +d9(0)*rhov_n(i,j,k)   &
                  + d9(1)*rhov_n(i,j,k+1)+d9(2)*rhov_n(i,j,k+2)+d9(3)*rhov_n(i,j,k+3)+d9(4)*rhov_n(i,j,k+4) &
                  + d9(1)*rhov_n(i,j,k-1)+d9(2)*rhov_n(i,j,k-2)+d9(3)*rhov_n(i,j,k-3)+d9(4)*rhov_n(i,j,k-4)
             Krhow(i,j,k)= Krhow(i,j,k)  +d9(0)*rhow_n(i,j,k)   &
                  + d9(1)*rhow_n(i,j,k+1)+d9(2)*rhow_n(i,j,k+2)+d9(3)*rhow_n(i,j,k+3)+d9(4)*rhow_n(i,j,k+4) &
                  + d9(1)*rhow_n(i,j,k-1)+d9(2)*rhow_n(i,j,k-2)+d9(3)*rhow_n(i,j,k-3)+d9(4)*rhow_n(i,j,k-4)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d9(0)*rhoe_n(i,j,k)   &
                  + d9(1)*rhoe_n(i,j,k+1)+d9(2)*rhoe_n(i,j,k+2)+d9(3)*rhoe_n(i,j,k+3)+d9(4)*rhoe_n(i,j,k+4) &
                  + d9(1)*rhoe_n(i,j,k-1)+d9(2)*rhoe_n(i,j,k-2)+d9(3)*rhoe_n(i,j,k-3)+d9(4)*rhoe_n(i,j,k-4)
          enddo
       enddo
    endif

    ! Interior points
    ! ---------------
    ! 10-th order centered filter
    do k=ndz,nfz
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  + d11(0)*rho_n(i,j,k) &
                  + d11(1)*( rho_n(i,j,k+1)+rho_n(i,j,k-1) )&
                  + d11(2)*( rho_n(i,j,k+2)+rho_n(i,j,k-2) )&
                  + d11(3)*( rho_n(i,j,k+3)+rho_n(i,j,k-3) )&
                  + d11(4)*( rho_n(i,j,k+4)+rho_n(i,j,k-4) )&
                  + d11(5)*( rho_n(i,j,k+5)+rho_n(i,j,k-5) )

             Krhou(i,j,k)= Krhou(i,j,k)  + d11(0)*rhou_n(i,j,k) &
                  + d11(1)*( rhou_n(i,j,k+1)+rhou_n(i,j,k-1) )&
                  + d11(2)*( rhou_n(i,j,k+2)+rhou_n(i,j,k-2) )&
                  + d11(3)*( rhou_n(i,j,k+3)+rhou_n(i,j,k-3) )&
                  + d11(4)*( rhou_n(i,j,k+4)+rhou_n(i,j,k-4) )&
                  + d11(5)*( rhou_n(i,j,k+5)+rhou_n(i,j,k-5) )

             Krhov(i,j,k)= Krhov(i,j,k)  + d11(0)*rhov_n(i,j,k) &
                  + d11(1)*( rhov_n(i,j,k+1)+rhov_n(i,j,k-1) )&
                  + d11(2)*( rhov_n(i,j,k+2)+rhov_n(i,j,k-2) )&
                  + d11(3)*( rhov_n(i,j,k+3)+rhov_n(i,j,k-3) )&
                  + d11(4)*( rhov_n(i,j,k+4)+rhov_n(i,j,k-4) )&
                  + d11(5)*( rhov_n(i,j,k+5)+rhov_n(i,j,k-5) )

             Krhow(i,j,k)= Krhow(i,j,k)  + d11(0)*rhow_n(i,j,k) &
                  + d11(1)*( rhow_n(i,j,k+1)+rhow_n(i,j,k-1) )&
                  + d11(2)*( rhow_n(i,j,k+2)+rhow_n(i,j,k-2) )&
                  + d11(3)*( rhow_n(i,j,k+3)+rhow_n(i,j,k-3) )&
                  + d11(4)*( rhow_n(i,j,k+4)+rhow_n(i,j,k-4) )&
                  + d11(5)*( rhow_n(i,j,k+5)+rhow_n(i,j,k-5) )

             Krhoe(i,j,k)= Krhoe(i,j,k)  + d11(0)*rhoe_n(i,j,k) &
                  + d11(1)*( rhoe_n(i,j,k+1)+rhoe_n(i,j,k-1) )&
                  + d11(2)*( rhoe_n(i,j,k+2)+rhoe_n(i,j,k-2) )&
                  + d11(3)*( rhoe_n(i,j,k+3)+rhoe_n(i,j,k-3) )&
                  + d11(4)*( rhoe_n(i,j,k+4)+rhoe_n(i,j,k-4) )&
                  + d11(5)*( rhoe_n(i,j,k+5)+rhoe_n(i,j,k-5) )
          enddo
       enddo
    enddo

    ! BC at kmax
    ! ----------
    if (is_boundary(3,2)) then
       ! 8-th order centered filter
       k=nz-4
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d9(0)*rho_n(i,j,k)   &
                  + d9(1)*rho_n(i,j,k+1)+d9(2)*rho_n(i,j,k+2)+d9(3)*rho_n(i,j,k+3)+d9(4)*rho_n(i,j,k+4) &
                  + d9(1)*rho_n(i,j,k-1)+d9(2)*rho_n(i,j,k-2)+d9(3)*rho_n(i,j,k-3)+d9(4)*rho_n(i,j,k-4)
             Krhou(i,j,k)= Krhou(i,j,k)  +d9(0)*rhou_n(i,j,k)   &
                  + d9(1)*rhou_n(i,j,k+1)+d9(2)*rhou_n(i,j,k+2)+d9(3)*rhou_n(i,j,k+3)+d9(4)*rhou_n(i,j,k+4) &
                  + d9(1)*rhou_n(i,j,k-1)+d9(2)*rhou_n(i,j,k-2)+d9(3)*rhou_n(i,j,k-3)+d9(4)*rhou_n(i,j,k-4)
             Krhov(i,j,k)= Krhov(i,j,k)  +d9(0)*rhov_n(i,j,k)   &
                  + d9(1)*rhov_n(i,j,k+1)+d9(2)*rhov_n(i,j,k+2)+d9(3)*rhov_n(i,j,k+3)+d9(4)*rhov_n(i,j,k+4) &
                  + d9(1)*rhov_n(i,j,k-1)+d9(2)*rhov_n(i,j,k-2)+d9(3)*rhov_n(i,j,k-3)+d9(4)*rhov_n(i,j,k-4)
             Krhow(i,j,k)= Krhow(i,j,k)  +d9(0)*rhow_n(i,j,k)   &
                  + d9(1)*rhow_n(i,j,k+1)+d9(2)*rhow_n(i,j,k+2)+d9(3)*rhow_n(i,j,k+3)+d9(4)*rhow_n(i,j,k+4) &
                  + d9(1)*rhow_n(i,j,k-1)+d9(2)*rhow_n(i,j,k-2)+d9(3)*rhow_n(i,j,k-3)+d9(4)*rhow_n(i,j,k-4)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d9(0)*rhoe_n(i,j,k)   &
                  + d9(1)*rhoe_n(i,j,k+1)+d9(2)*rhoe_n(i,j,k+2)+d9(3)*rhoe_n(i,j,k+3)+d9(4)*rhoe_n(i,j,k+4) &
                  + d9(1)*rhoe_n(i,j,k-1)+d9(2)*rhoe_n(i,j,k-2)+d9(3)*rhoe_n(i,j,k-3)+d9(4)*rhoe_n(i,j,k-4)
          enddo
       enddo
       ! 6-th order centered filter
       k=nz-3
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i,j,k+1)+d7(2)*rho_n(i,j,k+2)+d7(3)*rho_n(i,j,k+3) &
                  + d7(1)*rho_n(i,j,k-1)+d7(2)*rho_n(i,j,k-2)+d7(3)*rho_n(i,j,k-3)
             Krhou(i,j,k)= Krhou(i,j,k)  +d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i,j,k+1)+d7(2)*rhou_n(i,j,k+2)+d7(3)*rhou_n(i,j,k+3) &
                  + d7(1)*rhou_n(i,j,k-1)+d7(2)*rhou_n(i,j,k-2)+d7(3)*rhou_n(i,j,k-3)
             Krhov(i,j,k)= Krhov(i,j,k)  +d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i,j,k+1)+d7(2)*rhov_n(i,j,k+2)+d7(3)*rhov_n(i,j,k+3) &
                  + d7(1)*rhov_n(i,j,k-1)+d7(2)*rhov_n(i,j,k-2)+d7(3)*rhov_n(i,j,k-3)
             Krhow(i,j,k)= Krhow(i,j,k)  +d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i,j,k+1)+d7(2)*rhow_n(i,j,k+2)+d7(3)*rhow_n(i,j,k+3) &
                  + d7(1)*rhow_n(i,j,k-1)+d7(2)*rhow_n(i,j,k-2)+d7(3)*rhow_n(i,j,k-3)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i,j,k+1)+d7(2)*rhoe_n(i,j,k+2)+d7(3)*rhoe_n(i,j,k+3) &
                  + d7(1)*rhoe_n(i,j,k-1)+d7(2)*rhoe_n(i,j,k-2)+d7(3)*rhoe_n(i,j,k-3)
          enddo
       enddo
       ! 4-th order centered filter
       k=nz-2
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i,j,k+1)+d5(2)*rho_n(i,j,k+2) &
                  + d5(1)*rho_n(i,j,k-1)+d5(2)*rho_n(i,j,k-2)
             Krhou(i,j,k)= Krhou(i,j,k)  +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j,k+1)+d5(2)*rhou_n(i,j,k+2) &
                  + d5(1)*rhou_n(i,j,k-1)+d5(2)*rhou_n(i,j,k-2)
             Krhov(i,j,k)= Krhov(i,j,k)  +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j,k+1)+d5(2)*rhov_n(i,j,k+2) &
                  + d5(1)*rhov_n(i,j,k-1)+d5(2)*rhov_n(i,j,k-2)
             Krhow(i,j,k)= Krhow(i,j,k)  +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j,k+1)+d5(2)*rhow_n(i,j,k+2) &
                  + d5(1)*rhow_n(i,j,k-1)+d5(2)*rhow_n(i,j,k-2)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i,j,k+1)+d5(2)*rhoe_n(i,j,k+2) &
                  + d5(1)*rhoe_n(i,j,k-1)+d5(2)*rhoe_n(i,j,k-2)
          enddo
       enddo
       ! no filter for k=nz-1 & k=nz
    endif
        
  end subroutine filtering_11pts_

  !==============================================================================
  subroutine filtering_11pts
  !==============================================================================
    !> Apply filtering on 11-point stencils (+ boundaries)
    !> * old version without comput./comm. overlap *
  !==============================================================================
    use mod_flow
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    ! ---------------------------------------------------------------------------
    
    ! Filtering along i-direction
    ! ===========================
    
    ! BC at imin
    ! ----------
    if (is_boundary(1,1)) then
       ! no filter for i=1 & i=2 -> initialize array
       do j=1,ny
           Krho(1,j,:)=0.0_wp
          Krhou(1,j,:)=0.0_wp
          Krhov(1,j,:)=0.0_wp
          Krhow(1,j,:)=0.0_wp
          Krhoe(1,j,:)=0.0_wp
          
           Krho(2,j,:)=0.0_wp
          Krhou(2,j,:)=0.0_wp
          Krhov(2,j,:)=0.0_wp
          Krhow(2,j,:)=0.0_wp
          Krhoe(2,j,:)=0.0_wp
       enddo
       ! 4-th order centered filter
       i=3
       do k=1,nz
          do j=1,ny
             Krho(i,j,k) = d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i+1,j,k)+d5(2)*rho_n(i+2,j,k) &
                  + d5(1)*rho_n(i-1,j,k)+d5(2)*rho_n(i-2,j,k)
             Krhou(i,j,k)= d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i+1,j,k)+d5(2)*rhou_n(i+2,j,k) &
                  + d5(1)*rhou_n(i-1,j,k)+d5(2)*rhou_n(i-2,j,k)
             Krhov(i,j,k)= d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i+1,j,k)+d5(2)*rhov_n(i+2,j,k) &
                  + d5(1)*rhov_n(i-1,j,k)+d5(2)*rhov_n(i-2,j,k)
             Krhow(i,j,k)= d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i+1,j,k)+d5(2)*rhow_n(i+2,j,k) &
                  + d5(1)*rhow_n(i-1,j,k)+d5(2)*rhow_n(i-2,j,k)
             Krhoe(i,j,k)= d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i+1,j,k)+d5(2)*rhoe_n(i+2,j,k) &
                  + d5(1)*rhoe_n(i-1,j,k)+d5(2)*rhoe_n(i-2,j,k)
          enddo
       enddo
       ! 6-th order centered filter
       i=4
       do k=1,nz
          do j=1,ny
             Krho(i,j,k) = d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i+1,j,k)+d7(2)*rho_n(i+2,j,k)+d7(3)*rho_n(i+3,j,k) &
                  + d7(1)*rho_n(i-1,j,k)+d7(2)*rho_n(i-2,j,k)+d7(3)*rho_n(i-3,j,k)
             Krhou(i,j,k)= d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i+1,j,k)+d7(2)*rhou_n(i+2,j,k)+d7(3)*rhou_n(i+3,j,k) &
                  + d7(1)*rhou_n(i-1,j,k)+d7(2)*rhou_n(i-2,j,k)+d7(3)*rhou_n(i-3,j,k)
             Krhov(i,j,k)= d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i+1,j,k)+d7(2)*rhov_n(i+2,j,k)+d7(3)*rhov_n(i+3,j,k) &
                  + d7(1)*rhov_n(i-1,j,k)+d7(2)*rhov_n(i-2,j,k)+d7(3)*rhov_n(i-3,j,k)
             Krhow(i,j,k)= d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i+1,j,k)+d7(2)*rhow_n(i+2,j,k)+d7(3)*rhow_n(i+3,j,k) &
                  + d7(1)*rhow_n(i-1,j,k)+d7(2)*rhow_n(i-2,j,k)+d7(3)*rhow_n(i-3,j,k)
             Krhoe(i,j,k)= d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i+1,j,k)+d7(2)*rhoe_n(i+2,j,k)+d7(3)*rhoe_n(i+3,j,k) &
                  + d7(1)*rhoe_n(i-1,j,k)+d7(2)*rhoe_n(i-2,j,k)+d7(3)*rhoe_n(i-3,j,k)
          enddo
       enddo
       ! 8-th order centered filter
       i=5
       do k=1,nz
          do j=1,ny
             Krho(i,j,k) = d9(0)*rho_n(i,j,k)   &
                  + d9(1)*rho_n(i+1,j,k)+d9(2)*rho_n(i+2,j,k)+d9(3)*rho_n(i+3,j,k)+d9(4)*rho_n(i+4,j,k) &
                  + d9(1)*rho_n(i-1,j,k)+d9(2)*rho_n(i-2,j,k)+d9(3)*rho_n(i-3,j,k)+d9(4)*rho_n(i-4,j,k)
             Krhou(i,j,k)= d9(0)*rhou_n(i,j,k)   &
                  + d9(1)*rhou_n(i+1,j,k)+d9(2)*rhou_n(i+2,j,k)+d9(3)*rhou_n(i+3,j,k)+d9(4)*rhou_n(i+4,j,k) &
                  + d9(1)*rhou_n(i-1,j,k)+d9(2)*rhou_n(i-2,j,k)+d9(3)*rhou_n(i-3,j,k)+d9(4)*rhou_n(i-4,j,k)
             Krhov(i,j,k)= d9(0)*rhov_n(i,j,k)   &
                  + d9(1)*rhov_n(i+1,j,k)+d9(2)*rhov_n(i+2,j,k)+d9(3)*rhov_n(i+3,j,k)+d9(4)*rhov_n(i+4,j,k) &
                  + d9(1)*rhov_n(i-1,j,k)+d9(2)*rhov_n(i-2,j,k)+d9(3)*rhov_n(i-3,j,k)+d9(4)*rhov_n(i-4,j,k)
             Krhow(i,j,k)= d9(0)*rhow_n(i,j,k)   &
                  + d9(1)*rhow_n(i+1,j,k)+d9(2)*rhow_n(i+2,j,k)+d9(3)*rhow_n(i+3,j,k)+d9(4)*rhow_n(i+4,j,k) &
                  + d9(1)*rhow_n(i-1,j,k)+d9(2)*rhow_n(i-2,j,k)+d9(3)*rhow_n(i-3,j,k)+d9(4)*rhow_n(i-4,j,k)
             Krhoe(i,j,k)= d9(0)*rhoe_n(i,j,k)   &
                  + d9(1)*rhoe_n(i+1,j,k)+d9(2)*rhoe_n(i+2,j,k)+d9(3)*rhoe_n(i+3,j,k)+d9(4)*rhoe_n(i+4,j,k) &
                  + d9(1)*rhoe_n(i-1,j,k)+d9(2)*rhoe_n(i-2,j,k)+d9(3)*rhoe_n(i-3,j,k)+d9(4)*rhoe_n(i-4,j,k)
          enddo
       enddo
    endif

    ! Interior points
    ! ---------------
    ! 10-th order centered filter
    do k=1,nz
       do j=1,ny
          do i=ndx,nfx
             Krho(i,j,k) = d11(0)*rho_n(i,j,k) &
                  + d11(1)*( rho_n(i+1,j,k)+rho_n(i-1,j,k) )&
                  + d11(2)*( rho_n(i+2,j,k)+rho_n(i-2,j,k) )&
                  + d11(3)*( rho_n(i+3,j,k)+rho_n(i-3,j,k) )&
                  + d11(4)*( rho_n(i+4,j,k)+rho_n(i-4,j,k) )&
                  + d11(5)*( rho_n(i+5,j,k)+rho_n(i-5,j,k) )

             Krhou(i,j,k)= d11(0)*rhou_n(i,j,k) &
                  + d11(1)*( rhou_n(i+1,j,k)+rhou_n(i-1,j,k) )&
                  + d11(2)*( rhou_n(i+2,j,k)+rhou_n(i-2,j,k) )&
                  + d11(3)*( rhou_n(i+3,j,k)+rhou_n(i-3,j,k) )&
                  + d11(4)*( rhou_n(i+4,j,k)+rhou_n(i-4,j,k) )&
                  + d11(5)*( rhou_n(i+5,j,k)+rhou_n(i-5,j,k) )

             Krhov(i,j,k)= d11(0)*rhov_n(i,j,k) &
                  + d11(1)*( rhov_n(i+1,j,k)+rhov_n(i-1,j,k) )&
                  + d11(2)*( rhov_n(i+2,j,k)+rhov_n(i-2,j,k) )&
                  + d11(3)*( rhov_n(i+3,j,k)+rhov_n(i-3,j,k) )&
                  + d11(4)*( rhov_n(i+4,j,k)+rhov_n(i-4,j,k) )&
                  + d11(5)*( rhov_n(i+5,j,k)+rhov_n(i-5,j,k) )

             Krhow(i,j,k)= d11(0)*rhow_n(i,j,k) &
                  + d11(1)*( rhow_n(i+1,j,k)+rhow_n(i-1,j,k) )&
                  + d11(2)*( rhow_n(i+2,j,k)+rhow_n(i-2,j,k) )&
                  + d11(3)*( rhow_n(i+3,j,k)+rhow_n(i-3,j,k) )&
                  + d11(4)*( rhow_n(i+4,j,k)+rhow_n(i-4,j,k) )&
                  + d11(5)*( rhow_n(i+5,j,k)+rhow_n(i-5,j,k) )

             Krhoe(i,j,k)= d11(0)*rhoe_n(i,j,k) &
                  + d11(1)*( rhoe_n(i+1,j,k)+rhoe_n(i-1,j,k) )&
                  + d11(2)*( rhoe_n(i+2,j,k)+rhoe_n(i-2,j,k) )&
                  + d11(3)*( rhoe_n(i+3,j,k)+rhoe_n(i-3,j,k) )&
                  + d11(4)*( rhoe_n(i+4,j,k)+rhoe_n(i-4,j,k) )&
                  + d11(5)*( rhoe_n(i+5,j,k)+rhoe_n(i-5,j,k) )
          enddo
       enddo
    enddo

    ! BC at imax
    ! ----------
    if (is_boundary(1,2)) then
       ! 8-th order centered filter
       i=nx-4
       do k=1,nz
          do j=1,ny
             Krho(i,j,k) = d9(0)*rho_n(i,j,k)   &
                  + d9(1)*rho_n(i+1,j,k)+d9(2)*rho_n(i+2,j,k)+d9(3)*rho_n(i+3,j,k)+d9(4)*rho_n(i+4,j,k) &
                  + d9(1)*rho_n(i-1,j,k)+d9(2)*rho_n(i-2,j,k)+d9(3)*rho_n(i-3,j,k)+d9(4)*rho_n(i-4,j,k)
             Krhou(i,j,k)= d9(0)*rhou_n(i,j,k)   &
                  + d9(1)*rhou_n(i+1,j,k)+d9(2)*rhou_n(i+2,j,k)+d9(3)*rhou_n(i+3,j,k)+d9(4)*rhou_n(i+4,j,k) &
                  + d9(1)*rhou_n(i-1,j,k)+d9(2)*rhou_n(i-2,j,k)+d9(3)*rhou_n(i-3,j,k)+d9(4)*rhou_n(i-4,j,k)
             Krhov(i,j,k)= d9(0)*rhov_n(i,j,k)   &
                  + d9(1)*rhov_n(i+1,j,k)+d9(2)*rhov_n(i+2,j,k)+d9(3)*rhov_n(i+3,j,k)+d9(4)*rhov_n(i+4,j,k) &
                  + d9(1)*rhov_n(i-1,j,k)+d9(2)*rhov_n(i-2,j,k)+d9(3)*rhov_n(i-3,j,k)+d9(4)*rhov_n(i-4,j,k)
             Krhow(i,j,k)= d9(0)*rhow_n(i,j,k)   &
                  + d9(1)*rhow_n(i+1,j,k)+d9(2)*rhow_n(i+2,j,k)+d9(3)*rhow_n(i+3,j,k)+d9(4)*rhow_n(i+4,j,k) &
                  + d9(1)*rhow_n(i-1,j,k)+d9(2)*rhow_n(i-2,j,k)+d9(3)*rhow_n(i-3,j,k)+d9(4)*rhow_n(i-4,j,k)
             Krhoe(i,j,k)= d9(0)*rhoe_n(i,j,k)   &
                  + d9(1)*rhoe_n(i+1,j,k)+d9(2)*rhoe_n(i+2,j,k)+d9(3)*rhoe_n(i+3,j,k)+d9(4)*rhoe_n(i+4,j,k) &
                  + d9(1)*rhoe_n(i-1,j,k)+d9(2)*rhoe_n(i-2,j,k)+d9(3)*rhoe_n(i-3,j,k)+d9(4)*rhoe_n(i-4,j,k)
          enddo
       enddo
       ! 6-th order centered filter
       i=nx-3
       do k=1,nz
          do j=1,ny
             Krho(i,j,k) = d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i+1,j,k)+d7(2)*rho_n(i+2,j,k)+d7(3)*rho_n(i+3,j,k) &
                  + d7(1)*rho_n(i-1,j,k)+d7(2)*rho_n(i-2,j,k)+d7(3)*rho_n(i-3,j,k)
             Krhou(i,j,k)= d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i+1,j,k)+d7(2)*rhou_n(i+2,j,k)+d7(3)*rhou_n(i+3,j,k) &
                  + d7(1)*rhou_n(i-1,j,k)+d7(2)*rhou_n(i-2,j,k)+d7(3)*rhou_n(i-3,j,k)
             Krhov(i,j,k)= d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i+1,j,k)+d7(2)*rhov_n(i+2,j,k)+d7(3)*rhov_n(i+3,j,k) &
                  + d7(1)*rhov_n(i-1,j,k)+d7(2)*rhov_n(i-2,j,k)+d7(3)*rhov_n(i-3,j,k)
             Krhow(i,j,k)= d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i+1,j,k)+d7(2)*rhow_n(i+2,j,k)+d7(3)*rhow_n(i+3,j,k) &
                  + d7(1)*rhow_n(i-1,j,k)+d7(2)*rhow_n(i-2,j,k)+d7(3)*rhow_n(i-3,j,k)
             Krhoe(i,j,k)= d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i+1,j,k)+d7(2)*rhoe_n(i+2,j,k)+d7(3)*rhoe_n(i+3,j,k) &
                  + d7(1)*rhoe_n(i-1,j,k)+d7(2)*rhoe_n(i-2,j,k)+d7(3)*rhoe_n(i-3,j,k)
          enddo
       enddo
       ! 4-th order centered filter
       i=nx-2
       do k=1,nz
          do j=1,ny
             Krho(i,j,k) = d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i+1,j,k)+d5(2)*rho_n(i+2,j,k) &
                  + d5(1)*rho_n(i-1,j,k)+d5(2)*rho_n(i-2,j,k)
             Krhou(i,j,k)= d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i+1,j,k)+d5(2)*rhou_n(i+2,j,k) &
                  + d5(1)*rhou_n(i-1,j,k)+d5(2)*rhou_n(i-2,j,k)
             Krhov(i,j,k)= d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i+1,j,k)+d5(2)*rhov_n(i+2,j,k) &
                  + d5(1)*rhov_n(i-1,j,k)+d5(2)*rhov_n(i-2,j,k)
             Krhow(i,j,k)= d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i+1,j,k)+d5(2)*rhow_n(i+2,j,k) &
                  + d5(1)*rhow_n(i-1,j,k)+d5(2)*rhow_n(i-2,j,k)
             Krhoe(i,j,k)= d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i+1,j,k)+d5(2)*rhoe_n(i+2,j,k) &
                  + d5(1)*rhoe_n(i-1,j,k)+d5(2)*rhoe_n(i-2,j,k)
          enddo
       enddo
       ! no filter for i=nx-1 & i=nx -> initialize array
       do j=1,ny
           Krho(nx-1,j,:)=0.0_wp
          Krhou(nx-1,j,:)=0.0_wp
          Krhov(nx-1,j,:)=0.0_wp
          Krhow(nx-1,j,:)=0.0_wp
          Krhoe(nx-1,j,:)=0.0_wp
          
           Krho(nx,j,:)=0.0_wp
          Krhou(nx,j,:)=0.0_wp
          Krhov(nx,j,:)=0.0_wp
          Krhow(nx,j,:)=0.0_wp
          Krhoe(nx,j,:)=0.0_wp
       enddo
    endif
    
    ! Filtering along j-direction
    ! ===========================
    
    ! BC at jmin
    ! ----------
    if (is_boundary(2,1)) then
       ! no filter for j=1 & j=2
       
       ! 4-th order centered filter
       j=3
       do k=1,nz
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i,j+1,k)+d5(2)*rho_n(i,j+2,k) &
                  + d5(1)*rho_n(i,j-1,k)+d5(2)*rho_n(i,j-2,k)
             Krhou(i,j,k)= Krhou(i,j,k)  +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j+1,k)+d5(2)*rhou_n(i,j+2,k) &
                  + d5(1)*rhou_n(i,j-1,k)+d5(2)*rhou_n(i,j-2,k)
             Krhov(i,j,k)= Krhov(i,j,k)  +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j+1,k)+d5(2)*rhov_n(i,j+2,k) &
                  + d5(1)*rhov_n(i,j-1,k)+d5(2)*rhov_n(i,j-2,k)
             Krhow(i,j,k)= Krhow(i,j,k)  +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j+1,k)+d5(2)*rhow_n(i,j+2,k) &
                  + d5(1)*rhow_n(i,j-1,k)+d5(2)*rhow_n(i,j-2,k)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i,j+1,k)+d5(2)*rhoe_n(i,j+2,k) &
                  + d5(1)*rhoe_n(i,j-1,k)+d5(2)*rhoe_n(i,j-2,k)
          enddo
       enddo
       ! 6-th order centered filter
       j=4
       do k=1,nz
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i,j+1,k)+d7(2)*rho_n(i,j+2,k)+d7(3)*rho_n(i,j+3,k) &
                  + d7(1)*rho_n(i,j-1,k)+d7(2)*rho_n(i,j-2,k)+d7(3)*rho_n(i,j-3,k)
             Krhou(i,j,k)= Krhou(i,j,k)  +d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i,j+1,k)+d7(2)*rhou_n(i,j+2,k)+d7(3)*rhou_n(i,j+3,k) &
                  + d7(1)*rhou_n(i,j-1,k)+d7(2)*rhou_n(i,j-2,k)+d7(3)*rhou_n(i,j-3,k)
             Krhov(i,j,k)= Krhov(i,j,k)  +d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i,j+1,k)+d7(2)*rhov_n(i,j+2,k)+d7(3)*rhov_n(i,j+3,k) &
                  + d7(1)*rhov_n(i,j-1,k)+d7(2)*rhov_n(i,j-2,k)+d7(3)*rhov_n(i,j-3,k)
             Krhow(i,j,k)= Krhow(i,j,k)  +d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i,j+1,k)+d7(2)*rhow_n(i,j+2,k)+d7(3)*rhow_n(i,j+3,k) &
                  + d7(1)*rhow_n(i,j-1,k)+d7(2)*rhow_n(i,j-2,k)+d7(3)*rhow_n(i,j-3,k)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i,j+1,k)+d7(2)*rhoe_n(i,j+2,k)+d7(3)*rhoe_n(i,j+3,k) &
                  + d7(1)*rhoe_n(i,j-1,k)+d7(2)*rhoe_n(i,j-2,k)+d7(3)*rhoe_n(i,j-3,k)
          enddo
       enddo
       ! 8-th order centered filter
       j=5
       do k=1,nz
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d9(0)*rho_n(i,j,k)   &
                  + d9(1)*rho_n(i,j+1,k)+d9(2)*rho_n(i,j+2,k)+d9(3)*rho_n(i,j+3,k)+d9(4)*rho_n(i,j+4,k) &
                  + d9(1)*rho_n(i,j-1,k)+d9(2)*rho_n(i,j-2,k)+d9(3)*rho_n(i,j-3,k)+d9(4)*rho_n(i,j-4,k)
             Krhou(i,j,k)= Krhou(i,j,k)  +d9(0)*rhou_n(i,j,k)   &
                  + d9(1)*rhou_n(i,j+1,k)+d9(2)*rhou_n(i,j+2,k)+d9(3)*rhou_n(i,j+3,k)+d9(4)*rhou_n(i,j+4,k) &
                  + d9(1)*rhou_n(i,j-1,k)+d9(2)*rhou_n(i,j-2,k)+d9(3)*rhou_n(i,j-3,k)+d9(4)*rhou_n(i,j-4,k)
             Krhov(i,j,k)= Krhov(i,j,k)  +d9(0)*rhov_n(i,j,k)   &
                  + d9(1)*rhov_n(i,j+1,k)+d9(2)*rhov_n(i,j+2,k)+d9(3)*rhov_n(i,j+3,k)+d9(4)*rhov_n(i,j+4,k) &
                  + d9(1)*rhov_n(i,j-1,k)+d9(2)*rhov_n(i,j-2,k)+d9(3)*rhov_n(i,j-3,k)+d9(4)*rhov_n(i,j-4,k)
             Krhow(i,j,k)= Krhow(i,j,k)  +d9(0)*rhow_n(i,j,k)   &
                  + d9(1)*rhow_n(i,j+1,k)+d9(2)*rhow_n(i,j+2,k)+d9(3)*rhow_n(i,j+3,k)+d9(4)*rhow_n(i,j+4,k) &
                  + d9(1)*rhow_n(i,j-1,k)+d9(2)*rhow_n(i,j-2,k)+d9(3)*rhow_n(i,j-3,k)+d9(4)*rhow_n(i,j-4,k)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d9(0)*rhoe_n(i,j,k)   &
                  + d9(1)*rhoe_n(i,j+1,k)+d9(2)*rhoe_n(i,j+2,k)+d9(3)*rhoe_n(i,j+3,k)+d9(4)*rhoe_n(i,j+4,k) &
                  + d9(1)*rhoe_n(i,j-1,k)+d9(2)*rhoe_n(i,j-2,k)+d9(3)*rhoe_n(i,j-3,k)+d9(4)*rhoe_n(i,j-4,k)
          enddo
       enddo
    endif

    ! Interior points
    ! ---------------
    ! 10-th order centered filter
    do k=1,nz
       do j=ndy,nfy
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  + d11(0)*rho_n(i,j,k) &
                  + d11(1)*( rho_n(i,j+1,k)+rho_n(i,j-1,k) )&
                  + d11(2)*( rho_n(i,j+2,k)+rho_n(i,j-2,k) )&
                  + d11(3)*( rho_n(i,j+3,k)+rho_n(i,j-3,k) )&
                  + d11(4)*( rho_n(i,j+4,k)+rho_n(i,j-4,k) )&
                  + d11(5)*( rho_n(i,j+5,k)+rho_n(i,j-5,k) )

             Krhou(i,j,k)= Krhou(i,j,k)  + d11(0)*rhou_n(i,j,k) &
                  + d11(1)*( rhou_n(i,j+1,k)+rhou_n(i,j-1,k) )&
                  + d11(2)*( rhou_n(i,j+2,k)+rhou_n(i,j-2,k) )&
                  + d11(3)*( rhou_n(i,j+3,k)+rhou_n(i,j-3,k) )&
                  + d11(4)*( rhou_n(i,j+4,k)+rhou_n(i,j-4,k) )&
                  + d11(5)*( rhou_n(i,j+5,k)+rhou_n(i,j-5,k) )

             Krhov(i,j,k)= Krhov(i,j,k)  + d11(0)*rhov_n(i,j,k) &
                  + d11(1)*( rhov_n(i,j+1,k)+rhov_n(i,j-1,k) )&
                  + d11(2)*( rhov_n(i,j+2,k)+rhov_n(i,j-2,k) )&
                  + d11(3)*( rhov_n(i,j+3,k)+rhov_n(i,j-3,k) )&
                  + d11(4)*( rhov_n(i,j+4,k)+rhov_n(i,j-4,k) )&
                  + d11(5)*( rhov_n(i,j+5,k)+rhov_n(i,j-5,k) )

             Krhow(i,j,k)= Krhow(i,j,k)  + d11(0)*rhow_n(i,j,k) &
                  + d11(1)*( rhow_n(i,j+1,k)+rhow_n(i,j-1,k) )&
                  + d11(2)*( rhow_n(i,j+2,k)+rhow_n(i,j-2,k) )&
                  + d11(3)*( rhow_n(i,j+3,k)+rhow_n(i,j-3,k) )&
                  + d11(4)*( rhow_n(i,j+4,k)+rhow_n(i,j-4,k) )&
                  + d11(5)*( rhow_n(i,j+5,k)+rhow_n(i,j-5,k) )

             Krhoe(i,j,k)= Krhoe(i,j,k)  + d11(0)*rhoe_n(i,j,k) &
                  + d11(1)*( rhoe_n(i,j+1,k)+rhoe_n(i,j-1,k) )&
                  + d11(2)*( rhoe_n(i,j+2,k)+rhoe_n(i,j-2,k) )&
                  + d11(3)*( rhoe_n(i,j+3,k)+rhoe_n(i,j-3,k) )&
                  + d11(4)*( rhoe_n(i,j+4,k)+rhoe_n(i,j-4,k) )&
                  + d11(5)*( rhoe_n(i,j+5,k)+rhoe_n(i,j-5,k) )
          enddo
       enddo
    enddo

    ! BC at jmax
    ! ----------
    if (is_boundary(2,2)) then
       ! 8-th order centered filter
       j=ny-4
       do k=1,nz
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d9(0)*rho_n(i,j,k)   &
                  + d9(1)*rho_n(i,j+1,k)+d9(2)*rho_n(i,j+2,k)+d9(3)*rho_n(i,j+3,k)+d9(4)*rho_n(i,j+4,k) &
                  + d9(1)*rho_n(i,j-1,k)+d9(2)*rho_n(i,j-2,k)+d9(3)*rho_n(i,j-3,k)+d9(4)*rho_n(i,j-4,k)
             Krhou(i,j,k)= Krhou(i,j,k)  +d9(0)*rhou_n(i,j,k)   &
                  + d9(1)*rhou_n(i,j+1,k)+d9(2)*rhou_n(i,j+2,k)+d9(3)*rhou_n(i,j+3,k)+d9(4)*rhou_n(i,j+4,k) &
                  + d9(1)*rhou_n(i,j-1,k)+d9(2)*rhou_n(i,j-2,k)+d9(3)*rhou_n(i,j-3,k)+d9(4)*rhou_n(i,j-4,k)
             Krhov(i,j,k)= Krhov(i,j,k)  +d9(0)*rhov_n(i,j,k)   &
                  + d9(1)*rhov_n(i,j+1,k)+d9(2)*rhov_n(i,j+2,k)+d9(3)*rhov_n(i,j+3,k)+d9(4)*rhov_n(i,j+4,k) &
                  + d9(1)*rhov_n(i,j-1,k)+d9(2)*rhov_n(i,j-2,k)+d9(3)*rhov_n(i,j-3,k)+d9(4)*rhov_n(i,j-4,k)
             Krhow(i,j,k)= Krhow(i,j,k)  +d9(0)*rhow_n(i,j,k)   &
                  + d9(1)*rhow_n(i,j+1,k)+d9(2)*rhow_n(i,j+2,k)+d9(3)*rhow_n(i,j+3,k)+d9(4)*rhow_n(i,j+4,k) &
                  + d9(1)*rhow_n(i,j-1,k)+d9(2)*rhow_n(i,j-2,k)+d9(3)*rhow_n(i,j-3,k)+d9(4)*rhow_n(i,j-4,k)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d9(0)*rhoe_n(i,j,k)   &
                  + d9(1)*rhoe_n(i,j+1,k)+d9(2)*rhoe_n(i,j+2,k)+d9(3)*rhoe_n(i,j+3,k)+d9(4)*rhoe_n(i,j+4,k) &
                  + d9(1)*rhoe_n(i,j-1,k)+d9(2)*rhoe_n(i,j-2,k)+d9(3)*rhoe_n(i,j-3,k)+d9(4)*rhoe_n(i,j-4,k)
          enddo
       enddo
       ! 6-th order centered filter
       j=ny-3
       do k=1,nz
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i,j+1,k)+d7(2)*rho_n(i,j+2,k)+d7(3)*rho_n(i,j+3,k) &
                  + d7(1)*rho_n(i,j-1,k)+d7(2)*rho_n(i,j-2,k)+d7(3)*rho_n(i,j-3,k)
             Krhou(i,j,k)= Krhou(i,j,k)  +d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i,j+1,k)+d7(2)*rhou_n(i,j+2,k)+d7(3)*rhou_n(i,j+3,k) &
                  + d7(1)*rhou_n(i,j-1,k)+d7(2)*rhou_n(i,j-2,k)+d7(3)*rhou_n(i,j-3,k)
             Krhov(i,j,k)= Krhov(i,j,k)  +d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i,j+1,k)+d7(2)*rhov_n(i,j+2,k)+d7(3)*rhov_n(i,j+3,k) &
                  + d7(1)*rhov_n(i,j-1,k)+d7(2)*rhov_n(i,j-2,k)+d7(3)*rhov_n(i,j-3,k)
             Krhow(i,j,k)= Krhow(i,j,k)  +d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i,j+1,k)+d7(2)*rhow_n(i,j+2,k)+d7(3)*rhow_n(i,j+3,k) &
                  + d7(1)*rhow_n(i,j-1,k)+d7(2)*rhow_n(i,j-2,k)+d7(3)*rhow_n(i,j-3,k)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i,j+1,k)+d7(2)*rhoe_n(i,j+2,k)+d7(3)*rhoe_n(i,j+3,k) &
                  + d7(1)*rhoe_n(i,j-1,k)+d7(2)*rhoe_n(i,j-2,k)+d7(3)*rhoe_n(i,j-3,k)
          enddo
       enddo
       ! 4-th order centered filter
       j=ny-2
       do k=1,nz
          do i=1,nx
             Krho(i,j,k)= Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i,j+1,k)+d5(2)*rho_n(i,j+2,k) &
                  + d5(1)*rho_n(i,j-1,k)+d5(2)*rho_n(i,j-2,k)
             Krhou(i,j,k)= Krhou(i,j,k)  +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j+1,k)+d5(2)*rhou_n(i,j+2,k) &
                  + d5(1)*rhou_n(i,j-1,k)+d5(2)*rhou_n(i,j-2,k)
             Krhov(i,j,k)= Krhov(i,j,k)  +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j+1,k)+d5(2)*rhov_n(i,j+2,k) &
                  + d5(1)*rhov_n(i,j-1,k)+d5(2)*rhov_n(i,j-2,k)
             Krhow(i,j,k)= Krhow(i,j,k)  +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j+1,k)+d5(2)*rhow_n(i,j+2,k) &
                  + d5(1)*rhow_n(i,j-1,k)+d5(2)*rhow_n(i,j-2,k)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i,j+1,k)+d5(2)*rhoe_n(i,j+2,k) &
                  + d5(1)*rhoe_n(i,j-1,k)+d5(2)*rhoe_n(i,j-2,k)
          enddo
       enddo
       ! no filter for j=ny-1 & j=ny
    endif

    !****************
    if (is_2D) return
    !****************
    
    ! Filtering along k-direction
    ! ===========================
    
    ! BC at kmin
    ! ----------
    if (is_boundary(3,1)) then
       ! no filter for k=1 & k=2

       ! 4-th order centered filter
       k=3
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i,j,k+1)+d5(2)*rho_n(i,j,k+2) &
                  + d5(1)*rho_n(i,j,k-1)+d5(2)*rho_n(i,j,k-2)
             Krhou(i,j,k)= Krhou(i,j,k)  +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j,k+1)+d5(2)*rhou_n(i,j,k+2) &
                  + d5(1)*rhou_n(i,j,k-1)+d5(2)*rhou_n(i,j,k-2)
             Krhov(i,j,k)= Krhov(i,j,k)  +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j,k+1)+d5(2)*rhov_n(i,j,k+2) &
                  + d5(1)*rhov_n(i,j,k-1)+d5(2)*rhov_n(i,j,k-2)
             Krhow(i,j,k)= Krhow(i,j,k)  +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j,k+1)+d5(2)*rhow_n(i,j,k+2) &
                  + d5(1)*rhow_n(i,j,k-1)+d5(2)*rhow_n(i,j,k-2)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i,j,k+1)+d5(2)*rhoe_n(i,j,k+2) &
                  + d5(1)*rhoe_n(i,j,k-1)+d5(2)*rhoe_n(i,j,k-2)
          enddo
       enddo
       ! 6-th order centered filter
       k=4
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i,j,k+1)+d7(2)*rho_n(i,j,k+2)+d7(3)*rho_n(i,j,k+3) &
                  + d7(1)*rho_n(i,j,k-1)+d7(2)*rho_n(i,j,k-2)+d7(3)*rho_n(i,j,k-3)
             Krhou(i,j,k)= Krhou(i,j,k)  +d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i,j,k+1)+d7(2)*rhou_n(i,j,k+2)+d7(3)*rhou_n(i,j,k+3) &
                  + d7(1)*rhou_n(i,j,k-1)+d7(2)*rhou_n(i,j,k-2)+d7(3)*rhou_n(i,j,k-3)
             Krhov(i,j,k)= Krhov(i,j,k)  +d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i,j,k+1)+d7(2)*rhov_n(i,j,k+2)+d7(3)*rhov_n(i,j,k+3) &
                  + d7(1)*rhov_n(i,j,k-1)+d7(2)*rhov_n(i,j,k-2)+d7(3)*rhov_n(i,j,k-3)
             Krhow(i,j,k)= Krhow(i,j,k)  +d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i,j,k+1)+d7(2)*rhow_n(i,j,k+2)+d7(3)*rhow_n(i,j,k+3) &
                  + d7(1)*rhow_n(i,j,k-1)+d7(2)*rhow_n(i,j,k-2)+d7(3)*rhow_n(i,j,k-3)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i,j,k+1)+d7(2)*rhoe_n(i,j,k+2)+d7(3)*rhoe_n(i,j,k+3) &
                  + d7(1)*rhoe_n(i,j,k-1)+d7(2)*rhoe_n(i,j,k-2)+d7(3)*rhoe_n(i,j,k-3)
          enddo
       enddo
       ! 8-th order centered filter
       k=5
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d9(0)*rho_n(i,j,k)   &
                  + d9(1)*rho_n(i,j,k+1)+d9(2)*rho_n(i,j,k+2)+d9(3)*rho_n(i,j,k+3)+d9(4)*rho_n(i,j,k+4) &
                  + d9(1)*rho_n(i,j,k-1)+d9(2)*rho_n(i,j,k-2)+d9(3)*rho_n(i,j,k-3)+d9(4)*rho_n(i,j,k-4)
             Krhou(i,j,k)= Krhou(i,j,k)  +d9(0)*rhou_n(i,j,k)   &
                  + d9(1)*rhou_n(i,j,k+1)+d9(2)*rhou_n(i,j,k+2)+d9(3)*rhou_n(i,j,k+3)+d9(4)*rhou_n(i,j,k+4) &
                  + d9(1)*rhou_n(i,j,k-1)+d9(2)*rhou_n(i,j,k-2)+d9(3)*rhou_n(i,j,k-3)+d9(4)*rhou_n(i,j,k-4)
             Krhov(i,j,k)= Krhov(i,j,k)  +d9(0)*rhov_n(i,j,k)   &
                  + d9(1)*rhov_n(i,j,k+1)+d9(2)*rhov_n(i,j,k+2)+d9(3)*rhov_n(i,j,k+3)+d9(4)*rhov_n(i,j,k+4) &
                  + d9(1)*rhov_n(i,j,k-1)+d9(2)*rhov_n(i,j,k-2)+d9(3)*rhov_n(i,j,k-3)+d9(4)*rhov_n(i,j,k-4)
             Krhow(i,j,k)= Krhow(i,j,k)  +d9(0)*rhow_n(i,j,k)   &
                  + d9(1)*rhow_n(i,j,k+1)+d9(2)*rhow_n(i,j,k+2)+d9(3)*rhow_n(i,j,k+3)+d9(4)*rhow_n(i,j,k+4) &
                  + d9(1)*rhow_n(i,j,k-1)+d9(2)*rhow_n(i,j,k-2)+d9(3)*rhow_n(i,j,k-3)+d9(4)*rhow_n(i,j,k-4)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d9(0)*rhoe_n(i,j,k)   &
                  + d9(1)*rhoe_n(i,j,k+1)+d9(2)*rhoe_n(i,j,k+2)+d9(3)*rhoe_n(i,j,k+3)+d9(4)*rhoe_n(i,j,k+4) &
                  + d9(1)*rhoe_n(i,j,k-1)+d9(2)*rhoe_n(i,j,k-2)+d9(3)*rhoe_n(i,j,k-3)+d9(4)*rhoe_n(i,j,k-4)
          enddo
       enddo
    endif

    ! Interior points
    ! ---------------
    ! 10-th order centered filter
    do k=ndz,nfz
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  + d11(0)*rho_n(i,j,k) &
                  + d11(1)*( rho_n(i,j,k+1)+rho_n(i,j,k-1) )&
                  + d11(2)*( rho_n(i,j,k+2)+rho_n(i,j,k-2) )&
                  + d11(3)*( rho_n(i,j,k+3)+rho_n(i,j,k-3) )&
                  + d11(4)*( rho_n(i,j,k+4)+rho_n(i,j,k-4) )&
                  + d11(5)*( rho_n(i,j,k+5)+rho_n(i,j,k-5) )

             Krhou(i,j,k)= Krhou(i,j,k)  + d11(0)*rhou_n(i,j,k) &
                  + d11(1)*( rhou_n(i,j,k+1)+rhou_n(i,j,k-1) )&
                  + d11(2)*( rhou_n(i,j,k+2)+rhou_n(i,j,k-2) )&
                  + d11(3)*( rhou_n(i,j,k+3)+rhou_n(i,j,k-3) )&
                  + d11(4)*( rhou_n(i,j,k+4)+rhou_n(i,j,k-4) )&
                  + d11(5)*( rhou_n(i,j,k+5)+rhou_n(i,j,k-5) )

             Krhov(i,j,k)= Krhov(i,j,k)  + d11(0)*rhov_n(i,j,k) &
                  + d11(1)*( rhov_n(i,j,k+1)+rhov_n(i,j,k-1) )&
                  + d11(2)*( rhov_n(i,j,k+2)+rhov_n(i,j,k-2) )&
                  + d11(3)*( rhov_n(i,j,k+3)+rhov_n(i,j,k-3) )&
                  + d11(4)*( rhov_n(i,j,k+4)+rhov_n(i,j,k-4) )&
                  + d11(5)*( rhov_n(i,j,k+5)+rhov_n(i,j,k-5) )

             Krhow(i,j,k)= Krhow(i,j,k)  + d11(0)*rhow_n(i,j,k) &
                  + d11(1)*( rhow_n(i,j,k+1)+rhow_n(i,j,k-1) )&
                  + d11(2)*( rhow_n(i,j,k+2)+rhow_n(i,j,k-2) )&
                  + d11(3)*( rhow_n(i,j,k+3)+rhow_n(i,j,k-3) )&
                  + d11(4)*( rhow_n(i,j,k+4)+rhow_n(i,j,k-4) )&
                  + d11(5)*( rhow_n(i,j,k+5)+rhow_n(i,j,k-5) )

             Krhoe(i,j,k)= Krhoe(i,j,k)  + d11(0)*rhoe_n(i,j,k) &
                  + d11(1)*( rhoe_n(i,j,k+1)+rhoe_n(i,j,k-1) )&
                  + d11(2)*( rhoe_n(i,j,k+2)+rhoe_n(i,j,k-2) )&
                  + d11(3)*( rhoe_n(i,j,k+3)+rhoe_n(i,j,k-3) )&
                  + d11(4)*( rhoe_n(i,j,k+4)+rhoe_n(i,j,k-4) )&
                  + d11(5)*( rhoe_n(i,j,k+5)+rhoe_n(i,j,k-5) )
          enddo
       enddo
    enddo

    ! BC at kmax
    ! ----------
    if (is_boundary(3,2)) then
       ! 8-th order centered filter
       k=nz-4
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d9(0)*rho_n(i,j,k)   &
                  + d9(1)*rho_n(i,j,k+1)+d9(2)*rho_n(i,j,k+2)+d9(3)*rho_n(i,j,k+3)+d9(4)*rho_n(i,j,k+4) &
                  + d9(1)*rho_n(i,j,k-1)+d9(2)*rho_n(i,j,k-2)+d9(3)*rho_n(i,j,k-3)+d9(4)*rho_n(i,j,k-4)
             Krhou(i,j,k)= Krhou(i,j,k)  +d9(0)*rhou_n(i,j,k)   &
                  + d9(1)*rhou_n(i,j,k+1)+d9(2)*rhou_n(i,j,k+2)+d9(3)*rhou_n(i,j,k+3)+d9(4)*rhou_n(i,j,k+4) &
                  + d9(1)*rhou_n(i,j,k-1)+d9(2)*rhou_n(i,j,k-2)+d9(3)*rhou_n(i,j,k-3)+d9(4)*rhou_n(i,j,k-4)
             Krhov(i,j,k)= Krhov(i,j,k)  +d9(0)*rhov_n(i,j,k)   &
                  + d9(1)*rhov_n(i,j,k+1)+d9(2)*rhov_n(i,j,k+2)+d9(3)*rhov_n(i,j,k+3)+d9(4)*rhov_n(i,j,k+4) &
                  + d9(1)*rhov_n(i,j,k-1)+d9(2)*rhov_n(i,j,k-2)+d9(3)*rhov_n(i,j,k-3)+d9(4)*rhov_n(i,j,k-4)
             Krhow(i,j,k)= Krhow(i,j,k)  +d9(0)*rhow_n(i,j,k)   &
                  + d9(1)*rhow_n(i,j,k+1)+d9(2)*rhow_n(i,j,k+2)+d9(3)*rhow_n(i,j,k+3)+d9(4)*rhow_n(i,j,k+4) &
                  + d9(1)*rhow_n(i,j,k-1)+d9(2)*rhow_n(i,j,k-2)+d9(3)*rhow_n(i,j,k-3)+d9(4)*rhow_n(i,j,k-4)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d9(0)*rhoe_n(i,j,k)   &
                  + d9(1)*rhoe_n(i,j,k+1)+d9(2)*rhoe_n(i,j,k+2)+d9(3)*rhoe_n(i,j,k+3)+d9(4)*rhoe_n(i,j,k+4) &
                  + d9(1)*rhoe_n(i,j,k-1)+d9(2)*rhoe_n(i,j,k-2)+d9(3)*rhoe_n(i,j,k-3)+d9(4)*rhoe_n(i,j,k-4)
          enddo
       enddo
       ! 6-th order centered filter
       k=nz-3
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i,j,k+1)+d7(2)*rho_n(i,j,k+2)+d7(3)*rho_n(i,j,k+3) &
                  + d7(1)*rho_n(i,j,k-1)+d7(2)*rho_n(i,j,k-2)+d7(3)*rho_n(i,j,k-3)
             Krhou(i,j,k)= Krhou(i,j,k)  +d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i,j,k+1)+d7(2)*rhou_n(i,j,k+2)+d7(3)*rhou_n(i,j,k+3) &
                  + d7(1)*rhou_n(i,j,k-1)+d7(2)*rhou_n(i,j,k-2)+d7(3)*rhou_n(i,j,k-3)
             Krhov(i,j,k)= Krhov(i,j,k)  +d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i,j,k+1)+d7(2)*rhov_n(i,j,k+2)+d7(3)*rhov_n(i,j,k+3) &
                  + d7(1)*rhov_n(i,j,k-1)+d7(2)*rhov_n(i,j,k-2)+d7(3)*rhov_n(i,j,k-3)
             Krhow(i,j,k)= Krhow(i,j,k)  +d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i,j,k+1)+d7(2)*rhow_n(i,j,k+2)+d7(3)*rhow_n(i,j,k+3) &
                  + d7(1)*rhow_n(i,j,k-1)+d7(2)*rhow_n(i,j,k-2)+d7(3)*rhow_n(i,j,k-3)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i,j,k+1)+d7(2)*rhoe_n(i,j,k+2)+d7(3)*rhoe_n(i,j,k+3) &
                  + d7(1)*rhoe_n(i,j,k-1)+d7(2)*rhoe_n(i,j,k-2)+d7(3)*rhoe_n(i,j,k-3)
          enddo
       enddo
       ! 4-th order centered filter
       k=nz-2
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i,j,k+1)+d5(2)*rho_n(i,j,k+2) &
                  + d5(1)*rho_n(i,j,k-1)+d5(2)*rho_n(i,j,k-2)
             Krhou(i,j,k)= Krhou(i,j,k)  +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j,k+1)+d5(2)*rhou_n(i,j,k+2) &
                  + d5(1)*rhou_n(i,j,k-1)+d5(2)*rhou_n(i,j,k-2)
             Krhov(i,j,k)= Krhov(i,j,k)  +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j,k+1)+d5(2)*rhov_n(i,j,k+2) &
                  + d5(1)*rhov_n(i,j,k-1)+d5(2)*rhov_n(i,j,k-2)
             Krhow(i,j,k)= Krhow(i,j,k)  +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j,k+1)+d5(2)*rhow_n(i,j,k+2) &
                  + d5(1)*rhow_n(i,j,k-1)+d5(2)*rhow_n(i,j,k-2)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i,j,k+1)+d5(2)*rhoe_n(i,j,k+2) &
                  + d5(1)*rhoe_n(i,j,k-1)+d5(2)*rhoe_n(i,j,k-2)
          enddo
       enddo
       ! no filter for k=nz-1 & k=nz
    endif
        
  end subroutine filtering_11pts

  !==============================================================================
  subroutine filtering_9pts
  !==============================================================================
    !> Apply filtering on 9-point stencils (+ boundaries)
    !> * old version without comput./comm. overlap *
  !==============================================================================
    use mod_flow
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    ! ---------------------------------------------------------------------------

    ! Filtering along i-direction
    ! ===========================
    
    ! BC at imin
    ! ----------
    if (is_boundary(1,1)) then
       ! no filter for i=1 & i=2 -> initialize array
       do j=1,ny
           Krho(1,j,:)=0.0_wp
          Krhou(1,j,:)=0.0_wp
          Krhov(1,j,:)=0.0_wp
          Krhow(1,j,:)=0.0_wp
          Krhoe(1,j,:)=0.0_wp
          
           Krho(2,j,:)=0.0_wp
          Krhou(2,j,:)=0.0_wp
          Krhov(2,j,:)=0.0_wp
          Krhow(2,j,:)=0.0_wp
          Krhoe(2,j,:)=0.0_wp
       enddo
       ! 4-th order centered filter
       i=3
       do k=1,nz
          do j=1,ny
             Krho(i,j,k)  = d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i+1,j,k)+d5(2)*rho_n(i+2,j,k) &
                  + d5(1)*rho_n(i-1,j,k)+d5(2)*rho_n(i-2,j,k)
             Krhou(i,j,k) = d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i+1,j,k)+d5(2)*rhou_n(i+2,j,k) &
                  + d5(1)*rhou_n(i-1,j,k)+d5(2)*rhou_n(i-2,j,k)
             Krhov(i,j,k) = d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i+1,j,k)+d5(2)*rhov_n(i+2,j,k) &
                  + d5(1)*rhov_n(i-1,j,k)+d5(2)*rhov_n(i-2,j,k)
             Krhow(i,j,k) = d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i+1,j,k)+d5(2)*rhow_n(i+2,j,k) &
                  + d5(1)*rhow_n(i-1,j,k)+d5(2)*rhow_n(i-2,j,k)
             Krhoe(i,j,k) = d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i+1,j,k)+d5(2)*rhoe_n(i+2,j,k) &
                  + d5(1)*rhoe_n(i-1,j,k)+d5(2)*rhoe_n(i-2,j,k)
          enddo
       enddo
       ! 6-th order centered filter
       i=4
       do k=1,nz
          do j=1,ny
             Krho(i,j,k)  = d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i+1,j,k)+d7(2)*rho_n(i+2,j,k)+d7(3)*rho_n(i+3,j,k) &
                  + d7(1)*rho_n(i-1,j,k)+d7(2)*rho_n(i-2,j,k)+d7(3)*rho_n(i-3,j,k)
             Krhou(i,j,k) = d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i+1,j,k)+d7(2)*rhou_n(i+2,j,k)+d7(3)*rhou_n(i+3,j,k) &
                  + d7(1)*rhou_n(i-1,j,k)+d7(2)*rhou_n(i-2,j,k)+d7(3)*rhou_n(i-3,j,k)
             Krhov(i,j,k) = d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i+1,j,k)+d7(2)*rhov_n(i+2,j,k)+d7(3)*rhov_n(i+3,j,k) &
                  + d7(1)*rhov_n(i-1,j,k)+d7(2)*rhov_n(i-2,j,k)+d7(3)*rhov_n(i-3,j,k)
             Krhow(i,j,k) = d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i+1,j,k)+d7(2)*rhow_n(i+2,j,k)+d7(3)*rhow_n(i+3,j,k) &
                  + d7(1)*rhow_n(i-1,j,k)+d7(2)*rhow_n(i-2,j,k)+d7(3)*rhow_n(i-3,j,k)
             Krhoe(i,j,k) = d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i+1,j,k)+d7(2)*rhoe_n(i+2,j,k)+d7(3)*rhoe_n(i+3,j,k) &
                  + d7(1)*rhoe_n(i-1,j,k)+d7(2)*rhoe_n(i-2,j,k)+d7(3)*rhoe_n(i-3,j,k)
          enddo
       enddo
    endif

    ! Interior points
    ! ---------------
    ! 8-th order centered filter
    do k=1,nz
       do j=1,ny
          do i=ndx,nfx
             Krho(i,j,k)  = d9(0)*rho_n(i,j,k) &
                  + d9(1)*( rho_n(i+1,j,k)+rho_n(i-1,j,k) )&
                  + d9(2)*( rho_n(i+2,j,k)+rho_n(i-2,j,k) )&
                  + d9(3)*( rho_n(i+3,j,k)+rho_n(i-3,j,k) )&
                  + d9(4)*( rho_n(i+4,j,k)+rho_n(i-4,j,k) )

             Krhou(i,j,k) = d9(0)*rhou_n(i,j,k) &
                  + d9(1)*( rhou_n(i+1,j,k)+rhou_n(i-1,j,k) )&
                  + d9(2)*( rhou_n(i+2,j,k)+rhou_n(i-2,j,k) )&
                  + d9(3)*( rhou_n(i+3,j,k)+rhou_n(i-3,j,k) )&
                  + d9(4)*( rhou_n(i+4,j,k)+rhou_n(i-4,j,k) )

             Krhov(i,j,k) = d9(0)*rhov_n(i,j,k) &
                  + d9(1)*( rhov_n(i+1,j,k)+rhov_n(i-1,j,k) )&
                  + d9(2)*( rhov_n(i+2,j,k)+rhov_n(i-2,j,k) )&
                  + d9(3)*( rhov_n(i+3,j,k)+rhov_n(i-3,j,k) )&
                  + d9(4)*( rhov_n(i+4,j,k)+rhov_n(i-4,j,k) )

             Krhow(i,j,k) = d9(0)*rhow_n(i,j,k) &
                  + d9(1)*( rhow_n(i+1,j,k)+rhow_n(i-1,j,k) )&
                  + d9(2)*( rhow_n(i+2,j,k)+rhow_n(i-2,j,k) )&
                  + d9(3)*( rhow_n(i+3,j,k)+rhow_n(i-3,j,k) )&
                  + d9(4)*( rhow_n(i+4,j,k)+rhow_n(i-4,j,k) )

             Krhoe(i,j,k) = d9(0)*rhoe_n(i,j,k) &
                  + d9(1)*( rhoe_n(i+1,j,k)+rhoe_n(i-1,j,k) )&
                  + d9(2)*( rhoe_n(i+2,j,k)+rhoe_n(i-2,j,k) )&
                  + d9(3)*( rhoe_n(i+3,j,k)+rhoe_n(i-3,j,k) )&
                  + d9(4)*( rhoe_n(i+4,j,k)+rhoe_n(i-4,j,k) )
          enddo
       enddo
    enddo

    ! BC at imax
    ! ----------
    if (is_boundary(1,2)) then
       ! 6-th order centered filter
       i=nx-3
       do k=1,nz
          do j=1,ny
             Krho(i,j,k)  = d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i+1,j,k)+d7(2)*rho_n(i+2,j,k)+d7(3)*rho_n(i+3,j,k) &
                  + d7(1)*rho_n(i-1,j,k)+d7(2)*rho_n(i-2,j,k)+d7(3)*rho_n(i-3,j,k)
             Krhou(i,j,k) = d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i+1,j,k)+d7(2)*rhou_n(i+2,j,k)+d7(3)*rhou_n(i+3,j,k) &
                  + d7(1)*rhou_n(i-1,j,k)+d7(2)*rhou_n(i-2,j,k)+d7(3)*rhou_n(i-3,j,k)
             Krhov(i,j,k) = d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i+1,j,k)+d7(2)*rhov_n(i+2,j,k)+d7(3)*rhov_n(i+3,j,k) &
                  + d7(1)*rhov_n(i-1,j,k)+d7(2)*rhov_n(i-2,j,k)+d7(3)*rhov_n(i-3,j,k)
             Krhow(i,j,k) = d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i+1,j,k)+d7(2)*rhow_n(i+2,j,k)+d7(3)*rhow_n(i+3,j,k) &
                  + d7(1)*rhow_n(i-1,j,k)+d7(2)*rhow_n(i-2,j,k)+d7(3)*rhow_n(i-3,j,k)
             Krhoe(i,j,k) = d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i+1,j,k)+d7(2)*rhoe_n(i+2,j,k)+d7(3)*rhoe_n(i+3,j,k) &
                  + d7(1)*rhoe_n(i-1,j,k)+d7(2)*rhoe_n(i-2,j,k)+d7(3)*rhoe_n(i-3,j,k)
          enddo
       enddo
       ! 4-th order centered filter
       i=nx-2
       do k=1,nz
          do j=1,ny
             Krho(i,j,k)  = d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i+1,j,k)+d5(2)*rho_n(i+2,j,k) &
                  + d5(1)*rho_n(i-1,j,k)+d5(2)*rho_n(i-2,j,k)
             Krhou(i,j,k) = d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i+1,j,k)+d5(2)*rhou_n(i+2,j,k) &
                  + d5(1)*rhou_n(i-1,j,k)+d5(2)*rhou_n(i-2,j,k)
             Krhov(i,j,k) = d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i+1,j,k)+d5(2)*rhov_n(i+2,j,k) &
                  + d5(1)*rhov_n(i-1,j,k)+d5(2)*rhov_n(i-2,j,k)
             Krhow(i,j,k) = d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i+1,j,k)+d5(2)*rhow_n(i+2,j,k) &
                  + d5(1)*rhow_n(i-1,j,k)+d5(2)*rhow_n(i-2,j,k)
             Krhoe(i,j,k) = d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i+1,j,k)+d5(2)*rhoe_n(i+2,j,k) &
                  + d5(1)*rhoe_n(i-1,j,k)+d5(2)*rhoe_n(i-2,j,k)
          enddo
       enddo
       ! no filter for i=nx-1 & i=nx -> initialize array
       do j=1,ny
           Krho(nx-1,j,:)=0.0_wp
          Krhou(nx-1,j,:)=0.0_wp
          Krhov(nx-1,j,:)=0.0_wp
          Krhow(nx-1,j,:)=0.0_wp
          Krhoe(nx-1,j,:)=0.0_wp
          
           Krho(nx,j,:)=0.0_wp
          Krhou(nx,j,:)=0.0_wp
          Krhov(nx,j,:)=0.0_wp
          Krhow(nx,j,:)=0.0_wp
          Krhoe(nx,j,:)=0.0_wp
       enddo
    endif

    ! Filtering along j-direction
    ! ===========================
    
    ! BC at jmin
    ! ----------
    if (is_boundary(2,1)) then
       ! 4-th order centered filter
       j=3
       do k=1,nz
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i,j+1,k) +d5(2)*rho_n(i,j+2,k) &
                  + d5(1)*rho_n(i,j-1,k) +d5(2)*rho_n(i,j-2,k)
             Krhou(i,j,k) = Krhou(i,j,k) +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j+1,k)+d5(2)*rhou_n(i,j+2,k) &
                  + d5(1)*rhou_n(i,j-1,k)+d5(2)*rhou_n(i,j-2,k)
             Krhov(i,j,k) = Krhov(i,j,k) +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j+1,k)+d5(2)*rhov_n(i,j+2,k) &
                  + d5(1)*rhov_n(i,j-1,k)+d5(2)*rhov_n(i,j-2,k)
             Krhow(i,j,k) = Krhow(i,j,k) +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j+1,k)+d5(2)*rhow_n(i,j+2,k) &
                  + d5(1)*rhow_n(i,j-1,k)+d5(2)*rhow_n(i,j-2,k)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i,j+1,k)+d5(2)*rhoe_n(i,j+2,k) &
                  + d5(1)*rhoe_n(i,j-1,k)+d5(2)*rhoe_n(i,j-2,k)
          enddo
       enddo
       ! 6-th order centered filter
       j=4
       do k=1,nz
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  +d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i,j+1,k)+d7(2)*rho_n(i,j+2,k)+d7(3)*rho_n(i,j+3,k) &
                  + d7(1)*rho_n(i,j-1,k)+d7(2)*rho_n(i,j-2,k)+d7(3)*rho_n(i,j-3,k)
             Krhou(i,j,k) = Krhou(i,j,k)  +d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i,j+1,k)+d7(2)*rhou_n(i,j+2,k)+d7(3)*rhou_n(i,j+3,k) &
                  + d7(1)*rhou_n(i,j-1,k)+d7(2)*rhou_n(i,j-2,k)+d7(3)*rhou_n(i,j-3,k)
             Krhov(i,j,k) = Krhov(i,j,k)  +d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i,j+1,k)+d7(2)*rhov_n(i,j+2,k)+d7(3)*rhov_n(i,j+3,k) &
                  + d7(1)*rhov_n(i,j-1,k)+d7(2)*rhov_n(i,j-2,k)+d7(3)*rhov_n(i,j-3,k)
             Krhow(i,j,k) = Krhow(i,j,k)  +d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i,j+1,k)+d7(2)*rhow_n(i,j+2,k)+d7(3)*rhow_n(i,j+3,k) &
                  + d7(1)*rhow_n(i,j-1,k)+d7(2)*rhow_n(i,j-2,k)+d7(3)*rhow_n(i,j-3,k)
             Krhoe(i,j,k) = Krhoe(i,j,k)  +d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i,j+1,k)+d7(2)*rhoe_n(i,j+2,k)+d7(3)*rhoe_n(i,j+3,k) &
                  + d7(1)*rhoe_n(i,j-1,k)+d7(2)*rhoe_n(i,j-2,k)+d7(3)*rhoe_n(i,j-3,k)
          enddo
       enddo
    endif

    ! Interior points
    ! ---------------
    ! 8-th order centered filter
    do k=1,nz
       do j=ndy,nfy
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  + d9(0)*rho_n(i,j,k) &
                  + d9(1)*( rho_n(i,j+1,k)+rho_n(i,j-1,k) )&
                  + d9(2)*( rho_n(i,j+2,k)+rho_n(i,j-2,k) )&
                  + d9(3)*( rho_n(i,j+3,k)+rho_n(i,j-3,k) )&
                  + d9(4)*( rho_n(i,j+4,k)+rho_n(i,j-4,k) )

             Krhou(i,j,k) = Krhou(i,j,k)  + d9(0)*rhou_n(i,j,k) &
                  + d9(1)*( rhou_n(i,j+1,k)+rhou_n(i,j-1,k) )&
                  + d9(2)*( rhou_n(i,j+2,k)+rhou_n(i,j-2,k) )&
                  + d9(3)*( rhou_n(i,j+3,k)+rhou_n(i,j-3,k) )&
                  + d9(4)*( rhou_n(i,j+4,k)+rhou_n(i,j-4,k) )

             Krhov(i,j,k) = Krhov(i,j,k)  + d9(0)*rhov_n(i,j,k) &
                  + d9(1)*( rhov_n(i,j+1,k)+rhov_n(i,j-1,k) )&
                  + d9(2)*( rhov_n(i,j+2,k)+rhov_n(i,j-2,k) )&
                  + d9(3)*( rhov_n(i,j+3,k)+rhov_n(i,j-3,k) )&
                  + d9(4)*( rhov_n(i,j+4,k)+rhov_n(i,j-4,k) )

             Krhow(i,j,k) = Krhow(i,j,k)  + d9(0)*rhow_n(i,j,k) &
                  + d9(1)*( rhow_n(i,j+1,k)+rhow_n(i,j-1,k) )&
                  + d9(2)*( rhow_n(i,j+2,k)+rhow_n(i,j-2,k) )&
                  + d9(3)*( rhow_n(i,j+3,k)+rhow_n(i,j-3,k) )&
                  + d9(4)*( rhow_n(i,j+4,k)+rhow_n(i,j-4,k) )

             Krhoe(i,j,k) = Krhoe(i,j,k)  + d9(0)*rhoe_n(i,j,k) &
                  + d9(1)*( rhoe_n(i,j+1,k)+rhoe_n(i,j-1,k) )&
                  + d9(2)*( rhoe_n(i,j+2,k)+rhoe_n(i,j-2,k) )&
                  + d9(3)*( rhoe_n(i,j+3,k)+rhoe_n(i,j-3,k) )&
                  + d9(4)*( rhoe_n(i,j+4,k)+rhoe_n(i,j-4,k) )
          enddo
       enddo
    enddo

    ! BC at jmax
    ! ----------
    if (is_boundary(2,2)) then
       ! 6-th order centered filter
       j=ny-3
       do k=1,nz
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  +d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i,j+1,k)+d7(2)*rho_n(i,j+2,k)+d7(3)*rho_n(i,j+3,k) &
                  + d7(1)*rho_n(i,j-1,k)+d7(2)*rho_n(i,j-2,k)+d7(3)*rho_n(i,j-3,k)
             Krhou(i,j,k) = Krhou(i,j,k)  +d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i,j+1,k)+d7(2)*rhou_n(i,j+2,k)+d7(3)*rhou_n(i,j+3,k) &
                  + d7(1)*rhou_n(i,j-1,k)+d7(2)*rhou_n(i,j-2,k)+d7(3)*rhou_n(i,j-3,k)
             Krhov(i,j,k) = Krhov(i,j,k)  +d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i,j+1,k)+d7(2)*rhov_n(i,j+2,k)+d7(3)*rhov_n(i,j+3,k) &
                  + d7(1)*rhov_n(i,j-1,k)+d7(2)*rhov_n(i,j-2,k)+d7(3)*rhov_n(i,j-3,k)
             Krhow(i,j,k) = Krhow(i,j,k)  +d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i,j+1,k)+d7(2)*rhow_n(i,j+2,k)+d7(3)*rhow_n(i,j+3,k) &
                  + d7(1)*rhow_n(i,j-1,k)+d7(2)*rhow_n(i,j-2,k)+d7(3)*rhow_n(i,j-3,k)
             Krhoe(i,j,k) = Krhoe(i,j,k)  +d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i,j+1,k)+d7(2)*rhoe_n(i,j+2,k)+d7(3)*rhoe_n(i,j+3,k) &
                  + d7(1)*rhoe_n(i,j-1,k)+d7(2)*rhoe_n(i,j-2,k)+d7(3)*rhoe_n(i,j-3,k)
          enddo
       enddo
       ! 4-th order centered filter
       j=ny-2
       do k=1,nz
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i,j+1,k)+d5(2)*rho_n(i,j+2,k) &
                  + d5(1)*rho_n(i,j-1,k)+d5(2)*rho_n(i,j-2,k)
             Krhou(i,j,k) = Krhou(i,j,k)  +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j+1,k)+d5(2)*rhou_n(i,j+2,k) &
                  + d5(1)*rhou_n(i,j-1,k)+d5(2)*rhou_n(i,j-2,k)
             Krhov(i,j,k) = Krhov(i,j,k)  +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j+1,k)+d5(2)*rhov_n(i,j+2,k) &
                  + d5(1)*rhov_n(i,j-1,k)+d5(2)*rhov_n(i,j-2,k)
             Krhow(i,j,k) = Krhow(i,j,k)  +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j+1,k)+d5(2)*rhow_n(i,j+2,k) &
                  + d5(1)*rhow_n(i,j-1,k)+d5(2)*rhow_n(i,j-2,k)
             Krhoe(i,j,k) = Krhoe(i,j,k)  +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i,j+1,k)+d5(2)*rhoe_n(i,j+2,k) &
                  + d5(1)*rhoe_n(i,j-1,k)+d5(2)*rhoe_n(i,j-2,k)
          enddo
       enddo
    endif

    !****************
    if (is_2D) return
    !****************

    ! Filtering along k-direction
    ! ===========================
    
    ! BC at kmin
    ! ----------
    if (is_boundary(3,1)) then
       ! 4-th order centered filter
       k=3
       do j=1,ny
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i,j,k+1)+d5(2)*rho_n(i,j,k+2) &
                  + d5(1)*rho_n(i,j,k-1)+d5(2)*rho_n(i,j,k-2)
             Krhou(i,j,k) = Krhou(i,j,k)  +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j,k+1)+d5(2)*rhou_n(i,j,k+2) &
                  + d5(1)*rhou_n(i,j,k-1)+d5(2)*rhou_n(i,j,k-2)
             Krhov(i,j,k) = Krhov(i,j,k)  +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j,k+1)+d5(2)*rhov_n(i,j,k+2) &
                  + d5(1)*rhov_n(i,j,k-1)+d5(2)*rhov_n(i,j,k-2)
             Krhow(i,j,k) = Krhow(i,j,k)  +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j,k+1)+d5(2)*rhow_n(i,j,k+2) &
                  + d5(1)*rhow_n(i,j,k-1)+d5(2)*rhow_n(i,j,k-2)
             Krhoe(i,j,k) = Krhoe(i,j,k)  +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i,j,k+1)+d5(2)*rhoe_n(i,j,k+2) &
                  + d5(1)*rhoe_n(i,j,k-1)+d5(2)*rhoe_n(i,j,k-2)
          enddo
       enddo
       ! 6-th order centered filter
       k=4
       do j=1,ny
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  +d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i,j,k+1)+d7(2)*rho_n(i,j,k+2)+d7(3)*rho_n(i,j,k+3) &
                  + d7(1)*rho_n(i,j,k-1)+d7(2)*rho_n(i,j,k-2)+d7(3)*rho_n(i,j,k-3)
             Krhou(i,j,k) = Krhou(i,j,k)  +d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i,j,k+1)+d7(2)*rhou_n(i,j,k+2)+d7(3)*rhou_n(i,j,k+3) &
                  + d7(1)*rhou_n(i,j,k-1)+d7(2)*rhou_n(i,j,k-2)+d7(3)*rhou_n(i,j,k-3)
             Krhov(i,j,k) = Krhov(i,j,k)  +d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i,j,k+1)+d7(2)*rhov_n(i,j,k+2)+d7(3)*rhov_n(i,j,k+3) &
                  + d7(1)*rhov_n(i,j,k-1)+d7(2)*rhov_n(i,j,k-2)+d7(3)*rhov_n(i,j,k-3)
             Krhow(i,j,k) = Krhow(i,j,k)  +d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i,j,k+1)+d7(2)*rhow_n(i,j,k+2)+d7(3)*rhow_n(i,j,k+3) &
                  + d7(1)*rhow_n(i,j,k-1)+d7(2)*rhow_n(i,j,k-2)+d7(3)*rhow_n(i,j,k-3)
             Krhoe(i,j,k) = Krhoe(i,j,k)  +d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i,j,k+1)+d7(2)*rhoe_n(i,j,k+2)+d7(3)*rhoe_n(i,j,k+3) &
                  + d7(1)*rhoe_n(i,j,k-1)+d7(2)*rhoe_n(i,j,k-2)+d7(3)*rhoe_n(i,j,k-3)
          enddo
       enddo
    endif

    ! Interior points
    ! ---------------
    ! 8-th order centered filter
    do k=ndz,nfz
       do j=1,ny
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  + d9(0)*rho_n(i,j,k) &
                  + d9(1)*( rho_n(i,j,k+1)+rho_n(i,j,k-1) )&
                  + d9(2)*( rho_n(i,j,k+2)+rho_n(i,j,k-2) )&
                  + d9(3)*( rho_n(i,j,k+3)+rho_n(i,j,k-3) )&
                  + d9(4)*( rho_n(i,j,k+4)+rho_n(i,j,k-4) )

             Krhou(i,j,k) = Krhou(i,j,k)  + d9(0)*rhou_n(i,j,k) &
                  + d9(1)*( rhou_n(i,j,k+1)+rhou_n(i,j,k-1) )&
                  + d9(2)*( rhou_n(i,j,k+2)+rhou_n(i,j,k-2) )&
                  + d9(3)*( rhou_n(i,j,k+3)+rhou_n(i,j,k-3) )&
                  + d9(4)*( rhou_n(i,j,k+4)+rhou_n(i,j,k-4) )

             Krhov(i,j,k) = Krhov(i,j,k)  + d9(0)*rhov_n(i,j,k) &
                  + d9(1)*( rhov_n(i,j,k+1)+rhov_n(i,j,k-1) )&
                  + d9(2)*( rhov_n(i,j,k+2)+rhov_n(i,j,k-2) )&
                  + d9(3)*( rhov_n(i,j,k+3)+rhov_n(i,j,k-3) )&
                  + d9(4)*( rhov_n(i,j,k+4)+rhov_n(i,j,k-4) )

             Krhow(i,j,k) = Krhow(i,j,k)  + d9(0)*rhow_n(i,j,k) &
                  + d9(1)*( rhow_n(i,j,k+1)+rhow_n(i,j,k-1) )&
                  + d9(2)*( rhow_n(i,j,k+2)+rhow_n(i,j,k-2) )&
                  + d9(3)*( rhow_n(i,j,k+3)+rhow_n(i,j,k-3) )&
                  + d9(4)*( rhow_n(i,j,k+4)+rhow_n(i,j,k-4) )

             Krhoe(i,j,k) = Krhoe(i,j,k)  + d9(0)*rhoe_n(i,j,k) &
                  + d9(1)*( rhoe_n(i,j,k+1)+rhoe_n(i,j,k-1) )&
                  + d9(2)*( rhoe_n(i,j,k+2)+rhoe_n(i,j,k-2) )&
                  + d9(3)*( rhoe_n(i,j,k+3)+rhoe_n(i,j,k-3) )&
                  + d9(4)*( rhoe_n(i,j,k+4)+rhoe_n(i,j,k-4) )
          enddo
       enddo
    enddo

    ! BC at kmax
    ! ----------
    if (is_boundary(3,2)) then
       ! 6-th order centered filter
       k=nz-3
       do j=1,ny
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  +d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i,j,k+1)+d7(2)*rho_n(i,j,k+2)+d7(3)*rho_n(i,j,k+3) &
                  + d7(1)*rho_n(i,j,k-1)+d7(2)*rho_n(i,j,k-2)+d7(3)*rho_n(i,j,k-3)
             Krhou(i,j,k) = Krhou(i,j,k)  +d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i,j,k+1)+d7(2)*rhou_n(i,j,k+2)+d7(3)*rhou_n(i,j,k+3) &
                  + d7(1)*rhou_n(i,j,k-1)+d7(2)*rhou_n(i,j,k-2)+d7(3)*rhou_n(i,j,k-3)
             Krhov(i,j,k) = Krhov(i,j,k)  +d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i,j,k+1)+d7(2)*rhov_n(i,j,k+2)+d7(3)*rhov_n(i,j,k+3) &
                  + d7(1)*rhov_n(i,j,k-1)+d7(2)*rhov_n(i,j,k-2)+d7(3)*rhov_n(i,j,k-3)
             Krhow(i,j,k) = Krhow(i,j,k)  +d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i,j,k+1)+d7(2)*rhow_n(i,j,k+2)+d7(3)*rhow_n(i,j,k+3) &
                  + d7(1)*rhow_n(i,j,k-1)+d7(2)*rhow_n(i,j,k-2)+d7(3)*rhow_n(i,j,k-3)
             Krhoe(i,j,k) = Krhoe(i,j,k)  +d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i,j,k+1)+d7(2)*rhoe_n(i,j,k+2)+d7(3)*rhoe_n(i,j,k+3) &
                  + d7(1)*rhoe_n(i,j,k-1)+d7(2)*rhoe_n(i,j,k-2)+d7(3)*rhoe_n(i,j,k-3)
          enddo
       enddo
       ! 4-th order centered filter
       k=nz-2
       do j=1,ny
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i,j,k+1)+d5(2)*rho_n(i,j,k+2) &
                  + d5(1)*rho_n(i,j,k-1)+d5(2)*rho_n(i,j,k-2)
             Krhou(i,j,k) = Krhou(i,j,k)  +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j,k+1)+d5(2)*rhou_n(i,j,k+2) &
                  + d5(1)*rhou_n(i,j,k-1)+d5(2)*rhou_n(i,j,k-2)
             Krhov(i,j,k) = Krhov(i,j,k)  +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j,k+1)+d5(2)*rhov_n(i,j,k+2) &
                  + d5(1)*rhov_n(i,j,k-1)+d5(2)*rhov_n(i,j,k-2)
             Krhow(i,j,k) = Krhow(i,j,k)  +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j,k+1)+d5(2)*rhow_n(i,j,k+2) &
                  + d5(1)*rhow_n(i,j,k-1)+d5(2)*rhow_n(i,j,k-2)
             Krhoe(i,j,k) = Krhoe(i,j,k)  +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i,j,k+1)+d5(2)*rhoe_n(i,j,k+2) &
                  + d5(1)*rhoe_n(i,j,k-1)+d5(2)*rhoe_n(i,j,k-2)
          enddo
       enddo
    endif

  end subroutine filtering_9pts

  !==============================================================================
  subroutine filtering_7pts
  !==============================================================================
    !> Apply filtering on 7-point stencils (+ boundaries)
    !> * old version without comput./comm. overlap *
  !==============================================================================
    use mod_flow
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    ! ---------------------------------------------------------------------------

    ! Filtering along i-direction
    ! ===========================
    
    ! BC at imin
    ! ----------
    if (is_boundary(1,1)) then
       ! no filter for i=1 & i=2 -> initialize array
       do j=1,ny
           Krho(1,j,:)=0.0_wp
          Krhou(1,j,:)=0.0_wp
          Krhov(1,j,:)=0.0_wp
          Krhow(1,j,:)=0.0_wp
          Krhoe(1,j,:)=0.0_wp
          
           Krho(2,j,:)=0.0_wp
          Krhou(2,j,:)=0.0_wp
          Krhov(2,j,:)=0.0_wp
          Krhow(2,j,:)=0.0_wp
          Krhoe(2,j,:)=0.0_wp
       enddo
       ! 4-th order centered filter
       i=3
       do k=1,nz
          do j=1,ny
             Krho(i,j,k)  = d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i+1,j,k)+d5(2)*rho_n(i+2,j,k) &
                  + d5(1)*rho_n(i-1,j,k)+d5(2)*rho_n(i-2,j,k)
             Krhou(i,j,k) = d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i+1,j,k)+d5(2)*rhou_n(i+2,j,k) &
                  + d5(1)*rhou_n(i-1,j,k)+d5(2)*rhou_n(i-2,j,k)
             Krhov(i,j,k) = d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i+1,j,k)+d5(2)*rhov_n(i+2,j,k) &
                  + d5(1)*rhov_n(i-1,j,k)+d5(2)*rhov_n(i-2,j,k)
             Krhow(i,j,k) = d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i+1,j,k)+d5(2)*rhow_n(i+2,j,k) &
                  + d5(1)*rhow_n(i-1,j,k)+d5(2)*rhow_n(i-2,j,k)
             Krhoe(i,j,k) = d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i+1,j,k)+d5(2)*rhoe_n(i+2,j,k) &
                  + d5(1)*rhoe_n(i-1,j,k)+d5(2)*rhoe_n(i-2,j,k)
          enddo
       enddo
    endif

    ! Interior points
    ! ---------------
    ! 6-th order centered filter
    do k=1,nz
       do j=1,ny
          do i=ndx,nfx
             Krho(i,j,k)  = d7(0)*rho_n(i,j,k) &
                  + d7(1)*( rho_n(i+1,j,k)+rho_n(i-1,j,k) )&
                  + d7(2)*( rho_n(i+2,j,k)+rho_n(i-2,j,k) )&
                  + d7(3)*( rho_n(i+3,j,k)+rho_n(i-3,j,k) )

             Krhou(i,j,k) = d7(0)*rhou_n(i,j,k) &
                  + d7(1)*( rhou_n(i+1,j,k)+rhou_n(i-1,j,k) )&
                  + d7(2)*( rhou_n(i+2,j,k)+rhou_n(i-2,j,k) )&
                  + d7(3)*( rhou_n(i+3,j,k)+rhou_n(i-3,j,k) )

             Krhov(i,j,k) = d7(0)*rhov_n(i,j,k) &
                  + d7(1)*( rhov_n(i+1,j,k)+rhov_n(i-1,j,k) )&
                  + d7(2)*( rhov_n(i+2,j,k)+rhov_n(i-2,j,k) )&
                  + d7(3)*( rhov_n(i+3,j,k)+rhov_n(i-3,j,k) )

             Krhow(i,j,k) = d7(0)*rhow_n(i,j,k) &
                  + d7(1)*( rhow_n(i+1,j,k)+rhow_n(i-1,j,k) )&
                  + d7(2)*( rhow_n(i+2,j,k)+rhow_n(i-2,j,k) )&
                  + d7(3)*( rhow_n(i+3,j,k)+rhow_n(i-3,j,k) )

             Krhoe(i,j,k) = d7(0)*rhoe_n(i,j,k) &
                  + d7(1)*( rhoe_n(i+1,j,k)+rhoe_n(i-1,j,k) )&
                  + d7(2)*( rhoe_n(i+2,j,k)+rhoe_n(i-2,j,k) )&
                  + d7(3)*( rhoe_n(i+3,j,k)+rhoe_n(i-3,j,k) )
          enddo
       enddo
    enddo

    ! BC at imax
    ! ----------
    if (is_boundary(1,2)) then
       ! 4-th order centered filter
       i=nx-2
       do k=1,nz
          do j=1,ny
             Krho(i,j,k)  = d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i+1,j,k)+d5(2)*rho_n(i+2,j,k) &
                  + d5(1)*rho_n(i-1,j,k)+d5(2)*rho_n(i-2,j,k)
             Krhou(i,j,k) = d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i+1,j,k)+d5(2)*rhou_n(i+2,j,k) &
                  + d5(1)*rhou_n(i-1,j,k)+d5(2)*rhou_n(i-2,j,k)
             Krhov(i,j,k) = d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i+1,j,k)+d5(2)*rhov_n(i+2,j,k) &
                  + d5(1)*rhov_n(i-1,j,k)+d5(2)*rhov_n(i-2,j,k)
             Krhow(i,j,k) = d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i+1,j,k)+d5(2)*rhow_n(i+2,j,k) &
                  + d5(1)*rhow_n(i-1,j,k)+d5(2)*rhow_n(i-2,j,k)
             Krhoe(i,j,k) = d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i+1,j,k)+d5(2)*rhoe_n(i+2,j,k) &
                  + d5(1)*rhoe_n(i-1,j,k)+d5(2)*rhoe_n(i-2,j,k)
          enddo
       enddo
       ! no filter for i=nx-1 & i=nx -> initialize array
       do j=1,ny
           Krho(nx-1,j,:)=0.0_wp
          Krhou(nx-1,j,:)=0.0_wp
          Krhov(nx-1,j,:)=0.0_wp
          Krhow(nx-1,j,:)=0.0_wp
          Krhoe(nx-1,j,:)=0.0_wp
          
           Krho(nx,j,:)=0.0_wp
          Krhou(nx,j,:)=0.0_wp
          Krhov(nx,j,:)=0.0_wp
          Krhow(nx,j,:)=0.0_wp
          Krhoe(nx,j,:)=0.0_wp
       enddo
    endif

    ! Filtering along j-direction
    ! ===========================
    
    ! BC at jmin
    ! ----------
    if (is_boundary(2,1)) then
       ! 4-th order centered filter
       j=3
       do k=1,nz
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i,j+1,k)+d5(2)*rho_n(i,j+2,k) &
                  + d5(1)*rho_n(i,j-1,k)+d5(2)*rho_n(i,j-2,k)
             Krhou(i,j,k) = Krhou(i,j,k)  +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j+1,k)+d5(2)*rhou_n(i,j+2,k) &
                  + d5(1)*rhou_n(i,j-1,k)+d5(2)*rhou_n(i,j-2,k)
             Krhov(i,j,k) = Krhov(i,j,k)  +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j+1,k)+d5(2)*rhov_n(i,j+2,k) &
                  + d5(1)*rhov_n(i,j-1,k)+d5(2)*rhov_n(i,j-2,k)
             Krhow(i,j,k) = Krhow(i,j,k)  +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j+1,k)+d5(2)*rhow_n(i,j+2,k) &
                  + d5(1)*rhow_n(i,j-1,k)+d5(2)*rhow_n(i,j-2,k)
             Krhoe(i,j,k) = Krhoe(i,j,k)  +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i,j+1,k)+d5(2)*rhoe_n(i,j+2,k) &
                  + d5(1)*rhoe_n(i,j-1,k)+d5(2)*rhoe_n(i,j-2,k)
          enddo
       enddo
    endif

    ! Interior points
    ! ---------------
    ! 6-th order centered filter
    do k=1,nz
       do j=ndy,nfy
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  + d7(0)*rho_n(i,j,k) &
                  + d7(1)*( rho_n(i,j+1,k)+rho_n(i,j-1,k) )&
                  + d7(2)*( rho_n(i,j+2,k)+rho_n(i,j-2,k) )&
                  + d7(3)*( rho_n(i,j+3,k)+rho_n(i,j-3,k) )

             Krhou(i,j,k) = Krhou(i,j,k)  + d7(0)*rhou_n(i,j,k) &
                  + d7(1)*( rhou_n(i,j+1,k)+rhou_n(i,j-1,k) )&
                  + d7(2)*( rhou_n(i,j+2,k)+rhou_n(i,j-2,k) )&
                  + d7(3)*( rhou_n(i,j+3,k)+rhou_n(i,j-3,k) )

             Krhov(i,j,k) = Krhov(i,j,k)  + d7(0)*rhov_n(i,j,k) &
                  + d7(1)*( rhov_n(i,j+1,k)+rhov_n(i,j-1,k) )&
                  + d7(2)*( rhov_n(i,j+2,k)+rhov_n(i,j-2,k) )&
                  + d7(3)*( rhov_n(i,j+3,k)+rhov_n(i,j-3,k) )

             Krhow(i,j,k) = Krhow(i,j,k)  + d7(0)*rhow_n(i,j,k) &
                  + d7(1)*( rhow_n(i,j+1,k)+rhow_n(i,j-1,k) )&
                  + d7(2)*( rhow_n(i,j+2,k)+rhow_n(i,j-2,k) )&
                  + d7(3)*( rhow_n(i,j+3,k)+rhow_n(i,j-3,k) )

             Krhoe(i,j,k) = Krhoe(i,j,k)  + d7(0)*rhoe_n(i,j,k) &
                  + d7(1)*( rhoe_n(i,j+1,k)+rhoe_n(i,j-1,k) )&
                  + d7(2)*( rhoe_n(i,j+2,k)+rhoe_n(i,j-2,k) )&
                  + d7(3)*( rhoe_n(i,j+3,k)+rhoe_n(i,j-3,k) )
          enddo
       enddo
    enddo

    ! BC at jmax
    ! ----------
    if (is_boundary(2,2)) then
       ! 4-th order centered filter
       j=ny-2
       do k=1,nz
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i,j+1,k)+d5(2)*rho_n(i,j+2,k) &
                  + d5(1)*rho_n(i,j-1,k)+d5(2)*rho_n(i,j-2,k)
             Krhou(i,j,k) = Krhou(i,j,k)  +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j+1,k)+d5(2)*rhou_n(i,j+2,k) &
                  + d5(1)*rhou_n(i,j-1,k)+d5(2)*rhou_n(i,j-2,k)
             Krhov(i,j,k) = Krhov(i,j,k)  +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j+1,k)+d5(2)*rhov_n(i,j+2,k) &
                  + d5(1)*rhov_n(i,j-1,k)+d5(2)*rhov_n(i,j-2,k)
             Krhow(i,j,k) = Krhow(i,j,k)  +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j+1,k)+d5(2)*rhow_n(i,j+2,k) &
                  + d5(1)*rhow_n(i,j-1,k)+d5(2)*rhow_n(i,j-2,k)
             Krhoe(i,j,k) = Krhoe(i,j,k)  +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i,j+1,k)+d5(2)*rhoe_n(i,j+2,k) &
                  + d5(1)*rhoe_n(i,j-1,k)+d5(2)*rhoe_n(i,j-2,k)
          enddo
       enddo
    endif

    !****************
    if (is_2D) return
    !****************

    ! Filtering along k-direction
    ! ===========================
    
    ! BC at kmin
    ! ----------
    if (is_boundary(3,1)) then
       ! 4-th order centered filter
       k=3
       do j=1,ny
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i,j,k+1)+d5(2)*rho_n(i,j,k+2) &
                  + d5(1)*rho_n(i,j,k-1)+d5(2)*rho_n(i,j,k-2)
             Krhou(i,j,k) = Krhou(i,j,k)  +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j,k+1)+d5(2)*rhou_n(i,j,k+2) &
                  + d5(1)*rhou_n(i,j,k-1)+d5(2)*rhou_n(i,j,k-2)
             Krhov(i,j,k) = Krhov(i,j,k)  +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j,k+1)+d5(2)*rhov_n(i,j,k+2) &
                  + d5(1)*rhov_n(i,j,k-1)+d5(2)*rhov_n(i,j,k-2)
             Krhow(i,j,k) = Krhow(i,j,k)  +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j,k+1)+d5(2)*rhow_n(i,j,k+2) &
                  + d5(1)*rhow_n(i,j,k-1)+d5(2)*rhow_n(i,j,k-2)
             Krhoe(i,j,k) = Krhoe(i,j,k)  +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i,j,k+1)+d5(2)*rhoe_n(i,j,k+2) &
                  + d5(1)*rhoe_n(i,j,k-1)+d5(2)*rhoe_n(i,j,k-2)
          enddo
       enddo
    endif

    ! Interior points
    ! ---------------
    ! 6-th order centered filter
    do k=ndz,nfz
       do j=1,ny
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  + d7(0)*rho_n(i,j,k) &
                  + d7(1)*( rho_n(i,j,k+1)+rho_n(i,j,k-1) )&
                  + d7(2)*( rho_n(i,j,k+2)+rho_n(i,j,k-2) )&
                  + d7(3)*( rho_n(i,j,k+3)+rho_n(i,j,k-3) )

             Krhou(i,j,k) = Krhou(i,j,k)  + d7(0)*rhou_n(i,j,k) &
                  + d7(1)*( rhou_n(i,j,k+1)+rhou_n(i,j,k-1) )&
                  + d7(2)*( rhou_n(i,j,k+2)+rhou_n(i,j,k-2) )&
                  + d7(3)*( rhou_n(i,j,k+3)+rhou_n(i,j,k-3) )

             Krhov(i,j,k) = Krhov(i,j,k)  + d7(0)*rhov_n(i,j,k) &
                  + d7(1)*( rhov_n(i,j,k+1)+rhov_n(i,j,k-1) )&
                  + d7(2)*( rhov_n(i,j,k+2)+rhov_n(i,j,k-2) )&
                  + d7(3)*( rhov_n(i,j,k+3)+rhov_n(i,j,k-3) )

             Krhow(i,j,k) = Krhow(i,j,k)  + d7(0)*rhow_n(i,j,k) &
                  + d7(1)*( rhow_n(i,j,k+1)+rhow_n(i,j,k-1) )&
                  + d7(2)*( rhow_n(i,j,k+2)+rhow_n(i,j,k-2) )&
                  + d7(3)*( rhow_n(i,j,k+3)+rhow_n(i,j,k-3) )

             Krhoe(i,j,k) = Krhoe(i,j,k)  + d7(0)*rhoe_n(i,j,k) &
                  + d7(1)*( rhoe_n(i,j,k+1)+rhoe_n(i,j,k-1) )&
                  + d7(2)*( rhoe_n(i,j,k+2)+rhoe_n(i,j,k-2) )&
                  + d7(3)*( rhoe_n(i,j,k+3)+rhoe_n(i,j,k-3) )
          enddo
       enddo
    enddo

    ! BC at kmax
    ! ----------
    if (is_boundary(3,2)) then
       ! 4-th order centered filter
       k=nz-2
       do j=1,ny
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i,j,k+1)+d5(2)*rho_n(i,j,k+2) &
                  + d5(1)*rho_n(i,j,k-1)+d5(2)*rho_n(i,j,k-2)
             Krhou(i,j,k) = Krhou(i,j,k)  +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j,k+1)+d5(2)*rhou_n(i,j,k+2) &
                  + d5(1)*rhou_n(i,j,k-1)+d5(2)*rhou_n(i,j,k-2)
             Krhov(i,j,k) = Krhov(i,j,k)  +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j,k+1)+d5(2)*rhov_n(i,j,k+2) &
                  + d5(1)*rhov_n(i,j,k-1)+d5(2)*rhov_n(i,j,k-2)
             Krhow(i,j,k) = Krhow(i,j,k)  +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j,k+1)+d5(2)*rhow_n(i,j,k+2) &
                  + d5(1)*rhow_n(i,j,k-1)+d5(2)*rhow_n(i,j,k-2)
             Krhoe(i,j,k) = Krhoe(i,j,k)  +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i,j,k+1)+d5(2)*rhoe_n(i,j,k+2) &
                  + d5(1)*rhoe_n(i,j,k-1)+d5(2)*rhoe_n(i,j,k-2)
          enddo
       enddo
    endif

  end subroutine filtering_7pts

  !==============================================================================
  subroutine filtering_5pts
  !==============================================================================
    !> Apply filtering on 5-point stencils (+ boundaries)
    !> * old version without comput./comm. overlap *
  !==============================================================================
    use mod_flow
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    ! ---------------------------------------------------------------------------

    ! Filtering along i-direction
    ! ===========================
    
    ! BC at imin
    ! ----------
    if (is_boundary(1,1)) then
       ! no filter for i=1 & i=2 -> initialize array
       do j=1,ny
           Krho(1,j,:)=0.0_wp
          Krhou(1,j,:)=0.0_wp
          Krhov(1,j,:)=0.0_wp
          Krhow(1,j,:)=0.0_wp
          Krhoe(1,j,:)=0.0_wp
          
           Krho(2,j,:)=0.0_wp
          Krhou(2,j,:)=0.0_wp
          Krhov(2,j,:)=0.0_wp
          Krhow(2,j,:)=0.0_wp
          Krhoe(2,j,:)=0.0_wp
       enddo
    endif

    ! Interior points
    ! ---------------
    ! 4-th order centered filter
    do k=1,nz
       do j=1,ny
          do i=ndx,nfx
             Krho(i,j,k)  = d5(0)*rho_n(i,j,k) &
                  + d5(1)*( rho_n(i+1,j,k)+rho_n(i-1,j,k) )&
                  + d5(2)*( rho_n(i+2,j,k)+rho_n(i-2,j,k) )

             Krhou(i,j,k) = d5(0)*rhou_n(i,j,k) &
                  + d5(1)*( rhou_n(i+1,j,k)+rhou_n(i-1,j,k) )&
                  + d5(2)*( rhou_n(i+2,j,k)+rhou_n(i-2,j,k) )

             Krhov(i,j,k) = d5(0)*rhov_n(i,j,k) &
                  + d5(1)*( rhov_n(i+1,j,k)+rhov_n(i-1,j,k) )&
                  + d5(2)*( rhov_n(i+2,j,k)+rhov_n(i-2,j,k) )

             Krhow(i,j,k) = d5(0)*rhow_n(i,j,k) &
                  + d5(1)*( rhow_n(i+1,j,k)+rhow_n(i-1,j,k) )&
                  + d5(2)*( rhow_n(i+2,j,k)+rhow_n(i-2,j,k) )

             Krhoe(i,j,k) = d5(0)*rhoe_n(i,j,k) &
                  + d5(1)*( rhoe_n(i+1,j,k)+rhoe_n(i-1,j,k) )&
                  + d5(2)*( rhoe_n(i+2,j,k)+rhoe_n(i-2,j,k) )
          enddo
       enddo
    enddo

    ! BC at imax
    ! ----------
    if (is_boundary(1,2)) then
       ! no filter for i=nx-1 & i=nx -> initialize array
       do j=1,ny
           Krho(nx-1,j,:)=0.0_wp
          Krhou(nx-1,j,:)=0.0_wp
          Krhov(nx-1,j,:)=0.0_wp
          Krhow(nx-1,j,:)=0.0_wp
          Krhoe(nx-1,j,:)=0.0_wp
          
           Krho(nx,j,:)=0.0_wp
          Krhou(nx,j,:)=0.0_wp
          Krhov(nx,j,:)=0.0_wp
          Krhow(nx,j,:)=0.0_wp
          Krhoe(nx,j,:)=0.0_wp
       enddo
    endif

    ! Filtering along j-direction
    ! ===========================
    
    ! Interior points
    ! ---------------
    ! 4-th order centered filter
    do k=1,nz
       do j=ndy,nfy
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  + d5(0)*rho_n(i,j,k) &
                  + d5(1)*( rho_n(i,j+1,k)+rho_n(i,j-1,k) )&
                  + d5(2)*( rho_n(i,j+2,k)+rho_n(i,j-2,k) )

             Krhou(i,j,k) = Krhou(i,j,k)  + d5(0)*rhou_n(i,j,k) &
                  + d5(1)*( rhou_n(i,j+1,k)+rhou_n(i,j-1,k) )&
                  + d5(2)*( rhou_n(i,j+2,k)+rhou_n(i,j-2,k) )

             Krhov(i,j,k) = Krhov(i,j,k)  + d5(0)*rhov_n(i,j,k) &
                  + d5(1)*( rhov_n(i,j+1,k)+rhov_n(i,j-1,k) )&
                  + d5(2)*( rhov_n(i,j+2,k)+rhov_n(i,j-2,k) )

             Krhow(i,j,k) = Krhow(i,j,k)  + d5(0)*rhow_n(i,j,k) &
                  + d5(1)*( rhow_n(i,j+1,k)+rhow_n(i,j-1,k) )&
                  + d5(2)*( rhow_n(i,j+2,k)+rhow_n(i,j-2,k) )

             Krhoe(i,j,k) = Krhoe(i,j,k)  + d5(0)*rhoe_n(i,j,k) &
                  + d5(1)*( rhoe_n(i,j+1,k)+rhoe_n(i,j-1,k) )&
                  + d5(2)*( rhoe_n(i,j+2,k)+rhoe_n(i,j-2,k) )
          enddo
       enddo
    enddo

    !****************
    if (is_2D) return
    !****************

    ! Filtering along k-direction
    ! ===========================

    ! Interior points
    ! ---------------
    ! 4-th order centered filter
    do k=ndz,nfz
       do j=1,ny
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  + d5(0)*rho_n(i,j,k) &
                  + d5(1)*( rho_n(i,j,k+1)+rho_n(i,j,k-1) )&
                  + d5(2)*( rho_n(i,j,k+2)+rho_n(i,j,k-2) )

             Krhou(i,j,k) = Krhou(i,j,k)  + d5(0)*rhou_n(i,j,k) &
                  + d5(1)*( rhou_n(i,j,k+1)+rhou_n(i,j,k-1) )&
                  + d5(2)*( rhou_n(i,j,k+2)+rhou_n(i,j,k-2) )

             Krhov(i,j,k) = Krhov(i,j,k)  + d5(0)*rhov_n(i,j,k) &
                  + d5(1)*( rhov_n(i,j,k+1)+rhov_n(i,j,k-1) )&
                  + d5(2)*( rhov_n(i,j,k+2)+rhov_n(i,j,k-2) )

             Krhow(i,j,k) = Krhow(i,j,k)  + d5(0)*rhow_n(i,j,k) &
                  + d5(1)*( rhow_n(i,j,k+1)+rhow_n(i,j,k-1) )&
                  + d5(2)*( rhow_n(i,j,k+2)+rhow_n(i,j,k-2) )

             Krhoe(i,j,k) = Krhoe(i,j,k)  + d5(0)*rhoe_n(i,j,k) &
                  + d5(1)*( rhoe_n(i,j,k+1)+rhoe_n(i,j,k-1) )&
                  + d5(2)*( rhoe_n(i,j,k+2)+rhoe_n(i,j,k-2) )
          enddo
       enddo
    enddo

  end subroutine filtering_5pts

  !==============================================================================
  subroutine filtering_11pts_periodic
  !==============================================================================
    !> Apply filtering on 11-point stencils (PERIODIC VERSION)
    !> * old version without comput./comm. overlap *
  !==============================================================================
    use mod_flow
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    ! ---------------------------------------------------------------------------

    ! Filtering along i-direction
    ! ===========================
    do k=1,nz
       do j=1,ny
          do i=1,nx
             Krho(i,j,k)  = d11(0)*rho_n(i,j,k) &
                  + d11(1)*( rho_n(i+1,j,k)+rho_n(i-1,j,k) )&
                  + d11(2)*( rho_n(i+2,j,k)+rho_n(i-2,j,k) )&
                  + d11(3)*( rho_n(i+3,j,k)+rho_n(i-3,j,k) )&
                  + d11(4)*( rho_n(i+4,j,k)+rho_n(i-4,j,k) )&
                  + d11(5)*( rho_n(i+5,j,k)+rho_n(i-5,j,k) )

             Krhou(i,j,k) = d11(0)*rhou_n(i,j,k) &
                  + d11(1)*( rhou_n(i+1,j,k)+rhou_n(i-1,j,k) )&
                  + d11(2)*( rhou_n(i+2,j,k)+rhou_n(i-2,j,k) )&
                  + d11(3)*( rhou_n(i+3,j,k)+rhou_n(i-3,j,k) )&
                  + d11(4)*( rhou_n(i+4,j,k)+rhou_n(i-4,j,k) )&
                  + d11(5)*( rhou_n(i+5,j,k)+rhou_n(i-5,j,k) )

             Krhov(i,j,k) = d11(0)*rhov_n(i,j,k) &
                  + d11(1)*( rhov_n(i+1,j,k)+rhov_n(i-1,j,k) )&
                  + d11(2)*( rhov_n(i+2,j,k)+rhov_n(i-2,j,k) )&
                  + d11(3)*( rhov_n(i+3,j,k)+rhov_n(i-3,j,k) )&
                  + d11(4)*( rhov_n(i+4,j,k)+rhov_n(i-4,j,k) )&
                  + d11(5)*( rhov_n(i+5,j,k)+rhov_n(i-5,j,k) )

             Krhow(i,j,k) = d11(0)*rhow_n(i,j,k) &
                  + d11(1)*( rhow_n(i+1,j,k)+rhow_n(i-1,j,k) )&
                  + d11(2)*( rhow_n(i+2,j,k)+rhow_n(i-2,j,k) )&
                  + d11(3)*( rhow_n(i+3,j,k)+rhow_n(i-3,j,k) )&
                  + d11(4)*( rhow_n(i+4,j,k)+rhow_n(i-4,j,k) )&
                  + d11(5)*( rhow_n(i+5,j,k)+rhow_n(i-5,j,k) )

             Krhoe(i,j,k) = d11(0)*rhoe_n(i,j,k) &
                  + d11(1)*( rhoe_n(i+1,j,k)+rhoe_n(i-1,j,k) )&
                  + d11(2)*( rhoe_n(i+2,j,k)+rhoe_n(i-2,j,k) )&
                  + d11(3)*( rhoe_n(i+3,j,k)+rhoe_n(i-3,j,k) )&
                  + d11(4)*( rhoe_n(i+4,j,k)+rhoe_n(i-4,j,k) )&
                  + d11(5)*( rhoe_n(i+5,j,k)+rhoe_n(i-5,j,k) )
          enddo
       enddo
    enddo

    ! Filtering along j-direction
    ! ===========================
    do k=1,nz
       do j=1,ny
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  + d11(0)*rho_n(i,j,k) &
                  + d11(1)*( rho_n(i,j+1,k)+rho_n(i,j-1,k) )&
                  + d11(2)*( rho_n(i,j+2,k)+rho_n(i,j-2,k) )&
                  + d11(3)*( rho_n(i,j+3,k)+rho_n(i,j-3,k) )&
                  + d11(4)*( rho_n(i,j+4,k)+rho_n(i,j-4,k) )&
                  + d11(5)*( rho_n(i,j+5,k)+rho_n(i,j-5,k) )

             Krhou(i,j,k) = Krhou(i,j,k)  + d11(0)*rhou_n(i,j,k) &
                  + d11(1)*( rhou_n(i,j+1,k)+rhou_n(i,j-1,k) )&
                  + d11(2)*( rhou_n(i,j+2,k)+rhou_n(i,j-2,k) )&
                  + d11(3)*( rhou_n(i,j+3,k)+rhou_n(i,j-3,k) )&
                  + d11(4)*( rhou_n(i,j+4,k)+rhou_n(i,j-4,k) )&
                  + d11(5)*( rhou_n(i,j+5,k)+rhou_n(i,j-5,k) )

             Krhov(i,j,k) = Krhov(i,j,k)  + d11(0)*rhov_n(i,j,k) &
                  + d11(1)*( rhov_n(i,j+1,k)+rhov_n(i,j-1,k) )&
                  + d11(2)*( rhov_n(i,j+2,k)+rhov_n(i,j-2,k) )&
                  + d11(3)*( rhov_n(i,j+3,k)+rhov_n(i,j-3,k) )&
                  + d11(4)*( rhov_n(i,j+4,k)+rhov_n(i,j-4,k) )&
                  + d11(5)*( rhov_n(i,j+5,k)+rhov_n(i,j-5,k) )

             Krhow(i,j,k) = Krhow(i,j,k)  + d11(0)*rhow_n(i,j,k) &
                  + d11(1)*( rhow_n(i,j+1,k)+rhow_n(i,j-1,k) )&
                  + d11(2)*( rhow_n(i,j+2,k)+rhow_n(i,j-2,k) )&
                  + d11(3)*( rhow_n(i,j+3,k)+rhow_n(i,j-3,k) )&
                  + d11(4)*( rhow_n(i,j+4,k)+rhow_n(i,j-4,k) )&
                  + d11(5)*( rhow_n(i,j+5,k)+rhow_n(i,j-5,k) )

             Krhoe(i,j,k) = Krhoe(i,j,k)  + d11(0)*rhoe_n(i,j,k) &
                  + d11(1)*( rhoe_n(i,j+1,k)+rhoe_n(i,j-1,k) )&
                  + d11(2)*( rhoe_n(i,j+2,k)+rhoe_n(i,j-2,k) )&
                  + d11(3)*( rhoe_n(i,j+3,k)+rhoe_n(i,j-3,k) )&
                  + d11(4)*( rhoe_n(i,j+4,k)+rhoe_n(i,j-4,k) )&
                  + d11(5)*( rhoe_n(i,j+5,k)+rhoe_n(i,j-5,k) )
          enddo
       enddo
    enddo

    !****************
    if (is_2D) return
    !****************

    ! Filtering along k-direction
    ! ===========================
    do k=1,nz
       do j=1,ny
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  + d11(0)*rho_n(i,j,k) &
                  + d11(1)*( rho_n(i,j,k+1)+rho_n(i,j,k-1) )&
                  + d11(2)*( rho_n(i,j,k+2)+rho_n(i,j,k-2) )&
                  + d11(3)*( rho_n(i,j,k+3)+rho_n(i,j,k-3) )&
                  + d11(4)*( rho_n(i,j,k+4)+rho_n(i,j,k-4) )&
                  + d11(5)*( rho_n(i,j,k+5)+rho_n(i,j,k-5) )

             Krhou(i,j,k) = Krhou(i,j,k)  + d11(0)*rhou_n(i,j,k) &
                  + d11(1)*( rhou_n(i,j,k+1)+rhou_n(i,j,k-1) )&
                  + d11(2)*( rhou_n(i,j,k+2)+rhou_n(i,j,k-2) )&
                  + d11(3)*( rhou_n(i,j,k+3)+rhou_n(i,j,k-3) )&
                  + d11(4)*( rhou_n(i,j,k+4)+rhou_n(i,j,k-4) )&
                  + d11(5)*( rhou_n(i,j,k+5)+rhou_n(i,j,k-5) )

             Krhov(i,j,k) = Krhov(i,j,k)  + d11(0)*rhov_n(i,j,k) &
                  + d11(1)*( rhov_n(i,j,k+1)+rhov_n(i,j,k-1) )&
                  + d11(2)*( rhov_n(i,j,k+2)+rhov_n(i,j,k-2) )&
                  + d11(3)*( rhov_n(i,j,k+3)+rhov_n(i,j,k-3) )&
                  + d11(4)*( rhov_n(i,j,k+4)+rhov_n(i,j,k-4) )&
                  + d11(5)*( rhov_n(i,j,k+5)+rhov_n(i,j,k-5) )

             Krhow(i,j,k) = Krhow(i,j,k)  + d11(0)*rhow_n(i,j,k) &
                  + d11(1)*( rhow_n(i,j,k+1)+rhow_n(i,j,k-1) )&
                  + d11(2)*( rhow_n(i,j,k+2)+rhow_n(i,j,k-2) )&
                  + d11(3)*( rhow_n(i,j,k+3)+rhow_n(i,j,k-3) )&
                  + d11(4)*( rhow_n(i,j,k+4)+rhow_n(i,j,k-4) )&
                  + d11(5)*( rhow_n(i,j,k+5)+rhow_n(i,j,k-5) )

             Krhoe(i,j,k) = Krhoe(i,j,k)  + d11(0)*rhoe_n(i,j,k) &
                  + d11(1)*( rhoe_n(i,j,k+1)+rhoe_n(i,j,k-1) )&
                  + d11(2)*( rhoe_n(i,j,k+2)+rhoe_n(i,j,k-2) )&
                  + d11(3)*( rhoe_n(i,j,k+3)+rhoe_n(i,j,k-3) )&
                  + d11(4)*( rhoe_n(i,j,k+4)+rhoe_n(i,j,k-4) )&
                  + d11(5)*( rhoe_n(i,j,k+5)+rhoe_n(i,j,k-5) )
          enddo
       enddo
    enddo

  end subroutine filtering_11pts_periodic

end module mod_filtering

