! ===============================================================================
module mod_filtering_inc
! ===============================================================================
  !> Module for filtering variables
  !> Added to increments (for IRS smoothing)
! ===============================================================================
  use mod_flow
  use mod_filtering
  implicit none
  ! -----------------------------------------------------------------------------
  ! declaration of coefficients in mod_filtering.f90
  ! -----------------------------------------------------------------------------

contains

  !==============================================================================
  subroutine init_filtering_inc(stencil,is_DRP,X_sf)
  !==============================================================================
    !> Initialization of selective filters
    !> Added to increments (for IRS smoothing)
  !==============================================================================
    use mod_time
    implicit none
    ! ---------------------------------------------------------------------------
    ! Scheme parameters
    integer, intent(in)  :: stencil
    logical, intent(in)  :: is_DRP
    real(wp), intent(in) :: X_sf
    ! ---------------------------------------------------------------------------

    ! Fill selective filter coefficients
    ! ==================================
    call init_coeff_filters(stencil,is_DRP)

    ! Multiply coeff by filtering amplitude (chi_SF) and divide by deltat
    ! ===================================================================
    d3 =d3*X_sf/deltat
    d5 =d5*X_sf/deltat
    d7 =d7*X_sf/deltat
    d9 =d9*X_sf/deltat
    d11=d11*X_sf/deltat
    d13=d13*X_sf/deltat

    d15=d15*X_sf/deltat
    d24=d24*X_sf/deltat
    d28=d28*X_sf/deltat
    d37=d37*X_sf/deltat
    d46=d46*X_sf/deltat

  end subroutine init_filtering_inc

  !===============================================================================
  subroutine filtering_11pts_inc
  !===============================================================================
    !> Apply filtering on 11-point stencils (+ boundaries)
    !> Added to increments (for IRS smoothing)
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    ! ---------------------------------------------------------------------------

    ! Filtering along i-direction
    ! ===========================
    
    ! BC at imin
    ! ----------
    if (is_boundary(1,1)) then
       ! 4-th order centered filter
       i=3
       do k=1,nz
          do j=1,ny
             Krho(i,j,k)  = Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i+1,j,k) +d5(2)*rho_n(i+2,j,k) &
                  + d5(1)*rho_n(i-1,j,k) +d5(2)*rho_n(i-2,j,k)
             Krhou(i,j,k) = Krhou(i,j,k) +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i+1,j,k)+d5(2)*rhou_n(i+2,j,k) &
                  + d5(1)*rhou_n(i-1,j,k)+d5(2)*rhou_n(i-2,j,k)
             Krhov(i,j,k) = Krhov(i,j,k) +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i+1,j,k)+d5(2)*rhov_n(i+2,j,k) &
                  + d5(1)*rhov_n(i-1,j,k)+d5(2)*rhov_n(i-2,j,k)
             Krhow(i,j,k) = Krhow(i,j,k) +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i+1,j,k)+d5(2)*rhow_n(i+2,j,k) &
                  + d5(1)*rhow_n(i-1,j,k)+d5(2)*rhow_n(i-2,j,k)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i+1,j,k)+d5(2)*rhoe_n(i+2,j,k) &
                  + d5(1)*rhoe_n(i-1,j,k)+d5(2)*rhoe_n(i-2,j,k)
          enddo
       enddo
       ! 6-th order centered filter
       i=4
       do k=1,nz
          do j=1,ny
             Krho(i,j,k)  = Krho(i,j,k)  +d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i+1,j,k) +d7(2)*rho_n(i+2,j,k)+d7(3)*rho_n(i+3,j,k) &
                  + d7(1)*rho_n(i-1,j,k) +d7(2)*rho_n(i-2,j,k)+d7(3)*rho_n(i-3,j,k)
             Krhou(i,j,k) = Krhou(i,j,k) +d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i+1,j,k)+d7(2)*rhou_n(i+2,j,k)+d7(3)*rhou_n(i+3,j,k) &
                  + d7(1)*rhou_n(i-1,j,k)+d7(2)*rhou_n(i-2,j,k)+d7(3)*rhou_n(i-3,j,k)
             Krhov(i,j,k) = Krhov(i,j,k) +d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i+1,j,k)+d7(2)*rhov_n(i+2,j,k)+d7(3)*rhov_n(i+3,j,k) &
                  + d7(1)*rhov_n(i-1,j,k)+d7(2)*rhov_n(i-2,j,k)+d7(3)*rhov_n(i-3,j,k)
             Krhow(i,j,k) = Krhow(i,j,k) +d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i+1,j,k)+d7(2)*rhow_n(i+2,j,k)+d7(3)*rhow_n(i+3,j,k) &
                  + d7(1)*rhow_n(i-1,j,k)+d7(2)*rhow_n(i-2,j,k)+d7(3)*rhow_n(i-3,j,k)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i+1,j,k)+d7(2)*rhoe_n(i+2,j,k)+d7(3)*rhoe_n(i+3,j,k) &
                  + d7(1)*rhoe_n(i-1,j,k)+d7(2)*rhoe_n(i-2,j,k)+d7(3)*rhoe_n(i-3,j,k)
          enddo
       enddo
       ! 8-th order centered filter
       i=5
       do k=1,nz
          do j=1,ny
             Krho(i,j,k)  = Krho(i,j,k)  +d9(0)*rho_n(i,j,k)   &
                  + d9(1)*rho_n(i+1,j,k) +d9(2)*rho_n(i+2,j,k)+d9(3)*rho_n(i+3,j,k)+d9(4)*rho_n(i+4,j,k) &
                  + d9(1)*rho_n(i-1,j,k) +d9(2)*rho_n(i-2,j,k)+d9(3)*rho_n(i-3,j,k)+d9(4)*rho_n(i-4,j,k)
             Krhou(i,j,k) = Krhou(i,j,k) +d9(0)*rhou_n(i,j,k)   &
                  + d9(1)*rhou_n(i+1,j,k)+d9(2)*rhou_n(i+2,j,k)+d9(3)*rhou_n(i+3,j,k)+d9(4)*rhou_n(i+4,j,k) &
                  + d9(1)*rhou_n(i-1,j,k)+d9(2)*rhou_n(i-2,j,k)+d9(3)*rhou_n(i-3,j,k)+d9(4)*rhou_n(i-4,j,k)
             Krhov(i,j,k) = Krhov(i,j,k) +d9(0)*rhov_n(i,j,k)   &
                  + d9(1)*rhov_n(i+1,j,k)+d9(2)*rhov_n(i+2,j,k)+d9(3)*rhov_n(i+3,j,k)+d9(4)*rhov_n(i+4,j,k) &
                  + d9(1)*rhov_n(i-1,j,k)+d9(2)*rhov_n(i-2,j,k)+d9(3)*rhov_n(i-3,j,k)+d9(4)*rhov_n(i-4,j,k)
             Krhow(i,j,k) = Krhow(i,j,k) +d9(0)*rhow_n(i,j,k)   &
                  + d9(1)*rhow_n(i+1,j,k)+d9(2)*rhow_n(i+2,j,k)+d9(3)*rhow_n(i+3,j,k)+d9(4)*rhow_n(i+4,j,k) &
                  + d9(1)*rhow_n(i-1,j,k)+d9(2)*rhow_n(i-2,j,k)+d9(3)*rhow_n(i-3,j,k)+d9(4)*rhow_n(i-4,j,k)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d9(0)*rhoe_n(i,j,k)   &
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
             Krho(i,j,k)  = Krho(i,j,k)  + d11(0)*rho_n(i,j,k) &
                  + d11(1)*( rho_n(i+1,j,k)+rho_n(i-1,j,k) )&
                  + d11(2)*( rho_n(i+2,j,k)+rho_n(i-2,j,k) )&
                  + d11(3)*( rho_n(i+3,j,k)+rho_n(i-3,j,k) )&
                  + d11(4)*( rho_n(i+4,j,k)+rho_n(i-4,j,k) )&
                  + d11(5)*( rho_n(i+5,j,k)+rho_n(i-5,j,k) )

             Krhou(i,j,k) = Krhou(i,j,k)  + d11(0)*rhou_n(i,j,k) &
                  + d11(1)*( rhou_n(i+1,j,k)+rhou_n(i-1,j,k) )&
                  + d11(2)*( rhou_n(i+2,j,k)+rhou_n(i-2,j,k) )&
                  + d11(3)*( rhou_n(i+3,j,k)+rhou_n(i-3,j,k) )&
                  + d11(4)*( rhou_n(i+4,j,k)+rhou_n(i-4,j,k) )&
                  + d11(5)*( rhou_n(i+5,j,k)+rhou_n(i-5,j,k) )

             Krhov(i,j,k) = Krhov(i,j,k)  + d11(0)*rhov_n(i,j,k) &
                  + d11(1)*( rhov_n(i+1,j,k)+rhov_n(i-1,j,k) )&
                  + d11(2)*( rhov_n(i+2,j,k)+rhov_n(i-2,j,k) )&
                  + d11(3)*( rhov_n(i+3,j,k)+rhov_n(i-3,j,k) )&
                  + d11(4)*( rhov_n(i+4,j,k)+rhov_n(i-4,j,k) )&
                  + d11(5)*( rhov_n(i+5,j,k)+rhov_n(i-5,j,k) )

             Krhow(i,j,k) = Krhow(i,j,k)  + d11(0)*rhow_n(i,j,k) &
                  + d11(1)*( rhow_n(i+1,j,k)+rhow_n(i-1,j,k) )&
                  + d11(2)*( rhow_n(i+2,j,k)+rhow_n(i-2,j,k) )&
                  + d11(3)*( rhow_n(i+3,j,k)+rhow_n(i-3,j,k) )&
                  + d11(4)*( rhow_n(i+4,j,k)+rhow_n(i-4,j,k) )&
                  + d11(5)*( rhow_n(i+5,j,k)+rhow_n(i-5,j,k) )

             Krhoe(i,j,k) = Krhoe(i,j,k)  + d11(0)*rhoe_n(i,j,k) &
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
             Krho(i,j,k)  = Krho(i,j,k)  +d9(0)*rho_n(i,j,k)   &
                  + d9(1)*rho_n(i+1,j,k) +d9(2)*rho_n(i+2,j,k)+d9(3)*rho_n(i+3,j,k)+d9(4)*rho_n(i+4,j,k) &
                  + d9(1)*rho_n(i-1,j,k) +d9(2)*rho_n(i-2,j,k)+d9(3)*rho_n(i-3,j,k)+d9(4)*rho_n(i-4,j,k)
             Krhou(i,j,k) = Krhou(i,j,k) +d9(0)*rhou_n(i,j,k)   &
                  + d9(1)*rhou_n(i+1,j,k)+d9(2)*rhou_n(i+2,j,k)+d9(3)*rhou_n(i+3,j,k)+d9(4)*rhou_n(i+4,j,k) &
                  + d9(1)*rhou_n(i-1,j,k)+d9(2)*rhou_n(i-2,j,k)+d9(3)*rhou_n(i-3,j,k)+d9(4)*rhou_n(i-4,j,k)
             Krhov(i,j,k) = Krhov(i,j,k) +d9(0)*rhov_n(i,j,k)   &
                  + d9(1)*rhov_n(i+1,j,k)+d9(2)*rhov_n(i+2,j,k)+d9(3)*rhov_n(i+3,j,k)+d9(4)*rhov_n(i+4,j,k) &
                  + d9(1)*rhov_n(i-1,j,k)+d9(2)*rhov_n(i-2,j,k)+d9(3)*rhov_n(i-3,j,k)+d9(4)*rhov_n(i-4,j,k)
             Krhow(i,j,k) = Krhow(i,j,k) +d9(0)*rhow_n(i,j,k)   &
                  + d9(1)*rhow_n(i+1,j,k)+d9(2)*rhow_n(i+2,j,k)+d9(3)*rhow_n(i+3,j,k)+d9(4)*rhow_n(i+4,j,k) &
                  + d9(1)*rhow_n(i-1,j,k)+d9(2)*rhow_n(i-2,j,k)+d9(3)*rhow_n(i-3,j,k)+d9(4)*rhow_n(i-4,j,k)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d9(0)*rhoe_n(i,j,k)   &
                  + d9(1)*rhoe_n(i+1,j,k)+d9(2)*rhoe_n(i+2,j,k)+d9(3)*rhoe_n(i+3,j,k)+d9(4)*rhoe_n(i+4,j,k) &
                  + d9(1)*rhoe_n(i-1,j,k)+d9(2)*rhoe_n(i-2,j,k)+d9(3)*rhoe_n(i-3,j,k)+d9(4)*rhoe_n(i-4,j,k)
          enddo
       enddo
       ! 6-th order centered filter
       i=nx-3
       do k=1,nz
          do j=1,ny
             Krho(i,j,k)  = Krho(i,j,k)  +d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i+1,j,k) +d7(2)*rho_n(i+2,j,k)+d7(3)*rho_n(i+3,j,k) &
                  + d7(1)*rho_n(i-1,j,k) +d7(2)*rho_n(i-2,j,k)+d7(3)*rho_n(i-3,j,k)
             Krhou(i,j,k) = Krhou(i,j,k) +d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i+1,j,k)+d7(2)*rhou_n(i+2,j,k)+d7(3)*rhou_n(i+3,j,k) &
                  + d7(1)*rhou_n(i-1,j,k)+d7(2)*rhou_n(i-2,j,k)+d7(3)*rhou_n(i-3,j,k)
             Krhov(i,j,k) = Krhov(i,j,k) +d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i+1,j,k)+d7(2)*rhov_n(i+2,j,k)+d7(3)*rhov_n(i+3,j,k) &
                  + d7(1)*rhov_n(i-1,j,k)+d7(2)*rhov_n(i-2,j,k)+d7(3)*rhov_n(i-3,j,k)
             Krhow(i,j,k) = Krhow(i,j,k) +d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i+1,j,k)+d7(2)*rhow_n(i+2,j,k)+d7(3)*rhow_n(i+3,j,k) &
                  + d7(1)*rhow_n(i-1,j,k)+d7(2)*rhow_n(i-2,j,k)+d7(3)*rhow_n(i-3,j,k)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i+1,j,k)+d7(2)*rhoe_n(i+2,j,k)+d7(3)*rhoe_n(i+3,j,k) &
                  + d7(1)*rhoe_n(i-1,j,k)+d7(2)*rhoe_n(i-2,j,k)+d7(3)*rhoe_n(i-3,j,k)
          enddo
       enddo
       ! 4-th order centered filter
       i=nx-2
       do k=1,nz
          do j=1,ny
             Krho(i,j,k)  = Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i+1,j,k) +d5(2)*rho_n(i+2,j,k) &
                  + d5(1)*rho_n(i-1,j,k) +d5(2)*rho_n(i-2,j,k)
             Krhou(i,j,k) = Krhou(i,j,k) +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i+1,j,k)+d5(2)*rhou_n(i+2,j,k) &
                  + d5(1)*rhou_n(i-1,j,k)+d5(2)*rhou_n(i-2,j,k)
             Krhov(i,j,k) = Krhov(i,j,k) +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i+1,j,k)+d5(2)*rhov_n(i+2,j,k) &
                  + d5(1)*rhov_n(i-1,j,k)+d5(2)*rhov_n(i-2,j,k)
             Krhow(i,j,k) = Krhow(i,j,k) +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i+1,j,k)+d5(2)*rhow_n(i+2,j,k) &
                  + d5(1)*rhow_n(i-1,j,k)+d5(2)*rhow_n(i-2,j,k)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i+1,j,k)+d5(2)*rhoe_n(i+2,j,k) &
                  + d5(1)*rhoe_n(i-1,j,k)+d5(2)*rhoe_n(i-2,j,k)
          enddo
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
                  + d7(1)*rho_n(i,j+1,k) +d7(2)*rho_n(i,j+2,k)+d7(3)*rho_n(i,j+3,k) &
                  + d7(1)*rho_n(i,j-1,k) +d7(2)*rho_n(i,j-2,k)+d7(3)*rho_n(i,j-3,k)
             Krhou(i,j,k) = Krhou(i,j,k) +d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i,j+1,k)+d7(2)*rhou_n(i,j+2,k)+d7(3)*rhou_n(i,j+3,k) &
                  + d7(1)*rhou_n(i,j-1,k)+d7(2)*rhou_n(i,j-2,k)+d7(3)*rhou_n(i,j-3,k)
             Krhov(i,j,k) = Krhov(i,j,k) +d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i,j+1,k)+d7(2)*rhov_n(i,j+2,k)+d7(3)*rhov_n(i,j+3,k) &
                  + d7(1)*rhov_n(i,j-1,k)+d7(2)*rhov_n(i,j-2,k)+d7(3)*rhov_n(i,j-3,k)
             Krhow(i,j,k) = Krhow(i,j,k) +d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i,j+1,k)+d7(2)*rhow_n(i,j+2,k)+d7(3)*rhow_n(i,j+3,k) &
                  + d7(1)*rhow_n(i,j-1,k)+d7(2)*rhow_n(i,j-2,k)+d7(3)*rhow_n(i,j-3,k)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i,j+1,k)+d7(2)*rhoe_n(i,j+2,k)+d7(3)*rhoe_n(i,j+3,k) &
                  + d7(1)*rhoe_n(i,j-1,k)+d7(2)*rhoe_n(i,j-2,k)+d7(3)*rhoe_n(i,j-3,k)
          enddo
       enddo
       ! 8-th order centered filter
       j=5
       do k=1,nz
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  +d9(0)*rho_n(i,j,k)   &
                  + d9(1)*rho_n(i,j+1,k) +d9(2)*rho_n(i,j+2,k)+d9(3)*rho_n(i,j+3,k)+d9(4)*rho_n(i,j+4,k) &
                  + d9(1)*rho_n(i,j-1,k) +d9(2)*rho_n(i,j-2,k)+d9(3)*rho_n(i,j-3,k)+d9(4)*rho_n(i,j-4,k)
             Krhou(i,j,k) = Krhou(i,j,k) +d9(0)*rhou_n(i,j,k)   &
                  + d9(1)*rhou_n(i,j+1,k)+d9(2)*rhou_n(i,j+2,k)+d9(3)*rhou_n(i,j+3,k)+d9(4)*rhou_n(i,j+4,k) &
                  + d9(1)*rhou_n(i,j-1,k)+d9(2)*rhou_n(i,j-2,k)+d9(3)*rhou_n(i,j-3,k)+d9(4)*rhou_n(i,j-4,k)
             Krhov(i,j,k) = Krhov(i,j,k) +d9(0)*rhov_n(i,j,k)   &
                  + d9(1)*rhov_n(i,j+1,k)+d9(2)*rhov_n(i,j+2,k)+d9(3)*rhov_n(i,j+3,k)+d9(4)*rhov_n(i,j+4,k) &
                  + d9(1)*rhov_n(i,j-1,k)+d9(2)*rhov_n(i,j-2,k)+d9(3)*rhov_n(i,j-3,k)+d9(4)*rhov_n(i,j-4,k)
             Krhow(i,j,k) = Krhow(i,j,k) +d9(0)*rhow_n(i,j,k)   &
                  + d9(1)*rhow_n(i,j+1,k)+d9(2)*rhow_n(i,j+2,k)+d9(3)*rhow_n(i,j+3,k)+d9(4)*rhow_n(i,j+4,k) &
                  + d9(1)*rhow_n(i,j-1,k)+d9(2)*rhow_n(i,j-2,k)+d9(3)*rhow_n(i,j-3,k)+d9(4)*rhow_n(i,j-4,k)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d9(0)*rhoe_n(i,j,k)   &
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

    ! BC at jmax
    ! ----------
    if (is_boundary(2,2)) then
       ! 8-th order centered filter
       j=ny-4
       do k=1,nz
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  +d9(0)*rho_n(i,j,k)   &
                  + d9(1)*rho_n(i,j+1,k) +d9(2)*rho_n(i,j+2,k)+d9(3)*rho_n(i,j+3,k)+d9(4)*rho_n(i,j+4,k) &
                  + d9(1)*rho_n(i,j-1,k) +d9(2)*rho_n(i,j-2,k)+d9(3)*rho_n(i,j-3,k)+d9(4)*rho_n(i,j-4,k)
             Krhou(i,j,k) = Krhou(i,j,k) +d9(0)*rhou_n(i,j,k)   &
                  + d9(1)*rhou_n(i,j+1,k)+d9(2)*rhou_n(i,j+2,k)+d9(3)*rhou_n(i,j+3,k)+d9(4)*rhou_n(i,j+4,k) &
                  + d9(1)*rhou_n(i,j-1,k)+d9(2)*rhou_n(i,j-2,k)+d9(3)*rhou_n(i,j-3,k)+d9(4)*rhou_n(i,j-4,k)
             Krhov(i,j,k) = Krhov(i,j,k) +d9(0)*rhov_n(i,j,k)   &
                  + d9(1)*rhov_n(i,j+1,k)+d9(2)*rhov_n(i,j+2,k)+d9(3)*rhov_n(i,j+3,k)+d9(4)*rhov_n(i,j+4,k) &
                  + d9(1)*rhov_n(i,j-1,k)+d9(2)*rhov_n(i,j-2,k)+d9(3)*rhov_n(i,j-3,k)+d9(4)*rhov_n(i,j-4,k)
             Krhow(i,j,k) = Krhow(i,j,k) +d9(0)*rhow_n(i,j,k)   &
                  + d9(1)*rhow_n(i,j+1,k)+d9(2)*rhow_n(i,j+2,k)+d9(3)*rhow_n(i,j+3,k)+d9(4)*rhow_n(i,j+4,k) &
                  + d9(1)*rhow_n(i,j-1,k)+d9(2)*rhow_n(i,j-2,k)+d9(3)*rhow_n(i,j-3,k)+d9(4)*rhow_n(i,j-4,k)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d9(0)*rhoe_n(i,j,k)   &
                  + d9(1)*rhoe_n(i,j+1,k)+d9(2)*rhoe_n(i,j+2,k)+d9(3)*rhoe_n(i,j+3,k)+d9(4)*rhoe_n(i,j+4,k) &
                  + d9(1)*rhoe_n(i,j-1,k)+d9(2)*rhoe_n(i,j-2,k)+d9(3)*rhoe_n(i,j-3,k)+d9(4)*rhoe_n(i,j-4,k)
          enddo
       enddo
       ! 6-th order centered filter
       j=ny-3
       do k=1,nz
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  +d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i,j+1,k) +d7(2)*rho_n(i,j+2,k)+d7(3)*rho_n(i,j+3,k) &
                  + d7(1)*rho_n(i,j-1,k) +d7(2)*rho_n(i,j-2,k)+d7(3)*rho_n(i,j-3,k)
             Krhou(i,j,k) = Krhou(i,j,k) +d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i,j+1,k)+d7(2)*rhou_n(i,j+2,k)+d7(3)*rhou_n(i,j+3,k) &
                  + d7(1)*rhou_n(i,j-1,k)+d7(2)*rhou_n(i,j-2,k)+d7(3)*rhou_n(i,j-3,k)
             Krhov(i,j,k) = Krhov(i,j,k) +d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i,j+1,k)+d7(2)*rhov_n(i,j+2,k)+d7(3)*rhov_n(i,j+3,k) &
                  + d7(1)*rhov_n(i,j-1,k)+d7(2)*rhov_n(i,j-2,k)+d7(3)*rhov_n(i,j-3,k)
             Krhow(i,j,k) = Krhow(i,j,k) +d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i,j+1,k)+d7(2)*rhow_n(i,j+2,k)+d7(3)*rhow_n(i,j+3,k) &
                  + d7(1)*rhow_n(i,j-1,k)+d7(2)*rhow_n(i,j-2,k)+d7(3)*rhow_n(i,j-3,k)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i,j+1,k)+d7(2)*rhoe_n(i,j+2,k)+d7(3)*rhoe_n(i,j+3,k) &
                  + d7(1)*rhoe_n(i,j-1,k)+d7(2)*rhoe_n(i,j-2,k)+d7(3)*rhoe_n(i,j-3,k)
          enddo
       enddo
       ! 4-th order centered filter
       j=ny-2
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
                  + d5(1)*rho_n(i,j,k+1) +d5(2)*rho_n(i,j,k+2) &
                  + d5(1)*rho_n(i,j,k-1) +d5(2)*rho_n(i,j,k-2)
             Krhou(i,j,k) = Krhou(i,j,k) +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j,k+1)+d5(2)*rhou_n(i,j,k+2) &
                  + d5(1)*rhou_n(i,j,k-1)+d5(2)*rhou_n(i,j,k-2)
             Krhov(i,j,k) = Krhov(i,j,k) +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j,k+1)+d5(2)*rhov_n(i,j,k+2) &
                  + d5(1)*rhov_n(i,j,k-1)+d5(2)*rhov_n(i,j,k-2)
             Krhow(i,j,k) = Krhow(i,j,k) +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j,k+1)+d5(2)*rhow_n(i,j,k+2) &
                  + d5(1)*rhow_n(i,j,k-1)+d5(2)*rhow_n(i,j,k-2)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i,j,k+1)+d5(2)*rhoe_n(i,j,k+2) &
                  + d5(1)*rhoe_n(i,j,k-1)+d5(2)*rhoe_n(i,j,k-2)
          enddo
       enddo
       ! 6-th order centered filter
       k=4
       do j=1,ny
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  +d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i,j,k+1) +d7(2)*rho_n(i,j,k+2)+d7(3)*rho_n(i,j,k+3) &
                  + d7(1)*rho_n(i,j,k-1) +d7(2)*rho_n(i,j,k-2)+d7(3)*rho_n(i,j,k-3)
             Krhou(i,j,k) = Krhou(i,j,k) +d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i,j,k+1)+d7(2)*rhou_n(i,j,k+2)+d7(3)*rhou_n(i,j,k+3) &
                  + d7(1)*rhou_n(i,j,k-1)+d7(2)*rhou_n(i,j,k-2)+d7(3)*rhou_n(i,j,k-3)
             Krhov(i,j,k) = Krhov(i,j,k) +d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i,j,k+1)+d7(2)*rhov_n(i,j,k+2)+d7(3)*rhov_n(i,j,k+3) &
                  + d7(1)*rhov_n(i,j,k-1)+d7(2)*rhov_n(i,j,k-2)+d7(3)*rhov_n(i,j,k-3)
             Krhow(i,j,k) = Krhow(i,j,k) +d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i,j,k+1)+d7(2)*rhow_n(i,j,k+2)+d7(3)*rhow_n(i,j,k+3) &
                  + d7(1)*rhow_n(i,j,k-1)+d7(2)*rhow_n(i,j,k-2)+d7(3)*rhow_n(i,j,k-3)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i,j,k+1)+d7(2)*rhoe_n(i,j,k+2)+d7(3)*rhoe_n(i,j,k+3) &
                  + d7(1)*rhoe_n(i,j,k-1)+d7(2)*rhoe_n(i,j,k-2)+d7(3)*rhoe_n(i,j,k-3)
          enddo
       enddo
       ! 8-th order centered filter
       k=5
       do j=1,ny
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  +d9(0)*rho_n(i,j,k)   &
                  + d9(1)*rho_n(i,j,k+1) +d9(2)*rho_n(i,j,k+2)+d9(3)*rho_n(i,j,k+3)+d9(4)*rho_n(i,j,k+4) &
                  + d9(1)*rho_n(i,j,k-1) +d9(2)*rho_n(i,j,k-2)+d9(3)*rho_n(i,j,k-3)+d9(4)*rho_n(i,j,k-4)
             Krhou(i,j,k) = Krhou(i,j,k) +d9(0)*rhou_n(i,j,k)   &
                  + d9(1)*rhou_n(i,j,k+1)+d9(2)*rhou_n(i,j,k+2)+d9(3)*rhou_n(i,j,k+3)+d9(4)*rhou_n(i,j,k+4) &
                  + d9(1)*rhou_n(i,j,k-1)+d9(2)*rhou_n(i,j,k-2)+d9(3)*rhou_n(i,j,k-3)+d9(4)*rhou_n(i,j,k-4)
             Krhov(i,j,k) = Krhov(i,j,k) +d9(0)*rhov_n(i,j,k)   &
                  + d9(1)*rhov_n(i,j,k+1)+d9(2)*rhov_n(i,j,k+2)+d9(3)*rhov_n(i,j,k+3)+d9(4)*rhov_n(i,j,k+4) &
                  + d9(1)*rhov_n(i,j,k-1)+d9(2)*rhov_n(i,j,k-2)+d9(3)*rhov_n(i,j,k-3)+d9(4)*rhov_n(i,j,k-4)
             Krhow(i,j,k) = Krhow(i,j,k) +d9(0)*rhow_n(i,j,k)   &
                  + d9(1)*rhow_n(i,j,k+1)+d9(2)*rhow_n(i,j,k+2)+d9(3)*rhow_n(i,j,k+3)+d9(4)*rhow_n(i,j,k+4) &
                  + d9(1)*rhow_n(i,j,k-1)+d9(2)*rhow_n(i,j,k-2)+d9(3)*rhow_n(i,j,k-3)+d9(4)*rhow_n(i,j,k-4)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d9(0)*rhoe_n(i,j,k)   &
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

    ! BC at kmax
    ! ----------
    if (is_boundary(3,2)) then
       ! 8-th order centered filter
       k=nz-4
       do j=1,ny
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  +d9(0)*rho_n(i,j,k)   &
                  + d9(1)*rho_n(i,j,k+1) +d9(2)*rho_n(i,j,k+2)+d9(3)*rho_n(i,j,k+3)+d9(4)*rho_n(i,j,k+4) &
                  + d9(1)*rho_n(i,j,k-1) +d9(2)*rho_n(i,j,k-2)+d9(3)*rho_n(i,j,k-3)+d9(4)*rho_n(i,j,k-4)
             Krhou(i,j,k) = Krhou(i,j,k) +d9(0)*rhou_n(i,j,k)   &
                  + d9(1)*rhou_n(i,j,k+1)+d9(2)*rhou_n(i,j,k+2)+d9(3)*rhou_n(i,j,k+3)+d9(4)*rhou_n(i,j,k+4) &
                  + d9(1)*rhou_n(i,j,k-1)+d9(2)*rhou_n(i,j,k-2)+d9(3)*rhou_n(i,j,k-3)+d9(4)*rhou_n(i,j,k-4)
             Krhov(i,j,k) = Krhov(i,j,k) +d9(0)*rhov_n(i,j,k)   &
                  + d9(1)*rhov_n(i,j,k+1)+d9(2)*rhov_n(i,j,k+2)+d9(3)*rhov_n(i,j,k+3)+d9(4)*rhov_n(i,j,k+4) &
                  + d9(1)*rhov_n(i,j,k-1)+d9(2)*rhov_n(i,j,k-2)+d9(3)*rhov_n(i,j,k-3)+d9(4)*rhov_n(i,j,k-4)
             Krhow(i,j,k) = Krhow(i,j,k) +d9(0)*rhow_n(i,j,k)   &
                  + d9(1)*rhow_n(i,j,k+1)+d9(2)*rhow_n(i,j,k+2)+d9(3)*rhow_n(i,j,k+3)+d9(4)*rhow_n(i,j,k+4) &
                  + d9(1)*rhow_n(i,j,k-1)+d9(2)*rhow_n(i,j,k-2)+d9(3)*rhow_n(i,j,k-3)+d9(4)*rhow_n(i,j,k-4)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d9(0)*rhoe_n(i,j,k)   &
                  + d9(1)*rhoe_n(i,j,k+1)+d9(2)*rhoe_n(i,j,k+2)+d9(3)*rhoe_n(i,j,k+3)+d9(4)*rhoe_n(i,j,k+4) &
                  + d9(1)*rhoe_n(i,j,k-1)+d9(2)*rhoe_n(i,j,k-2)+d9(3)*rhoe_n(i,j,k-3)+d9(4)*rhoe_n(i,j,k-4)
          enddo
       enddo
       ! 6-th order centered filter
       k=nz-3
       do j=1,ny
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  +d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i,j,k+1) +d7(2)*rho_n(i,j,k+2)+d7(3)*rho_n(i,j,k+3) &
                  + d7(1)*rho_n(i,j,k-1) +d7(2)*rho_n(i,j,k-2)+d7(3)*rho_n(i,j,k-3)
             Krhou(i,j,k) = Krhou(i,j,k) +d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i,j,k+1)+d7(2)*rhou_n(i,j,k+2)+d7(3)*rhou_n(i,j,k+3) &
                  + d7(1)*rhou_n(i,j,k-1)+d7(2)*rhou_n(i,j,k-2)+d7(3)*rhou_n(i,j,k-3)
             Krhov(i,j,k) = Krhov(i,j,k) +d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i,j,k+1)+d7(2)*rhov_n(i,j,k+2)+d7(3)*rhov_n(i,j,k+3) &
                  + d7(1)*rhov_n(i,j,k-1)+d7(2)*rhov_n(i,j,k-2)+d7(3)*rhov_n(i,j,k-3)
             Krhow(i,j,k) = Krhow(i,j,k) +d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i,j,k+1)+d7(2)*rhow_n(i,j,k+2)+d7(3)*rhow_n(i,j,k+3) &
                  + d7(1)*rhow_n(i,j,k-1)+d7(2)*rhow_n(i,j,k-2)+d7(3)*rhow_n(i,j,k-3)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i,j,k+1)+d7(2)*rhoe_n(i,j,k+2)+d7(3)*rhoe_n(i,j,k+3) &
                  + d7(1)*rhoe_n(i,j,k-1)+d7(2)*rhoe_n(i,j,k-2)+d7(3)*rhoe_n(i,j,k-3)
          enddo
       enddo
       ! 4-th order centered filter
       k=nz-2
       do j=1,ny
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i,j,k+1) +d5(2)*rho_n(i,j,k+2) &
                  + d5(1)*rho_n(i,j,k-1) +d5(2)*rho_n(i,j,k-2)
             Krhou(i,j,k) = Krhou(i,j,k) +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j,k+1)+d5(2)*rhou_n(i,j,k+2) &
                  + d5(1)*rhou_n(i,j,k-1)+d5(2)*rhou_n(i,j,k-2)
             Krhov(i,j,k) = Krhov(i,j,k) +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j,k+1)+d5(2)*rhov_n(i,j,k+2) &
                  + d5(1)*rhov_n(i,j,k-1)+d5(2)*rhov_n(i,j,k-2)
             Krhow(i,j,k) = Krhow(i,j,k) +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j,k+1)+d5(2)*rhow_n(i,j,k+2) &
                  + d5(1)*rhow_n(i,j,k-1)+d5(2)*rhow_n(i,j,k-2)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i,j,k+1)+d5(2)*rhoe_n(i,j,k+2) &
                  + d5(1)*rhoe_n(i,j,k-1)+d5(2)*rhoe_n(i,j,k-2)
          enddo
       enddo
    endif

  end subroutine filtering_11pts_inc

  !===============================================================================
  subroutine filtering_9pts_inc
  !===============================================================================
    !> Apply filtering on 9-point stencils (+ boundaries)
    !> Added to increments (for IRS smoothing)
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    ! ---------------------------------------------------------------------------

    ! Filtering along i-direction
    ! ===========================
    
    ! BC at imin
    ! ----------
    if (is_boundary(1,1)) then
       ! 4-th order centered filter
       i=3
       do k=1,nz
          do j=1,ny
             Krho(i,j,k)  = Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i+1,j,k) +d5(2)*rho_n(i+2,j,k) &
                  + d5(1)*rho_n(i-1,j,k) +d5(2)*rho_n(i-2,j,k)
             Krhou(i,j,k) = Krhou(i,j,k) +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i+1,j,k)+d5(2)*rhou_n(i+2,j,k) &
                  + d5(1)*rhou_n(i-1,j,k)+d5(2)*rhou_n(i-2,j,k)
             Krhov(i,j,k) = Krhov(i,j,k) +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i+1,j,k)+d5(2)*rhov_n(i+2,j,k) &
                  + d5(1)*rhov_n(i-1,j,k)+d5(2)*rhov_n(i-2,j,k)
             Krhow(i,j,k) = Krhow(i,j,k) +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i+1,j,k)+d5(2)*rhow_n(i+2,j,k) &
                  + d5(1)*rhow_n(i-1,j,k)+d5(2)*rhow_n(i-2,j,k)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i+1,j,k)+d5(2)*rhoe_n(i+2,j,k) &
                  + d5(1)*rhoe_n(i-1,j,k)+d5(2)*rhoe_n(i-2,j,k)
          enddo
       enddo
       ! 6-th order centered filter
       i=4
       do k=1,nz
          do j=1,ny
             Krho(i,j,k)  = Krho(i,j,k)  +d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i+1,j,k) +d7(2)*rho_n(i+2,j,k)+d7(3)*rho_n(i+3,j,k) &
                  + d7(1)*rho_n(i-1,j,k) +d7(2)*rho_n(i-2,j,k)+d7(3)*rho_n(i-3,j,k)
             Krhou(i,j,k) = Krhou(i,j,k) +d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i+1,j,k)+d7(2)*rhou_n(i+2,j,k)+d7(3)*rhou_n(i+3,j,k) &
                  + d7(1)*rhou_n(i-1,j,k)+d7(2)*rhou_n(i-2,j,k)+d7(3)*rhou_n(i-3,j,k)
             Krhov(i,j,k) = Krhov(i,j,k) +d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i+1,j,k)+d7(2)*rhov_n(i+2,j,k)+d7(3)*rhov_n(i+3,j,k) &
                  + d7(1)*rhov_n(i-1,j,k)+d7(2)*rhov_n(i-2,j,k)+d7(3)*rhov_n(i-3,j,k)
             Krhow(i,j,k) = Krhow(i,j,k) +d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i+1,j,k)+d7(2)*rhow_n(i+2,j,k)+d7(3)*rhow_n(i+3,j,k) &
                  + d7(1)*rhow_n(i-1,j,k)+d7(2)*rhow_n(i-2,j,k)+d7(3)*rhow_n(i-3,j,k)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d7(0)*rhoe_n(i,j,k)   &
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
             Krho(i,j,k)  = Krho(i,j,k)  + d9(0)*rho_n(i,j,k) &
                  + d9(1)*( rho_n(i+1,j,k)+rho_n(i-1,j,k) )&
                  + d9(2)*( rho_n(i+2,j,k)+rho_n(i-2,j,k) )&
                  + d9(3)*( rho_n(i+3,j,k)+rho_n(i-3,j,k) )&
                  + d9(4)*( rho_n(i+4,j,k)+rho_n(i-4,j,k) )

             Krhou(i,j,k) = Krhou(i,j,k)  + d9(0)*rhou_n(i,j,k) &
                  + d9(1)*( rhou_n(i+1,j,k)+rhou_n(i-1,j,k) )&
                  + d9(2)*( rhou_n(i+2,j,k)+rhou_n(i-2,j,k) )&
                  + d9(3)*( rhou_n(i+3,j,k)+rhou_n(i-3,j,k) )&
                  + d9(4)*( rhou_n(i+4,j,k)+rhou_n(i-4,j,k) )

             Krhov(i,j,k) = Krhov(i,j,k)  + d9(0)*rhov_n(i,j,k) &
                  + d9(1)*( rhov_n(i+1,j,k)+rhov_n(i-1,j,k) )&
                  + d9(2)*( rhov_n(i+2,j,k)+rhov_n(i-2,j,k) )&
                  + d9(3)*( rhov_n(i+3,j,k)+rhov_n(i-3,j,k) )&
                  + d9(4)*( rhov_n(i+4,j,k)+rhov_n(i-4,j,k) )

             Krhow(i,j,k) = Krhow(i,j,k)  + d9(0)*rhow_n(i,j,k) &
                  + d9(1)*( rhow_n(i+1,j,k)+rhow_n(i-1,j,k) )&
                  + d9(2)*( rhow_n(i+2,j,k)+rhow_n(i-2,j,k) )&
                  + d9(3)*( rhow_n(i+3,j,k)+rhow_n(i-3,j,k) )&
                  + d9(4)*( rhow_n(i+4,j,k)+rhow_n(i-4,j,k) )

             Krhoe(i,j,k) = Krhoe(i,j,k)  + d9(0)*rhoe_n(i,j,k) &
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
             Krho(i,j,k)  = Krho(i,j,k)  +d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i+1,j,k) +d7(2)*rho_n(i+2,j,k)+d7(3)*rho_n(i+3,j,k) &
                  + d7(1)*rho_n(i-1,j,k) +d7(2)*rho_n(i-2,j,k)+d7(3)*rho_n(i-3,j,k)
             Krhou(i,j,k) = Krhou(i,j,k) +d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i+1,j,k)+d7(2)*rhou_n(i+2,j,k)+d7(3)*rhou_n(i+3,j,k) &
                  + d7(1)*rhou_n(i-1,j,k)+d7(2)*rhou_n(i-2,j,k)+d7(3)*rhou_n(i-3,j,k)
             Krhov(i,j,k) = Krhov(i,j,k) +d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i+1,j,k)+d7(2)*rhov_n(i+2,j,k)+d7(3)*rhov_n(i+3,j,k) &
                  + d7(1)*rhov_n(i-1,j,k)+d7(2)*rhov_n(i-2,j,k)+d7(3)*rhov_n(i-3,j,k)
             Krhow(i,j,k) = Krhow(i,j,k) +d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i+1,j,k)+d7(2)*rhow_n(i+2,j,k)+d7(3)*rhow_n(i+3,j,k) &
                  + d7(1)*rhow_n(i-1,j,k)+d7(2)*rhow_n(i-2,j,k)+d7(3)*rhow_n(i-3,j,k)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i+1,j,k)+d7(2)*rhoe_n(i+2,j,k)+d7(3)*rhoe_n(i+3,j,k) &
                  + d7(1)*rhoe_n(i-1,j,k)+d7(2)*rhoe_n(i-2,j,k)+d7(3)*rhoe_n(i-3,j,k)
          enddo
       enddo
       ! 4-th order centered filter
       i=nx-2
       do k=1,nz
          do j=1,ny
             Krho(i,j,k)  = Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i+1,j,k) +d5(2)*rho_n(i+2,j,k) &
                  + d5(1)*rho_n(i-1,j,k) +d5(2)*rho_n(i-2,j,k)
             Krhou(i,j,k) = Krhou(i,j,k) +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i+1,j,k)+d5(2)*rhou_n(i+2,j,k) &
                  + d5(1)*rhou_n(i-1,j,k)+d5(2)*rhou_n(i-2,j,k)
             Krhov(i,j,k) = Krhov(i,j,k) +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i+1,j,k)+d5(2)*rhov_n(i+2,j,k) &
                  + d5(1)*rhov_n(i-1,j,k)+d5(2)*rhov_n(i-2,j,k)
             Krhow(i,j,k) = Krhow(i,j,k) +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i+1,j,k)+d5(2)*rhow_n(i+2,j,k) &
                  + d5(1)*rhow_n(i-1,j,k)+d5(2)*rhow_n(i-2,j,k)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i+1,j,k)+d5(2)*rhoe_n(i+2,j,k) &
                  + d5(1)*rhoe_n(i-1,j,k)+d5(2)*rhoe_n(i-2,j,k)
          enddo
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
                  + d7(1)*rho_n(i,j+1,k) +d7(2)*rho_n(i,j+2,k)+d7(3)*rho_n(i,j+3,k) &
                  + d7(1)*rho_n(i,j-1,k) +d7(2)*rho_n(i,j-2,k)+d7(3)*rho_n(i,j-3,k)
             Krhou(i,j,k) = Krhou(i,j,k) +d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i,j+1,k)+d7(2)*rhou_n(i,j+2,k)+d7(3)*rhou_n(i,j+3,k) &
                  + d7(1)*rhou_n(i,j-1,k)+d7(2)*rhou_n(i,j-2,k)+d7(3)*rhou_n(i,j-3,k)
             Krhov(i,j,k) = Krhov(i,j,k) +d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i,j+1,k)+d7(2)*rhov_n(i,j+2,k)+d7(3)*rhov_n(i,j+3,k) &
                  + d7(1)*rhov_n(i,j-1,k)+d7(2)*rhov_n(i,j-2,k)+d7(3)*rhov_n(i,j-3,k)
             Krhow(i,j,k) = Krhow(i,j,k) +d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i,j+1,k)+d7(2)*rhow_n(i,j+2,k)+d7(3)*rhow_n(i,j+3,k) &
                  + d7(1)*rhow_n(i,j-1,k)+d7(2)*rhow_n(i,j-2,k)+d7(3)*rhow_n(i,j-3,k)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d7(0)*rhoe_n(i,j,k)   &
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
                  + d7(1)*rho_n(i,j+1,k) +d7(2)*rho_n(i,j+2,k)+d7(3)*rho_n(i,j+3,k) &
                  + d7(1)*rho_n(i,j-1,k) +d7(2)*rho_n(i,j-2,k)+d7(3)*rho_n(i,j-3,k)
             Krhou(i,j,k) = Krhou(i,j,k) +d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i,j+1,k)+d7(2)*rhou_n(i,j+2,k)+d7(3)*rhou_n(i,j+3,k) &
                  + d7(1)*rhou_n(i,j-1,k)+d7(2)*rhou_n(i,j-2,k)+d7(3)*rhou_n(i,j-3,k)
             Krhov(i,j,k) = Krhov(i,j,k) +d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i,j+1,k)+d7(2)*rhov_n(i,j+2,k)+d7(3)*rhov_n(i,j+3,k) &
                  + d7(1)*rhov_n(i,j-1,k)+d7(2)*rhov_n(i,j-2,k)+d7(3)*rhov_n(i,j-3,k)
             Krhow(i,j,k) = Krhow(i,j,k) +d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i,j+1,k)+d7(2)*rhow_n(i,j+2,k)+d7(3)*rhow_n(i,j+3,k) &
                  + d7(1)*rhow_n(i,j-1,k)+d7(2)*rhow_n(i,j-2,k)+d7(3)*rhow_n(i,j-3,k)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i,j+1,k)+d7(2)*rhoe_n(i,j+2,k)+d7(3)*rhoe_n(i,j+3,k) &
                  + d7(1)*rhoe_n(i,j-1,k)+d7(2)*rhoe_n(i,j-2,k)+d7(3)*rhoe_n(i,j-3,k)
          enddo
       enddo
       ! 4-th order centered filter
       j=ny-2
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
                  + d5(1)*rho_n(i,j,k+1) +d5(2)*rho_n(i,j,k+2) &
                  + d5(1)*rho_n(i,j,k-1) +d5(2)*rho_n(i,j,k-2)
             Krhou(i,j,k) = Krhou(i,j,k) +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j,k+1)+d5(2)*rhou_n(i,j,k+2) &
                  + d5(1)*rhou_n(i,j,k-1)+d5(2)*rhou_n(i,j,k-2)
             Krhov(i,j,k) = Krhov(i,j,k) +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j,k+1)+d5(2)*rhov_n(i,j,k+2) &
                  + d5(1)*rhov_n(i,j,k-1)+d5(2)*rhov_n(i,j,k-2)
             Krhow(i,j,k) = Krhow(i,j,k) +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j,k+1)+d5(2)*rhow_n(i,j,k+2) &
                  + d5(1)*rhow_n(i,j,k-1)+d5(2)*rhow_n(i,j,k-2)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i,j,k+1)+d5(2)*rhoe_n(i,j,k+2) &
                  + d5(1)*rhoe_n(i,j,k-1)+d5(2)*rhoe_n(i,j,k-2)
          enddo
       enddo
       ! 6-th order centered filter
       k=4
       do j=1,ny
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  +d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i,j,k+1) +d7(2)*rho_n(i,j,k+2)+d7(3)*rho_n(i,j,k+3) &
                  + d7(1)*rho_n(i,j,k-1) +d7(2)*rho_n(i,j,k-2)+d7(3)*rho_n(i,j,k-3)
             Krhou(i,j,k) = Krhou(i,j,k) +d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i,j,k+1)+d7(2)*rhou_n(i,j,k+2)+d7(3)*rhou_n(i,j,k+3) &
                  + d7(1)*rhou_n(i,j,k-1)+d7(2)*rhou_n(i,j,k-2)+d7(3)*rhou_n(i,j,k-3)
             Krhov(i,j,k) = Krhov(i,j,k) +d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i,j,k+1)+d7(2)*rhov_n(i,j,k+2)+d7(3)*rhov_n(i,j,k+3) &
                  + d7(1)*rhov_n(i,j,k-1)+d7(2)*rhov_n(i,j,k-2)+d7(3)*rhov_n(i,j,k-3)
             Krhow(i,j,k) = Krhow(i,j,k) +d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i,j,k+1)+d7(2)*rhow_n(i,j,k+2)+d7(3)*rhow_n(i,j,k+3) &
                  + d7(1)*rhow_n(i,j,k-1)+d7(2)*rhow_n(i,j,k-2)+d7(3)*rhow_n(i,j,k-3)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d7(0)*rhoe_n(i,j,k)   &
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
                  + d7(1)*rho_n(i,j,k+1) +d7(2)*rho_n(i,j,k+2)+d7(3)*rho_n(i,j,k+3) &
                  + d7(1)*rho_n(i,j,k-1) +d7(2)*rho_n(i,j,k-2)+d7(3)*rho_n(i,j,k-3)
             Krhou(i,j,k) = Krhou(i,j,k) +d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i,j,k+1)+d7(2)*rhou_n(i,j,k+2)+d7(3)*rhou_n(i,j,k+3) &
                  + d7(1)*rhou_n(i,j,k-1)+d7(2)*rhou_n(i,j,k-2)+d7(3)*rhou_n(i,j,k-3)
             Krhov(i,j,k) = Krhov(i,j,k) +d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i,j,k+1)+d7(2)*rhov_n(i,j,k+2)+d7(3)*rhov_n(i,j,k+3) &
                  + d7(1)*rhov_n(i,j,k-1)+d7(2)*rhov_n(i,j,k-2)+d7(3)*rhov_n(i,j,k-3)
             Krhow(i,j,k) = Krhow(i,j,k) +d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i,j,k+1)+d7(2)*rhow_n(i,j,k+2)+d7(3)*rhow_n(i,j,k+3) &
                  + d7(1)*rhow_n(i,j,k-1)+d7(2)*rhow_n(i,j,k-2)+d7(3)*rhow_n(i,j,k-3)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i,j,k+1)+d7(2)*rhoe_n(i,j,k+2)+d7(3)*rhoe_n(i,j,k+3) &
                  + d7(1)*rhoe_n(i,j,k-1)+d7(2)*rhoe_n(i,j,k-2)+d7(3)*rhoe_n(i,j,k-3)
          enddo
       enddo
       ! 4-th order centered filter
       k=nz-2
       do j=1,ny
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i,j,k+1) +d5(2)*rho_n(i,j,k+2) &
                  + d5(1)*rho_n(i,j,k-1) +d5(2)*rho_n(i,j,k-2)
             Krhou(i,j,k) = Krhou(i,j,k) +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j,k+1)+d5(2)*rhou_n(i,j,k+2) &
                  + d5(1)*rhou_n(i,j,k-1)+d5(2)*rhou_n(i,j,k-2)
             Krhov(i,j,k) = Krhov(i,j,k) +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j,k+1)+d5(2)*rhov_n(i,j,k+2) &
                  + d5(1)*rhov_n(i,j,k-1)+d5(2)*rhov_n(i,j,k-2)
             Krhow(i,j,k) = Krhow(i,j,k) +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j,k+1)+d5(2)*rhow_n(i,j,k+2) &
                  + d5(1)*rhow_n(i,j,k-1)+d5(2)*rhow_n(i,j,k-2)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i,j,k+1)+d5(2)*rhoe_n(i,j,k+2) &
                  + d5(1)*rhoe_n(i,j,k-1)+d5(2)*rhoe_n(i,j,k-2)
          enddo
       enddo
    endif

  end subroutine filtering_9pts_inc

  !===============================================================================
  subroutine filtering_7pts_inc
  !===============================================================================
    !> Apply filtering on 7-point stencils (+ boundaries)
    !> Added to increments (for IRS smoothing)
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    ! ---------------------------------------------------------------------------

    ! Filtering along i-direction
    ! ===========================
    
    ! BC at imin
    ! ----------
    if (is_boundary(1,1)) then
       ! 4-th order centered filter
       i=3
       do k=1,nz
          do j=1,ny
             Krho(i,j,k)  = Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i+1,j,k) +d5(2)*rho_n(i+2,j,k) &
                  + d5(1)*rho_n(i-1,j,k) +d5(2)*rho_n(i-2,j,k)
             Krhou(i,j,k) = Krhou(i,j,k) +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i+1,j,k)+d5(2)*rhou_n(i+2,j,k) &
                  + d5(1)*rhou_n(i-1,j,k)+d5(2)*rhou_n(i-2,j,k)
             Krhov(i,j,k) = Krhov(i,j,k) +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i+1,j,k)+d5(2)*rhov_n(i+2,j,k) &
                  + d5(1)*rhov_n(i-1,j,k)+d5(2)*rhov_n(i-2,j,k)
             Krhow(i,j,k) = Krhow(i,j,k) +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i+1,j,k)+d5(2)*rhow_n(i+2,j,k) &
                  + d5(1)*rhow_n(i-1,j,k)+d5(2)*rhow_n(i-2,j,k)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d5(0)*rhoe_n(i,j,k)   &
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
             Krho(i,j,k)  = Krho(i,j,k)  + d7(0)*rho_n(i,j,k) &
                  + d7(1)*( rho_n(i+1,j,k)+rho_n(i-1,j,k) )&
                  + d7(2)*( rho_n(i+2,j,k)+rho_n(i-2,j,k) )&
                  + d7(3)*( rho_n(i+3,j,k)+rho_n(i-3,j,k) )

             Krhou(i,j,k) = Krhou(i,j,k)  + d7(0)*rhou_n(i,j,k) &
                  + d7(1)*( rhou_n(i+1,j,k)+rhou_n(i-1,j,k) )&
                  + d7(2)*( rhou_n(i+2,j,k)+rhou_n(i-2,j,k) )&
                  + d7(3)*( rhou_n(i+3,j,k)+rhou_n(i-3,j,k) )

             Krhov(i,j,k) = Krhov(i,j,k)  + d7(0)*rhov_n(i,j,k) &
                  + d7(1)*( rhov_n(i+1,j,k)+rhov_n(i-1,j,k) )&
                  + d7(2)*( rhov_n(i+2,j,k)+rhov_n(i-2,j,k) )&
                  + d7(3)*( rhov_n(i+3,j,k)+rhov_n(i-3,j,k) )

             Krhow(i,j,k) = Krhow(i,j,k)  + d7(0)*rhow_n(i,j,k) &
                  + d7(1)*( rhow_n(i+1,j,k)+rhow_n(i-1,j,k) )&
                  + d7(2)*( rhow_n(i+2,j,k)+rhow_n(i-2,j,k) )&
                  + d7(3)*( rhow_n(i+3,j,k)+rhow_n(i-3,j,k) )

             Krhoe(i,j,k) = Krhoe(i,j,k)  + d7(0)*rhoe_n(i,j,k) &
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
             Krho(i,j,k)  = Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i+1,j,k) +d5(2)*rho_n(i+2,j,k) &
                  + d5(1)*rho_n(i-1,j,k) +d5(2)*rho_n(i-2,j,k)
             Krhou(i,j,k) = Krhou(i,j,k) +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i+1,j,k)+d5(2)*rhou_n(i+2,j,k) &
                  + d5(1)*rhou_n(i-1,j,k)+d5(2)*rhou_n(i-2,j,k)
             Krhov(i,j,k) = Krhov(i,j,k) +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i+1,j,k)+d5(2)*rhov_n(i+2,j,k) &
                  + d5(1)*rhov_n(i-1,j,k)+d5(2)*rhov_n(i-2,j,k)
             Krhow(i,j,k) = Krhow(i,j,k) +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i+1,j,k)+d5(2)*rhow_n(i+2,j,k) &
                  + d5(1)*rhow_n(i-1,j,k)+d5(2)*rhow_n(i-2,j,k)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i+1,j,k)+d5(2)*rhoe_n(i+2,j,k) &
                  + d5(1)*rhoe_n(i-1,j,k)+d5(2)*rhoe_n(i-2,j,k)
          enddo
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
                  + d5(1)*rho_n(i,j,k+1) +d5(2)*rho_n(i,j,k+2) &
                  + d5(1)*rho_n(i,j,k-1) +d5(2)*rho_n(i,j,k-2)
             Krhou(i,j,k) = Krhou(i,j,k) +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j,k+1)+d5(2)*rhou_n(i,j,k+2) &
                  + d5(1)*rhou_n(i,j,k-1)+d5(2)*rhou_n(i,j,k-2)
             Krhov(i,j,k) = Krhov(i,j,k) +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j,k+1)+d5(2)*rhov_n(i,j,k+2) &
                  + d5(1)*rhov_n(i,j,k-1)+d5(2)*rhov_n(i,j,k-2)
             Krhow(i,j,k) = Krhow(i,j,k) +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j,k+1)+d5(2)*rhow_n(i,j,k+2) &
                  + d5(1)*rhow_n(i,j,k-1)+d5(2)*rhow_n(i,j,k-2)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d5(0)*rhoe_n(i,j,k)   &
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
                  + d5(1)*rho_n(i,j,k+1) +d5(2)*rho_n(i,j,k+2) &
                  + d5(1)*rho_n(i,j,k-1) +d5(2)*rho_n(i,j,k-2)
             Krhou(i,j,k) = Krhou(i,j,k) +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j,k+1)+d5(2)*rhou_n(i,j,k+2) &
                  + d5(1)*rhou_n(i,j,k-1)+d5(2)*rhou_n(i,j,k-2)
             Krhov(i,j,k) = Krhov(i,j,k) +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j,k+1)+d5(2)*rhov_n(i,j,k+2) &
                  + d5(1)*rhov_n(i,j,k-1)+d5(2)*rhov_n(i,j,k-2)
             Krhow(i,j,k) = Krhow(i,j,k) +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j,k+1)+d5(2)*rhow_n(i,j,k+2) &
                  + d5(1)*rhow_n(i,j,k-1)+d5(2)*rhow_n(i,j,k-2)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i,j,k+1)+d5(2)*rhoe_n(i,j,k+2) &
                  + d5(1)*rhoe_n(i,j,k-1)+d5(2)*rhoe_n(i,j,k-2)
          enddo
       enddo
    endif

  end subroutine filtering_7pts_inc

  !===============================================================================
  subroutine filtering_5pts_inc
  !===============================================================================
    !> Apply filtering on 5-point stencils (+ boundaries)
    !> Added to increments (for IRS smoothing)
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    ! ---------------------------------------------------------------------------

    ! Filtering along i-direction
    ! ===========================

    ! Interior points
    ! ---------------
    ! 4-th order centered filter
    do k=1,nz
       do j=1,ny
          do i=ndx,nfx
             Krho(i,j,k)  = Krho(i,j,k)  + d5(0)*rho_n(i,j,k) &
                  + d5(1)*( rho_n(i+1,j,k)+rho_n(i-1,j,k) )&
                  + d5(2)*( rho_n(i+2,j,k)+rho_n(i-2,j,k) )

             Krhou(i,j,k) = Krhou(i,j,k)  + d5(0)*rhou_n(i,j,k) &
                  + d5(1)*( rhou_n(i+1,j,k)+rhou_n(i-1,j,k) )&
                  + d5(2)*( rhou_n(i+2,j,k)+rhou_n(i-2,j,k) )

             Krhov(i,j,k) = Krhov(i,j,k)  + d5(0)*rhov_n(i,j,k) &
                  + d5(1)*( rhov_n(i+1,j,k)+rhov_n(i-1,j,k) )&
                  + d5(2)*( rhov_n(i+2,j,k)+rhov_n(i-2,j,k) )

             Krhow(i,j,k) = Krhow(i,j,k)  + d5(0)*rhow_n(i,j,k) &
                  + d5(1)*( rhow_n(i+1,j,k)+rhow_n(i-1,j,k) )&
                  + d5(2)*( rhow_n(i+2,j,k)+rhow_n(i-2,j,k) )

             Krhoe(i,j,k) = Krhoe(i,j,k)  + d5(0)*rhoe_n(i,j,k) &
                  + d5(1)*( rhoe_n(i+1,j,k)+rhoe_n(i-1,j,k) )&
                  + d5(2)*( rhoe_n(i+2,j,k)+rhoe_n(i-2,j,k) )
          enddo
       enddo
    enddo

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

  end subroutine filtering_5pts_inc

end module mod_filtering_inc

