!=================================================================================
submodule (mod_flux_euler) smod_flux_euler_wedge_9pts_c
!=================================================================================
  !> author: XG
  !> date: January 2024
  !> 2.5D curvilinear version - compute Eulerian fluxes (inviscid part) with 9-pt stencil
  !> Add_on routine to treat double corners  (wedge)
!=================================================================================

contains

  !===============================================================================
  module subroutine flux_euler_w_imin_jmin_9pts_c
  !===============================================================================
    !> Derivatives of Euler fluxes for edge I_MIN/J_MIN close to two walls
    !> - 2.5D curvilinear coordinate - 9-point stencil -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: vc,rovc
    ! ----------------------------------------------------------------------------

    ! ~> First pass: 9-pt metrics & order reduction for derivatives
    ! =============================================================

    ! fluxes along ksi
    ! ----------------
    do k=ndz_e,nfz_e
       do j=1,4
          do i=-3,8
             ! contravariant velocity
             vc=uu(i,j,k)*y_eta(i,j)-vv(i,j,k)*x_eta(i,j)
             rovc=rho_n(i,j,k)*vc

             Frho(i,j,k)  = rovc
             Frhou(i,j,k) = rovc*uu(i,j,k)+prs(i,j,k)*y_eta(i,j)
             Frhov(i,j,k) = rovc*vv(i,j,k)-prs(i,j,k)*x_eta(i,j)
             Frhow(i,j,k) = rovc*ww(i,j,k)
             Frhoe(i,j,k) = (rhoe_n(i,j,k)+prs(i,j,k))*vc
          enddo
       enddo
    enddo

    ! derivatives along ksi
    ! ---------------------
    i=1
    do j=1,4
       do k=ndz_e,nfz_e
          Krho(i,j,k)  = ( Frho(i+1,j,k)- Frho(i,j,k))
          Krhou(i,j,k) = (Frhou(i+1,j,k)-Frhou(i,j,k))
          Krhov(i,j,k) = (Frhov(i+1,j,k)-Frhov(i,j,k))
          Krhow(i,j,k) = (Frhow(i+1,j,k)-Frhow(i,j,k))
          Krhoe(i,j,k) = (Frhoe(i+1,j,k)-Frhoe(i,j,k))
       enddo
    enddo

    i=2
    do j=1,4
       do k=ndz_e,nfz_e
          Krho(i,j,k)  = 0.5_wp*( Frho(i+1,j,k)- Frho(i-1,j,k))
          Krhou(i,j,k) = 0.5_wp*(Frhou(i+1,j,k)-Frhou(i-1,j,k))
          Krhov(i,j,k) = 0.5_wp*(Frhov(i+1,j,k)-Frhov(i-1,j,k))
          Krhow(i,j,k) = 0.5_wp*(Frhow(i+1,j,k)-Frhow(i-1,j,k))
          Krhoe(i,j,k) = 0.5_wp*(Frhoe(i+1,j,k)-Frhoe(i-1,j,k))
       enddo
    enddo

    i=3
    do j=1,4
       do k=ndz_e,nfz_e
          Krho(i,j,k)  = a5(1)*( Frho(i+1,j,k)- Frho(i-1,j,k)) &
                       + a5(2)*( Frho(i+2,j,k)- Frho(i-2,j,k))
          Krhou(i,j,k) = a5(1)*(Frhou(i+1,j,k)-Frhou(i-1,j,k)) &
                       + a5(2)*(Frhou(i+2,j,k)-Frhou(i-2,j,k))
          Krhov(i,j,k) = a5(1)*(Frhov(i+1,j,k)-Frhov(i-1,j,k)) &
                       + a5(2)*(Frhov(i+2,j,k)-Frhov(i-2,j,k))
          Krhow(i,j,k) = a5(1)*(Frhow(i+1,j,k)-Frhow(i-1,j,k)) &
                       + a5(2)*(Frhow(i+2,j,k)-Frhow(i-2,j,k))
          Krhoe(i,j,k) = a5(1)*(Frhoe(i+1,j,k)-Frhoe(i-1,j,k)) &
                       + a5(2)*(Frhoe(i+2,j,k)-Frhoe(i-2,j,k))
       enddo
    enddo

    i=4
    do j=1,4
       do k=ndz_e,nfz_e
          Krho(i,j,k)  = a7(1) * ( Frho(i+1,j,k)-Frho(i-1,j,k) ) &
                       + a7(2) * ( Frho(i+2,j,k)-Frho(i-2,j,k) ) &
                       + a7(3) * ( Frho(i+3,j,k)-Frho(i-3,j,k) )
          Krhou(i,j,k) = a7(1) * ( Frhou(i+1,j,k)-Frhou(i-1,j,k) ) &
                       + a7(2) * ( Frhou(i+2,j,k)-Frhou(i-2,j,k) ) &
                       + a7(3) * ( Frhou(i+3,j,k)-Frhou(i-3,j,k) )
          Krhov(i,j,k) = a7(1) * ( Frhov(i+1,j,k)-Frhov(i-1,j,k) ) &
                       + a7(2) * ( Frhov(i+2,j,k)-Frhov(i-2,j,k) ) &
                       + a7(3) * ( Frhov(i+3,j,k)-Frhov(i-3,j,k) )
          Krhow(i,j,k) = a7(1) * ( Frhow(i+1,j,k)-Frhow(i-1,j,k) ) &
                       + a7(2) * ( Frhow(i+2,j,k)-Frhow(i-2,j,k) ) &
                       + a7(3) * ( Frhow(i+3,j,k)-Frhow(i-3,j,k) )
          Krhoe(i,j,k) = a7(1) * ( Frhoe(i+1,j,k)-Frhoe(i-1,j,k) ) &
                       + a7(2) * ( Frhoe(i+2,j,k)-Frhoe(i-2,j,k) ) &
                       + a7(3) * ( Frhoe(i+3,j,k)-Frhoe(i-3,j,k) )
       enddo
    enddo

    ! modified fluxes along eta
    ! -------------------------
    do k=ndz_e,nfz_e
       do i=1,4
          do j=-3,8
             ! contravariant velocity
             vc=vv(i,j,k)*x_ksi(i,j)-uu(i,j,k)*y_ksi(i,j)
             rovc=rho_n(i,j,k)*vc

             Frho(i,j,k)  =  rovc
             Frhou(i,j,k) =  rovc*uu(i,j,k)-prs(i,j,k)*y_ksi(i,j)
             Frhov(i,j,k) =  rovc*vv(i,j,k)+prs(i,j,k)*x_ksi(i,j)
             Frhow(i,j,k) =  rovc*ww(i,j,k)
             Frhoe(i,j,k) = (rhoe_n(i,j,k)+prs(i,j,k))*vc
          enddo
       enddo
    enddo

    ! derivatives along eta
    ! ---------------------
    j=1
    do i=1,4
       do k=ndz_e,nfz_e  
          Krho(i,j,k)  = (Frho(i,j+1,k) -Frho(i,j,k)) + Krho(i,j,k)
          Krhou(i,j,k) = (Frhou(i,j+1,k)-Frhou(i,j,k)) + Krhou(i,j,k)
          Krhov(i,j,k) = (Frhov(i,j+1,k)-Frhov(i,j,k)) + Krhov(i,j,k)
          Krhow(i,j,k) = (Frhow(i,j+1,k)-Frhow(i,j,k)) + Krhow(i,j,k)
          Krhoe(i,j,k) = (Frhoe(i,j+1,k)-Frhoe(i,j,k)) + Krhoe(i,j,k)
       enddo
    enddo

    j=2
    do i=1,4
       do k=ndz_e,nfz_e  
          Krho(i,j,k)  = 0.5_wp*(Frho(i,j+1,k) -Frho(i,j-1,k)) + Krho(i,j,k)
          Krhou(i,j,k) = 0.5_wp*(Frhou(i,j+1,k)-Frhou(i,j-1,k)) + Krhou(i,j,k)
          Krhov(i,j,k) = 0.5_wp*(Frhov(i,j+1,k)-Frhov(i,j-1,k)) + Krhov(i,j,k)
          Krhow(i,j,k) = 0.5_wp*(Frhow(i,j+1,k)-Frhow(i,j-1,k)) + Krhow(i,j,k)
          Krhoe(i,j,k) = 0.5_wp*(Frhoe(i,j+1,k)-Frhoe(i,j-1,k)) + Krhoe(i,j,k)
       enddo
    enddo

    j=3
    do i=1,4
       do k=ndz_e,nfz_e
          Krho(i,j,k)  = a5(1)*(Frho(i,j+1,k)-Frho(i,j-1,k)) &
                       + a5(2)*(Frho(i,j+2,k)-Frho(i,j-2,k)) + Krho(i,j,k)
          Krhou(i,j,k) = a5(1)*(Frhou(i,j+1,k)-Frhou(i,j-1,k)) &
                       + a5(2)*(Frhou(i,j+2,k)-Frhou(i,j-2,k)) + Krhou(i,j,k)
          Krhov(i,j,k) = a5(1)*(Frhov(i,j+1,k)-Frhov(i,j-1,k)) &
                       + a5(2)*(Frhov(i,j+2,k)-Frhov(i,j-2,k)) + Krhov(i,j,k)
          Krhow(i,j,k) = a5(1)*(Frhow(i,j+1,k)-Frhow(i,j-1,k)) &
                       + a5(2)*(Frhow(i,j+2,k)-Frhow(i,j-2,k)) + Krhow(i,j,k)
          Krhoe(i,j,k) = a5(1)*(Frhoe(i,j+1,k)-Frhoe(i,j-1,k)) &
                       + a5(2)*(Frhoe(i,j+2,k)-Frhoe(i,j-2,k)) + Krhoe(i,j,k)
       enddo
    enddo

    j=4
    do i=1,4
       do k=ndz_e,nfz_e
          Krho(i,j,k)  = a7(1) * ( Frho(i,j+1,k)-Frho(i,j-1,k) ) &
                       + a7(2) * ( Frho(i,j+2,k)-Frho(i,j-2,k) ) &
                       + a7(3) * ( Frho(i,j+3,k)-Frho(i,j-3,k) ) + Krho(i,j,k)
          Krhou(i,j,k) = a7(1) * ( Frhou(i,j+1,k)-Frhou(i,j-1,k) ) &
                       + a7(2) * ( Frhou(i,j+2,k)-Frhou(i,j-2,k) ) &
                       + a7(3) * ( Frhou(i,j+3,k)-Frhou(i,j-3,k) ) + Krhou(i,j,k)
          Krhov(i,j,k) = a7(1) * ( Frhov(i,j+1,k)-Frhov(i,j-1,k) ) &
                       + a7(2) * ( Frhov(i,j+2,k)-Frhov(i,j-2,k) ) &
                       + a7(3) * ( Frhov(i,j+3,k)-Frhov(i,j-3,k) ) + Krhov(i,j,k)
          Krhow(i,j,k) = a7(1) * ( Frhow(i,j+1,k)-Frhow(i,j-1,k) ) &
                       + a7(2) * ( Frhow(i,j+2,k)-Frhow(i,j-2,k) ) &
                       + a7(3) * ( Frhow(i,j+3,k)-Frhow(i,j-3,k) ) + Krhow(i,j,k)
          Krhoe(i,j,k) = a7(1) * ( Frhoe(i,j+1,k)-Frhoe(i,j-1,k) ) &
                       + a7(2) * ( Frhoe(i,j+2,k)-Frhoe(i,j-2,k) ) &
                       + a7(3) * ( Frhoe(i,j+3,k)-Frhoe(i,j-3,k) ) + Krhoe(i,j,k)
       enddo
    enddo

    ! ~> Second pass: order reduction for metrics & 9-pt derivatives
    ! ==============================================================

    ! modified fluxes along ksi
    ! -------------------------
    do k=ndz_e,nfz_e
       do j=1,4
          do i=-3,8
             ! contravariant velocity
             vc=uu(i,j,k)*y_eta_imin_jmin(i,j)-vv(i,j,k)*x_eta_imin_jmin(i,j)
             rovc=rho_n(i,j,k)*vc

             Frho(i,j,k)  = rovc
             Frhou(i,j,k) = rovc*uu(i,j,k)+prs(i,j,k)*y_eta_imin_jmin(i,j)
             Frhov(i,j,k) = rovc*vv(i,j,k)-prs(i,j,k)*x_eta_imin_jmin(i,j)
             Frhow(i,j,k) = rovc*ww(i,j,k)
             Frhoe(i,j,k) = (rhoe_n(i,j,k)+prs(i,j,k))*vc
          enddo
       enddo
    enddo

    ! derivatives along ksi
    ! ---------------------
    do k=ndz_e,nfz_e
       do j=1,4
          do i=1,4
             Krho(i,j,k)  = a9(1) * ( Frho(i+1,j,k)-Frho(i-1,j,k) )   &
                          + a9(2) * ( Frho(i+2,j,k)-Frho(i-2,j,k) )   &
                          + a9(3) * ( Frho(i+3,j,k)-Frho(i-3,j,k) )   &
                          + a9(4) * ( Frho(i+4,j,k)-Frho(i-4,j,k) ) + Krho(i,j,k)

             Krhou(i,j,k) = a9(1) * ( Frhou(i+1,j,k)-Frhou(i-1,j,k) ) &
                          + a9(2) * ( Frhou(i+2,j,k)-Frhou(i-2,j,k) ) &
                          + a9(3) * ( Frhou(i+3,j,k)-Frhou(i-3,j,k) ) &
                          + a9(4) * ( Frhou(i+4,j,k)-Frhou(i-4,j,k) ) + Krhou(i,j,k)

             Krhov(i,j,k) = a9(1) * ( Frhov(i+1,j,k)-Frhov(i-1,j,k) ) &
                          + a9(2) * ( Frhov(i+2,j,k)-Frhov(i-2,j,k) ) &
                          + a9(3) * ( Frhov(i+3,j,k)-Frhov(i-3,j,k) ) &
                          + a9(4) * ( Frhov(i+4,j,k)-Frhov(i-4,j,k) ) + Krhov(i,j,k)

             Krhow(i,j,k) = a9(1) * ( Frhow(i+1,j,k)-Frhow(i-1,j,k) ) &
                          + a9(2) * ( Frhow(i+2,j,k)-Frhow(i-2,j,k) ) &
                          + a9(3) * ( Frhow(i+3,j,k)-Frhow(i-3,j,k) ) &
                          + a9(4) * ( Frhow(i+4,j,k)-Frhow(i-4,j,k) ) + Krhow(i,j,k)

             Krhoe(i,j,k) = a9(1) * ( Frhoe(i+1,j,k)-Frhoe(i-1,j,k) ) &
                          + a9(2) * ( Frhoe(i+2,j,k)-Frhoe(i-2,j,k) ) &
                          + a9(3) * ( Frhoe(i+3,j,k)-Frhoe(i-3,j,k) ) &
                          + a9(4) * ( Frhoe(i+4,j,k)-Frhoe(i-4,j,k) ) + Krhoe(i,j,k)
          enddo
       enddo
    enddo

    ! modified fluxes along eta
    ! -------------------------
    do k=ndz_e,nfz_e
       do i=1,4
          do j=-3,8
             ! contravariant velocity
             vc=vv(i,j,k)*x_ksi_imin_jmin(i,j)-uu(i,j,k)*y_ksi_imin_jmin(i,j)
             rovc=rho_n(i,j,k)*vc

             Frho(i,j,k)  =  rovc
             Frhou(i,j,k) =  rovc*uu(i,j,k)-prs(i,j,k)*y_ksi_imin_jmin(i,j)
             Frhov(i,j,k) =  rovc*vv(i,j,k)+prs(i,j,k)*x_ksi_imin_jmin(i,j)
             Frhow(i,j,k) =  rovc*ww(i,j,k)
             Frhoe(i,j,k) = (rhoe_n(i,j,k)+prs(i,j,k))*vc
          enddo
       enddo
    enddo

    ! derivatives along eta
    ! ---------------------
    do k=ndz_e,nfz_e
       do j=1,4
          do i=1,4
             Krho(i,j,k)  = a9(1) * ( Frho(i,j+1,k)-Frho(i,j-1,k) ) &
                          + a9(2) * ( Frho(i,j+2,k)-Frho(i,j-2,k) ) &
                          + a9(3) * ( Frho(i,j+3,k)-Frho(i,j-3,k) ) &
                          + a9(4) * ( Frho(i,j+4,k)-Frho(i,j-4,k) ) &
                          + Krho(i,j,k)

             Krhou(i,j,k) = a9(1) * ( Frhou(i,j+1,k)-Frhou(i,j-1,k) ) &
                          + a9(2) * ( Frhou(i,j+2,k)-Frhou(i,j-2,k) ) &
                          + a9(3) * ( Frhou(i,j+3,k)-Frhou(i,j-3,k) ) &
                          + a9(4) * ( Frhou(i,j+4,k)-Frhou(i,j-4,k) ) &
                          + Krhou(i,j,k)

             Krhov(i,j,k) = a9(1) * ( Frhov(i,j+1,k)-Frhov(i,j-1,k) ) &
                          + a9(2) * ( Frhov(i,j+2,k)-Frhov(i,j-2,k) ) &
                          + a9(3) * ( Frhov(i,j+3,k)-Frhov(i,j-3,k) ) &
                          + a9(4) * ( Frhov(i,j+4,k)-Frhov(i,j-4,k) ) &
                          + Krhov(i,j,k)

             Krhow(i,j,k) = a9(1) * ( Frhow(i,j+1,k)-Frhow(i,j-1,k) ) &
                          + a9(2) * ( Frhow(i,j+2,k)-Frhow(i,j-2,k) ) &
                          + a9(3) * ( Frhow(i,j+3,k)-Frhow(i,j-3,k) ) &
                          + a9(4) * ( Frhow(i,j+4,k)-Frhow(i,j-4,k) ) &
                          + Krhow(i,j,k)

             Krhoe(i,j,k) = a9(1) * ( Frhoe(i,j+1,k)-Frhoe(i,j-1,k) ) &
                          + a9(2) * ( Frhoe(i,j+2,k)-Frhoe(i,j-2,k) ) &
                          + a9(3) * ( Frhoe(i,j+3,k)-Frhoe(i,j-3,k) ) &
                          + a9(4) * ( Frhoe(i,j+4,k)-Frhoe(i,j-4,k) ) &
                          + Krhoe(i,j,k)
          enddo
       enddo
    enddo

    ! Average of the two derivative evaluations
    ! =========================================
    do k=ndz_e,nfz_e
       do j=1,4
          do i=1,4
             Krho(i,j,k)  = Krho(i,j,k)*0.5
             Krhou(i,j,k) = Krhou(i,j,k)*0.5
             Krhov(i,j,k) = Krhov(i,j,k)*0.5
             Krhow(i,j,k) = Krhow(i,j,k)*0.5
             Krhoe(i,j,k) = Krhoe(i,j,k)*0.5
          enddo
       enddo
    enddo

  end subroutine flux_euler_w_imin_jmin_9pts_c

  !===============================================================================
  module subroutine flux_euler_w_imin_jmax_9pts_c
  !===============================================================================
    !> Derivatives of Euler fluxes for edge I_MIN/J_MAX close to two walls
    !> - curvilinear coordinate - 9-point stencil -
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: vc,rovc
    ! ---------------------------------------------------------------------------

    ! ~> First pass: 9-pt metrics & order reduction for derivatives
    ! =============================================================

    ! fluxes along ksi
    ! ----------------
    do k=ndz_e,nfz_e
       do j=ny-3,ny
          do i=-3,8
             ! contravariant velocity
             vc=uu(i,j,k)*y_eta(i,j)-vv(i,j,k)*x_eta(i,j)
             rovc=rho_n(i,j,k)*vc

             Frho(i,j,k)  = rovc
             Frhou(i,j,k) = rovc*uu(i,j,k)+prs(i,j,k)*y_eta(i,j)
             Frhov(i,j,k) = rovc*vv(i,j,k)-prs(i,j,k)*x_eta(i,j)
             Frhow(i,j,k) = rovc*ww(i,j,k)
             Frhoe(i,j,k) = (rhoe_n(i,j,k)+prs(i,j,k))*vc
          enddo
       enddo
    enddo

    ! derivatives along ksi
    ! ---------------------
    i=1
    do j=ny-3,ny
       do k=ndz_e,nfz_e
          Krho(i,j,k)  = ( Frho(i+1,j,k)- Frho(i,j,k))
          Krhou(i,j,k) = (Frhou(i+1,j,k)-Frhou(i,j,k))
          Krhov(i,j,k) = (Frhov(i+1,j,k)-Frhov(i,j,k))
          Krhow(i,j,k) = (Frhow(i+1,j,k)-Frhow(i,j,k))
          Krhoe(i,j,k) = (Frhoe(i+1,j,k)-Frhoe(i,j,k))
       enddo
    enddo

    i=2
    do j=ny-3,ny
       do k=ndz_e,nfz_e
          Krho(i,j,k)  = 0.5_wp*( Frho(i+1,j,k)- Frho(i-1,j,k))
          Krhou(i,j,k) = 0.5_wp*(Frhou(i+1,j,k)-Frhou(i-1,j,k))
          Krhov(i,j,k) = 0.5_wp*(Frhov(i+1,j,k)-Frhov(i-1,j,k))
          Krhow(i,j,k) = 0.5_wp*(Frhow(i+1,j,k)-Frhow(i-1,j,k))
          Krhoe(i,j,k) = 0.5_wp*(Frhoe(i+1,j,k)-Frhoe(i-1,j,k))
       enddo
    enddo

    i=3
    do j=ny-3,ny
       do k=ndz_e,nfz_e
          Krho(i,j,k)  = a5(1)*( Frho(i+1,j,k)- Frho(i-1,j,k)) &
                       + a5(2)*( Frho(i+2,j,k)- Frho(i-2,j,k))
          Krhou(i,j,k) = a5(1)*(Frhou(i+1,j,k)-Frhou(i-1,j,k)) &
                       + a5(2)*(Frhou(i+2,j,k)-Frhou(i-2,j,k))
          Krhov(i,j,k) = a5(1)*(Frhov(i+1,j,k)-Frhov(i-1,j,k)) &
                       + a5(2)*(Frhov(i+2,j,k)-Frhov(i-2,j,k))
          Krhow(i,j,k) = a5(1)*(Frhow(i+1,j,k)-Frhow(i-1,j,k)) &
                       + a5(2)*(Frhow(i+2,j,k)-Frhow(i-2,j,k))
          Krhoe(i,j,k) = a5(1)*(Frhoe(i+1,j,k)-Frhoe(i-1,j,k)) &
                       + a5(2)*(Frhoe(i+2,j,k)-Frhoe(i-2,j,k))
       enddo
    enddo

    i=4
    do j=ny-3,ny
       do k=ndz_e,nfz_e
          Krho(i,j,k)  = a7(1) * ( Frho(i+1,j,k)-Frho(i-1,j,k) ) &
                       + a7(2) * ( Frho(i+2,j,k)-Frho(i-2,j,k) ) &
                       + a7(3) * ( Frho(i+3,j,k)-Frho(i-3,j,k) )
          Krhou(i,j,k) = a7(1) * ( Frhou(i+1,j,k)-Frhou(i-1,j,k) ) &
                       + a7(2) * ( Frhou(i+2,j,k)-Frhou(i-2,j,k) ) &
                       + a7(3) * ( Frhou(i+3,j,k)-Frhou(i-3,j,k) )
          Krhov(i,j,k) = a7(1) * ( Frhov(i+1,j,k)-Frhov(i-1,j,k) ) &
                       + a7(2) * ( Frhov(i+2,j,k)-Frhov(i-2,j,k) ) &
                       + a7(3) * ( Frhov(i+3,j,k)-Frhov(i-3,j,k) )
          Krhow(i,j,k) = a7(1) * ( Frhow(i+1,j,k)-Frhow(i-1,j,k) ) &
                       + a7(2) * ( Frhow(i+2,j,k)-Frhow(i-2,j,k) ) &
                       + a7(3) * ( Frhow(i+3,j,k)-Frhow(i-3,j,k) )
          Krhoe(i,j,k) = a7(1) * ( Frhoe(i+1,j,k)-Frhoe(i-1,j,k) ) &
                       + a7(2) * ( Frhoe(i+2,j,k)-Frhoe(i-2,j,k) ) &
                       + a7(3) * ( Frhoe(i+3,j,k)-Frhoe(i-3,j,k) )
       enddo
    enddo

    ! modified fluxes along eta
    ! -------------------------
    do k=ndz_e,nfz_e
       do i=1,4
          do j=ny-7,ny+4
             ! contravariant velocity
             vc=vv(i,j,k)*x_ksi(i,j)-uu(i,j,k)*y_ksi(i,j)
             rovc=rho_n(i,j,k)*vc

             Frho(i,j,k)  =  rovc
             Frhou(i,j,k) =  rovc*uu(i,j,k)-prs(i,j,k)*y_ksi(i,j)
             Frhov(i,j,k) =  rovc*vv(i,j,k)+prs(i,j,k)*x_ksi(i,j)
             Frhow(i,j,k) =  rovc*ww(i,j,k)
             Frhoe(i,j,k) = (rhoe_n(i,j,k)+prs(i,j,k))*vc
          enddo
       enddo
    enddo

    ! derivatives along eta
    ! ---------------------
    j=ny
    do i=1,4
       do k=ndz_e,nfz_e  
          Krho(i,j,k)  = (Frho(i,j,k) -Frho(i,j-1,k)) + Krho(i,j,k)
          Krhou(i,j,k) = (Frhou(i,j,k)-Frhou(i,j-1,k)) + Krhou(i,j,k)
          Krhov(i,j,k) = (Frhov(i,j,k)-Frhov(i,j-1,k)) + Krhov(i,j,k)
          Krhow(i,j,k) = (Frhow(i,j,k)-Frhow(i,j-1,k)) + Krhow(i,j,k)
          Krhoe(i,j,k) = (Frhoe(i,j,k)-Frhoe(i,j-1,k)) + Krhoe(i,j,k)
       enddo
    enddo

    j=ny-1
    do i=1,4
       do k=ndz_e,nfz_e  
          Krho(i,j,k)  = 0.5_wp*(Frho(i,j+1,k) -Frho(i,j-1,k)) + Krho(i,j,k)
          Krhou(i,j,k) = 0.5_wp*(Frhou(i,j+1,k)-Frhou(i,j-1,k)) + Krhou(i,j,k)
          Krhov(i,j,k) = 0.5_wp*(Frhov(i,j+1,k)-Frhov(i,j-1,k)) + Krhov(i,j,k)
          Krhow(i,j,k) = 0.5_wp*(Frhow(i,j+1,k)-Frhow(i,j-1,k)) + Krhow(i,j,k)
          Krhoe(i,j,k) = 0.5_wp*(Frhoe(i,j+1,k)-Frhoe(i,j-1,k)) + Krhoe(i,j,k)
       enddo
    enddo

    j=ny-2
    do i=1,4
       do k=ndz_e,nfz_e
          Krho(i,j,k)  = a5(1)*(Frho(i,j+1,k)-Frho(i,j-1,k)) &
                       + a5(2)*(Frho(i,j+2,k)-Frho(i,j-2,k)) + Krho(i,j,k)
          Krhou(i,j,k) = a5(1)*(Frhou(i,j+1,k)-Frhou(i,j-1,k)) &
                       + a5(2)*(Frhou(i,j+2,k)-Frhou(i,j-2,k)) + Krhou(i,j,k)
          Krhov(i,j,k) = a5(1)*(Frhov(i,j+1,k)-Frhov(i,j-1,k)) &
                       + a5(2)*(Frhov(i,j+2,k)-Frhov(i,j-2,k)) + Krhov(i,j,k)
          Krhow(i,j,k) = a5(1)*(Frhow(i,j+1,k)-Frhow(i,j-1,k)) &
                       + a5(2)*(Frhow(i,j+2,k)-Frhow(i,j-2,k)) + Krhow(i,j,k)
          Krhoe(i,j,k) = a5(1)*(Frhoe(i,j+1,k)-Frhoe(i,j-1,k)) &
                       + a5(2)*(Frhoe(i,j+2,k)-Frhoe(i,j-2,k)) + Krhoe(i,j,k)
       enddo
    enddo

    j=ny-3
    do i=1,4
       do k=ndz_e,nfz_e
          Krho(i,j,k)  = a7(1) * ( Frho(i,j+1,k)-Frho(i,j-1,k) ) &
                       + a7(2) * ( Frho(i,j+2,k)-Frho(i,j-2,k) ) &
                       + a7(3) * ( Frho(i,j+3,k)-Frho(i,j-3,k) ) + Krho(i,j,k)
          Krhou(i,j,k) = a7(1) * ( Frhou(i,j+1,k)-Frhou(i,j-1,k) ) &
                       + a7(2) * ( Frhou(i,j+2,k)-Frhou(i,j-2,k) ) &
                       + a7(3) * ( Frhou(i,j+3,k)-Frhou(i,j-3,k) ) + Krhou(i,j,k)
          Krhov(i,j,k) = a7(1) * ( Frhov(i,j+1,k)-Frhov(i,j-1,k) ) &
                       + a7(2) * ( Frhov(i,j+2,k)-Frhov(i,j-2,k) ) &
                       + a7(3) * ( Frhov(i,j+3,k)-Frhov(i,j-3,k) ) + Krhov(i,j,k)
          Krhow(i,j,k) = a7(1) * ( Frhow(i,j+1,k)-Frhow(i,j-1,k) ) &
                       + a7(2) * ( Frhow(i,j+2,k)-Frhow(i,j-2,k) ) &
                       + a7(3) * ( Frhow(i,j+3,k)-Frhow(i,j-3,k) ) + Krhow(i,j,k)
          Krhoe(i,j,k) = a7(1) * ( Frhoe(i,j+1,k)-Frhoe(i,j-1,k) ) &
                       + a7(2) * ( Frhoe(i,j+2,k)-Frhoe(i,j-2,k) ) &
                       + a7(3) * ( Frhoe(i,j+3,k)-Frhoe(i,j-3,k) ) + Krhoe(i,j,k)
       enddo
    enddo

    ! ~> Second pass: order reduction for metrics & 9-pt derivatives
    ! ==============================================================

    ! modified fluxes along ksi
    ! -------------------------
    do k=ndz_e,nfz_e
       do j=ny-3,ny
          do i=-3,8
             ! contravariant velocity
             vc=uu(i,j,k)*y_eta_imin_jmax(i,j)-vv(i,j,k)*x_eta_imin_jmax(i,j)
             rovc=rho_n(i,j,k)*vc

             Frho(i,j,k)  = rovc
             Frhou(i,j,k) = rovc*uu(i,j,k)+prs(i,j,k)*y_eta_imin_jmax(i,j)
             Frhov(i,j,k) = rovc*vv(i,j,k)-prs(i,j,k)*x_eta_imin_jmax(i,j)
             Frhow(i,j,k) = rovc*ww(i,j,k)
             Frhoe(i,j,k) = (rhoe_n(i,j,k)+prs(i,j,k))*vc
          enddo
       enddo
    enddo

    ! derivatives along ksi
    ! ---------------------
    do k=ndz_e,nfz_e
       do j=ny-3,ny
          do i=1,4
             Krho(i,j,k)  = a9(1) * ( Frho(i+1,j,k)-Frho(i-1,j,k) )   &
                          + a9(2) * ( Frho(i+2,j,k)-Frho(i-2,j,k) )   &
                          + a9(3) * ( Frho(i+3,j,k)-Frho(i-3,j,k) )   &
                          + a9(4) * ( Frho(i+4,j,k)-Frho(i-4,j,k) ) + Krho(i,j,k)

             Krhou(i,j,k) = a9(1) * ( Frhou(i+1,j,k)-Frhou(i-1,j,k) ) &
                          + a9(2) * ( Frhou(i+2,j,k)-Frhou(i-2,j,k) ) &
                          + a9(3) * ( Frhou(i+3,j,k)-Frhou(i-3,j,k) ) &
                          + a9(4) * ( Frhou(i+4,j,k)-Frhou(i-4,j,k) ) + Krhou(i,j,k)

             Krhov(i,j,k) = a9(1) * ( Frhov(i+1,j,k)-Frhov(i-1,j,k) ) &
                          + a9(2) * ( Frhov(i+2,j,k)-Frhov(i-2,j,k) ) &
                          + a9(3) * ( Frhov(i+3,j,k)-Frhov(i-3,j,k) ) &
                          + a9(4) * ( Frhov(i+4,j,k)-Frhov(i-4,j,k) ) + Krhov(i,j,k)

             Krhow(i,j,k) = a9(1) * ( Frhow(i+1,j,k)-Frhow(i-1,j,k) ) &
                          + a9(2) * ( Frhow(i+2,j,k)-Frhow(i-2,j,k) ) &
                          + a9(3) * ( Frhow(i+3,j,k)-Frhow(i-3,j,k) ) &
                          + a9(4) * ( Frhow(i+4,j,k)-Frhow(i-4,j,k) ) + Krhow(i,j,k)

             Krhoe(i,j,k) = a9(1) * ( Frhoe(i+1,j,k)-Frhoe(i-1,j,k) ) &
                          + a9(2) * ( Frhoe(i+2,j,k)-Frhoe(i-2,j,k) ) &
                          + a9(3) * ( Frhoe(i+3,j,k)-Frhoe(i-3,j,k) ) &
                          + a9(4) * ( Frhoe(i+4,j,k)-Frhoe(i-4,j,k) ) + Krhoe(i,j,k)
          enddo
       enddo
    enddo

    ! modified fluxes along eta
    ! -------------------------
    do k=ndz_e,nfz_e
       do i=1,4
          do j=ny-7,ny+4
             ! contravariant velocity
             vc=vv(i,j,k)*x_ksi_imin_jmax(i,j)-uu(i,j,k)*y_ksi_imin_jmax(i,j)
             rovc=rho_n(i,j,k)*vc

             Frho(i,j,k)  =  rovc
             Frhou(i,j,k) =  rovc*uu(i,j,k)-prs(i,j,k)*y_ksi_imin_jmax(i,j)
             Frhov(i,j,k) =  rovc*vv(i,j,k)+prs(i,j,k)*x_ksi_imin_jmax(i,j)
             Frhow(i,j,k) =  rovc*ww(i,j,k)
             Frhoe(i,j,k) = (rhoe_n(i,j,k)+prs(i,j,k))*vc
          enddo
       enddo
    enddo

    ! derivatives along eta
    ! ---------------------
    do k=ndz_e,nfz_e
       do j=ny-3,ny
          do i=1,4
             Krho(i,j,k)  = a9(1) * ( Frho(i,j+1,k)-Frho(i,j-1,k) ) &
                          + a9(2) * ( Frho(i,j+2,k)-Frho(i,j-2,k) ) &
                          + a9(3) * ( Frho(i,j+3,k)-Frho(i,j-3,k) ) &
                          + a9(4) * ( Frho(i,j+4,k)-Frho(i,j-4,k) ) &
                          + Krho(i,j,k)

             Krhou(i,j,k) = a9(1) * ( Frhou(i,j+1,k)-Frhou(i,j-1,k) ) &
                          + a9(2) * ( Frhou(i,j+2,k)-Frhou(i,j-2,k) ) &
                          + a9(3) * ( Frhou(i,j+3,k)-Frhou(i,j-3,k) ) &
                          + a9(4) * ( Frhou(i,j+4,k)-Frhou(i,j-4,k) ) &
                          + Krhou(i,j,k)

             Krhov(i,j,k) = a9(1) * ( Frhov(i,j+1,k)-Frhov(i,j-1,k) ) &
                          + a9(2) * ( Frhov(i,j+2,k)-Frhov(i,j-2,k) ) &
                          + a9(3) * ( Frhov(i,j+3,k)-Frhov(i,j-3,k) ) &
                          + a9(4) * ( Frhov(i,j+4,k)-Frhov(i,j-4,k) ) &
                          + Krhov(i,j,k)

             Krhow(i,j,k) = a9(1) * ( Frhow(i,j+1,k)-Frhow(i,j-1,k) ) &
                          + a9(2) * ( Frhow(i,j+2,k)-Frhow(i,j-2,k) ) &
                          + a9(3) * ( Frhow(i,j+3,k)-Frhow(i,j-3,k) ) &
                          + a9(4) * ( Frhow(i,j+4,k)-Frhow(i,j-4,k) ) &
                          + Krhow(i,j,k)

             Krhoe(i,j,k) = a9(1) * ( Frhoe(i,j+1,k)-Frhoe(i,j-1,k) ) &
                          + a9(2) * ( Frhoe(i,j+2,k)-Frhoe(i,j-2,k) ) &
                          + a9(3) * ( Frhoe(i,j+3,k)-Frhoe(i,j-3,k) ) &
                          + a9(4) * ( Frhoe(i,j+4,k)-Frhoe(i,j-4,k) ) &
                          + Krhoe(i,j,k)
          enddo
       enddo
    enddo

    ! Average of the two derivative evaluations
    ! =========================================
    do k=ndz_e,nfz_e
       do j=ny-3,ny
          do i=1,4
             Krho(i,j,k)  = Krho(i,j,k)*0.5
             Krhou(i,j,k) = Krhou(i,j,k)*0.5
             Krhov(i,j,k) = Krhov(i,j,k)*0.5
             Krhow(i,j,k) = Krhow(i,j,k)*0.5
             Krhoe(i,j,k) = Krhoe(i,j,k)*0.5
          enddo
       enddo
    enddo

  end subroutine flux_euler_w_imin_jmax_9pts_c

  !===============================================================================
  module subroutine flux_euler_w_imax_jmin_9pts_c
  !===============================================================================
    !> Derivatives of Euler fluxes for edge I_MAX/J_MIN close to two walls
    !> - curvilinear coordinate - 9-point stencil -
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: vc,rovc
    ! ---------------------------------------------------------------------------

    ! ~> First pass: 9-pt metrics & order reduction for derivatives
    ! =============================================================

    ! modified fluxes along ksi
    ! -------------------------
    do k=ndz_e,nfz_e
       do j=1,4
          do i=nx-7,nx+4
             ! contravariant velocity
             vc=uu(i,j,k)*y_eta(i,j)-vv(i,j,k)*x_eta(i,j)
             rovc=rho_n(i,j,k)*vc

             Frho(i,j,k)  = rovc
             Frhou(i,j,k) = rovc*uu(i,j,k)+prs(i,j,k)*y_eta(i,j)
             Frhov(i,j,k) = rovc*vv(i,j,k)-prs(i,j,k)*x_eta(i,j)
             Frhow(i,j,k) = rovc*ww(i,j,k)
             Frhoe(i,j,k) = (rhoe_n(i,j,k)+prs(i,j,k))*vc
          enddo
       enddo
    enddo

    ! derivatives along ksi
    ! ---------------------
    i=nx-3
    do j=1,4
       do k=ndz_e,nfz_e
          Krho(i,j,k)  = a7(1) * ( Frho(i+1,j,k)-Frho(i-1,j,k) ) &
                       + a7(2) * ( Frho(i+2,j,k)-Frho(i-2,j,k) ) &
                       + a7(3) * ( Frho(i+3,j,k)-Frho(i-3,j,k) )
          Krhou(i,j,k) = a7(1) * ( Frhou(i+1,j,k)-Frhou(i-1,j,k) ) &
                       + a7(2) * ( Frhou(i+2,j,k)-Frhou(i-2,j,k) ) &
                       + a7(3) * ( Frhou(i+3,j,k)-Frhou(i-3,j,k) )
          Krhov(i,j,k) = a7(1) * ( Frhov(i+1,j,k)-Frhov(i-1,j,k) ) &
                       + a7(2) * ( Frhov(i+2,j,k)-Frhov(i-2,j,k) ) &
                       + a7(3) * ( Frhov(i+3,j,k)-Frhov(i-3,j,k) )
          Krhow(i,j,k) = a7(1) * ( Frhow(i+1,j,k)-Frhow(i-1,j,k) ) &
                       + a7(2) * ( Frhow(i+2,j,k)-Frhow(i-2,j,k) ) &
                       + a7(3) * ( Frhow(i+3,j,k)-Frhow(i-3,j,k) )
          Krhoe(i,j,k) = a7(1) * ( Frhoe(i+1,j,k)-Frhoe(i-1,j,k) ) &
                       + a7(2) * ( Frhoe(i+2,j,k)-Frhoe(i-2,j,k) ) &
                       + a7(3) * ( Frhoe(i+3,j,k)-Frhoe(i-3,j,k) )
       enddo
    enddo

    i=nx-2
    do j=1,4
       do k=ndz_e,nfz_e
          Krho(i,j,k)  = a5(1)*( Frho(i+1,j,k)- Frho(i-1,j,k)) &
                       + a5(2)*( Frho(i+2,j,k)- Frho(i-2,j,k))
          Krhou(i,j,k) = a5(1)*(Frhou(i+1,j,k)-Frhou(i-1,j,k)) &
                       + a5(2)*(Frhou(i+2,j,k)-Frhou(i-2,j,k))
          Krhov(i,j,k) = a5(1)*(Frhov(i+1,j,k)-Frhov(i-1,j,k)) &
                       + a5(2)*(Frhov(i+2,j,k)-Frhov(i-2,j,k))
          Krhow(i,j,k) = a5(1)*(Frhow(i+1,j,k)-Frhow(i-1,j,k)) &
                       + a5(2)*(Frhow(i+2,j,k)-Frhow(i-2,j,k))
          Krhoe(i,j,k) = a5(1)*(Frhoe(i+1,j,k)-Frhoe(i-1,j,k)) &
                       + a5(2)*(Frhoe(i+2,j,k)-Frhoe(i-2,j,k))
       enddo
    enddo

    i=nx-1
    do j=1,4
       do k=ndz_e,nfz_e
          Krho(i,j,k)  = 0.5_wp*( Frho(i+1,j,k)- Frho(i-1,j,k))
          Krhou(i,j,k) = 0.5_wp*(Frhou(i+1,j,k)-Frhou(i-1,j,k))
          Krhov(i,j,k) = 0.5_wp*(Frhov(i+1,j,k)-Frhov(i-1,j,k))
          Krhow(i,j,k) = 0.5_wp*(Frhow(i+1,j,k)-Frhow(i-1,j,k))
          Krhoe(i,j,k) = 0.5_wp*(Frhoe(i+1,j,k)-Frhoe(i-1,j,k))
       enddo
    enddo

    i=nx
    do j=1,4
       do k=ndz_e,nfz_e
          Krho(i,j,k)  = ( Frho(i,j,k)- Frho(i-1,j,k))
          Krhou(i,j,k) = (Frhou(i,j,k)-Frhou(i-1,j,k))
          Krhov(i,j,k) = (Frhov(i,j,k)-Frhov(i-1,j,k))
          Krhow(i,j,k) = (Frhow(i,j,k)-Frhow(i-1,j,k))
          Krhoe(i,j,k) = (Frhoe(i,j,k)-Frhoe(i-1,j,k))
       enddo
    enddo

    ! modified fluxes along eta
    ! -------------------------
    do k=ndz_e,nfz_e
       do i=nx-3,nx
          do j=-3,8
             ! contravariant velocity
             vc=vv(i,j,k)*x_ksi(i,j)-uu(i,j,k)*y_ksi(i,j)
             rovc=rho_n(i,j,k)*vc

             Frho(i,j,k)  =  rovc
             Frhou(i,j,k) =  rovc*uu(i,j,k)-prs(i,j,k)*y_ksi(i,j)
             Frhov(i,j,k) =  rovc*vv(i,j,k)+prs(i,j,k)*x_ksi(i,j)
             Frhow(i,j,k) =  rovc*ww(i,j,k)
             Frhoe(i,j,k) = (rhoe_n(i,j,k)+prs(i,j,k))*vc
          enddo
       enddo
    enddo

    ! derivatives along eta
    ! ---------------------
    j=1
    do i=nx-3,nx
       do k=ndz_e,nfz_e  
          Krho(i,j,k)  = (Frho(i,j+1,k) -Frho(i,j,k)) + Krho(i,j,k)
          Krhou(i,j,k) = (Frhou(i,j+1,k)-Frhou(i,j,k)) + Krhou(i,j,k)
          Krhov(i,j,k) = (Frhov(i,j+1,k)-Frhov(i,j,k)) + Krhov(i,j,k)
          Krhow(i,j,k) = (Frhow(i,j+1,k)-Frhow(i,j,k)) + Krhow(i,j,k)
          Krhoe(i,j,k) = (Frhoe(i,j+1,k)-Frhoe(i,j,k)) + Krhoe(i,j,k)
       enddo
    enddo

    j=2
    do i=nx-3,nx
       do k=ndz_e,nfz_e  
          Krho(i,j,k)  = 0.5_wp*(Frho(i,j+1,k) -Frho(i,j-1,k)) + Krho(i,j,k)
          Krhou(i,j,k) = 0.5_wp*(Frhou(i,j+1,k)-Frhou(i,j-1,k)) + Krhou(i,j,k)
          Krhov(i,j,k) = 0.5_wp*(Frhov(i,j+1,k)-Frhov(i,j-1,k)) + Krhov(i,j,k)
          Krhow(i,j,k) = 0.5_wp*(Frhow(i,j+1,k)-Frhow(i,j-1,k)) + Krhow(i,j,k)
          Krhoe(i,j,k) = 0.5_wp*(Frhoe(i,j+1,k)-Frhoe(i,j-1,k)) + Krhoe(i,j,k)
       enddo
    enddo

    j=3
    do i=nx-3,nx
       do k=ndz_e,nfz_e
          Krho(i,j,k)  = a5(1)*(Frho(i,j+1,k)-Frho(i,j-1,k)) &
                       + a5(2)*(Frho(i,j+2,k)-Frho(i,j-2,k)) + Krho(i,j,k)
          Krhou(i,j,k) = a5(1)*(Frhou(i,j+1,k)-Frhou(i,j-1,k)) &
                       + a5(2)*(Frhou(i,j+2,k)-Frhou(i,j-2,k)) + Krhou(i,j,k)
          Krhov(i,j,k) = a5(1)*(Frhov(i,j+1,k)-Frhov(i,j-1,k)) &
                       + a5(2)*(Frhov(i,j+2,k)-Frhov(i,j-2,k)) + Krhov(i,j,k)
          Krhow(i,j,k) = a5(1)*(Frhow(i,j+1,k)-Frhow(i,j-1,k)) &
                       + a5(2)*(Frhow(i,j+2,k)-Frhow(i,j-2,k)) + Krhow(i,j,k)
          Krhoe(i,j,k) = a5(1)*(Frhoe(i,j+1,k)-Frhoe(i,j-1,k)) &
                       + a5(2)*(Frhoe(i,j+2,k)-Frhoe(i,j-2,k)) + Krhoe(i,j,k)
       enddo
    enddo

    j=4
    do i=nx-3,nx
       do k=ndz_e,nfz_e
          Krho(i,j,k)  = a7(1) * ( Frho(i,j+1,k)-Frho(i,j-1,k) ) &
                       + a7(2) * ( Frho(i,j+2,k)-Frho(i,j-2,k) ) &
                       + a7(3) * ( Frho(i,j+3,k)-Frho(i,j-3,k) ) + Krho(i,j,k)
          Krhou(i,j,k) = a7(1) * ( Frhou(i,j+1,k)-Frhou(i,j-1,k) ) &
                       + a7(2) * ( Frhou(i,j+2,k)-Frhou(i,j-2,k) ) &
                       + a7(3) * ( Frhou(i,j+3,k)-Frhou(i,j-3,k) ) + Krhou(i,j,k)
          Krhov(i,j,k) = a7(1) * ( Frhov(i,j+1,k)-Frhov(i,j-1,k) ) &
                       + a7(2) * ( Frhov(i,j+2,k)-Frhov(i,j-2,k) ) &
                       + a7(3) * ( Frhov(i,j+3,k)-Frhov(i,j-3,k) ) + Krhov(i,j,k)
          Krhow(i,j,k) = a7(1) * ( Frhow(i,j+1,k)-Frhow(i,j-1,k) ) &
                       + a7(2) * ( Frhow(i,j+2,k)-Frhow(i,j-2,k) ) &
                       + a7(3) * ( Frhow(i,j+3,k)-Frhow(i,j-3,k) ) + Krhow(i,j,k)
          Krhoe(i,j,k) = a7(1) * ( Frhoe(i,j+1,k)-Frhoe(i,j-1,k) ) &
                       + a7(2) * ( Frhoe(i,j+2,k)-Frhoe(i,j-2,k) ) &
                       + a7(3) * ( Frhoe(i,j+3,k)-Frhoe(i,j-3,k) ) + Krhoe(i,j,k)
       enddo
    enddo

    ! ~> Second pass: order reduction for metrics & 9-pt derivatives
    ! ==============================================================

    ! modified fluxes along ksi
    ! -------------------------
    do k=ndz_e,nfz_e
       do j=1,4
          do i=nx-7,nx+4
             ! contravariant velocity
             vc=uu(i,j,k)*y_eta_imax_jmin(i,j)-vv(i,j,k)*x_eta_imax_jmin(i,j)
             rovc=rho_n(i,j,k)*vc

             Frho(i,j,k)  = rovc
             Frhou(i,j,k) = rovc*uu(i,j,k)+prs(i,j,k)*y_eta_imax_jmin(i,j)
             Frhov(i,j,k) = rovc*vv(i,j,k)-prs(i,j,k)*x_eta_imax_jmin(i,j)
             Frhow(i,j,k) = rovc*ww(i,j,k)
             Frhoe(i,j,k) = (rhoe_n(i,j,k)+prs(i,j,k))*vc
          enddo
       enddo
    enddo

    ! derivatives along ksi
    ! ---------------------
    do k=ndz_e,nfz_e
       do j=1,4
          do i=nx-3,nx
             Krho(i,j,k)  = a9(1) * ( Frho(i+1,j,k)-Frho(i-1,j,k) )   &
                          + a9(2) * ( Frho(i+2,j,k)-Frho(i-2,j,k) )   &
                          + a9(3) * ( Frho(i+3,j,k)-Frho(i-3,j,k) )   &
                          + a9(4) * ( Frho(i+4,j,k)-Frho(i-4,j,k) ) + Krho(i,j,k)

             Krhou(i,j,k) = a9(1) * ( Frhou(i+1,j,k)-Frhou(i-1,j,k) ) &
                          + a9(2) * ( Frhou(i+2,j,k)-Frhou(i-2,j,k) ) &
                          + a9(3) * ( Frhou(i+3,j,k)-Frhou(i-3,j,k) ) &
                          + a9(4) * ( Frhou(i+4,j,k)-Frhou(i-4,j,k) ) + Krhou(i,j,k)

             Krhov(i,j,k) = a9(1) * ( Frhov(i+1,j,k)-Frhov(i-1,j,k) ) &
                          + a9(2) * ( Frhov(i+2,j,k)-Frhov(i-2,j,k) ) &
                          + a9(3) * ( Frhov(i+3,j,k)-Frhov(i-3,j,k) ) &
                          + a9(4) * ( Frhov(i+4,j,k)-Frhov(i-4,j,k) ) + Krhov(i,j,k)

             Krhow(i,j,k) = a9(1) * ( Frhow(i+1,j,k)-Frhow(i-1,j,k) ) &
                          + a9(2) * ( Frhow(i+2,j,k)-Frhow(i-2,j,k) ) &
                          + a9(3) * ( Frhow(i+3,j,k)-Frhow(i-3,j,k) ) &
                          + a9(4) * ( Frhow(i+4,j,k)-Frhow(i-4,j,k) ) + Krhow(i,j,k)

             Krhoe(i,j,k) = a9(1) * ( Frhoe(i+1,j,k)-Frhoe(i-1,j,k) ) &
                          + a9(2) * ( Frhoe(i+2,j,k)-Frhoe(i-2,j,k) ) &
                          + a9(3) * ( Frhoe(i+3,j,k)-Frhoe(i-3,j,k) ) &
                          + a9(4) * ( Frhoe(i+4,j,k)-Frhoe(i-4,j,k) ) + Krhoe(i,j,k)
          enddo
       enddo
    enddo

    ! modified fluxes along eta
    ! -------------------------
    do k=ndz_e,nfz_e
       do i=nx-3,nx
          do j=-3,8
             ! contravariant velocity
             vc=vv(i,j,k)*x_ksi_imax_jmin(i,j)-uu(i,j,k)*y_ksi_imax_jmin(i,j)
             rovc=rho_n(i,j,k)*vc

             Frho(i,j,k)  =  rovc
             Frhou(i,j,k) =  rovc*uu(i,j,k)-prs(i,j,k)*y_ksi_imax_jmin(i,j)
             Frhov(i,j,k) =  rovc*vv(i,j,k)+prs(i,j,k)*x_ksi_imax_jmin(i,j)
             Frhow(i,j,k) =  rovc*ww(i,j,k)
             Frhoe(i,j,k) = (rhoe_n(i,j,k)+prs(i,j,k))*vc
          enddo
       enddo
    enddo

    ! derivatives along eta
    ! ---------------------
    do k=ndz_e,nfz_e
       do j=1,4
          do i=nx-3,nx
             Krho(i,j,k)  = a9(1) * ( Frho(i,j+1,k)-Frho(i,j-1,k) ) &
                          + a9(2) * ( Frho(i,j+2,k)-Frho(i,j-2,k) ) &
                          + a9(3) * ( Frho(i,j+3,k)-Frho(i,j-3,k) ) &
                          + a9(4) * ( Frho(i,j+4,k)-Frho(i,j-4,k) ) &
                          + Krho(i,j,k)

             Krhou(i,j,k) = a9(1) * ( Frhou(i,j+1,k)-Frhou(i,j-1,k) ) &
                          + a9(2) * ( Frhou(i,j+2,k)-Frhou(i,j-2,k) ) &
                          + a9(3) * ( Frhou(i,j+3,k)-Frhou(i,j-3,k) ) &
                          + a9(4) * ( Frhou(i,j+4,k)-Frhou(i,j-4,k) ) &
                          + Krhou(i,j,k)

             Krhov(i,j,k) = a9(1) * ( Frhov(i,j+1,k)-Frhov(i,j-1,k) ) &
                          + a9(2) * ( Frhov(i,j+2,k)-Frhov(i,j-2,k) ) &
                          + a9(3) * ( Frhov(i,j+3,k)-Frhov(i,j-3,k) ) &
                          + a9(4) * ( Frhov(i,j+4,k)-Frhov(i,j-4,k) ) &
                          + Krhov(i,j,k)

             Krhow(i,j,k) = a9(1) * ( Frhow(i,j+1,k)-Frhow(i,j-1,k) ) &
                          + a9(2) * ( Frhow(i,j+2,k)-Frhow(i,j-2,k) ) &
                          + a9(3) * ( Frhow(i,j+3,k)-Frhow(i,j-3,k) ) &
                          + a9(4) * ( Frhow(i,j+4,k)-Frhow(i,j-4,k) ) &
                          + Krhow(i,j,k)

             Krhoe(i,j,k) = a9(1) * ( Frhoe(i,j+1,k)-Frhoe(i,j-1,k) ) &
                          + a9(2) * ( Frhoe(i,j+2,k)-Frhoe(i,j-2,k) ) &
                          + a9(3) * ( Frhoe(i,j+3,k)-Frhoe(i,j-3,k) ) &
                          + a9(4) * ( Frhoe(i,j+4,k)-Frhoe(i,j-4,k) ) &
                          + Krhoe(i,j,k)
          enddo
       enddo
    enddo

    ! Average of the two derivative evaluations
    ! =========================================
    do k=ndz_e,nfz_e
       do j=1,4
          do i=nx-3,nx
             Krho(i,j,k)  = Krho(i,j,k)*0.5
             Krhou(i,j,k) = Krhou(i,j,k)*0.5
             Krhov(i,j,k) = Krhov(i,j,k)*0.5
             Krhow(i,j,k) = Krhow(i,j,k)*0.5
             Krhoe(i,j,k) = Krhoe(i,j,k)*0.5
          enddo
       enddo
    enddo

  end subroutine flux_euler_w_imax_jmin_9pts_c

  !===============================================================================
  module subroutine flux_euler_w_imax_jmax_9pts_c
  !===============================================================================
    !> Derivatives of Euler fluxes for edge I_MAX/J_MAX close to two walls
    !> - curvilinear coordinate - 9-point stencil -
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: vc,rovc
    ! ---------------------------------------------------------------------------

    ! ~> First pass: 9-pt metrics & order reduction for derivatives
    ! =============================================================

    ! modified fluxes along ksi
    ! -------------------------
    do k=ndz_e,nfz_e
       do j=ny-3,ny
          do i=nx-7,nx+4
             ! contravariant velocity
             vc=uu(i,j,k)*y_eta(i,j)-vv(i,j,k)*x_eta(i,j)
             rovc=rho_n(i,j,k)*vc

             Frho(i,j,k)  = rovc
             Frhou(i,j,k) = rovc*uu(i,j,k)+prs(i,j,k)*y_eta(i,j)
             Frhov(i,j,k) = rovc*vv(i,j,k)-prs(i,j,k)*x_eta(i,j)
             Frhow(i,j,k) = rovc*ww(i,j,k)
             Frhoe(i,j,k) = (rhoe_n(i,j,k)+prs(i,j,k))*vc
          enddo
       enddo
    enddo

    ! derivatives along ksi
    ! ---------------------
    i=nx-3
    do j=ny-3,ny
       do k=ndz_e,nfz_e
          Krho(i,j,k)  = a7(1) * ( Frho(i+1,j,k)-Frho(i-1,j,k) ) &
                       + a7(2) * ( Frho(i+2,j,k)-Frho(i-2,j,k) ) &
                       + a7(3) * ( Frho(i+3,j,k)-Frho(i-3,j,k) )
          Krhou(i,j,k) = a7(1) * ( Frhou(i+1,j,k)-Frhou(i-1,j,k) ) &
                       + a7(2) * ( Frhou(i+2,j,k)-Frhou(i-2,j,k) ) &
                       + a7(3) * ( Frhou(i+3,j,k)-Frhou(i-3,j,k) )
          Krhov(i,j,k) = a7(1) * ( Frhov(i+1,j,k)-Frhov(i-1,j,k) ) &
                       + a7(2) * ( Frhov(i+2,j,k)-Frhov(i-2,j,k) ) &
                       + a7(3) * ( Frhov(i+3,j,k)-Frhov(i-3,j,k) )
          Krhow(i,j,k) = a7(1) * ( Frhow(i+1,j,k)-Frhow(i-1,j,k) ) &
                       + a7(2) * ( Frhow(i+2,j,k)-Frhow(i-2,j,k) ) &
                       + a7(3) * ( Frhow(i+3,j,k)-Frhow(i-3,j,k) )
          Krhoe(i,j,k) = a7(1) * ( Frhoe(i+1,j,k)-Frhoe(i-1,j,k) ) &
                       + a7(2) * ( Frhoe(i+2,j,k)-Frhoe(i-2,j,k) ) &
                       + a7(3) * ( Frhoe(i+3,j,k)-Frhoe(i-3,j,k) )
       enddo
    enddo

    i=nx-2
    do j=ny-3,ny
       do k=ndz_e,nfz_e
          Krho(i,j,k)  = a5(1)*( Frho(i+1,j,k)- Frho(i-1,j,k)) &
                       + a5(2)*( Frho(i+2,j,k)- Frho(i-2,j,k))
          Krhou(i,j,k) = a5(1)*(Frhou(i+1,j,k)-Frhou(i-1,j,k)) &
                       + a5(2)*(Frhou(i+2,j,k)-Frhou(i-2,j,k))
          Krhov(i,j,k) = a5(1)*(Frhov(i+1,j,k)-Frhov(i-1,j,k)) &
                       + a5(2)*(Frhov(i+2,j,k)-Frhov(i-2,j,k))
          Krhow(i,j,k) = a5(1)*(Frhow(i+1,j,k)-Frhow(i-1,j,k)) &
                       + a5(2)*(Frhow(i+2,j,k)-Frhow(i-2,j,k))
          Krhoe(i,j,k) = a5(1)*(Frhoe(i+1,j,k)-Frhoe(i-1,j,k)) &
                       + a5(2)*(Frhoe(i+2,j,k)-Frhoe(i-2,j,k))
       enddo
    enddo

    i=nx-1
    do j=ny-3,ny
       do k=ndz_e,nfz_e
          Krho(i,j,k)  = 0.5_wp*( Frho(i+1,j,k)- Frho(i-1,j,k))
          Krhou(i,j,k) = 0.5_wp*(Frhou(i+1,j,k)-Frhou(i-1,j,k))
          Krhov(i,j,k) = 0.5_wp*(Frhov(i+1,j,k)-Frhov(i-1,j,k))
          Krhow(i,j,k) = 0.5_wp*(Frhow(i+1,j,k)-Frhow(i-1,j,k))
          Krhoe(i,j,k) = 0.5_wp*(Frhoe(i+1,j,k)-Frhoe(i-1,j,k))
       enddo
    enddo

    i=nx
    do j=ny-3,ny
       do k=ndz_e,nfz_e
          Krho(i,j,k)  = ( Frho(i,j,k)- Frho(i-1,j,k))
          Krhou(i,j,k) = (Frhou(i,j,k)-Frhou(i-1,j,k))
          Krhov(i,j,k) = (Frhov(i,j,k)-Frhov(i-1,j,k))
          Krhow(i,j,k) = (Frhow(i,j,k)-Frhow(i-1,j,k))
          Krhoe(i,j,k) = (Frhoe(i,j,k)-Frhoe(i-1,j,k))
       enddo
    enddo

    ! modified fluxes along eta
    ! -------------------------
    do k=ndz_e,nfz_e
       do i=nx-3,nx
          do j=ny-7,ny+4
             ! contravariant velocity
             vc=vv(i,j,k)*x_ksi(i,j)-uu(i,j,k)*y_ksi(i,j)
             rovc=rho_n(i,j,k)*vc

             Frho(i,j,k)  =  rovc
             Frhou(i,j,k) =  rovc*uu(i,j,k)-prs(i,j,k)*y_ksi(i,j)
             Frhov(i,j,k) =  rovc*vv(i,j,k)+prs(i,j,k)*x_ksi(i,j)
             Frhow(i,j,k) =  rovc*ww(i,j,k)
             Frhoe(i,j,k) = (rhoe_n(i,j,k)+prs(i,j,k))*vc
          enddo
       enddo
    enddo

    ! derivatives along eta
    ! ---------------------
    j=ny
    do i=nx-3,nx
       do k=ndz_e,nfz_e  
          Krho(i,j,k)  = (Frho(i,j,k) -Frho(i,j-1,k)) + Krho(i,j,k)
          Krhou(i,j,k) = (Frhou(i,j,k)-Frhou(i,j-1,k)) + Krhou(i,j,k)
          Krhov(i,j,k) = (Frhov(i,j,k)-Frhov(i,j-1,k)) + Krhov(i,j,k)
          Krhow(i,j,k) = (Frhow(i,j,k)-Frhow(i,j-1,k)) + Krhow(i,j,k)
          Krhoe(i,j,k) = (Frhoe(i,j,k)-Frhoe(i,j-1,k)) + Krhoe(i,j,k)
       enddo
    enddo

    j=ny-1
    do i=nx-3,nx
       do k=ndz_e,nfz_e  
          Krho(i,j,k)  = 0.5_wp*(Frho(i,j+1,k) -Frho(i,j-1,k)) + Krho(i,j,k)
          Krhou(i,j,k) = 0.5_wp*(Frhou(i,j+1,k)-Frhou(i,j-1,k)) + Krhou(i,j,k)
          Krhov(i,j,k) = 0.5_wp*(Frhov(i,j+1,k)-Frhov(i,j-1,k)) + Krhov(i,j,k)
          Krhow(i,j,k) = 0.5_wp*(Frhow(i,j+1,k)-Frhow(i,j-1,k)) + Krhow(i,j,k)
          Krhoe(i,j,k) = 0.5_wp*(Frhoe(i,j+1,k)-Frhoe(i,j-1,k)) + Krhoe(i,j,k)
       enddo
    enddo

    j=ny-2
    do i=nx-3,nx
       do k=ndz_e,nfz_e
          Krho(i,j,k)  = a5(1)*(Frho(i,j+1,k)-Frho(i,j-1,k)) &
                       + a5(2)*(Frho(i,j+2,k)-Frho(i,j-2,k)) + Krho(i,j,k)
          Krhou(i,j,k) = a5(1)*(Frhou(i,j+1,k)-Frhou(i,j-1,k)) &
                       + a5(2)*(Frhou(i,j+2,k)-Frhou(i,j-2,k)) + Krhou(i,j,k)
          Krhov(i,j,k) = a5(1)*(Frhov(i,j+1,k)-Frhov(i,j-1,k)) &
                       + a5(2)*(Frhov(i,j+2,k)-Frhov(i,j-2,k)) + Krhov(i,j,k)
          Krhow(i,j,k) = a5(1)*(Frhow(i,j+1,k)-Frhow(i,j-1,k)) &
                       + a5(2)*(Frhow(i,j+2,k)-Frhow(i,j-2,k)) + Krhow(i,j,k)
          Krhoe(i,j,k) = a5(1)*(Frhoe(i,j+1,k)-Frhoe(i,j-1,k)) &
                       + a5(2)*(Frhoe(i,j+2,k)-Frhoe(i,j-2,k)) + Krhoe(i,j,k)
       enddo
    enddo

    j=ny-3
    do i=nx-3,nx
       do k=ndz_e,nfz_e
          Krho(i,j,k)  = a7(1) * ( Frho(i,j+1,k)-Frho(i,j-1,k) ) &
                       + a7(2) * ( Frho(i,j+2,k)-Frho(i,j-2,k) ) &
                       + a7(3) * ( Frho(i,j+3,k)-Frho(i,j-3,k) ) + Krho(i,j,k)
          Krhou(i,j,k) = a7(1) * ( Frhou(i,j+1,k)-Frhou(i,j-1,k) ) &
                       + a7(2) * ( Frhou(i,j+2,k)-Frhou(i,j-2,k) ) &
                       + a7(3) * ( Frhou(i,j+3,k)-Frhou(i,j-3,k) ) + Krhou(i,j,k)
          Krhov(i,j,k) = a7(1) * ( Frhov(i,j+1,k)-Frhov(i,j-1,k) ) &
                       + a7(2) * ( Frhov(i,j+2,k)-Frhov(i,j-2,k) ) &
                       + a7(3) * ( Frhov(i,j+3,k)-Frhov(i,j-3,k) ) + Krhov(i,j,k)
          Krhow(i,j,k) = a7(1) * ( Frhow(i,j+1,k)-Frhow(i,j-1,k) ) &
                       + a7(2) * ( Frhow(i,j+2,k)-Frhow(i,j-2,k) ) &
                       + a7(3) * ( Frhow(i,j+3,k)-Frhow(i,j-3,k) ) + Krhow(i,j,k)
          Krhoe(i,j,k) = a7(1) * ( Frhoe(i,j+1,k)-Frhoe(i,j-1,k) ) &
                       + a7(2) * ( Frhoe(i,j+2,k)-Frhoe(i,j-2,k) ) &
                       + a7(3) * ( Frhoe(i,j+3,k)-Frhoe(i,j-3,k) ) + Krhoe(i,j,k)
       enddo
    enddo

    ! ~> Second pass: order reduction for metrics & 9-pt derivatives
    ! ==============================================================

    ! modified fluxes along ksi
    ! -------------------------
    do k=ndz_e,nfz_e
       do j=ny-3,ny
          do i=nx-7,nx+4
             ! contravariant velocity
             vc=uu(i,j,k)*y_eta_imax_jmax(i,j)-vv(i,j,k)*x_eta_imax_jmax(i,j)
             rovc=rho_n(i,j,k)*vc

             Frho(i,j,k)  = rovc
             Frhou(i,j,k) = rovc*uu(i,j,k)+prs(i,j,k)*y_eta_imax_jmax(i,j)
             Frhov(i,j,k) = rovc*vv(i,j,k)-prs(i,j,k)*x_eta_imax_jmax(i,j)
             Frhow(i,j,k) = rovc*ww(i,j,k)
             Frhoe(i,j,k) = (rhoe_n(i,j,k)+prs(i,j,k))*vc
          enddo
       enddo
    enddo

    ! derivatives along ksi
    ! ---------------------
    do k=ndz_e,nfz_e
       do j=ny-3,ny
          do i=nx-3,nx
             Krho(i,j,k)  = a9(1) * ( Frho(i+1,j,k)-Frho(i-1,j,k) )   &
                          + a9(2) * ( Frho(i+2,j,k)-Frho(i-2,j,k) )   &
                          + a9(3) * ( Frho(i+3,j,k)-Frho(i-3,j,k) )   &
                          + a9(4) * ( Frho(i+4,j,k)-Frho(i-4,j,k) ) + Krho(i,j,k)

             Krhou(i,j,k) = a9(1) * ( Frhou(i+1,j,k)-Frhou(i-1,j,k) ) &
                          + a9(2) * ( Frhou(i+2,j,k)-Frhou(i-2,j,k) ) &
                          + a9(3) * ( Frhou(i+3,j,k)-Frhou(i-3,j,k) ) &
                          + a9(4) * ( Frhou(i+4,j,k)-Frhou(i-4,j,k) ) + Krhou(i,j,k)

             Krhov(i,j,k) = a9(1) * ( Frhov(i+1,j,k)-Frhov(i-1,j,k) ) &
                          + a9(2) * ( Frhov(i+2,j,k)-Frhov(i-2,j,k) ) &
                          + a9(3) * ( Frhov(i+3,j,k)-Frhov(i-3,j,k) ) &
                          + a9(4) * ( Frhov(i+4,j,k)-Frhov(i-4,j,k) ) + Krhov(i,j,k)

             Krhow(i,j,k) = a9(1) * ( Frhow(i+1,j,k)-Frhow(i-1,j,k) ) &
                          + a9(2) * ( Frhow(i+2,j,k)-Frhow(i-2,j,k) ) &
                          + a9(3) * ( Frhow(i+3,j,k)-Frhow(i-3,j,k) ) &
                          + a9(4) * ( Frhow(i+4,j,k)-Frhow(i-4,j,k) ) + Krhow(i,j,k)

             Krhoe(i,j,k) = a9(1) * ( Frhoe(i+1,j,k)-Frhoe(i-1,j,k) ) &
                          + a9(2) * ( Frhoe(i+2,j,k)-Frhoe(i-2,j,k) ) &
                          + a9(3) * ( Frhoe(i+3,j,k)-Frhoe(i-3,j,k) ) &
                          + a9(4) * ( Frhoe(i+4,j,k)-Frhoe(i-4,j,k) ) + Krhoe(i,j,k)
          enddo
       enddo
    enddo

    ! modified fluxes along eta
    ! -------------------------
    do k=ndz_e,nfz_e
       do i=nx-3,nx
          do j=ny-7,ny+4
             ! contravariant velocity
             vc=vv(i,j,k)*x_ksi_imax_jmax(i,j)-uu(i,j,k)*y_ksi_imax_jmax(i,j)
             rovc=rho_n(i,j,k)*vc

             Frho(i,j,k)  =  rovc
             Frhou(i,j,k) =  rovc*uu(i,j,k)-prs(i,j,k)*y_ksi_imax_jmax(i,j)
             Frhov(i,j,k) =  rovc*vv(i,j,k)+prs(i,j,k)*x_ksi_imax_jmax(i,j)
             Frhow(i,j,k) =  rovc*ww(i,j,k)
             Frhoe(i,j,k) = (rhoe_n(i,j,k)+prs(i,j,k))*vc
          enddo
       enddo
    enddo

    ! derivatives along eta
    ! ---------------------
    do k=ndz_e,nfz_e
       do j=ny-3,ny
          do i=nx-3,nx
             Krho(i,j,k)  = a9(1) * ( Frho(i,j+1,k)-Frho(i,j-1,k) ) &
                          + a9(2) * ( Frho(i,j+2,k)-Frho(i,j-2,k) ) &
                          + a9(3) * ( Frho(i,j+3,k)-Frho(i,j-3,k) ) &
                          + a9(4) * ( Frho(i,j+4,k)-Frho(i,j-4,k) ) &
                          + Krho(i,j,k)

             Krhou(i,j,k) = a9(1) * ( Frhou(i,j+1,k)-Frhou(i,j-1,k) ) &
                          + a9(2) * ( Frhou(i,j+2,k)-Frhou(i,j-2,k) ) &
                          + a9(3) * ( Frhou(i,j+3,k)-Frhou(i,j-3,k) ) &
                          + a9(4) * ( Frhou(i,j+4,k)-Frhou(i,j-4,k) ) &
                          + Krhou(i,j,k)

             Krhov(i,j,k) = a9(1) * ( Frhov(i,j+1,k)-Frhov(i,j-1,k) ) &
                          + a9(2) * ( Frhov(i,j+2,k)-Frhov(i,j-2,k) ) &
                          + a9(3) * ( Frhov(i,j+3,k)-Frhov(i,j-3,k) ) &
                          + a9(4) * ( Frhov(i,j+4,k)-Frhov(i,j-4,k) ) &
                          + Krhov(i,j,k)

             Krhow(i,j,k) = a9(1) * ( Frhow(i,j+1,k)-Frhow(i,j-1,k) ) &
                          + a9(2) * ( Frhow(i,j+2,k)-Frhow(i,j-2,k) ) &
                          + a9(3) * ( Frhow(i,j+3,k)-Frhow(i,j-3,k) ) &
                          + a9(4) * ( Frhow(i,j+4,k)-Frhow(i,j-4,k) ) &
                          + Krhow(i,j,k)

             Krhoe(i,j,k) = a9(1) * ( Frhoe(i,j+1,k)-Frhoe(i,j-1,k) ) &
                          + a9(2) * ( Frhoe(i,j+2,k)-Frhoe(i,j-2,k) ) &
                          + a9(3) * ( Frhoe(i,j+3,k)-Frhoe(i,j-3,k) ) &
                          + a9(4) * ( Frhoe(i,j+4,k)-Frhoe(i,j-4,k) ) &
                          + Krhoe(i,j,k)
          enddo
       enddo
    enddo

    ! Average of the two derivative evaluations
    ! =========================================
    do k=ndz_e,nfz_e
       do j=ny-3,ny
          do i=nx-3,nx
             Krho(i,j,k)  = Krho(i,j,k)*0.5
             Krhou(i,j,k) = Krhou(i,j,k)*0.5
             Krhov(i,j,k) = Krhov(i,j,k)*0.5
             Krhow(i,j,k) = Krhow(i,j,k)*0.5
             Krhoe(i,j,k) = Krhoe(i,j,k)*0.5
          enddo
       enddo
    enddo

  end subroutine flux_euler_w_imax_jmax_9pts_c

end submodule smod_flux_euler_wedge_9pts_c
