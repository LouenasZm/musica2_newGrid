!=================================================================================
submodule (mod_flux_euler) smod_flux_euler_wedge_11pts_SBP4_c
!=================================================================================
  !> author: XG
  !> date: January 2024
  !> 2.5D curvilinear version - compute Eulerian fluxes (inviscid part) with 11-pt stencil
  !> using SBP4 boundary schemes
  !> Add_on routine to treat double corners  (wedge)
!=================================================================================

contains

  !===============================================================================
  module subroutine flux_euler_w_imin_jmin_11pts_SBP4_c
  !===============================================================================
    !> Derivatives of Euler fluxes for edge I_MIN/J_MIN close to two walls
    !> - 2.5D curvilinear coordinate - 11-point stencil - Summation by Parts o4 -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: vc,rovc
    ! ----------------------------------------------------------------------------

    ! ~> First pass: 11-pt metrics & order reduction for derivatives
    ! ==============================================================

    ! fluxes along ksi
    ! ----------------
    do k=ndz_e,nfz_e
       do j=1,5
          do i=-4,10
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
    ! Point #1: stencil [o x x x]
    !i=1
    do j=1,5
       do k=ndz_e,nfz_e
          Krho(1,j,k) = as4p0(1)*Frho(1,j,k) +as4p0(2)*Frho(2,j,k) &
                      + as4p0(3)*Frho(3,j,k) +as4p0(4)*Frho(4,j,k)
          Krhou(1,j,k)= as4p0(1)*Frhou(1,j,k)+as4p0(2)*Frhou(2,j,k) &
                      + as4p0(3)*Frhou(3,j,k)+as4p0(4)*Frhou(4,j,k)
          Krhov(1,j,k)= as4p0(1)*Frhov(1,j,k)+as4p0(2)*Frhov(2,j,k) &
                      + as4p0(3)*Frhov(3,j,k)+as4p0(4)*Frhov(4,j,k)
          Krhow(1,j,k)= as4p0(1)*Frhow(1,j,k)+as4p0(2)*Frhow(2,j,k) &
                      + as4p0(3)*Frhow(3,j,k)+as4p0(4)*Frhow(4,j,k)
          Krhoe(1,j,k)= as4p0(1)*Frhoe(1,j,k)+as4p0(2)*Frhoe(2,j,k) &
                      + as4p0(3)*Frhoe(3,j,k)+as4p0(4)*Frhoe(4,j,k)
       enddo
    enddo

    ! Point #2: stencil [x o x x]
    i=2
    do j=1,5
       do k=ndz_e,nfz_e
          Krho(2,j,k) = as4p1(1)*Frho(1,j,k) +as4p1(2)*Frho(2,j,k) &
                      + as4p1(3)*Frho(3,j,k) +as4p1(4)*Frho(4,j,k)
          Krhou(2,j,k)= as4p1(1)*Frhou(1,j,k)+as4p1(2)*Frhou(2,j,k) &
                      + as4p1(3)*Frhou(3,j,k)+as4p1(4)*Frhou(4,j,k)
          Krhov(2,j,k)= as4p1(1)*Frhov(1,j,k)+as4p1(2)*Frhov(2,j,k) &
                      + as4p1(3)*Frhov(3,j,k)+as4p1(4)*Frhov(4,j,k)
          Krhow(2,j,k)= as4p1(1)*Frhow(1,j,k)+as4p1(2)*Frhow(2,j,k) &
                      + as4p1(3)*Frhow(3,j,k)+as4p1(4)*Frhow(4,j,k)
          Krhoe(2,j,k)= as4p1(1)*Frhoe(1,j,k)+as4p1(2)*Frhoe(2,j,k) &
                      + as4p1(3)*Frhoe(3,j,k)+as4p1(4)*Frhoe(4,j,k)
       enddo
    enddo

    ! Point #3: stencil [x x o x x]
    !i=3
    do j=1,5
       do k=ndz_e,nfz_e
          Krho(3,j,k) = as4p2(1)*Frho(1,j,k) +as4p2(2)*Frho(2,j,k) &
                      + as4p2(3)*Frho(3,j,k) +as4p2(4)*Frho(4,j,k) &
                      + as4p2(5)*Frho(5,j,k)
          Krhou(3,j,k)= as4p2(1)*Frhou(1,j,k)+as4p2(2)*Frhou(2,j,k) &
                      + as4p2(3)*Frhou(3,j,k)+as4p2(4)*Frhou(4,j,k) &
                      + as4p2(5)*Frhou(5,j,k)
          Krhov(3,j,k)= as4p2(1)*Frhov(1,j,k)+as4p2(2)*Frhov(2,j,k) &
                      + as4p2(3)*Frhov(3,j,k)+as4p2(4)*Frhov(4,j,k) &
                      + as4p2(5)*Frhov(5,j,k)
          Krhow(3,j,k)= as4p2(1)*Frhow(1,j,k)+as4p2(2)*Frhow(2,j,k) &
                      + as4p2(3)*Frhow(3,j,k)+as4p2(4)*Frhow(4,j,k) &
                      + as4p2(5)*Frhow(5,j,k)
          Krhoe(3,j,k)= as4p2(1)*Frhoe(1,j,k)+as4p2(2)*Frhoe(2,j,k) &
                      + as4p2(3)*Frhoe(3,j,k)+as4p2(4)*Frhoe(4,j,k) &
                      + as4p2(5)*Frhoe(5,j,k)
       enddo
    enddo

    ! Point #4: stencil [x x x o x x]
    !i=4
    do j=1,5
       do k=ndz_e,nfz_e
          Krho(4,j,k) = as4p3(1)*Frho(1,j,k) +as4p3(2)*Frho(2,j,k) &
                      + as4p3(3)*Frho(3,j,k) +as4p3(4)*Frho(4,j,k) &
                      + as4p3(5)*Frho(5,j,k) +as4p3(6)*Frho(6,j,k)
          Krhou(4,j,k)= as4p3(1)*Frhou(1,j,k)+as4p3(2)*Frhou(2,j,k) &
                      + as4p3(3)*Frhou(3,j,k)+as4p3(4)*Frhou(4,j,k) &
                      + as4p3(5)*Frhou(5,j,k)+as4p3(6)*Frhou(6,j,k)
          Krhov(4,j,k)= as4p3(1)*Frhov(1,j,k)+as4p3(2)*Frhov(2,j,k) &
                      + as4p3(3)*Frhov(3,j,k)+as4p3(4)*Frhov(4,j,k) &
                      + as4p3(5)*Frhov(5,j,k)+as4p3(6)*Frhov(6,j,k)
          Krhow(4,j,k)= as4p3(1)*Frhow(1,j,k)+as4p3(2)*Frhow(2,j,k) &
                      + as4p3(3)*Frhow(3,j,k)+as4p3(4)*Frhow(4,j,k) &
                      + as4p3(5)*Frhow(5,j,k)+as4p3(6)*Frhow(6,j,k)
          Krhoe(4,j,k)= as4p3(1)*Frhoe(1,j,k)+as4p3(2)*Frhoe(2,j,k) &
                      + as4p3(3)*Frhoe(3,j,k)+as4p3(4)*Frhoe(4,j,k) &
                      + as4p3(5)*Frhoe(5,j,k)+as4p3(6)*Frhoe(6,j,k)
       enddo
    enddo

    ! Eighth-order (9-pt stencil)
    i=5
    do j=1,5
       do k=ndz_e,nfz_e
          Krho(i,j,k)  = a9(1) * ( Frho(i+1,j,k)-Frho(i-1,j,k) ) &
                       + a9(2) * ( Frho(i+2,j,k)-Frho(i-2,j,k) ) &
                       + a9(3) * ( Frho(i+3,j,k)-Frho(i-3,j,k) ) &
                       + a9(4) * ( Frho(i+4,j,k)-Frho(i-4,j,k) )
          Krhou(i,j,k) = a9(1) * ( Frhou(i+1,j,k)-Frhou(i-1,j,k) ) &
                       + a9(2) * ( Frhou(i+2,j,k)-Frhou(i-2,j,k) ) &
                       + a9(3) * ( Frhou(i+3,j,k)-Frhou(i-3,j,k) ) &
                       + a9(4) * ( Frhou(i+4,j,k)-Frhou(i-4,j,k) )
          Krhov(i,j,k) = a9(1) * ( Frhov(i+1,j,k)-Frhov(i-1,j,k) ) &
                       + a9(2) * ( Frhov(i+2,j,k)-Frhov(i-2,j,k) ) &
                       + a9(3) * ( Frhov(i+3,j,k)-Frhov(i-3,j,k) ) &
                       + a9(4) * ( Frhov(i+4,j,k)-Frhov(i-4,j,k) )
          Krhow(i,j,k) = a9(1) * ( Frhow(i+1,j,k)-Frhow(i-1,j,k) ) &
                       + a9(2) * ( Frhow(i+2,j,k)-Frhow(i-2,j,k) ) &
                       + a9(3) * ( Frhow(i+3,j,k)-Frhow(i-3,j,k) ) &
                       + a9(4) * ( Frhow(i+4,j,k)-Frhow(i-4,j,k) )
          Krhoe(i,j,k) = a9(1) * ( Frhoe(i+1,j,k)-Frhoe(i-1,j,k) ) &
                       + a9(2) * ( Frhoe(i+2,j,k)-Frhoe(i-2,j,k) ) &
                       + a9(3) * ( Frhoe(i+3,j,k)-Frhoe(i-3,j,k) ) &
                       + a9(4) * ( Frhoe(i+4,j,k)-Frhoe(i-4,j,k) )
       enddo
    enddo

    ! modified fluxes along eta
    ! -------------------------
    do k=ndz_e,nfz_e
       do i=1,5
          do j=-4,10
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
    ! Point #1: stencil [o x x x]
    !j=1
    do i=1,5
       do k=ndz_e,nfz_e  
          Krho(i,1,k) = as4p0(1)*Frho(i,1,k) +as4p0(2)*Frho(i,2,k) &
                      + as4p0(3)*Frho(i,3,k) +as4p0(4)*Frho(i,4,k) + Krho(i,1,k)
          Krhou(i,1,k)= as4p0(1)*Frhou(i,1,k)+as4p0(2)*Frhou(i,2,k) &
                      + as4p0(3)*Frhou(i,3,k)+as4p0(4)*Frhou(i,4,k) + Krhou(i,1,k)
          Krhov(i,1,k)= as4p0(1)*Frhov(i,1,k)+as4p0(2)*Frhov(i,2,k) &
                      + as4p0(3)*Frhov(i,3,k)+as4p0(4)*Frhov(i,4,k) + Krhov(i,1,k)
          Krhow(i,1,k)= as4p0(1)*Frhow(i,1,k)+as4p0(2)*Frhow(i,2,k) &
                      + as4p0(3)*Frhow(i,3,k)+as4p0(4)*Frhow(i,4,k) + Krhow(i,1,k)
          Krhoe(i,1,k)= as4p0(1)*Frhoe(i,1,k)+as4p0(2)*Frhoe(i,2,k) &
                      + as4p0(3)*Frhoe(i,3,k)+as4p0(4)*Frhoe(i,4,k) + Krhoe(i,1,k)
       enddo
    enddo

    ! Point #2: stencil [x o x x]
    !j=2
    do i=1,5
       do k=ndz_e,nfz_e  
          Krho(i,2,k) = as4p1(1)*Frho(i,1,k) +as4p1(2)*Frho(i,2,k) &
                      + as4p1(3)*Frho(i,3,k) +as4p1(4)*Frho(i,4,k) + Krho(i,2,k)
          Krhou(i,2,k)= as4p1(1)*Frhou(i,1,k)+as4p1(2)*Frhou(i,2,k) &
                      + as4p1(3)*Frhou(i,3,k)+as4p1(4)*Frhou(i,4,k) + Krhou(i,2,k)
          Krhov(i,2,k)= as4p1(1)*Frhov(i,1,k)+as4p1(2)*Frhov(i,2,k) &
                      + as4p1(3)*Frhov(i,3,k)+as4p1(4)*Frhov(i,4,k) + Krhov(i,2,k)
          Krhow(i,2,k)= as4p1(1)*Frhow(i,1,k)+as4p1(2)*Frhow(i,2,k) &
                      + as4p1(3)*Frhow(i,3,k)+as4p1(4)*Frhow(i,4,k) + Krhow(i,2,k)
          Krhoe(i,2,k)= as4p1(1)*Frhoe(i,1,k)+as4p1(2)*Frhoe(i,2,k) &
                      + as4p1(3)*Frhoe(i,3,k)+as4p1(4)*Frhoe(i,4,k) + Krhoe(i,2,k)
       enddo
    enddo

    ! Point #3: stencil [x x o x x]
    !j=3
    do i=1,5
       do k=ndz_e,nfz_e
          Krho(i,3,k) = as4p2(1)*Frho(i,1,k) +as4p2(2)*Frho(i,2,k) &
                      + as4p2(3)*Frho(i,3,k) +as4p2(4)*Frho(i,4,k) &
                      + as4p2(5)*Frho(i,5,k) + Krho(i,3,k)
          Krhou(i,3,k)= as4p2(1)*Frhou(i,1,k)+as4p2(2)*Frhou(i,2,k) &
                      + as4p2(3)*Frhou(i,3,k)+as4p2(4)*Frhou(i,4,k) &
                      + as4p2(5)*Frhou(i,5,k) + Krhou(i,3,k)
          Krhov(i,3,k)= as4p2(1)*Frhov(i,1,k)+as4p2(2)*Frhov(i,2,k) &
                      + as4p2(3)*Frhov(i,3,k)+as4p2(4)*Frhov(i,4,k) &
                      + as4p2(5)*Frhov(i,5,k) + Krhov(i,3,k)
          Krhow(i,3,k)= as4p2(1)*Frhow(i,1,k)+as4p2(2)*Frhow(i,2,k) &
                      + as4p2(3)*Frhow(i,3,k)+as4p2(4)*Frhow(i,4,k) &
                      + as4p2(5)*Frhow(i,5,k) + Krhow(i,3,k)
          Krhoe(i,3,k)= as4p2(1)*Frhoe(i,1,k)+as4p2(2)*Frhoe(i,2,k) &
                      + as4p2(3)*Frhoe(i,3,k)+as4p2(4)*Frhoe(i,4,k) &
                      + as4p2(5)*Frhoe(i,5,k) + Krhoe(i,3,k)
       enddo
    enddo

    ! Point #4: stencil [x x x o x x]
    !j=4
    do i=1,5
       do k=ndz_e,nfz_e
          Krho(i,4,k) = as4p3(1)*Frho(i,1,k) +as4p3(2)*Frho(i,2,k) &
                      + as4p3(3)*Frho(i,3,k) +as4p3(4)*Frho(i,4,k) &
                      + as4p3(5)*Frho(i,5,k) +as4p3(6)*Frho(i,6,k) + Krho(i,4,k)
          Krhou(i,4,k)= as4p3(1)*Frhou(i,1,k)+as4p3(2)*Frhou(i,2,k) &
                      + as4p3(3)*Frhou(i,3,k)+as4p3(4)*Frhou(i,4,k) &
                      + as4p3(5)*Frhou(i,5,k)+as4p3(6)*Frhou(i,6,k) + Krhou(i,4,k)
          Krhov(i,4,k)= as4p3(1)*Frhov(i,1,k)+as4p3(2)*Frhov(i,2,k) &
                      + as4p3(3)*Frhov(i,3,k)+as4p3(4)*Frhov(i,4,k) &
                      + as4p3(5)*Frhov(i,5,k)+as4p3(6)*Frhov(i,6,k) + Krhov(i,4,k)
          Krhow(i,4,k)= as4p3(1)*Frhow(i,1,k)+as4p3(2)*Frhow(i,2,k) &
                      + as4p3(3)*Frhow(i,3,k)+as4p3(4)*Frhow(i,4,k) &
                      + as4p3(5)*Frhow(i,5,k)+as4p3(6)*Frhow(i,6,k) + Krhow(i,4,k)
          Krhoe(i,4,k)= as4p3(1)*Frhoe(i,1,k)+as4p3(2)*Frhoe(i,2,k) &
                      + as4p3(3)*Frhoe(i,3,k)+as4p3(4)*Frhoe(i,4,k) &
                      + as4p3(5)*Frhoe(i,5,k)+as4p3(6)*Frhoe(i,6,k) + Krhoe(i,4,k)
       enddo
    enddo

    ! Eighth-order (9-pt stencil)
    j=5
    do i=1,5
       do k=ndz_e,nfz_e
          Krho(i,j,k)  = a9(1) * ( Frho(i,j+1,k)-Frho(i,j-1,k) ) &
                       + a9(2) * ( Frho(i,j+2,k)-Frho(i,j-2,k) ) &
                       + a9(3) * ( Frho(i,j+3,k)-Frho(i,j-3,k) ) &
                       + a9(4) * ( Frho(i,j+4,k)-Frho(i,j-4,k) ) + Krho(i,j,k)      
          Krhou(i,j,k) = a9(1) * ( Frhou(i,j+1,k)-Frhou(i,j-1,k) ) &
                       + a9(2) * ( Frhou(i,j+2,k)-Frhou(i,j-2,k) ) &
                       + a9(3) * ( Frhou(i,j+3,k)-Frhou(i,j-3,k) ) &
                       + a9(4) * ( Frhou(i,j+4,k)-Frhou(i,j-4,k) ) + Krhou(i,j,k)
          Krhov(i,j,k) = a9(1) * ( Frhov(i,j+1,k)-Frhov(i,j-1,k) ) &
                       + a9(2) * ( Frhov(i,j+2,k)-Frhov(i,j-2,k) ) &
                       + a9(3) * ( Frhov(i,j+3,k)-Frhov(i,j-3,k) ) &
                       + a9(4) * ( Frhov(i,j+4,k)-Frhov(i,j-4,k) ) + Krhov(i,j,k)
          Krhow(i,j,k) = a9(1) * ( Frhow(i,j+1,k)-Frhow(i,j-1,k) ) &
                       + a9(2) * ( Frhow(i,j+2,k)-Frhow(i,j-2,k) ) &
                       + a9(3) * ( Frhow(i,j+3,k)-Frhow(i,j-3,k) ) &
                       + a9(4) * ( Frhow(i,j+4,k)-Frhow(i,j-4,k) ) + Krhow(i,j,k)
          Krhoe(i,j,k) = a9(1) * ( Frhoe(i,j+1,k)-Frhoe(i,j-1,k) ) &
                       + a9(2) * ( Frhoe(i,j+2,k)-Frhoe(i,j-2,k) ) &
                       + a9(3) * ( Frhoe(i,j+3,k)-Frhoe(i,j-3,k) ) &
                       + a9(4) * ( Frhoe(i,j+4,k)-Frhoe(i,j-4,k) ) + Krhoe(i,j,k)
       enddo
    enddo

    ! ~> Second pass: order reduction for metrics & 11-pt derivatives
    ! ===============================================================

    ! modified fluxes along ksi
    ! -------------------------
    do k=ndz_e,nfz_e
       do j=1,5
          do i=-4,10
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
       do j=1,5
          do i=1,5
             Krho(i,j,k)  = a11(1) * ( Frho(i+1,j,k)-Frho(i-1,j,k) )   &
                          + a11(2) * ( Frho(i+2,j,k)-Frho(i-2,j,k) )   &
                          + a11(3) * ( Frho(i+3,j,k)-Frho(i-3,j,k) )   &
                          + a11(4) * ( Frho(i+4,j,k)-Frho(i-4,j,k) )   &
                          + a11(5) * ( Frho(i+5,j,k)-Frho(i-5,j,k) ) + Krho(i,j,k)

             Krhou(i,j,k) = a11(1) * ( Frhou(i+1,j,k)-Frhou(i-1,j,k) ) &
                          + a11(2) * ( Frhou(i+2,j,k)-Frhou(i-2,j,k) ) &
                          + a11(3) * ( Frhou(i+3,j,k)-Frhou(i-3,j,k) ) &
                          + a11(4) * ( Frhou(i+4,j,k)-Frhou(i-4,j,k) ) &
                          + a11(5) * ( Frhou(i+5,j,k)-Frhou(i-5,j,k) ) + Krhou(i,j,k)

             Krhov(i,j,k) = a11(1) * ( Frhov(i+1,j,k)-Frhov(i-1,j,k) ) &
                          + a11(2) * ( Frhov(i+2,j,k)-Frhov(i-2,j,k) ) &
                          + a11(3) * ( Frhov(i+3,j,k)-Frhov(i-3,j,k) ) &
                          + a11(4) * ( Frhov(i+4,j,k)-Frhov(i-4,j,k) ) &
                          + a11(5) * ( Frhov(i+5,j,k)-Frhov(i-5,j,k) ) + Krhov(i,j,k)

             Krhow(i,j,k) = a11(1) * ( Frhow(i+1,j,k)-Frhow(i-1,j,k) ) &
                          + a11(2) * ( Frhow(i+2,j,k)-Frhow(i-2,j,k) ) &
                          + a11(3) * ( Frhow(i+3,j,k)-Frhow(i-3,j,k) ) &
                          + a11(4) * ( Frhow(i+4,j,k)-Frhow(i-4,j,k) ) &
                          + a11(5) * ( Frhow(i+5,j,k)-Frhow(i-5,j,k) ) + Krhow(i,j,k)

             Krhoe(i,j,k) = a11(1) * ( Frhoe(i+1,j,k)-Frhoe(i-1,j,k) ) &
                          + a11(2) * ( Frhoe(i+2,j,k)-Frhoe(i-2,j,k) ) &
                          + a11(3) * ( Frhoe(i+3,j,k)-Frhoe(i-3,j,k) ) &
                          + a11(4) * ( Frhoe(i+4,j,k)-Frhoe(i-4,j,k) ) &
                          + a11(5) * ( Frhoe(i+5,j,k)-Frhoe(i-5,j,k) ) + Krhoe(i,j,k)
          enddo
       enddo
    enddo

    ! modified fluxes along eta
    ! -------------------------
    do k=ndz_e,nfz_e
       do i=1,5
          do j=-4,10
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
       do j=1,5
          do i=1,5
             Krho(i,j,k)  = a11(1) * ( Frho(i,j+1,k)-Frho(i,j-1,k) ) &
                          + a11(2) * ( Frho(i,j+2,k)-Frho(i,j-2,k) ) &
                          + a11(3) * ( Frho(i,j+3,k)-Frho(i,j-3,k) ) &
                          + a11(4) * ( Frho(i,j+4,k)-Frho(i,j-4,k) ) &
                          + a11(5) * ( Frho(i,j+5,k)-Frho(i,j-5,k) ) &
                          + Krho(i,j,k)

             Krhou(i,j,k) = a11(1) * ( Frhou(i,j+1,k)-Frhou(i,j-1,k) ) &
                          + a11(2) * ( Frhou(i,j+2,k)-Frhou(i,j-2,k) ) &
                          + a11(3) * ( Frhou(i,j+3,k)-Frhou(i,j-3,k) ) &
                          + a11(4) * ( Frhou(i,j+4,k)-Frhou(i,j-4,k) ) &
                          + a11(5) * ( Frhou(i,j+5,k)-Frhou(i,j-5,k) ) &
                          + Krhou(i,j,k)

             Krhov(i,j,k) = a11(1) * ( Frhov(i,j+1,k)-Frhov(i,j-1,k) ) &
                          + a11(2) * ( Frhov(i,j+2,k)-Frhov(i,j-2,k) ) &
                          + a11(3) * ( Frhov(i,j+3,k)-Frhov(i,j-3,k) ) &
                          + a11(4) * ( Frhov(i,j+4,k)-Frhov(i,j-4,k) ) &
                          + a11(5) * ( Frhov(i,j+5,k)-Frhov(i,j-5,k) ) &
                          + Krhov(i,j,k)

             Krhow(i,j,k) = a11(1) * ( Frhow(i,j+1,k)-Frhow(i,j-1,k) ) &
                          + a11(2) * ( Frhow(i,j+2,k)-Frhow(i,j-2,k) ) &
                          + a11(3) * ( Frhow(i,j+3,k)-Frhow(i,j-3,k) ) &
                          + a11(4) * ( Frhow(i,j+4,k)-Frhow(i,j-4,k) ) &
                          + a11(5) * ( Frhow(i,j+5,k)-Frhow(i,j-5,k) ) &
                          + Krhow(i,j,k)

             Krhoe(i,j,k) = a11(1) * ( Frhoe(i,j+1,k)-Frhoe(i,j-1,k) ) &
                          + a11(2) * ( Frhoe(i,j+2,k)-Frhoe(i,j-2,k) ) &
                          + a11(3) * ( Frhoe(i,j+3,k)-Frhoe(i,j-3,k) ) &
                          + a11(4) * ( Frhoe(i,j+4,k)-Frhoe(i,j-4,k) ) &
                          + a11(5) * ( Frhoe(i,j+5,k)-Frhoe(i,j-5,k) ) &
                          + Krhoe(i,j,k)
          enddo
       enddo
    enddo

    ! Average of the two derivative evaluations
    ! =========================================
    do k=ndz_e,nfz_e
       do j=1,5
          do i=1,5
             Krho(i,j,k)  = Krho(i,j,k)*0.5_wp
             Krhou(i,j,k) = Krhou(i,j,k)*0.5_wp
             Krhov(i,j,k) = Krhov(i,j,k)*0.5_wp
             Krhow(i,j,k) = Krhow(i,j,k)*0.5_wp
             Krhoe(i,j,k) = Krhoe(i,j,k)*0.5_wp
          enddo
       enddo
    enddo

  end subroutine flux_euler_w_imin_jmin_11pts_SBP4_c

  !===============================================================================
  module subroutine flux_euler_w_imin_jmax_11pts_SBP4_c
  !===============================================================================
    !> Derivatives of Euler fluxes for edge I_MIN/J_MAX close to two walls
    !> - curvilinear coordinate - 11-point stencil - Summation by Parts o4 -
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: vc,rovc
    ! ---------------------------------------------------------------------------

    ! ~> First pass: 11-pt metrics & order reduction for derivatives
    ! ==============================================================

    ! fluxes along ksi
    ! ----------------
    do k=ndz_e,nfz_e
       do j=ny-4,ny
          do i=-4,10
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
    ! Point #1: stencil [o x x x]
    !i=1
    do j=ny-4,ny
       do k=ndz_e,nfz_e
          Krho(1,j,k) = as4p0(1)*Frho(1,j,k) +as4p0(2)*Frho(2,j,k) &
                      + as4p0(3)*Frho(3,j,k) +as4p0(4)*Frho(4,j,k)
          Krhou(1,j,k)= as4p0(1)*Frhou(1,j,k)+as4p0(2)*Frhou(2,j,k) &
                      + as4p0(3)*Frhou(3,j,k)+as4p0(4)*Frhou(4,j,k)
          Krhov(1,j,k)= as4p0(1)*Frhov(1,j,k)+as4p0(2)*Frhov(2,j,k) &
                      + as4p0(3)*Frhov(3,j,k)+as4p0(4)*Frhov(4,j,k)
          Krhow(1,j,k)= as4p0(1)*Frhow(1,j,k)+as4p0(2)*Frhow(2,j,k) &
                      + as4p0(3)*Frhow(3,j,k)+as4p0(4)*Frhow(4,j,k)
          Krhoe(1,j,k)= as4p0(1)*Frhoe(1,j,k)+as4p0(2)*Frhoe(2,j,k) &
                      + as4p0(3)*Frhoe(3,j,k)+as4p0(4)*Frhoe(4,j,k)
       enddo
    enddo

    ! Point #2: stencil [x o x x]
    !i=2
    do j=ny-4,ny
       do k=ndz_e,nfz_e
          Krho(2,j,k) = as4p1(1)*Frho(1,j,k) +as4p1(2)*Frho(2,j,k) &
                      + as4p1(3)*Frho(3,j,k) +as4p1(4)*Frho(4,j,k)
          Krhou(2,j,k)= as4p1(1)*Frhou(1,j,k)+as4p1(2)*Frhou(2,j,k) &
                      + as4p1(3)*Frhou(3,j,k)+as4p1(4)*Frhou(4,j,k)
          Krhov(2,j,k)= as4p1(1)*Frhov(1,j,k)+as4p1(2)*Frhov(2,j,k) &
                      + as4p1(3)*Frhov(3,j,k)+as4p1(4)*Frhov(4,j,k)
          Krhow(2,j,k)= as4p1(1)*Frhow(1,j,k)+as4p1(2)*Frhow(2,j,k) &
                      + as4p1(3)*Frhow(3,j,k)+as4p1(4)*Frhow(4,j,k)
          Krhoe(2,j,k)= as4p1(1)*Frhoe(1,j,k)+as4p1(2)*Frhoe(2,j,k) &
                      + as4p1(3)*Frhoe(3,j,k)+as4p1(4)*Frhoe(4,j,k)
       enddo
    enddo

    ! Point #3: stencil [x x o x x]
    !i=3
    do j=ny-4,ny
       do k=ndz_e,nfz_e
          Krho(3,j,k) = as4p2(1)*Frho(1,j,k) +as4p2(2)*Frho(2,j,k) &
                      + as4p2(3)*Frho(3,j,k) +as4p2(4)*Frho(4,j,k) &
                      + as4p2(5)*Frho(5,j,k)
          Krhou(3,j,k)= as4p2(1)*Frhou(1,j,k)+as4p2(2)*Frhou(2,j,k) &
                      + as4p2(3)*Frhou(3,j,k)+as4p2(4)*Frhou(4,j,k) &
                      + as4p2(5)*Frhou(5,j,k)
          Krhov(3,j,k)= as4p2(1)*Frhov(1,j,k)+as4p2(2)*Frhov(2,j,k) &
                      + as4p2(3)*Frhov(3,j,k)+as4p2(4)*Frhov(4,j,k) &
                      + as4p2(5)*Frhov(5,j,k)
          Krhow(3,j,k)= as4p2(1)*Frhow(1,j,k)+as4p2(2)*Frhow(2,j,k) &
                      + as4p2(3)*Frhow(3,j,k)+as4p2(4)*Frhow(4,j,k) &
                      + as4p2(5)*Frhow(5,j,k)
          Krhoe(3,j,k)= as4p2(1)*Frhoe(1,j,k)+as4p2(2)*Frhoe(2,j,k) &
                      + as4p2(3)*Frhoe(3,j,k)+as4p2(4)*Frhoe(4,j,k) &
                      + as4p2(5)*Frhoe(5,j,k)
       enddo
    enddo

    ! Point #4: stencil [x x x o x x]
    !i=4
    do j=ny-4,ny
       do k=ndz_e,nfz_e
          Krho(4,j,k) = as4p3(1)*Frho(1,j,k) +as4p3(2)*Frho(2,j,k) &
                      + as4p3(3)*Frho(3,j,k) +as4p3(4)*Frho(4,j,k) &
                      + as4p3(5)*Frho(5,j,k) +as4p3(6)*Frho(6,j,k)
          Krhou(4,j,k)= as4p3(1)*Frhou(1,j,k)+as4p3(2)*Frhou(2,j,k) &
                      + as4p3(3)*Frhou(3,j,k)+as4p3(4)*Frhou(4,j,k) &
                      + as4p3(5)*Frhou(5,j,k)+as4p3(6)*Frhou(6,j,k)
          Krhov(4,j,k)= as4p3(1)*Frhov(1,j,k)+as4p3(2)*Frhov(2,j,k) &
                      + as4p3(3)*Frhov(3,j,k)+as4p3(4)*Frhov(4,j,k) &
                      + as4p3(5)*Frhov(5,j,k)+as4p3(6)*Frhov(6,j,k)
          Krhow(4,j,k)= as4p3(1)*Frhow(1,j,k)+as4p3(2)*Frhow(2,j,k) &
                      + as4p3(3)*Frhow(3,j,k)+as4p3(4)*Frhow(4,j,k) &
                      + as4p3(5)*Frhow(5,j,k)+as4p3(6)*Frhow(6,j,k)
          Krhoe(4,j,k)= as4p3(1)*Frhoe(1,j,k)+as4p3(2)*Frhoe(2,j,k) &
                      + as4p3(3)*Frhoe(3,j,k)+as4p3(4)*Frhoe(4,j,k) &
                      + as4p3(5)*Frhoe(5,j,k)+as4p3(6)*Frhoe(6,j,k)
       enddo
    enddo

    ! Eighth-order (9-pt stencil)
    i=5
    do j=ny-4,ny
       do k=ndz_e,nfz_e
          Krho(i,j,k)  = a9(1) * ( Frho(i+1,j,k)-Frho(i-1,j,k) ) &
                       + a9(2) * ( Frho(i+2,j,k)-Frho(i-2,j,k) ) &
                       + a9(3) * ( Frho(i+3,j,k)-Frho(i-3,j,k) ) &
                       + a9(4) * ( Frho(i+4,j,k)-Frho(i-4,j,k) )
          Krhou(i,j,k) = a9(1) * ( Frhou(i+1,j,k)-Frhou(i-1,j,k) ) &
                       + a9(2) * ( Frhou(i+2,j,k)-Frhou(i-2,j,k) ) &
                       + a9(3) * ( Frhou(i+3,j,k)-Frhou(i-3,j,k) ) &
                       + a9(4) * ( Frhou(i+4,j,k)-Frhou(i-4,j,k) )
          Krhov(i,j,k) = a9(1) * ( Frhov(i+1,j,k)-Frhov(i-1,j,k) ) &
                       + a9(2) * ( Frhov(i+2,j,k)-Frhov(i-2,j,k) ) &
                       + a9(3) * ( Frhov(i+3,j,k)-Frhov(i-3,j,k) ) &
                       + a9(4) * ( Frhov(i+4,j,k)-Frhov(i-4,j,k) )
          Krhow(i,j,k) = a9(1) * ( Frhow(i+1,j,k)-Frhow(i-1,j,k) ) &
                       + a9(2) * ( Frhow(i+2,j,k)-Frhow(i-2,j,k) ) &
                       + a9(3) * ( Frhow(i+3,j,k)-Frhow(i-3,j,k) ) &
                       + a9(4) * ( Frhow(i+4,j,k)-Frhow(i-4,j,k) )
          Krhoe(i,j,k) = a9(1) * ( Frhoe(i+1,j,k)-Frhoe(i-1,j,k) ) &
                       + a9(2) * ( Frhoe(i+2,j,k)-Frhoe(i-2,j,k) ) &
                       + a9(3) * ( Frhoe(i+3,j,k)-Frhoe(i-3,j,k) ) &
                       + a9(4) * ( Frhoe(i+4,j,k)-Frhoe(i-4,j,k) )
       enddo
    enddo

    ! modified fluxes along eta
    ! -------------------------
    do k=ndz_e,nfz_e
       do i=1,5
          do j=ny-9,ny+5
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
    ! Point #1: stencil [o x x x]
    j=ny
    do i=1,5
       do k=ndz_e,nfz_e  
          Krho(i,j,k) = as4m0(1)*Frho(i,ny  ,k) +as4m0(2)*Frho(i,ny-1,k) &
                      + as4m0(3)*Frho(i,ny-2,k) +as4m0(4)*Frho(i,ny-3,k) + Krho(i,j,k)
          Krhou(i,j,k)= as4m0(1)*Frhou(i,ny  ,k)+as4m0(2)*Frhou(i,ny-1,k) &
                      + as4m0(3)*Frhou(i,ny-2,k)+as4m0(4)*Frhou(i,ny-3,k) + Krhou(i,j,k)
          Krhov(i,j,k)= as4m0(1)*Frhov(i,ny  ,k)+as4m0(2)*Frhov(i,ny-1,k) &
                      + as4m0(3)*Frhov(i,ny-2,k)+as4m0(4)*Frhov(i,ny-3,k) + Krhov(i,j,k)
          Krhow(i,j,k)= as4m0(1)*Frhow(i,ny  ,k)+as4m0(2)*Frhow(i,ny-1,k) &
                      + as4m0(3)*Frhow(i,ny-2,k)+as4m0(4)*Frhow(i,ny-3,k) + Krhow(i,j,k)
          Krhoe(i,j,k)= as4m0(1)*Frhoe(i,ny  ,k)+as4m0(2)*Frhoe(i,ny-1,k) &
                      + as4m0(3)*Frhoe(i,ny-2,k)+as4m0(4)*Frhoe(i,ny-3,k) + Krhoe(i,j,k)
       enddo
    enddo

    ! Point #2: stencil [x o x x]
    j=ny-1
    do i=1,5
       do k=ndz_e,nfz_e  
          Krho(i,j,k) = as4m1(1)*Frho(i,ny  ,k) +as4m1(2)*Frho(i,ny-1,k) &
                      + as4m1(3)*Frho(i,ny-2,k) +as4m1(4)*Frho(i,ny-3,k) + Krho(i,j,k)
          Krhou(i,j,k)= as4m1(1)*Frhou(i,ny  ,k)+as4m1(2)*Frhou(i,ny-1,k) &
                      + as4m1(3)*Frhou(i,ny-2,k)+as4m1(4)*Frhou(i,ny-3,k) + Krhou(i,j,k)
          Krhov(i,j,k)= as4m1(1)*Frhov(i,ny  ,k)+as4m1(2)*Frhov(i,ny-1,k) &
                      + as4m1(3)*Frhov(i,ny-2,k)+as4m1(4)*Frhov(i,ny-3,k) + Krhov(i,j,k)
          Krhow(i,j,k)= as4m1(1)*Frhow(i,ny  ,k)+as4m1(2)*Frhow(i,ny-1,k) &
                      + as4m1(3)*Frhow(i,ny-2,k)+as4m1(4)*Frhow(i,ny-3,k) + Krhow(i,j,k)
          Krhoe(i,j,k)= as4m1(1)*Frhoe(i,ny  ,k)+as4m1(2)*Frhoe(i,ny-1,k) &
                      + as4m1(3)*Frhoe(i,ny-2,k)+as4m1(4)*Frhoe(i,ny-3,k) + Krhoe(i,j,k)
       enddo
    enddo

    ! Point #3: stencil [x x o x x]
    j=ny-2
    do i=1,5
       do k=ndz_e,nfz_e
          Krho(i,j,k) = as4m2(1)*Frho(i,ny  ,k) +as4m2(2)*Frho(i,ny-1,k) &
                      + as4m2(3)*Frho(i,ny-2,k) +as4m2(4)*Frho(i,ny-3,k) &
                      + as4m2(5)*Frho(i,ny-4,k) + Krho(i,j,k)
          Krhou(i,j,k)= as4m2(1)*Frhou(i,ny  ,k)+as4m2(2)*Frhou(i,ny-1,k) &
                      + as4m2(3)*Frhou(i,ny-2,k)+as4m2(4)*Frhou(i,ny-3,k) &
                      + as4m2(5)*Frhou(i,ny-4,k)+ Krhou(i,j,k)
          Krhov(i,j,k)= as4m2(1)*Frhov(i,ny  ,k)+as4m2(2)*Frhov(i,ny-1,k) &
                      + as4m2(3)*Frhov(i,ny-2,k)+as4m2(4)*Frhov(i,ny-3,k) &
                      + as4m2(5)*Frhov(i,ny-4,k)+ Krhov(i,j,k)
          Krhow(i,j,k)= as4m2(1)*Frhow(i,ny  ,k)+as4m2(2)*Frhow(i,ny-1,k) &
                      + as4m2(3)*Frhow(i,ny-2,k)+as4m2(4)*Frhow(i,ny-3,k) &
                      + as4m2(5)*Frhow(i,ny-4,k)+ Krhow(i,j,k)
          Krhoe(i,j,k)= as4m2(1)*Frhoe(i,ny  ,k)+as4m2(2)*Frhoe(i,ny-1,k) &
                      + as4m2(3)*Frhoe(i,ny-2,k)+as4m2(4)*Frhoe(i,ny-3,k) &
                      + as4m2(5)*Frhoe(i,ny-4,k)+ Krhoe(i,j,k)
       enddo
    enddo

    ! Point #4: stencil [x x x o x x]
    j=ny-3
    do i=1,5
       do k=ndz_e,nfz_e
          Krho(i,j,k) = as4m3(1)*Frho(i,ny  ,k) +as4m3(2)*Frho(i,ny-1,k) &
                      + as4m3(3)*Frho(i,ny-2,k) +as4m3(4)*Frho(i,ny-3,k) &
                      + as4m3(5)*Frho(i,ny-4,k) +as4m3(6)*Frho(i,ny-5,k) + Krho(i,j,k)
          Krhou(i,j,k)= as4m3(1)*Frhou(i,ny  ,k)+as4m3(2)*Frhou(i,ny-1,k) &
                      + as4m3(3)*Frhou(i,ny-2,k)+as4m3(4)*Frhou(i,ny-3,k) &
                      + as4m3(5)*Frhou(i,ny-4,k)+as4m3(6)*Frhou(i,ny-5,k) + Krhou(i,j,k)
          Krhov(i,j,k)= as4m3(1)*Frhov(i,ny  ,k)+as4m3(2)*Frhov(i,ny-1,k) &
                      + as4m3(3)*Frhov(i,ny-2,k)+as4m3(4)*Frhov(i,ny-3,k) &
                      + as4m3(5)*Frhov(i,ny-4,k)+as4m3(6)*Frhov(i,ny-5,k) + Krhov(i,j,k)
          Krhow(i,j,k)= as4m3(1)*Frhow(i,ny  ,k)+as4m3(2)*Frhow(i,ny-1,k) &
                      + as4m3(3)*Frhow(i,ny-2,k)+as4m3(4)*Frhow(i,ny-3,k) &
                      + as4m3(5)*Frhow(i,ny-4,k)+as4m3(6)*Frhow(i,ny-5,k) + Krhow(i,j,k)
          Krhoe(i,j,k)= as4m3(1)*Frhoe(i,ny  ,k)+as4m3(2)*Frhoe(i,ny-1,k) &
                      + as4m3(3)*Frhoe(i,ny-2,k)+as4m3(4)*Frhoe(i,ny-3,k) &
                      + as4m3(5)*Frhoe(i,ny-4,k)+as4m3(6)*Frhoe(i,ny-5,k) + Krhoe(i,j,k)
       enddo
    enddo

    ! Eighth-order (9-pt stencil)
    j=ny-4
    do i=1,5
       do k=ndz_e,nfz_e
          Krho(i,j,k)  = a9(1) * ( Frho(i,j+1,k)-Frho(i,j-1,k) ) &
                       + a9(2) * ( Frho(i,j+2,k)-Frho(i,j-2,k) ) &
                       + a9(3) * ( Frho(i,j+3,k)-Frho(i,j-3,k) ) &
                       + a9(4) * ( Frho(i,j+4,k)-Frho(i,j-4,k) ) + Krho(i,j,k)      
          Krhou(i,j,k) = a9(1) * ( Frhou(i,j+1,k)-Frhou(i,j-1,k) ) &
                       + a9(2) * ( Frhou(i,j+2,k)-Frhou(i,j-2,k) ) &
                       + a9(3) * ( Frhou(i,j+3,k)-Frhou(i,j-3,k) ) &
                       + a9(4) * ( Frhou(i,j+4,k)-Frhou(i,j-4,k) ) + Krhou(i,j,k)
          Krhov(i,j,k) = a9(1) * ( Frhov(i,j+1,k)-Frhov(i,j-1,k) ) &
                       + a9(2) * ( Frhov(i,j+2,k)-Frhov(i,j-2,k) ) &
                       + a9(3) * ( Frhov(i,j+3,k)-Frhov(i,j-3,k) ) &
                       + a9(4) * ( Frhov(i,j+4,k)-Frhov(i,j-4,k) ) + Krhov(i,j,k)
          Krhow(i,j,k) = a9(1) * ( Frhow(i,j+1,k)-Frhow(i,j-1,k) ) &
                       + a9(2) * ( Frhow(i,j+2,k)-Frhow(i,j-2,k) ) &
                       + a9(3) * ( Frhow(i,j+3,k)-Frhow(i,j-3,k) ) &
                       + a9(4) * ( Frhow(i,j+4,k)-Frhow(i,j-4,k) ) + Krhow(i,j,k)
          Krhoe(i,j,k) = a9(1) * ( Frhoe(i,j+1,k)-Frhoe(i,j-1,k) ) &
                       + a9(2) * ( Frhoe(i,j+2,k)-Frhoe(i,j-2,k) ) &
                       + a9(3) * ( Frhoe(i,j+3,k)-Frhoe(i,j-3,k) ) &
                       + a9(4) * ( Frhoe(i,j+4,k)-Frhoe(i,j-4,k) ) + Krhoe(i,j,k)
       enddo
    enddo

    ! ~> Second pass: order reduction for metrics & 11-pt derivatives
    ! ===============================================================

    ! modified fluxes along ksi
    ! -------------------------
    do k=ndz_e,nfz_e
       do j=ny-4,ny
          do i=-4,10
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
       do j=ny-4,ny
          do i=1,5
             Krho(i,j,k)  = a11(1) * ( Frho(i+1,j,k)-Frho(i-1,j,k) )   &
                          + a11(2) * ( Frho(i+2,j,k)-Frho(i-2,j,k) )   &
                          + a11(3) * ( Frho(i+3,j,k)-Frho(i-3,j,k) )   &
                          + a11(4) * ( Frho(i+4,j,k)-Frho(i-4,j,k) )   &
                          + a11(5) * ( Frho(i+5,j,k)-Frho(i-5,j,k) ) + Krho(i,j,k)

             Krhou(i,j,k) = a11(1) * ( Frhou(i+1,j,k)-Frhou(i-1,j,k) ) &
                          + a11(2) * ( Frhou(i+2,j,k)-Frhou(i-2,j,k) ) &
                          + a11(3) * ( Frhou(i+3,j,k)-Frhou(i-3,j,k) ) &
                          + a11(4) * ( Frhou(i+4,j,k)-Frhou(i-4,j,k) ) &
                          + a11(5) * ( Frhou(i+5,j,k)-Frhou(i-5,j,k) ) + Krhou(i,j,k)

             Krhov(i,j,k) = a11(1) * ( Frhov(i+1,j,k)-Frhov(i-1,j,k) ) &
                          + a11(2) * ( Frhov(i+2,j,k)-Frhov(i-2,j,k) ) &
                          + a11(3) * ( Frhov(i+3,j,k)-Frhov(i-3,j,k) ) &
                          + a11(4) * ( Frhov(i+4,j,k)-Frhov(i-4,j,k) ) &
                          + a11(5) * ( Frhov(i+5,j,k)-Frhov(i-5,j,k) ) + Krhov(i,j,k)

             Krhow(i,j,k) = a11(1) * ( Frhow(i+1,j,k)-Frhow(i-1,j,k) ) &
                          + a11(2) * ( Frhow(i+2,j,k)-Frhow(i-2,j,k) ) &
                          + a11(3) * ( Frhow(i+3,j,k)-Frhow(i-3,j,k) ) &
                          + a11(4) * ( Frhow(i+4,j,k)-Frhow(i-4,j,k) ) &
                          + a11(5) * ( Frhow(i+5,j,k)-Frhow(i-5,j,k) ) + Krhow(i,j,k)

             Krhoe(i,j,k) = a11(1) * ( Frhoe(i+1,j,k)-Frhoe(i-1,j,k) ) &
                          + a11(2) * ( Frhoe(i+2,j,k)-Frhoe(i-2,j,k) ) &
                          + a11(3) * ( Frhoe(i+3,j,k)-Frhoe(i-3,j,k) ) &
                          + a11(4) * ( Frhoe(i+4,j,k)-Frhoe(i-4,j,k) ) &
                          + a11(5) * ( Frhoe(i+5,j,k)-Frhoe(i-5,j,k) ) + Krhoe(i,j,k)
          enddo
       enddo
    enddo

    ! modified fluxes along eta
    ! -------------------------
    do k=ndz_e,nfz_e
       do i=1,5
          do j=ny-9,ny+5
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
       do j=ny-4,ny
          do i=1,5
             Krho(i,j,k)  = a11(1) * ( Frho(i,j+1,k)-Frho(i,j-1,k) ) &
                          + a11(2) * ( Frho(i,j+2,k)-Frho(i,j-2,k) ) &
                          + a11(3) * ( Frho(i,j+3,k)-Frho(i,j-3,k) ) &
                          + a11(4) * ( Frho(i,j+4,k)-Frho(i,j-4,k) ) &
                          + a11(5) * ( Frho(i,j+5,k)-Frho(i,j-5,k) ) &
                          + Krho(i,j,k)

             Krhou(i,j,k) = a11(1) * ( Frhou(i,j+1,k)-Frhou(i,j-1,k) ) &
                          + a11(2) * ( Frhou(i,j+2,k)-Frhou(i,j-2,k) ) &
                          + a11(3) * ( Frhou(i,j+3,k)-Frhou(i,j-3,k) ) &
                          + a11(4) * ( Frhou(i,j+4,k)-Frhou(i,j-4,k) ) &
                          + a11(5) * ( Frhou(i,j+5,k)-Frhou(i,j-5,k) ) &
                          + Krhou(i,j,k)

             Krhov(i,j,k) = a11(1) * ( Frhov(i,j+1,k)-Frhov(i,j-1,k) ) &
                          + a11(2) * ( Frhov(i,j+2,k)-Frhov(i,j-2,k) ) &
                          + a11(3) * ( Frhov(i,j+3,k)-Frhov(i,j-3,k) ) &
                          + a11(4) * ( Frhov(i,j+4,k)-Frhov(i,j-4,k) ) &
                          + a11(5) * ( Frhov(i,j+5,k)-Frhov(i,j-5,k) ) &
                          + Krhov(i,j,k)

             Krhow(i,j,k) = a11(1) * ( Frhow(i,j+1,k)-Frhow(i,j-1,k) ) &
                          + a11(2) * ( Frhow(i,j+2,k)-Frhow(i,j-2,k) ) &
                          + a11(3) * ( Frhow(i,j+3,k)-Frhow(i,j-3,k) ) &
                          + a11(4) * ( Frhow(i,j+4,k)-Frhow(i,j-4,k) ) &
                          + a11(5) * ( Frhow(i,j+5,k)-Frhow(i,j-5,k) ) &
                          + Krhow(i,j,k)

             Krhoe(i,j,k) = a11(1) * ( Frhoe(i,j+1,k)-Frhoe(i,j-1,k) ) &
                          + a11(2) * ( Frhoe(i,j+2,k)-Frhoe(i,j-2,k) ) &
                          + a11(3) * ( Frhoe(i,j+3,k)-Frhoe(i,j-3,k) ) &
                          + a11(4) * ( Frhoe(i,j+4,k)-Frhoe(i,j-4,k) ) &
                          + a11(5) * ( Frhoe(i,j+5,k)-Frhoe(i,j-5,k) ) &
                          + Krhoe(i,j,k)
          enddo
       enddo
    enddo

    ! Average of the two derivative evaluations
    ! =========================================
    do k=ndz_e,nfz_e
       do j=ny-4,ny
          do i=1,5
             Krho(i,j,k)  = Krho(i,j,k)*0.5_wp
             Krhou(i,j,k) = Krhou(i,j,k)*0.5_wp
             Krhov(i,j,k) = Krhov(i,j,k)*0.5_wp
             Krhow(i,j,k) = Krhow(i,j,k)*0.5_wp
             Krhoe(i,j,k) = Krhoe(i,j,k)*0.5_wp
          enddo
       enddo
    enddo

  end subroutine flux_euler_w_imin_jmax_11pts_SBP4_c

  !===============================================================================
  module subroutine flux_euler_w_imax_jmin_11pts_SBP4_c
  !===============================================================================
    !> Derivatives of Euler fluxes for edge I_MAX/J_MIN close to two walls
    !> - curvilinear coordinate - 11-point stencil - Summation by Parts o4 -
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: vc,rovc
    ! ---------------------------------------------------------------------------

    ! ~> First pass: 11-pt metrics & order reduction for derivatives
    ! ==============================================================

    ! modified fluxes along ksi
    ! -------------------------
    do k=ndz_e,nfz_e
       do j=1,5
          do i=nx-9,nx+5
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
    ! Point #1: stencil [o x x x]
    i=nx
    do j=1,5
       do k=ndz_e,nfz_e
          Krho(i,j,k) = as4m0(1)*Frho(nx  ,j,k) +as4m0(2)*Frho(nx-1,j,k) &
                      + as4m0(3)*Frho(nx-2,j,k) +as4m0(4)*Frho(nx-3,j,k)
          Krhou(i,j,k)= as4m0(1)*Frhou(nx  ,j,k)+as4m0(2)*Frhou(nx-1,j,k) &
                      + as4m0(3)*Frhou(nx-2,j,k)+as4m0(4)*Frhou(nx-3,j,k)
          Krhov(i,j,k)= as4m0(1)*Frhov(nx  ,j,k)+as4m0(2)*Frhov(nx-1,j,k) &
                      + as4m0(3)*Frhov(nx-2,j,k)+as4m0(4)*Frhov(nx-3,j,k)
          Krhow(i,j,k)= as4m0(1)*Frhow(nx  ,j,k)+as4m0(2)*Frhow(nx-1,j,k) &
                      + as4m0(3)*Frhow(nx-2,j,k)+as4m0(4)*Frhow(nx-3,j,k)
          Krhoe(i,j,k)= as4m0(1)*Frhoe(nx  ,j,k)+as4m0(2)*Frhoe(nx-1,j,k) &
                      + as4m0(3)*Frhoe(nx-2,j,k)+as4m0(4)*Frhoe(nx-3,j,k)
       enddo
    enddo

    ! Point #2: stencil [x o x x]
    i=nx-1
    do j=1,5
       do k=ndz_e,nfz_e
          Krho(i,j,k) = as4m1(1)*Frho(nx  ,j,k) +as4m1(2)*Frho(nx-1,j,k) &
                      + as4m1(3)*Frho(nx-2,j,k) +as4m1(4)*Frho(nx-3,j,k)
          Krhou(i,j,k)= as4m1(1)*Frhou(nx  ,j,k)+as4m1(2)*Frhou(nx-1,j,k) &
                      + as4m1(3)*Frhou(nx-2,j,k)+as4m1(4)*Frhou(nx-3,j,k)
          Krhov(i,j,k)= as4m1(1)*Frhov(nx  ,j,k)+as4m1(2)*Frhov(nx-1,j,k) &
                      + as4m1(3)*Frhov(nx-2,j,k)+as4m1(4)*Frhov(nx-3,j,k)
          Krhow(i,j,k)= as4m1(1)*Frhow(nx  ,j,k)+as4m1(2)*Frhow(nx-1,j,k) &
                      + as4m1(3)*Frhow(nx-2,j,k)+as4m1(4)*Frhow(nx-3,j,k)
          Krhoe(i,j,k)= as4m1(1)*Frhoe(nx  ,j,k)+as4m1(2)*Frhoe(nx-1,j,k) &
                      + as4m1(3)*Frhoe(nx-2,j,k)+as4m1(4)*Frhoe(nx-3,j,k)
       enddo
    enddo

    ! Point #3: stencil [x x o x x]
    i=nx-2
    do j=1,5
       do k=ndz_e,nfz_e
          Krho(i,j,k) = as4m2(1)*Frho(nx  ,j,k) +as4m2(2)*Frho(nx-1,j,k) &
                      + as4m2(3)*Frho(nx-2,j,k) +as4m2(4)*Frho(nx-3,j,k) &
                      + as4m2(5)*Frho(nx-4,j,k)
          Krhou(i,j,k)= as4m2(1)*Frhou(nx  ,j,k)+as4m2(2)*Frhou(nx-1,j,k) &
                      + as4m2(3)*Frhou(nx-2,j,k)+as4m2(4)*Frhou(nx-3,j,k) &
                      + as4m2(5)*Frhou(nx-4,j,k)
          Krhov(i,j,k)= as4m2(1)*Frhov(nx  ,j,k)+as4m2(2)*Frhov(nx-1,j,k) &
                      + as4m2(3)*Frhov(nx-2,j,k)+as4m2(4)*Frhov(nx-3,j,k) &
                      + as4m2(5)*Frhov(nx-4,j,k)
          Krhow(i,j,k)= as4m2(1)*Frhow(nx  ,j,k)+as4m2(2)*Frhow(nx-1,j,k) &
                      + as4m2(3)*Frhow(nx-2,j,k)+as4m2(4)*Frhow(nx-3,j,k) &
                      + as4m2(5)*Frhow(nx-4,j,k)
          Krhoe(i,j,k)= as4m2(1)*Frhoe(nx  ,j,k)+as4m2(2)*Frhoe(nx-1,j,k) &
                      + as4m2(3)*Frhoe(nx-2,j,k)+as4m2(4)*Frhoe(nx-3,j,k) &
                      + as4m2(5)*Frhoe(nx-4,j,k)
       enddo
    enddo

    ! Point #4: stencil [x x x o x x]
    i=nx-3
    do j=1,5
       do k=ndz_e,nfz_e
          Krho(i,j,k) = as4m3(1)*Frho(nx  ,j,k) +as4m3(2)*Frho(nx-1,j,k) &
                      + as4m3(3)*Frho(nx-2,j,k) +as4m3(4)*Frho(nx-3,j,k) &
                      + as4m3(5)*Frho(nx-4,j,k) +as4m3(6)*Frho(nx-5,j,k)
          Krhou(i,j,k)= as4m3(1)*Frhou(nx  ,j,k)+as4m3(2)*Frhou(nx-1,j,k) &
                      + as4m3(3)*Frhou(nx-2,j,k)+as4m3(4)*Frhou(nx-3,j,k) &
                      + as4m3(5)*Frhou(nx-4,j,k)+as4m3(6)*Frhou(nx-5,j,k)
          Krhov(i,j,k)= as4m3(1)*Frhov(nx  ,j,k)+as4m3(2)*Frhov(nx-1,j,k) &
                      + as4m3(3)*Frhov(nx-2,j,k)+as4m3(4)*Frhov(nx-3,j,k) &
                      + as4m3(5)*Frhov(nx-4,j,k)+as4m3(6)*Frhov(nx-5,j,k)
          Krhow(i,j,k)= as4m3(1)*Frhow(nx  ,j,k)+as4m3(2)*Frhow(nx-1,j,k) &
                      + as4m3(3)*Frhow(nx-2,j,k)+as4m3(4)*Frhow(nx-3,j,k) &
                      + as4m3(5)*Frhow(nx-4,j,k)+as4m3(6)*Frhow(nx-5,j,k)
          Krhoe(i,j,k)= as4m3(1)*Frhoe(nx  ,j,k)+as4m3(2)*Frhoe(nx-1,j,k) &
                      + as4m3(3)*Frhoe(nx-2,j,k)+as4m3(4)*Frhoe(nx-3,j,k) &
                      + as4m3(5)*Frhoe(nx-4,j,k)+as4m3(6)*Frhoe(nx-5,j,k)
       enddo
    enddo

    ! Eighth-order (9-pt stencil)
    i=nx-4
    do j=1,5
       do k=ndz_e,nfz_e
          Krho(i,j,k)  = a9(1) * ( Frho(i+1,j,k)-Frho(i-1,j,k) ) &
                       + a9(2) * ( Frho(i+2,j,k)-Frho(i-2,j,k) ) &
                       + a9(3) * ( Frho(i+3,j,k)-Frho(i-3,j,k) ) &
                       + a9(4) * ( Frho(i+4,j,k)-Frho(i-4,j,k) )
          Krhou(i,j,k) = a9(1) * ( Frhou(i+1,j,k)-Frhou(i-1,j,k) ) &
                       + a9(2) * ( Frhou(i+2,j,k)-Frhou(i-2,j,k) ) &
                       + a9(3) * ( Frhou(i+3,j,k)-Frhou(i-3,j,k) ) &
                       + a9(4) * ( Frhou(i+4,j,k)-Frhou(i-4,j,k) )
          Krhov(i,j,k) = a9(1) * ( Frhov(i+1,j,k)-Frhov(i-1,j,k) ) &
                       + a9(2) * ( Frhov(i+2,j,k)-Frhov(i-2,j,k) ) &
                       + a9(3) * ( Frhov(i+3,j,k)-Frhov(i-3,j,k) ) &
                       + a9(4) * ( Frhov(i+4,j,k)-Frhov(i-4,j,k) )
          Krhow(i,j,k) = a9(1) * ( Frhow(i+1,j,k)-Frhow(i-1,j,k) ) &
                       + a9(2) * ( Frhow(i+2,j,k)-Frhow(i-2,j,k) ) &
                       + a9(3) * ( Frhow(i+3,j,k)-Frhow(i-3,j,k) ) &
                       + a9(4) * ( Frhow(i+4,j,k)-Frhow(i-4,j,k) )
          Krhoe(i,j,k) = a9(1) * ( Frhoe(i+1,j,k)-Frhoe(i-1,j,k) ) &
                       + a9(2) * ( Frhoe(i+2,j,k)-Frhoe(i-2,j,k) ) &
                       + a9(3) * ( Frhoe(i+3,j,k)-Frhoe(i-3,j,k) ) &
                       + a9(4) * ( Frhoe(i+4,j,k)-Frhoe(i-4,j,k) )
       enddo
    enddo

    ! modified fluxes along eta
    ! -------------------------
    do k=ndz_e,nfz_e
       do i=nx-4,nx
          do j=-4,10
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
    ! Point #1: stencil [o x x x]
    !j=1
    do i=nx-4,nx
       do k=ndz_e,nfz_e  
          Krho(i,1,k) = as4p0(1)*Frho(i,1,k) +as4p0(2)*Frho(i,2,k) &
                      + as4p0(3)*Frho(i,3,k) +as4p0(4)*Frho(i,4,k) + Krho(i,1,k)
          Krhou(i,1,k)= as4p0(1)*Frhou(i,1,k)+as4p0(2)*Frhou(i,2,k) &
                      + as4p0(3)*Frhou(i,3,k)+as4p0(4)*Frhou(i,4,k) + Krhou(i,1,k)
          Krhov(i,1,k)= as4p0(1)*Frhov(i,1,k)+as4p0(2)*Frhov(i,2,k) &
                      + as4p0(3)*Frhov(i,3,k)+as4p0(4)*Frhov(i,4,k) + Krhov(i,1,k)
          Krhow(i,1,k)= as4p0(1)*Frhow(i,1,k)+as4p0(2)*Frhow(i,2,k) &
                      + as4p0(3)*Frhow(i,3,k)+as4p0(4)*Frhow(i,4,k) + Krhow(i,1,k)
          Krhoe(i,1,k)= as4p0(1)*Frhoe(i,1,k)+as4p0(2)*Frhoe(i,2,k) &
                      + as4p0(3)*Frhoe(i,3,k)+as4p0(4)*Frhoe(i,4,k) + Krhoe(i,1,k)
       enddo
    enddo

    ! Point #2: stencil [x o x x]
    !j=2
    do i=nx-4,nx
       do k=ndz_e,nfz_e  
          Krho(i,2,k) = as4p1(1)*Frho(i,1,k) +as4p1(2)*Frho(i,2,k) &
                      + as4p1(3)*Frho(i,3,k) +as4p1(4)*Frho(i,4,k) + Krho(i,2,k)
          Krhou(i,2,k)= as4p1(1)*Frhou(i,1,k)+as4p1(2)*Frhou(i,2,k) &
                      + as4p1(3)*Frhou(i,3,k)+as4p1(4)*Frhou(i,4,k) + Krhou(i,2,k)
          Krhov(i,2,k)= as4p1(1)*Frhov(i,1,k)+as4p1(2)*Frhov(i,2,k) &
                      + as4p1(3)*Frhov(i,3,k)+as4p1(4)*Frhov(i,4,k) + Krhov(i,2,k)
          Krhow(i,2,k)= as4p1(1)*Frhow(i,1,k)+as4p1(2)*Frhow(i,2,k) &
                      + as4p1(3)*Frhow(i,3,k)+as4p1(4)*Frhow(i,4,k) + Krhow(i,2,k)
          Krhoe(i,2,k)= as4p1(1)*Frhoe(i,1,k)+as4p1(2)*Frhoe(i,2,k) &
                      + as4p1(3)*Frhoe(i,3,k)+as4p1(4)*Frhoe(i,4,k) + Krhoe(i,2,k)
       enddo
    enddo

    ! Point #3: stencil [x x o x x]
    !j=3
    do i=nx-4,nx
       do k=ndz_e,nfz_e
          Krho(i,3,k) = as4p2(1)*Frho(i,1,k) +as4p2(2)*Frho(i,2,k) &
                      + as4p2(3)*Frho(i,3,k) +as4p2(4)*Frho(i,4,k) &
                      + as4p2(5)*Frho(i,5,k) + Krho(i,3,k)
          Krhou(i,3,k)= as4p2(1)*Frhou(i,1,k)+as4p2(2)*Frhou(i,2,k) &
                      + as4p2(3)*Frhou(i,3,k)+as4p2(4)*Frhou(i,4,k) &
                      + as4p2(5)*Frhou(i,5,k) + Krhou(i,3,k)
          Krhov(i,3,k)= as4p2(1)*Frhov(i,1,k)+as4p2(2)*Frhov(i,2,k) &
                      + as4p2(3)*Frhov(i,3,k)+as4p2(4)*Frhov(i,4,k) &
                      + as4p2(5)*Frhov(i,5,k) + Krhov(i,3,k)
          Krhow(i,3,k)= as4p2(1)*Frhow(i,1,k)+as4p2(2)*Frhow(i,2,k) &
                      + as4p2(3)*Frhow(i,3,k)+as4p2(4)*Frhow(i,4,k) &
                      + as4p2(5)*Frhow(i,5,k) + Krhow(i,3,k)
          Krhoe(i,3,k)= as4p2(1)*Frhoe(i,1,k)+as4p2(2)*Frhoe(i,2,k) &
                      + as4p2(3)*Frhoe(i,3,k)+as4p2(4)*Frhoe(i,4,k) &
                      + as4p2(5)*Frhoe(i,5,k) + Krhoe(i,3,k)
       enddo
    enddo

    ! Point #4: stencil [x x x o x x]
    !j=4
    do i=nx-4,nx
       do k=ndz_e,nfz_e
          Krho(i,4,k) = as4p3(1)*Frho(i,1,k) +as4p3(2)*Frho(i,2,k) &
                      + as4p3(3)*Frho(i,3,k) +as4p3(4)*Frho(i,4,k) &
                      + as4p3(5)*Frho(i,5,k) +as4p3(6)*Frho(i,6,k) + Krho(i,4,k)
          Krhou(i,4,k)= as4p3(1)*Frhou(i,1,k)+as4p3(2)*Frhou(i,2,k) &
                      + as4p3(3)*Frhou(i,3,k)+as4p3(4)*Frhou(i,4,k) &
                      + as4p3(5)*Frhou(i,5,k)+as4p3(6)*Frhou(i,6,k) + Krhou(i,4,k)
          Krhov(i,4,k)= as4p3(1)*Frhov(i,1,k)+as4p3(2)*Frhov(i,2,k) &
                      + as4p3(3)*Frhov(i,3,k)+as4p3(4)*Frhov(i,4,k) &
                      + as4p3(5)*Frhov(i,5,k)+as4p3(6)*Frhov(i,6,k) + Krhov(i,4,k)
          Krhow(i,4,k)= as4p3(1)*Frhow(i,1,k)+as4p3(2)*Frhow(i,2,k) &
                      + as4p3(3)*Frhow(i,3,k)+as4p3(4)*Frhow(i,4,k) &
                      + as4p3(5)*Frhow(i,5,k)+as4p3(6)*Frhow(i,6,k) + Krhow(i,4,k)
          Krhoe(i,4,k)= as4p3(1)*Frhoe(i,1,k)+as4p3(2)*Frhoe(i,2,k) &
                      + as4p3(3)*Frhoe(i,3,k)+as4p3(4)*Frhoe(i,4,k) &
                      + as4p3(5)*Frhoe(i,5,k)+as4p3(6)*Frhoe(i,6,k) + Krhoe(i,4,k)
       enddo
    enddo

    ! Eighth-order (9-pt stencil)
    j=5
    do i=nx-4,nx
       do k=ndz_e,nfz_e
          Krho(i,j,k)  = a9(1) * ( Frho(i,j+1,k)-Frho(i,j-1,k) ) &
                       + a9(2) * ( Frho(i,j+2,k)-Frho(i,j-2,k) ) &
                       + a9(3) * ( Frho(i,j+3,k)-Frho(i,j-3,k) ) &
                       + a9(4) * ( Frho(i,j+4,k)-Frho(i,j-4,k) ) + Krho(i,j,k)      
          Krhou(i,j,k) = a9(1) * ( Frhou(i,j+1,k)-Frhou(i,j-1,k) ) &
                       + a9(2) * ( Frhou(i,j+2,k)-Frhou(i,j-2,k) ) &
                       + a9(3) * ( Frhou(i,j+3,k)-Frhou(i,j-3,k) ) &
                       + a9(4) * ( Frhou(i,j+4,k)-Frhou(i,j-4,k) ) + Krhou(i,j,k)
          Krhov(i,j,k) = a9(1) * ( Frhov(i,j+1,k)-Frhov(i,j-1,k) ) &
                       + a9(2) * ( Frhov(i,j+2,k)-Frhov(i,j-2,k) ) &
                       + a9(3) * ( Frhov(i,j+3,k)-Frhov(i,j-3,k) ) &
                       + a9(4) * ( Frhov(i,j+4,k)-Frhov(i,j-4,k) ) + Krhov(i,j,k)
          Krhow(i,j,k) = a9(1) * ( Frhow(i,j+1,k)-Frhow(i,j-1,k) ) &
                       + a9(2) * ( Frhow(i,j+2,k)-Frhow(i,j-2,k) ) &
                       + a9(3) * ( Frhow(i,j+3,k)-Frhow(i,j-3,k) ) &
                       + a9(4) * ( Frhow(i,j+4,k)-Frhow(i,j-4,k) ) + Krhow(i,j,k)
          Krhoe(i,j,k) = a9(1) * ( Frhoe(i,j+1,k)-Frhoe(i,j-1,k) ) &
                       + a9(2) * ( Frhoe(i,j+2,k)-Frhoe(i,j-2,k) ) &
                       + a9(3) * ( Frhoe(i,j+3,k)-Frhoe(i,j-3,k) ) &
                       + a9(4) * ( Frhoe(i,j+4,k)-Frhoe(i,j-4,k) ) + Krhoe(i,j,k)
       enddo
    enddo

    ! ~> Second pass: order reduction for metrics & 11-pt derivatives
    ! ===============================================================

    ! modified fluxes along ksi
    ! -------------------------
    do k=ndz_e,nfz_e
       do j=1,5
          do i=nx-9,nx+5
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
       do j=1,5
          do i=nx-4,nx
             Krho(i,j,k)  = a11(1) * ( Frho(i+1,j,k)-Frho(i-1,j,k) )   &
                          + a11(2) * ( Frho(i+2,j,k)-Frho(i-2,j,k) )   &
                          + a11(3) * ( Frho(i+3,j,k)-Frho(i-3,j,k) )   &
                          + a11(4) * ( Frho(i+4,j,k)-Frho(i-4,j,k) )   &
                          + a11(5) * ( Frho(i+5,j,k)-Frho(i-5,j,k) ) + Krho(i,j,k)

             Krhou(i,j,k) = a11(1) * ( Frhou(i+1,j,k)-Frhou(i-1,j,k) ) &
                          + a11(2) * ( Frhou(i+2,j,k)-Frhou(i-2,j,k) ) &
                          + a11(3) * ( Frhou(i+3,j,k)-Frhou(i-3,j,k) ) &
                          + a11(4) * ( Frhou(i+4,j,k)-Frhou(i-4,j,k) ) &
                          + a11(5) * ( Frhou(i+5,j,k)-Frhou(i-5,j,k) ) + Krhou(i,j,k)

             Krhov(i,j,k) = a11(1) * ( Frhov(i+1,j,k)-Frhov(i-1,j,k) ) &
                          + a11(2) * ( Frhov(i+2,j,k)-Frhov(i-2,j,k) ) &
                          + a11(3) * ( Frhov(i+3,j,k)-Frhov(i-3,j,k) ) &
                          + a11(4) * ( Frhov(i+4,j,k)-Frhov(i-4,j,k) ) &
                          + a11(5) * ( Frhov(i+5,j,k)-Frhov(i-5,j,k) ) + Krhov(i,j,k)

             Krhow(i,j,k) = a11(1) * ( Frhow(i+1,j,k)-Frhow(i-1,j,k) ) &
                          + a11(2) * ( Frhow(i+2,j,k)-Frhow(i-2,j,k) ) &
                          + a11(3) * ( Frhow(i+3,j,k)-Frhow(i-3,j,k) ) &
                          + a11(4) * ( Frhow(i+4,j,k)-Frhow(i-4,j,k) ) &
                          + a11(5) * ( Frhow(i+5,j,k)-Frhow(i-5,j,k) ) + Krhow(i,j,k)

             Krhoe(i,j,k) = a11(1) * ( Frhoe(i+1,j,k)-Frhoe(i-1,j,k) ) &
                          + a11(2) * ( Frhoe(i+2,j,k)-Frhoe(i-2,j,k) ) &
                          + a11(3) * ( Frhoe(i+3,j,k)-Frhoe(i-3,j,k) ) &
                          + a11(4) * ( Frhoe(i+4,j,k)-Frhoe(i-4,j,k) ) &
                          + a11(5) * ( Frhoe(i+5,j,k)-Frhoe(i-5,j,k) ) + Krhoe(i,j,k)
          enddo
       enddo
    enddo

    ! modified fluxes along eta
    ! -------------------------
    do k=ndz_e,nfz_e
       do i=nx-4,nx
          do j=-4,10
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
       do j=1,5
          do i=nx-4,nx
             Krho(i,j,k)  = a11(1) * ( Frho(i,j+1,k)-Frho(i,j-1,k) ) &
                          + a11(2) * ( Frho(i,j+2,k)-Frho(i,j-2,k) ) &
                          + a11(3) * ( Frho(i,j+3,k)-Frho(i,j-3,k) ) &
                          + a11(4) * ( Frho(i,j+4,k)-Frho(i,j-4,k) ) &
                          + a11(5) * ( Frho(i,j+5,k)-Frho(i,j-5,k) ) &
                          + Krho(i,j,k)

             Krhou(i,j,k) = a11(1) * ( Frhou(i,j+1,k)-Frhou(i,j-1,k) ) &
                          + a11(2) * ( Frhou(i,j+2,k)-Frhou(i,j-2,k) ) &
                          + a11(3) * ( Frhou(i,j+3,k)-Frhou(i,j-3,k) ) &
                          + a11(4) * ( Frhou(i,j+4,k)-Frhou(i,j-4,k) ) &
                          + a11(5) * ( Frhou(i,j+5,k)-Frhou(i,j-5,k) ) &
                          + Krhou(i,j,k)

             Krhov(i,j,k) = a11(1) * ( Frhov(i,j+1,k)-Frhov(i,j-1,k) ) &
                          + a11(2) * ( Frhov(i,j+2,k)-Frhov(i,j-2,k) ) &
                          + a11(3) * ( Frhov(i,j+3,k)-Frhov(i,j-3,k) ) &
                          + a11(4) * ( Frhov(i,j+4,k)-Frhov(i,j-4,k) ) &
                          + a11(5) * ( Frhov(i,j+5,k)-Frhov(i,j-5,k) ) &
                          + Krhov(i,j,k)

             Krhow(i,j,k) = a11(1) * ( Frhow(i,j+1,k)-Frhow(i,j-1,k) ) &
                          + a11(2) * ( Frhow(i,j+2,k)-Frhow(i,j-2,k) ) &
                          + a11(3) * ( Frhow(i,j+3,k)-Frhow(i,j-3,k) ) &
                          + a11(4) * ( Frhow(i,j+4,k)-Frhow(i,j-4,k) ) &
                          + a11(5) * ( Frhow(i,j+5,k)-Frhow(i,j-5,k) ) &
                          + Krhow(i,j,k)

             Krhoe(i,j,k) = a11(1) * ( Frhoe(i,j+1,k)-Frhoe(i,j-1,k) ) &
                          + a11(2) * ( Frhoe(i,j+2,k)-Frhoe(i,j-2,k) ) &
                          + a11(3) * ( Frhoe(i,j+3,k)-Frhoe(i,j-3,k) ) &
                          + a11(4) * ( Frhoe(i,j+4,k)-Frhoe(i,j-4,k) ) &
                          + a11(5) * ( Frhoe(i,j+5,k)-Frhoe(i,j-5,k) ) &
                          + Krhoe(i,j,k)
          enddo
       enddo
    enddo

    ! Average of the two derivative evaluations
    ! =========================================
    do k=ndz_e,nfz_e
       do j=1,5
          do i=nx-4,nx
             Krho(i,j,k)  = Krho(i,j,k)*0.5_wp
             Krhou(i,j,k) = Krhou(i,j,k)*0.5_wp
             Krhov(i,j,k) = Krhov(i,j,k)*0.5_wp
             Krhow(i,j,k) = Krhow(i,j,k)*0.5_wp
             Krhoe(i,j,k) = Krhoe(i,j,k)*0.5_wp
          enddo
       enddo
    enddo

  end subroutine flux_euler_w_imax_jmin_11pts_SBP4_c

  !===============================================================================
  module subroutine flux_euler_w_imax_jmax_11pts_SBP4_c
  !===============================================================================
    !> Derivatives of Euler fluxes for edge I_MAX/J_MAX close to two walls
    !> - curvilinear coordinate - 11-point stencil - Summation by Parts o4 -
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: vc,rovc
    ! ---------------------------------------------------------------------------

    ! ~> First pass: 11-pt metrics & order reduction for derivatives
    ! ==============================================================

    ! modified fluxes along ksi
    ! -------------------------
    do k=ndz_e,nfz_e
       do j=ny-4,ny
          do i=nx-9,nx+5
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
    ! Point #1: stencil [o x x x]
    i=nx
    do j=ny-4,ny
       do k=ndz_e,nfz_e
          Krho(i,j,k) = as4m0(1)*Frho(nx  ,j,k) +as4m0(2)*Frho(nx-1,j,k) &
                      + as4m0(3)*Frho(nx-2,j,k) +as4m0(4)*Frho(nx-3,j,k)
          Krhou(i,j,k)= as4m0(1)*Frhou(nx  ,j,k)+as4m0(2)*Frhou(nx-1,j,k) &
                      + as4m0(3)*Frhou(nx-2,j,k)+as4m0(4)*Frhou(nx-3,j,k)
          Krhov(i,j,k)= as4m0(1)*Frhov(nx  ,j,k)+as4m0(2)*Frhov(nx-1,j,k) &
                      + as4m0(3)*Frhov(nx-2,j,k)+as4m0(4)*Frhov(nx-3,j,k)
          Krhow(i,j,k)= as4m0(1)*Frhow(nx  ,j,k)+as4m0(2)*Frhow(nx-1,j,k) &
                      + as4m0(3)*Frhow(nx-2,j,k)+as4m0(4)*Frhow(nx-3,j,k)
          Krhoe(i,j,k)= as4m0(1)*Frhoe(nx  ,j,k)+as4m0(2)*Frhoe(nx-1,j,k) &
                      + as4m0(3)*Frhoe(nx-2,j,k)+as4m0(4)*Frhoe(nx-3,j,k)
       enddo
    enddo

    ! Point #2: stencil [x o x x]
    i=nx-1
    do j=ny-4,ny
       do k=ndz_e,nfz_e
          Krho(i,j,k) = as4m1(1)*Frho(nx  ,j,k) +as4m1(2)*Frho(nx-1,j,k) &
                      + as4m1(3)*Frho(nx-2,j,k) +as4m1(4)*Frho(nx-3,j,k)
          Krhou(i,j,k)= as4m1(1)*Frhou(nx  ,j,k)+as4m1(2)*Frhou(nx-1,j,k) &
                      + as4m1(3)*Frhou(nx-2,j,k)+as4m1(4)*Frhou(nx-3,j,k)
          Krhov(i,j,k)= as4m1(1)*Frhov(nx  ,j,k)+as4m1(2)*Frhov(nx-1,j,k) &
                      + as4m1(3)*Frhov(nx-2,j,k)+as4m1(4)*Frhov(nx-3,j,k)
          Krhow(i,j,k)= as4m1(1)*Frhow(nx  ,j,k)+as4m1(2)*Frhow(nx-1,j,k) &
                      + as4m1(3)*Frhow(nx-2,j,k)+as4m1(4)*Frhow(nx-3,j,k)
          Krhoe(i,j,k)= as4m1(1)*Frhoe(nx  ,j,k)+as4m1(2)*Frhoe(nx-1,j,k) &
                      + as4m1(3)*Frhoe(nx-2,j,k)+as4m1(4)*Frhoe(nx-3,j,k)
       enddo
    enddo

    ! Point #3: stencil [x x o x x]
    i=nx-2
    do j=ny-4,ny
       do k=ndz_e,nfz_e
          Krho(i,j,k) = as4m2(1)*Frho(nx  ,j,k) +as4m2(2)*Frho(nx-1,j,k) &
                      + as4m2(3)*Frho(nx-2,j,k) +as4m2(4)*Frho(nx-3,j,k) &
                      + as4m2(5)*Frho(nx-4,j,k)
          Krhou(i,j,k)= as4m2(1)*Frhou(nx  ,j,k)+as4m2(2)*Frhou(nx-1,j,k) &
                      + as4m2(3)*Frhou(nx-2,j,k)+as4m2(4)*Frhou(nx-3,j,k) &
                      + as4m2(5)*Frhou(nx-4,j,k)
          Krhov(i,j,k)= as4m2(1)*Frhov(nx  ,j,k)+as4m2(2)*Frhov(nx-1,j,k) &
                      + as4m2(3)*Frhov(nx-2,j,k)+as4m2(4)*Frhov(nx-3,j,k) &
                      + as4m2(5)*Frhov(nx-4,j,k)
          Krhow(i,j,k)= as4m2(1)*Frhow(nx  ,j,k)+as4m2(2)*Frhow(nx-1,j,k) &
                      + as4m2(3)*Frhow(nx-2,j,k)+as4m2(4)*Frhow(nx-3,j,k) &
                      + as4m2(5)*Frhow(nx-4,j,k)
          Krhoe(i,j,k)= as4m2(1)*Frhoe(nx  ,j,k)+as4m2(2)*Frhoe(nx-1,j,k) &
                      + as4m2(3)*Frhoe(nx-2,j,k)+as4m2(4)*Frhoe(nx-3,j,k) &
                      + as4m2(5)*Frhoe(nx-4,j,k)
       enddo
    enddo

    ! Point #4: stencil [x x x o x x]
    i=nx-3
    do j=ny-4,ny
       do k=ndz_e,nfz_e
          Krho(i,j,k) = as4m3(1)*Frho(nx  ,j,k) +as4m3(2)*Frho(nx-1,j,k) &
                      + as4m3(3)*Frho(nx-2,j,k) +as4m3(4)*Frho(nx-3,j,k) &
                      + as4m3(5)*Frho(nx-4,j,k) +as4m3(6)*Frho(nx-5,j,k)
          Krhou(i,j,k)= as4m3(1)*Frhou(nx  ,j,k)+as4m3(2)*Frhou(nx-1,j,k) &
                      + as4m3(3)*Frhou(nx-2,j,k)+as4m3(4)*Frhou(nx-3,j,k) &
                      + as4m3(5)*Frhou(nx-4,j,k)+as4m3(6)*Frhou(nx-5,j,k)
          Krhov(i,j,k)= as4m3(1)*Frhov(nx  ,j,k)+as4m3(2)*Frhov(nx-1,j,k) &
                      + as4m3(3)*Frhov(nx-2,j,k)+as4m3(4)*Frhov(nx-3,j,k) &
                      + as4m3(5)*Frhov(nx-4,j,k)+as4m3(6)*Frhov(nx-5,j,k)
          Krhow(i,j,k)= as4m3(1)*Frhow(nx  ,j,k)+as4m3(2)*Frhow(nx-1,j,k) &
                      + as4m3(3)*Frhow(nx-2,j,k)+as4m3(4)*Frhow(nx-3,j,k) &
                      + as4m3(5)*Frhow(nx-4,j,k)+as4m3(6)*Frhow(nx-5,j,k)
          Krhoe(i,j,k)= as4m3(1)*Frhoe(nx  ,j,k)+as4m3(2)*Frhoe(nx-1,j,k) &
                      + as4m3(3)*Frhoe(nx-2,j,k)+as4m3(4)*Frhoe(nx-3,j,k) &
                      + as4m3(5)*Frhoe(nx-4,j,k)+as4m3(6)*Frhoe(nx-5,j,k)
       enddo
    enddo

    ! Eighth-order (9-pt stencil)
    i=nx-4
    do j=ny-4,ny
       do k=ndz_e,nfz_e
          Krho(i,j,k)  = a9(1) * ( Frho(i+1,j,k)-Frho(i-1,j,k) ) &
                       + a9(2) * ( Frho(i+2,j,k)-Frho(i-2,j,k) ) &
                       + a9(3) * ( Frho(i+3,j,k)-Frho(i-3,j,k) ) &
                       + a9(4) * ( Frho(i+4,j,k)-Frho(i-4,j,k) )
          Krhou(i,j,k) = a9(1) * ( Frhou(i+1,j,k)-Frhou(i-1,j,k) ) &
                       + a9(2) * ( Frhou(i+2,j,k)-Frhou(i-2,j,k) ) &
                       + a9(3) * ( Frhou(i+3,j,k)-Frhou(i-3,j,k) ) &
                       + a9(4) * ( Frhou(i+4,j,k)-Frhou(i-4,j,k) )
          Krhov(i,j,k) = a9(1) * ( Frhov(i+1,j,k)-Frhov(i-1,j,k) ) &
                       + a9(2) * ( Frhov(i+2,j,k)-Frhov(i-2,j,k) ) &
                       + a9(3) * ( Frhov(i+3,j,k)-Frhov(i-3,j,k) ) &
                       + a9(4) * ( Frhov(i+4,j,k)-Frhov(i-4,j,k) )
          Krhow(i,j,k) = a9(1) * ( Frhow(i+1,j,k)-Frhow(i-1,j,k) ) &
                       + a9(2) * ( Frhow(i+2,j,k)-Frhow(i-2,j,k) ) &
                       + a9(3) * ( Frhow(i+3,j,k)-Frhow(i-3,j,k) ) &
                       + a9(4) * ( Frhow(i+4,j,k)-Frhow(i-4,j,k) )
          Krhoe(i,j,k) = a9(1) * ( Frhoe(i+1,j,k)-Frhoe(i-1,j,k) ) &
                       + a9(2) * ( Frhoe(i+2,j,k)-Frhoe(i-2,j,k) ) &
                       + a9(3) * ( Frhoe(i+3,j,k)-Frhoe(i-3,j,k) ) &
                       + a9(4) * ( Frhoe(i+4,j,k)-Frhoe(i-4,j,k) )
       enddo
    enddo

    ! modified fluxes along eta
    ! -------------------------
    do k=ndz_e,nfz_e
       do i=nx-4,nx
          do j=ny-9,ny+5
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
    ! Point #1: stencil [o x x x]
    j=ny
    do i=nx-4,nx
       do k=ndz_e,nfz_e  
          Krho(i,j,k) = as4m0(1)*Frho(i,ny  ,k) +as4m0(2)*Frho(i,ny-1,k) &
                      + as4m0(3)*Frho(i,ny-2,k) +as4m0(4)*Frho(i,ny-3,k) + Krho(i,j,k)
          Krhou(i,j,k)= as4m0(1)*Frhou(i,ny  ,k)+as4m0(2)*Frhou(i,ny-1,k) &
                      + as4m0(3)*Frhou(i,ny-2,k)+as4m0(4)*Frhou(i,ny-3,k) + Krhou(i,j,k)
          Krhov(i,j,k)= as4m0(1)*Frhov(i,ny  ,k)+as4m0(2)*Frhov(i,ny-1,k) &
                      + as4m0(3)*Frhov(i,ny-2,k)+as4m0(4)*Frhov(i,ny-3,k) + Krhov(i,j,k)
          Krhow(i,j,k)= as4m0(1)*Frhow(i,ny  ,k)+as4m0(2)*Frhow(i,ny-1,k) &
                      + as4m0(3)*Frhow(i,ny-2,k)+as4m0(4)*Frhow(i,ny-3,k) + Krhow(i,j,k)
          Krhoe(i,j,k)= as4m0(1)*Frhoe(i,ny  ,k)+as4m0(2)*Frhoe(i,ny-1,k) &
                      + as4m0(3)*Frhoe(i,ny-2,k)+as4m0(4)*Frhoe(i,ny-3,k) + Krhoe(i,j,k)
       enddo
    enddo

    ! Point #2: stencil [x o x x]
    j=ny-1
    do i=nx-4,nx
       do k=ndz_e,nfz_e  
          Krho(i,j,k) = as4m1(1)*Frho(i,ny  ,k) +as4m1(2)*Frho(i,ny-1,k) &
                      + as4m1(3)*Frho(i,ny-2,k) +as4m1(4)*Frho(i,ny-3,k) + Krho(i,j,k)
          Krhou(i,j,k)= as4m1(1)*Frhou(i,ny  ,k)+as4m1(2)*Frhou(i,ny-1,k) &
                      + as4m1(3)*Frhou(i,ny-2,k)+as4m1(4)*Frhou(i,ny-3,k) + Krhou(i,j,k)
          Krhov(i,j,k)= as4m1(1)*Frhov(i,ny  ,k)+as4m1(2)*Frhov(i,ny-1,k) &
                      + as4m1(3)*Frhov(i,ny-2,k)+as4m1(4)*Frhov(i,ny-3,k) + Krhov(i,j,k)
          Krhow(i,j,k)= as4m1(1)*Frhow(i,ny  ,k)+as4m1(2)*Frhow(i,ny-1,k) &
                      + as4m1(3)*Frhow(i,ny-2,k)+as4m1(4)*Frhow(i,ny-3,k) + Krhow(i,j,k)
          Krhoe(i,j,k)= as4m1(1)*Frhoe(i,ny  ,k)+as4m1(2)*Frhoe(i,ny-1,k) &
                      + as4m1(3)*Frhoe(i,ny-2,k)+as4m1(4)*Frhoe(i,ny-3,k) + Krhoe(i,j,k)
       enddo
    enddo

    ! Point #3: stencil [x x o x x]
    j=ny-2
    do i=nx-4,nx
       do k=ndz_e,nfz_e
          Krho(i,j,k) = as4m2(1)*Frho(i,ny  ,k) +as4m2(2)*Frho(i,ny-1,k) &
                      + as4m2(3)*Frho(i,ny-2,k) +as4m2(4)*Frho(i,ny-3,k) &
                      + as4m2(5)*Frho(i,ny-4,k) + Krho(i,j,k)
          Krhou(i,j,k)= as4m2(1)*Frhou(i,ny  ,k)+as4m2(2)*Frhou(i,ny-1,k) &
                      + as4m2(3)*Frhou(i,ny-2,k)+as4m2(4)*Frhou(i,ny-3,k) &
                      + as4m2(5)*Frhou(i,ny-4,k)+ Krhou(i,j,k)
          Krhov(i,j,k)= as4m2(1)*Frhov(i,ny  ,k)+as4m2(2)*Frhov(i,ny-1,k) &
                      + as4m2(3)*Frhov(i,ny-2,k)+as4m2(4)*Frhov(i,ny-3,k) &
                      + as4m2(5)*Frhov(i,ny-4,k)+ Krhov(i,j,k)
          Krhow(i,j,k)= as4m2(1)*Frhow(i,ny  ,k)+as4m2(2)*Frhow(i,ny-1,k) &
                      + as4m2(3)*Frhow(i,ny-2,k)+as4m2(4)*Frhow(i,ny-3,k) &
                      + as4m2(5)*Frhow(i,ny-4,k)+ Krhow(i,j,k)
          Krhoe(i,j,k)= as4m2(1)*Frhoe(i,ny  ,k)+as4m2(2)*Frhoe(i,ny-1,k) &
                      + as4m2(3)*Frhoe(i,ny-2,k)+as4m2(4)*Frhoe(i,ny-3,k) &
                      + as4m2(5)*Frhoe(i,ny-4,k)+ Krhoe(i,j,k)
       enddo
    enddo

    ! Point #4: stencil [x x x o x x]
    j=ny-3
    do i=nx-4,nx
       do k=ndz_e,nfz_e
          Krho(i,j,k) = as4m3(1)*Frho(i,ny  ,k) +as4m3(2)*Frho(i,ny-1,k) &
                      + as4m3(3)*Frho(i,ny-2,k) +as4m3(4)*Frho(i,ny-3,k) &
                      + as4m3(5)*Frho(i,ny-4,k) +as4m3(6)*Frho(i,ny-5,k) + Krho(i,j,k)
          Krhou(i,j,k)= as4m3(1)*Frhou(i,ny  ,k)+as4m3(2)*Frhou(i,ny-1,k) &
                      + as4m3(3)*Frhou(i,ny-2,k)+as4m3(4)*Frhou(i,ny-3,k) &
                      + as4m3(5)*Frhou(i,ny-4,k)+as4m3(6)*Frhou(i,ny-5,k) + Krhou(i,j,k)
          Krhov(i,j,k)= as4m3(1)*Frhov(i,ny  ,k)+as4m3(2)*Frhov(i,ny-1,k) &
                      + as4m3(3)*Frhov(i,ny-2,k)+as4m3(4)*Frhov(i,ny-3,k) &
                      + as4m3(5)*Frhov(i,ny-4,k)+as4m3(6)*Frhov(i,ny-5,k) + Krhov(i,j,k)
          Krhow(i,j,k)= as4m3(1)*Frhow(i,ny  ,k)+as4m3(2)*Frhow(i,ny-1,k) &
                      + as4m3(3)*Frhow(i,ny-2,k)+as4m3(4)*Frhow(i,ny-3,k) &
                      + as4m3(5)*Frhow(i,ny-4,k)+as4m3(6)*Frhow(i,ny-5,k) + Krhow(i,j,k)
          Krhoe(i,j,k)= as4m3(1)*Frhoe(i,ny  ,k)+as4m3(2)*Frhoe(i,ny-1,k) &
                      + as4m3(3)*Frhoe(i,ny-2,k)+as4m3(4)*Frhoe(i,ny-3,k) &
                      + as4m3(5)*Frhoe(i,ny-4,k)+as4m3(6)*Frhoe(i,ny-5,k) + Krhoe(i,j,k)
       enddo
    enddo

    ! Eighth-order (9-pt stencil)
    j=ny-4
    do i=nx-4,nx
       do k=ndz_e,nfz_e
          Krho(i,j,k)  = a9(1) * ( Frho(i,j+1,k)-Frho(i,j-1,k) ) &
                       + a9(2) * ( Frho(i,j+2,k)-Frho(i,j-2,k) ) &
                       + a9(3) * ( Frho(i,j+3,k)-Frho(i,j-3,k) ) &
                       + a9(4) * ( Frho(i,j+4,k)-Frho(i,j-4,k) ) + Krho(i,j,k)      
          Krhou(i,j,k) = a9(1) * ( Frhou(i,j+1,k)-Frhou(i,j-1,k) ) &
                       + a9(2) * ( Frhou(i,j+2,k)-Frhou(i,j-2,k) ) &
                       + a9(3) * ( Frhou(i,j+3,k)-Frhou(i,j-3,k) ) &
                       + a9(4) * ( Frhou(i,j+4,k)-Frhou(i,j-4,k) ) + Krhou(i,j,k)
          Krhov(i,j,k) = a9(1) * ( Frhov(i,j+1,k)-Frhov(i,j-1,k) ) &
                       + a9(2) * ( Frhov(i,j+2,k)-Frhov(i,j-2,k) ) &
                       + a9(3) * ( Frhov(i,j+3,k)-Frhov(i,j-3,k) ) &
                       + a9(4) * ( Frhov(i,j+4,k)-Frhov(i,j-4,k) ) + Krhov(i,j,k)
          Krhow(i,j,k) = a9(1) * ( Frhow(i,j+1,k)-Frhow(i,j-1,k) ) &
                       + a9(2) * ( Frhow(i,j+2,k)-Frhow(i,j-2,k) ) &
                       + a9(3) * ( Frhow(i,j+3,k)-Frhow(i,j-3,k) ) &
                       + a9(4) * ( Frhow(i,j+4,k)-Frhow(i,j-4,k) ) + Krhow(i,j,k)
          Krhoe(i,j,k) = a9(1) * ( Frhoe(i,j+1,k)-Frhoe(i,j-1,k) ) &
                       + a9(2) * ( Frhoe(i,j+2,k)-Frhoe(i,j-2,k) ) &
                       + a9(3) * ( Frhoe(i,j+3,k)-Frhoe(i,j-3,k) ) &
                       + a9(4) * ( Frhoe(i,j+4,k)-Frhoe(i,j-4,k) ) + Krhoe(i,j,k)
       enddo
    enddo

    ! ~> Second pass: order reduction for metrics & 11-pt derivatives
    ! ===============================================================

    ! modified fluxes along ksi
    ! -------------------------
    do k=ndz_e,nfz_e
       do j=ny-4,ny
          do i=nx-9,nx+5
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
       do j=ny-4,ny
          do i=nx-4,nx
             Krho(i,j,k)  = a11(1) * ( Frho(i+1,j,k)-Frho(i-1,j,k) )   &
                          + a11(2) * ( Frho(i+2,j,k)-Frho(i-2,j,k) )   &
                          + a11(3) * ( Frho(i+3,j,k)-Frho(i-3,j,k) )   &
                          + a11(4) * ( Frho(i+4,j,k)-Frho(i-4,j,k) )   &
                          + a11(5) * ( Frho(i+5,j,k)-Frho(i-5,j,k) ) + Krho(i,j,k)

             Krhou(i,j,k) = a11(1) * ( Frhou(i+1,j,k)-Frhou(i-1,j,k) ) &
                          + a11(2) * ( Frhou(i+2,j,k)-Frhou(i-2,j,k) ) &
                          + a11(3) * ( Frhou(i+3,j,k)-Frhou(i-3,j,k) ) &
                          + a11(4) * ( Frhou(i+4,j,k)-Frhou(i-4,j,k) ) &
                          + a11(5) * ( Frhou(i+5,j,k)-Frhou(i-5,j,k) ) + Krhou(i,j,k)

             Krhov(i,j,k) = a11(1) * ( Frhov(i+1,j,k)-Frhov(i-1,j,k) ) &
                          + a11(2) * ( Frhov(i+2,j,k)-Frhov(i-2,j,k) ) &
                          + a11(3) * ( Frhov(i+3,j,k)-Frhov(i-3,j,k) ) &
                          + a11(4) * ( Frhov(i+4,j,k)-Frhov(i-4,j,k) ) &
                          + a11(5) * ( Frhov(i+5,j,k)-Frhov(i-5,j,k) ) + Krhov(i,j,k)

             Krhow(i,j,k) = a11(1) * ( Frhow(i+1,j,k)-Frhow(i-1,j,k) ) &
                          + a11(2) * ( Frhow(i+2,j,k)-Frhow(i-2,j,k) ) &
                          + a11(3) * ( Frhow(i+3,j,k)-Frhow(i-3,j,k) ) &
                          + a11(4) * ( Frhow(i+4,j,k)-Frhow(i-4,j,k) ) &
                          + a11(5) * ( Frhow(i+5,j,k)-Frhow(i-5,j,k) ) + Krhow(i,j,k)

             Krhoe(i,j,k) = a11(1) * ( Frhoe(i+1,j,k)-Frhoe(i-1,j,k) ) &
                          + a11(2) * ( Frhoe(i+2,j,k)-Frhoe(i-2,j,k) ) &
                          + a11(3) * ( Frhoe(i+3,j,k)-Frhoe(i-3,j,k) ) &
                          + a11(4) * ( Frhoe(i+4,j,k)-Frhoe(i-4,j,k) ) &
                          + a11(5) * ( Frhoe(i+5,j,k)-Frhoe(i-5,j,k) ) + Krhoe(i,j,k)
          enddo
       enddo
    enddo

    ! modified fluxes along eta
    ! -------------------------
    do k=ndz_e,nfz_e
       do i=nx-4,nx
          do j=ny-9,ny+5
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
       do j=ny-4,ny
          do i=nx-4,nx
             Krho(i,j,k)  = a11(1) * ( Frho(i,j+1,k)-Frho(i,j-1,k) ) &
                          + a11(2) * ( Frho(i,j+2,k)-Frho(i,j-2,k) ) &
                          + a11(3) * ( Frho(i,j+3,k)-Frho(i,j-3,k) ) &
                          + a11(4) * ( Frho(i,j+4,k)-Frho(i,j-4,k) ) &
                          + a11(5) * ( Frho(i,j+5,k)-Frho(i,j-5,k) ) &
                          + Krho(i,j,k)

             Krhou(i,j,k) = a11(1) * ( Frhou(i,j+1,k)-Frhou(i,j-1,k) ) &
                          + a11(2) * ( Frhou(i,j+2,k)-Frhou(i,j-2,k) ) &
                          + a11(3) * ( Frhou(i,j+3,k)-Frhou(i,j-3,k) ) &
                          + a11(4) * ( Frhou(i,j+4,k)-Frhou(i,j-4,k) ) &
                          + a11(5) * ( Frhou(i,j+5,k)-Frhou(i,j-5,k) ) &
                          + Krhou(i,j,k)

             Krhov(i,j,k) = a11(1) * ( Frhov(i,j+1,k)-Frhov(i,j-1,k) ) &
                          + a11(2) * ( Frhov(i,j+2,k)-Frhov(i,j-2,k) ) &
                          + a11(3) * ( Frhov(i,j+3,k)-Frhov(i,j-3,k) ) &
                          + a11(4) * ( Frhov(i,j+4,k)-Frhov(i,j-4,k) ) &
                          + a11(5) * ( Frhov(i,j+5,k)-Frhov(i,j-5,k) ) &
                          + Krhov(i,j,k)

             Krhow(i,j,k) = a11(1) * ( Frhow(i,j+1,k)-Frhow(i,j-1,k) ) &
                          + a11(2) * ( Frhow(i,j+2,k)-Frhow(i,j-2,k) ) &
                          + a11(3) * ( Frhow(i,j+3,k)-Frhow(i,j-3,k) ) &
                          + a11(4) * ( Frhow(i,j+4,k)-Frhow(i,j-4,k) ) &
                          + a11(5) * ( Frhow(i,j+5,k)-Frhow(i,j-5,k) ) &
                          + Krhow(i,j,k)

             Krhoe(i,j,k) = a11(1) * ( Frhoe(i,j+1,k)-Frhoe(i,j-1,k) ) &
                          + a11(2) * ( Frhoe(i,j+2,k)-Frhoe(i,j-2,k) ) &
                          + a11(3) * ( Frhoe(i,j+3,k)-Frhoe(i,j-3,k) ) &
                          + a11(4) * ( Frhoe(i,j+4,k)-Frhoe(i,j-4,k) ) &
                          + a11(5) * ( Frhoe(i,j+5,k)-Frhoe(i,j-5,k) ) &
                          + Krhoe(i,j,k)
          enddo
       enddo
    enddo

    ! Average of the two derivative evaluations
    ! =========================================
    do k=ndz_e,nfz_e
       do j=ny-4,ny
          do i=nx-4,nx
             Krho(i,j,k)  = Krho(i,j,k)*0.5_wp
             Krhou(i,j,k) = Krhou(i,j,k)*0.5_wp
             Krhov(i,j,k) = Krhov(i,j,k)*0.5_wp
             Krhow(i,j,k) = Krhow(i,j,k)*0.5_wp
             Krhoe(i,j,k) = Krhoe(i,j,k)*0.5_wp
          enddo
       enddo
    enddo

  end subroutine flux_euler_w_imax_jmax_11pts_SBP4_c

end submodule smod_flux_euler_wedge_11pts_SBP4_c
