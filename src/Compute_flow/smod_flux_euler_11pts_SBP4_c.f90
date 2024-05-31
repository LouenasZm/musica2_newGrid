!=================================================================================
submodule (mod_flux_euler) smod_flux_euler_11pts_SBP4_c
!=================================================================================
  !> author: XG
  !> date: January 2024
  !> 2.5D curvilinear version - compute Eulerian fluxes (inviscid part) with 11-pt stencil
  !> using SBP4 boundary schemes
!=================================================================================

contains

  !==============================================================================
  module subroutine flux_euler_11pts_SBP4_c
  !==============================================================================
    !> Derivatives of Eulerian fluxes (inviscid part) - 11-point stencil -
    !> - 2.5D curvilinear version with SBP4 boundary schemes -
 !==============================================================================
    !use mod_fluid
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: vc,rovc
    !real(wp) :: dvdy
    ! ---------------------------------------------------------------------------

    ! Computation of inviscid curvilinear fluxes along ksi
    ! ====================================================
    do k=ndz_e,nfz_e
       do j=ndy_e,nfy_e
          do i=ndxt,nfxt
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

    ! Flux derivatives along ksi
    ! ==========================

    ! Wall BC at imin
    ! ---------------
    if (is_bc_wall2(1,1)) then
       !i=1
       do k=ndz_e,nfz_e
          do j=ndy_e,nfy_e
             Krho(1,j,k) = as4p0(1)* Frho(1,j,k)+as4p0(2)* Frho(2,j,k) &
                         + as4p0(3)* Frho(3,j,k)+as4p0(4)* Frho(4,j,k)
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
    endif

    if (is_bc_1pt(1,1)) then
       !i=2
       do k=ndz_e,nfz_e
          do j=ndy_e,nfy_e
             Krho(2,j,k) = as4p1(1)* Frho(1,j,k)+as4p1(2)* Frho(2,j,k) &
                         + as4p1(3)* Frho(3,j,k)+as4p1(4)* Frho(4,j,k)
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

       !i=3
       do k=ndz_e,nfz_e
          do j=ndy_e,nfy_e
             Krho(3,j,k) = as4p2(1)* Frho(1,j,k)+as4p2(2)* Frho(2,j,k) &
                         + as4p2(3)* Frho(3,j,k)+as4p2(4)* Frho(4,j,k) &
                         + as4p2(5)* Frho(5,j,k)
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

       !i=4
       do k=ndz_e,nfz_e
          do j=ndy_e,nfy_e
             Krho(4,j,k) = as4p3(1)* Frho(1,j,k)+as4p3(2)* Frho(2,j,k) &
                         + as4p3(3)* Frho(3,j,k)+as4p3(4)* Frho(4,j,k) &
                         + as4p3(5)* Frho(5,j,k)+as4p3(6)* Frho(6,j,k)
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

       i=5
       do k=ndz_e,nfz_e
          do j=ndy_e,nfy_e
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
    endif

    ! Interior points
    ! ---------------
    do k=ndz_e,nfz_e
       do j=ndy_e,nfy_e
          do i=ndx,nfx
             Krho(i,j,k)  = a11(1) * ( Frho(i+1,j,k)-Frho(i-1,j,k) )   &
                          + a11(2) * ( Frho(i+2,j,k)-Frho(i-2,j,k) )   &
                          + a11(3) * ( Frho(i+3,j,k)-Frho(i-3,j,k) )   &
                          + a11(4) * ( Frho(i+4,j,k)-Frho(i-4,j,k) )   &
                          + a11(5) * ( Frho(i+5,j,k)-Frho(i-5,j,k) )

             Krhou(i,j,k) = a11(1) * ( Frhou(i+1,j,k)-Frhou(i-1,j,k) ) &
                          + a11(2) * ( Frhou(i+2,j,k)-Frhou(i-2,j,k) ) &
                          + a11(3) * ( Frhou(i+3,j,k)-Frhou(i-3,j,k) ) &
                          + a11(4) * ( Frhou(i+4,j,k)-Frhou(i-4,j,k) ) &
                          + a11(5) * ( Frhou(i+5,j,k)-Frhou(i-5,j,k) )

             Krhov(i,j,k) = a11(1) * ( Frhov(i+1,j,k)-Frhov(i-1,j,k) ) &
                          + a11(2) * ( Frhov(i+2,j,k)-Frhov(i-2,j,k) ) &
                          + a11(3) * ( Frhov(i+3,j,k)-Frhov(i-3,j,k) ) &
                          + a11(4) * ( Frhov(i+4,j,k)-Frhov(i-4,j,k) ) &
                          + a11(5) * ( Frhov(i+5,j,k)-Frhov(i-5,j,k) )

             Krhow(i,j,k) = a11(1) * ( Frhow(i+1,j,k)-Frhow(i-1,j,k) ) &
                          + a11(2) * ( Frhow(i+2,j,k)-Frhow(i-2,j,k) ) &
                          + a11(3) * ( Frhow(i+3,j,k)-Frhow(i-3,j,k) ) &
                          + a11(4) * ( Frhow(i+4,j,k)-Frhow(i-4,j,k) ) &
                          + a11(5) * ( Frhow(i+5,j,k)-Frhow(i-5,j,k) )

             Krhoe(i,j,k) = a11(1) * ( Frhoe(i+1,j,k)-Frhoe(i-1,j,k) ) &
                          + a11(2) * ( Frhoe(i+2,j,k)-Frhoe(i-2,j,k) ) &
                          + a11(3) * ( Frhoe(i+3,j,k)-Frhoe(i-3,j,k) ) &
                          + a11(4) * ( Frhoe(i+4,j,k)-Frhoe(i-4,j,k) ) &
                          + a11(5) * ( Frhoe(i+5,j,k)-Frhoe(i-5,j,k) )
          enddo
       enddo
    enddo

    ! BC at imax
    ! ----------
    if (is_bc_1pt(1,2)) then
       i=nx-4
       do k=ndz_e,nfz_e
          do j=ndy_e,nfy_e
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

       i=nx-3
       do k=ndz_e,nfz_e
          do j=ndy_e,nfy_e
             Krho(i,j,k) = as4m3(1)* Frho(nx  ,j,k)+as4m3(2)* Frho(nx-1,j,k) &
                         + as4m3(3)* Frho(nx-2,j,k)+as4m3(4)* Frho(nx-3,j,k) &
                         + as4m3(5)* Frho(nx-4,j,k)+as4m3(6)* Frho(nx-5,j,k)
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

       i=nx-2
       do k=ndz_e,nfz_e
          do j=ndy_e,nfy_e
             Krho(i,j,k) = as4m2(1)* Frho(nx  ,j,k)+as4m2(2)* Frho(nx-1,j,k) &
                         + as4m2(3)* Frho(nx-2,j,k)+as4m2(4)* Frho(nx-3,j,k) &
                         + as4m2(5)* Frho(nx-4,j,k)
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

       i=nx-1
       do k=ndz_e,nfz_e
          do j=ndy_e,nfy_e
             Krho(i,j,k) = as4m1(1)* Frho(nx  ,j,k)+as4m1(2)* Frho(nx-1,j,k) &
                         + as4m1(3)* Frho(nx-2,j,k)+as4m1(4)* Frho(nx-3,j,k)
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
    endif

    ! Wall BC at imax
    ! ---------------
    if (is_bc_wall2(1,2)) then
       i=nx
       do k=ndz_e,nfz_e
          do j=ndy_e,nfy_e
             Krho(i,j,k) = as4m0(1)* Frho(nx  ,j,k)+as4m0(2)* Frho(nx-1,j,k) &
                         + as4m0(3)* Frho(nx-2,j,k)+as4m0(4)* Frho(nx-3,j,k)
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
    endif

    ! Computation of curvilinear fluxes along eta
    ! ===========================================
    do k=ndz_e,nfz_e
       do j=ndyt,nfyt
          do i=ndx_e,nfx_e
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

    ! Flux derivatives along eta
    ! ==========================
    
    ! Wall BC at jmin
    ! ---------------
    if (is_bc_wall2(2,1)) then
       !j=1
       do k=ndz_e,nfz_e
          do i=ndx_e,nfx_e          
             Krho(i,1,k) = as4p0(1)* Frho(i,1,k)+as4p0(2)* Frho(i,2,k) &
                         + as4p0(3)* Frho(i,3,k)+as4p0(4)* Frho(i,4,k) + Krho(i,1,k)
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
!!$       ! simplified
!!$       j=1
!!$       do k=ndz_e,nfz_e
!!$          do i=ndx_e,nfx_e          
!!$             dvdy= (as4p0(1)*vv(i,1,k)+as4p0(2)*vv(i,2,k) &
!!$                  + as4p0(3)*vv(i,3,k)+as4p0(4)*vv(i,4,k))*idy1_jmin
!!$             Krho(i,j,k)=rho_n(i,j,k)*dvdy + Krho(i,j,k)
!!$             Krhou(i,j,k)=rhou_n(i,j,k)*dvdy + Krhou(i,j,k)
!!$             Krhoe(i,j,k)=(gam*igm1*prs(i,j,k)+0.5_wp*rhou_n(i,j,k)*uu(i,j,k))*dvdy &
!!$                  + Krhoe(i,j,k)
!!$          enddo
!!$       enddo
    endif

    if (is_bc_1pt(2,1)) then
       !j=2
       do k=ndz_e,nfz_e
          do i=ndx_e,nfx_e          
             Krho(i,2,k) = as4p1(1)* Frho(i,1,k)+as4p1(2)* Frho(i,2,k) &
                         + as4p1(3)* Frho(i,3,k)+as4p1(4)* Frho(i,4,k) + Krho(i,2,k)
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

       !j=3
       do k=ndz_e,nfz_e
          do i=ndx_e,nfx_e
             Krho(i,3,k) = as4p2(1)* Frho(i,1,k)+as4p2(2)* Frho(i,2,k) &
                         + as4p2(3)* Frho(i,3,k)+as4p2(4)* Frho(i,4,k) &
                         + as4p2(5)* Frho(i,5,k) + Krho(i,3,k)
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

       !j=4
       do k=ndz_e,nfz_e
          do i=ndx_e,nfx_e
             Krho(i,4,k) = as4p3(1)* Frho(i,1,k)+as4p3(2)* Frho(i,2,k) &
                         + as4p3(3)* Frho(i,3,k)+as4p3(4)* Frho(i,4,k) &
                         + as4p3(5)* Frho(i,5,k)+as4p3(6)* Frho(i,6,k) + Krho(i,4,k)
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

       j=5
       do k=ndz_e,nfz_e
          do i=ndx_e,nfx_e
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
    endif

    ! Interior points
    ! ---------------
    do k=ndz_e,nfz_e
       do j=ndy,nfy
          do i=ndx_e,nfx_e
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

    ! BC at jmax
    ! ----------
    if (is_bc_1pt(2,2)) then
       j=ny-4
       do k=ndz_e,nfz_e
          do i=ndx_e,nfx_e
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

       j=ny-3
       do k=ndz_e,nfz_e
          do i=ndx_e,nfx_e
             Krho(i,j,k) = as4m3(1)* Frho(i,ny  ,k)+as4m3(2)* Frho(i,ny-1,k) &
                         + as4m3(3)* Frho(i,ny-2,k)+as4m3(4)* Frho(i,ny-3,k) &
                         + as4m3(5)* Frho(i,ny-4,k)+as4m3(6)* Frho(i,ny-5,k) + Krho(i,j,k)
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

       j=ny-2
       do k=ndz_e,nfz_e
          do i=ndx_e,nfx_e
             Krho(i,j,k) = as4m2(1)* Frho(i,ny  ,k)+as4m2(2)* Frho(i,ny-1,k) &
                         + as4m2(3)* Frho(i,ny-2,k)+as4m2(4)* Frho(i,ny-3,k) &
                         + as4m2(5)* Frho(i,ny-4,k) + Krho(i,j,k)
             Krhou(i,j,k)= as4m2(1)*Frhou(i,ny  ,k)+as4m2(2)*Frhou(i,ny-1,k) &
                         + as4m2(3)*Frhou(i,ny-2,k)+as4m2(4)*Frhou(i,ny-3,k) &
                         + as4m2(5)*Frhou(i,ny-4,k) + Krhou(i,j,k)
             Krhov(i,j,k)= as4m2(1)*Frhov(i,ny  ,k)+as4m2(2)*Frhov(i,ny-1,k) &
                         + as4m2(3)*Frhov(i,ny-2,k)+as4m2(4)*Frhov(i,ny-3,k) &
                         + as4m2(5)*Frhov(i,ny-4,k) + Krhov(i,j,k)
             Krhow(i,j,k)= as4m2(1)*Frhow(i,ny  ,k)+as4m2(2)*Frhow(i,ny-1,k) &
                         + as4m2(3)*Frhow(i,ny-2,k)+as4m2(4)*Frhow(i,ny-3,k) &
                         + as4m2(5)*Frhow(i,ny-4,k) + Krhow(i,j,k)
             Krhoe(i,j,k)= as4m2(1)*Frhoe(i,ny  ,k)+as4m2(2)*Frhoe(i,ny-1,k) &
                         + as4m2(3)*Frhoe(i,ny-2,k)+as4m2(4)*Frhoe(i,ny-3,k) &
                         + as4m2(5)*Frhoe(i,ny-4,k) + Krhoe(i,j,k)
          enddo
       enddo

       j=ny-1
       do k=ndz_e,nfz_e
          do i=ndx_e,nfx_e
             Krho(i,j,k) = as4m1(1)* Frho(i,ny  ,k)+as4m1(2)* Frho(i,ny-1,k) &
                         + as4m1(3)* Frho(i,ny-2,k)+as4m1(4)* Frho(i,ny-3,k) + Krho(i,j,k)
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
    endif

    ! Wall BC at jmax
    ! ---------------
    if (is_bc_wall2(2,2)) then
       j=ny
       do k=ndz_e,nfz_e
          do i=ndx_e,nfx_e
             Krho(i,j,k) = as4m0(1)* Frho(i,ny  ,k)+as4m0(2)* Frho(i,ny-1,k) &
                         + as4m0(3)* Frho(i,ny-2,k)+as4m0(4)* Frho(i,ny-3,k) + Krho(i,j,k)
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
    endif

    ! Compute derivatives of Euler fluxes for edges near two walls
    ! ============================================================
    if (BC_edge(1,1,1)%sort==2) &
         call flux_euler_w_imin_jmin_11pts_SBP4_c
    if (BC_edge(1,1,2)%sort==2) &
         call flux_euler_w_imin_jmax_11pts_SBP4_c
    if (BC_edge(1,2,1)%sort==2) &
         call flux_euler_w_imax_jmin_11pts_SBP4_c
    if (BC_edge(1,2,2)%sort==2) &
         call flux_euler_w_imax_jmax_11pts_SBP4_c

    ! Multiply by inverse Jacobian
    ! ============================
    do k=ndz_e,nfz_e
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Krho(i,j,k)=  Krho(i,j,k)*ijacob(i,j)
             Krhou(i,j,k)= Krhou(i,j,k)*ijacob(i,j)
             Krhov(i,j,k)= Krhov(i,j,k)*ijacob(i,j)
             Krhow(i,j,k)= Krhow(i,j,k)*ijacob(i,j)
             Krhoe(i,j,k)= Krhoe(i,j,k)*ijacob(i,j)
          enddo
       enddo
    enddo

    !****************
    if (is_2D) return
    !****************

    ! Computation of inviscid fluxes along z
    ! ======================================
    do k=ndzt,nfzt
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Frho(i,j,k)  =  rhow_n(i,j,k)
             Frhou(i,j,k) =  rhow_n(i,j,k)*uu(i,j,k)
             Frhov(i,j,k) =  rhow_n(i,j,k)*vv(i,j,k)
             Frhow(i,j,k) =  rhow_n(i,j,k)*ww(i,j,k) +prs(i,j,k)
             Frhoe(i,j,k) = (rhoe_n(i,j,k)+prs(i,j,k))*ww(i,j,k)
          enddo
       enddo
    enddo

    ! Flux derivatives along z
    ! ========================

    ! Wall BC at kmin
    ! ---------------
    if (is_bc_wall2(3,1)) then
       !k=1
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e          
             Krho(i,j,1) = (as4p0(1)* Frho(i,j,1)+as4p0(2)* Frho(i,j,2) &
                          + as4p0(3)* Frho(i,j,3)+as4p0(4)* Frho(i,j,4))*idz1_kmin + Krho(i,j,1)
             Krhou(i,j,1)= (as4p0(1)*Frhou(i,j,1)+as4p0(2)*Frhou(i,j,2) &
                          + as4p0(3)*Frhou(i,j,3)+as4p0(4)*Frhou(i,j,4))*idz1_kmin + Krhou(i,j,1)
             Krhov(i,j,1)= (as4p0(1)*Frhov(i,j,1)+as4p0(2)*Frhov(i,j,2) &
                          + as4p0(3)*Frhov(i,j,3)+as4p0(4)*Frhov(i,j,4))*idz1_kmin + Krhov(i,j,1)
             Krhow(i,j,1)= (as4p0(1)*Frhow(i,j,1)+as4p0(2)*Frhow(i,j,2) &
                          + as4p0(3)*Frhow(i,j,3)+as4p0(4)*Frhow(i,j,4))*idz1_kmin + Krhow(i,j,1)
             Krhoe(i,j,1)= (as4p0(1)*Frhoe(i,j,1)+as4p0(2)*Frhoe(i,j,2) &
                          + as4p0(3)*Frhoe(i,j,3)+as4p0(4)*Frhoe(i,j,4))*idz1_kmin + Krhoe(i,j,1)

          enddo
       enddo
    endif

    ! BC at kmin
    ! ----------
    if (is_bc_1pt(3,1)) then
       !k=2
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e          
             Krho(i,j,2) = (as4p1(1)* Frho(i,j,1)+as4p1(2)* Frho(i,j,2) &
                          + as4p1(3)* Frho(i,j,3)+as4p1(4)* Frho(i,j,4))*idz2_kmin + Krho(i,j,2)
             Krhou(i,j,2)= (as4p1(1)*Frhou(i,j,1)+as4p1(2)*Frhou(i,j,2) &
                          + as4p1(3)*Frhou(i,j,3)+as4p1(4)*Frhou(i,j,4))*idz2_kmin + Krhou(i,j,2)
             Krhov(i,j,2)= (as4p1(1)*Frhov(i,j,1)+as4p1(2)*Frhov(i,j,2) &
                          + as4p1(3)*Frhov(i,j,3)+as4p1(4)*Frhov(i,j,4))*idz2_kmin + Krhov(i,j,2)
             Krhow(i,j,2)= (as4p1(1)*Frhow(i,j,1)+as4p1(2)*Frhow(i,j,2) &
                          + as4p1(3)*Frhow(i,j,3)+as4p1(4)*Frhow(i,j,4))*idz2_kmin + Krhow(i,j,2)
             Krhoe(i,j,2)= (as4p1(1)*Frhoe(i,j,1)+as4p1(2)*Frhoe(i,j,2) &
                          + as4p1(3)*Frhoe(i,j,3)+as4p1(4)*Frhoe(i,j,4))*idz2_kmin + Krhoe(i,j,2)
          enddo
       enddo

       !k=3
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Krho(i,j,3) = (as4p2(1)* Frho(i,j,1)+as4p2(2)* Frho(i,j,2) &
                          + as4p2(3)* Frho(i,j,3)+as4p2(4)* Frho(i,j,4) &
                          + as4p2(5)* Frho(i,j,5))*idz4_kmin + Krho(i,j,3)
             Krhou(i,j,3)= (as4p2(1)*Frhou(i,j,1)+as4p2(2)*Frhou(i,j,2) &
                          + as4p2(3)*Frhou(i,j,3)+as4p2(4)*Frhou(i,j,4) &
                          + as4p2(5)*Frhou(i,j,5))*idz4_kmin + Krhou(i,j,3)
             Krhov(i,j,3)= (as4p2(1)*Frhov(i,j,1)+as4p2(2)*Frhov(i,j,2) &
                          + as4p2(3)*Frhov(i,j,3)+as4p2(4)*Frhov(i,j,4) &
                          + as4p2(5)*Frhov(i,j,5))*idz4_kmin + Krhov(i,j,3)
             Krhow(i,j,3)= (as4p2(1)*Frhow(i,j,1)+as4p2(2)*Frhow(i,j,2) &
                          + as4p2(3)*Frhow(i,j,3)+as4p2(4)*Frhow(i,j,4) &
                          + as4p2(5)*Frhow(i,j,5))*idz4_kmin + Krhow(i,j,3)
             Krhoe(i,j,3)= (as4p2(1)*Frhoe(i,j,1)+as4p2(2)*Frhoe(i,j,2) &
                          + as4p2(3)*Frhoe(i,j,3)+as4p2(4)*Frhoe(i,j,4) &
                          + as4p2(5)*Frhoe(i,j,5))*idz4_kmin + Krhoe(i,j,3)
          enddo
       enddo

       !k=4
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Krho(i,j,4) = (as4p3(1)* Frho(i,j,1)+as4p3(2)* Frho(i,j,2) &
                          + as4p3(3)* Frho(i,j,3)+as4p3(4)* Frho(i,j,4) &
                          + as4p3(5)* Frho(i,j,5)+as4p3(6)* Frho(i,j,6))*idz6_kmin + Krho(i,j,4)
             Krhou(i,j,4)= (as4p3(1)*Frhou(i,j,1)+as4p3(2)*Frhou(i,j,2) &
                          + as4p3(3)*Frhou(i,j,3)+as4p3(4)*Frhou(i,j,4) &
                          + as4p3(5)*Frhou(i,j,5)+as4p3(6)*Frhou(i,j,6))*idz6_kmin + Krhou(i,j,4)
             Krhov(i,j,4)= (as4p3(1)*Frhov(i,j,1)+as4p3(2)*Frhov(i,j,2) &
                          + as4p3(3)*Frhov(i,j,3)+as4p3(4)*Frhov(i,j,4) &
                          + as4p3(5)*Frhov(i,j,5)+as4p3(6)*Frhov(i,j,6))*idz6_kmin + Krhov(i,j,4)
             Krhow(i,j,4)= (as4p3(1)*Frhow(i,j,1)+as4p3(2)*Frhow(i,j,2) &
                          + as4p3(3)*Frhow(i,j,3)+as4p3(4)*Frhow(i,j,4) &
                          + as4p3(5)*Frhow(i,j,5)+as4p3(6)*Frhow(i,j,6))*idz6_kmin + Krhow(i,j,4)
             Krhoe(i,j,4)= (as4p3(1)*Frhoe(i,j,1)+as4p3(2)*Frhoe(i,j,2) &
                          + as4p3(3)*Frhoe(i,j,3)+as4p3(4)*Frhoe(i,j,4) &
                          + as4p3(5)*Frhoe(i,j,5)+as4p3(6)*Frhoe(i,j,6))*idz6_kmin + Krhoe(i,j,4)
          enddo
       enddo

       k=5
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Krho(i,j,k) = ( a9(1) * ( Frho(i,j,k+1)-Frho(i,j,k-1) ) &
                           + a9(2) * ( Frho(i,j,k+2)-Frho(i,j,k-2) ) &
                           + a9(3) * ( Frho(i,j,k+3)-Frho(i,j,k-3) ) &
                           + a9(4) * ( Frho(i,j,k+4)-Frho(i,j,k-4) ) )*idz8_kmin &
                           + Krho(i,j,k)
             Krhou(i,j,k)= ( a9(1) * ( Frhou(i,j,k+1)-Frhou(i,j,k-1) ) &
                           + a9(2) * ( Frhou(i,j,k+2)-Frhou(i,j,k-2) ) &
                           + a9(3) * ( Frhou(i,j,k+3)-Frhou(i,j,k-3) ) &
                           + a9(4) * ( Frhou(i,j,k+4)-Frhou(i,j,k-4) ) )*idz8_kmin &
                           + Krhou(i,j,k)
             Krhov(i,j,k)= ( a9(1) * ( Frhov(i,j,k+1)-Frhov(i,j,k-1) ) &
                           + a9(2) * ( Frhov(i,j,k+2)-Frhov(i,j,k-2) ) &
                           + a9(3) * ( Frhov(i,j,k+3)-Frhov(i,j,k-3) ) &
                           + a9(4) * ( Frhov(i,j,k+4)-Frhov(i,j,k-4) ) )*idz8_kmin &
                           + Krhov(i,j,k)
             Krhow(i,j,k)= ( a9(1) * ( Frhow(i,j,k+1)-Frhow(i,j,k-1) ) &
                           + a9(2) * ( Frhow(i,j,k+2)-Frhow(i,j,k-2) ) &
                           + a9(3) * ( Frhow(i,j,k+3)-Frhow(i,j,k-3) ) &
                           + a9(4) * ( Frhow(i,j,k+4)-Frhow(i,j,k-4) ) )*idz8_kmin &
                           + Krhow(i,j,k)
             Krhoe(i,j,k)= ( a9(1) * ( Frhoe(i,j,k+1)-Frhoe(i,j,k-1) ) &
                           + a9(2) * ( Frhoe(i,j,k+2)-Frhoe(i,j,k-2) ) &
                           + a9(3) * ( Frhoe(i,j,k+3)-Frhoe(i,j,k-3) ) &
                           + a9(4) * ( Frhoe(i,j,k+4)-Frhoe(i,j,k-4) ) )*idz8_kmin &
                           + Krhoe(i,j,k)
          enddo
       enddo
    endif

    ! Interior points
    ! ---------------
    do k=ndz,nfz
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Krho(i,j,k) = ( a11(1)* ( Frho(i,j,k+1)-Frho(i,j,k-1) ) &
                           + a11(2)* ( Frho(i,j,k+2)-Frho(i,j,k-2) ) &
                           + a11(3)* ( Frho(i,j,k+3)-Frho(i,j,k-3) ) &
                           + a11(4)* ( Frho(i,j,k+4)-Frho(i,j,k-4) ) &
                           + a11(5)* ( Frho(i,j,k+5)-Frho(i,j,k-5) ) )*idz(k) &
                           + Krho(i,j,k)

             Krhou(i,j,k)= ( a11(1)* ( Frhou(i,j,k+1)-Frhou(i,j,k-1) ) &
                           + a11(2)* ( Frhou(i,j,k+2)-Frhou(i,j,k-2) ) &
                           + a11(3)* ( Frhou(i,j,k+3)-Frhou(i,j,k-3) ) &
                           + a11(4)* ( Frhou(i,j,k+4)-Frhou(i,j,k-4) ) &
                           + a11(5)* ( Frhou(i,j,k+5)-Frhou(i,j,k-5) ) )*idz(k) &
                           + Krhou(i,j,k)

             Krhov(i,j,k)= ( a11(1)* ( Frhov(i,j,k+1)-Frhov(i,j,k-1) ) &
                           + a11(2)* ( Frhov(i,j,k+2)-Frhov(i,j,k-2) ) &
                           + a11(3)* ( Frhov(i,j,k+3)-Frhov(i,j,k-3) ) &
                           + a11(4)* ( Frhov(i,j,k+4)-Frhov(i,j,k-4) ) &
                           + a11(5)* ( Frhov(i,j,k+5)-Frhov(i,j,k-5) ) )*idz(k) &
                           + Krhov(i,j,k)

             Krhow(i,j,k)= ( a11(1)* ( Frhow(i,j,k+1)-Frhow(i,j,k-1) ) &
                           + a11(2)* ( Frhow(i,j,k+2)-Frhow(i,j,k-2) ) &
                           + a11(3)* ( Frhow(i,j,k+3)-Frhow(i,j,k-3) ) &
                           + a11(4)* ( Frhow(i,j,k+4)-Frhow(i,j,k-4) ) &
                           + a11(5)* ( Frhow(i,j,k+5)-Frhow(i,j,k-5) ) )*idz(k) &
                           + Krhow(i,j,k)

             Krhoe(i,j,k)= ( a11(1)* ( Frhoe(i,j,k+1)-Frhoe(i,j,k-1) ) &
                           + a11(2)* ( Frhoe(i,j,k+2)-Frhoe(i,j,k-2) ) &
                           + a11(3)* ( Frhoe(i,j,k+3)-Frhoe(i,j,k-3) ) &
                           + a11(4)* ( Frhoe(i,j,k+4)-Frhoe(i,j,k-4) ) &
                           + a11(5)* ( Frhoe(i,j,k+5)-Frhoe(i,j,k-5) ) )*idz(k) &
                           + Krhoe(i,j,k)
          enddo
       enddo
    enddo

    ! BC at kmax
    ! ----------
    if (is_bc_1pt(3,2)) then
       k=nz-4
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Krho(i,j,k) = ( a9(1) * ( Frho(i,j,k+1)-Frho(i,j,k-1) ) &
                           + a9(2) * ( Frho(i,j,k+2)-Frho(i,j,k-2) ) &
                           + a9(3) * ( Frho(i,j,k+3)-Frho(i,j,k-3) ) &
                           + a9(4) * ( Frho(i,j,k+4)-Frho(i,j,k-4) ) )*idz8_kmax &
                           + Krho(i,j,k)
             Krhou(i,j,k)= ( a9(1) * ( Frhou(i,j,k+1)-Frhou(i,j,k-1) ) &
                           + a9(2) * ( Frhou(i,j,k+2)-Frhou(i,j,k-2) ) &
                           + a9(3) * ( Frhou(i,j,k+3)-Frhou(i,j,k-3) ) &
                           + a9(4) * ( Frhou(i,j,k+4)-Frhou(i,j,k-4) ) )*idz8_kmax &
                           + Krhou(i,j,k)
             Krhov(i,j,k)= ( a9(1) * ( Frhov(i,j,k+1)-Frhov(i,j,k-1) ) &
                           + a9(2) * ( Frhov(i,j,k+2)-Frhov(i,j,k-2) ) &
                           + a9(3) * ( Frhov(i,j,k+3)-Frhov(i,j,k-3) ) &
                           + a9(4) * ( Frhov(i,j,k+4)-Frhov(i,j,k-4) ) )*idz8_kmax &
                           + Krhov(i,j,k)
             Krhow(i,j,k)= ( a9(1) * ( Frhow(i,j,k+1)-Frhow(i,j,k-1) ) &
                           + a9(2) * ( Frhow(i,j,k+2)-Frhow(i,j,k-2) ) &
                           + a9(3) * ( Frhow(i,j,k+3)-Frhow(i,j,k-3) ) &
                           + a9(4) * ( Frhow(i,j,k+4)-Frhow(i,j,k-4) ) )*idz8_kmax &
                           + Krhow(i,j,k)
             Krhoe(i,j,k)= ( a9(1) * ( Frhoe(i,j,k+1)-Frhoe(i,j,k-1) ) &
                           + a9(2) * ( Frhoe(i,j,k+2)-Frhoe(i,j,k-2) ) &
                           + a9(3) * ( Frhoe(i,j,k+3)-Frhoe(i,j,k-3) ) &
                           + a9(4) * ( Frhoe(i,j,k+4)-Frhoe(i,j,k-4) ) )*idz8_kmax &
                           + Krhoe(i,j,k)
          enddo
       enddo

       k=nz-3
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Krho(i,j,k) = (as4m3(1)* Frho(i,j,nz  )+as4m3(2)* Frho(i,j,nz-1) &
                          + as4m3(3)* Frho(i,j,nz-2)+as4m3(4)* Frho(i,j,nz-3) &
                          + as4m3(5)* Frho(i,j,nz-4)+as4m3(6)* Frho(i,j,nz-5))*idz6_kmax + Krho(i,j,k)
             Krhou(i,j,k)= (as4m3(1)*Frhou(i,j,nz  )+as4m3(2)*Frhou(i,j,nz-1) &
                          + as4m3(3)*Frhou(i,j,nz-2)+as4m3(4)*Frhou(i,j,nz-3) &
                          + as4m3(5)*Frhou(i,j,nz-4)+as4m3(6)*Frhou(i,j,nz-5))*idz6_kmax + Krhou(i,j,k)
             Krhov(i,j,k)= (as4m3(1)*Frhov(i,j,nz  )+as4m3(2)*Frhov(i,j,nz-1) &
                          + as4m3(3)*Frhov(i,j,nz-2)+as4m3(4)*Frhov(i,j,nz-3) &
                          + as4m3(5)*Frhov(i,j,nz-4)+as4m3(6)*Frhov(i,j,nz-5))*idz6_kmax + Krhov(i,j,k)
             Krhow(i,j,k)= (as4m3(1)*Frhow(i,j,nz  )+as4m3(2)*Frhow(i,j,nz-1) &
                          + as4m3(3)*Frhow(i,j,nz-2)+as4m3(4)*Frhow(i,j,nz-3) &
                          + as4m3(5)*Frhow(i,j,nz-4)+as4m3(6)*Frhow(i,j,nz-5))*idz6_kmax + Krhow(i,j,k)
             Krhoe(i,j,k)= (as4m3(1)*Frhoe(i,j,nz  )+as4m3(2)*Frhoe(i,j,nz-1) &
                          + as4m3(3)*Frhoe(i,j,nz-2)+as4m3(4)*Frhoe(i,j,nz-3) &
                          + as4m3(5)*Frhoe(i,j,nz-4)+as4m3(6)*Frhoe(i,j,nz-5))*idz6_kmax + Krhoe(i,j,k)
          enddo
       enddo

       k=nz-2
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Krho(i,j,k) = (as4m2(1)* Frho(i,j,nz  )+as4m2(2)* Frho(i,j,nz-1) &
                          + as4m2(3)* Frho(i,j,nz-2)+as4m2(4)* Frho(i,j,nz-3) &
                          + as4m2(5)* Frho(i,j,nz-4))*idz4_kmax + Krho(i,j,k)
             Krhou(i,j,k)= (as4m2(1)*Frhou(i,j,nz  )+as4m2(2)*Frhou(i,j,nz-1) &
                          + as4m2(3)*Frhou(i,j,nz-2)+as4m2(4)*Frhou(i,j,nz-3) &
                          + as4m2(5)*Frhou(i,j,nz-4))*idz4_kmax + Krhou(i,j,k)
             Krhov(i,j,k)= (as4m2(1)*Frhov(i,j,nz  )+as4m2(2)*Frhov(i,j,nz-1) &
                          + as4m2(3)*Frhov(i,j,nz-2)+as4m2(4)*Frhov(i,j,nz-3) &
                          + as4m2(5)*Frhov(i,j,nz-4))*idz4_kmax + Krhov(i,j,k)
             Krhow(i,j,k)= (as4m2(1)*Frhow(i,j,nz  )+as4m2(2)*Frhow(i,j,nz-1) &
                          + as4m2(3)*Frhow(i,j,nz-2)+as4m2(4)*Frhow(i,j,nz-3) &
                          + as4m2(5)*Frhow(i,j,nz-4))*idz4_kmax + Krhow(i,j,k)
             Krhoe(i,j,k)= (as4m2(1)*Frhoe(i,j,nz  )+as4m2(2)*Frhoe(i,j,nz-1) &
                          + as4m2(3)*Frhoe(i,j,nz-2)+as4m2(4)*Frhoe(i,j,nz-3) &
                          + as4m2(5)*Frhoe(i,j,nz-4))*idz4_kmax + Krhoe(i,j,k)
          enddo
       enddo

       k=nz-1
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Krho(i,j,k) = (as4m1(1)* Frho(i,j,nz  )+as4m1(2)* Frho(i,j,nz-1) &
                          + as4m1(3)* Frho(i,j,nz-2)+as4m1(4)* Frho(i,j,nz-3))*idz2_kmax + Krho(i,j,k)
             Krhou(i,j,k)= (as4m1(1)*Frhou(i,j,nz  )+as4m1(2)*Frhou(i,j,nz-1) &
                          + as4m1(3)*Frhou(i,j,nz-2)+as4m1(4)*Frhou(i,j,nz-3))*idz2_kmax + Krhou(i,j,k)
             Krhov(i,j,k)= (as4m1(1)*Frhov(i,j,nz  )+as4m1(2)*Frhov(i,j,nz-1) &
                          + as4m1(3)*Frhov(i,j,nz-2)+as4m1(4)*Frhov(i,j,nz-3))*idz2_kmax + Krhov(i,j,k)
             Krhov(i,j,k)= (as4m1(1)*Frhov(i,j,nz  )+as4m1(2)*Frhov(i,j,nz-1) &
                          + as4m1(3)*Frhov(i,j,nz-2)+as4m1(4)*Frhov(i,j,nz-3))*idz2_kmax + Krhov(i,j,k)
             Krhoe(i,j,k)= (as4m1(1)*Frhoe(i,j,nz  )+as4m1(2)*Frhoe(i,j,nz-1) &
                          + as4m1(3)*Frhoe(i,j,nz-2)+as4m1(4)*Frhoe(i,j,nz-3))*idz2_kmax + Krhoe(i,j,k)
          enddo
       enddo
    endif

    ! Wall BC at kmax
    ! ---------------
    if (is_bc_wall2(3,2)) then
       k=nz
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Krho(i,j,k) = (as4m0(1)* Frho(i,j,nz  )+as4m0(2)* Frho(i,j,nz-1) &
                          + as4m0(3)* Frho(i,j,nz-2)+as4m0(4)* Frho(i,j,nz-3))*idz1_kmax + Krho(i,j,k)
             Krhou(i,j,k)= (as4m0(1)*Frhou(i,j,nz  )+as4m0(2)*Frhou(i,j,nz-1) &
                          + as4m0(3)*Frhou(i,j,nz-2)+as4m0(4)*Frhou(i,j,nz-3))*idz1_kmax + Krhou(i,j,k)
             Krhov(i,j,k)= (as4m0(1)*Frhov(i,j,nz  )+as4m0(2)*Frhov(i,j,nz-1) &
                          + as4m0(3)*Frhov(i,j,nz-2)+as4m0(4)*Frhov(i,j,nz-3))*idz1_kmax + Krhov(i,j,k)
             Krhov(i,j,k)= (as4m0(1)*Frhov(i,j,nz  )+as4m0(2)*Frhov(i,j,nz-1) &
                          + as4m0(3)*Frhov(i,j,nz-2)+as4m0(4)*Frhov(i,j,nz-3))*idz1_kmax + Krhov(i,j,k)
             Krhoe(i,j,k)= (as4m0(1)*Frhoe(i,j,nz  )+as4m0(2)*Frhoe(i,j,nz-1) &
                          + as4m0(3)*Frhoe(i,j,nz-2)+as4m0(4)*Frhoe(i,j,nz-3))*idz1_kmax + Krhoe(i,j,k)
          enddo
       enddo
    endif
    
  end subroutine flux_euler_11pts_SBP4_c

end submodule smod_flux_euler_11pts_SBP4_c
