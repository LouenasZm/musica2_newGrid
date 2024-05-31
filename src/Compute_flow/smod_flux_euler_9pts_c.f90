!=================================================================================
submodule (mod_flux_euler) smod_flux_euler_9pts_c
!=================================================================================
  !> author: XG
  !> date: January 2024
  !> 2.5D curvilinear version - compute Eulerian fluxes (inviscid part) with 9-pt stencil
!=================================================================================

contains

  !==============================================================================
  module subroutine flux_euler_9pts_c
  !==============================================================================
    !> Derivatives of Eulerian fluxes (inviscid part) - 9-point stencil -
    !> - 2.5D curvilinear version -
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: vc,rovc
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

    ! Advancement of wall points
    ! --------------------------
    if (is_bc_wall2(1,1)) then
       !i=1
       do k=ndz_e,nfz_e
          do j=ndy_e,nfy_e
             Krho(1,j,k)= a02(1)*Frho(1,j,k) +a02(2)*Frho(2,j,k) +a02(3)*Frho(3,j,k)
             Krhou(1,j,k)=a02(1)*Frhou(1,j,k)+a02(2)*Frhou(2,j,k)+a02(3)*Frhou(3,j,k)
             Krhov(1,j,k)=a02(1)*Frhov(1,j,k)+a02(2)*Frhov(2,j,k)+a02(3)*Frhov(3,j,k)
             Krhow(1,j,k)=a02(1)*Frhow(1,j,k)+a02(2)*Frhow(2,j,k)+a02(3)*Frhow(3,j,k)
             Krhoe(1,j,k)=a02(1)*Frhoe(1,j,k)+a02(2)*Frhoe(2,j,k)+a02(3)*Frhoe(3,j,k)
          enddo
       enddo
    endif

    ! BC at imin
    ! ----------
    if (is_bc_1pt(1,1)) then
       i=2
       do k=ndz_e,nfz_e
          do j=ndy_e,nfy_e
             Krho(i,j,k)  = 0.5_wp*( Frho(i+1,j,k)- Frho(i-1,j,k))
             Krhou(i,j,k) = 0.5_wp*(Frhou(i+1,j,k)-Frhou(i-1,j,k))
             Krhov(i,j,k) = 0.5_wp*(Frhov(i+1,j,k)-Frhov(i-1,j,k))
             Krhow(i,j,k) = 0.5_wp*(Frhow(i+1,j,k)-Frhow(i-1,j,k))
             Krhoe(i,j,k) = 0.5_wp*(Frhoe(i+1,j,k)-Frhoe(i-1,j,k))
          enddo
       enddo

       i=3
       do k=ndz_e,nfz_e
          do j=ndy_e,nfy_e
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
       do k=ndz_e,nfz_e
          do j=ndy_e,nfy_e
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
    endif

    ! Interior points
    ! ---------------
    do k=ndz_e,nfz_e
       do j=ndy_e,nfy_e
          do i=ndx,nfx
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
    enddo

    ! BC at imax
    ! ----------
    if (is_bc_1pt(1,2)) then
       i=nx-3
       do k=ndz_e,nfz_e
          do j=ndy_e,nfy_e
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
       do k=ndz_e,nfz_e
          do j=ndy_e,nfy_e
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
       do k=ndz_e,nfz_e
          do j=ndy_e,nfy_e
             Krho(i,j,k)  = 0.5_wp*( Frho(i+1,j,k)- Frho(i-1,j,k))
             Krhou(i,j,k) = 0.5_wp*(Frhou(i+1,j,k)-Frhou(i-1,j,k))
             Krhov(i,j,k) = 0.5_wp*(Frhov(i+1,j,k)-Frhov(i-1,j,k))
             Krhow(i,j,k) = 0.5_wp*(Frhow(i+1,j,k)-Frhow(i-1,j,k))
             Krhoe(i,j,k) = 0.5_wp*(Frhoe(i+1,j,k)-Frhoe(i-1,j,k))
          enddo
       enddo
    endif

    ! Advancement of wall points
    ! --------------------------
    if (is_bc_wall2(1,2)) then
       i=nx
       do k=ndz_e,nfz_e
          do j=ndy_e,nfy_e
             Krho(i,j,k) =a20(3)*Frho(i-2,j,k) +a20(2)*Frho(i-1,j,k) +a20(1)*Frho(i,j,k)
             Krhou(i,j,k)=a20(3)*Frhou(i-2,j,k)+a20(2)*Frhou(i-1,j,k)+a20(1)*Frhou(i,j,k)
             Krhov(i,j,k)=a20(3)*Frhov(i-2,j,k)+a20(2)*Frhov(i-1,j,k)+a20(1)*Frhov(i,j,k)
             Krhow(i,j,k)=a20(3)*Frhow(i-2,j,k)+a20(2)*Frhow(i-1,j,k)+a20(1)*Frhow(i,j,k)
             Krhoe(i,j,k)=a20(3)*Frhoe(i-2,j,k)+a20(2)*Frhoe(i-1,j,k)+a20(1)*Frhoe(i,j,k)
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

    ! Advancement of wall points
    ! --------------------------
    if (is_bc_wall2(2,1)) then
       !j=1
       do k=ndz_e,nfz_e
          do i=ndx_e,nfx_e
             Krho(i,1,k)= a02(1)*Frho(i,1,k) +a02(2)*Frho(i,2,k) +a02(3)*Frho(i,3,k) + Krho(i,1,k)
             Krhou(i,1,k)=a02(1)*Frhou(i,1,k)+a02(2)*Frhou(i,2,k)+a02(3)*Frhou(i,3,k)+ Krhou(i,1,k)
             Krhov(i,1,k)=a02(1)*Frhov(i,1,k)+a02(2)*Frhov(i,2,k)+a02(3)*Frhov(i,3,k)+ Krhov(i,1,k)
             Krhow(i,1,k)=a02(1)*Frhow(i,1,k)+a02(2)*Frhow(i,2,k)+a02(3)*Frhow(i,3,k)+ Krhow(i,1,k)
             Krhoe(i,1,k)=a02(1)*Frhoe(i,1,k)+a02(2)*Frhoe(i,2,k)+a02(3)*Frhoe(i,3,k)+ Krhoe(i,1,k)
          enddo
       enddo
    endif

    !BC at jmin
    ! ---------
    if (is_bc_1pt(2,1)) then
       j=2
       do k=ndz_e,nfz_e
          do i=ndx_e,nfx_e
             Krho(i,j,k)  = 0.5_wp*(Frho(i,j+1,k) -Frho(i,j-1,k)) + Krho(i,j,k)
             Krhou(i,j,k) = 0.5_wp*(Frhou(i,j+1,k)-Frhou(i,j-1,k)) + Krhou(i,j,k)
             Krhov(i,j,k) = 0.5_wp*(Frhov(i,j+1,k)-Frhov(i,j-1,k)) + Krhov(i,j,k)
             Krhow(i,j,k) = 0.5_wp*(Frhow(i,j+1,k)-Frhow(i,j-1,k)) + Krhow(i,j,k)
             Krhoe(i,j,k) = 0.5_wp*(Frhoe(i,j+1,k)-Frhoe(i,j-1,k)) + Krhoe(i,j,k)
          enddo
       enddo

       j=3
       do k=ndz_e,nfz_e
          do i=ndx_e,nfx_e
             Krho(i,j,k) =  a5(1)*(Frho(i,j+1,k)-Frho(i,j-1,k)) &
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
       do k=ndz_e,nfz_e
          do i=ndx_e,nfx_e
             Krho(i,j,k) =  a7(1) * ( Frho(i,j+1,k)-Frho(i,j-1,k) ) &
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
    endif

    ! Interior points
    ! ---------------
    do k=ndz_e,nfz_e
       do j=ndy,nfy
          do i=ndx_e,nfx_e
             Krho(i,j,k) =  a9(1) * ( Frho(i,j+1,k)-Frho(i,j-1,k) ) &
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
    enddo

    ! BC at jmax
    ! ----------
    if (is_bc_1pt(2,2)) then
       j=ny-3
       do k=ndz_e,nfz_e
          do i=ndx_e,nfx_e
             Krho(i,j,k) =  a7(1) * ( Frho(i,j+1,k)-Frho(i,j-1,k) ) &
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

       j=ny-2
       do k=ndz_e,nfz_e
          do i=ndx_e,nfx_e
             Krho(i,j,k) =  a5(1)*(Frho(i,j+1,k)-Frho(i,j-1,k)) &
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

       j=ny-1
       do k=ndz_e,nfz_e
          do i=ndx_e,nfx_e
             Krho(i,j,k)  = 0.5_wp*(Frho(i,j+1,k) -Frho(i,j-1,k)) + Krho(i,j,k)
             Krhou(i,j,k) = 0.5_wp*(Frhou(i,j+1,k)-Frhou(i,j-1,k)) + Krhou(i,j,k)
             Krhov(i,j,k) = 0.5_wp*(Frhov(i,j+1,k)-Frhov(i,j-1,k)) + Krhov(i,j,k)
             Krhow(i,j,k) = 0.5_wp*(Frhow(i,j+1,k)-Frhow(i,j-1,k)) + Krhow(i,j,k)
             Krhoe(i,j,k) = 0.5_wp*(Frhoe(i,j+1,k)-Frhoe(i,j-1,k)) + Krhoe(i,j,k)
          enddo
       enddo
    endif

    ! Advancement of wall points
    ! --------------------------
    if (is_bc_wall2(2,2)) then
       j=ny
       do k=ndz_e,nfz_e
          do i=ndx_e,nfx_e
             Krho(i,j,k) =a20(3)*Frho(i,j-2,k) +a20(2)*Frho(i,j-1,k) +a20(1)*Frho(i,j,k) + Krho(i,j,k)
             Krhou(i,j,k)=a20(3)*Frhou(i,j-2,k)+a20(2)*Frhou(i,j-1,k)+a20(1)*Frhou(i,j,k)+ Krhou(i,j,k)
             Krhov(i,j,k)=a20(3)*Frhov(i,j-2,k)+a20(2)*Frhov(i,j-1,k)+a20(1)*Frhov(i,j,k)+ Krhov(i,j,k)
             Krhow(i,j,k)=a20(3)*Frhow(i,j-2,k)+a20(2)*Frhow(i,j-1,k)+a20(1)*Frhow(i,j,k)+ Krhow(i,j,k)
             Krhoe(i,j,k)=a20(3)*Frhoe(i,j-2,k)+a20(2)*Frhoe(i,j-1,k)+a20(1)*Frhoe(i,j,k)+ Krhoe(i,j,k)
          enddo
       enddo
    endif

    ! Compute derivatives of Euler fluxes for edges near two walls
    ! ============================================================
    if (BC_edge(1,1,1)%sort==2) call flux_euler_w_imin_jmin_9pts_c
    if (BC_edge(1,1,2)%sort==2) call flux_euler_w_imin_jmax_9pts_c
    if (BC_edge(1,2,1)%sort==2) call flux_euler_w_imax_jmin_9pts_c
    if (BC_edge(1,2,2)%sort==2) call flux_euler_w_imax_jmax_9pts_c

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

    ! Advancement of wall points
    ! --------------------------
    if (is_bc_wall2(3,1)) then
       !k=1
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Krho(i,j,1)= (a02(1)*Frho(i,j,1) +a02(2)*Frho(i,j,2) +a02(3)*Frho(i,j,3) )*idz1_kmin+ Krho(i,j,1)
             Krhou(i,j,1)=(a02(1)*Frhou(i,j,1)+a02(2)*Frhou(i,j,2)+a02(3)*Frhou(i,j,3))*idz1_kmin+ Krhou(i,j,1)
             Krhov(i,j,1)=(a02(1)*Frhov(i,j,1)+a02(2)*Frhov(i,j,2)+a02(3)*Frhov(i,j,3))*idz1_kmin+ Krhov(i,j,1)
             Krhow(i,j,1)=(a02(1)*Frhow(i,j,1)+a02(2)*Frhow(i,j,2)+a02(3)*Frhow(i,j,3))*idz1_kmin+ Krhow(i,j,1)
             Krhoe(i,j,1)=(a02(1)*Frhoe(i,j,1)+a02(2)*Frhoe(i,j,2)+a02(3)*Frhoe(i,j,3))*idz1_kmin+ Krhoe(i,j,1)
          enddo
       enddo
    endif

    ! BC at kmin
    ! ----------
    if (is_bc_1pt(3,1)) then
       k=2
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Krho(i,j,k)  =  Krho(i,j,k) + ( Frho(i,j,k+1)- Frho(i,j,k-1))*idz2_kmin
             Krhou(i,j,k) = Krhou(i,j,k) + (Frhou(i,j,k+1)-Frhou(i,j,k-1))*idz2_kmin
             Krhov(i,j,k) = Krhov(i,j,k) + (Frhov(i,j,k+1)-Frhov(i,j,k-1))*idz2_kmin
             Krhow(i,j,k) = Krhow(i,j,k) + (Frhow(i,j,k+1)-Frhow(i,j,k-1))*idz2_kmin
             Krhoe(i,j,k) = Krhoe(i,j,k) + (Frhoe(i,j,k+1)-Frhoe(i,j,k-1))*idz2_kmin
          enddo
       enddo

       k=3
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Krho(i,j,k) = ( a5(1)*( Frho(i,j,k+1)- Frho(i,j,k-1)) &
                           + a5(2)*( Frho(i,j,k+2)- Frho(i,j,k-2)) )*idz4_kmin + Krho(i,j,k)
             Krhou(i,j,k)= ( a5(1)*(Frhou(i,j,k+1)-Frhou(i,j,k-1)) &
                           + a5(2)*(Frhou(i,j,k+2)-Frhou(i,j,k-2)) )*idz4_kmin + Krhou(i,j,k)
             Krhov(i,j,k)= ( a5(1)*(Frhov(i,j,k+1)-Frhov(i,j,k-1)) &
                           + a5(2)*(Frhov(i,j,k+2)-Frhov(i,j,k-2)) )*idz4_kmin + Krhov(i,j,k)
             Krhow(i,j,k)= ( a5(1)*(Frhow(i,j,k+1)-Frhow(i,j,k-1)) &
                           + a5(2)*(Frhow(i,j,k+2)-Frhow(i,j,k-2)) )*idz4_kmin + Krhow(i,j,k)
             Krhoe(i,j,k)= ( a5(1)*(Frhoe(i,j,k+1)-Frhoe(i,j,k-1)) &
                           + a5(2)*(Frhoe(i,j,k+2)-Frhoe(i,j,k-2)) )*idz4_kmin + Krhoe(i,j,k)
          enddo
       enddo

       k=4
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Krho(i,j,k) = ( a7(1) * ( Frho(i,j,k+1)-Frho(i,j,k-1) ) &
                           + a7(2) * ( Frho(i,j,k+2)-Frho(i,j,k-2) ) &
                           + a7(3) * ( Frho(i,j,k+3)-Frho(i,j,k-3) ) )*idz6_kmin &
                           + Krho(i,j,k)
             Krhou(i,j,k)= ( a7(1) * ( Frhou(i,j,k+1)-Frhou(i,j,k-1) ) &
                           + a7(2) * ( Frhou(i,j,k+2)-Frhou(i,j,k-2) ) &
                           + a7(3) * ( Frhou(i,j,k+3)-Frhou(i,j,k-3) ) )*idz6_kmin &
                           + Krhou(i,j,k)
             Krhov(i,j,k)= ( a7(1) * ( Frhov(i,j,k+1)-Frhov(i,j,k-1) ) &
                           + a7(2) * ( Frhov(i,j,k+2)-Frhov(i,j,k-2) ) &
                           + a7(3) * ( Frhov(i,j,k+3)-Frhov(i,j,k-3) ) )*idz6_kmin &
                           + Krhov(i,j,k)
             Krhow(i,j,k)= ( a7(1) * ( Frhow(i,j,k+1)-Frhow(i,j,k-1) ) &
                           + a7(2) * ( Frhow(i,j,k+2)-Frhow(i,j,k-2) ) &
                           + a7(3) * ( Frhow(i,j,k+3)-Frhow(i,j,k-3) ) )*idz6_kmin &
                           + Krhow(i,j,k)
             Krhoe(i,j,k)= ( a7(1) * ( Frhoe(i,j,k+1)-Frhoe(i,j,k-1) ) &
                           + a7(2) * ( Frhoe(i,j,k+2)-Frhoe(i,j,k-2) ) &
                           + a7(3) * ( Frhoe(i,j,k+3)-Frhoe(i,j,k-3) ) )*idz6_kmin &
                           + Krhoe(i,j,k)
          enddo
       enddo
    endif

    ! Interior points
    ! ---------------
    do k=ndz,nfz
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Krho(i,j,k) = ( a9(1) * ( Frho(i,j,k+1)-Frho(i,j,k-1) ) &
                           + a9(2) * ( Frho(i,j,k+2)-Frho(i,j,k-2) ) &
                           + a9(3) * ( Frho(i,j,k+3)-Frho(i,j,k-3) ) &
                           + a9(4) * ( Frho(i,j,k+4)-Frho(i,j,k-4) ) )*idz(k) &
                           + Krho(i,j,k)
             Krhou(i,j,k)= ( a9(1) * ( Frhou(i,j,k+1)-Frhou(i,j,k-1) ) &
                           + a9(2) * ( Frhou(i,j,k+2)-Frhou(i,j,k-2) ) &
                           + a9(3) * ( Frhou(i,j,k+3)-Frhou(i,j,k-3) ) &
                           + a9(4) * ( Frhou(i,j,k+4)-Frhou(i,j,k-4) ) )*idz(k) &
                           + Krhou(i,j,k)
             Krhov(i,j,k)= ( a9(1) * ( Frhov(i,j,k+1)-Frhov(i,j,k-1) ) &
                           + a9(2) * ( Frhov(i,j,k+2)-Frhov(i,j,k-2) ) &
                           + a9(3) * ( Frhov(i,j,k+3)-Frhov(i,j,k-3) ) &
                           + a9(4) * ( Frhov(i,j,k+4)-Frhov(i,j,k-4) ) )*idz(k) &
                           + Krhov(i,j,k)
             Krhow(i,j,k)= ( a9(1) * ( Frhow(i,j,k+1)-Frhow(i,j,k-1) ) &
                           + a9(2) * ( Frhow(i,j,k+2)-Frhow(i,j,k-2) ) &
                           + a9(3) * ( Frhow(i,j,k+3)-Frhow(i,j,k-3) ) &
                           + a9(4) * ( Frhow(i,j,k+4)-Frhow(i,j,k-4) ) )*idz(k) &
                           + Krhow(i,j,k)
             Krhoe(i,j,k)= ( a9(1) * ( Frhoe(i,j,k+1)-Frhoe(i,j,k-1) ) &
                           + a9(2) * ( Frhoe(i,j,k+2)-Frhoe(i,j,k-2) ) &
                           + a9(3) * ( Frhoe(i,j,k+3)-Frhoe(i,j,k-3) ) &
                           + a9(4) * ( Frhoe(i,j,k+4)-Frhoe(i,j,k-4) ) )*idz(k) &
                           + Krhoe(i,j,k)
          enddo
       enddo
    enddo

    ! BC at kmax
    ! ----------
    if (is_bc_1pt(3,2)) then
       k=nz-3
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Krho(i,j,k) = ( a7(1) * ( Frho(i,j,k+1)-Frho(i,j,k-1) ) &
                           + a7(2) * ( Frho(i,j,k+2)-Frho(i,j,k-2) ) &
                           + a7(3) * ( Frho(i,j,k+3)-Frho(i,j,k-3) ) )*idz6_kmax &
                           + Krho(i,j,k)
             Krhou(i,j,k)= ( a7(1) * ( Frhou(i,j,k+1)-Frhou(i,j,k-1) ) &
                           + a7(2) * ( Frhou(i,j,k+2)-Frhou(i,j,k-2) ) &
                           + a7(3) * ( Frhou(i,j,k+3)-Frhou(i,j,k-3) ) )*idz6_kmax &
                           + Krhou(i,j,k)
             Krhov(i,j,k)= ( a7(1) * ( Frhov(i,j,k+1)-Frhov(i,j,k-1) ) &
                           + a7(2) * ( Frhov(i,j,k+2)-Frhov(i,j,k-2) ) &
                           + a7(3) * ( Frhov(i,j,k+3)-Frhov(i,j,k-3) ) )*idz6_kmax &
                           + Krhov(i,j,k)
             Krhow(i,j,k)= ( a7(1) * ( Frhow(i,j,k+1)-Frhow(i,j,k-1) ) &
                           + a7(2) * ( Frhow(i,j,k+2)-Frhow(i,j,k-2) ) &
                           + a7(3) * ( Frhow(i,j,k+3)-Frhow(i,j,k-3) ) )*idz6_kmax &
                           + Krhow(i,j,k)
             Krhoe(i,j,k)= ( a7(1) * ( Frhoe(i,j,k+1)-Frhoe(i,j,k-1) ) &
                           + a7(2) * ( Frhoe(i,j,k+2)-Frhoe(i,j,k-2) ) &
                           + a7(3) * ( Frhoe(i,j,k+3)-Frhoe(i,j,k-3) ) )*idz6_kmax &
                           + Krhoe(i,j,k)
          enddo
       enddo

       k=nz-2
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Krho(i,j,k) = ( a5(1)*( Frho(i,j,k+1)- Frho(i,j,k-1)) &
                           + a5(2)*( Frho(i,j,k+2)- Frho(i,j,k-2)) )*idz4_kmax + Krho(i,j,k)
             Krhou(i,j,k)= ( a5(1)*(Frhou(i,j,k+1)-Frhou(i,j,k-1)) &
                           + a5(2)*(Frhou(i,j,k+2)-Frhou(i,j,k-2)) )*idz4_kmax + Krhou(i,j,k)
             Krhov(i,j,k)= ( a5(1)*(Frhov(i,j,k+1)-Frhov(i,j,k-1)) &
                           + a5(2)*(Frhov(i,j,k+2)-Frhov(i,j,k-2)) )*idz4_kmax + Krhov(i,j,k)
             Krhow(i,j,k)= ( a5(1)*(Frhow(i,j,k+1)-Frhow(i,j,k-1)) &
                           + a5(2)*(Frhow(i,j,k+2)-Frhow(i,j,k-2)) )*idz4_kmax + Krhow(i,j,k)
             Krhoe(i,j,k)= ( a5(1)*(Frhoe(i,j,k+1)-Frhoe(i,j,k-1)) &
                           + a5(2)*(Frhoe(i,j,k+2)-Frhoe(i,j,k-2)) )*idz4_kmax + Krhoe(i,j,k)
          enddo
       enddo

       k=nz-1
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Krho(i,j,k)  = Krho(i,j,k)  + ( Frho(i,j,k+1)- Frho(i,j,k-1))*idz2_kmax
             Krhou(i,j,k) = Krhou(i,j,k) + (Frhou(i,j,k+1)-Frhou(i,j,k-1))*idz2_kmax
             Krhov(i,j,k) = Krhov(i,j,k) + (Frhov(i,j,k+1)-Frhov(i,j,k-1))*idz2_kmax
             Krhow(i,j,k) = Krhow(i,j,k) + (Frhow(i,j,k+1)-Frhow(i,j,k-1))*idz2_kmax
             Krhoe(i,j,k) = Krhoe(i,j,k) + (Frhoe(i,j,k+1)-Frhoe(i,j,k-1))*idz2_kmax
          enddo
       enddo
    endif

    ! Advancement of wall points
    ! --------------------------
    if (is_bc_wall2(3,2)) then
       k=nz
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Krho(i,j,k) =(a20(3)*Frho(i,j,k-2) +a20(2)*Frho(i,j,k-1) +a20(1)*Frho(i,j,k) )*idz1_kmax+ Krho(i,j,k)
             Krhou(i,j,k)=(a20(3)*Frhou(i,j,k-2)+a20(2)*Frhou(i,j,k-1)+a20(1)*Frhou(i,j,k))*idz1_kmax+ Krhou(i,j,k)
             Krhov(i,j,k)=(a20(3)*Frhov(i,j,k-2)+a20(2)*Frhov(i,j,k-1)+a20(1)*Frhov(i,j,k))*idz1_kmax+ Krhov(i,j,k)
             Krhow(i,j,k)=(a20(3)*Frhow(i,j,k-2)+a20(2)*Frhow(i,j,k-1)+a20(1)*Frhow(i,j,k))*idz1_kmax+ Krhow(i,j,k)
             Krhoe(i,j,k)=(a20(3)*Frhoe(i,j,k-2)+a20(2)*Frhoe(i,j,k-1)+a20(1)*Frhoe(i,j,k))*idz1_kmax+ Krhoe(i,j,k)
          enddo
       enddo
    endif

  end subroutine flux_euler_9pts_c

end submodule smod_flux_euler_9pts_c
