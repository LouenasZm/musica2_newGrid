!=================================================================================
submodule (mod_flux_euler) smod_flux_euler_5pts_SBP4_c
!=================================================================================
  !> author: XG
  !> date: January 2024
  !> 2.5D curvilinear version - compute Eulerian fluxes (inviscid part) with 5-pt stencil
  !> using SBP4 boundary schemes
!=================================================================================

contains

  !==============================================================================
  module subroutine flux_euler_5pts_SBP4_c
  !==============================================================================
    !> Derivatives of Eulerian fluxes (inviscid part) - 5-point stencil -
    !> - 2.5D curvilinear version with SBP4 boundary schemes -
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: vc,rovc
    ! ---------------------------------------------------------------------------

    ! Initialize flux derivative arrays at wall points
    ! ================================================
    ! Wall BC at imin
    ! ---------------
    !  if (is_bc_wall(1,1)) then
    !      Krho(1,ndy_e:nfy_e,ndz_e:nfz_e)=0.0_wp
    !     Krhou(1,ndy_e:nfy_e,ndz_e:nfz_e)=0.0_wp
    !     Krhov(1,ndy_e:nfy_e,ndz_e:nfz_e)=0.0_wp
    !     Krhow(1,ndy_e:nfy_e,ndz_e:nfz_e)=0.0_wp
    !     Krhoe(1,ndy_e:nfy_e,ndz_e:nfz_e)=0.0_wp
    !  endif
    !  ! Wall BC at imax
    !  ! ---------------
    !  if (is_bc_wall(1,2)) then
    !      Krho(nx,ndy_e:nfy_e,ndz_e:nfz_e)=0.0_wp
    !     Krhou(nx,ndy_e:nfy_e,ndz_e:nfz_e)=0.0_wp
    !     Krhov(nx,ndy_e:nfy_e,ndz_e:nfz_e)=0.0_wp
    !     Krhow(nx,ndy_e:nfy_e,ndz_e:nfz_e)=0.0_wp
    !     Krhoe(nx,ndy_e:nfy_e,ndz_e:nfz_e)=0.0_wp
    !  endif
!!$    ! Wall BC at jmin
!!$    ! ---------------
!!$    if (is_bc_wall(2,1)) then
!!$       Krho(ndx_e:nfx_e,1,ndz_e:nfz_e)=0.0_wp
!!$       Krhou(ndx_e:nfx_e,1,ndz_e:nfz_e)=0.0_wp
!!$       Krhov(ndx_e:nfx_e,1,ndz_e:nfz_e)=0.0_wp
!!$       Krhow(ndx_e:nfx_e,1,ndz_e:nfz_e)=0.0_wp
!!$       Krhoe(ndx_e:nfx_e,1,ndz_e:nfz_e)=0.0_wp
!!$    endif
!!$    ! Wall BC at jmax
!!$    ! ---------------
!!$    if (is_bc_wall(2,2)) then
!!$       Krho(ndx_e:nfx_e,ny,ndz_e:nfz_e)=0.0_wp
!!$       Krhou(ndx_e:nfx_e,ny,ndz_e:nfz_e)=0.0_wp
!!$       Krhov(ndx_e:nfx_e,ny,ndz_e:nfz_e)=0.0_wp
!!$       Krhow(ndx_e:nfx_e,ny,ndz_e:nfz_e)=0.0_wp
!!$       Krhoe(ndx_e:nfx_e,ny,ndz_e:nfz_e)=0.0_wp
!!$    endif
    ! Wall BC at kmin
    ! ---------------
    !  if (is_bc_wall(3,1)) then
    !      Krho(ndx_e:nfx_e,ndy_e:nfy_e,1)=0.0_wp
    !     Krhou(ndx_e:nfx_e,ndy_e:nfy_e,1)=0.0_wp
    !     Krhov(ndx_e:nfx_e,ndy_e:nfy_e,1)=0.0_wp
    !     Krhow(ndx_e:nfx_e,ndy_e:nfy_e,1)=0.0_wp
    !     Krhoe(ndx_e:nfx_e,ndy_e:nfy_e,1)=0.0_wp
    !  endif
    !  ! Wall BC at kmax
    !  ! ---------------
    !  if (is_bc_wall(3,2)) then
    !      Krho(ndx_e:nfx_e,ndy_e:nfy_e,nz)=0.0_wp
    !     Krhou(ndx_e:nfx_e,ndy_e:nfy_e,nz)=0.0_wp
    !     Krhov(ndx_e:nfx_e,ndy_e:nfy_e,nz)=0.0_wp
    !     Krhow(ndx_e:nfx_e,ndy_e:nfy_e,nz)=0.0_wp
    !     Krhoe(ndx_e:nfx_e,ndy_e:nfy_e,nz)=0.0_wp
    !  endif

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
    if (is_bc_1pt(1,1)) then
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
    endif

    ! Interior points
    ! ---------------
    do k=ndz_e,nfz_e
       do j=ndy_e,nfy_e
          do i=ndx,nfx
             Krho(i,j,k) = a5(1)*( Frho(i+1,j,k)- Frho(i-1,j,k)) &
                         + a5(2)*( Frho(i+2,j,k)- Frho(i-2,j,k))
             Krhou(i,j,k)= a5(1)*(Frhou(i+1,j,k)-Frhou(i-1,j,k)) &
                         + a5(2)*(Frhou(i+2,j,k)-Frhou(i-2,j,k))
             Krhov(i,j,k)= a5(1)*(Frhov(i+1,j,k)-Frhov(i-1,j,k)) &
                         + a5(2)*(Frhov(i+2,j,k)-Frhov(i-2,j,k))
             Krhow(i,j,k)= a5(1)*(Frhow(i+1,j,k)-Frhow(i-1,j,k)) &
                         + a5(2)*(Frhow(i+2,j,k)-Frhow(i-2,j,k))
             Krhoe(i,j,k)= a5(1)*(Frhoe(i+1,j,k)-Frhoe(i-1,j,k)) &
                         + a5(2)*(Frhoe(i+2,j,k)-Frhoe(i-2,j,k))
          enddo
       enddo
    enddo

    ! Wall BC at imax
    ! ---------------
    if (is_bc_1pt(1,2)) then
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

       !i=nx
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
    if (is_bc_1pt(2,1)) then
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
    endif

    ! Interior points
    ! ---------------
    do k=ndz_e,nfz_e
       do j=ndy,nfy
          do i=ndx_e,nfx_e
             Krho(i,j,k)=  a5(1)*(Frho(i,j+1,k)-Frho(i,j-1,k)) &
                         + a5(2)*(Frho(i,j+2,k)-Frho(i,j-2,k)) + Krho(i,j,k)
             Krhou(i,j,k)= a5(1)*(Frhou(i,j+1,k)-Frhou(i,j-1,k)) &
                         + a5(2)*(Frhou(i,j+2,k)-Frhou(i,j-2,k)) + Krhou(i,j,k)
             Krhov(i,j,k)= a5(1)*(Frhov(i,j+1,k)-Frhov(i,j-1,k)) &
                         + a5(2)*(Frhov(i,j+2,k)-Frhov(i,j-2,k)) + Krhov(i,j,k)
             Krhow(i,j,k)= a5(1)*(Frhow(i,j+1,k)-Frhow(i,j-1,k)) &
                         + a5(2)*(Frhow(i,j+2,k)-Frhow(i,j-2,k)) + Krhow(i,j,k)
             Krhoe(i,j,k)= a5(1)*(Frhoe(i,j+1,k)-Frhoe(i,j-1,k)) &
                         + a5(2)*(Frhoe(i,j+2,k)-Frhoe(i,j-2,k)) + Krhoe(i,j,k)
          enddo
       enddo
    enddo

    ! Wall BC at jmax
    ! ---------------
    if (is_bc_1pt(2,2)) then
       j=ny-1
       do k=ndz_e,nfz_e
          do i=ndx_e,nfx_e
             Krho(i,j,k) = as4m1(1)* Frho(i,ny  ,k)+as4m1(2)* Frho(i,ny-1,k) &
                         + as4m1(3)* Frho(i,ny-2,k)+as4m1(4)* Frho(i,ny-3,k) + Krho(i,j,k)
             Krhou(i,j,k)= as4m1(1)*Frhou(i,ny  ,k)+as4m1(2)*Frhou(i,ny-1,k) &
                         + as4m1(3)*Frhou(i,ny-2,k)+as4m1(4)*Frhou(i,ny-3,k) + Krhou(i,j,k)
             Krhov(i,j,k)= as4m1(1)*Frhov(i,ny  ,k)+as4m1(2)*Frhov(i,ny-1,k) &
                         + as4m1(3)*Frhov(i,ny-2,k)+as4m1(4)*Frhov(i,ny-3,k) + Krhov(i,j,k)
             Krhov(i,j,k)= as4m1(1)*Frhov(i,ny  ,k)+as4m1(2)*Frhov(i,ny-1,k) &
                         + as4m1(3)*Frhov(i,ny-2,k)+as4m1(4)*Frhov(i,ny-3,k) + Krhov(i,j,k)
             Krhoe(i,j,k)= as4m1(1)*Frhoe(i,ny  ,k)+as4m1(2)*Frhoe(i,ny-1,k) &
                         + as4m1(3)*Frhoe(i,ny-2,k)+as4m1(4)*Frhoe(i,ny-3,k) + Krhoe(i,j,k)
          enddo
       enddo

       j=ny
       do k=ndz_e,nfz_e
          do i=ndx_e,nfx_e
             Krho(i,j,k) = as4m0(1)* Frho(i,ny  ,k)+as4m0(2)* Frho(i,ny-1,k) &
                         + as4m0(3)* Frho(i,ny-2,k)+as4m0(4)* Frho(i,ny-3,k) + Krho(i,j,k)
             Krhou(i,j,k)= as4m0(1)*Frhou(i,ny  ,k)+as4m0(2)*Frhou(i,ny-1,k) &
                         + as4m0(3)*Frhou(i,ny-2,k)+as4m0(4)*Frhou(i,ny-3,k) + Krhou(i,j,k)
             Krhov(i,j,k)= as4m0(1)*Frhov(i,ny  ,k)+as4m0(2)*Frhov(i,ny-1,k) &
                         + as4m0(3)*Frhov(i,ny-2,k)+as4m0(4)*Frhov(i,ny-3,k) + Krhov(i,j,k)
             Krhov(i,j,k)= as4m0(1)*Frhov(i,ny  ,k)+as4m0(2)*Frhov(i,ny-1,k) &
                         + as4m0(3)*Frhov(i,ny-2,k)+as4m0(4)*Frhov(i,ny-3,k) + Krhov(i,j,k)
             Krhoe(i,j,k)= as4m0(1)*Frhoe(i,ny  ,k)+as4m0(2)*Frhoe(i,ny-1,k) &
                         + as4m0(3)*Frhoe(i,ny-2,k)+as4m0(4)*Frhoe(i,ny-3,k) + Krhoe(i,j,k)
          enddo
       enddo
    endif

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
    endif

    ! Interior points
    ! ---------------
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

    ! BC at kmax
    ! ----------
    if (is_bc_1pt(3,2)) then
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

  end subroutine flux_euler_5pts_SBP4_c

end submodule smod_flux_euler_5pts_SBP4_c
