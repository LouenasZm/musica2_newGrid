!=================================================================================
submodule (mod_gradient) smod_grad_vel_c3
!=================================================================================
  !> author: XG
  !> date: January 2024
  !> Full 3D curvilinear version - routines to compute gradients of velocity components
!=================================================================================

contains

  !===============================================================================
  module subroutine grad_vel_5pts_c3
  !===============================================================================
    !> Compute velocity and temperature derivatives
    !> 5-point stencil - 3D curvilinear version -
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: dvardx,dvardy,dvardz
    ! ---------------------------------------------------------------------------

    ! Derivatives along ksi (Caution: dux is du/dksi, etc ...)
    ! =====================
    if (is_boundary(1,1)) then
       i=1
       do k=1,nz
          do j=1,ny
             dux(i,j,k)= a04(1)*uu(i  ,j,k)+a04(2)*uu(i+1,j,k) &
                        +a04(3)*uu(i+2,j,k)+a04(4)*uu(i+3,j,k) &
                        +a04(5)*uu(i+4,j,k)
             dvx(i,j,k)= a04(1)*vv(i  ,j,k)+a04(2)*vv(i+1,j,k) &
                        +a04(3)*vv(i+2,j,k)+a04(4)*vv(i+3,j,k) &
                        +a04(5)*vv(i+4,j,k)
             dwx(i,j,k)= a04(1)*ww(i  ,j,k)+a04(2)*ww(i+1,j,k) &
                        +a04(3)*ww(i+2,j,k)+a04(4)*ww(i+3,j,k) &
                        +a04(5)*ww(i+4,j,k)
          enddo
       enddo

       i=2
       do k=1,nz
          do j=1,ny
             dux(i,j,k)= a13(1)*uu(i-1,j,k)+a13(2)*uu(i  ,j,k) &
                        +a13(3)*uu(i+1,j,k)+a13(4)*uu(i+2,j,k) &
                        +a13(5)*uu(i+3,j,k)
             dvx(i,j,k)= a13(1)*vv(i-1,j,k)+a13(2)*vv(i  ,j,k) &
                        +a13(3)*vv(i+1,j,k)+a13(4)*vv(i+2,j,k) &
                        +a13(5)*vv(i+3,j,k)
             dwx(i,j,k)= a13(1)*ww(i-1,j,k)+a13(2)*ww(i  ,j,k) &
                        +a13(3)*ww(i+1,j,k)+a13(4)*ww(i+2,j,k) &
                        +a13(5)*ww(i+3,j,k)
          enddo
       enddo
    endif

    do k=1,nz
       do j=1,ny
          do i=ndx_v3,nfx_v3
             dux(i,j,k) = a5(1)*( uu(i+1,j,k)- uu(i-1,j,k)) &
                        + a5(2)*( uu(i+2,j,k)- uu(i-2,j,k))
             dvx(i,j,k) = a5(1)*( vv(i+1,j,k)- vv(i-1,j,k)) &
                        + a5(2)*( vv(i+2,j,k)- vv(i-2,j,k))
             dwx(i,j,k) = a5(1)*( ww(i+1,j,k)- ww(i-1,j,k)) &
                        + a5(2)*( ww(i+2,j,k)- ww(i-2,j,k))
          enddo
       enddo
    enddo

    if (is_boundary(1,2)) then
       i=nx-1
       do k=1,nz
          do j=1,ny
             dux(i,j,k)= a31(5)*uu(i-3,j,k)+a31(4)*uu(i-2,j,k) &
                        +a31(3)*uu(i-1,j,k)+a31(2)*uu(i  ,j,k) &
                        +a31(1)*uu(i+1,j,k)
             dvx(i,j,k)= a31(5)*vv(i-3,j,k)+a31(4)*vv(i-2,j,k) &
                        +a31(3)*vv(i-1,j,k)+a31(2)*vv(i  ,j,k) &
                        +a31(1)*vv(i+1,j,k)
             dwx(i,j,k)= a31(5)*ww(i-3,j,k)+a31(4)*ww(i-2,j,k) &
                        +a31(3)*ww(i-1,j,k)+a31(2)*ww(i  ,j,k) &
                        +a31(1)*ww(i+1,j,k)
          enddo
       enddo

       i=nx
       do k=1,nz
          do j=1,ny
             dux(i,j,k)= a40(5)*uu(i-4,j,k)+a40(4)*uu(i-3,j,k) &
                        +a40(3)*uu(i-2,j,k)+a40(2)*uu(i-1,j,k) &
                        +a40(1)*uu(i  ,j,k)
             dvx(i,j,k)= a40(5)*vv(i-4,j,k)+a40(4)*vv(i-3,j,k) &
                        +a40(3)*vv(i-2,j,k)+a40(2)*vv(i-1,j,k) &
                        +a40(1)*vv(i  ,j,k)
             dwx(i,j,k)= a40(5)*ww(i-4,j,k)+a40(4)*ww(i-3,j,k) &
                        +a40(3)*ww(i-2,j,k)+a40(2)*ww(i-1,j,k) &
                        +a40(1)*ww(i  ,j,k)
          enddo
       enddo
    endif

    ! Derivatives along eta (Caution: duy is du/deta, etc ...)
    ! =====================
    if (is_boundary(2,1)) then
       j=1
       do k=1,nz
          do i=1,nx
             duy(i,j,k)= a04(1)*uu(i,j  ,k)+a04(2)*uu(i,j+1,k) &
                        +a04(3)*uu(i,j+2,k)+a04(4)*uu(i,j+3,k) &
                        +a04(5)*uu(i,j+4,k)
             dvy(i,j,k)= a04(1)*vv(i,j  ,k)+a04(2)*vv(i,j+1,k) &
                        +a04(3)*vv(i,j+2,k)+a04(4)*vv(i,j+3,k) &
                        +a04(5)*vv(i,j+4,k)
             dwy(i,j,k)= a04(1)*ww(i,j  ,k)+a04(2)*ww(i,j+1,k) &
                        +a04(3)*ww(i,j+2,k)+a04(4)*ww(i,j+3,k) &
                        +a04(5)*ww(i,j+4,k)
          enddo
       enddo

       j=2
       do k=1,nz
          do i=1,nx
             duy(i,j,k)= a13(1)*uu(i,j-1,k)+a13(2)*uu(i,j  ,k) &
                        +a13(3)*uu(i,j+1,k)+a13(4)*uu(i,j+2,k) &
                        +a13(5)*uu(i,j+3,k)
             dvy(i,j,k)= a13(1)*vv(i,j-1,k)+a13(2)*vv(i,j  ,k) &
                        +a13(3)*vv(i,j+1,k)+a13(4)*vv(i,j+2,k) &
                        +a13(5)*vv(i,j+3,k)
             dwy(i,j,k)= a13(1)*ww(i,j-1,k)+a13(2)*ww(i,j  ,k) &
                        +a13(3)*ww(i,j+1,k)+a13(4)*ww(i,j+2,k) &
                        +a13(5)*ww(i,j+3,k)
          enddo
       enddo
    endif

    do k=1,nz
       do j=ndy_v3,nfy_v3
          do i=1,nx
             duy(i,j,k) = a5(1)*( uu(i,j+1,k)- uu(i,j-1,k)) &
                        + a5(2)*( uu(i,j+2,k)- uu(i,j-2,k))
             dvy(i,j,k) = a5(1)*( vv(i,j+1,k)- vv(i,j-1,k)) &
                        + a5(2)*( vv(i,j+2,k)- vv(i,j-2,k))
             dwy(i,j,k) = a5(1)*( ww(i,j+1,k)- ww(i,j-1,k)) &
                        + a5(2)*( ww(i,j+2,k)- ww(i,j-2,k))
          enddo
       enddo
    enddo

    if (is_boundary(2,2)) then
       j=ny-1
       do k=1,nz
          do i=1,nx
             duy(i,j,k)= a31(5)*uu(i,j-3,k)+a31(4)*uu(i,j-2,k) &
                        +a31(3)*uu(i,j-1,k)+a31(2)*uu(i,j  ,k) &
                        +a31(1)*uu(i,j+1,k)
             dvy(i,j,k)= a31(5)*vv(i,j-3,k)+a31(4)*vv(i,j-2,k) &
                        +a31(3)*vv(i,j-1,k)+a31(2)*vv(i,j  ,k) &
                        +a31(1)*vv(i,j+1,k)
             dwy(i,j,k)= a31(5)*ww(i,j-3,k)+a31(4)*ww(i,j-2,k) &
                        +a31(3)*ww(i,j-1,k)+a31(2)*ww(i,j  ,k) &
                        +a31(1)*ww(i,j+1,k)
          enddo
       enddo

       j=ny
       do k=1,nz
          do i=1,nx
             duy(i,j,k)= a40(5)*uu(i,j-4,k)+a40(4)*uu(i,j-3,k) &
                        +a40(3)*uu(i,j-2,k)+a40(2)*uu(i,j-1,k) &
                        +a40(1)*uu(i,j  ,k)
             dvy(i,j,k)= a40(5)*vv(i,j-4,k)+a40(4)*vv(i,j-3,k) &
                        +a40(3)*vv(i,j-2,k)+a40(2)*vv(i,j-1,k) &
                        +a40(1)*vv(i,j  ,k)
             dwy(i,j,k)= a40(5)*ww(i,j-4,k)+a40(4)*ww(i,j-3,k) &
                        +a40(3)*ww(i,j-2,k)+a40(2)*ww(i,j-1,k) &
                        +a40(1)*ww(i,j  ,k)
          enddo
       enddo
    endif

    ! Derivatives along phi (Caution: duz is du/dphi, etc ...)
    ! =====================
    if (is_boundary(3,1)) then
       k=1
       do j=1,ny
          do i=1,nx
             duz(i,j,k)= a04(1)*uu(i,j,k  )+a04(2)*uu(i,j,k+1) &
                        +a04(3)*uu(i,j,k+2)+a04(4)*uu(i,j,k+3) &
                        +a04(5)*uu(i,j,k+4)
             dvz(i,j,k)= a04(1)*vv(i,j,k  )+a04(2)*vv(i,j,k+1) &
                        +a04(3)*vv(i,j,k+2)+a04(4)*vv(i,j,k+3) &
                        +a04(5)*vv(i,j,k+4)
             dwz(i,j,k)= a04(1)*ww(i,j,k  )+a04(2)*ww(i,j,k+1) &
                        +a04(3)*ww(i,j,k+2)+a04(4)*ww(i,j,k+3) &
                        +a04(5)*ww(i,j,k+4)
          enddo
       enddo

       k=2
       do j=1,ny
          do i=1,nx
             duz(i,j,k)= a13(1)*uu(i,j,k-1)+a13(2)*uu(i,j,k  ) &
                        +a13(3)*uu(i,j,k+1)+a13(4)*uu(i,j,k+2) &
                        +a13(5)*uu(i,j,k+3)
             dvz(i,j,k)= a13(1)*vv(i,j,k-1)+a13(2)*vv(i,j,k  ) &
                        +a13(3)*vv(i,j,k+1)+a13(4)*vv(i,j,k+2) &
                        +a13(5)*vv(i,j,k+3)
             dwz(i,j,k)= a13(1)*ww(i,j,k-1)+a13(2)*ww(i,j,k  ) &
                        +a13(3)*ww(i,j,k+1)+a13(4)*ww(i,j,k+2) &
                        +a13(5)*ww(i,j,k+3)
          enddo
       enddo
    endif

    do k=ndz_v3,nfz_v3
       do j=1,ny
          do i=1,nx
             duz(i,j,k)= a5(1)*( uu(i,j,k+1)- uu(i,j,k-1)) &
                       + a5(2)*( uu(i,j,k+2)- uu(i,j,k-2))
             dvz(i,j,k)= a5(1)*( vv(i,j,k+1)- vv(i,j,k-1)) &
                       + a5(2)*( vv(i,j,k+2)- vv(i,j,k-2))
             dwz(i,j,k)= a5(1)*( ww(i,j,k+1)- ww(i,j,k-1)) &
                       + a5(2)*( ww(i,j,k+2)- ww(i,j,k-2))
          enddo
       enddo
    enddo

    if (is_boundary(3,2)) then
       k=nz-1
       do j=1,ny
          do i=1,nx
             duz(i,j,k)= a31(5)*uu(i,j,k-3)+a31(4)*uu(i,j,k-2) &
                        +a31(3)*uu(i,j,k-1)+a31(2)*uu(i,j,k  ) &
                        +a31(1)*uu(i,j,k+1)
             dvz(i,j,k)= a31(5)*vv(i,j,k-3)+a31(4)*vv(i,j,k-2) &
                        +a31(3)*vv(i,j,k-1)+a31(2)*vv(i,j,k  ) &
                        +a31(1)*vv(i,j,k+1)
             dwz(i,j,k)= a31(5)*ww(i,j,k-3)+a31(4)*ww(i,j,k-2) &
                        +a31(3)*ww(i,j,k-1)+a31(2)*ww(i,j,k  ) &
                        +a31(1)*ww(i,j,k+1)
          enddo
       enddo

       k=nz
       do j=1,ny
          do i=1,nx
             duz(i,j,k)= a40(5)*uu(i,j,k-4)+a40(4)*uu(i,j,k-3) &
                        +a40(3)*uu(i,j,k-2)+a40(2)*uu(i,j,k-1) &
                        +a40(1)*uu(i,j,k  )
             dvz(i,j,k)= a40(5)*vv(i,j,k-4)+a40(4)*vv(i,j,k-3) &
                        +a40(3)*vv(i,j,k-2)+a40(2)*vv(i,j,k-1) &
                        +a40(1)*vv(i,j,k  )
             dwz(i,j,k)= a40(5)*ww(i,j,k-4)+a40(4)*ww(i,j,k-3) &
                        +a40(3)*ww(i,j,k-2)+a40(2)*ww(i,j,k-1) &
                        +a40(1)*ww(i,j,k  )
          enddo
       enddo
    endif

    ! Compute derivatives in Cartesian coordinates
    ! ============================================
    do k=1,nz
       do j=1,ny
          do i=1,nx
             ! u velocity (x-component)
             ! ------------------------
             dvardx=dux(i,j,k)*ksi_x_v(i,j,k)+duy(i,j,k)*eta_x_v(i,j,k)+duz(i,j,k)*phi_x_v(i,j,k)
             dvardy=dux(i,j,k)*ksi_y_v(i,j,k)+duy(i,j,k)*eta_y_v(i,j,k)+duz(i,j,k)*phi_y_v(i,j,k)
             dvardz=dux(i,j,k)*ksi_z_v(i,j,k)+duy(i,j,k)*eta_z_v(i,j,k)+duz(i,j,k)*phi_z_v(i,j,k)
             ! du/dx
             dux(i,j,k)=dvardx*ijacob3_v(i,j,k)
             ! du/dy
             duy(i,j,k)=dvardy*ijacob3_v(i,j,k)
             ! du/dz
             duz(i,j,k)=dvardz*ijacob3_v(i,j,k)

             ! v velocity (y-component)
             ! ------------------------
             dvardx=dvx(i,j,k)*ksi_x_v(i,j,k)+dvy(i,j,k)*eta_x_v(i,j,k)+dvz(i,j,k)*phi_x_v(i,j,k)
             dvardy=dvx(i,j,k)*ksi_y_v(i,j,k)+dvy(i,j,k)*eta_y_v(i,j,k)+dvz(i,j,k)*phi_y_v(i,j,k)
             dvardz=dvx(i,j,k)*ksi_z_v(i,j,k)+dvy(i,j,k)*eta_z_v(i,j,k)+dvz(i,j,k)*phi_z_v(i,j,k)
             ! dv/dx
             dvx(i,j,k)=dvardx*ijacob3_v(i,j,k)
             ! dv/dy
             dvy(i,j,k)=dvardy*ijacob3_v(i,j,k)
             ! dv/dz
             dvz(i,j,k)=dvardz*ijacob3_v(i,j,k)

             ! w velocity (z-component)
             ! ------------------------
             dvardx=dwx(i,j,k)*ksi_x_v(i,j,k)+dwy(i,j,k)*eta_x_v(i,j,k)+dwz(i,j,k)*phi_x_v(i,j,k)
             dvardy=dwx(i,j,k)*ksi_y_v(i,j,k)+dwy(i,j,k)*eta_y_v(i,j,k)+dwz(i,j,k)*phi_y_v(i,j,k)
             dvardz=dwx(i,j,k)*ksi_z_v(i,j,k)+dwy(i,j,k)*eta_z_v(i,j,k)+dwz(i,j,k)*phi_z_v(i,j,k)
             ! dw/dx
             dwx(i,j,k)=dvardx*ijacob3_v(i,j,k)
             ! dw/dy
             dwy(i,j,k)=dvardy*ijacob3_v(i,j,k)
             ! dw/dz
             dwz(i,j,k)=dvardz*ijacob3_v(i,j,k)
          enddo
       enddo
    enddo

  end subroutine grad_vel_5pts_c3

  !===============================================================================
  module subroutine grad_vel_3pts_c3
  !===============================================================================
    !> Compute velocity and temperature derivatives
    !> 3-point stencil - 3D curvilinear version -
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: dvardx,dvardy,dvardz
    ! ---------------------------------------------------------------------------

    ! Derivatives along ksi (Caution: dux is du/dksi, etc ...)
    ! =====================
    if (is_boundary(1,1)) then
       i=1
       do k=1,nz
          do j=1,ny
             dux(i,j,k)= a02(1)*uu(i  ,j,k)+a02(2)*uu(i+1,j,k) &
                        +a02(3)*uu(i+2,j,k)
             dvx(i,j,k)= a02(1)*vv(i  ,j,k)+a02(2)*vv(i+1,j,k) &
                        +a02(3)*vv(i+2,j,k)
             dwx(i,j,k)= a02(1)*ww(i  ,j,k)+a02(2)*ww(i+1,j,k) &
                        +a02(3)*ww(i+2,j,k)
          enddo
       enddo
    endif

    do k=1,nz
       do j=1,ny
          do i=ndx_v3,nfx_v3
             dux(i,j,k) = 0.5_wp*(uu(i+1,j,k)- uu(i-1,j,k))
             dvx(i,j,k) = 0.5_wp*(vv(i+1,j,k)- vv(i-1,j,k))
             dwx(i,j,k) = 0.5_wp*(ww(i+1,j,k)- ww(i-1,j,k))
          enddo
       enddo
    enddo

    if (is_boundary(1,2)) then
       i=nx
       do k=1,nz
          do j=1,ny
             dux(i,j,k)= a20(3)*uu(i-2,j,k)+a20(2)*uu(i-1,j,k) &
                        +a20(1)*uu(i  ,j,k)
             dvx(i,j,k)= a20(3)*vv(i-2,j,k)+a20(2)*vv(i-1,j,k) &
                        +a20(1)*vv(i  ,j,k)
             dwx(i,j,k)= a20(3)*ww(i-2,j,k)+a20(2)*ww(i-1,j,k) &
                        +a20(1)*ww(i  ,j,k)
          enddo
       enddo
    endif

    ! Derivatives along eta (Caution: duy is du/deta, etc ...)
    ! =====================
    if (is_boundary(2,1)) then
       j=1
       do k=1,nz
          do i=1,nx
             duy(i,j,k)= a02(1)*uu(i,j  ,k)+a02(2)*uu(i,j+1,k) &
                        +a02(3)*uu(i,j+2,k)
             dvy(i,j,k)= a02(1)*vv(i,j  ,k)+a02(2)*vv(i,j+1,k) &
                        +a02(3)*vv(i,j+2,k)
             dwy(i,j,k)= a02(1)*ww(i,j  ,k)+a02(2)*ww(i,j+1,k) &
                        +a02(3)*ww(i,j+2,k)
          enddo
       enddo
    endif

    do k=1,nz
       do j=ndy_v3,nfy_v3
          do i=1,nx
             duy(i,j,k) = 0.5_wp*(uu(i,j+1,k)- uu(i,j-1,k))
             dvy(i,j,k) = 0.5_wp*(vv(i,j+1,k)- vv(i,j-1,k))
             dwy(i,j,k) = 0.5_wp*(ww(i,j+1,k)- ww(i,j-1,k))
          enddo
       enddo
    enddo

    if (is_boundary(2,2)) then
       j=ny
       do k=1,nz
          do i=1,nx
             duy(i,j,k)= a20(3)*uu(i,j-2,k)+a20(2)*uu(i,j-1,k) &
                        +a20(1)*uu(i,j  ,k)
             dvy(i,j,k)= a20(3)*vv(i,j-2,k)+a20(2)*vv(i,j-1,k) &
                        +a20(1)*vv(i,j  ,k)
             dwy(i,j,k)= a20(3)*ww(i,j-2,k)+a20(2)*ww(i,j-1,k) &
                        +a20(1)*ww(i,j  ,k)
          enddo
       enddo
    endif

    ! Derivatives along phi (Caution: duz is du/dphi, etc ...)
    ! =====================
    if (is_boundary(3,1)) then
       k=1
       do j=1,ny
          do i=1,nx
             duz(i,j,k)= a02(1)*uu(i,j,k  )+a02(2)*uu(i,j,k+1) &
                        +a02(3)*uu(i,j,k+2)
             dvz(i,j,k)= a02(1)*vv(i,j,k  )+a02(2)*vv(i,j,k+1) &
                        +a02(3)*vv(i,j,k+2)
             dwz(i,j,k)= a02(1)*ww(i,j,k  )+a02(2)*ww(i,j,k+1) &
                        +a02(3)*ww(i,j,k+2)
          enddo
       enddo
    endif

    do k=ndz_v3,nfz_v3
       do j=1,ny
          do i=1,nx
             duz(i,j,k)= 0.5_wp*(uu(i,j,k+1)- uu(i,j,k-1))
             dvz(i,j,k)= 0.5_wp*(vv(i,j,k+1)- vv(i,j,k-1))
             dwz(i,j,k)= 0.5_wp*(ww(i,j,k+1)- ww(i,j,k-1))
          enddo
       enddo
    enddo

    if (is_boundary(3,2)) then
       k=nz
       do j=1,ny
          do i=1,nx
             duz(i,j,k)= a20(3)*uu(i,j,k-2)+a20(2)*uu(i,j,k-1) &
                        +a20(1)*uu(i,j,k  )
             dvz(i,j,k)= a20(3)*vv(i,j,k-2)+a20(2)*vv(i,j,k-1) &
                        +a20(1)*vv(i,j,k  )
             dwz(i,j,k)= a20(3)*ww(i,j,k-2)+a20(2)*ww(i,j,k-1) &
                        +a20(1)*ww(i,j,k  )
          enddo
       enddo
    endif

    ! Compute derivatives in Cartesian coordinates
    ! ============================================
    do k=1,nz
       do j=1,ny
          do i=1,nx
             ! u velocity (x-component)
             ! ------------------------
             dvardx=dux(i,j,k)*ksi_x_v(i,j,k)+duy(i,j,k)*eta_x_v(i,j,k)+duz(i,j,k)*phi_x_v(i,j,k)
             dvardy=dux(i,j,k)*ksi_y_v(i,j,k)+duy(i,j,k)*eta_y_v(i,j,k)+duz(i,j,k)*phi_y_v(i,j,k)
             dvardz=dux(i,j,k)*ksi_z_v(i,j,k)+duy(i,j,k)*eta_z_v(i,j,k)+duz(i,j,k)*phi_z_v(i,j,k)
             ! du/dx
             dux(i,j,k)=dvardx*ijacob3_v(i,j,k)
             ! du/dy
             duy(i,j,k)=dvardy*ijacob3_v(i,j,k)
             ! du/dz
             duz(i,j,k)=dvardz*ijacob3_v(i,j,k)

             ! v velocity (y-component)
             ! ------------------------
             dvardx=dvx(i,j,k)*ksi_x_v(i,j,k)+dvy(i,j,k)*eta_x_v(i,j,k)+dvz(i,j,k)*phi_x_v(i,j,k)
             dvardy=dvx(i,j,k)*ksi_y_v(i,j,k)+dvy(i,j,k)*eta_y_v(i,j,k)+dvz(i,j,k)*phi_y_v(i,j,k)
             dvardz=dvx(i,j,k)*ksi_z_v(i,j,k)+dvy(i,j,k)*eta_z_v(i,j,k)+dvz(i,j,k)*phi_z_v(i,j,k)
             ! dv/dx
             dvx(i,j,k)=dvardx*ijacob3_v(i,j,k)
             ! dv/dy
             dvy(i,j,k)=dvardy*ijacob3_v(i,j,k)
             ! dv/dz
             dvz(i,j,k)=dvardz*ijacob3_v(i,j,k)

             ! w velocity (z-component)
             ! ------------------------
             dvardx=dwx(i,j,k)*ksi_x_v(i,j,k)+dwy(i,j,k)*eta_x_v(i,j,k)+dwz(i,j,k)*phi_x_v(i,j,k)
             dvardy=dwx(i,j,k)*ksi_y_v(i,j,k)+dwy(i,j,k)*eta_y_v(i,j,k)+dwz(i,j,k)*phi_y_v(i,j,k)
             dvardz=dwx(i,j,k)*ksi_z_v(i,j,k)+dwy(i,j,k)*eta_z_v(i,j,k)+dwz(i,j,k)*phi_z_v(i,j,k)
             ! dw/dx
             dwx(i,j,k)=dvardx*ijacob3_v(i,j,k)
             ! dw/dy
             dwy(i,j,k)=dvardy*ijacob3_v(i,j,k)
             ! dw/dz
             dwz(i,j,k)=dvardz*ijacob3_v(i,j,k)
          enddo
       enddo
    enddo

  end subroutine grad_vel_3pts_c3

end submodule smod_grad_vel_c3
