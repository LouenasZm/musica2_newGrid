!=================================================================================
submodule (mod_gradient) smod_grad_vel
!=================================================================================
  !> author: XG
  !> date: January 2024
  !> Cartesian version - routines to compute gradients of velocity components
!=================================================================================

contains

  !===============================================================================
  module subroutine grad_vel_5pts
  !===============================================================================
    !> Compute velocity derivatives
    !> 5-point stencil - Cartesian version -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k
    ! ----------------------------------------------------------------------------

    ! Derivatives along x
    ! ===================
    if (is_boundary(1,1)) then
       i=1
       do k=1,nz
          do j=1,ny
             dux(i,j,k)= ( a04(1)*uu(i  ,j,k)+a04(2)*uu(i+1,j,k) &
                          +a04(3)*uu(i+2,j,k)+a04(4)*uu(i+3,j,k) &
                          +a04(5)*uu(i+4,j,k) )*idx_v(i)
             dvx(i,j,k)= ( a04(1)*vv(i  ,j,k)+a04(2)*vv(i+1,j,k) &
                          +a04(3)*vv(i+2,j,k)+a04(4)*vv(i+3,j,k) &
                          +a04(5)*vv(i+4,j,k) )*idx_v(i)
             dwx(i,j,k)= ( a04(1)*ww(i  ,j,k)+a04(2)*ww(i+1,j,k) &
                          +a04(3)*ww(i+2,j,k)+a04(4)*ww(i+3,j,k) &
                          +a04(5)*ww(i+4,j,k) )*idx_v(i)
          enddo
       enddo

       i=2
       do k=1,nz
          do j=1,ny
             dux(i,j,k)= ( a13(1)*uu(i-1,j,k)+a13(2)*uu(i  ,j,k) &
                          +a13(3)*uu(i+1,j,k)+a13(4)*uu(i+2,j,k) &
                          +a13(5)*uu(i+3,j,k) )*idx_v(i)
             dvx(i,j,k)= ( a13(1)*vv(i-1,j,k)+a13(2)*vv(i  ,j,k) &
                          +a13(3)*vv(i+1,j,k)+a13(4)*vv(i+2,j,k) &
                          +a13(5)*vv(i+3,j,k) )*idx_v(i)
             dwx(i,j,k)= ( a13(1)*ww(i-1,j,k)+a13(2)*ww(i  ,j,k) &
                          +a13(3)*ww(i+1,j,k)+a13(4)*ww(i+2,j,k) &
                          +a13(5)*ww(i+3,j,k) )*idx_v(i)
          enddo
       enddo
    endif

    do k=1,nz
       do j=1,ny
          do i=ndx_v3,nfx_v3
             dux(i,j,k)= ( a5(1)*( uu(i+1,j,k)- uu(i-1,j,k)) &
                         + a5(2)*( uu(i+2,j,k)- uu(i-2,j,k)) )*idx_v(i)
             dvx(i,j,k)= ( a5(1)*( vv(i+1,j,k)- vv(i-1,j,k)) &
                         + a5(2)*( vv(i+2,j,k)- vv(i-2,j,k)) )*idx_v(i)
             dwx(i,j,k)= ( a5(1)*( ww(i+1,j,k)- ww(i-1,j,k)) &
                         + a5(2)*( ww(i+2,j,k)- ww(i-2,j,k)) )*idx_v(i)
          enddo
       enddo
    enddo

    if (is_boundary(1,2)) then
       i=nx-1
       do k=1,nz
          do j=1,ny
             dux(i,j,k)= ( a31(5)*uu(i-3,j,k)+a31(4)*uu(i-2,j,k) &
                          +a31(3)*uu(i-1,j,k)+a31(2)*uu(i  ,j,k) &
                          +a31(1)*uu(i+1,j,k) )*idx_v(i)
             dvx(i,j,k)= ( a31(5)*vv(i-3,j,k)+a31(4)*vv(i-2,j,k) &
                          +a31(3)*vv(i-1,j,k)+a31(2)*vv(i  ,j,k) &
                          +a31(1)*vv(i+1,j,k) )*idx_v(i)
             dwx(i,j,k)= ( a31(5)*ww(i-3,j,k)+a31(4)*ww(i-2,j,k) &
                          +a31(3)*ww(i-1,j,k)+a31(2)*ww(i  ,j,k) &
                          +a31(1)*ww(i+1,j,k) )*idx_v(i)
          enddo
       enddo

       i=nx
       do k=1,nz
          do j=1,ny
             dux(i,j,k)= ( a40(5)*uu(i-4,j,k)+a40(4)*uu(i-3,j,k) &
                          +a40(3)*uu(i-2,j,k)+a40(2)*uu(i-1,j,k) &
                          +a40(1)*uu(i  ,j,k) )*idx_v(i)
             dvx(i,j,k)= ( a40(5)*vv(i-4,j,k)+a40(4)*vv(i-3,j,k) &
                          +a40(3)*vv(i-2,j,k)+a40(2)*vv(i-1,j,k) &
                          +a40(1)*vv(i  ,j,k) )*idx_v(i)
             dwx(i,j,k)= ( a40(5)*ww(i-4,j,k)+a40(4)*ww(i-3,j,k) &
                          +a40(3)*ww(i-2,j,k)+a40(2)*ww(i-1,j,k) &
                          +a40(1)*ww(i  ,j,k) )*idx_v(i)
          enddo
       enddo
    endif

    ! Derivatives along y
    ! ===================
    if (is_boundary(2,1)) then
       j=1
       do k=1,nz
          do i=1,nx
             duy(i,j,k)= ( a04(1)*uu(i,j  ,k)+a04(2)*uu(i,j+1,k) &
                          +a04(3)*uu(i,j+2,k)+a04(4)*uu(i,j+3,k) &
                          +a04(5)*uu(i,j+4,k) )*idy_v(j)
             dvy(i,j,k)= ( a04(1)*vv(i,j  ,k)+a04(2)*vv(i,j+1,k) &
                          +a04(3)*vv(i,j+2,k)+a04(4)*vv(i,j+3,k) &
                          +a04(5)*vv(i,j+4,k) )*idy_v(j)
             dwy(i,j,k)= ( a04(1)*ww(i,j  ,k)+a04(2)*ww(i,j+1,k) &
                          +a04(3)*ww(i,j+2,k)+a04(4)*ww(i,j+3,k) &
                          +a04(5)*ww(i,j+4,k) )*idy_v(j)
          enddo
       enddo

       j=2
       do k=1,nz
          do i=1,nx
             duy(i,j,k)= ( a13(1)*uu(i,j-1,k)+a13(2)*uu(i,j  ,k) &
                          +a13(3)*uu(i,j+1,k)+a13(4)*uu(i,j+2,k) &
                          +a13(5)*uu(i,j+3,k) )*idy_v(j)
             dvy(i,j,k)= ( a13(1)*vv(i,j-1,k)+a13(2)*vv(i,j  ,k) &
                          +a13(3)*vv(i,j+1,k)+a13(4)*vv(i,j+2,k) &
                          +a13(5)*vv(i,j+3,k) )*idy_v(j)
             dwy(i,j,k)= ( a13(1)*ww(i,j-1,k)+a13(2)*ww(i,j  ,k) &
                          +a13(3)*ww(i,j+1,k)+a13(4)*ww(i,j+2,k) &
                          +a13(5)*ww(i,j+3,k) )*idy_v(j)
          enddo
       enddo
    endif

    do k=1,nz
       do j=ndy_v3,nfy_v3
          do i=1,nx
             duy(i,j,k)= ( a5(1)*( uu(i,j+1,k)- uu(i,j-1,k)) &
                         + a5(2)*( uu(i,j+2,k)- uu(i,j-2,k)) )*idy_v(j)
             dvy(i,j,k)= ( a5(1)*( vv(i,j+1,k)- vv(i,j-1,k)) &
                         + a5(2)*( vv(i,j+2,k)- vv(i,j-2,k)) )*idy_v(j)
             dwy(i,j,k)= ( a5(1)*( ww(i,j+1,k)- ww(i,j-1,k)) &
                         + a5(2)*( ww(i,j+2,k)- ww(i,j-2,k)) )*idy_v(j)
          enddo
       enddo
    enddo

    if (is_boundary(2,2)) then
       j=ny-1
       do k=1,nz
          do i=1,nx
             duy(i,j,k)= ( a31(5)*uu(i,j-3,k)+a31(4)*uu(i,j-2,k) &
                          +a31(3)*uu(i,j-1,k)+a31(2)*uu(i,j  ,k) &
                          +a31(1)*uu(i,j+1,k) )*idy_v(j)
             dvy(i,j,k)= ( a31(5)*vv(i,j-3,k)+a31(4)*vv(i,j-2,k) &
                          +a31(3)*vv(i,j-1,k)+a31(2)*vv(i,j  ,k) &
                          +a31(1)*vv(i,j+1,k) )*idy_v(j)
             dwy(i,j,k)= ( a31(5)*ww(i,j-3,k)+a31(4)*ww(i,j-2,k) &
                          +a31(3)*ww(i,j-1,k)+a31(2)*ww(i,j  ,k) &
                          +a31(1)*ww(i,j+1,k) )*idy_v(j)
          enddo
       enddo

       j=ny
       do k=1,nz
          do i=1,nx
             duy(i,j,k)= ( a40(5)*uu(i,j-4,k)+a40(4)*uu(i,j-3,k) &
                          +a40(3)*uu(i,j-2,k)+a40(2)*uu(i,j-1,k) &
                          +a40(1)*uu(i,j  ,k) )*idy_v(j)
             dvy(i,j,k)= ( a40(5)*vv(i,j-4,k)+a40(4)*vv(i,j-3,k) &
                          +a40(3)*vv(i,j-2,k)+a40(2)*vv(i,j-1,k) &
                          +a40(1)*vv(i,j  ,k) )*idy_v(j)
             dwy(i,j,k)= ( a40(5)*ww(i,j-4,k)+a40(4)*ww(i,j-3,k) &
                          +a40(3)*ww(i,j-2,k)+a40(2)*ww(i,j-1,k) &
                          +a40(1)*ww(i,j  ,k) )*idy_v(j)
          enddo
       enddo
    endif

    !****************
    if (is_2D) return
    !****************

    ! Derivatives along z
    ! ===================
    if (is_boundary(3,1)) then
       k=1
       do j=1,ny
          do i=1,nx
             duz(i,j,k)= ( a04(1)*uu(i,j,k  )+a04(2)*uu(i,j,k+1) &
                          +a04(3)*uu(i,j,k+2)+a04(4)*uu(i,j,k+3) &
                          +a04(5)*uu(i,j,k+4) )*idz_v(k)
             dvz(i,j,k)= ( a04(1)*vv(i,j,k  )+a04(2)*vv(i,j,k+1) &
                          +a04(3)*vv(i,j,k+2)+a04(4)*vv(i,j,k+3) &
                          +a04(5)*vv(i,j,k+4) )*idz_v(k)
             dwz(i,j,k)= ( a04(1)*ww(i,j,k  )+a04(2)*ww(i,j,k+1) &
                          +a04(3)*ww(i,j,k+2)+a04(4)*ww(i,j,k+3) &
                          +a04(5)*ww(i,j,k+4) )*idz_v(k)
          enddo
       enddo

       k=2
       do j=1,ny
          do i=1,nx
             duz(i,j,k)= ( a13(1)*uu(i,j,k-1)+a13(2)*uu(i,j,k  ) &
                          +a13(3)*uu(i,j,k+1)+a13(4)*uu(i,j,k+2) &
                          +a13(5)*uu(i,j,k+3) )*idz_v(k)
             dvz(i,j,k)= ( a13(1)*vv(i,j,k-1)+a13(2)*vv(i,j,k  ) &
                          +a13(3)*vv(i,j,k+1)+a13(4)*vv(i,j,k+2) &
                          +a13(5)*vv(i,j,k+3) )*idz_v(k)
             dwz(i,j,k)= ( a13(1)*ww(i,j,k-1)+a13(2)*ww(i,j,k  ) &
                          +a13(3)*ww(i,j,k+1)+a13(4)*ww(i,j,k+2) &
                          +a13(5)*ww(i,j,k+3) )*idz_v(k)
          enddo
       enddo
    endif

    do k=ndz_v3,nfz_v3
       do j=1,ny
          do i=1,nx
             duz(i,j,k)= ( a5(1)*( uu(i,j,k+1)- uu(i,j,k-1)) &
                         + a5(2)*( uu(i,j,k+2)- uu(i,j,k-2)) )*idz_v(k)
             dvz(i,j,k)= ( a5(1)*( vv(i,j,k+1)- vv(i,j,k-1)) &
                         + a5(2)*( vv(i,j,k+2)- vv(i,j,k-2)) )*idz_v(k)
             dwz(i,j,k)= ( a5(1)*( ww(i,j,k+1)- ww(i,j,k-1)) &
                         + a5(2)*( ww(i,j,k+2)- ww(i,j,k-2)) )*idz_v(k)
          enddo
       enddo
    enddo

    if (is_boundary(3,2)) then
       k=nz-1
       do j=1,ny
          do i=1,nx
             duz(i,j,k)= ( a31(5)*uu(i,j,k-3)+a31(4)*uu(i,j,k-2) &
                          +a31(3)*uu(i,j,k-1)+a31(2)*uu(i,j,k  ) &
                          +a31(1)*uu(i,j,k+1) )*idz_v(k)
             dvz(i,j,k)= ( a31(5)*vv(i,j,k-3)+a31(4)*vv(i,j,k-2) &
                          +a31(3)*vv(i,j,k-1)+a31(2)*vv(i,j,k  ) &
                          +a31(1)*vv(i,j,k+1) )*idz_v(k)
             dwz(i,j,k)= ( a31(5)*ww(i,j,k-3)+a31(4)*ww(i,j,k-2) &
                          +a31(3)*ww(i,j,k-1)+a31(2)*ww(i,j,k  ) &
                          +a31(1)*ww(i,j,k+1) )*idz_v(k)
          enddo
       enddo

       k=nz
       do j=1,ny
          do i=1,nx
             duz(i,j,k)= ( a40(5)*uu(i,j,k-4)+a40(4)*uu(i,j,k-3) &
                          +a40(3)*uu(i,j,k-2)+a40(2)*uu(i,j,k-1) &
                          +a40(1)*uu(i,j,k  ) )*idz_v(k)
             dvz(i,j,k)= ( a40(5)*vv(i,j,k-4)+a40(4)*vv(i,j,k-3) &
                          +a40(3)*vv(i,j,k-2)+a40(2)*vv(i,j,k-1) &
                          +a40(1)*vv(i,j,k  ) )*idz_v(k)
             dwz(i,j,k)= ( a40(5)*ww(i,j,k-4)+a40(4)*ww(i,j,k-3) &
                          +a40(3)*ww(i,j,k-2)+a40(2)*ww(i,j,k-1) &
                          +a40(1)*ww(i,j,k  ) )*idz_v(k)
          enddo
       enddo
    endif

  end subroutine grad_vel_5pts

  !===============================================================================
  module subroutine grad_vel_3pts
  !===============================================================================
    !> Compute velocity derivatives
    !> 3-point stencil - Cartesian version -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k
    ! ----------------------------------------------------------------------------

    ! Derivatives along x
    ! ===================
    if (is_boundary(1,1)) then
       i=1
       do k=1,nz
          do j=1,ny
             dux(i,j,k)= ( a02(1)*uu(i  ,j,k)+a02(2)*uu(i+1,j,k) &
                          +a02(3)*uu(i+2,j,k))*idx_v(i)
             dvx(i,j,k)= ( a02(1)*vv(i  ,j,k)+a02(2)*vv(i+1,j,k) &
                          +a02(3)*vv(i+2,j,k))*idx_v(i)
             dwx(i,j,k)= ( a02(1)*ww(i  ,j,k)+a02(2)*ww(i+1,j,k) &
                          +a02(3)*ww(i+2,j,k))*idx_v(i)
          enddo
       enddo
    endif

    do k=1,nz
       do j=1,ny
          do i=ndx_v3,nfx_v3
             dux(i,j,k)= 0.5_wp*(uu(i+1,j,k)-uu(i-1,j,k))*idx_v(i)
             dvx(i,j,k)= 0.5_wp*(vv(i+1,j,k)-vv(i-1,j,k))*idx_v(i)
             dwx(i,j,k)= 0.5_wp*(ww(i+1,j,k)-ww(i-1,j,k))*idx_v(i)
          enddo
       enddo
    enddo

    if (is_boundary(1,2)) then
       i=nx
       do k=1,nz
          do j=1,ny
             dux(i,j,k)= ( a20(3)*uu(i-2,j,k)+a20(2)*uu(i-1,j,k) &
                          +a20(1)*uu(i  ,j,k) )*idx_v(i)
             dvx(i,j,k)= ( a20(3)*vv(i-2,j,k)+a20(2)*vv(i-1,j,k) &
                          +a20(1)*vv(i  ,j,k) )*idx_v(i)
             dwx(i,j,k)= ( a20(3)*ww(i-2,j,k)+a20(2)*ww(i-1,j,k) &
                          +a20(1)*ww(i  ,j,k) )*idx_v(i)
          enddo
       enddo
    endif

    ! Derivatives along y
    ! ===================
    if (is_boundary(2,1)) then
       j=1
       do k=1,nz
          do i=1,nx
             duy(i,j,k)= ( a02(1)*uu(i,j  ,k)+a02(2)*uu(i,j+1,k) &
                          +a02(3)*uu(i,j+2,k))*idy_v(j)
             dvy(i,j,k)= ( a02(1)*vv(i,j  ,k)+a02(2)*vv(i,j+1,k) &
                          +a02(3)*vv(i,j+2,k))*idy_v(j)
             dwy(i,j,k)= ( a02(1)*ww(i,j  ,k)+a02(2)*ww(i,j+1,k) &
                          +a02(3)*ww(i,j+2,k))*idy_v(j)
          enddo
       enddo
    endif

    do k=1,nz
       do j=ndy_v3,nfy_v3
          do i=1,nx
             duy(i,j,k)= 0.5_wp*(uu(i,j+1,k)-uu(i,j-1,k))*idy_v(j)
             dvy(i,j,k)= 0.5_wp*(vv(i,j+1,k)-vv(i,j-1,k))*idy_v(j)
             dwy(i,j,k)= 0.5_wp*(ww(i,j+1,k)-ww(i,j-1,k))*idy_v(j)
          enddo
       enddo
    enddo

    if (is_boundary(2,2)) then
       j=ny
       do k=1,nz
          do i=1,nx
             duy(i,j,k)= ( a20(3)*uu(i,j-2,k)+a20(2)*uu(i,j-1,k) &
                          +a20(1)*uu(i,j  ,k) )*idy_v(j)
             dvy(i,j,k)= ( a20(3)*vv(i,j-2,k)+a20(2)*vv(i,j-1,k) &
                          +a20(1)*vv(i,j  ,k) )*idy_v(j)
             dwy(i,j,k)= ( a20(3)*ww(i,j-2,k)+a20(2)*ww(i,j-1,k) &
                          +a20(1)*ww(i,j  ,k) )*idy_v(j)
          enddo
       enddo
    endif

    !****************
    if (is_2D) return
    !****************

    ! Derivatives along z
    ! ===================
    if (is_boundary(3,1)) then
       k=1
       do j=1,ny
          do i=1,nx
             duz(i,j,k)= ( a02(1)*uu(i,j,k  )+a02(2)*uu(i,j,k+1) &
                          +a02(3)*uu(i,j,k+2))*idz_v(k)
             dvz(i,j,k)= ( a02(1)*vv(i,j,k  )+a02(2)*vv(i,j,k+1) &
                          +a02(3)*vv(i,j,k+2))*idz_v(k)
             dwz(i,j,k)= ( a02(1)*ww(i,j,k  )+a02(2)*ww(i,j,k+1) &
                          +a02(3)*ww(i,j,k+2))*idz_v(k)
          enddo
       enddo
    endif

    do k=ndz_v3,nfz_v3
       do j=1,ny
          do i=1,nx
             duz(i,j,k)= 0.5_wp*(uu(i,j,k+1)-uu(i,j,k-1))*idz_v(k)
             dvz(i,j,k)= 0.5_wp*(vv(i,j,k+1)-vv(i,j,k-1))*idz_v(k)
             dwz(i,j,k)= 0.5_wp*(ww(i,j,k+1)-ww(i,j,k-1))*idz_v(k)
          enddo
       enddo
    enddo

    if (is_boundary(3,2)) then
       k=nz
       do j=1,ny
          do i=1,nx
             duz(i,j,k)= ( a20(3)*uu(i,j,k-2)+a20(2)*uu(i,j,k-1) &
                          +a20(1)*uu(i,j,k  ) )*idz_v(k)
             dvz(i,j,k)= ( a20(3)*vv(i,j,k-2)+a20(2)*vv(i,j,k-1) &
                          +a20(1)*vv(i,j,k  ) )*idz_v(k)
             dwz(i,j,k)= ( a20(3)*ww(i,j,k-2)+a20(2)*ww(i,j,k-1) &
                          +a20(1)*ww(i,j,k  ) )*idz_v(k)
          enddo
       enddo
    endif

  end subroutine grad_vel_3pts

end submodule smod_grad_vel
