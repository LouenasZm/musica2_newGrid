!=================================================================================
submodule (mod_flux_visc) smod_flux_visc_3pts
!=================================================================================
  !> author: XG
  !> date: January 2024
  !> Cartesian version - compute viscous fluxes with 3-pt stencil
!=================================================================================

contains

  !==============================================================================
  module subroutine flux_visc_3pts
  !==============================================================================
    !> Compute derivatives of viscous fluxes (3-point stencil - order 2)
    !> - Cartesian version -
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: tau11,tau22,tau33,tau12,tau13,tau23,trace,mu
    real(wp), parameter :: ONE_THIRD=1.0_wp/3.0_wp
    ! ---------------------------------------------------------------------------

    ! Compute viscous fluxes
    ! ======================
    do k=ndzt_v,nfzt_v
       do j=ndyt_v,nfyt_v
          do i=ndxt_v,nfxt_v

             ! compute S_ij
             tau11 = dux(i,j,k)
             tau22 = dvy(i,j,k)
             tau33 = dwz(i,j,k)
             tau12 = 0.5_wp*(duy(i,j,k)+dvx(i,j,k))
             tau13 = 0.5_wp*(duz(i,j,k)+dwx(i,j,k))
             tau23 = 0.5_wp*(dvz(i,j,k)+dwy(i,j,k))
             trace = ONE_THIRD*(tau11+tau22+tau33)

             ! compute -tau_ij
             mu =-2.0_wp*visc(i,j,k)
             tau11=mu*(tau11-trace)
             tau22=mu*(tau22-trace)
             tau33=mu*(tau33-trace) 
             tau12=mu*tau12
             tau13=mu*tau13
             tau23=mu*tau23

             ! viscous fluxes along x
             Frhou(i,j,k)= tau11
             Frhov(i,j,k)= tau12
             Frhow(i,j,k)= tau13
             Frhoe(i,j,k)= (uu(i,j,k)*Frhou(i,j,k)+vv(i,j,k)*Frhov(i,j,k)+ww(i,j,k)*Frhow(i,j,k)) &
                  -cok(i,j,k)*dTx(i,j,k)

             ! viscous fluxes along y
             Grhou(i,j,k)=tau12
             Grhov(i,j,k)=tau22
             Grhow(i,j,k)=tau23
             Grhoe(i,j,k)= (uu(i,j,k)*Grhou(i,j,k)+vv(i,j,k)*Grhov(i,j,k)+ww(i,j,k)*Grhow(i,j,k)) &
                  -cok(i,j,k)*dTy(i,j,k)

             ! viscous fluxes along z
             Hrhou(i,j,k)= tau13
             Hrhov(i,j,k)= tau23
             Hrhow(i,j,k)= tau33
             Hrhoe(i,j,k)= (uu(i,j,k)*Hrhou(i,j,k)+vv(i,j,k)*Hrhov(i,j,k)+ww(i,j,k)*Hrhow(i,j,k)) &
                  -cok(i,j,k)*dTz(i,j,k)
          enddo
       enddo
    enddo

    ! Compute derivatives of viscous fluxes along x
    ! =============================================
    if (is_bc_wall(1,1)) then
       i=1
       do k=ndz_v,nfz_v
          do j=ndy_v,nfy_v
             Krhou(i,j,k) = ( a02(1)*Frhou(i,j,k) + a02(2)*Frhou(i+1,j,k) &
                            + a02(3)*Frhou(i+2,j,k))*idx_v(i) + Krhou(i,j,k)
             Krhov(i,j,k) = ( a02(1)*Frhov(i,j,k) + a02(2)*Frhov(i+1,j,k) &
                            + a02(3)*Frhov(i+2,j,k))*idx_v(i) + Krhov(i,j,k)
             Krhow(i,j,k) = ( a02(1)*Frhow(i,j,k) + a02(2)*Frhow(i+1,j,k) &
                            + a02(3)*Frhow(i+2,j,k))*idx_v(i) + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a02(1)*Frhoe(i,j,k) + a02(2)*Frhoe(i+1,j,k) &
                            + a02(3)*Frhoe(i+2,j,k))*idx_v(i) + Krhoe(i,j,k)
          enddo
       enddo
    endif

    do k=ndz_v,nfz_v
       do j=ndy_v,nfy_v
          do i=ndx_vi,nfx_vi
             Krhou(i,j,k) = a3(1)*( Frhou(i+1,j,k)-Frhou(i-1,j,k) )*idx_v(i) + Krhou(i,j,k)
             Krhov(i,j,k) = a3(1)*( Frhov(i+1,j,k)-Frhov(i-1,j,k) )*idx_v(i) + Krhov(i,j,k)
             Krhow(i,j,k) = a3(1)*( Frhow(i+1,j,k)-Frhow(i-1,j,k) )*idx_v(i) + Krhow(i,j,k)
             Krhoe(i,j,k) = a3(1)*( Frhoe(i+1,j,k)-Frhoe(i-1,j,k) )*idx_v(i) + Krhoe(i,j,k)
          enddo
       enddo
    enddo

    if (is_bc_wall(1,2)) then
       i=nx
       do k=ndz_v,nfz_v
          do j=ndy_v,nfy_v
             Krhou(i,j,k) = ( a20(1)*Frhou(i,j,k) + a20(2)*Frhou(i-1,j,k) &
                            + a20(3)*Frhou(i-2,j,k))*idx_v(i) + Krhou(i,j,k)
             Krhov(i,j,k) = ( a20(1)*Frhov(i,j,k) + a20(2)*Frhov(i-1,j,k) &
                            + a20(3)*Frhov(i-2,j,k))*idx_v(i) + Krhov(i,j,k)
             Krhow(i,j,k) = ( a20(1)*Frhow(i,j,k) + a20(2)*Frhow(i-1,j,k) &
                            + a20(3)*Frhow(i-2,j,k))*idx_v(i) + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a20(1)*Frhoe(i,j,k) + a20(2)*Frhoe(i-1,j,k) &
                            + a20(3)*Frhoe(i-2,j,k))*idx_v(i) + Krhoe(i,j,k)
          enddo
       enddo
    endif

    ! Compute derivatives of viscous fluxes along y
    ! =============================================
    if (is_bc_wall(2,1)) then
       j=1
       do k=ndz_v,nfz_v
          do i=ndx_v,nfx_v
             Krhou(i,j,k) = ( a02(1)*Grhou(i,j,k) + a02(2)*Grhou(i,j+1,k) &
                            + a02(3)*Grhou(i,j+2,k))*idy_v(j) + Krhou(i,j,k)
             Krhov(i,j,k) = ( a02(1)*Grhov(i,j,k) + a02(2)*Grhov(i,j+1,k) &
                            + a02(3)*Grhov(i,j+2,k))*idy_v(j) + Krhov(i,j,k)
             Krhow(i,j,k) = ( a02(1)*Grhow(i,j,k) + a02(2)*Grhow(i,j+1,k) &
                            + a02(3)*Grhow(i,j+2,k))*idy_v(j) + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a02(1)*Grhoe(i,j,k) + a02(2)*Grhoe(i,j+1,k) &
                            + a02(3)*Grhoe(i,j+2,k))*idy_v(j) + Krhoe(i,j,k)
          enddo
       enddo
    endif

    do k=ndz_v,nfz_v
       do j=ndy_vi,nfy_vi
          do i=ndx_v,nfx_v
             Krhou(i,j,k) = a3(1)*( Grhou(i,j+1,k)-Grhou(i,j-1,k) )*idy_v(j) + Krhou(i,j,k)
             Krhov(i,j,k) = a3(1)*( Grhov(i,j+1,k)-Grhov(i,j-1,k) )*idy_v(j) + Krhov(i,j,k)
             Krhow(i,j,k) = a3(1)*( Grhow(i,j+1,k)-Grhow(i,j-1,k) )*idy_v(j) + Krhow(i,j,k)
             Krhoe(i,j,k) = a3(1)*( Grhoe(i,j+1,k)-Grhoe(i,j-1,k) )*idy_v(j) + Krhoe(i,j,k)
          enddo
       enddo
    enddo

    if (is_bc_wall(2,2)) then
       j=ny
       do k=ndz_v,nfz_v
          do i=ndx_v,nfx_v
             Krhou(i,j,k) = ( a20(1)*Grhou(i,j,k) + a20(2)*Grhou(i,j-1,k) &
                            + a20(3)*Grhou(i,j-2,k))*idy_v(j) + Krhou(i,j,k)
             Krhov(i,j,k) = ( a20(1)*Grhov(i,j,k) + a20(2)*Grhov(i,j-1,k) &
                            + a20(3)*Grhov(i,j-2,k))*idy_v(j) + Krhov(i,j,k)
             Krhow(i,j,k) = ( a20(1)*Grhow(i,j,k) + a20(2)*Grhow(i,j-1,k) &
                            + a20(3)*Grhow(i,j-2,k))*idy_v(j) + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a20(1)*Grhoe(i,j,k) + a20(2)*Grhoe(i,j-1,k) &
                            + a20(3)*Grhoe(i,j-2,k))*idy_v(j) + Krhoe(i,j,k)
          enddo
       enddo
    endif

    if (is_2D) return

    ! Compute derivatives of viscous fluxes along z
    ! =============================================
    if (is_bc_wall(3,1)) then
       k=1
       do j=ndy_v,nfy_v
          do i=ndx_v,nfx_v
             Krhou(i,j,k) = ( a02(1)*Hrhou(i,j,k  )+a02(2)*Hrhou(i,j,k+1) &
                            + a02(3)*Hrhou(i,j,k+2))*idz_v(k) + Krhou(i,j,k)
             Krhov(i,j,k) = ( a02(1)*Hrhov(i,j,k  )+a02(2)*Hrhov(i,j,k+1) &
                            + a02(3)*Hrhov(i,j,k+2))*idz_v(k) + Krhov(i,j,k)
             Krhow(i,j,k) = ( a02(1)*Hrhow(i,j,k  )+a02(2)*Hrhow(i,j,k+1) &
                            + a02(3)*Hrhow(i,j,k+2))*idz_v(k) + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a02(1)*Hrhoe(i,j,k  )+a02(2)*Hrhoe(i,j,k+1) &
                            + a02(3)*Hrhoe(i,j,k+2))*idz_v(k) + Krhoe(i,j,k)
          enddo
       enddo
    endif

    do k=ndz_vi,nfz_vi
       do j=ndy_v,nfy_v
          do i=ndx_v,nfx_v
             Krhou(i,j,k) = a3(1)*( Hrhou(i,j,k+1)-Hrhou(i,j,k-1) )*idz_v(k) + Krhou(i,j,k)
             Krhov(i,j,k) = a3(1)*( Hrhov(i,j,k+1)-Hrhov(i,j,k-1) )*idz_v(k) + Krhov(i,j,k)
             Krhow(i,j,k) = a3(1)*( Hrhow(i,j,k+1)-Hrhow(i,j,k-1) )*idz_v(k) + Krhow(i,j,k)
             Krhoe(i,j,k) = a3(1)*( Hrhoe(i,j,k+1)-Hrhoe(i,j,k-1) )*idz_v(k) + Krhoe(i,j,k)
          enddo
       enddo
    enddo

    if (is_bc_wall(3,2)) then
       k=nz
       do j=ndy_v,nfy_v
          do i=ndx_v,nfx_v
             Krhou(i,j,k) = ( a20(1)*Hrhou(i,j,k)  +a20(2)*Hrhou(i,j,k-1) &
                            + a20(3)*Hrhou(i,j,k-2))*idz_v(k) + Krhou(i,j,k)
             Krhov(i,j,k) = ( a20(1)*Hrhov(i,j,k)  +a20(2)*Hrhov(i,j,k-1) &
                            + a20(3)*Hrhov(i,j,k-2))*idz_v(k) + Krhov(i,j,k)
             Krhow(i,j,k) = ( a20(1)*Hrhow(i,j,k)  +a20(2)*Hrhow(i,j,k-1) &
                            + a20(3)*Hrhow(i,j,k-2))*idz_v(k) + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a20(1)*Hrhoe(i,j,k)  +a20(2)*Hrhoe(i,j,k-1) &
                            + a20(3)*Hrhoe(i,j,k-2))*idz_v(k) + Krhoe(i,j,k)
          enddo
       enddo
    endif

  end subroutine flux_visc_3pts

end submodule smod_flux_visc_3pts
