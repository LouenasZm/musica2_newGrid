!=================================================================================
submodule (mod_flux_visc) smod_flux_visc_5pts_c
!=================================================================================
  !> author: XG
  !> date: January 2024
  !> 2.5D curvilinear version - compute viscous fluxes with 5-pt stencil
!=================================================================================

contains

  !==============================================================================
  module subroutine flux_visc_5pts_c
  !==============================================================================
    !> Compute derivatives of viscous fluxes (5-point stencil - order 4)
    !> - 2.5D curvilinear version -
  !==============================================================================
    use mod_rans
    use mod_eos
    implicit none
    ! ---------------------------------------------------------------------------
    integer  :: i,j,k
    real(wp) :: tau11,tau22,tau33,tau12,tau13,tau23,trace,mu,cnu13
    ! ---------------------------------------------------------------------------
    real(wp), dimension(ndxt_v:nfxt_v,ndyt_v:nfyt_v,ndzt_v:nfzt_v) :: tauR11,tauR12,tauR13
    real(wp), dimension(ndxt_v:nfxt_v,ndyt_v:nfyt_v,ndzt_v:nfzt_v) :: tauR22,tauR23,tauR33,cokt
    ! ---------------------------------------------------------------------------

    ! Compute viscous fluxes
    ! ======================
    if (.not.is_RANS) then
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

                ! viscous fluxes along ksi
                Frhou(i,j,k)= tau11*y_eta_v(i,j)-tau12*x_eta_v(i,j)
                Frhov(i,j,k)= tau12*y_eta_v(i,j)-tau22*x_eta_v(i,j)
                Frhow(i,j,k)= tau13*y_eta_v(i,j)-tau23*x_eta_v(i,j)
                Frhoe(i,j,k)= (uu(i,j,k)*Frhou(i,j,k)+vv(i,j,k)*Frhov(i,j,k)+ww(i,j,k)*Frhow(i,j,k)) &
                              -cok(i,j,k)*(dTx(i,j,k)*y_eta_v(i,j)-dTy(i,j,k)*x_eta_v(i,j))

                ! viscous fluxes along eta
                Grhou(i,j,k)=-tau11*y_ksi_v(i,j)+tau12*x_ksi_v(i,j)
                Grhov(i,j,k)=-tau12*y_ksi_v(i,j)+tau22*x_ksi_v(i,j)
                Grhow(i,j,k)=-tau13*y_ksi_v(i,j)+tau23*x_ksi_v(i,j)
                Grhoe(i,j,k)= (uu(i,j,k)*Grhou(i,j,k)+vv(i,j,k)*Grhov(i,j,k)+ww(i,j,k)*Grhow(i,j,k)) &
                              -cok(i,j,k)*(-dTx(i,j,k)*y_ksi_v(i,j)+dTy(i,j,k)*x_ksi_v(i,j))

                ! viscous fluxes along z
                Hrhou(i,j,k)= tau13
                Hrhov(i,j,k)= tau23
                Hrhow(i,j,k)= tau33
                Hrhoe(i,j,k)= (uu(i,j,k)*Hrhou(i,j,k)+vv(i,j,k)*Hrhov(i,j,k)+ww(i,j,k)*Hrhow(i,j,k)) &
                              -cok(i,j,k)*dTz(i,j,k)
             enddo
          enddo
       enddo

    else
       
       if (model_RANS.eq.'SA') then
          cnu1 = 7.1_wp
          cnu13 = (7.1_wp)**3
          ! prrt = 1.0_wp

          do k=ndzt_v,nfzt_v
             do j=ndyt_v,nfyt_v
                do i=ndxt_v,nfxt_v
                   ! nutilde wall function
                   khi  = rho_n(i,j,k)*nutil_n(i,j,k)/visc(i,j,k)
                   fnu1 = khi**3/(khi**3+cnu13)

                   ! compute Reynolds stress
                   tauR11(i,j,k)  = 2*rho_n(i,j,k)*nutil_n(i,j,k)*fnu1* dux(i,j,k)
                   tauR22(i,j,k)  = 2*rho_n(i,j,k)*nutil_n(i,j,k)*fnu1* dvy(i,j,k)
                   tauR33(i,j,k)  = 2*rho_n(i,j,k)*nutil_n(i,j,k)*fnu1* dwz(i,j,k)
                   tauR12(i,j,k)  =   rho_n(i,j,k)*nutil_n(i,j,k)*fnu1*(duy(i,j,k)+dvx(i,j,k))
                   tauR13(i,j,k)  =   rho_n(i,j,k)*nutil_n(i,j,k)*fnu1*(duz(i,j,k)+dwx(i,j,k))
                   tauR23(i,j,k)  =   rho_n(i,j,k)*nutil_n(i,j,k)*fnu1*(dvz(i,j,k)+dwy(i,j,k))

                   ! turbulent heat transfer coefficient
                   cokt(i,j,k) = rho_n(i,j,k)*nutil_n(i,j,k)*fnu1*cpcalc_tro(Tmp(i,j,k),rho_n(i,j,k))!/prrt
                enddo
             enddo
          enddo

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

                   ! compute -(tau_ij+tauR_ij)
                   mu =-2.0_wp*visc(i,j,k)
                   tau11=mu*(tau11-trace)-tauR11(i,j,k)
                   tau22=mu*(tau22-trace)-tauR22(i,j,k)
                   tau33=mu*(tau33-trace)-tauR33(i,j,k)
                   tau12=mu*tau12-tauR12(i,j,k)
                   tau13=mu*tau13-tauR13(i,j,k)
                   tau23=mu*tau23-tauR23(i,j,k)

                   ! viscous fluxes along ksi
                   Frhou(i,j,k)= tau11*y_eta_v(i,j)-tau12*x_eta_v(i,j)
                   Frhov(i,j,k)= tau12*y_eta_v(i,j)-tau22*x_eta_v(i,j)
                   Frhow(i,j,k)= tau13*y_eta_v(i,j)-tau23*x_eta_v(i,j)
                   Frhoe(i,j,k)= (uu(i,j,k)*Frhou(i,j,k)+vv(i,j,k)*Frhov(i,j,k)+ww(i,j,k)*Frhow(i,j,k)) &
                               -(cok(i,j,k)+cokt(i,j,k))*(dTx(i,j,k)*y_eta_v(i,j)-dTy(i,j,k)*x_eta_v(i,j))

                   ! viscous fluxes along eta
                   Grhou(i,j,k)=-tau11*y_ksi_v(i,j)+tau12*x_ksi_v(i,j)
                   Grhov(i,j,k)=-tau12*y_ksi_v(i,j)+tau22*x_ksi_v(i,j)
                   Grhow(i,j,k)=-tau13*y_ksi_v(i,j)+tau23*x_ksi_v(i,j)
                   Grhoe(i,j,k)= (uu(i,j,k)*Grhou(i,j,k)+vv(i,j,k)*Grhov(i,j,k)+ww(i,j,k)*Grhow(i,j,k)) &
                               -(cok(i,j,k)+cokt(i,j,k))*(-dTx(i,j,k)*y_ksi_v(i,j)+dTy(i,j,k)*x_ksi_v(i,j))

                   ! viscous fluxes along z
                   Hrhou(i,j,k)= tau13
                   Hrhov(i,j,k)= tau23
                   Hrhow(i,j,k)= tau33
                   Hrhoe(i,j,k)= (uu(i,j,k)*Hrhou(i,j,k)+vv(i,j,k)*Hrhov(i,j,k)+ww(i,j,k)*Hrhow(i,j,k)) &
                               -(cok(i,j,k)+cokt(i,j,k))*dTz(i,j,k)
                enddo
             enddo
          enddo
       endif
    endif

    ! Compute derivatives of viscous fluxes along ksi
    ! ===============================================
    if (is_bc_wall(1,1)) then
       i=1
       do k=ndz_v,nfz_v
          do j=ndy_v,nfy_v
             Krhou(i,j,k) = ( a04(1)*Frhou(i  ,j,k)+a04(2)*Frhou(i+1,j,k) &
                            + a04(3)*Frhou(i+2,j,k)+a04(4)*Frhou(i+3,j,k) &
                            + a04(5)*Frhou(i+4,j,k) )*ijacob_v(i,j) + Krhou(i,j,k)
             Krhov(i,j,k) = ( a04(1)*Frhov(i  ,j,k)+a04(2)*Frhov(i+1,j,k) &
                            + a04(3)*Frhov(i+2,j,k)+a04(4)*Frhov(i+3,j,k) &
                            + a04(5)*Frhov(i+4,j,k) )*ijacob_v(i,j) + Krhov(i,j,k)
             Krhow(i,j,k) = ( a04(1)*Frhow(i  ,j,k)+a04(2)*Frhow(i+1,j,k) &
                            + a04(3)*Frhow(i+2,j,k)+a04(4)*Frhow(i+3,j,k) &
                            + a04(5)*Frhow(i+4,j,k) )*ijacob_v(i,j) + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a04(1)*Frhoe(i  ,j,k)+a04(2)*Frhoe(i+1,j,k) &
                            + a04(3)*Frhoe(i+2,j,k)+a04(4)*Frhoe(i+3,j,k) &
                            + a04(5)*Frhoe(i+4,j,k) )*ijacob_v(i,j) + Krhoe(i,j,k)
          enddo
       enddo
    endif
    if (is_bc_1pt(1,1)) then
       i=2
       do k=ndz_v,nfz_v
          do j=ndy_v,nfy_v
             Krhou(i,j,k) = ( a13(1)*Frhou(i-1,j,k)+a13(2)*Frhou(i,j,k) &
                            + a13(3)*Frhou(i+1,j,k)+a13(4)*Frhou(i+2,j,k) &
                            + a13(5)*Frhou(i+3,j,k) )*ijacob_v(i,j) + Krhou(i,j,k)
             Krhov(i,j,k) = ( a13(1)*Frhov(i-1,j,k)+a13(2)*Frhov(i,j,k) &
                            + a13(3)*Frhov(i+1,j,k)+a13(4)*Frhov(i+2,j,k) &
                            + a13(5)*Frhov(i+3,j,k) )*ijacob_v(i,j) + Krhov(i,j,k)
             Krhow(i,j,k) = ( a13(1)*Frhow(i-1,j,k)+a13(2)*Frhow(i,j,k) &
                            + a13(3)*Frhow(i+1,j,k)+a13(4)*Frhow(i+2,j,k) &
                            + a13(5)*Frhow(i+3,j,k) )*ijacob_v(i,j) + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a13(1)*Frhoe(i-1,j,k)+a13(2)*Frhoe(i,j,k) &
                            + a13(3)*Frhoe(i+1,j,k)+a13(4)*Frhoe(i+2,j,k) &
                            + a13(5)*Frhoe(i+3,j,k) )*ijacob_v(i,j) + Krhoe(i,j,k)
          enddo
       enddo
    endif

    do k=ndz_v,nfz_v
       do j=ndy_v,nfy_v
          do i=ndx_vi,nfx_vi
             Krhou(i,j,k) = ( a5(1)*( Frhou(i+1,j,k)-Frhou(i-1,j,k) ) &
                            + a5(2)*( Frhou(i+2,j,k)-Frhou(i-2,j,k) ) )*ijacob_v(i,j) + Krhou(i,j,k)
             Krhov(i,j,k) = ( a5(1)*( Frhov(i+1,j,k)-Frhov(i-1,j,k) ) &
                            + a5(2)*( Frhov(i+2,j,k)-Frhov(i-2,j,k) ) )*ijacob_v(i,j) + Krhov(i,j,k)
             Krhow(i,j,k) = ( a5(1)*( Frhow(i+1,j,k)-Frhow(i-1,j,k) ) &
                            + a5(2)*( Frhow(i+2,j,k)-Frhow(i-2,j,k) ) )*ijacob_v(i,j) + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a5(1)*( Frhoe(i+1,j,k)-Frhoe(i-1,j,k) ) &
                            + a5(2)*( Frhoe(i+2,j,k)-Frhoe(i-2,j,k) ) )*ijacob_v(i,j) + Krhoe(i,j,k)
          enddo
       enddo
    enddo

    if (is_bc_1pt(1,2)) then
       i=nx-1
       do k=ndz_v,nfz_v
          do j=ndy_v,nfy_v
             Krhou(i,j,k) = ( a31(1)*Frhou(i+1,j,k)+a31(2)*Frhou(i,j,k) &
                            + a31(3)*Frhou(i-1,j,k)+a31(4)*Frhou(i-2,j,k) &
                            + a31(5)*Frhou(i-3,j,k) )*ijacob_v(i,j) + Krhou(i,j,k)
             Krhov(i,j,k) = ( a31(1)*Frhov(i+1,j,k)+a31(2)*Frhov(i,j,k) &
                            + a31(3)*Frhov(i-1,j,k)+a31(4)*Frhov(i-2,j,k) &
                            + a31(5)*Frhov(i-3,j,k) )*ijacob_v(i,j) + Krhov(i,j,k)
             Krhow(i,j,k) = ( a31(1)*Frhow(i+1,j,k)+a31(2)*Frhow(i,j,k) &
                            + a31(3)*Frhow(i-1,j,k)+a31(4)*Frhow(i-2,j,k) &
                            + a31(5)*Frhow(i-3,j,k) )*ijacob_v(i,j) + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a31(1)*Frhoe(i+1,j,k)+a31(2)*Frhoe(i,j,k) &
                            + a31(3)*Frhoe(i-1,j,k)+a31(4)*Frhoe(i-2,j,k) &
                            + a31(5)*Frhoe(i-3,j,k) )*ijacob_v(i,j) + Krhoe(i,j,k)
          enddo
       enddo
    endif
    if (is_bc_wall(1,2)) then
       i=nx
       do k=ndz_v,nfz_v
          do j=ndy_v,nfy_v
             Krhou(i,j,k) = ( a40(1)*Frhou(i,j,k)  +a40(2)*Frhou(i-1,j,k) &
                            + a40(3)*Frhou(i-2,j,k)+a40(4)*Frhou(i-3,j,k) &
                            + a40(5)*Frhou(i-4,j,k) )*ijacob_v(i,j) + Krhou(i,j,k)
             Krhov(i,j,k) = ( a40(1)*Frhov(i,j,k)  +a40(2)*Frhov(i-1,j,k) &
                            + a40(3)*Frhov(i-2,j,k)+a40(4)*Frhov(i-3,j,k) &
                            + a40(5)*Frhov(i-4,j,k) )*ijacob_v(i,j) + Krhov(i,j,k)
             Krhow(i,j,k) = ( a40(1)*Frhow(i,j,k)  +a40(2)*Frhow(i-1,j,k) &
                            + a40(3)*Frhow(i-2,j,k)+a40(4)*Frhow(i-3,j,k) &
                            + a40(5)*Frhow(i-4,j,k) )*ijacob_v(i,j) + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a40(1)*Frhoe(i,j,k)  +a40(2)*Frhoe(i-1,j,k) &
                            + a40(3)*Frhoe(i-2,j,k)+a40(4)*Frhoe(i-3,j,k) &
                            + a40(5)*Frhoe(i-4,j,k) )*ijacob_v(i,j) + Krhoe(i,j,k)
          enddo
       enddo
    endif

    ! Compute derivatives of viscous fluxes along eta
    ! ===============================================
    if (is_bc_wall(2,1)) then
       j=1
       do k=ndz_v,nfz_v
          do i=ndx_v,nfx_v
             Krhou(i,j,k) = ( a04(1)*Grhou(i,j  ,k)+a04(2)*Grhou(i,j+1,k) &
                            + a04(3)*Grhou(i,j+2,k)+a04(4)*Grhou(i,j+3,k) &
                            + a04(5)*Grhou(i,j+4,k) )*ijacob_v(i,j) + Krhou(i,j,k)
             Krhov(i,j,k) = ( a04(1)*Grhov(i,j  ,k)+a04(2)*Grhov(i,j+1,k) &
                            + a04(3)*Grhov(i,j+2,k)+a04(4)*Grhov(i,j+3,k) &
                            + a04(5)*Grhov(i,j+4,k) )*ijacob_v(i,j) + Krhov(i,j,k)
             Krhow(i,j,k) = ( a04(1)*Grhow(i,j  ,k)+a04(2)*Grhow(i,j+1,k) &
                            + a04(3)*Grhow(i,j+2,k)+a04(4)*Grhow(i,j+3,k) &
                            + a04(5)*Grhow(i,j+4,k) )*ijacob_v(i,j) + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a04(1)*Grhoe(i,j  ,k)+a04(2)*Grhoe(i,j+1,k) &
                            + a04(3)*Grhoe(i,j+2,k)+a04(4)*Grhoe(i,j+3,k) &
                            + a04(5)*Grhoe(i,j+4,k) )*ijacob_v(i,j) + Krhoe(i,j,k)
          enddo
       enddo
    endif
    if (is_bc_1pt(2,1)) then
       j=2
       do k=ndz_v,nfz_v
          do i=ndx_v,nfx_v
             Krhou(i,j,k) = ( a13(1)*Grhou(i,j-1,k)+a13(2)*Grhou(i,j,k) &
                            + a13(3)*Grhou(i,j+1,k)+a13(4)*Grhou(i,j+2,k) &
                            + a13(5)*Grhou(i,j+3,k) )*ijacob_v(i,j) + Krhou(i,j,k)
             Krhov(i,j,k) = ( a13(1)*Grhov(i,j-1,k)+a13(2)*Grhov(i,j,k) &
                            + a13(3)*Grhov(i,j+1,k)+a13(4)*Grhov(i,j+2,k) &
                            + a13(5)*Grhov(i,j+3,k) )*ijacob_v(i,j) + Krhov(i,j,k)
             Krhow(i,j,k) = ( a13(1)*Grhow(i,j-1,k)+a13(2)*Grhow(i,j,k) &
                            + a13(3)*Grhow(i,j+1,k)+a13(4)*Grhow(i,j+2,k) &
                            + a13(5)*Grhow(i,j+3,k) )*ijacob_v(i,j) + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a13(1)*Grhoe(i,j-1,k)+a13(2)*Grhoe(i,j,k) &
                            + a13(3)*Grhoe(i,j+1,k)+a13(4)*Grhoe(i,j+2,k) &
                            + a13(5)*Grhoe(i,j+3,k) )*ijacob_v(i,j) + Krhoe(i,j,k)
          enddo
       enddo
    endif

    do k=ndz_v,nfz_v
       do j=ndy_vi,nfy_vi
          do i=ndx_v,nfx_v
             Krhou(i,j,k) = ( a5(1)*( Grhou(i,j+1,k)-Grhou(i,j-1,k)) &
                            + a5(2)*( Grhou(i,j+2,k)-Grhou(i,j-2,k)) )*ijacob_v(i,j) + Krhou(i,j,k)
             Krhov(i,j,k) = ( a5(1)*( Grhov(i,j+1,k)-Grhov(i,j-1,k)) &
                            + a5(2)*( Grhov(i,j+2,k)-Grhov(i,j-2,k)) )*ijacob_v(i,j) + Krhov(i,j,k)
             Krhow(i,j,k) = ( a5(1)*( Grhow(i,j+1,k)-Grhow(i,j-1,k)) &
                            + a5(2)*( Grhow(i,j+2,k)-Grhow(i,j-2,k)) )*ijacob_v(i,j) + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a5(1)*( Grhoe(i,j+1,k)-Grhoe(i,j-1,k)) &
                            + a5(2)*( Grhoe(i,j+2,k)-Grhoe(i,j-2,k)) )*ijacob_v(i,j) + Krhoe(i,j,k)
          enddo
       enddo
    enddo

    if (is_bc_1pt(2,2)) then
       j=ny-1
       do k=ndz_v,nfz_v
          do i=ndx_v,nfx_v
             Krhou(i,j,k) = ( a31(1)*Grhou(i,j+1,k)+a31(2)*Grhou(i,j,k) &
                            + a31(3)*Grhou(i,j-1,k)+a31(4)*Grhou(i,j-2,k) &
                            + a31(5)*Grhou(i,j-3,k) )*ijacob_v(i,j) + Krhou(i,j,k)
             Krhov(i,j,k) = ( a31(1)*Grhov(i,j+1,k)+a31(2)*Grhov(i,j,k) &
                            + a31(3)*Grhov(i,j-1,k)+a31(4)*Grhov(i,j-2,k) &
                            + a31(5)*Grhov(i,j-3,k) )*ijacob_v(i,j) + Krhov(i,j,k)
             Krhow(i,j,k) = ( a31(1)*Grhow(i,j+1,k)+a31(2)*Grhow(i,j,k) &
                            + a31(3)*Grhow(i,j-1,k)+a31(4)*Grhow(i,j-2,k) &
                            + a31(5)*Grhow(i,j-3,k) )*ijacob_v(i,j) + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a31(1)*Grhoe(i,j+1,k)+a31(2)*Grhoe(i,j,k) &
                            + a31(3)*Grhoe(i,j-1,k)+a31(4)*Grhoe(i,j-2,k) &
                            + a31(5)*Grhoe(i,j-3,k) )*ijacob_v(i,j) + Krhoe(i,j,k)
          enddo
       enddo
    endif
    if (is_bc_wall(2,2)) then
       j=ny
       do k=ndz_v,nfz_v
          do i=ndx_v,nfx_v
             Krhou(i,j,k) = ( a40(1)*Grhou(i,j,k)  +a40(2)*Grhou(i,j-1,k) &
                            + a40(3)*Grhou(i,j-2,k)+a40(4)*Grhou(i,j-3,k) &
                            + a40(5)*Grhou(i,j-4,k) )*ijacob_v(i,j) + Krhou(i,j,k)
             Krhov(i,j,k) = ( a40(1)*Grhov(i,j,k)  +a40(2)*Grhov(i,j-1,k) &
                            + a40(3)*Grhov(i,j-2,k)+a40(4)*Grhov(i,j-3,k) &
                            + a40(5)*Grhov(i,j-4,k) )*ijacob_v(i,j) + Krhov(i,j,k)
             Krhow(i,j,k) = ( a40(1)*Grhow(i,j,k)  +a40(2)*Grhow(i,j-1,k) &
                            + a40(3)*Grhow(i,j-2,k)+a40(4)*Grhow(i,j-3,k) &
                            + a40(5)*Grhow(i,j-4,k) )*ijacob_v(i,j) + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a40(1)*Grhoe(i,j,k)  +a40(2)*Grhoe(i,j-1,k) &
                            + a40(3)*Grhoe(i,j-2,k)+a40(4)*Grhoe(i,j-3,k) &
                            + a40(5)*Grhoe(i,j-4,k) )*ijacob_v(i,j) + Krhoe(i,j,k)
          enddo
       enddo
    endif

    !****************
    if (is_2D) return
    !****************

    ! Compute derivatives of viscous fluxes along z
    ! =============================================
    if (is_bc_wall(3,1)) then
       k=1
       do j=ndy_v,nfy_v
          do i=ndx_v,nfx_v
             Krhou(i,j,k) = ( a04(1)*Hrhou(i,j,k  )+a04(2)*Hrhou(i,j,k+1) &
                            + a04(3)*Hrhou(i,j,k+2)+a04(4)*Hrhou(i,j,k+3) &
                            + a04(5)*Hrhou(i,j,k+4) )*idz_v(k) + Krhou(i,j,k)
             Krhov(i,j,k) = ( a04(1)*Hrhov(i,j,k  )+a04(2)*Hrhov(i,j,k+1) &
                            + a04(3)*Hrhov(i,j,k+2)+a04(4)*Hrhov(i,j,k+3) &
                            + a04(5)*Hrhov(i,j,k+4) )*idz_v(k) + Krhov(i,j,k)
             Krhow(i,j,k) = ( a04(1)*Hrhow(i,j,k  )+a04(2)*Hrhow(i,j,k+1) &
                            + a04(3)*Hrhow(i,j,k+2)+a04(4)*Hrhow(i,j,k+3) &
                            + a04(5)*Hrhow(i,j,k+4) )*idz_v(k) + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a04(1)*Hrhoe(i,j,k  )+a04(2)*Hrhoe(i,j,k+1) &
                            + a04(3)*Hrhoe(i,j,k+2)+a04(4)*Hrhoe(i,j,k+3) &
                            + a04(5)*Hrhoe(i,j,k+4) )*idz_v(k) + Krhoe(i,j,k)
          enddo
       enddo
    endif
    if (is_bc_1pt(3,1)) then
       k=2
       do j=ndy_v,nfy_v
          do i=ndx_v,nfx_v
             Krhou(i,j,k) = ( a13(1)*Hrhou(i,j,k-1)+a13(2)*Hrhou(i,j,k) &
                            + a13(3)*Hrhou(i,j,k+1)+a13(4)*Hrhou(i,j,k+2) &
                            + a13(5)*Hrhou(i,j,k+3) )*idz_v(k) + Krhou(i,j,k)
             Krhov(i,j,k) = ( a13(1)*Hrhov(i,j,k-1)+a13(2)*Hrhov(i,j,k) &
                            + a13(3)*Hrhov(i,j,k+1)+a13(4)*Hrhov(i,j,k+2) &
                            + a13(5)*Hrhov(i,j,k+3) )*idz_v(k) + Krhov(i,j,k)
             Krhow(i,j,k) = ( a13(1)*Hrhow(i,j,k-1)+a13(2)*Hrhow(i,j,k) &
                            + a13(3)*Hrhow(i,j,k+1)+a13(4)*Hrhow(i,j,k+2) &
                            + a13(5)*Hrhow(i,j,k+3) )*idz_v(k) + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a13(1)*Hrhoe(i,j,k-1)+a13(2)*Hrhoe(i,j,k) &
                            + a13(3)*Hrhoe(i,j,k+1)+a13(4)*Hrhoe(i,j,k+2) &
                            + a13(5)*Hrhoe(i,j,k+3) )*idz_v(k) + Krhoe(i,j,k)
          enddo
       enddo
    endif

    do k=ndz_vi,nfz_vi
       do j=ndy_v,nfy_v
          do i=ndx_v,nfx_v
             Krhou(i,j,k) = ( a5(1)*( Hrhou(i,j,k+1)-Hrhou(i,j,k-1)) &
                            + a5(2)*( Hrhou(i,j,k+2)-Hrhou(i,j,k-2)) )*idz_v(k) + Krhou(i,j,k)
             Krhov(i,j,k) = ( a5(1)*( Hrhov(i,j,k+1)-Hrhov(i,j,k-1)) &
                            + a5(2)*( Hrhov(i,j,k+2)-Hrhov(i,j,k-2)) )*idz_v(k) + Krhov(i,j,k)
             Krhow(i,j,k) = ( a5(1)*( Hrhow(i,j,k+1)-Hrhow(i,j,k-1)) &
                            + a5(2)*( Hrhow(i,j,k+2)-Hrhow(i,j,k-2)) )*idz_v(k) + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a5(1)*( Hrhoe(i,j,k+1)-Hrhoe(i,j,k-1)) &
                            + a5(2)*( Hrhoe(i,j,k+2)-Hrhoe(i,j,k-2)) )*idz_v(k) + Krhoe(i,j,k)
          enddo
       enddo
    enddo

    if (is_bc_1pt(3,2)) then
       k=nz-1
       do j=ndy_v,nfy_v
          do i=ndx_v,nfx_v
             Krhou(i,j,k) = ( a31(1)*Hrhou(i,j,k+1)+a31(2)*Hrhou(i,j,k) &
                            + a31(3)*Hrhou(i,j,k-1)+a31(4)*Hrhou(i,j,k-2) &
                            + a31(5)*Hrhou(i,j,k-3) )*idz_v(k) + Krhou(i,j,k)
             Krhov(i,j,k) = ( a31(1)*Hrhov(i,j,k+1)+a31(2)*Hrhov(i,j,k) &
                            + a31(3)*Hrhov(i,j,k-1)+a31(4)*Hrhov(i,j,k-2) &
                            + a31(5)*Hrhov(i,j,k-3) )*idz_v(k) + Krhov(i,j,k)
             Krhow(i,j,k) = ( a31(1)*Hrhow(i,j,k+1)+a31(2)*Hrhow(i,j,k) &
                            + a31(3)*Hrhow(i,j,k-1)+a31(4)*Hrhow(i,j,k-2) &
                            + a31(5)*Hrhow(i,j,k-3) )*idz_v(k) + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a31(1)*Hrhoe(i,j,k+1)+a31(2)*Hrhoe(i,j,k) &
                            + a31(3)*Hrhoe(i,j,k-1)+a31(4)*Hrhoe(i,j,k-2) &
                            + a31(5)*Hrhoe(i,j,k-3) )*idz_v(k) + Krhoe(i,j,k)
          enddo
       enddo
    endif
    if (is_bc_wall(3,2)) then
       k=nz
       do j=ndy_v,nfy_v
          do i=ndx_v,nfx_v
             Krhou(i,j,k) = ( a40(1)*Hrhou(i,j,k)  +a40(2)*Hrhou(i,j,k-1) &
                            + a40(3)*Hrhou(i,j,k-2)+a40(4)*Hrhou(i,j,k-3) &
                            + a40(5)*Hrhou(i,j,k-4) )*idz_v(k) + Krhou(i,j,k)
             Krhov(i,j,k) = ( a40(1)*Hrhov(i,j,k)  +a40(2)*Hrhov(i,j,k-1) &
                            + a40(3)*Hrhov(i,j,k-2)+a40(4)*Hrhov(i,j,k-3) &
                            + a40(5)*Hrhov(i,j,k-4) )*idz_v(k) + Krhov(i,j,k)
             Krhow(i,j,k) = ( a40(1)*Hrhow(i,j,k)  +a40(2)*Hrhow(i,j,k-1) &
                            + a40(3)*Hrhow(i,j,k-2)+a40(4)*Hrhow(i,j,k-3) &
                            + a40(5)*Hrhow(i,j,k-4) )*idz_v(k) + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a40(1)*Hrhoe(i,j,k)  +a40(2)*Hrhoe(i,j,k-1) &
                            + a40(3)*Hrhoe(i,j,k-2)+a40(4)*Hrhoe(i,j,k-3) &
                            + a40(5)*Hrhoe(i,j,k-4) )*idz_v(k) + Krhoe(i,j,k)
          enddo
       enddo
    endif

  end subroutine flux_visc_5pts_c

end submodule smod_flux_visc_5pts_c
