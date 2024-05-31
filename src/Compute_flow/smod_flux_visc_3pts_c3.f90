!=================================================================================
submodule (mod_flux_visc) smod_flux_visc_3pts_c3
!=================================================================================
  !> author: XG
  !> date: January 2024
  !> 3D curvilinear version - compute viscous fluxes with 3-pt stencil
!=================================================================================

contains

  !==============================================================================
  module subroutine flux_visc_3pts_c3
  !==============================================================================
    !> Compute derivatives of viscous fluxes (3-point stencil - order 2)
    !> - Full 3D curvilinear version -
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
                Frhou(i,j,k)= tau11*ksi_x_v(i,j,k)+tau12*ksi_y_v(i,j,k)+tau13*ksi_z_v(i,j,k)
                Frhov(i,j,k)= tau12*ksi_x_v(i,j,k)+tau22*ksi_y_v(i,j,k)+tau23*ksi_z_v(i,j,k)
                Frhow(i,j,k)= tau13*ksi_x_v(i,j,k)+tau23*ksi_y_v(i,j,k)+tau33*ksi_z_v(i,j,k)
                Frhoe(i,j,k)= (uu(i,j,k)*Frhou(i,j,k)+vv(i,j,k)*Frhov(i,j,k)+ww(i,j,k)*Frhow(i,j,k)) &
                     -cok(i,j,k)*(dTx(i,j,k)*ksi_x_v(i,j,k)+dTy(i,j,k)*ksi_y_v(i,j,k)+dTz(i,j,k)*ksi_z_v(i,j,k))

                ! viscous fluxes along eta
                Grhou(i,j,k)= tau11*eta_x_v(i,j,k)+tau12*eta_y_v(i,j,k)+tau13*eta_z_v(i,j,k)
                Grhov(i,j,k)= tau12*eta_x_v(i,j,k)+tau22*eta_y_v(i,j,k)+tau23*eta_z_v(i,j,k)
                Grhow(i,j,k)= tau13*eta_x_v(i,j,k)+tau23*eta_y_v(i,j,k)+tau33*eta_z_v(i,j,k)
                Grhoe(i,j,k)= (uu(i,j,k)*Grhou(i,j,k)+vv(i,j,k)*Grhov(i,j,k)+ww(i,j,k)*Grhow(i,j,k)) &
                     -cok(i,j,k)*(dTx(i,j,k)*eta_x_v(i,j,k)+dTy(i,j,k)*eta_y_v(i,j,k)+dTz(i,j,k)*eta_z_v(i,j,k))

                ! viscous fluxes along phi
                Hrhou(i,j,k)= tau11*phi_x_v(i,j,k)+tau12*phi_y_v(i,j,k)+tau13*phi_z_v(i,j,k)
                Hrhov(i,j,k)= tau12*phi_x_v(i,j,k)+tau22*phi_y_v(i,j,k)+tau23*phi_z_v(i,j,k)
                Hrhow(i,j,k)= tau13*phi_x_v(i,j,k)+tau23*phi_y_v(i,j,k)+tau33*phi_z_v(i,j,k)
                Hrhoe(i,j,k)= (uu(i,j,k)*Hrhou(i,j,k)+vv(i,j,k)*Hrhov(i,j,k)+ww(i,j,k)*Hrhow(i,j,k)) &
                     -cok(i,j,k)*(dTx(i,j,k)*phi_x_v(i,j,k)+dTy(i,j,k)*phi_y_v(i,j,k)+dTz(i,j,k)*phi_z_v(i,j,k))
             enddo
          enddo
       enddo

    else

       ! if (model_RANS.eq.'SA') then
       !    cnu1 = 7.1_wp
       !    prrt = 1.0_wp

       !    do k=ndzt_v,nfzt_v
       !       do j=ndyt_v,nfyt_v
       !          do i=ndxt_v,nfxt_v

       !             ! nutilde wall function
       !             khi  = rho_n(i,j,k)*nutil_n(i,j,k)/visc(i,j,k)
       !             fnu1 = khi**3/(khi**3+cnu1**3)

       !             ! compute Reynolds stress
       !             tauR11  = 2*rho_n(i,j,k)*nutil_n(i,j,k)*fnu1* dux(i,j,k)
       !             tauR22  = 2*rho_n(i,j,k)*nutil_n(i,j,k)*fnu1* dvy(i,j,k)
       !             tauR33  = 2*rho_n(i,j,k)*nutil_n(i,j,k)*fnu1* dwz(i,j,k)
       !             tauR12  =   rho_n(i,j,k)*nutil_n(i,j,k)*fnu1*(duy(i,j,k)+dvx(i,j,k))
       !             tauR13  =   rho_n(i,j,k)*nutil_n(i,j,k)*fnu1*(duz(i,j,k)+dwx(i,j,k))
       !             tauR23  =   rho_n(i,j,k)*nutil_n(i,j,k)*fnu1*(dvz(i,j,k)+dwy(i,j,k))

       !             ! turbulent heat transfer coefficient
       !             cokt = rho_n(i,j,k)*nutil_n(i,j,k)*fnu1*cpcalc_tro(Tmp(i,j,k),rho_n(i,j,k))/prrt

       !             ! compute S_ij
       !             tau11 = dux(i,j,k)
       !             tau22 = dvy(i,j,k)
       !             tau33 = dwz(i,j,k)
       !             tau12 = 0.5_wp*(duy(i,j,k)+dvx(i,j,k))
       !             tau13 = 0.5_wp*(duz(i,j,k)+dwx(i,j,k))
       !             tau23 = 0.5_wp*(dvz(i,j,k)+dwy(i,j,k))
       !             trace = ONE_THIRD*(tau11+tau22+tau33)

       !             ! compute -(tau_ij+tauR_ij)
       !             mu =-2.0_wp*visc(i,j,k)
       !             tau11=mu*(tau11-trace)-tauR11
       !             tau22=mu*(tau22-trace)-tauR22
       !             tau33=mu*(tau33-trace)-tauR33
       !             tau12=mu*tau12-tauR12
       !             tau13=mu*tau13-tauR13
       !             tau23=mu*tau23-tauR23

       !             ! viscous fluxes along ksi
       !             Frhou(i,j,k)= tau11*ksi_x_v(i,j,k)+tau12*ksi_y_v(i,j,k)+tau13*ksi_z_v(i,j,k)
       !             Frhov(i,j,k)= tau12*ksi_x_v(i,j,k)+tau22*ksi_y_v(i,j,k)+tau23*ksi_z_v(i,j,k)
       !             Frhow(i,j,k)= tau13*ksi_x_v(i,j,k)+tau23*ksi_y_v(i,j,k)+tau33*ksi_z_v(i,j,k)
       !             Frhoe(i,j,k)= (uu(i,j,k)*Frhou(i,j,k)+vv(i,j,k)*Frhov(i,j,k)+ww(i,j,k)*Frhow(i,j,k)) &
       !                  -(cok(i,j,k)+cokt)*(dTx(i,j,k)*ksi_x_v(i,j,k)+dTy(i,j,k)*ksi_y_v(i,j,k)+dTz(i,j,k)*ksi_z_v(i,j,k))

       !             ! viscous fluxes along eta
       !             Grhou(i,j,k)= tau11*eta_x_v(i,j,k)+tau12*eta_y_v(i,j,k)+tau13*eta_z_v(i,j,k)
       !             Grhov(i,j,k)= tau12*eta_x_v(i,j,k)+tau22*eta_y_v(i,j,k)+tau23*eta_z_v(i,j,k)
       !             Grhow(i,j,k)= tau13*eta_x_v(i,j,k)+tau23*eta_y_v(i,j,k)+tau33*eta_z_v(i,j,k)
       !             Grhoe(i,j,k)= (uu(i,j,k)*Grhou(i,j,k)+vv(i,j,k)*Grhov(i,j,k)+ww(i,j,k)*Grhow(i,j,k)) &
       !                  -(cok(i,j,k)+cokt)*(dTx(i,j,k)*eta_x_v(i,j,k)+dTy(i,j,k)*eta_y_v(i,j,k)+dTz(i,j,k)*eta_z_v(i,j,k))

       !             ! viscous fluxes along phi
       !             Hrhou(i,j,k)= tau11*phi_x_v(i,j,k)+tau12*phi_y_v(i,j,k)+tau13*phi_z_v(i,j,k)
       !             Hrhov(i,j,k)= tau12*phi_x_v(i,j,k)+tau22*phi_y_v(i,j,k)+tau23*phi_z_v(i,j,k)
       !             Hrhow(i,j,k)= tau13*phi_x_v(i,j,k)+tau23*phi_y_v(i,j,k)+tau33*phi_z_v(i,j,k)
       !             Hrhoe(i,j,k)= (uu(i,j,k)*Hrhou(i,j,k)+vv(i,j,k)*Hrhov(i,j,k)+ww(i,j,k)*Hrhow(i,j,k)) &
       !                  -(cok(i,j,k)+cokt)*(dTx(i,j,k)*phi_x_v(i,j,k)+dTy(i,j,k)*phi_y_v(i,j,k)+dTz(i,j,k)*phi_z_v(i,j,k))
       !          enddo
       !       enddo
       !    enddo
       ! endif

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
                   Frhou(i,j,k)= tau11*ksi_x_v(i,j,k)+tau12*ksi_y_v(i,j,k)+tau13*ksi_z_v(i,j,k)
                   Frhov(i,j,k)= tau12*ksi_x_v(i,j,k)+tau22*ksi_y_v(i,j,k)+tau23*ksi_z_v(i,j,k)
                   Frhow(i,j,k)= tau13*ksi_x_v(i,j,k)+tau23*ksi_y_v(i,j,k)+tau33*ksi_z_v(i,j,k)
                   Frhoe(i,j,k)= (uu(i,j,k)*Frhou(i,j,k)+vv(i,j,k)*Frhov(i,j,k)+ww(i,j,k)*Frhow(i,j,k)) &
                        -(cok(i,j,k)+cokt(i,j,k))*(dTx(i,j,k)*ksi_x_v(i,j,k)+dTy(i,j,k)*ksi_y_v(i,j,k)+dTz(i,j,k)*ksi_z_v(i,j,k))

                   ! viscous fluxes along eta
                   Grhou(i,j,k)= tau11*eta_x_v(i,j,k)+tau12*eta_y_v(i,j,k)+tau13*eta_z_v(i,j,k)
                   Grhov(i,j,k)= tau12*eta_x_v(i,j,k)+tau22*eta_y_v(i,j,k)+tau23*eta_z_v(i,j,k)
                   Grhow(i,j,k)= tau13*eta_x_v(i,j,k)+tau23*eta_y_v(i,j,k)+tau33*eta_z_v(i,j,k)
                   Grhoe(i,j,k)= (uu(i,j,k)*Grhou(i,j,k)+vv(i,j,k)*Grhov(i,j,k)+ww(i,j,k)*Grhow(i,j,k)) &
                        -(cok(i,j,k)+cokt(i,j,k))*(dTx(i,j,k)*eta_x_v(i,j,k)+dTy(i,j,k)*eta_y_v(i,j,k)+dTz(i,j,k)*eta_z_v(i,j,k))

                   ! viscous fluxes along phi
                   Hrhou(i,j,k)= tau11*phi_x_v(i,j,k)+tau12*phi_y_v(i,j,k)+tau13*phi_z_v(i,j,k)
                   Hrhov(i,j,k)= tau12*phi_x_v(i,j,k)+tau22*phi_y_v(i,j,k)+tau23*phi_z_v(i,j,k)
                   Hrhow(i,j,k)= tau13*phi_x_v(i,j,k)+tau23*phi_y_v(i,j,k)+tau33*phi_z_v(i,j,k)
                   Hrhoe(i,j,k)= (uu(i,j,k)*Hrhou(i,j,k)+vv(i,j,k)*Hrhov(i,j,k)+ww(i,j,k)*Hrhow(i,j,k)) &
                        -(cok(i,j,k)+cokt(i,j,k))*(dTx(i,j,k)*phi_x_v(i,j,k)+dTy(i,j,k)*phi_y_v(i,j,k)+dTz(i,j,k)*phi_z_v(i,j,k))
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
             Krhou(i,j,k) = ( a02(1)*Frhou(i  ,j,k)+a02(2)*Frhou(i+1,j,k) &
                            + a02(3)*Frhou(i+2,j,k) )*ijacob3_v(i,j,k) + Krhou(i,j,k)
             Krhov(i,j,k) = ( a02(1)*Frhov(i  ,j,k)+a02(2)*Frhov(i+1,j,k) &
                            + a02(3)*Frhov(i+2,j,k) )*ijacob3_v(i,j,k) + Krhov(i,j,k)
             Krhow(i,j,k) = ( a02(1)*Frhow(i  ,j,k)+a02(2)*Frhow(i+1,j,k) &
                            + a02(3)*Frhow(i+2,j,k) )*ijacob3_v(i,j,k) + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a02(1)*Frhoe(i  ,j,k)+a02(2)*Frhoe(i+1,j,k) &
                            + a02(3)*Frhoe(i+2,j,k) )*ijacob3_v(i,j,k) + Krhoe(i,j,k)
          enddo
       enddo
    endif

    do k=ndz_v,nfz_v
       do j=ndy_v,nfy_v
          do i=ndx_vi,nfx_vi
             Krhou(i,j,k) = ( a3(1)*( Frhou(i+1,j,k)-Frhou(i-1,j,k) ) )*ijacob3_v(i,j,k) + Krhou(i,j,k)
             Krhov(i,j,k) = ( a3(1)*( Frhov(i+1,j,k)-Frhov(i-1,j,k) ) )*ijacob3_v(i,j,k) + Krhov(i,j,k)
             Krhow(i,j,k) = ( a3(1)*( Frhow(i+1,j,k)-Frhow(i-1,j,k) ) )*ijacob3_v(i,j,k) + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a3(1)*( Frhoe(i+1,j,k)-Frhoe(i-1,j,k) ) )*ijacob3_v(i,j,k) + Krhoe(i,j,k)
          enddo
       enddo
    enddo

    if (is_bc_wall(1,2)) then
       i=nx
       do k=ndz_v,nfz_v
          do j=ndy_v,nfy_v
             Krhou(i,j,k) = ( a20(1)*Frhou(i,j,k)  +a20(2)*Frhou(i-1,j,k) &
                            + a20(3)*Frhou(i-2,j,k) )*ijacob3_v(i,j,k) + Krhou(i,j,k)
             Krhov(i,j,k) = ( a20(1)*Frhov(i,j,k)  +a20(2)*Frhov(i-1,j,k) &
                            + a20(3)*Frhov(i-2,j,k) )*ijacob3_v(i,j,k) + Krhov(i,j,k)
             Krhow(i,j,k) = ( a20(1)*Frhow(i,j,k)  +a20(2)*Frhow(i-1,j,k) &
                            + a20(3)*Frhow(i-2,j,k) )*ijacob3_v(i,j,k) + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a20(1)*Frhoe(i,j,k)  +a20(2)*Frhoe(i-1,j,k) &
                            + a20(3)*Frhoe(i-2,j,k) )*ijacob3_v(i,j,k) + Krhoe(i,j,k)
          enddo
       enddo
    endif

    ! Compute derivatives of viscous fluxes along eta
    ! ===============================================
    if (is_bc_wall(2,1)) then
       j=1
       do k=ndz_v,nfz_v
          do i=ndx_v,nfx_v
             Krhou(i,j,k) = ( a02(1)*Grhou(i,j  ,k)+a02(2)*Grhou(i,j+1,k) &
                            + a02(3)*Grhou(i,j+2,k) )*ijacob3_v(i,j,k) + Krhou(i,j,k)
             Krhov(i,j,k) = ( a02(1)*Grhov(i,j  ,k)+a02(2)*Grhov(i,j+1,k) &
                            + a02(3)*Grhov(i,j+2,k) )*ijacob3_v(i,j,k) + Krhov(i,j,k)
             Krhow(i,j,k) = ( a02(1)*Grhow(i,j  ,k)+a02(2)*Grhow(i,j+1,k) &
                            + a02(3)*Grhow(i,j+2,k) )*ijacob3_v(i,j,k) + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a02(1)*Grhoe(i,j  ,k)+a02(2)*Grhoe(i,j+1,k) &
                            + a02(3)*Grhoe(i,j+2,k) )*ijacob3_v(i,j,k) + Krhoe(i,j,k)
          enddo
       enddo
    endif

    do k=ndz_v,nfz_v
       do j=ndy_vi,nfy_vi
          do i=ndx_v,nfx_v
             Krhou(i,j,k) = ( a3(1)*( Grhou(i,j+1,k)-Grhou(i,j-1,k)) )*ijacob3_v(i,j,k) + Krhou(i,j,k)
             Krhov(i,j,k) = ( a3(1)*( Grhov(i,j+1,k)-Grhov(i,j-1,k)) )*ijacob3_v(i,j,k) + Krhov(i,j,k)
             Krhow(i,j,k) = ( a3(1)*( Grhow(i,j+1,k)-Grhow(i,j-1,k)) )*ijacob3_v(i,j,k) + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a3(1)*( Grhoe(i,j+1,k)-Grhoe(i,j-1,k)) )*ijacob3_v(i,j,k) + Krhoe(i,j,k)
          enddo
       enddo
    enddo

    if (is_bc_wall(2,2)) then
       j=ny
       do k=ndz_v,nfz_v
          do i=ndx_v,nfx_v
             Krhou(i,j,k) = ( a20(1)*Grhou(i,j,k)  +a20(2)*Grhou(i,j-1,k) &
                            + a20(3)*Grhou(i,j-2,k) )*ijacob3_v(i,j,k) + Krhou(i,j,k)
             Krhov(i,j,k) = ( a20(1)*Grhov(i,j,k)  +a20(2)*Grhov(i,j-1,k) &
                            + a20(3)*Grhov(i,j-2,k) )*ijacob3_v(i,j,k) + Krhov(i,j,k)
             Krhow(i,j,k) = ( a20(1)*Grhow(i,j,k)  +a20(2)*Grhow(i,j-1,k) &
                            + a20(3)*Grhow(i,j-2,k) )*ijacob3_v(i,j,k) + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a20(1)*Grhoe(i,j,k)  +a20(2)*Grhoe(i,j-1,k) &
                            + a20(3)*Grhoe(i,j-2,k) )*ijacob3_v(i,j,k) + Krhoe(i,j,k)
          enddo
       enddo
    endif

    ! Compute derivatives of viscous fluxes along phi
    ! ===============================================
    if (is_bc_wall(3,1)) then
       k=1
       do j=ndy_v,nfy_v
          do i=ndx_v,nfx_v
             Krhou(i,j,k) = ( a02(1)*Hrhou(i,j,k  )+a02(2)*Hrhou(i,j,k+1) &
                            + a02(3)*Hrhou(i,j,k+2) )*ijacob3_v(i,j,k) + Krhou(i,j,k)
             Krhov(i,j,k) = ( a02(1)*Hrhov(i,j,k  )+a02(2)*Hrhov(i,j,k+1) &
                            + a02(3)*Hrhov(i,j,k+2) )*ijacob3_v(i,j,k) + Krhov(i,j,k)
             Krhow(i,j,k) = ( a02(1)*Hrhow(i,j,k  )+a02(2)*Hrhow(i,j,k+1) &
                            + a02(3)*Hrhow(i,j,k+2) )*ijacob3_v(i,j,k) + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a02(1)*Hrhoe(i,j,k  )+a02(2)*Hrhoe(i,j,k+1) &
                            + a02(3)*Hrhoe(i,j,k+2) )*ijacob3_v(i,j,k) + Krhoe(i,j,k)
          enddo
       enddo
    endif

    do k=ndz_vi,nfz_vi
       do j=ndy_v,nfy_v
          do i=ndx_v,nfx_v
             Krhou(i,j,k) = ( a3(1)*( Hrhou(i,j,k+1)-Hrhou(i,j,k-1)) )*ijacob3_v(i,j,k) + Krhou(i,j,k)
             Krhov(i,j,k) = ( a3(1)*( Hrhov(i,j,k+1)-Hrhov(i,j,k-1)) )*ijacob3_v(i,j,k) + Krhov(i,j,k)
             Krhow(i,j,k) = ( a3(1)*( Hrhow(i,j,k+1)-Hrhow(i,j,k-1)) )*ijacob3_v(i,j,k) + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a3(1)*( Hrhoe(i,j,k+1)-Hrhoe(i,j,k-1)) )*ijacob3_v(i,j,k) + Krhoe(i,j,k)
          enddo
       enddo
    enddo

    if (is_bc_wall(3,2)) then
       k=nz
       do j=ndy_v,nfy_v
          do i=ndx_v,nfx_v
             Krhou(i,j,k) = ( a20(1)*Hrhou(i,j,k)  +a20(2)*Hrhou(i,j,k-1) &
                            + a20(3)*Hrhou(i,j,k-2) )*ijacob3_v(i,j,k) + Krhou(i,j,k)
             Krhov(i,j,k) = ( a20(1)*Hrhov(i,j,k)  +a20(2)*Hrhov(i,j,k-1) &
                            + a20(3)*Hrhov(i,j,k-2) )*ijacob3_v(i,j,k) + Krhov(i,j,k)
             Krhow(i,j,k) = ( a20(1)*Hrhow(i,j,k)  +a20(2)*Hrhow(i,j,k-1) &
                            + a20(3)*Hrhow(i,j,k-2) )*ijacob3_v(i,j,k) + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a20(1)*Hrhoe(i,j,k)  +a20(2)*Hrhoe(i,j,k-1) &
                            + a20(3)*Hrhoe(i,j,k-2) )*ijacob3_v(i,j,k) + Krhoe(i,j,k)
          enddo
       enddo
    endif

  end subroutine flux_visc_3pts_c3

end submodule smod_flux_visc_3pts_c3
