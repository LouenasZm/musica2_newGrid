!==============================================================================
subroutine flux_visc_5pts_SM_wall
    !==============================================================================
      !> Compute derivatives of viscous fluxes (5-point stencil - order 4)
      !> with Smagorinsky model
      !> - Cartesian version -
    !==============================================================================
      use mod_coeff_deriv
      use mod_flow
      use mod_interface
      use mod_constant
      use mod_mpi
      use mod_wall_model
      use warnstop
      implicit none
      ! ---------------------------------------------------------------------------
      integer :: i,j,k
      real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2) :: dTx2,dTy2,dTz2,visco_sgs,k_sgs
      real(wp) :: deltac
      real(wp) :: tau11,tau22,tau33,tau12,tau13,tau23,trace,mu
      ! ---------------------------------------------------------------------------
      real(wp) :: Aplus,utau_mean,visc_mean,rho_mean
      real(wp), dimension(ny1:ny2) :: VD,VD2
    
      ! for a regular grid
      deltac = (idx(1)*idy(1)*idz(1))**(-1./3.)
    
      ! Initialization
      ! ==============
    
      dTx2 = 0.
      dTy2 = 0.
      dTz2 = 0.
    
      S11 = 0.
      S12 = 0.
      S13 = 0.
      S22 = 0.
      S23 = 0.
      S33 = 0.
    
      ! Compute Sij and grad(T)
      ! =======================
      ! Lower wall
      do k=ndzt_v,nfzt_v
         do j=ndyt_v,nfyt_v
            if(j <= wm_ind .or. j >=(nfyt_v-wm_ind))then
            do i=ndxt_v,nfxt_v
               ! compute S_ij
               S11(i,j,k) = dux(i,j,k)
               S22(i,j,k) = dvy(i,j,k)
               S33(i,j,k) = dwz(i,j,k)
               S12(i,j,k) = 0.5_wp*(duy(i,j,k)+dvx(i,j,k))
               S13(i,j,k) = 0.5_wp*(duz(i,j,k)+dwx(i,j,k))
               S23(i,j,k) = 0.5_wp*(dvz(i,j,k)+dwy(i,j,k))
            enddo
        endif
         enddo
      enddo

      ! Filtering of Sij
      ! ================
      if (is_filt_Sij) call  filtre_Sij_wall
    
      ! Turbulent viscosity
      ! ===================
      if (is_filt_nu) then
        do k=ndzt_v,nfzt_v
            do j=ndyt_v,nfyt_v
            if(j <= wm_ind .or. j >=(nfyt_v-wm_ind))then
               do i=ndxt_v,nfxt_v
                 !Calculation of turbulent viscosity
                  visco_sgs(i,j,k) = rho(i,j,k)*((Cs_SM*deltac)**2)*sqrt(2.*(S11f(i,j,k)**2+S22f(i,j,k)**2+S33f(i,j,k)**2 &
                                                  + 2.*(S12f(i,j,k)**2 +S13f(i,j,k)**2+S23f(i,j,k)**2)))
                  !Calculation of Taukk/2 = energie cin de sous maille
                  k_sgs(i,j,k) = Ci_SM * rho(i,j,k) * deltac**2 * (2.*(S11f(i,j,k)**2+S22f(i,j,k)**2+S33f(i,j,k)**2+2.*(S12f(i,j,k)**2+S13f(i,j,k)**2+S23f(i,j,k)**2)))
               enddo
            endif
            enddo
           
        enddo
      else
         do k=ndzt_v,nfzt_v
            ! Lower wall:
            do j=ndyt_v,nfyt_v
            if(j <= wm_ind .or. j >=(nfyt_v-wm_ind))then
               do i=ndxt_v,nfxt_v
                 !Calculation of turbulent viscosity
                  visco_sgs(i,j,k) = rho(i,j,k)*((Cs_SM*deltac)**2)*sqrt(2.*(S11(i,j,k)**2+S22(i,j,k)**2+S33(i,j,k)**2 &
                                                  + 2.*(S12(i,j,k)**2 +S13(i,j,k)**2+S23(i,j,k)**2)))
                  !Calculation of Taukk/2 = energie cin de sous maille
                  k_sgs(i,j,k) = Ci_SM * rho(i,j,k) * deltac**2 * (2.*(S11(i,j,k)**2+S22(i,j,k)**2+S33(i,j,k)**2+2.*(S12(i,j,k)**2+S13(i,j,k)**2+S23(i,j,k)**2)))
               enddo
            endif
            enddo
         enddo
      end if
    
      ! Van driest damping function /!\ Just for CHAN
      ! ===========================
      Aplus = 26.0_wp
      VD = 1.0_wp; VD2 = 1.0_wp
      if (is_bc_wall(2,1)) then
         utau_mean = SUM(utau_jmin(1:nx,1:nz))/(nz*nx)
         rho_mean = SUM(rho(1:nx,1,1:nz))/(nz*nx)
         visc_mean = SUM(visc(1:nx,1,1:nz))/(nz*nx)
         ! Lower wall:
         do j=wm_ind,nfyt_v-wm_ind
           VD(j) = 1.0_wp - exp(-(abs(y(j)-y(1))*utau_mean*rho_mean/visc_mean)/Aplus)
         enddo
         ! Upper wall: 
         do j=nfyt_v-wm_ind,nfyt_v
            VD(j) = 1.0_wp - exp(-(abs(y(j)-y(1))*utau_mean*rho_mean/visc_mean)/Aplus)
          enddo
      endif
      if (is_bc_wall(2,2)) then
         utau_mean = SUM(utau_jmax(1:nx,1:nz))/(nz*nx)
         rho_mean = SUM(rho(1:nx,ny,1:nz))/(nz*nx)
         visc_mean = SUM(visc(1:nx,ny,1:nz))/(nz*nx)
         ! Lower wall
         do j=wm_ind,nfyt_v-wm_ind-1
           VD2(j) = 1.0_wp - exp(-(abs(y(j)-y(ny))*utau_mean*rho_mean/visc_mean)/Aplus)
         enddo
         ! Upper wall 
         do j=nfyt_v-wm_ind,nfyt_v
            VD2(j) = 1.0_wp - exp(-(abs(y(j)-y(ny))*utau_mean*rho_mean/visc_mean)/Aplus)
          enddo
      endif
    
      if ((is_bc_wall(2,1)).or.(is_bc_wall(2,2))) then
         VD = VD*VD2
         do k=ndzt_v,nfzt_v
            ! Lower wall
            do j=wm_ind,nfyt_v-wm_ind-1
               do i=ndxt_v,nfxt_v
                  ! Van driest damping at the wall
                  visco_sgs(i,j,k) = visco_sgs(i,j,k)*VD(j)
               enddo
            enddo
            ! Upper wall
            do j=nfyt_v-wm_ind,nfyt_v
                do i=ndxt_v,nfxt_v
                   ! Van driest damping at the wall
                   visco_sgs(i,j,k) = visco_sgs(i,j,k)*VD(j)
                enddo
             enddo
         enddo
      endif
    
      if (is_filt_Sij) then
         do k=ndzt_v,nfzt_v
            do j=ndyt_v,nfyt_v
              if (j <= wm_ind .or. j >= (nfyt_v-wm_ind))then 
               do i=ndxt_v,nfxt_v
    
                  ! Calculation of viscoux fluxes
                  ! =============================
                  ! viscous fluxes along x
                  Frhou(i,j,k) = - visc(i,j,k) * 2 *( S11(i,j,k) - (S11(i,j,k) + S22(i,j,k) + S33(i,j,k))/3.)
                  Frhov(i,j,k) = - visc(i,j,k) * 2 * S12(i,j,k)
                  Frhow(i,j,k) = - visc(i,j,k) * 2 * S23(i,j,k)
                  Frhoe(i,j,k) = (uu(i,j,k)*Frhou(i,j,k) + vv(i,j,k)*Frhov(i,j,k) + ww(i,j,k)*Frhow(i,j,k)) &
                                  - cok(i,j,k)*dTx(i,j,k)
    
                  ! viscous fluxes along y
                  Grhou(i,j,k) = - visc(i,j,k) * 2.* S12(i,j,k)
                  Grhov(i,j,k) = - visc(i,j,k) * 2.*( S22(i,j,k) - (S11(i,j,k) + S22(i,j,k) + S33(i,j,k))/3.)
                  Grhow(i,j,k) = - visc(i,j,k) * 2.* S23(i,j,k)
                  Grhoe(i,j,k) = (uu(i,j,k)*Grhou(i,j,k) + vv(i,j,k)*Grhov(i,j,k) + ww(i,j,k)*Grhow(i,j,k)) &
                                  - cok(i,j,k)*dTy(i,j,k)
    
                  ! viscous fluxes along z
                  Hrhou(i,j,k) = - visc(i,j,k) * 2.* S13(i,j,k)
                  Hrhov(i,j,k) = - visc(i,j,k) * 2.* S23(i,j,k)
                  Hrhow(i,j,k) = - visc(i,j,k) * 2.*( S33(i,j,k) - (S11(i,j,k) + S22(i,j,k) + S33(i,j,k))/3.)
                  Hrhoe(i,j,k) = (uu(i,j,k)*Hrhou(i,j,k) + vv(i,j,k)*Hrhov(i,j,k) + ww(i,j,k)*Hrhow(i,j,k)) &
                                - cok(i,j,k)*dTz(i,j,k)
    
                  ! Calculation of filtered viscoux fluxes
                  ! ======================================
                  ! filtered viscous fluxes along x
                  Frhou(i,j,k) = Frhou(i,j,k) - visco_sgs(i,j,k) * 2 *( S11f(i,j,k) - (S11f(i,j,k)+S22f(i,j,k)+S33f(i,j,k))/3.) + 2*k_sgs(i,j,k)/3
                  Frhov(i,j,k) = Frhov(i,j,k) - visco_sgs(i,j,k) * 2 * S12f(i,j,k)
                  Frhow(i,j,k) = Frhow(i,j,k) - visco_sgs(i,j,k) * 2 * S23f(i,j,k)
                  Frhoe(i,j,k) = Frhoe(i,j,k) - visco_sgs(i,j,k) * 2.*( rhou(i,j,k)*(S11f(i,j,k) - (S11f(i,j,k)+S22f(i,j,k)+S33f(i,j,k))/3.) &
                                                                 + rhov(i,j,k)*S12f(i,j,k) + rhow(i,j,k)*S13f(i,j,k) ) /rho(i,j,k) &
                                                                 + 2*k_sgs(i,j,k)/3 * rhou(i,j,k)/rho(i,j,k)
    
                  ! filtered viscous fluxes along y
                  Grhou(i,j,k) = Grhou(i,j,k) - visco_sgs(i,j,k) * 2.*S12f(i,j,k)
                  Grhov(i,j,k) = Grhov(i,j,k) - visco_sgs(i,j,k) * 2.*(S22f(i,j,k) - (S11f(i,j,k)+S22f(i,j,k)+S33f(i,j,k))/3.) + 2*k_sgs(i,j,k)/3
                  Grhow(i,j,k) = Grhow(i,j,k) - visco_sgs(i,j,k) * 2.*S23f(i,j,k)
                  Grhoe(i,j,k) = Grhoe(i,j,k) - visco_sgs(i,j,k) * 2.*( rhou(i,j,k)*S12f(i,j,k) + rhov(i,j,k)*(S22f(i,j,k) - (S11f(i,j,k)+S22f(i,j,k) &
                                                                 +S33f(i,j,k))/3.) + rhow(i,j,k)*S23f(i,j,k) ) /rho(i,j,k) &
                                                                 + 2*k_sgs(i,j,k)/3 * rhov(i,j,k)/rho(i,j,k)
    
                  ! filtered viscous fluxes along z
                  Hrhou(i,j,k) = Hrhou(i,j,k) - visco_sgs(i,j,k) * 2.*S13f(i,j,k)
                  Hrhov(i,j,k) = Hrhov(i,j,k) - visco_sgs(i,j,k) * 2.*S23f(i,j,k)
                  Hrhow(i,j,k) = Hrhow(i,j,k) - visco_sgs(i,j,k) * 2.*(S33f(i,j,k) - (S11f(i,j,k)+S22f(i,j,k)+S33f(i,j,k))/3.) + 2*k_sgs(i,j,k)/3
                  Hrhoe(i,j,k) = Hrhoe(i,j,k) - visco_sgs(i,j,k) * 2.*( rhou(i,j,k)*S13f(i,j,k) + rhov(i,j,k)*S23f(i,j,k) + rhow(i,j,k)*(S33f(i,j,k) &
                                                                 - (S11f(i,j,k)+S22f(i,j,k)+S33f(i,j,k))/3.) ) /rho(i,j,k) &
                                                                 + 2*k_sgs(i,j,k)/3 * rhow(i,j,k)/rho(i,j,k)
    
               enddo
              endif
            enddo
         enddo
      else
         do k=ndzt_v,nfzt_v
            do j=ndyt_v,nfyt_v
              if (j <= wm_ind .or. j >= (nfyt_v-wm_ind))then 
               do i=ndxt_v,nfxt_v
    
                  ! Calculation of viscoux fluxes
                  ! =============================
                  ! viscous fluxes along x
                  Frhou(i,j,k) = - (visc(i,j,k) + visco_sgs(i,j,k)) * 2 *( S11(i,j,k) - (S11(i,j,k) + S22(i,j,k) + S33(i,j,k))/3.)
                  Frhov(i,j,k) = - (visc(i,j,k) + visco_sgs(i,j,k)) * 2 * S12(i,j,k)
                  Frhow(i,j,k) = - (visc(i,j,k) + visco_sgs(i,j,k)) * 2 * S23(i,j,k)
                  Frhoe(i,j,k) = (uu(i,j,k)*Frhou(i,j,k) + vv(i,j,k)*Frhov(i,j,k) + ww(i,j,k)*Frhow(i,j,k)) &
                                  - cok(i,j,k)*dTx(i,j,k)
    
                  ! viscous fluxes along y
                  Grhou(i,j,k) = - (visc(i,j,k) + visco_sgs(i,j,k)) * 2.* S12(i,j,k)
                  Grhov(i,j,k) = - (visc(i,j,k) + visco_sgs(i,j,k)) * 2.*( S22(i,j,k) - (S11(i,j,k) + S22(i,j,k) + S33(i,j,k))/3.)
                  Grhow(i,j,k) = - (visc(i,j,k) + visco_sgs(i,j,k)) * 2.* S23(i,j,k)
                  Grhoe(i,j,k) = (uu(i,j,k)*Grhou(i,j,k) + vv(i,j,k)*Grhov(i,j,k) + ww(i,j,k)*Grhow(i,j,k)) &
                                  - cok(i,j,k)*dTy(i,j,k)
    
                  ! viscous fluxes along z
                  Hrhou(i,j,k) = - (visc(i,j,k) + visco_sgs(i,j,k)) * 2.* S13(i,j,k)
                  Hrhov(i,j,k) = - (visc(i,j,k) + visco_sgs(i,j,k)) * 2.* S23(i,j,k)
                  Hrhow(i,j,k) = - (visc(i,j,k) + visco_sgs(i,j,k)) * 2.*( S33(i,j,k) - (S11(i,j,k) + S22(i,j,k) + S33(i,j,k))/3.)
                  Hrhoe(i,j,k) = (uu(i,j,k)*Hrhou(i,j,k) + vv(i,j,k)*Hrhov(i,j,k) + ww(i,j,k)*Hrhow(i,j,k)) &
                                - cok(i,j,k)*dTz(i,j,k)
    
               end do
            endif
            end do
         end do
      end if
    
      ! Compute viscous fluxes
        ! ======================
      do k=ndzt_v,nfzt_v
        do j=ndyt_v,nfyt_v
          if (j >= wm_ind .or. j <= (nfyt_v-wm_ind))then 
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
              Grhou(i,j,k)= tau12
              Grhov(i,j,k)= tau22
              Grhow(i,j,k)= tau23
              Grhoe(i,j,k)= (uu(i,j,k)*Grhou(i,j,k)+vv(i,j,k)*Grhov(i,j,k)+ww(i,j,k)*Grhow(i,j,k)) &
                   -cok(i,j,k)*dTy(i,j,k)
    
              ! viscous fluxes along z
              Hrhou(i,j,k)= tau13
              Hrhov(i,j,k)= tau23
              Hrhow(i,j,k)= tau33
              Hrhoe(i,j,k)= (uu(i,j,k)*Hrhou(i,j,k)+vv(i,j,k)*Hrhov(i,j,k)+ww(i,j,k)*Hrhow(i,j,k)) &
                   -cok(i,j,k)*dTz(i,j,k)
           enddo
          endif
        enddo
     enddo
    
    
      ! Wall modeling
      ! =============
      if (is_wall_model) then
         if (is_bc_wall(2,1)) call bc_wm_jmin
         if (is_bc_wall(2,2)) call bc_wm_jmax
      endif
    
    
      ! Compute derivatives of viscous fluxes along x
      ! =============================================
      if (is_bc_wall(1,1)) then
         i=1
         do k=ndz_v,nfz_v
            do j=ndy_v,nfy_v
               Krhou(i,j,k) = ( a04(1)*Frhou(i  ,j,k)+a04(2)*Frhou(i+1,j,k) &
                              + a04(3)*Frhou(i+2,j,k)+a04(4)*Frhou(i+3,j,k) &
                              + a04(5)*Frhou(i+4,j,k) )*idx_v(i) + Krhou(i,j,k)
               Krhov(i,j,k) = ( a04(1)*Frhov(i  ,j,k)+a04(2)*Frhov(i+1,j,k) &
                              + a04(3)*Frhov(i+2,j,k)+a04(4)*Frhov(i+3,j,k) &
                              + a04(5)*Frhov(i+4,j,k) )*idx_v(i) + Krhov(i,j,k)
               Krhow(i,j,k) = ( a04(1)*Frhow(i  ,j,k)+a04(2)*Frhow(i+1,j,k) &
                              + a04(3)*Frhow(i+2,j,k)+a04(4)*Frhow(i+3,j,k) &
                              + a04(5)*Frhow(i+4,j,k) )*idx_v(i) + Krhow(i,j,k)
               Krhoe(i,j,k) = ( a04(1)*Frhoe(i  ,j,k)+a04(2)*Frhoe(i+1,j,k) &
                              + a04(3)*Frhoe(i+2,j,k)+a04(4)*Frhoe(i+3,j,k) &
                              + a04(5)*Frhoe(i+4,j,k) )*idx_v(i) + Krhoe(i,j,k)
            enddo
         enddo
      endif
      if (is_bc_1pt(1,1)) then
         i=2
         do k=ndz_v,nfz_v
            do j=ndy_v,wm_ind-1
               Krhou(i,j,k) = ( a13(1)*Frhou(i-1,j,k)+a13(2)*Frhou(i,j,k) &
                              + a13(3)*Frhou(i+1,j,k)+a13(4)*Frhou(i+2,j,k) &
                              + a13(5)*Frhou(i+3,j,k) )*idx_v(i) + Krhou(i,j,k)
               Krhov(i,j,k) = ( a13(1)*Frhov(i-1,j,k)+a13(2)*Frhov(i,j,k) &
                              + a13(3)*Frhov(i+1,j,k)+a13(4)*Frhov(i+2,j,k) &
                              + a13(5)*Frhov(i+3,j,k) )*idx_v(i) + Krhov(i,j,k)
               Krhow(i,j,k) = ( a13(1)*Frhow(i-1,j,k)+a13(2)*Frhow(i,j,k) &
                              + a13(3)*Frhow(i+1,j,k)+a13(4)*Frhow(i+2,j,k) &
                              + a13(5)*Frhow(i+3,j,k) )*idx_v(i) + Krhow(i,j,k)
               Krhoe(i,j,k) = ( a13(1)*Frhoe(i-1,j,k)+a13(2)*Frhoe(i,j,k) &
                              + a13(3)*Frhoe(i+1,j,k)+a13(4)*Frhoe(i+2,j,k) &
                              + a13(5)*Frhoe(i+3,j,k) )*idx_v(i) + Krhoe(i,j,k)
            enddo
         enddo
      endif
    
      do k=ndz_v,nfz_v
         do j=ndy_v,wm_ind-1
            do i=ndx_vi,nfx_vi
               Krhou(i,j,k) = ( a5(1)*( Frhou(i+1,j,k)-Frhou(i-1,j,k) ) &
                              + a5(2)*( Frhou(i+2,j,k)-Frhou(i-2,j,k) ) )*idx_v(i) + Krhou(i,j,k)
               Krhov(i,j,k) = ( a5(1)*( Frhov(i+1,j,k)-Frhov(i-1,j,k) ) &
                              + a5(2)*( Frhov(i+2,j,k)-Frhov(i-2,j,k) ) )*idx_v(i) + Krhov(i,j,k)
               Krhow(i,j,k) = ( a5(1)*( Frhow(i+1,j,k)-Frhow(i-1,j,k) ) &
                              + a5(2)*( Frhow(i+2,j,k)-Frhow(i-2,j,k) ) )*idx_v(i) + Krhow(i,j,k)
               Krhoe(i,j,k) = ( a5(1)*( Frhoe(i+1,j,k)-Frhoe(i-1,j,k) ) &
                              + a5(2)*( Frhoe(i+2,j,k)-Frhoe(i-2,j,k) ) )*idx_v(i) + Krhoe(i,j,k)
            enddo
         enddo
      enddo
    
      if (is_bc_1pt(1,2)) then
         i=nx-1
         do k=ndz_v,nfz_v
            do j=ndy_v,wm_ind-1
               Krhou(i,j,k) = ( a31(1)*Frhou(i+1,j,k)+a31(2)*Frhou(i,j,k) &
                              + a31(3)*Frhou(i-1,j,k)+a31(4)*Frhou(i-2,j,k) &
                              + a31(5)*Frhou(i-3,j,k) )*idx_v(i) + Krhou(i,j,k)
               Krhov(i,j,k) = ( a31(1)*Frhov(i+1,j,k)+a31(2)*Frhov(i,j,k) &
                              + a31(3)*Frhov(i-1,j,k)+a31(4)*Frhov(i-2,j,k) &
                              + a31(5)*Frhov(i-3,j,k) )*idx_v(i) + Krhov(i,j,k)
               Krhow(i,j,k) = ( a31(1)*Frhow(i+1,j,k)+a31(2)*Frhow(i,j,k) &
                              + a31(3)*Frhow(i-1,j,k)+a31(4)*Frhow(i-2,j,k) &
                              + a31(5)*Frhow(i-3,j,k) )*idx_v(i) + Krhow(i,j,k)
               Krhoe(i,j,k) = ( a31(1)*Frhoe(i+1,j,k)+a31(2)*Frhoe(i,j,k) &
                              + a31(3)*Frhoe(i-1,j,k)+a31(4)*Frhoe(i-2,j,k) &
                              + a31(5)*Frhoe(i-3,j,k) )*idx_v(i) + Krhoe(i,j,k)
            enddo
         enddo
      endif
      if (is_bc_wall(1,2)) then
         i=nx
         do k=ndz_v,nfz_v
            do j=ndy_v,wm_ind-1
               Krhou(i,j,k) = ( a40(1)*Frhou(i,j,k)  +a40(2)*Frhou(i-1,j,k) &
                              + a40(3)*Frhou(i-2,j,k)+a40(4)*Frhou(i-3,j,k) &
                              + a40(5)*Frhou(i-4,j,k) )*idx_v(i) + Krhou(i,j,k)
               Krhov(i,j,k) = ( a40(1)*Frhov(i,j,k)  +a40(2)*Frhov(i-1,j,k) &
                              + a40(3)*Frhov(i-2,j,k)+a40(4)*Frhov(i-3,j,k) &
                              + a40(5)*Frhov(i-4,j,k) )*idx_v(i) + Krhov(i,j,k)
               Krhow(i,j,k) = ( a40(1)*Frhow(i,j,k)  +a40(2)*Frhow(i-1,j,k) &
                              + a40(3)*Frhow(i-2,j,k)+a40(4)*Frhow(i-3,j,k) &
                              + a40(5)*Frhow(i-4,j,k) )*idx_v(i) + Krhow(i,j,k)
               Krhoe(i,j,k) = ( a40(1)*Frhoe(i,j,k)  +a40(2)*Frhoe(i-1,j,k) &
                              + a40(3)*Frhoe(i-2,j,k)+a40(4)*Frhoe(i-3,j,k) &
                              + a40(5)*Frhoe(i-4,j,k) )*idx_v(i) + Krhoe(i,j,k)
            enddo
         enddo
      endif
    
      ! Compute derivatives of viscous fluxes along y
      ! =============================================
      if (is_bc_wall(2,1)) then
         j=1
         do k=ndz_v,nfz_v
            do i=ndx_v,nfx_v
               Krhou(i,j,k) = ( a04(1)*Grhou(i,j  ,k)+a04(2)*Grhou(i,j+1,k) &
                              + a04(3)*Grhou(i,j+2,k)+a04(4)*Grhou(i,j+3,k) &
                              + a04(5)*Grhou(i,j+4,k) )*idy_v(j) + Krhou(i,j,k)
               Krhov(i,j,k) = ( a04(1)*Grhov(i,j  ,k)+a04(2)*Grhov(i,j+1,k) &
                              + a04(3)*Grhov(i,j+2,k)+a04(4)*Grhov(i,j+3,k) &
                              + a04(5)*Grhov(i,j+4,k) )*idy_v(j) + Krhov(i,j,k)
               Krhow(i,j,k) = ( a04(1)*Grhow(i,j  ,k)+a04(2)*Grhow(i,j+1,k) &
                              + a04(3)*Grhow(i,j+2,k)+a04(4)*Grhow(i,j+3,k) &
                              + a04(5)*Grhow(i,j+4,k) )*idy_v(j) + Krhow(i,j,k)
               Krhoe(i,j,k) = ( a04(1)*Grhoe(i,j  ,k)+a04(2)*Grhoe(i,j+1,k) &
                              + a04(3)*Grhoe(i,j+2,k)+a04(4)*Grhoe(i,j+3,k) &
                              + a04(5)*Grhoe(i,j+4,k) )*idy_v(j) + Krhoe(i,j,k)
            enddo
         enddo
      endif
      if (is_bc_1pt(2,1)) then
         j=2
         do k=ndz_v,nfz_v
            do i=ndx_v,nfx_v
               Krhou(i,j,k) = ( a13(1)*Grhou(i,j-1,k)+a13(2)*Grhou(i,j,k) &
                              + a13(3)*Grhou(i,j+1,k)+a13(4)*Grhou(i,j+2,k) &
                              + a13(5)*Grhou(i,j+3,k) )*idy_v(j) + Krhou(i,j,k)
               Krhov(i,j,k) = ( a13(1)*Grhov(i,j-1,k)+a13(2)*Grhov(i,j,k) &
                              + a13(3)*Grhov(i,j+1,k)+a13(4)*Grhov(i,j+2,k) &
                              + a13(5)*Grhov(i,j+3,k) )*idy_v(j) + Krhov(i,j,k)
               Krhow(i,j,k) = ( a13(1)*Grhow(i,j-1,k)+a13(2)*Grhow(i,j,k) &
                              + a13(3)*Grhow(i,j+1,k)+a13(4)*Grhow(i,j+2,k) &
                              + a13(5)*Grhow(i,j+3,k) )*idy_v(j) + Krhow(i,j,k)
               Krhoe(i,j,k) = ( a13(1)*Grhoe(i,j-1,k)+a13(2)*Grhoe(i,j,k) &
                              + a13(3)*Grhoe(i,j+1,k)+a13(4)*Grhoe(i,j+2,k) &
                              + a13(5)*Grhoe(i,j+3,k) )*idy_v(j) + Krhoe(i,j,k)
            enddo
         enddo
      endif
    
      do k=ndz_v,nfz_v
         do j=ndy_vi,nfy_vi
            do i=ndx_v,nfx_v
               Krhou(i,j,k) = ( a5(1)*( Grhou(i,j+1,k)-Grhou(i,j-1,k)) &
                              + a5(2)*( Grhou(i,j+2,k)-Grhou(i,j-2,k)) )*idy_v(j) + Krhou(i,j,k)
               Krhov(i,j,k) = ( a5(1)*( Grhov(i,j+1,k)-Grhov(i,j-1,k)) &
                              + a5(2)*( Grhov(i,j+2,k)-Grhov(i,j-2,k)) )*idy_v(j) + Krhov(i,j,k)
               Krhow(i,j,k) = ( a5(1)*( Grhow(i,j+1,k)-Grhow(i,j-1,k)) &
                              + a5(2)*( Grhow(i,j+2,k)-Grhow(i,j-2,k)) )*idy_v(j) + Krhow(i,j,k)
               Krhoe(i,j,k) = ( a5(1)*( Grhoe(i,j+1,k)-Grhoe(i,j-1,k)) &
                              + a5(2)*( Grhoe(i,j+2,k)-Grhoe(i,j-2,k)) )*idy_v(j) + Krhoe(i,j,k)
            enddo
         enddo
      enddo
    
      if (is_bc_1pt(2,2)) then
         j=ny-1
         do k=ndz_v,nfz_v
            do i=ndx_v,nfx_v
               Krhou(i,j,k) = ( a31(1)*Grhou(i,j+1,k)+a31(2)*Grhou(i,j,k) &
                              + a31(3)*Grhou(i,j-1,k)+a31(4)*Grhou(i,j-2,k) &
                              + a31(5)*Grhou(i,j-3,k) )*idy_v(j) + Krhou(i,j,k)
               Krhov(i,j,k) = ( a31(1)*Grhov(i,j+1,k)+a31(2)*Grhov(i,j,k) &
                              + a31(3)*Grhov(i,j-1,k)+a31(4)*Grhov(i,j-2,k) &
                              + a31(5)*Grhov(i,j-3,k) )*idy_v(j) + Krhov(i,j,k)
               Krhow(i,j,k) = ( a31(1)*Grhow(i,j+1,k)+a31(2)*Grhow(i,j,k) &
                              + a31(3)*Grhow(i,j-1,k)+a31(4)*Grhow(i,j-2,k) &
                              + a31(5)*Grhow(i,j-3,k) )*idy_v(j) + Krhow(i,j,k)
               Krhoe(i,j,k) = ( a31(1)*Grhoe(i,j+1,k)+a31(2)*Grhoe(i,j,k) &
                              + a31(3)*Grhoe(i,j-1,k)+a31(4)*Grhoe(i,j-2,k) &
                              + a31(5)*Grhoe(i,j-3,k) )*idy_v(j) + Krhoe(i,j,k)
            enddo
         enddo
      endif
      if (is_bc_wall(2,2)) then
         j=ny
         do k=ndz_v,nfz_v
            do i=ndx_v,nfx_v
               Krhou(i,j,k) = ( a40(1)*Grhou(i,j,k)  +a40(2)*Grhou(i,j-1,k) &
                              + a40(3)*Grhou(i,j-2,k)+a40(4)*Grhou(i,j-3,k) &
                              + a40(5)*Grhou(i,j-4,k) )*idy_v(j) + Krhou(i,j,k)
               Krhov(i,j,k) = ( a40(1)*Grhov(i,j,k)  +a40(2)*Grhov(i,j-1,k) &
                              + a40(3)*Grhov(i,j-2,k)+a40(4)*Grhov(i,j-3,k) &
                              + a40(5)*Grhov(i,j-4,k) )*idy_v(j) + Krhov(i,j,k)
               Krhow(i,j,k) = ( a40(1)*Grhow(i,j,k)  +a40(2)*Grhow(i,j-1,k) &
                              + a40(3)*Grhow(i,j-2,k)+a40(4)*Grhow(i,j-3,k) &
                              + a40(5)*Grhow(i,j-4,k) )*idy_v(j) + Krhow(i,j,k)
               Krhoe(i,j,k) = ( a40(1)*Grhoe(i,j,k)  +a40(2)*Grhoe(i,j-1,k) &
                              + a40(3)*Grhoe(i,j-2,k)+a40(4)*Grhoe(i,j-3,k) &
                              + a40(5)*Grhoe(i,j-4,k) )*idy_v(j) + Krhoe(i,j,k)
            enddo
         enddo
      endif
    
      if (is_2D) return
    
      ! Compute derivatives of viscous fluxes along z
      ! =============================================
      if (is_bc_wall(3,1)) then
         k=1
         do j=ndy_v,wm_ind-1
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
         do j=ndy_v,wm_ind-1
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
         do j=ndy_v,wm_ind-1
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
         do j=ndy_v,wm_ind-1
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
         do j=ndy_v,wm_ind-1
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
    
    end subroutine flux_visc_5pts_SM_wall
    
    
    
    !==============================================================================
    subroutine filtre_Sij_wall
    !==============================================================================
      !> Filtering of Sij
      !> - Cartesian version -
    !==============================================================================
      use mod_coeff_deriv
      use mod_interface
      use mod_flow
      use mod_wall_model
      implicit none
      integer :: i,j,k
      real :: d11t(0:5)
      real, dimension(nx1:nx2,ny1:ny2,nz1:nz2) :: DS11,DS22,DS33,DS12,DS13,DS23
    
      call communication_(S11,S22,S33,S12,S13)
      call communication_(S11,S22,S33,S12,S23) ! quasi la meme chose 2 fois, a changer
    
      ! Initialisation
      ! ---------------
      S11f = 0.
      S22f = 0.
      S33f = 0.
      S12f = 0.
      S13f = 0.
      S23f = 0.
    
      DS11 = 0.
      DS22 = 0.
      DS33 = 0.
      DS12 = 0.
      DS13 = 0.
      DS23 = 0.
    
      ! Coupure en pi/3
      ! ---------------
      d11t(0)= 2./3.
      d11t(1)=-0.26775782
      d11t(2)=-0.12016956
      d11t(3)= 0.
      d11t(4)= 0.03683622
      d11t(5)= 0.01775782
    
      ! Filtrage en x sur 11 points
      ! ---------------------------
      do k=ndz_v,nfz_v
         do j=ndy_v,wm_ind-1
            do i=ndx_v,nfx_v
    
               DS11(i,j,k)= ( d11t(0) * S11(i,j,k)        &
                            + d11t(1) * ( S11(i+1,j,k)+S11(i-1,j,k) ) &
                            + d11t(2) * ( S11(i+2,j,k)+S11(i-2,j,k) ) &
                            + d11t(3) * ( S11(i+3,j,k)+S11(i-3,j,k) ) &
                            + d11t(4) * ( S11(i+4,j,k)+S11(i-4,j,k) ) &
                            + d11t(5) * ( S11(i+5,j,k)+S11(i-5,j,k) ) )
    
               DS22(i,j,k)= ( d11t(0) * S22(i,j,k)        &
                            + d11t(1) * ( S22(i+1,j,k)+S22(i-1,j,k) ) &
                            + d11t(2) * ( S22(i+2,j,k)+S22(i-2,j,k) ) &
                            + d11t(3) * ( S22(i+3,j,k)+S22(i-3,j,k) ) &
                            + d11t(4) * ( S22(i+4,j,k)+S22(i-4,j,k) )   &
                            + d11t(5) * ( S22(i+5,j,k)+S22(i-5,j,k) ) )
    
               DS33(i,j,k)= ( d11t(0) * S33(i,j,k)        &
                            + d11t(1) * ( S33(i+1,j,k)+S33(i-1,j,k) ) &
                            + d11t(2) * ( S33(i+2,j,k)+S33(i-2,j,k) ) &
                            + d11t(3) * ( S33(i+3,j,k)+S33(i-3,j,k) ) &
                            + d11t(4) * ( S33(i+4,j,k)+S33(i-4,j,k) ) &
                            + d11t(5) * ( S33(i+5,j,k)+S33(i-5,j,k) ) )
    
               DS12(i,j,k)= ( d11t(0) * S12(i,j,k)        &
                            + d11t(1) * ( S12(i+1,j,k)+S12(i-1,j,k) ) &
                            + d11t(2) * ( S12(i+2,j,k)+S12(i-2,j,k) ) &
                            + d11t(3) * ( S12(i+3,j,k)+S12(i-3,j,k) ) &
                            + d11t(4) * ( S12(i+4,j,k)+S12(i-4,j,k) ) &
                            + d11t(5) * ( S12(i+5,j,k)+S12(i-5,j,k) ) )
    
               DS13(i,j,k)= ( d11t(0) * S13(i,j,k)              &
                            + d11t(1) * ( S13(i+1,j,k)+S13(i-1,j,k) ) &
                            + d11t(2) * ( S13(i+2,j,k)+S13(i-2,j,k) ) &
                            + d11t(3) * ( S13(i+3,j,k)+S13(i-3,j,k) ) &
                            + d11t(4) * ( S13(i+4,j,k)+S13(i-4,j,k) ) &
                            + d11t(5) * ( S13(i+5,j,k)+S13(i-5,j,k) ) )
    
               DS23(i,j,k)= ( d11t(0) * S23(i,j,k)              &
                            + d11t(1) * ( S23(i+1,j,k)+S23(i-1,j,k) ) &
                            + d11t(2) * ( S23(i+2,j,k)+S23(i-2,j,k) ) &
                            + d11t(3) * ( S23(i+3,j,k)+S23(i-3,j,k) ) &
                            + d11t(4) * ( S23(i+4,j,k)+S23(i-4,j,k) ) &
                            + d11t(5) * ( S23(i+5,j,k)+S23(i-5,j,k) ) )
    
            enddo
         enddo
      enddo
    
      ! Filtrage des Sij
      ! ----------------
    
      do k=ndz_v,nfz_v
         do j=ndy_v,wm_ind-1
            do i=ndx_v,nfx_v
               S11f(i,j,k) = DS11(i,j,k)
               S22f(i,j,k) = DS22(i,j,k)
               S33f(i,j,k) = DS33(i,j,k)
               S12f(i,j,k) = DS12(i,j,k)
               S13f(i,j,k) = DS13(i,j,k)
               S23f(i,j,k) = DS23(i,j,k)
            enddo
         enddo
      enddo
    
    
      ! Filtrage en y sur 11 points
      ! ---------------------------
    
      do k=ndz_v,nfz_v
         do j=ndy_v,wm_ind-1
            do i=ndx_v,nfx_v
    
               DS11(i,j,k) =  ( d11t(0) * S11(i,j,k)  &
                              + d11t(1) * ( S11(i,j+1,k)+S11(i,j-1,k) ) &
                              + d11t(2) * ( S11(i,j+2,k)+S11(i,j-2,k) )     &
                              + d11t(3) * ( S11(i,j+3,k)+S11(i,j-3,k) )     &
                              + d11t(4) * ( S11(i,j+4,k)+S11(i,j-4,k) )     &
                              + d11t(5) * ( S11(i,j+5,k)+S11(i,j-5,k) ) )
    
               DS22(i,j,k)  = ( d11t(0) * S22(i,j,k)  &
                              + d11t(1) * ( S22(i,j+1,k)+S22(i,j-1,k) ) &
                              + d11t(2) * ( S22(i,j+2,k)+S22(i,j-2,k) )     &
                              + d11t(3) * ( S22(i,j+3,k)+S22(i,j-3,k) )     &
                              + d11t(4) * ( S22(i,j+4,k)+S22(i,j-4,k) )     &
                              + d11t(5) * ( S22(i,j+5,k)+S22(i,j-5,k) ) )
    
               DS33(i,j,k)  = ( d11t(0) * S33(i,j,k)  &
                              + d11t(1) * ( S33(i,j+1,k)+S33(i,j-1,k) ) &
                              + d11t(2) * ( S33(i,j+2,k)+S33(i,j-2,k) )     &
                              + d11t(3) * ( S33(i,j+3,k)+S33(i,j-3,k) )     &
                              + d11t(4) * ( S33(i,j+4,k)+S33(i,j-4,k) )     &
                              + d11t(5) * ( S33(i,j+5,k)+S33(i,j-5,k) ) )
    
               DS12(i,j,k)  = ( d11t(0) * S12(i,j,k)  &
                              + d11t(1) * ( S12(i,j+1,k)+S12(i,j-1,k) ) &
                              + d11t(2) * ( S12(i,j+2,k)+S12(i,j-2,k) )     &
                              + d11t(3) * ( S12(i,j+3,k)+S12(i,j-3,k) )     &
                              + d11t(4) * ( S12(i,j+4,k)+S12(i,j-4,k) )     &
                              + d11t(5) * ( S12(i,j+5,k)+S12(i,j-5,k) ) )
    
               DS13(i,j,k)  = ( d11t(0) * S13(i,j,k)  &
                              + d11t(1) * ( S13(i,j+1,k)+S13(i,j-1,k) ) &
                              + d11t(2) * ( S13(i,j+2,k)+S13(i,j-2,k) )     &
                              + d11t(3) * ( S13(i,j+3,k)+S13(i,j-3,k) )     &
                              + d11t(4) * ( S13(i,j+4,k)+S13(i,j-4,k) )     &
                              + d11t(5) * ( S13(i,j+5,k)+S13(i,j-5,k) ) )
    
               DS23(i,j,k)  = ( d11t(0) * S23(i,j,k)  &
                              + d11t(1) * ( S23(i,j+1,k)+S23(i,j-1,k) ) &
                              + d11t(2) * ( S23(i,j+2,k)+S23(i,j-2,k) )     &
                              + d11t(3) * ( S23(i,j+3,k)+S23(i,j-3,k) )     &
                              + d11t(4) * ( S23(i,j+4,k)+S23(i,j-4,k) )     &
                              + d11t(5) * ( S23(i,j+5,k)+S23(i,j-5,k) ) )
    
            enddo
         enddo
      enddo
    
      ! Filtrage des Sij
      ! ----------------
    
      do k=ndz_v,nfz_v
         do j=ndy_v,wm_ind-1
            do i=ndx_v,nfx_v
               S11f(i,j,k) = S11f(i,j,k) + DS11(i,j,k)
               S22f(i,j,k) = S22f(i,j,k) + DS22(i,j,k)
               S33f(i,j,k) = S33f(i,j,k) + DS33(i,j,k)
               S12f(i,j,k) = S12f(i,j,k) + DS12(i,j,k)
               S13f(i,j,k) = S13f(i,j,k) + DS13(i,j,k)
               S23f(i,j,k) = S23f(i,j,k) + DS23(i,j,k)
            enddo
         enddo
      enddo
    
    
      ! Filtrage en z sur 11 points
      ! ---------------------------
    
      do k=ndz_v,nfz_v
         do j=ndy_v,wm_ind-1
            do i=ndx_v,nfx_v
    
               DS11(i,j,k) =  ( d11t(0) * S11(i,j,k)  &
                              + d11t(1) * ( S11(i,j,k+1)+S11(i,j,k-1) ) &
                              + d11t(2) * ( S11(i,j,k+2)+S11(i,j,k-2) )     &
                              + d11t(3) * ( S11(i,j,k+3)+S11(i,j,k-3) )     &
                              + d11t(4) * ( S11(i,j,k+4)+S11(i,j,k-4) )     &
                              + d11t(5) * ( S11(i,j,k+5)+S11(i,j,k-5) ) )
    
               DS22(i,j,k)  = ( d11t(0) * S22(i,j,k)  &
                              + d11t(1) * ( S22(i,j,k+1)+S22(i,j,k-1) ) &
                              + d11t(2) * ( S22(i,j,k+2)+S22(i,j,k-2) )     &
                              + d11t(3) * ( S22(i,j,k+3)+S22(i,j,k-3) )     &
                              + d11t(4) * ( S22(i,j,k+4)+S22(i,j,k-4) )     &
                              + d11t(5) * ( S22(i,j,k+5)+S22(i,j,k-5) ) )
    
               DS33(i,j,k)  = ( d11t(0) * S33(i,j,k)  &
                              + d11t(1) * ( S33(i,j,k+1)+S33(i,j,k-1) ) &
                              + d11t(2) * ( S33(i,j,k+2)+S33(i,j,k-2) )     &
                              + d11t(3) * ( S33(i,j,k+3)+S33(i,j,k-3) )     &
                              + d11t(4) * ( S33(i,j,k+4)+S33(i,j,k-4) )     &
                              + d11t(5) * ( S33(i,j,k+5)+S33(i,j,k-5) ) )
    
               DS12(i,j,k)  = ( d11t(0) * S12(i,j,k)  &
                              + d11t(1) * ( S12(i,j,k+1)+S12(i,j,k-1) ) &
                              + d11t(2) * ( S12(i,j,k+2)+S12(i,j,k-2) )     &
                              + d11t(3) * ( S12(i,j,k+3)+S12(i,j,k-3) )     &
                              + d11t(4) * ( S12(i,j,k+4)+S12(i,j,k-4) )     &
                              + d11t(5) * ( S12(i,j,k+5)+S12(i,j,k-5) ) )
    
               DS13(i,j,k)  = ( d11t(0) * S13(i,j,k)  &
                              + d11t(1) * ( S13(i,j,k+1)+S13(i,j,k-1) ) &
                              + d11t(2) * ( S13(i,j,k+2)+S13(i,j,k-2) )     &
                              + d11t(3) * ( S13(i,j,k+3)+S13(i,j,k-3) )     &
                              + d11t(4) * ( S13(i,j,k+4)+S13(i,j,k-4) )     &
                              + d11t(5) * ( S13(i,j,k+5)+S13(i,j,k-5) ) )
    
               DS23(i,j,k)  = ( d11t(0) * S23(i,j,k)  &
                              + d11t(1) * ( S23(i,j,k+1)+S23(i,j,k-1) ) &
                              + d11t(2) * ( S23(i,j,k+2)+S23(i,j,k-2) )     &
                              + d11t(3) * ( S23(i,j,k+3)+S23(i,j,k-3) )     &
                              + d11t(4) * ( S23(i,j,k+4)+S23(i,j,k-4) )     &
                              + d11t(5) * ( S23(i,j,k+5)+S23(i,j,k-5) ) )
    
            enddo
         enddo
      enddo
    
      ! Filtrage des Sij
      ! ----------------
    
      do k=ndz_v,nfz_v
         do j=ndy_v,wm_ind-1
            do i=ndx_v,nfx_v
               S11f(i,j,k) = S11f(i,j,k) + DS11(i,j,k)
               S22f(i,j,k) = S22f(i,j,k) + DS22(i,j,k)
               S33f(i,j,k) = S33f(i,j,k) + DS33(i,j,k)
               S12f(i,j,k) = S12f(i,j,k) + DS12(i,j,k)
               S13f(i,j,k) = S13f(i,j,k) + DS13(i,j,k)
               S23f(i,j,k) = S23f(i,j,k) + DS23(i,j,k)
            enddo
         enddo
      enddo
    
end subroutine filtre_Sij_wall

