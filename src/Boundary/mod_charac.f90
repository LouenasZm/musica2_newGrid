!============================================================================
module mod_charac
!============================================================================
  !> author: XG
  !> date: February 2020
  !> Characteristic boundary conditions - Cartesian version -
  !> Att. restricted to perfect gas
!=============================================================================
  use mod_flow
  use mod_init_flow
  use mod_constant
  use mod_coeff_deriv
  use mod_mpi
  use mod_time
  use mod_eos
  implicit none
  !---------------------------------------------------------------------------
  real(wp) :: cp,av,cj ! thermodyn.
  ! real(wp) :: rt,ut,vt,wt,pt ! time derivatives
  real(wp) :: L1,L2,L3,L4,L5 ! charac. amplitudes
  !---------------------------------------------------------------------------

contains

  !===============================================================================
  subroutine bc_charac_imin
  !===============================================================================
    !> Characteristic BC: boundary condition at imin (left) - Cartesian version -
  !===============================================================================
    use mod_eigenmode
    use mod_RFM
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: j,k
    ! real(wp), dimension(nz) :: dry,dpy
    real(wp), dimension(ny,nz) :: drx,dpx
    real(wp), dimension(ny,nz) :: rt,ut,vt,wt,pt
    ! ----------------------------------------------------------------------------
    ! eigenmode disturbances
    !real(wp) :: um,vm,wm,rm,tm ! disturbances
    real(wp) :: dumt,dvmt,dwmt,drmt,dpmt ! time derivatives of disturbances
    ! RFM disturbances
    real(wp), dimension(ny,nz) :: ut_in,vt_in,wt_in
    ! ---------------------------------------------------------------------------

    ! Index of left boundary: i=1

    ! Compute derivatives along x
    ! ===========================
    if (is_curv) then
       do k=1,nz
          do j=1,ny
             dux(1,j,k)=( a04(1)*uu(1,j,k)+a04(2)*uu(2,j,k) &
                         +a04(3)*uu(3,j,k)+a04(4)*uu(4,j,k) &
                         +a04(5)*uu(5,j,k) )*y_eta_v(1,j)*ijacob_v(1,j)
             dvx(1,j,k)=( a04(1)*vv(1,j,k)+a04(2)*vv(2,j,k) &
                         +a04(3)*vv(3,j,k)+a04(4)*vv(4,j,k) &
                         +a04(5)*vv(5,j,k) )*y_eta_v(1,j)*ijacob_v(1,j)
             dwx(1,j,k)=( a04(1)*ww(1,j,k)+a04(2)*ww(2,j,k) &
                         +a04(3)*ww(3,j,k)+a04(4)*ww(4,j,k) &
                         +a04(5)*ww(5,j,k) )*y_eta_v(1,j)*ijacob_v(1,j)
             dpx(  j,k)=( a04(1)*prs(1,j,k)+a04(2)*prs(2,j,k) &
                         +a04(3)*prs(3,j,k)+a04(4)*prs(4,j,k) &
                         +a04(5)*prs(5,j,k) )*y_eta_v(1,j)*ijacob_v(1,j)
             drx(  j,k)=( a04(1)*rho_n(1,j,k)+a04(2)*rho_n(2,j,k) &
                         +a04(3)*rho_n(3,j,k)+a04(4)*rho_n(4,j,k) &
                         +a04(5)*rho_n(5,j,k) )*y_eta_v(1,j)*ijacob_v(1,j)
          enddo
       enddo
    else
       do k=1,nz
          do j=1,ny
             dux(1,j,k)=( a04(1)*uu(1,j,k)+a04(2)*uu(2,j,k) &
                         +a04(3)*uu(3,j,k)+a04(4)*uu(4,j,k) &
                         +a04(5)*uu(5,j,k) )*idx_v(1)
             dvx(1,j,k)=( a04(1)*vv(1,j,k)+a04(2)*vv(2,j,k) &
                         +a04(3)*vv(3,j,k)+a04(4)*vv(4,j,k) &
                         +a04(5)*vv(5,j,k) )*idx_v(1)
             dwx(1,j,k)= (a04(1)*ww(1,j,k)+a04(2)*ww(2,j,k) &
                         +a04(3)*ww(3,j,k)+a04(4)*ww(4,j,k) &
                         +a04(5)*ww(5,j,k) )*idx_v(1)
             dpx(  j,k)= (a04(1)*prs(1,j,k)+a04(2)*prs(2,j,k) &
                         +a04(3)*prs(3,j,k)+a04(4)*prs(4,j,k) &
                         +a04(5)*prs(5,j,k) )*idx_v(1)
             drx(  j,k)= (a04(1)*rho_n(1,j,k)+a04(2)*rho_n(2,j,k) &
                         +a04(3)*rho_n(3,j,k)+a04(4)*rho_n(4,j,k) &
                         +a04(5)*rho_n(5,j,k) )*idx_v(1)
          enddo
       enddo
    endif

    if (is_RFM) call disturb_inlet_RFM_charac(ut_in,vt_in,wt_in)
    
    do k=1,nz
       do j=1,ny
          ! Compute characteristic amplitudes
          ! =================================
          L5 = 0.0_wp
          L1 = (uu(1,j,k)-c_(1,j,k))*(dpx(j,k)-rho_n(1,j,k)*c_(1,j,k)*dux(1,j,k))
          if (uu(1,j,k).ge.c_(1,j,k)) L1 = 0.0_wp
          !L1 = 0.0_wp
          L2 = 0.0_wp
          L3 = 0.0_wp
          L4 = 0.0_wp

          !if (j==100) print *,L1,uu(1,j,k),c_(1,j,k)

          ! Time derivatives for primitive variables
          ! ========================================
          pt(j,k) = 0.5_wp*(L5+L1)
          rt(j,k) = -(pt(j,k)+L2)/(c_(1,j,k)**2)
          ut(j,k) = (L5-L1)/(2.0_wp*rho_n(1,j,k)*c_(1,j,k))
          vt(j,k) = L3
          wt(j,k) = L4
          
          ! compute alpha*(time derivatives of disturb)
          ! -------------------------------------------
          if (is_eigenmode) then ! add eigenmodes
             call eig_disturb1_dt(1,j,k,dumt,dvmt,dwmt,drmt,dpmt)
             pt(j,k) = pt(j,k) - dpmt
             rt(j,k) = rt(j,k) - drmt
             ut(j,k) = ut(j,k) - dumt
             vt(j,k) = vt(j,k) - dvmt
             wt(j,k) = wt(j,k) - dwmt
          endif
          
          if (is_RFM) then ! add eigenmodes
             ut(j,k) = ut(j,k) - ut_in(j,k)
             vt(j,k) = vt(j,k) - vt_in(j,k)
             wt(j,k) = wt(j,k) - wt_in(j,k)
          endif
          
       enddo
    enddo

    ! Update fluxes at each RK step
    ! =============================
    do k=1,nz
       do j=1,ny
          cp = cpcalc_tro(Tmp(1,j,k),rho_n(1,j,k))
          av = avcalc_tro(Tmp(1,j,k),rho_n(1,j,k))
          Krho(1,j,k)  = rt(j,k)
          Krhou(1,j,k) = uu(1,j,k)*rt(j,k)+rho_n(1,j,k)*ut(j,k)
          Krhov(1,j,k) = vv(1,j,k)*rt(j,k)+rho_n(1,j,k)*vt(j,k)
          Krhow(1,j,k) = ww(1,j,k)*rt(j,k)+rho_n(1,j,k)*wt(j,k)
          Krhoe(1,j,k) = cp/av*(pt(j,k)/c_(1,j,k)**2 - rt(j,k)) &
               + (rhoe_n(1,j,k)+prs(1,j,k))/rho_n(1,j,k)*rt(j,k) &
               + rho_n(1,j,k)*(uu(1,j,k)*ut(j,k)+vv(1,j,k)*vt(j,k)+ww(1,j,k)*wt(j,k))
       enddo
    enddo

!!$    ! Flux derivatives along y
!!$    ! ========================
!!$
!!$    if ((coord(2)==0).and.(BC_face(2,1)%sort==-3)) then
!!$       ! left-bottom edge
!!$       ! ================
!!$       j=1
!!$
!!$       ! Compute derivatives along y
!!$       ! ===========================
!!$       if (is_curv) then
!!$          do k=1,nz
!!$             duy(i,j,k)= ( a04(1)*uu(i,j  ,k)+a04(2)*uu(i,j+1,k) &
!!$                  +a04(3)*uu(i,j+2,k)+a04(4)*uu(i,j+3,k) &
!!$                  +a04(5)*uu(i,j+4,k) )*x_ksi_v(i,j)*ijacob_v(i,j)
!!$             dvy(i,j,k)= ( a04(1)*vv(i,j  ,k)+a04(2)*vv(i,j+1,k) &
!!$                  +a04(3)*vv(i,j+2,k)+a04(4)*vv(i,j+3,k) &
!!$                  +a04(5)*vv(i,j+4,k) )*x_ksi_v(i,j)*ijacob_v(i,j)
!!$             dwy(i,j,k)= ( a04(1)*ww(i,j  ,k)+a04(2)*ww(i,j+1,k) &
!!$                  +a04(3)*ww(i,j+2,k)+a04(4)*ww(i,j+3,k) &
!!$                  +a04(5)*ww(i,j+4,k) )*x_ksi_v(i,j)*ijacob_v(i,j)
!!$             dpy(k)    = ( a04(1)*prs(i,j  ,k)+a04(2)*prs(i,j+1,k) &
!!$                  +a04(3)*prs(i,j+2,k)+a04(4)*prs(i,j+3,k) &
!!$                  +a04(5)*prs(i,j+4,k) )*x_ksi_v(i,j)*ijacob_v(i,j)
!!$             dry(k)    = ( a04(1)*rho_n(i,j  ,k)+a04(2)*rho_n(i,j+1,k) &
!!$                  +a04(3)*rho_n(i,j+2,k)+a04(4)*rho_n(i,j+3,k) &
!!$                  +a04(5)*rho_n(i,j+4,k) )*x_ksi_v(i,j)*ijacob_v(i,j)
!!$          enddo
!!$       else
!!$          do k=1,nz
!!$             duy(i,j,k)= ( a04(1)*uu(i,j  ,k)+a04(2)*uu(i,j+1,k) &
!!$                  +a04(3)*uu(i,j+2,k)+a04(4)*uu(i,j+3,k) &
!!$                  +a04(5)*uu(i,j+4,k) )*idy_v(j)
!!$             dvy(i,j,k)= ( a04(1)*vv(i,j  ,k)+a04(2)*vv(i,j+1,k) &
!!$                  +a04(3)*vv(i,j+2,k)+a04(4)*vv(i,j+3,k) &
!!$                  +a04(5)*vv(i,j+4,k) )*idy_v(j)
!!$             dwy(i,j,k)= ( a04(1)*ww(i,j  ,k)+a04(2)*ww(i,j+1,k) &
!!$                  +a04(3)*ww(i,j+2,k)+a04(4)*ww(i,j+3,k) &
!!$                  +a04(5)*ww(i,j+4,k) )*idy_v(j)
!!$             dpy(k)    = ( a04(1)*prs(i,j  ,k)+a04(2)*prs(i,j+1,k) &
!!$                  +a04(3)*prs(i,j+2,k)+a04(4)*prs(i,j+3,k) &
!!$                  +a04(5)*prs(i,j+4,k) )*idy_v(j)
!!$             dry(k)    = ( a04(1)*rho_n(i,j  ,k)+a04(2)*rho_n(i,j+1,k) &
!!$                  +a04(3)*rho_n(i,j+2,k)+a04(4)*rho_n(i,j+3,k) &
!!$                  +a04(5)*rho_n(i,j+4,k) )*idy_v(j)
!!$          enddo
!!$       endif
!!$
!!$       ! Update sound velocity (necessary if first runge-kutta step ?)
!!$       ! =====================
!!$       do k=1,nz
!!$          c_(i,j,k) = (c2calc_tro(Tmp(i,j,k),rho_n(i,j,k)))**0.5
!!$       enddo
!!$
!!$       do k=1,nz
!!$          ! Compute characteristic amplitudes
!!$          ! =================================
!!$          L1 = (vv(i,j,k)-c_(i,j,k))*(dpy(k)-rho_n(i,j,k)*c_(i,j,k)*dvy(i,j,k))
!!$          ! L1 = (v0_sud(i,k)-c_(i,j,k))*(dpy(k)-rho_n(i,j,k)*c_(i,j,k)*dvy(i,j,k))
!!$          L5 = 0.0_wp
!!$
!!$          if (vv(i,j,k).ge.0.0_wp) then
!!$          ! if (v0_sud(i,k).ge.0.0_wp) then
!!$             L2 = 0.0_wp
!!$             L3 = 0.0_wp
!!$             L4 = 0.0_wp
!!$             if (vv(i,j,k).ge.c_(i,j,k)) L1 = 0.0_wp
!!$             ! if (v0_sud(i,k).ge.c_(i,j,k)) L1 = 0.0_wp
!!$          else
!!$             L2 = vv(i,j,k)*duy(i,j,k)
!!$             L3 = vv(i,j,k)*(dpy(k)-(c_(i,j,k)**2)*dry(k))
!!$             L4 = vv(i,j,k)*dwy(i,j,k)
!!$             if (vv(i,j,k).le.-c_(i,j,k)) L5=(vv(i,j,k)+c_(i,j,k))*(dpy(k)+rho_n(i,j,k)*c_(i,j,k)*dvy(i,j,k))
!!$             ! L2 = v0_sud(i,k)*duy(i,j,k)
!!$             ! L3 = v0_sud(i,k)*(dpy(k)-(c_(i,j,k)**2)*dry(k))
!!$             ! L4 = v0_sud(i,k)*dwy(i,j,k)
!!$             ! if (v0_sud(i,k).le.-c_(i,j,k)) L5=(v0_sud(i,k)+c_(i,j,k))*(dpy(k)+rho_n(i,j,k)*c_(i,j,k)*dvy(i,j,k))
!!$          endif
!!$
!!$          ! Time derivatives for primitive variables
!!$          ! ========================================
!!$          pt(j,k) = 0.5_wp*(L5+L1)
!!$          rt(j,k) = -(pt(j,k)+L3)/(c_(i,j,k)**2)
!!$          ut(j,k) = L2
!!$          vt(j,k) = (L5-L1)/(2.0_wp*rho_n(i,j,k)*c_(i,j,k))
!!$          wt(j,k) = L4
!!$       enddo
!!$
!!$       do k=1,nz
!!$          ! Update fluxes at each RK step
!!$          ! =============================
!!$          cp = cpcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
!!$          av = avcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
!!$          Krho(i,j,k)  = Krho(i,j,k)  + rt(j,k)
!!$          Krhou(i,j,k) = Krhou(i,j,k) + rho_n(i,j,k)*ut(j,k) + uu(i,j,k)*rt(j,k)
!!$          Krhov(i,j,k) = Krhov(i,j,k) + rho_n(i,j,k)*vt(j,k) + vv(i,j,k)*rt(j,k)
!!$          Krhow(i,j,k) = Krhow(i,j,k) + rho_n(i,j,k)*wt(j,k) + ww(i,j,k)*rt(j,k)
!!$          Krhoe(i,j,k) = Krhoe(i,j,k) +  cp/av*(pt(j,k)/c_(i,j,k)**2 - rt(j,k)) &
!!$               + (rhoe_n(i,j,k)+prs(i,j,k))/rho_n(i,j,k)*rt(j,k) &
!!$               + rho_n(i,j,k)*(uu(i,j,k)*ut(j,k)+vv(i,j,k)*vt(j,k)+ww(i,j,k)*wt(j,k))
!!$       enddo
!!$
!!$    endif
!!$
!!$    if ((coord(2)==ndomy-1).and.(BC_face(2,2)%sort==-3))then
!!$       ! left-top edge
!!$       ! =============
!!$       j=ny
!!$
!!$       ! Compute derivatives along y
!!$       ! ===========================
!!$       if (is_curv) then
!!$          do k=1,nz
!!$             duy(i,j,k)= ( a40(5)*uu(i,j-4,k)+a40(4)*uu(i,j-3,k) &
!!$                  +a40(3)*uu(i,j-2,k)+a40(2)*uu(i,j-1,k) &
!!$                  +a40(1)*uu(i,j  ,k) )*x_ksi_v(i,j)*ijacob_v(i,j)
!!$             dvy(i,j,k)= ( a40(5)*vv(i,j-4,k)+a40(4)*vv(i,j-3,k) &
!!$                  +a40(3)*vv(i,j-2,k)+a40(2)*vv(i,j-1,k) &
!!$                  +a40(1)*vv(i,j  ,k) )*x_ksi_v(i,j)*ijacob_v(i,j)
!!$             dwy(i,j,k)= ( a40(5)*ww(i,j-4,k)+a40(4)*ww(i,j-3,k) &
!!$                  +a40(3)*ww(i,j-2,k)+a40(2)*ww(i,j-1,k) &
!!$                  +a40(1)*ww(i,j  ,k) )*x_ksi_v(i,j)*ijacob_v(i,j)
!!$             dpy(k)    = ( a40(5)*prs(i,j-4,k)+a40(4)*prs(i,j-3,k) &
!!$                  +a40(3)*prs(i,j-2,k)+a40(2)*prs(i,j-1,k) &
!!$                  +a40(1)*prs(i,j  ,k) )*x_ksi_v(i,j)*ijacob_v(i,j)
!!$             dry(k)    = ( a40(5)*rho_n(i,j-4,k)+a40(4)*rho_n(i,j-3,k) &
!!$                  +a40(3)*rho_n(i,j-2,k)+a40(2)*rho_n(i,j-1,k) &
!!$                  +a40(1)*rho_n(i,j  ,k) )*x_ksi_v(i,j)*ijacob_v(i,j)
!!$          enddo
!!$       else
!!$          do k=1,nz
!!$             duy(i,j,k)= ( a40(5)*uu(i,j-4,k)+a40(4)*uu(i,j-3,k) &
!!$                  +a40(3)*uu(i,j-2,k)+a40(2)*uu(i,j-1,k) &
!!$                  +a40(1)*uu(i,j  ,k) )*idy_v(j)
!!$             dvy(i,j,k)= ( a40(5)*vv(i,j-4,k)+a40(4)*vv(i,j-3,k) &
!!$                  +a40(3)*vv(i,j-2,k)+a40(2)*vv(i,j-1,k) &
!!$                  +a40(1)*vv(i,j  ,k) )*idy_v(j)
!!$             dwy(i,j,k)= ( a40(5)*ww(i,j-4,k)+a40(4)*ww(i,j-3,k) &
!!$                  +a40(3)*ww(i,j-2,k)+a40(2)*ww(i,j-1,k) &
!!$                  +a40(1)*ww(i,j  ,k) )*idy_v(j)
!!$             dpy(k)    = ( a40(5)*prs(i,j-4,k)+a40(4)*prs(i,j-3,k) &
!!$                  +a40(3)*prs(i,j-2,k)+a40(2)*prs(i,j-1,k) &
!!$                  +a40(1)*prs(i,j  ,k) )*idy_v(j)
!!$             dry(k)    = ( a40(5)*rho_n(i,j-4,k)+a40(4)*rho_n(i,j-3,k) &
!!$                  +a40(3)*rho_n(i,j-2,k)+a40(2)*rho_n(i,j-1,k) &
!!$                  +a40(1)*rho_n(i,j  ,k) )*idy_v(j)
!!$          enddo
!!$       endif
!!$
!!$       ! Update sound velocity (necessary if first runge-kutta step ?)
!!$       ! =====================
!!$       do k=1,nz
!!$          c_(i,j,k) = (c2calc_tro(Tmp(i,j,k),rho_n(i,j,k)))**0.5
!!$       enddo
!!$
!!$       do k=1,nz
!!$          ! Compute characteristic amplitudes
!!$          ! =================================
!!$          L1 = 0.0_wp
!!$          L5 = (vv(i,j,k)+c_(i,j,k))*(dpy(k)+rho_n(i,j,k)*c_(i,j,k)*dvy(i,j,k))
!!$          if (vv(i,j,k).le.0.0_wp) then
!!$          ! if (v0_nord(i,k).le.0.0_wp) then
!!$             L2 = 0.0_wp
!!$             L3 = 0.0_wp
!!$             L4 = 0.0_wp
!!$             if (vv(i,j,k).le.-c_(i,j,k)) L5 = 0.0_wp
!!$          else
!!$             L2 = vv(i,j,k)*duy(i,j,k)
!!$             L3 = vv(i,j,k)*(dpy(k)-(c_(i,j,k)**2)*dry(k))
!!$             L4 = vv(i,j,k)*dwy(i,j,k)
!!$             if (vv(i,j,k).ge.c_(i,j,k)) L1=(vv(i,j,k)-c_(i,j,k))*(dpy(k)-rho_n(i,j,k)*c_(i,j,k)*dvy(i,j,k))
!!$          endif
!!$
!!$          ! Time derivatives for primitive variables
!!$          ! ========================================
!!$          pt(j,k) = 0.5_wp*(L5+L1)
!!$          rt(j,k) = -(pt(j,k)+L3)/(c_(i,j,k)**2)
!!$          ut(j,k) = L2
!!$          vt(j,k) = (L5-L1)/(2.0_wp*rho_n(i,j,k)*c_(i,j,k))
!!$          wt(j,k) = L4
!!$
!!$          ! Update fluxes at each RK step
!!$          ! =============================
!!$          cp = cpcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
!!$          av = avcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
!!$          Krho(i,j,k)  = Krho(i,j,k)  + rt(j,k)
!!$          Krhou(i,j,k) = Krhou(i,j,k) + rho_n(i,j,k)*ut(j,k) + uu(i,j,k)*rt(j,k)
!!$          Krhov(i,j,k) = Krhov(i,j,k) + rho_n(i,j,k)*vt(j,k) + vv(i,j,k)*rt(j,k)
!!$          Krhow(i,j,k) = Krhow(i,j,k) + rho_n(i,j,k)*wt(j,k) + ww(i,j,k)*rt(j,k)
!!$          Krhoe(i,j,k) = Krhoe(i,j,k) +  cp/av*(pt(j,k)/c_(i,j,k)**2 - rt(j,k)) &
!!$               + (rhoe_n(i,j,k)+prs(i,j,k))/rho_n(i,j,k)*rt(j,k) &
!!$               + rho_n(i,j,k)*(uu(i,j,k)*ut(j,k)+vv(i,j,k)*vt(j,k)+ww(i,j,k)*wt(j,k))
!!$       enddo
!!$    endif

  end subroutine bc_charac_imin

  !===============================================================================
  subroutine bc_charac_imax
  !===============================================================================
    !> Characteristic BC: boundary condition at imax (right) - Cartesian version -
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: j,k
    real(wp) :: cj,cj2,umc,upc
    ! real(wp), dimension(nz) :: dry,dpy
    real(wp), dimension(ny,nz) :: drx,dpx
    real(wp), dimension(ny,nz) :: rt,ut,vt,wt,pt
    ! ---------------------------------------------------------------------------

    ! Index of right boundary: i=nx

    ! Compute derivatives along x
    ! ===========================
    if (is_curv) then
       do k=1,nz
          do j=1,ny
             dux(nx,j,k)=( a40(5)*uu(nx-4,j,k)+a40(4)*uu(nx-3,j,k) &
                          +a40(3)*uu(nx-2,j,k)+a40(2)*uu(nx-1,j,k) &
                          +a40(1)*uu(nx  ,j,k) )*y_eta_v(nx,j)*ijacob_v(nx,j)
             dvx(nx,j,k)=( a40(5)*vv(nx-4,j,k)+a40(4)*vv(nx-3,j,k) &
                          +a40(3)*vv(nx-2,j,k)+a40(2)*vv(nx-1,j,k) &
                          +a40(1)*vv(nx  ,j,k) )*y_eta_v(nx,j)*ijacob_v(nx,j)
             dwx(nx,j,k)=( a40(5)*ww(nx-4,j,k)+a40(4)*ww(nx-3,j,k) &
                          +a40(3)*ww(nx-2,j,k)+a40(2)*ww(nx-1,j,k) &
                          +a40(1)*ww(nx  ,j,k) )*y_eta_v(nx,j)*ijacob_v(nx,j)
             dpx(   j,k)=( a40(5)*prs(nx-4,j,k)+a40(4)*prs(nx-3,j,k) &
                          +a40(3)*prs(nx-2,j,k)+a40(2)*prs(nx-1,j,k) &
                          +a40(1)*prs(nx  ,j,k) )*y_eta_v(nx,j)*ijacob_v(nx,j)
             drx(   j,k)=( a40(5)*rho_n(nx-4,j,k)+a40(4)*rho_n(nx-3,j,k) &
                          +a40(3)*rho_n(nx-2,j,k)+a40(2)*rho_n(nx-1,j,k) &
                          +a40(1)*rho_n(nx  ,j,k) )*y_eta_v(nx,j)*ijacob_v(nx,j)
          enddo
       enddo
    else
       do k=1,nz
          do j=1,ny
             dux(nx,j,k)=( a40(5)*uu(nx-4,j,k)+a40(4)*uu(nx-3,j,k) &
                          +a40(3)*uu(nx-2,j,k)+a40(2)*uu(nx-1,j,k) &
                          +a40(1)*uu(nx  ,j,k) )*idx_v(nx)
             dvx(nx,j,k)=( a40(5)*vv(nx-4,j,k)+a40(4)*vv(nx-3,j,k) &
                          +a40(3)*vv(nx-2,j,k)+a40(2)*vv(nx-1,j,k) &
                          +a40(1)*vv(nx  ,j,k) )*idx_v(nx)
             dwx(nx,j,k)=( a40(5)*ww(nx-4,j,k)+a40(4)*ww(nx-3,j,k) &
                          +a40(3)*ww(nx-2,j,k)+a40(2)*ww(nx-1,j,k) &
                          +a40(1)*ww(nx  ,j,k) )*idx_v(nx)
             dpx(   j,k)=( a40(5)*prs(nx-4,j,k)+a40(4)*prs(nx-3,j,k) &
                          +a40(3)*prs(nx-2,j,k)+a40(2)*prs(nx-1,j,k) &
                          +a40(1)*prs(nx  ,j,k) )*idx_v(nx)
             drx(   j,k)=( a40(5)*rho_n(nx-4,j,k)+a40(4)*rho_n(nx-3,j,k) &
                          +a40(3)*rho_n(nx-2,j,k)+a40(2)*rho_n(nx-1,j,k) &
                          +a40(1)*rho_n(nx  ,j,k) )*idx_v(nx)
          enddo
       enddo
    endif

    ! Update normal direction (LODI characteristics)
    ! =======================
    do k=1,nz
       do j=1,ny

          ! Characteristic velocities
          ! =========================
          cj=c_(nx,j,k)
          cj2=cj*cj
          upc=uu(nx,j,k)+cj
          umc=uu(nx,j,k)-cj

          ! Compute characteristic amplitudes
          ! =================================
          if (umc.le.0.0_wp) then
             !L1=0.0_wp
             !L1 = -0.28_wp*(1.0_wp-Mach**2)*cj/xg(ngx)*(prs(i,j,k)-p_ref)
             L1 = 0.1*(prs(nx,j,k)-0.528*p_ref)
          else
             L1=umc*(dpx(j,k)-rho_n(nx,j,k)*cj*dux(nx,j,k))
          endif
          if (upc.le.0.0_wp) then
             L5=0.0_wp
          else
             L5=upc*(dpx(j,k)+rho_n(nx,j,k)*cj*dux(nx,j,k))
          endif
          L2=uu(nx,j,k)*(dpx(j,k)-cj2*drx(j,k))
          L3=uu(nx,j,k)*dvx(nx,j,k)
          L4=uu(nx,j,k)*dwx(nx,j,k)

          ! Time derivatives for primitive variables
          ! ========================================
          pt(j,k)= 0.5_wp*(L5+L1)
          rt(j,k)= -(pt(j,k)+L2)/cj2
          ut(j,k)= 0.5_wp*(L5-L1)/(rho_n(nx,j,k)*cj)
          vt(j,k)= L3
          wt(j,k)= L4

          ! Update fluxes at each RK step
          ! =============================
          cp = cpcalc_tro(Tmp(nx,j,k),rho_n(nx,j,k))
          av = avcalc_tro(Tmp(nx,j,k),rho_n(nx,j,k))
          Krho(nx,j,k)  = rt(j,k)
          Krhou(nx,j,k) = rho_n(nx,j,k)*ut(j,k) + uu(nx,j,k)*rt(j,k)
          Krhov(nx,j,k) = rho_n(nx,j,k)*vt(j,k) + vv(nx,j,k)*rt(j,k)
          Krhow(nx,j,k) = rho_n(nx,j,k)*wt(j,k) + ww(nx,j,k)*rt(j,k)
          Krhoe(nx,j,k) = cp/av*(pt(j,k)/c_(nx,j,k)**2 - rt(j,k)) &
               + (rhoe_n(nx,j,k)+prs(nx,j,k))/rho_n(nx,j,k)*rt(j,k) &
               + rho_n(nx,j,k)*(uu(nx,j,k)*ut(j,k)+vv(nx,j,k)*vt(j,k)+ww(nx,j,k)*wt(j,k))
       enddo
    enddo

!!$    ! Flux derivatives along y
!!$    ! ========================
!!$
!!$    if ((coord(2)==0).and.(BC_face(2,1)%sort==-3)) then
!!$       ! right-bottom edge
!!$       ! =================
!!$       j=1
!!$
!!$       ! Compute derivatives along y
!!$       ! ===========================
!!$       if (is_curv) then
!!$          do k=1,nz
!!$             duy(i,j,k)= ( a04(1)*uu(i,j  ,k)+a04(2)*uu(i,j+1,k) &
!!$                  +a04(3)*uu(i,j+2,k)+a04(4)*uu(i,j+3,k) &
!!$                  +a04(5)*uu(i,j+4,k) )*x_ksi_v(i,j)*ijacob_v(i,j)
!!$             dvy(i,j,k)= ( a04(1)*vv(i,j  ,k)+a04(2)*vv(i,j+1,k) &
!!$                  +a04(3)*vv(i,j+2,k)+a04(4)*vv(i,j+3,k) &
!!$                  +a04(5)*vv(i,j+4,k) )*x_ksi_v(i,j)*ijacob_v(i,j)
!!$             dwy(i,j,k)= ( a04(1)*ww(i,j  ,k)+a04(2)*ww(i,j+1,k) &
!!$                  +a04(3)*ww(i,j+2,k)+a04(4)*ww(i,j+3,k) &
!!$                  +a04(5)*ww(i,j+4,k) )*x_ksi_v(i,j)*ijacob_v(i,j)
!!$             dpy(k)    = ( a04(1)*prs(i,j  ,k)+a04(2)*prs(i,j+1,k) &
!!$                  +a04(3)*prs(i,j+2,k)+a04(4)*prs(i,j+3,k) &
!!$                  +a04(5)*prs(i,j+4,k) )*x_ksi_v(i,j)*ijacob_v(i,j)
!!$             dry(k)    = ( a04(1)*rho_n(i,j  ,k)+a04(2)*rho_n(i,j+1,k) &
!!$                  +a04(3)*rho_n(i,j+2,k)+a04(4)*rho_n(i,j+3,k) &
!!$                  +a04(5)*rho_n(i,j+4,k) )*x_ksi_v(i,j)*ijacob_v(i,j)
!!$          enddo
!!$       else
!!$          do k=1,nz
!!$             duy(i,j,k)= ( a04(1)*uu(i,j  ,k)+a04(2)*uu(i,j+1,k) &
!!$                  +a04(3)*uu(i,j+2,k)+a04(4)*uu(i,j+3,k) &
!!$                  +a04(5)*uu(i,j+4,k) )*idy_v(j)
!!$             dvy(i,j,k)= ( a04(1)*vv(i,j  ,k)+a04(2)*vv(i,j+1,k) &
!!$                  +a04(3)*vv(i,j+2,k)+a04(4)*vv(i,j+3,k) &
!!$                  +a04(5)*vv(i,j+4,k) )*idy_v(j)
!!$             dwy(i,j,k)= ( a04(1)*ww(i,j  ,k)+a04(2)*ww(i,j+1,k) &
!!$                  +a04(3)*ww(i,j+2,k)+a04(4)*ww(i,j+3,k) &
!!$                  +a04(5)*ww(i,j+4,k) )*idy_v(j)
!!$             dpy(k)    = ( a04(1)*prs(i,j  ,k)+a04(2)*prs(i,j+1,k) &
!!$                  +a04(3)*prs(i,j+2,k)+a04(4)*prs(i,j+3,k) &
!!$                  +a04(5)*prs(i,j+4,k) )*idy_v(j)
!!$             dry(k)    = ( a04(1)*rho_n(i,j  ,k)+a04(2)*rho_n(i,j+1,k) &
!!$                  +a04(3)*rho_n(i,j+2,k)+a04(4)*rho_n(i,j+3,k) &
!!$                  +a04(5)*rho_n(i,j+4,k) )*idy_v(j)
!!$          enddo
!!$       endif
!!$
!!$       ! Update sound velocity (necessary if first runge-kutta step ?)
!!$       ! =====================
!!$       do k=1,nz
!!$          c_(i,j,k) = (c2calc_tro(Tmp(i,j,k),rho_n(i,j,k)))**0.5
!!$       enddo
!!$
!!$       do k=1,nz
!!$          ! Compute characteristic amplitudes
!!$          ! =================================
!!$          L1 = (vv(i,j,k)-c_(i,j,k))*(dpy(k)-rho_n(i,j,k)*c_(i,j,k)*dvy(i,j,k))
!!$          L5 = 0.0_wp
!!$          L1 = 0.0_wp
!!$          if (vv(i,j,k).ge.0.0_wp) then
!!$          ! if (v0_sud(i,k).ge.0.0_wp) then
!!$             L2 = 0.0_wp
!!$             L3 = 0.0_wp
!!$             L4 = 0.0_wp
!!$             if (vv(i,j,k).ge.c_(i,j,k)) L1 = 0.0_wp
!!$          else
!!$             L2 = vv(i,j,k)*duy(i,j,k)
!!$             L3 = vv(i,j,k)*(dpy(k)-(c_(i,j,k)**2)*dry(k))
!!$             L4 = vv(i,j,k)*dwy(i,j,k)
!!$             if (vv(i,j,k).le.-c_(i,j,k)) L5=(vv(i,j,k)+c_(i,j,k))*(dpy(k)+rho_n(i,j,k)*c_(i,j,k)*dvy(i,j,k))
!!$          endif
!!$
!!$          ! Time derivatives for primitive variables
!!$          ! ========================================
!!$          pt(j,k) = 0.5_wp*(L5+L1)
!!$          rt(j,k) = -(pt(j,k)+L3)/(c_(i,j,k)**2)
!!$          ut(j,k) = L2
!!$          vt(j,k) = (L5-L1)/(2.0_wp*rho_n(i,j,k)*c_(i,j,k))
!!$          wt(j,k) = L4
!!$       enddo
!!$
!!$       do k=1,nz
!!$          ! Update fluxes at each RK step
!!$          ! =============================
!!$          cp = cpcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
!!$          av = avcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
!!$          Krho(i,j,k)  = Krho(i,j,k)  + rt(j,k)
!!$          Krhou(i,j,k) = Krhou(i,j,k) + rho_n(i,j,k)*ut(j,k) + uu(i,j,k)*rt(j,k)
!!$          Krhov(i,j,k) = Krhov(i,j,k) + rho_n(i,j,k)*vt(j,k) + vv(i,j,k)*rt(j,k)
!!$          Krhow(i,j,k) = Krhow(i,j,k) + rho_n(i,j,k)*wt(j,k) + ww(i,j,k)*rt(j,k)
!!$          Krhoe(i,j,k) = Krhoe(i,j,k) +  cp/av*(pt(j,k)/c_(i,j,k)**2 - rt(j,k)) &
!!$               + (rhoe_n(i,j,k)+prs(i,j,k))/rho_n(i,j,k)*rt(j,k) &
!!$               + rho_n(i,j,k)*(uu(i,j,k)*ut(j,k)+vv(i,j,k)*vt(j,k)+ww(i,j,k)*wt(j,k))
!!$       enddo
!!$
!!$    endif
!!$
!!$    if ((coord(2)==ndomy-1).and.(BC_face(2,2)%sort==-3)) then
!!$       ! right-top edge
!!$       ! ==============
!!$       j=ny
!!$
!!$       ! Compute derivatives along y
!!$       ! ===========================
!!$       if (is_curv) then
!!$          do k=1,nz
!!$             duy(i,j,k)= ( a40(5)*uu(i,j-4,k)+a40(4)*uu(i,j-3,k) &
!!$                  +a40(3)*uu(i,j-2,k)+a40(2)*uu(i,j-1,k) &
!!$                  +a40(1)*uu(i,j  ,k) )*x_ksi_v(i,j)*ijacob_v(i,j)
!!$             dvy(i,j,k)= ( a40(5)*vv(i,j-4,k)+a40(4)*vv(i,j-3,k) &
!!$                  +a40(3)*vv(i,j-2,k)+a40(2)*vv(i,j-1,k) &
!!$                  +a40(1)*vv(i,j  ,k) )*x_ksi_v(i,j)*ijacob_v(i,j)
!!$             dwy(i,j,k)= ( a40(5)*ww(i,j-4,k)+a40(4)*ww(i,j-3,k) &
!!$                  +a40(3)*ww(i,j-2,k)+a40(2)*ww(i,j-1,k) &
!!$                  +a40(1)*ww(i,j  ,k) )*x_ksi_v(i,j)*ijacob_v(i,j)
!!$             dpy(k)    = ( a40(5)*prs(i,j-4,k)+a40(4)*prs(i,j-3,k) &
!!$                  +a40(3)*prs(i,j-2,k)+a40(2)*prs(i,j-1,k) &
!!$                  +a40(1)*prs(i,j  ,k) )*x_ksi_v(i,j)*ijacob_v(i,j)
!!$             dry(k)    = ( a40(5)*rho_n(i,j-4,k)+a40(4)*rho_n(i,j-3,k) &
!!$                  +a40(3)*rho_n(i,j-2,k)+a40(2)*rho_n(i,j-1,k) &
!!$                  +a40(1)*rho_n(i,j  ,k) )*x_ksi_v(i,j)*ijacob_v(i,j)
!!$          enddo
!!$       else
!!$          do k=1,nz
!!$             duy(i,j,k)= ( a40(5)*uu(i,j-4,k)+a40(4)*uu(i,j-3,k) &
!!$                  +a40(3)*uu(i,j-2,k)+a40(2)*uu(i,j-1,k) &
!!$                  +a40(1)*uu(i,j  ,k) )*idy_v(j)
!!$             dvy(i,j,k)= ( a40(5)*vv(i,j-4,k)+a40(4)*vv(i,j-3,k) &
!!$                  +a40(3)*vv(i,j-2,k)+a40(2)*vv(i,j-1,k) &
!!$                  +a40(1)*vv(i,j  ,k) )*idy_v(j)
!!$             dwy(i,j,k)= ( a40(5)*ww(i,j-4,k)+a40(4)*ww(i,j-3,k) &
!!$                  +a40(3)*ww(i,j-2,k)+a40(2)*ww(i,j-1,k) &
!!$                  +a40(1)*ww(i,j  ,k) )*idy_v(j)
!!$             dpy(k)    = ( a40(5)*prs(i,j-4,k)+a40(4)*prs(i,j-3,k) &
!!$                  +a40(3)*prs(i,j-2,k)+a40(2)*prs(i,j-1,k) &
!!$                  +a40(1)*prs(i,j  ,k) )*idy_v(j)
!!$             dry(k)    = ( a40(5)*rho_n(i,j-4,k)+a40(4)*rho_n(i,j-3,k) &
!!$                  +a40(3)*rho_n(i,j-2,k)+a40(2)*rho_n(i,j-1,k) &
!!$                  +a40(1)*rho_n(i,j  ,k) )*idy_v(j)
!!$          enddo
!!$       endif
!!$
!!$       ! Update sound velocity (necessary if first runge-kutta step ?)
!!$       ! =====================
!!$       do k=1,nz
!!$          c_(i,j,k) = (c2calc_tro(Tmp(i,j,k),rho_n(i,j,k)))**0.5
!!$       enddo
!!$
!!$       do k=1,nz
!!$          ! Compute characteristic amplitudes
!!$          ! =================================
!!$          L1 = 0.0_wp
!!$          L5 = (vv(i,j,k)+c_(i,j,k))*(dpy(k)+rho_n(i,j,k)*c_(i,j,k)*dvy(i,j,k))
!!$          if (vv(i,j,k).le.0.0_wp) then
!!$          ! if (v0_nord(i,k).le.0.0_wp) then
!!$             L2 = 0.0_wp
!!$             L3 = 0.0_wp
!!$             L4 = 0.0_wp
!!$             if (vv(i,j,k).le.-c_(i,j,k)) L5 = 0.0_wp
!!$          else
!!$             L2 = vv(i,j,k)*duy(i,j,k)
!!$             L3 = vv(i,j,k)*(dpy(k)-(c_(i,j,k)**2)*dry(k))
!!$             L4 = vv(i,j,k)*dwy(i,j,k)
!!$             if (vv(i,j,k).ge.c_(i,j,k)) L1 = (vv(i,j,k)-c_(i,j,k))*(dpy(k)-rho_n(i,j,k)*c_(i,j,k)*dvy(i,j,k))
!!$          endif
!!$
!!$          ! Time derivatives for primitive variables
!!$          ! ========================================
!!$          pt(j,k) = 0.5_wp*(L5+L1)
!!$          rt(j,k) = -(pt(j,k)+L3)/(c_(i,j,k)**2)
!!$          ut(j,k) = L2
!!$          vt(j,k) = (L5-L1)/(2.0_wp*rho_n(i,j,k)*c_(i,j,k))
!!$          wt(j,k) = L4
!!$       enddo
!!$
!!$       do k=1,nz
!!$          ! Update fluxes at each RK step
!!$          ! =============================
!!$          cp = cpcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
!!$          av = avcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
!!$          Krho(i,j,k)  = Krho(i,j,k)  + rt(j,k)
!!$          Krhou(i,j,k) = Krhou(i,j,k) + rho_n(i,j,k)*ut(j,k) + uu(i,j,k)*rt(j,k)
!!$          Krhov(i,j,k) = Krhov(i,j,k) + rho_n(i,j,k)*vt(j,k) + vv(i,j,k)*rt(j,k)
!!$          Krhow(i,j,k) = Krhow(i,j,k) + rho_n(i,j,k)*wt(j,k) + ww(i,j,k)*rt(j,k)
!!$          Krhoe(i,j,k) = Krhoe(i,j,k) +  cp/av*(pt(j,k)/c_(i,j,k)**2 - rt(j,k)) &
!!$               + (rhoe_n(i,j,k)+prs(i,j,k))/rho_n(i,j,k)*rt(j,k) &
!!$               + rho_n(i,j,k)*(uu(i,j,k)*ut(j,k)+vv(i,j,k)*vt(j,k)+ww(i,j,k)*wt(j,k))
!!$       enddo
!!$    endif
    ! ---------------------------------------------------------------------------

  end subroutine bc_charac_imax

  !===============================================================================
  subroutine bc_charac_jmin
  !===============================================================================
    !> Characteristic BC: boundary condition at jmin (bottom) - Cartesian version -
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer  :: i,j,k
    real(wp), dimension(nx,nz) :: rt,ut,vt,wt,pt
    real(wp), dimension(nx,nz) :: dry,dpy
    ! ---------------------------------------------------------------------------
    ! ! Temporary
    ! real(wp) :: inn

    ! ! Test : temporary
    ! inn = 1.0_wp/dble(10000)
    ! v0_sud(1:nx,1:nz)  = ((10000-1)*v0_sud(1:ny,1:nz) + vv(1:nx,1,1:nz))*inn

    ! Index of bottom boundary
    ! ========================
    j=1

    ! Wall condition applied (different from Tam&Dong)
    ! ======================
    ! Wall BC at imin
    ! ---------------
    if (is_bc_wall(1,1)) then
        Krho(1,j,ndz_e:nfz_e)=0.0_wp
       Krhou(1,j,ndz_e:nfz_e)=0.0_wp
       Krhov(1,j,ndz_e:nfz_e)=0.0_wp
       Krhow(1,j,ndz_e:nfz_e)=0.0_wp
       Krhoe(1,j,ndz_e:nfz_e)=0.0_wp
    endif
    ! Wall BC at imax
    ! ---------------
    if (is_bc_wall(1,2)) then
        Krho(nx,j,ndz_e:nfz_e)=0.0_wp
       Krhou(nx,j,ndz_e:nfz_e)=0.0_wp
       Krhov(nx,j,ndz_e:nfz_e)=0.0_wp
       Krhow(nx,j,ndz_e:nfz_e)=0.0_wp
       Krhoe(nx,j,ndz_e:nfz_e)=0.0_wp
    endif
    ! Wall BC at kmin
    ! ---------------
    if (is_bc_wall(3,1)) then
        Krho(ndx_e:nfx_e,j,1)=0.0_wp
       Krhou(ndx_e:nfx_e,j,1)=0.0_wp
       Krhov(ndx_e:nfx_e,j,1)=0.0_wp
       Krhow(ndx_e:nfx_e,j,1)=0.0_wp
       Krhoe(ndx_e:nfx_e,j,1)=0.0_wp
    endif
    ! Wall BC at kmax
    ! ---------------
    if (is_bc_wall(3,2)) then
        Krho(ndx_e:nfx_e,j,nz)=0.0_wp
       Krhou(ndx_e:nfx_e,j,nz)=0.0_wp
       Krhov(ndx_e:nfx_e,j,nz)=0.0_wp
       Krhow(ndx_e:nfx_e,j,nz)=0.0_wp
       Krhoe(ndx_e:nfx_e,j,nz)=0.0_wp
    endif

    ! Compute derivatives along y
    ! ===========================
    if (is_curv) then
       do k=1,nz
          ! do i=ndx3,nfx3
          ! do i=ndx_e,nfx_e
          do i=ndx_c,nfx_c
             duy(i,j,k)= ( a04(1)*uu(i,j  ,k)+a04(2)*uu(i,j+1,k) &
                  +a04(3)*uu(i,j+2,k)+a04(4)*uu(i,j+3,k) &
                  +a04(5)*uu(i,j+4,k) )*x_ksi_v(i,j)*ijacob_v(i,j)
             dvy(i,j,k)= ( a04(1)*vv(i,j  ,k)+a04(2)*vv(i,j+1,k) &
                  +a04(3)*vv(i,j+2,k)+a04(4)*vv(i,j+3,k) &
                  +a04(5)*vv(i,j+4,k) )*x_ksi_v(i,j)*ijacob_v(i,j)
             dwy(i,j,k)= ( a04(1)*ww(i,j  ,k)+a04(2)*ww(i,j+1,k) &
                  +a04(3)*ww(i,j+2,k)+a04(4)*ww(i,j+3,k) &
                  +a04(5)*ww(i,j+4,k) )*x_ksi_v(i,j)*ijacob_v(i,j)
             dpy(i,k)  = ( a04(1)*prs(i,j  ,k)+a04(2)*prs(i,j+1,k) &
                  +a04(3)*prs(i,j+2,k)+a04(4)*prs(i,j+3,k) &
                  +a04(5)*prs(i,j+4,k) )*x_ksi_v(i,j)*ijacob_v(i,j)
             dry(i,k)  = ( a04(1)*rho_n(i,j  ,k)+a04(2)*rho_n(i,j+1,k) &
                  +a04(3)*rho_n(i,j+2,k)+a04(4)*rho_n(i,j+3,k) &
                  +a04(5)*rho_n(i,j+4,k) )*x_ksi_v(i,j)*ijacob_v(i,j)
          enddo
       enddo
    else
       do k=1,nz
          ! do i=ndx3,nfx3
          ! do i=ndx_e,nfx_e
          do i=ndx_c,nfx_c
             duy(i,j,k)= ( a04(1)*uu(i,j  ,k)+a04(2)*uu(i,j+1,k) &
                  +a04(3)*uu(i,j+2,k)+a04(4)*uu(i,j+3,k) &
                  +a04(5)*uu(i,j+4,k) )*idy_v(j)
             dvy(i,j,k)= ( a04(1)*vv(i,j  ,k)+a04(2)*vv(i,j+1,k) &
                  +a04(3)*vv(i,j+2,k)+a04(4)*vv(i,j+3,k) &
                  +a04(5)*vv(i,j+4,k) )*idy_v(j)
             dwy(i,j,k)= ( a04(1)*ww(i,j  ,k)+a04(2)*ww(i,j+1,k) &
                  +a04(3)*ww(i,j+2,k)+a04(4)*ww(i,j+3,k) &
                  +a04(5)*ww(i,j+4,k) )*idy_v(j)
             dpy(i,k)  = ( a04(1)*prs(i,j  ,k)+a04(2)*prs(i,j+1,k) &
                  +a04(3)*prs(i,j+2,k)+a04(4)*prs(i,j+3,k) &
                  +a04(5)*prs(i,j+4,k) )*idy_v(j)
             dry(i,k)  = ( a04(1)*rho_n(i,j  ,k)+a04(2)*rho_n(i,j+1,k) &
                  +a04(3)*rho_n(i,j+2,k)+a04(4)*rho_n(i,j+3,k) &
                  +a04(5)*rho_n(i,j+4,k) )*idy_v(j)
          enddo
       enddo
    endif


    ! Update normal direction (LODI characteristics)
    ! =======================
    do k=1,nz
       do i=ndx_c,nfx_c

          ! Update sound velocity (necessary if first runge-kutta step ?)
          ! =====================
          c_(i,j,k) = (c2calc_tro(Tmp(i,j,k),rho_n(i,j,k)))**0.5

          ! Compute characteristic amplitudes
          ! =================================
          L1 = (vv(i,j,k)-c_(i,j,k))*(dpy(i,k)-rho_n(i,j,k)*c_(i,j,k)*dvy(i,j,k))
          L5 = 0.0_wp
          if (vv(i,j,k).ge.0.0_wp) then
          ! if (v0_sud(i,k).ge.0.0_wp) then
             L2 = 0.0_wp
             L3 = 0.0_wp
             L4 = 0.0_wp
             if (vv(i,j,k).ge.c_(i,j,k)) L1 = 0.0_wp
          else
             L2 = vv(i,j,k)*duy(i,j,k)
             L3 = vv(i,j,k)*(dpy(i,k)-(c_(i,j,k)**2)*dry(i,k))
             L4 = vv(i,j,k)*dwy(i,j,k)
             if (vv(i,j,k).le.-c_(i,j,k)) L5=(vv(i,j,k)+c_(i,j,k))*(dpy(i,k)+rho_n(i,j,k)*c_(i,j,k)*dvy(i,j,k))
          endif

          ! Time derivatives for primitive variables
          ! ========================================
          pt(i,k) = 0.5_wp*(L5+L1)
          rt(i,k) = -(pt(i,k)+L3)/(c_(i,j,k)**2)
          ut(i,k) = L2
          vt(i,k) = (L5-L1)/(2.0_wp*rho_n(i,j,k)*c_(i,j,k))
          wt(i,k) = L4

          ! Update fluxes at each RK step
          ! =============================
          cp = cpcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
          av = avcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
          Krho(i,j,k)  = rt(i,k)
          Krhou(i,j,k) = rho_n(i,j,k)*ut(i,k) + uu(i,j,k)*rt(i,k)
          Krhov(i,j,k) = rho_n(i,j,k)*vt(i,k) + vv(i,j,k)*rt(i,k)
          Krhow(i,j,k) = rho_n(i,j,k)*wt(i,k) + ww(i,j,k)*rt(i,k)
          Krhoe(i,j,k) = cp/av*(pt(i,k)/c_(i,j,k)**2 - rt(i,k)) &
               + (rhoe_n(i,j,k)+prs(i,j,k))/rho_n(i,j,k)*rt(i,k) &
               + rho_n(i,j,k)*(uu(i,j,k)*ut(i,k)+vv(i,j,k)*vt(i,k)+ww(i,j,k)*wt(i,k))
       enddo
    enddo

  end subroutine bc_charac_jmin

  !===============================================================================
  subroutine bc_charac_jmax
  !===============================================================================
    !> Characteristic BC: boundary condition at jmax (top) - Cartesian version -
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer  :: i,j,k
    real(wp), dimension(nx,nz) :: rt,ut,vt,wt,pt
    real(wp), dimension(nx,nz) :: dry,dpy
    ! ---------------------------------------------------------------------------
    ! ! Temporary
    ! real(wp) :: inn

    ! ! Test : temporary
    ! inn = 1.0_wp/dble(10000)
    ! v0_nord(1:nx,1:nz)  = ((10000-1)*v0_nord(1:ny,1:nz) + vv(1:nx,ny,1:nz))*inn

    ! Index of top boundary
    ! ======================
    j=ny

    ! Wall condition applied (different from Tam&Dong)
    ! ======================
    ! Wall BC at imin
    ! ---------------
    if (is_bc_wall(1,1)) then
        Krho(1,j,ndz_e:nfz_e)=0.0_wp
       Krhou(1,j,ndz_e:nfz_e)=0.0_wp
       Krhov(1,j,ndz_e:nfz_e)=0.0_wp
       Krhow(1,j,ndz_e:nfz_e)=0.0_wp
       Krhoe(1,j,ndz_e:nfz_e)=0.0_wp
    endif
    ! Wall BC at imax
    ! ---------------
    if (is_bc_wall(1,2)) then
        Krho(nx,j,ndz_e:nfz_e)=0.0_wp
       Krhou(nx,j,ndz_e:nfz_e)=0.0_wp
       Krhov(nx,j,ndz_e:nfz_e)=0.0_wp
       Krhow(nx,j,ndz_e:nfz_e)=0.0_wp
       Krhoe(nx,j,ndz_e:nfz_e)=0.0_wp
    endif
    ! Wall BC at kmin
    ! ---------------
    if (is_bc_wall(3,1)) then
        Krho(ndx_e:nfx_e,j,1)=0.0_wp
       Krhou(ndx_e:nfx_e,j,1)=0.0_wp
       Krhov(ndx_e:nfx_e,j,1)=0.0_wp
       Krhow(ndx_e:nfx_e,j,1)=0.0_wp
       Krhoe(ndx_e:nfx_e,j,1)=0.0_wp
    endif
    ! Wall BC at kmax
    ! ---------------
    if (is_bc_wall(3,2)) then
        Krho(ndx_e:nfx_e,j,nz)=0.0_wp
       Krhou(ndx_e:nfx_e,j,nz)=0.0_wp
       Krhov(ndx_e:nfx_e,j,nz)=0.0_wp
       Krhow(ndx_e:nfx_e,j,nz)=0.0_wp
       Krhoe(ndx_e:nfx_e,j,nz)=0.0_wp
    endif

    ! Compute derivatives along y
    ! ===========================
    if (is_curv) then
       do k=1,nz
          do i=ndx_c,nfx_c
             duy(i,j,k)= ( a40(5)*uu(i,j-4,k)+a40(4)*uu(i,j-3,k) &
                          +a40(3)*uu(i,j-2,k)+a40(2)*uu(i,j-1,k) &
                          +a40(1)*uu(i,j  ,k) )*x_ksi_v(i,j)*ijacob_v(i,j)
             dvy(i,j,k)= ( a40(5)*vv(i,j-4,k)+a40(4)*vv(i,j-3,k) &
                          +a40(3)*vv(i,j-2,k)+a40(2)*vv(i,j-1,k) &
                          +a40(1)*vv(i,j  ,k) )*x_ksi_v(i,j)*ijacob_v(i,j)
             dwy(i,j,k)= ( a40(5)*ww(i,j-4,k)+a40(4)*ww(i,j-3,k) &
                          +a40(3)*ww(i,j-2,k)+a40(2)*ww(i,j-1,k) &
                          +a40(1)*ww(i,j  ,k) )*x_ksi_v(i,j)*ijacob_v(i,j)
             dpy(i,k)  = ( a40(5)*prs(i,j-4,k)+a40(4)*prs(i,j-3,k) &
                          +a40(3)*prs(i,j-2,k)+a40(2)*prs(i,j-1,k) &
                          +a40(1)*prs(i,j  ,k) )*x_ksi_v(i,j)*ijacob_v(i,j)
             dry(i,k)  = ( a40(5)*rho_n(i,j-4,k)+a40(4)*rho_n(i,j-3,k) &
                          +a40(3)*rho_n(i,j-2,k)+a40(2)*rho_n(i,j-1,k) &
                          +a40(1)*rho_n(i,j  ,k) )*x_ksi_v(i,j)*ijacob_v(i,j)
          enddo
       enddo
    else
       do k=1,nz
          do i=ndx_c,nfx_c
             duy(i,j,k)= ( a40(5)*uu(i,j-4,k)+a40(4)*uu(i,j-3,k) &
                          +a40(3)*uu(i,j-2,k)+a40(2)*uu(i,j-1,k) &
                          +a40(1)*uu(i,j  ,k) )*idy_v(j)
             dvy(i,j,k)= ( a40(5)*vv(i,j-4,k)+a40(4)*vv(i,j-3,k) &
                          +a40(3)*vv(i,j-2,k)+a40(2)*vv(i,j-1,k) &
                          +a40(1)*vv(i,j  ,k) )*idy_v(j)
             dwy(i,j,k)= ( a40(5)*ww(i,j-4,k)+a40(4)*ww(i,j-3,k) &
                          +a40(3)*ww(i,j-2,k)+a40(2)*ww(i,j-1,k) &
                          +a40(1)*ww(i,j  ,k) )*idy_v(j)
             dpy(i,k)  = ( a40(5)*prs(i,j-4,k)+a40(4)*prs(i,j-3,k) &
                          +a40(3)*prs(i,j-2,k)+a40(2)*prs(i,j-1,k) &
                          +a40(1)*prs(i,j  ,k) )*idy_v(j)
             dry(i,k)  = ( a40(5)*rho_n(i,j-4,k)+a40(4)*rho_n(i,j-3,k) &
                          +a40(3)*rho_n(i,j-2,k)+a40(2)*rho_n(i,j-1,k) &
                          +a40(1)*rho_n(i,j  ,k) )*idy_v(j)
          enddo
       enddo
    endif

    ! Update normal direction (LODI characteristics)
    ! =======================
    do k=1,nz
        do i=ndx_c,nfx_c

          ! Update sound velocity (necessary if first runge-kutta step ?)
          ! =====================
          c_(i,j,k) = (c2calc_tro(Tmp(i,j,k),rho_n(i,j,k)))**0.5

          ! Compute characteristic amplitudes
          ! =================================
          L1 = 0.0_wp
          L5 = (vv(i,j,k)+c_(i,j,k))*(dpy(i,k)+rho_n(i,j,k)*c_(i,j,k)*dvy(i,j,k))
          if (vv(i,j,k).le.0.0_wp) then
          ! if (v0_nord(i,k).le.0.0_wp) then
             L2 = 0.0_wp
             L3 = 0.0_wp
             L4 = 0.0_wp
             if (vv(i,j,k).le.-c_(i,j,k)) L5=0.0_wp
          else
             L2 = vv(i,j,k)*duy(i,j,k)
             L3 = vv(i,j,k)*(dpy(i,k)-(c_(i,j,k)**2)*dry(i,k))
             L4 = vv(i,j,k)*dwy(i,j,k)
             if (vv(i,j,k).ge.c_(i,j,k)) L1=(vv(i,j,k)-c_(i,j,k))*(dpy(i,k)-rho_n(i,j,k)*c_(i,j,k)*dvy(i,j,k))
          endif

          ! Time derivatives for primitive variables
          ! ========================================
          pt(i,k) = 0.5_wp*(L5+L1)
          !rt(i,k) = -(pt(i,k)+L3)/(c_(i,j,k)**2)
          rt(i,k) = (pt(i,k)+L3)/(c_(i,j,k)**2)
          ut(i,k) = L2
          vt(i,k) = (L5-L1)/(2.0_wp*rho_n(i,j,k)*c_(i,j,k))
          wt(i,k) = L4

          ! Update fluxes at each RK step
          ! =============================
          cp = cpcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
          av = avcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
          Krho(i,j,k)  = rt(i,k)
          Krhou(i,j,k) = rho_n(i,j,k)*ut(i,k) + uu(i,j,k)*rt(i,k)
          Krhov(i,j,k) = rho_n(i,j,k)*vt(i,k) + vv(i,j,k)*rt(i,k)
          Krhow(i,j,k) = rho_n(i,j,k)*wt(i,k) + ww(i,j,k)*rt(i,k)
          Krhoe(i,j,k) = cp/av*(pt(i,k)/c_(i,j,k)**2 - rt(i,k)) &
               + (rhoe_n(i,j,k)+prs(i,j,k))/rho_n(i,j,k)*rt(i,k) &
               + rho_n(i,j,k)*(uu(i,j,k)*ut(i,k)+vv(i,j,k)*vt(i,k)+ww(i,j,k)*wt(i,k))
       enddo
    enddo

  end subroutine bc_charac_jmax

  !===============================================================================
  subroutine bc_charac_kmax
  !===============================================================================
    !> Characteristic BC: boundary condition at kmax (back) - Cartesian version -
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer  :: i,j,k
    ! ---------------------------------------------------------------------------

    ! Update normal direction (LODI characteristics)
    ! =======================
    k=1
    do j=1,ny
        do i=ndx_c,nfx_c
          Krho(i,j,k)  = 0.
          Krhou(i,j,k) = 0.
          Krhov(i,j,k) = 0.
          Krhow(i,j,k) = 0.
          Krhoe(i,j,k) = 0.
       enddo
    enddo

  end subroutine bc_charac_kmax

  !===============================================================================
  subroutine bc_supersonic_inlet
  !===============================================================================
    !> supersonic inlet: impose rhoin_ref,.. at imin (left)
  !===============================================================================
    use mod_eigenmode
    use mod_RFM
    use mod_time ! for: deltat,ck
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: j,k
    !-------------------------------------------------------------------------
    ! eigenmode disturbances
    real(wp) :: alpha,cp,c2,av
    real(wp) :: um,vm,wm,rm,tm ! disturbances
    real(wp) :: dumt,dvmt,dwmt,drmt,dpmt ! time derivatives of disturbances
    real(wp), dimension(ny,nz) :: ut_in,vt_in,wt_in
    ! ---------------------------------------------------------------------------

    ! time at RK stage
    ! ================
    alpha=ck(irk)*deltat

    if (is_eigenmode) then ! add eigenmodes

!!$       ! direct imposition of disturbances
!!$       ! =================================
!!$       do k=1,nz
!!$          do j=1,ny
!!$
!!$             ! add disturbances to primitives variables
!!$             ! ----------------------------------------
!!$             call eig_disturb1(1,j,k,um,vm,wm,rm,tm)
!!$             rho_n(1,j,k)= BC_face(1,1)%Uref(j,k,1) + rm*1.
!!$              uu(1,j,k)  = BC_face(1,1)%Uref(j,k,2) + um*1.
!!$              vv(1,j,k)  = BC_face(1,1)%Uref(j,k,3) + vm*1.
!!$              ww(1,j,k)  = BC_face(1,1)%Uref(j,k,4) + wm*1.
!!$             Tmp(1,j,k)  = BC_face(1,1)%Uref(j,k,6) + tm*1.
!!$
!!$             ! compute conservative variables
!!$             ! ------------------------------
!!$             rhou_n(1,j,k) = rho_n(1,j,k)*uu(1,j,k)
!!$             rhov_n(1,j,k) = rho_n(1,j,k)*vv(1,j,k)
!!$             rhow_n(1,j,k) = rho_n(1,j,k)*ww(1,j,k)
!!$             rhoe_n(1,j,k) = rho_n(1,j,k)*( ecalc_tro(Tmp(1,j,k),rho_n(1,j,k)) &
!!$                           + 0.5_wp*(uu(1,j,k)**2+vv(1,j,k)**2+ww(1,j,k)**2) )
!!$          enddo
!!$       enddo

       ! imposition of time derivatives of disturbances
       ! ==============================================
       do k=1,nz
          do j=1,ny

             ! compute alpha*(time derivatives of disturb)
             ! -------------------------------------------
             call eig_disturb1_dt(1,j,k,dumt,dvmt,dwmt,drmt,dpmt)

             ! update conservative variables
             ! -----------------------------
             rho_n(1,j,k)= BC_face(1,1)%Uref(j,k,1)
             rhou_n(1,j,k)= BC_face(1,1)%Uref(j,k,7) + alpha*rho_n(1,j,k)*dumt
             rhov_n(1,j,k)= BC_face(1,1)%Uref(j,k,8) + alpha*rho_n(1,j,k)*dvmt
             rhow_n(1,j,k)= BC_face(1,1)%Uref(j,k,9) + alpha*rho_n(1,j,k)*dwmt
             ! *****************************
             ! /!\ temporarily only for PFG
             ! *****************************
             rhoe_n(1,j,k)= BC_face(1,1)%Uref(j,k,10) + alpha*( dpmt*igm1 &
                  +rhou_n(1,j,k)*dumt+rhov_n(1,j,k)*dvmt+rhow_n(1,j,k)*dwmt )

!!$             drmt=0.
!!$             dpmt=0.
!!$
!!$             rho_n(1,j,k)= BC_face(1,1)%Uref(j,k,1) + alpha*drmt
!!$             rhou_n(1,j,k)= BC_face(1,1)%Uref(j,k,7) &
!!$                  + alpha*(BC_face(1,1)%Uref(j,k,1)*dumt+BC_face(1,1)%Uref(j,k,2)*drmt)
!!$             rhov_n(1,j,k)= BC_face(1,1)%Uref(j,k,8) &
!!$                  + alpha*(BC_face(1,1)%Uref(j,k,1)*dvmt+BC_face(1,1)%Uref(j,k,3)*drmt)
!!$             rhow_n(1,j,k)= BC_face(1,1)%Uref(j,k,9) &
!!$                  + alpha*(BC_face(1,1)%Uref(j,k,1)*dwmt+BC_face(1,1)%Uref(j,k,4)*drmt)
!!$             
!!$             cp = cpcalc_tro(Tmp(1,j,k),rho_n(1,j,k))
!!$             av = avcalc_tro(Tmp(1,j,k),rho_n(1,j,k))
!!$             c2 = c2calc_tro(Tmp(1,j,k),rho_n(1,j,k))
!!$          
!!$             rhoe_n(1,j,k)= BC_face(1,1)%Uref(j,k,10) + alpha*( cp/av*(dpmt/c2-drmt) &
!!$               + (BC_face(1,1)%Uref(j,k,10)+BC_face(1,1)%Uref(j,k,5))/BC_face(1,1)%Uref(j,k,1)*drmt &
!!$               + BC_face(1,1)%Uref(j,k,1)*(BC_face(1,1)%Uref(j,k,2)*dumt+BC_face(1,1)%Uref(j,k,3)*dvmt+BC_face(1,1)%Uref(j,k,4)*dwmt) )
          enddo
       enddo

    elseif (is_RFM) then ! add eigenmodes
       
       ! compute RFM disturbances
       ! ------------------------
       call disturb_inlet_RFM_charac(ut_in,vt_in,wt_in)
       !ut_in=0.0_wp
       !vt_in=0.0_wp
       !wt_in=0.0_wp

       do k=1,nz
          do j=1,ny

             ! update conservative variables
             ! -----------------------------
             rho_n(1,j,k)= BC_face(1,1)%Uref(j,k,1)
             rhou_n(1,j,k)= BC_face(1,1)%Uref(j,k,7) + alpha*rho_n(1,j,k)*ut_in(j,k)*20.
             rhov_n(1,j,k)= BC_face(1,1)%Uref(j,k,8) + alpha*rho_n(1,j,k)*vt_in(j,k)*20.
             rhow_n(1,j,k)= BC_face(1,1)%Uref(j,k,9) + alpha*rho_n(1,j,k)*wt_in(j,k)*20.
             ! *****************************
             ! /!\ temporarily only for PFG
             ! *****************************
             rhoe_n(1,j,k)= BC_face(1,1)%Uref(j,k,10) + alpha*( &
                  +rhou_n(1,j,k)*ut_in(j,k)+rhov_n(1,j,k)*vt_in(j,k)+rhow_n(1,j,k)*wt_in(j,k) )
          enddo
       enddo
       
!!$       ! direct imposition of disturbances
!!$       ! =================================
!!$       do k=1,nz
!!$          do j=1,ny-20
!!$
!!$             ! add disturbances to primitives variables
!!$             ! ----------------------------------------
!!$             rho_n(1,j,k)= BC_face(1,1)%Uref(j,k,1)
!!$              uu(1,j,k)  = BC_face(1,1)%Uref(j,k,2) + ut_in(j,k)
!!$              vv(1,j,k)  = BC_face(1,1)%Uref(j,k,3) + vt_in(j,k)
!!$              ww(1,j,k)  = BC_face(1,1)%Uref(j,k,4) + wt_in(j,k)
!!$             Tmp(1,j,k)  = BC_face(1,1)%Uref(j,k,6)
!!$
!!$             ! compute conservative variables
!!$             ! ------------------------------
!!$             rhou_n(1,j,k) = rho_n(1,j,k)*uu(1,j,k)
!!$             rhov_n(1,j,k) = rho_n(1,j,k)*vv(1,j,k)
!!$             rhow_n(1,j,k) = rho_n(1,j,k)*ww(1,j,k)
!!$             rhoe_n(1,j,k) = rho_n(1,j,k)*( ecalc_tro(Tmp(1,j,k),rho_n(1,j,k)) &
!!$                           + 0.5_wp*(uu(1,j,k)**2+vv(1,j,k)**2+ww(1,j,k)**2) )
!!$          enddo
!!$       enddo
    else

       do k=1,nz
          do j=1,ny
             rho_n(1,j,k)= BC_face(1,1)%Uref(j,k,1)
             rhou_n(1,j,k)= BC_face(1,1)%Uref(j,k,7)
             rhov_n(1,j,k)= BC_face(1,1)%Uref(j,k,8)
             rhow_n(1,j,k)= BC_face(1,1)%Uref(j,k,9)
             rhoe_n(1,j,k)= BC_face(1,1)%Uref(j,k,10)
          enddo
       enddo

    endif

  end subroutine bc_supersonic_inlet
  
end module mod_charac
