!===============================================================================
submodule (mod_TamDong2d) smod_TamDong2d_SBP4_edges
!===============================================================================
  !> author: XG
  !> date: February 2020 - modif January 2022
  !> Non-reflecting boundary conditions of Tam and Dong
  !> - 2D (periodic) Cartesian version - routines for edges
  !> /!\ dev only: version with SBP4 boundary schemes [not for regular use]
!=============================================================================== 

contains

  !===============================================================================
  module subroutine bc_TD2d_imin_jmin_SBP4
  !===============================================================================
    !> 2D Tam & Dong's BC: boundary condition at imin-jmin (edge 1,1,1 /left-bottom) - Cartesian version -
  !===============================================================================
    use mod_RFM
    use mod_eigenmode
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: cp,av,c2_
    real(wp), dimension(1:ngh+3,1:ngh+3,nz) :: rf,uf,vf,wf,pf
    real(wp), dimension(1:ngh,1:ngh,nz) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(1:ngh,1:ngh,nz) :: dux_,dvx_,dwx_,dpx_,drx_
    real(wp), dimension(1:ngh,1:ngh,nz) :: duy_,dvy_,dwy_,dpy_,dry_
    !-------------------------------------------------------------------------
    ! eigenmode disturbances
    real(wp) :: pt_in_,ut_in_,vt_in_,wt_in_,rt_in_
    real(wp) :: dp_in,du_in,dv_in,dw_in,dr_in
    real(wp), dimension(ngh,ny,nz) :: ut_in,vt_in,wt_in
    !-------------------------------------------------------------------------

    ! Pointers for polar coordinates
    ! ==============================
    ir=>BC_edge(1,1,1)%ir
    cosphi=>BC_edge(1,1,1)%cosphi
    sinphi=>BC_edge(1,1,1)%sinphi

    ! Compute fluctuations
    ! ====================
    do i=1,nghp3
       do j=1,nghp3
          do k=1,nz
             rf(i,j,k)=rho_n(i,j,k)-BC_face(1,1)%U0(i,j,k,1)
             uf(i,j,k)=   uu(i,j,k)-BC_face(1,1)%U0(i,j,k,2)
             vf(i,j,k)=   vv(i,j,k)-BC_face(1,1)%U0(i,j,k,3)
             wf(i,j,k)=   ww(i,j,k)-BC_face(1,1)%U0(i,j,k,4)
             pf(i,j,k)=  prs(i,j,k)-BC_face(1,1)%U0(i,j,k,5)
          enddo
       enddo
    enddo

    ! Compute group velocity vg
    ! =========================
    ! vg=u0*cos(phi)+v0*sin(phi) + sqrt(c0^2-w0^2-(u0*sin(phi)+v0*cos(phi))^2)
    do i=1,ngh
       do j=1,ngh
          do k=1,nz
             vg(i,j,k)= BC_face(1,1)%U0(i,j,k,2)*cosphi(i,j)+BC_face(1,1)%U0(i,j,k,3)*sinphi(i,j) &
                  + sqrt(BC_face(1,1)%U0(i,j,k,6)-BC_face(1,1)%U0(i,j,k,4)**2 &
                  -(BC_face(1,1)%U0(i,j,k,2)*sinphi(i,j)-BC_face(1,1)%U0(i,j,k,3)*cosphi(i,j))**2)
          enddo
       enddo
    enddo

    ! Non-centered derivatives in x-direction *cos(phi)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    i=1
    do j=1,ngh
       do k=1,nz
          dux_(i,j,k)= ( as4p0(1)*uf(1,j,k)+as4p0(2)*uf(2,j,k) &
                       + as4p0(3)*uf(3,j,k)+as4p0(4)*uf(4,j,k) ) *idx(i)*cosphi(i,j)
          dvx_(i,j,k)= ( as4p0(1)*vf(1,j,k)+as4p0(2)*vf(2,j,k) &
                       + as4p0(3)*vf(3,j,k)+as4p0(4)*vf(4,j,k) ) *idx(i)*cosphi(i,j)
          dwx_(i,j,k)= ( as4p0(1)*wf(1,j,k)+as4p0(2)*wf(2,j,k) &
                       + as4p0(3)*wf(3,j,k)+as4p0(4)*wf(4,j,k) ) *idx(i)*cosphi(i,j)
          dpx_(i,j,k)= ( as4p0(1)*pf(1,j,k)+as4p0(2)*pf(2,j,k) &
                       + as4p0(3)*pf(3,j,k)+as4p0(4)*pf(4,j,k) ) *idx(i)*cosphi(i,j)
          drx_(i,j,k)= ( as4p0(1)*rf(1,j,k)+as4p0(2)*rf(2,j,k) &
                       + as4p0(3)*rf(3,j,k)+as4p0(4)*rf(4,j,k) ) *idx(i)*cosphi(i,j)
       enddo
    enddo

    i=2
    do j=1,ngh
       do k=1,nz
          dux_(i,j,k)= ( as4p1(1)*uf(1,j,k)+as4p1(2)*uf(2,j,k) &
                       + as4p1(3)*uf(3,j,k)+as4p1(4)*uf(4,j,k) ) *idx(i)*cosphi(i,j)
          dvx_(i,j,k)= ( as4p1(1)*vf(1,j,k)+as4p1(2)*vf(2,j,k) &
                       + as4p1(3)*vf(3,j,k)+as4p1(4)*vf(4,j,k) ) *idx(i)*cosphi(i,j)
          dwx_(i,j,k)= ( as4p1(1)*wf(1,j,k)+as4p1(2)*wf(2,j,k) &
                       + as4p1(3)*wf(3,j,k)+as4p1(4)*wf(4,j,k) ) *idx(i)*cosphi(i,j)
          dpx_(i,j,k)= ( as4p1(1)*pf(1,j,k)+as4p1(2)*pf(2,j,k) &
                       + as4p1(3)*pf(3,j,k)+as4p1(4)*pf(4,j,k) ) *idx(i)*cosphi(i,j)
          drx_(i,j,k)= ( as4p1(1)*rf(1,j,k)+as4p1(2)*rf(2,j,k) &
                       + as4p1(3)*rf(3,j,k)+as4p1(4)*rf(4,j,k) ) *idx(i)*cosphi(i,j)
       enddo
    enddo

    i=3
    do j=1,ngh
       do k=1,nz
          dux_(i,j,k)= ( as4p2(1)*uf(1,j,k)+as4p2(2)*uf(2,j,k) &
                       + as4p2(3)*uf(3,j,k)+as4p2(4)*uf(4,j,k) &
                       + as4p2(5)*uf(5,j,k) ) *idx(i)*cosphi(i,j)
          dvx_(i,j,k)= ( as4p2(1)*vf(1,j,k)+as4p2(2)*vf(2,j,k) &
                       + as4p2(3)*vf(3,j,k)+as4p2(4)*vf(4,j,k) &
                       + as4p2(5)*vf(5,j,k) ) *idx(i)*cosphi(i,j)
          dwx_(i,j,k)= ( as4p2(1)*wf(1,j,k)+as4p2(2)*wf(2,j,k) &
                       + as4p2(3)*wf(3,j,k)+as4p2(4)*wf(4,j,k) &
                       + as4p2(5)*wf(5,j,k) ) *idx(i)*cosphi(i,j)
          dpx_(i,j,k)= ( as4p2(1)*pf(1,j,k)+as4p2(2)*pf(2,j,k) &
                       + as4p2(3)*pf(3,j,k)+as4p2(4)*pf(4,j,k) &
                       + as4p2(5)*pf(5,j,k) ) *idx(i)*cosphi(i,j)
          drx_(i,j,k)= ( as4p2(1)*rf(1,j,k)+as4p2(2)*rf(2,j,k) &
                       + as4p2(3)*rf(3,j,k)+as4p2(4)*rf(4,j,k) &
                       + as4p2(5)*rf(5,j,k) ) *idx(i)*cosphi(i,j)
       enddo
    enddo

    i=4
    do j=1,ngh
       do k=1,nz
          dux_(i,j,k)= ( as4p3(1)*uf(1,j,k)+as4p3(2)*uf(2,j,k) &
                       + as4p3(3)*uf(3,j,k)+as4p3(4)*uf(4,j,k) &
                       + as4p3(5)*uf(5,j,k)+as4p3(6)*uf(6,j,k) ) *idx(i)*cosphi(i,j)
          dvx_(i,j,k)= ( as4p3(1)*vf(1,j,k)+as4p3(2)*vf(2,j,k) &
                       + as4p3(3)*vf(3,j,k)+as4p3(4)*vf(4,j,k) &
                       + as4p3(5)*vf(5,j,k)+as4p3(6)*vf(6,j,k) ) *idx(i)*cosphi(i,j)
          dwx_(i,j,k)= ( as4p3(1)*wf(1,j,k)+as4p3(2)*wf(2,j,k) &
                       + as4p3(3)*wf(3,j,k)+as4p3(4)*wf(4,j,k) &
                       + as4p3(5)*wf(5,j,k)+as4p3(6)*wf(6,j,k) ) *idx(i)*cosphi(i,j)
          dpx_(i,j,k)= ( as4p3(1)*pf(1,j,k)+as4p3(2)*pf(2,j,k) &
                       + as4p3(3)*pf(3,j,k)+as4p3(4)*pf(4,j,k) &
                       + as4p3(5)*pf(5,j,k)+as4p3(6)*pf(6,j,k) ) *idx(i)*cosphi(i,j)
          drx_(i,j,k)= ( as4p3(1)*rf(1,j,k)+as4p3(2)*rf(2,j,k) &
                       + as4p3(3)*rf(3,j,k)+as4p3(4)*rf(4,j,k) &
                       + as4p3(5)*rf(5,j,k)+as4p3(6)*rf(6,j,k) ) *idx(i)*cosphi(i,j)
       enddo
    enddo

    do i=5,ngh
       do j=1,ngh
          do k=1,nz
             dux_(i,j,k)= ( a5(1)* (uf(i+1,j,k)-uf(i-1,j,k)) &
                          + a5(2)* (uf(i+2,j,k)-uf(i-2,j,k)) ) *idx(i)*cosphi(i,j)
             dvx_(i,j,k)= ( a5(1)* (vf(i+1,j,k)-vf(i-1,j,k)) &
                          + a5(2)* (vf(i+2,j,k)-vf(i-2,j,k)) ) *idx(i)*cosphi(i,j)
             dwx_(i,j,k)= ( a5(1)* (wf(i+1,j,k)-wf(i-1,j,k)) &
                          + a5(2)* (wf(i+2,j,k)-wf(i-2,j,k)) ) *idx(i)*cosphi(i,j)
             dpx_(i,j,k)= ( a5(1)* (pf(i+1,j,k)-pf(i-1,j,k)) &
                          + a5(2)* (pf(i+2,j,k)-pf(i-2,j,k)) ) *idx(i)*cosphi(i,j)
             drx_(i,j,k)= ( a5(1)* (rf(i+1,j,k)-rf(i-1,j,k)) &
                          + a5(2)* (rf(i+2,j,k)-rf(i-2,j,k)) ) *idx(i)*cosphi(i,j)
          enddo
       enddo
    enddo

    ! Non-centered derivatives in y-direction *sin(phi)
    ! =======================================
    j=1
    do i=1,ngh
       do k=1,nz
          duy_(i,j,k)= ( as4p0(1)*uf(i,1,k)+as4p0(2)*uf(i,2,k) &
                       + as4p0(3)*uf(i,3,k)+as4p0(4)*uf(i,4,k) ) *idy(j)*sinphi(i,j)
          dvy_(i,j,k)= ( as4p0(1)*vf(i,1,k)+as4p0(2)*vf(i,2,k) &
                       + as4p0(3)*vf(i,3,k)+as4p0(4)*vf(i,4,k) ) *idy(j)*sinphi(i,j)
          dwy_(i,j,k)= ( as4p0(1)*wf(i,1,k)+as4p0(2)*wf(i,2,k) &
                       + as4p0(3)*wf(i,3,k)+as4p0(4)*wf(i,4,k) ) *idy(j)*sinphi(i,j)
          dpy_(i,j,k)= ( as4p0(1)*pf(i,1,k)+as4p0(2)*pf(i,2,k) &
                       + as4p0(3)*pf(i,3,k)+as4p0(4)*pf(i,4,k) ) *idy(j)*sinphi(i,j)
          dry_(i,j,k)= ( as4p0(1)*rf(i,1,k)+as4p0(2)*rf(i,2,k) &
                       + as4p0(3)*rf(i,3,k)+as4p0(4)*rf(i,4,k) ) *idy(j)*sinphi(i,j)
       enddo
    enddo

    j=2
    do i=1,ngh
       do k=1,nz
          duy_(i,j,k)= ( as4p1(1)*uf(i,1,k)+as4p1(2)*uf(i,2,k) &
                       + as4p1(3)*uf(i,3,k)+as4p1(4)*uf(i,4,k) ) *idy(j)*sinphi(i,j)
          dvy_(i,j,k)= ( as4p1(1)*vf(i,1,k)+as4p1(2)*vf(i,2,k) &
                       + as4p1(3)*vf(i,3,k)+as4p1(4)*vf(i,4,k) ) *idy(j)*sinphi(i,j)
          dwy_(i,j,k)= ( as4p1(1)*wf(i,1,k)+as4p1(2)*wf(i,2,k) &
                       + as4p1(3)*wf(i,3,k)+as4p1(4)*wf(i,4,k) ) *idy(j)*sinphi(i,j)
          dpy_(i,j,k)= ( as4p1(1)*pf(i,1,k)+as4p1(2)*pf(i,2,k) &
                       + as4p1(3)*pf(i,3,k)+as4p1(4)*pf(i,4,k) ) *idy(j)*sinphi(i,j)
          dry_(i,j,k)= ( as4p1(1)*rf(i,1,k)+as4p1(2)*rf(i,2,k) &
                       + as4p1(3)*rf(i,3,k)+as4p1(4)*rf(i,4,k) ) *idy(j)*sinphi(i,j)
       enddo
    enddo

    j=3
    do i=1,ngh
       do k=1,nz
          duy_(i,j,k)= ( as4p2(1)*uf(i,1,k)+as4p2(2)*uf(i,2,k) &
                       + as4p2(3)*uf(i,3,k)+as4p2(4)*uf(i,4,k) &
                       + as4p2(5)*uf(i,5,k) ) *idy(j)*sinphi(i,j)
          dvy_(i,j,k)= ( as4p2(1)*vf(i,1,k)+as4p2(2)*vf(i,2,k) &
                       + as4p2(3)*vf(i,3,k)+as4p2(4)*vf(i,4,k) &
                       + as4p2(5)*vf(i,5,k) ) *idy(j)*sinphi(i,j)
          dwy_(i,j,k)= ( as4p2(1)*wf(i,1,k)+as4p2(2)*wf(i,2,k) &
                       + as4p2(3)*wf(i,3,k)+as4p2(4)*wf(i,4,k) &
                       + as4p2(5)*wf(i,5,k) ) *idy(j)*sinphi(i,j)
          dpy_(i,j,k)= ( as4p2(1)*pf(i,1,k)+as4p2(2)*pf(i,2,k) &
                       + as4p2(3)*pf(i,3,k)+as4p2(4)*pf(i,4,k) &
                       + as4p2(5)*pf(i,5,k) ) *idy(j)*sinphi(i,j)
          dry_(i,j,k)= ( as4p2(1)*rf(i,1,k)+as4p2(2)*rf(i,2,k) &
                       + as4p2(3)*rf(i,3,k)+as4p2(4)*rf(i,4,k) &
                       + as4p2(5)*rf(i,5,k) ) *idy(j)*sinphi(i,j)
       enddo
    enddo

    j=4
    do i=1,ngh
       do k=1,nz
          duy_(i,j,k)= ( as4p3(1)*uf(i,1,k)+as4p3(2)*uf(i,2,k) &
                       + as4p3(3)*uf(i,3,k)+as4p3(4)*uf(i,4,k) &
                       + as4p3(5)*uf(i,5,k)+as4p3(6)*uf(i,6,k) ) *idy(j)*sinphi(i,j)
          dvy_(i,j,k)= ( as4p3(1)*vf(i,1,k)+as4p3(2)*vf(i,2,k) &
                       + as4p3(3)*vf(i,3,k)+as4p3(4)*vf(i,4,k) &
                       + as4p3(5)*vf(i,5,k)+as4p3(6)*vf(i,6,k) ) *idy(j)*sinphi(i,j)
          dwy_(i,j,k)= ( as4p3(1)*wf(i,1,k)+as4p3(2)*wf(i,2,k) &
                       + as4p3(3)*wf(i,3,k)+as4p3(4)*wf(i,4,k) &
                       + as4p3(5)*wf(i,5,k)+as4p3(6)*wf(i,6,k) ) *idy(j)*sinphi(i,j)
          dpy_(i,j,k)= ( as4p3(1)*pf(i,1,k)+as4p3(2)*pf(i,2,k) &
                       + as4p3(3)*pf(i,3,k)+as4p3(4)*pf(i,4,k) &
                       + as4p3(5)*pf(i,5,k)+as4p3(6)*pf(i,6,k) ) *idy(j)*sinphi(i,j)
          dry_(i,j,k)= ( as4p3(1)*rf(i,1,k)+as4p3(2)*rf(i,2,k) &
                       + as4p3(3)*rf(i,3,k)+as4p3(4)*rf(i,4,k) &
                       + as4p3(5)*rf(i,5,k)+as4p3(6)*rf(i,6,k) ) *idy(j)*sinphi(i,j)
       enddo
    enddo

    do j=5,ngh
       do i=1,ngh
          do k=1,nz
             duy_(i,j,k)= ( a5(1)*(uf(i,j+1,k)-uf(i,j-1,k)) &
                          + a5(2)*(uf(i,j+2,k)-uf(i,j-2,k)) ) *idy(j)*sinphi(i,j)
             dvy_(i,j,k)= ( a5(1)*(vf(i,j+1,k)-vf(i,j-1,k)) &
                          + a5(2)*(vf(i,j+2,k)-vf(i,j-2,k)) ) *idy(j)*sinphi(i,j)
             dwy_(i,j,k)= ( a5(1)*(wf(i,j+1,k)-wf(i,j-1,k)) &
                          + a5(2)*(wf(i,j+2,k)-wf(i,j-2,k)) ) *idy(j)*sinphi(i,j)
             dpy_(i,j,k)= ( a5(1)*(pf(i,j+1,k)-pf(i,j-1,k)) &
                          + a5(2)*(pf(i,j+2,k)-pf(i,j-2,k)) ) *idy(j)*sinphi(i,j)
             dry_(i,j,k)= ( a5(1)*(rf(i,j+1,k)-rf(i,j-1,k)) &
                          + a5(2)*(rf(i,j+2,k)-rf(i,j-2,k)) ) *idy(j)*sinphi(i,j)
          enddo
       enddo
    enddo

    ! Compute rt,ut,vt,wt & pt
    ! ========================
    do i=1,ngh
       do j=1,ngh
          do k=1,nz
             pt(i,j,k) = vg(i,j,k)*(dpx_(i,j,k)+dpy_(i,j,k)+pf(i,j,k)*ir(i,j))
             ut(i,j,k) = vg(i,j,k)*(dux_(i,j,k)+duy_(i,j,k)+uf(i,j,k)*ir(i,j))
             vt(i,j,k) = vg(i,j,k)*(dvx_(i,j,k)+dvy_(i,j,k)+vf(i,j,k)*ir(i,j))
             wt(i,j,k) = vg(i,j,k)*(dwx_(i,j,k)+dwy_(i,j,k)+wf(i,j,k)*ir(i,j))
             rt(i,j,k) = vg(i,j,k)*(drx_(i,j,k)+dry_(i,j,k)+rf(i,j,k)*ir(i,j))
          enddo
       enddo
    enddo

    if (is_eigenmode) then
       do i=1,ngh
          do j=1,ngh
             do k=1,nz
                call eig_disturb2_imin_jmin(i,j,k,pt_in_,ut_in_,vt_in_,wt_in_,rt_in_,dp_in,du_in,dv_in,dw_in,dr_in)
                pt(i,j,k) = pt(i,j,k) - vg(i,j,k)*dp_in - pt_in_
                ut(i,j,k) = ut(i,j,k) - vg(i,j,k)*du_in - ut_in_
                vt(i,j,k) = vt(i,j,k) - vg(i,j,k)*dv_in - vt_in_
                wt(i,j,k) = wt(i,j,k) - vg(i,j,k)*dw_in - wt_in_
                rt(i,j,k) = rt(i,j,k) - vg(i,j,k)*dr_in - rt_in_
             enddo
          enddo
       enddo
    endif

    if (is_RFM) then
       call disturb_inlet_RFM_TamDong_imin_jmin(vg,ut_in,vt_in,wt_in)
       do i=1,ngh
          do j=1,ngh
             do k=1,nz
                ! ut_in = vg(i,j,k)*du_in + ut_in
                ut(i,j,k) = ut(i,j,k) - ut_in(i,j,k)
                vt(i,j,k) = vt(i,j,k) - vt_in(i,j,k)
                wt(i,j,k) = wt(i,j,k) - wt_in(i,j,k)
             enddo
          enddo
       enddo
    endif

    ! Update fluxes at each RK step
    ! =============================
    do k=1,nz
       do j=1,ngh
          do i=1,ngh
             cp=cpcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
             av=avcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
             c2_=c2calc_tro(Tmp(i,j,k),rho_n(i,j,k))

             Krho(i,j,k)  = rt(i,j,k)
             Krhou(i,j,k) = uu(i,j,k)*rt(i,j,k)+rho_n(i,j,k)*ut(i,j,k)
             Krhov(i,j,k) = vv(i,j,k)*rt(i,j,k)+rho_n(i,j,k)*vt(i,j,k)
             Krhow(i,j,k) = ww(i,j,k)*rt(i,j,k)+rho_n(i,j,k)*wt(i,j,k)
             Krhoe(i,j,k) = cp/av*(pt(i,j,k)/c2_-rt(i,j,k)) &
                  + (rhoe_n(i,j,k)+prs(i,j,k))/rho_n(i,j,k)*rt(i,j,k) &
                  + rho_n(i,j,k)*(uu(i,j,k)*ut(i,j,k)+vv(i,j,k)*vt(i,j,k)+ww(i,j,k)*wt(i,j,k))
          enddo
       enddo
    enddo

  end subroutine bc_TD2d_imin_jmin_SBP4

  !===============================================================================
  module subroutine bc_TD2d_imin_jmax_SBP4
  !===============================================================================
    !> 2D Tam & Dong's BC: boundary condition at imin-jmax (edge 1,1,2 /left-top) - Cartesian version -
  !===============================================================================
    use mod_RFM
    use mod_eigenmode
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l
    real(wp) :: cp,av,c2_
    real(wp), dimension(1:ngh+3,-2:ngh,nz) :: rf,uf,vf,wf,pf
    real(wp), dimension(1:ngh,1:ngh,nz) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(1:ngh,1:ngh,nz) :: dux_,dvx_,dwx_,dpx_,drx_
    real(wp), dimension(1:ngh,1:ngh,nz) :: duy_,dvy_,dwy_,dpy_,dry_
    !-------------------------------------------------------------------------
    ! eigenmode or RFM disturbances
    ! real(wp) :: pt_in_,ut_in_,vt_in_,wt_in_,rt_in_
    ! real(wp) :: dp_in,du_in,dv_in,dw_in,dr_in
    !-------------------------------------------------------------------------

    ! Pointers for polar coordinates
    ! ==============================
    ir=>BC_edge(1,1,2)%ir
    cosphi=>BC_edge(1,1,2)%cosphi
    sinphi=>BC_edge(1,1,2)%sinphi

    ! Compute fluctuations
    ! ====================
    do i=1,nghp3
       do j=nymngh-2,ny
          l=j-nymngh
          do k=1,nz
             rf(i,l,k)=rho_n(i,j,k)-BC_face(1,1)%U0(i,j,k,1)
             uf(i,l,k)=   uu(i,j,k)-BC_face(1,1)%U0(i,j,k,2)
             vf(i,l,k)=   vv(i,j,k)-BC_face(1,1)%U0(i,j,k,3)
             wf(i,l,k)=   ww(i,j,k)-BC_face(1,1)%U0(i,j,k,4)
             pf(i,l,k)=  prs(i,j,k)-BC_face(1,1)%U0(i,j,k,5)
          enddo
       enddo
    enddo

    ! Compute group velocity vg
    ! =========================
    ! vg=u0*cos(phi)+v0*sin(phi) + sqrt(c0^2-w0^2-(u0*sin(phi)+v0*cos(phi))^2)
    do i=1,ngh
       do j=nymnghp1,ny
          l=j-nymngh
          do k=1,nz
             vg(i,l,k)= BC_face(1,1)%U0(i,j,k,2)*cosphi(i,l)+BC_face(1,1)%U0(i,j,k,3)*sinphi(i,l) &
                  + sqrt(BC_face(1,1)%U0(i,j,k,6)-BC_face(1,1)%U0(i,j,k,4)**2 &
                  -(BC_face(1,1)%U0(i,j,k,2)*sinphi(i,l)-BC_face(1,1)%U0(i,j,k,3)*cosphi(i,l))**2)
          enddo
       enddo
    enddo

    ! Non-centered derivatives in x-direction *cos(phi)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    i=1
    do j=nymnghp1,ny
       l=j-nymngh
       do k=1,nz
          dux_(i,l,k)= ( as4p0(1)*uf(1,l,k)+as4p0(2)*uf(2,l,k) &
                       + as4p0(3)*uf(3,l,k)+as4p0(4)*uf(4,l,k) ) *idx(i)*cosphi(i,l)
          dvx_(i,l,k)= ( as4p0(1)*vf(1,l,k)+as4p0(2)*vf(2,l,k) &
                       + as4p0(3)*vf(3,l,k)+as4p0(4)*vf(4,l,k) ) *idx(i)*cosphi(i,l)
          dwx_(i,l,k)= ( as4p0(1)*wf(1,l,k)+as4p0(2)*wf(2,l,k) &
                       + as4p0(3)*wf(3,l,k)+as4p0(4)*wf(4,l,k) ) *idx(i)*cosphi(i,l)
          dpx_(i,l,k)= ( as4p0(1)*pf(1,l,k)+as4p0(2)*pf(2,l,k) &
                       + as4p0(3)*pf(3,l,k)+as4p0(4)*pf(4,l,k) ) *idx(i)*cosphi(i,l)
          drx_(i,l,k)= ( as4p0(1)*rf(1,l,k)+as4p0(2)*rf(2,l,k) &
                       + as4p0(3)*rf(3,l,k)+as4p0(4)*rf(4,l,k) ) *idx(i)*cosphi(i,l)
       enddo
    enddo

    i=2
    do j=nymnghp1,ny
       l=j-nymngh
       do k=1,nz
          dux_(i,l,k)= ( as4p1(1)*uf(1,l,k)+as4p1(2)*uf(2,l,k) &
                       + as4p1(3)*uf(3,l,k)+as4p1(4)*uf(4,l,k) ) *idx(i)*cosphi(i,l)
          dvx_(i,l,k)= ( as4p1(1)*vf(1,l,k)+as4p1(2)*vf(2,l,k) &
                       + as4p1(3)*vf(3,l,k)+as4p1(4)*vf(4,l,k) ) *idx(i)*cosphi(i,l)
          dwx_(i,l,k)= ( as4p1(1)*wf(1,l,k)+as4p1(2)*wf(2,l,k) &
                       + as4p1(3)*wf(3,l,k)+as4p1(4)*wf(4,l,k) ) *idx(i)*cosphi(i,l)
          dpx_(i,l,k)= ( as4p1(1)*pf(1,l,k)+as4p1(2)*pf(2,l,k) &
                       + as4p1(3)*pf(3,l,k)+as4p1(4)*pf(4,l,k) ) *idx(i)*cosphi(i,l)
          drx_(i,l,k)= ( as4p1(1)*rf(1,l,k)+as4p1(2)*rf(2,l,k) &
                       + as4p1(3)*rf(3,l,k)+as4p1(4)*rf(4,l,k) ) *idx(i)*cosphi(i,l)
       enddo
    enddo

    i=3
    do j=nymnghp1,ny
       l=j-nymngh
       do k=1,nz
          dux_(i,l,k)= ( as4p2(1)*uf(1,l,k)+as4p2(2)*uf(2,l,k) &
                       + as4p2(3)*uf(3,l,k)+as4p2(4)*uf(4,l,k) &
                       + as4p2(5)*uf(5,l,k) ) *idx(i)*cosphi(i,l)
          dvx_(i,l,k)= ( as4p2(1)*vf(1,l,k)+as4p2(2)*vf(2,l,k) &
                       + as4p2(3)*vf(3,l,k)+as4p2(4)*vf(4,l,k) &
                       + as4p2(5)*vf(5,l,k) ) *idx(i)*cosphi(i,l)
          dwx_(i,l,k)= ( as4p2(1)*wf(1,l,k)+as4p2(2)*wf(2,l,k) &
                       + as4p2(3)*wf(3,l,k)+as4p2(4)*wf(4,l,k) &
                       + as4p2(5)*wf(5,l,k) ) *idx(i)*cosphi(i,l)
          dpx_(i,l,k)= ( as4p2(1)*pf(1,l,k)+as4p2(2)*pf(2,l,k) &
                       + as4p2(3)*pf(3,l,k)+as4p2(4)*pf(4,l,k) &
                       + as4p2(5)*pf(5,l,k) ) *idx(i)*cosphi(i,l)
          drx_(i,l,k)= ( as4p2(1)*rf(1,l,k)+as4p2(2)*rf(2,l,k) &
                       + as4p2(3)*rf(3,l,k)+as4p2(4)*rf(4,l,k) &
                       + as4p2(5)*rf(5,l,k) ) *idx(i)*cosphi(i,l)
       enddo
    enddo

    i=4
    do j=nymnghp1,ny
       l=j-nymngh
       do k=1,nz
          dux_(i,l,k)= ( as4p3(1)*uf(1,l,k)+as4p3(2)*uf(2,l,k) &
                       + as4p3(3)*uf(3,l,k)+as4p3(4)*uf(4,l,k) &
                       + as4p3(5)*uf(5,l,k)+as4p3(6)*uf(6,l,k) ) *idx(i)*cosphi(i,l)
          dvx_(i,l,k)= ( as4p3(1)*vf(1,l,k)+as4p3(2)*vf(2,l,k) &
                       + as4p3(3)*vf(3,l,k)+as4p3(4)*vf(4,l,k) &
                       + as4p3(5)*vf(5,l,k)+as4p3(6)*vf(6,l,k) ) *idx(i)*cosphi(i,l)
          dwx_(i,l,k)= ( as4p3(1)*wf(1,l,k)+as4p3(2)*wf(2,l,k) &
                       + as4p3(3)*wf(3,l,k)+as4p3(4)*wf(4,l,k) &
                       + as4p3(5)*wf(5,l,k)+as4p3(6)*wf(6,l,k) ) *idx(i)*cosphi(i,l)
          dpx_(i,l,k)= ( as4p3(1)*pf(1,l,k)+as4p3(2)*pf(2,l,k) &
                       + as4p3(3)*pf(3,l,k)+as4p3(4)*pf(4,l,k) &
                       + as4p3(5)*pf(5,l,k)+as4p3(6)*pf(6,l,k) ) *idx(i)*cosphi(i,l)
          drx_(i,l,k)= ( as4p3(1)*rf(1,l,k)+as4p3(2)*rf(2,l,k) &
                       + as4p3(3)*rf(3,l,k)+as4p3(4)*rf(4,l,k) &
                       + as4p3(5)*rf(5,l,k)+as4p3(6)*rf(6,l,k) ) *idx(i)*cosphi(i,l)
       enddo
    enddo

    do i=5,ngh
       do j=nymnghp1,ny
          l=j-nymngh
          do k=1,nz
             dux_(i,l,k)= ( a5(1)* (uf(i+1,l,k)-uf(i-1,l,k)) &
                          + a5(2)* (uf(i+2,l,k)-uf(i-2,l,k)) ) *idx(i)*cosphi(i,l)
             dvx_(i,l,k)= ( a5(1)* (vf(i+1,l,k)-vf(i-1,l,k)) &
                          + a5(2)* (vf(i+2,l,k)-vf(i-2,l,k)) ) *idx(i)*cosphi(i,l)
             dwx_(i,l,k)= ( a5(1)* (wf(i+1,l,k)-wf(i-1,l,k)) &
                          + a5(2)* (wf(i+2,l,k)-wf(i-2,l,k)) ) *idx(i)*cosphi(i,l)
             dpx_(i,l,k)= ( a5(1)* (pf(i+1,l,k)-pf(i-1,l,k)) &
                          + a5(2)* (pf(i+2,l,k)-pf(i-2,l,k)) ) *idx(i)*cosphi(i,l)
             drx_(i,l,k)= ( a5(1)* (rf(i+1,l,k)-rf(i-1,l,k)) &
                          + a5(2)* (rf(i+2,l,k)-rf(i-2,l,k)) ) *idx(i)*cosphi(i,l)
          enddo
       enddo
    enddo

    ! Non-centered derivatives in y-direction *sin(phi)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    do j=nymnghp1,ny-4
       l=j-nymngh
       do i=1,ngh
          do k=1,nz
             duy_(i,l,k) = ( a5(1)*(uf(i,l+1,k)-uf(i,l-1,k)) &
                           + a5(2)*(uf(i,l+2,k)-uf(i,l-2,k)) ) *idy(j)*sinphi(i,l)
             dvy_(i,l,k) = ( a5(1)*(vf(i,l+1,k)-vf(i,l-1,k)) &
                           + a5(2)*(vf(i,l+2,k)-vf(i,l-2,k)) ) *idy(j)*sinphi(i,l)
             dwy_(i,l,k) = ( a5(1)*(wf(i,l+1,k)-wf(i,l-1,k)) &
                           + a5(2)*(wf(i,l+2,k)-wf(i,l-2,k)) ) *idy(j)*sinphi(i,l)
             dpy_(i,l,k) = ( a5(1)*(pf(i,l+1,k)-pf(i,l-1,k)) &
                           + a5(2)*(pf(i,l+2,k)-pf(i,l-2,k)) ) *idy(j)*sinphi(i,l)
             dry_(i,l,k) = ( a5(1)*(rf(i,l+1,k)-rf(i,l-1,k)) &
                           + a5(2)*(rf(i,l+2,k)-rf(i,l-2,k)) ) *idy(j)*sinphi(i,l)
          enddo
       enddo
    enddo

    j=ny-3
    l=j-nymngh
    do i=1,ngh
       do k=1,nz
          duy_(i,l,k) = ( as4m3(1)*uf(i,l+3,k)+as4m3(2)*uf(i,l+2,k) &
                        + as4m3(3)*uf(i,l+1,k)+as4m3(4)*uf(i,l  ,k) &
                        + as4m3(5)*uf(i,l-1,k)+as4m3(6)*uf(i,l-2,k) ) *idy(j)*sinphi(i,l)
          dvy_(i,l,k) = ( as4m3(1)*vf(i,l+3,k)+as4m3(2)*vf(i,l+2,k) &
                        + as4m3(3)*vf(i,l+1,k)+as4m3(4)*vf(i,l  ,k) &
                        + as4m3(5)*vf(i,l-1,k)+as4m3(6)*vf(i,l-2,k) ) *idy(j)*sinphi(i,l)
          dwy_(i,l,k) = ( as4m3(1)*wf(i,l+3,k)+as4m3(2)*wf(i,l+2,k) &
                        + as4m3(3)*wf(i,l+1,k)+as4m3(4)*wf(i,l  ,k) &
                        + as4m3(5)*wf(i,l-1,k)+as4m3(6)*wf(i,l-2,k) ) *idy(j)*sinphi(i,l)
          dpy_(i,l,k) = ( as4m3(1)*pf(i,l+3,k)+as4m3(2)*pf(i,l+2,k) &
                        + as4m3(3)*pf(i,l+1,k)+as4m3(4)*pf(i,l  ,k) &
                        + as4m3(5)*pf(i,l-1,k)+as4m3(6)*pf(i,l-2,k) ) *idy(j)*sinphi(i,l)
          dry_(i,l,k) = ( as4m3(1)*rf(i,l+3,k)+as4m3(2)*rf(i,l+2,k) &
                        + as4m3(3)*rf(i,l+1,k)+as4m3(4)*rf(i,l  ,k) &
                        + as4m3(5)*rf(i,l-1,k)+as4m3(6)*rf(i,l-2,k) ) *idy(j)*sinphi(i,l)
       enddo
    enddo

    j=ny-2
    l=j-nymngh
    do i=1,ngh
       do k=1,nz
          duy_(i,l,k) = ( as4m2(1)*uf(i,l+2,k)+as4m2(2)*uf(i,l+1,k) &
                        + as4m2(3)*uf(i,l  ,k)+as4m2(4)*uf(i,l-1,k) &
                        + as4m2(5)*uf(i,l-2,k) ) *idy(j)*sinphi(i,l)
          dvy_(i,l,k) = ( as4m2(1)*vf(i,l+2,k)+as4m2(2)*vf(i,l+1,k) &
                        + as4m2(3)*vf(i,l  ,k)+as4m2(4)*vf(i,l-1,k) &
                        + as4m2(5)*vf(i,l-2,k) ) *idy(j)*sinphi(i,l)
          dwy_(i,l,k) = ( as4m2(1)*wf(i,l+2,k)+as4m2(2)*wf(i,l+1,k) &
                        + as4m2(3)*wf(i,l  ,k)+as4m2(4)*wf(i,l-1,k) &
                        + as4m2(5)*wf(i,l-2,k) ) *idy(j)*sinphi(i,l)
          dpy_(i,l,k) = ( as4m2(1)*pf(i,l+2,k)+as4m2(2)*pf(i,l+1,k) &
                        + as4m2(3)*pf(i,l  ,k)+as4m2(4)*pf(i,l-1,k) &
                        + as4m2(5)*pf(i,l-2,k) ) *idy(j)*sinphi(i,l)
          dry_(i,l,k) = ( as4m2(1)*rf(i,l+2,k)+as4m2(2)*rf(i,l+1,k) &
                        + as4m2(3)*rf(i,l  ,k)+as4m2(4)*rf(i,l-1,k) &
                        + as4m2(5)*rf(i,l-2,k) ) *idy(j)*sinphi(i,l)
       enddo
    enddo

    j=ny-1
    l=j-nymngh
    do i=1,ngh
       do k=1,nz
          duy_(i,l,k) = ( as4m1(1)*uf(i,l+1,k)+as4m1(2)*uf(i,l  ,k) &
                        + as4m1(3)*uf(i,l-1,k)+as4m1(4)*uf(i,l-2,k) ) *idy(j)*sinphi(i,l)
          dvy_(i,l,k) = ( as4m1(1)*vf(i,l+1,k)+as4m1(2)*vf(i,l  ,k) &
                        + as4m1(3)*vf(i,l-1,k)+as4m1(4)*vf(i,l-2,k) ) *idy(j)*sinphi(i,l)
          dwy_(i,l,k) = ( as4m1(1)*wf(i,l+1,k)+as4m1(2)*wf(i,l  ,k) &
                        + as4m1(3)*wf(i,l-1,k)+as4m1(4)*wf(i,l-2,k) ) *idy(j)*sinphi(i,l)
          dpy_(i,l,k) = ( as4m1(1)*pf(i,l+1,k)+as4m1(2)*pf(i,l  ,k) &
                        + as4m1(3)*pf(i,l-1,k)+as4m1(4)*pf(i,l-2,k) ) *idy(j)*sinphi(i,l)
          dry_(i,l,k) = ( as4m1(1)*rf(i,l+1,k)+as4m1(2)*rf(i,l  ,k) &
                        + as4m1(3)*rf(i,l-1,k)+as4m1(4)*rf(i,l-2,k) ) *idy(j)*sinphi(i,l)
       enddo
    enddo

    j=ny
    l=j-nymngh
    do i=1,ngh
       do k=1,nz
          duy_(i,l,k) = ( as4m0(1)*uf(i,l  ,k)+as4m0(2)*uf(i,l-1,k) &
                        + as4m0(3)*uf(i,l-2,k)+as4m0(4)*uf(i,l-3,k) ) *idy(j)*sinphi(i,l)
          dvy_(i,l,k) = ( as4m0(1)*vf(i,l  ,k)+as4m0(2)*vf(i,l-1,k) &
                        + as4m0(3)*vf(i,l-2,k)+as4m0(4)*vf(i,l-3,k) ) *idy(j)*sinphi(i,l)
          dwy_(i,l,k) = ( as4m0(1)*wf(i,l  ,k)+as4m0(2)*wf(i,l-1,k) &
                        + as4m0(3)*wf(i,l-2,k)+as4m0(4)*wf(i,l-3,k) ) *idy(j)*sinphi(i,l)
          dpy_(i,l,k) = ( as4m0(1)*pf(i,l  ,k)+as4m0(2)*pf(i,l-1,k) &
                        + as4m0(3)*pf(i,l-2,k)+as4m0(4)*pf(i,l-3,k) ) *idy(j)*sinphi(i,l)
          dry_(i,l,k) = ( as4m0(1)*rf(i,l  ,k)+as4m0(2)*rf(i,l-1,k) &
                        + as4m0(3)*rf(i,l-2,k)+as4m0(4)*rf(i,l-3,k) ) *idy(j)*sinphi(i,l)
       enddo
    enddo

    ! Compute rt,ut,vt,wt & pt
    ! ========================
    do i=1,ngh
       do l=1,ngh
          do k=1,nz
             pt(i,l,k)= vg(i,l,k)*(dpx_(i,l,k)+ dpy_(i,l,k)+pf(i,l,k)*ir(i,l))
             ut(i,l,k)= vg(i,l,k)*(dux_(i,l,k)+ duy_(i,l,k)+uf(i,l,k)*ir(i,l))
             vt(i,l,k)= vg(i,l,k)*(dvx_(i,l,k)+ dvy_(i,l,k)+vf(i,l,k)*ir(i,l))
             wt(i,l,k)= vg(i,l,k)*(dwx_(i,l,k)+ dwy_(i,l,k)+wf(i,l,k)*ir(i,l))
             rt(i,l,k)= vg(i,l,k)*(drx_(i,l,k)+ dry_(i,l,k)+rf(i,l,k)*ir(i,l))
          enddo
       enddo
    enddo

!    if (is_eigenmode) then
!       do i=1,ngh
!          do j=nymnghp1,ny
!             l=j-nfy
!             l=j-nymngh !To be checked
!             do k=1,nz
!                call eig_disturb2_imin_jmin(i,j,k,pt_in_,ut_in_,vt_in_,wt_in_,rt_in_,dp_in,du_in,dv_in,dw_in,dr_in)
!                pt(i,l,k) = pt(i,l,k) - vg(i,j,k)*dp_in - pt_in_
!                ut(i,l,k) = ut(i,l,k) - vg(i,j,k)*du_in - ut_in_
!                vt(i,l,k) = vt(i,l,k) - vg(i,j,k)*dv_in - vt_in_
!                wt(i,l,k) = wt(i,l,k) - vg(i,j,k)*dw_in - wt_in_
!                rt(i,l,k) = rt(i,l,k) - vg(i,j,k)*dr_in - rt_in_
! !!$                pt(i,j,k) = - pt_in_
! !!$                ut(i,j,k) = - ut_in_
! !!$                vt(i,j,k) = - vt_in_
! !!$                wt(i,j,k) = - wt_in_
! !!$                rt(i,j,k) = - rt_in_
!             enddo
!          enddo
!       enddo
!    endif

!    if (is_RFM) then
!       do i=1,ngh
!          do j=nymnghp1,ny
!             l=j-nfy !To be checked
!             l=j-nymngh
!             do k=1,nz
!                call disturb_inlet_RFM_TamDong(i,j,k,ut_in,vt_in,wt_in,du_in,dv_in,dw_in)
!                ut(i,l,k) = ut(i,l,k) - vg(i,l,k)*du_in - ut_in
!                vt(i,l,k) = vt(i,l,k) - vg(i,l,k)*dv_in - vt_in
!                wt(i,l,k) = wt(i,l,k) - vg(i,l,k)*dw_in - wt_in
!                ! ut(i,j,k) = - ut_in
!                ! vt(i,j,k) = - vt_in
!                ! wt(i,j,k) = - wt_in
!             enddo
!          enddo
!       enddo
!    endif

    ! Update fluxes at each RK step
    ! =============================
    do k=1,nz
       do j=nymnghp1,ny
          l=j-nymngh
          do i=1,ngh
             cp=cpcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
             av=avcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
             c2_=c2calc_tro(Tmp(i,j,k),rho_n(i,j,k))

             Krho(i,j,k)  = rt(i,l,k)
             Krhou(i,j,k) = uu(i,j,k)*rt(i,l,k)+rho_n(i,j,k)*ut(i,l,k)
             Krhov(i,j,k) = vv(i,j,k)*rt(i,l,k)+rho_n(i,j,k)*vt(i,l,k)
             Krhow(i,j,k) = ww(i,j,k)*rt(i,l,k)+rho_n(i,j,k)*wt(i,l,k)
             Krhoe(i,j,k) = cp/av*(pt(i,l,k)/c2_-rt(i,l,k)) &
                  + (rhoe_n(i,j,k)+prs(i,j,k))/rho_n(i,j,k)*rt(i,l,k) &
                  + rho_n(i,j,k)*(uu(i,j,k)*ut(i,l,k)+vv(i,j,k)*vt(i,l,k)+ww(i,j,k)*wt(i,l,k))
          enddo
       enddo
    enddo

  end subroutine bc_TD2d_imin_jmax_SBP4

  !===============================================================================
  module subroutine bc_TD2d_imax_jmin_SBP4
  !===============================================================================
    !> 2D Tam & Dong's BC: boundary condition at imax-jmin (edge 1,2,1 /right-bottom) - Cartesian version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l
    real(wp) :: cp,av,c2_
    real(wp), dimension(-2:ngh,1:ngh+3,nz) :: rf,uf,vf,wf,pf
    real(wp), dimension(1:ngh,1:ngh,nz) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(1:ngh,1:ngh,nz) :: dux_,dvx_,dwx_,dpx_,drx_
    real(wp), dimension(1:ngh,1:ngh,nz) :: duy_,dvy_,dwy_,dpy_,dry_
    !-------------------------------------------------------------------------
   
    ! Pointers for polar coordinates
    ! ==============================
    ir=>BC_edge(1,2,1)%ir
    cosphi=>BC_edge(1,2,1)%cosphi
    sinphi=>BC_edge(1,2,1)%sinphi

    ! Compute fluctuations
    ! ====================
    do i=nxmngh-2,nx
       l=i-nxmngh
       do j=1,nghp3
          do k=1,nz
             rf(l,j,k)=rho_n(i,j,k)-BC_face(1,2)%U0(l,j,k,1)
             uf(l,j,k)=   uu(i,j,k)-BC_face(1,2)%U0(l,j,k,2)
             vf(l,j,k)=   vv(i,j,k)-BC_face(1,2)%U0(l,j,k,3)
             wf(l,j,k)=   ww(i,j,k)-BC_face(1,2)%U0(l,j,k,4)
             pf(l,j,k)=  prs(i,j,k)-BC_face(1,2)%U0(l,j,k,5)
          enddo
       enddo
    enddo

    ! Compute group velocity vg
    ! =========================
    ! vg=u0*cos(phi)+v0*sin(phi) + sqrt(c0^2-w0^2-(u0*sin(phi)+v0*cos(phi))^2)
    do l=1,ngh
       do j=1,ngh
          do k=1,nz
             vg(l,j,k)= BC_face(1,2)%U0(l,j,k,2)*cosphi(l,j)+BC_face(1,2)%U0(l,j,k,3)*sinphi(l,j) &
                  + sqrt(BC_face(1,2)%U0(l,j,k,6)-BC_face(1,2)%U0(l,j,k,4)**2 &
                  -(BC_face(1,2)%U0(l,j,k,2)*sinphi(l,j)-BC_face(1,2)%U0(l,j,k,3)*cosphi(l,j))**2)
          enddo
       enddo
    enddo

    ! Non-centered derivatives in x-direction *cos(phi)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    do i=nxmnghp1,nx-4
       l=i-nxmngh
       do j=1,ngh
          do k=1,nz
             dux_(l,j,k) = ( a5(1)* (uf(l+1,j,k)-uf(l-1,j,k)) &
                           + a5(2)* (uf(l+2,j,k)-uf(l-2,j,k)) ) *idx(i)*cosphi(l,j)
             dvx_(l,j,k) = ( a5(1)* (vf(l+1,j,k)-vf(l-1,j,k)) &
                           + a5(2)* (vf(l+2,j,k)-vf(l-2,j,k)) ) *idx(i)*cosphi(l,j)
             dwx_(l,j,k) = ( a5(1)* (wf(l+1,j,k)-wf(l-1,j,k)) &
                           + a5(2)* (wf(l+2,j,k)-wf(l-2,j,k)) ) *idx(i)*cosphi(l,j)
             dpx_(l,j,k) = ( a5(1)* (pf(l+1,j,k)-pf(l-1,j,k)) &
                           + a5(2)* (pf(l+2,j,k)-pf(l-2,j,k)) ) *idx(i)*cosphi(l,j)
             drx_(l,j,k) = ( a5(1)* (rf(l+1,j,k)-rf(l-1,j,k)) &
                           + a5(2)* (rf(l+2,j,k)-rf(l-2,j,k)) ) *idx(i)*cosphi(l,j)
          enddo
       enddo
    enddo
    
    i=nx-3
    l=i-nxmngh
    do j=1,ngh
       do k=1,nz
          dux_(l,j,k) = ( as4m3(1)*uf(l+3,j,k)+as4m3(2)*uf(l+2,j,k) &
                        + as4m3(3)*uf(l+1,j,k)+as4m3(4)*uf(l  ,j,k) &
                        + as4m3(5)*uf(l-1,j,k)+as4m3(6)*uf(l-2,j,k) ) *idx(i)*cosphi(l,j)
          dvx_(l,j,k) = ( as4m3(1)*vf(l+3,j,k)+as4m3(2)*vf(l+2,j,k) &
                        + as4m3(3)*vf(l+1,j,k)+as4m3(4)*vf(l  ,j,k) &
                        + as4m3(5)*vf(l-1,j,k)+as4m3(6)*vf(l-2,j,k) ) *idx(i)*cosphi(l,j)
          dwx_(l,j,k) = ( as4m3(1)*wf(l+3,j,k)+as4m3(2)*wf(l+2,j,k) &
                        + as4m3(3)*wf(l+1,j,k)+as4m3(4)*wf(l  ,j,k) &
                        + as4m3(5)*wf(l-1,j,k)+as4m3(6)*wf(l-2,j,k) ) *idx(i)*cosphi(l,j)
          dpx_(l,j,k) = ( as4m3(1)*pf(l+3,j,k)+as4m3(2)*pf(l+2,j,k) &
                        + as4m3(3)*pf(l+1,j,k)+as4m3(4)*pf(l  ,j,k) &
                        + as4m3(5)*pf(l-1,j,k)+as4m3(6)*pf(l-2,j,k) ) *idx(i)*cosphi(l,j)
          drx_(l,j,k) = ( as4m3(1)*rf(l+3,j,k)+as4m3(2)*rf(l+2,j,k) &
                        + as4m3(3)*rf(l+1,j,k)+as4m3(4)*rf(l  ,j,k) &
                        + as4m3(5)*rf(l-1,j,k)+as4m3(6)*rf(l-2,j,k) ) *idx(i)*cosphi(l,j)
       enddo
    enddo
    
    i=nx-2
    l=i-nxmngh
    do j=1,ngh
       do k=1,nz
          dux_(l,j,k) = ( as4m2(1)*uf(l+2,j,k)+as4m2(2)*uf(l+1,j,k) &
                        + as4m2(3)*uf(l  ,j,k)+as4m2(4)*uf(l-1,j,k) &
                        + as4m2(5)*uf(l-2,j,k) ) *idx(i)*cosphi(l,j)
          dvx_(l,j,k) = ( as4m2(1)*vf(l+2,j,k)+as4m2(2)*vf(l+1,j,k) &
                        + as4m2(3)*vf(l  ,j,k)+as4m2(4)*vf(l-1,j,k) &
                        + as4m2(5)*vf(l-2,j,k) ) *idx(i)*cosphi(l,j)
          dwx_(l,j,k) = ( as4m2(1)*wf(l+2,j,k)+as4m2(2)*wf(l+1,j,k) &
                        + as4m2(3)*wf(l  ,j,k)+as4m2(4)*wf(l-1,j,k) &
                        + as4m2(5)*wf(l-2,j,k) ) *idx(i)*cosphi(l,j)
          dpx_(l,j,k) = ( as4m2(1)*pf(l+2,j,k)+as4m2(2)*pf(l+1,j,k) &
                        + as4m2(3)*pf(l  ,j,k)+as4m2(4)*pf(l-1,j,k) &
                        + as4m2(5)*pf(l-2,j,k) ) *idx(i)*cosphi(l,j)
          drx_(l,j,k) = ( as4m2(1)*rf(l+2,j,k)+as4m2(2)*rf(l+1,j,k) &
                        + as4m2(3)*rf(l  ,j,k)+as4m2(4)*rf(l-1,j,k) &
                        + as4m2(5)*rf(l-2,j,k) ) *idx(i)*cosphi(l,j)
       enddo
    enddo
    
    i=nx-1
    l=i-nxmngh
    do j=1,ngh
       do k=1,nz
          dux_(l,j,k) = ( as4m1(1)*uf(l+1,j,k)+as4m1(2)*uf(l  ,j,k) &
                        + as4m1(3)*uf(l-1,j,k)+as4m1(4)*uf(l-2,j,k)  ) *idx(i)*cosphi(l,j)
          dvx_(l,j,k) = ( as4m1(1)*vf(l+1,j,k)+as4m1(2)*vf(l  ,j,k) &
                        + as4m1(3)*vf(l-1,j,k)+as4m1(4)*vf(l-2,j,k) ) *idx(i)*cosphi(l,j)
          dwx_(l,j,k) = ( as4m1(1)*wf(l+1,j,k)+as4m1(2)*wf(l  ,j,k) &
                        + as4m1(3)*wf(l-1,j,k)+as4m1(4)*wf(l-2,j,k) ) *idx(i)*cosphi(l,j)
          dpx_(l,j,k) = ( as4m1(1)*pf(l+1,j,k)+as4m1(2)*pf(l  ,j,k) &
                        + as4m1(3)*pf(l-1,j,k)+as4m1(4)*pf(l-2,j,k) ) *idx(i)*cosphi(l,j)
          drx_(l,j,k) = ( as4m1(1)*rf(l+1,j,k)+as4m1(2)*rf(l  ,j,k) &
                        + as4m1(3)*rf(l-1,j,k)+as4m1(4)*rf(l-2,j,k) ) *idx(i)*cosphi(l,j)
       enddo
    enddo
    
    i=nx
    l=i-nxmngh
    do j=1,ngh
       do k=1,nz
          dux_(l,j,k) = ( as4m0(1)*uf(l  ,j,k)+as4m0(2)*uf(l-1,j,k) &
                        + as4m0(3)*uf(l-2,j,k)+as4m0(4)*uf(l-3,j,k) ) *idx(i)*cosphi(l,j)
          dvx_(l,j,k) = ( as4m0(1)*vf(l  ,j,k)+as4m0(2)*vf(l-1,j,k) &
                        + as4m0(3)*vf(l-2,j,k)+as4m0(4)*vf(l-3,j,k) ) *idx(i)*cosphi(l,j)
          dwx_(l,j,k) = ( as4m0(1)*wf(l  ,j,k)+as4m0(2)*wf(l-1,j,k) &
                        + as4m0(3)*wf(l-2,j,k)+as4m0(4)*wf(l-3,j,k) ) *idx(i)*cosphi(l,j)
          dpx_(l,j,k) = ( as4m0(1)*pf(l  ,j,k)+as4m0(2)*pf(l-1,j,k) &
                        + as4m0(3)*pf(l-2,j,k)+as4m0(4)*pf(l-3,j,k) ) *idx(i)*cosphi(l,j)
          drx_(l,j,k) = ( as4m0(1)*rf(l  ,j,k)+as4m0(2)*rf(l-1,j,k) &
                        + as4m0(3)*rf(l-2,j,k)+as4m0(4)*rf(l-3,j,k) ) *idx(i)*cosphi(l,j)
       enddo
    enddo

    ! Non-centered derivatives in y-direction *sin(phi)
    ! =======================================
    j=1
    do i=nxmnghp1,nx
       l=i-nxmngh
       do k=1,nz
          duy_(l,j,k)= ( as4p0(1)*uf(l,1,k)+as4p0(2)*uf(l,2,k) &
                       + as4p0(3)*uf(l,3,k)+as4p0(4)*uf(l,4,k) ) *idy(j)*sinphi(l,j)
          dvy_(l,j,k)= ( as4p0(1)*vf(l,1,k)+as4p0(2)*vf(l,2,k) &
                       + as4p0(3)*vf(l,3,k)+as4p0(4)*vf(l,4,k) ) *idy(j)*sinphi(l,j)
          dwy_(l,j,k)= ( as4p0(1)*wf(l,1,k)+as4p0(2)*wf(l,2,k) &
                       + as4p0(3)*wf(l,3,k)+as4p0(4)*wf(l,4,k) ) *idy(j)*sinphi(l,j)
          dpy_(l,j,k)= ( as4p0(1)*pf(l,1,k)+as4p0(2)*pf(l,2,k) &
                       + as4p0(3)*pf(l,3,k)+as4p0(4)*pf(l,4,k) ) *idy(j)*sinphi(l,j)
          dry_(l,j,k)= ( as4p0(1)*rf(l,1,k)+as4p0(2)*rf(l,2,k) &
                       + as4p0(3)*rf(l,3,k)+as4p0(4)*rf(l,4,k) ) *idy(j)*sinphi(l,j)
       enddo
    enddo

    j=2
    do i=nxmnghp1,nx
       l=i-nxmngh
       do k=1,nz
          duy_(l,j,k)= ( as4p1(1)*uf(l,1,k)+as4p1(2)*uf(l,2,k) &
                       + as4p1(3)*uf(l,3,k)+as4p1(4)*uf(l,4,k)  ) *idy(j)*sinphi(l,j)
          dvy_(l,j,k)= ( as4p1(1)*vf(l,1,k)+as4p1(2)*vf(l,2,k) &
                       + as4p1(3)*vf(l,3,k)+as4p1(4)*vf(l,4,k) ) *idy(j)*sinphi(l,j)
          dwy_(l,j,k)= ( as4p1(1)*wf(l,1,k)+as4p1(2)*wf(l,2,k) &
                       + as4p1(3)*wf(l,3,k)+as4p1(4)*wf(l,4,k) ) *idy(j)*sinphi(l,j)
          dpy_(l,j,k)= ( as4p1(1)*pf(l,1,k)+as4p1(2)*pf(l,2,k) &
                       + as4p1(3)*pf(l,3,k)+as4p1(4)*pf(l,4,k) ) *idy(j)*sinphi(l,j)
          dry_(l,j,k)= ( as4p1(1)*rf(l,1,k)+as4p1(2)*rf(l,2,k) &
                       + as4p1(3)*rf(l,3,k)+as4p1(4)*rf(l,4,k) ) *idy(j)*sinphi(l,j)
       enddo
    enddo

    j=3
    do i=nxmnghp1,nx
       l=i-nxmngh
       do k=1,nz
          duy_(l,j,k)= ( as4p2(1)*uf(l,1,k)+as4p2(2)*uf(l,2,k) &
                       + as4p2(3)*uf(l,3,k)+as4p2(4)*uf(l,4,k) &
                       + as4p2(5)*uf(l,5,k) ) *idy(j)*sinphi(l,j)
          dvy_(l,j,k)= ( as4p2(1)*vf(l,1,k)+as4p2(2)*vf(l,2,k) &
                       + as4p2(3)*vf(l,3,k)+as4p2(4)*vf(l,4,k) &
                       + as4p2(5)*vf(l,5,k) ) *idy(j)*sinphi(l,j)
          dwy_(l,j,k)= ( as4p2(1)*wf(l,1,k)+as4p2(2)*wf(l,2,k) &
                       + as4p2(3)*wf(l,3,k)+as4p2(4)*wf(l,4,k) &
                       + as4p2(5)*wf(l,5,k) ) *idy(j)*sinphi(l,j)
          dpy_(l,j,k)= ( as4p2(1)*pf(l,1,k)+as4p2(2)*pf(l,2,k) &
                       + as4p2(3)*pf(l,3,k)+as4p2(4)*pf(l,4,k) &
                       + as4p2(5)*pf(l,5,k) ) *idy(j)*sinphi(l,j)
          dry_(l,j,k)= ( as4p2(1)*rf(l,1,k)+as4p2(2)*rf(l,2,k) &
                       + as4p2(3)*rf(l,3,k)+as4p2(4)*rf(l,4,k) &
                       + as4p2(5)*rf(l,5,k) ) *idy(j)*sinphi(l,j)
       enddo
    enddo

    j=4
    do i=nxmnghp1,nx
       l=i-nxmngh
       do k=1,nz
          duy_(l,j,k)= ( as4p3(1)*uf(l,1,k)+as4p3(2)*uf(l,2,k) &
                       + as4p3(3)*uf(l,3,k)+as4p3(4)*uf(l,4,k) &
                       + as4p3(5)*uf(l,5,k)+as4p3(6)*uf(l,6,k) ) *idy(j)*sinphi(l,j)
          dvy_(l,j,k)= ( as4p3(1)*vf(l,1,k)+as4p3(2)*vf(l,2,k) &
                       + as4p3(3)*vf(l,3,k)+as4p3(4)*vf(l,4,k) &
                       + as4p3(5)*vf(l,5,k)+as4p3(6)*vf(l,6,k) ) *idy(j)*sinphi(l,j)
          dwy_(l,j,k)= ( as4p3(1)*wf(l,1,k)+as4p3(2)*wf(l,2,k) &
                       + as4p3(3)*wf(l,3,k)+as4p3(4)*wf(l,4,k) &
                       + as4p3(5)*wf(l,5,k)+as4p3(6)*wf(l,6,k) ) *idy(j)*sinphi(l,j)
          dpy_(l,j,k)= ( as4p3(1)*pf(l,1,k)+as4p3(2)*pf(l,2,k) &
                       + as4p3(3)*pf(l,3,k)+as4p3(4)*pf(l,4,k) &
                       + as4p3(5)*pf(l,5,k)+as4p3(6)*pf(l,6,k) ) *idy(j)*sinphi(l,j)
          dry_(l,j,k)= ( as4p3(1)*rf(l,1,k)+as4p3(2)*rf(l,2,k) &
                       + as4p3(3)*rf(l,3,k)+as4p3(4)*rf(l,4,k) &
                       + as4p3(5)*rf(l,5,k)+as4p3(6)*rf(l,6,k) ) *idy(j)*sinphi(l,j)
       enddo
    enddo

    do j=5,ngh
       do i=nxmnghp1,nx
          l=i-nxmngh
          do k=1,nz
             duy_(l,j,k)= ( a5(1)*(uf(l,j+1,k)-uf(l,j-1,k)) &
                          + a5(2)*(uf(l,j+2,k)-uf(l,j-2,k)) ) *idy(j)*sinphi(l,j)
             dvy_(l,j,k)= ( a5(1)*(vf(l,j+1,k)-vf(l,j-1,k)) &
                          + a5(2)*(vf(l,j+2,k)-vf(l,j-2,k)) ) *idy(j)*sinphi(l,j)
             dwy_(l,j,k)= ( a5(1)*(wf(l,j+1,k)-wf(l,j-1,k)) &
                          + a5(2)*(wf(l,j+2,k)-wf(l,j-2,k)) ) *idy(j)*sinphi(l,j)
             dpy_(l,j,k)= ( a5(1)*(pf(l,j+1,k)-pf(l,j-1,k)) &
                          + a5(2)*(pf(l,j+2,k)-pf(l,j-2,k)) ) *idy(j)*sinphi(l,j)
             dry_(l,j,k)= ( a5(1)*(rf(l,j+1,k)-rf(l,j-1,k)) &
                          + a5(2)*(rf(l,j+2,k)-rf(l,j-2,k)) ) *idy(j)*sinphi(l,j)
          enddo
       enddo
    enddo

    ! Compute rt,ut,vt,wt & pt
    ! ========================
    do l=1,ngh
       do j=1,ngh
          do k=1,nz
             pt(l,j,k)= vg(l,j,k) * (dpx_(l,j,k)+dpy_(l,j,k)+pf(l,j,k)*ir(l,j))
             ut(l,j,k)= vg(l,j,k) * (dux_(l,j,k)+duy_(l,j,k)+uf(l,j,k)*ir(l,j))
             vt(l,j,k)= vg(l,j,k) * (dvx_(l,j,k)+dvy_(l,j,k)+vf(l,j,k)*ir(l,j))
             wt(l,j,k)= vg(l,j,k) * (dwx_(l,j,k)+dwy_(l,j,k)+wf(l,j,k)*ir(l,j))
             rt(l,j,k)= vg(l,j,k) * (drx_(l,j,k)+dry_(l,j,k)+rf(l,j,k)*ir(l,j))
          enddo
       enddo
    enddo

    ! Update fluxes at each RK step
    ! =============================
    do k=1,nz
       do i=nxmnghp1,nx
          l=i-nxmngh
          do j=1,ngh
             cp=cpcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
             av=avcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
             c2_=c2calc_tro(Tmp(i,j,k),rho_n(i,j,k))
          
             Krho(i,j,k)  = rt(l,j,k)
             Krhou(i,j,k) = uu(i,j,k)*rt(l,j,k)+rho_n(i,j,k)*ut(l,j,k)
             Krhov(i,j,k) = vv(i,j,k)*rt(l,j,k)+rho_n(i,j,k)*vt(l,j,k)
             Krhow(i,j,k) = ww(i,j,k)*rt(l,j,k)+rho_n(i,j,k)*wt(l,j,k)
             Krhoe(i,j,k) = cp/av*(pt(l,j,k)/c2_-rt(l,j,k)) &
                  + (rhoe_n(i,j,k)+prs(i,j,k))/rho_n(i,j,k)*rt(l,j,k) &
                  + rho_n(i,j,k)*(uu(i,j,k)*ut(l,j,k)+vv(i,j,k)*vt(l,j,k)+ww(i,j,k)*wt(l,j,k))
          enddo
       enddo
    enddo

  end subroutine bc_TD2d_imax_jmin_SBP4

  !===============================================================================
  module subroutine bc_TD2d_imax_jmax_SBP4
  !===============================================================================
    !> 2D Tam & Dong's BC: boundary condition at imax-jmax (edge 1,2,2 /right-top) - Cartesian version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l,m
    real(wp) :: cp,av,c2_
    real(wp), dimension(-2:ngh,-2:ngh,nz) :: rf,uf,vf,wf,pf
    real(wp), dimension(1:ngh,1:ngh,nz) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(1:ngh,1:ngh,nz) :: dux_,dvx_,dwx_,dpx_,drx_
    real(wp), dimension(1:ngh,1:ngh,nz) :: duy_,dvy_,dwy_,dpy_,dry_
    !-------------------------------------------------------------------------
    
    ! Pointers for polar coordinates
    ! ==============================
    ir=>BC_edge(1,2,2)%ir
    cosphi=>BC_edge(1,2,2)%cosphi
    sinphi=>BC_edge(1,2,2)%sinphi

    ! Compute fluctuations
    ! ====================
    do i=nxmngh-2,nx
       l=i-nxmngh
       do j=nymngh-2,ny
          m=j-nymngh
          do k=1,nz
             rf(l,m,k)=rho_n(i,j,k)-BC_face(2,2)%U0(i,m,k,1)
             uf(l,m,k)=   uu(i,j,k)-BC_face(2,2)%U0(i,m,k,2)
             vf(l,m,k)=   vv(i,j,k)-BC_face(2,2)%U0(i,m,k,3)
             wf(l,m,k)=   ww(i,j,k)-BC_face(2,2)%U0(i,m,k,4)
             pf(l,m,k)=  prs(i,j,k)-BC_face(2,2)%U0(i,m,k,5)
          enddo
       enddo
    enddo
    
    ! Compute group velocity vg
    ! =========================
    ! vg=u0*cos(phi)+v0*sin(phi) + sqrt(c0^2-w0^2-(u0*sin(phi)+v0*cos(phi))^2)
    do i=nxmnghp1,nx
       l=i-nxmngh
       do j=nymnghp1,ny
          m=j-nymngh
          do k=1,nz
             vg(l,m,k)= BC_face(2,2)%U0(i,m,k,2)*cosphi(l,m)+BC_face(2,2)%U0(i,m,k,3)*sinphi(l,m) &
                  + sqrt(BC_face(2,2)%U0(i,m,k,6)-BC_face(2,2)%U0(i,m,k,4)**2 &
                  -(BC_face(2,2)%U0(i,m,k,2)*sinphi(l,m)-BC_face(2,2)%U0(i,m,k,3)*cosphi(l,m))**2)
          enddo
       enddo
    enddo

    ! Non-centered derivatives in x-direction *cos(phi)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    do i=nxmnghp1,nx-4
       l=i-nxmngh
       do j=nymnghp1,ny
          m=j-nymngh
          do k=1,nz
             dux_(l,m,k) = ( a5(1)* (uf(l+1,m,k)-uf(l-1,m,k)) &
                           + a5(2)* (uf(l+2,m,k)-uf(l-2,m,k)) ) *idx(i)*cosphi(l,m)
             dvx_(l,m,k) = ( a5(1)* (vf(l+1,m,k)-vf(l-1,m,k)) &
                           + a5(2)* (vf(l+2,m,k)-vf(l-2,m,k)) ) *idx(i)*cosphi(l,m)
             dwx_(l,m,k) = ( a5(1)* (wf(l+1,m,k)-wf(l-1,m,k)) &
                           + a5(2)* (wf(l+2,m,k)-wf(l-2,m,k)) ) *idx(i)*cosphi(l,m)
             dpx_(l,m,k) = ( a5(1)* (pf(l+1,m,k)-pf(l-1,m,k)) &
                           + a5(2)* (pf(l+2,m,k)-pf(l-2,m,k)) ) *idx(i)*cosphi(l,m)
             drx_(l,m,k) = ( a5(1)* (rf(l+1,m,k)-rf(l-1,m,k)) &
                           + a5(2)* (rf(l+2,m,k)-rf(l-2,m,k)) ) *idx(i)*cosphi(l,m)
          enddo
       enddo
    enddo
    
    i=nx-3
    l=i-nxmngh
    do j=nymnghp1,ny
       m=j-nymngh
       do k=1,nz
          dux_(l,m,k) = ( as4m3(1)*uf(l+3,m,k)+as4m3(2)*uf(l+2,m,k) &
                        + as4m3(3)*uf(l+1,m,k)+as4m3(4)*uf(l  ,m,k) &
                        + as4m3(5)*uf(l-1,m,k)+as4m3(6)*uf(l-2,m,k) ) *idx(i)*cosphi(l,m)
          dvx_(l,m,k) = ( as4m3(1)*vf(l+3,m,k)+as4m3(2)*vf(l+2,m,k) &
                        + as4m3(3)*vf(l+1,m,k)+as4m3(4)*vf(l  ,m,k) &
                        + as4m3(5)*vf(l-1,m,k)+as4m3(6)*vf(l-2,m,k) ) *idx(i)*cosphi(l,m)
          dwx_(l,m,k) = ( as4m3(1)*wf(l+3,m,k)+as4m3(2)*wf(l+2,m,k) &
                        + as4m3(3)*wf(l+1,m,k)+as4m3(4)*wf(l  ,m,k) &
                        + as4m3(5)*wf(l-1,m,k)+as4m3(6)*wf(l-2,m,k) ) *idx(i)*cosphi(l,m)
          dpx_(l,m,k) = ( as4m3(1)*pf(l+3,m,k)+as4m3(2)*pf(l+2,m,k) &
                        + as4m3(3)*pf(l+1,m,k)+as4m3(4)*pf(l  ,m,k) &
                        + as4m3(5)*pf(l-1,m,k)+as4m3(6)*pf(l-2,m,k) ) *idx(i)*cosphi(l,m)
          drx_(l,m,k) = ( as4m3(1)*rf(l+3,m,k)+as4m3(2)*rf(l+2,m,k) &
                        + as4m3(3)*rf(l+1,m,k)+as4m3(4)*rf(l  ,m,k) &
                        + as4m3(5)*rf(l-1,m,k)+as4m3(6)*rf(l-2,m,k) ) *idx(i)*cosphi(l,m)
       enddo
    enddo
    
    i=nx-2
    l=i-nxmngh
    do j=nymnghp1,ny
       m=j-nymngh
       do k=1,nz
          dux_(l,m,k) = ( as4m2(1)*uf(l+2,m,k)+as4m2(2)*uf(l+1,m,k) &
                        + as4m2(3)*uf(l  ,m,k)+as4m2(4)*uf(l-1,m,k) &
                        + as4m2(5)*uf(l-2,m,k) ) *idx(i)*cosphi(l,m)
          dvx_(l,m,k) = ( as4m2(1)*vf(l+2,m,k)+as4m2(2)*vf(l+1,m,k) &
                        + as4m2(3)*vf(l  ,m,k)+as4m2(4)*vf(l-1,m,k) &
                        + as4m2(5)*vf(l-2,m,k) ) *idx(i)*cosphi(l,m)
          dwx_(l,m,k) = ( as4m2(1)*wf(l+2,m,k)+as4m2(2)*wf(l+1,m,k) &
                        + as4m2(3)*wf(l  ,m,k)+as4m2(4)*wf(l-1,m,k) &
                        + as4m2(5)*wf(l-2,m,k) ) *idx(i)*cosphi(l,m)
          dpx_(l,m,k) = ( as4m2(1)*pf(l+2,m,k)+as4m2(2)*pf(l+1,m,k) &
                        + as4m2(3)*pf(l  ,m,k)+as4m2(4)*pf(l-1,m,k) &
                        + as4m2(5)*pf(l-2,m,k) ) *idx(i)*cosphi(l,m)
          drx_(l,m,k) = ( as4m2(1)*rf(l+2,m,k)+as4m2(2)*rf(l+1,m,k) &
                        + as4m2(3)*rf(l  ,m,k)+as4m2(4)*rf(l-1,m,k) &
                        + as4m2(5)*rf(l-2,m,k) ) *idx(i)*cosphi(l,m)
       enddo
    enddo
    
    i=nx-1
    l=i-nxmngh
    do j=nymnghp1,ny
       m=j-nymngh
       do k=1,nz
          dux_(l,m,k) = ( as4m1(1)*uf(l+1,m,k)+as4m1(2)*uf(l  ,m,k) &
                        + as4m1(3)*uf(l-1,m,k)+as4m1(4)*uf(l-2,m,k) ) *idx(i)*cosphi(l,m)
          dvx_(l,m,k) = ( as4m1(1)*vf(l+1,m,k)+as4m1(2)*vf(l  ,m,k) &
                        + as4m1(3)*vf(l-1,m,k)+as4m1(4)*vf(l-2,m,k) ) *idx(i)*cosphi(l,m)
          dwx_(l,m,k) = ( as4m1(1)*wf(l+1,m,k)+as4m1(2)*wf(l  ,m,k) &
                        + as4m1(3)*wf(l-1,m,k)+as4m1(4)*wf(l-2,m,k) ) *idx(i)*cosphi(l,m)
          dpx_(l,m,k) = ( as4m1(1)*pf(l+1,m,k)+as4m1(2)*pf(l  ,m,k) &
                        + as4m1(3)*pf(l-1,m,k)+as4m1(4)*pf(l-2,m,k) ) *idx(i)*cosphi(l,m)
          drx_(l,m,k) = ( as4m1(1)*rf(l+1,m,k)+as4m1(2)*rf(l  ,m,k) &
                        + as4m1(3)*rf(l-1,m,k)+as4m1(4)*rf(l-2,m,k) ) *idx(i)*cosphi(l,m)
       enddo
    enddo
    
    i=nx
    l=i-nxmngh
    do j=nymnghp1,ny
       m=j-nymngh
       do k=1,nz
          dux_(l,m,k) = ( as4m0(1)*uf(l  ,m,k)+as4m0(2)*uf(l-1,m,k) &
                        + as4m0(3)*uf(l-2,m,k)+as4m0(4)*uf(l-3,m,k) ) *idx(i)*cosphi(l,m)
          dvx_(l,m,k) = ( as4m0(1)*vf(l  ,m,k)+as4m0(2)*vf(l-1,m,k) &
                        + as4m0(3)*vf(l-2,m,k)+as4m0(4)*vf(l-3,m,k) ) *idx(i)*cosphi(l,m)
          dwx_(l,m,k) = ( as4m0(1)*wf(l  ,m,k)+as4m0(2)*wf(l-1,m,k) &
                        + as4m0(3)*wf(l-2,m,k)+as4m0(4)*wf(l-3,m,k) ) *idx(i)*cosphi(l,m)
          dpx_(l,m,k) = ( as4m0(1)*pf(l  ,m,k)+as4m0(2)*pf(l-1,m,k) &
                        + as4m0(3)*pf(l-2,m,k)+as4m0(4)*pf(l-3,m,k) ) *idx(i)*cosphi(l,m)
          drx_(l,m,k) = ( as4m0(1)*rf(l  ,m,k)+as4m0(2)*rf(l-1,m,k) &
                        + as4m0(3)*rf(l-2,m,k)+as4m0(4)*rf(l-3,m,k) ) *idx(i)*cosphi(l,m)
       enddo
    enddo
 
    ! Non-centered derivatives in y-direction *sin(phi)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    do j=nymnghp1,ny-4
       m=j-nymngh
       do l=1,ngh
          do k=1,nz
             duy_(l,m,k) = ( a5(1)*(uf(l,m+1,k)-uf(l,m-1,k)) &
                           + a5(2)*(uf(l,m+2,k)-uf(l,m-2,k)) ) *idy(j)*sinphi(l,m)
             dvy_(l,m,k) = ( a5(1)*(vf(l,m+1,k)-vf(l,m-1,k)) &
                           + a5(2)*(vf(l,m+2,k)-vf(l,m-2,k)) ) *idy(j)*sinphi(l,m)
             dwy_(l,m,k) = ( a5(1)*(wf(l,m+1,k)-wf(l,m-1,k)) &
                           + a5(2)*(wf(l,m+2,k)-wf(l,m-2,k)) ) *idy(j)*sinphi(l,m)
             dpy_(l,m,k) = ( a5(1)*(pf(l,m+1,k)-pf(l,m-1,k)) &
                           + a5(2)*(pf(l,m+2,k)-pf(l,m-2,k)) ) *idy(j)*sinphi(l,m)
             dry_(l,m,k) = ( a5(1)*(rf(l,m+1,k)-rf(l,m-1,k)) &
                           + a5(2)*(rf(l,m+2,k)-rf(l,m-2,k)) ) *idy(j)*sinphi(l,m)
          enddo
       enddo
    enddo

    j=ny-3
    m=j-nymngh
    do l=1,ngh
       do k=1,nz
          duy_(l,m,k) = ( as4m3(1)*uf(l,m+3,k)+as4m3(2)*uf(l,m+2,k) &
                        + as4m3(3)*uf(l,m+1,k)+as4m3(4)*uf(l,m  ,k) &
                        + as4m3(5)*uf(l,m-1,k)+as4m3(6)*uf(l,m-2,k) ) *idy(j)*sinphi(l,m)
          dvy_(l,m,k) = ( as4m3(1)*vf(l,m+3,k)+as4m3(2)*vf(l,m+2,k) &
                        + as4m3(3)*vf(l,m+1,k)+as4m3(4)*vf(l,m  ,k) &
                        + as4m3(5)*vf(l,m-1,k)+as4m3(6)*vf(l,m-2,k) ) *idy(j)*sinphi(l,m)
          dwy_(l,m,k) = ( as4m3(1)*wf(l,m+3,k)+as4m3(2)*wf(l,m+2,k) &
                        + as4m3(3)*wf(l,m+1,k)+as4m3(4)*wf(l,m  ,k) &
                        + as4m3(5)*wf(l,m-1,k)+as4m3(6)*wf(l,m-2,k) ) *idy(j)*sinphi(l,m)
          dpy_(l,m,k) = ( as4m3(1)*pf(l,m+3,k)+as4m3(2)*pf(l,m+2,k) &
                        + as4m3(3)*pf(l,m+1,k)+as4m3(4)*pf(l,m  ,k) &
                        + as4m3(5)*pf(l,m-1,k)+as4m3(6)*pf(l,m-2,k) ) *idy(j)*sinphi(l,m)
          dry_(l,m,k) = ( as4m3(1)*rf(l,m+3,k)+as4m3(2)*rf(l,m+2,k) &
                        + as4m3(3)*rf(l,m+1,k)+as4m3(4)*rf(l,m  ,k) &
                        + as4m3(5)*rf(l,m-1,k)+as4m3(6)*rf(l,m-2,k) ) *idy(j)*sinphi(l,m)
       enddo
    enddo
    
    j=ny-2
    m=j-nymngh
    do l=1,ngh
       do k=1,nz
          duy_(l,m,k) = ( as4m2(1)*uf(l,m+2,k)+as4m2(2)*uf(l,m+1,k) &
                        + as4m2(3)*uf(l,m  ,k)+as4m2(4)*uf(l,m-1,k) &
                        + as4m2(5)*uf(l,m-2,k) ) *idy(j)*sinphi(l,m)
          dvy_(l,m,k) = ( as4m2(1)*vf(l,m+2,k)+as4m2(2)*vf(l,m+1,k) &
                        + as4m2(3)*vf(l,m  ,k)+as4m2(4)*vf(l,m-1,k) &
                        + as4m2(5)*vf(l,m-2,k) ) *idy(j)*sinphi(l,m)
          dwy_(l,m,k) = ( as4m2(1)*wf(l,m+2,k)+as4m2(2)*wf(l,m+1,k) &
                        + as4m2(3)*wf(l,m  ,k)+as4m2(4)*wf(l,m-1,k) &
                        + as4m2(5)*wf(l,m-2,k) ) *idy(j)*sinphi(l,m)
          dpy_(l,m,k) = ( as4m2(1)*pf(l,m+2,k)+as4m2(2)*pf(l,m+1,k) &
                        + as4m2(3)*pf(l,m  ,k)+as4m2(4)*pf(l,m-1,k) &
                        + as4m2(5)*pf(l,m-2,k) ) *idy(j)*sinphi(l,m)
          dry_(l,m,k) = ( as4m2(1)*rf(l,m+2,k)+as4m2(2)*rf(l,m+1,k) &
                        + as4m2(3)*rf(l,m  ,k)+as4m2(4)*rf(l,m-1,k) &
                        + as4m2(5)*rf(l,m-2,k) ) *idy(j)*sinphi(l,m)
       enddo
    enddo
    
    j=ny-1
    m=j-nymngh
    do l=1,ngh
       do k=1,nz
          duy_(l,m,k) = ( as4m1(1)*uf(l,m+1,k)+as4m1(2)*uf(l,m  ,k) &
                        + as4m1(3)*uf(l,m-1,k)+as4m1(4)*uf(l,m-2,k) ) *idy(j)*sinphi(l,m)
          dvy_(l,m,k) = ( as4m1(1)*vf(l,m+1,k)+as4m1(2)*vf(l,m  ,k) &
                        + as4m1(3)*vf(l,m-1,k)+as4m1(4)*vf(l,m-2,k) ) *idy(j)*sinphi(l,m)
          dwy_(l,m,k) = ( as4m1(1)*wf(l,m+1,k)+as4m1(2)*wf(l,m  ,k) &
                        + as4m1(3)*wf(l,m-1,k)+as4m1(4)*wf(l,m-2,k) ) *idy(j)*sinphi(l,m)
          dpy_(l,m,k) = ( as4m1(1)*pf(l,m+1,k)+as4m1(2)*pf(l,m  ,k) &
                        + as4m1(3)*pf(l,m-1,k)+as4m1(4)*pf(l,m-2,k) ) *idy(j)*sinphi(l,m)
          dry_(l,m,k) = ( as4m1(1)*rf(l,m+1,k)+as4m1(2)*rf(l,m  ,k) &
                        + as4m1(3)*rf(l,m-1,k)+as4m1(4)*rf(l,m-2,k) ) *idy(j)*sinphi(l,m)
       enddo
    enddo
    
    j=ny
    m=j-nymngh
    do l=1,ngh
       do k=1,nz
          duy_(l,m,k) = ( as4m0(1)*uf(l,m  ,k)+as4m0(2)*uf(l,m-1,k) &
                        + as4m0(3)*uf(l,m-2,k)+as4m0(4)*uf(l,m-3,k) ) *idy(j)*sinphi(l,m)
          dvy_(l,m,k) = ( as4m0(1)*vf(l,m  ,k)+as4m0(2)*vf(l,m-1,k) &
                        + as4m0(3)*vf(l,m-2,k)+as4m0(4)*vf(l,m-3,k) ) *idy(j)*sinphi(l,m)
          dwy_(l,m,k) = ( as4m0(1)*wf(l,m  ,k)+as4m0(2)*wf(l,m-1,k) &
                        + as4m0(3)*wf(l,m-2,k)+as4m0(4)*wf(l,m-3,k) ) *idy(j)*sinphi(l,m)
          dpy_(l,m,k) = ( as4m0(1)*pf(l,m  ,k)+as4m0(2)*pf(l,m-1,k) &
                        + as4m0(3)*pf(l,m-2,k)+as4m0(4)*pf(l,m-3,k) ) *idy(j)*sinphi(l,m)
          dry_(l,m,k) = ( as4m0(1)*rf(l,m  ,k)+as4m0(2)*rf(l,m-1,k) &
                        + as4m0(3)*rf(l,m-2,k)+as4m0(4)*rf(l,m-3,k) ) *idy(j)*sinphi(l,m)
       enddo
    enddo
              
    ! Compute rt,ut,vt,wt & pt
    ! ========================
    do l=1,ngh
       do m=1,ngh
          do k=1,nz
             pt(l,m,k) = vg(l,m,k) * (dpx_(l,m,k)+dpy_(l,m,k)+pf(l,m,k)*ir(l,m))
             ut(l,m,k) = vg(l,m,k) * (dux_(l,m,k)+duy_(l,m,k)+uf(l,m,k)*ir(l,m))
             vt(l,m,k) = vg(l,m,k) * (dvx_(l,m,k)+dvy_(l,m,k)+vf(l,m,k)*ir(l,m))
             wt(l,m,k) = vg(l,m,k) * (dwx_(l,m,k)+dwy_(l,m,k)+wf(l,m,k)*ir(l,m))
             rt(l,m,k) = vg(l,m,k) * (drx_(l,m,k)+dry_(l,m,k)+rf(l,m,k)*ir(l,m))
          enddo
       enddo
    enddo

    ! Update fluxes at each RK step
    ! =============================
    do k=1,nz
       do i=nxmnghp1,nx
          l=i-nxmngh
          do j=nymnghp1,ny
             cp=cpcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
             av=avcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
             c2_=c2calc_tro(Tmp(i,j,k),rho_n(i,j,k))
          
             m=j-nymngh
             Krho(i,j,k)  = rt(l,m,k)
             Krhou(i,j,k) = uu(i,j,k)*rt(l,m,k)+rho_n(i,j,k)*ut(l,m,k)
             Krhov(i,j,k) = vv(i,j,k)*rt(l,m,k)+rho_n(i,j,k)*vt(l,m,k)
             Krhow(i,j,k) = ww(i,j,k)*rt(l,m,k)+rho_n(i,j,k)*wt(l,m,k)
             Krhoe(i,j,k) = cp/av*(pt(l,m,k)/c2_-rt(l,m,k)) &
                  + (rhoe_n(i,j,k)+prs(i,j,k))/rho_n(i,j,k)*rt(l,m,k) &
                  + rho_n(i,j,k)*(uu(i,j,k)*ut(l,m,k)+vv(i,j,k)*vt(l,m,k)+ww(i,j,k)*wt(l,m,k))
          enddo
       enddo
    enddo

  end subroutine bc_TD2d_imax_jmax_SBP4

end submodule smod_TamDong2d_SBP4_edges
