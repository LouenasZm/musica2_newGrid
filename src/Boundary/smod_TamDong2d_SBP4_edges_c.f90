!===============================================================
submodule (mod_TamDong2d_c) smod_TamDong2d_SBP4_edges_c
!===============================================================
  !> author: XG
  !> date: February 2020 - modif January 2022
  !> Non-reflecting boundary conditions of Tam and Dong
  !> - 2D (periodic) curvilinear version - routines for edges
  !> /!\ dev only: version with SBP4 boundary schemes [not for regular use]
!===============================================================

contains

  !===============================================================================
  module subroutine bc_TD2d_imin_jmin_c_SBP4
  !===============================================================================
    !> 2D Tam & Dong's BC: boundary condition at imin-jmin (edge 1,1,1 /left-bottom) - curvilinear version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: cp,av,c2_
    real(wp), dimension(1:ngh+3,1:ngh+3,nz) :: rf,uf,vf,wf,pf
    real(wp), dimension(ngh,ngh,nz) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(ngh,ngh,nz) :: duksi,dvksi,dwksi,dpksi,drksi
    real(wp), dimension(ngh,ngh,nz) :: dueta,dveta,dweta,dpeta,dreta
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
    
    ! Non-centered derivatives in x-direction
    ! =======================================
    ! (Tam & Webb DRP schemes)
    i=1
    do j=1,ngh
       do k=1,nz
          duksi(i,j,k)= as4p0(1)*uf(1,j,k)+as4p0(2)*uf(2,j,k) &
                      + as4p0(3)*uf(3,j,k)+as4p0(4)*uf(4,j,k)
          dvksi(i,j,k)= as4p0(1)*vf(1,j,k)+as4p0(2)*vf(2,j,k) &
                      + as4p0(3)*vf(3,j,k)+as4p0(4)*vf(4,j,k)
          dwksi(i,j,k)= as4p0(1)*wf(1,j,k)+as4p0(2)*wf(2,j,k) &
                      + as4p0(3)*wf(3,j,k)+as4p0(4)*wf(4,j,k)
          dpksi(i,j,k)= as4p0(1)*pf(1,j,k)+as4p0(2)*pf(2,j,k) &
                      + as4p0(3)*pf(3,j,k)+as4p0(4)*pf(4,j,k)
          drksi(i,j,k)= as4p0(1)*rf(1,j,k)+as4p0(2)*rf(2,j,k) &
                      + as4p0(3)*rf(3,j,k)+as4p0(4)*rf(4,j,k)
       enddo
    enddo

    i=2
    do j=1,ngh
       do k=1,nz
          duksi(i,j,k)= as4p1(1)*uf(1,j,k)+as4p1(2)*uf(2,j,k) &
                      + as4p1(3)*uf(3,j,k)+as4p1(4)*uf(4,j,k)
          dvksi(i,j,k)= as4p1(1)*vf(1,j,k)+as4p1(2)*vf(2,j,k) &
                      + as4p1(3)*vf(3,j,k)+as4p1(4)*vf(4,j,k)
          dwksi(i,j,k)= as4p1(1)*wf(1,j,k)+as4p1(2)*wf(2,j,k) &
                      + as4p1(3)*wf(3,j,k)+as4p1(4)*wf(4,j,k)
          dpksi(i,j,k)= as4p1(1)*pf(1,j,k)+as4p1(2)*pf(2,j,k) &
                      + as4p1(3)*pf(3,j,k)+as4p1(4)*pf(4,j,k)
          drksi(i,j,k)= as4p1(1)*rf(1,j,k)+as4p1(2)*rf(2,j,k) &
                      + as4p1(3)*rf(3,j,k)+as4p1(4)*rf(4,j,k)
       enddo
    enddo

    i=3
    do j=1,ngh
       do k=1,nz
          duksi(i,j,k)= as4p2(1)*uf(1,j,k)+as4p2(2)*uf(2,j,k) &
                      + as4p2(3)*uf(3,j,k)+as4p2(4)*uf(4,j,k) &
                      + as4p2(5)*uf(5,j,k)
          dvksi(i,j,k)= as4p2(1)*vf(1,j,k)+as4p2(2)*vf(2,j,k) &
                      + as4p2(3)*vf(3,j,k)+as4p2(4)*vf(4,j,k) &
                      + as4p2(5)*vf(5,j,k)
          dwksi(i,j,k)= as4p2(1)*wf(1,j,k)+as4p2(2)*wf(2,j,k) &
                      + as4p2(3)*wf(3,j,k)+as4p2(4)*wf(4,j,k) &
                      + as4p2(5)*wf(5,j,k)
          dpksi(i,j,k)= as4p2(1)*pf(1,j,k)+as4p2(2)*pf(2,j,k) &
                      + as4p2(3)*pf(3,j,k)+as4p2(4)*pf(4,j,k) &
                      + as4p2(5)*pf(5,j,k)
          drksi(i,j,k)= as4p2(1)*rf(1,j,k)+as4p2(2)*rf(2,j,k) &
                      + as4p2(3)*rf(3,j,k)+as4p2(4)*rf(4,j,k) &
                      + as4p2(5)*rf(5,j,k)
       enddo
    enddo

    i=4
    do j=1,ngh
       do k=1,nz
          duksi(i,j,k)= as4p3(1)*uf(1,j,k)+as4p3(2)*uf(2,j,k) &
                      + as4p3(3)*uf(3,j,k)+as4p3(4)*uf(4,j,k) &
                      + as4p3(5)*uf(5,j,k)+as4p3(6)*uf(6,j,k)
          dvksi(i,j,k)= as4p3(1)*vf(1,j,k)+as4p3(2)*vf(2,j,k) &
                      + as4p3(3)*vf(3,j,k)+as4p3(4)*vf(4,j,k) &
                      + as4p3(5)*vf(5,j,k)+as4p3(6)*vf(6,j,k)
          dwksi(i,j,k)= as4p3(1)*wf(1,j,k)+as4p3(2)*wf(2,j,k) &
                      + as4p3(3)*wf(3,j,k)+as4p3(4)*wf(4,j,k) &
                      + as4p3(5)*wf(5,j,k)+as4p3(6)*wf(6,j,k)
          dpksi(i,j,k)= as4p3(1)*pf(1,j,k)+as4p3(2)*pf(2,j,k) &
                      + as4p3(3)*pf(3,j,k)+as4p3(4)*pf(4,j,k) &
                      + as4p3(5)*pf(5,j,k)+as4p3(6)*pf(6,j,k)
          drksi(i,j,k)= as4p3(1)*rf(1,j,k)+as4p3(2)*rf(2,j,k) &
                      + as4p3(3)*rf(3,j,k)+as4p3(4)*rf(4,j,k) &
                      + as4p3(5)*rf(5,j,k)+as4p3(6)*rf(6,j,k)
       enddo
    enddo

    do i=5,ngh
       do j=1,ngh
          do k=1,nz
             duksi(i,j,k)= a5(1)* (uf(i+1,j,k)-uf(i-1,j,k)) &
                         + a5(2)* (uf(i+2,j,k)-uf(i-2,j,k))
             dvksi(i,j,k)= a5(1)* (vf(i+1,j,k)-vf(i-1,j,k)) &
                         + a5(2)* (vf(i+2,j,k)-vf(i-2,j,k))
             dwksi(i,j,k)= a5(1)* (wf(i+1,j,k)-wf(i-1,j,k)) &
                         + a5(2)* (wf(i+2,j,k)-wf(i-2,j,k))
             dpksi(i,j,k)= a5(1)* (pf(i+1,j,k)-pf(i-1,j,k)) &
                         + a5(2)* (pf(i+2,j,k)-pf(i-2,j,k))
             drksi(i,j,k)= a5(1)* (rf(i+1,j,k)-rf(i-1,j,k)) &
                         + a5(2)* (rf(i+2,j,k)-rf(i-2,j,k))
          enddo
       enddo
    enddo

    ! Non-centered derivatives in y-direction
    ! =======================================
    j=1
    do i=1,ngh
       do k=1,nz
          dueta(i,j,k)= as4p0(1)*uf(i,1,k)+as4p0(2)*uf(i,2,k) &
                      + as4p0(3)*uf(i,3,k)+as4p0(4)*uf(i,4,k)
          dveta(i,j,k)= as4p0(1)*vf(i,1,k)+as4p0(2)*vf(i,2,k) &
                      + as4p0(3)*vf(i,3,k)+as4p0(4)*vf(i,4,k)
          dweta(i,j,k)= as4p0(1)*wf(i,1,k)+as4p0(2)*wf(i,2,k) &
                      + as4p0(3)*wf(i,3,k)+as4p0(4)*wf(i,4,k)
          dpeta(i,j,k)= as4p0(1)*pf(i,1,k)+as4p0(2)*pf(i,2,k) &
                      + as4p0(3)*pf(i,3,k)+as4p0(4)*pf(i,4,k)
          dreta(i,j,k)= as4p0(1)*rf(i,1,k)+as4p0(2)*rf(i,2,k) &
                      + as4p0(3)*rf(i,3,k)+as4p0(4)*rf(i,4,k)
       enddo
    enddo

    j=2
    do i=1,ngh
       do k=1,nz
          dueta(i,j,k)= as4p1(1)*uf(i,1,k)+as4p1(2)*uf(i,2,k) &
                      + as4p1(3)*uf(i,3,k)+as4p1(4)*uf(i,4,k)
          dveta(i,j,k)= as4p1(1)*vf(i,1,k)+as4p1(2)*vf(i,2,k) &
                      + as4p1(3)*vf(i,3,k)+as4p1(4)*vf(i,4,k)
          dweta(i,j,k)= as4p1(1)*wf(i,1,k)+as4p1(2)*wf(i,2,k) &
                      + as4p1(3)*wf(i,3,k)+as4p1(4)*wf(i,4,k)
          dpeta(i,j,k)= as4p1(1)*pf(i,1,k)+as4p1(2)*pf(i,2,k) &
                      + as4p1(3)*pf(i,3,k)+as4p1(4)*pf(i,4,k)
          dreta(i,j,k)= as4p1(1)*rf(i,1,k)+as4p1(2)*rf(i,2,k) &
                      + as4p1(3)*rf(i,3,k)+as4p1(4)*rf(i,4,k)
       enddo
    enddo

    j=3
    do i=1,ngh
       do k=1,nz
          dueta(i,j,k)= as4p2(1)*uf(i,1,k)+as4p2(2)*uf(i,2,k) &
                      + as4p2(3)*uf(i,3,k)+as4p2(4)*uf(i,4,k) &
                      + as4p2(5)*uf(i,5,k)
          dveta(i,j,k)= as4p2(1)*vf(i,1,k)+as4p2(2)*vf(i,2,k) &
                      + as4p2(3)*vf(i,3,k)+as4p2(4)*vf(i,4,k) &
                      + as4p2(5)*vf(i,5,k)
          dweta(i,j,k)= as4p2(1)*wf(i,1,k)+as4p2(2)*wf(i,2,k) &
                      + as4p2(3)*wf(i,3,k)+as4p2(4)*wf(i,4,k) &
                      + as4p2(5)*wf(i,5,k)
          dpeta(i,j,k)= as4p2(1)*pf(i,1,k)+as4p2(2)*pf(i,2,k) &
                      + as4p2(3)*pf(i,3,k)+as4p2(4)*pf(i,4,k) &
                      + as4p2(5)*pf(i,5,k)
          dreta(i,j,k)= as4p2(1)*rf(i,1,k)+as4p2(2)*rf(i,2,k) &
                      + as4p2(3)*rf(i,3,k)+as4p2(4)*rf(i,4,k) &
                      + as4p2(5)*rf(i,5,k)
       enddo
    enddo
    
    j=4
    do i=1,ngh
       do k=1,nz
          dueta(i,j,k)= as4p3(1)*uf(i,1,k)+as4p3(2)*uf(i,2,k) &
                      + as4p3(3)*uf(i,3,k)+as4p3(4)*uf(i,4,k) &
                      + as4p3(5)*uf(i,5,k)+as4p3(6)*uf(i,6,k)
          dveta(i,j,k)= as4p3(1)*vf(i,1,k)+as4p3(2)*vf(i,2,k) &
                      + as4p3(3)*vf(i,3,k)+as4p3(4)*vf(i,4,k) &
                      + as4p3(5)*vf(i,5,k)+as4p3(6)*vf(i,6,k)
          dweta(i,j,k)= as4p3(1)*wf(i,1,k)+as4p3(2)*wf(i,2,k) &
                      + as4p3(3)*wf(i,3,k)+as4p3(4)*wf(i,4,k) &
                      + as4p3(5)*wf(i,5,k)+as4p3(6)*wf(i,6,k)
          dpeta(i,j,k)= as4p3(1)*pf(i,1,k)+as4p3(2)*pf(i,2,k) &
                      + as4p3(3)*pf(i,3,k)+as4p3(4)*pf(i,4,k) &
                      + as4p3(5)*pf(i,5,k)+as4p3(6)*pf(i,6,k)
          dreta(i,j,k)= as4p3(1)*rf(i,1,k)+as4p3(2)*rf(i,2,k) &
                      + as4p3(3)*rf(i,3,k)+as4p3(4)*rf(i,4,k) &
                      + as4p3(5)*rf(i,5,k)+as4p3(6)*rf(i,6,k)
       enddo
    enddo
    
    do j=5,ngh
       do i=1,ngh
          do k=1,nz
             dueta(i,j,k)= a5(1)*(uf(i,j+1,k)-uf(i,j-1,k)) &
                         + a5(2)*(uf(i,j+2,k)-uf(i,j-2,k))
             dveta(i,j,k)= a5(1)*(vf(i,j+1,k)-vf(i,j-1,k)) &
                         + a5(2)*(vf(i,j+2,k)-vf(i,j-2,k))
             dweta(i,j,k)= a5(1)*(wf(i,j+1,k)-wf(i,j-1,k)) &
                         + a5(2)*(wf(i,j+2,k)-wf(i,j-2,k))
             dpeta(i,j,k)= a5(1)*(pf(i,j+1,k)-pf(i,j-1,k)) &
                         + a5(2)*(pf(i,j+2,k)-pf(i,j-2,k))
             dreta(i,j,k)= a5(1)*(rf(i,j+1,k)-rf(i,j-1,k)) &
                         + a5(2)*(rf(i,j+2,k)-rf(i,j-2,k))
          enddo
       enddo
    enddo

    ! Compute rt,ut,vt,wt & pt
    ! ========================
    do i=1,ngh
       do j=1,ngh
          do k=1,nz
             pt(i,j,k) = vg(i,j,k)*( ((dpksi(i,j,k)*y_eta(i,j)-dpeta(i,j,k)*y_ksi(i,j))*cosphi(i,j) &
                                     +(dpeta(i,j,k)*x_ksi(i,j)-dpksi(i,j,k)*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                     + pf(i,j,k)*ir(i,j) )
             ut(i,j,k) = vg(i,j,k)*( ((duksi(i,j,k)*y_eta(i,j)-dueta(i,j,k)*y_ksi(i,j))*cosphi(i,j) &
                                     +(dueta(i,j,k)*x_ksi(i,j)-duksi(i,j,k)*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                     + uf(i,j,k)*ir(i,j) )
             vt(i,j,k) = vg(i,j,k)*( ((dvksi(i,j,k)*y_eta(i,j)-dveta(i,j,k)*y_ksi(i,j))*cosphi(i,j) &
                                     +(dveta(i,j,k)*x_ksi(i,j)-dvksi(i,j,k)*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                     + vf(i,j,k)*ir(i,j) )
             wt(i,j,k) = vg(i,j,k)*( ((dwksi(i,j,k)*y_eta(i,j)-dweta(i,j,k)*y_ksi(i,j))*cosphi(i,j) &
                                     +(dweta(i,j,k)*x_ksi(i,j)-dwksi(i,j,k)*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                     + wf(i,j,k)*ir(i,j) )
             rt(i,j,k) = vg(i,j,k)*( ((drksi(i,j,k)*y_eta(i,j)-dreta(i,j,k)*y_ksi(i,j))*cosphi(i,j) &
                                     +(dreta(i,j,k)*x_ksi(i,j)-drksi(i,j,k)*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                     + rf(i,j,k)*ir(i,j) )
          enddo
       enddo
    enddo

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

  end subroutine bc_TD2d_imin_jmin_c_SBP4

  !===============================================================================
  module subroutine bc_TD2d_imin_jmax_c_SBP4
  !===============================================================================
    !> 2D Tam & Dong's BC: boundary condition at imin-jmax (edge 1,1,2 /left-top) - curvilinear version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l
    real(wp) :: cp,av,c2_
    real(wp), dimension(1:ngh+3,-2:ngh,nz) :: rf,uf,vf,wf,pf
    real(wp), dimension(ngh,ngh,nz)  :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(ngh,ngh,nz) :: duksi,dvksi,dwksi,dpksi,drksi
    real(wp), dimension(ngh,ngh,nz) :: dueta,dveta,dweta,dpeta,dreta
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

    ! Non-centered derivatives in x-direction
    ! =======================================
    ! (Tam & Webb DRP schemes)
    i=1
    do j=nymnghp1,ny
       l=j-nymngh
       do k=1,nz
          duksi(i,l,k)= as4p0(1)*uf(1,l,k)+as4p0(2)*uf(2,l,k) &
                      + as4p0(3)*uf(3,l,k)+as4p0(4)*uf(4,l,k)
          dvksi(i,l,k)= as4p0(1)*vf(1,l,k)+as4p0(2)*vf(2,l,k) &
                      + as4p0(3)*vf(3,l,k)+as4p0(4)*vf(4,l,k)
          dwksi(i,l,k)= as4p0(1)*wf(1,l,k)+as4p0(2)*wf(2,l,k) &
                      + as4p0(3)*wf(3,l,k)+as4p0(4)*wf(4,l,k)
          dpksi(i,l,k)= as4p0(1)*pf(1,l,k)+as4p0(2)*pf(2,l,k) &
                      + as4p0(3)*pf(3,l,k)+as4p0(4)*pf(4,l,k)
          drksi(i,l,k)= as4p0(1)*rf(1,l,k)+as4p0(2)*rf(2,l,k) &
                      + as4p0(3)*rf(3,l,k)+as4p0(4)*rf(4,l,k)
       enddo
    enddo

    i=2
    do j=nymnghp1,ny
       l=j-nymngh
       do k=1,nz
          duksi(i,l,k)= as4p1(1)*uf(1,l,k)+as4p1(2)*uf(2,l,k) &
                      + as4p1(3)*uf(3,l,k)+as4p1(4)*uf(4,l,k)
          dvksi(i,l,k)= as4p1(1)*vf(1,l,k)+as4p1(2)*vf(2,l,k) &
                      + as4p1(3)*vf(3,l,k)+as4p1(4)*vf(4,l,k)
          dwksi(i,l,k)= as4p1(1)*wf(1,l,k)+as4p1(2)*wf(2,l,k) &
                      + as4p1(3)*wf(3,l,k)+as4p1(4)*wf(4,l,k)
          dpksi(i,l,k)= as4p1(1)*pf(1,l,k)+as4p1(2)*pf(2,l,k) &
                      + as4p1(3)*pf(3,l,k)+as4p1(4)*pf(4,l,k)
          drksi(i,l,k)= as4p1(1)*rf(1,l,k)+as4p1(2)*rf(2,l,k) &
                      + as4p1(3)*rf(3,l,k)+as4p1(4)*rf(4,l,k)
       enddo
    enddo

    i=3
    do j=nymnghp1,ny
       l=j-nymngh
       do k=1,nz
          duksi(i,l,k)= as4p2(1)*uf(1,l,k)+as4p2(2)*uf(2,l,k) &
                      + as4p2(3)*uf(3,l,k)+as4p2(4)*uf(4,l,k) &
                      + as4p2(5)*uf(5,l,k)
          dvksi(i,l,k)= as4p2(1)*vf(1,l,k)+as4p2(2)*vf(2,l,k) &
                      + as4p2(3)*vf(3,l,k)+as4p2(4)*vf(4,l,k) &
                      + as4p2(5)*vf(5,l,k)
          dwksi(i,l,k)= as4p2(1)*wf(1,l,k)+as4p2(2)*wf(2,l,k) &
                      + as4p2(3)*wf(3,l,k)+as4p2(4)*wf(4,l,k) &
                      + as4p2(5)*wf(5,l,k)
          dpksi(i,l,k)= as4p2(1)*pf(1,l,k)+as4p2(2)*pf(2,l,k) &
                      + as4p2(3)*pf(3,l,k)+as4p2(4)*pf(4,l,k) &
                      + as4p2(5)*pf(5,l,k)
          drksi(i,l,k)= as4p2(1)*rf(1,l,k)+as4p2(2)*rf(2,l,k) &
                      + as4p2(3)*rf(3,l,k)+as4p2(4)*rf(4,l,k) &
                      + as4p2(5)*rf(5,l,k)
       enddo
    enddo

    i=4
    do j=nymnghp1,ny
       l=j-nymngh
       do k=1,nz
          duksi(i,l,k)= as4p3(1)*uf(1,l,k)+as4p3(2)*uf(2,l,k) &
                      + as4p3(3)*uf(3,l,k)+as4p3(4)*uf(4,l,k) &
                      + as4p3(5)*uf(5,l,k)+as4p3(6)*uf(6,l,k)
          dvksi(i,l,k)= as4p3(1)*vf(1,l,k)+as4p3(2)*vf(2,l,k) &
                      + as4p3(3)*vf(3,l,k)+as4p3(4)*vf(4,l,k) &
                      + as4p3(5)*vf(5,l,k)+as4p3(6)*vf(6,l,k)
          dwksi(i,l,k)= as4p3(1)*wf(1,l,k)+as4p3(2)*wf(2,l,k) &
                      + as4p3(3)*wf(3,l,k)+as4p3(4)*wf(4,l,k) &
                      + as4p3(5)*wf(5,l,k)+as4p3(6)*wf(6,l,k)
          dpksi(i,l,k)= as4p3(1)*pf(1,l,k)+as4p3(2)*pf(2,l,k) &
                      + as4p3(3)*pf(3,l,k)+as4p3(4)*pf(4,l,k) &
                      + as4p3(5)*pf(5,l,k)+as4p3(6)*pf(6,l,k)
          drksi(i,l,k)= as4p3(1)*rf(1,l,k)+as4p3(2)*rf(2,l,k) &
                      + as4p3(3)*rf(3,l,k)+as4p3(4)*rf(4,l,k) &
                      + as4p3(5)*rf(5,l,k)+as4p3(6)*rf(6,l,k)
       enddo
    enddo

    do i=5,ngh
       do j=nymnghp1,ny
          l=j-nymngh
          do k=1,nz
             duksi(i,l,k)= a5(1)* (uf(i+1,l,k)-uf(i-1,l,k)) &
                         + a5(2)* (uf(i+2,l,k)-uf(i-2,l,k))
             dvksi(i,l,k)= a5(1)* (vf(i+1,l,k)-vf(i-1,l,k)) &
                         + a5(2)* (vf(i+2,l,k)-vf(i-2,l,k))
             dwksi(i,l,k)= a5(1)* (wf(i+1,l,k)-wf(i-1,l,k)) &
                         + a5(2)* (wf(i+2,l,k)-wf(i-2,l,k))
             dpksi(i,l,k)= a5(1)* (pf(i+1,l,k)-pf(i-1,l,k)) &
                         + a5(2)* (pf(i+2,l,k)-pf(i-2,l,k))
             drksi(i,l,k)= a5(1)* (rf(i+1,l,k)-rf(i-1,l,k)) &
                         + a5(2)* (rf(i+2,l,k)-rf(i-2,l,k))
          enddo
       enddo
    enddo

    ! Non-centered derivatives in y-direction
    ! =======================================
    ! (Tam & Webb DRP schemes)
    do j=nymnghp1,ny-4
       l=j-nymngh
       do i=1,ngh
          do k=1,nz
             dueta(i,l,k) = a5(1)*(uf(i,l+1,k)-uf(i,l-1,k)) &
                          + a5(2)*(uf(i,l+2,k)-uf(i,l-2,k))
             dveta(i,l,k) = a5(1)*(vf(i,l+1,k)-vf(i,l-1,k)) &
                          + a5(2)*(vf(i,l+2,k)-vf(i,l-2,k))
             dweta(i,l,k) = a5(1)*(wf(i,l+1,k)-wf(i,l-1,k)) &
                          + a5(2)*(wf(i,l+2,k)-wf(i,l-2,k))
             dpeta(i,l,k) = a5(1)*(pf(i,l+1,k)-pf(i,l-1,k)) &
                          + a5(2)*(pf(i,l+2,k)-pf(i,l-2,k))
             dreta(i,l,k) = a5(1)*(rf(i,l+1,k)-rf(i,l-1,k)) &
                          + a5(2)*(rf(i,l+2,k)-rf(i,l-2,k))
          enddo
       enddo
    enddo

    j=ny-3
    l=j-nymngh
    do i=1,ngh
       do k=1,nz
          dueta(i,l,k) = as4m3(1)*uf(i,l+3,k)+as4m3(2)*uf(i,l+2,k) &
                       + as4m3(3)*uf(i,l+1,k)+as4m3(4)*uf(i,l  ,k) &
                       + as4m3(5)*uf(i,l-1,k)+as4m3(6)*uf(i,l-2,k)
          dveta(i,l,k) = as4m3(1)*vf(i,l+3,k)+as4m3(2)*vf(i,l+2,k) &
                       + as4m3(3)*vf(i,l+1,k)+as4m3(4)*vf(i,l  ,k) &
                       + as4m3(5)*vf(i,l-1,k)+as4m3(6)*vf(i,l-2,k)
          dweta(i,l,k) = as4m3(1)*wf(i,l+3,k)+as4m3(2)*wf(i,l+2,k) &
                       + as4m3(3)*wf(i,l+1,k)+as4m3(4)*wf(i,l  ,k) &
                       + as4m3(5)*wf(i,l-1,k)+as4m3(6)*wf(i,l-2,k)
          dpeta(i,l,k) = as4m3(1)*pf(i,l+3,k)+as4m3(2)*pf(i,l+2,k) &
                       + as4m3(3)*pf(i,l+1,k)+as4m3(4)*pf(i,l  ,k) &
                       + as4m3(5)*pf(i,l-1,k)+as4m3(6)*pf(i,l-2,k)
          dreta(i,l,k) = as4m3(1)*rf(i,l+3,k)+as4m3(2)*rf(i,l+2,k) &
                       + as4m3(3)*rf(i,l+1,k)+as4m3(4)*rf(i,l  ,k) &
                       + as4m3(5)*rf(i,l-1,k)+as4m3(6)*rf(i,l-2,k)
       enddo
    enddo
    
    j=ny-2
    l=j-nymngh
    do i=1,ngh
       do k=1,nz
          dueta(i,l,k) = as4m2(1)*uf(i,l+2,k)+as4m2(2)*uf(i,l+1,k) &
                       + as4m2(3)*uf(i,l  ,k)+as4m2(4)*uf(i,l-1,k) &
                       + as4m2(5)*uf(i,l-2,k)
          dveta(i,l,k) = as4m2(1)*vf(i,l+2,k)+as4m2(2)*vf(i,l+1,k) &
                       + as4m2(3)*vf(i,l  ,k)+as4m2(4)*vf(i,l-1,k) &
                       + as4m2(5)*vf(i,l-2,k)
          dweta(i,l,k) = as4m2(1)*wf(i,l+2,k)+as4m2(2)*wf(i,l+1,k) &
                       + as4m2(3)*wf(i,l  ,k)+as4m2(4)*wf(i,l-1,k) &
                       + as4m2(5)*wf(i,l-2,k)
          dpeta(i,l,k) = as4m2(1)*pf(i,l+2,k)+as4m2(2)*pf(i,l+1,k) &
                       + as4m2(3)*pf(i,l  ,k)+as4m2(4)*pf(i,l-1,k) &
                       + as4m2(5)*pf(i,l-2,k)
          dreta(i,l,k) = as4m2(1)*rf(i,l+2,k)+as4m2(2)*rf(i,l+1,k) &
                       + as4m2(3)*rf(i,l  ,k)+as4m2(4)*rf(i,l-1,k) &
                       + as4m2(5)*rf(i,l-2,k)
       enddo
    enddo
    
    j=ny-1
    l=j-nymngh
    do i=1,ngh
       do k=1,nz
          dueta(i,l,k) = as4m1(1)*uf(i,l+1,k)+as4m1(2)*uf(i,l  ,k) &
                       + as4m1(3)*uf(i,l-1,k)+as4m1(4)*uf(i,l-2,k)
          dveta(i,l,k) = as4m1(1)*vf(i,l+1,k)+as4m1(2)*vf(i,l  ,k) &
                       + as4m1(3)*vf(i,l-1,k)+as4m1(4)*vf(i,l-2,k)
          dweta(i,l,k) = as4m1(1)*wf(i,l+1,k)+as4m1(2)*wf(i,l  ,k) &
                       + as4m1(3)*wf(i,l-1,k)+as4m1(4)*wf(i,l-2,k)
          dpeta(i,l,k) = as4m1(1)*pf(i,l+1,k)+as4m1(2)*pf(i,l  ,k) &
                       + as4m1(3)*pf(i,l-1,k)+as4m1(4)*pf(i,l-2,k)
          dreta(i,l,k) = as4m1(1)*rf(i,l+1,k)+as4m1(2)*rf(i,l  ,k) &
                       + as4m1(3)*rf(i,l-1,k)+as4m1(4)*rf(i,l-2,k)
       enddo
    enddo
    
    j=ny
    l=j-nymngh
    do i=1,ngh
       do k=1,nz
          dueta(i,l,k) = as4m0(1)*uf(i,l  ,k)+as4m0(2)*uf(i,l-1,k) &
                       + as4m0(3)*uf(i,l-2,k)+as4m0(4)*uf(i,l-3,k)
          dveta(i,l,k) = as4m0(1)*vf(i,l  ,k)+as4m0(2)*vf(i,l-1,k) &
                       + as4m0(3)*vf(i,l-2,k)+as4m0(4)*vf(i,l-3,k)
          dweta(i,l,k) = as4m0(1)*wf(i,l  ,k)+as4m0(2)*wf(i,l-1,k) &
                       + as4m0(3)*wf(i,l-2,k)+as4m0(4)*wf(i,l-3,k)
          dpeta(i,l,k) = as4m0(1)*pf(i,l  ,k)+as4m0(2)*pf(i,l-1,k) &
                       + as4m0(3)*pf(i,l-2,k)+as4m0(4)*pf(i,l-3,k)
          dreta(i,l,k) = as4m0(1)*rf(i,l  ,k)+as4m0(2)*rf(i,l-1,k) &
                       + as4m0(3)*rf(i,l-2,k)+as4m0(4)*rf(i,l-3,k)
       enddo
    enddo
         
    ! Compute rt,ut,vt,wt & pt
    ! ========================
    do i=1,ngh
       do j=nymnghp1,ny
          l=j-nymngh
          do k=1,nz
             pt(i,l,k) = vg(i,l,k)*( ((dpksi(i,l,k)*y_eta(i,j)-dpeta(i,l,k)*y_ksi(i,j))*cosphi(i,l) &
                                     +(dpeta(i,l,k)*x_ksi(i,j)-dpksi(i,l,k)*x_eta(i,j))*sinphi(i,l))*ijacob(i,j) &
                                     + pf(i,l,k)*ir(i,l) )
             ut(i,l,k) = vg(i,l,k)*( ((duksi(i,l,k)*y_eta(i,j)-dueta(i,l,k)*y_ksi(i,j))*cosphi(i,l) &
                                     +(dueta(i,l,k)*x_ksi(i,j)-duksi(i,l,k)*x_eta(i,j))*sinphi(i,l))*ijacob(i,j) &
                                     + uf(i,l,k)*ir(i,l) )
             vt(i,l,k) = vg(i,l,k)*( ((dvksi(i,l,k)*y_eta(i,j)-dveta(i,l,k)*y_ksi(i,j))*cosphi(i,l) &
                                     +(dveta(i,l,k)*x_ksi(i,j)-dvksi(i,l,k)*x_eta(i,j))*sinphi(i,l))*ijacob(i,j) &
                                     + vf(i,l,k)*ir(i,l) )
             wt(i,l,k) = vg(i,l,k)*( ((dwksi(i,l,k)*y_eta(i,j)-dweta(i,l,k)*y_ksi(i,j))*cosphi(i,l) &
                                     +(dweta(i,l,k)*x_ksi(i,j)-dwksi(i,l,k)*x_eta(i,j))*sinphi(i,l))*ijacob(i,j) &
                                     + wf(i,l,k)*ir(i,l) )
             rt(i,l,k) = vg(i,l,k)*( ((drksi(i,l,k)*y_eta(i,j)-dreta(i,l,k)*y_ksi(i,j))*cosphi(i,l) &
                                     +(dreta(i,l,k)*x_ksi(i,j)-drksi(i,l,k)*x_eta(i,j))*sinphi(i,l))*ijacob(i,j) &
                                     + rf(i,l,k)*ir(i,l) )
          enddo
       enddo
    enddo

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

  end subroutine bc_TD2d_imin_jmax_c_SBP4

  !===============================================================================
  module subroutine bc_TD2d_imax_jmin_c_SBP4
  !===============================================================================
    !> 2D Tam & Dong's BC: boundary condition at imax-jmin (edge 1,2,1 /right-bottom) - curvilinear version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l
    real(wp) :: cp,av,c2_
    real(wp), dimension(-2:ngh,1:ngh+3,nz) :: rf,uf,vf,wf,pf
    real(wp), dimension(ngh,ngh,nz)  :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(ngh,ngh,nz)  :: duksi,dvksi,dwksi,dpksi,drksi
    real(wp), dimension(ngh,ngh,nz)  :: dueta,dveta,dweta,dpeta,dreta
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
!!$             rf(l,j,k)=rho_n(i,j,k)-BC_face(2,1)%U0(i,j,k,1)
!!$             uf(l,j,k)=   uu(i,j,k)-BC_face(2,1)%U0(i,j,k,2)
!!$             vf(l,j,k)=   vv(i,j,k)-BC_face(2,1)%U0(i,j,k,3)
!!$             wf(l,j,k)=   ww(i,j,k)-BC_face(2,1)%U0(i,j,k,4)
!!$             pf(l,j,k)=  prs(i,j,k)-BC_face(2,1)%U0(i,j,k,5)
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
    do i=nxmnghp1,nx
       l=i-nxmngh
       do j=1,ngh
          do k=1,nz
!!$             vg(l,j,k)= BC_face(2,1)%U0(i,j,k,2)*cosphi(l,j)+BC_face(2,1)%U0(i,j,k,3)*sinphi(l,j) &
!!$                  + sqrt(BC_face(2,1)%U0(i,j,k,6)-BC_face(2,1)%U0(i,j,k,4)**2 &
!!$                  -(BC_face(2,1)%U0(i,j,k,2)*sinphi(l,j)-BC_face(2,1)%U0(i,j,k,3)*cosphi(l,j))**2)
             vg(l,j,k)= BC_face(1,2)%U0(l,j,k,2)*cosphi(l,j)+BC_face(1,2)%U0(l,j,k,3)*sinphi(l,j) &
                  + sqrt(BC_face(1,2)%U0(l,j,k,6)-BC_face(1,2)%U0(l,j,k,4)**2 &
                  -(BC_face(1,2)%U0(l,j,k,2)*sinphi(l,j)-BC_face(1,2)%U0(l,j,k,3)*cosphi(l,j))**2)
          enddo
       enddo
    enddo

    ! Non-centered derivatives in x-direction
    ! =======================================
    ! (Tam & Webb DRP schemes)
    do i=nxmnghp1,nx-4
       l=i-nxmngh
       do j=1,ngh
          do k=1,nz
             duksi(l,j,k) = a5(1)* (uf(l+1,j,k)-uf(l-1,j,k)) &
                          + a5(2)* (uf(l+2,j,k)-uf(l-2,j,k))
             dvksi(l,j,k) = a5(1)* (vf(l+1,j,k)-vf(l-1,j,k)) &
                          + a5(2)* (vf(l+2,j,k)-vf(l-2,j,k))
             dwksi(l,j,k) = a5(1)* (wf(l+1,j,k)-wf(l-1,j,k)) &
                          + a5(2)* (wf(l+2,j,k)-wf(l-2,j,k))
             dpksi(l,j,k) = a5(1)* (pf(l+1,j,k)-pf(l-1,j,k)) &
                          + a5(2)* (pf(l+2,j,k)-pf(l-2,j,k))
             drksi(l,j,k) = a5(1)* (rf(l+1,j,k)-rf(l-1,j,k)) &
                          + a5(2)* (rf(l+2,j,k)-rf(l-2,j,k))
          enddo
       enddo
    enddo
    
    i=nx-3
    l=i-nxmngh
    do j=1,ngh
       do k=1,nz
          duksi(l,j,k) = as4m3(1)*uf(l+3,j,k)+as4m3(2)*uf(l+2,j,k) &
                       + as4m3(3)*uf(l+1,j,k)+as4m3(4)*uf(l  ,j,k) &
                       + as4m3(5)*uf(l-1,j,k)+as4m3(6)*uf(l-2,j,k)
          dvksi(l,j,k) = as4m3(1)*vf(l+3,j,k)+as4m3(2)*vf(l+2,j,k) &
                       + as4m3(3)*vf(l+1,j,k)+as4m3(4)*vf(l  ,j,k) &
                       + as4m3(5)*vf(l-1,j,k)+as4m3(6)*vf(l-2,j,k)
          dwksi(l,j,k) = as4m3(1)*wf(l+3,j,k)+as4m3(2)*wf(l+2,j,k) &
                       + as4m3(3)*wf(l+1,j,k)+as4m3(4)*wf(l  ,j,k) &
                       + as4m3(5)*wf(l-1,j,k)+as4m3(6)*wf(l-2,j,k)
          dpksi(l,j,k) = as4m3(1)*pf(l+3,j,k)+as4m3(2)*pf(l+2,j,k) &
                       + as4m3(3)*pf(l+1,j,k)+as4m3(4)*pf(l  ,j,k) &
                       + as4m3(5)*pf(l-1,j,k)+as4m3(6)*pf(l-2,j,k)
          drksi(l,j,k) = as4m3(1)*rf(l+3,j,k)+as4m3(2)*rf(l+2,j,k) &
                       + as4m3(3)*rf(l+1,j,k)+as4m3(4)*rf(l  ,j,k) &
                       + as4m3(5)*rf(l-1,j,k)+as4m3(6)*rf(l-2,j,k)
       enddo
    enddo
    
    i=nx-2
    l=i-nxmngh
    do j=1,ngh
       do k=1,nz
          duksi(l,j,k) = as4m2(1)*uf(l+2,j,k)+as4m2(2)*uf(l+1,j,k) &
                       + as4m2(3)*uf(l  ,j,k)+as4m2(4)*uf(l-1,j,k) &
                       + as4m2(5)*uf(l-2,j,k)
          dvksi(l,j,k) = as4m2(1)*vf(l+2,j,k)+as4m2(2)*vf(l+1,j,k) &
                       + as4m2(3)*vf(l  ,j,k)+as4m2(4)*vf(l-1,j,k) &
                       + as4m2(5)*vf(l-2,j,k)
          dwksi(l,j,k) = as4m2(1)*wf(l+2,j,k)+as4m2(2)*wf(l+1,j,k) &
                       + as4m2(3)*wf(l  ,j,k)+as4m2(4)*wf(l-1,j,k) &
                       + as4m2(5)*wf(l-2,j,k)
          dpksi(l,j,k) = as4m2(1)*pf(l+2,j,k)+as4m2(2)*pf(l+1,j,k) &
                       + as4m2(3)*pf(l  ,j,k)+as4m2(4)*pf(l-1,j,k) &
                       + as4m2(5)*pf(l-2,j,k)
          drksi(l,j,k) = as4m2(1)*rf(l+2,j,k)+as4m2(2)*rf(l+1,j,k) &
                       + as4m2(3)*rf(l  ,j,k)+as4m2(4)*rf(l-1,j,k) &
                       + as4m2(5)*rf(l-2,j,k)
       enddo
    enddo
    
    i=nx-1
    l=i-nxmngh
    do j=1,ngh
       do k=1,nz
          duksi(l,j,k) = as4m1(1)*uf(l+1,j,k)+as4m1(2)*uf(l  ,j,k) &
                       + as4m1(3)*uf(l-1,j,k)+as4m1(4)*uf(l-2,j,k)
          dvksi(l,j,k) = as4m1(1)*vf(l+1,j,k)+as4m1(2)*vf(l  ,j,k) &
                       + as4m1(3)*vf(l-1,j,k)+as4m1(4)*vf(l-2,j,k)
          dwksi(l,j,k) = as4m1(1)*wf(l+1,j,k)+as4m1(2)*wf(l  ,j,k) &
                       + as4m1(3)*wf(l-1,j,k)+as4m1(4)*wf(l-2,j,k)
          dpksi(l,j,k) = as4m1(1)*pf(l+1,j,k)+as4m1(2)*pf(l  ,j,k) &
                       + as4m1(3)*pf(l-1,j,k)+as4m1(4)*pf(l-2,j,k)
          drksi(l,j,k) = as4m1(1)*rf(l+1,j,k)+as4m1(2)*rf(l  ,j,k) &
                       + as4m1(3)*rf(l-1,j,k)+as4m1(4)*rf(l-2,j,k)
       enddo
    enddo
    
    i=nx
    l=i-nxmngh
    do j=1,ngh
       do k=1,nz
          duksi(l,j,k) = as4m0(1)*uf(l  ,j,k)+as4m0(2)*uf(l-1,j,k) &
                       + as4m0(3)*uf(l-2,j,k)+as4m0(4)*uf(l-3,j,k)
          dvksi(l,j,k) = as4m0(1)*vf(l  ,j,k)+as4m0(2)*vf(l-1,j,k) &
                       + as4m0(3)*vf(l-2,j,k)+as4m0(4)*vf(l-3,j,k)
          dwksi(l,j,k) = as4m0(1)*wf(l  ,j,k)+as4m0(2)*wf(l-1,j,k) &
                       + as4m0(3)*wf(l-2,j,k)+as4m0(4)*wf(l-3,j,k)
          dpksi(l,j,k) = as4m0(1)*pf(l  ,j,k)+as4m0(2)*pf(l-1,j,k) &
                       + as4m0(3)*pf(l-2,j,k)+as4m0(4)*pf(l-3,j,k)
          drksi(l,j,k) = as4m0(1)*rf(l  ,j,k)+as4m0(2)*rf(l-1,j,k) &
                       + as4m0(3)*rf(l-2,j,k)+as4m0(4)*rf(l-3,j,k)
       enddo
    enddo

    ! Non-centered derivatives in y-direction
    ! =======================================
    j=1
    do i=nxmnghp1,nx
       l=i-nxmngh
       do k=1,nz
          dueta(l,j,k)= as4p0(1)*uf(l,1,k)+as4p0(2)*uf(l,2,k) &
                      + as4p0(3)*uf(l,3,k)+as4p0(4)*uf(l,4,k)
          dveta(l,j,k)= as4p0(1)*vf(l,1,k)+as4p0(2)*vf(l,2,k) &
                      + as4p0(3)*vf(l,3,k)+as4p0(4)*vf(l,4,k)
          dweta(l,j,k)= as4p0(1)*wf(l,1,k)+as4p0(2)*wf(l,2,k) &
                      + as4p0(3)*wf(l,3,k)+as4p0(4)*wf(l,4,k)
          dpeta(l,j,k)= as4p0(1)*pf(l,1,k)+as4p0(2)*pf(l,2,k) &
                      + as4p0(3)*pf(l,3,k)+as4p0(4)*pf(l,4,k)
          dreta(l,j,k)= as4p0(1)*rf(l,1,k)+as4p0(2)*rf(l,2,k) &
                      + as4p0(3)*rf(l,3,k)+as4p0(4)*rf(l,4,k)
       enddo
    enddo

    j=2
    do i=nxmnghp1,nx
       l=i-nxmngh
       do k=1,nz
          dueta(l,j,k)= as4p1(1)*uf(l,1,k)+as4p1(2)*uf(l,2,k) &
                      + as4p1(3)*uf(l,3,k)+as4p1(4)*uf(l,4,k)
          dveta(l,j,k)= as4p1(1)*vf(l,1,k)+as4p1(2)*vf(l,2,k) &
                      + as4p1(3)*vf(l,3,k)+as4p1(4)*vf(l,4,k)
          dweta(l,j,k)= as4p1(1)*wf(l,1,k)+as4p1(2)*wf(l,2,k) &
                      + as4p1(3)*wf(l,3,k)+as4p1(4)*wf(l,4,k)
          dpeta(l,j,k)= as4p1(1)*pf(l,1,k)+as4p1(2)*pf(l,2,k) &
                      + as4p1(3)*pf(l,3,k)+as4p1(4)*pf(l,4,k)
          dreta(l,j,k)= as4p1(1)*rf(l,1,k)+as4p1(2)*rf(l,2,k) &
                      + as4p1(3)*rf(l,3,k)+as4p1(4)*rf(l,4,k)
       enddo
    enddo

    j=3
    do i=nxmnghp1,nx
       l=i-nxmngh
       do k=1,nz
          dueta(l,j,k)= as4p2(1)*uf(l,1,k)+as4p2(2)*uf(l,2,k) &
                      + as4p2(3)*uf(l,3,k)+as4p2(4)*uf(l,4,k) &
                      + as4p2(5)*uf(l,5,k)
          dveta(l,j,k)= as4p2(1)*vf(l,1,k)+as4p2(2)*vf(l,2,k) &
                      + as4p2(3)*vf(l,3,k)+as4p2(4)*vf(l,4,k) &
                      + as4p2(5)*vf(l,5,k)
          dweta(l,j,k)= as4p2(1)*wf(l,1,k)+as4p2(2)*wf(l,2,k) &
                      + as4p2(3)*wf(l,3,k)+as4p2(4)*wf(l,4,k) &
                      + as4p2(5)*wf(l,5,k)
          dpeta(l,j,k)= as4p2(1)*pf(l,1,k)+as4p2(2)*pf(l,2,k) &
                      + as4p2(3)*pf(l,3,k)+as4p2(4)*pf(l,4,k) &
                      + as4p2(5)*pf(l,5,k)
          dreta(l,j,k)= as4p2(1)*rf(l,1,k)+as4p2(2)*rf(l,2,k) &
                      + as4p2(3)*rf(l,3,k)+as4p2(4)*rf(l,4,k) &
                      + as4p2(5)*rf(l,5,k)
       enddo
    enddo

    j=4
    do i=nxmnghp1,nx
       l=i-nxmngh
       do k=1,nz
          dueta(l,j,k)= as4p3(1)*uf(l,1,k)+as4p3(2)*uf(l,2,k) &
                      + as4p3(3)*uf(l,3,k)+as4p3(4)*uf(l,4,k) &
                      + as4p3(5)*uf(l,5,k)+as4p3(6)*uf(l,6,k)
          dveta(l,j,k)= as4p3(1)*vf(l,1,k)+as4p3(2)*vf(l,2,k) &
                      + as4p3(3)*vf(l,3,k)+as4p3(4)*vf(l,4,k) &
                      + as4p3(5)*vf(l,5,k)+as4p3(6)*vf(l,6,k)
          dweta(l,j,k)= as4p3(1)*wf(l,1,k)+as4p3(2)*wf(l,2,k) &
                      + as4p3(3)*wf(l,3,k)+as4p3(4)*wf(l,4,k) &
                      + as4p3(5)*wf(l,5,k)+as4p3(6)*wf(l,6,k)
          dpeta(l,j,k)= as4p3(1)*pf(l,1,k)+as4p3(2)*pf(l,2,k) &
                      + as4p3(3)*pf(l,3,k)+as4p3(4)*pf(l,4,k) &
                      + as4p3(5)*pf(l,5,k)+as4p3(6)*pf(l,6,k)
          dreta(l,j,k)= as4p3(1)*rf(l,1,k)+as4p3(2)*rf(l,2,k) &
                      + as4p3(3)*rf(l,3,k)+as4p3(4)*rf(l,4,k) &
                      + as4p3(5)*rf(l,5,k)+as4p3(6)*rf(l,6,k)
       enddo
    enddo
    
    do j=5,ngh
       do i=nxmnghp1,nx
          l=i-nxmngh
          do k=1,nz
             dueta(l,j,k)= a5(1)*(uf(l,j+1,k)-uf(l,j-1,k)) &
                         + a5(2)*(uf(l,j+2,k)-uf(l,j-2,k))
             dveta(l,j,k)= a5(1)*(vf(l,j+1,k)-vf(l,j-1,k)) &
                         + a5(2)*(vf(l,j+2,k)-vf(l,j-2,k))
             dweta(l,j,k)= a5(1)*(wf(l,j+1,k)-wf(l,j-1,k)) &
                         + a5(2)*(wf(l,j+2,k)-wf(l,j-2,k))
             dpeta(l,j,k)= a5(1)*(pf(l,j+1,k)-pf(l,j-1,k)) &
                         + a5(2)*(pf(l,j+2,k)-pf(l,j-2,k))
             dreta(l,j,k)= a5(1)*(rf(l,j+1,k)-rf(l,j-1,k)) &
                         + a5(2)*(rf(l,j+2,k)-rf(l,j-2,k))
          enddo
       enddo
    enddo

    ! Compute rt,ut,vt,wt & pt
    ! ========================
    do i=nxmnghp1,nx
       l=i-nxmngh
       do j=1,ngh
          do k=1,nz
             pt(l,j,k) = vg(l,j,k)*( ((dpksi(l,j,k)*y_eta(i,j)-dpeta(l,j,k)*y_ksi(i,j))*cosphi(l,j) &
                                     +(dpeta(l,j,k)*x_ksi(i,j)-dpksi(l,j,k)*x_eta(i,j))*sinphi(l,j))*ijacob(i,j) &
                                     + pf(l,j,k)*ir(l,j) )
             ut(l,j,k) = vg(l,j,k)*( ((duksi(l,j,k)*y_eta(i,j)-dueta(l,j,k)*y_ksi(i,j))*cosphi(l,j) &
                                     +(dueta(l,j,k)*x_ksi(i,j)-duksi(l,j,k)*x_eta(i,j))*sinphi(l,j))*ijacob(i,j) &
                                     + uf(l,j,k)*ir(l,j) )
             vt(l,j,k) = vg(l,j,k)*( ((dvksi(l,j,k)*y_eta(i,j)-dveta(l,j,k)*y_ksi(i,j))*cosphi(l,j) &
                                     +(dveta(l,j,k)*x_ksi(i,j)-dvksi(l,j,k)*x_eta(i,j))*sinphi(l,j))*ijacob(i,j) &
                                     + vf(l,j,k)*ir(l,j) )
             wt(l,j,k) = vg(l,j,k)*( ((dwksi(l,j,k)*y_eta(i,j)-dweta(l,j,k)*y_ksi(i,j))*cosphi(l,j) &
                                     +(dweta(l,j,k)*x_ksi(i,j)-dwksi(l,j,k)*x_eta(i,j))*sinphi(l,j))*ijacob(i,j) &
                                     + wf(l,j,k)*ir(l,j) )
             rt(l,j,k) = vg(l,j,k)*( ((drksi(l,j,k)*y_eta(i,j)-dreta(l,j,k)*y_ksi(i,j))*cosphi(l,j) &
                                     +(dreta(l,j,k)*x_ksi(i,j)-drksi(l,j,k)*x_eta(i,j))*sinphi(l,j))*ijacob(i,j) &
                                     + rf(l,j,k)*ir(l,j) )
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

  end subroutine bc_TD2d_imax_jmin_c_SBP4

  !===============================================================================
  module subroutine bc_TD2d_imax_jmax_c_SBP4
  !===============================================================================
    !> 2D Tam & Dong's BC: boundary condition at imax-jmax (edge 1,2,2 /right-top) - curvilinear version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l,m
    real(wp) :: cp,av,c2_
    real(wp), dimension(-2:ngh,-2:ngh,nz) :: rf,uf,vf,wf,pf
    real(wp), dimension(ngh,ngh,nz)   :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(ngh,ngh,nz)   :: duksi,dvksi,dwksi,dpksi,drksi
    real(wp), dimension(ngh,ngh,nz)   :: dueta,dveta,dweta,dpeta,dreta
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

    ! Non-centered derivatives in x-direction
    ! =======================================
    ! (Tam & Webb DRP schemes)
    do i=nxmnghp1,nx-4
       l=i-nxmngh
       do j=nymnghp1,ny
          m=j-nymngh
          do k=1,nz
             duksi(l,m,k) = a5(1)* (uf(l+1,m,k)-uf(l-1,m,k)) &
                          + a5(2)* (uf(l+2,m,k)-uf(l-2,m,k))
             dvksi(l,m,k) = a5(1)* (vf(l+1,m,k)-vf(l-1,m,k)) &
                          + a5(2)* (vf(l+2,m,k)-vf(l-2,m,k))
             dwksi(l,m,k) = a5(1)* (wf(l+1,m,k)-wf(l-1,m,k)) &
                          + a5(2)* (wf(l+2,m,k)-wf(l-2,m,k))
             dpksi(l,m,k) = a5(1)* (pf(l+1,m,k)-pf(l-1,m,k)) &
                          + a5(2)* (pf(l+2,m,k)-pf(l-2,m,k))
             drksi(l,m,k) = a5(1)* (rf(l+1,m,k)-rf(l-1,m,k)) &
                          + a5(2)* (rf(l+2,m,k)-rf(l-2,m,k))
          enddo
       enddo
    enddo
    
    i=nx-3
    l=i-nxmngh
    do j=nymnghp1,ny
       m=j-nymngh
       do k=1,nz
          duksi(l,m,k) = as4m3(1)*uf(l+3,m,k)+as4m3(2)*uf(l+2,m,k) &
                       + as4m3(3)*uf(l+1,m,k)+as4m3(4)*uf(l  ,m,k) &
                       + as4m3(5)*uf(l-1,m,k)+as4m3(6)*uf(l-2,m,k)
          dvksi(l,m,k) = as4m3(1)*vf(l+3,m,k)+as4m3(2)*vf(l+2,m,k) &
                       + as4m3(3)*vf(l+1,m,k)+as4m3(4)*vf(l  ,m,k) &
                       + as4m3(5)*vf(l-1,m,k)+as4m3(6)*vf(l-2,m,k)
          dwksi(l,m,k) = as4m3(1)*wf(l+3,m,k)+as4m3(2)*wf(l+2,m,k) &
                       + as4m3(3)*wf(l+1,m,k)+as4m3(4)*wf(l  ,m,k) &
                       + as4m3(5)*wf(l-1,m,k)+as4m3(6)*wf(l-2,m,k)
          dpksi(l,m,k) = as4m3(1)*pf(l+3,m,k)+as4m3(2)*pf(l+2,m,k) &
                       + as4m3(3)*pf(l+1,m,k)+as4m3(4)*pf(l  ,m,k) &
                       + as4m3(5)*pf(l-1,m,k)+as4m3(6)*pf(l-2,m,k)
          drksi(l,m,k) = as4m3(1)*rf(l+3,m,k)+as4m3(2)*rf(l+2,m,k) &
                       + as4m3(3)*rf(l+1,m,k)+as4m3(4)*rf(l  ,m,k) &
                       + as4m3(5)*rf(l-1,m,k)+as4m3(6)*rf(l-2,m,k)
       enddo
    enddo
    
    i=nx-2
    l=i-nxmngh
    do j=nymnghp1,ny
       m=j-nymngh
       do k=1,nz
          duksi(l,m,k) = as4m2(1)*uf(l+2,m,k)+as4m2(2)*uf(l+1,m,k) &
                       + as4m2(3)*uf(l  ,m,k)+as4m2(4)*uf(l-1,m,k) &
                       + as4m2(5)*uf(l-2,m,k)
          dvksi(l,m,k) = as4m2(1)*vf(l+2,m,k)+as4m2(2)*vf(l+1,m,k) &
                       + as4m2(3)*vf(l  ,m,k)+as4m2(4)*vf(l-1,m,k) &
                       + as4m2(5)*vf(l-2,m,k)
          dwksi(l,m,k) = as4m2(1)*wf(l+2,m,k)+as4m2(2)*wf(l+1,m,k) &
                       + as4m2(3)*wf(l  ,m,k)+as4m2(4)*wf(l-1,m,k) &
                       + as4m2(5)*wf(l-2,m,k)
          dpksi(l,m,k) = as4m2(1)*pf(l+2,m,k)+as4m2(2)*pf(l+1,m,k) &
                       + as4m2(3)*pf(l  ,m,k)+as4m2(4)*pf(l-1,m,k) &
                       + as4m2(5)*pf(l-2,m,k)
          drksi(l,m,k) = as4m2(1)*rf(l+2,m,k)+as4m2(2)*rf(l+1,m,k) &
                       + as4m2(3)*rf(l  ,m,k)+as4m2(4)*rf(l-1,m,k) &
                       + as4m2(5)*rf(l-2,m,k)
       enddo
    enddo
    
    i=nx-1
    l=i-nxmngh
    do j=nymnghp1,ny
       m=j-nymngh
       do k=1,nz
          duksi(l,m,k) = as4m1(1)*uf(l+1,m,k)+as4m1(2)*uf(l  ,m,k) &
                       + as4m1(3)*uf(l-1,m,k)+as4m1(4)*uf(l-2,m,k)
          dvksi(l,m,k) = as4m1(1)*vf(l+1,m,k)+as4m1(2)*vf(l  ,m,k) &
                       + as4m1(3)*vf(l-1,m,k)+as4m1(4)*vf(l-2,m,k)
          dwksi(l,m,k) = as4m1(1)*wf(l+1,m,k)+as4m1(2)*wf(l  ,m,k) &
                       + as4m1(3)*wf(l-1,m,k)+as4m1(4)*wf(l-2,m,k)
          dpksi(l,m,k) = as4m1(1)*pf(l+1,m,k)+as4m1(2)*pf(l  ,m,k) &
                       + as4m1(3)*pf(l-1,m,k)+as4m1(4)*pf(l-2,m,k)
          drksi(l,m,k) = as4m1(1)*rf(l+1,m,k)+as4m1(2)*rf(l  ,m,k) &
                       + as4m1(3)*rf(l-1,m,k)+as4m1(4)*rf(l-2,m,k)
       enddo
    enddo
    
    i=nx
    l=i-nxmngh
    do j=nymnghp1,ny
       m=j-nymngh
       do k=1,nz
          duksi(l,m,k) = as4m0(1)*uf(l  ,m,k)+as4m0(2)*uf(l-1,m,k) &
                       + as4m0(3)*uf(l-2,m,k)+as4m0(4)*uf(l-3,m,k)
          dvksi(l,m,k) = as4m0(1)*vf(l  ,m,k)+as4m0(2)*vf(l-1,m,k) &
                       + as4m0(3)*vf(l-2,m,k)+as4m0(4)*vf(l-3,m,k)
          dwksi(l,m,k) = as4m0(1)*wf(l  ,m,k)+as4m0(2)*wf(l-1,m,k) &
                       + as4m0(3)*wf(l-2,m,k)+as4m0(4)*wf(l-3,m,k)
          dpksi(l,m,k) = as4m0(1)*pf(l  ,m,k)+as4m0(2)*pf(l-1,m,k) &
                       + as4m0(3)*pf(l-2,m,k)+as4m0(4)*pf(l-3,m,k)
          drksi(l,m,k) = as4m0(1)*rf(l  ,m,k)+as4m0(2)*rf(l-1,m,k) &
                       + as4m0(3)*rf(l-2,m,k)+as4m0(4)*rf(l-3,m,k)
       enddo
    enddo
 
    ! Non-centered derivatives in y-direction
    ! =======================================
    ! (Tam & Webb DRP schemes)
    !j=nymnghp1,ny-4
    do m=ngh-4,ngh-4
       do l=1,ngh
          do k=1,nz
             dueta(l,m,k) = a5(1)*(uf(l,m+1,k)-uf(l,m-1,k)) &
                          + a5(2)*(uf(l,m+2,k)-uf(l,m-2,k))
             dveta(l,m,k) = a5(1)*(vf(l,m+1,k)-vf(l,m-1,k)) &
                          + a5(2)*(vf(l,m+2,k)-vf(l,m-2,k))
             dweta(l,m,k) = a5(1)*(wf(l,m+1,k)-wf(l,m-1,k)) &
                          + a5(2)*(wf(l,m+2,k)-wf(l,m-2,k))
             dpeta(l,m,k) = a5(1)*(pf(l,m+1,k)-pf(l,m-1,k)) &
                          + a5(2)*(pf(l,m+2,k)-pf(l,m-2,k))
             dreta(l,m,k) = a5(1)*(rf(l,m+1,k)-rf(l,m-1,k)) &
                          + a5(2)*(rf(l,m+2,k)-rf(l,m-2,k))
          enddo
       enddo
    enddo

    !j=ny-3
    m=ngh-3
    do l=1,ngh
       do k=1,nz
          dueta(l,m,k) = as4m3(1)*uf(l,m+3,k)+as4m3(2)*uf(l,m+2,k) &
                       + as4m3(3)*uf(l,m+1,k)+as4m3(4)*uf(l,m  ,k) &
                       + as4m3(5)*uf(l,m-1,k)+as4m3(6)*uf(l,m-2,k)
          dveta(l,m,k) = as4m3(1)*vf(l,m+3,k)+as4m3(2)*vf(l,m+2,k) &
                       + as4m3(3)*vf(l,m+1,k)+as4m3(4)*vf(l,m  ,k) &
                       + as4m3(5)*vf(l,m-1,k)+as4m3(6)*vf(l,m-2,k)
          dweta(l,m,k) = as4m3(1)*wf(l,m+3,k)+as4m3(2)*wf(l,m+2,k) &
                       + as4m3(3)*wf(l,m+1,k)+as4m3(4)*wf(l,m  ,k) &
                       + as4m3(5)*wf(l,m-1,k)+as4m3(6)*wf(l,m-2,k)
          dpeta(l,m,k) = as4m3(1)*pf(l,m+3,k)+as4m3(2)*pf(l,m+2,k) &
                       + as4m3(3)*pf(l,m+1,k)+as4m3(4)*pf(l,m  ,k) &
                       + as4m3(5)*pf(l,m-1,k)+as4m3(6)*pf(l,m-2,k)
          dreta(l,m,k) = as4m3(1)*rf(l,m+3,k)+as4m3(2)*rf(l,m+2,k) &
                       + as4m3(3)*rf(l,m+1,k)+as4m3(4)*rf(l,m  ,k) &
                       + as4m3(5)*rf(l,m-1,k)+as4m3(6)*rf(l,m-2,k)
       enddo
    enddo
    
    !j=ny-2
    m=ngh-2
    do l=1,ngh
       do k=1,nz
          dueta(l,m,k) = as4m2(1)*uf(l,m+2,k)+as4m2(2)*uf(l,m+1,k) &
                       + as4m2(3)*uf(l,m  ,k)+as4m2(4)*uf(l,m-1,k) &
                       + as4m2(5)*uf(l,m-2,k)
          dveta(l,m,k) = as4m2(1)*vf(l,m+2,k)+as4m2(2)*vf(l,m+1,k) &
                       + as4m2(3)*vf(l,m  ,k)+as4m2(4)*vf(l,m-1,k) &
                       + as4m2(5)*vf(l,m-2,k)
          dweta(l,m,k) = as4m2(1)*wf(l,m+2,k)+as4m2(2)*wf(l,m+1,k) &
                       + as4m2(3)*wf(l,m  ,k)+as4m2(4)*wf(l,m-1,k) &
                       + as4m2(5)*wf(l,m-2,k)
          dpeta(l,m,k) = as4m2(1)*pf(l,m+2,k)+as4m2(2)*pf(l,m+1,k) &
                       + as4m2(3)*pf(l,m  ,k)+as4m2(4)*pf(l,m-1,k) &
                       + as4m2(5)*pf(l,m-2,k)
          dreta(l,m,k) = as4m2(1)*rf(l,m+2,k)+as4m2(2)*rf(l,m+1,k) &
                       + as4m2(3)*rf(l,m  ,k)+as4m2(4)*rf(l,m-1,k) &
                       + as4m2(5)*rf(l,m-2,k)
       enddo
    enddo
    
    !j=ny-1
    m=ngh-1
    do l=1,ngh
       do k=1,nz
          dueta(l,m,k) = as4m1(1)*uf(l,m+1,k)+as4m1(2)*uf(l,m  ,k) &
                       + as4m1(3)*uf(l,m-1,k)+as4m1(4)*uf(l,m-2,k)
          dveta(l,m,k) = as4m1(1)*vf(l,m+1,k)+as4m1(2)*vf(l,m  ,k) &
                       + as4m1(3)*vf(l,m-1,k)+as4m1(4)*vf(l,m-2,k)
          dweta(l,m,k) = as4m1(1)*wf(l,m+1,k)+as4m1(2)*wf(l,m  ,k) &
                       + as4m1(3)*wf(l,m-1,k)+as4m1(4)*wf(l,m-2,k)
          dpeta(l,m,k) = as4m1(1)*pf(l,m+1,k)+as4m1(2)*pf(l,m  ,k) &
                       + as4m1(3)*pf(l,m-1,k)+as4m1(4)*pf(l,m-2,k)
          dreta(l,m,k) = as4m1(1)*rf(l,m+1,k)+as4m1(2)*rf(l,m  ,k) &
                       + as4m1(3)*rf(l,m-1,k)+as4m1(4)*rf(l,m-2,k)
       enddo
    enddo
    
    !j=ny
    m=ngh
    do l=1,ngh
       do k=1,nz
          dueta(l,m,k) = as4m0(1)*uf(l,m  ,k)+as4m0(2)*uf(l,m-1,k) &
                       + as4m0(3)*uf(l,m-2,k)+as4m0(4)*uf(l,m-3,k)
          dveta(l,m,k) = as4m0(1)*vf(l,m  ,k)+as4m0(2)*vf(l,m-1,k) &
                       + as4m0(3)*vf(l,m-2,k)+as4m0(4)*vf(l,m-3,k)
          dweta(l,m,k) = as4m0(1)*wf(l,m  ,k)+as4m0(2)*wf(l,m-1,k) &
                       + as4m0(3)*wf(l,m-2,k)+as4m0(4)*wf(l,m-3,k)
          dpeta(l,m,k) = as4m0(1)*pf(l,m  ,k)+as4m0(2)*pf(l,m-1,k) &
                       + as4m0(3)*pf(l,m-2,k)+as4m0(4)*pf(l,m-3,k)
          dreta(l,m,k) = as4m0(1)*rf(l,m  ,k)+as4m0(2)*rf(l,m-1,k) &
                       + as4m0(3)*rf(l,m-2,k)+as4m0(4)*rf(l,m-3,k)
       enddo
    enddo
              
    ! Compute rt,ut,vt,wt & pt
    ! ========================
    do i=nxmnghp1,nx
       l=i-nxmngh
       do j=nymnghp1,ny
          m=j-nymngh
          do k=1,nz
             pt(l,m,k) = vg(l,m,k)*( ((dpksi(l,m,k)*y_eta(i,j)-dpeta(l,m,k)*y_ksi(i,j))*cosphi(l,m) &
                                     +(dpeta(l,m,k)*x_ksi(i,j)-dpksi(l,m,k)*x_eta(i,j))*sinphi(l,m))*ijacob(i,j) &
                                     + pf(l,m,k)*ir(l,m) )
             ut(l,m,k) = vg(l,m,k)*( ((duksi(l,m,k)*y_eta(i,j)-dueta(l,m,k)*y_ksi(i,j))*cosphi(l,m) &
                                     +(dueta(l,m,k)*x_ksi(i,j)-duksi(l,m,k)*x_eta(i,j))*sinphi(l,m))*ijacob(i,j) &
                                     + uf(l,m,k)*ir(l,m) )
             vt(l,m,k) = vg(l,m,k)*( ((dvksi(l,m,k)*y_eta(i,j)-dveta(l,m,k)*y_ksi(i,j))*cosphi(l,m) &
                                     +(dveta(l,m,k)*x_ksi(i,j)-dvksi(l,m,k)*x_eta(i,j))*sinphi(l,m))*ijacob(i,j) &
                                     + vf(l,m,k)*ir(l,m) )
             wt(l,m,k) = vg(l,m,k)*( ((dwksi(l,m,k)*y_eta(i,j)-dweta(l,m,k)*y_ksi(i,j))*cosphi(l,m) &
                                     +(dweta(l,m,k)*x_ksi(i,j)-dwksi(l,m,k)*x_eta(i,j))*sinphi(l,m))*ijacob(i,j) &
                                     + wf(l,m,k)*ir(l,m) )
             rt(l,m,k) = vg(l,m,k)*( ((drksi(l,m,k)*y_eta(i,j)-dreta(l,m,k)*y_ksi(i,j))*cosphi(l,m) &
                                     +(dreta(l,m,k)*x_ksi(i,j)-drksi(l,m,k)*x_eta(i,j))*sinphi(l,m))*ijacob(i,j) &
                                     + rf(l,m,k)*ir(l,m) )
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

  end subroutine bc_TD2d_imax_jmax_c_SBP4

end submodule smod_TamDong2d_SBP4_edges_c
