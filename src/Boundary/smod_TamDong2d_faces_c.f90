!===============================================================
submodule (mod_TamDong2d_c) smod_TamDong2d_faces_c
!===============================================================
  !> author: XG
  !> date: February 2020 - modif January 2022
  !> Non-reflecting boundary conditions of Tam and Dong
  !> - 2D (periodic) curvilinear version - routines for faces
!===============================================================

contains

  !===============================================================================
  module subroutine bc_TD2d_imin_c
  !===============================================================================
    !> 2D Tam & Dong's BC: boundary condition at imin (left) - curvilinear version -
  !===============================================================================
    use mod_RFM
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: cp,av,inn1,ntm1
    real(wp) :: duksi,dvksi,dwksi,dpksi,drksi
    real(wp) :: dueta,dveta,dweta,dpeta,dreta
    real(wp), dimension(1:ngh+3,ny1:ny2,nz) :: rf,uf,vf,wf,pf
    real(wp), dimension(ngh,ny,nz) :: rt,ut,vt,wt,pt,vg,c2_
    !-------------------------------------------------------------------------
    ! eigenmode or RFM disturbances
    ! real(wp) :: pt_in_,ut_in_,vt_in_,wt_in_,rt_in_
    ! real(wp) :: dp_in,du_in,dv_in,dw_in,dr_in
    real(wp), dimension(ngh,ny,nz) :: ut_in,vt_in,wt_in
    !-------------------------------------------------------------------------
    ! added for is_mean_ref
    integer :: n_moy
    real(wp) :: eps
    !-------------------------------------------------------------------------

    ! Pointers for polar coordinates
    ! ==============================
    ir=>BC_face(1,1)%ir
    cosphi=>BC_face(1,1)%cosphi
    sinphi=>BC_face(1,1)%sinphi
    
    ! Sound speed
    ! ===========
    do k=1,nz
       do j=1,ny
          do i=1,ngh
             c2_(i,j,k)=c2calc_tro(Tmp(i,j,k),rho_n(i,j,k))
          enddo
       enddo
    enddo

    ! Compute online time-averaged primitive variables
    ! ================================================
    if (irk==nrk) then
       inn1 = 1.0_wp/dble(ntotal)
       ntm1=dble(ntotal-1)
       if (BC_face(1,1)%is_mean_ref) then
          ! weight mean field with imposed reference field in inlet planes
          ! ~> weak imposition of fixed mean quantities
          eps=0.1_wp
          n_moy=100
          inn1=inn1*eps
          ! time-averaged primitive variables
          do k=1,nz
             do j=1,ny
                BC_face(1,1)%U0(1:nghp3,j,k,1)=(ntm1*BC_face(1,1)%U0(1:nghp3,j,k,1)+rho_n(1:nghp3,j,k))*inn1 &
                                                  + (1.0_wp-eps)*BC_face(1,1)%Uref(j,k,1)
                BC_face(1,1)%U0(1:nghp3,j,k,2)=(ntm1*BC_face(1,1)%U0(1:nghp3,j,k,2)+uu(1:nghp3,j,k))*inn1    &
                                               + (1.0_wp-eps)*BC_face(1,1)%Uref(j,k,2)
                BC_face(1,1)%U0(1:nghp3,j,k,3)=(ntm1*BC_face(1,1)%U0(1:nghp3,j,k,3)+vv(1:nghp3,j,k))*inn1    &
                                               + (1.0_wp-eps)*BC_face(1,1)%Uref(j,k,3)
                BC_face(1,1)%U0(1:nghp3,j,k,4)=(ntm1*BC_face(1,1)%U0(1:nghp3,j,k,4)+ww(1:nghp3,j,k))*inn1    &
                                               + (1.0_wp-eps)*BC_face(1,1)%Uref(j,k,4)
                BC_face(1,1)%U0(1:nghp3,j,k,5)=(ntm1*BC_face(1,1)%U0(1:nghp3,j,k,5)+prs(1:nghp3,j,k))*inn1   &
                                               + (1.0_wp-eps)*BC_face(1,1)%Uref(j,k,5)
             enddo
          enddo
          ! time-averaged sound speed squared
          do k=1,nz
             do j=1,ny
                BC_face(1,1)%U0(1:ngh,j,k,6)=(ntm1*BC_face(1,1)%U0(1:ngh,j,k,6)+c2_(1:ngh,j,k))*inn1 &
                                                  + (1.0_wp-eps)*BC_face(1,1)%Uref(j,k,6)
             enddo
          enddo
       else
          ! time-averaged primitive variables
          do k=1,nz
             BC_face(1,1)%U0(1:nghp3,:,k,1)=(ntm1*BC_face(1,1)%U0(1:nghp3,:,k,1)+rho_n(1:nghp3,:,k))*inn1
             BC_face(1,1)%U0(1:nghp3,:,k,2)=(ntm1*BC_face(1,1)%U0(1:nghp3,:,k,2)+uu(1:nghp3,:,k))*inn1
             BC_face(1,1)%U0(1:nghp3,:,k,3)=(ntm1*BC_face(1,1)%U0(1:nghp3,:,k,3)+vv(1:nghp3,:,k))*inn1
             BC_face(1,1)%U0(1:nghp3,:,k,4)=(ntm1*BC_face(1,1)%U0(1:nghp3,:,k,4)+ww(1:nghp3,:,k))*inn1
             BC_face(1,1)%U0(1:nghp3,:,k,5)=(ntm1*BC_face(1,1)%U0(1:nghp3,:,k,5)+prs(1:nghp3,:,k))*inn1
          enddo
          ! time-averaged sound speed squared
          do k=1,nz
             BC_face(1,1)%U0(1:ngh,1:ny,k,6)=(ntm1*BC_face(1,1)%U0(1:ngh,1:ny,k,6)+c2_(1:ngh,1:ny,k))*inn1
          enddo
       endif
    endif
    
    ! Compute fluctuations
    ! ====================
    do i=1,nghp3
       do j=ny1,ny2
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
       do j=ndy,nfy
          do k=1,nz
             ! vg(i,j,k)= BC_face(1,1)%U0(i,j,k,2)*cosphi(i,j)+BC_face(1,1)%U0(i,j,k,3)*sinphi(i,j) &
             !      + sqrt(BC_face(1,1)%U0(i,j,k,6)-BC_face(1,1)%U0(i,j,k,4)**2 &
             !      -(BC_face(1,1)%U0(i,j,k,2)*sinphi(i,j)-BC_face(1,1)%U0(i,j,k,3)*cosphi(i,j))**2)
             ! Protection of value with abs() if supersonic at outlet
             vg(i,j,k)= BC_face(1,1)%U0(i,j,k,2)*cosphi(i,j)+BC_face(1,1)%U0(i,j,k,3)*sinphi(i,j) &
                  + sqrt(abs(BC_face(1,1)%U0(i,j,k,6)-BC_face(1,1)%U0(i,j,k,4)**2 &
                  -(BC_face(1,1)%U0(i,j,k,2)*sinphi(i,j)-BC_face(1,1)%U0(i,j,k,3)*cosphi(i,j))**2))
          enddo
       enddo
    enddo

    ! Compute vg*[cos(phi)*dq/dx+sin(phi)*dq/dy+q/r]
    ! ==============================================
    ! (Tam & Webb DRP schemes)
    i=1
    do j=ndy,nfy
       do k=1,nz
          duksi = a06(1)*uf(1,j,k)+a06(2)*uf(2,j,k) &
                + a06(3)*uf(3,j,k)+a06(4)*uf(4,j,k) &
                + a06(5)*uf(5,j,k)+a06(6)*uf(6,j,k) &
                + a06(7)*uf(7,j,k)
          dvksi = a06(1)*vf(1,j,k)+a06(2)*vf(2,j,k) &
                + a06(3)*vf(3,j,k)+a06(4)*vf(4,j,k) &
                + a06(5)*vf(5,j,k)+a06(6)*vf(6,j,k) &
                + a06(7)*vf(7,j,k)
          dwksi = a06(1)*wf(1,j,k)+a06(2)*wf(2,j,k) &
                + a06(3)*wf(3,j,k)+a06(4)*wf(4,j,k) &
                + a06(5)*wf(5,j,k)+a06(6)*wf(6,j,k) &
                + a06(7)*wf(7,j,k)
          dpksi = a06(1)*pf(1,j,k)+a06(2)*pf(2,j,k) &
                + a06(3)*pf(3,j,k)+a06(4)*pf(4,j,k) &
                + a06(5)*pf(5,j,k)+a06(6)*pf(6,j,k) &
                + a06(7)*pf(7,j,k)
          drksi = a06(1)*rf(1,j,k)+a06(2)*rf(2,j,k) &
                + a06(3)*rf(3,j,k)+a06(4)*rf(4,j,k) &
                + a06(5)*rf(5,j,k)+a06(6)*rf(6,j,k) &
                + a06(7)*rf(7,j,k)

          dueta = a7(1)*( uf(i,j+1,k) - uf(i,j-1,k) ) + &
                  a7(2)*( uf(i,j+2,k) - uf(i,j-2,k) ) + &
                  a7(3)*( uf(i,j+3,k) - uf(i,j-3,k) )
          dveta = a7(1)*( vf(i,j+1,k) - vf(i,j-1,k) ) + &
                  a7(2)*( vf(i,j+2,k) - vf(i,j-2,k) ) + &
                  a7(3)*( vf(i,j+3,k) - vf(i,j-3,k) )
          dweta = a7(1)*( wf(i,j+1,k) - wf(i,j-1,k) ) + &
                  a7(2)*( wf(i,j+2,k) - wf(i,j-2,k) ) + &
                  a7(3)*( wf(i,j+3,k) - wf(i,j-3,k) )
          dpeta = a7(1)*( pf(i,j+1,k) - pf(i,j-1,k) ) + &
                  a7(2)*( pf(i,j+2,k) - pf(i,j-2,k) ) + &
                  a7(3)*( pf(i,j+3,k) - pf(i,j-3,k) )
          dreta = a7(1)*( rf(i,j+1,k) - rf(i,j-1,k) ) + &
                  a7(2)*( rf(i,j+2,k) - rf(i,j-2,k) ) + &
                  a7(3)*( rf(i,j+3,k) - rf(i,j-3,k) )
          
          pt(i,j,k) = vg(i,j,k)*( ((dpksi*y_eta(i,j)-dpeta*y_ksi(i,j))*cosphi(i,j) &
                                  +(dpeta*x_ksi(i,j)-dpksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                  + pf(i,j,k)*ir(i,j) )
          ut(i,j,k) = vg(i,j,k)*( ((duksi*y_eta(i,j)-dueta*y_ksi(i,j))*cosphi(i,j) &
                                  +(dueta*x_ksi(i,j)-duksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                  + uf(i,j,k)*ir(i,j) )
          vt(i,j,k) = vg(i,j,k)*( ((dvksi*y_eta(i,j)-dveta*y_ksi(i,j))*cosphi(i,j) &
                                  +(dveta*x_ksi(i,j)-dvksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                  + vf(i,j,k)*ir(i,j) )
          wt(i,j,k) = vg(i,j,k)*( ((dwksi*y_eta(i,j)-dweta*y_ksi(i,j))*cosphi(i,j) &
                                  +(dweta*x_ksi(i,j)-dwksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                  + wf(i,j,k)*ir(i,j) )
          rt(i,j,k) = vg(i,j,k)*( ((drksi*y_eta(i,j)-dreta*y_ksi(i,j))*cosphi(i,j) &
                                  +(dreta*x_ksi(i,j)-drksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                  + rf(i,j,k)*ir(i,j) )
       enddo
    enddo

    i=2
    do j=ndy,nfy
       do k=1,nz
          duksi = a15(1)*uf(1,j,k)+a15(2)*uf(2,j,k) &
                + a15(3)*uf(3,j,k)+a15(4)*uf(4,j,k) &
                + a15(5)*uf(5,j,k)+a15(6)*uf(6,j,k) &
                + a15(7)*uf(7,j,k)
          dvksi = a15(1)*vf(1,j,k)+a15(2)*vf(2,j,k) &
                + a15(3)*vf(3,j,k)+a15(4)*vf(4,j,k) &
                + a15(5)*vf(5,j,k)+a15(6)*vf(6,j,k) &
                + a15(7)*vf(7,j,k)
          dwksi = a15(1)*wf(1,j,k)+a15(2)*wf(2,j,k) &
                + a15(3)*wf(3,j,k)+a15(4)*wf(4,j,k) &
                + a15(5)*wf(5,j,k)+a15(6)*wf(6,j,k) &
                + a15(7)*wf(7,j,k)
          dpksi = a15(1)*pf(1,j,k)+a15(2)*pf(2,j,k) &
                + a15(3)*pf(3,j,k)+a15(4)*pf(4,j,k) &
                + a15(5)*pf(5,j,k)+a15(6)*pf(6,j,k) &
                + a15(7)*pf(7,j,k)
          drksi = a15(1)*rf(1,j,k)+a15(2)*rf(2,j,k) &
                + a15(3)*rf(3,j,k)+a15(4)*rf(4,j,k) &
                + a15(5)*rf(5,j,k)+a15(6)*rf(6,j,k) &
                + a15(7)*rf(7,j,k)

          dueta = a7(1)*( uf(i,j+1,k) - uf(i,j-1,k) ) + &
                  a7(2)*( uf(i,j+2,k) - uf(i,j-2,k) ) + &
                  a7(3)*( uf(i,j+3,k) - uf(i,j-3,k) )
          dveta = a7(1)*( vf(i,j+1,k) - vf(i,j-1,k) ) + &
                  a7(2)*( vf(i,j+2,k) - vf(i,j-2,k) ) + &
                  a7(3)*( vf(i,j+3,k) - vf(i,j-3,k) )
          dweta = a7(1)*( wf(i,j+1,k) - wf(i,j-1,k) ) + &
                  a7(2)*( wf(i,j+2,k) - wf(i,j-2,k) ) + &
                  a7(3)*( wf(i,j+3,k) - wf(i,j-3,k) )
          dpeta = a7(1)*( pf(i,j+1,k) - pf(i,j-1,k) ) + &
                  a7(2)*( pf(i,j+2,k) - pf(i,j-2,k) ) + &
                  a7(3)*( pf(i,j+3,k) - pf(i,j-3,k) )
          dreta = a7(1)*( rf(i,j+1,k) - rf(i,j-1,k) ) + &
                  a7(2)*( rf(i,j+2,k) - rf(i,j-2,k) ) + &
                  a7(3)*( rf(i,j+3,k) - rf(i,j-3,k) )

          pt(i,j,k) = vg(i,j,k)*( ((dpksi*y_eta(i,j)-dpeta*y_ksi(i,j))*cosphi(i,j) &
                                  +(dpeta*x_ksi(i,j)-dpksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                  + pf(i,j,k)*ir(i,j) )
          ut(i,j,k) = vg(i,j,k)*( ((duksi*y_eta(i,j)-dueta*y_ksi(i,j))*cosphi(i,j) &
                                  +(dueta*x_ksi(i,j)-duksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                  + uf(i,j,k)*ir(i,j) )
          vt(i,j,k) = vg(i,j,k)*( ((dvksi*y_eta(i,j)-dveta*y_ksi(i,j))*cosphi(i,j) &
                                  +(dveta*x_ksi(i,j)-dvksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                  + vf(i,j,k)*ir(i,j) )
          wt(i,j,k) = vg(i,j,k)*( ((dwksi*y_eta(i,j)-dweta*y_ksi(i,j))*cosphi(i,j) &
                                  +(dweta*x_ksi(i,j)-dwksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                  + wf(i,j,k)*ir(i,j) )
          rt(i,j,k) = vg(i,j,k)*( ((drksi*y_eta(i,j)-dreta*y_ksi(i,j))*cosphi(i,j) &
                                  +(dreta*x_ksi(i,j)-drksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                  + rf(i,j,k)*ir(i,j) )
       enddo
    enddo

    i=3
    do j=ndy,nfy
       do k=1,nz
          duksi = a24(1)*uf(1,j,k)+a24(2)*uf(2,j,k) &
                + a24(3)*uf(3,j,k)+a24(4)*uf(4,j,k) &
                + a24(5)*uf(5,j,k)+a24(6)*uf(6,j,k) &
                + a24(7)*uf(7,j,k)
          dvksi = a24(1)*vf(1,j,k)+a24(2)*vf(2,j,k) &
                + a24(3)*vf(3,j,k)+a24(4)*vf(4,j,k) &
                + a24(5)*vf(5,j,k)+a24(6)*vf(6,j,k) &
                + a24(7)*vf(7,j,k)
          dwksi = a24(1)*wf(1,j,k)+a24(2)*wf(2,j,k) &
                + a24(3)*wf(3,j,k)+a24(4)*wf(4,j,k) &
                + a24(5)*wf(5,j,k)+a24(6)*wf(6,j,k) &
                + a24(7)*wf(7,j,k)
          dpksi = a24(1)*pf(1,j,k)+a24(2)*pf(2,j,k) &
                + a24(3)*pf(3,j,k)+a24(4)*pf(4,j,k) &
                + a24(5)*pf(5,j,k)+a24(6)*pf(6,j,k) &
                + a24(7)*pf(7,j,k)
          drksi = a24(1)*rf(1,j,k)+a24(2)*rf(2,j,k) &
                + a24(3)*rf(3,j,k)+a24(4)*rf(4,j,k) &
                + a24(5)*rf(5,j,k)+a24(6)*rf(6,j,k) &
                + a24(7)*rf(7,j,k)

          dueta = a7(1)*( uf(i,j+1,k) - uf(i,j-1,k) ) + &
                  a7(2)*( uf(i,j+2,k) - uf(i,j-2,k) ) + &
                  a7(3)*( uf(i,j+3,k) - uf(i,j-3,k) )
          dveta = a7(1)*( vf(i,j+1,k) - vf(i,j-1,k) ) + &
                  a7(2)*( vf(i,j+2,k) - vf(i,j-2,k) ) + &
                  a7(3)*( vf(i,j+3,k) - vf(i,j-3,k) )
          dweta = a7(1)*( wf(i,j+1,k) - wf(i,j-1,k) ) + &
                  a7(2)*( wf(i,j+2,k) - wf(i,j-2,k) ) + &
                  a7(3)*( wf(i,j+3,k) - wf(i,j-3,k) )
          dpeta = a7(1)*( pf(i,j+1,k) - pf(i,j-1,k) ) + &
                  a7(2)*( pf(i,j+2,k) - pf(i,j-2,k) ) + &
                  a7(3)*( pf(i,j+3,k) - pf(i,j-3,k) )
          dreta = a7(1)*( rf(i,j+1,k) - rf(i,j-1,k) ) + &
                  a7(2)*( rf(i,j+2,k) - rf(i,j-2,k) ) + &
                  a7(3)*( rf(i,j+3,k) - rf(i,j-3,k) )

          pt(i,j,k) = vg(i,j,k)*( ((dpksi*y_eta(i,j)-dpeta*y_ksi(i,j))*cosphi(i,j) &
                                  +(dpeta*x_ksi(i,j)-dpksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                  + pf(i,j,k)*ir(i,j) )
          ut(i,j,k) = vg(i,j,k)*( ((duksi*y_eta(i,j)-dueta*y_ksi(i,j))*cosphi(i,j) &
                                  +(dueta*x_ksi(i,j)-duksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                  + uf(i,j,k)*ir(i,j) )
          vt(i,j,k) = vg(i,j,k)*( ((dvksi*y_eta(i,j)-dveta*y_ksi(i,j))*cosphi(i,j) &
                                  +(dveta*x_ksi(i,j)-dvksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                  + vf(i,j,k)*ir(i,j) )
          wt(i,j,k) = vg(i,j,k)*( ((dwksi*y_eta(i,j)-dweta*y_ksi(i,j))*cosphi(i,j) &
                                  +(dweta*x_ksi(i,j)-dwksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                  + wf(i,j,k)*ir(i,j) )
          rt(i,j,k) = vg(i,j,k)*( ((drksi*y_eta(i,j)-dreta*y_ksi(i,j))*cosphi(i,j) &
                                  +(dreta*x_ksi(i,j)-drksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                  + rf(i,j,k)*ir(i,j) )
       enddo
    enddo

    do i=4,ngh
       do j=ndy,nfy
          do k=1,nz
             duksi = a7(1)* (uf(i+1,j,k)-uf(i-1,j,k)) &
                   + a7(2)* (uf(i+2,j,k)-uf(i-2,j,k)) &
                   + a7(3)* (uf(i+3,j,k)-uf(i-3,j,k))
             dvksi = a7(1)* (vf(i+1,j,k)-vf(i-1,j,k)) &
                   + a7(2)* (vf(i+2,j,k)-vf(i-2,j,k)) &
                   + a7(3)* (vf(i+3,j,k)-vf(i-3,j,k))
             dwksi = a7(1)* (wf(i+1,j,k)-wf(i-1,j,k)) &
                   + a7(2)* (wf(i+2,j,k)-wf(i-2,j,k)) &
                   + a7(3)* (wf(i+3,j,k)-wf(i-3,j,k))
             dpksi = a7(1)* (pf(i+1,j,k)-pf(i-1,j,k)) &
                   + a7(2)* (pf(i+2,j,k)-pf(i-2,j,k)) &
                   + a7(3)* (pf(i+3,j,k)-pf(i-3,j,k))
             drksi = a7(1)* (rf(i+1,j,k)-rf(i-1,j,k)) &
                   + a7(2)* (rf(i+2,j,k)-rf(i-2,j,k)) &
                   + a7(3)* (rf(i+3,j,k)-rf(i-3,j,k))

             dueta = a7(1)*( uf(i,j+1,k) - uf(i,j-1,k) ) + &
                     a7(2)*( uf(i,j+2,k) - uf(i,j-2,k) ) + &
                     a7(3)*( uf(i,j+3,k) - uf(i,j-3,k) )
             dveta = a7(1)*( vf(i,j+1,k) - vf(i,j-1,k) ) + &
                     a7(2)*( vf(i,j+2,k) - vf(i,j-2,k) ) + &
                     a7(3)*( vf(i,j+3,k) - vf(i,j-3,k) )
             dweta = a7(1)*( wf(i,j+1,k) - wf(i,j-1,k) ) + &
                     a7(2)*( wf(i,j+2,k) - wf(i,j-2,k) ) + &
                     a7(3)*( wf(i,j+3,k) - wf(i,j-3,k) )
             dpeta = a7(1)*( pf(i,j+1,k) - pf(i,j-1,k) ) + &
                     a7(2)*( pf(i,j+2,k) - pf(i,j-2,k) ) + &
                     a7(3)*( pf(i,j+3,k) - pf(i,j-3,k) )
             dreta = a7(1)*( rf(i,j+1,k) - rf(i,j-1,k) ) + &
                     a7(2)*( rf(i,j+2,k) - rf(i,j-2,k) ) + &
                     a7(3)*( rf(i,j+3,k) - rf(i,j-3,k) )

             pt(i,j,k) = vg(i,j,k)*( ((dpksi*y_eta(i,j)-dpeta*y_ksi(i,j))*cosphi(i,j) &
                                     +(dpeta*x_ksi(i,j)-dpksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                     + pf(i,j,k)*ir(i,j) )
             ut(i,j,k) = vg(i,j,k)*( ((duksi*y_eta(i,j)-dueta*y_ksi(i,j))*cosphi(i,j) &
                                     +(dueta*x_ksi(i,j)-duksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                     + uf(i,j,k)*ir(i,j) )
             vt(i,j,k) = vg(i,j,k)*( ((dvksi*y_eta(i,j)-dveta*y_ksi(i,j))*cosphi(i,j) &
                                     +(dveta*x_ksi(i,j)-dvksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                     + vf(i,j,k)*ir(i,j) )
             wt(i,j,k) = vg(i,j,k)*( ((dwksi*y_eta(i,j)-dweta*y_ksi(i,j))*cosphi(i,j) &
                                     +(dweta*x_ksi(i,j)-dwksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                     + wf(i,j,k)*ir(i,j) )
             rt(i,j,k) = vg(i,j,k)*( ((drksi*y_eta(i,j)-dreta*y_ksi(i,j))*cosphi(i,j) &
                                     +(dreta*x_ksi(i,j)-drksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                     + rf(i,j,k)*ir(i,j) )
          enddo
       enddo
    enddo

    ! Superimpose Random Fourier Modes (RFM) synthetic turbulence
    ! ===========================================================
    ! ut_in is the temporal derivative from RHS T&T applied to inlet dist.
    if (is_RFM) then
       call disturb_inlet_RFM_TamDong_imin(vg,ut_in,vt_in,wt_in)
       do i=1,ngh
          do j=ndy,nfy
             do k=1,nz
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
       do j=ndy,nfy
          do i=1,ngh
             cp=cpcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
             av=avcalc_tro(Tmp(i,j,k),rho_n(i,j,k))

             Krho(i,j,k)  = rt(i,j,k)
             Krhou(i,j,k) = uu(i,j,k)*rt(i,j,k)+rho_n(i,j,k)*ut(i,j,k)
             Krhov(i,j,k) = vv(i,j,k)*rt(i,j,k)+rho_n(i,j,k)*vt(i,j,k)
             Krhow(i,j,k) = ww(i,j,k)*rt(i,j,k)+rho_n(i,j,k)*wt(i,j,k)
             Krhoe(i,j,k) = cp/av*(pt(i,j,k)/c2_(i,j,k)-rt(i,j,k)) &
                  + (rhoe_n(i,j,k)+prs(i,j,k))/rho_n(i,j,k)*rt(i,j,k) &
                  + rho_n(i,j,k)*(uu(i,j,k)*ut(i,j,k)+vv(i,j,k)*vt(i,j,k)+ww(i,j,k)*wt(i,j,k))
          enddo
       enddo
    enddo

  end subroutine bc_TD2d_imin_c

  !===============================================================================
  module subroutine bc_TD2d_imax_c
  !===============================================================================
    !> 2D Tam & Dong's BC: boundary condition at imax (right) - curvilinear version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l
    real(wp) :: cp,av,inn1,ntm1
    real(wp) :: duksi,dvksi,dwksi,dpksi,drksi
    real(wp) :: dueta,dveta,dweta,dpeta,dreta
    real(wp), dimension(-2:ngh,ny1:ny2,nz) :: rf,uf,vf,wf,pf
    real(wp), dimension(ngh,ny,nz) :: rt,ut,vt,wt,pt,vg,c2_
    !-------------------------------------------------------------------------
    ! added for is_mean_ref
    integer :: n_moy
    real(wp) :: eps
    !-------------------------------------------------------------------------

    ! Pointers for polar coordinates
    ! ==============================
    ir=>BC_face(1,2)%ir
    cosphi=>BC_face(1,2)%cosphi
    sinphi=>BC_face(1,2)%sinphi

    ! Sound speed
    ! ===========
    do k=1,nz
       do j=1,ny
          do i=nxmnghp1,nx
             l=i-nxmngh
             c2_(l,j,k)=c2calc_tro(Tmp(i,j,k),rho_n(i,j,k))
          enddo
       enddo
    enddo

    ! Compute online time-averaged primitive variables
    ! ================================================
    if (irk==nrk) then
       inn1 = 1.0_wp/dble(ntotal)
       ntm1=dble(ntotal-1)
       if (BC_face(1,2)%is_mean_ref) then
          ! weight mean field with imposed reference field in inlet planes
          ! ~> weak imposition of fixed mean quantities
          eps=0.1_wp
          n_moy=100
          inn1=inn1*eps
          ! time-averaged primitive variables
          do k=1,nz
             do j=1,ny
                BC_face(1,2)%U0(-2:ngh,j,k,1)=(ntm1*BC_face(1,2)%U0(-2:ngh,j,k,1)+rho_n(nxmngh-2:nx,j,k))*inn1 &
                                                  + (1.0_wp-eps)*BC_face(1,2)%Uref(j,k,1)

                BC_face(1,2)%U0(-2:ngh,j,k,2)=(ntm1*BC_face(1,2)%U0(-2:ngh,j,k,2)+uu(nxmngh-2:nx,j,k))*inn1    &
                                               + (1.0_wp-eps)*BC_face(1,2)%Uref(j,k,2)
                BC_face(1,2)%U0(-2:ngh,j,k,3)=(ntm1*BC_face(1,2)%U0(-2:ngh,j,k,3)+vv(nxmngh-2:nx,j,k))*inn1    &
                                               + (1.0_wp-eps)*BC_face(1,2)%Uref(j,k,3)
                BC_face(1,2)%U0(-2:ngh,j,k,4)=(ntm1*BC_face(1,2)%U0(-2:ngh,j,k,4)+ww(nxmngh-2:nx,j,k))*inn1    &
                                               + (1.0_wp-eps)*BC_face(1,2)%Uref(j,k,4)
                BC_face(1,2)%U0(-2:ngh,j,k,5)=(ntm1*BC_face(1,2)%U0(-2:ngh,j,k,5)+prs(nxmngh-2:nx,j,k))*inn1   &
                                               + (1.0_wp-eps)*BC_face(1,2)%Uref(j,k,5)
             enddo
          enddo
          ! time-averaged sound speed squared
          do k=1,nz
             do j=1,ny
                BC_face(1,2)%U0(1:ngh,j,k,6)=(ntm1*BC_face(1,2)%U0(1:ngh,j,k,6)+c2_(1:ngh,j,k))*inn1  &
                                                  + (1.0_wp-eps)*BC_face(1,2)%Uref(j,k,6)
             enddo
          enddo
       else
         ! time-averaged primitive variables
         do k=1,nz
            BC_face(1,2)%U0(-2:ngh,:,k,1)=(ntm1*BC_face(1,2)%U0(-2:ngh,:,k,1)+rho_n(nxmngh-2:nx,:,k))*inn1
            BC_face(1,2)%U0(-2:ngh,:,k,2)=(ntm1*BC_face(1,2)%U0(-2:ngh,:,k,2)+uu(nxmngh-2:nx,:,k))*inn1
            BC_face(1,2)%U0(-2:ngh,:,k,3)=(ntm1*BC_face(1,2)%U0(-2:ngh,:,k,3)+vv(nxmngh-2:nx,:,k))*inn1
            BC_face(1,2)%U0(-2:ngh,:,k,4)=(ntm1*BC_face(1,2)%U0(-2:ngh,:,k,4)+ww(nxmngh-2:nx,:,k))*inn1
            BC_face(1,2)%U0(-2:ngh,:,k,5)=(ntm1*BC_face(1,2)%U0(-2:ngh,:,k,5)+prs(nxmngh-2:nx,:,k))*inn1
         enddo
         ! time-averaged sound speed squared
         do k=1,nz
            BC_face(1,2)%U0(1:ngh,1:ny,k,6)=(ntm1*BC_face(1,2)%U0(1:ngh,1:ny,k,6)+c2_(1:ngh,1:ny,k))*inn1
         enddo
       endif
    endif

    ! Compute fluctuations
    ! ====================
    do i=nxmngh-2,nx
       l=i-nxmngh
       do j=ny1,ny2
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
       do j=ndy,nfy
          do k=1,nz
             ! vg(l,j,k)= BC_face(1,2)%U0(l,j,k,2)*cosphi(l,j)+BC_face(1,2)%U0(l,j,k,3)*sinphi(l,j) &
             !      + sqrt(BC_face(1,2)%U0(l,j,k,6)-BC_face(1,2)%U0(l,j,k,4)**2 &
             !      -(BC_face(1,2)%U0(l,j,k,2)*sinphi(l,j)-BC_face(1,2)%U0(l,j,k,3)*cosphi(l,j))**2)
             ! Protection of value with abs() if supersonic at outlet
             vg(l,j,k)= BC_face(1,2)%U0(l,j,k,2)*cosphi(l,j)+BC_face(1,2)%U0(l,j,k,3)*sinphi(l,j) &
                  + sqrt(abs(BC_face(1,2)%U0(l,j,k,6)-BC_face(1,2)%U0(l,j,k,4)**2 &
                  -(BC_face(1,2)%U0(l,j,k,2)*sinphi(l,j)-BC_face(1,2)%U0(l,j,k,3)*cosphi(l,j))**2))
          enddo
       enddo
    enddo

    ! Compute vg*[cos(phi)*dq/dx+sin(phi)*dq/dy+q/r]
    ! ==============================================
    ! (Tam & Webb DRP schemes)
    do i=nxmnghp1,nx-3
       l=i-nxmngh
       do j=ndy,nfy
          do k=1,nz
             duksi = a7(1)* (uf(l+1,j,k)-uf(l-1,j,k)) &
                   + a7(2)* (uf(l+2,j,k)-uf(l-2,j,k)) &
                   + a7(3)* (uf(l+3,j,k)-uf(l-3,j,k))
             dvksi = a7(1)* (vf(l+1,j,k)-vf(l-1,j,k)) &
                   + a7(2)* (vf(l+2,j,k)-vf(l-2,j,k)) &
                   + a7(3)* (vf(l+3,j,k)-vf(l-3,j,k))
             dwksi = a7(1)* (wf(l+1,j,k)-wf(l-1,j,k)) &
                   + a7(2)* (wf(l+2,j,k)-wf(l-2,j,k)) &
                   + a7(3)* (wf(l+3,j,k)-wf(l-3,j,k))
             dpksi = a7(1)* (pf(l+1,j,k)-pf(l-1,j,k)) &
                   + a7(2)* (pf(l+2,j,k)-pf(l-2,j,k)) &
                   + a7(3)* (pf(l+3,j,k)-pf(l-3,j,k))
             drksi = a7(1)* (rf(l+1,j,k)-rf(l-1,j,k)) &
                   + a7(2)* (rf(l+2,j,k)-rf(l-2,j,k)) &
                   + a7(3)* (rf(l+3,j,k)-rf(l-3,j,k))
             
             dueta = a7(1)*( uf(l,j+1,k) - uf(l,j-1,k) ) + &
                     a7(2)*( uf(l,j+2,k) - uf(l,j-2,k) ) + &
                     a7(3)*( uf(l,j+3,k) - uf(l,j-3,k) )
             dveta = a7(1)*( vf(l,j+1,k) - vf(l,j-1,k) ) + &
                     a7(2)*( vf(l,j+2,k) - vf(l,j-2,k) ) + &
                     a7(3)*( vf(l,j+3,k) - vf(l,j-3,k) )
             dweta = a7(1)*( wf(l,j+1,k) - wf(l,j-1,k) ) + &
                     a7(2)*( wf(l,j+2,k) - wf(l,j-2,k) ) + &
                     a7(3)*( wf(l,j+3,k) - wf(l,j-3,k) )
             dpeta = a7(1)*( pf(l,j+1,k) - pf(l,j-1,k) ) + &
                     a7(2)*( pf(l,j+2,k) - pf(l,j-2,k) ) + &
                     a7(3)*( pf(l,j+3,k) - pf(l,j-3,k) )
             dreta = a7(1)*( rf(l,j+1,k) - rf(l,j-1,k) ) + &
                     a7(2)*( rf(l,j+2,k) - rf(l,j-2,k) ) + &
                     a7(3)*( rf(l,j+3,k) - rf(l,j-3,k) )
          
             pt(l,j,k) = vg(l,j,k)*( ((dpksi*y_eta(i,j)-dpeta*y_ksi(i,j))*cosphi(l,j) &
                                     +(dpeta*x_ksi(i,j)-dpksi*x_eta(i,j))*sinphi(l,j))*ijacob(i,j) &
                                     + pf(l,j,k)*ir(l,j) )
             ut(l,j,k) = vg(l,j,k)*( ((duksi*y_eta(i,j)-dueta*y_ksi(i,j))*cosphi(l,j) &
                                     +(dueta*x_ksi(i,j)-duksi*x_eta(i,j))*sinphi(l,j))*ijacob(i,j) &
                                     + uf(l,j,k)*ir(l,j) )
             vt(l,j,k) = vg(l,j,k)*( ((dvksi*y_eta(i,j)-dveta*y_ksi(i,j))*cosphi(l,j) &
                                     +(dveta*x_ksi(i,j)-dvksi*x_eta(i,j))*sinphi(l,j))*ijacob(i,j) &
                                     + vf(l,j,k)*ir(l,j) )
             wt(l,j,k) = vg(l,j,k)*( ((dwksi*y_eta(i,j)-dweta*y_ksi(i,j))*cosphi(l,j) &
                                     +(dweta*x_ksi(i,j)-dwksi*x_eta(i,j))*sinphi(l,j))*ijacob(i,j) &
                                     + wf(l,j,k)*ir(l,j) )
             rt(l,j,k) = vg(l,j,k)*( ((drksi*y_eta(i,j)-dreta*y_ksi(i,j))*cosphi(l,j) &
                                     +(dreta*x_ksi(i,j)-drksi*x_eta(i,j))*sinphi(l,j))*ijacob(i,j) &
                                     + rf(l,j,k)*ir(l,j) )
          enddo
       enddo
    enddo
    
    i=nx-2
    l=i-nxmngh
    do j=ndy,nfy
       do k=1,nz
          duksi = a42(1)*uf(l+2,j,k)+a42(2)*uf(l+1,j,k) &
                + a42(3)*uf(l  ,j,k)+a42(4)*uf(l-1,j,k) &
                + a42(5)*uf(l-2,j,k)+a42(6)*uf(l-3,j,k) &
                + a42(7)*uf(l-4,j,k)
          dvksi = a42(1)*vf(l+2,j,k)+a42(2)*vf(l+1,j,k) &
                + a42(3)*vf(l  ,j,k)+a42(4)*vf(l-1,j,k) &
                + a42(5)*vf(l-2,j,k)+a42(6)*vf(l-3,j,k) &
                + a42(7)*vf(l-4,j,k)
          dwksi = a42(1)*wf(l+2,j,k)+a42(2)*wf(l+1,j,k) &
                + a42(3)*wf(l  ,j,k)+a42(4)*wf(l-1,j,k) &
                + a42(5)*wf(l-2,j,k)+a42(6)*wf(l-3,j,k) &
                + a42(7)*wf(l-4,j,k)
          dpksi = a42(1)*pf(l+2,j,k)+a42(2)*pf(l+1,j,k) &
                + a42(3)*pf(l  ,j,k)+a42(4)*pf(l-1,j,k) &
                + a42(5)*pf(l-2,j,k)+a42(6)*pf(l-3,j,k) &
                + a42(7)*pf(l-4,j,k)
          drksi = a42(1)*rf(l+2,j,k)+a42(2)*rf(l+1,j,k) &
                + a42(3)*rf(l  ,j,k)+a42(4)*rf(l-1,j,k) &
                + a42(5)*rf(l-2,j,k)+a42(6)*rf(l-3,j,k) &
                + a42(7)*rf(l-4,j,k)

          dueta = a7(1)*( uf(l,j+1,k) - uf(l,j-1,k) ) + &
                  a7(2)*( uf(l,j+2,k) - uf(l,j-2,k) ) + &
                  a7(3)*( uf(l,j+3,k) - uf(l,j-3,k) )
          dveta = a7(1)*( vf(l,j+1,k) - vf(l,j-1,k) ) + &
                  a7(2)*( vf(l,j+2,k) - vf(l,j-2,k) ) + &
                  a7(3)*( vf(l,j+3,k) - vf(l,j-3,k) )
          dweta = a7(1)*( wf(l,j+1,k) - wf(l,j-1,k) ) + &
                  a7(2)*( wf(l,j+2,k) - wf(l,j-2,k) ) + &
                  a7(3)*( wf(l,j+3,k) - wf(l,j-3,k) )
          dpeta = a7(1)*( pf(l,j+1,k) - pf(l,j-1,k) ) + &
                  a7(2)*( pf(l,j+2,k) - pf(l,j-2,k) ) + &
                  a7(3)*( pf(l,j+3,k) - pf(l,j-3,k) )
          dreta = a7(1)*( rf(l,j+1,k) - rf(l,j-1,k) ) + &
                  a7(2)*( rf(l,j+2,k) - rf(l,j-2,k) ) + &
                  a7(3)*( rf(l,j+3,k) - rf(l,j-3,k) )
          
          pt(l,j,k) = vg(l,j,k)*( ((dpksi*y_eta(i,j)-dpeta*y_ksi(i,j))*cosphi(l,j) &
                                  +(dpeta*x_ksi(i,j)-dpksi*x_eta(i,j))*sinphi(l,j))*ijacob(i,j) &
                                  + pf(l,j,k)*ir(l,j) )
          ut(l,j,k) = vg(l,j,k)*( ((duksi*y_eta(i,j)-dueta*y_ksi(i,j))*cosphi(l,j) &
                                  +(dueta*x_ksi(i,j)-duksi*x_eta(i,j))*sinphi(l,j))*ijacob(i,j) &
                                  + uf(l,j,k)*ir(l,j) )
          vt(l,j,k) = vg(l,j,k)*( ((dvksi*y_eta(i,j)-dveta*y_ksi(i,j))*cosphi(l,j) &
                                  +(dveta*x_ksi(i,j)-dvksi*x_eta(i,j))*sinphi(l,j))*ijacob(i,j) &
                                  + vf(l,j,k)*ir(l,j) )
          wt(l,j,k) = vg(l,j,k)*( ((dwksi*y_eta(i,j)-dweta*y_ksi(i,j))*cosphi(l,j) &
                                  +(dweta*x_ksi(i,j)-dwksi*x_eta(i,j))*sinphi(l,j))*ijacob(i,j) &
                                  + wf(l,j,k)*ir(l,j) )
          rt(l,j,k) = vg(l,j,k)*( ((drksi*y_eta(i,j)-dreta*y_ksi(i,j))*cosphi(l,j) &
                                  +(dreta*x_ksi(i,j)-drksi*x_eta(i,j))*sinphi(l,j))*ijacob(i,j) &
                                  + rf(l,j,k)*ir(l,j) )
       enddo
    enddo
    
    i=nx-1
    l=i-nxmngh
    do j=ndy,nfy
       do k=1,nz
          duksi = a51(1)*uf(l+1,j,k)+a51(2)*uf(l  ,j,k) &
                + a51(3)*uf(l-1,j,k)+a51(4)*uf(l-2,j,k) &
                + a51(5)*uf(l-3,j,k)+a51(6)*uf(l-4,j,k) &
                + a51(7)*uf(l-5,j,k)
          dvksi = a51(1)*vf(l+1,j,k)+a51(2)*vf(l  ,j,k) &
                + a51(3)*vf(l-1,j,k)+a51(4)*vf(l-2,j,k) &
                + a51(5)*vf(l-3,j,k)+a51(6)*vf(l-4,j,k) &
                + a51(7)*vf(l-5,j,k)
          dwksi = a51(1)*wf(l+1,j,k)+a51(2)*wf(l  ,j,k) &
                + a51(3)*wf(l-1,j,k)+a51(4)*wf(l-2,j,k) &
                + a51(5)*wf(l-3,j,k)+a51(6)*wf(l-4,j,k) &
                + a51(7)*wf(l-5,j,k)
          dpksi = a51(1)*pf(l+1,j,k)+a51(2)*pf(l  ,j,k) &
                + a51(3)*pf(l-1,j,k)+a51(4)*pf(l-2,j,k) &
                + a51(5)*pf(l-3,j,k)+a51(6)*pf(l-4,j,k) &
                + a51(7)*pf(l-5,j,k)
          drksi = a51(1)*rf(l+1,j,k)+a51(2)*rf(l  ,j,k) &
                + a51(3)*rf(l-1,j,k)+a51(4)*rf(l-2,j,k) &
                + a51(5)*rf(l-3,j,k)+a51(6)*rf(l-4,j,k) &
                + a51(7)*rf(l-5,j,k)

          dueta = a7(1)*( uf(l,j+1,k) - uf(l,j-1,k) ) + &
                  a7(2)*( uf(l,j+2,k) - uf(l,j-2,k) ) + &
                  a7(3)*( uf(l,j+3,k) - uf(l,j-3,k) )
          dveta = a7(1)*( vf(l,j+1,k) - vf(l,j-1,k) ) + &
                  a7(2)*( vf(l,j+2,k) - vf(l,j-2,k) ) + &
                  a7(3)*( vf(l,j+3,k) - vf(l,j-3,k) )
          dweta = a7(1)*( wf(l,j+1,k) - wf(l,j-1,k) ) + &
                  a7(2)*( wf(l,j+2,k) - wf(l,j-2,k) ) + &
                  a7(3)*( wf(l,j+3,k) - wf(l,j-3,k) )
          dpeta = a7(1)*( pf(l,j+1,k) - pf(l,j-1,k) ) + &
                  a7(2)*( pf(l,j+2,k) - pf(l,j-2,k) ) + &
                  a7(3)*( pf(l,j+3,k) - pf(l,j-3,k) )
          dreta = a7(1)*( rf(l,j+1,k) - rf(l,j-1,k) ) + &
                  a7(2)*( rf(l,j+2,k) - rf(l,j-2,k) ) + &
                  a7(3)*( rf(l,j+3,k) - rf(l,j-3,k) )
          
          pt(l,j,k) = vg(l,j,k)*( ((dpksi*y_eta(i,j)-dpeta*y_ksi(i,j))*cosphi(l,j) &
                                  +(dpeta*x_ksi(i,j)-dpksi*x_eta(i,j))*sinphi(l,j))*ijacob(i,j) &
                                  + pf(l,j,k)*ir(l,j) )
          ut(l,j,k) = vg(l,j,k)*( ((duksi*y_eta(i,j)-dueta*y_ksi(i,j))*cosphi(l,j) &
                                  +(dueta*x_ksi(i,j)-duksi*x_eta(i,j))*sinphi(l,j))*ijacob(i,j) &
                                  + uf(l,j,k)*ir(l,j) )
          vt(l,j,k) = vg(l,j,k)*( ((dvksi*y_eta(i,j)-dveta*y_ksi(i,j))*cosphi(l,j) &
                                  +(dveta*x_ksi(i,j)-dvksi*x_eta(i,j))*sinphi(l,j))*ijacob(i,j) &
                                  + vf(l,j,k)*ir(l,j) )
          wt(l,j,k) = vg(l,j,k)*( ((dwksi*y_eta(i,j)-dweta*y_ksi(i,j))*cosphi(l,j) &
                                  +(dweta*x_ksi(i,j)-dwksi*x_eta(i,j))*sinphi(l,j))*ijacob(i,j) &
                                  + wf(l,j,k)*ir(l,j) )
          rt(l,j,k) = vg(l,j,k)*( ((drksi*y_eta(i,j)-dreta*y_ksi(i,j))*cosphi(l,j) &
                                  +(dreta*x_ksi(i,j)-drksi*x_eta(i,j))*sinphi(l,j))*ijacob(i,j) &
                                  + rf(l,j,k)*ir(l,j) )
       enddo
    enddo
    
    i=nx
    l=i-nxmngh
    do j=ndy,nfy
       do k=1,nz
          duksi = a60(1)*uf(l  ,j,k)+a60(2)*uf(l-1,j,k) &
                + a60(3)*uf(l-2,j,k)+a60(4)*uf(l-3,j,k) &
                + a60(5)*uf(l-4,j,k)+a60(6)*uf(l-5,j,k) &
                + a60(7)*uf(l-6,j,k)
          dvksi = a60(1)*vf(l  ,j,k)+a60(2)*vf(l-1,j,k) &
                + a60(3)*vf(l-2,j,k)+a60(4)*vf(l-3,j,k) &
                + a60(5)*vf(l-4,j,k)+a60(6)*vf(l-5,j,k) &
                + a60(7)*vf(l-6,j,k)
          dwksi = a60(1)*wf(l  ,j,k)+a60(2)*wf(l-1,j,k) &
                + a60(3)*wf(l-2,j,k)+a60(4)*wf(l-3,j,k) &
                + a60(5)*wf(l-4,j,k)+a60(6)*wf(l-5,j,k) &
                + a60(7)*wf(l-6,j,k)
          dpksi = a60(1)*pf(l  ,j,k)+a60(2)*pf(l-1,j,k) &
                + a60(3)*pf(l-2,j,k)+a60(4)*pf(l-3,j,k) &
                + a60(5)*pf(l-4,j,k)+a60(6)*pf(l-5,j,k) &
                + a60(7)*pf(l-6,j,k)
          drksi = a60(1)*rf(l  ,j,k)+a60(2)*rf(l-1,j,k) &
                + a60(3)*rf(l-2,j,k)+a60(4)*rf(l-3,j,k) &
                + a60(5)*rf(l-4,j,k)+a60(6)*rf(l-5,j,k) &
                + a60(7)*rf(l-6,j,k)
          
          dueta = a7(1)*( uf(l,j+1,k) - uf(l,j-1,k) ) + &
                  a7(2)*( uf(l,j+2,k) - uf(l,j-2,k) ) + &
                  a7(3)*( uf(l,j+3,k) - uf(l,j-3,k) )
          dveta = a7(1)*( vf(l,j+1,k) - vf(l,j-1,k) ) + &
                  a7(2)*( vf(l,j+2,k) - vf(l,j-2,k) ) + &
                  a7(3)*( vf(l,j+3,k) - vf(l,j-3,k) )
          dweta = a7(1)*( wf(l,j+1,k) - wf(l,j-1,k) ) + &
                  a7(2)*( wf(l,j+2,k) - wf(l,j-2,k) ) + &
                  a7(3)*( wf(l,j+3,k) - wf(l,j-3,k) )
          dpeta = a7(1)*( pf(l,j+1,k) - pf(l,j-1,k) ) + &
                  a7(2)*( pf(l,j+2,k) - pf(l,j-2,k) ) + &
                  a7(3)*( pf(l,j+3,k) - pf(l,j-3,k) )
          dreta = a7(1)*( rf(l,j+1,k) - rf(l,j-1,k) ) + &
                  a7(2)*( rf(l,j+2,k) - rf(l,j-2,k) ) + &
                  a7(3)*( rf(l,j+3,k) - rf(l,j-3,k) )
          
          pt(l,j,k) = vg(l,j,k)*( ((dpksi*y_eta(i,j)-dpeta*y_ksi(i,j))*cosphi(l,j) &
                                  +(dpeta*x_ksi(i,j)-dpksi*x_eta(i,j))*sinphi(l,j))*ijacob(i,j) &
                                  + pf(l,j,k)*ir(l,j) )
          ut(l,j,k) = vg(l,j,k)*( ((duksi*y_eta(i,j)-dueta*y_ksi(i,j))*cosphi(l,j) &
                                  +(dueta*x_ksi(i,j)-duksi*x_eta(i,j))*sinphi(l,j))*ijacob(i,j) &
                                  + uf(l,j,k)*ir(l,j) )
          vt(l,j,k) = vg(l,j,k)*( ((dvksi*y_eta(i,j)-dveta*y_ksi(i,j))*cosphi(l,j) &
                                  +(dveta*x_ksi(i,j)-dvksi*x_eta(i,j))*sinphi(l,j))*ijacob(i,j) &
                                  + vf(l,j,k)*ir(l,j) )
          wt(l,j,k) = vg(l,j,k)*( ((dwksi*y_eta(i,j)-dweta*y_ksi(i,j))*cosphi(l,j) &
                                  +(dweta*x_ksi(i,j)-dwksi*x_eta(i,j))*sinphi(l,j))*ijacob(i,j) &
                                  + wf(l,j,k)*ir(l,j) )
          rt(l,j,k) = vg(l,j,k)*( ((drksi*y_eta(i,j)-dreta*y_ksi(i,j))*cosphi(l,j) &
                                  +(dreta*x_ksi(i,j)-drksi*x_eta(i,j))*sinphi(l,j))*ijacob(i,j) &
                                  + rf(l,j,k)*ir(l,j) )
       enddo
    enddo
    
    ! Update fluxes at each RK step
    ! =============================
    do k=1,nz
       do j=ndy,nfy
          do i=nxmnghp1,nx
             cp=cpcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
             av=avcalc_tro(Tmp(i,j,k),rho_n(i,j,k))

             l=i-nxmngh
             Krho(i,j,k)  = rt(l,j,k)
             Krhou(i,j,k) = uu(i,j,k)*rt(l,j,k)+rho_n(i,j,k)*ut(l,j,k)
             Krhov(i,j,k) = vv(i,j,k)*rt(l,j,k)+rho_n(i,j,k)*vt(l,j,k)
             Krhow(i,j,k) = ww(i,j,k)*rt(l,j,k)+rho_n(i,j,k)*wt(l,j,k)
             Krhoe(i,j,k) = cp/av*(pt(l,j,k)/c2_(l,j,k)-rt(l,j,k)) &
                  + (rhoe_n(i,j,k)+prs(i,j,k))/rho_n(i,j,k)*rt(l,j,k) &
                  + rho_n(i,j,k)*(uu(i,j,k)*ut(l,j,k)+vv(i,j,k)*vt(l,j,k)+ww(i,j,k)*wt(l,j,k))
          enddo
       enddo
    enddo

  end subroutine bc_TD2d_imax_c
  
  !===============================================================================
  module subroutine bc_TD2d_jmin_c
  !===============================================================================
    !> 2D Tam & Dong's BC: boundary condition at jmin (bottom) - curvilinear version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: cp,av,inn1,ntm1
    real(wp) :: duksi,dvksi,dwksi,dpksi,drksi
    real(wp) :: dueta,dveta,dweta,dpeta,dreta
    real(wp), dimension(nx1:nx2,1:ngh+3,nz) :: rf,uf,vf,wf,pf
    real(wp), dimension(nx,ngh,nz) :: rt,ut,vt,wt,pt,vg,c2_
    !-------------------------------------------------------------------------
    ! added for is_mean_ref
    integer :: n_moy
    real(wp) :: eps
    !-------------------------------------------------------------------------

    ! Pointers for polar coordinates
    ! ==============================
    ir=>BC_face(2,1)%ir
    cosphi=>BC_face(2,1)%cosphi
    sinphi=>BC_face(2,1)%sinphi

    ! Sound speed
    ! ===========
    do k=1,nz
       do i=1,nx
          do j=1,ngh
             c2_(i,j,k)=c2calc_tro(Tmp(i,j,k),rho_n(i,j,k))
          enddo
       enddo
    enddo

    ! Compute online time-averaged primitive variables
    ! ================================================
    if (irk==nrk) then
       inn1 = 1.0_wp/dble(ntotal)
       ntm1=dble(ntotal-1)
       if (BC_face(2,1)%is_mean_ref) then
          eps=0.1_wp
          n_moy=100
          inn1=inn1*eps
          ! time-averaged primitive variables
          do k=1,nz
             do i=1,nx
                BC_face(2,1)%U0(i,1:nghp3,k,1)=(ntm1*BC_face(2,1)%U0(i,1:nghp3,k,1)+rho_n(i,1:nghp3,k))*inn1 &
                                                  + (1.0_wp-eps)*BC_face(2,1)%Uref(i,k,1)
                BC_face(2,1)%U0(i,1:nghp3,k,2)=(ntm1*BC_face(2,1)%U0(i,1:nghp3,k,2)+uu(i,1:nghp3,k))*inn1 &
                                                  + (1.0_wp-eps)*BC_face(2,1)%Uref(i,k,2)
                BC_face(2,1)%U0(i,1:nghp3,k,3)=(ntm1*BC_face(2,1)%U0(i,1:nghp3,k,3)+vv(i,1:nghp3,k))*inn1 &
                                                  + (1.0_wp-eps)*BC_face(2,1)%Uref(i,k,3)
                BC_face(2,1)%U0(i,1:nghp3,k,4)=(ntm1*BC_face(2,1)%U0(i,1:nghp3,k,4)+ww(i,1:nghp3,k))*inn1 &
                                                  + (1.0_wp-eps)*BC_face(2,1)%Uref(i,k,4)
                BC_face(2,1)%U0(i,1:nghp3,k,5)=(ntm1*BC_face(2,1)%U0(i,1:nghp3,k,5)+prs(i,1:nghp3,k))*inn1 &
                                                  + (1.0_wp-eps)*BC_face(2,1)%Uref(i,k,5)
             enddo
          enddo
          ! time-averaged sound speed squared
          do k=1,nz
             do i=1,nx
                BC_face(2,1)%U0(i,1:ngh,k,6)=(ntm1*BC_face(2,1)%U0(i,1:ngh,k,6)+c2_(i,1:ngh,k))*inn1&
                                                  + (1.0_wp-eps)*BC_face(2,1)%Uref(i,k,6)
             enddo
          enddo
       else
          ! time-averaged primitive variables
          do k=1,nz
             BC_face(2,1)%U0(:,1:nghp3,k,1)=(ntm1*BC_face(2,1)%U0(:,1:nghp3,k,1)+rho_n(:,1:nghp3,k))*inn1
             BC_face(2,1)%U0(:,1:nghp3,k,2)=(ntm1*BC_face(2,1)%U0(:,1:nghp3,k,2)+uu(:,1:nghp3,k))*inn1
             BC_face(2,1)%U0(:,1:nghp3,k,3)=(ntm1*BC_face(2,1)%U0(:,1:nghp3,k,3)+vv(:,1:nghp3,k))*inn1
             BC_face(2,1)%U0(:,1:nghp3,k,4)=(ntm1*BC_face(2,1)%U0(:,1:nghp3,k,4)+ww(:,1:nghp3,k))*inn1
             BC_face(2,1)%U0(:,1:nghp3,k,5)=(ntm1*BC_face(2,1)%U0(:,1:nghp3,k,5)+prs(:,1:nghp3,k))*inn1
          enddo
          ! time-averaged sound speed squared
          do k=1,nz
             BC_face(2,1)%U0(1:nx,1:ngh,k,6)=(ntm1*BC_face(2,1)%U0(1:nx,1:ngh,k,6)+c2_(1:nx,1:ngh,k))*inn1
          enddo
       endif
    endif

    ! Compute fluctuations
    ! ====================
    do j=1,nghp3
       do i=nx1,nx2
          do k=1,nz
             rf(i,j,k)=rho_n(i,j,k)-BC_face(2,1)%U0(i,j,k,1)
             uf(i,j,k)=   uu(i,j,k)-BC_face(2,1)%U0(i,j,k,2)
             vf(i,j,k)=   vv(i,j,k)-BC_face(2,1)%U0(i,j,k,3)
             wf(i,j,k)=   ww(i,j,k)-BC_face(2,1)%U0(i,j,k,4)
             pf(i,j,k)=  prs(i,j,k)-BC_face(2,1)%U0(i,j,k,5)
          enddo
       enddo
    enddo

    ! Compute group velocity vg
    ! =========================
    ! vg=u0*cos(phi)+v0*sin(phi) + sqrt(c0^2-w0^2-(u0*sin(phi)+v0*cos(phi))^2)
    do j=1,ngh
       do i=ndx,nfx
          do k=1,nz
             ! vg(i,j,k)= BC_face(2,1)%U0(i,j,k,2)*cosphi(i,j)+BC_face(2,1)%U0(i,j,k,3)*sinphi(i,j) &
             !      + sqrt(BC_face(2,1)%U0(i,j,k,6)-BC_face(2,1)%U0(i,j,k,4)**2 &
             !      -(BC_face(2,1)%U0(i,j,k,2)*sinphi(i,j)-BC_face(2,1)%U0(i,j,k,3)*cosphi(i,j))**2)
             ! Protection of value with abs() if supersonic at outlet
             vg(i,j,k)= BC_face(2,1)%U0(i,j,k,2)*cosphi(i,j)+BC_face(2,1)%U0(i,j,k,3)*sinphi(i,j) &
                  + sqrt(abs(BC_face(2,1)%U0(i,j,k,6)-BC_face(2,1)%U0(i,j,k,4)**2 &
                  -(BC_face(2,1)%U0(i,j,k,2)*sinphi(i,j)-BC_face(2,1)%U0(i,j,k,3)*cosphi(i,j))**2))
          enddo
       enddo
    enddo

    ! Compute vg*[cos(phi)*dq/dx+sin(phi)*dq/dy+q/r]
    ! ==============================================
    ! (Tam & Webb DRP schemes)
    j=1
    do i=ndx,nfx
       do k=1,nz
          duksi = a7(1)*( uf(i+1,j,k) - uf(i-1,j,k) ) + &
                  a7(2)*( uf(i+2,j,k) - uf(i-2,j,k) ) + &
                  a7(3)*( uf(i+3,j,k) - uf(i-3,j,k) )
          dvksi = a7(1)*( vf(i+1,j,k) - vf(i-1,j,k) ) + &
                  a7(2)*( vf(i+2,j,k) - vf(i-2,j,k) ) + &
                  a7(3)*( vf(i+3,j,k) - vf(i-3,j,k) )
          dwksi = a7(1)*( wf(i+1,j,k) - wf(i-1,j,k) ) + &
                  a7(2)*( wf(i+2,j,k) - wf(i-2,j,k) ) + &
                  a7(3)*( wf(i+3,j,k) - wf(i-3,j,k) )
          dpksi = a7(1)*( pf(i+1,j,k) - pf(i-1,j,k) ) + &
                  a7(2)*( pf(i+2,j,k) - pf(i-2,j,k) ) + &
                  a7(3)*( pf(i+3,j,k) - pf(i-3,j,k) )
          drksi = a7(1)*( rf(i+1,j,k) - rf(i-1,j,k) ) + &
                  a7(2)*( rf(i+2,j,k) - rf(i-2,j,k) ) + &
                  a7(3)*( rf(i+3,j,k) - rf(i-3,j,k) )
          
          dueta = a06(1)*uf(i,1,k)+a06(2)*uf(i,2,k) &
                + a06(3)*uf(i,3,k)+a06(4)*uf(i,4,k) &
                + a06(5)*uf(i,5,k)+a06(6)*uf(i,6,k) &
                + a06(7)*uf(i,7,k)
          dveta = a06(1)*vf(i,1,k)+a06(2)*vf(i,2,k) &
                + a06(3)*vf(i,3,k)+a06(4)*vf(i,4,k) &
                + a06(5)*vf(i,5,k)+a06(6)*vf(i,6,k) &
                + a06(7)*vf(i,7,k)
          dweta = a06(1)*wf(i,1,k)+a06(2)*wf(i,2,k) &
                + a06(3)*wf(i,3,k)+a06(4)*wf(i,4,k) &
                + a06(5)*wf(i,5,k)+a06(6)*wf(i,6,k) &
                + a06(7)*wf(i,7,k)
          dpeta = a06(1)*pf(i,1,k)+a06(2)*pf(i,2,k) &
                + a06(3)*pf(i,3,k)+a06(4)*pf(i,4,k) &
                + a06(5)*pf(i,5,k)+a06(6)*pf(i,6,k) &
                + a06(7)*pf(i,7,k)
          dreta = a06(1)*rf(i,1,k)+a06(2)*rf(i,2,k) &
                + a06(3)*rf(i,3,k)+a06(4)*rf(i,4,k) &
                + a06(5)*rf(i,5,k)+a06(6)*rf(i,6,k) &
                + a06(7)*rf(i,7,k)

          pt(i,j,k) = vg(i,j,k)*( ((dpksi*y_eta(i,j)-dpeta*y_ksi(i,j))*cosphi(i,j) &
                                  +(dpeta*x_ksi(i,j)-dpksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                  + pf(i,j,k)*ir(i,j) )
          ut(i,j,k) = vg(i,j,k)*( ((duksi*y_eta(i,j)-dueta*y_ksi(i,j))*cosphi(i,j) &
                                  +(dueta*x_ksi(i,j)-duksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                  + uf(i,j,k)*ir(i,j) )
          vt(i,j,k) = vg(i,j,k)*( ((dvksi*y_eta(i,j)-dveta*y_ksi(i,j))*cosphi(i,j) &
                                  +(dveta*x_ksi(i,j)-dvksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                  + vf(i,j,k)*ir(i,j) )
          wt(i,j,k) = vg(i,j,k)*( ((dwksi*y_eta(i,j)-dweta*y_ksi(i,j))*cosphi(i,j) &
                                  +(dweta*x_ksi(i,j)-dwksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                  + wf(i,j,k)*ir(i,j) )
          rt(i,j,k) = vg(i,j,k)*( ((drksi*y_eta(i,j)-dreta*y_ksi(i,j))*cosphi(i,j) &
                                  +(dreta*x_ksi(i,j)-drksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                  + rf(i,j,k)*ir(i,j) )
       enddo
    enddo

    j=2
    do i=ndx,nfx
       do k=1,nz
          duksi = a7(1)*( uf(i+1,j,k) - uf(i-1,j,k) ) + &
                  a7(2)*( uf(i+2,j,k) - uf(i-2,j,k) ) + &
                  a7(3)*( uf(i+3,j,k) - uf(i-3,j,k) )
          dvksi = a7(1)*( vf(i+1,j,k) - vf(i-1,j,k) ) + &
                  a7(2)*( vf(i+2,j,k) - vf(i-2,j,k) ) + &
                  a7(3)*( vf(i+3,j,k) - vf(i-3,j,k) )
          dwksi = a7(1)*( wf(i+1,j,k) - wf(i-1,j,k) ) + &
                  a7(2)*( wf(i+2,j,k) - wf(i-2,j,k) ) + &
                  a7(3)*( wf(i+3,j,k) - wf(i-3,j,k) )
          dpksi = a7(1)*( pf(i+1,j,k) - pf(i-1,j,k) ) + &
                  a7(2)*( pf(i+2,j,k) - pf(i-2,j,k) ) + &
                  a7(3)*( pf(i+3,j,k) - pf(i-3,j,k) )
          drksi = a7(1)*( rf(i+1,j,k) - rf(i-1,j,k) ) + &
                  a7(2)*( rf(i+2,j,k) - rf(i-2,j,k) ) + &
                  a7(3)*( rf(i+3,j,k) - rf(i-3,j,k) )
          
          dueta = a15(1)*uf(i,1,k)+a15(2)*uf(i,2,k) &
                + a15(3)*uf(i,3,k)+a15(4)*uf(i,4,k) &
                + a15(5)*uf(i,5,k)+a15(6)*uf(i,6,k) &
                + a15(7)*uf(i,7,k)
          dveta = a15(1)*vf(i,1,k)+a15(2)*vf(i,2,k) &
                + a15(3)*vf(i,3,k)+a15(4)*vf(i,4,k) &
                + a15(5)*vf(i,5,k)+a15(6)*vf(i,6,k) &
                + a15(7)*vf(i,7,k)
          dweta = a15(1)*wf(i,1,k)+a15(2)*wf(i,2,k) &
                + a15(3)*wf(i,3,k)+a15(4)*wf(i,4,k) &
                + a15(5)*wf(i,5,k)+a15(6)*wf(i,6,k) &
                + a15(7)*wf(i,7,k)
          dpeta = a15(1)*pf(i,1,k)+a15(2)*pf(i,2,k) &
                + a15(3)*pf(i,3,k)+a15(4)*pf(i,4,k) &
                + a15(5)*pf(i,5,k)+a15(6)*pf(i,6,k) &
                + a15(7)*pf(i,7,k)
          dreta = a15(1)*rf(i,1,k)+a15(2)*rf(i,2,k) &
                + a15(3)*rf(i,3,k)+a15(4)*rf(i,4,k) &
                + a15(5)*rf(i,5,k)+a15(6)*rf(i,6,k) &
                + a15(7)*rf(i,7,k)

          pt(i,j,k) = vg(i,j,k)*( ((dpksi*y_eta(i,j)-dpeta*y_ksi(i,j))*cosphi(i,j) &
                                  +(dpeta*x_ksi(i,j)-dpksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                  + pf(i,j,k)*ir(i,j) )
          ut(i,j,k) = vg(i,j,k)*( ((duksi*y_eta(i,j)-dueta*y_ksi(i,j))*cosphi(i,j) &
                                  +(dueta*x_ksi(i,j)-duksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                  + uf(i,j,k)*ir(i,j) )
          vt(i,j,k) = vg(i,j,k)*( ((dvksi*y_eta(i,j)-dveta*y_ksi(i,j))*cosphi(i,j) &
                                  +(dveta*x_ksi(i,j)-dvksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                  + vf(i,j,k)*ir(i,j) )
          wt(i,j,k) = vg(i,j,k)*( ((dwksi*y_eta(i,j)-dweta*y_ksi(i,j))*cosphi(i,j) &
                                  +(dweta*x_ksi(i,j)-dwksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                  + wf(i,j,k)*ir(i,j) )
          rt(i,j,k) = vg(i,j,k)*( ((drksi*y_eta(i,j)-dreta*y_ksi(i,j))*cosphi(i,j) &
                                  +(dreta*x_ksi(i,j)-drksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                  + rf(i,j,k)*ir(i,j) )
       enddo
    enddo

    j=3
    do i=ndx,nfx
       do k=1,nz
          duksi = a7(1)*( uf(i+1,j,k) - uf(i-1,j,k) ) + &
                  a7(2)*( uf(i+2,j,k) - uf(i-2,j,k) ) + &
                  a7(3)*( uf(i+3,j,k) - uf(i-3,j,k) )
          dvksi = a7(1)*( vf(i+1,j,k) - vf(i-1,j,k) ) + &
                  a7(2)*( vf(i+2,j,k) - vf(i-2,j,k) ) + &
                  a7(3)*( vf(i+3,j,k) - vf(i-3,j,k) )
          dwksi = a7(1)*( wf(i+1,j,k) - wf(i-1,j,k) ) + &
                  a7(2)*( wf(i+2,j,k) - wf(i-2,j,k) ) + &
                  a7(3)*( wf(i+3,j,k) - wf(i-3,j,k) )
          dpksi = a7(1)*( pf(i+1,j,k) - pf(i-1,j,k) ) + &
                  a7(2)*( pf(i+2,j,k) - pf(i-2,j,k) ) + &
                  a7(3)*( pf(i+3,j,k) - pf(i-3,j,k) )
          drksi = a7(1)*( rf(i+1,j,k) - rf(i-1,j,k) ) + &
                  a7(2)*( rf(i+2,j,k) - rf(i-2,j,k) ) + &
                  a7(3)*( rf(i+3,j,k) - rf(i-3,j,k) )
          
          dueta = a24(1)*uf(i,1,k)+a24(2)*uf(i,2,k) &
                + a24(3)*uf(i,3,k)+a24(4)*uf(i,4,k) &
                + a24(5)*uf(i,5,k)+a24(6)*uf(i,6,k) &
                + a24(7)*uf(i,7,k)
          dveta = a24(1)*vf(i,1,k)+a24(2)*vf(i,2,k) &
                + a24(3)*vf(i,3,k)+a24(4)*vf(i,4,k) &
                + a24(5)*vf(i,5,k)+a24(6)*vf(i,6,k) &
                + a24(7)*vf(i,7,k)
          dweta = a24(1)*wf(i,1,k)+a24(2)*wf(i,2,k) &
                + a24(3)*wf(i,3,k)+a24(4)*wf(i,4,k) &
                + a24(5)*wf(i,5,k)+a24(6)*wf(i,6,k) &
                + a24(7)*wf(i,7,k)
          dpeta = a24(1)*pf(i,1,k)+a24(2)*pf(i,2,k) &
                + a24(3)*pf(i,3,k)+a24(4)*pf(i,4,k) &
                + a24(5)*pf(i,5,k)+a24(6)*pf(i,6,k) &
                + a24(7)*pf(i,7,k)
          dreta = a24(1)*rf(i,1,k)+a24(2)*rf(i,2,k) &
                + a24(3)*rf(i,3,k)+a24(4)*rf(i,4,k) &
                + a24(5)*rf(i,5,k)+a24(6)*rf(i,6,k) &
                + a24(7)*rf(i,7,k)

          pt(i,j,k) = vg(i,j,k)*( ((dpksi*y_eta(i,j)-dpeta*y_ksi(i,j))*cosphi(i,j) &
                                  +(dpeta*x_ksi(i,j)-dpksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                  + pf(i,j,k)*ir(i,j) )
          ut(i,j,k) = vg(i,j,k)*( ((duksi*y_eta(i,j)-dueta*y_ksi(i,j))*cosphi(i,j) &
                                  +(dueta*x_ksi(i,j)-duksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                  + uf(i,j,k)*ir(i,j) )
          vt(i,j,k) = vg(i,j,k)*( ((dvksi*y_eta(i,j)-dveta*y_ksi(i,j))*cosphi(i,j) &
                                  +(dveta*x_ksi(i,j)-dvksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                  + vf(i,j,k)*ir(i,j) )
          wt(i,j,k) = vg(i,j,k)*( ((dwksi*y_eta(i,j)-dweta*y_ksi(i,j))*cosphi(i,j) &
                                  +(dweta*x_ksi(i,j)-dwksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                  + wf(i,j,k)*ir(i,j) )
          rt(i,j,k) = vg(i,j,k)*( ((drksi*y_eta(i,j)-dreta*y_ksi(i,j))*cosphi(i,j) &
                                  +(dreta*x_ksi(i,j)-drksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                  + rf(i,j,k)*ir(i,j) )
       enddo
    enddo

    do j=4,ngh
       do i=ndx,nfx
          do k=1,nz
             duksi = a7(1)*( uf(i+1,j,k) - uf(i-1,j,k) ) + &
                     a7(2)*( uf(i+2,j,k) - uf(i-2,j,k) ) + &
                     a7(3)*( uf(i+3,j,k) - uf(i-3,j,k) )
             dvksi = a7(1)*( vf(i+1,j,k) - vf(i-1,j,k) ) + &
                     a7(2)*( vf(i+2,j,k) - vf(i-2,j,k) ) + &
                     a7(3)*( vf(i+3,j,k) - vf(i-3,j,k) )
             dwksi = a7(1)*( wf(i+1,j,k) - wf(i-1,j,k) ) + &
                     a7(2)*( wf(i+2,j,k) - wf(i-2,j,k) ) + &
                     a7(3)*( wf(i+3,j,k) - wf(i-3,j,k) )
             dpksi = a7(1)*( pf(i+1,j,k) - pf(i-1,j,k) ) + &
                     a7(2)*( pf(i+2,j,k) - pf(i-2,j,k) ) + &
                     a7(3)*( pf(i+3,j,k) - pf(i-3,j,k) )
             drksi = a7(1)*( rf(i+1,j,k) - rf(i-1,j,k) ) + &
                     a7(2)*( rf(i+2,j,k) - rf(i-2,j,k) ) + &
                     a7(3)*( rf(i+3,j,k) - rf(i-3,j,k) )

             dueta = a7(1)*(uf(i,j+1,k)-uf(i,j-1,k)) &
                   + a7(2)*(uf(i,j+2,k)-uf(i,j-2,k)) &
                   + a7(3)*(uf(i,j+3,k)-uf(i,j-3,k))
             dveta = a7(1)*(vf(i,j+1,k)-vf(i,j-1,k)) &
                   + a7(2)*(vf(i,j+2,k)-vf(i,j-2,k)) &
                   + a7(3)*(vf(i,j+3,k)-vf(i,j-3,k))
             dweta = a7(1)*(wf(i,j+1,k)-wf(i,j-1,k)) &
                   + a7(2)*(wf(i,j+2,k)-wf(i,j-2,k)) &
                   + a7(3)*(wf(i,j+3,k)-wf(i,j-3,k))
             dpeta = a7(1)*(pf(i,j+1,k)-pf(i,j-1,k)) &
                   + a7(2)*(pf(i,j+2,k)-pf(i,j-2,k)) &
                   + a7(3)*(pf(i,j+3,k)-pf(i,j-3,k))
             dreta = a7(1)*(rf(i,j+1,k)-rf(i,j-1,k)) &
                   + a7(2)*(rf(i,j+2,k)-rf(i,j-2,k)) &
                   + a7(3)*(rf(i,j+3,k)-rf(i,j-3,k))

             pt(i,j,k) = vg(i,j,k)*( ((dpksi*y_eta(i,j)-dpeta*y_ksi(i,j))*cosphi(i,j) &
                                     +(dpeta*x_ksi(i,j)-dpksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                     + pf(i,j,k)*ir(i,j) )
             ut(i,j,k) = vg(i,j,k)*( ((duksi*y_eta(i,j)-dueta*y_ksi(i,j))*cosphi(i,j) &
                                     +(dueta*x_ksi(i,j)-duksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                     + uf(i,j,k)*ir(i,j) )
             vt(i,j,k) = vg(i,j,k)*( ((dvksi*y_eta(i,j)-dveta*y_ksi(i,j))*cosphi(i,j) &
                                     +(dveta*x_ksi(i,j)-dvksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                     + vf(i,j,k)*ir(i,j) )
             wt(i,j,k) = vg(i,j,k)*( ((dwksi*y_eta(i,j)-dweta*y_ksi(i,j))*cosphi(i,j) &
                                     +(dweta*x_ksi(i,j)-dwksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                     + wf(i,j,k)*ir(i,j) )
             rt(i,j,k) = vg(i,j,k)*( ((drksi*y_eta(i,j)-dreta*y_ksi(i,j))*cosphi(i,j) &
                                     +(dreta*x_ksi(i,j)-drksi*x_eta(i,j))*sinphi(i,j))*ijacob(i,j) &
                                     + rf(i,j,k)*ir(i,j) )
          enddo
       enddo
    enddo

    ! Update fluxes at each RK step
    ! =============================
    do k=1,nz
       do j=1,ngh
          do i=ndx,nfx
             cp=cpcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
             av=avcalc_tro(Tmp(i,j,k),rho_n(i,j,k))

             Krho(i,j,k)  = rt(i,j,k)
             Krhou(i,j,k) = uu(i,j,k)*rt(i,j,k)+rho_n(i,j,k)*ut(i,j,k)
             Krhov(i,j,k) = vv(i,j,k)*rt(i,j,k)+rho_n(i,j,k)*vt(i,j,k)
             Krhow(i,j,k) = ww(i,j,k)*rt(i,j,k)+rho_n(i,j,k)*wt(i,j,k)
             Krhoe(i,j,k) = cp/av*(pt(i,j,k)/c2_(i,j,k)-rt(i,j,k)) &
                  + (rhoe_n(i,j,k)+prs(i,j,k))/rho_n(i,j,k)*rt(i,j,k) &
                  + rho_n(i,j,k)*(uu(i,j,k)*ut(i,j,k)+vv(i,j,k)*vt(i,j,k)+ww(i,j,k)*wt(i,j,k))
          enddo
       enddo
    enddo

  end subroutine bc_TD2d_jmin_c

  !===============================================================================
  module subroutine bc_TD2d_jmax_c
  !===============================================================================
    !> 2D Tam & Dong's BC: boundary condition at jmax (top) - curvilinear version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l
    real(wp) :: cp,av,inn1,ntm1
    real(wp) :: duksi,dvksi,dwksi,dpksi,drksi
    real(wp) :: dueta,dveta,dweta,dpeta,dreta
    real(wp), dimension(nx1:nx2,-2:ngh,nz) :: rf,uf,vf,wf,pf
    real(wp), dimension(nx,ngh,nz) :: rt,ut,vt,wt,pt,vg,c2_
    !-------------------------------------------------------------------------
    ! added for is_mean_ref
    integer :: n_moy
    real(wp) :: eps
    !-------------------------------------------------------------------------

    ! Pointers for polar coordinates
    ! ==============================
    ir=>BC_face(2,2)%ir
    cosphi=>BC_face(2,2)%cosphi
    sinphi=>BC_face(2,2)%sinphi

    ! Sound speed
    ! ===========
    do k=1,nz
       do i=1,nx
          do j=nymnghp1,ny
             l=j-nymngh
             c2_(i,l,k)=c2calc_tro(Tmp(i,j,k),rho_n(i,j,k))
          enddo
       enddo
    enddo

    ! Compute online time-averaged primitive variables
    ! ================================================
    if (irk==nrk) then
       inn1 = 1.0_wp/dble(ntotal)
       ntm1=dble(ntotal-1)
       if (BC_face(2,2)%is_mean_ref) then
          eps=0.1_wp
          n_moy=100
          inn1=inn1*eps
          ! time-averaged primitive variables
          do k=1,nz
             do i=1,nx
                BC_face(2,2)%U0(i,-2:ngh,k,1)=(ntm1*BC_face(2,2)%U0(i,-2:ngh,k,1)+rho_n(i,nymngh-2:ny,k))*inn1 &
                                                  + (1.0_wp-eps)*BC_face(2,2)%Uref(i,k,1)
                BC_face(2,2)%U0(i,-2:ngh,k,2)=(ntm1*BC_face(2,2)%U0(i,-2:ngh,k,2)+uu(i,nymngh-2:ny,k))*inn1 &
                                                  + (1.0_wp-eps)*BC_face(2,2)%Uref(i,k,2)
                BC_face(2,2)%U0(i,-2:ngh,k,3)=(ntm1*BC_face(2,2)%U0(i,-2:ngh,k,3)+vv(i,nymngh-2:ny,k))*inn1 &
                                                  + (1.0_wp-eps)*BC_face(2,2)%Uref(i,k,3)
                BC_face(2,2)%U0(i,-2:ngh,k,4)=(ntm1*BC_face(2,2)%U0(i,-2:ngh,k,4)+ww(i,nymngh-2:ny,k))*inn1 &
                                                  + (1.0_wp-eps)*BC_face(2,2)%Uref(i,k,4)
                BC_face(2,2)%U0(i,-2:ngh,k,5)=(ntm1*BC_face(2,2)%U0(i,-2:ngh,k,5)+prs(i,nymngh-2:ny,k))*inn1 &
                                                  + (1.0_wp-eps)*BC_face(2,2)%Uref(i,k,5)
             enddo
          enddo
          ! time-averaged sound speed squared
          do k=1,nz
             do i=1,nx
                BC_face(2,2)%U0(i,1:ngh,k,6)=(ntm1*BC_face(2,2)%U0(i,1:ngh,k,6)+c2_(i,1:ngh,k))*inn1&
                                                  + (1.0_wp-eps)*BC_face(2,2)%Uref(i,k,6)
             enddo
          enddo
       else
         ! time-averaged primitive variables
         do k=1,nz
            BC_face(2,2)%U0(:,-2:ngh,k,1)=(ntm1*BC_face(2,2)%U0(:,-2:ngh,k,1)+rho_n(:,nymngh-2:ny,k))*inn1
            BC_face(2,2)%U0(:,-2:ngh,k,2)=(ntm1*BC_face(2,2)%U0(:,-2:ngh,k,2)+uu(:,nymngh-2:ny,k))*inn1
            BC_face(2,2)%U0(:,-2:ngh,k,3)=(ntm1*BC_face(2,2)%U0(:,-2:ngh,k,3)+vv(:,nymngh-2:ny,k))*inn1
            BC_face(2,2)%U0(:,-2:ngh,k,4)=(ntm1*BC_face(2,2)%U0(:,-2:ngh,k,4)+ww(:,nymngh-2:ny,k))*inn1
            BC_face(2,2)%U0(:,-2:ngh,k,5)=(ntm1*BC_face(2,2)%U0(:,-2:ngh,k,5)+prs(:,nymngh-2:ny,k))*inn1
         enddo
         ! time-averaged sound speed squared
         do k=1,nz
            BC_face(2,2)%U0(1:nx,1:ngh,k,6)=(ntm1*BC_face(2,2)%U0(1:nx,1:ngh,k,6)+c2_(1:nx,1:ngh,k))*inn1
         enddo
       endif
    endif

    ! Compute fluctuations
    ! ====================
    do j=nymngh-2,ny
       l=j-nymngh
       do i=nx1,nx2
          do k=1,nz
             rf(i,l,k)=rho_n(i,j,k)-BC_face(2,2)%U0(i,l,k,1)
             uf(i,l,k)=   uu(i,j,k)-BC_face(2,2)%U0(i,l,k,2)
             vf(i,l,k)=   vv(i,j,k)-BC_face(2,2)%U0(i,l,k,3)
             wf(i,l,k)=   ww(i,j,k)-BC_face(2,2)%U0(i,l,k,4)
             pf(i,l,k)=  prs(i,j,k)-BC_face(2,2)%U0(i,l,k,5)
          enddo
       enddo
    enddo

    ! Compute group velocity vg
    ! =========================
    ! vg=u0*cos(phi)+v0*sin(phi) + sqrt(c0^2-w0^2-(u0*sin(phi)+v0*cos(phi))^2)
    do l=1,ngh
       do i=ndx,nfx
          do k=1,nz
             ! vg(i,l,k)= BC_face(2,2)%U0(i,l,k,2)*cosphi(i,l)+BC_face(2,2)%U0(i,l,k,3)*sinphi(i,l) &
             !      + sqrt(BC_face(2,2)%U0(i,l,k,6)-BC_face(2,2)%U0(i,l,k,4)**2 &
             !      -(BC_face(2,2)%U0(i,l,k,2)*sinphi(i,l)-BC_face(2,2)%U0(i,l,k,3)*cosphi(i,l))**2)
             ! Protection of value with abs() if supersonic at outlet
             vg(i,l,k)= BC_face(2,2)%U0(i,l,k,2)*cosphi(i,l)+BC_face(2,2)%U0(i,l,k,3)*sinphi(i,l) &
                  + sqrt(abs(BC_face(2,2)%U0(i,l,k,6)-BC_face(2,2)%U0(i,l,k,4)**2 &
                  -(BC_face(2,2)%U0(i,l,k,2)*sinphi(i,l)-BC_face(2,2)%U0(i,l,k,3)*cosphi(i,l))**2))
          enddo
       enddo
    enddo

    ! Compute vg*[cos(phi)*dq/dx+sin(phi)*dq/dy+q/r]
    ! ==============================================
    ! (Tam & Webb DRP schemes)
    do j=nymnghp1,ny-3
       l=j-nymngh
       do i=ndx,nfx
          do k=1,nz
             duksi = a7(1)*( uf(i+1,l,k) - uf(i-1,l,k) ) + &
                     a7(2)*( uf(i+2,l,k) - uf(i-2,l,k) ) + &
                     a7(3)*( uf(i+3,l,k) - uf(i-3,l,k) )
             dvksi = a7(1)*( vf(i+1,l,k) - vf(i-1,l,k) ) + &
                     a7(2)*( vf(i+2,l,k) - vf(i-2,l,k) ) + &
                     a7(3)*( vf(i+3,l,k) - vf(i-3,l,k) )
             dwksi = a7(1)*( wf(i+1,l,k) - wf(i-1,l,k) ) + &
                     a7(2)*( wf(i+2,l,k) - wf(i-2,l,k) ) + &
                     a7(3)*( wf(i+3,l,k) - wf(i-3,l,k) )
             dpksi = a7(1)*( pf(i+1,l,k) - pf(i-1,l,k) ) + &
                     a7(2)*( pf(i+2,l,k) - pf(i-2,l,k) ) + &
                     a7(3)*( pf(i+3,l,k) - pf(i-3,l,k) )
             drksi = a7(1)*( rf(i+1,l,k) - rf(i-1,l,k) ) + &
                     a7(2)*( rf(i+2,l,k) - rf(i-2,l,k) ) + &
                     a7(3)*( rf(i+3,l,k) - rf(i-3,l,k) )
             
             dueta = a7(1)*(uf(i,l+1,k)-uf(i,l-1,k)) &
                   + a7(2)*(uf(i,l+2,k)-uf(i,l-2,k)) &
                   + a7(3)*(uf(i,l+3,k)-uf(i,l-3,k))
             dveta = a7(1)*(vf(i,l+1,k)-vf(i,l-1,k)) &
                   + a7(2)*(vf(i,l+2,k)-vf(i,l-2,k)) &
                   + a7(3)*(vf(i,l+3,k)-vf(i,l-3,k))
             dweta = a7(1)*(wf(i,l+1,k)-wf(i,l-1,k)) &
                   + a7(2)*(wf(i,l+2,k)-wf(i,l-2,k)) &
                   + a7(3)*(wf(i,l+3,k)-wf(i,l-3,k))
             dpeta = a7(1)*(pf(i,l+1,k)-pf(i,l-1,k)) &
                   + a7(2)*(pf(i,l+2,k)-pf(i,l-2,k)) &
                   + a7(3)*(pf(i,l+3,k)-pf(i,l-3,k))
             dreta = a7(1)*(rf(i,l+1,k)-rf(i,l-1,k)) &
                   + a7(2)*(rf(i,l+2,k)-rf(i,l-2,k)) &
                   + a7(3)*(rf(i,l+3,k)-rf(i,l-3,k))
    
             pt(i,l,k) = vg(i,l,k)*( ((dpksi*y_eta(i,j)-dpeta*y_ksi(i,j))*cosphi(i,l) &
                                     +(dpeta*x_ksi(i,j)-dpksi*x_eta(i,j))*sinphi(i,l))*ijacob(i,j) &
                                     + pf(i,l,k)*ir(i,l) )
             ut(i,l,k) = vg(i,l,k)*( ((duksi*y_eta(i,j)-dueta*y_ksi(i,j))*cosphi(i,l) &
                                     +(dueta*x_ksi(i,j)-duksi*x_eta(i,j))*sinphi(i,l))*ijacob(i,j) &
                                     + uf(i,l,k)*ir(i,l) )
             vt(i,l,k) = vg(i,l,k)*( ((dvksi*y_eta(i,j)-dveta*y_ksi(i,j))*cosphi(i,l) &
                                     +(dveta*x_ksi(i,j)-dvksi*x_eta(i,j))*sinphi(i,l))*ijacob(i,j) &
                                     + vf(i,l,k)*ir(i,l) )
             wt(i,l,k) = vg(i,l,k)*( ((dwksi*y_eta(i,j)-dweta*y_ksi(i,j))*cosphi(i,l) &
                                     +(dweta*x_ksi(i,j)-dwksi*x_eta(i,j))*sinphi(i,l))*ijacob(i,j) &
                                     + wf(i,l,k)*ir(i,l) )
             rt(i,l,k) = vg(i,l,k)*( ((drksi*y_eta(i,j)-dreta*y_ksi(i,j))*cosphi(i,l) &
                                     +(dreta*x_ksi(i,j)-drksi*x_eta(i,j))*sinphi(i,l))*ijacob(i,j) &
                                     + rf(i,l,k)*ir(i,l) )
          enddo
       enddo
    enddo

    j=ny-2
    l=j-nymngh
    do i=ndx,nfx
       do k=1,nz
          duksi = a7(1)*( uf(i+1,l,k) - uf(i-1,l,k) ) + &
                  a7(2)*( uf(i+2,l,k) - uf(i-2,l,k) ) + &
                  a7(3)*( uf(i+3,l,k) - uf(i-3,l,k) )
          dvksi = a7(1)*( vf(i+1,l,k) - vf(i-1,l,k) ) + &
                  a7(2)*( vf(i+2,l,k) - vf(i-2,l,k) ) + &
                  a7(3)*( vf(i+3,l,k) - vf(i-3,l,k) )
          dwksi = a7(1)*( wf(i+1,l,k) - wf(i-1,l,k) ) + &
                  a7(2)*( wf(i+2,l,k) - wf(i-2,l,k) ) + &
                  a7(3)*( wf(i+3,l,k) - wf(i-3,l,k) )
          dpksi = a7(1)*( pf(i+1,l,k) - pf(i-1,l,k) ) + &
                  a7(2)*( pf(i+2,l,k) - pf(i-2,l,k) ) + &
                  a7(3)*( pf(i+3,l,k) - pf(i-3,l,k) )
          drksi = a7(1)*( rf(i+1,l,k) - rf(i-1,l,k) ) + &
                  a7(2)*( rf(i+2,l,k) - rf(i-2,l,k) ) + &
                  a7(3)*( rf(i+3,l,k) - rf(i-3,l,k) )
          
          dueta = a42(1)*uf(i,l+2,k)+a42(2)*uf(i,l+1,k) &
                + a42(3)*uf(i,l  ,k)+a42(4)*uf(i,l-1,k) &
                + a42(5)*uf(i,l-2,k)+a42(6)*uf(i,l-3,k) &
                + a42(7)*uf(i,l-4,k)
          dveta = a42(1)*vf(i,l+2,k)+a42(2)*vf(i,l+1,k) &
                + a42(3)*vf(i,l  ,k)+a42(4)*vf(i,l-1,k) &
                + a42(5)*vf(i,l-2,k)+a42(6)*vf(i,l-3,k) &
                + a42(7)*vf(i,l-4,k)
          dweta = a42(1)*wf(i,l+2,k)+a42(2)*wf(i,l+1,k) &
                + a42(3)*wf(i,l  ,k)+a42(4)*wf(i,l-1,k) &
                + a42(5)*wf(i,l-2,k)+a42(6)*wf(i,l-3,k) &
                + a42(7)*wf(i,l-4,k)
          dpeta = a42(1)*pf(i,l+2,k)+a42(2)*pf(i,l+1,k) &
                + a42(3)*pf(i,l  ,k)+a42(4)*pf(i,l-1,k) &
                + a42(5)*pf(i,l-2,k)+a42(6)*pf(i,l-3,k) &
                + a42(7)*pf(i,l-4,k)
          dreta = a42(1)*rf(i,l+2,k)+a42(2)*rf(i,l+1,k) &
                + a42(3)*rf(i,l  ,k)+a42(4)*rf(i,l-1,k) &
                + a42(5)*rf(i,l-2,k)+a42(6)*rf(i,l-3,k) &
                + a42(7)*rf(i,l-4,k)

          pt(i,l,k) = vg(i,l,k)*( ((dpksi*y_eta(i,j)-dpeta*y_ksi(i,j))*cosphi(i,l) &
                                  +(dpeta*x_ksi(i,j)-dpksi*x_eta(i,j))*sinphi(i,l))*ijacob(i,j) &
                                  + pf(i,l,k)*ir(i,l) )
          ut(i,l,k) = vg(i,l,k)*( ((duksi*y_eta(i,j)-dueta*y_ksi(i,j))*cosphi(i,l) &
                                  +(dueta*x_ksi(i,j)-duksi*x_eta(i,j))*sinphi(i,l))*ijacob(i,j) &
                                  + uf(i,l,k)*ir(i,l) )
          vt(i,l,k) = vg(i,l,k)*( ((dvksi*y_eta(i,j)-dveta*y_ksi(i,j))*cosphi(i,l) &
                                  +(dveta*x_ksi(i,j)-dvksi*x_eta(i,j))*sinphi(i,l))*ijacob(i,j) &
                                  + vf(i,l,k)*ir(i,l) )
          wt(i,l,k) = vg(i,l,k)*( ((dwksi*y_eta(i,j)-dweta*y_ksi(i,j))*cosphi(i,l) &
                                  +(dweta*x_ksi(i,j)-dwksi*x_eta(i,j))*sinphi(i,l))*ijacob(i,j) &
                                  + wf(i,l,k)*ir(i,l) )
          rt(i,l,k) = vg(i,l,k)*( ((drksi*y_eta(i,j)-dreta*y_ksi(i,j))*cosphi(i,l) &
                                  +(dreta*x_ksi(i,j)-drksi*x_eta(i,j))*sinphi(i,l))*ijacob(i,j) &
                                  + rf(i,l,k)*ir(i,l) )
       enddo
    enddo
    
    j=ny-1
    l=j-nymngh
    do i=ndx,nfx
       do k=1,nz
          duksi = a7(1)*( uf(i+1,l,k) - uf(i-1,l,k) ) + &
                  a7(2)*( uf(i+2,l,k) - uf(i-2,l,k) ) + &
                  a7(3)*( uf(i+3,l,k) - uf(i-3,l,k) )
          dvksi = a7(1)*( vf(i+1,l,k) - vf(i-1,l,k) ) + &
                  a7(2)*( vf(i+2,l,k) - vf(i-2,l,k) ) + &
                  a7(3)*( vf(i+3,l,k) - vf(i-3,l,k) )
          dwksi = a7(1)*( wf(i+1,l,k) - wf(i-1,l,k) ) + &
                  a7(2)*( wf(i+2,l,k) - wf(i-2,l,k) ) + &
                  a7(3)*( wf(i+3,l,k) - wf(i-3,l,k) )
          dpksi = a7(1)*( pf(i+1,l,k) - pf(i-1,l,k) ) + &
                  a7(2)*( pf(i+2,l,k) - pf(i-2,l,k) ) + &
                  a7(3)*( pf(i+3,l,k) - pf(i-3,l,k) )
          drksi = a7(1)*( rf(i+1,l,k) - rf(i-1,l,k) ) + &
                  a7(2)*( rf(i+2,l,k) - rf(i-2,l,k) ) + &
                  a7(3)*( rf(i+3,l,k) - rf(i-3,l,k) )
          
          dueta = a51(1)*uf(i,l+1,k)+a51(2)*uf(i,l  ,k) &
                + a51(3)*uf(i,l-1,k)+a51(4)*uf(i,l-2,k) &
                + a51(5)*uf(i,l-3,k)+a51(6)*uf(i,l-4,k) &
                + a51(7)*uf(i,l-5,k)
          dveta = a51(1)*vf(i,l+1,k)+a51(2)*vf(i,l  ,k) &
                + a51(3)*vf(i,l-1,k)+a51(4)*vf(i,l-2,k) &
                + a51(5)*vf(i,l-3,k)+a51(6)*vf(i,l-4,k) &
                + a51(7)*vf(i,l-5,k)
          dweta = a51(1)*wf(i,l+1,k)+a51(2)*wf(i,l  ,k) &
                + a51(3)*wf(i,l-1,k)+a51(4)*wf(i,l-2,k) &
                + a51(5)*wf(i,l-3,k)+a51(6)*wf(i,l-4,k) &
                + a51(7)*wf(i,l-5,k)
          dpeta = a51(1)*pf(i,l+1,k)+a51(2)*pf(i,l  ,k) &
                + a51(3)*pf(i,l-1,k)+a51(4)*pf(i,l-2,k) &
                + a51(5)*pf(i,l-3,k)+a51(6)*pf(i,l-4,k) &
                + a51(7)*pf(i,l-5,k)
          dreta = a51(1)*rf(i,l+1,k)+a51(2)*rf(i,l  ,k) &
                + a51(3)*rf(i,l-1,k)+a51(4)*rf(i,l-2,k) &
                + a51(5)*rf(i,l-3,k)+a51(6)*rf(i,l-4,k) &
                + a51(7)*rf(i,l-5,k)
          
          pt(i,l,k) = vg(i,l,k)*( ((dpksi*y_eta(i,j)-dpeta*y_ksi(i,j))*cosphi(i,l) &
                                  +(dpeta*x_ksi(i,j)-dpksi*x_eta(i,j))*sinphi(i,l))*ijacob(i,j) &
                                  + pf(i,l,k)*ir(i,l) )
          ut(i,l,k) = vg(i,l,k)*( ((duksi*y_eta(i,j)-dueta*y_ksi(i,j))*cosphi(i,l) &
                                  +(dueta*x_ksi(i,j)-duksi*x_eta(i,j))*sinphi(i,l))*ijacob(i,j) &
                                  + uf(i,l,k)*ir(i,l) )
          vt(i,l,k) = vg(i,l,k)*( ((dvksi*y_eta(i,j)-dveta*y_ksi(i,j))*cosphi(i,l) &
                                  +(dveta*x_ksi(i,j)-dvksi*x_eta(i,j))*sinphi(i,l))*ijacob(i,j) &
                                  + vf(i,l,k)*ir(i,l) )
          wt(i,l,k) = vg(i,l,k)*( ((dwksi*y_eta(i,j)-dweta*y_ksi(i,j))*cosphi(i,l) &
                                  +(dweta*x_ksi(i,j)-dwksi*x_eta(i,j))*sinphi(i,l))*ijacob(i,j) &
                                  + wf(i,l,k)*ir(i,l) )
          rt(i,l,k) = vg(i,l,k)*( ((drksi*y_eta(i,j)-dreta*y_ksi(i,j))*cosphi(i,l) &
                                  +(dreta*x_ksi(i,j)-drksi*x_eta(i,j))*sinphi(i,l))*ijacob(i,j) &
                                  + rf(i,l,k)*ir(i,l) )
       enddo
    enddo
    
    j=ny
    l=j-nymngh
    do i=ndx,nfx
       do k=1,nz
          duksi = a7(1)*( uf(i+1,l,k) - uf(i-1,l,k) ) + &
                  a7(2)*( uf(i+2,l,k) - uf(i-2,l,k) ) + &
                  a7(3)*( uf(i+3,l,k) - uf(i-3,l,k) )
          dvksi = a7(1)*( vf(i+1,l,k) - vf(i-1,l,k) ) + &
                  a7(2)*( vf(i+2,l,k) - vf(i-2,l,k) ) + &
                  a7(3)*( vf(i+3,l,k) - vf(i-3,l,k) )
          dwksi = a7(1)*( wf(i+1,l,k) - wf(i-1,l,k) ) + &
                  a7(2)*( wf(i+2,l,k) - wf(i-2,l,k) ) + &
                  a7(3)*( wf(i+3,l,k) - wf(i-3,l,k) )
          dpksi = a7(1)*( pf(i+1,l,k) - pf(i-1,l,k) ) + &
                  a7(2)*( pf(i+2,l,k) - pf(i-2,l,k) ) + &
                  a7(3)*( pf(i+3,l,k) - pf(i-3,l,k) )
          drksi = a7(1)*( rf(i+1,l,k) - rf(i-1,l,k) ) + &
                  a7(2)*( rf(i+2,l,k) - rf(i-2,l,k) ) + &
                  a7(3)*( rf(i+3,l,k) - rf(i-3,l,k) )
          
          dueta = a60(1)*uf(i,l  ,k)+a60(2)*uf(i,l-1,k) &
                + a60(3)*uf(i,l-2,k)+a60(4)*uf(i,l-3,k) &
                + a60(5)*uf(i,l-4,k)+a60(6)*uf(i,l-5,k) &
                + a60(7)*uf(i,l-6,k)
          dveta = a60(1)*vf(i,l  ,k)+a60(2)*vf(i,l-1,k) &
                + a60(3)*vf(i,l-2,k)+a60(4)*vf(i,l-3,k) &
                + a60(5)*vf(i,l-4,k)+a60(6)*vf(i,l-5,k) &
                + a60(7)*vf(i,l-6,k)
          dweta = a60(1)*wf(i,l  ,k)+a60(2)*wf(i,l-1,k) &
                + a60(3)*wf(i,l-2,k)+a60(4)*wf(i,l-3,k) &
                + a60(5)*wf(i,l-4,k)+a60(6)*wf(i,l-5,k) &
                + a60(7)*wf(i,l-6,k)
          dpeta = a60(1)*pf(i,l  ,k)+a60(2)*pf(i,l-1,k) &
                + a60(3)*pf(i,l-2,k)+a60(4)*pf(i,l-3,k) &
                + a60(5)*pf(i,l-4,k)+a60(6)*pf(i,l-5,k) &
                + a60(7)*pf(i,l-6,k)
          dreta = a60(1)*rf(i,l  ,k)+a60(2)*rf(i,l-1,k) &
                + a60(3)*rf(i,l-2,k)+a60(4)*rf(i,l-3,k) &
                + a60(5)*rf(i,l-4,k)+a60(6)*rf(i,l-5,k) &
                + a60(7)*rf(i,l-6,k)
          
          pt(i,l,k) = vg(i,l,k)*( ((dpksi*y_eta(i,j)-dpeta*y_ksi(i,j))*cosphi(i,l) &
                                  +(dpeta*x_ksi(i,j)-dpksi*x_eta(i,j))*sinphi(i,l))*ijacob(i,j) &
                                  + pf(i,l,k)*ir(i,l) )
          ut(i,l,k) = vg(i,l,k)*( ((duksi*y_eta(i,j)-dueta*y_ksi(i,j))*cosphi(i,l) &
                                  +(dueta*x_ksi(i,j)-duksi*x_eta(i,j))*sinphi(i,l))*ijacob(i,j) &
                                  + uf(i,l,k)*ir(i,l) )
          vt(i,l,k) = vg(i,l,k)*( ((dvksi*y_eta(i,j)-dveta*y_ksi(i,j))*cosphi(i,l) &
                                  +(dveta*x_ksi(i,j)-dvksi*x_eta(i,j))*sinphi(i,l))*ijacob(i,j) &
                                  + vf(i,l,k)*ir(i,l) )
          wt(i,l,k) = vg(i,l,k)*( ((dwksi*y_eta(i,j)-dweta*y_ksi(i,j))*cosphi(i,l) &
                                  +(dweta*x_ksi(i,j)-dwksi*x_eta(i,j))*sinphi(i,l))*ijacob(i,j) &
                                  + wf(i,l,k)*ir(i,l) )
          rt(i,l,k) = vg(i,l,k)*( ((drksi*y_eta(i,j)-dreta*y_ksi(i,j))*cosphi(i,l) &
                                  +(dreta*x_ksi(i,j)-drksi*x_eta(i,j))*sinphi(i,l))*ijacob(i,j) &
                                  + rf(i,l,k)*ir(i,l) )
       enddo
    enddo
    
    ! Update fluxes at each RK step
    ! =============================
    do k=1,nz
       do j=nymnghp1,ny
          l=j-nymngh
          do i=ndx,nfx
             cp=cpcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
             av=avcalc_tro(Tmp(i,j,k),rho_n(i,j,k))

             Krho(i,j,k)  = rt(i,l,k)
             Krhou(i,j,k) = uu(i,j,k)*rt(i,l,k)+rho_n(i,j,k)*ut(i,l,k)
             Krhov(i,j,k) = vv(i,j,k)*rt(i,l,k)+rho_n(i,j,k)*vt(i,l,k)
             Krhow(i,j,k) = ww(i,j,k)*rt(i,l,k)+rho_n(i,j,k)*wt(i,l,k)
             Krhoe(i,j,k) = cp/av*(pt(i,l,k)/c2_(i,l,k)-rt(i,l,k)) &
                  + (rhoe_n(i,j,k)+prs(i,j,k))/rho_n(i,j,k)*rt(i,l,k) &
                  + rho_n(i,j,k)*(uu(i,j,k)*ut(i,l,k)+vv(i,j,k)*vt(i,l,k)+ww(i,j,k)*wt(i,l,k))
          enddo
       enddo
    enddo

  end subroutine bc_TD2d_jmax_c

end submodule smod_TamDong2d_faces_c
