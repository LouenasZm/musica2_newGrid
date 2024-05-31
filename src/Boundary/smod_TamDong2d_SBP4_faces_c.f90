!===============================================================
submodule (mod_TamDong2d_c) smod_TamDong2d_SBP4_faces_c
!===============================================================
  !> author: XG
  !> date: February 2020 - modif January 2022
  !> Non-reflecting boundary conditions of Tam and Dong
  !> - 2D (periodic) curvilinear version - routines forr faces
  !> /!\ dev only: version with SBP4 boundary schemes [not for regular use]
!===============================================================

contains

  !===============================================================================
  module subroutine bc_TD2d_imin_c_SBP4
  !===============================================================================
    !> 2D Tam & Dong's BC: boundary condition at imin (left) - curvilinear version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: cp,av,inn1,ntm1
    real(wp) :: duksi,dvksi,dwksi,dpksi,drksi
    real(wp) :: dueta,dveta,dweta,dpeta,dreta
    real(wp), dimension(1:ngh+3,ny1:ny2,nz) :: rf,uf,vf,wf,pf
    real(wp), dimension(ngh,ny,nz) :: rt,ut,vt,wt,pt,vg,c2_
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
             vg(i,j,k)= BC_face(1,1)%U0(i,j,k,2)*cosphi(i,j)+BC_face(1,1)%U0(i,j,k,3)*sinphi(i,j) &
                  + sqrt(BC_face(1,1)%U0(i,j,k,6)-BC_face(1,1)%U0(i,j,k,4)**2 &
                  -(BC_face(1,1)%U0(i,j,k,2)*sinphi(i,j)-BC_face(1,1)%U0(i,j,k,3)*cosphi(i,j))**2)
          enddo
       enddo
    enddo

    ! Compute vg*[cos(phi)*dq/dx+sin(phi)*dq/dy+q/r]
    ! ==============================================
    ! (Tam & Webb DRP schemes)
    i=1
    do j=ndy,nfy
       do k=1,nz
          duksi = as4p0(1)*uf(1,j,k)+as4p0(2)*uf(2,j,k) &
                + as4p0(3)*uf(3,j,k)+as4p0(4)*uf(4,j,k)
          dvksi = as4p0(1)*vf(1,j,k)+as4p0(2)*vf(2,j,k) &
                + as4p0(3)*vf(3,j,k)+as4p0(4)*vf(4,j,k)
          dwksi = as4p0(1)*wf(1,j,k)+as4p0(2)*wf(2,j,k) &
                + as4p0(3)*wf(3,j,k)+as4p0(4)*wf(4,j,k)
          dpksi = as4p0(1)*pf(1,j,k)+as4p0(2)*pf(2,j,k) &
                + as4p0(3)*pf(3,j,k)+as4p0(4)*pf(4,j,k)
          drksi = as4p0(1)*rf(1,j,k)+as4p0(2)*rf(2,j,k) &
                + as4p0(3)*rf(3,j,k)+as4p0(4)*rf(4,j,k)

          dueta = a5(1)*( uf(i,j+1,k) - uf(i,j-1,k) ) + &
                  a5(2)*( uf(i,j+2,k) - uf(i,j-2,k) )
          dveta = a5(1)*( vf(i,j+1,k) - vf(i,j-1,k) ) + &
                  a5(2)*( vf(i,j+2,k) - vf(i,j-2,k) )
          dweta = a5(1)*( wf(i,j+1,k) - wf(i,j-1,k) ) + &
                  a5(2)*( wf(i,j+2,k) - wf(i,j-2,k) )
          dpeta = a5(1)*( pf(i,j+1,k) - pf(i,j-1,k) ) + &
                  a5(2)*( pf(i,j+2,k) - pf(i,j-2,k) )
          dreta = a5(1)*( rf(i,j+1,k) - rf(i,j-1,k) ) + &
                  a5(2)*( rf(i,j+2,k) - rf(i,j-2,k) )
          
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
          duksi = as4p1(1)*uf(1,j,k)+as4p1(2)*uf(2,j,k) &
                + as4p1(3)*uf(3,j,k)+as4p1(4)*uf(4,j,k) 
          dvksi = as4p1(1)*vf(1,j,k)+as4p1(2)*vf(2,j,k) &
                + as4p1(3)*vf(3,j,k)+as4p1(4)*vf(4,j,k)
          dwksi = as4p1(1)*wf(1,j,k)+as4p1(2)*wf(2,j,k) &
                + as4p1(3)*wf(3,j,k)+as4p1(4)*wf(4,j,k)
          dpksi = as4p1(1)*pf(1,j,k)+as4p1(2)*pf(2,j,k) &
                + as4p1(3)*pf(3,j,k)+as4p1(4)*pf(4,j,k)
          drksi = as4p1(1)*rf(1,j,k)+as4p1(2)*rf(2,j,k) &
                + as4p1(3)*rf(3,j,k)+as4p1(4)*rf(4,j,k)

          dueta = a5(1)*( uf(i,j+1,k) - uf(i,j-1,k) ) + &
                  a5(2)*( uf(i,j+2,k) - uf(i,j-2,k) )
          dveta = a5(1)*( vf(i,j+1,k) - vf(i,j-1,k) ) + &
                  a5(2)*( vf(i,j+2,k) - vf(i,j-2,k) )
          dweta = a5(1)*( wf(i,j+1,k) - wf(i,j-1,k) ) + &
                  a5(2)*( wf(i,j+2,k) - wf(i,j-2,k) )
          dpeta = a5(1)*( pf(i,j+1,k) - pf(i,j-1,k) ) + &
                  a5(2)*( pf(i,j+2,k) - pf(i,j-2,k) )
          dreta = a5(1)*( rf(i,j+1,k) - rf(i,j-1,k) ) + &
                  a5(2)*( rf(i,j+2,k) - rf(i,j-2,k) )

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
          duksi = as4p2(1)*uf(1,j,k)+as4p2(2)*uf(2,j,k) &
                + as4p2(3)*uf(3,j,k)+as4p2(4)*uf(4,j,k) &
                + as4p2(5)*uf(5,j,k)
          dvksi = as4p2(1)*vf(1,j,k)+as4p2(2)*vf(2,j,k) &
                + as4p2(3)*vf(3,j,k)+as4p2(4)*vf(4,j,k) &
                + as4p2(5)*vf(5,j,k)
          dwksi = as4p2(1)*wf(1,j,k)+as4p2(2)*wf(2,j,k) &
                + as4p2(3)*wf(3,j,k)+as4p2(4)*wf(4,j,k) &
                + as4p2(5)*wf(5,j,k)
          dpksi = as4p2(1)*pf(1,j,k)+as4p2(2)*pf(2,j,k) &
                + as4p2(3)*pf(3,j,k)+as4p2(4)*pf(4,j,k) &
                + as4p2(5)*pf(5,j,k)
          drksi = as4p2(1)*rf(1,j,k)+as4p2(2)*rf(2,j,k) &
                + as4p2(3)*rf(3,j,k)+as4p2(4)*rf(4,j,k) &
                + as4p2(5)*rf(5,j,k)

          dueta = a5(1)*( uf(i,j+1,k) - uf(i,j-1,k) ) + &
                  a5(2)*( uf(i,j+2,k) - uf(i,j-2,k) )
          dveta = a5(1)*( vf(i,j+1,k) - vf(i,j-1,k) ) + &
                  a5(2)*( vf(i,j+2,k) - vf(i,j-2,k) )
          dweta = a5(1)*( wf(i,j+1,k) - wf(i,j-1,k) ) + &
                  a5(2)*( wf(i,j+2,k) - wf(i,j-2,k) )
          dpeta = a5(1)*( pf(i,j+1,k) - pf(i,j-1,k) ) + &
                  a5(2)*( pf(i,j+2,k) - pf(i,j-2,k) )
          dreta = a5(1)*( rf(i,j+1,k) - rf(i,j-1,k) ) + &
                  a5(2)*( rf(i,j+2,k) - rf(i,j-2,k) )

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

    i=4
    do j=ndy,nfy
       do k=1,nz
          duksi = as4p3(1)*uf(1,j,k)+as4p3(2)*uf(2,j,k) &
                + as4p3(3)*uf(3,j,k)+as4p3(4)*uf(4,j,k) &
                + as4p3(5)*uf(5,j,k)+as4p3(6)*uf(6,j,k)
          dvksi = as4p3(1)*vf(1,j,k)+as4p3(2)*vf(2,j,k) &
                + as4p3(3)*vf(3,j,k)+as4p3(4)*vf(4,j,k) &
                + as4p3(5)*vf(5,j,k)+as4p3(6)*vf(6,j,k)
          dwksi = as4p3(1)*wf(1,j,k)+as4p3(2)*wf(2,j,k) &
                + as4p3(3)*wf(3,j,k)+as4p3(4)*wf(4,j,k) &
                + as4p3(5)*wf(5,j,k)+as4p3(6)*wf(6,j,k)
          dpksi = as4p3(1)*pf(1,j,k)+as4p3(2)*pf(2,j,k) &
                + as4p3(3)*pf(3,j,k)+as4p3(4)*pf(4,j,k) &
                + as4p3(5)*pf(5,j,k)+as4p3(6)*pf(6,j,k)
          drksi = as4p3(1)*rf(1,j,k)+as4p3(2)*rf(2,j,k) &
                + as4p3(3)*rf(3,j,k)+as4p3(4)*rf(4,j,k) &
                + as4p3(5)*rf(5,j,k)+as4p3(6)*rf(6,j,k)

          dueta = a5(1)*( uf(i,j+1,k) - uf(i,j-1,k) ) + &
                  a5(2)*( uf(i,j+2,k) - uf(i,j-2,k) )
          dveta = a5(1)*( vf(i,j+1,k) - vf(i,j-1,k) ) + &
                  a5(2)*( vf(i,j+2,k) - vf(i,j-2,k) )
          dweta = a5(1)*( wf(i,j+1,k) - wf(i,j-1,k) ) + &
                  a5(2)*( wf(i,j+2,k) - wf(i,j-2,k) )
          dpeta = a5(1)*( pf(i,j+1,k) - pf(i,j-1,k) ) + &
                  a5(2)*( pf(i,j+2,k) - pf(i,j-2,k) )
          dreta = a5(1)*( rf(i,j+1,k) - rf(i,j-1,k) ) + &
                  a5(2)*( rf(i,j+2,k) - rf(i,j-2,k) )

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

    do i=5,ngh
       do j=ndy,nfy
          do k=1,nz
             duksi = a5(1)* (uf(i+1,j,k)-uf(i-1,j,k)) &
                   + a5(2)* (uf(i+2,j,k)-uf(i-2,j,k))
             dvksi = a5(1)* (vf(i+1,j,k)-vf(i-1,j,k)) &
                   + a5(2)* (vf(i+2,j,k)-vf(i-2,j,k))
             dwksi = a5(1)* (wf(i+1,j,k)-wf(i-1,j,k)) &
                   + a5(2)* (wf(i+2,j,k)-wf(i-2,j,k))
             dpksi = a5(1)* (pf(i+1,j,k)-pf(i-1,j,k)) &
                   + a5(2)* (pf(i+2,j,k)-pf(i-2,j,k))
             drksi = a5(1)* (rf(i+1,j,k)-rf(i-1,j,k)) &
                   + a5(2)* (rf(i+2,j,k)-rf(i-2,j,k))

             dueta = a5(1)*( uf(i,j+1,k) - uf(i,j-1,k) ) + &
                     a5(2)*( uf(i,j+2,k) - uf(i,j-2,k) )
             dveta = a5(1)*( vf(i,j+1,k) - vf(i,j-1,k) ) + &
                     a5(2)*( vf(i,j+2,k) - vf(i,j-2,k) )
             dweta = a5(1)*( wf(i,j+1,k) - wf(i,j-1,k) ) + &
                     a5(2)*( wf(i,j+2,k) - wf(i,j-2,k) )
             dpeta = a5(1)*( pf(i,j+1,k) - pf(i,j-1,k) ) + &
                     a5(2)*( pf(i,j+2,k) - pf(i,j-2,k) )
             dreta = a5(1)*( rf(i,j+1,k) - rf(i,j-1,k) ) + &
                     a5(2)*( rf(i,j+2,k) - rf(i,j-2,k) )

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

  end subroutine bc_TD2d_imin_c_SBP4

  !===============================================================================
  module subroutine bc_TD2d_imax_c_SBP4
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
             vg(l,j,k)= BC_face(1,2)%U0(l,j,k,2)*cosphi(l,j)+BC_face(1,2)%U0(l,j,k,3)*sinphi(l,j) &
                  + sqrt(BC_face(1,2)%U0(l,j,k,6)-BC_face(1,2)%U0(l,j,k,4)**2 &
                  -(BC_face(1,2)%U0(l,j,k,2)*sinphi(l,j)-BC_face(1,2)%U0(l,j,k,3)*cosphi(l,j))**2)
          enddo
       enddo
    enddo

    ! Compute vg*[cos(phi)*dq/dx+sin(phi)*dq/dy+q/r]
    ! ==============================================
    ! (Tam & Webb DRP schemes)
    do i=nxmnghp1,nx-4
       l=i-nxmngh
       do j=ndy,nfy
          do k=1,nz
             duksi = a5(1)* (uf(l+1,j,k)-uf(l-1,j,k)) &
                   + a5(2)* (uf(l+2,j,k)-uf(l-2,j,k))
             dvksi = a5(1)* (vf(l+1,j,k)-vf(l-1,j,k)) &
                   + a5(2)* (vf(l+2,j,k)-vf(l-2,j,k))
             dwksi = a5(1)* (wf(l+1,j,k)-wf(l-1,j,k)) &
                   + a5(2)* (wf(l+2,j,k)-wf(l-2,j,k))
             dpksi = a5(1)* (pf(l+1,j,k)-pf(l-1,j,k)) &
                   + a5(2)* (pf(l+2,j,k)-pf(l-2,j,k))
             drksi = a5(1)* (rf(l+1,j,k)-rf(l-1,j,k)) &
                   + a5(2)* (rf(l+2,j,k)-rf(l-2,j,k))
             
             dueta = a5(1)*( uf(l,j+1,k) - uf(l,j-1,k) ) + &
                     a5(2)*( uf(l,j+2,k) - uf(l,j-2,k) )
             dveta = a5(1)*( vf(l,j+1,k) - vf(l,j-1,k) ) + &
                     a5(2)*( vf(l,j+2,k) - vf(l,j-2,k) )
             dweta = a5(1)*( wf(l,j+1,k) - wf(l,j-1,k) ) + &
                     a5(2)*( wf(l,j+2,k) - wf(l,j-2,k) )
             dpeta = a5(1)*( pf(l,j+1,k) - pf(l,j-1,k) ) + &
                     a5(2)*( pf(l,j+2,k) - pf(l,j-2,k) )
             dreta = a5(1)*( rf(l,j+1,k) - rf(l,j-1,k) ) + &
                     a5(2)*( rf(l,j+2,k) - rf(l,j-2,k) )
          
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
    
    i=nx-3
    l=i-nxmngh
    do j=ndy,nfy
       do k=1,nz
          duksi = as4m3(1)*uf(l+3,j,k)+as4m3(2)*uf(l+2,j,k) &
                + as4m3(3)*uf(l+1,j,k)+as4m3(4)*uf(l  ,j,k) &
                + as4m3(5)*uf(l-1,j,k)+as4m3(6)*uf(l-2,j,k)
          dvksi = as4m3(1)*vf(l+3,j,k)+as4m3(2)*vf(l+2,j,k) &
                + as4m3(3)*vf(l+1,j,k)+as4m3(4)*vf(l  ,j,k) &
                + as4m3(5)*vf(l-1,j,k)+as4m3(6)*vf(l-2,j,k)
          dwksi = as4m3(1)*wf(l+3,j,k)+as4m3(2)*wf(l+2,j,k) &
                + as4m3(3)*wf(l+1,j,k)+as4m3(4)*wf(l  ,j,k) &
                + as4m3(5)*wf(l-1,j,k)+as4m3(6)*wf(l-2,j,k)
          dpksi = as4m3(1)*pf(l+3,j,k)+as4m3(2)*pf(l+2,j,k) &
                + as4m3(3)*pf(l+1,j,k)+as4m3(4)*pf(l  ,j,k) &
                + as4m3(5)*pf(l-1,j,k)+as4m3(6)*pf(l-2,j,k)
          drksi = as4m3(1)*rf(l+3,j,k)+as4m3(2)*rf(l+2,j,k) &
                + as4m3(3)*rf(l+1,j,k)+as4m3(4)*rf(l  ,j,k) &
                + as4m3(5)*rf(l-1,j,k)+as4m3(6)*rf(l-2,j,k)

          dueta = a5(1)*( uf(l,j+1,k) - uf(l,j-1,k) ) + &
                  a5(2)*( uf(l,j+2,k) - uf(l,j-2,k) )
          dveta = a5(1)*( vf(l,j+1,k) - vf(l,j-1,k) ) + &
                  a5(2)*( vf(l,j+2,k) - vf(l,j-2,k) )
          dweta = a5(1)*( wf(l,j+1,k) - wf(l,j-1,k) ) + &
                  a5(2)*( wf(l,j+2,k) - wf(l,j-2,k) )
          dpeta = a5(1)*( pf(l,j+1,k) - pf(l,j-1,k) ) + &
                  a5(2)*( pf(l,j+2,k) - pf(l,j-2,k) )
          dreta = a5(1)*( rf(l,j+1,k) - rf(l,j-1,k) ) + &
                  a5(2)*( rf(l,j+2,k) - rf(l,j-2,k) )
          
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
    
    i=nx-2
    l=i-nxmngh
    do j=ndy,nfy
       do k=1,nz
          duksi = as4m2(1)*uf(l+2,j,k)+as4m2(2)*uf(l+1,j,k) &
                + as4m2(3)*uf(l  ,j,k)+as4m2(4)*uf(l-1,j,k) &
                + as4m2(5)*uf(l-2,j,k)
          dvksi = as4m2(1)*vf(l+2,j,k)+as4m2(2)*vf(l+1,j,k) &
                + as4m2(3)*vf(l  ,j,k)+as4m2(4)*vf(l-1,j,k) &
                + as4m2(5)*vf(l-2,j,k)
          dwksi = as4m2(1)*wf(l+2,j,k)+as4m2(2)*wf(l+1,j,k) &
                + as4m2(3)*wf(l  ,j,k)+as4m2(4)*wf(l-1,j,k) &
                + as4m2(5)*wf(l-2,j,k)
          dpksi = as4m2(1)*pf(l+2,j,k)+as4m2(2)*pf(l+1,j,k) &
                + as4m2(3)*pf(l  ,j,k)+as4m2(4)*pf(l-1,j,k) &
                + as4m2(5)*pf(l-2,j,k)
          drksi = as4m2(1)*rf(l+2,j,k)+as4m2(2)*rf(l+1,j,k) &
                + as4m2(3)*rf(l  ,j,k)+as4m2(4)*rf(l-1,j,k) &
                + as4m2(5)*rf(l-2,j,k)

          dueta = a5(1)*( uf(l,j+1,k) - uf(l,j-1,k) ) + &
                  a5(2)*( uf(l,j+2,k) - uf(l,j-2,k) )
          dveta = a5(1)*( vf(l,j+1,k) - vf(l,j-1,k) ) + &
                  a5(2)*( vf(l,j+2,k) - vf(l,j-2,k) )
          dweta = a5(1)*( wf(l,j+1,k) - wf(l,j-1,k) ) + &
                  a5(2)*( wf(l,j+2,k) - wf(l,j-2,k) )
          dpeta = a5(1)*( pf(l,j+1,k) - pf(l,j-1,k) ) + &
                  a5(2)*( pf(l,j+2,k) - pf(l,j-2,k) )
          dreta = a5(1)*( rf(l,j+1,k) - rf(l,j-1,k) ) + &
                  a5(2)*( rf(l,j+2,k) - rf(l,j-2,k) )
          
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
          duksi = as4m1(1)*uf(l+1,j,k)+as4m1(2)*uf(l  ,j,k) &
                + as4m1(3)*uf(l-1,j,k)+as4m1(4)*uf(l-2,j,k)
          dvksi = as4m1(1)*vf(l+1,j,k)+as4m1(2)*vf(l  ,j,k) &
                + as4m1(3)*vf(l-1,j,k)+as4m1(4)*vf(l-2,j,k)
          dwksi = as4m1(1)*wf(l+1,j,k)+as4m1(2)*wf(l  ,j,k) &
                + as4m1(3)*wf(l-1,j,k)+as4m1(4)*wf(l-2,j,k)
          dpksi = as4m1(1)*pf(l+1,j,k)+as4m1(2)*pf(l  ,j,k) &
                + as4m1(3)*pf(l-1,j,k)+as4m1(4)*pf(l-2,j,k)
          drksi = as4m1(1)*rf(l+1,j,k)+as4m1(2)*rf(l  ,j,k) &
                + as4m1(3)*rf(l-1,j,k)+as4m1(4)*rf(l-2,j,k)

          dueta = a5(1)*( uf(l,j+1,k) - uf(l,j-1,k) ) + &
                  a5(2)*( uf(l,j+2,k) - uf(l,j-2,k) )
          dveta = a5(1)*( vf(l,j+1,k) - vf(l,j-1,k) ) + &
                  a5(2)*( vf(l,j+2,k) - vf(l,j-2,k) )
          dweta = a5(1)*( wf(l,j+1,k) - wf(l,j-1,k) ) + &
                  a5(2)*( wf(l,j+2,k) - wf(l,j-2,k) )
          dpeta = a5(1)*( pf(l,j+1,k) - pf(l,j-1,k) ) + &
                  a5(2)*( pf(l,j+2,k) - pf(l,j-2,k) )
          dreta = a5(1)*( rf(l,j+1,k) - rf(l,j-1,k) ) + &
                  a5(2)*( rf(l,j+2,k) - rf(l,j-2,k) )
          
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
          duksi = as4m0(1)*uf(l  ,j,k)+as4m0(2)*uf(l-1,j,k) &
                + as4m0(3)*uf(l-2,j,k)+as4m0(4)*uf(l-3,j,k) 
          dvksi = as4m0(1)*vf(l  ,j,k)+as4m0(2)*vf(l-1,j,k) &
                + as4m0(3)*vf(l-2,j,k)+as4m0(4)*vf(l-3,j,k)
          dwksi = as4m0(1)*wf(l  ,j,k)+as4m0(2)*wf(l-1,j,k) &
                + as4m0(3)*wf(l-2,j,k)+as4m0(4)*wf(l-3,j,k)
          dpksi = as4m0(1)*pf(l  ,j,k)+as4m0(2)*pf(l-1,j,k) &
                + as4m0(3)*pf(l-2,j,k)+as4m0(4)*pf(l-3,j,k)
          drksi = as4m0(1)*rf(l  ,j,k)+as4m0(2)*rf(l-1,j,k) &
                + as4m0(3)*rf(l-2,j,k)+as4m0(4)*rf(l-3,j,k)
          
          dueta = a5(1)*( uf(l,j+1,k) - uf(l,j-1,k) ) + &
                  a5(2)*( uf(l,j+2,k) - uf(l,j-2,k) )
          dveta = a5(1)*( vf(l,j+1,k) - vf(l,j-1,k) ) + &
                  a5(2)*( vf(l,j+2,k) - vf(l,j-2,k) )
          dweta = a5(1)*( wf(l,j+1,k) - wf(l,j-1,k) ) + &
                  a5(2)*( wf(l,j+2,k) - wf(l,j-2,k) )
          dpeta = a5(1)*( pf(l,j+1,k) - pf(l,j-1,k) ) + &
                  a5(2)*( pf(l,j+2,k) - pf(l,j-2,k) )
          dreta = a5(1)*( rf(l,j+1,k) - rf(l,j-1,k) ) + &
                  a5(2)*( rf(l,j+2,k) - rf(l,j-2,k) )
          
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

  end subroutine bc_TD2d_imax_c_SBP4
  
  !===============================================================================
  module subroutine bc_TD2d_jmin_c_SBP4
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
             vg(i,j,k)= BC_face(2,1)%U0(i,j,k,2)*cosphi(i,j)+BC_face(2,1)%U0(i,j,k,3)*sinphi(i,j) &
                  + sqrt(BC_face(2,1)%U0(i,j,k,6)-BC_face(2,1)%U0(i,j,k,4)**2 &
                  -(BC_face(2,1)%U0(i,j,k,2)*sinphi(i,j)-BC_face(2,1)%U0(i,j,k,3)*cosphi(i,j))**2)
          enddo
       enddo
    enddo

    ! Compute vg*[cos(phi)*dq/dx+sin(phi)*dq/dy+q/r]
    ! ==============================================
    ! (Tam & Webb DRP schemes)
    j=1
    do i=ndx,nfx
       do k=1,nz
          duksi = a5(1)*( uf(i+1,j,k) - uf(i-1,j,k) ) + &
                  a5(2)*( uf(i+2,j,k) - uf(i-2,j,k) )
          dvksi = a5(1)*( vf(i+1,j,k) - vf(i-1,j,k) ) + &
                  a5(2)*( vf(i+2,j,k) - vf(i-2,j,k) )
          dwksi = a5(1)*( wf(i+1,j,k) - wf(i-1,j,k) ) + &
                  a5(2)*( wf(i+2,j,k) - wf(i-2,j,k) )
          dpksi = a5(1)*( pf(i+1,j,k) - pf(i-1,j,k) ) + &
                  a5(2)*( pf(i+2,j,k) - pf(i-2,j,k) )
          drksi = a5(1)*( rf(i+1,j,k) - rf(i-1,j,k) ) + &
                  a5(2)*( rf(i+2,j,k) - rf(i-2,j,k) )
          
          dueta = as4p0(1)*uf(i,1,k)+as4p0(2)*uf(i,2,k) &
                + as4p0(3)*uf(i,3,k)+as4p0(4)*uf(i,4,k)
          dveta = as4p0(1)*vf(i,1,k)+as4p0(2)*vf(i,2,k) &
                + as4p0(3)*vf(i,3,k)+as4p0(4)*vf(i,4,k)
          dweta = as4p0(1)*wf(i,1,k)+as4p0(2)*wf(i,2,k) &
                + as4p0(3)*wf(i,3,k)+as4p0(4)*wf(i,4,k)
          dpeta = as4p0(1)*pf(i,1,k)+as4p0(2)*pf(i,2,k) &
                + as4p0(3)*pf(i,3,k)+as4p0(4)*pf(i,4,k)
          dreta = as4p0(1)*rf(i,1,k)+as4p0(2)*rf(i,2,k) &
                + as4p0(3)*rf(i,3,k)+as4p0(4)*rf(i,4,k)

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
          duksi = a5(1)*( uf(i+1,j,k) - uf(i-1,j,k) ) + &
                  a5(2)*( uf(i+2,j,k) - uf(i-2,j,k) )
          dvksi = a5(1)*( vf(i+1,j,k) - vf(i-1,j,k) ) + &
                  a5(2)*( vf(i+2,j,k) - vf(i-2,j,k) )
          dwksi = a5(1)*( wf(i+1,j,k) - wf(i-1,j,k) ) + &
                  a5(2)*( wf(i+2,j,k) - wf(i-2,j,k) )
          dpksi = a5(1)*( pf(i+1,j,k) - pf(i-1,j,k) ) + &
                  a5(2)*( pf(i+2,j,k) - pf(i-2,j,k) )
          drksi = a5(1)*( rf(i+1,j,k) - rf(i-1,j,k) ) + &
                  a5(2)*( rf(i+2,j,k) - rf(i-2,j,k) )
          
          dueta = as4p1(1)*uf(i,1,k)+as4p1(2)*uf(i,2,k) &
                + as4p1(3)*uf(i,3,k)+as4p1(4)*uf(i,4,k)
          dveta = as4p1(1)*vf(i,1,k)+as4p1(2)*vf(i,2,k) &
                + as4p1(3)*vf(i,3,k)+as4p1(4)*vf(i,4,k)
          dweta = as4p1(1)*wf(i,1,k)+as4p1(2)*wf(i,2,k) &
                + as4p1(3)*wf(i,3,k)+as4p1(4)*wf(i,4,k)
          dpeta = as4p1(1)*pf(i,1,k)+as4p1(2)*pf(i,2,k) &
                + as4p1(3)*pf(i,3,k)+as4p1(4)*pf(i,4,k)
          dreta = as4p1(1)*rf(i,1,k)+as4p1(2)*rf(i,2,k) &
                + as4p1(3)*rf(i,3,k)+as4p1(4)*rf(i,4,k)

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
          duksi = a5(1)*( uf(i+1,j,k) - uf(i-1,j,k) ) + &
                  a5(2)*( uf(i+2,j,k) - uf(i-2,j,k) )
          dvksi = a5(1)*( vf(i+1,j,k) - vf(i-1,j,k) ) + &
                  a5(2)*( vf(i+2,j,k) - vf(i-2,j,k) )
          dwksi = a5(1)*( wf(i+1,j,k) - wf(i-1,j,k) ) + &
                  a5(2)*( wf(i+2,j,k) - wf(i-2,j,k) )
          dpksi = a5(1)*( pf(i+1,j,k) - pf(i-1,j,k) ) + &
                  a5(2)*( pf(i+2,j,k) - pf(i-2,j,k) )
          drksi = a5(1)*( rf(i+1,j,k) - rf(i-1,j,k) ) + &
                  a5(2)*( rf(i+2,j,k) - rf(i-2,j,k) )
          
          dueta = as4p2(1)*uf(i,1,k)+as4p2(2)*uf(i,2,k) &
                + as4p2(3)*uf(i,3,k)+as4p2(4)*uf(i,4,k) &
                + as4p2(5)*uf(i,5,k)
          dveta = as4p2(1)*vf(i,1,k)+as4p2(2)*vf(i,2,k) &
                + as4p2(3)*vf(i,3,k)+as4p2(4)*vf(i,4,k) &
                + as4p2(5)*vf(i,5,k)
          dweta = as4p2(1)*wf(i,1,k)+as4p2(2)*wf(i,2,k) &
                + as4p2(3)*wf(i,3,k)+as4p2(4)*wf(i,4,k) &
                + as4p2(5)*wf(i,5,k)
          dpeta = as4p2(1)*pf(i,1,k)+as4p2(2)*pf(i,2,k) &
                + as4p2(3)*pf(i,3,k)+as4p2(4)*pf(i,4,k) &
                + as4p2(5)*pf(i,5,k)
          dreta = as4p2(1)*rf(i,1,k)+as4p2(2)*rf(i,2,k) &
                + as4p2(3)*rf(i,3,k)+as4p2(4)*rf(i,4,k) &
                + as4p2(5)*rf(i,5,k)

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

    j=4
    do i=ndx,nfx
       do k=1,nz
          duksi = a5(1)*( uf(i+1,j,k) - uf(i-1,j,k) ) + &
                  a5(2)*( uf(i+2,j,k) - uf(i-2,j,k) )
          dvksi = a5(1)*( vf(i+1,j,k) - vf(i-1,j,k) ) + &
                  a5(2)*( vf(i+2,j,k) - vf(i-2,j,k) )
          dwksi = a5(1)*( wf(i+1,j,k) - wf(i-1,j,k) ) + &
                  a5(2)*( wf(i+2,j,k) - wf(i-2,j,k) )
          dpksi = a5(1)*( pf(i+1,j,k) - pf(i-1,j,k) ) + &
                  a5(2)*( pf(i+2,j,k) - pf(i-2,j,k) )
          drksi = a5(1)*( rf(i+1,j,k) - rf(i-1,j,k) ) + &
                  a5(2)*( rf(i+2,j,k) - rf(i-2,j,k) )
          
          dueta = as4p3(1)*uf(i,1,k)+as4p3(2)*uf(i,2,k) &
                + as4p3(3)*uf(i,3,k)+as4p3(4)*uf(i,4,k) &
                + as4p3(5)*uf(i,5,k)+as4p3(6)*uf(i,6,k)
          dveta = as4p3(1)*vf(i,1,k)+as4p3(2)*vf(i,2,k) &
                + as4p3(3)*vf(i,3,k)+as4p3(4)*vf(i,4,k) &
                + as4p3(5)*vf(i,5,k)+as4p3(6)*vf(i,6,k)
          dweta = as4p3(1)*wf(i,1,k)+as4p3(2)*wf(i,2,k) &
                + as4p3(3)*wf(i,3,k)+as4p3(4)*wf(i,4,k) &
                + as4p3(5)*wf(i,5,k)+as4p3(6)*wf(i,6,k)
          dpeta = as4p3(1)*pf(i,1,k)+as4p3(2)*pf(i,2,k) &
                + as4p3(3)*pf(i,3,k)+as4p3(4)*pf(i,4,k) &
                + as4p3(5)*pf(i,5,k)+as4p3(6)*pf(i,6,k)
          dreta = as4p3(1)*rf(i,1,k)+as4p3(2)*rf(i,2,k) &
                + as4p3(3)*rf(i,3,k)+as4p3(4)*rf(i,4,k) &
                + as4p3(5)*rf(i,5,k)+as4p3(6)*rf(i,6,k)

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

    do j=5,ngh
       do i=ndx,nfx
          do k=1,nz
             duksi = a5(1)*( uf(i+1,j,k) - uf(i-1,j,k) ) + &
                     a5(2)*( uf(i+2,j,k) - uf(i-2,j,k) )
             dvksi = a5(1)*( vf(i+1,j,k) - vf(i-1,j,k) ) + &
                     a5(2)*( vf(i+2,j,k) - vf(i-2,j,k) )
             dwksi = a5(1)*( wf(i+1,j,k) - wf(i-1,j,k) ) + &
                     a5(2)*( wf(i+2,j,k) - wf(i-2,j,k) )
             dpksi = a5(1)*( pf(i+1,j,k) - pf(i-1,j,k) ) + &
                     a5(2)*( pf(i+2,j,k) - pf(i-2,j,k) )
             drksi = a5(1)*( rf(i+1,j,k) - rf(i-1,j,k) ) + &
                     a5(2)*( rf(i+2,j,k) - rf(i-2,j,k) )

             dueta = a5(1)*(uf(i,j+1,k)-uf(i,j-1,k)) &
                   + a5(2)*(uf(i,j+2,k)-uf(i,j-2,k))
             dveta = a5(1)*(vf(i,j+1,k)-vf(i,j-1,k)) &
                   + a5(2)*(vf(i,j+2,k)-vf(i,j-2,k))
             dweta = a5(1)*(wf(i,j+1,k)-wf(i,j-1,k)) &
                   + a5(2)*(wf(i,j+2,k)-wf(i,j-2,k))
             dpeta = a5(1)*(pf(i,j+1,k)-pf(i,j-1,k)) &
                   + a5(2)*(pf(i,j+2,k)-pf(i,j-2,k))
             dreta = a5(1)*(rf(i,j+1,k)-rf(i,j-1,k)) &
                   + a5(2)*(rf(i,j+2,k)-rf(i,j-2,k))

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

  end subroutine bc_TD2d_jmin_c_SBP4

  !===============================================================================
  module subroutine bc_TD2d_jmax_c_SBP4
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
             vg(i,l,k)= BC_face(2,2)%U0(i,l,k,2)*cosphi(i,l)+BC_face(2,2)%U0(i,l,k,3)*sinphi(i,l) &
                  + sqrt(BC_face(2,2)%U0(i,l,k,6)-BC_face(2,2)%U0(i,l,k,4)**2 &
                  -(BC_face(2,2)%U0(i,l,k,2)*sinphi(i,l)-BC_face(2,2)%U0(i,l,k,3)*cosphi(i,l))**2)
          enddo
       enddo
    enddo

    ! Compute vg*[cos(phi)*dq/dx+sin(phi)*dq/dy+q/r]
    ! ==============================================
    ! (Tam & Webb DRP schemes)
    do j=nymnghp1,ny-4
       l=j-nymngh
       do i=ndx,nfx
          do k=1,nz
             duksi = a5(1)*( uf(i+1,l,k) - uf(i-1,l,k) ) + &
                     a5(2)*( uf(i+2,l,k) - uf(i-2,l,k) )
             dvksi = a5(1)*( vf(i+1,l,k) - vf(i-1,l,k) ) + &
                     a5(2)*( vf(i+2,l,k) - vf(i-2,l,k) )
             dwksi = a5(1)*( wf(i+1,l,k) - wf(i-1,l,k) ) + &
                     a5(2)*( wf(i+2,l,k) - wf(i-2,l,k) )
             dpksi = a5(1)*( pf(i+1,l,k) - pf(i-1,l,k) ) + &
                     a5(2)*( pf(i+2,l,k) - pf(i-2,l,k) )
             drksi = a5(1)*( rf(i+1,l,k) - rf(i-1,l,k) ) + &
                     a5(2)*( rf(i+2,l,k) - rf(i-2,l,k) )
             
             dueta = a5(1)*(uf(i,l+1,k)-uf(i,l-1,k)) &
                   + a5(2)*(uf(i,l+2,k)-uf(i,l-2,k))
             dveta = a5(1)*(vf(i,l+1,k)-vf(i,l-1,k)) &
                   + a5(2)*(vf(i,l+2,k)-vf(i,l-2,k))
             dweta = a5(1)*(wf(i,l+1,k)-wf(i,l-1,k)) &
                   + a5(2)*(wf(i,l+2,k)-wf(i,l-2,k))
             dpeta = a5(1)*(pf(i,l+1,k)-pf(i,l-1,k)) &
                   + a5(2)*(pf(i,l+2,k)-pf(i,l-2,k))
             dreta = a5(1)*(rf(i,l+1,k)-rf(i,l-1,k)) &
                   + a5(2)*(rf(i,l+2,k)-rf(i,l-2,k))
    
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

    j=ny-3
    l=j-nymngh
    do i=ndx,nfx
       do k=1,nz
          duksi = a5(1)*( uf(i+1,l,k) - uf(i-1,l,k) ) + &
                  a5(2)*( uf(i+2,l,k) - uf(i-2,l,k) )
          dvksi = a5(1)*( vf(i+1,l,k) - vf(i-1,l,k) ) + &
                  a5(2)*( vf(i+2,l,k) - vf(i-2,l,k) )
          dwksi = a5(1)*( wf(i+1,l,k) - wf(i-1,l,k) ) + &
                  a5(2)*( wf(i+2,l,k) - wf(i-2,l,k) )
          dpksi = a5(1)*( pf(i+1,l,k) - pf(i-1,l,k) ) + &
                  a5(2)*( pf(i+2,l,k) - pf(i-2,l,k) )
          drksi = a5(1)*( rf(i+1,l,k) - rf(i-1,l,k) ) + &
                  a5(2)*( rf(i+2,l,k) - rf(i-2,l,k) )
          
          dueta = as4m3(1)*uf(i,l+3,k)+as4m3(2)*uf(i,l+2,k) &
                + as4m3(3)*uf(i,l+1,k)+as4m3(4)*uf(i,l  ,k) &
                + as4m3(5)*uf(i,l-1,k)+as4m3(6)*uf(i,l-2,k)
          dveta = as4m3(1)*vf(i,l+3,k)+as4m3(2)*vf(i,l+2,k) &
                + as4m3(3)*vf(i,l+1,k)+as4m3(4)*vf(i,l  ,k) &
                + as4m3(5)*vf(i,l-1,k)+as4m3(6)*vf(i,l-2,k)
          dweta = as4m3(1)*wf(i,l+3,k)+as4m3(2)*wf(i,l+2,k) &
                + as4m3(3)*wf(i,l+1,k)+as4m3(4)*wf(i,l  ,k) &
                + as4m3(5)*wf(i,l-1,k)+as4m3(6)*wf(i,l-2,k)
          dpeta = as4m3(1)*pf(i,l+3,k)+as4m3(2)*pf(i,l+2,k) &
                + as4m3(3)*pf(i,l+1,k)+as4m3(4)*pf(i,l  ,k) &
                + as4m3(5)*pf(i,l-1,k)+as4m3(6)*pf(i,l-2,k)
          dreta = as4m3(1)*rf(i,l+3,k)+as4m3(2)*rf(i,l+2,k) &
                + as4m3(3)*rf(i,l+1,k)+as4m3(4)*rf(i,l  ,k) &
                + as4m3(5)*rf(i,l-1,k)+as4m3(6)*rf(i,l-2,k)

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
    
    j=ny-2
    l=j-nymngh
    do i=ndx,nfx
       do k=1,nz
          duksi = a5(1)*( uf(i+1,l,k) - uf(i-1,l,k) ) + &
                  a5(2)*( uf(i+2,l,k) - uf(i-2,l,k) )
          dvksi = a5(1)*( vf(i+1,l,k) - vf(i-1,l,k) ) + &
                  a5(2)*( vf(i+2,l,k) - vf(i-2,l,k) )
          dwksi = a5(1)*( wf(i+1,l,k) - wf(i-1,l,k) ) + &
                  a5(2)*( wf(i+2,l,k) - wf(i-2,l,k) )
          dpksi = a5(1)*( pf(i+1,l,k) - pf(i-1,l,k) ) + &
                  a5(2)*( pf(i+2,l,k) - pf(i-2,l,k) )
          drksi = a5(1)*( rf(i+1,l,k) - rf(i-1,l,k) ) + &
                  a5(2)*( rf(i+2,l,k) - rf(i-2,l,k) )
          
          dueta = as4m2(1)*uf(i,l+2,k)+as4m2(2)*uf(i,l+1,k) &
                + as4m2(3)*uf(i,l  ,k)+as4m2(4)*uf(i,l-1,k) &
                + as4m2(5)*uf(i,l-2,k)
          dveta = as4m2(1)*vf(i,l+2,k)+as4m2(2)*vf(i,l+1,k) &
                + as4m2(3)*vf(i,l  ,k)+as4m2(4)*vf(i,l-1,k) &
                + as4m2(5)*vf(i,l-2,k)
          dweta = as4m2(1)*wf(i,l+2,k)+as4m2(2)*wf(i,l+1,k) &
                + as4m2(3)*wf(i,l  ,k)+as4m2(4)*wf(i,l-1,k) &
                + as4m2(5)*wf(i,l-2,k)
          dpeta = as4m2(1)*pf(i,l+2,k)+as4m2(2)*pf(i,l+1,k) &
                + as4m2(3)*pf(i,l  ,k)+as4m2(4)*pf(i,l-1,k) &
                + as4m2(5)*pf(i,l-2,k)
          dreta = as4m2(1)*rf(i,l+2,k)+as4m2(2)*rf(i,l+1,k) &
                + as4m2(3)*rf(i,l  ,k)+as4m2(4)*rf(i,l-1,k) &
                + as4m2(5)*rf(i,l-2,k)

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
          duksi = a5(1)*( uf(i+1,l,k) - uf(i-1,l,k) ) + &
                  a5(2)*( uf(i+2,l,k) - uf(i-2,l,k) )
          dvksi = a5(1)*( vf(i+1,l,k) - vf(i-1,l,k) ) + &
                  a5(2)*( vf(i+2,l,k) - vf(i-2,l,k) )
          dwksi = a5(1)*( wf(i+1,l,k) - wf(i-1,l,k) ) + &
                  a5(2)*( wf(i+2,l,k) - wf(i-2,l,k) )
          dpksi = a5(1)*( pf(i+1,l,k) - pf(i-1,l,k) ) + &
                  a5(2)*( pf(i+2,l,k) - pf(i-2,l,k) )
          drksi = a5(1)*( rf(i+1,l,k) - rf(i-1,l,k) ) + &
                  a5(2)*( rf(i+2,l,k) - rf(i-2,l,k) )
          
          dueta = as4m1(1)*uf(i,l+1,k)+as4m1(2)*uf(i,l  ,k) &
                + as4m1(3)*uf(i,l-1,k)+as4m1(4)*uf(i,l-2,k)
          dveta = as4m1(1)*vf(i,l+1,k)+as4m1(2)*vf(i,l  ,k) &
                + as4m1(3)*vf(i,l-1,k)+as4m1(4)*vf(i,l-2,k)
          dweta = as4m1(1)*wf(i,l+1,k)+as4m1(2)*wf(i,l  ,k) &
                + as4m1(3)*wf(i,l-1,k)+as4m1(4)*wf(i,l-2,k)
          dpeta = as4m1(1)*pf(i,l+1,k)+as4m1(2)*pf(i,l  ,k) &
                + as4m1(3)*pf(i,l-1,k)+as4m1(4)*pf(i,l-2,k)
          dreta = as4m1(1)*rf(i,l+1,k)+as4m1(2)*rf(i,l  ,k) &
                + as4m1(3)*rf(i,l-1,k)+as4m1(4)*rf(i,l-2,k)
          
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
          duksi = a5(1)*( uf(i+1,l,k) - uf(i-1,l,k) ) + &
                  a5(2)*( uf(i+2,l,k) - uf(i-2,l,k) )
          dvksi = a5(1)*( vf(i+1,l,k) - vf(i-1,l,k) ) + &
                  a5(2)*( vf(i+2,l,k) - vf(i-2,l,k) )
          dwksi = a5(1)*( wf(i+1,l,k) - wf(i-1,l,k) ) + &
                  a5(2)*( wf(i+2,l,k) - wf(i-2,l,k) )
          dpksi = a5(1)*( pf(i+1,l,k) - pf(i-1,l,k) ) + &
                  a5(2)*( pf(i+2,l,k) - pf(i-2,l,k) )
          drksi = a5(1)*( rf(i+1,l,k) - rf(i-1,l,k) ) + &
                  a5(2)*( rf(i+2,l,k) - rf(i-2,l,k) )
          
          dueta = as4m0(1)*uf(i,l  ,k)+as4m0(2)*uf(i,l-1,k) &
                + as4m0(3)*uf(i,l-2,k)+as4m0(4)*uf(i,l-3,k)
          dveta = as4m0(1)*vf(i,l  ,k)+as4m0(2)*vf(i,l-1,k) &
                + as4m0(3)*vf(i,l-2,k)+as4m0(4)*vf(i,l-3,k)
          dweta = as4m0(1)*wf(i,l  ,k)+as4m0(2)*wf(i,l-1,k) &
                + as4m0(3)*wf(i,l-2,k)+as4m0(4)*wf(i,l-3,k)
          dpeta = as4m0(1)*pf(i,l  ,k)+as4m0(2)*pf(i,l-1,k) &
                + as4m0(3)*pf(i,l-2,k)+as4m0(4)*pf(i,l-3,k)
          dreta = as4m0(1)*rf(i,l  ,k)+as4m0(2)*rf(i,l-1,k) &
                + as4m0(3)*rf(i,l-2,k)+as4m0(4)*rf(i,l-3,k)
          
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

  end subroutine bc_TD2d_jmax_c_SBP4
  
end submodule smod_TamDong2d_SBP4_faces_c
