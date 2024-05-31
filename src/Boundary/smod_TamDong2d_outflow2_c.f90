!===============================================================
submodule (mod_TamDong2d_c) smod_TamDong2d_outflow2_c
!===============================================================
  !> author: XG
  !> date: February 2020 - modif January 2022
  !> Outflow boundary conditions of Tam and Dong
  !> - 2D (periodic) curvilinear version - outflow in dir.2 (i)
!===============================================================

contains

  !===============================================================================
  module subroutine bc_TD2d_outflow_jmin_c
  !===============================================================================
    !> 2D Tam & Dong's outflow BC: boundary condition at jmin (bottom) - curvilinear version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: cp,av,inn1,ntm1
    real(wp), dimension(nx1:nx2,1:ngh+3,nz1:nz2) :: rf,uf,vf,wf,pf
    real(wp), dimension(nx,ngh,nz) :: rt,ut,vt,wt,pt,vg,c2_
    real(wp) :: duksi,dvksi,dwksi,dpksi,drksi
    real(wp) :: dueta,dveta,dweta,dpeta,dreta
    real(wp) :: duz,dvz,dwz,dpz,drz,dpx,dpy
    !-------------------------------------------------------------------------

    ! Pointers for polar coordinates
    ! ==============================
    ir=>BC_face(2,1)%ir
    cosphi=>BC_face(2,1)%cosphi
    sinphi=>BC_face(2,1)%sinphi

    if (is_2D) then
       duz=0.0_wp
       dvz=0.0_wp
       dwz=0.0_wp
       dpz=0.0_wp
       drz=0.0_wp
    endif

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
       do k=nz1,nz2
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
          do k=nz1,nz2
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

          if (.not.is_2D) then
             duz = ( a7(1)* (uf(i,j,k+1)-uf(i,j,k-1)) &
                   + a7(2)* (uf(i,j,k+2)-uf(i,j,k-2)) + &
                   + a7(3)* (uf(i,j,k+3)-uf(i,j,k-3)) ) *idz(k)
             dvz = ( a7(1)* (vf(i,j,k+1)-vf(i,j,k-1)) &
                   + a7(2)* (vf(i,j,k+2)-vf(i,j,k-2)) + &
                   + a7(3)* (vf(i,j,k+3)-vf(i,j,k-3)) ) *idz(k)
             dwz = ( a7(1)* (wf(i,j,k+1)-wf(i,j,k-1)) &
                   + a7(2)* (wf(i,j,k+2)-wf(i,j,k-2)) + &
                   + a7(3)* (wf(i,j,k+3)-wf(i,j,k-3)) ) *idz(k)
             dpz = ( a7(1)* (pf(i,j,k+1)-pf(i,j,k-1)) &
                   + a7(2)* (pf(i,j,k+2)-pf(i,j,k-2)) + &
                   + a7(3)* (pf(i,j,k+3)-pf(i,j,k-3)) ) *idz(k)
             drz = ( a7(1)* (rf(i,j,k+1)-rf(i,j,k-1)) &
                   + a7(2)* (rf(i,j,k+2)-rf(i,j,k-2)) + &
                   + a7(3)* (rf(i,j,k+3)-rf(i,j,k-3)) ) *idz(k)
          endif

          dpx=(dpksi*y_eta(i,j)-dpeta*y_ksi(i,j))*ijacob(i,j)
          dpy=(dpeta*x_ksi(i,j)-dpksi*x_eta(i,j))*ijacob(i,j)

          pt(i,j,k)=vg(i,j,k)*(dpx*cosphi(i,j)+dpy*sinphi(i,j)+pf(i,j,k)*ir(i,j))

          ! ut=(u0.Grad)uf+dpx/rho0=u0*dux+v0*duy+w0*duz + dpx/rho0
          ut(i,j,k)=(BC_face(2,1)%U0(i,j,k,2)*(duksi*y_eta(i,j)-dueta*y_ksi(i,j)) &
                    +BC_face(2,1)%U0(i,j,k,3)*(dueta*x_ksi(i,j)-duksi*x_eta(i,j)))*ijacob(i,j) &
                    +BC_face(2,1)%U0(i,j,k,4)*duz &
                   + dpx/BC_face(2,1)%U0(i,j,k,1)
          ! vt=(u0.Grad)vf+dpy/rho0=u0*dvx+v0*dvy+w0*dvz + dpy/rho0
          vt(i,j,k)=(BC_face(2,1)%U0(i,j,k,2)*(dvksi*y_eta(i,j)-dveta*y_ksi(i,j)) &
                    +BC_face(2,1)%U0(i,j,k,3)*(dveta*x_ksi(i,j)-dvksi*x_eta(i,j)))*ijacob(i,j) &
                    +BC_face(2,1)%U0(i,j,k,4)*dvz &
                   + dpy/BC_face(2,1)%U0(i,j,k,1)
          ! wt=(u0.Grad)wf+dpz/rho0=u0*dwx+v0*dwy+w0*dwz + dpz/rho0
          wt(i,j,k)=(BC_face(2,1)%U0(i,j,k,2)*(dwksi*y_eta(i,j)-dweta*y_ksi(i,j)) &
                    +BC_face(2,1)%U0(i,j,k,3)*(dreta*x_ksi(i,j)-drksi*x_eta(i,j)))*ijacob(i,j) &
                    +BC_face(2,1)%U0(i,j,k,4)*dwz &
                   + dpz/BC_face(2,1)%U0(i,j,k,1)

          ! rhot=(u0.Grad)rhof-(pt+(u0.Grad)pf)/c0^2
          rt(i,j,k)=(BC_face(2,1)%U0(i,j,k,2)*(drksi*y_eta(i,j)-dreta*y_ksi(i,j)) &
                    +BC_face(2,1)%U0(i,j,k,3)*(dreta*x_ksi(i,j)-drksi*x_eta(i,j)))*ijacob(i,j) &
                    +BC_face(2,1)%U0(i,j,k,4)*drz &
                   -(pt(i,j,k) &
                    +BC_face(2,1)%U0(i,j,k,2)*dpx &
                    +BC_face(2,1)%U0(i,j,k,3)*dpy &
                    +BC_face(2,1)%U0(i,j,k,4)*dpz)/BC_face(2,1)%U0(i,j,k,6)
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

          if (.not.is_2D) then
             duz = ( a7(1)* (uf(i,j,k+1)-uf(i,j,k-1)) &
                   + a7(2)* (uf(i,j,k+2)-uf(i,j,k-2)) + &
                   + a7(3)* (uf(i,j,k+3)-uf(i,j,k-3)) ) *idz(k)
             dvz = ( a7(1)* (vf(i,j,k+1)-vf(i,j,k-1)) &
                   + a7(2)* (vf(i,j,k+2)-vf(i,j,k-2)) + &
                   + a7(3)* (vf(i,j,k+3)-vf(i,j,k-3)) ) *idz(k)
             dwz = ( a7(1)* (wf(i,j,k+1)-wf(i,j,k-1)) &
                   + a7(2)* (wf(i,j,k+2)-wf(i,j,k-2)) + &
                   + a7(3)* (wf(i,j,k+3)-wf(i,j,k-3)) ) *idz(k)
             dpz = ( a7(1)* (pf(i,j,k+1)-pf(i,j,k-1)) &
                   + a7(2)* (pf(i,j,k+2)-pf(i,j,k-2)) + &
                   + a7(3)* (pf(i,j,k+3)-pf(i,j,k-3)) ) *idz(k)
             drz = ( a7(1)* (rf(i,j,k+1)-rf(i,j,k-1)) &
                   + a7(2)* (rf(i,j,k+2)-rf(i,j,k-2)) + &
                   + a7(3)* (rf(i,j,k+3)-rf(i,j,k-3)) ) *idz(k)
          endif

          dpx=(dpksi*y_eta(i,j)-dpeta*y_ksi(i,j))*ijacob(i,j)
          dpy=(dpeta*x_ksi(i,j)-dpksi*x_eta(i,j))*ijacob(i,j)

          pt(i,j,k)=vg(i,j,k)*(dpx*cosphi(i,j)+dpy*sinphi(i,j)+pf(i,j,k)*ir(i,j))

          ! ut=(u0.Grad)uf+dpx/rho0=u0*dux+v0*duy+w0*duz + dpx/rho0
          ut(i,j,k)=(BC_face(2,1)%U0(i,j,k,2)*(duksi*y_eta(i,j)-dueta*y_ksi(i,j)) &
                    +BC_face(2,1)%U0(i,j,k,3)*(dueta*x_ksi(i,j)-duksi*x_eta(i,j)))*ijacob(i,j) &
                    +BC_face(2,1)%U0(i,j,k,4)*duz &
                   + dpx/BC_face(2,1)%U0(i,j,k,1)
          ! vt=(u0.Grad)vf+dpy/rho0=u0*dvx+v0*dvy+w0*dvz + dpy/rho0
          vt(i,j,k)=(BC_face(2,1)%U0(i,j,k,2)*(dvksi*y_eta(i,j)-dveta*y_ksi(i,j)) &
                    +BC_face(2,1)%U0(i,j,k,3)*(dveta*x_ksi(i,j)-dvksi*x_eta(i,j)))*ijacob(i,j) &
                    +BC_face(2,1)%U0(i,j,k,4)*dvz &
                   + dpy/BC_face(2,1)%U0(i,j,k,1)
          ! wt=(u0.Grad)wf+dpz/rho0=u0*dwx+v0*dwy+w0*dwz + dpz/rho0
          wt(i,j,k)=(BC_face(2,1)%U0(i,j,k,2)*(dwksi*y_eta(i,j)-dweta*y_ksi(i,j)) &
                    +BC_face(2,1)%U0(i,j,k,3)*(dreta*x_ksi(i,j)-drksi*x_eta(i,j)))*ijacob(i,j) &
                    +BC_face(2,1)%U0(i,j,k,4)*dwz &
                   + dpz/BC_face(2,1)%U0(i,j,k,1)

          ! rhot=(u0.Grad)rhof-(pt+(u0.Grad)pf)/c0^2
          rt(i,j,k)=(BC_face(2,1)%U0(i,j,k,2)*(drksi*y_eta(i,j)-dreta*y_ksi(i,j)) &
                    +BC_face(2,1)%U0(i,j,k,3)*(dreta*x_ksi(i,j)-drksi*x_eta(i,j)))*ijacob(i,j) &
                    +BC_face(2,1)%U0(i,j,k,4)*drz &
                   -(pt(i,j,k) &
                    +BC_face(2,1)%U0(i,j,k,2)*dpx &
                    +BC_face(2,1)%U0(i,j,k,3)*dpy &
                    +BC_face(2,1)%U0(i,j,k,4)*dpz)/BC_face(2,1)%U0(i,j,k,6)
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

          if (.not.is_2D) then
             duz = ( a7(1)* (uf(i,j,k+1)-uf(i,j,k-1)) &
                   + a7(2)* (uf(i,j,k+2)-uf(i,j,k-2)) + &
                   + a7(3)* (uf(i,j,k+3)-uf(i,j,k-3)) ) *idz(k)
             dvz = ( a7(1)* (vf(i,j,k+1)-vf(i,j,k-1)) &
                   + a7(2)* (vf(i,j,k+2)-vf(i,j,k-2)) + &
                   + a7(3)* (vf(i,j,k+3)-vf(i,j,k-3)) ) *idz(k)
             dwz = ( a7(1)* (wf(i,j,k+1)-wf(i,j,k-1)) &
                   + a7(2)* (wf(i,j,k+2)-wf(i,j,k-2)) + &
                   + a7(3)* (wf(i,j,k+3)-wf(i,j,k-3)) ) *idz(k)
             dpz = ( a7(1)* (pf(i,j,k+1)-pf(i,j,k-1)) &
                   + a7(2)* (pf(i,j,k+2)-pf(i,j,k-2)) + &
                   + a7(3)* (pf(i,j,k+3)-pf(i,j,k-3)) ) *idz(k)
             drz = ( a7(1)* (rf(i,j,k+1)-rf(i,j,k-1)) &
                   + a7(2)* (rf(i,j,k+2)-rf(i,j,k-2)) + &
                   + a7(3)* (rf(i,j,k+3)-rf(i,j,k-3)) ) *idz(k)
          endif

          dpx=(dpksi*y_eta(i,j)-dpeta*y_ksi(i,j))*ijacob(i,j)
          dpy=(dpeta*x_ksi(i,j)-dpksi*x_eta(i,j))*ijacob(i,j)

          pt(i,j,k)=vg(i,j,k)*(dpx*cosphi(i,j)+dpy*sinphi(i,j)+pf(i,j,k)*ir(i,j))

          ! ut=(u0.Grad)uf+dpx/rho0=u0*dux+v0*duy+w0*duz + dpx/rho0
          ut(i,j,k)=(BC_face(2,1)%U0(i,j,k,2)*(duksi*y_eta(i,j)-dueta*y_ksi(i,j)) &
                    +BC_face(2,1)%U0(i,j,k,3)*(dueta*x_ksi(i,j)-duksi*x_eta(i,j)))*ijacob(i,j) &
                    +BC_face(2,1)%U0(i,j,k,4)*duz &
                   + dpx/BC_face(2,1)%U0(i,j,k,1)
          ! vt=(u0.Grad)vf+dpy/rho0=u0*dvx+v0*dvy+w0*dvz + dpy/rho0
          vt(i,j,k)=(BC_face(2,1)%U0(i,j,k,2)*(dvksi*y_eta(i,j)-dveta*y_ksi(i,j)) &
                    +BC_face(2,1)%U0(i,j,k,3)*(dveta*x_ksi(i,j)-dvksi*x_eta(i,j)))*ijacob(i,j) &
                    +BC_face(2,1)%U0(i,j,k,4)*dvz &
                   + dpy/BC_face(2,1)%U0(i,j,k,1)
          ! wt=(u0.Grad)wf+dpz/rho0=u0*dwx+v0*dwy+w0*dwz + dpz/rho0
          wt(i,j,k)=(BC_face(2,1)%U0(i,j,k,2)*(dwksi*y_eta(i,j)-dweta*y_ksi(i,j)) &
                    +BC_face(2,1)%U0(i,j,k,3)*(dreta*x_ksi(i,j)-drksi*x_eta(i,j)))*ijacob(i,j) &
                    +BC_face(2,1)%U0(i,j,k,4)*dwz &
                   + dpz/BC_face(2,1)%U0(i,j,k,1)

          ! rhot=(u0.Grad)rhof-(pt+(u0.Grad)pf)/c0^2
          rt(i,j,k)=(BC_face(2,1)%U0(i,j,k,2)*(drksi*y_eta(i,j)-dreta*y_ksi(i,j)) &
                    +BC_face(2,1)%U0(i,j,k,3)*(dreta*x_ksi(i,j)-drksi*x_eta(i,j)))*ijacob(i,j) &
                    +BC_face(2,1)%U0(i,j,k,4)*drz &
                   -(pt(i,j,k) &
                    +BC_face(2,1)%U0(i,j,k,2)*dpx &
                    +BC_face(2,1)%U0(i,j,k,3)*dpy &
                    +BC_face(2,1)%U0(i,j,k,4)*dpz)/BC_face(2,1)%U0(i,j,k,6)
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

             if (.not.is_2D) then
                duz = ( a7(1)* (uf(i,j,k+1)-uf(i,j,k-1)) &
                      + a7(2)* (uf(i,j,k+2)-uf(i,j,k-2)) + &
                      + a7(3)* (uf(i,j,k+3)-uf(i,j,k-3)) ) *idz(k)
                dvz = ( a7(1)* (vf(i,j,k+1)-vf(i,j,k-1)) &
                      + a7(2)* (vf(i,j,k+2)-vf(i,j,k-2)) + &
                      + a7(3)* (vf(i,j,k+3)-vf(i,j,k-3)) ) *idz(k)
                dwz = ( a7(1)* (wf(i,j,k+1)-wf(i,j,k-1)) &
                      + a7(2)* (wf(i,j,k+2)-wf(i,j,k-2)) + &
                      + a7(3)* (wf(i,j,k+3)-wf(i,j,k-3)) ) *idz(k)
                dpz = ( a7(1)* (pf(i,j,k+1)-pf(i,j,k-1)) &
                      + a7(2)* (pf(i,j,k+2)-pf(i,j,k-2)) + &
                      + a7(3)* (pf(i,j,k+3)-pf(i,j,k-3)) ) *idz(k)
                drz = ( a7(1)* (rf(i,j,k+1)-rf(i,j,k-1)) &
                      + a7(2)* (rf(i,j,k+2)-rf(i,j,k-2)) + &
                      + a7(3)* (rf(i,j,k+3)-rf(i,j,k-3)) ) *idz(k)
             endif

             dpx=(dpksi*y_eta(i,j)-dpeta*y_ksi(i,j))*ijacob(i,j)
             dpy=(dpeta*x_ksi(i,j)-dpksi*x_eta(i,j))*ijacob(i,j)

             pt(i,j,k)=vg(i,j,k)*(dpx*cosphi(i,j)+dpy*sinphi(i,j)+pf(i,j,k)*ir(i,j))

             ! ut=(u0.Grad)uf+dpx/rho0=u0*dux+v0*duy+w0*duz + dpx/rho0
             ut(i,j,k)=(BC_face(2,1)%U0(i,j,k,2)*(duksi*y_eta(i,j)-dueta*y_ksi(i,j)) &
                       +BC_face(2,1)%U0(i,j,k,3)*(dueta*x_ksi(i,j)-duksi*x_eta(i,j)))*ijacob(i,j) &
                       +BC_face(2,1)%U0(i,j,k,4)*duz &
                      + dpx/BC_face(2,1)%U0(i,j,k,1)
             ! vt=(u0.Grad)vf+dpy/rho0=u0*dvx+v0*dvy+w0*dvz + dpy/rho0
             vt(i,j,k)=(BC_face(2,1)%U0(i,j,k,2)*(dvksi*y_eta(i,j)-dveta*y_ksi(i,j)) &
                       +BC_face(2,1)%U0(i,j,k,3)*(dveta*x_ksi(i,j)-dvksi*x_eta(i,j)))*ijacob(i,j) &
                       +BC_face(2,1)%U0(i,j,k,4)*dvz &
                      + dpy/BC_face(2,1)%U0(i,j,k,1)
             ! wt=(u0.Grad)wf+dpz/rho0=u0*dwx+v0*dwy+w0*dwz + dpz/rho0
             wt(i,j,k)=(BC_face(2,1)%U0(i,j,k,2)*(dwksi*y_eta(i,j)-dweta*y_ksi(i,j)) &
                       +BC_face(2,1)%U0(i,j,k,3)*(dreta*x_ksi(i,j)-drksi*x_eta(i,j)))*ijacob(i,j) &
                       +BC_face(2,1)%U0(i,j,k,4)*dwz &
                      + dpz/BC_face(2,1)%U0(i,j,k,1)

             ! rhot=(u0.Grad)rhof-(pt+(u0.Grad)pf)/c0^2
             rt(i,j,k)=(BC_face(2,1)%U0(i,j,k,2)*(drksi*y_eta(i,j)-dreta*y_ksi(i,j)) &
                       +BC_face(2,1)%U0(i,j,k,3)*(dreta*x_ksi(i,j)-drksi*x_eta(i,j)))*ijacob(i,j) &
                       +BC_face(2,1)%U0(i,j,k,4)*drz &
                      -(pt(i,j,k) &
                       +BC_face(2,1)%U0(i,j,k,2)*dpx &
                       +BC_face(2,1)%U0(i,j,k,3)*dpy &
                       +BC_face(2,1)%U0(i,j,k,4)*dpz)/BC_face(2,1)%U0(i,j,k,6)
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

  end subroutine bc_TD2d_outflow_jmin_c

  !===============================================================================
  module subroutine bc_TD2d_outflow_jmax_c
  !===============================================================================
    !> 2D Tam & Dong's outflow BC: boundary condition at jmax (top) - curvilinear version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l
    real(wp) :: cp,av,inn1,ntm1
    real(wp), dimension(nx1:nx2,-2:ngh,nz1:nz2) :: rf,uf,vf,wf,pf
    real(wp), dimension(nx,ngh,nz) :: rt,ut,vt,wt,pt,vg,c2_
    real(wp) :: duksi,dvksi,dwksi,dpksi,drksi
    real(wp) :: dueta,dveta,dweta,dpeta,dreta
    real(wp) :: duz,dvz,dwz,dpz,drz,dpx,dpy
    !-------------------------------------------------------------------------

    ! Pointers for polar coordinates
    ! ==============================
    ir=>BC_face(2,2)%ir
    cosphi=>BC_face(2,2)%cosphi
    sinphi=>BC_face(2,2)%sinphi

    if (is_2D) then
       duz=0.0_wp
       dvz=0.0_wp
       dwz=0.0_wp
       dpz=0.0_wp
       drz=0.0_wp
    endif

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
       do k=nz1,nz2
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
          do k=nz1,nz2
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

             if (.not.is_2D) then
                duz = ( a7(1)* (uf(i,l,k+1)-uf(i,l,k-1)) &
                      + a7(2)* (uf(i,l,k+2)-uf(i,l,k-2)) + &
                      + a7(3)* (uf(i,l,k+3)-uf(i,l,k-3)) ) *idz(k)
                dvz = ( a7(1)* (vf(i,l,k+1)-vf(i,l,k-1)) &
                      + a7(2)* (vf(i,l,k+2)-vf(i,l,k-2)) + &
                      + a7(3)* (vf(i,l,k+3)-vf(i,l,k-3)) ) *idz(k)
                dwz = ( a7(1)* (wf(i,l,k+1)-wf(i,l,k-1)) &
                      + a7(2)* (wf(i,l,k+2)-wf(i,l,k-2)) + &
                      + a7(3)* (wf(i,l,k+3)-wf(i,l,k-3)) ) *idz(k)
                dpz = ( a7(1)* (pf(i,l,k+1)-pf(i,l,k-1)) &
                      + a7(2)* (pf(i,l,k+2)-pf(i,l,k-2)) + &
                      + a7(3)* (pf(i,l,k+3)-pf(i,l,k-3)) ) *idz(k)
                drz = ( a7(1)* (rf(i,l,k+1)-rf(i,l,k-1)) &
                      + a7(2)* (rf(i,l,k+2)-rf(i,l,k-2)) + &
                      + a7(3)* (rf(i,l,k+3)-rf(i,l,k-3)) ) *idz(k)
             endif

             dpx=(dpksi*y_eta(i,j)-dpeta*y_ksi(i,j))*ijacob(i,j)
             dpy=(dpeta*x_ksi(i,j)-dpksi*x_eta(i,j))*ijacob(i,j)

             pt(i,l,k)=vg(i,l,k)*(dpx*cosphi(i,l)+dpy*sinphi(i,l)+pf(i,l,k)*ir(i,l))

             ! ut=(u0.Grad)uf+dpx/rho0=u0*dux+v0*duy+w0*duz + dpx/rho0
             ut(i,l,k)=(BC_face(2,2)%U0(i,l,k,2)*(duksi*y_eta(i,j)-dueta*y_ksi(i,j)) &
                       +BC_face(2,2)%U0(i,l,k,3)*(dueta*x_ksi(i,j)-duksi*x_eta(i,j)))*ijacob(i,j) &
                       +BC_face(2,2)%U0(i,l,k,4)*duz &
                      + dpx/BC_face(2,2)%U0(i,l,k,1)
             ! vt=(u0.Grad)vf+dpy/rho0=u0*dvx+v0*dvy+w0*dvz + dpy/rho0
             vt(i,l,k)=(BC_face(2,2)%U0(i,l,k,2)*(dvksi*y_eta(i,j)-dveta*y_ksi(i,j)) &
                       +BC_face(2,2)%U0(i,l,k,3)*(dveta*x_ksi(i,j)-dvksi*x_eta(i,j)))*ijacob(i,j) &
                       +BC_face(2,2)%U0(i,l,k,4)*dvz &
                      + dpy/BC_face(2,2)%U0(i,l,k,1)
             ! wt=(u0.Grad)wf+dpz/rho0=u0*dwx+v0*dwy+w0*dwz + dpz/rho0
             wt(i,l,k)=(BC_face(2,2)%U0(i,l,k,2)*(dwksi*y_eta(i,j)-dweta*y_ksi(i,j)) &
                       +BC_face(2,2)%U0(i,l,k,3)*(dreta*x_ksi(i,j)-drksi*x_eta(i,j)))*ijacob(i,j) &
                       +BC_face(2,2)%U0(i,l,k,4)*dwz &
                      + dpz/BC_face(2,2)%U0(i,l,k,1)

             ! rhot=(u0.Grad)rhof-(pt+(u0.Grad)pf)/c0^2
             rt(i,l,k)=(BC_face(2,2)%U0(i,l,k,2)*(drksi*y_eta(i,j)-dreta*y_ksi(i,j)) &
                       +BC_face(2,2)%U0(i,l,k,3)*(dreta*x_ksi(i,j)-drksi*x_eta(i,j)))*ijacob(i,j) &
                       +BC_face(2,2)%U0(i,l,k,4)*drz &
                      -(pt(i,l,k) &
                       +BC_face(2,2)%U0(i,l,k,2)*dpx &
                       +BC_face(2,2)%U0(i,l,k,3)*dpy &
                       +BC_face(2,2)%U0(i,l,k,4)*dpz)/BC_face(2,2)%U0(i,l,k,6)
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

          if (.not.is_2D) then
             duz = ( a7(1)* (uf(i,l,k+1)-uf(i,l,k-1)) &
                   + a7(2)* (uf(i,l,k+2)-uf(i,l,k-2)) + &
                   + a7(3)* (uf(i,l,k+3)-uf(i,l,k-3)) ) *idz(k)
             dvz = ( a7(1)* (vf(i,l,k+1)-vf(i,l,k-1)) &
                   + a7(2)* (vf(i,l,k+2)-vf(i,l,k-2)) + &
                   + a7(3)* (vf(i,l,k+3)-vf(i,l,k-3)) ) *idz(k)
             dwz = ( a7(1)* (wf(i,l,k+1)-wf(i,l,k-1)) &
                   + a7(2)* (wf(i,l,k+2)-wf(i,l,k-2)) + &
                   + a7(3)* (wf(i,l,k+3)-wf(i,l,k-3)) ) *idz(k)
             dpz = ( a7(1)* (pf(i,l,k+1)-pf(i,l,k-1)) &
                   + a7(2)* (pf(i,l,k+2)-pf(i,l,k-2)) + &
                   + a7(3)* (pf(i,l,k+3)-pf(i,l,k-3)) ) *idz(k)
             drz = ( a7(1)* (rf(i,l,k+1)-rf(i,l,k-1)) &
                   + a7(2)* (rf(i,l,k+2)-rf(i,l,k-2)) + &
                   + a7(3)* (rf(i,l,k+3)-rf(i,l,k-3)) ) *idz(k)
          endif

          dpx=(dpksi*y_eta(i,j)-dpeta*y_ksi(i,j))*ijacob(i,j)
          dpy=(dpeta*x_ksi(i,j)-dpksi*x_eta(i,j))*ijacob(i,j)

          pt(i,l,k)=vg(i,l,k)*(dpx*cosphi(i,l)+dpy*sinphi(i,l)+pf(i,l,k)*ir(i,l))

          ! ut=(u0.Grad)uf+dpx/rho0=u0*dux+v0*duy+w0*duz + dpx/rho0
          ut(i,l,k)=(BC_face(2,2)%U0(i,l,k,2)*(duksi*y_eta(i,j)-dueta*y_ksi(i,j)) &
                    +BC_face(2,2)%U0(i,l,k,3)*(dueta*x_ksi(i,j)-duksi*x_eta(i,j)))*ijacob(i,j) &
                    +BC_face(2,2)%U0(i,l,k,4)*duz &
                   + dpx/BC_face(2,2)%U0(i,l,k,1)
          ! vt=(u0.Grad)vf+dpy/rho0=u0*dvx+v0*dvy+w0*dvz + dpy/rho0
          vt(i,l,k)=(BC_face(2,2)%U0(i,l,k,2)*(dvksi*y_eta(i,j)-dveta*y_ksi(i,j)) &
                    +BC_face(2,2)%U0(i,l,k,3)*(dveta*x_ksi(i,j)-dvksi*x_eta(i,j)))*ijacob(i,j) &
                    +BC_face(2,2)%U0(i,l,k,4)*dvz &
                   + dpy/BC_face(2,2)%U0(i,l,k,1)
          ! wt=(u0.Grad)wf+dpz/rho0=u0*dwx+v0*dwy+w0*dwz + dpz/rho0
          wt(i,l,k)=(BC_face(2,2)%U0(i,l,k,2)*(dwksi*y_eta(i,j)-dweta*y_ksi(i,j)) &
                    +BC_face(2,2)%U0(i,l,k,3)*(dreta*x_ksi(i,j)-drksi*x_eta(i,j)))*ijacob(i,j) &
                    +BC_face(2,2)%U0(i,l,k,4)*dwz &
                   + dpz/BC_face(2,2)%U0(i,l,k,1)

          ! rhot=(u0.Grad)rhof-(pt+(u0.Grad)pf)/c0^2
          rt(i,l,k)=(BC_face(2,2)%U0(i,l,k,2)*(drksi*y_eta(i,j)-dreta*y_ksi(i,j)) &
                    +BC_face(2,2)%U0(i,l,k,3)*(dreta*x_ksi(i,j)-drksi*x_eta(i,j)))*ijacob(i,j) &
                    +BC_face(2,2)%U0(i,l,k,4)*drz &
                   -(pt(i,l,k) &
                    +BC_face(2,2)%U0(i,l,k,2)*dpx &
                    +BC_face(2,2)%U0(i,l,k,3)*dpy &
                    +BC_face(2,2)%U0(i,l,k,4)*dpz)/BC_face(2,2)%U0(i,l,k,6)
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
          
          if (.not.is_2D) then
             duz = ( a7(1)* (uf(i,l,k+1)-uf(i,l,k-1)) &
                   + a7(2)* (uf(i,l,k+2)-uf(i,l,k-2)) + &
                   + a7(3)* (uf(i,l,k+3)-uf(i,l,k-3)) ) *idz(k)
             dvz = ( a7(1)* (vf(i,l,k+1)-vf(i,l,k-1)) &
                   + a7(2)* (vf(i,l,k+2)-vf(i,l,k-2)) + &
                   + a7(3)* (vf(i,l,k+3)-vf(i,l,k-3)) ) *idz(k)
             dwz = ( a7(1)* (wf(i,l,k+1)-wf(i,l,k-1)) &
                   + a7(2)* (wf(i,l,k+2)-wf(i,l,k-2)) + &
                   + a7(3)* (wf(i,l,k+3)-wf(i,l,k-3)) ) *idz(k)
             dpz = ( a7(1)* (pf(i,l,k+1)-pf(i,l,k-1)) &
                   + a7(2)* (pf(i,l,k+2)-pf(i,l,k-2)) + &
                   + a7(3)* (pf(i,l,k+3)-pf(i,l,k-3)) ) *idz(k)
             drz = ( a7(1)* (rf(i,l,k+1)-rf(i,l,k-1)) &
                   + a7(2)* (rf(i,l,k+2)-rf(i,l,k-2)) + &
                   + a7(3)* (rf(i,l,k+3)-rf(i,l,k-3)) ) *idz(k)
          endif

          dpx=(dpksi*y_eta(i,j)-dpeta*y_ksi(i,j))*ijacob(i,j)
          dpy=(dpeta*x_ksi(i,j)-dpksi*x_eta(i,j))*ijacob(i,j)

          pt(i,l,k)=vg(i,l,k)*(dpx*cosphi(i,l)+dpy*sinphi(i,l)+pf(i,l,k)*ir(i,l))

          ! ut=(u0.Grad)uf+dpx/rho0=u0*dux+v0*duy+w0*duz + dpx/rho0
          ut(i,l,k)=(BC_face(2,2)%U0(i,l,k,2)*(duksi*y_eta(i,j)-dueta*y_ksi(i,j)) &
                    +BC_face(2,2)%U0(i,l,k,3)*(dueta*x_ksi(i,j)-duksi*x_eta(i,j)))*ijacob(i,j) &
                    +BC_face(2,2)%U0(i,l,k,4)*duz &
                   + dpx/BC_face(2,2)%U0(i,l,k,1)
          ! vt=(u0.Grad)vf+dpy/rho0=u0*dvx+v0*dvy+w0*dvz + dpy/rho0
          vt(i,l,k)=(BC_face(2,2)%U0(i,l,k,2)*(dvksi*y_eta(i,j)-dveta*y_ksi(i,j)) &
                    +BC_face(2,2)%U0(i,l,k,3)*(dveta*x_ksi(i,j)-dvksi*x_eta(i,j)))*ijacob(i,j) &
                    +BC_face(2,2)%U0(i,l,k,4)*dvz &
                   + dpy/BC_face(2,2)%U0(i,l,k,1)
          ! wt=(u0.Grad)wf+dpz/rho0=u0*dwx+v0*dwy+w0*dwz + dpz/rho0
          wt(i,l,k)=(BC_face(2,2)%U0(i,l,k,2)*(dwksi*y_eta(i,j)-dweta*y_ksi(i,j)) &
                    +BC_face(2,2)%U0(i,l,k,3)*(dreta*x_ksi(i,j)-drksi*x_eta(i,j)))*ijacob(i,j) &
                    +BC_face(2,2)%U0(i,l,k,4)*dwz &
                   + dpz/BC_face(2,2)%U0(i,l,k,1)

          ! rhot=(u0.Grad)rhof-(pt+(u0.Grad)pf)/c0^2
          rt(i,l,k)=(BC_face(2,2)%U0(i,l,k,2)*(drksi*y_eta(i,j)-dreta*y_ksi(i,j)) &
                    +BC_face(2,2)%U0(i,l,k,3)*(dreta*x_ksi(i,j)-drksi*x_eta(i,j)))*ijacob(i,j) &
                    +BC_face(2,2)%U0(i,l,k,4)*drz &
                   -(pt(i,l,k) &
                    +BC_face(2,2)%U0(i,l,k,2)*dpx &
                    +BC_face(2,2)%U0(i,l,k,3)*dpy &
                    +BC_face(2,2)%U0(i,l,k,4)*dpz)/BC_face(2,2)%U0(i,l,k,6)
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
          
          if (.not.is_2D) then
             duz = ( a7(1)* (uf(i,l,k+1)-uf(i,l,k-1)) &
                   + a7(2)* (uf(i,l,k+2)-uf(i,l,k-2)) + &
                   + a7(3)* (uf(i,l,k+3)-uf(i,l,k-3)) ) *idz(k)
             dvz = ( a7(1)* (vf(i,l,k+1)-vf(i,l,k-1)) &
                   + a7(2)* (vf(i,l,k+2)-vf(i,l,k-2)) + &
                   + a7(3)* (vf(i,l,k+3)-vf(i,l,k-3)) ) *idz(k)
             dwz = ( a7(1)* (wf(i,l,k+1)-wf(i,l,k-1)) &
                   + a7(2)* (wf(i,l,k+2)-wf(i,l,k-2)) + &
                   + a7(3)* (wf(i,l,k+3)-wf(i,l,k-3)) ) *idz(k)
             dpz = ( a7(1)* (pf(i,l,k+1)-pf(i,l,k-1)) &
                   + a7(2)* (pf(i,l,k+2)-pf(i,l,k-2)) + &
                   + a7(3)* (pf(i,l,k+3)-pf(i,l,k-3)) ) *idz(k)
             drz = ( a7(1)* (rf(i,l,k+1)-rf(i,l,k-1)) &
                   + a7(2)* (rf(i,l,k+2)-rf(i,l,k-2)) + &
                   + a7(3)* (rf(i,l,k+3)-rf(i,l,k-3)) ) *idz(k)
          endif

          dpx=(dpksi*y_eta(i,j)-dpeta*y_ksi(i,j))*ijacob(i,j)
          dpy=(dpeta*x_ksi(i,j)-dpksi*x_eta(i,j))*ijacob(i,j)

          pt(i,l,k)=vg(i,l,k)*(dpx*cosphi(i,l)+dpy*sinphi(i,l)+pf(i,l,k)*ir(i,l))

          ! ut=(u0.Grad)uf+dpx/rho0=u0*dux+v0*duy+w0*duz + dpx/rho0
          ut(i,l,k)=(BC_face(2,2)%U0(i,l,k,2)*(duksi*y_eta(i,j)-dueta*y_ksi(i,j)) &
                    +BC_face(2,2)%U0(i,l,k,3)*(dueta*x_ksi(i,j)-duksi*x_eta(i,j)))*ijacob(i,j) &
                    +BC_face(2,2)%U0(i,l,k,4)*duz &
                   + dpx/BC_face(2,2)%U0(i,l,k,1)
          ! vt=(u0.Grad)vf+dpy/rho0=u0*dvx+v0*dvy+w0*dvz + dpy/rho0
          vt(i,l,k)=(BC_face(2,2)%U0(i,l,k,2)*(dvksi*y_eta(i,j)-dveta*y_ksi(i,j)) &
                    +BC_face(2,2)%U0(i,l,k,3)*(dveta*x_ksi(i,j)-dvksi*x_eta(i,j)))*ijacob(i,j) &
                    +BC_face(2,2)%U0(i,l,k,4)*dvz &
                   + dpy/BC_face(2,2)%U0(i,l,k,1)
          ! wt=(u0.Grad)wf+dpz/rho0=u0*dwx+v0*dwy+w0*dwz + dpz/rho0
          wt(i,l,k)=(BC_face(2,2)%U0(i,l,k,2)*(dwksi*y_eta(i,j)-dweta*y_ksi(i,j)) &
                    +BC_face(2,2)%U0(i,l,k,3)*(dreta*x_ksi(i,j)-drksi*x_eta(i,j)))*ijacob(i,j) &
                    +BC_face(2,2)%U0(i,l,k,4)*dwz &
                   + dpz/BC_face(2,2)%U0(i,l,k,1)

          ! rhot=(u0.Grad)rhof-(pt+(u0.Grad)pf)/c0^2
          rt(i,l,k)=(BC_face(2,2)%U0(i,l,k,2)*(drksi*y_eta(i,j)-dreta*y_ksi(i,j)) &
                    +BC_face(2,2)%U0(i,l,k,3)*(dreta*x_ksi(i,j)-drksi*x_eta(i,j)))*ijacob(i,j) &
                    +BC_face(2,2)%U0(i,l,k,4)*drz &
                   -(pt(i,l,k) &
                    +BC_face(2,2)%U0(i,l,k,2)*dpx &
                    +BC_face(2,2)%U0(i,l,k,3)*dpy &
                    +BC_face(2,2)%U0(i,l,k,4)*dpz)/BC_face(2,2)%U0(i,l,k,6)
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

  end subroutine bc_TD2d_outflow_jmax_c

end submodule smod_TamDong2d_outflow2_c
