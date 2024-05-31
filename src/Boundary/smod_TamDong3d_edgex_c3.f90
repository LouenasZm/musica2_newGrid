!===============================================================================
submodule (mod_TamDong3d_c3) smod_TamDong3d_edgex_c3
!===============================================================================
  !> author: XG
  !> date: April 2023
  !> Non-reflecting boundary conditions of Tam and Dong
  !> - 3D curvilinear version - routines for edges along x
!=============================================================================== 

contains

  !===============================================================================
  module subroutine bc_TD3d_jmin_kmin_c3
  !===============================================================================
    !> 3D Tam & Dong's BC: boundary condition at jmin-kmin (edge 2,1,1 /bottom-front)
    !> - 3D curvilinear version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: cp,av,c2_
    real(wp), dimension(nx1:nx2,1:ngh+3,1:ngh+3) :: rf,uf,vf,wf,pf
    real(wp), dimension(nx,1:ngh,1:ngh) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(nx,1:ngh,1:ngh) :: dueta,dveta,dweta,dpeta,dreta
    real(wp), dimension(nx,1:ngh,1:ngh) :: duphi,dvphi,dwphi,dpphi,drphi
    real(wp) :: duksi,dvksi,dwksi,dpksi,drksi
    !-------------------------------------------------------------------------

    ! Pointers for spherical coordinates
    ! ==================================
    ir=>BC_edge(2,1,1)%i_r
    cosphi=>BC_edge(2,1,1)%cosp
    sinphi=>BC_edge(2,1,1)%sinp
    costeta=>BC_edge(2,1,1)%cost
    sinteta=>BC_edge(2,1,1)%sint
    costcosp=>BC_edge(2,1,1)%costcosp
    costsinp=>BC_edge(2,1,1)%costsinp
    sintcosp=>BC_edge(2,1,1)%sintcosp
    sintsinp=>BC_edge(2,1,1)%sintsinp
    
    ! Compute fluctuations
    ! ====================
    do k=1,nghp3
       do j=1,nghp3
          do i=nx1,nx2
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
    ! vg= u0*sin(teta)*cos(phi)+v0*sin(teta)*sin(phi)+w0*cos(teta)
    !   + sqrt{c0^2-[u0*cos(teta)*cos(phi)+v0*cos(teta)*sin(phi)-w0*sin(teta)]^2-[u0*sin(phi)-v0*cos(phi)]^2}
    do k=1,ngh
       do j=1,ngh
          do i=ndx,nfx
             vg(i,j,k)= BC_face(2,1)%U0(i,j,k,2)*sintcosp(i,j,k) +     &
               BC_face(2,1)%U0(i,j,k,3)*sintsinp(i,j,k) +     &
               BC_face(2,1)%U0(i,j,k,4)*costeta(i,j,k)          &
                      + sqrt( BC_face(2,1)%U0(i,j,k,6)-                &
                       ( BC_face(2,1)%U0(i,j,k,2)*costcosp(i,j,k) +    &
                   BC_face(2,1)%U0(i,j,k,3)*costsinp(i,j,k) -    &
                     BC_face(2,1)%U0(i,j,k,4)*sinteta(i,j,k) )**2- &
                  ( BC_face(2,1)%U0(i,j,k,2)*sinphi(i,j,k) -      &
                     BC_face(2,1)%U0(i,j,k,3)*cosphi(i,j,k) )**2 )
          enddo
       enddo
    enddo

    ! Non-centered derivatives in y-direction *sin(teta)*sin(phi)
    ! =======================================
    j=1
    do k=1,ngh
       do i=ndx,nfx
          dueta(i,j,k)= a06(1)*uf(i,1,k)+a06(2)*uf(i,2,k) &
                      + a06(3)*uf(i,3,k)+a06(4)*uf(i,4,k) &
                      + a06(5)*uf(i,5,k)+a06(6)*uf(i,6,k) &
                      + a06(7)*uf(i,7,k)
          dveta(i,j,k)= a06(1)*vf(i,1,k)+a06(2)*vf(i,2,k) &
                      + a06(3)*vf(i,3,k)+a06(4)*vf(i,4,k) &
                      + a06(5)*vf(i,5,k)+a06(6)*vf(i,6,k) &
                      + a06(7)*vf(i,7,k)
          dweta(i,j,k)= a06(1)*wf(i,1,k)+a06(2)*wf(i,2,k) &
                      + a06(3)*wf(i,3,k)+a06(4)*wf(i,4,k) &
                      + a06(5)*wf(i,5,k)+a06(6)*wf(i,6,k) &
                      + a06(7)*wf(i,7,k)
          dpeta(i,j,k)= a06(1)*pf(i,1,k)+a06(2)*pf(i,2,k) &
                      + a06(3)*pf(i,3,k)+a06(4)*pf(i,4,k) &
                      + a06(5)*pf(i,5,k)+a06(6)*pf(i,6,k) &
                      + a06(7)*pf(i,7,k)
          dreta(i,j,k)= a06(1)*rf(i,1,k)+a06(2)*rf(i,2,k) &
                      + a06(3)*rf(i,3,k)+a06(4)*rf(i,4,k) &
                      + a06(5)*rf(i,5,k)+a06(6)*rf(i,6,k) &
                      + a06(7)*rf(i,7,k)
       enddo
    enddo

    j=2
    do k=1,ngh
       do i=ndx,nfx
          dueta(i,j,k)= a15(1)*uf(i,1,k)+a15(2)*uf(i,2,k) &
                      + a15(3)*uf(i,3,k)+a15(4)*uf(i,4,k) &
                      + a15(5)*uf(i,5,k)+a15(6)*uf(i,6,k) &
                      + a15(7)*uf(i,7,k)
          dveta(i,j,k)= a15(1)*vf(i,1,k)+a15(2)*vf(i,2,k) &
                      + a15(3)*vf(i,3,k)+a15(4)*vf(i,4,k) &
                      + a15(5)*vf(i,5,k)+a15(6)*vf(i,6,k) &
                      + a15(7)*vf(i,7,k)
          dweta(i,j,k)= a15(1)*wf(i,1,k)+a15(2)*wf(i,2,k) &
                      + a15(3)*wf(i,3,k)+a15(4)*wf(i,4,k) &
                      + a15(5)*wf(i,5,k)+a15(6)*wf(i,6,k) &
                      + a15(7)*wf(i,7,k)
          dpeta(i,j,k)= a15(1)*pf(i,1,k)+a15(2)*pf(i,2,k) &
                      + a15(3)*pf(i,3,k)+a15(4)*pf(i,4,k) &
                      + a15(5)*pf(i,5,k)+a15(6)*pf(i,6,k) &
                      + a15(7)*pf(i,7,k)
          dreta(i,j,k)= a15(1)*rf(i,1,k)+a15(2)*rf(i,2,k) &
                      + a15(3)*rf(i,3,k)+a15(4)*rf(i,4,k) &
                      + a15(5)*rf(i,5,k)+a15(6)*rf(i,6,k) &
                      + a15(7)*rf(i,7,k)
       enddo
    enddo

    j=3
    do k=1,ngh
       do i=ndx,nfx
          dueta(i,j,k)= a24(1)*uf(i,1,k)+a24(2)*uf(i,2,k) &
                      + a24(3)*uf(i,3,k)+a24(4)*uf(i,4,k) &
                      + a24(5)*uf(i,5,k)+a24(6)*uf(i,6,k) &
                      + a24(7)*uf(i,7,k)
          dveta(i,j,k)= a24(1)*vf(i,1,k)+a24(2)*vf(i,2,k) &
                      + a24(3)*vf(i,3,k)+a24(4)*vf(i,4,k) &
                      + a24(5)*vf(i,5,k)+a24(6)*vf(i,6,k) &
                      + a24(7)*vf(i,7,k)
          dweta(i,j,k)= a24(1)*wf(i,1,k)+a24(2)*wf(i,2,k) &
                      + a24(3)*wf(i,3,k)+a24(4)*wf(i,4,k) &
                      + a24(5)*wf(i,5,k)+a24(6)*wf(i,6,k) &
                      + a24(7)*wf(i,7,k)
          dpeta(i,j,k)= a24(1)*pf(i,1,k)+a24(2)*pf(i,2,k) &
                      + a24(3)*pf(i,3,k)+a24(4)*pf(i,4,k) &
                      + a24(5)*pf(i,5,k)+a24(6)*pf(i,6,k) &
                      + a24(7)*pf(i,7,k)
          dreta(i,j,k)= a24(1)*rf(i,1,k)+a24(2)*rf(i,2,k) &
                      + a24(3)*rf(i,3,k)+a24(4)*rf(i,4,k) &
                      + a24(5)*rf(i,5,k)+a24(6)*rf(i,6,k) &
                      + a24(7)*rf(i,7,k)
       enddo
    enddo

    do j=4,ngh
       do k=1,ngh
          do i=ndx,nfx
             dueta(i,j,k)= a7(1)*(uf(i,j+1,k)-uf(i,j-1,k)) &
                         + a7(2)*(uf(i,j+2,k)-uf(i,j-2,k)) &
                         + a7(3)*(uf(i,j+3,k)-uf(i,j-3,k))
             dveta(i,j,k)= a7(1)*(vf(i,j+1,k)-vf(i,j-1,k)) &
                         + a7(2)*(vf(i,j+2,k)-vf(i,j-2,k)) &
                         + a7(3)*(vf(i,j+3,k)-vf(i,j-3,k))
             dweta(i,j,k)= a7(1)*(wf(i,j+1,k)-wf(i,j-1,k)) &
                         + a7(2)*(wf(i,j+2,k)-wf(i,j-2,k)) &
                         + a7(3)*(wf(i,j+3,k)-wf(i,j-3,k))
             dpeta(i,j,k)= a7(1)*(pf(i,j+1,k)-pf(i,j-1,k)) &
                         + a7(2)*(pf(i,j+2,k)-pf(i,j-2,k)) &
                         + a7(3)*(pf(i,j+3,k)-pf(i,j-3,k))
             dreta(i,j,k)= a7(1)*(rf(i,j+1,k)-rf(i,j-1,k)) &
                         + a7(2)*(rf(i,j+2,k)-rf(i,j-2,k)) &
                         + a7(3)*(rf(i,j+3,k)-rf(i,j-3,k))
          enddo
       enddo
    enddo

    ! Non-centered derivatives in z-direction *cos(teta)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    k=1
    do j=1,ngh
       do i=ndx,nfx
          duphi(i,j,k) = a06(1)*uf(i,j,1)+a06(2)*uf(i,j,2) &
                       + a06(3)*uf(i,j,3)+a06(4)*uf(i,j,4) &
                       + a06(5)*uf(i,j,5)+a06(6)*uf(i,j,6) &
                       + a06(7)*uf(i,j,7)
          dvphi(i,j,k) = a06(1)*vf(i,j,1)+a06(2)*vf(i,j,2) &
                       + a06(3)*vf(i,j,3)+a06(4)*vf(i,j,4) &
                       + a06(5)*vf(i,j,5)+a06(6)*vf(i,j,6) &
                       + a06(7)*vf(i,j,7)
          dwphi(i,j,k) = a06(1)*wf(i,j,1)+a06(2)*wf(i,j,2) &
                       + a06(3)*wf(i,j,3)+a06(4)*wf(i,j,4) &
                       + a06(5)*wf(i,j,5)+a06(6)*wf(i,j,6) &
                       + a06(7)*wf(i,j,7)
          dpphi(i,j,k) = a06(1)*pf(i,j,1)+a06(2)*pf(i,j,2) &
                       + a06(3)*pf(i,j,3)+a06(4)*pf(i,j,4) &
                       + a06(5)*pf(i,j,5)+a06(6)*pf(i,j,6) &
                       + a06(7)*pf(i,j,7)
          drphi(i,j,k) = a06(1)*rf(i,j,1)+a06(2)*rf(i,j,2) &
                       + a06(3)*rf(i,j,3)+a06(4)*rf(i,j,4) &
                       + a06(5)*rf(i,j,5)+a06(6)*rf(i,j,6) &
                       + a06(7)*rf(i,j,7)
       enddo
    enddo

    k=2
    do j=1,ngh
       do i=ndx,nfx
          duphi(i,j,k) = a15(1)*uf(i,j,1)+a15(2)*uf(i,j,2) &
                       + a15(3)*uf(i,j,3)+a15(4)*uf(i,j,4) &
                       + a15(5)*uf(i,j,5)+a15(6)*uf(i,j,6) &
                       + a15(7)*uf(i,j,7)
          dvphi(i,j,k) = a15(1)*vf(i,j,1)+a15(2)*vf(i,j,2) &
                       + a15(3)*vf(i,j,3)+a15(4)*vf(i,j,4) &
                       + a15(5)*vf(i,j,5)+a15(6)*vf(i,j,6) &
                       + a15(7)*vf(i,j,7)
          dwphi(i,j,k) = a15(1)*wf(i,j,1)+a15(2)*wf(i,j,2) &
                       + a15(3)*wf(i,j,3)+a15(4)*wf(i,j,4) &
                       + a15(5)*wf(i,j,5)+a15(6)*wf(i,j,6) &
                       + a15(7)*wf(i,j,7)
          dpphi(i,j,k) = a15(1)*pf(i,j,1)+a15(2)*pf(i,j,2) &
                       + a15(3)*pf(i,j,3)+a15(4)*pf(i,j,4) &
                       + a15(5)*pf(i,j,5)+a15(6)*pf(i,j,6) &
                       + a15(7)*pf(i,j,7)
          drphi(i,j,k) = a15(1)*rf(i,j,1)+a15(2)*rf(i,j,2) &
                       + a15(3)*rf(i,j,3)+a15(4)*rf(i,j,4) &
                       + a15(5)*rf(i,j,5)+a15(6)*rf(i,j,6) &
                       + a15(7)*rf(i,j,7)
       enddo
    enddo

    k=3
    do j=1,ngh
       do i=ndx,nfx
          duphi(i,j,k) = a24(1)*uf(i,j,1)+a24(2)*uf(i,j,2) &
                       + a24(3)*uf(i,j,3)+a24(4)*uf(i,j,4) &
                       + a24(5)*uf(i,j,5)+a24(6)*uf(i,j,6) &
                       + a24(7)*uf(i,j,7)
          dvphi(i,j,k) = a24(1)*vf(i,j,1)+a24(2)*vf(i,j,2) &
                       + a24(3)*vf(i,j,3)+a24(4)*vf(i,j,4) &
                       + a24(5)*vf(i,j,5)+a24(6)*vf(i,j,6) &
                       + a24(7)*vf(i,j,7)
          dwphi(i,j,k) = a24(1)*wf(i,j,1)+a24(2)*wf(i,j,2) &
                       + a24(3)*wf(i,j,3)+a24(4)*wf(i,j,4) &
                       + a24(5)*wf(i,j,5)+a24(6)*wf(i,j,6) &
                       + a24(7)*wf(i,j,7)
          dpphi(i,j,k) = a24(1)*pf(i,j,1)+a24(2)*pf(i,j,2) &
                       + a24(3)*pf(i,j,3)+a24(4)*pf(i,j,4) &
                       + a24(5)*pf(i,j,5)+a24(6)*pf(i,j,6) &
                       + a24(7)*pf(i,j,7)
          drphi(i,j,k) = a24(1)*rf(i,j,1)+a24(2)*rf(i,j,2) &
                       + a24(3)*rf(i,j,3)+a24(4)*rf(i,j,4) &
                       + a24(5)*rf(i,j,5)+a24(6)*rf(i,j,6) &
                       + a24(7)*rf(i,j,7)
       enddo
    enddo

    do k=4,ngh
       do j=1,ngh
          do i=ndx,nfx
             duphi(i,j,k)= a7(1)*(uf(i,j,k+1) - uf(i,j,k-1)) + &
                           a7(2)*(uf(i,j,k+2) - uf(i,j,k-2)) + &
                           a7(3)*(uf(i,j,k+3) - uf(i,j,k-3))
             dvphi(i,j,k)= a7(1)*(vf(i,j,k+1) - vf(i,j,k-1)) + &
                           a7(2)*(vf(i,j,k+2) - vf(i,j,k-2)) + &
                           a7(3)*(vf(i,j,k+3) - vf(i,j,k-3))
             dwphi(i,j,k)= a7(1)*(wf(i,j,k+1) - wf(i,j,k-1)) + &
                           a7(2)*(wf(i,j,k+2) - wf(i,j,k-2)) + &
                           a7(3)*(wf(i,j,k+3) - wf(i,j,k-3))
             dpphi(i,j,k)= a7(1)*(pf(i,j,k+1) - pf(i,j,k-1)) + &
                           a7(2)*(pf(i,j,k+2) - pf(i,j,k-2)) + &
                           a7(3)*(pf(i,j,k+3) - pf(i,j,k-3))
             drphi(i,j,k)= a7(1)*(rf(i,j,k+1) - rf(i,j,k-1)) + &
                           a7(2)*(rf(i,j,k+2) - rf(i,j,k-2)) + &
                           a7(3)*(rf(i,j,k+3) - rf(i,j,k-3))
          enddo
       enddo
    enddo

    ! Compute rt,ut,vt,wt & pt
    ! ========================
    do k=1,ngh
       do j=1,ngh
          do i=ndx,nfx
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

             pt(i,j,k) = vg(i,j,k)*( &
                       ((dpksi*ksi_x(i,j,k)+dpeta(i,j,k)*eta_x(i,j,k)+dpphi(i,j,k)*phi_x(i,j,k))*sintcosp(i,j,k) &
                       +(dpksi*ksi_y(i,j,k)+dpeta(i,j,k)*eta_y(i,j,k)+dpphi(i,j,k)*phi_y(i,j,k))*sintsinp(i,j,k) &
                       +(dpksi*ksi_z(i,j,k)+dpeta(i,j,k)*eta_z(i,j,k)+dpphi(i,j,k)*phi_z(i,j,k))*costeta(i,j,k)  &
                       )*ijacob3(i,j,k) + pf(i,j,k)*ir(i,j,k) )
             ut(i,j,k) = vg(i,j,k)*( &
                       ((duksi*ksi_x(i,j,k)+dueta(i,j,k)*eta_x(i,j,k)+duphi(i,j,k)*phi_x(i,j,k))*sintcosp(i,j,k) &
                       +(duksi*ksi_y(i,j,k)+dueta(i,j,k)*eta_y(i,j,k)+duphi(i,j,k)*phi_y(i,j,k))*sintsinp(i,j,k) &
                       +(duksi*ksi_z(i,j,k)+dueta(i,j,k)*eta_z(i,j,k)+duphi(i,j,k)*phi_z(i,j,k))*costeta(i,j,k)  &
                       )*ijacob3(i,j,k) + uf(i,j,k)*ir(i,j,k) )
             vt(i,j,k) = vg(i,j,k)*( &
                       ((dvksi*ksi_x(i,j,k)+dveta(i,j,k)*eta_x(i,j,k)+dvphi(i,j,k)*phi_x(i,j,k))*sintcosp(i,j,k) &
                       +(dvksi*ksi_y(i,j,k)+dveta(i,j,k)*eta_y(i,j,k)+dvphi(i,j,k)*phi_y(i,j,k))*sintsinp(i,j,k) &
                       +(dvksi*ksi_z(i,j,k)+dveta(i,j,k)*eta_z(i,j,k)+dvphi(i,j,k)*phi_z(i,j,k))*costeta(i,j,k)  &
                       )*ijacob3(i,j,k) + vf(i,j,k)*ir(i,j,k) )
             wt(i,j,k) = vg(i,j,k)*( &
                       ((dwksi*ksi_x(i,j,k)+dweta(i,j,k)*eta_x(i,j,k)+dwphi(i,j,k)*phi_x(i,j,k))*sintcosp(i,j,k) &
                       +(dwksi*ksi_y(i,j,k)+dweta(i,j,k)*eta_y(i,j,k)+dwphi(i,j,k)*phi_y(i,j,k))*sintsinp(i,j,k) &
                       +(dwksi*ksi_z(i,j,k)+dweta(i,j,k)*eta_z(i,j,k)+dwphi(i,j,k)*phi_z(i,j,k))*costeta(i,j,k)  &
                       )*ijacob3(i,j,k) + wf(i,j,k)*ir(i,j,k) )
             rt(i,j,k) = vg(i,j,k)*( &
                       ((drksi*ksi_x(i,j,k)+dreta(i,j,k)*eta_x(i,j,k)+drphi(i,j,k)*phi_x(i,j,k))*sintcosp(i,j,k) &
                       +(drksi*ksi_y(i,j,k)+dreta(i,j,k)*eta_y(i,j,k)+drphi(i,j,k)*phi_y(i,j,k))*sintsinp(i,j,k) &
                       +(drksi*ksi_z(i,j,k)+dreta(i,j,k)*eta_z(i,j,k)+drphi(i,j,k)*phi_z(i,j,k))*costeta(i,j,k)  &
                       )*ijacob3(i,j,k) + rf(i,j,k)*ir(i,j,k) )
          enddo
       enddo
    enddo

    ! Update fluxes at each RK step
    ! =============================
    do i=ndx,nfx
       do j=1,ngh
          do k=1,ngh
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

  end subroutine bc_TD3d_jmin_kmin_c3

  !===============================================================================
  module subroutine bc_TD3d_jmin_kmax_c3
  !===============================================================================
    !> 3D Tam & Dong's BC: boundary condition at jmin-kmax (edge 2,1,2 /bottom-back)
    !> - 3D curvilinear version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l
    real(wp) :: cp,av,c2_
    real(wp), dimension(nx1:nx2,1:ngh+3,-2:ngh) :: rf,uf,vf,wf,pf
    real(wp), dimension(nx,1:ngh,1:ngh) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(nx,1:ngh,1:ngh) :: dueta,dveta,dweta,dpeta,dreta
    real(wp), dimension(nx,1:ngh,1:ngh) :: duphi,dvphi,dwphi,dpphi,drphi
    real(wp) :: duksi,dvksi,dwksi,dpksi,drksi
    !-------------------------------------------------------------------------

    ! Pointers for spherical coordinates
    ! ==================================
    ir=>BC_edge(2,1,2)%i_r
    cosphi=>BC_edge(2,1,2)%cosp
    sinphi=>BC_edge(2,1,2)%sinp
    costeta=>BC_edge(2,1,2)%cost
    sinteta=>BC_edge(2,1,2)%sint
    costcosp=>BC_edge(2,1,2)%costcosp
    costsinp=>BC_edge(2,1,2)%costsinp
    sintcosp=>BC_edge(2,1,2)%sintcosp
    sintsinp=>BC_edge(2,1,2)%sintsinp

    ! Compute fluctuations
    ! ====================
    do j=1,nghp3
       do k=nzmngh-2,nz
          l=k-nzmngh
          do i=nx1,nx2
             rf(i,j,l)=rho_n(i,j,k)-BC_face(2,1)%U0(i,j,k,1)
             uf(i,j,l)=   uu(i,j,k)-BC_face(2,1)%U0(i,j,k,2)
             vf(i,j,l)=   vv(i,j,k)-BC_face(2,1)%U0(i,j,k,3)
             wf(i,j,l)=   ww(i,j,k)-BC_face(2,1)%U0(i,j,k,4)
             pf(i,j,l)=  prs(i,j,k)-BC_face(2,1)%U0(i,j,k,5)
          enddo
       enddo
    enddo

    ! Compute group velocity vg
    ! =========================
    ! vg= u0*sin(teta)*cos(phi)+v0*sin(teta)*sin(phi)+w0*cos(teta)
    !   + sqrt{c0^2-[u0*cos(teta)*cos(phi)+v0*cos(teta)*sin(phi)-w0*sin(teta)]^2-[u0*sin(phi)-v0*cos(phi)]^2}
    do j=1,ngh
       do k=nzmnghp1,nz
          l=k-nzmngh
          do i=ndx,nfx
             vg(i,j,l)= BC_face(2,1)%U0(i,j,k,2)*sintcosp(i,j,l) +     &
               BC_face(2,1)%U0(i,j,k,3)*sintsinp(i,j,l) +     &
               BC_face(2,1)%U0(i,j,k,4)*costeta(i,j,l)          &
                      + sqrt( BC_face(2,1)%U0(i,j,k,6)-                &
                       ( BC_face(2,1)%U0(i,j,k,2)*costcosp(i,j,l) +    &
                   BC_face(2,1)%U0(i,j,k,3)*costsinp(i,j,l) -    &
                     BC_face(2,1)%U0(i,j,k,4)*sinteta(i,j,l) )**2- &
                  ( BC_face(2,1)%U0(i,j,k,2)*sinphi(i,j,l) -      &
                     BC_face(2,1)%U0(i,j,k,3)*cosphi(i,j,l) )**2 )
          enddo
       enddo
    enddo

    ! Non-centered derivatives in y-direction *sin(teta)*sin(phi)
    ! =======================================
    j=1
    do l=1,ngh
       do i=ndx,nfx
          dueta(i,j,l)= a06(1)*uf(i,1,l)+a06(2)*uf(i,2,l) &
                      + a06(3)*uf(i,3,l)+a06(4)*uf(i,4,l) &
                      + a06(5)*uf(i,5,l)+a06(6)*uf(i,6,l) &
                      + a06(7)*uf(i,7,l)
          dveta(i,j,l)= a06(1)*vf(i,1,l)+a06(2)*vf(i,2,l) &
                      + a06(3)*vf(i,3,l)+a06(4)*vf(i,4,l) &
                      + a06(5)*vf(i,5,l)+a06(6)*vf(i,6,l) &
                      + a06(7)*vf(i,7,l)
          dweta(i,j,l)= a06(1)*wf(i,1,l)+a06(2)*wf(i,2,l) &
                      + a06(3)*wf(i,3,l)+a06(4)*wf(i,4,l) &
                      + a06(5)*wf(i,5,l)+a06(6)*wf(i,6,l) &
                      + a06(7)*wf(i,7,l)
          dpeta(i,j,l)= a06(1)*pf(i,1,l)+a06(2)*pf(i,2,l) &
                      + a06(3)*pf(i,3,l)+a06(4)*pf(i,4,l) &
                      + a06(5)*pf(i,5,l)+a06(6)*pf(i,6,l) &
                      + a06(7)*pf(i,7,l)
          dreta(i,j,l)= a06(1)*rf(i,1,l)+a06(2)*rf(i,2,l) &
                      + a06(3)*rf(i,3,l)+a06(4)*rf(i,4,l) &
                      + a06(5)*rf(i,5,l)+a06(6)*rf(i,6,l) &
                      + a06(7)*rf(i,7,l)
       enddo
    enddo

    j=2
    do l=1,ngh
       do i=ndx,nfx
          dueta(i,j,l)= a15(1)*uf(i,1,l)+a15(2)*uf(i,2,l) &
                      + a15(3)*uf(i,3,l)+a15(4)*uf(i,4,l) &
                      + a15(5)*uf(i,5,l)+a15(6)*uf(i,6,l) &
                      + a15(7)*uf(i,7,l)
          dveta(i,j,l)= a15(1)*vf(i,1,l)+a15(2)*vf(i,2,l) &
                      + a15(3)*vf(i,3,l)+a15(4)*vf(i,4,l) &
                      + a15(5)*vf(i,5,l)+a15(6)*vf(i,6,l) &
                      + a15(7)*vf(i,7,l)
          dweta(i,j,l)= a15(1)*wf(i,1,l)+a15(2)*wf(i,2,l) &
                      + a15(3)*wf(i,3,l)+a15(4)*wf(i,4,l) &
                      + a15(5)*wf(i,5,l)+a15(6)*wf(i,6,l) &
                      + a15(7)*wf(i,7,l)
          dpeta(i,j,l)= a15(1)*pf(i,1,l)+a15(2)*pf(i,2,l) &
                      + a15(3)*pf(i,3,l)+a15(4)*pf(i,4,l) &
                      + a15(5)*pf(i,5,l)+a15(6)*pf(i,6,l) &
                      + a15(7)*pf(i,7,l)
          dreta(i,j,l)= a15(1)*rf(i,1,l)+a15(2)*rf(i,2,l) &
                      + a15(3)*rf(i,3,l)+a15(4)*rf(i,4,l) &
                      + a15(5)*rf(i,5,l)+a15(6)*rf(i,6,l) &
                      + a15(7)*rf(i,7,l)
       enddo
    enddo

    j=3
    do l=1,ngh
       do i=ndx,nfx
          dueta(i,j,l)= a24(1)*uf(i,1,l)+a24(2)*uf(i,2,l) &
                      + a24(3)*uf(i,3,l)+a24(4)*uf(i,4,l) &
                      + a24(5)*uf(i,5,l)+a24(6)*uf(i,6,l) &
                      + a24(7)*uf(i,7,l)
          dveta(i,j,l)= a24(1)*vf(i,1,l)+a24(2)*vf(i,2,l) &
                      + a24(3)*vf(i,3,l)+a24(4)*vf(i,4,l) &
                      + a24(5)*vf(i,5,l)+a24(6)*vf(i,6,l) &
                      + a24(7)*vf(i,7,l)
          dweta(i,j,l)= a24(1)*wf(i,1,l)+a24(2)*wf(i,2,l) &
                      + a24(3)*wf(i,3,l)+a24(4)*wf(i,4,l) &
                      + a24(5)*wf(i,5,l)+a24(6)*wf(i,6,l) &
                      + a24(7)*wf(i,7,l)
          dpeta(i,j,l)= a24(1)*pf(i,1,l)+a24(2)*pf(i,2,l) &
                      + a24(3)*pf(i,3,l)+a24(4)*pf(i,4,l) &
                      + a24(5)*pf(i,5,l)+a24(6)*pf(i,6,l) &
                      + a24(7)*pf(i,7,l)
          dreta(i,j,l)= a24(1)*rf(i,1,l)+a24(2)*rf(i,2,l) &
                      + a24(3)*rf(i,3,l)+a24(4)*rf(i,4,l) &
                      + a24(5)*rf(i,5,l)+a24(6)*rf(i,6,l) &
                      + a24(7)*rf(i,7,l)
       enddo
    enddo

    do j=4,ngh
       do l=1,ngh
          do i=ndx,nfx
             dueta(i,j,l)= a7(1)*(uf(i,j+1,l)-uf(i,j-1,l)) &
                         + a7(2)*(uf(i,j+2,l)-uf(i,j-2,l)) &
                         + a7(3)*(uf(i,j+3,l)-uf(i,j-3,l))
             dveta(i,j,l)= a7(1)*(vf(i,j+1,l)-vf(i,j-1,l)) &
                         + a7(2)*(vf(i,j+2,l)-vf(i,j-2,l)) &
                         + a7(3)*(vf(i,j+3,l)-vf(i,j-3,l))
             dweta(i,j,l)= a7(1)*(wf(i,j+1,l)-wf(i,j-1,l)) &
                         + a7(2)*(wf(i,j+2,l)-wf(i,j-2,l)) &
                         + a7(3)*(wf(i,j+3,l)-wf(i,j-3,l))
             dpeta(i,j,l)= a7(1)*(pf(i,j+1,l)-pf(i,j-1,l)) &
                         + a7(2)*(pf(i,j+2,l)-pf(i,j-2,l)) &
                         + a7(3)*(pf(i,j+3,l)-pf(i,j-3,l))
             dreta(i,j,l)= a7(1)*(rf(i,j+1,l)-rf(i,j-1,l)) &
                         + a7(2)*(rf(i,j+2,l)-rf(i,j-2,l)) &
                         + a7(3)*(rf(i,j+3,l)-rf(i,j-3,l))
          enddo
       enddo
    enddo

    ! Non-centered derivatives in z-direction *cos(teta)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    do k=nzmnghp1,nz-3
       l=k-nzmngh
       do j=1,ngh
          do i=ndx,nfx
             duphi(i,j,l)= a7(1)*(uf(i,j,l+1) - uf(i,j,l-1)) + &
                           a7(2)*(uf(i,j,l+2) - uf(i,j,l-2)) + &
                           a7(3)*(uf(i,j,l+3) - uf(i,j,l-3))
             dvphi(i,j,l)= a7(1)*(vf(i,j,l+1) - vf(i,j,l-1)) + &
                           a7(2)*(vf(i,j,l+2) - vf(i,j,l-2)) + &
                           a7(3)*(vf(i,j,l+3) - vf(i,j,l-3))
             dwphi(i,j,l)= a7(1)*(wf(i,j,l+1) - wf(i,j,l-1)) + &
                           a7(2)*(wf(i,j,l+2) - wf(i,j,l-2)) + &
                           a7(3)*(wf(i,j,l+3) - wf(i,j,l-3))
             dpphi(i,j,l)= a7(1)*(pf(i,j,l+1) - pf(i,j,l-1)) + &
                           a7(2)*(pf(i,j,l+2) - pf(i,j,l-2)) + &
                           a7(3)*(pf(i,j,l+3) - pf(i,j,l-3))
             drphi(i,j,l)= a7(1)*(rf(i,j,l+1) - rf(i,j,l-1)) + &
                           a7(2)*(rf(i,j,l+2) - rf(i,j,l-2)) + &
                           a7(3)*(rf(i,j,l+3) - rf(i,j,l-3))
          enddo
       enddo
    enddo

    k=nz-2
    l=k-nzmngh
    do j=1,ngh
       do i=ndx,nfx
          duphi(i,j,l) = a42(1)*uf(i,j,l+2)+a42(2)*uf(i,j,l+1) &
                       + a42(3)*uf(i,j,l  )+a42(4)*uf(i,j,l-1) &
                       + a42(5)*uf(i,j,l-2)+a42(6)*uf(i,j,l-3) &
                       + a42(7)*uf(i,j,l-4)
          dvphi(i,j,l) = a42(1)*vf(i,j,l+2)+a42(2)*vf(i,j,l+1) &
                       + a42(3)*vf(i,j,l  )+a42(4)*vf(i,j,l-1) &
                       + a42(5)*vf(i,j,l-2)+a42(6)*vf(i,j,l-3) &
                       + a42(7)*vf(i,j,l-4)
          dwphi(i,j,l) = a42(1)*wf(i,j,l+2)+a42(2)*wf(i,j,l+1) &
                       + a42(3)*wf(i,j,l  )+a42(4)*wf(i,j,l-1) &
                       + a42(5)*wf(i,j,l-2)+a42(6)*wf(i,j,l-3) &
                       + a42(7)*wf(i,j,l-4)
          dpphi(i,j,l) = a42(1)*pf(i,j,l+2)+a42(2)*pf(i,j,l+1) &
                       + a42(3)*pf(i,j,l  )+a42(4)*pf(i,j,l-1) &
                       + a42(5)*pf(i,j,l-2)+a42(6)*pf(i,j,l-3) &
                       + a42(7)*pf(i,j,l-4)
          drphi(i,j,l) = a42(1)*rf(i,j,l+2)+a42(2)*rf(i,j,l+1) &
                       + a42(3)*rf(i,j,l  )+a42(4)*rf(i,j,l-1) &
                       + a42(5)*rf(i,j,l-2)+a42(6)*rf(i,j,l-3) &
                       + a42(7)*rf(i,j,l-4)
       enddo
    enddo

    k=nz-1
    l=k-nzmngh
    do j=1,ngh
       do i=ndx,nfx
          duphi(i,j,l) = a51(1)*uf(i,j,l+1)+a51(2)*uf(i,j,l  ) &
                       + a51(3)*uf(i,j,l-1)+a51(4)*uf(i,j,l-2) &
                       + a51(5)*uf(i,j,l-3)+a51(6)*uf(i,j,l-4) &
                       + a51(7)*uf(i,j,l-5)
          dvphi(i,j,l) = a51(1)*vf(i,j,l+1)+a51(2)*vf(i,j,l  ) &
                       + a51(3)*vf(i,j,l-1)+a51(4)*vf(i,j,l-2) &
                       + a51(5)*vf(i,j,l-3)+a51(6)*vf(i,j,l-4) &
                       + a51(7)*vf(i,j,l-5)
          dwphi(i,j,l) = a51(1)*wf(i,j,l+1)+a51(2)*wf(i,j,l  ) &
                       + a51(3)*wf(i,j,l-1)+a51(4)*wf(i,j,l-2) &
                       + a51(5)*wf(i,j,l-3)+a51(6)*wf(i,j,l-4) &
                       + a51(7)*wf(i,j,l-5)
          dpphi(i,j,l) = a51(1)*pf(i,j,l+1)+a51(2)*pf(i,j,l  ) &
                       + a51(3)*pf(i,j,l-1)+a51(4)*pf(i,j,l-2) &
                       + a51(5)*pf(i,j,l-3)+a51(6)*pf(i,j,l-4) &
                       + a51(7)*pf(i,j,l-5)
          drphi(i,j,l) = a51(1)*rf(i,j,l+1)+a51(2)*rf(i,j,l  ) &
                       + a51(3)*rf(i,j,l-1)+a51(4)*rf(i,j,l-2) &
                       + a51(5)*rf(i,j,l-3)+a51(6)*rf(i,j,l-4) &
                       + a51(7)*rf(i,j,l-5)
       enddo
    enddo

    k=nz
    l=k-nzmngh
    do j=1,ngh
       do i=ndx,nfx
          duphi(i,j,l) = a60(1)*uf(i,j,l  )+a60(2)*uf(i,j,l-1) &
                       + a60(3)*uf(i,j,l-2)+a60(4)*uf(i,j,l-3) &
                       + a60(5)*uf(i,j,l-4)+a60(6)*uf(i,j,l-5) &
                       + a60(7)*uf(i,j,l-6)
          dvphi(i,j,l) = a60(1)*vf(i,j,l  )+a60(2)*vf(i,j,l-1) &
                       + a60(3)*vf(i,j,l-2)+a60(4)*vf(i,j,l-3) &
                       + a60(5)*vf(i,j,l-4)+a60(6)*vf(i,j,l-5) &
                       + a60(7)*vf(i,j,l-6)
          dwphi(i,j,l) = a60(1)*wf(i,j,l  )+a60(2)*wf(i,j,l-1) &
                       + a60(3)*wf(i,j,l-2)+a60(4)*wf(i,j,l-3) &
                       + a60(5)*wf(i,j,l-4)+a60(6)*wf(i,j,l-5) &
                       + a60(7)*wf(i,j,l-6)
          dpphi(i,j,l) = a60(1)*pf(i,j,l  )+a60(2)*pf(i,j,l-1) &
                       + a60(3)*pf(i,j,l-2)+a60(4)*pf(i,j,l-3) &
                       + a60(5)*pf(i,j,l-4)+a60(6)*pf(i,j,l-5) &
                       + a60(7)*pf(i,j,l-6)
          drphi(i,j,l) = a60(1)*rf(i,j,l  )+a60(2)*rf(i,j,l-1) &
                       + a60(3)*rf(i,j,l-2)+a60(4)*rf(i,j,l-3) &
                       + a60(5)*rf(i,j,l-4)+a60(6)*rf(i,j,l-5) &
                       + a60(7)*rf(i,j,l-6)
       enddo
    enddo

    ! Compute rt,ut,vt,wt & pt
    ! ========================
    do i=ndx,nfx
       do j=1,ngh
          do k=nzmnghp1,nz
             l=k-nzmngh
             duksi = a7(1)*( uf(i+1,j,l) - uf(i-1,j,l) ) + &
                     a7(2)*( uf(i+2,j,l) - uf(i-2,j,l) ) + &
                     a7(3)*( uf(i+3,j,l) - uf(i-3,j,l) )
             dvksi = a7(1)*( vf(i+1,j,l) - vf(i-1,j,l) ) + &
                     a7(2)*( vf(i+2,j,l) - vf(i-2,j,l) ) + &
                     a7(3)*( vf(i+3,j,l) - vf(i-3,j,l) )
             dwksi = a7(1)*( wf(i+1,j,l) - wf(i-1,j,l) ) + &
                     a7(2)*( wf(i+2,j,l) - wf(i-2,j,l) ) + &
                     a7(3)*( wf(i+3,j,l) - wf(i-3,j,l) )
             dpksi = a7(1)*( pf(i+1,j,l) - pf(i-1,j,l) ) + &
                     a7(2)*( pf(i+2,j,l) - pf(i-2,j,l) ) + &
                     a7(3)*( pf(i+3,j,l) - pf(i-3,j,l) )
             drksi = a7(1)*( rf(i+1,j,l) - rf(i-1,j,l) ) + &
                     a7(2)*( rf(i+2,j,l) - rf(i-2,j,l) ) + &
                     a7(3)*( rf(i+3,j,l) - rf(i-3,j,l) )

             pt(i,j,l) = vg(i,j,l)*( &
                       ((dpksi*ksi_x(i,j,k)+dpeta(i,j,l)*eta_x(i,j,k)+dpphi(i,j,l)*phi_x(i,j,k))*sintcosp(i,j,l) &
                       +(dpksi*ksi_y(i,j,k)+dpeta(i,j,l)*eta_y(i,j,k)+dpphi(i,j,l)*phi_y(i,j,k))*sintsinp(i,j,l) &
                       +(dpksi*ksi_z(i,j,k)+dpeta(i,j,l)*eta_z(i,j,k)+dpphi(i,j,l)*phi_z(i,j,k))*costeta(i,j,l)  &
                       )*ijacob3(i,j,k) + pf(i,j,l)*ir(i,j,l) )
             ut(i,j,l) = vg(i,j,l)*( &
                       ((duksi*ksi_x(i,j,k)+dueta(i,j,l)*eta_x(i,j,k)+duphi(i,j,l)*phi_x(i,j,k))*sintcosp(i,j,l) &
                       +(duksi*ksi_y(i,j,k)+dueta(i,j,l)*eta_y(i,j,k)+duphi(i,j,l)*phi_y(i,j,k))*sintsinp(i,j,l) &
                       +(duksi*ksi_z(i,j,k)+dueta(i,j,l)*eta_z(i,j,k)+duphi(i,j,l)*phi_z(i,j,k))*costeta(i,j,l)  &
                       )*ijacob3(i,j,k) + uf(i,j,l)*ir(i,j,l) )
             vt(i,j,l) = vg(i,j,l)*( &
                       ((dvksi*ksi_x(i,j,k)+dveta(i,j,l)*eta_x(i,j,k)+dvphi(i,j,l)*phi_x(i,j,k))*sintcosp(i,j,l) &
                       +(dvksi*ksi_y(i,j,k)+dveta(i,j,l)*eta_y(i,j,k)+dvphi(i,j,l)*phi_y(i,j,k))*sintsinp(i,j,l) &
                       +(dvksi*ksi_z(i,j,k)+dveta(i,j,l)*eta_z(i,j,k)+dvphi(i,j,l)*phi_z(i,j,k))*costeta(i,j,l)  &
                       )*ijacob3(i,j,k) + vf(i,j,l)*ir(i,j,l) )
             wt(i,j,l) = vg(i,j,l)*( &
                       ((dwksi*ksi_x(i,j,k)+dweta(i,j,l)*eta_x(i,j,k)+dwphi(i,j,l)*phi_x(i,j,k))*sintcosp(i,j,l) &
                       +(dwksi*ksi_y(i,j,k)+dweta(i,j,l)*eta_y(i,j,k)+dwphi(i,j,l)*phi_y(i,j,k))*sintsinp(i,j,l) &
                       +(dwksi*ksi_z(i,j,k)+dweta(i,j,l)*eta_z(i,j,k)+dwphi(i,j,l)*phi_z(i,j,k))*costeta(i,j,l)  &
                       )*ijacob3(i,j,k) + wf(i,j,l)*ir(i,j,l) )
             rt(i,j,l) = vg(i,j,l)*( &
                       ((drksi*ksi_x(i,j,k)+dreta(i,j,l)*eta_x(i,j,k)+drphi(i,j,l)*phi_x(i,j,k))*sintcosp(i,j,l) &
                       +(drksi*ksi_y(i,j,k)+dreta(i,j,l)*eta_y(i,j,k)+drphi(i,j,l)*phi_y(i,j,k))*sintsinp(i,j,l) &
                       +(drksi*ksi_z(i,j,k)+dreta(i,j,l)*eta_z(i,j,k)+drphi(i,j,l)*phi_z(i,j,k))*costeta(i,j,l)  &
                       )*ijacob3(i,j,k) + rf(i,j,l)*ir(i,j,l) )
          enddo
       enddo
    enddo

    ! Update fluxes at each RK step
    ! =============================
    do i=ndx,nfx
       do j=1,ngh
          do k=nzmnghp1,nz
             l=k-nzmngh
             cp=cpcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
             av=avcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
             c2_=c2calc_tro(Tmp(i,j,k),rho_n(i,j,k))

             Krho(i,j,k)  = rt(i,j,l)
             Krhou(i,j,k) = uu(i,j,k)*rt(i,j,l)+rho_n(i,j,k)*ut(i,j,l)
             Krhov(i,j,k) = vv(i,j,k)*rt(i,j,l)+rho_n(i,j,k)*vt(i,j,l)
             Krhow(i,j,k) = ww(i,j,k)*rt(i,j,l)+rho_n(i,j,k)*wt(i,j,l)
             Krhoe(i,j,k) = cp/av*(pt(i,j,l)/c2_-rt(i,j,l)) &
                  + (rhoe_n(i,j,k)+prs(i,j,k))/rho_n(i,j,k)*rt(i,j,l) &
                  + rho_n(i,j,k)*(uu(i,j,k)*ut(i,j,l)+vv(i,j,k)*vt(i,j,l)+ww(i,j,k)*wt(i,j,l))
          enddo
       enddo
    enddo

  end subroutine bc_TD3d_jmin_kmax_c3

  !===============================================================================
  module subroutine bc_TD3d_jmax_kmin_c3
  !===============================================================================
    !> 3D Tam & Dong's BC: boundary condition at jmax-kmin (edge 2,2,1 /top-front)
    !> - 3D curvilinear version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l
    real(wp) :: cp,av,c2_
    real(wp), dimension(nx1:nx2,-2:ngh,1:ngh+3) :: rf,uf,vf,wf,pf
    real(wp), dimension(nx,1:ngh,1:ngh) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(nx,1:ngh,1:ngh) :: dueta,dveta,dweta,dpeta,dreta
    real(wp), dimension(nx,1:ngh,1:ngh) :: duphi,dvphi,dwphi,dpphi,drphi
    real(wp) :: duksi,dvksi,dwksi,dpksi,drksi
    !-------------------------------------------------------------------------

    ! Pointers for spherical coordinates
    ! ==================================
    ir=>BC_edge(2,2,1)%i_r
    cosphi=>BC_edge(2,2,1)%cosp
    sinphi=>BC_edge(2,2,1)%sinp
    costeta=>BC_edge(2,2,1)%cost
    sinteta=>BC_edge(2,2,1)%sint
    costcosp=>BC_edge(2,2,1)%costcosp
    costsinp=>BC_edge(2,2,1)%costsinp
    sintcosp=>BC_edge(2,2,1)%sintcosp
    sintsinp=>BC_edge(2,2,1)%sintsinp

    ! Compute fluctuations
    ! ====================
    do j=nymngh-2,ny
       l=j-nymngh
       do k=1,nghp3
          do i=nx1,nx2
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
    ! vg= u0*sin(teta)*cos(phi)+v0*sin(teta)*sin(phi)+w0*cos(teta)
    !   + sqrt{c0^2-[u0*cos(teta)*cos(phi)+v0*cos(teta)*sin(phi)-w0*sin(teta)]^2-[u0*sin(phi)-v0*cos(phi)]^2}
    do l=1,ngh
       do k=1,ngh
          do i=ndx,nfx
             vg(i,l,k)= BC_face(2,2)%U0(i,l,k,2)*sintcosp(i,l,k) +     &
               BC_face(2,2)%U0(i,l,k,3)*sintsinp(i,l,k) +     &
               BC_face(2,2)%U0(i,l,k,4)*costeta(i,l,k)          &
                      + sqrt( BC_face(2,2)%U0(i,l,k,6)-                &
                       ( BC_face(2,2)%U0(i,l,k,2)*costcosp(i,l,k) +    &
                   BC_face(2,2)%U0(i,l,k,3)*costsinp(i,l,k) -    &
                     BC_face(2,2)%U0(i,l,k,4)*sinteta(i,l,k) )**2- &
                  ( BC_face(2,2)%U0(i,l,k,2)*sinphi(i,l,k) -      &
                     BC_face(2,2)%U0(i,l,k,3)*cosphi(i,l,k) )**2 )
          enddo
       enddo
    enddo

    ! Non-centered derivatives in y-direction *sin(teta)*sin(phi)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    do j=nymnghp1,ny-3
       l=j-nymngh
       do k=1,ngh
          do i=ndx,nfx
             dueta(i,l,k) = a7(1)*(uf(i,l+1,k)-uf(i,l-1,k)) &
                          + a7(2)*(uf(i,l+2,k)-uf(i,l-2,k)) &
                          + a7(3)*(uf(i,l+3,k)-uf(i,l-3,k))
             dveta(i,l,k) = a7(1)*(vf(i,l+1,k)-vf(i,l-1,k)) &
                          + a7(2)*(vf(i,l+2,k)-vf(i,l-2,k)) &
                          + a7(3)*(vf(i,l+3,k)-vf(i,l-3,k))
             dweta(i,l,k) = a7(1)*(wf(i,l+1,k)-wf(i,l-1,k)) &
                          + a7(2)*(wf(i,l+2,k)-wf(i,l-2,k)) &
                          + a7(3)*(wf(i,l+3,k)-wf(i,l-3,k))
             dpeta(i,l,k) = a7(1)*(pf(i,l+1,k)-pf(i,l-1,k)) &
                          + a7(2)*(pf(i,l+2,k)-pf(i,l-2,k)) &
                          + a7(3)*(pf(i,l+3,k)-pf(i,l-3,k))
             dreta(i,l,k) = a7(1)*(rf(i,l+1,k)-rf(i,l-1,k)) &
                          + a7(2)*(rf(i,l+2,k)-rf(i,l-2,k)) &
                          + a7(3)*(rf(i,l+3,k)-rf(i,l-3,k))
          enddo
       enddo
    enddo

    j=ny-2
    l=j-nymngh
    do k=1,ngh
       do i=ndx,nfx
          dueta(i,l,k) = a42(1)*uf(i,l+2,k)+a42(2)*uf(i,l+1,k) &
                       + a42(3)*uf(i,l  ,k)+a42(4)*uf(i,l-1,k) &
                       + a42(5)*uf(i,l-2,k)+a42(6)*uf(i,l-3,k) &
                       + a42(7)*uf(i,l-4,k)
          dveta(i,l,k) = a42(1)*vf(i,l+2,k)+a42(2)*vf(i,l+1,k) &
                       + a42(3)*vf(i,l  ,k)+a42(4)*vf(i,l-1,k) &
                       + a42(5)*vf(i,l-2,k)+a42(6)*vf(i,l-3,k) &
                       + a42(7)*vf(i,l-4,k)
          dweta(i,l,k) = a42(1)*wf(i,l+2,k)+a42(2)*wf(i,l+1,k) &
                       + a42(3)*wf(i,l  ,k)+a42(4)*wf(i,l-1,k) &
                       + a42(5)*wf(i,l-2,k)+a42(6)*wf(i,l-3,k) &
                       + a42(7)*wf(i,l-4,k)
          dpeta(i,l,k) = a42(1)*pf(i,l+2,k)+a42(2)*pf(i,l+1,k) &
                       + a42(3)*pf(i,l  ,k)+a42(4)*pf(i,l-1,k) &
                       + a42(5)*pf(i,l-2,k)+a42(6)*pf(i,l-3,k) &
                       + a42(7)*pf(i,l-4,k)
          dreta(i,l,k) = a42(1)*rf(i,l+2,k)+a42(2)*rf(i,l+1,k) &
                       + a42(3)*rf(i,l  ,k)+a42(4)*rf(i,l-1,k) &
                       + a42(5)*rf(i,l-2,k)+a42(6)*rf(i,l-3,k) &
                       + a42(7)*rf(i,l-4,k)
       enddo
    enddo

    j=ny-1
    l=j-nymngh
    do k=1,ngh
       do i=ndx,nfx
          dueta(i,l,k) = a51(1)*uf(i,l+1,k)+a51(2)*uf(i,l  ,k) &
                       + a51(3)*uf(i,l-1,k)+a51(4)*uf(i,l-2,k) &
                       + a51(5)*uf(i,l-3,k)+a51(6)*uf(i,l-4,k) &
                       + a51(7)*uf(i,l-5,k)
          dveta(i,l,k) = a51(1)*vf(i,l+1,k)+a51(2)*vf(i,l  ,k) &
                       + a51(3)*vf(i,l-1,k)+a51(4)*vf(i,l-2,k) &
                       + a51(5)*vf(i,l-3,k)+a51(6)*vf(i,l-4,k) &
                       + a51(7)*vf(i,l-5,k)
          dweta(i,l,k) = a51(1)*wf(i,l+1,k)+a51(2)*wf(i,l  ,k) &
                       + a51(3)*wf(i,l-1,k)+a51(4)*wf(i,l-2,k) &
                       + a51(5)*wf(i,l-3,k)+a51(6)*wf(i,l-4,k) &
                       + a51(7)*wf(i,l-5,k)
          dpeta(i,l,k) = a51(1)*pf(i,l+1,k)+a51(2)*pf(i,l  ,k) &
                       + a51(3)*pf(i,l-1,k)+a51(4)*pf(i,l-2,k) &
                       + a51(5)*pf(i,l-3,k)+a51(6)*pf(i,l-4,k) &
                       + a51(7)*pf(i,l-5,k)
          dreta(i,l,k) = a51(1)*rf(i,l+1,k)+a51(2)*rf(i,l  ,k) &
                       + a51(3)*rf(i,l-1,k)+a51(4)*rf(i,l-2,k) &
                       + a51(5)*rf(i,l-3,k)+a51(6)*rf(i,l-4,k) &
                       + a51(7)*rf(i,l-5,k)
       enddo
    enddo

    j=ny
    l=j-nymngh
    do k=1,ngh
       do i=ndx,nfx
          dueta(i,l,k) = a60(1)*uf(i,l  ,k)+a60(2)*uf(i,l-1,k) &
                       + a60(3)*uf(i,l-2,k)+a60(4)*uf(i,l-3,k) &
                       + a60(5)*uf(i,l-4,k)+a60(6)*uf(i,l-5,k) &
                       + a60(7)*uf(i,l-6,k)
          dveta(i,l,k) = a60(1)*vf(i,l  ,k)+a60(2)*vf(i,l-1,k) &
                       + a60(3)*vf(i,l-2,k)+a60(4)*vf(i,l-3,k) &
                       + a60(5)*vf(i,l-4,k)+a60(6)*vf(i,l-5,k) &
                       + a60(7)*vf(i,l-6,k)
          dweta(i,l,k) = a60(1)*wf(i,l  ,k)+a60(2)*wf(i,l-1,k) &
                       + a60(3)*wf(i,l-2,k)+a60(4)*wf(i,l-3,k) &
                       + a60(5)*wf(i,l-4,k)+a60(6)*wf(i,l-5,k) &
                       + a60(7)*wf(i,l-6,k)
          dpeta(i,l,k) = a60(1)*pf(i,l  ,k)+a60(2)*pf(i,l-1,k) &
                       + a60(3)*pf(i,l-2,k)+a60(4)*pf(i,l-3,k) &
                       + a60(5)*pf(i,l-4,k)+a60(6)*pf(i,l-5,k) &
                       + a60(7)*pf(i,l-6,k)
          dreta(i,l,k) = a60(1)*rf(i,l  ,k)+a60(2)*rf(i,l-1,k) &
                       + a60(3)*rf(i,l-2,k)+a60(4)*rf(i,l-3,k) &
                       + a60(5)*rf(i,l-4,k)+a60(6)*rf(i,l-5,k) &
                       + a60(7)*rf(i,l-6,k)
       enddo
    enddo

    ! Non-centered derivatives in z-direction *cos(teta)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    k=1
    do l=1,ngh
       do i=ndx,nfx
          duphi(i,l,k) = a06(1)*uf(i,l,1)+a06(2)*uf(i,l,2) &
                       + a06(3)*uf(i,l,3)+a06(4)*uf(i,l,4) &
                       + a06(5)*uf(i,l,5)+a06(6)*uf(i,l,6) &
                       + a06(7)*uf(i,l,7)
          dvphi(i,l,k) = a06(1)*vf(i,l,1)+a06(2)*vf(i,l,2) &
                       + a06(3)*vf(i,l,3)+a06(4)*vf(i,l,4) &
                       + a06(5)*vf(i,l,5)+a06(6)*vf(i,l,6) &
                       + a06(7)*vf(i,l,7)
          dwphi(i,l,k) = a06(1)*wf(i,l,1)+a06(2)*wf(i,l,2) &
                       + a06(3)*wf(i,l,3)+a06(4)*wf(i,l,4) &
                       + a06(5)*wf(i,l,5)+a06(6)*wf(i,l,6) &
                       + a06(7)*wf(i,l,7)
          dpphi(i,l,k) = a06(1)*pf(i,l,1)+a06(2)*pf(i,l,2) &
                       + a06(3)*pf(i,l,3)+a06(4)*pf(i,l,4) &
                       + a06(5)*pf(i,l,5)+a06(6)*pf(i,l,6) &
                       + a06(7)*pf(i,l,7)
          drphi(i,l,k) = a06(1)*rf(i,l,1)+a06(2)*rf(i,l,2) &
                       + a06(3)*rf(i,l,3)+a06(4)*rf(i,l,4) &
                       + a06(5)*rf(i,l,5)+a06(6)*rf(i,l,6) &
                       + a06(7)*rf(i,l,7)
       enddo
    enddo

    k=2
    do l=1,ngh
       do i=ndx,nfx
          duphi(i,l,k) = a15(1)*uf(i,l,1)+a15(2)*uf(i,l,2) &
                       + a15(3)*uf(i,l,3)+a15(4)*uf(i,l,4) &
                       + a15(5)*uf(i,l,5)+a15(6)*uf(i,l,6) &
                       + a15(7)*uf(i,l,7)
          dvphi(i,l,k) = a15(1)*vf(i,l,1)+a15(2)*vf(i,l,2) &
                       + a15(3)*vf(i,l,3)+a15(4)*vf(i,l,4) &
                       + a15(5)*vf(i,l,5)+a15(6)*vf(i,l,6) &
                       + a15(7)*vf(i,l,7)
          dwphi(i,l,k) = a15(1)*wf(i,l,1)+a15(2)*wf(i,l,2) &
                       + a15(3)*wf(i,l,3)+a15(4)*wf(i,l,4) &
                       + a15(5)*wf(i,l,5)+a15(6)*wf(i,l,6) &
                       + a15(7)*wf(i,l,7)
          dpphi(i,l,k) = a15(1)*pf(i,l,1)+a15(2)*pf(i,l,2) &
                       + a15(3)*pf(i,l,3)+a15(4)*pf(i,l,4) &
                       + a15(5)*pf(i,l,5)+a15(6)*pf(i,l,6) &
                       + a15(7)*pf(i,l,7)
          drphi(i,l,k) = a15(1)*rf(i,l,1)+a15(2)*rf(i,l,2) &
                       + a15(3)*rf(i,l,3)+a15(4)*rf(i,l,4) &
                       + a15(5)*rf(i,l,5)+a15(6)*rf(i,l,6) &
                       + a15(7)*rf(i,l,7)
       enddo
    enddo

    k=3
    do l=1,ngh
       do i=ndx,nfx
          duphi(i,l,k) = a24(1)*uf(i,l,1)+a24(2)*uf(i,l,2) &
                       + a24(3)*uf(i,l,3)+a24(4)*uf(i,l,4) &
                       + a24(5)*uf(i,l,5)+a24(6)*uf(i,l,6) &
                       + a24(7)*uf(i,l,7)
          dvphi(i,l,k) = a24(1)*vf(i,l,1)+a24(2)*vf(i,l,2) &
                       + a24(3)*vf(i,l,3)+a24(4)*vf(i,l,4) &
                       + a24(5)*vf(i,l,5)+a24(6)*vf(i,l,6) &
                       + a24(7)*vf(i,l,7)
          dwphi(i,l,k) = a24(1)*wf(i,l,1)+a24(2)*wf(i,l,2) &
                       + a24(3)*wf(i,l,3)+a24(4)*wf(i,l,4) &
                       + a24(5)*wf(i,l,5)+a24(6)*wf(i,l,6) &
                       + a24(7)*wf(i,l,7)
          dpphi(i,l,k) = a24(1)*pf(i,l,1)+a24(2)*pf(i,l,2) &
                       + a24(3)*pf(i,l,3)+a24(4)*pf(i,l,4) &
                       + a24(5)*pf(i,l,5)+a24(6)*pf(i,l,6) &
                       + a24(7)*pf(i,l,7)
          drphi(i,l,k) = a24(1)*rf(i,l,1)+a24(2)*rf(i,l,2) &
                       + a24(3)*rf(i,l,3)+a24(4)*rf(i,l,4) &
                       + a24(5)*rf(i,l,5)+a24(6)*rf(i,l,6) &
                       + a24(7)*rf(i,l,7)
       enddo
    enddo

    do k=4,ngh
       do l=1,ngh
          do i=ndx,nfx
             duphi(i,l,k)= a7(1)*(uf(i,l,k+1) - uf(i,l,k-1)) + &
                           a7(2)*(uf(i,l,k+2) - uf(i,l,k-2)) + &
                           a7(3)*(uf(i,l,k+3) - uf(i,l,k-3))
             dvphi(i,l,k)= a7(1)*(vf(i,l,k+1) - vf(i,l,k-1)) + &
                           a7(2)*(vf(i,l,k+2) - vf(i,l,k-2)) + &
                           a7(3)*(vf(i,l,k+3) - vf(i,l,k-3))
             dwphi(i,l,k)= a7(1)*(wf(i,l,k+1) - wf(i,l,k-1)) + &
                           a7(2)*(wf(i,l,k+2) - wf(i,l,k-2)) + &
                           a7(3)*(wf(i,l,k+3) - wf(i,l,k-3))
             dpphi(i,l,k)= a7(1)*(pf(i,l,k+1) - pf(i,l,k-1)) + &
                           a7(2)*(pf(i,l,k+2) - pf(i,l,k-2)) + &
                           a7(3)*(pf(i,l,k+3) - pf(i,l,k-3))
             drphi(i,l,k)= a7(1)*(rf(i,l,k+1) - rf(i,l,k-1)) + &
                           a7(2)*(rf(i,l,k+2) - rf(i,l,k-2)) + &
                           a7(3)*(rf(i,l,k+3) - rf(i,l,k-3))
          enddo
       enddo
    enddo

    ! Compute rt,ut,vt,wt & pt
    ! ========================
    do i=ndx,nfx
       do j=nymnghp1,ny
          l=j-nymngh
          do k=1,ngh
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

             pt(i,l,k) = vg(i,l,k)*( &
                       ((dpksi*ksi_x(i,j,k)+dpeta(i,l,k)*eta_x(i,j,k)+dpphi(i,l,k)*phi_x(i,j,k))*sintcosp(i,l,k) &
                       +(dpksi*ksi_y(i,j,k)+dpeta(i,l,k)*eta_y(i,j,k)+dpphi(i,l,k)*phi_y(i,j,k))*sintsinp(i,l,k) &
                       +(dpksi*ksi_z(i,j,k)+dpeta(i,l,k)*eta_z(i,j,k)+dpphi(i,l,k)*phi_z(i,j,k))*costeta(i,l,k)  &
                       )*ijacob3(i,j,k) + pf(i,l,k)*ir(i,l,k) )
             ut(i,l,k) = vg(i,l,k)*( &
                       ((duksi*ksi_x(i,j,k)+dueta(i,l,k)*eta_x(i,j,k)+duphi(i,l,k)*phi_x(i,j,k))*sintcosp(i,l,k) &
                       +(duksi*ksi_y(i,j,k)+dueta(i,l,k)*eta_y(i,j,k)+duphi(i,l,k)*phi_y(i,j,k))*sintsinp(i,l,k) &
                       +(duksi*ksi_z(i,j,k)+dueta(i,l,k)*eta_z(i,j,k)+duphi(i,l,k)*phi_z(i,j,k))*costeta(i,l,k)  &
                       )*ijacob3(i,j,k) + uf(i,l,k)*ir(i,l,k) )
             vt(i,l,k) = vg(i,l,k)*( &
                       ((dvksi*ksi_x(i,j,k)+dveta(i,l,k)*eta_x(i,j,k)+dvphi(i,l,k)*phi_x(i,j,k))*sintcosp(i,l,k) &
                       +(dvksi*ksi_y(i,j,k)+dveta(i,l,k)*eta_y(i,j,k)+dvphi(i,l,k)*phi_y(i,j,k))*sintsinp(i,l,k) &
                       +(dvksi*ksi_z(i,j,k)+dveta(i,l,k)*eta_z(i,j,k)+dvphi(i,l,k)*phi_z(i,j,k))*costeta(i,l,k)  &
                       )*ijacob3(i,j,k) + vf(i,l,k)*ir(i,l,k) )
             wt(i,l,k) = vg(i,l,k)*( &
                       ((dwksi*ksi_x(i,j,k)+dweta(i,l,k)*eta_x(i,j,k)+dwphi(i,l,k)*phi_x(i,j,k))*sintcosp(i,l,k) &
                       +(dwksi*ksi_y(i,j,k)+dweta(i,l,k)*eta_y(i,j,k)+dwphi(i,l,k)*phi_y(i,j,k))*sintsinp(i,l,k) &
                       +(dwksi*ksi_z(i,j,k)+dweta(i,l,k)*eta_z(i,j,k)+dwphi(i,l,k)*phi_z(i,j,k))*costeta(i,l,k)  &
                       )*ijacob3(i,j,k) + wf(i,l,k)*ir(i,l,k) )
             rt(i,l,k) = vg(i,l,k)*( &
                       ((drksi*ksi_x(i,j,k)+dreta(i,l,k)*eta_x(i,j,k)+drphi(i,l,k)*phi_x(i,j,k))*sintcosp(i,l,k) &
                       +(drksi*ksi_y(i,j,k)+dreta(i,l,k)*eta_y(i,j,k)+drphi(i,l,k)*phi_y(i,j,k))*sintsinp(i,l,k) &
                       +(drksi*ksi_z(i,j,k)+dreta(i,l,k)*eta_z(i,j,k)+drphi(i,l,k)*phi_z(i,j,k))*costeta(i,l,k)  &
                       )*ijacob3(i,j,k) + rf(i,l,k)*ir(i,l,k) )
          enddo
       enddo
    enddo

    ! Update fluxes at each RK step
    ! =============================
    do i=ndx,nfx
       do j=nymnghp1,ny
          l=j-nymngh
          do k=1,ngh
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

  end subroutine bc_TD3d_jmax_kmin_c3

  !===============================================================================
  module subroutine bc_TD3d_jmax_kmax_c3
  !===============================================================================
    !> 3D Tam & Dong's BC: boundary condition at jmax-kmax (edge 2,2,2 /top-back)
    !> - 3D curvilinear version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l,m
    real(wp) :: cp,av,c2_
    real(wp), dimension(nx1:nx2,-2:ngh,-2:ngh) :: rf,uf,vf,wf,pf
    real(wp), dimension(nx,1:ngh,1:ngh) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(nx,1:ngh,1:ngh) :: dueta,dveta,dweta,dpeta,dreta
    real(wp), dimension(nx,1:ngh,1:ngh) :: duphi,dvphi,dwphi,dpphi,drphi
    real(wp) :: duksi,dvksi,dwksi,dpksi,drksi
    !-------------------------------------------------------------------------

    ! Pointers for spherical coordinates
    ! ==================================
    ir=>BC_edge(2,2,2)%i_r
    cosphi=>BC_edge(2,2,2)%cosp
    sinphi=>BC_edge(2,2,2)%sinp
    costeta=>BC_edge(2,2,2)%cost
    sinteta=>BC_edge(2,2,2)%sint
    costcosp=>BC_edge(2,2,2)%costcosp
    costsinp=>BC_edge(2,2,2)%costsinp
    sintcosp=>BC_edge(2,2,2)%sintcosp
    sintsinp=>BC_edge(2,2,2)%sintsinp

    ! Compute fluctuations
    ! ====================
    do k=nzmngh-2,nz
       m=k-nzmngh
       do j=nymngh-2,ny
          l=j-nymngh
          do i=nx1,nx2
             rf(i,l,m)=rho_n(i,j,k)-BC_face(2,2)%U0(i,l,k,1)
             uf(i,l,m)=   uu(i,j,k)-BC_face(2,2)%U0(i,l,k,2)
             vf(i,l,m)=   vv(i,j,k)-BC_face(2,2)%U0(i,l,k,3)
             wf(i,l,m)=   ww(i,j,k)-BC_face(2,2)%U0(i,l,k,4)
             pf(i,l,m)=  prs(i,j,k)-BC_face(2,2)%U0(i,l,k,5)
          enddo
       enddo
    enddo

    ! Compute group velocity vg
    ! =========================
    ! vg= u0*sin(teta)*cos(phi)+v0*sin(teta)*sin(phi)+w0*cos(teta)
    !   + sqrt{c0^2-[u0*cos(teta)*cos(phi)+v0*cos(teta)*sin(phi)-w0*sin(teta)]^2-[u0*sin(phi)-v0*cos(phi)]^2}
    do k=nzmnghp1,nz
       m=k-nzmngh
       do l=1,ngh
          do i=ndx,nfx
             vg(i,l,m)= BC_face(2,2)%U0(i,l,k,2)*sintcosp(i,l,m) +     &
               BC_face(2,2)%U0(i,l,k,3)*sintsinp(i,l,m) +     &
               BC_face(2,2)%U0(i,l,k,4)*costeta(i,l,m)          &
                      + sqrt( BC_face(2,2)%U0(i,l,k,6)-                &
                       ( BC_face(2,2)%U0(i,l,k,2)*costcosp(i,l,m) +    &
                   BC_face(2,2)%U0(i,l,k,3)*costsinp(i,l,m) -    &
                     BC_face(2,2)%U0(i,l,k,4)*sinteta(i,l,m) )**2- &
                  ( BC_face(2,2)%U0(i,l,k,2)*sinphi(i,l,m) -      &
                     BC_face(2,2)%U0(i,l,k,3)*cosphi(i,l,m) )**2 )
          enddo
       enddo
    enddo

    ! Non-centered derivatives in y-direction *sin(teta)*sin(phi)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    do j=nymnghp1,ny-3
       l=j-nymngh
       do m=1,ngh
          do i=ndx,nfx
             dueta(i,l,m) = a7(1)*(uf(i,l+1,m)-uf(i,l-1,m)) &
                          + a7(2)*(uf(i,l+2,m)-uf(i,l-2,m)) &
                          + a7(3)*(uf(i,l+3,m)-uf(i,l-3,m))
             dveta(i,l,m) = a7(1)*(vf(i,l+1,m)-vf(i,l-1,m)) &
                          + a7(2)*(vf(i,l+2,m)-vf(i,l-2,m)) &
                          + a7(3)*(vf(i,l+3,m)-vf(i,l-3,m))
             dweta(i,l,m) = a7(1)*(wf(i,l+1,m)-wf(i,l-1,m)) &
                          + a7(2)*(wf(i,l+2,m)-wf(i,l-2,m)) &
                          + a7(3)*(wf(i,l+3,m)-wf(i,l-3,m))
             dpeta(i,l,m) = a7(1)*(pf(i,l+1,m)-pf(i,l-1,m)) &
                          + a7(2)*(pf(i,l+2,m)-pf(i,l-2,m)) &
                          + a7(3)*(pf(i,l+3,m)-pf(i,l-3,m))
             dreta(i,l,m) = a7(1)*(rf(i,l+1,m)-rf(i,l-1,m)) &
                          + a7(2)*(rf(i,l+2,m)-rf(i,l-2,m)) &
                          + a7(3)*(rf(i,l+3,m)-rf(i,l-3,m))
          enddo
       enddo
    enddo

    j=ny-2
    l=j-nymngh
    do m=1,ngh
       do i=ndx,nfx
          dueta(i,l,m) = a42(1)*uf(i,l+2,m)+a42(2)*uf(i,l+1,m) &
                       + a42(3)*uf(i,l  ,m)+a42(4)*uf(i,l-1,m) &
                       + a42(5)*uf(i,l-2,m)+a42(6)*uf(i,l-3,m) &
                       + a42(7)*uf(i,l-4,m)
          dveta(i,l,m) = a42(1)*vf(i,l+2,m)+a42(2)*vf(i,l+1,m) &
                       + a42(3)*vf(i,l  ,m)+a42(4)*vf(i,l-1,m) &
                       + a42(5)*vf(i,l-2,m)+a42(6)*vf(i,l-3,m) &
                       + a42(7)*vf(i,l-4,m)
          dweta(i,l,m) = a42(1)*wf(i,l+2,m)+a42(2)*wf(i,l+1,m) &
                       + a42(3)*wf(i,l  ,m)+a42(4)*wf(i,l-1,m) &
                       + a42(5)*wf(i,l-2,m)+a42(6)*wf(i,l-3,m) &
                       + a42(7)*wf(i,l-4,m)
          dpeta(i,l,m) = a42(1)*pf(i,l+2,m)+a42(2)*pf(i,l+1,m) &
                       + a42(3)*pf(i,l  ,m)+a42(4)*pf(i,l-1,m) &
                       + a42(5)*pf(i,l-2,m)+a42(6)*pf(i,l-3,m) &
                       + a42(7)*pf(i,l-4,m)
          dreta(i,l,m) = a42(1)*rf(i,l+2,m)+a42(2)*rf(i,l+1,m) &
                       + a42(3)*rf(i,l  ,m)+a42(4)*rf(i,l-1,m) &
                       + a42(5)*rf(i,l-2,m)+a42(6)*rf(i,l-3,m) &
                       + a42(7)*rf(i,l-4,m)
       enddo
    enddo

    j=ny-1
    l=j-nymngh
    do m=1,ngh
       do i=ndx,nfx
          dueta(i,l,m) = a51(1)*uf(i,l+1,m)+a51(2)*uf(i,l  ,m) &
                       + a51(3)*uf(i,l-1,m)+a51(4)*uf(i,l-2,m) &
                       + a51(5)*uf(i,l-3,m)+a51(6)*uf(i,l-4,m) &
                       + a51(7)*uf(i,l-5,m)
          dveta(i,l,m) = a51(1)*vf(i,l+1,m)+a51(2)*vf(i,l  ,m) &
                       + a51(3)*vf(i,l-1,m)+a51(4)*vf(i,l-2,m) &
                       + a51(5)*vf(i,l-3,m)+a51(6)*vf(i,l-4,m) &
                       + a51(7)*vf(i,l-5,m)
          dweta(i,l,m) = a51(1)*wf(i,l+1,m)+a51(2)*wf(i,l  ,m) &
                       + a51(3)*wf(i,l-1,m)+a51(4)*wf(i,l-2,m) &
                       + a51(5)*wf(i,l-3,m)+a51(6)*wf(i,l-4,m) &
                       + a51(7)*wf(i,l-5,m)
          dpeta(i,l,m) = a51(1)*pf(i,l+1,m)+a51(2)*pf(i,l  ,m) &
                       + a51(3)*pf(i,l-1,m)+a51(4)*pf(i,l-2,m) &
                       + a51(5)*pf(i,l-3,m)+a51(6)*pf(i,l-4,m) &
                       + a51(7)*pf(i,l-5,m)
          dreta(i,l,m) = a51(1)*rf(i,l+1,m)+a51(2)*rf(i,l  ,m) &
                       + a51(3)*rf(i,l-1,m)+a51(4)*rf(i,l-2,m) &
                       + a51(5)*rf(i,l-3,m)+a51(6)*rf(i,l-4,m) &
                       + a51(7)*rf(i,l-5,m)
       enddo
    enddo

    j=ny
    l=j-nymngh
    do m=1,ngh
       do i=ndx,nfx
          dueta(i,l,m) = a60(1)*uf(i,l  ,m)+a60(2)*uf(i,l-1,m) &
                       + a60(3)*uf(i,l-2,m)+a60(4)*uf(i,l-3,m) &
                       + a60(5)*uf(i,l-4,m)+a60(6)*uf(i,l-5,m) &
                       + a60(7)*uf(i,l-6,m)
          dveta(i,l,m) = a60(1)*vf(i,l  ,m)+a60(2)*vf(i,l-1,m) &
                       + a60(3)*vf(i,l-2,m)+a60(4)*vf(i,l-3,m) &
                       + a60(5)*vf(i,l-4,m)+a60(6)*vf(i,l-5,m) &
                       + a60(7)*vf(i,l-6,m)
          dweta(i,l,m) = a60(1)*wf(i,l  ,m)+a60(2)*wf(i,l-1,m) &
                       + a60(3)*wf(i,l-2,m)+a60(4)*wf(i,l-3,m) &
                       + a60(5)*wf(i,l-4,m)+a60(6)*wf(i,l-5,m) &
                       + a60(7)*wf(i,l-6,m)
          dpeta(i,l,m) = a60(1)*pf(i,l  ,m)+a60(2)*pf(i,l-1,m) &
                       + a60(3)*pf(i,l-2,m)+a60(4)*pf(i,l-3,m) &
                       + a60(5)*pf(i,l-4,m)+a60(6)*pf(i,l-5,m) &
                       + a60(7)*pf(i,l-6,m)
          dreta(i,l,m) = a60(1)*rf(i,l  ,m)+a60(2)*rf(i,l-1,m) &
                       + a60(3)*rf(i,l-2,m)+a60(4)*rf(i,l-3,m) &
                       + a60(5)*rf(i,l-4,m)+a60(6)*rf(i,l-5,m) &
                       + a60(7)*rf(i,l-6,m)
       enddo
    enddo

    ! Non-centered derivatives in z-direction *cos(teta)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    do k=nzmnghp1,nz-3
       m=k-nzmngh
       do l=1,ngh
          do i=ndx,nfx
             duphi(i,l,m)= a7(1)*(uf(i,l,m+1) - uf(i,l,m-1)) + &
                           a7(2)*(uf(i,l,m+2) - uf(i,l,m-2)) + &
                           a7(3)*(uf(i,l,m+3) - uf(i,l,m-3))
             dvphi(i,l,m)= a7(1)*(vf(i,l,m+1) - vf(i,l,m-1)) + &
                           a7(2)*(vf(i,l,m+2) - vf(i,l,m-2)) + &
                           a7(3)*(vf(i,l,m+3) - vf(i,l,m-3))
             dwphi(i,l,m)= a7(1)*(wf(i,l,m+1) - wf(i,l,m-1)) + &
                           a7(2)*(wf(i,l,m+2) - wf(i,l,m-2)) + &
                           a7(3)*(wf(i,l,m+3) - wf(i,l,m-3))
             dpphi(i,l,m)= a7(1)*(pf(i,l,m+1) - pf(i,l,m-1)) + &
                           a7(2)*(pf(i,l,m+2) - pf(i,l,m-2)) + &
                           a7(3)*(pf(i,l,m+3) - pf(i,l,m-3))
             drphi(i,l,m)= a7(1)*(rf(i,l,m+1) - rf(i,l,m-1)) + &
                           a7(2)*(rf(i,l,m+2) - rf(i,l,m-2)) + &
                           a7(3)*(rf(i,l,m+3) - rf(i,l,m-3))
          enddo
       enddo
    enddo

    k=nz-2
    m=k-nzmngh
    do l=1,ngh
       do i=ndx,nfx
          duphi(i,l,m) = a42(1)*uf(i,l,m+2)+a42(2)*uf(i,l,m+1) &
                       + a42(3)*uf(i,l,m  )+a42(4)*uf(i,l,m-1) &
                       + a42(5)*uf(i,l,m-2)+a42(6)*uf(i,l,m-3) &
                       + a42(7)*uf(i,l,m-4)
          dvphi(i,l,m) = a42(1)*vf(i,l,m+2)+a42(2)*vf(i,l,m+1) &
                       + a42(3)*vf(i,l,m  )+a42(4)*vf(i,l,m-1) &
                       + a42(5)*vf(i,l,m-2)+a42(6)*vf(i,l,m-3) &
                       + a42(7)*vf(i,l,m-4)
          dwphi(i,l,m) = a42(1)*wf(i,l,m+2)+a42(2)*wf(i,l,m+1) &
                       + a42(3)*wf(i,l,m  )+a42(4)*wf(i,l,m-1) &
                       + a42(5)*wf(i,l,m-2)+a42(6)*wf(i,l,m-3) &
                       + a42(7)*wf(i,l,m-4)
          dpphi(i,l,m) = a42(1)*pf(i,l,m+2)+a42(2)*pf(i,l,m+1) &
                       + a42(3)*pf(i,l,m  )+a42(4)*pf(i,l,m-1) &
                       + a42(5)*pf(i,l,m-2)+a42(6)*pf(i,l,m-3) &
                       + a42(7)*pf(i,l,m-4)
          drphi(i,l,m) = a42(1)*rf(i,l,m+2)+a42(2)*rf(i,l,m+1) &
                       + a42(3)*rf(i,l,m  )+a42(4)*rf(i,l,m-1) &
                       + a42(5)*rf(i,l,m-2)+a42(6)*rf(i,l,m-3) &
                       + a42(7)*rf(i,l,m-4)
       enddo
    enddo

    k=nz-1
    m=k-nzmngh
    do l=1,ngh
       do i=ndx,nfx
          duphi(i,l,m) = a51(1)*uf(i,l,m+1)+a51(2)*uf(i,l,m  ) &
                       + a51(3)*uf(i,l,m-1)+a51(4)*uf(i,l,m-2) &
                       + a51(5)*uf(i,l,m-3)+a51(6)*uf(i,l,m-4) &
                       + a51(7)*uf(i,l,m-5)
          dvphi(i,l,m) = a51(1)*vf(i,l,m+1)+a51(2)*vf(i,l,m  ) &
                       + a51(3)*vf(i,l,m-1)+a51(4)*vf(i,l,m-2) &
                       + a51(5)*vf(i,l,m-3)+a51(6)*vf(i,l,m-4) &
                       + a51(7)*vf(i,l,m-5)
          dwphi(i,l,m) = a51(1)*wf(i,l,m+1)+a51(2)*wf(i,l,m  ) &
                       + a51(3)*wf(i,l,m-1)+a51(4)*wf(i,l,m-2) &
                       + a51(5)*wf(i,l,m-3)+a51(6)*wf(i,l,m-4) &
                       + a51(7)*wf(i,l,m-5)
          dpphi(i,l,m) = a51(1)*pf(i,l,m+1)+a51(2)*pf(i,l,m  ) &
                       + a51(3)*pf(i,l,m-1)+a51(4)*pf(i,l,m-2) &
                       + a51(5)*pf(i,l,m-3)+a51(6)*pf(i,l,m-4) &
                       + a51(7)*pf(i,l,m-5)
          drphi(i,l,m) = a51(1)*rf(i,l,m+1)+a51(2)*rf(i,l,m  ) &
                       + a51(3)*rf(i,l,m-1)+a51(4)*rf(i,l,m-2) &
                       + a51(5)*rf(i,l,m-3)+a51(6)*rf(i,l,m-4) &
                       + a51(7)*rf(i,l,m-5)
       enddo
    enddo

    k=nz
    m=k-nzmngh
    do l=1,ngh
       do i=ndx,nfx
          duphi(i,l,m) = a60(1)*uf(i,l,m  )+a60(2)*uf(i,l,m-1) &
                       + a60(3)*uf(i,l,m-2)+a60(4)*uf(i,l,m-3) &
                       + a60(5)*uf(i,l,m-4)+a60(6)*uf(i,l,m-5) &
                       + a60(7)*uf(i,l,m-6)
          dvphi(i,l,m) = a60(1)*vf(i,l,m  )+a60(2)*vf(i,l,m-1) &
                       + a60(3)*vf(i,l,m-2)+a60(4)*vf(i,l,m-3) &
                       + a60(5)*vf(i,l,m-4)+a60(6)*vf(i,l,m-5) &
                       + a60(7)*vf(i,l,m-6)
          dwphi(i,l,m) = a60(1)*wf(i,l,m  )+a60(2)*wf(i,l,m-1) &
                       + a60(3)*wf(i,l,m-2)+a60(4)*wf(i,l,m-3) &
                       + a60(5)*wf(i,l,m-4)+a60(6)*wf(i,l,m-5) &
                       + a60(7)*wf(i,l,m-6)
          dpphi(i,l,m) = a60(1)*pf(i,l,m  )+a60(2)*pf(i,l,m-1) &
                       + a60(3)*pf(i,l,m-2)+a60(4)*pf(i,l,m-3) &
                       + a60(5)*pf(i,l,m-4)+a60(6)*pf(i,l,m-5) &
                       + a60(7)*pf(i,l,m-6)
          drphi(i,l,m) = a60(1)*rf(i,l,m  )+a60(2)*rf(i,l,m-1) &
                       + a60(3)*rf(i,l,m-2)+a60(4)*rf(i,l,m-3) &
                       + a60(5)*rf(i,l,m-4)+a60(6)*rf(i,l,m-5) &
                       + a60(7)*rf(i,l,m-6)
       enddo
    enddo

    ! Compute rt,ut,vt,wt & pt
    ! ========================
    do i=ndx,nfx
       do j=nymnghp1,ny
          l=j-nymngh
          do k=nzmnghp1,nz
             m=k-nzmngh
             duksi = a7(1)*( uf(i+1,l,m) - uf(i-1,l,m) ) + &
                     a7(2)*( uf(i+2,l,m) - uf(i-2,l,m) ) + &
                     a7(3)*( uf(i+3,l,m) - uf(i-3,l,m) )
             dvksi = a7(1)*( vf(i+1,l,m) - vf(i-1,l,m) ) + &
                     a7(2)*( vf(i+2,l,m) - vf(i-2,l,m) ) + &
                     a7(3)*( vf(i+3,l,m) - vf(i-3,l,m) )
             dwksi = a7(1)*( wf(i+1,l,m) - wf(i-1,l,m) ) + &
                     a7(2)*( wf(i+2,l,m) - wf(i-2,l,m) ) + &
                     a7(3)*( wf(i+3,l,m) - wf(i-3,l,m) )
             dpksi = a7(1)*( pf(i+1,l,m) - pf(i-1,l,m) ) + &
                     a7(2)*( pf(i+2,l,m) - pf(i-2,l,m) ) + &
                     a7(3)*( pf(i+3,l,m) - pf(i-3,l,m) )
             drksi = a7(1)*( rf(i+1,l,m) - rf(i-1,l,m) ) + &
                     a7(2)*( rf(i+2,l,m) - rf(i-2,l,m) ) + &
                     a7(3)*( rf(i+3,l,m) - rf(i-3,l,m) )

             pt(i,l,m) = vg(i,l,m)*( &
                       ((dpksi*ksi_x(i,j,k)+dpeta(i,l,m)*eta_x(i,j,k)+dpphi(i,l,m)*phi_x(i,j,k))*sintcosp(i,l,m) &
                       +(dpksi*ksi_y(i,j,k)+dpeta(i,l,m)*eta_y(i,j,k)+dpphi(i,l,m)*phi_y(i,j,k))*sintsinp(i,l,m) &
                       +(dpksi*ksi_z(i,j,k)+dpeta(i,l,m)*eta_z(i,j,k)+dpphi(i,l,m)*phi_z(i,j,k))*costeta(i,l,m)  &
                       )*ijacob3(i,j,k) + pf(i,l,m)*ir(i,l,m) )
             ut(i,l,m) = vg(i,l,m)*( &
                       ((duksi*ksi_x(i,j,k)+dueta(i,l,m)*eta_x(i,j,k)+duphi(i,l,m)*phi_x(i,j,k))*sintcosp(i,l,m) &
                       +(duksi*ksi_y(i,j,k)+dueta(i,l,m)*eta_y(i,j,k)+duphi(i,l,m)*phi_y(i,j,k))*sintsinp(i,l,m) &
                       +(duksi*ksi_z(i,j,k)+dueta(i,l,m)*eta_z(i,j,k)+duphi(i,l,m)*phi_z(i,j,k))*costeta(i,l,m)  &
                       )*ijacob3(i,j,k) + uf(i,l,m)*ir(i,l,m) )
             vt(i,l,m) = vg(i,l,m)*( &
                       ((dvksi*ksi_x(i,j,k)+dveta(i,l,m)*eta_x(i,j,k)+dvphi(i,l,m)*phi_x(i,j,k))*sintcosp(i,l,m) &
                       +(dvksi*ksi_y(i,j,k)+dveta(i,l,m)*eta_y(i,j,k)+dvphi(i,l,m)*phi_y(i,j,k))*sintsinp(i,l,m) &
                       +(dvksi*ksi_z(i,j,k)+dveta(i,l,m)*eta_z(i,j,k)+dvphi(i,l,m)*phi_z(i,j,k))*costeta(i,l,m)  &
                       )*ijacob3(i,j,k) + vf(i,l,m)*ir(i,l,m) )
             wt(i,l,m) = vg(i,l,m)*( &
                       ((dwksi*ksi_x(i,j,k)+dweta(i,l,m)*eta_x(i,j,k)+dwphi(i,l,m)*phi_x(i,j,k))*sintcosp(i,l,m) &
                       +(dwksi*ksi_y(i,j,k)+dweta(i,l,m)*eta_y(i,j,k)+dwphi(i,l,m)*phi_y(i,j,k))*sintsinp(i,l,m) &
                       +(dwksi*ksi_z(i,j,k)+dweta(i,l,m)*eta_z(i,j,k)+dwphi(i,l,m)*phi_z(i,j,k))*costeta(i,l,m)  &
                       )*ijacob3(i,j,k) + wf(i,l,m)*ir(i,l,m) )
             rt(i,l,m) = vg(i,l,m)*( &
                       ((drksi*ksi_x(i,j,k)+dreta(i,l,m)*eta_x(i,j,k)+drphi(i,l,m)*phi_x(i,j,k))*sintcosp(i,l,m) &
                       +(drksi*ksi_y(i,j,k)+dreta(i,l,m)*eta_y(i,j,k)+drphi(i,l,m)*phi_y(i,j,k))*sintsinp(i,l,m) &
                       +(drksi*ksi_z(i,j,k)+dreta(i,l,m)*eta_z(i,j,k)+drphi(i,l,m)*phi_z(i,j,k))*costeta(i,l,m)  &
                       )*ijacob3(i,j,k) + rf(i,l,m)*ir(i,l,m) )
          enddo
       enddo
    enddo

    ! Update fluxes at each RK step
    ! =============================
    do i=ndx,nfx
       do j=nymnghp1,ny
          l=j-nymngh
          do k=nzmnghp1,nz
             m=k-nzmngh
             cp=cpcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
             av=avcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
             c2_=c2calc_tro(Tmp(i,j,k),rho_n(i,j,k))

             Krho(i,j,k)  = rt(i,l,m)
             Krhou(i,j,k) = uu(i,j,k)*rt(i,l,m)+rho_n(i,j,k)*ut(i,l,m)
             Krhov(i,j,k) = vv(i,j,k)*rt(i,l,m)+rho_n(i,j,k)*vt(i,l,m)
             Krhow(i,j,k) = ww(i,j,k)*rt(i,l,m)+rho_n(i,j,k)*wt(i,l,m)
             Krhoe(i,j,k) = cp/av*(pt(i,l,m)/c2_-rt(i,l,m)) &
                  + (rhoe_n(i,j,k)+prs(i,j,k))/rho_n(i,j,k)*rt(i,l,m) &
                  + rho_n(i,j,k)*(uu(i,j,k)*ut(i,l,m)+vv(i,j,k)*vt(i,l,m)+ww(i,j,k)*wt(i,l,m))
          enddo
       enddo
    enddo

  end subroutine bc_TD3d_jmax_kmax_c3

end submodule smod_TamDong3d_edgex_c3
