!===============================================================================
submodule (mod_TamDong3d_c3) smod_TamDong3d_faces_c3
!===============================================================================
  !> author: XG
  !> date: April 2023
  !> Non-reflecting boundary conditions of Tam and Dong
  !> - 3D curvilinear version - routines for faces
!=============================================================================== 

contains

  !===============================================================================
  module subroutine bc_TD3d_imin_c3
  !===============================================================================
    !> 3D Tam & Dong's BC: boundary condition at imin (left) - 3D curvilinear version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: cp,av,inn1,ntm1
    real(wp) :: duksi,dvksi,dwksi,dpksi,drksi
    real(wp) :: dueta,dveta,dweta,dpeta,dreta
    real(wp) :: duphi,dvphi,dwphi,dpphi,drphi
    real(wp), dimension(1:ngh+3,ny1:ny2,nz1:nz2) :: rf,uf,vf,wf,pf
    real(wp), dimension(ngh,ny,nz) :: rt,ut,vt,wt,pt,vg,c2_
    !-------------------------------------------------------------------------

    ! Pointers for spherical coordinates
    ! ==================================
    ir=>BC_face(1,1)%i_r
    cosphi=>BC_face(1,1)%cosp
    sinphi=>BC_face(1,1)%sinp
    costeta=>BC_face(1,1)%cost
    sinteta=>BC_face(1,1)%sint
    costcosp=>BC_face(1,1)%costcosp
    costsinp=>BC_face(1,1)%costsinp
    sintcosp=>BC_face(1,1)%sintcosp
    sintsinp=>BC_face(1,1)%sintsinp

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
       do k=nz1,nz2
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
          do k=nz1,nz2
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
    ! vg= u0*sin(teta)*cos(phi)+v0*sin(teta)*sin(phi)+w0*cos(teta)
    !   + sqrt{c0^2-[u0*cos(teta)*cos(phi)+v0*cos(teta)*sin(phi)-w0*sin(teta)]^2-[u0*sin(phi)-v0*cos(phi)]^2}
    do i=1,ngh
       do j=ndy,nfy
          do k=ndz,nfz
             vg(i,j,k)= BC_face(1,1)%U0(i,j,k,2)*sintcosp(i,j,k) +     &
               BC_face(1,1)%U0(i,j,k,3)*sintsinp(i,j,k) +     &
               BC_face(1,1)%U0(i,j,k,4)*costeta(i,j,k)          &
                      + sqrt( BC_face(1,1)%U0(i,j,k,6)-                &
                       ( BC_face(1,1)%U0(i,j,k,2)*costcosp(i,j,k) +    &
                   BC_face(1,1)%U0(i,j,k,3)*costsinp(i,j,k) -    &
                     BC_face(1,1)%U0(i,j,k,4)*sinteta(i,j,k) )**2- &
                  ( BC_face(1,1)%U0(i,j,k,2)*sinphi(i,j,k) -      &
                     BC_face(1,1)%U0(i,j,k,3)*cosphi(i,j,k) )**2 )
          enddo
       enddo
    enddo

    ! Compute vg*[sin(teta)*cos(phi)*dq/dx+sin(teta)*sin(phi)*dq/dy+cos(teta)*dq/dz+q/r]
    ! ==================================================================================
    ! (Tam & Webb DRP schemes)
    i=1
    do j=ndy,nfy
       do k=ndz,nfz
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

          duphi = a7(1)*( uf(i,j,k+1) - uf(i,j,k-1) ) + &
                  a7(2)*( uf(i,j,k+2) - uf(i,j,k-2) ) + &
                  a7(3)*( uf(i,j,k+3) - uf(i,j,k-3) )
          dvphi = a7(1)*( vf(i,j,k+1) - vf(i,j,k-1) ) + &
                  a7(2)*( vf(i,j,k+2) - vf(i,j,k-2) ) + &
                  a7(3)*( vf(i,j,k+3) - vf(i,j,k-3) )
          dwphi = a7(1)*( wf(i,j,k+1) - wf(i,j,k-1) ) + &
                  a7(2)*( wf(i,j,k+2) - wf(i,j,k-2) ) + &
                  a7(3)*( wf(i,j,k+3) - wf(i,j,k-3) )
          dpphi = a7(1)*( pf(i,j,k+1) - pf(i,j,k-1) ) + &
                  a7(2)*( pf(i,j,k+2) - pf(i,j,k-2) ) + &
                  a7(3)*( pf(i,j,k+3) - pf(i,j,k-3) )
          drphi = a7(1)*( rf(i,j,k+1) - rf(i,j,k-1) ) + &
                  a7(2)*( rf(i,j,k+2) - rf(i,j,k-2) ) + &
                  a7(3)*( rf(i,j,k+3) - rf(i,j,k-3) )

          pt(i,j,k) = vg(i,j,k)*( &
                    ((dpksi*ksi_x(i,j,k)+dpeta*eta_x(i,j,k)+dpphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(dpksi*ksi_y(i,j,k)+dpeta*eta_y(i,j,k)+dpphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(dpksi*ksi_z(i,j,k)+dpeta*eta_z(i,j,k)+dpphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + pf(i,j,k)*ir(i,j,k) )
          ut(i,j,k) = vg(i,j,k)*( &
                    ((duksi*ksi_x(i,j,k)+dueta*eta_x(i,j,k)+duphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(duksi*ksi_y(i,j,k)+dueta*eta_y(i,j,k)+duphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(duksi*ksi_z(i,j,k)+dueta*eta_z(i,j,k)+duphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + uf(i,j,k)*ir(i,j,k) )
          vt(i,j,k) = vg(i,j,k)*( &
                    ((dvksi*ksi_x(i,j,k)+dveta*eta_x(i,j,k)+dvphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(dvksi*ksi_y(i,j,k)+dveta*eta_y(i,j,k)+dvphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(dvksi*ksi_z(i,j,k)+dveta*eta_z(i,j,k)+dvphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + vf(i,j,k)*ir(i,j,k) )
          wt(i,j,k) = vg(i,j,k)*( &
                    ((dwksi*ksi_x(i,j,k)+dweta*eta_x(i,j,k)+dwphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(dwksi*ksi_y(i,j,k)+dweta*eta_y(i,j,k)+dwphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(dwksi*ksi_z(i,j,k)+dweta*eta_z(i,j,k)+dwphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + wf(i,j,k)*ir(i,j,k) )
          rt(i,j,k) = vg(i,j,k)*( &
                    ((drksi*ksi_x(i,j,k)+dreta*eta_x(i,j,k)+drphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(drksi*ksi_y(i,j,k)+dreta*eta_y(i,j,k)+drphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(drksi*ksi_z(i,j,k)+dreta*eta_z(i,j,k)+drphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + rf(i,j,k)*ir(i,j,k) )
       enddo
    enddo

    i=2
    do j=ndy,nfy
       do k=ndz,nfz
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

          duphi = a7(1)*( uf(i,j,k+1) - uf(i,j,k-1) ) + &
                  a7(2)*( uf(i,j,k+2) - uf(i,j,k-2) ) + &
                  a7(3)*( uf(i,j,k+3) - uf(i,j,k-3) )
          dvphi = a7(1)*( vf(i,j,k+1) - vf(i,j,k-1) ) + &
                  a7(2)*( vf(i,j,k+2) - vf(i,j,k-2) ) + &
                  a7(3)*( vf(i,j,k+3) - vf(i,j,k-3) )
          dwphi = a7(1)*( wf(i,j,k+1) - wf(i,j,k-1) ) + &
                  a7(2)*( wf(i,j,k+2) - wf(i,j,k-2) ) + &
                  a7(3)*( wf(i,j,k+3) - wf(i,j,k-3) )
          dpphi = a7(1)*( pf(i,j,k+1) - pf(i,j,k-1) ) + &
                  a7(2)*( pf(i,j,k+2) - pf(i,j,k-2) ) + &
                  a7(3)*( pf(i,j,k+3) - pf(i,j,k-3) )
          drphi = a7(1)*( rf(i,j,k+1) - rf(i,j,k-1) ) + &
                  a7(2)*( rf(i,j,k+2) - rf(i,j,k-2) ) + &
                  a7(3)*( rf(i,j,k+3) - rf(i,j,k-3) )

          pt(i,j,k) = vg(i,j,k)*( &
                    ((dpksi*ksi_x(i,j,k)+dpeta*eta_x(i,j,k)+dpphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(dpksi*ksi_y(i,j,k)+dpeta*eta_y(i,j,k)+dpphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(dpksi*ksi_z(i,j,k)+dpeta*eta_z(i,j,k)+dpphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + pf(i,j,k)*ir(i,j,k) )
          ut(i,j,k) = vg(i,j,k)*( &
                    ((duksi*ksi_x(i,j,k)+dueta*eta_x(i,j,k)+duphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(duksi*ksi_y(i,j,k)+dueta*eta_y(i,j,k)+duphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(duksi*ksi_z(i,j,k)+dueta*eta_z(i,j,k)+duphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + uf(i,j,k)*ir(i,j,k) )
          vt(i,j,k) = vg(i,j,k)*( &
                    ((dvksi*ksi_x(i,j,k)+dveta*eta_x(i,j,k)+dvphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(dvksi*ksi_y(i,j,k)+dveta*eta_y(i,j,k)+dvphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(dvksi*ksi_z(i,j,k)+dveta*eta_z(i,j,k)+dvphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + vf(i,j,k)*ir(i,j,k) )
          wt(i,j,k) = vg(i,j,k)*( &
                    ((dwksi*ksi_x(i,j,k)+dweta*eta_x(i,j,k)+dwphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(dwksi*ksi_y(i,j,k)+dweta*eta_y(i,j,k)+dwphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(dwksi*ksi_z(i,j,k)+dweta*eta_z(i,j,k)+dwphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + wf(i,j,k)*ir(i,j,k) )
          rt(i,j,k) = vg(i,j,k)*( &
                    ((drksi*ksi_x(i,j,k)+dreta*eta_x(i,j,k)+drphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(drksi*ksi_y(i,j,k)+dreta*eta_y(i,j,k)+drphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(drksi*ksi_z(i,j,k)+dreta*eta_z(i,j,k)+drphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + rf(i,j,k)*ir(i,j,k) )
       enddo
    enddo

    i=3
    do j=ndy,nfy
       do k=ndz,nfz
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

          duphi = a7(1)*( uf(i,j,k+1) - uf(i,j,k-1) ) + &
                  a7(2)*( uf(i,j,k+2) - uf(i,j,k-2) ) + &
                  a7(3)*( uf(i,j,k+3) - uf(i,j,k-3) )
          dvphi = a7(1)*( vf(i,j,k+1) - vf(i,j,k-1) ) + &
                  a7(2)*( vf(i,j,k+2) - vf(i,j,k-2) ) + &
                  a7(3)*( vf(i,j,k+3) - vf(i,j,k-3) )
          dwphi = a7(1)*( wf(i,j,k+1) - wf(i,j,k-1) ) + &
                  a7(2)*( wf(i,j,k+2) - wf(i,j,k-2) ) + &
                  a7(3)*( wf(i,j,k+3) - wf(i,j,k-3) )
          dpphi = a7(1)*( pf(i,j,k+1) - pf(i,j,k-1) ) + &
                  a7(2)*( pf(i,j,k+2) - pf(i,j,k-2) ) + &
                  a7(3)*( pf(i,j,k+3) - pf(i,j,k-3) )
          drphi = a7(1)*( rf(i,j,k+1) - rf(i,j,k-1) ) + &
                  a7(2)*( rf(i,j,k+2) - rf(i,j,k-2) ) + &
                  a7(3)*( rf(i,j,k+3) - rf(i,j,k-3) )

          pt(i,j,k) = vg(i,j,k)*( &
                    ((dpksi*ksi_x(i,j,k)+dpeta*eta_x(i,j,k)+dpphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(dpksi*ksi_y(i,j,k)+dpeta*eta_y(i,j,k)+dpphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(dpksi*ksi_z(i,j,k)+dpeta*eta_z(i,j,k)+dpphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + pf(i,j,k)*ir(i,j,k) )
          ut(i,j,k) = vg(i,j,k)*( &
                    ((duksi*ksi_x(i,j,k)+dueta*eta_x(i,j,k)+duphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(duksi*ksi_y(i,j,k)+dueta*eta_y(i,j,k)+duphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(duksi*ksi_z(i,j,k)+dueta*eta_z(i,j,k)+duphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + uf(i,j,k)*ir(i,j,k) )
          vt(i,j,k) = vg(i,j,k)*( &
                    ((dvksi*ksi_x(i,j,k)+dveta*eta_x(i,j,k)+dvphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(dvksi*ksi_y(i,j,k)+dveta*eta_y(i,j,k)+dvphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(dvksi*ksi_z(i,j,k)+dveta*eta_z(i,j,k)+dvphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + vf(i,j,k)*ir(i,j,k) )
          wt(i,j,k) = vg(i,j,k)*( &
                    ((dwksi*ksi_x(i,j,k)+dweta*eta_x(i,j,k)+dwphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(dwksi*ksi_y(i,j,k)+dweta*eta_y(i,j,k)+dwphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(dwksi*ksi_z(i,j,k)+dweta*eta_z(i,j,k)+dwphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + wf(i,j,k)*ir(i,j,k) )
          rt(i,j,k) = vg(i,j,k)*( &
                    ((drksi*ksi_x(i,j,k)+dreta*eta_x(i,j,k)+drphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(drksi*ksi_y(i,j,k)+dreta*eta_y(i,j,k)+drphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(drksi*ksi_z(i,j,k)+dreta*eta_z(i,j,k)+drphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + rf(i,j,k)*ir(i,j,k) )
       enddo
    enddo

    do i=4,ngh
       do j=ndy,nfy
          do k=ndz,nfz
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

             duphi = a7(1)*( uf(i,j,k+1) - uf(i,j,k-1) ) + &
                     a7(2)*( uf(i,j,k+2) - uf(i,j,k-2) ) + &
                     a7(3)*( uf(i,j,k+3) - uf(i,j,k-3) )
             dvphi = a7(1)*( vf(i,j,k+1) - vf(i,j,k-1) ) + &
                     a7(2)*( vf(i,j,k+2) - vf(i,j,k-2) ) + &
                     a7(3)*( vf(i,j,k+3) - vf(i,j,k-3) )
             dwphi = a7(1)*( wf(i,j,k+1) - wf(i,j,k-1) ) + &
                     a7(2)*( wf(i,j,k+2) - wf(i,j,k-2) ) + &
                     a7(3)*( wf(i,j,k+3) - wf(i,j,k-3) )
             dpphi = a7(1)*( pf(i,j,k+1) - pf(i,j,k-1) ) + &
                     a7(2)*( pf(i,j,k+2) - pf(i,j,k-2) ) + &
                     a7(3)*( pf(i,j,k+3) - pf(i,j,k-3) )
             drphi = a7(1)*( rf(i,j,k+1) - rf(i,j,k-1) ) + &
                     a7(2)*( rf(i,j,k+2) - rf(i,j,k-2) ) + &
                     a7(3)*( rf(i,j,k+3) - rf(i,j,k-3) )

             pt(i,j,k) = vg(i,j,k)*( &
                       ((dpksi*ksi_x(i,j,k)+dpeta*eta_x(i,j,k)+dpphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                       +(dpksi*ksi_y(i,j,k)+dpeta*eta_y(i,j,k)+dpphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                       +(dpksi*ksi_z(i,j,k)+dpeta*eta_z(i,j,k)+dpphi*phi_z(i,j,k))*costeta(i,j,k)  &
                       )*ijacob3(i,j,k) + pf(i,j,k)*ir(i,j,k) )
             ut(i,j,k) = vg(i,j,k)*( &
                       ((duksi*ksi_x(i,j,k)+dueta*eta_x(i,j,k)+duphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                       +(duksi*ksi_y(i,j,k)+dueta*eta_y(i,j,k)+duphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                       +(duksi*ksi_z(i,j,k)+dueta*eta_z(i,j,k)+duphi*phi_z(i,j,k))*costeta(i,j,k)  &
                       )*ijacob3(i,j,k) + uf(i,j,k)*ir(i,j,k) )
             vt(i,j,k) = vg(i,j,k)*( &
                       ((dvksi*ksi_x(i,j,k)+dveta*eta_x(i,j,k)+dvphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                       +(dvksi*ksi_y(i,j,k)+dveta*eta_y(i,j,k)+dvphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                       +(dvksi*ksi_z(i,j,k)+dveta*eta_z(i,j,k)+dvphi*phi_z(i,j,k))*costeta(i,j,k)  &
                       )*ijacob3(i,j,k) + vf(i,j,k)*ir(i,j,k) )
             wt(i,j,k) = vg(i,j,k)*( &
                       ((dwksi*ksi_x(i,j,k)+dweta*eta_x(i,j,k)+dwphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                       +(dwksi*ksi_y(i,j,k)+dweta*eta_y(i,j,k)+dwphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                       +(dwksi*ksi_z(i,j,k)+dweta*eta_z(i,j,k)+dwphi*phi_z(i,j,k))*costeta(i,j,k)  &
                       )*ijacob3(i,j,k) + wf(i,j,k)*ir(i,j,k) )
             rt(i,j,k) = vg(i,j,k)*( &
                       ((drksi*ksi_x(i,j,k)+dreta*eta_x(i,j,k)+drphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                       +(drksi*ksi_y(i,j,k)+dreta*eta_y(i,j,k)+drphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                       +(drksi*ksi_z(i,j,k)+dreta*eta_z(i,j,k)+drphi*phi_z(i,j,k))*costeta(i,j,k)  &
                       )*ijacob3(i,j,k) + rf(i,j,k)*ir(i,j,k) )
          enddo
       enddo
    enddo

    ! Update fluxes at each RK step
    ! =============================
    do k=ndz,nfz
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

  end subroutine bc_TD3d_imin_c3

  !===============================================================================
  module subroutine bc_TD3d_imax_c3
  !===============================================================================
    !> 3D Tam & Dong's BC: boundary condition at imax (right) - 3D curvilinear version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l
    real(wp) :: cp,av,inn1,ntm1
    real(wp) :: duksi,dvksi,dwksi,dpksi,drksi
    real(wp) :: dueta,dveta,dweta,dpeta,dreta
    real(wp) :: duphi,dvphi,dwphi,dpphi,drphi
    real(wp), dimension(-2:ngh,ny1:ny2,nz1:nz2) :: rf,uf,vf,wf,pf
    real(wp), dimension(ngh,ny,nz) :: rt,ut,vt,wt,pt,vg,c2_
    !-------------------------------------------------------------------------

    ! Pointers for spherical coordinates
    ! ==================================
    ir=>BC_face(1,2)%i_r
    cosphi=>BC_face(1,2)%cosp
    sinphi=>BC_face(1,2)%sinp
    costeta=>BC_face(1,2)%cost
    sinteta=>BC_face(1,2)%sint
    costcosp=>BC_face(1,2)%costcosp
    costsinp=>BC_face(1,2)%costsinp
    sintcosp=>BC_face(1,2)%sintcosp
    sintsinp=>BC_face(1,2)%sintsinp

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
       do k=nz1,nz2
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
          do k=nz1,nz2
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
    ! vg= u0*sin(teta)*cos(phi)+v0*sin(teta)*sin(phi)+w0*cos(teta)
    !   + sqrt{c0^2-[u0*cos(teta)*cos(phi)+v0*cos(teta)*sin(phi)-w0*sin(teta)]^2-[u0*sin(phi)-v0*cos(phi)]^2}
    do l=1,ngh
       do j=ndy,nfy
          do k=ndz,nfz
             vg(l,j,k)= BC_face(1,2)%U0(l,j,k,2)*sintcosp(l,j,k) +     &
               BC_face(1,2)%U0(l,j,k,3)*sintsinp(l,j,k) +     &
               BC_face(1,2)%U0(l,j,k,4)*costeta(l,j,k)          &
                      + sqrt( BC_face(1,2)%U0(l,j,k,6)-                &
                       ( BC_face(1,2)%U0(l,j,k,2)*costcosp(l,j,k) +    &
                   BC_face(1,2)%U0(l,j,k,3)*costsinp(l,j,k) -    &
                     BC_face(1,2)%U0(l,j,k,4)*sinteta(l,j,k) )**2- &
                  ( BC_face(1,2)%U0(l,j,k,2)*sinphi(l,j,k) -      &
                     BC_face(1,2)%U0(l,j,k,3)*cosphi(l,j,k) )**2 )
          enddo
       enddo
    enddo

    ! Compute vg*[sin(teta)*cos(phi)*dq/dx+sin(teta)*sin(phi)*dq/dy+cos(teta)*dq/dz+q/r]
    ! ==================================================================================
    ! (Tam & Webb DRP schemes)
    do i=nxmnghp1,nx-3
       l=i-nxmngh
       do j=ndy,nfy
          do k=ndz,nfz
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

             duphi = a7(1)*( uf(l,j,k+1) - uf(l,j,k-1) ) + &
                     a7(2)*( uf(l,j,k+2) - uf(l,j,k-2) ) + &
                     a7(3)*( uf(l,j,k+3) - uf(l,j,k-3) )
             dvphi = a7(1)*( vf(l,j,k+1) - vf(l,j,k-1) ) + &
                     a7(2)*( vf(l,j,k+2) - vf(l,j,k-2) ) + &
                     a7(3)*( vf(l,j,k+3) - vf(l,j,k-3) )
             dwphi = a7(1)*( wf(l,j,k+1) - wf(l,j,k-1) ) + &
                     a7(2)*( wf(l,j,k+2) - wf(l,j,k-2) ) + &
                     a7(3)*( wf(l,j,k+3) - wf(l,j,k-3) )
             dpphi = a7(1)*( pf(l,j,k+1) - pf(l,j,k-1) ) + &
                     a7(2)*( pf(l,j,k+2) - pf(l,j,k-2) ) + &
                     a7(3)*( pf(l,j,k+3) - pf(l,j,k-3) )
             drphi = a7(1)*( rf(l,j,k+1) - rf(l,j,k-1) ) + &
                     a7(2)*( rf(l,j,k+2) - rf(l,j,k-2) ) + &
                     a7(3)*( rf(l,j,k+3) - rf(l,j,k-3) )

             pt(l,j,k) = vg(l,j,k)*( &
                       ((dpksi*ksi_x(i,j,k)+dpeta*eta_x(i,j,k)+dpphi*phi_x(i,j,k))*sintcosp(l,j,k) &
                       +(dpksi*ksi_y(i,j,k)+dpeta*eta_y(i,j,k)+dpphi*phi_y(i,j,k))*sintsinp(l,j,k) &
                       +(dpksi*ksi_z(i,j,k)+dpeta*eta_z(i,j,k)+dpphi*phi_z(i,j,k))*costeta(l,j,k)  &
                       )*ijacob3(i,j,k) + pf(l,j,k)*ir(l,j,k) )
             ut(l,j,k) = vg(l,j,k)*( &
                       ((duksi*ksi_x(i,j,k)+dueta*eta_x(i,j,k)+duphi*phi_x(i,j,k))*sintcosp(l,j,k) &
                       +(duksi*ksi_y(i,j,k)+dueta*eta_y(i,j,k)+duphi*phi_y(i,j,k))*sintsinp(l,j,k) &
                       +(duksi*ksi_z(i,j,k)+dueta*eta_z(i,j,k)+duphi*phi_z(i,j,k))*costeta(l,j,k)  &
                       )*ijacob3(i,j,k) + uf(l,j,k)*ir(l,j,k) )
             vt(l,j,k) = vg(l,j,k)*( &
                       ((dvksi*ksi_x(i,j,k)+dveta*eta_x(i,j,k)+dvphi*phi_x(i,j,k))*sintcosp(l,j,k) &
                       +(dvksi*ksi_y(i,j,k)+dveta*eta_y(i,j,k)+dvphi*phi_y(i,j,k))*sintsinp(l,j,k) &
                       +(dvksi*ksi_z(i,j,k)+dveta*eta_z(i,j,k)+dvphi*phi_z(i,j,k))*costeta(l,j,k)  &
                       )*ijacob3(i,j,k) + vf(l,j,k)*ir(l,j,k) )
             wt(l,j,k) = vg(l,j,k)*( &
                       ((dwksi*ksi_x(i,j,k)+dweta*eta_x(i,j,k)+dwphi*phi_x(i,j,k))*sintcosp(l,j,k) &
                       +(dwksi*ksi_y(i,j,k)+dweta*eta_y(i,j,k)+dwphi*phi_y(i,j,k))*sintsinp(l,j,k) &
                       +(dwksi*ksi_z(i,j,k)+dweta*eta_z(i,j,k)+dwphi*phi_z(i,j,k))*costeta(l,j,k)  &
                       )*ijacob3(i,j,k) + wf(l,j,k)*ir(l,j,k) )
             rt(l,j,k) = vg(l,j,k)*( &
                       ((drksi*ksi_x(i,j,k)+dreta*eta_x(i,j,k)+drphi*phi_x(i,j,k))*sintcosp(l,j,k) &
                       +(drksi*ksi_y(i,j,k)+dreta*eta_y(i,j,k)+drphi*phi_y(i,j,k))*sintsinp(l,j,k) &
                       +(drksi*ksi_z(i,j,k)+dreta*eta_z(i,j,k)+drphi*phi_z(i,j,k))*costeta(l,j,k)  &
                       )*ijacob3(i,j,k) + rf(l,j,k)*ir(l,j,k) )
          enddo
       enddo
    enddo

    i=nx-2
    l=i-nxmngh
    do j=ndy,nfy
       do k=ndz,nfz
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

          duphi = a7(1)*( uf(l,j,k+1) - uf(l,j,k-1) ) + &
                  a7(2)*( uf(l,j,k+2) - uf(l,j,k-2) ) + &
                  a7(3)*( uf(l,j,k+3) - uf(l,j,k-3) )
          dvphi = a7(1)*( vf(l,j,k+1) - vf(l,j,k-1) ) + &
                  a7(2)*( vf(l,j,k+2) - vf(l,j,k-2) ) + &
                  a7(3)*( vf(l,j,k+3) - vf(l,j,k-3) )
          dwphi = a7(1)*( wf(l,j,k+1) - wf(l,j,k-1) ) + &
                  a7(2)*( wf(l,j,k+2) - wf(l,j,k-2) ) + &
                  a7(3)*( wf(l,j,k+3) - wf(l,j,k-3) )
          dpphi = a7(1)*( pf(l,j,k+1) - pf(l,j,k-1) ) + &
                  a7(2)*( pf(l,j,k+2) - pf(l,j,k-2) ) + &
                  a7(3)*( pf(l,j,k+3) - pf(l,j,k-3) )
          drphi = a7(1)*( rf(l,j,k+1) - rf(l,j,k-1) ) + &
                  a7(2)*( rf(l,j,k+2) - rf(l,j,k-2) ) + &
                  a7(3)*( rf(l,j,k+3) - rf(l,j,k-3) )

          pt(l,j,k) = vg(l,j,k)*( &
                    ((dpksi*ksi_x(i,j,k)+dpeta*eta_x(i,j,k)+dpphi*phi_x(i,j,k))*sintcosp(l,j,k) &
                    +(dpksi*ksi_y(i,j,k)+dpeta*eta_y(i,j,k)+dpphi*phi_y(i,j,k))*sintsinp(l,j,k) &
                    +(dpksi*ksi_z(i,j,k)+dpeta*eta_z(i,j,k)+dpphi*phi_z(i,j,k))*costeta(l,j,k)  &
                    )*ijacob3(i,j,k) + pf(l,j,k)*ir(l,j,k) )
          ut(l,j,k) = vg(l,j,k)*( &
                    ((duksi*ksi_x(i,j,k)+dueta*eta_x(i,j,k)+duphi*phi_x(i,j,k))*sintcosp(l,j,k) &
                    +(duksi*ksi_y(i,j,k)+dueta*eta_y(i,j,k)+duphi*phi_y(i,j,k))*sintsinp(l,j,k) &
                    +(duksi*ksi_z(i,j,k)+dueta*eta_z(i,j,k)+duphi*phi_z(i,j,k))*costeta(l,j,k)  &
                    )*ijacob3(i,j,k) + uf(l,j,k)*ir(l,j,k) )
          vt(l,j,k) = vg(l,j,k)*( &
                    ((dvksi*ksi_x(i,j,k)+dveta*eta_x(i,j,k)+dvphi*phi_x(i,j,k))*sintcosp(l,j,k) &
                    +(dvksi*ksi_y(i,j,k)+dveta*eta_y(i,j,k)+dvphi*phi_y(i,j,k))*sintsinp(l,j,k) &
                    +(dvksi*ksi_z(i,j,k)+dveta*eta_z(i,j,k)+dvphi*phi_z(i,j,k))*costeta(l,j,k)  &
                    )*ijacob3(i,j,k) + vf(l,j,k)*ir(l,j,k) )
          wt(l,j,k) = vg(l,j,k)*( &
                    ((dwksi*ksi_x(i,j,k)+dweta*eta_x(i,j,k)+dwphi*phi_x(i,j,k))*sintcosp(l,j,k) &
                    +(dwksi*ksi_y(i,j,k)+dweta*eta_y(i,j,k)+dwphi*phi_y(i,j,k))*sintsinp(l,j,k) &
                    +(dwksi*ksi_z(i,j,k)+dweta*eta_z(i,j,k)+dwphi*phi_z(i,j,k))*costeta(l,j,k)  &
                    )*ijacob3(i,j,k) + wf(l,j,k)*ir(l,j,k) )
          rt(l,j,k) = vg(l,j,k)*( &
                    ((drksi*ksi_x(i,j,k)+dreta*eta_x(i,j,k)+drphi*phi_x(i,j,k))*sintcosp(l,j,k) &
                    +(drksi*ksi_y(i,j,k)+dreta*eta_y(i,j,k)+drphi*phi_y(i,j,k))*sintsinp(l,j,k) &
                    +(drksi*ksi_z(i,j,k)+dreta*eta_z(i,j,k)+drphi*phi_z(i,j,k))*costeta(l,j,k)  &
                    )*ijacob3(i,j,k) + rf(l,j,k)*ir(l,j,k) )
       enddo
    enddo

    i=nx-1
    l=i-nxmngh
    do j=ndy,nfy
       do k=ndz,nfz
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

          duphi = a7(1)*( uf(l,j,k+1) - uf(l,j,k-1) ) + &
                  a7(2)*( uf(l,j,k+2) - uf(l,j,k-2) ) + &
                  a7(3)*( uf(l,j,k+3) - uf(l,j,k-3) )
          dvphi = a7(1)*( vf(l,j,k+1) - vf(l,j,k-1) ) + &
                  a7(2)*( vf(l,j,k+2) - vf(l,j,k-2) ) + &
                  a7(3)*( vf(l,j,k+3) - vf(l,j,k-3) )
          dwphi = a7(1)*( wf(l,j,k+1) - wf(l,j,k-1) ) + &
                  a7(2)*( wf(l,j,k+2) - wf(l,j,k-2) ) + &
                  a7(3)*( wf(l,j,k+3) - wf(l,j,k-3) )
          dpphi = a7(1)*( pf(l,j,k+1) - pf(l,j,k-1) ) + &
                  a7(2)*( pf(l,j,k+2) - pf(l,j,k-2) ) + &
                  a7(3)*( pf(l,j,k+3) - pf(l,j,k-3) )
          drphi = a7(1)*( rf(l,j,k+1) - rf(l,j,k-1) ) + &
                  a7(2)*( rf(l,j,k+2) - rf(l,j,k-2) ) + &
                  a7(3)*( rf(l,j,k+3) - rf(l,j,k-3) )

          pt(l,j,k) = vg(l,j,k)*( &
                    ((dpksi*ksi_x(i,j,k)+dpeta*eta_x(i,j,k)+dpphi*phi_x(i,j,k))*sintcosp(l,j,k) &
                    +(dpksi*ksi_y(i,j,k)+dpeta*eta_y(i,j,k)+dpphi*phi_y(i,j,k))*sintsinp(l,j,k) &
                    +(dpksi*ksi_z(i,j,k)+dpeta*eta_z(i,j,k)+dpphi*phi_z(i,j,k))*costeta(l,j,k)  &
                    )*ijacob3(i,j,k) + pf(l,j,k)*ir(l,j,k) )
          ut(l,j,k) = vg(l,j,k)*( &
                    ((duksi*ksi_x(i,j,k)+dueta*eta_x(i,j,k)+duphi*phi_x(i,j,k))*sintcosp(l,j,k) &
                    +(duksi*ksi_y(i,j,k)+dueta*eta_y(i,j,k)+duphi*phi_y(i,j,k))*sintsinp(l,j,k) &
                    +(duksi*ksi_z(i,j,k)+dueta*eta_z(i,j,k)+duphi*phi_z(i,j,k))*costeta(l,j,k)  &
                    )*ijacob3(i,j,k) + uf(l,j,k)*ir(l,j,k) )
          vt(l,j,k) = vg(l,j,k)*( &
                    ((dvksi*ksi_x(i,j,k)+dveta*eta_x(i,j,k)+dvphi*phi_x(i,j,k))*sintcosp(l,j,k) &
                    +(dvksi*ksi_y(i,j,k)+dveta*eta_y(i,j,k)+dvphi*phi_y(i,j,k))*sintsinp(l,j,k) &
                    +(dvksi*ksi_z(i,j,k)+dveta*eta_z(i,j,k)+dvphi*phi_z(i,j,k))*costeta(l,j,k)  &
                    )*ijacob3(i,j,k) + vf(l,j,k)*ir(l,j,k) )
          wt(l,j,k) = vg(l,j,k)*( &
                    ((dwksi*ksi_x(i,j,k)+dweta*eta_x(i,j,k)+dwphi*phi_x(i,j,k))*sintcosp(l,j,k) &
                    +(dwksi*ksi_y(i,j,k)+dweta*eta_y(i,j,k)+dwphi*phi_y(i,j,k))*sintsinp(l,j,k) &
                    +(dwksi*ksi_z(i,j,k)+dweta*eta_z(i,j,k)+dwphi*phi_z(i,j,k))*costeta(l,j,k)  &
                    )*ijacob3(i,j,k) + wf(l,j,k)*ir(l,j,k) )
          rt(l,j,k) = vg(l,j,k)*( &
                    ((drksi*ksi_x(i,j,k)+dreta*eta_x(i,j,k)+drphi*phi_x(i,j,k))*sintcosp(l,j,k) &
                    +(drksi*ksi_y(i,j,k)+dreta*eta_y(i,j,k)+drphi*phi_y(i,j,k))*sintsinp(l,j,k) &
                    +(drksi*ksi_z(i,j,k)+dreta*eta_z(i,j,k)+drphi*phi_z(i,j,k))*costeta(l,j,k)  &
                    )*ijacob3(i,j,k) + rf(l,j,k)*ir(l,j,k) )
       enddo
    enddo

    i=nx
    l=i-nxmngh
    do j=ndy,nfy
       do k=ndz,nfz
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

          duphi = a7(1)*( uf(l,j,k+1) - uf(l,j,k-1) ) + &
                  a7(2)*( uf(l,j,k+2) - uf(l,j,k-2) ) + &
                  a7(3)*( uf(l,j,k+3) - uf(l,j,k-3) )
          dvphi = a7(1)*( vf(l,j,k+1) - vf(l,j,k-1) ) + &
                  a7(2)*( vf(l,j,k+2) - vf(l,j,k-2) ) + &
                  a7(3)*( vf(l,j,k+3) - vf(l,j,k-3) )
          dwphi = a7(1)*( wf(l,j,k+1) - wf(l,j,k-1) ) + &
                  a7(2)*( wf(l,j,k+2) - wf(l,j,k-2) ) + &
                  a7(3)*( wf(l,j,k+3) - wf(l,j,k-3) )
          dpphi = a7(1)*( pf(l,j,k+1) - pf(l,j,k-1) ) + &
                  a7(2)*( pf(l,j,k+2) - pf(l,j,k-2) ) + &
                  a7(3)*( pf(l,j,k+3) - pf(l,j,k-3) )
          drphi = a7(1)*( rf(l,j,k+1) - rf(l,j,k-1) ) + &
                  a7(2)*( rf(l,j,k+2) - rf(l,j,k-2) ) + &
                  a7(3)*( rf(l,j,k+3) - rf(l,j,k-3) )

          pt(l,j,k) = vg(l,j,k)*( &
                    ((dpksi*ksi_x(i,j,k)+dpeta*eta_x(i,j,k)+dpphi*phi_x(i,j,k))*sintcosp(l,j,k) &
                    +(dpksi*ksi_y(i,j,k)+dpeta*eta_y(i,j,k)+dpphi*phi_y(i,j,k))*sintsinp(l,j,k) &
                    +(dpksi*ksi_z(i,j,k)+dpeta*eta_z(i,j,k)+dpphi*phi_z(i,j,k))*costeta(l,j,k)  &
                    )*ijacob3(i,j,k) + pf(l,j,k)*ir(l,j,k) )
          ut(l,j,k) = vg(l,j,k)*( &
                    ((duksi*ksi_x(i,j,k)+dueta*eta_x(i,j,k)+duphi*phi_x(i,j,k))*sintcosp(l,j,k) &
                    +(duksi*ksi_y(i,j,k)+dueta*eta_y(i,j,k)+duphi*phi_y(i,j,k))*sintsinp(l,j,k) &
                    +(duksi*ksi_z(i,j,k)+dueta*eta_z(i,j,k)+duphi*phi_z(i,j,k))*costeta(l,j,k)  &
                    )*ijacob3(i,j,k) + uf(l,j,k)*ir(l,j,k) )
          vt(l,j,k) = vg(l,j,k)*( &
                    ((dvksi*ksi_x(i,j,k)+dveta*eta_x(i,j,k)+dvphi*phi_x(i,j,k))*sintcosp(l,j,k) &
                    +(dvksi*ksi_y(i,j,k)+dveta*eta_y(i,j,k)+dvphi*phi_y(i,j,k))*sintsinp(l,j,k) &
                    +(dvksi*ksi_z(i,j,k)+dveta*eta_z(i,j,k)+dvphi*phi_z(i,j,k))*costeta(l,j,k)  &
                    )*ijacob3(i,j,k) + vf(l,j,k)*ir(l,j,k) )
          wt(l,j,k) = vg(l,j,k)*( &
                    ((dwksi*ksi_x(i,j,k)+dweta*eta_x(i,j,k)+dwphi*phi_x(i,j,k))*sintcosp(l,j,k) &
                    +(dwksi*ksi_y(i,j,k)+dweta*eta_y(i,j,k)+dwphi*phi_y(i,j,k))*sintsinp(l,j,k) &
                    +(dwksi*ksi_z(i,j,k)+dweta*eta_z(i,j,k)+dwphi*phi_z(i,j,k))*costeta(l,j,k)  &
                    )*ijacob3(i,j,k) + wf(l,j,k)*ir(l,j,k) )
          rt(l,j,k) = vg(l,j,k)*( &
                    ((drksi*ksi_x(i,j,k)+dreta*eta_x(i,j,k)+drphi*phi_x(i,j,k))*sintcosp(l,j,k) &
                    +(drksi*ksi_y(i,j,k)+dreta*eta_y(i,j,k)+drphi*phi_y(i,j,k))*sintsinp(l,j,k) &
                    +(drksi*ksi_z(i,j,k)+dreta*eta_z(i,j,k)+drphi*phi_z(i,j,k))*costeta(l,j,k)  &
                    )*ijacob3(i,j,k) + rf(l,j,k)*ir(l,j,k) )
       enddo
    enddo

    ! Update fluxes at each RK step
    ! =============================
    do k=ndz,nfz
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

  end subroutine bc_TD3d_imax_c3

  !===============================================================================
  module subroutine bc_TD3d_jmin_c3
  !===============================================================================
    !> 3D Tam & Dong's BC: boundary condition at jmin (bottom) - 3D curvilinear version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: cp,av,inn1,ntm1
    real(wp) :: duksi,dvksi,dwksi,dpksi,drksi
    real(wp) :: dueta,dveta,dweta,dpeta,dreta
    real(wp) :: duphi,dvphi,dwphi,dpphi,drphi
    real(wp), dimension(nx1:nx2,1:ngh+3,nz1:nz2) :: rf,uf,vf,wf,pf
    real(wp), dimension(nx,ngh,nz) :: rt,ut,vt,wt,pt,vg,c2_
    !-------------------------------------------------------------------------

    ! Pointers for spherical coordinates
    ! ==================================
    ir=>BC_face(2,1)%i_r
    cosphi=>BC_face(2,1)%cosp
    sinphi=>BC_face(2,1)%sinp
    costeta=>BC_face(2,1)%cost
    sinteta=>BC_face(2,1)%sint
    costcosp=>BC_face(2,1)%costcosp
    costsinp=>BC_face(2,1)%costsinp
    sintcosp=>BC_face(2,1)%sintcosp
    sintsinp=>BC_face(2,1)%sintsinp

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
    ! vg= u0*sin(teta)*cos(phi)+v0*sin(teta)*sin(phi)+w0*cos(teta)
    !   + sqrt{c0^2-[u0*cos(teta)*cos(phi)+v0*cos(teta)*sin(phi)-w0*sin(teta)]^2-[u0*sin(phi)-v0*cos(phi)]^2}
    do j=1,ngh
       do i=ndx,nfx
          do k=ndz,nfz
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

    ! Compute vg*[sin(teta)*cos(phi)*dq/dx+sin(teta)*sin(phi)*dq/dy+cos(teta)*dq/dz+q/r]
    ! ==================================================================================
    ! (Tam & Webb DRP schemes)
    j=1
    do i=ndx,nfx
       do k=ndz,nfz
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

          duphi = a7(1)*( uf(i,j,k+1) - uf(i,j,k-1) ) + &
                  a7(2)*( uf(i,j,k+2) - uf(i,j,k-2) ) + &
                  a7(3)*( uf(i,j,k+3) - uf(i,j,k-3) )
          dvphi = a7(1)*( vf(i,j,k+1) - vf(i,j,k-1) ) + &
                  a7(2)*( vf(i,j,k+2) - vf(i,j,k-2) ) + &
                  a7(3)*( vf(i,j,k+3) - vf(i,j,k-3) )
          dwphi = a7(1)*( wf(i,j,k+1) - wf(i,j,k-1) ) + &
                  a7(2)*( wf(i,j,k+2) - wf(i,j,k-2) ) + &
                  a7(3)*( wf(i,j,k+3) - wf(i,j,k-3) )
          dpphi = a7(1)*( pf(i,j,k+1) - pf(i,j,k-1) ) + &
                  a7(2)*( pf(i,j,k+2) - pf(i,j,k-2) ) + &
                  a7(3)*( pf(i,j,k+3) - pf(i,j,k-3) )
          drphi = a7(1)*( rf(i,j,k+1) - rf(i,j,k-1) ) + &
                  a7(2)*( rf(i,j,k+2) - rf(i,j,k-2) ) + &
                  a7(3)*( rf(i,j,k+3) - rf(i,j,k-3) )

          pt(i,j,k) = vg(i,j,k)*( &
                    ((dpksi*ksi_x(i,j,k)+dpeta*eta_x(i,j,k)+dpphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(dpksi*ksi_y(i,j,k)+dpeta*eta_y(i,j,k)+dpphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(dpksi*ksi_z(i,j,k)+dpeta*eta_z(i,j,k)+dpphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + pf(i,j,k)*ir(i,j,k) )
          ut(i,j,k) = vg(i,j,k)*( &
                    ((duksi*ksi_x(i,j,k)+dueta*eta_x(i,j,k)+duphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(duksi*ksi_y(i,j,k)+dueta*eta_y(i,j,k)+duphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(duksi*ksi_z(i,j,k)+dueta*eta_z(i,j,k)+duphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + uf(i,j,k)*ir(i,j,k) )
          vt(i,j,k) = vg(i,j,k)*( &
                    ((dvksi*ksi_x(i,j,k)+dveta*eta_x(i,j,k)+dvphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(dvksi*ksi_y(i,j,k)+dveta*eta_y(i,j,k)+dvphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(dvksi*ksi_z(i,j,k)+dveta*eta_z(i,j,k)+dvphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + vf(i,j,k)*ir(i,j,k) )
          wt(i,j,k) = vg(i,j,k)*( &
                    ((dwksi*ksi_x(i,j,k)+dweta*eta_x(i,j,k)+dwphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(dwksi*ksi_y(i,j,k)+dweta*eta_y(i,j,k)+dwphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(dwksi*ksi_z(i,j,k)+dweta*eta_z(i,j,k)+dwphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + wf(i,j,k)*ir(i,j,k) )
          rt(i,j,k) = vg(i,j,k)*( &
                    ((drksi*ksi_x(i,j,k)+dreta*eta_x(i,j,k)+drphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(drksi*ksi_y(i,j,k)+dreta*eta_y(i,j,k)+drphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(drksi*ksi_z(i,j,k)+dreta*eta_z(i,j,k)+drphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + rf(i,j,k)*ir(i,j,k) )
       enddo
    enddo

    j=2
    do i=ndx,nfx
       do k=ndz,nfz
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

          duphi = a7(1)*( uf(i,j,k+1) - uf(i,j,k-1) ) + &
                  a7(2)*( uf(i,j,k+2) - uf(i,j,k-2) ) + &
                  a7(3)*( uf(i,j,k+3) - uf(i,j,k-3) )
          dvphi = a7(1)*( vf(i,j,k+1) - vf(i,j,k-1) ) + &
                  a7(2)*( vf(i,j,k+2) - vf(i,j,k-2) ) + &
                  a7(3)*( vf(i,j,k+3) - vf(i,j,k-3) )
          dwphi = a7(1)*( wf(i,j,k+1) - wf(i,j,k-1) ) + &
                  a7(2)*( wf(i,j,k+2) - wf(i,j,k-2) ) + &
                  a7(3)*( wf(i,j,k+3) - wf(i,j,k-3) )
          dpphi = a7(1)*( pf(i,j,k+1) - pf(i,j,k-1) ) + &
                  a7(2)*( pf(i,j,k+2) - pf(i,j,k-2) ) + &
                  a7(3)*( pf(i,j,k+3) - pf(i,j,k-3) )
          drphi = a7(1)*( rf(i,j,k+1) - rf(i,j,k-1) ) + &
                  a7(2)*( rf(i,j,k+2) - rf(i,j,k-2) ) + &
                  a7(3)*( rf(i,j,k+3) - rf(i,j,k-3) )

          pt(i,j,k) = vg(i,j,k)*( &
                    ((dpksi*ksi_x(i,j,k)+dpeta*eta_x(i,j,k)+dpphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(dpksi*ksi_y(i,j,k)+dpeta*eta_y(i,j,k)+dpphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(dpksi*ksi_z(i,j,k)+dpeta*eta_z(i,j,k)+dpphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + pf(i,j,k)*ir(i,j,k) )
          ut(i,j,k) = vg(i,j,k)*( &
                    ((duksi*ksi_x(i,j,k)+dueta*eta_x(i,j,k)+duphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(duksi*ksi_y(i,j,k)+dueta*eta_y(i,j,k)+duphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(duksi*ksi_z(i,j,k)+dueta*eta_z(i,j,k)+duphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + uf(i,j,k)*ir(i,j,k) )
          vt(i,j,k) = vg(i,j,k)*( &
                    ((dvksi*ksi_x(i,j,k)+dveta*eta_x(i,j,k)+dvphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(dvksi*ksi_y(i,j,k)+dveta*eta_y(i,j,k)+dvphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(dvksi*ksi_z(i,j,k)+dveta*eta_z(i,j,k)+dvphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + vf(i,j,k)*ir(i,j,k) )
          wt(i,j,k) = vg(i,j,k)*( &
                    ((dwksi*ksi_x(i,j,k)+dweta*eta_x(i,j,k)+dwphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(dwksi*ksi_y(i,j,k)+dweta*eta_y(i,j,k)+dwphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(dwksi*ksi_z(i,j,k)+dweta*eta_z(i,j,k)+dwphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + wf(i,j,k)*ir(i,j,k) )
          rt(i,j,k) = vg(i,j,k)*( &
                    ((drksi*ksi_x(i,j,k)+dreta*eta_x(i,j,k)+drphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(drksi*ksi_y(i,j,k)+dreta*eta_y(i,j,k)+drphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(drksi*ksi_z(i,j,k)+dreta*eta_z(i,j,k)+drphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + rf(i,j,k)*ir(i,j,k) )
       enddo
    enddo

    j=3
    do i=ndx,nfx
       do k=ndz,nfz
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

          duphi = a7(1)*( uf(i,j,k+1) - uf(i,j,k-1) ) + &
                  a7(2)*( uf(i,j,k+2) - uf(i,j,k-2) ) + &
                  a7(3)*( uf(i,j,k+3) - uf(i,j,k-3) )
          dvphi = a7(1)*( vf(i,j,k+1) - vf(i,j,k-1) ) + &
                  a7(2)*( vf(i,j,k+2) - vf(i,j,k-2) ) + &
                  a7(3)*( vf(i,j,k+3) - vf(i,j,k-3) )
          dwphi = a7(1)*( wf(i,j,k+1) - wf(i,j,k-1) ) + &
                  a7(2)*( wf(i,j,k+2) - wf(i,j,k-2) ) + &
                  a7(3)*( wf(i,j,k+3) - wf(i,j,k-3) )
          dpphi = a7(1)*( pf(i,j,k+1) - pf(i,j,k-1) ) + &
                  a7(2)*( pf(i,j,k+2) - pf(i,j,k-2) ) + &
                  a7(3)*( pf(i,j,k+3) - pf(i,j,k-3) )
          drphi = a7(1)*( rf(i,j,k+1) - rf(i,j,k-1) ) + &
                  a7(2)*( rf(i,j,k+2) - rf(i,j,k-2) ) + &
                  a7(3)*( rf(i,j,k+3) - rf(i,j,k-3) )

          pt(i,j,k) = vg(i,j,k)*( &
                    ((dpksi*ksi_x(i,j,k)+dpeta*eta_x(i,j,k)+dpphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(dpksi*ksi_y(i,j,k)+dpeta*eta_y(i,j,k)+dpphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(dpksi*ksi_z(i,j,k)+dpeta*eta_z(i,j,k)+dpphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + pf(i,j,k)*ir(i,j,k) )
          ut(i,j,k) = vg(i,j,k)*( &
                    ((duksi*ksi_x(i,j,k)+dueta*eta_x(i,j,k)+duphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(duksi*ksi_y(i,j,k)+dueta*eta_y(i,j,k)+duphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(duksi*ksi_z(i,j,k)+dueta*eta_z(i,j,k)+duphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + uf(i,j,k)*ir(i,j,k) )
          vt(i,j,k) = vg(i,j,k)*( &
                    ((dvksi*ksi_x(i,j,k)+dveta*eta_x(i,j,k)+dvphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(dvksi*ksi_y(i,j,k)+dveta*eta_y(i,j,k)+dvphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(dvksi*ksi_z(i,j,k)+dveta*eta_z(i,j,k)+dvphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + vf(i,j,k)*ir(i,j,k) )
          wt(i,j,k) = vg(i,j,k)*( &
                    ((dwksi*ksi_x(i,j,k)+dweta*eta_x(i,j,k)+dwphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(dwksi*ksi_y(i,j,k)+dweta*eta_y(i,j,k)+dwphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(dwksi*ksi_z(i,j,k)+dweta*eta_z(i,j,k)+dwphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + wf(i,j,k)*ir(i,j,k) )
          rt(i,j,k) = vg(i,j,k)*( &
                    ((drksi*ksi_x(i,j,k)+dreta*eta_x(i,j,k)+drphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(drksi*ksi_y(i,j,k)+dreta*eta_y(i,j,k)+drphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(drksi*ksi_z(i,j,k)+dreta*eta_z(i,j,k)+drphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + rf(i,j,k)*ir(i,j,k) )
       enddo
    enddo

    do j=4,ngh
       do i=ndx,nfx
          do k=ndz,nfz
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

             duphi = a7(1)*( uf(i,j,k+1) - uf(i,j,k-1) ) + &
                     a7(2)*( uf(i,j,k+2) - uf(i,j,k-2) ) + &
                     a7(3)*( uf(i,j,k+3) - uf(i,j,k-3) )
             dvphi = a7(1)*( vf(i,j,k+1) - vf(i,j,k-1) ) + &
                     a7(2)*( vf(i,j,k+2) - vf(i,j,k-2) ) + &
                     a7(3)*( vf(i,j,k+3) - vf(i,j,k-3) )
             dwphi = a7(1)*( wf(i,j,k+1) - wf(i,j,k-1) ) + &
                     a7(2)*( wf(i,j,k+2) - wf(i,j,k-2) ) + &
                     a7(3)*( wf(i,j,k+3) - wf(i,j,k-3) )
             dpphi = a7(1)*( pf(i,j,k+1) - pf(i,j,k-1) ) + &
                     a7(2)*( pf(i,j,k+2) - pf(i,j,k-2) ) + &
                     a7(3)*( pf(i,j,k+3) - pf(i,j,k-3) )
             drphi = a7(1)*( rf(i,j,k+1) - rf(i,j,k-1) ) + &
                     a7(2)*( rf(i,j,k+2) - rf(i,j,k-2) ) + &
                     a7(3)*( rf(i,j,k+3) - rf(i,j,k-3) )

             pt(i,j,k) = vg(i,j,k)*( &
                       ((dpksi*ksi_x(i,j,k)+dpeta*eta_x(i,j,k)+dpphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                       +(dpksi*ksi_y(i,j,k)+dpeta*eta_y(i,j,k)+dpphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                       +(dpksi*ksi_z(i,j,k)+dpeta*eta_z(i,j,k)+dpphi*phi_z(i,j,k))*costeta(i,j,k)  &
                       )*ijacob3(i,j,k) + pf(i,j,k)*ir(i,j,k) )
             ut(i,j,k) = vg(i,j,k)*( &
                       ((duksi*ksi_x(i,j,k)+dueta*eta_x(i,j,k)+duphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                       +(duksi*ksi_y(i,j,k)+dueta*eta_y(i,j,k)+duphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                       +(duksi*ksi_z(i,j,k)+dueta*eta_z(i,j,k)+duphi*phi_z(i,j,k))*costeta(i,j,k)  &
                       )*ijacob3(i,j,k) + uf(i,j,k)*ir(i,j,k) )
             vt(i,j,k) = vg(i,j,k)*( &
                       ((dvksi*ksi_x(i,j,k)+dveta*eta_x(i,j,k)+dvphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                       +(dvksi*ksi_y(i,j,k)+dveta*eta_y(i,j,k)+dvphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                       +(dvksi*ksi_z(i,j,k)+dveta*eta_z(i,j,k)+dvphi*phi_z(i,j,k))*costeta(i,j,k)  &
                       )*ijacob3(i,j,k) + vf(i,j,k)*ir(i,j,k) )
             wt(i,j,k) = vg(i,j,k)*( &
                       ((dwksi*ksi_x(i,j,k)+dweta*eta_x(i,j,k)+dwphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                       +(dwksi*ksi_y(i,j,k)+dweta*eta_y(i,j,k)+dwphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                       +(dwksi*ksi_z(i,j,k)+dweta*eta_z(i,j,k)+dwphi*phi_z(i,j,k))*costeta(i,j,k)  &
                       )*ijacob3(i,j,k) + wf(i,j,k)*ir(i,j,k) )
             rt(i,j,k) = vg(i,j,k)*( &
                       ((drksi*ksi_x(i,j,k)+dreta*eta_x(i,j,k)+drphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                       +(drksi*ksi_y(i,j,k)+dreta*eta_y(i,j,k)+drphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                       +(drksi*ksi_z(i,j,k)+dreta*eta_z(i,j,k)+drphi*phi_z(i,j,k))*costeta(i,j,k)  &
                       )*ijacob3(i,j,k) + rf(i,j,k)*ir(i,j,k) )
          enddo
       enddo
    enddo

    ! Update fluxes at each RK step
    ! =============================
    do k=ndz,nfz
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

  end subroutine bc_TD3d_jmin_c3

  !===============================================================================
  module subroutine bc_TD3d_jmax_c3
  !===============================================================================
    !> 3D Tam & Dong's BC: boundary condition at jmax (top) - 3D curvilinear version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l
    real(wp) :: cp,av,inn1,ntm1
    real(wp) :: duksi,dvksi,dwksi,dpksi,drksi
    real(wp) :: dueta,dveta,dweta,dpeta,dreta
    real(wp) :: duphi,dvphi,dwphi,dpphi,drphi
    real(wp), dimension(nx1:nx2,-2:ngh,nz1:nz2) :: rf,uf,vf,wf,pf
    real(wp), dimension(nx,ngh,nz) :: rt,ut,vt,wt,pt,vg,c2_
    !-------------------------------------------------------------------------

    ! Pointers for spherical coordinates
    ! ==================================
    ir=>BC_face(2,2)%i_r
    cosphi=>BC_face(2,2)%cosp
    sinphi=>BC_face(2,2)%sinp
    costeta=>BC_face(2,2)%cost
    sinteta=>BC_face(2,2)%sint
    costcosp=>BC_face(2,2)%costcosp
    costsinp=>BC_face(2,2)%costsinp
    sintcosp=>BC_face(2,2)%sintcosp
    sintsinp=>BC_face(2,2)%sintsinp

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
    ! vg= u0*sin(teta)*cos(phi)+v0*sin(teta)*sin(phi)+w0*cos(teta)
    !   + sqrt{c0^2-[u0*cos(teta)*cos(phi)+v0*cos(teta)*sin(phi)-w0*sin(teta)]^2-[u0*sin(phi)-v0*cos(phi)]^2}
    do l=1,ngh
       do i=ndx,nfx
          do k=ndz,nfz
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

    ! Compute vg*[sin(teta)*cos(phi)*dq/dx+sin(teta)*sin(phi)*dq/dy+cos(teta)*dq/dz+q/r]
    ! ==================================================================================
    ! (Tam & Webb DRP schemes)
    do j=nymnghp1,ny-3
       l=j-nymngh
       do i=ndx,nfx
          do k=ndz,nfz
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

             duphi = a7(1)*( uf(i,l,k+1) - uf(i,l,k-1) ) + &
                     a7(2)*( uf(i,l,k+2) - uf(i,l,k-2) ) + &
                     a7(3)*( uf(i,l,k+3) - uf(i,l,k-3) )
             dvphi = a7(1)*( vf(i,l,k+1) - vf(i,l,k-1) ) + &
                     a7(2)*( vf(i,l,k+2) - vf(i,l,k-2) ) + &
                     a7(3)*( vf(i,l,k+3) - vf(i,l,k-3) )
             dwphi = a7(1)*( wf(i,l,k+1) - wf(i,l,k-1) ) + &
                     a7(2)*( wf(i,l,k+2) - wf(i,l,k-2) ) + &
                     a7(3)*( wf(i,l,k+3) - wf(i,l,k-3) )
             dpphi = a7(1)*( pf(i,l,k+1) - pf(i,l,k-1) ) + &
                     a7(2)*( pf(i,l,k+2) - pf(i,l,k-2) ) + &
                     a7(3)*( pf(i,l,k+3) - pf(i,l,k-3) )
             drphi = a7(1)*( rf(i,l,k+1) - rf(i,l,k-1) ) + &
                     a7(2)*( rf(i,l,k+2) - rf(i,l,k-2) ) + &
                     a7(3)*( rf(i,l,k+3) - rf(i,l,k-3) )

             pt(i,l,k) = vg(i,l,k)*( &
                       ((dpksi*ksi_x(i,j,k)+dpeta*eta_x(i,j,k)+dpphi*phi_x(i,j,k))*sintcosp(i,l,k) &
                       +(dpksi*ksi_y(i,j,k)+dpeta*eta_y(i,j,k)+dpphi*phi_y(i,j,k))*sintsinp(i,l,k) &
                       +(dpksi*ksi_z(i,j,k)+dpeta*eta_z(i,j,k)+dpphi*phi_z(i,j,k))*costeta(i,l,k)  &
                       )*ijacob3(i,j,k) + pf(i,l,k)*ir(i,l,k) )
             ut(i,l,k) = vg(i,l,k)*( &
                       ((duksi*ksi_x(i,j,k)+dueta*eta_x(i,j,k)+duphi*phi_x(i,j,k))*sintcosp(i,l,k) &
                       +(duksi*ksi_y(i,j,k)+dueta*eta_y(i,j,k)+duphi*phi_y(i,j,k))*sintsinp(i,l,k) &
                       +(duksi*ksi_z(i,j,k)+dueta*eta_z(i,j,k)+duphi*phi_z(i,j,k))*costeta(i,l,k)  &
                       )*ijacob3(i,j,k) + uf(i,l,k)*ir(i,l,k) )
             vt(i,l,k) = vg(i,l,k)*( &
                       ((dvksi*ksi_x(i,j,k)+dveta*eta_x(i,j,k)+dvphi*phi_x(i,j,k))*sintcosp(i,l,k) &
                       +(dvksi*ksi_y(i,j,k)+dveta*eta_y(i,j,k)+dvphi*phi_y(i,j,k))*sintsinp(i,l,k) &
                       +(dvksi*ksi_z(i,j,k)+dveta*eta_z(i,j,k)+dvphi*phi_z(i,j,k))*costeta(i,l,k)  &
                       )*ijacob3(i,j,k) + vf(i,l,k)*ir(i,l,k) )
             wt(i,l,k) = vg(i,l,k)*( &
                       ((dwksi*ksi_x(i,j,k)+dweta*eta_x(i,j,k)+dwphi*phi_x(i,j,k))*sintcosp(i,l,k) &
                       +(dwksi*ksi_y(i,j,k)+dweta*eta_y(i,j,k)+dwphi*phi_y(i,j,k))*sintsinp(i,l,k) &
                       +(dwksi*ksi_z(i,j,k)+dweta*eta_z(i,j,k)+dwphi*phi_z(i,j,k))*costeta(i,l,k)  &
                       )*ijacob3(i,j,k) + wf(i,l,k)*ir(i,l,k) )
             rt(i,l,k) = vg(i,l,k)*( &
                       ((drksi*ksi_x(i,j,k)+dreta*eta_x(i,j,k)+drphi*phi_x(i,j,k))*sintcosp(i,l,k) &
                       +(drksi*ksi_y(i,j,k)+dreta*eta_y(i,j,k)+drphi*phi_y(i,j,k))*sintsinp(i,l,k) &
                       +(drksi*ksi_z(i,j,k)+dreta*eta_z(i,j,k)+drphi*phi_z(i,j,k))*costeta(i,l,k)  &
                       )*ijacob3(i,j,k) + rf(i,l,k)*ir(i,l,k) )
          enddo
       enddo
    enddo

    j=ny-2
    l=j-nymngh
    do i=ndx,nfx
       do k=ndz,nfz
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

          duphi = a7(1)*( uf(i,l,k+1) - uf(i,l,k-1) ) + &
                  a7(2)*( uf(i,l,k+2) - uf(i,l,k-2) ) + &
                  a7(3)*( uf(i,l,k+3) - uf(i,l,k-3) )
          dvphi = a7(1)*( vf(i,l,k+1) - vf(i,l,k-1) ) + &
                  a7(2)*( vf(i,l,k+2) - vf(i,l,k-2) ) + &
                  a7(3)*( vf(i,l,k+3) - vf(i,l,k-3) )
          dwphi = a7(1)*( wf(i,l,k+1) - wf(i,l,k-1) ) + &
                  a7(2)*( wf(i,l,k+2) - wf(i,l,k-2) ) + &
                  a7(3)*( wf(i,l,k+3) - wf(i,l,k-3) )
          dpphi = a7(1)*( pf(i,l,k+1) - pf(i,l,k-1) ) + &
                  a7(2)*( pf(i,l,k+2) - pf(i,l,k-2) ) + &
                  a7(3)*( pf(i,l,k+3) - pf(i,l,k-3) )
          drphi = a7(1)*( rf(i,l,k+1) - rf(i,l,k-1) ) + &
                  a7(2)*( rf(i,l,k+2) - rf(i,l,k-2) ) + &
                  a7(3)*( rf(i,l,k+3) - rf(i,l,k-3) )

          pt(i,l,k) = vg(i,l,k)*( &
                    ((dpksi*ksi_x(i,j,k)+dpeta*eta_x(i,j,k)+dpphi*phi_x(i,j,k))*sintcosp(i,l,k) &
                    +(dpksi*ksi_y(i,j,k)+dpeta*eta_y(i,j,k)+dpphi*phi_y(i,j,k))*sintsinp(i,l,k) &
                    +(dpksi*ksi_z(i,j,k)+dpeta*eta_z(i,j,k)+dpphi*phi_z(i,j,k))*costeta(i,l,k)  &
                    )*ijacob3(i,j,k) + pf(i,l,k)*ir(i,l,k) )
          ut(i,l,k) = vg(i,l,k)*( &
                    ((duksi*ksi_x(i,j,k)+dueta*eta_x(i,j,k)+duphi*phi_x(i,j,k))*sintcosp(i,l,k) &
                    +(duksi*ksi_y(i,j,k)+dueta*eta_y(i,j,k)+duphi*phi_y(i,j,k))*sintsinp(i,l,k) &
                    +(duksi*ksi_z(i,j,k)+dueta*eta_z(i,j,k)+duphi*phi_z(i,j,k))*costeta(i,l,k)  &
                    )*ijacob3(i,j,k) + uf(i,l,k)*ir(i,l,k) )
          vt(i,l,k) = vg(i,l,k)*( &
                    ((dvksi*ksi_x(i,j,k)+dveta*eta_x(i,j,k)+dvphi*phi_x(i,j,k))*sintcosp(i,l,k) &
                    +(dvksi*ksi_y(i,j,k)+dveta*eta_y(i,j,k)+dvphi*phi_y(i,j,k))*sintsinp(i,l,k) &
                    +(dvksi*ksi_z(i,j,k)+dveta*eta_z(i,j,k)+dvphi*phi_z(i,j,k))*costeta(i,l,k)  &
                    )*ijacob3(i,j,k) + vf(i,l,k)*ir(i,l,k) )
          wt(i,l,k) = vg(i,l,k)*( &
                    ((dwksi*ksi_x(i,j,k)+dweta*eta_x(i,j,k)+dwphi*phi_x(i,j,k))*sintcosp(i,l,k) &
                    +(dwksi*ksi_y(i,j,k)+dweta*eta_y(i,j,k)+dwphi*phi_y(i,j,k))*sintsinp(i,l,k) &
                    +(dwksi*ksi_z(i,j,k)+dweta*eta_z(i,j,k)+dwphi*phi_z(i,j,k))*costeta(i,l,k)  &
                    )*ijacob3(i,j,k) + wf(i,l,k)*ir(i,l,k) )
          rt(i,l,k) = vg(i,l,k)*( &
                    ((drksi*ksi_x(i,j,k)+dreta*eta_x(i,j,k)+drphi*phi_x(i,j,k))*sintcosp(i,l,k) &
                    +(drksi*ksi_y(i,j,k)+dreta*eta_y(i,j,k)+drphi*phi_y(i,j,k))*sintsinp(i,l,k) &
                    +(drksi*ksi_z(i,j,k)+dreta*eta_z(i,j,k)+drphi*phi_z(i,j,k))*costeta(i,l,k)  &
                    )*ijacob3(i,j,k) + rf(i,l,k)*ir(i,l,k) )
       enddo
    enddo

    j=ny-1
    l=j-nymngh
    do i=ndx,nfx
       do k=ndz,nfz
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

          duphi = a7(1)*( uf(i,l,k+1) - uf(i,l,k-1) ) + &
                  a7(2)*( uf(i,l,k+2) - uf(i,l,k-2) ) + &
                  a7(3)*( uf(i,l,k+3) - uf(i,l,k-3) )
          dvphi = a7(1)*( vf(i,l,k+1) - vf(i,l,k-1) ) + &
                  a7(2)*( vf(i,l,k+2) - vf(i,l,k-2) ) + &
                  a7(3)*( vf(i,l,k+3) - vf(i,l,k-3) )
          dwphi = a7(1)*( wf(i,l,k+1) - wf(i,l,k-1) ) + &
                  a7(2)*( wf(i,l,k+2) - wf(i,l,k-2) ) + &
                  a7(3)*( wf(i,l,k+3) - wf(i,l,k-3) )
          dpphi = a7(1)*( pf(i,l,k+1) - pf(i,l,k-1) ) + &
                  a7(2)*( pf(i,l,k+2) - pf(i,l,k-2) ) + &
                  a7(3)*( pf(i,l,k+3) - pf(i,l,k-3) )
          drphi = a7(1)*( rf(i,l,k+1) - rf(i,l,k-1) ) + &
                  a7(2)*( rf(i,l,k+2) - rf(i,l,k-2) ) + &
                  a7(3)*( rf(i,l,k+3) - rf(i,l,k-3) )

          pt(i,l,k) = vg(i,l,k)*( &
                    ((dpksi*ksi_x(i,j,k)+dpeta*eta_x(i,j,k)+dpphi*phi_x(i,j,k))*sintcosp(i,l,k) &
                    +(dpksi*ksi_y(i,j,k)+dpeta*eta_y(i,j,k)+dpphi*phi_y(i,j,k))*sintsinp(i,l,k) &
                    +(dpksi*ksi_z(i,j,k)+dpeta*eta_z(i,j,k)+dpphi*phi_z(i,j,k))*costeta(i,l,k)  &
                    )*ijacob3(i,j,k) + pf(i,l,k)*ir(i,l,k) )
          ut(i,l,k) = vg(i,l,k)*( &
                    ((duksi*ksi_x(i,j,k)+dueta*eta_x(i,j,k)+duphi*phi_x(i,j,k))*sintcosp(i,l,k) &
                    +(duksi*ksi_y(i,j,k)+dueta*eta_y(i,j,k)+duphi*phi_y(i,j,k))*sintsinp(i,l,k) &
                    +(duksi*ksi_z(i,j,k)+dueta*eta_z(i,j,k)+duphi*phi_z(i,j,k))*costeta(i,l,k)  &
                    )*ijacob3(i,j,k) + uf(i,l,k)*ir(i,l,k) )
          vt(i,l,k) = vg(i,l,k)*( &
                    ((dvksi*ksi_x(i,j,k)+dveta*eta_x(i,j,k)+dvphi*phi_x(i,j,k))*sintcosp(i,l,k) &
                    +(dvksi*ksi_y(i,j,k)+dveta*eta_y(i,j,k)+dvphi*phi_y(i,j,k))*sintsinp(i,l,k) &
                    +(dvksi*ksi_z(i,j,k)+dveta*eta_z(i,j,k)+dvphi*phi_z(i,j,k))*costeta(i,l,k)  &
                    )*ijacob3(i,j,k) + vf(i,l,k)*ir(i,l,k) )
          wt(i,l,k) = vg(i,l,k)*( &
                    ((dwksi*ksi_x(i,j,k)+dweta*eta_x(i,j,k)+dwphi*phi_x(i,j,k))*sintcosp(i,l,k) &
                    +(dwksi*ksi_y(i,j,k)+dweta*eta_y(i,j,k)+dwphi*phi_y(i,j,k))*sintsinp(i,l,k) &
                    +(dwksi*ksi_z(i,j,k)+dweta*eta_z(i,j,k)+dwphi*phi_z(i,j,k))*costeta(i,l,k)  &
                    )*ijacob3(i,j,k) + wf(i,l,k)*ir(i,l,k) )
          rt(i,l,k) = vg(i,l,k)*( &
                    ((drksi*ksi_x(i,j,k)+dreta*eta_x(i,j,k)+drphi*phi_x(i,j,k))*sintcosp(i,l,k) &
                    +(drksi*ksi_y(i,j,k)+dreta*eta_y(i,j,k)+drphi*phi_y(i,j,k))*sintsinp(i,l,k) &
                    +(drksi*ksi_z(i,j,k)+dreta*eta_z(i,j,k)+drphi*phi_z(i,j,k))*costeta(i,l,k)  &
                    )*ijacob3(i,j,k) + rf(i,l,k)*ir(i,l,k) )
       enddo
    enddo

    j=ny
    l=j-nymngh
    do i=ndx,nfx
       do k=ndz,nfz
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

          duphi = a7(1)*( uf(i,l,k+1) - uf(i,l,k-1) ) + &
                  a7(2)*( uf(i,l,k+2) - uf(i,l,k-2) ) + &
                  a7(3)*( uf(i,l,k+3) - uf(i,l,k-3) )
          dvphi = a7(1)*( vf(i,l,k+1) - vf(i,l,k-1) ) + &
                  a7(2)*( vf(i,l,k+2) - vf(i,l,k-2) ) + &
                  a7(3)*( vf(i,l,k+3) - vf(i,l,k-3) )
          dwphi = a7(1)*( wf(i,l,k+1) - wf(i,l,k-1) ) + &
                  a7(2)*( wf(i,l,k+2) - wf(i,l,k-2) ) + &
                  a7(3)*( wf(i,l,k+3) - wf(i,l,k-3) )
          dpphi = a7(1)*( pf(i,l,k+1) - pf(i,l,k-1) ) + &
                  a7(2)*( pf(i,l,k+2) - pf(i,l,k-2) ) + &
                  a7(3)*( pf(i,l,k+3) - pf(i,l,k-3) )
          drphi = a7(1)*( rf(i,l,k+1) - rf(i,l,k-1) ) + &
                  a7(2)*( rf(i,l,k+2) - rf(i,l,k-2) ) + &
                  a7(3)*( rf(i,l,k+3) - rf(i,l,k-3) )

          pt(i,l,k) = vg(i,l,k)*( &
                    ((dpksi*ksi_x(i,j,k)+dpeta*eta_x(i,j,k)+dpphi*phi_x(i,j,k))*sintcosp(i,l,k) &
                    +(dpksi*ksi_y(i,j,k)+dpeta*eta_y(i,j,k)+dpphi*phi_y(i,j,k))*sintsinp(i,l,k) &
                    +(dpksi*ksi_z(i,j,k)+dpeta*eta_z(i,j,k)+dpphi*phi_z(i,j,k))*costeta(i,l,k)  &
                    )*ijacob3(i,j,k) + pf(i,l,k)*ir(i,l,k) )
          ut(i,l,k) = vg(i,l,k)*( &
                    ((duksi*ksi_x(i,j,k)+dueta*eta_x(i,j,k)+duphi*phi_x(i,j,k))*sintcosp(i,l,k) &
                    +(duksi*ksi_y(i,j,k)+dueta*eta_y(i,j,k)+duphi*phi_y(i,j,k))*sintsinp(i,l,k) &
                    +(duksi*ksi_z(i,j,k)+dueta*eta_z(i,j,k)+duphi*phi_z(i,j,k))*costeta(i,l,k)  &
                    )*ijacob3(i,j,k) + uf(i,l,k)*ir(i,l,k) )
          vt(i,l,k) = vg(i,l,k)*( &
                    ((dvksi*ksi_x(i,j,k)+dveta*eta_x(i,j,k)+dvphi*phi_x(i,j,k))*sintcosp(i,l,k) &
                    +(dvksi*ksi_y(i,j,k)+dveta*eta_y(i,j,k)+dvphi*phi_y(i,j,k))*sintsinp(i,l,k) &
                    +(dvksi*ksi_z(i,j,k)+dveta*eta_z(i,j,k)+dvphi*phi_z(i,j,k))*costeta(i,l,k)  &
                    )*ijacob3(i,j,k) + vf(i,l,k)*ir(i,l,k) )
          wt(i,l,k) = vg(i,l,k)*( &
                    ((dwksi*ksi_x(i,j,k)+dweta*eta_x(i,j,k)+dwphi*phi_x(i,j,k))*sintcosp(i,l,k) &
                    +(dwksi*ksi_y(i,j,k)+dweta*eta_y(i,j,k)+dwphi*phi_y(i,j,k))*sintsinp(i,l,k) &
                    +(dwksi*ksi_z(i,j,k)+dweta*eta_z(i,j,k)+dwphi*phi_z(i,j,k))*costeta(i,l,k)  &
                    )*ijacob3(i,j,k) + wf(i,l,k)*ir(i,l,k) )
          rt(i,l,k) = vg(i,l,k)*( &
                    ((drksi*ksi_x(i,j,k)+dreta*eta_x(i,j,k)+drphi*phi_x(i,j,k))*sintcosp(i,l,k) &
                    +(drksi*ksi_y(i,j,k)+dreta*eta_y(i,j,k)+drphi*phi_y(i,j,k))*sintsinp(i,l,k) &
                    +(drksi*ksi_z(i,j,k)+dreta*eta_z(i,j,k)+drphi*phi_z(i,j,k))*costeta(i,l,k)  &
                    )*ijacob3(i,j,k) + rf(i,l,k)*ir(i,l,k) )
       enddo
    enddo

    ! Update fluxes at each RK step
    ! =============================
    do k=ndz,nfz
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

  end subroutine bc_TD3d_jmax_c3

  !===============================================================================
  module subroutine bc_TD3d_kmin_c3
  !===============================================================================
    !> 3D Tam & Dong's BC: boundary condition at kmin (front) - 3D curvilinear version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: cp,av,inn1,ntm1
    real(wp) :: duksi,dvksi,dwksi,dpksi,drksi
    real(wp) :: dueta,dveta,dweta,dpeta,dreta
    real(wp) :: duphi,dvphi,dwphi,dpphi,drphi
    real(wp), dimension(nx1:nx2,ny1:ny2,1:ngh+3) :: rf,uf,vf,wf,pf
    real(wp), dimension(nx,ny,ngh) :: rt,ut,vt,wt,pt,vg,c2_
    !-------------------------------------------------------------------------

    ! Pointers for spherical coordinates
    ! ==================================
    ir=>BC_face(3,1)%i_r
    cosphi=>BC_face(3,1)%cosp
    sinphi=>BC_face(3,1)%sinp
    costeta=>BC_face(3,1)%cost
    sinteta=>BC_face(3,1)%sint
    costcosp=>BC_face(3,1)%costcosp
    costsinp=>BC_face(3,1)%costsinp
    sintcosp=>BC_face(3,1)%sintcosp
    sintsinp=>BC_face(3,1)%sintsinp

    ! Sound speed
    ! ===========
    do k=1,ngh
       do i=1,nx
          do j=1,ny
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
       do j=ny1,ny2
          BC_face(3,1)%U0(:,j,1:nghp3,1)=(ntm1*BC_face(3,1)%U0(:,j,1:nghp3,1)+rho_n(:,j,1:nghp3))*inn1
          BC_face(3,1)%U0(:,j,1:nghp3,2)=(ntm1*BC_face(3,1)%U0(:,j,1:nghp3,2)+uu(:,j,1:nghp3))*inn1
          BC_face(3,1)%U0(:,j,1:nghp3,3)=(ntm1*BC_face(3,1)%U0(:,j,1:nghp3,3)+vv(:,j,1:nghp3))*inn1
          BC_face(3,1)%U0(:,j,1:nghp3,4)=(ntm1*BC_face(3,1)%U0(:,j,1:nghp3,4)+ww(:,j,1:nghp3))*inn1
          BC_face(3,1)%U0(:,j,1:nghp3,5)=(ntm1*BC_face(3,1)%U0(:,j,1:nghp3,5)+prs(:,j,1:nghp3))*inn1
       enddo
       ! time-averaged sound speed squared
       do j=1,ny
          BC_face(3,1)%U0(1:nx,j,1:ngh,6)=(ntm1*BC_face(3,1)%U0(1:nx,j,1:ngh,6)+c2_(1:nx,j,1:ngh))*inn1
       enddo
    endif

    ! Compute fluctuations
    ! ====================
    do k=1,nghp3
       do i=nx1,nx2
          do j=ny1,ny2
             rf(i,j,k)=rho_n(i,j,k)-BC_face(3,1)%U0(i,j,k,1)
             uf(i,j,k)=   uu(i,j,k)-BC_face(3,1)%U0(i,j,k,2)
             vf(i,j,k)=   vv(i,j,k)-BC_face(3,1)%U0(i,j,k,3)
             wf(i,j,k)=   ww(i,j,k)-BC_face(3,1)%U0(i,j,k,4)
             pf(i,j,k)=  prs(i,j,k)-BC_face(3,1)%U0(i,j,k,5)
          enddo
       enddo
    enddo

    ! Compute group velocity vg
    ! =========================
    ! vg= u0*sin(teta)*cos(phi)+v0*sin(teta)*sin(phi)+w0*cos(teta)
    !   + sqrt{c0^2-[u0*cos(teta)*cos(phi)+v0*cos(teta)*sin(phi)-w0*sin(teta)]^2-[u0*sin(phi)-v0*cos(phi)]^2}
    do k=1,ngh
       do i=ndx,nfx
          do j=ndy,nfy
             vg(i,j,k)= BC_face(3,1)%U0(i,j,k,2)*sintcosp(i,j,k) +     &
               BC_face(3,1)%U0(i,j,k,3)*sintsinp(i,j,k) +     &
               BC_face(3,1)%U0(i,j,k,4)*costeta(i,j,k)          &
                      + sqrt( BC_face(3,1)%U0(i,j,k,6)-                &
                       ( BC_face(3,1)%U0(i,j,k,2)*costcosp(i,j,k) +    &
                   BC_face(3,1)%U0(i,j,k,3)*costsinp(i,j,k) -    &
                     BC_face(3,1)%U0(i,j,k,4)*sinteta(i,j,k) )**2- &
                  ( BC_face(3,1)%U0(i,j,k,2)*sinphi(i,j,k) -      &
                     BC_face(3,1)%U0(i,j,k,3)*cosphi(i,j,k) )**2 )
          enddo
       enddo
    enddo

    ! Compute vg*[sin(teta)*cos(phi)*dq/dx+sin(teta)*sin(phi)*dq/dy+cos(teta)*dq/dz+q/r]
    ! ==================================================================================
    ! (Tam & Webb DRP schemes)
    k=1
    do i=ndx,nfx
       do j=ndy,nfy
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

          duphi = a06(1)*uf(i,j,1)+a06(2)*uf(i,j,2) &
                + a06(3)*uf(i,j,3)+a06(4)*uf(i,j,4) &
                + a06(5)*uf(i,j,5)+a06(6)*uf(i,j,6) &
                + a06(7)*uf(i,j,7)
          dvphi = a06(1)*vf(i,j,1)+a06(2)*vf(i,j,2) &
                + a06(3)*vf(i,j,3)+a06(4)*vf(i,j,4) &
                + a06(5)*vf(i,j,5)+a06(6)*vf(i,j,6) &
                + a06(7)*vf(i,j,7)
          dwphi = a06(1)*wf(i,j,1)+a06(2)*wf(i,j,2) &
                + a06(3)*wf(i,j,3)+a06(4)*wf(i,j,4) &
                + a06(5)*wf(i,j,5)+a06(6)*wf(i,j,6) &
                + a06(7)*wf(i,j,7)
          dpphi = a06(1)*pf(i,j,1)+a06(2)*pf(i,j,2) &
                + a06(3)*pf(i,j,3)+a06(4)*pf(i,j,4) &
                + a06(5)*pf(i,j,5)+a06(6)*pf(i,j,6) &
                + a06(7)*pf(i,j,7)
          drphi = a06(1)*rf(i,j,1)+a06(2)*rf(i,j,2) &
                + a06(3)*rf(i,j,3)+a06(4)*rf(i,j,4) &
                + a06(5)*rf(i,j,5)+a06(6)*rf(i,j,6) &
                + a06(7)*rf(i,j,7)

          pt(i,j,k) = vg(i,j,k)*( &
                    ((dpksi*ksi_x(i,j,k)+dpeta*eta_x(i,j,k)+dpphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(dpksi*ksi_y(i,j,k)+dpeta*eta_y(i,j,k)+dpphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(dpksi*ksi_z(i,j,k)+dpeta*eta_z(i,j,k)+dpphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + pf(i,j,k)*ir(i,j,k) )
          ut(i,j,k) = vg(i,j,k)*( &
                    ((duksi*ksi_x(i,j,k)+dueta*eta_x(i,j,k)+duphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(duksi*ksi_y(i,j,k)+dueta*eta_y(i,j,k)+duphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(duksi*ksi_z(i,j,k)+dueta*eta_z(i,j,k)+duphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + uf(i,j,k)*ir(i,j,k) )
          vt(i,j,k) = vg(i,j,k)*( &
                    ((dvksi*ksi_x(i,j,k)+dveta*eta_x(i,j,k)+dvphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(dvksi*ksi_y(i,j,k)+dveta*eta_y(i,j,k)+dvphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(dvksi*ksi_z(i,j,k)+dveta*eta_z(i,j,k)+dvphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + vf(i,j,k)*ir(i,j,k) )
          wt(i,j,k) = vg(i,j,k)*( &
                    ((dwksi*ksi_x(i,j,k)+dweta*eta_x(i,j,k)+dwphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(dwksi*ksi_y(i,j,k)+dweta*eta_y(i,j,k)+dwphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(dwksi*ksi_z(i,j,k)+dweta*eta_z(i,j,k)+dwphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + wf(i,j,k)*ir(i,j,k) )
          rt(i,j,k) = vg(i,j,k)*( &
                    ((drksi*ksi_x(i,j,k)+dreta*eta_x(i,j,k)+drphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(drksi*ksi_y(i,j,k)+dreta*eta_y(i,j,k)+drphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(drksi*ksi_z(i,j,k)+dreta*eta_z(i,j,k)+drphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + rf(i,j,k)*ir(i,j,k) )
       enddo
    enddo

    k=2
    do i=ndx,nfx
       do j=ndy,nfy
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

          duphi = a15(1)*uf(i,j,1)+a15(2)*uf(i,j,2) &
                + a15(3)*uf(i,j,3)+a15(4)*uf(i,j,4) &
                + a15(5)*uf(i,j,5)+a15(6)*uf(i,j,6) &
                + a15(7)*uf(i,j,7)
          dvphi = a15(1)*vf(i,j,1)+a15(2)*vf(i,j,2) &
                + a15(3)*vf(i,j,3)+a15(4)*vf(i,j,4) &
                + a15(5)*vf(i,j,5)+a15(6)*vf(i,j,6) &
                + a15(7)*vf(i,j,7)
          dwphi = a15(1)*wf(i,j,1)+a15(2)*wf(i,j,2) &
                + a15(3)*wf(i,j,3)+a15(4)*wf(i,j,4) &
                + a15(5)*wf(i,j,5)+a15(6)*wf(i,j,6) &
                + a15(7)*wf(i,j,7)
          dpphi = a15(1)*pf(i,j,1)+a15(2)*pf(i,j,2) &
                + a15(3)*pf(i,j,3)+a15(4)*pf(i,j,4) &
                + a15(5)*pf(i,j,5)+a15(6)*pf(i,j,6) &
                + a15(7)*pf(i,j,7)
          drphi = a15(1)*rf(i,j,1)+a15(2)*rf(i,j,2) &
                + a15(3)*rf(i,j,3)+a15(4)*rf(i,j,4) &
                + a15(5)*rf(i,j,5)+a15(6)*rf(i,j,6) &
                + a15(7)*rf(i,j,7)

          pt(i,j,k) = vg(i,j,k)*( &
                    ((dpksi*ksi_x(i,j,k)+dpeta*eta_x(i,j,k)+dpphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(dpksi*ksi_y(i,j,k)+dpeta*eta_y(i,j,k)+dpphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(dpksi*ksi_z(i,j,k)+dpeta*eta_z(i,j,k)+dpphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + pf(i,j,k)*ir(i,j,k) )
          ut(i,j,k) = vg(i,j,k)*( &
                    ((duksi*ksi_x(i,j,k)+dueta*eta_x(i,j,k)+duphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(duksi*ksi_y(i,j,k)+dueta*eta_y(i,j,k)+duphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(duksi*ksi_z(i,j,k)+dueta*eta_z(i,j,k)+duphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + uf(i,j,k)*ir(i,j,k) )
          vt(i,j,k) = vg(i,j,k)*( &
                    ((dvksi*ksi_x(i,j,k)+dveta*eta_x(i,j,k)+dvphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(dvksi*ksi_y(i,j,k)+dveta*eta_y(i,j,k)+dvphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(dvksi*ksi_z(i,j,k)+dveta*eta_z(i,j,k)+dvphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + vf(i,j,k)*ir(i,j,k) )
          wt(i,j,k) = vg(i,j,k)*( &
                    ((dwksi*ksi_x(i,j,k)+dweta*eta_x(i,j,k)+dwphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(dwksi*ksi_y(i,j,k)+dweta*eta_y(i,j,k)+dwphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(dwksi*ksi_z(i,j,k)+dweta*eta_z(i,j,k)+dwphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + wf(i,j,k)*ir(i,j,k) )
          rt(i,j,k) = vg(i,j,k)*( &
                    ((drksi*ksi_x(i,j,k)+dreta*eta_x(i,j,k)+drphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(drksi*ksi_y(i,j,k)+dreta*eta_y(i,j,k)+drphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(drksi*ksi_z(i,j,k)+dreta*eta_z(i,j,k)+drphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + rf(i,j,k)*ir(i,j,k) )
       enddo
    enddo

    k=3
    do i=ndx,nfx
       do j=ndy,nfy
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

          duphi = a24(1)*uf(i,j,1)+a24(2)*uf(i,j,2) &
                + a24(3)*uf(i,j,3)+a24(4)*uf(i,j,4) &
                + a24(5)*uf(i,j,5)+a24(6)*uf(i,j,6) &
                + a24(7)*uf(i,j,7)
          dvphi = a24(1)*vf(i,j,1)+a24(2)*vf(i,j,2) &
                + a24(3)*vf(i,j,3)+a24(4)*vf(i,j,4) &
                + a24(5)*vf(i,j,5)+a24(6)*vf(i,j,6) &
                + a24(7)*vf(i,j,7)
          dwphi = a24(1)*wf(i,j,1)+a24(2)*wf(i,j,2) &
                + a24(3)*wf(i,j,3)+a24(4)*wf(i,j,4) &
                + a24(5)*wf(i,j,5)+a24(6)*wf(i,j,6) &
                + a24(7)*wf(i,j,7)
          dpphi = a24(1)*pf(i,j,1)+a24(2)*pf(i,j,2) &
                + a24(3)*pf(i,j,3)+a24(4)*pf(i,j,4) &
                + a24(5)*pf(i,j,5)+a24(6)*pf(i,j,6) &
                + a24(7)*pf(i,j,7)
          drphi = a24(1)*rf(i,j,1)+a24(2)*rf(i,j,2) &
                + a24(3)*rf(i,j,3)+a24(4)*rf(i,j,4) &
                + a24(5)*rf(i,j,5)+a24(6)*rf(i,j,6) &
                + a24(7)*rf(i,j,7)

          pt(i,j,k) = vg(i,j,k)*( &
                    ((dpksi*ksi_x(i,j,k)+dpeta*eta_x(i,j,k)+dpphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(dpksi*ksi_y(i,j,k)+dpeta*eta_y(i,j,k)+dpphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(dpksi*ksi_z(i,j,k)+dpeta*eta_z(i,j,k)+dpphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + pf(i,j,k)*ir(i,j,k) )
          ut(i,j,k) = vg(i,j,k)*( &
                    ((duksi*ksi_x(i,j,k)+dueta*eta_x(i,j,k)+duphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(duksi*ksi_y(i,j,k)+dueta*eta_y(i,j,k)+duphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(duksi*ksi_z(i,j,k)+dueta*eta_z(i,j,k)+duphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + uf(i,j,k)*ir(i,j,k) )
          vt(i,j,k) = vg(i,j,k)*( &
                    ((dvksi*ksi_x(i,j,k)+dveta*eta_x(i,j,k)+dvphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(dvksi*ksi_y(i,j,k)+dveta*eta_y(i,j,k)+dvphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(dvksi*ksi_z(i,j,k)+dveta*eta_z(i,j,k)+dvphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + vf(i,j,k)*ir(i,j,k) )
          wt(i,j,k) = vg(i,j,k)*( &
                    ((dwksi*ksi_x(i,j,k)+dweta*eta_x(i,j,k)+dwphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(dwksi*ksi_y(i,j,k)+dweta*eta_y(i,j,k)+dwphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(dwksi*ksi_z(i,j,k)+dweta*eta_z(i,j,k)+dwphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + wf(i,j,k)*ir(i,j,k) )
          rt(i,j,k) = vg(i,j,k)*( &
                    ((drksi*ksi_x(i,j,k)+dreta*eta_x(i,j,k)+drphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                    +(drksi*ksi_y(i,j,k)+dreta*eta_y(i,j,k)+drphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                    +(drksi*ksi_z(i,j,k)+dreta*eta_z(i,j,k)+drphi*phi_z(i,j,k))*costeta(i,j,k)  &
                    )*ijacob3(i,j,k) + rf(i,j,k)*ir(i,j,k) )
       enddo
    enddo

    do k=4,ngh
       do i=ndx,nfx
          do j=ndy,nfy
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

             duphi = a7(1)*( uf(i,j,k+1) - uf(i,j,k-1) ) + &
                     a7(2)*( uf(i,j,k+2) - uf(i,j,k-2) ) + &
                     a7(3)*( uf(i,j,k+3) - uf(i,j,k-3) )
             dvphi = a7(1)*( vf(i,j,k+1) - vf(i,j,k-1) ) + &
                     a7(2)*( vf(i,j,k+2) - vf(i,j,k-2) ) + &
                     a7(3)*( vf(i,j,k+3) - vf(i,j,k-3) )
             dwphi = a7(1)*( wf(i,j,k+1) - wf(i,j,k-1) ) + &
                     a7(2)*( wf(i,j,k+2) - wf(i,j,k-2) ) + &
                     a7(3)*( wf(i,j,k+3) - wf(i,j,k-3) )
             dpphi = a7(1)*( pf(i,j,k+1) - pf(i,j,k-1) ) + &
                     a7(2)*( pf(i,j,k+2) - pf(i,j,k-2) ) + &
                     a7(3)*( pf(i,j,k+3) - pf(i,j,k-3) )
             drphi = a7(1)*( rf(i,j,k+1) - rf(i,j,k-1) ) + &
                     a7(2)*( rf(i,j,k+2) - rf(i,j,k-2) ) + &
                     a7(3)*( rf(i,j,k+3) - rf(i,j,k-3) )

             pt(i,j,k) = vg(i,j,k)*( &
                       ((dpksi*ksi_x(i,j,k)+dpeta*eta_x(i,j,k)+dpphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                       +(dpksi*ksi_y(i,j,k)+dpeta*eta_y(i,j,k)+dpphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                       +(dpksi*ksi_z(i,j,k)+dpeta*eta_z(i,j,k)+dpphi*phi_z(i,j,k))*costeta(i,j,k)  &
                       )*ijacob3(i,j,k) + pf(i,j,k)*ir(i,j,k) )
             ut(i,j,k) = vg(i,j,k)*( &
                       ((duksi*ksi_x(i,j,k)+dueta*eta_x(i,j,k)+duphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                       +(duksi*ksi_y(i,j,k)+dueta*eta_y(i,j,k)+duphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                       +(duksi*ksi_z(i,j,k)+dueta*eta_z(i,j,k)+duphi*phi_z(i,j,k))*costeta(i,j,k)  &
                       )*ijacob3(i,j,k) + uf(i,j,k)*ir(i,j,k) )
             vt(i,j,k) = vg(i,j,k)*( &
                       ((dvksi*ksi_x(i,j,k)+dveta*eta_x(i,j,k)+dvphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                       +(dvksi*ksi_y(i,j,k)+dveta*eta_y(i,j,k)+dvphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                       +(dvksi*ksi_z(i,j,k)+dveta*eta_z(i,j,k)+dvphi*phi_z(i,j,k))*costeta(i,j,k)  &
                       )*ijacob3(i,j,k) + vf(i,j,k)*ir(i,j,k) )
             wt(i,j,k) = vg(i,j,k)*( &
                       ((dwksi*ksi_x(i,j,k)+dweta*eta_x(i,j,k)+dwphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                       +(dwksi*ksi_y(i,j,k)+dweta*eta_y(i,j,k)+dwphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                       +(dwksi*ksi_z(i,j,k)+dweta*eta_z(i,j,k)+dwphi*phi_z(i,j,k))*costeta(i,j,k)  &
                       )*ijacob3(i,j,k) + wf(i,j,k)*ir(i,j,k) )
             rt(i,j,k) = vg(i,j,k)*( &
                       ((drksi*ksi_x(i,j,k)+dreta*eta_x(i,j,k)+drphi*phi_x(i,j,k))*sintcosp(i,j,k) &
                       +(drksi*ksi_y(i,j,k)+dreta*eta_y(i,j,k)+drphi*phi_y(i,j,k))*sintsinp(i,j,k) &
                       +(drksi*ksi_z(i,j,k)+dreta*eta_z(i,j,k)+drphi*phi_z(i,j,k))*costeta(i,j,k)  &
                       )*ijacob3(i,j,k) + rf(i,j,k)*ir(i,j,k) )
          enddo
       enddo
    enddo

    ! Update fluxes at each RK step
    ! =============================
    do k=1,ngh
       do j=ndy,nfy
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

  end subroutine bc_TD3d_kmin_c3

  !===============================================================================
  module subroutine bc_TD3d_kmax_c3
  !===============================================================================
    !> 3D Tam & Dong's BC: boundary condition at kmax (back) - 3D curvilinear version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l
    real(wp) :: cp,av,inn1,ntm1
    real(wp) :: duksi,dvksi,dwksi,dpksi,drksi
    real(wp) :: dueta,dveta,dweta,dpeta,dreta
    real(wp) :: duphi,dvphi,dwphi,dpphi,drphi
    real(wp), dimension(nx1:nx2,ny1:ny2,-2:ngh) :: rf,uf,vf,wf,pf
    real(wp), dimension(nx,ny,ngh) :: rt,ut,vt,wt,pt,vg,c2_
    !-------------------------------------------------------------------------

    ! Pointers for spherical coordinates
    ! ==================================
    ir=>BC_face(3,2)%i_r
    cosphi=>BC_face(3,2)%cosp
    sinphi=>BC_face(3,2)%sinp
    costeta=>BC_face(3,2)%cost
    sinteta=>BC_face(3,2)%sint
    costcosp=>BC_face(3,2)%costcosp
    costsinp=>BC_face(3,2)%costsinp
    sintcosp=>BC_face(3,2)%sintcosp
    sintsinp=>BC_face(3,2)%sintsinp

    ! Sound speed
    ! ===========
    do k=nzmnghp1,nz
       l=k-nzmngh
       do i=1,nx
          do j=1,ny
             c2_(i,j,l)=c2calc_tro(Tmp(i,j,k),rho_n(i,j,k))
          enddo
       enddo
    enddo

    ! Compute online time-averaged primitive variables
    ! ================================================
    if (irk==nrk) then
       inn1 = 1.0_wp/dble(ntotal)
       ntm1=dble(ntotal-1)
       ! time-averaged primitive variables
       do j=ny1,ny2
          BC_face(3,2)%U0(:,j,-2:ngh,1)=(ntm1*BC_face(3,2)%U0(:,j,-2:ngh,1)+rho_n(:,j,nzmngh-2:nz))*inn1
          BC_face(3,2)%U0(:,j,-2:ngh,2)=(ntm1*BC_face(3,2)%U0(:,j,-2:ngh,2)+uu(:,j,nzmngh-2:nz))*inn1
          BC_face(3,2)%U0(:,j,-2:ngh,3)=(ntm1*BC_face(3,2)%U0(:,j,-2:ngh,3)+vv(:,j,nzmngh-2:nz))*inn1
          BC_face(3,2)%U0(:,j,-2:ngh,4)=(ntm1*BC_face(3,2)%U0(:,j,-2:ngh,4)+ww(:,j,nzmngh-2:nz))*inn1
          BC_face(3,2)%U0(:,j,-2:ngh,5)=(ntm1*BC_face(3,2)%U0(:,j,-2:ngh,5)+prs(:,j,nzmngh-2:nz))*inn1
       enddo
       ! time-averaged sound speed squared
       do j=1,ny
          BC_face(3,2)%U0(1:nx,j,1:ngh,6)=(ntm1*BC_face(3,2)%U0(1:nx,j,1:ngh,6)+c2_(1:nx,j,1:ngh))*inn1
       enddo
    endif

    ! Compute fluctuations
    ! ====================
    do k=nzmngh-2,nz
       l=k-nzmngh
       do i=nx1,nx2
          do j=ny1,ny2
             rf(i,j,l)=rho_n(i,j,k)-BC_face(3,2)%U0(i,j,l,1)
             uf(i,j,l)=   uu(i,j,k)-BC_face(3,2)%U0(i,j,l,2)
             vf(i,j,l)=   vv(i,j,k)-BC_face(3,2)%U0(i,j,l,3)
             wf(i,j,l)=   ww(i,j,k)-BC_face(3,2)%U0(i,j,l,4)
             pf(i,j,l)=  prs(i,j,k)-BC_face(3,2)%U0(i,j,l,5)
          enddo
       enddo
    enddo

    ! Compute group velocity vg
    ! =========================
    ! vg= u0*sin(teta)*cos(phi)+v0*sin(teta)*sin(phi)+w0*cos(teta)
    !   + sqrt{c0^2-[u0*cos(teta)*cos(phi)+v0*cos(teta)*sin(phi)-w0*sin(teta)]^2-[u0*sin(phi)-v0*cos(phi)]^2}
    do l=1,ngh
       do i=ndx,nfx
          do j=ndy,nfy
             vg(i,j,l)= BC_face(3,2)%U0(i,j,l,2)*sintcosp(i,j,l) +     &
               BC_face(3,2)%U0(i,j,l,3)*sintsinp(i,j,l) +     &
               BC_face(3,2)%U0(i,j,l,4)*costeta(i,j,l)          &
                      + sqrt( BC_face(3,2)%U0(i,j,l,6)-                &
                       ( BC_face(3,2)%U0(i,j,l,2)*costcosp(i,j,l) +    &
                   BC_face(3,2)%U0(i,j,l,3)*costsinp(i,j,l) -    &
                     BC_face(3,2)%U0(i,j,l,4)*sinteta(i,j,l) )**2- &
                  ( BC_face(3,2)%U0(i,j,l,2)*sinphi(i,j,l) -      &
                     BC_face(3,2)%U0(i,j,l,3)*cosphi(i,j,l) )**2 )
          enddo
       enddo
    enddo

    ! Compute vg*[sin(teta)*cos(phi)*dq/dx+sin(teta)*sin(phi)*dq/dy+cos(teta)*dq/dz+q/r]
    ! ==================================================================================
    ! (Tam & Webb DRP schemes)
    do k=nzmnghp1,nz-3
       l=k-nzmngh
       do i=ndx,nfx
          do j=ndy,nfy
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

             dueta = a7(1)*( uf(i,j+1,l) - uf(i,j-1,l) ) + &
                     a7(2)*( uf(i,j+2,l) - uf(i,j-2,l) ) + &
                     a7(3)*( uf(i,j+3,l) - uf(i,j-3,l) )
             dveta = a7(1)*( vf(i,j+1,l) - vf(i,j-1,l) ) + &
                     a7(2)*( vf(i,j+2,l) - vf(i,j-2,l) ) + &
                     a7(3)*( vf(i,j+3,l) - vf(i,j-3,l) )
             dweta = a7(1)*( wf(i,j+1,l) - wf(i,j-1,l) ) + &
                     a7(2)*( wf(i,j+2,l) - wf(i,j-2,l) ) + &
                     a7(3)*( wf(i,j+3,l) - wf(i,j-3,l) )
             dpeta = a7(1)*( pf(i,j+1,l) - pf(i,j-1,l) ) + &
                     a7(2)*( pf(i,j+2,l) - pf(i,j-2,l) ) + &
                     a7(3)*( pf(i,j+3,l) - pf(i,j-3,l) )
             dreta = a7(1)*( rf(i,j+1,l) - rf(i,j-1,l) ) + &
                     a7(2)*( rf(i,j+2,l) - rf(i,j-2,l) ) + &
                     a7(3)*( rf(i,j+3,l) - rf(i,j-3,l) )

             duphi = a7(1)*( uf(i,j,l+1) - uf(i,j,l-1) ) + &
                     a7(2)*( uf(i,j,l+2) - uf(i,j,l-2) ) + &
                     a7(3)*( uf(i,j,l+3) - uf(i,j,l-3) )
             dvphi = a7(1)*( vf(i,j,l+1) - vf(i,j,l-1) ) + &
                     a7(2)*( vf(i,j,l+2) - vf(i,j,l-2) ) + &
                     a7(3)*( vf(i,j,l+3) - vf(i,j,l-3) )
             dwphi = a7(1)*( wf(i,j,l+1) - wf(i,j,l-1) ) + &
                     a7(2)*( wf(i,j,l+2) - wf(i,j,l-2) ) + &
                     a7(3)*( wf(i,j,l+3) - wf(i,j,l-3) )
             dpphi = a7(1)*( pf(i,j,l+1) - pf(i,j,l-1) ) + &
                     a7(2)*( pf(i,j,l+2) - pf(i,j,l-2) ) + &
                     a7(3)*( pf(i,j,l+3) - pf(i,j,l-3) )
             drphi = a7(1)*( rf(i,j,l+1) - rf(i,j,l-1) ) + &
                     a7(2)*( rf(i,j,l+2) - rf(i,j,l-2) ) + &
                     a7(3)*( rf(i,j,l+3) - rf(i,j,l-3) )

             pt(i,j,l) = vg(i,j,l)*( &
                       ((dpksi*ksi_x(i,j,k)+dpeta*eta_x(i,j,k)+dpphi*phi_x(i,j,k))*sintcosp(i,j,l) &
                       +(dpksi*ksi_y(i,j,k)+dpeta*eta_y(i,j,k)+dpphi*phi_y(i,j,k))*sintsinp(i,j,l) &
                       +(dpksi*ksi_z(i,j,k)+dpeta*eta_z(i,j,k)+dpphi*phi_z(i,j,k))*costeta(i,j,l)  &
                       )*ijacob3(i,j,k) + pf(i,j,l)*ir(i,j,l) )
             ut(i,j,l) = vg(i,j,l)*( &
                       ((duksi*ksi_x(i,j,k)+dueta*eta_x(i,j,k)+duphi*phi_x(i,j,k))*sintcosp(i,j,l) &
                       +(duksi*ksi_y(i,j,k)+dueta*eta_y(i,j,k)+duphi*phi_y(i,j,k))*sintsinp(i,j,l) &
                       +(duksi*ksi_z(i,j,k)+dueta*eta_z(i,j,k)+duphi*phi_z(i,j,k))*costeta(i,j,l)  &
                       )*ijacob3(i,j,k) + uf(i,j,l)*ir(i,j,l) )
             vt(i,j,l) = vg(i,j,l)*( &
                       ((dvksi*ksi_x(i,j,k)+dveta*eta_x(i,j,k)+dvphi*phi_x(i,j,k))*sintcosp(i,j,l) &
                       +(dvksi*ksi_y(i,j,k)+dveta*eta_y(i,j,k)+dvphi*phi_y(i,j,k))*sintsinp(i,j,l) &
                       +(dvksi*ksi_z(i,j,k)+dveta*eta_z(i,j,k)+dvphi*phi_z(i,j,k))*costeta(i,j,l)  &
                       )*ijacob3(i,j,k) + vf(i,j,l)*ir(i,j,l) )
             wt(i,j,l) = vg(i,j,l)*( &
                       ((dwksi*ksi_x(i,j,k)+dweta*eta_x(i,j,k)+dwphi*phi_x(i,j,k))*sintcosp(i,j,l) &
                       +(dwksi*ksi_y(i,j,k)+dweta*eta_y(i,j,k)+dwphi*phi_y(i,j,k))*sintsinp(i,j,l) &
                       +(dwksi*ksi_z(i,j,k)+dweta*eta_z(i,j,k)+dwphi*phi_z(i,j,k))*costeta(i,j,l)  &
                       )*ijacob3(i,j,k) + wf(i,j,l)*ir(i,j,l) )
             rt(i,j,l) = vg(i,j,l)*( &
                       ((drksi*ksi_x(i,j,k)+dreta*eta_x(i,j,k)+drphi*phi_x(i,j,k))*sintcosp(i,j,l) &
                       +(drksi*ksi_y(i,j,k)+dreta*eta_y(i,j,k)+drphi*phi_y(i,j,k))*sintsinp(i,j,l) &
                       +(drksi*ksi_z(i,j,k)+dreta*eta_z(i,j,k)+drphi*phi_z(i,j,k))*costeta(i,j,l)  &
                       )*ijacob3(i,j,k) + rf(i,j,l)*ir(i,j,l) )
          enddo
       enddo
    enddo

    k=nz-2
    l=k-nzmngh
    do i=ndx,nfx
       do j=ndy,nfy
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

          dueta = a7(1)*( uf(i,j+1,l) - uf(i,j-1,l) ) + &
                  a7(2)*( uf(i,j+2,l) - uf(i,j-2,l) ) + &
                  a7(3)*( uf(i,j+3,l) - uf(i,j-3,l) )
          dveta = a7(1)*( vf(i,j+1,l) - vf(i,j-1,l) ) + &
                  a7(2)*( vf(i,j+2,l) - vf(i,j-2,l) ) + &
                  a7(3)*( vf(i,j+3,l) - vf(i,j-3,l) )
          dweta = a7(1)*( wf(i,j+1,l) - wf(i,j-1,l) ) + &
                  a7(2)*( wf(i,j+2,l) - wf(i,j-2,l) ) + &
                  a7(3)*( wf(i,j+3,l) - wf(i,j-3,l) )
          dpeta = a7(1)*( pf(i,j+1,l) - pf(i,j-1,l) ) + &
                  a7(2)*( pf(i,j+2,l) - pf(i,j-2,l) ) + &
                  a7(3)*( pf(i,j+3,l) - pf(i,j-3,l) )
          dreta = a7(1)*( rf(i,j+1,l) - rf(i,j-1,l) ) + &
                  a7(2)*( rf(i,j+2,l) - rf(i,j-2,l) ) + &
                  a7(3)*( rf(i,j+3,l) - rf(i,j-3,l) )

          duphi = a42(1)*uf(i,j,l+2)+a42(2)*uf(i,j,l+1) &
                + a42(3)*uf(i,j,l  )+a42(4)*uf(i,j,l-1) &
                + a42(5)*uf(i,j,l-2)+a42(6)*uf(i,j,l-3) &
                + a42(7)*uf(i,j,l-4)
          dvphi = a42(1)*vf(i,j,l+2)+a42(2)*vf(i,j,l+1) &
                + a42(3)*vf(i,j,l  )+a42(4)*vf(i,j,l-1) &
                + a42(5)*vf(i,j,l-2)+a42(6)*vf(i,j,l-3) &
                + a42(7)*vf(i,j,l-4)
          dwphi = a42(1)*wf(i,j,l+2)+a42(2)*wf(i,j,l+1) &
                + a42(3)*wf(i,j,l  )+a42(4)*wf(i,j,l-1) &
                + a42(5)*wf(i,j,l-2)+a42(6)*wf(i,j,l-3) &
                + a42(7)*wf(i,j,l-4)
          dpphi = a42(1)*pf(i,j,l+2)+a42(2)*pf(i,j,l+1) &
                + a42(3)*pf(i,j,l  )+a42(4)*pf(i,j,l-1) &
                + a42(5)*pf(i,j,l-2)+a42(6)*pf(i,j,l-3) &
                + a42(7)*pf(i,j,l-4)
          drphi = a42(1)*rf(i,j,l+2)+a42(2)*rf(i,j,l+1) &
                + a42(3)*rf(i,j,l  )+a42(4)*rf(i,j,l-1) &
                + a42(5)*rf(i,j,l-2)+a42(6)*rf(i,j,l-3) &
                + a42(7)*rf(i,j,l-4)

          pt(i,j,l) = vg(i,j,l)*( &
                    ((dpksi*ksi_x(i,j,k)+dpeta*eta_x(i,j,k)+dpphi*phi_x(i,j,k))*sintcosp(i,j,l) &
                    +(dpksi*ksi_y(i,j,k)+dpeta*eta_y(i,j,k)+dpphi*phi_y(i,j,k))*sintsinp(i,j,l) &
                    +(dpksi*ksi_z(i,j,k)+dpeta*eta_z(i,j,k)+dpphi*phi_z(i,j,k))*costeta(i,j,l)  &
                    )*ijacob3(i,j,k) + pf(i,j,l)*ir(i,j,l) )
          ut(i,j,l) = vg(i,j,l)*( &
                    ((duksi*ksi_x(i,j,k)+dueta*eta_x(i,j,k)+duphi*phi_x(i,j,k))*sintcosp(i,j,l) &
                    +(duksi*ksi_y(i,j,k)+dueta*eta_y(i,j,k)+duphi*phi_y(i,j,k))*sintsinp(i,j,l) &
                    +(duksi*ksi_z(i,j,k)+dueta*eta_z(i,j,k)+duphi*phi_z(i,j,k))*costeta(i,j,l)  &
                    )*ijacob3(i,j,k) + uf(i,j,l)*ir(i,j,l) )
          vt(i,j,l) = vg(i,j,l)*( &
                    ((dvksi*ksi_x(i,j,k)+dveta*eta_x(i,j,k)+dvphi*phi_x(i,j,k))*sintcosp(i,j,l) &
                    +(dvksi*ksi_y(i,j,k)+dveta*eta_y(i,j,k)+dvphi*phi_y(i,j,k))*sintsinp(i,j,l) &
                    +(dvksi*ksi_z(i,j,k)+dveta*eta_z(i,j,k)+dvphi*phi_z(i,j,k))*costeta(i,j,l)  &
                    )*ijacob3(i,j,k) + vf(i,j,l)*ir(i,j,l) )
          wt(i,j,l) = vg(i,j,l)*( &
                    ((dwksi*ksi_x(i,j,k)+dweta*eta_x(i,j,k)+dwphi*phi_x(i,j,k))*sintcosp(i,j,l) &
                    +(dwksi*ksi_y(i,j,k)+dweta*eta_y(i,j,k)+dwphi*phi_y(i,j,k))*sintsinp(i,j,l) &
                    +(dwksi*ksi_z(i,j,k)+dweta*eta_z(i,j,k)+dwphi*phi_z(i,j,k))*costeta(i,j,l)  &
                    )*ijacob3(i,j,k) + wf(i,j,l)*ir(i,j,l) )
          rt(i,j,l) = vg(i,j,l)*( &
                    ((drksi*ksi_x(i,j,k)+dreta*eta_x(i,j,k)+drphi*phi_x(i,j,k))*sintcosp(i,j,l) &
                    +(drksi*ksi_y(i,j,k)+dreta*eta_y(i,j,k)+drphi*phi_y(i,j,k))*sintsinp(i,j,l) &
                    +(drksi*ksi_z(i,j,k)+dreta*eta_z(i,j,k)+drphi*phi_z(i,j,k))*costeta(i,j,l)  &
                    )*ijacob3(i,j,k) + rf(i,j,l)*ir(i,j,l) )
       enddo
    enddo

    k=nz-1
    l=k-nzmngh
    do i=ndx,nfx
       do j=ndy,nfy
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

          dueta = a7(1)*( uf(i,j+1,l) - uf(i,j-1,l) ) + &
                  a7(2)*( uf(i,j+2,l) - uf(i,j-2,l) ) + &
                  a7(3)*( uf(i,j+3,l) - uf(i,j-3,l) )
          dveta = a7(1)*( vf(i,j+1,l) - vf(i,j-1,l) ) + &
                  a7(2)*( vf(i,j+2,l) - vf(i,j-2,l) ) + &
                  a7(3)*( vf(i,j+3,l) - vf(i,j-3,l) )
          dweta = a7(1)*( wf(i,j+1,l) - wf(i,j-1,l) ) + &
                  a7(2)*( wf(i,j+2,l) - wf(i,j-2,l) ) + &
                  a7(3)*( wf(i,j+3,l) - wf(i,j-3,l) )
          dpeta = a7(1)*( pf(i,j+1,l) - pf(i,j-1,l) ) + &
                  a7(2)*( pf(i,j+2,l) - pf(i,j-2,l) ) + &
                  a7(3)*( pf(i,j+3,l) - pf(i,j-3,l) )
          dreta = a7(1)*( rf(i,j+1,l) - rf(i,j-1,l) ) + &
                  a7(2)*( rf(i,j+2,l) - rf(i,j-2,l) ) + &
                  a7(3)*( rf(i,j+3,l) - rf(i,j-3,l) )

          duphi = a51(1)*uf(i,j,l+1)+a51(2)*uf(i,j,l  ) &
                + a51(3)*uf(i,j,l-1)+a51(4)*uf(i,j,l-2) &
                + a51(5)*uf(i,j,l-3)+a51(6)*uf(i,j,l-4) &
                + a51(7)*uf(i,j,l-5)
          dvphi = a51(1)*vf(i,j,l+1)+a51(2)*vf(i,j,l  ) &
                + a51(3)*vf(i,j,l-1)+a51(4)*vf(i,j,l-2) &
                + a51(5)*vf(i,j,l-3)+a51(6)*vf(i,j,l-4) &
                + a51(7)*vf(i,j,l-5)
          dwphi = a51(1)*wf(i,j,l+1)+a51(2)*wf(i,j,l  ) &
                + a51(3)*wf(i,j,l-1)+a51(4)*wf(i,j,l-2) &
                + a51(5)*wf(i,j,l-3)+a51(6)*wf(i,j,l-4) &
                + a51(7)*wf(i,j,l-5)
          dpphi = a51(1)*pf(i,j,l+1)+a51(2)*pf(i,j,l  ) &
                + a51(3)*pf(i,j,l-1)+a51(4)*pf(i,j,l-2) &
                + a51(5)*pf(i,j,l-3)+a51(6)*pf(i,j,l-4) &
                + a51(7)*pf(i,j,l-5)
          drphi = a51(1)*rf(i,j,l+1)+a51(2)*rf(i,j,l  ) &
                + a51(3)*rf(i,j,l-1)+a51(4)*rf(i,j,l-2) &
                + a51(5)*rf(i,j,l-3)+a51(6)*rf(i,j,l-4) &
                + a51(7)*rf(i,j,l-5)

          pt(i,j,l) = vg(i,j,l)*( &
                    ((dpksi*ksi_x(i,j,k)+dpeta*eta_x(i,j,k)+dpphi*phi_x(i,j,k))*sintcosp(i,j,l) &
                    +(dpksi*ksi_y(i,j,k)+dpeta*eta_y(i,j,k)+dpphi*phi_y(i,j,k))*sintsinp(i,j,l) &
                    +(dpksi*ksi_z(i,j,k)+dpeta*eta_z(i,j,k)+dpphi*phi_z(i,j,k))*costeta(i,j,l)  &
                    )*ijacob3(i,j,k) + pf(i,j,l)*ir(i,j,l) )
          ut(i,j,l) = vg(i,j,l)*( &
                    ((duksi*ksi_x(i,j,k)+dueta*eta_x(i,j,k)+duphi*phi_x(i,j,k))*sintcosp(i,j,l) &
                    +(duksi*ksi_y(i,j,k)+dueta*eta_y(i,j,k)+duphi*phi_y(i,j,k))*sintsinp(i,j,l) &
                    +(duksi*ksi_z(i,j,k)+dueta*eta_z(i,j,k)+duphi*phi_z(i,j,k))*costeta(i,j,l)  &
                    )*ijacob3(i,j,k) + uf(i,j,l)*ir(i,j,l) )
          vt(i,j,l) = vg(i,j,l)*( &
                    ((dvksi*ksi_x(i,j,k)+dveta*eta_x(i,j,k)+dvphi*phi_x(i,j,k))*sintcosp(i,j,l) &
                    +(dvksi*ksi_y(i,j,k)+dveta*eta_y(i,j,k)+dvphi*phi_y(i,j,k))*sintsinp(i,j,l) &
                    +(dvksi*ksi_z(i,j,k)+dveta*eta_z(i,j,k)+dvphi*phi_z(i,j,k))*costeta(i,j,l)  &
                    )*ijacob3(i,j,k) + vf(i,j,l)*ir(i,j,l) )
          wt(i,j,l) = vg(i,j,l)*( &
                    ((dwksi*ksi_x(i,j,k)+dweta*eta_x(i,j,k)+dwphi*phi_x(i,j,k))*sintcosp(i,j,l) &
                    +(dwksi*ksi_y(i,j,k)+dweta*eta_y(i,j,k)+dwphi*phi_y(i,j,k))*sintsinp(i,j,l) &
                    +(dwksi*ksi_z(i,j,k)+dweta*eta_z(i,j,k)+dwphi*phi_z(i,j,k))*costeta(i,j,l)  &
                    )*ijacob3(i,j,k) + wf(i,j,l)*ir(i,j,l) )
          rt(i,j,l) = vg(i,j,l)*( &
                    ((drksi*ksi_x(i,j,k)+dreta*eta_x(i,j,k)+drphi*phi_x(i,j,k))*sintcosp(i,j,l) &
                    +(drksi*ksi_y(i,j,k)+dreta*eta_y(i,j,k)+drphi*phi_y(i,j,k))*sintsinp(i,j,l) &
                    +(drksi*ksi_z(i,j,k)+dreta*eta_z(i,j,k)+drphi*phi_z(i,j,k))*costeta(i,j,l)  &
                    )*ijacob3(i,j,k) + rf(i,j,l)*ir(i,j,l) )
       enddo
    enddo

    k=nz
    l=k-nzmngh
    do i=ndx,nfx
       do j=ndy,nfy
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

          dueta = a7(1)*( uf(i,j+1,l) - uf(i,j-1,l) ) + &
                  a7(2)*( uf(i,j+2,l) - uf(i,j-2,l) ) + &
                  a7(3)*( uf(i,j+3,l) - uf(i,j-3,l) )
          dveta = a7(1)*( vf(i,j+1,l) - vf(i,j-1,l) ) + &
                  a7(2)*( vf(i,j+2,l) - vf(i,j-2,l) ) + &
                  a7(3)*( vf(i,j+3,l) - vf(i,j-3,l) )
          dweta = a7(1)*( wf(i,j+1,l) - wf(i,j-1,l) ) + &
                  a7(2)*( wf(i,j+2,l) - wf(i,j-2,l) ) + &
                  a7(3)*( wf(i,j+3,l) - wf(i,j-3,l) )
          dpeta = a7(1)*( pf(i,j+1,l) - pf(i,j-1,l) ) + &
                  a7(2)*( pf(i,j+2,l) - pf(i,j-2,l) ) + &
                  a7(3)*( pf(i,j+3,l) - pf(i,j-3,l) )
          dreta = a7(1)*( rf(i,j+1,l) - rf(i,j-1,l) ) + &
                  a7(2)*( rf(i,j+2,l) - rf(i,j-2,l) ) + &
                  a7(3)*( rf(i,j+3,l) - rf(i,j-3,l) )

          duphi = a60(1)*uf(i,j,l  )+a60(2)*uf(i,j,l-1) &
                + a60(3)*uf(i,j,l-2)+a60(4)*uf(i,j,l-3) &
                + a60(5)*uf(i,j,l-4)+a60(6)*uf(i,j,l-5) &
                + a60(7)*uf(i,j,l-6)
          dvphi = a60(1)*vf(i,j,l  )+a60(2)*vf(i,j,l-1) &
                + a60(3)*vf(i,j,l-2)+a60(4)*vf(i,j,l-3) &
                + a60(5)*vf(i,j,l-4)+a60(6)*vf(i,j,l-5) &
                + a60(7)*vf(i,j,l-6)
          dwphi = a60(1)*wf(i,j,l  )+a60(2)*wf(i,j,l-1) &
                + a60(3)*wf(i,j,l-2)+a60(4)*wf(i,j,l-3) &
                + a60(5)*wf(i,j,l-4)+a60(6)*wf(i,j,l-5) &
                + a60(7)*wf(i,j,l-6)
          dpphi = a60(1)*pf(i,j,l  )+a60(2)*pf(i,j,l-1) &
                + a60(3)*pf(i,j,l-2)+a60(4)*pf(i,j,l-3) &
                + a60(5)*pf(i,j,l-4)+a60(6)*pf(i,j,l-5) &
                + a60(7)*pf(i,j,l-6)
          drphi = a60(1)*rf(i,j,l  )+a60(2)*rf(i,j,l-1) &
                + a60(3)*rf(i,j,l-2)+a60(4)*rf(i,j,l-3) &
                + a60(5)*rf(i,j,l-4)+a60(6)*rf(i,j,l-5) &
                + a60(7)*rf(i,j,l-6)

          pt(i,j,l) = vg(i,j,l)*( &
                    ((dpksi*ksi_x(i,j,k)+dpeta*eta_x(i,j,k)+dpphi*phi_x(i,j,k))*sintcosp(i,j,l) &
                    +(dpksi*ksi_y(i,j,k)+dpeta*eta_y(i,j,k)+dpphi*phi_y(i,j,k))*sintsinp(i,j,l) &
                    +(dpksi*ksi_z(i,j,k)+dpeta*eta_z(i,j,k)+dpphi*phi_z(i,j,k))*costeta(i,j,l)  &
                    )*ijacob3(i,j,k) + pf(i,j,l)*ir(i,j,l) )
          ut(i,j,l) = vg(i,j,l)*( &
                    ((duksi*ksi_x(i,j,k)+dueta*eta_x(i,j,k)+duphi*phi_x(i,j,k))*sintcosp(i,j,l) &
                    +(duksi*ksi_y(i,j,k)+dueta*eta_y(i,j,k)+duphi*phi_y(i,j,k))*sintsinp(i,j,l) &
                    +(duksi*ksi_z(i,j,k)+dueta*eta_z(i,j,k)+duphi*phi_z(i,j,k))*costeta(i,j,l)  &
                    )*ijacob3(i,j,k) + uf(i,j,l)*ir(i,j,l) )
          vt(i,j,l) = vg(i,j,l)*( &
                    ((dvksi*ksi_x(i,j,k)+dveta*eta_x(i,j,k)+dvphi*phi_x(i,j,k))*sintcosp(i,j,l) &
                    +(dvksi*ksi_y(i,j,k)+dveta*eta_y(i,j,k)+dvphi*phi_y(i,j,k))*sintsinp(i,j,l) &
                    +(dvksi*ksi_z(i,j,k)+dveta*eta_z(i,j,k)+dvphi*phi_z(i,j,k))*costeta(i,j,l)  &
                    )*ijacob3(i,j,k) + vf(i,j,l)*ir(i,j,l) )
          wt(i,j,l) = vg(i,j,l)*( &
                    ((dwksi*ksi_x(i,j,k)+dweta*eta_x(i,j,k)+dwphi*phi_x(i,j,k))*sintcosp(i,j,l) &
                    +(dwksi*ksi_y(i,j,k)+dweta*eta_y(i,j,k)+dwphi*phi_y(i,j,k))*sintsinp(i,j,l) &
                    +(dwksi*ksi_z(i,j,k)+dweta*eta_z(i,j,k)+dwphi*phi_z(i,j,k))*costeta(i,j,l)  &
                    )*ijacob3(i,j,k) + wf(i,j,l)*ir(i,j,l) )
          rt(i,j,l) = vg(i,j,l)*( &
                    ((drksi*ksi_x(i,j,k)+dreta*eta_x(i,j,k)+drphi*phi_x(i,j,k))*sintcosp(i,j,l) &
                    +(drksi*ksi_y(i,j,k)+dreta*eta_y(i,j,k)+drphi*phi_y(i,j,k))*sintsinp(i,j,l) &
                    +(drksi*ksi_z(i,j,k)+dreta*eta_z(i,j,k)+drphi*phi_z(i,j,k))*costeta(i,j,l)  &
                    )*ijacob3(i,j,k) + rf(i,j,l)*ir(i,j,l) )
       enddo
    enddo

    ! Update fluxes at each RK step
    ! =============================
    do k=nzmnghp1,nz
       l=k-nzmngh
       do j=ndy,nfy
          do i=ndx,nfx
             cp=cpcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
             av=avcalc_tro(Tmp(i,j,k),rho_n(i,j,k))

             Krho(i,j,k)  = rt(i,j,l)
             Krhou(i,j,k) = uu(i,j,k)*rt(i,j,l)+rho_n(i,j,k)*ut(i,j,l)
             Krhov(i,j,k) = vv(i,j,k)*rt(i,j,l)+rho_n(i,j,k)*vt(i,j,l)
             Krhow(i,j,k) = ww(i,j,k)*rt(i,j,l)+rho_n(i,j,k)*wt(i,j,l)
             Krhoe(i,j,k) = cp/av*(pt(i,j,l)/c2_(i,j,l)-rt(i,j,l)) &
                  + (rhoe_n(i,j,k)+prs(i,j,k))/rho_n(i,j,k)*rt(i,j,l) &
                  + rho_n(i,j,k)*(uu(i,j,k)*ut(i,j,l)+vv(i,j,k)*vt(i,j,l)+ww(i,j,k)*wt(i,j,l))
          enddo
       enddo
    enddo

  end subroutine bc_TD3d_kmax_c3

end submodule smod_TamDong3d_faces_c3
