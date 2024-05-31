!===============================================================================
submodule (mod_TamDong3d_c3) smod_TamDong3d_corner_c3
!===============================================================================
  !> author: XG
  !> date: April 2023
  !> Non-reflecting boundary conditions of Tam and Dong
  !> - 3D curvilinear version - routines for coorners
!=============================================================================== 

contains

  !===============================================================================
  module subroutine bc_TD3d_imin_jmin_kmin_c3
  !===============================================================================
    !> 3D Tam & Dong's BC at imin-jmin-kmin (corner 1,1,1 /left-bottom-front)
    !> - 3D curvilinear version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: cp,av,c2_
    real(wp), dimension(1:ngh+3,1:ngh+3,1:ngh+3) :: rf,uf,vf,wf,pf
    real(wp), dimension(ngh,ngh,ngh) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(ngh,ngh,ngh) :: duksi,dvksi,dwksi,dpksi,drksi
    real(wp), dimension(ngh,ngh,ngh) :: dueta,dveta,dweta,dpeta,dreta
    real(wp), dimension(ngh,ngh,ngh) :: duphi,dvphi,dwphi,dpphi,drphi
    !-------------------------------------------------------------------------

    ! Pointers for spherical coordinates
    ! ==================================
    ir=>BC_corner(1,1,1)%i_r
    cosphi=>BC_corner(1,1,1)%cosp
    sinphi=>BC_corner(1,1,1)%sinp
    costeta=>BC_corner(1,1,1)%cost
    sinteta=>BC_corner(1,1,1)%sint
    costcosp=>BC_corner(1,1,1)%costcosp
    costsinp=>BC_corner(1,1,1)%costsinp
    sintcosp=>BC_corner(1,1,1)%sintcosp
    sintsinp=>BC_corner(1,1,1)%sintsinp
    
    ! Compute fluctuations
    ! ====================
    do i=1,nghp3
       do j=1,nghp3
          do k=1,nghp3
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
       do j=1,ngh
          do k=1,ngh
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

    ! Non-centered derivatives in x-direction *sin(teta)*cos(phi)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    i=1
    do j=1,ngh
       do k=1,ngh
          duksi(i,j,k)= a06(1)*uf(1,j,k)+a06(2)*uf(2,j,k) &
                      + a06(3)*uf(3,j,k)+a06(4)*uf(4,j,k) &
                      + a06(5)*uf(5,j,k)+a06(6)*uf(6,j,k) &
                      + a06(7)*uf(7,j,k)
          dvksi(i,j,k)= a06(1)*vf(1,j,k)+a06(2)*vf(2,j,k) &
                      + a06(3)*vf(3,j,k)+a06(4)*vf(4,j,k) &
                      + a06(5)*vf(5,j,k)+a06(6)*vf(6,j,k) &
                      + a06(7)*vf(7,j,k)
          dwksi(i,j,k)= a06(1)*wf(1,j,k)+a06(2)*wf(2,j,k) &
                      + a06(3)*wf(3,j,k)+a06(4)*wf(4,j,k) &
                      + a06(5)*wf(5,j,k)+a06(6)*wf(6,j,k) &
                      + a06(7)*wf(7,j,k)
          dpksi(i,j,k)= a06(1)*pf(1,j,k)+a06(2)*pf(2,j,k) &
                      + a06(3)*pf(3,j,k)+a06(4)*pf(4,j,k) &
                      + a06(5)*pf(5,j,k)+a06(6)*pf(6,j,k) &
                      + a06(7)*pf(7,j,k)
          drksi(i,j,k)= a06(1)*rf(1,j,k)+a06(2)*rf(2,j,k) &
                      + a06(3)*rf(3,j,k)+a06(4)*rf(4,j,k) &
                      + a06(5)*rf(5,j,k)+a06(6)*rf(6,j,k) &
                      + a06(7)*rf(7,j,k)
       enddo
    enddo

    i=2
    do j=1,ngh
       do k=1,ngh
          duksi(i,j,k)= a15(1)*uf(1,j,k)+a15(2)*uf(2,j,k) &
                      + a15(3)*uf(3,j,k)+a15(4)*uf(4,j,k) &
                      + a15(5)*uf(5,j,k)+a15(6)*uf(6,j,k) &
                      + a15(7)*uf(7,j,k)
          dvksi(i,j,k)= a15(1)*vf(1,j,k)+a15(2)*vf(2,j,k) &
                      + a15(3)*vf(3,j,k)+a15(4)*vf(4,j,k) &
                      + a15(5)*vf(5,j,k)+a15(6)*vf(6,j,k) &
                      + a15(7)*vf(7,j,k)
          dwksi(i,j,k)= a15(1)*wf(1,j,k)+a15(2)*wf(2,j,k) &
                      + a15(3)*wf(3,j,k)+a15(4)*wf(4,j,k) &
                      + a15(5)*wf(5,j,k)+a15(6)*wf(6,j,k) &
                      + a15(7)*wf(7,j,k)
          dpksi(i,j,k)= a15(1)*pf(1,j,k)+a15(2)*pf(2,j,k) &
                      + a15(3)*pf(3,j,k)+a15(4)*pf(4,j,k) &
                      + a15(5)*pf(5,j,k)+a15(6)*pf(6,j,k) &
                      + a15(7)*pf(7,j,k)
          drksi(i,j,k)= a15(1)*rf(1,j,k)+a15(2)*rf(2,j,k) &
                      + a15(3)*rf(3,j,k)+a15(4)*rf(4,j,k) &
                      + a15(5)*rf(5,j,k)+a15(6)*rf(6,j,k) &
                      + a15(7)*rf(7,j,k)
       enddo
    enddo

    i=3
    do j=1,ngh
       do k=1,ngh
          duksi(i,j,k)= a24(1)*uf(1,j,k)+a24(2)*uf(2,j,k) &
                      + a24(3)*uf(3,j,k)+a24(4)*uf(4,j,k) &
                      + a24(5)*uf(5,j,k)+a24(6)*uf(6,j,k) &
                      + a24(7)*uf(7,j,k)
          dvksi(i,j,k)= a24(1)*vf(1,j,k)+a24(2)*vf(2,j,k) &
                      + a24(3)*vf(3,j,k)+a24(4)*vf(4,j,k) &
                      + a24(5)*vf(5,j,k)+a24(6)*vf(6,j,k) &
                      + a24(7)*vf(7,j,k)
          dwksi(i,j,k)= a24(1)*wf(1,j,k)+a24(2)*wf(2,j,k) &
                      + a24(3)*wf(3,j,k)+a24(4)*wf(4,j,k) &
                      + a24(5)*wf(5,j,k)+a24(6)*wf(6,j,k) &
                      + a24(7)*wf(7,j,k)
          dpksi(i,j,k)= a24(1)*pf(1,j,k)+a24(2)*pf(2,j,k) &
                      + a24(3)*pf(3,j,k)+a24(4)*pf(4,j,k) &
                      + a24(5)*pf(5,j,k)+a24(6)*pf(6,j,k) &
                      + a24(7)*pf(7,j,k)
          drksi(i,j,k)= a24(1)*rf(1,j,k)+a24(2)*rf(2,j,k) &
                      + a24(3)*rf(3,j,k)+a24(4)*rf(4,j,k) &
                      + a24(5)*rf(5,j,k)+a24(6)*rf(6,j,k) &
                      + a24(7)*rf(7,j,k)
       enddo
    enddo

    do i=4,ngh
       do j=1,ngh
          do k=1,ngh
             duksi(i,j,k)= a7(1)* (uf(i+1,j,k)-uf(i-1,j,k)) &
                         + a7(2)* (uf(i+2,j,k)-uf(i-2,j,k)) &
                         + a7(3)* (uf(i+3,j,k)-uf(i-3,j,k))
             dvksi(i,j,k)= a7(1)* (vf(i+1,j,k)-vf(i-1,j,k)) &
                         + a7(2)* (vf(i+2,j,k)-vf(i-2,j,k)) &
                         + a7(3)* (vf(i+3,j,k)-vf(i-3,j,k))
             dwksi(i,j,k)= a7(1)* (wf(i+1,j,k)-wf(i-1,j,k)) &
                         + a7(2)* (wf(i+2,j,k)-wf(i-2,j,k)) &
                         + a7(3)* (wf(i+3,j,k)-wf(i-3,j,k))
             dpksi(i,j,k)= a7(1)* (pf(i+1,j,k)-pf(i-1,j,k)) &
                         + a7(2)* (pf(i+2,j,k)-pf(i-2,j,k)) &
                         + a7(3)* (pf(i+3,j,k)-pf(i-3,j,k))
             drksi(i,j,k)= a7(1)* (rf(i+1,j,k)-rf(i-1,j,k)) &
                         + a7(2)* (rf(i+2,j,k)-rf(i-2,j,k)) &
                         + a7(3)* (rf(i+3,j,k)-rf(i-3,j,k))
          enddo
       enddo
    enddo

    ! Non-centered derivatives in y-direction *sin(teta)*sin(phi)
    ! =======================================
    j=1
    do i=1,ngh
       do k=1,ngh
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
    do i=1,ngh
       do k=1,ngh
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
    do i=1,ngh
       do k=1,ngh
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
       do i=1,ngh
          do k=1,ngh
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
    do i=1,ngh
       do j=1,ngh
          duphi(i,j,k)=  a06(1)*uf(i,j,1)+a06(2)*uf(i,j,2) &
                       + a06(3)*uf(i,j,3)+a06(4)*uf(i,j,4) &
                       + a06(5)*uf(i,j,5)+a06(6)*uf(i,j,6) &
                       + a06(7)*uf(i,j,7)
          dvphi(i,j,k)=  a06(1)*vf(i,j,1)+a06(2)*vf(i,j,2) &
                       + a06(3)*vf(i,j,3)+a06(4)*vf(i,j,4) &
                       + a06(5)*vf(i,j,5)+a06(6)*vf(i,j,6) &
                       + a06(7)*vf(i,j,7)
          dwphi(i,j,k)=  a06(1)*wf(i,j,1)+a06(2)*wf(i,j,2) &
                       + a06(3)*wf(i,j,3)+a06(4)*wf(i,j,4) &
                       + a06(5)*wf(i,j,5)+a06(6)*wf(i,j,6) &
                       + a06(7)*wf(i,j,7)
          dpphi(i,j,k)=  a06(1)*pf(i,j,1)+a06(2)*pf(i,j,2) &
                       + a06(3)*pf(i,j,3)+a06(4)*pf(i,j,4) &
                       + a06(5)*pf(i,j,5)+a06(6)*pf(i,j,6) &
                       + a06(7)*pf(i,j,7)
          drphi(i,j,k)=  a06(1)*rf(i,j,1)+a06(2)*rf(i,j,2) &
                       + a06(3)*rf(i,j,3)+a06(4)*rf(i,j,4) &
                       + a06(5)*rf(i,j,5)+a06(6)*rf(i,j,6) &
                       + a06(7)*rf(i,j,7)
       enddo
    enddo

    k=2
    do i=1,ngh
       do j=1,ngh
          duphi(i,j,k)=  a15(1)*uf(i,j,1)+a15(2)*uf(i,j,2) &
                       + a15(3)*uf(i,j,3)+a15(4)*uf(i,j,4) &
                       + a15(5)*uf(i,j,5)+a15(6)*uf(i,j,6) &
                       + a15(7)*uf(i,j,7)
          dvphi(i,j,k)=  a15(1)*vf(i,j,1)+a15(2)*vf(i,j,2) &
                       + a15(3)*vf(i,j,3)+a15(4)*vf(i,j,4) &
                       + a15(5)*vf(i,j,5)+a15(6)*vf(i,j,6) &
                       + a15(7)*vf(i,j,7)
          dwphi(i,j,k)=  a15(1)*wf(i,j,1)+a15(2)*wf(i,j,2) &
                       + a15(3)*wf(i,j,3)+a15(4)*wf(i,j,4) &
                       + a15(5)*wf(i,j,5)+a15(6)*wf(i,j,6) &
                       + a15(7)*wf(i,j,7)
          dpphi(i,j,k)=  a15(1)*pf(i,j,1)+a15(2)*pf(i,j,2) &
                       + a15(3)*pf(i,j,3)+a15(4)*pf(i,j,4) &
                       + a15(5)*pf(i,j,5)+a15(6)*pf(i,j,6) &
                       + a15(7)*pf(i,j,7)
          drphi(i,j,k)=  a15(1)*rf(i,j,1)+a15(2)*rf(i,j,2) &
                       + a15(3)*rf(i,j,3)+a15(4)*rf(i,j,4) &
                       + a15(5)*rf(i,j,5)+a15(6)*rf(i,j,6) &
                       + a15(7)*rf(i,j,7)
       enddo
    enddo

    k=3
    do i=1,ngh
       do j=1,ngh
          duphi(i,j,k)=  a24(1)*uf(i,j,1)+a24(2)*uf(i,j,2) &
                       + a24(3)*uf(i,j,3)+a24(4)*uf(i,j,4) &
                       + a24(5)*uf(i,j,5)+a24(6)*uf(i,j,6) &
                       + a24(7)*uf(i,j,7)
          dvphi(i,j,k)=  a24(1)*vf(i,j,1)+a24(2)*vf(i,j,2) &
                       + a24(3)*vf(i,j,3)+a24(4)*vf(i,j,4) &
                       + a24(5)*vf(i,j,5)+a24(6)*vf(i,j,6) &
                       + a24(7)*vf(i,j,7)
          dwphi(i,j,k)=  a24(1)*wf(i,j,1)+a24(2)*wf(i,j,2) &
                       + a24(3)*wf(i,j,3)+a24(4)*wf(i,j,4) &
                       + a24(5)*wf(i,j,5)+a24(6)*wf(i,j,6) &
                       + a24(7)*wf(i,j,7)
          dpphi(i,j,k)=  a24(1)*pf(i,j,1)+a24(2)*pf(i,j,2) &
                       + a24(3)*pf(i,j,3)+a24(4)*pf(i,j,4) &
                       + a24(5)*pf(i,j,5)+a24(6)*pf(i,j,6) &
                       + a24(7)*pf(i,j,7)
          drphi(i,j,k)=  a24(1)*rf(i,j,1)+a24(2)*rf(i,j,2) &
                       + a24(3)*rf(i,j,3)+a24(4)*rf(i,j,4) &
                       + a24(5)*rf(i,j,5)+a24(6)*rf(i,j,6) &
                       + a24(7)*rf(i,j,7)
       enddo
    enddo

    do k=4,ngh
       do i=1,ngh
          do j=1,ngh
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
    do i=1,ngh
       do j=1,ngh
          do k=1,ngh
             pt(i,j,k) = vg(i,j,k)*( &
                       ((dpksi(i,j,k)*ksi_x(i,j,k)+dpeta(i,j,k)*eta_x(i,j,k)+dpphi(i,j,k)*phi_x(i,j,k))*sintcosp(i,j,k) &
                       +(dpksi(i,j,k)*ksi_y(i,j,k)+dpeta(i,j,k)*eta_y(i,j,k)+dpphi(i,j,k)*phi_y(i,j,k))*sintsinp(i,j,k) &
                       +(dpksi(i,j,k)*ksi_z(i,j,k)+dpeta(i,j,k)*eta_z(i,j,k)+dpphi(i,j,k)*phi_z(i,j,k))*costeta(i,j,k)  &
                       )*ijacob3(i,j,k) + pf(i,j,k)*ir(i,j,k) )
             ut(i,j,k) = vg(i,j,k)*( &
                       ((duksi(i,j,k)*ksi_x(i,j,k)+dueta(i,j,k)*eta_x(i,j,k)+duphi(i,j,k)*phi_x(i,j,k))*sintcosp(i,j,k) &
                       +(duksi(i,j,k)*ksi_y(i,j,k)+dueta(i,j,k)*eta_y(i,j,k)+duphi(i,j,k)*phi_y(i,j,k))*sintsinp(i,j,k) &
                       +(duksi(i,j,k)*ksi_z(i,j,k)+dueta(i,j,k)*eta_z(i,j,k)+duphi(i,j,k)*phi_z(i,j,k))*costeta(i,j,k)  &
                       )*ijacob3(i,j,k) + uf(i,j,k)*ir(i,j,k) )
             vt(i,j,k) = vg(i,j,k)*( &
                       ((dvksi(i,j,k)*ksi_x(i,j,k)+dveta(i,j,k)*eta_x(i,j,k)+dvphi(i,j,k)*phi_x(i,j,k))*sintcosp(i,j,k) &
                       +(dvksi(i,j,k)*ksi_y(i,j,k)+dveta(i,j,k)*eta_y(i,j,k)+dvphi(i,j,k)*phi_y(i,j,k))*sintsinp(i,j,k) &
                       +(dvksi(i,j,k)*ksi_z(i,j,k)+dveta(i,j,k)*eta_z(i,j,k)+dvphi(i,j,k)*phi_z(i,j,k))*costeta(i,j,k)  &
                       )*ijacob3(i,j,k) + vf(i,j,k)*ir(i,j,k) )
             wt(i,j,k) = vg(i,j,k)*( &
                       ((dwksi(i,j,k)*ksi_x(i,j,k)+dweta(i,j,k)*eta_x(i,j,k)+dwphi(i,j,k)*phi_x(i,j,k))*sintcosp(i,j,k) &
                       +(dwksi(i,j,k)*ksi_y(i,j,k)+dweta(i,j,k)*eta_y(i,j,k)+dwphi(i,j,k)*phi_y(i,j,k))*sintsinp(i,j,k) &
                       +(dwksi(i,j,k)*ksi_z(i,j,k)+dweta(i,j,k)*eta_z(i,j,k)+dwphi(i,j,k)*phi_z(i,j,k))*costeta(i,j,k)  &
                       )*ijacob3(i,j,k) + wf(i,j,k)*ir(i,j,k) )
             rt(i,j,k) = vg(i,j,k)*( &
                       ((drksi(i,j,k)*ksi_x(i,j,k)+dreta(i,j,k)*eta_x(i,j,k)+drphi(i,j,k)*phi_x(i,j,k))*sintcosp(i,j,k) &
                       +(drksi(i,j,k)*ksi_y(i,j,k)+dreta(i,j,k)*eta_y(i,j,k)+drphi(i,j,k)*phi_y(i,j,k))*sintsinp(i,j,k) &
                       +(drksi(i,j,k)*ksi_z(i,j,k)+dreta(i,j,k)*eta_z(i,j,k)+drphi(i,j,k)*phi_z(i,j,k))*costeta(i,j,k)  &
                       )*ijacob3(i,j,k) + rf(i,j,k)*ir(i,j,k) )
          enddo
       enddo
    enddo

    ! Update fluxes at each RK step
    ! =============================
    do k=1,ngh
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

  end subroutine bc_TD3d_imin_jmin_kmin_c3

  !===============================================================================
  module subroutine bc_TD3d_imin_jmin_kmax_c3
  !===============================================================================
    !> 3D Tam & Dong's BC at imin-jmin-kmax (corner 1,1,2 /left-bottom-back)
    !> - 3D curvilinear version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l
    real(wp) :: cp,av,c2_
    real(wp), dimension(1:ngh+3,1:ngh+3,-2:ngh) :: rf,uf,vf,wf,pf
    real(wp), dimension(ngh,ngh,ngh) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(ngh,ngh,ngh) :: duksi,dvksi,dwksi,dpksi,drksi
    real(wp), dimension(ngh,ngh,ngh) :: dueta,dveta,dweta,dpeta,dreta
    real(wp), dimension(ngh,ngh,ngh) :: duphi,dvphi,dwphi,dpphi,drphi
    !-------------------------------------------------------------------------

    ! Pointers for spherical coordinates
    ! ==================================
    ir=>BC_corner(1,1,2)%i_r
    cosphi=>BC_corner(1,1,2)%cosp
    sinphi=>BC_corner(1,1,2)%sinp
    costeta=>BC_corner(1,1,2)%cost
    sinteta=>BC_corner(1,1,2)%sint
    costcosp=>BC_corner(1,1,2)%costcosp
    costsinp=>BC_corner(1,1,2)%costsinp
    sintcosp=>BC_corner(1,1,2)%sintcosp
    sintsinp=>BC_corner(1,1,2)%sintsinp

    ! Compute fluctuations
    ! ====================
    do i=1,nghp3
       do j=1,nghp3
          do k=nzmngh-2,nz
             l=k-nzmngh
             rf(i,j,l)=rho_n(i,j,k)-BC_face(1,1)%U0(i,j,k,1)
             uf(i,j,l)=   uu(i,j,k)-BC_face(1,1)%U0(i,j,k,2)
             vf(i,j,l)=   vv(i,j,k)-BC_face(1,1)%U0(i,j,k,3)
             wf(i,j,l)=   ww(i,j,k)-BC_face(1,1)%U0(i,j,k,4)
             pf(i,j,l)=  prs(i,j,k)-BC_face(1,1)%U0(i,j,k,5)
          enddo
       enddo
    enddo

    ! Compute group velocity vg
    ! =========================
    ! vg= u0*sin(teta)*cos(phi)+v0*sin(teta)*sin(phi)+w0*cos(teta)
    !   + sqrt{c0^2-[u0*cos(teta)*cos(phi)+v0*cos(teta)*sin(phi)-w0*sin(teta)]^2-[u0*sin(phi)-v0*cos(phi)]^2}
    do i=1,ngh
       do j=1,ngh
          do k=nzmnghp1,nz
             l=k-nzmngh
             vg(i,j,l)= BC_face(1,1)%U0(i,j,k,2)*sintcosp(i,j,l) +     &
               BC_face(1,1)%U0(i,j,k,3)*sintsinp(i,j,l) +     &
               BC_face(1,1)%U0(i,j,k,4)*costeta(i,j,l)          &
                      + sqrt( BC_face(1,1)%U0(i,j,l,6)-                &
                       ( BC_face(1,1)%U0(i,j,k,2)*costcosp(i,j,l) +    &
                   BC_face(1,1)%U0(i,j,k,3)*costsinp(i,j,l) -    &
                     BC_face(1,1)%U0(i,j,k,4)*sinteta(i,j,l) )**2- &
                  ( BC_face(1,1)%U0(i,j,k,2)*sinphi(i,j,l) -      &
                     BC_face(1,1)%U0(i,j,k,3)*cosphi(i,j,l) )**2 )
          enddo
       enddo
    enddo

    ! Non-centered derivatives in x-direction *sin(teta)*cos(phi)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    i=1
    do j=1,ngh
       do l=1,ngh
          duksi(i,j,l)= a06(1)*uf(1,j,l)+a06(2)*uf(2,j,l) &
                      + a06(3)*uf(3,j,l)+a06(4)*uf(4,j,l) &
                      + a06(5)*uf(5,j,l)+a06(6)*uf(6,j,l) &
                      + a06(7)*uf(7,j,l)
          dvksi(i,j,l)= a06(1)*vf(1,j,l)+a06(2)*vf(2,j,l) &
                      + a06(3)*vf(3,j,l)+a06(4)*vf(4,j,l) &
                      + a06(5)*vf(5,j,l)+a06(6)*vf(6,j,l) &
                      + a06(7)*vf(7,j,l)
          dwksi(i,j,l)= a06(1)*wf(1,j,l)+a06(2)*wf(2,j,l) &
                      + a06(3)*wf(3,j,l)+a06(4)*wf(4,j,l) &
                      + a06(5)*wf(5,j,l)+a06(6)*wf(6,j,l) &
                      + a06(7)*wf(7,j,l)
          dpksi(i,j,l)= a06(1)*pf(1,j,l)+a06(2)*pf(2,j,l) &
                      + a06(3)*pf(3,j,l)+a06(4)*pf(4,j,l) &
                      + a06(5)*pf(5,j,l)+a06(6)*pf(6,j,l) &
                      + a06(7)*pf(7,j,l)
          drksi(i,j,l)= a06(1)*rf(1,j,l)+a06(2)*rf(2,j,l) &
                      + a06(3)*rf(3,j,l)+a06(4)*rf(4,j,l) &
                      + a06(5)*rf(5,j,l)+a06(6)*rf(6,j,l) &
                      + a06(7)*rf(7,j,l)
       enddo
    enddo

    i=2
    do j=1,ngh
       do l=1,ngh
          duksi(i,j,l)= a15(1)*uf(1,j,l)+a15(2)*uf(2,j,l) &
                      + a15(3)*uf(3,j,l)+a15(4)*uf(4,j,l) &
                      + a15(5)*uf(5,j,l)+a15(6)*uf(6,j,l) &
                      + a15(7)*uf(7,j,l)
          dvksi(i,j,l)= a15(1)*vf(1,j,l)+a15(2)*vf(2,j,l) &
                      + a15(3)*vf(3,j,l)+a15(4)*vf(4,j,l) &
                      + a15(5)*vf(5,j,l)+a15(6)*vf(6,j,l) &
                      + a15(7)*vf(7,j,l)
          dwksi(i,j,l)= a15(1)*wf(1,j,l)+a15(2)*wf(2,j,l) &
                      + a15(3)*wf(3,j,l)+a15(4)*wf(4,j,l) &
                      + a15(5)*wf(5,j,l)+a15(6)*wf(6,j,l) &
                      + a15(7)*wf(7,j,l)
          dpksi(i,j,l)= a15(1)*pf(1,j,l)+a15(2)*pf(2,j,l) &
                      + a15(3)*pf(3,j,l)+a15(4)*pf(4,j,l) &
                      + a15(5)*pf(5,j,l)+a15(6)*pf(6,j,l) &
                      + a15(7)*pf(7,j,l)
          drksi(i,j,l)= a15(1)*rf(1,j,l)+a15(2)*rf(2,j,l) &
                      + a15(3)*rf(3,j,l)+a15(4)*rf(4,j,l) &
                      + a15(5)*rf(5,j,l)+a15(6)*rf(6,j,l) &
                      + a15(7)*rf(7,j,l)
       enddo
    enddo

    i=3
    do j=1,ngh
       do l=1,ngh
          duksi(i,j,l)= a24(1)*uf(1,j,l)+a24(2)*uf(2,j,l) &
                      + a24(3)*uf(3,j,l)+a24(4)*uf(4,j,l) &
                      + a24(5)*uf(5,j,l)+a24(6)*uf(6,j,l) &
                      + a24(7)*uf(7,j,l)
          dvksi(i,j,l)= a24(1)*vf(1,j,l)+a24(2)*vf(2,j,l) &
                      + a24(3)*vf(3,j,l)+a24(4)*vf(4,j,l) &
                      + a24(5)*vf(5,j,l)+a24(6)*vf(6,j,l) &
                      + a24(7)*vf(7,j,l)
          dwksi(i,j,l)= a24(1)*wf(1,j,l)+a24(2)*wf(2,j,l) &
                      + a24(3)*wf(3,j,l)+a24(4)*wf(4,j,l) &
                      + a24(5)*wf(5,j,l)+a24(6)*wf(6,j,l) &
                      + a24(7)*wf(7,j,l)
          dpksi(i,j,l)= a24(1)*pf(1,j,l)+a24(2)*pf(2,j,l) &
                      + a24(3)*pf(3,j,l)+a24(4)*pf(4,j,l) &
                      + a24(5)*pf(5,j,l)+a24(6)*pf(6,j,l) &
                      + a24(7)*pf(7,j,l)
          drksi(i,j,l)= a24(1)*rf(1,j,l)+a24(2)*rf(2,j,l) &
                      + a24(3)*rf(3,j,l)+a24(4)*rf(4,j,l) &
                      + a24(5)*rf(5,j,l)+a24(6)*rf(6,j,l) &
                      + a24(7)*rf(7,j,l)
       enddo
    enddo

    do i=4,ngh
       do j=1,ngh
          do l=1,ngh
             duksi(i,j,l)= a7(1)* (uf(i+1,j,l)-uf(i-1,j,l)) &
                         + a7(2)* (uf(i+2,j,l)-uf(i-2,j,l)) &
                         + a7(3)* (uf(i+3,j,l)-uf(i-3,j,l))
             dvksi(i,j,l)= a7(1)* (vf(i+1,j,l)-vf(i-1,j,l)) &
                         + a7(2)* (vf(i+2,j,l)-vf(i-2,j,l)) &
                         + a7(3)* (vf(i+3,j,l)-vf(i-3,j,l))
             dwksi(i,j,l)= a7(1)* (wf(i+1,j,l)-wf(i-1,j,l)) &
                         + a7(2)* (wf(i+2,j,l)-wf(i-2,j,l)) &
                         + a7(3)* (wf(i+3,j,l)-wf(i-3,j,l))
             dpksi(i,j,l)= a7(1)* (pf(i+1,j,l)-pf(i-1,j,l)) &
                         + a7(2)* (pf(i+2,j,l)-pf(i-2,j,l)) &
                         + a7(3)* (pf(i+3,j,l)-pf(i-3,j,l))
             drksi(i,j,l)= a7(1)* (rf(i+1,j,l)-rf(i-1,j,l)) &
                         + a7(2)* (rf(i+2,j,l)-rf(i-2,j,l)) &
                         + a7(3)* (rf(i+3,j,l)-rf(i-3,j,l))
          enddo
       enddo
    enddo

    ! Non-centered derivatives in y-direction *sin(teta)*sin(phi)
    ! =======================================
    j=1
    do i=1,ngh
       do l=1,ngh
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
    do i=1,ngh
       do l=1,ngh
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
    do i=1,ngh
       do l=1,ngh
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
       do i=1,ngh
          do l=1,ngh
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
       do i=1,ngh
          do j=1,ngh
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
    do i=1,ngh
       do j=1,ngh
          duphi(i,j,l)=  a42(1)*uf(i,j,l+2)+a42(2)*uf(i,j,l+1) &
                       + a42(3)*uf(i,j,l  )+a42(4)*uf(i,j,l-1) &
                       + a42(5)*uf(i,j,l-2)+a42(6)*uf(i,j,l-3) &
                       + a42(7)*uf(i,j,l-4)
          dvphi(i,j,l)=  a42(1)*vf(i,j,l+2)+a42(2)*vf(i,j,l+1) &
                       + a42(3)*vf(i,j,l  )+a42(4)*vf(i,j,l-1) &
                       + a42(5)*vf(i,j,l-2)+a42(6)*vf(i,j,l-3) &
                       + a42(7)*vf(i,j,l-4)
          dwphi(i,j,l)=  a42(1)*wf(i,j,l+2)+a42(2)*wf(i,j,l+1) &
                       + a42(3)*wf(i,j,l  )+a42(4)*wf(i,j,l-1) &
                       + a42(5)*wf(i,j,l-2)+a42(6)*wf(i,j,l-3) &
                       + a42(7)*wf(i,j,l-4)
          dpphi(i,j,l)=  a42(1)*pf(i,j,l+2)+a42(2)*pf(i,j,l+1) &
                       + a42(3)*pf(i,j,l  )+a42(4)*pf(i,j,l-1) &
                       + a42(5)*pf(i,j,l-2)+a42(6)*pf(i,j,l-3) &
                       + a42(7)*pf(i,j,l-4)
          drphi(i,j,l)=  a42(1)*rf(i,j,l+2)+a42(2)*rf(i,j,l+1) &
                       + a42(3)*rf(i,j,l  )+a42(4)*rf(i,j,l-1) &
                       + a42(5)*rf(i,j,l-2)+a42(6)*rf(i,j,l-3) &
                       + a42(7)*rf(i,j,l-4)
       enddo
    enddo

    k=nz-1
    l=k-nzmngh
    do i=1,ngh
       do j=1,ngh
          duphi(i,j,l)=  a51(1)*uf(i,j,l+1)+a51(2)*uf(i,j,l  ) &
                       + a51(3)*uf(i,j,l-1)+a51(4)*uf(i,j,l-2) &
                       + a51(5)*uf(i,j,l-3)+a51(6)*uf(i,j,l-4) &
                       + a51(7)*uf(i,j,l-5)
          dvphi(i,j,l)=  a51(1)*vf(i,j,l+1)+a51(2)*vf(i,j,l  ) &
                       + a51(3)*vf(i,j,l-1)+a51(4)*vf(i,j,l-2) &
                       + a51(5)*vf(i,j,l-3)+a51(6)*vf(i,j,l-4) &
                       + a51(7)*vf(i,j,l-5)
          dwphi(i,j,l)=  a51(1)*wf(i,j,l+1)+a51(2)*wf(i,j,l  ) &
                       + a51(3)*wf(i,j,l-1)+a51(4)*wf(i,j,l-2) &
                       + a51(5)*wf(i,j,l-3)+a51(6)*wf(i,j,l-4) &
                       + a51(7)*wf(i,j,l-5)
          dpphi(i,j,l)=  a51(1)*pf(i,j,l+1)+a51(2)*pf(i,j,l  ) &
                       + a51(3)*pf(i,j,l-1)+a51(4)*pf(i,j,l-2) &
                       + a51(5)*pf(i,j,l-3)+a51(6)*pf(i,j,l-4) &
                       + a51(7)*pf(i,j,l-5)
          drphi(i,j,l)=  a51(1)*rf(i,j,l+1)+a51(2)*rf(i,j,l  ) &
                       + a51(3)*rf(i,j,l-1)+a51(4)*rf(i,j,l-2) &
                       + a51(5)*rf(i,j,l-3)+a51(6)*rf(i,j,l-4) &
                       + a51(7)*rf(i,j,l-5)
       enddo
    enddo

    k=nz
    l=k-nzmngh
    do i=1,ngh
       do j=1,ngh
          duphi(i,j,l)=  a60(1)*uf(i,j,l  )+a60(2)*uf(i,j,l-1) &
                       + a60(3)*uf(i,j,l-2)+a60(4)*uf(i,j,l-3) &
                       + a60(5)*uf(i,j,l-4)+a60(6)*uf(i,j,l-5) &
                       + a60(7)*uf(i,j,l-6)
          dvphi(i,j,l)=  a60(1)*vf(i,j,l  )+a60(2)*vf(i,j,l-1) &
                       + a60(3)*vf(i,j,l-2)+a60(4)*vf(i,j,l-3) &
                       + a60(5)*vf(i,j,l-4)+a60(6)*vf(i,j,l-5) &
                       + a60(7)*vf(i,j,l-6)
          dwphi(i,j,l)=  a60(1)*wf(i,j,l  )+a60(2)*wf(i,j,l-1) &
                       + a60(3)*wf(i,j,l-2)+a60(4)*wf(i,j,l-3) &
                       + a60(5)*wf(i,j,l-4)+a60(6)*wf(i,j,l-5) &
                       + a60(7)*wf(i,j,l-6)
          dpphi(i,j,l)=  a60(1)*pf(i,j,l  )+a60(2)*pf(i,j,l-1) &
                       + a60(3)*pf(i,j,l-2)+a60(4)*pf(i,j,l-3) &
                       + a60(5)*pf(i,j,l-4)+a60(6)*pf(i,j,l-5) &
                       + a60(7)*pf(i,j,l-6)
          drphi(i,j,l)=  a60(1)*rf(i,j,l  )+a60(2)*rf(i,j,l-1) &
                       + a60(3)*rf(i,j,l-2)+a60(4)*rf(i,j,l-3) &
                       + a60(5)*rf(i,j,l-4)+a60(6)*rf(i,j,l-5) &
                       + a60(7)*rf(i,j,l-6)
       enddo
    enddo

    ! Compute rt,ut,vt,wt & pt
    ! ========================
    do i=1,ngh
       do j=1,ngh
          do k=nzmnghp1,nz
             l=k-nzmngh
             pt(i,j,l) = vg(i,j,l)*( &
                       ((dpksi(i,j,l)*ksi_x(i,j,k)+dpeta(i,j,l)*eta_x(i,j,k)+dpphi(i,j,l)*phi_x(i,j,k))*sintcosp(i,j,l) &
                       +(dpksi(i,j,l)*ksi_y(i,j,k)+dpeta(i,j,l)*eta_y(i,j,k)+dpphi(i,j,l)*phi_y(i,j,k))*sintsinp(i,j,l) &
                       +(dpksi(i,j,l)*ksi_z(i,j,k)+dpeta(i,j,l)*eta_z(i,j,k)+dpphi(i,j,l)*phi_z(i,j,k))*costeta(i,j,l)  &
                       )*ijacob3(i,j,k) + pf(i,j,l)*ir(i,j,l) )
             ut(i,j,l) = vg(i,j,l)*( &
                       ((duksi(i,j,l)*ksi_x(i,j,k)+dueta(i,j,l)*eta_x(i,j,k)+duphi(i,j,l)*phi_x(i,j,k))*sintcosp(i,j,l) &
                       +(duksi(i,j,l)*ksi_y(i,j,k)+dueta(i,j,l)*eta_y(i,j,k)+duphi(i,j,l)*phi_y(i,j,k))*sintsinp(i,j,l) &
                       +(duksi(i,j,l)*ksi_z(i,j,k)+dueta(i,j,l)*eta_z(i,j,k)+duphi(i,j,l)*phi_z(i,j,k))*costeta(i,j,l)  &
                       )*ijacob3(i,j,k) + uf(i,j,l)*ir(i,j,l) )
             vt(i,j,l) = vg(i,j,l)*( &
                       ((dvksi(i,j,l)*ksi_x(i,j,k)+dveta(i,j,l)*eta_x(i,j,k)+dvphi(i,j,l)*phi_x(i,j,k))*sintcosp(i,j,l) &
                       +(dvksi(i,j,l)*ksi_y(i,j,k)+dveta(i,j,l)*eta_y(i,j,k)+dvphi(i,j,l)*phi_y(i,j,k))*sintsinp(i,j,l) &
                       +(dvksi(i,j,l)*ksi_z(i,j,k)+dveta(i,j,l)*eta_z(i,j,k)+dvphi(i,j,l)*phi_z(i,j,k))*costeta(i,j,l)  &
                       )*ijacob3(i,j,k) + vf(i,j,l)*ir(i,j,l) )
             wt(i,j,l) = vg(i,j,l)*( &
                       ((dwksi(i,j,l)*ksi_x(i,j,k)+dweta(i,j,l)*eta_x(i,j,k)+dwphi(i,j,l)*phi_x(i,j,k))*sintcosp(i,j,l) &
                       +(dwksi(i,j,l)*ksi_y(i,j,k)+dweta(i,j,l)*eta_y(i,j,k)+dwphi(i,j,l)*phi_y(i,j,k))*sintsinp(i,j,l) &
                       +(dwksi(i,j,l)*ksi_z(i,j,k)+dweta(i,j,l)*eta_z(i,j,k)+dwphi(i,j,l)*phi_z(i,j,k))*costeta(i,j,l)  &
                       )*ijacob3(i,j,k) + wf(i,j,l)*ir(i,j,l) )
             rt(i,j,l) = vg(i,j,l)*( &
                       ((drksi(i,j,l)*ksi_x(i,j,k)+dreta(i,j,l)*eta_x(i,j,k)+drphi(i,j,l)*phi_x(i,j,k))*sintcosp(i,j,l) &
                       +(drksi(i,j,l)*ksi_y(i,j,k)+dreta(i,j,l)*eta_y(i,j,k)+drphi(i,j,l)*phi_y(i,j,k))*sintsinp(i,j,l) &
                       +(drksi(i,j,l)*ksi_z(i,j,k)+dreta(i,j,l)*eta_z(i,j,k)+drphi(i,j,l)*phi_z(i,j,k))*costeta(i,j,l)  &
                       )*ijacob3(i,j,k) + rf(i,j,l)*ir(i,j,l) )
          enddo
       enddo
    enddo

    ! Update fluxes at each RK step
    ! =============================
    do k=nzmnghp1,nz
       l=k-nzmngh
       do j=1,ngh
          do i=1,ngh
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

  end subroutine bc_TD3d_imin_jmin_kmax_c3

  !===============================================================================
  module subroutine bc_TD3d_imin_jmax_kmin_c3
  !===============================================================================
    !> 3D Tam & Dong's BC at imin-jmax-kmin (corner 1,2,1 /left-top-front)
    !> - 3D curvilinear version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l
    real(wp) :: cp,av,c2_
    real(wp), dimension(1:ngh+3,-2:ngh,1:ngh+3) :: rf,uf,vf,wf,pf
    real(wp), dimension(ngh,ngh,ngh) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(ngh,ngh,ngh) :: duksi,dvksi,dwksi,dpksi,drksi
    real(wp), dimension(ngh,ngh,ngh) :: dueta,dveta,dweta,dpeta,dreta
    real(wp), dimension(ngh,ngh,ngh) :: duphi,dvphi,dwphi,dpphi,drphi
    !-------------------------------------------------------------------------

    ! Pointers for spherical coordinates
    ! ==================================
    ir=>BC_corner(1,2,1)%i_r
    cosphi=>BC_corner(1,2,1)%cosp
    sinphi=>BC_corner(1,2,1)%sinp
    costeta=>BC_corner(1,2,1)%cost
    sinteta=>BC_corner(1,2,1)%sint
    costcosp=>BC_corner(1,2,1)%costcosp
    costsinp=>BC_corner(1,2,1)%costsinp
    sintcosp=>BC_corner(1,2,1)%sintcosp
    sintsinp=>BC_corner(1,2,1)%sintsinp

    ! Compute fluctuations
    ! ====================
    do i=1,nghp3
       do j=nymngh-2,ny
          l=j-nymngh
          do k=1,nghp3
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
    ! vg= u0*sin(teta)*cos(phi)+v0*sin(teta)*sin(phi)+w0*cos(teta)
    !   + sqrt{c0^2-[u0*cos(teta)*cos(phi)+v0*cos(teta)*sin(phi)-w0*sin(teta)]^2-[u0*sin(phi)-v0*cos(phi)]^2}
    do i=1,ngh
       do j=nymnghp1,ny
          l=j-nymngh
          do k=1,ngh
             vg(i,l,k)= BC_face(1,1)%U0(i,j,k,2)*sintcosp(i,l,k) +     &
               BC_face(1,1)%U0(i,j,k,3)*sintsinp(i,l,k) +     &
               BC_face(1,1)%U0(i,j,k,4)*costeta(i,l,k)          &
                      + sqrt( BC_face(1,1)%U0(i,j,k,6)-                &
                       ( BC_face(1,1)%U0(i,j,k,2)*costcosp(i,l,k) +    &
                   BC_face(1,1)%U0(i,j,k,3)*costsinp(i,l,k) -    &
                     BC_face(1,1)%U0(i,j,k,4)*sinteta(i,l,k) )**2- &
                  ( BC_face(1,1)%U0(i,j,k,2)*sinphi(i,l,k) -      &
                     BC_face(1,1)%U0(i,j,k,3)*cosphi(i,l,k) )**2 )
          enddo
       enddo
    enddo

    ! Non-centered derivatives in x-direction *sin(teta)*cos(phi)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    i=1
    do j=nymnghp1,ny
       l=j-nymngh
       do k=1,ngh
          duksi(i,l,k)= a06(1)*uf(1,l,k)+a06(2)*uf(2,l,k) &
                      + a06(3)*uf(3,l,k)+a06(4)*uf(4,l,k) &
                      + a06(5)*uf(5,l,k)+a06(6)*uf(6,l,k) &
                      + a06(7)*uf(7,l,k)
          dvksi(i,l,k)= a06(1)*vf(1,l,k)+a06(2)*vf(2,l,k) &
                      + a06(3)*vf(3,l,k)+a06(4)*vf(4,l,k) &
                      + a06(5)*vf(5,l,k)+a06(6)*vf(6,l,k) &
                      + a06(7)*vf(7,l,k)
          dwksi(i,l,k)= a06(1)*wf(1,l,k)+a06(2)*wf(2,l,k) &
                      + a06(3)*wf(3,l,k)+a06(4)*wf(4,l,k) &
                      + a06(5)*wf(5,l,k)+a06(6)*wf(6,l,k) &
                      + a06(7)*wf(7,l,k)
          dpksi(i,l,k)= a06(1)*pf(1,l,k)+a06(2)*pf(2,l,k) &
                      + a06(3)*pf(3,l,k)+a06(4)*pf(4,l,k) &
                      + a06(5)*pf(5,l,k)+a06(6)*pf(6,l,k) &
                      + a06(7)*pf(7,l,k)
          drksi(i,l,k)= a06(1)*rf(1,l,k)+a06(2)*rf(2,l,k) &
                      + a06(3)*rf(3,l,k)+a06(4)*rf(4,l,k) &
                      + a06(5)*rf(5,l,k)+a06(6)*rf(6,l,k) &
                      + a06(7)*rf(7,l,k)
       enddo
    enddo

    i=2
    do j=nymnghp1,ny
       l=j-nymngh
       do k=1,ngh
          duksi(i,l,k)= a15(1)*uf(1,l,k)+a15(2)*uf(2,l,k) &
                      + a15(3)*uf(3,l,k)+a15(4)*uf(4,l,k) &
                      + a15(5)*uf(5,l,k)+a15(6)*uf(6,l,k) &
                      + a15(7)*uf(7,l,k)
          dvksi(i,l,k)= a15(1)*vf(1,l,k)+a15(2)*vf(2,l,k) &
                      + a15(3)*vf(3,l,k)+a15(4)*vf(4,l,k) &
                      + a15(5)*vf(5,l,k)+a15(6)*vf(6,l,k) &
                      + a15(7)*vf(7,l,k)
          dwksi(i,l,k)= a15(1)*wf(1,l,k)+a15(2)*wf(2,l,k) &
                      + a15(3)*wf(3,l,k)+a15(4)*wf(4,l,k) &
                      + a15(5)*wf(5,l,k)+a15(6)*wf(6,l,k) &
                      + a15(7)*wf(7,l,k)
          dpksi(i,l,k)= a15(1)*pf(1,l,k)+a15(2)*pf(2,l,k) &
                      + a15(3)*pf(3,l,k)+a15(4)*pf(4,l,k) &
                      + a15(5)*pf(5,l,k)+a15(6)*pf(6,l,k) &
                      + a15(7)*pf(7,l,k)
          drksi(i,l,k)= a15(1)*rf(1,l,k)+a15(2)*rf(2,l,k) &
                      + a15(3)*rf(3,l,k)+a15(4)*rf(4,l,k) &
                      + a15(5)*rf(5,l,k)+a15(6)*rf(6,l,k) &
                      + a15(7)*rf(7,l,k)
       enddo
    enddo

    i=3
    do j=nymnghp1,ny
       l=j-nymngh
       do k=1,ngh
          duksi(i,l,k)= a24(1)*uf(1,l,k)+a24(2)*uf(2,l,k) &
                      + a24(3)*uf(3,l,k)+a24(4)*uf(4,l,k) &
                      + a24(5)*uf(5,l,k)+a24(6)*uf(6,l,k) &
                      + a24(7)*uf(7,l,k)
          dvksi(i,l,k)= a24(1)*vf(1,l,k)+a24(2)*vf(2,l,k) &
                      + a24(3)*vf(3,l,k)+a24(4)*vf(4,l,k) &
                      + a24(5)*vf(5,l,k)+a24(6)*vf(6,l,k) &
                      + a24(7)*vf(7,l,k)
          dwksi(i,l,k)= a24(1)*wf(1,l,k)+a24(2)*wf(2,l,k) &
                      + a24(3)*wf(3,l,k)+a24(4)*wf(4,l,k) &
                      + a24(5)*wf(5,l,k)+a24(6)*wf(6,l,k) &
                      + a24(7)*wf(7,l,k)
          dpksi(i,l,k)= a24(1)*pf(1,l,k)+a24(2)*pf(2,l,k) &
                      + a24(3)*pf(3,l,k)+a24(4)*pf(4,l,k) &
                      + a24(5)*pf(5,l,k)+a24(6)*pf(6,l,k) &
                      + a24(7)*pf(7,l,k)
          drksi(i,l,k)= a24(1)*rf(1,l,k)+a24(2)*rf(2,l,k) &
                      + a24(3)*rf(3,l,k)+a24(4)*rf(4,l,k) &
                      + a24(5)*rf(5,l,k)+a24(6)*rf(6,l,k) &
                      + a24(7)*rf(7,l,k)
       enddo
    enddo

    do i=4,ngh
       do j=nymnghp1,ny
          l=j-nymngh
          do k=1,ngh
             duksi(i,l,k)= a7(1)* (uf(i+1,l,k)-uf(i-1,l,k)) &
                         + a7(2)* (uf(i+2,l,k)-uf(i-2,l,k)) &
                         + a7(3)* (uf(i+3,l,k)-uf(i-3,l,k))
             dvksi(i,l,k)= a7(1)* (vf(i+1,l,k)-vf(i-1,l,k)) &
                         + a7(2)* (vf(i+2,l,k)-vf(i-2,l,k)) &
                         + a7(3)* (vf(i+3,l,k)-vf(i-3,l,k))
             dwksi(i,l,k)= a7(1)* (wf(i+1,l,k)-wf(i-1,l,k)) &
                         + a7(2)* (wf(i+2,l,k)-wf(i-2,l,k)) &
                         + a7(3)* (wf(i+3,l,k)-wf(i-3,l,k))
             dpksi(i,l,k)= a7(1)* (pf(i+1,l,k)-pf(i-1,l,k)) &
                         + a7(2)* (pf(i+2,l,k)-pf(i-2,l,k)) &
                         + a7(3)* (pf(i+3,l,k)-pf(i-3,l,k))
             drksi(i,l,k)= a7(1)* (rf(i+1,l,k)-rf(i-1,l,k)) &
                         + a7(2)* (rf(i+2,l,k)-rf(i-2,l,k)) &
                         + a7(3)* (rf(i+3,l,k)-rf(i-3,l,k))
          enddo
       enddo
    enddo

    ! Non-centered derivatives in y-direction *sin(teta)*sin(phi)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    do j=nymnghp1,ny-3
       l=j-nymngh
       do i=1,ngh
          do k=1,ngh
             dueta(i,l,k)=  a7(1)*(uf(i,l+1,k)-uf(i,l-1,k)) &
                          + a7(2)*(uf(i,l+2,k)-uf(i,l-2,k)) &
                          + a7(3)*(uf(i,l+3,k)-uf(i,l-3,k))
             dveta(i,l,k)=  a7(1)*(vf(i,l+1,k)-vf(i,l-1,k)) &
                          + a7(2)*(vf(i,l+2,k)-vf(i,l-2,k)) &
                          + a7(3)*(vf(i,l+3,k)-vf(i,l-3,k))
             dweta(i,l,k)=  a7(1)*(wf(i,l+1,k)-wf(i,l-1,k)) &
                          + a7(2)*(wf(i,l+2,k)-wf(i,l-2,k)) &
                          + a7(3)*(wf(i,l+3,k)-wf(i,l-3,k))
             dpeta(i,l,k)=  a7(1)*(pf(i,l+1,k)-pf(i,l-1,k)) &
                          + a7(2)*(pf(i,l+2,k)-pf(i,l-2,k)) &
                          + a7(3)*(pf(i,l+3,k)-pf(i,l-3,k))
             dreta(i,l,k)=  a7(1)*(rf(i,l+1,k)-rf(i,l-1,k)) &
                          + a7(2)*(rf(i,l+2,k)-rf(i,l-2,k)) &
                          + a7(3)*(rf(i,l+3,k)-rf(i,l-3,k))
          enddo
       enddo
    enddo

    j=ny-2
    l=j-nymngh
    do i=1,ngh
       do k=1,ngh
          dueta(i,l,k)=  a42(1)*uf(i,l+2,k)+a42(2)*uf(i,l+1,k) &
                       + a42(3)*uf(i,l  ,k)+a42(4)*uf(i,l-1,k) &
                       + a42(5)*uf(i,l-2,k)+a42(6)*uf(i,l-3,k) &
                       + a42(7)*uf(i,l-4,k)
          dveta(i,l,k)=  a42(1)*vf(i,l+2,k)+a42(2)*vf(i,l+1,k) &
                       + a42(3)*vf(i,l  ,k)+a42(4)*vf(i,l-1,k) &
                       + a42(5)*vf(i,l-2,k)+a42(6)*vf(i,l-3,k) &
                       + a42(7)*vf(i,l-4,k)
          dweta(i,l,k)=  a42(1)*wf(i,l+2,k)+a42(2)*wf(i,l+1,k) &
                       + a42(3)*wf(i,l  ,k)+a42(4)*wf(i,l-1,k) &
                       + a42(5)*wf(i,l-2,k)+a42(6)*wf(i,l-3,k) &
                       + a42(7)*wf(i,l-4,k)
          dpeta(i,l,k)=  a42(1)*pf(i,l+2,k)+a42(2)*pf(i,l+1,k) &
                       + a42(3)*pf(i,l  ,k)+a42(4)*pf(i,l-1,k) &
                       + a42(5)*pf(i,l-2,k)+a42(6)*pf(i,l-3,k) &
                       + a42(7)*pf(i,l-4,k)
          dreta(i,l,k)=  a42(1)*rf(i,l+2,k)+a42(2)*rf(i,l+1,k) &
                       + a42(3)*rf(i,l  ,k)+a42(4)*rf(i,l-1,k) &
                       + a42(5)*rf(i,l-2,k)+a42(6)*rf(i,l-3,k) &
                       + a42(7)*rf(i,l-4,k)
       enddo
    enddo

    j=ny-1
    l=j-nymngh
    do i=1,ngh
       do k=1,ngh
          dueta(i,l,k)=  a51(1)*uf(i,l+1,k)+a51(2)*uf(i,l  ,k) &
                       + a51(3)*uf(i,l-1,k)+a51(4)*uf(i,l-2,k) &
                       + a51(5)*uf(i,l-3,k)+a51(6)*uf(i,l-4,k) &
                       + a51(7)*uf(i,l-5,k)
          dveta(i,l,k)=  a51(1)*vf(i,l+1,k)+a51(2)*vf(i,l  ,k) &
                       + a51(3)*vf(i,l-1,k)+a51(4)*vf(i,l-2,k) &
                       + a51(5)*vf(i,l-3,k)+a51(6)*vf(i,l-4,k) &
                       + a51(7)*vf(i,l-5,k)
          dweta(i,l,k)=  a51(1)*wf(i,l+1,k)+a51(2)*wf(i,l  ,k) &
                       + a51(3)*wf(i,l-1,k)+a51(4)*wf(i,l-2,k) &
                       + a51(5)*wf(i,l-3,k)+a51(6)*wf(i,l-4,k) &
                       + a51(7)*wf(i,l-5,k)
          dpeta(i,l,k)=  a51(1)*pf(i,l+1,k)+a51(2)*pf(i,l  ,k) &
                       + a51(3)*pf(i,l-1,k)+a51(4)*pf(i,l-2,k) &
                       + a51(5)*pf(i,l-3,k)+a51(6)*pf(i,l-4,k) &
                       + a51(7)*pf(i,l-5,k)
          dreta(i,l,k)=  a51(1)*rf(i,l+1,k)+a51(2)*rf(i,l  ,k) &
                       + a51(3)*rf(i,l-1,k)+a51(4)*rf(i,l-2,k) &
                       + a51(5)*rf(i,l-3,k)+a51(6)*rf(i,l-4,k) &
                       + a51(7)*rf(i,l-5,k)
       enddo
    enddo

    j=ny
    l=j-nymngh
    do i=1,ngh
       do k=1,ngh
          dueta(i,l,k)=  a60(1)*uf(i,l  ,k)+a60(2)*uf(i,l-1,k) &
                       + a60(3)*uf(i,l-2,k)+a60(4)*uf(i,l-3,k) &
                       + a60(5)*uf(i,l-4,k)+a60(6)*uf(i,l-5,k) &
                       + a60(7)*uf(i,l-6,k)
          dveta(i,l,k)=  a60(1)*vf(i,l  ,k)+a60(2)*vf(i,l-1,k) &
                       + a60(3)*vf(i,l-2,k)+a60(4)*vf(i,l-3,k) &
                       + a60(5)*vf(i,l-4,k)+a60(6)*vf(i,l-5,k) &
                       + a60(7)*vf(i,l-6,k)
          dweta(i,l,k)=  a60(1)*wf(i,l  ,k)+a60(2)*wf(i,l-1,k) &
                       + a60(3)*wf(i,l-2,k)+a60(4)*wf(i,l-3,k) &
                       + a60(5)*wf(i,l-4,k)+a60(6)*wf(i,l-5,k) &
                       + a60(7)*wf(i,l-6,k)
          dpeta(i,l,k)=  a60(1)*pf(i,l  ,k)+a60(2)*pf(i,l-1,k) &
                       + a60(3)*pf(i,l-2,k)+a60(4)*pf(i,l-3,k) &
                       + a60(5)*pf(i,l-4,k)+a60(6)*pf(i,l-5,k) &
                       + a60(7)*pf(i,l-6,k)
          dreta(i,l,k)=  a60(1)*rf(i,l  ,k)+a60(2)*rf(i,l-1,k) &
                       + a60(3)*rf(i,l-2,k)+a60(4)*rf(i,l-3,k) &
                       + a60(5)*rf(i,l-4,k)+a60(6)*rf(i,l-5,k) &
                       + a60(7)*rf(i,l-6,k)
       enddo
    enddo

    ! Non-centered derivatives in z-direction *cos(teta)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    k=1
    do i=1,ngh
       do l=1,ngh
          duphi(i,l,k)=  a06(1)*uf(i,l,1)+a06(2)*uf(i,l,2) &
                       + a06(3)*uf(i,l,3)+a06(4)*uf(i,l,4) &
                       + a06(5)*uf(i,l,5)+a06(6)*uf(i,l,6) &
                       + a06(7)*uf(i,l,7)
          dvphi(i,l,k)=  a06(1)*vf(i,l,1)+a06(2)*vf(i,l,2) &
                       + a06(3)*vf(i,l,3)+a06(4)*vf(i,l,4) &
                       + a06(5)*vf(i,l,5)+a06(6)*vf(i,l,6) &
                       + a06(7)*vf(i,l,7)
          dwphi(i,l,k)=  a06(1)*wf(i,l,1)+a06(2)*wf(i,l,2) &
                       + a06(3)*wf(i,l,3)+a06(4)*wf(i,l,4) &
                       + a06(5)*wf(i,l,5)+a06(6)*wf(i,l,6) &
                       + a06(7)*wf(i,l,7)
          dpphi(i,l,k)=  a06(1)*pf(i,l,1)+a06(2)*pf(i,l,2) &
                       + a06(3)*pf(i,l,3)+a06(4)*pf(i,l,4) &
                       + a06(5)*pf(i,l,5)+a06(6)*pf(i,l,6) &
                       + a06(7)*pf(i,l,7)
          drphi(i,l,k)=  a06(1)*rf(i,l,1)+a06(2)*rf(i,l,2) &
                       + a06(3)*rf(i,l,3)+a06(4)*rf(i,l,4) &
                       + a06(5)*rf(i,l,5)+a06(6)*rf(i,l,6) &
                       + a06(7)*rf(i,l,7)
       enddo
    enddo

    k=2
    do i=1,ngh
       do l=1,ngh
          duphi(i,l,k)=  a15(1)*uf(i,l,1)+a15(2)*uf(i,l,2) &
                       + a15(3)*uf(i,l,3)+a15(4)*uf(i,l,4) &
                       + a15(5)*uf(i,l,5)+a15(6)*uf(i,l,6) &
                       + a15(7)*uf(i,l,7)
          dvphi(i,l,k)=  a15(1)*vf(i,l,1)+a15(2)*vf(i,l,2) &
                       + a15(3)*vf(i,l,3)+a15(4)*vf(i,l,4) &
                       + a15(5)*vf(i,l,5)+a15(6)*vf(i,l,6) &
                       + a15(7)*vf(i,l,7)
          dwphi(i,l,k)=  a15(1)*wf(i,l,1)+a15(2)*wf(i,l,2) &
                       + a15(3)*wf(i,l,3)+a15(4)*wf(i,l,4) &
                       + a15(5)*wf(i,l,5)+a15(6)*wf(i,l,6) &
                       + a15(7)*wf(i,l,7)
          dpphi(i,l,k)=  a15(1)*pf(i,l,1)+a15(2)*pf(i,l,2) &
                       + a15(3)*pf(i,l,3)+a15(4)*pf(i,l,4) &
                       + a15(5)*pf(i,l,5)+a15(6)*pf(i,l,6) &
                       + a15(7)*pf(i,l,7)
          drphi(i,l,k)=  a15(1)*rf(i,l,1)+a15(2)*rf(i,l,2) &
                       + a15(3)*rf(i,l,3)+a15(4)*rf(i,l,4) &
                       + a15(5)*rf(i,l,5)+a15(6)*rf(i,l,6) &
                       + a15(7)*rf(i,l,7)
       enddo
    enddo

    k=3
    do i=1,ngh
       do l=1,ngh
          duphi(i,l,k)=  a24(1)*uf(i,l,1)+a24(2)*uf(i,l,2) &
                       + a24(3)*uf(i,l,3)+a24(4)*uf(i,l,4) &
                       + a24(5)*uf(i,l,5)+a24(6)*uf(i,l,6) &
                       + a24(7)*uf(i,l,7)
          dvphi(i,l,k)=  a24(1)*vf(i,l,1)+a24(2)*vf(i,l,2) &
                       + a24(3)*vf(i,l,3)+a24(4)*vf(i,l,4) &
                       + a24(5)*vf(i,l,5)+a24(6)*vf(i,l,6) &
                       + a24(7)*vf(i,l,7)
          dwphi(i,l,k)=  a24(1)*wf(i,l,1)+a24(2)*wf(i,l,2) &
                       + a24(3)*wf(i,l,3)+a24(4)*wf(i,l,4) &
                       + a24(5)*wf(i,l,5)+a24(6)*wf(i,l,6) &
                       + a24(7)*wf(i,l,7)
          dpphi(i,l,k)=  a24(1)*pf(i,l,1)+a24(2)*pf(i,l,2) &
                       + a24(3)*pf(i,l,3)+a24(4)*pf(i,l,4) &
                       + a24(5)*pf(i,l,5)+a24(6)*pf(i,l,6) &
                       + a24(7)*pf(i,l,7)
          drphi(i,l,k)=  a24(1)*rf(i,l,1)+a24(2)*rf(i,l,2) &
                       + a24(3)*rf(i,l,3)+a24(4)*rf(i,l,4) &
                       + a24(5)*rf(i,l,5)+a24(6)*rf(i,l,6) &
                       + a24(7)*rf(i,l,7)
       enddo
    enddo

    do k=4,ngh
       do i=1,ngh
          do l=1,ngh
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
    do i=1,ngh
       do j=nymnghp1,ny
          l=j-nymngh
          do k=1,ngh
             pt(i,l,k) = vg(i,l,k)*( &
                       ((dpksi(i,l,k)*ksi_x(i,j,k)+dpeta(i,l,k)*eta_x(i,j,k)+dpphi(i,l,k)*phi_x(i,j,k))*sintcosp(i,l,k) &
                       +(dpksi(i,l,k)*ksi_y(i,j,k)+dpeta(i,l,k)*eta_y(i,j,k)+dpphi(i,l,k)*phi_y(i,j,k))*sintsinp(i,l,k) &
                       +(dpksi(i,l,k)*ksi_z(i,j,k)+dpeta(i,l,k)*eta_z(i,j,k)+dpphi(i,l,k)*phi_z(i,j,k))*costeta(i,l,k)  &
                       )*ijacob3(i,j,k) + pf(i,l,k)*ir(i,l,k) )
             ut(i,l,k) = vg(i,l,k)*( &
                       ((duksi(i,l,k)*ksi_x(i,j,k)+dueta(i,l,k)*eta_x(i,j,k)+duphi(i,l,k)*phi_x(i,j,k))*sintcosp(i,l,k) &
                       +(duksi(i,l,k)*ksi_y(i,j,k)+dueta(i,l,k)*eta_y(i,j,k)+duphi(i,l,k)*phi_y(i,j,k))*sintsinp(i,l,k) &
                       +(duksi(i,l,k)*ksi_z(i,j,k)+dueta(i,l,k)*eta_z(i,j,k)+duphi(i,l,k)*phi_z(i,j,k))*costeta(i,l,k)  &
                       )*ijacob3(i,j,k) + uf(i,l,k)*ir(i,l,k) )
             vt(i,l,k) = vg(i,l,k)*( &
                       ((dvksi(i,l,k)*ksi_x(i,j,k)+dveta(i,l,k)*eta_x(i,j,k)+dvphi(i,l,k)*phi_x(i,j,k))*sintcosp(i,l,k) &
                       +(dvksi(i,l,k)*ksi_y(i,j,k)+dveta(i,l,k)*eta_y(i,j,k)+dvphi(i,l,k)*phi_y(i,j,k))*sintsinp(i,l,k) &
                       +(dvksi(i,l,k)*ksi_z(i,j,k)+dveta(i,l,k)*eta_z(i,j,k)+dvphi(i,l,k)*phi_z(i,j,k))*costeta(i,l,k)  &
                       )*ijacob3(i,j,k) + vf(i,l,k)*ir(i,l,k) )
             wt(i,l,k) = vg(i,l,k)*( &
                       ((dwksi(i,l,k)*ksi_x(i,j,k)+dweta(i,l,k)*eta_x(i,j,k)+dwphi(i,l,k)*phi_x(i,j,k))*sintcosp(i,l,k) &
                       +(dwksi(i,l,k)*ksi_y(i,j,k)+dweta(i,l,k)*eta_y(i,j,k)+dwphi(i,l,k)*phi_y(i,j,k))*sintsinp(i,l,k) &
                       +(dwksi(i,l,k)*ksi_z(i,j,k)+dweta(i,l,k)*eta_z(i,j,k)+dwphi(i,l,k)*phi_z(i,j,k))*costeta(i,l,k)  &
                       )*ijacob3(i,j,k) + wf(i,l,k)*ir(i,l,k) )
             rt(i,l,k) = vg(i,l,k)*( &
                       ((drksi(i,l,k)*ksi_x(i,j,k)+dreta(i,l,k)*eta_x(i,j,k)+drphi(i,l,k)*phi_x(i,j,k))*sintcosp(i,l,k) &
                       +(drksi(i,l,k)*ksi_y(i,j,k)+dreta(i,l,k)*eta_y(i,j,k)+drphi(i,l,k)*phi_y(i,j,k))*sintsinp(i,l,k) &
                       +(drksi(i,l,k)*ksi_z(i,j,k)+dreta(i,l,k)*eta_z(i,j,k)+drphi(i,l,k)*phi_z(i,j,k))*costeta(i,l,k)  &
                       )*ijacob3(i,j,k) + rf(i,l,k)*ir(i,l,k) )
          enddo
       enddo
    enddo

    ! Update fluxes at each RK step
    ! =============================
    do k=1,ngh
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

  end subroutine bc_TD3d_imin_jmax_kmin_c3

  !===============================================================================
  module subroutine bc_TD3d_imin_jmax_kmax_c3
  !===============================================================================
    !> 3D Tam & Dong's BC at imin-jmax-kmax (corner 1,2,2 /left-top-back)
    !> - 3D curvilinear version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l,m
    real(wp) :: cp,av,c2_
    real(wp), dimension(1:ngh+3,-2:ngh,-2:ngh) :: rf,uf,vf,wf,pf
    real(wp), dimension(ngh,ngh,ngh) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(ngh,ngh,ngh) :: duksi,dvksi,dwksi,dpksi,drksi
    real(wp), dimension(ngh,ngh,ngh) :: dueta,dveta,dweta,dpeta,dreta
    real(wp), dimension(ngh,ngh,ngh) :: duphi,dvphi,dwphi,dpphi,drphi
    !-------------------------------------------------------------------------

    ! Pointers for spherical coordinates
    ! ==================================
    ir=>BC_corner(1,2,2)%i_r
    cosphi=>BC_corner(1,2,2)%cosp
    sinphi=>BC_corner(1,2,2)%sinp
    costeta=>BC_corner(1,2,2)%cost
    sinteta=>BC_corner(1,2,2)%sint
    costcosp=>BC_corner(1,2,2)%costcosp
    costsinp=>BC_corner(1,2,2)%costsinp
    sintcosp=>BC_corner(1,2,2)%sintcosp
    sintsinp=>BC_corner(1,2,2)%sintsinp

    ! Compute fluctuations
    ! ====================
    do i=1,nghp3
       do j=nymngh-2,ny
          l=j-nymngh
          do k=nzmngh-2,nz
             m=k-nzmngh
             rf(i,l,m)=rho_n(i,j,k)-BC_face(1,1)%U0(i,j,k,1)
             uf(i,l,m)=   uu(i,j,k)-BC_face(1,1)%U0(i,j,k,2)
             vf(i,l,m)=   vv(i,j,k)-BC_face(1,1)%U0(i,j,k,3)
             wf(i,l,m)=   ww(i,j,k)-BC_face(1,1)%U0(i,j,k,4)
             pf(i,l,m)=  prs(i,j,k)-BC_face(1,1)%U0(i,j,k,5)
          enddo
       enddo
    enddo

    ! Compute group velocity vg
    ! =========================
    ! vg= u0*sin(teta)*cos(phi)+v0*sin(teta)*sin(phi)+w0*cos(teta)
    !   + sqrt{c0^2-[u0*cos(teta)*cos(phi)+v0*cos(teta)*sin(phi)-w0*sin(teta)]^2-[u0*sin(phi)-v0*cos(phi)]^2}
    do i=1,ngh
       do j=nymnghp1,ny
          l=j-nymngh
          do k=nzmnghp1,nz
             m=k-nzmngh
             vg(i,l,m)= BC_face(1,1)%U0(i,j,k,2)*sintcosp(i,l,m) +     &
               BC_face(1,1)%U0(i,j,k,3)*sintsinp(i,l,m) +     &
               BC_face(1,1)%U0(i,j,k,4)*costeta(i,l,m)          &
                      + sqrt( BC_face(1,1)%U0(i,j,k,6)-                &
                       ( BC_face(1,1)%U0(i,j,k,2)*costcosp(i,l,m) +    &
                   BC_face(1,1)%U0(i,j,k,3)*costsinp(i,l,m) -    &
                     BC_face(1,1)%U0(i,j,k,4)*sinteta(i,l,m) )**2- &
                  ( BC_face(1,1)%U0(i,j,k,2)*sinphi(i,l,m) -      &
                     BC_face(1,1)%U0(i,j,k,3)*cosphi(i,l,m) )**2 )
          enddo
       enddo
    enddo

    ! Non-centered derivatives in x-direction *sin(teta)*cos(phi)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    i=1
    do j=nymnghp1,ny
       l=j-nymngh
       do m=1,ngh
          duksi(i,l,m)= a06(1)*uf(1,l,m)+a06(2)*uf(2,l,m) &
                      + a06(3)*uf(3,l,m)+a06(4)*uf(4,l,m) &
                      + a06(5)*uf(5,l,m)+a06(6)*uf(6,l,m) &
                      + a06(7)*uf(7,l,m)
          dvksi(i,l,m)= a06(1)*vf(1,l,m)+a06(2)*vf(2,l,m) &
                      + a06(3)*vf(3,l,m)+a06(4)*vf(4,l,m) &
                      + a06(5)*vf(5,l,m)+a06(6)*vf(6,l,m) &
                      + a06(7)*vf(7,l,m)
          dwksi(i,l,m)= a06(1)*wf(1,l,m)+a06(2)*wf(2,l,m) &
                      + a06(3)*wf(3,l,m)+a06(4)*wf(4,l,m) &
                      + a06(5)*wf(5,l,m)+a06(6)*wf(6,l,m) &
                      + a06(7)*wf(7,l,m)
          dpksi(i,l,m)= a06(1)*pf(1,l,m)+a06(2)*pf(2,l,m) &
                      + a06(3)*pf(3,l,m)+a06(4)*pf(4,l,m) &
                      + a06(5)*pf(5,l,m)+a06(6)*pf(6,l,m) &
                      + a06(7)*pf(7,l,m)
          drksi(i,l,m)= a06(1)*rf(1,l,m)+a06(2)*rf(2,l,m) &
                      + a06(3)*rf(3,l,m)+a06(4)*rf(4,l,m) &
                      + a06(5)*rf(5,l,m)+a06(6)*rf(6,l,m) &
                      + a06(7)*rf(7,l,m)
       enddo
    enddo

    i=2
    do j=nymnghp1,ny
       l=j-nymngh
       do m=1,ngh
          duksi(i,l,m)= a15(1)*uf(1,l,m)+a15(2)*uf(2,l,m) &
                      + a15(3)*uf(3,l,m)+a15(4)*uf(4,l,m) &
                      + a15(5)*uf(5,l,m)+a15(6)*uf(6,l,m) &
                      + a15(7)*uf(7,l,m)
          dvksi(i,l,m)= a15(1)*vf(1,l,m)+a15(2)*vf(2,l,m) &
                      + a15(3)*vf(3,l,m)+a15(4)*vf(4,l,m) &
                      + a15(5)*vf(5,l,m)+a15(6)*vf(6,l,m) &
                      + a15(7)*vf(7,l,m)
          dwksi(i,l,m)= a15(1)*wf(1,l,m)+a15(2)*wf(2,l,m) &
                      + a15(3)*wf(3,l,m)+a15(4)*wf(4,l,m) &
                      + a15(5)*wf(5,l,m)+a15(6)*wf(6,l,m) &
                      + a15(7)*wf(7,l,m)
          dpksi(i,l,m)= a15(1)*pf(1,l,m)+a15(2)*pf(2,l,m) &
                      + a15(3)*pf(3,l,m)+a15(4)*pf(4,l,m) &
                      + a15(5)*pf(5,l,m)+a15(6)*pf(6,l,m) &
                      + a15(7)*pf(7,l,m)
          drksi(i,l,m)= a15(1)*rf(1,l,m)+a15(2)*rf(2,l,m) &
                      + a15(3)*rf(3,l,m)+a15(4)*rf(4,l,m) &
                      + a15(5)*rf(5,l,m)+a15(6)*rf(6,l,m) &
                      + a15(7)*rf(7,l,m)
       enddo
    enddo

    i=3
    do j=nymnghp1,ny
       l=j-nymngh
       do m=1,ngh
          duksi(i,l,m)= a24(1)*uf(1,l,m)+a24(2)*uf(2,l,m) &
                      + a24(3)*uf(3,l,m)+a24(4)*uf(4,l,m) &
                      + a24(5)*uf(5,l,m)+a24(6)*uf(6,l,m) &
                      + a24(7)*uf(7,l,m)
          dvksi(i,l,m)= a24(1)*vf(1,l,m)+a24(2)*vf(2,l,m) &
                      + a24(3)*vf(3,l,m)+a24(4)*vf(4,l,m) &
                      + a24(5)*vf(5,l,m)+a24(6)*vf(6,l,m) &
                      + a24(7)*vf(7,l,m)
          dwksi(i,l,m)= a24(1)*wf(1,l,m)+a24(2)*wf(2,l,m) &
                      + a24(3)*wf(3,l,m)+a24(4)*wf(4,l,m) &
                      + a24(5)*wf(5,l,m)+a24(6)*wf(6,l,m) &
                      + a24(7)*wf(7,l,m)
          dpksi(i,l,m)= a24(1)*pf(1,l,m)+a24(2)*pf(2,l,m) &
                      + a24(3)*pf(3,l,m)+a24(4)*pf(4,l,m) &
                      + a24(5)*pf(5,l,m)+a24(6)*pf(6,l,m) &
                      + a24(7)*pf(7,l,m)
          drksi(i,l,m)= a24(1)*rf(1,l,m)+a24(2)*rf(2,l,m) &
                      + a24(3)*rf(3,l,m)+a24(4)*rf(4,l,m) &
                      + a24(5)*rf(5,l,m)+a24(6)*rf(6,l,m) &
                      + a24(7)*rf(7,l,m)
       enddo
    enddo

    do i=4,ngh
       do j=nymnghp1,ny
          l=j-nymngh
          do m=1,ngh
             duksi(i,l,m)= a7(1)* (uf(i+1,l,m)-uf(i-1,l,m)) &
                         + a7(2)* (uf(i+2,l,m)-uf(i-2,l,m)) &
                         + a7(3)* (uf(i+3,l,m)-uf(i-3,l,m))
             dvksi(i,l,m)= a7(1)* (vf(i+1,l,m)-vf(i-1,l,m)) &
                         + a7(2)* (vf(i+2,l,m)-vf(i-2,l,m)) &
                         + a7(3)* (vf(i+3,l,m)-vf(i-3,l,m))
             dwksi(i,l,m)= a7(1)* (wf(i+1,l,m)-wf(i-1,l,m)) &
                         + a7(2)* (wf(i+2,l,m)-wf(i-2,l,m)) &
                         + a7(3)* (wf(i+3,l,m)-wf(i-3,l,m))
             dpksi(i,l,m)= a7(1)* (pf(i+1,l,m)-pf(i-1,l,m)) &
                         + a7(2)* (pf(i+2,l,m)-pf(i-2,l,m)) &
                         + a7(3)* (pf(i+3,l,m)-pf(i-3,l,m))
             drksi(i,l,m)= a7(1)* (rf(i+1,l,m)-rf(i-1,l,m)) &
                         + a7(2)* (rf(i+2,l,m)-rf(i-2,l,m)) &
                         + a7(3)* (rf(i+3,l,m)-rf(i-3,l,m))
          enddo
       enddo
    enddo

    ! Non-centered derivatives in y-direction *sin(teta)*sin(phi)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    do j=nymnghp1,ny-3
       l=j-nymngh
       do i=1,ngh
          do m=1,ngh
             dueta(i,l,m)=  a7(1)*(uf(i,l+1,m)-uf(i,l-1,m)) &
                          + a7(2)*(uf(i,l+2,m)-uf(i,l-2,m)) &
                          + a7(3)*(uf(i,l+3,m)-uf(i,l-3,m))
             dveta(i,l,m)=  a7(1)*(vf(i,l+1,m)-vf(i,l-1,m)) &
                          + a7(2)*(vf(i,l+2,m)-vf(i,l-2,m)) &
                          + a7(3)*(vf(i,l+3,m)-vf(i,l-3,m))
             dweta(i,l,m)=  a7(1)*(wf(i,l+1,m)-wf(i,l-1,m)) &
                          + a7(2)*(wf(i,l+2,m)-wf(i,l-2,m)) &
                          + a7(3)*(wf(i,l+3,m)-wf(i,l-3,m))
             dpeta(i,l,m)=  a7(1)*(pf(i,l+1,m)-pf(i,l-1,m)) &
                          + a7(2)*(pf(i,l+2,m)-pf(i,l-2,m)) &
                          + a7(3)*(pf(i,l+3,m)-pf(i,l-3,m))
             dreta(i,l,m)=  a7(1)*(rf(i,l+1,m)-rf(i,l-1,m)) &
                          + a7(2)*(rf(i,l+2,m)-rf(i,l-2,m)) &
                          + a7(3)*(rf(i,l+3,m)-rf(i,l-3,m))
          enddo
       enddo
    enddo

    j=ny-2
    l=j-nymngh
    do i=1,ngh
       do m=1,ngh
          dueta(i,l,m)=  a42(1)*uf(i,l+2,m)+a42(2)*uf(i,l+1,m) &
                       + a42(3)*uf(i,l  ,m)+a42(4)*uf(i,l-1,m) &
                       + a42(5)*uf(i,l-2,m)+a42(6)*uf(i,l-3,m) &
                       + a42(7)*uf(i,l-4,m)
          dveta(i,l,m)=  a42(1)*vf(i,l+2,m)+a42(2)*vf(i,l+1,m) &
                       + a42(3)*vf(i,l  ,m)+a42(4)*vf(i,l-1,m) &
                       + a42(5)*vf(i,l-2,m)+a42(6)*vf(i,l-3,m) &
                       + a42(7)*vf(i,l-4,m)
          dweta(i,l,m)=  a42(1)*wf(i,l+2,m)+a42(2)*wf(i,l+1,m) &
                       + a42(3)*wf(i,l  ,m)+a42(4)*wf(i,l-1,m) &
                       + a42(5)*wf(i,l-2,m)+a42(6)*wf(i,l-3,m) &
                       + a42(7)*wf(i,l-4,m)
          dpeta(i,l,m)=  a42(1)*pf(i,l+2,m)+a42(2)*pf(i,l+1,m) &
                       + a42(3)*pf(i,l  ,m)+a42(4)*pf(i,l-1,m) &
                       + a42(5)*pf(i,l-2,m)+a42(6)*pf(i,l-3,m) &
                       + a42(7)*pf(i,l-4,m)
          dreta(i,l,m)=  a42(1)*rf(i,l+2,m)+a42(2)*rf(i,l+1,m) &
                       + a42(3)*rf(i,l  ,m)+a42(4)*rf(i,l-1,m) &
                       + a42(5)*rf(i,l-2,m)+a42(6)*rf(i,l-3,m) &
                       + a42(7)*rf(i,l-4,m)
       enddo
    enddo

    j=ny-1
    l=j-nymngh
    do i=1,ngh
       do m=1,ngh
          dueta(i,l,m)=  a51(1)*uf(i,l+1,m)+a51(2)*uf(i,l  ,m) &
                       + a51(3)*uf(i,l-1,m)+a51(4)*uf(i,l-2,m) &
                       + a51(5)*uf(i,l-3,m)+a51(6)*uf(i,l-4,m) &
                       + a51(7)*uf(i,l-5,m)
          dveta(i,l,m)=  a51(1)*vf(i,l+1,m)+a51(2)*vf(i,l  ,m) &
                       + a51(3)*vf(i,l-1,m)+a51(4)*vf(i,l-2,m) &
                       + a51(5)*vf(i,l-3,m)+a51(6)*vf(i,l-4,m) &
                       + a51(7)*vf(i,l-5,m)
          dweta(i,l,m)=  a51(1)*wf(i,l+1,m)+a51(2)*wf(i,l  ,m) &
                       + a51(3)*wf(i,l-1,m)+a51(4)*wf(i,l-2,m) &
                       + a51(5)*wf(i,l-3,m)+a51(6)*wf(i,l-4,m) &
                       + a51(7)*wf(i,l-5,m)
          dpeta(i,l,m)=  a51(1)*pf(i,l+1,m)+a51(2)*pf(i,l  ,m) &
                       + a51(3)*pf(i,l-1,m)+a51(4)*pf(i,l-2,m) &
                       + a51(5)*pf(i,l-3,m)+a51(6)*pf(i,l-4,m) &
                       + a51(7)*pf(i,l-5,m)
          dreta(i,l,m)=  a51(1)*rf(i,l+1,m)+a51(2)*rf(i,l  ,m) &
                       + a51(3)*rf(i,l-1,m)+a51(4)*rf(i,l-2,m) &
                       + a51(5)*rf(i,l-3,m)+a51(6)*rf(i,l-4,m) &
                       + a51(7)*rf(i,l-5,m)
       enddo
    enddo

    j=ny
    l=j-nymngh
    do i=1,ngh
       do m=1,ngh
          dueta(i,l,m)=  a60(1)*uf(i,l  ,m)+a60(2)*uf(i,l-1,m) &
                       + a60(3)*uf(i,l-2,m)+a60(4)*uf(i,l-3,m) &
                       + a60(5)*uf(i,l-4,m)+a60(6)*uf(i,l-5,m) &
                       + a60(7)*uf(i,l-6,m)
          dveta(i,l,m)=  a60(1)*vf(i,l  ,m)+a60(2)*vf(i,l-1,m) &
                       + a60(3)*vf(i,l-2,m)+a60(4)*vf(i,l-3,m) &
                       + a60(5)*vf(i,l-4,m)+a60(6)*vf(i,l-5,m) &
                       + a60(7)*vf(i,l-6,m)
          dweta(i,l,m)=  a60(1)*wf(i,l  ,m)+a60(2)*wf(i,l-1,m) &
                       + a60(3)*wf(i,l-2,m)+a60(4)*wf(i,l-3,m) &
                       + a60(5)*wf(i,l-4,m)+a60(6)*wf(i,l-5,m) &
                       + a60(7)*wf(i,l-6,m)
          dpeta(i,l,m)=  a60(1)*pf(i,l  ,m)+a60(2)*pf(i,l-1,m) &
                       + a60(3)*pf(i,l-2,m)+a60(4)*pf(i,l-3,m) &
                       + a60(5)*pf(i,l-4,m)+a60(6)*pf(i,l-5,m) &
                       + a60(7)*pf(i,l-6,m)
          dreta(i,l,m)=  a60(1)*rf(i,l  ,m)+a60(2)*rf(i,l-1,m) &
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
       do i=1,ngh
          do l=1,ngh
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
    do i=1,ngh
       do l=1,ngh
          duphi(i,l,m)=  a42(1)*uf(i,l,m+2)+a42(2)*uf(i,l,m+1) &
                       + a42(3)*uf(i,l,m  )+a42(4)*uf(i,l,m-1) &
                       + a42(5)*uf(i,l,m-2)+a42(6)*uf(i,l,m-3) &
                       + a42(7)*uf(i,l,m-4)
          dvphi(i,l,m)=  a42(1)*vf(i,l,m+2)+a42(2)*vf(i,l,m+1) &
                       + a42(3)*vf(i,l,m  )+a42(4)*vf(i,l,m-1) &
                       + a42(5)*vf(i,l,m-2)+a42(6)*vf(i,l,m-3) &
                       + a42(7)*vf(i,l,m-4)
          dwphi(i,l,m)=  a42(1)*wf(i,l,m+2)+a42(2)*wf(i,l,m+1) &
                       + a42(3)*wf(i,l,m  )+a42(4)*wf(i,l,m-1) &
                       + a42(5)*wf(i,l,m-2)+a42(6)*wf(i,l,m-3) &
                       + a42(7)*wf(i,l,m-4)
          dpphi(i,l,m)=  a42(1)*pf(i,l,m+2)+a42(2)*pf(i,l,m+1) &
                       + a42(3)*pf(i,l,m  )+a42(4)*pf(i,l,m-1) &
                       + a42(5)*pf(i,l,m-2)+a42(6)*pf(i,l,m-3) &
                       + a42(7)*pf(i,l,m-4)
          drphi(i,l,m)=  a42(1)*rf(i,l,m+2)+a42(2)*rf(i,l,m+1) &
                       + a42(3)*rf(i,l,m  )+a42(4)*rf(i,l,m-1) &
                       + a42(5)*rf(i,l,m-2)+a42(6)*rf(i,l,m-3) &
                       + a42(7)*rf(i,l,m-4)
       enddo
    enddo

    k=nz-1
    m=k-nzmngh
    do i=1,ngh
       do l=1,ngh
          duphi(i,l,m)=  a51(1)*uf(i,l,m+1)+a51(2)*uf(i,l,m  ) &
                       + a51(3)*uf(i,l,m-1)+a51(4)*uf(i,l,m-2) &
                       + a51(5)*uf(i,l,m-3)+a51(6)*uf(i,l,m-4) &
                       + a51(7)*uf(i,l,m-5)
          dvphi(i,l,m)=  a51(1)*vf(i,l,m+1)+a51(2)*vf(i,l,m  ) &
                       + a51(3)*vf(i,l,m-1)+a51(4)*vf(i,l,m-2) &
                       + a51(5)*vf(i,l,m-3)+a51(6)*vf(i,l,m-4) &
                       + a51(7)*vf(i,l,m-5)
          dwphi(i,l,m)=  a51(1)*wf(i,l,m+1)+a51(2)*wf(i,l,m  ) &
                       + a51(3)*wf(i,l,m-1)+a51(4)*wf(i,l,m-2) &
                       + a51(5)*wf(i,l,m-3)+a51(6)*wf(i,l,m-4) &
                       + a51(7)*wf(i,l,m-5)
          dpphi(i,l,m)=  a51(1)*pf(i,l,m+1)+a51(2)*pf(i,l,m  ) &
                       + a51(3)*pf(i,l,m-1)+a51(4)*pf(i,l,m-2) &
                       + a51(5)*pf(i,l,m-3)+a51(6)*pf(i,l,m-4) &
                       + a51(7)*pf(i,l,m-5)
          drphi(i,l,m)=  a51(1)*rf(i,l,m+1)+a51(2)*rf(i,l,m  ) &
                       + a51(3)*rf(i,l,m-1)+a51(4)*rf(i,l,m-2) &
                       + a51(5)*rf(i,l,m-3)+a51(6)*rf(i,l,m-4) &
                       + a51(7)*rf(i,l,m-5)
       enddo
    enddo

    k=nz
    m=k-nzmngh
    do i=1,ngh
       do l=1,ngh
          duphi(i,l,m)=  a60(1)*uf(i,l,m  )+a60(2)*uf(i,l,m-1) &
                       + a60(3)*uf(i,l,m-2)+a60(4)*uf(i,l,m-3) &
                       + a60(5)*uf(i,l,m-4)+a60(6)*uf(i,l,m-5) &
                       + a60(7)*uf(i,l,m-6)
          dvphi(i,l,m)=  a60(1)*vf(i,l,m  )+a60(2)*vf(i,l,m-1) &
                       + a60(3)*vf(i,l,m-2)+a60(4)*vf(i,l,m-3) &
                       + a60(5)*vf(i,l,m-4)+a60(6)*vf(i,l,m-5) &
                       + a60(7)*vf(i,l,m-6)
          dwphi(i,l,m)=  a60(1)*wf(i,l,m  )+a60(2)*wf(i,l,m-1) &
                       + a60(3)*wf(i,l,m-2)+a60(4)*wf(i,l,m-3) &
                       + a60(5)*wf(i,l,m-4)+a60(6)*wf(i,l,m-5) &
                       + a60(7)*wf(i,l,m-6)
          dpphi(i,l,m)=  a60(1)*pf(i,l,m  )+a60(2)*pf(i,l,m-1) &
                       + a60(3)*pf(i,l,m-2)+a60(4)*pf(i,l,m-3) &
                       + a60(5)*pf(i,l,m-4)+a60(6)*pf(i,l,m-5) &
                       + a60(7)*pf(i,l,m-6)
          drphi(i,l,m)=  a60(1)*rf(i,l,m  )+a60(2)*rf(i,l,m-1) &
                       + a60(3)*rf(i,l,m-2)+a60(4)*rf(i,l,m-3) &
                       + a60(5)*rf(i,l,m-4)+a60(6)*rf(i,l,m-5) &
                       + a60(7)*rf(i,l,m-6)
       enddo
    enddo

    ! Compute rt,ut,vt,wt & pt
    ! ========================
    do i=1,ngh
       do j=nymnghp1,ny
          l=j-nymngh
          do k=nzmnghp1,nz
             m=k-nzmngh
             pt(i,l,m) = vg(i,l,m)*( &
                       ((dpksi(i,l,m)*ksi_x(i,j,k)+dpeta(i,l,m)*eta_x(i,j,k)+dpphi(i,l,m)*phi_x(i,j,k))*sintcosp(i,l,m) &
                       +(dpksi(i,l,m)*ksi_y(i,j,k)+dpeta(i,l,m)*eta_y(i,j,k)+dpphi(i,l,m)*phi_y(i,j,k))*sintsinp(i,l,m) &
                       +(dpksi(i,l,m)*ksi_z(i,j,k)+dpeta(i,l,m)*eta_z(i,j,k)+dpphi(i,l,m)*phi_z(i,j,k))*costeta(i,l,m)  &
                       )*ijacob3(i,j,k) + pf(i,l,m)*ir(i,l,m) )
             ut(i,l,m) = vg(i,l,m)*( &
                       ((duksi(i,l,m)*ksi_x(i,j,k)+dueta(i,l,m)*eta_x(i,j,k)+duphi(i,l,m)*phi_x(i,j,k))*sintcosp(i,l,m) &
                       +(duksi(i,l,m)*ksi_y(i,j,k)+dueta(i,l,m)*eta_y(i,j,k)+duphi(i,l,m)*phi_y(i,j,k))*sintsinp(i,l,m) &
                       +(duksi(i,l,m)*ksi_z(i,j,k)+dueta(i,l,m)*eta_z(i,j,k)+duphi(i,l,m)*phi_z(i,j,k))*costeta(i,l,m)  &
                       )*ijacob3(i,j,k) + uf(i,l,m)*ir(i,l,m) )
             vt(i,l,m) = vg(i,l,m)*( &
                       ((dvksi(i,l,m)*ksi_x(i,j,k)+dveta(i,l,m)*eta_x(i,j,k)+dvphi(i,l,m)*phi_x(i,j,k))*sintcosp(i,l,m) &
                       +(dvksi(i,l,m)*ksi_y(i,j,k)+dveta(i,l,m)*eta_y(i,j,k)+dvphi(i,l,m)*phi_y(i,j,k))*sintsinp(i,l,m) &
                       +(dvksi(i,l,m)*ksi_z(i,j,k)+dveta(i,l,m)*eta_z(i,j,k)+dvphi(i,l,m)*phi_z(i,j,k))*costeta(i,l,m)  &
                       )*ijacob3(i,j,k) + vf(i,l,m)*ir(i,l,m) )
             wt(i,l,m) = vg(i,l,m)*( &
                       ((dwksi(i,l,m)*ksi_x(i,j,k)+dweta(i,l,m)*eta_x(i,j,k)+dwphi(i,l,m)*phi_x(i,j,k))*sintcosp(i,l,m) &
                       +(dwksi(i,l,m)*ksi_y(i,j,k)+dweta(i,l,m)*eta_y(i,j,k)+dwphi(i,l,m)*phi_y(i,j,k))*sintsinp(i,l,m) &
                       +(dwksi(i,l,m)*ksi_z(i,j,k)+dweta(i,l,m)*eta_z(i,j,k)+dwphi(i,l,m)*phi_z(i,j,k))*costeta(i,l,m)  &
                       )*ijacob3(i,j,k) + wf(i,l,m)*ir(i,l,m) )
             rt(i,l,m) = vg(i,l,m)*( &
                       ((drksi(i,l,m)*ksi_x(i,j,k)+dreta(i,l,m)*eta_x(i,j,k)+drphi(i,l,m)*phi_x(i,j,k))*sintcosp(i,l,m) &
                       +(drksi(i,l,m)*ksi_y(i,j,k)+dreta(i,l,m)*eta_y(i,j,k)+drphi(i,l,m)*phi_y(i,j,k))*sintsinp(i,l,m) &
                       +(drksi(i,l,m)*ksi_z(i,j,k)+dreta(i,l,m)*eta_z(i,j,k)+drphi(i,l,m)*phi_z(i,j,k))*costeta(i,l,m)  &
                       )*ijacob3(i,j,k) + rf(i,l,m)*ir(i,l,m) )
          enddo
       enddo
    enddo

    ! Update fluxes at each RK step
    ! =============================
    do k=nzmnghp1,nz
       m=k-nzmngh
       do j=nymnghp1,ny
          l=j-nymngh
          do i=1,ngh
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

  end subroutine bc_TD3d_imin_jmax_kmax_c3

  !===============================================================================
  module subroutine bc_TD3d_imax_jmin_kmin_c3
  !===============================================================================
    !> 3D Tam & Dong's BC at imax-jmin-kmin (corner 2,1,1 /right-bottom-front)
    !> - 3D curvilinear version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l
    real(wp) :: cp,av,c2_
    real(wp), dimension(-2:ngh,1:ngh+3,1:ngh+3) :: rf,uf,vf,wf,pf
    real(wp), dimension(ngh,ngh,ngh) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(ngh,ngh,ngh) :: duksi,dvksi,dwksi,dpksi,drksi
    real(wp), dimension(ngh,ngh,ngh) :: dueta,dveta,dweta,dpeta,dreta
    real(wp), dimension(ngh,ngh,ngh) :: duphi,dvphi,dwphi,dpphi,drphi
    !-------------------------------------------------------------------------

    ! Pointers for spherical coordinates
    ! ==================================
    ir=>BC_corner(2,1,1)%i_r
    cosphi=>BC_corner(2,1,1)%cosp
    sinphi=>BC_corner(2,1,1)%sinp
    costeta=>BC_corner(2,1,1)%cost
    sinteta=>BC_corner(2,1,1)%sint
    costcosp=>BC_corner(2,1,1)%costcosp
    costsinp=>BC_corner(2,1,1)%costsinp
    sintcosp=>BC_corner(2,1,1)%sintcosp
    sintsinp=>BC_corner(2,1,1)%sintsinp

    ! Compute fluctuations
    ! ====================
    do i=nxmngh-2,nx
       l=i-nxmngh
       do j=1,nghp3
          do k=1,nghp3
             rf(l,j,k)=rho_n(i,j,k)-BC_face(2,1)%U0(i,j,k,1)
             uf(l,j,k)=   uu(i,j,k)-BC_face(2,1)%U0(i,j,k,2)
             vf(l,j,k)=   vv(i,j,k)-BC_face(2,1)%U0(i,j,k,3)
             wf(l,j,k)=   ww(i,j,k)-BC_face(2,1)%U0(i,j,k,4)
             pf(l,j,k)=  prs(i,j,k)-BC_face(2,1)%U0(i,j,k,5)
          enddo
       enddo
    enddo

    ! Compute group velocity vg
    ! =========================
    ! vg= u0*sin(teta)*cos(phi)+v0*sin(teta)*sin(phi)+w0*cos(teta)
    !   + sqrt{c0^2-[u0*cos(teta)*cos(phi)+v0*cos(teta)*sin(phi)-w0*sin(teta)]^2-[u0*sin(phi)-v0*cos(phi)]^2}
    do i=nxmnghp1,nx
       l=i-nxmngh
       do j=1,ngh
          do k=1,ngh
             vg(l,j,k)= BC_face(2,1)%U0(i,j,k,2)*sintcosp(l,j,k) +     &
               BC_face(2,1)%U0(i,j,k,3)*sintsinp(l,j,k) +     &
               BC_face(2,1)%U0(i,j,k,4)*costeta(l,j,k)          &
                      + sqrt( BC_face(2,1)%U0(i,j,k,6)-                &
                       ( BC_face(2,1)%U0(i,j,k,2)*costcosp(l,j,k) +    &
                   BC_face(2,1)%U0(i,j,k,3)*costsinp(l,j,k) -    &
                     BC_face(2,1)%U0(i,j,k,4)*sinteta(l,j,k) )**2- &
                  ( BC_face(2,1)%U0(i,j,k,2)*sinphi(l,j,k) -      &
                     BC_face(2,1)%U0(i,j,k,3)*cosphi(l,j,k) )**2 )
          enddo
       enddo
    enddo

    ! Non-centered derivatives in x-direction *sin(teta)*cos(phi)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    do i=nxmnghp1,nx-3
       l=i-nxmngh
       do j=1,ngh
          do k=1,ngh
             duksi(l,j,k)=  a7(1)* (uf(l+1,j,k)-uf(l-1,j,k)) &
                          + a7(2)* (uf(l+2,j,k)-uf(l-2,j,k)) &
                          + a7(3)* (uf(l+3,j,k)-uf(l-3,j,k))
             dvksi(l,j,k)=  a7(1)* (vf(l+1,j,k)-vf(l-1,j,k)) &
                          + a7(2)* (vf(l+2,j,k)-vf(l-2,j,k)) &
                          + a7(3)* (vf(l+3,j,k)-vf(l-3,j,k))
             dwksi(l,j,k)=  a7(1)* (wf(l+1,j,k)-wf(l-1,j,k)) &
                          + a7(2)* (wf(l+2,j,k)-wf(l-2,j,k)) &
                          + a7(3)* (wf(l+3,j,k)-wf(l-3,j,k))
             dpksi(l,j,k)=  a7(1)* (pf(l+1,j,k)-pf(l-1,j,k)) &
                          + a7(2)* (pf(l+2,j,k)-pf(l-2,j,k)) &
                          + a7(3)* (pf(l+3,j,k)-pf(l-3,j,k))
             drksi(l,j,k)=  a7(1)* (rf(l+1,j,k)-rf(l-1,j,k)) &
                          + a7(2)* (rf(l+2,j,k)-rf(l-2,j,k)) &
                          + a7(3)* (rf(l+3,j,k)-rf(l-3,j,k))
          enddo
       enddo
    enddo

    i=nx-2
    l=i-nxmngh
    do j=1,ngh
       do k=1,ngh
          duksi(l,j,k)=  a42(1)*uf(l+2,j,k)+a42(2)*uf(l+1,j,k) &
                       + a42(3)*uf(l  ,j,k)+a42(4)*uf(l-1,j,k) &
                       + a42(5)*uf(l-2,j,k)+a42(6)*uf(l-3,j,k) &
                       + a42(7)*uf(l-4,j,k)
          dvksi(l,j,k)=  a42(1)*vf(l+2,j,k)+a42(2)*vf(l+1,j,k) &
                       + a42(3)*vf(l  ,j,k)+a42(4)*vf(l-1,j,k) &
                       + a42(5)*vf(l-2,j,k)+a42(6)*vf(l-3,j,k) &
                       + a42(7)*vf(l-4,j,k)
          dwksi(l,j,k)=  a42(1)*wf(l+2,j,k)+a42(2)*wf(l+1,j,k) &
                       + a42(3)*wf(l  ,j,k)+a42(4)*wf(l-1,j,k) &
                       + a42(5)*wf(l-2,j,k)+a42(6)*wf(l-3,j,k) &
                       + a42(7)*wf(l-4,j,k)
          dpksi(l,j,k)=  a42(1)*pf(l+2,j,k)+a42(2)*pf(l+1,j,k) &
                       + a42(3)*pf(l  ,j,k)+a42(4)*pf(l-1,j,k) &
                       + a42(5)*pf(l-2,j,k)+a42(6)*pf(l-3,j,k) &
                       + a42(7)*pf(l-4,j,k)
          drksi(l,j,k)=  a42(1)*rf(l+2,j,k)+a42(2)*rf(l+1,j,k) &
                       + a42(3)*rf(l  ,j,k)+a42(4)*rf(l-1,j,k) &
                       + a42(5)*rf(l-2,j,k)+a42(6)*rf(l-3,j,k) &
                       + a42(7)*rf(l-4,j,k)
       enddo
    enddo

    i=nx-1
    l=i-nxmngh
    do j=1,ngh
       do k=1,ngh
          duksi(l,j,k)=  a51(1)*uf(l+1,j,k)+a51(2)*uf(l  ,j,k) &
                       + a51(3)*uf(l-1,j,k)+a51(4)*uf(l-2,j,k) &
                       + a51(5)*uf(l-3,j,k)+a51(6)*uf(l-4,j,k) &
                       + a51(7)*uf(l-5,j,k)
          dvksi(l,j,k)=  a51(1)*vf(l+1,j,k)+a51(2)*vf(l  ,j,k) &
                       + a51(3)*vf(l-1,j,k)+a51(4)*vf(l-2,j,k) &
                       + a51(5)*vf(l-3,j,k)+a51(6)*vf(l-4,j,k) &
                       + a51(7)*vf(l-5,j,k)
          dwksi(l,j,k)=  a51(1)*wf(l+1,j,k)+a51(2)*wf(l  ,j,k) &
                       + a51(3)*wf(l-1,j,k)+a51(4)*wf(l-2,j,k) &
                       + a51(5)*wf(l-3,j,k)+a51(6)*wf(l-4,j,k) &
                       + a51(7)*wf(l-5,j,k)
          dpksi(l,j,k)=  a51(1)*pf(l+1,j,k)+a51(2)*pf(l  ,j,k) &
                       + a51(3)*pf(l-1,j,k)+a51(4)*pf(l-2,j,k) &
                       + a51(5)*pf(l-3,j,k)+a51(6)*pf(l-4,j,k) &
                       + a51(7)*pf(l-5,j,k)
          drksi(l,j,k)=  a51(1)*rf(l+1,j,k)+a51(2)*rf(l  ,j,k) &
                       + a51(3)*rf(l-1,j,k)+a51(4)*rf(l-2,j,k) &
                       + a51(5)*rf(l-3,j,k)+a51(6)*rf(l-4,j,k) &
                       + a51(7)*rf(l-5,j,k)
       enddo
    enddo

    i=nx
    l=i-nxmngh
    do j=1,ngh
       do k=1,ngh
          duksi(l,j,k)=  a60(1)*uf(l  ,j,k)+a60(2)*uf(l-1,j,k) &
                       + a60(3)*uf(l-2,j,k)+a60(4)*uf(l-3,j,k) &
                       + a60(5)*uf(l-4,j,k)+a60(6)*uf(l-5,j,k) &
                       + a60(7)*uf(l-6,j,k)
          dvksi(l,j,k)=  a60(1)*vf(l  ,j,k)+a60(2)*vf(l-1,j,k) &
                       + a60(3)*vf(l-2,j,k)+a60(4)*vf(l-3,j,k) &
                       + a60(5)*vf(l-4,j,k)+a60(6)*vf(l-5,j,k) &
                       + a60(7)*vf(l-6,j,k)
          dwksi(l,j,k)=  a60(1)*wf(l  ,j,k)+a60(2)*wf(l-1,j,k) &
                       + a60(3)*wf(l-2,j,k)+a60(4)*wf(l-3,j,k) &
                       + a60(5)*wf(l-4,j,k)+a60(6)*wf(l-5,j,k) &
                       + a60(7)*wf(l-6,j,k)
          dpksi(l,j,k)=  a60(1)*pf(l  ,j,k)+a60(2)*pf(l-1,j,k) &
                       + a60(3)*pf(l-2,j,k)+a60(4)*pf(l-3,j,k) &
                       + a60(5)*pf(l-4,j,k)+a60(6)*pf(l-5,j,k) &
                       + a60(7)*pf(l-6,j,k)
          drksi(l,j,k)=  a60(1)*rf(l  ,j,k)+a60(2)*rf(l-1,j,k) &
                       + a60(3)*rf(l-2,j,k)+a60(4)*rf(l-3,j,k) &
                       + a60(5)*rf(l-4,j,k)+a60(6)*rf(l-5,j,k) &
                       + a60(7)*rf(l-6,j,k)
       enddo
    enddo

    ! Non-centered derivatives in y-direction *sin(teta)*sin(phi)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    j=1
    do l=1,ngh
       do k=1,ngh
          dueta(l,j,k)= a06(1)*uf(l,1,k)+a06(2)*uf(l,2,k) &
                      + a06(3)*uf(l,3,k)+a06(4)*uf(l,4,k) &
                      + a06(5)*uf(l,5,k)+a06(6)*uf(l,6,k) &
                      + a06(7)*uf(l,7,k)
          dveta(l,j,k)= a06(1)*vf(l,1,k)+a06(2)*vf(l,2,k) &
                      + a06(3)*vf(l,3,k)+a06(4)*vf(l,4,k) &
                      + a06(5)*vf(l,5,k)+a06(6)*vf(l,6,k) &
                      + a06(7)*vf(l,7,k)
          dweta(l,j,k)= a06(1)*wf(l,1,k)+a06(2)*wf(l,2,k) &
                      + a06(3)*wf(l,3,k)+a06(4)*wf(l,4,k) &
                      + a06(5)*wf(l,5,k)+a06(6)*wf(l,6,k) &
                      + a06(7)*wf(l,7,k)
          dpeta(l,j,k)= a06(1)*pf(l,1,k)+a06(2)*pf(l,2,k) &
                      + a06(3)*pf(l,3,k)+a06(4)*pf(l,4,k) &
                      + a06(5)*pf(l,5,k)+a06(6)*pf(l,6,k) &
                      + a06(7)*pf(l,7,k)
          dreta(l,j,k)= a06(1)*rf(l,1,k)+a06(2)*rf(l,2,k) &
                      + a06(3)*rf(l,3,k)+a06(4)*rf(l,4,k) &
                      + a06(5)*rf(l,5,k)+a06(6)*rf(l,6,k) &
                      + a06(7)*rf(l,7,k)
       enddo
    enddo

    j=2
    do l=1,ngh
       do k=1,ngh
          dueta(l,j,k)= a15(1)*uf(l,1,k)+a15(2)*uf(l,2,k) &
                      + a15(3)*uf(l,3,k)+a15(4)*uf(l,4,k) &
                      + a15(5)*uf(l,5,k)+a15(6)*uf(l,6,k) &
                      + a15(7)*uf(l,7,k)
          dveta(l,j,k)= a15(1)*vf(l,1,k)+a15(2)*vf(l,2,k) &
                      + a15(3)*vf(l,3,k)+a15(4)*vf(l,4,k) &
                      + a15(5)*vf(l,5,k)+a15(6)*vf(l,6,k) &
                      + a15(7)*vf(l,7,k)
          dweta(l,j,k)= a15(1)*wf(l,1,k)+a15(2)*wf(l,2,k) &
                      + a15(3)*wf(l,3,k)+a15(4)*wf(l,4,k) &
                      + a15(5)*wf(l,5,k)+a15(6)*wf(l,6,k) &
                      + a15(7)*wf(l,7,k)
          dpeta(l,j,k)= a15(1)*pf(l,1,k)+a15(2)*pf(l,2,k) &
                      + a15(3)*pf(l,3,k)+a15(4)*pf(l,4,k) &
                      + a15(5)*pf(l,5,k)+a15(6)*pf(l,6,k) &
                      + a15(7)*pf(l,7,k)
          dreta(l,j,k)= a15(1)*rf(l,1,k)+a15(2)*rf(l,2,k) &
                      + a15(3)*rf(l,3,k)+a15(4)*rf(l,4,k) &
                      + a15(5)*rf(l,5,k)+a15(6)*rf(l,6,k) &
                      + a15(7)*rf(l,7,k)
       enddo
    enddo

    j=3
    do l=1,ngh
       do k=1,ngh
          dueta(l,j,k)= a24(1)*uf(l,1,k)+a24(2)*uf(l,2,k) &
                      + a24(3)*uf(l,3,k)+a24(4)*uf(l,4,k) &
                      + a24(5)*uf(l,5,k)+a24(6)*uf(l,6,k) &
                      + a24(7)*uf(l,7,k)
          dveta(l,j,k)= a24(1)*vf(l,1,k)+a24(2)*vf(l,2,k) &
                      + a24(3)*vf(l,3,k)+a24(4)*vf(l,4,k) &
                      + a24(5)*vf(l,5,k)+a24(6)*vf(l,6,k) &
                      + a24(7)*vf(l,7,k)
          dweta(l,j,k)= a24(1)*wf(l,1,k)+a24(2)*wf(l,2,k) &
                      + a24(3)*wf(l,3,k)+a24(4)*wf(l,4,k) &
                      + a24(5)*wf(l,5,k)+a24(6)*wf(l,6,k) &
                      + a24(7)*wf(l,7,k)
          dpeta(l,j,k)= a24(1)*pf(l,1,k)+a24(2)*pf(l,2,k) &
                      + a24(3)*pf(l,3,k)+a24(4)*pf(l,4,k) &
                      + a24(5)*pf(l,5,k)+a24(6)*pf(l,6,k) &
                      + a24(7)*pf(l,7,k)
          dreta(l,j,k)= a24(1)*rf(l,1,k)+a24(2)*rf(l,2,k) &
                      + a24(3)*rf(l,3,k)+a24(4)*rf(l,4,k) &
                      + a24(5)*rf(l,5,k)+a24(6)*rf(l,6,k) &
                      + a24(7)*rf(l,7,k)
       enddo
    enddo

    do j=4,ngh
       do l=1,ngh
          do k=1,ngh
             dueta(l,j,k)= a7(1)*(uf(l,j+1,k)-uf(l,j-1,k)) &
                         + a7(2)*(uf(l,j+2,k)-uf(l,j-2,k)) &
                         + a7(3)*(uf(l,j+3,k)-uf(l,j-3,k))
             dveta(l,j,k)= a7(1)*(vf(l,j+1,k)-vf(l,j-1,k)) &
                         + a7(2)*(vf(l,j+2,k)-vf(l,j-2,k)) &
                         + a7(3)*(vf(l,j+3,k)-vf(l,j-3,k))
             dweta(l,j,k)= a7(1)*(wf(l,j+1,k)-wf(l,j-1,k)) &
                         + a7(2)*(wf(l,j+2,k)-wf(l,j-2,k)) &
                         + a7(3)*(wf(l,j+3,k)-wf(l,j-3,k))
             dpeta(l,j,k)= a7(1)*(pf(l,j+1,k)-pf(l,j-1,k)) &
                         + a7(2)*(pf(l,j+2,k)-pf(l,j-2,k)) &
                         + a7(3)*(pf(l,j+3,k)-pf(l,j-3,k))
             dreta(l,j,k)= a7(1)*(rf(l,j+1,k)-rf(l,j-1,k)) &
                         + a7(2)*(rf(l,j+2,k)-rf(l,j-2,k)) &
                         + a7(3)*(rf(l,j+3,k)-rf(l,j-3,k))
          enddo
       enddo
    enddo

    ! Non-centered derivatives in z-direction *cos(teta)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    k=1
    do l=1,ngh
       do j=1,ngh
          duphi(l,j,k)=  a06(1)*uf(l,j,1)+a06(2)*uf(l,j,2) &
                       + a06(3)*uf(l,j,3)+a06(4)*uf(l,j,4) &
                       + a06(5)*uf(l,j,5)+a06(6)*uf(l,j,6) &
                       + a06(7)*uf(l,j,7)
          dvphi(l,j,k)=  a06(1)*vf(l,j,1)+a06(2)*vf(l,j,2) &
                       + a06(3)*vf(l,j,3)+a06(4)*vf(l,j,4) &
                       + a06(5)*vf(l,j,5)+a06(6)*vf(l,j,6) &
                       + a06(7)*vf(l,j,7)
          dwphi(l,j,k)=  a06(1)*wf(l,j,1)+a06(2)*wf(l,j,2) &
                       + a06(3)*wf(l,j,3)+a06(4)*wf(l,j,4) &
                       + a06(5)*wf(l,j,5)+a06(6)*wf(l,j,6) &
                       + a06(7)*wf(l,j,7)
          dpphi(l,j,k)=  a06(1)*pf(l,j,1)+a06(2)*pf(l,j,2) &
                       + a06(3)*pf(l,j,3)+a06(4)*pf(l,j,4) &
                       + a06(5)*pf(l,j,5)+a06(6)*pf(l,j,6) &
                       + a06(7)*pf(l,j,7)
          drphi(l,j,k)=  a06(1)*rf(l,j,1)+a06(2)*rf(l,j,2) &
                       + a06(3)*rf(l,j,3)+a06(4)*rf(l,j,4) &
                       + a06(5)*rf(l,j,5)+a06(6)*rf(l,j,6) &
                       + a06(7)*rf(l,j,7)
       enddo
    enddo

    k=2
    do l=1,ngh
       do j=1,ngh
          duphi(l,j,k)=  a15(1)*uf(l,j,1)+a15(2)*uf(l,j,2) &
                       + a15(3)*uf(l,j,3)+a15(4)*uf(l,j,4) &
                       + a15(5)*uf(l,j,5)+a15(6)*uf(l,j,6) &
                       + a15(7)*uf(l,j,7)
          dvphi(l,j,k)=  a15(1)*vf(l,j,1)+a15(2)*vf(l,j,2) &
                       + a15(3)*vf(l,j,3)+a15(4)*vf(l,j,4) &
                       + a15(5)*vf(l,j,5)+a15(6)*vf(l,j,6) &
                       + a15(7)*vf(l,j,7)
          dwphi(l,j,k)=  a15(1)*wf(l,j,1)+a15(2)*wf(l,j,2) &
                       + a15(3)*wf(l,j,3)+a15(4)*wf(l,j,4) &
                       + a15(5)*wf(l,j,5)+a15(6)*wf(l,j,6) &
                       + a15(7)*wf(l,j,7)
          dpphi(l,j,k)=  a15(1)*pf(l,j,1)+a15(2)*pf(l,j,2) &
                       + a15(3)*pf(l,j,3)+a15(4)*pf(l,j,4) &
                       + a15(5)*pf(l,j,5)+a15(6)*pf(l,j,6) &
                       + a15(7)*pf(l,j,7)
          drphi(l,j,k)=  a15(1)*rf(l,j,1)+a15(2)*rf(l,j,2) &
                       + a15(3)*rf(l,j,3)+a15(4)*rf(l,j,4) &
                       + a15(5)*rf(l,j,5)+a15(6)*rf(l,j,6) &
                       + a15(7)*rf(l,j,7)
       enddo
    enddo

    k=3
    do l=1,ngh
       do j=1,ngh
          duphi(l,j,k)=  a24(1)*uf(l,j,1)+a24(2)*uf(l,j,2) &
                       + a24(3)*uf(l,j,3)+a24(4)*uf(l,j,4) &
                       + a24(5)*uf(l,j,5)+a24(6)*uf(l,j,6) &
                       + a24(7)*uf(l,j,7)
          dvphi(l,j,k)=  a24(1)*vf(l,j,1)+a24(2)*vf(l,j,2) &
                       + a24(3)*vf(l,j,3)+a24(4)*vf(l,j,4) &
                       + a24(5)*vf(l,j,5)+a24(6)*vf(l,j,6) &
                       + a24(7)*vf(l,j,7)
          dwphi(l,j,k)=  a24(1)*wf(l,j,1)+a24(2)*wf(l,j,2) &
                       + a24(3)*wf(l,j,3)+a24(4)*wf(l,j,4) &
                       + a24(5)*wf(l,j,5)+a24(6)*wf(l,j,6) &
                       + a24(7)*wf(l,j,7)
          dpphi(l,j,k)=  a24(1)*pf(l,j,1)+a24(2)*pf(l,j,2) &
                       + a24(3)*pf(l,j,3)+a24(4)*pf(l,j,4) &
                       + a24(5)*pf(l,j,5)+a24(6)*pf(l,j,6) &
                       + a24(7)*pf(l,j,7)
          drphi(l,j,k)=  a24(1)*rf(l,j,1)+a24(2)*rf(l,j,2) &
                       + a24(3)*rf(l,j,3)+a24(4)*rf(l,j,4) &
                       + a24(5)*rf(l,j,5)+a24(6)*rf(l,j,6) &
                       + a24(7)*rf(l,j,7)
       enddo
    enddo

    do k=4,ngh
       do l=1,ngh
          do j=1,ngh
             duphi(l,j,k)= a7(1)*(uf(l,j,k+1) - uf(l,j,k-1)) + &
                           a7(2)*(uf(l,j,k+2) - uf(l,j,k-2)) + &
                           a7(3)*(uf(l,j,k+3) - uf(l,j,k-3))
             dvphi(l,j,k)= a7(1)*(vf(l,j,k+1) - vf(l,j,k-1)) + &
                           a7(2)*(vf(l,j,k+2) - vf(l,j,k-2)) + &
                           a7(3)*(vf(l,j,k+3) - vf(l,j,k-3))
             dwphi(l,j,k)= a7(1)*(wf(l,j,k+1) - wf(l,j,k-1)) + &
                           a7(2)*(wf(l,j,k+2) - wf(l,j,k-2)) + &
                           a7(3)*(wf(l,j,k+3) - wf(l,j,k-3))
             dpphi(l,j,k)= a7(1)*(pf(l,j,k+1) - pf(l,j,k-1)) + &
                           a7(2)*(pf(l,j,k+2) - pf(l,j,k-2)) + &
                           a7(3)*(pf(l,j,k+3) - pf(l,j,k-3))
             drphi(l,j,k)= a7(1)*(rf(l,j,k+1) - rf(l,j,k-1)) + &
                           a7(2)*(rf(l,j,k+2) - rf(l,j,k-2)) + &
                           a7(3)*(rf(l,j,k+3) - rf(l,j,k-3))
          enddo
       enddo
    enddo

    ! Compute rt,ut,vt,wt & pt
    ! ========================
    do i=nxmnghp1,nx
       l=i-nxmngh
       do j=1,ngh
          do k=1,ngh
             pt(l,j,k) = vg(l,j,k)*( &
                       ((dpksi(l,j,k)*ksi_x(i,j,k)+dpeta(l,j,k)*eta_x(i,j,k)+dpphi(l,j,k)*phi_x(i,j,k))*sintcosp(l,j,k) &
                       +(dpksi(l,j,k)*ksi_y(i,j,k)+dpeta(l,j,k)*eta_y(i,j,k)+dpphi(l,j,k)*phi_y(i,j,k))*sintsinp(l,j,k) &
                       +(dpksi(l,j,k)*ksi_z(i,j,k)+dpeta(l,j,k)*eta_z(i,j,k)+dpphi(l,j,k)*phi_z(i,j,k))*costeta(l,j,k)  &
                       )*ijacob3(i,j,k) + pf(l,j,k)*ir(l,j,k) )
             ut(l,j,k) = vg(l,j,k)*( &
                       ((duksi(l,j,k)*ksi_x(i,j,k)+dueta(l,j,k)*eta_x(i,j,k)+duphi(l,j,k)*phi_x(i,j,k))*sintcosp(l,j,k) &
                       +(duksi(l,j,k)*ksi_y(i,j,k)+dueta(l,j,k)*eta_y(i,j,k)+duphi(l,j,k)*phi_y(i,j,k))*sintsinp(l,j,k) &
                       +(duksi(l,j,k)*ksi_z(i,j,k)+dueta(l,j,k)*eta_z(i,j,k)+duphi(l,j,k)*phi_z(i,j,k))*costeta(l,j,k)  &
                       )*ijacob3(i,j,k) + uf(l,j,k)*ir(l,j,k) )
             vt(l,j,k) = vg(l,j,k)*( &
                       ((dvksi(l,j,k)*ksi_x(i,j,k)+dveta(l,j,k)*eta_x(i,j,k)+dvphi(l,j,k)*phi_x(i,j,k))*sintcosp(l,j,k) &
                       +(dvksi(l,j,k)*ksi_y(i,j,k)+dveta(l,j,k)*eta_y(i,j,k)+dvphi(l,j,k)*phi_y(i,j,k))*sintsinp(l,j,k) &
                       +(dvksi(l,j,k)*ksi_z(i,j,k)+dveta(l,j,k)*eta_z(i,j,k)+dvphi(l,j,k)*phi_z(i,j,k))*costeta(l,j,k)  &
                       )*ijacob3(i,j,k) + vf(l,j,k)*ir(l,j,k) )
             wt(l,j,k) = vg(l,j,k)*( &
                       ((dwksi(l,j,k)*ksi_x(i,j,k)+dweta(l,j,k)*eta_x(i,j,k)+dwphi(l,j,k)*phi_x(i,j,k))*sintcosp(l,j,k) &
                       +(dwksi(l,j,k)*ksi_y(i,j,k)+dweta(l,j,k)*eta_y(i,j,k)+dwphi(l,j,k)*phi_y(i,j,k))*sintsinp(l,j,k) &
                       +(dwksi(l,j,k)*ksi_z(i,j,k)+dweta(l,j,k)*eta_z(i,j,k)+dwphi(l,j,k)*phi_z(i,j,k))*costeta(l,j,k)  &
                       )*ijacob3(i,j,k) + wf(l,j,k)*ir(l,j,k) )
             rt(l,j,k) = vg(l,j,k)*( &
                       ((drksi(l,j,k)*ksi_x(i,j,k)+dreta(l,j,k)*eta_x(i,j,k)+drphi(l,j,k)*phi_x(i,j,k))*sintcosp(l,j,k) &
                       +(drksi(l,j,k)*ksi_y(i,j,k)+dreta(l,j,k)*eta_y(i,j,k)+drphi(l,j,k)*phi_y(i,j,k))*sintsinp(l,j,k) &
                       +(drksi(l,j,k)*ksi_z(i,j,k)+dreta(l,j,k)*eta_z(i,j,k)+drphi(l,j,k)*phi_z(i,j,k))*costeta(l,j,k)  &
                       )*ijacob3(i,j,k) + rf(l,j,k)*ir(l,j,k) )
          enddo
       enddo
    enddo

    ! Update fluxes at each RK step
    ! =============================
    do k=1,ngh
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

  end subroutine bc_TD3d_imax_jmin_kmin_c3

  !===============================================================================
  module subroutine bc_TD3d_imax_jmin_kmax_c3
  !===============================================================================
    !> 3D Tam & Dong's BC at imax-jmin-kmax (corner 2,1,2 /right-bottom-back)
    !> - 3D curvilinear version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l,m
    real(wp) :: cp,av,c2_
    real(wp), dimension(-2:ngh,1:ngh+3,-2:ngh) :: rf,uf,vf,wf,pf
    real(wp), dimension(ngh,ngh,ngh) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(ngh,ngh,ngh) :: duksi,dvksi,dwksi,dpksi,drksi
    real(wp), dimension(ngh,ngh,ngh) :: dueta,dveta,dweta,dpeta,dreta
    real(wp), dimension(ngh,ngh,ngh) :: duphi,dvphi,dwphi,dpphi,drphi
    !-------------------------------------------------------------------------

    ! Pointers for spherical coordinates
    ! ==================================
    ir=>BC_corner(2,1,2)%i_r
    cosphi=>BC_corner(2,1,2)%cosp
    sinphi=>BC_corner(2,1,2)%sinp
    costeta=>BC_corner(2,1,2)%cost
    sinteta=>BC_corner(2,1,2)%sint
    costcosp=>BC_corner(2,1,2)%costcosp
    costsinp=>BC_corner(2,1,2)%costsinp
    sintcosp=>BC_corner(2,1,2)%sintcosp
    sintsinp=>BC_corner(2,1,2)%sintsinp

    ! Compute fluctuations
    ! ====================
    do i=nxmngh-2,nx
       l=i-nxmngh
       do j=1,nghp3
          do k=nzmngh-2,nz
             m=k-nzmngh
             rf(l,j,m)=rho_n(i,j,k)-BC_face(2,1)%U0(i,j,k,1)
             uf(l,j,m)=   uu(i,j,k)-BC_face(2,1)%U0(i,j,k,2)
             vf(l,j,m)=   vv(i,j,k)-BC_face(2,1)%U0(i,j,k,3)
             wf(l,j,m)=   ww(i,j,k)-BC_face(2,1)%U0(i,j,k,4)
             pf(l,j,m)=  prs(i,j,k)-BC_face(2,1)%U0(i,j,k,5)
          enddo
       enddo
    enddo

    ! Compute group velocity vg
    ! =========================
    ! vg= u0*sin(teta)*cos(phi)+v0*sin(teta)*sin(phi)+w0*cos(teta)
    !   + sqrt{c0^2-[u0*cos(teta)*cos(phi)+v0*cos(teta)*sin(phi)-w0*sin(teta)]^2-[u0*sin(phi)-v0*cos(phi)]^2}
    do i=nxmnghp1,nx
       l=i-nxmngh
       do j=1,ngh
          do k=nzmnghp1,nz
             m=k-nzmngh
             vg(l,j,m)= BC_face(2,1)%U0(i,j,k,2)*sintcosp(l,j,m) +     &
               BC_face(2,1)%U0(i,j,k,3)*sintsinp(l,j,m) +     &
               BC_face(2,1)%U0(i,j,k,4)*costeta(l,j,m)          &
                      + sqrt( BC_face(2,1)%U0(i,j,k,6)-                &
                       ( BC_face(2,1)%U0(i,j,k,2)*costcosp(l,j,m) +    &
                   BC_face(2,1)%U0(i,j,k,3)*costsinp(l,j,m) -    &
                     BC_face(2,1)%U0(i,j,k,4)*sinteta(l,j,m) )**2- &
                  ( BC_face(2,1)%U0(i,j,k,2)*sinphi(l,j,m) -      &
                     BC_face(2,1)%U0(i,j,k,3)*cosphi(l,j,m) )**2 )
          enddo
       enddo
    enddo

    ! Non-centered derivatives in x-direction *sin(teta)*cos(phi)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    do i=nxmnghp1,nx-3
       l=i-nxmngh
       do j=1,ngh
          do m=1,ngh
             duksi(l,j,m)=  a7(1)* (uf(l+1,j,m)-uf(l-1,j,m)) &
                          + a7(2)* (uf(l+2,j,m)-uf(l-2,j,m)) &
                          + a7(3)* (uf(l+3,j,m)-uf(l-3,j,m))
             dvksi(l,j,m)=  a7(1)* (vf(l+1,j,m)-vf(l-1,j,m)) &
                          + a7(2)* (vf(l+2,j,m)-vf(l-2,j,m)) &
                          + a7(3)* (vf(l+3,j,m)-vf(l-3,j,m))
             dwksi(l,j,m)=  a7(1)* (wf(l+1,j,m)-wf(l-1,j,m)) &
                          + a7(2)* (wf(l+2,j,m)-wf(l-2,j,m)) &
                          + a7(3)* (wf(l+3,j,m)-wf(l-3,j,m))
             dpksi(l,j,m)=  a7(1)* (pf(l+1,j,m)-pf(l-1,j,m)) &
                          + a7(2)* (pf(l+2,j,m)-pf(l-2,j,m)) &
                          + a7(3)* (pf(l+3,j,m)-pf(l-3,j,m))
             drksi(l,j,m)=  a7(1)* (rf(l+1,j,m)-rf(l-1,j,m)) &
                          + a7(2)* (rf(l+2,j,m)-rf(l-2,j,m)) &
                          + a7(3)* (rf(l+3,j,m)-rf(l-3,j,m))
          enddo
       enddo
    enddo

    i=nx-2
    l=i-nxmngh
    do j=1,ngh
       do m=1,ngh
          duksi(l,j,m)=  a42(1)*uf(l+2,j,m)+a42(2)*uf(l+1,j,m) &
                       + a42(3)*uf(l  ,j,m)+a42(4)*uf(l-1,j,m) &
                       + a42(5)*uf(l-2,j,m)+a42(6)*uf(l-3,j,m) &
                       + a42(7)*uf(l-4,j,m)
          dvksi(l,j,m)=  a42(1)*vf(l+2,j,m)+a42(2)*vf(l+1,j,m) &
                       + a42(3)*vf(l  ,j,m)+a42(4)*vf(l-1,j,m) &
                       + a42(5)*vf(l-2,j,m)+a42(6)*vf(l-3,j,m) &
                       + a42(7)*vf(l-4,j,m)
          dwksi(l,j,m)=  a42(1)*wf(l+2,j,m)+a42(2)*wf(l+1,j,m) &
                       + a42(3)*wf(l  ,j,m)+a42(4)*wf(l-1,j,m) &
                       + a42(5)*wf(l-2,j,m)+a42(6)*wf(l-3,j,m) &
                       + a42(7)*wf(l-4,j,m)
          dpksi(l,j,m)=  a42(1)*pf(l+2,j,m)+a42(2)*pf(l+1,j,m) &
                       + a42(3)*pf(l  ,j,m)+a42(4)*pf(l-1,j,m) &
                       + a42(5)*pf(l-2,j,m)+a42(6)*pf(l-3,j,m) &
                       + a42(7)*pf(l-4,j,m)
          drksi(l,j,m)=  a42(1)*rf(l+2,j,m)+a42(2)*rf(l+1,j,m) &
                       + a42(3)*rf(l  ,j,m)+a42(4)*rf(l-1,j,m) &
                       + a42(5)*rf(l-2,j,m)+a42(6)*rf(l-3,j,m) &
                       + a42(7)*rf(l-4,j,m)
       enddo
    enddo

    i=nx-1
    l=i-nxmngh
    do j=1,ngh
       do m=1,ngh
          duksi(l,j,m)=  a51(1)*uf(l+1,j,m)+a51(2)*uf(l  ,j,m) &
                       + a51(3)*uf(l-1,j,m)+a51(4)*uf(l-2,j,m) &
                       + a51(5)*uf(l-3,j,m)+a51(6)*uf(l-4,j,m) &
                       + a51(7)*uf(l-5,j,m)
          dvksi(l,j,m)=  a51(1)*vf(l+1,j,m)+a51(2)*vf(l  ,j,m) &
                       + a51(3)*vf(l-1,j,m)+a51(4)*vf(l-2,j,m) &
                       + a51(5)*vf(l-3,j,m)+a51(6)*vf(l-4,j,m) &
                       + a51(7)*vf(l-5,j,m)
          dwksi(l,j,m)=  a51(1)*wf(l+1,j,m)+a51(2)*wf(l  ,j,m) &
                       + a51(3)*wf(l-1,j,m)+a51(4)*wf(l-2,j,m) &
                       + a51(5)*wf(l-3,j,m)+a51(6)*wf(l-4,j,m) &
                       + a51(7)*wf(l-5,j,m)
          dpksi(l,j,m)=  a51(1)*pf(l+1,j,m)+a51(2)*pf(l  ,j,m) &
                       + a51(3)*pf(l-1,j,m)+a51(4)*pf(l-2,j,m) &
                       + a51(5)*pf(l-3,j,m)+a51(6)*pf(l-4,j,m) &
                       + a51(7)*pf(l-5,j,m)
          drksi(l,j,m)=  a51(1)*rf(l+1,j,m)+a51(2)*rf(l  ,j,m) &
                       + a51(3)*rf(l-1,j,m)+a51(4)*rf(l-2,j,m) &
                       + a51(5)*rf(l-3,j,m)+a51(6)*rf(l-4,j,m) &
                       + a51(7)*rf(l-5,j,m)
       enddo
    enddo

    i=nx
    l=i-nxmngh
    do j=1,ngh
       do m=1,ngh
          duksi(l,j,m)=  a60(1)*uf(l  ,j,m)+a60(2)*uf(l-1,j,m) &
                       + a60(3)*uf(l-2,j,m)+a60(4)*uf(l-3,j,m) &
                       + a60(5)*uf(l-4,j,m)+a60(6)*uf(l-5,j,m) &
                       + a60(7)*uf(l-6,j,m)
          dvksi(l,j,m)=  a60(1)*vf(l  ,j,m)+a60(2)*vf(l-1,j,m) &
                       + a60(3)*vf(l-2,j,m)+a60(4)*vf(l-3,j,m) &
                       + a60(5)*vf(l-4,j,m)+a60(6)*vf(l-5,j,m) &
                       + a60(7)*vf(l-6,j,m)
          dwksi(l,j,m)=  a60(1)*wf(l  ,j,m)+a60(2)*wf(l-1,j,m) &
                       + a60(3)*wf(l-2,j,m)+a60(4)*wf(l-3,j,m) &
                       + a60(5)*wf(l-4,j,m)+a60(6)*wf(l-5,j,m) &
                       + a60(7)*wf(l-6,j,m)
          dpksi(l,j,m)=  a60(1)*pf(l  ,j,m)+a60(2)*pf(l-1,j,m) &
                       + a60(3)*pf(l-2,j,m)+a60(4)*pf(l-3,j,m) &
                       + a60(5)*pf(l-4,j,m)+a60(6)*pf(l-5,j,m) &
                       + a60(7)*pf(l-6,j,m)
          drksi(l,j,m)=  a60(1)*rf(l  ,j,m)+a60(2)*rf(l-1,j,m) &
                       + a60(3)*rf(l-2,j,m)+a60(4)*rf(l-3,j,m) &
                       + a60(5)*rf(l-4,j,m)+a60(6)*rf(l-5,j,m) &
                       + a60(7)*rf(l-6,j,m)
       enddo
    enddo

    ! Non-centered derivatives in y-direction *sin(teta)*sin(phi)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    j=1
    do i=nxmnghp1,nx
       l=i-nxmngh
       do m=1,ngh
          dueta(l,j,m)= a06(1)*uf(l,1,m)+a06(2)*uf(l,2,m) &
                      + a06(3)*uf(l,3,m)+a06(4)*uf(l,4,m) &
                      + a06(5)*uf(l,5,m)+a06(6)*uf(l,6,m) &
                      + a06(7)*uf(l,7,m)
          dveta(l,j,m)= a06(1)*vf(l,1,m)+a06(2)*vf(l,2,m) &
                      + a06(3)*vf(l,3,m)+a06(4)*vf(l,4,m) &
                      + a06(5)*vf(l,5,m)+a06(6)*vf(l,6,m) &
                      + a06(7)*vf(l,7,m)
          dweta(l,j,m)= a06(1)*wf(l,1,m)+a06(2)*wf(l,2,m) &
                      + a06(3)*wf(l,3,m)+a06(4)*wf(l,4,m) &
                      + a06(5)*wf(l,5,m)+a06(6)*wf(l,6,m) &
                      + a06(7)*wf(l,7,m)
          dpeta(l,j,m)= a06(1)*pf(l,1,m)+a06(2)*pf(l,2,m) &
                      + a06(3)*pf(l,3,m)+a06(4)*pf(l,4,m) &
                      + a06(5)*pf(l,5,m)+a06(6)*pf(l,6,m) &
                      + a06(7)*pf(l,7,m)
          dreta(l,j,m)= a06(1)*rf(l,1,m)+a06(2)*rf(l,2,m) &
                      + a06(3)*rf(l,3,m)+a06(4)*rf(l,4,m) &
                      + a06(5)*rf(l,5,m)+a06(6)*rf(l,6,m) &
                      + a06(7)*rf(l,7,m)
       enddo
    enddo

    j=2
    do i=nxmnghp1,nx
       l=i-nxmngh
       do m=1,ngh
          dueta(l,j,m)= a15(1)*uf(l,1,m)+a15(2)*uf(l,2,m) &
                      + a15(3)*uf(l,3,m)+a15(4)*uf(l,4,m) &
                      + a15(5)*uf(l,5,m)+a15(6)*uf(l,6,m) &
                      + a15(7)*uf(l,7,m)
          dveta(l,j,m)= a15(1)*vf(l,1,m)+a15(2)*vf(l,2,m) &
                      + a15(3)*vf(l,3,m)+a15(4)*vf(l,4,m) &
                      + a15(5)*vf(l,5,m)+a15(6)*vf(l,6,m) &
                      + a15(7)*vf(l,7,m)
          dweta(l,j,m)= a15(1)*wf(l,1,m)+a15(2)*wf(l,2,m) &
                      + a15(3)*wf(l,3,m)+a15(4)*wf(l,4,m) &
                      + a15(5)*wf(l,5,m)+a15(6)*wf(l,6,m) &
                      + a15(7)*wf(l,7,m)
          dpeta(l,j,m)= a15(1)*pf(l,1,m)+a15(2)*pf(l,2,m) &
                      + a15(3)*pf(l,3,m)+a15(4)*pf(l,4,m) &
                      + a15(5)*pf(l,5,m)+a15(6)*pf(l,6,m) &
                      + a15(7)*pf(l,7,m)
          dreta(l,j,m)= a15(1)*rf(l,1,m)+a15(2)*rf(l,2,m) &
                      + a15(3)*rf(l,3,m)+a15(4)*rf(l,4,m) &
                      + a15(5)*rf(l,5,m)+a15(6)*rf(l,6,m) &
                      + a15(7)*rf(l,7,m)
       enddo
    enddo

    j=3
    do i=nxmnghp1,nx
       l=i-nxmngh
       do m=1,ngh
          dueta(l,j,m)= a24(1)*uf(l,1,m)+a24(2)*uf(l,2,m) &
                      + a24(3)*uf(l,3,m)+a24(4)*uf(l,4,m) &
                      + a24(5)*uf(l,5,m)+a24(6)*uf(l,6,m) &
                      + a24(7)*uf(l,7,m)
          dveta(l,j,m)= a24(1)*vf(l,1,m)+a24(2)*vf(l,2,m) &
                      + a24(3)*vf(l,3,m)+a24(4)*vf(l,4,m) &
                      + a24(5)*vf(l,5,m)+a24(6)*vf(l,6,m) &
                      + a24(7)*vf(l,7,m)
          dweta(l,j,m)= a24(1)*wf(l,1,m)+a24(2)*wf(l,2,m) &
                      + a24(3)*wf(l,3,m)+a24(4)*wf(l,4,m) &
                      + a24(5)*wf(l,5,m)+a24(6)*wf(l,6,m) &
                      + a24(7)*wf(l,7,m)
          dpeta(l,j,m)= a24(1)*pf(l,1,m)+a24(2)*pf(l,2,m) &
                      + a24(3)*pf(l,3,m)+a24(4)*pf(l,4,m) &
                      + a24(5)*pf(l,5,m)+a24(6)*pf(l,6,m) &
                      + a24(7)*pf(l,7,m)
          dreta(l,j,m)= a24(1)*rf(l,1,m)+a24(2)*rf(l,2,m) &
                      + a24(3)*rf(l,3,m)+a24(4)*rf(l,4,m) &
                      + a24(5)*rf(l,5,m)+a24(6)*rf(l,6,m) &
                      + a24(7)*rf(l,7,m)
       enddo
    enddo

    do j=4,ngh
       do i=nxmnghp1,nx
          l=i-nxmngh
          do m=1,ngh
             dueta(l,j,m)= a7(1)*(uf(l,j+1,m)-uf(l,j-1,m)) &
                         + a7(2)*(uf(l,j+2,m)-uf(l,j-2,m)) &
                         + a7(3)*(uf(l,j+3,m)-uf(l,j-3,m))
             dveta(l,j,m)= a7(1)*(vf(l,j+1,m)-vf(l,j-1,m)) &
                         + a7(2)*(vf(l,j+2,m)-vf(l,j-2,m)) &
                         + a7(3)*(vf(l,j+3,m)-vf(l,j-3,m))
             dweta(l,j,m)= a7(1)*(wf(l,j+1,m)-wf(l,j-1,m)) &
                         + a7(2)*(wf(l,j+2,m)-wf(l,j-2,m)) &
                         + a7(3)*(wf(l,j+3,m)-wf(l,j-3,m))
             dpeta(l,j,m)= a7(1)*(pf(l,j+1,m)-pf(l,j-1,m)) &
                         + a7(2)*(pf(l,j+2,m)-pf(l,j-2,m)) &
                         + a7(3)*(pf(l,j+3,m)-pf(l,j-3,m))
             dreta(l,j,m)= a7(1)*(rf(l,j+1,m)-rf(l,j-1,m)) &
                         + a7(2)*(rf(l,j+2,m)-rf(l,j-2,m)) &
                         + a7(3)*(rf(l,j+3,m)-rf(l,j-3,m))
          enddo
       enddo
    enddo

    ! Non-centered derivatives in z-direction *cos(teta)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    do k=nzmnghp1,nz-3
       m=k-nzmngh
       do l=1,ngh
          do j=1,ngh
             duphi(l,j,m)= a7(1)*(uf(l,j,m+1) - uf(l,j,m-1)) + &
                           a7(2)*(uf(l,j,m+2) - uf(l,j,m-2)) + &
                           a7(3)*(uf(l,j,m+3) - uf(l,j,m-3))
             dvphi(l,j,m)= a7(1)*(vf(l,j,m+1) - vf(l,j,m-1)) + &
                           a7(2)*(vf(l,j,m+2) - vf(l,j,m-2)) + &
                           a7(3)*(vf(l,j,m+3) - vf(l,j,m-3))
             dwphi(l,j,m)= a7(1)*(wf(l,j,m+1) - wf(l,j,m-1)) + &
                           a7(2)*(wf(l,j,m+2) - wf(l,j,m-2)) + &
                           a7(3)*(wf(l,j,m+3) - wf(l,j,m-3))
             dpphi(l,j,m)= a7(1)*(pf(l,j,m+1) - pf(l,j,m-1)) + &
                           a7(2)*(pf(l,j,m+2) - pf(l,j,m-2)) + &
                           a7(3)*(pf(l,j,m+3) - pf(l,j,m-3))
             drphi(l,j,m)= a7(1)*(rf(l,j,m+1) - rf(l,j,m-1)) + &
                           a7(2)*(rf(l,j,m+2) - rf(l,j,m-2)) + &
                           a7(3)*(rf(l,j,m+3) - rf(l,j,m-3))
          enddo
       enddo
    enddo

    k=nz-2
    m=k-nzmngh
    do l=1,ngh
       do j=1,ngh
          duphi(l,j,m)=  a42(1)*uf(l,j,m+2)+a42(2)*uf(l,j,m+1) &
                       + a42(3)*uf(l,j,m  )+a42(4)*uf(l,j,m-1) &
                       + a42(5)*uf(l,j,m-2)+a42(6)*uf(l,j,m-3) &
                       + a42(7)*uf(l,j,m-4)
          dvphi(l,j,m)=  a42(1)*vf(l,j,m+2)+a42(2)*vf(l,j,m+1) &
                       + a42(3)*vf(l,j,m  )+a42(4)*vf(l,j,m-1) &
                       + a42(5)*vf(l,j,m-2)+a42(6)*vf(l,j,m-3) &
                       + a42(7)*vf(l,j,m-4)
          dwphi(l,j,m)=  a42(1)*wf(l,j,m+2)+a42(2)*wf(l,j,m+1) &
                       + a42(3)*wf(l,j,m  )+a42(4)*wf(l,j,m-1) &
                       + a42(5)*wf(l,j,m-2)+a42(6)*wf(l,j,m-3) &
                       + a42(7)*wf(l,j,m-4)
          dpphi(l,j,m)=  a42(1)*pf(l,j,m+2)+a42(2)*pf(l,j,m+1) &
                       + a42(3)*pf(l,j,m  )+a42(4)*pf(l,j,m-1) &
                       + a42(5)*pf(l,j,m-2)+a42(6)*pf(l,j,m-3) &
                       + a42(7)*pf(l,j,m-4)
          drphi(l,j,m)=  a42(1)*rf(l,j,m+2)+a42(2)*rf(l,j,m+1) &
                       + a42(3)*rf(l,j,m  )+a42(4)*rf(l,j,m-1) &
                       + a42(5)*rf(l,j,m-2)+a42(6)*rf(l,j,m-3) &
                       + a42(7)*rf(l,j,m-4)
       enddo
    enddo

    k=nz-1
    m=k-nzmngh
    do l=1,ngh
       do j=1,ngh
          duphi(l,j,m)=  a51(1)*uf(l,j,m+1)+a51(2)*uf(l,j,m  ) &
                       + a51(3)*uf(l,j,m-1)+a51(4)*uf(l,j,m-2) &
                       + a51(5)*uf(l,j,m-3)+a51(6)*uf(l,j,m-4) &
                       + a51(7)*uf(l,j,m-5)
          dvphi(l,j,m)=  a51(1)*vf(l,j,m+1)+a51(2)*vf(l,j,m  ) &
                       + a51(3)*vf(l,j,m-1)+a51(4)*vf(l,j,m-2) &
                       + a51(5)*vf(l,j,m-3)+a51(6)*vf(l,j,m-4) &
                       + a51(7)*vf(l,j,m-5)
          dwphi(l,j,m)=  a51(1)*wf(l,j,m+1)+a51(2)*wf(l,j,m  ) &
                       + a51(3)*wf(l,j,m-1)+a51(4)*wf(l,j,m-2) &
                       + a51(5)*wf(l,j,m-3)+a51(6)*wf(l,j,m-4) &
                       + a51(7)*wf(l,j,m-5)
          dpphi(l,j,m)=  a51(1)*pf(l,j,m+1)+a51(2)*pf(l,j,m  ) &
                       + a51(3)*pf(l,j,m-1)+a51(4)*pf(l,j,m-2) &
                       + a51(5)*pf(l,j,m-3)+a51(6)*pf(l,j,m-4) &
                       + a51(7)*pf(l,j,m-5)
          drphi(l,j,m)=  a51(1)*rf(l,j,m+1)+a51(2)*rf(l,j,m  ) &
                       + a51(3)*rf(l,j,m-1)+a51(4)*rf(l,j,m-2) &
                       + a51(5)*rf(l,j,m-3)+a51(6)*rf(l,j,m-4) &
                       + a51(7)*rf(l,j,m-5)
       enddo
    enddo

    k=nz
    m=k-nzmngh
    do l=1,ngh
       do j=1,ngh
          duphi(l,j,m)=  a60(1)*uf(l,j,m  )+a60(2)*uf(l,j,m-1) &
                       + a60(3)*uf(l,j,m-2)+a60(4)*uf(l,j,m-3) &
                       + a60(5)*uf(l,j,m-4)+a60(6)*uf(l,j,m-5) &
                       + a60(7)*uf(l,j,m-6)
          dvphi(l,j,m)=  a60(1)*vf(l,j,m  )+a60(2)*vf(l,j,m-1) &
                       + a60(3)*vf(l,j,m-2)+a60(4)*vf(l,j,m-3) &
                       + a60(5)*vf(l,j,m-4)+a60(6)*vf(l,j,m-5) &
                       + a60(7)*vf(l,j,m-6)
          dwphi(l,j,m)=  a60(1)*wf(l,j,m  )+a60(2)*wf(l,j,m-1) &
                       + a60(3)*wf(l,j,m-2)+a60(4)*wf(l,j,m-3) &
                       + a60(5)*wf(l,j,m-4)+a60(6)*wf(l,j,m-5) &
                       + a60(7)*wf(l,j,m-6)
          dpphi(l,j,m)=  a60(1)*pf(l,j,m  )+a60(2)*pf(l,j,m-1) &
                       + a60(3)*pf(l,j,m-2)+a60(4)*pf(l,j,m-3) &
                       + a60(5)*pf(l,j,m-4)+a60(6)*pf(l,j,m-5) &
                       + a60(7)*pf(l,j,m-6)
          drphi(l,j,m)=  a60(1)*rf(l,j,m  )+a60(2)*rf(l,j,m-1) &
                       + a60(3)*rf(l,j,m-2)+a60(4)*rf(l,j,m-3) &
                       + a60(5)*rf(l,j,m-4)+a60(6)*rf(l,j,m-5) &
                       + a60(7)*rf(l,j,m-6)
       enddo
    enddo

    ! Compute rt,ut,vt,wt & pt
    ! ========================
    do i=nxmnghp1,nx
       l=i-nxmngh
       do j=1,ngh
          do k=nzmnghp1,nz
             m=k-nzmngh
             pt(l,j,m) = vg(l,j,m)*( &
                       ((dpksi(l,j,m)*ksi_x(i,j,k)+dpeta(l,j,m)*eta_x(i,j,k)+dpphi(l,j,m)*phi_x(i,j,k))*sintcosp(l,j,m) &
                       +(dpksi(l,j,m)*ksi_y(i,j,k)+dpeta(l,j,m)*eta_y(i,j,k)+dpphi(l,j,m)*phi_y(i,j,k))*sintsinp(l,j,m) &
                       +(dpksi(l,j,m)*ksi_z(i,j,k)+dpeta(l,j,m)*eta_z(i,j,k)+dpphi(l,j,m)*phi_z(i,j,k))*costeta(l,j,m)  &
                       )*ijacob3(i,j,k) + pf(l,j,m)*ir(l,j,m) )
             ut(l,j,m) = vg(l,j,m)*( &
                       ((duksi(l,j,m)*ksi_x(i,j,k)+dueta(l,j,m)*eta_x(i,j,k)+duphi(l,j,m)*phi_x(i,j,k))*sintcosp(l,j,m) &
                       +(duksi(l,j,m)*ksi_y(i,j,k)+dueta(l,j,m)*eta_y(i,j,k)+duphi(l,j,m)*phi_y(i,j,k))*sintsinp(l,j,m) &
                       +(duksi(l,j,m)*ksi_z(i,j,k)+dueta(l,j,m)*eta_z(i,j,k)+duphi(l,j,m)*phi_z(i,j,k))*costeta(l,j,m)  &
                       )*ijacob3(i,j,k) + uf(l,j,m)*ir(l,j,m) )
             vt(l,j,m) = vg(l,j,m)*( &
                       ((dvksi(l,j,m)*ksi_x(i,j,k)+dveta(l,j,m)*eta_x(i,j,k)+dvphi(l,j,m)*phi_x(i,j,k))*sintcosp(l,j,m) &
                       +(dvksi(l,j,m)*ksi_y(i,j,k)+dveta(l,j,m)*eta_y(i,j,k)+dvphi(l,j,m)*phi_y(i,j,k))*sintsinp(l,j,m) &
                       +(dvksi(l,j,m)*ksi_z(i,j,k)+dveta(l,j,m)*eta_z(i,j,k)+dvphi(l,j,m)*phi_z(i,j,k))*costeta(l,j,m)  &
                       )*ijacob3(i,j,k) + vf(l,j,m)*ir(l,j,m) )
             wt(l,j,m) = vg(l,j,m)*( &
                       ((dwksi(l,j,m)*ksi_x(i,j,k)+dweta(l,j,m)*eta_x(i,j,k)+dwphi(l,j,m)*phi_x(i,j,k))*sintcosp(l,j,m) &
                       +(dwksi(l,j,m)*ksi_y(i,j,k)+dweta(l,j,m)*eta_y(i,j,k)+dwphi(l,j,m)*phi_y(i,j,k))*sintsinp(l,j,m) &
                       +(dwksi(l,j,m)*ksi_z(i,j,k)+dweta(l,j,m)*eta_z(i,j,k)+dwphi(l,j,m)*phi_z(i,j,k))*costeta(l,j,m)  &
                       )*ijacob3(i,j,k) + wf(l,j,m)*ir(l,j,m) )
             rt(l,j,m) = vg(l,j,m)*( &
                       ((drksi(l,j,m)*ksi_x(i,j,k)+dreta(l,j,m)*eta_x(i,j,k)+drphi(l,j,m)*phi_x(i,j,k))*sintcosp(l,j,m) &
                       +(drksi(l,j,m)*ksi_y(i,j,k)+dreta(l,j,m)*eta_y(i,j,k)+drphi(l,j,m)*phi_y(i,j,k))*sintsinp(l,j,m) &
                       +(drksi(l,j,m)*ksi_z(i,j,k)+dreta(l,j,m)*eta_z(i,j,k)+drphi(l,j,m)*phi_z(i,j,k))*costeta(l,j,m)  &
                       )*ijacob3(i,j,k) + rf(l,j,m)*ir(l,j,m) )
          enddo
       enddo
    enddo

    ! Update fluxes at each RK step
    ! =============================
    do k=nzmnghp1,nz
       m=k-nzmngh
       do i=nxmnghp1,nx
          l=i-nxmngh
          do j=1,ngh
             cp=cpcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
             av=avcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
             c2_=c2calc_tro(Tmp(i,j,k),rho_n(i,j,k))

             Krho(i,j,k)  = rt(l,j,m)
             Krhou(i,j,k) = uu(i,j,k)*rt(l,j,m)+rho_n(i,j,k)*ut(l,j,m)
             Krhov(i,j,k) = vv(i,j,k)*rt(l,j,m)+rho_n(i,j,k)*vt(l,j,m)
             Krhow(i,j,k) = ww(i,j,k)*rt(l,j,m)+rho_n(i,j,k)*wt(l,j,m)
             Krhoe(i,j,k) = cp/av*(pt(l,j,m)/c2_-rt(l,j,m)) &
                  + (rhoe_n(i,j,k)+prs(i,j,k))/rho_n(i,j,k)*rt(l,j,m) &
                  + rho_n(i,j,k)*(uu(i,j,k)*ut(l,j,m)+vv(i,j,k)*vt(l,j,m)+ww(i,j,k)*wt(l,j,m))
          enddo
       enddo
    enddo

  end subroutine bc_TD3d_imax_jmin_kmax_c3

  !===============================================================================
  module subroutine bc_TD3d_imax_jmax_kmin_c3
  !===============================================================================
    !> 3D Tam & Dong's BC at imax-jmax-kmin (corner 2,2,1 /right-top-front)
    !> - 3D curvilinear version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l,m
    real(wp) :: cp,av,c2_
    real(wp), dimension(-2:ngh,-2:ngh,1:ngh+3) :: rf,uf,vf,wf,pf
    real(wp), dimension(ngh,ngh,ngh) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(ngh,ngh,ngh) :: duksi,dvksi,dwksi,dpksi,drksi
    real(wp), dimension(ngh,ngh,ngh) :: dueta,dveta,dweta,dpeta,dreta
    real(wp), dimension(ngh,ngh,ngh) :: duphi,dvphi,dwphi,dpphi,drphi
    !-------------------------------------------------------------------------

    ! Pointers for spherical coordinates
    ! ==================================
    ir=>BC_corner(2,2,1)%i_r
    cosphi=>BC_corner(2,2,1)%cosp
    sinphi=>BC_corner(2,2,1)%sinp
    costeta=>BC_corner(2,2,1)%cost
    sinteta=>BC_corner(2,2,1)%sint
    costcosp=>BC_corner(2,2,1)%costcosp
    costsinp=>BC_corner(2,2,1)%costsinp
    sintcosp=>BC_corner(2,2,1)%sintcosp
    sintsinp=>BC_corner(2,2,1)%sintsinp

    ! Compute fluctuations
    ! ====================
    do i=nxmngh-2,nx
       l=i-nxmngh
       do j=nymngh-2,ny
          m=j-nymngh
          do k=1,nghp3
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
    ! vg= u0*sin(teta)*cos(phi)+v0*sin(teta)*sin(phi)+w0*cos(teta)
    !   + sqrt{c0^2-[u0*cos(teta)*cos(phi)+v0*cos(teta)*sin(phi)-w0*sin(teta)]^2-[u0*sin(phi)-v0*cos(phi)]^2}
    do i=nxmnghp1,nx
       l=i-nxmngh
       do j=nymnghp1,ny
          m=j-nymngh
          do k=1,ngh
             vg(l,m,k)= BC_face(2,2)%U0(i,m,k,2)*sintcosp(l,m,k) +     &
               BC_face(2,2)%U0(i,m,k,3)*sintsinp(l,m,k) +     &
               BC_face(2,2)%U0(i,m,k,4)*costeta(l,m,k)          &
                      + sqrt( BC_face(2,2)%U0(i,m,k,6)-                &
                       ( BC_face(2,2)%U0(i,m,k,2)*costcosp(l,m,k) +    &
                   BC_face(2,2)%U0(i,m,k,3)*costsinp(l,m,k) -    &
                     BC_face(2,2)%U0(i,m,k,4)*sinteta(l,m,k) )**2- &
                  ( BC_face(2,2)%U0(i,m,k,2)*sinphi(l,m,k) -      &
                     BC_face(2,2)%U0(i,m,k,3)*cosphi(l,m,k) )**2 )
          enddo
       enddo
    enddo

    ! Non-centered derivatives in x-direction *sin(teta)*cos(phi)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    do i=nxmnghp1,nx-3
       l=i-nxmngh
       do j=nymnghp1,ny
          m=j-nymngh
          do k=1,ngh
             duksi(l,m,k)=  a7(1)* (uf(l+1,m,k)-uf(l-1,m,k)) &
                          + a7(2)* (uf(l+2,m,k)-uf(l-2,m,k)) &
                          + a7(3)* (uf(l+3,m,k)-uf(l-3,m,k))
             dvksi(l,m,k)=  a7(1)* (vf(l+1,m,k)-vf(l-1,m,k)) &
                          + a7(2)* (vf(l+2,m,k)-vf(l-2,m,k)) &
                          + a7(3)* (vf(l+3,m,k)-vf(l-3,m,k))
             dwksi(l,m,k)=  a7(1)* (wf(l+1,m,k)-wf(l-1,m,k)) &
                          + a7(2)* (wf(l+2,m,k)-wf(l-2,m,k)) &
                          + a7(3)* (wf(l+3,m,k)-wf(l-3,m,k))
             dpksi(l,m,k)=  a7(1)* (pf(l+1,m,k)-pf(l-1,m,k)) &
                          + a7(2)* (pf(l+2,m,k)-pf(l-2,m,k)) &
                          + a7(3)* (pf(l+3,m,k)-pf(l-3,m,k))
             drksi(l,m,k)=  a7(1)* (rf(l+1,m,k)-rf(l-1,m,k)) &
                          + a7(2)* (rf(l+2,m,k)-rf(l-2,m,k)) &
                          + a7(3)* (rf(l+3,m,k)-rf(l-3,m,k))
          enddo
       enddo
    enddo

    i=nx-2
    l=i-nxmngh
    do j=nymnghp1,ny
       m=j-nymngh
       do k=1,ngh
          duksi(l,m,k)=  a42(1)*uf(l+2,m,k)+a42(2)*uf(l+1,m,k) &
                       + a42(3)*uf(l  ,m,k)+a42(4)*uf(l-1,m,k) &
                       + a42(5)*uf(l-2,m,k)+a42(6)*uf(l-3,m,k) &
                       + a42(7)*uf(l-4,m,k)
          dvksi(l,m,k)=  a42(1)*vf(l+2,m,k)+a42(2)*vf(l+1,m,k) &
                       + a42(3)*vf(l  ,m,k)+a42(4)*vf(l-1,m,k) &
                       + a42(5)*vf(l-2,m,k)+a42(6)*vf(l-3,m,k) &
                       + a42(7)*vf(l-4,m,k)
          dwksi(l,m,k)=  a42(1)*wf(l+2,m,k)+a42(2)*wf(l+1,m,k) &
                       + a42(3)*wf(l  ,m,k)+a42(4)*wf(l-1,m,k) &
                       + a42(5)*wf(l-2,m,k)+a42(6)*wf(l-3,m,k) &
                       + a42(7)*wf(l-4,m,k)
          dpksi(l,m,k)=  a42(1)*pf(l+2,m,k)+a42(2)*pf(l+1,m,k) &
                       + a42(3)*pf(l  ,m,k)+a42(4)*pf(l-1,m,k) &
                       + a42(5)*pf(l-2,m,k)+a42(6)*pf(l-3,m,k) &
                       + a42(7)*pf(l-4,m,k)
          drksi(l,m,k)=  a42(1)*rf(l+2,m,k)+a42(2)*rf(l+1,m,k) &
                       + a42(3)*rf(l  ,m,k)+a42(4)*rf(l-1,m,k) &
                       + a42(5)*rf(l-2,m,k)+a42(6)*rf(l-3,m,k) &
                       + a42(7)*rf(l-4,m,k)
       enddo
    enddo

    i=nx-1
    l=i-nxmngh
    do j=nymnghp1,ny
       m=j-nymngh
       do k=1,ngh
          duksi(l,m,k)=  a51(1)*uf(l+1,m,k)+a51(2)*uf(l  ,m,k) &
                       + a51(3)*uf(l-1,m,k)+a51(4)*uf(l-2,m,k) &
                       + a51(5)*uf(l-3,m,k)+a51(6)*uf(l-4,m,k) &
                       + a51(7)*uf(l-5,m,k)
          dvksi(l,m,k)=  a51(1)*vf(l+1,m,k)+a51(2)*vf(l  ,m,k) &
                       + a51(3)*vf(l-1,m,k)+a51(4)*vf(l-2,m,k) &
                       + a51(5)*vf(l-3,m,k)+a51(6)*vf(l-4,m,k) &
                       + a51(7)*vf(l-5,m,k)
          dwksi(l,m,k)=  a51(1)*wf(l+1,m,k)+a51(2)*wf(l  ,m,k) &
                       + a51(3)*wf(l-1,m,k)+a51(4)*wf(l-2,m,k) &
                       + a51(5)*wf(l-3,m,k)+a51(6)*wf(l-4,m,k) &
                       + a51(7)*wf(l-5,m,k)
          dpksi(l,m,k)=  a51(1)*pf(l+1,m,k)+a51(2)*pf(l  ,m,k) &
                       + a51(3)*pf(l-1,m,k)+a51(4)*pf(l-2,m,k) &
                       + a51(5)*pf(l-3,m,k)+a51(6)*pf(l-4,m,k) &
                       + a51(7)*pf(l-5,m,k)
          drksi(l,m,k)=  a51(1)*rf(l+1,m,k)+a51(2)*rf(l  ,m,k) &
                       + a51(3)*rf(l-1,m,k)+a51(4)*rf(l-2,m,k) &
                       + a51(5)*rf(l-3,m,k)+a51(6)*rf(l-4,m,k) &
                       + a51(7)*rf(l-5,m,k)
       enddo
    enddo

    i=nx
    l=i-nxmngh
    do j=nymnghp1,ny
       m=j-nymngh
       do k=1,ngh
          duksi(l,m,k)=  a60(1)*uf(l  ,m,k)+a60(2)*uf(l-1,m,k) &
                       + a60(3)*uf(l-2,m,k)+a60(4)*uf(l-3,m,k) &
                       + a60(5)*uf(l-4,m,k)+a60(6)*uf(l-5,m,k) &
                       + a60(7)*uf(l-6,m,k)
          dvksi(l,m,k)=  a60(1)*vf(l  ,m,k)+a60(2)*vf(l-1,m,k) &
                       + a60(3)*vf(l-2,m,k)+a60(4)*vf(l-3,m,k) &
                       + a60(5)*vf(l-4,m,k)+a60(6)*vf(l-5,m,k) &
                       + a60(7)*vf(l-6,m,k)
          dwksi(l,m,k)=  a60(1)*wf(l  ,m,k)+a60(2)*wf(l-1,m,k) &
                       + a60(3)*wf(l-2,m,k)+a60(4)*wf(l-3,m,k) &
                       + a60(5)*wf(l-4,m,k)+a60(6)*wf(l-5,m,k) &
                       + a60(7)*wf(l-6,m,k)
          dpksi(l,m,k)=  a60(1)*pf(l  ,m,k)+a60(2)*pf(l-1,m,k) &
                       + a60(3)*pf(l-2,m,k)+a60(4)*pf(l-3,m,k) &
                       + a60(5)*pf(l-4,m,k)+a60(6)*pf(l-5,m,k) &
                       + a60(7)*pf(l-6,m,k)
          drksi(l,m,k)=  a60(1)*rf(l  ,m,k)+a60(2)*rf(l-1,m,k) &
                       + a60(3)*rf(l-2,m,k)+a60(4)*rf(l-3,m,k) &
                       + a60(5)*rf(l-4,m,k)+a60(6)*rf(l-5,m,k) &
                       + a60(7)*rf(l-6,m,k)
       enddo
    enddo

    ! Non-centered derivatives in y-direction *sin(teta)*sin(phi)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    do j=nymnghp1,ny-3
       m=j-nymngh
       do l=1,ngh
          do k=1,ngh
             dueta(l,m,k)=  a7(1)*(uf(l,m+1,k)-uf(l,m-1,k)) &
                          + a7(2)*(uf(l,m+2,k)-uf(l,m-2,k)) &
                          + a7(3)*(uf(l,m+3,k)-uf(l,m-3,k))
             dveta(l,m,k)=  a7(1)*(vf(l,m+1,k)-vf(l,m-1,k)) &
                          + a7(2)*(vf(l,m+2,k)-vf(l,m-2,k)) &
                          + a7(3)*(vf(l,m+3,k)-vf(l,m-3,k))
             dweta(l,m,k)=  a7(1)*(wf(l,m+1,k)-wf(l,m-1,k)) &
                          + a7(2)*(wf(l,m+2,k)-wf(l,m-2,k)) &
                          + a7(3)*(wf(l,m+3,k)-wf(l,m-3,k))
             dpeta(l,m,k)=  a7(1)*(pf(l,m+1,k)-pf(l,m-1,k)) &
                          + a7(2)*(pf(l,m+2,k)-pf(l,m-2,k)) &
                          + a7(3)*(pf(l,m+3,k)-pf(l,m-3,k))
             dreta(l,m,k)=  a7(1)*(rf(l,m+1,k)-rf(l,m-1,k)) &
                          + a7(2)*(rf(l,m+2,k)-rf(l,m-2,k)) &
                          + a7(3)*(rf(l,m+3,k)-rf(l,m-3,k))
          enddo
       enddo
    enddo

    j=ny-2
    m=j-nymngh
    do l=1,ngh
       do k=1,ngh
          dueta(l,m,k)=  a42(1)*uf(l,m+2,k)+a42(2)*uf(l,m+1,k) &
                       + a42(3)*uf(l,m  ,k)+a42(4)*uf(l,m-1,k) &
                       + a42(5)*uf(l,m-2,k)+a42(6)*uf(l,m-3,k) &
                       + a42(7)*uf(l,m-4,k)
          dveta(l,m,k)=  a42(1)*vf(l,m+2,k)+a42(2)*vf(l,m+1,k) &
                       + a42(3)*vf(l,m  ,k)+a42(4)*vf(l,m-1,k) &
                       + a42(5)*vf(l,m-2,k)+a42(6)*vf(l,m-3,k) &
                       + a42(7)*vf(l,m-4,k)
          dweta(l,m,k)=  a42(1)*wf(l,m+2,k)+a42(2)*wf(l,m+1,k) &
                       + a42(3)*wf(l,m  ,k)+a42(4)*wf(l,m-1,k) &
                       + a42(5)*wf(l,m-2,k)+a42(6)*wf(l,m-3,k) &
                       + a42(7)*wf(l,m-4,k)
          dpeta(l,m,k)=  a42(1)*pf(l,m+2,k)+a42(2)*pf(l,m+1,k) &
                       + a42(3)*pf(l,m  ,k)+a42(4)*pf(l,m-1,k) &
                       + a42(5)*pf(l,m-2,k)+a42(6)*pf(l,m-3,k) &
                       + a42(7)*pf(l,m-4,k)
          dreta(l,m,k)=  a42(1)*rf(l,m+2,k)+a42(2)*rf(l,m+1,k) &
                       + a42(3)*rf(l,m  ,k)+a42(4)*rf(l,m-1,k) &
                       + a42(5)*rf(l,m-2,k)+a42(6)*rf(l,m-3,k) &
                       + a42(7)*rf(l,m-4,k)
       enddo
    enddo

    j=ny-1
    m=j-nymngh
    do l=1,ngh
       do k=1,ngh
          dueta(l,m,k)=  a51(1)*uf(l,m+1,k)+a51(2)*uf(l,m  ,k) &
                       + a51(3)*uf(l,m-1,k)+a51(4)*uf(l,m-2,k) &
                       + a51(5)*uf(l,m-3,k)+a51(6)*uf(l,m-4,k) &
                       + a51(7)*uf(l,m-5,k)
          dveta(l,m,k)=  a51(1)*vf(l,m+1,k)+a51(2)*vf(l,m  ,k) &
                       + a51(3)*vf(l,m-1,k)+a51(4)*vf(l,m-2,k) &
                       + a51(5)*vf(l,m-3,k)+a51(6)*vf(l,m-4,k) &
                       + a51(7)*vf(l,m-5,k)
          dweta(l,m,k)=  a51(1)*wf(l,m+1,k)+a51(2)*wf(l,m  ,k) &
                       + a51(3)*wf(l,m-1,k)+a51(4)*wf(l,m-2,k) &
                       + a51(5)*wf(l,m-3,k)+a51(6)*wf(l,m-4,k) &
                       + a51(7)*wf(l,m-5,k)
          dpeta(l,m,k)=  a51(1)*pf(l,m+1,k)+a51(2)*pf(l,m  ,k) &
                       + a51(3)*pf(l,m-1,k)+a51(4)*pf(l,m-2,k) &
                       + a51(5)*pf(l,m-3,k)+a51(6)*pf(l,m-4,k) &
                       + a51(7)*pf(l,m-5,k)
          dreta(l,m,k)=  a51(1)*rf(l,m+1,k)+a51(2)*rf(l,m  ,k) &
                       + a51(3)*rf(l,m-1,k)+a51(4)*rf(l,m-2,k) &
                       + a51(5)*rf(l,m-3,k)+a51(6)*rf(l,m-4,k) &
                       + a51(7)*rf(l,m-5,k)
       enddo
    enddo

    j=ny
    m=j-nymngh
    do l=1,ngh
       do k=1,ngh
          dueta(l,m,k)=  a60(1)*uf(l,m  ,k)+a60(2)*uf(l,m-1,k) &
                       + a60(3)*uf(l,m-2,k)+a60(4)*uf(l,m-3,k) &
                       + a60(5)*uf(l,m-4,k)+a60(6)*uf(l,m-5,k) &
                       + a60(7)*uf(l,m-6,k)
          dveta(l,m,k)=  a60(1)*vf(l,m  ,k)+a60(2)*vf(l,m-1,k) &
                       + a60(3)*vf(l,m-2,k)+a60(4)*vf(l,m-3,k) &
                       + a60(5)*vf(l,m-4,k)+a60(6)*vf(l,m-5,k) &
                       + a60(7)*vf(l,m-6,k)
          dweta(l,m,k)=  a60(1)*wf(l,m  ,k)+a60(2)*wf(l,m-1,k) &
                       + a60(3)*wf(l,m-2,k)+a60(4)*wf(l,m-3,k) &
                       + a60(5)*wf(l,m-4,k)+a60(6)*wf(l,m-5,k) &
                       + a60(7)*wf(l,m-6,k)
          dpeta(l,m,k)=  a60(1)*pf(l,m  ,k)+a60(2)*pf(l,m-1,k) &
                       + a60(3)*pf(l,m-2,k)+a60(4)*pf(l,m-3,k) &
                       + a60(5)*pf(l,m-4,k)+a60(6)*pf(l,m-5,k) &
                       + a60(7)*pf(l,m-6,k)
          dreta(l,m,k)=  a60(1)*rf(l,m  ,k)+a60(2)*rf(l,m-1,k) &
                       + a60(3)*rf(l,m-2,k)+a60(4)*rf(l,m-3,k) &
                       + a60(5)*rf(l,m-4,k)+a60(6)*rf(l,m-5,k) &
                       + a60(7)*rf(l,m-6,k)
       enddo
    enddo

    ! Non-centered derivatives in z-direction *cos(teta)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    k=1
    do l=1,ngh
       do m=1,ngh
          duphi(l,m,k)=  a06(1)*uf(l,m,1)+a06(2)*uf(l,m,2) &
                       + a06(3)*uf(l,m,3)+a06(4)*uf(l,m,4) &
                       + a06(5)*uf(l,m,5)+a06(6)*uf(l,m,6) &
                       + a06(7)*uf(l,m,7)
          dvphi(l,m,k)=  a06(1)*vf(l,m,1)+a06(2)*vf(l,m,2) &
                       + a06(3)*vf(l,m,3)+a06(4)*vf(l,m,4) &
                       + a06(5)*vf(l,m,5)+a06(6)*vf(l,m,6) &
                       + a06(7)*vf(l,m,7)
          dwphi(l,m,k)=  a06(1)*wf(l,m,1)+a06(2)*wf(l,m,2) &
                       + a06(3)*wf(l,m,3)+a06(4)*wf(l,m,4) &
                       + a06(5)*wf(l,m,5)+a06(6)*wf(l,m,6) &
                       + a06(7)*wf(l,m,7)
          dpphi(l,m,k)=  a06(1)*pf(l,m,1)+a06(2)*pf(l,m,2) &
                       + a06(3)*pf(l,m,3)+a06(4)*pf(l,m,4) &
                       + a06(5)*pf(l,m,5)+a06(6)*pf(l,m,6) &
                       + a06(7)*pf(l,m,7)
          drphi(l,m,k)=  a06(1)*rf(l,m,1)+a06(2)*rf(l,m,2) &
                       + a06(3)*rf(l,m,3)+a06(4)*rf(l,m,4) &
                       + a06(5)*rf(l,m,5)+a06(6)*rf(l,m,6) &
                       + a06(7)*rf(l,m,7)
       enddo
    enddo

    k=2
    do l=1,ngh
       do m=1,ngh
          duphi(l,m,k)=  a15(1)*uf(l,m,1)+a15(2)*uf(l,m,2) &
                       + a15(3)*uf(l,m,3)+a15(4)*uf(l,m,4) &
                       + a15(5)*uf(l,m,5)+a15(6)*uf(l,m,6) &
                       + a15(7)*uf(l,m,7)
          dvphi(l,m,k)=  a15(1)*vf(l,m,1)+a15(2)*vf(l,m,2) &
                       + a15(3)*vf(l,m,3)+a15(4)*vf(l,m,4) &
                       + a15(5)*vf(l,m,5)+a15(6)*vf(l,m,6) &
                       + a15(7)*vf(l,m,7)
          dwphi(l,m,k)=  a15(1)*wf(l,m,1)+a15(2)*wf(l,m,2) &
                       + a15(3)*wf(l,m,3)+a15(4)*wf(l,m,4) &
                       + a15(5)*wf(l,m,5)+a15(6)*wf(l,m,6) &
                       + a15(7)*wf(l,m,7)
          dpphi(l,m,k)=  a15(1)*pf(l,m,1)+a15(2)*pf(l,m,2) &
                       + a15(3)*pf(l,m,3)+a15(4)*pf(l,m,4) &
                       + a15(5)*pf(l,m,5)+a15(6)*pf(l,m,6) &
                       + a15(7)*pf(l,m,7)
          drphi(l,m,k)=  a15(1)*rf(l,m,1)+a15(2)*rf(l,m,2) &
                       + a15(3)*rf(l,m,3)+a15(4)*rf(l,m,4) &
                       + a15(5)*rf(l,m,5)+a15(6)*rf(l,m,6) &
                       + a15(7)*rf(l,m,7)
       enddo
    enddo

    k=3
    do l=1,ngh
       do m=1,ngh
          duphi(l,m,k)=  a24(1)*uf(l,m,1)+a24(2)*uf(l,m,2) &
                       + a24(3)*uf(l,m,3)+a24(4)*uf(l,m,4) &
                       + a24(5)*uf(l,m,5)+a24(6)*uf(l,m,6) &
                       + a24(7)*uf(l,m,7)
          dvphi(l,m,k)=  a24(1)*vf(l,m,1)+a24(2)*vf(l,m,2) &
                       + a24(3)*vf(l,m,3)+a24(4)*vf(l,m,4) &
                       + a24(5)*vf(l,m,5)+a24(6)*vf(l,m,6) &
                       + a24(7)*vf(l,m,7)
          dwphi(l,m,k)=  a24(1)*wf(l,m,1)+a24(2)*wf(l,m,2) &
                       + a24(3)*wf(l,m,3)+a24(4)*wf(l,m,4) &
                       + a24(5)*wf(l,m,5)+a24(6)*wf(l,m,6) &
                       + a24(7)*wf(l,m,7)
          dpphi(l,m,k)=  a24(1)*pf(l,m,1)+a24(2)*pf(l,m,2) &
                       + a24(3)*pf(l,m,3)+a24(4)*pf(l,m,4) &
                       + a24(5)*pf(l,m,5)+a24(6)*pf(l,m,6) &
                       + a24(7)*pf(l,m,7)
          drphi(l,m,k)=  a24(1)*rf(l,m,1)+a24(2)*rf(l,m,2) &
                       + a24(3)*rf(l,m,3)+a24(4)*rf(l,m,4) &
                       + a24(5)*rf(l,m,5)+a24(6)*rf(l,m,6) &
                       + a24(7)*rf(l,m,7)
       enddo
    enddo

    do k=4,ngh
       do l=1,ngh
          do m=1,ngh
             duphi(l,m,k)= a7(1)*(uf(l,m,k+1) - uf(l,m,k-1)) + &
                           a7(2)*(uf(l,m,k+2) - uf(l,m,k-2)) + &
                           a7(3)*(uf(l,m,k+3) - uf(l,m,k-3))
             dvphi(l,m,k)= a7(1)*(vf(l,m,k+1) - vf(l,m,k-1)) + &
                           a7(2)*(vf(l,m,k+2) - vf(l,m,k-2)) + &
                           a7(3)*(vf(l,m,k+3) - vf(l,m,k-3))
             dwphi(l,m,k)= a7(1)*(wf(l,m,k+1) - wf(l,m,k-1)) + &
                           a7(2)*(wf(l,m,k+2) - wf(l,m,k-2)) + &
                           a7(3)*(wf(l,m,k+3) - wf(l,m,k-3))
             dpphi(l,m,k)= a7(1)*(pf(l,m,k+1) - pf(l,m,k-1)) + &
                           a7(2)*(pf(l,m,k+2) - pf(l,m,k-2)) + &
                           a7(3)*(pf(l,m,k+3) - pf(l,m,k-3))
             drphi(l,m,k)= a7(1)*(rf(l,m,k+1) - rf(l,m,k-1)) + &
                           a7(2)*(rf(l,m,k+2) - rf(l,m,k-2)) + &
                           a7(3)*(rf(l,m,k+3) - rf(l,m,k-3))
          enddo
       enddo
    enddo

    ! Compute rt,ut,vt,wt & pt
    ! ========================
    do i=nxmnghp1,nx
       l=i-nxmngh
       do j=nymnghp1,ny
          m=j-nymngh
          do k=1,ngh
             pt(l,m,k) = vg(l,m,k)*( &
                       ((dpksi(l,m,k)*ksi_x(i,j,k)+dpeta(l,m,k)*eta_x(i,j,k)+dpphi(l,m,k)*phi_x(i,j,k))*sintcosp(l,m,k) &
                       +(dpksi(l,m,k)*ksi_y(i,j,k)+dpeta(l,m,k)*eta_y(i,j,k)+dpphi(l,m,k)*phi_y(i,j,k))*sintsinp(l,m,k) &
                       +(dpksi(l,m,k)*ksi_z(i,j,k)+dpeta(l,m,k)*eta_z(i,j,k)+dpphi(l,m,k)*phi_z(i,j,k))*costeta(l,m,k)  &
                       )*ijacob3(i,j,k) + pf(l,m,k)*ir(l,m,k) )
             ut(l,m,k) = vg(l,m,k)*( &
                       ((duksi(l,m,k)*ksi_x(i,j,k)+dueta(l,m,k)*eta_x(i,j,k)+duphi(l,m,k)*phi_x(i,j,k))*sintcosp(l,m,k) &
                       +(duksi(l,m,k)*ksi_y(i,j,k)+dueta(l,m,k)*eta_y(i,j,k)+duphi(l,m,k)*phi_y(i,j,k))*sintsinp(l,m,k) &
                       +(duksi(l,m,k)*ksi_z(i,j,k)+dueta(l,m,k)*eta_z(i,j,k)+duphi(l,m,k)*phi_z(i,j,k))*costeta(l,m,k)  &
                       )*ijacob3(i,j,k) + uf(l,m,k)*ir(l,m,k) )
             vt(l,m,k) = vg(l,m,k)*( &
                       ((dvksi(l,m,k)*ksi_x(i,j,k)+dveta(l,m,k)*eta_x(i,j,k)+dvphi(l,m,k)*phi_x(i,j,k))*sintcosp(l,m,k) &
                       +(dvksi(l,m,k)*ksi_y(i,j,k)+dveta(l,m,k)*eta_y(i,j,k)+dvphi(l,m,k)*phi_y(i,j,k))*sintsinp(l,m,k) &
                       +(dvksi(l,m,k)*ksi_z(i,j,k)+dveta(l,m,k)*eta_z(i,j,k)+dvphi(l,m,k)*phi_z(i,j,k))*costeta(l,m,k)  &
                       )*ijacob3(i,j,k) + vf(l,m,k)*ir(l,m,k) )
             wt(l,m,k) = vg(l,m,k)*( &
                       ((dwksi(l,m,k)*ksi_x(i,j,k)+dweta(l,m,k)*eta_x(i,j,k)+dwphi(l,m,k)*phi_x(i,j,k))*sintcosp(l,m,k) &
                       +(dwksi(l,m,k)*ksi_y(i,j,k)+dweta(l,m,k)*eta_y(i,j,k)+dwphi(l,m,k)*phi_y(i,j,k))*sintsinp(l,m,k) &
                       +(dwksi(l,m,k)*ksi_z(i,j,k)+dweta(l,m,k)*eta_z(i,j,k)+dwphi(l,m,k)*phi_z(i,j,k))*costeta(l,m,k)  &
                       )*ijacob3(i,j,k) + wf(l,m,k)*ir(l,m,k) )
             rt(l,m,k) = vg(l,m,k)*( &
                       ((drksi(l,m,k)*ksi_x(i,j,k)+dreta(l,m,k)*eta_x(i,j,k)+drphi(l,m,k)*phi_x(i,j,k))*sintcosp(l,m,k) &
                       +(drksi(l,m,k)*ksi_y(i,j,k)+dreta(l,m,k)*eta_y(i,j,k)+drphi(l,m,k)*phi_y(i,j,k))*sintsinp(l,m,k) &
                       +(drksi(l,m,k)*ksi_z(i,j,k)+dreta(l,m,k)*eta_z(i,j,k)+drphi(l,m,k)*phi_z(i,j,k))*costeta(l,m,k)  &
                       )*ijacob3(i,j,k) + rf(l,m,k)*ir(l,m,k) )
          enddo
       enddo
    enddo

    ! Update fluxes at each RK step
    ! =============================
    do k=1,ngh
       do i=nxmnghp1,nx
          l=i-nxmngh
          do j=nymnghp1,ny
             m=j-nymngh
             cp=cpcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
             av=avcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
             c2_=c2calc_tro(Tmp(i,j,k),rho_n(i,j,k))

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

  end subroutine bc_TD3d_imax_jmax_kmin_c3

  !===============================================================================
  module subroutine bc_TD3d_imax_jmax_kmax_c3
  !===============================================================================
    !> 3D Tam & Dong's BC at imax-jmax-kmax (corner 2,2,2 /right-top-back)
    !> - 3D curvilinear version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l,m,n
    real(wp) :: cp,av,c2_
    real(wp), dimension(-2:ngh,-2:ngh,-2:ngh) :: rf,uf,vf,wf,pf
    real(wp), dimension(ngh,ngh,ngh) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(ngh,ngh,ngh) :: duksi,dvksi,dwksi,dpksi,drksi
    real(wp), dimension(ngh,ngh,ngh) :: dueta,dveta,dweta,dpeta,dreta
    real(wp), dimension(ngh,ngh,ngh) :: duphi,dvphi,dwphi,dpphi,drphi
    !-------------------------------------------------------------------------

    ! Pointers for spherical coordinates
    ! ==================================
    ir=>BC_corner(2,2,2)%i_r
    cosphi=>BC_corner(2,2,2)%cosp
    sinphi=>BC_corner(2,2,2)%sinp
    costeta=>BC_corner(2,2,2)%cost
    sinteta=>BC_corner(2,2,2)%sint
    costcosp=>BC_corner(2,2,2)%costcosp
    costsinp=>BC_corner(2,2,2)%costsinp
    sintcosp=>BC_corner(2,2,2)%sintcosp
    sintsinp=>BC_corner(2,2,2)%sintsinp

    ! Compute fluctuations
    ! ====================
    do i=nxmngh-2,nx
       l=i-nxmngh
       do j=nymngh-2,ny
          m=j-nymngh
          do k=nzmngh-2,nz
             n=k-nzmngh
             rf(l,m,n)=rho_n(i,j,k)-BC_face(2,2)%U0(i,m,k,1)
             uf(l,m,n)=   uu(i,j,k)-BC_face(2,2)%U0(i,m,k,2)
             vf(l,m,n)=   vv(i,j,k)-BC_face(2,2)%U0(i,m,k,3)
             wf(l,m,n)=   ww(i,j,k)-BC_face(2,2)%U0(i,m,k,4)
             pf(l,m,n)=  prs(i,j,k)-BC_face(2,2)%U0(i,m,k,5)
          enddo
       enddo
    enddo

    ! Compute group velocity vg
    ! =========================
    ! vg= u0*sin(teta)*cos(phi)+v0*sin(teta)*sin(phi)+w0*cos(teta)
    !   + sqrt{c0^2-[u0*cos(teta)*cos(phi)+v0*cos(teta)*sin(phi)-w0*sin(teta)]^2-[u0*sin(phi)-v0*cos(phi)]^2}
    do i=nxmnghp1,nx
       l=i-nxmngh
       do m=1,ngh
          do k=nzmnghp1,nz
             n=k-nzmngh
             vg(l,m,n)= BC_face(2,2)%U0(i,m,k,2)*sintcosp(l,m,n) +     &
               BC_face(2,2)%U0(i,m,k,3)*sintsinp(l,m,n) +     &
               BC_face(2,2)%U0(i,m,k,4)*costeta(l,m,n)          &
                      + sqrt( BC_face(2,2)%U0(i,m,n,6)-                &
                       ( BC_face(2,2)%U0(i,m,k,2)*costcosp(l,m,n) +    &
                   BC_face(2,2)%U0(i,m,k,3)*costsinp(l,m,n) -    &
                     BC_face(2,2)%U0(i,m,k,4)*sinteta(l,m,n) )**2- &
                  ( BC_face(2,2)%U0(i,m,k,2)*sinphi(l,m,n) -      &
                     BC_face(2,2)%U0(i,m,k,3)*cosphi(l,m,n) )**2 )
          enddo
       enddo
    enddo

    ! Non-centered derivatives in x-direction *sin(teta)*cos(phi)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    do i=nxmnghp1,nx-3
       l=i-nxmngh
       do m=1,ngh
          do n=1,ngh
             duksi(l,m,n)=  a7(1)* (uf(l+1,m,n)-uf(l-1,m,n)) &
                          + a7(2)* (uf(l+2,m,n)-uf(l-2,m,n)) &
                          + a7(3)* (uf(l+3,m,n)-uf(l-3,m,n))
             dvksi(l,m,n)=  a7(1)* (vf(l+1,m,n)-vf(l-1,m,n)) &
                          + a7(2)* (vf(l+2,m,n)-vf(l-2,m,n)) &
                          + a7(3)* (vf(l+3,m,n)-vf(l-3,m,n))
             dwksi(l,m,n)=  a7(1)* (wf(l+1,m,n)-wf(l-1,m,n)) &
                          + a7(2)* (wf(l+2,m,n)-wf(l-2,m,n)) &
                          + a7(3)* (wf(l+3,m,n)-wf(l-3,m,n))
             dpksi(l,m,n)=  a7(1)* (pf(l+1,m,n)-pf(l-1,m,n)) &
                          + a7(2)* (pf(l+2,m,n)-pf(l-2,m,n)) &
                          + a7(3)* (pf(l+3,m,n)-pf(l-3,m,n))
             drksi(l,m,n)=  a7(1)* (rf(l+1,m,n)-rf(l-1,m,n)) &
                          + a7(2)* (rf(l+2,m,n)-rf(l-2,m,n)) &
                          + a7(3)* (rf(l+3,m,n)-rf(l-3,m,n))
          enddo
       enddo
    enddo

    i=nx-2
    l=i-nxmngh
    do m=1,ngh
       do n=1,ngh
          duksi(l,m,n)=  a42(1)*uf(l+2,m,n)+a42(2)*uf(l+1,m,n) &
                       + a42(3)*uf(l  ,m,n)+a42(4)*uf(l-1,m,n) &
                       + a42(5)*uf(l-2,m,n)+a42(6)*uf(l-3,m,n) &
                       + a42(7)*uf(l-4,m,n)
          dvksi(l,m,n)=  a42(1)*vf(l+2,m,n)+a42(2)*vf(l+1,m,n) &
                       + a42(3)*vf(l  ,m,n)+a42(4)*vf(l-1,m,n) &
                       + a42(5)*vf(l-2,m,n)+a42(6)*vf(l-3,m,n) &
                       + a42(7)*vf(l-4,m,n)
          dwksi(l,m,n)=  a42(1)*wf(l+2,m,n)+a42(2)*wf(l+1,m,n) &
                       + a42(3)*wf(l  ,m,n)+a42(4)*wf(l-1,m,n) &
                       + a42(5)*wf(l-2,m,n)+a42(6)*wf(l-3,m,n) &
                       + a42(7)*wf(l-4,m,n)
          dpksi(l,m,n)=  a42(1)*pf(l+2,m,n)+a42(2)*pf(l+1,m,n) &
                       + a42(3)*pf(l  ,m,n)+a42(4)*pf(l-1,m,n) &
                       + a42(5)*pf(l-2,m,n)+a42(6)*pf(l-3,m,n) &
                       + a42(7)*pf(l-4,m,n)
          drksi(l,m,n)=  a42(1)*rf(l+2,m,n)+a42(2)*rf(l+1,m,n) &
                       + a42(3)*rf(l  ,m,n)+a42(4)*rf(l-1,m,n) &
                       + a42(5)*rf(l-2,m,n)+a42(6)*rf(l-3,m,n) &
                       + a42(7)*rf(l-4,m,n)
       enddo
    enddo

    i=nx-1
    l=i-nxmngh
    do m=1,ngh
       do n=1,ngh
          duksi(l,m,n)=  a51(1)*uf(l+1,m,n)+a51(2)*uf(l  ,m,n) &
                       + a51(3)*uf(l-1,m,n)+a51(4)*uf(l-2,m,n) &
                       + a51(5)*uf(l-3,m,n)+a51(6)*uf(l-4,m,n) &
                       + a51(7)*uf(l-5,m,n)
          dvksi(l,m,n)=  a51(1)*vf(l+1,m,n)+a51(2)*vf(l  ,m,n) &
                       + a51(3)*vf(l-1,m,n)+a51(4)*vf(l-2,m,n) &
                       + a51(5)*vf(l-3,m,n)+a51(6)*vf(l-4,m,n) &
                       + a51(7)*vf(l-5,m,n)
          dwksi(l,m,n)=  a51(1)*wf(l+1,m,n)+a51(2)*wf(l  ,m,n) &
                       + a51(3)*wf(l-1,m,n)+a51(4)*wf(l-2,m,n) &
                       + a51(5)*wf(l-3,m,n)+a51(6)*wf(l-4,m,n) &
                       + a51(7)*wf(l-5,m,n)
          dpksi(l,m,n)=  a51(1)*pf(l+1,m,n)+a51(2)*pf(l  ,m,n) &
                       + a51(3)*pf(l-1,m,n)+a51(4)*pf(l-2,m,n) &
                       + a51(5)*pf(l-3,m,n)+a51(6)*pf(l-4,m,n) &
                       + a51(7)*pf(l-5,m,n)
          drksi(l,m,n)=  a51(1)*rf(l+1,m,n)+a51(2)*rf(l  ,m,n) &
                       + a51(3)*rf(l-1,m,n)+a51(4)*rf(l-2,m,n) &
                       + a51(5)*rf(l-3,m,n)+a51(6)*rf(l-4,m,n) &
                       + a51(7)*rf(l-5,m,n)
       enddo
    enddo

    i=nx
    l=i-nxmngh
    do m=1,ngh
       do n=1,ngh
          duksi(l,m,n)=  a60(1)*uf(l  ,m,n)+a60(2)*uf(l-1,m,n) &
                       + a60(3)*uf(l-2,m,n)+a60(4)*uf(l-3,m,n) &
                       + a60(5)*uf(l-4,m,n)+a60(6)*uf(l-5,m,n) &
                       + a60(7)*uf(l-6,m,n)
          dvksi(l,m,n)=  a60(1)*vf(l  ,m,n)+a60(2)*vf(l-1,m,n) &
                       + a60(3)*vf(l-2,m,n)+a60(4)*vf(l-3,m,n) &
                       + a60(5)*vf(l-4,m,n)+a60(6)*vf(l-5,m,n) &
                       + a60(7)*vf(l-6,m,n)
          dwksi(l,m,n)=  a60(1)*wf(l  ,m,n)+a60(2)*wf(l-1,m,n) &
                       + a60(3)*wf(l-2,m,n)+a60(4)*wf(l-3,m,n) &
                       + a60(5)*wf(l-4,m,n)+a60(6)*wf(l-5,m,n) &
                       + a60(7)*wf(l-6,m,n)
          dpksi(l,m,n)=  a60(1)*pf(l  ,m,n)+a60(2)*pf(l-1,m,n) &
                       + a60(3)*pf(l-2,m,n)+a60(4)*pf(l-3,m,n) &
                       + a60(5)*pf(l-4,m,n)+a60(6)*pf(l-5,m,n) &
                       + a60(7)*pf(l-6,m,n)
          drksi(l,m,n)=  a60(1)*rf(l  ,m,n)+a60(2)*rf(l-1,m,n) &
                       + a60(3)*rf(l-2,m,n)+a60(4)*rf(l-3,m,n) &
                       + a60(5)*rf(l-4,m,n)+a60(6)*rf(l-5,m,n) &
                       + a60(7)*rf(l-6,m,n)
       enddo
    enddo

    ! Non-centered derivatives in y-direction *sin(teta)*sin(phi)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    do j=nymnghp1,ny-3
       m=j-nymngh
       do l=1,ngh
          do n=1,ngh
             dueta(l,m,n)=  a7(1)*(uf(l,m+1,n)-uf(l,m-1,n)) &
                          + a7(2)*(uf(l,m+2,n)-uf(l,m-2,n)) &
                          + a7(3)*(uf(l,m+3,n)-uf(l,m-3,n))
             dveta(l,m,n)=  a7(1)*(vf(l,m+1,n)-vf(l,m-1,n)) &
                          + a7(2)*(vf(l,m+2,n)-vf(l,m-2,n)) &
                          + a7(3)*(vf(l,m+3,n)-vf(l,m-3,n))
             dweta(l,m,n)=  a7(1)*(wf(l,m+1,n)-wf(l,m-1,n)) &
                          + a7(2)*(wf(l,m+2,n)-wf(l,m-2,n)) &
                          + a7(3)*(wf(l,m+3,n)-wf(l,m-3,n))
             dpeta(l,m,n)=  a7(1)*(pf(l,m+1,n)-pf(l,m-1,n)) &
                          + a7(2)*(pf(l,m+2,n)-pf(l,m-2,n)) &
                          + a7(3)*(pf(l,m+3,n)-pf(l,m-3,n))
             dreta(l,m,n)=  a7(1)*(rf(l,m+1,n)-rf(l,m-1,n)) &
                          + a7(2)*(rf(l,m+2,n)-rf(l,m-2,n)) &
                          + a7(3)*(rf(l,m+3,n)-rf(l,m-3,n))
          enddo
       enddo
    enddo

    j=ny-2
    m=j-nymngh
    do l=1,ngh
       do n=1,ngh
          dueta(l,m,n)=  a42(1)*uf(l,m+2,n)+a42(2)*uf(l,m+1,n) &
                       + a42(3)*uf(l,m  ,n)+a42(4)*uf(l,m-1,n) &
                       + a42(5)*uf(l,m-2,n)+a42(6)*uf(l,m-3,n) &
                       + a42(7)*uf(l,m-4,n)
          dveta(l,m,n)=  a42(1)*vf(l,m+2,n)+a42(2)*vf(l,m+1,n) &
                       + a42(3)*vf(l,m  ,n)+a42(4)*vf(l,m-1,n) &
                       + a42(5)*vf(l,m-2,n)+a42(6)*vf(l,m-3,n) &
                       + a42(7)*vf(l,m-4,n)
          dweta(l,m,n)=  a42(1)*wf(l,m+2,n)+a42(2)*wf(l,m+1,n) &
                       + a42(3)*wf(l,m  ,n)+a42(4)*wf(l,m-1,n) &
                       + a42(5)*wf(l,m-2,n)+a42(6)*wf(l,m-3,n) &
                       + a42(7)*wf(l,m-4,n)
          dpeta(l,m,n)=  a42(1)*pf(l,m+2,n)+a42(2)*pf(l,m+1,n) &
                       + a42(3)*pf(l,m  ,n)+a42(4)*pf(l,m-1,n) &
                       + a42(5)*pf(l,m-2,n)+a42(6)*pf(l,m-3,n) &
                       + a42(7)*pf(l,m-4,n)
          dreta(l,m,n)=  a42(1)*rf(l,m+2,n)+a42(2)*rf(l,m+1,n) &
                       + a42(3)*rf(l,m  ,n)+a42(4)*rf(l,m-1,n) &
                       + a42(5)*rf(l,m-2,n)+a42(6)*rf(l,m-3,n) &
                       + a42(7)*rf(l,m-4,n)
       enddo
    enddo

    j=ny-1
    m=j-nymngh
    do l=1,ngh
       do n=1,ngh
          dueta(l,m,n)=  a51(1)*uf(l,m+1,n)+a51(2)*uf(l,m  ,n) &
                       + a51(3)*uf(l,m-1,n)+a51(4)*uf(l,m-2,n) &
                       + a51(5)*uf(l,m-3,n)+a51(6)*uf(l,m-4,n) &
                       + a51(7)*uf(l,m-5,n)
          dveta(l,m,n)=  a51(1)*vf(l,m+1,n)+a51(2)*vf(l,m  ,n) &
                       + a51(3)*vf(l,m-1,n)+a51(4)*vf(l,m-2,n) &
                       + a51(5)*vf(l,m-3,n)+a51(6)*vf(l,m-4,n) &
                       + a51(7)*vf(l,m-5,n)
          dweta(l,m,n)=  a51(1)*wf(l,m+1,n)+a51(2)*wf(l,m  ,n) &
                       + a51(3)*wf(l,m-1,n)+a51(4)*wf(l,m-2,n) &
                       + a51(5)*wf(l,m-3,n)+a51(6)*wf(l,m-4,n) &
                       + a51(7)*wf(l,m-5,n)
          dpeta(l,m,n)=  a51(1)*pf(l,m+1,n)+a51(2)*pf(l,m  ,n) &
                       + a51(3)*pf(l,m-1,n)+a51(4)*pf(l,m-2,n) &
                       + a51(5)*pf(l,m-3,n)+a51(6)*pf(l,m-4,n) &
                       + a51(7)*pf(l,m-5,n)
          dreta(l,m,n)=  a51(1)*rf(l,m+1,n)+a51(2)*rf(l,m  ,n) &
                       + a51(3)*rf(l,m-1,n)+a51(4)*rf(l,m-2,n) &
                       + a51(5)*rf(l,m-3,n)+a51(6)*rf(l,m-4,n) &
                       + a51(7)*rf(l,m-5,n)
       enddo
    enddo

    j=ny
    m=j-nymngh
    do l=1,ngh
       do n=1,ngh
          dueta(l,m,n)=  a60(1)*uf(l,m  ,n)+a60(2)*uf(l,m-1,n) &
                       + a60(3)*uf(l,m-2,n)+a60(4)*uf(l,m-3,n) &
                       + a60(5)*uf(l,m-4,n)+a60(6)*uf(l,m-5,n) &
                       + a60(7)*uf(l,m-6,n)
          dveta(l,m,n)=  a60(1)*vf(l,m  ,n)+a60(2)*vf(l,m-1,n) &
                       + a60(3)*vf(l,m-2,n)+a60(4)*vf(l,m-3,n) &
                       + a60(5)*vf(l,m-4,n)+a60(6)*vf(l,m-5,n) &
                       + a60(7)*vf(l,m-6,n)
          dweta(l,m,n)=  a60(1)*wf(l,m  ,n)+a60(2)*wf(l,m-1,n) &
                       + a60(3)*wf(l,m-2,n)+a60(4)*wf(l,m-3,n) &
                       + a60(5)*wf(l,m-4,n)+a60(6)*wf(l,m-5,n) &
                       + a60(7)*wf(l,m-6,n)
          dpeta(l,m,n)=  a60(1)*pf(l,m  ,n)+a60(2)*pf(l,m-1,n) &
                       + a60(3)*pf(l,m-2,n)+a60(4)*pf(l,m-3,n) &
                       + a60(5)*pf(l,m-4,n)+a60(6)*pf(l,m-5,n) &
                       + a60(7)*pf(l,m-6,n)
          dreta(l,m,n)=  a60(1)*rf(l,m  ,n)+a60(2)*rf(l,m-1,n) &
                       + a60(3)*rf(l,m-2,n)+a60(4)*rf(l,m-3,n) &
                       + a60(5)*rf(l,m-4,n)+a60(6)*rf(l,m-5,n) &
                       + a60(7)*rf(l,m-6,n)
       enddo
    enddo

    ! Non-centered derivatives in z-direction *cos(teta)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    do k=nzmnghp1,nz-3
       n=k-nzmngh
       do l=1,ngh
          do m=1,ngh
             duphi(l,m,n)= a7(1)*(uf(l,m,n+1) - uf(l,m,n-1)) + &
                           a7(2)*(uf(l,m,n+2) - uf(l,m,n-2)) + &
                           a7(3)*(uf(l,m,n+3) - uf(l,m,n-3))
             dvphi(l,m,n)= a7(1)*(vf(l,m,n+1) - vf(l,m,n-1)) + &
                           a7(2)*(vf(l,m,n+2) - vf(l,m,n-2)) + &
                           a7(3)*(vf(l,m,n+3) - vf(l,m,n-3))
             dwphi(l,m,n)= a7(1)*(wf(l,m,n+1) - wf(l,m,n-1)) + &
                           a7(2)*(wf(l,m,n+2) - wf(l,m,n-2)) + &
                           a7(3)*(wf(l,m,n+3) - wf(l,m,n-3))
             dpphi(l,m,n)= a7(1)*(pf(l,m,n+1) - pf(l,m,n-1)) + &
                           a7(2)*(pf(l,m,n+2) - pf(l,m,n-2)) + &
                           a7(3)*(pf(l,m,n+3) - pf(l,m,n-3))
             drphi(l,m,n)= a7(1)*(rf(l,m,n+1) - rf(l,m,n-1)) + &
                           a7(2)*(rf(l,m,n+2) - rf(l,m,n-2)) + &
                           a7(3)*(rf(l,m,n+3) - rf(l,m,n-3))
          enddo
       enddo
    enddo

    k=nz-2
    n=k-nzmngh
    do l=1,ngh
       do m=1,ngh
          duphi(l,m,n)=  a42(1)*uf(l,m,n+2)+a42(2)*uf(l,m,n+1) &
                       + a42(3)*uf(l,m,n  )+a42(4)*uf(l,m,n-1) &
                       + a42(5)*uf(l,m,n-2)+a42(6)*uf(l,m,n-3) &
                       + a42(7)*uf(l,m,n-4)
          dvphi(l,m,n)=  a42(1)*vf(l,m,n+2)+a42(2)*vf(l,m,n+1) &
                       + a42(3)*vf(l,m,n  )+a42(4)*vf(l,m,n-1) &
                       + a42(5)*vf(l,m,n-2)+a42(6)*vf(l,m,n-3) &
                       + a42(7)*vf(l,m,n-4)
          dwphi(l,m,n)=  a42(1)*wf(l,m,n+2)+a42(2)*wf(l,m,n+1) &
                       + a42(3)*wf(l,m,n  )+a42(4)*wf(l,m,n-1) &
                       + a42(5)*wf(l,m,n-2)+a42(6)*wf(l,m,n-3) &
                       + a42(7)*wf(l,m,n-4)
          dpphi(l,m,n)=  a42(1)*pf(l,m,n+2)+a42(2)*pf(l,m,n+1) &
                       + a42(3)*pf(l,m,n  )+a42(4)*pf(l,m,n-1) &
                       + a42(5)*pf(l,m,n-2)+a42(6)*pf(l,m,n-3) &
                       + a42(7)*pf(l,m,n-4)
          drphi(l,m,n)=  a42(1)*rf(l,m,n+2)+a42(2)*rf(l,m,n+1) &
                       + a42(3)*rf(l,m,n  )+a42(4)*rf(l,m,n-1) &
                       + a42(5)*rf(l,m,n-2)+a42(6)*rf(l,m,n-3) &
                       + a42(7)*rf(l,m,n-4)
       enddo
    enddo

    k=nz-1
    n=k-nzmngh
    do l=1,ngh
       do m=1,ngh
          duphi(l,m,n)=  a51(1)*uf(l,m,n+1)+a51(2)*uf(l,m,n  ) &
                       + a51(3)*uf(l,m,n-1)+a51(4)*uf(l,m,n-2) &
                       + a51(5)*uf(l,m,n-3)+a51(6)*uf(l,m,n-4) &
                       + a51(7)*uf(l,m,n-5)
          dvphi(l,m,n)=  a51(1)*vf(l,m,n+1)+a51(2)*vf(l,m,n  ) &
                       + a51(3)*vf(l,m,n-1)+a51(4)*vf(l,m,n-2) &
                       + a51(5)*vf(l,m,n-3)+a51(6)*vf(l,m,n-4) &
                       + a51(7)*vf(l,m,n-5)
          dwphi(l,m,n)=  a51(1)*wf(l,m,n+1)+a51(2)*wf(l,m,n  ) &
                       + a51(3)*wf(l,m,n-1)+a51(4)*wf(l,m,n-2) &
                       + a51(5)*wf(l,m,n-3)+a51(6)*wf(l,m,n-4) &
                       + a51(7)*wf(l,m,n-5)
          dpphi(l,m,n)=  a51(1)*pf(l,m,n+1)+a51(2)*pf(l,m,n  ) &
                       + a51(3)*pf(l,m,n-1)+a51(4)*pf(l,m,n-2) &
                       + a51(5)*pf(l,m,n-3)+a51(6)*pf(l,m,n-4) &
                       + a51(7)*pf(l,m,n-5)
          drphi(l,m,n)=  a51(1)*rf(l,m,n+1)+a51(2)*rf(l,m,n  ) &
                       + a51(3)*rf(l,m,n-1)+a51(4)*rf(l,m,n-2) &
                       + a51(5)*rf(l,m,n-3)+a51(6)*rf(l,m,n-4) &
                       + a51(7)*rf(l,m,n-5)
       enddo
    enddo

    k=nz
    n=k-nzmngh
    do l=1,ngh
       do m=1,ngh
          duphi(l,m,n)=  a60(1)*uf(l,m,n  )+a60(2)*uf(l,m,n-1) &
                       + a60(3)*uf(l,m,n-2)+a60(4)*uf(l,m,n-3) &
                       + a60(5)*uf(l,m,n-4)+a60(6)*uf(l,m,n-5) &
                       + a60(7)*uf(l,m,n-6)
          dvphi(l,m,n)=  a60(1)*vf(l,m,n  )+a60(2)*vf(l,m,n-1) &
                       + a60(3)*vf(l,m,n-2)+a60(4)*vf(l,m,n-3) &
                       + a60(5)*vf(l,m,n-4)+a60(6)*vf(l,m,n-5) &
                       + a60(7)*vf(l,m,n-6)
          dwphi(l,m,n)=  a60(1)*wf(l,m,n  )+a60(2)*wf(l,m,n-1) &
                       + a60(3)*wf(l,m,n-2)+a60(4)*wf(l,m,n-3) &
                       + a60(5)*wf(l,m,n-4)+a60(6)*wf(l,m,n-5) &
                       + a60(7)*wf(l,m,n-6)
          dpphi(l,m,n)=  a60(1)*pf(l,m,n  )+a60(2)*pf(l,m,n-1) &
                       + a60(3)*pf(l,m,n-2)+a60(4)*pf(l,m,n-3) &
                       + a60(5)*pf(l,m,n-4)+a60(6)*pf(l,m,n-5) &
                       + a60(7)*pf(l,m,n-6)
          drphi(l,m,n)=  a60(1)*rf(l,m,n  )+a60(2)*rf(l,m,n-1) &
                       + a60(3)*rf(l,m,n-2)+a60(4)*rf(l,m,n-3) &
                       + a60(5)*rf(l,m,n-4)+a60(6)*rf(l,m,n-5) &
                       + a60(7)*rf(l,m,n-6)
       enddo
    enddo

    ! Compute rt,ut,vt,wt & pt
    ! ========================
    do i=nxmnghp1,nx
       l=i-nxmngh
       do m=1,ngh
          do k=nzmnghp1,nz
             n=k-nzmngh
             pt(l,m,n) = vg(l,m,n)*( &
                       ((dpksi(l,m,n)*ksi_x(i,j,k)+dpeta(l,m,n)*eta_x(i,j,k)+dpphi(l,m,n)*phi_x(i,j,k))*sintcosp(l,m,n) &
                       +(dpksi(l,m,n)*ksi_y(i,j,k)+dpeta(l,m,n)*eta_y(i,j,k)+dpphi(l,m,n)*phi_y(i,j,k))*sintsinp(l,m,n) &
                       +(dpksi(l,m,n)*ksi_z(i,j,k)+dpeta(l,m,n)*eta_z(i,j,k)+dpphi(l,m,n)*phi_z(i,j,k))*costeta(l,m,n)  &
                       )*ijacob3(i,j,k) + pf(l,m,n)*ir(l,m,n) )
             ut(l,m,n) = vg(l,m,n)*( &
                       ((duksi(l,m,n)*ksi_x(i,j,k)+dueta(l,m,n)*eta_x(i,j,k)+duphi(l,m,n)*phi_x(i,j,k))*sintcosp(l,m,n) &
                       +(duksi(l,m,n)*ksi_y(i,j,k)+dueta(l,m,n)*eta_y(i,j,k)+duphi(l,m,n)*phi_y(i,j,k))*sintsinp(l,m,n) &
                       +(duksi(l,m,n)*ksi_z(i,j,k)+dueta(l,m,n)*eta_z(i,j,k)+duphi(l,m,n)*phi_z(i,j,k))*costeta(l,m,n)  &
                       )*ijacob3(i,j,k) + uf(l,m,n)*ir(l,m,n) )
             vt(l,m,n) = vg(l,m,n)*( &
                       ((dvksi(l,m,n)*ksi_x(i,j,k)+dveta(l,m,n)*eta_x(i,j,k)+dvphi(l,m,n)*phi_x(i,j,k))*sintcosp(l,m,n) &
                       +(dvksi(l,m,n)*ksi_y(i,j,k)+dveta(l,m,n)*eta_y(i,j,k)+dvphi(l,m,n)*phi_y(i,j,k))*sintsinp(l,m,n) &
                       +(dvksi(l,m,n)*ksi_z(i,j,k)+dveta(l,m,n)*eta_z(i,j,k)+dvphi(l,m,n)*phi_z(i,j,k))*costeta(l,m,n)  &
                       )*ijacob3(i,j,k) + vf(l,m,n)*ir(l,m,n) )
             wt(l,m,n) = vg(l,m,n)*( &
                       ((dwksi(l,m,n)*ksi_x(i,j,k)+dweta(l,m,n)*eta_x(i,j,k)+dwphi(l,m,n)*phi_x(i,j,k))*sintcosp(l,m,n) &
                       +(dwksi(l,m,n)*ksi_y(i,j,k)+dweta(l,m,n)*eta_y(i,j,k)+dwphi(l,m,n)*phi_y(i,j,k))*sintsinp(l,m,n) &
                       +(dwksi(l,m,n)*ksi_z(i,j,k)+dweta(l,m,n)*eta_z(i,j,k)+dwphi(l,m,n)*phi_z(i,j,k))*costeta(l,m,n)  &
                       )*ijacob3(i,j,k) + wf(l,m,n)*ir(l,m,n) )
             rt(l,m,n) = vg(l,m,n)*( &
                       ((drksi(l,m,n)*ksi_x(i,j,k)+dreta(l,m,n)*eta_x(i,j,k)+drphi(l,m,n)*phi_x(i,j,k))*sintcosp(l,m,n) &
                       +(drksi(l,m,n)*ksi_y(i,j,k)+dreta(l,m,n)*eta_y(i,j,k)+drphi(l,m,n)*phi_y(i,j,k))*sintsinp(l,m,n) &
                       +(drksi(l,m,n)*ksi_z(i,j,k)+dreta(l,m,n)*eta_z(i,j,k)+drphi(l,m,n)*phi_z(i,j,k))*costeta(l,m,n)  &
                       )*ijacob3(i,j,k) + rf(l,m,n)*ir(l,m,n) )
         enddo
       enddo
    enddo

    ! Update fluxes at each RK step
    ! =============================
    do k=nzmnghp1,nz
       n=k-nzmngh
       do i=nxmnghp1,nx
          l=i-nxmngh
          do j=nymnghp1,ny
             m=j-nymngh
             cp=cpcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
             av=avcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
             c2_=c2calc_tro(Tmp(i,j,k),rho_n(i,j,k))

             Krho(i,j,k)  = rt(l,m,n)
             Krhou(i,j,k) = uu(i,j,k)*rt(l,m,n)+rho_n(i,j,k)*ut(l,m,n)
             Krhov(i,j,k) = vv(i,j,k)*rt(l,m,n)+rho_n(i,j,k)*vt(l,m,n)
             Krhow(i,j,k) = ww(i,j,k)*rt(l,m,n)+rho_n(i,j,k)*wt(l,m,n)
             Krhoe(i,j,k) = cp/av*(pt(l,m,n)/c2_-rt(l,m,n)) &
                  + (rhoe_n(i,j,k)+prs(i,j,k))/rho_n(i,j,k)*rt(l,m,n) &
                  + rho_n(i,j,k)*(uu(i,j,k)*ut(l,m,n)+vv(i,j,k)*vt(l,m,n)+ww(i,j,k)*wt(l,m,n))
          enddo
       enddo
    enddo

  end subroutine bc_TD3d_imax_jmax_kmax_c3

end submodule smod_TamDong3d_corner_c3
