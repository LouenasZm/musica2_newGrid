!===============================================================================
submodule (mod_TamDong3d) smod_TamDong3d_edgex
!===============================================================================
  !> author: XG
  !> date: February 2020 - modif January 2022
  !> Non-reflecting boundary conditions of Tam and Dong
  !> - 3D Cartesian version - edge along x
!=============================================================================== 

contains

  !===============================================================================
  module subroutine bc_TD3d_jmin_kmin
  !===============================================================================
    !> 3D Tam & Dong's BC: boundary condition at jmin-kmin (edge 2,1,1 /bottom-front)
    !> - Cartesian version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: cp,av,c2_
    real(wp), dimension(nx1:nx2,1:ngh+3,1:ngh+3) :: rf,uf,vf,wf,pf
    real(wp), dimension(nx,1:ngh,1:ngh) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(nx,1:ngh,1:ngh) :: duy,dvy,dwy,dpy,dry
    real(wp), dimension(nx,1:ngh,1:ngh) :: duz,dvz,dwz,dpz,drz
    real(wp) :: dux,dvx,dwx,dpx,drx
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
	           	BC_face(2,1)%U0(i,j,k,4)*costeta(i,j,k)	       &
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
    ! NOT(Tam & Webb DRP schemes) => order 2 (wall)
    j=1
    do k=1,ngh
       do i=ndx,nfx
          duy(i,j,k)=(uf(i,j+1,k)-uf(i,j,k))*idy1_jmin*sintsinp(i,j,k)
          dvy(i,j,k)=(vf(i,j+1,k)-vf(i,j,k))*idy1_jmin*sintsinp(i,j,k)
          dwy(i,j,k)=(wf(i,j+1,k)-wf(i,j,k))*idy1_jmin*sintsinp(i,j,k)
          dpy(i,j,k)=(pf(i,j+1,k)-pf(i,j,k))*idy1_jmin*sintsinp(i,j,k)
          dry(i,j,k)=(rf(i,j+1,k)-rf(i,j,k))*idy1_jmin*sintsinp(i,j,k)
       enddo
    enddo

    j=2
    do k=1,ngh
       do i=ndx,nfx
          duy(i,j,k)=a3(1)*(uf(i,j+1,k)-uf(i,j-1,k))*idy2_jmin*sintsinp(i,j,k)
          dvy(i,j,k)=a3(1)*(vf(i,j+1,k)-vf(i,j-1,k))*idy2_jmin*sintsinp(i,j,k)
          dwy(i,j,k)=a3(1)*(wf(i,j+1,k)-wf(i,j-1,k))*idy2_jmin*sintsinp(i,j,k)
          dpy(i,j,k)=a3(1)*(pf(i,j+1,k)-pf(i,j-1,k))*idy2_jmin*sintsinp(i,j,k)
          dry(i,j,k)=a3(1)*(rf(i,j+1,k)-rf(i,j-1,k))*idy2_jmin*sintsinp(i,j,k)
       enddo
    enddo

    j=3
    do k=1,ngh
       do i=ndx,nfx
          duy(i,j,k)= ( a5(1)*(uf(i,j+1,k)-uf(i,j-1,k)) &
                      + a5(2)*(uf(i,j+2,k)-uf(i,j-2,k)) ) *idy4_jmin*sintsinp(i,j,k)
          dvy(i,j,k)= ( a5(1)*(vf(i,j+1,k)-vf(i,j-1,k)) &
                      + a5(2)*(vf(i,j+2,k)-vf(i,j-2,k)) ) *idy4_jmin*sintsinp(i,j,k)
          dwy(i,j,k)= ( a5(1)*(wf(i,j+1,k)-wf(i,j-1,k)) &
                      + a5(2)*(wf(i,j+2,k)-wf(i,j-2,k)) ) *idy4_jmin*sintsinp(i,j,k)
          dpy(i,j,k)= ( a5(1)*(pf(i,j+1,k)-pf(i,j-1,k)) &
                      + a5(2)*(pf(i,j+2,k)-pf(i,j-2,k)) ) *idy4_jmin*sintsinp(i,j,k)
          dry(i,j,k)= ( a5(1)*(rf(i,j+1,k)-rf(i,j-1,k)) &
                      + a5(2)*(rf(i,j+2,k)-rf(i,j-2,k)) ) *idy4_jmin*sintsinp(i,j,k)
       enddo
    enddo
!!$
!!$    j=1
!!$    do k=1,ngh
!!$       do i=ndx,nfx
!!$          duy(i,j,k)= ( a06(1)*uf(i,1,k)+a06(2)*uf(i,2,k) &
!!$                      + a06(3)*uf(i,3,k)+a06(4)*uf(i,4,k) &
!!$                      + a06(5)*uf(i,5,k)+a06(6)*uf(i,6,k) &
!!$                      + a06(7)*uf(i,7,k) ) *idy(j)*sintsinp(i,j,k)
!!$          dvy(i,j,k)= ( a06(1)*vf(i,1,k)+a06(2)*vf(i,2,k) &
!!$                      + a06(3)*vf(i,3,k)+a06(4)*vf(i,4,k) &
!!$                      + a06(5)*vf(i,5,k)+a06(6)*vf(i,6,k) &
!!$                      + a06(7)*vf(i,7,k) ) *idy(j)*sintsinp(i,j,k)
!!$          dwy(i,j,k)= ( a06(1)*wf(i,1,k)+a06(2)*wf(i,2,k) &
!!$                      + a06(3)*wf(i,3,k)+a06(4)*wf(i,4,k) &
!!$                      + a06(5)*wf(i,5,k)+a06(6)*wf(i,6,k) &
!!$                      + a06(7)*wf(i,7,k) ) *idy(j)*sintsinp(i,j,k)
!!$          dpy(i,j,k)= ( a06(1)*pf(i,1,k)+a06(2)*pf(i,2,k) &
!!$                      + a06(3)*pf(i,3,k)+a06(4)*pf(i,4,k) &
!!$                      + a06(5)*pf(i,5,k)+a06(6)*pf(i,6,k) &
!!$                      + a06(7)*pf(i,7,k) ) *idy(j)*sintsinp(i,j,k)
!!$          dry(i,j,k)= ( a06(1)*rf(i,1,k)+a06(2)*rf(i,2,k) &
!!$                      + a06(3)*rf(i,3,k)+a06(4)*rf(i,4,k) &
!!$                      + a06(5)*rf(i,5,k)+a06(6)*rf(i,6,k) &
!!$                      + a06(7)*rf(i,7,k) ) *idy(j)*sintsinp(i,j,k)
!!$       enddo
!!$    enddo
!!$
!!$    j=2
!!$    do k=1,ngh
!!$       do i=ndx,nfx
!!$          duy(i,j,k)= ( a15(1)*uf(i,1,k)+a15(2)*uf(i,2,k) &
!!$                      + a15(3)*uf(i,3,k)+a15(4)*uf(i,4,k) &
!!$                      + a15(5)*uf(i,5,k)+a15(6)*uf(i,6,k) &
!!$                      + a15(7)*uf(i,7,k) ) *idy(j)*sintsinp(i,j,k)
!!$          dvy(i,j,k)= ( a15(1)*vf(i,1,k)+a15(2)*vf(i,2,k) &
!!$                      + a15(3)*vf(i,3,k)+a15(4)*vf(i,4,k) &
!!$                      + a15(5)*vf(i,5,k)+a15(6)*vf(i,6,k) &
!!$                      + a15(7)*vf(i,7,k) ) *idy(j)*sintsinp(i,j,k)
!!$          dwy(i,j,k)= ( a15(1)*wf(i,1,k)+a15(2)*wf(i,2,k) &
!!$                      + a15(3)*wf(i,3,k)+a15(4)*wf(i,4,k) &
!!$                      + a15(5)*wf(i,5,k)+a15(6)*wf(i,6,k) &
!!$                      + a15(7)*wf(i,7,k) ) *idy(j)*sintsinp(i,j,k)
!!$          dpy(i,j,k)= ( a15(1)*pf(i,1,k)+a15(2)*pf(i,2,k) &
!!$                      + a15(3)*pf(i,3,k)+a15(4)*pf(i,4,k) &
!!$                      + a15(5)*pf(i,5,k)+a15(6)*pf(i,6,k) &
!!$                      + a15(7)*pf(i,7,k) ) *idy(j)*sintsinp(i,j,k)
!!$          dry(i,j,k)= ( a15(1)*rf(i,1,k)+a15(2)*rf(i,2,k) &
!!$                      + a15(3)*rf(i,3,k)+a15(4)*rf(i,4,k) &
!!$                      + a15(5)*rf(i,5,k)+a15(6)*rf(i,6,k) &
!!$                      + a15(7)*rf(i,7,k) ) *idy(j)*sintsinp(i,j,k)
!!$       enddo
!!$    enddo
!!$
!!$    j=3
!!$    do k=1,ngh
!!$       do i=ndx,nfx
!!$          duy(i,j,k)= ( a24(1)*uf(i,1,k)+a24(2)*uf(i,2,k) &
!!$                      + a24(3)*uf(i,3,k)+a24(4)*uf(i,4,k) &
!!$                      + a24(5)*uf(i,5,k)+a24(6)*uf(i,6,k) &
!!$                      + a24(7)*uf(i,7,k) ) *idy(j)*sintsinp(i,j,k)
!!$          dvy(i,j,k)= ( a24(1)*vf(i,1,k)+a24(2)*vf(i,2,k) &
!!$                      + a24(3)*vf(i,3,k)+a24(4)*vf(i,4,k) &
!!$                      + a24(5)*vf(i,5,k)+a24(6)*vf(i,6,k) &
!!$                      + a24(7)*vf(i,7,k) ) *idy(j)*sintsinp(i,j,k)
!!$          dwy(i,j,k)= ( a24(1)*wf(i,1,k)+a24(2)*wf(i,2,k) &
!!$                      + a24(3)*wf(i,3,k)+a24(4)*wf(i,4,k) &
!!$                      + a24(5)*wf(i,5,k)+a24(6)*wf(i,6,k) &
!!$                      + a24(7)*wf(i,7,k) ) *idy(j)*sintsinp(i,j,k)
!!$          dpy(i,j,k)= ( a24(1)*pf(i,1,k)+a24(2)*pf(i,2,k) &
!!$                      + a24(3)*pf(i,3,k)+a24(4)*pf(i,4,k) &
!!$                      + a24(5)*pf(i,5,k)+a24(6)*pf(i,6,k) &
!!$                      + a24(7)*pf(i,7,k) ) *idy(j)*sintsinp(i,j,k)
!!$          dry(i,j,k)= ( a24(1)*rf(i,1,k)+a24(2)*rf(i,2,k) &
!!$                      + a24(3)*rf(i,3,k)+a24(4)*rf(i,4,k) &
!!$                      + a24(5)*rf(i,5,k)+a24(6)*rf(i,6,k) &
!!$                      + a24(7)*rf(i,7,k) ) *idy(j)*sintsinp(i,j,k)
!!$       enddo
!!$    enddo

    do j=4,ngh
       do k=1,ngh
          do i=ndx,nfx
             duy(i,j,k)= ( a7(1)*(uf(i,j+1,k)-uf(i,j-1,k)) &
                         + a7(2)*(uf(i,j+2,k)-uf(i,j-2,k)) &
                         + a7(3)*(uf(i,j+3,k)-uf(i,j-3,k)) ) *idy(j)*sintsinp(i,j,k)
             dvy(i,j,k)= ( a7(1)*(vf(i,j+1,k)-vf(i,j-1,k)) &
                         + a7(2)*(vf(i,j+2,k)-vf(i,j-2,k)) &
                         + a7(3)*(vf(i,j+3,k)-vf(i,j-3,k)) ) *idy(j)*sintsinp(i,j,k)
             dwy(i,j,k)= ( a7(1)*(wf(i,j+1,k)-wf(i,j-1,k)) &
                         + a7(2)*(wf(i,j+2,k)-wf(i,j-2,k)) &
                         + a7(3)*(wf(i,j+3,k)-wf(i,j-3,k)) ) *idy(j)*sintsinp(i,j,k)
             dpy(i,j,k)= ( a7(1)*(pf(i,j+1,k)-pf(i,j-1,k)) &
                         + a7(2)*(pf(i,j+2,k)-pf(i,j-2,k)) &
                         + a7(3)*(pf(i,j+3,k)-pf(i,j-3,k)) ) *idy(j)*sintsinp(i,j,k)
             dry(i,j,k)= ( a7(1)*(rf(i,j+1,k)-rf(i,j-1,k)) &
                         + a7(2)*(rf(i,j+2,k)-rf(i,j-2,k)) &
                         + a7(3)*(rf(i,j+3,k)-rf(i,j-3,k)) ) *idy(j)*sintsinp(i,j,k)
          enddo
       enddo
    enddo

    ! Non-centered derivatives in z-direction *cos(teta)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    k=1
    do j=1,ngh
       do i=ndx,nfx
          duz(i,j,k) = ( a06(1)*uf(i,j,1)+a06(2)*uf(i,j,2) &
                       + a06(3)*uf(i,j,3)+a06(4)*uf(i,j,4) &
                       + a06(5)*uf(i,j,5)+a06(6)*uf(i,j,6) &
                       + a06(7)*uf(i,j,7) ) *idz(k)*costeta(i,j,k)
          dvz(i,j,k) = ( a06(1)*vf(i,j,1)+a06(2)*vf(i,j,2) &
                       + a06(3)*vf(i,j,3)+a06(4)*vf(i,j,4) &
                       + a06(5)*vf(i,j,5)+a06(6)*vf(i,j,6) &
                       + a06(7)*vf(i,j,7) ) *idz(k)*costeta(i,j,k)
          dwz(i,j,k) = ( a06(1)*wf(i,j,1)+a06(2)*wf(i,j,2) &
                       + a06(3)*wf(i,j,3)+a06(4)*wf(i,j,4) &
                       + a06(5)*wf(i,j,5)+a06(6)*wf(i,j,6) &
                       + a06(7)*wf(i,j,7) ) *idz(k)*costeta(i,j,k)
          dpz(i,j,k) = ( a06(1)*pf(i,j,1)+a06(2)*pf(i,j,2) &
                       + a06(3)*pf(i,j,3)+a06(4)*pf(i,j,4) &
                       + a06(5)*pf(i,j,5)+a06(6)*pf(i,j,6) &
                       + a06(7)*pf(i,j,7) ) *idz(k)*costeta(i,j,k)
          drz(i,j,k) = ( a06(1)*rf(i,j,1)+a06(2)*rf(i,j,2) &
                       + a06(3)*rf(i,j,3)+a06(4)*rf(i,j,4) &
                       + a06(5)*rf(i,j,5)+a06(6)*rf(i,j,6) &
                       + a06(7)*rf(i,j,7) ) *idz(k)*costeta(i,j,k)
       enddo
    enddo

    k=2
    do j=1,ngh
       do i=ndx,nfx
          duz(i,j,k) = ( a15(1)*uf(i,j,1)+a15(2)*uf(i,j,2) &
                       + a15(3)*uf(i,j,3)+a15(4)*uf(i,j,4) &
                       + a15(5)*uf(i,j,5)+a15(6)*uf(i,j,6) &
                       + a15(7)*uf(i,j,7) ) *idz(k)*costeta(i,j,k)
          dvz(i,j,k) = ( a15(1)*vf(i,j,1)+a15(2)*vf(i,j,2) &
                       + a15(3)*vf(i,j,3)+a15(4)*vf(i,j,4) &
                       + a15(5)*vf(i,j,5)+a15(6)*vf(i,j,6) &
                       + a15(7)*vf(i,j,7) ) *idz(k)*costeta(i,j,k)
          dwz(i,j,k) = ( a15(1)*wf(i,j,1)+a15(2)*wf(i,j,2) &
                       + a15(3)*wf(i,j,3)+a15(4)*wf(i,j,4) &
                       + a15(5)*wf(i,j,5)+a15(6)*wf(i,j,6) &
                       + a15(7)*wf(i,j,7) ) *idz(k)*costeta(i,j,k)
          dpz(i,j,k) = ( a15(1)*pf(i,j,1)+a15(2)*pf(i,j,2) &
                       + a15(3)*pf(i,j,3)+a15(4)*pf(i,j,4) &
                       + a15(5)*pf(i,j,5)+a15(6)*pf(i,j,6) &
                       + a15(7)*pf(i,j,7) ) *idz(k)*costeta(i,j,k)
          drz(i,j,k) = ( a15(1)*rf(i,j,1)+a15(2)*rf(i,j,2) &
                       + a15(3)*rf(i,j,3)+a15(4)*rf(i,j,4) &
                       + a15(5)*rf(i,j,5)+a15(6)*rf(i,j,6) &
                       + a15(7)*rf(i,j,7) ) *idz(k)*costeta(i,j,k)
       enddo
    enddo

    k=3
    do j=1,ngh
       do i=ndx,nfx
          duz(i,j,k) = ( a24(1)*uf(i,j,1)+a24(2)*uf(i,j,2) &
                       + a24(3)*uf(i,j,3)+a24(4)*uf(i,j,4) &
                       + a24(5)*uf(i,j,5)+a24(6)*uf(i,j,6) &
                       + a24(7)*uf(i,j,7) ) *idz(k)*costeta(i,j,k)
          dvz(i,j,k) = ( a24(1)*vf(i,j,1)+a24(2)*vf(i,j,2) &
                       + a24(3)*vf(i,j,3)+a24(4)*vf(i,j,4) &
                       + a24(5)*vf(i,j,5)+a24(6)*vf(i,j,6) &
                       + a24(7)*vf(i,j,7) ) *idz(k)*costeta(i,j,k)
          dwz(i,j,k) = ( a24(1)*wf(i,j,1)+a24(2)*wf(i,j,2) &
                       + a24(3)*wf(i,j,3)+a24(4)*wf(i,j,4) &
                       + a24(5)*wf(i,j,5)+a24(6)*wf(i,j,6) &
                       + a24(7)*wf(i,j,7) ) *idz(k)*costeta(i,j,k)
          dpz(i,j,k) = ( a24(1)*pf(i,j,1)+a24(2)*pf(i,j,2) &
                       + a24(3)*pf(i,j,3)+a24(4)*pf(i,j,4) &
                       + a24(5)*pf(i,j,5)+a24(6)*pf(i,j,6) &
                       + a24(7)*pf(i,j,7) ) *idz(k)*costeta(i,j,k)
          drz(i,j,k) = ( a24(1)*rf(i,j,1)+a24(2)*rf(i,j,2) &
                       + a24(3)*rf(i,j,3)+a24(4)*rf(i,j,4) &
                       + a24(5)*rf(i,j,5)+a24(6)*rf(i,j,6) &
                       + a24(7)*rf(i,j,7) ) *idz(k)*costeta(i,j,k)
       enddo
    enddo

    do k=4,ngh
       do j=1,ngh
          do i=ndx,nfx
             duz(i,j,k) = (a7(1)*(uf(i,j,k+1) - uf(i,j,k-1)) + &
                           a7(2)*(uf(i,j,k+2) - uf(i,j,k-2)) + &
                           a7(3)*(uf(i,j,k+3) - uf(i,j,k-3)) ) *idz(k)*costeta(i,j,k)
             dvz(i,j,k) = (a7(1)*(vf(i,j,k+1) - vf(i,j,k-1)) + &
                           a7(2)*(vf(i,j,k+2) - vf(i,j,k-2)) + &
                           a7(3)*(vf(i,j,k+3) - vf(i,j,k-3)) ) *idz(k)*costeta(i,j,k)
             dwz(i,j,k) = (a7(1)*(wf(i,j,k+1) - wf(i,j,k-1)) + &
                           a7(2)*(wf(i,j,k+2) - wf(i,j,k-2)) + &
                           a7(3)*(wf(i,j,k+3) - wf(i,j,k-3)) ) *idz(k)*costeta(i,j,k)
             dpz(i,j,k) = (a7(1)*(pf(i,j,k+1) - pf(i,j,k-1)) + &
                           a7(2)*(pf(i,j,k+2) - pf(i,j,k-2)) + &
                           a7(3)*(pf(i,j,k+3) - pf(i,j,k-3)) ) *idz(k)*costeta(i,j,k)
             drz(i,j,k) = (a7(1)*(rf(i,j,k+1) - rf(i,j,k-1)) + &
                           a7(2)*(rf(i,j,k+2) - rf(i,j,k-2)) + &
                           a7(3)*(rf(i,j,k+3) - rf(i,j,k-3)) ) *idz(k)*costeta(i,j,k)
          enddo
       enddo
    enddo

    ! Compute rt,ut,vt,wt & pt
    ! ========================
    do k=1,ngh
       do j=1,ngh
          do i=ndx,nfx
             dux = ( a7(1)*( uf(i+1,j,k) - uf(i-1,j,k) ) + &
                     a7(2)*( uf(i+2,j,k) - uf(i-2,j,k) ) + &
                     a7(3)*( uf(i+3,j,k) - uf(i-3,j,k) ) ) *idx(i)*sintcosp(i,j,k)
             dvx = ( a7(1)*( vf(i+1,j,k) - vf(i-1,j,k) ) + &
                     a7(2)*( vf(i+2,j,k) - vf(i-2,j,k) ) + &
                     a7(3)*( vf(i+3,j,k) - vf(i-3,j,k) ) ) *idx(i)*sintcosp(i,j,k)
             dwx = ( a7(1)*( wf(i+1,j,k) - wf(i-1,j,k) ) + &
                     a7(2)*( wf(i+2,j,k) - wf(i-2,j,k) ) + &
                     a7(3)*( wf(i+3,j,k) - wf(i-3,j,k) ) ) *idx(i)*sintcosp(i,j,k)
             dpx = ( a7(1)*( pf(i+1,j,k) - pf(i-1,j,k) ) + &
                     a7(2)*( pf(i+2,j,k) - pf(i-2,j,k) ) + &
                     a7(3)*( pf(i+3,j,k) - pf(i-3,j,k) ) ) *idx(i)*sintcosp(i,j,k)
             drx = ( a7(1)*( rf(i+1,j,k) - rf(i-1,j,k) ) + &
                     a7(2)*( rf(i+2,j,k) - rf(i-2,j,k) ) + &
                     a7(3)*( rf(i+3,j,k) - rf(i-3,j,k) ) ) *idx(i)*sintcosp(i,j,k)

             pt(i,j,k) = vg(i,j,k)*(dpx+dpy(i,j,k)+dpz(i,j,k)+pf(i,j,k)*ir(i,j,k))
             ut(i,j,k) = vg(i,j,k)*(dux+duy(i,j,k)+duz(i,j,k)+uf(i,j,k)*ir(i,j,k))
             vt(i,j,k) = vg(i,j,k)*(dvx+dvy(i,j,k)+dvz(i,j,k)+vf(i,j,k)*ir(i,j,k))
             wt(i,j,k) = vg(i,j,k)*(dwx+dwy(i,j,k)+dwz(i,j,k)+wf(i,j,k)*ir(i,j,k))
             rt(i,j,k) = vg(i,j,k)*(drx+dry(i,j,k)+drz(i,j,k)+rf(i,j,k)*ir(i,j,k))
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

  end subroutine bc_TD3d_jmin_kmin

  !===============================================================================
  module subroutine bc_TD3d_jmin_kmax
  !===============================================================================
    !> 3D Tam & Dong's BC: boundary condition at jmin-kmax (edge 2,1,2 /bottom-back)
    !> - Cartesian version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l
    real(wp) :: cp,av,c2_
    real(wp), dimension(nx1:nx2,1:ngh+3,-2:ngh) :: rf,uf,vf,wf,pf
    real(wp), dimension(nx,1:ngh,1:ngh) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(nx,1:ngh,1:ngh) :: duy,dvy,dwy,dpy,dry
    real(wp), dimension(nx,1:ngh,1:ngh) :: duz,dvz,dwz,dpz,drz
    real(wp) :: dux,dvx,dwx,dpx,drx
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
	           	BC_face(2,1)%U0(i,j,k,4)*costeta(i,j,l)	       &
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
    ! (Tam & Webb DRP schemes)
    ! NOT(Tam & Webb DRP schemes) => order 2 (wall)
    j=1
    do l=1,ngh
       do i=ndx,nfx
          duy(i,j,l)=(uf(i,j+1,l)-uf(i,j,l))*idy1_jmin*sintsinp(i,j,l)
          dvy(i,j,l)=(vf(i,j+1,l)-vf(i,j,l))*idy1_jmin*sintsinp(i,j,l)
          dwy(i,j,l)=(wf(i,j+1,l)-wf(i,j,l))*idy1_jmin*sintsinp(i,j,l)
          dpy(i,j,l)=(pf(i,j+1,l)-pf(i,j,l))*idy1_jmin*sintsinp(i,j,l)
          dry(i,j,l)=(rf(i,j+1,l)-rf(i,j,l))*idy1_jmin*sintsinp(i,j,l)
       enddo
    enddo

    j=2
    do l=1,ngh
       do i=ndx,nfx
          duy(i,j,l)=a3(1)*(uf(i,j+1,l)-uf(i,j-1,l))*idy2_jmin*sintsinp(i,j,l)
          dvy(i,j,l)=a3(1)*(vf(i,j+1,l)-vf(i,j-1,l))*idy2_jmin*sintsinp(i,j,l)
          dwy(i,j,l)=a3(1)*(wf(i,j+1,l)-wf(i,j-1,l))*idy2_jmin*sintsinp(i,j,l)
          dpy(i,j,l)=a3(1)*(pf(i,j+1,l)-pf(i,j-1,l))*idy2_jmin*sintsinp(i,j,l)
          dry(i,j,l)=a3(1)*(rf(i,j+1,l)-rf(i,j-1,l))*idy2_jmin*sintsinp(i,j,l)
       enddo
    enddo

    j=3
    do l=1,ngh
       do i=ndx,nfx
          duy(i,j,l)= ( a5(1)*(uf(i,j+1,l)-uf(i,j-1,l)) &
                      + a5(2)*(uf(i,j+2,l)-uf(i,j-2,l)) ) *idy4_jmin*sintsinp(i,j,l)
          dvy(i,j,l)= ( a5(1)*(vf(i,j+1,l)-vf(i,j-1,l)) &
                      + a5(2)*(vf(i,j+2,l)-vf(i,j-2,l)) ) *idy4_jmin*sintsinp(i,j,l)
          dwy(i,j,l)= ( a5(1)*(wf(i,j+1,l)-wf(i,j-1,l)) &
                      + a5(2)*(wf(i,j+2,l)-wf(i,j-2,l)) ) *idy4_jmin*sintsinp(i,j,l)
          dpy(i,j,l)= ( a5(1)*(pf(i,j+1,l)-pf(i,j-1,l)) &
                      + a5(2)*(pf(i,j+2,l)-pf(i,j-2,l)) ) *idy4_jmin*sintsinp(i,j,l)
          dry(i,j,l)= ( a5(1)*(rf(i,j+1,l)-rf(i,j-1,l)) &
                      + a5(2)*(rf(i,j+2,l)-rf(i,j-2,l)) ) *idy4_jmin*sintsinp(i,j,l)
       enddo
    enddo
!!$
!!$    j=1
!!$    do l=1,ngh
!!$       do i=ndx,nfx
!!$          duy(i,j,l)= ( a06(1)*uf(i,1,l)+a06(2)*uf(i,2,l) &
!!$                      + a06(3)*uf(i,3,l)+a06(4)*uf(i,4,l) &
!!$                      + a06(5)*uf(i,5,l)+a06(6)*uf(i,6,l) &
!!$                      + a06(7)*uf(i,7,l) ) *idy(j)*sintsinp(i,j,l)
!!$          dvy(i,j,l)= ( a06(1)*vf(i,1,l)+a06(2)*vf(i,2,l) &
!!$                      + a06(3)*vf(i,3,l)+a06(4)*vf(i,4,l) &
!!$                      + a06(5)*vf(i,5,l)+a06(6)*vf(i,6,l) &
!!$                      + a06(7)*vf(i,7,l) ) *idy(j)*sintsinp(i,j,l)
!!$          dwy(i,j,l)= ( a06(1)*wf(i,1,l)+a06(2)*wf(i,2,l) &
!!$                      + a06(3)*wf(i,3,l)+a06(4)*wf(i,4,l) &
!!$                      + a06(5)*wf(i,5,l)+a06(6)*wf(i,6,l) &
!!$                      + a06(7)*wf(i,7,l) ) *idy(j)*sintsinp(i,j,l)
!!$          dpy(i,j,l)= ( a06(1)*pf(i,1,l)+a06(2)*pf(i,2,l) &
!!$                      + a06(3)*pf(i,3,l)+a06(4)*pf(i,4,l) &
!!$                      + a06(5)*pf(i,5,l)+a06(6)*pf(i,6,l) &
!!$                      + a06(7)*pf(i,7,l) ) *idy(j)*sintsinp(i,j,l)
!!$          dry(i,j,l)= ( a06(1)*rf(i,1,l)+a06(2)*rf(i,2,l) &
!!$                      + a06(3)*rf(i,3,l)+a06(4)*rf(i,4,l) &
!!$                      + a06(5)*rf(i,5,l)+a06(6)*rf(i,6,l) &
!!$                      + a06(7)*rf(i,7,l) ) *idy(j)*sintsinp(i,j,l)
!!$       enddo
!!$    enddo
!!$
!!$    j=2
!!$    do l=1,ngh
!!$       do i=ndx,nfx
!!$          duy(i,j,l)= ( a15(1)*uf(i,1,l)+a15(2)*uf(i,2,l) &
!!$                      + a15(3)*uf(i,3,l)+a15(4)*uf(i,4,l) &
!!$                      + a15(5)*uf(i,5,l)+a15(6)*uf(i,6,l) &
!!$                      + a15(7)*uf(i,7,l) ) *idy(j)*sintsinp(i,j,l)
!!$          dvy(i,j,l)= ( a15(1)*vf(i,1,l)+a15(2)*vf(i,2,l) &
!!$                      + a15(3)*vf(i,3,l)+a15(4)*vf(i,4,l) &
!!$                      + a15(5)*vf(i,5,l)+a15(6)*vf(i,6,l) &
!!$                      + a15(7)*vf(i,7,l) ) *idy(j)*sintsinp(i,j,l)
!!$          dwy(i,j,l)= ( a15(1)*wf(i,1,l)+a15(2)*wf(i,2,l) &
!!$                      + a15(3)*wf(i,3,l)+a15(4)*wf(i,4,l) &
!!$                      + a15(5)*wf(i,5,l)+a15(6)*wf(i,6,l) &
!!$                      + a15(7)*wf(i,7,l) ) *idy(j)*sintsinp(i,j,l)
!!$          dpy(i,j,l)= ( a15(1)*pf(i,1,l)+a15(2)*pf(i,2,l) &
!!$                      + a15(3)*pf(i,3,l)+a15(4)*pf(i,4,l) &
!!$                      + a15(5)*pf(i,5,l)+a15(6)*pf(i,6,l) &
!!$                      + a15(7)*pf(i,7,l) ) *idy(j)*sintsinp(i,j,l)
!!$          dry(i,j,l)= ( a15(1)*rf(i,1,l)+a15(2)*rf(i,2,l) &
!!$                      + a15(3)*rf(i,3,l)+a15(4)*rf(i,4,l) &
!!$                      + a15(5)*rf(i,5,l)+a15(6)*rf(i,6,l) &
!!$                      + a15(7)*rf(i,7,l) ) *idy(j)*sintsinp(i,j,l)
!!$       enddo
!!$    enddo
!!$
!!$    j=3
!!$    do l=1,ngh
!!$       do i=ndx,nfx
!!$          duy(i,j,l)= ( a24(1)*uf(i,1,l)+a24(2)*uf(i,2,l) &
!!$                      + a24(3)*uf(i,3,l)+a24(4)*uf(i,4,l) &
!!$                      + a24(5)*uf(i,5,l)+a24(6)*uf(i,6,l) &
!!$                      + a24(7)*uf(i,7,l) ) *idy(j)*sintsinp(i,j,l)
!!$          dvy(i,j,l)= ( a24(1)*vf(i,1,l)+a24(2)*vf(i,2,l) &
!!$                      + a24(3)*vf(i,3,l)+a24(4)*vf(i,4,l) &
!!$                      + a24(5)*vf(i,5,l)+a24(6)*vf(i,6,l) &
!!$                      + a24(7)*vf(i,7,l) ) *idy(j)*sintsinp(i,j,l)
!!$          dwy(i,j,l)= ( a24(1)*wf(i,1,l)+a24(2)*wf(i,2,l) &
!!$                      + a24(3)*wf(i,3,l)+a24(4)*wf(i,4,l) &
!!$                      + a24(5)*wf(i,5,l)+a24(6)*wf(i,6,l) &
!!$                      + a24(7)*wf(i,7,l) ) *idy(j)*sintsinp(i,j,l)
!!$          dpy(i,j,l)= ( a24(1)*pf(i,1,l)+a24(2)*pf(i,2,l) &
!!$                      + a24(3)*pf(i,3,l)+a24(4)*pf(i,4,l) &
!!$                      + a24(5)*pf(i,5,l)+a24(6)*pf(i,6,l) &
!!$                      + a24(7)*pf(i,7,l) ) *idy(j)*sintsinp(i,j,l)
!!$          dry(i,j,l)= ( a24(1)*rf(i,1,l)+a24(2)*rf(i,2,l) &
!!$                      + a24(3)*rf(i,3,l)+a24(4)*rf(i,4,l) &
!!$                      + a24(5)*rf(i,5,l)+a24(6)*rf(i,6,l) &
!!$                      + a24(7)*rf(i,7,l) ) *idy(j)*sintsinp(i,j,l)
!!$       enddo
!!$    enddo

    do j=4,ngh
       do l=1,ngh
          do i=ndx,nfx
             duy(i,j,l)= ( a7(1)*(uf(i,j+1,l)-uf(i,j-1,l)) &
                         + a7(2)*(uf(i,j+2,l)-uf(i,j-2,l)) &
                         + a7(3)*(uf(i,j+3,l)-uf(i,j-3,l)) ) *idy(j)*sintsinp(i,j,l)
             dvy(i,j,l)= ( a7(1)*(vf(i,j+1,l)-vf(i,j-1,l)) &
                         + a7(2)*(vf(i,j+2,l)-vf(i,j-2,l)) &
                         + a7(3)*(vf(i,j+3,l)-vf(i,j-3,l)) ) *idy(j)*sintsinp(i,j,l)
             dwy(i,j,l)= ( a7(1)*(wf(i,j+1,l)-wf(i,j-1,l)) &
                         + a7(2)*(wf(i,j+2,l)-wf(i,j-2,l)) &
                         + a7(3)*(wf(i,j+3,l)-wf(i,j-3,l)) ) *idy(j)*sintsinp(i,j,l)
             dpy(i,j,l)= ( a7(1)*(pf(i,j+1,l)-pf(i,j-1,l)) &
                         + a7(2)*(pf(i,j+2,l)-pf(i,j-2,l)) &
                         + a7(3)*(pf(i,j+3,l)-pf(i,j-3,l)) ) *idy(j)*sintsinp(i,j,l)
             dry(i,j,l)= ( a7(1)*(rf(i,j+1,l)-rf(i,j-1,l)) &
                         + a7(2)*(rf(i,j+2,l)-rf(i,j-2,l)) &
                         + a7(3)*(rf(i,j+3,l)-rf(i,j-3,l)) ) *idy(j)*sintsinp(i,j,l)
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
             duz(i,j,l) = (a7(1)*(uf(i,j,l+1) - uf(i,j,l-1)) + &
                           a7(2)*(uf(i,j,l+2) - uf(i,j,l-2)) + &
                           a7(3)*(uf(i,j,l+3) - uf(i,j,l-3)) ) *idz(k)*costeta(i,j,l)
             dvz(i,j,l) = (a7(1)*(vf(i,j,l+1) - vf(i,j,l-1)) + &
                           a7(2)*(vf(i,j,l+2) - vf(i,j,l-2)) + &
                           a7(3)*(vf(i,j,l+3) - vf(i,j,l-3)) ) *idz(k)*costeta(i,j,l)
             dwz(i,j,l) = (a7(1)*(wf(i,j,l+1) - wf(i,j,l-1)) + &
                           a7(2)*(wf(i,j,l+2) - wf(i,j,l-2)) + &
                           a7(3)*(wf(i,j,l+3) - wf(i,j,l-3)) ) *idz(k)*costeta(i,j,l)
             dpz(i,j,l) = (a7(1)*(pf(i,j,l+1) - pf(i,j,l-1)) + &
                           a7(2)*(pf(i,j,l+2) - pf(i,j,l-2)) + &
                           a7(3)*(pf(i,j,l+3) - pf(i,j,l-3)) ) *idz(k)*costeta(i,j,l)
             drz(i,j,l) = (a7(1)*(rf(i,j,l+1) - rf(i,j,l-1)) + &
                           a7(2)*(rf(i,j,l+2) - rf(i,j,l-2)) + &
                           a7(3)*(rf(i,j,l+3) - rf(i,j,l-3)) ) *idz(k)*costeta(i,j,l)
          enddo
       enddo
    enddo

    k=nz-2
    l=k-nzmngh
    do j=1,ngh
       do i=ndx,nfx
          duz(i,j,l) = ( a42(1)*uf(i,j,l+2)+a42(2)*uf(i,j,l+1) &
                       + a42(3)*uf(i,j,l  )+a42(4)*uf(i,j,l-1) &
                       + a42(5)*uf(i,j,l-2)+a42(6)*uf(i,j,l-3) &
                       + a42(7)*uf(i,j,l-4) ) *idz(k)*costeta(i,j,l)
          dvz(i,j,l) = ( a42(1)*vf(i,j,l+2)+a42(2)*vf(i,j,l+1) &
                       + a42(3)*vf(i,j,l  )+a42(4)*vf(i,j,l-1) &
                       + a42(5)*vf(i,j,l-2)+a42(6)*vf(i,j,l-3) &
                       + a42(7)*vf(i,j,l-4) ) *idz(k)*costeta(i,j,l)
          dwz(i,j,l) = ( a42(1)*wf(i,j,l+2)+a42(2)*wf(i,j,l+1) &
                       + a42(3)*wf(i,j,l  )+a42(4)*wf(i,j,l-1) &
                       + a42(5)*wf(i,j,l-2)+a42(6)*wf(i,j,l-3) &
                       + a42(7)*wf(i,j,l-4) ) *idz(k)*costeta(i,j,l)
          dpz(i,j,l) = ( a42(1)*pf(i,j,l+2)+a42(2)*pf(i,j,l+1) &
                       + a42(3)*pf(i,j,l  )+a42(4)*pf(i,j,l-1) &
                       + a42(5)*pf(i,j,l-2)+a42(6)*pf(i,j,l-3) &
                       + a42(7)*pf(i,j,l-4) ) *idz(k)*costeta(i,j,l)
          drz(i,j,l) = ( a42(1)*rf(i,j,l+2)+a42(2)*rf(i,j,l+1) &
                       + a42(3)*rf(i,j,l  )+a42(4)*rf(i,j,l-1) &
                       + a42(5)*rf(i,j,l-2)+a42(6)*rf(i,j,l-3) &
                       + a42(7)*rf(i,j,l-4) ) *idz(k)*costeta(i,j,l)
       enddo
    enddo

    k=nz-1
    l=k-nzmngh
    do j=1,ngh
       do i=ndx,nfx
          duz(i,j,l) = ( a51(1)*uf(i,j,l+1)+a51(2)*uf(i,j,l  ) &
                       + a51(3)*uf(i,j,l-1)+a51(4)*uf(i,j,l-2) &
                       + a51(5)*uf(i,j,l-3)+a51(6)*uf(i,j,l-4) &
                       + a51(7)*uf(i,j,l-5) ) *idz(k)*costeta(i,j,l)
          dvz(i,j,l) = ( a51(1)*vf(i,j,l+1)+a51(2)*vf(i,j,l  ) &
                       + a51(3)*vf(i,j,l-1)+a51(4)*vf(i,j,l-2) &
                       + a51(5)*vf(i,j,l-3)+a51(6)*vf(i,j,l-4) &
                       + a51(7)*vf(i,j,l-5) ) *idz(k)*costeta(i,j,l)
          dwz(i,j,l) = ( a51(1)*wf(i,j,l+1)+a51(2)*wf(i,j,l  ) &
                       + a51(3)*wf(i,j,l-1)+a51(4)*wf(i,j,l-2) &
                       + a51(5)*wf(i,j,l-3)+a51(6)*wf(i,j,l-4) &
                       + a51(7)*wf(i,j,l-5) ) *idz(k)*costeta(i,j,l)
          dpz(i,j,l) = ( a51(1)*pf(i,j,l+1)+a51(2)*pf(i,j,l  ) &
                       + a51(3)*pf(i,j,l-1)+a51(4)*pf(i,j,l-2) &
                       + a51(5)*pf(i,j,l-3)+a51(6)*pf(i,j,l-4) &
                       + a51(7)*pf(i,j,l-5) ) *idz(k)*costeta(i,j,l)
          drz(i,j,l) = ( a51(1)*rf(i,j,l+1)+a51(2)*rf(i,j,l  ) &
                       + a51(3)*rf(i,j,l-1)+a51(4)*rf(i,j,l-2) &
                       + a51(5)*rf(i,j,l-3)+a51(6)*rf(i,j,l-4) &
                       + a51(7)*rf(i,j,l-5) ) *idz(k)*costeta(i,j,l)
       enddo
    enddo

    k=nz
    l=k-nzmngh
    do j=1,ngh
       do i=ndx,nfx
          duz(i,j,l) = ( a60(1)*uf(i,j,l  )+a60(2)*uf(i,j,l-1) &
                       + a60(3)*uf(i,j,l-2)+a60(4)*uf(i,j,l-3) &
                       + a60(5)*uf(i,j,l-4)+a60(6)*uf(i,j,l-5) &
                       + a60(7)*uf(i,j,l-6) ) *idz(k)*costeta(i,j,l)
          dvz(i,j,l) = ( a60(1)*vf(i,j,l  )+a60(2)*vf(i,j,l-1) &
                       + a60(3)*vf(i,j,l-2)+a60(4)*vf(i,j,l-3) &
                       + a60(5)*vf(i,j,l-4)+a60(6)*vf(i,j,l-5) &
                       + a60(7)*vf(i,j,l-6) ) *idz(k)*costeta(i,j,l)
          dwz(i,j,l) = ( a60(1)*wf(i,j,l  )+a60(2)*wf(i,j,l-1) &
                       + a60(3)*wf(i,j,l-2)+a60(4)*wf(i,j,l-3) &
                       + a60(5)*wf(i,j,l-4)+a60(6)*wf(i,j,l-5) &
                       + a60(7)*wf(i,j,l-6) ) *idz(k)*costeta(i,j,l)
          dpz(i,j,l) = ( a60(1)*pf(i,j,l  )+a60(2)*pf(i,j,l-1) &
                       + a60(3)*pf(i,j,l-2)+a60(4)*pf(i,j,l-3) &
                       + a60(5)*pf(i,j,l-4)+a60(6)*pf(i,j,l-5) &
                       + a60(7)*pf(i,j,l-6) ) *idz(k)*costeta(i,j,l)
          drz(i,j,l) = ( a60(1)*rf(i,j,l  )+a60(2)*rf(i,j,l-1) &
                       + a60(3)*rf(i,j,l-2)+a60(4)*rf(i,j,l-3) &
                       + a60(5)*rf(i,j,l-4)+a60(6)*rf(i,j,l-5) &
                       + a60(7)*rf(i,j,l-6) ) *idz(k)*costeta(i,j,l)
       enddo
    enddo

    ! Compute rt,ut,vt,wt & pt
    ! ========================
    do l=1,ngh
       do j=1,ngh
          do i=ndx,nfx
             dux = ( a7(1)*( uf(i+1,j,l) - uf(i-1,j,l) ) + &
                     a7(2)*( uf(i+2,j,l) - uf(i-2,j,l) ) + &
                     a7(3)*( uf(i+3,j,l) - uf(i-3,j,l) ) ) *idx(i)*sintcosp(i,j,l)
             dvx = ( a7(1)*( vf(i+1,j,l) - vf(i-1,j,l) ) + &
                     a7(2)*( vf(i+2,j,l) - vf(i-2,j,l) ) + &
                     a7(3)*( vf(i+3,j,l) - vf(i-3,j,l) ) ) *idx(i)*sintcosp(i,j,l)
             dwx = ( a7(1)*( wf(i+1,j,l) - wf(i-1,j,l) ) + &
                     a7(2)*( wf(i+2,j,l) - wf(i-2,j,l) ) + &
                     a7(3)*( wf(i+3,j,l) - wf(i-3,j,l) ) ) *idx(i)*sintcosp(i,j,l)
             dpx = ( a7(1)*( pf(i+1,j,l) - pf(i-1,j,l) ) + &
                     a7(2)*( pf(i+2,j,l) - pf(i-2,j,l) ) + &
                     a7(3)*( pf(i+3,j,l) - pf(i-3,j,l) ) ) *idx(i)*sintcosp(i,j,l)
             drx = ( a7(1)*( rf(i+1,j,l) - rf(i-1,j,l) ) + &
                     a7(2)*( rf(i+2,j,l) - rf(i-2,j,l) ) + &
                     a7(3)*( rf(i+3,j,l) - rf(i-3,j,l) ) ) *idx(i)*sintcosp(i,j,l)

             pt(i,j,l) = vg(i,j,l)*(dpx+dpy(i,j,l)+dpz(i,j,l)+pf(i,j,l)*ir(i,j,l))
             ut(i,j,l) = vg(i,j,l)*(dux+duy(i,j,l)+duz(i,j,l)+uf(i,j,l)*ir(i,j,l))
             vt(i,j,l) = vg(i,j,l)*(dvx+dvy(i,j,l)+dvz(i,j,l)+vf(i,j,l)*ir(i,j,l))
             wt(i,j,l) = vg(i,j,l)*(dwx+dwy(i,j,l)+dwz(i,j,l)+wf(i,j,l)*ir(i,j,l))
             rt(i,j,l) = vg(i,j,l)*(drx+dry(i,j,l)+drz(i,j,l)+rf(i,j,l)*ir(i,j,l))
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

  end subroutine bc_TD3d_jmin_kmax

  !===============================================================================
  module subroutine bc_TD3d_jmax_kmin
  !===============================================================================
    !> 3D Tam & Dong's BC: boundary condition at jmax-kmin (edge 2,2,1 /top-front)
    !> - Cartesian version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l
    real(wp) :: cp,av,c2_
    real(wp), dimension(nx1:nx2,-2:ngh,1:ngh+3) :: rf,uf,vf,wf,pf
    real(wp), dimension(nx,1:ngh,1:ngh) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(nx,1:ngh,1:ngh) :: duy,dvy,dwy,dpy,dry
    real(wp), dimension(nx,1:ngh,1:ngh) :: duz,dvz,dwz,dpz,drz
    real(wp) :: dux,dvx,dwx,dpx,drx
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
	           	BC_face(2,2)%U0(i,l,k,4)*costeta(i,l,k)	       &
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
             duy(i,l,k) = ( a7(1)*(uf(i,l+1,k)-uf(i,l-1,k)) &
                          + a7(2)*(uf(i,l+2,k)-uf(i,l-2,k)) &
                          + a7(3)*(uf(i,l+3,k)-uf(i,l-3,k)) ) *idy(j)*sintsinp(i,l,k)
             dvy(i,l,k) = ( a7(1)*(vf(i,l+1,k)-vf(i,l-1,k)) &
                          + a7(2)*(vf(i,l+2,k)-vf(i,l-2,k)) &
                          + a7(3)*(vf(i,l+3,k)-vf(i,l-3,k)) ) *idy(j)*sintsinp(i,l,k)
             dwy(i,l,k) = ( a7(1)*(wf(i,l+1,k)-wf(i,l-1,k)) &
                          + a7(2)*(wf(i,l+2,k)-wf(i,l-2,k)) &
                          + a7(3)*(wf(i,l+3,k)-wf(i,l-3,k)) ) *idy(j)*sintsinp(i,l,k)
             dpy(i,l,k) = ( a7(1)*(pf(i,l+1,k)-pf(i,l-1,k)) &
                          + a7(2)*(pf(i,l+2,k)-pf(i,l-2,k)) &
                          + a7(3)*(pf(i,l+3,k)-pf(i,l-3,k)) ) *idy(j)*sintsinp(i,l,k)
             dry(i,l,k) = ( a7(1)*(rf(i,l+1,k)-rf(i,l-1,k)) &
                          + a7(2)*(rf(i,l+2,k)-rf(i,l-2,k)) &
                          + a7(3)*(rf(i,l+3,k)-rf(i,l-3,k)) ) *idy(j)*sintsinp(i,l,k)
          enddo
       enddo
    enddo

    j=ny-2
    l=j-nymngh
    do k=1,ngh
       do i=ndx,nfx
          duy(i,l,k) = ( a42(1)*uf(i,l+2,k)+a42(2)*uf(i,l+1,k) &
                       + a42(3)*uf(i,l  ,k)+a42(4)*uf(i,l-1,k) &
                       + a42(5)*uf(i,l-2,k)+a42(6)*uf(i,l-3,k) &
                       + a42(7)*uf(i,l-4,k) ) *idy(j)*sintsinp(i,l,k)
          dvy(i,l,k) = ( a42(1)*vf(i,l+2,k)+a42(2)*vf(i,l+1,k) &
                       + a42(3)*vf(i,l  ,k)+a42(4)*vf(i,l-1,k) &
                       + a42(5)*vf(i,l-2,k)+a42(6)*vf(i,l-3,k) &
                       + a42(7)*vf(i,l-4,k) ) *idy(j)*sintsinp(i,l,k)
          dwy(i,l,k) = ( a42(1)*wf(i,l+2,k)+a42(2)*wf(i,l+1,k) &
                       + a42(3)*wf(i,l  ,k)+a42(4)*wf(i,l-1,k) &
                       + a42(5)*wf(i,l-2,k)+a42(6)*wf(i,l-3,k) &
                       + a42(7)*wf(i,l-4,k) ) *idy(j)*sintsinp(i,l,k)
          dpy(i,l,k) = ( a42(1)*pf(i,l+2,k)+a42(2)*pf(i,l+1,k) &
                       + a42(3)*pf(i,l  ,k)+a42(4)*pf(i,l-1,k) &
                       + a42(5)*pf(i,l-2,k)+a42(6)*pf(i,l-3,k) &
                       + a42(7)*pf(i,l-4,k) ) *idy(j)*sintsinp(i,l,k)
          dry(i,l,k) = ( a42(1)*rf(i,l+2,k)+a42(2)*rf(i,l+1,k) &
                       + a42(3)*rf(i,l  ,k)+a42(4)*rf(i,l-1,k) &
                       + a42(5)*rf(i,l-2,k)+a42(6)*rf(i,l-3,k) &
                       + a42(7)*rf(i,l-4,k) ) *idy(j)*sintsinp(i,l,k)
       enddo
    enddo

    j=ny-1
    l=j-nymngh
    do k=1,ngh
       do i=ndx,nfx
          duy(i,l,k) = ( a51(1)*uf(i,l+1,k)+a51(2)*uf(i,l  ,k) &
                       + a51(3)*uf(i,l-1,k)+a51(4)*uf(i,l-2,k) &
                       + a51(5)*uf(i,l-3,k)+a51(6)*uf(i,l-4,k) &
                       + a51(7)*uf(i,l-5,k) ) *idy(j)*sintsinp(i,l,k)
          dvy(i,l,k) = ( a51(1)*vf(i,l+1,k)+a51(2)*vf(i,l  ,k) &
                       + a51(3)*vf(i,l-1,k)+a51(4)*vf(i,l-2,k) &
                       + a51(5)*vf(i,l-3,k)+a51(6)*vf(i,l-4,k) &
                       + a51(7)*vf(i,l-5,k) ) *idy(j)*sintsinp(i,l,k)
          dwy(i,l,k) = ( a51(1)*wf(i,l+1,k)+a51(2)*wf(i,l  ,k) &
                       + a51(3)*wf(i,l-1,k)+a51(4)*wf(i,l-2,k) &
                       + a51(5)*wf(i,l-3,k)+a51(6)*wf(i,l-4,k) &
                       + a51(7)*wf(i,l-5,k) ) *idy(j)*sintsinp(i,l,k)
          dpy(i,l,k) = ( a51(1)*pf(i,l+1,k)+a51(2)*pf(i,l  ,k) &
                       + a51(3)*pf(i,l-1,k)+a51(4)*pf(i,l-2,k) &
                       + a51(5)*pf(i,l-3,k)+a51(6)*pf(i,l-4,k) &
                       + a51(7)*pf(i,l-5,k) ) *idy(j)*sintsinp(i,l,k)
          dry(i,l,k) = ( a51(1)*rf(i,l+1,k)+a51(2)*rf(i,l  ,k) &
                       + a51(3)*rf(i,l-1,k)+a51(4)*rf(i,l-2,k) &
                       + a51(5)*rf(i,l-3,k)+a51(6)*rf(i,l-4,k) &
                       + a51(7)*rf(i,l-5,k) ) *idy(j)*sintsinp(i,l,k)
       enddo
    enddo

    j=ny
    l=j-nymngh
    do k=1,ngh
       do i=ndx,nfx
          duy(i,l,k) = ( a60(1)*uf(i,l  ,k)+a60(2)*uf(i,l-1,k) &
                       + a60(3)*uf(i,l-2,k)+a60(4)*uf(i,l-3,k) &
                       + a60(5)*uf(i,l-4,k)+a60(6)*uf(i,l-5,k) &
                       + a60(7)*uf(i,l-6,k) ) *idy(j)*sintsinp(i,l,k)
          dvy(i,l,k) = ( a60(1)*vf(i,l  ,k)+a60(2)*vf(i,l-1,k) &
                       + a60(3)*vf(i,l-2,k)+a60(4)*vf(i,l-3,k) &
                       + a60(5)*vf(i,l-4,k)+a60(6)*vf(i,l-5,k) &
                       + a60(7)*vf(i,l-6,k) ) *idy(j)*sintsinp(i,l,k)
          dwy(i,l,k) = ( a60(1)*wf(i,l  ,k)+a60(2)*wf(i,l-1,k) &
                       + a60(3)*wf(i,l-2,k)+a60(4)*wf(i,l-3,k) &
                       + a60(5)*wf(i,l-4,k)+a60(6)*wf(i,l-5,k) &
                       + a60(7)*wf(i,l-6,k) ) *idy(j)*sintsinp(i,l,k)
          dpy(i,l,k) = ( a60(1)*pf(i,l  ,k)+a60(2)*pf(i,l-1,k) &
                       + a60(3)*pf(i,l-2,k)+a60(4)*pf(i,l-3,k) &
                       + a60(5)*pf(i,l-4,k)+a60(6)*pf(i,l-5,k) &
                       + a60(7)*pf(i,l-6,k) ) *idy(j)*sintsinp(i,l,k)
          dry(i,l,k) = ( a60(1)*rf(i,l  ,k)+a60(2)*rf(i,l-1,k) &
                       + a60(3)*rf(i,l-2,k)+a60(4)*rf(i,l-3,k) &
                       + a60(5)*rf(i,l-4,k)+a60(6)*rf(i,l-5,k) &
                       + a60(7)*rf(i,l-6,k) ) *idy(j)*sintsinp(i,l,k)
       enddo
    enddo

    ! Non-centered derivatives in z-direction *cos(teta)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    k=1
    do l=1,ngh
       do i=ndx,nfx
          duz(i,l,k) = ( a06(1)*uf(i,l,1)+a06(2)*uf(i,l,2) &
                       + a06(3)*uf(i,l,3)+a06(4)*uf(i,l,4) &
                       + a06(5)*uf(i,l,5)+a06(6)*uf(i,l,6) &
                       + a06(7)*uf(i,l,7) ) *idz(k)*costeta(i,l,k)
          dvz(i,l,k) = ( a06(1)*vf(i,l,1)+a06(2)*vf(i,l,2) &
                       + a06(3)*vf(i,l,3)+a06(4)*vf(i,l,4) &
                       + a06(5)*vf(i,l,5)+a06(6)*vf(i,l,6) &
                       + a06(7)*vf(i,l,7) ) *idz(k)*costeta(i,l,k)
          dwz(i,l,k) = ( a06(1)*wf(i,l,1)+a06(2)*wf(i,l,2) &
                       + a06(3)*wf(i,l,3)+a06(4)*wf(i,l,4) &
                       + a06(5)*wf(i,l,5)+a06(6)*wf(i,l,6) &
                       + a06(7)*wf(i,l,7) ) *idz(k)*costeta(i,l,k)
          dpz(i,l,k) = ( a06(1)*pf(i,l,1)+a06(2)*pf(i,l,2) &
                       + a06(3)*pf(i,l,3)+a06(4)*pf(i,l,4) &
                       + a06(5)*pf(i,l,5)+a06(6)*pf(i,l,6) &
                       + a06(7)*pf(i,l,7) ) *idz(k)*costeta(i,l,k)
          drz(i,l,k) = ( a06(1)*rf(i,l,1)+a06(2)*rf(i,l,2) &
                       + a06(3)*rf(i,l,3)+a06(4)*rf(i,l,4) &
                       + a06(5)*rf(i,l,5)+a06(6)*rf(i,l,6) &
                       + a06(7)*rf(i,l,7) ) *idz(k)*costeta(i,l,k)
       enddo
    enddo

    k=2
    do l=1,ngh
       do i=ndx,nfx
          duz(i,l,k) = ( a15(1)*uf(i,l,1)+a15(2)*uf(i,l,2) &
                       + a15(3)*uf(i,l,3)+a15(4)*uf(i,l,4) &
                       + a15(5)*uf(i,l,5)+a15(6)*uf(i,l,6) &
                       + a15(7)*uf(i,l,7) ) *idz(k)*costeta(i,l,k)
          dvz(i,l,k) = ( a15(1)*vf(i,l,1)+a15(2)*vf(i,l,2) &
                       + a15(3)*vf(i,l,3)+a15(4)*vf(i,l,4) &
                       + a15(5)*vf(i,l,5)+a15(6)*vf(i,l,6) &
                       + a15(7)*vf(i,l,7) ) *idz(k)*costeta(i,l,k)
          dwz(i,l,k) = ( a15(1)*wf(i,l,1)+a15(2)*wf(i,l,2) &
                       + a15(3)*wf(i,l,3)+a15(4)*wf(i,l,4) &
                       + a15(5)*wf(i,l,5)+a15(6)*wf(i,l,6) &
                       + a15(7)*wf(i,l,7) ) *idz(k)*costeta(i,l,k)
          dpz(i,l,k) = ( a15(1)*pf(i,l,1)+a15(2)*pf(i,l,2) &
                       + a15(3)*pf(i,l,3)+a15(4)*pf(i,l,4) &
                       + a15(5)*pf(i,l,5)+a15(6)*pf(i,l,6) &
                       + a15(7)*pf(i,l,7) ) *idz(k)*costeta(i,l,k)
          drz(i,l,k) = ( a15(1)*rf(i,l,1)+a15(2)*rf(i,l,2) &
                       + a15(3)*rf(i,l,3)+a15(4)*rf(i,l,4) &
                       + a15(5)*rf(i,l,5)+a15(6)*rf(i,l,6) &
                       + a15(7)*rf(i,l,7) ) *idz(k)*costeta(i,l,k)
       enddo
    enddo

    k=3
    do l=1,ngh
       do i=ndx,nfx
          duz(i,l,k) = ( a24(1)*uf(i,l,1)+a24(2)*uf(i,l,2) &
                       + a24(3)*uf(i,l,3)+a24(4)*uf(i,l,4) &
                       + a24(5)*uf(i,l,5)+a24(6)*uf(i,l,6) &
                       + a24(7)*uf(i,l,7) ) *idz(k)*costeta(i,l,k)
          dvz(i,l,k) = ( a24(1)*vf(i,l,1)+a24(2)*vf(i,l,2) &
                       + a24(3)*vf(i,l,3)+a24(4)*vf(i,l,4) &
                       + a24(5)*vf(i,l,5)+a24(6)*vf(i,l,6) &
                       + a24(7)*vf(i,l,7) ) *idz(k)*costeta(i,l,k)
          dwz(i,l,k) = ( a24(1)*wf(i,l,1)+a24(2)*wf(i,l,2) &
                       + a24(3)*wf(i,l,3)+a24(4)*wf(i,l,4) &
                       + a24(5)*wf(i,l,5)+a24(6)*wf(i,l,6) &
                       + a24(7)*wf(i,l,7) ) *idz(k)*costeta(i,l,k)
          dpz(i,l,k) = ( a24(1)*pf(i,l,1)+a24(2)*pf(i,l,2) &
                       + a24(3)*pf(i,l,3)+a24(4)*pf(i,l,4) &
                       + a24(5)*pf(i,l,5)+a24(6)*pf(i,l,6) &
                       + a24(7)*pf(i,l,7) ) *idz(k)*costeta(i,l,k)
          drz(i,l,k) = ( a24(1)*rf(i,l,1)+a24(2)*rf(i,l,2) &
                       + a24(3)*rf(i,l,3)+a24(4)*rf(i,l,4) &
                       + a24(5)*rf(i,l,5)+a24(6)*rf(i,l,6) &
                       + a24(7)*rf(i,l,7) ) *idz(k)*costeta(i,l,k)
       enddo
    enddo

    do k=4,ngh
       do l=1,ngh
          do i=ndx,nfx
             duz(i,l,k) = (a7(1)*(uf(i,l,k+1) - uf(i,l,k-1)) + &
                           a7(2)*(uf(i,l,k+2) - uf(i,l,k-2)) + &
                           a7(3)*(uf(i,l,k+3) - uf(i,l,k-3)) ) *idz(k)*costeta(i,l,k)
             dvz(i,l,k) = (a7(1)*(vf(i,l,k+1) - vf(i,l,k-1)) + &
                           a7(2)*(vf(i,l,k+2) - vf(i,l,k-2)) + &
                           a7(3)*(vf(i,l,k+3) - vf(i,l,k-3)) ) *idz(k)*costeta(i,l,k)
             dwz(i,l,k) = (a7(1)*(wf(i,l,k+1) - wf(i,l,k-1)) + &
                           a7(2)*(wf(i,l,k+2) - wf(i,l,k-2)) + &
                           a7(3)*(wf(i,l,k+3) - wf(i,l,k-3)) ) *idz(k)*costeta(i,l,k)
             dpz(i,l,k) = (a7(1)*(pf(i,l,k+1) - pf(i,l,k-1)) + &
                           a7(2)*(pf(i,l,k+2) - pf(i,l,k-2)) + &
                           a7(3)*(pf(i,l,k+3) - pf(i,l,k-3)) ) *idz(k)*costeta(i,l,k)
             drz(i,l,k) = (a7(1)*(rf(i,l,k+1) - rf(i,l,k-1)) + &
                           a7(2)*(rf(i,l,k+2) - rf(i,l,k-2)) + &
                           a7(3)*(rf(i,l,k+3) - rf(i,l,k-3)) ) *idz(k)*costeta(i,l,k)
          enddo
       enddo
    enddo

    ! Compute rt,ut,vt,wt & pt
    ! ========================
    do k=1,ngh
       do l=1,ngh
          do i=ndx,nfx
             dux = ( a7(1)*( uf(i+1,l,k) - uf(i-1,l,k) ) + &
                     a7(2)*( uf(i+2,l,k) - uf(i-2,l,k) ) + &
                     a7(3)*( uf(i+3,l,k) - uf(i-3,l,k) ) ) *idx(i)*sintcosp(i,l,k)
             dvx = ( a7(1)*( vf(i+1,l,k) - vf(i-1,l,k) ) + &
                     a7(2)*( vf(i+2,l,k) - vf(i-2,l,k) ) + &
                     a7(3)*( vf(i+3,l,k) - vf(i-3,l,k) ) ) *idx(i)*sintcosp(i,l,k)
             dwx = ( a7(1)*( wf(i+1,l,k) - wf(i-1,l,k) ) + &
                     a7(2)*( wf(i+2,l,k) - wf(i-2,l,k) ) + &
                     a7(3)*( wf(i+3,l,k) - wf(i-3,l,k) ) ) *idx(i)*sintcosp(i,l,k)
             dpx = ( a7(1)*( pf(i+1,l,k) - pf(i-1,l,k) ) + &
                     a7(2)*( pf(i+2,l,k) - pf(i-2,l,k) ) + &
                     a7(3)*( pf(i+3,l,k) - pf(i-3,l,k) ) ) *idx(i)*sintcosp(i,l,k)
             drx = ( a7(1)*( rf(i+1,l,k) - rf(i-1,l,k) ) + &
                     a7(2)*( rf(i+2,l,k) - rf(i-2,l,k) ) + &
                     a7(3)*( rf(i+3,l,k) - rf(i-3,l,k) ) ) *idx(i)*sintcosp(i,l,k)

             pt(i,l,k) = vg(i,l,k)*(dpx+dpy(i,l,k)+dpz(i,l,k)+pf(i,l,k)*ir(i,l,k))
             ut(i,l,k) = vg(i,l,k)*(dux+duy(i,l,k)+duz(i,l,k)+uf(i,l,k)*ir(i,l,k))
             vt(i,l,k) = vg(i,l,k)*(dvx+dvy(i,l,k)+dvz(i,l,k)+vf(i,l,k)*ir(i,l,k))
             wt(i,l,k) = vg(i,l,k)*(dwx+dwy(i,l,k)+dwz(i,l,k)+wf(i,l,k)*ir(i,l,k))
             rt(i,l,k) = vg(i,l,k)*(drx+dry(i,l,k)+drz(i,l,k)+rf(i,l,k)*ir(i,l,k))
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

  end subroutine bc_TD3d_jmax_kmin

  !===============================================================================
  module subroutine bc_TD3d_jmax_kmax
  !===============================================================================
    !> 3D Tam & Dong's BC: boundary condition at jmax-kmax (edge 2,2,2 /top-back)
    !> - Cartesian version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l,m
    real(wp) :: cp,av,c2_
    real(wp), dimension(nx1:nx2,-2:ngh,-2:ngh) :: rf,uf,vf,wf,pf
    real(wp), dimension(nx,1:ngh,1:ngh) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(nx,1:ngh,1:ngh) :: duy,dvy,dwy,dpy,dry
    real(wp), dimension(nx,1:ngh,1:ngh) :: duz,dvz,dwz,dpz,drz
    real(wp) :: dux,dvx,dwx,dpx,drx
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
	           	BC_face(2,2)%U0(i,l,k,4)*costeta(i,l,m)	       &
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
             duy(i,l,m) = ( a7(1)*(uf(i,l+1,m)-uf(i,l-1,m)) &
                          + a7(2)*(uf(i,l+2,m)-uf(i,l-2,m)) &
                          + a7(3)*(uf(i,l+3,m)-uf(i,l-3,m)) ) *idy(j)*sintsinp(i,l,m)
             dvy(i,l,m) = ( a7(1)*(vf(i,l+1,m)-vf(i,l-1,m)) &
                          + a7(2)*(vf(i,l+2,m)-vf(i,l-2,m)) &
                          + a7(3)*(vf(i,l+3,m)-vf(i,l-3,m)) ) *idy(j)*sintsinp(i,l,m)
             dwy(i,l,m) = ( a7(1)*(wf(i,l+1,m)-wf(i,l-1,m)) &
                          + a7(2)*(wf(i,l+2,m)-wf(i,l-2,m)) &
                          + a7(3)*(wf(i,l+3,m)-wf(i,l-3,m)) ) *idy(j)*sintsinp(i,l,m)
             dpy(i,l,m) = ( a7(1)*(pf(i,l+1,m)-pf(i,l-1,m)) &
                          + a7(2)*(pf(i,l+2,m)-pf(i,l-2,m)) &
                          + a7(3)*(pf(i,l+3,m)-pf(i,l-3,m)) ) *idy(j)*sintsinp(i,l,m)
             dry(i,l,m) = ( a7(1)*(rf(i,l+1,m)-rf(i,l-1,m)) &
                          + a7(2)*(rf(i,l+2,m)-rf(i,l-2,m)) &
                          + a7(3)*(rf(i,l+3,m)-rf(i,l-3,m)) ) *idy(j)*sintsinp(i,l,m)
          enddo
       enddo
    enddo

    j=ny-2
    l=j-nymngh
    do m=1,ngh
       do i=ndx,nfx
          duy(i,l,m) = ( a42(1)*uf(i,l+2,m)+a42(2)*uf(i,l+1,m) &
                       + a42(3)*uf(i,l  ,m)+a42(4)*uf(i,l-1,m) &
                       + a42(5)*uf(i,l-2,m)+a42(6)*uf(i,l-3,m) &
                       + a42(7)*uf(i,l-4,m) ) *idy(j)*sintsinp(i,l,m)
          dvy(i,l,m) = ( a42(1)*vf(i,l+2,m)+a42(2)*vf(i,l+1,m) &
                       + a42(3)*vf(i,l  ,m)+a42(4)*vf(i,l-1,m) &
                       + a42(5)*vf(i,l-2,m)+a42(6)*vf(i,l-3,m) &
                       + a42(7)*vf(i,l-4,m) ) *idy(j)*sintsinp(i,l,m)
          dwy(i,l,m) = ( a42(1)*wf(i,l+2,m)+a42(2)*wf(i,l+1,m) &
                       + a42(3)*wf(i,l  ,m)+a42(4)*wf(i,l-1,m) &
                       + a42(5)*wf(i,l-2,m)+a42(6)*wf(i,l-3,m) &
                       + a42(7)*wf(i,l-4,m) ) *idy(j)*sintsinp(i,l,m)
          dpy(i,l,m) = ( a42(1)*pf(i,l+2,m)+a42(2)*pf(i,l+1,m) &
                       + a42(3)*pf(i,l  ,m)+a42(4)*pf(i,l-1,m) &
                       + a42(5)*pf(i,l-2,m)+a42(6)*pf(i,l-3,m) &
                       + a42(7)*pf(i,l-4,m) ) *idy(j)*sintsinp(i,l,m)
          dry(i,l,m) = ( a42(1)*rf(i,l+2,m)+a42(2)*rf(i,l+1,m) &
                       + a42(3)*rf(i,l  ,m)+a42(4)*rf(i,l-1,m) &
                       + a42(5)*rf(i,l-2,m)+a42(6)*rf(i,l-3,m) &
                       + a42(7)*rf(i,l-4,m) ) *idy(j)*sintsinp(i,l,m)
       enddo
    enddo
    
    j=ny-1
    l=j-nymngh
    do m=1,ngh
       do i=ndx,nfx
          duy(i,l,m) = ( a51(1)*uf(i,l+1,m)+a51(2)*uf(i,l  ,m) &
                       + a51(3)*uf(i,l-1,m)+a51(4)*uf(i,l-2,m) &
                       + a51(5)*uf(i,l-3,m)+a51(6)*uf(i,l-4,m) &
                       + a51(7)*uf(i,l-5,m) ) *idy(j)*sintsinp(i,l,m)
          dvy(i,l,m) = ( a51(1)*vf(i,l+1,m)+a51(2)*vf(i,l  ,m) &
                       + a51(3)*vf(i,l-1,m)+a51(4)*vf(i,l-2,m) &
                       + a51(5)*vf(i,l-3,m)+a51(6)*vf(i,l-4,m) &
                       + a51(7)*vf(i,l-5,m) ) *idy(j)*sintsinp(i,l,m)
          dwy(i,l,m) = ( a51(1)*wf(i,l+1,m)+a51(2)*wf(i,l  ,m) &
                       + a51(3)*wf(i,l-1,m)+a51(4)*wf(i,l-2,m) &
                       + a51(5)*wf(i,l-3,m)+a51(6)*wf(i,l-4,m) &
                       + a51(7)*wf(i,l-5,m) ) *idy(j)*sintsinp(i,l,m)
          dpy(i,l,m) = ( a51(1)*pf(i,l+1,m)+a51(2)*pf(i,l  ,m) &
                       + a51(3)*pf(i,l-1,m)+a51(4)*pf(i,l-2,m) &
                       + a51(5)*pf(i,l-3,m)+a51(6)*pf(i,l-4,m) &
                       + a51(7)*pf(i,l-5,m) ) *idy(j)*sintsinp(i,l,m)
          dry(i,l,m) = ( a51(1)*rf(i,l+1,m)+a51(2)*rf(i,l  ,m) &
                       + a51(3)*rf(i,l-1,m)+a51(4)*rf(i,l-2,m) &
                       + a51(5)*rf(i,l-3,m)+a51(6)*rf(i,l-4,m) &
                       + a51(7)*rf(i,l-5,m) ) *idy(j)*sintsinp(i,l,m)
       enddo
    enddo
    
    j=ny
    l=j-nymngh
    do m=1,ngh
       do i=ndx,nfx
          duy(i,l,m) = ( a60(1)*uf(i,l  ,m)+a60(2)*uf(i,l-1,m) &
                       + a60(3)*uf(i,l-2,m)+a60(4)*uf(i,l-3,m) &
                       + a60(5)*uf(i,l-4,m)+a60(6)*uf(i,l-5,m) &
                       + a60(7)*uf(i,l-6,m) ) *idy(j)*sintsinp(i,l,m)
          dvy(i,l,m) = ( a60(1)*vf(i,l  ,m)+a60(2)*vf(i,l-1,m) &
                       + a60(3)*vf(i,l-2,m)+a60(4)*vf(i,l-3,m) &
                       + a60(5)*vf(i,l-4,m)+a60(6)*vf(i,l-5,m) &
                       + a60(7)*vf(i,l-6,m) ) *idy(j)*sintsinp(i,l,m)
          dwy(i,l,m) = ( a60(1)*wf(i,l  ,m)+a60(2)*wf(i,l-1,m) &
                       + a60(3)*wf(i,l-2,m)+a60(4)*wf(i,l-3,m) &
                       + a60(5)*wf(i,l-4,m)+a60(6)*wf(i,l-5,m) &
                       + a60(7)*wf(i,l-6,m) ) *idy(j)*sintsinp(i,l,m)
          dpy(i,l,m) = ( a60(1)*pf(i,l  ,m)+a60(2)*pf(i,l-1,m) &
                       + a60(3)*pf(i,l-2,m)+a60(4)*pf(i,l-3,m) &
                       + a60(5)*pf(i,l-4,m)+a60(6)*pf(i,l-5,m) &
                       + a60(7)*pf(i,l-6,m) ) *idy(j)*sintsinp(i,l,m)
          dry(i,l,m) = ( a60(1)*rf(i,l  ,m)+a60(2)*rf(i,l-1,m) &
                       + a60(3)*rf(i,l-2,m)+a60(4)*rf(i,l-3,m) &
                       + a60(5)*rf(i,l-4,m)+a60(6)*rf(i,l-5,m) &
                       + a60(7)*rf(i,l-6,m) ) *idy(j)*sintsinp(i,l,m)
       enddo
    enddo
              
    ! Non-centered derivatives in z-direction *cos(teta)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    do k=nzmnghp1,nz-3
       m=k-nzmngh
       do l=1,ngh
          do i=ndx,nfx
             duz(i,l,m) = (a7(1)*(uf(i,l,m+1) - uf(i,l,m-1)) + &
                           a7(2)*(uf(i,l,m+2) - uf(i,l,m-2)) + &
                           a7(3)*(uf(i,l,m+3) - uf(i,l,m-3)) ) *idz(k)*costeta(i,l,m)
             dvz(i,l,m) = (a7(1)*(vf(i,l,m+1) - vf(i,l,m-1)) + &
                           a7(2)*(vf(i,l,m+2) - vf(i,l,m-2)) + &
                           a7(3)*(vf(i,l,m+3) - vf(i,l,m-3)) ) *idz(k)*costeta(i,l,m)
             dwz(i,l,m) = (a7(1)*(wf(i,l,m+1) - wf(i,l,m-1)) + &
                           a7(2)*(wf(i,l,m+2) - wf(i,l,m-2)) + &
                           a7(3)*(wf(i,l,m+3) - wf(i,l,m-3)) ) *idz(k)*costeta(i,l,m)
             dpz(i,l,m) = (a7(1)*(pf(i,l,m+1) - pf(i,l,m-1)) + &
                           a7(2)*(pf(i,l,m+2) - pf(i,l,m-2)) + &
                           a7(3)*(pf(i,l,m+3) - pf(i,l,m-3)) ) *idz(k)*costeta(i,l,m)
             drz(i,l,m) = (a7(1)*(rf(i,l,m+1) - rf(i,l,m-1)) + &
                           a7(2)*(rf(i,l,m+2) - rf(i,l,m-2)) + &
                           a7(3)*(rf(i,l,m+3) - rf(i,l,m-3)) ) *idz(k)*costeta(i,l,m)
          enddo
       enddo
    enddo

    k=nz-2
    m=k-nzmngh
    do l=1,ngh
       do i=ndx,nfx
          duz(i,l,m) = ( a42(1)*uf(i,l,m+2)+a42(2)*uf(i,l,m+1) &
                       + a42(3)*uf(i,l,m  )+a42(4)*uf(i,l,m-1) &
                       + a42(5)*uf(i,l,m-2)+a42(6)*uf(i,l,m-3) &
                       + a42(7)*uf(i,l,m-4) ) *idz(k)*costeta(i,l,m)
          dvz(i,l,m) = ( a42(1)*vf(i,l,m+2)+a42(2)*vf(i,l,m+1) &
                       + a42(3)*vf(i,l,m  )+a42(4)*vf(i,l,m-1) &
                       + a42(5)*vf(i,l,m-2)+a42(6)*vf(i,l,m-3) &
                       + a42(7)*vf(i,l,m-4) ) *idz(k)*costeta(i,l,m)
          dwz(i,l,m) = ( a42(1)*wf(i,l,m+2)+a42(2)*wf(i,l,m+1) &
                       + a42(3)*wf(i,l,m  )+a42(4)*wf(i,l,m-1) &
                       + a42(5)*wf(i,l,m-2)+a42(6)*wf(i,l,m-3) &
                       + a42(7)*wf(i,l,m-4) ) *idz(k)*costeta(i,l,m)
          dpz(i,l,m) = ( a42(1)*pf(i,l,m+2)+a42(2)*pf(i,l,m+1) &
                       + a42(3)*pf(i,l,m  )+a42(4)*pf(i,l,m-1) &
                       + a42(5)*pf(i,l,m-2)+a42(6)*pf(i,l,m-3) &
                       + a42(7)*pf(i,l,m-4) ) *idz(k)*costeta(i,l,m)
          drz(i,l,m) = ( a42(1)*rf(i,l,m+2)+a42(2)*rf(i,l,m+1) &
                       + a42(3)*rf(i,l,m  )+a42(4)*rf(i,l,m-1) &
                       + a42(5)*rf(i,l,m-2)+a42(6)*rf(i,l,m-3) &
                       + a42(7)*rf(i,l,m-4) ) *idz(k)*costeta(i,l,m)
       enddo
    enddo

    k=nz-1
    m=k-nzmngh
    do l=1,ngh
       do i=ndx,nfx
          duz(i,l,m) = ( a51(1)*uf(i,l,m+1)+a51(2)*uf(i,l,m  ) &
                       + a51(3)*uf(i,l,m-1)+a51(4)*uf(i,l,m-2) &
                       + a51(5)*uf(i,l,m-3)+a51(6)*uf(i,l,m-4) &
                       + a51(7)*uf(i,l,m-5) ) *idz(k)*costeta(i,l,m)
          dvz(i,l,m) = ( a51(1)*vf(i,l,m+1)+a51(2)*vf(i,l,m  ) &
                       + a51(3)*vf(i,l,m-1)+a51(4)*vf(i,l,m-2) &
                       + a51(5)*vf(i,l,m-3)+a51(6)*vf(i,l,m-4) &
                       + a51(7)*vf(i,l,m-5) ) *idz(k)*costeta(i,l,m)
          dwz(i,l,m) = ( a51(1)*wf(i,l,m+1)+a51(2)*wf(i,l,m  ) &
                       + a51(3)*wf(i,l,m-1)+a51(4)*wf(i,l,m-2) &
                       + a51(5)*wf(i,l,m-3)+a51(6)*wf(i,l,m-4) &
                       + a51(7)*wf(i,l,m-5) ) *idz(k)*costeta(i,l,m)
          dpz(i,l,m) = ( a51(1)*pf(i,l,m+1)+a51(2)*pf(i,l,m  ) &
                       + a51(3)*pf(i,l,m-1)+a51(4)*pf(i,l,m-2) &
                       + a51(5)*pf(i,l,m-3)+a51(6)*pf(i,l,m-4) &
                       + a51(7)*pf(i,l,m-5) ) *idz(k)*costeta(i,l,m)
          drz(i,l,m) = ( a51(1)*rf(i,l,m+1)+a51(2)*rf(i,l,m  ) &
                       + a51(3)*rf(i,l,m-1)+a51(4)*rf(i,l,m-2) &
                       + a51(5)*rf(i,l,m-3)+a51(6)*rf(i,l,m-4) &
                       + a51(7)*rf(i,l,m-5) ) *idz(k)*costeta(i,l,m)
       enddo
    enddo

    k=nz
    m=k-nzmngh
    do l=1,ngh
       do i=ndx,nfx
          duz(i,l,m) = ( a60(1)*uf(i,l,m  )+a60(2)*uf(i,l,m-1) &
                       + a60(3)*uf(i,l,m-2)+a60(4)*uf(i,l,m-3) &
                       + a60(5)*uf(i,l,m-4)+a60(6)*uf(i,l,m-5) &
                       + a60(7)*uf(i,l,m-6) ) *idz(k)*costeta(i,l,m)
          dvz(i,l,m) = ( a60(1)*vf(i,l,m  )+a60(2)*vf(i,l,m-1) &
                       + a60(3)*vf(i,l,m-2)+a60(4)*vf(i,l,m-3) &
                       + a60(5)*vf(i,l,m-4)+a60(6)*vf(i,l,m-5) &
                       + a60(7)*vf(i,l,m-6) ) *idz(k)*costeta(i,l,m)
          dwz(i,l,m) = ( a60(1)*wf(i,l,m  )+a60(2)*wf(i,l,m-1) &
                       + a60(3)*wf(i,l,m-2)+a60(4)*wf(i,l,m-3) &
                       + a60(5)*wf(i,l,m-4)+a60(6)*wf(i,l,m-5) &
                       + a60(7)*wf(i,l,m-6) ) *idz(k)*costeta(i,l,m)
          dpz(i,l,m) = ( a60(1)*pf(i,l,m  )+a60(2)*pf(i,l,m-1) &
                       + a60(3)*pf(i,l,m-2)+a60(4)*pf(i,l,m-3) &
                       + a60(5)*pf(i,l,m-4)+a60(6)*pf(i,l,m-5) &
                       + a60(7)*pf(i,l,m-6) ) *idz(k)*costeta(i,l,m)
          drz(i,l,m) = ( a60(1)*rf(i,l,m  )+a60(2)*rf(i,l,m-1) &
                       + a60(3)*rf(i,l,m-2)+a60(4)*rf(i,l,m-3) &
                       + a60(5)*rf(i,l,m-4)+a60(6)*rf(i,l,m-5) &
                       + a60(7)*rf(i,l,m-6) ) *idz(k)*costeta(i,l,m)
       enddo
    enddo

    ! Compute rt,ut,vt,wt & pt
    ! ========================
    do l=1,ngh
       do m=1,ngh
          do i=ndx,nfx
             dux = ( a7(1)*( uf(i+1,l,m) - uf(i-1,l,m) ) + &
                     a7(2)*( uf(i+2,l,m) - uf(i-2,l,m) ) + &
                     a7(3)*( uf(i+3,l,m) - uf(i-3,l,m) ) ) *idx(i)*sintcosp(i,l,m)
             dvx = ( a7(1)*( vf(i+1,l,m) - vf(i-1,l,m) ) + &
                     a7(2)*( vf(i+2,l,m) - vf(i-2,l,m) ) + &
                     a7(3)*( vf(i+3,l,m) - vf(i-3,l,m) ) ) *idx(i)*sintcosp(i,l,m)
             dwx = ( a7(1)*( wf(i+1,l,m) - wf(i-1,l,m) ) + &
                     a7(2)*( wf(i+2,l,m) - wf(i-2,l,m) ) + &
                     a7(3)*( wf(i+3,l,m) - wf(i-3,l,m) ) ) *idx(i)*sintcosp(i,l,m)
             dpx = ( a7(1)*( pf(i+1,l,m) - pf(i-1,l,m) ) + &
                     a7(2)*( pf(i+2,l,m) - pf(i-2,l,m) ) + &
                     a7(3)*( pf(i+3,l,m) - pf(i-3,l,m) ) ) *idx(i)*sintcosp(i,l,m)
             drx = ( a7(1)*( rf(i+1,l,m) - rf(i-1,l,m) ) + &
                     a7(2)*( rf(i+2,l,m) - rf(i-2,l,m) ) + &
                     a7(3)*( rf(i+3,l,m) - rf(i-3,l,m) ) ) *idx(i)*sintcosp(i,l,m)

             pt(i,l,m) = vg(i,l,m)*(dpx+dpy(i,l,m)+dpz(i,l,m)+pf(i,l,m)*ir(i,l,m))
             ut(i,l,m) = vg(i,l,m)*(dux+duy(i,l,m)+duz(i,l,m)+uf(i,l,m)*ir(i,l,m))
             vt(i,l,m) = vg(i,l,m)*(dvx+dvy(i,l,m)+dvz(i,l,m)+vf(i,l,m)*ir(i,l,m))
             wt(i,l,m) = vg(i,l,m)*(dwx+dwy(i,l,m)+dwz(i,l,m)+wf(i,l,m)*ir(i,l,m))
             rt(i,l,m) = vg(i,l,m)*(drx+dry(i,l,m)+drz(i,l,m)+rf(i,l,m)*ir(i,l,m))
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

  end subroutine bc_TD3d_jmax_kmax

end submodule smod_TamDong3d_edgex
