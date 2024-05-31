!===============================================================================
submodule (mod_TamDong3d) smod_TamDong3d_edgey
!===============================================================================
  !> author: XG
  !> date: February 2020 - modif January 2022
  !> Non-reflecting boundary conditions of Tam and Dong
  !> - 3D (periodic) Cartesian version - edges along y
!=============================================================================== 

contains

  !===============================================================================
  module subroutine bc_TD3d_kmin_imin
  !===============================================================================
    !> 3D Tam & Dong's BC: boundary condition at kmin-imin (edge 3,1,1 /front-left)
    !> - Cartesian version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: cp,av,c2_
    real(wp), dimension(1:ngh+3,ny1:ny2,1:ngh+3) :: rf,uf,vf,wf,pf
    real(wp), dimension(1:ngh,ny,1:ngh) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(1:ngh,ny,1:ngh) :: dux,dvx,dwx,dpx,drx
    real(wp), dimension(1:ngh,ny,1:ngh) :: duz,dvz,dwz,dpz,drz
    real(wp) :: duy,dvy,dwy,dpy,dry
    !-------------------------------------------------------------------------

    ! Pointers for spherical coordinates
    ! ==================================
    ir=>BC_edge(3,1,1)%i_r
    cosphi=>BC_edge(3,1,1)%cosp
    sinphi=>BC_edge(3,1,1)%sinp
    costeta=>BC_edge(3,1,1)%cost
    sinteta=>BC_edge(3,1,1)%sint
    costcosp=>BC_edge(3,1,1)%costcosp
    costsinp=>BC_edge(3,1,1)%costsinp
    sintcosp=>BC_edge(3,1,1)%sintcosp
    sintsinp=>BC_edge(3,1,1)%sintsinp
    
    ! Compute fluctuations
    ! ====================
    do i=1,nghp3
       do k=1,nghp3
          do j=ny1,ny2
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
       do k=1,ngh
          do j=ndy,nfy
             vg(i,j,k)= BC_face(1,1)%U0(i,j,k,2)*sintcosp(i,j,k) +     &
	           	BC_face(1,1)%U0(i,j,k,3)*sintsinp(i,j,k) +     &
	           	BC_face(1,1)%U0(i,j,k,4)*costeta(i,j,k)	       &
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
    do k=1,ngh
       do j=ndy,nfy
          dux(i,j,k)= ( a06(1)*uf(1,j,k)+a06(2)*uf(2,j,k) &
                      + a06(3)*uf(3,j,k)+a06(4)*uf(4,j,k) &
                      + a06(5)*uf(5,j,k)+a06(6)*uf(6,j,k) &
                      + a06(7)*uf(7,j,k) ) *idx(i)*sintcosp(i,j,k)
          dvx(i,j,k)= ( a06(1)*vf(1,j,k)+a06(2)*vf(2,j,k) &
                      + a06(3)*vf(3,j,k)+a06(4)*vf(4,j,k) &
                      + a06(5)*vf(5,j,k)+a06(6)*vf(6,j,k) &
                      + a06(7)*vf(7,j,k) ) *idx(i)*sintcosp(i,j,k)
          dwx(i,j,k)= ( a06(1)*wf(1,j,k)+a06(2)*wf(2,j,k) &
                      + a06(3)*wf(3,j,k)+a06(4)*wf(4,j,k) &
                      + a06(5)*wf(5,j,k)+a06(6)*wf(6,j,k) &
                      + a06(7)*wf(7,j,k) ) *idx(i)*sintcosp(i,j,k)
          dpx(i,j,k)= ( a06(1)*pf(1,j,k)+a06(2)*pf(2,j,k) &
                      + a06(3)*pf(3,j,k)+a06(4)*pf(4,j,k) &
                      + a06(5)*pf(5,j,k)+a06(6)*pf(6,j,k) &
                      + a06(7)*pf(7,j,k) ) *idx(i)*sintcosp(i,j,k)
          drx(i,j,k)= ( a06(1)*rf(1,j,k)+a06(2)*rf(2,j,k) &
                      + a06(3)*rf(3,j,k)+a06(4)*rf(4,j,k) &
                      + a06(5)*rf(5,j,k)+a06(6)*rf(6,j,k) &
                      + a06(7)*rf(7,j,k) ) *idx(i)*sintcosp(i,j,k)
       enddo
    enddo

    i=2
    do k=1,ngh
       do j=ndy,nfy
          dux(i,j,k)= ( a15(1)*uf(1,j,k)+a15(2)*uf(2,j,k) &
                      + a15(3)*uf(3,j,k)+a15(4)*uf(4,j,k) &
                      + a15(5)*uf(5,j,k)+a15(6)*uf(6,j,k) &
                      + a15(7)*uf(7,j,k) ) *idx(i)*sintcosp(i,j,k)
          dvx(i,j,k)= ( a15(1)*vf(1,j,k)+a15(2)*vf(2,j,k) &
                      + a15(3)*vf(3,j,k)+a15(4)*vf(4,j,k) &
                      + a15(5)*vf(5,j,k)+a15(6)*vf(6,j,k) &
                      + a15(7)*vf(7,j,k) ) *idx(i)*sintcosp(i,j,k)
          dwx(i,j,k)= ( a15(1)*wf(1,j,k)+a15(2)*wf(2,j,k) &
                      + a15(3)*wf(3,j,k)+a15(4)*wf(4,j,k) &
                      + a15(5)*wf(5,j,k)+a15(6)*wf(6,j,k) &
                      + a15(7)*wf(7,j,k) ) *idx(i)*sintcosp(i,j,k)
          dpx(i,j,k)= ( a15(1)*pf(1,j,k)+a15(2)*pf(2,j,k) &
                      + a15(3)*pf(3,j,k)+a15(4)*pf(4,j,k) &
                      + a15(5)*pf(5,j,k)+a15(6)*pf(6,j,k) &
                      + a15(7)*pf(7,j,k) ) *idx(i)*sintcosp(i,j,k)
          drx(i,j,k)= ( a15(1)*rf(1,j,k)+a15(2)*rf(2,j,k) &
                      + a15(3)*rf(3,j,k)+a15(4)*rf(4,j,k) &
                      + a15(5)*rf(5,j,k)+a15(6)*rf(6,j,k) &
                      + a15(7)*rf(7,j,k) ) *idx(i)*sintcosp(i,j,k)
       enddo
    enddo

    i=3
    do k=1,ngh
       do j=ndy,nfy
          dux(i,j,k)= ( a24(1)*uf(1,j,k)+a24(2)*uf(2,j,k) &
                      + a24(3)*uf(3,j,k)+a24(4)*uf(4,j,k) &
                      + a24(5)*uf(5,j,k)+a24(6)*uf(6,j,k) &
                      + a24(7)*uf(7,j,k) ) *idx(i)*sintcosp(i,j,k)
          dvx(i,j,k)= ( a24(1)*vf(1,j,k)+a24(2)*vf(2,j,k) &
                      + a24(3)*vf(3,j,k)+a24(4)*vf(4,j,k) &
                      + a24(5)*vf(5,j,k)+a24(6)*vf(6,j,k) &
                      + a24(7)*vf(7,j,k) ) *idx(i)*sintcosp(i,j,k)
          dwx(i,j,k)= ( a24(1)*wf(1,j,k)+a24(2)*wf(2,j,k) &
                      + a24(3)*wf(3,j,k)+a24(4)*wf(4,j,k) &
                      + a24(5)*wf(5,j,k)+a24(6)*wf(6,j,k) &
                      + a24(7)*wf(7,j,k) ) *idx(i)*sintcosp(i,j,k)
          dpx(i,j,k)= ( a24(1)*pf(1,j,k)+a24(2)*pf(2,j,k) &
                      + a24(3)*pf(3,j,k)+a24(4)*pf(4,j,k) &
                      + a24(5)*pf(5,j,k)+a24(6)*pf(6,j,k) &
                      + a24(7)*pf(7,j,k) ) *idx(i)*sintcosp(i,j,k)
          drx(i,j,k)= ( a24(1)*rf(1,j,k)+a24(2)*rf(2,j,k) &
                      + a24(3)*rf(3,j,k)+a24(4)*rf(4,j,k) &
                      + a24(5)*rf(5,j,k)+a24(6)*rf(6,j,k) &
                      + a24(7)*rf(7,j,k) ) *idx(i)*sintcosp(i,j,k)
       enddo
    enddo

    do i=4,ngh
       do k=1,ngh
          do j=ndy,nfy
             dux(i,j,k)= ( a7(1)* (uf(i+1,j,k)-uf(i-1,j,k)) &
                         + a7(2)* (uf(i+2,j,k)-uf(i-2,j,k)) &
                         + a7(3)* (uf(i+3,j,k)-uf(i-3,j,k)) ) *idx(i)*sintcosp(i,j,k)
             dvx(i,j,k)= ( a7(1)* (vf(i+1,j,k)-vf(i-1,j,k)) &
                         + a7(2)* (vf(i+2,j,k)-vf(i-2,j,k)) &
                         + a7(3)* (vf(i+3,j,k)-vf(i-3,j,k)) ) *idx(i)*sintcosp(i,j,k)
             dwx(i,j,k)= ( a7(1)* (wf(i+1,j,k)-wf(i-1,j,k)) &
                         + a7(2)* (wf(i+2,j,k)-wf(i-2,j,k)) &
                         + a7(3)* (wf(i+3,j,k)-wf(i-3,j,k)) ) *idx(i)*sintcosp(i,j,k)
             dpx(i,j,k)= ( a7(1)* (pf(i+1,j,k)-pf(i-1,j,k)) &
                         + a7(2)* (pf(i+2,j,k)-pf(i-2,j,k)) &
                         + a7(3)* (pf(i+3,j,k)-pf(i-3,j,k)) ) *idx(i)*sintcosp(i,j,k)
             drx(i,j,k)= ( a7(1)* (rf(i+1,j,k)-rf(i-1,j,k)) &
                         + a7(2)* (rf(i+2,j,k)-rf(i-2,j,k)) &
                         + a7(3)* (rf(i+3,j,k)-rf(i-3,j,k)) ) *idx(i)*sintcosp(i,j,k)
          enddo
       enddo
    enddo

    ! Non-centered derivatives in z-direction *cos(teta)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    k=1
    do i=1,ngh
       do j=ndy,nfy
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
    do i=1,ngh
       do j=ndy,nfy
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
    do i=1,ngh
       do j=ndy,nfy
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
       do i=1,ngh
          do j=ndy,nfy
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
       do i=1,ngh
          do j=ndy,nfy

             duy= ( a7(1)*(uf(i,j+1,k)-uf(i,j-1,k)) &
                  + a7(2)*(uf(i,j+2,k)-uf(i,j-2,k)) &
                  + a7(3)*(uf(i,j+3,k)-uf(i,j-3,k)) ) *idy(j)*sintsinp(i,j,k)
             dvy= ( a7(1)*(vf(i,j+1,k)-vf(i,j-1,k)) &
                  + a7(2)*(vf(i,j+2,k)-vf(i,j-2,k)) &
                  + a7(3)*(vf(i,j+3,k)-vf(i,j-3,k)) ) *idy(j)*sintsinp(i,j,k)
             dwy= ( a7(1)*(wf(i,j+1,k)-wf(i,j-1,k)) &
                  + a7(2)*(wf(i,j+2,k)-wf(i,j-2,k)) &
                  + a7(3)*(wf(i,j+3,k)-wf(i,j-3,k)) ) *idy(j)*sintsinp(i,j,k)
             dpy= ( a7(1)*(pf(i,j+1,k)-pf(i,j-1,k)) &
                  + a7(2)*(pf(i,j+2,k)-pf(i,j-2,k)) &
                  + a7(3)*(pf(i,j+3,k)-pf(i,j-3,k)) ) *idy(j)*sintsinp(i,j,k)
             dry= ( a7(1)*(rf(i,j+1,k)-rf(i,j-1,k)) &
                  + a7(2)*(rf(i,j+2,k)-rf(i,j-2,k)) &
                  + a7(3)*(rf(i,j+3,k)-rf(i,j-3,k)) ) *idy(j)*sintsinp(i,j,k)

             pt(i,j,k) = vg(i,j,k)*(dpx(i,j,k)+dpy+dpz(i,j,k)+pf(i,j,k)*ir(i,j,k))
             ut(i,j,k) = vg(i,j,k)*(dux(i,j,k)+duy+duz(i,j,k)+uf(i,j,k)*ir(i,j,k))
             vt(i,j,k) = vg(i,j,k)*(dvx(i,j,k)+dvy+dvz(i,j,k)+vf(i,j,k)*ir(i,j,k))
             wt(i,j,k) = vg(i,j,k)*(dwx(i,j,k)+dwy+dwz(i,j,k)+wf(i,j,k)*ir(i,j,k))
             rt(i,j,k) = vg(i,j,k)*(drx(i,j,k)+dry+drz(i,j,k)+rf(i,j,k)*ir(i,j,k))
          enddo
       enddo
    enddo

    ! Update fluxes at each RK step
    ! =============================
    do j=ndy,nfy
       do k=1,ngh
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

  end subroutine bc_TD3d_kmin_imin

  !===============================================================================
  module subroutine bc_TD3d_kmin_imax
  !===============================================================================
    !> 3D Tam & Dong's BC: boundary condition at kmin-imax (edge 3,1,2 /front-right)
    !> - Cartesian version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l
    real(wp) :: cp,av,c2_
    real(wp), dimension(-2:ngh,ny1:ny2,1:ngh+3) :: rf,uf,vf,wf,pf
    real(wp), dimension(1:ngh,ny,1:ngh) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(1:ngh,ny,1:ngh) :: dux,dvx,dwx,dpx,drx
    real(wp), dimension(1:ngh,ny,1:ngh) :: duz,dvz,dwz,dpz,drz
    real(wp) :: duy,dvy,dwy,dpy,dry
    !-------------------------------------------------------------------------

    ! Pointers for spherical coordinates
    ! ==================================
    ir=>BC_edge(3,1,2)%i_r
    cosphi=>BC_edge(3,1,2)%cosp
    sinphi=>BC_edge(3,1,2)%sinp
    costeta=>BC_edge(3,1,2)%cost
    sinteta=>BC_edge(3,1,2)%sint
    costcosp=>BC_edge(3,1,2)%costcosp
    costsinp=>BC_edge(3,1,2)%costsinp
    sintcosp=>BC_edge(3,1,2)%sintcosp
    sintsinp=>BC_edge(3,1,2)%sintsinp
    
    ! Compute fluctuations
    ! ====================
    do k=1,nghp3
       do i=nxmngh-2,nx
          l=i-nxmngh
          do j=ny1,ny2
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
    do k=1,ngh
       do l=1,ngh
          do j=ndy,nfy
             vg(l,j,k)= BC_face(1,2)%U0(l,j,k,2)*sintcosp(l,j,k) +     &
	           	BC_face(1,2)%U0(l,j,k,3)*sintsinp(l,j,k) +     &
	           	BC_face(1,2)%U0(l,j,k,4)*costeta(l,j,k)	       &
                      + sqrt( BC_face(1,2)%U0(l,j,k,6)-                & 				
                       ( BC_face(1,2)%U0(l,j,k,2)*costcosp(l,j,k) +    &				
                 	 BC_face(1,2)%U0(l,j,k,3)*costsinp(l,j,k) -    & 
          	         BC_face(1,2)%U0(l,j,k,4)*sinteta(l,j,k) )**2- &
	               ( BC_face(1,2)%U0(l,j,k,2)*sinphi(l,j,k) -      &
	      	         BC_face(1,2)%U0(l,j,k,3)*cosphi(l,j,k) )**2 )
          enddo
       enddo
    enddo

    ! Non-centered derivatives in x-direction *sin(teta)*cos(phi)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    do i=nxmnghp1,nx-3
       l=i-nxmngh
       do k=1,ngh
          do j=ndy,nfy
             dux(l,j,k) = ( a7(1)* (uf(l+1,j,k)-uf(l-1,j,k)) &
                          + a7(2)* (uf(l+2,j,k)-uf(l-2,j,k)) &
                          + a7(3)* (uf(l+3,j,k)-uf(l-3,j,k)) ) *idx(i)*sintcosp(l,j,k)
             dvx(l,j,k) = ( a7(1)* (vf(l+1,j,k)-vf(l-1,j,k)) &
                          + a7(2)* (vf(l+2,j,k)-vf(l-2,j,k)) &
                          + a7(3)* (vf(l+3,j,k)-vf(l-3,j,k)) ) *idx(i)*sintcosp(l,j,k)
             dwx(l,j,k) = ( a7(1)* (wf(l+1,j,k)-wf(l-1,j,k)) &
                          + a7(2)* (wf(l+2,j,k)-wf(l-2,j,k)) &
                          + a7(3)* (wf(l+3,j,k)-wf(l-3,j,k)) ) *idx(i)*sintcosp(l,j,k)
             dpx(l,j,k) = ( a7(1)* (pf(l+1,j,k)-pf(l-1,j,k)) &
                          + a7(2)* (pf(l+2,j,k)-pf(l-2,j,k)) &
                          + a7(3)* (pf(l+3,j,k)-pf(l-3,j,k)) ) *idx(i)*sintcosp(l,j,k)
             drx(l,j,k) = ( a7(1)* (rf(l+1,j,k)-rf(l-1,j,k)) &
                          + a7(2)* (rf(l+2,j,k)-rf(l-2,j,k)) &
                          + a7(3)* (rf(l+3,j,k)-rf(l-3,j,k)) ) *idx(i)*sintcosp(l,j,k)
          enddo
       enddo
    enddo
    
    i=nx-2
    l=i-nxmngh
    do k=1,ngh
       do j=ndy,nfy
          dux(l,j,k) = ( a42(1)*uf(l+2,j,k)+a42(2)*uf(l+1,j,k) &
                       + a42(3)*uf(l  ,j,k)+a42(4)*uf(l-1,j,k) &
                       + a42(5)*uf(l-2,j,k)+a42(6)*uf(l-3,j,k) &
                       + a42(7)*uf(l-4,j,k) ) *idx(i)*sintcosp(l,j,k)
          dvx(l,j,k) = ( a42(1)*vf(l+2,j,k)+a42(2)*vf(l+1,j,k) &
                       + a42(3)*vf(l  ,j,k)+a42(4)*vf(l-1,j,k) &
                       + a42(5)*vf(l-2,j,k)+a42(6)*vf(l-3,j,k) &
                       + a42(7)*vf(l-4,j,k) ) *idx(i)*sintcosp(l,j,k)
          dwx(l,j,k) = ( a42(1)*wf(l+2,j,k)+a42(2)*wf(l+1,j,k) &
                       + a42(3)*wf(l  ,j,k)+a42(4)*wf(l-1,j,k) &
                       + a42(5)*wf(l-2,j,k)+a42(6)*wf(l-3,j,k) &
                       + a42(7)*wf(l-4,j,k) ) *idx(i)*sintcosp(l,j,k)
          dpx(l,j,k) = ( a42(1)*pf(l+2,j,k)+a42(2)*pf(l+1,j,k) &
                       + a42(3)*pf(l  ,j,k)+a42(4)*pf(l-1,j,k) &
                       + a42(5)*pf(l-2,j,k)+a42(6)*pf(l-3,j,k) &
                       + a42(7)*pf(l-4,j,k) ) *idx(i)*sintcosp(l,j,k)
          drx(l,j,k) = ( a42(1)*rf(l+2,j,k)+a42(2)*rf(l+1,j,k) &
                       + a42(3)*rf(l  ,j,k)+a42(4)*rf(l-1,j,k) &
                       + a42(5)*rf(l-2,j,k)+a42(6)*rf(l-3,j,k) &
                       + a42(7)*rf(l-4,j,k) ) *idx(i)*sintcosp(l,j,k)
       enddo
    enddo
    
    i=nx-1
    l=i-nxmngh
    do k=1,ngh
       do j=ndy,nfy
          dux(l,j,k) = ( a51(1)*uf(l+1,j,k)+a51(2)*uf(l  ,j,k) &
                       + a51(3)*uf(l-1,j,k)+a51(4)*uf(l-2,j,k) &
                       + a51(5)*uf(l-3,j,k)+a51(6)*uf(l-4,j,k) &
                       + a51(7)*uf(l-5,j,k) ) *idx(i)*sintcosp(l,j,k)
          dvx(l,j,k) = ( a51(1)*vf(l+1,j,k)+a51(2)*vf(l  ,j,k) &
                       + a51(3)*vf(l-1,j,k)+a51(4)*vf(l-2,j,k) &
                       + a51(5)*vf(l-3,j,k)+a51(6)*vf(l-4,j,k) &
                       + a51(7)*vf(l-5,j,k) ) *idx(i)*sintcosp(l,j,k)
          dwx(l,j,k) = ( a51(1)*wf(l+1,j,k)+a51(2)*wf(l  ,j,k) &
                       + a51(3)*wf(l-1,j,k)+a51(4)*wf(l-2,j,k) &
                       + a51(5)*wf(l-3,j,k)+a51(6)*wf(l-4,j,k) &
                       + a51(7)*wf(l-5,j,k) ) *idx(i)*sintcosp(l,j,k)
          dpx(l,j,k) = ( a51(1)*pf(l+1,j,k)+a51(2)*pf(l  ,j,k) &
                       + a51(3)*pf(l-1,j,k)+a51(4)*pf(l-2,j,k) &
                       + a51(5)*pf(l-3,j,k)+a51(6)*pf(l-4,j,k) &
                       + a51(7)*pf(l-5,j,k) ) *idx(i)*sintcosp(l,j,k)
          drx(l,j,k) = ( a51(1)*rf(l+1,j,k)+a51(2)*rf(l  ,j,k) &
                       + a51(3)*rf(l-1,j,k)+a51(4)*rf(l-2,j,k) &
                       + a51(5)*rf(l-3,j,k)+a51(6)*rf(l-4,j,k) &
                       + a51(7)*rf(l-5,j,k) ) *idx(i)*sintcosp(l,j,k)
       enddo
    enddo
    
    i=nx
    l=i-nxmngh
    do k=1,ngh
       do j=ndy,nfy
          dux(l,j,k) = ( a60(1)*uf(l  ,j,k)+a60(2)*uf(l-1,j,k) &
                       + a60(3)*uf(l-2,j,k)+a60(4)*uf(l-3,j,k) &
                       + a60(5)*uf(l-4,j,k)+a60(6)*uf(l-5,j,k) &
                       + a60(7)*uf(l-6,j,k) ) *idx(i)*sintcosp(l,j,k)
          dvx(l,j,k) = ( a60(1)*vf(l  ,j,k)+a60(2)*vf(l-1,j,k) &
                       + a60(3)*vf(l-2,j,k)+a60(4)*vf(l-3,j,k) &
                       + a60(5)*vf(l-4,j,k)+a60(6)*vf(l-5,j,k) &
                       + a60(7)*vf(l-6,j,k) ) *idx(i)*sintcosp(l,j,k)
          dwx(l,j,k) = ( a60(1)*wf(l  ,j,k)+a60(2)*wf(l-1,j,k) &
                       + a60(3)*wf(l-2,j,k)+a60(4)*wf(l-3,j,k) &
                       + a60(5)*wf(l-4,j,k)+a60(6)*wf(l-5,j,k) &
                       + a60(7)*wf(l-6,j,k) ) *idx(i)*sintcosp(l,j,k)
          dpx(l,j,k) = ( a60(1)*pf(l  ,j,k)+a60(2)*pf(l-1,j,k) &
                       + a60(3)*pf(l-2,j,k)+a60(4)*pf(l-3,j,k) &
                       + a60(5)*pf(l-4,j,k)+a60(6)*pf(l-5,j,k) &
                       + a60(7)*pf(l-6,j,k) ) *idx(i)*sintcosp(l,j,k)
          drx(l,j,k) = ( a60(1)*rf(l  ,j,k)+a60(2)*rf(l-1,j,k) &
                       + a60(3)*rf(l-2,j,k)+a60(4)*rf(l-3,j,k) &
                       + a60(5)*rf(l-4,j,k)+a60(6)*rf(l-5,j,k) &
                       + a60(7)*rf(l-6,j,k) ) *idx(i)*sintcosp(l,j,k)
       enddo
    enddo
    
    ! Non-centered derivatives in z-direction *cos(teta)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    k=1
    do l=1,ngh
       do j=ndy,nfy
          duz(l,j,k) = ( a06(1)*uf(l,j,1)+a06(2)*uf(l,j,2) &
                       + a06(3)*uf(l,j,3)+a06(4)*uf(l,j,4) &
                       + a06(5)*uf(l,j,5)+a06(6)*uf(l,j,6) &
                       + a06(7)*uf(l,j,7) ) *idz(k)*costeta(l,j,k)
          dvz(l,j,k) = ( a06(1)*vf(l,j,1)+a06(2)*vf(l,j,2) &
                       + a06(3)*vf(l,j,3)+a06(4)*vf(l,j,4) &
                       + a06(5)*vf(l,j,5)+a06(6)*vf(l,j,6) &
                       + a06(7)*vf(l,j,7) ) *idz(k)*costeta(l,j,k)
          dwz(l,j,k) = ( a06(1)*wf(l,j,1)+a06(2)*wf(l,j,2) &
                       + a06(3)*wf(l,j,3)+a06(4)*wf(l,j,4) &
                       + a06(5)*wf(l,j,5)+a06(6)*wf(l,j,6) &
                       + a06(7)*wf(l,j,7) ) *idz(k)*costeta(l,j,k)
          dpz(l,j,k) = ( a06(1)*pf(l,j,1)+a06(2)*pf(l,j,2) &
                       + a06(3)*pf(l,j,3)+a06(4)*pf(l,j,4) &
                       + a06(5)*pf(l,j,5)+a06(6)*pf(l,j,6) &
                       + a06(7)*pf(l,j,7) ) *idz(k)*costeta(l,j,k)
          drz(l,j,k) = ( a06(1)*rf(l,j,1)+a06(2)*rf(l,j,2) &
                       + a06(3)*rf(l,j,3)+a06(4)*rf(l,j,4) &
                       + a06(5)*rf(l,j,5)+a06(6)*rf(l,j,6) &
                       + a06(7)*rf(l,j,7) ) *idz(k)*costeta(l,j,k)
       enddo
    enddo

    k=2
    do l=1,ngh
       do j=ndy,nfy
          duz(l,j,k) = ( a15(1)*uf(l,j,1)+a15(2)*uf(l,j,2) &
                       + a15(3)*uf(l,j,3)+a15(4)*uf(l,j,4) &
                       + a15(5)*uf(l,j,5)+a15(6)*uf(l,j,6) &
                       + a15(7)*uf(l,j,7) ) *idz(k)*costeta(l,j,k)
          dvz(l,j,k) = ( a15(1)*vf(l,j,1)+a15(2)*vf(l,j,2) &
                       + a15(3)*vf(l,j,3)+a15(4)*vf(l,j,4) &
                       + a15(5)*vf(l,j,5)+a15(6)*vf(l,j,6) &
                       + a15(7)*vf(l,j,7) ) *idz(k)*costeta(l,j,k)
          dwz(l,j,k) = ( a15(1)*wf(l,j,1)+a15(2)*wf(l,j,2) &
                       + a15(3)*wf(l,j,3)+a15(4)*wf(l,j,4) &
                       + a15(5)*wf(l,j,5)+a15(6)*wf(l,j,6) &
                       + a15(7)*wf(l,j,7) ) *idz(k)*costeta(l,j,k)
          dpz(l,j,k) = ( a15(1)*pf(l,j,1)+a15(2)*pf(l,j,2) &
                       + a15(3)*pf(l,j,3)+a15(4)*pf(l,j,4) &
                       + a15(5)*pf(l,j,5)+a15(6)*pf(l,j,6) &
                       + a15(7)*pf(l,j,7) ) *idz(k)*costeta(l,j,k)
          drz(l,j,k) = ( a15(1)*rf(l,j,1)+a15(2)*rf(l,j,2) &
                       + a15(3)*rf(l,j,3)+a15(4)*rf(l,j,4) &
                       + a15(5)*rf(l,j,5)+a15(6)*rf(l,j,6) &
                       + a15(7)*rf(l,j,7) ) *idz(k)*costeta(l,j,k)
       enddo
    enddo

    k=3
    do l=1,ngh
       do j=ndy,nfy
          duz(l,j,k) = ( a24(1)*uf(l,j,1)+a24(2)*uf(l,j,2) &
                       + a24(3)*uf(l,j,3)+a24(4)*uf(l,j,4) &
                       + a24(5)*uf(l,j,5)+a24(6)*uf(l,j,6) &
                       + a24(7)*uf(l,j,7) ) *idz(k)*costeta(l,j,k)
          dvz(l,j,k) = ( a24(1)*vf(l,j,1)+a24(2)*vf(l,j,2) &
                       + a24(3)*vf(l,j,3)+a24(4)*vf(l,j,4) &
                       + a24(5)*vf(l,j,5)+a24(6)*vf(l,j,6) &
                       + a24(7)*vf(l,j,7) ) *idz(k)*costeta(l,j,k)
          dwz(l,j,k) = ( a24(1)*wf(l,j,1)+a24(2)*wf(l,j,2) &
                       + a24(3)*wf(l,j,3)+a24(4)*wf(l,j,4) &
                       + a24(5)*wf(l,j,5)+a24(6)*wf(l,j,6) &
                       + a24(7)*wf(l,j,7) ) *idz(k)*costeta(l,j,k)
          dpz(l,j,k) = ( a24(1)*pf(l,j,1)+a24(2)*pf(l,j,2) &
                       + a24(3)*pf(l,j,3)+a24(4)*pf(l,j,4) &
                       + a24(5)*pf(l,j,5)+a24(6)*pf(l,j,6) &
                       + a24(7)*pf(l,j,7) ) *idz(k)*costeta(l,j,k)
          drz(l,j,k) = ( a24(1)*rf(l,j,1)+a24(2)*rf(l,j,2) &
                       + a24(3)*rf(l,j,3)+a24(4)*rf(l,j,4) &
                       + a24(5)*rf(l,j,5)+a24(6)*rf(l,j,6) &
                       + a24(7)*rf(l,j,7) ) *idz(k)*costeta(l,j,k)
       enddo
    enddo

    do k=4,ngh
       do l=1,ngh
          do j=ndy,nfy
             duz(l,j,k) = (a7(1)*(uf(l,j,k+1) - uf(l,j,k-1)) + &
                           a7(2)*(uf(l,j,k+2) - uf(l,j,k-2)) + &
                           a7(3)*(uf(l,j,k+3) - uf(l,j,k-3)) ) *idz(k)*costeta(l,j,k)
             dvz(l,j,k) = (a7(1)*(vf(l,j,k+1) - vf(l,j,k-1)) + &
                           a7(2)*(vf(l,j,k+2) - vf(l,j,k-2)) + &
                           a7(3)*(vf(l,j,k+3) - vf(l,j,k-3)) ) *idz(k)*costeta(l,j,k)
             dwz(l,j,k) = (a7(1)*(wf(l,j,k+1) - wf(l,j,k-1)) + &
                           a7(2)*(wf(l,j,k+2) - wf(l,j,k-2)) + &
                           a7(3)*(wf(l,j,k+3) - wf(l,j,k-3)) ) *idz(k)*costeta(l,j,k)
             dpz(l,j,k) = (a7(1)*(pf(l,j,k+1) - pf(l,j,k-1)) + &
                           a7(2)*(pf(l,j,k+2) - pf(l,j,k-2)) + &
                           a7(3)*(pf(l,j,k+3) - pf(l,j,k-3)) ) *idz(k)*costeta(l,j,k)
             drz(l,j,k) = (a7(1)*(rf(l,j,k+1) - rf(l,j,k-1)) + &
                           a7(2)*(rf(l,j,k+2) - rf(l,j,k-2)) + &
                           a7(3)*(rf(l,j,k+3) - rf(l,j,k-3)) ) *idz(k)*costeta(l,j,k)
          enddo
       enddo
    enddo

    ! Compute rt,ut,vt,wt & pt
    ! ========================
    do k=1,ngh
       do l=1,ngh
          do j=ndy,nfy
             duy= ( a7(1)*(uf(l,j+1,k)-uf(l,j-1,k)) &
                  + a7(2)*(uf(l,j+2,k)-uf(l,j-2,k)) &
                  + a7(3)*(uf(l,j+3,k)-uf(l,j-3,k)) ) *idy(j)*sintsinp(l,j,k)
             dvy= ( a7(1)*(vf(l,j+1,k)-vf(l,j-1,k)) &
                  + a7(2)*(vf(l,j+2,k)-vf(l,j-2,k)) &
                  + a7(3)*(vf(l,j+3,k)-vf(l,j-3,k)) ) *idy(j)*sintsinp(l,j,k)
             dwy= ( a7(1)*(wf(l,j+1,k)-wf(l,j-1,k)) &
                  + a7(2)*(wf(l,j+2,k)-wf(l,j-2,k)) &
                  + a7(3)*(wf(l,j+3,k)-wf(l,j-3,k)) ) *idy(j)*sintsinp(l,j,k)
             dpy= ( a7(1)*(pf(l,j+1,k)-pf(l,j-1,k)) &
                  + a7(2)*(pf(l,j+2,k)-pf(l,j-2,k)) &
                  + a7(3)*(pf(l,j+3,k)-pf(l,j-3,k)) ) *idy(j)*sintsinp(l,j,k)
             dry= ( a7(1)*(rf(l,j+1,k)-rf(l,j-1,k)) &
                  + a7(2)*(rf(l,j+2,k)-rf(l,j-2,k)) &
                  + a7(3)*(rf(l,j+3,k)-rf(l,j-3,k)) ) *idy(j)*sintsinp(l,j,k)

             pt(l,j,k) = vg(l,j,k)*(dpx(l,j,k)+dpy+dpz(l,j,k)+pf(l,j,k)*ir(l,j,k))
             ut(l,j,k) = vg(l,j,k)*(dux(l,j,k)+duy+duz(l,j,k)+uf(l,j,k)*ir(l,j,k))
             vt(l,j,k) = vg(l,j,k)*(dvx(l,j,k)+dvy+dvz(l,j,k)+vf(l,j,k)*ir(l,j,k))
             wt(l,j,k) = vg(l,j,k)*(dwx(l,j,k)+dwy+dwz(l,j,k)+wf(l,j,k)*ir(l,j,k))
             rt(l,j,k) = vg(l,j,k)*(drx(l,j,k)+dry+drz(l,j,k)+rf(l,j,k)*ir(l,j,k))
          enddo
       enddo
    enddo

    ! Update fluxes at each RK step
    ! =============================
    do j=ndy,nfy
       do i=nxmnghp1,nx
          l=i-nxmngh
          do k=1,ngh
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

  end subroutine bc_TD3d_kmin_imax

  !===============================================================================
  module subroutine bc_TD3d_kmax_imin
  !===============================================================================
    !> 3D Tam & Dong's BC: boundary condition at kmax-imin (edge 3,2,1 /back-left)
    !> - Cartesian version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l
    real(wp) :: cp,av,c2_
    real(wp), dimension(1:ngh+3,ny1:ny2,-2:ngh) :: rf,uf,vf,wf,pf
    real(wp), dimension(1:ngh,ny,1:ngh) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(1:ngh,ny,1:ngh) :: dux,dvx,dwx,dpx,drx
    real(wp), dimension(1:ngh,ny,1:ngh) :: duz,dvz,dwz,dpz,drz
    real(wp) :: duy,dvy,dwy,dpy,dry
    !-------------------------------------------------------------------------
   
    ! Pointers for spherical coordinates
    ! ==================================
    ir=>BC_edge(3,2,1)%i_r
    cosphi=>BC_edge(3,2,1)%cosp
    sinphi=>BC_edge(3,2,1)%sinp
    costeta=>BC_edge(3,2,1)%cost
    sinteta=>BC_edge(3,2,1)%sint
    costcosp=>BC_edge(3,2,1)%costcosp
    costsinp=>BC_edge(3,2,1)%costsinp
    sintcosp=>BC_edge(3,2,1)%sintcosp
    sintsinp=>BC_edge(3,2,1)%sintsinp
    
    ! Compute fluctuations
    ! ====================
    do k=nzmngh-2,nz
       l=k-nzmngh
       do i=1,nghp3
          do j=ny1,ny2
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
    do k=nzmnghp1,nz
       l=k-nzmngh
       do i=1,ngh
          do j=ndy,nfy
             vg(i,j,l)= BC_face(1,1)%U0(i,j,k,2)*sintcosp(i,j,l) +     &
	           	BC_face(1,1)%U0(i,j,k,3)*sintsinp(i,j,l) +     &
	           	BC_face(1,1)%U0(i,j,k,4)*costeta(i,j,l)	       &
                      + sqrt( BC_face(1,1)%U0(i,j,k,6)-                & 				
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
    do k=nzmnghp1,nz
       l=k-nzmngh
       do j=ndy,nfy
          dux(i,j,l)= ( a06(1)*uf(1,j,l)+a06(2)*uf(2,j,l) &
                      + a06(3)*uf(3,j,l)+a06(4)*uf(4,j,l) &
                      + a06(5)*uf(5,j,l)+a06(6)*uf(6,j,l) &
                      + a06(7)*uf(7,j,l) ) *idx(i)*sintcosp(i,j,l)
          dvx(i,j,l)= ( a06(1)*vf(1,j,l)+a06(2)*vf(2,j,l) &
                      + a06(3)*vf(3,j,l)+a06(4)*vf(4,j,l) &
                      + a06(5)*vf(5,j,l)+a06(6)*vf(6,j,l) &
                      + a06(7)*vf(7,j,l) ) *idx(i)*sintcosp(i,j,l)
          dwx(i,j,l)= ( a06(1)*wf(1,j,l)+a06(2)*wf(2,j,l) &
                      + a06(3)*wf(3,j,l)+a06(4)*wf(4,j,l) &
                      + a06(5)*wf(5,j,l)+a06(6)*wf(6,j,l) &
                      + a06(7)*wf(7,j,l) ) *idx(i)*sintcosp(i,j,l)
          dpx(i,j,l)= ( a06(1)*pf(1,j,l)+a06(2)*pf(2,j,l) &
                      + a06(3)*pf(3,j,l)+a06(4)*pf(4,j,l) &
                      + a06(5)*pf(5,j,l)+a06(6)*pf(6,j,l) &
                      + a06(7)*pf(7,j,l) ) *idx(i)*sintcosp(i,j,l)
          drx(i,j,l)= ( a06(1)*rf(1,j,l)+a06(2)*rf(2,j,l) &
                      + a06(3)*rf(3,j,l)+a06(4)*rf(4,j,l) &
                      + a06(5)*rf(5,j,l)+a06(6)*rf(6,j,l) &
                      + a06(7)*rf(7,j,l) ) *idx(i)*sintcosp(i,j,l)
       enddo
    enddo

    i=2
    do k=nzmnghp1,nz
       l=k-nzmngh
       do j=ndy,nfy
          dux(i,j,l)= ( a15(1)*uf(1,j,l)+a15(2)*uf(2,j,l) &
                      + a15(3)*uf(3,j,l)+a15(4)*uf(4,j,l) &
                      + a15(5)*uf(5,j,l)+a15(6)*uf(6,j,l) &
                      + a15(7)*uf(7,j,l) ) *idx(i)*sintcosp(i,j,l)
          dvx(i,j,l)= ( a15(1)*vf(1,j,l)+a15(2)*vf(2,j,l) &
                      + a15(3)*vf(3,j,l)+a15(4)*vf(4,j,l) &
                      + a15(5)*vf(5,j,l)+a15(6)*vf(6,j,l) &
                      + a15(7)*vf(7,j,l) ) *idx(i)*sintcosp(i,j,l)
          dwx(i,j,l)= ( a15(1)*wf(1,j,l)+a15(2)*wf(2,j,l) &
                      + a15(3)*wf(3,j,l)+a15(4)*wf(4,j,l) &
                      + a15(5)*wf(5,j,l)+a15(6)*wf(6,j,l) &
                      + a15(7)*wf(7,j,l) ) *idx(i)*sintcosp(i,j,l)
          dpx(i,j,l)= ( a15(1)*pf(1,j,l)+a15(2)*pf(2,j,l) &
                      + a15(3)*pf(3,j,l)+a15(4)*pf(4,j,l) &
                      + a15(5)*pf(5,j,l)+a15(6)*pf(6,j,l) &
                      + a15(7)*pf(7,j,l) ) *idx(i)*sintcosp(i,j,l)
          drx(i,j,l)= ( a15(1)*rf(1,j,l)+a15(2)*rf(2,j,l) &
                      + a15(3)*rf(3,j,l)+a15(4)*rf(4,j,l) &
                      + a15(5)*rf(5,j,l)+a15(6)*rf(6,j,l) &
                      + a15(7)*rf(7,j,l) ) *idx(i)*sintcosp(i,j,l)
       enddo
    enddo

    i=3
    do k=nzmnghp1,nz
       l=k-nzmngh
       do j=ndy,nfy
          dux(i,j,l)= ( a24(1)*uf(1,j,l)+a24(2)*uf(2,j,l) &
                      + a24(3)*uf(3,j,l)+a24(4)*uf(4,j,l) &
                      + a24(5)*uf(5,j,l)+a24(6)*uf(6,j,l) &
                      + a24(7)*uf(7,j,l) ) *idx(i)*sintcosp(i,j,l)
          dvx(i,j,l)= ( a24(1)*vf(1,j,l)+a24(2)*vf(2,j,l) &
                      + a24(3)*vf(3,j,l)+a24(4)*vf(4,j,l) &
                      + a24(5)*vf(5,j,l)+a24(6)*vf(6,j,l) &
                      + a24(7)*vf(7,j,l) ) *idx(i)*sintcosp(i,j,l)
          dwx(i,j,l)= ( a24(1)*wf(1,j,l)+a24(2)*wf(2,j,l) &
                      + a24(3)*wf(3,j,l)+a24(4)*wf(4,j,l) &
                      + a24(5)*wf(5,j,l)+a24(6)*wf(6,j,l) &
                      + a24(7)*wf(7,j,l) ) *idx(i)*sintcosp(i,j,l)
          dpx(i,j,l)= ( a24(1)*pf(1,j,l)+a24(2)*pf(2,j,l) &
                      + a24(3)*pf(3,j,l)+a24(4)*pf(4,j,l) &
                      + a24(5)*pf(5,j,l)+a24(6)*pf(6,j,l) &
                      + a24(7)*pf(7,j,l) ) *idx(i)*sintcosp(i,j,l)
          drx(i,j,l)= ( a24(1)*rf(1,j,l)+a24(2)*rf(2,j,l) &
                      + a24(3)*rf(3,j,l)+a24(4)*rf(4,j,l) &
                      + a24(5)*rf(5,j,l)+a24(6)*rf(6,j,l) &
                      + a24(7)*rf(7,j,l) ) *idx(i)*sintcosp(i,j,l)
       enddo
    enddo

    do i=4,ngh
       do k=nzmnghp1,nz
          l=k-nzmngh
          do j=ndy,nfy
             dux(i,j,l)= ( a7(1)* (uf(i+1,j,l)-uf(i-1,j,l)) &
                         + a7(2)* (uf(i+2,j,l)-uf(i-2,j,l)) &
                         + a7(3)* (uf(i+3,j,l)-uf(i-3,j,l)) ) *idx(i)*sintcosp(i,j,l)
             dvx(i,j,l)= ( a7(1)* (vf(i+1,j,l)-vf(i-1,j,l)) &
                         + a7(2)* (vf(i+2,j,l)-vf(i-2,j,l)) &
                         + a7(3)* (vf(i+3,j,l)-vf(i-3,j,l)) ) *idx(i)*sintcosp(i,j,l)
             dwx(i,j,l)= ( a7(1)* (wf(i+1,j,l)-wf(i-1,j,l)) &
                         + a7(2)* (wf(i+2,j,l)-wf(i-2,j,l)) &
                         + a7(3)* (wf(i+3,j,l)-wf(i-3,j,l)) ) *idx(i)*sintcosp(i,j,l)
             dpx(i,j,l)= ( a7(1)* (pf(i+1,j,l)-pf(i-1,j,l)) &
                         + a7(2)* (pf(i+2,j,l)-pf(i-2,j,l)) &
                         + a7(3)* (pf(i+3,j,l)-pf(i-3,j,l)) ) *idx(i)*sintcosp(i,j,l)
             drx(i,j,l)= ( a7(1)* (rf(i+1,j,l)-rf(i-1,j,l)) &
                         + a7(2)* (rf(i+2,j,l)-rf(i-2,j,l)) &
                         + a7(3)* (rf(i+3,j,l)-rf(i-3,j,l)) ) *idx(i)*sintcosp(i,j,l)
          enddo
       enddo
    enddo

    ! Non-centered derivatives in z-direction *cos(teta)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    do k=nzmnghp1,nz-3
       l=k-nzmngh
       do i=1,ngh
          do j=ndy,nfy
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
    do i=1,ngh
       do j=ndy,nfy
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
    do i=1,ngh
       do j=ndy,nfy
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
    do i=1,ngh
       do j=ndy,nfy
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
       do i=1,ngh
          do j=ndy,nfy
             duy= ( a7(1)*(uf(i,j+1,l)-uf(i,j-1,l)) &
                  + a7(2)*(uf(i,j+2,l)-uf(i,j-2,l)) &
                  + a7(3)*(uf(i,j+3,l)-uf(i,j-3,l)) ) *idy(j)*sintsinp(i,j,l)
             dvy= ( a7(1)*(vf(i,j+1,l)-vf(i,j-1,l)) &
                  + a7(2)*(vf(i,j+2,l)-vf(i,j-2,l)) &
                  + a7(3)*(vf(i,j+3,l)-vf(i,j-3,l)) ) *idy(j)*sintsinp(i,j,l)
             dwy= ( a7(1)*(wf(i,j+1,l)-wf(i,j-1,l)) &
                  + a7(2)*(wf(i,j+2,l)-wf(i,j-2,l)) &
                  + a7(3)*(wf(i,j+3,l)-wf(i,j-3,l)) ) *idy(j)*sintsinp(i,j,l)
             dpy= ( a7(1)*(pf(i,j+1,l)-pf(i,j-1,l)) &
                  + a7(2)*(pf(i,j+2,l)-pf(i,j-2,l)) &
                  + a7(3)*(pf(i,j+3,l)-pf(i,j-3,l)) ) *idy(j)*sintsinp(i,j,l)
             dry= ( a7(1)*(rf(i,j+1,l)-rf(i,j-1,l)) &
                  + a7(2)*(rf(i,j+2,l)-rf(i,j-2,l)) &
                  + a7(3)*(rf(i,j+3,l)-rf(i,j-3,l)) ) *idy(j)*sintsinp(i,j,l)

             pt(i,j,l) = vg(i,j,l)*(dpx(i,j,l)+dpy+dpz(i,j,l)+pf(i,j,l)*ir(i,j,l))
             ut(i,j,l) = vg(i,j,l)*(dux(i,j,l)+duy+duz(i,j,l)+uf(i,j,l)*ir(i,j,l))
             vt(i,j,l) = vg(i,j,l)*(dvx(i,j,l)+dvy+dvz(i,j,l)+vf(i,j,l)*ir(i,j,l))
             wt(i,j,l) = vg(i,j,l)*(dwx(i,j,l)+dwy+dwz(i,j,l)+wf(i,j,l)*ir(i,j,l))
             rt(i,j,l) = vg(i,j,l)*(drx(i,j,l)+dry+drz(i,j,l)+rf(i,j,l)*ir(i,j,l))
          enddo
       enddo
    enddo

    ! Update fluxes at each RK step
    ! =============================
    do j=ndy,nfy
       do k=nzmnghp1,nz
          l=k-nzmngh
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

  end subroutine bc_TD3d_kmax_imin

  !===============================================================================
  module subroutine bc_TD3d_kmax_imax
  !===============================================================================
    !> 3D Tam & Dong's BC: boundary condition at kmax-imax (edge 3,2,2 /back-right)
    !> - Cartesian version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l,m
    real(wp) :: cp,av,c2_
    real(wp), dimension(-2:ngh,ny1:ny2,-2:ngh) :: rf,uf,vf,wf,pf
    real(wp), dimension(1:ngh,ny,1:ngh) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(1:ngh,ny,1:ngh) :: dux,dvx,dwx,dpx,drx
    real(wp), dimension(1:ngh,ny,1:ngh) :: duz,dvz,dwz,dpz,drz
    real(wp) :: duy,dvy,dwy,dpy,dry
    !-------------------------------------------------------------------------
    
    ! Pointers for spherical coordinates
    ! ==================================
    ir=>BC_edge(3,2,2)%i_r
    cosphi=>BC_edge(3,2,2)%cosp
    sinphi=>BC_edge(3,2,2)%sinp
    costeta=>BC_edge(3,2,2)%cost
    sinteta=>BC_edge(3,2,2)%sint
    costcosp=>BC_edge(3,2,2)%costcosp
    costsinp=>BC_edge(3,2,2)%costsinp
    sintcosp=>BC_edge(3,2,2)%sintcosp
    sintsinp=>BC_edge(3,2,2)%sintsinp
    
    ! Compute fluctuations
    ! ====================
    do i=nxmngh-2,nx
       l=i-nxmngh
       do k=nzmngh-2,nz
          m=k-nzmngh
          do j=ny1,ny2
             rf(l,j,m)=rho_n(i,j,k)-BC_face(1,2)%U0(l,j,k,1)
             uf(l,j,m)=   uu(i,j,k)-BC_face(1,2)%U0(l,j,k,2)
             vf(l,j,m)=   vv(i,j,k)-BC_face(1,2)%U0(l,j,k,3)
             wf(l,j,m)=   ww(i,j,k)-BC_face(1,2)%U0(l,j,k,4)
             pf(l,j,m)=  prs(i,j,k)-BC_face(1,2)%U0(l,j,k,5)
          enddo
       enddo
    enddo
    
    ! Compute group velocity vg
    ! =========================
    ! vg= u0*sin(teta)*cos(phi)+v0*sin(teta)*sin(phi)+w0*cos(teta)
    !   + sqrt{c0^2-[u0*cos(teta)*cos(phi)+v0*cos(teta)*sin(phi)-w0*sin(teta)]^2-[u0*sin(phi)-v0*cos(phi)]^2}
    do l=1,ngh
       do k=nzmnghp1,nz
          m=k-nzmngh
          do j=ndy,nfy
             vg(l,j,m)= BC_face(1,2)%U0(l,j,k,2)*sintcosp(l,j,m) +     &
	           	BC_face(1,2)%U0(l,j,k,3)*sintsinp(l,j,m) +     &
	           	BC_face(1,2)%U0(l,j,k,4)*costeta(l,j,m)	       &
                      + sqrt( BC_face(1,2)%U0(l,j,k,6)-                & 				
                       ( BC_face(1,2)%U0(l,j,k,2)*costcosp(l,j,m) +    &				
                 	 BC_face(1,2)%U0(l,j,k,3)*costsinp(l,j,m) -    & 
          	         BC_face(1,2)%U0(l,j,k,4)*sinteta(l,j,m) )**2- &
	               ( BC_face(1,2)%U0(l,j,k,2)*sinphi(l,j,m) -      &
	      	         BC_face(1,2)%U0(l,j,k,3)*cosphi(l,j,m) )**2 )
          enddo
       enddo
    enddo

    ! Non-centered derivatives in x-direction *sin(teta)*cos(phi)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    do i=nxmnghp1,nx-3
       l=i-nxmngh
       do m=1,ngh
          do j=ndy,nfy
             dux(l,j,m) = ( a7(1)* (uf(l+1,j,m)-uf(l-1,j,m)) &
                          + a7(2)* (uf(l+2,j,m)-uf(l-2,j,m)) &
                          + a7(3)* (uf(l+3,j,m)-uf(l-3,j,m)) ) *idx(i)*sintcosp(l,j,m)
             dvx(l,j,m) = ( a7(1)* (vf(l+1,j,m)-vf(l-1,j,m)) &
                          + a7(2)* (vf(l+2,j,m)-vf(l-2,j,m)) &
                          + a7(3)* (vf(l+3,j,m)-vf(l-3,j,m)) ) *idx(i)*sintcosp(l,j,m)
             dwx(l,j,m) = ( a7(1)* (wf(l+1,j,m)-wf(l-1,j,m)) &
                          + a7(2)* (wf(l+2,j,m)-wf(l-2,j,m)) &
                          + a7(3)* (wf(l+3,j,m)-wf(l-3,j,m)) ) *idx(i)*sintcosp(l,j,m)
             dpx(l,j,m) = ( a7(1)* (pf(l+1,j,m)-pf(l-1,j,m)) &
                          + a7(2)* (pf(l+2,j,m)-pf(l-2,j,m)) &
                          + a7(3)* (pf(l+3,j,m)-pf(l-3,j,m)) ) *idx(i)*sintcosp(l,j,m)
             drx(l,j,m) = ( a7(1)* (rf(l+1,j,m)-rf(l-1,j,m)) &
                          + a7(2)* (rf(l+2,j,m)-rf(l-2,j,m)) &
                          + a7(3)* (rf(l+3,j,m)-rf(l-3,j,m)) ) *idx(i)*sintcosp(l,j,m)
          enddo
       enddo
    enddo
    
    i=nx-2
    l=i-nxmngh
    do m=1,ngh
       do j=ndy,nfy
          dux(l,j,m) = ( a42(1)*uf(l+2,j,m)+a42(2)*uf(l+1,j,m) &
                       + a42(3)*uf(l  ,j,m)+a42(4)*uf(l-1,j,m) &
                       + a42(5)*uf(l-2,j,m)+a42(6)*uf(l-3,j,m) &
                       + a42(7)*uf(l-4,j,m) ) *idx(i)*sintcosp(l,j,m)
          dvx(l,j,m) = ( a42(1)*vf(l+2,j,m)+a42(2)*vf(l+1,j,m) &
                       + a42(3)*vf(l  ,j,m)+a42(4)*vf(l-1,j,m) &
                       + a42(5)*vf(l-2,j,m)+a42(6)*vf(l-3,j,m) &
                       + a42(7)*vf(l-4,j,m) ) *idx(i)*sintcosp(l,j,m)
          dwx(l,j,m) = ( a42(1)*wf(l+2,j,m)+a42(2)*wf(l+1,j,m) &
                       + a42(3)*wf(l  ,j,m)+a42(4)*wf(l-1,j,m) &
                       + a42(5)*wf(l-2,j,m)+a42(6)*wf(l-3,j,m) &
                       + a42(7)*wf(l-4,j,m) ) *idx(i)*sintcosp(l,j,m)
          dpx(l,j,m) = ( a42(1)*pf(l+2,j,m)+a42(2)*pf(l+1,j,m) &
                       + a42(3)*pf(l  ,j,m)+a42(4)*pf(l-1,j,m) &
                       + a42(5)*pf(l-2,j,m)+a42(6)*pf(l-3,j,m) &
                       + a42(7)*pf(l-4,j,m) ) *idx(i)*sintcosp(l,j,m)
          drx(l,j,m) = ( a42(1)*rf(l+2,j,m)+a42(2)*rf(l+1,j,m) &
                       + a42(3)*rf(l  ,j,m)+a42(4)*rf(l-1,j,m) &
                       + a42(5)*rf(l-2,j,m)+a42(6)*rf(l-3,j,m) &
                       + a42(7)*rf(l-4,j,m) ) *idx(i)*sintcosp(l,j,m)
       enddo
    enddo
    
    i=nx-1
    l=i-nxmngh
    do m=1,ngh
       do j=ndy,nfy
          dux(l,j,m) = ( a51(1)*uf(l+1,j,m)+a51(2)*uf(l  ,j,m) &
                       + a51(3)*uf(l-1,j,m)+a51(4)*uf(l-2,j,m) &
                       + a51(5)*uf(l-3,j,m)+a51(6)*uf(l-4,j,m) &
                       + a51(7)*uf(l-5,j,m) ) *idx(i)*sintcosp(l,j,m)
          dvx(l,j,m) = ( a51(1)*vf(l+1,j,m)+a51(2)*vf(l  ,j,m) &
                       + a51(3)*vf(l-1,j,m)+a51(4)*vf(l-2,j,m) &
                       + a51(5)*vf(l-3,j,m)+a51(6)*vf(l-4,j,m) &
                       + a51(7)*vf(l-5,j,m) ) *idx(i)*sintcosp(l,j,m)
          dwx(l,j,m) = ( a51(1)*wf(l+1,j,m)+a51(2)*wf(l  ,j,m) &
                       + a51(3)*wf(l-1,j,m)+a51(4)*wf(l-2,j,m) &
                       + a51(5)*wf(l-3,j,m)+a51(6)*wf(l-4,j,m) &
                       + a51(7)*wf(l-5,j,m) ) *idx(i)*sintcosp(l,j,m)
          dpx(l,j,m) = ( a51(1)*pf(l+1,j,m)+a51(2)*pf(l  ,j,m) &
                       + a51(3)*pf(l-1,j,m)+a51(4)*pf(l-2,j,m) &
                       + a51(5)*pf(l-3,j,m)+a51(6)*pf(l-4,j,m) &
                       + a51(7)*pf(l-5,j,m) ) *idx(i)*sintcosp(l,j,m)
          drx(l,j,m) = ( a51(1)*rf(l+1,j,m)+a51(2)*rf(l  ,j,m) &
                       + a51(3)*rf(l-1,j,m)+a51(4)*rf(l-2,j,m) &
                       + a51(5)*rf(l-3,j,m)+a51(6)*rf(l-4,j,m) &
                       + a51(7)*rf(l-5,j,m) ) *idx(i)*sintcosp(l,j,m)
       enddo
    enddo
    
    i=nx
    l=i-nxmngh
    do m=1,ngh
       do j=ndy,nfy
          dux(l,j,m) = ( a60(1)*uf(l  ,j,m)+a60(2)*uf(l-1,j,m) &
                       + a60(3)*uf(l-2,j,m)+a60(4)*uf(l-3,j,m) &
                       + a60(5)*uf(l-4,j,m)+a60(6)*uf(l-5,j,m) &
                       + a60(7)*uf(l-6,j,m) ) *idx(i)*sintcosp(l,j,m)
          dvx(l,j,m) = ( a60(1)*vf(l  ,j,m)+a60(2)*vf(l-1,j,m) &
                       + a60(3)*vf(l-2,j,m)+a60(4)*vf(l-3,j,m) &
                       + a60(5)*vf(l-4,j,m)+a60(6)*vf(l-5,j,m) &
                       + a60(7)*vf(l-6,j,m) ) *idx(i)*sintcosp(l,j,m)
          dwx(l,j,m) = ( a60(1)*wf(l  ,j,m)+a60(2)*wf(l-1,j,m) &
                       + a60(3)*wf(l-2,j,m)+a60(4)*wf(l-3,j,m) &
                       + a60(5)*wf(l-4,j,m)+a60(6)*wf(l-5,j,m) &
                       + a60(7)*wf(l-6,j,m) ) *idx(i)*sintcosp(l,j,m)
          dpx(l,j,m) = ( a60(1)*pf(l  ,j,m)+a60(2)*pf(l-1,j,m) &
                       + a60(3)*pf(l-2,j,m)+a60(4)*pf(l-3,j,m) &
                       + a60(5)*pf(l-4,j,m)+a60(6)*pf(l-5,j,m) &
                       + a60(7)*pf(l-6,j,m) ) *idx(i)*sintcosp(l,j,m)
          drx(l,j,m) = ( a60(1)*rf(l  ,j,m)+a60(2)*rf(l-1,j,m) &
                       + a60(3)*rf(l-2,j,m)+a60(4)*rf(l-3,j,m) &
                       + a60(5)*rf(l-4,j,m)+a60(6)*rf(l-5,j,m) &
                       + a60(7)*rf(l-6,j,m) ) *idx(i)*sintcosp(l,j,m)
       enddo
    enddo
 
    ! Non-centered derivatives in z-direction *cos(teta)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    do k=nzmnghp1,nz-3
       m=k-nzmngh
       do l=1,ngh
          do j=ndy,nfy
             duz(l,j,m) = (a7(1)*(uf(l,j,m+1) - uf(l,j,m-1)) + &
                           a7(2)*(uf(l,j,m+2) - uf(l,j,m-2)) + &
                           a7(3)*(uf(l,j,m+3) - uf(l,j,m-3)) ) *idz(k)*costeta(l,j,m)
             dvz(l,j,m) = (a7(1)*(vf(l,j,m+1) - vf(l,j,m-1)) + &
                           a7(2)*(vf(l,j,m+2) - vf(l,j,m-2)) + &
                           a7(3)*(vf(l,j,m+3) - vf(l,j,m-3)) ) *idz(k)*costeta(l,j,m)
             dwz(l,j,m) = (a7(1)*(wf(l,j,m+1) - wf(l,j,m-1)) + &
                           a7(2)*(wf(l,j,m+2) - wf(l,j,m-2)) + &
                           a7(3)*(wf(l,j,m+3) - wf(l,j,m-3)) ) *idz(k)*costeta(l,j,m)
             dpz(l,j,m) = (a7(1)*(pf(l,j,m+1) - pf(l,j,m-1)) + &
                           a7(2)*(pf(l,j,m+2) - pf(l,j,m-2)) + &
                           a7(3)*(pf(l,j,m+3) - pf(l,j,m-3)) ) *idz(k)*costeta(l,j,m)
             drz(l,j,m) = (a7(1)*(rf(l,j,m+1) - rf(l,j,m-1)) + &
                           a7(2)*(rf(l,j,m+2) - rf(l,j,m-2)) + &
                           a7(3)*(rf(l,j,m+3) - rf(l,j,m-3)) ) *idz(k)*costeta(l,j,m)
          enddo
       enddo
    enddo

    k=nz-2
    m=k-nzmngh
    do l=1,ngh
       do j=ndy,nfy
          duz(l,j,m) = ( a42(1)*uf(l,j,m+2)+a42(2)*uf(l,j,m+1) &
                       + a42(3)*uf(l,j,m  )+a42(4)*uf(l,j,m-1) &
                       + a42(5)*uf(l,j,m-2)+a42(6)*uf(l,j,m-3) &
                       + a42(7)*uf(l,j,m-4) ) *idz(k)*costeta(l,j,m)
          dvz(l,j,m) = ( a42(1)*vf(l,j,m+2)+a42(2)*vf(l,j,m+1) &
                       + a42(3)*vf(l,j,m  )+a42(4)*vf(l,j,m-1) &
                       + a42(5)*vf(l,j,m-2)+a42(6)*vf(l,j,m-3) &
                       + a42(7)*vf(l,j,m-4) ) *idz(k)*costeta(l,j,m)
          dwz(l,j,m) = ( a42(1)*wf(l,j,m+2)+a42(2)*wf(l,j,m+1) &
                       + a42(3)*wf(l,j,m  )+a42(4)*wf(l,j,m-1) &
                       + a42(5)*wf(l,j,m-2)+a42(6)*wf(l,j,m-3) &
                       + a42(7)*wf(l,j,m-4) ) *idz(k)*costeta(l,j,m)
          dpz(l,j,m) = ( a42(1)*pf(l,j,m+2)+a42(2)*pf(l,j,m+1) &
                       + a42(3)*pf(l,j,m  )+a42(4)*pf(l,j,m-1) &
                       + a42(5)*pf(l,j,m-2)+a42(6)*pf(l,j,m-3) &
                       + a42(7)*pf(l,j,m-4) ) *idz(k)*costeta(l,j,m)
          drz(l,j,m) = ( a42(1)*rf(l,j,m+2)+a42(2)*rf(l,j,m+1) &
                       + a42(3)*rf(l,j,m  )+a42(4)*rf(l,j,m-1) &
                       + a42(5)*rf(l,j,m-2)+a42(6)*rf(l,j,m-3) &
                       + a42(7)*rf(l,j,m-4) ) *idz(k)*costeta(l,j,m)
       enddo
    enddo

    k=nz-1
    m=k-nzmngh
    do l=1,ngh
       do j=ndy,nfy
          duz(l,j,m) = ( a51(1)*uf(l,j,m+1)+a51(2)*uf(l,j,m  ) &
                       + a51(3)*uf(l,j,m-1)+a51(4)*uf(l,j,m-2) &
                       + a51(5)*uf(l,j,m-3)+a51(6)*uf(l,j,m-4) &
                       + a51(7)*uf(l,j,m-5) ) *idz(k)*costeta(l,j,m)
          dvz(l,j,m) = ( a51(1)*vf(l,j,m+1)+a51(2)*vf(l,j,m  ) &
                       + a51(3)*vf(l,j,m-1)+a51(4)*vf(l,j,m-2) &
                       + a51(5)*vf(l,j,m-3)+a51(6)*vf(l,j,m-4) &
                       + a51(7)*vf(l,j,m-5) ) *idz(k)*costeta(l,j,m)
          dwz(l,j,m) = ( a51(1)*wf(l,j,m+1)+a51(2)*wf(l,j,m  ) &
                       + a51(3)*wf(l,j,m-1)+a51(4)*wf(l,j,m-2) &
                       + a51(5)*wf(l,j,m-3)+a51(6)*wf(l,j,m-4) &
                       + a51(7)*wf(l,j,m-5) ) *idz(k)*costeta(l,j,m)
          dpz(l,j,m) = ( a51(1)*pf(l,j,m+1)+a51(2)*pf(l,j,m  ) &
                       + a51(3)*pf(l,j,m-1)+a51(4)*pf(l,j,m-2) &
                       + a51(5)*pf(l,j,m-3)+a51(6)*pf(l,j,m-4) &
                       + a51(7)*pf(l,j,m-5) ) *idz(k)*costeta(l,j,m)
          drz(l,j,m) = ( a51(1)*rf(l,j,m+1)+a51(2)*rf(l,j,m  ) &
                       + a51(3)*rf(l,j,m-1)+a51(4)*rf(l,j,m-2) &
                       + a51(5)*rf(l,j,m-3)+a51(6)*rf(l,j,m-4) &
                       + a51(7)*rf(l,j,m-5) ) *idz(k)*costeta(l,j,m)
       enddo
    enddo

    k=nz
    m=k-nzmngh
    do l=1,ngh
       do j=ndy,nfy
          duz(l,j,m) = ( a60(1)*uf(l,j,m  )+a60(2)*uf(l,j,m-1) &
                       + a60(3)*uf(l,j,m-2)+a60(4)*uf(l,j,m-3) &
                       + a60(5)*uf(l,j,m-4)+a60(6)*uf(l,j,m-5) &
                       + a60(7)*uf(l,j,m-6) ) *idz(k)*costeta(l,j,m)
          dvz(l,j,m) = ( a60(1)*vf(l,j,m  )+a60(2)*vf(l,j,m-1) &
                       + a60(3)*vf(l,j,m-2)+a60(4)*vf(l,j,m-3) &
                       + a60(5)*vf(l,j,m-4)+a60(6)*vf(l,j,m-5) &
                       + a60(7)*vf(l,j,m-6) ) *idz(k)*costeta(l,j,m)
          dwz(l,j,m) = ( a60(1)*wf(l,j,m  )+a60(2)*wf(l,j,m-1) &
                       + a60(3)*wf(l,j,m-2)+a60(4)*wf(l,j,m-3) &
                       + a60(5)*wf(l,j,m-4)+a60(6)*wf(l,j,m-5) &
                       + a60(7)*wf(l,j,m-6) ) *idz(k)*costeta(l,j,m)
          dpz(l,j,m) = ( a60(1)*pf(l,j,m  )+a60(2)*pf(l,j,m-1) &
                       + a60(3)*pf(l,j,m-2)+a60(4)*pf(l,j,m-3) &
                       + a60(5)*pf(l,j,m-4)+a60(6)*pf(l,j,m-5) &
                       + a60(7)*pf(l,j,m-6) ) *idz(k)*costeta(l,j,m)
          drz(l,j,m) = ( a60(1)*rf(l,j,m  )+a60(2)*rf(l,j,m-1) &
                       + a60(3)*rf(l,j,m-2)+a60(4)*rf(l,j,m-3) &
                       + a60(5)*rf(l,j,m-4)+a60(6)*rf(l,j,m-5) &
                       + a60(7)*rf(l,j,m-6) ) *idz(k)*costeta(l,j,m)
       enddo
    enddo

    ! Compute rt,ut,vt,wt & pt
    ! ========================
    do l=1,ngh
       do m=1,ngh
          do j=ndy,nfy
             duy= ( a7(1)*(uf(l,j+1,m)-uf(l,j-1,m)) &
                  + a7(2)*(uf(l,j+2,m)-uf(l,j-2,m)) &
                  + a7(3)*(uf(l,j+3,m)-uf(l,j-3,m)) ) *idy(j)*sintsinp(l,j,m)
             dvy= ( a7(1)*(vf(l,j+1,m)-vf(l,j-1,m)) &
                  + a7(2)*(vf(l,j+2,m)-vf(l,j-2,m)) &
                  + a7(3)*(vf(l,j+3,m)-vf(l,j-3,m)) ) *idy(j)*sintsinp(l,j,m)
             dwy= ( a7(1)*(wf(l,j+1,m)-wf(l,j-1,m)) &
                  + a7(2)*(wf(l,j+2,m)-wf(l,j-2,m)) &
                  + a7(3)*(wf(l,j+3,m)-wf(l,j-3,m)) ) *idy(j)*sintsinp(l,j,m)
             dpy= ( a7(1)*(pf(l,j+1,m)-pf(l,j-1,m)) &
                  + a7(2)*(pf(l,j+2,m)-pf(l,j-2,m)) &
                  + a7(3)*(pf(l,j+3,m)-pf(l,j-3,m)) ) *idy(j)*sintsinp(l,j,m)
             dry= ( a7(1)*(rf(l,j+1,m)-rf(l,j-1,m)) &
                  + a7(2)*(rf(l,j+2,m)-rf(l,j-2,m)) &
                  + a7(3)*(rf(l,j+3,m)-rf(l,j-3,m)) ) *idy(j)*sintsinp(l,j,m)

             pt(l,j,m) = vg(l,j,m)*(dpx(l,j,m)+dpy+dpz(l,j,m)+pf(l,j,m)*ir(l,j,m))
             ut(l,j,m) = vg(l,j,m)*(dux(l,j,m)+duy+duz(l,j,m)+uf(l,j,m)*ir(l,j,m))
             vt(l,j,m) = vg(l,j,m)*(dvx(l,j,m)+dvy+dvz(l,j,m)+vf(l,j,m)*ir(l,j,m))
             wt(l,j,m) = vg(l,j,m)*(dwx(l,j,m)+dwy+dwz(l,j,m)+wf(l,j,m)*ir(l,j,m))
             rt(l,j,m) = vg(l,j,m)*(drx(l,j,m)+dry+drz(l,j,m)+rf(l,j,m)*ir(l,j,m))
          enddo
       enddo
    enddo

    ! Update fluxes at each RK step
    ! =============================
    do j=ndy,nfy
       do i=nxmnghp1,nx
          l=i-nxmngh
          do k=nzmnghp1,nz
             m=k-nzmngh
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

  end subroutine bc_TD3d_kmax_imax

end submodule smod_TamDong3d_edgey
