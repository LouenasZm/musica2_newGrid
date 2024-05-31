!===============================================================================
submodule (mod_TamDong3d) smod_TamDong3d_corner
!===============================================================================
  !> author: XG
  !> date: February 2020 - modif January 2022
  !> Non-reflecting boundary conditions of Tam and Dong
  !> - 3D Cartesian version - Corners (edge intersection)
!=============================================================================== 

contains

  !===============================================================================
  module subroutine bc_TD3d_imin_jmin_kmin
  !===============================================================================
    !> 3D Tam & Dong's BC at imin-jmin-kmin (corner 1,1,1 /left-bottom-front)
    !> - Cartesian version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: cp,av,c2_
    real(wp), dimension(1:ngh+3,1:ngh+3,1:ngh+3) :: rf,uf,vf,wf,pf
    real(wp), dimension(ngh,ngh,ngh) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(ngh,ngh,ngh) :: dux,dvx,dwx,dpx,drx
    real(wp), dimension(ngh,ngh,ngh) :: duy,dvy,dwy,dpy,dry
    real(wp), dimension(ngh,ngh,ngh) :: duz,dvz,dwz,dpz,drz
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
    do j=1,ngh
       do k=1,ngh
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
    do j=1,ngh
       do k=1,ngh
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
    do j=1,ngh
       do k=1,ngh
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
       do j=1,ngh
          do k=1,ngh
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

    ! Non-centered derivatives in y-direction *sin(teta)*sin(phi)
    ! =======================================
    ! NOT(Tam & Webb DRP schemes) => order 2 (wall)
    j=1
    do i=1,ngh
       do k=1,ngh
          duy(i,j,k)=(uf(i,j+1,k)-uf(i,j,k))*idy1_jmin*sintsinp(i,j,k)
          dvy(i,j,k)=(vf(i,j+1,k)-vf(i,j,k))*idy1_jmin*sintsinp(i,j,k)
          dwy(i,j,k)=(wf(i,j+1,k)-wf(i,j,k))*idy1_jmin*sintsinp(i,j,k)
          dpy(i,j,k)=(pf(i,j+1,k)-pf(i,j,k))*idy1_jmin*sintsinp(i,j,k)
          dry(i,j,k)=(rf(i,j+1,k)-rf(i,j,k))*idy1_jmin*sintsinp(i,j,k)
       enddo
    enddo

    j=2
    do i=1,ngh
       do k=1,ngh
          duy(i,j,k)=a3(1)*(uf(i,j+1,k)-uf(i,j-1,k))*idy2_jmin*sintsinp(i,j,k)
          dvy(i,j,k)=a3(1)*(vf(i,j+1,k)-vf(i,j-1,k))*idy2_jmin*sintsinp(i,j,k)
          dwy(i,j,k)=a3(1)*(wf(i,j+1,k)-wf(i,j-1,k))*idy2_jmin*sintsinp(i,j,k)
          dpy(i,j,k)=a3(1)*(pf(i,j+1,k)-pf(i,j-1,k))*idy2_jmin*sintsinp(i,j,k)
          dry(i,j,k)=a3(1)*(rf(i,j+1,k)-rf(i,j-1,k))*idy2_jmin*sintsinp(i,j,k)
       enddo
    enddo

    j=3
    do i=1,ngh
       do k=1,ngh
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
!!$    do i=1,ngh
!!$       do k=1,ngh
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
!!$    do i=1,ngh
!!$       do k=1,ngh
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
!!$    do i=1,ngh
!!$       do k=1,ngh
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
       do i=1,ngh
          do k=1,ngh
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
    do i=1,ngh
       do j=1,ngh
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
       do j=1,ngh
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
       do j=1,ngh
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
          do j=1,ngh
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
    do i=1,ngh
       do j=1,ngh
          do k=1,ngh
             pt(i,j,k) = vg(i,j,k)*(dpx(i,j,k)+dpy(i,j,k)+dpz(i,j,k)+pf(i,j,k)*ir(i,j,k))
             ut(i,j,k) = vg(i,j,k)*(dux(i,j,k)+duy(i,j,k)+duz(i,j,k)+uf(i,j,k)*ir(i,j,k))
             vt(i,j,k) = vg(i,j,k)*(dvx(i,j,k)+dvy(i,j,k)+dvz(i,j,k)+vf(i,j,k)*ir(i,j,k))
             wt(i,j,k) = vg(i,j,k)*(dwx(i,j,k)+dwy(i,j,k)+dwz(i,j,k)+wf(i,j,k)*ir(i,j,k))
             rt(i,j,k) = vg(i,j,k)*(drx(i,j,k)+dry(i,j,k)+drz(i,j,k)+rf(i,j,k)*ir(i,j,k))
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

  end subroutine bc_TD3d_imin_jmin_kmin

  !===============================================================================
  module subroutine bc_TD3d_imin_jmin_kmax
  !===============================================================================
    !> 3D Tam & Dong's BC at imin-jmin-kmax (corner 1,1,2 /left-bottom-back)
    !> - Cartesian version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l
    real(wp) :: cp,av,c2_
    real(wp), dimension(1:ngh+3,1:ngh+3,-2:ngh) :: rf,uf,vf,wf,pf
    real(wp), dimension(ngh,ngh,ngh) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(ngh,ngh,ngh) :: dux,dvx,dwx,dpx,drx
    real(wp), dimension(ngh,ngh,ngh) :: duy,dvy,dwy,dpy,dry
    real(wp), dimension(ngh,ngh,ngh) :: duz,dvz,dwz,dpz,drz
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
	           	BC_face(1,1)%U0(i,j,k,4)*costeta(i,j,l)	       &
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
    do j=1,ngh
       do l=1,ngh
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
    do j=1,ngh
       do l=1,ngh
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
       do j=1,ngh
          do l=1,ngh
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

    ! Non-centered derivatives in y-direction *sin(teta)*sin(phi)
    ! =======================================
    ! NOT(Tam & Webb DRP schemes) => order 2 (wall)
    j=1
    do i=1,ngh
       do l=1,ngh
          duy(i,j,l)=(uf(i,j+1,l)-uf(i,j,l))*idy1_jmin*sintsinp(i,j,l)
          dvy(i,j,l)=(vf(i,j+1,l)-vf(i,j,l))*idy1_jmin*sintsinp(i,j,l)
          dwy(i,j,l)=(wf(i,j+1,l)-wf(i,j,l))*idy1_jmin*sintsinp(i,j,l)
          dpy(i,j,l)=(pf(i,j+1,l)-pf(i,j,l))*idy1_jmin*sintsinp(i,j,l)
          dry(i,j,l)=(rf(i,j+1,l)-rf(i,j,l))*idy1_jmin*sintsinp(i,j,l)
       enddo
    enddo

    j=2
    do i=1,ngh
       do l=1,ngh
          duy(i,j,l)=a3(1)*(uf(i,j+1,l)-uf(i,j-1,l))*idy2_jmin*sintsinp(i,j,l)
          dvy(i,j,l)=a3(1)*(vf(i,j+1,l)-vf(i,j-1,l))*idy2_jmin*sintsinp(i,j,l)
          dwy(i,j,l)=a3(1)*(wf(i,j+1,l)-wf(i,j-1,l))*idy2_jmin*sintsinp(i,j,l)
          dpy(i,j,l)=a3(1)*(pf(i,j+1,l)-pf(i,j-1,l))*idy2_jmin*sintsinp(i,j,l)
          dry(i,j,l)=a3(1)*(rf(i,j+1,l)-rf(i,j-1,l))*idy2_jmin*sintsinp(i,j,l)
       enddo
    enddo

    j=3
    do i=1,ngh
       do l=1,ngh
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
!!$    do i=1,ngh
!!$       do l=1,ngh
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
!!$    do i=1,ngh
!!$       do l=1,ngh
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
!!$    do i=1,ngh
!!$       do l=1,ngh
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
       do i=1,ngh
          do l=1,ngh
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
       do i=1,ngh
          do j=1,ngh
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
       do j=1,ngh
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
       do j=1,ngh
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
       do j=1,ngh
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
    do i=1,ngh
       do j=1,ngh
          do l=1,ngh
             pt(i,j,l) = vg(i,j,l)*(dpx(i,j,l)+dpy(i,j,l)+dpz(i,j,l)+pf(i,j,l)*ir(i,j,l))
             ut(i,j,l) = vg(i,j,l)*(dux(i,j,l)+duy(i,j,l)+duz(i,j,l)+uf(i,j,l)*ir(i,j,l))
             vt(i,j,l) = vg(i,j,l)*(dvx(i,j,l)+dvy(i,j,l)+dvz(i,j,l)+vf(i,j,l)*ir(i,j,l))
             wt(i,j,l) = vg(i,j,l)*(dwx(i,j,l)+dwy(i,j,l)+dwz(i,j,l)+wf(i,j,l)*ir(i,j,l))
             rt(i,j,l) = vg(i,j,l)*(drx(i,j,l)+dry(i,j,l)+drz(i,j,l)+rf(i,j,l)*ir(i,j,l))
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

  end subroutine bc_TD3d_imin_jmin_kmax

  !===============================================================================
  module subroutine bc_TD3d_imin_jmax_kmin
  !===============================================================================
    !> 3D Tam & Dong's BC at imin-jmax-kmin (corner 1,2,1 /left-top-front)
    !> - Cartesian version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l
    real(wp) :: cp,av,c2_
    real(wp), dimension(1:ngh+3,-2:ngh,1:ngh+3) :: rf,uf,vf,wf,pf
    real(wp), dimension(ngh,ngh,ngh) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(ngh,ngh,ngh) :: dux,dvx,dwx,dpx,drx
    real(wp), dimension(ngh,ngh,ngh) :: duy,dvy,dwy,dpy,dry
    real(wp), dimension(ngh,ngh,ngh) :: duz,dvz,dwz,dpz,drz
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
	           	BC_face(1,1)%U0(i,j,k,4)*costeta(i,l,k)	       &
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
          dux(i,l,k)= ( a06(1)*uf(1,l,k)+a06(2)*uf(2,l,k) &
                      + a06(3)*uf(3,l,k)+a06(4)*uf(4,l,k) &
                      + a06(5)*uf(5,l,k)+a06(6)*uf(6,l,k) &
                      + a06(7)*uf(7,l,k) ) *idx(i)*sintcosp(i,l,k)
          dvx(i,l,k)= ( a06(1)*vf(1,l,k)+a06(2)*vf(2,l,k) &
                      + a06(3)*vf(3,l,k)+a06(4)*vf(4,l,k) &
                      + a06(5)*vf(5,l,k)+a06(6)*vf(6,l,k) &
                      + a06(7)*vf(7,l,k) ) *idx(i)*sintcosp(i,l,k)
          dwx(i,l,k)= ( a06(1)*wf(1,l,k)+a06(2)*wf(2,l,k) &
                      + a06(3)*wf(3,l,k)+a06(4)*wf(4,l,k) &
                      + a06(5)*wf(5,l,k)+a06(6)*wf(6,l,k) &
                      + a06(7)*wf(7,l,k) ) *idx(i)*sintcosp(i,l,k)
          dpx(i,l,k)= ( a06(1)*pf(1,l,k)+a06(2)*pf(2,l,k) &
                      + a06(3)*pf(3,l,k)+a06(4)*pf(4,l,k) &
                      + a06(5)*pf(5,l,k)+a06(6)*pf(6,l,k) &
                      + a06(7)*pf(7,l,k) ) *idx(i)*sintcosp(i,l,k)
          drx(i,l,k)= ( a06(1)*rf(1,l,k)+a06(2)*rf(2,l,k) &
                      + a06(3)*rf(3,l,k)+a06(4)*rf(4,l,k) &
                      + a06(5)*rf(5,l,k)+a06(6)*rf(6,l,k) &
                      + a06(7)*rf(7,l,k) ) *idx(i)*sintcosp(i,l,k)
       enddo
    enddo

    i=2
    do j=nymnghp1,ny
       l=j-nymngh
       do k=1,ngh
          dux(i,l,k)= ( a15(1)*uf(1,l,k)+a15(2)*uf(2,l,k) &
                      + a15(3)*uf(3,l,k)+a15(4)*uf(4,l,k) &
                      + a15(5)*uf(5,l,k)+a15(6)*uf(6,l,k) &
                      + a15(7)*uf(7,l,k) ) *idx(i)*sintcosp(i,l,k)
          dvx(i,l,k)= ( a15(1)*vf(1,l,k)+a15(2)*vf(2,l,k) &
                      + a15(3)*vf(3,l,k)+a15(4)*vf(4,l,k) &
                      + a15(5)*vf(5,l,k)+a15(6)*vf(6,l,k) &
                      + a15(7)*vf(7,l,k) ) *idx(i)*sintcosp(i,l,k)
          dwx(i,l,k)= ( a15(1)*wf(1,l,k)+a15(2)*wf(2,l,k) &
                      + a15(3)*wf(3,l,k)+a15(4)*wf(4,l,k) &
                      + a15(5)*wf(5,l,k)+a15(6)*wf(6,l,k) &
                      + a15(7)*wf(7,l,k) ) *idx(i)*sintcosp(i,l,k)
          dpx(i,l,k)= ( a15(1)*pf(1,l,k)+a15(2)*pf(2,l,k) &
                      + a15(3)*pf(3,l,k)+a15(4)*pf(4,l,k) &
                      + a15(5)*pf(5,l,k)+a15(6)*pf(6,l,k) &
                      + a15(7)*pf(7,l,k) ) *idx(i)*sintcosp(i,l,k)
          drx(i,l,k)= ( a15(1)*rf(1,l,k)+a15(2)*rf(2,l,k) &
                      + a15(3)*rf(3,l,k)+a15(4)*rf(4,l,k) &
                      + a15(5)*rf(5,l,k)+a15(6)*rf(6,l,k) &
                      + a15(7)*rf(7,l,k) ) *idx(i)*sintcosp(i,l,k)
       enddo
    enddo

    i=3
    do j=nymnghp1,ny
       l=j-nymngh
       do k=1,ngh
          dux(i,l,k)= ( a24(1)*uf(1,l,k)+a24(2)*uf(2,l,k) &
                      + a24(3)*uf(3,l,k)+a24(4)*uf(4,l,k) &
                      + a24(5)*uf(5,l,k)+a24(6)*uf(6,l,k) &
                      + a24(7)*uf(7,l,k) ) *idx(i)*sintcosp(i,l,k)
          dvx(i,l,k)= ( a24(1)*vf(1,l,k)+a24(2)*vf(2,l,k) &
                      + a24(3)*vf(3,l,k)+a24(4)*vf(4,l,k) &
                      + a24(5)*vf(5,l,k)+a24(6)*vf(6,l,k) &
                      + a24(7)*vf(7,l,k) ) *idx(i)*sintcosp(i,l,k)
          dwx(i,l,k)= ( a24(1)*wf(1,l,k)+a24(2)*wf(2,l,k) &
                      + a24(3)*wf(3,l,k)+a24(4)*wf(4,l,k) &
                      + a24(5)*wf(5,l,k)+a24(6)*wf(6,l,k) &
                      + a24(7)*wf(7,l,k) ) *idx(i)*sintcosp(i,l,k)
          dpx(i,l,k)= ( a24(1)*pf(1,l,k)+a24(2)*pf(2,l,k) &
                      + a24(3)*pf(3,l,k)+a24(4)*pf(4,l,k) &
                      + a24(5)*pf(5,l,k)+a24(6)*pf(6,l,k) &
                      + a24(7)*pf(7,l,k) ) *idx(i)*sintcosp(i,l,k)
          drx(i,l,k)= ( a24(1)*rf(1,l,k)+a24(2)*rf(2,l,k) &
                      + a24(3)*rf(3,l,k)+a24(4)*rf(4,l,k) &
                      + a24(5)*rf(5,l,k)+a24(6)*rf(6,l,k) &
                      + a24(7)*rf(7,l,k) ) *idx(i)*sintcosp(i,l,k)
       enddo
    enddo

    do i=4,ngh
       do j=nymnghp1,ny
          l=j-nymngh
          do k=1,ngh
             dux(i,l,k)= ( a7(1)* (uf(i+1,l,k)-uf(i-1,l,k)) &
                         + a7(2)* (uf(i+2,l,k)-uf(i-2,l,k)) &
                         + a7(3)* (uf(i+3,l,k)-uf(i-3,l,k)) ) *idx(i)*sintcosp(i,l,k)
             dvx(i,l,k)= ( a7(1)* (vf(i+1,l,k)-vf(i-1,l,k)) &
                         + a7(2)* (vf(i+2,l,k)-vf(i-2,l,k)) &
                         + a7(3)* (vf(i+3,l,k)-vf(i-3,l,k)) ) *idx(i)*sintcosp(i,l,k)
             dwx(i,l,k)= ( a7(1)* (wf(i+1,l,k)-wf(i-1,l,k)) &
                         + a7(2)* (wf(i+2,l,k)-wf(i-2,l,k)) &
                         + a7(3)* (wf(i+3,l,k)-wf(i-3,l,k)) ) *idx(i)*sintcosp(i,l,k)
             dpx(i,l,k)= ( a7(1)* (pf(i+1,l,k)-pf(i-1,l,k)) &
                         + a7(2)* (pf(i+2,l,k)-pf(i-2,l,k)) &
                         + a7(3)* (pf(i+3,l,k)-pf(i-3,l,k)) ) *idx(i)*sintcosp(i,l,k)
             drx(i,l,k)= ( a7(1)* (rf(i+1,l,k)-rf(i-1,l,k)) &
                         + a7(2)* (rf(i+2,l,k)-rf(i-2,l,k)) &
                         + a7(3)* (rf(i+3,l,k)-rf(i-3,l,k)) ) *idx(i)*sintcosp(i,l,k)
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
    do i=1,ngh
       do k=1,ngh
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
    do i=1,ngh
       do k=1,ngh
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
    do i=1,ngh
       do k=1,ngh
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
    do i=1,ngh
       do l=1,ngh
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
    do i=1,ngh
       do l=1,ngh
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
    do i=1,ngh
       do l=1,ngh
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
       do i=1,ngh
          do l=1,ngh
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
    do i=1,ngh
       do l=1,ngh
          do k=1,ngh
             pt(i,l,k) = vg(i,l,k)*(dpx(i,l,k)+dpy(i,l,k)+dpz(i,l,k)+pf(i,l,k)*ir(i,l,k))
             ut(i,l,k) = vg(i,l,k)*(dux(i,l,k)+duy(i,l,k)+duz(i,l,k)+uf(i,l,k)*ir(i,l,k))
             vt(i,l,k) = vg(i,l,k)*(dvx(i,l,k)+dvy(i,l,k)+dvz(i,l,k)+vf(i,l,k)*ir(i,l,k))
             wt(i,l,k) = vg(i,l,k)*(dwx(i,l,k)+dwy(i,l,k)+dwz(i,l,k)+wf(i,l,k)*ir(i,l,k))
             rt(i,l,k) = vg(i,l,k)*(drx(i,l,k)+dry(i,l,k)+drz(i,l,k)+rf(i,l,k)*ir(i,l,k))
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

  end subroutine bc_TD3d_imin_jmax_kmin

  !===============================================================================
  module subroutine bc_TD3d_imin_jmax_kmax
  !===============================================================================
    !> 3D Tam & Dong's BC at imin-jmax-kmax (corner 1,2,2 /left-top-back)
    !> - Cartesian version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l,m
    real(wp) :: cp,av,c2_
    real(wp), dimension(1:ngh+3,-2:ngh,-2:ngh) :: rf,uf,vf,wf,pf
    real(wp), dimension(ngh,ngh,ngh) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(ngh,ngh,ngh) :: dux,dvx,dwx,dpx,drx
    real(wp), dimension(ngh,ngh,ngh) :: duy,dvy,dwy,dpy,dry
    real(wp), dimension(ngh,ngh,ngh) :: duz,dvz,dwz,dpz,drz
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
	           	BC_face(1,1)%U0(i,j,k,4)*costeta(i,l,m)	       &
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
          dux(i,l,m)= ( a06(1)*uf(1,l,m)+a06(2)*uf(2,l,m) &
                      + a06(3)*uf(3,l,m)+a06(4)*uf(4,l,m) &
                      + a06(5)*uf(5,l,m)+a06(6)*uf(6,l,m) &
                      + a06(7)*uf(7,l,m) ) *idx(i)*sintcosp(i,l,m)
          dvx(i,l,m)= ( a06(1)*vf(1,l,m)+a06(2)*vf(2,l,m) &
                      + a06(3)*vf(3,l,m)+a06(4)*vf(4,l,m) &
                      + a06(5)*vf(5,l,m)+a06(6)*vf(6,l,m) &
                      + a06(7)*vf(7,l,m) ) *idx(i)*sintcosp(i,l,m)
          dwx(i,l,m)= ( a06(1)*wf(1,l,m)+a06(2)*wf(2,l,m) &
                      + a06(3)*wf(3,l,m)+a06(4)*wf(4,l,m) &
                      + a06(5)*wf(5,l,m)+a06(6)*wf(6,l,m) &
                      + a06(7)*wf(7,l,m) ) *idx(i)*sintcosp(i,l,m)
          dpx(i,l,m)= ( a06(1)*pf(1,l,m)+a06(2)*pf(2,l,m) &
                      + a06(3)*pf(3,l,m)+a06(4)*pf(4,l,m) &
                      + a06(5)*pf(5,l,m)+a06(6)*pf(6,l,m) &
                      + a06(7)*pf(7,l,m) ) *idx(i)*sintcosp(i,l,m)
          drx(i,l,m)= ( a06(1)*rf(1,l,m)+a06(2)*rf(2,l,m) &
                      + a06(3)*rf(3,l,m)+a06(4)*rf(4,l,m) &
                      + a06(5)*rf(5,l,m)+a06(6)*rf(6,l,m) &
                      + a06(7)*rf(7,l,m) ) *idx(i)*sintcosp(i,l,m)
       enddo
    enddo

    i=2
    do j=nymnghp1,ny
       l=j-nymngh
       do m=1,ngh
          dux(i,l,m)= ( a15(1)*uf(1,l,m)+a15(2)*uf(2,l,m) &
                      + a15(3)*uf(3,l,m)+a15(4)*uf(4,l,m) &
                      + a15(5)*uf(5,l,m)+a15(6)*uf(6,l,m) &
                      + a15(7)*uf(7,l,m) ) *idx(i)*sintcosp(i,l,m)
          dvx(i,l,m)= ( a15(1)*vf(1,l,m)+a15(2)*vf(2,l,m) &
                      + a15(3)*vf(3,l,m)+a15(4)*vf(4,l,m) &
                      + a15(5)*vf(5,l,m)+a15(6)*vf(6,l,m) &
                      + a15(7)*vf(7,l,m) ) *idx(i)*sintcosp(i,l,m)
          dwx(i,l,m)= ( a15(1)*wf(1,l,m)+a15(2)*wf(2,l,m) &
                      + a15(3)*wf(3,l,m)+a15(4)*wf(4,l,m) &
                      + a15(5)*wf(5,l,m)+a15(6)*wf(6,l,m) &
                      + a15(7)*wf(7,l,m) ) *idx(i)*sintcosp(i,l,m)
          dpx(i,l,m)= ( a15(1)*pf(1,l,m)+a15(2)*pf(2,l,m) &
                      + a15(3)*pf(3,l,m)+a15(4)*pf(4,l,m) &
                      + a15(5)*pf(5,l,m)+a15(6)*pf(6,l,m) &
                      + a15(7)*pf(7,l,m) ) *idx(i)*sintcosp(i,l,m)
          drx(i,l,m)= ( a15(1)*rf(1,l,m)+a15(2)*rf(2,l,m) &
                      + a15(3)*rf(3,l,m)+a15(4)*rf(4,l,m) &
                      + a15(5)*rf(5,l,m)+a15(6)*rf(6,l,m) &
                      + a15(7)*rf(7,l,m) ) *idx(i)*sintcosp(i,l,m)
       enddo
    enddo

    i=3
    do j=nymnghp1,ny
       l=j-nymngh
       do m=1,ngh
          dux(i,l,m)= ( a24(1)*uf(1,l,m)+a24(2)*uf(2,l,m) &
                      + a24(3)*uf(3,l,m)+a24(4)*uf(4,l,m) &
                      + a24(5)*uf(5,l,m)+a24(6)*uf(6,l,m) &
                      + a24(7)*uf(7,l,m) ) *idx(i)*sintcosp(i,l,m)
          dvx(i,l,m)= ( a24(1)*vf(1,l,m)+a24(2)*vf(2,l,m) &
                      + a24(3)*vf(3,l,m)+a24(4)*vf(4,l,m) &
                      + a24(5)*vf(5,l,m)+a24(6)*vf(6,l,m) &
                      + a24(7)*vf(7,l,m) ) *idx(i)*sintcosp(i,l,m)
          dwx(i,l,m)= ( a24(1)*wf(1,l,m)+a24(2)*wf(2,l,m) &
                      + a24(3)*wf(3,l,m)+a24(4)*wf(4,l,m) &
                      + a24(5)*wf(5,l,m)+a24(6)*wf(6,l,m) &
                      + a24(7)*wf(7,l,m) ) *idx(i)*sintcosp(i,l,m)
          dpx(i,l,m)= ( a24(1)*pf(1,l,m)+a24(2)*pf(2,l,m) &
                      + a24(3)*pf(3,l,m)+a24(4)*pf(4,l,m) &
                      + a24(5)*pf(5,l,m)+a24(6)*pf(6,l,m) &
                      + a24(7)*pf(7,l,m) ) *idx(i)*sintcosp(i,l,m)
          drx(i,l,m)= ( a24(1)*rf(1,l,m)+a24(2)*rf(2,l,m) &
                      + a24(3)*rf(3,l,m)+a24(4)*rf(4,l,m) &
                      + a24(5)*rf(5,l,m)+a24(6)*rf(6,l,m) &
                      + a24(7)*rf(7,l,m) ) *idx(i)*sintcosp(i,l,m)
       enddo
    enddo

    do i=4,ngh
       do j=nymnghp1,ny
          l=j-nymngh
          do m=1,ngh
             dux(i,l,m)= ( a7(1)* (uf(i+1,l,m)-uf(i-1,l,m)) &
                         + a7(2)* (uf(i+2,l,m)-uf(i-2,l,m)) &
                         + a7(3)* (uf(i+3,l,m)-uf(i-3,l,m)) ) *idx(i)*sintcosp(i,l,m)
             dvx(i,l,m)= ( a7(1)* (vf(i+1,l,m)-vf(i-1,l,m)) &
                         + a7(2)* (vf(i+2,l,m)-vf(i-2,l,m)) &
                         + a7(3)* (vf(i+3,l,m)-vf(i-3,l,m)) ) *idx(i)*sintcosp(i,l,m)
             dwx(i,l,m)= ( a7(1)* (wf(i+1,l,m)-wf(i-1,l,m)) &
                         + a7(2)* (wf(i+2,l,m)-wf(i-2,l,m)) &
                         + a7(3)* (wf(i+3,l,m)-wf(i-3,l,m)) ) *idx(i)*sintcosp(i,l,m)
             dpx(i,l,m)= ( a7(1)* (pf(i+1,l,m)-pf(i-1,l,m)) &
                         + a7(2)* (pf(i+2,l,m)-pf(i-2,l,m)) &
                         + a7(3)* (pf(i+3,l,m)-pf(i-3,l,m)) ) *idx(i)*sintcosp(i,l,m)
             drx(i,l,m)= ( a7(1)* (rf(i+1,l,m)-rf(i-1,l,m)) &
                         + a7(2)* (rf(i+2,l,m)-rf(i-2,l,m)) &
                         + a7(3)* (rf(i+3,l,m)-rf(i-3,l,m)) ) *idx(i)*sintcosp(i,l,m)
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
    do i=1,ngh
       do m=1,ngh
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
    do i=1,ngh
       do m=1,ngh
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
    do i=1,ngh
       do m=1,ngh
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
       do i=1,ngh
          do l=1,ngh
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
    do i=1,ngh
       do l=1,ngh
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
    do i=1,ngh
       do l=1,ngh
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
    do i=1,ngh
       do l=1,ngh
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
    do i=1,ngh
       do l=1,ngh
          do m=1,ngh
             pt(i,l,m) = vg(i,l,m)*(dpx(i,l,m)+dpy(i,l,m)+dpz(i,l,m)+pf(i,l,m)*ir(i,l,m))
             ut(i,l,m) = vg(i,l,m)*(dux(i,l,m)+duy(i,l,m)+duz(i,l,m)+uf(i,l,m)*ir(i,l,m))
             vt(i,l,m) = vg(i,l,m)*(dvx(i,l,m)+dvy(i,l,m)+dvz(i,l,m)+vf(i,l,m)*ir(i,l,m))
             wt(i,l,m) = vg(i,l,m)*(dwx(i,l,m)+dwy(i,l,m)+dwz(i,l,m)+wf(i,l,m)*ir(i,l,m))
             rt(i,l,m) = vg(i,l,m)*(drx(i,l,m)+dry(i,l,m)+drz(i,l,m)+rf(i,l,m)*ir(i,l,m))
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

  end subroutine bc_TD3d_imin_jmax_kmax

  !===============================================================================
  module subroutine bc_TD3d_imax_jmin_kmin
  !===============================================================================
    !> 3D Tam & Dong's BC at imax-jmin-kmin (corner 2,1,1 /right-bottom-front)
    !> - Cartesian version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l
    real(wp) :: cp,av,c2_
    real(wp), dimension(-2:ngh,1:ngh+3,1:ngh+3) :: rf,uf,vf,wf,pf
    real(wp), dimension(ngh,ngh,ngh) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(ngh,ngh,ngh) :: dux,dvx,dwx,dpx,drx
    real(wp), dimension(ngh,ngh,ngh) :: duy,dvy,dwy,dpy,dry
    real(wp), dimension(ngh,ngh,ngh) :: duz,dvz,dwz,dpz,drz
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
	           	BC_face(2,1)%U0(i,j,k,4)*costeta(l,j,k)	       &
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
    do j=1,ngh
       do k=1,ngh
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
    do j=1,ngh
       do k=1,ngh
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
    do j=1,ngh
       do k=1,ngh
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

    ! Non-centered derivatives in y-direction *sin(teta)*sin(phi)
    ! =======================================
    ! NOT(Tam & Webb DRP schemes) => order 2 (wall)
    j=1
    do l=1,ngh
       do k=1,ngh
          duy(l,j,k)=(uf(l,j+1,k)-uf(l,j,k))*idy1_jmin*sintsinp(l,j,k)
          dvy(l,j,k)=(vf(l,j+1,k)-vf(l,j,k))*idy1_jmin*sintsinp(l,j,k)
          dwy(l,j,k)=(wf(l,j+1,k)-wf(l,j,k))*idy1_jmin*sintsinp(l,j,k)
          dpy(l,j,k)=(pf(l,j+1,k)-pf(l,j,k))*idy1_jmin*sintsinp(l,j,k)
          dry(l,j,k)=(rf(l,j+1,k)-rf(l,j,k))*idy1_jmin*sintsinp(l,j,k)
       enddo
    enddo

    j=2
    do l=1,ngh
       do k=1,ngh
          duy(l,j,k)=a3(1)*(uf(l,j+1,k)-uf(l,j-1,k))*idy2_jmin*sintsinp(l,j,k)
          dvy(l,j,k)=a3(1)*(vf(l,j+1,k)-vf(l,j-1,k))*idy2_jmin*sintsinp(l,j,k)
          dwy(l,j,k)=a3(1)*(wf(l,j+1,k)-wf(l,j-1,k))*idy2_jmin*sintsinp(l,j,k)
          dpy(l,j,k)=a3(1)*(pf(l,j+1,k)-pf(l,j-1,k))*idy2_jmin*sintsinp(l,j,k)
          dry(l,j,k)=a3(1)*(rf(l,j+1,k)-rf(l,j-1,k))*idy2_jmin*sintsinp(l,j,k)
       enddo
    enddo

    j=3
    do l=1,ngh
       do k=1,ngh
          duy(l,j,k)= ( a5(1)*(uf(l,j+1,k)-uf(l,j-1,k)) &
                      + a5(2)*(uf(l,j+2,k)-uf(l,j-2,k)) ) *idy4_jmin*sintsinp(l,j,k)
          dvy(l,j,k)= ( a5(1)*(vf(l,j+1,k)-vf(l,j-1,k)) &
                      + a5(2)*(vf(l,j+2,k)-vf(l,j-2,k)) ) *idy4_jmin*sintsinp(l,j,k)
          dwy(l,j,k)= ( a5(1)*(wf(l,j+1,k)-wf(l,j-1,k)) &
                      + a5(2)*(wf(l,j+2,k)-wf(l,j-2,k)) ) *idy4_jmin*sintsinp(l,j,k)
          dpy(l,j,k)= ( a5(1)*(pf(l,j+1,k)-pf(l,j-1,k)) &
                      + a5(2)*(pf(l,j+2,k)-pf(l,j-2,k)) ) *idy4_jmin*sintsinp(l,j,k)
          dry(l,j,k)= ( a5(1)*(rf(l,j+1,k)-rf(l,j-1,k)) &
                      + a5(2)*(rf(l,j+2,k)-rf(l,j-2,k)) ) *idy4_jmin*sintsinp(l,j,k)
       enddo
    enddo

    do j=4,ngh
       do l=1,ngh
          do k=1,ngh
             duy(l,j,k)= ( a7(1)*(uf(l,j+1,k)-uf(l,j-1,k)) &
                         + a7(2)*(uf(l,j+2,k)-uf(l,j-2,k)) &
                         + a7(3)*(uf(l,j+3,k)-uf(l,j-3,k)) ) *idy(j)*sintsinp(l,j,k)
             dvy(l,j,k)= ( a7(1)*(vf(l,j+1,k)-vf(l,j-1,k)) &
                         + a7(2)*(vf(l,j+2,k)-vf(l,j-2,k)) &
                         + a7(3)*(vf(l,j+3,k)-vf(l,j-3,k)) ) *idy(j)*sintsinp(l,j,k)
             dwy(l,j,k)= ( a7(1)*(wf(l,j+1,k)-wf(l,j-1,k)) &
                         + a7(2)*(wf(l,j+2,k)-wf(l,j-2,k)) &
                         + a7(3)*(wf(l,j+3,k)-wf(l,j-3,k)) ) *idy(j)*sintsinp(l,j,k)
             dpy(l,j,k)= ( a7(1)*(pf(l,j+1,k)-pf(l,j-1,k)) &
                         + a7(2)*(pf(l,j+2,k)-pf(l,j-2,k)) &
                         + a7(3)*(pf(l,j+3,k)-pf(l,j-3,k)) ) *idy(j)*sintsinp(l,j,k)
             dry(l,j,k)= ( a7(1)*(rf(l,j+1,k)-rf(l,j-1,k)) &
                         + a7(2)*(rf(l,j+2,k)-rf(l,j-2,k)) &
                         + a7(3)*(rf(l,j+3,k)-rf(l,j-3,k)) ) *idy(j)*sintsinp(l,j,k)
          enddo
       enddo
    enddo

    ! Non-centered derivatives in z-direction *cos(teta)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    k=1
    do l=1,ngh
       do j=1,ngh
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
       do j=1,ngh
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
       do j=1,ngh
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
          do j=1,ngh
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
    do l=1,ngh
       do j=1,ngh
          do k=1,ngh
             pt(l,j,k) = vg(l,j,k)*(dpx(l,j,k)+dpy(l,j,k)+dpz(l,j,k)+pf(l,j,k)*ir(l,j,k))
             ut(l,j,k) = vg(l,j,k)*(dux(l,j,k)+duy(l,j,k)+duz(l,j,k)+uf(l,j,k)*ir(l,j,k))
             vt(l,j,k) = vg(l,j,k)*(dvx(l,j,k)+dvy(l,j,k)+dvz(l,j,k)+vf(l,j,k)*ir(l,j,k))
             wt(l,j,k) = vg(l,j,k)*(dwx(l,j,k)+dwy(l,j,k)+dwz(l,j,k)+wf(l,j,k)*ir(l,j,k))
             rt(l,j,k) = vg(l,j,k)*(drx(l,j,k)+dry(l,j,k)+drz(l,j,k)+rf(l,j,k)*ir(l,j,k))
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

  end subroutine bc_TD3d_imax_jmin_kmin

  !===============================================================================
  module subroutine bc_TD3d_imax_jmin_kmax
  !===============================================================================
    !> 3D Tam & Dong's BC at imax-jmin-kmax (corner 2,1,2 /right-bottom-back)
    !> - Cartesian version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l,m
    real(wp) :: cp,av,c2_
    real(wp), dimension(-2:ngh,1:ngh+3,-2:ngh) :: rf,uf,vf,wf,pf
    real(wp), dimension(ngh,ngh,ngh) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(ngh,ngh,ngh) :: dux,dvx,dwx,dpx,drx
    real(wp), dimension(ngh,ngh,ngh) :: duy,dvy,dwy,dpy,dry
    real(wp), dimension(ngh,ngh,ngh) :: duz,dvz,dwz,dpz,drz
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
	           	BC_face(2,1)%U0(i,j,k,4)*costeta(l,j,m)	       &
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
    do j=1,ngh
       do m=1,ngh
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
    do j=1,ngh
       do m=1,ngh
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
    do j=1,ngh
       do m=1,ngh
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

    ! Non-centered derivatives in y-direction *sin(teta)*sin(phi)
    ! =======================================
    ! NOT(Tam & Webb DRP schemes) => order 2 (wall)
    j=1
    do i=nxmnghp1,nx
       l=i-nxmngh
       do m=1,ngh
          duy(l,j,m)=(uf(l,j+1,m)-uf(l,j,m))*idy1_jmin*sintsinp(l,j,m)
          dvy(l,j,m)=(vf(l,j+1,m)-vf(l,j,m))*idy1_jmin*sintsinp(l,j,m)
          dwy(l,j,m)=(wf(l,j+1,m)-wf(l,j,m))*idy1_jmin*sintsinp(l,j,m)
          dpy(l,j,m)=(pf(l,j+1,m)-pf(l,j,m))*idy1_jmin*sintsinp(l,j,m)
          dry(l,j,m)=(rf(l,j+1,m)-rf(l,j,m))*idy1_jmin*sintsinp(l,j,m)
       enddo
    enddo

    j=2
    do i=nxmnghp1,nx
       l=i-nxmngh
       do m=1,ngh
          duy(l,j,m)=a3(1)*(uf(l,j+1,m)-uf(l,j-1,m))*idy2_jmin*sintsinp(l,j,m)
          dvy(l,j,m)=a3(1)*(vf(l,j+1,m)-vf(l,j-1,m))*idy2_jmin*sintsinp(l,j,m)
          dwy(l,j,m)=a3(1)*(wf(l,j+1,m)-wf(l,j-1,m))*idy2_jmin*sintsinp(l,j,m)
          dpy(l,j,m)=a3(1)*(pf(l,j+1,m)-pf(l,j-1,m))*idy2_jmin*sintsinp(l,j,m)
          dry(l,j,m)=a3(1)*(rf(l,j+1,m)-rf(l,j-1,m))*idy2_jmin*sintsinp(l,j,m)
       enddo
    enddo

    j=3
    do i=nxmnghp1,nx
       l=i-nxmngh
       do m=1,ngh
          duy(l,j,m)= ( a5(1)*(uf(l,j+1,m)-uf(l,j-1,m)) &
                      + a5(2)*(uf(l,j+2,m)-uf(l,j-2,m)) ) *idy4_jmin*sintsinp(l,j,m)
          dvy(l,j,m)= ( a5(1)*(vf(l,j+1,m)-vf(l,j-1,m)) &
                      + a5(2)*(vf(l,j+2,m)-vf(l,j-2,m)) ) *idy4_jmin*sintsinp(l,j,m)
          dwy(l,j,m)= ( a5(1)*(wf(l,j+1,m)-wf(l,j-1,m)) &
                      + a5(2)*(wf(l,j+2,m)-wf(l,j-2,m)) ) *idy4_jmin*sintsinp(l,j,m)
          dpy(l,j,m)= ( a5(1)*(pf(l,j+1,m)-pf(l,j-1,m)) &
                      + a5(2)*(pf(l,j+2,m)-pf(l,j-2,m)) ) *idy4_jmin*sintsinp(l,j,m)
          dry(l,j,m)= ( a5(1)*(rf(l,j+1,m)-rf(l,j-1,m)) &
                      + a5(2)*(rf(l,j+2,m)-rf(l,j-2,m)) ) *idy4_jmin*sintsinp(l,j,m)
       enddo
    enddo

    do j=4,ngh
       do i=nxmnghp1,nx
          l=i-nxmngh
          do m=1,ngh
             duy(l,j,m)= ( a7(1)*(uf(l,j+1,m)-uf(l,j-1,m)) &
                         + a7(2)*(uf(l,j+2,m)-uf(l,j-2,m)) &
                         + a7(3)*(uf(l,j+3,m)-uf(l,j-3,m)) ) *idy(j)*sintsinp(l,j,m)
             dvy(l,j,m)= ( a7(1)*(vf(l,j+1,m)-vf(l,j-1,m)) &
                         + a7(2)*(vf(l,j+2,m)-vf(l,j-2,m)) &
                         + a7(3)*(vf(l,j+3,m)-vf(l,j-3,m)) ) *idy(j)*sintsinp(l,j,m)
             dwy(l,j,m)= ( a7(1)*(wf(l,j+1,m)-wf(l,j-1,m)) &
                         + a7(2)*(wf(l,j+2,m)-wf(l,j-2,m)) &
                         + a7(3)*(wf(l,j+3,m)-wf(l,j-3,m)) ) *idy(j)*sintsinp(l,j,m)
             dpy(l,j,m)= ( a7(1)*(pf(l,j+1,m)-pf(l,j-1,m)) &
                         + a7(2)*(pf(l,j+2,m)-pf(l,j-2,m)) &
                         + a7(3)*(pf(l,j+3,m)-pf(l,j-3,m)) ) *idy(j)*sintsinp(l,j,m)
             dry(l,j,m)= ( a7(1)*(rf(l,j+1,m)-rf(l,j-1,m)) &
                         + a7(2)*(rf(l,j+2,m)-rf(l,j-2,m)) &
                         + a7(3)*(rf(l,j+3,m)-rf(l,j-3,m)) ) *idy(j)*sintsinp(l,j,m)
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
       do j=1,ngh
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
       do j=1,ngh
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
       do j=1,ngh
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
       do j=1,ngh
          do m=1,ngh
             pt(l,j,m) = vg(l,j,m)*(dpx(l,j,m)+dpy(l,j,m)+dpz(l,j,m)+pf(l,j,m)*ir(l,j,m))
             ut(l,j,m) = vg(l,j,m)*(dux(l,j,m)+duy(l,j,m)+duz(l,j,m)+uf(l,j,m)*ir(l,j,m))
             vt(l,j,m) = vg(l,j,m)*(dvx(l,j,m)+dvy(l,j,m)+dvz(l,j,m)+vf(l,j,m)*ir(l,j,m))
             wt(l,j,m) = vg(l,j,m)*(dwx(l,j,m)+dwy(l,j,m)+dwz(l,j,m)+wf(l,j,m)*ir(l,j,m))
             rt(l,j,m) = vg(l,j,m)*(drx(l,j,m)+dry(l,j,m)+drz(l,j,m)+rf(l,j,m)*ir(l,j,m))
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

  end subroutine bc_TD3d_imax_jmin_kmax

  !===============================================================================
  module subroutine bc_TD3d_imax_jmax_kmin
  !===============================================================================
    !> 3D Tam & Dong's BC at imax-jmax-kmin (corner 2,2,1 /right-top-front)
    !> - Cartesian version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l,m
    real(wp) :: cp,av,c2_
    real(wp), dimension(-2:ngh,-2:ngh,1:ngh+3) :: rf,uf,vf,wf,pf
    real(wp), dimension(ngh,ngh,ngh) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(ngh,ngh,ngh) :: dux,dvx,dwx,dpx,drx
    real(wp), dimension(ngh,ngh,ngh) :: duy,dvy,dwy,dpy,dry
    real(wp), dimension(ngh,ngh,ngh) :: duz,dvz,dwz,dpz,drz
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
	           	BC_face(2,2)%U0(i,m,k,4)*costeta(l,m,k)	       &
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
             dux(l,m,k) = ( a7(1)* (uf(l+1,m,k)-uf(l-1,m,k)) &
                          + a7(2)* (uf(l+2,m,k)-uf(l-2,m,k)) &
                          + a7(3)* (uf(l+3,m,k)-uf(l-3,m,k)) ) *idx(i)*sintcosp(l,m,k)
             dvx(l,m,k) = ( a7(1)* (vf(l+1,m,k)-vf(l-1,m,k)) &
                          + a7(2)* (vf(l+2,m,k)-vf(l-2,m,k)) &
                          + a7(3)* (vf(l+3,m,k)-vf(l-3,m,k)) ) *idx(i)*sintcosp(l,m,k)
             dwx(l,m,k) = ( a7(1)* (wf(l+1,m,k)-wf(l-1,m,k)) &
                          + a7(2)* (wf(l+2,m,k)-wf(l-2,m,k)) &
                          + a7(3)* (wf(l+3,m,k)-wf(l-3,m,k)) ) *idx(i)*sintcosp(l,m,k)
             dpx(l,m,k) = ( a7(1)* (pf(l+1,m,k)-pf(l-1,m,k)) &
                          + a7(2)* (pf(l+2,m,k)-pf(l-2,m,k)) &
                          + a7(3)* (pf(l+3,m,k)-pf(l-3,m,k)) ) *idx(i)*sintcosp(l,m,k)
             drx(l,m,k) = ( a7(1)* (rf(l+1,m,k)-rf(l-1,m,k)) &
                          + a7(2)* (rf(l+2,m,k)-rf(l-2,m,k)) &
                          + a7(3)* (rf(l+3,m,k)-rf(l-3,m,k)) ) *idx(i)*sintcosp(l,m,k)
          enddo
       enddo
    enddo
    
    i=nx-2
    l=i-nxmngh
    do j=nymnghp1,ny
       m=j-nymngh
       do k=1,ngh
          dux(l,m,k) = ( a42(1)*uf(l+2,m,k)+a42(2)*uf(l+1,m,k) &
                       + a42(3)*uf(l  ,m,k)+a42(4)*uf(l-1,m,k) &
                       + a42(5)*uf(l-2,m,k)+a42(6)*uf(l-3,m,k) &
                       + a42(7)*uf(l-4,m,k) ) *idx(i)*sintcosp(l,m,k)
          dvx(l,m,k) = ( a42(1)*vf(l+2,m,k)+a42(2)*vf(l+1,m,k) &
                       + a42(3)*vf(l  ,m,k)+a42(4)*vf(l-1,m,k) &
                       + a42(5)*vf(l-2,m,k)+a42(6)*vf(l-3,m,k) &
                       + a42(7)*vf(l-4,m,k) ) *idx(i)*sintcosp(l,m,k)
          dwx(l,m,k) = ( a42(1)*wf(l+2,m,k)+a42(2)*wf(l+1,m,k) &
                       + a42(3)*wf(l  ,m,k)+a42(4)*wf(l-1,m,k) &
                       + a42(5)*wf(l-2,m,k)+a42(6)*wf(l-3,m,k) &
                       + a42(7)*wf(l-4,m,k) ) *idx(i)*sintcosp(l,m,k)
          dpx(l,m,k) = ( a42(1)*pf(l+2,m,k)+a42(2)*pf(l+1,m,k) &
                       + a42(3)*pf(l  ,m,k)+a42(4)*pf(l-1,m,k) &
                       + a42(5)*pf(l-2,m,k)+a42(6)*pf(l-3,m,k) &
                       + a42(7)*pf(l-4,m,k) ) *idx(i)*sintcosp(l,m,k)
          drx(l,m,k) = ( a42(1)*rf(l+2,m,k)+a42(2)*rf(l+1,m,k) &
                       + a42(3)*rf(l  ,m,k)+a42(4)*rf(l-1,m,k) &
                       + a42(5)*rf(l-2,m,k)+a42(6)*rf(l-3,m,k) &
                       + a42(7)*rf(l-4,m,k) ) *idx(i)*sintcosp(l,m,k)
       enddo
    enddo
    
    i=nx-1
    l=i-nxmngh
    do j=nymnghp1,ny
       m=j-nymngh
       do k=1,ngh
          dux(l,m,k) = ( a51(1)*uf(l+1,m,k)+a51(2)*uf(l  ,m,k) &
                       + a51(3)*uf(l-1,m,k)+a51(4)*uf(l-2,m,k) &
                       + a51(5)*uf(l-3,m,k)+a51(6)*uf(l-4,m,k) &
                       + a51(7)*uf(l-5,m,k) ) *idx(i)*sintcosp(l,m,k)
          dvx(l,m,k) = ( a51(1)*vf(l+1,m,k)+a51(2)*vf(l  ,m,k) &
                       + a51(3)*vf(l-1,m,k)+a51(4)*vf(l-2,m,k) &
                       + a51(5)*vf(l-3,m,k)+a51(6)*vf(l-4,m,k) &
                       + a51(7)*vf(l-5,m,k) ) *idx(i)*sintcosp(l,m,k)
          dwx(l,m,k) = ( a51(1)*wf(l+1,m,k)+a51(2)*wf(l  ,m,k) &
                       + a51(3)*wf(l-1,m,k)+a51(4)*wf(l-2,m,k) &
                       + a51(5)*wf(l-3,m,k)+a51(6)*wf(l-4,m,k) &
                       + a51(7)*wf(l-5,m,k) ) *idx(i)*sintcosp(l,m,k)
          dpx(l,m,k) = ( a51(1)*pf(l+1,m,k)+a51(2)*pf(l  ,m,k) &
                       + a51(3)*pf(l-1,m,k)+a51(4)*pf(l-2,m,k) &
                       + a51(5)*pf(l-3,m,k)+a51(6)*pf(l-4,m,k) &
                       + a51(7)*pf(l-5,m,k) ) *idx(i)*sintcosp(l,m,k)
          drx(l,m,k) = ( a51(1)*rf(l+1,m,k)+a51(2)*rf(l  ,m,k) &
                       + a51(3)*rf(l-1,m,k)+a51(4)*rf(l-2,m,k) &
                       + a51(5)*rf(l-3,m,k)+a51(6)*rf(l-4,m,k) &
                       + a51(7)*rf(l-5,m,k) ) *idx(i)*sintcosp(l,m,k)
       enddo
    enddo
    
    i=nx
    l=i-nxmngh
    do j=nymnghp1,ny
       m=j-nymngh
       do k=1,ngh
          dux(l,m,k) = ( a60(1)*uf(l  ,m,k)+a60(2)*uf(l-1,m,k) &
                       + a60(3)*uf(l-2,m,k)+a60(4)*uf(l-3,m,k) &
                       + a60(5)*uf(l-4,m,k)+a60(6)*uf(l-5,m,k) &
                       + a60(7)*uf(l-6,m,k) ) *idx(i)*sintcosp(l,m,k)
          dvx(l,m,k) = ( a60(1)*vf(l  ,m,k)+a60(2)*vf(l-1,m,k) &
                       + a60(3)*vf(l-2,m,k)+a60(4)*vf(l-3,m,k) &
                       + a60(5)*vf(l-4,m,k)+a60(6)*vf(l-5,m,k) &
                       + a60(7)*vf(l-6,m,k) ) *idx(i)*sintcosp(l,m,k)
          dwx(l,m,k) = ( a60(1)*wf(l  ,m,k)+a60(2)*wf(l-1,m,k) &
                       + a60(3)*wf(l-2,m,k)+a60(4)*wf(l-3,m,k) &
                       + a60(5)*wf(l-4,m,k)+a60(6)*wf(l-5,m,k) &
                       + a60(7)*wf(l-6,m,k) ) *idx(i)*sintcosp(l,m,k)
          dpx(l,m,k) = ( a60(1)*pf(l  ,m,k)+a60(2)*pf(l-1,m,k) &
                       + a60(3)*pf(l-2,m,k)+a60(4)*pf(l-3,m,k) &
                       + a60(5)*pf(l-4,m,k)+a60(6)*pf(l-5,m,k) &
                       + a60(7)*pf(l-6,m,k) ) *idx(i)*sintcosp(l,m,k)
          drx(l,m,k) = ( a60(1)*rf(l  ,m,k)+a60(2)*rf(l-1,m,k) &
                       + a60(3)*rf(l-2,m,k)+a60(4)*rf(l-3,m,k) &
                       + a60(5)*rf(l-4,m,k)+a60(6)*rf(l-5,m,k) &
                       + a60(7)*rf(l-6,m,k) ) *idx(i)*sintcosp(l,m,k)
       enddo
    enddo
 
    ! Non-centered derivatives in y-direction *sin(teta)*sin(phi)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    do j=nymnghp1,ny-3
       m=j-nymngh
       do l=1,ngh
          do k=1,ngh
             duy(l,m,k) = ( a7(1)*(uf(l,m+1,k)-uf(l,m-1,k)) &
                          + a7(2)*(uf(l,m+2,k)-uf(l,m-2,k)) &
                          + a7(3)*(uf(l,m+3,k)-uf(l,m-3,k)) ) *idy(j)*sintsinp(l,m,k)
             dvy(l,m,k) = ( a7(1)*(vf(l,m+1,k)-vf(l,m-1,k)) &
                          + a7(2)*(vf(l,m+2,k)-vf(l,m-2,k)) &
                          + a7(3)*(vf(l,m+3,k)-vf(l,m-3,k)) ) *idy(j)*sintsinp(l,m,k)
             dwy(l,m,k) = ( a7(1)*(wf(l,m+1,k)-wf(l,m-1,k)) &
                          + a7(2)*(wf(l,m+2,k)-wf(l,m-2,k)) &
                          + a7(3)*(wf(l,m+3,k)-wf(l,m-3,k)) ) *idy(j)*sintsinp(l,m,k)
             dpy(l,m,k) = ( a7(1)*(pf(l,m+1,k)-pf(l,m-1,k)) &
                          + a7(2)*(pf(l,m+2,k)-pf(l,m-2,k)) &
                          + a7(3)*(pf(l,m+3,k)-pf(l,m-3,k)) ) *idy(j)*sintsinp(l,m,k)
             dry(l,m,k) = ( a7(1)*(rf(l,m+1,k)-rf(l,m-1,k)) &
                          + a7(2)*(rf(l,m+2,k)-rf(l,m-2,k)) &
                          + a7(3)*(rf(l,m+3,k)-rf(l,m-3,k)) ) *idy(j)*sintsinp(l,m,k)
          enddo
       enddo
    enddo

    j=ny-2
    m=j-nymngh
    do l=1,ngh
       do k=1,ngh
          duy(l,m,k) = ( a42(1)*uf(l,m+2,k)+a42(2)*uf(l,m+1,k) &
                       + a42(3)*uf(l,m  ,k)+a42(4)*uf(l,m-1,k) &
                       + a42(5)*uf(l,m-2,k)+a42(6)*uf(l,m-3,k) &
                       + a42(7)*uf(l,m-4,k) ) *idy(j)*sintsinp(l,m,k)
          dvy(l,m,k) = ( a42(1)*vf(l,m+2,k)+a42(2)*vf(l,m+1,k) &
                       + a42(3)*vf(l,m  ,k)+a42(4)*vf(l,m-1,k) &
                       + a42(5)*vf(l,m-2,k)+a42(6)*vf(l,m-3,k) &
                       + a42(7)*vf(l,m-4,k) ) *idy(j)*sintsinp(l,m,k)
          dwy(l,m,k) = ( a42(1)*wf(l,m+2,k)+a42(2)*wf(l,m+1,k) &
                       + a42(3)*wf(l,m  ,k)+a42(4)*wf(l,m-1,k) &
                       + a42(5)*wf(l,m-2,k)+a42(6)*wf(l,m-3,k) &
                       + a42(7)*wf(l,m-4,k) ) *idy(j)*sintsinp(l,m,k)
          dpy(l,m,k) = ( a42(1)*pf(l,m+2,k)+a42(2)*pf(l,m+1,k) &
                       + a42(3)*pf(l,m  ,k)+a42(4)*pf(l,m-1,k) &
                       + a42(5)*pf(l,m-2,k)+a42(6)*pf(l,m-3,k) &
                       + a42(7)*pf(l,m-4,k) ) *idy(j)*sintsinp(l,m,k)
          dry(l,m,k) = ( a42(1)*rf(l,m+2,k)+a42(2)*rf(l,m+1,k) &
                       + a42(3)*rf(l,m  ,k)+a42(4)*rf(l,m-1,k) &
                       + a42(5)*rf(l,m-2,k)+a42(6)*rf(l,m-3,k) &
                       + a42(7)*rf(l,m-4,k) ) *idy(j)*sintsinp(l,m,k)
       enddo
    enddo
    
    j=ny-1
    m=j-nymngh
    do l=1,ngh
       do k=1,ngh
          duy(l,m,k) = ( a51(1)*uf(l,m+1,k)+a51(2)*uf(l,m  ,k) &
                       + a51(3)*uf(l,m-1,k)+a51(4)*uf(l,m-2,k) &
                       + a51(5)*uf(l,m-3,k)+a51(6)*uf(l,m-4,k) &
                       + a51(7)*uf(l,m-5,k) ) *idy(j)*sintsinp(l,m,k)
          dvy(l,m,k) = ( a51(1)*vf(l,m+1,k)+a51(2)*vf(l,m  ,k) &
                       + a51(3)*vf(l,m-1,k)+a51(4)*vf(l,m-2,k) &
                       + a51(5)*vf(l,m-3,k)+a51(6)*vf(l,m-4,k) &
                       + a51(7)*vf(l,m-5,k) ) *idy(j)*sintsinp(l,m,k)
          dwy(l,m,k) = ( a51(1)*wf(l,m+1,k)+a51(2)*wf(l,m  ,k) &
                       + a51(3)*wf(l,m-1,k)+a51(4)*wf(l,m-2,k) &
                       + a51(5)*wf(l,m-3,k)+a51(6)*wf(l,m-4,k) &
                       + a51(7)*wf(l,m-5,k) ) *idy(j)*sintsinp(l,m,k)
          dpy(l,m,k) = ( a51(1)*pf(l,m+1,k)+a51(2)*pf(l,m  ,k) &
                       + a51(3)*pf(l,m-1,k)+a51(4)*pf(l,m-2,k) &
                       + a51(5)*pf(l,m-3,k)+a51(6)*pf(l,m-4,k) &
                       + a51(7)*pf(l,m-5,k) ) *idy(j)*sintsinp(l,m,k)
          dry(l,m,k) = ( a51(1)*rf(l,m+1,k)+a51(2)*rf(l,m  ,k) &
                       + a51(3)*rf(l,m-1,k)+a51(4)*rf(l,m-2,k) &
                       + a51(5)*rf(l,m-3,k)+a51(6)*rf(l,m-4,k) &
                       + a51(7)*rf(l,m-5,k) ) *idy(j)*sintsinp(l,m,k)
       enddo
    enddo
    
    j=ny
    m=j-nymngh
    do l=1,ngh
       do k=1,ngh
          duy(l,m,k) = ( a60(1)*uf(l,m  ,k)+a60(2)*uf(l,m-1,k) &
                       + a60(3)*uf(l,m-2,k)+a60(4)*uf(l,m-3,k) &
                       + a60(5)*uf(l,m-4,k)+a60(6)*uf(l,m-5,k) &
                       + a60(7)*uf(l,m-6,k) ) *idy(j)*sintsinp(l,m,k)
          dvy(l,m,k) = ( a60(1)*vf(l,m  ,k)+a60(2)*vf(l,m-1,k) &
                       + a60(3)*vf(l,m-2,k)+a60(4)*vf(l,m-3,k) &
                       + a60(5)*vf(l,m-4,k)+a60(6)*vf(l,m-5,k) &
                       + a60(7)*vf(l,m-6,k) ) *idy(j)*sintsinp(l,m,k)
          dwy(l,m,k) = ( a60(1)*wf(l,m  ,k)+a60(2)*wf(l,m-1,k) &
                       + a60(3)*wf(l,m-2,k)+a60(4)*wf(l,m-3,k) &
                       + a60(5)*wf(l,m-4,k)+a60(6)*wf(l,m-5,k) &
                       + a60(7)*wf(l,m-6,k) ) *idy(j)*sintsinp(l,m,k)
          dpy(l,m,k) = ( a60(1)*pf(l,m  ,k)+a60(2)*pf(l,m-1,k) &
                       + a60(3)*pf(l,m-2,k)+a60(4)*pf(l,m-3,k) &
                       + a60(5)*pf(l,m-4,k)+a60(6)*pf(l,m-5,k) &
                       + a60(7)*pf(l,m-6,k) ) *idy(j)*sintsinp(l,m,k)
          dry(l,m,k) = ( a60(1)*rf(l,m  ,k)+a60(2)*rf(l,m-1,k) &
                       + a60(3)*rf(l,m-2,k)+a60(4)*rf(l,m-3,k) &
                       + a60(5)*rf(l,m-4,k)+a60(6)*rf(l,m-5,k) &
                       + a60(7)*rf(l,m-6,k) ) *idy(j)*sintsinp(l,m,k)
       enddo
    enddo
              
    ! Non-centered derivatives in z-direction *cos(teta)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    k=1
    do l=1,ngh
       do m=1,ngh
          duz(l,m,k) = ( a06(1)*uf(l,m,1)+a06(2)*uf(l,m,2) &
                       + a06(3)*uf(l,m,3)+a06(4)*uf(l,m,4) &
                       + a06(5)*uf(l,m,5)+a06(6)*uf(l,m,6) &
                       + a06(7)*uf(l,m,7) ) *idz(k)*costeta(l,m,k)
          dvz(l,m,k) = ( a06(1)*vf(l,m,1)+a06(2)*vf(l,m,2) &
                       + a06(3)*vf(l,m,3)+a06(4)*vf(l,m,4) &
                       + a06(5)*vf(l,m,5)+a06(6)*vf(l,m,6) &
                       + a06(7)*vf(l,m,7) ) *idz(k)*costeta(l,m,k)
          dwz(l,m,k) = ( a06(1)*wf(l,m,1)+a06(2)*wf(l,m,2) &
                       + a06(3)*wf(l,m,3)+a06(4)*wf(l,m,4) &
                       + a06(5)*wf(l,m,5)+a06(6)*wf(l,m,6) &
                       + a06(7)*wf(l,m,7) ) *idz(k)*costeta(l,m,k)
          dpz(l,m,k) = ( a06(1)*pf(l,m,1)+a06(2)*pf(l,m,2) &
                       + a06(3)*pf(l,m,3)+a06(4)*pf(l,m,4) &
                       + a06(5)*pf(l,m,5)+a06(6)*pf(l,m,6) &
                       + a06(7)*pf(l,m,7) ) *idz(k)*costeta(l,m,k)
          drz(l,m,k) = ( a06(1)*rf(l,m,1)+a06(2)*rf(l,m,2) &
                       + a06(3)*rf(l,m,3)+a06(4)*rf(l,m,4) &
                       + a06(5)*rf(l,m,5)+a06(6)*rf(l,m,6) &
                       + a06(7)*rf(l,m,7) ) *idz(k)*costeta(l,m,k)
       enddo
    enddo

    k=2
    do l=1,ngh
       do m=1,ngh
          duz(l,m,k) = ( a15(1)*uf(l,m,1)+a15(2)*uf(l,m,2) &
                       + a15(3)*uf(l,m,3)+a15(4)*uf(l,m,4) &
                       + a15(5)*uf(l,m,5)+a15(6)*uf(l,m,6) &
                       + a15(7)*uf(l,m,7) ) *idz(k)*costeta(l,m,k)
          dvz(l,m,k) = ( a15(1)*vf(l,m,1)+a15(2)*vf(l,m,2) &
                       + a15(3)*vf(l,m,3)+a15(4)*vf(l,m,4) &
                       + a15(5)*vf(l,m,5)+a15(6)*vf(l,m,6) &
                       + a15(7)*vf(l,m,7) ) *idz(k)*costeta(l,m,k)
          dwz(l,m,k) = ( a15(1)*wf(l,m,1)+a15(2)*wf(l,m,2) &
                       + a15(3)*wf(l,m,3)+a15(4)*wf(l,m,4) &
                       + a15(5)*wf(l,m,5)+a15(6)*wf(l,m,6) &
                       + a15(7)*wf(l,m,7) ) *idz(k)*costeta(l,m,k)
          dpz(l,m,k) = ( a15(1)*pf(l,m,1)+a15(2)*pf(l,m,2) &
                       + a15(3)*pf(l,m,3)+a15(4)*pf(l,m,4) &
                       + a15(5)*pf(l,m,5)+a15(6)*pf(l,m,6) &
                       + a15(7)*pf(l,m,7) ) *idz(k)*costeta(l,m,k)
          drz(l,m,k) = ( a15(1)*rf(l,m,1)+a15(2)*rf(l,m,2) &
                       + a15(3)*rf(l,m,3)+a15(4)*rf(l,m,4) &
                       + a15(5)*rf(l,m,5)+a15(6)*rf(l,m,6) &
                       + a15(7)*rf(l,m,7) ) *idz(k)*costeta(l,m,k)
       enddo
    enddo

    k=3
    do l=1,ngh
       do m=1,ngh
          duz(l,m,k) = ( a24(1)*uf(l,m,1)+a24(2)*uf(l,m,2) &
                       + a24(3)*uf(l,m,3)+a24(4)*uf(l,m,4) &
                       + a24(5)*uf(l,m,5)+a24(6)*uf(l,m,6) &
                       + a24(7)*uf(l,m,7) ) *idz(k)*costeta(l,m,k)
          dvz(l,m,k) = ( a24(1)*vf(l,m,1)+a24(2)*vf(l,m,2) &
                       + a24(3)*vf(l,m,3)+a24(4)*vf(l,m,4) &
                       + a24(5)*vf(l,m,5)+a24(6)*vf(l,m,6) &
                       + a24(7)*vf(l,m,7) ) *idz(k)*costeta(l,m,k)
          dwz(l,m,k) = ( a24(1)*wf(l,m,1)+a24(2)*wf(l,m,2) &
                       + a24(3)*wf(l,m,3)+a24(4)*wf(l,m,4) &
                       + a24(5)*wf(l,m,5)+a24(6)*wf(l,m,6) &
                       + a24(7)*wf(l,m,7) ) *idz(k)*costeta(l,m,k)
          dpz(l,m,k) = ( a24(1)*pf(l,m,1)+a24(2)*pf(l,m,2) &
                       + a24(3)*pf(l,m,3)+a24(4)*pf(l,m,4) &
                       + a24(5)*pf(l,m,5)+a24(6)*pf(l,m,6) &
                       + a24(7)*pf(l,m,7) ) *idz(k)*costeta(l,m,k)
          drz(l,m,k) = ( a24(1)*rf(l,m,1)+a24(2)*rf(l,m,2) &
                       + a24(3)*rf(l,m,3)+a24(4)*rf(l,m,4) &
                       + a24(5)*rf(l,m,5)+a24(6)*rf(l,m,6) &
                       + a24(7)*rf(l,m,7) ) *idz(k)*costeta(l,m,k)
       enddo
    enddo

    do k=4,ngh
       do l=1,ngh
          do m=1,ngh
             duz(l,m,k) = (a7(1)*(uf(l,m,k+1) - uf(l,m,k-1)) + &
                           a7(2)*(uf(l,m,k+2) - uf(l,m,k-2)) + &
                           a7(3)*(uf(l,m,k+3) - uf(l,m,k-3)) ) *idz(k)*costeta(l,m,k)
             dvz(l,m,k) = (a7(1)*(vf(l,m,k+1) - vf(l,m,k-1)) + &
                           a7(2)*(vf(l,m,k+2) - vf(l,m,k-2)) + &
                           a7(3)*(vf(l,m,k+3) - vf(l,m,k-3)) ) *idz(k)*costeta(l,m,k)
             dwz(l,m,k) = (a7(1)*(wf(l,m,k+1) - wf(l,m,k-1)) + &
                           a7(2)*(wf(l,m,k+2) - wf(l,m,k-2)) + &
                           a7(3)*(wf(l,m,k+3) - wf(l,m,k-3)) ) *idz(k)*costeta(l,m,k)
             dpz(l,m,k) = (a7(1)*(pf(l,m,k+1) - pf(l,m,k-1)) + &
                           a7(2)*(pf(l,m,k+2) - pf(l,m,k-2)) + &
                           a7(3)*(pf(l,m,k+3) - pf(l,m,k-3)) ) *idz(k)*costeta(l,m,k)
             drz(l,m,k) = (a7(1)*(rf(l,m,k+1) - rf(l,m,k-1)) + &
                           a7(2)*(rf(l,m,k+2) - rf(l,m,k-2)) + &
                           a7(3)*(rf(l,m,k+3) - rf(l,m,k-3)) ) *idz(k)*costeta(l,m,k)
          enddo
       enddo
    enddo

    ! Compute rt,ut,vt,wt & pt
    ! ========================
    do l=1,ngh
       do m=1,ngh
          do k=1,ngh
             pt(l,m,k) = vg(l,m,k)*(dpx(l,m,k)+dpy(l,m,k)+dpz(l,m,k)+pf(l,m,k)*ir(l,m,k))
             ut(l,m,k) = vg(l,m,k)*(dux(l,m,k)+duy(l,m,k)+duz(l,m,k)+uf(l,m,k)*ir(l,m,k))
             vt(l,m,k) = vg(l,m,k)*(dvx(l,m,k)+dvy(l,m,k)+dvz(l,m,k)+vf(l,m,k)*ir(l,m,k))
             wt(l,m,k) = vg(l,m,k)*(dwx(l,m,k)+dwy(l,m,k)+dwz(l,m,k)+wf(l,m,k)*ir(l,m,k))
             rt(l,m,k) = vg(l,m,k)*(drx(l,m,k)+dry(l,m,k)+drz(l,m,k)+rf(l,m,k)*ir(l,m,k))
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

  end subroutine bc_TD3d_imax_jmax_kmin

  !===============================================================================
  module subroutine bc_TD3d_imax_jmax_kmax
  !===============================================================================
    !> 3D Tam & Dong's BC at imax-jmax-kmax (corner 2,2,2 /right-top-back)
    !> - Cartesian version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l,m,n
    real(wp) :: cp,av,c2_
    real(wp), dimension(-2:ngh,-2:ngh,-2:ngh) :: rf,uf,vf,wf,pf
    real(wp), dimension(ngh,ngh,ngh) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(ngh,ngh,ngh) :: dux,dvx,dwx,dpx,drx
    real(wp), dimension(ngh,ngh,ngh) :: duy,dvy,dwy,dpy,dry
    real(wp), dimension(ngh,ngh,ngh) :: duz,dvz,dwz,dpz,drz
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
	           	BC_face(2,2)%U0(i,m,k,4)*costeta(l,m,n)	       &
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
             dux(l,m,n) = ( a7(1)* (uf(l+1,m,n)-uf(l-1,m,n)) &
                          + a7(2)* (uf(l+2,m,n)-uf(l-2,m,n)) &
                          + a7(3)* (uf(l+3,m,n)-uf(l-3,m,n)) ) *idx(i)*sintcosp(l,m,n)
             dvx(l,m,n) = ( a7(1)* (vf(l+1,m,n)-vf(l-1,m,n)) &
                          + a7(2)* (vf(l+2,m,n)-vf(l-2,m,n)) &
                          + a7(3)* (vf(l+3,m,n)-vf(l-3,m,n)) ) *idx(i)*sintcosp(l,m,n)
             dwx(l,m,n) = ( a7(1)* (wf(l+1,m,n)-wf(l-1,m,n)) &
                          + a7(2)* (wf(l+2,m,n)-wf(l-2,m,n)) &
                          + a7(3)* (wf(l+3,m,n)-wf(l-3,m,n)) ) *idx(i)*sintcosp(l,m,n)
             dpx(l,m,n) = ( a7(1)* (pf(l+1,m,n)-pf(l-1,m,n)) &
                          + a7(2)* (pf(l+2,m,n)-pf(l-2,m,n)) &
                          + a7(3)* (pf(l+3,m,n)-pf(l-3,m,n)) ) *idx(i)*sintcosp(l,m,n)
             drx(l,m,n) = ( a7(1)* (rf(l+1,m,n)-rf(l-1,m,n)) &
                          + a7(2)* (rf(l+2,m,n)-rf(l-2,m,n)) &
                          + a7(3)* (rf(l+3,m,n)-rf(l-3,m,n)) ) *idx(i)*sintcosp(l,m,n)
          enddo
       enddo
    enddo
    
    i=nx-2
    l=i-nxmngh
    do m=1,ngh
       do n=1,ngh
          dux(l,m,n) = ( a42(1)*uf(l+2,m,n)+a42(2)*uf(l+1,m,n) &
                       + a42(3)*uf(l  ,m,n)+a42(4)*uf(l-1,m,n) &
                       + a42(5)*uf(l-2,m,n)+a42(6)*uf(l-3,m,n) &
                       + a42(7)*uf(l-4,m,n) ) *idx(i)*sintcosp(l,m,n)
          dvx(l,m,n) = ( a42(1)*vf(l+2,m,n)+a42(2)*vf(l+1,m,n) &
                       + a42(3)*vf(l  ,m,n)+a42(4)*vf(l-1,m,n) &
                       + a42(5)*vf(l-2,m,n)+a42(6)*vf(l-3,m,n) &
                       + a42(7)*vf(l-4,m,n) ) *idx(i)*sintcosp(l,m,n)
          dwx(l,m,n) = ( a42(1)*wf(l+2,m,n)+a42(2)*wf(l+1,m,n) &
                       + a42(3)*wf(l  ,m,n)+a42(4)*wf(l-1,m,n) &
                       + a42(5)*wf(l-2,m,n)+a42(6)*wf(l-3,m,n) &
                       + a42(7)*wf(l-4,m,n) ) *idx(i)*sintcosp(l,m,n)
          dpx(l,m,n) = ( a42(1)*pf(l+2,m,n)+a42(2)*pf(l+1,m,n) &
                       + a42(3)*pf(l  ,m,n)+a42(4)*pf(l-1,m,n) &
                       + a42(5)*pf(l-2,m,n)+a42(6)*pf(l-3,m,n) &
                       + a42(7)*pf(l-4,m,n) ) *idx(i)*sintcosp(l,m,n)
          drx(l,m,n) = ( a42(1)*rf(l+2,m,n)+a42(2)*rf(l+1,m,n) &
                       + a42(3)*rf(l  ,m,n)+a42(4)*rf(l-1,m,n) &
                       + a42(5)*rf(l-2,m,n)+a42(6)*rf(l-3,m,n) &
                       + a42(7)*rf(l-4,m,n) ) *idx(i)*sintcosp(l,m,n)
       enddo
    enddo
    
    i=nx-1
    l=i-nxmngh
    do m=1,ngh
       do n=1,ngh
          dux(l,m,n) = ( a51(1)*uf(l+1,m,n)+a51(2)*uf(l  ,m,n) &
                       + a51(3)*uf(l-1,m,n)+a51(4)*uf(l-2,m,n) &
                       + a51(5)*uf(l-3,m,n)+a51(6)*uf(l-4,m,n) &
                       + a51(7)*uf(l-5,m,n) ) *idx(i)*sintcosp(l,m,n)
          dvx(l,m,n) = ( a51(1)*vf(l+1,m,n)+a51(2)*vf(l  ,m,n) &
                       + a51(3)*vf(l-1,m,n)+a51(4)*vf(l-2,m,n) &
                       + a51(5)*vf(l-3,m,n)+a51(6)*vf(l-4,m,n) &
                       + a51(7)*vf(l-5,m,n) ) *idx(i)*sintcosp(l,m,n)
          dwx(l,m,n) = ( a51(1)*wf(l+1,m,n)+a51(2)*wf(l  ,m,n) &
                       + a51(3)*wf(l-1,m,n)+a51(4)*wf(l-2,m,n) &
                       + a51(5)*wf(l-3,m,n)+a51(6)*wf(l-4,m,n) &
                       + a51(7)*wf(l-5,m,n) ) *idx(i)*sintcosp(l,m,n)
          dpx(l,m,n) = ( a51(1)*pf(l+1,m,n)+a51(2)*pf(l  ,m,n) &
                       + a51(3)*pf(l-1,m,n)+a51(4)*pf(l-2,m,n) &
                       + a51(5)*pf(l-3,m,n)+a51(6)*pf(l-4,m,n) &
                       + a51(7)*pf(l-5,m,n) ) *idx(i)*sintcosp(l,m,n)
          drx(l,m,n) = ( a51(1)*rf(l+1,m,n)+a51(2)*rf(l  ,m,n) &
                       + a51(3)*rf(l-1,m,n)+a51(4)*rf(l-2,m,n) &
                       + a51(5)*rf(l-3,m,n)+a51(6)*rf(l-4,m,n) &
                       + a51(7)*rf(l-5,m,n) ) *idx(i)*sintcosp(l,m,n)
       enddo
    enddo
    
    i=nx
    l=i-nxmngh
    do m=1,ngh
       do n=1,ngh
          dux(l,m,n) = ( a60(1)*uf(l  ,m,n)+a60(2)*uf(l-1,m,n) &
                       + a60(3)*uf(l-2,m,n)+a60(4)*uf(l-3,m,n) &
                       + a60(5)*uf(l-4,m,n)+a60(6)*uf(l-5,m,n) &
                       + a60(7)*uf(l-6,m,n) ) *idx(i)*sintcosp(l,m,n)
          dvx(l,m,n) = ( a60(1)*vf(l  ,m,n)+a60(2)*vf(l-1,m,n) &
                       + a60(3)*vf(l-2,m,n)+a60(4)*vf(l-3,m,n) &
                       + a60(5)*vf(l-4,m,n)+a60(6)*vf(l-5,m,n) &
                       + a60(7)*vf(l-6,m,n) ) *idx(i)*sintcosp(l,m,n)
          dwx(l,m,n) = ( a60(1)*wf(l  ,m,n)+a60(2)*wf(l-1,m,n) &
                       + a60(3)*wf(l-2,m,n)+a60(4)*wf(l-3,m,n) &
                       + a60(5)*wf(l-4,m,n)+a60(6)*wf(l-5,m,n) &
                       + a60(7)*wf(l-6,m,n) ) *idx(i)*sintcosp(l,m,n)
          dpx(l,m,n) = ( a60(1)*pf(l  ,m,n)+a60(2)*pf(l-1,m,n) &
                       + a60(3)*pf(l-2,m,n)+a60(4)*pf(l-3,m,n) &
                       + a60(5)*pf(l-4,m,n)+a60(6)*pf(l-5,m,n) &
                       + a60(7)*pf(l-6,m,n) ) *idx(i)*sintcosp(l,m,n)
          drx(l,m,n) = ( a60(1)*rf(l  ,m,n)+a60(2)*rf(l-1,m,n) &
                       + a60(3)*rf(l-2,m,n)+a60(4)*rf(l-3,m,n) &
                       + a60(5)*rf(l-4,m,n)+a60(6)*rf(l-5,m,n) &
                       + a60(7)*rf(l-6,m,n) ) *idx(i)*sintcosp(l,m,n)
       enddo
    enddo
 
    ! Non-centered derivatives in y-direction *sin(teta)*sin(phi)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    do j=nymnghp1,ny-3
       m=j-nymngh
       do l=1,ngh
          do n=1,ngh
             duy(l,m,n) = ( a7(1)*(uf(l,m+1,n)-uf(l,m-1,n)) &
                          + a7(2)*(uf(l,m+2,n)-uf(l,m-2,n)) &
                          + a7(3)*(uf(l,m+3,n)-uf(l,m-3,n)) ) *idy(j)*sintsinp(l,m,n)
             dvy(l,m,n) = ( a7(1)*(vf(l,m+1,n)-vf(l,m-1,n)) &
                          + a7(2)*(vf(l,m+2,n)-vf(l,m-2,n)) &
                          + a7(3)*(vf(l,m+3,n)-vf(l,m-3,n)) ) *idy(j)*sintsinp(l,m,n)
             dwy(l,m,n) = ( a7(1)*(wf(l,m+1,n)-wf(l,m-1,n)) &
                          + a7(2)*(wf(l,m+2,n)-wf(l,m-2,n)) &
                          + a7(3)*(wf(l,m+3,n)-wf(l,m-3,n)) ) *idy(j)*sintsinp(l,m,n)
             dpy(l,m,n) = ( a7(1)*(pf(l,m+1,n)-pf(l,m-1,n)) &
                          + a7(2)*(pf(l,m+2,n)-pf(l,m-2,n)) &
                          + a7(3)*(pf(l,m+3,n)-pf(l,m-3,n)) ) *idy(j)*sintsinp(l,m,n)
             dry(l,m,n) = ( a7(1)*(rf(l,m+1,n)-rf(l,m-1,n)) &
                          + a7(2)*(rf(l,m+2,n)-rf(l,m-2,n)) &
                          + a7(3)*(rf(l,m+3,n)-rf(l,m-3,n)) ) *idy(j)*sintsinp(l,m,n)
          enddo
       enddo
    enddo

    j=ny-2
    m=j-nymngh
    do l=1,ngh
       do n=1,ngh
          duy(l,m,n) = ( a42(1)*uf(l,m+2,n)+a42(2)*uf(l,m+1,n) &
                       + a42(3)*uf(l,m  ,n)+a42(4)*uf(l,m-1,n) &
                       + a42(5)*uf(l,m-2,n)+a42(6)*uf(l,m-3,n) &
                       + a42(7)*uf(l,m-4,n) ) *idy(j)*sintsinp(l,m,n)
          dvy(l,m,n) = ( a42(1)*vf(l,m+2,n)+a42(2)*vf(l,m+1,n) &
                       + a42(3)*vf(l,m  ,n)+a42(4)*vf(l,m-1,n) &
                       + a42(5)*vf(l,m-2,n)+a42(6)*vf(l,m-3,n) &
                       + a42(7)*vf(l,m-4,n) ) *idy(j)*sintsinp(l,m,n)
          dwy(l,m,n) = ( a42(1)*wf(l,m+2,n)+a42(2)*wf(l,m+1,n) &
                       + a42(3)*wf(l,m  ,n)+a42(4)*wf(l,m-1,n) &
                       + a42(5)*wf(l,m-2,n)+a42(6)*wf(l,m-3,n) &
                       + a42(7)*wf(l,m-4,n) ) *idy(j)*sintsinp(l,m,n)
          dpy(l,m,n) = ( a42(1)*pf(l,m+2,n)+a42(2)*pf(l,m+1,n) &
                       + a42(3)*pf(l,m  ,n)+a42(4)*pf(l,m-1,n) &
                       + a42(5)*pf(l,m-2,n)+a42(6)*pf(l,m-3,n) &
                       + a42(7)*pf(l,m-4,n) ) *idy(j)*sintsinp(l,m,n)
          dry(l,m,n) = ( a42(1)*rf(l,m+2,n)+a42(2)*rf(l,m+1,n) &
                       + a42(3)*rf(l,m  ,n)+a42(4)*rf(l,m-1,n) &
                       + a42(5)*rf(l,m-2,n)+a42(6)*rf(l,m-3,n) &
                       + a42(7)*rf(l,m-4,n) ) *idy(j)*sintsinp(l,m,n)
       enddo
    enddo
    
    j=ny-1
    m=j-nymngh
    do l=1,ngh
       do n=1,ngh
          duy(l,m,n) = ( a51(1)*uf(l,m+1,n)+a51(2)*uf(l,m  ,n) &
                       + a51(3)*uf(l,m-1,n)+a51(4)*uf(l,m-2,n) &
                       + a51(5)*uf(l,m-3,n)+a51(6)*uf(l,m-4,n) &
                       + a51(7)*uf(l,m-5,n) ) *idy(j)*sintsinp(l,m,n)
          dvy(l,m,n) = ( a51(1)*vf(l,m+1,n)+a51(2)*vf(l,m  ,n) &
                       + a51(3)*vf(l,m-1,n)+a51(4)*vf(l,m-2,n) &
                       + a51(5)*vf(l,m-3,n)+a51(6)*vf(l,m-4,n) &
                       + a51(7)*vf(l,m-5,n) ) *idy(j)*sintsinp(l,m,n)
          dwy(l,m,n) = ( a51(1)*wf(l,m+1,n)+a51(2)*wf(l,m  ,n) &
                       + a51(3)*wf(l,m-1,n)+a51(4)*wf(l,m-2,n) &
                       + a51(5)*wf(l,m-3,n)+a51(6)*wf(l,m-4,n) &
                       + a51(7)*wf(l,m-5,n) ) *idy(j)*sintsinp(l,m,n)
          dpy(l,m,n) = ( a51(1)*pf(l,m+1,n)+a51(2)*pf(l,m  ,n) &
                       + a51(3)*pf(l,m-1,n)+a51(4)*pf(l,m-2,n) &
                       + a51(5)*pf(l,m-3,n)+a51(6)*pf(l,m-4,n) &
                       + a51(7)*pf(l,m-5,n) ) *idy(j)*sintsinp(l,m,n)
          dry(l,m,n) = ( a51(1)*rf(l,m+1,n)+a51(2)*rf(l,m  ,n) &
                       + a51(3)*rf(l,m-1,n)+a51(4)*rf(l,m-2,n) &
                       + a51(5)*rf(l,m-3,n)+a51(6)*rf(l,m-4,n) &
                       + a51(7)*rf(l,m-5,n) ) *idy(j)*sintsinp(l,m,n)
       enddo
    enddo
    
    j=ny
    m=j-nymngh
    do l=1,ngh
       do n=1,ngh
          duy(l,m,n) = ( a60(1)*uf(l,m  ,n)+a60(2)*uf(l,m-1,n) &
                       + a60(3)*uf(l,m-2,n)+a60(4)*uf(l,m-3,n) &
                       + a60(5)*uf(l,m-4,n)+a60(6)*uf(l,m-5,n) &
                       + a60(7)*uf(l,m-6,n) ) *idy(j)*sintsinp(l,m,n)
          dvy(l,m,n) = ( a60(1)*vf(l,m  ,n)+a60(2)*vf(l,m-1,n) &
                       + a60(3)*vf(l,m-2,n)+a60(4)*vf(l,m-3,n) &
                       + a60(5)*vf(l,m-4,n)+a60(6)*vf(l,m-5,n) &
                       + a60(7)*vf(l,m-6,n) ) *idy(j)*sintsinp(l,m,n)
          dwy(l,m,n) = ( a60(1)*wf(l,m  ,n)+a60(2)*wf(l,m-1,n) &
                       + a60(3)*wf(l,m-2,n)+a60(4)*wf(l,m-3,n) &
                       + a60(5)*wf(l,m-4,n)+a60(6)*wf(l,m-5,n) &
                       + a60(7)*wf(l,m-6,n) ) *idy(j)*sintsinp(l,m,n)
          dpy(l,m,n) = ( a60(1)*pf(l,m  ,n)+a60(2)*pf(l,m-1,n) &
                       + a60(3)*pf(l,m-2,n)+a60(4)*pf(l,m-3,n) &
                       + a60(5)*pf(l,m-4,n)+a60(6)*pf(l,m-5,n) &
                       + a60(7)*pf(l,m-6,n) ) *idy(j)*sintsinp(l,m,n)
          dry(l,m,n) = ( a60(1)*rf(l,m  ,n)+a60(2)*rf(l,m-1,n) &
                       + a60(3)*rf(l,m-2,n)+a60(4)*rf(l,m-3,n) &
                       + a60(5)*rf(l,m-4,n)+a60(6)*rf(l,m-5,n) &
                       + a60(7)*rf(l,m-6,n) ) *idy(j)*sintsinp(l,m,n)
       enddo
    enddo
              
    ! Non-centered derivatives in z-direction *cos(teta)
    ! =======================================
    ! (Tam & Webb DRP schemes)
    do k=nzmnghp1,nz-3
       n=k-nzmngh
       do l=1,ngh
          do m=1,ngh
             duz(l,m,n) = (a7(1)*(uf(l,m,n+1) - uf(l,m,n-1)) + &
                           a7(2)*(uf(l,m,n+2) - uf(l,m,n-2)) + &
                           a7(3)*(uf(l,m,n+3) - uf(l,m,n-3)) ) *idz(k)*costeta(l,m,n)
             dvz(l,m,n) = (a7(1)*(vf(l,m,n+1) - vf(l,m,n-1)) + &
                           a7(2)*(vf(l,m,n+2) - vf(l,m,n-2)) + &
                           a7(3)*(vf(l,m,n+3) - vf(l,m,n-3)) ) *idz(k)*costeta(l,m,n)
             dwz(l,m,n) = (a7(1)*(wf(l,m,n+1) - wf(l,m,n-1)) + &
                           a7(2)*(wf(l,m,n+2) - wf(l,m,n-2)) + &
                           a7(3)*(wf(l,m,n+3) - wf(l,m,n-3)) ) *idz(k)*costeta(l,m,n)
             dpz(l,m,n) = (a7(1)*(pf(l,m,n+1) - pf(l,m,n-1)) + &
                           a7(2)*(pf(l,m,n+2) - pf(l,m,n-2)) + &
                           a7(3)*(pf(l,m,n+3) - pf(l,m,n-3)) ) *idz(k)*costeta(l,m,n)
             drz(l,m,n) = (a7(1)*(rf(l,m,n+1) - rf(l,m,n-1)) + &
                           a7(2)*(rf(l,m,n+2) - rf(l,m,n-2)) + &
                           a7(3)*(rf(l,m,n+3) - rf(l,m,n-3)) ) *idz(k)*costeta(l,m,n)
          enddo
       enddo
    enddo

    k=nz-2
    n=k-nzmngh
    do l=1,ngh
       do m=1,ngh
          duz(l,m,n) = ( a42(1)*uf(l,m,n+2)+a42(2)*uf(l,m,n+1) &
                       + a42(3)*uf(l,m,n  )+a42(4)*uf(l,m,n-1) &
                       + a42(5)*uf(l,m,n-2)+a42(6)*uf(l,m,n-3) &
                       + a42(7)*uf(l,m,n-4) ) *idz(k)*costeta(l,m,n)
          dvz(l,m,n) = ( a42(1)*vf(l,m,n+2)+a42(2)*vf(l,m,n+1) &
                       + a42(3)*vf(l,m,n  )+a42(4)*vf(l,m,n-1) &
                       + a42(5)*vf(l,m,n-2)+a42(6)*vf(l,m,n-3) &
                       + a42(7)*vf(l,m,n-4) ) *idz(k)*costeta(l,m,n)
          dwz(l,m,n) = ( a42(1)*wf(l,m,n+2)+a42(2)*wf(l,m,n+1) &
                       + a42(3)*wf(l,m,n  )+a42(4)*wf(l,m,n-1) &
                       + a42(5)*wf(l,m,n-2)+a42(6)*wf(l,m,n-3) &
                       + a42(7)*wf(l,m,n-4) ) *idz(k)*costeta(l,m,n)
          dpz(l,m,n) = ( a42(1)*pf(l,m,n+2)+a42(2)*pf(l,m,n+1) &
                       + a42(3)*pf(l,m,n  )+a42(4)*pf(l,m,n-1) &
                       + a42(5)*pf(l,m,n-2)+a42(6)*pf(l,m,n-3) &
                       + a42(7)*pf(l,m,n-4) ) *idz(k)*costeta(l,m,n)
          drz(l,m,n) = ( a42(1)*rf(l,m,n+2)+a42(2)*rf(l,m,n+1) &
                       + a42(3)*rf(l,m,n  )+a42(4)*rf(l,m,n-1) &
                       + a42(5)*rf(l,m,n-2)+a42(6)*rf(l,m,n-3) &
                       + a42(7)*rf(l,m,n-4) ) *idz(k)*costeta(l,m,n)
       enddo
    enddo

    k=nz-1
    n=k-nzmngh
    do l=1,ngh
       do m=1,ngh
          duz(l,m,n) = ( a51(1)*uf(l,m,n+1)+a51(2)*uf(l,m,n  ) &
                       + a51(3)*uf(l,m,n-1)+a51(4)*uf(l,m,n-2) &
                       + a51(5)*uf(l,m,n-3)+a51(6)*uf(l,m,n-4) &
                       + a51(7)*uf(l,m,n-5) ) *idz(k)*costeta(l,m,n)
          dvz(l,m,n) = ( a51(1)*vf(l,m,n+1)+a51(2)*vf(l,m,n  ) &
                       + a51(3)*vf(l,m,n-1)+a51(4)*vf(l,m,n-2) &
                       + a51(5)*vf(l,m,n-3)+a51(6)*vf(l,m,n-4) &
                       + a51(7)*vf(l,m,n-5) ) *idz(k)*costeta(l,m,n)
          dwz(l,m,n) = ( a51(1)*wf(l,m,n+1)+a51(2)*wf(l,m,n  ) &
                       + a51(3)*wf(l,m,n-1)+a51(4)*wf(l,m,n-2) &
                       + a51(5)*wf(l,m,n-3)+a51(6)*wf(l,m,n-4) &
                       + a51(7)*wf(l,m,n-5) ) *idz(k)*costeta(l,m,n)
          dpz(l,m,n) = ( a51(1)*pf(l,m,n+1)+a51(2)*pf(l,m,n  ) &
                       + a51(3)*pf(l,m,n-1)+a51(4)*pf(l,m,n-2) &
                       + a51(5)*pf(l,m,n-3)+a51(6)*pf(l,m,n-4) &
                       + a51(7)*pf(l,m,n-5) ) *idz(k)*costeta(l,m,n)
          drz(l,m,n) = ( a51(1)*rf(l,m,n+1)+a51(2)*rf(l,m,n  ) &
                       + a51(3)*rf(l,m,n-1)+a51(4)*rf(l,m,n-2) &
                       + a51(5)*rf(l,m,n-3)+a51(6)*rf(l,m,n-4) &
                       + a51(7)*rf(l,m,n-5) ) *idz(k)*costeta(l,m,n)
       enddo
    enddo

    k=nz
    n=k-nzmngh
    do l=1,ngh
       do m=1,ngh
          duz(l,m,n) = ( a60(1)*uf(l,m,n  )+a60(2)*uf(l,m,n-1) &
                       + a60(3)*uf(l,m,n-2)+a60(4)*uf(l,m,n-3) &
                       + a60(5)*uf(l,m,n-4)+a60(6)*uf(l,m,n-5) &
                       + a60(7)*uf(l,m,n-6) ) *idz(k)*costeta(l,m,n)
          dvz(l,m,n) = ( a60(1)*vf(l,m,n  )+a60(2)*vf(l,m,n-1) &
                       + a60(3)*vf(l,m,n-2)+a60(4)*vf(l,m,n-3) &
                       + a60(5)*vf(l,m,n-4)+a60(6)*vf(l,m,n-5) &
                       + a60(7)*vf(l,m,n-6) ) *idz(k)*costeta(l,m,n)
          dwz(l,m,n) = ( a60(1)*wf(l,m,n  )+a60(2)*wf(l,m,n-1) &
                       + a60(3)*wf(l,m,n-2)+a60(4)*wf(l,m,n-3) &
                       + a60(5)*wf(l,m,n-4)+a60(6)*wf(l,m,n-5) &
                       + a60(7)*wf(l,m,n-6) ) *idz(k)*costeta(l,m,n)
          dpz(l,m,n) = ( a60(1)*pf(l,m,n  )+a60(2)*pf(l,m,n-1) &
                       + a60(3)*pf(l,m,n-2)+a60(4)*pf(l,m,n-3) &
                       + a60(5)*pf(l,m,n-4)+a60(6)*pf(l,m,n-5) &
                       + a60(7)*pf(l,m,n-6) ) *idz(k)*costeta(l,m,n)
          drz(l,m,n) = ( a60(1)*rf(l,m,n  )+a60(2)*rf(l,m,n-1) &
                       + a60(3)*rf(l,m,n-2)+a60(4)*rf(l,m,n-3) &
                       + a60(5)*rf(l,m,n-4)+a60(6)*rf(l,m,n-5) &
                       + a60(7)*rf(l,m,n-6) ) *idz(k)*costeta(l,m,n)
       enddo
    enddo

    ! Compute rt,ut,vt,wt & pt
    ! ========================
    do l=1,ngh
       do m=1,ngh
          do n=1,ngh
             pt(l,m,n) = vg(l,m,n)*(dpx(l,m,n)+dpy(l,m,n)+dpz(l,m,n)+pf(l,m,n)*ir(l,m,n))
             ut(l,m,n) = vg(l,m,n)*(dux(l,m,n)+duy(l,m,n)+duz(l,m,n)+uf(l,m,n)*ir(l,m,n))
             vt(l,m,n) = vg(l,m,n)*(dvx(l,m,n)+dvy(l,m,n)+dvz(l,m,n)+vf(l,m,n)*ir(l,m,n))
             wt(l,m,n) = vg(l,m,n)*(dwx(l,m,n)+dwy(l,m,n)+dwz(l,m,n)+wf(l,m,n)*ir(l,m,n))
             rt(l,m,n) = vg(l,m,n)*(drx(l,m,n)+dry(l,m,n)+drz(l,m,n)+rf(l,m,n)*ir(l,m,n))             
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

  end subroutine bc_TD3d_imax_jmax_kmax

end submodule smod_TamDong3d_corner
