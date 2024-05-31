!===============================================================================
submodule (mod_TamDong3d) smod_TamDong3d_edgez
!===============================================================================
  !> author: XG
  !> date: February 2020 - modif January 2022
  !> Non-reflecting boundary conditions of Tam and Dong
  !> - 3D (periodic) Cartesian version - edges along z
!=============================================================================== 

contains

  !===============================================================================
  module subroutine bc_TD3d_imin_jmin
  !===============================================================================
    !> 3D Tam & Dong's BC: boundary condition at imin-jmin (edge 1,1,1 /left-bottom)
    !> - Cartesian version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: cp,av,c2_
    real(wp), dimension(1:ngh+3,1:ngh+3,nz1:nz2) :: rf,uf,vf,wf,pf
    real(wp), dimension(1:ngh,1:ngh,nz) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(1:ngh,1:ngh,nz) :: dux,dvx,dwx,dpx,drx
    real(wp), dimension(1:ngh,1:ngh,nz) :: duy,dvy,dwy,dpy,dry
    real(wp) :: duz,dvz,dwz,dpz,drz
    !-------------------------------------------------------------------------

    ! Pointers for spherical coordinates
    ! ==================================
    ir=>BC_edge(1,1,1)%i_r
    cosphi=>BC_edge(1,1,1)%cosp
    sinphi=>BC_edge(1,1,1)%sinp
    costeta=>BC_edge(1,1,1)%cost
    sinteta=>BC_edge(1,1,1)%sint
    costcosp=>BC_edge(1,1,1)%costcosp
    costsinp=>BC_edge(1,1,1)%costsinp
    sintcosp=>BC_edge(1,1,1)%sintcosp
    sintsinp=>BC_edge(1,1,1)%sintsinp
    
    ! Compute fluctuations
    ! ====================
    do i=1,nghp3
       do j=1,nghp3
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
       do j=1,ngh
          do k=ndz,nfz
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
       do k=ndz,nfz
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
       do k=ndz,nfz
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
       do k=ndz,nfz
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
          do k=ndz,nfz
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
       do k=ndz,nfz
          duy(i,j,k)=(uf(i,j+1,k)-uf(i,j,k))*idy1_jmin*sintsinp(i,j,k)
          dvy(i,j,k)=(vf(i,j+1,k)-vf(i,j,k))*idy1_jmin*sintsinp(i,j,k)
          dwy(i,j,k)=(wf(i,j+1,k)-wf(i,j,k))*idy1_jmin*sintsinp(i,j,k)
          dpy(i,j,k)=(pf(i,j+1,k)-pf(i,j,k))*idy1_jmin*sintsinp(i,j,k)
          dry(i,j,k)=(rf(i,j+1,k)-rf(i,j,k))*idy1_jmin*sintsinp(i,j,k)
       enddo
    enddo

    j=2
    do i=1,ngh
       do k=ndz,nfz
          duy(i,j,k)=a3(1)*(uf(i,j+1,k)-uf(i,j-1,k))*idy2_jmin*sintsinp(i,j,k)
          dvy(i,j,k)=a3(1)*(vf(i,j+1,k)-vf(i,j-1,k))*idy2_jmin*sintsinp(i,j,k)
          dwy(i,j,k)=a3(1)*(wf(i,j+1,k)-wf(i,j-1,k))*idy2_jmin*sintsinp(i,j,k)
          dpy(i,j,k)=a3(1)*(pf(i,j+1,k)-pf(i,j-1,k))*idy2_jmin*sintsinp(i,j,k)
          dry(i,j,k)=a3(1)*(rf(i,j+1,k)-rf(i,j-1,k))*idy2_jmin*sintsinp(i,j,k)
       enddo
    enddo

    j=3
    do i=1,ngh
       do k=ndz,nfz
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
!!$       do k=ndz,nfz
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
!!$       do k=ndz,nfz
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
!!$       do k=ndz,nfz
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
          do k=ndz,nfz
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

    ! Compute rt,ut,vt,wt & pt
    ! ========================
    do i=1,ngh
       do j=1,ngh
          do k=ndz,nfz

             duz = (a7(1)*( uf(i,j,k+1) - uf(i,j,k-1) ) + &
                    a7(2)*( uf(i,j,k+2) - uf(i,j,k-2) ) + &
                    a7(3)*( uf(i,j,k+3) - uf(i,j,k-3) ) ) *idz(k)*costeta(i,j,k)
             dvz = (a7(1)*( vf(i,j,k+1) - vf(i,j,k-1) ) + &
                    a7(2)*( vf(i,j,k+2) - vf(i,j,k-2) ) + &
                    a7(3)*( vf(i,j,k+3) - vf(i,j,k-3) ) ) *idz(k)*costeta(i,j,k)
             dwz = (a7(1)*( wf(i,j,k+1) - wf(i,j,k-1) ) + &
                    a7(2)*( wf(i,j,k+2) - wf(i,j,k-2) ) + &
                    a7(3)*( wf(i,j,k+3) - wf(i,j,k-3) ) ) *idz(k)*costeta(i,j,k)
             dpz = (a7(1)*( pf(i,j,k+1) - pf(i,j,k-1) ) + &
                    a7(2)*( pf(i,j,k+2) - pf(i,j,k-2) ) + &
                    a7(3)*( pf(i,j,k+3) - pf(i,j,k-3) ) ) *idz(k)*costeta(i,j,k)
             drz = (a7(1)*( rf(i,j,k+1) - rf(i,j,k-1) ) + &
                    a7(2)*( rf(i,j,k+2) - rf(i,j,k-2) ) + &
                    a7(3)*( rf(i,j,k+3) - rf(i,j,k-3) ) ) *idz(k)*costeta(i,j,k)

             pt(i,j,k) = vg(i,j,k)*(dpx(i,j,k)+dpy(i,j,k)+dpz+pf(i,j,k)*ir(i,j,k))
             ut(i,j,k) = vg(i,j,k)*(dux(i,j,k)+duy(i,j,k)+duz+uf(i,j,k)*ir(i,j,k))
             vt(i,j,k) = vg(i,j,k)*(dvx(i,j,k)+dvy(i,j,k)+dvz+vf(i,j,k)*ir(i,j,k))
             wt(i,j,k) = vg(i,j,k)*(dwx(i,j,k)+dwy(i,j,k)+dwz+wf(i,j,k)*ir(i,j,k))
             rt(i,j,k) = vg(i,j,k)*(drx(i,j,k)+dry(i,j,k)+drz+rf(i,j,k)*ir(i,j,k))
          enddo
       enddo
    enddo

    ! Update fluxes at each RK step
    ! =============================
    do k=ndz,nfz
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

  end subroutine bc_TD3d_imin_jmin

  !===============================================================================
  module subroutine bc_TD3d_imin_jmax
  !===============================================================================
    !> 3D Tam & Dong's BC: boundary condition at imin-jmax (edge 1,1,2 /left-top)
    !> - Cartesian version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l
    real(wp) :: cp,av,c2_
    real(wp), dimension(1:ngh+3,-2:ngh,nz1:nz2) :: rf,uf,vf,wf,pf
    real(wp), dimension(1:ngh,1:ngh,nz) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(1:ngh,1:ngh,nz) :: dux,dvx,dwx,dpx,drx
    real(wp), dimension(1:ngh,1:ngh,nz) :: duy,dvy,dwy,dpy,dry
    real(wp) :: duz,dvz,dwz,dpz,drz
    !-------------------------------------------------------------------------

    ! Pointers for spherical coordinates
    ! ==================================
    ir=>BC_edge(1,1,2)%i_r
    cosphi=>BC_edge(1,1,2)%cosp
    sinphi=>BC_edge(1,1,2)%sinp
    costeta=>BC_edge(1,1,2)%cost
    sinteta=>BC_edge(1,1,2)%sint
    costcosp=>BC_edge(1,1,2)%costcosp
    costsinp=>BC_edge(1,1,2)%costsinp
    sintcosp=>BC_edge(1,1,2)%sintcosp
    sintsinp=>BC_edge(1,1,2)%sintsinp
    
    ! Compute fluctuations
    ! ====================
    do i=1,nghp3
       do j=nymngh-2,ny
          l=j-nymngh
          do k=nz1,nz2
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
          do k=ndz,nfz
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
       do k=ndz,nfz
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
       do k=ndz,nfz
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
       do k=ndz,nfz
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
          do k=ndz,nfz
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
          do k=ndz,nfz
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
       do k=ndz,nfz
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
       do k=ndz,nfz
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
       do k=ndz,nfz
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

    ! Compute rt,ut,vt,wt & pt
    ! ========================
    do i=1,ngh
       do l=1,ngh
          do k=ndz,nfz
             duz = (a7(1)*( uf(i,l,k+1) - uf(i,l,k-1) ) + &
                    a7(2)*( uf(i,l,k+2) - uf(i,l,k-2) ) + &
                    a7(3)*( uf(i,l,k+3) - uf(i,l,k-3) ) ) *idz(k)*costeta(i,l,k)
             dvz = (a7(1)*( vf(i,l,k+1) - vf(i,l,k-1) ) + &
                    a7(2)*( vf(i,l,k+2) - vf(i,l,k-2) ) + &
                    a7(3)*( vf(i,l,k+3) - vf(i,l,k-3) ) ) *idz(k)*costeta(i,l,k)
             dwz = (a7(1)*( wf(i,l,k+1) - wf(i,l,k-1) ) + &
                    a7(2)*( wf(i,l,k+2) - wf(i,l,k-2) ) + &
                    a7(3)*( wf(i,l,k+3) - wf(i,l,k-3) ) ) *idz(k)*costeta(i,l,k)
             dpz = (a7(1)*( pf(i,l,k+1) - pf(i,l,k-1) ) + &
                    a7(2)*( pf(i,l,k+2) - pf(i,l,k-2) ) + &
                    a7(3)*( pf(i,l,k+3) - pf(i,l,k-3) ) ) *idz(k)*costeta(i,l,k)
             drz = (a7(1)*( rf(i,l,k+1) - rf(i,l,k-1) ) + &
                    a7(2)*( rf(i,l,k+2) - rf(i,l,k-2) ) + &
                    a7(3)*( rf(i,l,k+3) - rf(i,l,k-3) ) ) *idz(k)*costeta(i,l,k)

             pt(i,l,k) = vg(i,l,k)*(dpx(i,l,k)+dpy(i,l,k)+dpz+pf(i,l,k)*ir(i,l,k))
             ut(i,l,k) = vg(i,l,k)*(dux(i,l,k)+duy(i,l,k)+duz+uf(i,l,k)*ir(i,l,k))
             vt(i,l,k) = vg(i,l,k)*(dvx(i,l,k)+dvy(i,l,k)+dvz+vf(i,l,k)*ir(i,l,k))
             wt(i,l,k) = vg(i,l,k)*(dwx(i,l,k)+dwy(i,l,k)+dwz+wf(i,l,k)*ir(i,l,k))
             rt(i,l,k) = vg(i,l,k)*(drx(i,l,k)+dry(i,l,k)+drz+rf(i,l,k)*ir(i,l,k))
          enddo
       enddo
    enddo

    ! Update fluxes at each RK step
    ! =============================
    do k=ndz,nfz
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

  end subroutine bc_TD3d_imin_jmax

  !===============================================================================
  module subroutine bc_TD3d_imax_jmin
  !===============================================================================
    !> 3D Tam & Dong's BC: boundary condition at imax-jmin (edge 1,2,1 /right-bottom)
    !> - Cartesian version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l
    real(wp) :: cp,av,c2_
    real(wp), dimension(-2:ngh,1:ngh+3,nz1:nz2) :: rf,uf,vf,wf,pf
    real(wp), dimension(1:ngh,1:ngh,nz) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(1:ngh,1:ngh,nz) :: dux,dvx,dwx,dpx,drx
    real(wp), dimension(1:ngh,1:ngh,nz) :: duy,dvy,dwy,dpy,dry
    real(wp) :: duz,dvz,dwz,dpz,drz
    !-------------------------------------------------------------------------
   
    ! Pointers for spherical coordinates
    ! ==================================
    ir=>BC_edge(1,2,1)%i_r
    cosphi=>BC_edge(1,2,1)%cosp
    sinphi=>BC_edge(1,2,1)%sinp
    costeta=>BC_edge(1,2,1)%cost
    sinteta=>BC_edge(1,2,1)%sint
    costcosp=>BC_edge(1,2,1)%costcosp
    costsinp=>BC_edge(1,2,1)%costsinp
    sintcosp=>BC_edge(1,2,1)%sintcosp
    sintsinp=>BC_edge(1,2,1)%sintsinp
    
    ! Compute fluctuations
    ! ====================
    do i=nxmngh-2,nx
       l=i-nxmngh
       do j=1,nghp3
          do k=nz1,nz2
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
          do k=ndz,nfz
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
          do k=ndz,nfz
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
       do k=ndz,nfz
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
       do k=ndz,nfz
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
       do k=ndz,nfz
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
    do i=nxmnghp1,nx
       l=i-nxmngh
       do k=ndz,nfz
          duy(l,j,k)=(uf(l,j+1,k)-uf(l,j,k))*idy1_jmin*sintsinp(l,j,k)
          dvy(l,j,k)=(vf(l,j+1,k)-vf(l,j,k))*idy1_jmin*sintsinp(l,j,k)
          dwy(l,j,k)=(wf(l,j+1,k)-wf(l,j,k))*idy1_jmin*sintsinp(l,j,k)
          dpy(l,j,k)=(pf(l,j+1,k)-pf(l,j,k))*idy1_jmin*sintsinp(l,j,k)
          dry(l,j,k)=(rf(l,j+1,k)-rf(l,j,k))*idy1_jmin*sintsinp(l,j,k)
       enddo
    enddo

    j=2
    do i=nxmnghp1,nx
       l=i-nxmngh
       do k=ndz,nfz
          duy(l,j,k)=a3(1)*(uf(l,j+1,k)-uf(l,j-1,k))*idy2_jmin*sintsinp(l,j,k)
          dvy(l,j,k)=a3(1)*(vf(l,j+1,k)-vf(l,j-1,k))*idy2_jmin*sintsinp(l,j,k)
          dwy(l,j,k)=a3(1)*(wf(l,j+1,k)-wf(l,j-1,k))*idy2_jmin*sintsinp(l,j,k)
          dpy(l,j,k)=a3(1)*(pf(l,j+1,k)-pf(l,j-1,k))*idy2_jmin*sintsinp(l,j,k)
          dry(l,j,k)=a3(1)*(rf(l,j+1,k)-rf(l,j-1,k))*idy2_jmin*sintsinp(l,j,k)
       enddo
    enddo

    j=3
    do i=nxmnghp1,nx
       l=i-nxmngh
       do k=ndz,nfz
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
       do i=nxmnghp1,nx
          l=i-nxmngh
          do k=ndz,nfz
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

    ! Compute rt,ut,vt,wt & pt
    ! ========================
    do l=1,ngh
       do j=1,ngh
          do k=ndz,nfz
             duz = (a7(1)*( uf(l,j,k+1) - uf(l,j,k-1) ) + &
                    a7(2)*( uf(l,j,k+2) - uf(l,j,k-2) ) + &
                    a7(3)*( uf(l,j,k+3) - uf(l,j,k-3) ) ) *idz(k)*costeta(l,j,k)
             dvz = (a7(1)*( vf(l,j,k+1) - vf(l,j,k-1) ) + &
                    a7(2)*( vf(l,j,k+2) - vf(l,j,k-2) ) + &
                    a7(3)*( vf(l,j,k+3) - vf(l,j,k-3) ) ) *idz(k)*costeta(l,j,k)
             dwz = (a7(1)*( wf(l,j,k+1) - wf(l,j,k-1) ) + &
                    a7(2)*( wf(l,j,k+2) - wf(l,j,k-2) ) + &
                    a7(3)*( wf(l,j,k+3) - wf(l,j,k-3) ) ) *idz(k)*costeta(l,j,k)
             dpz = (a7(1)*( pf(l,j,k+1) - pf(l,j,k-1) ) + &
                    a7(2)*( pf(l,j,k+2) - pf(l,j,k-2) ) + &
                    a7(3)*( pf(l,j,k+3) - pf(l,j,k-3) ) ) *idz(k)*costeta(l,j,k)
             drz = (a7(1)*( rf(l,j,k+1) - rf(l,j,k-1) ) + &
                    a7(2)*( rf(l,j,k+2) - rf(l,j,k-2) ) + &
                    a7(3)*( rf(l,j,k+3) - rf(l,j,k-3) ) ) *idz(k)*costeta(l,j,k)

             pt(l,j,k) = vg(l,j,k)*(dpx(l,j,k)+dpy(l,j,k)+dpz+pf(l,j,k)*ir(l,j,k))
             ut(l,j,k) = vg(l,j,k)*(dux(l,j,k)+duy(l,j,k)+duz+uf(l,j,k)*ir(l,j,k))
             vt(l,j,k) = vg(l,j,k)*(dvx(l,j,k)+dvy(l,j,k)+dvz+vf(l,j,k)*ir(l,j,k))
             wt(l,j,k) = vg(l,j,k)*(dwx(l,j,k)+dwy(l,j,k)+dwz+wf(l,j,k)*ir(l,j,k))
             rt(l,j,k) = vg(l,j,k)*(drx(l,j,k)+dry(l,j,k)+drz+rf(l,j,k)*ir(l,j,k))
          enddo
       enddo
    enddo

    ! Update fluxes at each RK step
    ! =============================
    do k=ndz,nfz
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

  end subroutine bc_TD3d_imax_jmin

  !===============================================================================
  module subroutine bc_TD3d_imax_jmax
  !===============================================================================
    !> 3D Tam & Dong's BC: boundary condition at imax-jmax (edge 1,2,2 /right-top)
    !> - Cartesian version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l,m
    real(wp) :: cp,av,c2_
    real(wp), dimension(-2:ngh,-2:ngh,nz1:nz2) :: rf,uf,vf,wf,pf
    real(wp), dimension(1:ngh,1:ngh,nz) :: rt,ut,vt,wt,pt,vg
    real(wp), dimension(1:ngh,1:ngh,nz) :: dux,dvx,dwx,dpx,drx
    real(wp), dimension(1:ngh,1:ngh,nz) :: duy,dvy,dwy,dpy,dry
    real(wp) :: duz,dvz,dwz,dpz,drz
    !-------------------------------------------------------------------------
    
    ! Pointers for spherical coordinates
    ! ==================================
    ir=>BC_edge(1,2,2)%i_r
    cosphi=>BC_edge(1,2,2)%cosp
    sinphi=>BC_edge(1,2,2)%sinp
    costeta=>BC_edge(1,2,2)%cost
    sinteta=>BC_edge(1,2,2)%sint
    costcosp=>BC_edge(1,2,2)%costcosp
    costsinp=>BC_edge(1,2,2)%costsinp
    sintcosp=>BC_edge(1,2,2)%sintcosp
    sintsinp=>BC_edge(1,2,2)%sintsinp
    
    ! Compute fluctuations
    ! ====================
    do i=nxmngh-2,nx
       l=i-nxmngh
       do j=nymngh-2,ny
          m=j-nymngh
          do k=nz1,nz2
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
          do k=ndz,nfz
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
          do k=ndz,nfz
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
       do k=ndz,nfz
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
       do k=ndz,nfz
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
       do k=ndz,nfz
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
          do k=ndz,nfz
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
       do k=ndz,nfz
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
       do k=ndz,nfz
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
       do k=ndz,nfz
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
              
    ! Compute rt,ut,vt,wt & pt
    ! ========================
    do l=1,ngh
       do m=1,ngh
          do k=ndz,nfz
             duz = (a7(1)*( uf(l,m,k+1) - uf(l,m,k-1) ) + &
                    a7(2)*( uf(l,m,k+2) - uf(l,m,k-2) ) + &
                    a7(3)*( uf(l,m,k+3) - uf(l,m,k-3) ) ) *idz(k)*costeta(l,m,k)
             dvz = (a7(1)*( vf(l,m,k+1) - vf(l,m,k-1) ) + &
                    a7(2)*( vf(l,m,k+2) - vf(l,m,k-2) ) + &
                    a7(3)*( vf(l,m,k+3) - vf(l,m,k-3) ) ) *idz(k)*costeta(l,m,k)
             dwz = (a7(1)*( wf(l,m,k+1) - wf(l,m,k-1) ) + &
                    a7(2)*( wf(l,m,k+2) - wf(l,m,k-2) ) + &
                    a7(3)*( wf(l,m,k+3) - wf(l,m,k-3) ) ) *idz(k)*costeta(l,m,k)
             dpz = (a7(1)*( pf(l,m,k+1) - pf(l,m,k-1) ) + &
                    a7(2)*( pf(l,m,k+2) - pf(l,m,k-2) ) + &
                    a7(3)*( pf(l,m,k+3) - pf(l,m,k-3) ) ) *idz(k)*costeta(l,m,k)
             drz = (a7(1)*( rf(l,m,k+1) - rf(l,m,k-1) ) + &
                    a7(2)*( rf(l,m,k+2) - rf(l,m,k-2) ) + &
                    a7(3)*( rf(l,m,k+3) - rf(l,m,k-3) ) ) *idz(k)*costeta(l,m,k)

             pt(l,m,k) = vg(l,m,k)*(dpx(l,m,k)+dpy(l,m,k)+dpz+pf(l,m,k)*ir(l,m,k))
             ut(l,m,k) = vg(l,m,k)*(dux(l,m,k)+duy(l,m,k)+duz+uf(l,m,k)*ir(l,m,k))
             vt(l,m,k) = vg(l,m,k)*(dvx(l,m,k)+dvy(l,m,k)+dvz+vf(l,m,k)*ir(l,m,k))
             wt(l,m,k) = vg(l,m,k)*(dwx(l,m,k)+dwy(l,m,k)+dwz+wf(l,m,k)*ir(l,m,k))
             rt(l,m,k) = vg(l,m,k)*(drx(l,m,k)+dry(l,m,k)+drz+rf(l,m,k)*ir(l,m,k))
          enddo
       enddo
    enddo

    ! Update fluxes at each RK step
    ! =============================
    do k=ndz,nfz
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

  end subroutine bc_TD3d_imax_jmax

end submodule smod_TamDong3d_edgez
