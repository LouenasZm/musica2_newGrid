!===============================================================================
submodule (mod_TamDong2d) smod_TamDong2d_faces
!===============================================================================
  !> author: XG
  !> date: February 2020 - modif January 2022
  !> Non-reflecting boundary conditions of Tam and Dong
  !> - 2D (periodic) Cartesian version - routines for faces
!=============================================================================== 

contains

  !===============================================================================
  module subroutine bc_TD2d_imin
  !===============================================================================
    !> 2D Tam & Dong's BC: boundary condition at imin (left) - Cartesian version -
  !===============================================================================
    use mod_eigenmode
    use mod_RFM
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: cp,av,inn1,ntm1
    real(wp) :: du,dv,dw,dp,dr
    real(wp), dimension(1:ngh+3,ny1:ny2,nz) :: rf,uf,vf,wf,pf
    real(wp), dimension(ngh,ny,nz) :: rt,ut,vt,wt,pt,vg,c2_
    !-------------------------------------------------------------------------
    ! eigenmode or RFM disturbances
    real(wp) :: pt_in_,ut_in_,vt_in_,wt_in_,rt_in_
    real(wp) :: dp_in,du_in,dv_in,dw_in,dr_in
    real(wp), dimension(ngh,ny,nz) :: ut_in,vt_in,wt_in
    !-------------------------------------------------------------------------
    ! param vortex -> TEMPORARY, to be changed
    real(wp) :: exm,exm2,xa,ya,x0,ar,trk
    ! to directly impose u,v,w TO BE CHANGED
    ! real(wp) :: u_in,v_in,w_in
    real(wp) :: dutdx_in,dvtdx_in,dptdx_in,dutdy_in,dvtdy_in,dptdy_in
    real(wp) :: dutdt_in,dvtdt_in,dptdt_in,ampl
    !-------------------------------------------------------------------------
    ! added for is_mean_ref
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
          eps=0.1_wp
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
          du = ( a06(1)*uf(1,j,k)+a06(2)*uf(2,j,k) &
               + a06(3)*uf(3,j,k)+a06(4)*uf(4,j,k) &
               + a06(5)*uf(5,j,k)+a06(6)*uf(6,j,k) &
               + a06(7)*uf(7,j,k) ) *idx(i)*cosphi(i,j)
          dv = ( a06(1)*vf(1,j,k)+a06(2)*vf(2,j,k) &
               + a06(3)*vf(3,j,k)+a06(4)*vf(4,j,k) &
               + a06(5)*vf(5,j,k)+a06(6)*vf(6,j,k) &
               + a06(7)*vf(7,j,k) ) *idx(i)*cosphi(i,j)
          dw = ( a06(1)*wf(1,j,k)+a06(2)*wf(2,j,k) &
               + a06(3)*wf(3,j,k)+a06(4)*wf(4,j,k) &
               + a06(5)*wf(5,j,k)+a06(6)*wf(6,j,k) &
               + a06(7)*wf(7,j,k) ) *idx(i)*cosphi(i,j)
          dp = ( a06(1)*pf(1,j,k)+a06(2)*pf(2,j,k) &
               + a06(3)*pf(3,j,k)+a06(4)*pf(4,j,k) &
               + a06(5)*pf(5,j,k)+a06(6)*pf(6,j,k) &
               + a06(7)*pf(7,j,k) ) *idx(i)*cosphi(i,j)
          dr = ( a06(1)*rf(1,j,k)+a06(2)*rf(2,j,k) &
               + a06(3)*rf(3,j,k)+a06(4)*rf(4,j,k) &
               + a06(5)*rf(5,j,k)+a06(6)*rf(6,j,k) &
               + a06(7)*rf(7,j,k) ) *idx(i)*cosphi(i,j)

          du = du + (a7(1)*( uf(i,j+1,k) - uf(i,j-1,k) ) + &
                     a7(2)*( uf(i,j+2,k) - uf(i,j-2,k) ) + &
                     a7(3)*( uf(i,j+3,k) - uf(i,j-3,k) ) ) *idy(j)*sinphi(i,j)
          dv = dv + (a7(1)*( vf(i,j+1,k) - vf(i,j-1,k) ) + &
                     a7(2)*( vf(i,j+2,k) - vf(i,j-2,k) ) + &
                     a7(3)*( vf(i,j+3,k) - vf(i,j-3,k) ) ) *idy(j)*sinphi(i,j)
          dw = dw + (a7(1)*( wf(i,j+1,k) - wf(i,j-1,k) ) + &
                     a7(2)*( wf(i,j+2,k) - wf(i,j-2,k) ) + &
                     a7(3)*( wf(i,j+3,k) - wf(i,j-3,k) ) ) *idy(j)*sinphi(i,j)
          dp = dp + (a7(1)*( pf(i,j+1,k) - pf(i,j-1,k) ) + &
                     a7(2)*( pf(i,j+2,k) - pf(i,j-2,k) ) + &
                     a7(3)*( pf(i,j+3,k) - pf(i,j-3,k) ) ) *idy(j)*sinphi(i,j)
          dr = dr + (a7(1)*( rf(i,j+1,k) - rf(i,j-1,k) ) + &
                     a7(2)*( rf(i,j+2,k) - rf(i,j-2,k) ) + &
                     a7(3)*( rf(i,j+3,k) - rf(i,j-3,k) ) ) *idy(j)*sinphi(i,j)

          pt(i,j,k) = vg(i,j,k) * (dp+pf(i,j,k)*ir(i,j))
          ut(i,j,k) = vg(i,j,k) * (du+uf(i,j,k)*ir(i,j))
          vt(i,j,k) = vg(i,j,k) * (dv+vf(i,j,k)*ir(i,j))
          wt(i,j,k) = vg(i,j,k) * (dw+wf(i,j,k)*ir(i,j))
          rt(i,j,k) = vg(i,j,k) * (dr+rf(i,j,k)*ir(i,j))
       enddo
    enddo

    i=2
    do j=ndy,nfy
       do k=1,nz
          du = ( a15(1)*uf(1,j,k)+a15(2)*uf(2,j,k) &
               + a15(3)*uf(3,j,k)+a15(4)*uf(4,j,k) &
               + a15(5)*uf(5,j,k)+a15(6)*uf(6,j,k) &
               + a15(7)*uf(7,j,k) ) *idx(i)*cosphi(i,j)
          dv = ( a15(1)*vf(1,j,k)+a15(2)*vf(2,j,k) &
               + a15(3)*vf(3,j,k)+a15(4)*vf(4,j,k) &
               + a15(5)*vf(5,j,k)+a15(6)*vf(6,j,k) &
               + a15(7)*vf(7,j,k) ) *idx(i)*cosphi(i,j)
          dw = ( a15(1)*wf(1,j,k)+a15(2)*wf(2,j,k) &
               + a15(3)*wf(3,j,k)+a15(4)*wf(4,j,k) &
               + a15(5)*wf(5,j,k)+a15(6)*wf(6,j,k) &
               + a15(7)*wf(7,j,k) ) *idx(i)*cosphi(i,j)
          dp = ( a15(1)*pf(1,j,k)+a15(2)*pf(2,j,k) &
               + a15(3)*pf(3,j,k)+a15(4)*pf(4,j,k) &
               + a15(5)*pf(5,j,k)+a15(6)*pf(6,j,k) &
               + a15(7)*pf(7,j,k) ) *idx(i)*cosphi(i,j)
          dr = ( a15(1)*rf(1,j,k)+a15(2)*rf(2,j,k) &
               + a15(3)*rf(3,j,k)+a15(4)*rf(4,j,k) &
               + a15(5)*rf(5,j,k)+a15(6)*rf(6,j,k) &
               + a15(7)*rf(7,j,k) ) *idx(i)*cosphi(i,j)

          du = du + (a7(1)*( uf(i,j+1,k) - uf(i,j-1,k) ) + &
                     a7(2)*( uf(i,j+2,k) - uf(i,j-2,k) ) + &
                     a7(3)*( uf(i,j+3,k) - uf(i,j-3,k) ) ) *idy(j)*sinphi(i,j)
          dv = dv + (a7(1)*( vf(i,j+1,k) - vf(i,j-1,k) ) + &
                     a7(2)*( vf(i,j+2,k) - vf(i,j-2,k) ) + &
                     a7(3)*( vf(i,j+3,k) - vf(i,j-3,k) ) ) *idy(j)*sinphi(i,j)
          dw = dw + (a7(1)*( wf(i,j+1,k) - wf(i,j-1,k) ) + &
                     a7(2)*( wf(i,j+2,k) - wf(i,j-2,k) ) + &
                     a7(3)*( wf(i,j+3,k) - wf(i,j-3,k) ) ) *idy(j)*sinphi(i,j)
          dp = dp + (a7(1)*( pf(i,j+1,k) - pf(i,j-1,k) ) + &
                     a7(2)*( pf(i,j+2,k) - pf(i,j-2,k) ) + &
                     a7(3)*( pf(i,j+3,k) - pf(i,j-3,k) ) ) *idy(j)*sinphi(i,j)
          dr = dr + (a7(1)*( rf(i,j+1,k) - rf(i,j-1,k) ) + &
                     a7(2)*( rf(i,j+2,k) - rf(i,j-2,k) ) + &
                     a7(3)*( rf(i,j+3,k) - rf(i,j-3,k) ) ) *idy(j)*sinphi(i,j)

          pt(i,j,k) = vg(i,j,k) * (dp+pf(i,j,k)*ir(i,j))
          ut(i,j,k) = vg(i,j,k) * (du+uf(i,j,k)*ir(i,j))
          vt(i,j,k) = vg(i,j,k) * (dv+vf(i,j,k)*ir(i,j))
          wt(i,j,k) = vg(i,j,k) * (dw+wf(i,j,k)*ir(i,j))
          rt(i,j,k) = vg(i,j,k) * (dr+rf(i,j,k)*ir(i,j))
       enddo
    enddo

    i=3
    do j=ndy,nfy
       do k=1,nz
          du = ( a24(1)*uf(1,j,k)+a24(2)*uf(2,j,k) &
               + a24(3)*uf(3,j,k)+a24(4)*uf(4,j,k) &
               + a24(5)*uf(5,j,k)+a24(6)*uf(6,j,k) &
               + a24(7)*uf(7,j,k) ) *idx(i)*cosphi(i,j)
          dv = ( a24(1)*vf(1,j,k)+a24(2)*vf(2,j,k) &
               + a24(3)*vf(3,j,k)+a24(4)*vf(4,j,k) &
               + a24(5)*vf(5,j,k)+a24(6)*vf(6,j,k) &
               + a24(7)*vf(7,j,k) ) *idx(i)*cosphi(i,j)
          dw = ( a24(1)*wf(1,j,k)+a24(2)*wf(2,j,k) &
               + a24(3)*wf(3,j,k)+a24(4)*wf(4,j,k) &
               + a24(5)*wf(5,j,k)+a24(6)*wf(6,j,k) &
               + a24(7)*wf(7,j,k) ) *idx(i)*cosphi(i,j)
          dp = ( a24(1)*pf(1,j,k)+a24(2)*pf(2,j,k) &
               + a24(3)*pf(3,j,k)+a24(4)*pf(4,j,k) &
               + a24(5)*pf(5,j,k)+a24(6)*pf(6,j,k) &
               + a24(7)*pf(7,j,k) ) *idx(i)*cosphi(i,j)
          dr = ( a24(1)*rf(1,j,k)+a24(2)*rf(2,j,k) &
               + a24(3)*rf(3,j,k)+a24(4)*rf(4,j,k) &
               + a24(5)*rf(5,j,k)+a24(6)*rf(6,j,k) &
               + a24(7)*rf(7,j,k) ) *idx(i)*cosphi(i,j)

          du = du + (a7(1)*( uf(i,j+1,k) - uf(i,j-1,k) ) + &
                     a7(2)*( uf(i,j+2,k) - uf(i,j-2,k) ) + &
                     a7(3)*( uf(i,j+3,k) - uf(i,j-3,k) ) ) *idy(j)*sinphi(i,j)
          dv = dv + (a7(1)*( vf(i,j+1,k) - vf(i,j-1,k) ) + &
                     a7(2)*( vf(i,j+2,k) - vf(i,j-2,k) ) + &
                     a7(3)*( vf(i,j+3,k) - vf(i,j-3,k) ) ) *idy(j)*sinphi(i,j)
          dw = dw + (a7(1)*( wf(i,j+1,k) - wf(i,j-1,k) ) + &
                     a7(2)*( wf(i,j+2,k) - wf(i,j-2,k) ) + &
                     a7(3)*( wf(i,j+3,k) - wf(i,j-3,k) ) ) *idy(j)*sinphi(i,j)
          dp = dp + (a7(1)*( pf(i,j+1,k) - pf(i,j-1,k) ) + &
                     a7(2)*( pf(i,j+2,k) - pf(i,j-2,k) ) + &
                     a7(3)*( pf(i,j+3,k) - pf(i,j-3,k) ) ) *idy(j)*sinphi(i,j)
          dr = dr + (a7(1)*( rf(i,j+1,k) - rf(i,j-1,k) ) + &
                     a7(2)*( rf(i,j+2,k) - rf(i,j-2,k) ) + &
                     a7(3)*( rf(i,j+3,k) - rf(i,j-3,k) ) ) *idy(j)*sinphi(i,j)

          pt(i,j,k) = vg(i,j,k) * (dp+pf(i,j,k)*ir(i,j))
          ut(i,j,k) = vg(i,j,k) * (du+uf(i,j,k)*ir(i,j))
          vt(i,j,k) = vg(i,j,k) * (dv+vf(i,j,k)*ir(i,j))
          wt(i,j,k) = vg(i,j,k) * (dw+wf(i,j,k)*ir(i,j))
          rt(i,j,k) = vg(i,j,k) * (dr+rf(i,j,k)*ir(i,j))
       enddo
    enddo

    do i=4,ngh
       do j=ndy,nfy
          do k=1,nz
             du = ( a7(1)* (uf(i+1,j,k)-uf(i-1,j,k)) &
                  + a7(2)* (uf(i+2,j,k)-uf(i-2,j,k)) &
                  + a7(3)* (uf(i+3,j,k)-uf(i-3,j,k)) ) *idx(i)*cosphi(i,j)
             dv = ( a7(1)* (vf(i+1,j,k)-vf(i-1,j,k)) &
                  + a7(2)* (vf(i+2,j,k)-vf(i-2,j,k)) &
                  + a7(3)* (vf(i+3,j,k)-vf(i-3,j,k)) ) *idx(i)*cosphi(i,j)
             dw = ( a7(1)* (wf(i+1,j,k)-wf(i-1,j,k)) &
                  + a7(2)* (wf(i+2,j,k)-wf(i-2,j,k)) &
                  + a7(3)* (wf(i+3,j,k)-wf(i-3,j,k)) ) *idx(i)*cosphi(i,j)
             dp = ( a7(1)* (pf(i+1,j,k)-pf(i-1,j,k)) &
                  + a7(2)* (pf(i+2,j,k)-pf(i-2,j,k)) &
                  + a7(3)* (pf(i+3,j,k)-pf(i-3,j,k)) ) *idx(i)*cosphi(i,j)
             dr = ( a7(1)* (rf(i+1,j,k)-rf(i-1,j,k)) &
                  + a7(2)* (rf(i+2,j,k)-rf(i-2,j,k)) &
                  + a7(3)* (rf(i+3,j,k)-rf(i-3,j,k)) ) *idx(i)*cosphi(i,j)

             du = du + (a7(1)*( uf(i,j+1,k) - uf(i,j-1,k) ) + &
                        a7(2)*( uf(i,j+2,k) - uf(i,j-2,k) ) + &
                        a7(3)*( uf(i,j+3,k) - uf(i,j-3,k) ) ) *idy(j)*sinphi(i,j)
             dv = dv + (a7(1)*( vf(i,j+1,k) - vf(i,j-1,k) ) + &
                        a7(2)*( vf(i,j+2,k) - vf(i,j-2,k) ) + &
                        a7(3)*( vf(i,j+3,k) - vf(i,j-3,k) ) ) *idy(j)*sinphi(i,j)
             dw = dw + (a7(1)*( wf(i,j+1,k) - wf(i,j-1,k) ) + &
                        a7(2)*( wf(i,j+2,k) - wf(i,j-2,k) ) + &
                        a7(3)*( wf(i,j+3,k) - wf(i,j-3,k) ) ) *idy(j)*sinphi(i,j)
             dp = dp + (a7(1)*( pf(i,j+1,k) - pf(i,j-1,k) ) + &
                        a7(2)*( pf(i,j+2,k) - pf(i,j-2,k) ) + &
                        a7(3)*( pf(i,j+3,k) - pf(i,j-3,k) ) ) *idy(j)*sinphi(i,j)
             dr = dr + (a7(1)*( rf(i,j+1,k) - rf(i,j-1,k) ) + &
                        a7(2)*( rf(i,j+2,k) - rf(i,j-2,k) ) + &
                        a7(3)*( rf(i,j+3,k) - rf(i,j-3,k) ) ) *idy(j)*sinphi(i,j)

             pt(i,j,k) = vg(i,j,k)*(dp+pf(i,j,k)*ir(i,j))
             ut(i,j,k) = vg(i,j,k)*(du+uf(i,j,k)*ir(i,j))
             vt(i,j,k) = vg(i,j,k)*(dv+vf(i,j,k)*ir(i,j))
             wt(i,j,k) = vg(i,j,k)*(dw+wf(i,j,k)*ir(i,j))
             rt(i,j,k) = vg(i,j,k)*(dr+rf(i,j,k)*ir(i,j))
          enddo
       enddo
    enddo

    if (is_eigenmode) then
       do i=1,ngh
          do j=ndy,nfy
             do k=1,nz
                call eig_disturb2_imin(i,j,k,pt_in_,ut_in_,vt_in_,wt_in_,rt_in_,dp_in,du_in,dv_in,dw_in,dr_in)
                pt(i,j,k) = pt(i,j,k) - vg(i,j,k)*dp_in - pt_in_
                ut(i,j,k) = ut(i,j,k) - vg(i,j,k)*du_in - ut_in_
                vt(i,j,k) = vt(i,j,k) - vg(i,j,k)*dv_in - vt_in_
                wt(i,j,k) = wt(i,j,k) - vg(i,j,k)*dw_in - wt_in_
                rt(i,j,k) = rt(i,j,k) - vg(i,j,k)*dr_in - rt_in_
!!$                pt(i,j,k) = - pt_in_
!!$                ut(i,j,k) = - ut_in_
!!$                vt(i,j,k) = - vt_in_
!!$                wt(i,j,k) = - wt_in_
!!$                rt(i,j,k) = - rt_in_
             enddo
          enddo
       enddo
    endif

    if (is_RFM) then
       call disturb_inlet_RFM_TamDong_imin(vg,ut_in,vt_in,wt_in)
       do i=1,ngh
          do j=ndy,nfy
             do k=1,nz
                ! ut_in = vg(i,j,k)*du_in + ut_in
                ut(i,j,k) = ut(i,j,k) - ut_in(i,j,k)
                vt(i,j,k) = vt(i,j,k) - vt_in(i,j,k)
                wt(i,j,k) = wt(i,j,k) - wt_in(i,j,k)
             enddo
          enddo
       enddo
    endif


    ! if (is_RFM) then
    !    do i=1,ngh
    !       do j=ndy,nfy
    !          do k=1,nz
    !             ! ! Test 1
    !             ! if ((j+coord(2)*ny<=ngy-30).and.(j+coord(2)*ny>=30)) then
    !             !    if (is_2D) then
    !             !       call disturb_inlet_RFM_TamDong(i,j,k,ut_in,vt_in,wt_in,du_in,dv_in,dw_in)
    !             !       ut(i,j,k) = ut(i,j,k) - vg(i,j,k)*du_in - ut_in
    !             !       vt(i,j,k) = vt(i,j,k) - vg(i,j,k)*dv_in - vt_in
    !             !       wt(i,j,k) = wt(i,j,k) - vg(i,j,k)*dw_in - wt_in
    !             !    else if ((k+coord(3)*nz>=30).and.(k+coord(3)*nz<=ngz-30)) then
    !             !    ! else
    !             !       call disturb_inlet_RFM_TamDong(i,j,k,ut_in,vt_in,wt_in,du_in,dv_in,dw_in)
    !             !       ut(i,j,k) = ut(i,j,k) - vg(i,j,k)*du_in - ut_in
    !             !       vt(i,j,k) = vt(i,j,k) - vg(i,j,k)*dv_in - vt_in
    !             !       wt(i,j,k) = wt(i,j,k) - vg(i,j,k)*dw_in - wt_in
    !             !    else if ((k+coord(3)*nz>=40).and.(k+coord(3)*nz<50)) then
    !             !       call disturb_inlet_RFM_TamDong(i,j,k,ut_in,vt_in,wt_in,du_in,dv_in,dw_in)
    !             !       ut(i,j,k) = ut(i,j,k) - (vg(i,j,k)*du_in + ut_in)*(pi/2-atan(-abs(dble(k-49))))
    !             !       vt(i,j,k) = vt(i,j,k) - (vg(i,j,k)*dv_in + vt_in)*(pi/2-atan(-abs(dble(k-49))))
    !             !       wt(i,j,k) = wt(i,j,k) - (vg(i,j,k)*dw_in + wt_in)*(pi/2-atan(-abs(dble(k-49))))
    !             !    else if ((k+coord(3)*nz>ngz-50).and.(k+coord(3)*nz<=ngz-40)) then
    !             !       call disturb_inlet_RFM_TamDong(i,j,k,ut_in,vt_in,wt_in,du_in,dv_in,dw_in)
    !             !       ut(i,j,k) = ut(i,j,k) - (vg(i,j,k)*du_in + ut_in)*(pi/2-atan(-abs(dble(k-(ngz-49)))))
    !             !       vt(i,j,k) = vt(i,j,k) - (vg(i,j,k)*dv_in + vt_in)*(pi/2-atan(-abs(dble(k-(ngz-49)))))
    !             !       wt(i,j,k) = wt(i,j,k) - (vg(i,j,k)*dw_in + wt_in)*(pi/2-atan(-abs(dble(k-(ngz-49)))))
    !             !    endif
    !             ! else if ((j+coord(2)*ny>ngy-50).and.(j+coord(2)*ny<=ngy-40)) then
    !             !    if (is_2D) then
    !             !       call disturb_inlet_RFM_TamDong(i,j,k,ut_in,vt_in,wt_in,du_in,dv_in,dw_in)
    !             !       ut(i,j,k) = ut(i,j,k) - (vg(i,j,k)*du_in + ut_in)*(pi/2-atan(-abs(dble(j-(ngy-49)))))
    !             !       vt(i,j,k) = vt(i,j,k) - (vg(i,j,k)*dv_in + vt_in)*(pi/2-atan(-abs(dble(j-(ngy-49)))))
    !             !       wt(i,j,k) = wt(i,j,k) - (vg(i,j,k)*dw_in + wt_in)*(pi/2-atan(-abs(dble(j-(ngy-49)))))
    !             !    else if ((k+coord(3)*nz>=50).and.(k+coord(3)*nz<=ngz-50)) then
    !             !       call disturb_inlet_RFM_TamDong(i,j,k,ut_in,vt_in,wt_in,du_in,dv_in,dw_in)
    !             !       ut(i,j,k) = ut(i,j,k) - (vg(i,j,k)*du_in + ut_in)*(pi/2-atan(-abs(dble(j-(ngy-49)))))
    !             !       vt(i,j,k) = vt(i,j,k) - (vg(i,j,k)*dv_in + vt_in)*(pi/2-atan(-abs(dble(j-(ngy-49)))))
    !             !       wt(i,j,k) = wt(i,j,k) - (vg(i,j,k)*dw_in + wt_in)*(pi/2-atan(-abs(dble(j-(ngy-49)))))
    !             !    else if ((k+coord(3)*nz>=40).and.(k+coord(3)*nz<50)) then
    !             !       call disturb_inlet_RFM_TamDong(i,j,k,ut_in,vt_in,wt_in,du_in,dv_in,dw_in)
    !             !       ut(i,j,k) = ut(i,j,k) - (vg(i,j,k)*du_in + ut_in)*(pi/2-atan(-abs(dble(j-(ngy-49)))))*(pi/2-atan(-abs(dble(k-49))))
    !             !       vt(i,j,k) = vt(i,j,k) - (vg(i,j,k)*dv_in + vt_in)*(pi/2-atan(-abs(dble(j-(ngy-49)))))*(pi/2-atan(-abs(dble(k-49))))
    !             !       wt(i,j,k) = wt(i,j,k) - (vg(i,j,k)*dw_in + wt_in)*(pi/2-atan(-abs(dble(j-(ngy-49)))))*(pi/2-atan(-abs(dble(k-49))))
    !             !    else if ((k+coord(3)*nz>ngz-50).and.(k+coord(3)*nz<=ngz-40)) then
    !             !       call disturb_inlet_RFM_TamDong(i,j,k,ut_in,vt_in,wt_in,du_in,dv_in,dw_in)
    !             !       ut(i,j,k) = ut(i,j,k) - (vg(i,j,k)*du_in + ut_in)*(pi/2-atan(-abs(dble(j-(ngy-49)))))*(pi/2-atan(-abs(dble(k-(ngz-49)))))
    !             !       vt(i,j,k) = vt(i,j,k) - (vg(i,j,k)*dv_in + vt_in)*(pi/2-atan(-abs(dble(j-(ngy-49)))))*(pi/2-atan(-abs(dble(k-(ngz-49)))))
    !             !       wt(i,j,k) = wt(i,j,k) - (vg(i,j,k)*dw_in + wt_in)*(pi/2-atan(-abs(dble(j-(ngy-49)))))*(pi/2-atan(-abs(dble(k-(ngz-49)))))
    !             !    endif
    !             ! else if ((j+coord(2)*ny>=40).and.(j+coord(2)*ny<50)) then
    !             !    if (is_2D) then
    !             !       call disturb_inlet_RFM_TamDong(i,j,k,ut_in,vt_in,wt_in,du_in,dv_in,dw_in)
    !             !       ut(i,j,k) = ut(i,j,k) - (vg(i,j,k)*du_in + ut_in)*(pi/2-atan(-abs(dble(j-49))))
    !             !       vt(i,j,k) = vt(i,j,k) - (vg(i,j,k)*dv_in + vt_in)*(pi/2-atan(-abs(dble(j-49))))
    !             !       wt(i,j,k) = wt(i,j,k) - (vg(i,j,k)*dw_in + wt_in)*(pi/2-atan(-abs(dble(j-49))))
    !             !    else if ((k+coord(3)*nz>=50).and.(k+coord(3)*nz<=ngz-50)) then
    !             !       call disturb_inlet_RFM_TamDong(i,j,k,ut_in,vt_in,wt_in,du_in,dv_in,dw_in)
    !             !       ut(i,j,k) = ut(i,j,k) - (vg(i,j,k)*du_in + ut_in)*(pi/2-atan(-abs(dble(j-49))))
    !             !       vt(i,j,k) = vt(i,j,k) - (vg(i,j,k)*dv_in + vt_in)*(pi/2-atan(-abs(dble(j-49))))
    !             !       wt(i,j,k) = wt(i,j,k) - (vg(i,j,k)*dw_in + wt_in)*(pi/2-atan(-abs(dble(j-49))))
    !             !    else if ((k+coord(3)*nz>=40).and.(k+coord(3)*nz<50)) then
    !             !       call disturb_inlet_RFM_TamDong(i,j,k,ut_in,vt_in,wt_in,du_in,dv_in,dw_in)
    !             !       ut(i,j,k) = ut(i,j,k) - (vg(i,j,k)*du_in + ut_in)*(pi/2-atan(-abs(dble(j-49))))*(pi/2-atan(-abs(dble(k-49))))
    !             !       vt(i,j,k) = vt(i,j,k) - (vg(i,j,k)*dv_in + vt_in)*(pi/2-atan(-abs(dble(j-49))))*(pi/2-atan(-abs(dble(k-49))))
    !             !       wt(i,j,k) = wt(i,j,k) - (vg(i,j,k)*dw_in + wt_in)*(pi/2-atan(-abs(dble(j-49))))*(pi/2-atan(-abs(dble(k-49))))
    !             !    else if ((k+coord(3)*nz>ngz-50).and.(k+coord(3)*nz<=ngz-40)) then
    !             !       call disturb_inlet_RFM_TamDong(i,j,k,ut_in,vt_in,wt_in,du_in,dv_in,dw_in)
    !             !       ut(i,j,k) = ut(i,j,k) - (vg(i,j,k)*du_in + ut_in)*(pi/2-atan(-abs(dble(j-49))))*(pi/2-atan(-abs(dble(k-(ngz-49)))))
    !             !       vt(i,j,k) = vt(i,j,k) - (vg(i,j,k)*dv_in + vt_in)*(pi/2-atan(-abs(dble(j-49))))*(pi/2-atan(-abs(dble(k-(ngz-49)))))
    !             !       wt(i,j,k) = wt(i,j,k) - (vg(i,j,k)*dw_in + wt_in)*(pi/2-atan(-abs(dble(j-49))))*(pi/2-atan(-abs(dble(k-(ngz-49)))))
    !             !    endif
    !             ! endif

    !             ! ! Test 2
    !             ! if ((is_2D).and.(j+coord(2)*ny<=ngy-5).and.(j+coord(2)*ny>=5)) then
    !             !    call disturb_inlet_RFM_TamDong(i,j,k,ut_in,vt_in,wt_in,du_in,dv_in,dw_in)
    !             !    ut(i,j,k) = ut(i,j,k) - vg(i,j,k)*du_in - ut_in
    !             !    vt(i,j,k) = vt(i,j,k) - vg(i,j,k)*dv_in - vt_in
    !             !    wt(i,j,k) = wt(i,j,k) - vg(i,j,k)*dw_in - wt_in
    !             !    ! ut(i,j,k) = - ut_in
    !             !    ! vt(i,j,k) = - vt_in
    !             !    ! wt(i,j,k) = - wt_in
    !             ! else if ((j+coord(2)*ny<=ngy-5).and.(j+coord(2)*ny>=5).and. &
    !             !     (k+coord(3)*nz<=ngz-5).and.(k+coord(3)*nz>=5)) then
    !             !    call disturb_inlet_RFM_TamDong(i,j,k,ut_in,vt_in,wt_in,du_in,dv_in,dw_in)
    !             !    ut(i,j,k) = ut(i,j,k) - vg(i,j,k)*du_in - ut_in
    !             !    vt(i,j,k) = vt(i,j,k) - vg(i,j,k)*dv_in - vt_in
    !             !    wt(i,j,k) = wt(i,j,k) - vg(i,j,k)*dw_in - wt_in
    !             !    ! ut(i,j,k) = - ut_in
    !             !    ! vt(i,j,k) = - vt_in
    !             !    ! wt(i,j,k) = - wt_in
    !             ! endif

    !             call disturb_inlet_RFM_TamDong(i,j,k,ut_in,vt_in,wt_in,du_in,dv_in,dw_in)
    !             ut(i,j,k) = ut(i,j,k) - vg(i,j,k)*du_in - ut_in
    !             vt(i,j,k) = vt(i,j,k) - vg(i,j,k)*dv_in - vt_in
    !             wt(i,j,k) = wt(i,j,k) - vg(i,j,k)*dw_in - wt_in
    !             ! pt(i,j,k) = 0.0_wp
    !             ! ut(i,j,k) = - ut_in
    !             ! vt(i,j,k) = - vt_in
    !             ! wt(i,j,k) = - wt_in
    !             ! rt(i,j,k) = 0.0_wp

    !          enddo
    !       enddo
    !    enddo
    ! endif

!!$    if ((.not.is_vortex).and.(type_vortex.eq.-1)) then
!!$       ! Gaussian half-width
!!$       ar = log(2.0_wp)/(25.0_wp*deltay**2)
!!$       ! position of vortex center (+ virtual initial position)
!!$       trk = time + ck(irk)*deltat
!!$       x0 = -15.0_wp*deltay + U_ref*trk
!!$       ampl = 1.0_wp
!!$       do i=1,ngh
!!$          do j=ndy,nfy
!!$             do k=1,nz
!!$                xa=(x(i)-x0)
!!$                ya=(y(j)-yg(ngy/2))
!!$                exm=exp(-ar*(xa**2+ya**2))
!!$                exm2=exp(-2.0_wp*ar*(xa**2+ya**2))
!!$
!!$                dutdt_in = ampl*U_ref*2*ar*(xa*ya/deltay)*exm
!!$                dutdx_in = -ampl*2*ar*(xa*ya/deltay)*exm
!!$                dutdy_in = (ampl/deltay)*exm*(1-2*ar*ya**2)
!!$                du_in = ampl*(ya/deltay)*exm*ir(i,j) + dutdx_in*cosphi(i,j) + dutdy_in*sinphi(i,j)
!!$
!!$                dvtdt_in = (ampl*U_ref/deltay)*(1 - 2*ar*xa**2)*exm
!!$                dvtdx_in = (ampl/deltay)*(2*ar*xa**2 - 1)*exm
!!$                dvtdy_in = ampl*(xa*ya/deltay)*2*ar*exm
!!$                dv_in = -ampl*(xa/deltay)*exm*ir(i,j) + dvtdx_in*cosphi(i,j) + dvtdy_in*sinphi(i,j)
!!$
!!$                dptdt_in = -rho_ref*(ampl/deltay)**2 * U_ref*xa*exm2
!!$                dptdx_in = rho_ref*(ampl/deltay)**2 * xa*exm2
!!$                dptdy_in = rho_ref*(ampl/deltay)**2 * ya*exm2
!!$                dp_in = - (rho(i,j,k)*ampl**2)/(4*ar*deltay**2)*exm2*ir(i,j) + dptdx_in*cosphi(i,j) + dptdy_in*sinphi(i,j)
!!$
!!$                ! Case 1: with spatial and time derivatives
!!$                ut(i,j,k) = ut(i,j,k) - vg(i,j,k)*du_in - dutdt_in
!!$                vt(i,j,k) = vt(i,j,k) - vg(i,j,k)*dv_in - dvtdt_in
!!$                pt(i,j,k) = pt(i,j,k) - vg(i,j,k)*dp_in - dptdt_in
!!$                ! Case 2: with only time derivatives
!!$                ! ut(i,j,k) = ut(i,j,k) - dutdt_in
!!$                ! vt(i,j,k) = vt(i,j,k) - dvtdt_in
!!$                ! pt(i,j,k) = pt(i,j,k) - dptdt_in
!!$                ! Case 3: directly imposed
!!$                ! ut(i,j,k) = - dutdt_in
!!$                ! vt(i,j,k) = - dvtdt_in
!!$                ! wt(i,j,k) = - dptdt_in
!!$                ! pt(i,j,k) = 0.0_wp
!!$                ! rt(i,j,k) = 0.0_wp
!!$             enddo
!!$          enddo
!!$       enddo
!!$    endif

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

  end subroutine bc_TD2d_imin

  !===============================================================================
  module subroutine bc_TD2d_imax
  !===============================================================================
    !> 2D Tam & Dong's BC: boundary condition at imax (right) - Cartesian version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l
    real(wp) :: cp,av,inn1,ntm1
    real(wp) :: du,dv,dw,dp,dr
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
             du = ( a7(1)* (uf(l+1,j,k)-uf(l-1,j,k)) &
                  + a7(2)* (uf(l+2,j,k)-uf(l-2,j,k)) &
                  + a7(3)* (uf(l+3,j,k)-uf(l-3,j,k)) ) *idx(i)*cosphi(l,j)
             dv = ( a7(1)* (vf(l+1,j,k)-vf(l-1,j,k)) &
                  + a7(2)* (vf(l+2,j,k)-vf(l-2,j,k)) &
                  + a7(3)* (vf(l+3,j,k)-vf(l-3,j,k)) ) *idx(i)*cosphi(l,j)
             dw = ( a7(1)* (wf(l+1,j,k)-wf(l-1,j,k)) &
                  + a7(2)* (wf(l+2,j,k)-wf(l-2,j,k)) &
                  + a7(3)* (wf(l+3,j,k)-wf(l-3,j,k)) ) *idx(i)*cosphi(l,j)
             dp = ( a7(1)* (pf(l+1,j,k)-pf(l-1,j,k)) &
                  + a7(2)* (pf(l+2,j,k)-pf(l-2,j,k)) &
                  + a7(3)* (pf(l+3,j,k)-pf(l-3,j,k)) ) *idx(i)*cosphi(l,j)
             dr = ( a7(1)* (rf(l+1,j,k)-rf(l-1,j,k)) &
                  + a7(2)* (rf(l+2,j,k)-rf(l-2,j,k)) &
                  + a7(3)* (rf(l+3,j,k)-rf(l-3,j,k)) ) *idx(i)*cosphi(l,j)

             du = du + (a7(1)*( uf(l,j+1,k) - uf(l,j-1,k) ) + &
                        a7(2)*( uf(l,j+2,k) - uf(l,j-2,k) ) + &
                        a7(3)*( uf(l,j+3,k) - uf(l,j-3,k) ) ) *idy(j)*sinphi(l,j)
             dv = dv + (a7(1)*( vf(l,j+1,k) - vf(l,j-1,k) ) + &
                        a7(2)*( vf(l,j+2,k) - vf(l,j-2,k) ) + &
                        a7(3)*( vf(l,j+3,k) - vf(l,j-3,k) ) ) *idy(j)*sinphi(l,j)
             dw = dw + (a7(1)*( wf(l,j+1,k) - wf(l,j-1,k) ) + &
                        a7(2)*( wf(l,j+2,k) - wf(l,j-2,k) ) + &
                        a7(3)*( wf(l,j+3,k) - wf(l,j-3,k) ) ) *idy(j)*sinphi(l,j)
             dp = dp + (a7(1)*( pf(l,j+1,k) - pf(l,j-1,k) ) + &
                        a7(2)*( pf(l,j+2,k) - pf(l,j-2,k) ) + &
                        a7(3)*( pf(l,j+3,k) - pf(l,j-3,k) ) ) *idy(j)*sinphi(l,j)
             dr = dr + (a7(1)*( rf(l,j+1,k) - rf(l,j-1,k) ) + &
                        a7(2)*( rf(l,j+2,k) - rf(l,j-2,k) ) + &
                        a7(3)*( rf(l,j+3,k) - rf(l,j-3,k) ) ) *idy(j)*sinphi(l,j)

             pt(l,j,k) = vg(l,j,k)*(dp+pf(l,j,k)*ir(l,j))
             ut(l,j,k) = vg(l,j,k)*(du+uf(l,j,k)*ir(l,j))
             vt(l,j,k) = vg(l,j,k)*(dv+vf(l,j,k)*ir(l,j))
             wt(l,j,k) = vg(l,j,k)*(dw+wf(l,j,k)*ir(l,j))
             rt(l,j,k) = vg(l,j,k)*(dr+rf(l,j,k)*ir(l,j))
          enddo
       enddo
    enddo

    i=nx-2
    l=i-nxmngh
    do j=ndy,nfy
       do k=1,nz
          du = ( a42(1)*uf(l+2,j,k)+a42(2)*uf(l+1,j,k) &
               + a42(3)*uf(l  ,j,k)+a42(4)*uf(l-1,j,k) &
               + a42(5)*uf(l-2,j,k)+a42(6)*uf(l-3,j,k) &
               + a42(7)*uf(l-4,j,k) ) *idx(i)*cosphi(l,j)
          dv = ( a42(1)*vf(l+2,j,k)+a42(2)*vf(l+1,j,k) &
               + a42(3)*vf(l  ,j,k)+a42(4)*vf(l-1,j,k) &
               + a42(5)*vf(l-2,j,k)+a42(6)*vf(l-3,j,k) &
               + a42(7)*vf(l-4,j,k) ) *idx(i)*cosphi(l,j)
          dw = ( a42(1)*wf(l+2,j,k)+a42(2)*wf(l+1,j,k) &
               + a42(3)*wf(l  ,j,k)+a42(4)*wf(l-1,j,k) &
               + a42(5)*wf(l-2,j,k)+a42(6)*wf(l-3,j,k) &
               + a42(7)*wf(l-4,j,k) ) *idx(i)*cosphi(l,j)
          dp = ( a42(1)*pf(l+2,j,k)+a42(2)*pf(l+1,j,k) &
               + a42(3)*pf(l  ,j,k)+a42(4)*pf(l-1,j,k) &
               + a42(5)*pf(l-2,j,k)+a42(6)*pf(l-3,j,k) &
               + a42(7)*pf(l-4,j,k) ) *idx(i)*cosphi(l,j)
          dr = ( a42(1)*rf(l+2,j,k)+a42(2)*rf(l+1,j,k) &
               + a42(3)*rf(l  ,j,k)+a42(4)*rf(l-1,j,k) &
               + a42(5)*rf(l-2,j,k)+a42(6)*rf(l-3,j,k) &
               + a42(7)*rf(l-4,j,k) ) *idx(i)*cosphi(l,j)
          
          du = du + (a7(1)*( uf(l,j+1,k) - uf(l,j-1,k) ) + &
                     a7(2)*( uf(l,j+2,k) - uf(l,j-2,k) ) + &
                     a7(3)*( uf(l,j+3,k) - uf(l,j-3,k) ) ) *idy(j)*sinphi(l,j)
          dv = dv + (a7(1)*( vf(l,j+1,k) - vf(l,j-1,k) ) + &
                     a7(2)*( vf(l,j+2,k) - vf(l,j-2,k) ) + &
                     a7(3)*( vf(l,j+3,k) - vf(l,j-3,k) ) ) *idy(j)*sinphi(l,j)
          dw = dw + (a7(1)*( wf(l,j+1,k) - wf(l,j-1,k) ) + &
                     a7(2)*( wf(l,j+2,k) - wf(l,j-2,k) ) + &
                     a7(3)*( wf(l,j+3,k) - wf(l,j-3,k) ) ) *idy(j)*sinphi(l,j)
          dp = dp + (a7(1)*( pf(l,j+1,k) - pf(l,j-1,k) ) + &
                     a7(2)*( pf(l,j+2,k) - pf(l,j-2,k) ) + &
                     a7(3)*( pf(l,j+3,k) - pf(l,j-3,k) ) ) *idy(j)*sinphi(l,j)
          dr = dr + (a7(1)*( rf(l,j+1,k) - rf(l,j-1,k) ) + &
                     a7(2)*( rf(l,j+2,k) - rf(l,j-2,k) ) + &
                     a7(3)*( rf(l,j+3,k) - rf(l,j-3,k) ) ) *idy(j)*sinphi(l,j)

          pt(l,j,k) = vg(l,j,k)*(dp+pf(l,j,k)*ir(l,j))
          ut(l,j,k) = vg(l,j,k)*(du+uf(l,j,k)*ir(l,j))
          vt(l,j,k) = vg(l,j,k)*(dv+vf(l,j,k)*ir(l,j))
          wt(l,j,k) = vg(l,j,k)*(dw+wf(l,j,k)*ir(l,j))
          rt(l,j,k) = vg(l,j,k)*(dr+rf(l,j,k)*ir(l,j))
       enddo
    enddo

    i=nx-1
    l=i-nxmngh
    do j=ndy,nfy
       do k=1,nz
          du = ( a51(1)*uf(l+1,j,k)+a51(2)*uf(l  ,j,k) &
               + a51(3)*uf(l-1,j,k)+a51(4)*uf(l-2,j,k) &
               + a51(5)*uf(l-3,j,k)+a51(6)*uf(l-4,j,k) &
               + a51(7)*uf(l-5,j,k) ) *idx(i)*cosphi(l,j)
          dv = ( a51(1)*vf(l+1,j,k)+a51(2)*vf(l  ,j,k) &
               + a51(3)*vf(l-1,j,k)+a51(4)*vf(l-2,j,k) &
               + a51(5)*vf(l-3,j,k)+a51(6)*vf(l-4,j,k) &
               + a51(7)*vf(l-5,j,k) ) *idx(i)*cosphi(l,j)
          dw = ( a51(1)*wf(l+1,j,k)+a51(2)*wf(l  ,j,k) &
               + a51(3)*wf(l-1,j,k)+a51(4)*wf(l-2,j,k) &
               + a51(5)*wf(l-3,j,k)+a51(6)*wf(l-4,j,k) &
               + a51(7)*wf(l-5,j,k) ) *idx(i)*cosphi(l,j)
          dp = ( a51(1)*pf(l+1,j,k)+a51(2)*pf(l  ,j,k) &
               + a51(3)*pf(l-1,j,k)+a51(4)*pf(l-2,j,k) &
               + a51(5)*pf(l-3,j,k)+a51(6)*pf(l-4,j,k) &
               + a51(7)*pf(l-5,j,k) ) *idx(i)*cosphi(l,j)
          dr = ( a51(1)*rf(l+1,j,k)+a51(2)*rf(l  ,j,k) &
               + a51(3)*rf(l-1,j,k)+a51(4)*rf(l-2,j,k) &
               + a51(5)*rf(l-3,j,k)+a51(6)*rf(l-4,j,k) &
               + a51(7)*rf(l-5,j,k) ) *idx(i)*cosphi(l,j)

          du = du + (a7(1)*( uf(l,j+1,k) - uf(l,j-1,k) ) + &
                     a7(2)*( uf(l,j+2,k) - uf(l,j-2,k) ) + &
                     a7(3)*( uf(l,j+3,k) - uf(l,j-3,k) ) ) *idy(j)*sinphi(l,j)
          dv = dv + (a7(1)*( vf(l,j+1,k) - vf(l,j-1,k) ) + &
                     a7(2)*( vf(l,j+2,k) - vf(l,j-2,k) ) + &
                     a7(3)*( vf(l,j+3,k) - vf(l,j-3,k) ) ) *idy(j)*sinphi(l,j)
          dw = dw + (a7(1)*( wf(l,j+1,k) - wf(l,j-1,k) ) + &
                     a7(2)*( wf(l,j+2,k) - wf(l,j-2,k) ) + &
                     a7(3)*( wf(l,j+3,k) - wf(l,j-3,k) ) ) *idy(j)*sinphi(l,j)
          dp = dp + (a7(1)*( pf(l,j+1,k) - pf(l,j-1,k) ) + &
                     a7(2)*( pf(l,j+2,k) - pf(l,j-2,k) ) + &
                     a7(3)*( pf(l,j+3,k) - pf(l,j-3,k) ) ) *idy(j)*sinphi(l,j)
          dr = dr + (a7(1)*( rf(l,j+1,k) - rf(l,j-1,k) ) + &
                     a7(2)*( rf(l,j+2,k) - rf(l,j-2,k) ) + &
                     a7(3)*( rf(l,j+3,k) - rf(l,j-3,k) ) ) *idy(j)*sinphi(l,j)

          pt(l,j,k) = vg(l,j,k)*(dp+pf(l,j,k)*ir(l,j))
          ut(l,j,k) = vg(l,j,k)*(du+uf(l,j,k)*ir(l,j))
          vt(l,j,k) = vg(l,j,k)*(dv+vf(l,j,k)*ir(l,j))
          wt(l,j,k) = vg(l,j,k)*(dw+wf(l,j,k)*ir(l,j))
          rt(l,j,k) = vg(l,j,k)*(dr+rf(l,j,k)*ir(l,j))
       enddo
    enddo

    i=nx
    l=i-nxmngh
    do j=ndy,nfy
       do k=1,nz
          du = ( a60(1)*uf(l  ,j,k)+a60(2)*uf(l-1,j,k) &
               + a60(3)*uf(l-2,j,k)+a60(4)*uf(l-3,j,k) &
               + a60(5)*uf(l-4,j,k)+a60(6)*uf(l-5,j,k) &
               + a60(7)*uf(l-6,j,k) ) *idx(i)*cosphi(l,j)
          dv = ( a60(1)*vf(l  ,j,k)+a60(2)*vf(l-1,j,k) &
               + a60(3)*vf(l-2,j,k)+a60(4)*vf(l-3,j,k) &
               + a60(5)*vf(l-4,j,k)+a60(6)*vf(l-5,j,k) &
               + a60(7)*vf(l-6,j,k) ) *idx(i)*cosphi(l,j)
          dw = ( a60(1)*wf(l  ,j,k)+a60(2)*wf(l-1,j,k) &
               + a60(3)*wf(l-2,j,k)+a60(4)*wf(l-3,j,k) &
               + a60(5)*wf(l-4,j,k)+a60(6)*wf(l-5,j,k) &
               + a60(7)*wf(l-6,j,k) ) *idx(i)*cosphi(l,j)
          dp = ( a60(1)*pf(l  ,j,k)+a60(2)*pf(l-1,j,k) &
               + a60(3)*pf(l-2,j,k)+a60(4)*pf(l-3,j,k) &
               + a60(5)*pf(l-4,j,k)+a60(6)*pf(l-5,j,k) &
               + a60(7)*pf(l-6,j,k) ) *idx(i)*cosphi(l,j)
          dr = ( a60(1)*rf(l  ,j,k)+a60(2)*rf(l-1,j,k) &
               + a60(3)*rf(l-2,j,k)+a60(4)*rf(l-3,j,k) &
               + a60(5)*rf(l-4,j,k)+a60(6)*rf(l-5,j,k) &
               + a60(7)*rf(l-6,j,k) ) *idx(i)*cosphi(l,j)
          
          du = du + (a7(1)*( uf(l,j+1,k) - uf(l,j-1,k) ) + &
                     a7(2)*( uf(l,j+2,k) - uf(l,j-2,k) ) + &
                     a7(3)*( uf(l,j+3,k) - uf(l,j-3,k) ) ) *idy(j)*sinphi(l,j)
          dv = dv + (a7(1)*( vf(l,j+1,k) - vf(l,j-1,k) ) + &
                     a7(2)*( vf(l,j+2,k) - vf(l,j-2,k) ) + &
                     a7(3)*( vf(l,j+3,k) - vf(l,j-3,k) ) ) *idy(j)*sinphi(l,j)
          dw = dw + (a7(1)*( wf(l,j+1,k) - wf(l,j-1,k) ) + &
                     a7(2)*( wf(l,j+2,k) - wf(l,j-2,k) ) + &
                     a7(3)*( wf(l,j+3,k) - wf(l,j-3,k) ) ) *idy(j)*sinphi(l,j)
          dp = dp + (a7(1)*( pf(l,j+1,k) - pf(l,j-1,k) ) + &
                     a7(2)*( pf(l,j+2,k) - pf(l,j-2,k) ) + &
                     a7(3)*( pf(l,j+3,k) - pf(l,j-3,k) ) ) *idy(j)*sinphi(l,j)
          dr = dr + (a7(1)*( rf(l,j+1,k) - rf(l,j-1,k) ) + &
                     a7(2)*( rf(l,j+2,k) - rf(l,j-2,k) ) + &
                     a7(3)*( rf(l,j+3,k) - rf(l,j-3,k) ) ) *idy(j)*sinphi(l,j)

          pt(l,j,k) = vg(l,j,k)*(dp+pf(l,j,k)*ir(l,j))
          ut(l,j,k) = vg(l,j,k)*(du+uf(l,j,k)*ir(l,j))
          vt(l,j,k) = vg(l,j,k)*(dv+vf(l,j,k)*ir(l,j))
          wt(l,j,k) = vg(l,j,k)*(dw+wf(l,j,k)*ir(l,j))
          rt(l,j,k) = vg(l,j,k)*(dr+rf(l,j,k)*ir(l,j))
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

  end subroutine bc_TD2d_imax

  !===============================================================================
  module subroutine bc_TD2d_jmin
  !===============================================================================
    !> 2D Tam & Dong's BC: boundary condition at jmin (bottom) - Cartesian version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: cp,av,inn1,ntm1
    real(wp) :: du,dv,dw,dp,dr
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
          du = ( a7(1)*( uf(i+1,j,k) - uf(i-1,j,k) ) + &
                 a7(2)*( uf(i+2,j,k) - uf(i-2,j,k) ) + &
                 a7(3)*( uf(i+3,j,k) - uf(i-3,j,k) ) ) *idx(i)*cosphi(i,j)
          dv = ( a7(1)*( vf(i+1,j,k) - vf(i-1,j,k) ) + &
                 a7(2)*( vf(i+2,j,k) - vf(i-2,j,k) ) + &
                 a7(3)*( vf(i+3,j,k) - vf(i-3,j,k) ) ) *idx(i)*cosphi(i,j)
          dw = ( a7(1)*( wf(i+1,j,k) - wf(i-1,j,k) ) + &
                 a7(2)*( wf(i+2,j,k) - wf(i-2,j,k) ) + &
                 a7(3)*( wf(i+3,j,k) - wf(i-3,j,k) ) ) *idx(i)*cosphi(i,j)
          dp = ( a7(1)*( pf(i+1,j,k) - pf(i-1,j,k) ) + &
                 a7(2)*( pf(i+2,j,k) - pf(i-2,j,k) ) + &
                 a7(3)*( pf(i+3,j,k) - pf(i-3,j,k) ) ) *idx(i)*cosphi(i,j)
          dr = ( a7(1)*( rf(i+1,j,k) - rf(i-1,j,k) ) + &
                 a7(2)*( rf(i+2,j,k) - rf(i-2,j,k) ) + &
                 a7(3)*( rf(i+3,j,k) - rf(i-3,j,k) ) ) *idx(i)*cosphi(i,j)

          du = du + ( a06(1)*uf(i,1,k)+a06(2)*uf(i,2,k) &
                    + a06(3)*uf(i,3,k)+a06(4)*uf(i,4,k) &
                    + a06(5)*uf(i,5,k)+a06(6)*uf(i,6,k) &
                    + a06(7)*uf(i,7,k) ) *idy(j)*sinphi(i,j)
          dv = dv + ( a06(1)*vf(i,1,k)+a06(2)*vf(i,2,k) &
                    + a06(3)*vf(i,3,k)+a06(4)*vf(i,4,k) &
                    + a06(5)*vf(i,5,k)+a06(6)*vf(i,6,k) &
                    + a06(7)*vf(i,7,k) ) *idy(j)*sinphi(i,j)
          dw = dw + ( a06(1)*wf(i,1,k)+a06(2)*wf(i,2,k) &
                    + a06(3)*wf(i,3,k)+a06(4)*wf(i,4,k) &
                    + a06(5)*wf(i,5,k)+a06(6)*wf(i,6,k) &
                    + a06(7)*wf(i,7,k) ) *idy(j)*sinphi(i,j)
          dp = dp + ( a06(1)*pf(i,1,k)+a06(2)*pf(i,2,k) &
                    + a06(3)*pf(i,3,k)+a06(4)*pf(i,4,k) &
                    + a06(5)*pf(i,5,k)+a06(6)*pf(i,6,k) &
                    + a06(7)*pf(i,7,k) ) *idy(j)*sinphi(i,j)
          dr = dr + ( a06(1)*rf(i,1,k)+a06(2)*rf(i,2,k) &
                    + a06(3)*rf(i,3,k)+a06(4)*rf(i,4,k) &
                    + a06(5)*rf(i,5,k)+a06(6)*rf(i,6,k) &
                    + a06(7)*rf(i,7,k) ) *idy(j)*sinphi(i,j)

          pt(i,j,k) = vg(i,j,k) * (dp+pf(i,j,k)*ir(i,j))
          ut(i,j,k) = vg(i,j,k) * (du+uf(i,j,k)*ir(i,j))
          vt(i,j,k) = vg(i,j,k) * (dv+vf(i,j,k)*ir(i,j))
          wt(i,j,k) = vg(i,j,k) * (dw+wf(i,j,k)*ir(i,j))
          rt(i,j,k) = vg(i,j,k) * (dr+rf(i,j,k)*ir(i,j))
       enddo
    enddo

    j=2
    do i=ndx,nfx
       do k=1,nz
          du = ( a7(1)*( uf(i+1,j,k) - uf(i-1,j,k) ) + &
                 a7(2)*( uf(i+2,j,k) - uf(i-2,j,k) ) + &
                 a7(3)*( uf(i+3,j,k) - uf(i-3,j,k) ) ) *idx(i)*cosphi(i,j)
          dv = ( a7(1)*( vf(i+1,j,k) - vf(i-1,j,k) ) + &
                 a7(2)*( vf(i+2,j,k) - vf(i-2,j,k) ) + &
                 a7(3)*( vf(i+3,j,k) - vf(i-3,j,k) ) ) *idx(i)*cosphi(i,j)
          dw = ( a7(1)*( wf(i+1,j,k) - wf(i-1,j,k) ) + &
                 a7(2)*( wf(i+2,j,k) - wf(i-2,j,k) ) + &
                 a7(3)*( wf(i+3,j,k) - wf(i-3,j,k) ) ) *idx(i)*cosphi(i,j)
          dp = ( a7(1)*( pf(i+1,j,k) - pf(i-1,j,k) ) + &
                 a7(2)*( pf(i+2,j,k) - pf(i-2,j,k) ) + &
                 a7(3)*( pf(i+3,j,k) - pf(i-3,j,k) ) ) *idx(i)*cosphi(i,j)
          dr = ( a7(1)*( rf(i+1,j,k) - rf(i-1,j,k) ) + &
                 a7(2)*( rf(i+2,j,k) - rf(i-2,j,k) ) + &
                 a7(3)*( rf(i+3,j,k) - rf(i-3,j,k) ) ) *idx(i)*cosphi(i,j)

          du = du + ( a15(1)*uf(i,1,k)+a15(2)*uf(i,2,k) &
                    + a15(3)*uf(i,3,k)+a15(4)*uf(i,4,k) &
                    + a15(5)*uf(i,5,k)+a15(6)*uf(i,6,k) &
                    + a15(7)*uf(i,7,k) ) *idy(j)*sinphi(i,j)
          dv = dv + ( a15(1)*vf(i,1,k)+a15(2)*vf(i,2,k) &
                    + a15(3)*vf(i,3,k)+a15(4)*vf(i,4,k) &
                    + a15(5)*vf(i,5,k)+a15(6)*vf(i,6,k) &
                    + a15(7)*vf(i,7,k) ) *idy(j)*sinphi(i,j)
          dw = dw + ( a15(1)*wf(i,1,k)+a15(2)*wf(i,2,k) &
                    + a15(3)*wf(i,3,k)+a15(4)*wf(i,4,k) &
                    + a15(5)*wf(i,5,k)+a15(6)*wf(i,6,k) &
                    + a15(7)*wf(i,7,k) ) *idy(j)*sinphi(i,j)
          dp = dp + ( a15(1)*pf(i,1,k)+a15(2)*pf(i,2,k) &
                    + a15(3)*pf(i,3,k)+a15(4)*pf(i,4,k) &
                    + a15(5)*pf(i,5,k)+a15(6)*pf(i,6,k) &
                    + a15(7)*pf(i,7,k) ) *idy(j)*sinphi(i,j)
          dr = dr + ( a15(1)*rf(i,1,k)+a15(2)*rf(i,2,k) &
                    + a15(3)*rf(i,3,k)+a15(4)*rf(i,4,k) &
                    + a15(5)*rf(i,5,k)+a15(6)*rf(i,6,k) &
                    + a15(7)*rf(i,7,k) ) *idy(j)*sinphi(i,j)

          pt(i,j,k) = vg(i,j,k) * (dp+pf(i,j,k)*ir(i,j))
          ut(i,j,k) = vg(i,j,k) * (du+uf(i,j,k)*ir(i,j))
          vt(i,j,k) = vg(i,j,k) * (dv+vf(i,j,k)*ir(i,j))
          wt(i,j,k) = vg(i,j,k) * (dw+wf(i,j,k)*ir(i,j))
          rt(i,j,k) = vg(i,j,k) * (dr+rf(i,j,k)*ir(i,j))
       enddo
    enddo

    j=3
    do i=ndx,nfx
       do k=1,nz
          du = ( a7(1)*( uf(i+1,j,k) - uf(i-1,j,k) ) + &
                 a7(2)*( uf(i+2,j,k) - uf(i-2,j,k) ) + &
                 a7(3)*( uf(i+3,j,k) - uf(i-3,j,k) ) ) *idx(i)*cosphi(i,j)
          dv = ( a7(1)*( vf(i+1,j,k) - vf(i-1,j,k) ) + &
                 a7(2)*( vf(i+2,j,k) - vf(i-2,j,k) ) + &
                 a7(3)*( vf(i+3,j,k) - vf(i-3,j,k) ) ) *idx(i)*cosphi(i,j)
          dw = ( a7(1)*( wf(i+1,j,k) - wf(i-1,j,k) ) + &
                 a7(2)*( wf(i+2,j,k) - wf(i-2,j,k) ) + &
                 a7(3)*( wf(i+3,j,k) - wf(i-3,j,k) ) ) *idx(i)*cosphi(i,j)
          dp = ( a7(1)*( pf(i+1,j,k) - pf(i-1,j,k) ) + &
                 a7(2)*( pf(i+2,j,k) - pf(i-2,j,k) ) + &
                 a7(3)*( pf(i+3,j,k) - pf(i-3,j,k) ) ) *idx(i)*cosphi(i,j)
          dr = ( a7(1)*( rf(i+1,j,k) - rf(i-1,j,k) ) + &
                 a7(2)*( rf(i+2,j,k) - rf(i-2,j,k) ) + &
                 a7(3)*( rf(i+3,j,k) - rf(i-3,j,k) ) ) *idx(i)*cosphi(i,j)

          du = du + ( a24(1)*uf(i,1,k)+a24(2)*uf(i,2,k) &
                    + a24(3)*uf(i,3,k)+a24(4)*uf(i,4,k) &
                    + a24(5)*uf(i,5,k)+a24(6)*uf(i,6,k) &
                    + a24(7)*uf(i,7,k) ) *idy(j)*sinphi(i,j)
          dv = dv + ( a24(1)*vf(i,1,k)+a24(2)*vf(i,2,k) &
                    + a24(3)*vf(i,3,k)+a24(4)*vf(i,4,k) &
                    + a24(5)*vf(i,5,k)+a24(6)*vf(i,6,k) &
                    + a24(7)*vf(i,7,k) ) *idy(j)*sinphi(i,j)
          dw = dw + ( a24(1)*wf(i,1,k)+a24(2)*wf(i,2,k) &
                    + a24(3)*wf(i,3,k)+a24(4)*wf(i,4,k) &
                    + a24(5)*wf(i,5,k)+a24(6)*wf(i,6,k) &
                    + a24(7)*wf(i,7,k) ) *idy(j)*sinphi(i,j)
          dp = dp + ( a24(1)*pf(i,1,k)+a24(2)*pf(i,2,k) &
                    + a24(3)*pf(i,3,k)+a24(4)*pf(i,4,k) &
                    + a24(5)*pf(i,5,k)+a24(6)*pf(i,6,k) &
                    + a24(7)*pf(i,7,k) ) *idy(j)*sinphi(i,j)
          dr = dr + ( a24(1)*rf(i,1,k)+a24(2)*rf(i,2,k) &
                    + a24(3)*rf(i,3,k)+a24(4)*rf(i,4,k) &
                    + a24(5)*rf(i,5,k)+a24(6)*rf(i,6,k) &
                    + a24(7)*rf(i,7,k) ) *idy(j)*sinphi(i,j)

          pt(i,j,k) = vg(i,j,k) * (dp+pf(i,j,k)*ir(i,j))
          ut(i,j,k) = vg(i,j,k) * (du+uf(i,j,k)*ir(i,j))
          vt(i,j,k) = vg(i,j,k) * (dv+vf(i,j,k)*ir(i,j))
          wt(i,j,k) = vg(i,j,k) * (dw+wf(i,j,k)*ir(i,j))
          rt(i,j,k) = vg(i,j,k) * (dr+rf(i,j,k)*ir(i,j))
       enddo
    enddo

    do j=4,ngh
       do i=ndx,nfx
          do k=1,nz
             du = ( a7(1)*( uf(i+1,j,k) - uf(i-1,j,k) ) + &
                    a7(2)*( uf(i+2,j,k) - uf(i-2,j,k) ) + &
                    a7(3)*( uf(i+3,j,k) - uf(i-3,j,k) ) ) *idx(i)*cosphi(i,j)
             dv = ( a7(1)*( vf(i+1,j,k) - vf(i-1,j,k) ) + &
                    a7(2)*( vf(i+2,j,k) - vf(i-2,j,k) ) + &
                    a7(3)*( vf(i+3,j,k) - vf(i-3,j,k) ) ) *idx(i)*cosphi(i,j)
             dw = ( a7(1)*( wf(i+1,j,k) - wf(i-1,j,k) ) + &
                    a7(2)*( wf(i+2,j,k) - wf(i-2,j,k) ) + &
                    a7(3)*( wf(i+3,j,k) - wf(i-3,j,k) ) ) *idx(i)*cosphi(i,j)
             dp = ( a7(1)*( pf(i+1,j,k) - pf(i-1,j,k) ) + &
                    a7(2)*( pf(i+2,j,k) - pf(i-2,j,k) ) + &
                    a7(3)*( pf(i+3,j,k) - pf(i-3,j,k) ) ) *idx(i)*cosphi(i,j)
             dr = ( a7(1)*( rf(i+1,j,k) - rf(i-1,j,k) ) + &
                    a7(2)*( rf(i+2,j,k) - rf(i-2,j,k) ) + &
                    a7(3)*( rf(i+3,j,k) - rf(i-3,j,k) ) ) *idx(i)*cosphi(i,j)

             du = du + ( a7(1)*(uf(i,j+1,k)-uf(i,j-1,k)) &
                       + a7(2)*(uf(i,j+2,k)-uf(i,j-2,k)) &
                       + a7(3)*(uf(i,j+3,k)-uf(i,j-3,k)) ) *idy(j)*sinphi(i,j)
             dv = dv + ( a7(1)*(vf(i,j+1,k)-vf(i,j-1,k)) &
                       + a7(2)*(vf(i,j+2,k)-vf(i,j-2,k)) &
                       + a7(3)*(vf(i,j+3,k)-vf(i,j-3,k)) ) *idy(j)*sinphi(i,j)
             dw = dw + ( a7(1)*(wf(i,j+1,k)-wf(i,j-1,k)) &
                       + a7(2)*(wf(i,j+2,k)-wf(i,j-2,k)) &
                       + a7(3)*(wf(i,j+3,k)-wf(i,j-3,k)) ) *idy(j)*sinphi(i,j)
             dp = dp + ( a7(1)*(pf(i,j+1,k)-pf(i,j-1,k)) &
                       + a7(2)*(pf(i,j+2,k)-pf(i,j-2,k)) &
                       + a7(3)*(pf(i,j+3,k)-pf(i,j-3,k)) ) *idy(j)*sinphi(i,j)
             dr = dr + ( a7(1)*(rf(i,j+1,k)-rf(i,j-1,k)) &
                         + a7(2)*(rf(i,j+2,k)-rf(i,j-2,k)) &
                         + a7(3)*(rf(i,j+3,k)-rf(i,j-3,k)) ) *idy(j)*sinphi(i,j)

             pt(i,j,k) = vg(i,j,k) * (dp+pf(i,j,k)*ir(i,j))
             ut(i,j,k) = vg(i,j,k) * (du+uf(i,j,k)*ir(i,j))
             vt(i,j,k) = vg(i,j,k) * (dv+vf(i,j,k)*ir(i,j))
             wt(i,j,k) = vg(i,j,k) * (dw+wf(i,j,k)*ir(i,j))
             rt(i,j,k) = vg(i,j,k) * (dr+rf(i,j,k)*ir(i,j))
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

  end subroutine bc_TD2d_jmin

  !===============================================================================
  module subroutine bc_TD2d_jmax
  !===============================================================================
    !> 2D Tam & Dong's BC: boundary condition at jmax (top) - Cartesian version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l
    real(wp) :: cp,av,inn1,ntm1
    real(wp) :: du,dv,dw,dp,dr
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
             du = ( a7(1)*( uf(i+1,l,k) - uf(i-1,l,k) ) + &
                    a7(2)*( uf(i+2,l,k) - uf(i-2,l,k) ) + &
                    a7(3)*( uf(i+3,l,k) - uf(i-3,l,k) ) ) *idx(i)*cosphi(i,l)
             dv = ( a7(1)*( vf(i+1,l,k) - vf(i-1,l,k) ) + &
                    a7(2)*( vf(i+2,l,k) - vf(i-2,l,k) ) + &
                    a7(3)*( vf(i+3,l,k) - vf(i-3,l,k) ) ) *idx(i)*cosphi(i,l)
             dw = ( a7(1)*( wf(i+1,l,k) - wf(i-1,l,k) ) + &
                    a7(2)*( wf(i+2,l,k) - wf(i-2,l,k) ) + &
                    a7(3)*( wf(i+3,l,k) - wf(i-3,l,k) ) ) *idx(i)*cosphi(i,l)
             dp = ( a7(1)*( pf(i+1,l,k) - pf(i-1,l,k) ) + &
                    a7(2)*( pf(i+2,l,k) - pf(i-2,l,k) ) + &
                    a7(3)*( pf(i+3,l,k) - pf(i-3,l,k) ) ) *idx(i)*cosphi(i,l)
             dr = ( a7(1)*( rf(i+1,l,k) - rf(i-1,l,k) ) + &
                    a7(2)*( rf(i+2,l,k) - rf(i-2,l,k) ) + &
                    a7(3)*( rf(i+3,l,k) - rf(i-3,l,k) ) ) *idx(i)*cosphi(i,l)

             du = du + ( a7(1)*(uf(i,l+1,k)-uf(i,l-1,k)) &
                       + a7(2)*(uf(i,l+2,k)-uf(i,l-2,k)) &
                       + a7(3)*(uf(i,l+3,k)-uf(i,l-3,k)) ) *idy(j)*sinphi(i,l)
             dv = dv + ( a7(1)*(vf(i,l+1,k)-vf(i,l-1,k)) &
                       + a7(2)*(vf(i,l+2,k)-vf(i,l-2,k)) &
                       + a7(3)*(vf(i,l+3,k)-vf(i,l-3,k)) ) *idy(j)*sinphi(i,l)
             dw = dw + ( a7(1)*(wf(i,l+1,k)-wf(i,l-1,k)) &
                       + a7(2)*(wf(i,l+2,k)-wf(i,l-2,k)) &
                       + a7(3)*(wf(i,l+3,k)-wf(i,l-3,k)) ) *idy(j)*sinphi(i,l)
             dp = dp + ( a7(1)*(pf(i,l+1,k)-pf(i,l-1,k)) &
                       + a7(2)*(pf(i,l+2,k)-pf(i,l-2,k)) &
                       + a7(3)*(pf(i,l+3,k)-pf(i,l-3,k)) ) *idy(j)*sinphi(i,l)
             dr = dr + ( a7(1)*(rf(i,l+1,k)-rf(i,l-1,k)) &
                       + a7(2)*(rf(i,l+2,k)-rf(i,l-2,k)) &
                       + a7(3)*(rf(i,l+3,k)-rf(i,l-3,k)) ) *idy(j)*sinphi(i,l)

             pt(i,l,k) = vg(i,l,k) * (dp+pf(i,l,k)*ir(i,l))
             ut(i,l,k) = vg(i,l,k) * (du+uf(i,l,k)*ir(i,l))
             vt(i,l,k) = vg(i,l,k) * (dv+vf(i,l,k)*ir(i,l))
             wt(i,l,k) = vg(i,l,k) * (dw+wf(i,l,k)*ir(i,l))
             rt(i,l,k) = vg(i,l,k) * (dr+rf(i,l,k)*ir(i,l))
          enddo
       enddo
    enddo

    j=ny-2
    l=j-nymngh
    do i=ndx,nfx
       do k=1,nz
          du = ( a7(1)*( uf(i+1,l,k) - uf(i-1,l,k) ) + &
                 a7(2)*( uf(i+2,l,k) - uf(i-2,l,k) ) + &
                 a7(3)*( uf(i+3,l,k) - uf(i-3,l,k) ) ) *idx(i)*cosphi(i,l)
          dv = ( a7(1)*( vf(i+1,l,k) - vf(i-1,l,k) ) + &
                 a7(2)*( vf(i+2,l,k) - vf(i-2,l,k) ) + &
                 a7(3)*( vf(i+3,l,k) - vf(i-3,l,k) ) ) *idx(i)*cosphi(i,l)
          dw = ( a7(1)*( wf(i+1,l,k) - wf(i-1,l,k) ) + &
                 a7(2)*( wf(i+2,l,k) - wf(i-2,l,k) ) + &
                 a7(3)*( wf(i+3,l,k) - wf(i-3,l,k) ) ) *idx(i)*cosphi(i,l)
          dp = ( a7(1)*( pf(i+1,l,k) - pf(i-1,l,k) ) + &
                 a7(2)*( pf(i+2,l,k) - pf(i-2,l,k) ) + &
                 a7(3)*( pf(i+3,l,k) - pf(i-3,l,k) ) ) *idx(i)*cosphi(i,l)
          dr = ( a7(1)*( rf(i+1,l,k) - rf(i-1,l,k) ) + &
                 a7(2)*( rf(i+2,l,k) - rf(i-2,l,k) ) + &
                 a7(3)*( rf(i+3,l,k) - rf(i-3,l,k) ) ) *idx(i)*cosphi(i,l)

          du = du + ( a42(1)*uf(i,l+2,k)+a42(2)*uf(i,l+1,k) &
                    + a42(3)*uf(i,l  ,k)+a42(4)*uf(i,l-1,k) &
                    + a42(5)*uf(i,l-2,k)+a42(6)*uf(i,l-3,k) &
                    + a42(7)*uf(i,l-4,k) ) *idy(j)*sinphi(i,l)
          dv = dv + ( a42(1)*vf(i,l+2,k)+a42(2)*vf(i,l+1,k) &
                    + a42(3)*vf(i,l  ,k)+a42(4)*vf(i,l-1,k) &
                    + a42(5)*vf(i,l-2,k)+a42(6)*vf(i,l-3,k) &
                    + a42(7)*vf(i,l-4,k) ) *idy(j)*sinphi(i,l)
          dw = dw + ( a42(1)*wf(i,l+2,k)+a42(2)*wf(i,l+1,k) &
                    + a42(3)*wf(i,l  ,k)+a42(4)*wf(i,l-1,k) &
                    + a42(5)*wf(i,l-2,k)+a42(6)*wf(i,l-3,k) &
                    + a42(7)*wf(i,l-4,k) ) *idy(j)*sinphi(i,l)
          dp = dp + ( a42(1)*pf(i,l+2,k)+a42(2)*pf(i,l+1,k) &
                    + a42(3)*pf(i,l  ,k)+a42(4)*pf(i,l-1,k) &
                    + a42(5)*pf(i,l-2,k)+a42(6)*pf(i,l-3,k) &
                    + a42(7)*pf(i,l-4,k) ) *idy(j)*sinphi(i,l)
          dr = dr + ( a42(1)*rf(i,l+2,k)+a42(2)*rf(i,l+1,k) &
                    + a42(3)*rf(i,l  ,k)+a42(4)*rf(i,l-1,k) &
                    + a42(5)*rf(i,l-2,k)+a42(6)*rf(i,l-3,k) &
                    + a42(7)*rf(i,l-4,k) ) *idy(j)*sinphi(i,l)

          pt(i,l,k) = vg(i,l,k) * (dp+pf(i,l,k)*ir(i,l))
          ut(i,l,k) = vg(i,l,k) * (du+uf(i,l,k)*ir(i,l))
          vt(i,l,k) = vg(i,l,k) * (dv+vf(i,l,k)*ir(i,l))
          wt(i,l,k) = vg(i,l,k) * (dw+wf(i,l,k)*ir(i,l))
          rt(i,l,k) = vg(i,l,k) * (dr+rf(i,l,k)*ir(i,l))
       enddo
    enddo

    j=ny-1
    l=j-nymngh
    do i=ndx,nfx
       do k=1,nz
          du = ( a7(1)*( uf(i+1,l,k) - uf(i-1,l,k) ) + &
                 a7(2)*( uf(i+2,l,k) - uf(i-2,l,k) ) + &
                 a7(3)*( uf(i+3,l,k) - uf(i-3,l,k) ) ) *idx(i)*cosphi(i,l)
          dv = ( a7(1)*( vf(i+1,l,k) - vf(i-1,l,k) ) + &
                 a7(2)*( vf(i+2,l,k) - vf(i-2,l,k) ) + &
                 a7(3)*( vf(i+3,l,k) - vf(i-3,l,k) ) ) *idx(i)*cosphi(i,l)
          dw = ( a7(1)*( wf(i+1,l,k) - wf(i-1,l,k) ) + &
                 a7(2)*( wf(i+2,l,k) - wf(i-2,l,k) ) + &
                 a7(3)*( wf(i+3,l,k) - wf(i-3,l,k) ) ) *idx(i)*cosphi(i,l)
          dp = ( a7(1)*( pf(i+1,l,k) - pf(i-1,l,k) ) + &
                 a7(2)*( pf(i+2,l,k) - pf(i-2,l,k) ) + &
                 a7(3)*( pf(i+3,l,k) - pf(i-3,l,k) ) ) *idx(i)*cosphi(i,l)
          dr = ( a7(1)*( rf(i+1,l,k) - rf(i-1,l,k) ) + &
                 a7(2)*( rf(i+2,l,k) - rf(i-2,l,k) ) + &
                 a7(3)*( rf(i+3,l,k) - rf(i-3,l,k) ) ) *idx(i)*cosphi(i,l)

          du = du + ( a51(1)*uf(i,l+1,k)+a51(2)*uf(i,l  ,k) &
                    + a51(3)*uf(i,l-1,k)+a51(4)*uf(i,l-2,k) &
                    + a51(5)*uf(i,l-3,k)+a51(6)*uf(i,l-4,k) &
                    + a51(7)*uf(i,l-5,k) ) *idy(j)*sinphi(i,l)
          dv = dv + ( a51(1)*vf(i,l+1,k)+a51(2)*vf(i,l  ,k) &
                    + a51(3)*vf(i,l-1,k)+a51(4)*vf(i,l-2,k) &
                    + a51(5)*vf(i,l-3,k)+a51(6)*vf(i,l-4,k) &
                    + a51(7)*vf(i,l-5,k) ) *idy(j)*sinphi(i,l)
          dw = dw + ( a51(1)*wf(i,l+1,k)+a51(2)*wf(i,l  ,k) &
                    + a51(3)*wf(i,l-1,k)+a51(4)*wf(i,l-2,k) &
                    + a51(5)*wf(i,l-3,k)+a51(6)*wf(i,l-4,k) &
                    + a51(7)*wf(i,l-5,k) ) *idy(j)*sinphi(i,l)
          dp = dp + ( a51(1)*pf(i,l+1,k)+a51(2)*pf(i,l  ,k) &
                    + a51(3)*pf(i,l-1,k)+a51(4)*pf(i,l-2,k) &
                    + a51(5)*pf(i,l-3,k)+a51(6)*pf(i,l-4,k) &
                    + a51(7)*pf(i,l-5,k) ) *idy(j)*sinphi(i,l)
          dr = dr + ( a51(1)*rf(i,l+1,k)+a51(2)*rf(i,l  ,k) &
                    + a51(3)*rf(i,l-1,k)+a51(4)*rf(i,l-2,k) &
                    + a51(5)*rf(i,l-3,k)+a51(6)*rf(i,l-4,k) &
                    + a51(7)*rf(i,l-5,k) ) *idy(j)*sinphi(i,l)

          pt(i,l,k) = vg(i,l,k) * (dp+pf(i,l,k)*ir(i,l))
          ut(i,l,k) = vg(i,l,k) * (du+uf(i,l,k)*ir(i,l))
          vt(i,l,k) = vg(i,l,k) * (dv+vf(i,l,k)*ir(i,l))
          wt(i,l,k) = vg(i,l,k) * (dw+wf(i,l,k)*ir(i,l))
          rt(i,l,k) = vg(i,l,k) * (dr+rf(i,l,k)*ir(i,l))
       enddo
    enddo

    j=ny
    l=j-nymngh
    do i=ndx,nfx
       do k=1,nz
          du = ( a7(1)*( uf(i+1,l,k) - uf(i-1,l,k) ) + &
                 a7(2)*( uf(i+2,l,k) - uf(i-2,l,k) ) + &
                 a7(3)*( uf(i+3,l,k) - uf(i-3,l,k) ) ) *idx(i)*cosphi(i,l)
          dv = ( a7(1)*( vf(i+1,l,k) - vf(i-1,l,k) ) + &
                 a7(2)*( vf(i+2,l,k) - vf(i-2,l,k) ) + &
                 a7(3)*( vf(i+3,l,k) - vf(i-3,l,k) ) ) *idx(i)*cosphi(i,l)
          dw = ( a7(1)*( wf(i+1,l,k) - wf(i-1,l,k) ) + &
                 a7(2)*( wf(i+2,l,k) - wf(i-2,l,k) ) + &
                 a7(3)*( wf(i+3,l,k) - wf(i-3,l,k) ) ) *idx(i)*cosphi(i,l)
          dp = ( a7(1)*( pf(i+1,l,k) - pf(i-1,l,k) ) + &
                 a7(2)*( pf(i+2,l,k) - pf(i-2,l,k) ) + &
                 a7(3)*( pf(i+3,l,k) - pf(i-3,l,k) ) ) *idx(i)*cosphi(i,l)
          dr = ( a7(1)*( rf(i+1,l,k) - rf(i-1,l,k) ) + &
                 a7(2)*( rf(i+2,l,k) - rf(i-2,l,k) ) + &
                 a7(3)*( rf(i+3,l,k) - rf(i-3,l,k) ) ) *idx(i)*cosphi(i,l)

          du = du + ( a60(1)*uf(i,l  ,k)+a60(2)*uf(i,l-1,k) &
                    + a60(3)*uf(i,l-2,k)+a60(4)*uf(i,l-3,k) &
                    + a60(5)*uf(i,l-4,k)+a60(6)*uf(i,l-5,k) &
                    + a60(7)*uf(i,l-6,k) ) *idy(j)*sinphi(i,l)
          dv = dv + ( a60(1)*vf(i,l  ,k)+a60(2)*vf(i,l-1,k) &
                    + a60(3)*vf(i,l-2,k)+a60(4)*vf(i,l-3,k) &
                    + a60(5)*vf(i,l-4,k)+a60(6)*vf(i,l-5,k) &
                    + a60(7)*vf(i,l-6,k) ) *idy(j)*sinphi(i,l)
          dw = dw + ( a60(1)*wf(i,l  ,k)+a60(2)*wf(i,l-1,k) &
                    + a60(3)*wf(i,l-2,k)+a60(4)*wf(i,l-3,k) &
                    + a60(5)*wf(i,l-4,k)+a60(6)*wf(i,l-5,k) &
                    + a60(7)*wf(i,l-6,k) ) *idy(j)*sinphi(i,l)
          dp = dp + ( a60(1)*pf(i,l  ,k)+a60(2)*pf(i,l-1,k) &
                    + a60(3)*pf(i,l-2,k)+a60(4)*pf(i,l-3,k) &
                    + a60(5)*pf(i,l-4,k)+a60(6)*pf(i,l-5,k) &
                    + a60(7)*pf(i,l-6,k) ) *idy(j)*sinphi(i,l)
          dr = dr + ( a60(1)*rf(i,l  ,k)+a60(2)*rf(i,l-1,k) &
                    + a60(3)*rf(i,l-2,k)+a60(4)*rf(i,l-3,k) &
                    + a60(5)*rf(i,l-4,k)+a60(6)*rf(i,l-5,k) &
                    + a60(7)*rf(i,l-6,k) ) *idy(j)*sinphi(i,l)

          pt(i,l,k) = vg(i,l,k) * (dp+pf(i,l,k)*ir(i,l))
          ut(i,l,k) = vg(i,l,k) * (du+uf(i,l,k)*ir(i,l))
          vt(i,l,k) = vg(i,l,k) * (dv+vf(i,l,k)*ir(i,l))
          wt(i,l,k) = vg(i,l,k) * (dw+wf(i,l,k)*ir(i,l))
          rt(i,l,k) = vg(i,l,k) * (dr+rf(i,l,k)*ir(i,l))
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

  end subroutine bc_TD2d_jmax

end submodule smod_TamDong2d_faces
