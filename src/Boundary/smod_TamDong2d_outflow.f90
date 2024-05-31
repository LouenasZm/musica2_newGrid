!===============================================================================
submodule (mod_TamDong2d) smod_TamDong2d_outflow
!===============================================================================
  !> author: XG
  !> date: February 2020 - modif January 2022
  !> Outflow boundary conditions of Tam and Dong
  !> - 2D (periodic) Cartesian version -
  !> /!\ ~~~> restricted to imax (exit = right condition)
!=============================================================================== 

contains

  !===============================================================================
  module subroutine bc_TD2d_outflow_imax
  !===============================================================================
    !> 2D Tam & Dong's outflow BC: boundary condition at imax (right) - Cartesian version -
  !===============================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l
    real(wp) :: cp,av,inn1,ntm1
    real(wp), dimension(-2:ngh,ny1:ny2,nz1:nz2) :: rf,uf,vf,wf,pf
    real(wp), dimension(ngh,ny,nz) :: rt,ut,vt,wt,pt,vg,c2_
    real(wp) :: dux,dvx,dwx,dpx,drx
    real(wp) :: duy,dvy,dwy,dpy,dry
    real(wp) :: duz,dvz,dwz,dpz,drz
    !-------------------------------------------------------------------------

    ! Pointers for polar coordinates
    ! ==============================
    ir=>BC_face(1,2)%ir
    cosphi=>BC_face(1,2)%cosphi
    sinphi=>BC_face(1,2)%sinphi
    
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
    do i=nxmnghp1,nx-3
       l=i-nxmngh
       do j=ndy,nfy
          do k=1,nz
             dux = ( a7(1)* (uf(l+1,j,k)-uf(l-1,j,k)) &
                   + a7(2)* (uf(l+2,j,k)-uf(l-2,j,k)) &
                   + a7(3)* (uf(l+3,j,k)-uf(l-3,j,k)) ) *idx(i)
             dvx = ( a7(1)* (vf(l+1,j,k)-vf(l-1,j,k)) &
                   + a7(2)* (vf(l+2,j,k)-vf(l-2,j,k)) &
                   + a7(3)* (vf(l+3,j,k)-vf(l-3,j,k)) ) *idx(i)
             dwx = ( a7(1)* (wf(l+1,j,k)-wf(l-1,j,k)) &
                   + a7(2)* (wf(l+2,j,k)-wf(l-2,j,k)) &
                   + a7(3)* (wf(l+3,j,k)-wf(l-3,j,k)) ) *idx(i)
             dpx = ( a7(1)* (pf(l+1,j,k)-pf(l-1,j,k)) &
                   + a7(2)* (pf(l+2,j,k)-pf(l-2,j,k)) &
                   + a7(3)* (pf(l+3,j,k)-pf(l-3,j,k)) ) *idx(i)
             drx = ( a7(1)* (rf(l+1,j,k)-rf(l-1,j,k)) &
                   + a7(2)* (rf(l+2,j,k)-rf(l-2,j,k)) &
                   + a7(3)* (rf(l+3,j,k)-rf(l-3,j,k)) ) *idx(i)

             duy = ( a7(1)* (uf(l,j+1,k)-uf(l,j-1,k)) &
                   + a7(2)* (uf(l,j+2,k)-uf(l,j-2,k)) + &
                   + a7(3)* (uf(l,j+3,k)-uf(l,j-3,k)) ) *idy(j)
             dvy = ( a7(1)* (vf(l,j+1,k)-vf(l,j-1,k)) + &
                   + a7(2)* (vf(l,j+2,k)-vf(l,j-2,k)) + &
                   + a7(3)* (vf(l,j+3,k)-vf(l,j-3,k)) ) *idy(j)
             dwy = ( a7(1)* (wf(l,j+1,k)-wf(l,j-1,k)) + &
                   + a7(2)* (wf(l,j+2,k)-wf(l,j-2,k)) + &
                   + a7(3)* (wf(l,j+3,k)-wf(l,j-3,k)) ) *idy(j)
             dpy = ( a7(1)* (pf(l,j+1,k)-pf(l,j-1,k)) + &
                   + a7(2)* (pf(l,j+2,k)-pf(l,j-2,k)) + &
                   + a7(3)* (pf(l,j+3,k)-pf(l,j-3,k)) ) *idy(j)
             dry = ( a7(1)* (rf(l,j+1,k)-rf(l,j-1,k)) + &
                   + a7(2)* (rf(l,j+2,k)-rf(l,j-2,k)) + &
                   + a7(3)* (rf(l,j+3,k)-rf(l,j-3,k)) ) *idy(j)

             if (.not.is_2D) then
                duz = ( a7(1)* (uf(l,j,k+1)-uf(l,j,k-1)) &
                      + a7(2)* (uf(l,j,k+2)-uf(l,j,k-2)) + &
                      + a7(3)* (uf(l,j,k+3)-uf(l,j,k-3)) ) *idz(k)
                dvz = ( a7(1)* (vf(l,j,k+1)-vf(l,j,k-1)) &
                      + a7(2)* (vf(l,j,k+2)-vf(l,j,k-2)) + &
                      + a7(3)* (vf(l,j,k+3)-vf(l,j,k-3)) ) *idz(k)
                dwz = ( a7(1)* (wf(l,j,k+1)-wf(l,j,k-1)) &
                      + a7(2)* (wf(l,j,k+2)-wf(l,j,k-2)) + &
                      + a7(3)* (wf(l,j,k+3)-wf(l,j,k-3)) ) *idz(k)
                dpz = ( a7(1)* (pf(l,j,k+1)-pf(l,j,k-1)) &
                      + a7(2)* (pf(l,j,k+2)-pf(l,j,k-2)) + &
                      + a7(3)* (pf(l,j,k+3)-pf(l,j,k-3)) ) *idz(k)
                drz = ( a7(1)* (rf(l,j,k+1)-rf(l,j,k-1)) &
                      + a7(2)* (rf(l,j,k+2)-rf(l,j,k-2)) + &
                      + a7(3)* (rf(l,j,k+3)-rf(l,j,k-3)) ) *idz(k)
             endif

             pt(l,j,k)=vg(l,j,k)*(dpx*cosphi(l,j)+dpy*sinphi(l,j)+pf(l,j,k)*ir(l,j))

             ! ut=(u0.Grad)uf+dpx/rho0=u0*dux+v0*duy+w0*duz + dpx/rho0
             ut(l,j,k)=BC_face(1,2)%U0(l,j,k,2)*dux &
                      +BC_face(1,2)%U0(l,j,k,3)*duy &
                      +BC_face(1,2)%U0(l,j,k,4)*duz &
                      + dpx/BC_face(1,2)%U0(l,j,k,1)
             ! vt=(u0.Grad)vf+dpy/rho0=u0*dvx+v0*dvy+w0*dvz + dpy/rho0
             vt(l,j,k)=BC_face(1,2)%U0(l,j,k,2)*dvx &
                      +BC_face(1,2)%U0(l,j,k,3)*dvy &
                      +BC_face(1,2)%U0(l,j,k,4)*dvz &
                      + dpy/BC_face(1,2)%U0(l,j,k,1)
             ! wt=(u0.Grad)wf+dpz/rho0=u0*dwx+v0*dwy+w0*dwz + dpz/rho0
             wt(l,j,k)=BC_face(1,2)%U0(l,j,k,2)*dwx &
                      +BC_face(1,2)%U0(l,j,k,3)*dwy &
                      +BC_face(1,2)%U0(l,j,k,4)*dwz &
                      + dpz/BC_face(1,2)%U0(l,j,k,1)

             ! rhot=(u0.Grad)rhof-(pt+(u0.Grad)pf)/c0^2
             rt(l,j,k)=BC_face(1,2)%U0(l,j,k,2)*drx &
                      +BC_face(1,2)%U0(l,j,k,3)*dry &
                      +BC_face(1,2)%U0(l,j,k,4)*drz &
                      -(pt(l,j,k) &
                      +BC_face(1,2)%U0(l,j,k,2)*dpx &
                      +BC_face(1,2)%U0(l,j,k,3)*dpy &
                      +BC_face(1,2)%U0(l,j,k,4)*dpz)/BC_face(1,2)%U0(l,j,k,6)
          enddo
       enddo
    enddo

    i=nx-2
    l=i-nxmngh
    do j=ndy,nfy
       do k=1,nz
          dux = ( a42(1)*uf(l+2,j,k)+a42(2)*uf(l+1,j,k) &
                + a42(3)*uf(l  ,j,k)+a42(4)*uf(l-1,j,k) &
                + a42(5)*uf(l-2,j,k)+a42(6)*uf(l-3,j,k) &
                + a42(7)*uf(l-4,j,k) ) *idx(i)
          dvx = ( a42(1)*vf(l+2,j,k)+a42(2)*vf(l+1,j,k) &
                + a42(3)*vf(l  ,j,k)+a42(4)*vf(l-1,j,k) &
                + a42(5)*vf(l-2,j,k)+a42(6)*vf(l-3,j,k) &
                + a42(7)*vf(l-4,j,k) ) *idx(i)
          dwx = ( a42(1)*wf(l+2,j,k)+a42(2)*wf(l+1,j,k) &
                + a42(3)*wf(l  ,j,k)+a42(4)*wf(l-1,j,k) &
                + a42(5)*wf(l-2,j,k)+a42(6)*wf(l-3,j,k) &
                + a42(7)*wf(l-4,j,k) ) *idx(i)
          dpx = ( a42(1)*pf(l+2,j,k)+a42(2)*pf(l+1,j,k) &
                + a42(3)*pf(l  ,j,k)+a42(4)*pf(l-1,j,k) &
                + a42(5)*pf(l-2,j,k)+a42(6)*pf(l-3,j,k) &
                + a42(7)*pf(l-4,j,k) ) *idx(i)
          drx = ( a42(1)*rf(l+2,j,k)+a42(2)*rf(l+1,j,k) &
                + a42(3)*rf(l  ,j,k)+a42(4)*rf(l-1,j,k) &
                + a42(5)*rf(l-2,j,k)+a42(6)*rf(l-3,j,k) &
                + a42(7)*rf(l-4,j,k) ) *idx(i)
          
          duy = ( a7(1)* (uf(l,j+1,k)-uf(l,j-1,k)) &
                + a7(2)* (uf(l,j+2,k)-uf(l,j-2,k)) + &
                + a7(3)* (uf(l,j+3,k)-uf(l,j-3,k)) ) *idy(j)
          dvy = ( a7(1)* (vf(l,j+1,k)-vf(l,j-1,k)) + &
                + a7(2)* (vf(l,j+2,k)-vf(l,j-2,k)) + &
                + a7(3)* (vf(l,j+3,k)-vf(l,j-3,k)) ) *idy(j)
          dwy = ( a7(1)* (wf(l,j+1,k)-wf(l,j-1,k)) + &
                + a7(2)* (wf(l,j+2,k)-wf(l,j-2,k)) + &
                + a7(3)* (wf(l,j+3,k)-wf(l,j-3,k)) ) *idy(j)
          dpy = ( a7(1)* (pf(l,j+1,k)-pf(l,j-1,k)) + &
                + a7(2)* (pf(l,j+2,k)-pf(l,j-2,k)) + &
                + a7(3)* (pf(l,j+3,k)-pf(l,j-3,k)) ) *idy(j)
          dry = ( a7(1)* (rf(l,j+1,k)-rf(l,j-1,k)) + &
                + a7(2)* (rf(l,j+2,k)-rf(l,j-2,k)) + &
                + a7(3)* (rf(l,j+3,k)-rf(l,j-3,k)) ) *idy(j)

          if (.not.is_2D) then
             duz = ( a7(1)* (uf(l,j,k+1)-uf(l,j,k-1)) &
                   + a7(2)* (uf(l,j,k+2)-uf(l,j,k-2)) + &
                   + a7(3)* (uf(l,j,k+3)-uf(l,j,k-3)) ) *idz(k)
             dvz = ( a7(1)* (vf(l,j,k+1)-vf(l,j,k-1)) &
                   + a7(2)* (vf(l,j,k+2)-vf(l,j,k-2)) + &
                   + a7(3)* (vf(l,j,k+3)-vf(l,j,k-3)) ) *idz(k)
             dwz = ( a7(1)* (wf(l,j,k+1)-wf(l,j,k-1)) &
                   + a7(2)* (wf(l,j,k+2)-wf(l,j,k-2)) + &
                   + a7(3)* (wf(l,j,k+3)-wf(l,j,k-3)) ) *idz(k)
             dpz = ( a7(1)* (pf(l,j,k+1)-pf(l,j,k-1)) &
                   + a7(2)* (pf(l,j,k+2)-pf(l,j,k-2)) + &
                   + a7(3)* (pf(l,j,k+3)-pf(l,j,k-3)) ) *idz(k)
             drz = ( a7(1)* (rf(l,j,k+1)-rf(l,j,k-1)) &
                   + a7(2)* (rf(l,j,k+2)-rf(l,j,k-2)) + &
                   + a7(3)* (rf(l,j,k+3)-rf(l,j,k-3)) ) *idz(k)
          endif

          pt(l,j,k)=vg(l,j,k)*(dpx*cosphi(l,j)+dpy*sinphi(l,j)+pf(l,j,k)*ir(l,j))

          ! ut=(u0.Grad)uf+dpx/rho0=u0*dux+v0*duy+w0*duz + dpx/rho0
          ut(l,j,k)=BC_face(1,2)%U0(l,j,k,2)*dux &
                   +BC_face(1,2)%U0(l,j,k,3)*duy &
                   +BC_face(1,2)%U0(l,j,k,4)*duz &
                   + dpx/BC_face(1,2)%U0(l,j,k,1)
          ! vt=(u0.Grad)vf+dpy/rho0=u0*dvx+v0*dvy+w0*dvz + dpy/rho0
          vt(l,j,k)=BC_face(1,2)%U0(l,j,k,2)*dvx &
                   +BC_face(1,2)%U0(l,j,k,3)*dvy &
                   +BC_face(1,2)%U0(l,j,k,4)*dvz &
                   + dpy/BC_face(1,2)%U0(l,j,k,1)
          ! wt=(u0.Grad)wf+dpz/rho0=u0*dwx+v0*dwy+w0*dwz + dpz/rho0
          wt(l,j,k)=BC_face(1,2)%U0(l,j,k,2)*dwx &
                   +BC_face(1,2)%U0(l,j,k,3)*dwy &
                   +BC_face(1,2)%U0(l,j,k,4)*dwz &
                   + dpz/BC_face(1,2)%U0(l,j,k,1)

          ! rhot=(u0.Grad)rhof-(pt+(u0.Grad)pf)/c0^2
          rt(l,j,k)=BC_face(1,2)%U0(l,j,k,2)*drx &
                   +BC_face(1,2)%U0(l,j,k,3)*dry &
                   +BC_face(1,2)%U0(l,j,k,4)*drz &
                   -(pt(l,j,k) &
                   +BC_face(1,2)%U0(l,j,k,2)*dpx &
                   +BC_face(1,2)%U0(l,j,k,3)*dpy &
                   +BC_face(1,2)%U0(l,j,k,4)*dpz)/BC_face(1,2)%U0(l,j,k,6)
       enddo
    enddo

    i=nx-1
    l=i-nxmngh
    do j=ndy,nfy
       do k=1,nz
          dux = ( a51(1)*uf(l+1,j,k)+a51(2)*uf(l  ,j,k) &
                + a51(3)*uf(l-1,j,k)+a51(4)*uf(l-2,j,k) &
                + a51(5)*uf(l-3,j,k)+a51(6)*uf(l-4,j,k) &
                + a51(7)*uf(l-5,j,k) ) *idx(i)
          dvx = ( a51(1)*vf(l+1,j,k)+a51(2)*vf(l  ,j,k) &
                + a51(3)*vf(l-1,j,k)+a51(4)*vf(l-2,j,k) &
                + a51(5)*vf(l-3,j,k)+a51(6)*vf(l-4,j,k) &
                + a51(7)*vf(l-5,j,k) ) *idx(i)
          dwx = ( a51(1)*wf(l+1,j,k)+a51(2)*wf(l  ,j,k) &
                + a51(3)*wf(l-1,j,k)+a51(4)*wf(l-2,j,k) &
                + a51(5)*wf(l-3,j,k)+a51(6)*wf(l-4,j,k) &
                + a51(7)*wf(l-5,j,k) ) *idx(i)
          dpx = ( a51(1)*pf(l+1,j,k)+a51(2)*pf(l  ,j,k) &
                + a51(3)*pf(l-1,j,k)+a51(4)*pf(l-2,j,k) &
                + a51(5)*pf(l-3,j,k)+a51(6)*pf(l-4,j,k) &
                + a51(7)*pf(l-5,j,k) ) *idx(i)
          drx = ( a51(1)*rf(l+1,j,k)+a51(2)*rf(l  ,j,k) &
                + a51(3)*rf(l-1,j,k)+a51(4)*rf(l-2,j,k) &
                + a51(5)*rf(l-3,j,k)+a51(6)*rf(l-4,j,k) &
                + a51(7)*rf(l-5,j,k) ) *idx(i)

          duy = ( a7(1)* (uf(l,j+1,k)-uf(l,j-1,k)) &
                + a7(2)* (uf(l,j+2,k)-uf(l,j-2,k)) + &
                + a7(3)* (uf(l,j+3,k)-uf(l,j-3,k)) ) *idy(j)
          dvy = ( a7(1)* (vf(l,j+1,k)-vf(l,j-1,k)) + &
                + a7(2)* (vf(l,j+2,k)-vf(l,j-2,k)) + &
                + a7(3)* (vf(l,j+3,k)-vf(l,j-3,k)) ) *idy(j)
          dwy = ( a7(1)* (wf(l,j+1,k)-wf(l,j-1,k)) + &
                + a7(2)* (wf(l,j+2,k)-wf(l,j-2,k)) + &
                + a7(3)* (wf(l,j+3,k)-wf(l,j-3,k)) ) *idy(j)
          dpy = ( a7(1)* (pf(l,j+1,k)-pf(l,j-1,k)) + &
                + a7(2)* (pf(l,j+2,k)-pf(l,j-2,k)) + &
                + a7(3)* (pf(l,j+3,k)-pf(l,j-3,k)) ) *idy(j)
          dry = ( a7(1)* (rf(l,j+1,k)-rf(l,j-1,k)) + &
                + a7(2)* (rf(l,j+2,k)-rf(l,j-2,k)) + &
                + a7(3)* (rf(l,j+3,k)-rf(l,j-3,k)) ) *idy(j)

          if (.not.is_2D) then
             duz = ( a7(1)* (uf(l,j,k+1)-uf(l,j,k-1)) &
                   + a7(2)* (uf(l,j,k+2)-uf(l,j,k-2)) + &
                   + a7(3)* (uf(l,j,k+3)-uf(l,j,k-3)) ) *idz(k)
             dvz = ( a7(1)* (vf(l,j,k+1)-vf(l,j,k-1)) &
                   + a7(2)* (vf(l,j,k+2)-vf(l,j,k-2)) + &
                   + a7(3)* (vf(l,j,k+3)-vf(l,j,k-3)) ) *idz(k)
             dwz = ( a7(1)* (wf(l,j,k+1)-wf(l,j,k-1)) &
                   + a7(2)* (wf(l,j,k+2)-wf(l,j,k-2)) + &
                   + a7(3)* (wf(l,j,k+3)-wf(l,j,k-3)) ) *idz(k)
             dpz = ( a7(1)* (pf(l,j,k+1)-pf(l,j,k-1)) &
                   + a7(2)* (pf(l,j,k+2)-pf(l,j,k-2)) + &
                   + a7(3)* (pf(l,j,k+3)-pf(l,j,k-3)) ) *idz(k)
             drz = ( a7(1)* (rf(l,j,k+1)-rf(l,j,k-1)) &
                   + a7(2)* (rf(l,j,k+2)-rf(l,j,k-2)) + &
                   + a7(3)* (rf(l,j,k+3)-rf(l,j,k-3)) ) *idz(k)
          endif

          pt(l,j,k)=vg(l,j,k)*(dpx*cosphi(l,j)+dpy*sinphi(l,j)+pf(l,j,k)*ir(l,j))

          ! ut=(u0.Grad)uf+dpx/rho0=u0*dux+v0*duy+w0*duz + dpx/rho0
          ut(l,j,k)=BC_face(1,2)%U0(l,j,k,2)*dux &
                   +BC_face(1,2)%U0(l,j,k,3)*duy &
                   +BC_face(1,2)%U0(l,j,k,4)*duz &
                   + dpx/BC_face(1,2)%U0(l,j,k,1)
          ! vt=(u0.Grad)vf+dpy/rho0=u0*dvx+v0*dvy+w0*dvz + dpy/rho0
          vt(l,j,k)=BC_face(1,2)%U0(l,j,k,2)*dvx &
                   +BC_face(1,2)%U0(l,j,k,3)*dvy &
                   +BC_face(1,2)%U0(l,j,k,4)*dvz &
                   + dpy/BC_face(1,2)%U0(l,j,k,1)
          ! wt=(u0.Grad)wf+dpz/rho0=u0*dwx+v0*dwy+w0*dwz + dpz/rho0
          wt(l,j,k)=BC_face(1,2)%U0(l,j,k,2)*dwx &
                   +BC_face(1,2)%U0(l,j,k,3)*dwy &
                   +BC_face(1,2)%U0(l,j,k,4)*dwz &
                   + dpz/BC_face(1,2)%U0(l,j,k,1)

          ! rhot=(u0.Grad)rhof-(pt+(u0.Grad)pf)/c0^2
          rt(l,j,k)=BC_face(1,2)%U0(l,j,k,2)*drx &
                   +BC_face(1,2)%U0(l,j,k,3)*dry &
                   +BC_face(1,2)%U0(l,j,k,4)*drz &
                   -(pt(l,j,k) &
                   +BC_face(1,2)%U0(l,j,k,2)*dpx &
                   +BC_face(1,2)%U0(l,j,k,3)*dpy &
                   +BC_face(1,2)%U0(l,j,k,4)*dpz)/BC_face(1,2)%U0(l,j,k,6)
       enddo
    enddo

    i=nx
    l=i-nxmngh
    do j=ndy,nfy
       do k=1,nz
          dux = ( a60(1)*uf(l  ,j,k)+a60(2)*uf(l-1,j,k) &
                + a60(3)*uf(l-2,j,k)+a60(4)*uf(l-3,j,k) &
                + a60(5)*uf(l-4,j,k)+a60(6)*uf(l-5,j,k) &
                + a60(7)*uf(l-6,j,k) ) *idx(i)
          dvx = ( a60(1)*vf(l  ,j,k)+a60(2)*vf(l-1,j,k) &
                + a60(3)*vf(l-2,j,k)+a60(4)*vf(l-3,j,k) &
                + a60(5)*vf(l-4,j,k)+a60(6)*vf(l-5,j,k) &
                + a60(7)*vf(l-6,j,k) ) *idx(i)
          dwx = ( a60(1)*wf(l  ,j,k)+a60(2)*wf(l-1,j,k) &
                + a60(3)*wf(l-2,j,k)+a60(4)*wf(l-3,j,k) &
                + a60(5)*wf(l-4,j,k)+a60(6)*wf(l-5,j,k) &
                + a60(7)*wf(l-6,j,k) ) *idx(i)
          dpx = ( a60(1)*pf(l  ,j,k)+a60(2)*pf(l-1,j,k) &
                + a60(3)*pf(l-2,j,k)+a60(4)*pf(l-3,j,k) &
                + a60(5)*pf(l-4,j,k)+a60(6)*pf(l-5,j,k) &
                + a60(7)*pf(l-6,j,k) ) *idx(i)
          drx = ( a60(1)*rf(l  ,j,k)+a60(2)*rf(l-1,j,k) &
                + a60(3)*rf(l-2,j,k)+a60(4)*rf(l-3,j,k) &
                + a60(5)*rf(l-4,j,k)+a60(6)*rf(l-5,j,k) &
                + a60(7)*rf(l-6,j,k) ) *idx(i)
          
          duy = ( a7(1)* (uf(l,j+1,k)-uf(l,j-1,k)) &
                + a7(2)* (uf(l,j+2,k)-uf(l,j-2,k)) + &
                + a7(3)* (uf(l,j+3,k)-uf(l,j-3,k)) ) *idy(j)
          dvy = ( a7(1)* (vf(l,j+1,k)-vf(l,j-1,k)) + &
                + a7(2)* (vf(l,j+2,k)-vf(l,j-2,k)) + &
                + a7(3)* (vf(l,j+3,k)-vf(l,j-3,k)) ) *idy(j)
          dwy = ( a7(1)* (wf(l,j+1,k)-wf(l,j-1,k)) + &
                + a7(2)* (wf(l,j+2,k)-wf(l,j-2,k)) + &
                + a7(3)* (wf(l,j+3,k)-wf(l,j-3,k)) ) *idy(j)
          dpy = ( a7(1)* (pf(l,j+1,k)-pf(l,j-1,k)) + &
                + a7(2)* (pf(l,j+2,k)-pf(l,j-2,k)) + &
                + a7(3)* (pf(l,j+3,k)-pf(l,j-3,k)) ) *idy(j)
          dry = ( a7(1)* (rf(l,j+1,k)-rf(l,j-1,k)) + &
                + a7(2)* (rf(l,j+2,k)-rf(l,j-2,k)) + &
                + a7(3)* (rf(l,j+3,k)-rf(l,j-3,k)) ) *idy(j)

          if (.not.is_2D) then
             duz = ( a7(1)* (uf(l,j,k+1)-uf(l,j,k-1)) &
                   + a7(2)* (uf(l,j,k+2)-uf(l,j,k-2)) + &
                   + a7(3)* (uf(l,j,k+3)-uf(l,j,k-3)) ) *idz(k)
             dvz = ( a7(1)* (vf(l,j,k+1)-vf(l,j,k-1)) &
                   + a7(2)* (vf(l,j,k+2)-vf(l,j,k-2)) + &
                   + a7(3)* (vf(l,j,k+3)-vf(l,j,k-3)) ) *idz(k)
             dwz = ( a7(1)* (wf(l,j,k+1)-wf(l,j,k-1)) &
                   + a7(2)* (wf(l,j,k+2)-wf(l,j,k-2)) + &
                   + a7(3)* (wf(l,j,k+3)-wf(l,j,k-3)) ) *idz(k)
             dpz = ( a7(1)* (pf(l,j,k+1)-pf(l,j,k-1)) &
                   + a7(2)* (pf(l,j,k+2)-pf(l,j,k-2)) + &
                   + a7(3)* (pf(l,j,k+3)-pf(l,j,k-3)) ) *idz(k)
             drz = ( a7(1)* (rf(l,j,k+1)-rf(l,j,k-1)) &
                   + a7(2)* (rf(l,j,k+2)-rf(l,j,k-2)) + &
                   + a7(3)* (rf(l,j,k+3)-rf(l,j,k-3)) ) *idz(k)
          endif

          pt(l,j,k)=vg(l,j,k)*(dpx*cosphi(l,j)+dpy*sinphi(l,j)+pf(l,j,k)*ir(l,j))

          ! ut=(u0.Grad)uf+dpx/rho0=u0*dux+v0*duy+w0*duz + dpx/rho0
          ut(l,j,k)=BC_face(1,2)%U0(l,j,k,2)*dux &
                   +BC_face(1,2)%U0(l,j,k,3)*duy &
                   +BC_face(1,2)%U0(l,j,k,4)*duz &
                   + dpx/BC_face(1,2)%U0(l,j,k,1)
          ! vt=(u0.Grad)vf+dpy/rho0=u0*dvx+v0*dvy+w0*dvz + dpy/rho0
          vt(l,j,k)=BC_face(1,2)%U0(l,j,k,2)*dvx &
                   +BC_face(1,2)%U0(l,j,k,3)*dvy &
                   +BC_face(1,2)%U0(l,j,k,4)*dvz &
                   + dpy/BC_face(1,2)%U0(l,j,k,1)
          ! wt=(u0.Grad)wf+dpz/rho0=u0*dwx+v0*dwy+w0*dwz + dpz/rho0
          wt(l,j,k)=BC_face(1,2)%U0(l,j,k,2)*dwx &
                   +BC_face(1,2)%U0(l,j,k,3)*dwy &
                   +BC_face(1,2)%U0(l,j,k,4)*dwz &
                   + dpz/BC_face(1,2)%U0(l,j,k,1)

          ! rhot=(u0.Grad)rhof-(pt+(u0.Grad)pf)/c0^2
          rt(l,j,k)=BC_face(1,2)%U0(l,j,k,2)*drx &
                   +BC_face(1,2)%U0(l,j,k,3)*dry &
                   +BC_face(1,2)%U0(l,j,k,4)*drz &
                   -(pt(l,j,k) &
                   +BC_face(1,2)%U0(l,j,k,2)*dpx &
                   +BC_face(1,2)%U0(l,j,k,3)*dpy &
                   +BC_face(1,2)%U0(l,j,k,4)*dpz)/BC_face(1,2)%U0(l,j,k,6)
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

  end subroutine bc_TD2d_outflow_imax

end submodule smod_TamDong2d_outflow
