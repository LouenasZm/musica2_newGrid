!===============================================================================
submodule (mod_TamDong2d_1pt) smod_TamDong2d_1pt
!===============================================================================
  !> author: AB
  !> date: April 2022 (created from 5 points version)
  !> Non-reflecting boundary conditions of Tam and Dong
  !> - 2D / 2.5D (periodic) Cartesian version - version on 1 point
  !> /!\ NOT REGULAR - kept for dev (tests inflow conditions)
!===============================================================================

contains

  !==========================================================================================
  module subroutine bc_TD2d_1pt_imin
  !==========================================================================================
    !> 2D Tam & Dong's BC on 1 point: boundary condition at imin (left) - Cartesian version -
  !==========================================================================================
    use mod_eigenmode
    use mod_RFM
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: cp,av,inn1,ntm1
    real(wp) :: du,dv,dw,dp,dr,ideltax,ideltay
    real(wp), dimension(1:5,ny1:ny2,nz) :: rf,uf,vf,wf,pf
    real(wp), dimension(ny,nz) :: rt,ut,vt,wt,pt,vg,c2_
    !-------------------------------------------------------------------------
    ! eigenmode or RFM disturbances
    ! real(wp) :: pt_in,ut_in,vt_in,wt_in,rt_in
    ! real(wp) :: dp_in,du_in,dv_in,dw_in,dr_in
    real(wp), dimension(ny,nz) :: pt_in,ut_in,vt_in,wt_in,rt_in
    !-------------------------------------------------------------------------
    ! param vortex -> TEMPORARY, to be changed
    real(wp) :: exm,exm2,xa,ya,x0,ar,trk
    ! to directly impose u,v,w TO BE CHANGED
    ! real(wp) :: u_in,v_in,w_in
    real(wp) :: dutdx_in,dvtdx_in,dptdx_in,dutdy_in,dvtdy_in,dptdy_in
    real(wp) :: dutdt_in,dvtdt_in,dptdt_in,ampl
    real(wp) :: dp_in,du_in,dv_in
    !-------------------------------------------------------------------------
    ! added for is_mean_ref
    integer :: n_moy
    real(wp) :: eps
    !-------------------------------------------------------------------------

    ! Index of left boundary
    ! ======================
    i=1

    ! Pointers for polar coordinates
    ! ==============================
    ir=>BC_face(1,1)%ir
    cosphi=>BC_face(1,1)%cosphi
    sinphi=>BC_face(1,1)%sinphi

    ! Wall condition applied (different from Tam&Dong 5 points)
    ! ======================
    ! Wall BC at jmin
    ! ---------------
    if (is_bc_wall(2,1)) then
        Krho(i,1,ndz_e:nfz_e)=0.0_wp
       Krhou(i,1,ndz_e:nfz_e)=0.0_wp
       Krhov(i,1,ndz_e:nfz_e)=0.0_wp
       Krhow(i,1,ndz_e:nfz_e)=0.0_wp
       Krhoe(i,1,ndz_e:nfz_e)=0.0_wp
    endif
    ! Wall BC at jmax
    ! ---------------
    if (is_bc_wall(2,2)) then
        Krho(i,ny,ndz_e:nfz_e)=0.0_wp
       Krhou(i,ny,ndz_e:nfz_e)=0.0_wp
       Krhov(i,ny,ndz_e:nfz_e)=0.0_wp
       Krhow(i,ny,ndz_e:nfz_e)=0.0_wp
       Krhoe(i,ny,ndz_e:nfz_e)=0.0_wp
    endif
    ! Wall BC at kmin
    ! ---------------
    if (is_bc_wall(3,1)) then
        Krho(i,ndy_e:nfy_e,1)=0.0_wp
       Krhou(i,ndy_e:nfy_e,1)=0.0_wp
       Krhov(i,ndy_e:nfy_e,1)=0.0_wp
       Krhow(i,ndy_e:nfy_e,1)=0.0_wp
       Krhoe(i,ndy_e:nfy_e,1)=0.0_wp
    endif
    ! Wall BC at kmax
    ! ---------------
    if (is_bc_wall(3,2)) then
        Krho(i,ndy_e:nfy_e,nz)=0.0_wp
       Krhou(i,ndy_e:nfy_e,nz)=0.0_wp
       Krhov(i,ndy_e:nfy_e,nz)=0.0_wp
       Krhow(i,ndy_e:nfy_e,nz)=0.0_wp
       Krhoe(i,ndy_e:nfy_e,nz)=0.0_wp
    endif

    ! Sound speed
    ! ===========
    do k=1,nz
       do j=1,ny
          c2_(j,k)=c2calc_tro(Tmp(i,j,k),rho_n(i,j,k))
       enddo
    enddo

    ! Compute online time-averaged primitive variables
    ! ================================================
    if (irk==nrk) then
       inn1 = 1.0_wp/dble(ntotal)
       ntm1=dble(ntotal-1)
       if (BC_face(1,1)%is_mean_ref) then
          eps=0.1_wp
          n_moy=100
          inn1=inn1*eps
          ! time-averaged primitive variables
          do k=1,nz
             do j=ny1,ny2
                BC_face(1,1)%U0(1:5,j,k,1)=(ntm1*BC_face(1,1)%U0(1:5,j,k,1)+rho_n(1:5,j,k))*inn1 &
                                              + (1.0_wp-eps)*BC_face(1,1)%Uref(j,k,1)
                BC_face(1,1)%U0(1:5,j,k,2)=(ntm1*BC_face(1,1)%U0(1:5,j,k,2)+uu(1:5,j,k))*inn1    &
                                              + (1.0_wp-eps)*BC_face(1,1)%Uref(j,k,2)
                BC_face(1,1)%U0(1:5,j,k,3)=(ntm1*BC_face(1,1)%U0(1:5,j,k,3)+vv(1:5,j,k))*inn1    &
                                              + (1.0_wp-eps)*BC_face(1,1)%Uref(j,k,3)
                BC_face(1,1)%U0(1:5,j,k,4)=(ntm1*BC_face(1,1)%U0(1:5,j,k,4)+ww(1:5,j,k))*inn1    &
                                              + (1.0_wp-eps)*BC_face(1,1)%Uref(j,k,4)
                BC_face(1,1)%U0(1:5,j,k,5)=(ntm1*BC_face(1,1)%U0(1:5,j,k,5)+prs(1:5,j,k))*inn1   &
                                              + (1.0_wp-eps)*BC_face(1,1)%Uref(j,k,5)
             enddo
          enddo
          ! time-averaged sound speed squared
          do k=1,nz
             do j=1,ny
                BC_face(1,1)%U0(1,j,k,6)=(ntm1*BC_face(1,1)%U0(1,j,k,6)+c2_(j,k))*inn1 &
                                               + (1.0_wp-eps)*BC_face(1,1)%Uref(j,k,6)
             enddo
          enddo
       else
          ! time-averaged primitive variables
          do k=1,nz
             BC_face(1,1)%U0(1:5,:,k,1)=(ntm1*BC_face(1,1)%U0(1:5,:,k,1)+rho_n(1:5,:,k))*inn1
             BC_face(1,1)%U0(1:5,:,k,2)=(ntm1*BC_face(1,1)%U0(1:5,:,k,2)+uu(1:5,:,k))*inn1
             BC_face(1,1)%U0(1:5,:,k,3)=(ntm1*BC_face(1,1)%U0(1:5,:,k,3)+vv(1:5,:,k))*inn1
             BC_face(1,1)%U0(1:5,:,k,4)=(ntm1*BC_face(1,1)%U0(1:5,:,k,4)+ww(1:5,:,k))*inn1
             BC_face(1,1)%U0(1:5,:,k,5)=(ntm1*BC_face(1,1)%U0(1:5,:,k,5)+prs(1:5,:,k))*inn1
          enddo
          ! time-averaged sound speed squared
          do k=1,nz
             BC_face(1,1)%U0(1,1:ny,k,6)=(ntm1*BC_face(1,1)%U0(1,1:ny,k,6)+c2_(1:ny,k))*inn1
          enddo
       endif
    endif

    ! Compute fluctuations
    ! ====================
    do i=1,5
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
    i=1
    do j=ndy_td1,nfy_td1
       do k=1,nz
          vg(j,k)= BC_face(1,1)%U0(i,j,k,2)*cosphi(i,j)+BC_face(1,1)%U0(i,j,k,3)*sinphi(i,j) &
               + sqrt(BC_face(1,1)%U0(i,j,k,6)-BC_face(1,1)%U0(i,j,k,4)**2 &
               -(BC_face(1,1)%U0(i,j,k,2)*sinphi(i,j)-BC_face(1,1)%U0(i,j,k,3)*cosphi(i,j))**2)
       enddo
    enddo

    ! Compute vg*[cos(phi)*dq/dx+sin(phi)*dq/dy+q/r]
    ! ==============================================
    i=1
    ideltax = 1.0_wp/(x(i+1)-x(i))
    do j=ndy_td1,nfy_td1
       do k=1,nz
          ! (Tam & Webb DRP schemes)
          ! ------------------------
          ! du = ( a06(1)*uf(1,j,k)+a06(2)*uf(2,j,k) &
          !      + a06(3)*uf(3,j,k)+a06(4)*uf(4,j,k) &
          !      + a06(5)*uf(5,j,k)+a06(6)*uf(6,j,k) &
          !      + a06(7)*uf(7,j,k) ) *idx(i)*cosphi(i,j)
          ! dv = ( a06(1)*vf(1,j,k)+a06(2)*vf(2,j,k) &
          !      + a06(3)*vf(3,j,k)+a06(4)*vf(4,j,k) &
          !      + a06(5)*vf(5,j,k)+a06(6)*vf(6,j,k) &
          !      + a06(7)*vf(7,j,k) ) *idx(i)*cosphi(i,j)
          ! dw = ( a06(1)*wf(1,j,k)+a06(2)*wf(2,j,k) &
          !      + a06(3)*wf(3,j,k)+a06(4)*wf(4,j,k) &
          !      + a06(5)*wf(5,j,k)+a06(6)*wf(6,j,k) &
          !      + a06(7)*wf(7,j,k) ) *idx(i)*cosphi(i,j)
          ! dp = ( a06(1)*pf(1,j,k)+a06(2)*pf(2,j,k) &
          !      + a06(3)*pf(3,j,k)+a06(4)*pf(4,j,k) &
          !      + a06(5)*pf(5,j,k)+a06(6)*pf(6,j,k) &
          !      + a06(7)*pf(7,j,k) ) *idx(i)*cosphi(i,j)
          ! dr = ( a06(1)*rf(1,j,k)+a06(2)*rf(2,j,k) &
          !      + a06(3)*rf(3,j,k)+a06(4)*rf(4,j,k) &
          !      + a06(5)*rf(5,j,k)+a06(6)*rf(6,j,k) &
          !      + a06(7)*rf(7,j,k) ) *idx(i)*cosphi(i,j)

          ! du = du + (a7(1)*( uf(i,j+1,k) - uf(i,j-1,k) ) + &
          !            a7(2)*( uf(i,j+2,k) - uf(i,j-2,k) ) + &
          !            a7(3)*( uf(i,j+3,k) - uf(i,j-3,k) ) ) *idy(j)*sinphi(i,j)
          ! dv = dv + (a7(1)*( vf(i,j+1,k) - vf(i,j-1,k) ) + &
          !            a7(2)*( vf(i,j+2,k) - vf(i,j-2,k) ) + &
          !            a7(3)*( vf(i,j+3,k) - vf(i,j-3,k) ) ) *idy(j)*sinphi(i,j)
          ! dw = dw + (a7(1)*( wf(i,j+1,k) - wf(i,j-1,k) ) + &
          !            a7(2)*( wf(i,j+2,k) - wf(i,j-2,k) ) + &
          !            a7(3)*( wf(i,j+3,k) - wf(i,j-3,k) ) ) *idy(j)*sinphi(i,j)
          ! dp = dp + (a7(1)*( pf(i,j+1,k) - pf(i,j-1,k) ) + &
          !            a7(2)*( pf(i,j+2,k) - pf(i,j-2,k) ) + &
          !            a7(3)*( pf(i,j+3,k) - pf(i,j-3,k) ) ) *idy(j)*sinphi(i,j)
          ! dr = dr + (a7(1)*( rf(i,j+1,k) - rf(i,j-1,k) ) + &
          !            a7(2)*( rf(i,j+2,k) - rf(i,j-2,k) ) + &
          !            a7(3)*( rf(i,j+3,k) - rf(i,j-3,k) ) ) *idy(j)*sinphi(i,j)


          ! Standard schemes on 5-point stencil
          ! -----------------------------------
          ! order 4 (decentered 0-4)
          ! du = ( a04(1)*uf(i  ,j,k)+a04(2)*uf(i+1,j,k) &
          !      + a04(3)*uf(i+2,j,k)+a04(4)*uf(i+3,j,k) &
          !      + a04(5)*uf(i+4,j,k) )*idx(i)*cosphi(i,j)
          ! dv = ( a04(1)*vf(i  ,j,k)+a04(2)*vf(i+1,j,k) &
          !      + a04(3)*vf(i+2,j,k)+a04(4)*vf(i+3,j,k) &
          !      + a04(5)*vf(i+4,j,k) )*idx(i)*cosphi(i,j)
          ! dw = ( a04(1)*wf(i  ,j,k)+a04(2)*wf(i+1,j,k) &
          !      + a04(3)*wf(i+2,j,k)+a04(4)*wf(i+3,j,k) &
          !      + a04(5)*wf(i+4,j,k) )*idx(i)*cosphi(i,j)
          ! dp = ( a04(1)*pf(i  ,j,k)+a04(2)*pf(i+1,j,k) &
          !      + a04(3)*pf(i+2,j,k)+a04(4)*pf(i+3,j,k) &
          !      + a04(5)*pf(i+4,j,k) )*idx(i)*cosphi(i,j)
          ! dr = ( a04(1)*rf(i  ,j,k)+a04(2)*rf(i+1,j,k) &
          !      + a04(3)*rf(i+2,j,k)+a04(4)*rf(i+3,j,k) &
          !      + a04(5)*rf(i+4,j,k) )*idx(i)*cosphi(i,j)

          ! Centered order 4 in y-direction
          ! du = du + ( a5(1)*( uf(i,j+1,k)-uf(i,j-1,k)) &
          !           + a5(2)*( uf(i,j+2,k)-uf(i,j-2,k)) )*idy(j)*sinphi(i,j)
          ! dv = dv + ( a5(1)*( vf(i,j+1,k)-vf(i,j-1,k)) &
          !           + a5(2)*( vf(i,j+2,k)-vf(i,j-2,k)) )*idy(j)*sinphi(i,j)
          ! dw = dw + ( a5(1)*( wf(i,j+1,k)-wf(i,j-1,k)) &
          !           + a5(2)*( wf(i,j+2,k)-wf(i,j-2,k)) )*idy(j)*sinphi(i,j)
          ! dp = dp + ( a5(1)*( pf(i,j+1,k)-pf(i,j-1,k)) &
          !           + a5(2)*( pf(i,j+2,k)-pf(i,j-2,k)) )*idy(j)*sinphi(i,j)
          ! dr = dr + ( a5(1)*( rf(i,j+1,k)-rf(i,j-1,k)) &
          !           + a5(2)*( rf(i,j+2,k)-rf(i,j-2,k)) )*idy(j)*sinphi(i,j)

          ! Non-centered order 1 in x-direction (forward differences)
          du = (uf(i+1,j,k)-uf(i,j,k))*ideltax*cosphi(i,j)
          dv = (vf(i+1,j,k)-vf(i,j,k))*ideltax*cosphi(i,j)
          dw = (wf(i+1,j,k)-wf(i,j,k))*ideltax*cosphi(i,j)
          dp = (pf(i+1,j,k)-pf(i,j,k))*ideltax*cosphi(i,j)
          dr = (rf(i+1,j,k)-rf(i,j,k))*ideltax*cosphi(i,j)

          ! Centered order 4 in y-direction
          du = du + ( a5(1)*( uf(i,j+1,k)-uf(i,j-1,k)) &
                    + a5(2)*( uf(i,j+2,k)-uf(i,j-2,k)) )*idy(j)*sinphi(i,j)
          dv = dv + ( a5(1)*( vf(i,j+1,k)-vf(i,j-1,k)) &
                    + a5(2)*( vf(i,j+2,k)-vf(i,j-2,k)) )*idy(j)*sinphi(i,j)
          dw = dw + ( a5(1)*( wf(i,j+1,k)-wf(i,j-1,k)) &
                    + a5(2)*( wf(i,j+2,k)-wf(i,j-2,k)) )*idy(j)*sinphi(i,j)
          dp = dp + ( a5(1)*( pf(i,j+1,k)-pf(i,j-1,k)) &
                    + a5(2)*( pf(i,j+2,k)-pf(i,j-2,k)) )*idy(j)*sinphi(i,j)
          dr = dr + ( a5(1)*( rf(i,j+1,k)-rf(i,j-1,k)) &
                    + a5(2)*( rf(i,j+2,k)-rf(i,j-2,k)) )*idy(j)*sinphi(i,j)

          pt(j,k) = vg(j,k) * (dp + pf(i,j,k)*ir(i,j))
          ut(j,k) = vg(j,k) * (du + uf(i,j,k)*ir(i,j))
          vt(j,k) = vg(j,k) * (dv + vf(i,j,k)*ir(i,j))
          wt(j,k) = vg(j,k) * (dw + wf(i,j,k)*ir(i,j))
          rt(j,k) = vg(j,k) * (dr + rf(i,j,k)*ir(i,j))
       enddo
    enddo

!     if (is_eigenmode) then
!        do i=1,ngh
!           do j=ndy_td1,nfy_td1
!              do k=1,nz
!                 call eig_disturb2_imin(i,j,k,pt_in,ut_in,vt_in,wt_in,rt_in,dp_in,du_in,dv_in,dw_in,dr_in)
!                 pt(i,j,k) = pt(i,j,k) - vg(j,k)*dp_in - pt_in
!                 ut(i,j,k) = ut(i,j,k) - vg(j,k)*du_in - ut_in
!                 vt(i,j,k) = vt(i,j,k) - vg(j,k)*dv_in - vt_in
!                 wt(i,j,k) = wt(i,j,k) - vg(j,k)*dw_in - wt_in
!                 rt(i,j,k) = rt(i,j,k) - vg(j,k)*dr_in - rt_in
! !!$                pt(i,j,k) = - pt_in
! !!$                ut(i,j,k) = - ut_in
! !!$                vt(i,j,k) = - vt_in
! !!$                wt(i,j,k) = - wt_in
! !!$                rt(i,j,k) = - rt_in
!              enddo
!           enddo
!        enddo
!     endif

    ! Random Fourier Modes: not injected in the edge point
    if (is_RFM) then
       call disturb_inlet_RFM_TamDong1pt_imin(vg,ut_in,vt_in,wt_in)
       i=1
       do j=ndy_td1,nfy_td1
          do k=1,nz
             ! ut_in = vg(j,k)*du_in + ut_in
             ut(j,k) = ut(j,k) - ut_in(j,k)
             vt(j,k) = vt(j,k) - vt_in(j,k)
             wt(j,k) = wt(j,k) - wt_in(j,k)
          enddo
       enddo
    endif

!!$    if ((.not.is_vortex).and.(type_vortex.eq.-1)) then
!!$       ! Gaussian half-width
!!$       ar = log(2.0_wp)/(25.0_wp*deltay**2)
!!$       ! position of vortex center (+ virtual initial position)
!!$       trk = time + ck(irk)*deltat
!!$       x0 = -15.0_wp*deltay + U_ref*trk
!!$       ampl = 1.0_wp
!!$       i=1
!!$       do j=ndy_td1,nfy_td1
!!$          do k=1,nz
!!$             xa=(x(i)-x0)
!!$             ya=(y(j)-yg(ngy/2))
!!$             exm=exp(-ar*(xa**2+ya**2))
!!$             exm2=exp(-2.0_wp*ar*(xa**2+ya**2))
!!$
!!$             dutdt_in = ampl*U_ref*2*ar*(xa*ya/deltay)*exm
!!$             dutdx_in = -ampl*2*ar*(xa*ya/deltay)*exm
!!$             dutdy_in = (ampl/deltay)*exm*(1-2*ar*ya**2)
!!$             du_in = ampl*(ya/deltay)*exm*ir(i,j) + dutdx_in*cosphi(i,j) + dutdy_in*sinphi(i,j)
!!$
!!$             dvtdt_in = (ampl*U_ref/deltay)*(1 - 2*ar*xa**2)*exm
!!$             dvtdx_in = (ampl/deltay)*(2*ar*xa**2 - 1)*exm
!!$             dvtdy_in = ampl*(xa*ya/deltay)*2*ar*exm
!!$             dv_in = -ampl*(xa/deltay)*exm*ir(i,j) + dvtdx_in*cosphi(i,j) + dvtdy_in*sinphi(i,j)
!!$
!!$             dptdt_in = -rho_ref*(ampl/deltay)**2 * U_ref*xa*exm2
!!$             dptdx_in = rho_ref*(ampl/deltay)**2 * xa*exm2
!!$             dptdy_in = rho_ref*(ampl/deltay)**2 * ya*exm2
!!$             dp_in = - (rho(i,j,k)*ampl**2)/(4*ar*deltay**2)*exm2*ir(i,j) + dptdx_in*cosphi(i,j) + dptdy_in*sinphi(i,j)
!!$
!!$             ! Case 1: with spatial and time derivatives
!!$             ! ut(j,k) = ut(j,k) - vg(j,k)*du_in - dutdt_in
!!$             ! vt(j,k) = vt(j,k) - vg(j,k)*dv_in - dvtdt_in
!!$             ! pt(j,k) = pt(j,k) - vg(j,k)*dp_in - dptdt_in
!!$             ! Case 2: with only time derivatives
!!$             ! ut(j,k) = ut(j,k) - dutdt_in
!!$             ! vt(j,k) = vt(j,k) - dvtdt_in
!!$             ! pt(j,k) = pt(j,k) - dptdt_in
!!$             ! Case 3: directly imposed
!!$             ut(j,k) = - dutdt_in
!!$             vt(j,k) = - dvtdt_in
!!$             wt(j,k) = - dptdt_in
!!$             pt(j,k) = 0.0_wp
!!$             rt(j,k) = 0.0_wp
!!$          enddo
!!$       enddo
!!$    endif

    ! Boundary condition at imin-jmin (edge 1,1,1 /left-bottom)
    ! =========================================================
    if ((coord(2)==0).and.((BC_face(2,1)%sort==-11).or.(BC_face(2,1)%sort==0))) then
       i=1

       ! Pointers for polar coordinates
       ! ==============================
       ir=>BC_edge(1,1,1)%ir
       cosphi=>BC_edge(1,1,1)%cosphi
       sinphi=>BC_edge(1,1,1)%sinphi

       ! Compute fluctuations
       ! ====================
       ! already done before

       ! Compute group velocity vg
       ! =========================
       ! vg=u0*cos(phi)+v0*sin(phi) + sqrt(c0^2-w0^2-(u0*sin(phi)+v0*cos(phi))^2)
       do j=1,2
          do k=1,nz
             vg(j,k)= BC_face(1,1)%U0(i,j,k,2)*cosphi(i,j)+BC_face(1,1)%U0(i,j,k,3)*sinphi(i,j) &
                  + sqrt(BC_face(1,1)%U0(i,j,k,6)-BC_face(1,1)%U0(i,j,k,4)**2 &
                  -(BC_face(1,1)%U0(i,j,k,2)*sinphi(i,j)-BC_face(1,1)%U0(i,j,k,3)*cosphi(i,j))**2)
          enddo
       enddo

       ! Non-centered derivatives
       ! ========================
       if (BC_face(2,1)%sort==-11) then
          i=1;j=1
          ideltax = 1.0_wp/(x(i+1)-x(i))
          ideltay = 1.0_wp/(y(j+1)-y(j))
          do k=1,nz
             ! ! (Tam & Webb DRP schemes)
             ! ! Non-centered derivatives in x-direction *cos(phi)
             ! ! =======================================
             ! du = ( a06(1)*uf(1,j,k)+a06(2)*uf(2,j,k) &
             !      + a06(3)*uf(3,j,k)+a06(4)*uf(4,j,k) &
             !      + a06(5)*uf(5,j,k)+a06(6)*uf(6,j,k) &
             !      + a06(7)*uf(7,j,k) ) *idx(i)*cosphi(i,j)
             ! dv = ( a06(1)*vf(1,j,k)+a06(2)*vf(2,j,k) &
             !      + a06(3)*vf(3,j,k)+a06(4)*vf(4,j,k) &
             !      + a06(5)*vf(5,j,k)+a06(6)*vf(6,j,k) &
             !      + a06(7)*vf(7,j,k) ) *idx(i)*cosphi(i,j)
             ! dw = ( a06(1)*wf(1,j,k)+a06(2)*wf(2,j,k) &
             !      + a06(3)*wf(3,j,k)+a06(4)*wf(4,j,k) &
             !      + a06(5)*wf(5,j,k)+a06(6)*wf(6,j,k) &
             !      + a06(7)*wf(7,j,k) ) *idx(i)*cosphi(i,j)
             ! dp = ( a06(1)*pf(1,j,k)+a06(2)*pf(2,j,k) &
             !      + a06(3)*pf(3,j,k)+a06(4)*pf(4,j,k) &
             !      + a06(5)*pf(5,j,k)+a06(6)*pf(6,j,k) &
             !      + a06(7)*pf(7,j,k) ) *idx(i)*cosphi(i,j)
             ! dr = ( a06(1)*rf(1,j,k)+a06(2)*rf(2,j,k) &
             !      + a06(3)*rf(3,j,k)+a06(4)*rf(4,j,k) &
             !      + a06(5)*rf(5,j,k)+a06(6)*rf(6,j,k) &
             !      + a06(7)*rf(7,j,k) ) *idx(i)*cosphi(i,j)

             ! ! Non-centered derivatives in y-direction *sin(phi)
             ! ! =======================================
             ! du = du + ( a06(1)*uf(i,1,k)+a06(2)*uf(i,2,k) &
             !           + a06(3)*uf(i,3,k)+a06(4)*uf(i,4,k) &
             !           + a06(5)*uf(i,5,k)+a06(6)*uf(i,6,k) &
             !           + a06(7)*uf(i,7,k) ) *idy(j)*sinphi(i,j)
             ! dv = dv + ( a06(1)*vf(i,1,k)+a06(2)*vf(i,2,k) &
             !           + a06(3)*vf(i,3,k)+a06(4)*vf(i,4,k) &
             !           + a06(5)*vf(i,5,k)+a06(6)*vf(i,6,k) &
             !           + a06(7)*vf(i,7,k) ) *idy(j)*sinphi(i,j)
             ! dw = dw + ( a06(1)*wf(i,1,k)+a06(2)*wf(i,2,k) &
             !           + a06(3)*wf(i,3,k)+a06(4)*wf(i,4,k) &
             !           + a06(5)*wf(i,5,k)+a06(6)*wf(i,6,k) &
             !           + a06(7)*wf(i,7,k) ) *idy(j)*sinphi(i,j)
             ! dp = dp + ( a06(1)*pf(i,1,k)+a06(2)*pf(i,2,k) &
             !           + a06(3)*pf(i,3,k)+a06(4)*pf(i,4,k) &
             !           + a06(5)*pf(i,5,k)+a06(6)*pf(i,6,k) &
             !           + a06(7)*pf(i,7,k) ) *idy(j)*sinphi(i,j)
             ! dr = dr + ( a06(1)*rf(i,1,k)+a06(2)*rf(i,2,k) &
             !           + a06(3)*rf(i,3,k)+a06(4)*rf(i,4,k) &
             !           + a06(5)*rf(i,5,k)+a06(6)*rf(i,6,k) &
             !           + a06(7)*rf(i,7,k) ) *idy(j)*sinphi(i,j)


             ! ! Non-centered derivatives in x-direction *cos(phi)
             ! ! =======================================
             ! du = ( a04(1)*uf(i  ,j,k)+a04(2)*uf(i+1,j,k) &
             !      + a04(3)*uf(i+2,j,k)+a04(4)*uf(i+3,j,k) &
             !      + a04(5)*uf(i+4,j,k) )*idx(i)*cosphi(i,j)
             ! dv = ( a04(1)*vf(i  ,j,k)+a04(2)*vf(i+1,j,k) &
             !      + a04(3)*vf(i+2,j,k)+a04(4)*vf(i+3,j,k) &
             !      + a04(5)*vf(i+4,j,k) )*idx(i)*cosphi(i,j)
             ! dw = ( a04(1)*wf(i  ,j,k)+a04(2)*wf(i+1,j,k) &
             !      + a04(3)*wf(i+2,j,k)+a04(4)*wf(i+3,j,k) &
             !      + a04(5)*wf(i+4,j,k) )*idx(i)*cosphi(i,j)
             ! dp = ( a04(1)*pf(i  ,j,k)+a04(2)*pf(i+1,j,k) &
             !      + a04(3)*pf(i+2,j,k)+a04(4)*pf(i+3,j,k) &
             !      + a04(5)*pf(i+4,j,k) )*idx(i)*cosphi(i,j)
             ! dr = ( a04(1)*rf(i  ,j,k)+a04(2)*rf(i+1,j,k) &
             !      + a04(3)*rf(i+2,j,k)+a04(4)*rf(i+3,j,k) &
             !      + a04(5)*rf(i+4,j,k) )*idx(i)*cosphi(i,j)

             ! ! Non-centered derivatives in y-direction *sin(phi)
             ! ! =======================================
             ! du = du + ( a04(1)*uf(i,j  ,k)+a04(2)*uf(i,j+1,k) &
             !           + a04(3)*uf(i,j+2,k)+a04(4)*uf(i,j+3,k) &
             !           + a04(5)*uf(i,j+4,k) )*idy(j)*sinphi(i,j)
             ! dv = dv + ( a04(1)*vf(i,j  ,k)+a04(2)*vf(i,j+1,k) &
             !           + a04(3)*vf(i,j+2,k)+a04(4)*vf(i,j+3,k) &
             !           + a04(5)*vf(i,j+4,k) )*idy(j)*sinphi(i,j)
             ! dw = dw + ( a04(1)*wf(i,j  ,k)+a04(2)*wf(i,j+1,k) &
             !           + a04(3)*wf(i,j+2,k)+a04(4)*wf(i,j+3,k) &
             !           + a04(5)*wf(i,j+4,k) )*idy(j)*sinphi(i,j)
             ! dp = dp + ( a04(1)*pf(i,j  ,k)+a04(2)*pf(i,j+1,k) &
             !           + a04(3)*pf(i,j+2,k)+a04(4)*pf(i,j+3,k) &
             !           + a04(5)*pf(i,j+4,k) )*idy(j)*sinphi(i,j)
             ! dr = dr + ( a04(1)*rf(i,j  ,k)+a04(2)*rf(i,j+1,k) &
             !           + a04(3)*rf(i,j+2,k)+a04(4)*rf(i,j+3,k) &
             !           + a04(5)*rf(i,j+4,k) )*idy(j)*sinphi(i,j)

             ! Non-centered order 1 in x-direction (forward differences)
             du = (uf(i+1,j,k)-uf(i,j,k))*ideltax*cosphi(i,j)
             dv = (vf(i+1,j,k)-vf(i,j,k))*ideltax*cosphi(i,j)
             dw = (wf(i+1,j,k)-wf(i,j,k))*ideltax*cosphi(i,j)
             dp = (pf(i+1,j,k)-pf(i,j,k))*ideltax*cosphi(i,j)
             dr = (rf(i+1,j,k)-rf(i,j,k))*ideltax*cosphi(i,j)

             ! Non-centered order 1 in y-direction (forward differences)
             du = du + (uf(i,j+1,k) - uf(i,j,k))*ideltay*sinphi(i,j)
             dv = dv + (vf(i,j+1,k) - vf(i,j,k))*ideltay*sinphi(i,j)
             dw = dw + (wf(i,j+1,k) - wf(i,j,k))*ideltay*sinphi(i,j)
             dp = dp + (pf(i,j+1,k) - pf(i,j,k))*ideltay*sinphi(i,j)
             dr = dr + (rf(i,j+1,k) - rf(i,j,k))*ideltay*sinphi(i,j)

             pt(j,k) = vg(j,k)*(dp + pf(i,j,k)*ir(i,j))
             ut(j,k) = vg(j,k)*(du + uf(i,j,k)*ir(i,j))
             vt(j,k) = vg(j,k)*(dv + vf(i,j,k)*ir(i,j))
             wt(j,k) = vg(j,k)*(dw + wf(i,j,k)*ir(i,j))
             rt(j,k) = vg(j,k)*(dr + rf(i,j,k)*ir(i,j))


          enddo
       endif

       i=1;j=2
       ideltax = 1.0_wp/(x(i+1)-x(i))
       do k=1,nz
          ! ! Non-centered derivatives in x-direction *cos(phi)
          ! ! =======================================
          ! du = ( a04(1)*uf(i  ,j,k)+a04(2)*uf(i+1,j,k) &
          !      + a04(3)*uf(i+2,j,k)+a04(4)*uf(i+3,j,k) &
          !      + a04(5)*uf(i+4,j,k) )*idx(i)*cosphi(i,j)
          ! dv = ( a04(1)*vf(i  ,j,k)+a04(2)*vf(i+1,j,k) &
          !      + a04(3)*vf(i+2,j,k)+a04(4)*vf(i+3,j,k) &
          !      + a04(5)*vf(i+4,j,k) )*idx(i)*cosphi(i,j)
          ! dw = ( a04(1)*wf(i  ,j,k)+a04(2)*wf(i+1,j,k) &
          !      + a04(3)*wf(i+2,j,k)+a04(4)*wf(i+3,j,k) &
          !      + a04(5)*wf(i+4,j,k) )*idx(i)*cosphi(i,j)
          ! dp = ( a04(1)*pf(i  ,j,k)+a04(2)*pf(i+1,j,k) &
          !      + a04(3)*pf(i+2,j,k)+a04(4)*pf(i+3,j,k) &
          !      + a04(5)*pf(i+4,j,k) )*idx(i)*cosphi(i,j)
          ! dr = ( a04(1)*rf(i  ,j,k)+a04(2)*rf(i+1,j,k) &
          !      + a04(3)*rf(i+2,j,k)+a04(4)*rf(i+3,j,k) &
          !      + a04(5)*rf(i+4,j,k) )*idx(i)*cosphi(i,j)

          ! ! Non-centered derivatives in y-direction *sin(phi)
          ! ! =======================================
          ! du = du + ( a13(1)*uf(i,j-1,k)+a13(2)*uf(i,j,k) &
          !           + a13(3)*uf(i,j+1,k)+a13(4)*uf(i,j+2,k) &
          !           + a13(5)*uf(i,j+3,k) )*idy(j)*sinphi(i,j)
          ! dv = dv + ( a13(1)*vf(i,j-1,k)+a13(2)*vf(i,j,k) &
          !           + a13(3)*vf(i,j+1,k)+a13(4)*vf(i,j+2,k) &
          !           + a13(5)*vf(i,j+3,k) )*idy(j)*sinphi(i,j)
          ! dw = dw + ( a13(1)*wf(i,j-1,k)+a13(2)*wf(i,j,k) &
          !           + a13(3)*wf(i,j+1,k)+a13(4)*wf(i,j+2,k) &
          !           + a13(5)*wf(i,j+3,k) )*idy(j)*sinphi(i,j)
          ! dp = dp + ( a13(1)*pf(i,j-1,k)+a13(2)*pf(i,j,k) &
          !           + a13(3)*pf(i,j+1,k)+a13(4)*pf(i,j+2,k) &
          !           + a13(5)*pf(i,j+3,k) )*idy(j)*sinphi(i,j)
          ! dr = dr + ( a13(1)*rf(i,j-1,k)+a13(2)*rf(i,j,k) &
          !           + a13(3)*rf(i,j+1,k)+a13(4)*rf(i,j+2,k) &
          !           + a13(5)*rf(i,j+3,k) )*idy(j)*sinphi(i,j)

          ! Non-centered order 1 in x-direction (forward differences)
          du = (uf(i+1,j,k)-uf(i,j,k))*ideltax*cosphi(i,j)
          dv = (vf(i+1,j,k)-vf(i,j,k))*ideltax*cosphi(i,j)
          dw = (wf(i+1,j,k)-wf(i,j,k))*ideltax*cosphi(i,j)
          dp = (pf(i+1,j,k)-pf(i,j,k))*ideltax*cosphi(i,j)
          dr = (rf(i+1,j,k)-rf(i,j,k))*ideltax*cosphi(i,j)

          ! Standard schemes on 3-point stencil
          du =  du + a3(1)*( uf(i,j+1,k)- uf(i,j-1,k))*idy2_jmin*sinphi(i,j)
          dv =  dv + a3(1)*( vf(i,j+1,k)- vf(i,j-1,k))*idy2_jmin*sinphi(i,j)
          dw =  dw + a3(1)*( wf(i,j+1,k)- wf(i,j-1,k))*idy2_jmin*sinphi(i,j)
          dp =  dp + a3(1)*( pf(i,j+1,k)- pf(i,j-1,k))*idy2_jmin*sinphi(i,j)
          dr =  dr + a3(1)*( rf(i,j+1,k)- rf(i,j-1,k))*idy2_jmin*sinphi(i,j)

          pt(j,k) = vg(j,k)*(dp + pf(i,j,k)*ir(i,j))
          ut(j,k) = vg(j,k)*(du + uf(i,j,k)*ir(i,j))
          vt(j,k) = vg(j,k)*(dv + vf(i,j,k)*ir(i,j))
          wt(j,k) = vg(j,k)*(dw + wf(i,j,k)*ir(i,j))
          rt(j,k) = vg(j,k)*(dr + rf(i,j,k)*ir(i,j))
       enddo
    endif


    ! Boundary condition at imin-jmax (edge 1,1,2 /left-top)
    ! ======================================================
    if ((coord(2)==ndomy-1).and.((BC_face(2,2)%sort==-11).or.(BC_face(2,2)%sort==0))) then
       i=1

       ! Pointers for polar coordinates
       ! ==============================
       ir=>BC_edge(1,1,2)%ir
       cosphi=>BC_edge(1,1,2)%cosphi
       sinphi=>BC_edge(1,1,2)%sinphi

       ! Compute fluctuations
       ! ====================
       ! already done before

       ! Compute group velocity vg
       ! =========================
       ! vg=u0*cos(phi)+v0*sin(phi) + sqrt(c0^2-w0^2-(u0*sin(phi)+v0*cos(phi))^2)
       do j=ny-1,ny
          do k=1,nz
             vg(j,k)= BC_face(1,1)%U0(i,j,k,2)*cosphi(i,j-ny+2)+BC_face(1,1)%U0(i,j,k,3)*sinphi(i,j-ny+2) &
                  + sqrt(BC_face(1,1)%U0(i,j,k,6)-BC_face(1,1)%U0(i,j,k,4)**2 &
                  -(BC_face(1,1)%U0(i,j,k,2)*sinphi(i,j-ny+2)-BC_face(1,1)%U0(i,j,k,3)*cosphi(i,j-ny+2))**2)
          enddo
       enddo

       ! Non-centered derivatives
       ! ========================
       if (BC_face(2,2)%sort==-11) then
          i=1;j=ny
          ideltax = 1.0_wp/(x(i+1)-x(i))
          ideltay = 1.0_wp/(y(j)-y(j-1))
          do k=1,nz
             ! (Tam & Webb DRP schemes)
             ! ! Non-centered derivatives in x-direction *cos(phi)
             ! ! =======================================
             ! du = ( a06(1)*uf(1,l,k)+a06(2)*uf(2,l,k) &
             !      + a06(3)*uf(3,l,k)+a06(4)*uf(4,l,k) &
             !      + a06(5)*uf(5,l,k)+a06(6)*uf(6,l,k) &
             !      + a06(7)*uf(7,l,k) ) *idx(i)*cosphi(i,2)
             ! dv = ( a06(1)*vf(1,l,k)+a06(2)*vf(2,l,k) &
             !      + a06(3)*vf(3,l,k)+a06(4)*vf(4,l,k) &
             !      + a06(5)*vf(5,l,k)+a06(6)*vf(6,l,k) &
             !      + a06(7)*vf(7,l,k) ) *idx(i)*cosphi(i,2)
             ! dw = ( a06(1)*wf(1,l,k)+a06(2)*wf(2,l,k) &
             !      + a06(3)*wf(3,l,k)+a06(4)*wf(4,l,k) &
             !      + a06(5)*wf(5,l,k)+a06(6)*wf(6,l,k) &
             !      + a06(7)*wf(7,l,k) ) *idx(i)*cosphi(i,2)
             ! dp = ( a06(1)*pf(1,l,k)+a06(2)*pf(2,l,k) &
             !      + a06(3)*pf(3,l,k)+a06(4)*pf(4,l,k) &
             !      + a06(5)*pf(5,l,k)+a06(6)*pf(6,l,k) &
             !      + a06(7)*pf(7,l,k) ) *idx(i)*cosphi(i,2)
             ! dr = ( a06(1)*rf(1,l,k)+a06(2)*rf(2,l,k) &
             !      + a06(3)*rf(3,l,k)+a06(4)*rf(4,l,k) &
             !      + a06(5)*rf(5,l,k)+a06(6)*rf(6,l,k) &
             !      + a06(7)*rf(7,l,k) ) *idx(i)*cosphi(i,2)

             ! ! Non-centered derivatives in y-direction *sin(phi)
             ! ! =======================================
             ! du = du +( a60(1)*uf(i,l  ,k)+a60(2)*uf(i,l-1,k) &
             !          + a60(3)*uf(i,l-2,k)+a60(4)*uf(i,l-3,k) &
             !          + a60(5)*uf(i,l-4,k)+a60(6)*uf(i,l-5,k) &
             !          + a60(7)*uf(i,l-6,k) ) *idy(j)*sinphi(i,2)
             ! dv = dv +( a60(1)*vf(i,l  ,k)+a60(2)*vf(i,l-1,k) &
             !          + a60(3)*vf(i,l-2,k)+a60(4)*vf(i,l-3,k) &
             !          + a60(5)*vf(i,l-4,k)+a60(6)*vf(i,l-5,k) &
             !          + a60(7)*vf(i,l-6,k) ) *idy(j)*sinphi(i,2)
             ! dw = dw +( a60(1)*wf(i,l  ,k)+a60(2)*wf(i,l-1,k) &
             !          + a60(3)*wf(i,l-2,k)+a60(4)*wf(i,l-3,k) &
             !          + a60(5)*wf(i,l-4,k)+a60(6)*wf(i,l-5,k) &
             !          + a60(7)*wf(i,l-6,k) ) *idy(j)*sinphi(i,2)
             ! dp = dp +( a60(1)*pf(i,l  ,k)+a60(2)*pf(i,l-1,k) &
             !          + a60(3)*pf(i,l-2,k)+a60(4)*pf(i,l-3,k) &
             !          + a60(5)*pf(i,l-4,k)+a60(6)*pf(i,l-5,k) &
             !          + a60(7)*pf(i,l-6,k) ) *idy(j)*sinphi(i,2)
             ! dr = dr +( a60(1)*rf(i,l  ,k)+a60(2)*rf(i,l-1,k) &
             !          + a60(3)*rf(i,l-2,k)+a60(4)*rf(i,l-3,k) &
             !          + a60(5)*rf(i,l-4,k)+a60(6)*rf(i,l-5,k) &
             !          + a60(7)*rf(i,l-6,k) ) *idy(j)*sinphi(i,2)

             ! ! Non-centered derivatives in x-direction *cos(phi)
             ! ! =======================================
             ! du = ( a04(1)*uf(i  ,j,k)+a04(2)*uf(i+1,j,k) &
             !      + a04(3)*uf(i+2,j,k)+a04(4)*uf(i+3,j,k) &
             !      + a04(5)*uf(i+4,j,k) )*idx(i)*cosphi(i,2)
             ! dv = ( a04(1)*vf(i  ,j,k)+a04(2)*vf(i+1,j,k) &
             !      + a04(3)*vf(i+2,j,k)+a04(4)*vf(i+3,j,k) &
             !      + a04(5)*vf(i+4,j,k) )*idx(i)*cosphi(i,2)
             ! dw = ( a04(1)*wf(i  ,j,k)+a04(2)*wf(i+1,j,k) &
             !      + a04(3)*wf(i+2,j,k)+a04(4)*wf(i+3,j,k) &
             !      + a04(5)*wf(i+4,j,k) )*idx(i)*cosphi(i,2)
             ! dp = ( a04(1)*pf(i  ,j,k)+a04(2)*pf(i+1,j,k) &
             !      + a04(3)*pf(i+2,j,k)+a04(4)*pf(i+3,j,k) &
             !      + a04(5)*pf(i+4,j,k) )*idx(i)*cosphi(i,2)
             ! dr = ( a04(1)*rf(i  ,j,k)+a04(2)*rf(i+1,j,k) &
             !      + a04(3)*rf(i+2,j,k)+a04(4)*rf(i+3,j,k) &
             !      + a04(5)*rf(i+4,j,k) )*idx(i)*cosphi(i,2)

             ! ! Non-centered derivatives in y-direction *sin(phi)
             ! ! =======================================
             ! du = du +( a40(1)*uf(i,j,k)  +a40(2)*uf(i,j-1,k) &
             !          + a40(3)*uf(i,j-2,k)+a40(4)*uf(i,j-3,k) &
             !          + a40(5)*uf(i,j-4,k) )*idy(j)*sinphi(i,2)
             ! dv = dv +( a40(1)*vf(i,j,k)  +a40(2)*vf(i,j-1,k) &
             !          + a40(3)*vf(i,j-2,k)+a40(4)*vf(i,j-3,k) &
             !          + a40(5)*vf(i,j-4,k) )*idy(j)*sinphi(i,2)
             ! dw = dw +( a40(1)*wf(i,j,k)  +a40(2)*wf(i,j-1,k) &
             !          + a40(3)*wf(i,j-2,k)+a40(4)*wf(i,j-3,k) &
             !          + a40(5)*wf(i,j-4,k) )*idy(j)*sinphi(i,2)
             ! dp = dp +( a40(1)*pf(i,j,k)  +a40(2)*pf(i,j-1,k) &
             !          + a40(3)*pf(i,j-2,k)+a40(4)*pf(i,j-3,k) &
             !          + a40(5)*pf(i,j-4,k) )*idy(j)*sinphi(i,2)
             ! dr = dr +( a40(1)*rf(i,j,k)  +a40(2)*rf(i,j-1,k) &
             !          + a40(3)*rf(i,j-2,k)+a40(4)*rf(i,j-3,k) &
             !          + a40(5)*rf(i,j-4,k) )*idy(j)*sinphi(i,2)

             ! Non-centered order 1 in x-direction (forward differences)
             du = (uf(i+1,j,k)-uf(i,j,k))*ideltax*cosphi(i,2)
             dv = (vf(i+1,j,k)-vf(i,j,k))*ideltax*cosphi(i,2)
             dw = (wf(i+1,j,k)-wf(i,j,k))*ideltax*cosphi(i,2)
             dp = (pf(i+1,j,k)-pf(i,j,k))*ideltax*cosphi(i,2)
             dr = (rf(i+1,j,k)-rf(i,j,k))*ideltax*cosphi(i,2)

             ! Non-centered order 1 in y-direction (backward differences)
             du = du + (uf(i,j,k) - uf(i,j-1,k))*ideltay*sinphi(i,2)
             dv = dv + (vf(i,j,k) - vf(i,j-1,k))*ideltay*sinphi(i,2)
             dw = dw + (wf(i,j,k) - wf(i,j-1,k))*ideltay*sinphi(i,2)
             dp = dp + (pf(i,j,k) - pf(i,j-1,k))*ideltay*sinphi(i,2)
             dr = dr + (rf(i,j,k) - rf(i,j-1,k))*ideltay*sinphi(i,2)

             pt(j,k) = vg(j,k)*(dp + pf(i,j,k)*ir(i,2))
             ut(j,k) = vg(j,k)*(du + uf(i,j,k)*ir(i,2))
             vt(j,k) = vg(j,k)*(dv + vf(i,j,k)*ir(i,2))
             wt(j,k) = vg(j,k)*(dw + wf(i,j,k)*ir(i,2))
             rt(j,k) = vg(j,k)*(dr + rf(i,j,k)*ir(i,2))
          enddo
       endif

       i=1;j=ny-1
       ideltax = 1.0_wp/(x(i+1)-x(i))
       do k=1,nz
          ! ! Non-centered derivatives in x-direction *cos(phi)
          ! ! =======================================
          ! du = ( a04(1)*uf(i  ,j,k)+a04(2)*uf(i+1,j,k) &
          !      + a04(3)*uf(i+2,j,k)+a04(4)*uf(i+3,j,k) &
          !      + a04(5)*uf(i+4,j,k) )*idx(i)*cosphi(i,1)
          ! dv = ( a04(1)*vf(i  ,j,k)+a04(2)*vf(i+1,j,k) &
          !      + a04(3)*vf(i+2,j,k)+a04(4)*vf(i+3,j,k) &
          !      + a04(5)*vf(i+4,j,k) )*idx(i)*cosphi(i,1)
          ! dw = ( a04(1)*wf(i  ,j,k)+a04(2)*wf(i+1,j,k) &
          !      + a04(3)*wf(i+2,j,k)+a04(4)*wf(i+3,j,k) &
          !      + a04(5)*wf(i+4,j,k) )*idx(i)*cosphi(i,1)
          ! dp = ( a04(1)*pf(i  ,j,k)+a04(2)*pf(i+1,j,k) &
          !      + a04(3)*pf(i+2,j,k)+a04(4)*pf(i+3,j,k) &
          !      + a04(5)*pf(i+4,j,k) )*idx(i)*cosphi(i,1)
          ! dr = ( a04(1)*rf(i  ,j,k)+a04(2)*rf(i+1,j,k) &
          !      + a04(3)*rf(i+2,j,k)+a04(4)*rf(i+3,j,k) &
          !      + a04(5)*rf(i+4,j,k) )*idx(i)*cosphi(i,1)

          ! ! Non-centered derivatives in y-direction *sin(phi)
          ! ! =======================================
          ! du = du + ( a31(1)*uf(i,j+1,k)+a31(2)*uf(i,j,k) &
          !           + a31(3)*uf(i,j-1,k)+a31(4)*uf(i,j-2,k) &
          !           + a31(5)*uf(i,j-3,k) )*idy(j)*sinphi(i,1)
          ! dv = dv + ( a31(1)*vf(i,j+1,k)+a31(2)*vf(i,j,k) &
          !           + a31(3)*vf(i,j-1,k)+a31(4)*vf(i,j-2,k) &
          !           + a31(5)*vf(i,j-3,k) )*idy(j)*sinphi(i,1)
          ! dw = dw + ( a31(1)*wf(i,j+1,k)+a31(2)*wf(i,j,k) &
          !           + a31(3)*wf(i,j-1,k)+a31(4)*wf(i,j-2,k) &
          !           + a31(5)*wf(i,j-3,k) )*idy(j)*sinphi(i,1)
          ! dp = dp + ( a31(1)*pf(i,j+1,k)+a31(2)*pf(i,j,k) &
          !           + a31(3)*pf(i,j-1,k)+a31(4)*pf(i,j-2,k) &
          !           + a31(5)*pf(i,j-3,k) )*idy(j)*sinphi(i,1)
          ! dr = dr + ( a31(1)*rf(i,j+1,k)+a31(2)*rf(i,j,k) &
          !           + a31(3)*rf(i,j-1,k)+a31(4)*rf(i,j-2,k) &
          !           + a31(5)*rf(i,j-3,k) )*idy(j)*sinphi(i,1)

          ! Non-centered order 1 in x-direction (forward differences)
          du = (uf(i+1,j,k)-uf(i,j,k))*ideltax*cosphi(i,1)
          dv = (vf(i+1,j,k)-vf(i,j,k))*ideltax*cosphi(i,1)
          dw = (wf(i+1,j,k)-wf(i,j,k))*ideltax*cosphi(i,1)
          dp = (pf(i+1,j,k)-pf(i,j,k))*ideltax*cosphi(i,1)
          dr = (rf(i+1,j,k)-rf(i,j,k))*ideltax*cosphi(i,1)

          ! Standard schemes on 3-point stencil
          du =  du + a3(1)*( uf(i,j+1,k)- uf(i,j-1,k))*idy2_jmin*sinphi(i,1)
          dv =  dv + a3(1)*( vf(i,j+1,k)- vf(i,j-1,k))*idy2_jmin*sinphi(i,1)
          dw =  dw + a3(1)*( wf(i,j+1,k)- wf(i,j-1,k))*idy2_jmin*sinphi(i,1)
          dp =  dp + a3(1)*( pf(i,j+1,k)- pf(i,j-1,k))*idy2_jmin*sinphi(i,1)
          dr =  dr + a3(1)*( rf(i,j+1,k)- rf(i,j-1,k))*idy2_jmin*sinphi(i,1)

          pt(j,k) = vg(j,k)*(dp + pf(i,j,k)*ir(i,1))
          ut(j,k) = vg(j,k)*(du + uf(i,j,k)*ir(i,1))
          vt(j,k) = vg(j,k)*(dv + vf(i,j,k)*ir(i,1))
          wt(j,k) = vg(j,k)*(dw + wf(i,j,k)*ir(i,1))
          rt(j,k) = vg(j,k)*(dr + rf(i,j,k)*ir(i,1))
       enddo
    endif

    ! Update fluxes at each RK step
    ! =============================
    i=1
    ! pt=0;ut=0;vt=0;wt=0;rt=0;
    ! pt(2,1)=0;ut(2,1)=0;vt(2,1)=0;wt(2,1)=0;rt(2,1)=0
    ! pt(nx-1,1)=0;ut(nx-1,1)=0;vt(nx-1,1)=0;wt(nx-1,1)=0;rt(nx-1,1)=0
    do k=1,nz
       do j=ndy_td1m1,nfy_td1p1
          cp=cpcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
          av=avcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
          Krho(i,j,k)  = rt(j,k)
          Krhou(i,j,k) = uu(i,j,k)*rt(j,k)+rho_n(i,j,k)*ut(j,k)
          Krhov(i,j,k) = vv(i,j,k)*rt(j,k)+rho_n(i,j,k)*vt(j,k)
          Krhow(i,j,k) = ww(i,j,k)*rt(j,k)+rho_n(i,j,k)*wt(j,k)
          Krhoe(i,j,k) = cp/av*(pt(j,k)/c2_(j,k)-rt(j,k)) &
               + (rhoe_n(i,j,k)+prs(i,j,k))/rho_n(i,j,k)*rt(j,k) &
               + rho_n(i,j,k)*(uu(i,j,k)*ut(j,k)+vv(i,j,k)*vt(j,k)+ww(i,j,k)*wt(j,k))
       enddo
    enddo

  end subroutine bc_TD2d_1pt_imin

  !===========================================================================================
  module subroutine bc_TD2d_1pt_imax
  !===========================================================================================
    !> 2D Tam & Dong's BC on 1 point: boundary condition at imax (right) - Cartesian version -
  !===========================================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l,m
    real(wp) :: cp,av,inn1,ntm1,ideltax,ideltay
    real(wp) :: du,dv,dw,dp,dr
    real(wp), dimension(-3:1,ny1:ny2,nz) :: rf,uf,vf,wf,pf
    real(wp), dimension(ny,nz) :: rt,ut,vt,wt,pt,vg,c2_
    !-------------------------------------------------------------------------

    ! Index of right boundary
    ! =======================
    i=nx

    ! Pointers for polar coordinates
    ! ==============================
    ir=>BC_face(1,2)%ir
    cosphi=>BC_face(1,2)%cosphi
    sinphi=>BC_face(1,2)%sinphi

    ! Wall condition applied (different from Tam&Dong 5 points)
    ! ======================
    ! Wall BC at jmin
    ! ---------------
    if (is_bc_wall(2,1)) then
        Krho(i,1,ndz_e:nfz_e)=0.0_wp
       Krhou(i,1,ndz_e:nfz_e)=0.0_wp
       Krhov(i,1,ndz_e:nfz_e)=0.0_wp
       Krhow(i,1,ndz_e:nfz_e)=0.0_wp
       Krhoe(i,1,ndz_e:nfz_e)=0.0_wp
    endif
    ! Wall BC at jmax
    ! ---------------
    if (is_bc_wall(2,2)) then
        Krho(i,ny,ndz_e:nfz_e)=0.0_wp
       Krhou(i,ny,ndz_e:nfz_e)=0.0_wp
       Krhov(i,ny,ndz_e:nfz_e)=0.0_wp
       Krhow(i,ny,ndz_e:nfz_e)=0.0_wp
       Krhoe(i,ny,ndz_e:nfz_e)=0.0_wp
    endif
    ! Wall BC at kmin
    ! ---------------
    if (is_bc_wall(3,1)) then
        Krho(i,ndy_e:nfy_e,1)=0.0_wp
       Krhou(i,ndy_e:nfy_e,1)=0.0_wp
       Krhov(i,ndy_e:nfy_e,1)=0.0_wp
       Krhow(i,ndy_e:nfy_e,1)=0.0_wp
       Krhoe(i,ndy_e:nfy_e,1)=0.0_wp
    endif
    ! Wall BC at kmax
    ! ---------------
    if (is_bc_wall(3,2)) then
        Krho(i,ndy_e:nfy_e,nz)=0.0_wp
       Krhou(i,ndy_e:nfy_e,nz)=0.0_wp
       Krhov(i,ndy_e:nfy_e,nz)=0.0_wp
       Krhow(i,ndy_e:nfy_e,nz)=0.0_wp
       Krhoe(i,ndy_e:nfy_e,nz)=0.0_wp
    endif

    ! Sound speed
    ! ===========
    i=nx
    do k=1,nz
       do j=1,ny
          c2_(j,k)=c2calc_tro(Tmp(i,j,k),rho_n(i,j,k))
       enddo
    enddo

    ! Compute online time-averaged primitive variables
    ! ================================================
    if (irk==nrk) then
       inn1 = 1.0_wp/dble(ntotal)
       ntm1=dble(ntotal-1)
       ! time-averaged primitive variables
       do k=1,nz
          BC_face(1,2)%U0(-3:1,:,k,1)=(ntm1*BC_face(1,2)%U0(-3:1,:,k,1)+rho_n(nx-4:nx,:,k))*inn1
          BC_face(1,2)%U0(-3:1,:,k,2)=(ntm1*BC_face(1,2)%U0(-3:1,:,k,2)+   uu(nx-4:nx,:,k))*inn1
          BC_face(1,2)%U0(-3:1,:,k,3)=(ntm1*BC_face(1,2)%U0(-3:1,:,k,3)+   vv(nx-4:nx,:,k))*inn1
          BC_face(1,2)%U0(-3:1,:,k,4)=(ntm1*BC_face(1,2)%U0(-3:1,:,k,4)+   ww(nx-4:nx,:,k))*inn1
          BC_face(1,2)%U0(-3:1,:,k,5)=(ntm1*BC_face(1,2)%U0(-3:1,:,k,5)+  prs(nx-4:nx,:,k))*inn1
       enddo
       ! time-averaged sound speed squared
       i=nx
       l=i-(nx-1)
       do k=1,nz
          BC_face(1,2)%U0(l,1:ny,k,6)=(ntm1*BC_face(1,2)%U0(l,1:ny,k,6)+c2_(1:ny,k))*inn1
       enddo
    endif

    ! Compute fluctuations
    ! ====================
    do i=nx-4,nx
       l=i-(nx-1)
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
    i=nx
    l=i-(nx-1)
    do j=ndy_td1,nfy_td1
       do k=1,nz
          vg(j,k)= BC_face(1,2)%U0(l,j,k,2)*cosphi(1,j)+BC_face(1,2)%U0(l,j,k,3)*sinphi(1,j) &
               + sqrt(BC_face(1,2)%U0(l,j,k,6)-BC_face(1,2)%U0(l,j,k,4)**2 &
               -(BC_face(1,2)%U0(l,j,k,2)*sinphi(1,j)-BC_face(1,2)%U0(l,j,k,3)*cosphi(1,j))**2)
       enddo
    enddo

    ! Compute vg*[cos(phi)*dq/dx+sin(phi)*dq/dy+q/r]
    ! ==============================================
    i=nx
    l=i-(nx-1)
    ideltax = 1.0_wp/(x(i)-x(i-1))
    do j=ndy_td1,nfy_td1
       do k=1,nz
          ! (Tam & Webb DRP schemes)
          ! ------------------------
          ! du = ( a60(1)*uf(l  ,j,k)+a60(2)*uf(l-1,j,k) &
          !      + a60(3)*uf(l-2,j,k)+a60(4)*uf(l-3,j,k) &
          !      + a60(5)*uf(l-4,j,k)+a60(6)*uf(l-5,j,k) &
          !      + a60(7)*uf(l-6,j,k) ) *idx(i)*cosphi(1,j)
          ! dv = ( a60(1)*vf(l  ,j,k)+a60(2)*vf(l-1,j,k) &
          !      + a60(3)*vf(l-2,j,k)+a60(4)*vf(l-3,j,k) &
          !      + a60(5)*vf(l-4,j,k)+a60(6)*vf(l-5,j,k) &
          !      + a60(7)*vf(l-6,j,k) ) *idx(i)*cosphi(1,j)
          ! dw = ( a60(1)*wf(l  ,j,k)+a60(2)*wf(l-1,j,k) &
          !      + a60(3)*wf(l-2,j,k)+a60(4)*wf(l-3,j,k) &
          !      + a60(5)*wf(l-4,j,k)+a60(6)*wf(l-5,j,k) &
          !      + a60(7)*wf(l-6,j,k) ) *idx(i)*cosphi(1,j)
          ! dp = ( a60(1)*pf(l  ,j,k)+a60(2)*pf(l-1,j,k) &
          !      + a60(3)*pf(l-2,j,k)+a60(4)*pf(l-3,j,k) &
          !      + a60(5)*pf(l-4,j,k)+a60(6)*pf(l-5,j,k) &
          !      + a60(7)*pf(l-6,j,k) ) *idx(i)*cosphi(1,j)
          ! dr = ( a60(1)*rf(l  ,j,k)+a60(2)*rf(l-1,j,k) &
          !      + a60(3)*rf(l-2,j,k)+a60(4)*rf(l-3,j,k) &
          !      + a60(5)*rf(l-4,j,k)+a60(6)*rf(l-5,j,k) &
          !      + a60(7)*rf(l-6,j,k) ) *idx(i)*cosphi(1,j)

          ! du = du + (a7(1)*( uf(l,j+1,k) - uf(l,j-1,k) ) + &
          !            a7(2)*( uf(l,j+2,k) - uf(l,j-2,k) ) + &
          !            a7(3)*( uf(l,j+3,k) - uf(l,j-3,k) ) ) *idy(j)*sinphi(1,j)
          ! dv = dv + (a7(1)*( vf(l,j+1,k) - vf(l,j-1,k) ) + &
          !            a7(2)*( vf(l,j+2,k) - vf(l,j-2,k) ) + &
          !            a7(3)*( vf(l,j+3,k) - vf(l,j-3,k) ) ) *idy(j)*sinphi(1,j)
          ! dw = dw + (a7(1)*( wf(l,j+1,k) - wf(l,j-1,k) ) + &
          !            a7(2)*( wf(l,j+2,k) - wf(l,j-2,k) ) + &
          !            a7(3)*( wf(l,j+3,k) - wf(l,j-3,k) ) ) *idy(j)*sinphi(1,j)
          ! dp = dp + (a7(1)*( pf(l,j+1,k) - pf(l,j-1,k) ) + &
          !            a7(2)*( pf(l,j+2,k) - pf(l,j-2,k) ) + &
          !            a7(3)*( pf(l,j+3,k) - pf(l,j-3,k) ) ) *idy(j)*sinphi(1,j)
          ! dr = dr + (a7(1)*( rf(l,j+1,k) - rf(l,j-1,k) ) + &
          !            a7(2)*( rf(l,j+2,k) - rf(l,j-2,k) ) + &
          !            a7(3)*( rf(l,j+3,k) - rf(l,j-3,k) ) ) *idy(j)*sinphi(1,j)


          ! ! Standard schemes on 5-point stencil
          ! ! -----------------------------------
          ! ! order 4 (decentered 4-0)
          ! du = ( a40(1)*uf(l  ,j,k)+a40(2)*uf(l-1,j,k) &
          !      + a40(3)*uf(l-2,j,k)+a40(4)*uf(l-3,j,k) &
          !      + a40(5)*uf(l-4,j,k) )*idx(i)*cosphi(1,j)
          ! dv = ( a40(1)*vf(l  ,j,k)+a40(2)*vf(l-1,j,k) &
          !      + a40(3)*vf(l-2,j,k)+a40(4)*vf(l-3,j,k) &
          !      + a40(5)*vf(l-4,j,k) )*idx(i)*cosphi(1,j)
          ! dw = ( a40(1)*wf(l  ,j,k)+a40(2)*wf(l-1,j,k) &
          !      + a40(3)*wf(l-2,j,k)+a40(4)*wf(l-3,j,k) &
          !      + a40(5)*wf(l-4,j,k) )*idx(i)*cosphi(1,j)
          ! dp = ( a40(1)*pf(l  ,j,k)+a40(2)*pf(l-1,j,k) &
          !      + a40(3)*pf(l-2,j,k)+a40(4)*pf(l-3,j,k) &
          !      + a40(5)*pf(l-4,j,k) )*idx(i)*cosphi(1,j)
          ! dr = ( a40(1)*rf(l  ,j,k)+a40(2)*rf(l-1,j,k) &
          !      + a40(3)*rf(l-2,j,k)+a40(4)*rf(l-3,j,k) &
          !      + a40(5)*rf(l-4,j,k) )*idx(i)*cosphi(1,j)

          ! ! order 4 (centered)
          ! du = du + ( a5(1)*( uf(l,j+1,k)-uf(l,j-1,k)) &
          !           + a5(2)*( uf(l,j+2,k)-uf(l,j-2,k)) )*idy(j)*sinphi(1,j)
          ! dv = dv + ( a5(1)*( vf(l,j+1,k)-vf(l,j-1,k)) &
          !           + a5(2)*( vf(l,j+2,k)-vf(l,j-2,k)) )*idy(j)*sinphi(1,j)
          ! dw = dw + ( a5(1)*( wf(l,j+1,k)-wf(l,j-1,k)) &
          !           + a5(2)*( wf(l,j+2,k)-wf(l,j-2,k)) )*idy(j)*sinphi(1,j)
          ! dp = dp + ( a5(1)*( pf(l,j+1,k)-pf(l,j-1,k)) &
          !           + a5(2)*( pf(l,j+2,k)-pf(l,j-2,k)) )*idy(j)*sinphi(1,j)
          ! dr = dr + ( a5(1)*( rf(l,j+1,k)-rf(l,j-1,k)) &
          !           + a5(2)*( rf(l,j+2,k)-rf(l,j-2,k)) )*idy(j)*sinphi(1,j)

          ! Non-centered order 1 in x-direction (backward differences)
          du = (uf(l,j,k)-uf(l-1,j,k))*ideltax*cosphi(1,j)
          dv = (vf(l,j,k)-vf(l-1,j,k))*ideltax*cosphi(1,j)
          dw = (wf(l,j,k)-wf(l-1,j,k))*ideltax*cosphi(1,j)
          dp = (pf(l,j,k)-pf(l-1,j,k))*ideltax*cosphi(1,j)
          dr = (rf(l,j,k)-rf(l-1,j,k))*ideltax*cosphi(1,j)

          ! Centered order 4 in y-direction
          du = du + ( a5(1)*( uf(l,j+1,k)-uf(l,j-1,k)) &
                    + a5(2)*( uf(l,j+2,k)-uf(l,j-2,k)) )*idy(j)*sinphi(1,j)
          dv = dv + ( a5(1)*( vf(l,j+1,k)-vf(l,j-1,k)) &
                    + a5(2)*( vf(l,j+2,k)-vf(l,j-2,k)) )*idy(j)*sinphi(1,j)
          dw = dw + ( a5(1)*( wf(l,j+1,k)-wf(l,j-1,k)) &
                    + a5(2)*( wf(l,j+2,k)-wf(l,j-2,k)) )*idy(j)*sinphi(1,j)
          dp = dp + ( a5(1)*( pf(l,j+1,k)-pf(l,j-1,k)) &
                    + a5(2)*( pf(l,j+2,k)-pf(l,j-2,k)) )*idy(j)*sinphi(1,j)
          dr = dr + ( a5(1)*( rf(l,j+1,k)-rf(l,j-1,k)) &
                    + a5(2)*( rf(l,j+2,k)-rf(l,j-2,k)) )*idy(j)*sinphi(1,j)

          pt(j,k) = vg(j,k)*(dp + pf(l,j,k)*ir(1,j))
          ut(j,k) = vg(j,k)*(du + uf(l,j,k)*ir(1,j))
          vt(j,k) = vg(j,k)*(dv + vf(l,j,k)*ir(1,j))
          wt(j,k) = vg(j,k)*(dw + wf(l,j,k)*ir(1,j))
          rt(j,k) = vg(j,k)*(dr + rf(l,j,k)*ir(1,j))
       enddo
    enddo


    ! Boundary condition at imax-jmin (edge 1,2,1 /right-bottom)
    ! ==========================================================
    if ((coord(2)==0).and.((BC_face(2,1)%sort==-11).or.(BC_face(2,1)%sort==0))) then
       ! Pointers for polar coordinates
       ! ==============================
       ir=>BC_edge(1,2,1)%ir
       cosphi=>BC_edge(1,2,1)%cosphi
       sinphi=>BC_edge(1,2,1)%sinphi

       ! Compute fluctuations
       ! ====================
       ! already done before

       ! Compute group velocity vg
       ! =========================
       ! vg=u0*cos(phi)+v0*sin(phi) + sqrt(c0^2-w0^2-(u0*sin(phi)+v0*cos(phi))^2)
       i=nx
       do j=1,2
          do k=1,nz
             vg(j,k)= BC_face(1,2)%U0(l,j,k,2)*cosphi(i-nx+2,j)+BC_face(1,2)%U0(l,j,k,3)*sinphi(i-nx+2,j) &
                  + sqrt(BC_face(1,2)%U0(l,j,k,6)-BC_face(1,2)%U0(l,j,k,4)**2 &
                  -(BC_face(1,2)%U0(l,j,k,2)*sinphi(i-nx+2,j)-BC_face(1,2)%U0(l,j,k,3)*cosphi(i-nx+2,j))**2)
          enddo
       enddo

       ! Non-centered derivatives
       ! ========================
       ! (Tam & Webb DRP schemes)
       if (BC_face(2,1)%sort==-11) then
          i=nx;l=i-(nx-1);j=1
          ideltax = 1.0_wp/(x(i)-x(i-1))
          ideltay = 1.0_wp/(y(j+1)-y(j))
          do k=1,nz
             ! ! Non-centered derivatives in x-direction *cos(phi)
             ! ! =======================================
             ! du = ( a60(1)*uf(l  ,j,k)+a60(2)*uf(l-2,j,k) &
             !      + a60(3)*uf(l-2,j,k)+a60(4)*uf(l-3,j,k) &
             !      + a60(5)*uf(l-4,j,k)+a60(6)*uf(l-5,j,k) &
             !      + a60(7)*uf(l-6,j,k) ) *idx(i)*cosphi(2,j)
             ! dv = ( a60(1)*vf(l  ,j,k)+a60(2)*vf(l-2,j,k) &
             !      + a60(3)*vf(l-2,j,k)+a60(4)*vf(l-3,j,k) &
             !      + a60(5)*vf(l-4,j,k)+a60(6)*vf(l-5,j,k) &
             !      + a60(7)*vf(l-6,j,k) ) *idx(i)*cosphi(2,j)
             ! dw = ( a60(1)*wf(l  ,j,k)+a60(2)*wf(l-2,j,k) &
             !      + a60(3)*wf(l-2,j,k)+a60(4)*wf(l-3,j,k) &
             !      + a60(5)*wf(l-4,j,k)+a60(6)*wf(l-5,j,k) &
             !      + a60(7)*wf(l-6,j,k) ) *idx(i)*cosphi(2,j)
             ! dp = ( a60(1)*pf(l  ,j,k)+a60(2)*pf(l-2,j,k) &
             !      + a60(3)*pf(l-2,j,k)+a60(4)*pf(l-3,j,k) &
             !      + a60(5)*pf(l-4,j,k)+a60(6)*pf(l-5,j,k) &
             !      + a60(7)*pf(l-6,j,k) ) *idx(i)*cosphi(2,j)
             ! dr = ( a60(1)*rf(l  ,j,k)+a60(2)*rf(l-2,j,k) &
             !      + a60(3)*rf(l-2,j,k)+a60(4)*rf(l-3,j,k) &
             !      + a60(5)*rf(l-4,j,k)+a60(6)*rf(l-5,j,k) &
             !      + a60(7)*rf(l-6,j,k) ) *idx(i)*cosphi(2,j)

             ! ! Non-centered derivatives in y-direction *sin(phi)
             ! ! =======================================
             ! du = du + ( a06(1)*uf(l,1,k)+a06(2)*uf(l,2,k) &
             !           + a06(3)*uf(l,3,k)+a06(4)*uf(l,4,k) &
             !           + a06(5)*uf(l,5,k)+a06(6)*uf(l,6,k) &
             !           + a06(7)*uf(l,7,k) ) *idy(j)*sinphi(2,j)
             ! dv = dv + ( a06(1)*vf(l,1,k)+a06(2)*vf(l,2,k) &
             !           + a06(3)*vf(l,3,k)+a06(4)*vf(l,4,k) &
             !           + a06(5)*vf(l,5,k)+a06(6)*vf(l,6,k) &
             !           + a06(7)*vf(l,7,k) ) *idy(j)*sinphi(2,j)
             ! dw = dw + ( a06(1)*wf(l,1,k)+a06(2)*wf(l,2,k) &
             !           + a06(3)*wf(l,3,k)+a06(4)*wf(l,4,k) &
             !           + a06(5)*wf(l,5,k)+a06(6)*wf(l,6,k) &
             !           + a06(7)*wf(l,7,k) ) *idy(j)*sinphi(2,j)
             ! dp = dp + ( a06(1)*pf(l,1,k)+a06(2)*pf(l,2,k) &
             !           + a06(3)*pf(l,3,k)+a06(4)*pf(l,4,k) &
             !           + a06(5)*pf(l,5,k)+a06(6)*pf(l,6,k) &
             !           + a06(7)*pf(l,7,k) ) *idy(j)*sinphi(2,j)
             ! dr = dr + ( a06(1)*rf(l,1,k)+a06(2)*rf(l,2,k) &
             !           + a06(3)*rf(l,3,k)+a06(4)*rf(l,4,k) &
             !           + a06(5)*rf(l,5,k)+a06(6)*rf(l,6,k) &
             !           + a06(7)*rf(l,7,k) ) *idy(j)*sinphi(2,j)

             ! ! Non-centered derivatives in x-direction *cos(phi)
             ! ! =======================================
             ! du = ( a40(1)*uf(l  ,j,k)+a40(2)*uf(l-1,j,k) &
             !      + a40(3)*uf(l-2,j,k)+a40(4)*uf(l-3,j,k) &
             !      + a40(5)*uf(l-4,j,k) )*idx(i)*cosphi(2,j)
             ! dv = ( a40(1)*vf(l  ,j,k)+a40(2)*vf(l-1,j,k) &
             !      + a40(3)*vf(l-2,j,k)+a40(4)*vf(l-3,j,k) &
             !      + a40(5)*vf(l-4,j,k) )*idx(i)*cosphi(2,j)
             ! dw = ( a40(1)*wf(l  ,j,k)+a40(2)*wf(l-1,j,k) &
             !      + a40(3)*wf(l-2,j,k)+a40(4)*wf(l-3,j,k) &
             !      + a40(5)*wf(l-4,j,k) )*idx(i)*cosphi(2,j)
             ! dp = ( a40(1)*pf(l  ,j,k)+a40(2)*pf(l-1,j,k) &
             !      + a40(3)*pf(l-2,j,k)+a40(4)*pf(l-3,j,k) &
             !      + a40(5)*pf(l-4,j,k) )*idx(i)*cosphi(2,j)
             ! dr = ( a40(1)*rf(l  ,j,k)+a40(2)*rf(l-1,j,k) &
             !      + a40(3)*rf(l-2,j,k)+a40(4)*rf(l-3,j,k) &
             !      + a40(5)*rf(l-4,j,k) )*idx(i)*cosphi(2,j)

             ! ! Non-centered derivatives in y-direction *sin(phi)
             ! ! =======================================
             ! du = du + ( a04(1)*uf(l,j  ,k)+a04(2)*uf(l,j+1,k) &
             !           + a04(3)*uf(l,j+2,k)+a04(4)*uf(l,j+3,k) &
             !           + a04(5)*uf(l,j+4,k) )*idy(j)*sinphi(2,j)
             ! dv = dv + ( a04(1)*uf(l,j  ,k)+a04(2)*uf(l,j+1,k) &
             !           + a04(3)*uf(l,j+2,k)+a04(4)*uf(l,j+3,k) &
             !           + a04(5)*uf(l,j+4,k) )*idy(j)*sinphi(2,j)
             ! dw = dw + ( a04(1)*uf(l,j  ,k)+a04(2)*uf(l,j+1,k) &
             !           + a04(3)*uf(l,j+2,k)+a04(4)*uf(l,j+3,k) &
             !           + a04(5)*uf(l,j+4,k) )*idy(j)*sinphi(2,j)
             ! dp = dp + ( a04(1)*uf(l,j  ,k)+a04(2)*uf(l,j+1,k) &
             !           + a04(3)*uf(l,j+2,k)+a04(4)*uf(l,j+3,k) &
             !           + a04(5)*uf(l,j+4,k) )*idy(j)*sinphi(2,j)
             ! dr = dr + ( a04(1)*uf(l,j  ,k)+a04(2)*uf(l,j+1,k) &
             !           + a04(3)*uf(l,j+2,k)+a04(4)*uf(l,j+3,k) &
             !           + a04(5)*uf(l,j+4,k) )*idy(j)*sinphi(2,j)

             ! Non-centered order 1 in x-direction (backward differences)
             du = (uf(l,j,k)-uf(l-1,j,k))*ideltax*cosphi(2,j)
             dv = (vf(l,j,k)-vf(l-1,j,k))*ideltax*cosphi(2,j)
             dw = (wf(l,j,k)-wf(l-1,j,k))*ideltax*cosphi(2,j)
             dp = (pf(l,j,k)-pf(l-1,j,k))*ideltax*cosphi(2,j)
             dr = (rf(l,j,k)-rf(l-1,j,k))*ideltax*cosphi(2,j)

             ! Non-centered order 1 in y-direction (forward differences)
             du = du + (uf(l,j+1,k) - uf(l,j,k))*ideltay*sinphi(2,j)
             dv = dv + (vf(l,j+1,k) - vf(l,j,k))*ideltay*sinphi(2,j)
             dw = dw + (wf(l,j+1,k) - wf(l,j,k))*ideltay*sinphi(2,j)
             dp = dp + (pf(l,j+1,k) - pf(l,j,k))*ideltay*sinphi(2,j)
             dr = dr + (rf(l,j+1,k) - rf(l,j,k))*ideltay*sinphi(2,j)

             pt(j,k) = vg(j,k)*(dp + pf(l,j,k)*ir(2,j))
             ut(j,k) = vg(j,k)*(du + uf(l,j,k)*ir(2,j))
             vt(j,k) = vg(j,k)*(dv + vf(l,j,k)*ir(2,j))
             wt(j,k) = vg(j,k)*(dw + wf(l,j,k)*ir(2,j))
             rt(j,k) = vg(j,k)*(dr + rf(l,j,k)*ir(2,j))
          enddo
       endif

       i=nx;l=i-(nx-1);j=2
       ideltax = 1.0_wp/(x(i)-x(i-1))
       do k=1,nz
          ! ! Non-centered derivatives in x-direction *cos(phi)
          ! ! =======================================
          ! du = ( a40(1)*uf(l  ,j,k)+a40(2)*uf(l-1,j,k) &
          !      + a40(3)*uf(l-2,j,k)+a40(4)*uf(l-3,j,k) &
          !      + a40(5)*uf(l-4,j,k) )*idx(i)*cosphi(1,j)
          ! dv = ( a40(1)*vf(l  ,j,k)+a40(2)*vf(l-1,j,k) &
          !      + a40(3)*vf(l-2,j,k)+a40(4)*vf(l-3,j,k) &
          !      + a40(5)*vf(l-4,j,k) )*idx(i)*cosphi(1,j)
          ! dw = ( a40(1)*wf(l  ,j,k)+a40(2)*wf(l-1,j,k) &
          !      + a40(3)*wf(l-2,j,k)+a40(4)*wf(l-3,j,k) &
          !      + a40(5)*wf(l-4,j,k) )*idx(i)*cosphi(1,j)
          ! dp = ( a40(1)*pf(l  ,j,k)+a40(2)*pf(l-1,j,k) &
          !      + a40(3)*pf(l-2,j,k)+a40(4)*pf(l-3,j,k) &
          !      + a40(5)*pf(l-4,j,k) )*idx(i)*cosphi(1,j)
          ! dr = ( a40(1)*rf(l  ,j,k)+a40(2)*rf(l-1,j,k) &
          !      + a40(3)*rf(l-2,j,k)+a40(4)*rf(l-3,j,k) &
          !      + a40(5)*rf(l-4,j,k) )*idx(i)*cosphi(1,j)

          ! ! Non-centered derivatives in y-direction *sin(phi)
          ! ! =======================================
          ! du = du + ( a13(1)*uf(l,j-1,k)+a13(2)*uf(l,j,k) &
          !           + a13(3)*uf(l,j+1,k)+a13(4)*uf(l,j+2,k) &
          !           + a13(5)*uf(l,j+3,k) )*sinphi(1,j)
          ! dv = dv + ( a13(1)*vf(l,j-1,k)+a13(2)*vf(l,j,k) &
          !           + a13(3)*vf(l,j+1,k)+a13(4)*vf(l,j+2,k) &
          !           + a13(5)*vf(l,j+3,k) )*sinphi(1,j)
          ! dw = dw + ( a13(1)*wf(l,j-1,k)+a13(2)*wf(l,j,k) &
          !           + a13(3)*wf(l,j+1,k)+a13(4)*wf(l,j+2,k) &
          !           + a13(5)*wf(l,j+3,k) )*sinphi(1,j)
          ! dp = dp + ( a13(1)*pf(l,j-1,k)+a13(2)*pf(l,j,k) &
          !           + a13(3)*pf(l,j+1,k)+a13(4)*pf(l,j+2,k) &
          !           + a13(5)*pf(l,j+3,k) )*sinphi(1,j)
          ! dr = dr + ( a13(1)*rf(l,j-1,k)+a13(2)*rf(l,j,k) &
          !           + a13(3)*rf(l,j+1,k)+a13(4)*rf(l,j+2,k) &
          !           + a13(5)*rf(l,j+3,k) )*sinphi(1,j)

          ! Non-centered order 1 in x-direction (backward differences)
          du = (uf(l,j,k)-uf(l-1,j,k))*ideltax*cosphi(1,j)
          dv = (vf(l,j,k)-vf(l-1,j,k))*ideltax*cosphi(1,j)
          dw = (wf(l,j,k)-wf(l-1,j,k))*ideltax*cosphi(1,j)
          dp = (pf(l,j,k)-pf(l-1,j,k))*ideltax*cosphi(1,j)
          dr = (rf(l,j,k)-rf(l-1,j,k))*ideltax*cosphi(1,j)

          ! Standard schemes on 3-point stencil
          du =  du + a3(1)*( uf(l,j+1,k)- uf(l,j-1,k))*idy2_jmin*sinphi(1,j)
          dv =  dv + a3(1)*( vf(l,j+1,k)- vf(l,j-1,k))*idy2_jmin*sinphi(1,j)
          dw =  dw + a3(1)*( wf(l,j+1,k)- wf(l,j-1,k))*idy2_jmin*sinphi(1,j)
          dp =  dp + a3(1)*( pf(l,j+1,k)- pf(l,j-1,k))*idy2_jmin*sinphi(1,j)
          dr =  dr + a3(1)*( rf(l,j+1,k)- rf(l,j-1,k))*idy2_jmin*sinphi(1,j)

          pt(j,k) = vg(j,k)*(dp + pf(l,j,k)*ir(1,j))
          ut(j,k) = vg(j,k)*(du + uf(l,j,k)*ir(1,j))
          vt(j,k) = vg(j,k)*(dv + vf(l,j,k)*ir(1,j))
          wt(j,k) = vg(j,k)*(dw + wf(l,j,k)*ir(1,j))
          rt(j,k) = vg(j,k)*(dr + rf(l,j,k)*ir(1,j))
       enddo
    endif

    ! Boundary condition at imax-jmax (edge 1,2,2 /right-bottom)
    ! ==========================================================
    if ((coord(2)==ndomy-1).and.((BC_face(2,2)%sort==-11).or.(BC_face(2,2)%sort==0))) then
       ! Pointers for polar coordinates
       ! ==============================
       ir=>BC_edge(1,2,2)%ir
       cosphi=>BC_edge(1,2,2)%cosphi
       sinphi=>BC_edge(1,2,2)%sinphi

       ! Compute fluctuations
       ! ====================
       ! already done before

       ! Compute group velocity vg
       ! =========================
       ! vg=u0*cos(phi)+v0*sin(phi) + sqrt(c0^2-w0^2-(u0*sin(phi)+v0*cos(phi))^2)
       i=nx
       do j=ny-1,ny
          do k=1,nz
             vg(j,k)= BC_face(1,2)%U0(l,j,k,2)*cosphi(i-nx+2,j-ny+2)+BC_face(1,2)%U0(l,j,k,3)*sinphi(i-nx+2,j-ny+2) &
                  + sqrt(BC_face(1,2)%U0(l,j,k,6)-BC_face(1,2)%U0(l,j,k,4)**2 &
                  -(BC_face(1,2)%U0(l,j,k,2)*sinphi(i-nx+2,j-ny+2)-BC_face(1,2)%U0(l,j,k,3)*cosphi(i-nx+2,j-ny+2))**2)
          enddo
       enddo


       ! Non-centered derivatives
       ! ========================
       if (BC_face(2,2)%sort==-11) then
          i=nx;l=i-(nx-1);j=ny;m=j-(ny-1)
          ideltax = 1.0_wp/(x(i)-x(i-1))
          ideltay = 1.0_wp/(y(j)-y(j-1))
          do k=1,nz
             ! ! Non-centered derivatives in x-direction *cos(phi)
             ! ! =======================================
             ! du = ( a60(1)*uf(l  ,j,k)+a60(2)*uf(l-1,j,k) &
             !      + a60(3)*uf(l-2,j,k)+a60(4)*uf(l-3,j,k) &
             !      + a60(5)*uf(l-4,j,k)+a60(6)*uf(l-5,j,k) &
             !      + a60(7)*uf(l-6,j,k) ) *idx(i)*cosphi(2,2)
             ! dv = ( a60(1)*vf(l  ,j,k)+a60(2)*vf(l-1,j,k) &
             !      + a60(3)*vf(l-2,j,k)+a60(4)*vf(l-3,j,k) &
             !      + a60(5)*vf(l-4,j,k)+a60(6)*vf(l-5,j,k) &
             !      + a60(7)*vf(l-6,j,k) ) *idx(i)*cosphi(2,2)
             ! dw = ( a60(1)*wf(l  ,j,k)+a60(2)*wf(l-1,j,k) &
             !      + a60(3)*wf(l-2,j,k)+a60(4)*wf(l-3,j,k) &
             !      + a60(5)*wf(l-4,j,k)+a60(6)*wf(l-5,j,k) &
             !      + a60(7)*wf(l-6,j,k) ) *idx(i)*cosphi(2,2)
             ! dp = ( a60(1)*pf(l  ,j,k)+a60(2)*pf(l-1,j,k) &
             !      + a60(3)*pf(l-2,j,k)+a60(4)*pf(l-3,j,k) &
             !      + a60(5)*pf(l-4,j,k)+a60(6)*pf(l-5,j,k) &
             !      + a60(7)*pf(l-6,j,k) ) *idx(i)*cosphi(2,2)
             ! dr = ( a60(1)*rf(l  ,j,k)+a60(2)*rf(l-1,j,k) &
             !      + a60(3)*rf(l-2,j,k)+a60(4)*rf(l-3,j,k) &
             !      + a60(5)*rf(l-4,j,k)+a60(6)*rf(l-5,j,k) &
             !      + a60(7)*rf(l-6,j,k) ) *idx(i)*cosphi(2,2)

             ! ! Non-centered derivatives in y-direction *sin(phi)
             ! ! =======================================
             ! du = du + ( a60(1)*uf(l,j  ,k)+a60(2)*uf(l,j-1,k) &
             !           + a60(3)*uf(l,j-2,k)+a60(4)*uf(l,j-3,k) &
             !           + a60(5)*uf(l,j-4,k)+a60(6)*uf(l,j-5,k) &
             !           + a60(7)*uf(l,j-6,k) ) *idy(j)*sinphi(1,1)
             ! dv = dv + ( a60(1)*vf(l,j  ,k)+a60(2)*vf(l,j-1,k) &
             !           + a60(3)*vf(l,j-2,k)+a60(4)*vf(l,j-3,k) &
             !           + a60(5)*vf(l,j-4,k)+a60(6)*vf(l,j-5,k) &
             !           + a60(7)*vf(l,j-6,k) ) *idy(j)*sinphi(1,1)
             ! dw = dw + ( a60(1)*wf(l,j  ,k)+a60(2)*wf(l,j-1,k) &
             !           + a60(3)*wf(l,j-2,k)+a60(4)*wf(l,j-3,k) &
             !           + a60(5)*wf(l,j-4,k)+a60(6)*wf(l,j-5,k) &
             !           + a60(7)*wf(l,j-6,k) ) *idy(j)*sinphi(1,1)
             ! dp = dp + ( a60(1)*pf(l,j  ,k)+a60(2)*pf(l,j-1,k) &
             !           + a60(3)*pf(l,j-2,k)+a60(4)*pf(l,j-3,k) &
             !           + a60(5)*pf(l,j-4,k)+a60(6)*pf(l,j-5,k) &
             !           + a60(7)*pf(l,j-6,k) ) *idy(j)*sinphi(1,1)
             ! dr = dr + ( a60(1)*rf(l,j  ,k)+a60(2)*rf(l,j-1,k) &
             !           + a60(3)*rf(l,j-2,k)+a60(4)*rf(l,j-3,k) &
             !           + a60(5)*rf(l,j-4,k)+a60(6)*rf(l,j-5,k) &
             !           + a60(7)*rf(l,j-6,k) ) *idy(j)*sinphi(1,1)


             ! ! Non-centered derivatives in x-direction *cos(phi)
             ! ! =======================================
             ! du = ( a40(1)*uf(l  ,j,k)+a40(2)*uf(l-1,j,k) &
             !      + a40(3)*uf(l-2,j,k)+a40(4)*uf(l-3,j,k) &
             !      + a40(5)*uf(l-4,j,k) )*idx(i)*cosphi(2,2)
             ! dv = ( a40(1)*vf(l  ,j,k)+a40(2)*vf(l-1,j,k) &
             !      + a40(3)*vf(l-2,j,k)+a40(4)*vf(l-3,j,k) &
             !      + a40(5)*vf(l-4,j,k) )*idx(i)*cosphi(2,2)
             ! dw = ( a40(1)*wf(l  ,j,k)+a40(2)*wf(l-1,j,k) &
             !      + a40(3)*wf(l-2,j,k)+a40(4)*wf(l-3,j,k) &
             !      + a40(5)*wf(l-4,j,k) )*idx(i)*cosphi(2,2)
             ! dp = ( a40(1)*pf(l  ,j,k)+a40(2)*pf(l-1,j,k) &
             !      + a40(3)*pf(l-2,j,k)+a40(4)*pf(l-3,j,k) &
             !      + a40(5)*pf(l-4,j,k) )*idx(i)*cosphi(2,2)
             ! dr = ( a40(1)*rf(l  ,j,k)+a40(2)*rf(l-1,j,k) &
             !      + a40(3)*rf(l-2,j,k)+a40(4)*rf(l-3,j,k) &
             !      + a40(5)*rf(l-4,j,k) )*idx(i)*cosphi(2,2)

             ! ! Non-centered derivatives in y-direction *sin(phi)
             ! ! =======================================
             ! du = du +( a40(1)*uf(l,j,k)  +a40(2)*uf(l,j-1,k) &
             !          + a40(3)*uf(l,j-2,k)+a40(4)*uf(l,j-3,k) &
             !          + a40(5)*uf(l,j-4,k) )*idy(j)*sinphi(2,2)
             ! dv = dv +( a40(1)*vf(l,j,k)  +a40(2)*vf(l,j-1,k) &
             !          + a40(3)*vf(l,j-2,k)+a40(4)*vf(l,j-3,k) &
             !          + a40(5)*vf(l,j-4,k) )*idy(j)*sinphi(2,2)
             ! dw = dw +( a40(1)*wf(l,j,k)  +a40(2)*wf(l,j-1,k) &
             !          + a40(3)*wf(l,j-2,k)+a40(4)*wf(l,j-3,k) &
             !          + a40(5)*wf(l,j-4,k) )*idy(j)*sinphi(2,2)
             ! dp = dp +( a40(1)*pf(l,j,k)  +a40(2)*pf(l,j-1,k) &
             !          + a40(3)*pf(l,j-2,k)+a40(4)*pf(l,j-3,k) &
             !          + a40(5)*pf(l,j-4,k) )*idy(j)*sinphi(2,2)
             ! dr = dr +( a40(1)*rf(l,j,k)  +a40(2)*rf(l,j-1,k) &
             !          + a40(3)*rf(l,j-2,k)+a40(4)*rf(l,j-3,k) &
             !          + a40(5)*rf(l,j-4,k) )*idy(j)*sinphi(2,2)

             ! Non-centered order 1 in x-direction (backward differences)
             du = (uf(l,j,k)-uf(l-1,j,k))*ideltax*cosphi(2,2)
             dv = (vf(l,j,k)-vf(l-1,j,k))*ideltax*cosphi(2,2)
             dw = (wf(l,j,k)-wf(l-1,j,k))*ideltax*cosphi(2,2)
             dp = (pf(l,j,k)-pf(l-1,j,k))*ideltax*cosphi(2,2)
             dr = (rf(l,j,k)-rf(l-1,j,k))*ideltax*cosphi(2,2)

             ! Non-centered order 1 in y-direction (backward differences)
             du = du + (uf(l,j,k) - uf(l,j-1,k))*ideltay*sinphi(2,2)
             dv = dv + (vf(l,j,k) - vf(l,j-1,k))*ideltay*sinphi(2,2)
             dw = dw + (wf(l,j,k) - wf(l,j-1,k))*ideltay*sinphi(2,2)
             dp = dp + (pf(l,j,k) - pf(l,j-1,k))*ideltay*sinphi(2,2)
             dr = dr + (rf(l,j,k) - rf(l,j-1,k))*ideltay*sinphi(2,2)

             pt(j,k) = vg(j,k) * (dp + pf(l,j,k)*ir(2,2))
             ut(j,k) = vg(j,k) * (du + uf(l,j,k)*ir(2,2))
             vt(j,k) = vg(j,k) * (dv + vf(l,j,k)*ir(2,2))
             wt(j,k) = vg(j,k) * (dw + wf(l,j,k)*ir(2,2))
             rt(j,k) = vg(j,k) * (dr + rf(l,j,k)*ir(2,2))
          enddo
       endif

       i=nx;l=i-(nx-1);j=ny-1;m=j-(ny-1)
       do k=1,nz
          ! ! Non-centered derivatives in x-direction *cos(phi)
          ! ! =======================================
          ! du = ( a40(1)*uf(l  ,j,k)+a40(2)*uf(l-1,j,k) &
          !      + a40(3)*uf(l-2,j,k)+a40(4)*uf(l-3,j,k) &
          !      + a40(5)*uf(l-4,j,k) )*idx(i)*cosphi(2,1)
          ! dv = ( a40(1)*vf(l  ,j,k)+a40(2)*vf(l-1,j,k) &
          !      + a40(3)*vf(l-2,j,k)+a40(4)*vf(l-3,j,k) &
          !      + a40(5)*vf(l-4,j,k) )*idx(i)*cosphi(2,1)
          ! dw = ( a40(1)*wf(l  ,j,k)+a40(2)*wf(l-1,j,k) &
          !      + a40(3)*wf(l-2,j,k)+a40(4)*wf(l-3,j,k) &
          !      + a40(5)*wf(l-4,j,k) )*idx(i)*cosphi(2,1)
          ! dp = ( a40(1)*pf(l  ,j,k)+a40(2)*pf(l-1,j,k) &
          !      + a40(3)*pf(l-2,j,k)+a40(4)*pf(l-3,j,k) &
          !      + a40(5)*pf(l-4,j,k) )*idx(i)*cosphi(2,1)
          ! dr = ( a40(1)*rf(l  ,j,k)+a40(2)*rf(l-1,j,k) &
          !      + a40(3)*rf(l-2,j,k)+a40(4)*rf(l-3,j,k) &
          !      + a40(5)*rf(l-4,j,k) )*idx(i)*cosphi(2,1)

          ! ! Non-centered derivatives in y-direction *sin(phi)
          ! ! =======================================
          ! du = du + ( a31(1)*uf(l,j+1,k)+a31(2)*uf(l,j,k) &
          !           + a31(3)*uf(l,j-1,k)+a31(4)*uf(l,j-2,k) &
          !           + a31(5)*uf(l,j-3,k) )*idy(j)*sinphi(2,1)
          ! dv = dv + ( a31(1)*vf(l,j+1,k)+a31(2)*vf(l,j,k) &
          !           + a31(3)*vf(l,j-1,k)+a31(4)*vf(l,j-2,k) &
          !           + a31(5)*vf(l,j-3,k) )*idy(j)*sinphi(2,1)
          ! dw = dw + ( a31(1)*wf(l,j+1,k)+a31(2)*wf(l,j,k) &
          !           + a31(3)*wf(l,j-1,k)+a31(4)*wf(l,j-2,k) &
          !           + a31(5)*wf(l,j-3,k) )*idy(j)*sinphi(2,1)
          ! dp = dp + ( a31(1)*pf(l,j+1,k)+a31(2)*pf(l,j,k) &
          !           + a31(3)*pf(l,j-1,k)+a31(4)*pf(l,j-2,k) &
          !           + a31(5)*pf(l,j-3,k) )*idy(j)*sinphi(2,1)
          ! dr = dr + ( a31(1)*rf(l,j+1,k)+a31(2)*rf(l,j,k) &
          !           + a31(3)*rf(l,j-1,k)+a31(4)*rf(l,j-2,k) &
          !           + a31(5)*rf(l,j-3,k) )*idy(j)*sinphi(2,1)

          ! Non-centered order 1 in x-direction (backward differences)
          du = (uf(l,j,k)-uf(l-1,j,k))*ideltax*cosphi(2,1)
          dv = (vf(l,j,k)-vf(l-1,j,k))*ideltax*cosphi(2,1)
          dw = (wf(l,j,k)-wf(l-1,j,k))*ideltax*cosphi(2,1)
          dp = (pf(l,j,k)-pf(l-1,j,k))*ideltax*cosphi(2,1)
          dr = (rf(l,j,k)-rf(l-1,j,k))*ideltax*cosphi(2,1)

          ! Standard schemes on 3-point stencil
          du =  du + a3(1)*( uf(l,j+1,k)- uf(l,j-1,k))*idy2_jmin*sinphi(2,1)
          dv =  dv + a3(1)*( vf(l,j+1,k)- vf(l,j-1,k))*idy2_jmin*sinphi(2,1)
          dw =  dw + a3(1)*( wf(l,j+1,k)- wf(l,j-1,k))*idy2_jmin*sinphi(2,1)
          dp =  dp + a3(1)*( pf(l,j+1,k)- pf(l,j-1,k))*idy2_jmin*sinphi(2,1)
          dr =  dr + a3(1)*( rf(l,j+1,k)- rf(l,j-1,k))*idy2_jmin*sinphi(2,1)

          pt(j,k) = vg(j,k) * (dp + pf(l,j,k)*ir(2,1))
          ut(j,k) = vg(j,k) * (du + uf(l,j,k)*ir(2,1))
          vt(j,k) = vg(j,k) * (dv + vf(l,j,k)*ir(2,1))
          wt(j,k) = vg(j,k) * (dw + wf(l,j,k)*ir(2,1))
          rt(j,k) = vg(j,k) * (dr + rf(l,j,k)*ir(2,1))
       enddo
    endif

    ! Update fluxes at each RK step
    ! =============================
    i=nx
    do k=1,nz
       do j=ndy_td1m1,nfy_td1p1
          cp=cpcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
          av=avcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
          Krho(i,j,k)  = rt(j,k)
          Krhou(i,j,k) = uu(i,j,k)*rt(j,k)+rho_n(i,j,k)*ut(j,k)
          Krhov(i,j,k) = vv(i,j,k)*rt(j,k)+rho_n(i,j,k)*vt(j,k)
          Krhow(i,j,k) = ww(i,j,k)*rt(j,k)+rho_n(i,j,k)*wt(j,k)
          Krhoe(i,j,k) = cp/av*(pt(j,k)/c2_(j,k)-rt(j,k)) &
               + (rhoe_n(i,j,k)+prs(i,j,k))/rho_n(i,j,k)*rt(j,k) &
               + rho_n(i,j,k)*(uu(i,j,k)*ut(j,k)+vv(i,j,k)*vt(j,k)+ww(i,j,k)*wt(j,k))
       enddo
    enddo

  end subroutine bc_TD2d_1pt_imax

  !============================================================================================
  module subroutine bc_TD2d_1pt_jmin
  !============================================================================================
    !> 2D Tam & Dong's BC on 1 point: boundary condition at jmin (bottom) - Cartesian version -
  !============================================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: cp,av,inn1,ntm1
    real(wp) :: du,dv,dw,dp,dr
    real(wp), dimension(nx1:nx2,1:5,nz) :: rf,uf,vf,wf,pf
    real(wp), dimension(nx,nz) :: rt,ut,vt,wt,pt,vg,c2_
    !-------------------------------------------------------------------------

    ! Index of bottom boundary
    ! ========================
    j=1

    ! Pointers for polar coordinates
    ! ==============================
    ir=>BC_face(2,1)%ir
    cosphi=>BC_face(2,1)%cosphi
    sinphi=>BC_face(2,1)%sinphi

    ! Wall condition applied (different from Tam&Dong 5 points)
    ! ======================
    ! Wall BC at imin
    ! ---------------
    if (is_bc_wall(1,1)) then
        Krho(1,j,ndz_e:nfz_e)=0.0_wp
       Krhou(1,j,ndz_e:nfz_e)=0.0_wp
       Krhov(1,j,ndz_e:nfz_e)=0.0_wp
       Krhow(1,j,ndz_e:nfz_e)=0.0_wp
       Krhoe(1,j,ndz_e:nfz_e)=0.0_wp
    endif
    ! Wall BC at imax
    ! ---------------
    if (is_bc_wall(1,2)) then
        Krho(nx,j,ndz_e:nfz_e)=0.0_wp
       Krhou(nx,j,ndz_e:nfz_e)=0.0_wp
       Krhov(nx,j,ndz_e:nfz_e)=0.0_wp
       Krhow(nx,j,ndz_e:nfz_e)=0.0_wp
       Krhoe(nx,j,ndz_e:nfz_e)=0.0_wp
    endif
    ! Wall BC at kmin
    ! ---------------
    if (is_bc_wall(3,1)) then
        Krho(ndx_e:nfx_e,j,1)=0.0_wp
       Krhou(ndx_e:nfx_e,j,1)=0.0_wp
       Krhov(ndx_e:nfx_e,j,1)=0.0_wp
       Krhow(ndx_e:nfx_e,j,1)=0.0_wp
       Krhoe(ndx_e:nfx_e,j,1)=0.0_wp
    endif
    ! Wall BC at kmax
    ! ---------------
    if (is_bc_wall(3,2)) then
        Krho(ndx_e:nfx_e,j,nz)=0.0_wp
       Krhou(ndx_e:nfx_e,j,nz)=0.0_wp
       Krhov(ndx_e:nfx_e,j,nz)=0.0_wp
       Krhow(ndx_e:nfx_e,j,nz)=0.0_wp
       Krhoe(ndx_e:nfx_e,j,nz)=0.0_wp
    endif

    ! Sound speed
    ! ===========
    j=1
    do k=1,nz
       do i=1,nx
          c2_(i,k)=c2calc_tro(Tmp(i,j,k),rho_n(i,j,k))
       enddo
    enddo

    ! Compute online time-averaged primitive variables
    ! ================================================
    if (irk==nrk) then
       inn1 = 1.0_wp/dble(ntotal)
       ntm1=dble(ntotal-1)
       ! time-averaged primitive variables
       do k=1,nz
          BC_face(2,1)%U0(:,1:5,k,1)=(ntm1*BC_face(2,1)%U0(:,1:5,k,1)+rho_n(:,1:5,k))*inn1
          BC_face(2,1)%U0(:,1:5,k,2)=(ntm1*BC_face(2,1)%U0(:,1:5,k,2)+   uu(:,1:5,k))*inn1
          BC_face(2,1)%U0(:,1:5,k,3)=(ntm1*BC_face(2,1)%U0(:,1:5,k,3)+   vv(:,1:5,k))*inn1
          BC_face(2,1)%U0(:,1:5,k,4)=(ntm1*BC_face(2,1)%U0(:,1:5,k,4)+   ww(:,1:5,k))*inn1
          BC_face(2,1)%U0(:,1:5,k,5)=(ntm1*BC_face(2,1)%U0(:,1:5,k,5)+  prs(:,1:5,k))*inn1
       enddo
       ! time-averaged sound speed squared
       j=1
       do k=1,nz
          BC_face(2,1)%U0(1:nx,j,k,6)=(ntm1*BC_face(2,1)%U0(1:nx,j,k,6)+c2_(1:nx,k))*inn1
       enddo
    endif

    ! Compute fluctuations
    ! ====================
    do k=1,nz
       do j=1,5
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
    ! vg=u0*cos(phi)+v0*sin(phi) + sqrt(c0^2-w0^2-(u0*sin(phi)+v0*cos(phi))^2)
    j=1
    do k=1,nz
       do i=ndx_td1,nfx_td1
          vg(i,k)= BC_face(2,1)%U0(i,j,k,2)*cosphi(i,j)+BC_face(2,1)%U0(i,j,k,3)*sinphi(i,j) &
               + sqrt(BC_face(2,1)%U0(i,j,k,6)-BC_face(2,1)%U0(i,j,k,4)**2 &
               -(BC_face(2,1)%U0(i,j,k,2)*sinphi(i,j)-BC_face(2,1)%U0(i,j,k,3)*cosphi(i,j))**2)
       enddo
    enddo

    ! Compute vg*[cos(phi)*dq/dx+sin(phi)*dq/dy+q/r]
    ! ==============================================
    ! (Tam & Webb DRP schemes)
    j=1
    do k=1,nz
       do i=ndx_td1,nfx_td1
          ! du = ( a7(1)*( uf(i+1,j,k) - uf(i-1,j,k) ) + &
          !        a7(2)*( uf(i+2,j,k) - uf(i-2,j,k) ) + &
          !        a7(3)*( uf(i+3,j,k) - uf(i-3,j,k) ) ) *idx(i)*cosphi(i,j)
          ! dv = ( a7(1)*( vf(i+1,j,k) - vf(i-1,j,k) ) + &
          !        a7(2)*( vf(i+2,j,k) - vf(i-2,j,k) ) + &
          !        a7(3)*( vf(i+3,j,k) - vf(i-3,j,k) ) ) *idx(i)*cosphi(i,j)
          ! dw = ( a7(1)*( wf(i+1,j,k) - wf(i-1,j,k) ) + &
          !        a7(2)*( wf(i+2,j,k) - wf(i-2,j,k) ) + &
          !        a7(3)*( wf(i+3,j,k) - wf(i-3,j,k) ) ) *idx(i)*cosphi(i,j)
          ! dp = ( a7(1)*( pf(i+1,j,k) - pf(i-1,j,k) ) + &
          !        a7(2)*( pf(i+2,j,k) - pf(i-2,j,k) ) + &
          !        a7(3)*( pf(i+3,j,k) - pf(i-3,j,k) ) ) *idx(i)*cosphi(i,j)
          ! dr = ( a7(1)*( rf(i+1,j,k) - rf(i-1,j,k) ) + &
          !        a7(2)*( rf(i+2,j,k) - rf(i-2,j,k) ) + &
          !        a7(3)*( rf(i+3,j,k) - rf(i-3,j,k) ) ) *idx(i)*cosphi(i,j)

          ! du = du + ( a06(1)*uf(i,1,k)+a06(2)*uf(i,2,k) &
          !           + a06(3)*uf(i,3,k)+a06(4)*uf(i,4,k) &
          !           + a06(5)*uf(i,5,k)+a06(6)*uf(i,6,k) &
          !           + a06(7)*uf(i,7,k) ) *idy(j)*sinphi(i,j)
          ! dv = dv + ( a06(1)*vf(i,1,k)+a06(2)*vf(i,2,k) &
          !           + a06(3)*vf(i,3,k)+a06(4)*vf(i,4,k) &
          !           + a06(5)*vf(i,5,k)+a06(6)*vf(i,6,k) &
          !           + a06(7)*vf(i,7,k) ) *idy(j)*sinphi(i,j)
          ! dw = dw + ( a06(1)*wf(i,1,k)+a06(2)*wf(i,2,k) &
          !           + a06(3)*wf(i,3,k)+a06(4)*wf(i,4,k) &
          !           + a06(5)*wf(i,5,k)+a06(6)*wf(i,6,k) &
          !           + a06(7)*wf(i,7,k) ) *idy(j)*sinphi(i,j)
          ! dp = dp + ( a06(1)*pf(i,1,k)+a06(2)*pf(i,2,k) &
          !           + a06(3)*pf(i,3,k)+a06(4)*pf(i,4,k) &
          !           + a06(5)*pf(i,5,k)+a06(6)*pf(i,6,k) &
          !           + a06(7)*pf(i,7,k) ) *idy(j)*sinphi(i,j)
          ! dr = dr + ( a06(1)*rf(i,1,k)+a06(2)*rf(i,2,k) &
          !           + a06(3)*rf(i,3,k)+a06(4)*rf(i,4,k) &
          !           + a06(5)*rf(i,5,k)+a06(6)*rf(i,6,k) &
          !           + a06(7)*rf(i,7,k) ) *idy(j)*sinphi(i,j)


          du = ( a5(1)*( uf(i+1,j,k)-uf(i-1,j,k) ) &
               + a5(2)*( uf(i+2,j,k)-uf(i-2,j,k) ) )*idx(i)*cosphi(i,j)
          dv = ( a5(1)*( vf(i+1,j,k)-vf(i-1,j,k) ) &
               + a5(2)*( vf(i+2,j,k)-vf(i-2,j,k) ) )*idx(i)*cosphi(i,j)
          dw = ( a5(1)*( wf(i+1,j,k)-wf(i-1,j,k) ) &
               + a5(2)*( wf(i+2,j,k)-wf(i-2,j,k) ) )*idx(i)*cosphi(i,j)
          dp = ( a5(1)*( pf(i+1,j,k)-pf(i-1,j,k) ) &
               + a5(2)*( pf(i+2,j,k)-pf(i-2,j,k) ) )*idx(i)*cosphi(i,j)
          dr = ( a5(1)*( rf(i+1,j,k)-rf(i-1,j,k) ) &
               + a5(2)*( rf(i+2,j,k)-rf(i-2,j,k) ) )*idx(i)*cosphi(i,j)

          du = du + ( a04(1)*uf(i,j  ,k)+a04(2)*uf(i,j+1,k) &
                    + a04(3)*uf(i,j+2,k)+a04(4)*uf(i,j+3,k) &
                    + a04(5)*uf(i,j+4,k) )*idy(j)*sinphi(i,j)
          dv = dv + ( a04(1)*vf(i,j  ,k)+a04(2)*vf(i,j+1,k) &
                    + a04(3)*vf(i,j+2,k)+a04(4)*vf(i,j+3,k) &
                    + a04(5)*vf(i,j+4,k) )*idy(j)*sinphi(i,j)
          dw = dw + ( a04(1)*wf(i,j  ,k)+a04(2)*wf(i,j+1,k) &
                    + a04(3)*wf(i,j+2,k)+a04(4)*wf(i,j+3,k) &
                    + a04(5)*wf(i,j+4,k) )*idy(j)*sinphi(i,j)
          dp = dp + ( a04(1)*pf(i,j  ,k)+a04(2)*pf(i,j+1,k) &
                    + a04(3)*pf(i,j+2,k)+a04(4)*pf(i,j+3,k) &
                    + a04(5)*pf(i,j+4,k) )*idy(j)*sinphi(i,j)
          dr = dr + ( a04(1)*rf(i,j  ,k)+a04(2)*rf(i,j+1,k) &
                    + a04(3)*rf(i,j+2,k)+a04(4)*rf(i,j+3,k) &
                    + a04(5)*rf(i,j+4,k) )*idy(j)*sinphi(i,j)

          pt(i,k) = vg(i,k) * (dp+pf(i,j,k)*ir(i,j))
          ut(i,k) = vg(i,k) * (du+uf(i,j,k)*ir(i,j))
          vt(i,k) = vg(i,k) * (dv+vf(i,j,k)*ir(i,j))
          wt(i,k) = vg(i,k) * (dw+wf(i,j,k)*ir(i,j))
          rt(i,k) = vg(i,k) * (dr+rf(i,j,k)*ir(i,j))
       enddo
    enddo

    ! Boundary condition at imin-jmin (edge 1,1,1 /left-bottom)
    ! =========================================================
    if ((coord(1)==0).and.((BC_face(1,1)%sort==-11).or.(BC_face(1,1)%sort==0))) then
       ! Pointers for polar coordinates
       ! ==============================
       ir=>BC_edge(1,1,1)%ir
       cosphi=>BC_edge(1,1,1)%cosphi
       sinphi=>BC_edge(1,1,1)%sinphi

       ! Compute fluctuations
       ! ====================
       ! already done before

       i=2;j=1

       ! Compute group velocity vg
       ! =========================
       ! vg=u0*cos(phi)+v0*sin(phi) + sqrt(c0^2-w0^2-(u0*sin(phi)+v0*cos(phi))^2)
       do k=1,nz
          vg(i,k)= BC_face(2,1)%U0(i,j,k,2)*cosphi(i,j)+BC_face(2,1)%U0(i,j,k,3)*sinphi(i,j) &
               + sqrt(BC_face(2,1)%U0(i,j,k,6)-BC_face(2,1)%U0(i,j,k,4)**2 &
               -(BC_face(2,1)%U0(i,j,k,2)*sinphi(i,j)-BC_face(2,1)%U0(i,j,k,3)*cosphi(i,j))**2)
       enddo

       ! Non-centered derivatives
       ! ========================
       ! (Tam & Webb DRP schemes)
       do k=1,nz
          ! Non-centered derivatives in x-direction *cos(phi)
          ! =======================================
          du = ( a13(1)*uf(i-1,j,k)+a13(2)*uf(i,j,k) &
               + a13(3)*uf(i+1,j,k)+a13(4)*uf(i+2,j,k) &
               + a13(5)*uf(i+3,j,k) )*idx(i)*cosphi(i,j)
          dv = ( a13(1)*vf(i-1,j,k)+a13(2)*vf(i,j,k) &
               + a13(3)*vf(i+1,j,k)+a13(4)*vf(i+2,j,k) &
               + a13(5)*vf(i+3,j,k) )*idx(i)*cosphi(i,j)
          dw = ( a13(1)*wf(i-1,j,k)+a13(2)*wf(i,j,k) &
               + a13(3)*wf(i+1,j,k)+a13(4)*wf(i+2,j,k) &
               + a13(5)*wf(i+3,j,k) )*idx(i)*cosphi(i,j)
          dp = ( a13(1)*pf(i-1,j,k)+a13(2)*pf(i,j,k) &
               + a13(3)*pf(i+1,j,k)+a13(4)*pf(i+2,j,k) &
               + a13(5)*pf(i+3,j,k) )*idx(i)*cosphi(i,j)
          dr = ( a13(1)*rf(i-1,j,k)+a13(2)*rf(i,j,k) &
               + a13(3)*rf(i+1,j,k)+a13(4)*rf(i+2,j,k) &
               + a13(5)*rf(i+3,j,k) )*idx(i)*cosphi(i,j)

          ! Non-centered derivatives in y-direction *sin(phi)
          ! =======================================
          du = du + ( a04(1)*uf(i,j  ,k)+a04(2)*uf(i,j+1,k) &
                    + a04(3)*uf(i,j+2,k)+a04(4)*uf(i,j+3,k) &
                    + a04(5)*uf(i,j+4,k) )*idy(j)*sinphi(i,j)
          dv = dv + ( a04(1)*vf(i,j  ,k)+a04(2)*vf(i,j+1,k) &
                    + a04(3)*vf(i,j+2,k)+a04(4)*vf(i,j+3,k) &
                    + a04(5)*vf(i,j+4,k) )*idy(j)*sinphi(i,j)
          dw = dw + ( a04(1)*wf(i,j  ,k)+a04(2)*wf(i,j+1,k) &
                    + a04(3)*wf(i,j+2,k)+a04(4)*wf(i,j+3,k) &
                    + a04(5)*wf(i,j+4,k) )*idy(j)*sinphi(i,j)
          dp = dp + ( a04(1)*pf(i,j  ,k)+a04(2)*pf(i,j+1,k) &
                    + a04(3)*pf(i,j+2,k)+a04(4)*pf(i,j+3,k) &
                    + a04(5)*pf(i,j+4,k) )*idy(j)*sinphi(i,j)
          dr = dr + ( a04(1)*rf(i,j  ,k)+a04(2)*rf(i,j+1,k) &
                    + a04(3)*rf(i,j+2,k)+a04(4)*rf(i,j+3,k) &
                    + a04(5)*rf(i,j+4,k) )*idy(j)*sinphi(i,j)

          pt(i,k) = vg(i,k)*(dp + pf(i,j,k)*ir(i,j))
          ut(i,k) = vg(i,k)*(du + uf(i,j,k)*ir(i,j))
          vt(i,k) = vg(i,k)*(dv + vf(i,j,k)*ir(i,j))
          wt(i,k) = vg(i,k)*(dw + wf(i,j,k)*ir(i,j))
          rt(i,k) = vg(i,k)*(dr + rf(i,j,k)*ir(i,j))
       enddo
    endif


    ! Boundary condition at imax-jmin (edge 1,2,1 /right-bottom)
    ! =========================================================
    if ((coord(2)==0).and.((BC_face(1,2)%sort==-11).or.(BC_face(1,2)%sort==0))) then
       ! Pointers for polar coordinates
       ! ==============================
       ir=>BC_edge(1,2,1)%ir
       cosphi=>BC_edge(1,2,1)%cosphi
       sinphi=>BC_edge(1,2,1)%sinphi

       ! Compute fluctuations
       ! ====================
       ! already done before

       i=nx-1;j=1

       ! Compute group velocity vg
       ! =========================
       ! vg=u0*cos(phi)+v0*sin(phi) + sqrt(c0^2-w0^2-(u0*sin(phi)+v0*cos(phi))^2)
       do k=1,nz
          vg(i,k)= BC_face(2,1)%U0(i,j,k,2)*cosphi(i-nx+2,j)+BC_face(2,1)%U0(i,j,k,3)*sinphi(i-nx+2,j) &
               + sqrt(BC_face(2,1)%U0(i,j,k,6)-BC_face(2,1)%U0(i,j,k,4)**2 &
               -(BC_face(2,1)%U0(i,j,k,2)*sinphi(i-nx+2,j)-BC_face(2,1)%U0(i,j,k,3)*cosphi(i-nx+2,j))**2)
       enddo

       ! Non-centered derivatives
       ! ========================
       ! (Tam & Webb DRP schemes)
       do k=1,nz
          ! Non-centered derivatives in x-direction *cos(phi)
          ! =======================================
          du = ( a31(1)*uf(i+1,j,k)+a31(2)*uf(i,j,k) &
               + a31(3)*uf(i-1,j,k)+a31(4)*uf(i-2,j,k) &
               + a31(5)*uf(i-3,j,k) )*idx(i)*cosphi(i-nx+2,j)
          dv = ( a31(1)*vf(i+1,j,k)+a31(2)*vf(i,j,k) &
               + a31(3)*vf(i-1,j,k)+a31(4)*vf(i-2,j,k) &
               + a31(5)*vf(i-3,j,k) )*idx(i)*cosphi(i-nx+2,j)
          dw = ( a31(1)*wf(i+1,j,k)+a31(2)*wf(i,j,k) &
               + a31(3)*wf(i-1,j,k)+a31(4)*wf(i-2,j,k) &
               + a31(5)*wf(i-3,j,k) )*idx(i)*cosphi(i-nx+2,j)
          dp = ( a31(1)*pf(i+1,j,k)+a31(2)*pf(i,j,k) &
               + a31(3)*pf(i-1,j,k)+a31(4)*pf(i-2,j,k) &
               + a31(5)*pf(i-3,j,k) )*idx(i)*cosphi(i-nx+2,j)
          dr = ( a31(1)*rf(i+1,j,k)+a31(2)*rf(i,j,k) &
               + a31(3)*rf(i-1,j,k)+a31(4)*rf(i-2,j,k) &
               + a31(5)*rf(i-3,j,k) )*idx(i)*cosphi(i-nx+2,j)

          ! Non-centered derivatives in y-direction *sin(phi)
          ! =======================================
          du = du + ( a04(1)*uf(i,j  ,k)+a04(2)*uf(i,j+1,k) &
                    + a04(3)*uf(i,j+2,k)+a04(4)*uf(i,j+3,k) &
                    + a04(5)*uf(i,j+4,k) )*idy(j)*sinphi(i-nx+2,j)
          dv = dv + ( a04(1)*vf(i,j  ,k)+a04(2)*vf(i,j+1,k) &
                    + a04(3)*vf(i,j+2,k)+a04(4)*vf(i,j+3,k) &
                    + a04(5)*vf(i,j+4,k) )*idy(j)*sinphi(i-nx+2,j)
          dw = dw + ( a04(1)*wf(i,j  ,k)+a04(2)*wf(i,j+1,k) &
                    + a04(3)*wf(i,j+2,k)+a04(4)*wf(i,j+3,k) &
                    + a04(5)*wf(i,j+4,k) )*idy(j)*sinphi(i-nx+2,j)
          dp = dp + ( a04(1)*pf(i,j  ,k)+a04(2)*pf(i,j+1,k) &
                    + a04(3)*pf(i,j+2,k)+a04(4)*pf(i,j+3,k) &
                    + a04(5)*pf(i,j+4,k) )*idy(j)*sinphi(i-nx+2,j)
          dr = dr + ( a04(1)*rf(i,j  ,k)+a04(2)*rf(i,j+1,k) &
                    + a04(3)*rf(i,j+2,k)+a04(4)*rf(i,j+3,k) &
                    + a04(5)*rf(i,j+4,k) )*idy(j)*sinphi(i-nx+2,j)

          pt(i,k) = vg(i,k)*(dp + pf(i,j,k)*ir(i-nx+2,j))
          ut(i,k) = vg(i,k)*(du + uf(i,j,k)*ir(i-nx+2,j))
          vt(i,k) = vg(i,k)*(dv + vf(i,j,k)*ir(i-nx+2,j))
          wt(i,k) = vg(i,k)*(dw + wf(i,j,k)*ir(i-nx+2,j))
          rt(i,k) = vg(i,k)*(dr + rf(i,j,k)*ir(i-nx+2,j))
       enddo
    endif

    ! Update fluxes at each RK step
    ! =============================
    j=1
    do k=1,nz
       do i=ndx_td1m1,nfx_td1p1
          cp=cpcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
          av=avcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
          Krho(i,j,k)  = rt(i,k)
          Krhou(i,j,k) = uu(i,j,k)*rt(i,k)+rho_n(i,j,k)*ut(i,k)
          Krhov(i,j,k) = vv(i,j,k)*rt(i,k)+rho_n(i,j,k)*vt(i,k)
          Krhow(i,j,k) = ww(i,j,k)*rt(i,k)+rho_n(i,j,k)*wt(i,k)
          Krhoe(i,j,k) = cp/av*(pt(i,k)/c2_(i,k)-rt(i,k)) &
               + (rhoe_n(i,j,k)+prs(i,j,k))/rho_n(i,j,k)*rt(i,k) &
               + rho_n(i,j,k)*(uu(i,j,k)*ut(i,k)+vv(i,j,k)*vt(i,k)+ww(i,j,k)*wt(i,k))
       enddo
    enddo

  end subroutine bc_TD2d_1pt_jmin

  !=========================================================================================
  module subroutine bc_TD2d_1pt_jmax
  !=========================================================================================
    !> 2D Tam & Dong's BC on 1 point: boundary condition at jmax (top) - Cartesian version -
  !=========================================================================================
    implicit none
    !-------------------------------------------------------------------------
    integer :: i,j,k,l
    real(wp) :: cp,av,inn1,ntm1
    real(wp) :: du,dv,dw,dp,dr,ideltay
    real(wp), dimension(nx1:nx2,-3:1,nz) :: rf,uf,vf,wf,pf
    real(wp), dimension(nx,nz) :: rt,ut,vt,wt,pt,vg,c2_
    !-------------------------------------------------------------------------

    ! Index of top boundary
    ! ======================
    j=ny

    ! Pointers for polar coordinates
    ! ==============================
    ir=>BC_face(2,2)%ir
    cosphi=>BC_face(2,2)%cosphi
    sinphi=>BC_face(2,2)%sinphi

    ! Wall condition applied (different from Tam&Dong 5 points)
    ! ======================
    ! Wall BC at imin
    ! ---------------
    if (is_bc_wall(1,1)) then
        Krho(1,j,ndz_e:nfz_e)=0.0_wp
       Krhou(1,j,ndz_e:nfz_e)=0.0_wp
       Krhov(1,j,ndz_e:nfz_e)=0.0_wp
       Krhow(1,j,ndz_e:nfz_e)=0.0_wp
       Krhoe(1,j,ndz_e:nfz_e)=0.0_wp
    endif
    ! Wall BC at imax
    ! ---------------
    if (is_bc_wall(1,2)) then
        Krho(nx,j,ndz_e:nfz_e)=0.0_wp
       Krhou(nx,j,ndz_e:nfz_e)=0.0_wp
       Krhov(nx,j,ndz_e:nfz_e)=0.0_wp
       Krhow(nx,j,ndz_e:nfz_e)=0.0_wp
       Krhoe(nx,j,ndz_e:nfz_e)=0.0_wp
    endif
    ! Wall BC at kmin
    ! ---------------
    if (is_bc_wall(3,1)) then
        Krho(ndx_e:nfx_e,j,1)=0.0_wp
       Krhou(ndx_e:nfx_e,j,1)=0.0_wp
       Krhov(ndx_e:nfx_e,j,1)=0.0_wp
       Krhow(ndx_e:nfx_e,j,1)=0.0_wp
       Krhoe(ndx_e:nfx_e,j,1)=0.0_wp
    endif
    ! Wall BC at kmax
    ! ---------------
    if (is_bc_wall(3,2)) then
        Krho(ndx_e:nfx_e,j,nz)=0.0_wp
       Krhou(ndx_e:nfx_e,j,nz)=0.0_wp
       Krhov(ndx_e:nfx_e,j,nz)=0.0_wp
       Krhow(ndx_e:nfx_e,j,nz)=0.0_wp
       Krhoe(ndx_e:nfx_e,j,nz)=0.0_wp
    endif

    ! Sound speed
    ! ===========
    j=ny
    do k=1,nz
       do i=1,nx
          c2_(i,k)=c2calc_tro(Tmp(i,j,k),rho_n(i,j,k))
       enddo
    enddo

    ! Compute online time-averaged primitive variables
    ! ================================================
    if (irk==nrk) then
       inn1 = 1.0_wp/dble(ntotal)
       ntm1=dble(ntotal-1)
       ! time-averaged primitive variables
       do k=1,nz
          BC_face(2,2)%U0(:,-3:1,k,1)=(ntm1*BC_face(2,2)%U0(:,-3:1,k,1)+rho_n(:,ny-4:ny,k))*inn1
          BC_face(2,2)%U0(:,-3:1,k,2)=(ntm1*BC_face(2,2)%U0(:,-3:1,k,2)+   uu(:,ny-4:ny,k))*inn1
          BC_face(2,2)%U0(:,-3:1,k,3)=(ntm1*BC_face(2,2)%U0(:,-3:1,k,3)+   vv(:,ny-4:ny,k))*inn1
          BC_face(2,2)%U0(:,-3:1,k,4)=(ntm1*BC_face(2,2)%U0(:,-3:1,k,4)+   ww(:,ny-4:ny,k))*inn1
          BC_face(2,2)%U0(:,-3:1,k,5)=(ntm1*BC_face(2,2)%U0(:,-3:1,k,5)+  prs(:,ny-4:ny,k))*inn1
       enddo
       ! time-averaged sound speed squared
       do k=1,nz
          BC_face(2,2)%U0(1:nx,1,k,6)=(ntm1*BC_face(2,2)%U0(1:nx,1,k,6)+c2_(1:nx,k))*inn1
       enddo
    endif

    ! Compute fluctuations
    ! ====================
    do j=ny-4,ny
       l=j-(ny-1)
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
    j=ny
    l=j-(ny-1)
    do i=ndx_td1,nfx_td1
       do k=1,nz
          vg(i,k)= BC_face(2,2)%U0(i,l,k,2)*cosphi(i,1)+BC_face(2,2)%U0(i,l,k,3)*sinphi(i,1) &
               + sqrt(BC_face(2,2)%U0(i,l,k,6)-BC_face(2,2)%U0(i,l,k,4)**2 &
               -(BC_face(2,2)%U0(i,l,k,2)*sinphi(i,1)-BC_face(2,2)%U0(i,l,k,3)*cosphi(i,1))**2)
       enddo
    enddo

    ! Compute vg*[cos(phi)*dq/dx+sin(phi)*dq/dy+q/r]
    ! ==============================================
    ! (Tam & Webb DRP schemes)
    j=ny;l=j-(ny-1)
    ideltay = 1.0_wp/(y(j)-y(j-1))
    do i=ndx_td1,nfx_td1
       do k=1,nz
          ! du = ( a7(1)*( uf(i+1,l,k) - uf(i-1,l,k) ) + &
          !        a7(2)*( uf(i+2,l,k) - uf(i-2,l,k) ) + &
          !        a7(3)*( uf(i+3,l,k) - uf(i-3,l,k) ) ) *idx(i)*cosphi(i,1)
          ! dv = ( a7(1)*( vf(i+1,l,k) - vf(i-1,l,k) ) + &
          !        a7(2)*( vf(i+2,l,k) - vf(i-2,l,k) ) + &
          !        a7(3)*( vf(i+3,l,k) - vf(i-3,l,k) ) ) *idx(i)*cosphi(i,1)
          ! dw = ( a7(1)*( wf(i+1,l,k) - wf(i-1,l,k) ) + &
          !        a7(2)*( wf(i+2,l,k) - wf(i-2,l,k) ) + &
          !        a7(3)*( wf(i+3,l,k) - wf(i-3,l,k) ) ) *idx(i)*cosphi(i,1)
          ! dp = ( a7(1)*( pf(i+1,l,k) - pf(i-1,l,k) ) + &
          !        a7(2)*( pf(i+2,l,k) - pf(i-2,l,k) ) + &
          !        a7(3)*( pf(i+3,l,k) - pf(i-3,l,k) ) ) *idx(i)*cosphi(i,1)
          ! dr = ( a7(1)*( rf(i+1,l,k) - rf(i-1,l,k) ) + &
          !        a7(2)*( rf(i+2,l,k) - rf(i-2,l,k) ) + &
          !        a7(3)*( rf(i+3,l,k) - rf(i-3,l,k) ) ) *idx(i)*cosphi(i,1)

          ! du = du + ( a60(1)*uf(i,l  ,k)+a60(2)*uf(i,l-1,k) &
          !           + a60(3)*uf(i,l-2,k)+a60(4)*uf(i,l-3,k) &
          !           + a60(5)*uf(i,l-4,k)+a60(6)*uf(i,l-5,k) &
          !           + a60(7)*uf(i,l-6,k) ) *idy(j)*sinphi(i,1)
          ! dv = dv + ( a60(1)*vf(i,l  ,k)+a60(2)*vf(i,l-1,k) &
          !           + a60(3)*vf(i,l-2,k)+a60(4)*vf(i,l-3,k) &
          !           + a60(5)*vf(i,l-4,k)+a60(6)*vf(i,l-5,k) &
          !           + a60(7)*vf(i,l-6,k) ) *idy(j)*sinphi(i,1)
          ! dw = dw + ( a60(1)*wf(i,l  ,k)+a60(2)*wf(i,l-1,k) &
          !           + a60(3)*wf(i,l-2,k)+a60(4)*wf(i,l-3,k) &
          !           + a60(5)*wf(i,l-4,k)+a60(6)*wf(i,l-5,k) &
          !           + a60(7)*wf(i,l-6,k) ) *idy(j)*sinphi(i,1)
          ! dp = dp + ( a60(1)*pf(i,l  ,k)+a60(2)*pf(i,l-1,k) &
          !           + a60(3)*pf(i,l-2,k)+a60(4)*pf(i,l-3,k) &
          !           + a60(5)*pf(i,l-4,k)+a60(6)*pf(i,l-5,k) &
          !           + a60(7)*pf(i,l-6,k) ) *idy(j)*sinphi(i,1)
          ! dr = dr + ( a60(1)*rf(i,l  ,k)+a60(2)*rf(i,l-1,k) &
          !           + a60(3)*rf(i,l-2,k)+a60(4)*rf(i,l-3,k) &
          !           + a60(5)*rf(i,l-4,k)+a60(6)*rf(i,l-5,k) &
          !           + a60(7)*rf(i,l-6,k) ) *idy(j)*sinphi(i,1)


          ! ! Centered order 4 in x-direction
          ! du = ( a5(1)*( uf(i+1,l,k)-uf(i-1,l,k) ) &
          !      + a5(2)*( uf(i+2,l,k)-uf(i-2,l,k) ) )*idx(i)*cosphi(i,1)
          ! dv = ( a5(1)*( vf(i+1,l,k)-vf(i-1,l,k) ) &
          !      + a5(2)*( vf(i+2,l,k)-vf(i-2,l,k) ) )*idx(i)*cosphi(i,1)
          ! dw = ( a5(1)*( wf(i+1,l,k)-wf(i-1,l,k) ) &
          !      + a5(2)*( wf(i+2,l,k)-wf(i-2,l,k) ) )*idx(i)*cosphi(i,1)
          ! dp = ( a5(1)*( pf(i+1,l,k)-pf(i-1,l,k) ) &
          !      + a5(2)*( pf(i+2,l,k)-pf(i-2,l,k) ) )*idx(i)*cosphi(i,1)
          ! dr = ( a5(1)*( rf(i+1,l,k)-rf(i-1,l,k) ) &
          !      + a5(2)*( rf(i+2,l,k)-rf(i-2,l,k) ) )*idx(i)*cosphi(i,1)

          ! ! order 4 (decentered 4-0)
          ! du = du + ( a40(1)*uf(i,l,k)  +a40(2)*uf(i,l-1,k) &
          !           + a40(3)*uf(i,l-2,k)+a40(4)*uf(i,l-3,k) &
          !           + a40(5)*uf(i,l-4,k) )*idy(j)*sinphi(i,1)
          ! dv = dv + ( a40(1)*vf(i,l,k)  +a40(2)*vf(i,l-1,k) &
          !           + a40(3)*vf(i,l-2,k)+a40(4)*vf(i,l-3,k) &
          !           + a40(5)*vf(i,l-4,k) )*idy(j)*sinphi(i,1)
          ! dw = dw + ( a40(1)*wf(i,l,k)  +a40(2)*wf(i,l-1,k) &
          !           + a40(3)*wf(i,l-2,k)+a40(4)*wf(i,l-3,k) &
          !           + a40(5)*wf(i,l-4,k) )*idy(j)*sinphi(i,1)
          ! dp = dp + ( a40(1)*pf(i,l,k)  +a40(2)*pf(i,l-1,k) &
          !           + a40(3)*pf(i,l-2,k)+a40(4)*pf(i,l-3,k) &
          !           + a40(5)*pf(i,l-4,k) )*idy(j)*sinphi(i,1)
          ! dr = dr + ( a40(1)*rf(i,l,k)  +a40(2)*rf(i,l-1,k) &
          !           + a40(3)*rf(i,l-2,k)+a40(4)*rf(i,l-3,k) &
          !           + a40(5)*rf(i,l-4,k) )*idy(j)*sinphi(i,1)

          ! Centered order 4 in x-direction
          du = ( a5(1)*( uf(i+1,l,k)-uf(i-1,l,k) ) &
               + a5(2)*( uf(i+2,l,k)-uf(i-2,l,k) ) )*idx(i)*cosphi(i,1)
          dv = ( a5(1)*( vf(i+1,l,k)-vf(i-1,l,k) ) &
               + a5(2)*( vf(i+2,l,k)-vf(i-2,l,k) ) )*idx(i)*cosphi(i,1)
          dw = ( a5(1)*( wf(i+1,l,k)-wf(i-1,l,k) ) &
               + a5(2)*( wf(i+2,l,k)-wf(i-2,l,k) ) )*idx(i)*cosphi(i,1)
          dp = ( a5(1)*( pf(i+1,l,k)-pf(i-1,l,k) ) &
               + a5(2)*( pf(i+2,l,k)-pf(i-2,l,k) ) )*idx(i)*cosphi(i,1)
          dr = ( a5(1)*( rf(i+1,l,k)-rf(i-1,l,k) ) &
               + a5(2)*( rf(i+2,l,k)-rf(i-2,l,k) ) )*idx(i)*cosphi(i,1)

          ! Non-centered order 1 in y-direction (backward differences)
          du = du + (uf(i,l,k)-uf(i,l-1,k))*ideltay*sinphi(i,1)
          dv = dv + (vf(i,l,k)-vf(i,l-1,k))*ideltay*sinphi(i,1)
          dw = dw + (wf(i,l,k)-wf(i,l-1,k))*ideltay*sinphi(i,1)
          dp = dp + (pf(i,l,k)-pf(i,l-1,k))*ideltay*sinphi(i,1)
          dr = dr + (rf(i,l,k)-rf(i,l-1,k))*ideltay*sinphi(i,1)

          pt(i,k) = vg(i,k) * (dp+pf(i,l,k)*ir(i,1))
          ut(i,k) = vg(i,k) * (du+uf(i,l,k)*ir(i,1))
          vt(i,k) = vg(i,k) * (dv+vf(i,l,k)*ir(i,1))
          wt(i,k) = vg(i,k) * (dw+wf(i,l,k)*ir(i,1))
          rt(i,k) = vg(i,k) * (dr+rf(i,l,k)*ir(i,1))
       enddo
    enddo

    ! Boundary condition at imin-jmax (edge 1,1,2 /left-top)
    ! =========================================================
    if ((coord(1)==0).and.((BC_face(1,1)%sort==-11).or.(BC_face(1,1)%sort==0))) then
       ! Pointers for polar coordinates
       ! ==============================
       ir=>BC_edge(1,1,2)%ir
       cosphi=>BC_edge(1,1,2)%cosphi
       sinphi=>BC_edge(1,1,2)%sinphi

       ! Compute fluctuations
       ! ====================
       ! already done before

       i=2;j=ny;l=j-(ny-1)

       ! Compute group velocity vg
       ! =========================
       ! vg=u0*cos(phi)+v0*sin(phi) + sqrt(c0^2-w0^2-(u0*sin(phi)+v0*cos(phi))^2)
       do k=1,nz
          vg(i,k)= BC_face(2,2)%U0(i,l,k,2)*cosphi(i,2)+BC_face(2,2)%U0(i,l,k,3)*sinphi(i,2) &
               + sqrt(BC_face(2,2)%U0(i,l,k,6)-BC_face(2,2)%U0(i,l,k,4)**2 &
               -(BC_face(2,2)%U0(i,l,k,2)*sinphi(i,2)-BC_face(2,2)%U0(i,l,k,3)*cosphi(i,2))**2)
       enddo

       ideltay = 1.0_wp/(y(j)-y(j-1))
       do k=1,nz
          ! Non-centered derivatives (Tam & Webb DRP schemes)
          ! =================================================
          ! ! Non-centered derivatives in x-direction *cos(phi)
          ! ! =======================================
          ! du = ( a13(1)*uf(i-1,l,k)+a13(2)*uf(i  ,l,k) &
          !      + a13(3)*uf(i+1,l,k)+a13(4)*uf(i+2,l,k) &
          !      + a13(5)*uf(i+3,l,k) )*idx(i)*cosphi(i,2)
          ! dv = ( a13(1)*vf(i-1,l,k)+a13(2)*vf(i  ,l,k) &
          !      + a13(3)*vf(i+1,l,k)+a13(4)*vf(i+2,l,k) &
          !      + a13(5)*vf(i+3,l,k) )*idx(i)*cosphi(i,2)
          ! dw = ( a13(1)*wf(i-1,l,k)+a13(2)*wf(i  ,l,k) &
          !      + a13(3)*wf(i+1,l,k)+a13(4)*wf(i+2,l,k) &
          !      + a13(5)*wf(i+3,l,k) )*idx(i)*cosphi(i,2)
          ! dp = ( a13(1)*pf(i-1,l,k)+a13(2)*pf(i  ,l,k) &
          !      + a13(3)*pf(i+1,l,k)+a13(4)*pf(i+2,l,k) &
          !      + a13(5)*pf(i+3,l,k) )*idx(i)*cosphi(i,2)
          ! dr = ( a13(1)*rf(i-1,l,k)+a13(2)*rf(i  ,l,k) &
          !      + a13(3)*rf(i+1,l,k)+a13(4)*rf(i+2,l,k) &
          !      + a13(5)*rf(i+3,l,k) )*idx(i)*cosphi(i,2)

          ! ! Non-centered derivatives in y-direction *sin(phi)
          ! ! =======================================
          ! du = du + ( a40(1)*uf(i,l,k)  +a40(2)*uf(i,l-1,k) &
          !           + a40(3)*uf(i,l-2,k)+a40(4)*uf(i,l-3,k) &
          !           + a40(5)*uf(i,l-4,k) )*idy(j)*sinphi(i,2)
          ! dv = dv + ( a40(1)*vf(i,l,k)  +a40(2)*vf(i,l-1,k) &
          !           + a40(3)*vf(i,l-2,k)+a40(4)*vf(i,l-3,k) &
          !           + a40(5)*vf(i,l-4,k) )*idy(j)*sinphi(i,2)
          ! dw = dw + ( a40(1)*wf(i,l,k)  +a40(2)*wf(i,l-1,k) &
          !           + a40(3)*wf(i,l-2,k)+a40(4)*wf(i,l-3,k) &
          !           + a40(5)*wf(i,l-4,k) )*idy(j)*sinphi(i,2)
          ! dp = dp + ( a40(1)*pf(i,l,k)  +a40(2)*pf(i,l-1,k) &
          !           + a40(3)*pf(i,l-2,k)+a40(4)*pf(i,l-3,k) &
          !           + a40(5)*pf(i,l-4,k) )*idy(j)*sinphi(i,2)
          ! dr = dr + ( a40(1)*rf(i,l,k)  +a40(2)*rf(i,l-1,k) &
          !           + a40(3)*rf(i,l-2,k)+a40(4)*rf(i,l-3,k) &
          !           + a40(5)*rf(i,l-4,k) )*idy(j)*sinphi(i,2)

          ! Standard schemes on 3-point stencil
          du =  a3(1)*( uf(i+1,l,k)- uf(i-1,l,k))*idx2_imin*cosphi(i,2)
          dv =  a3(1)*( vf(i+1,l,k)- vf(i-1,l,k))*idx2_imin*cosphi(i,2)
          dw =  a3(1)*( wf(i+1,l,k)- wf(i-1,l,k))*idx2_imin*cosphi(i,2)
          dp =  a3(1)*( pf(i+1,l,k)- pf(i-1,l,k))*idx2_imin*cosphi(i,2)
          dr =  a3(1)*( rf(i+1,l,k)- rf(i-1,l,k))*idx2_imin*cosphi(i,2)

          ! Non-centered order 1 in y-direction (backward differences)
          du = du + (uf(i,l,k)-uf(i,l-1,k))*ideltay*sinphi(i,2)
          dv = dv + (vf(i,l,k)-vf(i,l-1,k))*ideltay*sinphi(i,2)
          dw = dw + (wf(i,l,k)-wf(i,l-1,k))*ideltay*sinphi(i,2)
          dp = dp + (pf(i,l,k)-pf(i,l-1,k))*ideltay*sinphi(i,2)
          dr = dr + (rf(i,l,k)-rf(i,l-1,k))*ideltay*sinphi(i,2)

          pt(i,k) = vg(i,k)*(dp + pf(i,l,k)*ir(i,2))
          ut(i,k) = vg(i,k)*(du + uf(i,l,k)*ir(i,2))
          vt(i,k) = vg(i,k)*(dv + vf(i,l,k)*ir(i,2))
          wt(i,k) = vg(i,k)*(dw + wf(i,l,k)*ir(i,2))
          rt(i,k) = vg(i,k)*(dr + rf(i,l,k)*ir(i,2))
       enddo
    endif

    ! Boundary condition at imax-jmax (edge 1,2,2 /right-top)
    ! =======================================================
    if ((coord(2)==ndomy-1).and.((BC_face(1,2)%sort==-11).or.(BC_face(1,2)%sort==0))) then
       ! Pointers for polar coordinates
       ! ==============================
       ir=>BC_edge(1,2,2)%ir
       cosphi=>BC_edge(1,2,2)%cosphi
       sinphi=>BC_edge(1,2,2)%sinphi

       ! Compute fluctuations
       ! ====================
       ! already done before

       i=nx-1;j=ny;l=j-(ny-1)

       ! Compute group velocity vg
       ! =========================
       ! vg=u0*cos(phi)+v0*sin(phi) + sqrt(c0^2-w0^2-(u0*sin(phi)+v0*cos(phi))^2)
       do k=1,nz
          vg(i,k)= BC_face(2,2)%U0(i,l,k,2)*cosphi(1,2)+BC_face(2,2)%U0(i,l,k,3)*sinphi(1,2) &
               + sqrt(BC_face(2,2)%U0(i,l,k,6)-BC_face(2,2)%U0(i,l,k,4)**2 &
               -(BC_face(2,2)%U0(i,l,k,2)*sinphi(1,2)-BC_face(2,2)%U0(i,l,k,3)*cosphi(1,2))**2)
       enddo

       i=nx-1;j=ny;l=j-(ny-1)
       ideltay = 1.0_wp/(y(j)-y(j-1))
       do k=1,nz
          ! Non-centered derivatives (Tam & Webb DRP schemes)
          ! =================================================
          ! ! Non-centered derivatives in x-direction *cos(phi)
          ! ! =======================================
          ! du = ( a31(1)*uf(i+1,l,k)+a31(2)*uf(i  ,l,k) &
          !      + a31(3)*uf(i-1,l,k)+a31(4)*uf(i-2,l,k) &
          !      + a31(5)*uf(i-3,l,k) )*idx(i)*cosphi(1,2)
          ! dv = ( a31(1)*vf(i+1,l,k)+a31(2)*vf(i  ,l,k) &
          !      + a31(3)*vf(i-1,l,k)+a31(4)*vf(i-2,l,k) &
          !      + a31(5)*vf(i-3,l,k) )*idx(i)*cosphi(1,2)
          ! dw = ( a31(1)*wf(i+1,l,k)+a31(2)*wf(i  ,l,k) &
          !      + a31(3)*wf(i-1,l,k)+a31(4)*wf(i-2,l,k) &
          !      + a31(5)*wf(i-3,l,k) )*idx(i)*cosphi(1,2)
          ! dp = ( a31(1)*pf(i+1,l,k)+a31(2)*pf(i  ,l,k) &
          !      + a31(3)*pf(i-1,l,k)+a31(4)*pf(i-2,l,k) &
          !      + a31(5)*pf(i-3,l,k) )*idx(i)*cosphi(1,2)
          ! dr = ( a31(1)*rf(i+1,l,k)+a31(2)*rf(i  ,l,k) &
          !      + a31(3)*rf(i-1,l,k)+a31(4)*rf(i-2,l,k) &
          !      + a31(5)*rf(i-3,l,k) )*idx(i)*cosphi(1,2)

          ! ! Non-centered derivatives in y-direction *sin(phi)
          ! ! =======================================
          ! du = du + ( a40(1)*uf(i,l,k)  +a40(2)*uf(i,l-1,k) &
          !           + a40(3)*uf(i,l-2,k)+a40(4)*uf(i,l-3,k) &
          !           + a40(5)*uf(i,l-4,k) )*idy(j)*sinphi(1,2)
          ! dv = dv + ( a40(1)*vf(i,l,k)  +a40(2)*vf(i,l-1,k) &
          !           + a40(3)*vf(i,l-2,k)+a40(4)*vf(i,l-3,k) &
          !           + a40(5)*vf(i,l-4,k) )*idy(j)*sinphi(1,2)
          ! dw = dw + ( a40(1)*wf(i,l,k)  +a40(2)*wf(i,l-1,k) &
          !           + a40(3)*wf(i,l-2,k)+a40(4)*wf(i,l-3,k) &
          !           + a40(5)*wf(i,l-4,k) )*idy(j)*sinphi(1,2)
          ! dp = dp + ( a40(1)*pf(i,l,k)  +a40(2)*pf(i,l-1,k) &
          !           + a40(3)*pf(i,l-2,k)+a40(4)*pf(i,l-3,k) &
          !           + a40(5)*pf(i,l-4,k) )*idy(j)*sinphi(1,2)
          ! dr = dr + ( a40(1)*rf(i,l,k)  +a40(2)*rf(i,l-1,k) &
          !           + a40(3)*rf(i,l-2,k)+a40(4)*rf(i,l-3,k) &
          !           + a40(5)*rf(i,l-4,k) )*idy(j)*sinphi(1,2)

          ! Standard schemes on 3-point stencil
          du =  a3(1)*( uf(i+1,l,k)- uf(i-1,l,k))*idx2_imin*cosphi(1,2)
          dv =  a3(1)*( vf(i+1,l,k)- vf(i-1,l,k))*idx2_imin*cosphi(1,2)
          dw =  a3(1)*( wf(i+1,l,k)- wf(i-1,l,k))*idx2_imin*cosphi(1,2)
          dp =  a3(1)*( pf(i+1,l,k)- pf(i-1,l,k))*idx2_imin*cosphi(1,2)
          dr =  a3(1)*( rf(i+1,l,k)- rf(i-1,l,k))*idx2_imin*cosphi(1,2)

          ! Non-centered order 1 in y-direction (backward differences)
          du = du + (uf(i,l,k)-uf(i,l-1,k))*ideltay*sinphi(1,2)
          dv = dv + (vf(i,l,k)-vf(i,l-1,k))*ideltay*sinphi(1,2)
          dw = dw + (wf(i,l,k)-wf(i,l-1,k))*ideltay*sinphi(1,2)
          dp = dp + (pf(i,l,k)-pf(i,l-1,k))*ideltay*sinphi(1,2)
          dr = dr + (rf(i,l,k)-rf(i,l-1,k))*ideltay*sinphi(1,2)

          pt(i,k) = vg(i,k)*(dp + pf(i,l,k)*ir(1,2))
          ut(i,k) = vg(i,k)*(du + uf(i,l,k)*ir(1,2))
          vt(i,k) = vg(i,k)*(dv + vf(i,l,k)*ir(1,2))
          wt(i,k) = vg(i,k)*(dw + wf(i,l,k)*ir(1,2))
          rt(i,k) = vg(i,k)*(dr + rf(i,l,k)*ir(1,2))
       enddo
    endif

    ! Update fluxes at each RK step
    ! =============================
    j=ny
    l=j-(ny-1)
    do k=1,nz
       do i=ndx_td1m1,nfx_td1p1
          cp=cpcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
          av=avcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
          Krho(i,j,k)  = rt(i,k)
          Krhou(i,j,k) = uu(i,j,k)*rt(i,k)+rho_n(i,j,k)*ut(i,k)
          Krhov(i,j,k) = vv(i,j,k)*rt(i,k)+rho_n(i,j,k)*vt(i,k)
          Krhow(i,j,k) = ww(i,j,k)*rt(i,k)+rho_n(i,j,k)*wt(i,k)
          Krhoe(i,j,k) = cp/av*(pt(i,k)/c2_(i,k)-rt(i,k)) &
               + (rhoe_n(i,j,k)+prs(i,j,k))/rho_n(i,j,k)*rt(i,k) &
               + rho_n(i,j,k)*(uu(i,j,k)*ut(i,k)+vv(i,j,k)*vt(i,k)+ww(i,j,k)*wt(i,k))
       enddo
    enddo

  end subroutine bc_TD2d_1pt_jmax

end submodule smod_TamDong2d_1pt
