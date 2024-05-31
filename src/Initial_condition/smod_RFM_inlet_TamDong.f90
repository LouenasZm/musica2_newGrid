!===============================================================================
submodule (mod_RFM) smod_RFM_TamDong
!===============================================================================
  !> author: AB & XG
  !> date: February 2024
  !> Generation of Random Fourier Modes (RFM)
  !> Injection in Tam & Dong's BCs [only for imin]
!=============================================================================== 

contains

  !===============================================================================
  module subroutine disturb_inlet_RFM_TamDong_imin(vg,ut_in,vt_in,wt_in)
  !===============================================================================
    !> Impose RFM disturbances at inlet (for Tam & Dong's BC)
  !===============================================================================
    use mod_time
    implicit none
    ! ----------------------------------------------------------------------------
    ! group velocity
    real(wp), dimension(ngh,ny,nz), intent(in) :: vg
    ! RFM disturbances
    real(wp), dimension(ngh,ny,nz), intent(out) :: ut_in,vt_in,wt_in
    ! ----------------------------------------------------------------------------
    ! local var
    integer :: i,j,k,nm
    real(wp) :: du_in,dv_in,dw_in
    real(wp) :: trk,kx,ky
    real(wp) :: phase,cphase,sphase
    real(wp) :: ut,vt,wt
    real(wp) :: dutdx,dvtdx,dwtdx,dutdy,dvtdy,dwtdy
    real(wp) :: xs,ys,zs ! scaled cooordinates
    ! scaled velocity components
    real(wp) :: us,vs,ws
    real(wp) :: dusdx,dvsdx,dwsdx
    real(wp) :: dusdy,dvsdy,dwsdy
    real(wp) :: dusdt,dvsdt,dwsdt
    ! ----------------------------------------------------------------------------
    ! TBL injection
    integer :: ndeb_j1,nend_j1,jbc
    integer, dimension(2) :: ndeb_j2,nend_j2

    trk = time + ck(irk)*deltat

    ! Determination of loop indices if TBL imin & imax
    if (is_RFM_FST) then
       if ((is_RFM_TBL(1)).and.(is_RFM_TBL(2))) then
          ndeb_j1=ndy_RFM-coord(2)*ny; nend_j1=nfy_RFM-coord(2)*ny ! Indices for FST
          ndeb_j2(1)=ndy; nend_j2(1)=ndeb_j1-1 ! Indices for TBL imin
          ndeb_j2(2)=nend_j1+1; nend_j2(2)=nfy  ! Indices for TBL imax
       ! Determination if TBL imin
       elseif (is_RFM_TBL(1)) then
          ndeb_j1=min(ndy_RFM-coord(2)*ny,nfy); nend_j1=nfy ! Indices for FST
          ndeb_j2(1)=ndy; nend_j2(1)=ndeb_j1-1 ! Indices for TBL imin
       ! Determination if TBL imax
       elseif (is_RFM_TBL(2)) then
          ndeb_j1=ndy; nend_j1=max(nfy_RFM-coord(2)*ny,ndy) ! Indices for FST
          ndeb_j2(2)=nend_j1+1; nend_j2(2)=nfy  ! Indices for TBL imax
       else
          ndeb_j1=ndy; nend_j1=nfy
       endif
    else
       ndeb_j2=ndy; nend_j2=nfy
    endif

    ut_in=0.0_wp; vt_in=0.0_wp; wt_in=0.0_wp

    ! -------------------
    ! I/ Injection of FST
    ! -------------------
    if (is_RFM_FST) then
       ! Construct turbulent stochastic field
       ! ------------------------------------
       if (anisotropy=='S') then
          if (is_curv3) then
             do i=1,ngh
                do j=ndeb_j1,nend_j1
                   do k=1,nz
                      ! scaled coordinates: x_scaled=R^T x
                      xs = (xc3(i,j,k)-u_m*trk)*vr(1,1,j) + (yc3(i,j,k)-v_m*trk)*vr(2,1,j) + (zc3(i,j,k)-w_m*trk)*vr(3,1,j)
                      ys = (xc3(i,j,k)-u_m*trk)*vr(1,2,j) + (yc3(i,j,k)-v_m*trk)*vr(2,2,j) + (zc3(i,j,k)-w_m*trk)*vr(3,2,j)
                      zs = (xc3(i,j,k)-u_m*trk)*vr(1,3,j) + (yc3(i,j,k)-v_m*trk)*vr(2,3,j) + (zc3(i,j,k)-w_m*trk)*vr(3,3,j)

                      ! summation over modes for scaled velocity components
                      us = 0.0_wp
                      vs = 0.0_wp
                      ws = 0.0_wp
                      dusdx = 0.0_wp
                      dvsdx = 0.0_wp
                      dwsdx = 0.0_wp
                      dusdy = 0.0_wp
                      dvsdy = 0.0_wp
                      dwsdy = 0.0_wp
                      dusdt = 0.0_wp
                      dvsdt = 0.0_wp
                      dwsdt = 0.0_wp
                      do nm=1,Nmode
                         phase = xk1s(nm,j)*xs+xk2s(nm,j)*ys+xk3s(nm,j)*zs + omn(nm)*trk + psi(nm)
                         cphase = cos(phase)
                         sphase = sin(phase)
                         us = us + cphase*sigma1s(nm,j)
                         vs = vs + cphase*sigma2s(nm,j)
                         ws = ws + cphase*sigma3s(nm,j)

                         kx = xk1s(nm,j)*vr(1,1,j) + xk2s(nm,j)*vr(1,2,j) + xk3s(nm,j)*vr(1,3,j)
                         dusdx = dusdx - kx*sphase*sigma1s(nm,j)
                         dvsdx = dvsdx - kx*sphase*sigma2s(nm,j)
                         dwsdx = dwsdx - kx*sphase*sigma3s(nm,j)

                         ky = xk1s(nm,j)*vr(2,1,j) + xk2s(nm,j)*vr(2,2,j) + xk3s(nm,j)*vr(2,3,j)
                         dusdy = dusdy-ky*sphase*sigma1s(nm,j)
                         dvsdy = dvsdy-ky*sphase*sigma2s(nm,j)
                         dwsdy = dwsdy-ky*sphase*sigma3s(nm,j)

                         dusdt = dusdt - omn(nm)*sphase*sigma1s(nm,j)
                         dvsdt = dvsdt - omn(nm)*sphase*sigma2s(nm,j)
                         dwsdt = dwsdt - omn(nm)*sphase*sigma3s(nm,j)
                      enddo

                      ! retrieve anisotropic field from scaled velocity field: u= R u_scaled
                      ut=us*vr(1,1,j)+vs*vr(1,2,j)+ws*vr(1,3,j)
                      vt=us*vr(2,1,j)+vs*vr(2,2,j)+ws*vr(2,3,j)
                      wt=us*vr(3,1,j)+vs*vr(3,2,j)+ws*vr(3,3,j)

                      dutdx=dusdx*vr(1,1,j)+dvsdx*vr(1,2,j)+dwsdx*vr(1,3,j)
                      dvtdx=dusdx*vr(2,1,j)+dvsdx*vr(2,2,j)+dwsdx*vr(2,3,j)
                      dwtdx=dusdx*vr(3,1,j)+dvsdx*vr(3,2,j)+dwsdx*vr(3,3,j)

                      dutdy=dusdy*vr(1,1,j)+dvsdy*vr(1,2,j)+dwsdy*vr(1,3,j)
                      dvtdy=dusdy*vr(2,1,j)+dvsdy*vr(2,2,j)+dwsdy*vr(2,3,j)
                      dwtdy=dusdy*vr(3,1,j)+dvsdy*vr(3,2,j)+dwsdy*vr(3,3,j)

                      du_in=ut*BC_face(1,1)%ir(i,j)+dutdx*BC_face(1,1)%cosphi(i,j)+dutdy*BC_face(1,1)%sinphi(i,j)
                      dv_in=vt*BC_face(1,1)%ir(i,j)+dvtdx*BC_face(1,1)%cosphi(i,j)+dvtdy*BC_face(1,1)%sinphi(i,j)
                      dw_in=wt*BC_face(1,1)%ir(i,j)+dwtdx*BC_face(1,1)%cosphi(i,j)+dwtdy*BC_face(1,1)%sinphi(i,j)

                      ut_in(i,j,k) = dusdt*vr(1,1,j) + dvsdt*vr(1,2,j) + dwsdt*vr(1,3,j) + vg(i,j,k)*du_in
                      vt_in(i,j,k) = dusdt*vr(2,1,j) + dvsdt*vr(2,2,j) + dwsdt*vr(2,3,j) + vg(i,j,k)*dv_in
                      wt_in(i,j,k) = dusdt*vr(3,1,j) + dvsdt*vr(3,2,j) + dwsdt*vr(3,3,j) + vg(i,j,k)*dw_in
                   enddo
                enddo
             enddo
          elseif (is_curv) then
             do i=1,ngh
                do j=ndeb_j1,nend_j1
                   do k=1,nz
                      ! scaled coordinates: x_scaled=R^T x
                      xs = (xc(i,j)-u_m*trk)*vr(1,1,j) + (yc(i,j)-v_m*trk)*vr(2,1,j) + (z(k)-w_m*trk)*vr(3,1,j)
                      ys = (xc(i,j)-u_m*trk)*vr(1,2,j) + (yc(i,j)-v_m*trk)*vr(2,2,j) + (z(k)-w_m*trk)*vr(3,2,j)
                      zs = (xc(i,j)-u_m*trk)*vr(1,3,j) + (yc(i,j)-v_m*trk)*vr(2,3,j) + (z(k)-w_m*trk)*vr(3,3,j)

                      ! summation over modes for scaled velocity components
                      us = 0.0_wp
                      vs = 0.0_wp
                      ws = 0.0_wp
                      dusdx = 0.0_wp
                      dvsdx = 0.0_wp
                      dwsdx = 0.0_wp
                      dusdy = 0.0_wp
                      dvsdy = 0.0_wp
                      dwsdy = 0.0_wp
                      dusdt = 0.0_wp
                      dvsdt = 0.0_wp
                      dwsdt = 0.0_wp
                      do nm=1,Nmode
                         phase = xk1s(nm,j)*xs+xk2s(nm,j)*ys+xk3s(nm,j)*zs + omn(nm)*trk + psi(nm)
                         cphase = cos(phase)
                         sphase = sin(phase)
                         us = us + cphase*sigma1s(nm,j)
                         vs = vs + cphase*sigma2s(nm,j)
                         ws = ws + cphase*sigma3s(nm,j)

                         kx = xk1s(nm,j)*vr(1,1,j) + xk2s(nm,j)*vr(1,2,j) + xk3s(nm,j)*vr(1,3,j)
                         dusdx = dusdx - kx*sphase*sigma1s(nm,j)
                         dvsdx = dvsdx - kx*sphase*sigma2s(nm,j)
                         dwsdx = dwsdx - kx*sphase*sigma3s(nm,j)

                         ky = xk1s(nm,j)*vr(2,1,j) + xk2s(nm,j)*vr(2,2,j) + xk3s(nm,j)*vr(2,3,j)
                         dusdy = dusdy-ky*sphase*sigma1s(nm,j)
                         dvsdy = dvsdy-ky*sphase*sigma2s(nm,j)
                         dwsdy = dwsdy-ky*sphase*sigma3s(nm,j)

                         dusdt = dusdt - omn(nm)*sphase*sigma1s(nm,j)
                         dvsdt = dvsdt - omn(nm)*sphase*sigma2s(nm,j)
                         dwsdt = dwsdt - omn(nm)*sphase*sigma3s(nm,j)
                      enddo

                      ! retrieve anisotropic field from scaled velocity field: u= R u_scaled
                      ut=us*vr(1,1,j)+vs*vr(1,2,j)+ws*vr(1,3,j)
                      vt=us*vr(2,1,j)+vs*vr(2,2,j)+ws*vr(2,3,j)
                      wt=us*vr(3,1,j)+vs*vr(3,2,j)+ws*vr(3,3,j)

                      dutdx=dusdx*vr(1,1,j)+dvsdx*vr(1,2,j)+dwsdx*vr(1,3,j)
                      dvtdx=dusdx*vr(2,1,j)+dvsdx*vr(2,2,j)+dwsdx*vr(2,3,j)
                      dwtdx=dusdx*vr(3,1,j)+dvsdx*vr(3,2,j)+dwsdx*vr(3,3,j)

                      dutdy=dusdy*vr(1,1,j)+dvsdy*vr(1,2,j)+dwsdy*vr(1,3,j)
                      dvtdy=dusdy*vr(2,1,j)+dvsdy*vr(2,2,j)+dwsdy*vr(2,3,j)
                      dwtdy=dusdy*vr(3,1,j)+dvsdy*vr(3,2,j)+dwsdy*vr(3,3,j)

                      du_in=ut*BC_face(1,1)%ir(i,j)+dutdx*BC_face(1,1)%cosphi(i,j)+dutdy*BC_face(1,1)%sinphi(i,j)
                      dv_in=vt*BC_face(1,1)%ir(i,j)+dvtdx*BC_face(1,1)%cosphi(i,j)+dvtdy*BC_face(1,1)%sinphi(i,j)
                      dw_in=wt*BC_face(1,1)%ir(i,j)+dwtdx*BC_face(1,1)%cosphi(i,j)+dwtdy*BC_face(1,1)%sinphi(i,j)

                      ut_in(i,j,k) = dusdt*vr(1,1,j) + dvsdt*vr(1,2,j) + dwsdt*vr(1,3,j) + vg(i,j,k)*du_in
                      vt_in(i,j,k) = dusdt*vr(2,1,j) + dvsdt*vr(2,2,j) + dwsdt*vr(2,3,j) + vg(i,j,k)*dv_in
                      wt_in(i,j,k) = dusdt*vr(3,1,j) + dvsdt*vr(3,2,j) + dwsdt*vr(3,3,j) + vg(i,j,k)*dw_in
                   enddo
                enddo
             enddo
          else
             do i=1,ngh
                do j=ndeb_j1,nend_j1
                   do k=1,nz
                      ! scaled coordinates: x_scaled=R^T x
                      xs = (x(i)-u_m*trk)*vr(1,1,j) + (y(j)-v_m*trk)*vr(2,1,j) + (z(k)-w_m*trk)*vr(3,1,j)
                      ys = (x(i)-u_m*trk)*vr(1,2,j) + (y(j)-v_m*trk)*vr(2,2,j) + (z(k)-w_m*trk)*vr(3,2,j)
                      zs = (x(i)-u_m*trk)*vr(1,3,j) + (y(j)-v_m*trk)*vr(2,3,j) + (z(k)-w_m*trk)*vr(3,3,j)

                      ! summation over modes for scaled velocity components
                      us = 0.0_wp
                      vs = 0.0_wp
                      ws = 0.0_wp
                      dusdx = 0.0_wp
                      dvsdx = 0.0_wp
                      dwsdx = 0.0_wp
                      dusdy = 0.0_wp
                      dvsdy = 0.0_wp
                      dwsdy = 0.0_wp
                      dusdt = 0.0_wp
                      dvsdt = 0.0_wp
                      dwsdt = 0.0_wp
                      do nm=1,Nmode
                         phase = xk1s(nm,j)*xs+xk2s(nm,j)*ys+xk3s(nm,j)*zs + omn(nm)*trk + psi(nm)
                         cphase = cos(phase)
                         sphase = sin(phase)
                         us = us + cphase*sigma1s(nm,j)
                         vs = vs + cphase*sigma2s(nm,j)
                         ws = ws + cphase*sigma3s(nm,j)

                         kx = xk1s(nm,j)*vr(1,1,j) + xk2s(nm,j)*vr(1,2,j) + xk3s(nm,j)*vr(1,3,j)
                         dusdx = dusdx - kx*sphase*sigma1s(nm,j)
                         dvsdx = dvsdx - kx*sphase*sigma2s(nm,j)
                         dwsdx = dwsdx - kx*sphase*sigma3s(nm,j)

                         ky = xk1s(nm,j)*vr(2,1,j) + xk2s(nm,j)*vr(2,2,j) + xk3s(nm,j)*vr(2,3,j)
                         dusdy = dusdy-ky*sphase*sigma1s(nm,j)
                         dvsdy = dvsdy-ky*sphase*sigma2s(nm,j)
                         dwsdy = dwsdy-ky*sphase*sigma3s(nm,j)

                         dusdt = dusdt - omn(nm)*sphase*sigma1s(nm,j)
                         dvsdt = dvsdt - omn(nm)*sphase*sigma2s(nm,j)
                         dwsdt = dwsdt - omn(nm)*sphase*sigma3s(nm,j)
                      enddo

                      ! retrieve anisotropic field from scaled velocity field: u= R u_scaled
                      ut=us*vr(1,1,j)+vs*vr(1,2,j)+ws*vr(1,3,j)
                      vt=us*vr(2,1,j)+vs*vr(2,2,j)+ws*vr(2,3,j)
                      wt=us*vr(3,1,j)+vs*vr(3,2,j)+ws*vr(3,3,j)

                      dutdx=dusdx*vr(1,1,j)+dvsdx*vr(1,2,j)+dwsdx*vr(1,3,j)
                      dvtdx=dusdx*vr(2,1,j)+dvsdx*vr(2,2,j)+dwsdx*vr(2,3,j)
                      dwtdx=dusdx*vr(3,1,j)+dvsdx*vr(3,2,j)+dwsdx*vr(3,3,j)

                      dutdy=dusdy*vr(1,1,j)+dvsdy*vr(1,2,j)+dwsdy*vr(1,3,j)
                      dvtdy=dusdy*vr(2,1,j)+dvsdy*vr(2,2,j)+dwsdy*vr(2,3,j)
                      dwtdy=dusdy*vr(3,1,j)+dvsdy*vr(3,2,j)+dwsdy*vr(3,3,j)

                      du_in=ut*BC_face(1,1)%ir(i,j)+dutdx*BC_face(1,1)%cosphi(i,j)+dutdy*BC_face(1,1)%sinphi(i,j)
                      dv_in=vt*BC_face(1,1)%ir(i,j)+dvtdx*BC_face(1,1)%cosphi(i,j)+dvtdy*BC_face(1,1)%sinphi(i,j)
                      dw_in=wt*BC_face(1,1)%ir(i,j)+dwtdx*BC_face(1,1)%cosphi(i,j)+dwtdy*BC_face(1,1)%sinphi(i,j)

                      ut_in(i,j,k) = dusdt*vr(1,1,j) + dvsdt*vr(1,2,j) + dwsdt*vr(1,3,j) + vg(i,j,k)*du_in
                      vt_in(i,j,k) = dusdt*vr(2,1,j) + dvsdt*vr(2,2,j) + dwsdt*vr(2,3,j) + vg(i,j,k)*dv_in
                      wt_in(i,j,k) = dusdt*vr(3,1,j) + dvsdt*vr(3,2,j) + dwsdt*vr(3,3,j) + vg(i,j,k)*dw_in
                   enddo
                enddo
             enddo
          endif
       else
          ut_in=0.0_wp
          vt_in=0.0_wp
          wt_in=0.0_wp
          if (is_curv3) then
             do i=1,ngh
                do j=ndeb_j1,nend_j1
                   do k=1,nz
                      ! summation over modes
                      ut=0.0_wp
                      vt=0.0_wp
                      wt=0.0_wp
                      dutdx=0.0_wp
                      dvtdx=0.0_wp
                      dwtdx=0.0_wp
                      dutdy=0.0_wp
                      dvtdy=0.0_wp
                      dwtdy=0.0_wp
                      do nm=1,Nmode
                         phase = xk1(nm)*(xc3(i,j,k)-u_m*trk) + xk2(nm)*(yc3(i,j,k)-v_m*trk) + xk3(nm)*(zc3(i,j,k)-w_m*trk) + omn(nm)*trk + psi(nm)
                         cphase=cos(phase)
                         sphase=sin(phase)

                         ut=ut+cphase*sigma1(nm)
                         vt=vt+cphase*sigma2(nm)
                         wt=wt+cphase*sigma3(nm)

                         dutdx=dutdx-xk1(nm)*sphase*sigma1(nm)
                         dvtdx=dvtdx-xk1(nm)*sphase*sigma2(nm)
                         dwtdx=dwtdx-xk1(nm)*sphase*sigma3(nm)

                         dutdy=dutdy-xk2(nm)*sphase*sigma1(nm)
                         dvtdy=dvtdy-xk2(nm)*sphase*sigma2(nm)
                         dwtdy=dwtdy-xk2(nm)*sphase*sigma3(nm)

                         ut_in(i,j,k) = ut_in(i,j,k) - (omn(nm) - xk1(nm)*u_m - xk2(nm)*v_m - xk3(nm)*w_m)*sphase*sigma1(nm)
                         vt_in(i,j,k) = vt_in(i,j,k) - (omn(nm) - xk1(nm)*u_m - xk2(nm)*v_m - xk3(nm)*w_m)*sphase*sigma2(nm)
                         wt_in(i,j,k) = wt_in(i,j,k) - (omn(nm) - xk1(nm)*u_m - xk2(nm)*v_m - xk3(nm)*w_m)*sphase*sigma3(nm)
                      enddo
                      du_in = ut*BC_face(1,1)%ir(i,j)+dutdx*BC_face(1,1)%cosphi(i,j)+dutdy*BC_face(1,1)%sinphi(i,j)
                      dv_in = vt*BC_face(1,1)%ir(i,j)+dvtdx*BC_face(1,1)%cosphi(i,j)+dvtdy*BC_face(1,1)%sinphi(i,j)
                      dw_in = wt*BC_face(1,1)%ir(i,j)+dwtdx*BC_face(1,1)%cosphi(i,j)+dwtdy*BC_face(1,1)%sinphi(i,j)

                      ut_in(i,j,k) = ut_in(i,j,k) + vg(i,j,k)*du_in
                      vt_in(i,j,k) = vt_in(i,j,k) + vg(i,j,k)*dv_in
                      wt_in(i,j,k) = wt_in(i,j,k) + vg(i,j,k)*dw_in
                   enddo
                enddo
             enddo
          else if (is_curv) then
             do i=1,ngh
                do j=ndeb_j1,nend_j1
                   do k=1,nz
                      ! summation over modes
                      ut=0.0_wp
                      vt=0.0_wp
                      wt=0.0_wp
                      dutdx=0.0_wp
                      dvtdx=0.0_wp
                      dwtdx=0.0_wp
                      dutdy=0.0_wp
                      dvtdy=0.0_wp
                      dwtdy=0.0_wp
                      do nm=1,Nmode
                         phase = xk1(nm)*(xc(i,j)-u_m*trk) + xk2(nm)*(yc(i,j)-v_m*trk) + xk3(nm)*(z(k)-w_m*trk) + omn(nm)*trk + psi(nm)
                         cphase=cos(phase)
                         sphase=sin(phase)

                         ut=ut+cphase*sigma1(nm)
                         vt=vt+cphase*sigma2(nm)
                         wt=wt+cphase*sigma3(nm)

                         dutdx=dutdx-xk1(nm)*sphase*sigma1(nm)
                         dvtdx=dvtdx-xk1(nm)*sphase*sigma2(nm)
                         dwtdx=dwtdx-xk1(nm)*sphase*sigma3(nm)

                         dutdy=dutdy-xk2(nm)*sphase*sigma1(nm)
                         dvtdy=dvtdy-xk2(nm)*sphase*sigma2(nm)
                         dwtdy=dwtdy-xk2(nm)*sphase*sigma3(nm)

                         ut_in(i,j,k) = ut_in(i,j,k) - (omn(nm) - xk1(nm)*u_m - xk2(nm)*v_m - xk3(nm)*w_m)*sphase*sigma1(nm)
                         vt_in(i,j,k) = vt_in(i,j,k) - (omn(nm) - xk1(nm)*u_m - xk2(nm)*v_m - xk3(nm)*w_m)*sphase*sigma2(nm)
                         wt_in(i,j,k) = wt_in(i,j,k) - (omn(nm) - xk1(nm)*u_m - xk2(nm)*v_m - xk3(nm)*w_m)*sphase*sigma3(nm)
                      enddo
                      du_in = ut*BC_face(1,1)%ir(i,j)+dutdx*BC_face(1,1)%cosphi(i,j)+dutdy*BC_face(1,1)%sinphi(i,j)
                      dv_in = vt*BC_face(1,1)%ir(i,j)+dvtdx*BC_face(1,1)%cosphi(i,j)+dvtdy*BC_face(1,1)%sinphi(i,j)
                      dw_in = wt*BC_face(1,1)%ir(i,j)+dwtdx*BC_face(1,1)%cosphi(i,j)+dwtdy*BC_face(1,1)%sinphi(i,j)

                      ut_in(i,j,k) = ut_in(i,j,k) + vg(i,j,k)*du_in
                      vt_in(i,j,k) = vt_in(i,j,k) + vg(i,j,k)*dv_in
                      wt_in(i,j,k) = wt_in(i,j,k) + vg(i,j,k)*dw_in
                   enddo
                enddo
             enddo
          else
             do i=1,ngh
                do j=ndeb_j1,nend_j1
                   do k=1,nz
                      ! summation over modes
                      ut=0.0_wp
                      vt=0.0_wp
                      wt=0.0_wp
                      dutdx=0.0_wp
                      dvtdx=0.0_wp
                      dwtdx=0.0_wp
                      dutdy=0.0_wp
                      dvtdy=0.0_wp
                      dwtdy=0.0_wp
                      do nm=1,Nmode
                         phase=xk1(nm)*(x(i)-u_m*trk) + xk2(nm)*(y(j)-v_m*trk) + xk3(nm)*(z(k)-w_m*trk) + omn(nm)*trk + psi(nm)
                         cphase=cos(phase)
                         sphase=sin(phase)

                         ut=ut+cphase*sigma1(nm)
                         vt=vt+cphase*sigma2(nm)
                         wt=wt+cphase*sigma3(nm)

                         dutdx=dutdx-xk1(nm)*sphase*sigma1(nm)
                         dvtdx=dvtdx-xk1(nm)*sphase*sigma2(nm)
                         dwtdx=dwtdx-xk1(nm)*sphase*sigma3(nm)

                         dutdy=dutdy-xk2(nm)*sphase*sigma1(nm)
                         dvtdy=dvtdy-xk2(nm)*sphase*sigma2(nm)
                         dwtdy=dwtdy-xk2(nm)*sphase*sigma3(nm)

                         ut_in(i,j,k) = ut_in(i,j,k) - (omn(nm)-xk1(nm)*u_m-xk2(nm)*v_m-xk3(nm)*w_m)*sphase*sigma1(nm)
                         vt_in(i,j,k) = vt_in(i,j,k) - (omn(nm)-xk1(nm)*u_m-xk2(nm)*v_m-xk3(nm)*w_m)*sphase*sigma2(nm)
                         wt_in(i,j,k) = wt_in(i,j,k) - (omn(nm)-xk1(nm)*u_m-xk2(nm)*v_m-xk3(nm)*w_m)*sphase*sigma3(nm)
                      enddo
                      du_in = ut*BC_face(1,1)%ir(i,j)+dutdx*BC_face(1,1)%cosphi(i,j)+dutdy*BC_face(1,1)%sinphi(i,j)
                      dv_in = vt*BC_face(1,1)%ir(i,j)+dvtdx*BC_face(1,1)%cosphi(i,j)+dvtdy*BC_face(1,1)%sinphi(i,j)
                      dw_in = wt*BC_face(1,1)%ir(i,j)+dwtdx*BC_face(1,1)%cosphi(i,j)+dwtdy*BC_face(1,1)%sinphi(i,j)

                      ut_in(i,j,k) = ut_in(i,j,k) + vg(i,j,k)*du_in
                      vt_in(i,j,k) = vt_in(i,j,k) + vg(i,j,k)*dv_in
                      wt_in(i,j,k) = wt_in(i,j,k) + vg(i,j,k)*dw_in
                   enddo
                enddo
             enddo
          endif
       endif

       ! Apply weighting function
       ! ========================
       if (anisotropy=='W') then
          do i=1,ngh
             do j=ndeb_j1,nend_j1
                do k=1,nz
                   ! multiply by rms profiles
                   ut_in(i,j,k) = ut_in(i,j,k)*u_rms(j)
                   vt_in(i,j,k) = vt_in(i,j,k)*v_rms(j)
                   wt_in(i,j,k) = wt_in(i,j,k)*w_rms(j)
                enddo
             enddo
          enddo
       endif

       ! Apply damping coefficients
       ! ==========================
       if (is_curv3) then
          do i=1,ngh
             do j=ndeb_j1,nend_j1
                do k=1,nz
                   ut_in(i,j,k) = ut_in(i,j,k)*damping_coeff3(j,k)
                   vt_in(i,j,k) = vt_in(i,j,k)*damping_coeff3(j,k)
                   wt_in(i,j,k) = wt_in(i,j,k)*damping_coeff3(j,k)
                enddo
             enddo
          enddo
       else
          do i=1,ngh
             do j=ndeb_j1,nend_j1
                do k=1,nz
                   ut_in(i,j,k) = ut_in(i,j,k)*damping_coeff(j)
                   vt_in(i,j,k) = vt_in(i,j,k)*damping_coeff(j)
                   wt_in(i,j,k) = wt_in(i,j,k)*damping_coeff(j)
                enddo
             enddo
          enddo
       endif
    endif

    ! --------------------
    ! II/ Injection of TBL
    ! --------------------
    do jbc=1,2
       if (is_RFM_TBL(jbc)) then
          if (is_curv3) then
             call mpistop("is_curv3 not implemented yet ",1)
          else if (is_curv) then
             call mpistop("is_curv not implemented yet ",1)
          else
             do i=1,ngh
                do j=ndeb_j2(jbc),nend_j2(jbc)
                   do k=1,nz
                      ! scaled coordinates: x_scaled=R^T x
                      xs = (x(i)-u_m_TBL*trk)*vr_TBL(1,1,j,jbc) + (y(j)-v_m_TBL*trk)*vr_TBL(2,1,j,jbc) + (z(k)-w_m_TBL*trk)*vr_TBL(3,1,j,jbc)
                      ys = (x(i)-u_m_TBL*trk)*vr_TBL(1,2,j,jbc) + (y(j)-v_m_TBL*trk)*vr_TBL(2,2,j,jbc) + (z(k)-w_m_TBL*trk)*vr_TBL(3,2,j,jbc)
                      zs = (x(i)-u_m_TBL*trk)*vr_TBL(1,3,j,jbc) + (y(j)-v_m_TBL*trk)*vr_TBL(2,3,j,jbc) + (z(k)-w_m_TBL*trk)*vr_TBL(3,3,j,jbc)

                      ! summation over modes for scaled velocity components
                      us = 0.0_wp
                      vs = 0.0_wp
                      ws = 0.0_wp
                      dusdx = 0.0_wp
                      dvsdx = 0.0_wp
                      dwsdx = 0.0_wp
                      dusdy = 0.0_wp
                      dvsdy = 0.0_wp
                      dwsdy = 0.0_wp
                      dusdt = 0.0_wp
                      dvsdt = 0.0_wp
                      dwsdt = 0.0_wp
                      do nm=1,Nmode
                         phase = xk1s_TBL(nm,j,jbc)*xs + xk2s_TBL(nm,j,jbc)*ys + xk3s_TBL(nm,j,jbc)*zs + psi_TBL(nm,jbc) ! + omn(nm)*trk
                         cphase = cos(phase)
                         sphase = sin(phase)
                         us = us + cphase*sigma1s_TBL(nm,j,jbc)
                         vs = vs + cphase*sigma2s_TBL(nm,j,jbc)
                         ws = ws + cphase*sigma3s_TBL(nm,j,jbc)

                         kx = xk1s_TBL(nm,j,jbc)*vr_TBL(1,1,j,jbc) + xk2s_TBL(nm,j,jbc)*vr_TBL(1,2,j,jbc) + xk3s_TBL(nm,j,jbc)*vr_TBL(1,3,j,jbc)
                         dusdx = dusdx - kx*sphase*sigma1s_TBL(nm,j,jbc)
                         dvsdx = dvsdx - kx*sphase*sigma2s_TBL(nm,j,jbc)
                         dwsdx = dwsdx - kx*sphase*sigma3s_TBL(nm,j,jbc)

                         ky = xk1s_TBL(nm,j,jbc)*vr_TBL(2,1,j,jbc) + xk2s_TBL(nm,j,jbc)*vr_TBL(2,2,j,jbc) + xk3s_TBL(nm,j,jbc)*vr_TBL(2,3,j,jbc)
                         dusdy = dusdy-ky*sphase*sigma1s_TBL(nm,j,jbc)
                         dvsdy = dvsdy-ky*sphase*sigma2s_TBL(nm,j,jbc)
                         dwsdy = dwsdy-ky*sphase*sigma3s_TBL(nm,j,jbc)

                         dusdt = dusdt !- omn(nm)*sphase*sigma1s_TBL(nm,j,jbc)
                         dvsdt = dvsdt !- omn(nm)*sphase*sigma2s_TBL(nm,j,jbc)
                         dwsdt = dwsdt !- omn(nm)*sphase*sigma3s_TBL(nm,j,jbc)
                      enddo

                      ! retrieve anisotropic field from scaled velocity field: u= R u_scaled
                      ut = us*vr_TBL(1,1,j,jbc) + vs*vr_TBL(1,2,j,jbc) + ws*vr_TBL(1,3,j,jbc)
                      vt = us*vr_TBL(2,1,j,jbc) + vs*vr_TBL(2,2,j,jbc) + ws*vr_TBL(2,3,j,jbc)
                      wt = us*vr_TBL(3,1,j,jbc) + vs*vr_TBL(3,2,j,jbc) + ws*vr_TBL(3,3,j,jbc)

                      dutdx = dusdx*vr_TBL(1,1,j,jbc) + dvsdx*vr_TBL(1,2,j,jbc) + dwsdx*vr_TBL(1,3,j,jbc)
                      dvtdx = dusdx*vr_TBL(2,1,j,jbc) + dvsdx*vr_TBL(2,2,j,jbc) + dwsdx*vr_TBL(2,3,j,jbc)
                      dwtdx = dusdx*vr_TBL(3,1,j,jbc) + dvsdx*vr_TBL(3,2,j,jbc) + dwsdx*vr_TBL(3,3,j,jbc)

                      dutdy = dusdy*vr_TBL(1,1,j,jbc) + dvsdy*vr_TBL(1,2,j,jbc) + dwsdy*vr_TBL(1,3,j,jbc)
                      dvtdy = dusdy*vr_TBL(2,1,j,jbc) + dvsdy*vr_TBL(2,2,j,jbc) + dwsdy*vr_TBL(2,3,j,jbc)
                      dwtdy = dusdy*vr_TBL(3,1,j,jbc) + dvsdy*vr_TBL(3,2,j,jbc) + dwsdy*vr_TBL(3,3,j,jbc)

                      du_in = ut*BC_face(1,1)%ir(i,j) + dutdx*BC_face(1,1)%cosphi(i,j) + dutdy*BC_face(1,1)%sinphi(i,j)
                      dv_in = vt*BC_face(1,1)%ir(i,j) + dvtdx*BC_face(1,1)%cosphi(i,j) + dvtdy*BC_face(1,1)%sinphi(i,j)
                      dw_in = wt*BC_face(1,1)%ir(i,j) + dwtdx*BC_face(1,1)%cosphi(i,j) + dwtdy*BC_face(1,1)%sinphi(i,j)

                      ut_in(i,j,k) = dusdt*vr_TBL(1,1,j,jbc) + dvsdt*vr_TBL(1,2,j,jbc) + dwsdt*vr_TBL(1,3,j,jbc) + vg(i,j,k)*du_in
                      vt_in(i,j,k) = dusdt*vr_TBL(2,1,j,jbc) + dvsdt*vr_TBL(2,2,j,jbc) + dwsdt*vr_TBL(2,3,j,jbc) + vg(i,j,k)*dv_in
                      wt_in(i,j,k) = dusdt*vr_TBL(3,1,j,jbc) + dvsdt*vr_TBL(3,2,j,jbc) + dwsdt*vr_TBL(3,3,j,jbc) + vg(i,j,k)*dw_in
                   enddo
                enddo
             enddo
          endif
       endif
    enddo

  end subroutine disturb_inlet_RFM_TamDong_imin

  !===============================================================================
  module subroutine disturb_inlet_RFM_TamDong_imin_jmin(vg,ut_in,vt_in,wt_in)
  !===============================================================================
    !> Impose RFM disturbances at inlet (for Tam & Dong's BC)
  !===============================================================================
    use mod_time
    implicit none
    ! ----------------------------------------------------------------------------
    ! group velocity
    real(wp), dimension(ngh,ngh,nz), intent(in) :: vg
    ! RFM disturbances
    real(wp), dimension(ngh,ngh,nz), intent(out) :: ut_in,vt_in,wt_in
    ! ----------------------------------------------------------------------------
    ! local var
    integer :: i,j,k,nm
    real(wp) :: du_in,dv_in,dw_in
    real(wp) :: trk,kx,ky,fy
    real(wp) :: phase,cphase,sphase
    real(wp) :: ut,vt,wt
    real(wp) :: dutdx,dvtdx,dwtdx,dutdy,dvtdy,dwtdy
    real(wp) :: xs,ys,zs ! scaled cooordinates
    ! scaled velocity components
    real(wp) :: us,vs,ws
    real(wp) :: dusdx,dvsdx,dwsdx
    real(wp) :: dusdy,dvsdy,dwsdy
    real(wp) :: dusdt,dvsdt,dwsdt
    ! ----------------------------------------------------------------------------

    trk = time + ck(irk)*deltat

    ut_in=0.0_wp; vt_in=0.0_wp; wt_in=0.0_wp

    ! -------------------
    ! I/ Injection of FST
    ! -------------------
    if (is_RFM_FST) then
       ! Construct turbulent stochastic field
       ! ------------------------------------
       if (anisotropy=='S') then
          if (is_curv3) then
             do i=1,ngh
                do j=1,ngh
                   do k=1,nz
                      ! scaled coordinates: x_scaled=R^T x
                      xs = (xc3(i,j,k)-u_m*trk)*vr(1,1,j) + (yc3(i,j,k)-v_m*trk)*vr(2,1,j) + (zc3(i,j,k)-w_m*trk)*vr(3,1,j)
                      ys = (xc3(i,j,k)-u_m*trk)*vr(1,2,j) + (yc3(i,j,k)-v_m*trk)*vr(2,2,j) + (zc3(i,j,k)-w_m*trk)*vr(3,2,j)
                      zs = (xc3(i,j,k)-u_m*trk)*vr(1,3,j) + (yc3(i,j,k)-v_m*trk)*vr(2,3,j) + (zc3(i,j,k)-w_m*trk)*vr(3,3,j)

                      ! summation over modes for scaled velocity components
                      us = 0.0_wp
                      vs = 0.0_wp
                      ws = 0.0_wp
                      dusdx = 0.0_wp
                      dvsdx = 0.0_wp
                      dwsdx = 0.0_wp
                      dusdy = 0.0_wp
                      dvsdy = 0.0_wp
                      dwsdy = 0.0_wp
                      dusdt = 0.0_wp
                      dvsdt = 0.0_wp
                      dwsdt = 0.0_wp
                      do nm=1,Nmode
                         phase = xk1s(nm,j)*xs+xk2s(nm,j)*ys+xk3s(nm,j)*zs + omn(nm)*trk + psi(nm)
                         cphase = cos(phase)
                         sphase = sin(phase)
                         us = us + cphase*sigma1s(nm,j)
                         vs = vs + cphase*sigma2s(nm,j)
                         ws = ws + cphase*sigma3s(nm,j)

                         kx = xk1s(nm,j)*vr(1,1,j) + xk2s(nm,j)*vr(1,2,j) + xk3s(nm,j)*vr(1,3,j)
                         dusdx = dusdx - kx*sphase*sigma1s(nm,j)
                         dvsdx = dvsdx - kx*sphase*sigma2s(nm,j)
                         dwsdx = dwsdx - kx*sphase*sigma3s(nm,j)

                         ky = xk1s(nm,j)*vr(2,1,j) + xk2s(nm,j)*vr(2,2,j) + xk3s(nm,j)*vr(2,3,j)
                         dusdy = dusdy-ky*sphase*sigma1s(nm,j)
                         dvsdy = dvsdy-ky*sphase*sigma2s(nm,j)
                         dwsdy = dwsdy-ky*sphase*sigma3s(nm,j)

                         dusdt = dusdt - omn(nm)*sphase*sigma1s(nm,j)
                         dvsdt = dvsdt - omn(nm)*sphase*sigma2s(nm,j)
                         dwsdt = dwsdt - omn(nm)*sphase*sigma3s(nm,j)
                      enddo

                      ! retrieve anisotropic field from scaled velocity field: u= R u_scaled
                      ut=us*vr(1,1,j)+vs*vr(1,2,j)+ws*vr(1,3,j)
                      vt=us*vr(2,1,j)+vs*vr(2,2,j)+ws*vr(2,3,j)
                      wt=us*vr(3,1,j)+vs*vr(3,2,j)+ws*vr(3,3,j)

                      dutdx=dusdx*vr(1,1,j)+dvsdx*vr(1,2,j)+dwsdx*vr(1,3,j)
                      dvtdx=dusdx*vr(2,1,j)+dvsdx*vr(2,2,j)+dwsdx*vr(2,3,j)
                      dwtdx=dusdx*vr(3,1,j)+dvsdx*vr(3,2,j)+dwsdx*vr(3,3,j)

                      dutdy=dusdy*vr(1,1,j)+dvsdy*vr(1,2,j)+dwsdy*vr(1,3,j)
                      dvtdy=dusdy*vr(2,1,j)+dvsdy*vr(2,2,j)+dwsdy*vr(2,3,j)
                      dwtdy=dusdy*vr(3,1,j)+dvsdy*vr(3,2,j)+dwsdy*vr(3,3,j)

                      du_in=ut*BC_edge(1,1,1)%ir(i,j)+dutdx*BC_edge(1,1,1)%cosphi(i,j)+dutdy*BC_edge(1,1,1)%sinphi(i,j)
                      dv_in=vt*BC_edge(1,1,1)%ir(i,j)+dvtdx*BC_edge(1,1,1)%cosphi(i,j)+dvtdy*BC_edge(1,1,1)%sinphi(i,j)
                      dw_in=wt*BC_edge(1,1,1)%ir(i,j)+dwtdx*BC_edge(1,1,1)%cosphi(i,j)+dwtdy*BC_edge(1,1,1)%sinphi(i,j)

                      ut_in(i,j,k) = dusdt*vr(1,1,j) + dvsdt*vr(1,2,j) + dwsdt*vr(1,3,j) + vg(i,j,k)*du_in
                      vt_in(i,j,k) = dusdt*vr(2,1,j) + dvsdt*vr(2,2,j) + dwsdt*vr(2,3,j) + vg(i,j,k)*dv_in
                      wt_in(i,j,k) = dusdt*vr(3,1,j) + dvsdt*vr(3,2,j) + dwsdt*vr(3,3,j) + vg(i,j,k)*dw_in
                   enddo
                enddo
             enddo
          else if (is_curv) then
             do i=1,ngh
                do j=1,ngh
                   do k=1,nz
                      ! scaled coordinates: x_scaled=R^T x
                      xs = (xc(i,j)-u_m*trk)*vr(1,1,j) + (yc(i,j)-v_m*trk)*vr(2,1,j) + (z(k)-w_m*trk)*vr(3,1,j)
                      ys = (xc(i,j)-u_m*trk)*vr(1,2,j) + (yc(i,j)-v_m*trk)*vr(2,2,j) + (z(k)-w_m*trk)*vr(3,2,j)
                      zs = (xc(i,j)-u_m*trk)*vr(1,3,j) + (yc(i,j)-v_m*trk)*vr(2,3,j) + (z(k)-w_m*trk)*vr(3,3,j)

                      ! summation over modes for scaled velocity components
                      us = 0.0_wp
                      vs = 0.0_wp
                      ws = 0.0_wp
                      dusdx = 0.0_wp
                      dvsdx = 0.0_wp
                      dwsdx = 0.0_wp
                      dusdy = 0.0_wp
                      dvsdy = 0.0_wp
                      dwsdy = 0.0_wp
                      dusdt = 0.0_wp
                      dvsdt = 0.0_wp
                      dwsdt = 0.0_wp
                      do nm=1,Nmode
                         phase = xk1s(nm,j)*xs+xk2s(nm,j)*ys+xk3s(nm,j)*zs + omn(nm)*trk + psi(nm)
                         cphase = cos(phase)
                         sphase = sin(phase)
                         us = us + cphase*sigma1s(nm,j)
                         vs = vs + cphase*sigma2s(nm,j)
                         ws = ws + cphase*sigma3s(nm,j)

                         kx = xk1s(nm,j)*vr(1,1,j) + xk2s(nm,j)*vr(1,2,j) + xk3s(nm,j)*vr(1,3,j)
                         dusdx = dusdx - kx*sphase*sigma1s(nm,j)
                         dvsdx = dvsdx - kx*sphase*sigma2s(nm,j)
                         dwsdx = dwsdx - kx*sphase*sigma3s(nm,j)

                         ky = xk1s(nm,j)*vr(2,1,j) + xk2s(nm,j)*vr(2,2,j) + xk3s(nm,j)*vr(2,3,j)
                         dusdy = dusdy-ky*sphase*sigma1s(nm,j)
                         dvsdy = dvsdy-ky*sphase*sigma2s(nm,j)
                         dwsdy = dwsdy-ky*sphase*sigma3s(nm,j)

                         dusdt = dusdt - omn(nm)*sphase*sigma1s(nm,j)
                         dvsdt = dvsdt - omn(nm)*sphase*sigma2s(nm,j)
                         dwsdt = dwsdt - omn(nm)*sphase*sigma3s(nm,j)
                      enddo

                      ! retrieve anisotropic field from scaled velocity field: u= R u_scaled
                      ut=us*vr(1,1,j)+vs*vr(1,2,j)+ws*vr(1,3,j)
                      vt=us*vr(2,1,j)+vs*vr(2,2,j)+ws*vr(2,3,j)
                      wt=us*vr(3,1,j)+vs*vr(3,2,j)+ws*vr(3,3,j)

                      dutdx=dusdx*vr(1,1,j)+dvsdx*vr(1,2,j)+dwsdx*vr(1,3,j)
                      dvtdx=dusdx*vr(2,1,j)+dvsdx*vr(2,2,j)+dwsdx*vr(2,3,j)
                      dwtdx=dusdx*vr(3,1,j)+dvsdx*vr(3,2,j)+dwsdx*vr(3,3,j)

                      dutdy=dusdy*vr(1,1,j)+dvsdy*vr(1,2,j)+dwsdy*vr(1,3,j)
                      dvtdy=dusdy*vr(2,1,j)+dvsdy*vr(2,2,j)+dwsdy*vr(2,3,j)
                      dwtdy=dusdy*vr(3,1,j)+dvsdy*vr(3,2,j)+dwsdy*vr(3,3,j)

                      du_in=ut*BC_edge(1,1,1)%ir(i,j)+dutdx*BC_edge(1,1,1)%cosphi(i,j)+dutdy*BC_edge(1,1,1)%sinphi(i,j)
                      dv_in=vt*BC_edge(1,1,1)%ir(i,j)+dvtdx*BC_edge(1,1,1)%cosphi(i,j)+dvtdy*BC_edge(1,1,1)%sinphi(i,j)
                      dw_in=wt*BC_edge(1,1,1)%ir(i,j)+dwtdx*BC_edge(1,1,1)%cosphi(i,j)+dwtdy*BC_edge(1,1,1)%sinphi(i,j)

                      ut_in(i,j,k) = dusdt*vr(1,1,j) + dvsdt*vr(1,2,j) + dwsdt*vr(1,3,j) + vg(i,j,k)*du_in
                      vt_in(i,j,k) = dusdt*vr(2,1,j) + dvsdt*vr(2,2,j) + dwsdt*vr(2,3,j) + vg(i,j,k)*dv_in
                      wt_in(i,j,k) = dusdt*vr(3,1,j) + dvsdt*vr(3,2,j) + dwsdt*vr(3,3,j) + vg(i,j,k)*dw_in
                   enddo
                enddo
             enddo
          else
             do i=1,ngh
                do j=1,ngh
                   do k=1,nz
                      ! scaled coordinates: x_scaled=R^T x
                      xs = (x(i)-u_m*trk)*vr(1,1,j) + (y(j)-v_m*trk)*vr(2,1,j) + (z(k)-w_m*trk)*vr(3,1,j)
                      ys = (x(i)-u_m*trk)*vr(1,2,j) + (y(j)-v_m*trk)*vr(2,2,j) + (z(k)-w_m*trk)*vr(3,2,j)
                      zs = (x(i)-u_m*trk)*vr(1,3,j) + (y(j)-v_m*trk)*vr(2,3,j) + (z(k)-w_m*trk)*vr(3,3,j)

                      ! summation over modes for scaled velocity components
                      us = 0.0_wp
                      vs = 0.0_wp
                      ws = 0.0_wp
                      dusdx = 0.0_wp
                      dvsdx = 0.0_wp
                      dwsdx = 0.0_wp
                      dusdy = 0.0_wp
                      dvsdy = 0.0_wp
                      dwsdy = 0.0_wp
                      dusdt = 0.0_wp
                      dvsdt = 0.0_wp
                      dwsdt = 0.0_wp
                      do nm=1,Nmode
                         phase = xk1s(nm,j)*xs+xk2s(nm,j)*ys+xk3s(nm,j)*zs + omn(nm)*trk + psi(nm)
                         cphase = cos(phase)
                         sphase = sin(phase)
                         us = us + cphase*sigma1s(nm,j)
                         vs = vs + cphase*sigma2s(nm,j)
                         ws = ws + cphase*sigma3s(nm,j)

                         kx = xk1s(nm,j)*vr(1,1,j) + xk2s(nm,j)*vr(1,2,j) + xk3s(nm,j)*vr(1,3,j)
                         dusdx = dusdx - kx*sphase*sigma1s(nm,j)
                         dvsdx = dvsdx - kx*sphase*sigma2s(nm,j)
                         dwsdx = dwsdx - kx*sphase*sigma3s(nm,j)

                         ky = xk1s(nm,j)*vr(2,1,j) + xk2s(nm,j)*vr(2,2,j) + xk3s(nm,j)*vr(2,3,j)
                         dusdy = dusdy-ky*sphase*sigma1s(nm,j)
                         dvsdy = dvsdy-ky*sphase*sigma2s(nm,j)
                         dwsdy = dwsdy-ky*sphase*sigma3s(nm,j)

                         dusdt = dusdt - omn(nm)*sphase*sigma1s(nm,j)
                         dvsdt = dvsdt - omn(nm)*sphase*sigma2s(nm,j)
                         dwsdt = dwsdt - omn(nm)*sphase*sigma3s(nm,j)
                      enddo

                      ! retrieve anisotropic field from scaled velocity field: u= R u_scaled
                      ut=us*vr(1,1,j)+vs*vr(1,2,j)+ws*vr(1,3,j)
                      vt=us*vr(2,1,j)+vs*vr(2,2,j)+ws*vr(2,3,j)
                      wt=us*vr(3,1,j)+vs*vr(3,2,j)+ws*vr(3,3,j)

                      dutdx=dusdx*vr(1,1,j)+dvsdx*vr(1,2,j)+dwsdx*vr(1,3,j)
                      dvtdx=dusdx*vr(2,1,j)+dvsdx*vr(2,2,j)+dwsdx*vr(2,3,j)
                      dwtdx=dusdx*vr(3,1,j)+dvsdx*vr(3,2,j)+dwsdx*vr(3,3,j)

                      dutdy=dusdy*vr(1,1,j)+dvsdy*vr(1,2,j)+dwsdy*vr(1,3,j)
                      dvtdy=dusdy*vr(2,1,j)+dvsdy*vr(2,2,j)+dwsdy*vr(2,3,j)
                      dwtdy=dusdy*vr(3,1,j)+dvsdy*vr(3,2,j)+dwsdy*vr(3,3,j)

                      du_in=ut*BC_edge(1,1,1)%ir(i,j)+dutdx*BC_edge(1,1,1)%cosphi(i,j)+dutdy*BC_edge(1,1,1)%sinphi(i,j)
                      dv_in=vt*BC_edge(1,1,1)%ir(i,j)+dvtdx*BC_edge(1,1,1)%cosphi(i,j)+dvtdy*BC_edge(1,1,1)%sinphi(i,j)
                      dw_in=wt*BC_edge(1,1,1)%ir(i,j)+dwtdx*BC_edge(1,1,1)%cosphi(i,j)+dwtdy*BC_edge(1,1,1)%sinphi(i,j)

                      ut_in(i,j,k) = dusdt*vr(1,1,j) + dvsdt*vr(1,2,j) + dwsdt*vr(1,3,j) + vg(i,j,k)*du_in
                      vt_in(i,j,k) = dusdt*vr(2,1,j) + dvsdt*vr(2,2,j) + dwsdt*vr(2,3,j) + vg(i,j,k)*dv_in
                      wt_in(i,j,k) = dusdt*vr(3,1,j) + dvsdt*vr(3,2,j) + dwsdt*vr(3,3,j) + vg(i,j,k)*dw_in
                   enddo
                enddo
             enddo
          endif
       else
          if (is_curv3) then
             do i=1,ngh
                do j=1,ngh
                   do k=1,nz
                      ! summation over modes
                      ut=0.0_wp
                      vt=0.0_wp
                      wt=0.0_wp
                      dutdx=0.0_wp
                      dvtdx=0.0_wp
                      dwtdx=0.0_wp
                      dutdy=0.0_wp
                      dvtdy=0.0_wp
                      dwtdy=0.0_wp
                      do nm=1,Nmode
                         phase = xk1(nm)*(xc3(i,j,k)-u_m*trk) + xk2(nm)*(yc3(i,j,k)-v_m*trk) + xk3(nm)*(zc3(i,j,k)-w_m*trk) + omn(nm)*trk + psi(nm)
                         cphase=cos(phase)
                         sphase=sin(phase)
                         ut=ut+cphase*sigma1(nm)
                         vt=vt+cphase*sigma2(nm)
                         wt=wt+cphase*sigma3(nm)

                         dutdx=dutdx-xk1(nm)*sphase*sigma1(nm)
                         dvtdx=dvtdx-xk1(nm)*sphase*sigma2(nm)
                         dwtdx=dwtdx-xk1(nm)*sphase*sigma3(nm)

                         dutdy=dutdy-xk2(nm)*sphase*sigma1(nm)
                         dvtdy=dvtdy-xk2(nm)*sphase*sigma2(nm)
                         dwtdy=dwtdy-xk2(nm)*sphase*sigma3(nm)

                         ut_in(i,j,k) = ut_in(i,j,k) - (omn(nm) - xk1(nm)*u_m - xk2(nm)*v_m - xk3(nm)*w_m)*sphase*sigma1(nm)
                         vt_in(i,j,k) = vt_in(i,j,k) - (omn(nm) - xk1(nm)*u_m - xk2(nm)*v_m - xk3(nm)*w_m)*sphase*sigma2(nm)
                         wt_in(i,j,k) = wt_in(i,j,k) - (omn(nm) - xk1(nm)*u_m - xk2(nm)*v_m - xk3(nm)*w_m)*sphase*sigma3(nm)
                      enddo
                      du_in = ut*BC_edge(1,1,1)%ir(i,j)+dutdx*BC_edge(1,1,1)%cosphi(i,j)+dutdy*BC_edge(1,1,1)%sinphi(i,j)
                      dv_in = vt*BC_edge(1,1,1)%ir(i,j)+dvtdx*BC_edge(1,1,1)%cosphi(i,j)+dvtdy*BC_edge(1,1,1)%sinphi(i,j)
                      dw_in = wt*BC_edge(1,1,1)%ir(i,j)+dwtdx*BC_edge(1,1,1)%cosphi(i,j)+dwtdy*BC_edge(1,1,1)%sinphi(i,j)

                      ut_in(i,j,k) = ut_in(i,j,k) + vg(i,j,k)*du_in
                      vt_in(i,j,k) = vt_in(i,j,k) + vg(i,j,k)*dv_in
                      wt_in(i,j,k) = wt_in(i,j,k) + vg(i,j,k)*dw_in
                   enddo
                enddo
             enddo
          else if (is_curv) then
             do i=1,ngh
                do j=1,ngh
                   do k=1,nz
                      ! summation over modes
                      ut=0.0_wp
                      vt=0.0_wp
                      wt=0.0_wp
                      dutdx=0.0_wp
                      dvtdx=0.0_wp
                      dwtdx=0.0_wp
                      dutdy=0.0_wp
                      dvtdy=0.0_wp
                      dwtdy=0.0_wp
                      do nm=1,Nmode
                         phase = xk1(nm)*(xc(i,j)-u_m*trk) + xk2(nm)*(yc(i,j)-v_m*trk) + xk3(nm)*(z(k)-w_m*trk) + omn(nm)*trk + psi(nm)
                         cphase=cos(phase)
                         sphase=sin(phase)
                         ut=ut+cphase*sigma1(nm)
                         vt=vt+cphase*sigma2(nm)
                         wt=wt+cphase*sigma3(nm)

                         dutdx=dutdx-xk1(nm)*sphase*sigma1(nm)
                         dvtdx=dvtdx-xk1(nm)*sphase*sigma2(nm)
                         dwtdx=dwtdx-xk1(nm)*sphase*sigma3(nm)

                         dutdy=dutdy-xk2(nm)*sphase*sigma1(nm)
                         dvtdy=dvtdy-xk2(nm)*sphase*sigma2(nm)
                         dwtdy=dwtdy-xk2(nm)*sphase*sigma3(nm)

                         ut_in(i,j,k) = ut_in(i,j,k) - (omn(nm) - xk1(nm)*u_m - xk2(nm)*v_m - xk3(nm)*w_m)*sphase*sigma1(nm)
                         vt_in(i,j,k) = vt_in(i,j,k) - (omn(nm) - xk1(nm)*u_m - xk2(nm)*v_m - xk3(nm)*w_m)*sphase*sigma2(nm)
                         wt_in(i,j,k) = wt_in(i,j,k) - (omn(nm) - xk1(nm)*u_m - xk2(nm)*v_m - xk3(nm)*w_m)*sphase*sigma3(nm)
                      enddo
                      du_in = ut*BC_edge(1,1,1)%ir(i,j)+dutdx*BC_edge(1,1,1)%cosphi(i,j)+dutdy*BC_edge(1,1,1)%sinphi(i,j)
                      dv_in = vt*BC_edge(1,1,1)%ir(i,j)+dvtdx*BC_edge(1,1,1)%cosphi(i,j)+dvtdy*BC_edge(1,1,1)%sinphi(i,j)
                      dw_in = wt*BC_edge(1,1,1)%ir(i,j)+dwtdx*BC_edge(1,1,1)%cosphi(i,j)+dwtdy*BC_edge(1,1,1)%sinphi(i,j)

                      ut_in(i,j,k) = ut_in(i,j,k) + vg(i,j,k)*du_in
                      vt_in(i,j,k) = vt_in(i,j,k) + vg(i,j,k)*dv_in
                      wt_in(i,j,k) = wt_in(i,j,k) + vg(i,j,k)*dw_in
                   enddo
                enddo
             enddo
          else
             do i=1,ngh
                do j=1,ngh
                   do k=1,nz
                      ! summation over modes
                      ut=0.0_wp
                      vt=0.0_wp
                      wt=0.0_wp
                      dutdx=0.0_wp
                      dvtdx=0.0_wp
                      dwtdx=0.0_wp
                      dutdy=0.0_wp
                      dvtdy=0.0_wp
                      dwtdy=0.0_wp
                      do nm=1,Nmode
                         phase=xk1(nm)*(x(i)-u_m*trk) + xk2(nm)*(y(j)-v_m*trk) + xk3(nm)*(z(k)-w_m*trk) + omn(nm)*trk + psi(nm)
                         cphase=cos(phase)
                         sphase=sin(phase)
                         ut=ut+cphase*sigma1(nm)
                         vt=vt+cphase*sigma2(nm)
                         wt=wt+cphase*sigma3(nm)

                         dutdx=dutdx-xk1(nm)*sphase*sigma1(nm)
                         dvtdx=dvtdx-xk1(nm)*sphase*sigma2(nm)
                         dwtdx=dwtdx-xk1(nm)*sphase*sigma3(nm)

                         dutdy=dutdy-xk2(nm)*sphase*sigma1(nm)
                         dvtdy=dvtdy-xk2(nm)*sphase*sigma2(nm)
                         dwtdy=dwtdy-xk2(nm)*sphase*sigma3(nm)

                         ut_in(i,j,k) = ut_in(i,j,k) - (omn(nm) - xk1(nm)*u_m - xk2(nm)*v_m - xk3(nm)*w_m)*sphase*sigma1(nm)
                         vt_in(i,j,k) = vt_in(i,j,k) - (omn(nm) - xk1(nm)*u_m - xk2(nm)*v_m - xk3(nm)*w_m)*sphase*sigma2(nm)
                         wt_in(i,j,k) = wt_in(i,j,k) - (omn(nm) - xk1(nm)*u_m - xk2(nm)*v_m - xk3(nm)*w_m)*sphase*sigma3(nm)
                      enddo
                      du_in = ut*BC_edge(1,1,1)%ir(i,j)+dutdx*BC_edge(1,1,1)%cosphi(i,j)+dutdy*BC_edge(1,1,1)%sinphi(i,j)
                      dv_in = vt*BC_edge(1,1,1)%ir(i,j)+dvtdx*BC_edge(1,1,1)%cosphi(i,j)+dvtdy*BC_edge(1,1,1)%sinphi(i,j)
                      dw_in = wt*BC_edge(1,1,1)%ir(i,j)+dwtdx*BC_edge(1,1,1)%cosphi(i,j)+dwtdy*BC_edge(1,1,1)%sinphi(i,j)

                      ut_in(i,j,k) = ut_in(i,j,k) + vg(i,j,k)*du_in
                      vt_in(i,j,k) = vt_in(i,j,k) + vg(i,j,k)*dv_in
                      wt_in(i,j,k) = wt_in(i,j,k) + vg(i,j,k)*dw_in
                   enddo
                enddo
             enddo
          endif
       endif

       ! Apply weighting function
       ! ========================
       if (anisotropy=='W') then
          do i=1,ngh
             do j=1,ngh
                fy=(y(j)-y(1))/(y(ngh+1)-y(1))
                do k=1,nz
                   ! add linear damping close to the wall
                   ut_in(i,j,k) = ut_in(i,j,k)*u_rms(j)*fy
                   vt_in(i,j,k) = vt_in(i,j,k)*v_rms(j)*fy
                   wt_in(i,j,k) = wt_in(i,j,k)*w_rms(j)*fy
                enddo
             enddo
          enddo
       endif

       ! Apply damping coefficients
       ! ==========================
       if (is_curv3) then
          do i=1,ngh
             do j=1,ngh
                do k=1,nz
                   ut_in(i,j,k) = ut_in(i,j,k)*damping_coeff3(j,k)
                   vt_in(i,j,k) = vt_in(i,j,k)*damping_coeff3(j,k)
                   wt_in(i,j,k) = wt_in(i,j,k)*damping_coeff3(j,k)
                enddo
             enddo
          enddo
       else
          do i=1,ngh
             do j=1,ngh
                do k=1,nz
                   ut_in(i,j,k) = ut_in(i,j,k)*damping_coeff(j)
                   vt_in(i,j,k) = vt_in(i,j,k)*damping_coeff(j)
                   wt_in(i,j,k) = wt_in(i,j,k)*damping_coeff(j)
                enddo
             enddo
          enddo
       endif
    endif

  end subroutine disturb_inlet_RFM_TamDong_imin_jmin

  !===============================================================================
  module subroutine disturb_inlet_RFM_TamDong_imin_jmax(vg,ut_in,vt_in,wt_in)
  !===============================================================================
    !> Impose RFM disturbances at inlet (for Tam & Dong's BC)
  !===============================================================================
    use mod_time
    implicit none
    ! ----------------------------------------------------------------------------
    ! group velocity
    real(wp), dimension(ngh,ngh,nz), intent(in) :: vg
    ! RFM disturbances
    real(wp), dimension(ngh,ngh,nz), intent(out) :: ut_in,vt_in,wt_in
    ! ----------------------------------------------------------------------------
    ! local var
    integer :: i,j,k,l,nm
    real(wp) :: du_in,dv_in,dw_in
    real(wp) :: trk,kx,ky,fy
    real(wp) :: phase,cphase,sphase
    real(wp) :: ut,vt,wt
    real(wp) :: dutdx,dvtdx,dwtdx,dutdy,dvtdy,dwtdy
    real(wp) :: xs,ys,zs ! scaled cooordinates
    ! scaled velocity components
    real(wp) :: us,vs,ws
    real(wp) :: dusdx,dvsdx,dwsdx
    real(wp) :: dusdy,dvsdy,dwsdy
    real(wp) :: dusdt,dvsdt,dwsdt
    ! ----------------------------------------------------------------------------

    trk = time + ck(irk)*deltat

    ut_in=0.0_wp; vt_in=0.0_wp; wt_in=0.0_wp

    ! -------------------
    ! I/ Injection of FST
    ! -------------------
    if (is_RFM_FST) then
       ! Construct turbulent stochastic field
       ! ------------------------------------
       if (anisotropy=='S') then
          if (is_curv3) then
             do i=1,ngh
                do j=ny-ngh+1,ny
                   l=j-(ny-ngh)
                   do k=1,nz
                      ! scaled coordinates: x_scaled=R^T x
                      xs = (xc3(i,j,k)-u_m*trk)*vr(1,1,j) + (yc3(i,j,k)-v_m*trk)*vr(2,1,j) + (zc3(i,j,k)-w_m*trk)*vr(3,1,j)
                      ys = (xc3(i,j,k)-u_m*trk)*vr(1,2,j) + (yc3(i,j,k)-v_m*trk)*vr(2,2,j) + (zc3(i,j,k)-w_m*trk)*vr(3,2,j)
                      zs = (xc3(i,j,k)-u_m*trk)*vr(1,3,j) + (yc3(i,j,k)-v_m*trk)*vr(2,3,j) + (zc3(i,j,k)-w_m*trk)*vr(3,3,j)

                      ! summation over modes for scaled velocity components
                      us = 0.0_wp
                      vs = 0.0_wp
                      ws = 0.0_wp
                      dusdx = 0.0_wp
                      dvsdx = 0.0_wp
                      dwsdx = 0.0_wp
                      dusdy = 0.0_wp
                      dvsdy = 0.0_wp
                      dwsdy = 0.0_wp
                      dusdt = 0.0_wp
                      dvsdt = 0.0_wp
                      dwsdt = 0.0_wp
                      do nm=1,Nmode
                         phase = xk1s(nm,j)*xs+xk2s(nm,j)*ys+xk3s(nm,j)*zs + omn(nm)*trk + psi(nm)
                         cphase = cos(phase)
                         sphase = sin(phase)
                         us = us + cphase*sigma1s(nm,j)
                         vs = vs + cphase*sigma2s(nm,j)
                         ws = ws + cphase*sigma3s(nm,j)

                         kx = xk1s(nm,j)*vr(1,1,j) + xk2s(nm,j)*vr(1,2,j) + xk3s(nm,j)*vr(1,3,j)
                         dusdx = dusdx - kx*sphase*sigma1s(nm,j)
                         dvsdx = dvsdx - kx*sphase*sigma2s(nm,j)
                         dwsdx = dwsdx - kx*sphase*sigma3s(nm,j)

                         ky = xk1s(nm,j)*vr(2,1,j) + xk2s(nm,j)*vr(2,2,j) + xk3s(nm,j)*vr(2,3,j)
                         dusdy = dusdy-ky*sphase*sigma1s(nm,j)
                         dvsdy = dvsdy-ky*sphase*sigma2s(nm,j)
                         dwsdy = dwsdy-ky*sphase*sigma3s(nm,j)

                         dusdt = dusdt - omn(nm)*sphase*sigma1s(nm,j)
                         dvsdt = dvsdt - omn(nm)*sphase*sigma2s(nm,j)
                         dwsdt = dwsdt - omn(nm)*sphase*sigma3s(nm,j)
                      enddo

                      ! retrieve anisotropic field from scaled velocity field: u= R u_scaled
                      ut=us*vr(1,1,j)+vs*vr(1,2,j)+ws*vr(1,3,j)
                      vt=us*vr(2,1,j)+vs*vr(2,2,j)+ws*vr(2,3,j)
                      wt=us*vr(3,1,j)+vs*vr(3,2,j)+ws*vr(3,3,j)

                      dutdx=dusdx*vr(1,1,j)+dvsdx*vr(1,2,j)+dwsdx*vr(1,3,j)
                      dvtdx=dusdx*vr(2,1,j)+dvsdx*vr(2,2,j)+dwsdx*vr(2,3,j)
                      dwtdx=dusdx*vr(3,1,j)+dvsdx*vr(3,2,j)+dwsdx*vr(3,3,j)

                      dutdy=dusdy*vr(1,1,j)+dvsdy*vr(1,2,j)+dwsdy*vr(1,3,j)
                      dvtdy=dusdy*vr(2,1,j)+dvsdy*vr(2,2,j)+dwsdy*vr(2,3,j)
                      dwtdy=dusdy*vr(3,1,j)+dvsdy*vr(3,2,j)+dwsdy*vr(3,3,j)

                      du_in=ut*BC_edge(1,1,2)%ir(i,l)+dutdx*BC_edge(1,1,2)%cosphi(i,l)+dutdy*BC_edge(1,1,2)%sinphi(i,l)
                      dv_in=vt*BC_edge(1,1,2)%ir(i,l)+dvtdx*BC_edge(1,1,2)%cosphi(i,l)+dvtdy*BC_edge(1,1,2)%sinphi(i,l)
                      dw_in=wt*BC_edge(1,1,2)%ir(i,l)+dwtdx*BC_edge(1,1,2)%cosphi(i,l)+dwtdy*BC_edge(1,1,2)%sinphi(i,l)

                      ut_in(i,l,k) = dusdt*vr(1,1,j) + dvsdt*vr(1,2,j) + dwsdt*vr(1,3,j) + vg(i,l,k)*du_in
                      vt_in(i,l,k) = dusdt*vr(2,1,j) + dvsdt*vr(2,2,j) + dwsdt*vr(2,3,j) + vg(i,l,k)*dv_in
                      wt_in(i,l,k) = dusdt*vr(3,1,j) + dvsdt*vr(3,2,j) + dwsdt*vr(3,3,j) + vg(i,l,k)*dw_in
                   enddo
                enddo
             enddo
          else if (is_curv) then
             do i=1,ngh
                do j=ny-ngh+1,ny
                   l=j-(ny-ngh)
                   do k=1,nz
                      ! scaled coordinates: x_scaled=R^T x
                      xs = (xc(i,j)-u_m*trk)*vr(1,1,j) + (yc(i,j)-v_m*trk)*vr(2,1,j) + (z(k)-w_m*trk)*vr(3,1,j)
                      ys = (xc(i,j)-u_m*trk)*vr(1,2,j) + (yc(i,j)-v_m*trk)*vr(2,2,j) + (z(k)-w_m*trk)*vr(3,2,j)
                      zs = (xc(i,j)-u_m*trk)*vr(1,3,j) + (yc(i,j)-v_m*trk)*vr(2,3,j) + (z(k)-w_m*trk)*vr(3,3,j)

                      ! summation over modes for scaled velocity components
                      us = 0.0_wp
                      vs = 0.0_wp
                      ws = 0.0_wp
                      dusdx = 0.0_wp
                      dvsdx = 0.0_wp
                      dwsdx = 0.0_wp
                      dusdy = 0.0_wp
                      dvsdy = 0.0_wp
                      dwsdy = 0.0_wp
                      dusdt = 0.0_wp
                      dvsdt = 0.0_wp
                      dwsdt = 0.0_wp
                      do nm=1,Nmode
                         phase = xk1s(nm,j)*xs+xk2s(nm,j)*ys+xk3s(nm,j)*zs + omn(nm)*trk + psi(nm)
                         cphase = cos(phase)
                         sphase = sin(phase)
                         us = us + cphase*sigma1s(nm,j)
                         vs = vs + cphase*sigma2s(nm,j)
                         ws = ws + cphase*sigma3s(nm,j)

                         kx = xk1s(nm,j)*vr(1,1,j) + xk2s(nm,j)*vr(1,2,j) + xk3s(nm,j)*vr(1,3,j)
                         dusdx = dusdx - kx*sphase*sigma1s(nm,j)
                         dvsdx = dvsdx - kx*sphase*sigma2s(nm,j)
                         dwsdx = dwsdx - kx*sphase*sigma3s(nm,j)

                         ky = xk1s(nm,j)*vr(2,1,j) + xk2s(nm,j)*vr(2,2,j) + xk3s(nm,j)*vr(2,3,j)
                         dusdy = dusdy-ky*sphase*sigma1s(nm,j)
                         dvsdy = dvsdy-ky*sphase*sigma2s(nm,j)
                         dwsdy = dwsdy-ky*sphase*sigma3s(nm,j)

                         dusdt = dusdt - omn(nm)*sphase*sigma1s(nm,j)
                         dvsdt = dvsdt - omn(nm)*sphase*sigma2s(nm,j)
                         dwsdt = dwsdt - omn(nm)*sphase*sigma3s(nm,j)
                      enddo

                      ! retrieve anisotropic field from scaled velocity field: u= R u_scaled
                      ut=us*vr(1,1,j)+vs*vr(1,2,j)+ws*vr(1,3,j)
                      vt=us*vr(2,1,j)+vs*vr(2,2,j)+ws*vr(2,3,j)
                      wt=us*vr(3,1,j)+vs*vr(3,2,j)+ws*vr(3,3,j)

                      dutdx=dusdx*vr(1,1,j)+dvsdx*vr(1,2,j)+dwsdx*vr(1,3,j)
                      dvtdx=dusdx*vr(2,1,j)+dvsdx*vr(2,2,j)+dwsdx*vr(2,3,j)
                      dwtdx=dusdx*vr(3,1,j)+dvsdx*vr(3,2,j)+dwsdx*vr(3,3,j)

                      dutdy=dusdy*vr(1,1,j)+dvsdy*vr(1,2,j)+dwsdy*vr(1,3,j)
                      dvtdy=dusdy*vr(2,1,j)+dvsdy*vr(2,2,j)+dwsdy*vr(2,3,j)
                      dwtdy=dusdy*vr(3,1,j)+dvsdy*vr(3,2,j)+dwsdy*vr(3,3,j)

                      du_in=ut*BC_edge(1,1,2)%ir(i,l)+dutdx*BC_edge(1,1,2)%cosphi(i,l)+dutdy*BC_edge(1,1,2)%sinphi(i,l)
                      dv_in=vt*BC_edge(1,1,2)%ir(i,l)+dvtdx*BC_edge(1,1,2)%cosphi(i,l)+dvtdy*BC_edge(1,1,2)%sinphi(i,l)
                      dw_in=wt*BC_edge(1,1,2)%ir(i,l)+dwtdx*BC_edge(1,1,2)%cosphi(i,l)+dwtdy*BC_edge(1,1,2)%sinphi(i,l)

                      ut_in(i,l,k) = dusdt*vr(1,1,j) + dvsdt*vr(1,2,j) + dwsdt*vr(1,3,j) + vg(i,l,k)*du_in
                      vt_in(i,l,k) = dusdt*vr(2,1,j) + dvsdt*vr(2,2,j) + dwsdt*vr(2,3,j) + vg(i,l,k)*dv_in
                      wt_in(i,l,k) = dusdt*vr(3,1,j) + dvsdt*vr(3,2,j) + dwsdt*vr(3,3,j) + vg(i,l,k)*dw_in
                   enddo
                enddo
             enddo
          else
             do i=1,ngh
                do j=ny-ngh+1,ny
                   l=j-(ny-ngh)
                   do k=1,nz
                      ! scaled coordinates: x_scaled=R^T x
                      xs = (x(i)-u_m*trk)*vr(1,1,j) + (y(j)-v_m*trk)*vr(2,1,j) + (z(k)-w_m*trk)*vr(3,1,j)
                      ys = (x(i)-u_m*trk)*vr(1,2,j) + (y(j)-v_m*trk)*vr(2,2,j) + (z(k)-w_m*trk)*vr(3,2,j)
                      zs = (x(i)-u_m*trk)*vr(1,3,j) + (y(j)-v_m*trk)*vr(2,3,j) + (z(k)-w_m*trk)*vr(3,3,j)

                      ! summation over modes for scaled velocity components
                      us = 0.0_wp
                      vs = 0.0_wp
                      ws = 0.0_wp
                      dusdx = 0.0_wp
                      dvsdx = 0.0_wp
                      dwsdx = 0.0_wp
                      dusdy = 0.0_wp
                      dvsdy = 0.0_wp
                      dwsdy = 0.0_wp
                      dusdt = 0.0_wp
                      dvsdt = 0.0_wp
                      dwsdt = 0.0_wp
                      do nm=1,Nmode
                         phase = xk1s(nm,j)*xs+xk2s(nm,j)*ys+xk3s(nm,j)*zs + omn(nm)*trk + psi(nm)
                         cphase = cos(phase)
                         sphase = sin(phase)
                         us = us + cphase*sigma1s(nm,j)
                         vs = vs + cphase*sigma2s(nm,j)
                         ws = ws + cphase*sigma3s(nm,j)

                         kx = xk1s(nm,j)*vr(1,1,j) + xk2s(nm,j)*vr(1,2,j) + xk3s(nm,j)*vr(1,3,j)
                         dusdx = dusdx - kx*sphase*sigma1s(nm,j)
                         dvsdx = dvsdx - kx*sphase*sigma2s(nm,j)
                         dwsdx = dwsdx - kx*sphase*sigma3s(nm,j)

                         ky = xk1s(nm,j)*vr(2,1,j) + xk2s(nm,j)*vr(2,2,j) + xk3s(nm,j)*vr(2,3,j)
                         dusdy = dusdy-ky*sphase*sigma1s(nm,j)
                         dvsdy = dvsdy-ky*sphase*sigma2s(nm,j)
                         dwsdy = dwsdy-ky*sphase*sigma3s(nm,j)

                         dusdt = dusdt - omn(nm)*sphase*sigma1s(nm,j)
                         dvsdt = dvsdt - omn(nm)*sphase*sigma2s(nm,j)
                         dwsdt = dwsdt - omn(nm)*sphase*sigma3s(nm,j)
                      enddo

                      ! retrieve anisotropic field from scaled velocity field: u= R u_scaled
                      ut=us*vr(1,1,j)+vs*vr(1,2,j)+ws*vr(1,3,j)
                      vt=us*vr(2,1,j)+vs*vr(2,2,j)+ws*vr(2,3,j)
                      wt=us*vr(3,1,j)+vs*vr(3,2,j)+ws*vr(3,3,j)

                      dutdx=dusdx*vr(1,1,j)+dvsdx*vr(1,2,j)+dwsdx*vr(1,3,j)
                      dvtdx=dusdx*vr(2,1,j)+dvsdx*vr(2,2,j)+dwsdx*vr(2,3,j)
                      dwtdx=dusdx*vr(3,1,j)+dvsdx*vr(3,2,j)+dwsdx*vr(3,3,j)

                      dutdy=dusdy*vr(1,1,j)+dvsdy*vr(1,2,j)+dwsdy*vr(1,3,j)
                      dvtdy=dusdy*vr(2,1,j)+dvsdy*vr(2,2,j)+dwsdy*vr(2,3,j)
                      dwtdy=dusdy*vr(3,1,j)+dvsdy*vr(3,2,j)+dwsdy*vr(3,3,j)

                      du_in=ut*BC_edge(1,1,2)%ir(i,l)+dutdx*BC_edge(1,1,2)%cosphi(i,l)+dutdy*BC_edge(1,1,2)%sinphi(i,l)
                      dv_in=vt*BC_edge(1,1,2)%ir(i,l)+dvtdx*BC_edge(1,1,2)%cosphi(i,l)+dvtdy*BC_edge(1,1,2)%sinphi(i,l)
                      dw_in=wt*BC_edge(1,1,2)%ir(i,l)+dwtdx*BC_edge(1,1,2)%cosphi(i,l)+dwtdy*BC_edge(1,1,2)%sinphi(i,l)

                      ut_in(i,l,k) = dusdt*vr(1,1,j) + dvsdt*vr(1,2,j) + dwsdt*vr(1,3,j) + vg(i,l,k)*du_in
                      vt_in(i,l,k) = dusdt*vr(2,1,j) + dvsdt*vr(2,2,j) + dwsdt*vr(2,3,j) + vg(i,l,k)*dv_in
                      wt_in(i,l,k) = dusdt*vr(3,1,j) + dvsdt*vr(3,2,j) + dwsdt*vr(3,3,j) + vg(i,l,k)*dw_in
                   enddo
                enddo
             enddo
          endif
       else
          ut_in=0.0_wp
          vt_in=0.0_wp
          wt_in=0.0_wp
          if (is_curv3) then
             do i=1,ngh
                do j=ny-ngh+1,ny
                   l=j-(ny-ngh)
                   do k=1,nz
                      ! summation over modes
                      ut=0.0_wp
                      vt=0.0_wp
                      wt=0.0_wp
                      dutdx=0.0_wp
                      dvtdx=0.0_wp
                      dwtdx=0.0_wp
                      dutdy=0.0_wp
                      dvtdy=0.0_wp
                      dwtdy=0.0_wp
                      do nm=1,Nmode
                         phase = xk1(nm)*(xc3(i,j,k)-u_m*trk) + xk2(nm)*(yc3(i,j,k)-v_m*trk) + xk3(nm)*(zc3(i,j,k)-w_m*trk) + omn(nm)*trk + psi(nm)
                         cphase=cos(phase)
                         sphase=sin(phase)
                         ut=ut+cphase*sigma1(nm)
                         vt=vt+cphase*sigma2(nm)
                         wt=wt+cphase*sigma3(nm)

                         dutdx=dutdx-xk1(nm)*sphase*sigma1(nm)
                         dvtdx=dvtdx-xk1(nm)*sphase*sigma2(nm)
                         dwtdx=dwtdx-xk1(nm)*sphase*sigma3(nm)

                         dutdy=dutdy-xk2(nm)*sphase*sigma1(nm)
                         dvtdy=dvtdy-xk2(nm)*sphase*sigma2(nm)
                         dwtdy=dwtdy-xk2(nm)*sphase*sigma3(nm)

                         ut_in(i,l,k) = ut_in(i,l,k) - (omn(nm)-xk1(nm)*u_m-xk2(nm)*v_m-xk3(nm)*w_m)*sphase*sigma1(nm)
                         vt_in(i,l,k) = vt_in(i,l,k) - (omn(nm)-xk1(nm)*u_m-xk2(nm)*v_m-xk3(nm)*w_m)*sphase*sigma2(nm)
                         wt_in(i,l,k) = wt_in(i,l,k) - (omn(nm)-xk1(nm)*u_m-xk2(nm)*v_m-xk3(nm)*w_m)*sphase*sigma3(nm)
                      enddo
                      du_in = ut*BC_edge(1,1,2)%ir(i,l)+dutdx*BC_edge(1,1,2)%cosphi(i,l)+dutdy*BC_edge(1,1,2)%sinphi(i,l)
                      dv_in = vt*BC_edge(1,1,2)%ir(i,l)+dvtdx*BC_edge(1,1,2)%cosphi(i,l)+dvtdy*BC_edge(1,1,2)%sinphi(i,l)
                      dw_in = wt*BC_edge(1,1,2)%ir(i,l)+dwtdx*BC_edge(1,1,2)%cosphi(i,l)+dwtdy*BC_edge(1,1,2)%sinphi(i,l)

                      ut_in(i,l,k) = ut_in(i,l,k) + vg(i,l,k)*du_in
                      vt_in(i,l,k) = vt_in(i,l,k) + vg(i,l,k)*dv_in
                      wt_in(i,l,k) = wt_in(i,l,k) + vg(i,l,k)*dw_in
                   enddo
                enddo
             enddo
          else if (is_curv) then
             do i=1,ngh
                do j=ny-ngh+1,ny
                   l=j-(ny-ngh)
                   do k=1,nz
                      ! summation over modes
                      ut=0.0_wp
                      vt=0.0_wp
                      wt=0.0_wp
                      dutdx=0.0_wp
                      dvtdx=0.0_wp
                      dwtdx=0.0_wp
                      dutdy=0.0_wp
                      dvtdy=0.0_wp
                      dwtdy=0.0_wp
                      do nm=1,Nmode
                         phase = xk1(nm)*(xc(i,j)-u_m*trk) + xk2(nm)*(yc(i,j)-v_m*trk) + xk3(nm)*(z(k)-w_m*trk) + omn(nm)*trk + psi(nm)
                         cphase=cos(phase)
                         sphase=sin(phase)
                         ut=ut+cphase*sigma1(nm)
                         vt=vt+cphase*sigma2(nm)
                         wt=wt+cphase*sigma3(nm)

                         dutdx=dutdx-xk1(nm)*sphase*sigma1(nm)
                         dvtdx=dvtdx-xk1(nm)*sphase*sigma2(nm)
                         dwtdx=dwtdx-xk1(nm)*sphase*sigma3(nm)

                         dutdy=dutdy-xk2(nm)*sphase*sigma1(nm)
                         dvtdy=dvtdy-xk2(nm)*sphase*sigma2(nm)
                         dwtdy=dwtdy-xk2(nm)*sphase*sigma3(nm)

                         ut_in(i,l,k) = ut_in(i,l,k) - (omn(nm)-xk1(nm)*u_m-xk2(nm)*v_m-xk3(nm)*w_m)*sphase*sigma1(nm)
                         vt_in(i,l,k) = vt_in(i,l,k) - (omn(nm)-xk1(nm)*u_m-xk2(nm)*v_m-xk3(nm)*w_m)*sphase*sigma2(nm)
                         wt_in(i,l,k) = wt_in(i,l,k) - (omn(nm)-xk1(nm)*u_m-xk2(nm)*v_m-xk3(nm)*w_m)*sphase*sigma3(nm)
                      enddo
                      du_in = ut*BC_edge(1,1,2)%ir(i,l)+dutdx*BC_edge(1,1,2)%cosphi(i,l)+dutdy*BC_edge(1,1,2)%sinphi(i,l)
                      dv_in = vt*BC_edge(1,1,2)%ir(i,l)+dvtdx*BC_edge(1,1,2)%cosphi(i,l)+dvtdy*BC_edge(1,1,2)%sinphi(i,l)
                      dw_in = wt*BC_edge(1,1,2)%ir(i,l)+dwtdx*BC_edge(1,1,2)%cosphi(i,l)+dwtdy*BC_edge(1,1,2)%sinphi(i,l)

                      ut_in(i,l,k) = ut_in(i,l,k) + vg(i,l,k)*du_in
                      vt_in(i,l,k) = vt_in(i,l,k) + vg(i,l,k)*dv_in
                      wt_in(i,l,k) = wt_in(i,l,k) + vg(i,l,k)*dw_in
                   enddo
                enddo
             enddo
          else
             do i=1,ngh
                do j=ny-ngh+1,ny
                   l=j-(ny-ngh)
                   do k=1,nz
                      ! summation over modes
                      ut=0.0_wp
                      vt=0.0_wp
                      wt=0.0_wp
                      dutdx=0.0_wp
                      dvtdx=0.0_wp
                      dwtdx=0.0_wp
                      dutdy=0.0_wp
                      dvtdy=0.0_wp
                      dwtdy=0.0_wp
                      do nm=1,Nmode
                         phase=xk1(nm)*(x(i)-u_m*trk) + xk2(nm)*(y(j)-v_m*trk) + xk3(nm)*(z(k)-w_m*trk) + omn(nm)*trk + psi(nm)
                         cphase=cos(phase)
                         sphase=sin(phase)
                         ut=ut+cphase*sigma1(nm)
                         vt=vt+cphase*sigma2(nm)
                         wt=wt+cphase*sigma3(nm)

                         dutdx=dutdx-xk1(nm)*sphase*sigma1(nm)
                         dvtdx=dvtdx-xk1(nm)*sphase*sigma2(nm)
                         dwtdx=dwtdx-xk1(nm)*sphase*sigma3(nm)

                         dutdy=dutdy-xk2(nm)*sphase*sigma1(nm)
                         dvtdy=dvtdy-xk2(nm)*sphase*sigma2(nm)
                         dwtdy=dwtdy-xk2(nm)*sphase*sigma3(nm)

                         ut_in(i,l,k) = ut_in(i,l,k) - (omn(nm)-xk1(nm)*u_m-xk2(nm)*v_m-xk3(nm)*w_m)*sphase*sigma1(nm)
                         vt_in(i,l,k) = vt_in(i,l,k) - (omn(nm)-xk1(nm)*u_m-xk2(nm)*v_m-xk3(nm)*w_m)*sphase*sigma2(nm)
                         wt_in(i,l,k) = wt_in(i,l,k) - (omn(nm)-xk1(nm)*u_m-xk2(nm)*v_m-xk3(nm)*w_m)*sphase*sigma3(nm)
                      enddo
                      du_in = ut*BC_edge(1,1,2)%ir(i,l)+dutdx*BC_edge(1,1,2)%cosphi(i,l)+dutdy*BC_edge(1,1,2)%sinphi(i,l)
                      dv_in = vt*BC_edge(1,1,2)%ir(i,l)+dvtdx*BC_edge(1,1,2)%cosphi(i,l)+dvtdy*BC_edge(1,1,2)%sinphi(i,l)
                      dw_in = wt*BC_edge(1,1,2)%ir(i,l)+dwtdx*BC_edge(1,1,2)%cosphi(i,l)+dwtdy*BC_edge(1,1,2)%sinphi(i,l)

                      ut_in(i,l,k) = ut_in(i,l,k) + vg(i,l,k)*du_in
                      vt_in(i,l,k) = vt_in(i,l,k) + vg(i,l,k)*dv_in
                      wt_in(i,l,k) = wt_in(i,l,k) + vg(i,l,k)*dw_in
                   enddo
                enddo
             enddo
          endif
       endif

       ! Apply weighting function
       ! ========================
       if (anisotropy=='W') then
          if (is_bc_wall(2,2)) then
             do i=1,ngh
                do j=ny-ngh+1,ny
                   l=j-(ny-ngh)
                   fy=(y(j)-y(1))/(y(ngh+1)-y(1))
                   do k=1,nz
                      ! add linear damping close to the wall
                      ut_in(i,l,k) = ut_in(i,l,k)*u_rms(j)*fy
                      vt_in(i,l,k) = vt_in(i,l,k)*v_rms(j)*fy
                      wt_in(i,l,k) = wt_in(i,l,k)*w_rms(j)*fy
                   enddo
                enddo
             enddo
          else
             do i=1,ngh
                do j=ny-ngh+1,ny
                   l=j-(ny-ngh)
                   do k=1,nz
                      ! multiply by rms profiles
                      ut_in(i,l,k) = ut_in(i,l,k)*u_rms(j)
                      vt_in(i,l,k) = vt_in(i,l,k)*v_rms(j)
                      wt_in(i,l,k) = wt_in(i,l,k)*w_rms(j)
                   enddo
                enddo
             enddo
          endif
       endif

       ! Apply damping coefficients
       ! ==========================
       if (is_curv3) then
          do i=1,ngh
             do j=ny-ngh+1,ny
                l=j-(ny-ngh)
                do k=1,nz
                   ut_in(i,l,k) = ut_in(i,l,k)*damping_coeff3(j,k)
                   vt_in(i,l,k) = vt_in(i,l,k)*damping_coeff3(j,k)
                   wt_in(i,l,k) = wt_in(i,l,k)*damping_coeff3(j,k)
                enddo
             enddo
          enddo
       else
          do i=1,ngh
             do j=ny-ngh+1,ny
                l=j-(ny-ngh)
                do k=1,nz
                   ut_in(i,l,k) = ut_in(i,l,k)*damping_coeff(j)
                   vt_in(i,l,k) = vt_in(i,l,k)*damping_coeff(j)
                   wt_in(i,l,k) = wt_in(i,l,k)*damping_coeff(j)
                enddo
             enddo
          enddo
       endif
    endif

  end subroutine disturb_inlet_RFM_TamDong_imin_jmax

  !===============================================================================
  module subroutine disturb_inlet_RFM_TamDong1pt_imin(vg,ut_in,vt_in,wt_in)
  !===============================================================================
    !> Impose RFM disturbances at inlet (for Tam & Dong's BC)
  !===============================================================================
    use mod_time
    implicit none
    ! ----------------------------------------------------------------------------
    ! group velocity
    real(wp), dimension(ny,nz), intent(in) :: vg
    ! RFM disturbances
    real(wp), dimension(ny,nz), intent(out) :: ut_in,vt_in,wt_in
    ! ----------------------------------------------------------------------------
    ! local var
    integer :: i,j,k,nm
    real(wp) :: du_in,dv_in,dw_in
    real(wp) :: trk,kx,ky
    real(wp) :: phase,cphase,sphase
    real(wp) :: ut,vt,wt
    real(wp) :: dutdx,dvtdx,dwtdx,dutdy,dvtdy,dwtdy
    real(wp) :: xs,ys,zs ! scaled cooordinates
    ! scaled velocity components
    real(wp) :: us,vs,ws
    real(wp) :: dusdx,dvsdx,dwsdx
    real(wp) :: dusdy,dvsdy,dwsdy
    real(wp) :: dusdt,dvsdt,dwsdt
    ! ----------------------------------------------------------------------------

    ! Index of left boundary
    ! ======================
    i=1

    trk = time + ck(irk)*deltat

    ! Construct turbulent stochastic field
    ! ------------------------------------
    if (anisotropy=='S') then
       if (is_curv3) then
          do j=ndy_td1,nfy_td1
             do k=1,nz
                ! scaled coordinates: x_scaled=R^T x
                xs = (xc3(i,j,k)-u_m*trk)*vr(1,1,j) + (yc3(i,j,k)-v_m*trk)*vr(2,1,j) + (zc3(i,j,k)-w_m*trk)*vr(3,1,j)
                ys = (xc3(i,j,k)-u_m*trk)*vr(1,2,j) + (yc3(i,j,k)-v_m*trk)*vr(2,2,j) + (zc3(i,j,k)-w_m*trk)*vr(3,2,j)
                zs = (xc3(i,j,k)-u_m*trk)*vr(1,3,j) + (yc3(i,j,k)-v_m*trk)*vr(2,3,j) + (zc3(i,j,k)-w_m*trk)*vr(3,3,j)

                ! summation over modes for scaled velocity components
                us = 0.0_wp
                vs = 0.0_wp
                ws = 0.0_wp
                dusdx = 0.0_wp
                dvsdx = 0.0_wp
                dwsdx = 0.0_wp
                dusdy = 0.0_wp
                dvsdy = 0.0_wp
                dwsdy = 0.0_wp
                dusdt = 0.0_wp
                dvsdt = 0.0_wp
                dwsdt = 0.0_wp
                do nm=1,Nmode
                   phase = xk1s(nm,j)*xs+xk2s(nm,j)*ys+xk3s(nm,j)*zs + omn(nm)*trk + psi(nm)
                   cphase = cos(phase)
                   sphase = sin(phase)
                   us = us + cphase*sigma1s(nm,j)
                   vs = vs + cphase*sigma2s(nm,j)
                   ws = ws + cphase*sigma3s(nm,j)

                   kx = xk1s(nm,j)*vr(1,1,j) + xk2s(nm,j)*vr(1,2,j) + xk3s(nm,j)*vr(1,3,j)
                   dusdx = dusdx - kx*sphase*sigma1s(nm,j)
                   dvsdx = dvsdx - kx*sphase*sigma2s(nm,j)
                   dwsdx = dwsdx - kx*sphase*sigma3s(nm,j)

                   ky = xk1s(nm,j)*vr(2,1,j) + xk2s(nm,j)*vr(2,2,j) + xk3s(nm,j)*vr(2,3,j)
                   dusdy = dusdy-ky*sphase*sigma1s(nm,j)
                   dvsdy = dvsdy-ky*sphase*sigma2s(nm,j)
                   dwsdy = dwsdy-ky*sphase*sigma3s(nm,j)

                   dusdt = dusdt - omn(nm)*sphase*sigma1s(nm,j)
                   dvsdt = dvsdt - omn(nm)*sphase*sigma2s(nm,j)
                   dwsdt = dwsdt - omn(nm)*sphase*sigma3s(nm,j)
                enddo

                ! retrieve anisotropic field from scaled velocity field: u= R u_scaled
                ut=us*vr(1,1,j)+vs*vr(1,2,j)+ws*vr(1,3,j)
                vt=us*vr(2,1,j)+vs*vr(2,2,j)+ws*vr(2,3,j)
                wt=us*vr(3,1,j)+vs*vr(3,2,j)+ws*vr(3,3,j)

                dutdx=dusdx*vr(1,1,j)+dvsdx*vr(1,2,j)+dwsdx*vr(1,3,j)
                dvtdx=dusdx*vr(2,1,j)+dvsdx*vr(2,2,j)+dwsdx*vr(2,3,j)
                dwtdx=dusdx*vr(3,1,j)+dvsdx*vr(3,2,j)+dwsdx*vr(3,3,j)

                dutdy=dusdy*vr(1,1,j)+dvsdy*vr(1,2,j)+dwsdy*vr(1,3,j)
                dvtdy=dusdy*vr(2,1,j)+dvsdy*vr(2,2,j)+dwsdy*vr(2,3,j)
                dwtdy=dusdy*vr(3,1,j)+dvsdy*vr(3,2,j)+dwsdy*vr(3,3,j)

                du_in=ut*BC_face(1,1)%ir(i,j)+dutdx*BC_face(1,1)%cosphi(i,j)+dutdy*BC_face(1,1)%sinphi(i,j)
                dv_in=vt*BC_face(1,1)%ir(i,j)+dvtdx*BC_face(1,1)%cosphi(i,j)+dvtdy*BC_face(1,1)%sinphi(i,j)
                dw_in=wt*BC_face(1,1)%ir(i,j)+dwtdx*BC_face(1,1)%cosphi(i,j)+dwtdy*BC_face(1,1)%sinphi(i,j)

                ut_in(j,k) = dusdt*vr(1,1,j) + dvsdt*vr(1,2,j) + dwsdt*vr(1,3,j) + vg(j,k)*du_in
                vt_in(j,k) = dusdt*vr(2,1,j) + dvsdt*vr(2,2,j) + dwsdt*vr(2,3,j) + vg(j,k)*dv_in
                wt_in(j,k) = dusdt*vr(3,1,j) + dvsdt*vr(3,2,j) + dwsdt*vr(3,3,j) + vg(j,k)*dw_in
             enddo
          enddo
       else if (is_curv) then
          do j=ndy_td1,nfy_td1
             do k=1,nz
                ! scaled coordinates: x_scaled=R^T x
                xs = (xc(i,j)-u_m*trk)*vr(1,1,j) + (yc(i,j)-v_m*trk)*vr(2,1,j) + (z(k)-w_m*trk)*vr(3,1,j)
                ys = (xc(i,j)-u_m*trk)*vr(1,2,j) + (yc(i,j)-v_m*trk)*vr(2,2,j) + (z(k)-w_m*trk)*vr(3,2,j)
                zs = (xc(i,j)-u_m*trk)*vr(1,3,j) + (yc(i,j)-v_m*trk)*vr(2,3,j) + (z(k)-w_m*trk)*vr(3,3,j)

                ! summation over modes for scaled velocity components
                us = 0.0_wp
                vs = 0.0_wp
                ws = 0.0_wp
                dusdx = 0.0_wp
                dvsdx = 0.0_wp
                dwsdx = 0.0_wp
                dusdy = 0.0_wp
                dvsdy = 0.0_wp
                dwsdy = 0.0_wp
                dusdt = 0.0_wp
                dvsdt = 0.0_wp
                dwsdt = 0.0_wp
                do nm=1,Nmode
                   phase = xk1s(nm,j)*xs+xk2s(nm,j)*ys+xk3s(nm,j)*zs + omn(nm)*trk + psi(nm)
                   cphase = cos(phase)
                   sphase = sin(phase)
                   us = us + cphase*sigma1s(nm,j)
                   vs = vs + cphase*sigma2s(nm,j)
                   ws = ws + cphase*sigma3s(nm,j)

                   kx = xk1s(nm,j)*vr(1,1,j) + xk2s(nm,j)*vr(1,2,j) + xk3s(nm,j)*vr(1,3,j)
                   dusdx = dusdx - kx*sphase*sigma1s(nm,j)
                   dvsdx = dvsdx - kx*sphase*sigma2s(nm,j)
                   dwsdx = dwsdx - kx*sphase*sigma3s(nm,j)

                   ky = xk1s(nm,j)*vr(2,1,j) + xk2s(nm,j)*vr(2,2,j) + xk3s(nm,j)*vr(2,3,j)
                   dusdy = dusdy-ky*sphase*sigma1s(nm,j)
                   dvsdy = dvsdy-ky*sphase*sigma2s(nm,j)
                   dwsdy = dwsdy-ky*sphase*sigma3s(nm,j)

                   dusdt = dusdt - omn(nm)*sphase*sigma1s(nm,j)
                   dvsdt = dvsdt - omn(nm)*sphase*sigma2s(nm,j)
                   dwsdt = dwsdt - omn(nm)*sphase*sigma3s(nm,j)
                enddo

                ! retrieve anisotropic field from scaled velocity field: u= R u_scaled
                ut=us*vr(1,1,j)+vs*vr(1,2,j)+ws*vr(1,3,j)
                vt=us*vr(2,1,j)+vs*vr(2,2,j)+ws*vr(2,3,j)
                wt=us*vr(3,1,j)+vs*vr(3,2,j)+ws*vr(3,3,j)

                dutdx=dusdx*vr(1,1,j)+dvsdx*vr(1,2,j)+dwsdx*vr(1,3,j)
                dvtdx=dusdx*vr(2,1,j)+dvsdx*vr(2,2,j)+dwsdx*vr(2,3,j)
                dwtdx=dusdx*vr(3,1,j)+dvsdx*vr(3,2,j)+dwsdx*vr(3,3,j)

                dutdy=dusdy*vr(1,1,j)+dvsdy*vr(1,2,j)+dwsdy*vr(1,3,j)
                dvtdy=dusdy*vr(2,1,j)+dvsdy*vr(2,2,j)+dwsdy*vr(2,3,j)
                dwtdy=dusdy*vr(3,1,j)+dvsdy*vr(3,2,j)+dwsdy*vr(3,3,j)

                du_in=ut*BC_face(1,1)%ir(i,j)+dutdx*BC_face(1,1)%cosphi(i,j)+dutdy*BC_face(1,1)%sinphi(i,j)
                dv_in=vt*BC_face(1,1)%ir(i,j)+dvtdx*BC_face(1,1)%cosphi(i,j)+dvtdy*BC_face(1,1)%sinphi(i,j)
                dw_in=wt*BC_face(1,1)%ir(i,j)+dwtdx*BC_face(1,1)%cosphi(i,j)+dwtdy*BC_face(1,1)%sinphi(i,j)

                ut_in(j,k) = dusdt*vr(1,1,j) + dvsdt*vr(1,2,j) + dwsdt*vr(1,3,j) + vg(j,k)*du_in
                vt_in(j,k) = dusdt*vr(2,1,j) + dvsdt*vr(2,2,j) + dwsdt*vr(2,3,j) + vg(j,k)*dv_in
                wt_in(j,k) = dusdt*vr(3,1,j) + dvsdt*vr(3,2,j) + dwsdt*vr(3,3,j) + vg(j,k)*dw_in
             enddo
          enddo
       else
          do j=ndy_td1,nfy_td1
             do k=1,nz
                ! scaled coordinates: x_scaled=R^T x
                xs = (x(i)-u_m*trk)*vr(1,1,j) + (y(j)-v_m*trk)*vr(2,1,j) + (z(k)-w_m*trk)*vr(3,1,j)
                ys = (x(i)-u_m*trk)*vr(1,2,j) + (y(j)-v_m*trk)*vr(2,2,j) + (z(k)-w_m*trk)*vr(3,2,j)
                zs = (x(i)-u_m*trk)*vr(1,3,j) + (y(j)-v_m*trk)*vr(2,3,j) + (z(k)-w_m*trk)*vr(3,3,j)

                ! summation over modes for scaled velocity components
                us = 0.0_wp
                vs = 0.0_wp
                ws = 0.0_wp
                dusdx = 0.0_wp
                dvsdx = 0.0_wp
                dwsdx = 0.0_wp
                dusdy = 0.0_wp
                dvsdy = 0.0_wp
                dwsdy = 0.0_wp
                dusdt = 0.0_wp
                dvsdt = 0.0_wp
                dwsdt = 0.0_wp
                do nm=1,Nmode
                   phase = xk1s(nm,j)*xs+xk2s(nm,j)*ys+xk3s(nm,j)*zs + omn(nm)*trk + psi(nm)
                   cphase = cos(phase)
                   sphase = sin(phase)
                   us = us + cphase*sigma1s(nm,j)
                   vs = vs + cphase*sigma2s(nm,j)
                   ws = ws + cphase*sigma3s(nm,j)

                   kx = xk1s(nm,j)*vr(1,1,j) + xk2s(nm,j)*vr(1,2,j) + xk3s(nm,j)*vr(1,3,j)
                   dusdx = dusdx - kx*sphase*sigma1s(nm,j)
                   dvsdx = dvsdx - kx*sphase*sigma2s(nm,j)
                   dwsdx = dwsdx - kx*sphase*sigma3s(nm,j)

                   ky = xk1s(nm,j)*vr(2,1,j) + xk2s(nm,j)*vr(2,2,j) + xk3s(nm,j)*vr(2,3,j)
                   dusdy = dusdy - ky*sphase*sigma1s(nm,j)
                   dvsdy = dvsdy - ky*sphase*sigma2s(nm,j)
                   dwsdy = dwsdy - ky*sphase*sigma3s(nm,j)

                   dusdt = dusdt - omn(nm)*sphase*sigma1s(nm,j)
                   dvsdt = dvsdt - omn(nm)*sphase*sigma2s(nm,j)
                   dwsdt = dwsdt - omn(nm)*sphase*sigma3s(nm,j)
                enddo

                ! retrieve anisotropic field from scaled velocity field: u= R u_scaled
                ut=us*vr(1,1,j)+vs*vr(1,2,j)+ws*vr(1,3,j)
                vt=us*vr(2,1,j)+vs*vr(2,2,j)+ws*vr(2,3,j)
                wt=us*vr(3,1,j)+vs*vr(3,2,j)+ws*vr(3,3,j)

                dutdx=dusdx*vr(1,1,j)+dvsdx*vr(1,2,j)+dwsdx*vr(1,3,j)
                dvtdx=dusdx*vr(2,1,j)+dvsdx*vr(2,2,j)+dwsdx*vr(2,3,j)
                dwtdx=dusdx*vr(3,1,j)+dvsdx*vr(3,2,j)+dwsdx*vr(3,3,j)

                dutdy=dusdy*vr(1,1,j)+dvsdy*vr(1,2,j)+dwsdy*vr(1,3,j)
                dvtdy=dusdy*vr(2,1,j)+dvsdy*vr(2,2,j)+dwsdy*vr(2,3,j)
                dwtdy=dusdy*vr(3,1,j)+dvsdy*vr(3,2,j)+dwsdy*vr(3,3,j)

                du_in=ut*BC_face(1,1)%ir(i,j)+dutdx*BC_face(1,1)%cosphi(i,j)+dutdy*BC_face(1,1)%sinphi(i,j)
                dv_in=vt*BC_face(1,1)%ir(i,j)+dvtdx*BC_face(1,1)%cosphi(i,j)+dvtdy*BC_face(1,1)%sinphi(i,j)
                dw_in=wt*BC_face(1,1)%ir(i,j)+dwtdx*BC_face(1,1)%cosphi(i,j)+dwtdy*BC_face(1,1)%sinphi(i,j)

                ut_in(j,k) = dusdt*vr(1,1,j) + dvsdt*vr(1,2,j) + dwsdt*vr(1,3,j) + vg(j,k)*du_in
                vt_in(j,k) = dusdt*vr(2,1,j) + dvsdt*vr(2,2,j) + dwsdt*vr(2,3,j) + vg(j,k)*dv_in
                wt_in(j,k) = dusdt*vr(3,1,j) + dvsdt*vr(3,2,j) + dwsdt*vr(3,3,j) + vg(j,k)*dw_in
             enddo
          enddo
       endif
    else
       ut_in=0.0_wp
       vt_in=0.0_wp
       wt_in=0.0_wp
       if (is_curv3) then
          do j=ndy_td1,nfy_td1
             do k=1,nz
                ! summation over modes
                ut=0.0_wp
                vt=0.0_wp
                wt=0.0_wp
                dutdx=0.0_wp
                dvtdx=0.0_wp
                dwtdx=0.0_wp
                dutdy=0.0_wp
                dvtdy=0.0_wp
                dwtdy=0.0_wp
                do nm=1,Nmode
                   phase = xk1(nm)*(xc3(i,j,k)-u_m*trk) + xk2(nm)*(yc3(i,j,k)-v_m*trk) + xk3(nm)*(zc3(i,j,k)-w_m*trk) + omn(nm)*trk + psi(nm)
                   cphase=cos(phase)
                   sphase=sin(phase)

                   ut = ut + cphase*sigma1(nm)
                   vt = vt + cphase*sigma2(nm)
                   wt = wt + cphase*sigma3(nm)

                   dutdy = dutdy - xk2(nm)*sphase*sigma1(nm)
                   dvtdy = dvtdy - xk2(nm)*sphase*sigma2(nm)
                   dwtdy = dwtdy - xk2(nm)*sphase*sigma3(nm)

                   ut_in(j,k) = ut_in(j,k) - (omn(nm)-xk1(nm)*u_m-xk2(nm)*v_m-xk3(nm)*w_m)*sphase*sigma1(nm)
                   vt_in(j,k) = vt_in(j,k) - (omn(nm)-xk1(nm)*u_m-xk2(nm)*v_m-xk3(nm)*w_m)*sphase*sigma2(nm)
                   wt_in(j,k) = wt_in(j,k) - (omn(nm)-xk1(nm)*u_m-xk2(nm)*v_m-xk3(nm)*w_m)*sphase*sigma3(nm)
                enddo
                du_in = ut*BC_face(1,1)%ir(i,j)+dutdx*BC_face(1,1)%cosphi(i,j)+dutdy*BC_face(1,1)%sinphi(i,j)
                dv_in = vt*BC_face(1,1)%ir(i,j)+dvtdx*BC_face(1,1)%cosphi(i,j)+dvtdy*BC_face(1,1)%sinphi(i,j)
                dw_in = wt*BC_face(1,1)%ir(i,j)+dwtdx*BC_face(1,1)%cosphi(i,j)+dwtdy*BC_face(1,1)%sinphi(i,j)

                ut_in(j,k) = ut_in(j,k) + vg(j,k)*du_in
                vt_in(j,k) = vt_in(j,k) + vg(j,k)*dv_in
                wt_in(j,k) = wt_in(j,k) + vg(j,k)*dw_in
             enddo
          enddo
       else if (is_curv) then
          do j=ndy_td1,nfy_td1
             do k=1,nz
                ! summation over modes
                ut=0.0_wp
                vt=0.0_wp
                wt=0.0_wp
                dutdx=0.0_wp
                dvtdx=0.0_wp
                dwtdx=0.0_wp
                dutdy=0.0_wp
                dvtdy=0.0_wp
                dwtdy=0.0_wp
                do nm=1,Nmode
                   phase = xk1(nm)*(xc(i,j)-u_m*trk) + xk2(nm)*(yc(i,j)-v_m*trk) + xk3(nm)*(z(k)-w_m*trk) + omn(nm)*trk + psi(nm)
                   cphase=cos(phase)
                   sphase=sin(phase)

                   ut = ut + cphase*sigma1(nm)
                   vt = vt + cphase*sigma2(nm)
                   wt = wt + cphase*sigma3(nm)

                   dutdy = dutdy - xk2(nm)*sphase*sigma1(nm)
                   dvtdy = dvtdy - xk2(nm)*sphase*sigma2(nm)
                   dwtdy = dwtdy - xk2(nm)*sphase*sigma3(nm)

                   ut_in(j,k) = ut_in(j,k) - (omn(nm)-xk1(nm)*u_m-xk2(nm)*v_m-xk3(nm)*w_m)*sphase*sigma1(nm)
                   vt_in(j,k) = vt_in(j,k) - (omn(nm)-xk1(nm)*u_m-xk2(nm)*v_m-xk3(nm)*w_m)*sphase*sigma2(nm)
                   wt_in(j,k) = wt_in(j,k) - (omn(nm)-xk1(nm)*u_m-xk2(nm)*v_m-xk3(nm)*w_m)*sphase*sigma3(nm)
                enddo
                du_in = ut*BC_face(1,1)%ir(i,j)+dutdx*BC_face(1,1)%cosphi(i,j)+dutdy*BC_face(1,1)%sinphi(i,j)
                dv_in = vt*BC_face(1,1)%ir(i,j)+dvtdx*BC_face(1,1)%cosphi(i,j)+dvtdy*BC_face(1,1)%sinphi(i,j)
                dw_in = wt*BC_face(1,1)%ir(i,j)+dwtdx*BC_face(1,1)%cosphi(i,j)+dwtdy*BC_face(1,1)%sinphi(i,j)

                ut_in(j,k) = ut_in(j,k) + vg(j,k)*du_in
                vt_in(j,k) = vt_in(j,k) + vg(j,k)*dv_in
                wt_in(j,k) = wt_in(j,k) + vg(j,k)*dw_in
             enddo
          enddo
       else
          do j=ndy_td1,nfy_td1
             do k=1,nz
                ! summation over modes
                ut=0.0_wp
                vt=0.0_wp
                wt=0.0_wp
                dutdx=0.0_wp
                dvtdx=0.0_wp
                dwtdx=0.0_wp
                dutdy=0.0_wp
                dvtdy=0.0_wp
                dwtdy=0.0_wp
                do nm=1,Nmode
                   phase = xk1(nm)*(x(i)-u_m*trk) + xk2(nm)*(y(j)-v_m*trk) + xk3(nm)*(z(k)-w_m*trk) + omn(nm)*trk + psi(nm)
                   cphase=cos(phase)
                   sphase=sin(phase)

                   ut = ut + cphase*sigma1(nm)
                   vt = vt + cphase*sigma2(nm)
                   wt = wt + cphase*sigma3(nm)

                   dutdy = dutdy - xk2(nm)*sphase*sigma1(nm)
                   dvtdy = dvtdy - xk2(nm)*sphase*sigma2(nm)
                   dwtdy = dwtdy - xk2(nm)*sphase*sigma3(nm)

                   ut_in(j,k) = ut_in(j,k) - (omn(nm)-xk1(nm)*u_m-xk2(nm)*v_m-xk3(nm)*w_m)*sphase*sigma1(nm)
                   vt_in(j,k) = vt_in(j,k) - (omn(nm)-xk1(nm)*u_m-xk2(nm)*v_m-xk3(nm)*w_m)*sphase*sigma2(nm)
                   wt_in(j,k) = wt_in(j,k) - (omn(nm)-xk1(nm)*u_m-xk2(nm)*v_m-xk3(nm)*w_m)*sphase*sigma3(nm)
                enddo
                du_in = ut*BC_face(1,1)%ir(i,j)+dutdx*BC_face(1,1)%cosphi(i,j)+dutdy*BC_face(1,1)%sinphi(i,j)
                dv_in = vt*BC_face(1,1)%ir(i,j)+dvtdx*BC_face(1,1)%cosphi(i,j)+dvtdy*BC_face(1,1)%sinphi(i,j)
                dw_in = wt*BC_face(1,1)%ir(i,j)+dwtdx*BC_face(1,1)%cosphi(i,j)+dwtdy*BC_face(1,1)%sinphi(i,j)

                ut_in(j,k) = ut_in(j,k) + vg(j,k)*du_in
                vt_in(j,k) = vt_in(j,k) + vg(j,k)*dv_in
                wt_in(j,k) = wt_in(j,k) + vg(j,k)*dw_in
             enddo
          enddo
       endif
    endif

    ! Apply weighting function
    ! ========================
    if (anisotropy=='W') then
       do i=1,ngh
          do j=ndy_td1,nfy_td1
             do k=1,nz
                ! multiply by rms profiles
                ut_in(j,k) = ut_in(j,k)*u_rms(j)
                vt_in(j,k) = vt_in(j,k)*v_rms(j)
                wt_in(j,k) = wt_in(j,k)*w_rms(j)
             enddo
          enddo
       enddo
    endif

    ! Apply damping coefficients
    ! ==========================
    if (is_curv3) then
       do i=1,ngh
          do j=ndy_td1,nfy_td1
             do k=1,nz
                ut_in(j,k) = ut_in(j,k)*damping_coeff3(j,k)
                vt_in(j,k) = vt_in(j,k)*damping_coeff3(j,k)
                wt_in(j,k) = wt_in(j,k)*damping_coeff3(j,k)
             enddo
          enddo
       enddo
    else
       do i=1,ngh
          do j=ndy_td1,nfy_td1
             do k=1,nz
                ut_in(j,k) = ut_in(j,k)*damping_coeff(j)
                vt_in(j,k) = vt_in(j,k)*damping_coeff(j)
                wt_in(j,k) = wt_in(j,k)*damping_coeff(j)
             enddo
          enddo
       enddo
    endif

  end subroutine disturb_inlet_RFM_TamDong1pt_imin

end submodule smod_RFM_TamDong
