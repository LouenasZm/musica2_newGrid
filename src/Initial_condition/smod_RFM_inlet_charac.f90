!===============================================================================
submodule (mod_RFM) smod_RFM_charac
!===============================================================================
  !> author: AB & XG
  !> date: February 2024
  !> Generation of Random Fourier Modes (RFM)
  !> Injection in characteristic BCs [only for imin]
!=============================================================================== 

contains

  !===============================================================================
  module subroutine disturb_inlet_RFM_charac(ut_in,vt_in,wt_in)
  !===============================================================================
    !> Impose RFM disturbances at inlet (for charac BC)
  !===============================================================================
    use mod_time
    implicit none
    ! ----------------------------------------------------------------------------
    ! RFM disturbances
    real(wp), dimension(ny,nz), intent(out) :: ut_in,vt_in,wt_in
    ! ----------------------------------------------------------------------------
    ! local var
    integer :: i,j,k,nm
    real(wp) :: trk,fy
    real(wp) :: phase,dsphase
    real(wp) :: xs,ys,zs ! scaled cooordinates
    ! scaled velocity components
    real(wp) :: dusdt,dvsdt,dwsdt
    ! ----------------------------------------------------------------------------
    ! TBL injection
    integer :: ndeb_j1,nend_j1,jbc
    integer, dimension(2) :: ndeb_j2,nend_j2

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

    i=1
    
    ! physical time at current RK stage
    ! ---------------------------------
    trk=time+ck(irk)*deltat

    ! -------------------
    ! I/ Injection of FST
    ! -------------------
    if (is_RFM_FST) then
       ! Construct turbulent stochastic field
       ! ------------------------------------
       if (anisotropy=='S') then
          if (is_curv3) then
             do k=1,nz
                do j=ndy_c,nfy_c
                   ! scaled coordinates: x_scaled=R^T x
                   xs=(xc3(i,j,k)-u_m*trk)*vr(1,1,j)+(yc3(i,j,k)-v_m*trk)*vr(2,1,j)+(zc3(i,j,k)-w_m*trk)*vr(3,1,j)
                   ys=(xc3(i,j,k)-u_m*trk)*vr(1,2,j)+(yc3(i,j,k)-v_m*trk)*vr(2,2,j)+(zc3(i,j,k)-w_m*trk)*vr(3,2,j)
                   zs=(xc3(i,j,k)-u_m*trk)*vr(1,3,j)+(yc3(i,j,k)-v_m*trk)*vr(2,3,j)+(zc3(i,j,k)-w_m*trk)*vr(3,3,j)

                   ! summation over modes for scaled velocity components
                   dusdt= 0.0_wp
                   dvsdt= 0.0_wp
                   dwsdt= 0.0_wp
                   do nm=1,Nmode
                      phase = xk1s(nm,j)*xs+xk2s(nm,j)*ys+xk3s(nm,j)*zs+omn(nm)*trk + psi(nm)
                      dsphase=(omn(nm)-xk1(nm)*u_m-xk2(nm)*v_m-xk3(nm)*w_m)*sin(phase)

                      dusdt= dusdt - dsphase*sigma1s(nm,j)
                      dvsdt= dvsdt - dsphase*sigma2s(nm,j)
                      dwsdt= dwsdt - dsphase*sigma3s(nm,j)
                   enddo

                   ! retrieve anisotropic field from scaled velocity field: u= R u_scaled
                   ut_in(j,k)= dusdt*vr(1,1,j) + dvsdt*vr(1,2,j) + dwsdt*vr(1,3,j)
                   vt_in(j,k)= dusdt*vr(2,1,j) + dvsdt*vr(2,2,j) + dwsdt*vr(2,3,j)
                   wt_in(j,k)= dusdt*vr(3,1,j) + dvsdt*vr(3,2,j) + dwsdt*vr(3,3,j)
                enddo
             enddo
          elseif (is_curv) then
             do k=1,nz
                do j=ndy_c,nfy_c
                   ! scaled coordinates: x_scaled=R^T x
                   xs= (xc(i,j)-u_m*trk)*vr(1,1,j)+(yc(i,j)-v_m*trk)*vr(2,1,j)+(z(k)-w_m*trk)*vr(3,1,j)
                   ys= (xc(i,j)-u_m*trk)*vr(1,2,j)+(yc(i,j)-v_m*trk)*vr(2,2,j)+(z(k)-w_m*trk)*vr(3,2,j)
                   zs= (xc(i,j)-u_m*trk)*vr(1,3,j)+(yc(i,j)-v_m*trk)*vr(2,3,j)+(z(k)-w_m*trk)*vr(3,3,j)

                   ! summation over modes for scaled velocity components
                   dusdt= 0.0_wp
                   dvsdt= 0.0_wp
                   dwsdt= 0.0_wp
                   do nm=1,Nmode
                      phase = xk1s(nm,j)*xs+xk2s(nm,j)*ys+xk3s(nm,j)*zs+omn(nm)*trk + psi(nm)
                      dsphase=(omn(nm)-xk1(nm)*u_m-xk2(nm)*v_m-xk3(nm)*w_m)*sin(phase)

                      dusdt= dusdt - dsphase*sigma1s(nm,j)
                      dvsdt= dvsdt - dsphase*sigma2s(nm,j)
                      dwsdt= dwsdt - dsphase*sigma3s(nm,j)
                   enddo

                   ! retrieve anisotropic field from scaled velocity field: u= R u_scaled
                   ut_in(j,k)= dusdt*vr(1,1,j) + dvsdt*vr(1,2,j) + dwsdt*vr(1,3,j)
                   vt_in(j,k)= dusdt*vr(2,1,j) + dvsdt*vr(2,2,j) + dwsdt*vr(2,3,j)
                   wt_in(j,k)= dusdt*vr(3,1,j) + dvsdt*vr(3,2,j) + dwsdt*vr(3,3,j)
                enddo
             enddo
          else
             do k=1,nz
                do j=ndy_c,nfy_c
                   ! scaled coordinates: x_scaled=R^T x
                   xs= (x(i)-u_m*trk)*vr(1,1,j)+(y(j)-v_m*trk)*vr(2,1,j)+(z(k)-w_m*trk)*vr(3,1,j)
                   ys= (x(i)-u_m*trk)*vr(1,2,j)+(y(j)-v_m*trk)*vr(2,2,j)+(z(k)-w_m*trk)*vr(3,2,j)
                   zs= (x(i)-u_m*trk)*vr(1,3,j)+(y(j)-v_m*trk)*vr(2,3,j)+(z(k)-w_m*trk)*vr(3,3,j)

                   ! summation over modes for scaled velocity components
                   dusdt= 0.0_wp
                   dvsdt= 0.0_wp
                   dwsdt= 0.0_wp
                   do nm=1,Nmode
                      phase = xk1s(nm,j)*xs + xk2s(nm,j)*ys + xk3s(nm,j)*zs + omn(nm)*trk + psi(nm)
                      dsphase=(omn(nm)-xk1(nm)*u_m-xk2(nm)*v_m-xk3(nm)*w_m)*sin(phase)

                      dusdt= dusdt - dsphase*sigma1s(nm,j)
                      dvsdt= dvsdt - dsphase*sigma2s(nm,j)
                      dwsdt= dwsdt - dsphase*sigma3s(nm,j)
                   enddo

                   ! retrieve anisotropic field from scaled velocity field: u= R u_scaled
                   ut_in(j,k)= dusdt*vr(1,1,j) + dvsdt*vr(1,2,j) + dwsdt*vr(1,3,j)
                   vt_in(j,k)= dusdt*vr(2,1,j) + dvsdt*vr(2,2,j) + dwsdt*vr(2,3,j)
                   wt_in(j,k)= dusdt*vr(3,1,j) + dvsdt*vr(3,2,j) + dwsdt*vr(3,3,j)
                enddo
             enddo
          endif
       else
          ut_in=0.0_wp
          vt_in=0.0_wp
          wt_in=0.0_wp
          if (is_curv3) then
             do k=1,nz
                do j=ndy_c,nfy_c
                   ! summation over modes
                   do nm=1,Nmode
                      phase = xk1(nm)*(xc3(i,j,k)-u_m*trk) &
                            + xk2(nm)*(yc3(i,j,k)-v_m*trk) &
                            + xk3(nm)*(zc3(i,j,k)-w_m*trk) + omn(nm)*trk + psi(nm)
                      dsphase=(omn(nm)-xk1(nm)*u_m-xk2(nm)*v_m-xk3(nm)*w_m)*sin(phase)

                      ut_in(j,k)= ut_in(j,k)-dsphase*sigma1(nm)
                      vt_in(j,k)= vt_in(j,k)-dsphase*sigma2(nm)
                      wt_in(j,k)= wt_in(j,k)-dsphase*sigma3(nm)
                   enddo
                enddo
             enddo
          elseif (is_curv) then
             do k=1,nz
                do j=ndy_c,nfy_c
                   ! summation over modes
                   do nm=1,Nmode
                      phase = xk1(nm)*(xc(i,j)-u_m*trk) &
                            + xk2(nm)*(yc(i,j)-v_m*trk) &
                            + xk3(nm)*(z(k)-w_m*trk) + omn(nm)*trk + psi(nm)
                      dsphase=(omn(nm)-xk1(nm)*u_m-xk2(nm)*v_m-xk3(nm)*w_m)*sin(phase)

                      ut_in(j,k)= ut_in(j,k)-dsphase*sigma1(nm)
                      vt_in(j,k)= vt_in(j,k)-dsphase*sigma2(nm)
                      wt_in(j,k)= wt_in(j,k)-dsphase*sigma3(nm)
                   enddo
                enddo
             enddo
          else
             do k=1,nz
                do j=ndy_c,nfy_c
                   ! summation over modes
                   do nm=1,Nmode
                      phase = xk1(nm)*(x(i)-u_m*trk) &
                            + xk2(nm)*(y(j)-v_m*trk) &
                            + xk3(nm)*(z(k)-w_m*trk) + omn(nm)*trk + psi(nm)
                      dsphase=(omn(nm)-xk1(nm)*u_m-xk2(nm)*v_m-xk3(nm)*w_m)*sin(phase)

                      ut_in(j,k)= ut_in(j,k) - dsphase*sigma1(nm)
                      vt_in(j,k)= vt_in(j,k) - dsphase*sigma2(nm)
                      wt_in(j,k)= wt_in(j,k) - dsphase*sigma3(nm)
                   enddo
                enddo
             enddo
          endif
       endif

       ! Apply weighting function
       ! =======================
       if ((anisotropy=='W').and.(coord(2)==0)) then
          ! add linear damping close to the wall
          do k=1,nz
             fy= (y(j)-y(1))/(y(6)-y(1))
             do j=ndy_c,5
                ut_in(j,k)= ut_in(j,k)*u_rms(j)*fy
                vt_in(j,k)= vt_in(j,k)*v_rms(j)*fy
                wt_in(j,k)= wt_in(j,k)*w_rms(j)*fy
             enddo
          enddo
          ! multiply by rms profiles
          do k=1,nz
             do j=ndy_c,nfy_c
                ut_in(j,k)= ut_in(j,k)*u_rms(j)
                vt_in(j,k)= vt_in(j,k)*v_rms(j)
                wt_in(j,k)= wt_in(j,k)*w_rms(j)
             enddo
          enddo
       elseif (anisotropy=='W') then
          ! multiply by rms profiles
          do k=1,nz
             do j=ndy_c,nfy_c
                ut_in(j,k)= ut_in(j,k)*u_rms(j)
                vt_in(j,k)= vt_in(j,k)*v_rms(j)
                wt_in(j,k)= wt_in(j,k)*w_rms(j)
             enddo
          enddo
       endif
       
       ! Apply damping coefficients
       ! ==========================
       if (is_curv3) then
          do k=1,nz
             do j=ndy_c,nfy_c
                ut_in(j,k)= ut_in(j,k)*damping_coeff3(j,k)
                vt_in(j,k)= vt_in(j,k)*damping_coeff3(j,k)
                wt_in(j,k)= wt_in(j,k)*damping_coeff3(j,k)
             enddo
          enddo
       else
          do k=1,nz
             do j=ndy_c,nfy_c
                ut_in(j,k)= ut_in(j,k)*damping_coeff(j)
                vt_in(j,k)= vt_in(j,k)*damping_coeff(j)
                wt_in(j,k)= wt_in(j,k)*damping_coeff(j)
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
             do k=1,nz
                do j=ndeb_j2(jbc),nend_j2(jbc)
                   ! scaled coordinates: x_scaled=R^T x
                   xs= (x(i)-u_m_TBL*trk)*vr_TBL(1,1,j,jbc) &
                      +(y(j)-v_m_TBL*trk)*vr_TBL(2,1,j,jbc) &
                      +(z(k)-w_m_TBL*trk)*vr_TBL(3,1,j,jbc)
                   ys= (x(i)-u_m_TBL*trk)*vr_TBL(1,2,j,jbc) &
                      +(y(j)-v_m_TBL*trk)*vr_TBL(2,2,j,jbc) &
                      +(z(k)-w_m_TBL*trk)*vr_TBL(3,2,j,jbc)
                   zs= (x(i)-u_m_TBL*trk)*vr_TBL(1,3,j,jbc) &
                      +(y(j)-v_m_TBL*trk)*vr_TBL(2,3,j,jbc) &
                      +(z(k)-w_m_TBL*trk)*vr_TBL(3,3,j,jbc)

                   ! summation over modes for scaled velocity components
                   dusdt= 0.0_wp
                   dvsdt= 0.0_wp
                   dwsdt= 0.0_wp
                   do nm=1,Nmode
                      phase= xk1s_TBL(nm,j,jbc)*xs&
                            +xk2s_TBL(nm,j,jbc)*ys&
                            +xk3s_TBL(nm,j,jbc)*zs+psi_TBL(nm,jbc)+omn_TBL(nm,jbc)*trk
                      dsphase=(omn_TBL(nm,jbc)-xk1s_TBL(nm,j,jbc)*u_m &
                                              -xk2s_TBL(nm,j,jbc)*v_m &
                                              -xk3s_TBL(nm,j,jbc)*w_m)*sin(phase)
!!$                      dsphase=cos(phase)

                      dusdt= dusdt - dsphase*sigma1s_TBL(nm,j,jbc)
                      dvsdt= dvsdt - dsphase*sigma2s_TBL(nm,j,jbc)
                      dwsdt= dwsdt - dsphase*sigma3s_TBL(nm,j,jbc)
                   enddo

                   ! retrieve anisotropic field from scaled velocity field: u= R u_scaled
                   ut_in(j,k)= dusdt*vr_TBL(1,1,j,jbc)+dvsdt*vr_TBL(1,2,j,jbc)+dwsdt*vr_TBL(1,3,j,jbc)
                   vt_in(j,k)= dusdt*vr_TBL(2,1,j,jbc)+dvsdt*vr_TBL(2,2,j,jbc)+dwsdt*vr_TBL(2,3,j,jbc)
                   wt_in(j,k)= dusdt*vr_TBL(3,1,j,jbc)+dvsdt*vr_TBL(3,2,j,jbc)+dwsdt*vr_TBL(3,3,j,jbc)
                enddo
             enddo
          endif
       endif
    enddo
    
  end subroutine disturb_inlet_RFM_charac

  !===============================================================================
  module subroutine disturb_inlet_RFM_turb(u_in,v_in,w_in)
  !===============================================================================
    !> Impose RFM disturbances at inlet (for turb BC)
  !===============================================================================
    use mod_time
    implicit none
    ! ----------------------------------------------------------------------------
    ! RFM disturbances
    real(wp), dimension(ny,nz), intent(out) :: u_in,v_in,w_in
    ! ----------------------------------------------------------------------------
    ! local var
    integer :: i,j,k,nm
    real(wp) :: trk
    real(wp) :: phase,cphase
    ! ----------------------------------------------------------------------------

    i= 1
    trk= time + ck(irk)*deltat

    ! Construct turbulent stochastic field
    ! ------------------------------------
    if (anisotropy=='S') then
       call mpistop('Anisotropy not implemented for turb BC in disturb_inlet_RFM_turb (mod_RFM.f90)!',1)
    else
       u_in=0.0_wp
       v_in=0.0_wp
       w_in=0.0_wp
       if (is_curv3) then
          do k=1,nz
             do j=ndy_c,nfy_c
                ! summation over modes
                do nm=1,Nmode
                   phase= xk1(nm)*(xc3(i,j,k)-u_m*trk) + xk2(nm)*(yc3(i,j,k)-v_m*trk) + xk3(nm)*(zc3(i,j,k)-w_m*trk) + omn(nm)*trk + psi(nm)
                   cphase=cos(phase)

                   u_in(j,k)= u_in(j,k) + cphase*sigma1(nm)
                   v_in(j,k)= v_in(j,k) + cphase*sigma2(nm)
                   w_in(j,k)= w_in(j,k) + cphase*sigma3(nm)
                enddo
             enddo
          enddo
       elseif (is_curv) then
          do k=1,nz
             do j=ndy_c,nfy_c
                ! summation over modes
                do nm=1,Nmode
                   phase= xk1(nm)*(xc(i,j)-u_m*trk) + xk2(nm)*(yc(i,j)-v_m*trk) + xk3(nm)*(z(k)-w_m*trk) + omn(nm)*trk + psi(nm)
                   cphase=cos(phase)

                   u_in(j,k)= u_in(j,k) + cphase*sigma1(nm)
                   v_in(j,k)= v_in(j,k) + cphase*sigma2(nm)
                   w_in(j,k)= w_in(j,k) + cphase*sigma3(nm)
                enddo
             enddo
          enddo
       else
          do k=1,nz
             do j=ndy_c,nfy_c
                ! summation over modes
                do nm=1,Nmode
                   phase=xk1(nm)*(x(i)-u_m*trk) + xk2(nm)*(y(j)-v_m*trk) + xk3(nm)*(z(k)-w_m*trk) + omn(nm)*trk + psi(nm)
                   cphase=cos(phase)

                   u_in(j,k)= u_in(j,k) + cphase*sigma1(nm)
                   v_in(j,k)= v_in(j,k) + cphase*sigma2(nm)
                   w_in(j,k)= w_in(j,k) + cphase*sigma3(nm)
                enddo
             enddo
          enddo
       endif
    endif

    ! Apply damping coefficients
    ! ==========================
    if (is_curv3) then
       do k=1,nz
          do j=ndy_c,nfy_c
             u_in(j,k)= u_in(j,k)*damping_coeff3(j,k)
             v_in(j,k)= v_in(j,k)*damping_coeff3(j,k)
             w_in(j,k)= w_in(j,k)*damping_coeff3(j,k)
          enddo
       enddo
    else
       do k=1,nz
          do j=ndy_c,nfy_c
             u_in(j,k)= u_in(j,k)*damping_coeff(j)
             v_in(j,k)= v_in(j,k)*damping_coeff(j)
             w_in(j,k)= w_in(j,k)*damping_coeff(j)
          enddo
       enddo
    endif

  end subroutine disturb_inlet_RFM_turb

  ! !===============================================================================
  ! module subroutine disturb_inlet_RFM_Charac(dpx,ut_in,vt_in,wt_in,pt_in,rt_in)
  ! !===============================================================================
  !   !> Impose RFM disturbances at inlet (for charac BC)
  ! !===============================================================================
  !   use mod_time
  !   implicit none
  !   ! ----------------------------------------------------------------------------
  !   integer :: i,j,k
  !   ! RFM disturbances
  !   real(wp), dimension(ny,nz), intent(in)  :: dpx
  !   real(wp), dimension(ny,nz), intent(out) :: ut_in,vt_in,wt_in,pt_in,rt_in
  !   ! ----------------------------------------------------------------------------
  !   ! local var
  !   integer  :: nm
  !   real(wp) :: trk,fy,L1
  !   real(wp) :: phase,sphase,cphase
  !   ! real(wp) :: ut,vt,wt
  !   real(wp) :: u_in,dutdx
  !   real(wp) :: xs,ys,zs ! scaled cooordinates
  !   ! scaled velocity components
  !   real(wp) :: dusdt,dvsdt,dwsdt
  !   ! ----------------------------------------------------------------------------

  !   i= 1
  !   trk= time + ck(irk)*deltat

  !   ! Construct turbulent stochastic field
  !   ! ------------------------------------
  !   if (anisotropy=='S') then
  !      do k=1,nz
  !         do j=ndy_c,nfy_c
  !            ! scaled coordinates: x_scaled=R^T x
  !            !xs= (x(i)-0.*u_m*t)*vr(1,1,j)+y(j)*vr(2,1,j)+z(k)*vr(3,1,j)
  !            !ys= (x(i)-0.*u_m*t)*vr(1,2,j)+y(j)*vr(2,2,j)+z(k)*vr(3,2,j)
  !            !zs= (x(i)-0.*u_m*t)*vr(1,3,j)+y(j)*vr(2,3,j)+z(k)*vr(3,3,j)
  !            xs= x(i)*vr(1,1,j)+y(j)*vr(2,1,j)+z(k)*vr(3,1,j)
  !            ys= x(i)*vr(1,2,j)+y(j)*vr(2,2,j)+z(k)*vr(3,2,j)
  !            zs= x(i)*vr(1,3,j)+y(j)*vr(2,3,j)+z(k)*vr(3,3,j)
  !            ! summation over modes for scaled velocity components
  !            dusdt=0.0_wp
  !            dvsdt=0.0_wp
  !            dwsdt=0.0_wp
  !            do nm=1,Nmode
  !               phase=xk1s(nm,j)*xs+xk2s(nm,j)*ys+xk3s(nm,j)*zs + omn(nm)*trk + psi(nm)
  !               sphase=sin(phase)

  !               dusdt=dusdt-omn(nm)*sphase*sigma1(nm)
  !               dvsdt=dvsdt-omn(nm)*sphase*sigma2(nm)
  !               dwsdt=dwsdt-omn(nm)*sphase*sigma3(nm)
  !            enddo
  !            ! retrieve anisotropic field from scaled velocity field: u= R u_scaled
  !            ut_in(j,k)= dusdt*vr(1,1,j) + dvsdt*vr(1,2,j) + dwsdt*vr(1,3,j)
  !            vt_in(j,k)= dusdt*vr(2,1,j) + dvsdt*vr(2,2,j) + dwsdt*vr(2,3,j)
  !            wt_in(j,k)= dusdt*vr(3,1,j) + dvsdt*vr(3,2,j) + dwsdt*vr(3,3,j)
  !         enddo
  !      enddo
  !   else
  !      ut_in=0.0_wp
  !      vt_in=0.0_wp
  !      wt_in=0.0_wp
  !      do k=1,nz
  !         do j=ndy_c,nfy_c
  !            L1= 0.0_wp
  !            u_in= 0.0_wp
  !            dutdx= 0.0_wp
  !            ! summation over modes
  !            do nm=1,Nmode
  !               phase=xk1(nm)*(x(i)-u_m*trk)+xk2(nm)*y(j)+xk3(nm)*z(k) + omn(nm)*trk + psi(nm)
  !               ! phase=xk1(nm)*(x(i)-u_ref*trk)+xk2(nm)*y(j)+xk3(nm)*z(k) + omn(nm)*trk + psi(nm)
  !               cphase=cos(phase)
  !               sphase=sin(phase)


  !               ! Compute sum(u_pertu), sum(d(u_pertu)/dx)
  !               u_in= u_in + cphase*sigma1(nm)
  !               dutdx= dutdx - xk1(nm)*sphase*sigma1(nm)

  !               ut_in(j,k)= ut_in(j,k) - (omn(nm)-xk1(nm)*u_ref)*sphase*sigma1(nm)
  !               vt_in(j,k)= vt_in(j,k) - (omn(nm)-xk1(nm)*u_ref)*sphase*sigma2(nm)
  !               wt_in(j,k)= wt_in(j,k) - (omn(nm)-xk1(nm)*u_ref)*sphase*sigma3(nm)
  !            enddo

  !            ! Compute characteristic amplitudes
  !            ! =================================
  !            ! L1= (uu-u_pertu - c_)*(dpx - rho*c_*(dux - d(u_pertu)/dx))     (d(p_pertu)/dx unknown)
  !            L1= (uu(i,j,k)-u_in-c_(i,j,k))*(dpx(j,k)-rho_n(i,j,k)*c_(i,j,k)*(dux(i,j,k)-dutdx))

  !            ! Time derivatives for primitive variables
  !            ! ========================================
  !            ut_in(j,k)= ut_in(j,k) + L1/(2*rho_n(i,j,k)*c_(i,j,k))
  !            pt_in(j,k)= -0.5_wp*L1
  !            rt_in(j,k)= 0.5_wp*L1/(c_(i,j,k)**2)
  !         enddo
  !      enddo
  !   endif

  !   ! Apply weighting function
  !   ! ========================
  !   if ((anisotropy=='W').and.(coord(2)==0)) then
  !      ! add linear damping close to the wall
  !      do k=1,nz
  !         fy= (y(j)-y(1))/(y(6)-y(1))
  !         do j=ndy_c,5
  !            ut_in(j,k)= ut_in(j,k)*u_rms(j)*fy
  !            vt_in(j,k)= vt_in(j,k)*v_rms(j)*fy
  !            wt_in(j,k)= wt_in(j,k)*w_rms(j)*fy
  !         enddo
  !      enddo
  !      ! multiply by rms profiles
  !      do k=1,nz
  !         do j=ndy_c,nfy_c
  !            ut_in(j,k)= ut_in(j,k)*u_rms(j)
  !            vt_in(j,k)= vt_in(j,k)*v_rms(j)
  !            wt_in(j,k)= wt_in(j,k)*w_rms(j)
  !         enddo
  !      enddo
  !   else if (anisotropy=='W') then
  !      ! multiply by rms profiles
  !      do k=1,nz
  !         do j=ndy_c,nfy_c
  !            ut_in(j,k)= ut_in(j,k)*u_rms(j)
  !            vt_in(j,k)= vt_in(j,k)*v_rms(j)
  !            wt_in(j,k)= wt_in(j,k)*w_rms(j)
  !         enddo
  !      enddo
  !   endif

  !   ! Apply damping coefficients
  !   ! ==========================
  !   do k=1,nz
  !      do j=ndy_c,nfy_c
  !         ut_in(j,k)= ut_in(j,k)*damping_coeff(j)
  !         vt_in(j,k)= vt_in(j,k)*damping_coeff(j)
  !         wt_in(j,k)= wt_in(j,k)*damping_coeff(j)
  !      enddo
  !   enddo

  ! end subroutine disturb_inlet_RFM_Charac

end submodule smod_RFM_charac
