!===============================================================================
submodule (mod_RFM) smod_RFM_init_vel
!===============================================================================
  !> author: AB & XG
  !> date: February 2024
  !> Generation of Random Fourier Modes (RFM)
  !> Initialization of a field with RFM & check RFM routine
!=============================================================================== 

contains

  !===============================================================================
  module subroutine init_vel_RFM
  !===============================================================================
    !> Initialization of a velocity field with RFM modes
  !===============================================================================
    use mod_time
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k,nm,nys,jbc
    real(wp) :: phase,trk
    real(wp) :: xs,ys,zs ! scaled cooordinates
    real(wp) :: us,vs,ws ! scaled velocity components
    ! ----------------------------------------------------------------------------
    real(wp), dimension(:,:,:), allocatable :: uu_temp,vv_temp,ww_temp
    ! ----------------------------------------------------------------------------

    if (iproc.eq.0) print *,"~> Initialisation of field with RFM"

    ! RFM for initialization only on the first processor (except if CHAN)
    ! --------------------------------------------------
    if ((coord(1).ne.0).and.(.not.CHAN)) return

    ! Allocation & Initialization
    ! ---------------------------
    allocate(uu_temp(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(vv_temp(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(ww_temp(nx1:nx2,ny1:ny2,nz1:nz2))
    uu_temp=0.0_wp
    vv_temp=0.0_wp
    ww_temp=0.0_wp

    ! -------------------
    ! I/ Injection of FST
    ! -------------------
    if (is_RFM_FST) then
       ! Construct turbulent stochastic field
       ! ------------------------------------
       if (anisotropy=='S') then
          call mpistop('init_vel_RFM in anisotropy needs to be updated ! ',1)
          do i=1,nx
             trk = -(x(i)-xg_in)/u_ref + ntotal*deltat
             do j=1,ny
                do k=1,nz
                   ! scaled coordinates: x_scaled=R^T x
                   if (is_curv) then
                      xs=xc(i,j)*vr(1,1,j)+yc(i,j)*vr(2,1,j)+z(k)*vr(3,1,j)
                      ys=xc(i,j)*vr(1,2,j)+yc(i,j)*vr(2,2,j)+z(k)*vr(3,2,j)
                      zs=xc(i,j)*vr(1,3,j)+yc(i,j)*vr(2,3,j)+z(k)*vr(3,3,j)
                   else
                      xs=x(i)*vr(1,1,j)+y(j)*vr(2,1,j)+z(k)*vr(3,1,j)
                      ys=x(i)*vr(1,2,j)+y(j)*vr(2,2,j)+z(k)*vr(3,2,j)
                      zs=x(i)*vr(1,3,j)+y(j)*vr(2,3,j)+z(k)*vr(3,3,j)
                   endif
                   ! summation over modes for scaled velocity components
                   us=0.0_wp
                   vs=0.0_wp
                   ws=0.0_wp
                   do nm=1,Nmode
                      phase = xk1s(nm,j)*xs + xk2s(nm,j)*ys + xk3s(nm,j)*zs + omn(nm)*trk + psi(nm)
                      phase=cos(phase)
                      us=us+phase*sigma1s(nm,j)
                      vs=vs+phase*sigma2s(nm,j)
                      ws=ws+phase*sigma3s(nm,j)
                   enddo
                   ! retrieve anisotropic field from scaled velocity field: u= R u_scaled
                   uu_temp(i,j,k) = us*vr(1,1,j) + vs*vr(1,2,j) + ws*vr(1,3,j)
                   vv_temp(i,j,k) = us*vr(2,1,j) + vs*vr(2,2,j) + ws*vr(2,3,j)
                   ww_temp(i,j,k) = us*vr(3,1,j) + vs*vr(3,2,j) + ws*vr(3,3,j)
                enddo
             enddo
          enddo
       else
          if (is_curv3) then
             do i=1,nx
                do j=1,ny
                   do k=1,nz
                      trk = -(xc3(i,j,k)-xg_in)/u_m + ntotal*deltat
                      ! summation over modes
                      do nm=1,Nmode
                         phase = xk1(nm)*(xc3(i,j,k)-u_m*ntotal*deltat) &
                               + xk2(nm)*(yc3(i,j,k)-v_m*ntotal*deltat) &
                               + xk3(nm)*(zc3(i,j,k)-w_m*ntotal*deltat) + omn(nm)*trk + psi(nm)
                         phase = cos(phase)

                         uu_temp(i,j,k) = uu_temp(i,j,k) + phase*sigma1(nm)
                         vv_temp(i,j,k) = vv_temp(i,j,k) + phase*sigma2(nm)
                         ww_temp(i,j,k) = ww_temp(i,j,k) + phase*sigma3(nm)
                      enddo
                   enddo
                enddo
             enddo
          else if (is_curv) then
             do i=1,nx
                do j=1,ny
                   do k=1,nz
                      trk = -(xc(i,j)-xg_in)/u_m + ntotal*deltat
                      ! summation over modes
                      do nm=1,Nmode
                         phase = xk1(nm)*(xc(i,j)-u_m*ntotal*deltat) &
                               + xk2(nm)*(yc(i,j)-v_m*ntotal*deltat) &
                               + xk3(nm)*(z(k)-w_m*ntotal*deltat) + omn(nm)*trk + psi(nm)
                         phase = cos(phase)

                         uu_temp(i,j,k) = uu_temp(i,j,k) + phase*sigma1(nm)
                         vv_temp(i,j,k) = vv_temp(i,j,k) + phase*sigma2(nm)
                         ww_temp(i,j,k) = ww_temp(i,j,k) + phase*sigma3(nm)
                      enddo
                   enddo
                enddo
             enddo
          else
             do i=1,nx
                do j=1,ny
                   do k=1,nz
                      trk = -(x(i)-xg_in)/u_m + ntotal*deltat
                      ! summation over modes
                      do nm=1,Nmode
                         phase = xk1(nm)*(x(i)-u_m*ntotal*deltat) &
                               + xk2(nm)*(y(j)-v_m*ntotal*deltat) &
                               + xk3(nm)*(z(k)-w_m*ntotal*deltat) + omn(nm)*trk + psi(nm)
                         phase = cos(phase)

                         uu_temp(i,j,k) = uu_temp(i,j,k) + phase*sigma1(nm)
                         vv_temp(i,j,k) = vv_temp(i,j,k) + phase*sigma2(nm)
                         ww_temp(i,j,k) = ww_temp(i,j,k) + phase*sigma3(nm)
                      enddo
                   enddo
                enddo
             enddo
          endif
       endif

       ! Apply weighting function
       ! ========================
       if (anisotropy=='W') then
          if (coord(2)==0) then
             nys=21
             ! linear damping close to the wall (20 first points)
             do i=1,nx
                do j=1,20
                   do k=1,nz
                      uu_temp(i,j,k)=uu_temp(i,j,k)*u_rms(j)*(y(j)-y(1))/(y(20)-y(1))
                      vv_temp(i,j,k)=vv_temp(i,j,k)*v_rms(j)*(y(j)-y(1))/(y(20)-y(1))
                      ww_temp(i,j,k)=ww_temp(i,j,k)*w_rms(j)*(y(j)-y(1))/(y(20)-y(1))
                   enddo
                enddo
             enddo
          else
             nys=1
          endif
          ! multiply by rms profiles
          do i=1,nx
             do j=nys,ny
                do k=1,nz
                   uu_temp(i,j,k)=uu_temp(i,j,k)*u_rms(j)
                   vv_temp(i,j,k)=vv_temp(i,j,k)*v_rms(j)
                   ww_temp(i,j,k)=ww_temp(i,j,k)*w_rms(j)
                enddo
             enddo
          enddo
       endif

       ! Apply damping coefficients
       ! ==========================
       if (is_curv3) then
          do i=1,nx
             do j=ndy_d,nfy_d
                do k=1,nz
                   uu_temp(i,j,k) = uu_temp(i,j,k)*damping_coeff3(j,k)
                   vv_temp(i,j,k) = vv_temp(i,j,k)*damping_coeff3(j,k)
                   ww_temp(i,j,k) = ww_temp(i,j,k)*damping_coeff3(j,k)
                enddo
             enddo
          enddo
       else
          do i=1,nx
             do j=ndy_d,nfy_d
                do k=1,nz
                   uu_temp(i,j,k) = uu_temp(i,j,k)*damping_coeff(j)
                   vv_temp(i,j,k) = vv_temp(i,j,k)*damping_coeff(j)
                   ww_temp(i,j,k) = ww_temp(i,j,k)*damping_coeff(j)
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
             call mpistop("init_vel_RFM in is_curv3 for TBL not implemented yet ",1)
          else if (is_curv) then
             call mpistop("init_vel_RFM in is_curv for TBL not implemented yet ",1)
          else
             do i=1,nx
                trk = -(x(i)-xg_in)/u_ref + ntotal*deltat ! XG WHY u_ref?? is it um +vm,wm
                do k=1,nz
                   do j=1,ny
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
                      us=0.0_wp
                      vs=0.0_wp
                      ws=0.0_wp
                      do nm=1,Nmode
                         phase= xk1s_TBL(nm,j,jbc)*xs&
                               +xk2s_TBL(nm,j,jbc)*ys&
                               +xk3s_TBL(nm,j,jbc)*zs+psi_TBL(nm,jbc)+omn_TBL(nm,jbc)*trk
                         phase=cos(phase)
                         us=us+phase*sigma1s_TBL(nm,j,jbc)
                         vs=vs+phase*sigma2s_TBL(nm,j,jbc)
                         ws=ws+phase*sigma3s_TBL(nm,j,jbc)
                      enddo
                      ! retrieve anisotropic field from scaled velocity field: u= R u_scaled
                      uu_temp(i,j,k)= us*vr_TBL(1,1,j,jbc)+vs*vr_TBL(1,2,j,jbc)+ws*vr_TBL(1,3,j,jbc)
                      vv_temp(i,j,k)= us*vr_TBL(2,1,j,jbc)+vs*vr_TBL(2,2,j,jbc)+ws*vr_TBL(2,3,j,jbc)
                      ww_temp(i,j,k)= us*vr_TBL(3,1,j,jbc)+vs*vr_TBL(3,2,j,jbc)+ws*vr_TBL(3,3,j,jbc)
                   enddo
                enddo
             enddo
             !call mpistop("init_vel_RFM in is_cartesian for TBL not implemented yet ",1)
          endif
       endif
    enddo

    ! Add mean streamwise velocity field
    ! ==================================
    do i=1,nx
       do j=ndy_d,nfy_d
          do k=1,nz
             uu(i,j,k) = uu(i,j,k) + uu_temp(i,j,k)
             vv(i,j,k) = vv(i,j,k) + vv_temp(i,j,k)
             ww(i,j,k) = ww(i,j,k) + ww_temp(i,j,k)
          enddo
       enddo
    enddo

    if (is_2D) ww=0.0_wp

    deallocate(uu_temp,vv_temp,ww_temp)

    ! Update of variables
    ! ===================
    rhou = rho*uu
    rhov = rho*vv
    rhow = rho*ww
    do k=ndzt,nfzt
       do j=ndyt,nfyt
          do i=ndxt,nfxt
             rhoe(i,j,k)= rho(i,j,k)*(ecalc_tro(Tmp(i,j,k),rho(i,j,k)) &
                        + 0.5_wp*(uu(i,j,k)**2+vv(i,j,k)**2+ww(i,j,k)**2))
          enddo
       enddo
    enddo

  end subroutine init_vel_RFM

  !===============================================================================
  module subroutine compute_RFM_check
  !===============================================================================
    !> compute stochastic turbulent field to check convergence of rms profiles
  !===============================================================================
    use mod_time
    use mod_mpi_part
    use mod_io_snapshots
    use mod_block
    use mod_utils
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k,nm,nys,wr
    real(wp) :: t,phase
    real(wp) :: xs,ys,zs ! scaled cooordinates
    real(wp) :: us,vs,ws ! scaled velocity components
    real(wp), dimension(nob_inp,ny,nz) :: ut,vt,wt
    ! ----------------------------------------------------------------------------
    ! compute rms profiles
    real(wp) :: nr,nr1
    real(wp), dimension(nob_inp,ny,nz) :: urt,vrt,wrt,uvrt,uwrt,vwrt
    real(wp), dimension(ny) :: urms,vrms,wrms,uvrms,uwrms,vwrms
    ! ----------------------------------------------------------------------------

    ! Initialization
    ! --------------
    ut = 0.0_wp; vt = 0.0_wp; wt = 0.0_wp

    ! Check convergence
    ! -----------------
    if (iproc==0) then
       open(49,file='turbul.bin',form='unformatted',status='unknown')
       rewind(49)
       write(49) ngy
       write(49) ngz
       write(49) (yg(j),j=1,ngy)
       write(49) (zg(k),k=1,ngz)
    endif

    ! local timestep
    ! --------------
    ! CFL=0.7
    t=0.0_wp
    deltat=CFL*deltay/(c_ref+U_ref)*1.0_wp ! TO BE CHECKED *100 to accelerate mean
    deltat=3.697607637574558E-008 ! TO BE CHECKED *100 to accelerate mean
    if (iproc.eq.0) print *,3.697607637574558E-008,deltat

    ! initialization
    ! --------------
    urt=0.0_wp
    vrt=0.0_wp
    wrt=0.0_wp
    uvrt=0.0_wp
    uwrt=0.0_wp
    vwrt=0.0_wp
    ntime=0.0_wp

    ! initialisation of the planes
    ! ----------------------------
    if (.not.is_2D) call init_RFM_planes

    ! Temporal loop (to test convergence of rms fluctuations)
    ! -------------
    do while (ntime<nmax) ! [use nmax in param.ini file]

       ! Temporal increment
       ! ------------------
       ntime=ntime+1
       t=t+deltat
       if ((iproc==0).and.(mod(ntime,100)==0)) then
          write(6,*) 'Iteration number:', ntime, 'Time:',t
       endif

       ! Construct turbulent stochastic field
       ! ------------------------------------
       if (anisotropy=='S') then
          if ((is_curv3).or.(is_curv)) call mpistop("compute_RFM_check not implemented for anisotropic turbulence in curvilinear grid",1)
          do j=1,ny
             do k=1,nz
                do i=1,nob_inp
                   ! scaled coordinates: x_scaled=R^T x
                   xs = (x(i)-u_m*t)*vr(1,1,j) + (y(j)-v_m*t)*vr(2,1,j) + (z(k)-w_m*t)*vr(3,1,j)
                   ys = (x(i)-u_m*t)*vr(1,2,j) + (y(j)-v_m*t)*vr(2,2,j) + (z(k)-w_m*t)*vr(3,2,j)
                   zs = (x(i)-u_m*t)*vr(1,3,j) + (y(j)-v_m*t)*vr(2,3,j) + (z(k)-w_m*t)*vr(3,3,j)

                   ! summation over modes for scaled velocity components
                   us=0.0_wp
                   vs=0.0_wp
                   ws=0.0_wp
                   do nm=1,Nmode
                      phase=xk1s(nm,j)*xs+xk2s(nm,j)*ys+xk3s(nm,j)*zs + omn(nm)*t + psi(nm)
                      phase=cos(phase)
                      us=us+phase*sigma1s(nm,j)
                      vs=vs+phase*sigma2s(nm,j)
                      ws=ws+phase*sigma3s(nm,j)
                   enddo
                   ! retrieve anisotropic field from scaled velocity field: u= R u_scaled
                   ut(i,j,k)=us*vr(1,1,j)+vs*vr(1,2,j)+ws*vr(1,3,j)
                   vt(i,j,k)=us*vr(2,1,j)+vs*vr(2,2,j)+ws*vr(2,3,j)
                   wt(i,j,k)=us*vr(3,1,j)+vs*vr(3,2,j)+ws*vr(3,3,j)
                enddo
             enddo
          enddo
       else
          if (is_curv3) then
             do k=1,nz
                do j=1,ny
                   do i=1,nob_inp
                      ! summation over modes
                      do nm=1,Nmode
                         phase = xk1(nm)*(xc3(i,j,k)-u_m*t) + xk2(nm)*(yc3(i,j,k)-v_m*t) + xk3(nm)*(zc3(i,j,k)-w_m*t) + omn(nm)*t + psi(nm)
                         phase = cos(phase)

                         ut(i,j,k) = ut(i,j,k) + phase*sigma1(nm)
                         vt(i,j,k) = vt(i,j,k) + phase*sigma2(nm)
                         wt(i,j,k) = wt(i,j,k) + phase*sigma3(nm)
                      enddo
                   enddo
                enddo
             enddo
          else if (is_curv) then
             do k=1,nz
                do j=1,ny
                   do i=1,nob_inp
                      ! summation over modes
                      do nm=1,Nmode
                         phase= xk1(nm)*(xc(i,j)-u_m*t)+xk2(nm)*(yc(i,j)-v_m*t)+xk3(nm)*(z(k)-w_m*t)+omn(nm)*t+psi(nm)
                         phase= cos(phase)

                         ut(i,j,k) = ut(i,j,k) + phase*sigma1(nm)
                         vt(i,j,k) = vt(i,j,k) + phase*sigma2(nm)
                         wt(i,j,k) = wt(i,j,k) + phase*sigma3(nm)
                      enddo
                   enddo
                enddo
             enddo
          else
             do k=1,nz
                do j=1,ny
                   do i=1,nob_inp
                      ! summation over modes
                      do nm=1,Nmode
                         phase=xk1(nm)*(x(i)-u_m*t)+xk2(nm)*(y(j)-v_m*t)+xk3(nm)*(z(k)-w_m*t)+omn(nm)*t+psi(nm)
                         phase = cos(phase)

                         ut(i,j,k) = ut(i,j,k) + phase*sigma1(nm)
                         vt(i,j,k) = vt(i,j,k) + phase*sigma2(nm)
                         wt(i,j,k) = wt(i,j,k) + phase*sigma3(nm)
                      enddo
                   enddo
                enddo
             enddo
          endif
       endif

       ! Apply weighting function
       ! ========================
       if (anisotropy=='W') then
          if (coord(2)==0) then
             nys=21
             ! linear damping close to the wall (20 first points)
             do i=1,nob_inp
                do j=1,nys-1
                   ut(i,j,:) = ut(i,j,:)*u_rms(j)*(y(j)-y(1))/(y(20)-y(1))
                   vt(i,j,:) = vt(i,j,:)*v_rms(j)*(y(j)-y(1))/(y(20)-y(1))
                   wt(i,j,:) = wt(i,j,:)*w_rms(j)*(y(j)-y(1))/(y(20)-y(1))
                enddo
             enddo
          else
             nys=1
          endif
          ! multiply by rms profiles
          do j=nys,ny
             ut(:,j,:)=ut(:,j,:)*u_rms(j)
             vt(:,j,:)=vt(:,j,:)*v_rms(j)
             wt(:,j,:)=wt(:,j,:)*w_rms(j)
          enddo
       endif

       ! Apply damping coefficients
       ! ==========================
       if (is_curv3) then
          do k=1,nz
             do j=1,ny
                do i=1,nob_inp
                   ut(i,j,k) = ut(i,j,k)*damping_coeff3(j,k)
                   vt(i,j,k) = vt(i,j,k)*damping_coeff3(j,k)
                   wt(i,j,k) = wt(i,j,k)*damping_coeff3(j,k)
                enddo
             enddo
          enddo
       else
          do k=1,nz
             do j=1,ny
                do i=1,nob_inp
                   ut(i,j,k) = ut(i,j,k)*damping_coeff(j)
                   vt(i,j,k) = vt(i,j,k)*damping_coeff(j)
                   wt(i,j,k) = wt(i,j,k)*damping_coeff(j)
                enddo
             enddo
          enddo
       endif

       ! Write fluctuations in file (Att. not MPI)
       ! --------------------------
       if (nproc==1) then
          if (ntime==nmax.or.ntime==nmax-nmax/4.or.ntime==nmax-nmax/2.or.ntime==nmax-3*nmax/4) then
             write(6,*) '~> write binary file at iteration ', ntime
             write(49) ((ut(1,j,k),j=1,ny),k=1,nz)
             write(49) ((vt(1,j,k),j=1,ny),k=1,nz)
             write(49) ((wt(1,j,k),j=1,ny),k=1,nz)
          endif
       endif

       ! Write fluctuations of inlet planes
       ! ----------------------------------
       if ((mod(ntime,freq_plane).eq.0).and.(.not.is_2D)) then
         wr = 0
          do i=1,nob_inp
             if (inlet_ip(i).ne.0)  then
                if (wr.eq.0) then
                   if (iproc==0) write(6,*) '~> write inlet planes files at iteration ', ntime
                   uvar(1:nob_inp,:,:,1) = ut(:,:,:)
                   uvar(1:nob_inp,:,:,2) = vt(:,:,:)
                   uvar(1:nob_inp,:,:,3) = wt(:,:,:)
                   wr = 1
                endif
                call write_snapshot(inlet_ip(i))
             endif
          end do
      endif



       ! Temporal mean of Reynolds stresses
       ! ----------------------------------
       nr =dble(ntime)
       nr1=dble(ntime-1)
       urt =(nr1*urt +ut*ut)/nr
       vrt =(nr1*vrt +vt*vt)/nr
       wrt =(nr1*wrt +wt*wt)/nr
       uvrt=(nr1*uvrt+ut*vt)/nr
       uwrt=(nr1*uwrt+ut*wt)/nr
       vwrt=(nr1*vwrt+vt*wt)/nr

    enddo

    ! Spanwise mean of Reynolds stresses
    ! ==================================
    if (iproc==0) print *,"Calculation of spanwise Reynold stresses"
    urms=0.0_wp
    vrms=0.0_wp
    wrms=0.0_wp
    uvrms=0.0_wp
    uwrms=0.0_wp
    vwrms=0.0_wp
    do k=1,nz
       urms=urms+urt(1,:,k)
       vrms=vrms+vrt(1,:,k)
       wrms=wrms+wrt(1,:,k)
       uvrms=uvrms+uvrt(1,:,k)
       uwrms=uwrms+uwrt(1,:,k)
       vwrms=vwrms+vwrt(1,:,k)
    enddo
    call MPI_ALLREDUCE(MPI_IN_PLACE, urms,ny,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXY,info)
    call MPI_ALLREDUCE(MPI_IN_PLACE, vrms,ny,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXY,info)
    call MPI_ALLREDUCE(MPI_IN_PLACE, wrms,ny,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXY,info)
    call MPI_ALLREDUCE(MPI_IN_PLACE,uvrms,ny,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXY,info)
    call MPI_ALLREDUCE(MPI_IN_PLACE,uwrms,ny,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXY,info)
    call MPI_ALLREDUCE(MPI_IN_PLACE,vwrms,ny,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXY,info)
    urms=sqrt(urms/ngz)
    vrms=sqrt(vrms/ngz)
    wrms=sqrt(wrms/ngz)
    uvrms=uvrms/ngz
    uwrms=uwrms/ngz
    vwrms=vwrms/ngz

    if (iproc==0) print *,"End of calculation of spanwise Reynold stresses"

    ! MPI reconstruction and writing of profiles
    ! ==========================================
    call write_RFM_rms(urms/u_ref,0)
    call write_RFM_rms(vrms/u_ref,1)
    call write_RFM_rms(wrms/u_ref,1)
    call write_RFM_rms(uvrms/u_ref**2,1)
    call write_RFM_rms(uwrms/u_ref**2,1)
    call write_RFM_rms(vwrms/u_ref**2,1)

    ! if (iproc==0) print *,"Writting of mpi files for spanwise reynold stresses disabled"
    if (iproc==0) print *,"Writting of mpi files for spanwise reynold stresses activated"
        
  end subroutine compute_RFM_check

  !===============================================================================
  module subroutine write_RFM_rms(varms,ind)
  !===============================================================================
    !> check RFM: routine to write rms profile with MPI
  !===============================================================================
    use mod_mpi_part
    implicit none
    ! ----------------------------------------------------------------------------
    integer, intent(in) :: ind ! if 0 init communicator array iy
    real(wp), dimension(ny) :: varms
    ! ----------------------------------------------------------------------------
    ! compute rms profiles
    integer :: ip,j
    logical, dimension(0:nproc-1) :: iy
    real(wp), dimension(ngy) :: varmsg
    ! ----------------------------------------------------------------------------
    
    ! procs along y line
    ! ------------------
    if (ind==0) then
       iy=.false.
       if ((coord(1)==0).and.(coord(3)==0)) iy(iproc)=.true.
       ! share indicator iy
       do ip=0,nproc-1
          call MPI_BCAST(iy(ip),1,MPI_INTEGER,ip,COMM_global,info)
       enddo
    endif

    ! reconstruct and write varms
    ! ---------------------------
    
    if (iproc.ne.0) then
       if (iy(iproc)) call MPI_SEND(varms,ny,MPI_DOUBLE_PRECISION,0,tag,COMM_global,info)
    else
       varmsg(1:ny)=varms
       ! receive from other procs in the line
       do ip=1,nproc-1
          if (iy(ip)) call MPI_RECV(varmsg(coordy(ip)),ny,MPI_DOUBLE_PRECISION, &
                                    ip,tag,COMM_global,status,info)
       enddo
       ! write in file
       write(49) (varmsg(j),j=1,ngy)
    endif

  end subroutine write_RFM_rms

end submodule smod_RFM_init_vel
