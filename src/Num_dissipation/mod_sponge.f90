!=================================================================================
module mod_sponge
!=================================================================================
  !> Module for sponge zone
  ! New version XG 27/08/21
  ! restricted to 2D (periodic) configuration  
!=================================================================================
  use mod_constant
  use mod_mpi
  use mod_flow
  implicit none
  !-------------------------------------------------------------------------------
  ! Global damping coefficient ~> specified in param.ini
  real(wp) :: alpha_sz
  ! indicator to apply sponge zone for each proc
  logical, dimension(:), allocatable :: isponge
  ! local damping coefficient (2D version)
  real(wp), dimension(:,:), allocatable :: alpha1,alpha2
  ! local damping coefficient (3D version)
  real(wp), dimension(:,:,:), allocatable :: alpha3
  !-------------------------------------------------------------------------------

contains

  !===============================================================================
  subroutine init_sponge
  !===============================================================================
    !> Sponge zone initialization (only for 2D or z-periodic)
    !===============================================================================
    use mod_block
    use mod_time
    use warnstop
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: nbl,i,j,ierr
    integer :: is1,is2,d1_is,d2_is,js1,js2,d1_js,d2_js
    real(wp) :: beta
    real(wp), dimension(ngx) :: alphaxg
    real(wp), dimension(ngy) :: alphayg
    real(wp), dimension(ngx,ngy) :: alphaxgc
    real(wp), dimension(ngx,ngy) :: alphaygc
    ! ----------------------------------------------------------------------------

    if (iproc==0) write(*,*) repeat('=',70)
    if (iproc==0) print *,'Init sponge zone'
    ierr=0
    
    ! Allocate isponge for all procs
    ! ==============================
    allocate(isponge(0:nproc-1))
    isponge=.false.
    
    ! Number of current block
    ! =======================
    nbl=nob(iproc)
    
    ! check compatibility of sponge zone indices
    ! ==========================================
    if (bl(nbl)%is_sponge) then
       ! i-direction
       is1=bl(nbl)%is1
       is2=bl(nbl)%is2
       d1_is=bl(nbl)%d_is
       if (is1==1) d1_is=0
       d2_is=bl(nbl)%d_is
       if (is2==ngx) d2_is=0
       ! j-direction
       js1=bl(nbl)%js1
       js2=bl(nbl)%js2
       d1_js=bl(nbl)%d_js
       if (js1==1) d1_js=0
       d2_js=bl(nbl)%d_js
       if (js2==ngy) d2_js=0

       ! check i-direction
       if (is2>ngx) ierr=1
       if ((is2<ngx).and.(.not.is_curv)) ierr=2
       if ((is1+d1_is)>(is2-d2_is)) ierr=3
       ! check j-direction
       if (js2>ngy) ierr=4
       if ((js1+d1_js)>(js2-d2_js)) ierr=5
    endif

    ! print error message
    ! ===================
    select case (ierr)
    case (1)
       call mpistop('bad definition of sponge zone: is2>ngx',0)
    case (2)
       call mpistop('bad definition of sponge zone: is2 should be ngx (exit at right boundary)',0)
    case (3)
       call mpistop('bad definition of sponge zone: d_is too large',0)
    case (4)
       call mpistop('bad definition of sponge zone: js2>ngy',0)
    case (5)
       call mpistop('bad definition of sponge zone: d_js too large',0)
    end select
        
    if (bl(nbl)%is_sponge) then

       ! indices defining the sponge zone
       ! ================================
       ! i-direction
       is1=bl(nbl)%is1
       is2=bl(nbl)%is2
       d1_is=bl(nbl)%d_is
       if (is1==1) d1_is=0
       d2_is=bl(nbl)%d_is
       if (is2==ngy) d2_is=0
       ! j-direction
       js1=bl(nbl)%js1
       js2=bl(nbl)%js2
       d1_js=bl(nbl)%d_js
       if (js1==1) d1_js=0
       d2_js=bl(nbl)%d_js
       if (js2==ngy) d2_js=0

       ! nominal amplitude of Laplacian filter
       ! =====================================
       alpha_sz=0.01_wp

       ! damping exponent
       ! ================
       beta=1.5_wp

       ! indicator to know if domain 'iproc' contains part of the sponge zone
       ! =====================================================================
       isponge(iproc)=.true.
       if (coord(1)*nx+1 >is2) isponge(iproc)=(isponge(iproc)).and.(.false.)
       if (coord(1)*nx+nx<is1) isponge(iproc)=(isponge(iproc)).and.(.false.)
       if (coord(2)*ny+1 >js2) isponge(iproc)=(isponge(iproc)).and.(.false.)
       if (coord(2)*ny+ny<js1) isponge(iproc)=(isponge(iproc)).and.(.false.)

       ! profile of the sponge zone along x (global grid)
       ! ==================================

       if (is_curv) then
          ! for curvilinear grids
          ! ---------------------
          alphaxgc=0.0_wp

          do i=is1,is1+d1_is-1
             alphaxgc(i,1:ngy)=(dble(i-is1)/dble(d1_is))**beta
          enddo

          do i=is1+d1_is,is2
             alphaxgc(i,1:ngy)=1.0_wp
          enddo
          
          do i=is2-d2_is+1,is2
             alphaxgc(i,1:ngy)= (dble(is2-i)/dble(d2_is))**beta
          enddo
       else
          ! for Cartesian grids
          ! -------------------
          alphaxg=0.0_wp

          ! Based on indices of grid
          do i=is1,is1+d1_is-1
             alphaxg(i)=(dble(i-is1)/dble(d1_is))**beta
          enddo

          do i=is1+d1_is,is2-d2_is
             alphaxg(i)=1.0_wp
          enddo

          do i=is2-d2_is+1,is2
             alphaxg(i)= (dble(is2-i)/dble(d2_is))**beta
          enddo

          ! ! Based on global grid
          ! do i=is1,is1+d1_is-1
          !    alphaxg(i)=((xg(i)-xg(is1))/(xg(is1+d1_is-1)-xg(is1)))**beta
          ! enddo

          ! do i=is1+d1_is,is2-d2_is
          !    alphaxg(i)=1.0_wp
          ! enddo

          ! do i=is2-d2_is+1,is2
          !    alphaxg(i)=((xg(i)-xg(is2))/(xg(is2-d2_is+1)-xg(is2)))**beta
          ! enddo
       endif
    
       ! profile of the sponge zone along y (global grid)
       ! ==================================
       
       if (is_curv) then
          ! for curvilinear grids
          ! ---------------------
          alphaygc=0.0_wp

          do j=js1,js1+d1_js-1
             alphaygc(1:ngx,j)= (dble(j-js1)/dble(d1_js))**beta
          enddo

          do j=js1+d1_js,js2-d2_js
             alphaygc(1:ngx,j)= 1.0_wp
          enddo
          
          do j=js2-d2_js+1,js2
             alphaygc(1:ngx,j)= (dble(js2-j)/dble(d2_js))**beta
          enddo
       else
          ! for Cartesian grids
          ! -------------------
          alphayg=0.0_wp

          ! Based on indices of grid
          do j=js1,js1+d1_js-1
             alphayg(j)= (dble(j-js1)/dble(d1_js))**beta
          enddo

          do j=js1+d1_js,js2-d2_js
             alphayg(j)= 1.0_wp
          enddo

          do j=js2-d2_js+1,js2
             alphayg(j)= (dble(js2-j)/dble(d2_js))**beta
          enddo

          ! ! Based on global grid
          ! do j=js1,js1+d1_js-1
          !    alphayg(j)= ((yg(j)-yg(js1))/(yg(js1+d1_js-1)-yg(js1)))**beta
          ! enddo

          ! do j=js1+d1_js,js2-d2_js
          !    alphayg(j)= 1.0_wp
          ! enddo

          ! do j=js2-d2_js+1,js2
          !    alphayg(j)= ((yg(j)-yg(js2))/(yg(js2-d2_js+1)-yg(js2)))**beta
          ! enddo
       endif

       ! local damping coefficient alpha1 (local grid)
       ! ================================
       if (isponge(iproc)) then
          allocate(alpha1(nx,ny))
          allocate(alpha2(nx,ny))
          if (is_curv) then
             do j=1,ny
                do i=1,nx
                   alpha1(i,j) = alphaxgc(i+coord(1)*nx,j+coord(2)*ny)*alphaygc(i+coord(1)*nx,j+coord(2)*ny)
                   alpha2(i,j) = alphaygc(i+coord(1)*nx,j+coord(2)*ny)
                enddo
             enddo
          else
             do j=1,ny
                do i=1,nx
                   alpha1(i,j) = alphaxg(i+coord(1)*nx)*alphayg(j+coord(2)*ny)
                   alpha2(i,j) = alphayg(j+coord(2)*ny)
                enddo
             enddo
          endif
          if (is_dissip_in_increments) then
             ! /!\ Change sign because added to increment
             alpha1=-alpha_sz*alpha1/4.0_wp/deltat
          else
             alpha1=alpha_sz*alpha1/4.0_wp
          endif
       endif
             
       ! user-defined output for sponge zone coefficient
       ! ===============================================
       !if (isponge(iproc)) then
       !   uvar(:,:,1,1)=-alpha1*deltat*4.0_wp/alpha_sz
       !   !uvar(:,:,1,1)=alpha1
       !   uvar(:,:,1,2)=alpha2
       !endif
       
    endif
    
    ! indicate proc with sponge zone
    ! ==============================
    !if (isponge(iproc)) print *,'isponge for proc',iproc
    
    if (iproc==0) print *,'End sponge zone'
    if (iproc==0) write(*,*) repeat('=',70)
    
  end subroutine init_sponge

  !===============================================================================
  subroutine apply_sponge
  !===============================================================================
    !> Apply sponge zone
    ! only periodic BC along z
  !===============================================================================
    use mod_eos
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2) :: rhof,rhouf,rhovf,rhowf,rhoef
    ! ----------------------------------------------------------------------------

    ! Laplacian filter in x
    ! =====================
    ! Compute fluctuations
    rhof  = rho_n  !- rho0
    rhouf = rhou_n !- rho0*u0 
    rhovf = rhov_n !- rho0*v0 
    rhowf = rhow_n !- rho0*w0 
    rhoef = rhoe_n !- (p0*igm1+0.5_wp*rho0*(u0**2+v0**2+w0**2))

    do k=1,nz
       do j=1,ny
          do i=ndx3,nfx3
             Krho(i,j,k) = Krho(i,j,k)+alpha1(i,j)*(-2.0_wp* rhof(i,j,k)+ rhof(i+1,j,k)+ rhof(i-1,j,k))
             Krhou(i,j,k)=Krhou(i,j,k)+alpha1(i,j)*(-2.0_wp*rhouf(i,j,k)+rhouf(i+1,j,k)+rhouf(i-1,j,k))
             Krhov(i,j,k)=Krhov(i,j,k)+alpha1(i,j)*(-2.0_wp*rhovf(i,j,k)+rhovf(i+1,j,k)+rhovf(i-1,j,k))
             Krhow(i,j,k)=Krhow(i,j,k)+alpha1(i,j)*(-2.0_wp*rhowf(i,j,k)+rhowf(i+1,j,k)+rhowf(i-1,j,k))
             Krhoe(i,j,k)=Krhoe(i,j,k)+alpha1(i,j)*(-2.0_wp*rhoef(i,j,k)+rhoef(i+1,j,k)+rhoef(i-1,j,k))
          enddo
       enddo
    enddo

    ! Laplacian filter in y
    ! =====================
    ! Compute fluctuations
    rhof  = rho_n  !- rho0
    rhouf = rhou_n !- rho0*u0 
    rhovf = rhov_n !- rho0*v0 
    rhowf = rhow_n !- rho0*w0 
    rhoef = rhoe_n !- (p0*igm1+0.5_wp*rho0*(u0**2+v0**2+w0**2))

    do k=1,nz
       do j=ndy3,nfy3
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)+alpha1(i,j)*(-2.0_wp* rhof(i,j,k)+ rhof(i,j+1,k)+ rhof(i,j-1,k))
             Krhou(i,j,k)=Krhou(i,j,k)+alpha1(i,j)*(-2.0_wp*rhouf(i,j,k)+rhouf(i,j+1,k)+rhouf(i,j-1,k))
             Krhov(i,j,k)=Krhov(i,j,k)+alpha1(i,j)*(-2.0_wp*rhovf(i,j,k)+rhovf(i,j+1,k)+rhovf(i,j-1,k))
             Krhow(i,j,k)=Krhow(i,j,k)+alpha1(i,j)*(-2.0_wp*rhowf(i,j,k)+rhowf(i,j+1,k)+rhowf(i,j-1,k))
             Krhoe(i,j,k)=Krhoe(i,j,k)+alpha1(i,j)*(-2.0_wp*rhoef(i,j,k)+rhoef(i,j+1,k)+rhoef(i,j-1,k))
          enddo
       enddo
    enddo
    
    if (is_2d) return
   
    ! Laplacian filter in z
    ! =====================
    ! Compute fluctuations
    rhof  = rho_n  !- rho0
    rhouf = rhou_n !- rho0*u0 
    rhovf = rhov_n !- rho0*v0 
    rhowf = rhow_n !- rho0*w0 
    rhoef = rhoe_n !- (p0*igm1+0.5_wp*rho0*(u0**2+v0**2+w0**2))

    do k=ndz3,nfz3
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)+alpha1(i,j)*(-2.0_wp* rhof(i,j,k)+ rhof(i,j,k+1)+ rhof(i,j,k-1))
             Krhou(i,j,k)=Krhou(i,j,k)+alpha1(i,j)*(-2.0_wp*rhouf(i,j,k)+rhouf(i,j,k+1)+rhouf(i,j,k-1))
             Krhov(i,j,k)=Krhov(i,j,k)+alpha1(i,j)*(-2.0_wp*rhovf(i,j,k)+rhovf(i,j,k+1)+rhovf(i,j,k-1))
             Krhow(i,j,k)=Krhow(i,j,k)+alpha1(i,j)*(-2.0_wp*rhowf(i,j,k)+rhowf(i,j,k+1)+rhowf(i,j,k-1))
             Krhoe(i,j,k)=Krhoe(i,j,k)+alpha1(i,j)*(-2.0_wp*rhoef(i,j,k)+rhoef(i,j,k+1)+rhoef(i,j,k-1))
          enddo
       enddo
    enddo

  end subroutine apply_sponge

  !===============================================================================
  subroutine init_sponge_3d
  !===============================================================================
    !> Sponge zone initialization (only for 2D or z-periodic)
    !===============================================================================
    use mod_block
    use mod_time
    use warnstop
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: nbl,i,j,k,ierr
    integer :: is1,is2,d1_is,d2_is
    integer :: js1,js2,d1_js,d2_js
    integer :: ks1,ks2,d1_ks,d2_ks
    real(wp) :: beta
    real(wp), dimension(:), allocatable :: alphaxg,alphayg,alphazg
    real(wp), dimension(:,:), allocatable :: alphaxgc,alphaygc
    real(wp), dimension(:,:,:), allocatable :: alphaxgc3,alphaygc3,alphazgc3
    ! ----------------------------------------------------------------------------

    if (iproc==0) write(*,*) repeat('=',70)
    if (iproc==0) print *,'Init sponge zone (3D version)'
    ierr=0

    ! Allocate isponge for all procs
    ! ==============================
    allocate(isponge(0:nproc-1))
    isponge=.false.

    ! Number of current block
    ! =======================
    nbl=nob(iproc)

    ! check compatibility of sponge zone indices
    ! ==========================================
    if (bl(nbl)%is_sponge) then
       ! i-direction
       is1=bl(nbl)%is1
       is2=bl(nbl)%is2
       d1_is=bl(nbl)%d_is
       if (is1==1) d1_is=0
       d2_is=bl(nbl)%d_is
       if (is2==ngx) d2_is=0
       ! j-direction
       js1=bl(nbl)%js1
       js2=bl(nbl)%js2
       d1_js=bl(nbl)%d_js
       if (js1==1) d1_js=0
       d2_js=bl(nbl)%d_js
       if (js2==ngy) d2_js=0
       ! k-direction
       ks1=bl(nbl)%ks1
       ks2=bl(nbl)%ks2
       d1_ks=bl(nbl)%d_ks
       if (ks1==1) d1_ks=0
       d2_ks=bl(nbl)%d_ks
       if (ks2==ngz) d2_ks=0

       ! check i-direction
       if (is2>ngx) ierr=1
       if ((is2<ngx).and.(.not.is_curv).and.(.not.is_curv3)) ierr=2
       if ((is1+d1_is)>(is2-d2_is)) ierr=3
       ! check j-direction
       if (js2>ngy) ierr=4
       if ((js1+d1_js)>(js2-d2_js)) ierr=5
       ! check k-direction
       if (ks2>ngz) ierr=6
       if ((ks1+d1_ks)>(ks2-d2_ks)) ierr=7
    endif

    ! print error message
    ! ===================
    select case (ierr)
    case (1)
       call mpistop('bad definition of sponge zone: is2>ngx',0)
    case (2)
       call mpistop('bad definition of sponge zone: is2 should be ngx (exit at right boundary)',0)
    case (3)
       call mpistop('bad definition of sponge zone: d_is too large',0)
    case (4)
       call mpistop('bad definition of sponge zone: js2>ngy',0)
    case (5)
       call mpistop('bad definition of sponge zone: d_js too large',0)
    case (6)
       call mpistop('bad definition of sponge zone: ks2>ngz',0)
    case (7)
       call mpistop('bad definition of sponge zone: d_ks too large',0)
    end select

    if (bl(nbl)%is_sponge) then

       ! indices defining the sponge zone
       ! ================================
       ! i-direction
       is1=bl(nbl)%is1
       is2=bl(nbl)%is2
       d1_is=bl(nbl)%d_is
       if (is1==1) d1_is=0
       d2_is=bl(nbl)%d_is
       if (is2==ngy) d2_is=0
       ! j-direction
       js1=bl(nbl)%js1
       js2=bl(nbl)%js2
       d1_js=bl(nbl)%d_js
       if (js1==1) d1_js=0
       d2_js=bl(nbl)%d_js
       if (js2==ngy) d2_js=0
       ! k-direction
       ks1=bl(nbl)%ks1
       ks2=bl(nbl)%ks2
       d1_ks=bl(nbl)%d_ks
       if (ks1==1) d1_ks=0
       d2_ks=bl(nbl)%d_ks
       if (ks2==ngz) d2_ks=0

       ! nominal amplitude of Laplacian filter
       ! =====================================
       !alpha_sz=0.005_wp
       alpha_sz=0.01_wp

       ! damping exponent
       ! ================
       beta=1.5_wp

       ! indicator to know if domain 'iproc' contains part of the sponge zone
       ! =====================================================================
       isponge(iproc)=.true.
       if (coord(1)*nx+1 >is2) isponge(iproc)=(isponge(iproc)).and.(.false.)
       if (coord(1)*nx+nx<is1) isponge(iproc)=(isponge(iproc)).and.(.false.)
       if (coord(2)*ny+1 >js2) isponge(iproc)=(isponge(iproc)).and.(.false.)
       if (coord(2)*ny+ny<js1) isponge(iproc)=(isponge(iproc)).and.(.false.)
       if (coord(3)*nz+1 >ks2) isponge(iproc)=(isponge(iproc)).and.(.false.)
       if (coord(3)*nz+nz<ks1) isponge(iproc)=(isponge(iproc)).and.(.false.)

       ! profile of the sponge zone along x (global grid)
       ! ==================================

       if (is_curv3) then
          ! for 3D curvilinear grids
          ! ------------------------
          allocate(alphaxgc3(ngx,ngy,ngz))
          alphaxgc3=0.0_wp

          do i=is1,is1+d1_is-1
             alphaxgc3(i,1:ngy,1:ngz)=(dble(i-is1)/dble(d1_is))**beta
          enddo

          do i=is1+d1_is,is2
             alphaxgc3(i,1:ngy,1:ngz)=1.0_wp
          enddo

          do i=is2-d2_is+1,is2
             alphaxgc3(i,1:ngy,1:ngz)= (dble(is2-i)/dble(d2_is))**beta
          enddo
       else
          if (is_curv) then
             ! for 2D curvilinear grids
             ! ------------------------
             allocate(alphaxgc(ngx,ngy))
             alphaxgc=0.0_wp

             do i=is1,is1+d1_is-1
                alphaxgc(i,1:ngy)=(dble(i-is1)/dble(d1_is))**beta
             enddo

             do i=is1+d1_is,is2
                alphaxgc(i,1:ngy)=1.0_wp
             enddo

             do i=is2-d2_is+1,is2
                alphaxgc(i,1:ngy)= (dble(is2-i)/dble(d2_is))**beta
             enddo
          else
             ! for Cartesian grids
             ! -------------------
             allocate(alphaxg(ngx))
             alphaxg=0.0_wp

             do i=is1,is1+d1_is-1
                alphaxg(i)=((xg(i)-xg(is1))/(xg(is1+d1_is-1)-xg(is1)))**beta
             enddo

             do i=is1+d1_is,is2-d2_is
                alphaxg(i)=1.0_wp
             enddo

             do i=is2-d2_is+1,is2
                alphaxg(i)=((xg(i)-xg(is2))/(xg(is2-d2_is+1)-xg(is2)))**beta
             enddo
          endif
       endif

       ! profile of the sponge zone along y (global grid)
       ! ==================================

       if (is_curv3) then
          ! for 2D curvilinear grids
          ! ------------------------
          allocate(alphaygc3(ngx,ngy,ngz))
          alphaygc3=0.0_wp

          do j=js1,js1+d1_js-1
             alphaygc3(1:ngx,j,1:ngz)= (dble(j-js1)/dble(d1_js))**beta
          enddo

          do j=js1+d1_js,js2-d2_js
             alphaygc3(1:ngx,j,1:ngz)= 1.0_wp
          enddo

          do j=js2-d2_js+1,js2
             alphaygc3(1:ngx,j,1:ngz)= (dble(js2-j)/dble(d2_js))**beta
          enddo
       else
          if (is_curv) then
             ! for 2D curvilinear grids
             ! ------------------------
             allocate(alphaygc(ngx,ngy))
             alphaygc=0.0_wp

             do j=js1,js1+d1_js-1
                alphaygc(1:ngx,j)= (dble(j-js1)/dble(d1_js))**beta
             enddo

             do j=js1+d1_js,js2-d2_js
                alphaygc(1:ngx,j)= 1.0_wp
             enddo

             do j=js2-d2_js+1,js2
                alphaygc(1:ngx,j)= (dble(js2-j)/dble(d2_js))**beta
             enddo
          else
             ! for Cartesian grids
             ! -------------------
             allocate(alphayg(ngy))
             alphayg=0.0_wp

             do j=js1,js1+d1_js-1
                alphayg(j)= ((yg(j)-yg(js1))/(yg(js1+d1_js-1)-yg(js1)))**beta
             enddo

             do j=js1+d1_js,js2-d2_js
                alphayg(j)= 1.0_wp
             enddo

             do j=js2-d2_js+1,js2
                alphayg(j)= ((yg(j)-yg(js2))/(yg(js2-d2_js+1)-yg(js2)))**beta
             enddo
          endif
       endif

       ! profile of the sponge zone along z (global grid)
       ! ==================================

       if (is_curv3) then
          ! for 2D curvilinear grids
          ! ------------------------
          allocate(alphazgc3(ngx,ngy,ngz))
          alphazgc3=0.0_wp

          do k=ks1,ks1+d1_ks-1
             alphazgc3(1:ngx,1:ngy,k)= (dble(k-ks1)/dble(d1_ks))**beta
          enddo

          do k=ks1+d1_ks,ks2-d2_ks
             alphazgc3(1:ngx,1:ngy,k)= 1.0_wp
          enddo

          do k=ks2-d2_ks+1,ks2
             alphazgc3(1:ngx,1:ngy,k)= (dble(ks2-k)/dble(d2_ks))**beta
          enddo
       else
          ! for Cartesian grids
          ! -------------------
          allocate(alphazg(ngz))
          alphazg=0.0_wp

          do k=ks1,ks1+d1_ks-1
             alphazg(k)= ((zg(k)-zg(ks1))/(zg(ks1+d1_ks-1)-zg(ks1)))**beta
          enddo

          do k=ks1+d1_ks,ks2-d2_ks
             alphazg(k)= 1.0_wp
          enddo

          do k=ks2-d2_ks+1,ks2
             alphazg(k)= ((zg(k)-zg(ks2))/(zg(ks2-d2_ks+1)-zg(ks2)))**beta
          enddo
       endif

       ! local damping coefficient alpha1 (local grid)
       ! ================================
       if (isponge(iproc)) then
          allocate(alpha3(nx,ny,nz))

          if (is_curv3) then
             do k=1,nz
                do j=1,ny
                   do i=1,nx
                      alpha3(i,j,k) = alphaxgc3(i+coord(1)*nx,j+coord(2)*ny,k+coord(3)*nz) &
                                     *alphaygc3(i+coord(1)*nx,j+coord(2)*ny,k+coord(3)*nz) &
                                     *alphazgc3(i+coord(1)*nx,j+coord(2)*ny,k+coord(3)*nz)
                   enddo
                enddo
             enddo
             deallocate(alphaxgc3,alphaygc3,alphazgc3)
          else
             if (is_curv) then
                do k=1,nz
                   do j=1,ny
                      do i=1,nx
                         alpha3(i,j,k) = alphaxgc(i+coord(1)*nx,j+coord(2)*ny) &
                                        *alphaygc(i+coord(1)*nx,j+coord(2)*ny)*alphazg(k+coord(3)*nz)
                      enddo
                   enddo
                enddo
                deallocate(alphaxgc,alphaygc,alphazg)
             else
                do k=1,nz
                   do j=1,ny
                      do i=1,nx
                         alpha3(i,j,k) = alphaxg(i+coord(1)*nx)*alphayg(j+coord(2)*ny)*alphazg(k+coord(3)*nz)
                      enddo
                   enddo
                enddo
                deallocate(alphaxg,alphayg,alphazg)
             endif
          endif

          if (is_dissip_in_increments) then
             ! /!\ Change sign because added to increment
             alpha3=-alpha_sz*alpha3/4.0_wp/deltat
          else
             alpha3=alpha_sz*alpha3/4.0_wp
          endif
       endif

       ! user-defined output for sponge zone coefficient
       ! ===============================================
       !if (isponge(iproc)) then
       !   uvar(:,:,1,1)=-alpha3*deltat*4.0_wp/alpha_sz
       !   !uvar(:,:,1,1)=alpha3
       !endif

    endif

    ! indicate proc with sponge zone
    ! ==============================
    !if (isponge(iproc)) print *,'isponge for proc',iproc

    if (iproc==0) print *,'End sponge zone'
    if (iproc==0) write(*,*) repeat('=',70)

  end subroutine init_sponge_3d

  !===============================================================================
  subroutine apply_sponge_3d
  !===============================================================================
    !> Apply sponge zone
    ! only periodic BC along z
  !===============================================================================
    use mod_eos
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2) :: rhof,rhouf,rhovf,rhowf,rhoef
    ! ----------------------------------------------------------------------------

    ! Laplacian filter in x
    ! =====================
    ! Compute fluctuations
    rhof  = rho_n  !- rho0
    rhouf = rhou_n !- rho0*u0
    rhovf = rhov_n !- rho0*v0
    rhowf = rhow_n !- rho0*w0
    rhoef = rhoe_n !- (p0*igm1+0.5_wp*rho0*(u0**2+v0**2+w0**2))

    do k=1,nz
       do j=1,ny
          do i=ndx3,nfx3
             Krho(i,j,k) = Krho(i,j,k)+alpha3(i,j,k)*(-2.0_wp* rhof(i,j,k)+ rhof(i+1,j,k)+ rhof(i-1,j,k))
             Krhou(i,j,k)=Krhou(i,j,k)+alpha3(i,j,k)*(-2.0_wp*rhouf(i,j,k)+rhouf(i+1,j,k)+rhouf(i-1,j,k))
             Krhov(i,j,k)=Krhov(i,j,k)+alpha3(i,j,k)*(-2.0_wp*rhovf(i,j,k)+rhovf(i+1,j,k)+rhovf(i-1,j,k))
             Krhow(i,j,k)=Krhow(i,j,k)+alpha3(i,j,k)*(-2.0_wp*rhowf(i,j,k)+rhowf(i+1,j,k)+rhowf(i-1,j,k))
             Krhoe(i,j,k)=Krhoe(i,j,k)+alpha3(i,j,k)*(-2.0_wp*rhoef(i,j,k)+rhoef(i+1,j,k)+rhoef(i-1,j,k))
          enddo
       enddo
    enddo

    ! Laplacian filter in y
    ! =====================
    ! Compute fluctuations
    rhof  = rho_n  !- rho0
    rhouf = rhou_n !- rho0*u0
    rhovf = rhov_n !- rho0*v0
    rhowf = rhow_n !- rho0*w0
    rhoef = rhoe_n !- (p0*igm1+0.5_wp*rho0*(u0**2+v0**2+w0**2))

    do k=1,nz
       do j=ndy3,nfy3
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)+alpha3(i,j,k)*(-2.0_wp* rhof(i,j,k)+ rhof(i,j+1,k)+ rhof(i,j-1,k))
             Krhou(i,j,k)=Krhou(i,j,k)+alpha3(i,j,k)*(-2.0_wp*rhouf(i,j,k)+rhouf(i,j+1,k)+rhouf(i,j-1,k))
             Krhov(i,j,k)=Krhov(i,j,k)+alpha3(i,j,k)*(-2.0_wp*rhovf(i,j,k)+rhovf(i,j+1,k)+rhovf(i,j-1,k))
             Krhow(i,j,k)=Krhow(i,j,k)+alpha3(i,j,k)*(-2.0_wp*rhowf(i,j,k)+rhowf(i,j+1,k)+rhowf(i,j-1,k))
             Krhoe(i,j,k)=Krhoe(i,j,k)+alpha3(i,j,k)*(-2.0_wp*rhoef(i,j,k)+rhoef(i,j+1,k)+rhoef(i,j-1,k))
          enddo
       enddo
    enddo

    if (is_2d) return

    ! Laplacian filter in z
    ! =====================
    ! Compute fluctuations
    rhof  = rho_n  !- rho0
    rhouf = rhou_n !- rho0*u0
    rhovf = rhov_n !- rho0*v0
    rhowf = rhow_n !- rho0*w0
    rhoef = rhoe_n !- (p0*igm1+0.5_wp*rho0*(u0**2+v0**2+w0**2))

    do k=ndz3,nfz3
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)+alpha3(i,j,k)*(-2.0_wp* rhof(i,j,k)+ rhof(i,j,k+1)+ rhof(i,j,k-1))
             Krhou(i,j,k)=Krhou(i,j,k)+alpha3(i,j,k)*(-2.0_wp*rhouf(i,j,k)+rhouf(i,j,k+1)+rhouf(i,j,k-1))
             Krhov(i,j,k)=Krhov(i,j,k)+alpha3(i,j,k)*(-2.0_wp*rhovf(i,j,k)+rhovf(i,j,k+1)+rhovf(i,j,k-1))
             Krhow(i,j,k)=Krhow(i,j,k)+alpha3(i,j,k)*(-2.0_wp*rhowf(i,j,k)+rhowf(i,j,k+1)+rhowf(i,j,k-1))
             Krhoe(i,j,k)=Krhoe(i,j,k)+alpha3(i,j,k)*(-2.0_wp*rhoef(i,j,k)+rhoef(i,j,k+1)+rhoef(i,j,k-1))
          enddo
       enddo
    enddo

  end subroutine apply_sponge_3d

end module mod_sponge
