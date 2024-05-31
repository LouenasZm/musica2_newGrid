!=================================================================================
submodule (mod_grid_metrics_c3) smod_grid_metrics3_c3
!=================================================================================
  !> Module to compute full 3D curvilinear metrics
  !> with Geometric Conservation Law (GCL)
!=================================================================================
  
contains

  !===============================================================================
  module subroutine grid_metrics_ijacob_3d
  !===============================================================================
    !> Compute Jacobians of metrics (full 3D version)
    !> - stencil -ngh:+ngh for inviscid fluxes -
  !===============================================================================
    use mod_mpi ! for iproc in print screen
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    ! ---------------------------------------------------------------------------

    ! Allocate Jacobian array
    ! -----------------------
    allocate(ijacob3(0:nfx1,0:nfy1,0:nfz1))

    ! Jacobian of metrics transformation
    ! ==================================
    ijacob3=1.0_wp

    do k=1,nz
       do j=1,ny
          do i=1,nx
             ijacob3(i,j,k)= xdksi(i,j,k)*(ydeta(i,j,k)*zdphi(i,j,k)-ydphi(i,j,k)*zdeta(i,j,k)) &
                           - xdeta(i,j,k)*(ydksi(i,j,k)*zdphi(i,j,k)-ydphi(i,j,k)*zdksi(i,j,k)) &
                           + xdphi(i,j,k)*(ydksi(i,j,k)*zdeta(i,j,k)-ydeta(i,j,k)*zdksi(i,j,k))
             if (ijacob3(i,j,k)==0) print *,'jac',iproc,i,j,k
          enddo
       enddo
    enddo

    ! Define Jacobian at ghost cells 0 & nx+1 for conservative formulation
    ! =======================================
    ! used for example in mod_art_visc

    ! at imin: i=0
    if (ndx1==1) then
       ijacob3(0,:,:)=ijacob3(1,:,:)
    else
       i=ndx1
       do k=1,nz
          do j=1,ny
             ijacob3(i,j,k)= xdksi(i,j,k)*(ydeta(i,j,k)*zdphi(i,j,k)-ydphi(i,j,k)*zdeta(i,j,k)) &
                           - xdeta(i,j,k)*(ydksi(i,j,k)*zdphi(i,j,k)-ydphi(i,j,k)*zdksi(i,j,k)) &
                           + xdphi(i,j,k)*(ydksi(i,j,k)*zdeta(i,j,k)-ydeta(i,j,k)*zdksi(i,j,k))
             if (ijacob3(i,j,k)==0) ijacob3(i,j,k)=1.0e-20_wp
          enddo
       enddo
    endif

    ! at imax: i=nx+1
    if (nfx1.ne.nx) then
       i=nfx1
       do k=1,nz
          do j=1,ny
             ijacob3(i,j,k)= xdksi(i,j,k)*(ydeta(i,j,k)*zdphi(i,j,k)-ydphi(i,j,k)*zdeta(i,j,k)) &
                           - xdeta(i,j,k)*(ydksi(i,j,k)*zdphi(i,j,k)-ydphi(i,j,k)*zdksi(i,j,k)) &
                           + xdphi(i,j,k)*(ydksi(i,j,k)*zdeta(i,j,k)-ydeta(i,j,k)*zdksi(i,j,k))
             if (ijacob3(i,j,k)==0) ijacob3(i,j,k)=1.0e-20_wp
          enddo
       enddo
    endif

    ! at jmin: j=0
    if (ndy1==1) then
       ijacob3(:,0,:)=ijacob3(:,1,:)
    else
       j=ndy1
       do k=1,nz
          do i=1,nx
             ijacob3(i,j,k)= xdksi(i,j,k)*(ydeta(i,j,k)*zdphi(i,j,k)-ydphi(i,j,k)*zdeta(i,j,k)) &
                           - xdeta(i,j,k)*(ydksi(i,j,k)*zdphi(i,j,k)-ydphi(i,j,k)*zdksi(i,j,k)) &
                           + xdphi(i,j,k)*(ydksi(i,j,k)*zdeta(i,j,k)-ydeta(i,j,k)*zdksi(i,j,k))
             if (ijacob3(i,j,k)==0) ijacob3(i,j,k)=1.0e-20_wp
          enddo
       enddo
    endif

    ! at jmax: j=ny+1
    if (nfy1.ne.ny) then
       j=nfy1
       do k=1,nz
          do i=1,nx
             ijacob3(i,j,k)= xdksi(i,j,k)*(ydeta(i,j,k)*zdphi(i,j,k)-ydphi(i,j,k)*zdeta(i,j,k)) &
                           - xdeta(i,j,k)*(ydksi(i,j,k)*zdphi(i,j,k)-ydphi(i,j,k)*zdksi(i,j,k)) &
                           + xdphi(i,j,k)*(ydksi(i,j,k)*zdeta(i,j,k)-ydeta(i,j,k)*zdksi(i,j,k))
             if (ijacob3(i,j,k)==0) ijacob3(i,j,k)=1.0e-20_wp
          enddo
       enddo
    endif
    
    ! at kmin: k=0
    if (ndz1==1) then
       ijacob3(:,:,0)=ijacob3(:,:,1)
    else
       k=ndz1
       do j=1,ny
          do i=1,nx
             ijacob3(i,j,k)= xdksi(i,j,k)*(ydeta(i,j,k)*zdphi(i,j,k)-ydphi(i,j,k)*zdeta(i,j,k)) &
                           - xdeta(i,j,k)*(ydksi(i,j,k)*zdphi(i,j,k)-ydphi(i,j,k)*zdksi(i,j,k)) &
                           + xdphi(i,j,k)*(ydksi(i,j,k)*zdeta(i,j,k)-ydeta(i,j,k)*zdksi(i,j,k))
             if (ijacob3(i,j,k)==0) ijacob3(i,j,k)=1.0e-20_wp
          enddo
       enddo
    endif

    ! at jmax: k=nz+1
    if (nfz1.ne.nz) then
       k=nfz1
       do j=1,ny
          do i=1,nx
             ijacob3(i,j,k)= xdksi(i,j,k)*(ydeta(i,j,k)*zdphi(i,j,k)-ydphi(i,j,k)*zdeta(i,j,k)) &
                           - xdeta(i,j,k)*(ydksi(i,j,k)*zdphi(i,j,k)-ydphi(i,j,k)*zdksi(i,j,k)) &
                           + xdphi(i,j,k)*(ydksi(i,j,k)*zdeta(i,j,k)-ydeta(i,j,k)*zdksi(i,j,k))
             if (ijacob3(i,j,k)==0) ijacob3(i,j,k)=1.0e-20_wp
          enddo
       enddo
    endif

    ! Inverse of Jacob3ian
    ! ====================
    ijacob3=1.0_wp/ijacob3
    
  end subroutine grid_metrics_ijacob_3d

  !===============================================================================
  module subroutine grid_metrics_gradients_3d
  !===============================================================================
    !> Compute norms of metrics gradients (full 3D version)
    !> - stencil -ngh:+ngh for inviscid fluxes -
  !===============================================================================
    use mod_grid
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    ! ---------------------------------------------------------------------------

    ! Allocate arrays
    ! ---------------
    allocate(g3_ksi(0:nfx1,0:nfy1,0:nfz1))
    allocate(g3_eta(0:nfx1,0:nfy1,0:nfz1))
    allocate(g3_phi(0:nfx1,0:nfy1,0:nfz1))

    ! Norm of metrics gradients
    ! =========================
    ! sqrt[grad(ksi).grad(ksi)]=sqrt(ksi_x^2+ksi_y^2+ksi_z^2)/|J|
    ! sqrt[grad(eta).grad(eta)]=sqrt(eta_x^2+eta_y^2+eta_z^2)/|J|
    ! sqrt[grad(eta).grad(eta)]=sqrt(eta_x^2+eta_y^2+eta_z^2)/|J|   
    do k=1,nz
       do j=1,ny
          do i=1,nx
             g3_ksi(i,j,k)=sqrt(ksi_x(i,j,k)**2+ksi_y(i,j,k)**2+ksi_z(i,j,k)**2)*abs(ijacob3(i,j,k))
             g3_eta(i,j,k)=sqrt(eta_x(i,j,k)**2+eta_y(i,j,k)**2+eta_z(i,j,k)**2)*abs(ijacob3(i,j,k))
             g3_phi(i,j,k)=sqrt(phi_x(i,j,k)**2+phi_y(i,j,k)**2+phi_z(i,j,k)**2)*abs(ijacob3(i,j,k))
          enddo
       enddo
    enddo

    ! Define gradients at ghost cells 0 & nx+1 for conservative formulation
    ! ========================================
    ! used for example in mod_art_visc

    ! at imin: i=0
    if (ndx1==1) then
       g3_ksi(0,:,:)=g3_ksi(1,:,:)
       g3_eta(0,:,:)=g3_eta(1,:,:)
       g3_phi(0,:,:)=g3_phi(1,:,:)
    else
       i=ndx1
       do k=1,nz
          do j=1,ny
             g3_ksi(i,j,k)=sqrt(ksi_x(i,j,k)**2+ksi_y(i,j,k)**2+ksi_z(i,j,k)**2)*abs(ijacob3(i,j,k))
             g3_eta(i,j,k)=sqrt(eta_x(i,j,k)**2+eta_y(i,j,k)**2+eta_z(i,j,k)**2)*abs(ijacob3(i,j,k))
             g3_phi(i,j,k)=sqrt(phi_x(i,j,k)**2+phi_y(i,j,k)**2+phi_z(i,j,k)**2)*abs(ijacob3(i,j,k))
          enddo
       enddo
    endif

    ! at imax: i=nx+1
    if (nfx1.ne.nx) then
       i=nfx1
       do k=1,nz
          do j=1,ny
             g3_ksi(i,j,k)=sqrt(ksi_x(i,j,k)**2+ksi_y(i,j,k)**2+ksi_z(i,j,k)**2)*abs(ijacob3(i,j,k))
             g3_eta(i,j,k)=sqrt(eta_x(i,j,k)**2+eta_y(i,j,k)**2+eta_z(i,j,k)**2)*abs(ijacob3(i,j,k))
             g3_phi(i,j,k)=sqrt(phi_x(i,j,k)**2+phi_y(i,j,k)**2+phi_z(i,j,k)**2)*abs(ijacob3(i,j,k))
          enddo
       enddo
    endif

    ! at jmin: j=0
    if (ndy1==1) then
       g3_ksi(:,0,:)=g3_ksi(:,1,:)
       g3_eta(:,0,:)=g3_eta(:,1,:)
       g3_phi(:,0,:)=g3_phi(:,1,:)
    else
       j=ndy1
       do k=1,nz
          do i=1,nx
             g3_ksi(i,j,k)=sqrt(ksi_x(i,j,k)**2+ksi_y(i,j,k)**2+ksi_z(i,j,k)**2)*abs(ijacob3(i,j,k))
             g3_eta(i,j,k)=sqrt(eta_x(i,j,k)**2+eta_y(i,j,k)**2+eta_z(i,j,k)**2)*abs(ijacob3(i,j,k))
             g3_phi(i,j,k)=sqrt(phi_x(i,j,k)**2+phi_y(i,j,k)**2+phi_z(i,j,k)**2)*abs(ijacob3(i,j,k))
          enddo
       enddo
    endif

    ! at jmax: j=ny+1
    if (nfy1.ne.ny) then
       j=nfy1
       do k=1,nz
          do i=1,nx
             g3_ksi(i,j,k)=sqrt(ksi_x(i,j,k)**2+ksi_y(i,j,k)**2+ksi_z(i,j,k)**2)*abs(ijacob3(i,j,k))
             g3_eta(i,j,k)=sqrt(eta_x(i,j,k)**2+eta_y(i,j,k)**2+eta_z(i,j,k)**2)*abs(ijacob3(i,j,k))
             g3_phi(i,j,k)=sqrt(phi_x(i,j,k)**2+phi_y(i,j,k)**2+phi_z(i,j,k)**2)*abs(ijacob3(i,j,k))
          enddo
       enddo
    endif
    
    ! at kmin: k=0
    if (ndz1==1) then
       g3_ksi(:,:,0)=g3_ksi(:,:,1)
       g3_eta(:,:,0)=g3_eta(:,:,1)
       g3_phi(:,:,0)=g3_phi(:,:,1)
    else
       k=ndz1
       do j=1,ny
          do i=1,nx
             g3_ksi(i,j,k)=sqrt(ksi_x(i,j,k)**2+ksi_y(i,j,k)**2+ksi_z(i,j,k)**2)*abs(ijacob3(i,j,k))
             g3_eta(i,j,k)=sqrt(eta_x(i,j,k)**2+eta_y(i,j,k)**2+eta_z(i,j,k)**2)*abs(ijacob3(i,j,k))
             g3_phi(i,j,k)=sqrt(phi_x(i,j,k)**2+phi_y(i,j,k)**2+phi_z(i,j,k)**2)*abs(ijacob3(i,j,k))
          enddo
       enddo
    endif

    ! at jmax: k=nz+1
    if (nfz1.ne.nz) then
       k=nfz1
       do j=1,ny
          do i=1,nx
             g3_ksi(i,j,k)=sqrt(ksi_x(i,j,k)**2+ksi_y(i,j,k)**2+ksi_z(i,j,k)**2)*abs(ijacob3(i,j,k))
             g3_eta(i,j,k)=sqrt(eta_x(i,j,k)**2+eta_y(i,j,k)**2+eta_z(i,j,k)**2)*abs(ijacob3(i,j,k))
             g3_phi(i,j,k)=sqrt(phi_x(i,j,k)**2+phi_y(i,j,k)**2+phi_z(i,j,k)**2)*abs(ijacob3(i,j,k))
          enddo
       enddo
    endif

  end subroutine grid_metrics_gradients_3d

  !===============================================================================
  module subroutine grid_metrics_ijacob_3d_v
  !===============================================================================
    !> Compute Jacobians of metrics (full 3D version)
    !> - stencil -ngh_v:+ngh_v for viscous fluxes -
  !===============================================================================
    !use warnstop
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    ! ---------------------------------------------------------------------------

    ! Allocate Jacobian array
    ! -----------------------
    allocate(ijacob3_v(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))

    ! Jacobian of metrics transformation
    ! ==================================
    ijacob3_v=1.0_wp

    do k=1,nz
       do j=1,ny
          do i=1,nx
             ijacob3_v(i,j,k)= xdksi(i,j,k)*(ydeta(i,j,k)*zdphi(i,j,k)-ydphi(i,j,k)*zdeta(i,j,k)) &
                             - xdeta(i,j,k)*(ydksi(i,j,k)*zdphi(i,j,k)-ydphi(i,j,k)*zdksi(i,j,k)) &
                             + xdphi(i,j,k)*(ydksi(i,j,k)*zdeta(i,j,k)-ydeta(i,j,k)*zdksi(i,j,k))
             if (ijacob3_v(i,j,k)==0) then
                ijacob3_v(i,j,k)=1.0e-20_wp
             !   print *,'jac_v',iproc,i,j,k
             !   call mpistop('pb Jacobian',0)
             endif
          enddo
       enddo
    enddo

    ! Define Jacobian at ghost cells for double derivative ?? useful ??
    ! ==============================
    
    ! at imin: i=0
    if (ndx_v1.ne.1) then
       i=ndx_v1
       do k=1,nz
          do j=1,ny
             ijacob3_v(i,j,k)= xdksi(i,j,k)*(ydeta(i,j,k)*zdphi(i,j,k)-ydphi(i,j,k)*zdeta(i,j,k)) &
                             - xdeta(i,j,k)*(ydksi(i,j,k)*zdphi(i,j,k)-ydphi(i,j,k)*zdksi(i,j,k)) &
                             + xdphi(i,j,k)*(ydksi(i,j,k)*zdeta(i,j,k)-ydeta(i,j,k)*zdksi(i,j,k))
             if (ijacob3_v(i,j,k)==0) ijacob3_v(i,j,k)=1.0e-20_wp
          enddo
       enddo
    endif

    ! at imax: i=nx+1
    if (nfx_v1.ne.nx) then
       i=nfx_v1
       do k=1,nz
          do j=1,ny
             ijacob3_v(i,j,k)= xdksi(i,j,k)*(ydeta(i,j,k)*zdphi(i,j,k)-ydphi(i,j,k)*zdeta(i,j,k)) &
                             - xdeta(i,j,k)*(ydksi(i,j,k)*zdphi(i,j,k)-ydphi(i,j,k)*zdksi(i,j,k)) &
                             + xdphi(i,j,k)*(ydksi(i,j,k)*zdeta(i,j,k)-ydeta(i,j,k)*zdksi(i,j,k))
             if (ijacob3_v(i,j,k)==0) ijacob3_v(i,j,k)=1.0e-20_wp
          enddo
       enddo
    endif

    ! at jmin: j=0
    if (ndy_v1.ne.1) then
       j=ndy_v1
       do k=1,nz
          do i=1,nx
             ijacob3_v(i,j,k)= xdksi(i,j,k)*(ydeta(i,j,k)*zdphi(i,j,k)-ydphi(i,j,k)*zdeta(i,j,k)) &
                             - xdeta(i,j,k)*(ydksi(i,j,k)*zdphi(i,j,k)-ydphi(i,j,k)*zdksi(i,j,k)) &
                             + xdphi(i,j,k)*(ydksi(i,j,k)*zdeta(i,j,k)-ydeta(i,j,k)*zdksi(i,j,k))
             if (ijacob3_v(i,j,k)==0) ijacob3_v(i,j,k)=1.0e-20_wp
          enddo
       enddo
    endif

    ! at jmax: j=ny+1
    if (nfy_v1.ne.ny) then
       j=nfy_v1
       do k=1,nz
          do i=1,nx
             ijacob3_v(i,j,k)= xdksi(i,j,k)*(ydeta(i,j,k)*zdphi(i,j,k)-ydphi(i,j,k)*zdeta(i,j,k)) &
                             - xdeta(i,j,k)*(ydksi(i,j,k)*zdphi(i,j,k)-ydphi(i,j,k)*zdksi(i,j,k)) &
                             + xdphi(i,j,k)*(ydksi(i,j,k)*zdeta(i,j,k)-ydeta(i,j,k)*zdksi(i,j,k))
             if (ijacob3_v(i,j,k)==0) ijacob3_v(i,j,k)=1.0e-20_wp
          enddo
       enddo
    endif
    
    ! at kmin: k=0
    if (ndz_v1.ne.1) then
       k=ndz_v1
       do j=1,ny
          do i=1,nx
             ijacob3_v(i,j,k)= xdksi(i,j,k)*(ydeta(i,j,k)*zdphi(i,j,k)-ydphi(i,j,k)*zdeta(i,j,k)) &
                             - xdeta(i,j,k)*(ydksi(i,j,k)*zdphi(i,j,k)-ydphi(i,j,k)*zdksi(i,j,k)) &
                             + xdphi(i,j,k)*(ydksi(i,j,k)*zdeta(i,j,k)-ydeta(i,j,k)*zdksi(i,j,k))
             if (ijacob3_v(i,j,k)==0) ijacob3_v(i,j,k)=1.0e-20_wp
          enddo
       enddo
    endif

    ! at jmax: k=nz+1
    if (nfz_v1.ne.nz) then
       k=nfz_v1
       do j=1,ny
          do i=1,nx
             ijacob3_v(i,j,k)= xdksi(i,j,k)*(ydeta(i,j,k)*zdphi(i,j,k)-ydphi(i,j,k)*zdeta(i,j,k)) &
                             - xdeta(i,j,k)*(ydksi(i,j,k)*zdphi(i,j,k)-ydphi(i,j,k)*zdksi(i,j,k)) &
                             + xdphi(i,j,k)*(ydksi(i,j,k)*zdeta(i,j,k)-ydeta(i,j,k)*zdksi(i,j,k))
             if (ijacob3_v(i,j,k)==0) ijacob3_v(i,j,k)=1.0e-20_wp
          enddo
       enddo
    endif

    ! Inverse of Jacobian
    ! ===================
    ijacob3_v=1.0_wp/ijacob3_v

  end subroutine grid_metrics_ijacob_3d_v
  
  !===============================================================================
  module subroutine grid_metrics_ijacob_3d2_v
  !===============================================================================
    !> Compute Jacobians of metrics (full 3D version)
    !> - stencil -ngh_v:+ngh_v for viscous fluxes -
  !===============================================================================
    !use warnstop
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    ! ---------------------------------------------------------------------------

    ! Allocate Jacobian array
    ! -----------------------
    allocate(ijacob3_v(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))

    ! Jacobian of metrics transformation
    ! ==================================
    ijacob3_v=1.0_wp

    do k=1,nz
       do j=1,ny
          do i=1,nx
             ijacob3_v(i,j,k)= ksi_x_v(i,j,k)*(eta_y_v(i,j,k)*phi_z_v(i,j,k)-phi_y_v(i,j,k)*eta_z_v(i,j,k)) &
                             - eta_x_v(i,j,k)*(ksi_y_v(i,j,k)*phi_z_v(i,j,k)-phi_y_v(i,j,k)*ksi_z_v(i,j,k)) &
                             + phi_x_v(i,j,k)*(ksi_y_v(i,j,k)*eta_z_v(i,j,k)-eta_y_v(i,j,k)*ksi_z_v(i,j,k))
          enddo
       enddo
    enddo

    ! Define Jacobian at ghost cells for double derivative ?? useful ??
    ! ==============================
    
    ! at imin: i=0
    if (ndx_v1.ne.1) then
       i=ndx_v1
       do k=1,nz
          do j=1,ny
             ijacob3_v(i,j,k)= ksi_x_v(i,j,k)*(eta_y_v(i,j,k)*phi_z_v(i,j,k)-phi_y_v(i,j,k)*eta_z_v(i,j,k)) &
                             - eta_x_v(i,j,k)*(ksi_y_v(i,j,k)*phi_z_v(i,j,k)-phi_y_v(i,j,k)*ksi_z_v(i,j,k)) &
                             + phi_x_v(i,j,k)*(ksi_y_v(i,j,k)*eta_z_v(i,j,k)-eta_y_v(i,j,k)*ksi_z_v(i,j,k))
          enddo
       enddo
    endif

    ! at imax: i=nx+1
    if (nfx_v1.ne.nx) then
       i=nfx_v1
       do k=1,nz
          do j=1,ny
             ijacob3_v(i,j,k)= ksi_x_v(i,j,k)*(eta_y_v(i,j,k)*phi_z_v(i,j,k)-phi_y_v(i,j,k)*eta_z_v(i,j,k)) &
                             - eta_x_v(i,j,k)*(ksi_y_v(i,j,k)*phi_z_v(i,j,k)-phi_y_v(i,j,k)*ksi_z_v(i,j,k)) &
                             + phi_x_v(i,j,k)*(ksi_y_v(i,j,k)*eta_z_v(i,j,k)-eta_y_v(i,j,k)*ksi_z_v(i,j,k))
          enddo
       enddo
    endif

    ! at jmin: j=0
    if (ndy_v1.ne.1) then
       j=ndy_v1
       do k=1,nz
          do i=1,nx
             ijacob3_v(i,j,k)= ksi_x_v(i,j,k)*(eta_y_v(i,j,k)*phi_z_v(i,j,k)-phi_y_v(i,j,k)*eta_z_v(i,j,k)) &
                             - eta_x_v(i,j,k)*(ksi_y_v(i,j,k)*phi_z_v(i,j,k)-phi_y_v(i,j,k)*ksi_z_v(i,j,k)) &
                             + phi_x_v(i,j,k)*(ksi_y_v(i,j,k)*eta_z_v(i,j,k)-eta_y_v(i,j,k)*ksi_z_v(i,j,k))
          enddo
       enddo
    endif

    ! at jmax: j=ny+1
    if (nfy_v1.ne.ny) then
       j=nfy_v1
       do k=1,nz
          do i=1,nx
             ijacob3_v(i,j,k)= ksi_x_v(i,j,k)*(eta_y_v(i,j,k)*phi_z_v(i,j,k)-phi_y_v(i,j,k)*eta_z_v(i,j,k)) &
                             - eta_x_v(i,j,k)*(ksi_y_v(i,j,k)*phi_z_v(i,j,k)-phi_y_v(i,j,k)*ksi_z_v(i,j,k)) &
                             + phi_x_v(i,j,k)*(ksi_y_v(i,j,k)*eta_z_v(i,j,k)-eta_y_v(i,j,k)*ksi_z_v(i,j,k))
          enddo
       enddo
    endif
    
    ! at kmin: k=0
    if (ndz_v1.ne.1) then
       k=ndz_v1
       do j=1,ny
          do i=1,nx
             ijacob3_v(i,j,k)= ksi_x_v(i,j,k)*(eta_y_v(i,j,k)*phi_z_v(i,j,k)-phi_y_v(i,j,k)*eta_z_v(i,j,k)) &
                             - eta_x_v(i,j,k)*(ksi_y_v(i,j,k)*phi_z_v(i,j,k)-phi_y_v(i,j,k)*ksi_z_v(i,j,k)) &
                             + phi_x_v(i,j,k)*(ksi_y_v(i,j,k)*eta_z_v(i,j,k)-eta_y_v(i,j,k)*ksi_z_v(i,j,k))
          enddo
       enddo
    endif

    ! at jmax: k=nz+1
    if (nfz_v1.ne.nz) then
       k=nfz_v1
       do j=1,ny
          do i=1,nx
             ijacob3_v(i,j,k)= ksi_x_v(i,j,k)*(eta_y_v(i,j,k)*phi_z_v(i,j,k)-phi_y_v(i,j,k)*eta_z_v(i,j,k)) &
                             - eta_x_v(i,j,k)*(ksi_y_v(i,j,k)*phi_z_v(i,j,k)-phi_y_v(i,j,k)*ksi_z_v(i,j,k)) &
                             + phi_x_v(i,j,k)*(ksi_y_v(i,j,k)*eta_z_v(i,j,k)-eta_y_v(i,j,k)*ksi_z_v(i,j,k))
          enddo
       enddo
    endif

    ! Inverse of Jacob3ian
    ! ====================
    ijacob3_v=-1.0_wp/sqrt(abs(ijacob3_v))

  end subroutine grid_metrics_ijacob_3d2_v
  
end submodule smod_grid_metrics3_c3
