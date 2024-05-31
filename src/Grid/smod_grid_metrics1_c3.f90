!=================================================================================
submodule (mod_grid_metrics_c3) smod_grid_metrics1_c3
!=================================================================================

contains

  !===============================================================================
  module subroutine grid_metrics_3d
  !===============================================================================
    !> Compute curvilinear metrics along phi  (full 3D version)
    !> - stencil -ngh:+ngh for inviscid fluxes -
    !===============================================================================
    use mod_mpi ! for: iproc,COMM_global,info
    implicit none
    ! ---------------------------------------------------------------------------
    logical :: isGCL
    real(wp), dimension(:,:,:), allocatable :: xdksi_y,ydksi_z,zdksi_x
    real(wp), dimension(:,:,:), allocatable :: xdeta_y,ydeta_z,zdeta_x
    real(wp), dimension(:,:,:), allocatable :: xdphi_y,ydphi_z,zdphi_x
    real(wp), dimension(:,:,:), allocatable :: ksi_x1,ksi_y1,ksi_z1
    real(wp), dimension(:,:,:), allocatable :: eta_x1,eta_y1,eta_z1
    real(wp), dimension(:,:,:), allocatable :: phi_x1,phi_y1,phi_z1    
    real(wp), dimension(:,:,:), allocatable :: ksi_x2,ksi_y2,ksi_z2
    real(wp), dimension(:,:,:), allocatable :: eta_x2,eta_y2,eta_z2
    real(wp), dimension(:,:,:), allocatable :: phi_x2,phi_y2,phi_z2   
    ! ---------------------------------------------------------------------------
  
    isGCL=.true.
    !isGCL=.false.

    ! Grid derivatives
    ! ================
    ! stored in intermediate 3D array, deleted after 'grid_normals' calculation
    
    ! Allocation & initialization of first derivatives of the grid
    ! ------------------------------------------------------------
    allocate(xdksi(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(ydksi(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(zdksi(nx1:nx2,ny1:ny2,nz1:nz2))
    xdksi=0.0_wp
    ydksi=0.0_wp
    zdksi=0.0_wp
    allocate(xdeta(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(ydeta(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(zdeta(nx1:nx2,ny1:ny2,nz1:nz2))
    xdeta=0.0_wp
    ydeta=0.0_wp
    zdeta=0.0_wp
    allocate(xdphi(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(ydphi(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(zdphi(nx1:nx2,ny1:ny2,nz1:nz2))
    xdphi=0.0_wp
    ydphi=0.0_wp
    zdphi=0.0_wp
    
    ! Reconstruct ghost cells of degenerate edges only for i-direction
    ! ----------------------------------------------------------------
    !call correct_edges_i

    ! Compute first derivatives of the grid in i-direction
    ! ----------------------------------------------------
    call derivative_ksi(xc3,xdksi,yc3,ydksi,zc3,zdksi)

    ! Reconstruct ghost cells of degenerate edges only for j-direction
    ! ----------------------------------------------------------------
    !call correct_edges_j

    ! Compute first derivatives of the grid in j-direction
    ! ----------------------------------------------------
    call derivative_eta(xc3,xdeta,yc3,ydeta,zc3,zdeta)
    
    ! Compute first derivatives of the grid in k-direction
    ! ----------------------------------------------------
    call derivative_phi(xc3,xdphi,yc3,ydphi,zc3,zdphi)
    
    ! Communications of first derivatives
    ! -----------------------------------
    call comm_metrics_3d(xdksi,xdeta,xdphi)
    call comm_metrics_3d(ydksi,ydeta,ydphi)
    call comm_metrics_3d(zdksi,zdeta,zdphi)
    
    call MPI_BARRIER(COMM_global,info)

    ! Correct metrics derivatives due to swap and reverse in neighbors
    ! ----------------------------------------------------------------
    call correct_deriv_sign(xdksi,xdeta,xdphi)
    call correct_deriv_sign(ydksi,ydeta,ydphi)
    call correct_deriv_sign(zdksi,zdeta,zdphi)
        
    ! Initializations
    ! ---------------
    allocate(ksi_x(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(ksi_y(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(ksi_z(nx1:nx2,ny1:ny2,nz1:nz2))
    ksi_x=0.0_wp
    ksi_y=0.0_wp
    ksi_z=0.0_wp
    allocate(eta_x(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(eta_y(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(eta_z(nx1:nx2,ny1:ny2,nz1:nz2))
    eta_x=0.0_wp
    eta_y=0.0_wp
    eta_z=0.0_wp
    allocate(phi_x(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(phi_y(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(phi_z(nx1:nx2,ny1:ny2,nz1:nz2))
    phi_x=0.0_wp
    phi_y=0.0_wp
    phi_z=0.0_wp
    
    allocate(ksi_x1(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(ksi_y1(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(ksi_z1(nx1:nx2,ny1:ny2,nz1:nz2))
    ksi_x1=0.0_wp
    ksi_y1=0.0_wp
    ksi_z1=0.0_wp
    allocate(eta_x1(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(eta_y1(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(eta_z1(nx1:nx2,ny1:ny2,nz1:nz2))
    eta_x1=0.0_wp
    eta_y1=0.0_wp
    eta_z1=0.0_wp
    allocate(phi_x1(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(phi_y1(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(phi_z1(nx1:nx2,ny1:ny2,nz1:nz2))    
    phi_x1=0.0_wp
    phi_y1=0.0_wp
    phi_z1=0.0_wp
    
    allocate(ksi_x2(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(ksi_y2(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(ksi_z2(nx1:nx2,ny1:ny2,nz1:nz2))
    ksi_x2=0.0_wp
    ksi_y2=0.0_wp
    ksi_z2=0.0_wp
    allocate(eta_x2(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(eta_y2(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(eta_z2(nx1:nx2,ny1:ny2,nz1:nz2))
    eta_x2=0.0_wp
    eta_y2=0.0_wp
    eta_z2=0.0_wp
    allocate(phi_x2(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(phi_y2(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(phi_z2(nx1:nx2,ny1:ny2,nz1:nz2))
    phi_x2=0.0_wp
    phi_y2=0.0_wp
    phi_z2=0.0_wp
   
    ! 3D curvilinear metrics
    ! ======================

    if (isGCL) then

       if (iproc==0) print *,'~> compute conservative metrics for inviscid fluxes...'
       if (iproc==0) print *,'   (GCL of Thomas & Lombard)'
       ! GCL along x: ksi_x=(y_eta*z)_phi-(y_phi*z)_eta, ...
       ! GCL along y: ksi_y=(z_eta*x)_phi-(z_phi*x)_eta, ...
       ! GCL along z: ksi_z=(x_eta*y)_phi-(x_phi*y)_eta, ...

       ! Define intermediate arrays
       ! --------------------------
       allocate(xdksi_y(nx1:nx2,ny1:ny2,nz1:nz2))
       allocate(ydksi_z(nx1:nx2,ny1:ny2,nz1:nz2))
       allocate(zdksi_x(nx1:nx2,ny1:ny2,nz1:nz2))
       xdksi_y=xdksi*yc3
       ydksi_z=ydksi*zc3
       zdksi_x=zdksi*xc3

       allocate(xdeta_y(nx1:nx2,ny1:ny2,nz1:nz2))
       allocate(ydeta_z(nx1:nx2,ny1:ny2,nz1:nz2))
       allocate(zdeta_x(nx1:nx2,ny1:ny2,nz1:nz2))
       xdeta_y=xdeta*yc3
       ydeta_z=ydeta*zc3
       zdeta_x=zdeta*xc3
       
       allocate(xdphi_y(nx1:nx2,ny1:ny2,nz1:nz2))
       allocate(ydphi_z(nx1:nx2,ny1:ny2,nz1:nz2))
       allocate(zdphi_x(nx1:nx2,ny1:ny2,nz1:nz2))
       xdphi_y=xdphi*yc3
       ydphi_z=ydphi*zc3
       zdphi_x=zdphi*xc3
                             
       ! compute first part of metrics
       ! -----------------------------
       ! along x: ksi_x=(y_eta*z)_phi, ...
       ! along y: ksi_y=(z_eta*x)_phi, ...
       ! along z: ksi_z=(x_eta*y)_phi, ...
       call derivative_phi(ydeta_z,ksi_x1,zdeta_x,ksi_y1,xdeta_y,ksi_z1)
       call derivative_ksi(ydphi_z,eta_x1,zdphi_x,eta_y1,xdphi_y,eta_z1)
       call derivative_eta(ydksi_z,phi_x1,zdksi_x,phi_y1,xdksi_y,phi_z1)
          
       ! compute second part of metrics
       ! ------------------------------
       ! along x: ksi_x=(y_phi*z)_eta, ...
       ! along y: ksi_y=(z_phi*x)_eta, ...
       ! along z: ksi_z=(x_phi*y)_eta, ...
       call derivative_eta(ydphi_z,ksi_x2,zdphi_x,ksi_y2,xdphi_y,ksi_z2)
       call derivative_phi(ydksi_z,eta_x2,zdksi_x,eta_y2,xdksi_y,eta_z2)
       call derivative_ksi(ydeta_z,phi_x2,zdeta_x,phi_y2,xdeta_y,phi_z2)

       ! communicate metrics
       ! -------------------
       call comm_metrics_3de(ksi_x1,eta_x1,phi_x1,ksi_x2,eta_x2,phi_x2)
       call comm_metrics_3de(ksi_y1,eta_y1,phi_y1,ksi_y2,eta_y2,phi_y2)
       call comm_metrics_3de(ksi_z1,eta_z1,phi_z1,ksi_z2,eta_z2,phi_z2)
       call MPI_BARRIER(COMM_global,info)
       
       ! correct signs due to reverse
       ! ----------------------------
       call correct_dderiv_sign(ksi_x1,eta_x1,phi_x1,ksi_x2,eta_x2,phi_x2)
       call correct_dderiv_sign(ksi_y1,eta_y1,phi_y1,ksi_y2,eta_y2,phi_y2)
       call correct_dderiv_sign(ksi_z1,eta_z1,phi_z1,ksi_z2,eta_z2,phi_z2)
       
       ! form final metrics
       ! ------------------
       ksi_x = ksi_x1-ksi_x2
       ksi_y = ksi_y1-ksi_y2
       ksi_z = ksi_z1-ksi_z2
       eta_x = eta_x1-eta_x2
       eta_y = eta_y1-eta_y2
       eta_z = eta_z1-eta_z2
       phi_x = phi_x1-phi_x2
       phi_y = phi_y1-phi_y2
       phi_z = phi_z1-phi_z2
       
       ! Free memory
       ! -----------
       deallocate(xdksi_y,ydksi_z,zdksi_x)
       deallocate(xdeta_y,ydeta_z,zdeta_x)
       deallocate(xdphi_y,ydphi_z,zdphi_x)

       deallocate(ksi_x1,ksi_y1,ksi_z1)
       deallocate(eta_x1,eta_y1,eta_z1)
       deallocate(phi_x1,phi_y1,phi_z1)   
       deallocate(ksi_x2,ksi_y2,ksi_z2)
       deallocate(eta_x2,eta_y2,eta_z2)
       deallocate(phi_x2,phi_y2,phi_z2)  
    else

       if (iproc==0) print *,'~> compute non-conservative metrics for inviscid fluxes...'

       ksi_x = ydeta*zdphi-ydphi*zdeta
       ksi_y = zdeta*xdphi-zdphi*xdeta
       ksi_z = xdeta*ydphi-xdphi*ydeta

       eta_x = ydphi*zdksi-ydksi*zdphi
       eta_y = zdphi*xdksi-zdksi*xdphi
       eta_z = xdphi*ydksi-xdksi*ydphi

       phi_x = ydksi*zdeta-ydeta*zdksi
       phi_y = zdksi*xdeta-zdeta*xdksi
       phi_z = xdksi*ydeta-xdeta*ydksi
       
   endif

    ! inverse Jacobian
    ! -----------------
    ! for inviscid fluxes (ngh ghost cells)
    call grid_metrics_ijacob_3d
       
    ! 3D curvilinear metrics for viscous fluxes (ngh_v ghost cells)
    ! =============================================================
    if ((iorder_visc.ne.0).or.(idepart==POST_PROCESSING)) then
       
       ! Allocation & initialization of first derivatives of the grid
       ! ------------------------------------------------------------
       deallocate(xdksi,ydksi,zdksi)
       allocate(xdksi(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
       allocate(ydksi(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
       allocate(zdksi(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
       xdksi=0.0_wp
       ydksi=0.0_wp
       zdksi=0.0_wp
       deallocate(xdeta,ydeta,zdeta)
       allocate(xdeta(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
       allocate(ydeta(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
       allocate(zdeta(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
       xdeta=0.0_wp
       ydeta=0.0_wp
       zdeta=0.0_wp
       deallocate(xdphi,ydphi,zdphi)
       allocate(xdphi(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
       allocate(ydphi(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
       allocate(zdphi(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
       xdphi=0.0_wp
       ydphi=0.0_wp
       zdphi=0.0_wp

       ! Compute first derivatives of the grid
       ! -------------------------------------
       call derivative_ksi_v(xc3(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v),xdksi, &
                             yc3(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v),ydksi, &
                             zc3(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v),zdksi)
       call derivative_eta_v(xc3(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v),xdeta, &
                             yc3(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v),ydeta, &
                             zc3(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v),zdeta)
       call derivative_phi_v(xc3(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v),xdphi, &
                             yc3(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v),ydphi, &
                             zc3(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v),zdphi)

       ! Communications of first derivatives
       ! -----------------------------------
       call comm_metrics_3d_v(xdksi,xdeta,xdphi)
       call comm_metrics_3d_v(ydksi,ydeta,ydphi)
       call comm_metrics_3d_v(zdksi,zdeta,zdphi)
       call MPI_BARRIER(COMM_global,info)
       
       ! Correct metrics derivatives due to swap and reverse in neighbors
       ! ----------------------------------------------------------------
       call correct_deriv_sign_v(xdksi,xdeta,xdphi)
       call correct_deriv_sign_v(ydksi,ydeta,ydphi)
       call correct_deriv_sign_v(zdksi,zdeta,zdphi)

       ! Initializations
       ! ---------------
       allocate(ksi_x_v(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
       allocate(ksi_y_v(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
       allocate(ksi_z_v(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
       ksi_x_v=0.0_wp
       ksi_y_v=0.0_wp
       ksi_z_v=0.0_wp
       allocate(eta_x_v(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
       allocate(eta_y_v(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
       allocate(eta_z_v(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
       eta_x_v=0.0_wp
       eta_y_v=0.0_wp
       eta_z_v=0.0_wp
       allocate(phi_x_v(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
       allocate(phi_y_v(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
       allocate(phi_z_v(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
       phi_x_v=0.0_wp
       phi_y_v=0.0_wp
       phi_z_v=0.0_wp

       !isGCL=.false.

       if (isGCL) then
          if (iproc==0) print *,'~> compute conservative metrics for viscous fluxes...'
          ! GCL along x: ksi_x=(y_eta*z)_phi-(y_phi*z)_eta, ...
          ! GCL along y: ksi_y=(z_eta*x)_phi-(z_phi*x)_eta, ...
          ! GCL along z: ksi_z=(x_eta*y)_phi-(x_phi*y)_eta, ...

          ! Define intermediate arrays
          ! --------------------------
          allocate(xdksi_y(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
          allocate(ydksi_z(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
          allocate(zdksi_x(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
          xdksi_y=xdksi*yc3(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v)
          ydksi_z=ydksi*zc3(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v)
          zdksi_x=zdksi*xc3(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v)

          allocate(xdeta_y(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
          allocate(ydeta_z(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
          allocate(zdeta_x(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
          xdeta_y=xdeta*yc3(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v)
          ydeta_z=ydeta*zc3(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v)
          zdeta_x=zdeta*xc3(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v)

          allocate(xdphi_y(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
          allocate(ydphi_z(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
          allocate(zdphi_x(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
          xdphi_y=xdphi*yc3(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v)
          ydphi_z=ydphi*zc3(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v)
          zdphi_x=zdphi*xc3(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v)

          allocate(ksi_x1(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
          allocate(ksi_y1(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
          allocate(ksi_z1(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
          ksi_x1=0.0_wp
          ksi_y1=0.0_wp
          ksi_z1=0.0_wp
          allocate(eta_x1(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
          allocate(eta_y1(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
          allocate(eta_z1(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
          eta_x1=0.0_wp
          eta_y1=0.0_wp
          eta_z1=0.0_wp
          allocate(phi_x1(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
          allocate(phi_y1(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
          allocate(phi_z1(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))    
          phi_x1=0.0_wp
          phi_y1=0.0_wp
          phi_z1=0.0_wp

          allocate(ksi_x2(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
          allocate(ksi_y2(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
          allocate(ksi_z2(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
          ksi_x2=0.0_wp
          ksi_y2=0.0_wp
          ksi_z2=0.0_wp
          allocate(eta_x2(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
          allocate(eta_y2(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
          allocate(eta_z2(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
          eta_x2=0.0_wp
          eta_y2=0.0_wp
          eta_z2=0.0_wp
          allocate(phi_x2(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
          allocate(phi_y2(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
          allocate(phi_z2(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
          phi_x2=0.0_wp
          phi_y2=0.0_wp
          phi_z2=0.0_wp

          ! compute first part of metrics
          ! -----------------------------
          ! along x: ksi_x=(y_eta*z)_phi, ...
          ! along y: ksi_y=(z_eta*x)_phi, ...
          ! along z: ksi_z=(x_eta*y)_phi, ...
          call derivative_phi_v(ydeta_z,ksi_x1,zdeta_x,ksi_y1,xdeta_y,ksi_z1)
          call derivative_ksi_v(ydphi_z,eta_x1,zdphi_x,eta_y1,xdphi_y,eta_z1)
          call derivative_eta_v(ydksi_z,phi_x1,zdksi_x,phi_y1,xdksi_y,phi_z1)

          ! compute second part of metrics
          ! ------------------------------
          ! along x: ksi_x=(y_phi*z)_eta, ...
          ! along y: ksi_y=(z_phi*x)_eta, ...
          ! along z: ksi_z=(x_phi*y)_eta, ...
          call derivative_eta_v(ydphi_z,ksi_x2,zdphi_x,ksi_y2,xdphi_y,ksi_z2)
          call derivative_phi_v(ydksi_z,eta_x2,zdksi_x,eta_y2,xdksi_y,eta_z2)
          call derivative_ksi_v(ydeta_z,phi_x2,zdeta_x,phi_y2,xdeta_y,phi_z2)

          ! communicate metrics
          ! -------------------
          call comm_metrics_3de_v(ksi_x1,eta_x1,phi_x1,ksi_x2,eta_x2,phi_x2)
          call comm_metrics_3de_v(ksi_y1,eta_y1,phi_y1,ksi_y2,eta_y2,phi_y2)
          call comm_metrics_3de_v(ksi_z1,eta_z1,phi_z1,ksi_z2,eta_z2,phi_z2)
          call MPI_BARRIER(COMM_global,info)

          ! correct signs due to reverse
          ! ----------------------------
          call correct_dderiv_sign_v(ksi_x1,eta_x1,phi_x1,ksi_x2,eta_x2,phi_x2)
          call correct_dderiv_sign_v(ksi_y1,eta_y1,phi_y1,ksi_y2,eta_y2,phi_y2)
          call correct_dderiv_sign_v(ksi_z1,eta_z1,phi_z1,ksi_z2,eta_z2,phi_z2)

          ! form final metrics
          ! ------------------
          ksi_x_v = ksi_x1-ksi_x2
          ksi_y_v = ksi_y1-ksi_y2
          ksi_z_v = ksi_z1-ksi_z2
          eta_x_v = eta_x1-eta_x2
          eta_y_v = eta_y1-eta_y2
          eta_z_v = eta_z1-eta_z2
          phi_x_v = phi_x1-phi_x2
          phi_y_v = phi_y1-phi_y2
          phi_z_v = phi_z1-phi_z2
       
          ! Free memory
          ! -----------
          deallocate(xdksi_y,ydksi_z,zdksi_x)
          deallocate(xdeta_y,ydeta_z,zdeta_x)
          deallocate(xdphi_y,ydphi_z,zdphi_x)

          deallocate(ksi_x1,ksi_y1,ksi_z1)
          deallocate(eta_x1,eta_y1,eta_z1)
          deallocate(phi_x1,phi_y1,phi_z1)   
          deallocate(ksi_x2,ksi_y2,ksi_z2)
          deallocate(eta_x2,eta_y2,eta_z2)
          deallocate(phi_x2,phi_y2,phi_z2)  
       else

          if (iproc==0) print *,'~> compute non-conservative metrics for viscous fluxes...'

          ksi_x_v = ydeta*zdphi-ydphi*zdeta
          ksi_y_v = zdeta*xdphi-zdphi*xdeta
          ksi_z_v = xdeta*ydphi-xdphi*ydeta

          eta_x_v = ydphi*zdksi-ydksi*zdphi
          eta_y_v = zdphi*xdksi-zdksi*xdphi
          eta_z_v = xdphi*ydksi-xdksi*ydphi

          phi_x_v = ydksi*zdeta-ydeta*zdksi
          phi_y_v = zdksi*xdeta-zdeta*xdksi
          phi_z_v = xdksi*ydeta-xdeta*ydksi
          
       endif

       ! inverse Jacobian
       ! -----------------
       ! for viscous fluxes (ngh_v ghost cells)
       call grid_metrics_ijacob_3d_v

    endif

  end subroutine grid_metrics_3d

end submodule smod_grid_metrics1_c3
