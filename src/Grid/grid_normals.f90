!==============================================================================
subroutine grid_normals
!==============================================================================
  !> Computation of wall normals and BC parameters
  !> - curvilinear grid -
!==============================================================================
  use mod_flow ! for grid,BC, ...
  use mod_grid_metrics_c3 ! for xdksi,ydksi,zdksi,xdeta,ydeta,zdeta,xdphi,ydphi,zdphi
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j,k
  real(wp) :: n2_eta,n2_ksi,n2_phi
  real(wp) :: gksi,geta,gphi,eeta_ephi,eksi_ephi,eksi_eeta,fe
  ! sign correction to have inwards wall normals
  real(wp) :: din,di(3)
  ! ---------------------------------------------------------------------------

  if (is_curv3) then

     ! Determination of normals for BC imin (ksi=cste)
     ! ====================================
     if (BC_face(1,1)%sort.le.0) then     

        ! allocate arrays for BC normals
        ! ------------------------------
        allocate(nxn_imin(ny,nz),nyn_imin(ny,nz),nzn_imin(ny,nz))
        allocate(getagksi3_imin(ny,nz),gphigksi3_imin(ny,nz))

        ! boundary normals along ksi
        ! --------------------------
        do k=1,nz
           do j=1,ny
              ! J²*||grad(ksi)||²
              n2_ksi=ksi_x(1,j,k)**2+ksi_y(1,j,k)**2+ksi_z(1,j,k)**2
              ! BC normal components: J*grad(ksi)/J*||grad(ksi)||
              nxn_imin(j,k) = ksi_x(1,j,k)/sqrt(n2_ksi)
              nyn_imin(j,k) = ksi_y(1,j,k)/sqrt(n2_ksi)
              nzn_imin(j,k) = ksi_z(1,j,k)/sqrt(n2_ksi)
              ! compute getagksi3 for BC: -grad(eta).grad(ksi)/||grad(ksi)||²
              getagksi3_imin(j,k)=eta_x(1,j,k)*ksi_x(1,j,k)+eta_y(1,j,k)*ksi_y(1,j,k)+eta_z(1,j,k)*ksi_z(1,j,k)
              getagksi3_imin(j,k)=-getagksi3_imin(j,k)/n2_ksi
              if (abs(getagksi3_imin(j,k))<1.e-16_wp) getagksi3_imin(j,k)=0.0_wp
              ! compute gphigksi3 for BC: -grad(phi).grad(ksi)/||grad(ksi)||²
              gphigksi3_imin(j,k)=phi_x(1,j,k)*ksi_x(1,j,k)+phi_y(1,j,k)*ksi_y(1,j,k)+phi_z(1,j,k)*ksi_z(1,j,k)
              gphigksi3_imin(j,k)=-gphigksi3_imin(j,k)/n2_ksi
              if (abs(gphigksi3_imin(j,k))<1.e-16_wp) gphigksi3_imin(j,k)=0.0_wp
           enddo
        enddo

        if (BC_face(1,1)%sort==0) then

           ! Projection on wall tangent
           ! ==========================
           ! 1/ to enforce adiabaticity in grad_T_*_c3.f90
           ! 2/ to suppress wall-normal velocity for Eulerian walls in mod_bc_wall_2_c3.f90
           
           ! allocations
           ! -----------
           allocate(tx_eta_imin(ny1_v:ny2_v,nz1_v:nz2_v),ty_eta_imin(ny1_v:ny2_v,nz1_v:nz2_v),tz_eta_imin(ny1_v:ny2_v,nz1_v:nz2_v))
           allocate(tx_phi_imin(ny1_v:ny2_v,nz1_v:nz2_v),ty_phi_imin(ny1_v:ny2_v,nz1_v:nz2_v),tz_phi_imin(ny1_v:ny2_v,nz1_v:nz2_v))
           allocate(txn_eta_imin(ny1_v:ny2_v,nz1_v:nz2_v),tyn_eta_imin(ny1_v:ny2_v,nz1_v:nz2_v),tzn_eta_imin(ny1_v:ny2_v,nz1_v:nz2_v))
           allocate(txn_phi_imin(ny1_v:ny2_v,nz1_v:nz2_v),tyn_phi_imin(ny1_v:ny2_v,nz1_v:nz2_v),tzn_phi_imin(ny1_v:ny2_v,nz1_v:nz2_v))

           ! tangent components & projection on tangent
           ! ------------------------------------------
           do k=ndz_v1,nfz_v1
              do j=ndy_v1,nfy_v1
                 geta=sqrt(xdeta(1,j,k)**2+ydeta(1,j,k)**2+zdeta(1,j,k)**2)
                 if (geta==0.0_wp) geta=1.0_wp
                 gphi=sqrt(xdphi(1,j,k)**2+ydphi(1,j,k)**2+zdphi(1,j,k)**2)
                 if (gphi==0.0_wp) gphi=1.0_wp
                 eeta_ephi=xdeta(1,j,k)*xdphi(1,j,k)+ydeta(1,j,k)*ydphi(1,j,k)+zdeta(1,j,k)*zdphi(1,j,k)
                 eeta_ephi=eeta_ephi/geta/gphi
                 fe=1.0_wp/(1.0_wp-eeta_ephi**2)                
                 tx_eta_imin(j,k)=fe*(xdeta(1,j,k)/geta-eeta_ephi*xdphi(1,j,k)/gphi)
                 ty_eta_imin(j,k)=fe*(ydeta(1,j,k)/geta-eeta_ephi*ydphi(1,j,k)/gphi)
                 tz_eta_imin(j,k)=fe*(zdeta(1,j,k)/geta-eeta_ephi*zdphi(1,j,k)/gphi)
                 tx_phi_imin(j,k)=fe*(xdphi(1,j,k)/gphi-eeta_ephi*xdeta(1,j,k)/geta)
                 ty_phi_imin(j,k)=fe*(ydphi(1,j,k)/gphi-eeta_ephi*ydeta(1,j,k)/geta)
                 tz_phi_imin(j,k)=fe*(zdphi(1,j,k)/gphi-eeta_ephi*zdeta(1,j,k)/geta)
                 txn_eta_imin(j,k)=xdeta(1,j,k)/geta
                 tyn_eta_imin(j,k)=ydeta(1,j,k)/geta
                 tzn_eta_imin(j,k)=zdeta(1,j,k)/geta
                 txn_phi_imin(j,k)=xdphi(1,j,k)/gphi
                 tyn_phi_imin(j,k)=ydphi(1,j,k)/gphi
                 tzn_phi_imin(j,k)=zdphi(1,j,k)/gphi
              enddo
           enddo
           
           ! Oriented surface elements
           ! =========================
           ! 1/ to compute drag and lift forces
           
           ! allocate
           ! --------
           allocate(nxnds_imin(ny,nz),nynds_imin(ny,nz),nznds_imin(ny,nz),inksi_imin(ny,nz))

           ! determine correct_sign of wall normals
           ! --------------------------------------
           ! /!\ direction depends on (i,j) orientation
           !     NEEDS to set INWARDS normals for inlet

           ! vector pointing inwards (between first and second cells)
           di(1)=xc3(2,1,1)-xc3(1,1,1)
           di(2)=yc3(2,1,1)-yc3(1,1,1)
           di(3)=zc3(2,1,1)-zc3(1,1,1)
           ! scalar product with wall normal
           din=di(1)*nxn_imin(1,1)+di(2)*nyn_imin(1,1)+di(3)*nzn_imin(1,1)
           ! sign
           sgn_imin=-din/abs(din)

           ! oriented surface elements
           do k=1,nz
              do j=1,ny
                 inksi_imin(j,k)=1.0_wp/sqrt(ksi_x(1,j,k)**2+ksi_y(1,j,k)**2+ksi_z(1,j,k)**2)
                 nxnds_imin(j,k) = ksi_x(1,j,k)*sgn_imin
                 nynds_imin(j,k) = ksi_y(1,j,k)*sgn_imin
                 nznds_imin(j,k) = ksi_z(1,j,k)*sgn_imin
              enddo
           enddo

        endif

     endif

     ! Determination of normals for BC imax (ksi=cste)
     ! ====================================
     if (BC_face(1,2)%sort.le.0) then     

        ! allocate arrays for BC normals
        ! ------------------------------
        allocate(nxn_imax(ny,nz),nyn_imax(ny,nz),nzn_imax(ny,nz))
        allocate(getagksi3_imax(ny,nz),gphigksi3_imax(ny,nz))

        ! boundary along ksi
        ! ------------------
        do k=1,nz
           do j=1,ny
              ! J²*||grad(ksi)||²
              n2_ksi=ksi_x(nx,j,k)**2+ksi_y(nx,j,k)**2+ksi_z(nx,j,k)**2
              ! BC normal components: J*grad(ksi)/J*||grad(ksi)||
              nxn_imax(j,k) = ksi_x(nx,j,k)/sqrt(n2_ksi)
              nyn_imax(j,k) = ksi_y(nx,j,k)/sqrt(n2_ksi)
              nzn_imax(j,k) = ksi_z(nx,j,k)/sqrt(n2_ksi)
              ! compute getagksi3 for BC: -grad(eta).grad(ksi)/||grad(ksi)||²
              getagksi3_imax(j,k)=eta_x(nx,j,k)*ksi_x(nx,j,k)+eta_y(nx,j,k)*ksi_y(nx,j,k)+eta_z(nx,j,k)*ksi_z(nx,j,k)
              getagksi3_imax(j,k)=-getagksi3_imax(j,k)/n2_ksi
              if (abs(getagksi3_imax(j,k))<1.e-16_wp) getagksi3_imax(j,k)=0.0_wp
              ! compute gphigksi3 for BC: -grad(phi).grad(ksi)/||grad(ksi)||²
              gphigksi3_imax(j,k)=phi_x(nx,j,k)*ksi_x(nx,j,k)+phi_y(nx,j,k)*ksi_y(nx,j,k)+phi_z(nx,j,k)*ksi_z(nx,j,k)
              gphigksi3_imax(j,k)=-gphigksi3_imax(j,k)/n2_ksi
              if (abs(gphigksi3_imax(j,k))<1.e-16_wp) gphigksi3_imax(j,k)=0.0_wp
           enddo
        enddo
        
        if (BC_face(1,2)%sort==0) then

           ! Projection on wall tangent
           ! ==========================
           ! 1/ to enforce adiabaticity in grad_T_*_c3.f90
           ! 2/ to suppress wall-normal velocity for Eulerian walls in mod_bc_wall_2_c3.f90
           
           ! allocations
           ! -----------
           allocate(tx_eta_imax(ny1_v:ny2_v,nz1_v:nz2_v),ty_eta_imax(ny1_v:ny2_v,nz1_v:nz2_v),tz_eta_imax(ny1_v:ny2_v,nz1_v:nz2_v))
           allocate(tx_phi_imax(ny1_v:ny2_v,nz1_v:nz2_v),ty_phi_imax(ny1_v:ny2_v,nz1_v:nz2_v),tz_phi_imax(ny1_v:ny2_v,nz1_v:nz2_v))
           allocate(txn_eta_imax(ny1_v:ny2_v,nz1_v:nz2_v),tyn_eta_imax(ny1_v:ny2_v,nz1_v:nz2_v),tzn_eta_imax(ny1_v:ny2_v,nz1_v:nz2_v))
           allocate(txn_phi_imax(ny1_v:ny2_v,nz1_v:nz2_v),tyn_phi_imax(ny1_v:ny2_v,nz1_v:nz2_v),tzn_phi_imax(ny1_v:ny2_v,nz1_v:nz2_v))

           ! tangent components & projection on tangent
           ! ------------------------------------------
           do k=ndz_v1,nfz_v1
              do j=ndy_v1,nfy_v1
                 geta=sqrt(xdeta(nx,j,k)**2+ydeta(nx,j,k)**2+zdeta(nx,j,k)**2)
                 if (geta==0.0_wp) geta=1.0_wp
                 gphi=sqrt(xdphi(nx,j,k)**2+ydphi(nx,j,k)**2+zdphi(nx,j,k)**2)
                 if (gphi==0.0_wp) gphi=1.0_wp
                 eeta_ephi=xdeta(nx,j,k)*xdphi(nx,j,k)+ydeta(nx,j,k)*ydphi(nx,j,k)+zdeta(nx,j,k)*zdphi(nx,j,k)
                 eeta_ephi=eeta_ephi/geta/gphi
                 fe=1.0_wp/(1.0_wp-eeta_ephi**2)                
                 tx_eta_imax(j,k)=fe*(xdeta(nx,j,k)/geta-eeta_ephi*xdphi(nx,j,k)/gphi)
                 ty_eta_imax(j,k)=fe*(ydeta(nx,j,k)/geta-eeta_ephi*ydphi(nx,j,k)/gphi)
                 tz_eta_imax(j,k)=fe*(zdeta(nx,j,k)/geta-eeta_ephi*zdphi(nx,j,k)/gphi)
                 tx_phi_imax(j,k)=fe*(xdphi(nx,j,k)/gphi-eeta_ephi*xdeta(nx,j,k)/geta)
                 ty_phi_imax(j,k)=fe*(ydphi(nx,j,k)/gphi-eeta_ephi*ydeta(nx,j,k)/geta)
                 tz_phi_imax(j,k)=fe*(zdphi(nx,j,k)/gphi-eeta_ephi*zdeta(nx,j,k)/geta)
                 txn_eta_imax(j,k)=xdeta(nx,j,k)/geta
                 tyn_eta_imax(j,k)=ydeta(nx,j,k)/geta
                 tzn_eta_imax(j,k)=zdeta(nx,j,k)/geta
                 txn_phi_imax(j,k)=xdphi(nx,j,k)/gphi
                 tyn_phi_imax(j,k)=ydphi(nx,j,k)/gphi
                 tzn_phi_imax(j,k)=zdphi(nx,j,k)/gphi
              enddo
           enddo
           
           ! Oriented surface elements
           ! =========================
           ! 1/ to compute drag and lift forces
           
           ! allocate
           ! --------
           allocate(nxnds_imax(ny,nz),nynds_imax(ny,nz),nznds_imax(ny,nz),inksi_imax(ny,nz))

           ! determine correct_sign of wall normals
           ! --------------------------------------
           ! /!\ direction depends on (i,j) orientation
           !     NEEDS to set INWARDS normals for inlet

           ! vector pointing inwards (between first and second cells)
           di(1)=xc3(nx,1,1)-xc3(nx-1,1,1)
           di(2)=yc3(nx,1,1)-yc3(nx-1,1,1)
           di(3)=zc3(nx,1,1)-zc3(nx-1,1,1)
           ! scalar product with wall normal
           din=di(1)*nxn_imax(1,1)+di(2)*nyn_imax(1,1)+di(3)*nzn_imax(1,1)
           ! sign
           sgn_imax=din/abs(din)

           ! oriented surface elements
           do k=1,nz
              do j=1,ny
                 inksi_imax(j,k)=1.0_wp/sqrt(ksi_x(nx,j,k)**2+ksi_y(nx,j,k)**2+ksi_z(nx,j,k)**2)
                 nxnds_imax(j,k) = ksi_x(nx,j,k)*sgn_imax
                 nynds_imax(j,k) = ksi_y(nx,j,k)*sgn_imax
                 nznds_imax(j,k) = ksi_z(nx,j,k)*sgn_imax
              enddo
           enddo

        endif

     endif

     ! Determination of normals for BC jmin (eta=cste)
     ! ====================================
     if (BC_face(2,1)%sort.le.0) then     

        ! allocate arrays for BC normals
        ! ------------------------------
        allocate(nxn_jmin(nx,nz),nyn_jmin(nx,nz),nzn_jmin(nx,nz))
        allocate(gksigeta3_jmin(nx,nz),gphigeta3_jmin(nx,nz))

        ! boundary along eta
        ! ------------------
        do k=1,nz
           do i=1,nx
              ! J²*||grad(eta)||²
              n2_eta=eta_x(i,1,k)**2+eta_y(i,1,k)**2+eta_z(i,1,k)**2
              ! BC normal components: J*grad(eta)/J*||grad(eta)||
              nxn_jmin(i,k) = eta_x(i,1,k)/sqrt(n2_eta)
              nyn_jmin(i,k) = eta_y(i,1,k)/sqrt(n2_eta)
              nzn_jmin(i,k) = eta_z(i,1,k)/sqrt(n2_eta)
              ! compute gksigeta3 for BC: -grad(ksi).grad(eta)/||grad(eta)||²
              gksigeta3_jmin(i,k)=ksi_x(i,1,k)*eta_x(i,1,k)+ksi_y(i,1,k)*eta_y(i,1,k)+ksi_z(i,1,k)*eta_z(i,1,k)
              gksigeta3_jmin(i,k)=-gksigeta3_jmin(i,k)/n2_eta
              if (abs(gksigeta3_jmin(i,k))<1.e-16_wp) gksigeta3_jmin(i,k)=0.0_wp
              ! compute gphigeta3 for BC: -grad(phi).grad(eta)/||grad(eta)||²
              gphigeta3_jmin(i,k)=phi_x(i,1,k)*eta_x(i,1,k)+phi_y(i,1,k)*eta_y(i,1,k)+phi_z(i,1,k)*eta_z(i,1,k)
              gphigeta3_jmin(i,k)=-gphigeta3_jmin(i,k)/n2_eta
              if (abs(gphigeta3_jmin(i,k))<1.e-16_wp) gphigeta3_jmin(i,k)=0.0_wp
           enddo
        enddo

        if (BC_face(2,1)%sort==0) then

           ! Projection on wall tangent
           ! ==========================
           ! 1/ to enforce adiabaticity in grad_T_*_c3.f90
           ! 2/ to suppress wall-normal velocity for Eulerian walls in mod_bc_wall_2_c3.f90
           
           ! allocations
           ! -----------
           allocate(tx_ksi_jmin(nx1_v:nx2_v,nz1_v:nz2_v),ty_ksi_jmin(nx1_v:nx2_v,nz1_v:nz2_v),tz_ksi_jmin(nx1_v:nx2_v,nz1_v:nz2_v))
           allocate(tx_phi_jmin(nx1_v:nx2_v,nz1_v:nz2_v),ty_phi_jmin(nx1_v:nx2_v,nz1_v:nz2_v),tz_phi_jmin(nx1_v:nx2_v,nz1_v:nz2_v))
           allocate(txn_ksi_jmin(nx1_v:nx2_v,nz1_v:nz2_v),tyn_ksi_jmin(nx1_v:nx2_v,nz1_v:nz2_v),tzn_ksi_jmin(nx1_v:nx2_v,nz1_v:nz2_v))
           allocate(txn_phi_jmin(nx1_v:nx2_v,nz1_v:nz2_v),tyn_phi_jmin(nx1_v:nx2_v,nz1_v:nz2_v),tzn_phi_jmin(nx1_v:nx2_v,nz1_v:nz2_v))

           ! tangent components & projection on tangent
           ! ------------------------------------------
           do k=ndz_v1,nfz_v1
              do i=ndx_v1,nfx_v1
                 gksi=sqrt(xdksi(i,1,k)**2+ydksi(i,1,k)**2+zdksi(i,1,k)**2)
                 if (gksi==0.0_wp) gksi=1.0_wp
                 gphi=sqrt(xdphi(i,1,k)**2+ydphi(i,1,k)**2+zdphi(i,1,k)**2)
                 if (gphi==0.0_wp) gphi=1.0_wp
                 eksi_ephi=xdksi(i,1,k)*xdphi(i,1,k)+ydksi(i,1,k)*ydphi(i,1,k)+zdksi(i,1,k)*zdphi(i,1,k)
                 eksi_ephi=eksi_ephi/gksi/gphi
                 fe=1.0_wp/(1.0_wp-eksi_ephi**2)                
                 tx_ksi_jmin(i,k)=fe*(xdksi(i,1,k)/gksi-eksi_ephi*xdphi(i,1,k)/gphi)
                 ty_ksi_jmin(i,k)=fe*(ydksi(i,1,k)/gksi-eksi_ephi*ydphi(i,1,k)/gphi)
                 tz_ksi_jmin(i,k)=fe*(zdksi(i,1,k)/gksi-eksi_ephi*zdphi(i,1,k)/gphi)
                 tx_phi_jmin(i,k)=fe*(xdphi(i,1,k)/gphi-eksi_ephi*xdksi(i,1,k)/gksi)
                 ty_phi_jmin(i,k)=fe*(ydphi(i,1,k)/gphi-eksi_ephi*ydksi(i,1,k)/gksi)
                 tz_phi_jmin(i,k)=fe*(zdphi(i,1,k)/gphi-eksi_ephi*zdksi(i,1,k)/gksi)
                 txn_ksi_jmin(i,k)=xdksi(i,1,k)/gksi
                 tyn_ksi_jmin(i,k)=ydksi(i,1,k)/gksi
                 tzn_ksi_jmin(i,k)=zdksi(i,1,k)/gksi
                 txn_phi_jmin(i,k)=xdphi(i,1,k)/gphi
                 tyn_phi_jmin(i,k)=ydphi(i,1,k)/gphi
                 tzn_phi_jmin(i,k)=zdphi(i,1,k)/gphi
              enddo
           enddo
           
           ! Oriented surface elements
           ! =========================
           ! 1/ to compute drag and lift forces
           
           ! allocate
           ! --------
           allocate(nxnds_jmin(nx,nz),nynds_jmin(nx,nz),nznds_jmin(nx,nz),ineta_jmin(nx,nz))

           ! determine correct_sign of wall normals
           ! --------------------------------------
           ! /!\ direction depends on (i,j) orientation
           !     NEEDS to set INWARDS normals for inlet

           ! vector pointing inwards (between first and second cells)
           di(1)=xc3(1,2,1)-xc3(1,1,1)
           di(2)=yc3(1,2,1)-yc3(1,1,1)
           di(3)=zc3(1,2,1)-zc3(1,1,1)
           ! scalar product with wall normal
           din=di(1)*nxn_jmin(1,1)+di(2)*nyn_jmin(1,1)+di(3)*nzn_jmin(1,1)
           ! sign
           sgn_jmin=-din/abs(din)

           ! oriented surface elements
           do k=1,nz
              do i=1,nx
                 ineta_jmin(i,k)=1.0_wp/sqrt(eta_x(i,1,k)**2+eta_y(i,1,k)**2+eta_z(i,1,k)**2)
                 nxnds_jmin(i,k) = eta_x(i,1,k)*sgn_jmin
                 nynds_jmin(i,k) = eta_y(i,1,k)*sgn_jmin
                 nznds_jmin(i,k) = eta_z(i,1,k)*sgn_jmin
              enddo
           enddo

        endif

     endif

     ! Determination of normals for BC jmax (eta=cste)
     ! ====================================
     if (BC_face(2,2)%sort.le.0) then     

        ! allocate arrays for BC normals
        ! ------------------------------
        allocate(nxn_jmax(nx,nz),nyn_jmax(nx,nz),nzn_jmax(nx,nz))
        allocate(gksigeta3_jmax(nx,nz),gphigeta3_jmax(nx,nz))

        ! boundary along eta
        ! ------------------
        do k=1,nz
           do i=1,nx
              ! J²*||grad(eta)||²
              n2_eta=eta_x(i,ny,k)**2+eta_y(i,ny,k)**2+eta_z(i,ny,k)**2
              ! BC normal components: J*grad(eta)/J*||grad(eta)||
              nxn_jmax(i,k) = eta_x(i,ny,k)/sqrt(n2_eta)
              nyn_jmax(i,k) = eta_y(i,ny,k)/sqrt(n2_eta)
              nzn_jmax(i,k) = eta_z(i,ny,k)/sqrt(n2_eta)
              ! compute gksigeta3 for BC: -grad(ksi).grad(eta)/||grad(eta)||²
              gksigeta3_jmax(i,k)=ksi_x(i,ny,k)*eta_x(i,ny,k)+ksi_y(i,ny,k)*eta_y(i,ny,k)+ksi_z(i,ny,k)*eta_z(i,ny,k)
              gksigeta3_jmax(i,k)=-gksigeta3_jmax(i,k)/n2_eta
              if (abs(gksigeta3_jmax(i,k))<1.e-16_wp) gksigeta3_jmax(i,k)=0.0_wp
              ! compute gphigeta3 for BC: -grad(phi).grad(eta)/||grad(eta)||²
              gphigeta3_jmax(i,k)=phi_x(i,ny,k)*eta_x(i,ny,k)+phi_y(i,ny,k)*eta_y(i,ny,k)+phi_z(i,ny,k)*eta_z(i,ny,k)
              gphigeta3_jmax(i,k)=-gphigeta3_jmax(i,k)/n2_eta
              if (abs(gphigeta3_jmax(i,k))<1.e-16_wp) gphigeta3_jmax(i,k)=0.0_wp
           enddo
        enddo

        if (BC_face(2,2)%sort==0) then

           ! Projection on wall tangent
           ! ==========================
           ! 1/ to enforce adiabaticity in grad_T_*_c3.f90
           ! 2/ to suppress wall-normal velocity for Eulerian walls in mod_bc_wall_2_c3.f90
           
           ! allocations
           ! -----------
           allocate(tx_ksi_jmax(nx1_v:nx2_v,nz1_v:nz2_v),ty_ksi_jmax(nx1_v:nx2_v,nz1_v:nz2_v),tz_ksi_jmax(nx1_v:nx2_v,nz1_v:nz2_v))
           allocate(tx_phi_jmax(nx1_v:nx2_v,nz1_v:nz2_v),ty_phi_jmax(nx1_v:nx2_v,nz1_v:nz2_v),tz_phi_jmax(nx1_v:nx2_v,nz1_v:nz2_v))
           allocate(txn_ksi_jmax(nx1_v:nx2_v,nz1_v:nz2_v),tyn_ksi_jmax(nx1_v:nx2_v,nz1_v:nz2_v),tzn_ksi_jmax(nx1_v:nx2_v,nz1_v:nz2_v))
           allocate(txn_phi_jmax(nx1_v:nx2_v,nz1_v:nz2_v),tyn_phi_jmax(nx1_v:nx2_v,nz1_v:nz2_v),tzn_phi_jmax(nx1_v:nx2_v,nz1_v:nz2_v))

           ! tangent components & projection on tangent
           ! ------------------------------------------
           do k=ndz_v1,nfz_v1
              do i=ndx_v1,nfx_v1
                 gksi=sqrt(xdksi(i,ny,k)**2+ydksi(i,ny,k)**2+zdksi(i,ny,k)**2)
                 if (gksi==0.0_wp) gksi=1.0_wp
                 gphi=sqrt(xdphi(i,ny,k)**2+ydphi(i,ny,k)**2+zdphi(i,ny,k)**2)
                 if (gphi==0.0_wp) gphi=1.0_wp
                 eksi_ephi=xdksi(i,ny,k)*xdphi(i,ny,k)+ydksi(i,ny,k)*ydphi(i,ny,k)+zdksi(i,ny,k)*zdphi(i,ny,k)
                 eksi_ephi=eksi_ephi/gksi/gphi
                 fe=1.0_wp/(1.0_wp-eksi_ephi**2)                
                 tx_ksi_jmax(i,k)=fe*(xdksi(i,ny,k)/gksi-eksi_ephi*xdphi(i,ny,k)/gphi)
                 ty_ksi_jmax(i,k)=fe*(ydksi(i,ny,k)/gksi-eksi_ephi*ydphi(i,ny,k)/gphi)
                 tz_ksi_jmax(i,k)=fe*(zdksi(i,ny,k)/gksi-eksi_ephi*zdphi(i,ny,k)/gphi)
                 tx_phi_jmax(i,k)=fe*(xdphi(i,ny,k)/gphi-eksi_ephi*xdksi(i,ny,k)/gksi)
                 ty_phi_jmax(i,k)=fe*(ydphi(i,ny,k)/gphi-eksi_ephi*ydksi(i,ny,k)/gksi)
                 tz_phi_jmax(i,k)=fe*(zdphi(i,ny,k)/gphi-eksi_ephi*zdksi(i,ny,k)/gksi)
                 txn_ksi_jmax(i,k)=xdksi(i,ny,k)/gksi
                 tyn_ksi_jmax(i,k)=ydksi(i,ny,k)/gksi
                 tzn_ksi_jmax(i,k)=zdksi(i,ny,k)/gksi
                 txn_phi_jmax(i,k)=xdphi(i,ny,k)/gphi
                 tyn_phi_jmax(i,k)=ydphi(i,ny,k)/gphi
                 tzn_phi_jmax(i,k)=zdphi(i,ny,k)/gphi
              enddo
           enddo
           
           ! Oriented surface elements
           ! =========================
           ! 1/ to compute drag and lift forces
           
           ! allocate
           ! --------
           allocate(nxnds_jmax(nx,nz),nynds_jmax(nx,nz),nznds_jmax(nx,nz),ineta_jmax(nx,nz))

           ! determine correct_sign of wall normals
           ! --------------------------------------
           ! /!\ direction depends on (i,j) orientation
           !     NEEDS to set INWARDS normals for inlet

           ! vector pointing inwards (between first and second cells)
           di(1)=xc3(1,ny,1)-xc3(1,ny-1,1)
           di(2)=yc3(1,ny,1)-yc3(1,ny-1,1)
           di(3)=zc3(1,ny,1)-zc3(1,ny-1,1)
           ! scalar product with wall normal
           din=di(1)*nxn_jmax(1,1)+di(2)*nyn_jmax(1,1)+di(3)*nzn_jmax(1,1)
           ! sign
           sgn_jmax=din/abs(din)

           ! oriented surface elements
           do k=1,nz
              do i=1,nx
                 ineta_jmax(i,k)=1.0_wp/sqrt(eta_x(i,ny,k)**2+eta_y(i,ny,k)**2+eta_z(i,ny,k)**2)
                 nxnds_jmax(i,k) = eta_x(i,ny,k)*sgn_jmax
                 nynds_jmax(i,k) = eta_y(i,ny,k)*sgn_jmax
                 nznds_jmax(i,k) = eta_z(i,ny,k)*sgn_jmax
              enddo
           enddo

        endif

     endif
     
     ! Determination of normals for BC kmin (ksi=cste)
     ! ====================================
     if (BC_face(3,1)%sort.le.0) then     

        ! allocate arrays for BC normals
        ! ------------------------------
        allocate(nxn_kmin(nx,ny),nyn_kmin(nx,ny),nzn_kmin(nx,ny))
        allocate(gksigphi3_kmin(nx,ny),getagphi3_kmin(nx,ny))

        ! boundary normals along ksi
        ! --------------------------
        do j=1,ny
           do i=1,nx
              ! J²*||grad(ksi)||²
              n2_phi=phi_x(i,j,1)**2+phi_y(i,j,1)**2+phi_z(i,j,1)**2
              ! BC normal components: J*grad(phi)/J*||grad(phi)||
              nxn_kmin(i,j) = phi_x(i,j,1)/sqrt(n2_phi)
              nyn_kmin(i,j) = phi_y(i,j,1)/sqrt(n2_phi)
              nzn_kmin(i,j) = phi_z(i,j,1)/sqrt(n2_phi)
              ! compute gksigphi3 for BC: -grad(ksi).grad(phi)/||grad(phi)||²
              gksigphi3_kmin(i,j)=ksi_x(i,j,1)*phi_x(i,j,1)+ksi_y(i,j,1)*phi_y(i,j,1)+ksi_z(i,j,1)*phi_z(i,j,1)
              gksigphi3_kmin(i,j)=-gksigphi3_kmin(i,j)/n2_phi
              if (abs(gksigphi3_kmin(i,j))<1.e-16_wp) gksigphi3_kmin(i,j)=0.0_wp
              ! compute getagphi3 for BC: -grad(eta).grad(phi)/||grad(phi)||²
              getagphi3_kmin(i,j)=eta_x(i,j,1)*phi_x(i,j,1)+eta_y(i,j,1)*phi_y(i,j,1)+eta_z(i,j,1)*phi_z(i,j,1)
              getagphi3_kmin(i,j)=-getagphi3_kmin(i,j)/n2_phi
              if (abs(getagphi3_kmin(i,j))<1.e-16_wp) getagphi3_kmin(i,j)=0.0_wp
           enddo
        enddo

        if (BC_face(3,1)%sort==0) then

           ! Projection on wall tangent
           ! ==========================
           ! 1/ to enforce adiabaticity in grad_T_*_c3.f90
           ! 2/ to suppress wall-normal velocity for Eulerian walls in mod_bc_wall_2_c3.f90
           
           ! allocations
           ! -----------
           allocate(tx_ksi_kmin(nx1_v:nx2_v,ny1_v:ny2_v),ty_ksi_kmin(nx1_v:nx2_v,ny1_v:ny2_v),tz_ksi_kmin(nx1_v:nx2_v,ny1_v:ny2_v))
           allocate(tx_eta_kmin(nx1_v:nx2_v,ny1_v:ny2_v),ty_eta_kmin(nx1_v:nx2_v,ny1_v:ny2_v),tz_eta_kmin(nx1_v:nx2_v,ny1_v:ny2_v))
           allocate(txn_ksi_kmin(nx1_v:nx2_v,ny1_v:ny2_v),tyn_ksi_kmin(nx1_v:nx2_v,ny1_v:ny2_v),tzn_ksi_kmin(nx1_v:nx2_v,ny1_v:ny2_v))
           allocate(txn_eta_kmin(nx1_v:nx2_v,ny1_v:ny2_v),tyn_eta_kmin(nx1_v:nx2_v,ny1_v:ny2_v),tzn_eta_kmin(nx1_v:nx2_v,ny1_v:ny2_v))

           ! tangent components & projection on tangent
           ! ------------------------------------------
           do j=ndy_v1,nfy_v1
              do i=ndx_v1,nfx_v1
                 gksi=sqrt(xdksi(i,j,1)**2+ydksi(i,j,1)**2+zdksi(i,j,1)**2)
                 if (gksi==0.0_wp) gksi=1.0_wp
                 geta=sqrt(xdeta(i,j,1)**2+ydeta(i,j,1)**2+zdeta(i,j,1)**2)
                 if (geta==0.0_wp) geta=1.0_wp
                 eksi_eeta=xdeta(i,j,1)*xdksi(i,j,1)+ydeta(i,j,1)*ydksi(i,j,1)+zdeta(i,j,1)*zdksi(i,j,1)
                 eksi_eeta=eksi_eeta/gksi/geta
                 fe=1.0_wp/(1.0_wp-eksi_eeta**2)                
                 tx_ksi_kmin(i,j)=fe*(xdksi(i,j,1)/gksi-eksi_eeta*xdeta(i,j,1)/geta)
                 ty_ksi_kmin(i,j)=fe*(ydksi(i,j,1)/gksi-eksi_eeta*ydeta(i,j,1)/geta)
                 tz_ksi_kmin(i,j)=fe*(zdksi(i,j,1)/gksi-eksi_eeta*zdeta(i,j,1)/geta)
                 tx_eta_kmin(i,j)=fe*(xdeta(i,j,1)/geta-eksi_eeta*xdksi(i,j,1)/gksi)
                 ty_eta_kmin(i,j)=fe*(ydeta(i,j,1)/geta-eksi_eeta*ydksi(i,j,1)/gksi)
                 tz_eta_kmin(i,j)=fe*(zdeta(i,j,1)/geta-eksi_eeta*zdksi(i,j,1)/gksi)
                 txn_ksi_kmin(i,j)=xdksi(i,j,1)/gksi
                 tyn_ksi_kmin(i,j)=ydksi(i,j,1)/gksi
                 tzn_ksi_kmin(i,j)=zdksi(i,j,1)/gksi
                 txn_eta_kmin(i,j)=xdeta(i,j,1)/geta
                 tyn_eta_kmin(i,j)=ydeta(i,j,1)/geta
                 tzn_eta_kmin(i,j)=zdeta(i,j,1)/geta
              enddo
           enddo
           
           ! Oriented surface elements
           ! =========================
           ! 1/ to compute drag and lift forces
           
           ! allocate
           ! --------
           allocate(nxnds_kmin(nx,ny),nynds_kmin(nx,ny),nznds_kmin(nx,ny),inphi_kmin(nx,ny))

           ! determine correct_sign of wall normals
           ! --------------------------------------
           ! /!\ direction depends on (i,j) orientation
           !     NEEDS to set INWARDS normals for inlet

           ! vector pointing inwards (between first and second cells)
           di(1)=xc3(1,1,2)-xc3(1,1,1)
           di(2)=yc3(1,1,2)-yc3(1,1,1)
           di(3)=zc3(1,1,2)-zc3(1,1,1)
           ! scalar product with wall normal
           din=di(1)*nxn_kmin(1,1)+di(2)*nyn_kmin(1,1)+di(3)*nzn_kmin(1,1)
           ! sign
           sgn_kmin=-din/abs(din)

           ! oriented surface elements
           do j=1,ny
              do i=1,nx
                 inphi_kmin(i,j)=1.0_wp/sqrt(phi_x(i,j,1)**2+phi_y(i,j,1)**2+phi_z(i,j,1)**2)
                 nxnds_kmin(i,j) = phi_x(i,j,1)*sgn_kmin
                 nynds_kmin(i,j) = phi_y(i,j,1)*sgn_kmin
                 nznds_kmin(i,j) = phi_z(i,j,1)*sgn_kmin
              enddo
           enddo

        endif

     endif

     ! Determination of normals for BC kmax (ksi=cste)
     ! ====================================
     if (BC_face(3,2)%sort.le.0) then     

        ! allocate arrays for BC normals
        ! ------------------------------
        allocate(nxn_kmax(nx,ny),nyn_kmax(nx,ny),nzn_kmax(nx,ny))
        allocate(gksigphi3_kmax(nx,ny),getagphi3_kmax(nx,ny))

        ! boundary normals along ksi
        ! --------------------------
        do j=1,ny
           do i=1,nx
              ! J²*||grad(ksi)||²
              n2_phi=phi_x(i,j,nz)**2+phi_y(i,j,nz)**2+phi_z(i,j,nz)**2
              ! BC normal components: J*grad(phi)/J*||grad(phi)||
              nxn_kmax(i,j) = phi_x(i,j,nz)/sqrt(n2_phi)
              nyn_kmax(i,j) = phi_y(i,j,nz)/sqrt(n2_phi)
              nzn_kmax(i,j) = phi_z(i,j,nz)/sqrt(n2_phi)
              ! compute gksigphi3 for BC: -grad(ksi).grad(phi)/||grad(phi)||²
              gksigphi3_kmax(i,j)=ksi_x(i,j,nz)*phi_x(i,j,nz)+ksi_y(i,j,nz)*phi_y(i,j,nz)+ksi_z(i,j,nz)*phi_z(i,j,nz)
              gksigphi3_kmax(i,j)=-gksigphi3_kmax(i,j)/n2_phi
              if (abs(gksigphi3_kmax(i,j))<1.e-16_wp) gksigphi3_kmax(i,j)=0.0_wp
              ! compute getagphi3 for BC: -grad(eta).grad(phi)/||grad(phi)||²
              getagphi3_kmax(i,j)=eta_x(i,j,nz)*phi_x(i,j,nz)+eta_y(i,j,nz)*phi_y(i,j,nz)+eta_z(i,j,nz)*phi_z(i,j,nz)
              getagphi3_kmax(i,j)=-getagphi3_kmax(i,j)/n2_phi
              if (abs(getagphi3_kmax(i,j))<1.e-16_wp) getagphi3_kmax(i,j)=0.0_wp
           enddo
        enddo

        if (BC_face(3,2)%sort==0) then

           ! Projection on wall tangent
           ! ==========================
           ! 1/ to enforce adiabaticity in grad_T_*_c3.f90
           ! 2/ to suppress wall-normal velocity for Eulerian walls in mod_bc_wall_2_c3.f90
           
           ! allocations
           ! -----------
           allocate(tx_ksi_kmax(nx1_v:nx2_v,ny1_v:ny2_v),ty_ksi_kmax(nx1_v:nx2_v,ny1_v:ny2_v),tz_ksi_kmax(nx1_v:nx2_v,ny1_v:ny2_v))
           allocate(tx_eta_kmax(nx1_v:nx2_v,ny1_v:ny2_v),ty_eta_kmax(nx1_v:nx2_v,ny1_v:ny2_v),tz_eta_kmax(nx1_v:nx2_v,ny1_v:ny2_v))
           allocate(txn_ksi_kmax(nx1_v:nx2_v,ny1_v:ny2_v),tyn_ksi_kmax(nx1_v:nx2_v,ny1_v:ny2_v),tzn_ksi_kmax(nx1_v:nx2_v,ny1_v:ny2_v))
           allocate(txn_eta_kmax(nx1_v:nx2_v,ny1_v:ny2_v),tyn_eta_kmax(nx1_v:nx2_v,ny1_v:ny2_v),tzn_eta_kmax(nx1_v:nx2_v,ny1_v:ny2_v))

           ! tangent components & projection on tangent
           ! ------------------------------------------
           do j=ndy_v1,nfy_v1
              do i=ndx_v1,nfx_v1
                 gksi=sqrt(xdksi(i,j,nz)**2+ydksi(i,j,nz)**2+zdksi(i,j,nz)**2)
                 if (gksi==0.0_wp) gksi=1.0_wp
                 geta=sqrt(xdeta(i,j,nz)**2+ydeta(i,j,nz)**2+zdeta(i,j,nz)**2)
                 if (geta==0.0_wp) geta=1.0_wp
                 eksi_eeta=xdeta(i,j,nz)*xdksi(i,j,nz)+ydeta(i,j,nz)*ydksi(i,j,nz)+zdeta(i,j,nz)*zdksi(i,j,nz)
                 eksi_eeta=eksi_eeta/gksi/geta
                 fe=1.0_wp/(1.0_wp-eksi_eeta**2)                
                 tx_ksi_kmax(i,j)=fe*(xdksi(i,j,nz)/gksi-eksi_eeta*xdeta(i,j,nz)/geta)
                 ty_ksi_kmax(i,j)=fe*(ydksi(i,j,nz)/gksi-eksi_eeta*ydeta(i,j,nz)/geta)
                 tz_ksi_kmax(i,j)=fe*(zdksi(i,j,nz)/gksi-eksi_eeta*zdeta(i,j,nz)/geta)
                 tx_eta_kmax(i,j)=fe*(xdeta(i,j,nz)/geta-eksi_eeta*xdksi(i,j,nz)/gksi)
                 ty_eta_kmax(i,j)=fe*(ydeta(i,j,nz)/geta-eksi_eeta*ydksi(i,j,nz)/gksi)
                 tz_eta_kmax(i,j)=fe*(zdeta(i,j,nz)/geta-eksi_eeta*zdksi(i,j,nz)/gksi)
                 txn_ksi_kmax(i,j)=xdksi(i,j,nz)/gksi
                 tyn_ksi_kmax(i,j)=ydksi(i,j,nz)/gksi
                 tzn_ksi_kmax(i,j)=zdksi(i,j,nz)/gksi
                 txn_eta_kmax(i,j)=xdeta(i,j,nz)/geta
                 tyn_eta_kmax(i,j)=ydeta(i,j,nz)/geta
                 tzn_eta_kmax(i,j)=zdeta(i,j,nz)/geta
              enddo
           enddo
           
           ! Oriented surface elements
           ! =========================
           ! 1/ to compute drag and lift forces
           
           ! allocate
           ! --------
           allocate(nxnds_kmax(nx,ny),nynds_kmax(nx,ny),nznds_kmax(nx,ny),inphi_kmax(nx,ny))

           ! determine correct_sign of wall normals
           ! --------------------------------------
           ! /!\ direction depends on (i,j) orientation
           !     NEEDS to set INWARDS normals for inlet

           ! vector pointing inwards (between first and second cells)
           di(1)=xc3(1,1,nz)-xc3(1,1,nz-1)
           di(2)=yc3(1,1,nz)-yc3(1,1,nz-1)
           di(3)=zc3(1,1,nz)-zc3(1,1,nz-1)
           ! scalar product with wall normal
           din=di(1)*nxn_kmax(1,1)+di(2)*nyn_kmax(1,1)+di(3)*nzn_kmax(1,1)
           ! sign
           sgn_kmax=din/abs(din)

           ! oriented surface elements
           do j=1,ny
              do i=1,nx
                 inphi_kmax(i,j)=1.0_wp/sqrt(phi_x(i,j,nz)**2+phi_y(i,j,nz)**2+phi_z(i,j,nz)**2)
                 nxnds_kmax(i,j) = phi_x(i,j,nz)*sgn_kmax
                 nynds_kmax(i,j) = phi_y(i,j,nz)*sgn_kmax
                 nznds_kmax(i,j) = phi_z(i,j,nz)*sgn_kmax
              enddo
           enddo

        endif

     endif

     ! Free intermediate arrays
     ! ------------------------
     deallocate(xdksi,ydksi,zdksi)
     deallocate(xdeta,ydeta,zdeta)
     deallocate(xdphi,ydphi,zdphi)

       ! **************
  else ! 2D curvilinear
       ! **************

     ! Determination of normals for BC imin (ksi=cste)
     ! ====================================
     if (BC_face(1,1)%sort.le.0) then     

        ! allocate arrays for BC normals
        ! ------------------------------
        allocate(nxn_imin(ny,nz),nyn_imin(ny,nz),nzn_imin(ny,nz))
        allocate(gksigeta_imin(ny))

        ! boundary normals along ksi
        ! --------------------------
        do k=1,nz
           do j=1,ny
              ! J²*||grad(ksi)||²
              n2_ksi=y_eta(1,j)**2+x_eta(1,j)**2
              ! BC normal components: J*grad(ksi)/J*||grad(ksi)||
              nxn_imin(j,k) =  y_eta(1,j)/sqrt(n2_ksi)
              nyn_imin(j,k) = -x_eta(1,j)/sqrt(n2_ksi)
              ! compute gksigeta for BC: -grad(eta).grad(ksi)/||grad(ksi)||²
              gksigeta_imin(j)=x_ksi(1,j)*x_eta(1,j)+y_ksi(1,j)*y_eta(1,j)
              gksigeta_imin(j)=-gksigeta_imin(j)/n2_ksi
              if (abs(gksigeta_imin(j))<1.e-16_wp) gksigeta_imin(j)=0.0_wp
           enddo
        enddo
        nzn_imin=0.0_wp

        if (BC_face(1,1)%sort==0) then

           ! compute components of n.dl along wall BC along wall BC
           ! ------------------------------------------------------
           allocate(nxndl_imin(1:ny),nyndl_imin(1:ny))

           ! determine correct_sign of wall normals
           ! --------------------------------------
           ! /!\ direction depends on (i,j) orientation
           !     NEEDS to set INWARDS normals for inlet

           ! vector pointing inwards (between first and second cells)
           di(1)=xc(2,1)-xc(1,1)
           di(2)=yc(2,1)-yc(1,1)
           ! scalar product with wall normal
           din=di(1)*nxn_imin(1,1)+di(2)*nyn_imin(1,1)
           ! sign
           sgn_imin=-din/abs(din)

           nxndl_imin=0.0_wp
           nyndl_imin=0.0_wp
           do j=1,ny
              nxndl_imin(j) =  y_eta(1,j)*sgn_imin
              nyndl_imin(j) = -x_eta(1,j)*sgn_imin
           enddo
           
           ! projection on wall tangent (to enforce adiabaticity)
           ! ----------------------------------------------------
           allocate(txeta_imin(ny1_v:ny2_v),tyeta_imin(ny1_v:ny2_v))

           do j=ndy_v1,nfy_v1
              txeta_imin(j)=x_eta_v(1,j)/(y_eta_v(1,j)**2+x_eta_v(1,j)**2)
              tyeta_imin(j)=y_eta_v(1,j)/(y_eta_v(1,j)**2+x_eta_v(1,j)**2)
           enddo

        endif

     endif

     ! Determination of normals for BC imax (ksi=cste)
     ! ====================================
     if (BC_face(1,2)%sort.le.0) then     

        ! allocate arrays for BC normals
        ! ------------------------------
        allocate(nxn_imax(ny,nz),nyn_imax(ny,nz),nzn_imax(ny,nz))
        allocate(gksigeta_imax(ny))

        ! boundary along ksi
        ! ------------------
        do k=1,nz
           do j=1,ny
              ! J²*||grad(ksi)||²
              n2_ksi=y_eta(nx,j)**2+x_eta(nx,j)**2
              ! BC normal components: J*grad(ksi)/J*||grad(ksi)||
              nxn_imax(j,k) =  y_eta(nx,j)/sqrt(n2_ksi)
              nyn_imax(j,k) = -x_eta(nx,j)/sqrt(n2_ksi)
              ! compute gksigeta for BC: -grad(eta).grad(ksi)/||grad(ksi)||²
              gksigeta_imax(j)=x_ksi(nx,j)*x_eta(nx,j)+y_ksi(nx,j)*y_eta(nx,j)
              gksigeta_imax(j)=-gksigeta_imax(j)/n2_ksi
              if (abs(gksigeta_imax(j))<1.e-16_wp) gksigeta_imax(j)=0.0_wp
           enddo
        enddo
        nzn_imax=0.0_wp

        if (BC_face(1,2)%sort==0) then     

           ! compute components of n.dl along wall BC
           ! ----------------------------------------
           allocate(nxndl_imax(ny),nyndl_imax(ny))

           ! determine correct_sign of wall normals
           ! --------------------------------------
           ! /!\ direction depends on (i,j) orientation
           !     NEEDS to set INWARDS normals for inlet

           ! vector pointing inwards (between first and second cells)
           di(1)=xc(nx,1)-xc(nx-1,1)
           di(2)=yc(nx,1)-yc(nx-1,1)
           ! scalar product with wall normal
           din=di(1)*nxn_imax(1,1)+di(2)*nyn_imax(1,1)
           ! sign
           sgn_imax=din/abs(din)

           nxndl_imax=0.0_wp
           nyndl_imax=0.0_wp
           do j=1,ny
              nxndl_imax(j) =  y_eta(nx,j)*sgn_imax
              nyndl_imax(j) = -x_eta(nx,j)*sgn_imax
           enddo
           
           ! projection on wall tangent (to enforce adiabaticity)
           ! ----------------------------------------------------
           allocate(txeta_imax(ny1_v:ny2_v),tyeta_imax(ny1_v:ny2_v))

           do j=ndy_v1,nfy_v1
              txeta_imax(j)=x_eta_v(nx,j)/(y_eta_v(nx,j)**2+x_eta_v(nx,j)**2)
              tyeta_imax(j)=y_eta_v(nx,j)/(y_eta_v(nx,j)**2+x_eta_v(nx,j)**2)
           enddo

        endif

     endif

     ! Determination of normals for BC jmin (eta=cste)
     ! ====================================
     if (BC_face(2,1)%sort.le.0) then     

        ! allocate arrays for BC normals
        ! ------------------------------
        allocate(nxn_jmin(nx,nz),nyn_jmin(nx,nz),nzn_jmin(nx,nz))
        allocate(gksigeta_jmin(nx))

        ! boundary along eta
        ! ------------------
        do k=1,nz
           do i=1,nx
              ! J²*||grad(eta)||²
              n2_eta=y_ksi(i,1)**2+x_ksi(i,1)**2
              ! BC normal components: J*grad(eta)/J*||grad(eta)||
              nxn_jmin(i,k) = -y_ksi(i,1)/sqrt(n2_eta)
              nyn_jmin(i,k) =  x_ksi(i,1)/sqrt(n2_eta)
              ! compute gksigeta for BC: -grad(ksi).grad(eta)/||grad(eta)||²
              gksigeta_jmin(i)=x_ksi(i,1)*x_eta(i,1)+y_ksi(i,1)*y_eta(i,1)
              gksigeta_jmin(i)=-gksigeta_jmin(i)/n2_eta
              if (abs(gksigeta_jmin(i))<1.e-16_wp) gksigeta_jmin(i)=0.0_wp
           enddo
        enddo
        nzn_jmin=0.0_wp

        if (BC_face(2,1)%sort==0) then

           ! compute components of n.dl along wall BC
           ! ----------------------------------------
           allocate(nxndl_jmin(nx),nyndl_jmin(nx))

           ! determine correct_sign of wall normals
           ! --------------------------------------
           ! /!\ direction depends on (i,j) orientation
           !     NEEDS to set INWARDS normals for inlet

           ! vector pointing inwards (between first and second cells)
           di(1)=xc(1,2)-xc(1,1)
           di(2)=yc(1,2)-yc(1,1)
           ! scalar product with wall normal
           din=di(1)*nxn_jmin(1,1)+di(2)*nyn_jmin(1,1)
           ! sign
           sgn_jmin=-din/abs(din)

           nxndl_jmin=0.0_wp
           nyndl_jmin=0.0_wp
           do i=1,nx
              nxndl_jmin(i) = -y_ksi(i,1)*sgn_jmin
              nyndl_jmin(i) =  x_ksi(i,1)*sgn_jmin
           enddo

           ! projection on wall tangent (to enforce adiabaticity)
           ! ----------------------------------------------------
           allocate(txksi_jmin(nx1_v:nx2_v),tyksi_jmin(nx1_v:nx2_v))

           do i=ndx_v1,nfx_v1
              txksi_jmin(i)=x_ksi_v(i,1)/(y_ksi_v(i,1)**2+x_ksi_v(i,1)**2)
              tyksi_jmin(i)=y_ksi_v(i,1)/(y_ksi_v(i,1)**2+x_ksi_v(i,1)**2)
           enddo

        endif

     endif

     ! Determination of normals for BC jmax (eta=cste)
     ! ====================================
     if (BC_face(2,2)%sort.le.0) then     

        ! allocate arrays for BC normals
        ! ------------------------------
        allocate(nxn_jmax(nx,nz),nyn_jmax(nx,nz),nzn_jmax(nx,nz))
        allocate(gksigeta_jmax(nx))

        ! boundary along eta
        ! ------------------
        do k=1,nz
           do i=1,nx
              ! J²*||grad(eta)||²
              n2_eta=y_ksi(i,ny)**2+x_ksi(i,ny)**2
              ! BC normal components: J*grad(eta)/J*||grad(eta)||
              nxn_jmax(i,k) = -y_ksi(i,ny)/sqrt(n2_eta)
              nyn_jmax(i,k) =  x_ksi(i,ny)/sqrt(n2_eta)
              ! compute gksigeta for BC: -grad(ksi).grad(eta)/||grad(eta)||²
              gksigeta_jmax(i)=x_ksi(i,ny)*x_eta(i,ny)+y_ksi(i,ny)*y_eta(i,ny)
              gksigeta_jmax(i)=-gksigeta_jmax(i)/n2_eta
              if (abs(gksigeta_jmax(i))<1.e-16_wp) gksigeta_jmax(i)=0.0_wp
           enddo
        enddo
        nzn_jmax=0.0_wp

        if (BC_face(2,2)%sort==0) then

           ! compute components of n.dl along wall BC
           ! ----------------------------------------
           allocate(nxndl_jmax(nx),nyndl_jmax(nx))

           ! determine correct_sign of wall normals
           ! --------------------------------------
           ! /!\ direction depends on (i,j) orientation
           !     NEEDS to set INWARDS normals for inlet

           ! vector pointing inwards (between first and second cells)
           di(1)=xc(1,ny)-xc(1,ny-1)
           di(2)=yc(1,ny)-yc(1,ny-1)
           ! scalar product with wall normal
           din=di(1)*nxn_jmax(1,1)+di(2)*nyn_jmax(1,1)
           ! sign
           sgn_jmax=din/abs(din)

           nxndl_jmax=0.0_wp
           nyndl_jmax=0.0_wp
           do i=1,nx
              nxndl_jmax(i) = -y_ksi(i,ny)*sgn_jmax
              nyndl_jmax(i) =  x_ksi(i,ny)*sgn_jmax
           enddo
           
           ! projection on wall tangent (to enforce adiabaticity)
           ! ----------------------------------------------------
           allocate(txksi_jmax(nx1_v:nx2_v),tyksi_jmax(nx1_v:nx2_v))

           do i=ndx_v1,nfx_v1
              txksi_jmax(i)=x_ksi_v(i,ny)/(y_ksi_v(i,ny)**2+x_ksi_v(i,ny)**2)
              tyksi_jmax(i)=y_ksi_v(i,ny)/(y_ksi_v(i,ny)**2+x_ksi_v(i,ny)**2)
           enddo
           
        endif

     endif

  endif

end subroutine grid_normals

