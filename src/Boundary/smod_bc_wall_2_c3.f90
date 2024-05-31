!================================================================================
submodule (mod_bc_wall) smod_bc_wall_2_c3
!================================================================================
  !> Submodule to apply wall Boundary Conditions on conservative variables
  !> - 3D curvilinear version - (advancement of rho on boundary)
!================================================================================

contains

  !==============================================================================
  module subroutine bc_wall_imin_slip_c3
  !==============================================================================
    !> Apply wall Boundary Conditions at imin: slip wall
    !> - 3D curvilinear version - (advancement of rho on boundary)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: j,k
    real(wp) :: rhovtan_eta,rhovtan_phi
    ! ---------------------------------------------------------------------------

    ! Impose wall conditions on conservative variables
    ! ================================================
    do k=1,nz
       do j=1,ny
          ! define tangential velocity components
          rhovtan_eta= rhou_n(1,j,k)*tx_eta_imin(j,k) &
                     + rhov_n(1,j,k)*ty_eta_imin(j,k) &
                     + rhow_n(1,j,k)*tz_eta_imin(j,k)
          rhovtan_phi= rhou_n(1,j,k)*tx_phi_imin(j,k) &
                     + rhov_n(1,j,k)*ty_phi_imin(j,k) &
                     + rhow_n(1,j,k)*tz_phi_imin(j,k)
          
          ! slip condition: nullity of normal velocity component
          rhou_n(1,j,k)=rhovtan_eta*txn_eta_imin(j,k)+rhovtan_phi*txn_phi_imin(j,k)
          rhov_n(1,j,k)=rhovtan_eta*tyn_eta_imin(j,k)+rhovtan_phi*tyn_phi_imin(j,k)
          rhow_n(1,j,k)=rhovtan_eta*tzn_eta_imin(j,k)+rhovtan_phi*tzn_phi_imin(j,k)
       enddo
    enddo

  end subroutine bc_wall_imin_slip_c3

  !==============================================================================
  module subroutine bc_wall_imin_adiabatic_c3
  !==============================================================================
    !> Apply wall Boundary Conditions at imin: no slip + adiabatic
    !> - 3D curvilinear version - (advancement of rho on boundary)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: j,k
    ! ---------------------------------------------------------------------------

    ! Impose wall conditions on conservative variables
    ! ================================================
    do k=1,nz
       do j=1,ny
          ! no slip conditions
          rhou_n(1,j,k)=0.0_wp
          rhov_n(1,j,k)=0.0_wp
          rhow_n(1,j,k)=0.0_wp
       enddo
    enddo
    
    do k=ndz_imin,nfz_imin
       do j=ndy_imin,nfy_imin
          ! dT/dn=0
          T_wall= onethird*( 4.0_wp*Tmp(2,j,k)-Tmp(3,j,k)         &
                + (Tmp(1,j+1,k)-Tmp(1,j-1,k))*getagksi3_imin(j,k) &
                + (Tmp(1,j,k+1)-Tmp(1,j,k-1))*gphigksi3_imin(j,k) )
          
          rhoe_n(1,j,k)=rho_n(1,j,k)*ecalc_tro(T_wall,rho_n(1,j,k))
       enddo
    enddo

    ! edge with other BC at jmin
    ! --------------------------
    if (ndy_imin==2) then
       j=1
       do k=ndz_imin,nfz_imin
          ! dT/dn=0
          T_wall= onethird*( 4.0_wp*Tmp(2,j,k)-Tmp(3,j,k)        &
                +2.0_wp*(Tmp(1,j+1,k)-Tmp(1,j,k))*getagksi3_imin(j,k) &
                +(Tmp(1,j,k+1)-Tmp(1,j,k-1))*gphigksi3_imin(j,k) )

          rhoe_n(1,j,k)=rho_n(1,j,k)*ecalc_tro(T_wall,rho_n(1,j,k))
       enddo
    endif

    ! edge with other BC at jmax
    ! --------------------------
    if (ndy_imin==ny-1) then
       j=ny
       do k=ndz_imin,nfz_imin
          ! dT/dn=0
          T_wall= onethird*( 4.0_wp*Tmp(2,j,k)-Tmp(3,j,k)        &
                +2.0_wp*(Tmp(1,j,k)-Tmp(1,j-1,k))*getagksi3_imin(j,k) &
                +(Tmp(1,j,k+1)-Tmp(1,j,k-1))*gphigksi3_imin(j,k) )

          rhoe_n(1,j,k)=rho_n(1,j,k)*ecalc_tro(T_wall,rho_n(1,j,k))
       enddo
    endif

    ! edge with other BC at kmin
    ! --------------------------
    if (ndz_imin==2) then
       k=1
       do j=ndy_imin,nfy_imin
          ! dT/dn=0
          T_wall= onethird*( 4.0_wp*Tmp(2,j,k)-Tmp(3,j,k)        &
                +(Tmp(1,j+1,k)-Tmp(1,j-1,k))*getagksi3_imin(j,k) &
                +2.0_wp*(Tmp(1,j,k+1)-Tmp(1,j,k))*gphigksi3_imin(j,k) )

          rhoe_n(1,j,k)=rho_n(1,j,k)*ecalc_tro(T_wall,rho_n(1,j,k))
       enddo
    endif

    ! edge with other BC at kmax
    ! --------------------------
    if (ndz_imin==nz-1) then
       k=nz
       do j=ndy_imin,nfy_imin
          ! dT/dn=0
          T_wall= onethird*( 4.0_wp*Tmp(2,j,k)-Tmp(3,j,k)        &
                +(Tmp(1,j+1,k)-Tmp(1,j-1,k))*getagksi3_imin(j,k) &
                +2.0_wp*(Tmp(1,j,k)-Tmp(1,j,k-1))*gphigksi3_imin(j,k) )

          rhoe_n(1,j,k)=rho_n(1,j,k)*ecalc_tro(T_wall,rho_n(1,j,k))
       enddo
    endif

    ! corner with other BC at jmin-kmin
    ! ---------------------------------
    if ((ndy_imin==2).and.(ndz_imin==2)) then
       j=1
       k=1
       ! dT/dn=0
       T_wall= onethird*( 4.0_wp*Tmp(2,j,k)-Tmp(3,j,k)        &
             +2.0_wp*(Tmp(1,j+1,k)-Tmp(1,j,k))*getagksi3_imin(j,k) &
             +2.0_wp*(Tmp(1,j,k+1)-Tmp(1,j,k))*gphigksi3_imin(j,k) )

       rhoe_n(1,j,k)=rho_n(1,j,k)*ecalc_tro(T_wall,rho_n(1,j,k))
    endif
    
    ! corner with other BC at jmin-kmax
    ! ---------------------------------
    if ((ndy_imin==2).and.(ndz_imin==nz-1)) then
       j=1
       k=nz
       ! dT/dn=0
       T_wall= onethird*( 4.0_wp*Tmp(2,j,k)-Tmp(3,j,k)        &
             +2.0_wp*(Tmp(1,j+1,k)-Tmp(1,j,k))*getagksi3_imin(j,k) &
             +2.0_wp*(Tmp(1,j,k)-Tmp(1,j,k-1))*gphigksi3_imin(j,k) )

       rhoe_n(1,j,k)=rho_n(1,j,k)*ecalc_tro(T_wall,rho_n(1,j,k))
    endif
    
    ! corner with other BC at jmax-kmin
    ! ---------------------------------
    if ((ndy_imin==ny-1).and.(ndz_imin==2)) then
       j=ny
       k=1
       ! dT/dn=0
       T_wall= onethird*( 4.0_wp*Tmp(2,j,k)-Tmp(3,j,k)        &
             +2.0_wp*(Tmp(1,j,k)-Tmp(1,j-1,k))*getagksi3_imin(j,k) &
             +2.0_wp*(Tmp(1,j,k+1)-Tmp(1,j,k))*gphigksi3_imin(j,k) )

       rhoe_n(1,j,k)=rho_n(1,j,k)*ecalc_tro(T_wall,rho_n(1,j,k))
    endif

    ! corner with other BC at jmax-kmax
    ! ---------------------------------
    if ((ndy_imin==ny-1).and.(ndz_imin==nz-1)) then
       j=ny
       k=nz
       ! dT/dn=0
       T_wall= onethird*( 4.0_wp*Tmp(2,j,k)-Tmp(3,j,k)        &
             +2.0_wp*(Tmp(1,j,k)-Tmp(1,j-1,k))*getagksi3_imin(j,k) &
             +2.0_wp*(Tmp(1,j,k)-Tmp(1,j,k-1))*gphigksi3_imin(j,k) )

       rhoe_n(1,j,k)=rho_n(1,j,k)*ecalc_tro(T_wall,rho_n(1,j,k))
    endif
    
  end subroutine bc_wall_imin_adiabatic_c3

  !==============================================================================
  module subroutine bc_wall_imax_slip_c3
  !==============================================================================
    !> Apply wall Boundary Conditions at imax: slip wall
    !> - 3D curvilinear version - (advancement of rho on boundary)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: j,k
    real(wp) :: rhovtan_eta,rhovtan_phi
    ! ---------------------------------------------------------------------------

    ! Impose wall conditions on conservative variables
    ! ================================================
    do k=1,nz
       do j=1,ny
          ! define tangential velocity components
          rhovtan_eta= rhou_n(nx,j,k)*tx_eta_imax(j,k) &
                     + rhov_n(nx,j,k)*ty_eta_imax(j,k) &
                     + rhow_n(nx,j,k)*tz_eta_imax(j,k)
          rhovtan_phi= rhou_n(nx,j,k)*tx_phi_imax(j,k) &
                     + rhov_n(nx,j,k)*ty_phi_imax(j,k) &
                     + rhow_n(nx,j,k)*tz_phi_imax(j,k)
          
          ! slip condition: nullity of normal velocity component
          rhou_n(nx,j,k)=rhovtan_eta*txn_eta_imax(j,k)+rhovtan_phi*txn_phi_imax(j,k)
          rhov_n(nx,j,k)=rhovtan_eta*tyn_eta_imax(j,k)+rhovtan_phi*tyn_phi_imax(j,k)
          rhow_n(nx,j,k)=rhovtan_eta*tzn_eta_imax(j,k)+rhovtan_phi*tzn_phi_imax(j,k)
       enddo
    enddo

  end subroutine bc_wall_imax_slip_c3

  !==============================================================================
  module subroutine bc_wall_imax_adiabatic_c3
  !==============================================================================
    !> Apply wall Boundary Conditions at imax: no slip + adiabatic
    !> - 3D curvilinear version - (advancement of rho on boundary)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: j,k
    ! ---------------------------------------------------------------------------

    ! Impose wall conditions on conservative variables
    ! ================================================
    do k=1,nz
       do j=1,ny
          ! no slip conditions
          rhou_n(nx,j,k)=0.0_wp
          rhov_n(nx,j,k)=0.0_wp
          rhow_n(nx,j,k)=0.0_wp
       enddo
    enddo

    do k=ndz_imax,nfz_imax
       do j=ndy_imax,nfy_imax
          ! dT/dn=0
          T_wall= onethird*( 4.0_wp*Tmp(nx-1,j,k)-Tmp(nx-2,j,k)     &
                + (Tmp(nx,j+1,k)-Tmp(nx,j-1,k))*getagksi3_imax(j,k) &
                + (Tmp(nx,j,k+1)-Tmp(nx,j,k-1))*gphigksi3_imax(j,k) )
          
          rhoe_n(nx,j,k)=rho_n(nx,j,k)*ecalc_tro(T_wall,rho_n(nx,j,k))
       enddo
    enddo

    ! edge with other BC at jmin
    ! --------------------------
    if (ndy_imax==2) then
       j=1
       do k=ndz_imax,nfz_imax
          ! dT/dn=0
          T_wall= onethird*( 4.0_wp*Tmp(nx-1,j,k)-Tmp(nx-2,j,k)    &
                +2.0_wp*(Tmp(nx,j+1,k)-Tmp(nx,j,k))*getagksi3_imax(j,k) &
                +(Tmp(nx,j,k+1)-Tmp(nx,j,k-1))*gphigksi3_imax(j,k) )

          rhoe_n(nx,j,k)=rho_n(nx,j,k)*ecalc_tro(T_wall,rho_n(nx,j,k))
       enddo
    endif

    ! edge with other BC at jmax
    ! --------------------------
    if (ndy_imax==ny-1) then
       j=ny
       do k=ndz_imax,nfz_imax
          ! dT/dn=0
          T_wall= onethird*( 4.0_wp*Tmp(nx-1,j,k)-Tmp(nx-2,j,k)    &
                +2.0_wp*(Tmp(nx,j,k)-Tmp(nx,j-1,k))*getagksi3_imax(j,k) &
                +(Tmp(nx,j,k+1)-Tmp(nx,j,k-1))*gphigksi3_imax(j,k) )

          rhoe_n(nx,j,k)=rho_n(nx,j,k)*ecalc_tro(T_wall,rho_n(nx,j,k))
       enddo
    endif

    ! edge with other BC at kmin
    ! --------------------------
    if (ndz_imax==2) then
       k=1
       do j=ndy_imax,nfy_imax
          ! dT/dn=0
          T_wall= onethird*( 4.0_wp*Tmp(nx-1,j,k)-Tmp(nx-2,j,k)    &
                +(Tmp(nx,j+1,k)-Tmp(nx,j-1,k))*getagksi3_imax(j,k) &
                +2.0_wp*(Tmp(nx,j,k+1)-Tmp(nx,j,k))*gphigksi3_imax(j,k) )

          rhoe_n(nx,j,k)=rho_n(nx,j,k)*ecalc_tro(T_wall,rho_n(nx,j,k))
       enddo
    endif

    ! edge with other BC at kmax
    ! --------------------------
    if (ndz_imax==nz-1) then
       k=nz
       do j=ndy_imax,nfy_imax
          ! dT/dn=0
          T_wall= onethird*( 4.0_wp*Tmp(nx-1,j,k)-Tmp(nx-2,j,k)    &
                +(Tmp(nx,j+1,k)-Tmp(nx,j-1,k))*getagksi3_imax(j,k) &
                +2.0_wp*(Tmp(nx,j,k)-Tmp(nx,j,k-1))*gphigksi3_imax(j,k) )

          rhoe_n(nx,j,k)=rho_n(nx,j,k)*ecalc_tro(T_wall,rho_n(nx,j,k))
       enddo
    endif

    ! corner with other BC at jmin-kmin
    ! ---------------------------------
    if ((ndy_imax==2).and.(ndz_imax==2)) then
       j=1
       k=1
       ! dT/dn=0
       T_wall= onethird*( 4.0_wp*Tmp(nx-1,j,k)-Tmp(nx-2,j,k)    &
             +2.0_wp*(Tmp(nx,j+1,k)-Tmp(nx,j,k))*getagksi3_imax(j,k) &
             +2.0_wp*(Tmp(nx,j,k+1)-Tmp(nx,j,k))*gphigksi3_imax(j,k) )

       rhoe_n(nx,j,k)=rho_n(nx,j,k)*ecalc_tro(T_wall,rho_n(nx,j,k))
    endif
    
    ! corner with other BC at jmin-kmax
    ! ---------------------------------
    if ((ndy_imax==2).and.(ndz_imax==nz-1)) then
       j=1
       k=nz
       ! dT/dn=0
       T_wall= onethird*( 4.0_wp*Tmp(nx-1,j,k)-Tmp(nx-2,j,k)    &
             +2.0_wp*(Tmp(nx,j+1,k)-Tmp(nx,j,k))*getagksi3_imax(j,k) &
             +2.0_wp*(Tmp(nx,j,k)-Tmp(nx,j,k-1))*gphigksi3_imax(j,k) )

       rhoe_n(nx,j,k)=rho_n(nx,j,k)*ecalc_tro(T_wall,rho_n(nx,j,k))
    endif
    
    ! corner with other BC at jmax-kmin
    ! ---------------------------------
    if ((ndy_imax==ny-1).and.(ndz_imax==2)) then
       j=ny
       k=1
       ! dT/dn=0
       T_wall= onethird*( 4.0_wp*Tmp(nx-1,j,k)-Tmp(nx-2,j,k)    &
             +2.0_wp*(Tmp(nx,j,k)-Tmp(nx,j-1,k))*getagksi3_imax(j,k) &
             +2.0_wp*(Tmp(nx,j,k+1)-Tmp(nx,j,k))*gphigksi3_imax(j,k) )

       rhoe_n(nx,j,k)=rho_n(nx,j,k)*ecalc_tro(T_wall,rho_n(nx,j,k))
    endif

    ! corner with other BC at jmax-kmax
    ! ---------------------------------
    if ((ndy_imax==ny-1).and.(ndz_imax==nz-1)) then
       j=ny
       k=nz
       ! dT/dn=0
       T_wall= onethird*( 4.0_wp*Tmp(nx-1,j,k)-Tmp(nx-2,j,k)    &
             +2.0_wp*(Tmp(nx,j,k)-Tmp(nx,j-1,k))*getagksi3_imax(j,k) &
             +2.0_wp*(Tmp(nx,j,k)-Tmp(nx,j,k-1))*gphigksi3_imax(j,k) )

       rhoe_n(nx,j,k)=rho_n(nx,j,k)*ecalc_tro(T_wall,rho_n(nx,j,k))
    endif
    
  end subroutine bc_wall_imax_adiabatic_c3

  !==============================================================================
  module subroutine bc_wall_jmin_slip_c3
  !==============================================================================
    !> Apply wall Boundary Conditions at jmin: slip wall
    !> - 3D curvilinear version - (advancement of rho on boundary)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,k
    real(wp) :: rhovtan_ksi,rhovtan_phi
    ! ---------------------------------------------------------------------------

    ! Impose wall conditions on conservative variables
    ! ================================================
    do k=1,nz
       do i=1,nx
          ! define tangential velocity components
          rhovtan_ksi= rhou_n(i,1,k)*tx_ksi_jmin(i,k) &
                     + rhov_n(i,1,k)*ty_ksi_jmin(i,k) &
                     + rhow_n(i,1,k)*tz_ksi_jmin(i,k)
          rhovtan_phi= rhou_n(i,1,k)*tx_phi_jmin(i,k) &
                     + rhov_n(i,1,k)*ty_phi_jmin(i,k) &
                     + rhow_n(i,1,k)*tz_phi_jmin(i,k)
          
          ! slip condition: nullity of normal velocity component
          rhou_n(i,1,k)=rhovtan_ksi*txn_ksi_jmin(i,k)+rhovtan_phi*txn_phi_jmin(i,k)
          rhov_n(i,1,k)=rhovtan_ksi*tyn_ksi_jmin(i,k)+rhovtan_phi*tyn_phi_jmin(i,k)
          rhow_n(i,1,k)=rhovtan_ksi*tzn_ksi_jmin(i,k)+rhovtan_phi*tzn_phi_jmin(i,k)
       enddo
    enddo

  end subroutine bc_wall_jmin_slip_c3

  !==============================================================================
  module subroutine bc_wall_jmin_adiabatic_c3
  !==============================================================================
    !> Apply wall Boundary Conditions at jmin: no slip + adiabatic
    !> - 3D curvilinear version - (advancement of rho on boundary)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,k
    ! ---------------------------------------------------------------------------

    ! Impose wall conditions on conservative variables
    ! ================================================
    do k=1,nz
       do i=1,nx
          ! no slip conditions
          rhou_n(i,1,k)=0.0_wp
          rhov_n(i,1,k)=0.0_wp
          rhow_n(i,1,k)=0.0_wp
       enddo
    enddo

    do k=ndz_jmin,nfz_jmin
       do i=ndx_jmin,nfx_jmin
          ! dT/dn=0
          T_wall= onethird*( 4.0_wp*Tmp(i,2,k)-Tmp(i,3,k)         &
                + (Tmp(i+1,1,k)-Tmp(i-1,1,k))*gksigeta3_jmin(i,k) &
                + (Tmp(i,1,k+1)-Tmp(i,1,k-1))*gphigeta3_jmin(i,k) )

          rhoe_n(i,1,k)=rho_n(i,1,k)*ecalc_tro(T_wall,rho_n(i,1,k))
       enddo
    enddo

    ! edge with other BC at imin
    ! --------------------------
    if (ndx_jmin==2) then
       i=1
       do k=ndz_jmin,nfz_jmin
          ! dT/dn=0
          T_wall= onethird*( 4.0_wp*Tmp(i,2,k)-Tmp(i,3,k)        &
                +2.0_wp*(Tmp(i+1,1,k)-Tmp(i,1,k))*gksigeta3_jmin(i,k) &
                +(Tmp(i,1,k+1)-Tmp(i,1,k-1))*gphigeta3_jmin(i,k) )

          rhoe_n(i,1,k)=rho_n(i,1,k)*ecalc_tro(T_wall,rho_n(i,1,k))
       enddo
    endif

    ! edge with other BC at imax
    ! --------------------------
    if (ndx_jmin==nx-1) then
       i=nx
       do k=ndz_jmin,nfz_jmin
          ! dT/dn=0
          T_wall= onethird*( 4.0_wp*Tmp(i,2,k)-Tmp(i,3,k)        &
                +2.0_wp*(Tmp(i,1,k)-Tmp(i-1,1,k))*gksigeta3_jmin(i,k) &
                +(Tmp(i,1,k+1)-Tmp(i,1,k-1))*gphigeta3_jmin(i,k) )

          rhoe_n(i,1,k)=rho_n(i,1,k)*ecalc_tro(T_wall,rho_n(i,1,k))
       enddo
    endif

    ! edge with other BC at kmin
    ! --------------------------
    if (ndz_jmin==2) then
       k=1
       do i=ndx_jmin,nfx_jmin
          ! dT/dn=0
          T_wall= onethird*( 4.0_wp*Tmp(i,2,k)-Tmp(i,3,k)        &
                +(Tmp(i+1,1,k)-Tmp(i-1,1,k))*gksigeta3_jmin(i,k) &
                +2.0_wp*(Tmp(i,1,k+1)-Tmp(i,1,k))*gphigeta3_jmin(i,k) )

          rhoe_n(i,1,k)=rho_n(i,1,k)*ecalc_tro(T_wall,rho_n(i,1,k))
       enddo
    endif

    ! edge with other BC at kmax
    ! --------------------------
    if (ndz_jmin==nz-1) then
       k=nz
       do i=ndx_jmin,nfx_jmin
          ! dT/dn=0
          T_wall= onethird*( 4.0_wp*Tmp(i,2,k)-Tmp(i,3,k)        &
                +(Tmp(i+1,1,k)-Tmp(i-1,1,k))*gksigeta3_jmin(i,k) &
                +2.0_wp*(Tmp(i,1,k)-Tmp(i,1,k-1))*gphigeta3_jmin(i,k) )

          rhoe_n(i,1,k)=rho_n(i,1,k)*ecalc_tro(T_wall,rho_n(i,1,k))
       enddo
    endif

    ! corner with other BC at imin-kmin
    ! ---------------------------------
    if ((ndx_jmin==2).and.(ndz_jmin==2)) then
       i=1
       k=1
       ! dT/dn=0
       T_wall= onethird*( 4.0_wp*Tmp(i,2,k)-Tmp(i,3,k)        &
             +2.0_wp*(Tmp(i+1,1,k)-Tmp(i,1,k))*gksigeta3_jmin(i,k) &
             +2.0_wp*(Tmp(i,1,k+1)-Tmp(i,1,k))*gphigeta3_jmin(i,k) )

       rhoe_n(i,1,k)=rho_n(i,1,k)*ecalc_tro(T_wall,rho_n(i,1,k))
    endif
    
    ! corner with other BC at imin-kmax
    ! ---------------------------------
    if ((ndx_jmin==2).and.(ndz_jmin==nz-1)) then
       i=1
       k=nz
       ! dT/dn=0
       T_wall= onethird*( 4.0_wp*Tmp(i,2,k)-Tmp(i,3,k)        &
             +2.0_wp*(Tmp(i+1,1,k)-Tmp(i,1,k))*gksigeta3_jmin(i,k) &
             +2.0_wp*(Tmp(i,1,k)-Tmp(i,1,k-1))*gphigeta3_jmin(i,k) )

       rhoe_n(i,1,k)=rho_n(i,1,k)*ecalc_tro(T_wall,rho_n(i,1,k))
    endif
    
    ! corner with other BC at imax-kmin
    ! ---------------------------------
    if ((ndx_jmin==nx-1).and.(ndz_jmin==2)) then
       i=nx
       k=1
       ! dT/dn=0
       T_wall= onethird*( 4.0_wp*Tmp(i,2,k)-Tmp(i,3,k)        &
             +2.0_wp*(Tmp(i,1,k)-Tmp(i-1,1,k))*gksigeta3_jmin(i,k) &
             +2.0_wp*(Tmp(i,1,k+1)-Tmp(i,1,k))*gphigeta3_jmin(i,k) )

       rhoe_n(i,1,k)=rho_n(i,1,k)*ecalc_tro(T_wall,rho_n(i,1,k))
    endif

    ! corner with other BC at imax-kmax
    ! ---------------------------------
    if ((ndx_jmin==nx-1).and.(ndz_jmin==nz-1)) then
       i=nx
       k=nz
       ! dT/dn=0
       T_wall= onethird*( 4.0_wp*Tmp(i,2,k)-Tmp(i,3,k)        &
             +2.0_wp*(Tmp(i,1,k)-Tmp(i-1,1,k))*gksigeta3_jmin(i,k) &
             +2.0_wp*(Tmp(i,1,k)-Tmp(i,1,k-1))*gphigeta3_jmin(i,k) )

       rhoe_n(i,1,k)=rho_n(i,1,k)*ecalc_tro(T_wall,rho_n(i,1,k))
    endif
    
  end subroutine bc_wall_jmin_adiabatic_c3

  !==============================================================================
  module subroutine bc_wall_jmax_slip_c3
  !==============================================================================
    !> Apply wall Boundary Conditions at jmax: slip wall
    !> - 3D curvilinear version - (advancement of rho on boundary)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,k
    real(wp) :: rhovtan_ksi,rhovtan_phi
    ! ---------------------------------------------------------------------------

    ! Impose wall conditions on conservative variables
    ! ================================================
    do k=1,nz
       do i=1,nx
          ! define tangential velocity components
          rhovtan_ksi= rhou_n(i,ny,k)*tx_ksi_jmax(i,k) &
                     + rhov_n(i,ny,k)*ty_ksi_jmax(i,k) &
                     + rhow_n(i,ny,k)*tz_ksi_jmax(i,k)
          rhovtan_phi= rhou_n(i,ny,k)*tx_phi_jmax(i,k) &
                     + rhov_n(i,ny,k)*ty_phi_jmax(i,k) &
                     + rhow_n(i,ny,k)*tz_phi_jmax(i,k)
          
          ! slip condition: nullity of normal velocity component
          rhou_n(i,ny,k)=rhovtan_ksi*txn_ksi_jmax(i,k)+rhovtan_phi*txn_phi_jmax(i,k)
          rhov_n(i,ny,k)=rhovtan_ksi*tyn_ksi_jmax(i,k)+rhovtan_phi*tyn_phi_jmax(i,k)
          rhow_n(i,ny,k)=rhovtan_ksi*tzn_ksi_jmax(i,k)+rhovtan_phi*tzn_phi_jmax(i,k)
       enddo
    enddo

  end subroutine bc_wall_jmax_slip_c3

  !==============================================================================
  module subroutine bc_wall_jmax_adiabatic_c3
  !==============================================================================
    !> Apply wall Boundary Conditions at jmax: no slip + adiabatic
    !> - 3D curvilinear version - (advancement of rho on boundary)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,k
    ! ---------------------------------------------------------------------------

    ! Impose wall conditions on conservative variables
    ! ================================================
    do k=1,nz
       do i=1,nx
          ! no slip conditions
          rhou_n(i,ny,k)=0.0_wp
          rhov_n(i,ny,k)=0.0_wp
          rhow_n(i,ny,k)=0.0_wp
       enddo
    enddo

    do k=ndz_jmax,nfz_jmax
       do i=ndx_jmax,nfx_jmax
          ! dT/dn=0
          T_wall=onethird*( 4.0_wp*Tmp(i,ny-1,k)-Tmp(i,ny-2,k)     &
                +(Tmp(i+1,ny,k)-Tmp(i-1,ny,k))*gksigeta3_jmax(i,k) &
                +(Tmp(i,ny,k+1)-Tmp(i,ny,k-1))*gphigeta3_jmax(i,k) )

          rhoe_n(i,ny,k)=rho_n(i,ny,k)*ecalc_tro(T_wall,rho_n(i,ny,k))
       enddo
    enddo

    ! edge with other BC at imin
    ! --------------------------
    if (ndx_jmax==2) then
       i=1
       do k=ndz_jmax,nfz_jmax
          ! dT/dn=0
          T_wall= onethird*( 4.0_wp*Tmp(i,ny-1,k)-Tmp(i,ny-2,k)    &
                +2.0_wp*(Tmp(i+1,ny,k)-Tmp(i,ny,k))*gksigeta3_jmax(i,k) &
                +(Tmp(i,ny,k+1)-Tmp(i,ny,k-1))*gphigeta3_jmax(i,k) )

          rhoe_n(i,ny,k)=rho_n(i,ny,k)*ecalc_tro(T_wall,rho_n(i,ny,k))
       enddo
    endif

    ! edge with other BC at imax
    ! --------------------------
    if (ndx_jmax==nx-1) then
       i=nx
       do k=ndz_jmax,nfz_jmax
          ! dT/dn=0
          T_wall= onethird*( 4.0_wp*Tmp(i,ny-1,k)-Tmp(i,ny-2,k)    &
                +2.0_wp*(Tmp(i,ny,k)-Tmp(i-1,ny,k))*gksigeta3_jmax(i,k) &
                +(Tmp(i,ny,k+1)-Tmp(i,ny,k-1))*gphigeta3_jmax(i,k) )

          rhoe_n(i,ny,k)=rho_n(i,ny,k)*ecalc_tro(T_wall,rho_n(i,ny,k))
       enddo
    endif

    ! edge with other BC at kmin
    ! --------------------------
    if (ndz_jmax==2) then
       k=1
       do i=ndx_jmax,nfx_jmax
          ! dT/dn=0
          T_wall= onethird*( 4.0_wp*Tmp(i,ny-1,k)-Tmp(i,ny-2,k)    &
                +(Tmp(i+1,ny,k)-Tmp(i-1,ny,k))*gksigeta3_jmax(i,k) &
                +2.0_wp*(Tmp(i,ny,k+1)-Tmp(i,ny,k))*gphigeta3_jmax(i,k) )

          rhoe_n(i,ny,k)=rho_n(i,ny,k)*ecalc_tro(T_wall,rho_n(i,ny,k))
       enddo
    endif

    ! edge with other BC at kmax
    ! --------------------------
    if (ndz_jmax==nz-1) then
       k=nz
       do i=ndx_jmax,nfx_jmax
          ! dT/dn=0
          T_wall= onethird*( 4.0_wp*Tmp(i,ny-1,k)-Tmp(i,ny-2,k)    &
                +(Tmp(i+1,ny,k)-Tmp(i-1,ny,k))*gksigeta3_jmax(i,k) &
                +2.0_wp*(Tmp(i,ny,k)-Tmp(i,ny,k-1))*gphigeta3_jmax(i,k) )

          rhoe_n(i,ny,k)=rho_n(i,ny,k)*ecalc_tro(T_wall,rho_n(i,ny,k))
       enddo
    endif

    ! corner with other BC at imin-kmin
    ! ---------------------------------
    if ((ndx_jmax==2).and.(ndz_jmax==2)) then
       i=1
       k=1
       ! dT/dn=0
       T_wall= onethird*( 4.0_wp*Tmp(i,ny-1,k)-Tmp(i,ny-2,k)    &
             +2.0_wp*(Tmp(i+1,ny,k)-Tmp(i,ny,k))*gksigeta3_jmax(i,k) &
             +2.0_wp*(Tmp(i,ny,k+1)-Tmp(i,ny,k))*gphigeta3_jmax(i,k) )

       rhoe_n(i,ny,k)=rho_n(i,ny,k)*ecalc_tro(T_wall,rho_n(i,ny,k))
    endif
    
    ! corner with other BC at imin-kmax
    ! ---------------------------------
    if ((ndx_jmax==2).and.(ndz_jmax==nz-1)) then
       i=1
       k=nz
       ! dT/dn=0
       T_wall= onethird*( 4.0_wp*Tmp(i,ny-1,k)-Tmp(i,ny-2,k)    &
             +2.0_wp*(Tmp(i+1,ny,k)-Tmp(i,ny,k))*gksigeta3_jmax(i,k) &
             +2.0_wp*(Tmp(i,ny,k)-Tmp(i,ny,k-1))*gphigeta3_jmax(i,k) )

       rhoe_n(i,ny,k)=rho_n(i,ny,k)*ecalc_tro(T_wall,rho_n(i,ny,k))
    endif
    
    ! corner with other BC at imax-kmin
    ! ---------------------------------
    if ((ndx_jmax==nx-1).and.(ndz_jmax==2)) then
       i=nx
       k=1
       ! dT/dn=0
       T_wall= onethird*( 4.0_wp*Tmp(i,ny-1,k)-Tmp(i,ny-2,k)    &
             +2.0_wp*(Tmp(i,ny,k)-Tmp(i-1,ny,k))*gksigeta3_jmax(i,k) &
             +2.0_wp*(Tmp(i,ny,k+1)-Tmp(i,ny,k))*gphigeta3_jmax(i,k) )

       rhoe_n(i,ny,k)=rho_n(i,ny,k)*ecalc_tro(T_wall,rho_n(i,ny,k))
    endif

    ! corner with other BC at imax-kmax
    ! ---------------------------------
    if ((ndx_jmax==nx-1).and.(ndz_jmax==nz-1)) then
       i=nx
       k=nz
       ! dT/dn=0
       T_wall= onethird*( 4.0_wp*Tmp(i,ny-1,k)-Tmp(i,ny-2,k)    &
             +2.0_wp*(Tmp(i,ny,k)-Tmp(i-1,ny,k))*gksigeta3_jmax(i,k) &
             +2.0_wp*(Tmp(i,ny,k)-Tmp(i,ny,k-1))*gphigeta3_jmax(i,k) )

       rhoe_n(i,ny,k)=rho_n(i,ny,k)*ecalc_tro(T_wall,rho_n(i,ny,k))
    endif
    
  end subroutine bc_wall_jmax_adiabatic_c3

  !==============================================================================
  module subroutine bc_wall_kmin_slip_c3
  !==============================================================================
    !> Apply wall Boundary Conditions at kmin: slip wall
    !> - 3D curvilinear version - (advancement of rho on boundary)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j
    real(wp) :: rhovtan_ksi,rhovtan_eta
    ! ---------------------------------------------------------------------------

    ! Impose wall conditions on conservative variables
    ! ================================================
    do j=1,ny
       do i=1,nx
          ! define tangential velocity components
          rhovtan_ksi= rhou_n(i,j,1)*tx_ksi_kmin(i,j) &
                     + rhov_n(i,j,1)*ty_ksi_kmin(i,j) &
                     + rhow_n(i,j,1)*tz_ksi_kmin(i,j)
          rhovtan_eta= rhou_n(i,j,1)*tx_eta_kmin(i,j) &
                     + rhov_n(i,j,1)*ty_eta_kmin(i,j) &
                     + rhow_n(i,j,1)*tz_eta_kmin(i,j)
          
          ! slip condition: nullity of normal velocity component
          rhou_n(i,j,1)=rhovtan_ksi*txn_ksi_kmin(i,j)+rhovtan_eta*txn_eta_kmin(i,j)
          rhov_n(i,j,1)=rhovtan_ksi*tyn_ksi_kmin(i,j)+rhovtan_eta*tyn_eta_kmin(i,j)
          rhow_n(i,j,1)=rhovtan_ksi*tzn_ksi_kmin(i,j)+rhovtan_eta*tzn_eta_kmin(i,j)
       enddo
    enddo

  end subroutine bc_wall_kmin_slip_c3

  !==============================================================================
  module subroutine bc_wall_kmin_adiabatic_c3
  !==============================================================================
    !> Apply wall Boundary Conditions at kmin: no slip + adiabatic
    !> - 3D curvilinear version - (advancement of rho on boundary)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j
    ! ---------------------------------------------------------------------------

    ! Impose wall conditions on conservative variables
    ! ================================================
    do j=1,ny
       do i=1,nx
          ! no slip conditions
          rhou_n(i,j,1)=0.0_wp
          rhov_n(i,j,1)=0.0_wp
          rhow_n(i,j,1)=0.0_wp
       enddo
    enddo

    do j=ndy_kmin,nfy_kmin
       do i=ndx_kmin,nfx_kmin
          ! dT/dn=0
          T_wall=onethird*( 4.0_wp*Tmp(i,j,2)-Tmp(i,j,3)         &
                +(Tmp(i+1,j,1)-Tmp(i-1,j,1))*gksigphi3_kmin(i,j) &
                +(Tmp(i,j+1,1)-Tmp(i,j-1,1))*getagphi3_kmin(i,j) )

          rhoe_n(i,j,1)=rho_n(i,j,1)*ecalc_tro(T_wall,rho_n(i,j,1))
       enddo
    enddo

    ! edge with other BC at imin
    ! --------------------------
    if (ndx_kmin==2) then
       i=1
       do j=ndy_kmin,nfy_kmin
          ! dT/dn=0
          T_wall= onethird*( 4.0_wp*Tmp(i,j,2)-Tmp(i,j,3)        &
                +2.0_wp*(Tmp(i+1,j,1)-Tmp(i,j,1))*gksigphi3_kmin(i,j) &
                +(Tmp(i,j+1,1)-Tmp(i,j-1,1))*getagphi3_kmin(i,j) )

          rhoe_n(i,j,1)=rho_n(i,j,1)*ecalc_tro(T_wall,rho_n(i,j,1))
       enddo
    endif

    ! edge with other BC at imax
    ! --------------------------
    if (ndx_kmin==nx-1) then
       i=nx
       do j=ndy_kmin,nfy_kmin
          ! dT/dn=0
          T_wall= onethird*( 4.0_wp*Tmp(i,j,2)-Tmp(i,j,3)        &
                +2.0_wp*(Tmp(i,j,1)-Tmp(i-1,j,1))*gksigphi3_kmin(i,j) &
                +(Tmp(i,j+1,1)-Tmp(i,j-1,1))*getagphi3_kmin(i,j) )

          rhoe_n(i,j,1)=rho_n(i,j,1)*ecalc_tro(T_wall,rho_n(i,j,1))
       enddo
    endif

    ! edge with other BC at jmin
    ! --------------------------
    if (ndy_kmin==2) then
       j=1
       do i=ndx_kmin,nfx_kmin
          ! dT/dn=0
          T_wall= onethird*( 4.0_wp*Tmp(i,j,2)-Tmp(i,j,3)        &
                +(Tmp(i+1,j,1)-Tmp(i-1,j,1))*gksigphi3_kmin(i,j) &
                +2.0_wp*(Tmp(i,j+1,1)-Tmp(i,j,1))*getagphi3_kmin(i,j) )

          rhoe_n(i,j,1)=rho_n(i,j,1)*ecalc_tro(T_wall,rho_n(i,j,1))
       enddo
    endif

    ! edge with other BC at jmax
    ! --------------------------
    if (ndy_kmin==ny-1) then
       j=ny
       do i=ndx_kmin,nfx_kmin
          ! dT/dn=0
          T_wall= onethird*( 4.0_wp*Tmp(i,j,2)-Tmp(i,j,3)        &
                +(Tmp(i+1,j,1)-Tmp(i-1,j,1))*gksigphi3_kmin(i,j) &
                +2.0_wp*(Tmp(i,j,1)-Tmp(i,j-1,1))*getagphi3_kmin(i,j) )

          rhoe_n(i,j,1)=rho_n(i,j,1)*ecalc_tro(T_wall,rho_n(i,j,1))
       enddo
    endif

    ! corner with other BC at imin-jmin
    ! ---------------------------------
    if ((ndx_kmin==2).and.(ndy_kmin==2)) then
       i=1
       j=1
       ! dT/dn=0
       T_wall= onethird*( 4.0_wp*Tmp(i,j,2)-Tmp(i,j,3)        &
             +2.0_wp*(Tmp(i+1,j,1)-Tmp(i,j,1))*gksigphi3_kmin(i,j) &
             +2.0_wp*(Tmp(i,j+1,1)-Tmp(i,j,1))*getagphi3_kmin(i,j) )

       rhoe_n(i,j,1)=rho_n(i,j,1)*ecalc_tro(T_wall,rho_n(i,j,1))
    endif
    
    ! corner with other BC at imin-jmax
    ! ---------------------------------
    if ((ndx_kmin==2).and.(ndy_kmin==ny-1)) then
       i=1
       j=ny
       ! dT/dn=0
       T_wall= onethird*( 4.0_wp*Tmp(i,j,2)-Tmp(i,j,3)        &
             +2.0_wp*(Tmp(i+1,j,1)-Tmp(i,j,1))*gksigphi3_kmin(i,j) &
             +2.0_wp*(Tmp(i,j,1)-Tmp(i,j-1,1))*getagphi3_kmin(i,j) )

       rhoe_n(i,j,1)=rho_n(i,j,1)*ecalc_tro(T_wall,rho_n(i,j,1))
    endif

    ! corner with other BC at imax-jmin
    ! ---------------------------------
    if ((ndx_kmin==nx-1).and.(ndy_kmin==2)) then
       i=nx
       j=1
       ! dT/dn=0
       T_wall= onethird*( 4.0_wp*Tmp(i,j,2)-Tmp(i,j,3)        &
             +2.0_wp*(Tmp(i,j,1)-Tmp(i-1,j,1))*gksigphi3_kmin(i,j) &
             +2.0_wp*(Tmp(i,j+1,1)-Tmp(i,j,1))*getagphi3_kmin(i,j) )

       rhoe_n(i,j,1)=rho_n(i,j,1)*ecalc_tro(T_wall,rho_n(i,j,1))
    endif
        
    ! corner with other BC at imax-jmax
    ! ---------------------------------
    if ((ndx_kmin==nx-1).and.(ndy_kmin==ny-1)) then
       i=nx
       j=ny
       ! dT/dn=0
       T_wall= onethird*( 4.0_wp*Tmp(i,j,2)-Tmp(i,j,3)        &
             +2.0_wp*(Tmp(i,j,1)-Tmp(i-1,j,1))*gksigphi3_kmin(i,j) &
             +2.0_wp*(Tmp(i,j,1)-Tmp(i,j-1,1))*getagphi3_kmin(i,j) )

       rhoe_n(i,j,1)=rho_n(i,j,1)*ecalc_tro(T_wall,rho_n(i,j,1))
    endif

  end subroutine bc_wall_kmin_adiabatic_c3

  !==============================================================================
  module subroutine bc_wall_kmax_slip_c3
  !==============================================================================
    !> Apply wall Boundary Conditions at kmax: slip wall
    !> - 3D curvilinear version - (advancement of rho on boundary)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j
    real(wp) :: rhovtan_ksi,rhovtan_eta
    ! ---------------------------------------------------------------------------

    ! Impose wall conditions on conservative variables
    ! ================================================
    do j=1,ny
       do i=1,nx
          ! define tangential velocity components
          rhovtan_ksi= rhou_n(i,j,nz)*tx_ksi_kmax(i,j) &
                     + rhov_n(i,j,nz)*ty_ksi_kmax(i,j) &
                     + rhow_n(i,j,nz)*tz_ksi_kmax(i,j)
          rhovtan_eta= rhou_n(i,j,nz)*tx_eta_kmax(i,j) &
                     + rhov_n(i,j,nz)*ty_eta_kmax(i,j) &
                     + rhow_n(i,j,nz)*tz_eta_kmax(i,j)
          
          ! slip condition: nullity of normal velocity component
          rhou_n(i,j,nz)=rhovtan_ksi*txn_ksi_kmax(i,j)+rhovtan_eta*txn_eta_kmax(i,j)
          rhov_n(i,j,nz)=rhovtan_ksi*tyn_ksi_kmax(i,j)+rhovtan_eta*tyn_eta_kmax(i,j)
          rhow_n(i,j,nz)=rhovtan_ksi*tzn_ksi_kmax(i,j)+rhovtan_eta*tzn_eta_kmax(i,j)
       enddo
    enddo

  end subroutine bc_wall_kmax_slip_c3

  !==============================================================================
  module subroutine bc_wall_kmax_adiabatic_c3
  !==============================================================================
    !> Apply wall Boundary Conditions at kmax: no slip + adiabatic
    !> - 3D curvilinear version - (advancement of rho on boundary)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j
    ! ---------------------------------------------------------------------------

    ! Impose wall conditions on conservative variables
    ! ================================================
    do j=1,ny
       do i=1,nx
          ! no slip conditions
          rhou_n(i,j,nz)=0.0_wp
          rhov_n(i,j,nz)=0.0_wp
          rhow_n(i,j,nz)=0.0_wp
       enddo
    enddo

    do j=ndy_kmax,nfy_kmax
       do i=ndx_kmax,nfx_kmax
          ! dT/dn=0
          T_wall=onethird*( 4.0_wp*Tmp(i,j,2)-Tmp(i,j,3)         &
                +(Tmp(i+1,j,nz)-Tmp(i-1,j,nz))*gksigphi3_kmax(i,j) &
                +(Tmp(i,j+1,nz)-Tmp(i,j-1,nz))*getagphi3_kmax(i,j) )

          rhoe_n(i,j,nz)=rho_n(i,j,nz)*ecalc_tro(T_wall,rho_n(i,j,nz))
       enddo
    enddo

    ! edge with other BC at imin
    ! --------------------------
    if (ndx_kmax==2) then
       i=1
       do j=ndy_kmax,nfy_kmax
          ! dT/dn=0
          T_wall= onethird*( 4.0_wp*Tmp(i,j,nz-1)-Tmp(i,j,nz-2)    &
                +2.0_wp*(Tmp(i+1,j,nz)-Tmp(i,j,nz))*gksigphi3_kmax(i,j) &
                +(Tmp(i,j+1,nz)-Tmp(i,j-1,nz))*getagphi3_kmax(i,j) )

          rhoe_n(i,j,nz)=rho_n(i,j,nz)*ecalc_tro(T_wall,rho_n(i,j,nz))
       enddo
    endif

    ! edge with other BC at imax
    ! --------------------------
    if (ndx_kmax==nx-1) then
       i=nx
       do j=ndy_kmax,nfy_kmax
          ! dT/dn=0
          T_wall= onethird*( 4.0_wp*Tmp(i,j,nz-1)-Tmp(i,j,nz-2)    &
                +2.0_wp*(Tmp(i,j,nz)-Tmp(i-1,j,nz))*gksigphi3_kmax(i,j) &
                +(Tmp(i,j+1,nz)-Tmp(i,j-1,nz))*getagphi3_kmax(i,j) )

          rhoe_n(i,j,nz)=rho_n(i,j,nz)*ecalc_tro(T_wall,rho_n(i,j,nz))
       enddo
    endif

    ! edge with other BC at jmin
    ! --------------------------
    if (ndy_kmax==2) then
       j=1
       do i=ndx_kmax,nfx_kmax
          ! dT/dn=0
          T_wall= onethird*( 4.0_wp*Tmp(i,j,nz-1)-Tmp(i,j,nz-2)    &
                +(Tmp(i+1,j,nz)-Tmp(i-1,j,nz))*gksigphi3_kmax(i,j) &
                +2.0_wp*(Tmp(i,j+1,nz)-Tmp(i,j,nz))*getagphi3_kmax(i,j) )

          rhoe_n(i,j,nz)=rho_n(i,j,nz)*ecalc_tro(T_wall,rho_n(i,j,nz))
       enddo
    endif

    ! edge with other BC at jmax
    ! --------------------------
    if (ndy_kmax==ny-1) then
       j=ny
       do i=ndx_kmax,nfx_kmax
          ! dT/dn=0
          T_wall= onethird*( 4.0_wp*Tmp(i,j,nz-1)-Tmp(i,j,nz-2)    &
                +(Tmp(i+1,j,nz)-Tmp(i-1,j,nz))*gksigphi3_kmax(i,j) &
                +2.0_wp*(Tmp(i,j,nz)-Tmp(i,j-1,nz))*getagphi3_kmax(i,j) )

          rhoe_n(i,j,nz)=rho_n(i,j,nz)*ecalc_tro(T_wall,rho_n(i,j,nz))
       enddo
    endif

    ! corner with other BC at imin-jmin
    ! ---------------------------------
    if ((ndx_kmax==2).and.(ndy_kmax==2)) then
       i=1
       j=1
       ! dT/dn=0
       T_wall= onethird*( 4.0_wp*Tmp(i,j,nz-1)-Tmp(i,j,nz-2)    &
             +2.0_wp*(Tmp(i+1,j,nz)-Tmp(i,j,nz))*gksigphi3_kmax(i,j) &
             +2.0_wp*(Tmp(i,j+1,nz)-Tmp(i,j,nz))*getagphi3_kmax(i,j) )

       rhoe_n(i,j,nz)=rho_n(i,j,nz)*ecalc_tro(T_wall,rho_n(i,j,nz))
    endif
    
    ! corner with other BC at imin-jmax
    ! ---------------------------------
    if ((ndx_kmax==2).and.(ndy_kmax==ny-1)) then
       i=1
       j=ny
       ! dT/dn=0
       T_wall= onethird*( 4.0_wp*Tmp(i,j,nz-1)-Tmp(i,j,nz-2)    &
             +2.0_wp*(Tmp(i+1,j,nz)-Tmp(i,j,nz))*gksigphi3_kmax(i,j) &
             +2.0_wp*(Tmp(i,j,nz)-Tmp(i,j-1,nz))*getagphi3_kmax(i,j) )

       rhoe_n(i,j,nz)=rho_n(i,j,nz)*ecalc_tro(T_wall,rho_n(i,j,nz))
    endif

    ! corner with other BC at imax-jmin
    ! ---------------------------------
    if ((ndx_kmax==nx-1).and.(ndy_kmax==2)) then
       i=nx
       j=1
       ! dT/dn=0
       T_wall= onethird*( 4.0_wp*Tmp(i,j,nz-1)-Tmp(i,j,nz-2)    &
             +2.0_wp*(Tmp(i,j,nz)-Tmp(i-1,j,nz))*gksigphi3_kmax(i,j) &
             +2.0_wp*(Tmp(i,j+1,nz)-Tmp(i,j,nz))*getagphi3_kmax(i,j) )

       rhoe_n(i,j,nz)=rho_n(i,j,nz)*ecalc_tro(T_wall,rho_n(i,j,nz))
    endif
        
    ! corner with other BC at imax-jmax
    ! ---------------------------------
    if ((ndx_kmax==nx-1).and.(ndy_kmax==ny-1)) then
       i=nx
       j=ny
       ! dT/dn=0
       T_wall= onethird*( 4.0_wp*Tmp(i,j,nz-1)-Tmp(i,j,nz-2)    &
             +2.0_wp*(Tmp(i,j,nz)-Tmp(i-1,j,nz))*gksigphi3_kmax(i,j) &
             +2.0_wp*(Tmp(i,j,nz)-Tmp(i,j-1,nz))*getagphi3_kmax(i,j) )

       rhoe_n(i,j,nz)=rho_n(i,j,nz)*ecalc_tro(T_wall,rho_n(i,j,nz))
    endif

  end subroutine bc_wall_kmax_adiabatic_c3

end submodule smod_bc_wall_2_c3
