!================================================================================
submodule (mod_bc_wall) smod_bc_wall_2_c
!================================================================================
  !> Submodule to apply wall Boundary Conditions on conservative variables
  !> - curvilinear version - (advancement of rho on boundary)
!================================================================================

contains

  !==============================================================================
  module subroutine bc_wall_imin_slip_c
  !==============================================================================
    !> Apply wall Boundary Conditions at imin: slip wall
    !> - curvilinear version - (advancement of rho on boundary)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: j,k
    real(wp) :: rhovtan
    ! ---------------------------------------------------------------------------

    ! Impose wall conditions on conservative variables
    ! ================================================
    ! -> slip condition: nullity of normal velocity component
    do k=1,nz
       do j=1,ny
          ! compute tangent velocity component
          rhovtan=rhou_n(1,j,k)*nyn_imin(j,k)-rhov_n(1,j,k)*nxn_imin(j,k)
          ! suppress normal component by projecting only tangent component
          rhou_n(1,j,k)= rhovtan*nyn_imin(j,k)
          rhov_n(1,j,k)=-rhovtan*nxn_imin(j,k)
       enddo
    enddo

  end subroutine bc_wall_imin_slip_c

  !==============================================================================
  module subroutine bc_wall_imin_adiabatic_c
  !==============================================================================
    !> Apply wall Boundary Conditions at imin: no slip + adiabatic
    !> - curvilinear version - (advancement of rho on boundary)
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
    
    do k=1,nz
       do j=ndy_imin,nfy_imin
          ! dT/dn=0
          T_wall= onethird*(4.0_wp*Tmp(2,j,k)-Tmp(3,j,k) &
                + (Tmp(1,j+1,k)-Tmp(1,j-1,k))*gksigeta_imin(j))

          rhoe_n(1,j,k)=rho_n(1,j,k)*ecalc_tro(T_wall,rho_n(1,j,k))
       enddo
    enddo

    ! corner with other BC at jmin
    ! ----------------------------
    if (ndy_imin==2) then
       j=1
       do k=1,nz
          ! dT/dn=0
          T_wall= onethird*(4.0_wp*Tmp(2,j,k)-Tmp(3,j,k) &
                + 2.0_wp*(Tmp(1,j+1,k)-Tmp(1,j,k))*gksigeta_imin(j))

          rhoe_n(1,j,k)=rho_n(1,j,k)*ecalc_tro(T_wall,rho_n(1,j,k))
       enddo
    endif

    ! corner with other BC at jmax
    ! ----------------------------
    if (nfy_imin==ny-1) then
       j=ny
       do k=1,nz
          ! dT/dn=0
          T_wall= onethird*(4.0_wp*Tmp(2,j,k)-Tmp(3,j,k) &
                + 2.0_wp*(Tmp(1,j,k)-Tmp(1,j-1,k))*gksigeta_imin(j))

          rhoe_n(1,j,k)=rho_n(1,j,k)*ecalc_tro(T_wall,rho_n(1,j,k))
       enddo
    endif

  end subroutine bc_wall_imin_adiabatic_c

  !==============================================================================
  module subroutine bc_wall_imax_slip_c
  !==============================================================================
    !> Apply wall Boundary Conditions at imax: slip wall
    !> - curvilinear version - (advancement of rho on boundary)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: j,k
    real(wp) :: rhovtan
    ! ---------------------------------------------------------------------------

    ! Impose wall conditions on conservative variables
    ! ================================================
    ! -> slip condition: nullity of normal velocity component
    do k=1,nz
       do j=1,ny
          ! compute tangent velocity component
          rhovtan=rhou_n(nx,j,k)*nyn_imax(j,k)-rhov_n(nx,j,k)*nxn_imax(j,k)
          ! suppress normal component by projecting only tangent component
          rhou_n(nx,j,k)= rhovtan*nyn_imax(j,k)
          rhov_n(nx,j,k)=-rhovtan*nxn_imax(j,k)
       enddo
    enddo

  end subroutine bc_wall_imax_slip_c

  !==============================================================================
  module subroutine bc_wall_imax_adiabatic_c
  !==============================================================================
    !> Apply wall Boundary Conditions at imax: no slip + adiabatic
    !> - curvilinear version - (advancement of rho on boundary)
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

    do k=1,nz
       do j=ndy_imax,nfy_imax
          ! dT/dn=0
          T_wall= onethird*(4.0_wp*Tmp(nx-1,j,k)-Tmp(nx-2,j,k) &
                + (Tmp(nx,j+1,k)-Tmp(nx,j-1,k))*gksigeta_imax(j))

          rhoe_n(nx,j,k)=rho_n(nx,j,k)*ecalc_tro(T_wall,rho_n(nx,j,k))
       enddo
    enddo

    ! corner with other BC at jmin
    ! ----------------------------
    if (ndy_imax==2) then
       j=1
       do k=1,nz
          ! dT/dn=0
          T_wall= onethird*(4.0_wp*Tmp(nx-1,j,k)-Tmp(nx-2,j,k) &
                + 2.0_wp*(Tmp(nx,j+1,k)-Tmp(nx,j,k))*gksigeta_imax(j))

          rhoe_n(nx,j,k)=rho_n(nx,j,k)*ecalc_tro(T_wall,rho_n(nx,j,k))
       enddo
    endif

    ! corner with other BC at jmax
    ! ----------------------------
    if (nfy_imax==ny-1) then
       j=ny
       do k=1,nz
          ! dT/dn=0
          T_wall= onethird*(4.0_wp*Tmp(nx-1,j,k)-Tmp(nx-2,j,k) &
                + 2.0_wp*(Tmp(nx,j,k)-Tmp(nx,j-1,k))*gksigeta_imax(j))

          rhoe_n(nx,j,k)=rho_n(nx,j,k)*ecalc_tro(T_wall,rho_n(nx,j,k))
       enddo
    endif
    
  end subroutine bc_wall_imax_adiabatic_c

  !==============================================================================
  module subroutine bc_wall_jmin_slip_c
  !==============================================================================
    !> Apply wall Boundary Conditions at jmin: slip wall
    !> - curvilinear version - (advancement of rho on boundary)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,k
    real(wp) :: rhovtan
    ! ---------------------------------------------------------------------------

    ! Impose wall conditions on conservative variables
    ! ================================================
    ! -> slip condition: nullity of normal velocity component
    do k=1,nz
       do i=1,nx
          ! compute tangent velocity component
          rhovtan=rhou_n(i,1,k)*nyn_jmin(i,k)-rhov_n(i,1,k)*nxn_jmin(i,k)
          ! suppress normal component by projecting only tangent component
          rhou_n(i,1,k)= rhovtan*nyn_jmin(i,k)
          rhov_n(i,1,k)=-rhovtan*nxn_jmin(i,k)
       enddo
    enddo

  end subroutine bc_wall_jmin_slip_c

  !==============================================================================
  module subroutine bc_wall_jmin_adiabatic_c
  !==============================================================================
    !> Apply wall Boundary Conditions at jmin: no slip + adiabatic
    !> - curvilinear version - (advancement of rho on boundary)
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

    do k=1,nz
       do i=ndx_jmin,nfx_jmin
          ! dT/dn=0
          T_wall= onethird*(4.0_wp*Tmp(i,2,k)-Tmp(i,3,k) &
                + (Tmp(i+1,1,k)-Tmp(i-1,1,k))*gksigeta_jmin(i))

          rhoe_n(i,1,k)=rho_n(i,1,k)*ecalc_tro(T_wall,rho_n(i,1,k))
       enddo
    enddo

    ! corner with other BC at imin
    ! ----------------------------
    if (ndx_jmin==2) then
       i=1
       do k=1,nz
          ! dT/dn=0
          T_wall= onethird*(4.0_wp*Tmp(i,2,k)-Tmp(i,3,k) &
                + 2.0_wp*(Tmp(i+1,1,k)-Tmp(i,1,k))*gksigeta_jmin(i))

          rhoe_n(i,1,k)=rho_n(i,1,k)*ecalc_tro(T_wall,rho_n(i,1,k))
       enddo
    endif

    ! corner with other BC at imax
    ! ----------------------------
    if (nfx_jmin==nx-1) then
       i=nx
       do k=1,nz
          ! dT/dn=0
          T_wall= onethird*(4.0_wp*Tmp(i,2,k)-Tmp(i,3,k) &
                + 2.0_wp*(Tmp(i,1,k)-Tmp(i-1,1,k))*gksigeta_jmin(i))

          rhoe_n(i,1,k)=rho_n(i,1,k)*ecalc_tro(T_wall,rho_n(i,1,k))
       enddo
    endif

  end subroutine bc_wall_jmin_adiabatic_c

  !==============================================================================
  module subroutine bc_wall_jmax_slip_c
  !==============================================================================
    !> Apply wall Boundary Conditions at jmax: slip wall
    !> - curvilinear version - (advancement of rho on boundary)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,k
    real(wp) :: rhovtan
    ! ---------------------------------------------------------------------------

    ! Impose wall conditions on conservative variables
    ! ================================================
    ! -> slip condition: nullity of normal velocity component
    do k=1,nz
       do i=1,nx
          ! compute tangent velocity component
          rhovtan=rhou_n(i,ny,k)*nyn_jmax(i,k)-rhov_n(i,ny,k)*nxn_jmax(i,k)
          ! suppress normal component by projecting only tangent component
          rhou_n(i,ny,k)= rhovtan*nyn_jmax(i,k)
          rhov_n(i,ny,k)=-rhovtan*nxn_jmax(i,k)
       enddo
    enddo

  end subroutine bc_wall_jmax_slip_c

  !==============================================================================
  module subroutine bc_wall_jmax_adiabatic_c
  !==============================================================================
    !> Apply wall Boundary Conditions at jmax: no slip + adiabatic
    !> - curvilinear version - (advancement of rho on boundary)
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

    do k=1,nz
       do i=ndx_jmax,nfx_jmax
          ! dT/dn=0
          T_wall= onethird*(4.0_wp*Tmp(i,ny-1,k)-Tmp(i,ny-2,k) &
                + (Tmp(i+1,ny,k)-Tmp(i-1,ny,k))*gksigeta_jmax(i))

          rhoe_n(i,ny,k)=rho_n(i,ny,k)*ecalc_tro(T_wall,rho_n(i,ny,k))
       enddo
    enddo

    ! corner with other BC at imin
    ! ----------------------------
    if (ndx_jmax==2) then
       i=1
       do k=1,nz
          ! dT/dn=0
          T_wall= onethird*(4.0_wp*Tmp(i,ny-1,k)-Tmp(i,ny-2,k) &
                + 2.0_wp*(Tmp(i+1,ny,k)-Tmp(i,ny,k))*gksigeta_jmax(i))

          rhoe_n(i,ny,k)=rho_n(i,ny,k)*ecalc_tro(T_wall,rho_n(i,ny,k))
       enddo
    endif

    ! corner with other BC at imax
    ! ----------------------------
    if (nfx_jmax==nx-1) then
       i=nx
       do k=1,nz
          ! dT/dn=0
          T_wall= onethird*(4.0_wp*Tmp(i,ny-1,k)-Tmp(i,ny-2,k) &
                + 2.0_wp*(Tmp(i,ny,k)-Tmp(i-1,ny,k))*gksigeta_jmax(i))

          rhoe_n(i,ny,k)=rho_n(i,ny,k)*ecalc_tro(T_wall,rho_n(i,ny,k))
       enddo
    endif

  end subroutine bc_wall_jmax_adiabatic_c

end submodule smod_bc_wall_2_c
