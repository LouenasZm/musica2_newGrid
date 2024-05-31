!================================================================================
submodule (mod_bc_wall) smod_bc_wall_c
!================================================================================
  !> Submodule to apply wall Boundary Conditions on conservative variables
  !> - curvilinear version - (all variables are imposed)
!================================================================================

contains

  !==============================================================================
  module subroutine bc_wall_imin_adiabatic_dpdn_c
  !==============================================================================
    !> Apply wall Boundary Conditions at imin: adiabatic + dp/dn=0
    !> - curvilinear version - (all variables are imposed)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: j,k
    real(wp) :: ro_wall,p_wall
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
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(2,j,k)-prs(3,j,k) &
                + (prs(1,j+1,k)-prs(1,j-1,k))*gksigeta_imin(j))
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(2,j,k))

          ! conservative variables
          rho_n(1,j,k) =ro_wall
          rhoe_n(1,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(2,j,k)-prs(3,j,k) &
                + 2.0_wp*(prs(1,j+1,k)-prs(1,j,k))*gksigeta_imin(j))
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(2,j,k))

          ! conservative variables
          rho_n(1,j,k) =ro_wall
          rhoe_n(1,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(2,j,k)-prs(3,j,k) &
                + 2.0_wp*(prs(1,j,k)-prs(1,j-1,k))*gksigeta_imin(j))
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(2,j,k))

          ! conservative variables
          rho_n(1,j,k) =ro_wall
          rhoe_n(1,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! Set increments to zero
    ! ======================
    do k=ndz_e,nfz_e
       Krho(1,ndy_e:nfy_e,k)=0.0_wp
       Krhou(1,ndy_e:nfy_e,k)=0.0_wp
       Krhov(1,ndy_e:nfy_e,k)=0.0_wp
       Krhow(1,ndy_e:nfy_e,k)=0.0_wp
       Krhoe(1,ndy_e:nfy_e,k)=0.0_wp
    enddo
    
  end subroutine bc_wall_imin_adiabatic_dpdn_c

  !==============================================================================
  module subroutine bc_wall_imin_isotherm_dpdn_c
  !==============================================================================
    !> Apply wall Boundary Conditions at imin: isothermal + dp/dn=0
    !> - curvilinear version - (all variables are imposed)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: j,k
    real(wp) :: ro_wall,p_wall
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
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(2,j,k)-prs(3,j,k) &
                + (prs(1,j+1,k)-prs(1,j-1,k))*gksigeta_imin(j))
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(2,j,k))

          ! conservative variables
          rho_n(1,j,k) =ro_wall
          rhoe_n(1,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    enddo

    ! corner with other BC at jmin
    ! ----------------------------
    if (ndy_imin==2) then
       j=1
       do k=1,nz
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(2,j,k)-prs(3,j,k) &
                + 2.0_wp*(prs(1,j+1,k)-prs(1,j,k))*gksigeta_imin(j))
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(2,j,k))

          ! conservative variables
          rho_n(1,j,k) =ro_wall
          rhoe_n(1,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! corner with other BC at jmax
    ! ----------------------------
    if (nfy_imin==ny-1) then
       j=ny
       do k=1,nz
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(2,j,k)-prs(3,j,k) &
                + 2.0_wp*(prs(1,j,k)-prs(1,j-1,k))*gksigeta_imin(j))
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(2,j,k))

          ! conservative variables
          rho_n(1,j,k) =ro_wall
          rhoe_n(1,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! Set increments to zero
    ! ======================
    do k=ndz_e,nfz_e
       Krho(1,ndy_e:nfy_e,k)=0.0_wp
       Krhou(1,ndy_e:nfy_e,k)=0.0_wp
       Krhov(1,ndy_e:nfy_e,k)=0.0_wp
       Krhow(1,ndy_e:nfy_e,k)=0.0_wp
       Krhoe(1,ndy_e:nfy_e,k)=0.0_wp
    enddo
    
  end subroutine bc_wall_imin_isotherm_dpdn_c

  !==============================================================================
  module subroutine bc_wall_imax_adiabatic_dpdn_c
  !==============================================================================
    !> Apply wall Boundary Conditions at imax: adiabatic + dp/dn=0
    !> - curvilinear version - (all variables are imposed)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: j,k
    real(wp) :: ro_wall,p_wall
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
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(nx-1,j,k)-prs(nx-2,j,k) &
                + (prs(nx,j+1,k)-prs(nx,j-1,k))*gksigeta_imax(j))
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(nx-1,j,k))

          ! conservative variables
          rho_n(nx,j,k) =ro_wall
          rhoe_n(nx,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(nx-1,j,k)-prs(nx-2,j,k) &
                + 2.0_wp*(prs(nx,j+1,k)-prs(nx,j,k))*gksigeta_imax(j))
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(nx-1,j,k))

          ! conservative variables
          rho_n(nx,j,k) =ro_wall
          rhoe_n(nx,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(nx-1,j,k)-prs(nx-2,j,k) &
                + 2.0_wp*(prs(nx,j,k)-prs(nx,j-1,k))*gksigeta_imax(j))
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(nx-1,j,k))

          ! conservative variables
          rho_n(nx,j,k) =ro_wall
          rhoe_n(nx,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! Set increments to zero
    ! ======================
    do k=ndz_e,nfz_e
       Krho(nx,ndy_e:nfy_e,k)=0.0_wp
       Krhou(nx,ndy_e:nfy_e,k)=0.0_wp
       Krhov(nx,ndy_e:nfy_e,k)=0.0_wp
       Krhow(nx,ndy_e:nfy_e,k)=0.0_wp
       Krhoe(nx,ndy_e:nfy_e,k)=0.0_wp
    enddo

  end subroutine bc_wall_imax_adiabatic_dpdn_c

  !==============================================================================
  module subroutine bc_wall_imax_isotherm_dpdn_c
  !==============================================================================
    !> Apply wall Boundary Conditions at imax: isothermal + dp/dn=0
    !> - curvilinear version - (all variables are imposed)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: j,k
    real(wp) :: ro_wall,p_wall
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
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(nx-1,j,k)-prs(nx-2,j,k) &
                + (prs(nx,j+1,k)-prs(nx,j-1,k))*gksigeta_imax(j))
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(nx-1,j,k))

          ! conservative variables
          rho_n(nx,j,k) =ro_wall
         rhoe_n(nx,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    enddo

    ! corner with other BC at jmin
    ! ----------------------------
    if (ndy_imax==2) then
       j=1
       do k=1,nz
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(nx-1,j,k)-prs(nx-2,j,k) &
                + 2.0_wp*(prs(nx,j+1,k)-prs(nx,j,k))*gksigeta_imax(j))
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(nx-1,j,k))

          ! conservative variables
          rho_n(nx,j,k) =ro_wall
          rhoe_n(nx,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! corner with other BC at jmax
    ! ----------------------------
    if (nfy_imax==ny-1) then
       j=ny
       do k=1,nz
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(nx-1,j,k)-prs(nx-2,j,k) &
                + 2.0_wp*(prs(nx,j,k)-prs(nx,j-1,k))*gksigeta_imax(j))
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(nx-1,j,k))

          ! conservative variables
          rho_n(nx,j,k) =ro_wall
          rhoe_n(nx,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! Set increments to zero
    ! ======================
    do k=ndz_e,nfz_e
       Krho(nx,ndy_e:nfy_e,k)=0.0_wp
       Krhou(nx,ndy_e:nfy_e,k)=0.0_wp
       Krhov(nx,ndy_e:nfy_e,k)=0.0_wp
       Krhow(nx,ndy_e:nfy_e,k)=0.0_wp
       Krhoe(nx,ndy_e:nfy_e,k)=0.0_wp
    enddo

  end subroutine bc_wall_imax_isotherm_dpdn_c

  !==============================================================================
  module subroutine bc_wall_jmin_adiabatic_dpdn_c
  !==============================================================================
    !> Apply wall Boundary Conditions at jmin: adiabatic + dp/dn=0
    !> - curvilinear version - (all variables are imposed)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,k
    real(wp) :: ro_wall,p_wall
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
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(i,2,k)-prs(i,3,k) &
                + (prs(i+1,1,k)-prs(i-1,1,k))*gksigeta_jmin(i))
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,2,k))

          ! conservative variables
          rho_n(i,1,k) =ro_wall
          rhoe_n(i,1,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(i,2,k)-prs(i,3,k) &
                + 2.0_wp*(prs(i+1,1,k)-prs(i,1,k))*gksigeta_jmin(i))
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,2,k))

          ! conservative variables
          rho_n(i,1,k) =ro_wall
          rhoe_n(i,1,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(i,2,k)-prs(i,3,k) &
                + 2.0_wp*(prs(i,1,k)-prs(i-1,1,k))*gksigeta_jmin(i))
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,2,k))

          ! conservative variables
          rho_n(i,1,k) =ro_wall
          rhoe_n(i,1,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! Set increments to zero
    ! ======================
    do k=ndz_e,nfz_e
       Krho(ndx_e:nfx_e,1,k)=0.0_wp
       Krhou(ndx_e:nfx_e,1,k)=0.0_wp
       Krhov(ndx_e:nfx_e,1,k)=0.0_wp
       Krhow(ndx_e:nfx_e,1,k)=0.0_wp
       Krhoe(ndx_e:nfx_e,1,k)=0.0_wp
    enddo
    
  end subroutine bc_wall_jmin_adiabatic_dpdn_c

  !==============================================================================
  module subroutine bc_wall_jmin_isotherm_dpdn_c
  !==============================================================================
    !> Apply wall Boundary Conditions at jmin: isothermal + dp/dn=0
    !> - curvilinear version - (all variables are imposed)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,k
    real(wp) :: ro_wall,p_wall
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
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(i,2,k)-prs(i,3,k) &
                + (prs(i+1,1,k)-prs(i-1,1,k))*gksigeta_jmin(i))
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,2,k))

          ! conservative variables
          rho_n(i,1,k) =ro_wall
          rhoe_n(i,1,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    enddo

    ! corner with other BC at imin
    ! ----------------------------
    if (ndx_jmin==2) then
       i=1
       do k=1,nz
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(i,2,k)-prs(i,3,k) &
                + 2.0_wp*(prs(i+1,1,k)-prs(i,1,k))*gksigeta_jmin(i))
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,2,k))

          ! conservative variables
          rho_n(i,1,k) =ro_wall
          rhoe_n(i,1,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! corner with other BC at imax
    ! ----------------------------
    if (nfx_jmin==nx-1) then
       i=nx
       do k=1,nz
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(i,2,k)-prs(i,3,k) &
                + 2.0_wp*(prs(i,1,k)-prs(i-1,1,k))*gksigeta_jmin(i))
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,2,k))

          ! conservative variables
          rho_n(i,1,k) =ro_wall
          rhoe_n(i,1,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! Set increments to zero
    ! ======================
    do k=ndz_e,nfz_e
       Krho(ndx_e:nfx_e,1,k)=0.0_wp
       Krhou(ndx_e:nfx_e,1,k)=0.0_wp
       Krhov(ndx_e:nfx_e,1,k)=0.0_wp
       Krhow(ndx_e:nfx_e,1,k)=0.0_wp
       Krhoe(ndx_e:nfx_e,1,k)=0.0_wp
    enddo
    
  end subroutine bc_wall_jmin_isotherm_dpdn_c

  !==============================================================================
  module subroutine bc_wall_jmax_adiabatic_dpdn_c
  !==============================================================================
    !> Apply wall Boundary Conditions at jmax: adiabatic + dp/dn=0
    !> - curvilinear version - (all variables are imposed)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,k
    real(wp) :: ro_wall,p_wall
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
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(i,ny-1,k)-prs(i,ny-2,k) &
                + (prs(i+1,ny,k)-prs(i-1,ny,k))*gksigeta_jmax(i))
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,ny-1,k))

          ! conservative variables
          rho_n(i,ny,k) =ro_wall
          rhoe_n(i,ny,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(i,ny-1,k)-prs(i,ny-2,k) &
                + 2.0_wp*(prs(i+1,ny,k)-prs(i,ny,k))*gksigeta_jmax(i))
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,ny-1,k))

          ! conservative variables
          rho_n(i,ny,k) =ro_wall
          rhoe_n(i,ny,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(i,ny-1,k)-prs(i,ny-2,k) &
                + 2.0_wp*(prs(i,ny,k)-prs(i-1,ny,k))*gksigeta_jmax(i))
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,ny-1,k))

          ! conservative variables
          rho_n(i,ny,k) =ro_wall
          rhoe_n(i,ny,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! Set increments to zero
    ! ======================
    do k=ndz_e,nfz_e
       Krho(ndx_e:nfx_e,ny,k)=0.0_wp
       Krhou(ndx_e:nfx_e,ny,k)=0.0_wp
       Krhov(ndx_e:nfx_e,ny,k)=0.0_wp
       Krhow(ndx_e:nfx_e,ny,k)=0.0_wp
       Krhoe(ndx_e:nfx_e,ny,k)=0.0_wp
    enddo
    
  end subroutine bc_wall_jmax_adiabatic_dpdn_c

  !==============================================================================
  module subroutine bc_wall_jmax_isotherm_dpdn_c
  !==============================================================================
    !> Apply wall Boundary Conditions at jmax: isothermal + dp/dn=0
    !> - curvilinear version - (all variables are imposed)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,k
    real(wp) :: ro_wall,p_wall
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
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(i,ny-1,k)-prs(i,ny-2,k) &
                + (prs(i+1,ny,k)-prs(i-1,ny,k))*gksigeta_jmax(i))
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,ny-1,k))

          ! conservative variables
          rho_n(i,ny,k) =ro_wall
          rhoe_n(i,ny,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    enddo

    ! corner with other BC at imin
    ! ----------------------------
    if (ndx_jmax==2) then
       i=1
       do k=1,nz
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(i,ny-1,k)-prs(i,ny-2,k) &
                + 2.0_wp*(prs(i+1,ny,k)-prs(i,ny,k))*gksigeta_jmax(i))
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,ny-1,k))

          ! conservative variables
          rho_n(i,ny,k) =ro_wall
          rhoe_n(i,ny,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! corner with other BC at imax
    ! ----------------------------
    if (nfx_jmax==nx-1) then
       i=nx
       do k=1,nz
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(i,ny-1,k)-prs(i,ny-2,k) &
                + 2.0_wp*(prs(i,ny,k)-prs(i-1,ny,k))*gksigeta_jmax(i))
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,ny-1,k))

          ! conservative variables
          rho_n(i,ny,k) =ro_wall
          rhoe_n(i,ny,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! Set increments to zero
    ! ======================
    do k=ndz_e,nfz_e
       Krho(ndx_e:nfx_e,ny,k)=0.0_wp
       Krhou(ndx_e:nfx_e,ny,k)=0.0_wp
       Krhov(ndx_e:nfx_e,ny,k)=0.0_wp
       Krhow(ndx_e:nfx_e,ny,k)=0.0_wp
       Krhoe(ndx_e:nfx_e,ny,k)=0.0_wp
    enddo
    
  end subroutine bc_wall_jmax_isotherm_dpdn_c

end submodule smod_bc_wall_c
