!================================================================================
submodule (mod_bc_wall) smod_bc_wall
!================================================================================
  !> Submodule to apply wall Boundary Conditions on conservative variables
  !> - Cartesian version - (all variables are imposed)
!================================================================================

contains

  !==============================================================================
  module subroutine bc_wall_imin_adiabatic_dpdn
  !==============================================================================
    !> Apply wall Boundary Conditions at imin: adiabatic + dp/dn=0
    !> - Cartesian version - (all variables are imposed)
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
          
          ! dT/dn=0
          T_wall=cextrp2(1,1)*Tmp(2,j,k)-cextrp3(1,1)*Tmp(3,j,k)
          ! dp/dn=0
          p_wall=cextrp2(1,1)*prs(2,j,k)-cextrp3(1,1)*prs(3,j,k)
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(2,j,k))
          
          ! conservative variables
          rho_n(1,j,k) =ro_wall
          rhoe_n(1,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    enddo

    ! Set increments to zero
    ! ======================
    do k=ndz_e,nfz_e
       Krho(1,ndy_e:nfy_e,k)=0.0_wp
       Krhou(1,ndy_e:nfy_e,k)=0.0_wp
       Krhov(1,ndy_e:nfy_e,k)=0.0_wp
       Krhow(1,ndy_e:nfy_e,k)=0.0_wp
       Krhoe(1,ndy_e:nfy_e,k)=0.0_wp
    enddo
    
  end subroutine bc_wall_imin_adiabatic_dpdn

  !==============================================================================
  module subroutine bc_wall_imin_isotherm_dpdn
  !==============================================================================
    !> Apply wall Boundary Conditions at imin: isothermal + dp/dn=0
    !> - Cartesian version - (all variables are imposed)
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
          
          ! dp/dn=0
          p_wall=cextrp2(1,1)*prs(2,j,k)-cextrp3(1,1)*prs(3,j,k)
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(2,j,k))

          ! conservative variables
          rho_n(1,j,k) =ro_wall
          rhoe_n(1,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    enddo

    ! Set increments to zero
    ! ======================
    do k=ndz_e,nfz_e
       Krho(1,ndy_e:nfy_e,k)=0.0_wp
       Krhou(1,ndy_e:nfy_e,k)=0.0_wp
       Krhov(1,ndy_e:nfy_e,k)=0.0_wp
       Krhow(1,ndy_e:nfy_e,k)=0.0_wp
       Krhoe(1,ndy_e:nfy_e,k)=0.0_wp
    enddo
    
  end subroutine bc_wall_imin_isotherm_dpdn

  !==============================================================================
  module subroutine bc_wall_imax_adiabatic_dpdn
  !==============================================================================
    !> Apply wall Boundary Conditions at imax: adiabatic + dp/dn=0
    !> - Cartesian version - (all variables are imposed)
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
          
          ! dT/dn=0
          T_wall=cextrp2(1,2)*Tmp(nx-1,j,k)-cextrp3(1,2)*Tmp(nx-2,j,k)
          ! dp/dn=0
          p_wall=cextrp2(1,2)*prs(nx-1,j,k)-cextrp3(1,2)*prs(nx-2,j,k)
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(nx-1,j,k))

          ! conservative variables
          rho_n(nx,j,k) =ro_wall
          rhoe_n(nx,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    enddo

    ! Set increments to zero
    ! ======================
    do k=ndz_e,nfz_e
       Krho(nx,ndy_e:nfy_e,k)=0.0_wp
       Krhou(nx,ndy_e:nfy_e,k)=0.0_wp
       Krhov(nx,ndy_e:nfy_e,k)=0.0_wp
       Krhow(nx,ndy_e:nfy_e,k)=0.0_wp
       Krhoe(nx,ndy_e:nfy_e,k)=0.0_wp
    enddo
    
  end subroutine bc_wall_imax_adiabatic_dpdn

  !==============================================================================
  module subroutine bc_wall_imax_isotherm_dpdn
  !==============================================================================
    !> Apply wall Boundary Conditions at imax: isothermal + dp/dn=0
    !> - Cartesian version - (all variables are imposed)
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
          
          ! dp/dn=0
          p_wall=cextrp2(1,2)*prs(nx-1,j,k)-cextrp3(1,2)*prs(nx-2,j,k)
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(nx-1,j,k))

          ! conservative variables
          rho_n(nx,j,k) =ro_wall
          rhoe_n(nx,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    enddo

    ! Set increments to zero
    ! ======================
    do k=ndz_e,nfz_e
       Krho(nx,ndy_e:nfy_e,k)=0.0_wp
       Krhou(nx,ndy_e:nfy_e,k)=0.0_wp
       Krhov(nx,ndy_e:nfy_e,k)=0.0_wp
       Krhow(nx,ndy_e:nfy_e,k)=0.0_wp
       Krhoe(nx,ndy_e:nfy_e,k)=0.0_wp
    enddo
    
  end subroutine bc_wall_imax_isotherm_dpdn

  !==============================================================================
  module subroutine bc_wall_jmin_adiabatic_dpdn
  !==============================================================================
    !> Apply wall Boundary Conditions at jmin: adiabatic + dp/dn=0
    !> - Cartesian version - (all variables are imposed)
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
          
          ! dT/dn=0
          T_wall=cextrp2(2,1)*Tmp(i,2,k)-cextrp3(2,1)*Tmp(i,3,k)
          ! dp/dn=0
          p_wall=cextrp2(2,1)*prs(i,2,k)-cextrp3(2,1)*prs(i,3,k)
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,2,k))

          ! conservative variables
          rho_n(i,1,k) =ro_wall
          rhoe_n(i,1,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    enddo

    ! Set increments to zero
    ! ======================
    do k=ndz_e,nfz_e
       Krho(ndx_e:nfx_e,1,k)=0.0_wp
       Krhou(ndx_e:nfx_e,1,k)=0.0_wp
       Krhov(ndx_e:nfx_e,1,k)=0.0_wp
       Krhow(ndx_e:nfx_e,1,k)=0.0_wp
       Krhoe(ndx_e:nfx_e,1,k)=0.0_wp
    enddo
    
  end subroutine bc_wall_jmin_adiabatic_dpdn

  !==============================================================================
  module subroutine bc_wall_jmin_isotherm_dpdn
  !==============================================================================
    !> Apply wall Boundary Conditions at jmin: isothermal + dp/dn=0
    !> - Cartesian version - (all variables are imposed)
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

          ! dp/dn=0
          p_wall=cextrp2(2,1)*prs(i,2,k)-cextrp3(2,1)*prs(i,3,k)
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,2,k))

          ! conservative variables
          rho_n(i,1,k) =ro_wall
          rhoe_n(i,1,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    enddo

    ! Set increments to zero
    ! ======================
    do k=ndz_e,nfz_e
       Krho(ndx_e:nfx_e,1,k)=0.0_wp
       Krhou(ndx_e:nfx_e,1,k)=0.0_wp
       Krhov(ndx_e:nfx_e,1,k)=0.0_wp
       Krhow(ndx_e:nfx_e,1,k)=0.0_wp
       Krhoe(ndx_e:nfx_e,1,k)=0.0_wp
    enddo
    
  end subroutine bc_wall_jmin_isotherm_dpdn

  !==============================================================================
  module subroutine bc_wall_jmax_adiabatic_dpdn
  !==============================================================================
    !> Apply wall Boundary Conditions at jmax: adiabatic + dp/dn=0
    !> - Cartesian version - (all variables are imposed)
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

          ! dT/dn=0
          T_wall=cextrp2(2,2)*Tmp(i,ny-1,k)-cextrp3(2,2)*Tmp(i,ny-2,k)
          ! dp/dn=0
          p_wall=cextrp2(2,2)*prs(i,ny-1,k)-cextrp3(2,2)*prs(i,ny-2,k)
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,ny-1,k))

          ! conservative variables
          rho_n(i,ny,k) =ro_wall
          rhoe_n(i,ny,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    enddo

    ! Set increments to zero
    ! ======================
    do k=ndz_e,nfz_e
       Krho(ndx_e:nfx_e,ny,k)=0.0_wp
       Krhou(ndx_e:nfx_e,ny,k)=0.0_wp
       Krhov(ndx_e:nfx_e,ny,k)=0.0_wp
       Krhow(ndx_e:nfx_e,ny,k)=0.0_wp
       Krhoe(ndx_e:nfx_e,ny,k)=0.0_wp
    enddo
    
  end subroutine bc_wall_jmax_adiabatic_dpdn

  !==============================================================================
  module subroutine bc_wall_jmax_isotherm_dpdn
  !==============================================================================
    !> Apply wall Boundary Conditions at jmax: isothermal + dp/dn=0
    !> - Cartesian version - (all variables are imposed)
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
          
          ! dp/dn=0
          p_wall=cextrp2(2,2)*prs(i,ny-1,k)-cextrp3(2,2)*prs(i,ny-2,k)
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,ny-1,k))

          ! conservative variables
          rho_n(i,ny,k) =ro_wall
          rhoe_n(i,ny,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    enddo

    ! Set increments to zero
    ! ======================
    do k=ndz_e,nfz_e
       Krho(ndx_e:nfx_e,ny,k)=0.0_wp
       Krhou(ndx_e:nfx_e,ny,k)=0.0_wp
       Krhov(ndx_e:nfx_e,ny,k)=0.0_wp
       Krhow(ndx_e:nfx_e,ny,k)=0.0_wp
       Krhoe(ndx_e:nfx_e,ny,k)=0.0_wp
    enddo
    
  end subroutine bc_wall_jmax_isotherm_dpdn

  !==============================================================================
  module subroutine bc_wall_kmin_adiabatic_dpdn
  !==============================================================================
    !> Apply wall Boundary Conditions at kmin: adiabatic + dp/dn=0
    !> - Cartesian version - (all variables are imposed)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j
    real(wp) :: ro_wall,p_wall
    ! ---------------------------------------------------------------------------

    ! Impose wall conditions on conservative variables
    ! ================================================
    do j=1,ny
       do i=1,nx
          ! no slip conditions
          rhou_n(i,j,1)=0.0_wp
          rhov_n(i,j,1)=0.0_wp
          rhow_n(i,j,1)=0.0_wp

          ! dT/dn=0
          T_wall=cextrp2(3,1)*Tmp(i,j,2)-cextrp3(3,1)*Tmp(i,j,3)
          ! dp/dn=0
          p_wall=cextrp2(3,1)*prs(i,j,2)-cextrp3(3,1)*prs(i,j,3)
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,2))

          ! conservative variables
          rho_n(i,j,1) =ro_wall
          rhoe_n(i,j,1)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    enddo

    ! Set increments to zero
    ! ======================
    do i=ndx_e,nfx_e
       Krho(i,ndy_e:nfy_e,1)=0.0_wp
       Krhou(i,ndy_e:nfy_e,1)=0.0_wp
       Krhov(i,ndy_e:nfy_e,1)=0.0_wp
       Krhow(i,ndy_e:nfy_e,1)=0.0_wp
       Krhoe(i,ndy_e:nfy_e,1)=0.0_wp
    enddo
    
  end subroutine bc_wall_kmin_adiabatic_dpdn

  !==============================================================================
  module subroutine bc_wall_kmin_isotherm_dpdn
  !==============================================================================
    !> Apply wall Boundary Conditions at kmin: isothermal + dp/dn=0
    !> - Cartesian version - (all variables are imposed)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j
    real(wp) :: ro_wall,p_wall
    ! ---------------------------------------------------------------------------

    ! Impose wall conditions on conservative variables
    ! ================================================
    do j=1,ny
       do i=1,nx
          ! no slip conditions
          rhou_n(i,j,1)=0.0_wp
          rhov_n(i,j,1)=0.0_wp
          rhow_n(i,j,1)=0.0_wp

          ! dp/dn=0
          p_wall=cextrp2(3,1)*prs(i,j,2)-cextrp3(3,1)*prs(i,j,3)
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,2))

          ! conservative variables
          rho_n(i,j,1) =ro_wall
          rhoe_n(i,j,1)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    enddo

    ! Set increments to zero
    ! ======================
    do i=ndx_e,nfx_e
       Krho(i,ndy_e:nfy_e,1)=0.0_wp
       Krhou(i,ndy_e:nfy_e,1)=0.0_wp
       Krhov(i,ndy_e:nfy_e,1)=0.0_wp
       Krhow(i,ndy_e:nfy_e,1)=0.0_wp
       Krhoe(i,ndy_e:nfy_e,1)=0.0_wp
    enddo
    
  end subroutine bc_wall_kmin_isotherm_dpdn

  !==============================================================================
  module subroutine bc_wall_kmax_adiabatic_dpdn
  !==============================================================================
    !> Apply wall Boundary Conditions at kmax: adiabatic + dp/dn=0
    !> - Cartesian version - (all variables are imposed)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j
    real(wp) :: ro_wall,p_wall
    ! ---------------------------------------------------------------------------

    ! Impose wall conditions on conservative variables
    ! ================================================
    do j=1,ny
       do i=1,nx
          ! no slip conditions
          rhou_n(i,j,nz)=0.0_wp
          rhov_n(i,j,nz)=0.0_wp
          rhow_n(i,j,nz)=0.0_wp

          ! dT/dn=0
          T_wall=cextrp2(3,2)*Tmp(i,j,nz-1)-cextrp3(3,2)*Tmp(i,j,nz-2)
          ! dp/dn=0
          p_wall=cextrp2(3,2)*prs(i,j,nz-1)-cextrp3(3,2)*prs(i,j,nz-2)
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,nz-1))

          ! conservative variables
          rho_n(i,j,nz) =ro_wall
          rhoe_n(i,j,nz)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    enddo

    ! Set increments to zero
    ! ======================
    do i=ndx_e,nfx_e
       Krho(i,ndy_e:nfy_e,nz)=0.0_wp
       Krhou(i,ndy_e:nfy_e,nz)=0.0_wp
       Krhov(i,ndy_e:nfy_e,nz)=0.0_wp
       Krhow(i,ndy_e:nfy_e,nz)=0.0_wp
       Krhoe(i,ndy_e:nfy_e,nz)=0.0_wp
    enddo
    
  end subroutine bc_wall_kmax_adiabatic_dpdn

  !==============================================================================
  module subroutine bc_wall_kmax_isotherm_dpdn
  !==============================================================================
    !> Apply wall Boundary Conditions at jmax: isothermal + dp/dn=0
    !> - Cartesian version - (all variables are imposed)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j
    real(wp) :: ro_wall,p_wall
    ! ---------------------------------------------------------------------------

    ! Impose wall conditions on conservative variables
    ! ================================================
    do j=1,ny
       do i=1,nx
          ! no slip conditions
          rhou_n(i,j,nz)=0.0_wp
          rhov_n(i,j,nz)=0.0_wp
          rhow_n(i,j,nz)=0.0_wp

          ! dp/dn=0
          p_wall=cextrp2(3,2)*prs(i,j,nz-1)-cextrp3(3,2)*prs(i,j,nz-2)
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,nz-1))

          ! conservative variables
          rho_n(i,j,nz) =ro_wall
          rhoe_n(i,j,nz)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    enddo

    ! Set increments to zero
    ! ======================
    do i=ndx_e,nfx_e
       Krho(i,ndy_e:nfy_e,nz)=0.0_wp
       Krhou(i,ndy_e:nfy_e,nz)=0.0_wp
       Krhov(i,ndy_e:nfy_e,nz)=0.0_wp
       Krhow(i,ndy_e:nfy_e,nz)=0.0_wp
       Krhoe(i,ndy_e:nfy_e,nz)=0.0_wp
    enddo
    
  end subroutine bc_wall_kmax_isotherm_dpdn

end submodule smod_bc_wall
