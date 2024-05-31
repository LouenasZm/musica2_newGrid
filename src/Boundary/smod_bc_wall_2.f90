!================================================================================
submodule (mod_bc_wall) smod_bc_wall_2
!================================================================================
  !> Submodule to apply wall Boundary Conditions on conservative variables
  !> - Cartesian version - (advancement of rho on boundary)
!================================================================================

contains

  !==============================================================================
  module subroutine bc_wall_imin_slip
  !==============================================================================
    !> Apply wall Boundary Conditions at imin: slip wall
    !> - Cartesian version - (advancement of rho on boundary)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: j,k
    ! ---------------------------------------------------------------------------

    ! Impose wall conditions on conservative variables
    ! ================================================
    do k=1,nz
       do j=1,ny
          ! slip condition: nullity of normal velocity component
          rhou_n(1,j,k)=0.0_wp
       enddo
    enddo

  end subroutine bc_wall_imin_slip

  !==============================================================================
  module subroutine bc_wall_imin_adiabatic
  !==============================================================================
    !> Apply wall Boundary Conditions at imin: no slip + adiabatic
    !> - Cartesian version - (advancement of rho on boundary)
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
          
          ! dT/dn=0
          T_wall=cextrp2(1,1)*Tmp(2,j,k)-cextrp3(1,1)*Tmp(3,j,k)
          rhoe_n(1,j,k)=rho_n(1,j,k)*ecalc_tro(T_wall,rho_n(1,j,k))
       enddo
    enddo

  end subroutine bc_wall_imin_adiabatic

  !==============================================================================
  module subroutine bc_wall_imin_isotherm
  !==============================================================================
    !> Apply wall Boundary Conditions at imin: no slip + isothermal
    !> - Cartesian version - (advancement of rho on boundary)
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
          
          ! T=cste
          rhoe_n(1,j,k)=rho_n(1,j,k)*ecalc_tro(T_wall,rho_n(1,j,k))
       enddo
    enddo

  end subroutine bc_wall_imin_isotherm

  !==============================================================================
  module subroutine bc_wall_imax_slip
  !==============================================================================
    !> Apply wall Boundary Conditions at imax: slip wall
    !> - Cartesian version - (advancement of rho on boundary)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: j,k
    ! ---------------------------------------------------------------------------

    ! Impose wall conditions on conservative variables
    ! ================================================
    do k=1,nz
       do j=1,ny
          ! slip condition: nullity of normal velocity component
          rhou_n(nx,j,k)=0.0_wp
       enddo
    enddo

  end subroutine bc_wall_imax_slip

  !==============================================================================
  module subroutine bc_wall_imax_adiabatic
  !==============================================================================
    !> Apply wall Boundary Conditions at imax: no slip + adiabatic
    !> - Cartesian version - (advancement of rho on boundary)
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
          
          ! dT/dn=0
          T_wall=cextrp2(1,2)*Tmp(nx-1,j,k)-cextrp3(1,2)*Tmp(nx-2,j,k)
          rhoe_n(nx,j,k)=rho_n(nx,j,k)*ecalc_tro(T_wall,rho_n(nx,j,k))
       enddo
    enddo

  end subroutine bc_wall_imax_adiabatic

  !==============================================================================
  module subroutine bc_wall_imax_isotherm
  !==============================================================================
    !> Apply wall Boundary Conditions at imax: no slip + isothermal
    !> - Cartesian version - (advancement of rho on boundary)
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
          
          ! T=cste
          rhoe_n(nx,j,k)=rho_n(nx,j,k)*ecalc_tro(T_wall,rho_n(nx,j,k))
       enddo
    enddo

  end subroutine bc_wall_imax_isotherm

  !==============================================================================
  module subroutine bc_wall_jmin_slip
  !==============================================================================
    !> Apply wall Boundary Conditions at jmin: slip wall
    !> - Cartesian version - (advancement of rho on boundary)
  !==============================================================================
    use mod_constant
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,k
    ! ---------------------------------------------------------------------------

    ! Impose wall conditions on conservative variables
    ! ================================================
    do k=1,nz
       do i=1,nx
          ! slip condition: nullity of normal velocity component
          rhov_n(i,1,k)=0.0_wp
       enddo
    enddo
   
  end subroutine bc_wall_jmin_slip

  !==============================================================================
  module subroutine bc_wall_jmin_adiabatic
  !==============================================================================
    !> Apply wall Boundary Conditions at jmin: no slip + adiabatic
    !> - Cartesian version - (advancement of rho on boundary)
  !==============================================================================
    use mod_constant
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
          
          ! dT/dn=0
          T_wall=cextrp2(2,1)*Tmp(i,2,k)-cextrp3(2,1)*Tmp(i,3,k)
          rhoe_n(i,1,k)=rho_n(i,1,k)*ecalc_tro(T_wall,rho_n(i,1,k))
       enddo
    enddo
   
  end subroutine bc_wall_jmin_adiabatic

  !==============================================================================
  module subroutine bc_wall_jmin_isotherm
  !==============================================================================
    !> Apply wall Boundary Conditions at jmin: no slip + isothermal
    !> - Cartesian version - (advancement of rho on boundary)
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
          
          ! T=cste
          rhoe_n(i,1,k)=rho_n(i,1,k)*ecalc_tro(T_wall,rho_n(i,1,k))
       enddo
    enddo
   
  end subroutine bc_wall_jmin_isotherm

  !==============================================================================
  module subroutine bc_wall_jmax_slip
  !==============================================================================
    !> Apply wall Boundary Conditions at jmax: slip wall
    !> - Cartesian version - (advancement of rho on boundary)
  !==============================================================================
    use mod_constant
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,k
    ! ---------------------------------------------------------------------------

    ! Impose wall conditions on conservative variables
    ! ================================================
    do k=1,nz
       do i=1,nx
          ! slip condition: nullity of normal velocity component
          rhov_n(i,ny,k)=0.0_wp
       enddo
    enddo
   
  end subroutine bc_wall_jmax_slip

  !==============================================================================
  module subroutine bc_wall_jmax_adiabatic
  !==============================================================================
    !> Apply wall Boundary Conditions at jmax: no slip + adiabatic
    !> - Cartesian version - (advancement of rho on boundary)
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
          
          ! dT/dn=0
          T_wall=cextrp2(2,2)*Tmp(i,ny-1,k)-cextrp3(2,2)*Tmp(i,ny-2,k)
          rhoe_n(i,ny,k)=rho_n(i,ny,k)*ecalc_tro(T_wall,rho_n(i,ny,k))
       enddo
    enddo

  end subroutine bc_wall_jmax_adiabatic

  !==============================================================================
  module subroutine bc_wall_jmax_isotherm
  !==============================================================================
    !> Apply wall Boundary Conditions at jmax: no slip + isothermal
    !> - Cartesian version - (advancement of rho on boundary)
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
          
          ! T=cste
          rhoe_n(i,ny,k)=rho_n(i,ny,k)*ecalc_tro(T_wall,rho_n(i,ny,k))
       enddo
    enddo

  end subroutine bc_wall_jmax_isotherm

  !==============================================================================
  module subroutine bc_wall_kmin_slip
  !==============================================================================
    !> Apply wall Boundary Conditions at kmin: slip wall
    !> - Cartesian version - (advancement of rho on boundary)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j
    ! ---------------------------------------------------------------------------

    ! Impose wall conditions on conservative variables
    ! ================================================
    do j=1,ny
       do i=1,nx
          ! slip condition: nullity of normal velocity component
          rhow_n(i,j,1)=0.0_wp
       enddo
    enddo

  end subroutine bc_wall_kmin_slip

  !==============================================================================
  module subroutine bc_wall_kmin_adiabatic
  !==============================================================================
    !> Apply wall Boundary Conditions at kmin: no slip + adiabatic
    !> - Cartesian version - (advancement of rho on boundary)
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
          
          ! dT/dn=0
          T_wall=cextrp2(3,1)*Tmp(i,j,2)-cextrp3(3,1)*Tmp(i,j,3)
          rhoe_n(i,j,1)=rho_n(i,j,1)*ecalc_tro(T_wall,rho_n(i,j,1))
       enddo
    enddo

  end subroutine bc_wall_kmin_adiabatic

  !==============================================================================
  module subroutine bc_wall_kmin_isotherm
  !==============================================================================
    !> Apply wall Boundary Conditions at kmin: no slip + isothermal
    !> - Cartesian version - (advancement of rho on boundary)
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
          
          ! T=cste
          rhoe_n(i,j,1)=rho_n(i,j,1)*ecalc_tro(T_wall,rho_n(i,j,1))
       enddo
    enddo

  end subroutine bc_wall_kmin_isotherm

  !==============================================================================
  module subroutine bc_wall_kmax_slip
  !==============================================================================
    !> Apply wall Boundary Conditions at kmax: slip wall
    !> - Cartesian version - (advancement of rho on boundary)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j
    ! ---------------------------------------------------------------------------

    ! Impose wall conditions on conservative variables
    ! ================================================
    do j=1,ny
       do i=1,nx
          ! slip condition: nullity of normal velocity component
          rhow_n(i,j,nz)=0.0_wp
       enddo
    enddo

  end subroutine bc_wall_kmax_slip

  !==============================================================================
  module subroutine bc_wall_kmax_adiabatic
  !==============================================================================
    !> Apply wall Boundary Conditions at kmax: no slip + adiabatic
    !> - Cartesian version - (advancement of rho on boundary)
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
          
          ! dT/dn=0
          T_wall=cextrp2(3,2)*Tmp(i,j,nz-1)-cextrp3(3,2)*Tmp(i,j,nz-2)
          rhoe_n(i,j,nz)=rho_n(i,j,nz)*ecalc_tro(T_wall,rho_n(i,j,nz))
       enddo
    enddo

  end subroutine bc_wall_kmax_adiabatic

  !==============================================================================
  module subroutine bc_wall_kmax_isotherm
  !==============================================================================
    !> Apply wall Boundary Conditions at jmax: no slip + isothermal
    !> - Cartesian version - (advancement of rho on boundary)
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
          
          ! T=cste
          rhoe_n(i,j,nz)=rho_n(i,j,nz)*ecalc_tro(T_wall,rho_n(i,j,nz))
       enddo
    enddo

  end subroutine bc_wall_kmax_isotherm

end submodule smod_bc_wall_2
