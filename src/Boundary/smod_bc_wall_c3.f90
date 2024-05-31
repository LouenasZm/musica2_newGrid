!================================================================================
submodule (mod_bc_wall) smod_bc_wall_c3
!================================================================================
  !> Submodule to apply wall Boundary Conditions - 3D curvilinear version -
  !> - 3D curvilinear version - (all variables are imposed)
!================================================================================
contains

  !==============================================================================
  module subroutine bc_wall_imin_adiabatic_dpdn_c3
  !==============================================================================
    !> Apply wall Boundary Conditions at imin: adiabatic + dp/dn=0
    !> - 3D curvilinear version - (all variables are imposed)
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
          
    do k=ndz_imin,nfz_imin
       do j=ndy_imin,nfy_imin
          ! dT/dn=0
          T_wall= onethird*( 4.0_wp*Tmp(2,j,k)-Tmp(3,j,k)        &
                +(Tmp(1,j+1,k)-Tmp(1,j-1,k))*getagksi3_imin(j,k) &
                +(Tmp(1,j,k+1)-Tmp(1,j,k-1))*gphigksi3_imin(j,k) )
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(2,j,k)-prs(3,j,k)        &
                +(prs(1,j+1,k)-prs(1,j-1,k))*getagksi3_imin(j,k) &
                +(prs(1,j,k+1)-prs(1,j,k-1))*gphigksi3_imin(j,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(2,j,k) )

          ! conservative variables
          rho_n(1,j,k) =ro_wall
          rhoe_n(1,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(2,j,k)-prs(3,j,k)        &
                +2.0_wp*(prs(1,j+1,k)-prs(1,j,k))*getagksi3_imin(j,k) &
                +(prs(1,j,k+1)-prs(1,j,k-1))*gphigksi3_imin(j,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(2,j,k) )

          ! conservative variables
          rho_n(1,j,k) =ro_wall
          rhoe_n(1,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(2,j,k)-prs(3,j,k)        &
                +2.0_wp*(prs(1,j,k)-prs(1,j-1,k))*getagksi3_imin(j,k) &
                +(prs(1,j,k+1)-prs(1,j,k-1))*gphigksi3_imin(j,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(2,j,k) )

          ! conservative variables
          rho_n(1,j,k) =ro_wall
          rhoe_n(1,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(2,j,k)-prs(3,j,k)        &
                +(prs(1,j+1,k)-prs(1,j-1,k))*getagksi3_imin(j,k) &
                +2.0_wp*(prs(1,j,k+1)-prs(1,j,k))*gphigksi3_imin(j,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(2,j,k) )

          ! conservative variables
          rho_n(1,j,k) =ro_wall
          rhoe_n(1,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(2,j,k)-prs(3,j,k)        &
                +(prs(1,j+1,k)-prs(1,j-1,k))*getagksi3_imin(j,k) &
                +2.0_wp*(prs(1,j,k)-prs(1,j,k-1))*gphigksi3_imin(j,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(2,j,k) )

          ! conservative variables
          rho_n(1,j,k) =ro_wall
          rhoe_n(1,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(2,j,k)-prs(3,j,k)        &
             +2.0_wp*(prs(1,j+1,k)-prs(1,j,k))*getagksi3_imin(j,k) &
             +2.0_wp*(prs(1,j,k+1)-prs(1,j,k))*gphigksi3_imin(j,k) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(2,j,k) )

       ! conservative variables
       rho_n(1,j,k) =ro_wall
       rhoe_n(1,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(2,j,k)-prs(3,j,k)        &
             +2.0_wp*(prs(1,j+1,k)-prs(1,j,k))*getagksi3_imin(j,k) &
             +2.0_wp*(prs(1,j,k)-prs(1,j,k-1))*gphigksi3_imin(j,k) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(2,j,k) )

       ! conservative variables
       rho_n(1,j,k) =ro_wall
       rhoe_n(1,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(2,j,k)-prs(3,j,k)        &
             +2.0_wp*(prs(1,j,k)-prs(1,j-1,k))*getagksi3_imin(j,k) &
             +2.0_wp*(prs(1,j,k+1)-prs(1,j,k))*gphigksi3_imin(j,k) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(2,j,k) )

       ! conservative variables
       rho_n(1,j,k) =ro_wall
       rhoe_n(1,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(2,j,k)-prs(3,j,k)        &
             +2.0_wp*(prs(1,j,k)-prs(1,j-1,k))*getagksi3_imin(j,k) &
             +2.0_wp*(prs(1,j,k)-prs(1,j,k-1))*gphigksi3_imin(j,k) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(2,j,k) )

       ! conservative variables
       rho_n(1,j,k) =ro_wall
       rhoe_n(1,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
    
  end subroutine bc_wall_imin_adiabatic_dpdn_c3

  !==============================================================================
  module subroutine bc_wall_imin_isotherm_dpdn_c3
  !==============================================================================
    !> Apply wall Boundary Conditions at imin: isothermal + dp/dn=0
    !> - 3D curvilinear version - (all variables are imposed)
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
          
    do k=ndz_imin,nfz_imin
       do j=ndy_imin,nfy_imin
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(2,j,k)-prs(3,j,k)        &
                +(prs(1,j+1,k)-prs(1,j-1,k))*getagksi3_imin(j,k) &
                +(prs(1,j,k+1)-prs(1,j,k-1))*gphigksi3_imin(j,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(2,j,k) )

          ! conservative variables
          rho_n(1,j,k) =ro_wall
          rhoe_n(1,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    enddo

    ! edge with other BC at jmin
    ! --------------------------
    if (ndy_imin==2) then
       j=1
       do k=ndz_imin,nfz_imin
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(2,j,k)-prs(3,j,k)        &
                +2.0_wp*(prs(1,j+1,k)-prs(1,j,k))*getagksi3_imin(j,k) &
                +(prs(1,j,k+1)-prs(1,j,k-1))*gphigksi3_imin(j,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(2,j,k) )

          ! conservative variables
          rho_n(1,j,k) =ro_wall
          rhoe_n(1,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! edge with other BC at jmax
    ! --------------------------
    if (ndy_imin==ny-1) then
       j=ny
       do k=ndz_imin,nfz_imin
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(2,j,k)-prs(3,j,k)        &
                +2.0_wp*(prs(1,j,k)-prs(1,j-1,k))*getagksi3_imin(j,k) &
                +(prs(1,j,k+1)-prs(1,j,k-1))*gphigksi3_imin(j,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(2,j,k) )

          ! conservative variables
          rho_n(1,j,k) =ro_wall
          rhoe_n(1,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! edge with other BC at kmin
    ! --------------------------
    if (ndz_imin==2) then
       k=1
       do j=ndy_imin,nfy_imin
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(2,j,k)-prs(3,j,k)        &
                +(prs(1,j+1,k)-prs(1,j-1,k))*getagksi3_imin(j,k) &
                +2.0_wp*(prs(1,j,k+1)-prs(1,j,k))*gphigksi3_imin(j,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(2,j,k) )

          ! conservative variables
          rho_n(1,j,k) =ro_wall
          rhoe_n(1,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! edge with other BC at kmax
    ! --------------------------
    if (ndz_imin==nz-1) then
       k=nz
       do j=ndy_imin,nfy_imin
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(2,j,k)-prs(3,j,k)        &
                +(prs(1,j+1,k)-prs(1,j-1,k))*getagksi3_imin(j,k) &
                +2.0_wp*(prs(1,j,k)-prs(1,j,k-1))*gphigksi3_imin(j,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(2,j,k) )

          ! conservative variables
          rho_n(1,j,k) =ro_wall
          rhoe_n(1,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! corner with other BC at jmin-kmin
    ! ---------------------------------
    if ((ndy_imin==2).and.(ndz_imin==2)) then
       j=1
       k=1
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(2,j,k)-prs(3,j,k)        &
             +2.0_wp*(prs(1,j+1,k)-prs(1,j,k))*getagksi3_imin(j,k) &
             +2.0_wp*(prs(1,j,k+1)-prs(1,j,k))*gphigksi3_imin(j,k) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(2,j,k) )

       ! conservative variables
       rho_n(1,j,k) =ro_wall
       rhoe_n(1,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
    endif
    
    ! corner with other BC at jmin-kmax
    ! ---------------------------------
    if ((ndy_imin==2).and.(ndz_imin==nz-1)) then
       j=1
       k=nz
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(2,j,k)-prs(3,j,k)        &
             +2.0_wp*(prs(1,j+1,k)-prs(1,j,k))*getagksi3_imin(j,k) &
             +2.0_wp*(prs(1,j,k)-prs(1,j,k-1))*gphigksi3_imin(j,k) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(2,j,k) )

       ! conservative variables
       rho_n(1,j,k) =ro_wall
       rhoe_n(1,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
    endif
    
    ! corner with other BC at jmax-kmin
    ! ---------------------------------
    if ((ndy_imin==ny-1).and.(ndz_imin==2)) then
       j=ny
       k=1
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(2,j,k)-prs(3,j,k)        &
             +2.0_wp*(prs(1,j,k)-prs(1,j-1,k))*getagksi3_imin(j,k) &
             +2.0_wp*(prs(1,j,k+1)-prs(1,j,k))*gphigksi3_imin(j,k) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(2,j,k) )

       ! conservative variables
       rho_n(1,j,k) =ro_wall
       rhoe_n(1,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
    endif

    ! corner with other BC at jmax-kmax
    ! ---------------------------------
    if ((ndy_imin==ny-1).and.(ndz_imin==nz-1)) then
       j=ny
       k=nz
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(2,j,k)-prs(3,j,k)        &
             +2.0_wp*(prs(1,j,k)-prs(1,j-1,k))*getagksi3_imin(j,k) &
             +2.0_wp*(prs(1,j,k)-prs(1,j,k-1))*gphigksi3_imin(j,k) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(2,j,k) )

       ! conservative variables
       rho_n(1,j,k) =ro_wall
       rhoe_n(1,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
    
  end subroutine bc_wall_imin_isotherm_dpdn_c3

  !==============================================================================
  module subroutine bc_wall_imax_adiabatic_dpdn_c3
  !==============================================================================
    !> Apply wall Boundary Conditions at imax: adiabatic + dp/dn=0
    !> - 3D curvilinear version - (all variables are imposed)
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

    do k=ndz_imax,nfz_imax
       do j=ndy_imax,nfy_imax
          ! dT/dn=0
          T_wall= onethird*( 4.0_wp*Tmp(nx-1,j,k)-Tmp(nx-2,j,k)    &
                +(Tmp(nx,j+1,k)-Tmp(nx,j-1,k))*getagksi3_imax(j,k) &
                +(Tmp(nx,j,k+1)-Tmp(nx,j,k-1))*gphigksi3_imax(j,k) )
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(nx-1,j,k)-prs(nx-2,j,k)    &
                +(prs(nx,j+1,k)-prs(nx,j-1,k))*getagksi3_imax(j,k) &
                +(prs(nx,j,k+1)-prs(nx,j,k-1))*gphigksi3_imax(j,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(nx-1,j,k) )

          ! conservative variables
          rho_n(nx,j,k) =ro_wall
          rhoe_n(nx,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(nx-1,j,k)-prs(nx-2,j,k)    &
                +2.0_wp*(prs(nx,j+1,k)-prs(nx,j,k))*getagksi3_imax(j,k) &
                +(prs(nx,j,k+1)-prs(nx,j,k-1))*gphigksi3_imax(j,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(nx-1,j,k) )

          ! conservative variables
          rho_n(nx,j,k) =ro_wall
          rhoe_n(nx,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(nx-1,j,k)-prs(nx-2,j,k)    &
                +2.0_wp*(prs(nx,j,k)-prs(nx,j-1,k))*getagksi3_imax(j,k) &
                +(prs(nx,j,k+1)-prs(nx,j,k-1))*gphigksi3_imax(j,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(nx-1,j,k) )

          ! conservative variables
          rho_n(nx,j,k) =ro_wall
          rhoe_n(nx,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(nx-1,j,k)-prs(nx-2,j,k)    &
                +(prs(nx,j+1,k)-prs(nx,j-1,k))*getagksi3_imax(j,k) &
                +2.0_wp*(prs(nx,j,k+1)-prs(nx,j,k))*gphigksi3_imax(j,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(nx-1,j,k) )

          ! conservative variables
          rho_n(nx,j,k) =ro_wall
          rhoe_n(nx,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(nx-1,j,k)-prs(nx-2,j,k)    &
                +(prs(nx,j+1,k)-prs(nx,j-1,k))*getagksi3_imax(j,k) &
                +2.0_wp*(prs(nx,j,k)-prs(nx,j,k-1))*gphigksi3_imax(j,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(nx-1,j,k) )

          ! conservative variables
          rho_n(nx,j,k) =ro_wall
          rhoe_n(nx,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(nx-1,j,k)-prs(nx-2,j,k)    &
             +2.0_wp*(prs(nx,j+1,k)-prs(nx,j,k))*getagksi3_imax(j,k) &
             +2.0_wp*(prs(nx,j,k+1)-prs(nx,j,k))*gphigksi3_imax(j,k) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(nx-1,j,k) )

       ! conservative variables
       rho_n(nx,j,k) =ro_wall
       rhoe_n(nx,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(nx-1,j,k)-prs(nx-2,j,k)    &
             +2.0_wp*(prs(nx,j+1,k)-prs(nx,j,k))*getagksi3_imax(j,k) &
             +2.0_wp*(prs(nx,j,k)-prs(nx,j,k-1))*gphigksi3_imax(j,k) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(nx-1,j,k) )

       ! conservative variables
       rho_n(nx,j,k) =ro_wall
       rhoe_n(nx,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(nx-1,j,k)-prs(nx-2,j,k)    &
             +2.0_wp*(prs(nx,j,k)-prs(nx,j-1,k))*getagksi3_imax(j,k) &
             +2.0_wp*(prs(nx,j,k+1)-prs(nx,j,k))*gphigksi3_imax(j,k) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(nx-1,j,k) )

       ! conservative variables
       rho_n(nx,j,k) =ro_wall
       rhoe_n(nx,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(nx-1,j,k)-prs(nx-2,j,k)    &
             +2.0_wp*(prs(nx,j,k)-prs(nx,j-1,k))*getagksi3_imax(j,k) &
             +2.0_wp*(prs(nx,j,k)-prs(nx,j,k-1))*gphigksi3_imax(j,k) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(nx-1,j,k) )

       ! conservative variables
       rho_n(nx,j,k) =ro_wall
       rhoe_n(nx,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
    
  end subroutine bc_wall_imax_adiabatic_dpdn_c3

  !==============================================================================
  module subroutine bc_wall_imax_isotherm_dpdn_c3
  !==============================================================================
    !> Apply wall Boundary Conditions at imax: isothermal + dp/dn=0
    !> - 3D curvilinear version - (all variables are imposed)
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

    do k=ndz_imax,nfz_imax
       do j=ndy_imax,nfy_imax
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(nx-1,j,k)-prs(nx-2,j,k)    &
                +(prs(nx,j+1,k)-prs(nx,j-1,k))*getagksi3_imax(j,k) &
                +(prs(nx,j,k+1)-prs(nx,j,k-1))*gphigksi3_imax(j,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(nx-1,j,k) )

          ! conservative variables
          rho_n(nx,j,k) =ro_wall
          rhoe_n(nx,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    enddo

    ! edge with other BC at jmin
    ! --------------------------
    if (ndy_imax==2) then
       j=1
       do k=ndz_imax,nfz_imax
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(nx-1,j,k)-prs(nx-2,j,k)    &
                +2.0_wp*(prs(nx,j+1,k)-prs(nx,j,k))*getagksi3_imax(j,k) &
                +(prs(nx,j,k+1)-prs(nx,j,k-1))*gphigksi3_imax(j,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(nx-1,j,k) )

          ! conservative variables
          rho_n(nx,j,k) =ro_wall
          rhoe_n(nx,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! edge with other BC at jmax
    ! --------------------------
    if (ndy_imax==ny-1) then
       j=ny
       do k=ndz_imax,nfz_imax
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(nx-1,j,k)-prs(nx-2,j,k)    &
                +2.0_wp*(prs(nx,j,k)-prs(nx,j-1,k))*getagksi3_imax(j,k) &
                +(prs(nx,j,k+1)-prs(nx,j,k-1))*gphigksi3_imax(j,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(nx-1,j,k) )

          ! conservative variables
          rho_n(nx,j,k) =ro_wall
          rhoe_n(nx,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! edge with other BC at kmin
    ! --------------------------
    if (ndz_imax==2) then
       k=1
       do j=ndy_imax,nfy_imax
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(nx-1,j,k)-prs(nx-2,j,k)    &
                +(prs(nx,j+1,k)-prs(nx,j-1,k))*getagksi3_imax(j,k) &
                +2.0_wp*(prs(nx,j,k+1)-prs(nx,j,k))*gphigksi3_imax(j,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(nx-1,j,k) )

          ! conservative variables
          rho_n(nx,j,k) =ro_wall
          rhoe_n(nx,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! edge with other BC at kmax
    ! --------------------------
    if (ndz_imax==nz-1) then
       k=nz
       do j=ndy_imax,nfy_imax
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(nx-1,j,k)-prs(nx-2,j,k)    &
                +(prs(nx,j+1,k)-prs(nx,j-1,k))*getagksi3_imax(j,k) &
                +2.0_wp*(prs(nx,j,k)-prs(nx,j,k-1))*gphigksi3_imax(j,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(nx-1,j,k) )

          ! conservative variables
          rho_n(nx,j,k) =ro_wall
          rhoe_n(nx,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! corner with other BC at jmin-kmin
    ! ---------------------------------
    if ((ndy_imax==2).and.(ndz_imax==2)) then
       j=1
       k=1
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(nx-1,j,k)-prs(nx-2,j,k)    &
             +2.0_wp*(prs(nx,j+1,k)-prs(nx,j,k))*getagksi3_imax(j,k) &
             +2.0_wp*(prs(nx,j,k+1)-prs(nx,j,k))*gphigksi3_imax(j,k) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(nx-1,j,k) )

       ! conservative variables
       rho_n(nx,j,k) =ro_wall
       rhoe_n(nx,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
    endif
    
    ! corner with other BC at jmin-kmax
    ! ---------------------------------
    if ((ndy_imax==2).and.(ndz_imax==nz-1)) then
       j=1
       k=nz
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(nx-1,j,k)-prs(nx-2,j,k)    &
             +2.0_wp*(prs(nx,j+1,k)-prs(nx,j,k))*getagksi3_imax(j,k) &
             +2.0_wp*(prs(nx,j,k)-prs(nx,j,k-1))*gphigksi3_imax(j,k) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(nx-1,j,k) )

       ! conservative variables
       rho_n(nx,j,k) =ro_wall
       rhoe_n(nx,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
    endif
    
    ! corner with other BC at jmax-kmin
    ! ---------------------------------
    if ((ndy_imax==ny-1).and.(ndz_imax==2)) then
       j=ny
       k=1
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(nx-1,j,k)-prs(nx-2,j,k)    &
             +2.0_wp*(prs(nx,j,k)-prs(nx,j-1,k))*getagksi3_imax(j,k) &
             +2.0_wp*(prs(nx,j,k+1)-prs(nx,j,k))*gphigksi3_imax(j,k) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(nx-1,j,k) )

       ! conservative variables
       rho_n(nx,j,k) =ro_wall
       rhoe_n(nx,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
    endif

    ! corner with other BC at jmax-kmax
    ! ---------------------------------
    if ((ndy_imax==ny-1).and.(ndz_imax==nz-1)) then
       j=ny
       k=nz
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(nx-1,j,k)-prs(nx-2,j,k)    &
             +2.0_wp*(prs(nx,j,k)-prs(nx,j-1,k))*getagksi3_imax(j,k) &
             +2.0_wp*(prs(nx,j,k)-prs(nx,j,k-1))*gphigksi3_imax(j,k) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(nx-1,j,k) )

       ! conservative variables
       rho_n(nx,j,k) =ro_wall
       rhoe_n(nx,j,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
    
  end subroutine bc_wall_imax_isotherm_dpdn_c3

  !==============================================================================
  module subroutine bc_wall_jmin_adiabatic_dpdn_c3
  !==============================================================================
    !> Apply wall Boundary Conditions at jmin: adiabatic + dp/dn=0
    !> - 3D curvilinear version - (all variables are imposed)
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

    do k=ndz_jmin,nfz_jmin
       do i=ndx_jmin,nfx_jmin
          ! dT/dn=0
          T_wall= onethird*( 4.0_wp*Tmp(i,2,k)-Tmp(i,3,k)        &
                +(Tmp(i+1,1,k)-Tmp(i-1,1,k))*gksigeta3_jmin(i,k) &
                +(Tmp(i,1,k+1)-Tmp(i,1,k-1))*gphigeta3_jmin(i,k) )
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(i,2,k)-prs(i,3,k)         &
                +(prs(i+1,1,k)-prs(i-1,1,k))*gksigeta3_jmin(i,k) &
                +(prs(i,1,k+1)-prs(i,1,k-1))*gphigeta3_jmin(i,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,2,k) )

          ! conservative variables
          rho_n(i,1,k) =ro_wall
          rhoe_n(i,1,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(i,2,k)-prs(i,3,k)        &
                +2.0_wp*(prs(i+1,1,k)-prs(i,1,k))*gksigeta3_jmin(i,k) &
                +(prs(i,1,k+1)-prs(i,1,k-1))*gphigeta3_jmin(i,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,2,k))

          ! conservative variables
          rho_n(i,1,k) =ro_wall
          rhoe_n(i,1,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(i,2,k)-prs(i,3,k)        &
                +2.0_wp*(prs(i,1,k)-prs(i-1,1,k))*gksigeta3_jmin(i,k) &
                +(prs(i,1,k+1)-prs(i,1,k-1))*gphigeta3_jmin(i,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,2,k))

          ! conservative variables
          rho_n(i,1,k) =ro_wall
          rhoe_n(i,1,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(i,2,k)-prs(i,3,k)        &
                +(prs(i+1,1,k)-prs(i-1,1,k))*gksigeta3_jmin(i,k) &
                +2.0_wp*(prs(i,1,k+1)-prs(i,1,k))*gphigeta3_jmin(i,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,2,k))

          ! conservative variables
          rho_n(i,1,k) =ro_wall
          rhoe_n(i,1,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(i,2,k)-prs(i,3,k)        &
                +(prs(i+1,1,k)-prs(i-1,1,k))*gksigeta3_jmin(i,k) &
                +2.0_wp*(prs(i,1,k)-prs(i,1,k-1))*gphigeta3_jmin(i,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,2,k))

          ! conservative variables
          rho_n(i,1,k) =ro_wall
          rhoe_n(i,1,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(i,2,k)-prs(i,3,k)        &
             +2.0_wp*(prs(i+1,1,k)-prs(i,1,k))*gksigeta3_jmin(i,k) &
             +2.0_wp*(prs(i,1,k+1)-prs(i,1,k))*gphigeta3_jmin(i,k) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,2,k))

       ! conservative variables
       rho_n(i,1,k) =ro_wall
       rhoe_n(i,1,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(i,2,k)-prs(i,3,k)        &
             +2.0_wp*(prs(i+1,1,k)-prs(i,1,k))*gksigeta3_jmin(i,k) &
             +2.0_wp*(prs(i,1,k)-prs(i,1,k-1))*gphigeta3_jmin(i,k) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,2,k))

       ! conservative variables
       rho_n(i,1,k) =ro_wall
       rhoe_n(i,1,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(i,2,k)-prs(i,3,k)        &
             +2.0_wp*(prs(i,1,k)-prs(i-1,1,k))*gksigeta3_jmin(i,k) &
             +2.0_wp*(prs(i,1,k+1)-prs(i,1,k))*gphigeta3_jmin(i,k) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,2,k))

       ! conservative variables
       rho_n(i,1,k) =ro_wall
       rhoe_n(i,1,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(i,2,k)-prs(i,3,k)        &
             +2.0_wp*(prs(i,1,k)-prs(i-1,1,k))*gksigeta3_jmin(i,k) &
             +2.0_wp*(prs(i,1,k)-prs(i,1,k-1))*gphigeta3_jmin(i,k) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,2,k))

       ! conservative variables
       rho_n(i,1,k) =ro_wall
       rhoe_n(i,1,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
    
  end subroutine bc_wall_jmin_adiabatic_dpdn_c3

  !==============================================================================
  module subroutine bc_wall_jmin_isotherm_dpdn_c3
  !==============================================================================
    !> Apply wall Boundary Conditions at jmin: isothermal + dp/dn=0
    !> - 3D curvilinear version - (all variables are imposed)
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

    do k=ndz_jmin,nfz_jmin
       do i=ndx_jmin,nfx_jmin
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(i,2,k)-prs(i,3,k)        &
                +(prs(i+1,1,k)-prs(i-1,1,k))*gksigeta3_jmin(i,k) &
                +(prs(i,1,k+1)-prs(i,1,k-1))*gphigeta3_jmin(i,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,2,k))

          ! conservative variables
          rho_n(i,1,k) =ro_wall
          rhoe_n(i,1,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    enddo

    ! edge with other BC at imin
    ! --------------------------
    if (ndx_jmin==2) then
       i=1
       do k=ndz_jmin,nfz_jmin
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(i,2,k)-prs(i,3,k)        &
                +2.0_wp*(prs(i+1,1,k)-prs(i,1,k))*gksigeta3_jmin(i,k) &
                +(prs(i,1,k+1)-prs(i,1,k-1))*gphigeta3_jmin(i,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,2,k))

          ! conservative variables
          rho_n(i,1,k) =ro_wall
          rhoe_n(i,1,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! edge with other BC at imax
    ! --------------------------
    if (ndx_jmin==nx-1) then
       i=nx
       do k=ndz_jmin,nfz_jmin
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(i,2,k)-prs(i,3,k)        &
                +2.0_wp*(prs(i,1,k)-prs(i-1,1,k))*gksigeta3_jmin(i,k) &
                +(prs(i,1,k+1)-prs(i,1,k-1))*gphigeta3_jmin(i,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,2,k))

          ! conservative variables
          rho_n(i,1,k) =ro_wall
          rhoe_n(i,1,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! edge with other BC at kmin
    ! --------------------------
    if (ndz_jmin==2) then
       k=1
       do i=ndx_jmin,nfx_jmin
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(i,2,k)-prs(i,3,k)        &
                +(prs(i+1,1,k)-prs(i-1,1,k))*gksigeta3_jmin(i,k) &
                +2.0_wp*(prs(i,1,k+1)-prs(i,1,k))*gphigeta3_jmin(i,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,2,k))

          ! conservative variables
          rho_n(i,1,k) =ro_wall
          rhoe_n(i,1,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! edge with other BC at kmax
    ! --------------------------
    if (ndz_jmin==nz-1) then
       k=nz
       do i=ndx_jmin,nfx_jmin
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(i,2,k)-prs(i,3,k)        &
                +(prs(i+1,1,k)-prs(i-1,1,k))*gksigeta3_jmin(i,k) &
                +2.0_wp*(prs(i,1,k)-prs(i,1,k-1))*gphigeta3_jmin(i,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,2,k))

          ! conservative variables
          rho_n(i,1,k) =ro_wall
          rhoe_n(i,1,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! corner with other BC at imin-kmin
    ! ---------------------------------
    if ((ndx_jmin==2).and.(ndz_jmin==2)) then
       i=1
       k=1
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(i,2,k)-prs(i,3,k)        &
             +2.0_wp*(prs(i+1,1,k)-prs(i,1,k))*gksigeta3_jmin(i,k) &
             +2.0_wp*(prs(i,1,k+1)-prs(i,1,k))*gphigeta3_jmin(i,k) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,2,k))

       ! conservative variables
       rho_n(i,1,k) =ro_wall
       rhoe_n(i,1,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
    endif
    
    ! corner with other BC at imin-kmax
    ! ---------------------------------
    if ((ndx_jmin==2).and.(ndz_jmin==nz-1)) then
       i=1
       k=nz
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(i,2,k)-prs(i,3,k)        &
             +2.0_wp*(prs(i+1,1,k)-prs(i,1,k))*gksigeta3_jmin(i,k) &
             +2.0_wp*(prs(i,1,k)-prs(i,1,k-1))*gphigeta3_jmin(i,k) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,2,k))

       ! conservative variables
       rho_n(i,1,k) =ro_wall
       rhoe_n(i,1,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
    endif
    
    ! corner with other BC at imax-kmin
    ! ---------------------------------
    if ((ndx_jmin==nx-1).and.(ndz_jmin==2)) then
       i=nx
       k=1
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(i,2,k)-prs(i,3,k)        &
             +2.0_wp*(prs(i,1,k)-prs(i-1,1,k))*gksigeta3_jmin(i,k) &
             +2.0_wp*(prs(i,1,k+1)-prs(i,1,k))*gphigeta3_jmin(i,k) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,2,k))

       ! conservative variables
       rho_n(i,1,k) =ro_wall
       rhoe_n(i,1,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
    endif

    ! corner with other BC at imax-kmax
    ! ---------------------------------
    if ((ndx_jmin==nx-1).and.(ndz_jmin==nz-1)) then
       i=nx
       k=nz
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(i,2,k)-prs(i,3,k)        &
             +2.0_wp*(prs(i,1,k)-prs(i-1,1,k))*gksigeta3_jmin(i,k) &
             +2.0_wp*(prs(i,1,k)-prs(i,1,k-1))*gphigeta3_jmin(i,k) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,2,k))

       ! conservative variables
       rho_n(i,1,k) =ro_wall
       rhoe_n(i,1,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
    
  end subroutine bc_wall_jmin_isotherm_dpdn_c3

  !==============================================================================
  module subroutine bc_wall_jmax_adiabatic_dpdn_c3
  !==============================================================================
    !> Apply wall Boundary Conditions at jmax: adiabatic + dp/dn=0
    !> - 3D curvilinear version - (all variables are imposed)
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

    do k=ndz_jmax,nfz_jmax
       do i=ndx_jmax,nfx_jmax
          ! dT/dn=0
          T_wall= onethird*( 4.0_wp*Tmp(i,ny-1,k)-Tmp(i,ny-2,k)    &
                +(Tmp(i+1,ny,k)-Tmp(i-1,ny,k))*gksigeta3_jmax(i,k) &
                +(Tmp(i,ny,k+1)-Tmp(i,ny,k-1))*gphigeta3_jmax(i,k) )
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(i,ny-1,k)-prs(i,ny-2,k)    &
                +(prs(i+1,ny,k)-prs(i-1,ny,k))*gksigeta3_jmax(i,k) &
                +(prs(i,ny,k+1)-prs(i,ny,k-1))*gphigeta3_jmax(i,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,ny-1,k) )

          ! conservative variables
          rho_n(i,ny,k) =ro_wall
          rhoe_n(i,ny,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(i,ny-1,k)-prs(i,ny-2,k)    &
                +2.0_wp*(prs(i+1,ny,k)-prs(i,ny,k))*gksigeta3_jmax(i,k) &
                +(prs(i,ny,k+1)-prs(i,ny,k-1))*gphigeta3_jmax(i,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,ny-1,k) )

          ! conservative variables
          rho_n(i,ny,k) =ro_wall
          rhoe_n(i,ny,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(i,ny-1,k)-prs(i,ny-2,k)    &
                +2.0_wp*(prs(i,ny,k)-prs(i-1,ny,k))*gksigeta3_jmax(i,k) &
                +(prs(i,ny,k+1)-prs(i,ny,k-1))*gphigeta3_jmax(i,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,ny-1,k) )

          ! conservative variables
          rho_n(i,ny,k) =ro_wall
          rhoe_n(i,ny,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(i,ny-1,k)-prs(i,ny-2,k)    &
                +(prs(i+1,ny,k)-prs(i-1,ny,k))*gksigeta3_jmax(i,k) &
                +2.0_wp*(prs(i,ny,k+1)-prs(i,ny,k))*gphigeta3_jmax(i,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,ny-1,k) )

          ! conservative variables
          rho_n(i,ny,k) =ro_wall
          rhoe_n(i,ny,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(i,ny-1,k)-prs(i,ny-2,k)    &
                +(prs(i+1,ny,k)-prs(i-1,ny,k))*gksigeta3_jmax(i,k) &
                +2.0_wp*(prs(i,ny,k)-prs(i,ny,k-1))*gphigeta3_jmax(i,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,ny-1,k) )

          ! conservative variables
          rho_n(i,ny,k) =ro_wall
          rhoe_n(i,ny,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(i,ny-1,k)-prs(i,ny-2,k)    &
             +2.0_wp*(prs(i+1,ny,k)-prs(i,ny,k))*gksigeta3_jmax(i,k) &
             +2.0_wp*(prs(i,ny,k+1)-prs(i,ny,k))*gphigeta3_jmax(i,k) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,ny-1,k) )

       ! conservative variables
       rho_n(i,ny,k) =ro_wall
       rhoe_n(i,ny,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(i,ny-1,k)-prs(i,ny-2,k)    &
             +2.0_wp*(prs(i+1,ny,k)-prs(i,ny,k))*gksigeta3_jmax(i,k) &
             +2.0_wp*(prs(i,ny,k)-prs(i,ny,k-1))*gphigeta3_jmax(i,k) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,ny-1,k) )

       ! conservative variables
       rho_n(i,ny,k) =ro_wall
       rhoe_n(i,ny,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(i,ny-1,k)-prs(i,ny-2,k)    &
             +2.0_wp*(prs(i,ny,k)-prs(i-1,ny,k))*gksigeta3_jmax(i,k) &
             +2.0_wp*(prs(i,ny,k+1)-prs(i,ny,k))*gphigeta3_jmax(i,k) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,ny-1,k) )

       ! conservative variables
       rho_n(i,ny,k) =ro_wall
       rhoe_n(i,ny,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(i,ny-1,k)-prs(i,ny-2,k)    &
             +2.0_wp*(prs(i,ny,k)-prs(i-1,ny,k))*gksigeta3_jmax(i,k) &
             +2.0_wp*(prs(i,ny,k)-prs(i,ny,k-1))*gphigeta3_jmax(i,k) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,ny-1,k) )

       ! conservative variables
       rho_n(i,ny,k) =ro_wall
       rhoe_n(i,ny,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
    
  end subroutine bc_wall_jmax_adiabatic_dpdn_c3

  !==============================================================================
  module subroutine bc_wall_jmax_isotherm_dpdn_c3
  !==============================================================================
    !> Apply wall Boundary Conditions at jmax: isothermal + dp/dn=0
    !> - 3D curvilinear version - (all variables are imposed)
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

    do k=ndz_jmax,nfz_jmax
       do i=ndx_jmax,nfx_jmax
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(i,ny-1,k)-prs(i,ny-2,k)    &
                +(prs(i+1,ny,k)-prs(i-1,ny,k))*gksigeta3_jmax(i,k) &
                +(prs(i,ny,k+1)-prs(i,ny,k-1))*gphigeta3_jmax(i,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,ny-1,k) )

          ! conservative variables
          rho_n(i,ny,k) =ro_wall
          rhoe_n(i,ny,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    enddo

    ! edge with other BC at imin
    ! --------------------------
    if (ndx_jmax==2) then
       i=1
       do k=ndz_jmax,nfz_jmax
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(i,ny-1,k)-prs(i,ny-2,k)    &
                +2.0_wp*(prs(i+1,ny,k)-prs(i,ny,k))*gksigeta3_jmax(i,k) &
                +(prs(i,ny,k+1)-prs(i,ny,k-1))*gphigeta3_jmax(i,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,ny-1,k) )

          ! conservative variables
          rho_n(i,ny,k) =ro_wall
          rhoe_n(i,ny,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! edge with other BC at imax
    ! --------------------------
    if (ndx_jmax==nx-1) then
       i=nx
       do k=ndz_jmax,nfz_jmax
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(i,ny-1,k)-prs(i,ny-2,k)    &
                +2.0_wp*(prs(i,ny,k)-prs(i-1,ny,k))*gksigeta3_jmax(i,k) &
                +(prs(i,ny,k+1)-prs(i,ny,k-1))*gphigeta3_jmax(i,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,ny-1,k) )

          ! conservative variables
          rho_n(i,ny,k) =ro_wall
          rhoe_n(i,ny,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! edge with other BC at kmin
    ! --------------------------
    if (ndz_jmax==2) then
       k=1
       do i=ndx_jmax,nfx_jmax
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(i,ny-1,k)-prs(i,ny-2,k)    &
                +(prs(i+1,ny,k)-prs(i-1,ny,k))*gksigeta3_jmax(i,k) &
                +2.0_wp*(prs(i,ny,k+1)-prs(i,ny,k))*gphigeta3_jmax(i,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,ny-1,k) )

          ! conservative variables
          rho_n(i,ny,k) =ro_wall
          rhoe_n(i,ny,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! edge with other BC at kmax
    ! --------------------------
    if (ndz_jmax==nz-1) then
       k=nz
       do i=ndx_jmax,nfx_jmax
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(i,ny-1,k)-prs(i,ny-2,k)    &
                +(prs(i+1,ny,k)-prs(i-1,ny,k))*gksigeta3_jmax(i,k) &
                +2.0_wp*(prs(i,ny,k)-prs(i,ny,k-1))*gphigeta3_jmax(i,k) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,ny-1,k) )

          ! conservative variables
          rho_n(i,ny,k) =ro_wall
          rhoe_n(i,ny,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! corner with other BC at imin-kmin
    ! ---------------------------------
    if ((ndx_jmax==2).and.(ndz_jmax==2)) then
       i=1
       k=1
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(i,ny-1,k)-prs(i,ny-2,k)    &
             +2.0_wp*(prs(i+1,ny,k)-prs(i,ny,k))*gksigeta3_jmax(i,k) &
             +2.0_wp*(prs(i,ny,k+1)-prs(i,ny,k))*gphigeta3_jmax(i,k) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,ny-1,k) )

       ! conservative variables
       rho_n(i,ny,k) =ro_wall
       rhoe_n(i,ny,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
    endif
    
    ! corner with other BC at imin-kmax
    ! ---------------------------------
    if ((ndx_jmax==2).and.(ndz_jmax==nz-1)) then
       i=1
       k=nz
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(i,ny-1,k)-prs(i,ny-2,k)    &
             +2.0_wp*(prs(i+1,ny,k)-prs(i,ny,k))*gksigeta3_jmax(i,k) &
             +2.0_wp*(prs(i,ny,k)-prs(i,ny,k-1))*gphigeta3_jmax(i,k) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,ny-1,k) )

       ! conservative variables
       rho_n(i,ny,k) =ro_wall
       rhoe_n(i,ny,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
    endif
    
    ! corner with other BC at imax-kmin
    ! ---------------------------------
    if ((ndx_jmax==nx-1).and.(ndz_jmax==2)) then
       i=nx
       k=1
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(i,ny-1,k)-prs(i,ny-2,k)    &
             +2.0_wp*(prs(i,ny,k)-prs(i-1,ny,k))*gksigeta3_jmax(i,k) &
             +2.0_wp*(prs(i,ny,k+1)-prs(i,ny,k))*gphigeta3_jmax(i,k) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,ny-1,k) )

       ! conservative variables
       rho_n(i,ny,k) =ro_wall
       rhoe_n(i,ny,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
    endif

    ! corner with other BC at imax-kmax
    ! ---------------------------------
    if ((ndx_jmax==nx-1).and.(ndz_jmax==nz-1)) then
       i=nx
       k=nz
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(i,ny-1,k)-prs(i,ny-2,k)    &
             +2.0_wp*(prs(i,ny,k)-prs(i-1,ny,k))*gksigeta3_jmax(i,k) &
             +2.0_wp*(prs(i,ny,k)-prs(i,ny,k-1))*gphigeta3_jmax(i,k) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,ny-1,k) )

       ! conservative variables
       rho_n(i,ny,k) =ro_wall
       rhoe_n(i,ny,k)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
    
  end subroutine bc_wall_jmax_isotherm_dpdn_c3

  !==============================================================================
  module subroutine bc_wall_kmin_adiabatic_dpdn_c3
  !==============================================================================
    !> Apply wall Boundary Conditions at kmin: adiabatic + dp/dn=0
    !> - 3D curvilinear version - (all variables are imposed)
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
       enddo
    enddo

    do j=ndy_kmin,nfy_kmin
       do i=ndx_kmin,nfx_kmin
          ! dT/dn=0
          T_wall= onethird*( 4.0_wp*Tmp(i,j,2)-Tmp(i,j,3)        &
                +(Tmp(i+1,j,1)-Tmp(i-1,j,1))*gksigphi3_kmin(i,j) &
                +(Tmp(i,j+1,1)-Tmp(i,j-1,1))*getagphi3_kmin(i,j) )
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(i,j,2)-prs(i,j,3)         &
                +(prs(i+1,j,1)-prs(i-1,j,1))*gksigphi3_kmin(i,j) &
                +(prs(i,j+1,1)-prs(i,j-1,1))*getagphi3_kmin(i,j) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,2) )

          ! conservative variables
          rho_n(i,j,1) =ro_wall
          rhoe_n(i,j,1)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(i,j,2)-prs(i,j,3)         &
                +2.0_wp*(prs(i+1,j,1)-prs(i,j,1))*gksigphi3_kmin(i,j) &
                +(prs(i,j+1,1)-prs(i,j-1,1))*getagphi3_kmin(i,j) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,2) )

          ! conservative variables
          rho_n(i,j,1) =ro_wall
          rhoe_n(i,j,1)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(i,j,2)-prs(i,j,3)         &
                +2.0_wp*(prs(i,j,1)-prs(i-1,j,1))*gksigphi3_kmin(i,j) &
                +(prs(i,j+1,1)-prs(i,j-1,1))*getagphi3_kmin(i,j) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,2) )

          ! conservative variables
          rho_n(i,j,1) =ro_wall
          rhoe_n(i,j,1)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(i,j,2)-prs(i,j,3)         &
                +(prs(i+1,j,1)-prs(i-1,j,1))*gksigphi3_kmin(i,j) &
                +2.0_wp*(prs(i,j+1,1)-prs(i,j,1))*getagphi3_kmin(i,j) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,2) )

          ! conservative variables
          rho_n(i,j,1) =ro_wall
          rhoe_n(i,j,1)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(i,j,2)-prs(i,j,3)         &
                +(prs(i+1,j,1)-prs(i-1,j,1))*gksigphi3_kmin(i,j) &
                +2.0_wp*(prs(i,j,1)-prs(i,j-1,1))*getagphi3_kmin(i,j) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,2) )

          ! conservative variables
          rho_n(i,j,1) =ro_wall
          rhoe_n(i,j,1)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
       ! dp/dn=0
       p_wall= onethird*(4.0_wp*prs(i,j,2)-prs(i,j,3)         &
             +2.0_wp*(prs(i+1,j,1)-prs(i,j,1))*gksigphi3_kmin(i,j) &
             +2.0_wp*(prs(i,j+1,1)-prs(i,j,1))*getagphi3_kmin(i,j) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,2) )

       ! conservative variables
       rho_n(i,j,1) =ro_wall
       rhoe_n(i,j,1)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
       ! dp/dn=0
       p_wall= onethird*(4.0_wp*prs(i,j,2)-prs(i,j,3)         &
             +2.0_wp*(prs(i+1,j,1)-prs(i,j,1))*gksigphi3_kmin(i,j) &
             +2.0_wp*(prs(i,j,1)-prs(i,j-1,1))*getagphi3_kmin(i,j) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,2) )

       ! conservative variables
       rho_n(i,j,1) =ro_wall
       rhoe_n(i,j,1)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
       ! dp/dn=0
       p_wall= onethird*(4.0_wp*prs(i,j,2)-prs(i,j,3)         &
             +2.0_wp*(prs(i,j,1)-prs(i-1,j,1))*gksigphi3_kmin(i,j) &
             +2.0_wp*(prs(i,j+1,1)-prs(i,j,1))*getagphi3_kmin(i,j) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,2) )

       ! conservative variables
       rho_n(i,j,1) =ro_wall
       rhoe_n(i,j,1)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
       ! dp/dn=0
       p_wall= onethird*(4.0_wp*prs(i,j,2)-prs(i,j,3)         &
             +2.0_wp*(prs(i,j,1)-prs(i-1,j,1))*gksigphi3_kmin(i,j) &
             +2.0_wp*(prs(i,j,1)-prs(i,j-1,1))*getagphi3_kmin(i,j) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,2) )

       ! conservative variables
       rho_n(i,j,1) =ro_wall
       rhoe_n(i,j,1)=ro_wall*ecalc_tro(T_wall,ro_wall)
    endif

    ! Set increments to zero
    ! ======================
    do i=ndx_e,nfx_e
       Krho(i,ndy_e:nfy_e,1)=0.0_wp
       Krhou(i,ndy_e:nfy_e,1)=0.0_wp
       Krhov(i,ndy_e:nfy_e,1)=0.0_wp
       Krhow(i,ndy_e:nfy_e,1)=0.0_wp
       Krhoe(i,ndy_e:nfy_e,1)=0.0_wp
    enddo
    
  end subroutine bc_wall_kmin_adiabatic_dpdn_c3

  !==============================================================================
  module subroutine bc_wall_kmin_isotherm_dpdn_c3
  !==============================================================================
    !> Apply wall Boundary Conditions at kmin: isothermal + dp/dn=0
    !> - 3D curvilinear version - (all variables are imposed)
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
       enddo
    enddo

    do j=ndy_kmin,nfy_kmin
       do i=ndx_kmin,nfx_kmin
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(i,j,2)-prs(i,j,3)        &
                +(prs(i+1,j,1)-prs(i-1,j,1))*gksigphi3_kmin(i,j) &
                +(prs(i,j+1,1)-prs(i,j-1,1))*getagphi3_kmin(i,j) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,2))

          ! conservative variables
          rho_n(i,j,1) =ro_wall
          rhoe_n(i,j,1)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    enddo

    ! edge with other BC at imin
    ! --------------------------
    if (ndx_kmin==2) then
       i=1
       do j=ndy_kmin,nfy_kmin
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(i,j,2)-prs(i,j,3)         &
                +2.0_wp*(prs(i+1,j,1)-prs(i,j,1))*gksigphi3_kmin(i,j) &
                +(prs(i,j+1,1)-prs(i,j-1,1))*getagphi3_kmin(i,j) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,2) )

          ! conservative variables
          rho_n(i,j,1) =ro_wall
          rhoe_n(i,j,1)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! edge with other BC at imax
    ! --------------------------
    if (ndx_kmin==nx-1) then
       i=nx
       do j=ndy_kmin,nfy_kmin
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(i,j,2)-prs(i,j,3)         &
                +2.0_wp*(prs(i,j,1)-prs(i-1,j,1))*gksigphi3_kmin(i,j) &
                +(prs(i,j+1,1)-prs(i,j-1,1))*getagphi3_kmin(i,j) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,2) )

          ! conservative variables
          rho_n(i,j,1) =ro_wall
          rhoe_n(i,j,1)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! edge with other BC at jmin
    ! --------------------------
    if (ndy_kmin==2) then
       j=1
       do i=ndx_kmin,nfx_kmin
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(i,j,2)-prs(i,j,3)         &
                +(prs(i+1,j,1)-prs(i-1,j,1))*gksigphi3_kmin(i,j) &
                +2.0_wp*(prs(i,j+1,1)-prs(i,j,1))*getagphi3_kmin(i,j) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,2) )

          ! conservative variables
          rho_n(i,j,1) =ro_wall
          rhoe_n(i,j,1)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! edge with other BC at jmax
    ! --------------------------
    if (ndy_kmin==ny-1) then
       j=ny
       do i=ndx_kmin,nfx_kmin
          ! dp/dn=0
          p_wall= onethird*(4.0_wp*prs(i,j,2)-prs(i,j,3)         &
                +(prs(i+1,j,1)-prs(i-1,j,1))*gksigphi3_kmin(i,j) &
                +2.0_wp*(prs(i,j,1)-prs(i,j-1,1))*getagphi3_kmin(i,j) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,2) )

          ! conservative variables
          rho_n(i,j,1) =ro_wall
          rhoe_n(i,j,1)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! corner with other BC at imin-jmin
    ! ---------------------------------
    if ((ndx_kmin==2).and.(ndy_kmin==2)) then
       i=1
       j=1
       ! dp/dn=0
       p_wall= onethird*(4.0_wp*prs(i,j,2)-prs(i,j,3)         &
             +2.0_wp*(prs(i+1,j,1)-prs(i,j,1))*gksigphi3_kmin(i,j) &
             +2.0_wp*(prs(i,j+1,1)-prs(i,j,1))*getagphi3_kmin(i,j) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,2) )

       ! conservative variables
       rho_n(i,j,1) =ro_wall
       rhoe_n(i,j,1)=ro_wall*ecalc_tro(T_wall,ro_wall)
    endif
    
    ! corner with other BC at imin-jmax
    ! ---------------------------------
    if ((ndx_kmin==2).and.(ndy_kmin==ny-1)) then
       i=1
       j=ny
       ! dp/dn=0
       p_wall= onethird*(4.0_wp*prs(i,j,2)-prs(i,j,3)         &
             +2.0_wp*(prs(i+1,j,1)-prs(i,j,1))*gksigphi3_kmin(i,j) &
             +2.0_wp*(prs(i,j,1)-prs(i,j-1,1))*getagphi3_kmin(i,j) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,2) )

       ! conservative variables
       rho_n(i,j,1) =ro_wall
       rhoe_n(i,j,1)=ro_wall*ecalc_tro(T_wall,ro_wall)
    endif

    ! corner with other BC at imax-jmin
    ! ---------------------------------
    if ((ndx_kmin==nx-1).and.(ndy_kmin==2)) then
       i=nx
       j=1
       ! dp/dn=0
       p_wall= onethird*(4.0_wp*prs(i,j,2)-prs(i,j,3)         &
             +2.0_wp*(prs(i,j,1)-prs(i-1,j,1))*gksigphi3_kmin(i,j) &
             +2.0_wp*(prs(i,j+1,1)-prs(i,j,1))*getagphi3_kmin(i,j) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,2) )

       ! conservative variables
       rho_n(i,j,1) =ro_wall
       rhoe_n(i,j,1)=ro_wall*ecalc_tro(T_wall,ro_wall)
    endif
        
    ! corner with other BC at imax-jmax
    ! ---------------------------------
    if ((ndx_kmin==nx-1).and.(ndy_kmin==ny-1)) then
       i=nx
       j=ny
       ! dp/dn=0
       p_wall= onethird*(4.0_wp*prs(i,j,2)-prs(i,j,3)         &
             +2.0_wp*(prs(i,j,1)-prs(i-1,j,1))*gksigphi3_kmin(i,j) &
             +2.0_wp*(prs(i,j,1)-prs(i,j-1,1))*getagphi3_kmin(i,j) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,2) )

       ! conservative variables
       rho_n(i,j,1) =ro_wall
       rhoe_n(i,j,1)=ro_wall*ecalc_tro(T_wall,ro_wall)
    endif

    ! Set increments to zero
    ! ======================
    do i=ndx_e,nfx_e
       Krho(i,ndy_e:nfy_e,1)=0.0_wp
       Krhou(i,ndy_e:nfy_e,1)=0.0_wp
       Krhov(i,ndy_e:nfy_e,1)=0.0_wp
       Krhow(i,ndy_e:nfy_e,1)=0.0_wp
       Krhoe(i,ndy_e:nfy_e,1)=0.0_wp
    enddo
    
  end subroutine bc_wall_kmin_isotherm_dpdn_c3

  !==============================================================================
  module subroutine bc_wall_kmax_adiabatic_dpdn_c3
  !==============================================================================
    !> Apply wall Boundary Conditions at kmax: adiabatic + dp/dn=0
    !> - 3D curvilinear version - (all variables are imposed)
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
       enddo
    enddo

    do j=ndy_kmin,nfy_kmin
       do i=ndx_kmax,nfx_kmax
          ! dT/dn=0
          T_wall= onethird*( 4.0_wp*Tmp(i,j,nz-1)-Tmp(i,j,nz-2)    &
                +(Tmp(i+1,j,nz)-Tmp(i-1,j,nz))*gksigphi3_kmax(i,j) &
                +(Tmp(i,j+1,nz)-Tmp(i,j-1,nz))*getagphi3_kmax(i,j) )
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(i,j,nz-1)-prs(i,j,nz-2)    &
                +(prs(i+1,j,nz)-prs(i-1,j,nz))*gksigphi3_kmax(i,j) &
                +(prs(i,j+1,nz)-prs(i,j-1,nz))*getagphi3_kmax(i,j) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,nz-1) )

          ! conservative variables
          rho_n(i,j,nz) =ro_wall
          rhoe_n(i,j,nz)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(i,j,nz-1)-prs(i,j,nz-2)    &
                +2.0_wp*(prs(i+1,j,nz)-prs(i,j,nz))*gksigphi3_kmax(i,j) &
                +(prs(i,j+1,nz)-prs(i,j-1,nz))*getagphi3_kmax(i,j) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,nz-1) )

          ! conservative variables
          rho_n(i,j,nz) =ro_wall
          rhoe_n(i,j,nz)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(i,j,nz-1)-prs(i,j,nz-2)    &
                +2.0_wp*(prs(i,j,nz)-prs(i-1,j,nz))*gksigphi3_kmax(i,j) &
                +(prs(i,j+1,nz)-prs(i,j-1,nz))*getagphi3_kmax(i,j) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,nz-1) )

          ! conservative variables
          rho_n(i,j,nz) =ro_wall
          rhoe_n(i,j,nz)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(i,j,nz-1)-prs(i,j,nz-2)    &
                +(prs(i+1,j,nz)-prs(i-1,j,nz))*gksigphi3_kmax(i,j) &
                +2.0_wp*(prs(i,j+1,nz)-prs(i,j,nz))*getagphi3_kmax(i,j) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,nz-1) )

          ! conservative variables
          rho_n(i,j,nz) =ro_wall
          rhoe_n(i,j,nz)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(i,j,nz-1)-prs(i,j,nz-2)    &
                +(prs(i+1,j,nz)-prs(i-1,j,nz))*gksigphi3_kmax(i,j) &
                +2.0_wp*(prs(i,j,nz)-prs(i,j-1,nz))*getagphi3_kmax(i,j) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,nz-1) )

          ! conservative variables
          rho_n(i,j,nz) =ro_wall
          rhoe_n(i,j,nz)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(i,j,nz-1)-prs(i,j,nz-2)    &
             +2.0_wp*(prs(i+1,j,nz)-prs(i,j,nz))*gksigphi3_kmax(i,j) &
             +2.0_wp*(prs(i,j+1,nz)-prs(i,j,nz))*getagphi3_kmax(i,j) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,nz-1) )

       ! conservative variables
       rho_n(i,j,nz) =ro_wall
       rhoe_n(i,j,nz)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(i,j,nz-1)-prs(i,j,nz-2)    &
             +2.0_wp*(prs(i+1,j,nz)-prs(i,j,nz))*gksigphi3_kmax(i,j) &
             +2.0_wp*(prs(i,j,nz)-prs(i,j-1,nz))*getagphi3_kmax(i,j) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,nz-1) )

       ! conservative variables
       rho_n(i,j,nz) =ro_wall
       rhoe_n(i,j,nz)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(i,j,nz-1)-prs(i,j,nz-2)    &
             +2.0_wp*(prs(i,j,nz)-prs(i-1,j,nz))*gksigphi3_kmax(i,j) &
             +2.0_wp*(prs(i,j+1,nz)-prs(i,j,nz))*getagphi3_kmax(i,j) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,nz-1) )

       ! conservative variables
       rho_n(i,j,nz) =ro_wall
       rhoe_n(i,j,nz)=ro_wall*ecalc_tro(T_wall,ro_wall)
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
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(i,j,nz-1)-prs(i,j,nz-2)    &
             +2.0_wp*(prs(i,j,nz)-prs(i-1,j,nz))*gksigphi3_kmax(i,j) &
             +2.0_wp*(prs(i,j,nz)-prs(i,j-1,nz))*getagphi3_kmax(i,j) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,nz-1) )

       ! conservative variables
       rho_n(i,j,nz) =ro_wall
       rhoe_n(i,j,nz)=ro_wall*ecalc_tro(T_wall,ro_wall)
    endif

    ! Set increments to zero
    ! ======================
    do i=ndx_e,nfx_e
       Krho(i,ndy_e:nfy_e,nz)=0.0_wp
       Krhou(i,ndy_e:nfy_e,nz)=0.0_wp
       Krhov(i,ndy_e:nfy_e,nz)=0.0_wp
       Krhow(i,ndy_e:nfy_e,nz)=0.0_wp
       Krhoe(i,ndy_e:nfy_e,nz)=0.0_wp
    enddo
    
  end subroutine bc_wall_kmax_adiabatic_dpdn_c3

  !==============================================================================
  module subroutine bc_wall_kmax_isotherm_dpdn_c3
  !==============================================================================
    !> Apply wall Boundary Conditions at kmax: isothermal + dp/dn=0
    !> - 3D curvilinear version - (all variables are imposed)
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
       enddo
    enddo

    do j=ndy_kmin,nfy_kmin
       do i=ndx_kmax,nfx_kmax
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(i,j,nz-1)-prs(i,j,nz-2)    &
                +(prs(i+1,j,nz)-prs(i-1,j,nz))*gksigphi3_kmax(i,j) &
                +(prs(i,j+1,nz)-prs(i,j-1,nz))*getagphi3_kmax(i,j) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,nz-1) )

          ! conservative variables
          rho_n(i,j,nz) =ro_wall
          rhoe_n(i,j,nz)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    enddo

    ! edge with other BC at imin
    ! --------------------------
    if (ndx_kmax==2) then
       i=1
       do j=ndy_kmax,nfy_kmax
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(i,j,nz-1)-prs(i,j,nz-2)    &
                +2.0_wp*(prs(i+1,j,nz)-prs(i,j,nz))*gksigphi3_kmax(i,j) &
                +(prs(i,j+1,nz)-prs(i,j-1,nz))*getagphi3_kmax(i,j) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,nz-1) )

          ! conservative variables
          rho_n(i,j,nz) =ro_wall
          rhoe_n(i,j,nz)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! edge with other BC at imax
    ! --------------------------
    if (ndx_kmax==nx-1) then
       i=nx
       do j=ndy_kmax,nfy_kmax
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(i,j,nz-1)-prs(i,j,nz-2)    &
                +2.0_wp*(prs(i,j,nz)-prs(i-1,j,nz))*gksigphi3_kmax(i,j) &
                +(prs(i,j+1,nz)-prs(i,j-1,nz))*getagphi3_kmax(i,j) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,nz-1) )

          ! conservative variables
          rho_n(i,j,nz) =ro_wall
          rhoe_n(i,j,nz)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! edge with other BC at jmin
    ! --------------------------
    if (ndy_kmax==2) then
       j=1
       do i=ndx_kmax,nfx_kmax
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(i,j,nz-1)-prs(i,j,nz-2)    &
                +(prs(i+1,j,nz)-prs(i-1,j,nz))*gksigphi3_kmax(i,j) &
                +2.0_wp*(prs(i,j+1,nz)-prs(i,j,nz))*getagphi3_kmax(i,j) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,nz-1) )

          ! conservative variables
          rho_n(i,j,nz) =ro_wall
          rhoe_n(i,j,nz)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! edge with other BC at jmax
    ! --------------------------
    if (ndy_kmax==ny-1) then
       j=ny
       do i=ndx_kmax,nfx_kmax
          ! dp/dn=0
          p_wall= onethird*( 4.0_wp*prs(i,j,nz-1)-prs(i,j,nz-2)    &
                +(prs(i+1,j,nz)-prs(i-1,j,nz))*gksigphi3_kmax(i,j) &
                +2.0_wp*(prs(i,j,nz)-prs(i,j-1,nz))*getagphi3_kmax(i,j) )
          ! calc rho
          ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,nz-1) )

          ! conservative variables
          rho_n(i,j,nz) =ro_wall
          rhoe_n(i,j,nz)=ro_wall*ecalc_tro(T_wall,ro_wall)
       enddo
    endif

    ! corner with other BC at imin-jmin
    ! ---------------------------------
    if ((ndx_kmax==2).and.(ndy_kmax==2)) then
       i=1
       j=1
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(i,j,nz-1)-prs(i,j,nz-2)    &
             +2.0_wp*(prs(i+1,j,nz)-prs(i,j,nz))*gksigphi3_kmax(i,j) &
             +2.0_wp*(prs(i,j+1,nz)-prs(i,j,nz))*getagphi3_kmax(i,j) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,nz-1) )

       ! conservative variables
       rho_n(i,j,nz) =ro_wall
       rhoe_n(i,j,nz)=ro_wall*ecalc_tro(T_wall,ro_wall)
    endif
    
    ! corner with other BC at imin-jmax
    ! ---------------------------------
    if ((ndx_kmax==2).and.(ndy_kmax==ny-1)) then
       i=1
       j=ny
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(i,j,nz-1)-prs(i,j,nz-2)    &
             +2.0_wp*(prs(i+1,j,nz)-prs(i,j,nz))*gksigphi3_kmax(i,j) &
             +2.0_wp*(prs(i,j,nz)-prs(i,j-1,nz))*getagphi3_kmax(i,j) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,nz-1) )

       ! conservative variables
       rho_n(i,j,nz) =ro_wall
       rhoe_n(i,j,nz)=ro_wall*ecalc_tro(T_wall,ro_wall)
    endif

    ! corner with other BC at imax-jmin
    ! ---------------------------------
    if ((ndx_kmax==nx-1).and.(ndy_kmax==2)) then
       i=nx
       j=1
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(i,j,nz-1)-prs(i,j,nz-2)    &
             +2.0_wp*(prs(i,j,nz)-prs(i-1,j,nz))*gksigphi3_kmax(i,j) &
             +2.0_wp*(prs(i,j+1,nz)-prs(i,j,nz))*getagphi3_kmax(i,j) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,nz-1) )

       ! conservative variables
       rho_n(i,j,nz) =ro_wall
       rhoe_n(i,j,nz)=ro_wall*ecalc_tro(T_wall,ro_wall)
    endif
        
    ! corner with other BC at imax-jmax
    ! ---------------------------------
    if ((ndx_kmax==nx-1).and.(ndy_kmax==ny-1)) then
       i=nx
       j=ny
       ! dp/dn=0
       p_wall= onethird*( 4.0_wp*prs(i,j,nz-1)-prs(i,j,nz-2)    &
             +2.0_wp*(prs(i,j,nz)-prs(i-1,j,nz))*gksigphi3_kmax(i,j) &
             +2.0_wp*(prs(i,j,nz)-prs(i,j-1,nz))*getagphi3_kmax(i,j) )
       ! calc rho
       ro_wall=rocalc_pt(p_wall,T_wall,rho_n(i,j,nz-1) )

       ! conservative variables
       rho_n(i,j,nz) =ro_wall
       rhoe_n(i,j,nz)=ro_wall*ecalc_tro(T_wall,ro_wall)
    endif

    ! Set increments to zero
    ! ======================
    do i=ndx_e,nfx_e
       Krho(i,ndy_e:nfy_e,nz)=0.0_wp
       Krhou(i,ndy_e:nfy_e,nz)=0.0_wp
       Krhov(i,ndy_e:nfy_e,nz)=0.0_wp
       Krhow(i,ndy_e:nfy_e,nz)=0.0_wp
       Krhoe(i,ndy_e:nfy_e,nz)=0.0_wp
    enddo
    
  end subroutine bc_wall_kmax_isotherm_dpdn_c3

end submodule smod_bc_wall_c3
