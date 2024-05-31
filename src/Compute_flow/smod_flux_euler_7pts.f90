!=================================================================================
submodule (mod_flux_euler) smod_flux_euler_7pts
!=================================================================================
  !> author: XG
  !> date: January 2024
  !> Cartesian version - compute Eulerian fluxes (inviscid part) with 7-pt stencil
!=================================================================================

contains

  !==============================================================================
  module subroutine flux_euler_7pts
  !==============================================================================
    !> Derivatives of Eulerian fluxes (inviscid part) - 7-point stencil -
    !> - Cartesian version -
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    ! ---------------------------------------------------------------------------

    ! Computation of inviscid fluxes along x
    ! ======================================
    do k=ndz_e,nfz_e
       do j=ndy_e,nfy_e
          do i=ndxt,nfxt
             Frho(i,j,k) = rhou_n(i,j,k)
             Frhou(i,j,k)= rhou_n(i,j,k)*uu(i,j,k) +prs(i,j,k)
             Frhov(i,j,k)= rhou_n(i,j,k)*vv(i,j,k)
             Frhow(i,j,k)= rhou_n(i,j,k)*ww(i,j,k)
             Frhoe(i,j,k)=(rhoe_n(i,j,k)+prs(i,j,k))*uu(i,j,k)
          enddo
       enddo
    enddo

    ! Flux derivatives along x
    ! ========================
    
    ! Advancement of wall points
    ! --------------------------
    if (is_bc_wall2(1,1)) then
       !i=1
       do k=ndz_e,nfz_e
          do j=ndy_e,nfy_e
             Krho(1,j,k)= (a02(1)*Frho(1,j,k) +a02(2)*Frho(2,j,k) +a02(3)*Frho(3,j,k) )*idx1_imin
             Krhou(1,j,k)=(a02(1)*Frhou(1,j,k)+a02(2)*Frhou(2,j,k)+a02(3)*Frhou(3,j,k))*idx1_imin
             Krhov(1,j,k)=(a02(1)*Frhov(1,j,k)+a02(2)*Frhov(2,j,k)+a02(3)*Frhov(3,j,k))*idx1_imin
             Krhow(1,j,k)=(a02(1)*Frhow(1,j,k)+a02(2)*Frhow(2,j,k)+a02(3)*Frhow(3,j,k))*idx1_imin
             Krhoe(1,j,k)=(a02(1)*Frhoe(1,j,k)+a02(2)*Frhoe(2,j,k)+a02(3)*Frhoe(3,j,k))*idx1_imin
          enddo
       enddo
    endif

    ! BC at imin
    ! ----------
    if (is_bc_1pt(1,1)) then  
       i=2
       do k=ndz_e,nfz_e
          do j=ndy_e,nfy_e
             Krho(i,j,k)  = a3(1)*( Frho(i+1,j,k)- Frho(i-1,j,k))*idx2_imin
             Krhou(i,j,k) = a3(1)*(Frhou(i+1,j,k)-Frhou(i-1,j,k))*idx2_imin
             Krhov(i,j,k) = a3(1)*(Frhov(i+1,j,k)-Frhov(i-1,j,k))*idx2_imin
             Krhow(i,j,k) = a3(1)*(Frhow(i+1,j,k)-Frhow(i-1,j,k))*idx2_imin
             Krhoe(i,j,k) = a3(1)*(Frhoe(i+1,j,k)-Frhoe(i-1,j,k))*idx2_imin
          enddo
       enddo

       i=3
       do k=ndz_e,nfz_e
          do j=ndy_e,nfy_e
             Krho(i,j,k)  = ( a5(1)*( Frho(i+1,j,k)- Frho(i-1,j,k)) &
                            + a5(2)*( Frho(i+2,j,k)- Frho(i-2,j,k)) )*idx4_imin
             Krhou(i,j,k) = ( a5(1)*(Frhou(i+1,j,k)-Frhou(i-1,j,k)) &
                            + a5(2)*(Frhou(i+2,j,k)-Frhou(i-2,j,k)) )*idx4_imin
             Krhov(i,j,k) = ( a5(1)*(Frhov(i+1,j,k)-Frhov(i-1,j,k)) &
                            + a5(2)*(Frhov(i+2,j,k)-Frhov(i-2,j,k)) )*idx4_imin
             Krhow(i,j,k) = ( a5(1)*(Frhow(i+1,j,k)-Frhow(i-1,j,k)) &
                            + a5(2)*(Frhow(i+2,j,k)-Frhow(i-2,j,k)) )*idx4_imin
             Krhoe(i,j,k) = ( a5(1)*(Frhoe(i+1,j,k)-Frhoe(i-1,j,k)) &
                            + a5(2)*(Frhoe(i+2,j,k)-Frhoe(i-2,j,k)) )*idx4_imin
          enddo
       enddo
    endif

    ! Interior points
    ! ---------------
    do k=ndz_e,nfz_e
       do j=ndy_e,nfy_e
          do i=ndx,nfx

             Krho(i,j,k)  = ( a7(1) * ( Frho(i+1,j,k)-Frho(i-1,j,k) )   &
                            + a7(2) * ( Frho(i+2,j,k)-Frho(i-2,j,k) )   &
                            + a7(3) * ( Frho(i+3,j,k)-Frho(i-3,j,k) ) ) *idx(i)

             Krhou(i,j,k) = ( a7(1) * ( Frhou(i+1,j,k)-Frhou(i-1,j,k) ) &
                            + a7(2) * ( Frhou(i+2,j,k)-Frhou(i-2,j,k) ) &
                            + a7(3) * ( Frhou(i+3,j,k)-Frhou(i-3,j,k) ) ) *idx(i)

             Krhov(i,j,k) = ( a7(1) * ( Frhov(i+1,j,k)-Frhov(i-1,j,k) ) &
                            + a7(2) * ( Frhov(i+2,j,k)-Frhov(i-2,j,k) ) &
                            + a7(3) * ( Frhov(i+3,j,k)-Frhov(i-3,j,k) ) ) *idx(i)

             Krhow(i,j,k) = ( a7(1) * ( Frhow(i+1,j,k)-Frhow(i-1,j,k) ) &
                            + a7(2) * ( Frhow(i+2,j,k)-Frhow(i-2,j,k) ) &
                            + a7(3) * ( Frhow(i+3,j,k)-Frhow(i-3,j,k) ) ) *idx(i)

             Krhoe(i,j,k) = ( a7(1) * ( Frhoe(i+1,j,k)-Frhoe(i-1,j,k) ) &
                            + a7(2) * ( Frhoe(i+2,j,k)-Frhoe(i-2,j,k) ) &
                            + a7(3) * ( Frhoe(i+3,j,k)-Frhoe(i-3,j,k) ) ) *idx(i)
          enddo
       enddo
    enddo

    ! BC at imax
    ! ----------
    if (is_bc_1pt(1,2)) then
       i=nx-2
       do k=ndz_e,nfz_e
          do j=ndy_e,nfy_e
             Krho(i,j,k)  = ( a5(1)*( Frho(i+1,j,k)- Frho(i-1,j,k)) &
                            + a5(2)*( Frho(i+2,j,k)- Frho(i-2,j,k)) )*idx4_imax
             Krhou(i,j,k) = ( a5(1)*(Frhou(i+1,j,k)-Frhou(i-1,j,k)) &
                            + a5(2)*(Frhou(i+2,j,k)-Frhou(i-2,j,k)) )*idx4_imax
             Krhov(i,j,k) = ( a5(1)*(Frhov(i+1,j,k)-Frhov(i-1,j,k)) &
                            + a5(2)*(Frhov(i+2,j,k)-Frhov(i-2,j,k)) )*idx4_imax
             Krhow(i,j,k) = ( a5(1)*(Frhow(i+1,j,k)-Frhow(i-1,j,k)) &
                            + a5(2)*(Frhow(i+2,j,k)-Frhow(i-2,j,k)) )*idx4_imax
             Krhoe(i,j,k) = ( a5(1)*(Frhoe(i+1,j,k)-Frhoe(i-1,j,k)) &
                            + a5(2)*(Frhoe(i+2,j,k)-Frhoe(i-2,j,k)) )*idx4_imax
          enddo
       enddo

       i=nx-1
       do k=ndz_e,nfz_e
          do j=ndy_e,nfy_e
             Krho(i,j,k)  = a3(1)*( Frho(i+1,j,k)- Frho(i-1,j,k))*idx2_imax
             Krhou(i,j,k) = a3(1)*(Frhou(i+1,j,k)-Frhou(i-1,j,k))*idx2_imax
             Krhov(i,j,k) = a3(1)*(Frhov(i+1,j,k)-Frhov(i-1,j,k))*idx2_imax
             Krhow(i,j,k) = a3(1)*(Frhow(i+1,j,k)-Frhow(i-1,j,k))*idx2_imax
             Krhoe(i,j,k) = a3(1)*(Frhoe(i+1,j,k)-Frhoe(i-1,j,k))*idx2_imax
          enddo
       enddo
    endif

    ! Advancement of wall points
    ! --------------------------
    if (is_bc_wall2(1,2)) then
       i=nx
       do k=ndz_e,nfz_e
          do j=ndy_e,nfy_e
             Krho(i,j,k) =(a20(3)*Frho(i-2,j,k) +a20(2)*Frho(i-1,j,k) +a20(1)*Frho(i,j,k) )*idx1_imax
             Krhou(i,j,k)=(a20(3)*Frhou(i-2,j,k)+a20(2)*Frhou(i-1,j,k)+a20(1)*Frhou(i,j,k))*idx1_imax
             Krhov(i,j,k)=(a20(3)*Frhov(i-2,j,k)+a20(2)*Frhov(i-1,j,k)+a20(1)*Frhov(i,j,k))*idx1_imax
             Krhow(i,j,k)=(a20(3)*Frhow(i-2,j,k)+a20(2)*Frhow(i-1,j,k)+a20(1)*Frhow(i,j,k))*idx1_imax
             Krhoe(i,j,k)=(a20(3)*Frhoe(i-2,j,k)+a20(2)*Frhoe(i-1,j,k)+a20(1)*Frhoe(i,j,k))*idx1_imax
          enddo
       enddo
    endif

    ! Computation of inviscid fluxes along y
    ! ======================================
    do k=ndz_e,nfz_e
       do j=ndyt,nfyt
          do i=ndx_e,nfx_e
             Frho(i,j,k) =  rhov_n(i,j,k)
             Frhou(i,j,k)=  rhov_n(i,j,k)*uu(i,j,k)
             Frhov(i,j,k)=  rhov_n(i,j,k)*vv(i,j,k) +prs(i,j,k)
             Frhow(i,j,k)=  rhov_n(i,j,k)*ww(i,j,k)
             Frhoe(i,j,k)= (rhoe_n(i,j,k)+prs(i,j,k))*vv(i,j,k)
          enddo
       enddo
    enddo

    ! Flux derivatives along y
    ! ========================

    ! Advancement of wall points
    ! --------------------------
    if (is_bc_wall2(2,1)) then
       !j=1
       do k=ndz_e,nfz_e
          do i=ndx_e,nfx_e
             Krho(i,1,k)= (a02(1)*Frho(i,1,k) +a02(2)*Frho(i,2,k) +a02(3)*Frho(i,3,k) )*idy1_jmin+ Krho(i,1,k)
             Krhou(i,1,k)=(a02(1)*Frhou(i,1,k)+a02(2)*Frhou(i,2,k)+a02(3)*Frhou(i,3,k))*idy1_jmin+ Krhou(i,1,k)
             Krhov(i,1,k)=(a02(1)*Frhov(i,1,k)+a02(2)*Frhov(i,2,k)+a02(3)*Frhov(i,3,k))*idy1_jmin+ Krhov(i,1,k)
             Krhow(i,1,k)=(a02(1)*Frhow(i,1,k)+a02(2)*Frhow(i,2,k)+a02(3)*Frhow(i,3,k))*idy1_jmin+ Krhow(i,1,k)
             Krhoe(i,1,k)=(a02(1)*Frhoe(i,1,k)+a02(2)*Frhoe(i,2,k)+a02(3)*Frhoe(i,3,k))*idy1_jmin+ Krhoe(i,1,k)
          enddo
       enddo
    endif

    ! BC at jmin
    ! ----------
    if (is_bc_1pt(2,1)) then     
       j=2
       do k=ndz_e,nfz_e
          do i=ndx_e,nfx_e
             Krho(i,j,k)  =  Krho(i,j,k) + a3(1)*( Frho(i,j+1,k)- Frho(i,j-1,k))*idy2_jmin
             Krhou(i,j,k) = Krhou(i,j,k) + a3(1)*(Frhou(i,j+1,k)-Frhou(i,j-1,k))*idy2_jmin
             Krhov(i,j,k) = Krhov(i,j,k) + a3(1)*(Frhov(i,j+1,k)-Frhov(i,j-1,k))*idy2_jmin
             Krhow(i,j,k) = Krhow(i,j,k) + a3(1)*(Frhow(i,j+1,k)-Frhow(i,j-1,k))*idy2_jmin
             Krhoe(i,j,k) = Krhoe(i,j,k) + a3(1)*(Frhoe(i,j+1,k)-Frhoe(i,j-1,k))*idy2_jmin
          enddo
       enddo

       j=3
       do k=ndz_e,nfz_e
          do i=ndx_e,nfx_e
             Krho(i,j,k)  = ( a5(1)*(Frho(i,j+1,k) - Frho(i,j-1,k)) &
                            + a5(2)*(Frho(i,j+2,k) - Frho(i,j-2,k)) )*idy4_jmin + Krho(i,j,k)
             Krhou(i,j,k) = ( a5(1)*(Frhou(i,j+1,k)-Frhou(i,j-1,k)) &
                            + a5(2)*(Frhou(i,j+2,k)-Frhou(i,j-2,k)) )*idy4_jmin + Krhou(i,j,k)
             Krhov(i,j,k) = ( a5(1)*(Frhov(i,j+1,k)-Frhov(i,j-1,k)) &
                            + a5(2)*(Frhov(i,j+2,k)-Frhov(i,j-2,k)) )*idy4_jmin + Krhov(i,j,k)
             Krhow(i,j,k) = ( a5(1)*(Frhow(i,j+1,k)-Frhow(i,j-1,k)) &
                            + a5(2)*(Frhow(i,j+2,k)-Frhow(i,j-2,k)) )*idy4_jmin + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a5(1)*(Frhoe(i,j+1,k)-Frhoe(i,j-1,k)) &
                            + a5(2)*(Frhoe(i,j+2,k)-Frhoe(i,j-2,k)) )*idy4_jmin + Krhoe(i,j,k)
          enddo
       enddo
    endif

    ! Interior points
    ! ---------------
    do k=ndz_e,nfz_e
       do j=ndy,nfy
          do i=ndx_e,nfx_e
             Krho(i,j,k)  = ( a7(1) * ( Frho(i,j+1,k)-Frho(i,j-1,k) ) &
                            + a7(2) * ( Frho(i,j+2,k)-Frho(i,j-2,k) ) &
                            + a7(3) * ( Frho(i,j+3,k)-Frho(i,j-3,k) ) ) *idy(j) &
                            + Krho(i,j,k)

             Krhou(i,j,k) = ( a7(1) * ( Frhou(i,j+1,k)-Frhou(i,j-1,k) ) &
                            + a7(2) * ( Frhou(i,j+2,k)-Frhou(i,j-2,k) ) &
                            + a7(3) * ( Frhou(i,j+3,k)-Frhou(i,j-3,k) ) ) *idy(j) &
                            + Krhou(i,j,k)

             Krhov(i,j,k) = ( a7(1) * ( Frhov(i,j+1,k)-Frhov(i,j-1,k) ) &
                            + a7(2) * ( Frhov(i,j+2,k)-Frhov(i,j-2,k) ) &
                            + a7(3) * ( Frhov(i,j+3,k)-Frhov(i,j-3,k) ) ) *idy(j) &
                            + Krhov(i,j,k)

             Krhow(i,j,k) = ( a7(1) * ( Frhow(i,j+1,k)-Frhow(i,j-1,k) ) &
                            + a7(2) * ( Frhow(i,j+2,k)-Frhow(i,j-2,k) ) &
                            + a7(3) * ( Frhow(i,j+3,k)-Frhow(i,j-3,k) ) ) *idy(j) &
                            + Krhow(i,j,k)

             Krhoe(i,j,k) = ( a7(1) * ( Frhoe(i,j+1,k)-Frhoe(i,j-1,k) ) &
                            + a7(2) * ( Frhoe(i,j+2,k)-Frhoe(i,j-2,k) ) &
                            + a7(3) * ( Frhoe(i,j+3,k)-Frhoe(i,j-3,k) ) ) *idy(j) &
                            + Krhoe(i,j,k)
          enddo
       enddo
    enddo

    ! BC at jmax
    ! ----------
    if (is_bc_1pt(2,2)) then
       j=ny-2
       do k=ndz_e,nfz_e
          do i=ndx_e,nfx_e
             Krho(i,j,k)  = ( a5(1)*(Frho(i,j+1,k) - Frho(i,j-1,k))  &
                            + a5(2)*(Frho(i,j+2,k) - Frho(i,j-2,k)) )*idy4_jmax + Krho(i,j,k)
             Krhou(i,j,k) = ( a5(1)*(Frhou(i,j+1,k)-Frhou(i,j-1,k))  &
                            + a5(2)*(Frhou(i,j+2,k)-Frhou(i,j-2,k)) )*idy4_jmax + Krhou(i,j,k)
             Krhov(i,j,k) = ( a5(1)*(Frhov(i,j+1,k)-Frhov(i,j-1,k))  &
                            + a5(2)*(Frhov(i,j+2,k)-Frhov(i,j-2,k)) )*idy4_jmax + Krhov(i,j,k)
             Krhow(i,j,k) = ( a5(1)*(Frhow(i,j+1,k)-Frhow(i,j-1,k))  &
                            + a5(2)*(Frhow(i,j+2,k)-Frhow(i,j-2,k)) )*idy4_jmax + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a5(1)*(Frhoe(i,j+1,k)-Frhoe(i,j-1,k))  &
                            + a5(2)*(Frhoe(i,j+2,k)-Frhoe(i,j-2,k)) )*idy4_jmax + Krhoe(i,j,k)
          enddo
       enddo

       j=ny-1
       do k=ndz_e,nfz_e
          do i=ndx_e,nfx_e
             Krho(i,j,k)  =  Krho(i,j,k) + a3(1)*( Frho(i,j+1,k)- Frho(i,j-1,k))*idy2_jmax
             Krhou(i,j,k) = Krhou(i,j,k) + a3(1)*(Frhou(i,j+1,k)-Frhou(i,j-1,k))*idy2_jmax
             Krhov(i,j,k) = Krhov(i,j,k) + a3(1)*(Frhov(i,j+1,k)-Frhov(i,j-1,k))*idy2_jmax
             Krhow(i,j,k) = Krhow(i,j,k) + a3(1)*(Frhow(i,j+1,k)-Frhow(i,j-1,k))*idy2_jmax
             Krhoe(i,j,k) = Krhoe(i,j,k) + a3(1)*(Frhoe(i,j+1,k)-Frhoe(i,j-1,k))*idy2_jmax
          enddo
       enddo
    endif

    ! Advancement of wall points
    ! --------------------------
    if (is_bc_wall2(2,2)) then
       j=ny
       do k=ndz_e,nfz_e
          do i=ndx_e,nfx_e
             Krho(i,j,k) =(a20(3)*Frho(i,j-2,k) +a20(2)*Frho(i,j-1,k) +a20(1)*Frho(i,j,k) )*idy1_jmax+ Krho(i,j,k)
             Krhou(i,j,k)=(a20(3)*Frhou(i,j-2,k)+a20(2)*Frhou(i,j-1,k)+a20(1)*Frhou(i,j,k))*idy1_jmax+ Krhou(i,j,k)
             Krhov(i,j,k)=(a20(3)*Frhov(i,j-2,k)+a20(2)*Frhov(i,j-1,k)+a20(1)*Frhov(i,j,k))*idy1_jmax+ Krhov(i,j,k)
             Krhow(i,j,k)=(a20(3)*Frhow(i,j-2,k)+a20(2)*Frhow(i,j-1,k)+a20(1)*Frhow(i,j,k))*idy1_jmax+ Krhow(i,j,k)
             Krhoe(i,j,k)=(a20(3)*Frhoe(i,j-2,k)+a20(2)*Frhoe(i,j-1,k)+a20(1)*Frhoe(i,j,k))*idy1_jmax+ Krhoe(i,j,k)
          enddo
       enddo
    endif

    !****************
    if (is_2D) return
    !****************

    ! Computation of inviscid fluxes along z
    ! ======================================
    do k=ndzt,nfzt
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Frho(i,j,k) =  rhow_n(i,j,k)
             Frhou(i,j,k)=  rhow_n(i,j,k)*uu(i,j,k)
             Frhov(i,j,k)=  rhow_n(i,j,k)*vv(i,j,k)
             Frhow(i,j,k)=  rhow_n(i,j,k)*ww(i,j,k) +prs(i,j,k)
             Frhoe(i,j,k)= (rhoe_n(i,j,k)+prs(i,j,k))*ww(i,j,k)
          enddo
       enddo
    enddo

    ! Flux derivatives along z
    ! ========================

    ! Advancement of wall points
    ! --------------------------
    if (is_bc_wall2(3,1)) then
       !k=1
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Krho(i,j,1)= (a02(1)*Frho(i,j,1) +a02(2)*Frho(i,j,2) +a02(3)*Frho(i,j,3) )*idz1_kmin+ Krho(i,j,1)
             Krhou(i,j,1)=(a02(1)*Frhou(i,j,1)+a02(2)*Frhou(i,j,2)+a02(3)*Frhou(i,j,3))*idz1_kmin+ Krhou(i,j,1)
             Krhov(i,j,1)=(a02(1)*Frhov(i,j,1)+a02(2)*Frhov(i,j,2)+a02(3)*Frhov(i,j,3))*idz1_kmin+ Krhov(i,j,1)
             Krhow(i,j,1)=(a02(1)*Frhow(i,j,1)+a02(2)*Frhow(i,j,2)+a02(3)*Frhow(i,j,3))*idz1_kmin+ Krhow(i,j,1)
             Krhoe(i,j,1)=(a02(1)*Frhoe(i,j,1)+a02(2)*Frhoe(i,j,2)+a02(3)*Frhoe(i,j,3))*idz1_kmin+ Krhoe(i,j,1)
          enddo
       enddo
    endif

    ! BC at kmin
    ! ----------
    if (is_bc_1pt(3,1)) then
       k=2
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Krho(i,j,k)  =  Krho(i,j,k) + ( Frho(i,j,k+1)- Frho(i,j,k-1))*idz2_kmin
             Krhou(i,j,k) = Krhou(i,j,k) + (Frhou(i,j,k+1)-Frhou(i,j,k-1))*idz2_kmin
             Krhov(i,j,k) = Krhov(i,j,k) + (Frhov(i,j,k+1)-Frhov(i,j,k-1))*idz2_kmin
             Krhow(i,j,k) = Krhow(i,j,k) + (Frhow(i,j,k+1)-Frhow(i,j,k-1))*idz2_kmin
             Krhoe(i,j,k) = Krhoe(i,j,k) + (Frhoe(i,j,k+1)-Frhoe(i,j,k-1))*idz2_kmin
          enddo
       enddo

       k=3
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Krho(i,j,k)  = ( a5(1)*( Frho(i,j,k+1)- Frho(i,j,k-1)) &
                            + a5(2)*( Frho(i,j,k+2)- Frho(i,j,k-2)) )*idz4_kmin + Krho(i,j,k)
             Krhou(i,j,k) = ( a5(1)*(Frhou(i,j,k+1)-Frhou(i,j,k-1)) &
                            + a5(2)*(Frhou(i,j,k+2)-Frhou(i,j,k-2)) )*idz4_kmin + Krhou(i,j,k)
             Krhov(i,j,k) = ( a5(1)*(Frhov(i,j,k+1)-Frhov(i,j,k-1)) &
                            + a5(2)*(Frhov(i,j,k+2)-Frhov(i,j,k-2)) )*idz4_kmin + Krhov(i,j,k)
             Krhow(i,j,k) = ( a5(1)*(Frhow(i,j,k+1)-Frhow(i,j,k-1)) &
                            + a5(2)*(Frhow(i,j,k+2)-Frhow(i,j,k-2)) )*idz4_kmin + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a5(1)*(Frhoe(i,j,k+1)-Frhoe(i,j,k-1)) &
                            + a5(2)*(Frhoe(i,j,k+2)-Frhoe(i,j,k-2)) )*idz4_kmin + Krhoe(i,j,k)
          enddo
       enddo
    endif

    ! Interior points
    ! ---------------
    do k=ndz,nfz
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Krho(i,j,k)  = ( a7(1)* ( Frho(i,j,k+1)-Frho(i,j,k-1) ) &
                            + a7(2)* ( Frho(i,j,k+2)-Frho(i,j,k-2) ) &
                            + a7(3)* ( Frho(i,j,k+3)-Frho(i,j,k-3) ) ) *idz(k) &
                            + Krho(i,j,k)

             Krhou(i,j,k) = ( a7(1)* ( Frhou(i,j,k+1)-Frhou(i,j,k-1) ) &
                            + a7(2)* ( Frhou(i,j,k+2)-Frhou(i,j,k-2) ) &
                            + a7(3)* ( Frhou(i,j,k+3)-Frhou(i,j,k-3) ) ) *idz(k) &
                            + Krhou(i,j,k)

             Krhov(i,j,k) = ( a7(1)* ( Frhov(i,j,k+1)-Frhov(i,j,k-1) ) &
                            + a7(2)* ( Frhov(i,j,k+2)-Frhov(i,j,k-2) ) &
                            + a7(3)* ( Frhov(i,j,k+3)-Frhov(i,j,k-3) ) ) *idz(k) &
                            + Krhov(i,j,k)

             Krhow(i,j,k) = ( a7(1)* ( Frhow(i,j,k+1)-Frhow(i,j,k-1) ) &
                            + a7(2)* ( Frhow(i,j,k+2)-Frhow(i,j,k-2) ) &
                            + a7(3)* ( Frhow(i,j,k+3)-Frhow(i,j,k-3) ) ) *idz(k) &
                            + Krhow(i,j,k)

             Krhoe(i,j,k) = ( a7(1)* ( Frhoe(i,j,k+1)-Frhoe(i,j,k-1) ) &
                            + a7(2)* ( Frhoe(i,j,k+2)-Frhoe(i,j,k-2) ) &
                            + a7(3)* ( Frhoe(i,j,k+3)-Frhoe(i,j,k-3) ) ) *idz(k) &
                            + Krhoe(i,j,k)
          enddo
       enddo
    enddo

    ! BC at kmax
    ! ----------
    if (is_bc_1pt(3,2)) then
       k=nz-2
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Krho(i,j,k)  = ( a5(1)*( Frho(i,j,k+1)- Frho(i,j,k-1)) &
                            + a5(2)*( Frho(i,j,k+2)- Frho(i,j,k-2)) )*idz4_kmax + Krho(i,j,k)
             Krhou(i,j,k) = ( a5(1)*(Frhou(i,j,k+1)-Frhou(i,j,k-1)) &
                            + a5(2)*(Frhou(i,j,k+2)-Frhou(i,j,k-2)) )*idz4_kmax + Krhou(i,j,k)
             Krhov(i,j,k) = ( a5(1)*(Frhov(i,j,k+1)-Frhov(i,j,k-1)) &
                            + a5(2)*(Frhov(i,j,k+2)-Frhov(i,j,k-2)) )*idz4_kmax + Krhov(i,j,k)
             Krhow(i,j,k) = ( a5(1)*(Frhow(i,j,k+1)-Frhow(i,j,k-1)) &
                            + a5(2)*(Frhow(i,j,k+2)-Frhow(i,j,k-2)) )*idz4_kmax + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a5(1)*(Frhoe(i,j,k+1)-Frhoe(i,j,k-1)) &
                            + a5(2)*(Frhoe(i,j,k+2)-Frhoe(i,j,k-2)) )*idz4_kmax + Krhoe(i,j,k)
          enddo
       enddo

       k=nz-1
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Krho(i,j,k)  = Krho(i,j,k)  + ( Frho(i,j,k+1)- Frho(i,j,k-1))*idz2_kmax
             Krhou(i,j,k) = Krhou(i,j,k) + (Frhou(i,j,k+1)-Frhou(i,j,k-1))*idz2_kmax
             Krhov(i,j,k) = Krhov(i,j,k) + (Frhov(i,j,k+1)-Frhov(i,j,k-1))*idz2_kmax
             Krhow(i,j,k) = Krhow(i,j,k) + (Frhow(i,j,k+1)-Frhow(i,j,k-1))*idz2_kmax
             Krhoe(i,j,k) = Krhoe(i,j,k) + (Frhoe(i,j,k+1)-Frhoe(i,j,k-1))*idz2_kmax
          enddo
       enddo
    endif

    ! Advancement of wall points
    ! --------------------------
    if (is_bc_wall2(3,2)) then
       k=nz
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Krho(i,j,k) =(a20(3)*Frho(i,j,k-2) +a20(2)*Frho(i,j,k-1) +a20(1)*Frho(i,j,k) )*idz1_kmax+ Krho(i,j,k)
             Krhou(i,j,k)=(a20(3)*Frhou(i,j,k-2)+a20(2)*Frhou(i,j,k-1)+a20(1)*Frhou(i,j,k))*idz1_kmax+ Krhou(i,j,k)
             Krhov(i,j,k)=(a20(3)*Frhov(i,j,k-2)+a20(2)*Frhov(i,j,k-1)+a20(1)*Frhov(i,j,k))*idz1_kmax+ Krhov(i,j,k)
             Krhow(i,j,k)=(a20(3)*Frhow(i,j,k-2)+a20(2)*Frhow(i,j,k-1)+a20(1)*Frhow(i,j,k))*idz1_kmax+ Krhow(i,j,k)
             Krhoe(i,j,k)=(a20(3)*Frhoe(i,j,k-2)+a20(2)*Frhoe(i,j,k-1)+a20(1)*Frhoe(i,j,k))*idz1_kmax+ Krhoe(i,j,k)
          enddo
       enddo
    endif

  end subroutine flux_euler_7pts

end submodule smod_flux_euler_7pts
