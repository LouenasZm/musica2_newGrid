!=================================================================================
submodule (mod_flux_euler) smod_flux_euler_11pts
!=================================================================================
  !> author: XG
  !> date: January 2024
  !> Cartesian version - compute Eulerian fluxes (inviscid part) with 11-pt stencil
!=================================================================================

contains

  !==============================================================================
  module subroutine flux_euler_11pts
  !==============================================================================
    !> Derivatives of Eulerian fluxes (inviscid part) - 11-point stencil -
    !> - Cartesian version -
  !==============================================================================
    !use mod_fluid
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    ! ---------------------------------------------------------------------------
    !! real(wp) :: drdy,dudy,dvdy,dpdy,vq ! XG DEV
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

       i=4
       do k=ndz_e,nfz_e
          do j=ndy_e,nfy_e
             Krho(i,j,k)  = ( a7(1) * ( Frho(i+1,j,k)-Frho(i-1,j,k) ) &
                            + a7(2) * ( Frho(i+2,j,k)-Frho(i-2,j,k) ) &
                            + a7(3) * ( Frho(i+3,j,k)-Frho(i-3,j,k) ) )*idx6_imin
             Krhou(i,j,k) = ( a7(1) * ( Frhou(i+1,j,k)-Frhou(i-1,j,k) ) &
                            + a7(2) * ( Frhou(i+2,j,k)-Frhou(i-2,j,k) ) &
                            + a7(3) * ( Frhou(i+3,j,k)-Frhou(i-3,j,k) ) )*idx6_imin
             Krhov(i,j,k) = ( a7(1) * ( Frhov(i+1,j,k)-Frhov(i-1,j,k) ) &
                            + a7(2) * ( Frhov(i+2,j,k)-Frhov(i-2,j,k) ) &
                            + a7(3) * ( Frhov(i+3,j,k)-Frhov(i-3,j,k) ) )*idx6_imin
             Krhow(i,j,k) = ( a7(1) * ( Frhow(i+1,j,k)-Frhow(i-1,j,k) ) &
                            + a7(2) * ( Frhow(i+2,j,k)-Frhow(i-2,j,k) ) &
                            + a7(3) * ( Frhow(i+3,j,k)-Frhow(i-3,j,k) ) )*idx6_imin
             Krhoe(i,j,k) = ( a7(1) * ( Frhoe(i+1,j,k)-Frhoe(i-1,j,k) ) &
                            + a7(2) * ( Frhoe(i+2,j,k)-Frhoe(i-2,j,k) ) &
                            + a7(3) * ( Frhoe(i+3,j,k)-Frhoe(i-3,j,k) ) )*idx6_imin
          enddo
       enddo

       i=5
       do k=ndz_e,nfz_e
          do j=ndy_e,nfy_e
             Krho(i,j,k)  = ( a9(1) * ( Frho(i+1,j,k)-Frho(i-1,j,k) ) &
                            + a9(2) * ( Frho(i+2,j,k)-Frho(i-2,j,k) ) &
                            + a9(3) * ( Frho(i+3,j,k)-Frho(i-3,j,k) ) &
                            + a9(4) * ( Frho(i+4,j,k)-Frho(i-4,j,k) ) )*idx8_imin
             Krhou(i,j,k) = ( a9(1) * ( Frhou(i+1,j,k)-Frhou(i-1,j,k) ) &
                            + a9(2) * ( Frhou(i+2,j,k)-Frhou(i-2,j,k) ) &
                            + a9(3) * ( Frhou(i+3,j,k)-Frhou(i-3,j,k) ) &
                            + a9(4) * ( Frhou(i+4,j,k)-Frhou(i-4,j,k) ) )*idx8_imin
             Krhov(i,j,k) = ( a9(1) * ( Frhov(i+1,j,k)-Frhov(i-1,j,k) ) &
                            + a9(2) * ( Frhov(i+2,j,k)-Frhov(i-2,j,k) ) &
                            + a9(3) * ( Frhov(i+3,j,k)-Frhov(i-3,j,k) ) &
                            + a9(4) * ( Frhov(i+4,j,k)-Frhov(i-4,j,k) ) )*idx8_imin
             Krhow(i,j,k) = ( a9(1) * ( Frhow(i+1,j,k)-Frhow(i-1,j,k) ) &
                            + a9(2) * ( Frhow(i+2,j,k)-Frhow(i-2,j,k) ) &
                            + a9(3) * ( Frhow(i+3,j,k)-Frhow(i-3,j,k) ) &
                            + a9(4) * ( Frhow(i+4,j,k)-Frhow(i-4,j,k) ) )*idx8_imin
             Krhoe(i,j,k) = ( a9(1) * ( Frhoe(i+1,j,k)-Frhoe(i-1,j,k) ) &
                            + a9(2) * ( Frhoe(i+2,j,k)-Frhoe(i-2,j,k) ) &
                            + a9(3) * ( Frhoe(i+3,j,k)-Frhoe(i-3,j,k) ) &
                            + a9(4) * ( Frhoe(i+4,j,k)-Frhoe(i-4,j,k) ) )*idx8_imin
          enddo
       enddo
    endif

    ! Interior points
    ! ---------------
    do k=ndz_e,nfz_e
       do j=ndy_e,nfy_e
          do i=ndx,nfx

             Krho(i,j,k)  = ( a11(1) * ( Frho(i+1,j,k)-Frho(i-1,j,k) )   &
                            + a11(2) * ( Frho(i+2,j,k)-Frho(i-2,j,k) )   &
                            + a11(3) * ( Frho(i+3,j,k)-Frho(i-3,j,k) )   &
                            + a11(4) * ( Frho(i+4,j,k)-Frho(i-4,j,k) )   &
                            + a11(5) * ( Frho(i+5,j,k)-Frho(i-5,j,k) ) )*idx(i)

             Krhou(i,j,k) = ( a11(1) * ( Frhou(i+1,j,k)-Frhou(i-1,j,k) ) &
                            + a11(2) * ( Frhou(i+2,j,k)-Frhou(i-2,j,k) ) &
                            + a11(3) * ( Frhou(i+3,j,k)-Frhou(i-3,j,k) ) &
                            + a11(4) * ( Frhou(i+4,j,k)-Frhou(i-4,j,k) ) &
                            + a11(5) * ( Frhou(i+5,j,k)-Frhou(i-5,j,k) ) )*idx(i)

             Krhov(i,j,k) = ( a11(1) * ( Frhov(i+1,j,k)-Frhov(i-1,j,k) ) &
                            + a11(2) * ( Frhov(i+2,j,k)-Frhov(i-2,j,k) ) &
                            + a11(3) * ( Frhov(i+3,j,k)-Frhov(i-3,j,k) ) &
                            + a11(4) * ( Frhov(i+4,j,k)-Frhov(i-4,j,k) ) &
                            + a11(5) * ( Frhov(i+5,j,k)-Frhov(i-5,j,k) ) )*idx(i)

             Krhow(i,j,k) = ( a11(1) * ( Frhow(i+1,j,k)-Frhow(i-1,j,k) ) &
                            + a11(2) * ( Frhow(i+2,j,k)-Frhow(i-2,j,k) ) &
                            + a11(3) * ( Frhow(i+3,j,k)-Frhow(i-3,j,k) ) &
                            + a11(4) * ( Frhow(i+4,j,k)-Frhow(i-4,j,k) ) &
                            + a11(5) * ( Frhow(i+5,j,k)-Frhow(i-5,j,k) ) )*idx(i)

             Krhoe(i,j,k) = ( a11(1) * ( Frhoe(i+1,j,k)-Frhoe(i-1,j,k) ) &
                            + a11(2) * ( Frhoe(i+2,j,k)-Frhoe(i-2,j,k) ) &
                            + a11(3) * ( Frhoe(i+3,j,k)-Frhoe(i-3,j,k) ) &
                            + a11(4) * ( Frhoe(i+4,j,k)-Frhoe(i-4,j,k) ) &
                            + a11(5) * ( Frhoe(i+5,j,k)-Frhoe(i-5,j,k) ) )*idx(i)
          enddo
       enddo
    enddo

    ! BC at imax
    ! ----------
    if (is_bc_1pt(1,2)) then
       i=nx-4
       do k=ndz_e,nfz_e
          do j=ndy_e,nfy_e
             Krho(i,j,k)  = ( a9(1) * ( Frho(i+1,j,k)-Frho(i-1,j,k) ) &
                            + a9(2) * ( Frho(i+2,j,k)-Frho(i-2,j,k) ) &
                            + a9(3) * ( Frho(i+3,j,k)-Frho(i-3,j,k) ) &
                            + a9(4) * ( Frho(i+4,j,k)-Frho(i-4,j,k) ) )*idx8_imax
             Krhou(i,j,k) = ( a9(1) * ( Frhou(i+1,j,k)-Frhou(i-1,j,k) ) &
                            + a9(2) * ( Frhou(i+2,j,k)-Frhou(i-2,j,k) ) &
                            + a9(3) * ( Frhou(i+3,j,k)-Frhou(i-3,j,k) ) &
                            + a9(4) * ( Frhou(i+4,j,k)-Frhou(i-4,j,k) ) )*idx8_imax
             Krhov(i,j,k) = ( a9(1) * ( Frhov(i+1,j,k)-Frhov(i-1,j,k) ) &
                            + a9(2) * ( Frhov(i+2,j,k)-Frhov(i-2,j,k) ) &
                            + a9(3) * ( Frhov(i+3,j,k)-Frhov(i-3,j,k) ) &
                            + a9(4) * ( Frhov(i+4,j,k)-Frhov(i-4,j,k) ) )*idx8_imax
             Krhow(i,j,k) = ( a9(1) * ( Frhow(i+1,j,k)-Frhow(i-1,j,k) ) &
                            + a9(2) * ( Frhow(i+2,j,k)-Frhow(i-2,j,k) ) &
                            + a9(3) * ( Frhow(i+3,j,k)-Frhow(i-3,j,k) ) &
                            + a9(4) * ( Frhow(i+4,j,k)-Frhow(i-4,j,k) ) )*idx8_imax
             Krhoe(i,j,k) = ( a9(1) * ( Frhoe(i+1,j,k)-Frhoe(i-1,j,k) ) &
                            + a9(2) * ( Frhoe(i+2,j,k)-Frhoe(i-2,j,k) ) &
                            + a9(3) * ( Frhoe(i+3,j,k)-Frhoe(i-3,j,k) ) &
                            + a9(4) * ( Frhoe(i+4,j,k)-Frhoe(i-4,j,k) ) )*idx8_imax
          enddo
       enddo

       i=nx-3
       do k=ndz_e,nfz_e
          do j=ndy_e,nfy_e
             Krho(i,j,k)  = ( a7(1) * ( Frho(i+1,j,k)-Frho(i-1,j,k) ) &
                            + a7(2) * ( Frho(i+2,j,k)-Frho(i-2,j,k) ) &
                            + a7(3) * ( Frho(i+3,j,k)-Frho(i-3,j,k) ) )*idx6_imax
             Krhou(i,j,k) = ( a7(1) * ( Frhou(i+1,j,k)-Frhou(i-1,j,k) ) &
                            + a7(2) * ( Frhou(i+2,j,k)-Frhou(i-2,j,k) ) &
                            + a7(3) * ( Frhou(i+3,j,k)-Frhou(i-3,j,k) ) )*idx6_imax
             Krhov(i,j,k) = ( a7(1) * ( Frhov(i+1,j,k)-Frhov(i-1,j,k) ) &
                            + a7(2) * ( Frhov(i+2,j,k)-Frhov(i-2,j,k) ) &
                            + a7(3) * ( Frhov(i+3,j,k)-Frhov(i-3,j,k) ) )*idx6_imax
             Krhow(i,j,k) = ( a7(1) * ( Frhow(i+1,j,k)-Frhow(i-1,j,k) ) &
                            + a7(2) * ( Frhow(i+2,j,k)-Frhow(i-2,j,k) ) &
                            + a7(3) * ( Frhow(i+3,j,k)-Frhow(i-3,j,k) ) )*idx6_imax
             Krhoe(i,j,k) = ( a7(1) * ( Frhoe(i+1,j,k)-Frhoe(i-1,j,k) ) &
                            + a7(2) * ( Frhoe(i+2,j,k)-Frhoe(i-2,j,k) ) &
                            + a7(3) * ( Frhoe(i+3,j,k)-Frhoe(i-3,j,k) ) )*idx6_imax
          enddo
       enddo

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
       ! XG DEV
!!$       ! simplified
!!$       j=1
!!$       do k=ndz_e,nfz_e
!!$          do i=ndx_e,nfx_e          
!!$             dvdy= (as4p0(1)*vv(i,1,k)+as4p0(2)*vv(i,2,k) &
!!$                  + as4p0(3)*vv(i,3,k)+as4p0(4)*vv(i,4,k))*idy1_jmin
!!$
!!$             !!          dvdy=(a04(2)*vv(i,j+1,k)+a04(3)*vv(i,j+2,k) &
!!$             !!               +a04(4)*vv(i,j+3,k)+a04(5)*vv(i,j+4,k))*idy1_jmin
!!$
!!$             !!           dvdy= (a06(1)*vv(i,1,k)+a06(2)*vv(i,2,k) &
!!$             !!                + a06(3)*vv(i,3,k)+a06(4)*vv(i,4,k) &
!!$             !!                + a06(5)*vv(i,5,k)+a06(6)*vv(i,6,k) &
!!$             !!                + a06(7)*vv(i,7,k))*idy1_jmin
!!$             Krho(i,j,k)=rho_n(i,j,k)*dvdy + Krho(i,j,k)
!!$             Krhou(i,j,k)=rhou_n(i,j,k)*dvdy + Krhou(i,j,k)
!!$             Krhoe(i,j,k)=(gam*igm1*prs(i,j,k)+0.5_wp*rhou_n(i,j,k)*uu(i,j,k))*dvdy &
!!$                  + Krhoe(i,j,k)
!!$          enddo
!!$       enddo

       j=2
       do k=ndz_e,nfz_e
          do i=ndx_e,nfx_e
             ! XG DEV
!!$             Krho(i,j,k)= (a15(1)*Frho(i,1,k)+a15(2)*Frho(i,2,k) &
!!$                         + a15(3)*Frho(i,3,k)+a15(4)*Frho(i,4,k) &
!!$                         + a15(5)*Frho(i,5,k)+a15(6)*Frho(i,6,k) &
!!$                         + a15(7)*Frho(i,7,k))*idy2_jmin+ Krho(i,j,k)
!!$             Krhou(i,j,k)= (a15(1)*Frhou(i,1,k)+a15(2)*Frhou(i,2,k) &
!!$                         + a15(3)*Frhou(i,3,k)+a15(4)*Frhou(i,4,k) &
!!$                         + a15(5)*Frhou(i,5,k)+a15(6)*Frhou(i,6,k) &
!!$                         + a15(7)*Frhou(i,7,k))*idy2_jmin+ Krhou(i,j,k)
!!$             Krhov(i,j,k)= (a15(1)*Frhov(i,1,k)+a15(2)*Frhov(i,2,k) &
!!$                         + a15(3)*Frhov(i,3,k)+a15(4)*Frhov(i,4,k) &
!!$                         + a15(5)*Frhov(i,5,k)+a15(6)*Frhov(i,6,k) &
!!$                         + a15(7)*Frhov(i,7,k))*idy2_jmin+ Krhov(i,j,k)
!!$             Krhow(i,j,k)= (a15(1)*Frhow(i,1,k)+a15(2)*Frhow(i,2,k) &
!!$                         + a15(3)*Frhow(i,3,k)+a15(4)*Frhow(i,4,k) &
!!$                         + a15(5)*Frhow(i,5,k)+a15(6)*Frhow(i,6,k) &
!!$                         + a15(7)*Frhow(i,7,k))*idy2_jmin+ Krhow(i,j,k)
!!$             Krhoe(i,j,k)= (a15(1)*Frhoe(i,1,k)+a15(2)*Frhoe(i,2,k) &
!!$                         + a15(3)*Frhoe(i,3,k)+a15(4)*Frhoe(i,4,k) &
!!$                         + a15(5)*Frhoe(i,5,k)+a15(6)*Frhoe(i,6,k) &
!!$                         + a15(7)*Frhoe(i,7,k))*idy2_jmin+ Krhoe(i,j,k)
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
             ! XG DEV
!!$             Krho(i,j,k)= (a24(1)*Frho(i,1,k)+a24(2)*Frho(i,2,k) &
!!$                         + a24(3)*Frho(i,3,k)+a24(4)*Frho(i,4,k) &
!!$                         + a24(5)*Frho(i,5,k)+a24(6)*Frho(i,6,k) &
!!$                         + a24(7)*Frho(i,7,k))*idy4_jmin + Krho(i,j,k)
!!$             Krhou(i,j,k)= (a24(1)*Frhou(i,1,k)+a24(2)*Frhou(i,2,k) &
!!$                         + a24(3)*Frhou(i,3,k)+a24(4)*Frhou(i,4,k) &
!!$                         + a24(5)*Frhou(i,5,k)+a24(6)*Frhou(i,6,k) &
!!$                         + a24(7)*Frhou(i,7,k))*idy4_jmin + Krhou(i,j,k)
!!$             Krhov(i,j,k)= (a24(1)*Frhov(i,1,k)+a24(2)*Frhov(i,2,k) &
!!$                         + a24(3)*Frhov(i,3,k)+a24(4)*Frhov(i,4,k) &
!!$                         + a24(5)*Frhov(i,5,k)+a24(6)*Frhov(i,6,k) &
!!$                         + a24(7)*Frhov(i,7,k))*idy4_jmin + Krhov(i,j,k)
!!$             Krhow(i,j,k)= (a24(1)*Frhow(i,1,k)+a24(2)*Frhow(i,2,k) &
!!$                         + a24(3)*Frhow(i,3,k)+a24(4)*Frhow(i,4,k) &
!!$                         + a24(5)*Frhow(i,5,k)+a24(6)*Frhow(i,6,k) &
!!$                         + a24(7)*Frhow(i,7,k))*idy4_jmin + Krhow(i,j,k)
!!$             Krhoe(i,j,k)= (a24(1)*Frhoe(i,1,k)+a24(2)*Frhoe(i,2,k) &
!!$                         + a24(3)*Frhoe(i,3,k)+a24(4)*Frhoe(i,4,k) &
!!$                         + a24(5)*Frhoe(i,5,k)+a24(6)*Frhoe(i,6,k) &
!!$                         + a24(7)*Frhoe(i,7,k))*idy4_jmin + Krhoe(i,j,k)
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

       j=4
       do k=ndz_e,nfz_e
          do i=ndx_e,nfx_e
             Krho(i,j,k)  = ( a7(1) * ( Frho(i,j+1,k)-Frho(i,j-1,k) ) &
                            + a7(2) * ( Frho(i,j+2,k)-Frho(i,j-2,k) ) &
                            + a7(3) * ( Frho(i,j+3,k)-Frho(i,j-3,k) ) )*idy6_jmin &
                            + Krho(i,j,k)
             Krhou(i,j,k) = ( a7(1) * ( Frhou(i,j+1,k)-Frhou(i,j-1,k) ) &
                            + a7(2) * ( Frhou(i,j+2,k)-Frhou(i,j-2,k) ) &
                            + a7(3) * ( Frhou(i,j+3,k)-Frhou(i,j-3,k) ) )*idy6_jmin &
                            + Krhou(i,j,k)
             Krhov(i,j,k) = ( a7(1) * ( Frhov(i,j+1,k)-Frhov(i,j-1,k) ) &
                            + a7(2) * ( Frhov(i,j+2,k)-Frhov(i,j-2,k) ) &
                            + a7(3) * ( Frhov(i,j+3,k)-Frhov(i,j-3,k) ) )*idy6_jmin &
                            + Krhov(i,j,k)
             Krhow(i,j,k) = ( a7(1) * ( Frhow(i,j+1,k)-Frhow(i,j-1,k) ) &
                            + a7(2) * ( Frhow(i,j+2,k)-Frhow(i,j-2,k) ) &
                            + a7(3) * ( Frhow(i,j+3,k)-Frhow(i,j-3,k) ) )*idy6_jmin &
                            + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a7(1) * ( Frhoe(i,j+1,k)-Frhoe(i,j-1,k) ) &
                            + a7(2) * ( Frhoe(i,j+2,k)-Frhoe(i,j-2,k) ) &
                            + a7(3) * ( Frhoe(i,j+3,k)-Frhoe(i,j-3,k) ) )*idy6_jmin &
                            + Krhoe(i,j,k)
          enddo
       enddo

       j=5
       do k=ndz_e,nfz_e
          do i=ndx_e,nfx_e
             Krho(i,j,k)  = ( a9(1) * ( Frho(i,j+1,k)-Frho(i,j-1,k) ) &
                            + a9(2) * ( Frho(i,j+2,k)-Frho(i,j-2,k) ) &
                            + a9(3) * ( Frho(i,j+3,k)-Frho(i,j-3,k) ) &
                            + a9(4) * ( Frho(i,j+4,k)-Frho(i,j-4,k) ) )*idy8_jmin &
                            + Krho(i,j,k)
             Krhou(i,j,k) = ( a9(1) * ( Frhou(i,j+1,k)-Frhou(i,j-1,k) ) &
                            + a9(2) * ( Frhou(i,j+2,k)-Frhou(i,j-2,k) ) &
                            + a9(3) * ( Frhou(i,j+3,k)-Frhou(i,j-3,k) ) &
                            + a9(4) * ( Frhou(i,j+4,k)-Frhou(i,j-4,k) ) )*idy8_jmin &
                            + Krhou(i,j,k)
             Krhov(i,j,k) = ( a9(1) * ( Frhov(i,j+1,k)-Frhov(i,j-1,k) ) &
                            + a9(2) * ( Frhov(i,j+2,k)-Frhov(i,j-2,k) ) &
                            + a9(3) * ( Frhov(i,j+3,k)-Frhov(i,j-3,k) ) &
                            + a9(4) * ( Frhov(i,j+4,k)-Frhov(i,j-4,k) ) )*idy8_jmin &
                            + Krhov(i,j,k)
             Krhow(i,j,k) = ( a9(1) * ( Frhow(i,j+1,k)-Frhow(i,j-1,k) ) &
                            + a9(2) * ( Frhow(i,j+2,k)-Frhow(i,j-2,k) ) &
                            + a9(3) * ( Frhow(i,j+3,k)-Frhow(i,j-3,k) ) &
                            + a9(4) * ( Frhow(i,j+4,k)-Frhow(i,j-4,k) ) )*idy8_jmin &
                            + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a9(1) * ( Frhoe(i,j+1,k)-Frhoe(i,j-1,k) ) &
                            + a9(2) * ( Frhoe(i,j+2,k)-Frhoe(i,j-2,k) ) &
                            + a9(3) * ( Frhoe(i,j+3,k)-Frhoe(i,j-3,k) ) &
                            + a9(4) * ( Frhoe(i,j+4,k)-Frhoe(i,j-4,k) ) )*idy8_jmin &
                            + Krhoe(i,j,k)
          enddo
       enddo
    endif

    ! Interior points
    ! ---------------
    do k=ndz_e,nfz_e
       do j=ndy,nfy
          do i=ndx_e,nfx_e
             Krho(i,j,k)  = ( a11(1) * ( Frho(i,j+1,k)-Frho(i,j-1,k) ) &
                            + a11(2) * ( Frho(i,j+2,k)-Frho(i,j-2,k) ) &
                            + a11(3) * ( Frho(i,j+3,k)-Frho(i,j-3,k) ) &
                            + a11(4) * ( Frho(i,j+4,k)-Frho(i,j-4,k) ) &
                            + a11(5) * ( Frho(i,j+5,k)-Frho(i,j-5,k) ) )*idy(j) &
                            + Krho(i,j,k)

             Krhou(i,j,k) = ( a11(1) * ( Frhou(i,j+1,k)-Frhou(i,j-1,k) ) &
                            + a11(2) * ( Frhou(i,j+2,k)-Frhou(i,j-2,k) ) &
                            + a11(3) * ( Frhou(i,j+3,k)-Frhou(i,j-3,k) ) &
                            + a11(4) * ( Frhou(i,j+4,k)-Frhou(i,j-4,k) ) &
                            + a11(5) * ( Frhou(i,j+5,k)-Frhou(i,j-5,k) ) )*idy(j) &
                            + Krhou(i,j,k)

             Krhov(i,j,k) = ( a11(1) * ( Frhov(i,j+1,k)-Frhov(i,j-1,k) ) &
                            + a11(2) * ( Frhov(i,j+2,k)-Frhov(i,j-2,k) ) &
                            + a11(3) * ( Frhov(i,j+3,k)-Frhov(i,j-3,k) ) &
                            + a11(4) * ( Frhov(i,j+4,k)-Frhov(i,j-4,k) ) &
                            + a11(5) * ( Frhov(i,j+5,k)-Frhov(i,j-5,k) ) )*idy(j) &
                            + Krhov(i,j,k)

             Krhow(i,j,k) = ( a11(1) * ( Frhow(i,j+1,k)-Frhow(i,j-1,k) ) &
                            + a11(2) * ( Frhow(i,j+2,k)-Frhow(i,j-2,k) ) &
                            + a11(3) * ( Frhow(i,j+3,k)-Frhow(i,j-3,k) ) &
                            + a11(4) * ( Frhow(i,j+4,k)-Frhow(i,j-4,k) ) &
                            + a11(5) * ( Frhow(i,j+5,k)-Frhow(i,j-5,k) ) )*idy(j) &
                            + Krhow(i,j,k)

             Krhoe(i,j,k) = ( a11(1) * ( Frhoe(i,j+1,k)-Frhoe(i,j-1,k) ) &
                            + a11(2) * ( Frhoe(i,j+2,k)-Frhoe(i,j-2,k) ) &
                            + a11(3) * ( Frhoe(i,j+3,k)-Frhoe(i,j-3,k) ) &
                            + a11(4) * ( Frhoe(i,j+4,k)-Frhoe(i,j-4,k) ) &
                            + a11(5) * ( Frhoe(i,j+5,k)-Frhoe(i,j-5,k) ) )*idy(j) &
                            + Krhoe(i,j,k)
          enddo
       enddo
    enddo

    ! BC at jmax
    ! ----------
    if (is_bc_1pt(2,2)) then
       j=ny-4
       do k=ndz_e,nfz_e
          do i=ndx_e,nfx_e
             Krho(i,j,k)  = ( a9(1) * ( Frho(i,j+1,k)-Frho(i,j-1,k) ) &
                            + a9(2) * ( Frho(i,j+2,k)-Frho(i,j-2,k) ) &
                            + a9(3) * ( Frho(i,j+3,k)-Frho(i,j-3,k) ) &
                            + a9(4) * ( Frho(i,j+4,k)-Frho(i,j-4,k) ) )*idy8_jmax &
                            + Krho(i,j,k)
             Krhou(i,j,k) = ( a9(1) * ( Frhou(i,j+1,k)-Frhou(i,j-1,k) ) &
                            + a9(2) * ( Frhou(i,j+2,k)-Frhou(i,j-2,k) ) &
                            + a9(3) * ( Frhou(i,j+3,k)-Frhou(i,j-3,k) ) &
                            + a9(4) * ( Frhou(i,j+4,k)-Frhou(i,j-4,k) ) )*idy8_jmax &
                            + Krhou(i,j,k)
             Krhov(i,j,k) = ( a9(1) * ( Frhov(i,j+1,k)-Frhov(i,j-1,k) ) &
                            + a9(2) * ( Frhov(i,j+2,k)-Frhov(i,j-2,k) ) &
                            + a9(3) * ( Frhov(i,j+3,k)-Frhov(i,j-3,k) ) &
                            + a9(4) * ( Frhov(i,j+4,k)-Frhov(i,j-4,k) ) )*idy8_jmax &
                            + Krhov(i,j,k)
             Krhow(i,j,k) = ( a9(1) * ( Frhow(i,j+1,k)-Frhow(i,j-1,k) ) &
                            + a9(2) * ( Frhow(i,j+2,k)-Frhow(i,j-2,k) ) &
                            + a9(3) * ( Frhow(i,j+3,k)-Frhow(i,j-3,k) ) &
                            + a9(4) * ( Frhow(i,j+4,k)-Frhow(i,j-4,k) ) )*idy8_jmax &
                            + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a9(1) * ( Frhoe(i,j+1,k)-Frhoe(i,j-1,k) ) &
                            + a9(2) * ( Frhoe(i,j+2,k)-Frhoe(i,j-2,k) ) &
                            + a9(3) * ( Frhoe(i,j+3,k)-Frhoe(i,j-3,k) ) &
                            + a9(4) * ( Frhoe(i,j+4,k)-Frhoe(i,j-4,k) ) )*idy8_jmax &
                            + Krhoe(i,j,k)
          enddo
       enddo

       j=ny-3
       do k=ndz_e,nfz_e
          do i=ndx_e,nfx_e
             Krho(i,j,k)  = ( a7(1) * ( Frho(i,j+1,k)-Frho(i,j-1,k) ) &
                            + a7(2) * ( Frho(i,j+2,k)-Frho(i,j-2,k) ) &
                            + a7(3) * ( Frho(i,j+3,k)-Frho(i,j-3,k) ) )*idy6_jmax &
                            + Krho(i,j,k)
             Krhou(i,j,k) = ( a7(1) * ( Frhou(i,j+1,k)-Frhou(i,j-1,k) ) &
                            + a7(2) * ( Frhou(i,j+2,k)-Frhou(i,j-2,k) ) &
                            + a7(3) * ( Frhou(i,j+3,k)-Frhou(i,j-3,k) ) )*idy6_jmax &
                            + Krhou(i,j,k)
             Krhov(i,j,k) = ( a7(1) * ( Frhov(i,j+1,k)-Frhov(i,j-1,k) ) &
                            + a7(2) * ( Frhov(i,j+2,k)-Frhov(i,j-2,k) ) &
                            + a7(3) * ( Frhov(i,j+3,k)-Frhov(i,j-3,k) ) )*idy6_jmax &
                            + Krhov(i,j,k)
             Krhow(i,j,k) = ( a7(1) * ( Frhow(i,j+1,k)-Frhow(i,j-1,k) ) &
                            + a7(2) * ( Frhow(i,j+2,k)-Frhow(i,j-2,k) ) &
                            + a7(3) * ( Frhow(i,j+3,k)-Frhow(i,j-3,k) ) )*idy6_jmax &
                            + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a7(1) * ( Frhoe(i,j+1,k)-Frhoe(i,j-1,k) ) &
                            + a7(2) * ( Frhoe(i,j+2,k)-Frhoe(i,j-2,k) ) &
                            + a7(3) * ( Frhoe(i,j+3,k)-Frhoe(i,j-3,k) ) )*idy6_jmax &
                            + Krhoe(i,j,k)
          enddo
       enddo

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

       k=4
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Krho(i,j,k)  = ( a7(1) * ( Frho(i,j,k+1)-Frho(i,j,k-1) ) &
                            + a7(2) * ( Frho(i,j,k+2)-Frho(i,j,k-2) ) &
                            + a7(3) * ( Frho(i,j,k+3)-Frho(i,j,k-3) ) )*idz6_kmin &
                            + Krho(i,j,k)
             Krhou(i,j,k) = ( a7(1) * ( Frhou(i,j,k+1)-Frhou(i,j,k-1) ) &
                            + a7(2) * ( Frhou(i,j,k+2)-Frhou(i,j,k-2) ) &
                            + a7(3) * ( Frhou(i,j,k+3)-Frhou(i,j,k-3) ) )*idz6_kmin &
                            + Krhou(i,j,k)
             Krhov(i,j,k) = ( a7(1) * ( Frhov(i,j,k+1)-Frhov(i,j,k-1) ) &
                            + a7(2) * ( Frhov(i,j,k+2)-Frhov(i,j,k-2) ) &
                            + a7(3) * ( Frhov(i,j,k+3)-Frhov(i,j,k-3) ) )*idz6_kmin &
                            + Krhov(i,j,k)
             Krhow(i,j,k) = ( a7(1) * ( Frhow(i,j,k+1)-Frhow(i,j,k-1) ) &
                            + a7(2) * ( Frhow(i,j,k+2)-Frhow(i,j,k-2) ) &
                            + a7(3) * ( Frhow(i,j,k+3)-Frhow(i,j,k-3) ) )*idz6_kmin &
                            + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a7(1) * ( Frhoe(i,j,k+1)-Frhoe(i,j,k-1) ) &
                            + a7(2) * ( Frhoe(i,j,k+2)-Frhoe(i,j,k-2) ) &
                            + a7(3) * ( Frhoe(i,j,k+3)-Frhoe(i,j,k-3) ) )*idz6_kmin &
                            + Krhoe(i,j,k)
          enddo
       enddo

       k=5
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Krho(i,j,k)  = ( a9(1) * ( Frho(i,j,k+1)-Frho(i,j,k-1) ) &
                            + a9(2) * ( Frho(i,j,k+2)-Frho(i,j,k-2) ) &
                            + a9(3) * ( Frho(i,j,k+3)-Frho(i,j,k-3) ) &
                            + a9(4) * ( Frho(i,j,k+4)-Frho(i,j,k-4) ) )*idz8_kmin &
                            + Krho(i,j,k)
             Krhou(i,j,k) = ( a9(1) * ( Frhou(i,j,k+1)-Frhou(i,j,k-1) ) &
                            + a9(2) * ( Frhou(i,j,k+2)-Frhou(i,j,k-2) ) &
                            + a9(3) * ( Frhou(i,j,k+3)-Frhou(i,j,k-3) ) &
                            + a9(4) * ( Frhou(i,j,k+4)-Frhou(i,j,k-4) ) )*idz8_kmin &
                            + Krhou(i,j,k)
             Krhov(i,j,k) = ( a9(1) * ( Frhov(i,j,k+1)-Frhov(i,j,k-1) ) &
                            + a9(2) * ( Frhov(i,j,k+2)-Frhov(i,j,k-2) ) &
                            + a9(3) * ( Frhov(i,j,k+3)-Frhov(i,j,k-3) ) &
                            + a9(4) * ( Frhov(i,j,k+4)-Frhov(i,j,k-4) ) )*idz8_kmin &
                            + Krhov(i,j,k)
             Krhow(i,j,k) = ( a9(1) * ( Frhow(i,j,k+1)-Frhow(i,j,k-1) ) &
                            + a9(2) * ( Frhow(i,j,k+2)-Frhow(i,j,k-2) ) &
                            + a9(3) * ( Frhow(i,j,k+3)-Frhow(i,j,k-3) ) &
                            + a9(4) * ( Frhow(i,j,k+4)-Frhow(i,j,k-4) ) )*idz8_kmin &
                            + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a9(1) * ( Frhoe(i,j,k+1)-Frhoe(i,j,k-1) ) &
                            + a9(2) * ( Frhoe(i,j,k+2)-Frhoe(i,j,k-2) ) &
                            + a9(3) * ( Frhoe(i,j,k+3)-Frhoe(i,j,k-3) ) &
                            + a9(4) * ( Frhoe(i,j,k+4)-Frhoe(i,j,k-4) ) )*idz8_kmin &
                            + Krhoe(i,j,k)
          enddo
       enddo
    endif

    ! Interior points
    ! ---------------
    do k=ndz,nfz
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Krho(i,j,k)  = ( a11(1)* ( Frho(i,j,k+1)-Frho(i,j,k-1) ) &
                            + a11(2)* ( Frho(i,j,k+2)-Frho(i,j,k-2) ) &
                            + a11(3)* ( Frho(i,j,k+3)-Frho(i,j,k-3) ) &
                            + a11(4)* ( Frho(i,j,k+4)-Frho(i,j,k-4) ) &
                            + a11(5)* ( Frho(i,j,k+5)-Frho(i,j,k-5) ) )*idz(k) &
                            + Krho(i,j,k)

             Krhou(i,j,k) = ( a11(1)* ( Frhou(i,j,k+1)-Frhou(i,j,k-1) ) &
                            + a11(2)* ( Frhou(i,j,k+2)-Frhou(i,j,k-2) ) &
                            + a11(3)* ( Frhou(i,j,k+3)-Frhou(i,j,k-3) ) &
                            + a11(4)* ( Frhou(i,j,k+4)-Frhou(i,j,k-4) ) &
                            + a11(5)* ( Frhou(i,j,k+5)-Frhou(i,j,k-5) ) )*idz(k) &
                            + Krhou(i,j,k)

             Krhov(i,j,k) = ( a11(1)* ( Frhov(i,j,k+1)-Frhov(i,j,k-1) ) &
                            + a11(2)* ( Frhov(i,j,k+2)-Frhov(i,j,k-2) ) &
                            + a11(3)* ( Frhov(i,j,k+3)-Frhov(i,j,k-3) ) &
                            + a11(4)* ( Frhov(i,j,k+4)-Frhov(i,j,k-4) ) &
                            + a11(5)* ( Frhov(i,j,k+5)-Frhov(i,j,k-5) ) )*idz(k) &
                            + Krhov(i,j,k)

             Krhow(i,j,k) = ( a11(1)* ( Frhow(i,j,k+1)-Frhow(i,j,k-1) ) &
                            + a11(2)* ( Frhow(i,j,k+2)-Frhow(i,j,k-2) ) &
                            + a11(3)* ( Frhow(i,j,k+3)-Frhow(i,j,k-3) ) &
                            + a11(4)* ( Frhow(i,j,k+4)-Frhow(i,j,k-4) ) &
                            + a11(5)* ( Frhow(i,j,k+5)-Frhow(i,j,k-5) ) )*idz(k) &
                            + Krhow(i,j,k)

             Krhoe(i,j,k) = ( a11(1)* ( Frhoe(i,j,k+1)-Frhoe(i,j,k-1) ) &
                            + a11(2)* ( Frhoe(i,j,k+2)-Frhoe(i,j,k-2) ) &
                            + a11(3)* ( Frhoe(i,j,k+3)-Frhoe(i,j,k-3) ) &
                            + a11(4)* ( Frhoe(i,j,k+4)-Frhoe(i,j,k-4) ) &
                            + a11(5)* ( Frhoe(i,j,k+5)-Frhoe(i,j,k-5) ) )*idz(k) &
                            + Krhoe(i,j,k)
          enddo
       enddo
    enddo

    ! BC at kmax
    ! ----------
    if (is_bc_1pt(3,2)) then
       k=nz-4
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Krho(i,j,k)  = ( a9(1) * ( Frho(i,j,k+1)-Frho(i,j,k-1) ) &
                            + a9(2) * ( Frho(i,j,k+2)-Frho(i,j,k-2) ) &
                            + a9(3) * ( Frho(i,j,k+3)-Frho(i,j,k-3) ) &
                            + a9(4) * ( Frho(i,j,k+4)-Frho(i,j,k-4) ) )*idz8_kmax &
                            + Krho(i,j,k)
             Krhou(i,j,k) = ( a9(1) * ( Frhou(i,j,k+1)-Frhou(i,j,k-1) ) &
                            + a9(2) * ( Frhou(i,j,k+2)-Frhou(i,j,k-2) ) &
                            + a9(3) * ( Frhou(i,j,k+3)-Frhou(i,j,k-3) ) &
                            + a9(4) * ( Frhou(i,j,k+4)-Frhou(i,j,k-4) ) )*idz8_kmax &
                            + Krhou(i,j,k)
             Krhov(i,j,k) = ( a9(1) * ( Frhov(i,j,k+1)-Frhov(i,j,k-1) ) &
                            + a9(2) * ( Frhov(i,j,k+2)-Frhov(i,j,k-2) ) &
                            + a9(3) * ( Frhov(i,j,k+3)-Frhov(i,j,k-3) ) &
                            + a9(4) * ( Frhov(i,j,k+4)-Frhov(i,j,k-4) ) )*idz8_kmax &
                            + Krhov(i,j,k)
             Krhow(i,j,k) = ( a9(1) * ( Frhow(i,j,k+1)-Frhow(i,j,k-1) ) &
                            + a9(2) * ( Frhow(i,j,k+2)-Frhow(i,j,k-2) ) &
                            + a9(3) * ( Frhow(i,j,k+3)-Frhow(i,j,k-3) ) &
                            + a9(4) * ( Frhow(i,j,k+4)-Frhow(i,j,k-4) ) )*idz8_kmax &
                            + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a9(1) * ( Frhoe(i,j,k+1)-Frhoe(i,j,k-1) ) &
                            + a9(2) * ( Frhoe(i,j,k+2)-Frhoe(i,j,k-2) ) &
                            + a9(3) * ( Frhoe(i,j,k+3)-Frhoe(i,j,k-3) ) &
                            + a9(4) * ( Frhoe(i,j,k+4)-Frhoe(i,j,k-4) ) )*idz8_kmax &
                            + Krhoe(i,j,k)
          enddo
       enddo

       k=nz-3
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Krho(i,j,k)  = ( a7(1) * ( Frho(i,j,k+1)-Frho(i,j,k-1) ) &
                            + a7(2) * ( Frho(i,j,k+2)-Frho(i,j,k-2) ) &
                            + a7(3) * ( Frho(i,j,k+3)-Frho(i,j,k-3) ) )*idz6_kmax &
                            + Krho(i,j,k)
             Krhou(i,j,k) = ( a7(1) * ( Frhou(i,j,k+1)-Frhou(i,j,k-1) ) &
                            + a7(2) * ( Frhou(i,j,k+2)-Frhou(i,j,k-2) ) &
                            + a7(3) * ( Frhou(i,j,k+3)-Frhou(i,j,k-3) ) )*idz6_kmax &
                            + Krhou(i,j,k)
             Krhov(i,j,k) = ( a7(1) * ( Frhov(i,j,k+1)-Frhov(i,j,k-1) ) &
                            + a7(2) * ( Frhov(i,j,k+2)-Frhov(i,j,k-2) ) &
                            + a7(3) * ( Frhov(i,j,k+3)-Frhov(i,j,k-3) ) )*idz6_kmax &
                            + Krhov(i,j,k)
             Krhow(i,j,k) = ( a7(1) * ( Frhow(i,j,k+1)-Frhow(i,j,k-1) ) &
                            + a7(2) * ( Frhow(i,j,k+2)-Frhow(i,j,k-2) ) &
                            + a7(3) * ( Frhow(i,j,k+3)-Frhow(i,j,k-3) ) )*idz6_kmax &
                            + Krhow(i,j,k)
             Krhoe(i,j,k) = ( a7(1) * ( Frhoe(i,j,k+1)-Frhoe(i,j,k-1) ) &
                            + a7(2) * ( Frhoe(i,j,k+2)-Frhoe(i,j,k-2) ) &
                            + a7(3) * ( Frhoe(i,j,k+3)-Frhoe(i,j,k-3) ) )*idz6_kmax &
                            + Krhoe(i,j,k)
          enddo
       enddo

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

  end subroutine flux_euler_11pts

end submodule smod_flux_euler_11pts
