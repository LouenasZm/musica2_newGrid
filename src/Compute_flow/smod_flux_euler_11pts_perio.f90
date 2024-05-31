!=================================================================================
submodule (mod_flux_euler) smod_flux_euler_11pts_perio
!=================================================================================
  !> author: XG
  !> date: January 2024
  !> Cartesian version - compute Eulerian fluxes (inviscid part) with 11-pt stencil
  !> Simplified periodic version ** FOR TESTS ONLY **
!=================================================================================

contains

  !==============================================================================
  subroutine flux_euler_11pts_p
  !==============================================================================
    !> Derivatives of Eulerian fluxes (inviscid part) - 11-point stencil -
    !> - Cartesian version - Simplified periodic version
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
    do k=ndz_e,nfz_e
       do j=ndy_e,nfy_e
          do i=ndx,nfx

             Krho(i,j,k)  = ( a11(1) * ( Frho(i+1,j,k)-Frho(i-1,j,k) )   &
                            + a11(2) * ( Frho(i+2,j,k)-Frho(i-2,j,k) )   &
                            + a11(3) * ( Frho(i+3,j,k)-Frho(i-3,j,k) )   &
                            + a11(4) * ( Frho(i+4,j,k)-Frho(i-4,j,k) )   &
                            + a11(5) * ( Frho(i+5,j,k)-Frho(i-5,j,k) ) ) *idx(i)

             Krhou(i,j,k) = ( a11(1) * ( Frhou(i+1,j,k)-Frhou(i-1,j,k) ) &
                            + a11(2) * ( Frhou(i+2,j,k)-Frhou(i-2,j,k) ) &
                            + a11(3) * ( Frhou(i+3,j,k)-Frhou(i-3,j,k) ) &
                            + a11(4) * ( Frhou(i+4,j,k)-Frhou(i-4,j,k) ) &
                            + a11(5) * ( Frhou(i+5,j,k)-Frhou(i-5,j,k) ) ) *idx(i)

             Krhov(i,j,k) = ( a11(1) * ( Frhov(i+1,j,k)-Frhov(i-1,j,k) ) &
                            + a11(2) * ( Frhov(i+2,j,k)-Frhov(i-2,j,k) ) &
                            + a11(3) * ( Frhov(i+3,j,k)-Frhov(i-3,j,k) ) &
                            + a11(4) * ( Frhov(i+4,j,k)-Frhov(i-4,j,k) ) &
                            + a11(5) * ( Frhov(i+5,j,k)-Frhov(i-5,j,k) ) ) *idx(i)

             Krhow(i,j,k) = ( a11(1) * ( Frhow(i+1,j,k)-Frhow(i-1,j,k) ) &
                            + a11(2) * ( Frhow(i+2,j,k)-Frhow(i-2,j,k) ) &
                            + a11(3) * ( Frhow(i+3,j,k)-Frhow(i-3,j,k) ) &
                            + a11(4) * ( Frhow(i+4,j,k)-Frhow(i-4,j,k) ) &
                            + a11(5) * ( Frhow(i+5,j,k)-Frhow(i-5,j,k) ) ) *idx(i)

             Krhoe(i,j,k) = ( a11(1) * ( Frhoe(i+1,j,k)-Frhoe(i-1,j,k) ) &
                            + a11(2) * ( Frhoe(i+2,j,k)-Frhoe(i-2,j,k) ) &
                            + a11(3) * ( Frhoe(i+3,j,k)-Frhoe(i-3,j,k) ) &
                            + a11(4) * ( Frhoe(i+4,j,k)-Frhoe(i-4,j,k) ) &
                            + a11(5) * ( Frhoe(i+5,j,k)-Frhoe(i-5,j,k) ) ) *idx(i)
          enddo
       enddo
    enddo

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
    do k=ndz_e,nfz_e
       do j=ndy,nfy
          do i=ndx_e,nfx_e
             Krho(i,j,k)  = ( a11(1) * ( Frho(i,j+1,k)-Frho(i,j-1,k) ) &
                            + a11(2) * ( Frho(i,j+2,k)-Frho(i,j-2,k) ) &
                            + a11(3) * ( Frho(i,j+3,k)-Frho(i,j-3,k) ) &
                            + a11(4) * ( Frho(i,j+4,k)-Frho(i,j-4,k) ) &
                            + a11(5) * ( Frho(i,j+5,k)-Frho(i,j-5,k) ) ) *idy(j) &
                            + Krho(i,j,k)

             Krhou(i,j,k) = ( a11(1) * ( Frhou(i,j+1,k)-Frhou(i,j-1,k) ) &
                            + a11(2) * ( Frhou(i,j+2,k)-Frhou(i,j-2,k) ) &
                            + a11(3) * ( Frhou(i,j+3,k)-Frhou(i,j-3,k) ) &
                            + a11(4) * ( Frhou(i,j+4,k)-Frhou(i,j-4,k) ) &
                            + a11(5) * ( Frhou(i,j+5,k)-Frhou(i,j-5,k) ) ) *idy(j) &
                            + Krhou(i,j,k)

             Krhov(i,j,k) = ( a11(1) * ( Frhov(i,j+1,k)-Frhov(i,j-1,k) ) &
                            + a11(2) * ( Frhov(i,j+2,k)-Frhov(i,j-2,k) ) &
                            + a11(3) * ( Frhov(i,j+3,k)-Frhov(i,j-3,k) ) &
                            + a11(4) * ( Frhov(i,j+4,k)-Frhov(i,j-4,k) ) &
                            + a11(5) * ( Frhov(i,j+5,k)-Frhov(i,j-5,k) ) ) *idy(j) &
                            + Krhov(i,j,k)

             Krhow(i,j,k) = ( a11(1) * ( Frhow(i,j+1,k)-Frhow(i,j-1,k) ) &
                            + a11(2) * ( Frhow(i,j+2,k)-Frhow(i,j-2,k) ) &
                            + a11(3) * ( Frhow(i,j+3,k)-Frhow(i,j-3,k) ) &
                            + a11(4) * ( Frhow(i,j+4,k)-Frhow(i,j-4,k) ) &
                            + a11(5) * ( Frhow(i,j+5,k)-Frhow(i,j-5,k) ) ) *idy(j) &
                            + Krhow(i,j,k)

             Krhoe(i,j,k) = ( a11(1) * ( Frhoe(i,j+1,k)-Frhoe(i,j-1,k) ) &
                            + a11(2) * ( Frhoe(i,j+2,k)-Frhoe(i,j-2,k) ) &
                            + a11(3) * ( Frhoe(i,j+3,k)-Frhoe(i,j-3,k) ) &
                            + a11(4) * ( Frhoe(i,j+4,k)-Frhoe(i,j-4,k) ) &
                            + a11(5) * ( Frhoe(i,j+5,k)-Frhoe(i,j-5,k) ) ) *idy(j) &
                            + Krhoe(i,j,k)
          enddo
       enddo
    enddo


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
    do k=ndz,nfz
       do j=ndy_e,nfy_e
          do i=ndx_e,nfx_e
             Krho(i,j,k)  = ( a11(1)* ( Frho(i,j,k+1)-Frho(i,j,k-1) ) &
                            + a11(2)* ( Frho(i,j,k+2)-Frho(i,j,k-2) ) &
                            + a11(3)* ( Frho(i,j,k+3)-Frho(i,j,k-3) ) &
                            + a11(4)* ( Frho(i,j,k+4)-Frho(i,j,k-4) ) &
                            + a11(5)* ( Frho(i,j,k+5)-Frho(i,j,k-5) ) ) *idz(k) &
                            + Krho(i,j,k)

             Krhou(i,j,k) = ( a11(1)* ( Frhou(i,j,k+1)-Frhou(i,j,k-1) ) &
                            + a11(2)* ( Frhou(i,j,k+2)-Frhou(i,j,k-2) ) &
                            + a11(3)* ( Frhou(i,j,k+3)-Frhou(i,j,k-3) ) &
                            + a11(4)* ( Frhou(i,j,k+4)-Frhou(i,j,k-4) ) &
                            + a11(5)* ( Frhou(i,j,k+5)-Frhou(i,j,k-5) ) ) *idz(k) &
                            + Krhou(i,j,k)

             Krhov(i,j,k) = ( a11(1)* ( Frhov(i,j,k+1)-Frhov(i,j,k-1) ) &
                            + a11(2)* ( Frhov(i,j,k+2)-Frhov(i,j,k-2) ) &
                            + a11(3)* ( Frhov(i,j,k+3)-Frhov(i,j,k-3) ) &
                            + a11(4)* ( Frhov(i,j,k+4)-Frhov(i,j,k-4) ) &
                            + a11(5)* ( Frhov(i,j,k+5)-Frhov(i,j,k-5) ) ) *idz(k) &
                            + Krhov(i,j,k)

             Krhow(i,j,k) = ( a11(1)* ( Frhow(i,j,k+1)-Frhow(i,j,k-1) ) &
                            + a11(2)* ( Frhow(i,j,k+2)-Frhow(i,j,k-2) ) &
                            + a11(3)* ( Frhow(i,j,k+3)-Frhow(i,j,k-3) ) &
                            + a11(4)* ( Frhow(i,j,k+4)-Frhow(i,j,k-4) ) &
                            + a11(5)* ( Frhow(i,j,k+5)-Frhow(i,j,k-5) ) ) *idz(k) &
                            + Krhow(i,j,k)

             Krhoe(i,j,k) = ( a11(1)* ( Frhoe(i,j,k+1)-Frhoe(i,j,k-1) ) &
                            + a11(2)* ( Frhoe(i,j,k+2)-Frhoe(i,j,k-2) ) &
                            + a11(3)* ( Frhoe(i,j,k+3)-Frhoe(i,j,k-3) ) &
                            + a11(4)* ( Frhoe(i,j,k+4)-Frhoe(i,j,k-4) ) &
                            + a11(5)* ( Frhoe(i,j,k+5)-Frhoe(i,j,k-5) ) ) *idz(k) &
                            + Krhoe(i,j,k)
          enddo
       enddo
    enddo

  end subroutine flux_euler_11pts_p

end submodule smod_flux_euler_11pts_perio
