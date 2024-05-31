  !===============================================================================
  subroutine filtering_11pts_overlap
  !===============================================================================
    !> Apply filtering on 11-point stencils (+ boundaries)
    !> with computation-communication overlap
  !===============================================================================
    use mod_flow
    use mod_filtering
    use mod_comm1
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    ! ---------------------------------------------------------------------------

    ! Filtering along i-direction
    ! ===========================

    ! Perform one-sided communications in i-direction
    ! -----------------------------------------------
    ! open target windows
    call window_comm1
    ! communicate var
    call communication_1_i
    
    ! Interior points
    ! ---------------
    ! 10-th order centered filter
    do k=1,nz
       do j=1,ny
          do i=6,nx-5
             Krho(i,j,k) = d11(0)*rho_n(i,j,k) &
                  + d11(1)*( rho_n(i+1,j,k)+rho_n(i-1,j,k) )&
                  + d11(2)*( rho_n(i+2,j,k)+rho_n(i-2,j,k) )&
                  + d11(3)*( rho_n(i+3,j,k)+rho_n(i-3,j,k) )&
                  + d11(4)*( rho_n(i+4,j,k)+rho_n(i-4,j,k) )&
                  + d11(5)*( rho_n(i+5,j,k)+rho_n(i-5,j,k) )

             Krhou(i,j,k)= d11(0)*rhou_n(i,j,k) &
                  + d11(1)*( rhou_n(i+1,j,k)+rhou_n(i-1,j,k) )&
                  + d11(2)*( rhou_n(i+2,j,k)+rhou_n(i-2,j,k) )&
                  + d11(3)*( rhou_n(i+3,j,k)+rhou_n(i-3,j,k) )&
                  + d11(4)*( rhou_n(i+4,j,k)+rhou_n(i-4,j,k) )&
                  + d11(5)*( rhou_n(i+5,j,k)+rhou_n(i-5,j,k) )

             Krhov(i,j,k)= d11(0)*rhov_n(i,j,k) &
                  + d11(1)*( rhov_n(i+1,j,k)+rhov_n(i-1,j,k) )&
                  + d11(2)*( rhov_n(i+2,j,k)+rhov_n(i-2,j,k) )&
                  + d11(3)*( rhov_n(i+3,j,k)+rhov_n(i-3,j,k) )&
                  + d11(4)*( rhov_n(i+4,j,k)+rhov_n(i-4,j,k) )&
                  + d11(5)*( rhov_n(i+5,j,k)+rhov_n(i-5,j,k) )

             Krhow(i,j,k)= d11(0)*rhow_n(i,j,k) &
                  + d11(1)*( rhow_n(i+1,j,k)+rhow_n(i-1,j,k) )&
                  + d11(2)*( rhow_n(i+2,j,k)+rhow_n(i-2,j,k) )&
                  + d11(3)*( rhow_n(i+3,j,k)+rhow_n(i-3,j,k) )&
                  + d11(4)*( rhow_n(i+4,j,k)+rhow_n(i-4,j,k) )&
                  + d11(5)*( rhow_n(i+5,j,k)+rhow_n(i-5,j,k) )

             Krhoe(i,j,k)= d11(0)*rhoe_n(i,j,k) &
                  + d11(1)*( rhoe_n(i+1,j,k)+rhoe_n(i-1,j,k) )&
                  + d11(2)*( rhoe_n(i+2,j,k)+rhoe_n(i-2,j,k) )&
                  + d11(3)*( rhoe_n(i+3,j,k)+rhoe_n(i-3,j,k) )&
                  + d11(4)*( rhoe_n(i+4,j,k)+rhoe_n(i-4,j,k) )&
                  + d11(5)*( rhoe_n(i+5,j,k)+rhoe_n(i-5,j,k) )
          enddo
       enddo
    enddo
    
    ! Check that communications are fullfilled
    ! ----------------------------------------
    call window_comm1
    
    ! BC at imin
    ! ----------
    if (is_boundary(1,1)) then
       ! no filter for i=1 & i=2 -> initialize array
       do j=1,ny
           Krho(1,j,:)=0.0_wp
          Krhou(1,j,:)=0.0_wp
          Krhov(1,j,:)=0.0_wp
          Krhow(1,j,:)=0.0_wp
          Krhoe(1,j,:)=0.0_wp
          
           Krho(2,j,:)=0.0_wp
          Krhou(2,j,:)=0.0_wp
          Krhov(2,j,:)=0.0_wp
          Krhow(2,j,:)=0.0_wp
          Krhoe(2,j,:)=0.0_wp
       enddo
       ! 4-th order centered filter
       i=3
       do k=1,nz
          do j=1,ny
             Krho(i,j,k) = d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i+1,j,k)+d5(2)*rho_n(i+2,j,k) &
                  + d5(1)*rho_n(i-1,j,k)+d5(2)*rho_n(i-2,j,k)
             Krhou(i,j,k)= d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i+1,j,k)+d5(2)*rhou_n(i+2,j,k) &
                  + d5(1)*rhou_n(i-1,j,k)+d5(2)*rhou_n(i-2,j,k)
             Krhov(i,j,k)= d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i+1,j,k)+d5(2)*rhov_n(i+2,j,k) &
                  + d5(1)*rhov_n(i-1,j,k)+d5(2)*rhov_n(i-2,j,k)
             Krhow(i,j,k)= d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i+1,j,k)+d5(2)*rhow_n(i+2,j,k) &
                  + d5(1)*rhow_n(i-1,j,k)+d5(2)*rhow_n(i-2,j,k)
             Krhoe(i,j,k)= d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i+1,j,k)+d5(2)*rhoe_n(i+2,j,k) &
                  + d5(1)*rhoe_n(i-1,j,k)+d5(2)*rhoe_n(i-2,j,k)
          enddo
       enddo
       ! 6-th order centered filter
       i=4
       do k=1,nz
          do j=1,ny
             Krho(i,j,k) = d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i+1,j,k)+d7(2)*rho_n(i+2,j,k)+d7(3)*rho_n(i+3,j,k) &
                  + d7(1)*rho_n(i-1,j,k)+d7(2)*rho_n(i-2,j,k)+d7(3)*rho_n(i-3,j,k)
             Krhou(i,j,k)= d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i+1,j,k)+d7(2)*rhou_n(i+2,j,k)+d7(3)*rhou_n(i+3,j,k) &
                  + d7(1)*rhou_n(i-1,j,k)+d7(2)*rhou_n(i-2,j,k)+d7(3)*rhou_n(i-3,j,k)
             Krhov(i,j,k)= d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i+1,j,k)+d7(2)*rhov_n(i+2,j,k)+d7(3)*rhov_n(i+3,j,k) &
                  + d7(1)*rhov_n(i-1,j,k)+d7(2)*rhov_n(i-2,j,k)+d7(3)*rhov_n(i-3,j,k)
             Krhow(i,j,k)= d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i+1,j,k)+d7(2)*rhow_n(i+2,j,k)+d7(3)*rhow_n(i+3,j,k) &
                  + d7(1)*rhow_n(i-1,j,k)+d7(2)*rhow_n(i-2,j,k)+d7(3)*rhow_n(i-3,j,k)
             Krhoe(i,j,k)= d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i+1,j,k)+d7(2)*rhoe_n(i+2,j,k)+d7(3)*rhoe_n(i+3,j,k) &
                  + d7(1)*rhoe_n(i-1,j,k)+d7(2)*rhoe_n(i-2,j,k)+d7(3)*rhoe_n(i-3,j,k)
          enddo
       enddo
       ! 8-th order centered filter
       i=5
       do k=1,nz
          do j=1,ny
             Krho(i,j,k) = d9(0)*rho_n(i,j,k)   &
                  + d9(1)*rho_n(i+1,j,k)+d9(2)*rho_n(i+2,j,k)+d9(3)*rho_n(i+3,j,k)+d9(4)*rho_n(i+4,j,k) &
                  + d9(1)*rho_n(i-1,j,k)+d9(2)*rho_n(i-2,j,k)+d9(3)*rho_n(i-3,j,k)+d9(4)*rho_n(i-4,j,k)
             Krhou(i,j,k)= d9(0)*rhou_n(i,j,k)   &
                  + d9(1)*rhou_n(i+1,j,k)+d9(2)*rhou_n(i+2,j,k)+d9(3)*rhou_n(i+3,j,k)+d9(4)*rhou_n(i+4,j,k) &
                  + d9(1)*rhou_n(i-1,j,k)+d9(2)*rhou_n(i-2,j,k)+d9(3)*rhou_n(i-3,j,k)+d9(4)*rhou_n(i-4,j,k)
             Krhov(i,j,k)= d9(0)*rhov_n(i,j,k)   &
                  + d9(1)*rhov_n(i+1,j,k)+d9(2)*rhov_n(i+2,j,k)+d9(3)*rhov_n(i+3,j,k)+d9(4)*rhov_n(i+4,j,k) &
                  + d9(1)*rhov_n(i-1,j,k)+d9(2)*rhov_n(i-2,j,k)+d9(3)*rhov_n(i-3,j,k)+d9(4)*rhov_n(i-4,j,k)
             Krhow(i,j,k)= d9(0)*rhow_n(i,j,k)   &
                  + d9(1)*rhow_n(i+1,j,k)+d9(2)*rhow_n(i+2,j,k)+d9(3)*rhow_n(i+3,j,k)+d9(4)*rhow_n(i+4,j,k) &
                  + d9(1)*rhow_n(i-1,j,k)+d9(2)*rhow_n(i-2,j,k)+d9(3)*rhow_n(i-3,j,k)+d9(4)*rhow_n(i-4,j,k)
             Krhoe(i,j,k)= d9(0)*rhoe_n(i,j,k)   &
                  + d9(1)*rhoe_n(i+1,j,k)+d9(2)*rhoe_n(i+2,j,k)+d9(3)*rhoe_n(i+3,j,k)+d9(4)*rhoe_n(i+4,j,k) &
                  + d9(1)*rhoe_n(i-1,j,k)+d9(2)*rhoe_n(i-2,j,k)+d9(3)*rhoe_n(i-3,j,k)+d9(4)*rhoe_n(i-4,j,k)
          enddo
       enddo
    else
       ! 10-th order centered filter
       do k=1,nz
          do j=1,ny
             do i=1,5
                Krho(i,j,k) = d11(0)*rho_n(i,j,k) &
                     + d11(1)*( rho_n(i+1,j,k)+rho_n(i-1,j,k) )&
                     + d11(2)*( rho_n(i+2,j,k)+rho_n(i-2,j,k) )&
                     + d11(3)*( rho_n(i+3,j,k)+rho_n(i-3,j,k) )&
                     + d11(4)*( rho_n(i+4,j,k)+rho_n(i-4,j,k) )&
                     + d11(5)*( rho_n(i+5,j,k)+rho_n(i-5,j,k) )

                Krhou(i,j,k)= d11(0)*rhou_n(i,j,k) &
                     + d11(1)*( rhou_n(i+1,j,k)+rhou_n(i-1,j,k) )&
                     + d11(2)*( rhou_n(i+2,j,k)+rhou_n(i-2,j,k) )&
                     + d11(3)*( rhou_n(i+3,j,k)+rhou_n(i-3,j,k) )&
                     + d11(4)*( rhou_n(i+4,j,k)+rhou_n(i-4,j,k) )&
                     + d11(5)*( rhou_n(i+5,j,k)+rhou_n(i-5,j,k) )

                Krhov(i,j,k)= d11(0)*rhov_n(i,j,k) &
                     + d11(1)*( rhov_n(i+1,j,k)+rhov_n(i-1,j,k) )&
                     + d11(2)*( rhov_n(i+2,j,k)+rhov_n(i-2,j,k) )&
                     + d11(3)*( rhov_n(i+3,j,k)+rhov_n(i-3,j,k) )&
                     + d11(4)*( rhov_n(i+4,j,k)+rhov_n(i-4,j,k) )&
                     + d11(5)*( rhov_n(i+5,j,k)+rhov_n(i-5,j,k) )

                Krhow(i,j,k)= d11(0)*rhow_n(i,j,k) &
                     + d11(1)*( rhow_n(i+1,j,k)+rhow_n(i-1,j,k) )&
                     + d11(2)*( rhow_n(i+2,j,k)+rhow_n(i-2,j,k) )&
                     + d11(3)*( rhow_n(i+3,j,k)+rhow_n(i-3,j,k) )&
                     + d11(4)*( rhow_n(i+4,j,k)+rhow_n(i-4,j,k) )&
                     + d11(5)*( rhow_n(i+5,j,k)+rhow_n(i-5,j,k) )

                Krhoe(i,j,k)= d11(0)*rhoe_n(i,j,k) &
                     + d11(1)*( rhoe_n(i+1,j,k)+rhoe_n(i-1,j,k) )&
                     + d11(2)*( rhoe_n(i+2,j,k)+rhoe_n(i-2,j,k) )&
                     + d11(3)*( rhoe_n(i+3,j,k)+rhoe_n(i-3,j,k) )&
                     + d11(4)*( rhoe_n(i+4,j,k)+rhoe_n(i-4,j,k) )&
                     + d11(5)*( rhoe_n(i+5,j,k)+rhoe_n(i-5,j,k) )
             enddo
          enddo
       enddo
    endif

    ! BC at imax
    ! ----------
    if (is_boundary(1,2)) then
       ! 8-th order centered filter
       i=nx-4
       do k=1,nz
          do j=1,ny
             Krho(i,j,k) = d9(0)*rho_n(i,j,k)   &
                  + d9(1)*rho_n(i+1,j,k)+d9(2)*rho_n(i+2,j,k)+d9(3)*rho_n(i+3,j,k)+d9(4)*rho_n(i+4,j,k) &
                  + d9(1)*rho_n(i-1,j,k)+d9(2)*rho_n(i-2,j,k)+d9(3)*rho_n(i-3,j,k)+d9(4)*rho_n(i-4,j,k)
             Krhou(i,j,k)= d9(0)*rhou_n(i,j,k)   &
                  + d9(1)*rhou_n(i+1,j,k)+d9(2)*rhou_n(i+2,j,k)+d9(3)*rhou_n(i+3,j,k)+d9(4)*rhou_n(i+4,j,k) &
                  + d9(1)*rhou_n(i-1,j,k)+d9(2)*rhou_n(i-2,j,k)+d9(3)*rhou_n(i-3,j,k)+d9(4)*rhou_n(i-4,j,k)
             Krhov(i,j,k)= d9(0)*rhov_n(i,j,k)   &
                  + d9(1)*rhov_n(i+1,j,k)+d9(2)*rhov_n(i+2,j,k)+d9(3)*rhov_n(i+3,j,k)+d9(4)*rhov_n(i+4,j,k) &
                  + d9(1)*rhov_n(i-1,j,k)+d9(2)*rhov_n(i-2,j,k)+d9(3)*rhov_n(i-3,j,k)+d9(4)*rhov_n(i-4,j,k)
             Krhow(i,j,k)= d9(0)*rhow_n(i,j,k)   &
                  + d9(1)*rhow_n(i+1,j,k)+d9(2)*rhow_n(i+2,j,k)+d9(3)*rhow_n(i+3,j,k)+d9(4)*rhow_n(i+4,j,k) &
                  + d9(1)*rhow_n(i-1,j,k)+d9(2)*rhow_n(i-2,j,k)+d9(3)*rhow_n(i-3,j,k)+d9(4)*rhow_n(i-4,j,k)
             Krhoe(i,j,k)= d9(0)*rhoe_n(i,j,k)   &
                  + d9(1)*rhoe_n(i+1,j,k)+d9(2)*rhoe_n(i+2,j,k)+d9(3)*rhoe_n(i+3,j,k)+d9(4)*rhoe_n(i+4,j,k) &
                  + d9(1)*rhoe_n(i-1,j,k)+d9(2)*rhoe_n(i-2,j,k)+d9(3)*rhoe_n(i-3,j,k)+d9(4)*rhoe_n(i-4,j,k)
          enddo
       enddo
       ! 6-th order centered filter
       i=nx-3
       do k=1,nz
          do j=1,ny
             Krho(i,j,k) = d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i+1,j,k)+d7(2)*rho_n(i+2,j,k)+d7(3)*rho_n(i+3,j,k) &
                  + d7(1)*rho_n(i-1,j,k)+d7(2)*rho_n(i-2,j,k)+d7(3)*rho_n(i-3,j,k)
             Krhou(i,j,k)= d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i+1,j,k)+d7(2)*rhou_n(i+2,j,k)+d7(3)*rhou_n(i+3,j,k) &
                  + d7(1)*rhou_n(i-1,j,k)+d7(2)*rhou_n(i-2,j,k)+d7(3)*rhou_n(i-3,j,k)
             Krhov(i,j,k)= d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i+1,j,k)+d7(2)*rhov_n(i+2,j,k)+d7(3)*rhov_n(i+3,j,k) &
                  + d7(1)*rhov_n(i-1,j,k)+d7(2)*rhov_n(i-2,j,k)+d7(3)*rhov_n(i-3,j,k)
             Krhow(i,j,k)= d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i+1,j,k)+d7(2)*rhow_n(i+2,j,k)+d7(3)*rhow_n(i+3,j,k) &
                  + d7(1)*rhow_n(i-1,j,k)+d7(2)*rhow_n(i-2,j,k)+d7(3)*rhow_n(i-3,j,k)
             Krhoe(i,j,k)= d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i+1,j,k)+d7(2)*rhoe_n(i+2,j,k)+d7(3)*rhoe_n(i+3,j,k) &
                  + d7(1)*rhoe_n(i-1,j,k)+d7(2)*rhoe_n(i-2,j,k)+d7(3)*rhoe_n(i-3,j,k)
          enddo
       enddo
       ! 4-th order centered filter
       i=nx-2
       do k=1,nz
          do j=1,ny
             Krho(i,j,k) = d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i+1,j,k)+d5(2)*rho_n(i+2,j,k) &
                  + d5(1)*rho_n(i-1,j,k)+d5(2)*rho_n(i-2,j,k)
             Krhou(i,j,k)= d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i+1,j,k)+d5(2)*rhou_n(i+2,j,k) &
                  + d5(1)*rhou_n(i-1,j,k)+d5(2)*rhou_n(i-2,j,k)
             Krhov(i,j,k)= d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i+1,j,k)+d5(2)*rhov_n(i+2,j,k) &
                  + d5(1)*rhov_n(i-1,j,k)+d5(2)*rhov_n(i-2,j,k)
             Krhow(i,j,k)= d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i+1,j,k)+d5(2)*rhow_n(i+2,j,k) &
                  + d5(1)*rhow_n(i-1,j,k)+d5(2)*rhow_n(i-2,j,k)
             Krhoe(i,j,k)= d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i+1,j,k)+d5(2)*rhoe_n(i+2,j,k) &
                  + d5(1)*rhoe_n(i-1,j,k)+d5(2)*rhoe_n(i-2,j,k)
          enddo
       enddo
       ! no filter for i=nx-1 & i=nx -> initialize array
       do j=1,ny
           Krho(nx-1,j,:)=0.0_wp
          Krhou(nx-1,j,:)=0.0_wp
          Krhov(nx-1,j,:)=0.0_wp
          Krhow(nx-1,j,:)=0.0_wp
          Krhoe(nx-1,j,:)=0.0_wp
          
           Krho(nx,j,:)=0.0_wp
          Krhou(nx,j,:)=0.0_wp
          Krhov(nx,j,:)=0.0_wp
          Krhow(nx,j,:)=0.0_wp
          Krhoe(nx,j,:)=0.0_wp
       enddo
    else
       ! 10-th order centered filter
       do k=1,nz
          do j=1,ny
             do i=nx-4,nx
                Krho(i,j,k) = d11(0)*rho_n(i,j,k) &
                     + d11(1)*( rho_n(i+1,j,k)+rho_n(i-1,j,k) )&
                     + d11(2)*( rho_n(i+2,j,k)+rho_n(i-2,j,k) )&
                     + d11(3)*( rho_n(i+3,j,k)+rho_n(i-3,j,k) )&
                     + d11(4)*( rho_n(i+4,j,k)+rho_n(i-4,j,k) )&
                     + d11(5)*( rho_n(i+5,j,k)+rho_n(i-5,j,k) )

                Krhou(i,j,k)= d11(0)*rhou_n(i,j,k) &
                     + d11(1)*( rhou_n(i+1,j,k)+rhou_n(i-1,j,k) )&
                     + d11(2)*( rhou_n(i+2,j,k)+rhou_n(i-2,j,k) )&
                     + d11(3)*( rhou_n(i+3,j,k)+rhou_n(i-3,j,k) )&
                     + d11(4)*( rhou_n(i+4,j,k)+rhou_n(i-4,j,k) )&
                     + d11(5)*( rhou_n(i+5,j,k)+rhou_n(i-5,j,k) )

                Krhov(i,j,k)= d11(0)*rhov_n(i,j,k) &
                     + d11(1)*( rhov_n(i+1,j,k)+rhov_n(i-1,j,k) )&
                     + d11(2)*( rhov_n(i+2,j,k)+rhov_n(i-2,j,k) )&
                     + d11(3)*( rhov_n(i+3,j,k)+rhov_n(i-3,j,k) )&
                     + d11(4)*( rhov_n(i+4,j,k)+rhov_n(i-4,j,k) )&
                     + d11(5)*( rhov_n(i+5,j,k)+rhov_n(i-5,j,k) )

                Krhow(i,j,k)= d11(0)*rhow_n(i,j,k) &
                     + d11(1)*( rhow_n(i+1,j,k)+rhow_n(i-1,j,k) )&
                     + d11(2)*( rhow_n(i+2,j,k)+rhow_n(i-2,j,k) )&
                     + d11(3)*( rhow_n(i+3,j,k)+rhow_n(i-3,j,k) )&
                     + d11(4)*( rhow_n(i+4,j,k)+rhow_n(i-4,j,k) )&
                     + d11(5)*( rhow_n(i+5,j,k)+rhow_n(i-5,j,k) )

                Krhoe(i,j,k)= d11(0)*rhoe_n(i,j,k) &
                     + d11(1)*( rhoe_n(i+1,j,k)+rhoe_n(i-1,j,k) )&
                     + d11(2)*( rhoe_n(i+2,j,k)+rhoe_n(i-2,j,k) )&
                     + d11(3)*( rhoe_n(i+3,j,k)+rhoe_n(i-3,j,k) )&
                     + d11(4)*( rhoe_n(i+4,j,k)+rhoe_n(i-4,j,k) )&
                     + d11(5)*( rhoe_n(i+5,j,k)+rhoe_n(i-5,j,k) )
             enddo
          enddo
       enddo
    endif
    
    ! Filtering along j-direction
    ! ===========================
    
    ! Perform one-sided communications in j-direction
    ! -----------------------------------------------
    ! open target windows
    call window_comm1
    ! communicate var
    call communication_1_j
    
    ! Interior points
    ! ---------------
    ! 10-th order centered filter
    do k=1,nz
       do j=6,ny-5
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  + d11(0)*rho_n(i,j,k) &
                  + d11(1)*( rho_n(i,j+1,k)+rho_n(i,j-1,k) )&
                  + d11(2)*( rho_n(i,j+2,k)+rho_n(i,j-2,k) )&
                  + d11(3)*( rho_n(i,j+3,k)+rho_n(i,j-3,k) )&
                  + d11(4)*( rho_n(i,j+4,k)+rho_n(i,j-4,k) )&
                  + d11(5)*( rho_n(i,j+5,k)+rho_n(i,j-5,k) )

             Krhou(i,j,k)= Krhou(i,j,k)  + d11(0)*rhou_n(i,j,k) &
                  + d11(1)*( rhou_n(i,j+1,k)+rhou_n(i,j-1,k) )&
                  + d11(2)*( rhou_n(i,j+2,k)+rhou_n(i,j-2,k) )&
                  + d11(3)*( rhou_n(i,j+3,k)+rhou_n(i,j-3,k) )&
                  + d11(4)*( rhou_n(i,j+4,k)+rhou_n(i,j-4,k) )&
                  + d11(5)*( rhou_n(i,j+5,k)+rhou_n(i,j-5,k) )

             Krhov(i,j,k)= Krhov(i,j,k)  + d11(0)*rhov_n(i,j,k) &
                  + d11(1)*( rhov_n(i,j+1,k)+rhov_n(i,j-1,k) )&
                  + d11(2)*( rhov_n(i,j+2,k)+rhov_n(i,j-2,k) )&
                  + d11(3)*( rhov_n(i,j+3,k)+rhov_n(i,j-3,k) )&
                  + d11(4)*( rhov_n(i,j+4,k)+rhov_n(i,j-4,k) )&
                  + d11(5)*( rhov_n(i,j+5,k)+rhov_n(i,j-5,k) )

             Krhow(i,j,k)= Krhow(i,j,k)  + d11(0)*rhow_n(i,j,k) &
                  + d11(1)*( rhow_n(i,j+1,k)+rhow_n(i,j-1,k) )&
                  + d11(2)*( rhow_n(i,j+2,k)+rhow_n(i,j-2,k) )&
                  + d11(3)*( rhow_n(i,j+3,k)+rhow_n(i,j-3,k) )&
                  + d11(4)*( rhow_n(i,j+4,k)+rhow_n(i,j-4,k) )&
                  + d11(5)*( rhow_n(i,j+5,k)+rhow_n(i,j-5,k) )

             Krhoe(i,j,k)= Krhoe(i,j,k)  + d11(0)*rhoe_n(i,j,k) &
                  + d11(1)*( rhoe_n(i,j+1,k)+rhoe_n(i,j-1,k) )&
                  + d11(2)*( rhoe_n(i,j+2,k)+rhoe_n(i,j-2,k) )&
                  + d11(3)*( rhoe_n(i,j+3,k)+rhoe_n(i,j-3,k) )&
                  + d11(4)*( rhoe_n(i,j+4,k)+rhoe_n(i,j-4,k) )&
                  + d11(5)*( rhoe_n(i,j+5,k)+rhoe_n(i,j-5,k) )
          enddo
       enddo
    enddo
    
    ! Check that communications are fullfilled
    ! ----------------------------------------
    call window_comm1

    ! BC at jmin
    ! ----------
    if (is_boundary(2,1)) then
       ! no filter for j=1 & j=2
       
       ! 4-th order centered filter
       j=3
       do k=1,nz
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i,j+1,k)+d5(2)*rho_n(i,j+2,k) &
                  + d5(1)*rho_n(i,j-1,k)+d5(2)*rho_n(i,j-2,k)
             Krhou(i,j,k)= Krhou(i,j,k)  +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j+1,k)+d5(2)*rhou_n(i,j+2,k) &
                  + d5(1)*rhou_n(i,j-1,k)+d5(2)*rhou_n(i,j-2,k)
             Krhov(i,j,k)= Krhov(i,j,k)  +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j+1,k)+d5(2)*rhov_n(i,j+2,k) &
                  + d5(1)*rhov_n(i,j-1,k)+d5(2)*rhov_n(i,j-2,k)
             Krhow(i,j,k)= Krhow(i,j,k)  +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j+1,k)+d5(2)*rhow_n(i,j+2,k) &
                  + d5(1)*rhow_n(i,j-1,k)+d5(2)*rhow_n(i,j-2,k)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i,j+1,k)+d5(2)*rhoe_n(i,j+2,k) &
                  + d5(1)*rhoe_n(i,j-1,k)+d5(2)*rhoe_n(i,j-2,k)
          enddo
       enddo
       ! 6-th order centered filter
       j=4
       do k=1,nz
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i,j+1,k)+d7(2)*rho_n(i,j+2,k)+d7(3)*rho_n(i,j+3,k) &
                  + d7(1)*rho_n(i,j-1,k)+d7(2)*rho_n(i,j-2,k)+d7(3)*rho_n(i,j-3,k)
             Krhou(i,j,k)= Krhou(i,j,k)  +d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i,j+1,k)+d7(2)*rhou_n(i,j+2,k)+d7(3)*rhou_n(i,j+3,k) &
                  + d7(1)*rhou_n(i,j-1,k)+d7(2)*rhou_n(i,j-2,k)+d7(3)*rhou_n(i,j-3,k)
             Krhov(i,j,k)= Krhov(i,j,k)  +d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i,j+1,k)+d7(2)*rhov_n(i,j+2,k)+d7(3)*rhov_n(i,j+3,k) &
                  + d7(1)*rhov_n(i,j-1,k)+d7(2)*rhov_n(i,j-2,k)+d7(3)*rhov_n(i,j-3,k)
             Krhow(i,j,k)= Krhow(i,j,k)  +d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i,j+1,k)+d7(2)*rhow_n(i,j+2,k)+d7(3)*rhow_n(i,j+3,k) &
                  + d7(1)*rhow_n(i,j-1,k)+d7(2)*rhow_n(i,j-2,k)+d7(3)*rhow_n(i,j-3,k)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i,j+1,k)+d7(2)*rhoe_n(i,j+2,k)+d7(3)*rhoe_n(i,j+3,k) &
                  + d7(1)*rhoe_n(i,j-1,k)+d7(2)*rhoe_n(i,j-2,k)+d7(3)*rhoe_n(i,j-3,k)
          enddo
       enddo
       ! 8-th order centered filter
       j=5
       do k=1,nz
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d9(0)*rho_n(i,j,k)   &
                  + d9(1)*rho_n(i,j+1,k)+d9(2)*rho_n(i,j+2,k)+d9(3)*rho_n(i,j+3,k)+d9(4)*rho_n(i,j+4,k) &
                  + d9(1)*rho_n(i,j-1,k)+d9(2)*rho_n(i,j-2,k)+d9(3)*rho_n(i,j-3,k)+d9(4)*rho_n(i,j-4,k)
             Krhou(i,j,k)= Krhou(i,j,k)  +d9(0)*rhou_n(i,j,k)   &
                  + d9(1)*rhou_n(i,j+1,k)+d9(2)*rhou_n(i,j+2,k)+d9(3)*rhou_n(i,j+3,k)+d9(4)*rhou_n(i,j+4,k) &
                  + d9(1)*rhou_n(i,j-1,k)+d9(2)*rhou_n(i,j-2,k)+d9(3)*rhou_n(i,j-3,k)+d9(4)*rhou_n(i,j-4,k)
             Krhov(i,j,k)= Krhov(i,j,k)  +d9(0)*rhov_n(i,j,k)   &
                  + d9(1)*rhov_n(i,j+1,k)+d9(2)*rhov_n(i,j+2,k)+d9(3)*rhov_n(i,j+3,k)+d9(4)*rhov_n(i,j+4,k) &
                  + d9(1)*rhov_n(i,j-1,k)+d9(2)*rhov_n(i,j-2,k)+d9(3)*rhov_n(i,j-3,k)+d9(4)*rhov_n(i,j-4,k)
             Krhow(i,j,k)= Krhow(i,j,k)  +d9(0)*rhow_n(i,j,k)   &
                  + d9(1)*rhow_n(i,j+1,k)+d9(2)*rhow_n(i,j+2,k)+d9(3)*rhow_n(i,j+3,k)+d9(4)*rhow_n(i,j+4,k) &
                  + d9(1)*rhow_n(i,j-1,k)+d9(2)*rhow_n(i,j-2,k)+d9(3)*rhow_n(i,j-3,k)+d9(4)*rhow_n(i,j-4,k)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d9(0)*rhoe_n(i,j,k)   &
                  + d9(1)*rhoe_n(i,j+1,k)+d9(2)*rhoe_n(i,j+2,k)+d9(3)*rhoe_n(i,j+3,k)+d9(4)*rhoe_n(i,j+4,k) &
                  + d9(1)*rhoe_n(i,j-1,k)+d9(2)*rhoe_n(i,j-2,k)+d9(3)*rhoe_n(i,j-3,k)+d9(4)*rhoe_n(i,j-4,k)
          enddo
       enddo
    else
       ! 10-th order centered filter
       do k=1,nz
          do j=1,5
             do i=1,nx
                Krho(i,j,k) = Krho(i,j,k)  + d11(0)*rho_n(i,j,k) &
                     + d11(1)*( rho_n(i,j+1,k)+rho_n(i,j-1,k) )&
                     + d11(2)*( rho_n(i,j+2,k)+rho_n(i,j-2,k) )&
                     + d11(3)*( rho_n(i,j+3,k)+rho_n(i,j-3,k) )&
                     + d11(4)*( rho_n(i,j+4,k)+rho_n(i,j-4,k) )&
                     + d11(5)*( rho_n(i,j+5,k)+rho_n(i,j-5,k) )

                Krhou(i,j,k)= Krhou(i,j,k)  + d11(0)*rhou_n(i,j,k) &
                     + d11(1)*( rhou_n(i,j+1,k)+rhou_n(i,j-1,k) )&
                     + d11(2)*( rhou_n(i,j+2,k)+rhou_n(i,j-2,k) )&
                     + d11(3)*( rhou_n(i,j+3,k)+rhou_n(i,j-3,k) )&
                     + d11(4)*( rhou_n(i,j+4,k)+rhou_n(i,j-4,k) )&
                     + d11(5)*( rhou_n(i,j+5,k)+rhou_n(i,j-5,k) )

                Krhov(i,j,k)= Krhov(i,j,k)  + d11(0)*rhov_n(i,j,k) &
                     + d11(1)*( rhov_n(i,j+1,k)+rhov_n(i,j-1,k) )&
                     + d11(2)*( rhov_n(i,j+2,k)+rhov_n(i,j-2,k) )&
                     + d11(3)*( rhov_n(i,j+3,k)+rhov_n(i,j-3,k) )&
                     + d11(4)*( rhov_n(i,j+4,k)+rhov_n(i,j-4,k) )&
                     + d11(5)*( rhov_n(i,j+5,k)+rhov_n(i,j-5,k) )

                Krhow(i,j,k)= Krhow(i,j,k)  + d11(0)*rhow_n(i,j,k) &
                     + d11(1)*( rhow_n(i,j+1,k)+rhow_n(i,j-1,k) )&
                     + d11(2)*( rhow_n(i,j+2,k)+rhow_n(i,j-2,k) )&
                     + d11(3)*( rhow_n(i,j+3,k)+rhow_n(i,j-3,k) )&
                     + d11(4)*( rhow_n(i,j+4,k)+rhow_n(i,j-4,k) )&
                     + d11(5)*( rhow_n(i,j+5,k)+rhow_n(i,j-5,k) )

                Krhoe(i,j,k)= Krhoe(i,j,k)  + d11(0)*rhoe_n(i,j,k) &
                     + d11(1)*( rhoe_n(i,j+1,k)+rhoe_n(i,j-1,k) )&
                     + d11(2)*( rhoe_n(i,j+2,k)+rhoe_n(i,j-2,k) )&
                     + d11(3)*( rhoe_n(i,j+3,k)+rhoe_n(i,j-3,k) )&
                     + d11(4)*( rhoe_n(i,j+4,k)+rhoe_n(i,j-4,k) )&
                     + d11(5)*( rhoe_n(i,j+5,k)+rhoe_n(i,j-5,k) )
             enddo
          enddo
       enddo
    endif

    ! BC at jmax
    ! ----------
    if (is_boundary(2,2)) then
       ! 8-th order centered filter
       j=ny-4
       do k=1,nz
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d9(0)*rho_n(i,j,k)   &
                  + d9(1)*rho_n(i,j+1,k)+d9(2)*rho_n(i,j+2,k)+d9(3)*rho_n(i,j+3,k)+d9(4)*rho_n(i,j+4,k) &
                  + d9(1)*rho_n(i,j-1,k)+d9(2)*rho_n(i,j-2,k)+d9(3)*rho_n(i,j-3,k)+d9(4)*rho_n(i,j-4,k)
             Krhou(i,j,k)= Krhou(i,j,k)  +d9(0)*rhou_n(i,j,k)   &
                  + d9(1)*rhou_n(i,j+1,k)+d9(2)*rhou_n(i,j+2,k)+d9(3)*rhou_n(i,j+3,k)+d9(4)*rhou_n(i,j+4,k) &
                  + d9(1)*rhou_n(i,j-1,k)+d9(2)*rhou_n(i,j-2,k)+d9(3)*rhou_n(i,j-3,k)+d9(4)*rhou_n(i,j-4,k)
             Krhov(i,j,k)= Krhov(i,j,k)  +d9(0)*rhov_n(i,j,k)   &
                  + d9(1)*rhov_n(i,j+1,k)+d9(2)*rhov_n(i,j+2,k)+d9(3)*rhov_n(i,j+3,k)+d9(4)*rhov_n(i,j+4,k) &
                  + d9(1)*rhov_n(i,j-1,k)+d9(2)*rhov_n(i,j-2,k)+d9(3)*rhov_n(i,j-3,k)+d9(4)*rhov_n(i,j-4,k)
             Krhow(i,j,k)= Krhow(i,j,k)  +d9(0)*rhow_n(i,j,k)   &
                  + d9(1)*rhow_n(i,j+1,k)+d9(2)*rhow_n(i,j+2,k)+d9(3)*rhow_n(i,j+3,k)+d9(4)*rhow_n(i,j+4,k) &
                  + d9(1)*rhow_n(i,j-1,k)+d9(2)*rhow_n(i,j-2,k)+d9(3)*rhow_n(i,j-3,k)+d9(4)*rhow_n(i,j-4,k)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d9(0)*rhoe_n(i,j,k)   &
                  + d9(1)*rhoe_n(i,j+1,k)+d9(2)*rhoe_n(i,j+2,k)+d9(3)*rhoe_n(i,j+3,k)+d9(4)*rhoe_n(i,j+4,k) &
                  + d9(1)*rhoe_n(i,j-1,k)+d9(2)*rhoe_n(i,j-2,k)+d9(3)*rhoe_n(i,j-3,k)+d9(4)*rhoe_n(i,j-4,k)
          enddo
       enddo
       ! 6-th order centered filter
       j=ny-3
       do k=1,nz
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i,j+1,k)+d7(2)*rho_n(i,j+2,k)+d7(3)*rho_n(i,j+3,k) &
                  + d7(1)*rho_n(i,j-1,k)+d7(2)*rho_n(i,j-2,k)+d7(3)*rho_n(i,j-3,k)
             Krhou(i,j,k)= Krhou(i,j,k)  +d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i,j+1,k)+d7(2)*rhou_n(i,j+2,k)+d7(3)*rhou_n(i,j+3,k) &
                  + d7(1)*rhou_n(i,j-1,k)+d7(2)*rhou_n(i,j-2,k)+d7(3)*rhou_n(i,j-3,k)
             Krhov(i,j,k)= Krhov(i,j,k)  +d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i,j+1,k)+d7(2)*rhov_n(i,j+2,k)+d7(3)*rhov_n(i,j+3,k) &
                  + d7(1)*rhov_n(i,j-1,k)+d7(2)*rhov_n(i,j-2,k)+d7(3)*rhov_n(i,j-3,k)
             Krhow(i,j,k)= Krhow(i,j,k)  +d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i,j+1,k)+d7(2)*rhow_n(i,j+2,k)+d7(3)*rhow_n(i,j+3,k) &
                  + d7(1)*rhow_n(i,j-1,k)+d7(2)*rhow_n(i,j-2,k)+d7(3)*rhow_n(i,j-3,k)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i,j+1,k)+d7(2)*rhoe_n(i,j+2,k)+d7(3)*rhoe_n(i,j+3,k) &
                  + d7(1)*rhoe_n(i,j-1,k)+d7(2)*rhoe_n(i,j-2,k)+d7(3)*rhoe_n(i,j-3,k)
          enddo
       enddo
       ! 4-th order centered filter
       j=ny-2
       do k=1,nz
          do i=1,nx
             Krho(i,j,k)= Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i,j+1,k)+d5(2)*rho_n(i,j+2,k) &
                  + d5(1)*rho_n(i,j-1,k)+d5(2)*rho_n(i,j-2,k)
             Krhou(i,j,k)= Krhou(i,j,k)  +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j+1,k)+d5(2)*rhou_n(i,j+2,k) &
                  + d5(1)*rhou_n(i,j-1,k)+d5(2)*rhou_n(i,j-2,k)
             Krhov(i,j,k)= Krhov(i,j,k)  +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j+1,k)+d5(2)*rhov_n(i,j+2,k) &
                  + d5(1)*rhov_n(i,j-1,k)+d5(2)*rhov_n(i,j-2,k)
             Krhow(i,j,k)= Krhow(i,j,k)  +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j+1,k)+d5(2)*rhow_n(i,j+2,k) &
                  + d5(1)*rhow_n(i,j-1,k)+d5(2)*rhow_n(i,j-2,k)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i,j+1,k)+d5(2)*rhoe_n(i,j+2,k) &
                  + d5(1)*rhoe_n(i,j-1,k)+d5(2)*rhoe_n(i,j-2,k)
          enddo
       enddo
       ! no filter for j=ny-1 & j=ny
    else
       ! 10-th order centered filter
       do k=1,nz
          do j=ny-4,ny
             do i=1,nx
                Krho(i,j,k) = Krho(i,j,k)  + d11(0)*rho_n(i,j,k) &
                     + d11(1)*( rho_n(i,j+1,k)+rho_n(i,j-1,k) )&
                     + d11(2)*( rho_n(i,j+2,k)+rho_n(i,j-2,k) )&
                     + d11(3)*( rho_n(i,j+3,k)+rho_n(i,j-3,k) )&
                     + d11(4)*( rho_n(i,j+4,k)+rho_n(i,j-4,k) )&
                     + d11(5)*( rho_n(i,j+5,k)+rho_n(i,j-5,k) )

                Krhou(i,j,k)= Krhou(i,j,k)  + d11(0)*rhou_n(i,j,k) &
                     + d11(1)*( rhou_n(i,j+1,k)+rhou_n(i,j-1,k) )&
                     + d11(2)*( rhou_n(i,j+2,k)+rhou_n(i,j-2,k) )&
                     + d11(3)*( rhou_n(i,j+3,k)+rhou_n(i,j-3,k) )&
                     + d11(4)*( rhou_n(i,j+4,k)+rhou_n(i,j-4,k) )&
                     + d11(5)*( rhou_n(i,j+5,k)+rhou_n(i,j-5,k) )

                Krhov(i,j,k)= Krhov(i,j,k)  + d11(0)*rhov_n(i,j,k) &
                     + d11(1)*( rhov_n(i,j+1,k)+rhov_n(i,j-1,k) )&
                     + d11(2)*( rhov_n(i,j+2,k)+rhov_n(i,j-2,k) )&
                     + d11(3)*( rhov_n(i,j+3,k)+rhov_n(i,j-3,k) )&
                     + d11(4)*( rhov_n(i,j+4,k)+rhov_n(i,j-4,k) )&
                     + d11(5)*( rhov_n(i,j+5,k)+rhov_n(i,j-5,k) )

                Krhow(i,j,k)= Krhow(i,j,k)  + d11(0)*rhow_n(i,j,k) &
                     + d11(1)*( rhow_n(i,j+1,k)+rhow_n(i,j-1,k) )&
                     + d11(2)*( rhow_n(i,j+2,k)+rhow_n(i,j-2,k) )&
                     + d11(3)*( rhow_n(i,j+3,k)+rhow_n(i,j-3,k) )&
                     + d11(4)*( rhow_n(i,j+4,k)+rhow_n(i,j-4,k) )&
                     + d11(5)*( rhow_n(i,j+5,k)+rhow_n(i,j-5,k) )

                Krhoe(i,j,k)= Krhoe(i,j,k)  + d11(0)*rhoe_n(i,j,k) &
                     + d11(1)*( rhoe_n(i,j+1,k)+rhoe_n(i,j-1,k) )&
                     + d11(2)*( rhoe_n(i,j+2,k)+rhoe_n(i,j-2,k) )&
                     + d11(3)*( rhoe_n(i,j+3,k)+rhoe_n(i,j-3,k) )&
                     + d11(4)*( rhoe_n(i,j+4,k)+rhoe_n(i,j-4,k) )&
                     + d11(5)*( rhoe_n(i,j+5,k)+rhoe_n(i,j-5,k) )
             enddo
          enddo
       enddo
    endif

    !****************
    if (is_2D) return
    !****************
    
    ! Filtering along k-direction
    ! ===========================

    ! Perform one-sided communications in k-direction
    ! -----------------------------------------------
    ! open target windows
    call window_comm1
    ! communicate var
    call communication_1_k
    
    ! Interior points
    ! ---------------
    ! 10-th order centered filter
    do k=6,nz-5
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  + d11(0)*rho_n(i,j,k) &
                  + d11(1)*( rho_n(i,j,k+1)+rho_n(i,j,k-1) )&
                  + d11(2)*( rho_n(i,j,k+2)+rho_n(i,j,k-2) )&
                  + d11(3)*( rho_n(i,j,k+3)+rho_n(i,j,k-3) )&
                  + d11(4)*( rho_n(i,j,k+4)+rho_n(i,j,k-4) )&
                  + d11(5)*( rho_n(i,j,k+5)+rho_n(i,j,k-5) )

             Krhou(i,j,k)= Krhou(i,j,k)  + d11(0)*rhou_n(i,j,k) &
                  + d11(1)*( rhou_n(i,j,k+1)+rhou_n(i,j,k-1) )&
                  + d11(2)*( rhou_n(i,j,k+2)+rhou_n(i,j,k-2) )&
                  + d11(3)*( rhou_n(i,j,k+3)+rhou_n(i,j,k-3) )&
                  + d11(4)*( rhou_n(i,j,k+4)+rhou_n(i,j,k-4) )&
                  + d11(5)*( rhou_n(i,j,k+5)+rhou_n(i,j,k-5) )

             Krhov(i,j,k)= Krhov(i,j,k)  + d11(0)*rhov_n(i,j,k) &
                  + d11(1)*( rhov_n(i,j,k+1)+rhov_n(i,j,k-1) )&
                  + d11(2)*( rhov_n(i,j,k+2)+rhov_n(i,j,k-2) )&
                  + d11(3)*( rhov_n(i,j,k+3)+rhov_n(i,j,k-3) )&
                  + d11(4)*( rhov_n(i,j,k+4)+rhov_n(i,j,k-4) )&
                  + d11(5)*( rhov_n(i,j,k+5)+rhov_n(i,j,k-5) )

             Krhow(i,j,k)= Krhow(i,j,k)  + d11(0)*rhow_n(i,j,k) &
                  + d11(1)*( rhow_n(i,j,k+1)+rhow_n(i,j,k-1) )&
                  + d11(2)*( rhow_n(i,j,k+2)+rhow_n(i,j,k-2) )&
                  + d11(3)*( rhow_n(i,j,k+3)+rhow_n(i,j,k-3) )&
                  + d11(4)*( rhow_n(i,j,k+4)+rhow_n(i,j,k-4) )&
                  + d11(5)*( rhow_n(i,j,k+5)+rhow_n(i,j,k-5) )

             Krhoe(i,j,k)= Krhoe(i,j,k)  + d11(0)*rhoe_n(i,j,k) &
                  + d11(1)*( rhoe_n(i,j,k+1)+rhoe_n(i,j,k-1) )&
                  + d11(2)*( rhoe_n(i,j,k+2)+rhoe_n(i,j,k-2) )&
                  + d11(3)*( rhoe_n(i,j,k+3)+rhoe_n(i,j,k-3) )&
                  + d11(4)*( rhoe_n(i,j,k+4)+rhoe_n(i,j,k-4) )&
                  + d11(5)*( rhoe_n(i,j,k+5)+rhoe_n(i,j,k-5) )
          enddo
       enddo
    enddo
    
    ! Check that communications are fullfilled
    ! ----------------------------------------
    call window_comm1

    ! BC at kmin
    ! ----------
    if (is_boundary(3,1)) then
       ! no filter for k=1 & k=2

       ! 4-th order centered filter
       k=3
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i,j,k+1)+d5(2)*rho_n(i,j,k+2) &
                  + d5(1)*rho_n(i,j,k-1)+d5(2)*rho_n(i,j,k-2)
             Krhou(i,j,k)= Krhou(i,j,k)  +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j,k+1)+d5(2)*rhou_n(i,j,k+2) &
                  + d5(1)*rhou_n(i,j,k-1)+d5(2)*rhou_n(i,j,k-2)
             Krhov(i,j,k)= Krhov(i,j,k)  +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j,k+1)+d5(2)*rhov_n(i,j,k+2) &
                  + d5(1)*rhov_n(i,j,k-1)+d5(2)*rhov_n(i,j,k-2)
             Krhow(i,j,k)= Krhow(i,j,k)  +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j,k+1)+d5(2)*rhow_n(i,j,k+2) &
                  + d5(1)*rhow_n(i,j,k-1)+d5(2)*rhow_n(i,j,k-2)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i,j,k+1)+d5(2)*rhoe_n(i,j,k+2) &
                  + d5(1)*rhoe_n(i,j,k-1)+d5(2)*rhoe_n(i,j,k-2)
          enddo
       enddo
       ! 6-th order centered filter
       k=4
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i,j,k+1)+d7(2)*rho_n(i,j,k+2)+d7(3)*rho_n(i,j,k+3) &
                  + d7(1)*rho_n(i,j,k-1)+d7(2)*rho_n(i,j,k-2)+d7(3)*rho_n(i,j,k-3)
             Krhou(i,j,k)= Krhou(i,j,k)  +d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i,j,k+1)+d7(2)*rhou_n(i,j,k+2)+d7(3)*rhou_n(i,j,k+3) &
                  + d7(1)*rhou_n(i,j,k-1)+d7(2)*rhou_n(i,j,k-2)+d7(3)*rhou_n(i,j,k-3)
             Krhov(i,j,k)= Krhov(i,j,k)  +d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i,j,k+1)+d7(2)*rhov_n(i,j,k+2)+d7(3)*rhov_n(i,j,k+3) &
                  + d7(1)*rhov_n(i,j,k-1)+d7(2)*rhov_n(i,j,k-2)+d7(3)*rhov_n(i,j,k-3)
             Krhow(i,j,k)= Krhow(i,j,k)  +d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i,j,k+1)+d7(2)*rhow_n(i,j,k+2)+d7(3)*rhow_n(i,j,k+3) &
                  + d7(1)*rhow_n(i,j,k-1)+d7(2)*rhow_n(i,j,k-2)+d7(3)*rhow_n(i,j,k-3)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i,j,k+1)+d7(2)*rhoe_n(i,j,k+2)+d7(3)*rhoe_n(i,j,k+3) &
                  + d7(1)*rhoe_n(i,j,k-1)+d7(2)*rhoe_n(i,j,k-2)+d7(3)*rhoe_n(i,j,k-3)
          enddo
       enddo
       ! 8-th order centered filter
       k=5
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d9(0)*rho_n(i,j,k)   &
                  + d9(1)*rho_n(i,j,k+1)+d9(2)*rho_n(i,j,k+2)+d9(3)*rho_n(i,j,k+3)+d9(4)*rho_n(i,j,k+4) &
                  + d9(1)*rho_n(i,j,k-1)+d9(2)*rho_n(i,j,k-2)+d9(3)*rho_n(i,j,k-3)+d9(4)*rho_n(i,j,k-4)
             Krhou(i,j,k)= Krhou(i,j,k)  +d9(0)*rhou_n(i,j,k)   &
                  + d9(1)*rhou_n(i,j,k+1)+d9(2)*rhou_n(i,j,k+2)+d9(3)*rhou_n(i,j,k+3)+d9(4)*rhou_n(i,j,k+4) &
                  + d9(1)*rhou_n(i,j,k-1)+d9(2)*rhou_n(i,j,k-2)+d9(3)*rhou_n(i,j,k-3)+d9(4)*rhou_n(i,j,k-4)
             Krhov(i,j,k)= Krhov(i,j,k)  +d9(0)*rhov_n(i,j,k)   &
                  + d9(1)*rhov_n(i,j,k+1)+d9(2)*rhov_n(i,j,k+2)+d9(3)*rhov_n(i,j,k+3)+d9(4)*rhov_n(i,j,k+4) &
                  + d9(1)*rhov_n(i,j,k-1)+d9(2)*rhov_n(i,j,k-2)+d9(3)*rhov_n(i,j,k-3)+d9(4)*rhov_n(i,j,k-4)
             Krhow(i,j,k)= Krhow(i,j,k)  +d9(0)*rhow_n(i,j,k)   &
                  + d9(1)*rhow_n(i,j,k+1)+d9(2)*rhow_n(i,j,k+2)+d9(3)*rhow_n(i,j,k+3)+d9(4)*rhow_n(i,j,k+4) &
                  + d9(1)*rhow_n(i,j,k-1)+d9(2)*rhow_n(i,j,k-2)+d9(3)*rhow_n(i,j,k-3)+d9(4)*rhow_n(i,j,k-4)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d9(0)*rhoe_n(i,j,k)   &
                  + d9(1)*rhoe_n(i,j,k+1)+d9(2)*rhoe_n(i,j,k+2)+d9(3)*rhoe_n(i,j,k+3)+d9(4)*rhoe_n(i,j,k+4) &
                  + d9(1)*rhoe_n(i,j,k-1)+d9(2)*rhoe_n(i,j,k-2)+d9(3)*rhoe_n(i,j,k-3)+d9(4)*rhoe_n(i,j,k-4)
          enddo
       enddo
    else
       ! 10-th order centered filter
       do k=1,5
          do j=1,ny
             do i=1,nx
                Krho(i,j,k) = Krho(i,j,k)  + d11(0)*rho_n(i,j,k) &
                     + d11(1)*( rho_n(i,j,k+1)+rho_n(i,j,k-1) )&
                     + d11(2)*( rho_n(i,j,k+2)+rho_n(i,j,k-2) )&
                     + d11(3)*( rho_n(i,j,k+3)+rho_n(i,j,k-3) )&
                     + d11(4)*( rho_n(i,j,k+4)+rho_n(i,j,k-4) )&
                     + d11(5)*( rho_n(i,j,k+5)+rho_n(i,j,k-5) )

                Krhou(i,j,k)= Krhou(i,j,k)  + d11(0)*rhou_n(i,j,k) &
                     + d11(1)*( rhou_n(i,j,k+1)+rhou_n(i,j,k-1) )&
                     + d11(2)*( rhou_n(i,j,k+2)+rhou_n(i,j,k-2) )&
                     + d11(3)*( rhou_n(i,j,k+3)+rhou_n(i,j,k-3) )&
                     + d11(4)*( rhou_n(i,j,k+4)+rhou_n(i,j,k-4) )&
                     + d11(5)*( rhou_n(i,j,k+5)+rhou_n(i,j,k-5) )

                Krhov(i,j,k)= Krhov(i,j,k)  + d11(0)*rhov_n(i,j,k) &
                     + d11(1)*( rhov_n(i,j,k+1)+rhov_n(i,j,k-1) )&
                     + d11(2)*( rhov_n(i,j,k+2)+rhov_n(i,j,k-2) )&
                     + d11(3)*( rhov_n(i,j,k+3)+rhov_n(i,j,k-3) )&
                     + d11(4)*( rhov_n(i,j,k+4)+rhov_n(i,j,k-4) )&
                     + d11(5)*( rhov_n(i,j,k+5)+rhov_n(i,j,k-5) )

                Krhow(i,j,k)= Krhow(i,j,k)  + d11(0)*rhow_n(i,j,k) &
                     + d11(1)*( rhow_n(i,j,k+1)+rhow_n(i,j,k-1) )&
                     + d11(2)*( rhow_n(i,j,k+2)+rhow_n(i,j,k-2) )&
                     + d11(3)*( rhow_n(i,j,k+3)+rhow_n(i,j,k-3) )&
                     + d11(4)*( rhow_n(i,j,k+4)+rhow_n(i,j,k-4) )&
                     + d11(5)*( rhow_n(i,j,k+5)+rhow_n(i,j,k-5) )

                Krhoe(i,j,k)= Krhoe(i,j,k)  + d11(0)*rhoe_n(i,j,k) &
                     + d11(1)*( rhoe_n(i,j,k+1)+rhoe_n(i,j,k-1) )&
                     + d11(2)*( rhoe_n(i,j,k+2)+rhoe_n(i,j,k-2) )&
                     + d11(3)*( rhoe_n(i,j,k+3)+rhoe_n(i,j,k-3) )&
                     + d11(4)*( rhoe_n(i,j,k+4)+rhoe_n(i,j,k-4) )&
                     + d11(5)*( rhoe_n(i,j,k+5)+rhoe_n(i,j,k-5) )
             enddo
          enddo
       enddo
    endif

    ! BC at kmax
    ! ----------
    if (is_boundary(3,2)) then
       ! 8-th order centered filter
       k=nz-4
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d9(0)*rho_n(i,j,k)   &
                  + d9(1)*rho_n(i,j,k+1)+d9(2)*rho_n(i,j,k+2)+d9(3)*rho_n(i,j,k+3)+d9(4)*rho_n(i,j,k+4) &
                  + d9(1)*rho_n(i,j,k-1)+d9(2)*rho_n(i,j,k-2)+d9(3)*rho_n(i,j,k-3)+d9(4)*rho_n(i,j,k-4)
             Krhou(i,j,k)= Krhou(i,j,k)  +d9(0)*rhou_n(i,j,k)   &
                  + d9(1)*rhou_n(i,j,k+1)+d9(2)*rhou_n(i,j,k+2)+d9(3)*rhou_n(i,j,k+3)+d9(4)*rhou_n(i,j,k+4) &
                  + d9(1)*rhou_n(i,j,k-1)+d9(2)*rhou_n(i,j,k-2)+d9(3)*rhou_n(i,j,k-3)+d9(4)*rhou_n(i,j,k-4)
             Krhov(i,j,k)= Krhov(i,j,k)  +d9(0)*rhov_n(i,j,k)   &
                  + d9(1)*rhov_n(i,j,k+1)+d9(2)*rhov_n(i,j,k+2)+d9(3)*rhov_n(i,j,k+3)+d9(4)*rhov_n(i,j,k+4) &
                  + d9(1)*rhov_n(i,j,k-1)+d9(2)*rhov_n(i,j,k-2)+d9(3)*rhov_n(i,j,k-3)+d9(4)*rhov_n(i,j,k-4)
             Krhow(i,j,k)= Krhow(i,j,k)  +d9(0)*rhow_n(i,j,k)   &
                  + d9(1)*rhow_n(i,j,k+1)+d9(2)*rhow_n(i,j,k+2)+d9(3)*rhow_n(i,j,k+3)+d9(4)*rhow_n(i,j,k+4) &
                  + d9(1)*rhow_n(i,j,k-1)+d9(2)*rhow_n(i,j,k-2)+d9(3)*rhow_n(i,j,k-3)+d9(4)*rhow_n(i,j,k-4)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d9(0)*rhoe_n(i,j,k)   &
                  + d9(1)*rhoe_n(i,j,k+1)+d9(2)*rhoe_n(i,j,k+2)+d9(3)*rhoe_n(i,j,k+3)+d9(4)*rhoe_n(i,j,k+4) &
                  + d9(1)*rhoe_n(i,j,k-1)+d9(2)*rhoe_n(i,j,k-2)+d9(3)*rhoe_n(i,j,k-3)+d9(4)*rhoe_n(i,j,k-4)
          enddo
       enddo
       ! 6-th order centered filter
       k=nz-3
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i,j,k+1)+d7(2)*rho_n(i,j,k+2)+d7(3)*rho_n(i,j,k+3) &
                  + d7(1)*rho_n(i,j,k-1)+d7(2)*rho_n(i,j,k-2)+d7(3)*rho_n(i,j,k-3)
             Krhou(i,j,k)= Krhou(i,j,k)  +d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i,j,k+1)+d7(2)*rhou_n(i,j,k+2)+d7(3)*rhou_n(i,j,k+3) &
                  + d7(1)*rhou_n(i,j,k-1)+d7(2)*rhou_n(i,j,k-2)+d7(3)*rhou_n(i,j,k-3)
             Krhov(i,j,k)= Krhov(i,j,k)  +d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i,j,k+1)+d7(2)*rhov_n(i,j,k+2)+d7(3)*rhov_n(i,j,k+3) &
                  + d7(1)*rhov_n(i,j,k-1)+d7(2)*rhov_n(i,j,k-2)+d7(3)*rhov_n(i,j,k-3)
             Krhow(i,j,k)= Krhow(i,j,k)  +d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i,j,k+1)+d7(2)*rhow_n(i,j,k+2)+d7(3)*rhow_n(i,j,k+3) &
                  + d7(1)*rhow_n(i,j,k-1)+d7(2)*rhow_n(i,j,k-2)+d7(3)*rhow_n(i,j,k-3)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i,j,k+1)+d7(2)*rhoe_n(i,j,k+2)+d7(3)*rhoe_n(i,j,k+3) &
                  + d7(1)*rhoe_n(i,j,k-1)+d7(2)*rhoe_n(i,j,k-2)+d7(3)*rhoe_n(i,j,k-3)
          enddo
       enddo
       ! 4-th order centered filter
       k=nz-2
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i,j,k+1)+d5(2)*rho_n(i,j,k+2) &
                  + d5(1)*rho_n(i,j,k-1)+d5(2)*rho_n(i,j,k-2)
             Krhou(i,j,k)= Krhou(i,j,k)  +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j,k+1)+d5(2)*rhou_n(i,j,k+2) &
                  + d5(1)*rhou_n(i,j,k-1)+d5(2)*rhou_n(i,j,k-2)
             Krhov(i,j,k)= Krhov(i,j,k)  +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j,k+1)+d5(2)*rhov_n(i,j,k+2) &
                  + d5(1)*rhov_n(i,j,k-1)+d5(2)*rhov_n(i,j,k-2)
             Krhow(i,j,k)= Krhow(i,j,k)  +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j,k+1)+d5(2)*rhow_n(i,j,k+2) &
                  + d5(1)*rhow_n(i,j,k-1)+d5(2)*rhow_n(i,j,k-2)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i,j,k+1)+d5(2)*rhoe_n(i,j,k+2) &
                  + d5(1)*rhoe_n(i,j,k-1)+d5(2)*rhoe_n(i,j,k-2)
          enddo
       enddo
       ! no filter for k=nz-1 & k=nz
    else
       ! 10-th order centered filter
       do k=nz-4,nz
          do j=1,ny
             do i=1,nx
                Krho(i,j,k) = Krho(i,j,k)  + d11(0)*rho_n(i,j,k) &
                     + d11(1)*( rho_n(i,j,k+1)+rho_n(i,j,k-1) )&
                     + d11(2)*( rho_n(i,j,k+2)+rho_n(i,j,k-2) )&
                     + d11(3)*( rho_n(i,j,k+3)+rho_n(i,j,k-3) )&
                     + d11(4)*( rho_n(i,j,k+4)+rho_n(i,j,k-4) )&
                     + d11(5)*( rho_n(i,j,k+5)+rho_n(i,j,k-5) )

                Krhou(i,j,k)= Krhou(i,j,k)  + d11(0)*rhou_n(i,j,k) &
                     + d11(1)*( rhou_n(i,j,k+1)+rhou_n(i,j,k-1) )&
                     + d11(2)*( rhou_n(i,j,k+2)+rhou_n(i,j,k-2) )&
                     + d11(3)*( rhou_n(i,j,k+3)+rhou_n(i,j,k-3) )&
                     + d11(4)*( rhou_n(i,j,k+4)+rhou_n(i,j,k-4) )&
                     + d11(5)*( rhou_n(i,j,k+5)+rhou_n(i,j,k-5) )

                Krhov(i,j,k)= Krhov(i,j,k)  + d11(0)*rhov_n(i,j,k) &
                     + d11(1)*( rhov_n(i,j,k+1)+rhov_n(i,j,k-1) )&
                     + d11(2)*( rhov_n(i,j,k+2)+rhov_n(i,j,k-2) )&
                     + d11(3)*( rhov_n(i,j,k+3)+rhov_n(i,j,k-3) )&
                     + d11(4)*( rhov_n(i,j,k+4)+rhov_n(i,j,k-4) )&
                     + d11(5)*( rhov_n(i,j,k+5)+rhov_n(i,j,k-5) )

                Krhow(i,j,k)= Krhow(i,j,k)  + d11(0)*rhow_n(i,j,k) &
                     + d11(1)*( rhow_n(i,j,k+1)+rhow_n(i,j,k-1) )&
                     + d11(2)*( rhow_n(i,j,k+2)+rhow_n(i,j,k-2) )&
                     + d11(3)*( rhow_n(i,j,k+3)+rhow_n(i,j,k-3) )&
                     + d11(4)*( rhow_n(i,j,k+4)+rhow_n(i,j,k-4) )&
                     + d11(5)*( rhow_n(i,j,k+5)+rhow_n(i,j,k-5) )

                Krhoe(i,j,k)= Krhoe(i,j,k)  + d11(0)*rhoe_n(i,j,k) &
                     + d11(1)*( rhoe_n(i,j,k+1)+rhoe_n(i,j,k-1) )&
                     + d11(2)*( rhoe_n(i,j,k+2)+rhoe_n(i,j,k-2) )&
                     + d11(3)*( rhoe_n(i,j,k+3)+rhoe_n(i,j,k-3) )&
                     + d11(4)*( rhoe_n(i,j,k+4)+rhoe_n(i,j,k-4) )&
                     + d11(5)*( rhoe_n(i,j,k+5)+rhoe_n(i,j,k-5) )
             enddo
          enddo
       enddo
    endif
        
  end subroutine filtering_11pts_overlap
  
  !===============================================================================
  subroutine filtering_9pts_overlap
  !===============================================================================
    !> Apply filtering on 9-point stencils (+ boundaries)
    !> with computation-communication overlap
  !===============================================================================
    use mod_flow
    use mod_filtering
    use mod_comm1
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    ! ---------------------------------------------------------------------------

    ! Filtering along i-direction
    ! ===========================

    ! Perform one-sided communications in i-direction
    ! -----------------------------------------------
    ! open target windows
    call window_comm1
    ! communicate var
    call communication_1_i
    
    ! Interior points
    ! ---------------
    ! 8-th order centered filter
    do k=1,nz
       do j=1,ny
          do i=5,nx-4
             Krho(i,j,k) = d9(0)*rho_n(i,j,k) &
                  + d9(1)*( rho_n(i+1,j,k)+rho_n(i-1,j,k) )&
                  + d9(2)*( rho_n(i+2,j,k)+rho_n(i-2,j,k) )&
                  + d9(3)*( rho_n(i+3,j,k)+rho_n(i-3,j,k) )&
                  + d9(4)*( rho_n(i+4,j,k)+rho_n(i-4,j,k) )

             Krhou(i,j,k)= d9(0)*rhou_n(i,j,k) &
                  + d9(1)*( rhou_n(i+1,j,k)+rhou_n(i-1,j,k) )&
                  + d9(2)*( rhou_n(i+2,j,k)+rhou_n(i-2,j,k) )&
                  + d9(3)*( rhou_n(i+3,j,k)+rhou_n(i-3,j,k) )&
                  + d9(4)*( rhou_n(i+4,j,k)+rhou_n(i-4,j,k) )

             Krhov(i,j,k)= d9(0)*rhov_n(i,j,k) &
                  + d9(1)*( rhov_n(i+1,j,k)+rhov_n(i-1,j,k) )&
                  + d9(2)*( rhov_n(i+2,j,k)+rhov_n(i-2,j,k) )&
                  + d9(3)*( rhov_n(i+3,j,k)+rhov_n(i-3,j,k) )&
                  + d9(4)*( rhov_n(i+4,j,k)+rhov_n(i-4,j,k) )

             Krhow(i,j,k)= d9(0)*rhow_n(i,j,k) &
                  + d9(1)*( rhow_n(i+1,j,k)+rhow_n(i-1,j,k) )&
                  + d9(2)*( rhow_n(i+2,j,k)+rhow_n(i-2,j,k) )&
                  + d9(3)*( rhow_n(i+3,j,k)+rhow_n(i-3,j,k) )&
                  + d9(4)*( rhow_n(i+4,j,k)+rhow_n(i-4,j,k) )

             Krhoe(i,j,k)= d9(0)*rhoe_n(i,j,k) &
                  + d9(1)*( rhoe_n(i+1,j,k)+rhoe_n(i-1,j,k) )&
                  + d9(2)*( rhoe_n(i+2,j,k)+rhoe_n(i-2,j,k) )&
                  + d9(3)*( rhoe_n(i+3,j,k)+rhoe_n(i-3,j,k) )&
                  + d9(4)*( rhoe_n(i+4,j,k)+rhoe_n(i-4,j,k) )
          enddo
       enddo
    enddo
    
    ! Check that communications are fullfilled
    ! ----------------------------------------
    call window_comm1
    
    ! BC at imin
    ! ----------
    if (is_boundary(1,1)) then
       ! no filter for i=1 & i=2 -> initialize array
       do j=1,ny
           Krho(1,j,:)=0.0_wp
          Krhou(1,j,:)=0.0_wp
          Krhov(1,j,:)=0.0_wp
          Krhow(1,j,:)=0.0_wp
          Krhoe(1,j,:)=0.0_wp
          
           Krho(2,j,:)=0.0_wp
          Krhou(2,j,:)=0.0_wp
          Krhov(2,j,:)=0.0_wp
          Krhow(2,j,:)=0.0_wp
          Krhoe(2,j,:)=0.0_wp
       enddo
       ! 4-th order centered filter
       i=3
       do k=1,nz
          do j=1,ny
             Krho(i,j,k) = d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i+1,j,k)+d5(2)*rho_n(i+2,j,k) &
                  + d5(1)*rho_n(i-1,j,k)+d5(2)*rho_n(i-2,j,k)
             Krhou(i,j,k)= d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i+1,j,k)+d5(2)*rhou_n(i+2,j,k) &
                  + d5(1)*rhou_n(i-1,j,k)+d5(2)*rhou_n(i-2,j,k)
             Krhov(i,j,k)= d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i+1,j,k)+d5(2)*rhov_n(i+2,j,k) &
                  + d5(1)*rhov_n(i-1,j,k)+d5(2)*rhov_n(i-2,j,k)
             Krhow(i,j,k)= d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i+1,j,k)+d5(2)*rhow_n(i+2,j,k) &
                  + d5(1)*rhow_n(i-1,j,k)+d5(2)*rhow_n(i-2,j,k)
             Krhoe(i,j,k)= d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i+1,j,k)+d5(2)*rhoe_n(i+2,j,k) &
                  + d5(1)*rhoe_n(i-1,j,k)+d5(2)*rhoe_n(i-2,j,k)
          enddo
       enddo
       ! 6-th order centered filter
       i=4
       do k=1,nz
          do j=1,ny
             Krho(i,j,k) = d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i+1,j,k)+d7(2)*rho_n(i+2,j,k)+d7(3)*rho_n(i+3,j,k) &
                  + d7(1)*rho_n(i-1,j,k)+d7(2)*rho_n(i-2,j,k)+d7(3)*rho_n(i-3,j,k)
             Krhou(i,j,k)= d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i+1,j,k)+d7(2)*rhou_n(i+2,j,k)+d7(3)*rhou_n(i+3,j,k) &
                  + d7(1)*rhou_n(i-1,j,k)+d7(2)*rhou_n(i-2,j,k)+d7(3)*rhou_n(i-3,j,k)
             Krhov(i,j,k)= d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i+1,j,k)+d7(2)*rhov_n(i+2,j,k)+d7(3)*rhov_n(i+3,j,k) &
                  + d7(1)*rhov_n(i-1,j,k)+d7(2)*rhov_n(i-2,j,k)+d7(3)*rhov_n(i-3,j,k)
             Krhow(i,j,k)= d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i+1,j,k)+d7(2)*rhow_n(i+2,j,k)+d7(3)*rhow_n(i+3,j,k) &
                  + d7(1)*rhow_n(i-1,j,k)+d7(2)*rhow_n(i-2,j,k)+d7(3)*rhow_n(i-3,j,k)
             Krhoe(i,j,k)= d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i+1,j,k)+d7(2)*rhoe_n(i+2,j,k)+d7(3)*rhoe_n(i+3,j,k) &
                  + d7(1)*rhoe_n(i-1,j,k)+d7(2)*rhoe_n(i-2,j,k)+d7(3)*rhoe_n(i-3,j,k)
          enddo
       enddo
    else
       ! 8-th order centered filter
       do k=1,nz
          do j=1,ny
             do i=1,4
                Krho(i,j,k) = d9(0)*rho_n(i,j,k) &
                     + d9(1)*( rho_n(i+1,j,k)+rho_n(i-1,j,k) )&
                     + d9(2)*( rho_n(i+2,j,k)+rho_n(i-2,j,k) )&
                     + d9(3)*( rho_n(i+3,j,k)+rho_n(i-3,j,k) )&
                     + d9(4)*( rho_n(i+4,j,k)+rho_n(i-4,j,k) )

                Krhou(i,j,k)= d9(0)*rhou_n(i,j,k) &
                     + d9(1)*( rhou_n(i+1,j,k)+rhou_n(i-1,j,k) )&
                     + d9(2)*( rhou_n(i+2,j,k)+rhou_n(i-2,j,k) )&
                     + d9(3)*( rhou_n(i+3,j,k)+rhou_n(i-3,j,k) )&
                     + d9(4)*( rhou_n(i+4,j,k)+rhou_n(i-4,j,k) )

                Krhov(i,j,k)= d9(0)*rhov_n(i,j,k) &
                     + d9(1)*( rhov_n(i+1,j,k)+rhov_n(i-1,j,k) )&
                     + d9(2)*( rhov_n(i+2,j,k)+rhov_n(i-2,j,k) )&
                     + d9(3)*( rhov_n(i+3,j,k)+rhov_n(i-3,j,k) )&
                     + d9(4)*( rhov_n(i+4,j,k)+rhov_n(i-4,j,k) )

                Krhow(i,j,k)= d9(0)*rhow_n(i,j,k) &
                     + d9(1)*( rhow_n(i+1,j,k)+rhow_n(i-1,j,k) )&
                     + d9(2)*( rhow_n(i+2,j,k)+rhow_n(i-2,j,k) )&
                     + d9(3)*( rhow_n(i+3,j,k)+rhow_n(i-3,j,k) )&
                     + d9(4)*( rhow_n(i+4,j,k)+rhow_n(i-4,j,k) )

                Krhoe(i,j,k)= d9(0)*rhoe_n(i,j,k) &
                     + d9(1)*( rhoe_n(i+1,j,k)+rhoe_n(i-1,j,k) )&
                     + d9(2)*( rhoe_n(i+2,j,k)+rhoe_n(i-2,j,k) )&
                     + d9(3)*( rhoe_n(i+3,j,k)+rhoe_n(i-3,j,k) )&
                     + d9(4)*( rhoe_n(i+4,j,k)+rhoe_n(i-4,j,k) )
             enddo
          enddo
       enddo
    endif

    ! BC at imax
    ! ----------
    if (is_boundary(1,2)) then
       ! 6-th order centered filter
       i=nx-3
       do k=1,nz
          do j=1,ny
             Krho(i,j,k) = d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i+1,j,k)+d7(2)*rho_n(i+2,j,k)+d7(3)*rho_n(i+3,j,k) &
                  + d7(1)*rho_n(i-1,j,k)+d7(2)*rho_n(i-2,j,k)+d7(3)*rho_n(i-3,j,k)
             Krhou(i,j,k)= d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i+1,j,k)+d7(2)*rhou_n(i+2,j,k)+d7(3)*rhou_n(i+3,j,k) &
                  + d7(1)*rhou_n(i-1,j,k)+d7(2)*rhou_n(i-2,j,k)+d7(3)*rhou_n(i-3,j,k)
             Krhov(i,j,k)= d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i+1,j,k)+d7(2)*rhov_n(i+2,j,k)+d7(3)*rhov_n(i+3,j,k) &
                  + d7(1)*rhov_n(i-1,j,k)+d7(2)*rhov_n(i-2,j,k)+d7(3)*rhov_n(i-3,j,k)
             Krhow(i,j,k)= d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i+1,j,k)+d7(2)*rhow_n(i+2,j,k)+d7(3)*rhow_n(i+3,j,k) &
                  + d7(1)*rhow_n(i-1,j,k)+d7(2)*rhow_n(i-2,j,k)+d7(3)*rhow_n(i-3,j,k)
             Krhoe(i,j,k)= d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i+1,j,k)+d7(2)*rhoe_n(i+2,j,k)+d7(3)*rhoe_n(i+3,j,k) &
                  + d7(1)*rhoe_n(i-1,j,k)+d7(2)*rhoe_n(i-2,j,k)+d7(3)*rhoe_n(i-3,j,k)
          enddo
       enddo
       ! 4-th order centered filter
       i=nx-2
       do k=1,nz
          do j=1,ny
             Krho(i,j,k) = d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i+1,j,k)+d5(2)*rho_n(i+2,j,k) &
                  + d5(1)*rho_n(i-1,j,k)+d5(2)*rho_n(i-2,j,k)
             Krhou(i,j,k)= d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i+1,j,k)+d5(2)*rhou_n(i+2,j,k) &
                  + d5(1)*rhou_n(i-1,j,k)+d5(2)*rhou_n(i-2,j,k)
             Krhov(i,j,k)= d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i+1,j,k)+d5(2)*rhov_n(i+2,j,k) &
                  + d5(1)*rhov_n(i-1,j,k)+d5(2)*rhov_n(i-2,j,k)
             Krhow(i,j,k)= d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i+1,j,k)+d5(2)*rhow_n(i+2,j,k) &
                  + d5(1)*rhow_n(i-1,j,k)+d5(2)*rhow_n(i-2,j,k)
             Krhoe(i,j,k)= d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i+1,j,k)+d5(2)*rhoe_n(i+2,j,k) &
                  + d5(1)*rhoe_n(i-1,j,k)+d5(2)*rhoe_n(i-2,j,k)
          enddo
       enddo
       ! no filter for i=nx-1 & i=nx -> initialize array
       do j=1,ny
           Krho(nx-1,j,:)=0.0_wp
          Krhou(nx-1,j,:)=0.0_wp
          Krhov(nx-1,j,:)=0.0_wp
          Krhow(nx-1,j,:)=0.0_wp
          Krhoe(nx-1,j,:)=0.0_wp
          
           Krho(nx,j,:)=0.0_wp
          Krhou(nx,j,:)=0.0_wp
          Krhov(nx,j,:)=0.0_wp
          Krhow(nx,j,:)=0.0_wp
          Krhoe(nx,j,:)=0.0_wp
       enddo
    else
       ! 8-th order centered filter
       do k=1,nz
          do j=1,ny
             do i=nx-3,nx
                Krho(i,j,k) = d9(0)*rho_n(i,j,k) &
                     + d9(1)*( rho_n(i+1,j,k)+rho_n(i-1,j,k) )&
                     + d9(2)*( rho_n(i+2,j,k)+rho_n(i-2,j,k) )&
                     + d9(3)*( rho_n(i+3,j,k)+rho_n(i-3,j,k) )&
                     + d9(4)*( rho_n(i+4,j,k)+rho_n(i-4,j,k) )

                Krhou(i,j,k)= d9(0)*rhou_n(i,j,k) &
                     + d9(1)*( rhou_n(i+1,j,k)+rhou_n(i-1,j,k) )&
                     + d9(2)*( rhou_n(i+2,j,k)+rhou_n(i-2,j,k) )&
                     + d9(3)*( rhou_n(i+3,j,k)+rhou_n(i-3,j,k) )&
                     + d9(4)*( rhou_n(i+4,j,k)+rhou_n(i-4,j,k) )

                Krhov(i,j,k)= d9(0)*rhov_n(i,j,k) &
                     + d9(1)*( rhov_n(i+1,j,k)+rhov_n(i-1,j,k) )&
                     + d9(2)*( rhov_n(i+2,j,k)+rhov_n(i-2,j,k) )&
                     + d9(3)*( rhov_n(i+3,j,k)+rhov_n(i-3,j,k) )&
                     + d9(4)*( rhov_n(i+4,j,k)+rhov_n(i-4,j,k) )

                Krhow(i,j,k)= d9(0)*rhow_n(i,j,k) &
                     + d9(1)*( rhow_n(i+1,j,k)+rhow_n(i-1,j,k) )&
                     + d9(2)*( rhow_n(i+2,j,k)+rhow_n(i-2,j,k) )&
                     + d9(3)*( rhow_n(i+3,j,k)+rhow_n(i-3,j,k) )&
                     + d9(4)*( rhow_n(i+4,j,k)+rhow_n(i-4,j,k) )

                Krhoe(i,j,k)= d9(0)*rhoe_n(i,j,k) &
                     + d9(1)*( rhoe_n(i+1,j,k)+rhoe_n(i-1,j,k) )&
                     + d9(2)*( rhoe_n(i+2,j,k)+rhoe_n(i-2,j,k) )&
                     + d9(3)*( rhoe_n(i+3,j,k)+rhoe_n(i-3,j,k) )&
                     + d9(4)*( rhoe_n(i+4,j,k)+rhoe_n(i-4,j,k) )
             enddo
          enddo
       enddo
    endif
    
    ! Filtering along j-direction
    ! ===========================
    
    ! Perform one-sided communications in j-direction
    ! -----------------------------------------------
    ! open target windows
    call window_comm1
    ! communicate var
    call communication_1_j
    
    ! Interior points
    ! ---------------
    ! 8-th order centered filter
    do k=1,nz
       do j=5,ny-4
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  + d9(0)*rho_n(i,j,k) &
                  + d9(1)*( rho_n(i,j+1,k)+rho_n(i,j-1,k) )&
                  + d9(2)*( rho_n(i,j+2,k)+rho_n(i,j-2,k) )&
                  + d9(3)*( rho_n(i,j+3,k)+rho_n(i,j-3,k) )&
                  + d9(4)*( rho_n(i,j+4,k)+rho_n(i,j-4,k) )

             Krhou(i,j,k)= Krhou(i,j,k)  + d9(0)*rhou_n(i,j,k) &
                  + d9(1)*( rhou_n(i,j+1,k)+rhou_n(i,j-1,k) )&
                  + d9(2)*( rhou_n(i,j+2,k)+rhou_n(i,j-2,k) )&
                  + d9(3)*( rhou_n(i,j+3,k)+rhou_n(i,j-3,k) )&
                  + d9(4)*( rhou_n(i,j+4,k)+rhou_n(i,j-4,k) )

             Krhov(i,j,k)= Krhov(i,j,k)  + d9(0)*rhov_n(i,j,k) &
                  + d9(1)*( rhov_n(i,j+1,k)+rhov_n(i,j-1,k) )&
                  + d9(2)*( rhov_n(i,j+2,k)+rhov_n(i,j-2,k) )&
                  + d9(3)*( rhov_n(i,j+3,k)+rhov_n(i,j-3,k) )&
                  + d9(4)*( rhov_n(i,j+4,k)+rhov_n(i,j-4,k) )

             Krhow(i,j,k)= Krhow(i,j,k)  + d9(0)*rhow_n(i,j,k) &
                  + d9(1)*( rhow_n(i,j+1,k)+rhow_n(i,j-1,k) )&
                  + d9(2)*( rhow_n(i,j+2,k)+rhow_n(i,j-2,k) )&
                  + d9(3)*( rhow_n(i,j+3,k)+rhow_n(i,j-3,k) )&
                  + d9(4)*( rhow_n(i,j+4,k)+rhow_n(i,j-4,k) )

             Krhoe(i,j,k)= Krhoe(i,j,k)  + d9(0)*rhoe_n(i,j,k) &
                  + d9(1)*( rhoe_n(i,j+1,k)+rhoe_n(i,j-1,k) )&
                  + d9(2)*( rhoe_n(i,j+2,k)+rhoe_n(i,j-2,k) )&
                  + d9(3)*( rhoe_n(i,j+3,k)+rhoe_n(i,j-3,k) )&
                  + d9(4)*( rhoe_n(i,j+4,k)+rhoe_n(i,j-4,k) )
          enddo
       enddo
    enddo
    
    ! Check that communications are fullfilled
    ! ----------------------------------------
    call window_comm1

    ! BC at jmin
    ! ----------
    if (is_boundary(2,1)) then
       ! no filter for j=1 & j=2
       
       ! 4-th order centered filter
       j=3
       do k=1,nz
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i,j+1,k)+d5(2)*rho_n(i,j+2,k) &
                  + d5(1)*rho_n(i,j-1,k)+d5(2)*rho_n(i,j-2,k)
             Krhou(i,j,k)= Krhou(i,j,k)  +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j+1,k)+d5(2)*rhou_n(i,j+2,k) &
                  + d5(1)*rhou_n(i,j-1,k)+d5(2)*rhou_n(i,j-2,k)
             Krhov(i,j,k)= Krhov(i,j,k)  +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j+1,k)+d5(2)*rhov_n(i,j+2,k) &
                  + d5(1)*rhov_n(i,j-1,k)+d5(2)*rhov_n(i,j-2,k)
             Krhow(i,j,k)= Krhow(i,j,k)  +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j+1,k)+d5(2)*rhow_n(i,j+2,k) &
                  + d5(1)*rhow_n(i,j-1,k)+d5(2)*rhow_n(i,j-2,k)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i,j+1,k)+d5(2)*rhoe_n(i,j+2,k) &
                  + d5(1)*rhoe_n(i,j-1,k)+d5(2)*rhoe_n(i,j-2,k)
          enddo
       enddo
       ! 6-th order centered filter
       j=4
       do k=1,nz
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i,j+1,k)+d7(2)*rho_n(i,j+2,k)+d7(3)*rho_n(i,j+3,k) &
                  + d7(1)*rho_n(i,j-1,k)+d7(2)*rho_n(i,j-2,k)+d7(3)*rho_n(i,j-3,k)
             Krhou(i,j,k)= Krhou(i,j,k)  +d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i,j+1,k)+d7(2)*rhou_n(i,j+2,k)+d7(3)*rhou_n(i,j+3,k) &
                  + d7(1)*rhou_n(i,j-1,k)+d7(2)*rhou_n(i,j-2,k)+d7(3)*rhou_n(i,j-3,k)
             Krhov(i,j,k)= Krhov(i,j,k)  +d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i,j+1,k)+d7(2)*rhov_n(i,j+2,k)+d7(3)*rhov_n(i,j+3,k) &
                  + d7(1)*rhov_n(i,j-1,k)+d7(2)*rhov_n(i,j-2,k)+d7(3)*rhov_n(i,j-3,k)
             Krhow(i,j,k)= Krhow(i,j,k)  +d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i,j+1,k)+d7(2)*rhow_n(i,j+2,k)+d7(3)*rhow_n(i,j+3,k) &
                  + d7(1)*rhow_n(i,j-1,k)+d7(2)*rhow_n(i,j-2,k)+d7(3)*rhow_n(i,j-3,k)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i,j+1,k)+d7(2)*rhoe_n(i,j+2,k)+d7(3)*rhoe_n(i,j+3,k) &
                  + d7(1)*rhoe_n(i,j-1,k)+d7(2)*rhoe_n(i,j-2,k)+d7(3)*rhoe_n(i,j-3,k)
          enddo
       enddo
    else
       ! 8-th order centered filter
       do k=1,nz
          do j=1,4
             do i=1,nx
                Krho(i,j,k) = Krho(i,j,k)  + d9(0)*rho_n(i,j,k) &
                     + d9(1)*( rho_n(i,j+1,k)+rho_n(i,j-1,k) )&
                     + d9(2)*( rho_n(i,j+2,k)+rho_n(i,j-2,k) )&
                     + d9(3)*( rho_n(i,j+3,k)+rho_n(i,j-3,k) )&
                     + d9(4)*( rho_n(i,j+4,k)+rho_n(i,j-4,k) )

                Krhou(i,j,k)= Krhou(i,j,k)  + d9(0)*rhou_n(i,j,k) &
                     + d9(1)*( rhou_n(i,j+1,k)+rhou_n(i,j-1,k) )&
                     + d9(2)*( rhou_n(i,j+2,k)+rhou_n(i,j-2,k) )&
                     + d9(3)*( rhou_n(i,j+3,k)+rhou_n(i,j-3,k) )&
                     + d9(4)*( rhou_n(i,j+4,k)+rhou_n(i,j-4,k) )

                Krhov(i,j,k)= Krhov(i,j,k)  + d9(0)*rhov_n(i,j,k) &
                     + d9(1)*( rhov_n(i,j+1,k)+rhov_n(i,j-1,k) )&
                     + d9(2)*( rhov_n(i,j+2,k)+rhov_n(i,j-2,k) )&
                     + d9(3)*( rhov_n(i,j+3,k)+rhov_n(i,j-3,k) )&
                     + d9(4)*( rhov_n(i,j+4,k)+rhov_n(i,j-4,k) )

                Krhow(i,j,k)= Krhow(i,j,k)  + d9(0)*rhow_n(i,j,k) &
                     + d9(1)*( rhow_n(i,j+1,k)+rhow_n(i,j-1,k) )&
                     + d9(2)*( rhow_n(i,j+2,k)+rhow_n(i,j-2,k) )&
                     + d9(3)*( rhow_n(i,j+3,k)+rhow_n(i,j-3,k) )&
                     + d9(4)*( rhow_n(i,j+4,k)+rhow_n(i,j-4,k) )

                Krhoe(i,j,k)= Krhoe(i,j,k)  + d9(0)*rhoe_n(i,j,k) &
                     + d9(1)*( rhoe_n(i,j+1,k)+rhoe_n(i,j-1,k) )&
                     + d9(2)*( rhoe_n(i,j+2,k)+rhoe_n(i,j-2,k) )&
                     + d9(3)*( rhoe_n(i,j+3,k)+rhoe_n(i,j-3,k) )&
                     + d9(4)*( rhoe_n(i,j+4,k)+rhoe_n(i,j-4,k) )
             enddo
          enddo
       enddo
    endif

    ! BC at jmax
    ! ----------
    if (is_boundary(2,2)) then
       ! 6-th order centered filter
       j=ny-3
       do k=1,nz
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i,j+1,k)+d7(2)*rho_n(i,j+2,k)+d7(3)*rho_n(i,j+3,k) &
                  + d7(1)*rho_n(i,j-1,k)+d7(2)*rho_n(i,j-2,k)+d7(3)*rho_n(i,j-3,k)
             Krhou(i,j,k)= Krhou(i,j,k)  +d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i,j+1,k)+d7(2)*rhou_n(i,j+2,k)+d7(3)*rhou_n(i,j+3,k) &
                  + d7(1)*rhou_n(i,j-1,k)+d7(2)*rhou_n(i,j-2,k)+d7(3)*rhou_n(i,j-3,k)
             Krhov(i,j,k)= Krhov(i,j,k)  +d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i,j+1,k)+d7(2)*rhov_n(i,j+2,k)+d7(3)*rhov_n(i,j+3,k) &
                  + d7(1)*rhov_n(i,j-1,k)+d7(2)*rhov_n(i,j-2,k)+d7(3)*rhov_n(i,j-3,k)
             Krhow(i,j,k)= Krhow(i,j,k)  +d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i,j+1,k)+d7(2)*rhow_n(i,j+2,k)+d7(3)*rhow_n(i,j+3,k) &
                  + d7(1)*rhow_n(i,j-1,k)+d7(2)*rhow_n(i,j-2,k)+d7(3)*rhow_n(i,j-3,k)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i,j+1,k)+d7(2)*rhoe_n(i,j+2,k)+d7(3)*rhoe_n(i,j+3,k) &
                  + d7(1)*rhoe_n(i,j-1,k)+d7(2)*rhoe_n(i,j-2,k)+d7(3)*rhoe_n(i,j-3,k)
          enddo
       enddo
       ! 4-th order centered filter
       j=ny-2
       do k=1,nz
          do i=1,nx
             Krho(i,j,k)= Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i,j+1,k)+d5(2)*rho_n(i,j+2,k) &
                  + d5(1)*rho_n(i,j-1,k)+d5(2)*rho_n(i,j-2,k)
             Krhou(i,j,k)= Krhou(i,j,k)  +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j+1,k)+d5(2)*rhou_n(i,j+2,k) &
                  + d5(1)*rhou_n(i,j-1,k)+d5(2)*rhou_n(i,j-2,k)
             Krhov(i,j,k)= Krhov(i,j,k)  +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j+1,k)+d5(2)*rhov_n(i,j+2,k) &
                  + d5(1)*rhov_n(i,j-1,k)+d5(2)*rhov_n(i,j-2,k)
             Krhow(i,j,k)= Krhow(i,j,k)  +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j+1,k)+d5(2)*rhow_n(i,j+2,k) &
                  + d5(1)*rhow_n(i,j-1,k)+d5(2)*rhow_n(i,j-2,k)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i,j+1,k)+d5(2)*rhoe_n(i,j+2,k) &
                  + d5(1)*rhoe_n(i,j-1,k)+d5(2)*rhoe_n(i,j-2,k)
          enddo
       enddo
       ! no filter for j=ny-1 & j=ny
    else
       ! 8-th order centered filter
       do k=1,nz
          do j=ny-3,ny
             do i=1,nx
                Krho(i,j,k) = Krho(i,j,k)  + d9(0)*rho_n(i,j,k) &
                     + d9(1)*( rho_n(i,j+1,k)+rho_n(i,j-1,k) )&
                     + d9(2)*( rho_n(i,j+2,k)+rho_n(i,j-2,k) )&
                     + d9(3)*( rho_n(i,j+3,k)+rho_n(i,j-3,k) )&
                     + d9(4)*( rho_n(i,j+4,k)+rho_n(i,j-4,k) )

                Krhou(i,j,k)= Krhou(i,j,k)  + d9(0)*rhou_n(i,j,k) &
                     + d9(1)*( rhou_n(i,j+1,k)+rhou_n(i,j-1,k) )&
                     + d9(2)*( rhou_n(i,j+2,k)+rhou_n(i,j-2,k) )&
                     + d9(3)*( rhou_n(i,j+3,k)+rhou_n(i,j-3,k) )&
                     + d9(4)*( rhou_n(i,j+4,k)+rhou_n(i,j-4,k) )

                Krhov(i,j,k)= Krhov(i,j,k)  + d9(0)*rhov_n(i,j,k) &
                     + d9(1)*( rhov_n(i,j+1,k)+rhov_n(i,j-1,k) )&
                     + d9(2)*( rhov_n(i,j+2,k)+rhov_n(i,j-2,k) )&
                     + d9(3)*( rhov_n(i,j+3,k)+rhov_n(i,j-3,k) )&
                     + d9(4)*( rhov_n(i,j+4,k)+rhov_n(i,j-4,k) )

                Krhow(i,j,k)= Krhow(i,j,k)  + d9(0)*rhow_n(i,j,k) &
                     + d9(1)*( rhow_n(i,j+1,k)+rhow_n(i,j-1,k) )&
                     + d9(2)*( rhow_n(i,j+2,k)+rhow_n(i,j-2,k) )&
                     + d9(3)*( rhow_n(i,j+3,k)+rhow_n(i,j-3,k) )&
                     + d9(4)*( rhow_n(i,j+4,k)+rhow_n(i,j-4,k) )

                Krhoe(i,j,k)= Krhoe(i,j,k)  + d9(0)*rhoe_n(i,j,k) &
                     + d9(1)*( rhoe_n(i,j+1,k)+rhoe_n(i,j-1,k) )&
                     + d9(2)*( rhoe_n(i,j+2,k)+rhoe_n(i,j-2,k) )&
                     + d9(3)*( rhoe_n(i,j+3,k)+rhoe_n(i,j-3,k) )&
                     + d9(4)*( rhoe_n(i,j+4,k)+rhoe_n(i,j-4,k) )
             enddo
          enddo
       enddo
    endif

    !****************
    if (is_2D) return
    !****************
    
    ! Filtering along k-direction
    ! ===========================

    ! Perform one-sided communications in k-direction
    ! -----------------------------------------------
    ! open target windows
    call window_comm1
    ! communicate var
    call communication_1_k
    
    ! Interior points
    ! ---------------
    ! 8-th order centered filter
    do k=5,nz-4
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  + d9(0)*rho_n(i,j,k) &
                  + d9(1)*( rho_n(i,j,k+1)+rho_n(i,j,k-1) )&
                  + d9(2)*( rho_n(i,j,k+2)+rho_n(i,j,k-2) )&
                  + d9(3)*( rho_n(i,j,k+3)+rho_n(i,j,k-3) )&
                  + d9(4)*( rho_n(i,j,k+4)+rho_n(i,j,k-4) )

             Krhou(i,j,k)= Krhou(i,j,k)  + d9(0)*rhou_n(i,j,k) &
                  + d9(1)*( rhou_n(i,j,k+1)+rhou_n(i,j,k-1) )&
                  + d9(2)*( rhou_n(i,j,k+2)+rhou_n(i,j,k-2) )&
                  + d9(3)*( rhou_n(i,j,k+3)+rhou_n(i,j,k-3) )&
                  + d9(4)*( rhou_n(i,j,k+4)+rhou_n(i,j,k-4) )

             Krhov(i,j,k)= Krhov(i,j,k)  + d9(0)*rhov_n(i,j,k) &
                  + d9(1)*( rhov_n(i,j,k+1)+rhov_n(i,j,k-1) )&
                  + d9(2)*( rhov_n(i,j,k+2)+rhov_n(i,j,k-2) )&
                  + d9(3)*( rhov_n(i,j,k+3)+rhov_n(i,j,k-3) )&
                  + d9(4)*( rhov_n(i,j,k+4)+rhov_n(i,j,k-4) )

             Krhow(i,j,k)= Krhow(i,j,k)  + d9(0)*rhow_n(i,j,k) &
                  + d9(1)*( rhow_n(i,j,k+1)+rhow_n(i,j,k-1) )&
                  + d9(2)*( rhow_n(i,j,k+2)+rhow_n(i,j,k-2) )&
                  + d9(3)*( rhow_n(i,j,k+3)+rhow_n(i,j,k-3) )&
                  + d9(4)*( rhow_n(i,j,k+4)+rhow_n(i,j,k-4) )

             Krhoe(i,j,k)= Krhoe(i,j,k)  + d9(0)*rhoe_n(i,j,k) &
                  + d9(1)*( rhoe_n(i,j,k+1)+rhoe_n(i,j,k-1) )&
                  + d9(2)*( rhoe_n(i,j,k+2)+rhoe_n(i,j,k-2) )&
                  + d9(3)*( rhoe_n(i,j,k+3)+rhoe_n(i,j,k-3) )&
                  + d9(4)*( rhoe_n(i,j,k+4)+rhoe_n(i,j,k-4) )
          enddo
       enddo
    enddo
    
    ! Check that communications are fullfilled
    ! ----------------------------------------
    call window_comm1

    ! BC at kmin
    ! ----------
    if (is_boundary(3,1)) then
       ! no filter for k=1 & k=2

       ! 4-th order centered filter
       k=3
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i,j,k+1)+d5(2)*rho_n(i,j,k+2) &
                  + d5(1)*rho_n(i,j,k-1)+d5(2)*rho_n(i,j,k-2)
             Krhou(i,j,k)= Krhou(i,j,k)  +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j,k+1)+d5(2)*rhou_n(i,j,k+2) &
                  + d5(1)*rhou_n(i,j,k-1)+d5(2)*rhou_n(i,j,k-2)
             Krhov(i,j,k)= Krhov(i,j,k)  +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j,k+1)+d5(2)*rhov_n(i,j,k+2) &
                  + d5(1)*rhov_n(i,j,k-1)+d5(2)*rhov_n(i,j,k-2)
             Krhow(i,j,k)= Krhow(i,j,k)  +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j,k+1)+d5(2)*rhow_n(i,j,k+2) &
                  + d5(1)*rhow_n(i,j,k-1)+d5(2)*rhow_n(i,j,k-2)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i,j,k+1)+d5(2)*rhoe_n(i,j,k+2) &
                  + d5(1)*rhoe_n(i,j,k-1)+d5(2)*rhoe_n(i,j,k-2)
          enddo
       enddo
       ! 6-th order centered filter
       k=4
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i,j,k+1)+d7(2)*rho_n(i,j,k+2)+d7(3)*rho_n(i,j,k+3) &
                  + d7(1)*rho_n(i,j,k-1)+d7(2)*rho_n(i,j,k-2)+d7(3)*rho_n(i,j,k-3)
             Krhou(i,j,k)= Krhou(i,j,k)  +d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i,j,k+1)+d7(2)*rhou_n(i,j,k+2)+d7(3)*rhou_n(i,j,k+3) &
                  + d7(1)*rhou_n(i,j,k-1)+d7(2)*rhou_n(i,j,k-2)+d7(3)*rhou_n(i,j,k-3)
             Krhov(i,j,k)= Krhov(i,j,k)  +d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i,j,k+1)+d7(2)*rhov_n(i,j,k+2)+d7(3)*rhov_n(i,j,k+3) &
                  + d7(1)*rhov_n(i,j,k-1)+d7(2)*rhov_n(i,j,k-2)+d7(3)*rhov_n(i,j,k-3)
             Krhow(i,j,k)= Krhow(i,j,k)  +d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i,j,k+1)+d7(2)*rhow_n(i,j,k+2)+d7(3)*rhow_n(i,j,k+3) &
                  + d7(1)*rhow_n(i,j,k-1)+d7(2)*rhow_n(i,j,k-2)+d7(3)*rhow_n(i,j,k-3)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i,j,k+1)+d7(2)*rhoe_n(i,j,k+2)+d7(3)*rhoe_n(i,j,k+3) &
                  + d7(1)*rhoe_n(i,j,k-1)+d7(2)*rhoe_n(i,j,k-2)+d7(3)*rhoe_n(i,j,k-3)
          enddo
       enddo
    else
       ! 8-th order centered filter
       do k=1,4
          do j=1,ny
             do i=1,nx
                Krho(i,j,k) = Krho(i,j,k)  + d9(0)*rho_n(i,j,k) &
                     + d9(1)*( rho_n(i,j,k+1)+rho_n(i,j,k-1) )&
                     + d9(2)*( rho_n(i,j,k+2)+rho_n(i,j,k-2) )&
                     + d9(3)*( rho_n(i,j,k+3)+rho_n(i,j,k-3) )&
                     + d9(4)*( rho_n(i,j,k+4)+rho_n(i,j,k-4) )

                Krhou(i,j,k)= Krhou(i,j,k)  + d9(0)*rhou_n(i,j,k) &
                     + d9(1)*( rhou_n(i,j,k+1)+rhou_n(i,j,k-1) )&
                     + d9(2)*( rhou_n(i,j,k+2)+rhou_n(i,j,k-2) )&
                     + d9(3)*( rhou_n(i,j,k+3)+rhou_n(i,j,k-3) )&
                     + d9(4)*( rhou_n(i,j,k+4)+rhou_n(i,j,k-4) )

                Krhov(i,j,k)= Krhov(i,j,k)  + d9(0)*rhov_n(i,j,k) &
                     + d9(1)*( rhov_n(i,j,k+1)+rhov_n(i,j,k-1) )&
                     + d9(2)*( rhov_n(i,j,k+2)+rhov_n(i,j,k-2) )&
                     + d9(3)*( rhov_n(i,j,k+3)+rhov_n(i,j,k-3) )&
                     + d9(4)*( rhov_n(i,j,k+4)+rhov_n(i,j,k-4) )

                Krhow(i,j,k)= Krhow(i,j,k)  + d9(0)*rhow_n(i,j,k) &
                     + d9(1)*( rhow_n(i,j,k+1)+rhow_n(i,j,k-1) )&
                     + d9(2)*( rhow_n(i,j,k+2)+rhow_n(i,j,k-2) )&
                     + d9(3)*( rhow_n(i,j,k+3)+rhow_n(i,j,k-3) )&
                     + d9(4)*( rhow_n(i,j,k+4)+rhow_n(i,j,k-4) )

                Krhoe(i,j,k)= Krhoe(i,j,k)  + d9(0)*rhoe_n(i,j,k) &
                     + d9(1)*( rhoe_n(i,j,k+1)+rhoe_n(i,j,k-1) )&
                     + d9(2)*( rhoe_n(i,j,k+2)+rhoe_n(i,j,k-2) )&
                     + d9(3)*( rhoe_n(i,j,k+3)+rhoe_n(i,j,k-3) )&
                     + d9(4)*( rhoe_n(i,j,k+4)+rhoe_n(i,j,k-4) )
             enddo
          enddo
       enddo
    endif

    ! BC at kmax
    ! ----------
    if (is_boundary(3,2)) then
       ! 6-th order centered filter
       k=nz-3
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d7(0)*rho_n(i,j,k)   &
                  + d7(1)*rho_n(i,j,k+1)+d7(2)*rho_n(i,j,k+2)+d7(3)*rho_n(i,j,k+3) &
                  + d7(1)*rho_n(i,j,k-1)+d7(2)*rho_n(i,j,k-2)+d7(3)*rho_n(i,j,k-3)
             Krhou(i,j,k)= Krhou(i,j,k)  +d7(0)*rhou_n(i,j,k)   &
                  + d7(1)*rhou_n(i,j,k+1)+d7(2)*rhou_n(i,j,k+2)+d7(3)*rhou_n(i,j,k+3) &
                  + d7(1)*rhou_n(i,j,k-1)+d7(2)*rhou_n(i,j,k-2)+d7(3)*rhou_n(i,j,k-3)
             Krhov(i,j,k)= Krhov(i,j,k)  +d7(0)*rhov_n(i,j,k)   &
                  + d7(1)*rhov_n(i,j,k+1)+d7(2)*rhov_n(i,j,k+2)+d7(3)*rhov_n(i,j,k+3) &
                  + d7(1)*rhov_n(i,j,k-1)+d7(2)*rhov_n(i,j,k-2)+d7(3)*rhov_n(i,j,k-3)
             Krhow(i,j,k)= Krhow(i,j,k)  +d7(0)*rhow_n(i,j,k)   &
                  + d7(1)*rhow_n(i,j,k+1)+d7(2)*rhow_n(i,j,k+2)+d7(3)*rhow_n(i,j,k+3) &
                  + d7(1)*rhow_n(i,j,k-1)+d7(2)*rhow_n(i,j,k-2)+d7(3)*rhow_n(i,j,k-3)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d7(0)*rhoe_n(i,j,k)   &
                  + d7(1)*rhoe_n(i,j,k+1)+d7(2)*rhoe_n(i,j,k+2)+d7(3)*rhoe_n(i,j,k+3) &
                  + d7(1)*rhoe_n(i,j,k-1)+d7(2)*rhoe_n(i,j,k-2)+d7(3)*rhoe_n(i,j,k-3)
          enddo
       enddo
       ! 4-th order centered filter
       k=nz-2
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i,j,k+1)+d5(2)*rho_n(i,j,k+2) &
                  + d5(1)*rho_n(i,j,k-1)+d5(2)*rho_n(i,j,k-2)
             Krhou(i,j,k)= Krhou(i,j,k)  +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j,k+1)+d5(2)*rhou_n(i,j,k+2) &
                  + d5(1)*rhou_n(i,j,k-1)+d5(2)*rhou_n(i,j,k-2)
             Krhov(i,j,k)= Krhov(i,j,k)  +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j,k+1)+d5(2)*rhov_n(i,j,k+2) &
                  + d5(1)*rhov_n(i,j,k-1)+d5(2)*rhov_n(i,j,k-2)
             Krhow(i,j,k)= Krhow(i,j,k)  +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j,k+1)+d5(2)*rhow_n(i,j,k+2) &
                  + d5(1)*rhow_n(i,j,k-1)+d5(2)*rhow_n(i,j,k-2)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i,j,k+1)+d5(2)*rhoe_n(i,j,k+2) &
                  + d5(1)*rhoe_n(i,j,k-1)+d5(2)*rhoe_n(i,j,k-2)
          enddo
       enddo
       ! no filter for k=nz-1 & k=nz
    else
       ! 8-th order centered filter
       do k=nz-3,nz
          do j=1,ny
             do i=1,nx
                Krho(i,j,k) = Krho(i,j,k)  + d9(0)*rho_n(i,j,k) &
                     + d9(1)*( rho_n(i,j,k+1)+rho_n(i,j,k-1) )&
                     + d9(2)*( rho_n(i,j,k+2)+rho_n(i,j,k-2) )&
                     + d9(3)*( rho_n(i,j,k+3)+rho_n(i,j,k-3) )&
                     + d9(4)*( rho_n(i,j,k+4)+rho_n(i,j,k-4) )

                Krhou(i,j,k)= Krhou(i,j,k)  + d9(0)*rhou_n(i,j,k) &
                     + d9(1)*( rhou_n(i,j,k+1)+rhou_n(i,j,k-1) )&
                     + d9(2)*( rhou_n(i,j,k+2)+rhou_n(i,j,k-2) )&
                     + d9(3)*( rhou_n(i,j,k+3)+rhou_n(i,j,k-3) )&
                     + d9(4)*( rhou_n(i,j,k+4)+rhou_n(i,j,k-4) )

                Krhov(i,j,k)= Krhov(i,j,k)  + d9(0)*rhov_n(i,j,k) &
                     + d9(1)*( rhov_n(i,j,k+1)+rhov_n(i,j,k-1) )&
                     + d9(2)*( rhov_n(i,j,k+2)+rhov_n(i,j,k-2) )&
                     + d9(3)*( rhov_n(i,j,k+3)+rhov_n(i,j,k-3) )&
                     + d9(4)*( rhov_n(i,j,k+4)+rhov_n(i,j,k-4) )

                Krhow(i,j,k)= Krhow(i,j,k)  + d9(0)*rhow_n(i,j,k) &
                     + d9(1)*( rhow_n(i,j,k+1)+rhow_n(i,j,k-1) )&
                     + d9(2)*( rhow_n(i,j,k+2)+rhow_n(i,j,k-2) )&
                     + d9(3)*( rhow_n(i,j,k+3)+rhow_n(i,j,k-3) )&
                     + d9(4)*( rhow_n(i,j,k+4)+rhow_n(i,j,k-4) )

                Krhoe(i,j,k)= Krhoe(i,j,k)  + d9(0)*rhoe_n(i,j,k) &
                     + d9(1)*( rhoe_n(i,j,k+1)+rhoe_n(i,j,k-1) )&
                     + d9(2)*( rhoe_n(i,j,k+2)+rhoe_n(i,j,k-2) )&
                     + d9(3)*( rhoe_n(i,j,k+3)+rhoe_n(i,j,k-3) )&
                     + d9(4)*( rhoe_n(i,j,k+4)+rhoe_n(i,j,k-4) )
             enddo
          enddo
       enddo
    endif
        
  end subroutine filtering_9pts_overlap
  
  !===============================================================================
  subroutine filtering_7pts_overlap
  !===============================================================================
    !> Apply filtering on 7-point stencils (+ boundaries)
    !> with computation-communication overlap
  !===============================================================================
    use mod_flow
    use mod_filtering
    use mod_comm1
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    ! ---------------------------------------------------------------------------

    ! Filtering along i-direction
    ! ===========================

    ! Perform one-sided communications in i-direction
    ! -----------------------------------------------
    ! open target windows
    call window_comm1
    ! communicate var
    call communication_1_i
    
    ! Interior points
    ! ---------------
    ! 6-th order centered filter
    do k=1,nz
       do j=1,ny
          do i=4,nx-3
             Krho(i,j,k) = d7(0)*rho_n(i,j,k) &
                  + d7(1)*( rho_n(i+1,j,k)+rho_n(i-1,j,k) )&
                  + d7(2)*( rho_n(i+2,j,k)+rho_n(i-2,j,k) )&
                  + d7(3)*( rho_n(i+3,j,k)+rho_n(i-3,j,k) )

             Krhou(i,j,k)= d7(0)*rhou_n(i,j,k) &
                  + d7(1)*( rhou_n(i+1,j,k)+rhou_n(i-1,j,k) )&
                  + d7(2)*( rhou_n(i+2,j,k)+rhou_n(i-2,j,k) )&
                  + d7(3)*( rhou_n(i+3,j,k)+rhou_n(i-3,j,k) )

             Krhov(i,j,k)= d7(0)*rhov_n(i,j,k) &
                  + d7(1)*( rhov_n(i+1,j,k)+rhov_n(i-1,j,k) )&
                  + d7(2)*( rhov_n(i+2,j,k)+rhov_n(i-2,j,k) )&
                  + d7(3)*( rhov_n(i+3,j,k)+rhov_n(i-3,j,k) )

             Krhow(i,j,k)= d7(0)*rhow_n(i,j,k) &
                  + d7(1)*( rhow_n(i+1,j,k)+rhow_n(i-1,j,k) )&
                  + d7(2)*( rhow_n(i+2,j,k)+rhow_n(i-2,j,k) )&
                  + d7(3)*( rhow_n(i+3,j,k)+rhow_n(i-3,j,k) )

             Krhoe(i,j,k)= d7(0)*rhoe_n(i,j,k) &
                  + d7(1)*( rhoe_n(i+1,j,k)+rhoe_n(i-1,j,k) )&
                  + d7(2)*( rhoe_n(i+2,j,k)+rhoe_n(i-2,j,k) )&
                  + d7(3)*( rhoe_n(i+3,j,k)+rhoe_n(i-3,j,k) )
          enddo
       enddo
    enddo
    
    ! Check that communications are fullfilled
    ! ----------------------------------------
    call window_comm1
    
    ! BC at imin
    ! ----------
    if (is_boundary(1,1)) then
       ! no filter for i=1 & i=2 -> initialize array
       do j=1,ny
           Krho(1,j,:)=0.0_wp
          Krhou(1,j,:)=0.0_wp
          Krhov(1,j,:)=0.0_wp
          Krhow(1,j,:)=0.0_wp
          Krhoe(1,j,:)=0.0_wp
          
           Krho(2,j,:)=0.0_wp
          Krhou(2,j,:)=0.0_wp
          Krhov(2,j,:)=0.0_wp
          Krhow(2,j,:)=0.0_wp
          Krhoe(2,j,:)=0.0_wp
       enddo
       ! 4-th order centered filter
       i=3
       do k=1,nz
          do j=1,ny
             Krho(i,j,k) = d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i+1,j,k)+d5(2)*rho_n(i+2,j,k) &
                  + d5(1)*rho_n(i-1,j,k)+d5(2)*rho_n(i-2,j,k)
             Krhou(i,j,k)= d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i+1,j,k)+d5(2)*rhou_n(i+2,j,k) &
                  + d5(1)*rhou_n(i-1,j,k)+d5(2)*rhou_n(i-2,j,k)
             Krhov(i,j,k)= d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i+1,j,k)+d5(2)*rhov_n(i+2,j,k) &
                  + d5(1)*rhov_n(i-1,j,k)+d5(2)*rhov_n(i-2,j,k)
             Krhow(i,j,k)= d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i+1,j,k)+d5(2)*rhow_n(i+2,j,k) &
                  + d5(1)*rhow_n(i-1,j,k)+d5(2)*rhow_n(i-2,j,k)
             Krhoe(i,j,k)= d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i+1,j,k)+d5(2)*rhoe_n(i+2,j,k) &
                  + d5(1)*rhoe_n(i-1,j,k)+d5(2)*rhoe_n(i-2,j,k)
          enddo
       enddo
    else
       ! 6-th order centered filter
       do k=1,nz
          do j=1,ny
             do i=1,3
                Krho(i,j,k) = d7(0)*rho_n(i,j,k) &
                     + d7(1)*( rho_n(i+1,j,k)+rho_n(i-1,j,k) )&
                     + d7(2)*( rho_n(i+2,j,k)+rho_n(i-2,j,k) )&
                     + d7(3)*( rho_n(i+3,j,k)+rho_n(i-3,j,k) )

                Krhou(i,j,k)= d7(0)*rhou_n(i,j,k) &
                     + d7(1)*( rhou_n(i+1,j,k)+rhou_n(i-1,j,k) )&
                     + d7(2)*( rhou_n(i+2,j,k)+rhou_n(i-2,j,k) )&
                     + d7(3)*( rhou_n(i+3,j,k)+rhou_n(i-3,j,k) )

                Krhov(i,j,k)= d7(0)*rhov_n(i,j,k) &
                     + d7(1)*( rhov_n(i+1,j,k)+rhov_n(i-1,j,k) )&
                     + d7(2)*( rhov_n(i+2,j,k)+rhov_n(i-2,j,k) )&
                     + d7(3)*( rhov_n(i+3,j,k)+rhov_n(i-3,j,k) )

                Krhow(i,j,k)= d7(0)*rhow_n(i,j,k) &
                     + d7(1)*( rhow_n(i+1,j,k)+rhow_n(i-1,j,k) )&
                     + d7(2)*( rhow_n(i+2,j,k)+rhow_n(i-2,j,k) )&
                     + d7(3)*( rhow_n(i+3,j,k)+rhow_n(i-3,j,k) )

                Krhoe(i,j,k)= d7(0)*rhoe_n(i,j,k) &
                     + d7(1)*( rhoe_n(i+1,j,k)+rhoe_n(i-1,j,k) )&
                     + d7(2)*( rhoe_n(i+2,j,k)+rhoe_n(i-2,j,k) )&
                     + d7(3)*( rhoe_n(i+3,j,k)+rhoe_n(i-3,j,k) )
             enddo
          enddo
       enddo
    endif

    ! BC at imax
    ! ----------
    if (is_boundary(1,2)) then
       ! 4-th order centered filter
       i=nx-2
       do k=1,nz
          do j=1,ny
             Krho(i,j,k) = d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i+1,j,k)+d5(2)*rho_n(i+2,j,k) &
                  + d5(1)*rho_n(i-1,j,k)+d5(2)*rho_n(i-2,j,k)
             Krhou(i,j,k)= d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i+1,j,k)+d5(2)*rhou_n(i+2,j,k) &
                  + d5(1)*rhou_n(i-1,j,k)+d5(2)*rhou_n(i-2,j,k)
             Krhov(i,j,k)= d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i+1,j,k)+d5(2)*rhov_n(i+2,j,k) &
                  + d5(1)*rhov_n(i-1,j,k)+d5(2)*rhov_n(i-2,j,k)
             Krhow(i,j,k)= d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i+1,j,k)+d5(2)*rhow_n(i+2,j,k) &
                  + d5(1)*rhow_n(i-1,j,k)+d5(2)*rhow_n(i-2,j,k)
             Krhoe(i,j,k)= d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i+1,j,k)+d5(2)*rhoe_n(i+2,j,k) &
                  + d5(1)*rhoe_n(i-1,j,k)+d5(2)*rhoe_n(i-2,j,k)
          enddo
       enddo
       ! no filter for i=nx-1 & i=nx -> initialize array
       do j=1,ny
           Krho(nx-1,j,:)=0.0_wp
          Krhou(nx-1,j,:)=0.0_wp
          Krhov(nx-1,j,:)=0.0_wp
          Krhow(nx-1,j,:)=0.0_wp
          Krhoe(nx-1,j,:)=0.0_wp
          
           Krho(nx,j,:)=0.0_wp
          Krhou(nx,j,:)=0.0_wp
          Krhov(nx,j,:)=0.0_wp
          Krhow(nx,j,:)=0.0_wp
          Krhoe(nx,j,:)=0.0_wp
       enddo
    else
       ! 6-th order centered filter
       do k=1,nz
          do j=1,ny
             do i=nx-2,nx
                Krho(i,j,k) = d7(0)*rho_n(i,j,k) &
                     + d7(1)*( rho_n(i+1,j,k)+rho_n(i-1,j,k) )&
                     + d7(2)*( rho_n(i+2,j,k)+rho_n(i-2,j,k) )&
                     + d7(3)*( rho_n(i+3,j,k)+rho_n(i-3,j,k) )

                Krhou(i,j,k)= d7(0)*rhou_n(i,j,k) &
                     + d7(1)*( rhou_n(i+1,j,k)+rhou_n(i-1,j,k) )&
                     + d7(2)*( rhou_n(i+2,j,k)+rhou_n(i-2,j,k) )&
                     + d7(3)*( rhou_n(i+3,j,k)+rhou_n(i-3,j,k) )

                Krhov(i,j,k)= d7(0)*rhov_n(i,j,k) &
                     + d7(1)*( rhov_n(i+1,j,k)+rhov_n(i-1,j,k) )&
                     + d7(2)*( rhov_n(i+2,j,k)+rhov_n(i-2,j,k) )&
                     + d7(3)*( rhov_n(i+3,j,k)+rhov_n(i-3,j,k) )

                Krhow(i,j,k)= d7(0)*rhow_n(i,j,k) &
                     + d7(1)*( rhow_n(i+1,j,k)+rhow_n(i-1,j,k) )&
                     + d7(2)*( rhow_n(i+2,j,k)+rhow_n(i-2,j,k) )&
                     + d7(3)*( rhow_n(i+3,j,k)+rhow_n(i-3,j,k) )

                Krhoe(i,j,k)= d7(0)*rhoe_n(i,j,k) &
                     + d7(1)*( rhoe_n(i+1,j,k)+rhoe_n(i-1,j,k) )&
                     + d7(2)*( rhoe_n(i+2,j,k)+rhoe_n(i-2,j,k) )&
                     + d7(3)*( rhoe_n(i+3,j,k)+rhoe_n(i-3,j,k) )
             enddo
          enddo
       enddo
    endif
    
    ! Filtering along j-direction
    ! ===========================
    
    ! Perform one-sided communications in j-direction
    ! -----------------------------------------------
    ! open target windows
    call window_comm1
    ! communicate var
    call communication_1_j
    
    ! Interior points
    ! ---------------
    ! 6-th order centered filter
    do k=1,nz
       do j=4,ny-3
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  + d7(0)*rho_n(i,j,k) &
                  + d7(1)*( rho_n(i,j+1,k)+rho_n(i,j-1,k) )&
                  + d7(2)*( rho_n(i,j+2,k)+rho_n(i,j-2,k) )&
                  + d7(3)*( rho_n(i,j+3,k)+rho_n(i,j-3,k) )

             Krhou(i,j,k)= Krhou(i,j,k)  + d7(0)*rhou_n(i,j,k) &
                  + d7(1)*( rhou_n(i,j+1,k)+rhou_n(i,j-1,k) )&
                  + d7(2)*( rhou_n(i,j+2,k)+rhou_n(i,j-2,k) )&
                  + d7(3)*( rhou_n(i,j+3,k)+rhou_n(i,j-3,k) )

             Krhov(i,j,k)= Krhov(i,j,k)  + d7(0)*rhov_n(i,j,k) &
                  + d7(1)*( rhov_n(i,j+1,k)+rhov_n(i,j-1,k) )&
                  + d7(2)*( rhov_n(i,j+2,k)+rhov_n(i,j-2,k) )&
                  + d7(3)*( rhov_n(i,j+3,k)+rhov_n(i,j-3,k) )

             Krhow(i,j,k)= Krhow(i,j,k)  + d7(0)*rhow_n(i,j,k) &
                  + d7(1)*( rhow_n(i,j+1,k)+rhow_n(i,j-1,k) )&
                  + d7(2)*( rhow_n(i,j+2,k)+rhow_n(i,j-2,k) )&
                  + d7(3)*( rhow_n(i,j+3,k)+rhow_n(i,j-3,k) )

             Krhoe(i,j,k)= Krhoe(i,j,k)  + d7(0)*rhoe_n(i,j,k) &
                  + d7(1)*( rhoe_n(i,j+1,k)+rhoe_n(i,j-1,k) )&
                  + d7(2)*( rhoe_n(i,j+2,k)+rhoe_n(i,j-2,k) )&
                  + d7(3)*( rhoe_n(i,j+3,k)+rhoe_n(i,j-3,k) )
          enddo
       enddo
    enddo
    
    ! Check that communications are fullfilled
    ! ----------------------------------------
    call window_comm1

    ! BC at jmin
    ! ----------
    if (is_boundary(2,1)) then
       ! no filter for j=1 & j=2
       
       ! 4-th order centered filter
       j=3
       do k=1,nz
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i,j+1,k)+d5(2)*rho_n(i,j+2,k) &
                  + d5(1)*rho_n(i,j-1,k)+d5(2)*rho_n(i,j-2,k)
             Krhou(i,j,k)= Krhou(i,j,k)  +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j+1,k)+d5(2)*rhou_n(i,j+2,k) &
                  + d5(1)*rhou_n(i,j-1,k)+d5(2)*rhou_n(i,j-2,k)
             Krhov(i,j,k)= Krhov(i,j,k)  +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j+1,k)+d5(2)*rhov_n(i,j+2,k) &
                  + d5(1)*rhov_n(i,j-1,k)+d5(2)*rhov_n(i,j-2,k)
             Krhow(i,j,k)= Krhow(i,j,k)  +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j+1,k)+d5(2)*rhow_n(i,j+2,k) &
                  + d5(1)*rhow_n(i,j-1,k)+d5(2)*rhow_n(i,j-2,k)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i,j+1,k)+d5(2)*rhoe_n(i,j+2,k) &
                  + d5(1)*rhoe_n(i,j-1,k)+d5(2)*rhoe_n(i,j-2,k)
          enddo
       enddo
    else
       ! 6-th order centered filter
       do k=1,nz
          do j=1,3
             do i=1,nx
                Krho(i,j,k) = Krho(i,j,k)  + d7(0)*rho_n(i,j,k) &
                     + d7(1)*( rho_n(i,j+1,k)+rho_n(i,j-1,k) )&
                     + d7(2)*( rho_n(i,j+2,k)+rho_n(i,j-2,k) )&
                     + d7(3)*( rho_n(i,j+3,k)+rho_n(i,j-3,k) )

                Krhou(i,j,k)= Krhou(i,j,k)  + d7(0)*rhou_n(i,j,k) &
                     + d7(1)*( rhou_n(i,j+1,k)+rhou_n(i,j-1,k) )&
                     + d7(2)*( rhou_n(i,j+2,k)+rhou_n(i,j-2,k) )&
                     + d7(3)*( rhou_n(i,j+3,k)+rhou_n(i,j-3,k) )

                Krhov(i,j,k)= Krhov(i,j,k)  + d7(0)*rhov_n(i,j,k) &
                     + d7(1)*( rhov_n(i,j+1,k)+rhov_n(i,j-1,k) )&
                     + d7(2)*( rhov_n(i,j+2,k)+rhov_n(i,j-2,k) )&
                     + d7(3)*( rhov_n(i,j+3,k)+rhov_n(i,j-3,k) )

                Krhow(i,j,k)= Krhow(i,j,k)  + d7(0)*rhow_n(i,j,k) &
                     + d7(1)*( rhow_n(i,j+1,k)+rhow_n(i,j-1,k) )&
                     + d7(2)*( rhow_n(i,j+2,k)+rhow_n(i,j-2,k) )&
                     + d7(3)*( rhow_n(i,j+3,k)+rhow_n(i,j-3,k) )

                Krhoe(i,j,k)= Krhoe(i,j,k)  + d7(0)*rhoe_n(i,j,k) &
                     + d7(1)*( rhoe_n(i,j+1,k)+rhoe_n(i,j-1,k) )&
                     + d7(2)*( rhoe_n(i,j+2,k)+rhoe_n(i,j-2,k) )&
                     + d7(3)*( rhoe_n(i,j+3,k)+rhoe_n(i,j-3,k) )
             enddo
          enddo
       enddo
    endif

    ! BC at jmax
    ! ----------
    if (is_boundary(2,2)) then
       ! 4-th order centered filter
       j=ny-2
       do k=1,nz
          do i=1,nx
             Krho(i,j,k)= Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i,j+1,k)+d5(2)*rho_n(i,j+2,k) &
                  + d5(1)*rho_n(i,j-1,k)+d5(2)*rho_n(i,j-2,k)
             Krhou(i,j,k)= Krhou(i,j,k)  +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j+1,k)+d5(2)*rhou_n(i,j+2,k) &
                  + d5(1)*rhou_n(i,j-1,k)+d5(2)*rhou_n(i,j-2,k)
             Krhov(i,j,k)= Krhov(i,j,k)  +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j+1,k)+d5(2)*rhov_n(i,j+2,k) &
                  + d5(1)*rhov_n(i,j-1,k)+d5(2)*rhov_n(i,j-2,k)
             Krhow(i,j,k)= Krhow(i,j,k)  +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j+1,k)+d5(2)*rhow_n(i,j+2,k) &
                  + d5(1)*rhow_n(i,j-1,k)+d5(2)*rhow_n(i,j-2,k)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i,j+1,k)+d5(2)*rhoe_n(i,j+2,k) &
                  + d5(1)*rhoe_n(i,j-1,k)+d5(2)*rhoe_n(i,j-2,k)
          enddo
       enddo
       ! no filter for j=ny-1 & j=ny
    else
       ! 6-th order centered filter
       do k=1,nz
          do j=ny-2,ny
             do i=1,nx
                Krho(i,j,k) = Krho(i,j,k)  + d7(0)*rho_n(i,j,k) &
                     + d7(1)*( rho_n(i,j+1,k)+rho_n(i,j-1,k) )&
                     + d7(2)*( rho_n(i,j+2,k)+rho_n(i,j-2,k) )&
                     + d7(3)*( rho_n(i,j+3,k)+rho_n(i,j-3,k) )

                Krhou(i,j,k)= Krhou(i,j,k)  + d7(0)*rhou_n(i,j,k) &
                     + d7(1)*( rhou_n(i,j+1,k)+rhou_n(i,j-1,k) )&
                     + d7(2)*( rhou_n(i,j+2,k)+rhou_n(i,j-2,k) )&
                     + d7(3)*( rhou_n(i,j+3,k)+rhou_n(i,j-3,k) )

                Krhov(i,j,k)= Krhov(i,j,k)  + d7(0)*rhov_n(i,j,k) &
                     + d7(1)*( rhov_n(i,j+1,k)+rhov_n(i,j-1,k) )&
                     + d7(2)*( rhov_n(i,j+2,k)+rhov_n(i,j-2,k) )&
                     + d7(3)*( rhov_n(i,j+3,k)+rhov_n(i,j-3,k) )

                Krhow(i,j,k)= Krhow(i,j,k)  + d7(0)*rhow_n(i,j,k) &
                     + d7(1)*( rhow_n(i,j+1,k)+rhow_n(i,j-1,k) )&
                     + d7(2)*( rhow_n(i,j+2,k)+rhow_n(i,j-2,k) )&
                     + d7(3)*( rhow_n(i,j+3,k)+rhow_n(i,j-3,k) )

                Krhoe(i,j,k)= Krhoe(i,j,k)  + d7(0)*rhoe_n(i,j,k) &
                     + d7(1)*( rhoe_n(i,j+1,k)+rhoe_n(i,j-1,k) )&
                     + d7(2)*( rhoe_n(i,j+2,k)+rhoe_n(i,j-2,k) )&
                     + d7(3)*( rhoe_n(i,j+3,k)+rhoe_n(i,j-3,k) )
             enddo
          enddo
       enddo
    endif

    !****************
    if (is_2D) return
    !****************
    
    ! Filtering along k-direction
    ! ===========================

    ! Perform one-sided communications in k-direction
    ! -----------------------------------------------
    ! open target windows
    call window_comm1
    ! communicate var
    call communication_1_k
    
    ! Interior points
    ! ---------------
    ! 6-th order centered filter
    do k=4,nz-3
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  + d7(0)*rho_n(i,j,k) &
                  + d7(1)*( rho_n(i,j,k+1)+rho_n(i,j,k-1) )&
                  + d7(2)*( rho_n(i,j,k+2)+rho_n(i,j,k-2) )&
                  + d7(3)*( rho_n(i,j,k+3)+rho_n(i,j,k-3) )

             Krhou(i,j,k)= Krhou(i,j,k)  + d7(0)*rhou_n(i,j,k) &
                  + d7(1)*( rhou_n(i,j,k+1)+rhou_n(i,j,k-1) )&
                  + d7(2)*( rhou_n(i,j,k+2)+rhou_n(i,j,k-2) )&
                  + d7(3)*( rhou_n(i,j,k+3)+rhou_n(i,j,k-3) )

             Krhov(i,j,k)= Krhov(i,j,k)  + d7(0)*rhov_n(i,j,k) &
                  + d7(1)*( rhov_n(i,j,k+1)+rhov_n(i,j,k-1) )&
                  + d7(2)*( rhov_n(i,j,k+2)+rhov_n(i,j,k-2) )&
                  + d7(3)*( rhov_n(i,j,k+3)+rhov_n(i,j,k-3) )

             Krhow(i,j,k)= Krhow(i,j,k)  + d7(0)*rhow_n(i,j,k) &
                  + d7(1)*( rhow_n(i,j,k+1)+rhow_n(i,j,k-1) )&
                  + d7(2)*( rhow_n(i,j,k+2)+rhow_n(i,j,k-2) )&
                  + d7(3)*( rhow_n(i,j,k+3)+rhow_n(i,j,k-3) )

             Krhoe(i,j,k)= Krhoe(i,j,k)  + d7(0)*rhoe_n(i,j,k) &
                  + d7(1)*( rhoe_n(i,j,k+1)+rhoe_n(i,j,k-1) )&
                  + d7(2)*( rhoe_n(i,j,k+2)+rhoe_n(i,j,k-2) )&
                  + d7(3)*( rhoe_n(i,j,k+3)+rhoe_n(i,j,k-3) )
          enddo
       enddo
    enddo
    
    ! Check that communications are fullfilled
    ! ----------------------------------------
    call window_comm1

    ! BC at kmin
    ! ----------
    if (is_boundary(3,1)) then
       ! no filter for k=1 & k=2

       ! 4-th order centered filter
       k=3
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i,j,k+1)+d5(2)*rho_n(i,j,k+2) &
                  + d5(1)*rho_n(i,j,k-1)+d5(2)*rho_n(i,j,k-2)
             Krhou(i,j,k)= Krhou(i,j,k)  +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j,k+1)+d5(2)*rhou_n(i,j,k+2) &
                  + d5(1)*rhou_n(i,j,k-1)+d5(2)*rhou_n(i,j,k-2)
             Krhov(i,j,k)= Krhov(i,j,k)  +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j,k+1)+d5(2)*rhov_n(i,j,k+2) &
                  + d5(1)*rhov_n(i,j,k-1)+d5(2)*rhov_n(i,j,k-2)
             Krhow(i,j,k)= Krhow(i,j,k)  +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j,k+1)+d5(2)*rhow_n(i,j,k+2) &
                  + d5(1)*rhow_n(i,j,k-1)+d5(2)*rhow_n(i,j,k-2)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i,j,k+1)+d5(2)*rhoe_n(i,j,k+2) &
                  + d5(1)*rhoe_n(i,j,k-1)+d5(2)*rhoe_n(i,j,k-2)
          enddo
       enddo
    else
       ! 6-th order centered filter
       do k=1,3
          do j=1,ny
             do i=1,nx
                Krho(i,j,k) = Krho(i,j,k)  + d7(0)*rho_n(i,j,k) &
                     + d7(1)*( rho_n(i,j,k+1)+rho_n(i,j,k-1) )&
                     + d7(2)*( rho_n(i,j,k+2)+rho_n(i,j,k-2) )&
                     + d7(3)*( rho_n(i,j,k+3)+rho_n(i,j,k-3) )

                Krhou(i,j,k)= Krhou(i,j,k)  + d7(0)*rhou_n(i,j,k) &
                     + d7(1)*( rhou_n(i,j,k+1)+rhou_n(i,j,k-1) )&
                     + d7(2)*( rhou_n(i,j,k+2)+rhou_n(i,j,k-2) )&
                     + d7(3)*( rhou_n(i,j,k+3)+rhou_n(i,j,k-3) )

                Krhov(i,j,k)= Krhov(i,j,k)  + d7(0)*rhov_n(i,j,k) &
                     + d7(1)*( rhov_n(i,j,k+1)+rhov_n(i,j,k-1) )&
                     + d7(2)*( rhov_n(i,j,k+2)+rhov_n(i,j,k-2) )&
                     + d7(3)*( rhov_n(i,j,k+3)+rhov_n(i,j,k-3) )

                Krhow(i,j,k)= Krhow(i,j,k)  + d7(0)*rhow_n(i,j,k) &
                     + d7(1)*( rhow_n(i,j,k+1)+rhow_n(i,j,k-1) )&
                     + d7(2)*( rhow_n(i,j,k+2)+rhow_n(i,j,k-2) )&
                     + d7(3)*( rhow_n(i,j,k+3)+rhow_n(i,j,k-3) )

                Krhoe(i,j,k)= Krhoe(i,j,k)  + d7(0)*rhoe_n(i,j,k) &
                     + d7(1)*( rhoe_n(i,j,k+1)+rhoe_n(i,j,k-1) )&
                     + d7(2)*( rhoe_n(i,j,k+2)+rhoe_n(i,j,k-2) )&
                     + d7(3)*( rhoe_n(i,j,k+3)+rhoe_n(i,j,k-3) )
             enddo
          enddo
       enddo
    endif

    ! BC at kmax
    ! ----------
    if (is_boundary(3,2)) then
       ! 4-th order centered filter
       k=nz-2
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  +d5(0)*rho_n(i,j,k)   &
                  + d5(1)*rho_n(i,j,k+1)+d5(2)*rho_n(i,j,k+2) &
                  + d5(1)*rho_n(i,j,k-1)+d5(2)*rho_n(i,j,k-2)
             Krhou(i,j,k)= Krhou(i,j,k)  +d5(0)*rhou_n(i,j,k)   &
                  + d5(1)*rhou_n(i,j,k+1)+d5(2)*rhou_n(i,j,k+2) &
                  + d5(1)*rhou_n(i,j,k-1)+d5(2)*rhou_n(i,j,k-2)
             Krhov(i,j,k)= Krhov(i,j,k)  +d5(0)*rhov_n(i,j,k)   &
                  + d5(1)*rhov_n(i,j,k+1)+d5(2)*rhov_n(i,j,k+2) &
                  + d5(1)*rhov_n(i,j,k-1)+d5(2)*rhov_n(i,j,k-2)
             Krhow(i,j,k)= Krhow(i,j,k)  +d5(0)*rhow_n(i,j,k)   &
                  + d5(1)*rhow_n(i,j,k+1)+d5(2)*rhow_n(i,j,k+2) &
                  + d5(1)*rhow_n(i,j,k-1)+d5(2)*rhow_n(i,j,k-2)
             Krhoe(i,j,k)= Krhoe(i,j,k)  +d5(0)*rhoe_n(i,j,k)   &
                  + d5(1)*rhoe_n(i,j,k+1)+d5(2)*rhoe_n(i,j,k+2) &
                  + d5(1)*rhoe_n(i,j,k-1)+d5(2)*rhoe_n(i,j,k-2)
          enddo
       enddo
       ! no filter for k=nz-1 & k=nz
    else
       ! 6-th order centered filter
       do k=nz-2,nz
          do j=1,ny
             do i=1,nx
                Krho(i,j,k) = Krho(i,j,k)  + d7(0)*rho_n(i,j,k) &
                     + d7(1)*( rho_n(i,j,k+1)+rho_n(i,j,k-1) )&
                     + d7(2)*( rho_n(i,j,k+2)+rho_n(i,j,k-2) )&
                     + d7(3)*( rho_n(i,j,k+3)+rho_n(i,j,k-3) )

                Krhou(i,j,k)= Krhou(i,j,k)  + d7(0)*rhou_n(i,j,k) &
                     + d7(1)*( rhou_n(i,j,k+1)+rhou_n(i,j,k-1) )&
                     + d7(2)*( rhou_n(i,j,k+2)+rhou_n(i,j,k-2) )&
                     + d7(3)*( rhou_n(i,j,k+3)+rhou_n(i,j,k-3) )

                Krhov(i,j,k)= Krhov(i,j,k)  + d7(0)*rhov_n(i,j,k) &
                     + d7(1)*( rhov_n(i,j,k+1)+rhov_n(i,j,k-1) )&
                     + d7(2)*( rhov_n(i,j,k+2)+rhov_n(i,j,k-2) )&
                     + d7(3)*( rhov_n(i,j,k+3)+rhov_n(i,j,k-3) )

                Krhow(i,j,k)= Krhow(i,j,k)  + d7(0)*rhow_n(i,j,k) &
                     + d7(1)*( rhow_n(i,j,k+1)+rhow_n(i,j,k-1) )&
                     + d7(2)*( rhow_n(i,j,k+2)+rhow_n(i,j,k-2) )&
                     + d7(3)*( rhow_n(i,j,k+3)+rhow_n(i,j,k-3) )

                Krhoe(i,j,k)= Krhoe(i,j,k)  + d7(0)*rhoe_n(i,j,k) &
                     + d7(1)*( rhoe_n(i,j,k+1)+rhoe_n(i,j,k-1) )&
                     + d7(2)*( rhoe_n(i,j,k+2)+rhoe_n(i,j,k-2) )&
                     + d7(3)*( rhoe_n(i,j,k+3)+rhoe_n(i,j,k-3) )
             enddo
          enddo
       enddo
    endif
        
  end subroutine filtering_7pts_overlap
  
  !===============================================================================
  subroutine filtering_5pts_overlap
  !===============================================================================
    !> Apply filtering on 5-point stencils (+ boundaries)
    !> with computation-communication overlap
  !===============================================================================
    use mod_flow
    use mod_filtering
    use mod_comm1
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    ! ---------------------------------------------------------------------------

    ! Filtering along i-direction
    ! ===========================

    ! Perform one-sided communications in i-direction
    ! -----------------------------------------------
    ! open target windows
    call window_comm1
    ! communicate var
    call communication_1_i
    
    ! Interior points
    ! ---------------
    ! 4-th order centered filter
    do k=1,nz
       do j=1,ny
          do i=3,nx-2
             Krho(i,j,k) = d5(0)*rho_n(i,j,k) &
                  + d5(1)*( rho_n(i+1,j,k)+rho_n(i-1,j,k) )&
                  + d5(2)*( rho_n(i+2,j,k)+rho_n(i-2,j,k) )

             Krhou(i,j,k)= d5(0)*rhou_n(i,j,k) &
                  + d5(1)*( rhou_n(i+1,j,k)+rhou_n(i-1,j,k) )&
                  + d5(2)*( rhou_n(i+2,j,k)+rhou_n(i-2,j,k) )

             Krhov(i,j,k)= d5(0)*rhov_n(i,j,k) &
                  + d5(1)*( rhov_n(i+1,j,k)+rhov_n(i-1,j,k) )&
                  + d5(2)*( rhov_n(i+2,j,k)+rhov_n(i-2,j,k) )

             Krhow(i,j,k)= d5(0)*rhow_n(i,j,k) &
                  + d5(1)*( rhow_n(i+1,j,k)+rhow_n(i-1,j,k) )&
                  + d5(2)*( rhow_n(i+2,j,k)+rhow_n(i-2,j,k) )

             Krhoe(i,j,k)= d5(0)*rhoe_n(i,j,k) &
                  + d5(1)*( rhoe_n(i+1,j,k)+rhoe_n(i-1,j,k) )&
                  + d5(2)*( rhoe_n(i+2,j,k)+rhoe_n(i-2,j,k) )
          enddo
       enddo
    enddo
    
    ! Check that communications are fullfilled
    ! ----------------------------------------
    call window_comm1
    
    ! BC at imin
    ! ----------
    if (is_boundary(1,1)) then
       ! no filter for i=1 & i=2 -> initialize array
       do j=1,ny
           Krho(1,j,:)=0.0_wp
          Krhou(1,j,:)=0.0_wp
          Krhov(1,j,:)=0.0_wp
          Krhow(1,j,:)=0.0_wp
          Krhoe(1,j,:)=0.0_wp
          
           Krho(2,j,:)=0.0_wp
          Krhou(2,j,:)=0.0_wp
          Krhov(2,j,:)=0.0_wp
          Krhow(2,j,:)=0.0_wp
          Krhoe(2,j,:)=0.0_wp
       enddo
    else
       ! 4-th order centered filter
       do k=1,nz
          do j=1,ny
             do i=1,2
                Krho(i,j,k) = d5(0)*rho_n(i,j,k) &
                     + d5(1)*( rho_n(i+1,j,k)+rho_n(i-1,j,k) )&
                     + d5(2)*( rho_n(i+2,j,k)+rho_n(i-2,j,k) )

                Krhou(i,j,k)= d5(0)*rhou_n(i,j,k) &
                     + d5(1)*( rhou_n(i+1,j,k)+rhou_n(i-1,j,k) )&
                     + d5(2)*( rhou_n(i+2,j,k)+rhou_n(i-2,j,k) )

                Krhov(i,j,k)= d5(0)*rhov_n(i,j,k) &
                     + d5(1)*( rhov_n(i+1,j,k)+rhov_n(i-1,j,k) )&
                     + d5(2)*( rhov_n(i+2,j,k)+rhov_n(i-2,j,k) )

                Krhow(i,j,k)= d5(0)*rhow_n(i,j,k) &
                     + d5(1)*( rhow_n(i+1,j,k)+rhow_n(i-1,j,k) )&
                     + d5(2)*( rhow_n(i+2,j,k)+rhow_n(i-2,j,k) )

                Krhoe(i,j,k)= d5(0)*rhoe_n(i,j,k) &
                     + d5(1)*( rhoe_n(i+1,j,k)+rhoe_n(i-1,j,k) )&
                     + d5(2)*( rhoe_n(i+2,j,k)+rhoe_n(i-2,j,k) )
             enddo
          enddo
       enddo
    endif

    ! BC at imax
    ! ----------
    if (is_boundary(1,2)) then
       ! no filter for i=nx-1 & i=nx -> initialize array
       do j=1,ny
           Krho(nx-1,j,:)=0.0_wp
          Krhou(nx-1,j,:)=0.0_wp
          Krhov(nx-1,j,:)=0.0_wp
          Krhow(nx-1,j,:)=0.0_wp
          Krhoe(nx-1,j,:)=0.0_wp
          
           Krho(nx,j,:)=0.0_wp
          Krhou(nx,j,:)=0.0_wp
          Krhov(nx,j,:)=0.0_wp
          Krhow(nx,j,:)=0.0_wp
          Krhoe(nx,j,:)=0.0_wp
       enddo
    else
       ! 4-th order centered filter
       do k=1,nz
          do j=1,ny
             do i=nx-1,nx
                Krho(i,j,k) = d5(0)*rho_n(i,j,k) &
                     + d5(1)*( rho_n(i+1,j,k)+rho_n(i-1,j,k) )&
                     + d5(2)*( rho_n(i+2,j,k)+rho_n(i-2,j,k) )

                Krhou(i,j,k)= d5(0)*rhou_n(i,j,k) &
                     + d5(1)*( rhou_n(i+1,j,k)+rhou_n(i-1,j,k) )&
                     + d5(2)*( rhou_n(i+2,j,k)+rhou_n(i-2,j,k) )

                Krhov(i,j,k)= d5(0)*rhov_n(i,j,k) &
                     + d5(1)*( rhov_n(i+1,j,k)+rhov_n(i-1,j,k) )&
                     + d5(2)*( rhov_n(i+2,j,k)+rhov_n(i-2,j,k) )

                Krhow(i,j,k)= d5(0)*rhow_n(i,j,k) &
                     + d5(1)*( rhow_n(i+1,j,k)+rhow_n(i-1,j,k) )&
                     + d5(2)*( rhow_n(i+2,j,k)+rhow_n(i-2,j,k) )

                Krhoe(i,j,k)= d5(0)*rhoe_n(i,j,k) &
                     + d5(1)*( rhoe_n(i+1,j,k)+rhoe_n(i-1,j,k) )&
                     + d5(2)*( rhoe_n(i+2,j,k)+rhoe_n(i-2,j,k) )
             enddo
          enddo
       enddo
    endif
    
    ! Filtering along j-direction
    ! ===========================
    
    ! Perform one-sided communications in j-direction
    ! -----------------------------------------------
    ! open target windows
    call window_comm1
    ! communicate var
    call communication_1_j
    
    ! Interior points
    ! ---------------
    ! 4-th order centered filter
    do k=1,nz
       do j=3,ny-2
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  + d5(0)*rho_n(i,j,k) &
                  + d5(1)*( rho_n(i,j+1,k)+rho_n(i,j-1,k) )&
                  + d5(2)*( rho_n(i,j+2,k)+rho_n(i,j-2,k) )

             Krhou(i,j,k)= Krhou(i,j,k)  + d5(0)*rhou_n(i,j,k) &
                  + d5(1)*( rhou_n(i,j+1,k)+rhou_n(i,j-1,k) )&
                  + d5(2)*( rhou_n(i,j+2,k)+rhou_n(i,j-2,k) )

             Krhov(i,j,k)= Krhov(i,j,k)  + d5(0)*rhov_n(i,j,k) &
                  + d5(1)*( rhov_n(i,j+1,k)+rhov_n(i,j-1,k) )&
                  + d5(2)*( rhov_n(i,j+2,k)+rhov_n(i,j-2,k) )

             Krhow(i,j,k)= Krhow(i,j,k)  + d5(0)*rhow_n(i,j,k) &
                  + d5(1)*( rhow_n(i,j+1,k)+rhow_n(i,j-1,k) )&
                  + d5(2)*( rhow_n(i,j+2,k)+rhow_n(i,j-2,k) )

             Krhoe(i,j,k)= Krhoe(i,j,k)  + d5(0)*rhoe_n(i,j,k) &
                  + d5(1)*( rhoe_n(i,j+1,k)+rhoe_n(i,j-1,k) )&
                  + d5(2)*( rhoe_n(i,j+2,k)+rhoe_n(i,j-2,k) )
          enddo
       enddo
    enddo
    
    ! Check that communications are fullfilled
    ! ----------------------------------------
    call window_comm1

    ! BC at jmin
    ! ----------
    if (.not.is_boundary(2,1)) then
       ! 4-th order centered filter
       do k=1,nz
          do j=1,2
             do i=1,nx
                Krho(i,j,k) = Krho(i,j,k)  + d5(0)*rho_n(i,j,k) &
                     + d5(1)*( rho_n(i,j+1,k)+rho_n(i,j-1,k) )&
                     + d5(2)*( rho_n(i,j+2,k)+rho_n(i,j-2,k) )

                Krhou(i,j,k)= Krhou(i,j,k)  + d5(0)*rhou_n(i,j,k) &
                     + d5(1)*( rhou_n(i,j+1,k)+rhou_n(i,j-1,k) )&
                     + d5(2)*( rhou_n(i,j+2,k)+rhou_n(i,j-2,k) )

                Krhov(i,j,k)= Krhov(i,j,k)  + d5(0)*rhov_n(i,j,k) &
                     + d5(1)*( rhov_n(i,j+1,k)+rhov_n(i,j-1,k) )&
                     + d5(2)*( rhov_n(i,j+2,k)+rhov_n(i,j-2,k) )

                Krhow(i,j,k)= Krhow(i,j,k)  + d5(0)*rhow_n(i,j,k) &
                     + d5(1)*( rhow_n(i,j+1,k)+rhow_n(i,j-1,k) )&
                     + d5(2)*( rhow_n(i,j+2,k)+rhow_n(i,j-2,k) )

                Krhoe(i,j,k)= Krhoe(i,j,k)  + d5(0)*rhoe_n(i,j,k) &
                     + d5(1)*( rhoe_n(i,j+1,k)+rhoe_n(i,j-1,k) )&
                     + d5(2)*( rhoe_n(i,j+2,k)+rhoe_n(i,j-2,k) )
             enddo
          enddo
       enddo
    endif

    ! BC at jmax
    ! ----------
    if (.not.is_boundary(2,2)) then
       ! 4-th order centered filter
       do k=1,nz
          do j=ny-1,ny
             do i=1,nx
                Krho(i,j,k) = Krho(i,j,k)  + d5(0)*rho_n(i,j,k) &
                     + d5(1)*( rho_n(i,j+1,k)+rho_n(i,j-1,k) )&
                     + d5(2)*( rho_n(i,j+2,k)+rho_n(i,j-2,k) )

                Krhou(i,j,k)= Krhou(i,j,k)  + d5(0)*rhou_n(i,j,k) &
                     + d5(1)*( rhou_n(i,j+1,k)+rhou_n(i,j-1,k) )&
                     + d5(2)*( rhou_n(i,j+2,k)+rhou_n(i,j-2,k) )

                Krhov(i,j,k)= Krhov(i,j,k)  + d5(0)*rhov_n(i,j,k) &
                     + d5(1)*( rhov_n(i,j+1,k)+rhov_n(i,j-1,k) )&
                     + d5(2)*( rhov_n(i,j+2,k)+rhov_n(i,j-2,k) )

                Krhow(i,j,k)= Krhow(i,j,k)  + d5(0)*rhow_n(i,j,k) &
                     + d5(1)*( rhow_n(i,j+1,k)+rhow_n(i,j-1,k) )&
                     + d5(2)*( rhow_n(i,j+2,k)+rhow_n(i,j-2,k) )

                Krhoe(i,j,k)= Krhoe(i,j,k)  + d5(0)*rhoe_n(i,j,k) &
                     + d5(1)*( rhoe_n(i,j+1,k)+rhoe_n(i,j-1,k) )&
                     + d5(2)*( rhoe_n(i,j+2,k)+rhoe_n(i,j-2,k) )
             enddo
          enddo
       enddo
    endif

    !****************
    if (is_2D) return
    !****************
    
    ! Filtering along k-direction
    ! ===========================

    ! Perform one-sided communications in k-direction
    ! -----------------------------------------------
    ! open target windows
    call window_comm1
    ! communicate var
    call communication_1_k
    
    ! Interior points
    ! ---------------
    ! 4-th order centered filter
    do k=3,nz-2
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  + d5(0)*rho_n(i,j,k) &
                  + d5(1)*( rho_n(i,j,k+1)+rho_n(i,j,k-1) )&
                  + d5(2)*( rho_n(i,j,k+2)+rho_n(i,j,k-2) )

             Krhou(i,j,k)= Krhou(i,j,k)  + d5(0)*rhou_n(i,j,k) &
                  + d5(1)*( rhou_n(i,j,k+1)+rhou_n(i,j,k-1) )&
                  + d5(2)*( rhou_n(i,j,k+2)+rhou_n(i,j,k-2) )

             Krhov(i,j,k)= Krhov(i,j,k)  + d5(0)*rhov_n(i,j,k) &
                  + d5(1)*( rhov_n(i,j,k+1)+rhov_n(i,j,k-1) )&
                  + d5(2)*( rhov_n(i,j,k+2)+rhov_n(i,j,k-2) )

             Krhow(i,j,k)= Krhow(i,j,k)  + d5(0)*rhow_n(i,j,k) &
                  + d5(1)*( rhow_n(i,j,k+1)+rhow_n(i,j,k-1) )&
                  + d5(2)*( rhow_n(i,j,k+2)+rhow_n(i,j,k-2) )

             Krhoe(i,j,k)= Krhoe(i,j,k)  + d5(0)*rhoe_n(i,j,k) &
                  + d5(1)*( rhoe_n(i,j,k+1)+rhoe_n(i,j,k-1) )&
                  + d5(2)*( rhoe_n(i,j,k+2)+rhoe_n(i,j,k-2) )
          enddo
       enddo
    enddo
    
    ! Check that communications are fullfilled
    ! ----------------------------------------
    call window_comm1

    ! BC at kmin
    ! ----------
    if (.not.is_boundary(3,1)) then
       ! 4-th order centered filter
       do k=1,2
          do j=1,ny
             do i=1,nx
                Krho(i,j,k) = Krho(i,j,k)  + d5(0)*rho_n(i,j,k) &
                     + d5(1)*( rho_n(i,j,k+1)+rho_n(i,j,k-1) )&
                     + d5(2)*( rho_n(i,j,k+2)+rho_n(i,j,k-2) )

                Krhou(i,j,k)= Krhou(i,j,k)  + d5(0)*rhou_n(i,j,k) &
                     + d5(1)*( rhou_n(i,j,k+1)+rhou_n(i,j,k-1) )&
                     + d5(2)*( rhou_n(i,j,k+2)+rhou_n(i,j,k-2) )

                Krhov(i,j,k)= Krhov(i,j,k)  + d5(0)*rhov_n(i,j,k) &
                     + d5(1)*( rhov_n(i,j,k+1)+rhov_n(i,j,k-1) )&
                     + d5(2)*( rhov_n(i,j,k+2)+rhov_n(i,j,k-2) )

                Krhow(i,j,k)= Krhow(i,j,k)  + d5(0)*rhow_n(i,j,k) &
                     + d5(1)*( rhow_n(i,j,k+1)+rhow_n(i,j,k-1) )&
                     + d5(2)*( rhow_n(i,j,k+2)+rhow_n(i,j,k-2) )

                Krhoe(i,j,k)= Krhoe(i,j,k)  + d5(0)*rhoe_n(i,j,k) &
                     + d5(1)*( rhoe_n(i,j,k+1)+rhoe_n(i,j,k-1) )&
                     + d5(2)*( rhoe_n(i,j,k+2)+rhoe_n(i,j,k-2) )
             enddo
          enddo
       enddo
    endif

    ! BC at kmax
    ! ----------
    if (.not.is_boundary(3,2)) then
       ! 4-th order centered filter
       do k=nz-1,nz
          do j=1,ny
             do i=1,nx
                Krho(i,j,k) = Krho(i,j,k)  + d5(0)*rho_n(i,j,k) &
                     + d5(1)*( rho_n(i,j,k+1)+rho_n(i,j,k-1) )&
                     + d5(2)*( rho_n(i,j,k+2)+rho_n(i,j,k-2) )

                Krhou(i,j,k)= Krhou(i,j,k)  + d5(0)*rhou_n(i,j,k) &
                     + d5(1)*( rhou_n(i,j,k+1)+rhou_n(i,j,k-1) )&
                     + d5(2)*( rhou_n(i,j,k+2)+rhou_n(i,j,k-2) )

                Krhov(i,j,k)= Krhov(i,j,k)  + d5(0)*rhov_n(i,j,k) &
                     + d5(1)*( rhov_n(i,j,k+1)+rhov_n(i,j,k-1) )&
                     + d5(2)*( rhov_n(i,j,k+2)+rhov_n(i,j,k-2) )

                Krhow(i,j,k)= Krhow(i,j,k)  + d5(0)*rhow_n(i,j,k) &
                     + d5(1)*( rhow_n(i,j,k+1)+rhow_n(i,j,k-1) )&
                     + d5(2)*( rhow_n(i,j,k+2)+rhow_n(i,j,k-2) )

                Krhoe(i,j,k)= Krhoe(i,j,k)  + d5(0)*rhoe_n(i,j,k) &
                     + d5(1)*( rhoe_n(i,j,k+1)+rhoe_n(i,j,k-1) )&
                     + d5(2)*( rhoe_n(i,j,k+2)+rhoe_n(i,j,k-2) )
             enddo
          enddo
       enddo
    endif
        
  end subroutine filtering_5pts_overlap
  
  !===============================================================================
  subroutine filtering_3pts_overlap
  !===============================================================================
    !> Apply filtering on 3-point stencils (+ boundaries)
    !> with computation-communication overlap
  !===============================================================================
    use mod_flow
    use mod_filtering
    use mod_comm1
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    ! ---------------------------------------------------------------------------

    ! Filtering along i-direction
    ! ===========================

    ! Perform one-sided communications in i-direction
    ! -----------------------------------------------
    ! open target windows
    call window_comm1
    ! communicate var
    call communication_1_i
    
    ! Interior points
    ! ---------------
    ! 2-nd order centered filter
    do k=1,nz
       do j=1,ny
          do i=2,nx-1
             Krho(i,j,k) = d3(0)*rho_n(i,j,k) &
                  + d3(1)*( rho_n(i+1,j,k)+rho_n(i-1,j,k) )

             Krhou(i,j,k)= d3(0)*rhou_n(i,j,k) &
                  + d3(1)*( rhou_n(i+1,j,k)+rhou_n(i-1,j,k) )

             Krhov(i,j,k)= d3(0)*rhov_n(i,j,k) &
                  + d3(1)*( rhov_n(i+1,j,k)+rhov_n(i-1,j,k) )

             Krhow(i,j,k)= d3(0)*rhow_n(i,j,k) &
                  + d3(1)*( rhow_n(i+1,j,k)+rhow_n(i-1,j,k) )

             Krhoe(i,j,k)= d3(0)*rhoe_n(i,j,k) &
                  + d3(1)*( rhoe_n(i+1,j,k)+rhoe_n(i-1,j,k) )
          enddo
       enddo
    enddo
    
    ! Check that communications are fullfilled
    ! ----------------------------------------
    call window_comm1
    
    ! BC at imin
    ! ----------
    if (is_boundary(1,1)) then
       ! no filter for i=1 & i=2 -> initialize array
       do j=1,ny
           Krho(1,j,:)=0.0_wp
          Krhou(1,j,:)=0.0_wp
          Krhov(1,j,:)=0.0_wp
          Krhow(1,j,:)=0.0_wp
          Krhoe(1,j,:)=0.0_wp
       enddo
    else
       ! 2-nd order centered filter
       i=1
       do k=1,nz
          do j=1,ny
             Krho(i,j,k) = d3(0)*rho_n(i,j,k) &
                  + d3(1)*( rho_n(i+1,j,k)+rho_n(i-1,j,k) )

             Krhou(i,j,k)= d3(0)*rhou_n(i,j,k) &
                  + d3(1)*( rhou_n(i+1,j,k)+rhou_n(i-1,j,k) )

             Krhov(i,j,k)= d3(0)*rhov_n(i,j,k) &
                  + d3(1)*( rhov_n(i+1,j,k)+rhov_n(i-1,j,k) )

             Krhow(i,j,k)= d3(0)*rhow_n(i,j,k) &
                  + d3(1)*( rhow_n(i+1,j,k)+rhow_n(i-1,j,k) )

             Krhoe(i,j,k)= d3(0)*rhoe_n(i,j,k) &
                  + d3(1)*( rhoe_n(i+1,j,k)+rhoe_n(i-1,j,k) )
          enddo
       enddo
    endif

    ! BC at imax
    ! ----------
    if (is_boundary(1,2)) then
       ! no filter for i=nx-1 & i=nx -> initialize array
       do j=1,ny
           Krho(nx,j,:)=0.0_wp
          Krhou(nx,j,:)=0.0_wp
          Krhov(nx,j,:)=0.0_wp
          Krhow(nx,j,:)=0.0_wp
          Krhoe(nx,j,:)=0.0_wp
       enddo
    else
       ! 2-nd order centered filter
       i=nx
       do k=1,nz
          do j=1,ny
             Krho(i,j,k) = d3(0)*rho_n(i,j,k) &
                  + d3(1)*( rho_n(i+1,j,k)+rho_n(i-1,j,k) )

             Krhou(i,j,k)= d3(0)*rhou_n(i,j,k) &
                  + d3(1)*( rhou_n(i+1,j,k)+rhou_n(i-1,j,k) )

             Krhov(i,j,k)= d3(0)*rhov_n(i,j,k) &
                  + d3(1)*( rhov_n(i+1,j,k)+rhov_n(i-1,j,k) )

             Krhow(i,j,k)= d3(0)*rhow_n(i,j,k) &
                  + d3(1)*( rhow_n(i+1,j,k)+rhow_n(i-1,j,k) )

             Krhoe(i,j,k)= d3(0)*rhoe_n(i,j,k) &
                  + d3(1)*( rhoe_n(i+1,j,k)+rhoe_n(i-1,j,k) )
          enddo
       enddo
    endif
    
    ! Filtering along j-direction
    ! ===========================
    
    ! Perform one-sided communications in j-direction
    ! -----------------------------------------------
    ! open target windows
    call window_comm1
    ! communicate var
    call communication_1_j
    
    ! Interior points
    ! ---------------
    ! 2-nd order centered filter
    do k=1,nz
       do j=2,ny-1
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  + d3(0)*rho_n(i,j,k) &
                  + d3(1)*( rho_n(i,j+1,k)+rho_n(i,j-1,k) )

             Krhou(i,j,k)= Krhou(i,j,k)  + d3(0)*rhou_n(i,j,k) &
                  + d3(1)*( rhou_n(i,j+1,k)+rhou_n(i,j-1,k) )

             Krhov(i,j,k)= Krhov(i,j,k)  + d3(0)*rhov_n(i,j,k) &
                  + d3(1)*( rhov_n(i,j+1,k)+rhov_n(i,j-1,k) )

             Krhow(i,j,k)= Krhow(i,j,k)  + d3(0)*rhow_n(i,j,k) &
                  + d3(1)*( rhow_n(i,j+1,k)+rhow_n(i,j-1,k) )

             Krhoe(i,j,k)= Krhoe(i,j,k)  + d3(0)*rhoe_n(i,j,k) &
                  + d3(1)*( rhoe_n(i,j+1,k)+rhoe_n(i,j-1,k) )
          enddo
       enddo
    enddo
    
    ! Check that communications are fullfilled
    ! ----------------------------------------
    call window_comm1

    ! BC at jmin
    ! ----------
    if (.not.is_boundary(2,1)) then
       ! 2-nd order centered filter
       j=1
       do k=1,nz
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  + d3(0)*rho_n(i,j,k) &
                  + d3(1)*( rho_n(i,j+1,k)+rho_n(i,j-1,k) )

             Krhou(i,j,k)= Krhou(i,j,k)  + d3(0)*rhou_n(i,j,k) &
                  + d3(1)*( rhou_n(i,j+1,k)+rhou_n(i,j-1,k) )

             Krhov(i,j,k)= Krhov(i,j,k)  + d3(0)*rhov_n(i,j,k) &
                  + d3(1)*( rhov_n(i,j+1,k)+rhov_n(i,j-1,k) )

             Krhow(i,j,k)= Krhow(i,j,k)  + d3(0)*rhow_n(i,j,k) &
                  + d3(1)*( rhow_n(i,j+1,k)+rhow_n(i,j-1,k) )

             Krhoe(i,j,k)= Krhoe(i,j,k)  + d3(0)*rhoe_n(i,j,k) &
                  + d3(1)*( rhoe_n(i,j+1,k)+rhoe_n(i,j-1,k) )
          enddo
       enddo
    endif

    ! BC at jmax
    ! ----------
    if (.not.is_boundary(2,2)) then
       ! 2-nd order centered filter
       j=ny
       do k=1,nz
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  + d3(0)*rho_n(i,j,k) &
                  + d3(1)*( rho_n(i,j+1,k)+rho_n(i,j-1,k) )

             Krhou(i,j,k)= Krhou(i,j,k)  + d3(0)*rhou_n(i,j,k) &
                  + d3(1)*( rhou_n(i,j+1,k)+rhou_n(i,j-1,k) )

             Krhov(i,j,k)= Krhov(i,j,k)  + d3(0)*rhov_n(i,j,k) &
                  + d3(1)*( rhov_n(i,j+1,k)+rhov_n(i,j-1,k) )

             Krhow(i,j,k)= Krhow(i,j,k)  + d3(0)*rhow_n(i,j,k) &
                  + d3(1)*( rhow_n(i,j+1,k)+rhow_n(i,j-1,k) )

             Krhoe(i,j,k)= Krhoe(i,j,k)  + d3(0)*rhoe_n(i,j,k) &
                  + d3(1)*( rhoe_n(i,j+1,k)+rhoe_n(i,j-1,k) )
          enddo
       enddo
    endif

    !****************
    if (is_2D) return
    !****************
    
    ! Filtering along k-direction
    ! ===========================

    ! Perform one-sided communications in k-direction
    ! -----------------------------------------------
    ! open target windows
    call window_comm1
    ! communicate var
    call communication_1_k
    
    ! Interior points
    ! ---------------
    ! 2-nd order centered filter
    do k=2,nz-1
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  + d3(0)*rho_n(i,j,k) &
                  + d3(1)*( rho_n(i,j,k+1)+rho_n(i,j,k-1) )

             Krhou(i,j,k)= Krhou(i,j,k)  + d3(0)*rhou_n(i,j,k) &
                  + d3(1)*( rhou_n(i,j,k+1)+rhou_n(i,j,k-1) )

             Krhov(i,j,k)= Krhov(i,j,k)  + d3(0)*rhov_n(i,j,k) &
                  + d3(1)*( rhov_n(i,j,k+1)+rhov_n(i,j,k-1) )

             Krhow(i,j,k)= Krhow(i,j,k)  + d3(0)*rhow_n(i,j,k) &
                  + d3(1)*( rhow_n(i,j,k+1)+rhow_n(i,j,k-1) )

             Krhoe(i,j,k)= Krhoe(i,j,k)  + d3(0)*rhoe_n(i,j,k) &
                  + d3(1)*( rhoe_n(i,j,k+1)+rhoe_n(i,j,k-1) )
          enddo
       enddo
    enddo
    
    ! Check that communications are fullfilled
    ! ----------------------------------------
    call window_comm1

    ! BC at kmin
    ! ----------
    if (.not.is_boundary(3,1)) then
       ! 2-nd order centered filter
       k=1
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  + d3(0)*rho_n(i,j,k) &
                  + d3(1)*( rho_n(i,j,k+1)+rho_n(i,j,k-1) )

             Krhou(i,j,k)= Krhou(i,j,k)  + d3(0)*rhou_n(i,j,k) &
                  + d3(1)*( rhou_n(i,j,k+1)+rhou_n(i,j,k-1) )

             Krhov(i,j,k)= Krhov(i,j,k)  + d3(0)*rhov_n(i,j,k) &
                  + d3(1)*( rhov_n(i,j,k+1)+rhov_n(i,j,k-1) )

             Krhow(i,j,k)= Krhow(i,j,k)  + d3(0)*rhow_n(i,j,k) &
                  + d3(1)*( rhow_n(i,j,k+1)+rhow_n(i,j,k-1) )

             Krhoe(i,j,k)= Krhoe(i,j,k)  + d3(0)*rhoe_n(i,j,k) &
                  + d3(1)*( rhoe_n(i,j,k+1)+rhoe_n(i,j,k-1) )
          enddo
       enddo
    endif

    ! BC at kmax
    ! ----------
    if (.not.is_boundary(3,2)) then
       ! 2-nd order centered filter
       k=nz
       do j=1,ny
          do i=1,nx
             Krho(i,j,k) = Krho(i,j,k)  + d3(0)*rho_n(i,j,k) &
                  + d3(1)*( rho_n(i,j,k+1)+rho_n(i,j,k-1) )

             Krhou(i,j,k)= Krhou(i,j,k)  + d3(0)*rhou_n(i,j,k) &
                  + d3(1)*( rhou_n(i,j,k+1)+rhou_n(i,j,k-1) )

             Krhov(i,j,k)= Krhov(i,j,k)  + d3(0)*rhov_n(i,j,k) &
                  + d3(1)*( rhov_n(i,j,k+1)+rhov_n(i,j,k-1) )

             Krhow(i,j,k)= Krhow(i,j,k)  + d3(0)*rhow_n(i,j,k) &
                  + d3(1)*( rhow_n(i,j,k+1)+rhow_n(i,j,k-1) )

             Krhoe(i,j,k)= Krhoe(i,j,k)  + d3(0)*rhoe_n(i,j,k) &
                  + d3(1)*( rhoe_n(i,j,k+1)+rhoe_n(i,j,k-1) )
          enddo
       enddo
    endif
        
  end subroutine filtering_3pts_overlap
