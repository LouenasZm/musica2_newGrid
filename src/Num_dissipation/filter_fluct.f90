  !===============================================================================
  subroutine filtering_fluct_11pts
  !===============================================================================
    !> Apply filtering on 11-point stencils (+ boundaries)
  !===============================================================================
    use mod_flow
    use mod_flow0
    use mod_filtering
    use mod_eos
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2) :: rhof,rhouf,rhovf,rhowf,rhoef
    ! ---------------------------------------------------------------------------
    
    rhof  = rho_n  - rho0
    rhouf = rhou_n - rho0*u0 
    rhovf = rhov_n - rho0*v0 
    rhowf = rhow_n - rho0*w0 
    rhoef = rhoe_n - ( p0/(gam-1)+0.5*rho0*(u0**2+v0**2+w0**2) )

    ! Filtering along x
    ! =================
    ! BC at imin
    ! ----------
    if (is_boundary(1,1)) then
!!$       ! 2-th order centered filter
!!$       i=2
!!$       do k=1,nz
!!$          do j=1,ny
!!$             Krho(i,j,k)  = Krho(i,j,k)  +d3(0)*rhof(i,j,k)  &
!!$                  + d3(1)*rhof(i+1,j,k)+d3(1)*rhof(i-1,j,k)
!!$             Krhou(i,j,k) = Krhou(i,j,k) +d3(0)*rhouf(i,j,k) &
!!$                  + d3(1)*rhouf(i+1,j,k)+d3(1)*rhouf(i-1,j,k)
!!$             Krhov(i,j,k) = Krhov(i,j,k) +d3(0)*rhovf(i,j,k) &
!!$                  + d3(1)*rhovf(i+1,j,k)+d3(1)*rhovf(i-1,j,k)
!!$             Krhow(i,j,k) = Krhow(i,j,k) +d3(0)*rhowf(i,j,k) &
!!$                  + d3(1)*rhowf(i+1,j,k)+d3(1)*rhowf(i-1,j,k)
!!$             Krhoe(i,j,k) = Krhoe(i,j,k) +d3(0)*rhoef(i,j,k) &
!!$                  + d3(1)*rhoef(i+1,j,k)+d3(1)*rhoef(i-1,j,k)
!!$          enddo
!!$       enddo
!!$       ! 4-th order centered filter
!!$       i=3
!!$       do k=1,nz
!!$          do j=1,ny
!!$             Krho(i,j,k)  = Krho(i,j,k)  +d5(0)*rhof(i,j,k)   &
!!$                  + d5(1)*rhof(i+1,j,k)+d5(2)*rhof(i+2,j,k) &
!!$                  + d5(1)*rhof(i-1,j,k)+d5(2)*rhof(i-2,j,k)
!!$             Krhou(i,j,k) = Krhou(i,j,k) +d5(0)*rhouf(i,j,k)   &
!!$                  + d5(1)*rhouf(i+1,j,k)+d5(2)*rhouf(i+2,j,k) &
!!$                  + d5(1)*rhouf(i-1,j,k)+d5(2)*rhouf(i-2,j,k)
!!$             Krhov(i,j,k) = Krhov(i,j,k) +d5(0)*rhovf(i,j,k)   &
!!$                  + d5(1)*rhovf(i+1,j,k)+d5(2)*rhovf(i+2,j,k) &
!!$                  + d5(1)*rhovf(i-1,j,k)+d5(2)*rhovf(i-2,j,k)
!!$             Krhow(i,j,k) = Krhow(i,j,k) +d5(0)*rhowf(i,j,k)   &
!!$                  + d5(1)*rhowf(i+1,j,k)+d5(2)*rhowf(i+2,j,k) &
!!$                  + d5(1)*rhowf(i-1,j,k)+d5(2)*rhowf(i-2,j,k)
!!$             Krhoe(i,j,k) = Krhoe(i,j,k) +d5(0)*rhoef(i,j,k)   &
!!$                  + d5(1)*rhoef(i+1,j,k)+d5(2)*rhoef(i+2,j,k) &
!!$                  + d5(1)*rhoef(i-1,j,k)+d5(2)*rhoef(i-2,j,k)
!!$          enddo
!!$       enddo
       ! 6-th order centered filter
       i=4
       do k=1,nz
          do j=1,ny
             Krho(i,j,k)  = Krho(i,j,k)  +d7(0)*rhof(i,j,k)   &
                  + d7(1)*rhof(i+1,j,k)+d7(2)*rhof(i+2,j,k)+d7(3)*rhof(i+3,j,k) &
                  + d7(1)*rhof(i-1,j,k)+d7(2)*rhof(i-2,j,k)+d7(3)*rhof(i-3,j,k)
             Krhou(i,j,k) = Krhou(i,j,k) +d7(0)*rhouf(i,j,k)   &
                  + d7(1)*rhouf(i+1,j,k)+d7(2)*rhouf(i+2,j,k)+d7(3)*rhouf(i+3,j,k) &
                  + d7(1)*rhouf(i-1,j,k)+d7(2)*rhouf(i-2,j,k)+d7(3)*rhouf(i-3,j,k)
             Krhov(i,j,k) = Krhov(i,j,k) +d7(0)*rhovf(i,j,k)   &
                  + d7(1)*rhovf(i+1,j,k)+d7(2)*rhovf(i+2,j,k)+d7(3)*rhovf(i+3,j,k) &
                  + d7(1)*rhovf(i-1,j,k)+d7(2)*rhovf(i-2,j,k)+d7(3)*rhovf(i-3,j,k)
             Krhow(i,j,k) = Krhow(i,j,k) +d7(0)*rhowf(i,j,k)   &
                  + d7(1)*rhowf(i+1,j,k)+d7(2)*rhowf(i+2,j,k)+d7(3)*rhowf(i+3,j,k) &
                  + d7(1)*rhowf(i-1,j,k)+d7(2)*rhowf(i-2,j,k)+d7(3)*rhowf(i-3,j,k)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d7(0)*rhoef(i,j,k)   &
                  + d7(1)*rhoef(i+1,j,k)+d7(2)*rhoef(i+2,j,k)+d7(3)*rhoef(i+3,j,k) &
                  + d7(1)*rhoef(i-1,j,k)+d7(2)*rhoef(i-2,j,k)+d7(3)*rhoef(i-3,j,k)
          enddo
       enddo
       ! 8-th order centered filter
       i=5
       do k=1,nz
          do j=1,ny
             Krho(i,j,k)  = Krho(i,j,k)  +d9(0)*rhof(i,j,k)   &
                  + d9(1)*rhof(i+1,j,k)+d9(2)*rhof(i+2,j,k)+d9(3)*rhof(i+3,j,k)+d9(4)*rhof(i+4,j,k) &
                  + d9(1)*rhof(i-1,j,k)+d9(2)*rhof(i-2,j,k)+d9(3)*rhof(i-3,j,k)+d9(4)*rhof(i-4,j,k)
             Krhou(i,j,k) = Krhou(i,j,k) +d9(0)*rhouf(i,j,k)   &
                  + d9(1)*rhouf(i+1,j,k)+d9(2)*rhouf(i+2,j,k)+d9(3)*rhouf(i+3,j,k)+d9(4)*rhouf(i+4,j,k) &
                  + d9(1)*rhouf(i-1,j,k)+d9(2)*rhouf(i-2,j,k)+d9(3)*rhouf(i-3,j,k)+d9(4)*rhouf(i-4,j,k)
             Krhov(i,j,k) = Krhov(i,j,k) +d9(0)*rhovf(i,j,k)   &
                  + d9(1)*rhovf(i+1,j,k)+d9(2)*rhovf(i+2,j,k)+d9(3)*rhovf(i+3,j,k)+d9(4)*rhovf(i+4,j,k) &
                  + d9(1)*rhovf(i-1,j,k)+d9(2)*rhovf(i-2,j,k)+d9(3)*rhovf(i-3,j,k)+d9(4)*rhovf(i-4,j,k)
             Krhow(i,j,k) = Krhow(i,j,k) +d9(0)*rhowf(i,j,k)   &
                  + d9(1)*rhowf(i+1,j,k)+d9(2)*rhowf(i+2,j,k)+d9(3)*rhowf(i+3,j,k)+d9(4)*rhowf(i+4,j,k) &
                  + d9(1)*rhowf(i-1,j,k)+d9(2)*rhowf(i-2,j,k)+d9(3)*rhowf(i-3,j,k)+d9(4)*rhowf(i-4,j,k)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d9(0)*rhoef(i,j,k)   &
                  + d9(1)*rhoef(i+1,j,k)+d9(2)*rhoef(i+2,j,k)+d9(3)*rhoef(i+3,j,k)+d9(4)*rhoef(i+4,j,k) &
                  + d9(1)*rhoef(i-1,j,k)+d9(2)*rhoef(i-2,j,k)+d9(3)*rhoef(i-3,j,k)+d9(4)*rhoef(i-4,j,k)
          enddo
       enddo
    endif

    ! Interior points
    ! ---------------
    ! 10-th order centered filter
    do k=1,nz
       do j=1,ny
          do i=ndx,nfx
             Krho(i,j,k)  = Krho(i,j,k)  + d11(0)*rhof(i,j,k) &
                  + d11(1)*( rhof(i+1,j,k)+rhof(i-1,j,k) )&
                  + d11(2)*( rhof(i+2,j,k)+rhof(i-2,j,k) )&
                  + d11(3)*( rhof(i+3,j,k)+rhof(i-3,j,k) )&
                  + d11(4)*( rhof(i+4,j,k)+rhof(i-4,j,k) )&
                  + d11(5)*( rhof(i+5,j,k)+rhof(i-5,j,k) )

             Krhou(i,j,k) = Krhou(i,j,k)  + d11(0)*rhouf(i,j,k) &
                  + d11(1)*( rhouf(i+1,j,k)+rhouf(i-1,j,k) )&
                  + d11(2)*( rhouf(i+2,j,k)+rhouf(i-2,j,k) )&
                  + d11(3)*( rhouf(i+3,j,k)+rhouf(i-3,j,k) )&
                  + d11(4)*( rhouf(i+4,j,k)+rhouf(i-4,j,k) )&
                  + d11(5)*( rhouf(i+5,j,k)+rhouf(i-5,j,k) )

             Krhov(i,j,k) = Krhov(i,j,k)  + d11(0)*rhovf(i,j,k) &
                  + d11(1)*( rhovf(i+1,j,k)+rhovf(i-1,j,k) )&
                  + d11(2)*( rhovf(i+2,j,k)+rhovf(i-2,j,k) )&
                  + d11(3)*( rhovf(i+3,j,k)+rhovf(i-3,j,k) )&
                  + d11(4)*( rhovf(i+4,j,k)+rhovf(i-4,j,k) )&
                  + d11(5)*( rhovf(i+5,j,k)+rhovf(i-5,j,k) )

             Krhow(i,j,k) = Krhow(i,j,k)  + d11(0)*rhowf(i,j,k) &
                  + d11(1)*( rhowf(i+1,j,k)+rhowf(i-1,j,k) )&
                  + d11(2)*( rhowf(i+2,j,k)+rhowf(i-2,j,k) )&
                  + d11(3)*( rhowf(i+3,j,k)+rhowf(i-3,j,k) )&
                  + d11(4)*( rhowf(i+4,j,k)+rhowf(i-4,j,k) )&
                  + d11(5)*( rhowf(i+5,j,k)+rhowf(i-5,j,k) )

             Krhoe(i,j,k) = Krhoe(i,j,k)  + d11(0)*rhoef(i,j,k) &
                  + d11(1)*( rhoef(i+1,j,k)+rhoef(i-1,j,k) )&
                  + d11(2)*( rhoef(i+2,j,k)+rhoef(i-2,j,k) )&
                  + d11(3)*( rhoef(i+3,j,k)+rhoef(i-3,j,k) )&
                  + d11(4)*( rhoef(i+4,j,k)+rhoef(i-4,j,k) )&
                  + d11(5)*( rhoef(i+5,j,k)+rhoef(i-5,j,k) )
          enddo
       enddo
    enddo

    ! BC at imax
    ! ----------
    if (is_boundary(1,2)) then
       ! 8-th order centered filter
       i=nx-4
       do k=1,nz
          do j=1,ny
             Krho(i,j,k)  = Krho(i,j,k)  +d9(0)*rhof(i,j,k)   &
                  + d9(1)*rhof(i+1,j,k)+d9(2)*rhof(i+2,j,k)+d9(3)*rhof(i+3,j,k)+d9(4)*rhof(i+4,j,k) &
                  + d9(1)*rhof(i-1,j,k)+d9(2)*rhof(i-2,j,k)+d9(3)*rhof(i-3,j,k)+d9(4)*rhof(i-4,j,k)
             Krhou(i,j,k) = Krhou(i,j,k) +d9(0)*rhouf(i,j,k)   &
                  + d9(1)*rhouf(i+1,j,k)+d9(2)*rhouf(i+2,j,k)+d9(3)*rhouf(i+3,j,k)+d9(4)*rhouf(i+4,j,k) &
                  + d9(1)*rhouf(i-1,j,k)+d9(2)*rhouf(i-2,j,k)+d9(3)*rhouf(i-3,j,k)+d9(4)*rhouf(i-4,j,k)
             Krhov(i,j,k) = Krhov(i,j,k) +d9(0)*rhovf(i,j,k)   &
                  + d9(1)*rhovf(i+1,j,k)+d9(2)*rhovf(i+2,j,k)+d9(3)*rhovf(i+3,j,k)+d9(4)*rhovf(i+4,j,k) &
                  + d9(1)*rhovf(i-1,j,k)+d9(2)*rhovf(i-2,j,k)+d9(3)*rhovf(i-3,j,k)+d9(4)*rhovf(i-4,j,k)
             Krhow(i,j,k) = Krhow(i,j,k) +d9(0)*rhowf(i,j,k)   &
                  + d9(1)*rhowf(i+1,j,k)+d9(2)*rhowf(i+2,j,k)+d9(3)*rhowf(i+3,j,k)+d9(4)*rhowf(i+4,j,k) &
                  + d9(1)*rhowf(i-1,j,k)+d9(2)*rhowf(i-2,j,k)+d9(3)*rhowf(i-3,j,k)+d9(4)*rhowf(i-4,j,k)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d9(0)*rhoef(i,j,k)   &
                  + d9(1)*rhoef(i+1,j,k)+d9(2)*rhoef(i+2,j,k)+d9(3)*rhoef(i+3,j,k)+d9(4)*rhoef(i+4,j,k) &
                  + d9(1)*rhoef(i-1,j,k)+d9(2)*rhoef(i-2,j,k)+d9(3)*rhoef(i-3,j,k)+d9(4)*rhoef(i-4,j,k)
          enddo
       enddo
       ! 6-th order centered filter
       i=nx-3
       do k=1,nz
          do j=1,ny
             Krho(i,j,k)  = Krho(i,j,k)  +d7(0)*rhof(i,j,k)   &
                  + d7(1)*rhof(i+1,j,k)+d7(2)*rhof(i+2,j,k)+d7(3)*rhof(i+3,j,k) &
                  + d7(1)*rhof(i-1,j,k)+d7(2)*rhof(i-2,j,k)+d7(3)*rhof(i-3,j,k)
             Krhou(i,j,k) = Krhou(i,j,k) +d7(0)*rhouf(i,j,k)   &
                  + d7(1)*rhouf(i+1,j,k)+d7(2)*rhouf(i+2,j,k)+d7(3)*rhouf(i+3,j,k) &
                  + d7(1)*rhouf(i-1,j,k)+d7(2)*rhouf(i-2,j,k)+d7(3)*rhouf(i-3,j,k)
             Krhov(i,j,k) = Krhov(i,j,k) +d7(0)*rhovf(i,j,k)   &
                  + d7(1)*rhovf(i+1,j,k)+d7(2)*rhovf(i+2,j,k)+d7(3)*rhovf(i+3,j,k) &
                  + d7(1)*rhovf(i-1,j,k)+d7(2)*rhovf(i-2,j,k)+d7(3)*rhovf(i-3,j,k)
             Krhow(i,j,k) = Krhow(i,j,k) +d7(0)*rhowf(i,j,k)   &
                  + d7(1)*rhowf(i+1,j,k)+d7(2)*rhowf(i+2,j,k)+d7(3)*rhowf(i+3,j,k) &
                  + d7(1)*rhowf(i-1,j,k)+d7(2)*rhowf(i-2,j,k)+d7(3)*rhowf(i-3,j,k)
             Krhoe(i,j,k) = Krhoe(i,j,k) +d7(0)*rhoef(i,j,k)   &
                  + d7(1)*rhoef(i+1,j,k)+d7(2)*rhoef(i+2,j,k)+d7(3)*rhoef(i+3,j,k) &
                  + d7(1)*rhoef(i-1,j,k)+d7(2)*rhoef(i-2,j,k)+d7(3)*rhoef(i-3,j,k)
          enddo
       enddo
!!$       ! 4-th order centered filter
!!$       i=nx-2
!!$       do k=1,nz
!!$          do j=1,ny
!!$             Krho(i,j,k)  = Krho(i,j,k)  +d5(0)*rhof(i,j,k)   &
!!$                  + d5(1)*rhof(i+1,j,k)+d5(2)*rhof(i+2,j,k) &
!!$                  + d5(1)*rhof(i-1,j,k)+d5(2)*rhof(i-2,j,k)
!!$             Krhou(i,j,k) = Krhou(i,j,k) +d5(0)*rhouf(i,j,k)   &
!!$                  + d5(1)*rhouf(i+1,j,k)+d5(2)*rhouf(i+2,j,k) &
!!$                  + d5(1)*rhouf(i-1,j,k)+d5(2)*rhouf(i-2,j,k)
!!$             Krhov(i,j,k) = Krhov(i,j,k) +d5(0)*rhovf(i,j,k)   &
!!$                  + d5(1)*rhovf(i+1,j,k)+d5(2)*rhovf(i+2,j,k) &
!!$                  + d5(1)*rhovf(i-1,j,k)+d5(2)*rhovf(i-2,j,k)
!!$             Krhow(i,j,k) = Krhow(i,j,k) +d5(0)*rhowf(i,j,k)   &
!!$                  + d5(1)*rhowf(i+1,j,k)+d5(2)*rhowf(i+2,j,k) &
!!$                  + d5(1)*rhowf(i-1,j,k)+d5(2)*rhowf(i-2,j,k)
!!$             Krhoe(i,j,k) = Krhoe(i,j,k) +d5(0)*rhoef(i,j,k)   &
!!$                  + d5(1)*rhoef(i+1,j,k)+d5(2)*rhoef(i+2,j,k) &
!!$                  + d5(1)*rhoef(i-1,j,k)+d5(2)*rhoef(i-2,j,k)
!!$          enddo
!!$       enddo
!!$       ! 2-th order centered filter
!!$       i=nx-1
!!$       do k=1,nz
!!$          do j=1,ny
!!$             Krho(i,j,k)  = Krho(i,j,k)  +d3(0)*rhof(i,j,k)  &
!!$                  + d3(1)*rhof(i+1,j,k)+d3(1)*rhof(i-1,j,k)
!!$             Krhou(i,j,k) = Krhou(i,j,k) +d3(0)*rhouf(i,j,k) &
!!$                  + d3(1)*rhouf(i+1,j,k)+d3(1)*rhouf(i-1,j,k)
!!$             Krhov(i,j,k) = Krhov(i,j,k) +d3(0)*rhovf(i,j,k) &
!!$                  + d3(1)*rhovf(i+1,j,k)+d3(1)*rhovf(i-1,j,k)
!!$             Krhow(i,j,k) = Krhow(i,j,k) +d3(0)*rhowf(i,j,k) &
!!$                  + d3(1)*rhowf(i+1,j,k)+d3(1)*rhowf(i-1,j,k)
!!$             Krhoe(i,j,k) = Krhoe(i,j,k) +d3(0)*rhoef(i,j,k) &
!!$                  + d3(1)*rhoef(i+1,j,k)+d3(1)*rhoef(i-1,j,k)
!!$          enddo
!!$       enddo
    endif

    ! Filtering along y
    ! =================
    ! BC at jmin
    ! ----------
    if (is_boundary(2,1)) then
!!$       ! non-centered filtering on stencil 1-5
!!$       !j=2
!!$       !do k=1,nz
!!$       !   do i=1,nx
!!$       !      Krho(i,j,k)  = Krho(i,j,k) +d15(1)*rhof(i,j-1,k) &
!!$       !           + d15(2)*rhof(i,j,k)   +d15(3)*rhof(i,j+1,k) &
!!$       !           + d15(4)*rhof(i,j+2,k) +d15(5)*rhof(i,j+3,k) &
!!$       !           + d15(6)*rhof(i,j+4,k) +d15(7)*rhof(i,j+5,k)
!!$       !      Krhou(i,j,k) = Krhou(i,j,k) +d15(1)*rhouf(i,j-1,k) &
!!$       !           + d15(2)*rhouf(i,j,k)  + d15(3)*rhouf(i,j+1,k) &
!!$       !           + d15(4)*rhouf(i,j+2,k)+ d15(5)*rhouf(i,j+3,k) &
!!$       !           + d15(6)*rhouf(i,j+4,k)+ d15(7)*rhouf(i,j+5,k)
!!$       !      Krhov(i,j,k) = Krhov(i,j,k) +d15(1)*rhovf(i,j-1,k) &
!!$       !           + d15(2)*rhovf(i,j,k)  + d15(3)*rhovf(i,j+1,k) &
!!$       !           + d15(4)*rhovf(i,j+2,k)+ d15(5)*rhovf(i,j+3,k) &
!!$       !           + d15(6)*rhovf(i,j+4,k)+ d15(7)*rhovf(i,j+5,k)
!!$       !      Krhow(i,j,k) = Krhow(i,j,k) +d15(1)*rhowf(i,j-1,k) &
!!$       !           + d15(2)*rhowf(i,j,k)  + d15(3)*rhowf(i,j+1,k) &
!!$       !           + d15(4)*rhowf(i,j+2,k)+ d15(5)*rhowf(i,j+3,k) &
!!$       !           + d15(6)*rhowf(i,j+4,k)+ d15(7)*rhowf(i,j+5,k)
!!$       !      Krhoe(i,j,k) = Krhoe(i,j,k) +d15(1)*rhoef(i,j-1,k) &
!!$       !           + d15(2)*rhoef(i,j,k)  + d15(3)*rhoef(i,j+1,k) &
!!$       !           + d15(4)*rhoef(i,j+2,k)+ d15(5)*rhoef(i,j+3,k) &
!!$       !           + d15(6)*rhoef(i,j+4,k)+ d15(7)*rhoef(i,j+5,k)
!!$       !   enddo
!!$       !enddo
!!$       ! non-centered filtering on stencil 2-8
!!$       j=3
!!$       do k=1,nz
!!$          do i=1,nx
!!$             Krho(i,j,k)  = Krho(i,j,k) + d28(1)*rhof(i,j-2,k) &
!!$                  + d28(2)*rhof(i,j-1,k) + d28(3)*rhof(i,j,k)   &
!!$                  + d28(4)*rhof(i,j+1,k) + d28(5)*rhof(i,j+2,k) &
!!$                  + d28(6)*rhof(i,j+3,k) + d28(7)*rhof(i,j+4,k) &
!!$                  + d28(8)*rhof(i,j+5,k) + d28(9)*rhof(i,j+6,k) &
!!$                  + d28(10)*rhof(i,j+7,k)+ d28(11)*rhof(i,j+8,k)
!!$             Krhou(i,j,k) = Krhou(i,j,k) +d28(1)*rhouf(i,j-2,k) &
!!$                  + d28(2)*rhouf(i,j-1,k) +d28(3)*rhouf(i,j,k)   &
!!$                  + d28(4)*rhouf(i,j+1,k) +d28(5)*rhouf(i,j+2,k) &
!!$                  + d28(6)*rhouf(i,j+3,k) +d28(7)*rhouf(i,j+4,k) &
!!$                  + d28(8)*rhouf(i,j+5,k) +d28(9)*rhouf(i,j+6,k) &
!!$                  + d28(10)*rhouf(i,j+7,k)+d28(11)*rhouf(i,j+8,k)
!!$             Krhov(i,j,k) = Krhov(i,j,k) +d28(1)*rhovf(i,j-2,k) &
!!$                  + d28(2)*rhovf(i,j-1,k) +d28(3)*rhovf(i,j,k)   &
!!$                  + d28(4)*rhovf(i,j+1,k) +d28(5)*rhovf(i,j+2,k) &
!!$                  + d28(6)*rhovf(i,j+3,k) +d28(7)*rhovf(i,j+4,k) &
!!$                  + d28(8)*rhovf(i,j+5,k) +d28(9)*rhovf(i,j+6,k) &
!!$                  + d28(10)*rhovf(i,j+7,k)+d28(11)*rhovf(i,j+8,k)
!!$             Krhow(i,j,k) = Krhow(i,j,k) +d28(1)*rhowf(i,j-2,k) &
!!$                  + d28(2)*rhowf(i,j-1,k) +d28(3)*rhowf(i,j,k)   &
!!$                  + d28(4)*rhowf(i,j+1,k) +d28(5)*rhowf(i,j+2,k) &
!!$                  + d28(6)*rhowf(i,j+3,k) +d28(7)*rhowf(i,j+4,k) &
!!$                  + d28(8)*rhowf(i,j+5,k) +d28(9)*rhowf(i,j+6,k) &
!!$                  + d28(10)*rhowf(i,j+7,k)+d28(11)*rhowf(i,j+8,k)
!!$             Krhoe(i,j,k) = Krhoe(i,j,k) +d28(1)*rhoef(i,j-2,k) &
!!$                  + d28(2)*rhoef(i,j-1,k) +d28(3)*rhoef(i,j,k)   &
!!$                  + d28(4)*rhoef(i,j+1,k) +d28(5)*rhoef(i,j+2,k) &
!!$                  + d28(6)*rhoef(i,j+3,k) +d28(7)*rhoef(i,j+4,k) &
!!$                  + d28(8)*rhoef(i,j+5,k) +d28(9)*rhoef(i,j+6,k) &
!!$                  + d28(10)*rhoef(i,j+7,k)+d28(11)*rhoef(i,j+8,k)
!!$          enddo
!!$       enddo
!!$       ! non-centered filtering on stencil 3-7
!!$       j=4
!!$       do k=1,nz
!!$          do i=1,nx
!!$             Krho(i,j,k)  = Krho(i,j,k) + d37(1)*rhof(i,j-3,k) &
!!$                  + d37(2)*rhof(i,j-2,k) + d37(3)*rhof(i,j-1,k) &
!!$                  + d37(4)*rhof(i,j,k)   + d37(5)*rhof(i,j+1,k) &
!!$                  + d37(6)*rhof(i,j+2,k) + d37(7)*rhof(i,j+3,k) &
!!$                  + d37(8)*rhof(i,j+4,k) + d37(9)*rhof(i,j+5,k) &
!!$                  + d37(10)*rhof(i,j+6,k)+ d37(11)*rhof(i,j+7,k)
!!$             Krhou(i,j,k) = Krhou(i,j,k) +d37(1)*rhouf(i,j-3,k) &
!!$                  + d37(2)*rhouf(i,j-2,k) +d37(3)*rhouf(i,j-1,k) &
!!$                  + d37(4)*rhouf(i,j,k)   +d37(5)*rhouf(i,j+1,k) &
!!$                  + d37(6)*rhouf(i,j+2,k) +d37(7)*rhouf(i,j+3,k) &
!!$                  + d37(8)*rhouf(i,j+4,k) +d37(9)*rhouf(i,j+5,k) &
!!$                  + d37(10)*rhouf(i,j+6,k)+d37(11)*rhouf(i,j+7,k)
!!$             Krhov(i,j,k) = Krhov(i,j,k) +d37(1)*rhovf(i,j-3,k) &
!!$                  + d37(2)*rhovf(i,j-2,k) +d37(3)*rhovf(i,j-1,k) &
!!$                  + d37(4)*rhovf(i,j,k)   +d37(5)*rhovf(i,j+1,k) &
!!$                  + d37(6)*rhovf(i,j+2,k) +d37(7)*rhovf(i,j+3,k) &
!!$                  + d37(8)*rhovf(i,j+4,k) +d37(9)*rhovf(i,j+5,k) &
!!$                  + d37(10)*rhovf(i,j+6,k)+d37(11)*rhovf(i,j+7,k)
!!$             Krhow(i,j,k) = Krhow(i,j,k) +d37(1)*rhowf(i,j-3,k) &
!!$                  + d37(2)*rhowf(i,j-2,k) +d37(3)*rhowf(i,j-1,k) &
!!$                  + d37(4)*rhowf(i,j,k)   +d37(5)*rhowf(i,j+1,k) &
!!$                  + d37(6)*rhowf(i,j+2,k) +d37(7)*rhowf(i,j+3,k) &
!!$                  + d37(8)*rhowf(i,j+4,k) +d37(9)*rhowf(i,j+5,k) &
!!$                  + d37(10)*rhowf(i,j+6,k)+d37(11)*rhowf(i,j+7,k)
!!$             Krhoe(i,j,k) = Krhoe(i,j,k) +d37(1)*rhoef(i,j-3,k) &
!!$                  + d37(2)*rhoef(i,j-2,k) +d37(3)*rhoef(i,j-1,k) &
!!$                  + d37(4)*rhoef(i,j,k)   +d37(5)*rhoef(i,j+1,k) &
!!$                  + d37(6)*rhoef(i,j+2,k) +d37(7)*rhoef(i,j+3,k) &
!!$                  + d37(8)*rhoef(i,j+4,k) +d37(9)*rhoef(i,j+5,k) &
!!$                  + d37(10)*rhoef(i,j+6,k)+d37(11)*rhoef(i,j+7,k)
!!$          enddo
!!$       enddo
!!$       ! non-centered filtering on stencil 4-6
!!$       j=5
!!$       do k=1,nz
!!$          do i=1,nx
!!$             Krho(i,j,k)  = Krho(i,j,k) + d46(1)*rhof(i,j-4,k) &
!!$                  + d46(2)*rhof(i,j-3,k) + d46(3)*rhof(i,j-2,k) &
!!$                  + d46(4)*rhof(i,j-1,k) + d46(5)*rhof(i,j,k)   &
!!$                  + d46(6)*rhof(i,j+1,k) + d46(7)*rhof(i,j+2,k) &
!!$                  + d46(8)*rhof(i,j+3,k) + d46(9)*rhof(i,j+4,k) &
!!$                  + d46(10)*rhof(i,j+5,k)+ d46(11)*rhof(i,j+6,k)
!!$             Krhou(i,j,k) = Krhou(i,j,k) +d46(1)*rhouf(i,j-4,k) &
!!$                  + d46(2)*rhouf(i,j-3,k) +d46(3)*rhouf(i,j-2,k) &
!!$                  + d46(4)*rhouf(i,j-1,k) +d46(5)*rhouf(i,j,k)   &
!!$                  + d46(6)*rhouf(i,j+1,k) +d46(7)*rhouf(i,j+2,k) &
!!$                  + d46(8)*rhouf(i,j+3,k) +d46(9)*rhouf(i,j+4,k) &
!!$                  + d46(10)*rhouf(i,j+5,k)+d46(11)*rhouf(i,j+6,k)
!!$             Krhov(i,j,k) = Krhov(i,j,k) +d46(1)*rhovf(i,j-4,k) &
!!$                  + d46(2)*rhovf(i,j-3,k) +d46(3)*rhovf(i,j-2,k) &
!!$                  + d46(4)*rhovf(i,j-1,k) +d46(5)*rhovf(i,j,k)   &
!!$                  + d46(6)*rhovf(i,j+1,k) +d46(7)*rhovf(i,j+2,k) &
!!$                  + d46(8)*rhovf(i,j+3,k) +d46(9)*rhovf(i,j+4,k) &
!!$                  + d46(10)*rhovf(i,j+5,k)+d46(11)*rhovf(i,j+6,k)
!!$             Krhow(i,j,k) = Krhow(i,j,k) +d46(1)*rhowf(i,j-4,k) &
!!$                  + d46(2)*rhowf(i,j-3,k) +d46(3)*rhowf(i,j-2,k) &
!!$                  + d46(4)*rhowf(i,j-1,k) +d46(5)*rhowf(i,j,k)   &
!!$                  + d46(6)*rhowf(i,j+1,k) +d46(7)*rhowf(i,j+2,k) &
!!$                  + d46(8)*rhowf(i,j+3,k) +d46(9)*rhowf(i,j+4,k) &
!!$                  + d46(10)*rhowf(i,j+5,k)+d46(11)*rhowf(i,j+6,k)
!!$             Krhoe(i,j,k) = Krhoe(i,j,k) +d46(1)*rhoef(i,j-4,k) &
!!$                  + d46(2)*rhoef(i,j-3,k) +d46(3)*rhoef(i,j-2,k) &
!!$                  + d46(4)*rhoef(i,j-1,k) +d46(5)*rhoef(i,j,k)   &
!!$                  + d46(6)*rhoef(i,j+1,k) +d46(7)*rhoef(i,j+2,k) &
!!$                  + d46(8)*rhoef(i,j+3,k) +d46(9)*rhoef(i,j+4,k) &
!!$                  + d46(10)*rhoef(i,j+5,k)+d46(11)*rhoef(i,j+6,k)
!!$          enddo
!!$       enddo

       ! 2-th order centered filter
       !j=2
       !do k=1,nz
       !   do i=1,nx
       !      Krho(i,j,k)  = Krho(i,j,k)  +d3(0)*rhof(i,j,k)   &
       !           + d3(1)*rhof(i,j+1,k)+d3(1)*rhof(i,j-1,k)
       !      Krhou(i,j,k) = Krhou(i,j,k)  +d3(0)*rhouf(i,j,k)   &
       !           + d3(1)*rhouf(i,j+1,k)+d3(1)*rhouf(i,j-1,k)
       !      Krhov(i,j,k) = Krhov(i,j,k)  +d3(0)*rhovf(i,j,k)   &
       !           + d3(1)*rhovf(i,j+1,k)+d3(1)*rhovf(i,j-1,k)
       !      Krhow(i,j,k) = Krhow(i,j,k)  +d3(0)*rhowf(i,j,k)   &
       !           + d3(1)*rhowf(i,j+1,k)+d3(1)*rhowf(i,j-1,k)
       !      Krhoe(i,j,k) = Krhoe(i,j,k)  +d3(0)*rhoef(i,j,k)   &
       !           + d3(1)*rhoef(i,j+1,k)+d3(1)*rhoef(i,j-1,k)
       !   enddo
       !enddo
       ! 4-th order centered filter
       j=3
       do k=1,nz
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  +d5(0)*rhof(i,j,k)   &
                  + d5(1)*rhof(i,j+1,k)+d5(2)*rhof(i,j+2,k) &
                  + d5(1)*rhof(i,j-1,k)+d5(2)*rhof(i,j-2,k)
             Krhou(i,j,k) = Krhou(i,j,k)  +d5(0)*rhouf(i,j,k)   &
                  + d5(1)*rhouf(i,j+1,k)+d5(2)*rhouf(i,j+2,k) &
                  + d5(1)*rhouf(i,j-1,k)+d5(2)*rhouf(i,j-2,k)
             Krhov(i,j,k) = Krhov(i,j,k)  +d5(0)*rhovf(i,j,k)   &
                  + d5(1)*rhovf(i,j+1,k)+d5(2)*rhovf(i,j+2,k) &
                  + d5(1)*rhovf(i,j-1,k)+d5(2)*rhovf(i,j-2,k)
             Krhow(i,j,k) = Krhow(i,j,k)  +d5(0)*rhowf(i,j,k)   &
                  + d5(1)*rhowf(i,j+1,k)+d5(2)*rhowf(i,j+2,k) &
                  + d5(1)*rhowf(i,j-1,k)+d5(2)*rhowf(i,j-2,k)
             Krhoe(i,j,k) = Krhoe(i,j,k)  +d5(0)*rhoef(i,j,k)   &
                  + d5(1)*rhoef(i,j+1,k)+d5(2)*rhoef(i,j+2,k) &
                  + d5(1)*rhoef(i,j-1,k)+d5(2)*rhoef(i,j-2,k)
          enddo
       enddo
       ! 6-th order centered filter
       j=4
       do k=1,nz
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  +d7(0)*rhof(i,j,k)   &
                  + d7(1)*rhof(i,j+1,k)+d7(2)*rhof(i,j+2,k)+d7(3)*rhof(i,j+3,k) &
                  + d7(1)*rhof(i,j-1,k)+d7(2)*rhof(i,j-2,k)+d7(3)*rhof(i,j-3,k)
             Krhou(i,j,k) = Krhou(i,j,k)  +d7(0)*rhouf(i,j,k)   &
                  + d7(1)*rhouf(i,j+1,k)+d7(2)*rhouf(i,j+2,k)+d7(3)*rhouf(i,j+3,k) &
                  + d7(1)*rhouf(i,j-1,k)+d7(2)*rhouf(i,j-2,k)+d7(3)*rhouf(i,j-3,k)
             Krhov(i,j,k) = Krhov(i,j,k)  +d7(0)*rhovf(i,j,k)   &
                  + d7(1)*rhovf(i,j+1,k)+d7(2)*rhovf(i,j+2,k)+d7(3)*rhovf(i,j+3,k) &
                  + d7(1)*rhovf(i,j-1,k)+d7(2)*rhovf(i,j-2,k)+d7(3)*rhovf(i,j-3,k)
             Krhow(i,j,k) = Krhow(i,j,k)  +d7(0)*rhowf(i,j,k)   &
                  + d7(1)*rhowf(i,j+1,k)+d7(2)*rhowf(i,j+2,k)+d7(3)*rhowf(i,j+3,k) &
                  + d7(1)*rhowf(i,j-1,k)+d7(2)*rhowf(i,j-2,k)+d7(3)*rhowf(i,j-3,k)
             Krhoe(i,j,k) = Krhoe(i,j,k)  +d7(0)*rhoef(i,j,k)   &
                  + d7(1)*rhoef(i,j+1,k)+d7(2)*rhoef(i,j+2,k)+d7(3)*rhoef(i,j+3,k) &
                  + d7(1)*rhoef(i,j-1,k)+d7(2)*rhoef(i,j-2,k)+d7(3)*rhoef(i,j-3,k)
          enddo
       enddo
       ! 8-th order centered filter
       j=5
       do k=1,nz
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  +d9(0)*rhof(i,j,k)   &
                  + d9(1)*rhof(i,j+1,k)+d9(2)*rhof(i,j+2,k)+d9(3)*rhof(i,j+3,k)+d9(4)*rhof(i,j+4,k) &
                  + d9(1)*rhof(i,j-1,k)+d9(2)*rhof(i,j-2,k)+d9(3)*rhof(i,j-3,k)+d9(4)*rhof(i,j-4,k)
             Krhou(i,j,k) = Krhou(i,j,k)  +d9(0)*rhouf(i,j,k)   &
                  + d9(1)*rhouf(i,j+1,k)+d9(2)*rhouf(i,j+2,k)+d9(3)*rhouf(i,j+3,k)+d9(4)*rhouf(i,j+4,k) &
                  + d9(1)*rhouf(i,j-1,k)+d9(2)*rhouf(i,j-2,k)+d9(3)*rhouf(i,j-3,k)+d9(4)*rhouf(i,j-4,k)
             Krhov(i,j,k) = Krhov(i,j,k)  +d9(0)*rhovf(i,j,k)   &
                  + d9(1)*rhovf(i,j+1,k)+d9(2)*rhovf(i,j+2,k)+d9(3)*rhovf(i,j+3,k)+d9(4)*rhovf(i,j+4,k) &
                  + d9(1)*rhovf(i,j-1,k)+d9(2)*rhovf(i,j-2,k)+d9(3)*rhovf(i,j-3,k)+d9(4)*rhovf(i,j-4,k)
             Krhow(i,j,k) = Krhow(i,j,k)  +d9(0)*rhowf(i,j,k)   &
                  + d9(1)*rhowf(i,j+1,k)+d9(2)*rhowf(i,j+2,k)+d9(3)*rhowf(i,j+3,k)+d9(4)*rhowf(i,j+4,k) &
                  + d9(1)*rhowf(i,j-1,k)+d9(2)*rhowf(i,j-2,k)+d9(3)*rhowf(i,j-3,k)+d9(4)*rhowf(i,j-4,k)
             Krhoe(i,j,k) = Krhoe(i,j,k)  +d9(0)*rhoef(i,j,k)   &
                  + d9(1)*rhoef(i,j+1,k)+d9(2)*rhoef(i,j+2,k)+d9(3)*rhoef(i,j+3,k)+d9(4)*rhoef(i,j+4,k) &
                  + d9(1)*rhoef(i,j-1,k)+d9(2)*rhoef(i,j-2,k)+d9(3)*rhoef(i,j-3,k)+d9(4)*rhoef(i,j-4,k)
          enddo
       enddo
    endif

    ! Interior points
    ! ---------------
    ! 10-th order centered filter
    do k=1,nz
       do j=ndy,nfy
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  + d11(0)*rhof(i,j,k) &
                  + d11(1)*( rhof(i,j+1,k)+rhof(i,j-1,k) )&
                  + d11(2)*( rhof(i,j+2,k)+rhof(i,j-2,k) )&
                  + d11(3)*( rhof(i,j+3,k)+rhof(i,j-3,k) )&
                  + d11(4)*( rhof(i,j+4,k)+rhof(i,j-4,k) )&
                  + d11(5)*( rhof(i,j+5,k)+rhof(i,j-5,k) )

             Krhou(i,j,k) = Krhou(i,j,k)  + d11(0)*rhouf(i,j,k) &
                  + d11(1)*( rhouf(i,j+1,k)+rhouf(i,j-1,k) )&
                  + d11(2)*( rhouf(i,j+2,k)+rhouf(i,j-2,k) )&
                  + d11(3)*( rhouf(i,j+3,k)+rhouf(i,j-3,k) )&
                  + d11(4)*( rhouf(i,j+4,k)+rhouf(i,j-4,k) )&
                  + d11(5)*( rhouf(i,j+5,k)+rhouf(i,j-5,k) )

             Krhov(i,j,k) = Krhov(i,j,k)  + d11(0)*rhovf(i,j,k) &
                  + d11(1)*( rhovf(i,j+1,k)+rhovf(i,j-1,k) )&
                  + d11(2)*( rhovf(i,j+2,k)+rhovf(i,j-2,k) )&
                  + d11(3)*( rhovf(i,j+3,k)+rhovf(i,j-3,k) )&
                  + d11(4)*( rhovf(i,j+4,k)+rhovf(i,j-4,k) )&
                  + d11(5)*( rhovf(i,j+5,k)+rhovf(i,j-5,k) )

             Krhow(i,j,k) = Krhow(i,j,k)  + d11(0)*rhowf(i,j,k) &
                  + d11(1)*( rhowf(i,j+1,k)+rhowf(i,j-1,k) )&
                  + d11(2)*( rhowf(i,j+2,k)+rhowf(i,j-2,k) )&
                  + d11(3)*( rhowf(i,j+3,k)+rhowf(i,j-3,k) )&
                  + d11(4)*( rhowf(i,j+4,k)+rhowf(i,j-4,k) )&
                  + d11(5)*( rhowf(i,j+5,k)+rhowf(i,j-5,k) )

             Krhoe(i,j,k) = Krhoe(i,j,k)  + d11(0)*rhoef(i,j,k) &
                  + d11(1)*( rhoef(i,j+1,k)+rhoef(i,j-1,k) )&
                  + d11(2)*( rhoef(i,j+2,k)+rhoef(i,j-2,k) )&
                  + d11(3)*( rhoef(i,j+3,k)+rhoef(i,j-3,k) )&
                  + d11(4)*( rhoef(i,j+4,k)+rhoef(i,j-4,k) )&
                  + d11(5)*( rhoef(i,j+5,k)+rhoef(i,j-5,k) )
          enddo
       enddo
    enddo

    ! BC at jmax
    ! ----------
    if (is_boundary(2,2)) then
       ! 8-th order centered filter
       j=ny-4
       do k=1,nz
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  +d9(0)*rhof(i,j,k)   &
                  + d9(1)*rhof(i,j+1,k)+d9(2)*rhof(i,j+2,k)+d9(3)*rhof(i,j+3,k)+d9(4)*rhof(i,j+4,k) &
                  + d9(1)*rhof(i,j-1,k)+d9(2)*rhof(i,j-2,k)+d9(3)*rhof(i,j-3,k)+d9(4)*rhof(i,j-4,k)
             Krhou(i,j,k) = Krhou(i,j,k)  +d9(0)*rhouf(i,j,k)   &
                  + d9(1)*rhouf(i,j+1,k)+d9(2)*rhouf(i,j+2,k)+d9(3)*rhouf(i,j+3,k)+d9(4)*rhouf(i,j+4,k) &
                  + d9(1)*rhouf(i,j-1,k)+d9(2)*rhouf(i,j-2,k)+d9(3)*rhouf(i,j-3,k)+d9(4)*rhouf(i,j-4,k)
             Krhov(i,j,k) = Krhov(i,j,k)  +d9(0)*rhovf(i,j,k)   &
                  + d9(1)*rhovf(i,j+1,k)+d9(2)*rhovf(i,j+2,k)+d9(3)*rhovf(i,j+3,k)+d9(4)*rhovf(i,j+4,k) &
                  + d9(1)*rhovf(i,j-1,k)+d9(2)*rhovf(i,j-2,k)+d9(3)*rhovf(i,j-3,k)+d9(4)*rhovf(i,j-4,k)
             Krhow(i,j,k) = Krhow(i,j,k)  +d9(0)*rhowf(i,j,k)   &
                  + d9(1)*rhowf(i,j+1,k)+d9(2)*rhowf(i,j+2,k)+d9(3)*rhowf(i,j+3,k)+d9(4)*rhowf(i,j+4,k) &
                  + d9(1)*rhowf(i,j-1,k)+d9(2)*rhowf(i,j-2,k)+d9(3)*rhowf(i,j-3,k)+d9(4)*rhowf(i,j-4,k)
             Krhoe(i,j,k) = Krhoe(i,j,k)  +d9(0)*rhoef(i,j,k)   &
                  + d9(1)*rhoef(i,j+1,k)+d9(2)*rhoef(i,j+2,k)+d9(3)*rhoef(i,j+3,k)+d9(4)*rhoef(i,j+4,k) &
                  + d9(1)*rhoef(i,j-1,k)+d9(2)*rhoef(i,j-2,k)+d9(3)*rhoef(i,j-3,k)+d9(4)*rhoef(i,j-4,k)
          enddo
       enddo
       ! 6-th order centered filter
       j=ny-3
       do k=1,nz
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  +d7(0)*rhof(i,j,k)   &
                  + d7(1)*rhof(i,j+1,k)+d7(2)*rhof(i,j+2,k)+d7(3)*rhof(i,j+3,k) &
                  + d7(1)*rhof(i,j-1,k)+d7(2)*rhof(i,j-2,k)+d7(3)*rhof(i,j-3,k)
             Krhou(i,j,k) = Krhou(i,j,k)  +d7(0)*rhouf(i,j,k)   &
                  + d7(1)*rhouf(i,j+1,k)+d7(2)*rhouf(i,j+2,k)+d7(3)*rhouf(i,j+3,k) &
                  + d7(1)*rhouf(i,j-1,k)+d7(2)*rhouf(i,j-2,k)+d7(3)*rhouf(i,j-3,k)
             Krhov(i,j,k) = Krhov(i,j,k)  +d7(0)*rhovf(i,j,k)   &
                  + d7(1)*rhovf(i,j+1,k)+d7(2)*rhovf(i,j+2,k)+d7(3)*rhovf(i,j+3,k) &
                  + d7(1)*rhovf(i,j-1,k)+d7(2)*rhovf(i,j-2,k)+d7(3)*rhovf(i,j-3,k)
             Krhow(i,j,k) = Krhow(i,j,k)  +d7(0)*rhowf(i,j,k)   &
                  + d7(1)*rhowf(i,j+1,k)+d7(2)*rhowf(i,j+2,k)+d7(3)*rhowf(i,j+3,k) &
                  + d7(1)*rhowf(i,j-1,k)+d7(2)*rhowf(i,j-2,k)+d7(3)*rhowf(i,j-3,k)
             Krhoe(i,j,k) = Krhoe(i,j,k)  +d7(0)*rhoef(i,j,k)   &
                  + d7(1)*rhoef(i,j+1,k)+d7(2)*rhoef(i,j+2,k)+d7(3)*rhoef(i,j+3,k) &
                  + d7(1)*rhoef(i,j-1,k)+d7(2)*rhoef(i,j-2,k)+d7(3)*rhoef(i,j-3,k)
          enddo
       enddo
       ! 4-th order centered filter
       j=ny-2
       do k=1,nz
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  +d5(0)*rhof(i,j,k)   &
                  + d5(1)*rhof(i,j+1,k)+d5(2)*rhof(i,j+2,k) &
                  + d5(1)*rhof(i,j-1,k)+d5(2)*rhof(i,j-2,k)
             Krhou(i,j,k) = Krhou(i,j,k)  +d5(0)*rhouf(i,j,k)   &
                  + d5(1)*rhouf(i,j+1,k)+d5(2)*rhouf(i,j+2,k) &
                  + d5(1)*rhouf(i,j-1,k)+d5(2)*rhouf(i,j-2,k)
             Krhov(i,j,k) = Krhov(i,j,k)  +d5(0)*rhovf(i,j,k)   &
                  + d5(1)*rhovf(i,j+1,k)+d5(2)*rhovf(i,j+2,k) &
                  + d5(1)*rhovf(i,j-1,k)+d5(2)*rhovf(i,j-2,k)
             Krhow(i,j,k) = Krhow(i,j,k)  +d5(0)*rhowf(i,j,k)   &
                  + d5(1)*rhowf(i,j+1,k)+d5(2)*rhowf(i,j+2,k) &
                  + d5(1)*rhowf(i,j-1,k)+d5(2)*rhowf(i,j-2,k)
             Krhoe(i,j,k) = Krhoe(i,j,k)  +d5(0)*rhoef(i,j,k)   &
                  + d5(1)*rhoef(i,j+1,k)+d5(2)*rhoef(i,j+2,k) &
                  + d5(1)*rhoef(i,j-1,k)+d5(2)*rhoef(i,j-2,k)
          enddo
       enddo
!!$       ! 2-th order centered filter
!!$       j=ny-1
!!$       do k=1,nz
!!$          do i=1,nx
!!$             Krho(i,j,k)  = Krho(i,j,k)  +d3(0)*rhof(i,j,k)   &
!!$                  + d3(1)*rhof(i,j+1,k)+d3(1)*rhof(i,j-1,k)
!!$             Krhou(i,j,k) = Krhou(i,j,k)  +d3(0)*rhouf(i,j,k)   &
!!$                  + d3(1)*rhouf(i,j+1,k)+d3(1)*rhouf(i,j-1,k)
!!$             Krhov(i,j,k) = Krhov(i,j,k)  +d3(0)*rhovf(i,j,k)   &
!!$                  + d3(1)*rhovf(i,j+1,k)+d3(1)*rhovf(i,j-1,k)
!!$             Krhow(i,j,k) = Krhow(i,j,k)  +d3(0)*rhowf(i,j,k)   &
!!$                  + d3(1)*rhowf(i,j+1,k)+d3(1)*rhowf(i,j-1,k)
!!$             Krhoe(i,j,k) = Krhoe(i,j,k)  +d3(0)*rhoef(i,j,k)   &
!!$                  + d3(1)*rhoef(i,j+1,k)+d3(1)*rhoef(i,j-1,k)
!!$          enddo
!!$       enddo
    endif

    !****************
    if (is_2D) return
    !****************

    ! Filtering along z
    ! =================
    ! BC at kmin
    ! ----------
    if (is_boundary(3,1)) then
!!$       ! 2-th order centered filter
!!$       k=2
!!$       do j=1,ny
!!$          do i=1,nx
!!$             Krho(i,j,k)  = Krho(i,j,k)  +d3(0)*rhof(i,j,k)   &
!!$                  + d3(1)*rhof(i,j,k+1)+d3(1)*rhof(i,j,k-1)
!!$             Krhou(i,j,k) = Krhou(i,j,k)  +d3(0)*rhouf(i,j,k)   &
!!$                  + d3(1)*rhouf(i,j,k+1)+d3(1)*rhouf(i,j,k-1)
!!$             Krhov(i,j,k) = Krhov(i,j,k)  +d3(0)*rhovf(i,j,k)   &
!!$                  + d3(1)*rhovf(i,j,k+1)+d3(1)*rhovf(i,j,k-1)
!!$             Krhow(i,j,k) = Krhow(i,j,k)  +d3(0)*rhowf(i,j,k)   &
!!$                  + d3(1)*rhowf(i,j,k+1)+d3(1)*rhowf(i,j,k-1)
!!$             Krhoe(i,j,k) = Krhoe(i,j,k)  +d3(0)*rhoef(i,j,k)   &
!!$                  + d3(1)*rhoef(i,j,k+1)+d3(1)*rhoef(i,j,k-1)
!!$          enddo
!!$       enddo
!!$       ! 4-th order centered filter
!!$       k=3
!!$       do j=1,ny
!!$          do i=1,nx
!!$             Krho(i,j,k)  = Krho(i,j,k)  +d5(0)*rhof(i,j,k)   &
!!$                  + d5(1)*rhof(i,j,k+1)+d5(2)*rhof(i,j,k+2) &
!!$                  + d5(1)*rhof(i,j,k-1)+d5(2)*rhof(i,j,k-2)
!!$             Krhou(i,j,k) = Krhou(i,j,k)  +d5(0)*rhouf(i,j,k)   &
!!$                  + d5(1)*rhouf(i,j,k+1)+d5(2)*rhouf(i,j,k+2) &
!!$                  + d5(1)*rhouf(i,j,k-1)+d5(2)*rhouf(i,j,k-2)
!!$             Krhov(i,j,k) = Krhov(i,j,k)  +d5(0)*rhovf(i,j,k)   &
!!$                  + d5(1)*rhovf(i,j,k+1)+d5(2)*rhovf(i,j,k+2) &
!!$                  + d5(1)*rhovf(i,j,k-1)+d5(2)*rhovf(i,j,k-2)
!!$             Krhow(i,j,k) = Krhow(i,j,k)  +d5(0)*rhowf(i,j,k)   &
!!$                  + d5(1)*rhowf(i,j,k+1)+d5(2)*rhowf(i,j,k+2) &
!!$                  + d5(1)*rhowf(i,j,k-1)+d5(2)*rhowf(i,j,k-2)
!!$             Krhoe(i,j,k) = Krhoe(i,j,k)  +d5(0)*rhoef(i,j,k)   &
!!$                  + d5(1)*rhoef(i,j,k+1)+d5(2)*rhoef(i,j,k+2) &
!!$                  + d5(1)*rhoef(i,j,k-1)+d5(2)*rhoef(i,j,k-2)
!!$          enddo
!!$       enddo
       ! 6-th order centered filter
       k=4
       do j=1,ny
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  +d7(0)*rhof(i,j,k)   &
                  + d7(1)*rhof(i,j,k+1)+d7(2)*rhof(i,j,k+2)+d7(3)*rhof(i,j,k+3) &
                  + d7(1)*rhof(i,j,k-1)+d7(2)*rhof(i,j,k-2)+d7(3)*rhof(i,j,k-3)
             Krhou(i,j,k) = Krhou(i,j,k)  +d7(0)*rhouf(i,j,k)   &
                  + d7(1)*rhouf(i,j,k+1)+d7(2)*rhouf(i,j,k+2)+d7(3)*rhouf(i,j,k+3) &
                  + d7(1)*rhouf(i,j,k-1)+d7(2)*rhouf(i,j,k-2)+d7(3)*rhouf(i,j,k-3)
             Krhov(i,j,k) = Krhov(i,j,k)  +d7(0)*rhovf(i,j,k)   &
                  + d7(1)*rhovf(i,j,k+1)+d7(2)*rhovf(i,j,k+2)+d7(3)*rhovf(i,j,k+3) &
                  + d7(1)*rhovf(i,j,k-1)+d7(2)*rhovf(i,j,k-2)+d7(3)*rhovf(i,j,k-3)
             Krhow(i,j,k) = Krhow(i,j,k)  +d7(0)*rhowf(i,j,k)   &
                  + d7(1)*rhowf(i,j,k+1)+d7(2)*rhowf(i,j,k+2)+d7(3)*rhowf(i,j,k+3) &
                  + d7(1)*rhowf(i,j,k-1)+d7(2)*rhowf(i,j,k-2)+d7(3)*rhowf(i,j,k-3)
             Krhoe(i,j,k) = Krhoe(i,j,k)  +d7(0)*rhoef(i,j,k)   &
                  + d7(1)*rhoef(i,j,k+1)+d7(2)*rhoef(i,j,k+2)+d7(3)*rhoef(i,j,k+3) &
                  + d7(1)*rhoef(i,j,k-1)+d7(2)*rhoef(i,j,k-2)+d7(3)*rhoef(i,j,k-3)
          enddo
       enddo
       ! 8-th order centered filter
       k=5
       do j=1,ny
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  +d9(0)*rhof(i,j,k)   &
                  + d9(1)*rhof(i,j,k+1)+d9(2)*rhof(i,j,k+2)+d9(3)*rhof(i,j,k+3)+d9(4)*rhof(i,j,k+4) &
                  + d9(1)*rhof(i,j,k-1)+d9(2)*rhof(i,j,k-2)+d9(3)*rhof(i,j,k-3)+d9(4)*rhof(i,j,k-4)
             Krhou(i,j,k) = Krhou(i,j,k)  +d9(0)*rhouf(i,j,k)   &
                  + d9(1)*rhouf(i,j,k+1)+d9(2)*rhouf(i,j,k+2)+d9(3)*rhouf(i,j,k+3)+d9(4)*rhouf(i,j,k+4) &
                  + d9(1)*rhouf(i,j,k-1)+d9(2)*rhouf(i,j,k-2)+d9(3)*rhouf(i,j,k-3)+d9(4)*rhouf(i,j,k-4)
             Krhov(i,j,k) = Krhov(i,j,k)  +d9(0)*rhovf(i,j,k)   &
                  + d9(1)*rhovf(i,j,k+1)+d9(2)*rhovf(i,j,k+2)+d9(3)*rhovf(i,j,k+3)+d9(4)*rhovf(i,j,k+4) &
                  + d9(1)*rhovf(i,j,k-1)+d9(2)*rhovf(i,j,k-2)+d9(3)*rhovf(i,j,k-3)+d9(4)*rhovf(i,j,k-4)
             Krhow(i,j,k) = Krhow(i,j,k)  +d9(0)*rhowf(i,j,k)   &
                  + d9(1)*rhowf(i,j,k+1)+d9(2)*rhowf(i,j,k+2)+d9(3)*rhowf(i,j,k+3)+d9(4)*rhowf(i,j,k+4) &
                  + d9(1)*rhowf(i,j,k-1)+d9(2)*rhowf(i,j,k-2)+d9(3)*rhowf(i,j,k-3)+d9(4)*rhowf(i,j,k-4)
             Krhoe(i,j,k) = Krhoe(i,j,k)  +d9(0)*rhoef(i,j,k)   &
                  + d9(1)*rhoef(i,j,k+1)+d9(2)*rhoef(i,j,k+2)+d9(3)*rhoef(i,j,k+3)+d9(4)*rhoef(i,j,k+4) &
                  + d9(1)*rhoef(i,j,k-1)+d9(2)*rhoef(i,j,k-2)+d9(3)*rhoef(i,j,k-3)+d9(4)*rhoef(i,j,k-4)
          enddo
       enddo
    endif

    ! Interior points
    ! ---------------
    do k=ndz,nfz
       do j=1,ny
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  + d11(0)*rhof(i,j,k) &
                  + d11(1)*( rhof(i,j,k+1)+rhof(i,j,k-1) )&
                  + d11(2)*( rhof(i,j,k+2)+rhof(i,j,k-2) )&
                  + d11(3)*( rhof(i,j,k+3)+rhof(i,j,k-3) )&
                  + d11(4)*( rhof(i,j,k+4)+rhof(i,j,k-4) )&
                  + d11(5)*( rhof(i,j,k+5)+rhof(i,j,k-5) )

             Krhou(i,j,k) = Krhou(i,j,k)  + d11(0)*rhouf(i,j,k) &
                  + d11(1)*( rhouf(i,j,k+1)+rhouf(i,j,k-1) )&
                  + d11(2)*( rhouf(i,j,k+2)+rhouf(i,j,k-2) )&
                  + d11(3)*( rhouf(i,j,k+3)+rhouf(i,j,k-3) )&
                  + d11(4)*( rhouf(i,j,k+4)+rhouf(i,j,k-4) )&
                  + d11(5)*( rhouf(i,j,k+5)+rhouf(i,j,k-5) )

             Krhov(i,j,k) = Krhov(i,j,k)  + d11(0)*rhovf(i,j,k) &
                  + d11(1)*( rhovf(i,j,k+1)+rhovf(i,j,k-1) )&
                  + d11(2)*( rhovf(i,j,k+2)+rhovf(i,j,k-2) )&
                  + d11(3)*( rhovf(i,j,k+3)+rhovf(i,j,k-3) )&
                  + d11(4)*( rhovf(i,j,k+4)+rhovf(i,j,k-4) )&
                  + d11(5)*( rhovf(i,j,k+5)+rhovf(i,j,k-5) )

             Krhow(i,j,k) = Krhow(i,j,k)  + d11(0)*rhowf(i,j,k) &
                  + d11(1)*( rhowf(i,j,k+1)+rhowf(i,j,k-1) )&
                  + d11(2)*( rhowf(i,j,k+2)+rhowf(i,j,k-2) )&
                  + d11(3)*( rhowf(i,j,k+3)+rhowf(i,j,k-3) )&
                  + d11(4)*( rhowf(i,j,k+4)+rhowf(i,j,k-4) )&
                  + d11(5)*( rhowf(i,j,k+5)+rhowf(i,j,k-5) )

             Krhoe(i,j,k) = Krhoe(i,j,k)  + d11(0)*rhoef(i,j,k) &
                  + d11(1)*( rhoef(i,j,k+1)+rhoef(i,j,k-1) )&
                  + d11(2)*( rhoef(i,j,k+2)+rhoef(i,j,k-2) )&
                  + d11(3)*( rhoef(i,j,k+3)+rhoef(i,j,k-3) )&
                  + d11(4)*( rhoef(i,j,k+4)+rhoef(i,j,k-4) )&
                  + d11(5)*( rhoef(i,j,k+5)+rhoef(i,j,k-5) )
          enddo
       enddo
    enddo

    ! BC at kmax
    ! ----------
    if (is_boundary(3,2)) then
       ! 8-th order centered filter
       k=nz-4
       do j=1,ny
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  +d9(0)*rhof(i,j,k)   &
                  + d9(1)*rhof(i,j,k+1)+d9(2)*rhof(i,j,k+2)+d9(3)*rhof(i,j,k+3)+d9(4)*rhof(i,j,k+4) &
                  + d9(1)*rhof(i,j,k-1)+d9(2)*rhof(i,j,k-2)+d9(3)*rhof(i,j,k-3)+d9(4)*rhof(i,j,k-4)
             Krhou(i,j,k) = Krhou(i,j,k)  +d9(0)*rhouf(i,j,k)   &
                  + d9(1)*rhouf(i,j,k+1)+d9(2)*rhouf(i,j,k+2)+d9(3)*rhouf(i,j,k+3)+d9(4)*rhouf(i,j,k+4) &
                  + d9(1)*rhouf(i,j,k-1)+d9(2)*rhouf(i,j,k-2)+d9(3)*rhouf(i,j,k-3)+d9(4)*rhouf(i,j,k-4)
             Krhov(i,j,k) = Krhov(i,j,k)  +d9(0)*rhovf(i,j,k)   &
                  + d9(1)*rhovf(i,j,k+1)+d9(2)*rhovf(i,j,k+2)+d9(3)*rhovf(i,j,k+3)+d9(4)*rhovf(i,j,k+4) &
                  + d9(1)*rhovf(i,j,k-1)+d9(2)*rhovf(i,j,k-2)+d9(3)*rhovf(i,j,k-3)+d9(4)*rhovf(i,j,k-4)
             Krhow(i,j,k) = Krhow(i,j,k)  +d9(0)*rhowf(i,j,k)   &
                  + d9(1)*rhowf(i,j,k+1)+d9(2)*rhowf(i,j,k+2)+d9(3)*rhowf(i,j,k+3)+d9(4)*rhowf(i,j,k+4) &
                  + d9(1)*rhowf(i,j,k-1)+d9(2)*rhowf(i,j,k-2)+d9(3)*rhowf(i,j,k-3)+d9(4)*rhowf(i,j,k-4)
             Krhoe(i,j,k) = Krhoe(i,j,k)  +d9(0)*rhoef(i,j,k)   &
                  + d9(1)*rhoef(i,j,k+1)+d9(2)*rhoef(i,j,k+2)+d9(3)*rhoef(i,j,k+3)+d9(4)*rhoef(i,j,k+4) &
                  + d9(1)*rhoef(i,j,k-1)+d9(2)*rhoef(i,j,k-2)+d9(3)*rhoef(i,j,k-3)+d9(4)*rhoef(i,j,k-4)
          enddo
       enddo
       ! 6-th order centered filter
       k=nz-3
       do j=1,ny
          do i=1,nx
             Krho(i,j,k)  = Krho(i,j,k)  +d7(0)*rhof(i,j,k)   &
                  + d7(1)*rhof(i,j,k+1)+d7(2)*rhof(i,j,k+2)+d7(3)*rhof(i,j,k+3) &
                  + d7(1)*rhof(i,j,k-1)+d7(2)*rhof(i,j,k-2)+d7(3)*rhof(i,j,k-3)
             Krhou(i,j,k) = Krhou(i,j,k)  +d7(0)*rhouf(i,j,k)   &
                  + d7(1)*rhouf(i,j,k+1)+d7(2)*rhouf(i,j,k+2)+d7(3)*rhouf(i,j,k+3) &
                  + d7(1)*rhouf(i,j,k-1)+d7(2)*rhouf(i,j,k-2)+d7(3)*rhouf(i,j,k-3)
             Krhov(i,j,k) = Krhov(i,j,k)  +d7(0)*rhovf(i,j,k)   &
                  + d7(1)*rhovf(i,j,k+1)+d7(2)*rhovf(i,j,k+2)+d7(3)*rhovf(i,j,k+3) &
                  + d7(1)*rhovf(i,j,k-1)+d7(2)*rhovf(i,j,k-2)+d7(3)*rhovf(i,j,k-3)
             Krhow(i,j,k) = Krhow(i,j,k)  +d7(0)*rhowf(i,j,k)   &
                  + d7(1)*rhowf(i,j,k+1)+d7(2)*rhowf(i,j,k+2)+d7(3)*rhowf(i,j,k+3) &
                  + d7(1)*rhowf(i,j,k-1)+d7(2)*rhowf(i,j,k-2)+d7(3)*rhowf(i,j,k-3)
             Krhoe(i,j,k) = Krhoe(i,j,k)  +d7(0)*rhoef(i,j,k)   &
                  + d7(1)*rhoef(i,j,k+1)+d7(2)*rhoef(i,j,k+2)+d7(3)*rhoef(i,j,k+3) &
                  + d7(1)*rhoef(i,j,k-1)+d7(2)*rhoef(i,j,k-2)+d7(3)*rhoef(i,j,k-3)
          enddo
       enddo
!!$       ! 4-th order centered filter
!!$       k=nz-2
!!$       do j=1,ny
!!$          do i=1,nx
!!$             Krho(i,j,k)  = Krho(i,j,k)  +d5(0)*rhof(i,j,k)   &
!!$                  + d5(1)*rhof(i,j,k+1)+d5(2)*rhof(i,j,k+2) &
!!$                  + d5(1)*rhof(i,j,k-1)+d5(2)*rhof(i,j,k-2)
!!$             Krhou(i,j,k) = Krhou(i,j,k)  +d5(0)*rhouf(i,j,k)   &
!!$                  + d5(1)*rhouf(i,j,k+1)+d5(2)*rhouf(i,j,k+2) &
!!$                  + d5(1)*rhouf(i,j,k-1)+d5(2)*rhouf(i,j,k-2)
!!$             Krhov(i,j,k) = Krhov(i,j,k)  +d5(0)*rhovf(i,j,k)   &
!!$                  + d5(1)*rhovf(i,j,k+1)+d5(2)*rhovf(i,j,k+2) &
!!$                  + d5(1)*rhovf(i,j,k-1)+d5(2)*rhovf(i,j,k-2)
!!$             Krhow(i,j,k) = Krhow(i,j,k)  +d5(0)*rhowf(i,j,k)   &
!!$                  + d5(1)*rhowf(i,j,k+1)+d5(2)*rhowf(i,j,k+2) &
!!$                  + d5(1)*rhowf(i,j,k-1)+d5(2)*rhowf(i,j,k-2)
!!$             Krhoe(i,j,k) = Krhoe(i,j,k)  +d5(0)*rhoef(i,j,k)   &
!!$                  + d5(1)*rhoef(i,j,k+1)+d5(2)*rhoef(i,j,k+2) &
!!$                  + d5(1)*rhoef(i,j,k-1)+d5(2)*rhoef(i,j,k-2)
!!$          enddo
!!$       enddo
!!$       ! 2-th order centered filter
!!$       k=nz-1
!!$       do j=1,ny
!!$          do i=1,nx
!!$             Krho(i,j,k)  = Krho(i,j,k)  +d3(0)*rhof(i,j,k)   &
!!$                  + d3(1)*rhof(i,j,k+1)+d3(1)*rhof(i,j,k-1)
!!$             Krhou(i,j,k) = Krhou(i,j,k)  +d3(0)*rhouf(i,j,k)   &
!!$                  + d3(1)*rhouf(i,j,k+1)+d3(1)*rhouf(i,j,k-1)
!!$             Krhov(i,j,k) = Krhov(i,j,k)  +d3(0)*rhovf(i,j,k)   &
!!$                  + d3(1)*rhovf(i,j,k+1)+d3(1)*rhovf(i,j,k-1)
!!$             Krhow(i,j,k) = Krhow(i,j,k)  +d3(0)*rhowf(i,j,k)   &
!!$                  + d3(1)*rhowf(i,j,k+1)+d3(1)*rhowf(i,j,k-1)
!!$             Krhoe(i,j,k) = Krhoe(i,j,k)  +d3(0)*rhoef(i,j,k)   &
!!$                  + d3(1)*rhoef(i,j,k+1)+d3(1)*rhoef(i,j,k-1)
!!$          enddo
!!$       enddo
    endif

  end subroutine filtering_fluct_11pts
