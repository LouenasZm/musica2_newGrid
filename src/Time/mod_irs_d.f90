!=================================================================================
module mod_irs_d
!=================================================================================
  !> author: Aurelien Bienner, modifs XG 
  !> date: 2021
  !> Module for Implicit Residual Smoothing parallelized with ghost cells
  !> * fusion of Cartesian & curvilinear version *
!=================================================================================
  use mod_constant ! <- for is_wall2
  implicit none
  ! ------------------------------------------------------------------------------
  ! ------------------------------------------------------------------------------

contains

  !===============================================================================
  subroutine irs_ngh_d
  !===============================================================================
    !> 4th-order Implicit Residual Smoothing (IRS4) * per direction *
    !> Parallel inversion of pentadiagonal matrix using ghost points
  !===============================================================================
    use mod_interface ! <- procedure pointer interfaces
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    call apply_irs_1
    call apply_irs_2
    call apply_irs_3

  end subroutine irs_ngh_d

  !===============================================================================
  subroutine irs2_ngh_d
  !===============================================================================
    !> 2nd-order Implicit Residual Smoothing (IRS4) * per direction *
    !> Parallel inversion of pentadiagonal matrix using ghost points
  !===============================================================================
    use mod_time ! <- for is_irs
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    if (is_irs_i) call irs2_ngh_i

    if (is_irs_j) call irs2_ngh_j

    if (is_irs_k) call irs2_ngh_k

  end subroutine irs2_ngh_d
   
  !===============================================================================
  subroutine irs4_ngh_d
  !===============================================================================
    !> 4th-order Implicit Residual Smoothing (IRS4) * per direction *
    !> Parallel inversion of pentadiagonal matrix using ghost points
  !===============================================================================
    use mod_time ! <- for is_irs
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    if (is_irs_i) call irs4_ngh_i

    if (is_irs_j) call irs4_ngh_j

    if (is_irs_k) call irs4_ngh_k

  end subroutine irs4_ngh_d

  !===============================================================================
  subroutine irs2_ngh_i
  !===============================================================================
    !> 2nd-order Implicit Residual Smoothing (IRS2) - Implicitation of i-direction
    !> Parallel inversion of tridiagonal matrix using ghost points
  !===============================================================================
    use mod_time   ! <- for deltat,is_irs
    use mod_flow   ! <- for velocities,increments,metrics
    use mod_comm1  ! <- for communication1_inc & CFL
    use mod_interface  ! <- for communication_inc
    implicit none
    ! ----------------------------------------------------------------------------
    ! indices
    integer :: i,j,k
    integer :: i1,i2,i3
    integer :: i1min,i1max,i2min,i2max,i3min,i3max
    ! contravariant velocities & spectral radius
    real(wp) :: vcmh,vcph,rspecmh,rspecph
    ! tridiagonal in i-direction
    ! --------------------------
    ! matrix elements a,d,c & right-hand side (RHS)
    real(wp), dimension(ny,nx1_irs:nx2_irs) :: a_i,d_i,c_i
    real(wp), dimension(ny,nx1_irs:nx2_irs) :: RHS_i,RHSu_i,RHSv_i,RHSw_i,RHSe_i
    real(wp), dimension(ny) :: temp_i ! work array
    ! ----------------------------------------------------------------------------

    !  i3 : i;   i1 : k;   i2 : j
    i1min=1;   i2min=1;   i3min=ndx-ngh_irs(1)
    i1max=nz;  i2max=ny;  i3max=nfx+ngh_irs(2)
    if ((is_boundary(1,1)).and.(i3min<1))  i3min=1
    if ((is_boundary(1,2)).and.(i3max>nx)) i3max=nx

    ! If wall, do not take into account the wall point
    ! ================================================
    if (.not.is_wall2) then
       if (is_bc_wall(2,1)) i2min=2
       if (is_bc_wall(2,2)) i2max=ny-1
       if (is_bc_wall(3,1)) i1min=2
       if (is_bc_wall(3,2)) i1max=nz-1
!!$    else ! for new wall conditions, set some increments to zero
!!$       if (is_bc_wall(1,1)) then
!!$          do i2=i2min,i2max
!!$             Krhou(1,i2,:)=0.0_wp
!!$             Krhov(1,i2,:)=0.0_wp
!!$             Krhow(1,i2,:)=0.0_wp
!!$             Krhoe(1,i2,:)=0.0_wp
!!$          enddo
!!$       endif
!!$       if (is_bc_wall(1,2)) then
!!$          do i2=i2min,i2max
!!$             Krhou(nx,i2,:)=0.0_wp
!!$             Krhov(nx,i2,:)=0.0_wp
!!$             Krhow(nx,i2,:)=0.0_wp
!!$             Krhoe(nx,i2,:)=0.0_wp
!!$          enddo
!!$       endif
!!$       if (is_bc_wall(2,1)) then
!!$          do i3=i3min,i3max
!!$             Krhou(i3,1,:)=0.0_wp
!!$             Krhov(i3,1,:)=0.0_wp
!!$             Krhow(i3,1,:)=0.0_wp
!!$             Krhoe(i3,1,:)=0.0_wp
!!$          enddo
!!$       endif
!!$       if (is_bc_wall(2,2)) then
!!$          do i3=i3min,i3max
!!$             Krhou(i3,ny,:)=0.0_wp
!!$             Krhov(i3,ny,:)=0.0_wp
!!$             Krhow(i3,ny,:)=0.0_wp
!!$             Krhoe(i3,ny,:)=0.0_wp
!!$          enddo
!!$       endif
!!$       if (is_bc_wall(3,1)) then
!!$          do i3=i3min,i3max
!!$             Krhou(i3,:,1)=0.0_wp
!!$             Krhov(i3,:,1)=0.0_wp
!!$             Krhow(i3,:,1)=0.0_wp
!!$             Krhoe(i3,:,1)=0.0_wp
!!$          enddo
!!$       endif
!!$       if (is_bc_wall(3,2)) then
!!$          do i3=i3min,i3max
!!$             Krhou(i3,:,nz)=0.0_wp
!!$             Krhov(i3,:,nz)=0.0_wp
!!$             Krhow(i3,:,nz)=0.0_wp
!!$             Krhoe(i3,:,nz)=0.0_wp
!!$          enddo
!!$       endif
    endif
 
    ! One-sided communications of increments in i-direction
    ! =====================================================
!!$    ! open target windows
!!$    call window_comm1_inc
!!$    ! communicate var
!!$    call communication_1_inc_i
!!$    call communication1_inc_i

    ! Calculation of the spectral radius
    ! ==================================
    if (is_dtlocal) then
       if (is_curv3) then
          do k=1,nz
             do j=1,ny
                do i=1,nx
                   vcmh=sqrt((uu(i-1,j,k)*ksi_x(i-1,j,k)+vv(i-1,j,k)*ksi_y(i-1,j,k)+ww(i-1,j,k)*ksi_z(i-1,j,k))**2 &
                            +(uu(i-1,j,k)*eta_x(i-1,j,k)+vv(i-1,j,k)*eta_y(i-1,j,k)+ww(i-1,j,k)*eta_z(i-1,j,k))**2 &
                            +(uu(i-1,j,k)*phi_x(i-1,j,k)+vv(i-1,j,k)*phi_y(i-1,j,k)+ww(i-1,j,k)*phi_z(i-1,j,k))**2 &
                            )*ijacob3(i-1,j,k)
                   vcph=sqrt((uu(i,j,k)*ksi_x(i,j,k)+vv(i,j,k)*ksi_y(i,j,k)+ww(i,j,k)*ksi_z(i,j,k))**2 &
                            +(uu(i,j,k)*eta_x(i,j,k)+vv(i,j,k)*eta_y(i,j,k)+ww(i,j,k)*eta_z(i,j,k))**2 &
                            +(uu(i,j,k)*phi_x(i,j,k)+vv(i,j,k)*phi_y(i,j,k)+ww(i,j,k)*phi_z(i,j,k))**2 &
                            )*ijacob3(i,j,k)
                   rspecmh=abs(vcmh)+c_(i-1,j,k)*sqrt(g3_ksi(i-1,j,k)**2+g3_eta(i-1,j,k)**2+g3_phi(i-1,j,k)**2)
                   rspecph=abs(vcph)+c_(i  ,j,k)*sqrt(g3_ksi(i  ,j,k)**2+g3_eta(i  ,j,k)**2+g3_phi(i  ,j,k)**2)
                   cfl_l(i,j,k)=0.5_wp*(rspecmh+rspecph)*dt_local(i,j,k)
                enddo
             enddo
          enddo
       else
          if (is_curv) then
             do k=1,nz
                do j=1,ny
                   do i=1,nx
                      vcmh=sqrt((uu(i-1,j,k)*y_eta(i-1,j)-vv(i-1,j,k)*x_eta(i-1,j))**2 &
                               +(vv(i-1,j,k)*x_ksi(i-1,j)-uu(i,j-1,k)*y_ksi(i-1,j))**2)*ijacob(i-1,j)
                      vcph=sqrt((uu(i,j  ,k)*y_eta(i,j  )-vv(i,j  ,k)*x_eta(i,j  ))**2 &
                               +(vv(i,j  ,k)*x_ksi(i,j  )-uu(i,j  ,k)*y_ksi(i,j  ))**2)*ijacob(i,j)
                      rspecmh=abs(vcmh)+c_(i-1,j,k)*sqrt(g_ksi(i-1,j)**2+g_eta(i-1,j)**2)
                      rspecph=abs(vcph)+c_(i  ,j,k)*sqrt(g_ksi(i  ,j)**2+g_eta(i  ,j)**2)
                      cfl_l(i,j,k)=0.5_wp*(rspecmh+rspecph)*dt_local(i,j,k)
                   enddo
                enddo
             enddo
          else
             do k=1,nz
                do j=1,ny
                   do i=1,nx
                      rspecmh=(sqrt(uu(i-1,j,k)**2+vv(i-1,j,k)**2)+c_(i,j-1,k))*idy(j-1)
                      rspecph=(sqrt(uu(i  ,j,k)**2+vv(i  ,j,k)**2)+c_(i,j  ,k))*idy(j)
                      cfl_l(i,j,k)=0.5_wp*(rspecmh+rspecph)*dt_local(i,j,k)
                   enddo
                enddo
             enddo
          endif
       endif
    else
       if (is_curv3) then
          do k=1,nz
             do j=1,ny
                do i=1,nx
                   vcmh=sqrt((uu(i-1,j,k)*ksi_x(i-1,j,k)+vv(i-1,j,k)*ksi_y(i-1,j,k)+ww(i-1,j,k)*ksi_z(i-1,j,k))**2 &
                            +(uu(i-1,j,k)*eta_x(i-1,j,k)+vv(i-1,j,k)*eta_y(i-1,j,k)+ww(i-1,j,k)*eta_z(i-1,j,k))**2 &
                            +(uu(i-1,j,k)*phi_x(i-1,j,k)+vv(i-1,j,k)*phi_y(i-1,j,k)+ww(i-1,j,k)*phi_z(i-1,j,k))**2 &
                            )*ijacob3(i-1,j,k)
                   vcph=sqrt((uu(i,j,k)*ksi_x(i,j,k)+vv(i,j,k)*ksi_y(i,j,k)+ww(i,j,k)*ksi_z(i,j,k))**2 &
                            +(uu(i,j,k)*eta_x(i,j,k)+vv(i,j,k)*eta_y(i,j,k)+ww(i,j,k)*eta_z(i,j,k))**2 &
                            +(uu(i,j,k)*phi_x(i,j,k)+vv(i,j,k)*phi_y(i,j,k)+ww(i,j,k)*phi_z(i,j,k))**2 &
                            )*ijacob3(i,j,k)
                   rspecmh=abs(vcmh)+c_(i-1,j,k)*sqrt(g3_ksi(i-1,j,k)**2+g3_eta(i-1,j,k)**2+g3_phi(i-1,j,k)**2)
                   rspecph=abs(vcph)+c_(i  ,j,k)*sqrt(g3_ksi(i  ,j,k)**2+g3_eta(i  ,j,k)**2+g3_phi(i  ,j,k)**2)
                   cfl_l(i,j,k)=0.5_wp*(rspecmh+rspecph)*deltat
                enddo
             enddo
          enddo
       else
          if (is_curv) then
             do k=1,nz
                do j=1,ny
                   do i=1,nx
                      vcmh=sqrt((uu(i-1,j,k)*y_eta(i-1,j)-vv(i-1,j,k)*x_eta(i-1,j))**2 &
                               +(vv(i-1,j,k)*x_ksi(i-1,j)-uu(i,j-1,k)*y_ksi(i-1,j))**2)*ijacob(i-1,j)
                      vcph=sqrt((uu(i,j  ,k)*y_eta(i,j  )-vv(i,j  ,k)*x_eta(i,j  ))**2 &
                               +(vv(i,j  ,k)*x_ksi(i,j  )-uu(i,j  ,k)*y_ksi(i,j  ))**2)*ijacob(i,j)
                      rspecmh=abs(vcmh)+c_(i-1,j,k)*sqrt(g_ksi(i-1,j)**2+g_eta(i-1,j)**2)
                      rspecph=abs(vcph)+c_(i  ,j,k)*sqrt(g_ksi(i  ,j)**2+g_eta(i  ,j)**2)
                      cfl_l(i,j,k)=0.5_wp*(rspecmh+rspecph)*deltat
                   enddo
                enddo
             enddo
          else
             do k=1,nz
                do j=1,ny
                   do i=1,nx
                      rspecmh=(abs(uu(i-1,j,k))+c_(i-1,j,k))*idx(i-1)
                      rspecph=(abs(uu(i  ,j,k))+c_(i  ,j,k))*idx(i)
                      cfl_l(i,j,k)=0.5_wp*(rspecmh+rspecph)*deltat
                   enddo
                enddo
             enddo
          endif
       endif
    endif

!!$    ! Check that communications are fullfilled
!!$    ! ----------------------------------------
!!$    call window_comm1_inc

!!$    ! One-sided communications of CFL in i-direction
!!$    ! ==============================================
!!$    ! open target windows
!!$    call window_comm1_cfl
!!$    ! communicate var
!!$    call communication_1_cfl_i
!!$
!!$    ! Check that communications are fullfilled
!!$    ! ----------------------------------------
!!$    call window_comm1_cfl

    call communication_inc(Krho,Krhou,Krhov,Krhow,Krhoe,cfl_l)

    do i1=i1min,i1max

       ! Fill tridiagonal matrix elements & RHS
       ! ======================================
       ! Construction of the tridiagonal matrix
       ! The linear systems is Mx=y with :
       !     _M=(m_ij) for i,j in {1,..,n} and with :
       !                 m_ii=d_i for i in {1,..,n}
       !                 m_ii-1=c_i-1 for i in {2,..,n}
       !                 m_ii+1=a_i for i in {1,..,n-1}
       !     _ y=RHS=dwi
       !     _d and y are of dimensions n, a and c of dimensions n-1,
       i3=i3min
       do i2=i2min,i2max
          a_i(i2,i3)=-theta_irs1*cfl_l(i3+1,i2,i1)
          d_i(i2,i3)=1.0_wp-a_i(i2,i3)
       enddo

       i3=i3max
       do i2=i2min,i2max
          c_i(i2,i3-1)=-theta_irs1*cfl_l(i3,i2,i1)
          d_i(i2,i3)  =1.0_wp-c_i(i2,i3-1)
       enddo

       do i3=i3min+1,i3max-1
          do i2=i2min, i2max
             c_i(i2,i3-1)=-theta_irs2*cfl_l(i3,i2,i1)**2
             a_i(i2,i3)  =-theta_irs2*cfl_l(i3+1,i2,i1)**2
             d_i(i2,i3)  =1.0_wp-(c_i(i2,i3-1)+a_i(i2,i3))
          enddo
       enddo

       ! Filling of the RHS
       do i3=i3min,i3max
          do i2=i2min,i2max
             RHS_i(i2,i3)= Krho(i3,i2,i1)
             RHSu_i(i2,i3)=Krhou(i3,i2,i1)
             RHSv_i(i2,i3)=Krhov(i3,i2,i1)
             RHSw_i(i2,i3)=Krhow(i3,i2,i1)
             RHSe_i(i2,i3)=Krhoe(i3,i2,i1)
          enddo
       enddo

       ! Resolution of the tridiagonal linear system
       ! ===========================================
       ! Thomas' algorithm: Forward elimination phase
       do i2=i2min,i2max
          temp_i(i2)=1.0_wp/d_i(i2,i3min)
          a_i(i2,i3min)=a_i(i2,i3min)*temp_i(i2)
          RHS_i(i2,i3min)= RHS_i(i2,i3min)*temp_i(i2)
          RHSu_i(i2,i3min)=RHSu_i(i2,i3min)*temp_i(i2)
          RHSv_i(i2,i3min)=RHSv_i(i2,i3min)*temp_i(i2)
          RHSw_i(i2,i3min)=RHSw_i(i2,i3min)*temp_i(i2)
          RHSe_i(i2,i3min)=RHSe_i(i2,i3min)*temp_i(i2)
       enddo

       do i3=i3min+1,i3max-1
          do i2=i2min,i2max
             temp_i(i2)=1.0_wp/(d_i(i2,i3)-c_i(i2,i3-1)*a_i(i2,i3-1))
             a_i(i2,i3)=a_i(i2,i3)*temp_i(i2)
             RHS_i(i2,i3)=( RHS_i(i2,i3)-c_i(i2,i3-1)* RHS_i(i2,i3-1))*temp_i(i2)
             RHSu_i(i2,i3)=(RHSu_i(i2,i3)-c_i(i2,i3-1)*RHSu_i(i2,i3-1))*temp_i(i2)
             RHSv_i(i2,i3)=(RHSv_i(i2,i3)-c_i(i2,i3-1)*RHSv_i(i2,i3-1))*temp_i(i2)
             RHSw_i(i2,i3)=(RHSw_i(i2,i3)-c_i(i2,i3-1)*RHSw_i(i2,i3-1))*temp_i(i2)
             RHSe_i(i2,i3)=(RHSe_i(i2,i3)-c_i(i2,i3-1)*RHSe_i(i2,i3-1))*temp_i(i2)
          enddo
       enddo

       ! Backward substitution phase
       do i2=i2min,i2max
          temp_i(i2)=1.0_wp/(d_i(i2,i3max)-c_i(i2,i3max-1)*a_i(i2,i3max-1))
          RHS_i(i2,i3max)=( RHS_i(i2,i3max)-c_i(i2,i3max-1)* RHS_i(i2,i3max-1))*temp_i(i2)
          RHSu_i(i2,i3max)=(RHSu_i(i2,i3max)-c_i(i2,i3max-1)*RHSu_i(i2,i3max-1))*temp_i(i2)
          RHSv_i(i2,i3max)=(RHSv_i(i2,i3max)-c_i(i2,i3max-1)*RHSv_i(i2,i3max-1))*temp_i(i2)
          RHSw_i(i2,i3max)=(RHSw_i(i2,i3max)-c_i(i2,i3max-1)*RHSw_i(i2,i3max-1))*temp_i(i2)
          RHSe_i(i2,i3max)=(RHSe_i(i2,i3max)-c_i(i2,i3max-1)*RHSe_i(i2,i3max-1))*temp_i(i2)
       enddo

       do i3=i3max-1,i3min,-1
          do i2=i2min,i2max
             RHS_i(i2,i3)= RHS_i(i2,i3)-a_i(i2,i3)* RHS_i(i2,i3+1)
             RHSu_i(i2,i3)=RHSu_i(i2,i3)-a_i(i2,i3)*RHSu_i(i2,i3+1)
             RHSv_i(i2,i3)=RHSv_i(i2,i3)-a_i(i2,i3)*RHSv_i(i2,i3+1)
             RHSw_i(i2,i3)=RHSw_i(i2,i3)-a_i(i2,i3)*RHSw_i(i2,i3+1)
             RHSe_i(i2,i3)=RHSe_i(i2,i3)-a_i(i2,i3)*RHSe_i(i2,i3+1)
          enddo
       enddo

       ! Update with the solution
       do i3=i3min,i3max
          do i2=i2min,i2max
             Krho(i3,i2,i1)= RHS_i(i2,i3)
             Krhou(i3,i2,i1)=RHSu_i(i2,i3)
             Krhov(i3,i2,i1)=RHSv_i(i2,i3)
             Krhow(i3,i2,i1)=RHSw_i(i2,i3)
             Krhoe(i3,i2,i1)=RHSe_i(i2,i3)
          enddo
       enddo

    enddo

  end subroutine irs2_ngh_i
  
  !===============================================================================
  subroutine irs2_ngh_j
  !===============================================================================
    !> 2nd-order Implicit Residual Smoothing (IRS2) - Implicitation of j-direction
    !> Parallel inversion of tridiagonal matrix using ghost points
  !===============================================================================
    use mod_time   ! <- for deltat,is_irs
    use mod_flow   ! <- for velocities,increments,metrics
    use mod_comm1  ! <- for communication1_inc & CFL
    use mod_interface  ! <- for communication_inc
    implicit none
    ! ----------------------------------------------------------------------------
    ! indices
    integer :: i,j,k
    integer :: i1,i2,i3
    integer :: i1min,i1max,i2min,i2max,i3min,i3max
    ! contravariant velocities & spectral radius
    real(wp) :: vcmh,vcph,rspecmh,rspecph
    ! tridiagonal in j-direction
    ! --------------------------
    ! matrix elements a,d,c
    real(wp), dimension(nx,ny1_irs:ny2_irs) :: a_j,d_j,c_j    
    real(wp), dimension(nx) :: temp_j ! work array
    ! ----------------------------------------------------------------------------

    !  i3 : j;   i1 : k;   i2 : i
    i1min=1;   i2min=1;   i3min=ndy-ngh_irs(3)
    i1max=nz;  i2max=nx;  i3max=nfy+ngh_irs(4)
    if ((is_boundary(2,1)).and.(i3min<1))  i3min=1
    if ((is_boundary(2,2)).and.(i3max>ny)) i3max=ny

    ! If wall, do not take into account the wall point
    ! ================================================
    if (.not.is_wall2) then
       if (is_bc_wall(1,1)) i2min=2
       if (is_bc_wall(1,2)) i2max=nx-1
       if (is_bc_wall(3,1)) i1min=2
       if (is_bc_wall(3,2)) i1max=nz-1
       if (irk==nrk) then
          if (is_bc_wall(1,1)) then
             do i3=i3min,i3max
                Krhou(1,i3,:)=0.0_wp
                Krhov(1,i3,:)=0.0_wp
                Krhow(1,i3,:)=0.0_wp
                Krhoe(1,i3,:)=0.0_wp
             enddo
          endif
          if (is_bc_wall(1,2)) then
             do i3=i3min,i3max
                Krhou(nx,i3,:)=0.0_wp
                Krhov(nx,i3,:)=0.0_wp
                Krhow(nx,i3,:)=0.0_wp
                Krhoe(nx,i3,:)=0.0_wp
             enddo
          endif
          if (is_bc_wall(2,1)) then
             do i2=i2min,i2max
                Krhou(i2,1,:)=0.0_wp
                Krhov(i2,1,:)=0.0_wp
                Krhow(i2,1,:)=0.0_wp
                Krhoe(i2,1,:)=0.0_wp
             enddo
          endif
          if (is_bc_wall(2,2)) then
             do i2=i2min,i2max
                Krhou(i2,ny,:)=0.0_wp
                Krhov(i2,ny,:)=0.0_wp
                Krhow(i2,ny,:)=0.0_wp
                Krhoe(i2,ny,:)=0.0_wp
             enddo
          endif
          if (is_bc_wall(3,1)) then
             do i2=i2min,i2max
                Krhou(i2,:,1)=0.0_wp
                Krhov(i2,:,1)=0.0_wp
                Krhow(i2,:,1)=0.0_wp
                Krhoe(i2,:,1)=0.0_wp
             enddo
          endif
          if (is_bc_wall(3,2)) then
             do i2=i2min,i2max
                Krhou(i2,:,nz)=0.0_wp
                Krhov(i2,:,nz)=0.0_wp
                Krhow(i2,:,nz)=0.0_wp
                Krhoe(i2,:,nz)=0.0_wp
             enddo
          endif
       endif
!!$    else ! for new wall conditions, set some increments to zero
!!$       if (is_bc_wall(1,1)) then
!!$          do i3=i3min,i3max
!!$             Krhou(1,i3,:)=0.0_wp
!!$             Krhov(1,i3,:)=0.0_wp
!!$             Krhow(1,i3,:)=0.0_wp
!!$             Krhoe(1,i3,:)=0.0_wp
!!$          enddo
!!$       endif
!!$       if (is_bc_wall(1,2)) then
!!$          do i3=i3min,i3max
!!$             Krhou(nx,i3,:)=0.0_wp
!!$             Krhov(nx,i3,:)=0.0_wp
!!$             Krhow(nx,i3,:)=0.0_wp
!!$             Krhoe(nx,i3,:)=0.0_wp
!!$          enddo
!!$       endif
!!$       if (is_bc_wall(2,1)) then
!!$          do i2=i2min,i2max
!!$             Krhou(i2,1,:)=0.0_wp
!!$             Krhov(i2,1,:)=0.0_wp
!!$             Krhow(i2,1,:)=0.0_wp
!!$             Krhoe(i2,1,:)=0.0_wp
!!$          enddo
!!$       endif
!!$       if (is_bc_wall(2,2)) then
!!$          do i2=i2min,i2max
!!$             Krhou(i2,ny,:)=0.0_wp
!!$             Krhov(i2,ny,:)=0.0_wp
!!$             Krhow(i2,ny,:)=0.0_wp
!!$             Krhoe(i2,ny,:)=0.0_wp
!!$          enddo
!!$       endif
!!$       if (is_bc_wall(3,1)) then
!!$          do i2=i2min,i2max
!!$             Krhou(i2,:,1)=0.0_wp
!!$             Krhov(i2,:,1)=0.0_wp
!!$             Krhow(i2,:,1)=0.0_wp
!!$             Krhoe(i2,:,1)=0.0_wp
!!$          enddo
!!$       endif
!!$       if (is_bc_wall(3,2)) then
!!$          do i2=i2min,i2max
!!$             Krhou(i2,:,nz)=0.0_wp
!!$             Krhov(i2,:,nz)=0.0_wp
!!$             Krhow(i2,:,nz)=0.0_wp
!!$             Krhoe(i2,:,nz)=0.0_wp
!!$          enddo
!!$       endif
    endif       

    ! One-sided communications of increments in j-direction
    ! =====================================================
!!$    ! open target windows
!!$    call window_comm1_inc
!!$    ! communicate var
!!$    call communication_1_inc_j
!!$    call communication1_inc_j

    ! Calculation of the spectral radius
    ! ==================================
    if (is_dtlocal) then
       if (is_curv3) then
          do k=1,nz
             do j=1,ny
                do i=1,nx
                   vcmh=sqrt((uu(i,j-1,k)*ksi_x(i,j-1,k)+vv(i,j-1,k)*ksi_y(i,j-1,k)+ww(i,j-1,k)*ksi_z(i,j-1,k))**2 &
                            +(uu(i,j-1,k)*eta_x(i,j-1,k)+vv(i,j-1,k)*eta_y(i,j-1,k)+ww(i,j-1,k)*eta_z(i,j-1,k))**2 &
                            +(uu(i,j-1,k)*phi_x(i,j-1,k)+vv(i,j-1,k)*phi_y(i,j-1,k)+ww(i,j-1,k)*phi_z(i,j-1,k))**2 &
                            )*ijacob3(i,j-1,k)
                   vcph=sqrt((uu(i,j,k)*ksi_x(i,j,k)+vv(i,j,k)*ksi_y(i,j,k)+ww(i,j,k)*ksi_z(i,j,k))**2 &
                            +(uu(i,j,k)*eta_x(i,j,k)+vv(i,j,k)*eta_y(i,j,k)+ww(i,j,k)*eta_z(i,j,k))**2 &
                            +(uu(i,j,k)*phi_x(i,j,k)+vv(i,j,k)*phi_y(i,j,k)+ww(i,j,k)*phi_z(i,j,k))**2 &
                            )*ijacob3(i,j,k)
                   rspecmh=abs(vcmh)+c_(i,j-1,k)*sqrt(g3_ksi(i,j-1,k)**2+g3_eta(i,j-1,k)**2+g3_phi(i,j-1,k)**2)
                   rspecph=abs(vcph)+c_(i,j  ,k)*sqrt(g3_ksi(i,j  ,k)**2+g3_eta(i,j  ,k)**2+g3_phi(i,j  ,k)**2)
                   cfl_l(i,j,k)=0.5_wp*(rspecmh+rspecph)*dt_local(i,j,k)
                enddo
             enddo
          enddo
       else
          if (is_curv) then
             do k=1,nz
                do j=1,ny
                   do i=1,nx
                      vcmh=sqrt((uu(i,j-1,k)*y_eta(i,j-1)-vv(i,j-1,k)*x_eta(i,j-1))**2 &
                               +(vv(i,j-1,k)*x_ksi(i,j-1)-uu(i,j-1,k)*y_ksi(i,j-1))**2)*ijacob(i,j-1)
                      vcph=sqrt((uu(i,j  ,k)*y_eta(i,j  )-vv(i,j  ,k)*x_eta(i,j  ))**2 &
                               +(vv(i,j  ,k)*x_ksi(i,j  )-uu(i,j  ,k)*y_ksi(i,j  ))**2)*ijacob(i,j  )
                      rspecmh=abs(vcmh)+c_(i,j-1,k)*sqrt(g_ksi(i,j-1)**2+g_eta(i,j-1)**2)
                      rspecph=abs(vcph)+c_(i,j  ,k)*sqrt(g_ksi(i,j  )**2+g_eta(i,j  )**2)
                      cfl_l(i,j,k)=0.5_wp*(rspecmh+rspecph)*dt_local(i,j,k)
                   enddo
                enddo
             enddo
          else
             do k=1,nz
                do j=1,ny
                   do i=1,nx
                      rspecmh=(sqrt(uu(i,j-1,k)**2+vv(i,j-1,k)**2)+c_(i,j-1,k))*idy(j-1)
                      rspecph=(sqrt(uu(i,j  ,k)**2+vv(i,j  ,k)**2)+c_(i,j  ,k))*idy(j)
                      cfl_l(i,j,k)=0.5_wp*(rspecmh+rspecph)*dt_local(i,j,k)
                   enddo
                enddo
             enddo
          endif
       endif
    else
       if (is_curv3) then
          do k=1,nz
             do j=1,ny
                do i=1,nx
                   vcmh=sqrt((uu(i,j-1,k)*ksi_x(i,j-1,k)+vv(i,j-1,k)*ksi_y(i,j-1,k)+ww(i,j-1,k)*ksi_z(i,j-1,k))**2 &
                            +(uu(i,j-1,k)*eta_x(i,j-1,k)+vv(i,j-1,k)*eta_y(i,j-1,k)+ww(i,j-1,k)*eta_z(i,j-1,k))**2 &
                            +(uu(i,j-1,k)*phi_x(i,j-1,k)+vv(i,j-1,k)*phi_y(i,j-1,k)+ww(i,j-1,k)*phi_z(i,j-1,k))**2 &
                            )*ijacob3(i,j-1,k)
                   vcph=sqrt((uu(i,j,k)*ksi_x(i,j,k)+vv(i,j,k)*ksi_y(i,j,k)+ww(i,j,k)*ksi_z(i,j,k))**2 &
                            +(uu(i,j,k)*eta_x(i,j,k)+vv(i,j,k)*eta_y(i,j,k)+ww(i,j,k)*eta_z(i,j,k))**2 &
                            +(uu(i,j,k)*phi_x(i,j,k)+vv(i,j,k)*phi_y(i,j,k)+ww(i,j,k)*phi_z(i,j,k))**2 &
                            )*ijacob3(i,j,k)
                   rspecmh=abs(vcmh)+c_(i,j-1,k)*sqrt(g3_ksi(i,j-1,k)**2+g3_eta(i,j-1,k)**2+g3_phi(i,j-1,k)**2)
                   rspecph=abs(vcph)+c_(i,j  ,k)*sqrt(g3_ksi(i,j  ,k)**2+g3_eta(i,j  ,k)**2+g3_phi(i,j  ,k)**2)
                   cfl_l(i,j,k)=0.5_wp*(rspecmh+rspecph)*deltat
                enddo
             enddo
          enddo
       else
          if (is_curv) then
             do k=1,nz
                do j=1,ny
                   do i=1,nx
                      vcmh=sqrt((uu(i,j-1,k)*y_eta(i,j-1)-vv(i,j-1,k)*x_eta(i,j-1))**2 &
                               +(vv(i,j-1,k)*x_ksi(i,j-1)-uu(i,j-1,k)*y_ksi(i,j-1))**2)*ijacob(i,j-1)
                      vcph=sqrt((uu(i,j  ,k)*y_eta(i,j  )-vv(i,j  ,k)*x_eta(i,j  ))**2 &
                               +(vv(i,j  ,k)*x_ksi(i,j  )-uu(i,j  ,k)*y_ksi(i,j  ))**2)*ijacob(i,j  )
                      rspecmh=abs(vcmh)+c_(i,j-1,k)*sqrt(g_ksi(i,j-1)**2+g_eta(i,j-1)**2)
                      rspecph=abs(vcph)+c_(i,j  ,k)*sqrt(g_ksi(i,j  )**2+g_eta(i,j  )**2)
                      cfl_l(i,j,k)=0.5_wp*(rspecmh+rspecph)*deltat
                   enddo
                enddo
             enddo
          else
             do k=1,nz
                do j=1,ny
                   do i=1,nx
                      rspecmh=(abs(vv(i,j-1,k))+c_(i,j-1,k))*idy(j-1)
                      rspecph=(abs(vv(i,j  ,k))+c_(i,j  ,k))*idy(j)
                      cfl_l(i,j,k)=0.5_wp*(rspecmh+rspecph)*deltat
                   enddo
                enddo
             enddo
          endif
       endif
    endif

!!$    ! Check that communications are fullfilled
!!$    ! ----------------------------------------
!!$    call window_comm1_inc

!!$    ! One-sided communications of CFL in i-direction
!!$    ! ==============================================
!!$    ! open target windows
!!$    call window_comm1_cfl
!!$    ! communicate var
!!$    call communication_1_cfl_j
!!$
!!$    ! Check that communications are fullfilled
!!$    ! ----------------------------------------
!!$    call window_comm1_cfl

    call communication_inc(Krho,Krhou,Krhov,Krhow,Krhoe,cfl_l)

    do i1=i1min,i1max

       ! Fill tridiagonal matrix elements (RHS=Krho)
       ! ================================
       i3=i3min
       do i2=i2min,i2max
          a_j(i2,i3)=-theta_irs1*cfl_l(i2,i3+1,i1)
          d_j(i2,i3)=1.0_wp-a_j(i2,i3)
       enddo

       i3=i3max
       do i2=i2min,i2max
          c_j(i2,i3-1)=-theta_irs1*cfl_l(i2,i3,i1)
          d_j(i2,i3)  =1.0_wp-c_j(i2,i3-1)
       enddo

       do i3=i3min+1,i3max-1
          do i2=i2min,i2max
             c_j(i2,i3-1)=-theta_irs2*cfl_l(i2,i3,i1)**2
             a_j(i2,i3)  =-theta_irs2*cfl_l(i2,i3+1,i1)**2
             d_j(i2,i3)  =1.0_wp-(c_j(i2,i3-1)+a_j(i2,i3))
          enddo
       enddo

       ! Resolution of the tridiagonal linear system
       ! ===========================================
       ! Thomas' algorithm: Forward elimination phase
       do i2=i2min,i2max
          temp_j(i2)=1.0_wp/d_j(i2,i3min)
          a_j(i2,i3min)=a_j(i2,i3min)*temp_j(i2)
          Krho(i2,i3min,i1)= Krho(i2,i3min,i1)*temp_j(i2)
          Krhou(i2,i3min,i1)=Krhou(i2,i3min,i1)*temp_j(i2)
          Krhov(i2,i3min,i1)=Krhov(i2,i3min,i1)*temp_j(i2)
          Krhow(i2,i3min,i1)=Krhow(i2,i3min,i1)*temp_j(i2)
          Krhoe(i2,i3min,i1)=Krhoe(i2,i3min,i1)*temp_j(i2)
       enddo

       do i3=i3min+1,i3max-1
          do i2=i2min,i2max
             temp_j(i2)=1.0_wp/(d_j(i2,i3)-c_j(i2,i3-1)*a_j(i2,i3-1))
             a_j(i2,i3)=a_j(i2,i3)*temp_j(i2)
             Krho(i2,i3,i1)=( Krho(i2,i3,i1)-c_j(i2,i3-1)* Krho(i2,i3-1,i1))*temp_j(i2)
             Krhou(i2,i3,i1)=(Krhou(i2,i3,i1)-c_j(i2,i3-1)*Krhou(i2,i3-1,i1))*temp_j(i2)
             Krhov(i2,i3,i1)=(Krhov(i2,i3,i1)-c_j(i2,i3-1)*Krhov(i2,i3-1,i1))*temp_j(i2)
             Krhow(i2,i3,i1)=(Krhow(i2,i3,i1)-c_j(i2,i3-1)*Krhow(i2,i3-1,i1))*temp_j(i2)
             Krhoe(i2,i3,i1)=(Krhoe(i2,i3,i1)-c_j(i2,i3-1)*Krhoe(i2,i3-1,i1))*temp_j(i2)
          enddo
       enddo

       ! Backward substitution phase
       do i2=i2min,i2max
          temp_j(i2)=1.0_wp/(d_j(i2,i3max)-c_j(i2,i3max-1)*a_j(i2,i3max-1))
          Krho(i2,i3max,i1)=( Krho(i2,i3max,i1)-c_j(i2,i3max-1)* Krho(i2,i3max-1,i1))*temp_j(i2)
          Krhou(i2,i3max,i1)=(Krhou(i2,i3max,i1)-c_j(i2,i3max-1)*Krhou(i2,i3max-1,i1))*temp_j(i2)
          Krhov(i2,i3max,i1)=(Krhov(i2,i3max,i1)-c_j(i2,i3max-1)*Krhov(i2,i3max-1,i1))*temp_j(i2)
          Krhow(i2,i3max,i1)=(Krhow(i2,i3max,i1)-c_j(i2,i3max-1)*Krhow(i2,i3max-1,i1))*temp_j(i2)
          Krhoe(i2,i3max,i1)=(Krhoe(i2,i3max,i1)-c_j(i2,i3max-1)*Krhoe(i2,i3max-1,i1))*temp_j(i2)
       enddo

       do i3=i3max-1,i3min,-1
          do i2=i2min,i2max
             Krho(i2,i3,i1)= Krho(i2,i3,i1)-a_j(i2,i3)* Krho(i2,i3+1,i1)
             Krhou(i2,i3,i1)=Krhou(i2,i3,i1)-a_j(i2,i3)*Krhou(i2,i3+1,i1)
             Krhov(i2,i3,i1)=Krhov(i2,i3,i1)-a_j(i2,i3)*Krhov(i2,i3+1,i1)
             Krhow(i2,i3,i1)=Krhow(i2,i3,i1)-a_j(i2,i3)*Krhow(i2,i3+1,i1)
             Krhoe(i2,i3,i1)=Krhoe(i2,i3,i1)-a_j(i2,i3)*Krhoe(i2,i3+1,i1)
          enddo
       enddo
       
    enddo

  end subroutine irs2_ngh_j
  
  !===============================================================================
  subroutine irs2_ngh_k
  !===============================================================================
    !> 2nd-order Implicit Residual Smoothing (IRS2) - Implicitation of k-direction
    !> Parallel inversion of tridiagonal matrix using ghost points
  !===============================================================================
    use mod_time   ! <- for deltat,is_irs
    use mod_flow   ! <- for velocities,increments,metrics
    use mod_comm1  ! <- for communication1_inc & CFL
    use mod_comm  ! <- for communication_inc_k
    use mod_interface  ! <- for communication_inc
    implicit none
    ! ----------------------------------------------------------------------------
    ! indices
    integer :: i,j,k
    integer :: i1,i2,i3
    integer :: i1min,i1max,i2min,i2max,i3min,i3max
    ! contravariant velocities & spectral radius
    real(wp) :: vcmh,vcph,rspecmh,rspecph
    ! tridiagonal in k-direction
    ! --------------------------
    ! matrix elements a,d,c & right-hand side (RHS)
    real(wp), dimension(nx,nz1_irs:nz2_irs) :: a_k,d_k,c_k
    real(wp), dimension(nx,nz1_irs:nz2_irs) :: RHS_k,RHSu_k,RHSv_k,RHSw_k,RHSe_k
    real(wp), dimension(nx) :: temp_k ! work array
    ! ----------------------------------------------------------------------------

    !  i3 : k;   i1 : j;   i2 : i
    i1min=1;   i2min=1;   i3min=ndz-ngh_irs(5)
    i1max=ny;  i2max=nx;  i3max=nfz+ngh_irs(6)
    if ((is_boundary(3,1)).and.(i3min<1)) i3min=1
    if ((is_boundary(3,2)).and.(i3max>nz)) i3max=nz

    ! If wall, do not take into account the wall point
    ! ================================================
    if (.not.is_wall2) then
       if (is_bc_wall(1,1)) i2min=2
       if (is_bc_wall(1,2)) i2max=nx-1
       if (is_bc_wall(2,1)) i1min=2
       if (is_bc_wall(2,2)) i1max=ny-1
    endif

    ! One-sided communications of increments in k-direction
    ! =====================================================
!!$    ! open target windows
!!$    call window_comm1_inc
!!$    ! communicate var
!!$    call communication_1_inc_k
    !call communication1_inc_k

    ! Calculation of the local CFL number at ind-1/2
    ! ==============================================
    if (is_dtlocal) then
       if (is_curv3) then
          do k=1,nz
             do j=1,ny
                do i=1,nx
                   vcmh=sqrt((uu(i,j,k-1)*ksi_x(i,j,k-1)+vv(i,j,k-1)*ksi_y(i,j,k-1)+ww(i,j,k-1)*ksi_z(i,j,k-1))**2 &
                            +(uu(i,j,k-1)*eta_x(i,j,k-1)+vv(i,j,k-1)*eta_y(i,j,k-1)+ww(i,j,k-1)*eta_z(i,j,k-1))**2 &
                            +(uu(i,j,k-1)*phi_x(i,j,k-1)+vv(i,j,k-1)*phi_y(i,j,k-1)+ww(i,j,k-1)*phi_z(i,j,k-1))**2 &
                            )*ijacob3(i,j,k-1)
                   vcph=sqrt((uu(i,j,k)*ksi_x(i,j,k)+vv(i,j,k)*ksi_y(i,j,k)+ww(i,j,k)*ksi_z(i,j,k))**2 &
                            +(uu(i,j,k)*eta_x(i,j,k)+vv(i,j,k)*eta_y(i,j,k)+ww(i,j,k)*eta_z(i,j,k))**2 &
                            +(uu(i,j,k)*phi_x(i,j,k)+vv(i,j,k)*phi_y(i,j,k)+ww(i,j,k)*phi_z(i,j,k))**2 &
                            )*ijacob3(i,j,k)
                   rspecmh=abs(vcmh)+c_(i,j,k-1)*sqrt(g3_ksi(i,j,k-1)**2+g3_eta(i,j,k-1)**2+g3_phi(i,j,k-1)**2)
                   rspecph=abs(vcph)+c_(i,j,k  )*sqrt(g3_ksi(i,j,k  )**2+g3_eta(i,j,k  )**2+g3_phi(i,j,k  )**2)
                   cfl_l(i,j,k)=0.5_wp*(rspecmh+rspecph)*dt_local(i,j,k)
                enddo
             enddo
          enddo
       else
          do k=1,nz
             do j=1,ny
                do i=1,nx
                   rspecmh=(abs(ww(i,j,k-1))+c_(i,j,k-1))*idz(k-1)
                   rspecph=(abs(ww(i,j,k  ))+c_(i,j,k  ))*idz(k  )
                   cfl_l(i,j,k)=0.5_wp*(rspecmh+rspecph)*dt_local(i,j,k)
                enddo
             enddo
          enddo
       endif
    else
       if (is_curv3) then
          do k=1,nz
             do j=1,ny
                do i=1,nx
                   vcmh=sqrt((uu(i,j,k-1)*ksi_x(i,j,k-1)+vv(i,j,k-1)*ksi_y(i,j,k-1)+ww(i,j,k-1)*ksi_z(i,j,k-1))**2 &
                            +(uu(i,j,k-1)*eta_x(i,j,k-1)+vv(i,j,k-1)*eta_y(i,j,k-1)+ww(i,j,k-1)*eta_z(i,j,k-1))**2 &
                            +(uu(i,j,k-1)*phi_x(i,j,k-1)+vv(i,j,k-1)*phi_y(i,j,k-1)+ww(i,j,k-1)*phi_z(i,j,k-1))**2 &
                            )*ijacob3(i,j,k-1)
                   vcph=sqrt((uu(i,j,k)*ksi_x(i,j,k)+vv(i,j,k)*ksi_y(i,j,k)+ww(i,j,k)*ksi_z(i,j,k))**2 &
                            +(uu(i,j,k)*eta_x(i,j,k)+vv(i,j,k)*eta_y(i,j,k)+ww(i,j,k)*eta_z(i,j,k))**2 &
                            +(uu(i,j,k)*phi_x(i,j,k)+vv(i,j,k)*phi_y(i,j,k)+ww(i,j,k)*phi_z(i,j,k))**2 &
                            )*ijacob3(i,j,k)
                   rspecmh=abs(vcmh)+c_(i,j,k-1)*sqrt(g3_ksi(i,j,k-1)**2+g3_eta(i,j,k-1)**2+g3_phi(i,j,k-1)**2)
                   rspecph=abs(vcph)+c_(i,j,k  )*sqrt(g3_ksi(i,j,k  )**2+g3_eta(i,j,k  )**2+g3_phi(i,j,k  )**2)
                   cfl_l(i,j,k)=0.5_wp*(rspecmh+rspecph)*deltat
                enddo
             enddo
          enddo
       else
          do k=1,nz
             do j=1,ny
                do i=1,nx
                   rspecmh=(abs(ww(i,j,k-1))+c_(i,j,k-1))*idz(k-1)
                   rspecph=(abs(ww(i,j,k  ))+c_(i,j,k  ))*idz(k  )
                   cfl_l(i,j,k)=0.5_wp*(rspecmh+rspecph)*deltat
                enddo
             enddo
          enddo
       endif
    endif

!!$    ! Check that communications are fullfilled
!!$    ! ----------------------------------------
!!$    call window_comm1_inc

    ! One-sided communications of CFL in i-direction
    ! ==============================================
    ! open target windows
!    call window_comm1_cfl
    ! communicate var
!    call communication_1_cfl_k

    ! Check that communications are fullfilled
    ! ----------------------------------------
!    call window_comm1_cfl

    ! Communication of ghost points
    ! =============================
    !call communication_inc_k(Krho,Krhou,Krhov,Krhow,Krhoe,cfl_l)
    call communication_inc(Krho,Krhou,Krhov,Krhow,Krhoe,cfl_l)

    do i1=i1min,i1max
       
       ! Fill tridiagonal matrix elements & RHS
       ! ======================================
       i3=i3min
       do i2=i2min,i2max
          a_k(i2,i3)=-theta_irs1*cfl_l(i2,i1,i3+1)
          d_k(i2,i3)=1.0_wp-a_k(i2,i3)
       enddo

       i3=i3max
       do i2=i2min,i2max
          c_k(i2,i3-1)=-theta_irs1*cfl_l(i2,i1,i3)
          d_k(i2,i3)  =1.0_wp-c_k(i2,i3-1)
       enddo

       do i3=i3min+1,i3max-1
          do i2=i2min,i2max
             c_k(i2,i3-1)=-theta_irs2*cfl_l(i2,i1,i3)**2
             a_k(i2,i3)  =-theta_irs2*cfl_l(i2,i1,i3+1)**2
             d_k(i2,i3)  =1.0_wp-(c_k(i2,i3-1)+a_k(i2,i3))
          enddo
       enddo

       ! Filling of the RHS
       do i3=i3min,i3max
          do i2=i2min,i2max
             RHS_k(i2,i3)= Krho(i2,i1,i3)
             RHSu_k(i2,i3)=Krhou(i2,i1,i3)
             RHSv_k(i2,i3)=Krhov(i2,i1,i3)
             RHSw_k(i2,i3)=Krhow(i2,i1,i3)
             RHSe_k(i2,i3)=Krhoe(i2,i1,i3)
          enddo
       enddo

       ! Resolution of the tridiagonal linear system
       ! ===========================================
       ! Thomas' algorithm-Forward elimination phase
       do i2=i2min,i2max
          temp_k(i2)=1.0_wp/d_k(i2,i3min)
          a_k(i2,i3min)=a_k(i2,i3min)*temp_k(i2)
          RHS_k(i2,i3min)= RHS_k(i2,i3min)*temp_k(i2)
          RHSu_k(i2,i3min)=RHSu_k(i2,i3min)*temp_k(i2)
          RHSv_k(i2,i3min)=RHSv_k(i2,i3min)*temp_k(i2)
          RHSw_k(i2,i3min)=RHSw_k(i2,i3min)*temp_k(i2)
          RHSe_k(i2,i3min)=RHSe_k(i2,i3min)*temp_k(i2)
       enddo

       do i3=i3min+1,i3max-1
          do i2=i2min,i2max
             temp_k(i2)=1.0_wp/(d_k(i2,i3)-c_k(i2,i3-1)*a_k(i2,i3-1))
             a_k(i2,i3)=a_k(i2,i3)*temp_k(i2)
             RHS_k(i2,i3)=( RHS_k(i2,i3)-c_k(i2,i3-1)* RHS_k(i2,i3-1))*temp_k(i2)
             RHSu_k(i2,i3)=(RHSu_k(i2,i3)-c_k(i2,i3-1)*RHSu_k(i2,i3-1))*temp_k(i2)
             RHSv_k(i2,i3)=(RHSv_k(i2,i3)-c_k(i2,i3-1)*RHSv_k(i2,i3-1))*temp_k(i2)
             RHSw_k(i2,i3)=(RHSw_k(i2,i3)-c_k(i2,i3-1)*RHSw_k(i2,i3-1))*temp_k(i2)
             RHSe_k(i2,i3)=(RHSe_k(i2,i3)-c_k(i2,i3-1)*RHSe_k(i2,i3-1))*temp_k(i2)
          enddo
       enddo

       ! Backward substitution phase
       do i2=i2min,i2max
          temp_k(i2)=1.0_wp/(d_k(i2,i3max)-c_k(i2,i3max-1)*a_k(i2,i3max-1))
          RHS_k(i2,i3max)=( RHS_k(i2,i3max)-c_k(i2,i3max-1)* RHS_k(i2,i3max-1))*temp_k(i2)
          RHSu_k(i2,i3max)=(RHSu_k(i2,i3max)-c_k(i2,i3max-1)*RHSu_k(i2,i3max-1))*temp_k(i2)
          RHSv_k(i2,i3max)=(RHSv_k(i2,i3max)-c_k(i2,i3max-1)*RHSv_k(i2,i3max-1))*temp_k(i2)
          RHSw_k(i2,i3max)=(RHSw_k(i2,i3max)-c_k(i2,i3max-1)*RHSw_k(i2,i3max-1))*temp_k(i2)
          RHSe_k(i2,i3max)=(RHSe_k(i2,i3max)-c_k(i2,i3max-1)*RHSe_k(i2,i3max-1))*temp_k(i2)
       enddo

       do i3=i3max-1,i3min,-1
          do i2=i2min,i2max
             RHS_k(i2,i3)= RHS_k(i2,i3)-a_k(i2,i3)* RHS_k(i2,i3+1)
             RHSu_k(i2,i3)=RHSu_k(i2,i3)-a_k(i2,i3)*RHSu_k(i2,i3+1)
             RHSv_k(i2,i3)=RHSv_k(i2,i3)-a_k(i2,i3)*RHSv_k(i2,i3+1)
             RHSw_k(i2,i3)=RHSw_k(i2,i3)-a_k(i2,i3)*RHSw_k(i2,i3+1)
             RHSe_k(i2,i3)=RHSe_k(i2,i3)-a_k(i2,i3)*RHSe_k(i2,i3+1)
          enddo
       enddo

       ! Update with the solution
       do i3=i3min,i3max
          do i2=i2min,i2max
             Krho(i2,i1,i3)= RHS_k(i2,i3)
             Krhou(i2,i1,i3)=RHSu_k(i2,i3)
             Krhov(i2,i1,i3)=RHSv_k(i2,i3)
             Krhow(i2,i1,i3)=RHSw_k(i2,i3)
             Krhoe(i2,i1,i3)=RHSe_k(i2,i3)
          enddo
       enddo
       
    enddo

  end subroutine irs2_ngh_k
  
  !===============================================================================
  subroutine irs4_ngh_i
  !===============================================================================
    !> 4th-order Implicit Residual Smoothing (IRS4) - Implicitation of i-direction
    !> Parallel inversion of pentadiagonal matrix using ghost points
  !===============================================================================
    use mod_time   ! <- for deltat,is_irs
    use mod_flow   ! <- for velocities,increments,metrics
    use mod_comm1  ! <- for communication1_inc & cfl
    use mod_interface  ! <- for communication_inc
    implicit none
    ! ----------------------------------------------------------------------------
    ! indices
    integer :: i,j,k
    integer :: i1,i2,i3
    integer :: i1min,i1max,i2min,i2max,i3min,i3max
    ! contravariant velocities & spectral radius
    real(wp) :: vcmh,vcph,rspecmh,rspecph
    ! pentadiagonal in i-direction
    ! ----------------------------
    ! matrix elements b,a,d,c,e & right-hand side (RHS)
    real(wp), dimension(ny,nx1_irs:nx2_irs) :: b_i,a_i,d_i,c_i,e_i
    real(wp), dimension(ny,nx1_irs:nx2_irs) :: RHS_i,RHSu_i,RHSv_i,RHSw_i,RHSe_i
    real(wp), dimension(ny) :: temp_i1,temp_i2 ! work arrays
    ! ----------------------------------------------------------------------------
    real(wp) :: theta
    
    !  i3 : i;   i1 : k;   i2 : j
    i1min=1;   i2min=1;   i3min=ndx-ngh_irs(1)
    i1max=nz;  i2max=ny;  i3max=nfx+ngh_irs(2)
    if ((is_boundary(1,1)).and.(i3min<1))  i3min=1
    if ((is_boundary(1,2)).and.(i3max>nx)) i3max=nx

    ! If wall, do not take into account the wall point
    ! ================================================
    if (.not.is_wall2) then
       if (is_bc_wall(2,1)) i2min=2
       if (is_bc_wall(2,2)) i2max=ny-1
       if (is_bc_wall(3,1)) i1min=2
       if (is_bc_wall(3,2)) i1max=nz-1
       if (irk==nrk) then
          if (is_bc_wall(1,1)) then
             do i2=i2min,i2max
                Krhou(1,i2,:)=0.0_wp
                Krhov(1,i2,:)=0.0_wp
                Krhow(1,i2,:)=0.0_wp
                Krhoe(1,i2,:)=0.0_wp
             enddo
          endif
          if (is_bc_wall(1,2)) then
             do i2=i2min,i2max
                Krhou(nx,i2,:)=0.0_wp
                Krhov(nx,i2,:)=0.0_wp
                Krhow(nx,i2,:)=0.0_wp
                Krhoe(nx,i2,:)=0.0_wp
             enddo
          endif
          if (is_bc_wall(2,1)) then
             do i3=i3min,i3max
                Krhou(i3,1,:)=0.0_wp
                Krhov(i3,1,:)=0.0_wp
                Krhow(i3,1,:)=0.0_wp
                Krhoe(i3,1,:)=0.0_wp
             enddo
          endif
          if (is_bc_wall(2,2)) then
             do i3=i3min,i3max
                Krhou(i3,ny,:)=0.0_wp
                Krhov(i3,ny,:)=0.0_wp
                Krhow(i3,ny,:)=0.0_wp
                Krhoe(i3,ny,:)=0.0_wp
             enddo
          endif
          if (is_bc_wall(3,1)) then
             do i3=i3min,i3max
                Krhou(i3,:,1)=0.0_wp
                Krhov(i3,:,1)=0.0_wp
                Krhow(i3,:,1)=0.0_wp
                Krhoe(i3,:,1)=0.0_wp
             enddo
          endif
          if (is_bc_wall(3,2)) then
             do i3=i3min,i3max
                Krhou(i3,:,nz)=0.0_wp
                Krhov(i3,:,nz)=0.0_wp
                Krhow(i3,:,nz)=0.0_wp
                Krhoe(i3,:,nz)=0.0_wp
             enddo
          endif
       endif
    else ! for new wall conditions, set some increments to zero
       if (is_bc_wall(1,1)) then
          do i2=i2min,i2max
             Krhou(1,i2,:)=0.0_wp
             Krhov(1,i2,:)=0.0_wp
             Krhow(1,i2,:)=0.0_wp
             Krhoe(1,i2,:)=0.0_wp
          enddo
       endif
       if (is_bc_wall(1,2)) then
          do i2=i2min,i2max
             Krhou(nx,i2,:)=0.0_wp
             Krhov(nx,i2,:)=0.0_wp
             Krhow(nx,i2,:)=0.0_wp
             Krhoe(nx,i2,:)=0.0_wp
          enddo
       endif
       if (is_bc_wall(2,1)) then
          do i3=i3min,i3max
             Krhou(i3,1,:)=0.0_wp
             Krhov(i3,1,:)=0.0_wp
             Krhow(i3,1,:)=0.0_wp
             Krhoe(i3,1,:)=0.0_wp
          enddo
       endif
       if (is_bc_wall(2,2)) then
          do i3=i3min,i3max
             Krhou(i3,ny,:)=0.0_wp
             Krhov(i3,ny,:)=0.0_wp
             Krhow(i3,ny,:)=0.0_wp
             Krhoe(i3,ny,:)=0.0_wp
          enddo
       endif
       if (is_bc_wall(3,1)) then
          do i3=i3min,i3max
             Krhou(i3,:,1)=0.0_wp
             Krhov(i3,:,1)=0.0_wp
             Krhow(i3,:,1)=0.0_wp
             Krhoe(i3,:,1)=0.0_wp
          enddo
       endif
       if (is_bc_wall(3,2)) then
          do i3=i3min,i3max
             Krhou(i3,:,nz)=0.0_wp
             Krhov(i3,:,nz)=0.0_wp
             Krhow(i3,:,nz)=0.0_wp
             Krhoe(i3,:,nz)=0.0_wp
          enddo
       endif
    endif

    ! One-sided communications of increments in i-direction
    ! =====================================================
!!$    ! open target windows
!!$    call window_comm1_inc
!!$    ! communicate var
!!$    !call communication_1_inc_i
!!$    call communication_1_inc_ij

!!$    call window_comm1_inc(Krho,Krhou,Krhov,Krhow,Krhoe)
!!$    call communication_1_inc_ij
!!$    call window_comm1_inc(Krho,Krhou,Krhov,Krhow,Krhoe)

    !call communication1_inc_ij
!!$  call communication1_inc_i

!!$    ! start check
!!$    call MPI_WIN_FENCE(0,win_cfl, info)
!!$
!!$    ! start commnunications
!!$    call commun1_inc_ij(cfl_l,win_cfl)
!!$
!!$    ! check end of commnunication
!!$    call MPI_WIN_FENCE(0,win_cfl, info)
   
    ! Calculation of the local CFL number at ind-1/2
    ! ==============================================
    if (is_dtlocal) then
       if (is_curv3) then
          do k=1,nz
             do j=1,ny
                do i=1,nx
                   vcmh=sqrt((uu(i-1,j,k)*ksi_x(i-1,j,k)+vv(i-1,j,k)*ksi_y(i-1,j,k)+ww(i-1,j,k)*ksi_z(i-1,j,k))**2 &
                            +(uu(i-1,j,k)*eta_x(i-1,j,k)+vv(i-1,j,k)*eta_y(i-1,j,k)+ww(i-1,j,k)*eta_z(i-1,j,k))**2 &
                            +(uu(i-1,j,k)*phi_x(i-1,j,k)+vv(i-1,j,k)*phi_y(i-1,j,k)+ww(i-1,j,k)*phi_z(i-1,j,k))**2 &
                            )*ijacob3(i-1,j,k)
                   vcph=sqrt((uu(i,j,k)*ksi_x(i,j,k)+vv(i,j,k)*ksi_y(i,j,k)+ww(i,j,k)*ksi_z(i,j,k))**2 &
                            +(uu(i,j,k)*eta_x(i,j,k)+vv(i,j,k)*eta_y(i,j,k)+ww(i,j,k)*eta_z(i,j,k))**2 &
                            +(uu(i,j,k)*phi_x(i,j,k)+vv(i,j,k)*phi_y(i,j,k)+ww(i,j,k)*phi_z(i,j,k))**2 &
                            )*ijacob3(i,j,k)
                   rspecmh=abs(vcmh)+c_(i-1,j,k)*sqrt(g3_ksi(i-1,j,k)**2+g3_eta(i-1,j,k)**2+g3_phi(i-1,j,k)**2)
                   rspecph=abs(vcph)+c_(i  ,j,k)*sqrt(g3_ksi(i  ,j,k)**2+g3_eta(i  ,j,k)**2+g3_phi(i  ,j,k)**2)
                   cfl_l(i,j,k)=0.5_wp*(rspecmh+rspecph)*dt_local(i,j,k)
                enddo
             enddo
          enddo
       else
          if (is_curv) then
             do k=1,nz
                do j=1,ny
                   do i=1,nx
                      vcmh=sqrt((uu(i-1,j,k)*y_eta(i-1,j)-vv(i-1,j,k)*x_eta(i-1,j))**2 &
                               +(vv(i-1,j,k)*x_ksi(i-1,j)-uu(i-1,j,k)*y_ksi(i-1,j))**2)*ijacob(i-1,j)
                      vcph=sqrt((uu(i,j  ,k)*y_eta(i,j  )-vv(i,j  ,k)*x_eta(i,j  ))**2 &
                               +(vv(i,j  ,k)*x_ksi(i,j  )-uu(i,j  ,k)*y_ksi(i,j  ))**2)*ijacob(i,j)
                      rspecmh=abs(vcmh)+c_(i-1,j,k)*sqrt(g_ksi(i-1,j)**2+g_eta(i-1,j)**2)
                      rspecph=abs(vcph)+c_(i  ,j,k)*sqrt(g_ksi(i  ,j)**2+g_eta(i  ,j)**2)
                      cfl_l(i,j,k)=0.5_wp*(rspecmh+rspecph)*dt_local(i,j,k)
                   enddo
                enddo
             enddo
          else
             do k=1,nz
                do j=1,ny
                   do i=1,nx
                      rspecmh=(abs(uu(i-1,j,k))+c_(i-1,j,k))*idx(i-1)
                      rspecph=(abs(uu(i  ,j,k))+c_(i  ,j,k))*idx(i)
                      cfl_l(i,j,k)=0.5_wp*(rspecmh+rspecph)*dt_local(i,j,k)
                   enddo
                enddo
             enddo
          endif
       endif
    else
       if (is_curv3) then
          do k=1,nz
             do j=1,ny
                do i=1,nx
                   vcmh=sqrt((uu(i-1,j,k)*ksi_x(i-1,j,k)+vv(i-1,j,k)*ksi_y(i-1,j,k)+ww(i-1,j,k)*ksi_z(i-1,j,k))**2 &
                            +(uu(i-1,j,k)*eta_x(i-1,j,k)+vv(i-1,j,k)*eta_y(i-1,j,k)+ww(i-1,j,k)*eta_z(i-1,j,k))**2 &
                            +(uu(i-1,j,k)*phi_x(i-1,j,k)+vv(i-1,j,k)*phi_y(i-1,j,k)+ww(i-1,j,k)*phi_z(i-1,j,k))**2 &
                            )*ijacob3(i-1,j,k)
                   vcph=sqrt((uu(i,j,k)*ksi_x(i,j,k)+vv(i,j,k)*ksi_y(i,j,k)+ww(i,j,k)*ksi_z(i,j,k))**2 &
                            +(uu(i,j,k)*eta_x(i,j,k)+vv(i,j,k)*eta_y(i,j,k)+ww(i,j,k)*eta_z(i,j,k))**2 &
                            +(uu(i,j,k)*phi_x(i,j,k)+vv(i,j,k)*phi_y(i,j,k)+ww(i,j,k)*phi_z(i,j,k))**2 &
                            )*ijacob3(i,j,k)
                   rspecmh=abs(vcmh)+c_(i-1,j,k)*sqrt(g3_ksi(i-1,j,k)**2+g3_eta(i-1,j,k)**2+g3_phi(i-1,j,k)**2)
                   rspecph=abs(vcph)+c_(i  ,j,k)*sqrt(g3_ksi(i  ,j,k)**2+g3_eta(i  ,j,k)**2+g3_phi(i  ,j,k)**2)
                   cfl_l(i,j,k)=0.5_wp*(rspecmh+rspecph)*deltat
                enddo
             enddo
          enddo
       else
          if (is_curv) then
             do k=1,nz
                do j=1,ny
                   do i=1,nx
                      !vcmh=(uu(i-1,j,k)*y_eta(i-1,j)-vv(i-1,j,k)*x_eta(i-1,j))*ijacob(i-1,j)
                      !vcph=(uu(i  ,j,k)*y_eta(i  ,j)-vv(i  ,j,k)*x_eta(i  ,j))*ijacob(i  ,j)
                      !rspecmh=abs(vcmh)+c_(i-1,j,k)*g_ksi(i-1,j)
                      !rspecph=abs(vcph)+c_(i  ,j,k)*g_ksi(i  ,j)
                      vcmh=sqrt((uu(i-1,j,k)*y_eta(i-1,j)-vv(i-1,j,k)*x_eta(i-1,j))**2 &
                               +(vv(i-1,j,k)*x_ksi(i-1,j)-uu(i-1,j,k)*y_ksi(i-1,j))**2)*ijacob(i-1,j)
                      vcph=sqrt((uu(i,j  ,k)*y_eta(i,j  )-vv(i,j  ,k)*x_eta(i,j  ))**2 &
                               +(vv(i,j  ,k)*x_ksi(i,j  )-uu(i,j  ,k)*y_ksi(i,j  ))**2)*ijacob(i,j)
!!$                      vcmh=max(uu(i-1,j,k)*y_eta(i-1,j)-vv(i-1,j,k)*x_eta(i-1,j), &
!!$                               vv(i-1,j,k)*x_ksi(i-1,j)-uu(i-1,j,k)*y_ksi(i-1,j))*ijacob(i-1,j)
!!$                      vcph=max(uu(i,j  ,k)*y_eta(i,j  )-vv(i,j  ,k)*x_eta(i,j  ), &
!!$                               vv(i,j  ,k)*x_ksi(i,j  )-uu(i,j  ,k)*y_ksi(i,j  ))*ijacob(i,j)
                      rspecmh=abs(vcmh)+c_(i-1,j,k)*sqrt(g_ksi(i-1,j)**2+g_eta(i-1,j)**2)
                      rspecph=abs(vcph)+c_(i  ,j,k)*sqrt(g_ksi(i  ,j)**2+g_eta(i  ,j)**2)
!!$                      rspecmh=max(abs(uu(i-1,j,k)*y_eta(i-1,j)-vv(i-1,j,k)*x_eta(i-1,j)*ijacob(i-1,j))+c_(i-1,j,k)*g_ksi(i-1,j),&
!!$                                  abs(vv(i-1,j,k)*x_ksi(i-1,j)-uu(i-1,j,k)*y_ksi(i-1,j)*ijacob(i-1,j))+c_(i-1,j,k)*g_eta(i-1,j))
!!$                      rspecph=max(abs(uu(i  ,j,k)*y_eta(i  ,j)-vv(i  ,j,k)*x_eta(i  ,j)*ijacob(i  ,j))+c_(i  ,j,k)*g_ksi(i  ,j),&
!!$                                  abs(vv(i  ,j,k)*x_ksi(i  ,j)-uu(i  ,j,k)*y_ksi(i  ,j)*ijacob(i  ,j))+c_(i  ,j,k)*g_eta(i  ,j))
!!$                      rspecmh=max(abs(uu(i-1,j,k))+c_(i-1,j,k),&
!!$                                  abs(vv(i-1,j,k))+c_(i-1,j,k))*ijacob(i-1,j)
!!$                      rspecph=max(abs(uu(i  ,j,k))+c_(i  ,j,k),&
!!$                                  abs(vv(i  ,j,k))+c_(i  ,j,k))*ijacob(i  ,j)
                      cfl_l(i,j,k)=0.5_wp*(rspecmh+rspecph)*deltat
!!$                cfl_l(i,j,k)=rspecph*deltat
                   enddo
                enddo
             enddo
          else
             do k=1,nz
                do j=1,ny
                   do i=1,nx
                      rspecmh=(abs(uu(i-1,j,k))+c_(i-1,j,k))*idx(i-1)
                      rspecph=(abs(uu(i  ,j,k))+c_(i  ,j,k))*idx(i)
                      cfl_l(i,j,k)=0.5_wp*(rspecmh+rspecph)*deltat
                   enddo
                enddo
             enddo
          endif
       endif
    endif

!!$    ! Check that communications are fullfilled
!!$    ! ----------------------------------------
!!$    call window_comm1_inc

    ! One-sided communications of CFL in i-direction
    ! ==============================================
!!$    ! open target windows
!!$    call window_comm1_cfl
!!$    ! communicate var
!!$    call communication_1_cfl_ij
!!$    !call communication_1_cfl_i
!!$
!!$    ! Check that communications are fullfilled
!!$    ! ----------------------------------------
!!$    call window_comm1_cfl
    
    !call communication1_cfl_i
    call communication_inc(Krho,Krhou,Krhov,Krhow,Krhoe,cfl_l)
    
    do i1=i1min,i1max

       ! Construction of the pentadiagonal matrix
       ! ========================================
       ! The linear systems is Mx=y with :
       !     _M=(m_ij) for i,j in {1,..,n} and with :
       !                 m_ii=d_i for i in {1,..,n}
       !                 m_ii-2=e_i-2 for i in {3,..,n}
       !                 m_ii-1=c_i-1 for i in {2,..,n}
       !                 m_ii+2=b_i for i in {1,..,n-2}
       !                 m_ii+1=a_i for i in {1,..,n-1}
       !     _ y=RHS=dwi
       !     _d and y are of dimensions n, a and c of dimensions n-1,
       !           e and b of dimensions n-2
       i3=i3min
       do i2=i2min,i2max
          b_i(i2,i3)=0.0_wp
          a_i(i2,i3)=-theta_irs1*cfl_l(i3+1,i2,i1)
          d_i(i2,i3)=1.0_wp-a_i(i2,i3)
       enddo

       i3=i3max
       do i2=i2min,i2max
          c_i(i2,i3-1)=-theta_irs1*cfl_l(i3,i2,i1)
          d_i(i2,i3)  =1.0_wp-c_i(i2,i3-1)
          e_i(i2,i3-2)=0.0_wp
       enddo

       i3=i3min+1
       do i2=i2min,i2max
          b_i(i2,i3)  =0.0_wp
          a_i(i2,i3)  =-theta_irs2*cfl_l(i3+1,i2,i1)**2
          c_i(i2,i3-1)=-theta_irs2*cfl_l(i3,i2,i1)**2
          d_i(i2,i3)  =1.0_wp-(a_i(i2,i3)+c_i(i2,i3-1))
       enddo

       i3=i3max-1
       do i2=i2min,i2max
          a_i(i2,i3)  =-theta_irs2*cfl_l(i3+1,i2,i1)**2
          c_i(i2,i3-1)=-theta_irs2*cfl_l(i3,i2,i1)**2
          d_i(i2,i3)  =1.0_wp-(a_i(i2,i3)+c_i(i2,i3-1))
          e_i(i2,i3-2)=0.0_wp
       enddo

       do i3=i3min+2,i3max-2
          do i2=i2min,i2max
             !theta=1.2*theta_irs4*(1.0_wp-exp(-(cfl_l(i3,i2,i1)**8)/100.0_wp))/cfl_l(i3,i2,i1)**0.12_wp
             theta=theta_irs4
             e_i(i2,i3-2)=theta*cfl_l(i3,i2,i1)**4
             b_i(i2,i3)  =theta*cfl_l(i3+1,i2,i1)**4
             c_i(i2,i3-1)=-3.0_wp*e_i(i2,i3-2)-b_i(i2,i3)
             a_i(i2,i3)  =-3.0_wp*b_i(i2,i3)-e_i(i2,i3-2)
             d_i(i2,i3)  =1.0_wp+3.0_wp*(e_i(i2,i3-2)+b_i(i2,i3))
          enddo
       enddo

       ! Filling of the RHS
       do i3=i3min,i3max
          do i2=i2min,i2max
             RHS_i(i2,i3)= Krho(i3,i2,i1)
             RHSu_i(i2,i3)=Krhou(i3,i2,i1)
             RHSv_i(i2,i3)=Krhov(i3,i2,i1)
             RHSw_i(i2,i3)=Krhow(i3,i2,i1)
             RHSe_i(i2,i3)=Krhoe(i3,i2,i1)
          enddo
       enddo
       
       ! Resolution of the pentadiagonal linear system
       ! =============================================
       i3=i3min
       do i2=i2min,i2max
          temp_i1(i2)=1.0_wp/d_i(i2,i3)
          a_i(i2,i3)=a_i(i2,i3)*temp_i1(i2)
          b_i(i2,i3)=b_i(i2,i3)*temp_i1(i2)
          RHS_i(i2,i3)= RHS_i(i2,i3)*temp_i1(i2)
          RHSu_i(i2,i3)=RHSu_i(i2,i3)*temp_i1(i2)
          RHSv_i(i2,i3)=RHSv_i(i2,i3)*temp_i1(i2)
          RHSw_i(i2,i3)=RHSw_i(i2,i3)*temp_i1(i2)
          RHSe_i(i2,i3)=RHSe_i(i2,i3)*temp_i1(i2)
       enddo

       i3=i3min+1
       do i2=i2min,i2max
          temp_i2(i2)=c_i(i2,i3-1)
          temp_i1(i2)=1.0_wp/(d_i(i2,i3)-a_i(i2,i3-1)*temp_i2(i2))
          a_i(i2,i3)=(a_i(i2,i3)-b_i(i2,i3-1)*temp_i2(i2))*temp_i1(i2)
          b_i(i2,i3)=b_i(i2,i3)*temp_i1(i2)
          RHS_i(i2,i3)=( RHS_i(i2,i3)- RHS_i(i2,i3-1)*temp_i2(i2))*temp_i1(i2)
          RHSu_i(i2,i3)=(RHSu_i(i2,i3)-RHSu_i(i2,i3-1)*temp_i2(i2))*temp_i1(i2)
          RHSv_i(i2,i3)=(RHSv_i(i2,i3)-RHSv_i(i2,i3-1)*temp_i2(i2))*temp_i1(i2)
          RHSw_i(i2,i3)=(RHSw_i(i2,i3)-RHSw_i(i2,i3-1)*temp_i2(i2))*temp_i1(i2)
          RHSe_i(i2,i3)=(RHSe_i(i2,i3)-RHSe_i(i2,i3-1)*temp_i2(i2))*temp_i1(i2)
       enddo

       ! Step 5: for i=3,...,n-2
       do i3=i3min+2,i3max-2
          do i2=i2min,i2max
             temp_i2(i2)=c_i(i2,i3-1)-a_i(i2,i3-2)*e_i(i2,i3-2)
             temp_i1(i2)=1.0_wp/(d_i(i2,i3)-b_i(i2,i3-2)*e_i(i2,i3-2)-a_i(i2,i3-1)*temp_i2(i2))
             a_i(i2,i3)=(a_i(i2,i3)-b_i(i2,i3-1)*temp_i2(i2))*temp_i1(i2)
             b_i(i2,i3)=b_i(i2,i3)*temp_i1(i2)
             RHS_i(i2,i3)=( RHS_i(i2,i3)- RHS_i(i2,i3-2)*e_i(i2,i3-2)- RHS_i(i2,i3-1)*temp_i2(i2))*temp_i1(i2)
             RHSu_i(i2,i3)=(RHSu_i(i2,i3)-RHSu_i(i2,i3-2)*e_i(i2,i3-2)-RHSu_i(i2,i3-1)*temp_i2(i2))*temp_i1(i2)
             RHSv_i(i2,i3)=(RHSv_i(i2,i3)-RHSv_i(i2,i3-2)*e_i(i2,i3-2)-RHSv_i(i2,i3-1)*temp_i2(i2))*temp_i1(i2)
             RHSw_i(i2,i3)=(RHSw_i(i2,i3)-RHSw_i(i2,i3-2)*e_i(i2,i3-2)-RHSw_i(i2,i3-1)*temp_i2(i2))*temp_i1(i2)
             RHSe_i(i2,i3)=(RHSe_i(i2,i3)-RHSe_i(i2,i3-2)*e_i(i2,i3-2)-RHSe_i(i2,i3-1)*temp_i2(i2))*temp_i1(i2)
          enddo
       enddo

       ! Step 5: for i=n-1, n
       i3=i3max-1
       do i2=i2min,i2max
          temp_i2(i2)=c_i(i2,i3-1)-a_i(i2,i3-2)*e_i(i2,i3-2)
          temp_i1(i2)=1.0_wp/(d_i(i2,i3)-b_i(i2,i3-2)*e_i(i2,i3-2)-a_i(i2,i3-1)*temp_i2(i2))
          a_i(i2,i3)=(a_i(i2,i3)-b_i(i2,i3-1)*temp_i2(i2))*temp_i1(i2)
          RHS_i(i2,i3)=( RHS_i(i2,i3)- RHS_i(i2,i3-2)*e_i(i2,i3-2)- RHS_i(i2,i3-1)*temp_i2(i2))*temp_i1(i2)
          RHSu_i(i2,i3)=(RHSu_i(i2,i3)-RHSu_i(i2,i3-2)*e_i(i2,i3-2)-RHSu_i(i2,i3-1)*temp_i2(i2))*temp_i1(i2)
          RHSv_i(i2,i3)=(RHSv_i(i2,i3)-RHSv_i(i2,i3-2)*e_i(i2,i3-2)-RHSv_i(i2,i3-1)*temp_i2(i2))*temp_i1(i2)
          RHSw_i(i2,i3)=(RHSw_i(i2,i3)-RHSw_i(i2,i3-2)*e_i(i2,i3-2)-RHSw_i(i2,i3-1)*temp_i2(i2))*temp_i1(i2)
          RHSe_i(i2,i3)=(RHSe_i(i2,i3)-RHSe_i(i2,i3-2)*e_i(i2,i3-2)-RHSe_i(i2,i3-1)*temp_i2(i2))*temp_i1(i2)
       enddo

       i3=i3max
       do i2=i2min,i2max
          temp_i2(i2)=c_i(i2,i3-1)-a_i(i2,i3-2)*e_i(i2,i3-2)
          temp_i1(i2)=1.0_wp/(d_i(i2,i3)-b_i(i2,i3-2)*e_i(i2,i3-2)-a_i(i2,i3-1)*temp_i2(i2))
          RHS_i(i2,i3)=( RHS_i(i2,i3)- RHS_i(i2,i3-2)*e_i(i2,i3-2)- RHS_i(i2,i3-1)*temp_i2(i2))*temp_i1(i2)
          RHSu_i(i2,i3)=(RHSu_i(i2,i3)-RHSu_i(i2,i3-2)*e_i(i2,i3-2)-RHSu_i(i2,i3-1)*temp_i2(i2))*temp_i1(i2)
          RHSv_i(i2,i3)=(RHSv_i(i2,i3)-RHSv_i(i2,i3-2)*e_i(i2,i3-2)-RHSv_i(i2,i3-1)*temp_i2(i2))*temp_i1(i2)
          RHSw_i(i2,i3)=(RHSw_i(i2,i3)-RHSw_i(i2,i3-2)*e_i(i2,i3-2)-RHSw_i(i2,i3-1)*temp_i2(i2))*temp_i1(i2)
          RHSe_i(i2,i3)=(RHSe_i(i2,i3)-RHSe_i(i2,i3-2)*e_i(i2,i3-2)-RHSe_i(i2,i3-1)*temp_i2(i2))*temp_i1(i2)
       enddo


       ! Step 6: computation of the solution vector
       i3=i3max-1
       do i2=i2min,i2max
          RHS_i(i2,i3)= RHS_i(i2,i3)-a_i(i2,i3)* RHS_i(i2,i3+1)
          RHSu_i(i2,i3)=RHSu_i(i2,i3)-a_i(i2,i3)*RHSu_i(i2,i3+1)
          RHSv_i(i2,i3)=RHSv_i(i2,i3)-a_i(i2,i3)*RHSv_i(i2,i3+1)
          RHSw_i(i2,i3)=RHSw_i(i2,i3)-a_i(i2,i3)*RHSw_i(i2,i3+1)
          RHSe_i(i2,i3)=RHSe_i(i2,i3)-a_i(i2,i3)*RHSe_i(i2,i3+1)
       enddo

       do i3=i3max-2,i3min,-1
          do i2=i2min,i2max
             RHS_i(i2,i3)= RHS_i(i2,i3)-a_i(i2,i3)* RHS_i(i2,i3+1)-b_i(i2,i3)* RHS_i(i2,i3+2)
             RHSu_i(i2,i3)=RHSu_i(i2,i3)-a_i(i2,i3)*RHSu_i(i2,i3+1)-b_i(i2,i3)*RHSu_i(i2,i3+2)
             RHSv_i(i2,i3)=RHSv_i(i2,i3)-a_i(i2,i3)*RHSv_i(i2,i3+1)-b_i(i2,i3)*RHSv_i(i2,i3+2)
             RHSw_i(i2,i3)=RHSw_i(i2,i3)-a_i(i2,i3)*RHSw_i(i2,i3+1)-b_i(i2,i3)*RHSw_i(i2,i3+2)
             RHSe_i(i2,i3)=RHSe_i(i2,i3)-a_i(i2,i3)*RHSe_i(i2,i3+1)-b_i(i2,i3)*RHSe_i(i2,i3+2)
          enddo
       enddo

       ! Update with the solution
       do i3=i3min,i3max
          do i2=i2min,i2max
             Krho(i3,i2,i1)= RHS_i(i2,i3)
             Krhou(i3,i2,i1)=RHSu_i(i2,i3)
             Krhov(i3,i2,i1)=RHSv_i(i2,i3)
             Krhow(i3,i2,i1)=RHSw_i(i2,i3)
             Krhoe(i3,i2,i1)=RHSe_i(i2,i3)
          enddo
       enddo
       
    enddo

  end subroutine irs4_ngh_i
  
  !===============================================================================
  subroutine irs4_ngh_j
  !===============================================================================
    !> 4th-order Implicit Residual Smoothing (IRS4) - Implicitation of j-direction
    !> Parallel inversion of pentadiagonal matrix using ghost points
  !===============================================================================
    use mod_time   ! <- for deltat,is_irs
    use mod_flow   ! <- for velocities,increments,metrics
    use mod_comm1  ! <- for communication1_inc & cfl
    use mod_interface  ! <- for communication_inc
    implicit none
    ! ----------------------------------------------------------------------------
    ! indices
    integer :: i,j,k
    integer :: i1,i2,i3
    integer :: i1min,i1max,i2min,i2max,i3min,i3max
    ! contravariant velocities & spectral radius
    real(wp) :: vcmh,vcph,rspecmh,rspecph
    ! pentadiagonal in j-direction
    ! ----------------------------
    ! matrix elements b,a,d,c,e
    real(wp), dimension(nx,ny1_irs:ny2_irs) :: b_j,a_j,d_j,c_j,e_j
    real(wp), dimension(nx) :: temp_j1,temp_j2 ! work arrays
    ! ----------------------------------------------------------------------------
    real(wp) :: theta

    !  i3 : j;   i1 : k;   i2 : i
    i1min=1;   i2min=1;   i3min=ndy-ngh_irs(3)
    i1max=nz;  i2max=nx;  i3max=nfy+ngh_irs(4)
    if ((is_boundary(2,1)).and.(i3min<1))  i3min=1
    if ((is_boundary(2,2)).and.(i3max>ny)) i3max=ny

    ! If wall, do not take into account the wall point
    ! ================================================
    if (.not.is_wall2) then
       if (is_bc_wall(1,1)) i2min=2
       if (is_bc_wall(1,2)) i2max=nx-1
       if (is_bc_wall(3,1)) i1min=2
       if (is_bc_wall(3,2)) i1max=nz-1
       if (irk==nrk) then
          if (is_bc_wall(1,1)) then
             do i3=i3min,i3max
                Krhou(1,i3,:)=0.0_wp
                Krhov(1,i3,:)=0.0_wp
                Krhow(1,i3,:)=0.0_wp
                Krhoe(1,i3,:)=0.0_wp
             enddo
          endif
          if (is_bc_wall(1,2)) then
             do i3=i3min,i3max
                Krhou(nx,i3,:)=0.0_wp
                Krhov(nx,i3,:)=0.0_wp
                Krhow(nx,i3,:)=0.0_wp
                Krhoe(nx,i3,:)=0.0_wp
             enddo
          endif
          if (is_bc_wall(2,1)) then
             do i2=i2min,i2max
                Krhou(i2,1,:)=0.0_wp
                Krhov(i2,1,:)=0.0_wp
                Krhow(i2,1,:)=0.0_wp
                Krhoe(i2,1,:)=0.0_wp
             enddo
          endif
          if (is_bc_wall(2,2)) then
             do i2=i2min,i2max
                Krhou(i2,ny,:)=0.0_wp
                Krhov(i2,ny,:)=0.0_wp
                Krhow(i2,ny,:)=0.0_wp
                Krhoe(i2,ny,:)=0.0_wp
             enddo
          endif
          if (is_bc_wall(3,1)) then
             do i2=i2min,i2max
                Krhou(i2,:,1)=0.0_wp
                Krhov(i2,:,1)=0.0_wp
                Krhow(i2,:,1)=0.0_wp
                Krhoe(i2,:,1)=0.0_wp
             enddo
          endif
          if (is_bc_wall(3,2)) then
             do i2=i2min,i2max
                Krhou(i2,:,nz)=0.0_wp
                Krhov(i2,:,nz)=0.0_wp
                Krhow(i2,:,nz)=0.0_wp
                Krhoe(i2,:,nz)=0.0_wp
             enddo
          endif
       endif
    else ! for new wall conditions, set some increments to zero
       if (is_bc_wall(1,1)) then
          do i3=i3min,i3max
             Krhou(1,i3,:)=0.0_wp
             Krhov(1,i3,:)=0.0_wp
             Krhow(1,i3,:)=0.0_wp
             Krhoe(1,i3,:)=0.0_wp
          enddo
       endif
       if (is_bc_wall(1,2)) then
          do i3=i3min,i3max
             Krhou(nx,i3,:)=0.0_wp
             Krhov(nx,i3,:)=0.0_wp
             Krhow(nx,i3,:)=0.0_wp
             Krhoe(nx,i3,:)=0.0_wp
          enddo
       endif
       if (is_bc_wall(2,1)) then
          do i2=i2min,i2max
             Krhou(i2,1,:)=0.0_wp
             Krhov(i2,1,:)=0.0_wp
             Krhow(i2,1,:)=0.0_wp
             Krhoe(i2,1,:)=0.0_wp
          enddo
       endif
       if (is_bc_wall(2,2)) then
          do i2=i2min,i2max
             Krhou(i2,ny,:)=0.0_wp
             Krhov(i2,ny,:)=0.0_wp
             Krhow(i2,ny,:)=0.0_wp
             Krhoe(i2,ny,:)=0.0_wp
          enddo
       endif
       if (is_bc_wall(3,1)) then
          do i2=i2min,i2max
             Krhou(i2,:,1)=0.0_wp
             Krhov(i2,:,1)=0.0_wp
             Krhow(i2,:,1)=0.0_wp
             Krhoe(i2,:,1)=0.0_wp
          enddo
       endif
       if (is_bc_wall(3,2)) then
          do i2=i2min,i2max
             Krhou(i2,:,nz)=0.0_wp
             Krhov(i2,:,nz)=0.0_wp
             Krhow(i2,:,nz)=0.0_wp
             Krhoe(i2,:,nz)=0.0_wp
          enddo
       endif
    endif       

    ! One-sided communications of increments in j-direction
    ! =====================================================
!!$    ! open target windows
!!$    call window_comm1_inc
!!$    ! communicate var
!!$    !call communication_1_inc_j
!!$    call communication_1_inc_ij
    !call communication1_inc_ij

    ! Calculation of the local CFL number at ind-1/2
    ! ==============================================
    if (is_dtlocal) then
       if (is_curv3) then
          do k=1,nz
             do j=1,ny
                do i=1,nx
                   vcmh=sqrt((uu(i,j-1,k)*ksi_x(i,j-1,k)+vv(i,j-1,k)*ksi_y(i,j-1,k)+ww(i,j-1,k)*ksi_z(i,j-1,k))**2 &
                            +(uu(i,j-1,k)*eta_x(i,j-1,k)+vv(i,j-1,k)*eta_y(i,j-1,k)+ww(i,j-1,k)*eta_z(i,j-1,k))**2 &
                            +(uu(i,j-1,k)*phi_x(i,j-1,k)+vv(i,j-1,k)*phi_y(i,j-1,k)+ww(i,j-1,k)*phi_z(i,j-1,k))**2 &
                            )*ijacob3(i,j-1,k)
                   vcph=sqrt((uu(i,j,k)*ksi_x(i,j,k)+vv(i,j,k)*ksi_y(i,j,k)+ww(i,j,k)*ksi_z(i,j,k))**2 &
                            +(uu(i,j,k)*eta_x(i,j,k)+vv(i,j,k)*eta_y(i,j,k)+ww(i,j,k)*eta_z(i,j,k))**2 &
                            +(uu(i,j,k)*phi_x(i,j,k)+vv(i,j,k)*phi_y(i,j,k)+ww(i,j,k)*phi_z(i,j,k))**2 &
                            )*ijacob3(i,j,k)
                   rspecmh=abs(vcmh)+c_(i,j-1,k)*sqrt(g3_ksi(i,j-1,k)**2+g3_eta(i,j-1,k)**2+g3_phi(i,j-1,k)**2)
                   rspecph=abs(vcph)+c_(i,j  ,k)*sqrt(g3_ksi(i,j  ,k)**2+g3_eta(i,j  ,k)**2+g3_phi(i,j  ,k)**2)
                   cfl_l(i,j,k)=0.5_wp*(rspecmh+rspecph)*dt_local(i,j,k)
                enddo
             enddo
          enddo
       else
          if (is_curv) then
             do k=1,nz
                do j=1,ny
                   do i=1,nx
                      vcmh=sqrt((uu(i,j-1,k)*y_eta(i,j-1)-vv(i,j-1,k)*x_eta(i,j-1))**2 &
                               +(vv(i,j-1,k)*x_ksi(i,j-1)-uu(i,j-1,k)*y_ksi(i,j-1))**2)*ijacob(i,j-1)
                      vcph=sqrt((uu(i,j  ,k)*y_eta(i,j  )-vv(i,j  ,k)*x_eta(i,j  ))**2 &
                               +(vv(i,j  ,k)*x_ksi(i,j  )-uu(i,j  ,k)*y_ksi(i,j  ))**2)*ijacob(i,j  )
                      rspecmh=abs(vcmh)+c_(i,j-1,k)*sqrt(g_ksi(i,j-1)**2+g_eta(i,j-1)**2)
                      rspecph=abs(vcph)+c_(i,j  ,k)*sqrt(g_ksi(i,j  )**2+g_eta(i,j  )**2)
                      cfl_l(i,j,k)=0.5_wp*(rspecmh+rspecph)*dt_local(i,j,k)
                   enddo
                enddo
             enddo
          else
             do k=1,nz
                do j=1,ny
                   do i=1,nx
                      rspecmh=(abs(vv(i,j-1,k))+c_(i,j-1,k))*idy(j-1)
                      rspecph=(abs(vv(i,j  ,k))+c_(i,j  ,k))*idy(j)
                      cfl_l(i,j,k)=0.5_wp*(rspecmh+rspecph)*dt_local(i,j,k)
                   enddo
                enddo
             enddo
          endif
       endif
    else
       if (is_curv3) then
          do k=1,nz
             do j=1,ny
                do i=1,nx
                   vcmh=sqrt((uu(i,j-1,k)*ksi_x(i,j-1,k)+vv(i,j-1,k)*ksi_y(i,j-1,k)+ww(i,j-1,k)*ksi_z(i,j-1,k))**2 &
                            +(uu(i,j-1,k)*eta_x(i,j-1,k)+vv(i,j-1,k)*eta_y(i,j-1,k)+ww(i,j-1,k)*eta_z(i,j-1,k))**2 &
                            +(uu(i,j-1,k)*phi_x(i,j-1,k)+vv(i,j-1,k)*phi_y(i,j-1,k)+ww(i,j-1,k)*phi_z(i,j-1,k))**2 &
                            )*ijacob3(i,j-1,k)
                   vcph=sqrt((uu(i,j,k)*ksi_x(i,j,k)+vv(i,j,k)*ksi_y(i,j,k)+ww(i,j,k)*ksi_z(i,j,k))**2 &
                            +(uu(i,j,k)*eta_x(i,j,k)+vv(i,j,k)*eta_y(i,j,k)+ww(i,j,k)*eta_z(i,j,k))**2 &
                            +(uu(i,j,k)*phi_x(i,j,k)+vv(i,j,k)*phi_y(i,j,k)+ww(i,j,k)*phi_z(i,j,k))**2 &
                            )*ijacob3(i,j,k)
                   rspecmh=abs(vcmh)+c_(i,j-1,k)*sqrt(g3_ksi(i,j-1,k)**2+g3_eta(i,j-1,k)**2+g3_phi(i,j-1,k)**2)
                   rspecph=abs(vcph)+c_(i,j  ,k)*sqrt(g3_ksi(i,j  ,k)**2+g3_eta(i,j  ,k)**2+g3_phi(i,j  ,k)**2)
                   cfl_l(i,j,k)=0.5_wp*(rspecmh+rspecph)*deltat
                enddo
             enddo
          enddo
       else
          if (is_curv) then
             do k=1,nz
                do j=1,ny
                   do i=1,nx
                      !vcmh=(vv(i,j-1,k)*x_ksi(i,j-1)-uu(i,j-1,k)*y_ksi(i,j-1))*ijacob(i,j-1)
                      !vcph=(vv(i,j  ,k)*x_ksi(i,j  )-uu(i,j  ,k)*y_ksi(i,j  ))*ijacob(i,j  )
                      !rspecmh=abs(vcmh)+c_(i,j-1,k)*g_eta(i,j-1)
                      !rspecph=abs(vcph)+c_(i,j  ,k)*g_eta(i,j  )
                      vcmh=sqrt((uu(i,j-1,k)*y_eta(i,j-1)-vv(i,j-1,k)*x_eta(i,j-1))**2 &
                               +(vv(i,j-1,k)*x_ksi(i,j-1)-uu(i,j-1,k)*y_ksi(i,j-1))**2)*ijacob(i,j-1)
                      vcph=sqrt((uu(i,j  ,k)*y_eta(i,j  )-vv(i,j  ,k)*x_eta(i,j  ))**2 &
                               +(vv(i,j  ,k)*x_ksi(i,j  )-uu(i,j  ,k)*y_ksi(i,j  ))**2)*ijacob(i,j)
!!$                vcmh=max(uu(i,j-1,k)*y_eta(i,j-1)-vv(i,j-1,k)*x_eta(i,j-1), &
!!$                         vv(i,j-1,k)*x_ksi(i,j-1)-uu(i,j-1,k)*y_ksi(i,j-1))*ijacob(i,j-1)
!!$                vcph=max(uu(i,j  ,k)*y_eta(i,j  )-vv(i,j  ,k)*x_eta(i,j  ), &
!!$                         vv(i,j  ,k)*x_ksi(i,j  )-uu(i,j  ,k)*y_ksi(i,j  ))*ijacob(i,j)
                      rspecmh=abs(vcmh)+c_(i,j-1,k)*sqrt(g_ksi(i,j-1)**2+g_eta(i,j-1)**2)
                      rspecph=abs(vcph)+c_(i,j  ,k)*sqrt(g_ksi(i,j  )**2+g_eta(i,j  )**2)
!!$                rspecmh=max(abs(uu(i,j-1,k)*y_eta(i,j-1)-vv(i,j-1,k)*x_eta(i,j-1)*ijacob(i,j-1))+c_(i,j-1,k)*g_ksi(i,j-1),&
!!$                            abs(vv(i,j-1,k)*x_ksi(i,j-1)-uu(i,j-1,k)*y_ksi(i,j-1)*ijacob(i,j-1))+c_(i,j-1,k)*g_eta(i,j-1))
!!$                rspecph=max(abs(uu(i,j  ,k)*y_eta(i,j  )-vv(i,j  ,k)*x_eta(i,j  )*ijacob(i,j  ))+c_(i,j  ,k)*g_ksi(i,j  ),&
!!$                            abs(vv(i,j  ,k)*x_ksi(i,j  )-uu(i,j  ,k)*y_ksi(i,j  )*ijacob(i,j  ))+c_(i,j  ,k)*g_eta(i,j  ))
!!$                rspecmh=max(abs(uu(i,j-1,k))+c_(i,j-1,k),&
!!$                            abs(vv(i,j-1,k))+c_(i,j-1,k))*ijacob(i,j-1)
!!$                rspecph=max(abs(uu(i,j  ,k))+c_(i,j  ,k),&
!!$                            abs(vv(i,j  ,k))+c_(i,j  ,k))*ijacob(i,j  )
                      cfl_l(i,j,k)=0.5_wp*(rspecmh+rspecph)*deltat
!!$                cfl_l(i,j,k)=rspecph*deltat
                   enddo
                enddo
             enddo
          else
             do k=1,nz
                do j=1,ny
                   do i=1,nx
                      rspecmh=(abs(vv(i,j-1,k))+c_(i,j-1,k))*idy(j-1)
                      rspecph=(abs(vv(i,j  ,k))+c_(i,j  ,k))*idy(j)
                      cfl_l(i,j,k)=0.5_wp*(rspecmh+rspecph)*deltat
                   enddo
                enddo
             enddo
          endif
       endif
    endif

!!$    ! Check that communications are fullfilled
!!$    ! ----------------------------------------
!!$    call window_comm1_inc

    ! Communication of ghost points
    ! =============================
    !call communication1_cfl_j
    !call communication1_cfl_j
    call communication_inc(Krho,Krhou,Krhov,Krhow,Krhoe,cfl_l)

    do i1=i1min,i1max

       ! Construction of the pentadiagonal matrix (RHS=Krho)
       ! ========================================
       i3=i3min
       do i2=i2min,i2max
          b_j(i2,i3)=0.0_wp
          a_j(i2,i3)=-theta_irs1*cfl_l(i2,i3+1,i1)
          d_j(i2,i3)=1.0_wp-a_j(i2,i3)
       enddo

       i3=i3max
       do i2=i2min,i2max
          c_j(i2,i3-1)=-theta_irs1*cfl_l(i2,i3,i1)
          d_j(i2,i3)  =1.0_wp-c_j(i2,i3-1)
          e_j(i2,i3-2)=0.0_wp
       enddo

       i3=i3min+1
       do i2=i2min,i2max
          b_j(i2,i3)  =0.0_wp
          a_j(i2,i3)  =-theta_irs2*cfl_l(i2,i3+1,i1)**2
          c_j(i2,i3-1)=-theta_irs2*cfl_l(i2,i3,i1)**2
          d_j(i2,i3)  =1.0_wp-(a_j(i2,i3)+c_j(i2,i3-1))
       enddo

       i3=i3max-1
       do i2=i2min,i2max
          a_j(i2,i3)  =-theta_irs2*cfl_l(i2,i3+1,i1)**2
          c_j(i2,i3-1)=-theta_irs2*cfl_l(i2,i3,i1)**2
          d_j(i2,i3)  =1.0_wp-(a_j(i2,i3)+c_j(i2,i3-1))
          e_j(i2,i3-2)=0.0_wp
       enddo

       do i3=i3min+2,i3max-2
          do i2=i2min,i2max
             !theta=1.2*theta_irs4*(1.0_wp-exp(-(cfl_l(i3,i2,i1)**8)/100.0_wp))/cfl_l(i3,i2,i1)**0.12_wp
             theta=theta_irs4
             e_j(i2,i3-2)=theta*cfl_l(i2,i3,i1)**4
             b_j(i2,i3)  =theta*cfl_l(i2,i3+1,i1)**4
             c_j(i2,i3-1)=-3.0_wp*e_j(i2,i3-2)-b_j(i2,i3)
             a_j(i2,i3)  =-3.0_wp*b_j(i2,i3)-e_j(i2,i3-2)
             d_j(i2,i3)  =1.0_wp+3.0_wp*(e_j(i2,i3-2)+b_j(i2,i3))
          enddo
       enddo

       ! Resolution of the pentadiagonal linear system
       i3=i3min
       do i2=i2min,i2max
          temp_j1(i2)=1.0_wp/d_j(i2,i3)
          a_j(i2,i3)=a_j(i2,i3)*temp_j1(i2)
          b_j(i2,i3)=b_j(i2,i3)*temp_j1(i2)
          Krho(i2,i3,i1)= Krho(i2,i3,i1)*temp_j1(i2)
          Krhou(i2,i3,i1)=Krhou(i2,i3,i1)*temp_j1(i2)
          Krhov(i2,i3,i1)=Krhov(i2,i3,i1)*temp_j1(i2)
          Krhow(i2,i3,i1)=Krhow(i2,i3,i1)*temp_j1(i2)
          Krhoe(i2,i3,i1)=Krhoe(i2,i3,i1)*temp_j1(i2)
       enddo

       i3=i3min+1
       do i2=i2min,i2max
          temp_j2(i2)=c_j(i2,i3-1)
          temp_j1(i2)=1.0_wp/(d_j(i2,i3)-a_j(i2,i3-1)*temp_j2(i2))
          a_j(i2,i3)=(a_j(i2,i3)-b_j(i2,i3-1)*temp_j2(i2))*temp_j1(i2)
          b_j(i2,i3)=b_j(i2,i3)*temp_j1(i2)
          Krho(i2,i3,i1)=( Krho(i2,i3,i1)- Krho(i2,i3-1,i1)*temp_j2(i2))*temp_j1(i2)
          Krhou(i2,i3,i1)=(Krhou(i2,i3,i1)-Krhou(i2,i3-1,i1)*temp_j2(i2))*temp_j1(i2)
          Krhov(i2,i3,i1)=(Krhov(i2,i3,i1)-Krhov(i2,i3-1,i1)*temp_j2(i2))*temp_j1(i2)
          Krhow(i2,i3,i1)=(Krhow(i2,i3,i1)-Krhow(i2,i3-1,i1)*temp_j2(i2))*temp_j1(i2)
          Krhoe(i2,i3,i1)=(Krhoe(i2,i3,i1)-Krhoe(i2,i3-1,i1)*temp_j2(i2))*temp_j1(i2)
       enddo

       ! Step 5: for i=3,...,n-2
       do i3=i3min+2,i3max-2
          do i2=i2min,i2max
             temp_j2(i2)=c_j(i2,i3-1)-a_j(i2,i3-2)*e_j(i2,i3-2)
             temp_j1(i2)=1.0_wp/(d_j(i2,i3)-b_j(i2,i3-2)*e_j(i2,i3-2) &
                  -a_j(i2,i3-1)*temp_j2(i2))
             a_j(i2,i3)=(a_j(i2,i3)-b_j(i2,i3-1)*temp_j2(i2))*temp_j1(i2)
             b_j(i2,i3)=b_j(i2,i3)*temp_j1(i2)
             Krho(i2,i3,i1) =( Krho(i2,i3,i1)- Krho(i2,i3-2,i1)*e_j(i2,i3-2)- Krho(i2,i3-1,i1)*temp_j2(i2))*temp_j1(i2)
             Krhou(i2,i3,i1)=(Krhou(i2,i3,i1)-Krhou(i2,i3-2,i1)*e_j(i2,i3-2)-Krhou(i2,i3-1,i1)*temp_j2(i2))*temp_j1(i2)
             Krhov(i2,i3,i1)=(Krhov(i2,i3,i1)-Krhov(i2,i3-2,i1)*e_j(i2,i3-2)-Krhov(i2,i3-1,i1)*temp_j2(i2))*temp_j1(i2)
             Krhow(i2,i3,i1)=(Krhow(i2,i3,i1)-Krhow(i2,i3-2,i1)*e_j(i2,i3-2)-Krhow(i2,i3-1,i1)*temp_j2(i2))*temp_j1(i2)
             Krhoe(i2,i3,i1)=(Krhoe(i2,i3,i1)-Krhoe(i2,i3-2,i1)*e_j(i2,i3-2)-Krhoe(i2,i3-1,i1)*temp_j2(i2))*temp_j1(i2)
          enddo
       enddo

       ! Step 5: for i=n-1, n
       i3=i3max-1
       do i2=i2min,i2max
          temp_j2(i2)=c_j(i2,i3-1)-a_j(i2,i3-2)*e_j(i2,i3-2)
          temp_j1(i2)=1.0_wp/(d_j(i2,i3)-b_j(i2,i3-2)*e_j(i2,i3-2)-a_j(i2,i3-1)*temp_j2(i2))
          a_j(i2,i3)=(a_j(i2,i3)-b_j(i2,i3-1)*temp_j2(i2))*temp_j1(i2)
          Krho(i2,i3,i1) =( Krho(i2,i3,i1)- Krho(i2,i3-2,i1)*e_j(i2,i3-2)- Krho(i2,i3-1,i1)*temp_j2(i2))*temp_j1(i2)
          Krhou(i2,i3,i1)=(Krhou(i2,i3,i1)-Krhou(i2,i3-2,i1)*e_j(i2,i3-2)-Krhou(i2,i3-1,i1)*temp_j2(i2))*temp_j1(i2)
          Krhov(i2,i3,i1)=(Krhov(i2,i3,i1)-Krhov(i2,i3-2,i1)*e_j(i2,i3-2)-Krhov(i2,i3-1,i1)*temp_j2(i2))*temp_j1(i2)
          Krhow(i2,i3,i1)=(Krhow(i2,i3,i1)-Krhow(i2,i3-2,i1)*e_j(i2,i3-2)-Krhow(i2,i3-1,i1)*temp_j2(i2))*temp_j1(i2)
          Krhoe(i2,i3,i1)=(Krhoe(i2,i3,i1)-Krhoe(i2,i3-2,i1)*e_j(i2,i3-2)-Krhoe(i2,i3-1,i1)*temp_j2(i2))*temp_j1(i2)
       enddo

       i3=i3max
       do i2=i2min,i2max
          temp_j2(i2)=c_j(i2,i3-1)-a_j(i2,i3-2)*e_j(i2,i3-2)
          temp_j1(i2)=1.0_wp/(d_j(i2,i3)-b_j(i2,i3-2)*e_j(i2,i3-2)-a_j(i2,i3-1)*temp_j2(i2))
          Krho(i2,i3,i1) =( Krho(i2,i3,i1)- Krho(i2,i3-2,i1)*e_j(i2,i3-2)- Krho(i2,i3-1,i1)*temp_j2(i2))*temp_j1(i2)
          Krhou(i2,i3,i1)=(Krhou(i2,i3,i1)-Krhou(i2,i3-2,i1)*e_j(i2,i3-2)-Krhou(i2,i3-1,i1)*temp_j2(i2))*temp_j1(i2)
          Krhov(i2,i3,i1)=(Krhov(i2,i3,i1)-Krhov(i2,i3-2,i1)*e_j(i2,i3-2)-Krhov(i2,i3-1,i1)*temp_j2(i2))*temp_j1(i2)
          Krhow(i2,i3,i1)=(Krhow(i2,i3,i1)-Krhow(i2,i3-2,i1)*e_j(i2,i3-2)-Krhow(i2,i3-1,i1)*temp_j2(i2))*temp_j1(i2)
          Krhoe(i2,i3,i1)=(Krhoe(i2,i3,i1)-Krhoe(i2,i3-2,i1)*e_j(i2,i3-2)-Krhoe(i2,i3-1,i1)*temp_j2(i2))*temp_j1(i2)
       enddo


       ! Step 6: computation of the solution vector
       i3=i3max-1
       do i2=i2min,i2max
          Krho(i2,i3,i1)= Krho(i2,i3,i1)-a_j(i2,i3)* Krho(i2,i3+1,i1)
          Krhou(i2,i3,i1)=Krhou(i2,i3,i1)-a_j(i2,i3)*Krhou(i2,i3+1,i1)
          Krhov(i2,i3,i1)=Krhov(i2,i3,i1)-a_j(i2,i3)*Krhov(i2,i3+1,i1)
          Krhow(i2,i3,i1)=Krhow(i2,i3,i1)-a_j(i2,i3)*Krhow(i2,i3+1,i1)
          Krhoe(i2,i3,i1)=Krhoe(i2,i3,i1)-a_j(i2,i3)*Krhoe(i2,i3+1,i1)
       enddo

       do i3=i3max-2,i3min,-1
          do i2=i2min,i2max
             Krho(i2,i3,i1)= Krho(i2,i3,i1)-a_j(i2,i3)* Krho(i2,i3+1,i1)-b_j(i2,i3)* Krho(i2,i3+2,i1)
             Krhou(i2,i3,i1)=Krhou(i2,i3,i1)-a_j(i2,i3)*Krhou(i2,i3+1,i1)-b_j(i2,i3)*Krhou(i2,i3+2,i1)
             Krhov(i2,i3,i1)=Krhov(i2,i3,i1)-a_j(i2,i3)*Krhov(i2,i3+1,i1)-b_j(i2,i3)*Krhov(i2,i3+2,i1)
             Krhow(i2,i3,i1)=Krhow(i2,i3,i1)-a_j(i2,i3)*Krhow(i2,i3+1,i1)-b_j(i2,i3)*Krhow(i2,i3+2,i1)
             Krhoe(i2,i3,i1)=Krhoe(i2,i3,i1)-a_j(i2,i3)*Krhoe(i2,i3+1,i1)-b_j(i2,i3)*Krhoe(i2,i3+2,i1)
          enddo
       enddo
       
    enddo

  end subroutine irs4_ngh_j
  
  !===============================================================================
  subroutine irs4_ngh_k
  !===============================================================================
    !> 4th-order Implicit Residual Smoothing (IRS4) - Implicitation of k-direction
    !> Parallel inversion of pentadiagonal matrix using ghost points
  !===============================================================================
    use mod_time   ! <- for deltat,is_irs
    use mod_flow   ! <- for velocities,increments,metrics
    use mod_comm1  ! <- for communication1_inc & cfl
    use mod_comm   ! <- for communication_inc_k
    use mod_interface  ! <- for communication_inc
    implicit none
    ! ----------------------------------------------------------------------------
    ! indices
    integer :: i,j,k
    integer :: i1,i2,i3
    integer :: i1min,i1max,i2min,i2max,i3min,i3max
    ! contravariant velocities & spectral radius
    real(wp) :: vcmh,vcph,rspecmh,rspecph
    ! pentadiagonal in k-direction
    ! ----------------------------
    ! matrix elements b,a,d,c,e & right-hand side (RHS)
    real(wp), dimension(nx,nz1_irs:nz2_irs) :: b_k,a_k,d_k,c_k,e_k
    real(wp), dimension(nx,nz1_irs:nz2_irs) :: RHS_k,RHSu_k,RHSv_k,RHSw_k,RHSe_k
    real(wp), dimension(nx) :: temp_k1,temp_k2 ! work arrays
    ! ----------------------------------------------------------------------------

    !  i3 : k;   i1 : i;   i2 : j
    i1min=1;   i2min=1;   i3min=ndz-ngh_irs(5)
    i1max=ny;  i2max=nx;  i3max=nfz+ngh_irs(6)
    if ((is_boundary(3,1)).and.(i3min<1)) i3min=1
    if ((is_boundary(3,2)).and.(i3max>nz)) i3max=nz

    ! If wall, do not take into account the wall point
    ! ================================================
    if (.not.is_wall2) then
       if (is_bc_wall(1,1)) i2min=2
       if (is_bc_wall(1,2)) i2max=nx-1
       if (is_bc_wall(2,1)) i1min=2
       if (is_bc_wall(2,2)) i1max=ny-1
       if (irk==nrk) then
          if (is_bc_wall(1,1)) then
             do i1=i1min,i1max
                Krhou(1,i1,:)=0.0_wp
                Krhov(1,i1,:)=0.0_wp
                Krhow(1,i1,:)=0.0_wp
                Krhoe(1,i1,:)=0.0_wp
             enddo
          endif
          if (is_bc_wall(1,2)) then
             do i1=i1min,i1max
                Krhou(nx,i1,:)=0.0_wp
                Krhov(nx,i1,:)=0.0_wp
                Krhow(nx,i1,:)=0.0_wp
                Krhoe(nx,i1,:)=0.0_wp
             enddo
          endif
          if (is_bc_wall(2,1)) then
             do i2=i2min,i2max
                Krhou(i2,1,:)=0.0_wp
                Krhov(i2,1,:)=0.0_wp
                Krhow(i2,1,:)=0.0_wp
                Krhoe(i2,1,:)=0.0_wp
             enddo
          endif
          if (is_bc_wall(2,2)) then
             do i2=i2min,i2max
                Krhou(i2,ny,:)=0.0_wp
                Krhov(i2,ny,:)=0.0_wp
                Krhow(i2,ny,:)=0.0_wp
                Krhoe(i2,ny,:)=0.0_wp
             enddo
          endif
          if (is_bc_wall(3,1)) then
             do i2=i2min,i2max
                Krhou(i2,:,1)=0.0_wp
                Krhov(i2,:,1)=0.0_wp
                Krhow(i2,:,1)=0.0_wp
                Krhoe(i2,:,1)=0.0_wp
             enddo
          endif
          if (is_bc_wall(3,2)) then
             do i2=i2min,i2max
                Krhou(i2,:,nz)=0.0_wp
                Krhov(i2,:,nz)=0.0_wp
                Krhow(i2,:,nz)=0.0_wp
                Krhoe(i2,:,nz)=0.0_wp
             enddo
          endif
       endif
    else ! for new wall conditions, set some increments to zero
       if (is_bc_wall(1,1)) then
          do i1=i1min,i1max
             Krhou(1,i1,:)=0.0_wp
             Krhov(1,i1,:)=0.0_wp
             Krhow(1,i1,:)=0.0_wp
             Krhoe(1,i1,:)=0.0_wp
          enddo
       endif
       if (is_bc_wall(1,2)) then
          do i1=i1min,i1max
             Krhou(nx,i1,:)=0.0_wp
             Krhov(nx,i1,:)=0.0_wp
             Krhow(nx,i1,:)=0.0_wp
             Krhoe(nx,i1,:)=0.0_wp
          enddo
       endif
       if (is_bc_wall(2,1)) then
          do i2=i2min,i2max
             Krhou(i2,1,:)=0.0_wp
             Krhov(i2,1,:)=0.0_wp
             Krhow(i2,1,:)=0.0_wp
             Krhoe(i2,1,:)=0.0_wp
          enddo
       endif
       if (is_bc_wall(2,2)) then
          do i2=i2min,i2max
             Krhou(i2,ny,:)=0.0_wp
             Krhov(i2,ny,:)=0.0_wp
             Krhow(i2,ny,:)=0.0_wp
             Krhoe(i2,ny,:)=0.0_wp
          enddo
       endif
       if (is_bc_wall(3,1)) then
          do i2=i2min,i2max
             Krhou(i2,:,1)=0.0_wp
             Krhov(i2,:,1)=0.0_wp
             Krhow(i2,:,1)=0.0_wp
             Krhoe(i2,:,1)=0.0_wp
          enddo
       endif
       if (is_bc_wall(3,2)) then
          do i2=i2min,i2max
             Krhou(i2,:,nz)=0.0_wp
             Krhov(i2,:,nz)=0.0_wp
             Krhow(i2,:,nz)=0.0_wp
             Krhoe(i2,:,nz)=0.0_wp
          enddo
       endif
    endif

    ! One-sided communications of increments in k-direction
    ! =====================================================
!!$    ! open target windows
!!$    call window_comm1_inc
!!$    ! communicate var
!!$    call communication_1_inc_k
!    call communication1_inc_k

    ! Calculation of the local CFL number at ind-1/2
    ! ==============================================
    if (is_dtlocal) then
       if (is_curv3) then
          do k=1,nz
             do j=1,ny
                do i=1,nx
                   vcmh=sqrt((uu(i,j,k-1)*ksi_x(i,j,k-1)+vv(i,j,k-1)*ksi_y(i,j,k-1)+ww(i,j,k-1)*ksi_z(i,j,k-1))**2 &
                            +(uu(i,j,k-1)*eta_x(i,j,k-1)+vv(i,j,k-1)*eta_y(i,j,k-1)+ww(i,j,k-1)*eta_z(i,j,k-1))**2 &
                            +(uu(i,j,k-1)*phi_x(i,j,k-1)+vv(i,j,k-1)*phi_y(i,j,k-1)+ww(i,j,k-1)*phi_z(i,j,k-1))**2 &
                            )*ijacob3(i,j,k-1)
                   vcph=sqrt((uu(i,j,k)*ksi_x(i,j,k)+vv(i,j,k)*ksi_y(i,j,k)+ww(i,j,k)*ksi_z(i,j,k))**2 &
                            +(uu(i,j,k)*eta_x(i,j,k)+vv(i,j,k)*eta_y(i,j,k)+ww(i,j,k)*eta_z(i,j,k))**2 &
                            +(uu(i,j,k)*phi_x(i,j,k)+vv(i,j,k)*phi_y(i,j,k)+ww(i,j,k)*phi_z(i,j,k))**2 &
                            )*ijacob3(i,j,k)
                   rspecmh=abs(vcmh)+c_(i,j,k-1)*sqrt(g3_ksi(i,j,k-1)**2+g3_eta(i,j,k-1)**2+g3_phi(i,j,k-1)**2)
                   rspecph=abs(vcph)+c_(i,j,k  )*sqrt(g3_ksi(i,j,k  )**2+g3_eta(i,j,k  )**2+g3_phi(i,j,k  )**2)
                   cfl_l(i,j,k)=0.5_wp*(rspecmh+rspecph)*dt_local(i,j,k)
                enddo
             enddo
          enddo
       else
          do k=1,nz
             do j=1,ny
                do i=1,nx
                   rspecmh=(abs(ww(i,j,k-1))+c_(i,j,k-1))*idz(k-1)
                   rspecph=(abs(ww(i,j,k  ))+c_(i,j,k  ))*idz(k  )
                   cfl_l(i,j,k)=0.5_wp*(rspecmh+rspecph)*dt_local(i,j,k)
                enddo
             enddo
          enddo
       endif
    else
       if (is_curv3) then
          do k=1,nz
             do j=1,ny
                do i=1,nx
                   vcmh=sqrt((uu(i,j,k-1)*ksi_x(i,j,k-1)+vv(i,j,k-1)*ksi_y(i,j,k-1)+ww(i,j,k-1)*ksi_z(i,j,k-1))**2 &
                            +(uu(i,j,k-1)*eta_x(i,j,k-1)+vv(i,j,k-1)*eta_y(i,j,k-1)+ww(i,j,k-1)*eta_z(i,j,k-1))**2 &
                            +(uu(i,j,k-1)*phi_x(i,j,k-1)+vv(i,j,k-1)*phi_y(i,j,k-1)+ww(i,j,k-1)*phi_z(i,j,k-1))**2 &
                            )*ijacob3(i,j,k-1)
                   vcph=sqrt((uu(i,j,k)*ksi_x(i,j,k)+vv(i,j,k)*ksi_y(i,j,k)+ww(i,j,k)*ksi_z(i,j,k))**2 &
                            +(uu(i,j,k)*eta_x(i,j,k)+vv(i,j,k)*eta_y(i,j,k)+ww(i,j,k)*eta_z(i,j,k))**2 &
                            +(uu(i,j,k)*phi_x(i,j,k)+vv(i,j,k)*phi_y(i,j,k)+ww(i,j,k)*phi_z(i,j,k))**2 &
                            )*ijacob3(i,j,k)
                   rspecmh=abs(vcmh)+c_(i,j,k-1)*sqrt(g3_ksi(i,j,k-1)**2+g3_eta(i,j,k-1)**2+g3_phi(i,j,k-1)**2)
                   rspecph=abs(vcph)+c_(i,j,k  )*sqrt(g3_ksi(i,j,k  )**2+g3_eta(i,j,k  )**2+g3_phi(i,j,k  )**2)
                   cfl_l(i,j,k)=0.5_wp*(rspecmh+rspecph)*deltat
                enddo
             enddo
          enddo
       else
          do k=1,nz
             do j=1,ny
                do i=1,nx
                   rspecmh=(abs(ww(i,j,k-1))+c_(i,j,k-1))*idz(k-1)
                   rspecph=(abs(ww(i,j,k  ))+c_(i,j,k  ))*idz(k  )
                   cfl_l(i,j,k)=0.5_wp*(rspecmh+rspecph)*deltat
                enddo
             enddo
          enddo
       endif
    endif

!!$    ! Check that communications are fullfilled
!!$    ! ----------------------------------------
!!$    call window_comm1_inc

    ! One-sided communications of CFL in i-direction
    ! ==============================================
    ! open target windows
!    call window_comm1_cfl
    ! communicate var
!    call communication_1_cfl_k

    ! Check that communications are fullfilled
    ! ----------------------------------------
!    call window_comm1_cfl
    
    ! Communication of ghost points
    ! =============================
    !call communication_inc_k(Krho,Krhou,Krhov,Krhow,Krhoe,cfl_l)
    call communication_inc(Krho,Krhou,Krhov,Krhow,Krhoe,cfl_l)

    do i1=i1min,i1max

       ! Construction of the pentadiagonal matrix
       ! ========================================
       i3=i3min
       do i2=i2min,i2max
          b_k(i2,i3)=0.0_wp
          a_k(i2,i3)=-theta_irs1*cfl_l(i2,i1,i3+1)
          d_k(i2,i3)=1.0_wp-a_k(i2,i3)
       enddo

       i3=i3max
       do i2=i2min,i2max
          c_k(i2,i3-1)=-theta_irs1*cfl_l(i2,i1,i3)
          d_k(i2,i3)  =1.0_wp-c_k(i2,i3-1)
          e_k(i2,i3-2)=0.0_wp
       enddo

       i3=i3min+1
       do i2=i2min,i2max
          b_k(i2,i3)  =0.0_wp
          a_k(i2,i3)  =-theta_irs2*cfl_l(i2,i1,i3+1)**2
          c_k(i2,i3-1)=-theta_irs2*cfl_l(i2,i1,i3  )**2
          d_k(i2,i3)  =1.0_wp-(a_k(i2,i3)+c_k(i2,i3-1))
       enddo

       i3=i3max-1
       do i2=i2min,i2max
          a_k(i2,i3)  =-theta_irs2*cfl_l(i2,i1,i3+1)**2
          c_k(i2,i3-1)=-theta_irs2*cfl_l(i2,i1,i3  )**2
          d_k(i2,i3)  =1.0_wp-(a_k(i2,i3)+c_k(i2,i3-1))
          e_k(i2,i3-2)=0.0_wp
       enddo

       do i3=i3min+2,i3max-2
          do i2=i2min,i2max
             e_k(i2,i3-2)=theta_irs4*cfl_l(i2,i1,i3  )**4
             b_k(i2,i3)  =theta_irs4*cfl_l(i2,i1,i3+1)**4
             c_k(i2,i3-1)=-3.0_wp*e_k(i2,i3-2)-b_k(i2,i3)
             a_k(i2,i3)  =-3.0_wp*b_k(i2,i3)-e_k(i2,i3-2)
             d_k(i2,i3)  =1.0_wp+3.0_wp*(e_k(i2,i3-2)+b_k(i2,i3))
          enddo
       enddo

       ! Filling of the RHS
       do i3=i3min,i3max
          do i2=i2min,i2max
             RHS_k(i2,i3)= Krho(i2,i1,i3)
             RHSu_k(i2,i3)=Krhou(i2,i1,i3)
             RHSv_k(i2,i3)=Krhov(i2,i1,i3)
             RHSw_k(i2,i3)=Krhow(i2,i1,i3)
             RHSe_k(i2,i3)=Krhoe(i2,i1,i3)
          enddo
       enddo

       ! Resolution of the pentadiagonal linear system
       ! =============================================
       i3=i3min
       do i2=i2min,i2max
          temp_k1(i2)=1.0_wp/d_k(i2,i3)
          a_k(i2,i3)=a_k(i2,i3)*temp_k1(i2)
          b_k(i2,i3)=b_k(i2,i3)*temp_k1(i2)
          RHS_k(i2,i3)= RHS_k(i2,i3)*temp_k1(i2)
          RHSu_k(i2,i3)=RHSu_k(i2,i3)*temp_k1(i2)
          RHSv_k(i2,i3)=RHSv_k(i2,i3)*temp_k1(i2)
          RHSw_k(i2,i3)=RHSw_k(i2,i3)*temp_k1(i2)
          RHSe_k(i2,i3)=RHSe_k(i2,i3)*temp_k1(i2)
       enddo

       i3=i3min+1
       do i2=i2min,i2max
          temp_k2(i2)=c_k(i2,i3-1)
          temp_k1(i2)=1.0_wp/(d_k(i2,i3)-a_k(i2,i3-1)*temp_k2(i2))
          a_k(i2,i3)=(a_k(i2,i3)-b_k(i2,i3-1)*temp_k2(i2))*temp_k1(i2)
          b_k(i2,i3)=b_k(i2,i3)*temp_k1(i2)
          RHS_k(i2,i3)=( RHS_k(i2,i3)- RHS_k(i2,i3-1)*temp_k2(i2))*temp_k1(i2)
          RHSu_k(i2,i3)=(RHSu_k(i2,i3)-RHSu_k(i2,i3-1)*temp_k2(i2))*temp_k1(i2)
          RHSv_k(i2,i3)=(RHSv_k(i2,i3)-RHSv_k(i2,i3-1)*temp_k2(i2))*temp_k1(i2)
          RHSw_k(i2,i3)=(RHSw_k(i2,i3)-RHSw_k(i2,i3-1)*temp_k2(i2))*temp_k1(i2)
          RHSe_k(i2,i3)=(RHSe_k(i2,i3)-RHSe_k(i2,i3-1)*temp_k2(i2))*temp_k1(i2)
       enddo

       ! Step 5: for i=3,...,n-2
       do i3=i3min+2,i3max-2
          do i2=i2min,i2max
             temp_k2(i2)=c_k(i2,i3-1)-a_k(i2,i3-2)*e_k(i2,i3-2)
             temp_k1(i2)=1.0_wp/(d_k(i2,i3)-b_k(i2,i3-2)*e_k(i2,i3-2)-a_k(i2,i3-1)*temp_k2(i2))
             a_k(i2,i3)=(a_k(i2,i3)-b_k(i2,i3-1)*temp_k2(i2))*temp_k1(i2)
             b_k(i2,i3)=b_k(i2,i3)*temp_k1(i2)
             RHS_k(i2,i3) =( RHS_k(i2,i3)- RHS_k(i2,i3-2)*e_k(i2,i3-2)- RHS_k(i2,i3-1)*temp_k2(i2))*temp_k1(i2)
             RHSu_k(i2,i3)=(RHSu_k(i2,i3)-RHSu_k(i2,i3-2)*e_k(i2,i3-2)-RHSu_k(i2,i3-1)*temp_k2(i2))*temp_k1(i2)
             RHSv_k(i2,i3)=(RHSv_k(i2,i3)-RHSv_k(i2,i3-2)*e_k(i2,i3-2)-RHSv_k(i2,i3-1)*temp_k2(i2))*temp_k1(i2)
             RHSw_k(i2,i3)=(RHSw_k(i2,i3)-RHSw_k(i2,i3-2)*e_k(i2,i3-2)-RHSw_k(i2,i3-1)*temp_k2(i2))*temp_k1(i2)
             RHSe_k(i2,i3)=(RHSe_k(i2,i3)-RHSe_k(i2,i3-2)*e_k(i2,i3-2)-RHSe_k(i2,i3-1)*temp_k2(i2))*temp_k1(i2)
          enddo
       enddo

       ! Step 5: for i=n-1, n
       i3=i3max-1
       do i2=i2min,i2max
          temp_k2(i2)=c_k(i2,i3-1)-a_k(i2,i3-2)*e_k(i2,i3-2)
          temp_k1(i2)=1.0_wp/(d_k(i2,i3)-b_k(i2,i3-2)*e_k(i2,i3-2)-a_k(i2,i3-1)*temp_k2(i2))
          a_k(i2,i3)=(a_k(i2,i3)-b_k(i2,i3-1)*temp_k2(i2))*temp_k1(i2)
          RHS_k(i2,i3)=( RHS_k(i2,i3)- RHS_k(i2,i3-2)*e_k(i2,i3-2)- RHS_k(i2,i3-1)*temp_k2(i2))*temp_k1(i2)
          RHSu_k(i2,i3)=(RHSu_k(i2,i3)-RHSu_k(i2,i3-2)*e_k(i2,i3-2)-RHSu_k(i2,i3-1)*temp_k2(i2))*temp_k1(i2)
          RHSv_k(i2,i3)=(RHSv_k(i2,i3)-RHSv_k(i2,i3-2)*e_k(i2,i3-2)-RHSv_k(i2,i3-1)*temp_k2(i2))*temp_k1(i2)
          RHSw_k(i2,i3)=(RHSw_k(i2,i3)-RHSw_k(i2,i3-2)*e_k(i2,i3-2)-RHSw_k(i2,i3-1)*temp_k2(i2))*temp_k1(i2)
          RHSe_k(i2,i3)=(RHSe_k(i2,i3)-RHSe_k(i2,i3-2)*e_k(i2,i3-2)-RHSe_k(i2,i3-1)*temp_k2(i2))*temp_k1(i2)
       enddo

       i3=i3max
       do i2=i2min,i2max
          temp_k2(i2)=c_k(i2,i3-1)-a_k(i2,i3-2)*e_k(i2,i3-2)
          temp_k1(i2)=1.0_wp/(d_k(i2,i3)-b_k(i2,i3-2)*e_k(i2,i3-2)-a_k(i2,i3-1)*temp_k2(i2))
          RHS_k(i2,i3)=( RHS_k(i2,i3)- RHS_k(i2,i3-2)*e_k(i2,i3-2)- RHS_k(i2,i3-1)*temp_k2(i2))*temp_k1(i2)
          RHSu_k(i2,i3)=(RHSu_k(i2,i3)-RHSu_k(i2,i3-2)*e_k(i2,i3-2)-RHSu_k(i2,i3-1)*temp_k2(i2))*temp_k1(i2)
          RHSv_k(i2,i3)=(RHSv_k(i2,i3)-RHSv_k(i2,i3-2)*e_k(i2,i3-2)-RHSv_k(i2,i3-1)*temp_k2(i2))*temp_k1(i2)
          RHSw_k(i2,i3)=(RHSw_k(i2,i3)-RHSw_k(i2,i3-2)*e_k(i2,i3-2)-RHSw_k(i2,i3-1)*temp_k2(i2))*temp_k1(i2)
          RHSe_k(i2,i3)=(RHSe_k(i2,i3)-RHSe_k(i2,i3-2)*e_k(i2,i3-2)-RHSe_k(i2,i3-1)*temp_k2(i2))*temp_k1(i2)
       enddo


       ! Step 6: computation of the solution vector
       i3=i3max-1
       do i2=i2min,i2max
          RHS_k(i2,i3)= RHS_k(i2,i3)-a_k(i2,i3)* RHS_k(i2,i3+1)
          RHSu_k(i2,i3)=RHSu_k(i2,i3)-a_k(i2,i3)*RHSu_k(i2,i3+1)
          RHSv_k(i2,i3)=RHSv_k(i2,i3)-a_k(i2,i3)*RHSv_k(i2,i3+1)
          RHSw_k(i2,i3)=RHSw_k(i2,i3)-a_k(i2,i3)*RHSw_k(i2,i3+1)
          RHSe_k(i2,i3)=RHSe_k(i2,i3)-a_k(i2,i3)*RHSe_k(i2,i3+1)
       enddo

       do i3=i3max-2,i3min,-1
          do i2=i2min,i2max
             RHS_k(i2,i3)= RHS_k(i2,i3)-a_k(i2,i3)* RHS_k(i2,i3+1)-b_k(i2,i3)* RHS_k(i2,i3+2)
             RHSu_k(i2,i3)=RHSu_k(i2,i3)-a_k(i2,i3)*RHSu_k(i2,i3+1)-b_k(i2,i3)*RHSu_k(i2,i3+2)
             RHSv_k(i2,i3)=RHSv_k(i2,i3)-a_k(i2,i3)*RHSv_k(i2,i3+1)-b_k(i2,i3)*RHSv_k(i2,i3+2)
             RHSw_k(i2,i3)=RHSw_k(i2,i3)-a_k(i2,i3)*RHSw_k(i2,i3+1)-b_k(i2,i3)*RHSw_k(i2,i3+2)
             RHSe_k(i2,i3)=RHSe_k(i2,i3)-a_k(i2,i3)*RHSe_k(i2,i3+1)-b_k(i2,i3)*RHSe_k(i2,i3+2)
          enddo
       enddo

       ! Update with the solution
       do i3=i3min,i3max
          do i2=i2min,i2max
             Krho(i2,i1,i3)= RHS_k(i2,i3)
             Krhou(i2,i1,i3)=RHSu_k(i2,i3)
             Krhov(i2,i1,i3)=RHSv_k(i2,i3)
             Krhow(i2,i1,i3)=RHSw_k(i2,i3)
             Krhoe(i2,i1,i3)=RHSe_k(i2,i3)
          enddo
       enddo
       
    enddo

  end subroutine irs4_ngh_k
  
end module mod_irs_d
