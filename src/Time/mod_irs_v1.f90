!=================================================================================
module mod_irs_v1
!=================================================================================
  !> author: Aurelien Bienner, modifs XG 
  !> date: 2021
  !> Module for Implicit Residual Smoothing, parallelized with ghost cells
!=================================================================================
  implicit none
  ! ------------------------------------------------------------------------------
  ! ------------------------------------------------------------------------------

contains

  !===============================================================================
  subroutine irs2_ngh_v1
  !===============================================================================
    !> 2nd-order Implicit Residual Smoothing (IRS2)
    !> Parallel inversion of tridiagonal matrix using ghost points
  !===============================================================================
    use mod_time      ! <- for deltat,is_irs
    use mod_flow      ! <- for velocities,increments,metrics
    use mod_interface ! <- for communication
    implicit none
    ! ----------------------------------------------------------------------------
    ! indices
    integer :: i1,i2,i3
    integer :: i1min,i1max,i2min,i2max,i3min,i3max
    ! spectral radius
    real(wp) :: rspecmh,rspecph
    real(wp), dimension(ndx-ngh:nfx+ngh,ndy-ngh:nfy+ngh,ndz-ngh:nfz+ngh) :: r_spec2
    ! tridiagonal in i-direction
    ! --------------------------
    ! matrix elements a,d,c
    real(wp), dimension(1:ny,ndx-ngh:nfx+ngh)   :: d_i
    real(wp), dimension(1:ny,ndx-ngh:nfx+ngh-1) :: a_i,c_i
    ! right-hand side (RHS)
    real(wp), dimension(1:ny,ndx-ngh:nfx+ngh) :: RHS_i,RHSu_i,RHSv_i,RHSw_i,RHSe_i
    ! work arrays
    real(wp), dimension(ny) :: temp_i
    ! tridiagonal in j-direction
    ! --------------------------
    ! matrix elements a,d,c
    real(wp), dimension(1:nx,ndy-ngh:nfy+ngh)   :: d_j
    real(wp), dimension(1:nx,ndy-ngh:nfy+ngh-1) :: a_j,c_j
    ! work arrays
    real(wp), dimension(nx) :: temp_j
    ! tridiagonal in k-direction
    ! --------------------------
    ! matrix elements a,d,c
    real(wp), dimension(1:nx,ndz-ngh:nfz+ngh)   :: d_k
    real(wp), dimension(1:nx,ndz-ngh:nfz+ngh-1) :: a_k,c_k
    ! right-hand side (RHS)
    real(wp), dimension(1:nx,ndz-ngh:nfz+ngh) :: RHS_k,RHSu_k,RHSv_k,RHSw_k,RHSe_k
    ! work arrays
    real(wp), dimension(nx) :: temp_k
    ! ----------------------------------------------------------------------------
 
    !----------------------------------------
    ! Implicitation of i-direction
    !----------------------------------------
    if (is_irs_i) then
       !  i3 : i;   i1 : k;   i2 : j
       i1min=1;   i2min=1;   i3min=ndx-ngh
       i1max=nz;  i2max=ny;  i3max=nfx+ngh

       ! If wall, do not take into account the wall point
       ! ================================================
       if (is_bc_wall(2,1)) i2min=2
       if (is_bc_wall(2,2)) i2max=ny-1
       if (is_bc_wall(3,1)) i1min=2
       if (is_bc_wall(3,2)) i1max=nz-1

       ! Communication of ghost points
       ! =============================
       call communication_(Krho,Krhou,Krhov,Krhow,Krhoe)

       ! Calculation of the spectral radius at ind-1/2
       ! =============================================
       do i1=i1min,i1max
          do i2=i2min,i2max
             do i3=i3min+1,i3max
                rspecmh=(abs(uu(i3-1,i2,i1))+c_(i3-1,i2,i1))*idx(i3-1)
                rspecph=(abs(uu(i3  ,i2,i1))+c_(i3  ,i2,i1))*idx(i3)
                r_spec2(i3,i2,i1)=(0.5_wp*(rspecmh+rspecph)*deltat)**2
             enddo
          enddo
       enddo

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
             a_i(i2,i3)=-theta_irs1*r_spec2(i3+1,i2,i1)**0.5
             d_i(i2,i3)=1.0_wp-a_i(i2,i3)
          enddo

          i3=i3max
          do i2=i2min,i2max
             c_i(i2,i3-1)=-theta_irs1*r_spec2(i3,i2,i1)**0.5
             d_i(i2,i3)  =1.0_wp-c_i(i2,i3-1)
          enddo

          do i3=i3min+1,i3max-1
             do i2=i2min, i2max
                c_i(i2,i3-1)=-theta_irs2*r_spec2(i3,i2,i1)
                a_i(i2,i3)  =-theta_irs2*r_spec2(i3+1,i2,i1)
                d_i(i2,i3)  =1.0_wp-(c_i(i2,i3-1)+a_i(i2,i3))
             enddo
          enddo

          ! Filling of the RHS
          do i3=i3min,i3max
             do i2=i2min,i2max
                RHS_i(i2,i3)=  Krho(i3,i2,i1)
                RHSu_i(i2,i3)= Krhou(i3,i2,i1)
                RHSv_i(i2,i3)= Krhov(i3,i2,i1)
                RHSw_i(i2,i3)= Krhow(i3,i2,i1)
                RHSe_i(i2,i3)= Krhoe(i3,i2,i1)
             enddo
          enddo

          ! Resolution of the tridiagonal linear system
          ! ===========================================
          ! Thomas' algorithm: Forward elimination phase
          do i2=i2min,i2max
             temp_i(i2)=1.0_wp/d_i(i2,i3min)
             a_i(i2,i3min)=a_i(i2,i3min)*temp_i(i2)
             RHS_i(i2,i3min) = RHS_i(i2,i3min) *temp_i(i2)
             RHSu_i(i2,i3min)= RHSu_i(i2,i3min)*temp_i(i2)
             RHSv_i(i2,i3min)= RHSv_i(i2,i3min)*temp_i(i2)
             RHSw_i(i2,i3min)= RHSw_i(i2,i3min)*temp_i(i2)
             RHSe_i(i2,i3min)= RHSe_i(i2,i3min)*temp_i(i2)
          enddo

          do i3=i3min+1,i3max-1
             do i2=i2min,i2max
                temp_i(i2)=1.0_wp/(d_i(i2,i3)-c_i(i2,i3-1)*a_i(i2,i3-1))
                a_i(i2,i3)=a_i(i2,i3)*temp_i(i2)
                RHS_i(i2,i3) =( RHS_i(i2,i3)-c_i(i2,i3-1)* RHS_i(i2,i3-1))*temp_i(i2)
                RHSu_i(i2,i3)=(RHSu_i(i2,i3)-c_i(i2,i3-1)*RHSu_i(i2,i3-1))*temp_i(i2)
                RHSv_i(i2,i3)=(RHSv_i(i2,i3)-c_i(i2,i3-1)*RHSv_i(i2,i3-1))*temp_i(i2)
                RHSw_i(i2,i3)=(RHSw_i(i2,i3)-c_i(i2,i3-1)*RHSw_i(i2,i3-1))*temp_i(i2)
                RHSe_i(i2,i3)=(RHSe_i(i2,i3)-c_i(i2,i3-1)*RHSe_i(i2,i3-1))*temp_i(i2)
             enddo
          enddo

          ! Backward substitution phase
          do i2=i2min,i2max
             temp_i(i2)=1.0_wp/(d_i(i2,i3max)-c_i(i2,i3max-1)*a_i(i2,i3max-1))
             RHS_i(i2,i3max) =( RHS_i(i2,i3max)-c_i(i2,i3max-1)* RHS_i(i2,i3max-1))*temp_i(i2)
             RHSu_i(i2,i3max)=(RHSu_i(i2,i3max)-c_i(i2,i3max-1)*RHSu_i(i2,i3max-1))*temp_i(i2)
             RHSv_i(i2,i3max)=(RHSv_i(i2,i3max)-c_i(i2,i3max-1)*RHSv_i(i2,i3max-1))*temp_i(i2)
             RHSw_i(i2,i3max)=(RHSw_i(i2,i3max)-c_i(i2,i3max-1)*RHSw_i(i2,i3max-1))*temp_i(i2)
             RHSe_i(i2,i3max)=(RHSe_i(i2,i3max)-c_i(i2,i3max-1)*RHSe_i(i2,i3max-1))*temp_i(i2)
          enddo

          do i3=i3max-1,i3min,-1
             do i2=i2min,i2max
                RHS_i(i2,i3) = RHS_i(i2,i3)-a_i(i2,i3)* RHS_i(i2,i3+1)
                RHSu_i(i2,i3)=RHSu_i(i2,i3)-a_i(i2,i3)*RHSu_i(i2,i3+1)
                RHSv_i(i2,i3)=RHSv_i(i2,i3)-a_i(i2,i3)*RHSv_i(i2,i3+1)
                RHSw_i(i2,i3)=RHSw_i(i2,i3)-a_i(i2,i3)*RHSw_i(i2,i3+1)
                RHSe_i(i2,i3)=RHSe_i(i2,i3)-a_i(i2,i3)*RHSe_i(i2,i3+1)
             enddo
          enddo

          ! Update with the solution
          do i3=i3min,i3max
             do i2=i2min,i2max
                Krho(i3,i2,i1) =  RHS_i(i2,i3)
                Krhou(i3,i2,i1)= RHSu_i(i2,i3)
                Krhov(i3,i2,i1)= RHSv_i(i2,i3)
                Krhow(i3,i2,i1)= RHSw_i(i2,i3)
                Krhoe(i3,i2,i1)= RHSe_i(i2,i3)
             enddo
          enddo
       enddo

    endif

    !----------------------------------------
    ! Implicitation of j-direction
    !----------------------------------------
    if (is_irs_j) then
       !  i3 : j;   i1 : k;   i2 : i
       i1min=1;   i2min=1;   i3min=ndy-ngh
       i1max=nz;  i2max=nx;  i3max=nfy+ngh

       ! If wall, do not take into account the wall point
       ! ================================================
       if (is_bc_wall(1,1)) i2min=2
       if (is_bc_wall(1,2)) i2max=nx-1
       if (is_bc_wall(3,1)) i1min=2
       if (is_bc_wall(3,2)) i1max=nz-1

       ! Communication of ghost points
       ! =============================
       call communication_(Krho,Krhou,Krhov,Krhow,Krhoe)

       ! Calculation of the spectral radius at ind-1/2
       ! =============================================
       do i1=i1min,i1max
          do i3=i3min+1,i3max
             do i2=i2min,i2max
                rspecmh=(abs(vv(i2,i3-1,i1))+c_(i2,i3-1,i1))*idy(i3-1)
                rspecph=(abs(vv(i2,i3  ,i1))+c_(i2,i3  ,i1))*idy(i3)
                r_spec2(i2,i3,i1)=(0.5_wp*(rspecmh+rspecph)*deltat)**2
             enddo
          enddo
       enddo
          
       do i1=i1min,i1max
          
          ! Fill tridiagonal matrix elements (RHS=Krho)
          ! ================================          
          i3=i3min
          do i2=i2min,i2max
             a_j(i2,i3)=-theta_irs1*(r_spec2(i2,i3+1,i1))**0.5
             d_j(i2,i3)=1.0_wp-a_j(i2,i3)
          enddo

          i3=i3max
          do i2=i2min,i2max
             c_j(i2,i3-1)=-theta_irs1*(r_spec2(i2,i3,i1))**0.5
             d_j(i2,i3)  =1.0_wp-c_j(i2,i3-1)
          enddo

          do i3=i3min+1,i3max-1
             do i2=i2min,i2max
                c_j(i2,i3-1)=-theta_irs2*r_spec2(i2,i3,i1)
                a_j(i2,i3)  =-theta_irs2*r_spec2(i2,i3+1,i1)
                d_j(i2,i3)  =1.0_wp-(c_j(i2,i3-1)+a_j(i2,i3))
             enddo
          enddo

          ! Resolution of the tridiagonal linear system
          ! ===========================================
          ! Thomas' algorithm: Forward elimination phase
          do i2=i2min,i2max
             temp_j(i2)=1.0_wp/d_j(i2,i3min)
             a_j(i2,i3min)=a_j(i2,i3min)*temp_j(i2)
             Krho(i2,i3min,i1)=  Krho(i2,i3min,i1) *temp_j(i2)
             Krhou(i2,i3min,i1)= Krhou(i2,i3min,i1)*temp_j(i2)
             Krhov(i2,i3min,i1)= Krhov(i2,i3min,i1)*temp_j(i2)
             Krhow(i2,i3min,i1)= Krhow(i2,i3min,i1)*temp_j(i2)
             Krhoe(i2,i3min,i1)= Krhoe(i2,i3min,i1)*temp_j(i2)
          enddo

          do i3=i3min+1,i3max-1
             do i2=i2min,i2max
                temp_j(i2)=1.0_wp/(d_j(i2,i3)-c_j(i2,i3-1)*a_j(i2,i3-1))
                a_j(i2,i3)=a_j(i2,i3)*temp_j(i2)
                Krho(i2,i3,i1) =( Krho(i2,i3,i1)-c_j(i2,i3-1)* Krho(i2,i3-1,i1))*temp_j(i2)
                Krhou(i2,i3,i1)=(Krhou(i2,i3,i1)-c_j(i2,i3-1)*Krhou(i2,i3-1,i1))*temp_j(i2)
                Krhov(i2,i3,i1)=(Krhov(i2,i3,i1)-c_j(i2,i3-1)*Krhov(i2,i3-1,i1))*temp_j(i2)
                Krhow(i2,i3,i1)=(Krhow(i2,i3,i1)-c_j(i2,i3-1)*Krhow(i2,i3-1,i1))*temp_j(i2)
                Krhoe(i2,i3,i1)=(Krhoe(i2,i3,i1)-c_j(i2,i3-1)*Krhoe(i2,i3-1,i1))*temp_j(i2)
             enddo
          enddo

          ! Backward substitution phase
          do i2=i2min,i2max
             temp_j(i2)=1.0_wp/(d_j(i2,i3max)-c_j(i2,i3max-1)*a_j(i2,i3max-1))
             Krho(i2,i3max,i1) =( Krho(i2,i3max,i1)-c_j(i2,i3max-1)* Krho(i2,i3max-1,i1))*temp_j(i2)
             Krhou(i2,i3max,i1)=(Krhou(i2,i3max,i1)-c_j(i2,i3max-1)*Krhou(i2,i3max-1,i1))*temp_j(i2)
             Krhov(i2,i3max,i1)=(Krhov(i2,i3max,i1)-c_j(i2,i3max-1)*Krhov(i2,i3max-1,i1))*temp_j(i2)
             Krhow(i2,i3max,i1)=(Krhow(i2,i3max,i1)-c_j(i2,i3max-1)*Krhow(i2,i3max-1,i1))*temp_j(i2)
             Krhoe(i2,i3max,i1)=(Krhoe(i2,i3max,i1)-c_j(i2,i3max-1)*Krhoe(i2,i3max-1,i1))*temp_j(i2)
          enddo

          do i3=i3max-1,i3min,-1
             do i2=i2min,i2max
                Krho(i2,i3,i1) = Krho(i2,i3,i1)-a_j(i2,i3)* Krho(i2,i3+1,i1)
                Krhou(i2,i3,i1)=Krhou(i2,i3,i1)-a_j(i2,i3)*Krhou(i2,i3+1,i1)
                Krhov(i2,i3,i1)=Krhov(i2,i3,i1)-a_j(i2,i3)*Krhov(i2,i3+1,i1)
                Krhow(i2,i3,i1)=Krhow(i2,i3,i1)-a_j(i2,i3)*Krhow(i2,i3+1,i1)
                Krhoe(i2,i3,i1)=Krhoe(i2,i3,i1)-a_j(i2,i3)*Krhoe(i2,i3+1,i1)
             enddo
          enddo
       enddo

    endif

    !----------------------------------------
    ! Implicitation of k-direction
    !----------------------------------------
    if (is_irs_k) then
       !  i3 : k;   i1 : j;   i2 : i
       i1min=1;   i2min=1;   i3min=ndz-ngh
       i1max=ny;  i2max=nx;  i3max=nfz+ngh

       ! If wall, do not take into account the wall point
       ! ================================================
       if (is_bc_wall(1,1)) i2min=2
       if (is_bc_wall(1,2)) i2max=nx-1
       if (is_bc_wall(2,1)) i1min=2
       if (is_bc_wall(2,2)) i1max=ny-1

       ! Communication of ghost points
       ! =============================
       call communication_(Krho,Krhou,Krhov,Krhow,Krhoe)

       ! Calculation of the spectral radius at ind-1/2
       ! =============================================
       do i3=i3min+1,i3max
          do i1=i1min,i1max
             do i2=i2min,i2max
                rspecmh=(abs(ww(i2,i1,i3-1))+c_(i2,i1,i3-1))*idz(i3-1)
                rspecph=(abs(ww(i2,i1,i3  ))+c_(i2,i1,i3  ))*idz(i3)
                r_spec2(i2,i1,i3)=(0.5_wp*(rspecmh+rspecph)*deltat)**2
             enddo
          enddo
       enddo

       ! Fill tridiagonal matrix elements & RHS
       ! ======================================
       do i1=i1min,i1max
          i3=i3min
          do i2=i2min,i2max
             a_k(i2,i3)=-theta_irs1*(r_spec2(i2,i1,i3+1))**0.5
             d_k(i2,i3)=1.0_wp-a_k(i2,i3)
          enddo

          i3=i3max
          do i2=i2min,i2max
             c_k(i2,i3-1)=-theta_irs1*(r_spec2(i2,i1,i3))**0.5
             d_k(i2,i3)  =1.0_wp-c_k(i2,i3-1)
          enddo

          do i3=i3min+1,i3max-1
             do i2=i2min, i2max
                c_k(i2,i3-1)=-theta_irs2*r_spec2(i2,i1,i3)
                a_k(i2,i3)  =-theta_irs2*r_spec2(i2,i1,i3+1)
                d_k(i2,i3)  =1.0_wp-(c_k(i2,i3-1)+a_k(i2,i3))
             enddo
          enddo

          ! Filling of the RHS
          do i3=i3min,i3max
             do i2=i2min,i2max
                RHS_k(i2,i3) = Krho(i2,i1,i3)
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
             RHS_k(i2,i3min) = RHS_k(i2,i3min)*temp_k(i2)
             RHSu_k(i2,i3min)=RHSu_k(i2,i3min)*temp_k(i2)
             RHSv_k(i2,i3min)=RHSv_k(i2,i3min)*temp_k(i2)
             RHSw_k(i2,i3min)=RHSw_k(i2,i3min)*temp_k(i2)
             RHSe_k(i2,i3min)=RHSe_k(i2,i3min)*temp_k(i2)
          enddo

          do i3=i3min+1,i3max-1
             do i2=i2min,i2max
                temp_k(i2)=1.0_wp/(d_k(i2,i3)-c_k(i2,i3-1)*a_k(i2,i3-1))
                a_k(i2,i3)=a_k(i2,i3)*temp_k(i2)
                RHS_k(i2,i3) =( RHS_k(i2,i3)-c_k(i2,i3-1)* RHS_k(i2,i3-1))*temp_k(i2)
                RHSu_k(i2,i3)=(RHSu_k(i2,i3)-c_k(i2,i3-1)*RHSu_k(i2,i3-1))*temp_k(i2)
                RHSv_k(i2,i3)=(RHSv_k(i2,i3)-c_k(i2,i3-1)*RHSv_k(i2,i3-1))*temp_k(i2)
                RHSw_k(i2,i3)=(RHSw_k(i2,i3)-c_k(i2,i3-1)*RHSw_k(i2,i3-1))*temp_k(i2)
                RHSe_k(i2,i3)=(RHSe_k(i2,i3)-c_k(i2,i3-1)*RHSe_k(i2,i3-1))*temp_k(i2)
             enddo
          enddo

          ! Backward substitution phase
          do i2=i2min,i2max
             temp_k(i2)=1.0_wp/(d_k(i2,i3max)-c_k(i2,i3max-1)*a_k(i2,i3max-1))
             RHS_k(i2,i3max) =( RHS_k(i2,i3max)-c_k(i2,i3max-1)* RHS_k(i2,i3max-1))*temp_k(i2)
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

    endif
  end subroutine irs2_ngh_v1

  !===============================================================================
  subroutine irs4_ngh_v1
  !===============================================================================
    !> 4th-order Implicit Residual Smoothing (IRS4)
    !> Parallel inversion of pentadiagonal matrix using ghost points
  !===============================================================================
    use mod_time      ! <- for deltat,is_irs
    use mod_flow      ! <- for velocities,increments,metrics
    use mod_interface ! <- for communication
    implicit none
    ! ----------------------------------------------------------------------------
    ! indices
    integer :: i1,i2,i3
    integer :: i1min,i1max,i2min,i2max,i3min,i3max
    ! spectral radius
    real(wp) :: rspecmh,rspecph
    real(wp), dimension(ndx-ngh:nfx+ngh,ndy-ngh:nfy+ngh,ndz-ngh:nfz+ngh) :: r_spec4
    ! pentadiagonal in i-direction
    ! ----------------------------
    ! matrix elements a,b,c,d,e
    real(wp), dimension(1:ny,ndx-ngh:nfx+ngh)   :: d_i
    real(wp), dimension(1:ny,ndx-ngh:nfx+ngh-1) :: a_i,c_i
    real(wp), dimension(1:ny,ndx-ngh:nfx+ngh-2) :: b_i,e_i
    ! right-hand side (RHS)
    real(wp), dimension(1:ny,ndx-ngh:nfx+ngh) :: RHS_i,RHSu_i,RHSv_i,RHSw_i,RHSe_i
    ! work arrays
    real(wp), dimension(ny) :: temp_i1,temp_i2
    ! pentadiagonal in j-direction
    ! ----------------------------
    ! matrix elements a,b,c,d,e
    real(wp), dimension(1:nx,ndy-ngh:nfy+ngh)   :: d_j
    real(wp), dimension(1:nx,ndy-ngh:nfy+ngh-1) :: a_j,c_j
    real(wp), dimension(1:nx,ndy-ngh:nfy+ngh-2) :: b_j,e_j
    ! work arrays
    real(wp), dimension(nx) :: temp_j1,temp_j2
    ! pentadiagonal in k-direction
    ! ----------------------------
    ! matrix elements a,b,c,d,e
    real(wp), dimension(1:nx,ndz-ngh:nfz+ngh)   :: d_k
    real(wp), dimension(1:nx,ndz-ngh:nfz+ngh-1) :: a_k,c_k
    real(wp), dimension(1:nx,ndz-ngh:nfz+ngh-2) :: b_k,e_k
    ! right-hand side (RHS)
    real(wp), dimension(1:nx,ndz-ngh:nfz+ngh) :: RHS_k,RHSu_k,RHSv_k,RHSw_k,RHSe_k
    ! work arrays
    real(wp), dimension(nx) :: temp_k1,temp_k2
    ! ----------------------------------------------------------------------------

    !----------------------------------------
    ! Implicitation of i-direction
    !----------------------------------------
    if (is_irs_i) then
       !  i3 : i;   i1 : k;   i2 : j
       i1min=1;   i2min=1;   i3min=ndx-ngh
       i1max=nz;  i2max=ny;  i3max=nfx+ngh

       ! If wall, do not take into account the wall point
       ! ================================================
       if (is_bc_wall(2,1)) i2min=2
       if (is_bc_wall(2,2)) i2max=ny-1
       if (is_bc_wall(3,1)) i1min=2
       if (is_bc_wall(3,2)) i1max=nz-1

       ! Communication of ghost points
       ! =============================
       call communication_(Krho,Krhou,Krhov,Krhow,Krhoe)

       ! Calculation of the spectral radius at ind-1/2
       ! =============================================
       do i1=i1min,i1max
          do i2=i2min,i2max
             do i3=i3min+1,i3max
                rspecmh=(abs(uu(i3-1,i2,i1))+c_(i3-1,i2,i1))*idx(i3-1)
                rspecph=(abs(uu(i3  ,i2,i1))+c_(i3  ,i2,i1))*idx(i3)
                r_spec4(i3,i2,i1)=(0.5_wp*(rspecmh+rspecph)*deltat)**4
             enddo
          enddo
       enddo

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
             a_i(i2,i3)=-theta_irs1*r_spec4(i3+1,i2,i1)**0.25
             d_i(i2,i3)=1.0_wp-a_i(i2,i3)
          enddo

          i3=i3max
          do i2=i2min,i2max
             c_i(i2,i3-1)=-theta_irs1*r_spec4(i3,i2,i1)**0.25
             d_i(i2,i3)  =1.0_wp-c_i(i2,i3-1)
             e_i(i2,i3-2)=0.0_wp
          enddo

          i3=i3min+1
          do i2=i2min,i2max
             b_i(i2,i3)  =0.0_wp
             a_i(i2,i3)  =-theta_irs2*(r_spec4(i3+1,i2,i1))**0.5
             c_i(i2,i3-1)=-theta_irs2*(r_spec4(i3,i2,i1))**0.5
             d_i(i2,i3)  =1.0_wp-(a_i(i2,i3)+c_i(i2,i3-1))
          enddo

          i3=i3max-1
          do i2=i2min,i2max
             a_i(i2,i3)  =-theta_irs2*(r_spec4(i3+1,i2,i1))**0.5
             c_i(i2,i3-1)=-theta_irs2*(r_spec4(i3,i2,i1))**0.5
             d_i(i2,i3)  =1.0_wp-(a_i(i2,i3)+c_i(i2,i3-1))
             e_i(i2,i3-2)=0.0_wp
          enddo

          do i3=i3min+2,i3max-2
             do i2=i2min,i2max
                e_i(i2,i3-2)=theta_irs4*r_spec4(i3,i2,i1)
                b_i(i2,i3)  =theta_irs4*r_spec4(i3+1,i2,i1)
                
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

    endif

    !----------------------------------------
    ! Implicitation of j-direction
    !----------------------------------------
    if (is_irs_j) then
       !  i3 : j;   i1 : k;   i2 : i
       i1min=1;   i2min=1;   i3min=ndy-ngh
       i1max=nz;  i2max=nx;  i3max=nfy+ngh

       ! If wall, do not take into account the wall point
       ! ================================================
       if (is_bc_wall(1,1)) i2min=2
       if (is_bc_wall(1,2)) i2max=nx-1
       if (is_bc_wall(3,1)) i1min=2
       if (is_bc_wall(3,2)) i1max=nz-1

       ! Communication of ghost points
       ! =============================
       call communication_(Krho,Krhou,Krhov,Krhow,Krhoe)

       ! Calculation of the spectral radius at ind-1/2
       ! =============================================
       do i1=i1min,i1max
          do i3=i3min+1,i3max
             do i2=i2min,i2max
                rspecmh=(abs(vv(i2,i3-1,i1))+c_(i2,i3-1,i1))*idy(i3-1)
                rspecph=(abs(vv(i2,i3  ,i1))+c_(i2,i3  ,i1))*idy(i3)
                r_spec4(i2,i3,i1)=(0.5_wp*(rspecmh+rspecph)*deltat)**4
             enddo
          enddo
       enddo

       do i1=i1min,i1max
          
          ! Construction of the pentadiagonal matrix (RHS=Krho)
          ! ========================================
          i3=i3min
          do i2=i2min,i2max
             b_j(i2,i3)=0.0_wp
             a_j(i2,i3)=-theta_irs1*r_spec4(i2,i3+1,i1)**0.25
             d_j(i2,i3)=1.0_wp-a_j(i2,i3)
          enddo

          i3=i3max
          do i2=i2min,i2max
             c_j(i2,i3-1)=-theta_irs1*r_spec4(i2,i3,i1)**0.25
             d_j(i2,i3)  =1.0_wp-c_j(i2,i3-1)
             e_j(i2,i3-2)=0.0_wp
          enddo

          i3=i3min+1
          do i2=i2min,i2max
             b_j(i2,i3)  =0.0_wp
             a_j(i2,i3)  =-theta_irs2*(r_spec4(i2,i3+1,i1))**0.5
             c_j(i2,i3-1)=-theta_irs2*(r_spec4(i2,i3,i1))**0.5
             d_j(i2,i3)  =1.0_wp-(a_j(i2,i3)+c_j(i2,i3-1))
          enddo

          i3=i3max-1
          do i2=i2min,i2max
             a_j(i2,i3)  =-theta_irs2*(r_spec4(i2,i3+1,i1))**0.5
             c_j(i2,i3-1)=-theta_irs2*(r_spec4(i2,i3,i1))**0.5
             d_j(i2,i3)  =1.0_wp-(a_j(i2,i3)+c_j(i2,i3-1))
             e_j(i2,i3-2)=0.0_wp
          enddo

          do i3=i3min+2,i3max-2
             do i2=i2min,i2max
                e_j(i2,i3-2)=theta_irs4*r_spec4(i2,i3,i1)
                b_j(i2,i3)  =theta_irs4*r_spec4(i2,i3+1,i1)
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
          do i2=i2min, i2max
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
             temp_j1(i2)=1.0_wp/(d_j(i2,i3)-b_j(i2,i3-2)*e_j(i2,i3-2) &
                                  -a_j(i2,i3-1)*temp_j2(i2))
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

    endif

    !----------------------------------------
    ! Implicitation of k-direction
    !----------------------------------------
    if (is_irs_k) then
       !  i3 : k;   i1 : i;   i2 : j
       i1min=1;   i2min=1;   i3min=ndz-ngh
       i1max=ny;  i2max=nx;  i3max=nfz+ngh

       ! If wall, do not take into account the wall point
       ! ================================================
       if (is_bc_wall(1,1)) i2min=2
       if (is_bc_wall(1,2)) i2max=nx-1
       if (is_bc_wall(2,1)) i1min=2
       if (is_bc_wall(2,2)) i1max=ny-1

       ! Communication of ghost points
       ! =============================
       call communication_(Krho,Krhou,Krhov,Krhow,Krhoe)

       ! Calculation of the spectral radius at ind-1/2
       ! =============================================
       do i3=i3min+1,i3max
          do i1=i1min,i1max
             do i2=i2min,i2max
                rspecmh=(abs(ww(i2,i1,i3-1)**2)+c_(i2,i1,i3-1))*idz(i3-1)
                rspecph=(abs(ww(i2,i1,i3  )**2)+c_(i2,i1,i3  ))*idz(i3)
                r_spec4(i2,i1,i3)=(0.5_wp*(rspecmh+rspecph)*deltat)**4
             enddo
          enddo
       enddo

       do i1=i1min,i1max
          
          ! Construction of the pentadiagonal matrix
          ! ========================================
          i3=i3min
          do i2=i2min,i2max
             b_k(i2,i3)=0.0_wp
             a_k(i2,i3)=-theta_irs1*(r_spec4(i2,i1,i3+1))**0.25
             d_k(i2,i3)=1.0_wp-a_k(i2,i3)
          enddo

          i3=i3max
          do i2=i2min,i2max
             c_k(i2,i3-1)=-theta_irs1*(r_spec4(i2,i1,i3))**0.25
             d_k(i2,i3)  =1.0_wp-c_k(i2,i3-1)
             e_k(i2,i3-2)=0.0_wp
          enddo

          i3=i3min+1
          do i2=i2min,i2max
             b_k(i2,i3)  =0.0_wp
             a_k(i2,i3)  =-theta_irs2*(r_spec4(i2,i1,i3+1))**0.5
             c_k(i2,i3-1)=-theta_irs2*(r_spec4(i2,i1,i3))**0.5
             d_k(i2,i3)  =1.0_wp-(a_k(i2,i3)+c_k(i2,i3-1))
          enddo

          i3=i3max-1
          do i2=i2min,i2max
             a_k(i2,i3)  =-theta_irs2*(r_spec4(i2,i1,i3+1))**0.5
             c_k(i2,i3-1)=-theta_irs2*(r_spec4(i2,i1,i3))**0.5
             d_k(i2,i3)  =1.0_wp-(a_k(i2,i3)+c_k(i2,i3-1))
             e_k(i2,i3-2)=0.0_wp
          enddo

          do i3=i3min+2,i3max-2
             do i2=i2min,i2max
                e_k(i2,i3-2)=theta_irs4*r_spec4(i2,i1,i3)
                b_k(i2,i3)  =theta_irs4*r_spec4(i2,i1,i3+1)
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

    endif

  end subroutine irs4_ngh_v1
  
end module mod_irs_v1
