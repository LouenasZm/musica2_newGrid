!===============================================================================
subroutine init_irs_c_v0
!===============================================================================
  !> Initialisation of the Implicit Residual Smoothing
!===============================================================================
   use mod_constant
   use mod_time
   use mod_interface
   use mod_mpi_part
   use mod_grid
   use mod_flow
   use warnstop
   implicit none

   if (iproc==0) print *, 'Initialisation of IRS...'

   ! Safety checks
   ! =============
   if (is_2d) is_irs_k = .false.

   ! Parallelisation with ghost points only implemented for IRS2 or IRS4
   if (type_para.ne.1 .and. type_para.ne.3) then
      call mpistop('Only ghost points for IRS2-IRS4 and Pascal for IRS2 implemented. Shutting down...',0)
   else if (iirs.ne.2 .and. iirs.ne.4) then
      call mpistop('Only ghost points for IRS2-IRS4 and Pascal for IRS2 implemented. Shutting down...',0)
   end if

   ! Definition of the necessary metrics for IRS
   ! -------------------------------------------
   call init_grid_irs_c_v0

   ! Parallelisation with ghost points
   ! =================================
   if (type_para==1) then
      if (iproc==0) print *,' Parallelisation of IRS with ghost points'
      
      if (iproc==0) print *,' Number of ghost points used for IRS : ', ngh

   ! Parallelisation with pascal for IRS2
   ! ====================================
   else if (type_para==3) then
      if (iproc==0) print *,'Parallelisation of IRS with Pascal'
      ! number of RHS
      NRHS = 5

      if (is_irs_i) allocate(DLp_ksi(1:ny,1:nx,NRHS),Dp_ksi(1:ny,1:nx,NRHS),DUp_ksi(1:ny,1:nx,NRHS))
      if (is_irs_j) allocate(DLp_eta(1:nx,1:ny,NRHS),Dp_eta(1:nx,1:ny,NRHS),DUp_eta(1:nx,1:ny,NRHS))
      if (is_irs_k) allocate(DLp_z(1:nx,1:nz,NRHS),Dp_z(1:nx,1:nz,NRHS),DUp_z(1:nx,1:nz,NRHS))

      ! Create PaScaL_TDMA plans for several tridiagonal systems.
      ! For x direction
      if (is_irs_i) call PaScaL_TDMA_plan_many_create(p_many_x,ny,iprocyz,nprocyz,COMMYZ)
      ! For y direction
      if (is_irs_j) call PaScaL_TDMA_plan_many_create(p_many_y,nx,iprocxz,nprocxz,COMMXZ)
      ! For z direction
      if (.not.is_2d .and. is_irs_k) call PaScaL_TDMA_plan_many_create(p_many_z,nx,iprocxy,nprocxy,COMMXY)

   ! Wrong parameter for type of parallelisation
   ! ===========================================
   else
      call mpistop('Definition of type_para for irs incorrect. Only ghost points for IRS2-IRS4 and Pascal &
                   for IRS2 implemented. Shutting down...',0)
   end if

end subroutine init_irs_c_v0

!===============================================================================
subroutine init_grid_irs_c_v0
!===============================================================================
  !> Initialisation of metrics arrays for IRS, based on ijacob and idz
!===============================================================================
   use mod_bc
   use mod_grid
   use mod_constant
   use mod_time
   use mod_mpi
   implicit none
   integer :: i,j

   ! Allocation for idz_irs and ijacob_irs
   ! -------------------------------------
   allocate(ijacob_irs(1-ngh:nx+ngh,1-ngh:ny+ngh))
   if (.not.is_2d) allocate(idz_irs(1-ngh:nz+ngh))


   ! Jacobian of metrics transformation
   ! ----------------------------------
   ijacob_irs=1.0_wp
   do i=ndx-ngh,nfx+ngh
      do j=1,ny
         ijacob_irs(i,j) = x_ksi(i,j)*y_eta(i,j) - x_eta(i,j)*y_ksi(i,j)
      end do
   end do

   do i=1,nx
      do j=ndy-ngh,ndy-1
         ijacob_irs(i,j) = x_ksi(i,j)*y_eta(i,j) - x_eta(i,j)*y_ksi(i,j)
      end do
   end do

   do i=1,nx
      do j=nfy+1,nfy+ngh
         ijacob_irs(i,j) = x_ksi(i,j)*y_eta(i,j) - x_eta(i,j)*y_ksi(i,j)
      end do
   end do

   ! Inverse of Jacobian
   ! -------------------
   ijacob_irs = 1.0_wp/ijacob_irs


   ! Definition of idz_irs
   ! ---------------------
   if (.not.is_2d) then
      ! Attribution for interior points
      ! -------------------------------
      idz_irs(1:nz) = idz(1:nz)

      ! Communication
      ! -------------
      call grid_comm_irs_c_v0(idz_irs,nz,3)

      if (iirs==4)then
         idz_irs(ndz+1-ngh:nfz+ngh) = 0.5_wp * (idz_irs(ndz-ngh:nfz-1+ngh) + idz_irs(ndz+1-ngh:nfz+ngh))
      end if
   end if

end subroutine init_grid_irs_c_v0

!===============================================================================
subroutine grid_comm_irs_c_v0(dx_,nx_,dir)
!===============================================================================
  !> Communicate Cartesian metrics to fill ghost pts (for viscous metrics only)
!===============================================================================
  use mod_mpi_part
  use mod_grid
  use mod_constant
  implicit none
  ! ----------------------------------------------------------------------------
  integer, intent(in) :: dir  ! direction
  integer, intent(in) :: nx_  ! size
  real(wp), dimension(1-ngh:nx_+ngh), intent(inout) :: dx_ ! metrics
  ! ---------------------------------------------------------------------------
  integer :: n1,n2   ! neighbors
  integer :: type_dx ! MPI type
  ! ----------------------------------------------------------------------------

  ! Define neighbors depending on dimension direction
  ! -------------------------------------------------
  n1 = 2*dir
  n2 = 2*dir-1

  ! MPI type for 1D exchange of size ngh_
  ! -------------------------------------
  call MPI_TYPE_VECTOR(ngh,1,1,MPI_DOUBLE_PRECISION,type_dx,info)
  call MPI_TYPE_COMMIT(type_dx,info)

  ! MPI SENDRECV
  ! ------------
  ! Send to neighbor 1 and reception from neighbor 2
  call MPI_SENDRECV(dx_(nx_-ngh+1),1,type_dx,neighbor(n1),tag &
                   ,dx_(   -ngh+1),1,type_dx,neighbor(n2),tag,COMM_global,status,info)
  ! Send to neighbor 2 and reception from neighbor 1
  call MPI_SENDRECV(dx_(1)    ,1,type_dx,neighbor(n2),tag &
                   ,dx_(nx_+1),1,type_dx,neighbor(n1),tag,COMM_global,status,info)

end subroutine grid_comm_irs_c_v0

!!$!===============================================================================
!!$subroutine irs_routine_c_v0
!!$!===============================================================================
!!$  !> Application of the Implicit Residual Smoothing (IRS)
!!$!===============================================================================
!!$   use mod_interface
!!$   implicit none
!!$
!!$   call update_var_in_dw
!!$
!!$   call irs_solver
!!$
!!$   call update_var_of_dw
!!$
!!$end subroutine irs_routine_c_v0

!===============================================================================
subroutine irs2_ngh_c_v0
!===============================================================================
  !> Application of the ghost points parallelisation for IRS2 in curvilinear
!===============================================================================
   use mod_flow
   use mod_constant
   use mod_time
   use mod_mpi
   use mod_interface
   use warnstop
   implicit none
   integer          :: ind1,ind2,ind3
   integer :: ind1min,ind1max,ind2min,ind2max,ind3min,ind3max
   real(wp),dimension(ndx-ngh:nfx+ngh,ndy-ngh:nfy+ngh,1:nz)   :: vc_ksi, vc_eta
   real(wp),dimension(ndx-ngh:nfx+ngh,ndy-ngh:nfy+ngh)   :: n2_ksi, n2_eta
   real(wp),dimension(1:ny,ndx-ngh:nfx+ngh)   :: array_dksi
   real(wp),dimension(1:ny,ndx-ngh:nfx+ngh-1) :: array_aksi,array_cksi
   real(wp),dimension(1:ny,ndx-ngh:nfx+ngh)   :: RHS_ksi,RHSu_ksi,RHSv_ksi,RHSw_ksi,RHSe_ksi
   real(wp),dimension(1:ny)                           :: temp_ksi
   real(wp),dimension(1:nx,ndy-ngh:nfy+ngh)   :: array_deta
   real(wp),dimension(1:nx,ndy-ngh:nfy+ngh-1) :: array_aeta,array_ceta
   real(wp),dimension(1:nx)                           :: temp_eta
   real(wp),dimension(1:nx,ndz-ngh:nfz+ngh)   :: array_dz
   real(wp),dimension(1:nx,ndz-ngh:nfz+ngh-1) :: array_az,array_cz
   real(wp),dimension(1:nx,ndz-ngh:nfz+ngh)   :: RHS_z,RHSu_z,RHSv_z,RHSw_z,RHSe_z
   real(wp),dimension(1:nx)                           :: temp_z
   real(wp),dimension(ndx-ngh:nfx+ngh,ndy-ngh:nfy+ngh,ndz-ngh:nfz+ngh)  :: r_spec2
   real(wp) :: coef,rspecmh,rspecph

   !----------------------------------------
   ! Implicitation of direction ksi
   !----------------------------------------
   if (is_irs_i) then
!      ind3 : i;     ind1 : k;     ind2 : j
      ind1min = 1;   ind2min = 1;   ind3min = ndx - ngh
      ind1max = nz;  ind2max = ny;  ind3max = nfx + ngh

      ! If wall, do not take into account the wall point
      ! ================================================
      if (BC_face(2,1)%sort==0) ind2min = 2
      if (BC_face(2,2)%sort==0) ind2max = ny - 1
      if (BC_face(3,1)%sort==0) ind1min = 2
      if (BC_face(3,2)%sort==0) ind1max = nz - 1

      ! Communication of ghost points
      ! =============================
      ! call communication_EW(Krho,Krhou,Krhov,Krhow,Krhoe)
      call communication_(Krho,Krhou,Krhov,Krhow,Krhoe)

      ! Calculation of contra-variant speed and norm of ksi
      ! ===================================================
      do ind2=ind2min,ind2max
         do ind3=ind3min,ind3max
            n2_ksi(ind3,ind2) = (y_eta(ind3,ind2)**2 + x_eta(ind3,ind2)**2) * ijacob_irs(ind3,ind2)**2
         end do
      end do

      do ind1=ind1min,ind1max
         do ind2=ind2min,ind2max
            do ind3=ind3min,ind3max
               vc_ksi(ind3,ind2,ind1) = (uu(ind3,ind2,ind1)*y_eta(ind3,ind2)-vv(ind3,ind2,ind1)*x_eta(ind3,ind2)) * ijacob_irs(ind3,ind2)
            end do
         end do
      end do

      ! Calculation of the spectral radius at ind-1/2
      ! =============================================
      do ind1=ind1min,ind1max
         do ind2=ind2min,ind2max
            do ind3=ind3min+1,ind3max
               rspecmh = sqrt(vc_ksi(ind3-1,ind2,ind1)**2) + c_(ind3-1,ind2,ind1)*sqrt(n2_ksi(ind3-1,ind2))
               rspecph = sqrt(vc_ksi(ind3,ind2,ind1)**2)   + c_(ind3,ind2,ind1)*sqrt(n2_ksi(ind3,ind2))
               r_spec2(ind3,ind2,ind1) = (0.5 * (rspecmh + rspecph))**2
            enddo
         enddo
      enddo

      ! ! ===========
      ! !    Test
      ! ! ===========
      ! ! A ENLEVER
      ! ! Calculation of contra-variant speed and norm of eta
      ! ! ===================================================
      ! do ind2=ind2min,ind2max
      !    do ind3=ind3min,ind3max
      !       n2_eta(ind3,ind2) = (y_ksi(ind3,ind2)**2 + x_ksi(ind3,ind2)**2) * ijacob_irs(ind3,ind2)**2
      !    end do
      ! end do

      ! do ind1=ind1min,ind1max
      !    do ind2=ind2min,ind2max
      !       do ind3=ind3min,ind3max
      !          vc_eta(ind3,ind2,ind1) = (vv(ind3,ind2,ind1)*x_ksi(ind3,ind2) - uu(ind3,ind2,ind1)*y_ksi(ind3,ind2)) * ijacob_irs(ind3,ind2)
      !       end do
      !    end do
      ! end do

      ! ! Calculation of the spectral radius at ind-1/2
      ! ! =============================================
      ! do ind1=ind1min,ind1max
      !    do ind2=ind2min,ind2max
      !       do ind3=ind3min+1,ind3max
      !          ! rspecmh = sqrt(vc_ksi(ind3-1,ind2,ind1)**2) + sqrt(vc_eta(ind3-1,ind2,ind1)**2) + c_(ind3-1,ind2,ind1)*(sqrt(n2_ksi(ind3-1,ind2)) + sqrt(n2_eta(ind3-1,ind2)))
      !          ! rspecph = sqrt(vc_ksi(ind3,ind2,ind1)**2)   + sqrt(vc_eta(ind3,ind2,ind1)**2)   + c_(ind3,ind2,ind1)*(  sqrt(n2_ksi(ind3,ind2))   + sqrt(n2_eta(ind3,ind2)))
      !          rspecmh = sqrt(vc_ksi(ind3-1,ind2,ind1)**2 + vc_eta(ind3-1,ind2,ind1)**2) + c_(ind3-1,ind2,ind1)*(sqrt(n2_ksi(ind3-1,ind2)) + sqrt(n2_eta(ind3-1,ind2)))
      !          rspecph = sqrt(vc_ksi(ind3,ind2,ind1)**2   + vc_eta(ind3,ind2,ind1)**2)   + c_(ind3,ind2,ind1)*(  sqrt(n2_ksi(ind3,ind2))   + sqrt(n2_eta(ind3,ind2)))
      !          ! rspecmh = sqrt(vc_ksi(ind3-1,ind2,ind1)**2) + sqrt(vc_eta(ind3-1,ind2,ind1)**2) + c_(ind3-1,ind2,ind1)*ijacob_irs(ind3-1,ind2)
      !          ! rspecph = sqrt(vc_ksi(ind3,ind2,ind1)**2)   + sqrt(vc_eta(ind3,ind2,ind1)**2)   + c_(ind3,ind2,ind1)*ijacob_irs(ind3,ind2)
      !          ! rspecmh = sqrt(vc_ksi(ind3-1,ind2,ind1)**2 + vc_eta(ind3-1,ind2,ind1)**2) + c_(ind3-1,ind2,ind1)*(sqrt(n2_ksi(ind3-1,ind2) + n2_eta(ind3-1,ind2)))
      !          ! rspecph = sqrt(vc_ksi(ind3,ind2,ind1)**2   + vc_eta(ind3,ind2,ind1)**2)   + c_(ind3,ind2,ind1)*(  sqrt(n2_ksi(ind3,ind2)   + n2_eta(ind3,ind2)))
      !          r_spec2(ind3,ind2,ind1) = (0.5 * (rspecmh + rspecph))**2
      !       enddo
      !    enddo
      ! enddo

      ! ! =========
      ! ! fin test
      ! ! =========

      do ind1 = ind1min, ind1max
         ! Construction of the tridiagonal matrix
         ! The linear systems is Mx = y with :
         !     _M = (m_ij) for i,j in {1,..,n} and with :
         !                 m_ii = d_i for i in {1,..,n}
         !                 m_ii-1 = c_i-1 for i in {2,..,n}
         !                 m_ii+1 = a_i for i in {1,..,n-1}
         !     _ y = RHS = dwi
         !     _d and y are of dimensions n, a and c of dimensions n-1,
         ind3=ind3min
         do ind2 = ind2min, ind2max
            array_aksi(ind2,ind3)  = - theta_irs1 * deltat * (r_spec2(ind3+1,ind2,ind1))**0.5
            array_dksi(ind2,ind3)  = 1. - array_aksi(ind2,ind3)
         end do

         ind3=ind3max
         do ind2 = ind2min, ind2max
            array_cksi(ind2,ind3-1)  = - theta_irs1 * deltat  * (r_spec2(ind3,ind2,ind1))**0.5
            array_dksi(ind2,ind3)    = 1. - array_cksi(ind2,ind3-1)
         end do

         coef = - theta_irs2 * deltat**2
         do ind3 = ind3min+1,ind3max-1
            do ind2 = ind2min, ind2max
               array_cksi(ind2,ind3-1)= coef * r_spec2(ind3,ind2,ind1)
               array_aksi(ind2,ind3)  = coef * r_spec2(ind3+1,ind2,ind1)
               array_dksi(ind2,ind3)  = 1 - (array_cksi(ind2,ind3-1) + array_aksi(ind2,ind3))
            enddo
         enddo

         ! Filling of the RHS
         do ind3 = ind3min,ind3max
            do ind2 = ind2min,ind2max
                RHS_ksi(ind2,ind3) =  Krho(ind3,ind2,ind1)
               RHSu_ksi(ind2,ind3) = Krhou(ind3,ind2,ind1)
               RHSv_ksi(ind2,ind3) = Krhov(ind3,ind2,ind1)
               RHSw_ksi(ind2,ind3) = Krhow(ind3,ind2,ind1)
               RHSe_ksi(ind2,ind3) = Krhoe(ind3,ind2,ind1)
            end do
         end do

         ! Resolution of the tridiagonal linear system
         ! ===========================================
         ! Thomas' algorithm - Forward elimination phase
         do ind2 = ind2min, ind2max
            temp_ksi(ind2) = 1. / array_dksi(ind2,ind3min)
            array_aksi(ind2,ind3min) = array_aksi(ind2,ind3min) * temp_ksi(ind2)
             RHS_ksi(ind2,ind3min) =  RHS_ksi(ind2,ind3min) * temp_ksi(ind2)
            RHSu_ksi(ind2,ind3min) = RHSu_ksi(ind2,ind3min) * temp_ksi(ind2)
            RHSv_ksi(ind2,ind3min) = RHSv_ksi(ind2,ind3min) * temp_ksi(ind2)
            RHSw_ksi(ind2,ind3min) = RHSw_ksi(ind2,ind3min) * temp_ksi(ind2)
            RHSe_ksi(ind2,ind3min) = RHSe_ksi(ind2,ind3min) * temp_ksi(ind2)
         end do

         do ind3 = ind3min+1, ind3max-1
            do ind2 = ind2min, ind2max
               temp_ksi(ind2) = 1. / (array_dksi(ind2,ind3) - array_cksi(ind2,ind3-1)*array_aksi(ind2,ind3-1))
               array_aksi(ind2,ind3)   = array_aksi(ind2,ind3) * temp_ksi(ind2)
                RHS_ksi(ind2,ind3) = ( RHS_ksi(ind2,ind3) - array_cksi(ind2,ind3-1)* RHS_ksi(ind2,ind3-1)) * temp_ksi(ind2)
               RHSu_ksi(ind2,ind3) = (RHSu_ksi(ind2,ind3) - array_cksi(ind2,ind3-1)*RHSu_ksi(ind2,ind3-1)) * temp_ksi(ind2)
               RHSv_ksi(ind2,ind3) = (RHSv_ksi(ind2,ind3) - array_cksi(ind2,ind3-1)*RHSv_ksi(ind2,ind3-1)) * temp_ksi(ind2)
               RHSw_ksi(ind2,ind3) = (RHSw_ksi(ind2,ind3) - array_cksi(ind2,ind3-1)*RHSw_ksi(ind2,ind3-1)) * temp_ksi(ind2)
               RHSe_ksi(ind2,ind3) = (RHSe_ksi(ind2,ind3) - array_cksi(ind2,ind3-1)*RHSe_ksi(ind2,ind3-1)) * temp_ksi(ind2)
            end do
         end do

         ! Backward substitution phase
         do ind2 = ind2min, ind2max
            temp_ksi(ind2) = 1. / (array_dksi(ind2,ind3max) - array_cksi(ind2,ind3max-1)*array_aksi(ind2,ind3max-1))
             RHS_ksi(ind2,ind3max) = ( RHS_ksi(ind2,ind3max) - array_cksi(ind2,ind3max-1)* RHS_ksi(ind2,ind3max-1))*temp_ksi(ind2)
            RHSu_ksi(ind2,ind3max) = (RHSu_ksi(ind2,ind3max) - array_cksi(ind2,ind3max-1)*RHSu_ksi(ind2,ind3max-1))*temp_ksi(ind2)
            RHSv_ksi(ind2,ind3max) = (RHSv_ksi(ind2,ind3max) - array_cksi(ind2,ind3max-1)*RHSv_ksi(ind2,ind3max-1))*temp_ksi(ind2)
            RHSw_ksi(ind2,ind3max) = (RHSw_ksi(ind2,ind3max) - array_cksi(ind2,ind3max-1)*RHSw_ksi(ind2,ind3max-1))*temp_ksi(ind2)
            RHSe_ksi(ind2,ind3max) = (RHSe_ksi(ind2,ind3max) - array_cksi(ind2,ind3max-1)*RHSe_ksi(ind2,ind3max-1))*temp_ksi(ind2)
         end do

         do ind3 = ind3max-1, ind3min, -1
            do ind2 = ind2min, ind2max
                RHS_ksi(ind2,ind3) =  RHS_ksi(ind2,ind3) - array_aksi(ind2,ind3)* RHS_ksi(ind2,ind3+1)
               RHSu_ksi(ind2,ind3) = RHSu_ksi(ind2,ind3) - array_aksi(ind2,ind3)*RHSu_ksi(ind2,ind3+1)
               RHSv_ksi(ind2,ind3) = RHSv_ksi(ind2,ind3) - array_aksi(ind2,ind3)*RHSv_ksi(ind2,ind3+1)
               RHSw_ksi(ind2,ind3) = RHSw_ksi(ind2,ind3) - array_aksi(ind2,ind3)*RHSw_ksi(ind2,ind3+1)
               RHSe_ksi(ind2,ind3) = RHSe_ksi(ind2,ind3) - array_aksi(ind2,ind3)*RHSe_ksi(ind2,ind3+1)
            end do
         enddo

         ! Update with the solution
         do ind3 = ind3min,ind3max
            do ind2 = ind2min,ind2max
                Krho(ind3,ind2,ind1) =  RHS_ksi(ind2,ind3)
               Krhou(ind3,ind2,ind1) = RHSu_ksi(ind2,ind3)
               Krhov(ind3,ind2,ind1) = RHSv_ksi(ind2,ind3)
               Krhow(ind3,ind2,ind1) = RHSw_ksi(ind2,ind3)
               Krhoe(ind3,ind2,ind1) = RHSe_ksi(ind2,ind3)
            end do
         end do
      enddo

   end if

   !----------------------------------------
   ! Implicitation of direction eta
   !----------------------------------------
   if (is_irs_j) then
!      ind3 : j;     ind1 : k;     ind2 : i
      ind1min = 1;   ind2min = 1;   ind3min = ndy - ngh
      ind1max = nz;  ind2max = nx;  ind3max = nfy + ngh

      ! If wall, do not take into account the wall point
      ! ================================================
      if (BC_face(1,1)%sort==0) ind2min = 2
      if (BC_face(1,2)%sort==0) ind2max = nx - 1
      if (BC_face(3,1)%sort==0) ind1min = 2
      if (BC_face(3,2)%sort==0) ind1max = nz - 1

      ! Communication of ghost points
      ! =============================
      ! call communication_NS(Krho,Krhou,Krhov,Krhow,Krhoe)
      call communication_(Krho,Krhou,Krhov,Krhow,Krhoe)

      ! Calculation of contra-variant speed and norm of eta
      ! ===================================================
      do ind3=ind3min,ind3max
         do ind2=ind2min,ind2max
            n2_eta(ind2,ind3) = (y_ksi(ind2,ind3)**2 + x_ksi(ind2,ind3)**2) * ijacob_irs(ind2,ind3)**2
         end do
      end do

      do ind1=ind1min,ind1max
         do ind3=ind3min,ind3max
            do ind2=ind2min,ind2max
               vc_eta(ind2,ind3,ind1) = (vv(ind2,ind3,ind1)*x_ksi(ind2,ind3) - uu(ind2,ind3,ind1)*y_ksi(ind2,ind3)) * ijacob_irs(ind2,ind3)
            end do
         end do
      end do

      ! Calculation of the spectral radius at ind-1/2
      ! =============================================
      do ind1=ind1min,ind1max
         do ind3=ind3min+1,ind3max
            do ind2=ind2min,ind2max
               rspecmh = sqrt(vc_eta(ind2,ind3-1,ind1)**2) + c_(ind2,ind3-1,ind1)*sqrt(n2_eta(ind2,ind3-1))
               rspecph = sqrt(vc_eta(ind2,ind3,ind1)**2) + c_(ind2,ind3,ind1)*sqrt(n2_eta(ind2,ind3))
               r_spec2(ind2,ind3,ind1) = (0.5 * (rspecmh + rspecph))**2
            enddo
         enddo
      enddo

      ! ! ===========
      ! !    Test
      ! ! ===========
      ! ! A ENLEVER
      ! ! Calculation of contra-variant speed and norm of ksi
      ! ! ===================================================
      ! do ind3=ind3min,ind3max
      !    do ind2=ind2min,ind2max
      !       n2_ksi(ind2,ind3) = (y_eta(ind2,ind3)**2 + x_eta(ind2,ind3)**2) * ijacob_irs(ind2,ind3)**2
      !    end do
      ! end do

      ! do ind1=ind1min,ind1max
      !    do ind2=ind2min,ind2max
      !       do ind3=ind3min,ind3max
      !          vc_ksi(ind2,ind3,ind1) = (uu(ind2,ind3,ind1)*y_eta(ind2,ind3)-vv(ind2,ind3,ind1)*x_eta(ind2,ind3)) * ijacob_irs(ind2,ind3)
      !       end do
      !    end do
      ! end do

      ! ! Calculation of the spectral radius at ind-1/2
      ! ! =============================================
      ! do ind1=ind1min,ind1max
      !    do ind3=ind3min+1,ind3max
      !       do ind2=ind2min,ind2max
      !          ! rspecmh = sqrt(vc_eta(ind2,ind3-1,ind1)**2) + sqrt(vc_ksi(ind2,ind3-1,ind1)**2) + c_(ind2,ind3-1,ind1)*(sqrt(n2_eta(ind2,ind3-1)) + sqrt(n2_ksi(ind2,ind3-1)))
      !          ! rspecph = sqrt(vc_eta(ind2,ind3,ind1)**2)   + sqrt(vc_ksi(ind2,ind3,ind1)**2)   + c_(ind2,ind3,ind1)*(  sqrt(n2_eta(ind2,ind3)) + sqrt(n2_ksi(ind2,ind3)))
      !          rspecmh = sqrt(vc_eta(ind2,ind3-1,ind1)**2 + vc_ksi(ind2,ind3-1,ind1)**2) + c_(ind2,ind3-1,ind1)*(sqrt(n2_eta(ind2,ind3-1)) + sqrt(n2_ksi(ind2,ind3-1)))
      !          rspecph = sqrt(vc_eta(ind2,ind3,ind1)**2   + vc_ksi(ind2,ind3,ind1)**2)   + c_(ind2,ind3,ind1)*(  sqrt(n2_eta(ind2,ind3))   + sqrt(n2_ksi(ind2,ind3)))
      !          ! rspecmh = sqrt(vc_eta(ind2,ind3-1,ind1)**2) + sqrt(vc_ksi(ind2,ind3-1,ind1)**2) + c_(ind2,ind3-1,ind1)*ijacob_irs(ind2,ind3-1)
      !          ! rspecph = sqrt(vc_eta(ind2,ind3,ind1)**2)   + sqrt(vc_ksi(ind2,ind3,ind1)**2)   + c_(ind2,ind3,ind1)*ijacob_irs(ind2,ind3)
      !          ! rspecmh = sqrt(vc_eta(ind2,ind3-1,ind1)**2 + vc_ksi(ind2,ind3-1,ind1)**2) + c_(ind2,ind3-1,ind1)*(sqrt(n2_eta(ind2,ind3-1) + n2_ksi(ind2,ind3-1)))
      !          ! rspecph = sqrt(vc_eta(ind2,ind3,ind1)**2   + vc_ksi(ind2,ind3,ind1)**2)   + c_(ind2,ind3,ind1)*(  sqrt(n2_eta(ind2,ind3)   + n2_ksi(ind2,ind3)))
      !          r_spec2(ind2,ind3,ind1) = (0.5 * (rspecmh + rspecph))**2
      !       enddo
      !    enddo
      ! enddo

      ! ! =========
      ! ! fin test
      ! ! =========

      do ind1 = ind1min, ind1max
         ! Construction of the tridiagonal matrix
         ind3=ind3min
         do ind2 = ind2min, ind2max
            array_aeta(ind2,ind3)  = - theta_irs1 * deltat * (r_spec2(ind2,ind3+1,ind1))**0.5
            array_deta(ind2,ind3)  = 1. - array_aeta(ind2,ind3)
         end do

         ind3=ind3max
         do ind2 = ind2min, ind2max
            array_ceta(ind2,ind3-1)  = - theta_irs1 * deltat  * (r_spec2(ind2,ind3,ind1))**0.5
            array_deta(ind2,ind3)    = 1. - array_ceta(ind2,ind3-1)
         end do

         coef = - theta_irs2 * deltat**2
         do ind3 = ind3min+1,ind3max-1
            do ind2 = ind2min, ind2max
               array_ceta(ind2,ind3-1)= coef * r_spec2(ind2,ind3,ind1)
               array_aeta(ind2,ind3)  = coef * r_spec2(ind2,ind3+1,ind1)
               array_deta(ind2,ind3)  = 1 - (array_ceta(ind2,ind3-1) + array_aeta(ind2,ind3))
            enddo
         enddo

         ! Resolution of the tridiagonal linear system
         ! ===========================================
         ! Thomas' algorithm - Forward elimination phase
         do ind2 = ind2min, ind2max
            temp_eta(ind2) = 1. / array_deta(ind2,ind3min)
            array_aeta(ind2,ind3min) = array_aeta(ind2,ind3min) * temp_eta(ind2)
             Krho(ind2,ind3min,ind1) =  Krho(ind2,ind3min,ind1) * temp_eta(ind2)
            Krhou(ind2,ind3min,ind1) = Krhou(ind2,ind3min,ind1) * temp_eta(ind2)
            Krhov(ind2,ind3min,ind1) = Krhov(ind2,ind3min,ind1) * temp_eta(ind2)
            Krhow(ind2,ind3min,ind1) = Krhow(ind2,ind3min,ind1) * temp_eta(ind2)
            Krhoe(ind2,ind3min,ind1) = Krhoe(ind2,ind3min,ind1) * temp_eta(ind2)
         end do

         do ind3 = ind3min+1, ind3max-1
            do ind2 = ind2min, ind2max
               temp_eta(ind2) = 1. / (array_deta(ind2,ind3) - array_ceta(ind2,ind3-1)*array_aeta(ind2,ind3-1))
               array_aeta(ind2,ind3)   = array_aeta(ind2,ind3) * temp_eta(ind2)
                Krho(ind2,ind3,ind1) = ( Krho(ind2,ind3,ind1) - array_ceta(ind2,ind3-1)* Krho(ind2,ind3-1,ind1)) * temp_eta(ind2)
               Krhou(ind2,ind3,ind1) = (Krhou(ind2,ind3,ind1) - array_ceta(ind2,ind3-1)*Krhou(ind2,ind3-1,ind1)) * temp_eta(ind2)
               Krhov(ind2,ind3,ind1) = (Krhov(ind2,ind3,ind1) - array_ceta(ind2,ind3-1)*Krhov(ind2,ind3-1,ind1)) * temp_eta(ind2)
               Krhow(ind2,ind3,ind1) = (Krhow(ind2,ind3,ind1) - array_ceta(ind2,ind3-1)*Krhow(ind2,ind3-1,ind1)) * temp_eta(ind2)
               Krhoe(ind2,ind3,ind1) = (Krhoe(ind2,ind3,ind1) - array_ceta(ind2,ind3-1)*Krhoe(ind2,ind3-1,ind1)) * temp_eta(ind2)
            end do
         end do

         ! Backward substitution phase
         do ind2 = ind2min, ind2max
            temp_eta(ind2) = 1. / (array_deta(ind2,ind3max) - array_ceta(ind2,ind3max-1)*array_aeta(ind2,ind3max-1))
             Krho(ind2,ind3max,ind1) = ( Krho(ind2,ind3max,ind1) - array_ceta(ind2,ind3max-1)* Krho(ind2,ind3max-1,ind1))*temp_eta(ind2)
            Krhou(ind2,ind3max,ind1) = (Krhou(ind2,ind3max,ind1) - array_ceta(ind2,ind3max-1)*Krhou(ind2,ind3max-1,ind1))*temp_eta(ind2)
            Krhov(ind2,ind3max,ind1) = (Krhov(ind2,ind3max,ind1) - array_ceta(ind2,ind3max-1)*Krhov(ind2,ind3max-1,ind1))*temp_eta(ind2)
            Krhow(ind2,ind3max,ind1) = (Krhow(ind2,ind3max,ind1) - array_ceta(ind2,ind3max-1)*Krhow(ind2,ind3max-1,ind1))*temp_eta(ind2)
            Krhoe(ind2,ind3max,ind1) = (Krhoe(ind2,ind3max,ind1) - array_ceta(ind2,ind3max-1)*Krhoe(ind2,ind3max-1,ind1))*temp_eta(ind2)
         end do

         do ind3 = ind3max-1, ind3min, -1
            do ind2 = ind2min, ind2max
                Krho(ind2,ind3,ind1) =  Krho(ind2,ind3,ind1) - array_aeta(ind2,ind3)* Krho(ind2,ind3+1,ind1)
               Krhou(ind2,ind3,ind1) = Krhou(ind2,ind3,ind1) - array_aeta(ind2,ind3)*Krhou(ind2,ind3+1,ind1)
               Krhov(ind2,ind3,ind1) = Krhov(ind2,ind3,ind1) - array_aeta(ind2,ind3)*Krhov(ind2,ind3+1,ind1)
               Krhow(ind2,ind3,ind1) = Krhow(ind2,ind3,ind1) - array_aeta(ind2,ind3)*Krhow(ind2,ind3+1,ind1)
               Krhoe(ind2,ind3,ind1) = Krhoe(ind2,ind3,ind1) - array_aeta(ind2,ind3)*Krhoe(ind2,ind3+1,ind1)
            end do
         enddo
      enddo

   end if

   !----------------------------------------
   ! Implicitation of direction z
   !----------------------------------------
   if (is_irs_k .and. .not. is_2d) then
      !      ind3 : k;     ind1 : j;     ind2 : i
      ind1min = 1;   ind2min = 1;   ind3min = ndz - ngh
      ind1max = ny;  ind2max = nx;  ind3max = nfz + ngh

      ! If wall, do not take into account the wall point
      ! ================================================
      if (BC_face(1,1)%sort==0) ind2min = 2
      if (BC_face(1,2)%sort==0) ind2max = nx - 1
      if (BC_face(2,1)%sort==0) ind1min = 2
      if (BC_face(2,2)%sort==0) ind1max = ny - 1

      ! Communication of ghost points
      ! =============================
      ! call communication_FB(Krho,Krhou,Krhov,Krhow,Krhoe)
      call communication_(Krho,Krhou,Krhov,Krhow,Krhoe)

      ! Calculation of the spectral radius at ind-1/2
      ! =============================================
      do ind3=ind3min+1,ind3max
         do ind1=ind1min,ind1max
            do ind2=ind2min,ind2max
               rspecmh = sqrt(ww(ind2,ind1,ind3-1)**2) + c_(ind2,ind1,ind3-1)
               rspecph = sqrt(ww(ind2,ind1,ind3)**2)   + c_(ind2,ind1,ind3)
               r_spec2(ind2,ind1,ind3) = (0.5 * (rspecmh + rspecph))**2
            enddo
         enddo
      enddo

      do ind1 = ind1min, ind1max
         ind3=ind3min
         do ind2 = ind2min, ind2max
            array_az(ind2,ind3)  = - theta_irs1 * deltat * idz_irs(ind3)* (r_spec2(ind2,ind1,ind3+1))**0.5
            array_dz(ind2,ind3)  = 1. - array_az(ind2,ind3)
         end do

         ind3=ind3max
         do ind2 = ind2min, ind2max
            array_cz(ind2,ind3-1)  = - theta_irs1 * deltat * idz_irs(ind3) * (r_spec2(ind2,ind1,ind3))**0.5
            array_dz(ind2,ind3)    = 1. - array_cz(ind2,ind3-1)
         end do

         do ind3 = ind3min+1,ind3max-1
            do ind2 = ind2min, ind2max
               temp_z(ind2)         = - theta_irs2 * (deltat * idz_irs(ind3))**2
               array_cz(ind2,ind3-1)= temp_z(ind2) * r_spec2(ind2,ind1,ind3)
               array_az(ind2,ind3)  = temp_z(ind2) * r_spec2(ind2,ind1,ind3+1)
               array_dz(ind2,ind3)  = 1 - (array_cz(ind2,ind3-1) + array_az(ind2,ind3))
            enddo
         enddo

         ! Filling of the RHS
         do ind3 = ind3min,ind3max
            do ind2 = ind2min,ind2max
                RHS_z(ind2,ind3) =  Krho(ind2,ind1,ind3)
               RHSu_z(ind2,ind3) = Krhou(ind2,ind1,ind3)
               RHSv_z(ind2,ind3) = Krhov(ind2,ind1,ind3)
               RHSw_z(ind2,ind3) = Krhow(ind2,ind1,ind3)
               RHSe_z(ind2,ind3) = Krhoe(ind2,ind1,ind3)
            end do
         end do

         ! Resolution of the tridiagonal linear system
         ! ===========================================
         ! Thomas' algorithm - Forward elimination phase
         do ind2 = ind2min, ind2max
            temp_z(ind2) = 1. / array_dz(ind2,ind3min)
            array_az(ind2,ind3min) = array_az(ind2,ind3min) * temp_z(ind2)
             RHS_z(ind2,ind3min) =  RHS_z(ind2,ind3min) * temp_z(ind2)
            RHSu_z(ind2,ind3min) = RHSu_z(ind2,ind3min) * temp_z(ind2)
            RHSv_z(ind2,ind3min) = RHSv_z(ind2,ind3min) * temp_z(ind2)
            RHSw_z(ind2,ind3min) = RHSw_z(ind2,ind3min) * temp_z(ind2)
            RHSe_z(ind2,ind3min) = RHSe_z(ind2,ind3min) * temp_z(ind2)
         end do

         do ind3 = ind3min+1, ind3max-1
            do ind2 = ind2min, ind2max
               temp_z(ind2) = 1. / (array_dz(ind2,ind3) - array_cz(ind2,ind3-1)*array_az(ind2,ind3-1))
               array_az(ind2,ind3)   = array_az(ind2,ind3) * temp_z(ind2)
                RHS_z(ind2,ind3) = ( RHS_z(ind2,ind3) - array_cz(ind2,ind3-1)* RHS_z(ind2,ind3-1)) * temp_z(ind2)
               RHSu_z(ind2,ind3) = (RHSu_z(ind2,ind3) - array_cz(ind2,ind3-1)*RHSu_z(ind2,ind3-1)) * temp_z(ind2)
               RHSv_z(ind2,ind3) = (RHSv_z(ind2,ind3) - array_cz(ind2,ind3-1)*RHSv_z(ind2,ind3-1)) * temp_z(ind2)
               RHSw_z(ind2,ind3) = (RHSw_z(ind2,ind3) - array_cz(ind2,ind3-1)*RHSw_z(ind2,ind3-1)) * temp_z(ind2)
               RHSe_z(ind2,ind3) = (RHSe_z(ind2,ind3) - array_cz(ind2,ind3-1)*RHSe_z(ind2,ind3-1)) * temp_z(ind2)
            end do
         end do

         ! Backward substitution phase
         do ind2 = ind2min, ind2max
            temp_z(ind2) = 1. / (array_dz(ind2,ind3max) - array_cz(ind2,ind3max-1)*array_az(ind2,ind3max-1))
             RHS_z(ind2,ind3max) = ( RHS_z(ind2,ind3max) - array_cz(ind2,ind3max-1)* RHS_z(ind2,ind3max-1))*temp_z(ind2)
            RHSu_z(ind2,ind3max) = (RHSu_z(ind2,ind3max) - array_cz(ind2,ind3max-1)*RHSu_z(ind2,ind3max-1))*temp_z(ind2)
            RHSv_z(ind2,ind3max) = (RHSv_z(ind2,ind3max) - array_cz(ind2,ind3max-1)*RHSv_z(ind2,ind3max-1))*temp_z(ind2)
            RHSw_z(ind2,ind3max) = (RHSw_z(ind2,ind3max) - array_cz(ind2,ind3max-1)*RHSw_z(ind2,ind3max-1))*temp_z(ind2)
            RHSe_z(ind2,ind3max) = (RHSe_z(ind2,ind3max) - array_cz(ind2,ind3max-1)*RHSe_z(ind2,ind3max-1))*temp_z(ind2)
         end do

         do ind3 = ind3max-1, ind3min, -1
            do ind2 = ind2min, ind2max
                RHS_z(ind2,ind3) =  RHS_z(ind2,ind3) - array_az(ind2,ind3)* RHS_z(ind2,ind3+1)
               RHSu_z(ind2,ind3) = RHSu_z(ind2,ind3) - array_az(ind2,ind3)*RHSu_z(ind2,ind3+1)
               RHSv_z(ind2,ind3) = RHSv_z(ind2,ind3) - array_az(ind2,ind3)*RHSv_z(ind2,ind3+1)
               RHSw_z(ind2,ind3) = RHSw_z(ind2,ind3) - array_az(ind2,ind3)*RHSw_z(ind2,ind3+1)
               RHSe_z(ind2,ind3) = RHSe_z(ind2,ind3) - array_az(ind2,ind3)*RHSe_z(ind2,ind3+1)
            end do
         enddo

         ! Update with the solution
         do ind3 = ind3min,ind3max
            do ind2 = ind2min,ind2max
                Krho(ind2,ind1,ind3) =  RHS_z(ind2,ind3)
               Krhou(ind2,ind1,ind3) = RHSu_z(ind2,ind3)
               Krhov(ind2,ind1,ind3) = RHSv_z(ind2,ind3)
               Krhow(ind2,ind1,ind3) = RHSw_z(ind2,ind3)
               Krhoe(ind2,ind1,ind3) = RHSe_z(ind2,ind3)
            end do
         end do
      enddo

   end if
end subroutine irs2_ngh_c_v0

!===============================================================================
subroutine irs4_ngh_c_v0
!===============================================================================
  !> Application of the ghost points parallelisation for IRS4
!===============================================================================
   use mod_flow
   use mod_constant
   use mod_time
   use mod_mpi
   use mod_interface
   use warnstop
   implicit none
   integer          :: ind1,ind2,ind3
   integer :: ind1min,ind1max,ind2min,ind2max,ind3min,ind3max
   real(wp),dimension(ndx-ngh:nfx+ngh,ndy-ngh:nfy+ngh,1:nz)   :: vc_ksi, vc_eta
   real(wp),dimension(ndx-ngh:nfx+ngh,ndy-ngh:nfy+ngh)   :: n2_ksi, n2_eta
   real(wp),dimension(1:ny,ndx-ngh:nfx+ngh)   :: array_dksi
   real(wp),dimension(1:ny,ndx-ngh:nfx+ngh-1) :: array_aksi,array_cksi
   real(wp),dimension(1:ny,ndx-ngh:nfx+ngh-2) :: array_bksi,array_eksi
   real(wp),dimension(1:ny,ndx-ngh:nfx+ngh)   :: RHS_ksi,RHSu_ksi,RHSv_ksi,RHSw_ksi,RHSe_ksi
   real(wp),dimension(1:ny)                           :: temp_ksi1,temp_ksi2
   real(wp),dimension(1:nx,ndy-ngh:nfy+ngh)   :: array_deta
   real(wp),dimension(1:nx,ndy-ngh:nfy+ngh-1) :: array_aeta,array_ceta
   real(wp),dimension(1:nx,ndy-ngh:nfy+ngh-2) :: array_by,array_ey
   real(wp),dimension(1:nx)                           :: temp_eta1,temp_eta2
   real(wp),dimension(1:nx,ndz-ngh:nfz+ngh)   :: array_dz
   real(wp),dimension(1:nx,ndz-ngh:nfz+ngh-1) :: array_az,array_cz
   real(wp),dimension(1:nx,ndz-ngh:nfz+ngh-2) :: array_bz,array_ez
   real(wp),dimension(1:nx,ndz-ngh:nfz+ngh)   :: RHS_z,RHSu_z,RHSv_z,RHSw_z,RHSe_z
   real(wp),dimension(1:nx)                           :: temp_z1,temp_z2
   real(wp),dimension(ndx-ngh:nfx+ngh,ndy-ngh:nfy+ngh,ndz-ngh:nfz+ngh)  :: r_spec4 !,r_spec4_2
   real(wp) :: coef,coefmh,coefph,rspecmh,rspecph !,rspecmh2,rspecph2,rspecmh3,rspecph3
   ! ! TO BE CHANGED
   ! real(wp),dimension(ndy-ngh:nfy+ngh)   :: coef_j

   !----------------------------------------
   ! Implicitation of direction ksi
   !----------------------------------------
   if (is_irs_i) then
!      ind3 : i;     ind1 : k;     ind2 : j
      ind1min = 1;   ind2min = 1;   ind3min = ndx - ngh
      ind1max = nz;  ind2max = ny;  ind3max = nfx + ngh

      ! If wall, do not take into account the wall point
      ! ================================================
      if (BC_face(2,1)%sort==0) ind2min = 2
      if (BC_face(2,2)%sort==0) ind2max = ny - 1
      if (BC_face(3,1)%sort==0) ind1min = 2
      if (BC_face(3,2)%sort==0) ind1max = nz - 1

      ! Communication of ghost points
      ! =============================
      ! call communication_EW(Krho,Krhou,Krhov,Krhow,Krhoe)
      call communication_(Krho,Krhou,Krhov,Krhow,Krhoe)

      ! Calculation of contra-variant speed and norm of ksi
      ! ===================================================
      do ind2=ind2min,ind2max
         do ind3=ind3min,ind3max
            n2_ksi(ind3,ind2) = (y_eta(ind3,ind2)**2 + x_eta(ind3,ind2)**2) * ijacob_irs(ind3,ind2)**2
         end do
      end do

      do ind1=ind1min,ind1max
         do ind2=ind2min,ind2max
            do ind3=ind3min,ind3max
               vc_ksi(ind3,ind2,ind1) = (uu(ind3,ind2,ind1)*y_eta(ind3,ind2)-vv(ind3,ind2,ind1)*x_eta(ind3,ind2)) * ijacob_irs(ind3,ind2)
            end do
         end do
      end do

!!$      ! Calculation of the spectral radius at ind-1/2
!!$      ! =============================================
!!$      do ind1=ind1min,ind1max
!!$         do ind2=ind2min,ind2max
!!$            do ind3=ind3min+1,ind3max
!!$               rspecmh = sqrt(vc_ksi(ind3-1,ind2,ind1)**2) + c_(ind3-1,ind2,ind1)*sqrt(n2_ksi(ind3-1,ind2))
!!$               rspecph = sqrt(vc_ksi(ind3,ind2,ind1)**2)   + c_(ind3,ind2,ind1)*sqrt(n2_ksi(ind3,ind2))
!!$               r_spec4(ind3,ind2,ind1) = (0.5 * (rspecmh + rspecph))**4
!!$            enddo
!!$         enddo
!!$      enddo

      ! ===========
      !    Test
      ! ===========
      ! A ENLEVER
      ! Calculation of contra-variant speed and norm of eta
      ! ===================================================
      do ind2=ind2min,ind2max
         do ind3=ind3min,ind3max
            n2_eta(ind3,ind2) = (y_ksi(ind3,ind2)**2 + x_ksi(ind3,ind2)**2) * ijacob_irs(ind3,ind2)**2
         end do
      end do

      ! Calculation of the spectral radius at ind-1/2
      ! =============================================
      do ind1=ind1min,ind1max
         do ind2=ind2min,ind2max
            do ind3=ind3min,ind3max
               vc_eta(ind3,ind2,ind1) = (vv(ind3,ind2,ind1)*x_ksi(ind3,ind2) - uu(ind3,ind2,ind1)*y_ksi(ind3,ind2)) * ijacob_irs(ind3,ind2)
            end do
         end do
      end do

      do ind1=ind1min,ind1max
         do ind2=ind2min,ind2max
            do ind3=ind3min+1,ind3max
               !rspecmh = sqrt(vc_ksi(ind3-1,ind2,ind1)**2) + sqrt(vc_eta(ind3-1,ind2,ind1)**2) + c_(ind3-1,ind2,ind1)*(sqrt(n2_ksi(ind3-1,ind2)) + sqrt(n2_eta(ind3-1,ind2)))
               !rspecph = sqrt(vc_ksi(ind3,ind2,ind1)**2)   + sqrt(vc_eta(ind3,ind2,ind1)**2)   + c_(ind3,ind2,ind1)*(  sqrt(n2_ksi(ind3,ind2))   + sqrt(n2_eta(ind3,ind2)))
               !rspecmh = sqrt(vc_ksi(ind3-1,ind2,ind1)**2 + vc_eta(ind3-1,ind2,ind1)**2) + c_(ind3-1,ind2,ind1)*(sqrt(n2_ksi(ind3-1,ind2)) + sqrt(n2_eta(ind3-1,ind2)))
               !rspecph = sqrt(vc_ksi(ind3,ind2,ind1)**2   + vc_eta(ind3,ind2,ind1)**2)   + c_(ind3,ind2,ind1)*(  sqrt(n2_ksi(ind3,ind2))   + sqrt(n2_eta(ind3,ind2)))
               !rspecmh = sqrt(vc_ksi(ind3-1,ind2,ind1)**2) + sqrt(vc_eta(ind3-1,ind2,ind1)**2) + c_(ind3-1,ind2,ind1)*ijacob_irs(ind3-1,ind2)
               !rspecph = sqrt(vc_ksi(ind3,ind2,ind1)**2)   + sqrt(vc_eta(ind3,ind2,ind1)**2)   + c_(ind3,ind2,ind1)*ijacob_irs(ind3,ind2)
               rspecmh = sqrt(vc_ksi(ind3-1,ind2,ind1)**2 + vc_eta(ind3-1,ind2,ind1)**2) + c_(ind3-1,ind2,ind1)*(sqrt(n2_ksi(ind3-1,ind2) + n2_eta(ind3-1,ind2)))
               rspecph = sqrt(vc_ksi(ind3,ind2,ind1)**2   + vc_eta(ind3,ind2,ind1)**2)   + c_(ind3,ind2,ind1)*(  sqrt(n2_ksi(ind3,ind2)   + n2_eta(ind3,ind2)))
               r_spec4(ind3,ind2,ind1) = (0.5 * (rspecmh + rspecph))**4
               ! r_spec4(ind3,ind2,ind1) = MAX(rspecmh,rspecph)**4
            enddo
         enddo
      enddo

      ! Façon Martinelli
      ! ----------------
      ! do ind1=ind1min,ind1max
      !    do ind2=ind2min+1,ind2max-1
      !       do ind3=ind3min+1,ind3max
      !          ! rspecmh = sqrt(vc_ksi(ind3,ind2,ind1)**2 + vc_eta(ind3,ind2,ind1)**2) + c_(ind3,ind2,ind1)*(sqrt(n2_ksi(ind3,ind2) + n2_eta(ind3,ind2)))
      !          ! rspecph = sqrt(vc_ksi(ind3+1,ind2,ind1)**2   + vc_eta(ind3+1,ind2,ind1)**2)   + c_(ind3+1,ind2,ind1)*(  sqrt(n2_ksi(ind3+1,ind2)   + n2_eta(ind3+1,ind2)))
      !          rspecmh2 = 0.5_wp*( &
      !                   sqrt(vc_ksi(ind3-1,ind2-1,ind1)**2 + vc_eta(ind3-1,ind2-1,ind1)**2) + c_(ind3-1,ind2-1,ind1)*(sqrt(n2_ksi(ind3-1,ind2-1) + n2_eta(ind3-1,ind2-1))) &
      !                 + sqrt(vc_ksi(ind3-1,ind2,ind1)**2   + vc_eta(ind3-1,ind2,ind1)**2)   + c_(ind3-1,ind2,ind1)*  (sqrt(n2_ksi(ind3-1,ind2)   + n2_eta(ind3-1,ind2))))
      !          rspecph2 = 0.5_wp*( &
      !                   sqrt(vc_ksi(ind3-1,ind2,ind1)**2   + vc_eta(ind3-1,ind2,ind1)**2)   + c_(ind3-1,ind2,ind1)*  (sqrt(n2_ksi(ind3-1,ind2)   + n2_eta(ind3-1,ind2))) &
      !                 + sqrt(vc_ksi(ind3-1,ind2+1,ind1)**2 + vc_eta(ind3-1,ind2+1,ind1)**2) + c_(ind3-1,ind2+1,ind1)*(sqrt(n2_ksi(ind3-1,ind2+1) + n2_eta(ind3-1,ind2+1))))
      !          rspecmh3 = 0.5_wp*( &
      !                   sqrt(vc_ksi(ind3,ind2-1,ind1)**2 + vc_eta(ind3,ind2-1,ind1)**2) + c_(ind3,ind2-1,ind1)*(sqrt(n2_ksi(ind3,ind2-1) + n2_eta(ind3,ind2-1))) &
      !                 + sqrt(vc_ksi(ind3,ind2,ind1)**2   + vc_eta(ind3,ind2,ind1)**2)   + c_(ind3,ind2,ind1)*  (sqrt(n2_ksi(ind3,ind2)   + n2_eta(ind3,ind2))))
      !          rspecph3 = 0.5_wp*( &
      !                   sqrt(vc_ksi(ind3,ind2,ind1)**2 + vc_eta(ind3,ind2,ind1)**2) + c_(ind3,ind2,ind1)*(sqrt(n2_ksi(ind3,ind2) + n2_eta(ind3,ind2))) &
      !                 + sqrt(vc_ksi(ind3,ind2+1,ind1)**2 + vc_eta(ind3,ind2+1,ind1)**2) + c_(ind3,ind2+1,ind1)*(sqrt(n2_ksi(ind3,ind2+1) + n2_eta(ind3,ind2+1))))

      !          r_spec4_2(ind3,ind2,ind1) = (r_spec4(ind3,ind2,ind1)**0.25_wp + 0.2_wp*MAX(rspecmh2,rspecph2,rspecmh3,rspecph3)**0.7_wp)**4
      !       enddo
      !    enddo
      ! enddo

      ! do ind1=ind1min,ind1max
      !    ind2 = ind2min
      !    do ind3=ind3min+1,ind3max
      !       r_spec4_2(ind3,ind2,ind1) = r_spec4(ind3,ind2,ind1)
      !    enddo
      !    ind2 = ind2max
      !    do ind3=ind3min+1,ind3max
      !       r_spec4_2(ind3,ind2,ind1) = r_spec4(ind3,ind2,ind1)
      !    enddo
      ! enddo
      ! =========
      ! fin test
      ! =========

      do ind1 = ind1min, ind1max
         ! Construction of the pentadiagonal matrix
         ! The linear systems is Mx = y with :
         !     _M = (m_ij) for i,j in {1,..,n} and with :
         !                 m_ii = d_i for i in {1,..,n}
         !                 m_ii-2 = e_i-2 for i in {3,..,n}
         !                 m_ii-1 = c_i-1 for i in {2,..,n}
         !                 m_ii+2 = b_i for i in {1,..,n-2}
         !                 m_ii+1 = a_i for i in {1,..,n-1}
         !     _ y = RHS = dwi
         !     _d and y are of dimensions n, a and c of dimensions n-1,
         !           e and b of dimensions n-2
         ind3=ind3min
         do ind2 = ind2min, ind2max
            coef     = theta_irs1 * deltat
            array_bksi(ind2,ind3)  = 0.
            array_aksi(ind2,ind3)  = - coef * (r_spec4(ind3+1,ind2,ind1))**0.25
            array_dksi(ind2,ind3)  = 1. + coef * (r_spec4(ind3+1,ind2,ind1))**0.25
         end do

         ind3=ind3max
         do ind2 = ind2min, ind2max
            coef     = theta_irs1 * deltat
            array_dksi(ind2,ind3)    = 1. + coef * (r_spec4(ind3,ind2,ind1))**0.25
            array_cksi(ind2,ind3-1)  = - coef * (r_spec4(ind3,ind2,ind1))**0.25
            array_eksi(ind2,ind3-2)  = 0.
         end do

         ind3=ind3min+1
         do ind2 = ind2min, ind2max
            coef     = - theta_irs2 * deltat**2
            array_bksi(ind2,ind3)    = 0.
            array_aksi(ind2,ind3)    = coef * (r_spec4(ind3+1,ind2,ind1))**0.5
            array_cksi(ind2,ind3-1)  = coef * (r_spec4(ind3,ind2,ind1))**0.5
            array_dksi(ind2,ind3)    = 1. - (array_aksi(ind2,ind3) + array_cksi(ind2,ind3-1))
         end do

         ind3=ind3max-1
         do ind2 = ind2min, ind2max
            coef     = - theta_irs2 * deltat**2
            array_aksi(ind2,ind3)    = coef * (r_spec4(ind3+1,ind2,ind1))**0.5
            array_cksi(ind2,ind3-1)  = coef * (r_spec4(ind3,ind2,ind1))**0.5
            array_dksi(ind2,ind3)    = 1. - (array_aksi(ind2,ind3) + array_cksi(ind2,ind3-1))
            array_eksi(ind2,ind3-2)  = 0.
         end do

         do ind3 = ind3min+2,ind3max-2
            do ind2 = ind2min, ind2max
               coef     = theta_irs4 * deltat**4
               array_eksi(ind2,ind3-2)  = coef * r_spec4(ind3,ind2,ind1)
               array_bksi(ind2,ind3)    = coef * r_spec4(ind3+1,ind2,ind1)
               array_cksi(ind2,ind3-1)  =  - 3.0_wp * array_eksi(ind2,ind3-2) - array_bksi(ind2,ind3)
               array_aksi(ind2,ind3)    =  - 3.0_wp * array_bksi(ind2,ind3) - array_eksi(ind2,ind3-2)
               array_dksi(ind2,ind3)    = 1.0_wp + 3.0_wp * (array_eksi(ind2,ind3-2) + array_bksi(ind2,ind3))
            end do
         enddo

         ! Filling of the RHS
         do ind3 = ind3min,ind3max
            do ind2 = ind2min,ind2max
                RHS_ksi(ind2,ind3) =  Krho(ind3,ind2,ind1)
               RHSu_ksi(ind2,ind3) = Krhou(ind3,ind2,ind1)
               RHSv_ksi(ind2,ind3) = Krhov(ind3,ind2,ind1)
               RHSw_ksi(ind2,ind3) = Krhow(ind3,ind2,ind1)
               RHSe_ksi(ind2,ind3) = Krhoe(ind3,ind2,ind1)
            end do
         end do

         ! Resolution of the pentadiagonal linear system
         ! =============================================
         ind3 = ind3min
         do ind2 = ind2min, ind2max
            temp_ksi1(ind2) = 1 / array_dksi(ind2,ind3)
            array_aksi(ind2,ind3) = array_aksi(ind2,ind3) * temp_ksi1(ind2)
            array_bksi(ind2,ind3) = array_bksi(ind2,ind3) * temp_ksi1(ind2)
             RHS_ksi(ind2,ind3)    =  RHS_ksi(ind2,ind3) * temp_ksi1(ind2)
            RHSu_ksi(ind2,ind3)    = RHSu_ksi(ind2,ind3) * temp_ksi1(ind2)
            RHSv_ksi(ind2,ind3)    = RHSv_ksi(ind2,ind3) * temp_ksi1(ind2)
            RHSw_ksi(ind2,ind3)    = RHSw_ksi(ind2,ind3) * temp_ksi1(ind2)
            RHSe_ksi(ind2,ind3)    = RHSe_ksi(ind2,ind3) * temp_ksi1(ind2)
         end do

         ind3 = ind3min + 1
         do ind2 = ind2min, ind2max
            temp_ksi2(ind2) = array_cksi(ind2,ind3-1)
            temp_ksi1(ind2) = 1. / (array_dksi(ind2,ind3) - array_aksi(ind2,ind3-1)*temp_ksi2(ind2))
            array_aksi(ind2,ind3) = (array_aksi(ind2,ind3) - array_bksi(ind2,ind3-1)*temp_ksi2(ind2)) * temp_ksi1(ind2)
            array_bksi(ind2,ind3) = array_bksi(ind2,ind3) * temp_ksi1(ind2)
             RHS_ksi(ind2,ind3) = ( RHS_ksi(ind2,ind3) -  RHS_ksi(ind2,ind3-1)*temp_ksi2(ind2)) * temp_ksi1(ind2)
            RHSu_ksi(ind2,ind3) = (RHSu_ksi(ind2,ind3) - RHSu_ksi(ind2,ind3-1)*temp_ksi2(ind2)) * temp_ksi1(ind2)
            RHSv_ksi(ind2,ind3) = (RHSv_ksi(ind2,ind3) - RHSv_ksi(ind2,ind3-1)*temp_ksi2(ind2)) * temp_ksi1(ind2)
            RHSw_ksi(ind2,ind3) = (RHSw_ksi(ind2,ind3) - RHSw_ksi(ind2,ind3-1)*temp_ksi2(ind2)) * temp_ksi1(ind2)
            RHSe_ksi(ind2,ind3) = (RHSe_ksi(ind2,ind3) - RHSe_ksi(ind2,ind3-1)*temp_ksi2(ind2)) * temp_ksi1(ind2)
         end do

         ! Step 5: for i = 3,...,n-2
         do ind3 = ind3min+2, ind3max-2
            do ind2 = ind2min, ind2max
               temp_ksi2(ind2) = array_cksi(ind2,ind3-1) - array_aksi(ind2,ind3-2)*array_eksi(ind2,ind3-2)
               temp_ksi1(ind2) = 1. / (array_dksi(ind2,ind3) - array_bksi(ind2,ind3-2)*array_eksi(ind2,ind3-2) - array_aksi(ind2,ind3-1)*temp_ksi2(ind2))
               array_aksi(ind2,ind3) = (array_aksi(ind2,ind3) - array_bksi(ind2,ind3-1)*temp_ksi2(ind2)) * temp_ksi1(ind2)
               array_bksi(ind2,ind3) = array_bksi(ind2,ind3) * temp_ksi1(ind2)
                RHS_ksi(ind2,ind3) = ( RHS_ksi(ind2,ind3) -  RHS_ksi(ind2,ind3-2)*array_eksi(ind2,ind3-2) -  RHS_ksi(ind2,ind3-1)*temp_ksi2(ind2)) * temp_ksi1(ind2)
               RHSu_ksi(ind2,ind3) = (RHSu_ksi(ind2,ind3) - RHSu_ksi(ind2,ind3-2)*array_eksi(ind2,ind3-2) - RHSu_ksi(ind2,ind3-1)*temp_ksi2(ind2)) * temp_ksi1(ind2)
               RHSv_ksi(ind2,ind3) = (RHSv_ksi(ind2,ind3) - RHSv_ksi(ind2,ind3-2)*array_eksi(ind2,ind3-2) - RHSv_ksi(ind2,ind3-1)*temp_ksi2(ind2)) * temp_ksi1(ind2)
               RHSw_ksi(ind2,ind3) = (RHSw_ksi(ind2,ind3) - RHSw_ksi(ind2,ind3-2)*array_eksi(ind2,ind3-2) - RHSw_ksi(ind2,ind3-1)*temp_ksi2(ind2)) * temp_ksi1(ind2)
               RHSe_ksi(ind2,ind3) = (RHSe_ksi(ind2,ind3) - RHSe_ksi(ind2,ind3-2)*array_eksi(ind2,ind3-2) - RHSe_ksi(ind2,ind3-1)*temp_ksi2(ind2)) * temp_ksi1(ind2)
            end do
         end do

         ! Step 5: for i = n-1, n
         ind3 = ind3max - 1
         do ind2 = ind2min, ind2max
            temp_ksi2(ind2) = array_cksi(ind2,ind3-1) - array_aksi(ind2,ind3-2)*array_eksi(ind2,ind3-2)
            temp_ksi1(ind2) = 1. / (array_dksi(ind2,ind3) - array_bksi(ind2,ind3-2)*array_eksi(ind2,ind3-2) - array_aksi(ind2,ind3-1)*temp_ksi2(ind2))
            array_aksi(ind2,ind3) = (array_aksi(ind2,ind3) - array_bksi(ind2,ind3-1)*temp_ksi2(ind2)) * temp_ksi1(ind2)
             RHS_ksi(ind2,ind3) = ( RHS_ksi(ind2,ind3) -  RHS_ksi(ind2,ind3-2)*array_eksi(ind2,ind3-2) -  RHS_ksi(ind2,ind3-1)*temp_ksi2(ind2)) * temp_ksi1(ind2)
            RHSu_ksi(ind2,ind3) = (RHSu_ksi(ind2,ind3) - RHSu_ksi(ind2,ind3-2)*array_eksi(ind2,ind3-2) - RHSu_ksi(ind2,ind3-1)*temp_ksi2(ind2)) * temp_ksi1(ind2)
            RHSv_ksi(ind2,ind3) = (RHSv_ksi(ind2,ind3) - RHSv_ksi(ind2,ind3-2)*array_eksi(ind2,ind3-2) - RHSv_ksi(ind2,ind3-1)*temp_ksi2(ind2)) * temp_ksi1(ind2)
            RHSw_ksi(ind2,ind3) = (RHSw_ksi(ind2,ind3) - RHSw_ksi(ind2,ind3-2)*array_eksi(ind2,ind3-2) - RHSw_ksi(ind2,ind3-1)*temp_ksi2(ind2)) * temp_ksi1(ind2)
            RHSe_ksi(ind2,ind3) = (RHSe_ksi(ind2,ind3) - RHSe_ksi(ind2,ind3-2)*array_eksi(ind2,ind3-2) - RHSe_ksi(ind2,ind3-1)*temp_ksi2(ind2)) * temp_ksi1(ind2)
         end do

         ind3 = ind3max
         do ind2 = ind2min, ind2max
            temp_ksi2(ind2) = array_cksi(ind2,ind3-1) - array_aksi(ind2,ind3-2)*array_eksi(ind2,ind3-2)
            temp_ksi1(ind2) = 1. / (array_dksi(ind2,ind3) - array_bksi(ind2,ind3-2)*array_eksi(ind2,ind3-2) - array_aksi(ind2,ind3-1)*temp_ksi2(ind2))
             RHS_ksi(ind2,ind3) = ( RHS_ksi(ind2,ind3) -  RHS_ksi(ind2,ind3-2)*array_eksi(ind2,ind3-2) -  RHS_ksi(ind2,ind3-1)*temp_ksi2(ind2)) * temp_ksi1(ind2)
            RHSu_ksi(ind2,ind3) = (RHSu_ksi(ind2,ind3) - RHSu_ksi(ind2,ind3-2)*array_eksi(ind2,ind3-2) - RHSu_ksi(ind2,ind3-1)*temp_ksi2(ind2)) * temp_ksi1(ind2)
            RHSv_ksi(ind2,ind3) = (RHSv_ksi(ind2,ind3) - RHSv_ksi(ind2,ind3-2)*array_eksi(ind2,ind3-2) - RHSv_ksi(ind2,ind3-1)*temp_ksi2(ind2)) * temp_ksi1(ind2)
            RHSw_ksi(ind2,ind3) = (RHSw_ksi(ind2,ind3) - RHSw_ksi(ind2,ind3-2)*array_eksi(ind2,ind3-2) - RHSw_ksi(ind2,ind3-1)*temp_ksi2(ind2)) * temp_ksi1(ind2)
            RHSe_ksi(ind2,ind3) = (RHSe_ksi(ind2,ind3) - RHSe_ksi(ind2,ind3-2)*array_eksi(ind2,ind3-2) - RHSe_ksi(ind2,ind3-1)*temp_ksi2(ind2)) * temp_ksi1(ind2)
         end do


         ! Step 6: computation of the solution vector
         ind3 = ind3max - 1
         do ind2 = ind2min, ind2max
             RHS_ksi(ind2,ind3) =  RHS_ksi(ind2,ind3) - array_aksi(ind2,ind3)* RHS_ksi(ind2,ind3+1)
            RHSu_ksi(ind2,ind3) = RHSu_ksi(ind2,ind3) - array_aksi(ind2,ind3)*RHSu_ksi(ind2,ind3+1)
            RHSv_ksi(ind2,ind3) = RHSv_ksi(ind2,ind3) - array_aksi(ind2,ind3)*RHSv_ksi(ind2,ind3+1)
            RHSw_ksi(ind2,ind3) = RHSw_ksi(ind2,ind3) - array_aksi(ind2,ind3)*RHSw_ksi(ind2,ind3+1)
            RHSe_ksi(ind2,ind3) = RHSe_ksi(ind2,ind3) - array_aksi(ind2,ind3)*RHSe_ksi(ind2,ind3+1)
         end do

         do ind3 = ind3max-2,ind3min,-1
            do ind2 = ind2min, ind2max
                RHS_ksi(ind2,ind3) =  RHS_ksi(ind2,ind3) - array_aksi(ind2,ind3)* RHS_ksi(ind2,ind3+1)  - array_bksi(ind2,ind3)* RHS_ksi(ind2,ind3+2)
               RHSu_ksi(ind2,ind3) = RHSu_ksi(ind2,ind3) - array_aksi(ind2,ind3)*RHSu_ksi(ind2,ind3+1)  - array_bksi(ind2,ind3)*RHSu_ksi(ind2,ind3+2)
               RHSv_ksi(ind2,ind3) = RHSv_ksi(ind2,ind3) - array_aksi(ind2,ind3)*RHSv_ksi(ind2,ind3+1)  - array_bksi(ind2,ind3)*RHSv_ksi(ind2,ind3+2)
               RHSw_ksi(ind2,ind3) = RHSw_ksi(ind2,ind3) - array_aksi(ind2,ind3)*RHSw_ksi(ind2,ind3+1)  - array_bksi(ind2,ind3)*RHSw_ksi(ind2,ind3+2)
               RHSe_ksi(ind2,ind3) = RHSe_ksi(ind2,ind3) - array_aksi(ind2,ind3)*RHSe_ksi(ind2,ind3+1)  - array_bksi(ind2,ind3)*RHSe_ksi(ind2,ind3+2)
            end do
         end do

         ! Update with the solution
         do ind3 = ind3min,ind3max
            do ind2 = ind2min,ind2max
                Krho(ind3,ind2,ind1) =  RHS_ksi(ind2,ind3)
               Krhou(ind3,ind2,ind1) = RHSu_ksi(ind2,ind3)
               Krhov(ind3,ind2,ind1) = RHSv_ksi(ind2,ind3)
               Krhow(ind3,ind2,ind1) = RHSw_ksi(ind2,ind3)
               Krhoe(ind3,ind2,ind1) = RHSe_ksi(ind2,ind3)
            end do
         end do
      enddo

   end if

   !----------------------------------------
   ! Implicitation of direction eta
   !----------------------------------------
   if (is_irs_j) then
!      ind3 : j;     ind1 : k;     ind2 : i
      ind1min = 1;   ind2min = 1;   ind3min = ndy - ngh
      ind1max = nz;  ind2max = nx;  ind3max = nfy + ngh

      ! If wall, do not take into account the wall point
      ! ================================================
      if (BC_face(1,1)%sort==0) ind2min = 2
      if (BC_face(1,2)%sort==0) ind2max = nx - 1
      if (BC_face(3,1)%sort==0) ind1min = 2
      if (BC_face(3,2)%sort==0) ind1max = nz - 1

      ! Determination of coefficients
      ! =============================
      ! coef_j = 0.0_wp
      ! do ind3=ind3min+5,ind3max-5
      !    coef_j(ind3) = theta_irs4 * deltat**4
      ! enddo
      ! if (BC_face(2,1)%sort==0) then
      !    ! coef_j(ind3min)   = theta_irs1 * deltat
      !    coef_j(ind3min)   = 0.0_wp
      !    coef_j(ind3min+1) = - 1.0*theta_irs2 * deltat**2
      !    coef_j(ind3min+2) = 1.0*theta_irs4 * deltat**4
      !    coef_j(ind3min+3) = 1.0*theta_irs4 * deltat**4
      !    coef_j(ind3min+4) = 1.0*theta_irs4 * deltat**4
      ! else
      !    coef_j(ind3min)   = 1*theta_irs1 * deltat
      !    coef_j(ind3min+1) = - 1.0*theta_irs2 * deltat**2
      !    coef_j(ind3min+2) = 1.0*theta_irs4 * deltat**4
      !    coef_j(ind3min+3) = 1.0*theta_irs4 * deltat**4
      !    coef_j(ind3min+4) = 1.0*theta_irs4 * deltat**4
      ! endif
      ! if (BC_face(2,2)%sort==0) then
      !    ! coef_j(ind3max)   = theta_irs1 * deltat
      !    coef_j(ind3max)   = 0.0_wp
      !    coef_j(ind3max-1) = - theta_irs2 * deltat**2
      !    coef_j(ind3max-2) = theta_irs4 * deltat**4
      !    coef_j(ind3max-3) = theta_irs4 * deltat**4
      !    coef_j(ind3max-4) = theta_irs4 * deltat**4
      ! else
      !    coef_j(ind3max)   = 1.0*theta_irs1 * deltat
      !    coef_j(ind3max-1) = - 1.0*theta_irs2 * deltat**2
      !    coef_j(ind3max-2) = 1.0*theta_irs4 * deltat**4
      !    coef_j(ind3max-3) = 1.0*theta_irs4 * deltat**4
      !    coef_j(ind3max-4) = 1.0*theta_irs4 * deltat**4
      ! endif

      ! Communication of ghost points
      ! =============================
      ! call communication_NS(Krho,Krhou,Krhov,Krhow,Krhoe)
      call communication_(Krho,Krhou,Krhov,Krhow,Krhoe)


      ! Calculation of contra-variant speed and norm of eta
      ! ===================================================
      do ind3=ind3min,ind3max
         do ind2=ind2min,ind2max
            n2_eta(ind2,ind3) = (y_ksi(ind2,ind3)**2 + x_ksi(ind2,ind3)**2) * ijacob_irs(ind2,ind3)**2
         end do
      end do

      do ind1=ind1min,ind1max
         do ind3=ind3min,ind3max
            do ind2=ind2min,ind2max
               vc_eta(ind2,ind3,ind1) = (vv(ind2,ind3,ind1)*x_ksi(ind2,ind3) - uu(ind2,ind3,ind1)*y_ksi(ind2,ind3)) * ijacob_irs(ind2,ind3)
            end do
         end do
      end do

!!$      ! Calculation of the spectral radius at ind-1/2
!!$      ! =============================================
!!$      do ind1=ind1min,ind1max
!!$         do ind3=ind3min+1,ind3max
!!$            do ind2=ind2min,ind2max
!!$               rspecmh = sqrt(vc_eta(ind2,ind3-1,ind1)**2) + c_(ind2,ind3-1,ind1)*sqrt(n2_eta(ind2,ind3-1))
!!$               rspecph = sqrt(vc_eta(ind2,ind3,ind1)**2)   + c_(ind2,ind3,ind1)*sqrt(n2_eta(ind2,ind3))
!!$               r_spec4(ind2,ind3,ind1) = (0.5 * (rspecmh + rspecph))**4
!!$            enddo
!!$         enddo
!!$      enddo

      ! ===========
      !    Test
      ! ===========
      ! A ENLEVER
      ! Calculation of contra-variant speed and norm of ksi
      ! ===================================================
      do ind3=ind3min,ind3max
         do ind2=ind2min,ind2max
            n2_ksi(ind2,ind3) = (y_eta(ind2,ind3)**2 + x_eta(ind2,ind3)**2) * ijacob_irs(ind2,ind3)**2
         end do
      end do

      do ind1=ind1min,ind1max
         do ind2=ind2min,ind2max
            do ind3=ind3min,ind3max
               vc_ksi(ind2,ind3,ind1) = (uu(ind2,ind3,ind1)*y_eta(ind2,ind3)-vv(ind2,ind3,ind1)*x_eta(ind2,ind3)) * ijacob_irs(ind2,ind3)
            end do
         end do
      end do

      ! Calculation of the spectral radius at ind-1/2
      ! =============================================
      do ind1=ind1min,ind1max
         do ind3=ind3min+1,ind3max
            do ind2=ind2min,ind2max
               !rspecmh = sqrt(vc_eta(ind2,ind3-1,ind1)**2) + sqrt(vc_ksi(ind2,ind3-1,ind1)**2) + c_(ind2,ind3-1,ind1)*(sqrt(n2_eta(ind2,ind3-1)) + sqrt(n2_ksi(ind2,ind3-1)))
               !rspecph = sqrt(vc_eta(ind2,ind3,ind1)**2)   + sqrt(vc_ksi(ind2,ind3,ind1)**2)   + c_(ind2,ind3,ind1)*(  sqrt(n2_eta(ind2,ind3)) + sqrt(n2_ksi(ind2,ind3)))
               !rspecmh = sqrt(vc_eta(ind2,ind3-1,ind1)**2 + vc_ksi(ind2,ind3-1,ind1)**2) + c_(ind2,ind3-1,ind1)*(sqrt(n2_eta(ind2,ind3-1)) + sqrt(n2_ksi(ind2,ind3-1)))
               !rspecph = sqrt(vc_eta(ind2,ind3,ind1)**2   + vc_ksi(ind2,ind3,ind1)**2)   + c_(ind2,ind3,ind1)*(  sqrt(n2_eta(ind2,ind3))   + sqrt(n2_ksi(ind2,ind3)))
               !rspecmh = sqrt(vc_eta(ind2,ind3-1,ind1)**2) + sqrt(vc_ksi(ind2,ind3-1,ind1)**2) + c_(ind2,ind3-1,ind1)*ijacob_irs(ind2,ind3-1)
               !rspecph = sqrt(vc_eta(ind2,ind3,ind1)**2)   + sqrt(vc_ksi(ind2,ind3,ind1)**2)   + c_(ind2,ind3,ind1)*ijacob_irs(ind2,ind3)
               rspecmh = sqrt(vc_eta(ind2,ind3-1,ind1)**2 + vc_ksi(ind2,ind3-1,ind1)**2) + c_(ind2,ind3-1,ind1)*(sqrt(n2_eta(ind2,ind3-1) + n2_ksi(ind2,ind3-1)))
               rspecph = sqrt(vc_eta(ind2,ind3,ind1)**2   + vc_ksi(ind2,ind3,ind1)**2)   + c_(ind2,ind3,ind1)*(  sqrt(n2_eta(ind2,ind3)   + n2_ksi(ind2,ind3)))
               r_spec4(ind2,ind3,ind1) = (0.5 * (rspecmh + rspecph))**4
               ! r_spec4(ind2,ind3,ind1) = MAX(rspecmh,rspecph)**4
            enddo
         enddo
      enddo

      ! Façon Martinelli
      ! ----------------
      ! do ind1=ind1min,ind1max
      !    do ind3=ind3min+1,ind3max
      !       do ind2=ind2min+1,ind2max-1
      !          rspecmh2 = 0.5_wp*(&
      !                     sqrt(vc_eta(ind2-1,ind3-1,ind1)**2 + vc_ksi(ind2-1,ind3-1,ind1)**2) + c_(ind2-1,ind3-1,ind1)*(sqrt(n2_eta(ind2-1,ind3-1) + n2_ksi(ind2-1,ind3-1))) &
      !                   + sqrt(vc_eta(ind2,  ind3-1,ind1)**2 + vc_ksi(ind2,  ind3-1,ind1)**2) + c_(ind2,  ind3-1,ind1)*(sqrt(n2_eta(ind2,  ind3-1) + n2_ksi(ind2,  ind3-1))))
      !          rspecph2 = 0.5_wp*(&
      !                     sqrt(vc_eta(ind2+1,ind3-1,ind1)**2 + vc_ksi(ind2+1,ind3-1,ind1)**2) + c_(ind2+1,ind3-1,ind1)*(sqrt(n2_eta(ind2+1,ind3-1) + n2_ksi(ind2+1,ind3-1))) &
      !                   + sqrt(vc_eta(ind2,  ind3-1,ind1)**2 + vc_ksi(ind2,  ind3-1,ind1)**2) + c_(ind2,  ind3-1,ind1)*(sqrt(n2_eta(ind2,  ind3-1) + n2_ksi(ind2,  ind3-1))))
      !          rspecmh3 = 0.5_wp*(&
      !                     sqrt(vc_eta(ind2-1,ind3,ind1)**2 + vc_ksi(ind2-1,ind3,ind1)**2) + c_(ind2-1,ind3,ind1)*(sqrt(n2_eta(ind2-1,ind3) + n2_ksi(ind2-1,ind3))) &
      !                   + sqrt(vc_eta(ind2,  ind3,ind1)**2 + vc_ksi(ind2,  ind3,ind1)**2) + c_(ind2,  ind3,ind1)*(sqrt(n2_eta(ind2,  ind3) + n2_ksi(ind2,  ind3))))
      !          rspecph3 = 0.5_wp*(&
      !                     sqrt(vc_eta(ind2+1,ind3,ind1)**2 + vc_ksi(ind2+1,ind3,ind1)**2) + c_(ind2+1,ind3,ind1)*(sqrt(n2_eta(ind2+1,ind3) + n2_ksi(ind2+1,ind3))) &
      !                   + sqrt(vc_eta(ind2,  ind3,ind1)**2 + vc_ksi(ind2,  ind3,ind1)**2) + c_(ind2,  ind3,ind1)*(sqrt(n2_eta(ind2,  ind3) + n2_ksi(ind2,  ind3))))
      !          r_spec4_2(ind2,ind3,ind1) = (r_spec4(ind2,ind3,ind1)**0.25_wp + 0.2_wp*MAX(rspecmh2,rspecph2)**0.2_wp)**4
      !       enddo
      !    enddo
      ! enddo

      ! do ind1=ind1min,ind1max
      !    ind2 = ind2min
      !    do ind3=ind3min+1,ind3max
      !       r_spec4_2(ind2,ind3,ind1) = r_spec4(ind2,ind3,ind1)
      !    enddo
      !    ind2 = ind2max
      !    do ind3=ind3min+1,ind3max
      !       r_spec4_2(ind2,ind3,ind1) = r_spec4(ind2,ind3,ind1)
      !    enddo
      ! enddo

      ! =========
      ! fin test
      ! =========

      do ind1 = ind1min, ind1max
         ! Construction of the pentadiagonal matrix
         ! ========================================
         ind3=ind3min
         do ind2 = ind2min, ind2max
            coef     = theta_irs1 * deltat
            array_by(ind2,ind3)  = 0.
            array_aeta(ind2,ind3)  = - coef * (r_spec4(ind2,ind3+1,ind1))**0.25
            array_deta(ind2,ind3)  = 1. + coef * (r_spec4(ind2,ind3+1,ind1))**0.25
         end do

         ind3=ind3max
         do ind2 = ind2min, ind2max
            coef     = theta_irs1 * deltat
            array_deta(ind2,ind3)    = 1. + coef * (r_spec4(ind2,ind3,ind1))**0.25
            array_ceta(ind2,ind3-1)  = - coef * (r_spec4(ind2,ind3,ind1))**0.25
            array_ey(ind2,ind3-2)  = 0.
         end do

         ind3=ind3min+1
         do ind2 = ind2min, ind2max
            coef     = - theta_irs2 * deltat**2
            array_by(ind2,ind3)    = 0.
            array_aeta(ind2,ind3)    = coef * (r_spec4(ind2,ind3+1,ind1))**0.5
            array_ceta(ind2,ind3-1)  = coef * (r_spec4(ind2,ind3,ind1))**0.5
            array_deta(ind2,ind3)    = 1. - (array_aeta(ind2,ind3) + array_ceta(ind2,ind3-1))
         end do

         ind3=ind3max-1
         do ind2 = ind2min, ind2max
            coef     = - theta_irs2 * deltat**2
            array_aeta(ind2,ind3)    = coef * (r_spec4(ind2,ind3+1,ind1))**0.5
            array_ceta(ind2,ind3-1)  = coef * (r_spec4(ind2,ind3,ind1))**0.5
            array_deta(ind2,ind3)    = 1. - (array_aeta(ind2,ind3) + array_ceta(ind2,ind3-1))
            array_ey(ind2,ind3-2)  = 0.
         end do

         do ind3 = ind3min+2,ind3max-2
            do ind2 = ind2min, ind2max
               coef     = theta_irs4 * deltat**4
               array_ey(ind2,ind3-2)   = coef * r_spec4(ind2,ind3,ind1)
               array_by(ind2,ind3)     = coef * r_spec4(ind2,ind3+1,ind1)
               array_ceta(ind2,ind3-1) =  - 3.0_wp * array_ey(ind2,ind3-2) - array_by(ind2,ind3)
               array_aeta(ind2,ind3)   =  - 3.0_wp * array_by(ind2,ind3) - array_ey(ind2,ind3-2)
               array_deta(ind2,ind3)   = 1.0_wp + 3.0_wp * (array_ey(ind2,ind3-2) + array_by(ind2,ind3))
            end do
         enddo

         ! Resolution of the pentadiagonal linear system
         ind3 = ind3min
         do ind2 = ind2min, ind2max
            temp_eta1(ind2) = 1 / array_deta(ind2,ind3)
            array_aeta(ind2,ind3) = array_aeta(ind2,ind3) * temp_eta1(ind2)
            array_by(ind2,ind3) = array_by(ind2,ind3) * temp_eta1(ind2)
             Krho(ind2,ind3,ind1)    =  Krho(ind2,ind3,ind1) * temp_eta1(ind2)
            Krhou(ind2,ind3,ind1)    = Krhou(ind2,ind3,ind1) * temp_eta1(ind2)
            Krhov(ind2,ind3,ind1)    = Krhov(ind2,ind3,ind1) * temp_eta1(ind2)
            Krhow(ind2,ind3,ind1)    = Krhow(ind2,ind3,ind1) * temp_eta1(ind2)
            Krhoe(ind2,ind3,ind1)    = Krhoe(ind2,ind3,ind1) * temp_eta1(ind2)
         end do

         ind3 = ind3min + 1
         do ind2 = ind2min, ind2max
            temp_eta2(ind2) = array_ceta(ind2,ind3-1)
            temp_eta1(ind2) = 1. / (array_deta(ind2,ind3) - array_aeta(ind2,ind3-1)*temp_eta2(ind2))
            array_aeta(ind2,ind3) = (array_aeta(ind2,ind3) - array_by(ind2,ind3-1)*temp_eta2(ind2)) * temp_eta1(ind2)
            array_by(ind2,ind3) = array_by(ind2,ind3) * temp_eta1(ind2)
             Krho(ind2,ind3,ind1) = ( Krho(ind2,ind3,ind1) -  Krho(ind2,ind3-1,ind1)*temp_eta2(ind2)) * temp_eta1(ind2)
            Krhou(ind2,ind3,ind1) = (Krhou(ind2,ind3,ind1) - Krhou(ind2,ind3-1,ind1)*temp_eta2(ind2)) * temp_eta1(ind2)
            Krhov(ind2,ind3,ind1) = (Krhov(ind2,ind3,ind1) - Krhov(ind2,ind3-1,ind1)*temp_eta2(ind2)) * temp_eta1(ind2)
            Krhow(ind2,ind3,ind1) = (Krhow(ind2,ind3,ind1) - Krhow(ind2,ind3-1,ind1)*temp_eta2(ind2)) * temp_eta1(ind2)
            Krhoe(ind2,ind3,ind1) = (Krhoe(ind2,ind3,ind1) - Krhoe(ind2,ind3-1,ind1)*temp_eta2(ind2)) * temp_eta1(ind2)
         end do

         ! Step 5: for i = 3,...,n-2
         do ind3 = ind3min+2, ind3max-2
            do ind2 = ind2min, ind2max
               temp_eta2(ind2) = array_ceta(ind2,ind3-1) - array_aeta(ind2,ind3-2)*array_ey(ind2,ind3-2)
               temp_eta1(ind2) = 1. / (array_deta(ind2,ind3) - array_by(ind2,ind3-2)*array_ey(ind2,ind3-2) - array_aeta(ind2,ind3-1)*temp_eta2(ind2))
               array_aeta(ind2,ind3) = (array_aeta(ind2,ind3) - array_by(ind2,ind3-1)*temp_eta2(ind2)) * temp_eta1(ind2)
               array_by(ind2,ind3) = array_by(ind2,ind3) * temp_eta1(ind2)
                Krho(ind2,ind3,ind1) = ( Krho(ind2,ind3,ind1) -  Krho(ind2,ind3-2,ind1)*array_ey(ind2,ind3-2) -  Krho(ind2,ind3-1,ind1)*temp_eta2(ind2)) * temp_eta1(ind2)
               Krhou(ind2,ind3,ind1) = (Krhou(ind2,ind3,ind1) - Krhou(ind2,ind3-2,ind1)*array_ey(ind2,ind3-2) - Krhou(ind2,ind3-1,ind1)*temp_eta2(ind2)) * temp_eta1(ind2)
               Krhov(ind2,ind3,ind1) = (Krhov(ind2,ind3,ind1) - Krhov(ind2,ind3-2,ind1)*array_ey(ind2,ind3-2) - Krhov(ind2,ind3-1,ind1)*temp_eta2(ind2)) * temp_eta1(ind2)
               Krhow(ind2,ind3,ind1) = (Krhow(ind2,ind3,ind1) - Krhow(ind2,ind3-2,ind1)*array_ey(ind2,ind3-2) - Krhow(ind2,ind3-1,ind1)*temp_eta2(ind2)) * temp_eta1(ind2)
               Krhoe(ind2,ind3,ind1) = (Krhoe(ind2,ind3,ind1) - Krhoe(ind2,ind3-2,ind1)*array_ey(ind2,ind3-2) - Krhoe(ind2,ind3-1,ind1)*temp_eta2(ind2)) * temp_eta1(ind2)
            end do
         end do

         ! Step 5: for i = n-1, n
         ind3 = ind3max - 1
         do ind2 = ind2min, ind2max
            temp_eta2(ind2) = array_ceta(ind2,ind3-1) - array_aeta(ind2,ind3-2)*array_ey(ind2,ind3-2)
            temp_eta1(ind2) = 1. / (array_deta(ind2,ind3) - array_by(ind2,ind3-2)*array_ey(ind2,ind3-2) - array_aeta(ind2,ind3-1)*temp_eta2(ind2))
            array_aeta(ind2,ind3) = (array_aeta(ind2,ind3) - array_by(ind2,ind3-1)*temp_eta2(ind2)) * temp_eta1(ind2)
             Krho(ind2,ind3,ind1) = ( Krho(ind2,ind3,ind1) -  Krho(ind2,ind3-2,ind1)*array_ey(ind2,ind3-2) -  Krho(ind2,ind3-1,ind1)*temp_eta2(ind2)) * temp_eta1(ind2)
            Krhou(ind2,ind3,ind1) = (Krhou(ind2,ind3,ind1) - Krhou(ind2,ind3-2,ind1)*array_ey(ind2,ind3-2) - Krhou(ind2,ind3-1,ind1)*temp_eta2(ind2)) * temp_eta1(ind2)
            Krhov(ind2,ind3,ind1) = (Krhov(ind2,ind3,ind1) - Krhov(ind2,ind3-2,ind1)*array_ey(ind2,ind3-2) - Krhov(ind2,ind3-1,ind1)*temp_eta2(ind2)) * temp_eta1(ind2)
            Krhow(ind2,ind3,ind1) = (Krhow(ind2,ind3,ind1) - Krhow(ind2,ind3-2,ind1)*array_ey(ind2,ind3-2) - Krhow(ind2,ind3-1,ind1)*temp_eta2(ind2)) * temp_eta1(ind2)
            Krhoe(ind2,ind3,ind1) = (Krhoe(ind2,ind3,ind1) - Krhoe(ind2,ind3-2,ind1)*array_ey(ind2,ind3-2) - Krhoe(ind2,ind3-1,ind1)*temp_eta2(ind2)) * temp_eta1(ind2)
         end do

         ind3 = ind3max
         do ind2 = ind2min, ind2max
            temp_eta2(ind2) = array_ceta(ind2,ind3-1) - array_aeta(ind2,ind3-2)*array_ey(ind2,ind3-2)
            temp_eta1(ind2) = 1. / (array_deta(ind2,ind3) - array_by(ind2,ind3-2)*array_ey(ind2,ind3-2) - array_aeta(ind2,ind3-1)*temp_eta2(ind2))
             Krho(ind2,ind3,ind1) = ( Krho(ind2,ind3,ind1) -  Krho(ind2,ind3-2,ind1)*array_ey(ind2,ind3-2) -  Krho(ind2,ind3-1,ind1)*temp_eta2(ind2)) * temp_eta1(ind2)
            Krhou(ind2,ind3,ind1) = (Krhou(ind2,ind3,ind1) - Krhou(ind2,ind3-2,ind1)*array_ey(ind2,ind3-2) - Krhou(ind2,ind3-1,ind1)*temp_eta2(ind2)) * temp_eta1(ind2)
            Krhov(ind2,ind3,ind1) = (Krhov(ind2,ind3,ind1) - Krhov(ind2,ind3-2,ind1)*array_ey(ind2,ind3-2) - Krhov(ind2,ind3-1,ind1)*temp_eta2(ind2)) * temp_eta1(ind2)
            Krhow(ind2,ind3,ind1) = (Krhow(ind2,ind3,ind1) - Krhow(ind2,ind3-2,ind1)*array_ey(ind2,ind3-2) - Krhow(ind2,ind3-1,ind1)*temp_eta2(ind2)) * temp_eta1(ind2)
            Krhoe(ind2,ind3,ind1) = (Krhoe(ind2,ind3,ind1) - Krhoe(ind2,ind3-2,ind1)*array_ey(ind2,ind3-2) - Krhoe(ind2,ind3-1,ind1)*temp_eta2(ind2)) * temp_eta1(ind2)
         end do


         ! Step 6: computation of the solution vector
         ind3 = ind3max - 1
         do ind2 = ind2min, ind2max
             Krho(ind2,ind3,ind1) =  Krho(ind2,ind3,ind1) - array_aeta(ind2,ind3)* Krho(ind2,ind3+1,ind1)
            Krhou(ind2,ind3,ind1) = Krhou(ind2,ind3,ind1) - array_aeta(ind2,ind3)*Krhou(ind2,ind3+1,ind1)
            Krhov(ind2,ind3,ind1) = Krhov(ind2,ind3,ind1) - array_aeta(ind2,ind3)*Krhov(ind2,ind3+1,ind1)
            Krhow(ind2,ind3,ind1) = Krhow(ind2,ind3,ind1) - array_aeta(ind2,ind3)*Krhow(ind2,ind3+1,ind1)
            Krhoe(ind2,ind3,ind1) = Krhoe(ind2,ind3,ind1) - array_aeta(ind2,ind3)*Krhoe(ind2,ind3+1,ind1)
         end do

         do ind3 = ind3max-2,ind3min,-1
            do ind2 = ind2min, ind2max
                Krho(ind2,ind3,ind1) =  Krho(ind2,ind3,ind1) - array_aeta(ind2,ind3)* Krho(ind2,ind3+1,ind1)  - array_by(ind2,ind3)* Krho(ind2,ind3+2,ind1)
               Krhou(ind2,ind3,ind1) = Krhou(ind2,ind3,ind1) - array_aeta(ind2,ind3)*Krhou(ind2,ind3+1,ind1)  - array_by(ind2,ind3)*Krhou(ind2,ind3+2,ind1)
               Krhov(ind2,ind3,ind1) = Krhov(ind2,ind3,ind1) - array_aeta(ind2,ind3)*Krhov(ind2,ind3+1,ind1)  - array_by(ind2,ind3)*Krhov(ind2,ind3+2,ind1)
               Krhow(ind2,ind3,ind1) = Krhow(ind2,ind3,ind1) - array_aeta(ind2,ind3)*Krhow(ind2,ind3+1,ind1)  - array_by(ind2,ind3)*Krhow(ind2,ind3+2,ind1)
               Krhoe(ind2,ind3,ind1) = Krhoe(ind2,ind3,ind1) - array_aeta(ind2,ind3)*Krhoe(ind2,ind3+1,ind1)  - array_by(ind2,ind3)*Krhoe(ind2,ind3+2,ind1)
            end do
         end do
      enddo

   end if

   !----------------------------------------
   ! Implicitation of direction z
   !----------------------------------------
   if (is_irs_k .and. .not. is_2d) then
      !      ind3 : k;     ind1 : i;     ind2 : j
      ind1min = 1;   ind2min = 1;   ind3min = ndz - ngh
      ind1max = ny;  ind2max = nx;  ind3max = nfz + ngh

      ! If wall, do not take into account the wall point
      ! ================================================
      if (BC_face(1,1)%sort==0) ind2min = 2
      if (BC_face(1,2)%sort==0) ind2max = nx - 1
      if (BC_face(2,1)%sort==0) ind1min = 2
      if (BC_face(2,2)%sort==0) ind1max = ny - 1

      ! Communication of ghost points
      ! =============================
      ! call communication_FB(Krho,Krhou,Krhov,Krhow,Krhoe)
      call communication_(Krho,Krhou,Krhov,Krhow,Krhoe)

      ! Calculation of the spectral radius at ind-1/2
      ! =============================================
      do ind3=ind3min+1,ind3max
         do ind1=ind1min,ind1max
            do ind2=ind2min,ind2max
               rspecmh = sqrt(ww(ind2,ind1,ind3-1)**2)   + c_(ind2,ind1,ind3-1)
               rspecph = sqrt(ww(ind2,ind1,ind3)**2)     + c_(ind2,ind1,ind3)
               r_spec4(ind2,ind1,ind3) = (0.5 * (rspecmh + rspecph))**4
            enddo
         enddo
      enddo

      do ind1 = ind1min, ind1max
         ! Construction of the pentadiagonal matrix
         ! ========================================
         ind3=ind3min
         do ind2 = ind2min, ind2max
            coef     = theta_irs1 * deltat * idz_irs(ind3)
            array_bz(ind2,ind3)  = 0.
            array_az(ind2,ind3)  = - coef * (r_spec4(ind2,ind1,ind3+1))**0.25
            array_dz(ind2,ind3)  = 1. + coef * (r_spec4(ind2,ind1,ind3+1))**0.25
         end do

         ind3=ind3max
         do ind2 = ind2min, ind2max
            coef     = theta_irs1 * deltat * idz_irs(ind3)
            array_dz(ind2,ind3)    = 1. + coef * (r_spec4(ind2,ind1,ind3))**0.25
            array_cz(ind2,ind3-1)  = - coef * (r_spec4(ind2,ind1,ind3))**0.25
            array_ez(ind2,ind3-2)  = 0.
         end do

         ind3=ind3min+1
         do ind2 = ind2min, ind2max
            coef     = - theta_irs2 * (deltat * idz_irs(ind3))**2
            array_bz(ind2,ind3)    = 0.
            array_az(ind2,ind3)    = coef * (r_spec4(ind2,ind1,ind3+1))**0.5
            array_cz(ind2,ind3-1)  = coef * (r_spec4(ind2,ind1,ind3))**0.5
            array_dz(ind2,ind3)    = 1. - (array_az(ind2,ind3) + array_cz(ind2,ind3-1))
         end do

         ind3=ind3max-1
         do ind2 = ind2min, ind2max
            coef     = - theta_irs2 * (deltat * idz_irs(ind3))**2
            array_az(ind2,ind3)    = coef * (r_spec4(ind2,ind1,ind3+1))**0.5
            array_cz(ind2,ind3-1)  = coef * (r_spec4(ind2,ind1,ind3))**0.5
            array_dz(ind2,ind3)    = 1. - (array_az(ind2,ind3) + array_cz(ind2,ind3-1))
            array_ez(ind2,ind3-2)  = 0.
         end do

         do ind3 = ind3min+2,ind3max-2
            do ind2 = ind2min, ind2max
               coefmh   = deltat * idz_irs(ind3)
               coefph   = deltat * idz_irs(ind3+1)
               coef     = theta_irs4 * deltat * 0.5_wp*(idz_irs(ind3)+idz_irs(ind3+1))

               array_ez(ind2,ind3-2)  = coef * coefmh**3 * r_spec4(ind2,ind1,ind3)
               array_bz(ind2,ind3)    = coef * coefph**3 * r_spec4(ind2,ind1,ind3+1)
               array_cz(ind2,ind3-1)  =  - 3.0_wp * array_ez(ind2,ind3-2) - array_bz(ind2,ind3)
               array_az(ind2,ind3)    =  - 3.0_wp * array_bz(ind2,ind3) - array_ez(ind2,ind3-2)
               array_dz(ind2,ind3)    = 1.0_wp + 3.0_wp * (array_ez(ind2,ind3-2) + array_bz(ind2,ind3))
            enddo
         end do


         ! Filling of the RHS
         do ind3 = ind3min,ind3max
            do ind2 = ind2min,ind2max
                RHS_z(ind2,ind3) =  Krho(ind2,ind1,ind3)
               RHSu_z(ind2,ind3) = Krhou(ind2,ind1,ind3)
               RHSv_z(ind2,ind3) = Krhov(ind2,ind1,ind3)
               RHSw_z(ind2,ind3) = Krhow(ind2,ind1,ind3)
               RHSe_z(ind2,ind3) = Krhoe(ind2,ind1,ind3)
            end do
         end do

         ! Resolution of the pentadiagonal linear system
         ! =============================================
         ind3 = ind3min
         do ind2 = ind2min, ind2max
            temp_z1(ind2) = 1 / array_dz(ind2,ind3)
            array_az(ind2,ind3) = array_az(ind2,ind3) * temp_z1(ind2)
            array_bz(ind2,ind3) = array_bz(ind2,ind3) * temp_z1(ind2)
             RHS_z(ind2,ind3)    =  RHS_z(ind2,ind3) * temp_z1(ind2)
            RHSu_z(ind2,ind3)    = RHSu_z(ind2,ind3) * temp_z1(ind2)
            RHSv_z(ind2,ind3)    = RHSv_z(ind2,ind3) * temp_z1(ind2)
            RHSw_z(ind2,ind3)    = RHSw_z(ind2,ind3) * temp_z1(ind2)
            RHSe_z(ind2,ind3)    = RHSe_z(ind2,ind3) * temp_z1(ind2)
         end do

         ind3 = ind3min + 1
         do ind2 = ind2min, ind2max
            temp_z2(ind2) = array_cz(ind2,ind3-1)
            temp_z1(ind2) = 1. / (array_dz(ind2,ind3) - array_az(ind2,ind3-1)*temp_z2(ind2))
            array_az(ind2,ind3) = (array_az(ind2,ind3) - array_bz(ind2,ind3-1)*temp_z2(ind2)) * temp_z1(ind2)
            array_bz(ind2,ind3) = array_bz(ind2,ind3) * temp_z1(ind2)
             RHS_z(ind2,ind3) = ( RHS_z(ind2,ind3) -  RHS_z(ind2,ind3-1)*temp_z2(ind2)) * temp_z1(ind2)
            RHSu_z(ind2,ind3) = (RHSu_z(ind2,ind3) - RHSu_z(ind2,ind3-1)*temp_z2(ind2)) * temp_z1(ind2)
            RHSv_z(ind2,ind3) = (RHSv_z(ind2,ind3) - RHSv_z(ind2,ind3-1)*temp_z2(ind2)) * temp_z1(ind2)
            RHSw_z(ind2,ind3) = (RHSw_z(ind2,ind3) - RHSw_z(ind2,ind3-1)*temp_z2(ind2)) * temp_z1(ind2)
            RHSe_z(ind2,ind3) = (RHSe_z(ind2,ind3) - RHSe_z(ind2,ind3-1)*temp_z2(ind2)) * temp_z1(ind2)
         end do

         ! Step 5: for i = 3,...,n-2
         do ind3 = ind3min+2, ind3max-2
            do ind2 = ind2min, ind2max
               temp_z2(ind2) = array_cz(ind2,ind3-1) - array_az(ind2,ind3-2)*array_ez(ind2,ind3-2)
               temp_z1(ind2) = 1. / (array_dz(ind2,ind3) - array_bz(ind2,ind3-2)*array_ez(ind2,ind3-2) - array_az(ind2,ind3-1)*temp_z2(ind2))
               array_az(ind2,ind3) = (array_az(ind2,ind3) - array_bz(ind2,ind3-1)*temp_z2(ind2)) * temp_z1(ind2)
               array_bz(ind2,ind3) = array_bz(ind2,ind3) * temp_z1(ind2)
                RHS_z(ind2,ind3) = ( RHS_z(ind2,ind3) -  RHS_z(ind2,ind3-2)*array_ez(ind2,ind3-2) -  RHS_z(ind2,ind3-1)*temp_z2(ind2)) * temp_z1(ind2)
               RHSu_z(ind2,ind3) = (RHSu_z(ind2,ind3) - RHSu_z(ind2,ind3-2)*array_ez(ind2,ind3-2) - RHSu_z(ind2,ind3-1)*temp_z2(ind2)) * temp_z1(ind2)
               RHSv_z(ind2,ind3) = (RHSv_z(ind2,ind3) - RHSv_z(ind2,ind3-2)*array_ez(ind2,ind3-2) - RHSv_z(ind2,ind3-1)*temp_z2(ind2)) * temp_z1(ind2)
               RHSw_z(ind2,ind3) = (RHSw_z(ind2,ind3) - RHSw_z(ind2,ind3-2)*array_ez(ind2,ind3-2) - RHSw_z(ind2,ind3-1)*temp_z2(ind2)) * temp_z1(ind2)
               RHSe_z(ind2,ind3) = (RHSe_z(ind2,ind3) - RHSe_z(ind2,ind3-2)*array_ez(ind2,ind3-2) - RHSe_z(ind2,ind3-1)*temp_z2(ind2)) * temp_z1(ind2)
            end do
         end do

         ! Step 5: for i = n-1, n
         ind3 = ind3max - 1
         do ind2 = ind2min, ind2max
            temp_z2(ind2) = array_cz(ind2,ind3-1) - array_az(ind2,ind3-2)*array_ez(ind2,ind3-2)
            temp_z1(ind2) = 1. / (array_dz(ind2,ind3) - array_bz(ind2,ind3-2)*array_ez(ind2,ind3-2) - array_az(ind2,ind3-1)*temp_z2(ind2))
            array_az(ind2,ind3) = (array_az(ind2,ind3) - array_bz(ind2,ind3-1)*temp_z2(ind2)) * temp_z1(ind2)
             RHS_z(ind2,ind3) = ( RHS_z(ind2,ind3) -  RHS_z(ind2,ind3-2)*array_ez(ind2,ind3-2) -  RHS_z(ind2,ind3-1)*temp_z2(ind2)) * temp_z1(ind2)
            RHSu_z(ind2,ind3) = (RHSu_z(ind2,ind3) - RHSu_z(ind2,ind3-2)*array_ez(ind2,ind3-2) - RHSu_z(ind2,ind3-1)*temp_z2(ind2)) * temp_z1(ind2)
            RHSv_z(ind2,ind3) = (RHSv_z(ind2,ind3) - RHSv_z(ind2,ind3-2)*array_ez(ind2,ind3-2) - RHSv_z(ind2,ind3-1)*temp_z2(ind2)) * temp_z1(ind2)
            RHSw_z(ind2,ind3) = (RHSw_z(ind2,ind3) - RHSw_z(ind2,ind3-2)*array_ez(ind2,ind3-2) - RHSw_z(ind2,ind3-1)*temp_z2(ind2)) * temp_z1(ind2)
            RHSe_z(ind2,ind3) = (RHSe_z(ind2,ind3) - RHSe_z(ind2,ind3-2)*array_ez(ind2,ind3-2) - RHSe_z(ind2,ind3-1)*temp_z2(ind2)) * temp_z1(ind2)
         end do

         ind3 = ind3max
         do ind2 = ind2min, ind2max
            temp_z2(ind2) = array_cz(ind2,ind3-1) - array_az(ind2,ind3-2)*array_ez(ind2,ind3-2)
            temp_z1(ind2) = 1. / (array_dz(ind2,ind3) - array_bz(ind2,ind3-2)*array_ez(ind2,ind3-2) - array_az(ind2,ind3-1)*temp_z2(ind2))
             RHS_z(ind2,ind3) = ( RHS_z(ind2,ind3) -  RHS_z(ind2,ind3-2)*array_ez(ind2,ind3-2) -  RHS_z(ind2,ind3-1)*temp_z2(ind2)) * temp_z1(ind2)
            RHSu_z(ind2,ind3) = (RHSu_z(ind2,ind3) - RHSu_z(ind2,ind3-2)*array_ez(ind2,ind3-2) - RHSu_z(ind2,ind3-1)*temp_z2(ind2)) * temp_z1(ind2)
            RHSv_z(ind2,ind3) = (RHSv_z(ind2,ind3) - RHSv_z(ind2,ind3-2)*array_ez(ind2,ind3-2) - RHSv_z(ind2,ind3-1)*temp_z2(ind2)) * temp_z1(ind2)
            RHSw_z(ind2,ind3) = (RHSw_z(ind2,ind3) - RHSw_z(ind2,ind3-2)*array_ez(ind2,ind3-2) - RHSw_z(ind2,ind3-1)*temp_z2(ind2)) * temp_z1(ind2)
            RHSe_z(ind2,ind3) = (RHSe_z(ind2,ind3) - RHSe_z(ind2,ind3-2)*array_ez(ind2,ind3-2) - RHSe_z(ind2,ind3-1)*temp_z2(ind2)) * temp_z1(ind2)
         end do


         ! Step 6: computation of the solution vector
         ind3 = ind3max - 1
         do ind2 = ind2min, ind2max
             RHS_z(ind2,ind3) =  RHS_z(ind2,ind3) - array_az(ind2,ind3)* RHS_z(ind2,ind3+1)
            RHSu_z(ind2,ind3) = RHSu_z(ind2,ind3) - array_az(ind2,ind3)*RHSu_z(ind2,ind3+1)
            RHSv_z(ind2,ind3) = RHSv_z(ind2,ind3) - array_az(ind2,ind3)*RHSv_z(ind2,ind3+1)
            RHSw_z(ind2,ind3) = RHSw_z(ind2,ind3) - array_az(ind2,ind3)*RHSw_z(ind2,ind3+1)
            RHSe_z(ind2,ind3) = RHSe_z(ind2,ind3) - array_az(ind2,ind3)*RHSe_z(ind2,ind3+1)
         end do

         do ind3 = ind3max-2,ind3min,-1
            do ind2 = ind2min, ind2max
                RHS_z(ind2,ind3) =  RHS_z(ind2,ind3) - array_az(ind2,ind3)* RHS_z(ind2,ind3+1)  - array_bz(ind2,ind3)* RHS_z(ind2,ind3+2)
               RHSu_z(ind2,ind3) = RHSu_z(ind2,ind3) - array_az(ind2,ind3)*RHSu_z(ind2,ind3+1)  - array_bz(ind2,ind3)*RHSu_z(ind2,ind3+2)
               RHSv_z(ind2,ind3) = RHSv_z(ind2,ind3) - array_az(ind2,ind3)*RHSv_z(ind2,ind3+1)  - array_bz(ind2,ind3)*RHSv_z(ind2,ind3+2)
               RHSw_z(ind2,ind3) = RHSw_z(ind2,ind3) - array_az(ind2,ind3)*RHSw_z(ind2,ind3+1)  - array_bz(ind2,ind3)*RHSw_z(ind2,ind3+2)
               RHSe_z(ind2,ind3) = RHSe_z(ind2,ind3) - array_az(ind2,ind3)*RHSe_z(ind2,ind3+1)  - array_bz(ind2,ind3)*RHSe_z(ind2,ind3+2)
            end do
         end do

         ! Update with the solution
         do ind3 = ind3min,ind3max
            do ind2 = ind2min,ind2max
                Krho(ind2,ind1,ind3) =  RHS_z(ind2,ind3)
               Krhou(ind2,ind1,ind3) = RHSu_z(ind2,ind3)
               Krhov(ind2,ind1,ind3) = RHSv_z(ind2,ind3)
               Krhow(ind2,ind1,ind3) = RHSw_z(ind2,ind3)
               Krhoe(ind2,ind1,ind3) = RHSe_z(ind2,ind3)
            end do
         end do
      enddo

   end if

end subroutine irs4_ngh_c_v0

!===============================================================================
subroutine irs_pascal_c_v0
!===============================================================================
!> Application of the tridiagonal pascal parallelisation
!===============================================================================
   use mod_flow
   use mod_block
   use mod_constant
   use mod_time
   use mod_mpi_part
   use warnstop
   implicit none
   integer          :: ind1,ind2,ind3,irhs
   integer :: ind1min,ind1max,ind2min,ind2max,ind3min,ind3max
   real(wp),dimension(0:nx+1,0:ny+1,0:nz+1)  :: vc_ksi,vc_eta
   real(wp),dimension(0:nx+1,0:ny+1)         :: n2_ksi,n2_eta
   real(wp), dimension(1:nx+1,1:ny+1,1:nz+1) :: r_spec2
   real(wp), dimension(1:ny,1:nx,1:NRHS)     :: RHS_ksi
   real(wp) :: coef,rspecmh,rspecph

   !----------------------------------------
   ! Implicitation of direction ksi
   !----------------------------------------
   if (is_irs_i) then
!      ind3 : i;     ind1 : k;     ind2 : j
      ind1min = 1;   ind2min = 1;   ind3min = 1
      ind1max = nz;  ind2max = ny;  ind3max = nx

      ! If wall, do not take into account the wall point
      ! ================================================
      if (BC_face(2,1)%sort==0) ind2min = 2
      if (BC_face(2,2)%sort==0) ind2max = ny - 1
      if (BC_face(3,1)%sort==0) ind1min = 2
      if (BC_face(3,2)%sort==0) ind1max = nz - 1

      ! Calculation of contra-variant speed and norm of ksi
      ! ===================================================
      do ind2=ind2min,ind2max
         do ind3=ind3min-1,ind3max+1
            n2_ksi(ind3,ind2) = (y_eta(ind3,ind2)**2 + x_eta(ind3,ind2)**2) * ijacob_irs(ind3,ind2)**2
         end do
      end do

      do ind1=ind1min,ind1max
         do ind2=ind2min,ind2max
            do ind3=ind3min-1,ind3max+1
               vc_ksi(ind3,ind2,ind1) = (uu(ind3,ind2,ind1)*y_eta(ind3,ind2)-vv(ind3,ind2,ind1)*x_eta(ind3,ind2)) * ijacob_irs(ind3,ind2)
            end do
         end do
      end do

      ! Calculation of the spectral radius at ind-1/2
      ! =============================================
      do ind1=ind1min,ind1max
         do ind2=ind2min,ind2max
            do ind3=ind3min,ind3max+1
               rspecmh = sqrt(vc_ksi(ind3-1,ind2,ind1)**2) + c_(ind3-1,ind2,ind1)*sqrt(n2_ksi(ind3-1,ind2))
               rspecph = sqrt(vc_ksi(ind3,ind2,ind1)**2)   + c_(ind3,ind2,ind1)*sqrt(n2_ksi(ind3,ind2))
               r_spec2(ind3,ind2,ind1) = (0.5 * (rspecmh + rspecph))**2
            enddo
         enddo
      enddo


      do ind1 = ind1min, ind1max
         ! Preparation of the arrays and matrix
         ! ====================================
         coef = - theta_irs2 * deltat**2
         do irhs = 1,NRHS
            do ind2 = ind2min, ind2max
               do ind3 = ind3min,ind3max
                  ! Lower diagonal
                  DLp_ksi(ind2,ind3,irhs) = coef * r_spec2(ind3,ind2,ind1)
                  ! Upper diagonal
                  DUp_ksi(ind2,ind3,irhs) = coef * r_spec2(ind3+1,ind2,ind1)
                  ! Main diagonal
                  Dp_ksi(ind2,ind3,irhs)  = 1  - (DLp_ksi(ind2,ind3,irhs) + DUp_ksi(ind2,ind3,irhs))
               enddo
            end do
         end do

         do ind2 = ind2min, ind2max
            do ind3 = ind3min,ind3max
               RHS_ksi(ind2,ind3,1)   = Krho(ind3,ind2,ind1)
               RHS_ksi(ind2,ind3,2)   = Krhou(ind3,ind2,ind1)
               RHS_ksi(ind2,ind3,3)   = Krhov(ind3,ind2,ind1)
               RHS_ksi(ind2,ind3,4)   = Krhow(ind3,ind2,ind1)
               RHS_ksi(ind2,ind3,5)   = Krhoe(ind3,ind2,ind1)
            end do
         end do

         if (BC_face(1,1)%sort<=0) then
            coef = - theta_irs1 * deltat
            ind3 = ind3min
            do ind2 = ind2min, ind2max
               ! Lower diagonal
               DLp_ksi(ind2,ind3,irhs) = 0.
               ! Upper diagonal
               DUp_ksi(ind2,ind3,irhs) = coef  * (r_spec2(ind3+1,ind2,ind1))**0.5
               ! Main diagonal
               Dp_ksi(ind2,ind3,irhs)  = 1  - DUp_ksi(ind2,ind3,irhs)
            end do
         end if
         if (BC_face(1,2)%sort<=0) then
            coef = - theta_irs1 * deltat
            ind3 = ind3max
            do ind2 = ind2min, ind2max
               ! Lower diagonal
               DLp_ksi(ind2,ind3,irhs) = coef * (r_spec2(ind3,ind2,ind1))**0.5
               ! Upper diagonal
               DUp_ksi(ind2,ind3,irhs) = 0.
               ! Main diagonal
               Dp_ksi(ind2,ind3,irhs)  = 1 - DLp_ksi(ind2,ind3,irhs)
            end do
         end if


         ! Resolution of the tridiagonal linear system
         ! ===========================================
         if (bl(1)%BC(1)==1) then
            do irhs=1,NRHS
               call PaScaL_TDMA_many_solve_cycle(p_many_x,DLp_ksi(:,:,irhs),Dp_ksi(:,:,irhs),DUp_ksi(:,:,irhs),RHS_ksi(:,:,irhs),ny,nx)
            end do
         else
            do irhs=1,NRHS
               call PaScaL_TDMA_many_solve(p_many_x,DLp_ksi(:,:,irhs),Dp_ksi(:,:,irhs),DUp_ksi(:,:,irhs),RHS_ksi(:,:,irhs),ny,nx)
            end do
         end if

         do ind2 = ind2min, ind2max
            do ind3 = ind3min,ind3max
               Krho(ind3,ind2,ind1)  = RHS_ksi(ind2,ind3,1)
               Krhou(ind3,ind2,ind1) = RHS_ksi(ind2,ind3,2)
               Krhov(ind3,ind2,ind1) = RHS_ksi(ind2,ind3,3)
               Krhow(ind3,ind2,ind1) = RHS_ksi(ind2,ind3,4)
               Krhoe(ind3,ind2,ind1) = RHS_ksi(ind2,ind3,5)
            end do
         end do
      enddo

   end if

   !----------------------------------------
   ! Implicitation of direction eta
   !----------------------------------------
   if (is_irs_j) then
!      ind3 : j;     ind1 : k;     ind2 : i
      ind1min = 1;   ind2min = 1;   ind3min = 1
      ind1max = nz;  ind2max = nx;  ind3max = ny

      ! If wall, do not take into account the wall point
      ! ================================================
      if (BC_face(1,1)%sort==0) ind2min = 2
      if (BC_face(1,2)%sort==0) ind2max = nx - 1
      if (BC_face(3,1)%sort==0) ind1min = 2
      if (BC_face(3,2)%sort==0) ind1max = nz - 1

      ! Calculation of contra-variant speed and norm of eta
      ! ===================================================
      do ind3=ind3min-1,ind3max+1
         do ind2=ind2min,ind2max
            n2_eta(ind2,ind3) = (y_ksi(ind2,ind3)**2 + x_ksi(ind2,ind3)**2) * ijacob_irs(ind2,ind3)**2
         end do
      end do

      do ind1=ind1min,ind1max
         do ind3=ind3min-1,ind3max+1
            do ind2=ind2min,ind2max
               vc_eta(ind2,ind3,ind1) = (vv(ind2,ind3,ind1)*x_ksi(ind2,ind3) - uu(ind2,ind3,ind1)*y_ksi(ind2,ind3)) * ijacob_irs(ind2,ind3)
            end do
         end do
      end do

      ! Calculation of the spectral radius at ind-1/2
      ! =============================================
      do ind1=ind1min,ind1max
         do ind3=ind3min,ind3max+1
            do ind2=ind2min,ind2max
               rspecmh = sqrt(vc_eta(ind2,ind3-1,ind1)**2) + c_(ind2,ind3-1,ind1)*sqrt(n2_eta(ind2,ind3-1))
               rspecph = sqrt(vc_eta(ind2,ind3,ind1)**2)   + c_(ind2,ind3,ind1)*  sqrt(n2_eta(ind2,ind3))
               r_spec2(ind2,ind3,ind1) = (0.5 * (rspecmh + rspecph))**2
            enddo
         enddo
      enddo


      do ind1 = ind1min, ind1max
         ! Preparation of the arrays and matrix
         ! ====================================
         coef = - theta_irs2 * deltat**2
         do irhs = 1,NRHS
            do ind2 = ind2min, ind2max
               do ind3 = ind3min,ind3max
                  ! Lower diagonal
                  DLp_eta(ind2,ind3,irhs) = coef * r_spec2(ind2,ind3,ind1)
                  ! Upper diagonal
                  DUp_eta(ind2,ind3,irhs) = coef * r_spec2(ind2,ind3+1,ind1)
                  ! Main diagonal
                  Dp_eta(ind2,ind3,irhs)  = 1  - (DLp_eta(ind2,ind3,irhs) + DUp_eta(ind2,ind3,irhs))
               enddo
            end do
         end do

         if (BC_face(2,1)%sort<=0) then
            coef = - theta_irs1 * deltat
            ind3 = ind3min
            do irhs=1,NRHS
               do ind2 = ind2min, ind2max
                  ! Lower diagonal
                  DLp_eta(ind2,ind3,irhs) = 0.
                  ! Upper diagonal
                  DUp_eta(ind2,ind3,irhs) = coef * (r_spec2(ind2,ind3+1,ind1))**0.5
                  ! Main diagonal
                  Dp_eta(ind2,ind3,irhs)  = 1  - DUp_eta(ind2,ind3,irhs)
               end do
            end do
         end if
         if (BC_face(2,2)%sort<=0) then
            coef = - theta_irs1 * deltat
            ind3 = ind3max
            do irhs=1,NRHS
               do ind2 = ind2min, ind2max
                  ! Lower diagonal
                  DLp_eta(ind2,ind3,irhs) = coef * (r_spec2(ind2,ind3,ind1))**0.5
                  ! Upper diagonal
                  DUp_eta(ind2,ind3,irhs) = 0.
                  ! Main diagonal
                  Dp_eta(ind2,ind3,irhs)  = 1 - DLp_eta(ind2,ind3,irhs)
               end do
            end do
         end if

         ! Resolution of the tridiagonal linear system
         ! ===========================================
         if (bl(1)%BC(3)==1) then
            call PaScaL_TDMA_many_solve_cycle(p_many_y,DLp_eta(:,:,1),Dp_eta(:,:,1),DUp_eta(:,:,1), Krho(1:nx,1:ny,ind1),nx,ny)
            call PaScaL_TDMA_many_solve_cycle(p_many_y,DLp_eta(:,:,2),Dp_eta(:,:,2),DUp_eta(:,:,2),Krhou(1:nx,1:ny,ind1),nx,ny)
            call PaScaL_TDMA_many_solve_cycle(p_many_y,DLp_eta(:,:,3),Dp_eta(:,:,3),DUp_eta(:,:,3),Krhov(1:nx,1:ny,ind1),nx,ny)
            call PaScaL_TDMA_many_solve_cycle(p_many_y,DLp_eta(:,:,4),Dp_eta(:,:,4),DUp_eta(:,:,4),Krhow(1:nx,1:ny,ind1),nx,ny)
            call PaScaL_TDMA_many_solve_cycle(p_many_y,DLp_eta(:,:,5),Dp_eta(:,:,5),DUp_eta(:,:,5),Krhoe(1:nx,1:ny,ind1),nx,ny)
         else
            call PaScaL_TDMA_many_solve(p_many_y,DLp_eta(:,:,1),Dp_eta(:,:,1),DUp_eta(:,:,1), Krho(1:nx,1:ny,ind1),nx,ny)
            call PaScaL_TDMA_many_solve(p_many_y,DLp_eta(:,:,2),Dp_eta(:,:,2),DUp_eta(:,:,2),Krhou(1:nx,1:ny,ind1),nx,ny)
            call PaScaL_TDMA_many_solve(p_many_y,DLp_eta(:,:,3),Dp_eta(:,:,3),DUp_eta(:,:,3),Krhov(1:nx,1:ny,ind1),nx,ny)
            call PaScaL_TDMA_many_solve(p_many_y,DLp_eta(:,:,4),Dp_eta(:,:,4),DUp_eta(:,:,4),Krhow(1:nx,1:ny,ind1),nx,ny)
            call PaScaL_TDMA_many_solve(p_many_y,DLp_eta(:,:,5),Dp_eta(:,:,5),DUp_eta(:,:,5),Krhoe(1:nx,1:ny,ind1),nx,ny)
         end if
      enddo

   end if

   !----------------------------------------
   ! Implicitation of direction z
   !----------------------------------------
   if (is_irs_k .and. .not. is_2d) then
!      ind3 : k;     ind1 : j;     ind2 : i
      ind1min = 1;   ind2min = 1;   ind3min = 1
      ind1max = ny;  ind2max = nx;  ind3max = nz

      ! If wall, do not take into account the wall point
      ! ================================================
      if (BC_face(1,1)%sort==0) ind2min = 2
      if (BC_face(1,2)%sort==0) ind2max = nx - 1
      if (BC_face(2,1)%sort==0) ind1min = 2
      if (BC_face(2,2)%sort==0) ind1max = ny - 1

      ! Calculation of the spectral radius at ind-1/2
      ! =============================================
      do ind3=ind3min,ind3max+1
         do ind1=ind1min,ind1max
            do ind2=ind2min,ind2max
               rspecmh = sqrt(ww(ind2,ind1,ind3-1)**2) + c_(ind2,ind1,ind3-1)
               rspecph = sqrt(ww(ind2,ind1,ind3)**2)   + c_(ind2,ind1,ind3)
               r_spec2(ind2,ind1,ind3) = (0.5 * (rspecmh + rspecph))**2
            enddo
         enddo
      enddo

      do ind1 = ind1min, ind1max
         ! Preparation of the arrays and matrix
         ! ====================================
         coef = - theta_irs2 * deltat**2
         do irhs = 1,NRHS
            do ind3 = ind3min,ind3max
               do ind2 = ind2min,ind2max
                  ! Lower diagonal
                  DLp_z(ind2,ind3,irhs) = coef * r_spec2(ind2,ind1,ind3)
                  ! Upper diagonal
                  DUp_z(ind2,ind3,irhs) = coef * r_spec2(ind2,ind1,ind3+1)
                  ! Main diagonal
                  Dp_z(ind2,ind3,irhs)  = 1  - (DLp_z(ind2,ind3,irhs) + DUp_z(ind2,ind3,irhs))
               enddo
            end do
         end do

         if (BC_face(3,1)%sort<=0) then
            coef = - theta_irs1 * deltat
            ind3 = ind3min
            do irhs=1,NRHS
               do ind2 = ind2min, ind2max
                  ! Lower diagonal
                  DLp_z(ind2,ind3,irhs) = 0.
                  ! Upper diagonal
                  DUp_z(ind2,ind3,irhs) = coef * (r_spec2(ind2,ind1,ind3+1))**0.5
                  ! Main diagonal
                  Dp_z(ind2,ind3,irhs)  = 1  - DUp_z(ind2,ind3,irhs)
               end do
            end do
         end if
         if (BC_face(3,2)%sort<=0) then
            coef = - theta_irs1 * deltat
            ind3 = ind3max
            do irhs=1,NRHS
               do ind2 = ind2min, ind2max
                  ! Lower diagonal
                  DLp_z(ind2,ind3,irhs) = coef * (r_spec2(ind2,ind1,ind3))**0.5
                  ! Upper diagonal
                  DUp_z(ind2,ind3,irhs) = 0.
                  ! Main diagonal
                  Dp_z(ind2,ind3,irhs)  = 1 - DLp_z(ind2,ind3,irhs)
               end do
            end do
         end if

         ! Resolution of the tridiagonal linear system
         ! ===========================================
         if (bl(1)%BC(5)==1) then
            call PaScaL_TDMA_many_solve_cycle(p_many_z,DLp_z(:,:,1),Dp_z(:,:,1),DUp_z(:,:,1), Krho(1:nx,ind1,1:nz),nx,nz)
            call PaScaL_TDMA_many_solve_cycle(p_many_z,DLp_z(:,:,2),Dp_z(:,:,2),DUp_z(:,:,2),Krhou(1:nx,ind1,1:nz),nx,nz)
            call PaScaL_TDMA_many_solve_cycle(p_many_z,DLp_z(:,:,3),Dp_z(:,:,3),DUp_z(:,:,3),Krhov(1:nx,ind1,1:nz),nx,nz)
            call PaScaL_TDMA_many_solve_cycle(p_many_z,DLp_z(:,:,4),Dp_z(:,:,4),DUp_z(:,:,4),Krhow(1:nx,ind1,1:nz),nx,nz)
            call PaScaL_TDMA_many_solve_cycle(p_many_z,DLp_z(:,:,5),Dp_z(:,:,5),DUp_z(:,:,5),Krhoe(1:nx,ind1,1:nz),nx,nz)
         else
            call PaScaL_TDMA_many_solve(p_many_z,DLp_z(:,:,1),Dp_z(:,:,1),DUp_z(:,:,1), Krho(1:nx,ind1,1:nz),nx,nz)
            call PaScaL_TDMA_many_solve(p_many_z,DLp_z(:,:,2),Dp_z(:,:,2),DUp_z(:,:,2),Krhou(1:nx,ind1,1:nz),nx,nz)
            call PaScaL_TDMA_many_solve(p_many_z,DLp_z(:,:,3),Dp_z(:,:,3),DUp_z(:,:,3),Krhov(1:nx,ind1,1:nz),nx,nz)
            call PaScaL_TDMA_many_solve(p_many_z,DLp_z(:,:,4),Dp_z(:,:,4),DUp_z(:,:,4),Krhow(1:nx,ind1,1:nz),nx,nz)
            call PaScaL_TDMA_many_solve(p_many_z,DLp_z(:,:,5),Dp_z(:,:,5),DUp_z(:,:,5),Krhoe(1:nx,ind1,1:nz),nx,nz)
         end if
      enddo

   end if

end subroutine irs_pascal_c_v0
