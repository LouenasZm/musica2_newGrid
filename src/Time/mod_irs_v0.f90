!===============================================================================
subroutine init_irs_v0
!===============================================================================
  !> Initialisation of the Implicit Residual Smoothing
!===============================================================================
  use mod_constant
  use mod_time
  use mod_interface
  use mod_mpi_part
  use mod_grid
  use warnstop
  implicit none

  if (iproc==0) print *, 'Initialisation of IRS...'

  ! Safety checks
  ! =============
  if (is_2d) is_irs_k = .false.

  ! Parallelisation with pascal tridiagonal algo only for IRS2
  if (type_para==3 .and. iirs/=2) type_para = 2

  ! Parallelisation with ghost points only implemented for IRS2 or IRS4
  if (type_para==1) then
     if (iirs/=2 .and. iirs/=4 .and. iirs/=6) then
        type_para = 2
        if (iproc==0) print *,'Ghost points parallelisation only implemented for IRS2, IRS4 or IRS6',&
             'type_para changed to 2 : Scalapack band algo.'
     endif
  endif

  ! Parallelisation with ghost points
  ! =================================
  if (type_para==1) then
     if (iproc==0) print *,' Parallelisation of IRS with ghost points'

     if (iproc==0) print *,' Number of ghost points used for IRS : ', ngh

     ! Parallelisation with scalapack
     ! ==============================
  else if (type_para==2) then
     if (iproc==0) print *,'Parallelisation of IRS with SCALAPACK'
     call init_ictxt_v0

     ! Submatrix A(1:N, JA:JA+N-1): full matrix
     JA=1 ! ~> A(1:N,1:N)

     ! Parallelisation tridiag. for IRS2
     ! =================================
     if (iirs==2) then
        if (is_irs_i) then
           allocate(DL_x(nx),D_x(nx),DU_x(nx))

           ! Matrix descriptor
           ! =================
           DESCA_x(1) = 501
           DESCA_x(2) = ictxt(1)
           DESCA_x(3) = ngx
           DESCA_x(4) = nx
           DESCA_x(5) = 0
           DESCA_x(6) = nx

           ! RHS descriptor
           ! ==============
           DESCB_x(1) = 502
           DESCB_x(2) = ictxt(1)
           DESCB_x(3) = ngx
           DESCB_x(4) = nx
           DESCB_x(5) = 0
           DESCB_x(6) = nx

           ! Auxiliary array
           ! ===============
           LAF_x = 12*ndomx + 3*nx
           allocate(AF_x(LAF_x))
        endif

        if (is_irs_j) then
           allocate(DL_y(ny),D_y(ny),DU_y(ny))

           DESCA_y(1) = 501
           DESCA_y(2) = ictxt(2)
           DESCA_y(3) = ngy
           DESCA_y(4) = ny
           DESCA_y(5) = 0
           DESCA_y(6) = ny

           DESCB_y(1) = 502
           DESCB_y(2) = ictxt(2)
           DESCB_y(3) = ngy
           DESCB_y(4) = ny
           DESCB_y(5) = 0
           DESCB_y(6) = ny

           LAF_y = 12*ndomy + 3*ny
           allocate(AF_y(LAF_y))
        endif

        if (is_irs_k) then
           allocate(DL_z(nz),D_z(nz),DU_z(nz))

           DESCA_z(1) = 501
           DESCA_z(2) = ictxt(3)
           DESCA_z(3) = ngz
           DESCA_z(4) = nz
           DESCA_z(5) = 0
           DESCA_z(6) = nz

           DESCB_z(1) = 502
           DESCB_z(2) = ictxt(3)
           DESCB_z(3) = ngz
           DESCB_z(4) = nz
           DESCB_z(5) = 0
           DESCB_z(6) = nz

           LAF_z = 12*ndomz + 3*nz
           allocate(AF_z(LAF_z))
        endif

        ! Fill RHS
        ! ========
        ! number of RHS
        NRHS = 5
        ! Subvector of B: full matrix
        IB = 1 ! ~> B(1:N,NRHS)

        ! Work array
        ! ==========
        LWORK = 10*MAX(ndomx,ndomy,ndomz) + 4*NRHS
        allocate(WORK(LWORK))


        ! Parallelisation band for IRS4, 6, ...
        ! ====================================
     else
        if (iirs==4) bw = 2 ! matrix bandwidth
        if (iirs==6) bw = 3
        if (iirs==8) bw = 4


        if (is_irs_i) then
           allocate(A_x(2*bw+1,nx))
           A_x = 0.

           ! Matrix descriptor
           ! =================
           DESCA_x(1) = 1
           DESCA_x(2) = ictxt(1)
           DESCA_x(3) = ngx
           DESCA_x(4) = ngx
           DESCA_x(5) = nx
           DESCA_x(6) = nx
           DESCA_x(7) = 0
           DESCA_x(8) = 0
           DESCA_x(9) = 2*bw+1

           ! RHS descriptor
           ! ==============
           DESCB_x(1) = 502
           DESCB_x(2) = ictxt(1)
           DESCB_x(3) = ngx
           DESCB_x(4) = nx
           DESCB_x(5) = 0
           DESCB_x(6) = nx

           ! Auxiliary array
           ! ===============
           LAF_x = 2*nx*bw + 6*bw**2
           allocate(AF_x(LAF_x))
        endif

        if (is_irs_j) then
           allocate(A_y(2*bw+1,ny))
           A_y = 0.

           DESCA_y(1) = 1
           DESCA_y(2) = ictxt(2)
           DESCA_y(3) = ngy
           DESCA_y(4) = ngy
           DESCA_y(5) = ny
           DESCA_y(6) = ny
           DESCA_y(7) = 0
           DESCA_y(8) = 0
           DESCA_y(9) = 2*bw+1

           DESCB_y(1) = 502
           DESCB_y(2) = ictxt(2)
           DESCB_y(3) = ngy
           DESCB_y(4) = ny
           DESCB_y(5) = 0
           DESCB_y(6) = ny

           LAF_y = 2*ny*bw + 6*bw**2
           allocate(AF_y(LAF_y))
        endif

        if (is_irs_k) then
           allocate(A_z(2*bw+1,nz))
           A_z = 0.

           DESCA_z(1) = 1
           DESCA_z(2) = 0
           DESCA_z(3) = ngz
           DESCA_z(4) = ngz
           DESCA_z(5) = nz
           DESCA_z(6) = nz
           DESCA_z(7) = 0
           DESCA_z(8) = 0
           DESCA_z(9) = 2*bw+1

           DESCB_z(1) = 502
           DESCB_z(2) = 0
           DESCB_z(3) = ngz
           DESCB_z(4) = nz
           DESCB_z(5) = 0
           DESCB_z(6) = nz

           LAF_z = 2*nz*bw + 6*bw**2
           allocate(AF_z(LAF_z))
        endif

        ! number of RHS
        NRHS = 5
        ! Subvector of B: full matrix
        IB = 1 ! ~> B(1:N,NRHS)

        ! Work array
        ! ==========
        LWORK = bw*max(bw,NRHS)
        allocate(WORK(LWORK))

     endif


     ! Parallelisation with pascal for IRS2
     ! ====================================
  else if (type_para==3) then
     if (iproc==0) print *,'Parallelisation of IRS with Pascal'
     ! number of RHS
     NRHS = 5

     if (is_irs_i) allocate(DLp_x(1:ny,1:nx,NRHS),Dp_x(1:ny,1:nx,NRHS),DUp_x(1:ny,1:nx,NRHS))
     if (is_irs_j) allocate(DLp_y(1:nx,1:ny,NRHS),Dp_y(1:nx,1:ny,NRHS),DUp_y(1:nx,1:ny,NRHS))
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
     call mpistop('Definition of type_para for irs incorrect. Shutting down...',0)
  endif

  ! Creation of the arrays idx, idy, idz
  ! =======================================
  if (iirs==4) call init_grid_irs_v0
  if (iirs==2 .and. type_para==1) call init_grid_irs_v0
  if (iirs==6 .and. type_para==1) call init_grid_irs_v0

end subroutine init_irs_v0

!===============================================================================
subroutine init_ictxt_v0
!===============================================================================
  !> Initialisation of the contexts (BLACS) for each direction - Scalapack
!===============================================================================
   use mod_mpi_part
   use mod_constant
   use mod_time
   use mod_ngh
   implicit none
   integer, dimension(1,ndomx)   :: map_x
   integer, dimension(1,ndomy)   :: map_y
   integer, dimension(1,ndomz)   :: map_z
   integer  :: ind,incr

   ! Creation of Blacs context for direction x
   ! =========================================
   incr = 1
   do ind=0,nproc-1
      if (coord(2).eq.coord2(ind) .and. coord(3).eq.coord3(ind)) then
         map_x(1,incr) = ind
         incr = incr + 1
      endif
   enddo

   call BLACS_GET(0, 0, ictxt(1))
   call BLACS_GRIDMAP(ictxt(1),map_x,1,1,ndomx)
   DESCA_x(2) = ictxt(1)
   DESCB_x(2) = ictxt(1)

   ! Creation of Blacs context for direction y
   ! =========================================
   incr = 1
   do ind=0,nproc-1
      if (coord(1).eq.coord1(ind) .and. coord(3).eq.coord3(ind)) then
         map_y(1,incr) = ind
         incr = incr + 1
      endif
   enddo
   call BLACS_GET(0, 0, ictxt(2))
   call BLACS_GRIDMAP(ictxt(2),map_y,1,1,ndomy)
   DESCA_y(2) = ictxt(2)
   DESCB_y(2) = ictxt(2)

   ! Creation of Blacs context for direction z
   ! =========================================
   if (.not. is_2d) then
      incr = 1
      do ind=0,nproc-1
         if (coord(1).eq.coord1(ind) .and. coord(2).eq.coord2(ind)) then
            map_z(1,incr) = ind
            incr = incr + 1
         endif
      enddo
      call BLACS_GET(0, 0, ictxt(3))
      call BLACS_GRIDMAP(ictxt(3),map_z,1,1,ndomz)
      DESCA_z(2) = ictxt(3)
      DESCB_z(2) = ictxt(3)
   endif

end subroutine init_ictxt_v0

!===============================================================================
subroutine init_grid_irs_v0
!===============================================================================
  !> Initialisation of metrics arrays for IRS, based on idx-y-z
!===============================================================================
  use mod_bc
  use mod_grid
  use mod_constant
  use mod_time
  use mod_mpi
  implicit none

  ! Allocation for idx_irs, idy_irs and idz_irs arrays
  ! --------------------------------------------------
  if (.not.is_2d) allocate(idx_irs(1-ngh:nx+ngh),idy_irs(1-ngh:ny+ngh),idz_irs(1-ngh:nz+ngh))
  if (is_2d)      allocate(idx_irs(1-ngh:nx+ngh),idy_irs(1-ngh:ny+ngh))

  ! Attribution for interior points
  ! -------------------------------
  idx_irs(1:nx) = idx(1:nx)
  idy_irs(1:ny) = idy(1:ny)
  if (.not.is_2d) idz_irs(1:nz) = idz(1:nz)

  call grid_comm_irs_v0(idx_irs,nx,1)
  call grid_comm_irs_v0(idy_irs,ny,2)
  if (.not.is_2d) call grid_comm_irs_v0(idz_irs,nz,3)

  if (iirs==4)then
     idx_irs(ndx+1-ngh:nfx+ngh) = 0.5_wp * (idx_irs(ndx-ngh:nfx-1+ngh) + idx_irs(ndx+1-ngh:nfx+ngh))
     idy_irs(ndy+1-ngh:nfy+ngh) = 0.5_wp * (idy_irs(ndy-ngh:nfy-1+ngh) + idy_irs(ndy+1-ngh:nfy+ngh))
     if (.not. is_2d) idz_irs(ndz+1-ngh:nfz+ngh) = 0.5_wp * (idz_irs(ndz-ngh:nfz-1+ngh) + idz_irs(ndz+1-ngh:nfz+ngh))
  endif

end subroutine init_grid_irs_v0

!===============================================================================
subroutine grid_comm_irs_v0(dx_,nx_,dir)
!===============================================================================
  !> Communicate Cartesian metrics to fill ghost pts (for viscous metrics only)
!===============================================================================
  use mod_mpi_part
  use mod_grid
  use mod_constant
  use mod_time
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

end subroutine grid_comm_irs_v0

!!$!===============================================================================
!!$subroutine irs_routine_v0
!!$!===============================================================================
!!$  !> Application of the Implicit Residual Smoothing (IRS)
!!$!===============================================================================
!!$  use mod_interface
!!$  implicit none
!!$
!!$  call update_var_in_dw
!!$
!!$  call irs_solver
!!$
!!$  call update_var_of_dw
!!$
!!$end subroutine irs_routine_v0

!===============================================================================
subroutine irs2_ngh_v0
!===============================================================================
  !> Application of the ghost points parallelisation for IRS2
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
   real(wp),dimension(1:ny,ndx-ngh:nfx+ngh)   :: array_dx
   real(wp),dimension(1:ny,ndx-ngh:nfx+ngh-1) :: array_ax,array_cx
   real(wp),dimension(1:ny,ndx-ngh:nfx+ngh)   :: RHS_x,RHSu_x,RHSv_x,RHSw_x,RHSe_x
   real(wp),dimension(1:ny)                           :: temp_x
   real(wp),dimension(1:nx,ndy-ngh:nfy+ngh)   :: array_dy
   real(wp),dimension(1:nx,ndy-ngh:nfy+ngh-1) :: array_ay,array_cy
   real(wp),dimension(1:nx)                           :: temp_y
   real(wp),dimension(1:nx,ndz-ngh:nfz+ngh)   :: array_dz
   real(wp),dimension(1:nx,ndz-ngh:nfz+ngh-1) :: array_az,array_cz
   real(wp),dimension(1:nx,ndz-ngh:nfz+ngh)   :: RHS_z,RHSu_z,RHSv_z,RHSw_z,RHSe_z
   real(wp),dimension(1:nx)                           :: temp_z
   real(wp),dimension(ndx-ngh:nfx+ngh,ndy-ngh:nfy+ngh,ndz-ngh:nfz+ngh)  :: r_spec2
   real(wp) :: rspecmh,rspecph
   external :: ls_irs2_gh

   !----------------------------------------
   ! Implicitation of direction x
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

      ! Calculation of the spectral radius at ind-1/2
      ! =============================================
      if (.not. is_2d) then
         do ind1=ind1min,ind1max
            do ind2=ind2min,ind2max
               do ind3=ind3min+1,ind3max
                  rspecmh =  sqrt(uu(ind3-1,ind2,ind1)**2) + c_(ind3-1,ind2,ind1)
                  rspecph = sqrt(uu(ind3,ind2,ind1)**2) + c_(ind3,ind2,ind1)
                  r_spec2(ind3,ind2,ind1) = (0.5 * (rspecmh + rspecph))**2
               enddo
            enddo
         enddo
      else
         do ind2=ind2min,ind2max
            do ind3=ind3min+1,ind3max
               rspecmh = sqrt(uu(ind3-1,ind2,1)**2) + c_(ind3-1,ind2,1)
               rspecph = sqrt(uu(ind3,ind2,1)**2) + c_(ind3,ind2,1)
               r_spec2(ind3,ind2,1) = (0.5 * (rspecmh + rspecph))**2
            enddo
         enddo
      endif

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
            array_ax(ind2,ind3)  = - theta_irs1 * deltat * idx_irs(ind3)* (r_spec2(ind3+1,ind2,ind1))**0.5
            array_dx(ind2,ind3)  = 1. - array_ax(ind2,ind3)
         enddo

         ind3=ind3max
         do ind2 = ind2min, ind2max
            array_cx(ind2,ind3-1)  = - theta_irs1 * deltat * idx_irs(ind3) * (r_spec2(ind3,ind2,ind1))**0.5
            array_dx(ind2,ind3)    = 1. - array_cx(ind2,ind3-1)
         enddo

         do ind3 = ind3min+1,ind3max-1
            do ind2 = ind2min, ind2max
               temp_x(ind2)         = - theta_irs2 * (deltat * idx_irs(ind3))**2
               array_cx(ind2,ind3-1)= temp_x(ind2) * r_spec2(ind3,ind2,ind1)
               array_ax(ind2,ind3)  = temp_x(ind2) * r_spec2(ind3+1,ind2,ind1)
               array_dx(ind2,ind3)  = 1 - (array_cx(ind2,ind3-1) + array_ax(ind2,ind3))
            enddo
         enddo

         ! Filling of the RHS
         do ind3 = ind3min,ind3max
            do ind2 = ind2min,ind2max
                RHS_x(ind2,ind3) =  Krho(ind3,ind2,ind1)
               RHSu_x(ind2,ind3) = Krhou(ind3,ind2,ind1)
               RHSv_x(ind2,ind3) = Krhov(ind3,ind2,ind1)
               RHSw_x(ind2,ind3) = Krhow(ind3,ind2,ind1)
               RHSe_x(ind2,ind3) = Krhoe(ind3,ind2,ind1)
            enddo
         enddo

         ! Resolution of the tridiagonal linear system
         ! ===========================================
         ! Thomas' algorithm - Forward elimination phase
         do ind2 = ind2min, ind2max
            temp_x(ind2) = 1. / array_dx(ind2,ind3min)
            array_ax(ind2,ind3min) = array_ax(ind2,ind3min) * temp_x(ind2)
             RHS_x(ind2,ind3min) =  RHS_x(ind2,ind3min) * temp_x(ind2)
            RHSu_x(ind2,ind3min) = RHSu_x(ind2,ind3min) * temp_x(ind2)
            RHSv_x(ind2,ind3min) = RHSv_x(ind2,ind3min) * temp_x(ind2)
            RHSw_x(ind2,ind3min) = RHSw_x(ind2,ind3min) * temp_x(ind2)
            RHSe_x(ind2,ind3min) = RHSe_x(ind2,ind3min) * temp_x(ind2)
         enddo

         do ind3 = ind3min+1, ind3max-1
            do ind2 = ind2min, ind2max
               temp_x(ind2) = 1. / (array_dx(ind2,ind3) - array_cx(ind2,ind3-1)*array_ax(ind2,ind3-1))
               array_ax(ind2,ind3)   = array_ax(ind2,ind3) * temp_x(ind2)
                RHS_x(ind2,ind3) = ( RHS_x(ind2,ind3) - array_cx(ind2,ind3-1)* RHS_x(ind2,ind3-1)) * temp_x(ind2)
               RHSu_x(ind2,ind3) = (RHSu_x(ind2,ind3) - array_cx(ind2,ind3-1)*RHSu_x(ind2,ind3-1)) * temp_x(ind2)
               RHSv_x(ind2,ind3) = (RHSv_x(ind2,ind3) - array_cx(ind2,ind3-1)*RHSv_x(ind2,ind3-1)) * temp_x(ind2)
               RHSw_x(ind2,ind3) = (RHSw_x(ind2,ind3) - array_cx(ind2,ind3-1)*RHSw_x(ind2,ind3-1)) * temp_x(ind2)
               RHSe_x(ind2,ind3) = (RHSe_x(ind2,ind3) - array_cx(ind2,ind3-1)*RHSe_x(ind2,ind3-1)) * temp_x(ind2)
            enddo
         enddo

         ! Backward substitution phase
         do ind2 = ind2min, ind2max
            temp_x(ind2) = 1. / (array_dx(ind2,ind3max) - array_cx(ind2,ind3max-1)*array_ax(ind2,ind3max-1))
             RHS_x(ind2,ind3max) = ( RHS_x(ind2,ind3max) - array_cx(ind2,ind3max-1)* RHS_x(ind2,ind3max-1))*temp_x(ind2)
            RHSu_x(ind2,ind3max) = (RHSu_x(ind2,ind3max) - array_cx(ind2,ind3max-1)*RHSu_x(ind2,ind3max-1))*temp_x(ind2)
            RHSv_x(ind2,ind3max) = (RHSv_x(ind2,ind3max) - array_cx(ind2,ind3max-1)*RHSv_x(ind2,ind3max-1))*temp_x(ind2)
            RHSw_x(ind2,ind3max) = (RHSw_x(ind2,ind3max) - array_cx(ind2,ind3max-1)*RHSw_x(ind2,ind3max-1))*temp_x(ind2)
            RHSe_x(ind2,ind3max) = (RHSe_x(ind2,ind3max) - array_cx(ind2,ind3max-1)*RHSe_x(ind2,ind3max-1))*temp_x(ind2)
         enddo

         do ind3 = ind3max-1, ind3min, -1
            do ind2 = ind2min, ind2max
                RHS_x(ind2,ind3) =  RHS_x(ind2,ind3) - array_ax(ind2,ind3)* RHS_x(ind2,ind3+1)
               RHSu_x(ind2,ind3) = RHSu_x(ind2,ind3) - array_ax(ind2,ind3)*RHSu_x(ind2,ind3+1)
               RHSv_x(ind2,ind3) = RHSv_x(ind2,ind3) - array_ax(ind2,ind3)*RHSv_x(ind2,ind3+1)
               RHSw_x(ind2,ind3) = RHSw_x(ind2,ind3) - array_ax(ind2,ind3)*RHSw_x(ind2,ind3+1)
               RHSe_x(ind2,ind3) = RHSe_x(ind2,ind3) - array_ax(ind2,ind3)*RHSe_x(ind2,ind3+1)
            enddo
         enddo

         ! Update with the solution
         do ind3 = ind3min,ind3max
            do ind2 = ind2min,ind2max
                Krho(ind3,ind2,ind1) =  RHS_x(ind2,ind3)
               Krhou(ind3,ind2,ind1) = RHSu_x(ind2,ind3)
               Krhov(ind3,ind2,ind1) = RHSv_x(ind2,ind3)
               Krhow(ind3,ind2,ind1) = RHSw_x(ind2,ind3)
               Krhoe(ind3,ind2,ind1) = RHSe_x(ind2,ind3)
            enddo
         enddo
      enddo

   endif

   !----------------------------------------
   ! Implicitation of direction y
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

      ! Calculation of the spectral radius at ind-1/2
      ! =============================================
      if (.not. is_2d) then
         do ind1=ind1min,ind1max
            do ind3=ind3min+1,ind3max
               do ind2=ind2min,ind2max
                  rspecmh = sqrt(vv(ind2,ind3-1,ind1)**2) + c_(ind2,ind3-1,ind1)
                  rspecph = sqrt(vv(ind2,ind3,ind1)**2)   + c_(ind2,ind3,ind1)
                  r_spec2(ind2,ind3,ind1) = (0.5 * (rspecmh + rspecph))**2
               enddo
            enddo
         enddo
      else
         do ind3=ind3min+1,ind3max
            do ind2=ind2min,ind2max
               rspecmh = sqrt(vv(ind2,ind3-1,1)**2) + c_(ind2,ind3-1,1)
               rspecph = sqrt(vv(ind2,ind3,1)**2)   + c_(ind2,ind3,1)
               r_spec2(ind2,ind3,1) = (0.5 * (rspecmh + rspecph))**2
            enddo
         enddo
      endif

      do ind1 = ind1min, ind1max
         ind3=ind3min
         do ind2 = ind2min, ind2max
            array_ay(ind2,ind3)  = - theta_irs1 * deltat * idy_irs(ind3)* (r_spec2(ind2,ind3+1,ind1))**0.5
            array_dy(ind2,ind3)  = 1. - array_ay(ind2,ind3)
         enddo

         ind3=ind3max
         do ind2 = ind2min, ind2max
            array_cy(ind2,ind3-1)  = - theta_irs1 * deltat * idy_irs(ind3) * (r_spec2(ind2,ind3,ind1))**0.5
            array_dy(ind2,ind3)    = 1. - array_cy(ind2,ind3-1)
         enddo

         do ind3 = ind3min+1,ind3max-1
            do ind2 = ind2min, ind2max
               temp_y(ind2)         = - theta_irs2 * (deltat * idy_irs(ind3))**2
               array_cy(ind2,ind3-1)= temp_y(ind2) * r_spec2(ind2,ind3,ind1)
               array_ay(ind2,ind3)  = temp_y(ind2) * r_spec2(ind2,ind3+1,ind1)
               array_dy(ind2,ind3)  = 1 - (array_cy(ind2,ind3-1) + array_ay(ind2,ind3))
            enddo
         enddo

         ! Resolution of the tridiagonal linear system
         ! ===========================================
         ! Thomas' algorithm - Forward elimination phase
         do ind2 = ind2min, ind2max
            temp_y(ind2) = 1. / array_dy(ind2,ind3min)
            array_ay(ind2,ind3min) = array_ay(ind2,ind3min) * temp_y(ind2)
             Krho(ind2,ind3min,ind1) =  Krho(ind2,ind3min,ind1) * temp_y(ind2)
            Krhou(ind2,ind3min,ind1) = Krhou(ind2,ind3min,ind1) * temp_y(ind2)
            Krhov(ind2,ind3min,ind1) = Krhov(ind2,ind3min,ind1) * temp_y(ind2)
            Krhow(ind2,ind3min,ind1) = Krhow(ind2,ind3min,ind1) * temp_y(ind2)
            Krhoe(ind2,ind3min,ind1) = Krhoe(ind2,ind3min,ind1) * temp_y(ind2)
         enddo

         do ind3 = ind3min+1, ind3max-1
            do ind2 = ind2min, ind2max
               temp_y(ind2) = 1. / (array_dy(ind2,ind3) - array_cy(ind2,ind3-1)*array_ay(ind2,ind3-1))
               array_ay(ind2,ind3)   = array_ay(ind2,ind3) * temp_y(ind2)
                Krho(ind2,ind3,ind1) = ( Krho(ind2,ind3,ind1) - array_cy(ind2,ind3-1)* Krho(ind2,ind3-1,ind1)) * temp_y(ind2)
               Krhou(ind2,ind3,ind1) = (Krhou(ind2,ind3,ind1) - array_cy(ind2,ind3-1)*Krhou(ind2,ind3-1,ind1)) * temp_y(ind2)
               Krhov(ind2,ind3,ind1) = (Krhov(ind2,ind3,ind1) - array_cy(ind2,ind3-1)*Krhov(ind2,ind3-1,ind1)) * temp_y(ind2)
               Krhow(ind2,ind3,ind1) = (Krhow(ind2,ind3,ind1) - array_cy(ind2,ind3-1)*Krhow(ind2,ind3-1,ind1)) * temp_y(ind2)
               Krhoe(ind2,ind3,ind1) = (Krhoe(ind2,ind3,ind1) - array_cy(ind2,ind3-1)*Krhoe(ind2,ind3-1,ind1)) * temp_y(ind2)
            enddo
         enddo

         ! Backward substitution phase
         do ind2 = ind2min, ind2max
            temp_y(ind2) = 1. / (array_dy(ind2,ind3max) - array_cy(ind2,ind3max-1)*array_ay(ind2,ind3max-1))
             Krho(ind2,ind3max,ind1) = ( Krho(ind2,ind3max,ind1) - array_cy(ind2,ind3max-1)* Krho(ind2,ind3max-1,ind1))*temp_y(ind2)
            Krhou(ind2,ind3max,ind1) = (Krhou(ind2,ind3max,ind1) - array_cy(ind2,ind3max-1)*Krhou(ind2,ind3max-1,ind1))*temp_y(ind2)
            Krhov(ind2,ind3max,ind1) = (Krhov(ind2,ind3max,ind1) - array_cy(ind2,ind3max-1)*Krhov(ind2,ind3max-1,ind1))*temp_y(ind2)
            Krhow(ind2,ind3max,ind1) = (Krhow(ind2,ind3max,ind1) - array_cy(ind2,ind3max-1)*Krhow(ind2,ind3max-1,ind1))*temp_y(ind2)
            Krhoe(ind2,ind3max,ind1) = (Krhoe(ind2,ind3max,ind1) - array_cy(ind2,ind3max-1)*Krhoe(ind2,ind3max-1,ind1))*temp_y(ind2)
         enddo

         do ind3 = ind3max-1, ind3min, -1
            do ind2 = ind2min, ind2max
                Krho(ind2,ind3,ind1) =  Krho(ind2,ind3,ind1) - array_ay(ind2,ind3)* Krho(ind2,ind3+1,ind1)
               Krhou(ind2,ind3,ind1) = Krhou(ind2,ind3,ind1) - array_ay(ind2,ind3)*Krhou(ind2,ind3+1,ind1)
               Krhov(ind2,ind3,ind1) = Krhov(ind2,ind3,ind1) - array_ay(ind2,ind3)*Krhov(ind2,ind3+1,ind1)
               Krhow(ind2,ind3,ind1) = Krhow(ind2,ind3,ind1) - array_ay(ind2,ind3)*Krhow(ind2,ind3+1,ind1)
               Krhoe(ind2,ind3,ind1) = Krhoe(ind2,ind3,ind1) - array_ay(ind2,ind3)*Krhoe(ind2,ind3+1,ind1)
            enddo
         enddo
      enddo

   endif

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
         enddo

         ind3=ind3max
         do ind2 = ind2min, ind2max
            array_cz(ind2,ind3-1)  = - theta_irs1 * deltat * idz_irs(ind3) * (r_spec2(ind2,ind1,ind3))**0.5
            array_dz(ind2,ind3)    = 1. - array_cz(ind2,ind3-1)
         enddo

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
            enddo
         enddo

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
         enddo

         do ind3 = ind3min+1, ind3max-1
            do ind2 = ind2min, ind2max
               temp_z(ind2) = 1. / (array_dz(ind2,ind3) - array_cz(ind2,ind3-1)*array_az(ind2,ind3-1))
               array_az(ind2,ind3)   = array_az(ind2,ind3) * temp_z(ind2)
                RHS_z(ind2,ind3) = ( RHS_z(ind2,ind3) - array_cz(ind2,ind3-1)* RHS_z(ind2,ind3-1)) * temp_z(ind2)
               RHSu_z(ind2,ind3) = (RHSu_z(ind2,ind3) - array_cz(ind2,ind3-1)*RHSu_z(ind2,ind3-1)) * temp_z(ind2)
               RHSv_z(ind2,ind3) = (RHSv_z(ind2,ind3) - array_cz(ind2,ind3-1)*RHSv_z(ind2,ind3-1)) * temp_z(ind2)
               RHSw_z(ind2,ind3) = (RHSw_z(ind2,ind3) - array_cz(ind2,ind3-1)*RHSw_z(ind2,ind3-1)) * temp_z(ind2)
               RHSe_z(ind2,ind3) = (RHSe_z(ind2,ind3) - array_cz(ind2,ind3-1)*RHSe_z(ind2,ind3-1)) * temp_z(ind2)
            enddo
         enddo

         ! Backward substitution phase
         do ind2 = ind2min, ind2max
            temp_z(ind2) = 1. / (array_dz(ind2,ind3max) - array_cz(ind2,ind3max-1)*array_az(ind2,ind3max-1))
             RHS_z(ind2,ind3max) = ( RHS_z(ind2,ind3max) - array_cz(ind2,ind3max-1)* RHS_z(ind2,ind3max-1))*temp_z(ind2)
            RHSu_z(ind2,ind3max) = (RHSu_z(ind2,ind3max) - array_cz(ind2,ind3max-1)*RHSu_z(ind2,ind3max-1))*temp_z(ind2)
            RHSv_z(ind2,ind3max) = (RHSv_z(ind2,ind3max) - array_cz(ind2,ind3max-1)*RHSv_z(ind2,ind3max-1))*temp_z(ind2)
            RHSw_z(ind2,ind3max) = (RHSw_z(ind2,ind3max) - array_cz(ind2,ind3max-1)*RHSw_z(ind2,ind3max-1))*temp_z(ind2)
            RHSe_z(ind2,ind3max) = (RHSe_z(ind2,ind3max) - array_cz(ind2,ind3max-1)*RHSe_z(ind2,ind3max-1))*temp_z(ind2)
         enddo

         do ind3 = ind3max-1, ind3min, -1
            do ind2 = ind2min, ind2max
                RHS_z(ind2,ind3) =  RHS_z(ind2,ind3) - array_az(ind2,ind3)* RHS_z(ind2,ind3+1)
               RHSu_z(ind2,ind3) = RHSu_z(ind2,ind3) - array_az(ind2,ind3)*RHSu_z(ind2,ind3+1)
               RHSv_z(ind2,ind3) = RHSv_z(ind2,ind3) - array_az(ind2,ind3)*RHSv_z(ind2,ind3+1)
               RHSw_z(ind2,ind3) = RHSw_z(ind2,ind3) - array_az(ind2,ind3)*RHSw_z(ind2,ind3+1)
               RHSe_z(ind2,ind3) = RHSe_z(ind2,ind3) - array_az(ind2,ind3)*RHSe_z(ind2,ind3+1)
            enddo
         enddo

         ! Update with the solution
         do ind3 = ind3min,ind3max
            do ind2 = ind2min,ind2max
                Krho(ind2,ind1,ind3) =  RHS_z(ind2,ind3)
               Krhou(ind2,ind1,ind3) = RHSu_z(ind2,ind3)
               Krhov(ind2,ind1,ind3) = RHSv_z(ind2,ind3)
               Krhow(ind2,ind1,ind3) = RHSw_z(ind2,ind3)
               Krhoe(ind2,ind1,ind3) = RHSe_z(ind2,ind3)
            enddo
         enddo
      enddo

   endif
end subroutine irs2_ngh_v0

!===============================================================================
subroutine irs4_ngh_v0
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
   real(wp),dimension(1:ny,ndx-ngh:nfx+ngh)   :: array_dx
   real(wp),dimension(1:ny,ndx-ngh:nfx+ngh-1) :: array_ax,array_cx
   real(wp),dimension(1:ny,ndx-ngh:nfx+ngh-2) :: array_bx,array_ex
   real(wp),dimension(1:ny,ndx-ngh:nfx+ngh)   :: RHS_x,RHSu_x,RHSv_x,RHSw_x,RHSe_x
   real(wp),dimension(1:ny)                           :: temp_x1,temp_x2
   real(wp),dimension(1:nx,ndy-ngh:nfy+ngh)   :: array_dy
   real(wp),dimension(1:nx,ndy-ngh:nfy+ngh-1) :: array_ay,array_cy
   real(wp),dimension(1:nx,ndy-ngh:nfy+ngh-2) :: array_by,array_ey
   real(wp),dimension(1:nx)                           :: temp_y1,temp_y2
   real(wp),dimension(1:nx,ndz-ngh:nfz+ngh)   :: array_dz
   real(wp),dimension(1:nx,ndz-ngh:nfz+ngh-1) :: array_az,array_cz
   real(wp),dimension(1:nx,ndz-ngh:nfz+ngh-2) :: array_bz,array_ez
   real(wp),dimension(1:nx,ndz-ngh:nfz+ngh)   :: RHS_z,RHSu_z,RHSv_z,RHSw_z,RHSe_z
   real(wp),dimension(1:nx)                           :: temp_z1,temp_z2
   real(wp),dimension(ndx-ngh:nfx+ngh,ndy-ngh:nfy+ngh,ndz-ngh:nfz+ngh)  :: r_spec4
   real(wp) :: coef,coefmh,coefph,rspecmh,rspecph
   external :: ls_irs4_gh

   !----------------------------------------
   ! Implicitation of direction x
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

      ! Calculation of the spectral radius at ind-1/2
      ! =============================================
      if (.not. is_2d) then
         do ind1=ind1min,ind1max
            do ind2=ind2min,ind2max
               do ind3=ind3min+1,ind3max
                  rspecmh = sqrt(uu(ind3-1,ind2,ind1)**2)   + c_(ind3-1,ind2,ind1)
                  rspecph = sqrt(uu(ind3,ind2,ind1)**2)     + c_(ind3,ind2,ind1)
                  r_spec4(ind3,ind2,ind1) = (0.5 * (rspecmh + rspecph))**4
               enddo
            enddo
         enddo
      else
         do ind2=ind2min,ind2max
            do ind3=ind3min+1,ind3max
               rspecmh = sqrt(uu(ind3-1,ind2,1)**2)   + c_(ind3-1,ind2,1)
               rspecph = sqrt(uu(ind3,ind2,1)**2)     + c_(ind3,ind2,1)
               r_spec4(ind3,ind2,1) = (0.5 * (rspecmh + rspecph))**4
            enddo
         enddo
      endif

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
            coef     = theta_irs1 * deltat * idx_irs(ind3)
            array_bx(ind2,ind3)  = 0.
            array_ax(ind2,ind3)  = - coef * (r_spec4(ind3+1,ind2,ind1))**0.25
            array_dx(ind2,ind3)  = 1. + coef * (r_spec4(ind3+1,ind2,ind1))**0.25
         enddo

         ind3=ind3max
         do ind2 = ind2min, ind2max
            coef     = theta_irs1 * deltat * idx_irs(ind3)
            array_dx(ind2,ind3)    = 1. + coef * (r_spec4(ind3,ind2,ind1))**0.25
            array_cx(ind2,ind3-1)  = - coef * (r_spec4(ind3,ind2,ind1))**0.25
            array_ex(ind2,ind3-2)  = 0.
         enddo

         ind3=ind3min+1
         do ind2 = ind2min, ind2max
            coef     = - theta_irs2 * (deltat * idx_irs(ind3))**2
            array_bx(ind2,ind3)    = 0.
            array_ax(ind2,ind3)    = coef * (r_spec4(ind3+1,ind2,ind1))**0.5
            array_cx(ind2,ind3-1)  = coef * (r_spec4(ind3,ind2,ind1))**0.5
            array_dx(ind2,ind3)    = 1. - (array_ax(ind2,ind3) + array_cx(ind2,ind3-1))
         enddo

         ind3=ind3max-1
         do ind2 = ind2min, ind2max
            coef     = - theta_irs2 * (deltat * idx_irs(ind3))**2
            array_ax(ind2,ind3)    = coef * (r_spec4(ind3+1,ind2,ind1))**0.5
            array_cx(ind2,ind3-1)  = coef * (r_spec4(ind3,ind2,ind1))**0.5
            array_dx(ind2,ind3)    = 1. - (array_ax(ind2,ind3) + array_cx(ind2,ind3-1))
            array_ex(ind2,ind3-2)  = 0.
         enddo

         do ind3 = ind3min+2,ind3max-2
            do ind2 = ind2min, ind2max
               coefmh   = deltat * idx_irs(ind3)
               coefph   = deltat * idx_irs(ind3+1)

               coef     = theta_irs4 * deltat * 0.5_wp*(idx_irs(ind3)+idx_irs(ind3+1))

               array_ex(ind2,ind3-2)  = coef * coefmh**3 * r_spec4(ind3,ind2,ind1)
               array_bx(ind2,ind3)    = coef * coefph**3 * r_spec4(ind3+1,ind2,ind1)
               array_cx(ind2,ind3-1)  =  - 3.0_wp * array_ex(ind2,ind3-2) - array_bx(ind2,ind3)
               array_ax(ind2,ind3)    =  - 3.0_wp * array_bx(ind2,ind3) - array_ex(ind2,ind3-2)
               array_dx(ind2,ind3)    = 1.0_wp + 3.0_wp * (array_ex(ind2,ind3-2) + array_bx(ind2,ind3))
            enddo
         enddo

         ! Filling of the RHS
         do ind3 = ind3min,ind3max
            do ind2 = ind2min,ind2max
                RHS_x(ind2,ind3) =  Krho(ind3,ind2,ind1)
               RHSu_x(ind2,ind3) = Krhou(ind3,ind2,ind1)
               RHSv_x(ind2,ind3) = Krhov(ind3,ind2,ind1)
               RHSw_x(ind2,ind3) = Krhow(ind3,ind2,ind1)
               RHSe_x(ind2,ind3) = Krhoe(ind3,ind2,ind1)
            enddo
         enddo

         ! Resolution of the pentadiagonal linear system
         ! =============================================
         ind3 = ind3min
         do ind2 = ind2min, ind2max
            temp_x1(ind2) = 1 / array_dx(ind2,ind3)
            array_ax(ind2,ind3) = array_ax(ind2,ind3) * temp_x1(ind2)
            array_bx(ind2,ind3) = array_bx(ind2,ind3) * temp_x1(ind2)
             RHS_x(ind2,ind3)    =  RHS_x(ind2,ind3) * temp_x1(ind2)
            RHSu_x(ind2,ind3)    = RHSu_x(ind2,ind3) * temp_x1(ind2)
            RHSv_x(ind2,ind3)    = RHSv_x(ind2,ind3) * temp_x1(ind2)
            RHSw_x(ind2,ind3)    = RHSw_x(ind2,ind3) * temp_x1(ind2)
            RHSe_x(ind2,ind3)    = RHSe_x(ind2,ind3) * temp_x1(ind2)
         enddo

         ind3 = ind3min + 1
         do ind2 = ind2min, ind2max
            temp_x2(ind2) = array_cx(ind2,ind3-1)
            temp_x1(ind2) = 1. / (array_dx(ind2,ind3) - array_ax(ind2,ind3-1)*temp_x2(ind2))
            array_ax(ind2,ind3) = (array_ax(ind2,ind3) - array_bx(ind2,ind3-1)*temp_x2(ind2)) * temp_x1(ind2)
            array_bx(ind2,ind3) = array_bx(ind2,ind3) * temp_x1(ind2)
             RHS_x(ind2,ind3) = ( RHS_x(ind2,ind3) -  RHS_x(ind2,ind3-1)*temp_x2(ind2)) * temp_x1(ind2)
            RHSu_x(ind2,ind3) = (RHSu_x(ind2,ind3) - RHSu_x(ind2,ind3-1)*temp_x2(ind2)) * temp_x1(ind2)
            RHSv_x(ind2,ind3) = (RHSv_x(ind2,ind3) - RHSv_x(ind2,ind3-1)*temp_x2(ind2)) * temp_x1(ind2)
            RHSw_x(ind2,ind3) = (RHSw_x(ind2,ind3) - RHSw_x(ind2,ind3-1)*temp_x2(ind2)) * temp_x1(ind2)
            RHSe_x(ind2,ind3) = (RHSe_x(ind2,ind3) - RHSe_x(ind2,ind3-1)*temp_x2(ind2)) * temp_x1(ind2)
         enddo

         ! Step 5: for i = 3,...,n-2
         do ind3 = ind3min+2, ind3max-2
            do ind2 = ind2min, ind2max
               temp_x2(ind2) = array_cx(ind2,ind3-1) - array_ax(ind2,ind3-2)*array_ex(ind2,ind3-2)
               temp_x1(ind2) = 1. / (array_dx(ind2,ind3) - array_bx(ind2,ind3-2)*array_ex(ind2,ind3-2) - array_ax(ind2,ind3-1)*temp_x2(ind2))
               array_ax(ind2,ind3) = (array_ax(ind2,ind3) - array_bx(ind2,ind3-1)*temp_x2(ind2)) * temp_x1(ind2)
               array_bx(ind2,ind3) = array_bx(ind2,ind3) * temp_x1(ind2)
                RHS_x(ind2,ind3) = ( RHS_x(ind2,ind3) -  RHS_x(ind2,ind3-2)*array_ex(ind2,ind3-2) -  RHS_x(ind2,ind3-1)*temp_x2(ind2)) * temp_x1(ind2)
               RHSu_x(ind2,ind3) = (RHSu_x(ind2,ind3) - RHSu_x(ind2,ind3-2)*array_ex(ind2,ind3-2) - RHSu_x(ind2,ind3-1)*temp_x2(ind2)) * temp_x1(ind2)
               RHSv_x(ind2,ind3) = (RHSv_x(ind2,ind3) - RHSv_x(ind2,ind3-2)*array_ex(ind2,ind3-2) - RHSv_x(ind2,ind3-1)*temp_x2(ind2)) * temp_x1(ind2)
               RHSw_x(ind2,ind3) = (RHSw_x(ind2,ind3) - RHSw_x(ind2,ind3-2)*array_ex(ind2,ind3-2) - RHSw_x(ind2,ind3-1)*temp_x2(ind2)) * temp_x1(ind2)
               RHSe_x(ind2,ind3) = (RHSe_x(ind2,ind3) - RHSe_x(ind2,ind3-2)*array_ex(ind2,ind3-2) - RHSe_x(ind2,ind3-1)*temp_x2(ind2)) * temp_x1(ind2)
            enddo
         enddo

         ! Step 5: for i = n-1, n
         ind3 = ind3max - 1
         do ind2 = ind2min, ind2max
            temp_x2(ind2) = array_cx(ind2,ind3-1) - array_ax(ind2,ind3-2)*array_ex(ind2,ind3-2)
            temp_x1(ind2) = 1. / (array_dx(ind2,ind3) - array_bx(ind2,ind3-2)*array_ex(ind2,ind3-2) - array_ax(ind2,ind3-1)*temp_x2(ind2))
            array_ax(ind2,ind3) = (array_ax(ind2,ind3) - array_bx(ind2,ind3-1)*temp_x2(ind2)) * temp_x1(ind2)
             RHS_x(ind2,ind3) = ( RHS_x(ind2,ind3) -  RHS_x(ind2,ind3-2)*array_ex(ind2,ind3-2) -  RHS_x(ind2,ind3-1)*temp_x2(ind2)) * temp_x1(ind2)
            RHSu_x(ind2,ind3) = (RHSu_x(ind2,ind3) - RHSu_x(ind2,ind3-2)*array_ex(ind2,ind3-2) - RHSu_x(ind2,ind3-1)*temp_x2(ind2)) * temp_x1(ind2)
            RHSv_x(ind2,ind3) = (RHSv_x(ind2,ind3) - RHSv_x(ind2,ind3-2)*array_ex(ind2,ind3-2) - RHSv_x(ind2,ind3-1)*temp_x2(ind2)) * temp_x1(ind2)
            RHSw_x(ind2,ind3) = (RHSw_x(ind2,ind3) - RHSw_x(ind2,ind3-2)*array_ex(ind2,ind3-2) - RHSw_x(ind2,ind3-1)*temp_x2(ind2)) * temp_x1(ind2)
            RHSe_x(ind2,ind3) = (RHSe_x(ind2,ind3) - RHSe_x(ind2,ind3-2)*array_ex(ind2,ind3-2) - RHSe_x(ind2,ind3-1)*temp_x2(ind2)) * temp_x1(ind2)
         enddo

         ind3 = ind3max
         do ind2 = ind2min, ind2max
            temp_x2(ind2) = array_cx(ind2,ind3-1) - array_ax(ind2,ind3-2)*array_ex(ind2,ind3-2)
            temp_x1(ind2) = 1. / (array_dx(ind2,ind3) - array_bx(ind2,ind3-2)*array_ex(ind2,ind3-2) - array_ax(ind2,ind3-1)*temp_x2(ind2))
             RHS_x(ind2,ind3) = ( RHS_x(ind2,ind3) -  RHS_x(ind2,ind3-2)*array_ex(ind2,ind3-2) -  RHS_x(ind2,ind3-1)*temp_x2(ind2)) * temp_x1(ind2)
            RHSu_x(ind2,ind3) = (RHSu_x(ind2,ind3) - RHSu_x(ind2,ind3-2)*array_ex(ind2,ind3-2) - RHSu_x(ind2,ind3-1)*temp_x2(ind2)) * temp_x1(ind2)
            RHSv_x(ind2,ind3) = (RHSv_x(ind2,ind3) - RHSv_x(ind2,ind3-2)*array_ex(ind2,ind3-2) - RHSv_x(ind2,ind3-1)*temp_x2(ind2)) * temp_x1(ind2)
            RHSw_x(ind2,ind3) = (RHSw_x(ind2,ind3) - RHSw_x(ind2,ind3-2)*array_ex(ind2,ind3-2) - RHSw_x(ind2,ind3-1)*temp_x2(ind2)) * temp_x1(ind2)
            RHSe_x(ind2,ind3) = (RHSe_x(ind2,ind3) - RHSe_x(ind2,ind3-2)*array_ex(ind2,ind3-2) - RHSe_x(ind2,ind3-1)*temp_x2(ind2)) * temp_x1(ind2)
         enddo


         ! Step 6: computation of the solution vector
         ind3 = ind3max - 1
         do ind2 = ind2min, ind2max
             RHS_x(ind2,ind3) =  RHS_x(ind2,ind3) - array_ax(ind2,ind3)* RHS_x(ind2,ind3+1)
            RHSu_x(ind2,ind3) = RHSu_x(ind2,ind3) - array_ax(ind2,ind3)*RHSu_x(ind2,ind3+1)
            RHSv_x(ind2,ind3) = RHSv_x(ind2,ind3) - array_ax(ind2,ind3)*RHSv_x(ind2,ind3+1)
            RHSw_x(ind2,ind3) = RHSw_x(ind2,ind3) - array_ax(ind2,ind3)*RHSw_x(ind2,ind3+1)
            RHSe_x(ind2,ind3) = RHSe_x(ind2,ind3) - array_ax(ind2,ind3)*RHSe_x(ind2,ind3+1)
         enddo

         do ind3 = ind3max-2,ind3min,-1
            do ind2 = ind2min, ind2max
                RHS_x(ind2,ind3) =  RHS_x(ind2,ind3) - array_ax(ind2,ind3)* RHS_x(ind2,ind3+1)  - array_bx(ind2,ind3)* RHS_x(ind2,ind3+2)
               RHSu_x(ind2,ind3) = RHSu_x(ind2,ind3) - array_ax(ind2,ind3)*RHSu_x(ind2,ind3+1)  - array_bx(ind2,ind3)*RHSu_x(ind2,ind3+2)
               RHSv_x(ind2,ind3) = RHSv_x(ind2,ind3) - array_ax(ind2,ind3)*RHSv_x(ind2,ind3+1)  - array_bx(ind2,ind3)*RHSv_x(ind2,ind3+2)
               RHSw_x(ind2,ind3) = RHSw_x(ind2,ind3) - array_ax(ind2,ind3)*RHSw_x(ind2,ind3+1)  - array_bx(ind2,ind3)*RHSw_x(ind2,ind3+2)
               RHSe_x(ind2,ind3) = RHSe_x(ind2,ind3) - array_ax(ind2,ind3)*RHSe_x(ind2,ind3+1)  - array_bx(ind2,ind3)*RHSe_x(ind2,ind3+2)
            enddo
         enddo

         ! Update with the solution
         do ind3 = ind3min,ind3max
            do ind2 = ind2min,ind2max
                Krho(ind3,ind2,ind1) =  RHS_x(ind2,ind3)
               Krhou(ind3,ind2,ind1) = RHSu_x(ind2,ind3)
               Krhov(ind3,ind2,ind1) = RHSv_x(ind2,ind3)
               Krhow(ind3,ind2,ind1) = RHSw_x(ind2,ind3)
               Krhoe(ind3,ind2,ind1) = RHSe_x(ind2,ind3)
            enddo
         enddo
      enddo

   endif

   !----------------------------------------
   ! Implicitation of direction y
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

      ! Calculation of the spectral radius at ind-1/2
      ! =============================================
      if (.not. is_2d) then
         do ind1=ind1min,ind1max
            do ind3=ind3min+1,ind3max
               do ind2=ind2min,ind2max
                  rspecmh = sqrt(vv(ind2,ind3-1,ind1)**2)   + c_(ind2,ind3-1,ind1)
                  rspecph = sqrt(vv(ind2,ind3,ind1)**2)     + c_(ind2,ind3,ind1)
                  r_spec4(ind2,ind3,ind1) = (0.5 * (rspecmh + rspecph))**4
               enddo
            enddo
         enddo
      else
         do ind3=ind3min+1,ind3max
            do ind2=ind2min,ind2max
               rspecmh = sqrt(vv(ind2,ind3-1,1)**2)   + c_(ind2,ind3-1,1)
               rspecph = sqrt(vv(ind2,ind3,1)**2)     + c_(ind2,ind3,1)
               r_spec4(ind2,ind3,1) = (0.5 * (rspecmh + rspecph))**4
            enddo
         enddo
      endif

      do ind1 = ind1min, ind1max
         ! Construction of the pentadiagonal matrix
         ! ========================================
         ind3=ind3min
         do ind2 = ind2min, ind2max
            coef     = theta_irs1 * deltat * idy_irs(ind3)
            array_by(ind2,ind3)  = 0.
            array_ay(ind2,ind3)  = - coef * (r_spec4(ind2,ind3+1,ind1))**0.25
            array_dy(ind2,ind3)  = 1. + coef * (r_spec4(ind2,ind3+1,ind1))**0.25
         enddo

         ind3=ind3max
         do ind2 = ind2min, ind2max
            coef     = theta_irs1 * deltat * idy_irs(ind3)
            array_dy(ind2,ind3)    = 1. + coef * (r_spec4(ind2,ind3,ind1))**0.25
            array_cy(ind2,ind3-1)  = - coef * (r_spec4(ind2,ind3,ind1))**0.25
            array_ey(ind2,ind3-2)  = 0.
         enddo

         ind3=ind3min+1
         do ind2 = ind2min, ind2max
            coef     = - theta_irs2 * (deltat * idy_irs(ind3))**2
            array_by(ind2,ind3)    = 0.
            array_ay(ind2,ind3)    = coef * (r_spec4(ind2,ind3+1,ind1))**0.5
            array_cy(ind2,ind3-1)  = coef * (r_spec4(ind2,ind3,ind1))**0.5
            array_dy(ind2,ind3)    = 1. - (array_ay(ind2,ind3) + array_cy(ind2,ind3-1))
         enddo

         ind3=ind3max-1
         do ind2 = ind2min, ind2max
            coef     = - theta_irs2 * (deltat * idy_irs(ind3))**2
            array_ay(ind2,ind3)    = coef * (r_spec4(ind2,ind3+1,ind1))**0.5
            array_cy(ind2,ind3-1)  = coef * (r_spec4(ind2,ind3,ind1))**0.5
            array_dy(ind2,ind3)    = 1. - (array_ay(ind2,ind3) + array_cy(ind2,ind3-1))
            array_ey(ind2,ind3-2)  = 0.
         enddo

         do ind3 = ind3min+2,ind3max-2
            do ind2 = ind2min, ind2max
               coefmh   = deltat * idy_irs(ind3)
               coefph   = deltat * idy_irs(ind3+1)

               coef     = theta_irs4 * deltat * 0.5_wp*(idy_irs(ind3)+idy_irs(ind3+1))

               array_ey(ind2,ind3-2)  = coef * coefmh**3 * r_spec4(ind2,ind3,ind1)
               array_by(ind2,ind3)    = coef * coefph**3 * r_spec4(ind2,ind3+1,ind1)
               array_cy(ind2,ind3-1)  =  - 3.0_wp * array_ey(ind2,ind3-2) - array_by(ind2,ind3)
               array_ay(ind2,ind3)    =  - 3.0_wp * array_by(ind2,ind3) - array_ey(ind2,ind3-2)
               array_dy(ind2,ind3)    = 1.0_wp + 3.0_wp * (array_ey(ind2,ind3-2) + array_by(ind2,ind3))
            enddo
         enddo

         ! Resolution of the pentadiagonal linear system
         ind3 = ind3min
         do ind2 = ind2min, ind2max
            temp_y1(ind2) = 1 / array_dy(ind2,ind3)
            array_ay(ind2,ind3) = array_ay(ind2,ind3) * temp_y1(ind2)
            array_by(ind2,ind3) = array_by(ind2,ind3) * temp_y1(ind2)
             Krho(ind2,ind3,ind1)    =  Krho(ind2,ind3,ind1) * temp_y1(ind2)
            Krhou(ind2,ind3,ind1)    = Krhou(ind2,ind3,ind1) * temp_y1(ind2)
            Krhov(ind2,ind3,ind1)    = Krhov(ind2,ind3,ind1) * temp_y1(ind2)
            Krhow(ind2,ind3,ind1)    = Krhow(ind2,ind3,ind1) * temp_y1(ind2)
            Krhoe(ind2,ind3,ind1)    = Krhoe(ind2,ind3,ind1) * temp_y1(ind2)
         enddo

         ind3 = ind3min + 1
         do ind2 = ind2min, ind2max
            temp_y2(ind2) = array_cy(ind2,ind3-1)
            temp_y1(ind2) = 1. / (array_dy(ind2,ind3) - array_ay(ind2,ind3-1)*temp_y2(ind2))
            array_ay(ind2,ind3) = (array_ay(ind2,ind3) - array_by(ind2,ind3-1)*temp_y2(ind2)) * temp_y1(ind2)
            array_by(ind2,ind3) = array_by(ind2,ind3) * temp_y1(ind2)
             Krho(ind2,ind3,ind1) = ( Krho(ind2,ind3,ind1) -  Krho(ind2,ind3-1,ind1)*temp_y2(ind2)) * temp_y1(ind2)
            Krhou(ind2,ind3,ind1) = (Krhou(ind2,ind3,ind1) - Krhou(ind2,ind3-1,ind1)*temp_y2(ind2)) * temp_y1(ind2)
            Krhov(ind2,ind3,ind1) = (Krhov(ind2,ind3,ind1) - Krhov(ind2,ind3-1,ind1)*temp_y2(ind2)) * temp_y1(ind2)
            Krhow(ind2,ind3,ind1) = (Krhow(ind2,ind3,ind1) - Krhow(ind2,ind3-1,ind1)*temp_y2(ind2)) * temp_y1(ind2)
            Krhoe(ind2,ind3,ind1) = (Krhoe(ind2,ind3,ind1) - Krhoe(ind2,ind3-1,ind1)*temp_y2(ind2)) * temp_y1(ind2)
         enddo

         ! Step 5: for i = 3,...,n-2
         do ind3 = ind3min+2, ind3max-2
            do ind2 = ind2min, ind2max
               temp_y2(ind2) = array_cy(ind2,ind3-1) - array_ay(ind2,ind3-2)*array_ey(ind2,ind3-2)
               temp_y1(ind2) = 1. / (array_dy(ind2,ind3) - array_by(ind2,ind3-2)*array_ey(ind2,ind3-2) - array_ay(ind2,ind3-1)*temp_y2(ind2))
               array_ay(ind2,ind3) = (array_ay(ind2,ind3) - array_by(ind2,ind3-1)*temp_y2(ind2)) * temp_y1(ind2)
               array_by(ind2,ind3) = array_by(ind2,ind3) * temp_y1(ind2)
                Krho(ind2,ind3,ind1) = ( Krho(ind2,ind3,ind1) -  Krho(ind2,ind3-2,ind1)*array_ey(ind2,ind3-2) -  Krho(ind2,ind3-1,ind1)*temp_y2(ind2)) * temp_y1(ind2)
               Krhou(ind2,ind3,ind1) = (Krhou(ind2,ind3,ind1) - Krhou(ind2,ind3-2,ind1)*array_ey(ind2,ind3-2) - Krhou(ind2,ind3-1,ind1)*temp_y2(ind2)) * temp_y1(ind2)
               Krhov(ind2,ind3,ind1) = (Krhov(ind2,ind3,ind1) - Krhov(ind2,ind3-2,ind1)*array_ey(ind2,ind3-2) - Krhov(ind2,ind3-1,ind1)*temp_y2(ind2)) * temp_y1(ind2)
               Krhow(ind2,ind3,ind1) = (Krhow(ind2,ind3,ind1) - Krhow(ind2,ind3-2,ind1)*array_ey(ind2,ind3-2) - Krhow(ind2,ind3-1,ind1)*temp_y2(ind2)) * temp_y1(ind2)
               Krhoe(ind2,ind3,ind1) = (Krhoe(ind2,ind3,ind1) - Krhoe(ind2,ind3-2,ind1)*array_ey(ind2,ind3-2) - Krhoe(ind2,ind3-1,ind1)*temp_y2(ind2)) * temp_y1(ind2)
            enddo
         enddo

         ! Step 5: for i = n-1, n
         ind3 = ind3max - 1
         do ind2 = ind2min, ind2max
            temp_y2(ind2) = array_cy(ind2,ind3-1) - array_ay(ind2,ind3-2)*array_ey(ind2,ind3-2)
            temp_y1(ind2) = 1. / (array_dy(ind2,ind3) - array_by(ind2,ind3-2)*array_ey(ind2,ind3-2) - array_ay(ind2,ind3-1)*temp_y2(ind2))
            array_ay(ind2,ind3) = (array_ay(ind2,ind3) - array_by(ind2,ind3-1)*temp_y2(ind2)) * temp_y1(ind2)
             Krho(ind2,ind3,ind1) = ( Krho(ind2,ind3,ind1) -  Krho(ind2,ind3-2,ind1)*array_ey(ind2,ind3-2) -  Krho(ind2,ind3-1,ind1)*temp_y2(ind2)) * temp_y1(ind2)
            Krhou(ind2,ind3,ind1) = (Krhou(ind2,ind3,ind1) - Krhou(ind2,ind3-2,ind1)*array_ey(ind2,ind3-2) - Krhou(ind2,ind3-1,ind1)*temp_y2(ind2)) * temp_y1(ind2)
            Krhov(ind2,ind3,ind1) = (Krhov(ind2,ind3,ind1) - Krhov(ind2,ind3-2,ind1)*array_ey(ind2,ind3-2) - Krhov(ind2,ind3-1,ind1)*temp_y2(ind2)) * temp_y1(ind2)
            Krhow(ind2,ind3,ind1) = (Krhow(ind2,ind3,ind1) - Krhow(ind2,ind3-2,ind1)*array_ey(ind2,ind3-2) - Krhow(ind2,ind3-1,ind1)*temp_y2(ind2)) * temp_y1(ind2)
            Krhoe(ind2,ind3,ind1) = (Krhoe(ind2,ind3,ind1) - Krhoe(ind2,ind3-2,ind1)*array_ey(ind2,ind3-2) - Krhoe(ind2,ind3-1,ind1)*temp_y2(ind2)) * temp_y1(ind2)
         enddo

         ind3 = ind3max
         do ind2 = ind2min, ind2max
            temp_y2(ind2) = array_cy(ind2,ind3-1) - array_ay(ind2,ind3-2)*array_ey(ind2,ind3-2)
            temp_y1(ind2) = 1. / (array_dy(ind2,ind3) - array_by(ind2,ind3-2)*array_ey(ind2,ind3-2) - array_ay(ind2,ind3-1)*temp_y2(ind2))
             Krho(ind2,ind3,ind1) = ( Krho(ind2,ind3,ind1) -  Krho(ind2,ind3-2,ind1)*array_ey(ind2,ind3-2) -  Krho(ind2,ind3-1,ind1)*temp_y2(ind2)) * temp_y1(ind2)
            Krhou(ind2,ind3,ind1) = (Krhou(ind2,ind3,ind1) - Krhou(ind2,ind3-2,ind1)*array_ey(ind2,ind3-2) - Krhou(ind2,ind3-1,ind1)*temp_y2(ind2)) * temp_y1(ind2)
            Krhov(ind2,ind3,ind1) = (Krhov(ind2,ind3,ind1) - Krhov(ind2,ind3-2,ind1)*array_ey(ind2,ind3-2) - Krhov(ind2,ind3-1,ind1)*temp_y2(ind2)) * temp_y1(ind2)
            Krhow(ind2,ind3,ind1) = (Krhow(ind2,ind3,ind1) - Krhow(ind2,ind3-2,ind1)*array_ey(ind2,ind3-2) - Krhow(ind2,ind3-1,ind1)*temp_y2(ind2)) * temp_y1(ind2)
            Krhoe(ind2,ind3,ind1) = (Krhoe(ind2,ind3,ind1) - Krhoe(ind2,ind3-2,ind1)*array_ey(ind2,ind3-2) - Krhoe(ind2,ind3-1,ind1)*temp_y2(ind2)) * temp_y1(ind2)
         enddo


         ! Step 6: computation of the solution vector
         ind3 = ind3max - 1
         do ind2 = ind2min, ind2max
             Krho(ind2,ind3,ind1) =  Krho(ind2,ind3,ind1) - array_ay(ind2,ind3)* Krho(ind2,ind3+1,ind1)
            Krhou(ind2,ind3,ind1) = Krhou(ind2,ind3,ind1) - array_ay(ind2,ind3)*Krhou(ind2,ind3+1,ind1)
            Krhov(ind2,ind3,ind1) = Krhov(ind2,ind3,ind1) - array_ay(ind2,ind3)*Krhov(ind2,ind3+1,ind1)
            Krhow(ind2,ind3,ind1) = Krhow(ind2,ind3,ind1) - array_ay(ind2,ind3)*Krhow(ind2,ind3+1,ind1)
            Krhoe(ind2,ind3,ind1) = Krhoe(ind2,ind3,ind1) - array_ay(ind2,ind3)*Krhoe(ind2,ind3+1,ind1)
         enddo

         do ind3 = ind3max-2,ind3min,-1
            do ind2 = ind2min, ind2max
                Krho(ind2,ind3,ind1) =  Krho(ind2,ind3,ind1) - array_ay(ind2,ind3)* Krho(ind2,ind3+1,ind1)  - array_by(ind2,ind3)* Krho(ind2,ind3+2,ind1)
               Krhou(ind2,ind3,ind1) = Krhou(ind2,ind3,ind1) - array_ay(ind2,ind3)*Krhou(ind2,ind3+1,ind1)  - array_by(ind2,ind3)*Krhou(ind2,ind3+2,ind1)
               Krhov(ind2,ind3,ind1) = Krhov(ind2,ind3,ind1) - array_ay(ind2,ind3)*Krhov(ind2,ind3+1,ind1)  - array_by(ind2,ind3)*Krhov(ind2,ind3+2,ind1)
               Krhow(ind2,ind3,ind1) = Krhow(ind2,ind3,ind1) - array_ay(ind2,ind3)*Krhow(ind2,ind3+1,ind1)  - array_by(ind2,ind3)*Krhow(ind2,ind3+2,ind1)
               Krhoe(ind2,ind3,ind1) = Krhoe(ind2,ind3,ind1) - array_ay(ind2,ind3)*Krhoe(ind2,ind3+1,ind1)  - array_by(ind2,ind3)*Krhoe(ind2,ind3+2,ind1)
            enddo
         enddo
      enddo

   endif

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
         enddo

         ind3=ind3max
         do ind2 = ind2min, ind2max
            coef     = theta_irs1 * deltat * idz_irs(ind3)
            array_dz(ind2,ind3)    = 1. + coef * (r_spec4(ind2,ind1,ind3))**0.25
            array_cz(ind2,ind3-1)  = - coef * (r_spec4(ind2,ind1,ind3))**0.25
            array_ez(ind2,ind3-2)  = 0.
         enddo

         ind3=ind3min+1
         do ind2 = ind2min, ind2max
            coef     = - theta_irs2 * (deltat * idz_irs(ind3))**2
            array_bz(ind2,ind3)    = 0.
            array_az(ind2,ind3)    = coef * (r_spec4(ind2,ind1,ind3+1))**0.5
            array_cz(ind2,ind3-1)  = coef * (r_spec4(ind2,ind1,ind3))**0.5
            array_dz(ind2,ind3)    = 1. - (array_az(ind2,ind3) + array_cz(ind2,ind3-1))
         enddo

         ind3=ind3max-1
         do ind2 = ind2min, ind2max
            coef     = - theta_irs2 * (deltat * idz_irs(ind3))**2
            array_az(ind2,ind3)    = coef * (r_spec4(ind2,ind1,ind3+1))**0.5
            array_cz(ind2,ind3-1)  = coef * (r_spec4(ind2,ind1,ind3))**0.5
            array_dz(ind2,ind3)    = 1. - (array_az(ind2,ind3) + array_cz(ind2,ind3-1))
            array_ez(ind2,ind3-2)  = 0.
         enddo

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
         enddo


         ! Filling of the RHS
         do ind3 = ind3min,ind3max
            do ind2 = ind2min,ind2max
                RHS_z(ind2,ind3) =  Krho(ind2,ind1,ind3)
               RHSu_z(ind2,ind3) = Krhou(ind2,ind1,ind3)
               RHSv_z(ind2,ind3) = Krhov(ind2,ind1,ind3)
               RHSw_z(ind2,ind3) = Krhow(ind2,ind1,ind3)
               RHSe_z(ind2,ind3) = Krhoe(ind2,ind1,ind3)
            enddo
         enddo

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
         enddo

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
         enddo

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
            enddo
         enddo

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
         enddo

         ind3 = ind3max
         do ind2 = ind2min, ind2max
            temp_z2(ind2) = array_cz(ind2,ind3-1) - array_az(ind2,ind3-2)*array_ez(ind2,ind3-2)
            temp_z1(ind2) = 1. / (array_dz(ind2,ind3) - array_bz(ind2,ind3-2)*array_ez(ind2,ind3-2) - array_az(ind2,ind3-1)*temp_z2(ind2))
             RHS_z(ind2,ind3) = ( RHS_z(ind2,ind3) -  RHS_z(ind2,ind3-2)*array_ez(ind2,ind3-2) -  RHS_z(ind2,ind3-1)*temp_z2(ind2)) * temp_z1(ind2)
            RHSu_z(ind2,ind3) = (RHSu_z(ind2,ind3) - RHSu_z(ind2,ind3-2)*array_ez(ind2,ind3-2) - RHSu_z(ind2,ind3-1)*temp_z2(ind2)) * temp_z1(ind2)
            RHSv_z(ind2,ind3) = (RHSv_z(ind2,ind3) - RHSv_z(ind2,ind3-2)*array_ez(ind2,ind3-2) - RHSv_z(ind2,ind3-1)*temp_z2(ind2)) * temp_z1(ind2)
            RHSw_z(ind2,ind3) = (RHSw_z(ind2,ind3) - RHSw_z(ind2,ind3-2)*array_ez(ind2,ind3-2) - RHSw_z(ind2,ind3-1)*temp_z2(ind2)) * temp_z1(ind2)
            RHSe_z(ind2,ind3) = (RHSe_z(ind2,ind3) - RHSe_z(ind2,ind3-2)*array_ez(ind2,ind3-2) - RHSe_z(ind2,ind3-1)*temp_z2(ind2)) * temp_z1(ind2)
         enddo


         ! Step 6: computation of the solution vector
         ind3 = ind3max - 1
         do ind2 = ind2min, ind2max
             RHS_z(ind2,ind3) =  RHS_z(ind2,ind3) - array_az(ind2,ind3)* RHS_z(ind2,ind3+1)
            RHSu_z(ind2,ind3) = RHSu_z(ind2,ind3) - array_az(ind2,ind3)*RHSu_z(ind2,ind3+1)
            RHSv_z(ind2,ind3) = RHSv_z(ind2,ind3) - array_az(ind2,ind3)*RHSv_z(ind2,ind3+1)
            RHSw_z(ind2,ind3) = RHSw_z(ind2,ind3) - array_az(ind2,ind3)*RHSw_z(ind2,ind3+1)
            RHSe_z(ind2,ind3) = RHSe_z(ind2,ind3) - array_az(ind2,ind3)*RHSe_z(ind2,ind3+1)
         enddo

         do ind3 = ind3max-2,ind3min,-1
            do ind2 = ind2min, ind2max
                RHS_z(ind2,ind3) =  RHS_z(ind2,ind3) - array_az(ind2,ind3)* RHS_z(ind2,ind3+1)  - array_bz(ind2,ind3)* RHS_z(ind2,ind3+2)
               RHSu_z(ind2,ind3) = RHSu_z(ind2,ind3) - array_az(ind2,ind3)*RHSu_z(ind2,ind3+1)  - array_bz(ind2,ind3)*RHSu_z(ind2,ind3+2)
               RHSv_z(ind2,ind3) = RHSv_z(ind2,ind3) - array_az(ind2,ind3)*RHSv_z(ind2,ind3+1)  - array_bz(ind2,ind3)*RHSv_z(ind2,ind3+2)
               RHSw_z(ind2,ind3) = RHSw_z(ind2,ind3) - array_az(ind2,ind3)*RHSw_z(ind2,ind3+1)  - array_bz(ind2,ind3)*RHSw_z(ind2,ind3+2)
               RHSe_z(ind2,ind3) = RHSe_z(ind2,ind3) - array_az(ind2,ind3)*RHSe_z(ind2,ind3+1)  - array_bz(ind2,ind3)*RHSe_z(ind2,ind3+2)
            enddo
         enddo

         ! Update with the solution
         do ind3 = ind3min,ind3max
            do ind2 = ind2min,ind2max
                Krho(ind2,ind1,ind3) =  RHS_z(ind2,ind3)
               Krhou(ind2,ind1,ind3) = RHSu_z(ind2,ind3)
               Krhov(ind2,ind1,ind3) = RHSv_z(ind2,ind3)
               Krhow(ind2,ind1,ind3) = RHSw_z(ind2,ind3)
               Krhoe(ind2,ind1,ind3) = RHSe_z(ind2,ind3)
            enddo
         enddo
      enddo

   endif
   
end subroutine irs4_ngh_v0

!===============================================================================
subroutine irs6_ngh_v0
!===============================================================================
!> Application of the ghost points parallelisation for IRS6
!===============================================================================
   use mod_flow
   use mod_constant
   use mod_time
   use mod_mpi
   use mod_interface
   use warnstop
   implicit none
   integer  :: ind1,ind2,ind3
   integer  :: ind1min,ind1max,ind2min,ind2max,ind3min,ind3max
   real(wp),dimension(ndx-ngh:nfx+ngh,ndy-ngh:nfy+ngh,ndz-ngh:nfz+ngh) :: r_spec6
   real(wp),dimension(ndx-ngh:nfx+ngh)     :: array_dx
   real(wp),dimension(ndx-ngh:nfx+ngh-1)   :: array_ax,array_ex
   real(wp),dimension(ndx-ngh:nfx+ngh-2)   :: array_bx,array_fx
   real(wp),dimension(ndx-ngh:nfx+ngh-3)   :: array_cx,array_gx
   real(wp),dimension(ndy-ngh:nfy+ngh)     :: array_dy
   real(wp),dimension(ndy-ngh:nfy+ngh-1)   :: array_ay,array_ey
   real(wp),dimension(ndy-ngh:nfy+ngh-2)   :: array_by,array_fy
   real(wp),dimension(ndy-ngh:nfy+ngh-3)   :: array_cy,array_gy
   real(wp),dimension(ndz-ngh:nfz+ngh)     :: array_dz
   real(wp),dimension(ndz-ngh:nfz+ngh-1)   :: array_az,array_ez
   real(wp),dimension(ndz-ngh:nfz+ngh-2)   :: array_bz,array_fz
   real(wp),dimension(ndz-ngh:nfz+ngh-3)   :: array_cz,array_gz
   real(wp) :: coef,coefmh,coefph,rspecmh,rspecph,unsur6

   unsur6 = 1./6.

   !----------------------------------------
   ! Implicitation of direction x
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

      ! Calculation of the spectral radius at ind-1/2
      ! =============================================
      if (.not. is_2d) then
         do ind1=ind1min,ind1max
            do ind2=ind2min,ind2max
               do ind3=ind3min+1,ind3max
                  rspecmh = sqrt(uu(ind3-1,ind2,ind1)**2)   + c_(ind3-1,ind2,ind1)
                  rspecph = sqrt(uu(ind3,ind2,ind1)**2)     + c_(ind3,ind2,ind1)
                  r_spec6(ind3,ind2,ind1) = (0.5 * (rspecmh + rspecph))**6
               enddo
            enddo
         enddo
      else
         do ind2=ind2min,ind2max
            do ind3=ind3min+1,ind3max
               rspecmh = sqrt(uu(ind3-1,ind2,1)**2)   + c_(ind3-1,ind2,1)
               rspecph = sqrt(uu(ind3,ind2,1)**2)     + c_(ind3,ind2,1)
               r_spec6(ind3,ind2,1) = (0.5 * (rspecmh + rspecph))**6
            enddo
         enddo
      endif

      do ind1 = ind1min, ind1max
         do ind2 = ind2min, ind2max
            ! Preparation of the arrays and matrix
            ! ====================================
            ! The linear systems is Mx = y with :
            !     _M = (m_ij) for i,j in {1,..,n} and with :
            !                 m_ii = d_i for i in {1,..,n}
            !                 m_ii+1 = a_i for i in {1,..,n-1}
            !                 m_ii+2 = b_i for i in {1,..,n-2}
            !                 m_ii+3 = c_i for i in {1,..,n-3}
            !                 m_ii-1 = e_i-1 for i in {2,..,n}
            !                 m_ii-2 = f_i-2 for i in {3,..,n}
            !                 m_ii-3 = g_i-3 for i in {4,..,n}
            do ind3 = ind3min+3,ind3max-3
               coef     = theta_irs6 * (deltat * idx_irs(ind3))**6
               coefmh   = coef * r_spec6(ind3,ind2,ind1)
               coefph   = coef * r_spec6(ind3+1,ind2,ind1)

               ! Third lower diagonal
               array_gx(ind3-3)  = -coefmh
               ! Third upper diagonal
               array_cx(ind3)    = -coefph
               ! Second lower diagonal
               array_fx(ind3-2)  = 5.0_wp*coefmh + coefph
               ! Second upper diagonal
               array_bx(ind3)    = coefmh + 5.0_wp*coefph
               ! First lower diagonal
               array_ex(ind3-1)  = -10.0_wp*coefmh - 5.0_wp*coefph
               ! First upper diagonal
               array_ax(ind3)    = -5.0_wp*coefmh - 10.0_wp*coefph
               ! Main diagonal
               array_dx(ind3)    = 10.0_wp*(coefmh + coefph) + 1.0_wp
            enddo

!            ind3 = ind3min
!            array_cx(ind3) = 0.0_wp
!            array_bx(ind3)  = 0.0_wp
!            array_ax(ind3)  = 0.0_wp
!            array_dx(ind3)  = 1.0_wp
!
!            ind3 = ind3max
!            array_dx(ind3)    = 1.0_wp
!            array_ex(ind3-1)  = 0.0_wp
!            array_fx(ind3-2)  = 0.0_wp
!            array_gx(ind3-3)  = 0.0_wp
!
!            ind3 = ind3min + 1
!            array_cx(ind3) = 0.0_wp
!            array_bx(ind3)  = 0.0_wp
!            array_ax(ind3)    = 0.0_wp
!            array_ex(ind3-1)  = 0.0_wp
!            array_dx(ind3)    = 1.0_wp
!
!            ind3 = ind3max - 1
!            array_ax(ind3)    = 0.0_wp
!            array_ex(ind3-1)  = 0.0_wp
!            array_dx(ind3)    = 1.0_wp
!            array_fx(ind3-2)  = 0.0_wp
!            array_gx(ind3-3)  = 0.0_wp
!
!            ind3 = ind3min + 2
!            array_cx(ind3)    = 0.0_wp
!            array_fx(ind3-2)  = 0.0_wp
!            array_bx(ind3)    = 0.0_wp
!            array_ex(ind3-1)  = 0.0_wp
!            array_ax(ind3)    = 0.0_wp
!            array_dx(ind3)    = 1.0_wp
!
!            ind3 = ind3max - 2
!            array_gx(ind3-3)  = 0.0_wp
!            array_fx(ind3-2)  = 0.0_wp
!            array_bx(ind3)    = 0.0_wp
!            array_ex(ind3-1)  = 0.0_wp
!            array_ax(ind3)    = 0.0_wp
!            array_dx(ind3)    = 1.0_wp

            ! IRS1
            ind3 = ind3min
            coef     = theta_irs1 * deltat * idx_irs(ind3)
            array_cx(ind3) = 0.0_wp
            array_bx(ind3)  = 0.0_wp
            array_ax(ind3)  = - coef * (r_spec6(ind3+1,ind2,ind1))**(unsur6)
            array_dx(ind3)  = 1.0_wp + coef * (r_spec6(ind3+1,ind2,ind1))**(unsur6)

            ind3 = ind3max
            coef     = theta_irs1 * deltat * idx_irs(ind3)
            array_dx(ind3)    = 1.0_wp + coef * (r_spec6(ind3,ind2,ind1))**(unsur6)
            array_ex(ind3-1)  = - coef * (r_spec6(ind3,ind2,ind1))**(unsur6)
            array_fx(ind3-2)  = 0.0_wp
            array_gx(ind3-3)  = 0.0_wp

            ! IRS2
            ind3 = ind3min + 1
            coef     = - theta_irs2 * (deltat * idx_irs(ind3))**2
            array_cx(ind3) = 0.0_wp
            array_bx(ind3)  = 0.0_wp
            array_ax(ind3)    = coef * (r_spec6(ind3+1,ind2,ind1))**(unsur6*2)
            array_ex(ind3-1)  = coef * (r_spec6(ind3,ind2,ind1))**(unsur6*2)
            array_dx(ind3)    = 1.0_wp - (array_ax(ind3) + array_ex(ind3-1))

            ind3 = ind3max - 1
            coef     = - theta_irs2 * (deltat * idx_irs(ind3))**2
            array_ax(ind3)    = coef * (r_spec6(ind3+1,ind2,ind1))**(unsur6*2)
            array_ex(ind3-1)  = coef * (r_spec6(ind3,ind2,ind1))**(unsur6*2)
            array_dx(ind3)    = 1.0_wp - (array_ax(ind3) + array_ex(ind3-1))
            array_fx(ind3-2)  = 0.0_wp
            array_gx(ind3-3)  = 0.0_wp


            ! IRS4
            ind3 = ind3min + 2
            coef     = theta_irs4 * (deltat * idx_irs(ind3))**4
            array_cx(ind3)    = 0.0_wp
            array_fx(ind3-2)  = coef * r_spec6(ind3,ind2,ind1)**(4*unsur6)
            array_bx(ind3)    = coef * r_spec6(ind3+1,ind2,ind1)**(4*unsur6)
            array_ex(ind3-1)  =  - 3.0_wp * array_fx(ind3-2) - array_bx(ind3)
            array_ax(ind3)    =  - 3.0_wp * array_bx(ind3) - array_fx(ind3-2)
            array_dx(ind3)    = 1.0_wp + 3.0_wp * (array_fx(ind3-2) + array_bx(ind3))

            ind3 = ind3max - 2
            coef     = theta_irs4 * (deltat * idx_irs(ind3))**4
            array_gx(ind3-3)  = 0.0_wp
            array_fx(ind3-2)  = coef * r_spec6(ind3,ind2,ind1)**(4*unsur6)
            array_bx(ind3)    = coef * r_spec6(ind3+1,ind2,ind1)**(4*unsur6)
            array_ex(ind3-1)  =  - 3.0_wp * array_fx(ind3-2) - array_bx(ind3)
            array_ax(ind3)    =  - 3.0_wp * array_bx(ind3) - array_fx(ind3-2)
            array_dx(ind3)    = 1.0_wp + 3.0_wp * (array_fx(ind3-2) + array_bx(ind3))

            ! Resolution of the heptadiagonal linear system
            ! =============================================
            call ls_irs6_gh(array_ax,array_bx,array_cx,array_dx,array_ex,array_fx,array_gx,Krho(ind3min:ind3max,ind2,ind1),Krhou(ind3min:ind3max,ind2,ind1),&
                                          Krhov(ind3min:ind3max,ind2,ind1),Krhow(ind3min:ind3max,ind2,ind1),Krhoe(ind3min:ind3max,ind2,ind1),ind3min,ind3max)

         enddo
      enddo
   endif

   !----------------------------------------
   ! Implicitation of direction y
   !----------------------------------------
   if (is_irs_j) then
!      ind3 : j;     ind1 : k;     ind2 : i
      ind1min = 1;   ind2min = 1;   ind3min = ndy - ngh
      ind1max = nz;  ind2max = nx;  ind3max = nfy + ngh

      ! Communication of ghost points
      ! =============================
      ! call communication_NS(Krho,Krhou,Krhov,Krhow,Krhoe)
      call communication_(Krho,Krhou,Krhov,Krhow,Krhoe)

      ! If wall, do not take into account the wall point
      ! ================================================
      if (BC_face(1,1)%sort==0) ind2min = 2
      if (BC_face(1,2)%sort==0) ind2max = nx - 1
      if (BC_face(3,1)%sort==0) ind1min = 2
      if (BC_face(3,2)%sort==0) ind1max = nz - 1

      ! Calculation of the spectral radius at ind-1/2
      ! =============================================
      if (.not. is_2d) then
         do ind1=ind1min,ind1max
            do ind3=ind3min+1,ind3max
               do ind2=ind2min,ind2max
                  rspecmh = sqrt(vv(ind2,ind3-1,ind1)**2)   + c_(ind2,ind3-1,ind1)
                  rspecph = sqrt(vv(ind2,ind3,ind1)**2)     + c_(ind2,ind3,ind1)
                  r_spec6(ind2,ind3,ind1) = (0.5 * (rspecmh + rspecph))**6
               enddo
            enddo
         enddo
      else
         do ind3=ind3min+1,ind3max
            do ind2=ind2min,ind2max
               rspecmh = sqrt(vv(ind2,ind3-1,1)**2)   + c_(ind2,ind3-1,1)
               rspecph = sqrt(vv(ind2,ind3,1)**2)     + c_(ind2,ind3,1)
               r_spec6(ind2,ind3,1) = (0.5 * (rspecmh + rspecph))**6
            enddo
         enddo
      endif

      do ind1 = ind1min, ind1max
         do ind2 = ind2min, ind2max
            ! Preparation of the arrays and matrix
            ! ====================================
            do ind3 = ind3min+3,ind3max-3
               coef     = theta_irs6 * (deltat * idy_irs(ind3))**6
               coefmh   = coef * r_spec6(ind2,ind3,ind1)
               coefph   = coef * r_spec6(ind2,ind3+1,ind1)

               ! Third lower diagonal
               array_gy(ind3-3)  = -coefmh
               ! Third upper diagonal
               array_cy(ind3)    = -coefph
               ! Second lower diagonal
               array_fy(ind3-2)  = 5.0_wp*coefmh + coefph
               ! Second upper diagonal
               array_by(ind3)    = coefmh + 5.0_wp*coefph
               ! First lower diagonal
               array_ey(ind3-1)  = -10.0_wp*coefmh - 5.0_wp*coefph
               ! First upper diagonal
               array_ay(ind3)    = -5.0_wp*coefmh - 10.0_wp*coefph
               ! Main diagonal
               array_dy(ind3)    = 10.0_wp*(coefmh + coefph) + 1.0_wp
            enddo

!            ind3 = ind3min
!            array_cy(ind3) = 0.0_wp
!            array_by(ind3)  = 0.0_wp
!            array_ay(ind3)  = 0.0_wp
!            array_dy(ind3)  = 1.0_wp
!
!            ind3 = ind3max
!            array_dy(ind3)    = 1.0_wp
!            array_ey(ind3-1)  = 0.0_wp
!            array_fy(ind3-2)  = 0.0_wp
!            array_gy(ind3-3)  = 0.0_wp
!
!            ind3 = ind3min + 1
!            array_cy(ind3) = 0.0_wp
!            array_by(ind3)  = 0.0_wp
!            array_ay(ind3)    = 0.0_wp
!            array_ey(ind3-1)  = 0.0_wp
!            array_dy(ind3)    = 1.0_wp
!
!            ind3 = ind3max - 1
!            array_ay(ind3)    = 0.0_wp
!            array_ey(ind3-1)  = 0.0_wp
!            array_dy(ind3)    = 1.0_wp
!            array_fy(ind3-2)  = 0.0_wp
!            array_gy(ind3-3)  = 0.0_wp
!
!            ind3 = ind3min + 2
!            array_cy(ind3)    = 0.0_wp
!            array_fy(ind3-2)  = 0.0_wp
!            array_by(ind3)    = 0.0_wp
!            array_ey(ind3-1)  = 0.0_wp
!            array_ay(ind3)    = 0.0_wp
!            array_dy(ind3)    = 1.0_wp
!
!            ind3 = ind3max - 2
!            array_gy(ind3-3)  = 0.0_wp
!            array_fy(ind3-2)  = 0.0_wp
!            array_by(ind3)    = 0.0_wp
!            array_ey(ind3-1)  = 0.0_wp
!            array_ay(ind3)    = 0.0_wp
!            array_dy(ind3)    = 1.0_wp

            ! IRS1
            ind3 = ind3min
            coef     = theta_irs1 * deltat * idy_irs(ind3)
            array_cy(ind3) = 0.0_wp
            array_by(ind3)  = 0.0_wp
            array_ay(ind3)  = - coef * (r_spec6(ind2,ind3+1,ind1))**(unsur6)
            array_dy(ind3)  = 1.0_wp + coef * (r_spec6(ind2,ind3+1,ind1))**(unsur6)

            ind3 = ind3max
            coef     = theta_irs1 * deltat * idy_irs(ind3)
            array_dy(ind3)    = 1.0_wp + coef * (r_spec6(ind2,ind3,ind1))**(unsur6)
            array_ey(ind3-1)  = - coef * (r_spec6(ind2,ind3,ind1))**(unsur6)
            array_fy(ind3-2)  = 0.0_wp
            array_gy(ind3-3)  = 0.0_wp

            ! IRS2
            ind3 = ind3min + 1
            coef     = - theta_irs2 * (deltat * idy_irs(ind3))**2
            array_cy(ind3)  = 0.0_wp
            array_by(ind3)  = 0.0_wp
            array_ay(ind3)    = coef * (r_spec6(ind2,ind3+1,ind1))**(unsur6*2)
            array_ey(ind3-1)  = coef * (r_spec6(ind2,ind3,ind1))**(unsur6*2)
            array_dy(ind3)    = 1. - (array_ay(ind3) + array_ey(ind3-1))

            ind3 = ind3max - 1
            coef     = - theta_irs2 * (deltat * idy_irs(ind3))**2
            array_ay(ind3)    = coef * (r_spec6(ind2,ind3+1,ind1))**(unsur6*2)
            array_ey(ind3-1)  = coef * (r_spec6(ind2,ind3,ind1))**(unsur6*2)
            array_dy(ind3)    = 1.0_wp - (array_ay(ind3) + array_ey(ind3-1))
            array_fy(ind3-2)  = 0.0_wp
            array_gy(ind3-3)  = 0.0_wp

            ! IRS4
            ind3 = ind3min + 2
            coef     = theta_irs4 * (deltat * idy_irs(ind3))**4
            array_cy(ind3)    = 0.0_wp
            array_fy(ind3-2)  = coef * r_spec6(ind2,ind3,ind1)**(4*unsur6)
            array_by(ind3)    = coef * r_spec6(ind2,ind3+1,ind1)**(4*unsur6)
            array_ey(ind3-1)  =  - 3.0_wp * array_fy(ind3-2) - array_by(ind3)
            array_ay(ind3)    =  - 3.0_wp * array_by(ind3) - array_fy(ind3-2)
            array_dy(ind3)    = 1.0_wp + 3.0_wp * (array_fy(ind3-2) + array_by(ind3))

            ind3 = ind3max - 2
            coef     = theta_irs4 * (deltat * idy_irs(ind3))**4
            array_gy(ind3-3)  = 0.0_wp
            array_fy(ind3-2)  = coef * r_spec6(ind2,ind3,ind1)**(4*unsur6)
            array_by(ind3)    = coef * r_spec6(ind2,ind3+1,ind1)**(4*unsur6)
            array_ey(ind3-1)  =  - 3.0_wp * array_fy(ind3-2) - array_by(ind3)
            array_ay(ind3)    =  - 3.0_wp * array_by(ind3) - array_fy(ind3-2)
            array_dy(ind3)    = 1.0_wp + 3.0_wp * (array_fy(ind3-2) + array_by(ind3))

            ! Resolution of the heptadiagonal linear system
            ! =============================================
            call ls_irs6_gh(array_ay,array_by,array_cy,array_dy,array_ey,array_fy,array_gy,Krho(ind2,ind3min:ind3max,ind1),Krhou(ind2,ind3min:ind3max,ind1),&
                              Krhov(ind2,ind3min:ind3max,ind1),Krhow(ind2,ind3min:ind3max,ind1),Krhoe(ind2,ind3min:ind3max,ind1),ind3min,ind3max)
         enddo
      enddo
   endif

   !----------------------------------------
   ! Implicitation of direction z
   !----------------------------------------
   if (is_irs_k) then
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
               r_spec6(ind2,ind1,ind3) = (0.5 * (rspecmh + rspecph))**6
            enddo
         enddo
      enddo

      do ind1 = ind1min, ind1max
         do ind2 = ind2min, ind2max
            ! Preparation of the arrays and matrix
            ! ====================================
            do ind3 = ind3min+3,ind3max-3
               coef     = theta_irs6 * (deltat * idz_irs(ind3))**6
               coefmh   = coef * r_spec6(ind2,ind1,ind3)
               coefph   = coef * r_spec6(ind2,ind1,ind3+1)

               ! Third lower diagonal
               array_gz(ind3-3)  = -coefmh
               ! Third upper diagonal
               array_cz(ind3)    = -coefph
               ! Second lower diagonal
               array_fz(ind3-2)  = 5.0_wp*coefmh + coefph
               ! Second upper diagonal
               array_bz(ind3)    = coefmh + 5.0_wp*coefph
               ! First lower diagonal
               array_ez(ind3-1)  = -10.0_wp*coefmh - 5.0_wp*coefph
               ! First upper diagonal
               array_az(ind3)    = -5.0_wp*coefmh - 10.0_wp*coefph
               ! Main diagonal
               array_dz(ind3)    = 10.0_wp*(coefmh + coefph) + 1.0_wp
            enddo

!            ind3 = ind3min
!            array_cz(ind3) = 0.0_wp
!            array_bz(ind3)  = 0.0_wp
!            array_az(ind3)  = 0.0_wp
!            array_dz(ind3)  = 1.0_wp
!
!            ind3 = ind3max
!            array_dz(ind3)    = 1.0_wp
!            array_ez(ind3-1)  = 0.0_wp
!            array_fz(ind3-2)  = 0.0_wp
!            array_gz(ind3-3)  = 0.0_wp
!
!            ind3 = ind3min + 1
!            array_cz(ind3) = 0.0_wp
!            array_bz(ind3)  = 0.0_wp
!            array_az(ind3)    = 0.0_wp
!            array_ez(ind3-1)  = 0.0_wp
!            array_dz(ind3)    = 1.0_wp
!
!            ind3 = ind3max - 1
!            array_az(ind3)    = 0.0_wp
!            array_ez(ind3-1)  = 0.0_wp
!            array_dz(ind3)    = 1.0_wp
!            array_fz(ind3-2)  = 0.0_wp
!            array_gz(ind3-3)  = 0.0_wp
!
!            ind3 = ind3min + 2
!            array_cz(ind3)    = 0.0_wp
!            array_fz(ind3-2)  = 0.0_wp
!            array_bz(ind3)    = 0.0_wp
!            array_ez(ind3-1)  = 0.0_wp
!            array_az(ind3)    = 0.0_wp
!            array_dz(ind3)    = 1.0_wp
!
!            ind3 = ind3max - 2
!            array_gz(ind3-3)  = 0.0_wp
!            array_fz(ind3-2)  = 0.0_wp
!            array_bz(ind3)    = 0.0_wp
!            array_ez(ind3-1)  = 0.0_wp
!            array_az(ind3)    = 0.0_wp
!            array_dz(ind3)    = 1.0_wp

            ! IRS1
            ind3 = ind3min
            coef     = theta_irs1 * deltat * idz_irs(ind3)
            array_cz(ind3) = 0.0_wp
            array_bz(ind3)  = 0.0_wp
            array_az(ind3)  = - coef * (r_spec6(ind2,ind1,ind3+1))**(unsur6)
            array_dz(ind3)  = 1.0_wp + coef * (r_spec6(ind2,ind1,ind3+1))**(unsur6)

            ind3 = ind3max
            coef     = theta_irs1 * deltat * idz_irs(ind3)
            array_dz(ind3)    = 1.0_wp + coef * (r_spec6(ind2,ind1,ind3))**(unsur6)
            array_ez(ind3-1)  = - coef * (r_spec6(ind2,ind1,ind3))**(unsur6)
            array_fz(ind3-2)  = 0.0_wp
            array_gz(ind3-3)  = 0.0_wp

            ! IRS2
            ind3 = ind3min + 1
            coef     = - theta_irs2 * (deltat * idz_irs(ind3))**2
            array_cz(ind3)  = 0.0_wp
            array_bz(ind3)  = 0.0_wp
            array_az(ind3)    = coef * (r_spec6(ind2,ind1,ind3+1))**(unsur6*2)
            array_ez(ind3-1)  = coef * (r_spec6(ind2,ind1,ind3))**(unsur6*2)
            array_dz(ind3)    = 1. - (array_az(ind3) + array_ez(ind3-1))

            ind3 = ind3max - 1
            coef     = - theta_irs2 * (deltat * idz_irs(ind3))**2
            array_az(ind3)    = coef * (r_spec6(ind2,ind1,ind3+1))**(unsur6*2)
            array_ez(ind3-1)  = coef * (r_spec6(ind2,ind1,ind3))**(unsur6*2)
            array_dz(ind3)    = 1.0_wp - (array_az(ind3) + array_ez(ind3-1))
            array_fz(ind3-2)  = 0.0_wp
            array_gz(ind3-3)  = 0.0_wp

            ! IRS4
            ind3 = ind3min + 2
            coef     = theta_irs4 * (deltat * idz_irs(ind3))**4
            array_cz(ind3)    = 0.0_wp
            array_fz(ind3-2)  = coef * r_spec6(ind2,ind1,ind3)**(4*unsur6)
            array_bz(ind3)    = coef * r_spec6(ind2,ind1,ind3+1)**(4*unsur6)
            array_ez(ind3-1)  =  - 3.0_wp * array_fz(ind3-2) - array_bz(ind3)
            array_az(ind3)    =  - 3.0_wp * array_bz(ind3) - array_fz(ind3-2)
            array_dz(ind3)    = 1.0_wp + 3.0_wp * (array_fz(ind3-2) + array_bz(ind3))

            ind3 = ind3max - 2
            coef     = theta_irs4 * (deltat * idz_irs(ind3))**4
            array_gz(ind3-3)  = 0.0_wp
            array_fz(ind3-2)  = coef * r_spec6(ind2,ind1,ind3)**(4*unsur6)
            array_bz(ind3)    = coef * r_spec6(ind2,ind1,ind3+1)**(4*unsur6)
            array_ez(ind3-1)  =  - 3.0_wp * array_fz(ind3-2) - array_bz(ind3)
            array_az(ind3)    =  - 3.0_wp * array_bz(ind3) - array_fz(ind3-2)
            array_dz(ind3)    = 1.0_wp + 3.0_wp * (array_fz(ind3-2) + array_bz(ind3))

            ! Resolution of the heptadiagonal linear system
            ! =============================================
            call ls_irs6_gh(array_az,array_bz,array_cz,array_dz,array_ez,array_fz,array_gz,Krho(ind2,ind1,ind3min:ind3max),Krhou(ind2,ind1,ind3min:ind3max),&
                              Krhov(ind2,ind1,ind3min:ind3max),Krhow(ind2,ind1,ind3min:ind3max),Krhoe(ind2,ind1,ind3min:ind3max),ind3min,ind3max)
         enddo
      enddo
   endif

end subroutine irs6_ngh_v0

!===============================================================================
subroutine irs4_sca_v0
!===============================================================================
!> Preparation of matrix for the band matrix scalapack parallelisation for IRS4
!===============================================================================
   use mod_flow
   use mod_constant
   use mod_time
   use mod_mpi
   use warnstop
   implicit none
   integer  :: ind1,ind2,ind3
   integer  :: ind1min,ind1max,ind2min,ind2max,ind3min,ind3max
   real(wp), dimension(1:nx+1,1:ny+1,1:nz+1)    :: r_spec4
   real(wp), dimension(1:nx,1:ny,1:nz,5)        :: RHS
   real(wp) :: coef,coefmh,coefph,rspecmh,rspecph

   ! Remplissage de la matrice RHS
   !     /!\ Temporaire /!\
   ! =============================
   RHS(1:nx,1:ny,1:nz,1)   = Krho(1:nx,1:ny,1:nz)
   RHS(1:nx,1:ny,1:nz,2)   = Krhou(1:nx,1:ny,1:nz)
   RHS(1:nx,1:ny,1:nz,3)   = Krhov(1:nx,1:ny,1:nz)
   RHS(1:nx,1:ny,1:nz,4)   = Krhow(1:nx,1:ny,1:nz)
   RHS(1:nx,1:ny,1:nz,5)   = Krhoe(1:nx,1:ny,1:nz)

   !----------------------------------------
   ! Implicitation of direction x
   !----------------------------------------
   if (is_irs_i) then
!      ind3 : i;     ind1 : k;     ind2 : j
      ind1min = 1;   ind2min = 1;   ind3min = 1
      ind1max = nz;  ind2max = ny;  ind3max = nx

      ! Calculation of the spectral radius at ind-1/2
      ! =============================================
      if (.not. is_2d) then
         do ind1=ind1min,ind1max
            do ind2=ind2min,ind2max
               do ind3=ind3min,ind3max+1
                  rspecmh = sqrt(uu(ind3-1,ind2,ind1)**2)   + c_(ind3-1,ind2,ind1)
                  rspecph = sqrt(uu(ind3,ind2,ind1)**2)     + c_(ind3,ind2,ind1)
                  r_spec4(ind3,ind2,ind1) = (0.5 * (rspecmh + rspecph))**4
               enddo
            enddo
         enddo
      else
         do ind2=ind2min,ind2max
            do ind3=ind3min,ind3max+1
               rspecmh = sqrt(uu(ind3-1,ind2,1)**2)   + c_(ind3-1,ind2,1)
               rspecph = sqrt(uu(ind3,ind2,1)**2)     + c_(ind3,ind2,1)
               r_spec4(ind3,ind2,1) = (0.5 * (rspecmh + rspecph))**4
            enddo
         enddo
      endif

      do ind1 = ind1min, ind1max
         do ind2 = ind2min, ind2max
            ! Preparation of the arrays and matrix
            ! ====================================
            do ind3 = ind3min,ind3max
               coefmh   = deltat * idx_irs(ind3)
               coefph   = deltat * idx_irs(ind3+1)

               coef     = theta_irs4 * deltat * 0.5_wp*(idx_irs(ind3)+idx_irs(ind3+1))

               ! Second lower diagonal
               A_x(1,ind3)  = coef * coefmh**3 * r_spec4(ind3,ind2,ind1)
               ! Second upper diagonal
               A_x(5,ind3)  = coef * coefph**3 * r_spec4(ind3+1,ind2,ind1)
               ! First lower diagonal
               A_x(2,ind3)  = - 3.0_wp * A_x(1,ind3) - A_x(5,ind3)
               ! First upper diagonal
               A_x(4,ind3)  = - 3.0_wp * A_x(5,ind3) - A_x(1,ind3)
               ! Main diagonal
               A_x(3,ind3) = 1.0_wp + 3.0_wp * (A_x(1,ind3) + A_x(5,ind3))
            enddo

!            if (is_boundary(1,1)) then
!               ! A COMPLETER
!               !??
!            endif
!            if (is_boundary(1,2)) then
!               ! A COMPLETER
!               !??
!            endif

            ! Resolution of the pentadiagonal linear system
            ! =============================================
            call pddbtrf(ngx,bw,bw,A_x,JA,DESCA_x,AF_x,LAF_x,WORK,LWORK,info)
            call pddbtrs('T',ngx,bw,bw,NRHS,A_x,JA,DESCA_x,RHS(:,ind2,ind1,:),IB,DESCB_x,AF_x,LAF_x,WORK,LWORK,INFO)

         enddo
      enddo
   endif

   !----------------------------------------
   ! Implicitation of direction y
   !----------------------------------------
   if (is_irs_j) then
!      ind3 : j;     ind1 : k;     ind2 : i
      ind1min = 1;   ind2min = 1;   ind3min = 1
      ind1max = nz;  ind2max = nx;  ind3max = ny

      ! Calculation of the spectral radius at ind-1/2
      ! =============================================
      if (.not. is_2d) then
         do ind1=ind1min,ind1max
            do ind2=ind2min,ind2max
               do ind3=ind3min,ind3max+1
                  rspecmh = sqrt(vv(ind2,ind3-1,ind1)**2)     + c_(ind2,ind3-1,ind1)
                  rspecph = sqrt(vv(ind2,ind3,ind1)**2)     + c_(ind2,ind3,ind1)
                  r_spec4(ind2,ind3,ind1) = (0.5 * (rspecmh + rspecph))**4
               enddo
            enddo
         enddo
      else
         do ind2=ind2min,ind2max
            do ind3=ind3min,ind3max+1
               rspecmh = sqrt(vv(ind2,ind3-1,1)**2)   + c_(ind2,ind3-1,1)
               rspecph = sqrt(vv(ind2,ind3,1)**2)     + c_(ind2,ind3,1)
               r_spec4(ind2,ind3,1) = (0.5 * (rspecmh + rspecph))**4
            enddo
         enddo
      endif

      do ind1 = ind1min, ind1max
         do ind2 = ind2min, ind2max
            ! Preparation of the arrays and matrix
            ! ====================================
            do ind3 = ind3min,ind3max
               coefmh   = deltat * idy_irs(ind3)
               coefph   = deltat * idy_irs(ind3+1)

               coef     = theta_irs4 * deltat * 0.5_wp*(idy_irs(ind3)+idy_irs(ind3+1))

               ! Second lower diagonal
               A_y(1,ind3)  = coef * coefmh**3 * r_spec4(ind2,ind3,ind1)
               ! Second upper diagonal
               A_y(5,ind3)  = coef * coefph**3 * r_spec4(ind2,ind3+1,ind1)
               ! First lower diagonal
               A_y(2,ind3)  = - 3.0_wp * A_y(1,ind3) - A_y(5,ind3)
               ! First upper diagonal
               A_y(4,ind3)  = - 3.0_wp * A_y(5,ind3) - A_y(1,ind3)
               ! Main diagonal
               A_y(3,ind3) = 1.0_wp + 3.0_wp * (A_y(1,ind3) + A_y(5,ind3))
            enddo

!            if (is_boundary(1,1)) then
!               ! A COMPLETER
!               !??
!            endif
!            if (is_boundary(1,2)) then
!               ! A COMPLETER
!               !??
!            endif

            ! Resolution of the pentadiagonal linear system
            ! =============================================
            call pddbtrf(ngy,bw,bw,A_y,JA,DESCA_y,AF_y,LAF_y,WORK,LWORK,info)
            call pddbtrs('T',ngy,bw,bw,NRHS,A_y,JA,DESCA_y,RHS(ind2,:,ind1,:),IB,DESCB_y,AF_y,LAF_y,WORK,LWORK,INFO)

         enddo
      enddo
   endif

   !----------------------------------------
   ! Implicitation of direction z
   !----------------------------------------
   if (is_irs_k) then
      call mpistop('Writing of the implicitation of the z direction not ready yet. Shutting down...',0)
   endif


   ! Remplissage des flux a partir du RHS
   !           /!\ Temporaire /!\
   ! ====================================
   Krho(1:nx,1:ny,1:nz)    = RHS(1:nx,1:ny,1:nz,1)
   Krhou(1:nx,1:ny,1:nz)   = RHS(1:nx,1:ny,1:nz,2)
   Krhov(1:nx,1:ny,1:nz)   = RHS(1:nx,1:ny,1:nz,3)
   Krhow(1:nx,1:ny,1:nz)   = RHS(1:nx,1:ny,1:nz,4)
   Krhoe(1:nx,1:ny,1:nz)   = RHS(1:nx,1:ny,1:nz,5)

end subroutine irs4_sca_v0

!===============================================================================
subroutine irs6_sca_v0
!===============================================================================
!> Preparation of matrix for the band matrix scalapack parallelisation for IRS6
!===============================================================================
   use mod_flow
   use mod_constant
   use mod_time
   use mod_mpi
   use warnstop
   implicit none
   integer  :: ind1,ind2,ind3
   integer  :: ind1min,ind1max,ind2min,ind2max,ind3min,ind3max
   real(wp), dimension(1:nx+1,1:ny+1,1:nz+1)    :: r_spec6
   real(wp), dimension(1:nx,1:ny,1:nz,5)        :: RHS
   real(wp) :: coef,coefmh,coefph,rspecmh,rspecph

   ! Remplissage de la matrice RHS
   !     /!\ Temporaire /!\
   ! =============================
   RHS(1:nx,1:ny,1:nz,1)   = Krho(1:nx,1:ny,1:nz)
   RHS(1:nx,1:ny,1:nz,2)   = Krhou(1:nx,1:ny,1:nz)
   RHS(1:nx,1:ny,1:nz,3)   = Krhov(1:nx,1:ny,1:nz)
   RHS(1:nx,1:ny,1:nz,4)   = Krhow(1:nx,1:ny,1:nz)
   RHS(1:nx,1:ny,1:nz,5)   = Krhoe(1:nx,1:ny,1:nz)

   !----------------------------------------
   ! Implicitation of direction x
   !----------------------------------------
   if (is_irs_i) then
!      ind3 : i;     ind1 : k;     ind2 : j
      ind1min = 1;   ind2min = 1;   ind3min = 1
      ind1max = nz;  ind2max = ny;  ind3max = nx

      ! Calculation of the spectral radius at ind-1/2
      ! =============================================
      if (.not. is_2d) then
         do ind1=ind1min,ind1max
            do ind2=ind2min,ind2max
               do ind3=ind3min,ind3max+1
                  rspecmh = sqrt(uu(ind3-1,ind2,ind1)**2)   + c_(ind3-1,ind2,ind1)
                  rspecph = sqrt(uu(ind3,ind2,ind1)**2)     + c_(ind3,ind2,ind1)
                  r_spec6(ind3,ind2,ind1) = (0.5 * (rspecmh + rspecph))**6
               enddo
            enddo
         enddo
      else
         do ind2=ind2min,ind2max
            do ind3=ind3min,ind3max+1
               rspecmh = sqrt(uu(ind3-1,ind2,1)**2)   + c_(ind3-1,ind2,1)
               rspecph = sqrt(uu(ind3,ind2,1)**2)     + c_(ind3,ind2,1)
               r_spec6(ind3,ind2,1) = (0.5 * (rspecmh + rspecph))**6
            enddo
         enddo
      endif

      do ind1 = ind1min, ind1max
         do ind2 = ind2min, ind2max
            ! Preparation of the arrays and matrix
            ! ====================================
            do ind3 = ind3min,ind3max
               coef     = theta_irs6 * (deltat * idx(ind3))**6
               coefmh   = coef * r_spec6(ind3,ind2,ind1)
               coefph   = coef * r_spec6(ind3+1,ind2,ind1)

               ! Third lower diagonal
               A_x(1,ind3) = -coefmh
               ! Third upper diagonal
               A_x(7,ind3) = -coefph
               ! Second lower diagonal
               A_x(2,ind3) = 5*coefmh + coefph
               ! Second upper diagonal
               A_x(6,ind3) = coefmh + 5*coefph
               ! First lower diagonal
               A_x(3,ind3) = -10*coefmh - 5*coefph
               ! First upper diagonal
               A_x(5,ind3) = -5*coefmh - 10*coefph
               ! Main diagonal
               A_x(4,ind3) = 10*(coefmh + coefph) + 1
            enddo

!            if (is_boundary(1,1)) then
!               ! A COMPLETER
!               !??
!            endif
!            if (is_boundary(1,2)) then
!               ! A COMPLETER
!               !??
!            endif

            ! Resolution of the heptadiagonal linear system
            ! =============================================
            call pddbtrf(ngx,bw,bw,A_x,JA,DESCA_x,AF_x,LAF_x,WORK,LWORK,info)
            call pddbtrs('T',ngx,bw,bw,NRHS,A_x,JA,DESCA_x,RHS(:,ind2,ind1,:),IB,DESCB_x,AF_x,LAF_x,WORK,LWORK,INFO)

         enddo
      enddo
   endif

   !----------------------------------------
   ! Implicitation of direction y
   !----------------------------------------
   if (is_irs_j) then
!      ind3 : j;     ind1 : k;     ind2 : i
      ind1min = 1;   ind2min = 1;   ind3min = 1
      ind1max = nz;  ind2max = nx;  ind3max = ny

      ! Calculation of the spectral radius at ind-1/2
      ! =============================================
      if (.not. is_2d) then
         do ind1=ind1min,ind1max
            do ind2=ind2min,ind2max
               do ind3=ind3min,ind3max+1
                  rspecmh = sqrt(vv(ind2,ind3-1,ind1)**2)   + c_(ind2,ind3-1,ind1)
                  rspecph = sqrt(vv(ind2,ind3,ind1)**2)     + c_(ind2,ind3,ind1)
                  r_spec6(ind2,ind3,ind1) = (0.5 * (rspecmh + rspecph))**6
               enddo
            enddo
         enddo
      else
         do ind2=ind2min,ind2max
            do ind3=ind3min,ind3max+1
               rspecmh = sqrt(vv(ind2,ind3-1,1)**2)   + c_(ind2,ind3-1,1)
               rspecph = sqrt(vv(ind2,ind3,1)**2)     + c_(ind2,ind3,1)
               r_spec6(ind2,ind3,1) = (0.5 * (rspecmh + rspecph))**6
            enddo
         enddo
      endif

      do ind1 = ind1min, ind1max
         do ind2 = ind2min, ind2max
            ! Preparation of the arrays and matrix
            ! ====================================
            do ind3 = ind3min,ind3max
               coef     = theta_irs6 * (deltat * idy(ind3))**6
               coefmh   = coef * r_spec6(ind2,ind3,ind1)
               coefph   = coef * r_spec6(ind2,ind3+1,ind1)

               ! Third lower diagonal
               A_y(1,ind3) = -coefmh
               ! Third upper diagonal
               A_y(7,ind3) = -coefph
               ! Second lower diagonal
               A_y(2,ind3) = 5*coefmh + coefph
               ! Second upper diagonal
               A_y(6,ind3) = coefmh + 5*coefph
               ! First lower diagonal
               A_y(3,ind3) = -10*coefmh - 5*coefph
               ! First upper diagonal
               A_y(5,ind3) = -5*coefmh - 10*coefph
               ! Main diagonal
               A_y(4,ind3) = 10*coefmh + 10*coefph + 1
            enddo

!            if (is_boundary(1,1)) then
!               ! A COMPLETER
!               !??
!            endif
!            if (is_boundary(1,2)) then
!               ! A COMPLETER
!               !??
!            endif

            ! Resolution of the heptadiagonal linear system
            ! =============================================
            call pddbtrf(ngy,bw,bw,A_y,JA,DESCA_y,AF_y,LAF_y,WORK,LWORK,info)
            call pddbtrs('T',ngy,bw,bw,NRHS,A_y,JA,DESCA_y,RHS(ind2,:,ind1,:),IB,DESCB_y,AF_y,LAF_y,WORK,LWORK,INFO)

         enddo
      enddo
   endif

   !----------------------------------------
   ! Implicitation of direction z
   !----------------------------------------
   if (is_irs_k) then
      call mpistop('Writing of the implicitation of the z direction not ready yet. Shutting down...',0)
   endif


   ! Remplissage des flux a partir du RHS
   !           /!\ Temporaire /!\
   ! ====================================
   Krho(1:nx,1:ny,1:nz)    = RHS(1:nx,1:ny,1:nz,1)
   Krhou(1:nx,1:ny,1:nz)   = RHS(1:nx,1:ny,1:nz,2)
   Krhov(1:nx,1:ny,1:nz)   = RHS(1:nx,1:ny,1:nz,3)
   Krhow(1:nx,1:ny,1:nz)   = RHS(1:nx,1:ny,1:nz,4)
   Krhoe(1:nx,1:ny,1:nz)   = RHS(1:nx,1:ny,1:nz,5)

end subroutine irs6_sca_v0

!===============================================================================
subroutine irs8_sca_v0
!===============================================================================
!> Preparation of matrix for the band matrix scalapack parallelisation for IRS8
!===============================================================================
   use mod_flow
   use mod_constant
   use mod_time
   use mod_mpi
   use warnstop
   implicit none
   integer  :: ind1,ind2,ind3
   integer  :: ind1min,ind1max,ind2min,ind2max,ind3min,ind3max
   real(wp), dimension(1:nx+1,1:ny+1,1:nz+1)    :: r_spec8
   real(wp), dimension(1:nx,1:ny,1:nz,5)        :: RHS
   real(wp) :: coef,coefmh,coefph,rspecmh,rspecph

   ! Remplissage de la matrice RHS
   !     /!\ Temporaire /!\
   ! =============================
   RHS(1:nx,1:ny,1:nz,1)   = Krho(1:nx,1:ny,1:nz)
   RHS(1:nx,1:ny,1:nz,2)   = Krhou(1:nx,1:ny,1:nz)
   RHS(1:nx,1:ny,1:nz,3)   = Krhov(1:nx,1:ny,1:nz)
   RHS(1:nx,1:ny,1:nz,4)   = Krhow(1:nx,1:ny,1:nz)
   RHS(1:nx,1:ny,1:nz,5)   = Krhoe(1:nx,1:ny,1:nz)

   !----------------------------------------
   ! Implicitation of direction x
   !----------------------------------------
   if (is_irs_i) then
!      ind3 : i;     ind1 : k;     ind2 : j
      ind1min = 1;   ind2min = 1;   ind3min = 1
      ind1max = nz;  ind2max = ny;  ind3max = nx

      ! Calculation of the spectral radius at ind-1/2
      ! =============================================
      if (.not. is_2d) then
         do ind1=ind1min,ind1max
            do ind2=ind2min,ind2max
               do ind3=ind3min,ind3max+1
                  rspecmh = sqrt(uu(ind3-1,ind2,ind1)**2)   + c_(ind3-1,ind2,ind1)
                  rspecph = sqrt(uu(ind3,ind2,ind1)**2)     + c_(ind3,ind2,ind1)
                  r_spec8(ind3,ind2,ind1) = (0.5 * (rspecmh + rspecph))**8
               enddo
            enddo
         enddo
      else
         do ind2=ind2min,ind2max
            do ind3=ind3min,ind3max+1
               rspecmh = sqrt(uu(ind3-1,ind2,1)**2)   + c_(ind3-1,ind2,1)
               rspecph = sqrt(uu(ind3,ind2,1)**2)     + c_(ind3,ind2,1)
               r_spec8(ind3,ind2,1) = (0.5 * (rspecmh + rspecph))**8
            enddo
         enddo
      endif

      do ind1 = ind1min, ind1max
         do ind2 = ind2min, ind2max
            ! Preparation of the arrays and matrix
            ! ====================================
            do ind3 = ind3min,ind3max
               coef     = theta_irs8 * (deltat * idx(ind3))**8
               coefmh   = coef * r_spec8(ind3,ind2,ind1)
               coefph   = coef * r_spec8(ind3+1,ind2,ind1)

               ! Fourth lower diagonal
               A_x(1,ind3) = coefmh
               ! Fourth upper diagonal
               A_x(9,ind3) = coefph
               ! Third lower diagonal
               A_x(2,ind3) = -7*coefmh - coefph
               ! Third upper diagonal
               A_x(8,ind3) = -coefmh - 7*coefph
               ! Second lower diagonal
               A_x(3,ind3) = 21*coefmh + 7*coefph
               ! Second upper diagonal
               A_x(7,ind3) = 7*coefmh + 21*coefph
               ! First lower diagonal
               A_x(4,ind3) = -35*coefmh - 21*coefph
               ! First upper diagonal
               A_x(6,ind3) = -21*coefmh - 35*coefph
               ! Main diagonal
               A_x(5,ind3) = 35*(coefmh + coefph) + 1
            enddo

!            if (is_boundary(1,1)) then
!               ! A COMPLETER
!               !??
!            endif
!            if (is_boundary(1,2)) then
!               ! A COMPLETER
!               !??
!            endif

            ! Resolution of the enneadiagonal linear system
            ! =============================================
            call pddbtrf(ngx,bw,bw,A_x,JA,DESCA_x,AF_x,LAF_x,WORK,LWORK,info)
            call pddbtrs('T',ngx,bw,bw,NRHS,A_x,JA,DESCA_x,RHS(:,ind2,ind1,:),IB,DESCB_x,AF_x,LAF_x,WORK,LWORK,INFO)

         enddo
      enddo
   endif

   !----------------------------------------
   ! Implicitation of direction y
   !----------------------------------------
   if (is_irs_j) then
!      ind3 : j;     ind1 : k;     ind2 : i
      ind1min = 1;   ind2min = 1;   ind3min = 1
      ind1max = nz;  ind2max = nx;  ind3max = ny

      ! Calculation of the spectral radius at ind-1/2
      ! =============================================
      if (.not. is_2d) then
         do ind1=ind1min,ind1max
            do ind2=ind2min,ind2max
               do ind3=ind3min,ind3max+1
                  rspecmh = sqrt(vv(ind2,ind3-1,ind1)**2)   + c_(ind2,ind3-1,ind1)
                  rspecph = sqrt(vv(ind2,ind3,ind1)**2)     + c_(ind2,ind3,ind1)
                  r_spec8(ind2,ind3,ind1) = (0.5 * (rspecmh + rspecph))**8
               enddo
            enddo
         enddo
      else
         do ind2=ind2min,ind2max
            do ind3=ind3min,ind3max+1
               rspecmh = sqrt(vv(ind2,ind3-1,1)**2)   + c_(ind2,ind3-1,1)
               rspecph = sqrt(vv(ind2,ind3,1)**2)     + c_(ind2,ind3,1)
               r_spec8(ind2,ind3,1) = (0.5 * (rspecmh + rspecph))**8
            enddo
         enddo
      endif

      do ind1 = ind1min, ind1max
         do ind2 = ind2min, ind2max
            ! Preparation of the arrays and matrix
            ! ====================================
            do ind3 = ind3min,ind3max
               coef     = theta_irs8 * (deltat * idx(ind3))**8
               coefmh   = coef * r_spec8(ind2,ind3,ind1)
               coefph   = coef * r_spec8(ind2,ind3+1,ind1)

               ! Fourth lower diagonal
               A_y(1,ind3) = coefmh
               ! Fourth upper diagonal
               A_y(9,ind3) = coefph
               ! Third lower diagonal
               A_y(2,ind3) = -7*coefmh - coefph
               ! Third upper diagonal
               A_y(8,ind3) = -coefmh - 7*coefph
               ! Second lower diagonal
               A_y(3,ind3) = 21*coefmh + 7*coefph
               ! Second upper diagonal
               A_y(7,ind3) = 7*coefmh + 21*coefph
               ! First lower diagonal
               A_y(4,ind3) = -35*coefmh - 21*coefph
               ! First upper diagonal
               A_y(6,ind3) = -21*coefmh - 35*coefph
               ! Main diagonal
               A_y(5,ind3) = 35*(coefmh + coefph) + 1
            enddo

!            if (is_boundary(1,1)) then
!               ! A COMPLETER
!               !??
!            endif
!            if (is_boundary(1,2)) then
!               ! A COMPLETER
!               !??
!            endif

            ! Resolution of the enneadiagonal linear system
            ! =============================================
            call pddbtrf(ngy,bw,bw,A_y,JA,DESCA_y,AF_y,LAF_y,WORK,LWORK,info)
            call pddbtrs('T',ngy,bw,bw,NRHS,A_y,JA,DESCA_y,RHS(ind2,:,ind1,:),IB,DESCB_y,AF_y,LAF_y,WORK,LWORK,INFO)

         enddo
      enddo
   endif

   !----------------------------------------
   ! Implicitation of direction z
   !----------------------------------------
   if (is_irs_k) then
      call mpistop('Writing of the implicitation of the z direction not ready yet. Shutting down...',0)
   endif


   ! Remplissage des flux a partir du RHS
   !           /!\ Temporaire /!\
   ! ====================================
   Krho(1:nx,1:ny,1:nz)    = RHS(1:nx,1:ny,1:nz,1)
   Krhou(1:nx,1:ny,1:nz)   = RHS(1:nx,1:ny,1:nz,2)
   Krhov(1:nx,1:ny,1:nz)   = RHS(1:nx,1:ny,1:nz,3)
   Krhow(1:nx,1:ny,1:nz)   = RHS(1:nx,1:ny,1:nz,4)
   Krhoe(1:nx,1:ny,1:nz)   = RHS(1:nx,1:ny,1:nz,5)

end subroutine irs8_sca_v0

!===============================================================================
subroutine irs_sca_tri_v0
!===============================================================================
!> Application of the tridiagonal scalapack parallelisation
!===============================================================================
   use mod_flow
   use mod_constant
   use mod_time
   use mod_mpi
   use warnstop
   implicit none
   integer          :: ind1,ind2,ind3
   integer :: ind1min,ind1max,ind2min,ind2max,ind3min,ind3max
   real(wp), dimension(1:nx+1,1:ny+1,1:nz+1)    :: r_spec2
   real(wp), dimension(1:nx,1:ny,1:nz,5)        :: RHS
   real(wp) :: coef,rspecmh,rspecph

   ! Remplissage de la matrice RHS
   !     /!\ Temporaire /!\
   ! =============================
   RHS(1:nx,1:ny,1:nz,1)   = Krho(1:nx,1:ny,1:nz)
   RHS(1:nx,1:ny,1:nz,2)   = Krhou(1:nx,1:ny,1:nz)
   RHS(1:nx,1:ny,1:nz,3)   = Krhov(1:nx,1:ny,1:nz)
   RHS(1:nx,1:ny,1:nz,4)   = Krhow(1:nx,1:ny,1:nz)
   RHS(1:nx,1:ny,1:nz,5)   = Krhoe(1:nx,1:ny,1:nz)

   !----------------------------------------
   ! Implicitation of direction x
   !----------------------------------------
   if (is_irs_i) then
!      ind3 : i;     ind1 : k;     ind2 : j
      ind1min = 1;   ind2min = 1;   ind3min = 1
      ind1max = nz;  ind2max = ny;  ind3max = nx

      ! Calculation of the spectral radius at ind-1/2
      ! =============================================
      if (.not. is_2d) then
         do ind1=ind1min,ind1max
            do ind2=ind2min,ind2max
               do ind3=ind3min,ind3max+1
                  rspecmh = sqrt(uu(ind3-1,ind2,ind1)**2) + c_(ind3-1,ind2,ind1)
                  rspecph = sqrt(uu(ind3,ind2,ind1)**2) + c_(ind3,ind2,ind1)
                  r_spec2(ind3,ind2,ind1) = (0.5 * (rspecmh + rspecph))**2
               enddo
            enddo
         enddo
      else
         do ind2=ind2min,ind2max
            do ind3=ind3min,ind3max+1
               rspecmh = sqrt(uu(ind3-1,ind2,1)**2) + c_(ind3-1,ind2,1)
               rspecph = sqrt(uu(ind3,ind2,1)**2) + c_(ind3,ind2,1)
               r_spec2(ind3,ind2,1) = (0.5 * (rspecmh + rspecph))**2
            enddo
         enddo
      endif

      do ind1 = ind1min, ind1max
         do ind2 = ind2min, ind2max
            ! Preparation of the arrays and matrix
            ! ====================================
            do ind3 = ind3min,ind3max
               coef = - theta_irs2 * (deltat * idx(ind3))**2

               ! Lower diagonal
               DL_x(ind3)  = coef * r_spec2(ind3,ind2,ind1)
               ! Upper diagonal
               DU_x(ind3)  = coef * r_spec2(ind3+1,ind2,ind1)
               ! Main diagonal
               D_x(ind3) = 1  - (DL_x(ind3) + DU_x(ind3))
            enddo

!            if (is_boundary(1,1)) then
!               ! A COMPLETER
!               !??
!            endif
!            if (is_boundary(1,2)) then
!               ! A COMPLETER
!               !??
!            endif

            ! Resolution of the tridiagonal linear system
            ! ===========================================
            call pddttrf(ngx,DL_x,D_x,DU_x,JA,DESCA_x,AF_x,LAF_x,WORK,LWORK,info)
            call pddttrs('N',ngx,NRHS,DL_x,D_x,DU_x,JA,DESCA_x,RHS(:,ind2,ind1,:),IB,DESCB_x,AF_x,LAF_x,WORK,LWORK,INFO)

         enddo
      enddo
   endif

   !----------------------------------------
   ! Implicitation of direction y
   !----------------------------------------
   if (is_irs_j) then
!      ind3 : j;     ind1 : k;     ind2 : i
      ind1min = 1;   ind2min = 1;   ind3min = 1
      ind1max = nz;  ind2max = nx;  ind3max = ny

      ! Calculation of the spectral radius at ind-1/2
      ! =============================================
      if (.not. is_2d) then
         do ind1=ind1min,ind1max
            do ind2=ind2min,ind2max
               do ind3=ind3min,ind3max+1
                  rspecmh = sqrt(vv(ind2,ind3-1,ind1)**2) + c_(ind2,ind3-1,ind1)
                  rspecph = sqrt(vv(ind2,ind3,ind1)**2) + c_(ind2,ind3,ind1)
                  r_spec2(ind2,ind3,ind1) = (0.5 * (rspecmh + rspecph))**2
               enddo
            enddo
         enddo
      else
         do ind2=ind2min,ind2max
            do ind3=ind3min,ind3max+1
               rspecmh = sqrt(vv(ind2,ind3-1,1)**2) + c_(ind2,ind3-1,1)
               rspecph = sqrt(vv(ind2,ind3,1)**2) + c_(ind2,ind3,1)
               r_spec2(ind2,ind3,1) = (0.5 * (rspecmh + rspecph))**2
            enddo
         enddo
      endif

      do ind1 = ind1min, ind1max
         do ind2 = ind2min, ind2max
            ! Preparation of the arrays and matrix
            ! ====================================
            do ind3 = ind3min,ind3max
               coef = - theta_irs2 * (deltat * idy(ind3))**2

               ! Lower diagonal
               DL_y(ind3)  = coef * r_spec2(ind2,ind3,ind1)
               ! Upper diagonal
               DU_y(ind3)  = coef * r_spec2(ind2,ind3+1,ind1)
               ! Main diagonal
               D_y(ind3) = 1  - (DL_y(ind3) + DU_y(ind3))
            enddo

!            if (is_boundary(1,1)) then
!               ! A COMPLETER
!               !??
!            endif
!            if (is_boundary(1,2)) then
!               ! A COMPLETER
!               !??
!            endif

            ! Resolution of the tridiagonal linear system
            ! ===========================================
            call pddttrf(ngy,DL_y,D_y,DU_y,JA,DESCA_y,AF_y,LAF_y,WORK,LWORK,info)
            call pddttrs('N',ngy,NRHS,DL_y,D_y,DU_y,JA,DESCA_y,RHS(ind2,:,ind1,:),IB,DESCB_y,AF_y,LAF_y,WORK,LWORK,INFO)
         enddo
      enddo
   endif

   !----------------------------------------
   ! Implicitation of direction z
   !----------------------------------------
   if (is_irs_k .and. .not. is_2d) then
      call mpistop('Writing of the implicitation of the z direction not ready yet. Shutting down...',0)
   endif

   ! Remplissage des flux a partir du RHS
   !           /!\ Temporaire /!\
   ! ====================================
   Krho(1:nx,1:ny,1:nz)    = RHS(1:nx,1:ny,1:nz,1)
   Krhou(1:nx,1:ny,1:nz)   = RHS(1:nx,1:ny,1:nz,2)
   Krhov(1:nx,1:ny,1:nz)   = RHS(1:nx,1:ny,1:nz,3)
   Krhow(1:nx,1:ny,1:nz)   = RHS(1:nx,1:ny,1:nz,4)
   Krhoe(1:nx,1:ny,1:nz)   = RHS(1:nx,1:ny,1:nz,5)

end subroutine irs_sca_tri_v0

!===============================================================================
subroutine irs_pascal_v0
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
   real(wp), dimension(1:nx+1,1:ny+1,1:nz+1) :: r_spec2
   real(wp), dimension(1:ny,1:nx,1:NRHS)     :: RHS_x
   real(wp) :: coef,rspecmh,rspecph

   !----------------------------------------
   ! Implicitation of direction x
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

      ! Calculation of the spectral radius at ind-1/2
      ! =============================================
      if (.not. is_2d) then
         do ind1=ind1min,ind1max
            do ind2=ind2min,ind2max
               do ind3=ind3min,ind3max+1
                  rspecmh = sqrt(uu(ind3-1,ind2,ind1)**2) + c_(ind3-1,ind2,ind1)
                  rspecph = sqrt(uu(ind3,ind2,ind1)**2) + c_(ind3,ind2,ind1)
                  r_spec2(ind3,ind2,ind1) = (0.5 * (rspecmh + rspecph))**2
               enddo
            enddo
         enddo
      else
         do ind2=ind2min,ind2max
            do ind3=ind3min,ind3max+1
               rspecmh = sqrt(uu(ind3-1,ind2,1)**2) + c_(ind3-1,ind2,1)
               rspecph = sqrt(uu(ind3,ind2,1)**2) + c_(ind3,ind2,1)
               r_spec2(ind3,ind2,1) = (0.5 * (rspecmh + rspecph))**2
            enddo
         enddo
      endif


      do ind1 = ind1min, ind1max
         ! Preparation of the arrays and matrix
         ! ====================================
         do irhs = 1,NRHS
            do ind2 = ind2min, ind2max
               do ind3 = ind3min,ind3max
                  coef = - theta_irs2 * (deltat * idx(ind3))**2

                  ! Lower diagonal
                  DLp_x(ind2,ind3,irhs) = coef * r_spec2(ind3,ind2,ind1)
                  ! Upper diagonal
                  DUp_x(ind2,ind3,irhs) = coef * r_spec2(ind3+1,ind2,ind1)
                  ! Main diagonal
                  Dp_x(ind2,ind3,irhs)  = 1  - (DLp_x(ind2,ind3,irhs) + DUp_x(ind2,ind3,irhs))
               enddo
            enddo
         enddo

         do ind2 = ind2min, ind2max
            do ind3 = ind3min,ind3max
               RHS_x(ind2,ind3,1)   = Krho(ind3,ind2,ind1)
               RHS_x(ind2,ind3,2)   = Krhou(ind3,ind2,ind1)
               RHS_x(ind2,ind3,3)   = Krhov(ind3,ind2,ind1)
               RHS_x(ind2,ind3,4)   = Krhow(ind3,ind2,ind1)
               RHS_x(ind2,ind3,5)   = Krhoe(ind3,ind2,ind1)
            enddo
         enddo

         if (BC_face(1,1)%sort<=0) then
            ind3 = ind3min
            do ind2 = ind2min, ind2max
               ! Lower diagonal
               DLp_x(ind2,ind3,irhs) = 0.
               ! Upper diagonal
               DUp_x(ind2,ind3,irhs) = - theta_irs1 * deltat * idx(ind3)* (r_spec2(ind3+1,ind2,ind1))**0.5
               ! Main diagonal
               Dp_x(ind2,ind3,irhs)  = 1  - DUp_x(ind2,ind3,irhs)
            enddo
         endif
         if (BC_face(1,2)%sort<=0) then
            ind3 = ind3max
            do ind2 = ind2min, ind2max
               ! Lower diagonal
               DLp_x(ind2,ind3,irhs) = - theta_irs1 * deltat * idx(ind3)* (r_spec2(ind3,ind2,ind1))**0.5
               ! Upper diagonal
               DUp_x(ind2,ind3,irhs) = 0.
               ! Main diagonal
               Dp_x(ind2,ind3,irhs)  = 1 - DLp_x(ind2,ind3,irhs)
            enddo
         endif

         ! Resolution of the tridiagonal linear system
         ! ===========================================
         if (bl(1)%BC(1)==1) then
            do irhs=1,NRHS
               call PaScaL_TDMA_many_solve_cycle(p_many_x,DLp_x(:,:,irhs),Dp_x(:,:,irhs),DUp_x(:,:,irhs),RHS_x(:,:,irhs),ny,nx)
            enddo
         else
            do irhs=1,NRHS
               call PaScaL_TDMA_many_solve(p_many_x,DLp_x(:,:,irhs),Dp_x(:,:,irhs),DUp_x(:,:,irhs),RHS_x(:,:,irhs),ny,nx)
            enddo
         endif

         do ind2 = ind2min, ind2max
            do ind3 = ind3min,ind3max
               Krho(ind3,ind2,ind1)  = RHS_x(ind2,ind3,1)
               Krhou(ind3,ind2,ind1) = RHS_x(ind2,ind3,2)
               Krhov(ind3,ind2,ind1) = RHS_x(ind2,ind3,3)
               Krhow(ind3,ind2,ind1) = RHS_x(ind2,ind3,4)
               Krhoe(ind3,ind2,ind1) = RHS_x(ind2,ind3,5)
            enddo
         enddo
      enddo

   endif

   !----------------------------------------
   ! Implicitation of direction y
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

      ! Calculation of the spectral radius at ind-1/2
      ! =============================================
      if (.not. is_2d) then
         do ind1=ind1min,ind1max
            do ind2=ind2min,ind2max
               do ind3=ind3min,ind3max+1
                  rspecmh = sqrt(vv(ind2,ind3-1,ind1)**2) + c_(ind2,ind3-1,ind1)
                  rspecph = sqrt(vv(ind2,ind3,ind1)**2) + c_(ind2,ind3,ind1)
                  r_spec2(ind2,ind3,ind1) = (0.5 * (rspecmh + rspecph))**2
               enddo
            enddo
         enddo
      else
         do ind2=ind2min,ind2max
            do ind3=ind3min,ind3max+1
               rspecmh = sqrt(vv(ind2,ind3-1,1)**2) + c_(ind2,ind3-1,1)
               rspecph = sqrt(vv(ind2,ind3,1)**2) + c_(ind2,ind3,1)
               r_spec2(ind2,ind3,1) = (0.5 * (rspecmh + rspecph))**2
            enddo
         enddo
      endif

      do ind1 = ind1min, ind1max
         ! Preparation of the arrays and matrix
         ! ====================================
         do irhs = 1,NRHS
            do ind2 = ind2min, ind2max
               do ind3 = ind3min,ind3max
                  coef = - theta_irs2 * (deltat * idy(ind3))**2

                  ! Lower diagonal
                  DLp_y(ind2,ind3,irhs) = coef * r_spec2(ind2,ind3,ind1)
                  ! Upper diagonal
                  DUp_y(ind2,ind3,irhs) = coef * r_spec2(ind2,ind3+1,ind1)
                  ! Main diagonal
                  Dp_y(ind2,ind3,irhs)  = 1  - (DLp_y(ind2,ind3,irhs) + DUp_y(ind2,ind3,irhs))
               enddo
            enddo
         enddo

         if (BC_face(2,1)%sort<=0) then
            ind3 = ind3min
            do irhs=1,NRHS
               do ind2 = ind2min, ind2max
                  ! Lower diagonal
                  DLp_y(ind2,ind3,irhs) = 0.
                  ! Upper diagonal
                  DUp_y(ind2,ind3,irhs) = - theta_irs1 * deltat * idy(ind3)* (r_spec2(ind2,ind3+1,ind1))**0.5
                  ! Main diagonal
                  Dp_y(ind2,ind3,irhs)  = 1  - DUp_y(ind2,ind3,irhs)
               enddo
            enddo
         endif
         if (BC_face(2,2)%sort<=0) then
            ind3 = ind3max
            do irhs=1,NRHS
               do ind2 = ind2min, ind2max
                  ! Lower diagonal
                  DLp_y(ind2,ind3,irhs) = - theta_irs1 * deltat * idy(ind3)* (r_spec2(ind2,ind3,ind1))**0.5
                  ! Upper diagonal
                  DUp_y(ind2,ind3,irhs) = 0.
                  ! Main diagonal
                  Dp_y(ind2,ind3,irhs)  = 1 - DLp_y(ind2,ind3,irhs)
               enddo
            enddo
         endif

         ! Resolution of the tridiagonal linear system
         ! ===========================================
         if (bl(1)%BC(3)==1) then
            call PaScaL_TDMA_many_solve_cycle(p_many_y,DLp_y(:,:,1),Dp_y(:,:,1),DUp_y(:,:,1), Krho(1:nx,1:ny,ind1),nx,ny)
            call PaScaL_TDMA_many_solve_cycle(p_many_y,DLp_y(:,:,2),Dp_y(:,:,2),DUp_y(:,:,2),Krhou(1:nx,1:ny,ind1),nx,ny)
            call PaScaL_TDMA_many_solve_cycle(p_many_y,DLp_y(:,:,3),Dp_y(:,:,3),DUp_y(:,:,3),Krhov(1:nx,1:ny,ind1),nx,ny)
            call PaScaL_TDMA_many_solve_cycle(p_many_y,DLp_y(:,:,4),Dp_y(:,:,4),DUp_y(:,:,4),Krhow(1:nx,1:ny,ind1),nx,ny)
            call PaScaL_TDMA_many_solve_cycle(p_many_y,DLp_y(:,:,5),Dp_y(:,:,5),DUp_y(:,:,5),Krhoe(1:nx,1:ny,ind1),nx,ny)
         else
            call PaScaL_TDMA_many_solve(p_many_y,DLp_y(:,:,1),Dp_y(:,:,1),DUp_y(:,:,1), Krho(1:nx,1:ny,ind1),nx,ny)
            call PaScaL_TDMA_many_solve(p_many_y,DLp_y(:,:,2),Dp_y(:,:,2),DUp_y(:,:,2),Krhou(1:nx,1:ny,ind1),nx,ny)
            call PaScaL_TDMA_many_solve(p_many_y,DLp_y(:,:,3),Dp_y(:,:,3),DUp_y(:,:,3),Krhov(1:nx,1:ny,ind1),nx,ny)
            call PaScaL_TDMA_many_solve(p_many_y,DLp_y(:,:,4),Dp_y(:,:,4),DUp_y(:,:,4),Krhow(1:nx,1:ny,ind1),nx,ny)
            call PaScaL_TDMA_many_solve(p_many_y,DLp_y(:,:,5),Dp_y(:,:,5),DUp_y(:,:,5),Krhoe(1:nx,1:ny,ind1),nx,ny)
         endif
      enddo

   endif

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
         do irhs = 1,NRHS
            do ind3 = ind3min,ind3max
               do ind2 = ind2min,ind2max
                  coef = - theta_irs2 * (deltat * idz(ind3))**2

                  ! Lower diagonal
                  DLp_z(ind2,ind3,irhs) = coef * r_spec2(ind2,ind1,ind3)
                  ! Upper diagonal
                  DUp_z(ind2,ind3,irhs) = coef * r_spec2(ind2,ind1,ind3+1)
                  ! Main diagonal
                  Dp_z(ind2,ind3,irhs)  = 1  - (DLp_z(ind2,ind3,irhs) + DUp_z(ind2,ind3,irhs))
               enddo
            enddo
         enddo

         if (BC_face(3,1)%sort<=0) then
            ind3 = ind3min
            do irhs=1,NRHS
               do ind2 = ind2min, ind2max
                  ! Lower diagonal
                  DLp_z(ind2,ind3,irhs) = 0.
                  ! Upper diagonal
                  DUp_z(ind2,ind3,irhs) = - theta_irs1 * deltat * idz(ind3)* (r_spec2(ind2,ind1,ind3+1))**0.5
                  ! Main diagonal
                  Dp_z(ind2,ind3,irhs)  = 1  - DUp_z(ind2,ind3,irhs)
               enddo
            enddo
         endif
         if (BC_face(3,2)%sort<=0) then
            ind3 = ind3max
            do irhs=1,NRHS
               do ind2 = ind2min, ind2max
                  ! Lower diagonal
                  DLp_z(ind2,ind3,irhs) = - theta_irs1 * deltat * idz(ind3)* (r_spec2(ind2,ind1,ind3))**0.5
                  ! Upper diagonal
                  DUp_z(ind2,ind3,irhs) = 0.
                  ! Main diagonal
                  Dp_z(ind2,ind3,irhs)  = 1 - DLp_z(ind2,ind3,irhs)
               enddo
            enddo
         endif

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
         endif
      enddo

   endif

end subroutine irs_pascal_v0

! !===============================================================================
! subroutine ls_irs2_gh(array_a,array_c,array_d,array_y1,array_y2,array_y3,array_y4,array_y5,ind3min,ind3max)
! !===============================================================================
! !> Resolution of the tridiagonal system associated to IRS2 with Thomas' algo
! !===============================================================================
!    use precision
!    implicit none
!    real(wp),dimension(ind3min:ind3max)    :: array_d
!    real(wp),dimension(ind3min:ind3max-1)  :: array_a,array_c
!    real(wp),dimension(ind3min:ind3max)    :: array_y1, array_y2, array_y3, array_y4, array_y5
!    real(wp) :: temp
!    integer  :: ind3,ind3min,ind3max

!    ! Thomas' algorithm
!    ! Forward elimination phase
!    array_a(ind3min)  = array_a(ind3min) / array_d(ind3min)
!    array_y1(ind3min) = array_y1(ind3min) / array_d(ind3min)
!    array_y2(ind3min) = array_y2(ind3min) / array_d(ind3min)
!    array_y3(ind3min) = array_y3(ind3min) / array_d(ind3min)
!    array_y4(ind3min) = array_y4(ind3min) / array_d(ind3min)
!    array_y5(ind3min) = array_y5(ind3min) / array_d(ind3min)

!    do ind3 = ind3min+1, ind3max-1
!       temp = 1. / (array_d(ind3) - array_c(ind3-1)*array_a(ind3-1))
!       array_a(ind3)   = array_a(ind3) * temp
!       array_y1(ind3)   = (array_y1(ind3) - array_c(ind3-1)*array_y1(ind3-1)) * temp
!       array_y2(ind3)   = (array_y2(ind3) - array_c(ind3-1)*array_y2(ind3-1)) * temp
!       array_y3(ind3)   = (array_y3(ind3) - array_c(ind3-1)*array_y3(ind3-1)) * temp
!       array_y4(ind3)   = (array_y4(ind3) - array_c(ind3-1)*array_y4(ind3-1)) * temp
!       array_y5(ind3)   = (array_y5(ind3) - array_c(ind3-1)*array_y5(ind3-1)) * temp
!    enddo

!    ! Backward substitution phase
!    temp = 1. / (array_d(ind3max) - array_c(ind3max-1)*array_a(ind3max-1))
!    array_y1(ind3max) = (array_y1(ind3max) - array_c(ind3max-1)*array_y1(ind3max-1)) * temp
!    array_y2(ind3max) = (array_y2(ind3max) - array_c(ind3max-1)*array_y2(ind3max-1)) * temp
!    array_y3(ind3max) = (array_y3(ind3max) - array_c(ind3max-1)*array_y3(ind3max-1)) * temp
!    array_y4(ind3max) = (array_y4(ind3max) - array_c(ind3max-1)*array_y4(ind3max-1)) * temp
!    array_y5(ind3max) = (array_y5(ind3max) - array_c(ind3max-1)*array_y5(ind3max-1)) * temp

!    do ind3 = ind3max-1, ind3min, -1
!       array_y1(ind3) = array_y1(ind3) - array_a(ind3)*array_y1(ind3+1)
!       array_y2(ind3) = array_y2(ind3) - array_a(ind3)*array_y2(ind3+1)
!       array_y3(ind3) = array_y3(ind3) - array_a(ind3)*array_y3(ind3+1)
!       array_y4(ind3) = array_y4(ind3) - array_a(ind3)*array_y4(ind3+1)
!       array_y5(ind3) = array_y5(ind3) - array_a(ind3)*array_y5(ind3+1)
!    enddo

! end subroutine ls_irs2_gh


! !===============================================================================
! subroutine ls_irs4_gh(array_a,array_b,array_c,array_d,array_e,array_y1,array_y2,array_y3,array_y4,array_y5,ind3min,ind3max)
! !===============================================================================
! !> Resolution of the pentadiagonal linear system associated to IRS
! !===============================================================================
!    use precision
!    implicit none
!    real(wp),dimension(ind3min:ind3max)    :: array_d
!    real(wp),dimension(ind3min:ind3max-1)  :: array_a,array_c
!    real(wp),dimension(ind3min:ind3max-2)  :: array_b,array_e
!    real(wp),dimension(ind3min:ind3max)    :: array_y1, array_y2, array_y3, array_y4, array_y5
!    real(wp) :: gam_,imu
!    integer  :: ind3,ind3min,ind3max

!    ! Resolution of the pentadiagonal linear system
!    ! Algo from http://www.hindawi.com/journals/mpe/2015/232456/
!    ! Step 3 & 4: for i = 1,2
!    ind3 = ind3min
!    imu = 1 / array_d(ind3)
!    array_a(ind3) = array_a(ind3) * imu
!    array_b(ind3) = array_b(ind3) * imu
!    array_y1(ind3)    = array_y1(ind3) * imu
!    array_y2(ind3)    = array_y2(ind3) * imu
!    array_y3(ind3)    = array_y3(ind3) * imu
!    array_y4(ind3)    = array_y4(ind3) * imu
!    array_y5(ind3)    = array_y5(ind3) * imu

!    ind3 = ind3min + 1
!    gam_ = array_c(ind3-1)
!    imu = 1. / (array_d(ind3) - array_a(ind3-1)*gam_)
!    array_a(ind3) = (array_a(ind3) - array_b(ind3-1)*gam_) * imu
!    array_b(ind3) = array_b(ind3) * imu
!    array_y1(ind3) = (array_y1(ind3) - array_y1(ind3-1)*gam_) * imu
!    array_y2(ind3) = (array_y2(ind3) - array_y2(ind3-1)*gam_) * imu
!    array_y3(ind3) = (array_y3(ind3) - array_y3(ind3-1)*gam_) * imu
!    array_y4(ind3) = (array_y4(ind3) - array_y4(ind3-1)*gam_) * imu
!    array_y5(ind3) = (array_y5(ind3) - array_y5(ind3-1)*gam_) * imu

!    ! Step 5: for i = 3,...,n-2
!    do ind3 = ind3min+2, ind3max-2
!       gam_ = array_c(ind3-1) - array_a(ind3-2)*array_e(ind3-2)
!       imu = 1. / (array_d(ind3) - array_b(ind3-2)*array_e(ind3-2) - array_a(ind3-1)*gam_)
!       array_a(ind3) = (array_a(ind3) - array_b(ind3-1)*gam_) * imu
!       array_b(ind3) = array_b(ind3) * imu
!       array_y1(ind3) = (array_y1(ind3) - array_y1(ind3-2)*array_e(ind3-2) - array_y1(ind3-1)*gam_) * imu
!       array_y2(ind3) = (array_y2(ind3) - array_y2(ind3-2)*array_e(ind3-2) - array_y2(ind3-1)*gam_) * imu
!       array_y3(ind3) = (array_y3(ind3) - array_y3(ind3-2)*array_e(ind3-2) - array_y3(ind3-1)*gam_) * imu
!       array_y4(ind3) = (array_y4(ind3) - array_y4(ind3-2)*array_e(ind3-2) - array_y4(ind3-1)*gam_) * imu
!       array_y5(ind3) = (array_y5(ind3) - array_y5(ind3-2)*array_e(ind3-2) - array_y5(ind3-1)*gam_) * imu
!    enddo

!    ! Step 5: for i = n-1, n
!    ind3 = ind3max - 1
!    gam_ = array_c(ind3-1) - array_a(ind3-2)*array_e(ind3-2)
!    imu = 1. / (array_d(ind3) - array_b(ind3-2)*array_e(ind3-2) - array_a(ind3-1)*gam_)
!    array_a(ind3) = (array_a(ind3) - array_b(ind3-1)*gam_) * imu
!    array_y1(ind3) = (array_y1(ind3) - array_y1(ind3-2)*array_e(ind3-2) - array_y1(ind3-1)*gam_) * imu
!    array_y2(ind3) = (array_y2(ind3) - array_y2(ind3-2)*array_e(ind3-2) - array_y2(ind3-1)*gam_) * imu
!    array_y3(ind3) = (array_y3(ind3) - array_y3(ind3-2)*array_e(ind3-2) - array_y3(ind3-1)*gam_) * imu
!    array_y4(ind3) = (array_y4(ind3) - array_y4(ind3-2)*array_e(ind3-2) - array_y4(ind3-1)*gam_) * imu
!    array_y5(ind3) = (array_y5(ind3) - array_y5(ind3-2)*array_e(ind3-2) - array_y5(ind3-1)*gam_) * imu

!    ind3 = ind3max
!    gam_ = array_c(ind3-1) - array_a(ind3-2)*array_e(ind3-2)
!    imu = 1. / (array_d(ind3) - array_b(ind3-2)*array_e(ind3-2) - array_a(ind3-1)*gam_)
!    array_y1(ind3) = (array_y1(ind3) - array_y1(ind3-2)*array_e(ind3-2) - array_y1(ind3-1)*gam_) * imu
!    array_y2(ind3) = (array_y2(ind3) - array_y2(ind3-2)*array_e(ind3-2) - array_y2(ind3-1)*gam_) * imu
!    array_y3(ind3) = (array_y3(ind3) - array_y3(ind3-2)*array_e(ind3-2) - array_y3(ind3-1)*gam_) * imu
!    array_y4(ind3) = (array_y4(ind3) - array_y4(ind3-2)*array_e(ind3-2) - array_y4(ind3-1)*gam_) * imu
!    array_y5(ind3) = (array_y5(ind3) - array_y5(ind3-2)*array_e(ind3-2) - array_y5(ind3-1)*gam_) * imu


!    ! Step 6: computation of the solution vector
!    ind3 = ind3max - 1
!    array_y1(ind3) = array_y1(ind3) - array_a(ind3)*array_y1(ind3+1)
!    array_y2(ind3) = array_y2(ind3) - array_a(ind3)*array_y2(ind3+1)
!    array_y3(ind3) = array_y3(ind3) - array_a(ind3)*array_y3(ind3+1)
!    array_y4(ind3) = array_y4(ind3) - array_a(ind3)*array_y4(ind3+1)
!    array_y5(ind3) = array_y5(ind3) - array_a(ind3)*array_y5(ind3+1)

!    do ind3 = ind3max-2,ind3min,-1
!       array_y1(ind3) = array_y1(ind3) - array_a(ind3)*array_y1(ind3+1)  - array_b(ind3)*array_y1(ind3+2)
!       array_y2(ind3) = array_y2(ind3) - array_a(ind3)*array_y2(ind3+1)  - array_b(ind3)*array_y2(ind3+2)
!       array_y3(ind3) = array_y3(ind3) - array_a(ind3)*array_y3(ind3+1)  - array_b(ind3)*array_y3(ind3+2)
!       array_y4(ind3) = array_y4(ind3) - array_a(ind3)*array_y4(ind3+1)  - array_b(ind3)*array_y4(ind3+2)
!       array_y5(ind3) = array_y5(ind3) - array_a(ind3)*array_y5(ind3+1)  - array_b(ind3)*array_y5(ind3+2)
!    enddo
! end subroutine ls_irs4_gh


!===============================================================================
subroutine ls_irs6_gh(array_a,array_b,array_c,array_d,array_e,array_f,array_g,array_y1,array_y2,array_y3,array_y4,array_y5,ind3min,ind3max)
!===============================================================================
!> Resolution of the heptadiagonal linear system associated to IRS6
!===============================================================================
   implicit none
   double precision, dimension(ind3min:ind3max)    :: array_d
   double precision, dimension(ind3min:ind3max-1)  :: array_a,array_e
   double precision, dimension(ind3min:ind3max-2)  :: array_b,array_f
   double precision, dimension(ind3min:ind3max-3)  :: array_c,array_g
   double precision, dimension(ind3min:ind3max)    :: array_y1, array_y2, array_y3, array_y4, array_y5
   double precision :: temp
   integer  :: ind3,ind3min,ind3max

   ! Resolution of the pentadiagonal linear system
   ! Algo from https://www.researchgate.net/publication/47820097_A_New_Algorithm_for_General_Cyclic_
   !           Heptadiagonal_Linear_Systems_Using_Sherman-Morrison-Woodbury_formula
   ! Step 1 & step 5
   ind3 = ind3min
   temp = 1. / array_d(ind3)
   array_e(ind3) = array_e(ind3) * temp
   array_f(ind3) = array_f(ind3) * temp

   ind3 = ind3min + 1
   array_d(ind3)  = array_d(ind3) - array_e(ind3-1)*array_a(ind3-1)
   array_a(ind3)  = array_a(ind3) - array_e(ind3-1)*array_b(ind3-1)
   array_y1(ind3) = array_y1(ind3) - array_e(ind3-1)*array_y1(ind3-1)
   array_y2(ind3) = array_y2(ind3) - array_e(ind3-1)*array_y2(ind3-1)
   array_y3(ind3) = array_y3(ind3) - array_e(ind3-1)*array_y3(ind3-1)
   array_y4(ind3) = array_y4(ind3) - array_e(ind3-1)*array_y4(ind3-1)
   array_y5(ind3) = array_y5(ind3) - array_e(ind3-1)*array_y5(ind3-1)
   array_e(ind3)  = (array_e(ind3) - array_f(ind3-1)*array_a(ind3-1)) / array_d(ind3)

   ind3 = ind3min + 2
   array_d(ind3)  = array_d(ind3) - array_f(ind3-2)*array_b(ind3) - array_e(ind3-1)*array_a(ind3-1)
   array_y1(ind3)  = array_y1(ind3) - array_f(ind3-2)*array_y1(ind3-2) - array_e(ind3-1)*array_y1(ind3-1)
   array_y2(ind3)  = array_y2(ind3) - array_f(ind3-2)*array_y2(ind3-2) - array_e(ind3-1)*array_y2(ind3-1)
   array_y3(ind3)  = array_y3(ind3) - array_f(ind3-2)*array_y3(ind3-2) - array_e(ind3-1)*array_y3(ind3-1)
   array_y4(ind3)  = array_y4(ind3) - array_f(ind3-2)*array_y4(ind3-2) - array_e(ind3-1)*array_y4(ind3-1)
   array_y5(ind3)  = array_y5(ind3) - array_f(ind3-2)*array_y5(ind3-2) - array_e(ind3-1)*array_y5(ind3-1)

   ! Step 2 & step 5
   do ind3 = ind3min+3,ind3max
      temp = 1. / array_d(ind3-3)
      array_f(ind3-2) = (array_f(ind3-2) - array_g(ind3-3)*array_a(ind3-3)*temp) / array_d(ind3-2)
      array_e(ind3-1) = (array_e(ind3-1) - array_g(ind3-3)*array_b(ind3-3)*temp - array_f(ind3-2)*array_a(ind3-2)) / array_d(ind3-1)
      array_b(ind3-2) = array_b(ind3-2) - array_e(ind3-3)*array_c(ind3-3)
      array_a(ind3-1) = array_a(ind3-1) - array_e(ind3-2)*array_b(ind3-2) - array_f(ind3-3)*array_c(ind3-3)
      array_d(ind3)   = array_d(ind3) - array_g(ind3-3)*array_c(ind3-3)*temp - array_f(ind3-2)*array_b(ind3-2) - array_e(ind3-1)*array_a(ind3-1)
      array_y1(ind3)  = array_y1(ind3) - array_g(ind3-3)*array_y1(ind3-3)*temp - array_f(ind3-2)*array_y1(ind3-2) - array_e(ind3-1)*array_y1(ind3-1)
      array_y2(ind3)  = array_y2(ind3) - array_g(ind3-3)*array_y2(ind3-3)*temp - array_f(ind3-2)*array_y2(ind3-2) - array_e(ind3-1)*array_y2(ind3-1)
      array_y3(ind3)  = array_y3(ind3) - array_g(ind3-3)*array_y3(ind3-3)*temp - array_f(ind3-2)*array_y3(ind3-2) - array_e(ind3-1)*array_y3(ind3-1)
      array_y4(ind3)  = array_y4(ind3) - array_g(ind3-3)*array_y4(ind3-3)*temp - array_f(ind3-2)*array_y4(ind3-2) - array_e(ind3-1)*array_y4(ind3-1)
      array_y5(ind3)  = array_y5(ind3) - array_g(ind3-3)*array_y5(ind3-3)*temp - array_f(ind3-2)*array_y5(ind3-2) - array_e(ind3-1)*array_y5(ind3-1)
   enddo

   ! Step 6
   array_y1(ind3max)   = array_y1(ind3max) / array_d(ind3max)
   array_y2(ind3max)   = array_y2(ind3max) / array_d(ind3max)
   array_y3(ind3max)   = array_y3(ind3max) / array_d(ind3max)
   array_y4(ind3max)   = array_y4(ind3max) / array_d(ind3max)
   array_y5(ind3max)   = array_y5(ind3max) / array_d(ind3max)

   array_y1(ind3max-1) = (array_y1(ind3max-1) - array_a(ind3max-1)*array_y1(ind3max)) / array_d(ind3max-1)
   array_y2(ind3max-1) = (array_y2(ind3max-1) - array_a(ind3max-1)*array_y2(ind3max)) / array_d(ind3max-1)
   array_y3(ind3max-1) = (array_y3(ind3max-1) - array_a(ind3max-1)*array_y3(ind3max)) / array_d(ind3max-1)
   array_y4(ind3max-1) = (array_y4(ind3max-1) - array_a(ind3max-1)*array_y4(ind3max)) / array_d(ind3max-1)
   array_y5(ind3max-1) = (array_y5(ind3max-1) - array_a(ind3max-1)*array_y5(ind3max)) / array_d(ind3max-1)

   array_y1(ind3max-2) = (array_y1(ind3max-2) - array_a(ind3max-2)*array_y1(ind3max-1) - array_b(ind3max-2)*array_y1(ind3max)) / array_d(ind3max-2)
   array_y2(ind3max-2) = (array_y2(ind3max-2) - array_a(ind3max-2)*array_y2(ind3max-1) - array_b(ind3max-2)*array_y2(ind3max)) / array_d(ind3max-2)
   array_y3(ind3max-2) = (array_y3(ind3max-2) - array_a(ind3max-2)*array_y3(ind3max-1) - array_b(ind3max-2)*array_y3(ind3max)) / array_d(ind3max-2)
   array_y4(ind3max-2) = (array_y4(ind3max-2) - array_a(ind3max-2)*array_y4(ind3max-1) - array_b(ind3max-2)*array_y4(ind3max)) / array_d(ind3max-2)
   array_y5(ind3max-2) = (array_y5(ind3max-2) - array_a(ind3max-2)*array_y5(ind3max-1) - array_b(ind3max-2)*array_y5(ind3max)) / array_d(ind3max-2)

   do ind3 = ind3max-3,ind3min,-1
      array_y1(ind3) = (array_y1(ind3) - array_a(ind3)*array_y1(ind3+1) - array_b(ind3)*array_y1(ind3+2) - array_c(ind3)*array_y1(ind3+3)) / array_d(ind3)
      array_y2(ind3) = (array_y2(ind3) - array_a(ind3)*array_y2(ind3+1) - array_b(ind3)*array_y2(ind3+2) - array_c(ind3)*array_y2(ind3+3)) / array_d(ind3)
      array_y3(ind3) = (array_y3(ind3) - array_a(ind3)*array_y3(ind3+1) - array_b(ind3)*array_y3(ind3+2) - array_c(ind3)*array_y3(ind3+3)) / array_d(ind3)
      array_y4(ind3) = (array_y4(ind3) - array_a(ind3)*array_y4(ind3+1) - array_b(ind3)*array_y4(ind3+2) - array_c(ind3)*array_y4(ind3+3)) / array_d(ind3)
      array_y5(ind3) = (array_y5(ind3) - array_a(ind3)*array_y5(ind3+1) - array_b(ind3)*array_y5(ind3+2) - array_c(ind3)*array_y5(ind3+3)) / array_d(ind3)
   enddo

end subroutine ls_irs6_gh





