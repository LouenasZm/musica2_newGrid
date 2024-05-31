!===============================================================================
subroutine grid_comm_metrics_cart(dx_,nx_,dir,ngh_)
!===============================================================================
  !> Communicate Cartesian metrics to fill ghost pts
  !> (simple blocking sendrecv)
!===============================================================================
  use mod_mpi_part
  use mod_grid
  implicit none
  ! ----------------------------------------------------------------------------
  integer, intent(in) :: dir  ! direction
  integer, intent(in) :: nx_  ! size
  integer, intent(in) :: ngh_ ! number of ghost points
  real(wp), dimension(1-ngh_:nx_+ngh_), intent(inout) :: dx_ ! metrics
  ! ---------------------------------------------------------------------------
  integer :: n1,n2   ! neighbors
  integer :: type_dx ! MPI type
  ! ----------------------------------------------------------------------------
  
  ! Fill ghost points (ngh_)
  ! ========================
  
  ! Define neighbors depending on dimension direction
  ! -------------------------------------------------
  n1=2*dir
  n2=2*dir-1
  
  ! MPI type for 1D exchange of size ngh_
  ! -------------------------------------
  call MPI_TYPE_VECTOR(ngh_,1,1,MPI_DOUBLE_PRECISION,type_dx,info)
  call MPI_TYPE_COMMIT(type_dx,info)

  ! MPI SENDRECV
  ! ------------
    
  if(is_adjoint_block) then ! for adjoint block interfaces
    
    ! Send to neighbor 1 and reception from neighbor 2
    call MPI_SENDRECV(dx_(nx_-ngh   ),1,type_dx,neighbor(n1),tag &
                     ,dx_(   -ngh_+1),1,type_dx,neighbor(n2),tag,COMM_global,status,info)  
    
    ! Send to neighbor 2 and reception from neighbor 1
    call MPI_SENDRECV(dx_(2)    ,1,type_dx,neighbor(n2),tag &
                     ,dx_(nx_+1),1,type_dx,neighbor(n1),tag,COMM_global,status,info)
    
  else
    
    ! Send to neighbor 1 and reception from neighbor 2
    call MPI_SENDRECV(dx_(nx_-ngh_+1),1,type_dx,neighbor(n1),tag &
                     ,dx_(   -ngh_+1),1,type_dx,neighbor(n2),tag,COMM_global,status,info)  
    ! Send to neighbor 2 and reception from neighbor 1
    call MPI_SENDRECV(dx_(1)    ,1,type_dx,neighbor(n2),tag &
                     ,dx_(nx_+1),1,type_dx,neighbor(n1),tag,COMM_global,status,info)
                     
  endif
    
end subroutine grid_comm_metrics_cart

!===============================================================================
subroutine grid_comm_metrics_curv
!===============================================================================
  !> Communicate curvilinear metrics to fill ghost points
!===============================================================================
  use mod_mpi_types_two_sided
  use mod_grid
  implicit none
  ! ----------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------
  
  ! Communications of metrics for inviscid fluxes (ngh ghost points)
  ! =============================================
  call comm_metrics(x_ksi,x_eta)
  call comm_metrics(y_ksi,y_eta)
  call MPI_BARRIER(COMM_global,info)

  ! Communications of metrics for viscous fluxes (2*ngh_v ghost points)
  ! ============================================
  call comm_metrics_v(x_ksi_v,x_eta_v)
  call comm_metrics_v(y_ksi_v,y_eta_v)
  call MPI_BARRIER(COMM_global,info)

  ! Delete MPI types for metrics communications
  call free_mpi_types_metrics

contains

  !=========================================================================
  subroutine comm_metrics(var1,var2)
  !=========================================================================
    !> Routine for communications of curvilinear metrics [inviscid]
  !=========================================================================
    implicit none
    ! ----------------------------------------------------------------------
    real(wp), dimension(nx1:nx2,ny1:ny2), intent(inout) :: var1,var2
    ! ----------------------------------------------------------------------
    integer :: n
    integer, parameter :: size2d=2*4
    integer, dimension(size2d) :: request2d
    integer, dimension(MPI_STATUS_SIZE,size2d) :: status2d
    ! ----------------------------------------------------------------------

    n=1
    call MPI_ISEND(var1(iis(n),ijs(n)),1,type_met(n), &
                   neighbor(n),tags(n),COMM_global,request2d(1),info)
    call MPI_IRECV(var1(iir(n),ijr(n)),1,type_met(n), &
                   neighbor(n),tagr(n),COMM_global,request2d(2),info)

    n=2
    call MPI_ISEND(var1(iis(n),ijs(n)),1,type_met(n), &
                   neighbor(n),tags(n),COMM_global,request2d(3),info)
    call MPI_IRECV(var1(iir(n),ijr(n)),1,type_met(n), &
                   neighbor(n),tagr(n),COMM_global,request2d(4),info)

    n=3
    call MPI_ISEND(var2(iis(n),ijs(n)),1,type_met(n), &
                   neighbor(n),tags(n),COMM_global,request2d(5),info)
    call MPI_IRECV(var2(iir(n),ijr(n)),1,type_met(n), &
                   neighbor(n),tagr(n),COMM_global,request2d(6),info)

    n=4
    call MPI_ISEND(var2(iis(n),ijs(n)),1,type_met(n), &
                   neighbor(n),tags(n),COMM_global,request2d(7),info)
    call MPI_IRECV(var2(iir(n),ijr(n)),1,type_met(n), &
                   neighbor(n),tagr(n),COMM_global,request2d(8),info)

    call MPI_WAITALL(size2d,request2d,status2d,info)
    
  end subroutine comm_metrics
  !=========================================================================

  !=========================================================================
  subroutine comm_metrics_v(var1,var2)
  !=========================================================================
    !> Routine for communications of curvilinear metrics [viscous]
  !=========================================================================
    implicit none
    ! ----------------------------------------------------------------------
    real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v), intent(inout) :: var1,var2
    ! ----------------------------------------------------------------------
    integer :: n
    integer, parameter :: size2d=2*4
    integer, dimension(size2d) :: request2d
    integer, dimension(MPI_STATUS_SIZE,size2d) :: status2d
    ! ----------------------------------------------------------------------

    n=1
    call MPI_ISEND(var1(iis_v(n),ijs_v(n)),1,type_met_v(n), &
                   neighbor(n),tags(n),COMM_global,request2d(1),info)
    call MPI_IRECV(var1(iir_v(n),ijr_v(n)),1,type_met_v(n), &
                   neighbor(n),tagr(n),COMM_global,request2d(2),info)

    n=2
    call MPI_ISEND(var1(iis_v(n),ijs_v(n)),1,type_met_v(n), &
                   neighbor(n),tags(n),COMM_global,request2d(3),info)
    call MPI_IRECV(var1(iir_v(n),ijr_v(n)),1,type_met_v(n), &
                   neighbor(n),tagr(n),COMM_global,request2d(4),info)

    n=3
    call MPI_ISEND(var2(iis_v(n),ijs_v(n)),1,type_met_v(n), &
                   neighbor(n),tags(n),COMM_global,request2d(5),info)
    call MPI_IRECV(var2(iir_v(n),ijr_v(n)),1,type_met_v(n), &
                   neighbor(n),tagr(n),COMM_global,request2d(6),info)

    n=4
    call MPI_ISEND(var2(iis_v(n),ijs_v(n)),1,type_met_v(n), &
                   neighbor(n),tags(n),COMM_global,request2d(7),info)
    call MPI_IRECV(var2(iir_v(n),ijr_v(n)),1,type_met_v(n), &
                   neighbor(n),tagr(n),COMM_global,request2d(8),info)

    call MPI_WAITALL(size2d,request2d,status2d,info)

  end subroutine comm_metrics_v
  !=========================================================================

end subroutine grid_comm_metrics_curv

!==============================================================================
subroutine comm_metrics_3d(var1,var2,var3)
!==============================================================================
  !> Routine for communications of 3D curvilinear metrics
!==============================================================================
  use mod_mpi_types_two_sided
  use mod_grid
  implicit none
  ! ---------------------------------------------------------------------------
  real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(inout) :: var1,var2,var3
  ! ---------------------------------------------------------------------------
  integer :: n
  integer, parameter :: size3d=2*6
  integer, dimension(size3d) :: request3d
  integer, dimension(MPI_STATUS_SIZE,size3d) :: status3d
  ! ---------------------------------------------------------------------------

  n=1
  call MPI_ISEND(var1(iis(n),ijs(n),iks(n)),1,type_face(n), &
                 neighbor(n),tags(n),COMM_global,request3d(1),info)
  call MPI_IRECV(var1(iir(n),ijr(n),ikr(n)),1,type_face(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(2),info)

  n=2
  call MPI_ISEND(var1(iis(n),ijs(n),iks(n)),1,type_face(n), &
                 neighbor(n),tags(n),COMM_global,request3d(3),info)
  call MPI_IRECV(var1(iir(n),ijr(n),ikr(n)),1,type_face(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(4),info)

  n=3
  call MPI_ISEND(var2(iis(n),ijs(n),iks(n)),1,type_face(n), &
                 neighbor(n),tags(n),COMM_global,request3d(5),info)
  call MPI_IRECV(var2(iir(n),ijr(n),ikr(n)),1,type_face(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(6),info)

  n=4
  call MPI_ISEND(var2(iis(n),ijs(n),iks(n)),1,type_face(n), &
                 neighbor(n),tags(n),COMM_global,request3d(7),info)
  call MPI_IRECV(var2(iir(n),ijr(n),ikr(n)),1,type_face(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(8),info)

  n=5
  call MPI_ISEND(var3(iis(n),ijs(n),iks(n)),1,type_face(n), &
                 neighbor(n),tags(n),COMM_global,request3d(9),info)
  call MPI_IRECV(var3(iir(n),ijr(n),ikr(n)),1,type_face(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(10),info)

  n=6
  call MPI_ISEND(var3(iis(n),ijs(n),iks(n)),1,type_face(n), &
                 neighbor(n),tags(n),COMM_global,request3d(11),info)
  call MPI_IRECV(var3(iir(n),ijr(n),ikr(n)),1,type_face(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(12),info)

  call MPI_WAITALL(size3d,request3d,status3d,info)

end subroutine comm_metrics_3d
!==============================================================================

!==============================================================================
subroutine comm_metrics_3de(ksi_1,eta_1,phi_1,ksi_2,eta_2,phi_2)
!==============================================================================
  !> Routine for communications of 3D curvilinear metrics
!==============================================================================
  use mod_mpi_types_two_sided
  use mod_grid
  implicit none
  ! ---------------------------------------------------------------------------
  real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(inout) :: ksi_1,eta_1,phi_1
  real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(inout) :: ksi_2,eta_2,phi_2
  ! ---------------------------------------------------------------------------
  integer :: n
  integer, parameter :: size3d=2*6*2
  integer, dimension(size3d) :: request3d
  integer, dimension(MPI_STATUS_SIZE,size3d) :: status3d
  ! ---------------------------------------------------------------------------

  n=1
  if ((is_swapij2(n)).or.(is_swapik2(n))) then
     call MPI_ISEND(phi_2(iis(n),ijs(n),iks(n)),1,type_face(n), &
                    neighbor(n),tags(n),COMM_global,request3d(1),info)
  else
     call MPI_ISEND(eta_1(iis(n),ijs(n),iks(n)),1,type_face(n), &
                    neighbor(n),tags(n),COMM_global,request3d(1),info)
  endif
  call MPI_IRECV(eta_1(iir(n),ijr(n),ikr(n)),1,type_face(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(2),info)

  n=2
  if ((is_swapij2(n)).or.(is_swapik2(n))) then
     call MPI_ISEND(phi_2(iis(n),ijs(n),iks(n)),1,type_face(n), &
                    neighbor(n),tags(n),COMM_global,request3d(3),info)
  else
     call MPI_ISEND(eta_1(iis(n),ijs(n),iks(n)),1,type_face(n), &
                    neighbor(n),tags(n),COMM_global,request3d(3),info)
  endif
  call MPI_IRECV(eta_1(iir(n),ijr(n),ikr(n)),1,type_face(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(4),info)

  n=3
  if ((is_swapij2(n)).or.(is_swapjk2(n))) then
     call MPI_ISEND(ksi_2(iis(n),ijs(n),iks(n)),1,type_face(n), &
                    neighbor(n),tags(n),COMM_global,request3d(5),info)
  else
     call MPI_ISEND(phi_1(iis(n),ijs(n),iks(n)),1,type_face(n), &
                    neighbor(n),tags(n),COMM_global,request3d(5),info)
  endif
  call MPI_IRECV(phi_1(iir(n),ijr(n),ikr(n)),1,type_face(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(6),info)

  n=4
  if ((is_swapij2(n)).or.(is_swapjk2(n))) then
     call MPI_ISEND(ksi_2(iis(n),ijs(n),iks(n)),1,type_face(n), &
                    neighbor(n),tags(n),COMM_global,request3d(7),info)
  else
     call MPI_ISEND(phi_1(iis(n),ijs(n),iks(n)),1,type_face(n), &
                    neighbor(n),tags(n),COMM_global,request3d(7),info)
  endif
  call MPI_IRECV(phi_1(iir(n),ijr(n),ikr(n)),1,type_face(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(8),info)

  n=5
  if ((is_swapik2(n)).or.(is_swapjk2(n))) then
     call MPI_ISEND(eta_2(iis(n),ijs(n),iks(n)),1,type_face(n), &
                    neighbor(n),tags(n),COMM_global,request3d(9),info)
  else
     call MPI_ISEND(ksi_1(iis(n),ijs(n),iks(n)),1,type_face(n), &
                    neighbor(n),tags(n),COMM_global,request3d(9),info)
  endif
  call MPI_IRECV(ksi_1(iir(n),ijr(n),ikr(n)),1,type_face(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(10),info)

  n=6
  if ((is_swapik2(n)).or.(is_swapjk2(n))) then
     call MPI_ISEND(eta_2(iis(n),ijs(n),iks(n)),1,type_face(n), &
                    neighbor(n),tags(n),COMM_global,request3d(11),info)
  else
     call MPI_ISEND(ksi_1(iis(n),ijs(n),iks(n)),1,type_face(n), &
                    neighbor(n),tags(n),COMM_global,request3d(11),info)
  endif
  call MPI_IRECV(ksi_1(iir(n),ijr(n),ikr(n)),1,type_face(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(12),info)

  n=1
  if ((is_swapij2(n)).or.(is_swapik2(n))) then
     call MPI_ISEND(eta_1(iis(n),ijs(n),iks(n)),1,type_face(n), &
                    neighbor(n),tags(n),COMM_global,request3d(13),info)
  else
     call MPI_ISEND(phi_2(iis(n),ijs(n),iks(n)),1,type_face(n), &
                    neighbor(n),tags(n),COMM_global,request3d(13),info)
  endif
  call MPI_IRECV(phi_2(iir(n),ijr(n),ikr(n)),1,type_face(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(14),info)

  n=2
  if ((is_swapij2(n)).or.(is_swapik2(n))) then
     call MPI_ISEND(eta_1(iis(n),ijs(n),iks(n)),1,type_face(n), &
                    neighbor(n),tags(n),COMM_global,request3d(15),info)
  else
     call MPI_ISEND(phi_2(iis(n),ijs(n),iks(n)),1,type_face(n), &
                    neighbor(n),tags(n),COMM_global,request3d(15),info)
  endif
  call MPI_IRECV(phi_2(iir(n),ijr(n),ikr(n)),1,type_face(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(16),info)

  n=3
  if ((is_swapij2(n)).or.(is_swapjk2(n))) then
     call MPI_ISEND(phi_1(iis(n),ijs(n),iks(n)),1,type_face(n), &
                    neighbor(n),tags(n),COMM_global,request3d(17),info)
  else
     call MPI_ISEND(ksi_2(iis(n),ijs(n),iks(n)),1,type_face(n), &
                    neighbor(n),tags(n),COMM_global,request3d(17),info)
  endif
  call MPI_IRECV(ksi_2(iir(n),ijr(n),ikr(n)),1,type_face(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(18),info)

  n=4
  if ((is_swapij2(n)).or.(is_swapjk2(n))) then
     call MPI_ISEND(phi_1(iis(n),ijs(n),iks(n)),1,type_face(n), &
                    neighbor(n),tags(n),COMM_global,request3d(19),info)
  else
     call MPI_ISEND(ksi_2(iis(n),ijs(n),iks(n)),1,type_face(n), &
                    neighbor(n),tags(n),COMM_global,request3d(19),info)
  endif
  call MPI_IRECV(ksi_2(iir(n),ijr(n),ikr(n)),1,type_face(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(20),info)

  n=5
  if ((is_swapik2(n)).or.(is_swapjk2(n))) then
     call MPI_ISEND(ksi_1(iis(n),ijs(n),iks(n)),1,type_face(n), &
                    neighbor(n),tags(n),COMM_global,request3d(21),info)
  else
     call MPI_ISEND(eta_2(iis(n),ijs(n),iks(n)),1,type_face(n), &
                    neighbor(n),tags(n),COMM_global,request3d(21),info)
  endif
  call MPI_IRECV(eta_2(iir(n),ijr(n),ikr(n)),1,type_face(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(22),info)

  n=6
  if ((is_swapik2(n)).or.(is_swapjk2(n))) then
     call MPI_ISEND(ksi_1(iis(n),ijs(n),iks(n)),1,type_face(n), &
                    neighbor(n),tags(n),COMM_global,request3d(23),info)
  else
     call MPI_ISEND(eta_2(iis(n),ijs(n),iks(n)),1,type_face(n), &
                    neighbor(n),tags(n),COMM_global,request3d(23),info)
  endif
  call MPI_IRECV(eta_2(iir(n),ijr(n),ikr(n)),1,type_face(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(24),info)

  call MPI_WAITALL(size3d,request3d,status3d,info)

end subroutine comm_metrics_3de
!==============================================================================

!==============================================================================
subroutine comm_metrics_3d_v(var1,var2,var3)
!==============================================================================
  !> Routine for communications of curvilinear metrics
!==============================================================================
  use mod_mpi_types_two_sided
  use mod_grid
  implicit none
  ! ---------------------------------------------------------------------------
  real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v), intent(inout) :: var1,var2,var3
  ! ---------------------------------------------------------------------------
  integer :: n
  integer, parameter :: size3d=2*6
  integer, dimension(size3d) :: request3d
  integer, dimension(MPI_STATUS_SIZE,size3d) :: status3d
  ! ---------------------------------------------------------------------------

  n=1
  call MPI_ISEND(var1(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                 neighbor(n),tags(n),COMM_global,request3d(1),info)
  call MPI_IRECV(var1(iir_v(n),ijr_v(n),ikr_v(n)),1,type_face_v(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(2),info)

  n=2
  call MPI_ISEND(var1(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                 neighbor(n),tags(n),COMM_global,request3d(3),info)
  call MPI_IRECV(var1(iir_v(n),ijr_v(n),ikr_v(n)),1,type_face_v(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(4),info)

  n=3
  call MPI_ISEND(var2(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                 neighbor(n),tags(n),COMM_global,request3d(5),info)
  call MPI_IRECV(var2(iir_v(n),ijr_v(n),ikr_v(n)),1,type_face_v(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(6),info)

  n=4
  call MPI_ISEND(var2(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                 neighbor(n),tags(n),COMM_global,request3d(7),info)
  call MPI_IRECV(var2(iir_v(n),ijr_v(n),ikr_v(n)),1,type_face_v(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(8),info)

  n=5
  call MPI_ISEND(var3(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                 neighbor(n),tags(n),COMM_global,request3d(9),info)
  call MPI_IRECV(var3(iir_v(n),ijr_v(n),ikr_v(n)),1,type_face_v(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(10),info)

  n=6
  call MPI_ISEND(var3(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                 neighbor(n),tags(n),COMM_global,request3d(11),info)
  call MPI_IRECV(var3(iir_v(n),ijr_v(n),ikr_v(n)),1,type_face_v(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(12),info)

  call MPI_WAITALL(size3d,request3d,status3d,info)

end subroutine comm_metrics_3d_v

!==============================================================================
subroutine comm_metrics_3de_v(ksi_1,eta_1,phi_1,ksi_2,eta_2,phi_2)
!==============================================================================
  !> Routine for communications of 3D curvilinear metrics along eta
!==============================================================================
  use mod_mpi_types_two_sided
  use mod_grid
  implicit none
  ! ---------------------------------------------------------------------------
  real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v), intent(inout) :: ksi_1,eta_1,phi_1
  real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v), intent(inout) :: ksi_2,eta_2,phi_2
  ! ---------------------------------------------------------------------------
  integer :: n
  integer, parameter :: size3d=2*6*2
  integer, dimension(size3d) :: request3d
  integer, dimension(MPI_STATUS_SIZE,size3d) :: status3d
  ! ---------------------------------------------------------------------------

  n=1
  if ((is_swapij2(n)).or.(is_swapik2(n))) then
     call MPI_ISEND(phi_2(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                    neighbor(n),tags(n),COMM_global,request3d(1),info)
  else
     call MPI_ISEND(eta_1(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                    neighbor(n),tags(n),COMM_global,request3d(1),info)
  endif
  call MPI_IRECV(eta_1(iir_v(n),ijr_v(n),ikr_v(n)),1,type_face_v(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(2),info)

  n=2
  if ((is_swapij2(n)).or.(is_swapik2(n))) then
     call MPI_ISEND(phi_2(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                    neighbor(n),tags(n),COMM_global,request3d(3),info)
  else
     call MPI_ISEND(eta_1(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                    neighbor(n),tags(n),COMM_global,request3d(3),info)
  endif
  call MPI_IRECV(eta_1(iir_v(n),ijr_v(n),ikr_v(n)),1,type_face_v(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(4),info)

  n=3
  if ((is_swapij2(n)).or.(is_swapjk2(n))) then
     call MPI_ISEND(ksi_2(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                    neighbor(n),tags(n),COMM_global,request3d(5),info)
  else
     call MPI_ISEND(phi_1(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                    neighbor(n),tags(n),COMM_global,request3d(5),info)
  endif
  call MPI_IRECV(phi_1(iir_v(n),ijr_v(n),ikr_v(n)),1,type_face_v(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(6),info)

  n=4
  if ((is_swapij2(n)).or.(is_swapjk2(n))) then
     call MPI_ISEND(ksi_2(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                    neighbor(n),tags(n),COMM_global,request3d(7),info)
  else
     call MPI_ISEND(phi_1(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                    neighbor(n),tags(n),COMM_global,request3d(7),info)
  endif
  call MPI_IRECV(phi_1(iir_v(n),ijr_v(n),ikr_v(n)),1,type_face_v(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(8),info)

  n=5
  if ((is_swapik2(n)).or.(is_swapjk2(n))) then
     call MPI_ISEND(eta_2(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                    neighbor(n),tags(n),COMM_global,request3d(9),info)
  else
     call MPI_ISEND(ksi_1(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                    neighbor(n),tags(n),COMM_global,request3d(9),info)
  endif
  call MPI_IRECV(ksi_1(iir_v(n),ijr_v(n),ikr_v(n)),1,type_face_v(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(10),info)

  n=6
  if ((is_swapik2(n)).or.(is_swapjk2(n))) then
     call MPI_ISEND(eta_2(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                    neighbor(n),tags(n),COMM_global,request3d(11),info)
  else
     call MPI_ISEND(ksi_1(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                    neighbor(n),tags(n),COMM_global,request3d(11),info)
  endif
  call MPI_IRECV(ksi_1(iir_v(n),ijr_v(n),ikr_v(n)),1,type_face_v(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(12),info)

  n=1
  if ((is_swapij2(n)).or.(is_swapik2(n))) then
     call MPI_ISEND(eta_1(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                    neighbor(n),tags(n),COMM_global,request3d(13),info)
  else
     call MPI_ISEND(phi_2(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                    neighbor(n),tags(n),COMM_global,request3d(13),info)
  endif
  call MPI_IRECV(phi_2(iir_v(n),ijr_v(n),ikr_v(n)),1,type_face_v(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(14),info)

  n=2
  if ((is_swapij2(n)).or.(is_swapik2(n))) then
     call MPI_ISEND(eta_1(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                    neighbor(n),tags(n),COMM_global,request3d(15),info)
  else
     call MPI_ISEND(phi_2(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                    neighbor(n),tags(n),COMM_global,request3d(15),info)
  endif
  call MPI_IRECV(phi_2(iir_v(n),ijr_v(n),ikr_v(n)),1,type_face_v(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(16),info)

  n=3
  if ((is_swapij2(n)).or.(is_swapjk2(n))) then
     call MPI_ISEND(phi_1(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                    neighbor(n),tags(n),COMM_global,request3d(17),info)
  else
     call MPI_ISEND(ksi_2(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                    neighbor(n),tags(n),COMM_global,request3d(17),info)
  endif
  call MPI_IRECV(ksi_2(iir_v(n),ijr_v(n),ikr_v(n)),1,type_face_v(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(18),info)

  n=4
  if ((is_swapij2(n)).or.(is_swapjk2(n))) then
     call MPI_ISEND(phi_1(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                    neighbor(n),tags(n),COMM_global,request3d(19),info)
  else
     call MPI_ISEND(ksi_2(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                    neighbor(n),tags(n),COMM_global,request3d(19),info)
  endif
  call MPI_IRECV(ksi_2(iir_v(n),ijr_v(n),ikr_v(n)),1,type_face_v(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(20),info)

  n=5
  if ((is_swapik2(n)).or.(is_swapjk2(n))) then
     call MPI_ISEND(ksi_1(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                    neighbor(n),tags(n),COMM_global,request3d(21),info)
  else
     call MPI_ISEND(eta_2(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                    neighbor(n),tags(n),COMM_global,request3d(21),info)
  endif
  call MPI_IRECV(eta_2(iir_v(n),ijr_v(n),ikr_v(n)),1,type_face_v(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(22),info)

  n=6
  if ((is_swapik2(n)).or.(is_swapjk2(n))) then
     call MPI_ISEND(ksi_1(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                    neighbor(n),tags(n),COMM_global,request3d(23),info)
  else
     call MPI_ISEND(eta_2(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                    neighbor(n),tags(n),COMM_global,request3d(23),info)
  endif
  call MPI_IRECV(eta_2(iir_v(n),ijr_v(n),ikr_v(n)),1,type_face_v(n), &
                 neighbor(n),tagr(n),COMM_global,request3d(24),info)

  call MPI_WAITALL(size3d,request3d,status3d,info)

end subroutine comm_metrics_3de_v
!==============================================================================

!===============================================================================
subroutine grid_comm_metrics_curv_old
!===============================================================================
  !> Communicate curvilinear metrics to fill ghost points ** OLD VERSION **
!===============================================================================
  use mod_mpi_types_two_sided_old
  use mod_grid
  implicit none
  ! ----------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------

  ! Communications of metrics for inviscid fluxes (ngh ghost points)
  ! =============================================
  call comm_metrics_old(x_ksi,x_eta)
  call comm_metrics_old(y_ksi,y_eta)
  call MPI_BARRIER(COMM_global,info)

  ! Communications of metrics for viscous fluxes (2*ngh_v ghost points)
  ! ============================================
  call comm_metrics_v_old(x_ksi_v,x_eta_v)
  call comm_metrics_v_old(y_ksi_v,y_eta_v)
  call MPI_BARRIER(COMM_global,info)

  ! Delete MPI types for metrics communications
  call free_mpi_types_metrics_old

contains

  !=========================================================================
  subroutine comm_metrics_old(var1,var2)
  !=========================================================================
    !> Routine for communications of curvilinear metrics [inviscid] ** OLD VERSION **
  !=========================================================================
    implicit none
    ! ----------------------------------------------------------------------
    real(wp), dimension(nx1:nx2,ny1:ny2), intent(inout) :: var1,var2
    ! ----------------------------------------------------------------------
    integer, parameter :: size2d=2*4
    integer, dimension(size2d) :: request2d
    integer, dimension(MPI_STATUS_SIZE,size2d) :: status2d
    ! ----------------------------------------------------------------------

    if (is_adjoint_block) then ! for adjoint block interfaces
    
      ! Send to neighbor W and reception from neighbor W
      call MPI_ISEND(var1(inWs,ipW),1,type_mW,neighbor(nW),tags(nW),COMM_global,request2d(1),info)
      call MPI_IRECV(var1(inWr,ipW),1,type_mW,neighbor(nW),tagr(nW),COMM_global,request2d(2),info)
      
      ! Send to neighbor E and reception from neighbor E
      call MPI_ISEND(var1(inEs,ipE),1,type_mE,neighbor(nE),tags(nE),COMM_global,request2d(3),info)
      call MPI_IRECV(var1(inEr,ipE),1,type_mE,neighbor(nE),tagr(nE),COMM_global,request2d(4),info)
      
      ! Send to neighbor S and reception from neighbor S
      call MPI_ISEND(var2(ipS,inSs),1,type_mS,neighbor(nS),tags(nS),COMM_global,request2d(5),info)
      call MPI_IRECV(var2(ipS,inSr),1,type_mS,neighbor(nS),tagr(nS),COMM_global,request2d(6),info)
      
      ! Send to neighbor N and reception from neighbor N
      call MPI_ISEND(var2(ipN,inNs),1,type_mN,neighbor(nN),tags(nN),COMM_global,request2d(7),info)
      call MPI_IRECV(var2(ipN,inNr),1,type_mN,neighbor(nN),tagr(nN),COMM_global,request2d(8),info)
      
    else
      
      ! Send to neighbor W and reception from neighbor W
      call MPI_ISEND(var1(inW    ,ipW),1,type_mW,neighbor(nW),tags(nW),COMM_global,request2d(1),info)
      call MPI_IRECV(var1(inW-ngh,ipW),1,type_mW,neighbor(nW),tagr(nW),COMM_global,request2d(2),info)
   
      ! Send to neighbor E and reception from neighbor E
      call MPI_ISEND(var1(inE-ngh,ipE),1,type_mE,neighbor(nE),tags(nE),COMM_global,request2d(3),info)
      call MPI_IRECV(var1(inE    ,ipE),1,type_mE,neighbor(nE),tagr(nE),COMM_global,request2d(4),info)
   
      ! Send to neighbor S and reception from neighbor S
      call MPI_ISEND(var2(ipS,inS    ),1,type_mS,neighbor(nS),tags(nS),COMM_global,request2d(5),info)
      call MPI_IRECV(var2(ipS,inS-ngh),1,type_mS,neighbor(nS),tagr(nS),COMM_global,request2d(6),info)
   
      ! Send to neighbor N and reception from neighbor N
      call MPI_ISEND(var2(ipN,inN-ngh),1,type_mN,neighbor(nN),tags(nN),COMM_global,request2d(7),info)
      call MPI_IRECV(var2(ipN,inN    ),1,type_mN,neighbor(nN),tagr(nN),COMM_global,request2d(8),info)
    
    endif

    call MPI_WAITALL(size2d,request2d,status2d,info)
    
  end subroutine comm_metrics_old
  !=========================================================================

  !=========================================================================
  subroutine comm_metrics_v_old(var1,var2)
  !=========================================================================
    !> Routine for communications of curvilinear metrics [viscous] ** OLD VERSION **
  !=========================================================================
    implicit none
    ! ----------------------------------------------------------------------
    real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v), intent(inout) :: var1,var2
    ! ----------------------------------------------------------------------
    integer, parameter :: size2d=2*4
    integer, dimension(size2d) :: request2d
    integer, dimension(MPI_STATUS_SIZE,size2d) :: status2d
    ! ----------------------------------------------------------------------
    
    if (is_adjoint_block) then ! for adjoint block interfaces
    
      ! Send to neighbor W and reception from neighbor W
      call MPI_ISEND(var1(inWs_v,ipW),1,type_mW_v,neighbor(nW),tags(nW),COMM_global,request2d(1),info)
      call MPI_IRECV(var1(inWr_v,ipW),1,type_mW_v,neighbor(nW),tagr(nW),COMM_global,request2d(2),info)
      
      ! Send to neighbor E and reception from neighbor E
      call MPI_ISEND(var1(inEs_v,ipE),1,type_mE_v,neighbor(nE),tags(nE),COMM_global,request2d(3),info)
      call MPI_IRECV(var1(inEr_v,ipE),1,type_mE_v,neighbor(nE),tagr(nE),COMM_global,request2d(4),info)
      
      ! Send to neighbor S and reception from neighbor S
      call MPI_ISEND(var2(ipS,inSs_v),1,type_mS_v,neighbor(nS),tags(nS),COMM_global,request2d(5),info)
      call MPI_IRECV(var2(ipS,inSr_v),1,type_mS_v,neighbor(nS),tagr(nS),COMM_global,request2d(6),info)
      
      ! Send to neighbor N and reception from neighbor N
      call MPI_ISEND(var2(ipN,inNs_v),1,type_mN_v,neighbor(nN),tags(nN),COMM_global,request2d(7),info)
      call MPI_IRECV(var2(ipN,inNr_v),1,type_mN_v,neighbor(nN),tagr(nN),COMM_global,request2d(8),info)
      
    else
      
      ! Send to neighbor W and reception from neighbor W
      call MPI_ISEND(var1(inW_v      ,ipW),1,type_mW_v,neighbor(nW),tags(nW),COMM_global,request2d(1),info)
      call MPI_IRECV(var1(inW_v-ngh_v,ipW),1,type_mW_v,neighbor(nW),tagr(nW),COMM_global,request2d(2),info)
    
      ! Send to neighbor E and reception from neighbor E
      call MPI_ISEND(var1(inE_v-ngh_v,ipE),1,type_mE_v,neighbor(nE),tags(nE),COMM_global,request2d(3),info)
      call MPI_IRECV(var1(inE_v      ,ipE),1,type_mE_v,neighbor(nE),tagr(nE),COMM_global,request2d(4),info)
    
      ! Send to neighbor S and reception from neighbor S
      call MPI_ISEND(var2(ipS,inS_v      ),1,type_mS_v,neighbor(nS),tags(nS),COMM_global,request2d(5),info)
      call MPI_IRECV(var2(ipS,inS_v-ngh_v),1,type_mS_v,neighbor(nS),tagr(nS),COMM_global,request2d(6),info)
    
      ! Send to neighbor N and reception from neighbor N
      call MPI_ISEND(var2(ipN,inN_v-ngh_v),1,type_mN_v,neighbor(nN),tags(nN),COMM_global,request2d(7),info)
      call MPI_IRECV(var2(ipN,inN_v      ),1,type_mN_v,neighbor(nN),tagr(nN),COMM_global,request2d(8),info)
    
    endif

    call MPI_WAITALL(size2d,request2d,status2d,info)
    
  end subroutine comm_metrics_v_old
  !=========================================================================

end subroutine grid_comm_metrics_curv_old

