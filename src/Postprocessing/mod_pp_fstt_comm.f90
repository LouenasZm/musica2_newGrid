!==============================================================================
module mod_pp_fstt_comm
!==============================================================================
  !> Module for communication in FSTT PP
!==============================================================================
  use warnstop
  use mod_pp_var
  use mod_mpi_part
  use mod_grid_directions
  implicit none
  ! ---------------------------------------------------------------------------
  ! 3D version, for 3 var
  integer, parameter :: size_3d_3v=12*3 ! 2(send/recv)*6(neighbors)*3(variables)
  integer, dimension(size_3d_3v) :: request_3d_3v
  integer, dimension(MPI_STATUS_SIZE,size_3d_3v) :: status_3d_3v
  ! for 3 vars, for interpolated grid
  integer, parameter :: size_3d_pp=4*3 ! 2(send/recv)*2(neighbors)*3(variables)
  integer, dimension(size_3d_pp) :: request_3d_pp
  integer, dimension(MPI_STATUS_SIZE,size_3d_pp) :: status_3d_pp
  ! for 1 var, for interpolated grid
  integer, parameter :: size_3d_pp2=4*1 ! 2(send/recv)*2(neighbors)*1(variables)
  integer, dimension(size_3d_pp2) :: request_3d_pp2
  integer, dimension(MPI_STATUS_SIZE,size_3d_pp2) :: status_3d_pp2
  ! type face (ngh_pp layers of cells) for streaks direction
  integer :: type_faceW_interp,type_faceE_interp
  integer :: type_faceF_interp,type_faceB_interp
  integer :: type_mW_interp,type_mE_interp ! for metrics
  ! indices to start communications in the direction [p]arallel to the face
  integer :: ipW_pp,ipE_pp
  ! indices to start communications in the direction [n]ormal to the face
  integer :: inW_pp,inE_pp
  ! ---------------------------------------------------------------------------

contains

  !============================================================================
  subroutine mpi_types_comm_fstt
  !============================================================================
    !> author: AB
    !> date: May 2022
    !> Subroutine for FSTT PP communication initialisation
    !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: nex,ney ! extended size
    integer :: sign_i,sign_j
    integer :: memjump ! jump in memory between block starts
    ! strides in memory between block starts
    integer(kind=MPI_ADDRESS_KIND) :: stride,stride2
    integer :: type_base,type_edge
    ! -------------------------------------------------------------------------

    ! Extended sizes (+ ghost cells)
    ! ------------------------------
    nex=nil_interp+2*ngh_pp
    ney=nj_interp

    ! stride between edges along x or y
    ! ---------------------------------
    stride2=nex*ney*sizeofreal

    ! MPI-type construction for imin face "along y" (W: west)
    ! -------------------------------------------------------
    ! for imin/W, parallel dir is j and normal dir is i
    inW_pp=1
    ipW_pp=1
    sign_i=1
    sign_j=1
    ! if rev parallel -> +j => -j
    if (is_rev(1,1)) then
       sign_j=-1
       ipW_pp=nj_interp
    endif
    ! if rev  normal  -> +i => -i
    if (is_rev(1,2)) then
       sign_i=-1
       inW_pp=ngh_pp
    endif
    if (is_swapij(1)) then
       ! if swap invert role of i and j
       memjump=sign_j*nex  ! j
       stride =sign_i*sizeofreal ! i
       call MPI_TYPE_VECTOR(nj_interp,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ngh_pp,1,stride,type_base,type_mW_interp,info)
    else
       ! default
       memjump=sign_i ! i
       stride =sign_j*nex*sizeofreal ! j
       call MPI_TYPE_VECTOR(ngh_pp,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(nj_interp,1,stride,type_base,type_mW_interp,info)
    endif
    call MPI_TYPE_COMMIT(type_mW_interp,info)
    call MPI_TYPE_CREATE_HVECTOR(nk_interp,1,stride2,type_mW_interp,type_faceW_interp,info)
    call MPI_TYPE_COMMIT(type_faceW_interp,info)

    ! MPI-type construction for imax face "along y" (E: east)
    ! -------------------------------------------------------
    ! for imax/E, parallel dir is j and normal dir is i
    inE_pp=nil_interp+1
    ipE_pp=1
    sign_i=1
    sign_j=1
    ! if rev parallel -> +j => -j
    if (is_rev(2,1)) then
       sign_j=-1
       ipE_pp=nj_interp
    endif
    ! if rev  normal  -> +i => -i
    if (is_rev(2,2)) then
       sign_i=-1
       inE_pp=nil_interp+ngh_pp
    endif
    if (is_swapij(2)) then
       ! if swap invert role of i and j
       memjump=sign_j*nex  ! j
       stride =sign_i*sizeofreal ! i
       call MPI_TYPE_VECTOR(nj_interp,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ngh_pp,1,stride,type_base,type_mE_interp,info)
    else
       ! default
       memjump=sign_i ! i
       stride =sign_j*nex*sizeofreal ! j
       call MPI_TYPE_VECTOR(ngh_pp,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(nj_interp,1,stride,type_base,type_mE_interp,info)
    endif
    call MPI_TYPE_COMMIT(type_mE_interp,info)
    call MPI_TYPE_CREATE_HVECTOR(nk_interp,1,stride2,type_mE_interp,type_faceE_interp,info)
    call MPI_TYPE_COMMIT(type_faceE_interp,info)


    stride=nex*sizeofreal
    call MPI_TYPE_VECTOR(nil_interp+2*ngh_pp,1,1,MPI_DOUBLE_PRECISION,type_base,info)

    ! MPI-type construction for kmin face "along z" (F: forward)
    ! --------------------------------------------------------
    call MPI_TYPE_CREATE_HVECTOR(nj_interp,1, stride,type_base,type_edge,info)
    call MPI_TYPE_CREATE_HVECTOR(nkernel_k,1,stride2,type_edge,type_faceF_interp,info)
    call MPI_TYPE_COMMIT(type_faceF_interp,info)

    ! MPI-type construction for kmax face "along z" (B: backward)
    ! -----------------------------------------------------------
    call MPI_TYPE_CREATE_HVECTOR(nj_interp,1, stride,type_base,type_edge,info)
    call MPI_TYPE_CREATE_HVECTOR(nkernel_k,1,stride2,type_edge,type_faceB_interp,info)
    call MPI_TYPE_COMMIT(type_faceB_interp,info)

  end subroutine mpi_types_comm_fstt

  !===============================================================================
  subroutine commun3d_interp_WE(var,i)
  !===============================================================================
    !> EW communications using non-blocking ISEND/IRECV (faces)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,k
    real(wp), dimension(1-ngh_pp:nil_interp+ngh_pp,1:nj_interp,1-nkernel_k:nk_interp+nkernel_k), intent(inout) :: var
    ! ----------------------------------------------------------------------------

    ! Counter: 2(send/recv)*2(neighbors)=12*(variable number)
    k=4*(i-1)

    ! Send to neighbor W and reception from neighbor W
    call MPI_ISEND(var(inW_pp       ,ipW_pp,1),1,type_faceW_interp,neighbor(nW),tags(nW),COMM_global,request_3d_pp(k+1),info)
    call MPI_IRECV(var(inW_pp-ngh_pp,ipW_pp,1),1,type_faceW_interp,neighbor(nW),tagr(nW),COMM_global,request_3d_pp(k+2),info)

    ! Send to neighbor E and reception from neighbor E
    call MPI_ISEND(var(inE_pp-ngh_pp,ipE_pp,1),1,type_faceE_interp,neighbor(nE),tags(nE),COMM_global,request_3d_pp(k+3),info)
    call MPI_IRECV(var(inE_pp       ,ipE_pp,1),1,type_faceE_interp,neighbor(nE),tagr(nE),COMM_global,request_3d_pp(k+4),info)

  end subroutine commun3d_interp_WE


!===============================================================================
  subroutine commun3d_interp_FB(var,i)
  !===============================================================================
    !> EW communications using non-blocking ISEND/IRECV (faces + edges)
    !> edges ~> with communication neighbor F and B
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,k
    real(wp), dimension(1-ngh_pp:nil_interp+ngh_pp,1:nj_interp,1-nkernel_k:nk_interp+nkernel_k), intent(inout) :: var
    ! ----------------------------------------------------------------------------

    ! Counter: 2(send/recv)*2(neighbors)=12*(variable number)
    k=4*(i-1)

    ! Send to neighbor F and reception from neighbor F
    call MPI_ISEND(var(1-ngh_pp,1,          1),1,type_faceF_interp,neighbor(nF),tags(nF),COMM_global,request_3d_pp(k+1),info)
    call MPI_IRECV(var(1-ngh_pp,1,1-nkernel_k),1,type_faceF_interp,neighbor(nF),tagr(nF),COMM_global,request_3d_pp(k+2),info)

    ! Send to neighbor B and reception from neighbor B
    call MPI_ISEND(var(1-ngh_pp,1,nk_interp-nkernel_k+1),1,type_faceB_interp,neighbor(nB),tags(nB),COMM_global,request_3d_pp(k+3),info)
    call MPI_IRECV(var(1-ngh_pp,1,          nk_interp+1),1,type_faceB_interp,neighbor(nB),tagr(nB),COMM_global,request_3d_pp(k+4),info)

  end subroutine commun3d_interp_FB

  !===============================================================================
  subroutine communication_fstt3(var1,var2,var3)
  !===============================================================================
    !> Communication of values for interpolated field for 3 var
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), dimension(1-ngh_pp:nil_interp+ngh_pp,1:nj_interp,1-nkernel_k:nk_interp+nkernel_k), intent(inout) :: var1,var2,var3
    ! ----------------------------------------------------------------------------

    call commun3d_interp_WE(var1,1)
    call commun3d_interp_WE(var2,2)
    call commun3d_interp_WE(var3,3)

    call MPI_WAITALL(size_3d_pp,request_3d_pp,status_3d_pp,info)

    call commun3d_interp_FB(var1,1)
    call commun3d_interp_FB(var2,2)
    call commun3d_interp_FB(var3,3)

    call MPI_WAITALL(size_3d_pp,request_3d_pp,status_3d_pp,info)



  end subroutine communication_fstt3

  !===============================================================================
  subroutine communication_fstt1(var)
  !===============================================================================
    !> Communication of values for interpolated field for 1 var (faces + edges)
    !> edges ~> with communication neighbor F and B
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: k
    real(wp), dimension(1-ngh_pp:nil_interp+ngh_pp,1:nj_interp,1-nkernel_k:nk_interp+nkernel_k), intent(inout) :: var
    ! ----------------------------------------------------------------------------

    ! Counter: 2(send/recv)*6(neighbors)=12*(variable number)
    k=0

    ! Send to neighbor W and reception from neighbor W
    call MPI_ISEND(var(inW_pp    ,ipW_pp,1),1,type_faceW_interp,neighbor(nW),tags(nW),COMM_global,request_3d_pp2(k+1),info)
    call MPI_IRECV(var(inW_pp-ngh_pp,ipW_pp,1),1,type_faceW_interp,neighbor(nW),tagr(nW),COMM_global,request_3d_pp2(k+2),info)

    ! Send to neighbor E and reception from neighbor E
    call MPI_ISEND(var(inE_pp-ngh_pp,ipE_pp,1),1,type_faceE_interp,neighbor(nE),tags(nE),COMM_global,request_3d_pp2(k+3),info)
    call MPI_IRECV(var(inE_pp       ,ipE_pp,1),1,type_faceE_interp,neighbor(nE),tagr(nE),COMM_global,request_3d_pp2(k+4),info)

    call MPI_WAITALL(size_3d_pp2,request_3d_pp2,status_3d_pp2,info)

    ! Send to neighbor F and reception from neighbor F
    call MPI_ISEND(var(1-ngh_pp,1,          1),1,type_faceF_interp,neighbor(nF),tags(nF),COMM_global,request_3d_pp2(k+1),info)
    call MPI_IRECV(var(1-ngh_pp,1,1-nkernel_k),1,type_faceF_interp,neighbor(nF),tagr(nF),COMM_global,request_3d_pp2(k+2),info)

    ! Send to neighbor B and reception from neighbor B
    call MPI_ISEND(var(1-ngh_pp,1,nk_interp-nkernel_k+1),1,type_faceB_interp,neighbor(nB),tags(nB),COMM_global,request_3d_pp2(k+3),info)
    call MPI_IRECV(var(1-ngh_pp,1,          nk_interp+1),1,type_faceB_interp,neighbor(nB),tagr(nB),COMM_global,request_3d_pp2(k+4),info)

    call MPI_WAITALL(size_3d_pp2,request_3d_pp2,status_3d_pp2,info)


  end subroutine communication_fstt1

  ! !===============================================================================
  ! subroutine communication_3d_3var(var1,var2,var3)
  ! !===============================================================================
  !   !> Call non-blocking communications in 3D (faces only)
  ! !===============================================================================
  !   use mod_comm
  !   implicit none
  !   ! ----------------------------------------------------------------------------
  !   real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(inout) :: var1,var2,var3
  !   ! ----------------------------------------------------------------------------

  !   call commun3d(var1,1)
  !   call commun3d(var2,2)
  !   call commun3d(var3,3)

  !   call MPI_WAITALL(size_3d_3v,request_3d_3v,status_3d_3v,info)

  ! end subroutine communication_3d_3var

end module mod_pp_fstt_comm
