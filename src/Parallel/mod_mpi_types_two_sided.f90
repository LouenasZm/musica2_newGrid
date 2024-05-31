!==================================================================================
module mod_mpi_types_two_sided
!==================================================================================
  !> Module to create MPI types
!==================================================================================
  use mod_mpi_part
  use mod_grid_directions
  implicit none
  ! -------------------------------------------------------------------------------
  ! ngh layers of cells [flux_euler & cie]
  ! ===================
  ! send/recv indices
  integer, dimension(6) :: iis,ijs,iks,iir,ijr,ikr
  ! MPI types for communications
  integer, dimension(6) :: type_face
  ! MPI types for metrics
  integer, dimension(6) :: type_met
  ! -------------------------------------------------------------------------------
  ! ngh_v layers of cells [flux_visc & cie]
  ! =====================
  ! send/recv indices
  integer, dimension(6) :: iis_v,ijs_v,iks_v,iir_v,ijr_v,ikr_v
  ! MPI types for communications
  integer, dimension(6) :: type_face_v
  ! MPI types for metrics
  integer, dimension(6) :: type_met_v
  ! -------------------------------------------------------------------------------
  ! nghirs layers of cells [increments for IRS]
  ! ======================
  ! send/recv indices
  integer, dimension(6) :: iis_i,ijs_i,iks_i,iir_i,ijr_i,ikr_i
  ! MPI types for communications
  integer, dimension(6) :: type_face_i
  ! -------------------------------------------------------------------------------
  ! ngh layers of cells [flux_euler & cie] -> extended
  ! ===================
  ! send/recv indices
  integer, dimension(6) :: iis_e,ijs_e,iks_e,iir_e,ijr_e,ikr_e
  ! MPI types for communications
  integer, dimension(6) :: type_face_e
  ! MPI types for communications
  integer, dimension(6) :: type_face_e2d

contains

  !================================================================================
  subroutine mpi_types_comm(ngh_,type)
  !================================================================================
    !> Definition of MPI types for communications
    !> type face (ngh layers of cells) for each direction [flux_euler & cie]
  !================================================================================
    implicit none
    ! -----------------------------------------------------------------------------
    integer, intent(in) :: ngh_
    character, intent(in) :: type
    ! -----------------------------------------------------------------------------
    integer :: n
    integer :: n1,n2,n3,n_tmp
    integer :: nex,ney ! extended sizes
    integer :: sign_i,sign_j,sign_k
    ! jump in memory per direction
    integer :: memjump1,memjump2,memjump3,memj_tmp
    ! strides in memory
    integer(kind=MPI_ADDRESS_KIND) :: stride,stride2
    integer :: type_base,type_edge
    ! -----------------------------------------------------------------------------

    ! Extended sizes (+ ghost cells)
    ! ------------------------------
    nex=nx+2*ngh_
    ney=ny+2*ngh_

    do n=1,6
       ! imin (n=1), imax (n=2)
       ! jmin (n=3), jmax (n=4)
       ! kmin (n=5), kmax (n=6)

       if (type=='r') then
          call determine_mpi_face(n,iis(n),ijs(n),iks(n),iir(n),ijr(n),ikr(n), &
                                  sign_i,sign_j,sign_k,ngh_)
       elseif (type=='v') then
          call determine_mpi_face(n,iis_v(n),ijs_v(n),iks_v(n),iir_v(n),ijr_v(n),ikr_v(n), &
                                  sign_i,sign_j,sign_k,ngh_)
       elseif (type=='i') then
          call determine_mpi_face(n,iis_i(n),ijs_i(n),iks_i(n),iir_i(n),ijr_i(n),ikr_i(n), &
                                  sign_i,sign_j,sign_k,ngh_)
       endif

       ! default values
       ! --------------
       ! memory jumps with reverse sign corrections
       memjump1=sign_i ! i
       memjump2=sign_j*nex ! j
       memjump3=sign_k*nex*ney ! k
       ! faces imin/imax
       if ((n==1).or.(n==2)) then
          n1=ngh_
          n2=ny
          n3=nz
       ! faces jmin/jmax
       elseif ((n==3).or.(n==4)) then
          n1=nx
          n2=ngh_
          n3=nz
       ! faces kmin/kmax
       elseif ((n==5).or.(n==6)) then
          n1=nx
          n2=ny
          n3=ngh_
       endif

       ! if swap_ij -> invert role of i and j
       ! ------------------------------------
       if (is_swapij(n)) then
          memj_tmp=memjump1
          memjump1=memjump2
          memjump2=memj_tmp
          n_tmp=n1
          n1=n2
          n2=n_tmp
       endif

       ! if swap_ik -> invert role of i and k
       ! ------------------------------------
       if (is_swapik(n)) then
          memj_tmp=memjump1
          memjump1=memjump3
          memjump3=memj_tmp
          n_tmp=n1
          n1=n3
          n3=n_tmp
       endif

       ! if swap_jk -> invert role of j and k
       ! ------------------------------------
       if (is_swapjk(n)) then
          memj_tmp=memjump2
          memjump2=memjump3
          memjump3=memj_tmp
          n_tmp=n2
          n2=n3
          n3=n_tmp
       endif

       ! if swap_ijk -> permute role of i,j,k
       ! ------------------------------------
       if (is_swapijk(n)) then
          memj_tmp=memjump1
          memjump1=memjump3
          memjump3=memjump2
          memjump2=memj_tmp
          n_tmp=n1
          n1=n3
          n3=n2
          n2=n_tmp
       endif

       ! MPI-type construction
       ! ---------------------

       stride=memjump2*sizeofreal
       stride2=memjump3*sizeofreal

       call MPI_TYPE_VECTOR(n1,1,memjump1,MPI_DOUBLE_PRECISION,type_base,info)

       ! ngh layers of cells [flux_euler & cie]
       ! ===================
       if (type=='r') then
          call MPI_TYPE_CREATE_HVECTOR(n2,1,stride,type_base,type_met(n),info)
          ! commit 2D type for metrics
          call MPI_TYPE_COMMIT(type_met(n),info)
          ! extend in 3D and commit
          call MPI_TYPE_CREATE_HVECTOR(n3,1,stride2,type_met(n),type_face(n),info)
          call MPI_TYPE_COMMIT(type_face(n),info)

       ! ngh_v layers of cells [flux_visc & cie]
       ! =====================
       elseif (type=='v') then
          call MPI_TYPE_CREATE_HVECTOR(n2,1,stride,type_base,type_met_v(n),info)
          ! commit 2D type for metrics
          call MPI_TYPE_COMMIT(type_met_v(n),info)
          ! extend in 3D and commit
          call MPI_TYPE_CREATE_HVECTOR(n3,1,stride2,type_met_v(n),type_face_v(n),info)
          call MPI_TYPE_COMMIT(type_face_v(n),info)

       ! nghirs layers of cells [increments for IRS]
       ! ======================
       elseif (type=='i') then
          call MPI_TYPE_CREATE_HVECTOR(n2,1,stride,type_base,type_edge,info)
          ! extend in 3D and commit
          call MPI_TYPE_CREATE_HVECTOR(n3,1,stride2,type_edge,type_face_i(n),info)
          call MPI_TYPE_COMMIT(type_face_i(n),info)
       endif

    enddo

  end subroutine mpi_types_comm

  !================================================================================
  subroutine mpi_types_comm_ex(ngh_)
  !================================================================================
    !> Definition of MPI types for communications
    !> type face (ngh layers of cells) for each direction [flux_euler & cie]
  !================================================================================
    implicit none
    ! -----------------------------------------------------------------------------
    integer, intent(in) :: ngh_
    ! -----------------------------------------------------------------------------
    integer :: n
    integer :: n1,n2,n3,n_tmp
    integer :: nex,ney,nez ! extended sizes
    integer :: sign_i,sign_j,sign_k
    ! jump in memory per direction
    integer :: memjump1,memjump2,memjump3,memj_tmp
    ! strides in memory
    integer(kind=MPI_ADDRESS_KIND) :: stride,stride2
    integer :: type_base
    ! -----------------------------------------------------------------------------

    ! Extended sizes (+ ghost cells)
    ! ------------------------------
    nex=nx+2*ngh_
    ney=ny+2*ngh_
    nez=nz+2*ngh_

    ! indices
    ! iis(n),ijs(n),iks(n),iir(n),ijr(n),ikr(n)

    do n=1,6
       ! imin (n=1), imax (n=2)
       ! jmin (n=3), jmax (n=4)
       ! kmin (n=5), kmax (n=6)

       ! if (ngh_==ngh) then
          call determine_mpi_face_e(n,iis_e(n),ijs_e(n),iks_e(n),iir_e(n),ijr_e(n),ikr_e(n), &
                                    sign_i,sign_j,sign_k,ngh_)
       ! endif

       ! default values
       ! --------------
       ! memory jumps with reverse sign corrections
       memjump1=sign_i ! i
       memjump2=sign_j*nex ! j
       memjump3=sign_k*nex*ney ! k
       ! faces imin/imax
       if ((n==1).or.(n==2)) then
          n1=ngh_
          n2=ney
          n3=nez
       ! faces jmin/jmax
       elseif ((n==3).or.(n==4)) then
          n1=nex
          n2=ngh_
          n3=nez
       ! faces kmin/kmax
       elseif ((n==5).or.(n==6)) then
          n1=nex
          n2=ney
          n3=ngh_
       endif

       ! if swap_ij -> invert role of i and j
       ! ------------------------------------
       if (is_swapij(n)) then
          memj_tmp=memjump1
          memjump1=memjump2
          memjump2=memj_tmp
          n_tmp=n1
          n1=n2
          n2=n_tmp
       endif

       ! if swap_ik -> invert role of i and k
       ! ------------------------------------
       if (is_swapik(n)) then
          memj_tmp=memjump1
          memjump1=memjump3
          memjump3=memj_tmp
          n_tmp=n1
          n1=n3
          n3=n_tmp
       endif

       ! if swap_jk -> invert role of j and k
       ! ------------------------------------
       if (is_swapjk(n)) then
          memj_tmp=memjump2
          memjump2=memjump3
          memjump3=memj_tmp
          n_tmp=n2
          n2=n3
          n3=n_tmp
       endif

       ! if swap_ijk -> permute role of i,j,k
       ! ------------------------------------
       if (is_swapijk(n)) then
          memj_tmp=memjump1
          memjump1=memjump3
          memjump3=memjump2
          memjump2=memj_tmp
          n_tmp=n1
          n1=n3
          n3=n2
          n2=n_tmp
       endif

       ! MPI-type construction
       ! ---------------------

       stride=memjump2*sizeofreal
       stride2=memjump3*sizeofreal

       call MPI_TYPE_VECTOR(n1,1,memjump1,MPI_DOUBLE_PRECISION,type_base,info)

       ! ngh layers of cells [flux_euler & cie]
       ! ===================
       ! if (ngh_==ngh) then
          call MPI_TYPE_CREATE_HVECTOR(n2,1,stride,type_base,type_face_e2d(n),info)
          call MPI_TYPE_COMMIT(type_face_e2d(n),info)
          ! extend in 3D and commit
          call MPI_TYPE_CREATE_HVECTOR(n3,1,stride2,type_face_e2d(n),type_face_e(n),info)
          call MPI_TYPE_COMMIT(type_face_e(n),info)
       ! endif

    enddo

  end subroutine mpi_types_comm_ex

  !================================================================================
  subroutine determine_mpi_face(n_f,ii_send,ij_send,ik_send, &
                                    ii_recv,ij_recv,ik_recv, &
                                    sign_i,sign_j,sign_k,ngh_)
  !================================================================================
    !> Determine face characteristics for MPI exchange
    !> -> send/recv indices and sign for reverse directions
  !================================================================================
    implicit none
    ! -----------------------------------------------------------------------------
    integer, intent(in) :: n_f,ngh_
    integer, intent(out) :: sign_i,sign_j,sign_k
    integer, intent(out) :: ii_send,ij_send,ik_send
    integer, intent(out) :: ii_recv,ij_recv,ik_recv
    ! -----------------------------------------------------------------------------
    integer :: dir,n_n,n_p1,n_p2
    integer :: first_proc,last_proc
    integer :: in_send,in_recv
    integer :: ip1_send,ip1_recv,ip2_send,ip2_recv
    integer :: sign_n,sign_p1,sign_p2
    ! -----------------------------------------------------------------------------

    ! Determine direction: 1=imin,imax; 2=jmin,jmax; 3=kmin,kmax;
    ! ====================
    if ((n_f==1).or.(n_f==2)) dir=1
    if ((n_f==3).or.(n_f==4)) dir=2
    if ((n_f==5).or.(n_f==6)) dir=3

    ! Initialize face characteristics
    ! ===============================
    first_proc=0

    if (dir==1) then
       ! i is the normal direction
       n_n=nx
       last_proc=ndomx-1
       ! j,k are parallel directions
       n_p1=ny
       n_p2=nz
    elseif (dir==2) then
       ! j is the normal direction
       n_n=ny
       last_proc=ndomy-1
       ! i,k are parallel directions
       n_p1=nx
       n_p2=nz
    elseif (dir==3) then
       ! k is the normal direction
       n_n=nz
       last_proc=ndomz-1
       ! i,j are parallel directions
       n_p1=nx
       n_p2=ny
    endif

    ! default values for send/recv indices in normal direction
    if (mod(n_f,2)==0) then ! for imax jmax kmax
       in_send=n_n+1-ngh_
       in_recv=n_n+1
    else ! for imin jmin kmin
       in_send=1
       in_recv=1-ngh_
    endif

    ! default values for send/recv indices in parallel directions
    ip1_send=1
    ip1_recv=1
    ip2_send=1
    ip2_recv=1

    ! default values for signs
    sign_n =1
    sign_p1=1
    sign_p2=1

    ! Correction at block interface for adjoint version
    ! -------------------------------------------------
    if (is_adjoint_block) then
       if (mod(n_f,2)==0) then ! for imax jmax kmax
          if (neighbor(n_f).ne.MPI_PROC_NULL) then
             ! block interface
             if (nob(neighbor(n_f)).ne.nob(iproc)) in_send=in_send-1
             ! particular case of periodicity (O-grid for instance)
             if ((nob(neighbor(n_f))==nob(iproc)).and.(coord(dir)==last_proc)) in_send=in_send-1
          endif
       else ! for imin jmin kmin
          if (neighbor(n_f).ne.MPI_PROC_NULL) then
             ! block interface
             if (nob(neighbor(n_f)).ne.nob(iproc)) in_send=in_send+1
             ! particular case of periodicity (O-grid for instance)
             if ((nob(neighbor(n_f))==nob(iproc)).and.(coord(dir)==first_proc)) in_send=in_send+1
          endif
       endif
    endif

    ! if rev parallel 1
    ! -----------------
    if (is_rev(n_f,1)) then
       ! -> parallel direction #1 starts from end
       ip1_send=n_p1
       ip1_recv=n_p1
       ! -> parallel direction #1 is described in reverse sense
       sign_p1=-1
    endif

    ! if rev parallel 2
    ! -----------------
    if (is_rev(n_f,3)) then
       ! -> parallel direction #2 starts from end
       ip2_send=n_p2
       ip2_recv=n_p2
       ! -> parallel direction #2 is described in reverse sense
       sign_p2=-1
    endif

    ! if rev  normal
    ! --------------
    if (is_rev(n_f,2)) then
       if (mod(n_f,2)==0) then ! for imax jmax kmax
          if (is_adjoint_block) then
             in_send=n_n-1
          else
             !! nx+1-ngh_ -> nx+1-ngh_ +ngh_-1
             in_send=n_n
          endif
          !! nx+1 -> nx+1 +ngh_-1
          in_recv=n_n+ngh_
       else ! for imin jmin kmin
          if (is_adjoint_block) then
             in_send=ngh_+1
          else
             !! 1 -> 1+ngh_-1
             in_send=ngh_
          endif
          !! 1-ngh_ -> 1-ngh_+ngh_-1
          in_recv=0
       endif

       ! -> normal direction is described in reverse sense
       sign_n=-1
    endif

    ! Attribute face characteristics
    ! ==============================
    if (dir==1) then
       ! i is the normal direction
       ii_send=in_send
       ii_recv=in_recv
       sign_i=sign_n
       ! j,k are parallel directions
       ij_send=ip1_send
       ij_recv=ip1_recv
       sign_j=sign_p1
       ik_send=ip2_send
       ik_recv=ip2_recv
       sign_k=sign_p2
    elseif (dir==2) then
       ! j is the normal direction
       ij_send=in_send
       ij_recv=in_recv
       sign_j=sign_n
       ! i,k are parallel directions
       ii_send=ip1_send
       ii_recv=ip1_recv
       sign_i=sign_p1
       ik_send=ip2_send
       ik_recv=ip2_recv
       sign_k=sign_p2
    elseif (dir==3) then
       ! k is the normal direction
       ik_send=in_send
       ik_recv=in_recv
       sign_k=sign_n
       ! i,j are parallel directions
       ii_send=ip1_send
       ii_recv=ip1_recv
       sign_i=sign_p1
       ij_send=ip2_send
       ij_recv=ip2_recv
       sign_j=sign_p2
    endif

  end subroutine determine_mpi_face

  !================================================================================
  subroutine determine_mpi_face_e(n_f,ii_send,ij_send,ik_send, &
                                      ii_recv,ij_recv,ik_recv, &
                                      sign_i,sign_j,sign_k,ngh_)
  !================================================================================
    !> Determine face characteristics for MPI exchange
    !> -> send/recv indices and sign for reverse directions
  !================================================================================
    implicit none
    ! -----------------------------------------------------------------------------
    integer, intent(in) :: n_f,ngh_
    integer, intent(out) :: sign_i,sign_j,sign_k
    integer, intent(out) :: ii_send,ij_send,ik_send
    integer, intent(out) :: ii_recv,ij_recv,ik_recv
    ! -----------------------------------------------------------------------------
    integer :: dir,n_n,n_p1,n_p2
    integer :: first_proc,last_proc
    integer :: in_send,in_recv
    integer :: ip1_send,ip1_recv,ip2_send,ip2_recv
    integer :: sign_n,sign_p1,sign_p2
    ! -----------------------------------------------------------------------------

    ! Determine direction: 1=imin,imax; 2=jmin,jmax; 3=kmin,kmax;
    ! ====================
    if ((n_f==1).or.(n_f==2)) dir=1
    if ((n_f==3).or.(n_f==4)) dir=2
    if ((n_f==5).or.(n_f==6)) dir=3

    ! Initialize face characteristics
    ! ===============================
    first_proc=0

    if (dir==1) then
       ! i is the normal direction
       n_n=nx
       last_proc=ndomx-1
       ! j,k are parallel directions
       n_p1=ny+ngh_
       n_p2=nz+ngh_
    elseif (dir==2) then
       ! j is the normal direction
       n_n=ny
       last_proc=ndomy-1
       ! i,k are parallel directions
       n_p1=nx+ngh_
       n_p2=nz+ngh_
    elseif (dir==3) then
       ! k is the normal direction
       n_n=nz
       last_proc=ndomz-1
       ! i,j are parallel directions
       n_p1=nx+ngh_
       n_p2=ny+ngh_
    endif

    ! default values for send/recv indices in normal direction
    if (mod(n_f,2)==0) then ! for imax jmax kmax
       in_send=n_n+1-ngh_
       in_recv=n_n+1
    else ! for imin jmin kmin
       in_send=1
       in_recv=1-ngh_
    endif

    ! default values for send/recv indices in parallel directions
    ip1_send=1-ngh_
    ip1_recv=1-ngh_
    ip2_send=1-ngh_
    ip2_recv=1-ngh_

    ! default values for signs
    sign_n =1
    sign_p1=1
    sign_p2=1

    ! Correction at block interface for adjoint version
    ! -------------------------------------------------
    if (is_adjoint_block) then
       if (mod(n_f,2)==0) then ! for imax jmax kmax
          if (neighbor(n_f).ne.MPI_PROC_NULL) then
             ! block interface
             if (nob(neighbor(n_f)).ne.nob(iproc)) in_send=in_send-1
             ! particular case of periodicity (O-grid for instance)
             if ((nob(neighbor(n_f))==nob(iproc)).and.(coord(dir)==last_proc)) in_send=in_send-1
          endif
       else ! for imin jmin kmin
          if (neighbor(n_f).ne.MPI_PROC_NULL) then
             ! block interface
             if (nob(neighbor(n_f)).ne.nob(iproc)) in_send=in_send+1
             ! particular case of periodicity (O-grid for instance)
             if ((nob(neighbor(n_f))==nob(iproc)).and.(coord(dir)==first_proc)) in_send=in_send+1
          endif
       endif
    endif

    ! if rev parallel 1
    ! -----------------
    if (is_rev(n_f,1)) then
       ! -> parallel direction #1 starts from end
       ip1_send=n_p1
       ip1_recv=n_p1
       ! -> parallel direction #1 is described in reverse sense
       sign_p1=-1
    endif

    ! if rev parallel 2
    ! -----------------
    if (is_rev(n_f,3)) then
       ! -> parallel direction #2 starts from end
       ip2_send=n_p2
       ip2_recv=n_p2
       ! -> parallel direction #2 is described in reverse sense
       sign_p2=-1
    endif

    ! if rev  normal
    ! --------------
    if (is_rev(n_f,2)) then
       if (mod(n_f,2)==0) then ! for imax jmax kmax
          if (is_adjoint_block) then
             in_send=n_n-1
          else
             !! nx+1-ngh_ -> nx+1-ngh_ +ngh_-1
             in_send=n_n
          endif
          !! nx+1 -> nx+1 +ngh_-1
          in_recv=n_n+ngh_
       else ! for imin jmin kmin
          if (is_adjoint_block) then
             in_send=ngh_+1
          else
             !! 1 -> 1+ngh_-1
             in_send=ngh_
          endif
          !! 1-ngh_ -> 1-ngh_+ngh_-1
          in_recv=0
       endif

       ! -> normal direction is described in reverse sense
       sign_n=-1
    endif

    ! Attribute face characteristics
    ! ==============================
    if (dir==1) then
       ! i is the normal direction
       ii_send=in_send
       ii_recv=in_recv
       sign_i=sign_n
       ! j,k are parallel directions
       ij_send=ip1_send
       ij_recv=ip1_recv
       sign_j=sign_p1
       ik_send=ip2_send
       ik_recv=ip2_recv
       sign_k=sign_p2
    elseif (dir==2) then
       ! j is the normal direction
       ij_send=in_send
       ij_recv=in_recv
       sign_j=sign_n
       ! i,k are parallel directions
       ii_send=ip1_send
       ii_recv=ip1_recv
       sign_i=sign_p1
       ik_send=ip2_send
       ik_recv=ip2_recv
       sign_k=sign_p2
    elseif (dir==3) then
       ! k is the normal direction
       ik_send=in_send
       ik_recv=in_recv
       sign_k=sign_n
       ! i,j are parallel directions
       ii_send=ip1_send
       ii_recv=ip1_recv
       sign_i=sign_p1
       ij_send=ip2_send
       ij_recv=ip2_recv
       sign_j=sign_p2
    endif

  end subroutine determine_mpi_face_e

  !===============================================================
  subroutine free_mpi_types_comm
  !===============================================================
    !> Free MPI types for communications of block faces
  !===============================================================
    implicit none
    ! ------------------------------------------------------------
    integer :: n
    ! ------------------------------------------------------------

    ! MPI types for communications of conservative variables
    ! [inviscid fluxes]
    do n=1,6
       call MPI_TYPE_FREE(type_face(n),info)
    enddo

    ! MPI types for communications of velocity derivatives
    ! [viscous fluxes]
    do n=1,6
       call MPI_TYPE_FREE(type_face_v(n),info)
    enddo

  end subroutine free_mpi_types_comm
  
  !===============================================================
  subroutine free_mpi_types_comm_ex
  !===============================================================
    !> Free MPI types for communications of extended faces
  !===============================================================
    implicit none
    ! ------------------------------------------------------------
    integer :: n
    ! ------------------------------------------------------------
    
    ! MPI types for communications of conservative variables
    ! [inviscid fluxes]
    do n=1,6
       call MPI_TYPE_FREE(type_face_e(n),info)
    enddo

  end subroutine free_mpi_types_comm_ex
  
  !===============================================================
  subroutine free_mpi_types_comm_inc
  !===============================================================
    !> Free MPI types for communications of increments (IRS)
  !===============================================================
    implicit none
    ! ------------------------------------------------------------
    integer :: n
    ! ------------------------------------------------------------

    ! MPI types for communications of conservative variables
    do n=1,6
       call MPI_TYPE_FREE(type_face_i(n),info)
    enddo

  end subroutine free_mpi_types_comm_inc
  
  !===============================================================
  subroutine free_mpi_types_metrics
  !===============================================================
    !> Free MPI types for communications of metrics
  !===============================================================
    implicit none
    ! ------------------------------------------------------------
    integer :: n
    ! ------------------------------------------------------------

    ! MPI types for communications of 2D curvilinear metrics
    ! [inviscid fluxes]
    do n=1,6
       call MPI_TYPE_FREE(type_met(n),info)
    enddo

    ! MPI types for communications of 2D curvilinear metrics
    ! [viscous fluxes]
    do n=1,6
       call MPI_TYPE_FREE(type_met_v(n),info)
    enddo

  end subroutine free_mpi_types_metrics

end module mod_mpi_types_two_sided
