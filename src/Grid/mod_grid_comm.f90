!=================================================================================
module mod_grid_comm
!=================================================================================
  !> Module to define subroutines for grid communications
!=================================================================================
  use mod_ngh
  use mod_block
  use mod_mpi_part
  implicit none
  ! ----------------------------------------------------------------------------
  integer :: nbl
  ! ----------------------------------------------------------------------------
  ! MPI types for grid communications
  integer :: type_x
  integer :: type_x_facex,type_x_facey,type_x_facez
  integer :: type_x_faceW,type_x_faceE
  integer :: type_x_faceS,type_x_faceN
  integer :: type_x_faceF,type_x_faceB
  ! indices to start communications in the direction [p]arallel to the block face
  integer :: ipW_bl,ipE_bl,ipS_bl,ipN_bl
  ! indices to start communications in the direction [n]ormal to the block face
  integer :: inW_bl,inE_bl,inS_bl,inN_bl
  ! for the cases of adjoint block interfaces
  ! indices to start send routines in the direction [n]ormal to the block face
  integer :: inWs_bl,inEs_bl,inSs_bl,inNs_bl
  ! indices to start receive routines in the direction [n]ormal to the block face
  integer :: inWr_bl,inEr_bl,inSr_bl,inNr_bl
  ! ----------------------------------------------------------------------------
  ! send/recv indices
  integer, dimension(6) :: iis_bl,ijs_bl,iks_bl,iir_bl,ijr_bl,ikr_bl
  ! ----------------------------------------------------------------------------
  ! type face grid (ngh layers of cells) for each direction [flux_euler & cie]
  integer, dimension(6) :: type_x_face ! for communications
  ! ----------------------------------------------------------------------------
contains
  
  !===============================================================================
  subroutine grid_comm_type_cart
  !===============================================================================
    !> Definition of MPI types for Cartesian grid communications
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    ! MPI type for 1D exchange of size ngh
    ! ------------------------------------
    call MPI_TYPE_VECTOR(ngh,1,1,MPI_DOUBLE_PRECISION,type_x,info)
    call MPI_TYPE_COMMIT(type_x,info)

  end subroutine grid_comm_type_cart

  !===============================================================================
  subroutine grid_comm_type_curv_2(nx_,ny_,nz_)
  !===============================================================================
    !> Definition of MPI types for interblock grid communications
    !> * 2D/3D Curvilinear version *
  !===============================================================================
    use mod_grid_directions
    use warnstop
    implicit none
    ! ----------------------------------------------------------------------------
    integer, intent(in) :: nx_,ny_,nz_  ! sizes
    ! ----------------------------------------------------------------------------
    integer :: n
    integer :: n1,n2,n3,n_tmp
    integer :: nex,ney,nez ! extended sizes
    integer :: sign_i,sign_j,sign_k
    ! jump in memory per direction
    integer :: memjump1,memjump2,memjump3,memj_tmp
    ! strides in memory
    integer(kind=MPI_ADDRESS_KIND) :: stride,stride2
    integer :: type_base,type_edge
    ! ----------------------------------------------------------------------------
    integer i,k
    integer :: type_test

    ! Extended sizes (+ ghost cells)
    ! ------------------------------
    nex=nx_+2*ngh
    ney=ny_+2*ngh
    nez=nz_+2*ngh

    ! MPI types for intrablock grid comm
    ! ----------------------------------
    stride=nex*sizeofreal
    stride2=nex*ney*sizeofreal

    if (is_curv3) then
       ! Construction of type "face along y" (west or east)
       call MPI_TYPE_VECTOR(ngh,1,1,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ney,1,stride,type_base,type_edge,info)
       call MPI_TYPE_CREATE_HVECTOR(nez,1,stride2,type_edge,type_x_facey,info)
       call MPI_TYPE_COMMIT(type_x_facey,info)

       ! Construction of type "face along x" (north or south)
       call MPI_TYPE_VECTOR(nex,1,1,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride,type_base,type_edge,info)
       call MPI_TYPE_CREATE_HVECTOR(nez,1,stride2,type_edge,type_x_facex,info)
       call MPI_TYPE_COMMIT(type_x_facex,info)

       ! Construction of type "face along z" (forward or backward)
       call MPI_TYPE_CREATE_HVECTOR(ney,1, stride,type_base,type_edge,info)
       call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride2,type_edge,type_x_facez,info)
       call MPI_TYPE_COMMIT(type_x_facez,info)
    else
       ! Construction of type "face along y" (west or east)
       call MPI_TYPE_VECTOR(ngh,1,1,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ny_+2*ngh,1,stride,type_base,type_x_facey,info)
       call MPI_TYPE_COMMIT(type_x_facey,info)

       ! Construction of type "face along x" (north or south)
       call MPI_TYPE_VECTOR(nx_+2*ngh,1,1,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride,type_base,type_x_facex,info)
       call MPI_TYPE_COMMIT(type_x_facex,info)
    endif

    ! MPI types for interblock grid comm
    ! ----------------------------------
    do n=1,6

       call determine_mpi_facex(n,iis_bl(n),ijs_bl(n),iks_bl(n), &
                                  iir_bl(n),ijr_bl(n),ikr_bl(n), &
                                  sign_i,sign_j,sign_k)

       !print *,sign_i,sign_j,sign_k,'proc',iproc

       ! default values
       ! --------------
       ! memory jumps with reverse sign corrections
       memjump1=sign_i ! i
       memjump2=sign_j*nex ! j
       memjump3=sign_k*nex*ney ! k
       ! faces imin/imax
       if ((n==1).or.(n==2)) then
          n1=ngh
          n2=ney
          n3=nez
       ! faces jmin/jmax
       elseif ((n==3).or.(n==4)) then
          n1=nex
          n2=ngh
          n3=nez
       ! faces kmin/kmax
       elseif ((n==5).or.(n==6)) then
          n1=nex
          n2=ney
          n3=ngh
       endif

       ! if swap_ij -> invert role of i and j
       ! ------------------------------------
       if (is_swapij_bl(n)) then
          !print *,'swap ij proc',iproc, 'face',n
          !print *,n1,n2,nx_,ny_,nz_
          !print *,memjump1,memjump2,memjump3
          !print *,'----------'
          memj_tmp=memjump1
          memjump1=memjump2
          memjump2=memj_tmp
          n_tmp=n1
          n1=n2
          n2=n_tmp
          !print *,n1,n2
          !print *,memjump1,memjump2,memjump3
       endif

       ! if swap_ik -> invert role of i and k
       ! ------------------------------------
       if (is_swapik_bl(n)) then
          !print *,'swap ik proc',iproc
          !print *,n1,n3,nx_,ny_,nz_
          !print *,memjump1,memjump2,memjump3
          !print *,'----------'
          memj_tmp=memjump1
          memjump1=memjump3
          memjump3=memj_tmp
          n_tmp=n1
          n1=n3
          n3=n_tmp
          !print *,n1,n3
          !print *,memjump1,memjump2,memjump3
       endif

       ! if swap_jk -> invert role of j and k
       ! ------------------------------------
       if (is_swapjk_bl(n)) then
          memj_tmp=memjump2
          memjump2=memjump3
          memjump3=memj_tmp
          n_tmp=n2
          n2=n3
          n3=n_tmp
       endif

       ! if swap_ijk -> permute role of i,j,k
       ! ------------------------------------
       if (is_swapijk_bl(n)) then
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

       if (is_curv3) then
          call MPI_TYPE_CREATE_HVECTOR(n2,1,stride,type_base,type_edge,info)
          ! extend in 3D and commit 3D type for grid
          call MPI_TYPE_CREATE_HVECTOR(n3,1,stride2,type_edge,type_x_face(n),info)
          call MPI_TYPE_COMMIT(type_x_face(n),info)
       else
          call MPI_TYPE_CREATE_HVECTOR(n2,1,stride,type_base,type_x_face(n),info)
          ! commit 2D type for grid
          call MPI_TYPE_COMMIT(type_x_face(n),info)
       endif

    enddo

!!$    if (iproc==0) then
!!$       print *,'proc',iproc
!!$       do i=1,5
!!$          print *,i,xgc3(i,1,1:5)
!!$       enddo
!!$       ! create type test
!!$       !call MPI_TYPE_VECTOR(5,1,1,MPI_DOUBLE_PRECISION,type_base,info)
!!$       !call MPI_TYPE_CREATE_HVECTOR(1,1,nex*sizeofreal,type_base,type_edge,info)
!!$       !call MPI_TYPE_CREATE_HVECTOR(5,1,nex*ney*sizeofreal,type_edge,type_test,info)
!!$       !call MPI_TYPE_COMMIT(type_test,info)
!!$    endif
!!$    if (iproc==1) then
!!$       print *,'proc',iproc
!!$       do i=1,5
!!$          print *,i,xgc3(i,ny_,1:5)
!!$       enddo
!!$       ! create type test
!!$       !call MPI_TYPE_VECTOR(5,1,1,MPI_DOUBLE_PRECISION,type_base,info)
!!$       !call MPI_TYPE_CREATE_HVECTOR(1,1,nex*sizeofreal,type_base,type_edge,info)
!!$       !call MPI_TYPE_CREATE_HVECTOR(5,1,nex*ney*sizeofreal,type_edge,type_test,info)
!!$       !call MPI_TYPE_COMMIT(type_test,info)
!!$    endif
!!$    call mpistop('here',0)

  contains

    !===========================================================================
    subroutine determine_mpi_facex(n_f,ii_send,ij_send,ik_send, &
                                       ii_recv,ij_recv,ik_recv, &
                                       sign_i,sign_j,sign_k)
    !===========================================================================
      !> Determine face characteristics for MPI exchange
      !> -> send/recv indices and sign for reverse directions
    !===========================================================================
      implicit none
      ! ------------------------------------------------------------------------
      integer, intent(in) :: n_f
      integer, intent(out) :: sign_i,sign_j,sign_k
      integer, intent(out) :: ii_send,ij_send,ik_send
      integer, intent(out) :: ii_recv,ij_recv,ik_recv
      ! ------------------------------------------------------------------------
      integer :: dir,n_n,n_p1,n_p2
      integer :: first_proc,last_proc
      integer :: in_send,in_recv
      integer :: ip1_send,ip1_recv,ip2_send,ip2_recv
      integer :: sign_n,sign_p1,sign_p2
      ! ------------------------------------------------------------------------

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
         n_n=nx_
         last_proc=ndomx-1
         ! j,k are parallel directions
         n_p1=ny_+ngh
         n_p2=nz_+ngh
      elseif (dir==2) then
         ! j is the normal direction
         n_n=ny_
         last_proc=ndomy-1
         ! i,k are parallel directions
         n_p1=nx_+ngh
         n_p2=nz_+ngh
      elseif (dir==3) then
         ! k is the normal direction
         n_n=nz_
         last_proc=ndomz-1
         ! i,j are parallel directions
         n_p1=nx_+ngh
         n_p2=ny_+ngh
      endif

      ! default values for send/recv indices in normal direction
      if (mod(n_f,2)==0) then ! for imax jmax kmax
         in_send=n_n+1-ngh
         in_recv=n_n+1
      else ! for imin jmin kmin
         in_send=1
         in_recv=1-ngh
      endif

      ! default values for send/recv indices in parallel directions
      ip1_send=1-ngh
      ip1_recv=1-ngh
      ip2_send=1-ngh
      ip2_recv=1-ngh

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
      if (is_rev_bl(n_f,1)) then
         ! -> parallel direction #1 starts from end
         ip1_send=n_p1
         ip1_recv=n_p1
         ! -> parallel direction #1 is described in reverse sense
         sign_p1=-1
      endif

      ! if rev parallel 2
      ! -----------------
      if (is_rev_bl(n_f,3)) then
         ! -> parallel direction #2 starts from end
         ip2_send=n_p2
         ip2_recv=n_p2
         ! -> parallel direction #2 is described in reverse sense
         sign_p2=-1
      endif

      ! if rev  normal
      ! --------------
      if (is_rev_bl(n_f,2)) then
         if (mod(n_f,2)==0) then ! for imax jmax kmax
            if (is_adjoint_block) then
               in_send=n_n-1
            else
               !! nx_+1-ngh -> nx_+1-ngh +ngh-1
               in_send=n_n
            endif
            !! nx_+1 -> nx_+1 +ngh-1
            in_recv=n_n+ngh
         else ! for imin jmin kmin
            if (is_adjoint_block) then
               in_send=ngh+1
            else
               !! 1 -> 1+ngh-1
               in_send=ngh
            endif
            !! 1-ngh -> 1-ngh+ngh-1
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

    end subroutine determine_mpi_facex

  end subroutine grid_comm_type_curv_2

  !===============================================================================
  subroutine grid_comm_type_curv(nx_,ny_)
  !===============================================================================
    !> Definition of MPI types for interblock grid communications
    !> * 2D Curvilinear version *
  !===============================================================================
    use mod_grid_directions
    implicit none
    ! ----------------------------------------------------------------------------
    integer, intent(in) :: nx_,ny_  ! sizes
    ! ----------------------------------------------------------------------------
    integer :: nex,ney ! extended sizes
    integer :: sign_i,sign_j
    integer :: memjump ! jump in memory between block starts
    integer(kind=MPI_ADDRESS_KIND) :: stride ! stride in memory between block starts
    integer :: type_base
    ! ----------------------------------------------------------------------------
    
    ! MPI types for intrablock grid comm
    ! ----------------------------------
    stride=(nx_+2*ngh)*sizeofreal

    ! Construction of type "face along x" (north or south)
    call MPI_TYPE_VECTOR(nx_+2*ngh,1,1,MPI_DOUBLE_PRECISION,type_base,info)
    call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride,type_base,type_x_facex,info)
    call MPI_TYPE_COMMIT(type_x_facex,info)

    ! Construction of type "face along y" (west or east)
    call MPI_TYPE_VECTOR(ngh,1,1,MPI_DOUBLE_PRECISION,type_base,info)
    call MPI_TYPE_CREATE_HVECTOR(ny_+2*ngh,1,stride,type_base,type_x_facey,info)
    call MPI_TYPE_COMMIT(type_x_facey,info)  
    
    ! Extended sizes (+ ghost cells)
    ! ------------------------------
    nex=nx_+2*ngh
    ney=ny_+2*ngh

    ! MPI-type construction for imin face "along y" (W: west)
    ! -------------------------------------------------------
    ! for imin/W, parallel dir is j and normal dir is i
    inW_bl=1
    inWs_bl=2
    inwr_bl=1-ngh
    ipW_bl=1-ngh  
    sign_i=1
    sign_j=1
    ! if rev parallel -> +j => -j
    if (is_rev_bl(1,1)) then
       sign_j=-1
       ipW_bl=ny_+ngh
    endif
    ! if rev  normal  -> +i => -i
    if (is_rev_bl(1,2)) then
       sign_i=-1
       inW_bl=ngh
       inWs_bl=ngh+1
       inwr_bl=0
    endif
    if (is_swapij_bl(1)) then
       ! if swap invert role of i and j
       memjump=sign_j*nex  ! j
       stride =sign_i*sizeofreal ! i
       call MPI_TYPE_VECTOR(ney,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride,type_base,type_x_faceW,info)
    else
       ! default
       memjump=sign_i ! i
       stride =sign_j*nex*sizeofreal ! j
       call MPI_TYPE_VECTOR(ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ney,1,stride,type_base,type_x_faceW,info)
    endif
    call MPI_TYPE_COMMIT(type_x_faceW,info)

    ! MPI-type construction for imax face "along y" (E: east)
    ! -------------------------------------------------------
    ! for imax/E, parallel dir is j and normal dir is i
    inE_bl=nx_+1
    inEs_bl=nx_-ngh
    inEr_bl=nx_+1
    ipE_bl=1-ngh
    sign_i=1
    sign_j=1
    ! if rev parallel -> +j => -j
    if (is_rev_bl(2,1)) then
       sign_j=-1
       ipE_bl=ny_+ngh
    endif
    ! if rev  normal  -> +i => -i
    if (is_rev_bl(2,2)) then
       sign_i=-1
       inE_bl=nx_+ngh
       inEs_bl=nx_-1
       inEr_bl=nx_+ngh
    endif
    if (is_swapij_bl(2)) then
       ! if swap invert role of i and j
       memjump=sign_j*nex  ! j
       stride =sign_i*sizeofreal ! i
       call MPI_TYPE_VECTOR(ney,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride,type_base,type_x_faceE,info)
    else
       ! default
       memjump=sign_i ! i
       stride =sign_j*nex*sizeofreal ! j
       call MPI_TYPE_VECTOR(ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ney,1,stride,type_base,type_x_faceE,info)
    endif
    call MPI_TYPE_COMMIT(type_x_faceE,info)

    ! MPI-type construction for jmin face "along x" (S: south)
    ! --------------------------------------------------------
    ! for jmin/S, parallel dir is i and normal dir is j
    ipS_bl=1-ngh
    inS_bl=1
    inSs_bl=2
    inSr_bl=1-ngh
    sign_i=1
    sign_j=1
    ! if rev parallel -> +i => -i
    if (is_rev_bl(3,1)) then
       sign_i=-1
       ipS_bl=nx_+ngh
    endif
    ! if rev  normal  -> +j => -j
    if (is_rev_bl(3,2)) then
       sign_j=-1
       inS_bl=ngh
       inSs_bl=ngh+1
       inSr_bl=0
    endif
    if (is_swapij_bl(3)) then
       ! if swap invert role of i and j
       memjump=sign_j*nex  ! j
       stride =sign_i*sizeofreal ! i
       call MPI_TYPE_VECTOR(ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(nex,1,stride,type_base,type_x_faceS,info)
    else
       ! default
       memjump=sign_i ! i
       stride =sign_j*nex*sizeofreal ! j
       call MPI_TYPE_VECTOR(nex,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride,type_base,type_x_faceS,info)
    endif
    call MPI_TYPE_COMMIT(type_x_faceS,info)

    ! MPI-type construction for jmax face "along x" (N: north)
    ! --------------------------------------------------------
    ! for jmax/N, parallel dir is i and normal dir is j
    ipN_bl=1-ngh
    inN_bl=ny_+1
    inNs_bl=ny_-ngh
    inNr_bl=ny_+1
    sign_i=1
    sign_j=1
    ! if rev parallel -> +i => -i
    if (is_rev_bl(4,1)) then
       sign_i=-1
       ipN_bl=nx_+ngh
    endif
    ! if rev  normal  -> +j => -j
    if (is_rev_bl(4,2)) then
       sign_j=-1
       inN_bl=ny_+ngh
       inNs_bl=ny_-1
       inNr_bl=ny_+ngh
    endif
    if (is_swapij_bl(4)) then
       ! if swap invert role of i and j
       memjump=sign_j*nex  ! j
       stride =sign_i*sizeofreal ! i
       call MPI_TYPE_VECTOR(ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(nex,1,stride,type_base,type_x_faceN,info)
    else
       ! default
       memjump=sign_i ! i
       stride =sign_j*nex*sizeofreal ! j
       call MPI_TYPE_VECTOR(nex,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride,type_base,type_x_faceN,info)
    endif
    call MPI_TYPE_COMMIT(type_x_faceN,info)
  
  end subroutine grid_comm_type_curv

  !===============================================================================
  subroutine grid_comm_type_curv3(nx_,ny_,nz_)
  !===============================================================================
    !> Definition of MPI types for interblock grid communications
    !> * 3D Curvilinear version *
  !===============================================================================
    use mod_grid_directions
    implicit none
    ! ----------------------------------------------------------------------------
    integer, intent(in) :: nx_,ny_,nz_  ! sizes
    ! ----------------------------------------------------------------------------
    integer :: nex,ney,nez ! extended sizes
    integer :: sign_i,sign_j
    integer :: memjump ! jump in memory between block starts
    ! strides in memory between block starts
    integer(kind=MPI_ADDRESS_KIND) :: stride,stride2
    integer :: type_base,type_edge
    ! ----------------------------------------------------------------------------

    ! Extended sizes (+ ghost cells)
    ! ------------------------------
    nex=nx_+2*ngh
    ney=ny_+2*ngh
    nez=nz_+2*ngh

    ! MPI types for intrablock grid comm
    ! ----------------------------------
    stride=nex*sizeofreal
    stride2=nex*ney*sizeofreal

    ! Construction of type "face along y" (west or east)
    call MPI_TYPE_VECTOR(ngh,1,1,MPI_DOUBLE_PRECISION,type_base,info)
    call MPI_TYPE_CREATE_HVECTOR(ney,1,stride,type_base,type_edge,info)
    call MPI_TYPE_CREATE_HVECTOR(nez,1,stride2,type_edge,type_x_facey,info)
    call MPI_TYPE_COMMIT(type_x_facey,info)

    ! Construction of type "face along x" (north or south)
    call MPI_TYPE_VECTOR(nex,1,1,MPI_DOUBLE_PRECISION,type_base,info)
    call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride,type_base,type_edge,info)
    call MPI_TYPE_CREATE_HVECTOR(nez,1,stride2,type_edge,type_x_facex,info)
    call MPI_TYPE_COMMIT(type_x_facex,info)

    ! Construction of type "face along z" (forward or backward)
    call MPI_TYPE_CREATE_HVECTOR(ney,1, stride,type_base,type_edge,info)
    call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride2,type_edge,type_x_facez,info)
    call MPI_TYPE_COMMIT(type_x_facez,info)

    ! MPI-type construction for imin face "along y" (W: west)
    ! -------------------------------------------------------
    ! for imin/W, parallel dir is j and normal dir is i
    inW_bl=1
    inWs_bl=2
    inwr_bl=1-ngh
    ipW_bl=1-ngh
    sign_i=1
    sign_j=1
    ! if rev parallel -> +j => -j
    if (is_rev_bl(1,1)) then
       sign_j=-1
       ipW_bl=ny_+ngh
    endif
    ! if rev  normal  -> +i => -i
    if (is_rev_bl(1,2)) then
       sign_i=-1
       inW_bl=ngh
       inWs_bl=ngh+1
       inwr_bl=0
    endif
    if (is_swapij_bl(1)) then
       ! if swap invert role of i and j
       memjump=sign_j*nex  ! j
       stride =sign_i*sizeofreal ! i
       call MPI_TYPE_VECTOR(ney,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride,type_base,type_edge,info)
    else
       ! default
       memjump=sign_i ! i
       stride =sign_j*nex*sizeofreal ! j
       call MPI_TYPE_VECTOR(ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ney,1,stride,type_base,type_edge,info)
    endif
    call MPI_TYPE_CREATE_HVECTOR(nez,1,stride2,type_edge,type_x_faceW,info)
    call MPI_TYPE_COMMIT(type_x_faceW,info)

    ! MPI-type construction for imax face "along y" (E: east)
    ! -------------------------------------------------------
    ! for imax/E, parallel dir is j and normal dir is i
    inE_bl=nx_+1
    inEs_bl=nx_-ngh
    inEr_bl=nx_+1
    ipE_bl=1-ngh
    sign_i=1
    sign_j=1
    ! if rev parallel -> +j => -j
    if (is_rev_bl(2,1)) then
       sign_j=-1
       ipE_bl=ny_+ngh
    endif
    ! if rev  normal  -> +i => -i
    if (is_rev_bl(2,2)) then
       sign_i=-1
       inE_bl=nx_+ngh
       inEs_bl=nx_-1
       inEr_bl=nx_+ngh
    endif
    if (is_swapij_bl(2)) then
       ! if swap invert role of i and j
       memjump=sign_j*nex  ! j
       stride =sign_i*sizeofreal ! i
       call MPI_TYPE_VECTOR(ney,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride,type_base,type_edge,info)
    else
       ! default
       memjump=sign_i ! i
       stride =sign_j*nex*sizeofreal ! j
       call MPI_TYPE_VECTOR(ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ney,1,stride,type_base,type_edge,info)
    endif
    call MPI_TYPE_CREATE_HVECTOR(nez,1,stride2,type_edge,type_x_faceE,info)
    call MPI_TYPE_COMMIT(type_x_faceE,info)

    ! MPI-type construction for jmin face "along x" (S: south)
    ! --------------------------------------------------------
    ! for jmin/S, parallel dir is i and normal dir is j
    ipS_bl=1-ngh
    inS_bl=1
    inSs_bl=2
    inSr_bl=1-ngh
    sign_i=1
    sign_j=1
    ! if rev parallel -> +i => -i
    if (is_rev_bl(3,1)) then
       sign_i=-1
       ipS_bl=nx_+ngh
    endif
    ! if rev  normal  -> +j => -j
    if (is_rev_bl(3,2)) then
       sign_j=-1
       inS_bl=ngh
       inSs_bl=ngh+1
       inSr_bl=0
    endif
    if (is_swapij_bl(3)) then
       ! if swap invert role of i and j
       memjump=sign_j*nex  ! j
       stride =sign_i*sizeofreal ! i
       call MPI_TYPE_VECTOR(ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(nex,1,stride,type_base,type_edge,info)
    else
       ! default
       memjump=sign_i ! i
       stride =sign_j*nex*sizeofreal ! j
       call MPI_TYPE_VECTOR(nex,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride,type_base,type_edge,info)
    endif
    call MPI_TYPE_CREATE_HVECTOR(nez,1,stride2,type_edge,type_x_faceS,info)
    call MPI_TYPE_COMMIT(type_x_faceS,info)

    ! MPI-type construction for jmax face "along x" (N: north)
    ! --------------------------------------------------------
    ! for jmax/N, parallel dir is i and normal dir is j
    ipN_bl=1-ngh
    inN_bl=ny_+1
    inNs_bl=ny_-ngh
    inNr_bl=ny_+1
    sign_i=1
    sign_j=1
    ! if rev parallel -> +i => -i
    if (is_rev_bl(4,1)) then
       sign_i=-1
       ipN_bl=nx_+ngh
    endif
    ! if rev  normal  -> +j => -j
    if (is_rev_bl(4,2)) then
       sign_j=-1
       inN_bl=ny_+ngh
       inEs_bl=ny_-1
       inEr_bl=ny_+ngh
    endif
    if (is_swapij_bl(4)) then
       ! if swap invert role of i and j
       memjump=sign_j*nex  ! j
       stride =sign_i*sizeofreal ! i
       call MPI_TYPE_VECTOR(ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(nex,1,stride,type_base,type_edge,info)
    else
       ! default
       memjump=sign_i ! i
       stride =sign_j*nex*sizeofreal ! j
       call MPI_TYPE_VECTOR(nex,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride,type_base,type_edge,info)
    endif
    call MPI_TYPE_CREATE_HVECTOR(nez,1,stride2,type_edge,type_x_faceN,info)
    call MPI_TYPE_COMMIT(type_x_faceN,info)

    stride=nex*sizeofreal
    call MPI_TYPE_VECTOR(nex,1,1,MPI_DOUBLE_PRECISION,type_base,info)

    ! MPI-type construction for kmin face "along z" (F: forward)
    ! --------------------------------------------------------
    call MPI_TYPE_CREATE_HVECTOR(ney,1, stride,type_base,type_edge,info)
    call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride2,type_edge,type_x_faceF,info)
    call MPI_TYPE_COMMIT(type_x_faceF,info)

    ! MPI-type construction for kmax face "along z" (B: backward)
    ! -----------------------------------------------------------
    call MPI_TYPE_CREATE_HVECTOR(ney,1, stride,type_base,type_edge,info)
    call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride2,type_edge,type_x_faceB,info)
    call MPI_TYPE_COMMIT(type_x_faceB,info)

  end subroutine grid_comm_type_curv3

  !===============================================================================
  subroutine free_grid_comm_type_cart
  !===============================================================================
    !> Free MPI types for grid communications - Cartesian version
  !===============================================================================
   implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------
    
    call MPI_TYPE_FREE(type_x,info)
    
  end subroutine free_grid_comm_type_cart

  !===============================================================================
  subroutine free_grid_comm_type_curv
  !===============================================================================
    !> Free MPI types for grid communications - 2D curvilinear version
  !===============================================================================
   implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------
    
    call MPI_TYPE_FREE(type_x_facex,info)
    call MPI_TYPE_FREE(type_x_facey,info)
    call MPI_TYPE_FREE(type_x_faceW,info)
    call MPI_TYPE_FREE(type_x_faceE,info)
    call MPI_TYPE_FREE(type_x_faceS,info)
    call MPI_TYPE_FREE(type_x_faceN,info)
    
  end subroutine free_grid_comm_type_curv

  !===============================================================================
  subroutine free_grid_comm_type_curv3
  !===============================================================================
    !> Free MPI types for grid communications - 3D curvilinear version
  !===============================================================================
   implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    call MPI_TYPE_FREE(type_x_facex,info)
    call MPI_TYPE_FREE(type_x_facey,info)
    call MPI_TYPE_FREE(type_x_facez,info)
    call MPI_TYPE_FREE(type_x_faceW,info)
    call MPI_TYPE_FREE(type_x_faceE,info)
    call MPI_TYPE_FREE(type_x_faceS,info)
    call MPI_TYPE_FREE(type_x_faceN,info)
    call MPI_TYPE_FREE(type_x_faceF,info)
    call MPI_TYPE_FREE(type_x_faceB,info)

  end subroutine free_grid_comm_type_curv3

  !===============================================================================
  subroutine grid_comm_interblock_cart(x_,nx_,dim)
  !===============================================================================
    !> Routine for interblock communications of global grid
    !> * Cartesian version *
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: dim  ! direction
    integer :: nx_  ! size
    real(wp), dimension(1-ngh:nx_+ngh), intent(inout)  :: x_  ! grid
    ! ----------------------------------------------------------------------------
    integer :: n1,n2 ! neighbors
    ! ----------------------------------------------------------------------------
    
    ! Define neighbors depending on dimension
    ! ---------------------------------------
    n1=2*dim
    n2=2*dim-1
    
    ! MPI SENDRECV
    ! ------------
    
    if(is_adjoint_block) then ! for adjoint block interfaces
    
      ! Send to neighbor 1 and reception from neighbor 2
      call MPI_SENDRECV(x_(nx_-ngh),1,type_x,neighbor_bl(n1),tag &
                       ,x_(-ngh+1) ,1,type_x,neighbor_bl(n2),tag,COMM_interblock,status,info) 
                       
      ! Send to neighbor 2 and reception from neighbor 1
      call MPI_SENDRECV(x_(2)    ,1,type_x,neighbor_bl(n2),tag &
                       ,x_(nx_+1),1,type_x,neighbor_bl(n1),tag,COMM_interblock,status,info)
                     
    else
    
      ! Send to neighbor 1 and reception from neighbor 2
      call MPI_SENDRECV(x_(nx_-ngh+1),1,type_x,neighbor_bl(n1),tag &
                       ,x_(-ngh+1)   ,1,type_x,neighbor_bl(n2),tag,COMM_interblock,status,info)  
      ! Send to neighbor 2 and reception from neighbor 1
      call MPI_SENDRECV(x_(1)    ,1,type_x,neighbor_bl(n2),tag &
                       ,x_(nx_+1),1,type_x,neighbor_bl(n1),tag,COMM_interblock,status,info)
        
    endif

  end subroutine grid_comm_interblock_cart

  !===============================================================================
  subroutine grid_comm_intrablock_cart(x_,nx_)
  !===============================================================================
    !> Routine for intrablock communications of global grid
    !> * Cartesian version *
  !===============================================================================
   implicit none
    ! ----------------------------------------------------------------------------
    integer :: nx_  ! size
    real(wp), dimension(1-ngh:nx_+ngh), intent(inout)  :: x_  ! grid
    ! ----------------------------------------------------------------------------
    integer :: ip
    ! ----------------------------------------------------------------------------

    nbl=nob(iproc)

    ! imin ghost points
    ! -----------------
    if (iproc.eq.iproc_leader(nbl)) then
       do ip=1,bl(nbl)%nproc-1
          call MPI_SEND(x_(-ngh+1),1,type_x,ip,tag,COMM_intrablock,info)
       enddo
    else
       call MPI_RECV(x_(-ngh+1),1,type_x,0,tag,COMM_intrablock,status,info)
    endif
    call MPI_BARRIER(COMM_global,info)

    ! imax ghost points
    ! -----------------
    if (iproc.eq.iproc_leader(nbl)) then
       do ip=1,bl(nbl)%nproc-1
          call MPI_SEND(x_(nx_+1),1,type_x,ip,tag,COMM_intrablock,info)
       enddo
    else
       call MPI_RECV(x_(nx_+1),1,type_x,0,tag,COMM_intrablock,status,info)
    endif
    call MPI_BARRIER(COMM_global,info)

  end subroutine grid_comm_intrablock_cart

  !===============================================================================
  subroutine grid_comm_interblock_curv(x_,nx_,ny_)
  !===============================================================================
    !> Routine for interblock communications of global grid
    !> * 2D Curvilinear version *
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer, intent(in) :: nx_,ny_  ! size
    real(wp), dimension(1-ngh:nx_+ngh,1-ngh:ny_+ngh), intent(inout) :: x_
    ! ---------------------------------------------------------------------------
    integer, parameter :: size_bl=8
    integer, dimension(size_bl) :: request_bl
    integer, dimension(MPI_STATUS_SIZE,size_bl) :: status_bl
    ! ---------------------------------------------------------------------------
    
    if (is_adjoint_block) then ! for adjoint block interfaces

       ! Send to neighbor W and reception from neighbor W
       call MPI_ISEND(x_(inWs_bl,ipW_bl),1,type_x_faceW,neighbor_bl(nW),tags_bl(nW),COMM_interblock,request_bl(1),info)
       call MPI_IRECV(x_(inWr_bl,ipW_bl),1,type_x_faceW,neighbor_bl(nW),tagr_bl(nW),COMM_interblock,request_bl(2),info)

       ! Send to neighbor E and reception from neighbor E
       call MPI_ISEND(x_(inEs_bl,ipE_bl),1,type_x_faceE,neighbor_bl(nE),tags_bl(nE),COMM_interblock,request_bl(3),info)
       call MPI_IRECV(x_(inEr_bl,ipE_bl),1,type_x_faceE,neighbor_bl(nE),tagr_bl(nE),COMM_interblock,request_bl(4),info)

       ! Send to neighbor S and reception from neighbor S
       call MPI_ISEND(x_(ipS_bl,inSs_bl),1,type_x_faceS,neighbor_bl(nS),tags_bl(nS),COMM_interblock,request_bl(5),info)
       call MPI_IRECV(x_(ipS_bl,inSr_bl),1,type_x_faceS,neighbor_bl(nS),tagr_bl(nS),COMM_interblock,request_bl(6),info)

       ! Send to neighbor N and reception from neighbor N
       call MPI_ISEND(x_(ipN_bl,inNs_bl),1,type_x_faceN,neighbor_bl(nN),tags_bl(nN),COMM_interblock,request_bl(7),info)
       call MPI_IRECV(x_(ipN_bl,inNr_bl),1,type_x_faceN,neighbor_bl(nN),tagr_bl(nN),COMM_interblock,request_bl(8),info)

    else

       ! Send to neighbor W and reception from neighbor W
       call MPI_ISEND(x_(inW_bl    ,ipW_bl),1,type_x_faceW,neighbor_bl(nW),tags_bl(nW),COMM_interblock,request_bl(1),info)
       call MPI_IRECV(x_(inW_bl-ngh,ipW_bl),1,type_x_faceW,neighbor_bl(nW),tagr_bl(nW),COMM_interblock,request_bl(2),info)

       ! Send to neighbor E and reception from neighbor E
       call MPI_ISEND(x_(inE_bl-ngh,ipE_bl),1,type_x_faceE,neighbor_bl(nE),tags_bl(nE),COMM_interblock,request_bl(3),info)
       call MPI_IRECV(x_(inE_bl    ,ipE_bl),1,type_x_faceE,neighbor_bl(nE),tagr_bl(nE),COMM_interblock,request_bl(4),info)

       ! Send to neighbor S and reception from neighbor S
       call MPI_ISEND(x_(ipS_bl,inS_bl    ),1,type_x_faceS,neighbor_bl(nS),tags_bl(nS),COMM_interblock,request_bl(5),info)
       call MPI_IRECV(x_(ipS_bl,inS_bl-ngh),1,type_x_faceS,neighbor_bl(nS),tagr_bl(nS),COMM_interblock,request_bl(6),info)

       ! Send to neighbor N and reception from neighbor N
       call MPI_ISEND(x_(ipN_bl,inN_bl-ngh),1,type_x_faceN,neighbor_bl(nN),tags_bl(nN),COMM_interblock,request_bl(7),info)
       call MPI_IRECV(x_(ipN_bl,inN_bl    ),1,type_x_faceN,neighbor_bl(nN),tagr_bl(nN),COMM_interblock,request_bl(8),info)

    endif
    
    call MPI_WAITALL(size_bl,request_bl,status_bl,info)
    
  end subroutine grid_comm_interblock_curv

  !===============================================================================
  subroutine grid_comm_interblock_curv_2(x_,nx_,ny_)
  !===============================================================================
    !> Routine for interblock communications of global grid
    !> * 2D Curvilinear version *
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer, intent(in) :: nx_,ny_  ! size
    real(wp), dimension(1-ngh:nx_+ngh,1-ngh:ny_+ngh), intent(inout) :: x_
    ! ---------------------------------------------------------------------------
    integer :: n
    integer, parameter :: size_bl=8
    integer, dimension(size_bl) :: request_bl
    integer, dimension(MPI_STATUS_SIZE,size_bl) :: status_bl
    ! ---------------------------------------------------------------------------

    do n=1,4
       call MPI_ISEND(x_(iis_bl(n),ijs_bl(n)),1,type_x_face(n), &
                      neighbor_bl(n),tags_bl(n),COMM_interblock,request_bl(2*n-1),info)
       call MPI_IRECV(x_(iir_bl(n),ijr_bl(n)),1,type_x_face(n), &
                      neighbor_bl(n),tagr_bl(n),COMM_interblock,request_bl(2*n),info)
    enddo

!!$    if (is_adjoint_block) then ! for adjoint block interfaces
!!$
!!$       ! Send to neighbor W and reception from neighbor W
!!$       call MPI_ISEND(x_(inWs_bl,ipW_bl),1,type_x_faceW,neighbor_bl(nW),tags_bl(nW),COMM_interblock,request_bl(1),info)
!!$       call MPI_IRECV(x_(inWr_bl,ipW_bl),1,type_x_faceW,neighbor_bl(nW),tagr_bl(nW),COMM_interblock,request_bl(2),info)
!!$
!!$       ! Send to neighbor E and reception from neighbor E
!!$       call MPI_ISEND(x_(inEs_bl,ipE_bl),1,type_x_faceE,neighbor_bl(nE),tags_bl(nE),COMM_interblock,request_bl(3),info)
!!$       call MPI_IRECV(x_(inEr_bl,ipE_bl),1,type_x_faceE,neighbor_bl(nE),tagr_bl(nE),COMM_interblock,request_bl(4),info)
!!$
!!$       ! Send to neighbor S and reception from neighbor S
!!$       call MPI_ISEND(x_(ipS_bl,inSs_bl),1,type_x_faceS,neighbor_bl(nS),tags_bl(nS),COMM_interblock,request_bl(5),info)
!!$       call MPI_IRECV(x_(ipS_bl,inSr_bl),1,type_x_faceS,neighbor_bl(nS),tagr_bl(nS),COMM_interblock,request_bl(6),info)
!!$
!!$       ! Send to neighbor N and reception from neighbor N
!!$       call MPI_ISEND(x_(ipN_bl,inNs_bl),1,type_x_faceN,neighbor_bl(nN),tags_bl(nN),COMM_interblock,request_bl(7),info)
!!$       call MPI_IRECV(x_(ipN_bl,inNr_bl),1,type_x_faceN,neighbor_bl(nN),tagr_bl(nN),COMM_interblock,request_bl(8),info)
!!$
!!$    else
!!$
!!$       ! Send to neighbor W and reception from neighbor W
!!$       call MPI_ISEND(x_(inW_bl    ,ipW_bl),1,type_x_faceW,neighbor_bl(nW),tags_bl(nW),COMM_interblock,request_bl(1),info)
!!$       call MPI_IRECV(x_(inW_bl-ngh,ipW_bl),1,type_x_faceW,neighbor_bl(nW),tagr_bl(nW),COMM_interblock,request_bl(2),info)
!!$
!!$       ! Send to neighbor E and reception from neighbor E
!!$       call MPI_ISEND(x_(inE_bl-ngh,ipE_bl),1,type_x_faceE,neighbor_bl(nE),tags_bl(nE),COMM_interblock,request_bl(3),info)
!!$       call MPI_IRECV(x_(inE_bl    ,ipE_bl),1,type_x_faceE,neighbor_bl(nE),tagr_bl(nE),COMM_interblock,request_bl(4),info)
!!$
!!$       ! Send to neighbor S and reception from neighbor S
!!$       call MPI_ISEND(x_(ipS_bl,inS_bl    ),1,type_x_faceS,neighbor_bl(nS),tags_bl(nS),COMM_interblock,request_bl(5),info)
!!$       call MPI_IRECV(x_(ipS_bl,inS_bl-ngh),1,type_x_faceS,neighbor_bl(nS),tagr_bl(nS),COMM_interblock,request_bl(6),info)
!!$
!!$       ! Send to neighbor N and reception from neighbor N
!!$       call MPI_ISEND(x_(ipN_bl,inN_bl-ngh),1,type_x_faceN,neighbor_bl(nN),tags_bl(nN),COMM_interblock,request_bl(7),info)
!!$       call MPI_IRECV(x_(ipN_bl,inN_bl    ),1,type_x_faceN,neighbor_bl(nN),tagr_bl(nN),COMM_interblock,request_bl(8),info)
!!$
!!$    endif

    call MPI_WAITALL(size_bl,request_bl,status_bl,info)

  end subroutine grid_comm_interblock_curv_2

  !===============================================================================
  subroutine grid_comm_intrablock_curv(x_,nx_,ny_)
  !===============================================================================
    !> Routine for intrablock communications of global grid
    !> * 2D Curvilinear version *
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer, intent(in) :: nx_,ny_  ! size
    real(wp), dimension(1-ngh:nx_+ngh,1-ngh:ny_+ngh), intent(inout)  :: x_  ! grid
    ! ----------------------------------------------------------------------------
    integer :: ip
    ! ----------------------------------------------------------------------------
 
    nbl=nob(iproc)

    ! imin ghost points
    ! -----------------
    if (iproc.eq.iproc_leader(nbl)) then       
       do ip=1,bl(nbl)%nproc-1
          call MPI_SEND(x_(-ngh+1,-ngh+1),1,type_x_facey,ip,tag,COMM_intrablock,info)
       enddo
    else
       call MPI_RECV(x_(-ngh+1,-ngh+1),1,type_x_facey,0,tag,COMM_intrablock,status,info)
    endif
    call MPI_BARRIER(COMM_global,info)

    ! imax ghost points
    ! -----------------
    if (iproc.eq.iproc_leader(nbl)) then
       do ip=1,bl(nbl)%nproc-1
          call MPI_SEND(x_(nx_+1,-ngh+1),1,type_x_facey,ip,tag,COMM_intrablock,info)
       enddo
    else
       call MPI_RECV(x_(nx_+1,-ngh+1),1,type_x_facey,0,tag,COMM_intrablock,status,info)
    endif
    call MPI_BARRIER(COMM_global,info)

    ! jmin ghost points
    ! -----------------
    if (iproc.eq.iproc_leader(nbl)) then       
       do ip=1,bl(nbl)%nproc-1
          call MPI_SEND(x_(-ngh+1,-ngh+1),1,type_x_facex,ip,tag,COMM_intrablock,info)
       enddo
    else
       call MPI_RECV(x_(-ngh+1,-ngh+1),1,type_x_facex,0,tag,COMM_intrablock,status,info)
    endif
    call MPI_BARRIER(COMM_global,info)

    ! jmax ghost points
    ! -----------------
    if (iproc.eq.iproc_leader(nbl)) then
       do ip=1,bl(nbl)%nproc-1
          call MPI_SEND(x_(-ngh+1,ny_+1),1,type_x_facex,ip,tag,COMM_intrablock,info)
       enddo
    else
       call MPI_RECV(x_(-ngh+1,ny_+1),1,type_x_facex,0,tag,COMM_intrablock,status,info)
    endif
    call MPI_BARRIER(COMM_global,info)

  end subroutine grid_comm_intrablock_curv

  !===============================================================================
  subroutine grid_comm_interblock_curv3(x_,nx_,ny_,nz_)
  !===============================================================================
    !> Routine for interblock communications of global grid
    !> * 3D Curvilinear version *
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer, intent(in) :: nx_,ny_,nz_  ! size
    real(wp), dimension(1-ngh:nx_+ngh,1-ngh:ny_+ngh,1-ngh:nz_+ngh), intent(inout) :: x_
    ! ---------------------------------------------------------------------------
    !integer, parameter :: size_bl=12
    integer, parameter :: size_bl=8
    integer, dimension(size_bl) :: request_bl
    integer, dimension(MPI_STATUS_SIZE,size_bl) :: status_bl
    ! ---------------------------------------------------------------------------

    if (is_adjoint_block) then ! for adjoint block interfaces

       ! Send to neighbor W and reception from neighbor W
       call MPI_ISEND(x_(inWs_bl,ipW_bl,1-ngh),1,type_x_faceW,neighbor_bl(nW),tags_bl(nW),COMM_interblock,request_bl(1),info)
       call MPI_IRECV(x_(inWr_bl,ipW_bl,1-ngh),1,type_x_faceW,neighbor_bl(nW),tagr_bl(nW),COMM_interblock,request_bl(2),info)

       ! Send to neighbor E and reception from neighbor E
       call MPI_ISEND(x_(inEs_bl,ipE_bl,1-ngh),1,type_x_faceE,neighbor_bl(nE),tags_bl(nE),COMM_interblock,request_bl(3),info)
       call MPI_IRECV(x_(inEr_bl,ipE_bl,1-ngh),1,type_x_faceE,neighbor_bl(nE),tagr_bl(nE),COMM_interblock,request_bl(4),info)

       ! Send to neighbor S and reception from neighbor S
       call MPI_ISEND(x_(ipS_bl,inSs_bl,1-ngh),1,type_x_faceS,neighbor_bl(nS),tags_bl(nS),COMM_interblock,request_bl(5),info)
       call MPI_IRECV(x_(ipS_bl,inSr_bl,1-ngh),1,type_x_faceS,neighbor_bl(nS),tagr_bl(nS),COMM_interblock,request_bl(6),info)

       ! Send to neighbor N and reception from neighbor N
       call MPI_ISEND(x_(ipN_bl,inNs_bl,1-ngh),1,type_x_faceN,neighbor_bl(nN),tags_bl(nN),COMM_interblock,request_bl(7),info)
       call MPI_IRECV(x_(ipN_bl,inNr_bl,1-ngh),1,type_x_faceN,neighbor_bl(nN),tagr_bl(nN),COMM_interblock,request_bl(8),info)

    else

!!$       ! Send to neighbor W and reception from neighbor W
!!$       call MPI_ISEND(x_(inW_bl    ,ipW_bl,1-ngh),1,type_x_faceW,neighbor_bl(nW),tags_bl(nW),COMM_interblock,request_bl(1),info)
!!$       call MPI_IRECV(x_(inW_bl-ngh,ipW_bl,1-ngh),1,type_x_faceW,neighbor_bl(nW),tagr_bl(nW),COMM_interblock,request_bl(2),info)
!!$
!!$       ! Send to neighbor E and reception from neighbor E
!!$       call MPI_ISEND(x_(inE_bl-ngh,ipE_bl,1-ngh),1,type_x_faceE,neighbor_bl(nE),tags_bl(nE),COMM_interblock,request_bl(3),info)
!!$       call MPI_IRECV(x_(inE_bl    ,ipE_bl,1-ngh),1,type_x_faceE,neighbor_bl(nE),tagr_bl(nE),COMM_interblock,request_bl(4),info)
!!$
!!$       ! Send to neighbor S and reception from neighbor S
!!$       call MPI_ISEND(x_(ipS_bl,inS_bl    ,1-ngh),1,type_x_faceS,neighbor_bl(nS),tags_bl(nS),COMM_interblock,request_bl(5),info)
!!$       call MPI_IRECV(x_(ipS_bl,inS_bl-ngh,1-ngh),1,type_x_faceS,neighbor_bl(nS),tagr_bl(nS),COMM_interblock,request_bl(6),info)
!!$
!!$       ! Send to neighbor N and reception from neighbor N
!!$       call MPI_ISEND(x_(ipN_bl,inN_bl-ngh,1-ngh),1,type_x_faceN,neighbor_bl(nN),tags_bl(nN),COMM_interblock,request_bl(7),info)
!!$       call MPI_IRECV(x_(ipN_bl,inN_bl    ,1-ngh),1,type_x_faceN,neighbor_bl(nN),tagr_bl(nN),COMM_interblock,request_bl(8),info)

       ! Send to neighbor S and reception from neighbor S
       call MPI_ISEND(x_(ipS_bl,inS_bl    ,1-ngh),1,type_x_faceS,neighbor_bl(nS),tags_bl(nS),COMM_interblock,request_bl(1),info)
       call MPI_IRECV(x_(ipS_bl,inS_bl-ngh,1-ngh),1,type_x_faceS,neighbor_bl(nS),tagr_bl(nS),COMM_interblock,request_bl(2),info)

       ! Send to neighbor N and reception from neighbor N
       call MPI_ISEND(x_(ipN_bl,inN_bl-ngh,1-ngh),1,type_x_faceN,neighbor_bl(nN),tags_bl(nN),COMM_interblock,request_bl(3),info)
       call MPI_IRECV(x_(ipN_bl,inN_bl    ,1-ngh),1,type_x_faceN,neighbor_bl(nN),tagr_bl(nN),COMM_interblock,request_bl(4),info)

!!$       ! Send to neighbor S and reception from neighbor S
!!$       call MPI_ISEND(x_(ipS_bl,inS_bl    ,1),1,type_x_faceS,neighbor_bl(nS),tags_bl(nS),COMM_interblock,request_bl(1),info)
!!$       call MPI_IRECV(x_(ipS_bl,inS_bl-ngh,1),1,type_x_faceS,neighbor_bl(nS),tagr_bl(nS),COMM_interblock,request_bl(2),info)
!!$
!!$       ! Send to neighbor N and reception from neighbor N
!!$       call MPI_ISEND(x_(ipN_bl,inN_bl-ngh,1),1,type_x_faceN,neighbor_bl(nN),tags_bl(nN),COMM_interblock,request_bl(3),info)
!!$       call MPI_IRECV(x_(ipN_bl,inN_bl    ,1),1,type_x_faceN,neighbor_bl(nN),tagr_bl(nN),COMM_interblock,request_bl(4),info)

!!$       ! Send to neighbor W and reception from neighbor W
!!$       call MPI_ISEND(x_(inW_bl    ,ipW_bl,1),1,type_x_faceW,neighbor_bl(nW),tags_bl(nW),COMM_interblock,request_bl(1),info)
!!$       call MPI_IRECV(x_(inW_bl-ngh,ipW_bl,1),1,type_x_faceW,neighbor_bl(nW),tagr_bl(nW),COMM_interblock,request_bl(2),info)
!!$
!!$       ! Send to neighbor E and reception from neighbor E
!!$       call MPI_ISEND(x_(inE_bl-ngh,ipE_bl,1),1,type_x_faceE,neighbor_bl(nE),tags_bl(nE),COMM_interblock,request_bl(3),info)
!!$       call MPI_IRECV(x_(inE_bl    ,ipE_bl,1),1,type_x_faceE,neighbor_bl(nE),tagr_bl(nE),COMM_interblock,request_bl(4),info)
!!$
!!$       ! Send to neighbor S and reception from neighbor S
!!$       call MPI_ISEND(x_(ipS_bl,inS_bl    ,1),1,type_x_faceS,neighbor_bl(nS),tags_bl(nS),COMM_interblock,request_bl(5),info)
!!$       call MPI_IRECV(x_(ipS_bl,inS_bl-ngh,1),1,type_x_faceS,neighbor_bl(nS),tagr_bl(nS),COMM_interblock,request_bl(6),info)
!!$
!!$       ! Send to neighbor N and reception from neighbor N
!!$       call MPI_ISEND(x_(ipN_bl,inN_bl-ngh,1),1,type_x_faceN,neighbor_bl(nN),tags_bl(nN),COMM_interblock,request_bl(7),info)
!!$       call MPI_IRECV(x_(ipN_bl,inN_bl    ,1),1,type_x_faceN,neighbor_bl(nN),tagr_bl(nN),COMM_interblock,request_bl(8),info)

    endif

!!$    ! Send to neighbor F and reception from neighbor F
!!$    call MPI_ISEND(x_(1-ngh,1-ngh,    1),1,type_x_faceF,neighbor_bl(nF),tags_bl(nF),COMM_interblock,request_bl(9),info)
!!$    call MPI_IRECV(x_(1-ngh,1-ngh,1-ngh),1,type_x_faceF,neighbor_bl(nF),tagr_bl(nF),COMM_interblock,request_bl(10),info)
!!$
!!$    ! Send to neighbor B and reception from neighbor B
!!$    call MPI_ISEND(x_(1-ngh,1-ngh,nz_-ngh+1),1,type_x_faceB,neighbor_bl(nB),tags_bl(nB),COMM_interblock,request_bl(11),info)
!!$    call MPI_IRECV(x_(1-ngh,1-ngh,    nz_+1),1,type_x_faceB,neighbor_bl(nB),tagr_bl(nB),COMM_interblock,request_bl(12),info)

    ! Send to neighbor F and reception from neighbor F
    call MPI_ISEND(x_(1-ngh,1-ngh,    1),1,type_x_faceF,neighbor_bl(nF),tags_bl(nF),COMM_interblock,request_bl(5),info)
    call MPI_IRECV(x_(1-ngh,1-ngh,1-ngh),1,type_x_faceF,neighbor_bl(nF),tagr_bl(nF),COMM_interblock,request_bl(6),info)

    ! Send to neighbor B and reception from neighbor B
    call MPI_ISEND(x_(1-ngh,1-ngh,nz_-ngh+1),1,type_x_faceB,neighbor_bl(nB),tags_bl(nB),COMM_interblock,request_bl(7),info)
    call MPI_IRECV(x_(1-ngh,1-ngh,    nz_+1),1,type_x_faceB,neighbor_bl(nB),tagr_bl(nB),COMM_interblock,request_bl(8),info)

!!$    ! Send to neighbor F and reception from neighbor F
!!$    call MPI_ISEND(x_(1,1,    1),1,type_x_faceF,neighbor_bl(nF),tags_bl(nF),COMM_interblock,request_bl(5),info)
!!$    call MPI_IRECV(x_(1,1,1-ngh),1,type_x_faceF,neighbor_bl(nF),tagr_bl(nF),COMM_interblock,request_bl(6),info)
!!$
!!$    ! Send to neighbor B and reception from neighbor B
!!$    call MPI_ISEND(x_(1,1,nz_-ngh+1),1,type_x_faceB,neighbor_bl(nB),tags_bl(nB),COMM_interblock,request_bl(7),info)
!!$    call MPI_IRECV(x_(1,1,    nz_+1),1,type_x_faceB,neighbor_bl(nB),tagr_bl(nB),COMM_interblock,request_bl(8),info)

!!$    ! Send to neighbor F and reception from neighbor F
!!$    call MPI_ISEND(x_(1,1,    1),1,type_x_faceF,neighbor_bl(nF),tags_bl(nF),COMM_interblock,request_bl(9),info)
!!$    call MPI_IRECV(x_(1,1,1-ngh),1,type_x_faceF,neighbor_bl(nF),tagr_bl(nF),COMM_interblock,request_bl(10),info)
!!$
!!$    ! Send to neighbor B and reception from neighbor B
!!$    call MPI_ISEND(x_(1,1,nz_-ngh+1),1,type_x_faceB,neighbor_bl(nB),tags_bl(nB),COMM_interblock,request_bl(11),info)
!!$    call MPI_IRECV(x_(1,1,    nz_+1),1,type_x_faceB,neighbor_bl(nB),tagr_bl(nB),COMM_interblock,request_bl(12),info)

    call MPI_WAITALL(size_bl,request_bl,status_bl,info)

  end subroutine grid_comm_interblock_curv3

  !===============================================================================
  subroutine grid_comm_interblock_curv3_2(x_,nx_,ny_,nz_)
  !===============================================================================
    !> Routine for interblock communications of global grid
    !> * 3D Curvilinear version *
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer, intent(in) :: nx_,ny_,nz_  ! size
    real(wp), dimension(1-ngh:nx_+ngh,1-ngh:ny_+ngh,1-ngh:nz_+ngh), intent(inout) :: x_
    ! ---------------------------------------------------------------------------
    integer :: n
    integer, parameter :: size_bl=12
    integer, dimension(size_bl) :: request_bl
    integer, dimension(MPI_STATUS_SIZE,size_bl) :: status_bl
    ! ---------------------------------------------------------------------------

!!$    if (iproc==0) then
!!$       n=3
!!$       call MPI_ISEND(x_(iis_bl(n),ijs_bl(n),iks_bl(n)),1,type_x_face(n), &
!!$                      neighbor_bl(n),tags_bl(n),COMM_interblock,request_bl(1),info)
!!$       call MPI_IRECV(x_(iir_bl(n),ijr_bl(n),ikr_bl(n)),1,type_x_face(n), &
!!$                      neighbor_bl(n),tagr_bl(n),COMM_interblock,request_bl(2),info)
!!$    endif
!!$
!!$    if (iproc==1) then
!!$       n=4
!!$       call MPI_ISEND(x_(iis_bl(n),ijs_bl(n),iks_bl(n)),1,type_x_face(n), &
!!$                      neighbor_bl(n),tags_bl(n),COMM_interblock,request_bl(1),info)
!!$       call MPI_IRECV(x_(iir_bl(n),ijr_bl(n),ikr_bl(n)),1,type_x_face(n), &
!!$                      neighbor_bl(n),tagr_bl(n),COMM_interblock,request_bl(2),info)
!!$    endif

    do n=1,6

!!$       if ((neighbor(n)>=0)) then
!!$          print *,iproc,'ind send',iis_bl(n),ijs_bl(n),iks_bl(n),tags_bl(n)
!!$          print *,iproc,'ind recv',iir_bl(n),ijr_bl(n),ikr_bl(n),tagr_bl(n)
!!$       endif

       call MPI_ISEND(x_(iis_bl(n),ijs_bl(n),iks_bl(n)),1,type_x_face(n), &
                      neighbor_bl(n),tags_bl(n),COMM_interblock,request_bl(2*n-1),info)
       call MPI_IRECV(x_(iir_bl(n),ijr_bl(n),ikr_bl(n)),1,type_x_face(n), &
                      neighbor_bl(n),tagr_bl(n),COMM_interblock,request_bl(2*n),info)
    enddo


!!$    if (is_adjoint_block) then ! for adjoint block interfaces
!!$
!!$       ! Send to neighbor W and reception from neighbor W
!!$       call MPI_ISEND(x_(inWs_bl,ipW_bl,1-ngh),1,type_x_faceW,neighbor_bl(nW),tags_bl(nW),COMM_interblock,request_bl(1),info)
!!$       call MPI_IRECV(x_(inWr_bl,ipW_bl,1-ngh),1,type_x_faceW,neighbor_bl(nW),tagr_bl(nW),COMM_interblock,request_bl(2),info)
!!$
!!$       ! Send to neighbor E and reception from neighbor E
!!$       call MPI_ISEND(x_(inEs_bl,ipE_bl,1-ngh),1,type_x_faceE,neighbor_bl(nE),tags_bl(nE),COMM_interblock,request_bl(3),info)
!!$       call MPI_IRECV(x_(inEr_bl,ipE_bl,1-ngh),1,type_x_faceE,neighbor_bl(nE),tagr_bl(nE),COMM_interblock,request_bl(4),info)
!!$
!!$       ! Send to neighbor S and reception from neighbor S
!!$       call MPI_ISEND(x_(ipS_bl,inSs_bl,1-ngh),1,type_x_faceS,neighbor_bl(nS),tags_bl(nS),COMM_interblock,request_bl(5),info)
!!$       call MPI_IRECV(x_(ipS_bl,inSr_bl,1-ngh),1,type_x_faceS,neighbor_bl(nS),tagr_bl(nS),COMM_interblock,request_bl(6),info)
!!$
!!$       ! Send to neighbor N and reception from neighbor N
!!$       call MPI_ISEND(x_(ipN_bl,inNs_bl,1-ngh),1,type_x_faceN,neighbor_bl(nN),tags_bl(nN),COMM_interblock,request_bl(7),info)
!!$       call MPI_IRECV(x_(ipN_bl,inNr_bl,1-ngh),1,type_x_faceN,neighbor_bl(nN),tagr_bl(nN),COMM_interblock,request_bl(8),info)
!!$
!!$    else
!!$
!!$       ! Send to neighbor W and reception from neighbor W
!!$       call MPI_ISEND(x_(inW_bl    ,ipW_bl,1-ngh),1,type_x_faceW,neighbor_bl(nW),tags_bl(nW),COMM_interblock,request_bl(1),info)
!!$       call MPI_IRECV(x_(inW_bl-ngh,ipW_bl,1-ngh),1,type_x_faceW,neighbor_bl(nW),tagr_bl(nW),COMM_interblock,request_bl(2),info)
!!$
!!$       ! Send to neighbor E and reception from neighbor E
!!$       call MPI_ISEND(x_(inE_bl-ngh,ipE_bl,1-ngh),1,type_x_faceE,neighbor_bl(nE),tags_bl(nE),COMM_interblock,request_bl(3),info)
!!$       call MPI_IRECV(x_(inE_bl    ,ipE_bl,1-ngh),1,type_x_faceE,neighbor_bl(nE),tagr_bl(nE),COMM_interblock,request_bl(4),info)
!!$
!!$       ! Send to neighbor S and reception from neighbor S
!!$       call MPI_ISEND(x_(ipS_bl,inS_bl    ,1-ngh),1,type_x_faceS,neighbor_bl(nS),tags_bl(nS),COMM_interblock,request_bl(5),info)
!!$       call MPI_IRECV(x_(ipS_bl,inS_bl-ngh,1-ngh),1,type_x_faceS,neighbor_bl(nS),tagr_bl(nS),COMM_interblock,request_bl(6),info)
!!$
!!$       ! Send to neighbor N and reception from neighbor N
!!$       call MPI_ISEND(x_(ipN_bl,inN_bl-ngh,1-ngh),1,type_x_faceN,neighbor_bl(nN),tags_bl(nN),COMM_interblock,request_bl(7),info)
!!$       call MPI_IRECV(x_(ipN_bl,inN_bl    ,1-ngh),1,type_x_faceN,neighbor_bl(nN),tagr_bl(nN),COMM_interblock,request_bl(8),info)
!!$
!!$    endif
!!$
!!$    ! Send to neighbor F and reception from neighbor F
!!$    call MPI_ISEND(x_(1-ngh,1-ngh,    1),1,type_x_faceF,neighbor_bl(nF),tags_bl(nF),COMM_interblock,request_bl(9),info)
!!$    call MPI_IRECV(x_(1-ngh,1-ngh,1-ngh),1,type_x_faceF,neighbor_bl(nF),tagr_bl(nF),COMM_interblock,request_bl(10),info)
!!$
!!$    ! Send to neighbor B and reception from neighbor B
!!$    call MPI_ISEND(x_(1-ngh,1-ngh,nz_-ngh+1),1,type_x_faceB,neighbor_bl(nB),tags_bl(nB),COMM_interblock,request_bl(11),info)
!!$    call MPI_IRECV(x_(1-ngh,1-ngh,    nz_+1),1,type_x_faceB,neighbor_bl(nB),tagr_bl(nB),COMM_interblock,request_bl(12),info)

    call MPI_WAITALL(size_bl,request_bl,status_bl,info)

  end subroutine grid_comm_interblock_curv3_2

  !===============================================================================
  subroutine grid_comm_intrablock_curv3(x_,nx_,ny_,nz_)
  !===============================================================================
    !> Routine for intrablock communications of global grid
    !> * 3D Curvilinear version *
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer, intent(in) :: nx_,ny_,nz_  ! size
    real(wp), dimension(1-ngh:nx_+ngh,1-ngh:ny_+ngh,1-ngh:nz_+ngh), intent(inout)  :: x_  ! grid
    ! ----------------------------------------------------------------------------
    integer :: ip
    ! ----------------------------------------------------------------------------

    nbl=nob(iproc)

    ! imin ghost points
    ! -----------------
    if (iproc.eq.iproc_leader(nbl)) then
       do ip=1,bl(nbl)%nproc-1
          call MPI_SEND(x_(-ngh+1,-ngh+1,-ngh+1),1,type_x_facey,ip,tag,COMM_intrablock,info)
       enddo
    else
       call MPI_RECV(x_(-ngh+1,-ngh+1,-ngh+1),1,type_x_facey,0,tag,COMM_intrablock,status,info)
    endif
    call MPI_BARRIER(COMM_global,info)

    ! imax ghost points
    ! -----------------
    if (iproc.eq.iproc_leader(nbl)) then
       do ip=1,bl(nbl)%nproc-1
          call MPI_SEND(x_(nx_+1,-ngh+1,-ngh+1),1,type_x_facey,ip,tag,COMM_intrablock,info)
       enddo
    else
       call MPI_RECV(x_(nx_+1,-ngh+1,-ngh+1),1,type_x_facey,0,tag,COMM_intrablock,status,info)
    endif
    call MPI_BARRIER(COMM_global,info)

    ! jmin ghost points
    ! -----------------
    if (iproc.eq.iproc_leader(nbl)) then
       do ip=1,bl(nbl)%nproc-1
          call MPI_SEND(x_(-ngh+1,-ngh+1,-ngh+1),1,type_x_facex,ip,tag,COMM_intrablock,info)
       enddo
    else
       call MPI_RECV(x_(-ngh+1,-ngh+1,-ngh+1),1,type_x_facex,0,tag,COMM_intrablock,status,info)
    endif
    call MPI_BARRIER(COMM_global,info)

    ! jmax ghost points
    ! -----------------
    if (iproc.eq.iproc_leader(nbl)) then
       do ip=1,bl(nbl)%nproc-1
          call MPI_SEND(x_(-ngh+1,ny_+1,-ngh+1),1,type_x_facex,ip,tag,COMM_intrablock,info)
       enddo
    else
       call MPI_RECV(x_(-ngh+1,ny_+1,-ngh+1),1,type_x_facex,0,tag,COMM_intrablock,status,info)
    endif
    call MPI_BARRIER(COMM_global,info)

    ! kmin ghost points
    ! -----------------
    if (iproc.eq.iproc_leader(nbl)) then
       do ip=1,bl(nbl)%nproc-1
          call MPI_SEND(x_(-ngh+1,-ngh+1,-ngh+1),1,type_x_facez,ip,tag,COMM_intrablock,info)
       enddo
    else
       call MPI_RECV(x_(-ngh+1,-ngh+1,-ngh+1),1,type_x_facez,0,tag,COMM_intrablock,status,info)
    endif
    call MPI_BARRIER(COMM_global,info)

    ! kmax ghost points
    ! -----------------
    if (iproc.eq.iproc_leader(nbl)) then
       do ip=1,bl(nbl)%nproc-1
          call MPI_SEND(x_(-ngh+1,-ngh+1,nz_+1),1,type_x_facez,ip,tag,COMM_intrablock,info)
       enddo
    else
       call MPI_RECV(x_(-ngh+1,-ngh+1,nz_+1),1,type_x_facez,0,tag,COMM_intrablock,status,info)
    endif
    call MPI_BARRIER(COMM_global,info)

  end subroutine grid_comm_intrablock_curv3

end module mod_grid_comm
