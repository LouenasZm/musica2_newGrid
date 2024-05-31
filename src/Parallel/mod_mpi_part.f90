!===============================================================================
module mod_mpi_part
!===============================================================================
  !> Module for MPI parallelization
  !> 1. partitioning 
  !> 2. intrablock/interblock comm for grid partioning
  !> 3. sub-communicators for slicing
  !> 4. definition of neighbors for communications
  !> 5. MPI types for communications
  !> 6. further MPI constants for communications
!===============================================================================
  use mod_mpi
  implicit none
  ! ----------------------------------------------------------------------------
  !  1. partitioning
  !     ============
  logical, dimension(3) :: periods
  integer, dimension(:), allocatable :: coordx,coordy,coordz ! <- used in IO NOT DEFINED??????????????????,
  integer, dimension(:), allocatable :: coord1,coord2,coord3
  ! ----------------------------------------------------------------------------
  !  2. intrablock/interblock communicators for grid partioning
  !     =======================================================
  !     [used by Grid/grid_comm.f90]
  integer :: COMM_interblock,COMM_intrablock
  ! receive [tags] and send [tags] tags for interblock communications
  integer, dimension(6) :: tagr_bl,tags_bl
  ! interblock neighborhood
  integer, dimension(6) :: neighbor_bl
  ! detection of double corners for reconstruction block edges
  logical, dimension(4) :: dble_corner
  ! ----------------------------------------------------------------------------
  !  3. sub-communicators for slicing [in IO routines]
  !     =============================
  integer :: COMMX,COMMY,COMMZ,nprocz
  integer :: iprocxy,nprocxy,COMMXY
  integer :: iprocxz,nprocxz,COMMXZ
  integer :: iprocyz,nprocyz,COMMYZ
!!$  integer, dimension(:,:), allocatable :: ss_comm
!!$  logical, dimension(:), allocatable :: ixy,iyz
  ! ----------------------------------------------------------------------------
  !  4. definition of neighbors for communications
  !     ==========================================
  ! neighborhood for faces
  integer, parameter :: nW=1,nE=2,nS=3,nN=4,nF=5,nB=6
  integer, dimension(6) :: neighbor
  ! neighborhood for z-edges TO BE CHANGED [no more used for double derivatives BUT appears in Interpolation/mod_flow_o]
  integer, parameter :: nNW=1,nNE=2,nSE=3,nSW=4
  integer, dimension(4) :: neighbor2
  ! ----------------------------------------------------------------------------
  !  5. for MPI communications
  !     ======================
  ! receive [tags] and send [tags] tags for communications
  integer, dimension(6) :: tagr,tags
  
  ! indices to start communications in the direction [p]arallel to the face
  integer :: ipWi,ipEi,ipSi,ipNi
  ! indices to start communications in the direction [n]ormal to the face
  integer :: inWi,inEi,inSi,inNi
  ! indices to start send routines in the direction [n]ormal to the face
  integer :: inWsi,inEsi,inSsi,inNsi
  ! indices to start receive routines in the direction [n]ormal to the face
  integer :: inWri,inEri,inSri,inNri

  ! indices to start communications in the direction [p]arallel to the face
  integer :: ipW,ipE,ipS,ipN
  ! indices to start communications in the direction [n]ormal to the face
  integer :: inW,inE,inS,inN
  ! indices to start communications in the direction [n]ormal to the face (viscous terms)
  integer :: inW_v,inE_v,inS_v,inN_v
  !--------for the cases of adjoint block interfaces
  ! indices to start send routines in the direction [n]ormal to the face
  integer :: inWs,inEs,inSs,inNs
  ! indices to start receive routines in the direction [n]ormal to the face
  integer :: inWr,inEr,inSr,inNr
  ! indices to start send routines in the direction [n]ormal to the face (viscous terms)
  integer :: inWs_v,inEs_v,inSs_v,inNs_v
  ! indices to start receive routines in the direction [n]ormal to the face (viscous terms)
  integer :: inWr_v,inEr_v,inSr_v,inNr_v
  !-----------------------------------------------------------------------------
  ! indices to start communications in the direction [p]arallel for face+edges
  integer :: ipW_e,ipE_e,ipS_e,ipN_e
  ! indices to start communications in the direction [n]ormal for face+edges
  integer :: inW_e,inE_e,inS_e,inN_e
  
  ! indices to start send routines in the direction [n]ormal to the face [one-sided; variables]
  integer :: in1s,in2s,in3s,in4s
  ! indices to start send routines in the direction [n]ormal to the face [one-sided; increments]
  integer :: ii1s,ii2s,ii3s,ii4s,ii5s,ii6s
  ! ----------------------------------------------------------------------------
  !  6. further MPI constants for communications
  !     ========================================
  !     [used by Grid/grid_comm.f90,grid_comm_metrics.f90]
  !     [used by Parallel/mod_comm...]
  integer :: sizeofreal ! size of MPI real numbers
  integer, parameter :: tag=100 ! generic tag
  integer, dimension(MPI_STATUS_SIZE) :: status
  ! ----------------------------------------------------------------------------
  !  7. definiton of double corner type in case of adjoint block interfaces   NOTE: only for 5-point corners
  !     ===================================================================
  logical :: is_adjoint_block
  type double_corner
    integer :: dc_index(2),dc_blcs(5),face_index(5,6),face_blcs(5,2)
    logical :: dc_exist,face_exist(5)
!     real(wp), dimension(:,:), allocatable :: dc_flow
  end type double_corner
  type(double_corner), dimension(:), allocatable :: dble_crnr
  integer :: dc_blocks(4,5),dc_blocks_tmp(4,5),dc_tmp(5),no_dble_crnr
  integer, dimension(:,:,:), allocatable :: dc_blocks_gb_tmp
  integer, dimension(:,:), allocatable :: dc_blocks_gb,dc_COMM
  integer, dimension(:), allocatable :: dc_side,dc_side_tmp
  logical, dimension(:), allocatable :: dc_mask
  ! ----------------------------------------------------------------------------
  !  8. switch for one-sided or two-sided communication routines
  !     ========================================================
  logical :: is_two_sided_comm
  ! ----------------------------------------------------------------------------
  
contains

  !===============================================================
  subroutine mpi_set_dim
  !===============================================================
    !> Definition of MPI partitioning
  !===============================================================
    use mod_block
    use mod_constant ! <- for is_2d ??
    use mod_flow
    use warnstop
    implicit none
    ! ------------------------------------------------------------
    integer :: nbl
    integer :: cpt1,cpt2
    ! ------------------------------------------------------------
    ! partitioning
    integer :: ndom_tot ! number of proc per block
    ! ------------------------------------------------------------

    if (iproc==0) print *,'MPI Partitioning: set dimensions'
    !if (iproc==0) print *,'----------------'
    
    ! Check number of procs (block version)
    ! =====================
    ! sum numbers of proc per block
    ndom_tot=0
    do nbl=1,nbloc
       ndom_tot=ndom_tot+bl(nbl)%nproc
    enddo
    ! check if the total number of proc is the same as nproc (current number from mpirun -n nproc ...)
    if (ndom_tot.ne.nproc) then
       if (iproc==0) print 98,ndom_tot
       call mpistop('Change nproc or partitioning. Shutting down...', 0)
    endif
98  format(1x,'BAD PARTITIONING, the number of procs must be:',i6)

    ! Determination of the block number for each process
    ! ==================================================
    ! Global array nob is the block number for iproc [from 1 to nbloc]
    allocate(nob(0:nproc-1))
    ! determine block number
    cpt1=0
    do nbl=1,nbloc
       cpt2=cpt1+bl(nbl)%nproc-1
       nob(cpt1:cpt2)=nbl
       cpt1=cpt2+1
    enddo
    !print 99,iproc,nob(iproc)
99 format(1x,'proc',i4,' belongs to block',i4)

    ! Determination of the rank of the leader process per block
    ! =========================================================
    allocate(iproc_leader(nbloc))
    iproc_leader(1)=0
    do nbl=2,nbloc
       iproc_leader(nbl)=iproc_leader(nbl-1)+bl(nbl-1)%nproc
    enddo

    ! Current block
    ! =============
    nbl=nob(iproc)

    ! Compute local ni,nj,nk for each blocks
    ! ======================================
    if (mod(bl(nbl)%ni,bl(nbl)%ndomi).ne.0) then
       print 100,nbl
       call mpistop('WARNING: not possible (ni must be divisible by ndomi)', 0)
    endif
    if (mod(bl(nbl)%nj,bl(nbl)%ndomj).ne.0) then
       print 100,nbl
       call mpistop('WARNING: not possible (nj must be divisible by ndomj)', 0)
    endif
    if (mod(bl(nbl)%nk,bl(nbl)%ndomk).ne.0) then
       print 100,nbl
       call mpistop('WARNING: not possible (nk must be divisible by ndomk)', 0)
    endif
100  format(1x,'Problem for Block #',i3)

    ! Set global grid sizes from block definitions
    ! ============================================
    ngx=bl(nbl)%ni
    ngy=bl(nbl)%nj
    ngz=bl(nbl)%nk
    
    ! if nk=1 then call 2-D solver
    ! ----------------------------
    if (ngz==1) is_2d=.true.
    if (is_2d) then
       dim=2
    else
       dim=3
    endif

    ! Set number of proc per direction from block definitions
    ! =======================================================
    ndomx=bl(nbl)%ndomi
    ndomy=bl(nbl)%ndomj
    ndomz=bl(nbl)%ndomk

    ! Set local grid sizes after partitioning
    ! =======================================
    nx=ngx/ndomx
    ny=ngy/ndomy
    nz=ngz/ndomz

    ! Set local grid sizes extended to ghost points (inviscid)
    ! =============================================
    nx1= 1-ngh
    nx2=nx+ngh
    ny1= 1-ngh
    ny2=ny+ngh
    nz1= 1-ngh
    nz2=nz+ngh
    
    if (is_2d) then
       nz=1
       nz1=1
       nz2=1
    endif

    ! Set local grid sizes extended to ghost points (viscous)
    ! =============================================
    nx1_v= 1-ngh_v
    nx2_v=nx+ngh_v
    ny1_v= 1-ngh_v
    ny2_v=ny+ngh_v
    nz1_v= 1-ngh_v
    nz2_v=nz+ngh_v
    
    if (is_2d) then
       nz1_v=1
       nz2_v=1
    endif

    ! Print informations at screen
    ! ============================
    if (iproc==0) then
       do nbl=1,nbloc
          print 96,nbl,bl(nbl)%ni,bl(nbl)%nj,bl(nbl)%nk,bl(nbl)%ndomi,bl(nbl)%ndomj,bl(nbl)%ndomk
          print *,'==========='
          print 97,bl(nbl)%BC
       enddo
       print *,'==========='
       write(6,fmt='(1x,''size of MPI domains per block: ni='',i5,'', nj='',i5,'', nk='',i5)') nx,ny,nz
    endif
96  format(1x,'Block #',i3,' (',i5,' x',i5,' x',i5,' points) on (',i3,' x',i3,' x',i3,' procs)')
97  format(1x,'~> boundary conditions at Imin:',i3,' / Imax:',i3,' / Jmin:',i3,' / Jmax:',i3,' / Kmin:',i3,' / Kmax:',i3)

  end subroutine mpi_set_dim
    
  !===============================================================
  subroutine mpi_intrablock
  !===============================================================
    !> Definition of MPI partitioning
  !===============================================================
    use mod_block
    use mod_grid
    implicit none
    ! ------------------------------------------------------------
    integer :: nbl
    ! ------------------------------------------------------------
    ! communicator splitting
    integer :: color
    integer :: COMM_bl ! sub-comm after splitting and before MPI_CART    
    ! ------------------------------------------------------------
    ! MPI_CART intrablock subcomm
    integer :: dims(3)
    logical :: reorder
    ! ------------------------------------------------------------

    if (iproc==0) print *,'MPI Partitioning: intrablock communicator'
    
    ! Creation of sub-communicators for blocks
    ! ========================================
    nbl=nob(iproc)
    color=nbl-1
    call MPI_COMM_SPLIT(COMM_global,color,iproc,COMM_bl,info)
    
    ! Create MPI_CART intrablock subcommunicators
    ! ===========================================

    ! Fill dimensions
    dims(1)=bl(nbl)%ndomi
    dims(2)=bl(nbl)%ndomj
    dims(3)=bl(nbl)%ndomk  

    ! Determine periodicity for intrablock communications
    periods(1) = .false.
    periods(2) = .false.
    periods(3) = .false.
    if ((bl(nbl)%BC(1)>0).and.(bl(nbl)%BC(2)==nbl)) &
         periods(1) = .true.
    if ((bl(nbl)%BC(3)>0).and.(bl(nbl)%BC(4)==nbl)) &
         periods(2) = .true.
    if ((.not.is_2D).and. &
        (bl(nbl)%BC(5)>0).and.(bl(nbl)%BC(6)==nbl)) &
         periods(3) = .true.

    ! Reorder rank in intrablock
    reorder=.true.

    ! Create communicator
    call MPI_CART_CREATE(COMM_bl,3,dims,periods,reorder,COMM_intrablock,info)

    ! Coordinates of proc in MPI_CART grid
    call MPI_CART_GET(COMM_intrablock,3,dims,periods,coord,info)
    !print 101,iproc,coord(1),coord(2),coord(3)
101 format(1x,'iproc',i4,', coord in subcomm (',i4,',',i4,',',i4,')')

    call MPI_BARRIER(COMM_global,info)

  end subroutine mpi_intrablock

  !===============================================================
  subroutine mpi_connect
  !===============================================================
    !> Definition of MPI partitioning
  !===============================================================
    use mod_block
    use mod_constant
    use mod_flow
    use mod_bc_periodicity ! for: init periodicity
    use mod_grid_directions
    use mod_grid_directions_3d
    use warnstop
    implicit none
    ! ------------------------------------------------------------
    integer :: i,ip,nbl,j,k,n
    ! ------------------------------------------------------------
    ! communicator splitting
    integer :: color
    ! integer :: COMM_bl ! sub-comm after splitting and before MPI_CART
    ! ------------------------------------------------------------
    ! MPI_CART intrablock subcomm
    ! integer :: dims(3)
    ! logical :: reorder
    ! ------------------------------------------------------------
    ! gather MPI intrablock infos to search neighbors
    integer :: iproc_cart,n_shift,nbl_n
    integer, dimension(:), allocatable :: rank
    integer :: n_W,n_E,n_S,n_N,n_F,n_B
    integer, dimension(:), allocatable :: neighborW,neighborE,neighborS,neighborN
    integer :: i_r ! for interblock rank
    integer, dimension(:), allocatable :: i_rank
    integer :: ndi_n,ndj_n,ndk_n
    integer :: coord_1,coord_2,coord_3
    ! ------------------------------------------------------------
    ! search double corners
    integer :: nD1,nD2
    ! ------------------------------------------------------------

    if (iproc==0) print *,'MPI Partitioning: connectivity'
    if (iproc==0) print *,'----------------'

    ! Creation of sub-communicators for blocks
    ! ========================================
    nbl=nob(iproc)

!!$    color=nbl-1
!!$    call MPI_COMM_SPLIT(COMM_global,color,iproc,COMM_bl,info)
!!$
!!$    ! Create MPI_CART intrablock subcommunicators
!!$    ! ===========================================
!!$
!!$    ! Fill dimensions
!!$    dims(1)=bl(nbl)%ndomi
!!$    dims(2)=bl(nbl)%ndomj
!!$    dims(3)=bl(nbl)%ndomk
!!$
!!$    ! Determine periodicity for intrablock communications
!!$    periods(1) = .false.
!!$    periods(2) = .false.
!!$    periods(3) = .false.
!!$    if ((bl(nbl)%BC(1)>0).and.(bl(nbl)%BC(2)==nbl)) &
!!$         periods(1) = .true.
!!$    if ((bl(nbl)%BC(3)>0).and.(bl(nbl)%BC(4)==nbl)) &
!!$         periods(2) = .true.
!!$    if ((.not.is_2D).and. &
!!$        (bl(nbl)%BC(5)>0).and.(bl(nbl)%BC(6)==nbl)) &
!!$         periods(3) = .true.
!!$
!!$    ! Reorder rank in intrablock
!!$    reorder=.true.
!!$
!!$    ! Create communicator
!!$    call MPI_CART_CREATE(COMM_bl,3,dims,periods,reorder,COMM_intrablock,info)

    ! No proc in MPI_CART grid
    call MPI_COMM_RANK(COMM_intrablock,iproc_cart,info)

!!$    ! Coordinates of proc in MPI_CART grid
!!$    call MPI_CART_GET(COMM_intrablock,3,dims,periods,coord,info)
!!$    print 101,iproc,coord(1),coord(2),coord(3)
!!$101 format(1x,'iproc',i4,', coord in subcomm (',i4,',',i4,',',i4,')')

    call MPI_BARRIER(COMM_global,info)

    ! Gather MPI_CART partitioning on all procs
    ! =========================================
    allocate(rank(0:nproc-1))
    allocate(coord1(0:nproc-1),coord2(0:nproc-1),coord3(0:nproc-1))
    n_shift=0
    if (nbl>1) then
       do i=1,nbl-1
          n_shift=n_shift+bl(i)%nproc
       enddo
    endif
    ! gather global ranks in the MPI_CART communicator on all procs
    call MPI_ALLGATHER(iproc_cart+n_shift,1,MPI_INTEGER,rank,1,MPI_INTEGER,COMM_global,info)
    ! gather coord in the MPI_CART communicator
    call MPI_ALLGATHER(coord(1),1,MPI_INTEGER,coord1,1,MPI_INTEGER,COMM_global,info)
    call MPI_ALLGATHER(coord(2),1,MPI_INTEGER,coord2,1,MPI_INTEGER,COMM_global,info)
    call MPI_ALLGATHER(coord(3),1,MPI_INTEGER,coord3,1,MPI_INTEGER,COMM_global,info)
    call MPI_BARRIER(COMM_global,info)
   
    ! Search North, South, West, East, Front and Back neighbors
    ! =========================================================
    neighbor=MPI_PROC_NULL
    call MPI_CART_SHIFT(COMM_intrablock,0,1,n_W,n_E,info)
    if (n_W.ne.MPI_PROC_NULL) neighbor(nW)=n_W+n_shift
    if (n_E.ne.MPI_PROC_NULL) neighbor(nE)=n_E+n_shift
    call MPI_CART_SHIFT(COMM_intrablock,1,1,n_S,n_N,info)
    if (n_S.ne.MPI_PROC_NULL) neighbor(nS)=n_S+n_shift
    if (n_N.ne.MPI_PROC_NULL) neighbor(nN)=n_N+n_shift
    call MPI_CART_SHIFT(COMM_intrablock,2,1,n_F,n_B,info)
    if (n_F.ne.MPI_PROC_NULL) neighbor(nF)=n_F+n_shift
    if (n_B.ne.MPI_PROC_NULL) neighbor(nB)=n_B+n_shift
    !print 102,iproc,neighbor(1),neighbor(2),neighbor(3),neighbor(4),neighbor(5),neighbor(6)
102 format(1x,'proc',i4,': neighbor ',i4,i4,i4,i4,i4,i4)

    ! Init periodicity Boundary Conditions
    ! ====================================
    call init_periodicity

    ! Search swap and reverse directions
    ! ==================================
    if (is_curv3) then
       call grid_directions_3d(neighbor)
    else
       call grid_directions(neighbor(1:4))
    endif

    ! Define interblock neighbors
    ! ===========================
    
    ! Search neighbor of imin-face in neighboring block
    ! -------------------------------------------------
    if ((coord(1)==0).and.(bl(nbl)%BC(1)>0)) then
       nbl_n=bl(nbl)%BC(1)
       ndi_n=bl(nbl_n)%ndomi-1
       ndj_n=bl(nbl_n)%ndomj-1       
       ndk_n=bl(nbl_n)%ndomk-1

       do ip=bl(nbl_n)%proc_min,bl(nbl_n)%proc_max
          if (is_swapij2(1)) then
            if (is_rev2(1,2)) then
                coord_2=0
             else
                coord_2=ndj_n
             endif
             if (is_rev2(1,1)) then
                coord_1=ndi_n-coord(2)
             else
                coord_1=coord(2)
             endif
              if (is_rev2(1,3)) then
                coord_3=ndk_n-coord(3)
             else
                coord_3=coord(3)
             endif
          elseif (is_swapik2(1)) then
             if (is_rev2(1,2)) then
                coord_3=0
             else
                coord_3=ndk_n
             endif
             if (is_rev2(1,1)) then
                coord_2=ndj_n-coord(2)
             else
                coord_2=coord(2)
             endif
             if (is_rev2(1,3)) then
                coord_1=ndi_n-coord(3)
             else
                coord_1=coord(3)
             endif
          elseif (is_swapjk2(1)) then
             if (is_rev2(1,2)) then
                coord_1=0
             else
                coord_1=ndi_n
             endif
             if (is_rev2(1,1)) then
                coord_2=ndk_n-coord(3)
             else
                coord_2=coord(3)
             endif
             if (is_rev2(1,3)) then
                coord_3=ndj_n-coord(2)
             else
                coord_3=coord(2)
             endif
          elseif (is_swapijk2(1)) then !??? 2 possibilités
             if (is_rev2(1,1)) then
                coord_2=ndj_n-coord(2)
             else
                coord_2=coord(2)
             endif
             if (is_rev2(1,2)) then
                coord_1=0
             else
                coord_1=ndi_n
             endif
             if (is_rev2(1,3)) then
                coord_3=ndk_n-coord(3)
             else
                coord_3=coord(3)
             endif
          else
             if (is_rev2(1,2)) then
                coord_1=0
             else
                coord_1=ndi_n
             endif
             if (is_rev2(1,1)) then
                coord_2=ndj_n-coord(2)
             else
                coord_2=coord(2)
             endif
             if (is_rev2(1,3)) then
                coord_3=ndk_n-coord(3)
             else
                coord_3=coord(3)
             endif
          endif
          if ((coord1(ip)==coord_1).and. &
              (coord2(ip)==coord_2).and.(coord3(ip)==coord_3)) then
             neighbor(nW)=rank(ip)
             !print *,'proc',iproc,' (swapij',is_swapij2(1),'rev',is_rev2(1,:),') neighbor W after correct',neighbor(nW)
          endif
       enddo
    endif

    ! Search neighbor of imax-face in neighboring block
    ! -------------------------------------------------
    if ((coord(1)==ndomx-1).and.(bl(nbl)%BC(2)>0)) then
       nbl_n=bl(nbl)%BC(2)
       ndi_n=bl(nbl_n)%ndomi-1
       ndj_n=bl(nbl_n)%ndomj-1       
       ndk_n=bl(nbl_n)%ndomk-1

       do ip=bl(nbl_n)%proc_min,bl(nbl_n)%proc_max
          if (is_swapij2(2)) then
             if (is_rev2(2,1)) then
                coord_1=ndi_n-coord(2)
             else
                coord_1=coord(2)
             endif
             if (is_rev2(2,2)) then
                coord_2=ndj_n
             else
                coord_2=0
             endif
             if (is_rev2(2,3)) then
                coord_3=ndk_n-coord(3)
             else
                coord_3=coord(3)
             endif
          elseif (is_swapik2(2)) then
             if (is_rev2(2,2)) then
                coord_3=ndk_n
             else
                coord_3=0
             endif
             if (is_rev2(2,1)) then
                coord_2=ndj_n-coord(2)
             else
                coord_2=coord(2)
             endif
             if (is_rev2(2,3)) then
                coord_1=ndi_n-coord(3)
             else
                coord_1=coord(3)
             endif
          elseif (is_swapjk2(2)) then
             if (is_rev2(2,2)) then
                coord_1=ndi_n
             else
                coord_1=0
             endif
             if (is_rev2(2,1)) then
                coord_2=ndk_n-coord(3)
             else
                coord_2=coord(3)
             endif
             if (is_rev2(2,3)) then
                coord_3=ndj_n-coord(2)
             else
                coord_3=coord(2)
             endif
         elseif (is_swapijk2(2)) then ! ????? 2 possibilities
             if (is_rev2(2,1)) then
                coord_2=ndj_n-coord(2)
             else
                coord_2=coord(2)
             endif
             if (is_rev2(2,2)) then
                coord_1=ndi_n
             else
                coord_1=0
             endif
             if (is_rev2(2,3)) then
                coord_3=ndk_n-coord(3)
             else
                coord_3=coord(3)
             endif
          else
             if (is_rev2(2,1)) then
                coord_2=ndj_n-coord(2)
             else
                coord_2=coord(2)
             endif
             if (is_rev2(2,2)) then
                coord_1=ndi_n
            else
                coord_1=0
             endif
             if (is_rev2(2,3)) then
                coord_3=ndk_n-coord(3)
             else
                coord_3=coord(3)
             endif
          endif
          if ((coord1(ip)==coord_1).and. &
              (coord2(ip)==coord_2).and.(coord3(ip)==coord_3)) then
             neighbor(2)=rank(ip)
             !print *,'proc',iproc,' (swapij',is_swapij2(2),'rev',is_rev2(2,:),') neighbor E after correct',neighbor(nE)
          endif
       enddo
    endif
    
    ! Search neighbor of jmin-face in neighboring block
    ! -------------------------------------------------
    if ((coord(2)==0).and.(bl(nbl)%BC(3)>0)) then
       nbl_n=bl(nbl)%BC(3)
       ndi_n=bl(nbl_n)%ndomi-1
       ndj_n=bl(nbl_n)%ndomj-1       
       ndk_n=bl(nbl_n)%ndomk-1

       do ip=bl(nbl_n)%proc_min,bl(nbl_n)%proc_max
          if (is_swapij2(3)) then
             if (is_rev2(3,1)) then
                coord_2=ndj_n-coord(1)
             else
                coord_2=coord(1)
             endif
             if (is_rev2(3,2)) then
                coord_1=0
             else
                coord_1=ndi_n
             endif
             if (is_rev2(3,3)) then
                coord_3=ndk_n-coord(3)
             else
                coord_3=coord(3)
             endif
          elseif (is_swapik2(3)) then
             if (is_rev2(3,1)) then
                coord_1=ndk_n-coord(3)
             else
                coord_1=coord(3)
             endif
             if (is_rev2(3,2)) then
                coord_2=0
             else
                coord_2=ndj_n
             endif
             if (is_rev2(3,3)) then
                coord_3=ndi_n-coord(1)
             else
                coord_3=coord(1)
             endif
          elseif (is_swapjk2(3)) then
             if (is_rev2(3,1)) then
                coord_1=ndi_n-coord(1)
             else
                coord_1=coord(1)
             endif
             if (is_rev2(3,2)) then
                coord_3=0
             else
                coord_3=ndk_n
             endif
             if (is_rev2(3,3)) then
                coord_2=ndj_n-coord(3)
             else
                coord_2=coord(3)
             endif
          elseif (is_swapijk2(3)) then
             if (is_rev2(3,1)) then
                coord_3=ndk_n-coord(3)
             else
                coord_3=coord(3)
             endif
             if (is_rev2(3,2)) then
                coord_1=0
             else
                coord_1=ndi_n
             endif
             if (is_rev2(3,3)) then
                coord_2=ndk_n-coord(2)
             else
                coord_2=coord(2)
             endif
          else
             if (is_rev2(3,1)) then
                coord_1=ndi_n-coord(1)
             else
                coord_1=coord(1)
             endif
             if (is_rev2(3,2)) then
                coord_2=0
             else
                coord_2=ndj_n
             endif
             if (is_rev2(3,3)) then
                coord_3=ndk_n-coord(3)
             else
                coord_3=coord(3)
             endif
          endif
          if ((coord1(ip)==coord_1).and. &
              (coord2(ip)==coord_2).and.(coord3(ip)==coord_3)) then
             neighbor(3)=rank(ip)
             !print *,'proc',iproc,' (swapij',is_swapij2(3),'rev',is_rev2(3,:),') neighbor S after correct:',neighbor(nS)
          endif
       enddo
    endif

    ! Search neighbor of jmax-face in neighboring block
    ! -------------------------------------------------
    if ((coord(2)==ndomy-1).and.(bl(nbl)%BC(4)>0)) then
       nbl_n=bl(nbl)%BC(4)
       ndi_n=bl(nbl_n)%ndomi-1
       ndj_n=bl(nbl_n)%ndomj-1       
       ndk_n=bl(nbl_n)%ndomk-1

!!$       if (iproc==12) print *,'I am here, voisin',nbl_n
!!$       if (iproc==12) print *,'dim voisin',ndi_n,ndj_n,ndk_n
!!$       if (iproc==12) print *,'proc voisin',bl(nbl_n)%proc_min,bl(nbl_n)%proc_max
!!$       if (iproc==12) print *,'swapij',is_swapij2(4),'swapik',is_swapik2(4),'swapjk',is_swapjk2(4)
!!$       if (iproc==12) print *,'rev',is_rev2(4,:)

       do ip=bl(nbl_n)%proc_min,bl(nbl_n)%proc_max
          if (is_swapij2(4)) then
             if (is_rev2(4,1)) then
                coord_2=ndj_n-coord(1)
             else
                coord_2=coord(1)
             endif
             if (is_rev2(4,2)) then
                coord_1=ndi_n
             else
                coord_1=0
             endif
             if (is_rev2(4,3)) then
                coord_3=ndk_n-coord(3)
             else
                coord_3=coord(3)
             endif
          elseif (is_swapik2(4)) then
             if (is_rev2(4,1)) then
                coord_1=ndk_n-coord(3)
             else
                coord_1=coord(3)
             endif
             if (is_rev2(4,2)) then
                coord_2=ndj_n
             else
                coord_2=0
             endif
             if (is_rev2(4,3)) then
                coord_3=ndi_n-coord(1)
             else
                coord_3=coord(1)
             endif
          elseif (is_swapjk2(4)) then
             if (is_rev2(4,1)) then
                coord_1=ndi_n-coord(1)
             else
                coord_1=coord(1)
             endif
             if (is_rev2(4,2)) then
                coord_3=ndk_n
             else
                coord_3=0
             endif
             if (is_rev2(4,3)) then
                coord_2=ndj_n-coord(3)
             else
                coord_2=coord(3)
             endif
          elseif (is_swapijk2(4)) then
             if (is_rev2(4,1)) then
                coord_3=ndk_n-coord(3)
             else
                coord_3=coord(3)
             endif
             if (is_rev2(4,2)) then
                coord_1=ndi_n
             else
                coord_1=0
             endif
             if (is_rev2(4,3)) then
                coord_2=ndk_n-coord(2)
             else
                coord_2=coord(2)
             endif
          else
             if (is_rev2(4,1)) then
                coord_1=ndi_n-coord(1)
             else
                coord_1=coord(1)
             endif
             if (is_rev2(4,2)) then
                coord_2=ndj_n
             else
                coord_2=0
             endif
             if (is_rev2(4,3)) then
                coord_3=ndk_n-coord(3)
             else
                coord_3=coord(3)
             endif
          endif
          if ((coord1(ip)==coord_1).and. &
              (coord2(ip)==coord_2).and.(coord3(ip)==coord_3)) then
             neighbor(4)=rank(ip)
             !print *,'proc',iproc,'swapij',is_swapij2(4),'swapik',is_swapik2(4),'swapjk',is_swapjk2(4),'rev',is_rev2(4,:),') neighbor jmax after correct',neighbor(4)
          endif
       enddo
    endif

    ! Search neighbor of kmin-face in neighboring block
    ! -------------------------------------------------
    if ((coord(3)==0).and.(bl(nbl)%BC(5)>0)) then
       nbl_n=bl(nbl)%BC(5)
       ndi_n=bl(nbl_n)%ndomi-1
       ndj_n=bl(nbl_n)%ndomj-1
       ndk_n=bl(nbl_n)%ndomk-1

       do ip=bl(nbl_n)%proc_min,bl(nbl_n)%proc_max

          if (is_swapij2(5)) then ! TO BE DONE necessary ???
             if (is_rev2(5,1)) then
                coord_1=ndi_n-coord(2)
             else
                coord_1=coord(2)
             endif
             if (is_rev2(5,3)) then
                coord_2=ndj_n-coord(1)
             else
                coord_2=coord(1)
             endif
             if (is_rev2(5,2)) then
                coord_3=0
             else
                coord_3=ndk_n
             endif
          elseif (is_swapik2(5)) then
             if (is_rev2(5,1)) then
                coord_3=ndk_n-coord(1)
             else
                coord_3=coord(1)
             endif
             if (is_rev2(5,3)) then
                coord_2=ndj_n-coord(2)
             else
                coord_2=coord(2)
             endif
             if (is_rev2(5,2)) then
                coord_1=0
             else
                coord_1=ndi_n
             endif
          elseif (is_swapjk2(5)) then
             if (is_rev2(5,1)) then
                coord_1=ndi_n-coord(1)
             else
                coord_1=coord(1)
             endif
             if (is_rev2(5,3)) then
                coord_3=ndk_n-coord(2)
             else
                coord_3=coord(2)
             endif
             if (is_rev2(5,2)) then
                coord_2=0
             else
                coord_2=ndj_n
             endif
          elseif (is_swapijk2(5)) then ! TO BE DONE
             if (is_rev2(5,1)) then
                coord_1=ndi_n-coord(1)
             else
                coord_1=coord(1)
             endif
             if (is_rev2(5,3)) then
                coord_2=ndj_n-coord(2)
             else
                coord_2=coord(2)
             endif
             if (is_rev2(5,2)) then
                coord_3=0
             else
                coord_3=ndk_n
             endif
          else
             if (is_rev2(5,1)) then
                coord_1=ndi_n-coord(1)
             else
                coord_1=coord(1)
             endif
             if (is_rev2(5,3)) then
                coord_2=ndj_n-coord(2)
             else
                coord_2=coord(2)
             endif
             if (is_rev2(5,2)) then
                coord_3=0
             else
                coord_3=ndk_n
             endif
          endif
          if ((coord1(ip)==coord_1).and. &
              (coord2(ip)==coord_2).and.(coord3(ip)==coord_3)) then
             neighbor(5)=rank(ip)
          endif
       enddo
    endif
    
    ! Search neighbor of kmax-face in neighboring block
    ! -------------------------------------------------
    if ((coord(3)==ndomz-1).and.(bl(nbl)%BC(6)>0)) then
       nbl_n=bl(nbl)%BC(6)
       ndi_n=bl(nbl_n)%ndomi-1
       ndj_n=bl(nbl_n)%ndomj-1
       ndk_n=bl(nbl_n)%ndomk-1

       do ip=bl(nbl_n)%proc_min,bl(nbl_n)%proc_max

          if (is_swapij2(6)) then ! TO BE DONE necessary ???
             if (is_rev2(6,1)) then
                coord_1=ndi_n-coord(2)
             else
                coord_1=coord(2)
             endif
             if (is_rev2(6,3)) then
                coord_2=ndj_n-coord(1)
             else
                coord_2=coord(1)
             endif
             if (is_rev2(6,2)) then
                coord_3=ndk_n
             else
                coord_3=0
             endif
          elseif (is_swapik2(6)) then
             if (is_rev2(6,1)) then
                coord_3=ndk_n-coord(1)
             else
                coord_3=coord(1)
             endif
             if (is_rev2(6,3)) then
                coord_2=ndj_n-coord(2)
             else
                coord_2=coord(2)
             endif
             if (is_rev2(6,2)) then
                coord_1=ndi_n
             else
                coord_1=0
             endif
          elseif (is_swapjk2(6)) then
             if (is_rev2(6,1)) then
                coord_1=ndi_n-coord(1)
             else
                coord_1=coord(1)
             endif
             if (is_rev2(6,3)) then
                coord_3=ndk_n-coord(2)
             else
                coord_3=coord(2)
             endif
             if (is_rev2(6,2)) then
                coord_2=ndj_n
             else
                coord_2=0
             endif
          elseif (is_swapijk2(6)) then ! TO BE DONE
             if (is_rev2(6,1)) then
                coord_1=ndi_n-coord(1)
             else
                coord_1=coord(1)
             endif
             if (is_rev2(6,3)) then
                coord_2=ndj_n-coord(2)
             else
                coord_2=coord(2)
             endif
             if (is_rev2(6,2)) then
                coord_3=ndk_n
             else
                coord_3=0
             endif
          else
             if (is_rev2(6,1)) then
                coord_1=ndi_n-coord(1)
             else
                coord_1=coord(1)
             endif
             if (is_rev2(6,3)) then
                coord_2=ndj_n-coord(2)
             else
                coord_2=coord(2)
             endif
             if (is_rev2(6,2)) then
                coord_3=ndk_n
             else
                coord_3=0
             endif
          endif
          !coord_1=coord(1)
          !coord_2=coord(2)
          !coord_3=0

          if ((coord1(ip)==coord_1).and. &
              (coord2(ip)==coord_2).and.(coord3(ip)==coord_3)) then
             neighbor(6)=rank(ip)
          endif
       enddo
    endif

    !print 107,iproc,neighbor(1),neighbor(2),neighbor(3),neighbor(4),neighbor(5),neighbor(6)
107 format(1x,'proc',i4,': neighbor after correct ',i4,i4,i4,i4,i4,i4)

    !call mpistop('connectivity',0)

    ! Define tags for Send/Recv communications
    ! ========================================
    ! Send tags: tags
    ! ---------
    do i=1,6
       tags(i)=(iproc+1)*10+i
    enddo
    ! Receive tags: tagr
    ! ------------

    ! Pb periodicite LS59 9 blocs Euler
    if (TURB) then
       if ((nob(iproc)==1).or.(nob(iproc)==2).or.(nob(iproc)==8).or.(nob(iproc)==9)) then
       !if ((nob(iproc)==1).or.(nob(iproc)==2).or.(nob(iproc)==3).or.(nob(iproc)==4)) then
       !if ((iproc==0).or.(iproc==1).or.(iproc==7).or.(iproc==8)) then
          ndir(3)=4 ! by default 3/S exchange with 4/N
          ndir(4)=3 ! by default 4/N exchange with 3/S
       endif
    endif

!!$    ! Pb periodicite LS59 9 blocs Euler
!!$    if ((iproc==0).or.(iproc==3).or.(iproc==32).or.(iproc==35)) then
!!$       ndir(3)=4 ! by default 3/S exchange with 4/N
!!$       ndir(4)=3 ! by default 4/N exchange with 3/S
!!$    endif

!!$    ! Pb periodicite
!!$    ndir(5)=6 ! by default 3/S exchange with 4/N
!!$    ndir(6)=5 ! by default 4/N exchange with 3/S

!!$    if (iproc==4) then
!!$       ndir(3)=4
!!$       ndir(4)=3
!!$    else
!!$       ndir(5)=6
!!$       ndir(6)=5
!!$    endif

!!$       if (iproc==4) ndir(3)=6
!!$       if (iproc==5) ndir(6)=3

    do i=1,6
       ip=-1
       if (neighbor(i)>=0) ip=neighbor(i)
       tagr(i)=(ip+1)*10+ndir(i)
    enddo
    
    ! print 108,iproc,tags(1),tagr(1),tags(2),tagr(2),tags(3),tagr(3),tags(4),tagr(4)
    ! print 112,iproc,tags(5),tagr(5),tags(6),tagr(6)
108 format(1x,'iproc',i4,', tags/tagr:',i6,'/',i6,' ',i6,'/',i6,' ',i6,'/',i6,' ',i6,'/',i6)
112 format(1x,'iproc',i4,', tags/tagr k:',i6,'/',i6,' ',i6,'/',i6)

    ! call mpistop('here',0)

    ! Creation of interblock communicator
    ! ===================================
    ! (useful to communicate global grids)
    ! -----------------------------------   
    nbl=nob(iproc)
    if (iproc.eq.iproc_leader(nbl)) then
       color=0
    else
       color=1
    endif
    call MPI_COMM_SPLIT(COMM_global,color,iproc,COMM_interblock,info)

    ! Gather rank of all procs in subcommunicator COMM_interblock
    ! -----------------------------------------------------------
    call MPI_COMM_RANK(COMM_interblock,i_r,info)
    allocate(i_rank(0:nproc-1))
    call MPI_ALLGATHER(i_r,1,MPI_INTEGER,i_rank,1,MPI_INTEGER,COMM_global,info)

    call MPI_BARRIER(COMM_global,info)
    
    ! Neighborhood for interblock communications
    ! ------------------------------------------   
    neighbor_bl=MPI_PROC_NULL
    if (iproc.eq.iproc_leader(nbl)) then
       do i=1,6
          if (bl(nbl)%BC(i)>0) neighbor_bl(i)=i_rank(iproc_leader(bl(nbl)%BC(i)))
       enddo
       if (is_2D) neighbor_bl(5:6)=MPI_PROC_NULL

       !print 109,iproc,i_rank(iproc), &
       !     neighbor_bl(1),neighbor_bl(2),neighbor_bl(3),neighbor_bl(4),neighbor_bl(5),neighbor_bl(6)
    endif
109 format(1x,'proc',i4,' (',i4,')',': neighbor_bl',i4,i4,i4,i4,i4,i4)

    ! Detect blocks with double corners (at block level)
    ! =================================
    ! arises only for curvilinear grids with multiple points (e.g. 5-point corner)
    
    dble_corner=.false.
    dc_blocks=-1

    ! search only once per block, i.e. for leader procs:
    if ((is_curv.or.is_curv3) .and.(iproc.eq.iproc_leader(nob(iproc)))) then

       ! Share W,E,S,N block-neighbors (neighbor_bl) on all procs
       ! --------------------------------------------------------
       allocate(neighborW(0:nbloc-1),neighborE(0:nbloc-1))
       allocate(neighborS(0:nbloc-1),neighborN(0:nbloc-1))
       call MPI_ALLGATHER(neighbor_bl(nW),1,MPI_INTEGER,neighborW,1,MPI_INTEGER,COMM_interblock,info)
       call MPI_ALLGATHER(neighbor_bl(nE),1,MPI_INTEGER,neighborE,1,MPI_INTEGER,COMM_interblock,info)
       call MPI_ALLGATHER(neighbor_bl(nS),1,MPI_INTEGER,neighborS,1,MPI_INTEGER,COMM_interblock,info)
       call MPI_ALLGATHER(neighbor_bl(nN),1,MPI_INTEGER,neighborN,1,MPI_INTEGER,COMM_interblock,info)

       ! corner NW (1)
       ! -------------
       if ((neighbor_bl(nN).ne.MPI_PROC_NULL).and.(neighbor_bl(nW).ne.MPI_PROC_NULL)) then
          ! determine neighbor W of neighbor N
          if (is_swapij2_bl(nN)) then
             if (is_rev2_bl(nN,1)) then
                nD1=neighborN(neighbor_bl(nN))
             else
                nD1=neighborS(neighbor_bl(nN))
             endif
          else
             if (is_rev2_bl(nN,1)) then
                nD1=neighborE(neighbor_bl(nN))
             else
                nD1=neighborW(neighbor_bl(nN))
             endif
          endif
          ! determine neighbor N of neighbor W
          if (is_swapij2_bl(nW)) then
             if (is_rev2_bl(nW,1)) then
                nD2=neighborW(neighbor_bl(nW))
             else
                nD2=neighborE(neighbor_bl(nW))
             endif
          else
             if (is_rev2_bl(nW,1)) then
                nD2=neighborS(neighbor_bl(nW))
             else
                nD2=neighborN(neighbor_bl(nW))
             endif
          endif
          ! check if both determination leads different results
          if (nD1.ne.nD2) dble_corner(nNW)=.true.
          if ((nD1<0).and.(nD2<0)) dble_corner(nNW)=.true.
          !print *,'iproc',iproc,'W(N)',nD1,'N(W)',nD2
          
          ! store the rank of blocks sharing the double corner NW
          if (dble_corner(nNW)) dc_blocks(1,:) = (/i_rank(iproc),nD1,nD2,neighbor_bl(nN),neighbor_bl(nW)/)
       endif

       ! corner NE (2)
       ! -------------
       if ((neighbor_bl(nN).ne.MPI_PROC_NULL).and.(neighbor_bl(nE).ne.MPI_PROC_NULL)) then
          ! determine neighbor E of neighbor N
          if (is_swapij2_bl(nN)) then
             if (is_rev2_bl(nN,1)) then
                nD1=neighborS(neighbor_bl(nN))
             else
                nD1=neighborN(neighbor_bl(nN))
             endif
          else
             if (is_rev2_bl(nN,1)) then
                nD1=neighborW(neighbor_bl(nN))
             else
                nD1=neighborE(neighbor_bl(nN))
             endif
          endif
          ! determine neighbor N of neighbor E
          if (is_swapij2_bl(nE)) then
             if (is_rev2_bl(nE,1)) then
                nD2=neighborW(neighbor_bl(nE))
             else
                nD2=neighborE(neighbor_bl(nE))
             endif
          else
             if (is_rev2_bl(nE,1)) then
                nD2=neighborS(neighbor_bl(nE))
             else
                nD2=neighborN(neighbor_bl(nE))
             endif
          endif
          if (nD1.ne.nD2) dble_corner(nNE)=.true.
          if ((nD1<0).and.(nD2<0)) dble_corner(nNE)=.true.
          !print *,'iproc',iproc,'E(N)',nD1,'N(E)',nD2
          
          ! store the rank of blocks sharing the double corner NE
          if (dble_corner(nNE)) dc_blocks(2,:) = (/i_rank(iproc),nD1,nD2,neighbor_bl(nN),neighbor_bl(nE)/)
       endif

       ! corner SE (3)
       ! -------------
       if ((neighbor_bl(nS).ne.MPI_PROC_NULL).and.(neighbor_bl(nE).ne.MPI_PROC_NULL)) then
          ! determine neighbor E of neighbor S
          if (is_swapij2_bl(nS)) then
             if (is_rev2_bl(nS,1)) then
                nD1=neighborS(neighbor_bl(nS))
             else
                nD1=neighborN(neighbor_bl(nS))
             endif
          else
             if (is_rev2_bl(nS,1)) then
                nD1=neighborW(neighbor_bl(nS))
             else
                nD1=neighborE(neighbor_bl(nS))
             endif
          endif
          ! determine neighbor S of neighbor E
          if (is_swapij2_bl(nE)) then
             if (is_rev2_bl(nE,1)) then
                nD2=neighborE(neighbor_bl(nE))
             else
                nD2=neighborW(neighbor_bl(nE))
             endif
          else
             if (is_rev2_bl(nE,1)) then
                nD2=neighborN(neighbor_bl(nE))
             else
                nD2=neighborS(neighbor_bl(nE))
             endif
          endif
          if (nD1.ne.nD2) dble_corner(nSE)=.true.
          if ((nD1<0).and.(nD2<0)) dble_corner(nSE)=.true.
          !print *,'iproc',iproc,'E(S)',nD1,'S(E)',nD2
          
          ! store the rank of blocks sharing the double corner SE
          if (dble_corner(nSE)) dc_blocks(3,:) = (/i_rank(iproc),nD1,nD2,neighbor_bl(nS),neighbor_bl(nE)/)
       endif

       ! corner SW (4)
       ! -------------
       if ((neighbor_bl(nS).ne.MPI_PROC_NULL).and.(neighbor_bl(nW).ne.MPI_PROC_NULL)) then
          ! determine neighbor W of neighbor S
          if (is_swapij2_bl(nS)) then
             if (is_rev2_bl(nS,1)) then
                nD1=neighborN(neighbor_bl(nS))
             else
                nD1=neighborS(neighbor_bl(nS))
             endif
          else
             if (is_rev2_bl(nS,1)) then
                nD1=neighborE(neighbor_bl(nS))
             else
                nD1=neighborW(neighbor_bl(nS))
             endif
          endif
          ! determine neighbor S of neighbor W
          if (is_swapij2_bl(nW)) then
             if (is_rev2_bl(nW,1)) then
                nD2=neighborE(neighbor_bl(nW))
             else
                nD2=neighborW(neighbor_bl(nW))
             endif
          else
             if (is_rev2_bl(nW,1)) then
                nD2=neighborN(neighbor_bl(nW))
             else
                nD2=neighborS(neighbor_bl(nW))
             endif
          endif
          if (nD1.ne.nD2) dble_corner(nSW)=.true.
          if ((nD1<0).and.(nD2<0)) dble_corner(nSW)=.true.
          !print *,'iproc',iproc,'W(S)',nD1,'S(W)',nD2
          
          ! store the rank of blocks sharing the double corner SW
          if (dble_corner(nSW)) dc_blocks(4,:) = (/i_rank(iproc),nD1,nD2,neighbor_bl(nS),neighbor_bl(nW)/)
       endif

       ! free temporary shared neighbors
       if (.not.is_adjoint_block) deallocate(neighborW,neighborE,neighborS,neighborN)
       
       !------------------------------------------------------------------
       ! double corner arrangements for adjoint block case (starting here)
       !------------------------------------------------------------------      
       if (is_adjoint_block) then

          ! order (sort increasing) the ranks of blocks sharing the same double corner
          allocate(dc_mask(5))

          do j=1,4 ! each face
             dc_mask=.true.
             do i=1,5
                dc_tmp(i)=minval(dc_blocks(j,:),dc_mask)
                dc_mask(minloc(dc_blocks(j,:),dc_mask))=.false.
             enddo
             dc_blocks(j,:) = dc_tmp(:)
          enddo

          deallocate(dc_mask)

          ! dc_blocks rearrange in increasing order

          allocate(dc_mask(nbloc*4),dc_blocks_gb(nbloc*4,5),dc_blocks_gb_tmp(4,5,0:nbloc-1))

          ! gather all double corner block ranks
          call MPI_ALLGATHER(dc_blocks,20,MPI_INTEGER,dc_blocks_gb_tmp,20,MPI_INTEGER,COMM_interblock,info)

          ! rearrange global dc_blocks (first index is for all faces and blocks)
          k=0
          do j=0,nbloc-1
             do i=1,4
                k=k+1
                dc_blocks_gb(k,:)=dc_blocks_gb_tmp(i,:,j)
             enddo
          enddo

          deallocate(dc_blocks_gb_tmp)

          dc_mask=.true.

          ! eliminate the dublicate rows by logical masking
          do i=1,nbloc*4
             if(dc_blocks_gb(i,1)==-1) then
                dc_mask(i)=.false.
             else
                do j=i+1,nbloc*4
                   if (all(dc_blocks_gb(i,:)==dc_blocks_gb(j,:)).and.dc_mask(j)) dc_mask(j)=.false.
                enddo
             endif
          enddo

          no_dble_crnr = count(dc_mask)

       endif

    endif

    if (is_curv.and.is_adjoint_block.and.(.false.)) then

       ! broadcast the number of double corners and allocate "dble_crnr" derived type
       call MPI_BCAST(no_dble_crnr,1,MPI_INTEGER,iproc_leader(1),COMM_intrablock,info)

       ! allocate and initialize "dble_crnr" derived type
       allocate(dble_crnr(no_dble_crnr),dc_side(no_dble_crnr),dc_side_tmp(no_dble_crnr))
       dc_side = -1
       do i=1,no_dble_crnr
          dble_crnr(i)%dc_exist = .false.
          dble_crnr(i)%dc_index(:) = -1
          dble_crnr(i)%dc_blcs(:) = -1
          dble_crnr(i)%face_exist(:) = .false.
          do j=1,5
             dble_crnr(i)%face_index(j,:) = -1
          enddo
       enddo

       ! determine which side the double corners belong to (NW=1, NE=2, SE=3, or SW=4)
       if (iproc==iproc_leader(nob(iproc))) then
          i=0
          do j=1,nbloc*4
             if (dc_mask(j)) then
                i=i+1
                do k=1,4
                   if (all(dc_blocks(k,:)==dc_blocks_gb(j,:))) dc_side(i)=k
                enddo
             endif
          enddo
          deallocate(dc_blocks_gb,dc_mask)
       endif

       ! fill "dble_crnr" of each procs with the corner point info
       do i=1,nbloc
          
          dc_side_tmp = dc_side
          dc_blocks_tmp = dc_blocks

!!$          print *,'bloc',i,'dc_side',dc_side,'dc_blocks',dc_blocks
!!$          print *,'bloc',i,'alloc',allocated(dc_side),'size',size(dc_side)
          
          if (no_dble_crnr>0) call MPI_BCAST(dc_side_tmp,4,MPI_INTEGER,iproc_leader(i),COMM_global,info)
          call MPI_BCAST(dc_blocks_tmp,20,MPI_INTEGER,iproc_leader(i),COMM_global,info)
          
          if (i==nob(iproc)) then
             dc_side = dc_side_tmp
             dc_blocks = dc_blocks_tmp
             do j=1,no_dble_crnr
                if (dc_side(j)/=-1) then
                   if((dc_side(j)==1).and.(neighbor(nN).ne.MPI_PROC_NULL).and.(neighbor(nW).ne.MPI_PROC_NULL)) then ! NW
                      if((nob(neighbor(nN))/=nob(iproc)).and.(nob(neighbor(nW))/=nob(iproc))) then
                         dble_crnr(j)%dc_index = (/1 ,ny/)  ! (i,j)
                         dble_crnr(j)%dc_exist = .true.
                         dble_crnr(j)%dc_blcs(:) = dc_blocks(1,:)
                      endif
                   elseif((dc_side(j)==2).and.(neighbor(nN).ne.MPI_PROC_NULL).and.(neighbor(nE).ne.MPI_PROC_NULL)) then ! NE
                      if((nob(neighbor(nN))/=nob(iproc)).and.(nob(neighbor(nE))/=nob(iproc))) then
                         dble_crnr(j)%dc_index = (/nx,ny/) 
                         dble_crnr(j)%dc_exist = .true.
                         dble_crnr(j)%dc_blcs(:) = dc_blocks(2,:)
                      endif
                   elseif((dc_side(j)==3).and.(neighbor(nS).ne.MPI_PROC_NULL).and.(neighbor(nE).ne.MPI_PROC_NULL)) then ! SE
                      if((nob(neighbor(nS))/=nob(iproc)).and.(nob(neighbor(nE))/=nob(iproc))) then
                         dble_crnr(j)%dc_index = (/nx,1 /)
                         dble_crnr(j)%dc_exist = .true.
                         dble_crnr(j)%dc_blcs(:) = dc_blocks(3,:)
                      endif
                   elseif((dc_side(j)==4).and.(neighbor(nS).ne.MPI_PROC_NULL).and.(neighbor(nW).ne.MPI_PROC_NULL)) then ! SW
                      if((nob(neighbor(nS))/=nob(iproc)).and.(nob(neighbor(nW))/=nob(iproc))) then
                         dble_crnr(j)%dc_index = (/1 ,1 /)
                         dble_crnr(j)%dc_exist = .true.
                         dble_crnr(j)%dc_blcs(:) = dc_blocks(4,:)
                      endif
                   endif
                endif
             enddo
          endif
       enddo

!!$       if (iproc==3) then
!!$          print *,'======================='
!!$          do j=1,2
!!$             print *,iproc,'corner',j
!!$             print *,iproc,dble_crnr(j)%dc_index
!!$             print *,iproc,dble_crnr(j)%dc_exist
!!$             print *,iproc,dble_crnr(j)%dc_blcs
!!$          enddo
!!$          print *,'======================='
!!$       endif

       ! determine info of 5 interfaces connected to the same double corner
       if(iproc/=iproc_leader(nob(iproc))) then
          allocate(neighborW(0:nbloc-1),neighborE(0:nbloc-1))
          allocate(neighborS(0:nbloc-1),neighborN(0:nbloc-1))
       endif

       call MPI_BCAST(neighborW,nbloc,MPI_INTEGER,iproc_leader(1),COMM_global,info)
       call MPI_BCAST(neighborE,nbloc,MPI_INTEGER,iproc_leader(1),COMM_global,info)
       call MPI_BCAST(neighborS,nbloc,MPI_INTEGER,iproc_leader(1),COMM_global,info)
       call MPI_BCAST(neighborN,nbloc,MPI_INTEGER,iproc_leader(1),COMM_global,info)

       ! determine the interface ID and the neighbor block
       do i=1,no_dble_crnr
          if (dble_crnr(i)%dc_exist) then
             k = 1
             do j=1,5
                n = 1
                do while((j+n)<=5)
                   if (dble_crnr(i)%dc_blcs(j)>=0) then
                      if (neighborW(dble_crnr(i)%dc_blcs(j)).ne.MPI_PROC_NULL) then
                         if (neighborW(dble_crnr(i)%dc_blcs(j))==dble_crnr(i)%dc_blcs(j+n)) then
                            dble_crnr(i)%face_blcs(k,1) = dble_crnr(i)%dc_blcs(j)
                            dble_crnr(i)%face_blcs(k,2) = dble_crnr(i)%dc_blcs(j+n)
                            if(any((nob(iproc)-1)==dble_crnr(i)%face_blcs(k,:))) &
                                 dble_crnr(i)%face_exist(k) = .true.
                            k = k+1
                         endif
                      endif
                      if (neighborE(dble_crnr(i)%dc_blcs(j)).ne.MPI_PROC_NULL) then
                         if (neighborE(dble_crnr(i)%dc_blcs(j))==dble_crnr(i)%dc_blcs(j+n)) then
                            dble_crnr(i)%face_blcs(k,1) = dble_crnr(i)%dc_blcs(j)
                            dble_crnr(i)%face_blcs(k,2) = dble_crnr(i)%dc_blcs(j+n)
                            if(any((nob(iproc)-1)==dble_crnr(i)%face_blcs(k,:))) &
                                 dble_crnr(i)%face_exist(k) = .true.
                            k = k+1
                         endif
                      endif
                      if (neighborS(dble_crnr(i)%dc_blcs(j)).ne.MPI_PROC_NULL) then
                         if (neighborS(dble_crnr(i)%dc_blcs(j))==dble_crnr(i)%dc_blcs(j+n)) then
                            dble_crnr(i)%face_blcs(k,1) = dble_crnr(i)%dc_blcs(j)
                            dble_crnr(i)%face_blcs(k,2) = dble_crnr(i)%dc_blcs(j+n)
                            if(any((nob(iproc)-1)==dble_crnr(i)%face_blcs(k,:))) &
                                 dble_crnr(i)%face_exist(k) = .true.
                            k = k+1
                         endif
                      endif
                      if (neighborN(dble_crnr(i)%dc_blcs(j)).ne.MPI_PROC_NULL) then
                         if (neighborN(dble_crnr(i)%dc_blcs(j))==dble_crnr(i)%dc_blcs(j+n)) then
                            dble_crnr(i)%face_blcs(k,1) = dble_crnr(i)%dc_blcs(j)
                            dble_crnr(i)%face_blcs(k,2) = dble_crnr(i)%dc_blcs(j+n)
                            if(any((nob(iproc)-1)==dble_crnr(i)%face_blcs(k,:))) &
                                 dble_crnr(i)%face_exist(k) = .true.
                            k = k+1
                         endif
                      endif
                   endif
                   n = n+1
                enddo
             enddo

             do j=1,5
                if(dble_crnr(i)%face_exist(j)) then
                   if (neighbor(nW).ne.MPI_PROC_NULL) then
                      if (any((nob(neighbor(nW))-1)==dble_crnr(i)%face_blcs(j,:)).and.(nob(iproc)-1)/=(nob(neighbor(nW))-1)) then
                         dble_crnr(i)%face_index(j,1:3) = (/1,1,1/)
                         if(dble_crnr(i)%dc_index(2)==1)  dble_crnr(i)%face_index(j,4:6) = (/1,5,1/)
                         if(dble_crnr(i)%dc_index(2)==ny) dble_crnr(i)%face_index(j,4:6) = (/ny,ny-4,-1/)
                      endif
                   endif
                   if (neighbor(nE).ne.MPI_PROC_NULL) then
                      if (any((nob(neighbor(nE))-1)==dble_crnr(i)%face_blcs(j,:)).and.(nob(iproc)-1)/=(nob(neighbor(nE))-1)) then
                         dble_crnr(i)%face_index(j,1:3) = (/nx,nx,1/)
                         if(dble_crnr(i)%dc_index(2)==1)  dble_crnr(i)%face_index(j,4:6) = (/1,5,1/)
                         if(dble_crnr(i)%dc_index(2)==ny) dble_crnr(i)%face_index(j,4:6) = (/ny,ny-4,-1/)
                      endif
                   endif
                   if (neighbor(nS).ne.MPI_PROC_NULL) then
                      if (any((nob(neighbor(nS))-1)==dble_crnr(i)%face_blcs(j,:)).and.(nob(iproc)-1)/=(nob(neighbor(nS))-1)) then
                         dble_crnr(i)%face_index(j,4:6) = (/1,1,1/)
                         if(dble_crnr(i)%dc_index(1)==1)  dble_crnr(i)%face_index(j,1:3) = (/1,5,1/)
                         if(dble_crnr(i)%dc_index(1)==nx) dble_crnr(i)%face_index(j,1:3) = (/nx,nx-4,-1/)
                      endif
                   endif
                   if (neighbor(nN).ne.MPI_PROC_NULL) then
                      if (any((nob(neighbor(nN))-1)==dble_crnr(i)%face_blcs(j,:)).and.(nob(iproc)-1)/=(nob(neighbor(nN))-1)) then
                         dble_crnr(i)%face_index(j,4:6) = (/ny,ny,1/)
                         if(dble_crnr(i)%dc_index(1)==1)  dble_crnr(i)%face_index(j,1:3) = (/1,5,1/)
                         if(dble_crnr(i)%dc_index(1)==nx) dble_crnr(i)%face_index(j,1:3) = (/nx,nx-4,-1/)
                      endif
                   endif
                endif
             enddo
          endif
       enddo

       deallocate(neighborW,neighborE,neighborS,neighborN)

    endif
    !----------------------------------------------------------------
    ! double corner arrangements for adjoint block case (ending here)
    !----------------------------------------------------------------

    dble_corner=.true.    ! This line should be uncommented in case of 4-point corner in the vortex convection problem

    ! Define tags for Send/Recv interblock communications
    ! ===================================================
    if (iproc.eq.iproc_leader(nob(iproc))) then

       ! Send tags: tags_bl
       ! ------------------
       do i=1,6
          tags_bl(i)=(nob(iproc))*10+i
       enddo
       
       if (TURB) then
          ! Pb periodicite LS59 9 blocs Euler
          !if ((iproc==0).or.(iproc==1).or.(iproc==7).or.(iproc==8)) then
          !if ((nob(iproc)==1).or.(nob(iproc)==2).or.(nob(iproc)==3).or.(nob(iproc)==4)) then
          if ((nob(iproc)==1).or.(nob(iproc)==2).or.(nob(iproc)==8).or.(nob(iproc)==9)) then
             ndir_bl(3)=4 ! by default 3/S exchange with 4/N
             ndir_bl(4)=3 ! by default 4/N exchange with 3/S
          endif
       endif

!!$       ! Pb periodicite LS59 9 blocs Euler
!!$       if ((nob(iproc)==1).or.(nob(iproc)==2).or.(nob(iproc)==8).or.(nob(iproc)==9)) then
!!$          ndir_bl(3)=4 ! by default 3/S exchange with 4/N
!!$          ndir_bl(4)=3 ! by default 4/N exchange with 3/S
!!$       endif

!!$       ndir_bl(5)=6 ! by default 3/S exchange with 4/N
!!$       ndir_bl(6)=5 ! by default 4/N exchange with 3/S

!!$       if (iproc==4) ndir_bl(3)=6
!!$       if (iproc==5) ndir_bl(6)=3

       ! Receive tags: tagr_bl
       ! ---------------------
       do i=1,6
          ip=-1
          if (neighbor_bl(i)>=0) ip=neighbor_bl(i)
          tagr_bl(i)=(ip+1)*10+ndir_bl(i)
       enddo

       ! print 110,iproc,tags_bl(1),tagr_bl(1),tags_bl(2),tagr_bl(2),tags_bl(3),tagr_bl(3),tags_bl(4),tagr_bl(4)
       ! print 111,iproc,tags_bl(5),tagr_bl(5),tags_bl(6),tagr_bl(6)
110 format(1x,'iproc',i4,', tags_bl/tagr_bl:',i6,'/',i6,' ',i6,'/',i6,' ',i6,'/',i6,' ',i6,'/',i6)
111 format(1x,'iproc',i4,', tags_bl/tagr_bl k:',i6,'/',i6,' ',i6,'/',i6)
       
    endif

    ! call mpistop('check',0)

    ! Determine 'sizeofreal'
    ! ======================
    call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION,sizeofreal,info)
    
    ! Start and end indices in the global grid *** TO BE CHANGED for multiblock ***
    ! ========================================
    allocate(coordx(0:nproc-1),coordy(0:nproc-1),coordz(0:nproc-1))
    coordx(iproc)=coord(1)*nx+1
    coordy(iproc)=coord(2)*ny+1
    coordz(iproc)=coord(3)*nz+1
    ! Send information to other all procs
    do ip=0,nproc-1
       call MPI_BCAST(coordx(ip),1,MPI_INTEGER,ip,COMM_global,info)
    enddo
    do ip=0,nproc-1
       call MPI_BCAST(coordy(ip),1,MPI_INTEGER,ip,COMM_global,info)
    enddo
    do ip=0,nproc-1
       call MPI_BCAST(coordz(ip),1,MPI_INTEGER,ip,COMM_global,info)
    enddo
    
    ! Create sub-communicators for slicing
    ! ====================================
    call sub_comm

  end subroutine mpi_connect

  !===============================================================
  subroutine sub_comm
  !===============================================================
    !> Create sub-communicators for slicing
  !===============================================================
    use mod_constant
    use mod_grid
    implicit none
    ! ------------------------------------------------------------
    integer :: color
    ! ------------------------------------------------------------

    ! create communicator of all processors with the same x-coordinate: COMMX
    ! =======================================================================
    color = coord(1)
    call MPI_COMM_SPLIT(COMM_intrablock,color,iproc,COMMX,info)
    
    ! create communicator of all processors with the same y-coordinate: COMMY
    ! =======================================================================
    color = coord(2)
    call MPI_COMM_SPLIT(COMM_intrablock,color,iproc,COMMY,info)
    
    ! create communicator of all processors with the same z-coordinate: COMMZ
    ! =======================================================================
    if (is_2D) then
       COMMZ = COMM_intrablock
    else
       color = coord(3)
       call MPI_COMM_SPLIT(COMM_intrablock,color,iproc,COMMZ,info)
       call MPI_COMM_SIZE(COMMZ,nprocz,info)
    endif

    ! create communicator with same xy-coordinate: COMMXY
    ! ===================================================
    color = coord(2)
    call MPI_COMM_SPLIT(COMMX,color,iproc,COMMXY,info)
    call MPI_COMM_RANK(COMMXY,iprocxy,info)
    call MPI_COMM_SIZE(COMMXY,nprocxy,info)

    ! create communicator with same yz-coordinate: COMMYZ
    ! ===================================================
    color = coord(2)
    call MPI_COMM_SPLIT(COMMZ,color,iproc,COMMYZ,info)
    call MPI_COMM_RANK(COMMYZ,iprocyz,info)
    call MPI_COMM_SIZE(COMMYZ,nprocyz,info)

    ! create communicator with same xz-coordinate: COMMXZ
    ! ===================================================
    color = coord(1)
    call MPI_COMM_SPLIT(COMMZ,color,iproc,COMMXZ,info)
    call MPI_COMM_RANK(COMMXZ,iprocxz,info)
    call MPI_COMM_SIZE(COMMXZ,nprocxz,info)

!!$    ! Create sub-communicators along z-lines
!!$    ! ======================================
!!$    ! (used for averaging in spanwise dimension in post-processing)
!!$    allocate(ss_comm(0:ndomx-1,0:ndomy-1))
!!$    color=coord(1)+coord(2)*ndomx
!!$    call MPI_COMM_SPLIT(COMM_global,color,iproc,ss_comm(coord(1),coord(2)),info)   
!!$    !print 100,iproc,color
100 format(1x,'proc #',i4,' my color is',i3)

!!$    ! create communicator with the same double corner: dc_COMM (used for adjoint block cases)
!!$    ! =======================================================================
!!$    allocate(dc_COMM(no_dble_crnr,6))
!!$    do i=1,no_dble_crnr
!!$       color=0
!!$       if (dble_crnr(i)%dc_exist) color=i
!!$       !print *,'corner',i,'proc',iproc,'color',color
!!$       call MPI_COMM_SPLIT(COMM_global,color,iproc,dc_COMM(i,6),info) ! corner sub-communicator
!!$       color=0
!!$       do j=1,5
!!$          if (dble_crnr(i)%face_exist(j)) color=j
!!$          call MPI_COMM_SPLIT(COMM_global,color,iproc,dc_COMM(i,j),info) ! interface sub-communicator
!!$       enddo
!!$    enddo

  end subroutine sub_comm

end module mod_mpi_part
