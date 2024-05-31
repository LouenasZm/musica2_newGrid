!==============================================================================================
module mod_io_snapshots
!==============================================================================================
  !> authors: XG & AB (adapted from Luca)
  !> date: June 2023
  !> Module for all snapshot outputs ~> volumes, planes, lines [init/write/read]
!==============================================================================================
  use mod_io
  use mod_pp_vorticity
  implicit none
  ! -------------------------------------------------------------------------------------------
  ! SNAPSHOT TYPE defining its attributes
  ! -------------------------------------------------------------------------------------------
  type snapshot_type
     integer :: ind_i1, ind_i2                      ! snapshot position: i indexes
     integer :: ind_j1, ind_j2                      ! snapshot position: j indexes
     integer :: ind_k1, ind_k2                      ! snapshot position: k indexes
     integer :: nvar                                ! number of variables
     character(6), dimension(:), allocatable :: var ! name of variables
     logical :: has_snap                            ! boolean, true if proc has snapshot
     logical :: has_mz                              ! boolean, true if proc averages in z
     logical :: stamp                               ! indicator to write snapshot with timestamp
     logical :: is_check                            ! is_check: if T, suppressed at launch of simu.
     integer :: freq                                ! writting frequency
     integer :: type                                ! 0: point, 1: line, 2: plane, 3: vol
     integer :: type_ind                            ! position regarding the type (ex: plane num. 5)
     character(len=8) :: type_name='snapshot'       ! name of the type of data for filename
     type(teciotype)  :: tectype                    ! tecio type with I/O parameters
     ! for planes
     integer :: normal                              ! plane normal direction
     integer :: index                               ! plane position index
     integer, dimension(:), allocatable :: jvar     ! variable j-index (jva)
     ! for lines
     integer :: dir                                 ! line direction
  end type snapshot_type
  ! -------------------------------------------------------------------------------------------
  integer :: nsnapshots     ! number of snapshots
  integer :: npoints,nlines,nplanes,nvolumes     ! number of points, lines, planes & volumes
  integer, dimension(:), allocatable :: ipt_2_isn,ili_2_isn,ipl_2_isn,ivl_2_isn
  ! -------------------------------------------------------------------------------------------
  type(snapshot_type) :: volume_stats,plane_stats,line_stats,wall_stats
  type(snapshot_type), dimension(:), allocatable :: snapshots
  ! -------------------------------------------------------------------------------------------

contains

  !==============================================================================================
  subroutine init_io_snapshots
  !==============================================================================================
    !> authors: XG & AB
    !> date: June 2023
    !> Initialize all components needed for snapshot outputs
  !==============================================================================================
    use mod_utils
    use mod_block
    implicit none
    ! -------------------------------------------------------------------------------------------
    integer :: nbl,isn,m,cpt
    integer :: nvar_udf
    integer :: ipt,ili,ipl,ivl
    ! -------------------------------------------------------------------------------------------
    logical :: iexist
    character(:), allocatable :: snapshotfile
    ! -------------------------------------------------------------------------------------------

    ! Number of current block
    ! =======================
    nbl=nob(iproc)

    ! Initialization of indices position
    ! ==================================
    npoints=0; nlines=0; nplanes=0; nvolumes=0

    if (bl(nbl)%nsnapshot.gt.0) then
       ! Initialisation of subsnapshot
       ! ===========================
       if ((coord(1)==0).and.(coord(2)==0).and.(coord(3)==0)) print *,"init snapshots outputs"

       ! Allocate snapshots for a block
       ! ==============================
       nsnapshots=bl(nbl)%nsnapshot
       allocate(snapshots(nsnapshots))

       ! Determine maximum number of user-defined variables
       ! ==================================================
       nvar_udf=0
       do isn=1,nsnapshots
          cpt=0
          ! number of user-defined variables for the snapshot
          do m=1,bl(nbl)%snapshot(isn)%nvar
             if (trim(bl(nbl)%snapshot(isn)%var(m))=='udf') cpt=cpt+1
          enddo
          nvar_udf=max(nvar_udf,cpt)
       enddo
       ! allocate array of user-defined variables [uvar in mod_flow]
       if (nvar_udf>0) then
          allocate(uvar(nx,ny,nz,nvar_udf))
          uvar=0.0_wp
       endif

       !print *,'uvar allocated',nvar_udf

       ! Initialize snapshots
       ! ====================
       do isn=1,nsnapshots
          ! snapshot indexes
          snapshots(isn)%ind_i1=bl(nbl)%snapshot(isn)%ind_i1
          snapshots(isn)%ind_i2=bl(nbl)%snapshot(isn)%ind_i2
          snapshots(isn)%ind_j1=bl(nbl)%snapshot(isn)%ind_j1
          snapshots(isn)%ind_j2=bl(nbl)%snapshot(isn)%ind_j2
          snapshots(isn)%ind_k1=bl(nbl)%snapshot(isn)%ind_k1
          snapshots(isn)%ind_k2=bl(nbl)%snapshot(isn)%ind_k2
          ! snapshot variables (number & names)
          snapshots(isn)%nvar=bl(nbl)%snapshot(isn)%nvar
          allocate(snapshots(isn)%var(snapshots(isn)%nvar))
          snapshots(isn)%var=bl(nbl)%snapshot(isn)%var(1:bl(nbl)%snapshot(isn)%nvar)

          ! check if variables are averaged in z [has_mz attribute]
          if ((trim(snapshots(isn)%var(1))=='umz').or. &
              (trim(snapshots(isn)%var(1))=='Mmz').or. &
              (trim(snapshots(isn)%var(1))=='pmz').or. &
              (trim(snapshots(isn)%var(1))=='rhomz').or. &
              (trim(snapshots(isn)%var(1))=='droxmz').or. &
              (trim(snapshots(isn)%var(1))=='droymz').or. &
              (trim(snapshots(isn)%var(1))=='voxymz')) then
             snapshots(isn)%has_mz=.true.
             ! check
             if (snapshots(isn)%ind_k1.ne.snapshots(isn)%ind_k2) &
                  call mpistop('to create snapshot with z averaging, index k must be constant',0)
          else
             snapshots(isn)%has_mz=.false.
          endif

          ! snapshot timestamp
          snapshots(isn)%stamp=is_timestamp

          ! Initialization of teciotype for each snapshot
          ! ==========================================
          call init_snapshot(snapshots(isn))

          ! snapshot name
          write(snapshots(isn)%tectype%zonename,'(A5,I3.3)') 'snapshot',isn
          snapshots(isn)%tectype%strandID=isn

          ! snapshot frequency output
          ! -------------------------
          if (bl(nbl)%snapshot(isn)%freq.eq.0) then
             snapshots(isn)%is_check = .false.
             if (snapshots(isn)%type.eq.0) snapshots(isn)%freq = freq_point
             if (snapshots(isn)%type.eq.1) snapshots(isn)%freq = freq_line
             if (snapshots(isn)%type.eq.2) snapshots(isn)%freq = freq_plane
             if (snapshots(isn)%type.eq.3) snapshots(isn)%freq = freq_volume
          else if (bl(nbl)%snapshot(isn)%freq.lt.0) then
             snapshots(isn)%is_check = .true.
             snapshots(isn)%freq = - nmax/bl(nbl)%snapshot(isn)%freq
             if (snapshots(isn)%freq.eq.0) snapshots(isn)%freq=nmax
          else
             snapshots(isn)%is_check = .false.
             snapshots(isn)%freq = bl(nbl)%snapshot(isn)%freq
          endif

          ! Part added to allow appending in the same file
          ! ==============================================
          if (.not.snapshots(isn)%stamp) then

             ! logical is_app (is append_file)
             snapshots(isn)%tectype%is_app=.true.
             ! define offset: number of points * number of variables * 8 bytes
             snapshots(isn)%tectype%disp=snapshots(isn)%tectype%ngx &
                                    *snapshots(isn)%tectype%ngy &
                                    *snapshots(isn)%tectype%ngz
             snapshots(isn)%tectype%disp=snapshots(isn)%tectype%disp*8*snapshots(isn)%nvar
             ! inform about restart mode
             if ((idepart.eq.FROM_FILE).or.(idepart.eq.POST_PROCESSING)) then
                snapshots(isn)%tectype%restart=.true.
             else
                snapshots(isn)%tectype%restart=.false.
             endif
             ! if ((coord(1).eq.0).and.(coord(2).eq.0).and.(coord(3).eq.0)) print *,'disp to be checked',nbl,snapshots(isn)%tectype%disp
          endif
       enddo

       ! Save link between snapshot num. and vol,... + attribution of new numbers
       ! ------------------------------------------------------------------------
       allocate(ipt_2_isn(1:npoints))
       allocate(ili_2_isn(1:nlines))
       allocate(ipl_2_isn(1:nplanes))
       allocate(ivl_2_isn(1:nvolumes))
       ipt=0;ili=0;ipl=0;ivl=0
       do isn=1,nsnapshots
          if (snapshots(isn)%type.eq.0) then
             ipt = ipt + 1
             ipt_2_isn(ipt) = isn
             write(snapshots(isn)%tectype%zonename,'(A5,I3.3)') 'point',ipt
          else if (snapshots(isn)%type.eq.1) then
             ili = ili + 1
             ili_2_isn(ili) = isn
             write(snapshots(isn)%tectype%zonename,'(A5,I3.3)') 'line',ili
          else if (snapshots(isn)%type.eq.2) then
             ipl = ipl + 1
             ipl_2_isn(ipl) = isn
             write(snapshots(isn)%tectype%zonename,'(A5,I3.3)') 'plane',ipl
          else if (snapshots(isn)%type.eq.3) then
             ivl = ivl + 1
             ivl_2_isn(ivl) = isn
             write(snapshots(isn)%tectype%zonename,'(A5,I3.3)') 'volume',ivl
          endif
       enddo

       ! If is_check, remove file if exists and is_timestamp=F
       ! -----------------------------------------------------
       if (idepart.ne.4) then          ! Not done if post-processing
          if ((coord(1)==0).and.(coord(2)==0).and.(coord(3)==0)) then
             do isn=1,nsnapshots
                if (((snapshots(isn)%is_check).and.(.not.snapshots(isn)%stamp)).or. &
                    (idepart.eq.FROM_SCRATCH)) then
                   ! name of the snapshot
                   ! --------------------
                   if (.not.snapshots(isn)%tectype%is_IOtec_write) then
                      snapshotfile=trim(snapshots(isn)%type_name)//'_'//snapshots(isn)%tectype%zonename(6:8)//'_bl'//trim(numchar(nob(iproc)))//'.bin'
                   else
                      snapshotfile=trim(snapshots(isn)%type_name)//'_'//snapshots(isn)%tectype%zonename(6:8)//'_bl'//trim(numchar(nob(iproc)))//'.plt'
                   endif
                   inquire(file=trim(snapshotfile),exist=iexist)
                   if (iexist) call system('rm '//trim(snapshotfile))
                endif
             enddo
          endif
       endif

    endif

    !call mpistop('stop in init_io_snapshots',0)

    ! Initialization of planes or volumes for stats output TO BE CHANGED
    ! ====================================================
    if (STBL.or.CYL.or.SHIT.or.ACT.or.TURB.or.TE) then
       if (is_curv3) then
          if (TURB) then
             ! Turbine blade with rough wall [particular case]
             ! -> 3 outputs: planes xy, volumes & wall planes xz
             
             ! 1. Prepare a xy-plane
             ! ---------------------
             plane_stats%ind_i1=1; plane_stats%ind_i2=ngx
             plane_stats%ind_j1=1; plane_stats%ind_j2=ngy
             plane_stats%ind_k1=1; plane_stats%ind_k2=1
             plane_stats%normal=3
             plane_stats%index=1
             ! plane variables (number & names)
             plane_stats%nvar=1
             allocate(plane_stats%var(plane_stats%nvar))
             plane_stats%var='mean'
             ! plane timestamp
             plane_stats%stamp=is_timestamp

             write(plane_stats%tectype%zonename,'(A5)') 'stats'
             plane_stats%tectype%strandID = -2

             call init_snapshot(plane_stats)

             ! 2. Prepare a xyz-volume
             ! -----------------------
             ! volume indexes
             volume_stats%ind_i1=1; volume_stats%ind_i2=ngx
             volume_stats%ind_j1=1; volume_stats%ind_j2=ngy
             volume_stats%ind_k1=1; volume_stats%ind_k2=ngz
             ! volume variables (number & names)
             volume_stats%nvar=1
             allocate(volume_stats%var(volume_stats%nvar))
             volume_stats%var='mean'
             ! volume timestamp
             volume_stats%stamp=is_timestamp

             write(volume_stats%tectype%zonename,'(A5)') 'stats'
             volume_stats%tectype%strandID = -2

             call init_snapshot(volume_stats)
            
             ! 3.Prepare a xz-plane (only for jmin-walls)
             ! --------------------
             wall_stats%ind_i1=1; wall_stats%ind_i2=ngx
             wall_stats%ind_j1=1; wall_stats%ind_j2=1
             wall_stats%ind_k1=1; wall_stats%ind_k2=ngz
             wall_stats%normal=2
             wall_stats%index=1
             ! plane variables (number & names)
             wall_stats%nvar=1
             allocate(wall_stats%var(wall_stats%nvar))
             wall_stats%var='mean'
             ! plane timestamp
             wall_stats%stamp=is_timestamp

             write(wall_stats%tectype%zonename,'(A5)') 'stats'
             wall_stats%tectype%strandID = -2

             call init_snapshot(wall_stats)
          else            
             ! Prepare a xyz-volume
             ! --------------------
             ! volume indexes
             volume_stats%ind_i1=1; volume_stats%ind_i2=ngx
             volume_stats%ind_j1=1; volume_stats%ind_j2=ngy
             volume_stats%ind_k1=1; volume_stats%ind_k2=ngz
             ! volume variables (number & names)
             volume_stats%nvar=1
             allocate(volume_stats%var(volume_stats%nvar))
             volume_stats%var='mean'
             ! volume timestamp
             volume_stats%stamp=is_timestamp

             write(volume_stats%tectype%zonename,'(A5)') 'stats'
             volume_stats%tectype%strandID = -2

             call init_snapshot(volume_stats)
          endif
       else
          ! Prepare a xy-plane
          ! ------------------
          plane_stats%ind_i1=1; plane_stats%ind_i2=ngx
          plane_stats%ind_j1=1; plane_stats%ind_j2=ngy
          plane_stats%ind_k1=1; plane_stats%ind_k2=1
          plane_stats%normal=3
          plane_stats%index=1
          ! plane variables (number & names)
          plane_stats%nvar=1
          allocate(plane_stats%var(plane_stats%nvar))
          plane_stats%var='mean'
          ! plane timestamp
          plane_stats%stamp=is_timestamp

          write(plane_stats%tectype%zonename,'(A5)') 'stats'
          plane_stats%tectype%strandID = -2

          call init_snapshot(plane_stats)
       endif

    endif

  end subroutine init_io_snapshots

  !==============================================================================================
  subroutine init_snapshot(snapshot)
  !==============================================================================================
    !> authors: XG & AB
    !> date: June 2023
    !> Initialize the sub-communicator an I/O-type for snapshot outputs
  !==============================================================================================
    use mod_utils
    use mod_block
    use mod_constant
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Input/Output: plane type
    type(snapshot_type), intent(inout) :: snapshot
    ! -------------------------------------------------------------------------------------------
    logical :: find_offset
    integer :: nbl
    integer :: i,snapnprc
    integer :: GROUP_intrablock,GROUP_snapshot
    integer :: isnapshot,iproc_has_snapshot
    integer, dimension(3) :: coord_snap,coord_offset
    integer, dimension(:), allocatable :: iproc_has_snapshot_glo,snapshot_ranks
    ! -------------------------------------------------------------------------------------------

    ! Number of current block
    ! =======================
    nbl=nob(iproc)

    ! Define number of points and number of processors
    ! ================================================
    ! Verification if coherent limits given
    ! -------------------------------------
    if ((snapshot%ind_i1).gt.(snapshot%ind_i2)) call mpistop('  ~> Problem for sub-snapshot definition: i1 greater than i2',0)
    if ((snapshot%ind_j1).gt.(snapshot%ind_j2)) call mpistop('  ~> Problem for sub-snapshot definition: j1 greater than j2',0)
    if ((snapshot%ind_k1).gt.(snapshot%ind_k2)) call mpistop('  ~> Problem for sub-snapshot definition: k1 greater than k2',0)


    ! Check the snapshot type
    ! =======================
    ! Initialisation of specific attributes
    snapshot%dir=0; snapshot%normal=0; snapshot%index=0

    ! Check if point
    ! --------------
    if (((snapshot%ind_i1).eq.(snapshot%ind_i2)).and.((snapshot%ind_j1).eq.(snapshot%ind_j2))&
         .and.((snapshot%ind_k1).eq.(snapshot%ind_k2))) then
       snapshot%type = 0
    ! Check if line
    ! -------------
    else if (((snapshot%ind_i1).eq.(snapshot%ind_i2)).and.((snapshot%ind_j1).eq.(snapshot%ind_j2))) then
       snapshot%type = 1
       snapshot%dir  = 3
    else if (((snapshot%ind_j1).eq.(snapshot%ind_j2)).and.((snapshot%ind_k1).eq.(snapshot%ind_k2))) then
       snapshot%type = 1
       snapshot%dir  = 1
    else if (((snapshot%ind_i1).eq.(snapshot%ind_i2)).and.((snapshot%ind_k1).eq.(snapshot%ind_k2))) then
       snapshot%type = 1
       snapshot%dir  = 2
    ! Check if plane
    ! --------------
    else if ((snapshot%ind_i1).eq.(snapshot%ind_i2)) then
       snapshot%type = 2
       snapshot%normal = 1
       snapshot%index = snapshot%ind_i1
    else if ((snapshot%ind_j1).eq.(snapshot%ind_j2)) then
       snapshot%type = 2
       snapshot%normal = 2
       snapshot%index = snapshot%ind_j1
    else if ((snapshot%ind_k1).eq.(snapshot%ind_k2)) then
       snapshot%type = 2
       snapshot%normal = 3
       snapshot%index = snapshot%ind_k1
    else
    ! else, is volume
    ! ---------------
       snapshot%type = 3
    endif

    ! Specify name of data type and its position
    ! ------------------------------------------
    if (snapshot%type.eq.0) then
       snapshot%type_name = "point"
       npoints = npoints + 1
       snapshot%type_ind = npoints
    else if (snapshot%type.eq.1) then
       snapshot%type_name = "line"
       nlines = nlines + 1
       snapshot%type_ind = nlines
    else if (snapshot%type.eq.2) then
       snapshot%type_name = "plane"
       nplanes = nplanes + 1
       snapshot%type_ind = nplanes
    else if (snapshot%type.eq.3) then
       snapshot%type_name = "volume"
       nvolumes = nvolumes + 1
       snapshot%type_ind = nvolumes
    endif

    ! Check on dimensions
    ! -------------------
    ! Dimensions of snapshot taken to correspond to a MPI interface (interface between 2 processors)
    if ((float(snapshot%ind_i1-1)/nx.ne.float((snapshot%ind_i1-1)/nx)).and.(snapshot%ind_i1.lt.1)&
        .and.((snapshot%ind_i1).ne.(snapshot%ind_i2))) then
       snapshot%ind_i1 = max(((snapshot%ind_i1-1)/nx+1)*nx,1)
       if (iproc==0) print 101,snapshot%tectype%strandID,nbl,snapshot%ind_i1
    endif
101 format(1x,'~> Sub-snapshot ',i3,' of block ',i3,' : I1 not a imin index processor, value put to',i4)
    if (float(snapshot%ind_i2)/nx.ne.float((snapshot%ind_i2)/nx)&
        .and.((snapshot%ind_i1).ne.(snapshot%ind_i2))) then
       snapshot%ind_i2 = min(((snapshot%ind_i2)/nx+1)*nx,ngx)
       if (iproc==0) print 102,snapshot%tectype%strandID,nbl,snapshot%ind_i2
    endif
102 format(1x,'~> Sub-snapshot ',i3,' of block ',i3,' : I2 not a imax index processor, value put to',i4)
    if ((float(snapshot%ind_j1-1)/ny.ne.float((snapshot%ind_j1-1)/ny)).and.(snapshot%ind_j1.lt.1)&
        .and.((snapshot%ind_j1).ne.(snapshot%ind_j2))) then
       snapshot%ind_j1 = max(((snapshot%ind_j1-1)/ny+1)*ny,1)
       if (iproc==0) print 103,snapshot%tectype%strandID,nbl,snapshot%ind_j1
    endif
103 format(1x,'~> Sub-snapshot ',i3,' of block ',i3,' : J1 not a jmin index processor, value put to',i4)
    if (float(snapshot%ind_j2)/ny.ne.float((snapshot%ind_j2)/ny)&
        .and.((snapshot%ind_j1).ne.(snapshot%ind_j2))) then
       snapshot%ind_j2 = min(((snapshot%ind_j2)/ny+1)*ny,ngy)
       if (iproc==0) print 104,snapshot%tectype%strandID,nbl,snapshot%ind_j2
    endif
104 format(1x,'~> Sub-snapshot ',i3,' of block ',i3,' : J2 not a jmax index processor, value put to',i4)
    if ((float(snapshot%ind_k1-1)/nz.ne.float((snapshot%ind_k1-1)/nz)).and.(snapshot%ind_k1.lt.1)&
        .and.((snapshot%ind_k1).ne.(snapshot%ind_k2))) then
       snapshot%ind_k1 = max(((snapshot%ind_k1-1)/nz+1)*nz,1)
       if (iproc==0) print 105,snapshot%tectype%strandID,nbl,snapshot%ind_k1
    endif
105 format(1x,'~> Sub-snapshot ',i3,' of block ',i3,' : K1 not a kmin index processor, value put to',i4)
    if (float(snapshot%ind_k2)/nz.ne.float((snapshot%ind_k2)/nz)&
        .and.((snapshot%ind_k1).ne.(snapshot%ind_k2))) then
       snapshot%ind_k2 = min(((snapshot%ind_k2)/nz+1)*nz,ngz)
       if (iproc==0) print 106,snapshot%tectype%strandID,nbl,snapshot%ind_k2
    endif
106 format(1x,'~> Sub-snapshot ',i3,' of block ',i3,' : K2 not a kmax index processor, value put to',i4)

    ! Set the correct indexes
    ! -----------------------
    snapshot%tectype%ngx = snapshot%ind_i2 - snapshot%ind_i1 + 1
    snapshot%tectype%ngy = snapshot%ind_j2 - snapshot%ind_j1 + 1
    snapshot%tectype%ngz = snapshot%ind_k2 - snapshot%ind_k1 + 1

    ! Create snapshot communicators
    ! =============================

    ! Initialize proc list for the snapshots
    ! --------------------------------------
    iproc_has_snapshot=-1

    ! Determine procs belonging to a given snapshot
    ! ---------------------------------------------
    if (((coord(1)*nx.lt.snapshot%ind_i1.and.(coord(1)+1)*nx.ge.snapshot%ind_i1) .or. &
         (coord(1)*nx.lt.snapshot%ind_i2.and.(coord(1)+1)*nx.ge.snapshot%ind_i2) .or. &
         (coord(1)*nx.gt.snapshot%ind_i1.and. coord(1)   *nx.lt.snapshot%ind_i2)).and. &
        ((coord(2)*ny.lt.snapshot%ind_j1.and.(coord(2)+1)*ny.ge.snapshot%ind_j1) .or. &
         (coord(2)*ny.lt.snapshot%ind_j2.and.(coord(2)+1)*ny.ge.snapshot%ind_j2) .or.&
         (coord(2)*ny.gt.snapshot%ind_j1.and. coord(2)   *ny.lt.snapshot%ind_j2)).and.&
        ((coord(3)*nz.lt.snapshot%ind_k1.and.(coord(3)+1)*nz.ge.snapshot%ind_k1) .or. &
         (coord(3)*nz.lt.snapshot%ind_k1.and.(coord(3)+1)*nz.ge.snapshot%ind_k1) .or. &
         (coord(3)*nz.gt.snapshot%ind_k1.and. coord(3)   *nz.lt.snapshot%ind_k2))) then
         iproc_has_snapshot = 1
    endif

    ! Save of if this proc has this snapshot
    ! --------------------------------------
    if (iproc_has_snapshot.eq.-1) then
       snapshot%has_snap = .false.
    else
       snapshot%has_snap = .true.
    endif

    ! Communicate full list
    ! ---------------------
    allocate(iproc_has_snapshot_glo(bl(nbl)%nproc))
    iproc_has_snapshot_glo= -9999
    call MPI_ALLGATHER(iproc_has_snapshot,1,MPI_INTEGER,iproc_has_snapshot_glo, &
                       1,MPI_INTEGER,COMM_intrablock,info)

    ! Allocate snapshot list
    ! ----------------------
    ! snapnprc = 0
    ! do i=1,bl(nbl)%nproc
    !    if (iproc_has_snapshot_glo(i).gt.0) then
    !       snapnprc=snapnprc+1
    !    endif
    ! enddo
    snapnprc = 0; find_offset = .true.; coord_offset = 0; coord_snap = 0
    do i=1,bl(nbl)%nproc
       if (iproc_has_snapshot_glo(i).gt.0) then
          snapnprc=snapnprc+1
          if (find_offset) then
             find_offset = .false.
             coord_offset(1) = (i-1)/(ndomz*ndomy)
             coord_offset(2) = ((i-1)-coord_offset(1)*ndomz*ndomy)/ndomz
             coord_offset(3) = ((i-1)-coord_offset(1)*ndomz*ndomy-coord_offset(2)*ndomz)
          endif
       endif
    enddo

    do i=1,3
       coord_snap(i) = coord(i) - coord_offset(i)
    enddo

    ! Fill snapshot list
    ! ------------------
    allocate(snapshot_ranks(snapnprc))
    isnapshot= 1
    do i=1,bl(nbl)%nproc
       if (iproc_has_snapshot_glo(i).gt.0) then
          snapshot_ranks(isnapshot)= i-1
          isnapshot=isnapshot+1
       endif
    enddo

    ! Create sub-communicator
    ! -----------------------
    ! Get world group
    call MPI_COMM_GROUP(COMM_intrablock,GROUP_intrablock,info)
    ! Create plane group
    call MPI_GROUP_INCL(GROUP_intrablock,snapnprc,snapshot_ranks,GROUP_snapshot,info)
    ! Destroy world group
    call MPI_GROUP_FREE(GROUP_intrablock,info)
    ! Create communicator
    call MPI_COMM_CREATE(COMM_intrablock,GROUP_snapshot,snapshot%tectype%MPI_COMM,info)
    ! Destroy plane group
    call MPI_GROUP_FREE(GROUP_snapshot,info)

    ! free allocated memory
    deallocate(snapshot_ranks,iproc_has_snapshot_glo)

    ! Initialize plane local indices
    ! ==============================

    ! default initializations
    ! -----------------------
    snapshot%tectype%nx1 = 1;  snapshot%tectype%nx2 = nx
    snapshot%tectype%ny1 = 1;  snapshot%tectype%ny2 = ny
    snapshot%tectype%nz1 = 1;  snapshot%tectype%nz2 = nz

    ! If point, 3 constant indices
    if (snapshot%type.eq.0) then
       snapshot%tectype%nx1 = max(1,min(nx,snapshot%ind_i1-coord(1)*nx))
       snapshot%tectype%nx2 = snapshot%tectype%nx1
       snapshot%tectype%ny1 = max(1,min(ny,snapshot%ind_j1-coord(2)*ny))
       snapshot%tectype%ny2 = snapshot%tectype%ny1
       snapshot%tectype%nz1 = max(1,min(nz,snapshot%ind_k1-coord(3)*nz))
       snapshot%tectype%nz2 = snapshot%tectype%nz1
    ! If lines, indices of 2 constants directions put to index value
    else if (snapshot%type.eq.1) then
       if (snapshot%dir.eq.1) then
          snapshot%tectype%ny1 = max(1,min(ny,snapshot%ind_j1-coord(2)*ny))
          snapshot%tectype%ny2 = snapshot%tectype%ny1
          snapshot%tectype%nz1 = max(1,min(nz,snapshot%ind_k1-coord(3)*nz))
          snapshot%tectype%nz2 = snapshot%tectype%nz1
       else if (snapshot%dir.eq.2) then
          snapshot%tectype%nx1 = max(1,min(nx,snapshot%ind_i1-coord(1)*nx))
          snapshot%tectype%nx2 = snapshot%tectype%nx1
          snapshot%tectype%nz1 = max(1,min(nz,snapshot%ind_k1-coord(3)*nz))
          snapshot%tectype%nz2 = snapshot%tectype%nz1
       else if (snapshot%dir.eq.3) then
          snapshot%tectype%nx1 = max(1,min(nx,snapshot%ind_i1-coord(1)*nx))
          snapshot%tectype%nx2 = snapshot%tectype%nx1
          snapshot%tectype%ny1 = max(1,min(ny,snapshot%ind_j1-coord(2)*ny))
          snapshot%tectype%ny2 = snapshot%tectype%ny1
       endif
    ! If planes, indices of normal direction put to index value
    else if (snapshot%type.eq.2) then
       if (snapshot%normal.eq.1) then
          snapshot%tectype%nx1 = max(1,min(nx,snapshot%index-coord(1)*nx))
          snapshot%tectype%nx2 = snapshot%tectype%nx1
       else if (snapshot%normal.eq.2) then
          snapshot%tectype%ny1 = max(1,min(ny,snapshot%index-coord(2)*ny))
          snapshot%tectype%ny2 = snapshot%tectype%ny1
       else if (snapshot%normal.eq.3) then
          snapshot%tectype%nz1 = max(1,min(nz,snapshot%index-coord(3)*nz))
          snapshot%tectype%nz2 = snapshot%tectype%nz1
       endif
    endif

    ! Initialize I/O type
    ! ===================
    if (idepart==POST_PROCESSING) then
       if (iproc_has_snapshot.eq.1) then
          call mod_io_init(snapshot%tectype%ngx,snapshot%tectype%ngy,snapshot%tectype%ngz, &
               snapshot%tectype%nx2-snapshot%tectype%nx1+1, &
               snapshot%tectype%ny2-snapshot%tectype%ny1+1, &
               snapshot%tectype%nz2-snapshot%tectype%nz1+1, &
               snapshot%tectype%nx2-snapshot%tectype%nx1+1, &
               snapshot%tectype%ny2-snapshot%tectype%ny1+1, &
               snapshot%tectype%nz2-snapshot%tectype%nz1+1, &
               ngh,3,coord_snap,is_IOtec_read,is_IOtec_write,snapshot%tectype)
       endif
    else
       if (iproc_has_snapshot.eq.1) then
          call mod_io_init(snapshot%tectype%ngx,snapshot%tectype%ngy,snapshot%tectype%ngz, &
               snapshot%tectype%nx2-snapshot%tectype%nx1+1, &
               snapshot%tectype%ny2-snapshot%tectype%ny1+1, &
               snapshot%tectype%nz2-snapshot%tectype%nz1+1, &
               snapshot%tectype%nx2-snapshot%tectype%nx1+1, &
               snapshot%tectype%ny2-snapshot%tectype%ny1+1, &
               snapshot%tectype%nz2-snapshot%tectype%nz1+1, &
               ngh,3,coord_snap,.false.,.false.,snapshot%tectype)
       endif
    endif

  end subroutine init_snapshot

  !==============================================================================================
  subroutine write_grid_snapshot(isn)
  !==============================================================================================
    !> authors: XG & AB
    !> date: June 2023
    !> Write GRID for snapshot (only for tecplot format)
  !==============================================================================================
    use mod_utils
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Input: index of snapshot to be written
    integer, intent(in) :: isn
    ! -------------------------------------------------------------------------------------------
    integer :: i,j,k,m,filetype,ndata
    integer :: nx1l,nx2l,ny1l,ny2l,nz1l,nz2l
    !character(len=23) :: gridfile='./snapshot_xxx_grid.plt'
    character(:), allocatable :: gridfile
    logical :: iexist
    ! -------------------------------------------------------------------------------------------

    ! Prepare number and name of data: x,y,z coordinates
    ! ===============================
    ndata=3
    allocate(dataname(ndata),varlist(ndata))
    dataname(1:ndata) = (/'X','Y','Z'/)

    ! for snapshots isn
    ! =================

    ! Check if the snapshot belong to this proc
    ! -----------------------------------------
    if (.not.(snapshots(isn)%has_snap)) return

    ! name of the snapshot: gridfile "snapshot_xxx_grid.plt"
    ! -------------------
    !gridfile(12:14)=snapshots(isn)%tectype%zonename(6:8)
    gridfile='grid'//trim(snapshots(isn)%type_name)// &
             '_'//snapshots(isn)%tectype%zonename(6:8)// &
             '_bl'//trim(numchar(nob(iproc)))//'.plt'

    ! check if grid file already exists
    ! ---------------------------------
    inquire(file=gridfile,exist=iexist)

    ! if not: write GRID
    ! ------------------
    if (.not.iexist) then

       ! snapshot indices
       nx1l= snapshots(isn)%tectype%nx1
       nx2l= snapshots(isn)%tectype%nx2
       ny1l= snapshots(isn)%tectype%ny1
       ny2l= snapshots(isn)%tectype%ny2
       nz1l= snapshots(isn)%tectype%nz1
       nz2l= snapshots(isn)%tectype%nz2

       ! allocate and initialize dummy variable for datas
       allocate(dummy(nx1l:nx2l,ny1l:ny2l,nz1l:nz2l,ndata))

       ! fill data x,y,z
       if (is_curv3) then ! curvilinear grid
          do k=nz1l,nz2l
             do j=ny1l,ny2l
                do i=nx1l,nx2l
                   dummy(i,j,k,1)= xc3(i,j,k)
                   dummy(i,j,k,2)= yc3(i,j,k)
                   dummy(i,j,k,3)= zc3(i,j,k)
                enddo
             enddo
          enddo
       elseif (is_curv) then ! curvilinear grid
          do k=nz1l,nz2l
             do j=ny1l,ny2l
                do i=nx1l,nx2l
                   dummy(i,j,k,1)= xc(i,j)
                   dummy(i,j,k,2)= yc(i,j)
                   dummy(i,j,k,3)= z(k)
                enddo
             enddo
          enddo
       else ! Cartesian grid
          do k=nz1l,nz2l
             do j=ny1l,ny2l
                do i=nx1l,nx2l
                   dummy(i,j,k,1)= x(i)
                   dummy(i,j,k,2)= y(j)
                   dummy(i,j,k,3)= z(k)
                enddo
             enddo
          enddo
       endif

       ! fill varlist (pointer on datas)
       do m=1,ndata
          varlist(m)%data=>dummy(:,:,:,m)
          varlist(m)%name= dataname(m)
       enddo

       ! write data operation
       if (iproc.eq.0) write(*,*) 'Writing grid snapshot... '//gridfile
       filetype=1 ! [0:full; 1:only grid; 2:only solution]
       call write_tec(gridfile,filetype,tstar,snapshots(isn)%tectype)

       ! free temporary data
       deallocate(dummy)

    endif ! end if test iexist

    ! free pointer varlist & dataname
    ! ===============================
    deallocate(varlist,dataname)

  end subroutine write_grid_snapshot

  !==============================================================================================
  subroutine write_snapshot(isn)
  !==============================================================================================
    !> authors: XG & AB
    !> date: June 2023
    !> Write variables for snapshot # isn
  !==============================================================================================
    use mod_eos
    use mod_utils
    use mod_deriv
    use mod_deriv_c
    use mod_deriv_c3
    use mod_rans
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Input: index of snapshot to be written
    integer, intent(in) :: isn
    ! -------------------------------------------------------------------------------------------
    ! Local variables
    integer :: i,j,k,m,filetype,ndata,cpt
    integer :: nx1sn,nx2sn,ny1sn,ny2sn,nz1sn,nz2sn
    ! -------------------------------------------------------------------------------------------
    character(:), allocatable :: snapshotfile
    ! -------------------------------------------------------------------------------------------
    ! RANS variables
    real(wp) :: nu,cnu13
    real(wp), dimension(nx,ny) :: mz
    ! real(wp) :: nu,khi,fnu1
    ! integer  :: cnu1
    ! -------------------------------------------------------------------------------------------

    ! Check if the snapshot belong to this proc
    ! -----------------------------------------
    if ((.not.(snapshots(isn)%has_snap)).and.(.not.(snapshots(isn)%has_mz))) return

    if (snapshots(isn)%has_snap) then

       ! Determine snapshot name
       ! =======================
       if (snapshots(isn)%stamp) then
          ! name of the snapshot with timestamp
          ! --------------------------------
          if (.not.snapshots(isn)%tectype%is_IOtec_write) then
             snapshotfile=trim(snapshots(isn)%type_name)// &
                               '_'//snapshots(isn)%tectype%zonename(6:8)// &
                               '_bl'//trim(numchar(nob(iproc)))// &
                               '_'//filestamp//'.bin'
          else
             snapshotfile=trim(snapshots(isn)%type_name)// &
                               '_'//snapshots(isn)%tectype%zonename(6:8)// &
                               '_bl'//trim(numchar(nob(iproc)))// &
                               '_'//filestamp//'.plt'
          endif
       else
          ! name of the snapshot without timestamp
          ! -----------------------------------
          if (.not.snapshots(isn)%tectype%is_IOtec_write) then
             snapshotfile=trim(snapshots(isn)%type_name)// &
                               '_'//snapshots(isn)%tectype%zonename(6:8)// &
                               '_bl'//trim(numchar(nob(iproc)))//'.bin'
          else
             snapshotfile=trim(snapshots(isn)%type_name)// &
                               '_'//snapshots(isn)%tectype%zonename(6:8)// &
                               '_bl'//trim(numchar(nob(iproc)))//'.plt'
          endif
       endif

       ! Prepare number and name of data
       ! ===============================
       ndata=snapshots(isn)%nvar
       allocate(varlist(ndata),dataname(ndata))
       do m=1,ndata
          dataname(m)=trim(snapshots(isn)%var(m))
       enddo

       ! Snapshot indices
       ! =============
       nx1sn= snapshots(isn)%tectype%nx1
       nx2sn= snapshots(isn)%tectype%nx2
       ny1sn= snapshots(isn)%tectype%ny1
       ny2sn= snapshots(isn)%tectype%ny2
       nz1sn= snapshots(isn)%tectype%nz1
       nz2sn= snapshots(isn)%tectype%nz2

       ! allocate dummy variable for datas
       ! ---------------------------------
       allocate(dummy(nx1sn:nx2sn,ny1sn:ny2sn,nz1sn:nz2sn,ndata))

       ! Fill data in temporary file dummy
       ! =================================
       ! List of possible variables: prs,uu,vv,ww,rho,Tmp,div,Mach,Gamma,Frhov,Grhow,udf

       ! Initialize counter for user-defined (udf) variables
       cpt=1

       ! Spalart-Allmaras constants
       cnu1 = 7.1_wp
       cnu13=cnu1**3

       do m=1,ndata
          select case (trim(snapshots(isn)%var(m)))
          case('prs') ! pressure
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      dummy(i,j,k,m)=prs(i,j,k)
                   enddo
                enddo
             enddo
          case('mut') ! Mu_t
             if (is_RANS) then
                do k=nz1sn,nz2sn
                   do j=ny1sn,ny2sn
                      do i=nx1sn,nx2sn
                         nu   = visc(i,j,k)/rho(i,j,k)
                         khi  = nutil(i,j,k)/nu
                         fnu1 = khi**3/(khi**3+cnu13)
                         dummy(i,j,k,m)=rho(i,j,k)*nutil(i,j,k)*fnu1
                         !dummy(i,j,k,m)=Sterm(i,j,k)
                      enddo
                   enddo
                enddo
             else
                dummy(:,:,:,m)=0.0_wp
             endif
          case('lDES') ! lengthscale RANS
             if (is_RANS) then
                do k=nz1sn,nz2sn
                   do j=ny1sn,ny2sn
                      do i=nx1sn,nx2sn
                         dummy(i,j,k,m)=lengthscale(i,j,k)
                      enddo
                   enddo
                enddo
             else
                dummy(:,:,:,m)=0.0_wp
             endif
          case( 'uu') ! streamwise velocity
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      dummy(i,j,k,m)= uu(i,j,k)
                   enddo
                enddo
             enddo
          case( 'vv') ! crossflow velocity
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      dummy(i,j,k,m)= vv(i,j,k)
                   enddo
                enddo
             enddo
          case( 'ww') ! spanwise velocity
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      dummy(i,j,k,m)= ww(i,j,k)
                   enddo
                enddo
             enddo
          case('rho') ! density
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      dummy(i,j,k,m)=rho(i,j,k)
                   enddo
                enddo
             enddo
          case('Tmp') ! temperature
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      dummy(i,j,k,m)=Tmp(i,j,k)
                   enddo
                enddo
             enddo
          case('div') ! divergence
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      dummy(i,j,k,m)=dux(i,j,k)+dvy(i,j,k)+dwz(i,j,k)
                   enddo
                enddo
             enddo
          case('duy') ! divergence
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      dummy(i,j,k,m)=duy(i,j,k)
                   enddo
                enddo
             enddo
          case('duz') ! divergence
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      dummy(i,j,k,m)=duz(i,j,k)
                   enddo
                enddo
             enddo
          case('dwx') ! spanwise velocity gradient in x direction
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      dummy(i,j,k,m)=dwx(i,j,k)
                   enddo
                enddo
             enddo
          case('dwy') ! spanwise velocity gradient in y direction
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      dummy(i,j,k,m)=dwy(i,j,k)
                   enddo
                enddo
             enddo
          case('drhoi') ! density gradient in x direction
             if (is_curv3) then
                call deriv_x_5pts_c3(rho(nx1:nx2,ny1:ny2,nz1:nz2), &
                     dummy(nx1sn:nx2sn,ny1sn:ny2sn,nz1sn:nz2sn,m), &
                     nx1sn,nx2sn,ny1sn,ny2sn,nz1sn,nz2sn)
             elseif (is_curv) then
                call deriv_x_5pts_c(rho(nx1:nx2,ny1:ny2,nz1:nz2), &
                     dummy(nx1sn:nx2sn,ny1sn:ny2sn,nz1sn:nz2sn,m), &
                     nx1sn,nx2sn,ny1sn,ny2sn,nz1sn,nz2sn)
             else
                call mpistop('Output drhoi not implemented yet !',0)
             endif
          case('drhoj') ! density gradient in y direction
             if (is_curv3) then
                call deriv_y_5pts_c3(rho(nx1:nx2,ny1:ny2,nz1:nz2), &
                     dummy(nx1sn:nx2sn,ny1sn:ny2sn,nz1sn:nz2sn,m), &
                     nx1sn,nx2sn,ny1sn,ny2sn,nz1sn,nz2sn)
             elseif (is_curv) then
                call deriv_y_5pts_c(rho(nx1:nx2,ny1:ny2,nz1:nz2), &
                     dummy(nx1sn:nx2sn,ny1sn:ny2sn,nz1sn:nz2sn,m), &
                     nx1sn,nx2sn,ny1sn,ny2sn,nz1sn,nz2sn)
             else
                call mpistop('Output drhoj not implemented yet !',0)
             endif
          case('Mach') ! local Mach number
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      dummy(i,j,k,m)=sqrt(uu(i,j,k)**2+vv(i,j,k)**2+ww(i,j,k)**2)/c_(i,j,k)
                   enddo
                enddo
             enddo
          case('Gamma') ! fundamental derivative of gas dynamics
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      dummy(i,j,k,m)=gcalc_tro(Tmp(i,j,k),rho(i,j,k))
                   enddo
                enddo
             enddo
          case('Frhov') ! Frhov
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      dummy(i,j,k,m)=Frhov(i,j,k)
                   enddo
                enddo
             enddo
          case('Grhow') ! Grhow
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      dummy(i,j,k,m)=Grhow(i,j,k)
                   enddo
                enddo
             enddo
          case('s') ! entropy
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      dummy(i,j,k,m)=scalc_tro(Tmp(i,j,k),rho(i,j,k))
                   enddo
                enddo
             enddo
          case('h') ! enthalpy
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      dummy(i,j,k,m)=ecalc_pro(prs(i,j,k),rho(i,j,k),Tmp(i,j,k))+prs(i,j,k)/rho(i,j,k)
                   enddo
                enddo
             enddo
          case('cfl_i') ! cfl_i
             if (is_curv) then
                do k=nz1sn,nz2sn
                   do j=ny1sn,ny2sn
                      do i=nx1sn,nx2sn
                         ! Calculation of contravariant vitesse ksi + c_sound in grid (ksi,eta)
                         dummy(i,j,k,m) = sqrt(((uu(i,j,k)*y_eta(i,j)-vv(i,j,k)*x_eta(i,j)) * ijacob(i,j))**2)&
                              + c_(i,j,k)*sqrt((y_eta(i,j)**2 + x_eta(i,j)**2) * ijacob(i,j)**2)
                         ! Calculation of cfl_ksi
                         dummy(i,j,k,m) = dummy(i,j,k,m)*deltat
                      enddo
                   enddo
                enddo
             else
                do k=nz1sn,nz2sn
                   do j=ny1sn,ny2sn
                      do i=nx1sn,nx2sn
                         dummy(i,j,k,m)=(abs(uu(i,j,k)) + c_(i,j,k))*deltat*idx(i)
                      enddo
                   enddo
                enddo
             endif
          case('cfl_j') ! cfl_j
             if (is_curv) then
                do k=nz1sn,nz2sn
                   do j=ny1sn,ny2sn
                      do i=nx1sn,nx2sn
                         ! Calculation of contravariant vitesse eta + c_sound in grid (ksi,eta)
                         dummy(i,j,k,m) = sqrt(((vv(i,j,k)*x_ksi(i,j) - uu(i,j,k)*y_ksi(i,j)) * ijacob(i,j))**2) &
                              + c_(i,j,k)*sqrt((y_ksi(i,j)**2 + x_ksi(i,j)**2) * ijacob(i,j)**2)
                         ! Calculation of cfl_eta
                         dummy(i,j,k,m) = dummy(i,j,k,m)*deltat
                      enddo
                   enddo
                enddo
             else
                do k=nz1sn,nz2sn
                   do j=ny1sn,ny2sn
                      do i=nx1sn,nx2sn
                         dummy(i,j,k,m)=(abs(vv(i,j,k)) + c_(i,j,k))*deltat*idy(j)
                      enddo
                   enddo
                enddo
             endif
          case('cfl_ij') ! cfl_ij
             if (is_curv) then
                do k=nz1sn,nz2sn
                   do j=ny1sn,ny2sn
                      do i=nx1sn,nx2sn
                         ! Calculation of contravariant vitesse eta + c_sound in grid (ksi,eta)
                         dummy(i,j,k,m) = sqrt(((vv(i,j,k)*x_ksi(i,j) - uu(i,j,k)*y_ksi(i,j) + &
                                                 uu(i,j,k)*y_eta(i,j) - vv(i,j,k)*x_eta(i,j)) * ijacob(i,j))**2) &
                              + c_(i,j,k)*sqrt((y_ksi(i,j)**2 + x_ksi(i,j)**2 + y_eta(i,j)**2 + x_eta(i,j)**2) * ijacob(i,j)**2)
                         ! Calculation of cfl_eta
                         dummy(i,j,k,m) = dummy(i,j,k,m)*deltat
                      enddo
                   enddo
                enddo
             else
                call mpistop('Output cfl_ij not implemented yet in cartesian !',1)
             endif
          case('cfr_ij') ! cfl_ij
             if (is_curv) then
                do k=nz1sn,nz2sn
                   do j=ny1sn,ny2sn
                      do i=nx1sn,nx2sn
                         ! Calculation of contravariant vitesse eta + c_sound in grid (ksi,eta)
                         dummy(i,j,k,m) = cok(i,j,k)/(rho(i,j,k)*cvcalc_tro(Tmp(i,j,k),rho(i,j,k)))*(g_ksi(i,j)**2+g_eta(i,j)**2)
                         ! Calculation of cfl_ij
                         dummy(i,j,k,m) = dummy(i,j,k,m)*deltat
                      enddo
                   enddo
                enddo
             else
                call mpistop('Output cfl_ij not implemented yet in cartesian or is_curv3 !',1)
             endif
          case( 'QQ') ! Q vorticity criterion
             if (.not.allocated(QQ)) then
                print *,'comp QQ'
                allocate(QQ(nx,ny,nz))
                call compute_Q_criterion
             endif
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      dummy(i,j,k,m)= QQ(i,j,k)
                   enddo
                enddo
             enddo
             deallocate(QQ)
          case( 'vort') ! vorticity norm
             allocate(vort(nx,ny,nz))
             call compute_vorticity3d
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      dummy(i,j,k,m)= vort(i,j,k)
                   enddo
                enddo
             enddo
             deallocate(vort)
          case( 'vortxy') ! vorticity norm
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      dummy(i,j,k,m)=dvx(i,j,k)-duy(i,j,k)
                   enddo
                enddo
             enddo
          case( 'umz') ! streamwise velocity averaged along z
             ! average along z for current proc
             mz=0.0_wp
             do k=1,nz
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      mz(i,j)=mz(i,j)+uu(i,j,k)
                   enddo
                enddo
             enddo
             mz=mz/dble(nz)
             ! average along z for all procs in z
             !call MPI_REDUCE(mz,mzg,nx*ny,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMMXY,info)
             call MPI_ALLREDUCE(MPI_IN_PLACE,mz,nx*ny,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXY,info)
             dummy(:,:,nz1sn,m)=mz/dble(ndomz)
          case( 'Mmz') ! Mach number averaged along z
             ! average along z for current proc
             mz=0.0_wp
             do k=1,nz
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      mz(i,j)=mz(i,j)+sqrt(uu(i,j,k)**2+vv(i,j,k)**2+ww(i,j,k)**2)/c_(i,j,k)
                   enddo
                enddo
             enddo
             mz=mz/dble(nz)
             ! average along z for all procs in z
             !call MPI_REDUCE(mz,mzg,nx*ny,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMMXY,info)
             call MPI_ALLREDUCE(MPI_IN_PLACE,mz,nx*ny,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXY,info)
             dummy(:,:,nz1sn,m)=mz/dble(ndomz)
          case( 'pmz') ! pressure averaged along z
             ! average along z for current proc
             mz=0.0_wp
             do k=1,nz
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      mz(i,j)=mz(i,j)+prs(i,j,k)
                   enddo
                enddo
             enddo
             mz=mz/dble(nz)
             ! average along z for all procs in z
             !call MPI_REDUCE(mz,mzg,nx*ny,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMMXY,info)
             call MPI_ALLREDUCE(MPI_IN_PLACE,mz,nx*ny,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXY,info)
             dummy(:,:,nz1sn,m)=mz/dble(ndomz)
          case( 'rhomz') ! density averaged along z
             ! average along z for current proc
             mz=0.0_wp
             do k=1,nz
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      mz(i,j)=mz(i,j)+rho(i,j,k)
                   enddo
                enddo
             enddo
             mz=mz/dble(nz)
             ! average along z for all procs in z
             !call MPI_REDUCE(mz,mzg,nx*ny,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMMXY,info)
             call MPI_ALLREDUCE(MPI_IN_PLACE,mz,nx*ny,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXY,info)
             dummy(:,:,nz1sn,m)=mz/dble(ndomz)
          case( 'droxmz') ! horizontal density gradient averaged along z
             ! average along z for current proc
             if (is_curv3) then
                call deriv_x_mz_5pts_c3(rho(nx1:nx2,ny1:ny2,nz1:nz2), &
                     mz,nx1sn,nx2sn,ny1sn,ny2sn,1,nz)
             elseif (is_curv) then
                call deriv_x_mz_5pts_c(rho(nx1:nx2,ny1:ny2,nz1:nz2), &
                     mz,nx1sn,nx2sn,ny1sn,ny2sn,1,nz)
             else
                call mpistop('Output drhoi not implemented yet !',0)
             endif
             ! average along z for all procs in z
             !call MPI_REDUCE(mz,mzg,nx*ny,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMMXY,info)
             call MPI_ALLREDUCE(MPI_IN_PLACE,mz,nx*ny,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXY,info)
             dummy(:,:,nz1sn,m)=mz/dble(ndomz)
          case( 'droymz') ! vertical density gradient averaged along z
             ! average along z for current proc
             if (is_curv3) then
                call deriv_y_mz_5pts_c3(rho(nx1:nx2,ny1:ny2,nz1:nz2), &
                     mz,nx1sn,nx2sn,ny1sn,ny2sn,1,nz)
             elseif (is_curv) then
                call deriv_y_mz_5pts_c(rho(nx1:nx2,ny1:ny2,nz1:nz2), &
                     mz,nx1sn,nx2sn,ny1sn,ny2sn,1,nz)
             else
                call mpistop('Output drhoj not implemented yet !',0)
             endif
             ! average along z for all procs in z
             !call MPI_REDUCE(mz,mzg,nx*ny,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMMXY,info)
             call MPI_ALLREDUCE(MPI_IN_PLACE,mz,nx*ny,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXY,info)
             dummy(:,:,nz1sn,m)=mz/dble(ndomz)
          case( 'voxymz') ! xy-vorticity averaged along z                ! average along z for current proc
             mz=0.0_wp
             do k=1,nz
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      mz(i,j)=mz(i,j)+dvx(i,j,k)-duy(i,j,k)
                   enddo
                enddo
             enddo
             mz=mz/dble(nz)
             ! average along z for all procs in z
             !call MPI_REDUCE(mz,mzg,nx*ny,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMMXY,info)
             call MPI_ALLREDUCE(MPI_IN_PLACE,mz,nx*ny,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXY,info)
             dummy(:,:,nz1sn,m)=mz/dble(ndomz)
          case('udf') ! user-defined variable
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      dummy(i,j,k,m)=uvar(i,j,k,cpt)
                   enddo
                enddo
             enddo
             cpt=cpt+1
          case default
             call mpistop('variable '//trim(snapshots(isn)%var(m)) &
                  //' not listed in write_snapshot [mod_io_snapshots.f90]',0)
          end select
       enddo

       ! Fill varlist (pointer on datas)
       ! ===============================
       do m=1,ndata
          varlist(m)%data=>dummy(:,:,:,m)
          varlist(m)%name= dataname(m)
       enddo

       ! Write data [write_tec in mod_tecplot.f90]
       ! ==========
       if ((iproc.eq.0).and.(isn.eq.1)) print *,'Writing snapshot... '//snapshotfile
       filetype=2 ! [0:full; 1:only grid; 2:only solution]
       call write_tec(snapshotfile,filetype,tstar,snapshots(isn)%tectype)

       ! Free temporary memory
       ! =====================
       deallocate(varlist,dataname,dummy)

    else

       ! snapshot indices
       ! ----------------
       nx1sn= 1
       nx2sn= nx
       ny1sn= 1
       ny2sn= ny
       nz1sn= 1
       nz2sn= 1
       ndata=snapshots(isn)%nvar

       ! Fill data in temporary file dummy
       ! =================================
       do m=1,ndata
          select case (trim(snapshots(isn)%var(m)))
          case( 'umz') ! streamwise velocity averaged along z
             ! average along z for current proc
             mz=0.0_wp
             do k=1,nz
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      mz(i,j)=mz(i,j)+uu(i,j,k)
                   enddo
                enddo
             enddo
             mz=mz/dble(nz)
          case( 'Mmz') ! Mach number averaged along z
             ! average along z for current proc
             mz=0.0_wp
             do k=1,nz
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      mz(i,j)=mz(i,j)+sqrt(uu(i,j,k)**2+vv(i,j,k)**2+ww(i,j,k)**2)/c_(i,j,k)
                   enddo
                enddo
             enddo
             mz=mz/dble(nz)
          case( 'pmz') ! pressure averaged along z
             ! average along z for current proc
             mz=0.0_wp
             do k=1,nz
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      mz(i,j)=mz(i,j)+prs(i,j,k)
                   enddo
                enddo
             enddo
             mz=mz/dble(nz)
          case( 'rhomz') ! density averaged along z
             ! average along z for current proc
             mz=0.0_wp
             do k=1,nz
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      mz(i,j)=mz(i,j)+rho(i,j,k)
                   enddo
                enddo
             enddo
             mz=mz/dble(nz)
          case( 'droxmz') ! horizontal density gradient averaged along z
             ! average along z for current proc
             if (is_curv3) then
                call deriv_x_mz_5pts_c3(rho(nx1:nx2,ny1:ny2,nz1:nz2), &
                     mz,nx1sn,nx2sn,ny1sn,ny2sn,1,nz)
             elseif (is_curv) then
                call deriv_x_mz_5pts_c(rho(nx1:nx2,ny1:ny2,nz1:nz2), &
                     mz,nx1sn,nx2sn,ny1sn,ny2sn,1,nz)
             else
                call mpistop('Output drhoi not implemented yet !',0)
             endif
          case( 'droymz') ! vertical density gradient averaged along z
             ! average along z for current proc
             if (is_curv3) then
                call deriv_y_mz_5pts_c3(rho(nx1:nx2,ny1:ny2,nz1:nz2), &
                     mz,nx1sn,nx2sn,ny1sn,ny2sn,1,nz)
             elseif (is_curv) then
                call deriv_y_mz_5pts_c(rho(nx1:nx2,ny1:ny2,nz1:nz2), &
                     mz,nx1sn,nx2sn,ny1sn,ny2sn,1,nz)
             else
                call mpistop('Output drhoj not implemented yet !',0)
             endif
          case( 'voxymz') ! xy-vorticity averaged along z
             ! average along z for current proc
             mz=0.0_wp
             do k=1,nz
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      mz(i,j)=mz(i,j)+dvx(i,j,k)-duy(i,j,k)
                   enddo
                enddo
             enddo
             mz=mz/dble(nz)
          end select
          
          ! Average along z for all procs in z
          ! ----------------------------------
          call MPI_ALLREDUCE(MPI_IN_PLACE,mz,nx*ny,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXY,info)
          !call MPI_REDUCE(mz,mzg,nx*ny,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMMXY,info)
       enddo

    endif
    
  end subroutine write_snapshot

  !==============================================================================================
  subroutine read_snapshot(isn,ibl2read,dirDATA)
  !==============================================================================================
    !> authors: XG & AB
    !> date: June 2023
    !> Read variables for snapshot # isn in directory dirDATA
  !==============================================================================================
    use mod_utils
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Inputs: index of snapshot to be written
    integer, intent(in) :: isn,ibl2read        ! index of snapshot to be written
    character(len=60), intent(in) :: dirDATA   ! directory where snapshot is stored
    ! -------------------------------------------------------------------------------------------
    ! Local variables
    integer :: i,j,k,m,ndata,cpt
    integer :: nx1sn,nx2sn,ny1sn,ny2sn,nz1sn,nz2sn
    ! -------------------------------------------------------------------------------------------
    character(:), allocatable :: snapshotfile
    ! -------------------------------------------------------------------------------------------

    ! Check if the snapshot belong to this proc
    ! -----------------------------------------
    if (.not.(snapshots(isn)%has_snap)) return

    ! Determine snapshot name
    ! =======================
    if (snapshots(isn)%stamp) then
       ! name of the snapshot with timestamp
       ! --------------------------------
       if (.not.snapshots(isn)%tectype%is_IOtec_read) then
          snapshotfile=trim(snapshots(isn)%type_name)//'_'//snapshots(isn)%tectype%zonename(6:8)//'_bl_'//filestamp//trim(numchar(ibl2read))//'.bin'
       else
          snapshotfile=trim(snapshots(isn)%type_name)//'_'//snapshots(isn)%tectype%zonename(6:8)//'_bl_'//filestamp//trim(numchar(ibl2read))//'.plt'
       endif
    else
       ! name of the snapshot without timestamp
       ! -----------------------------------
       if (.not.snapshots(isn)%tectype%is_IOtec_read) then
          snapshotfile=trim(snapshots(isn)%type_name)//'_'//snapshots(isn)%tectype%zonename(6:8)//'_bl'//trim(numchar(ibl2read))//'.bin'
       else
          snapshotfile=trim(snapshots(isn)%type_name)//'_'//snapshots(isn)%tectype%zonename(6:8)//'_bl'//trim(numchar(ibl2read))//'.plt'
       endif
    endif

    ! Prepare number and name of data
    ! ===============================
    ndata=snapshots(isn)%nvar
    allocate(varlist(ndata),dataname(ndata))
    do m=1,ndata
       dataname(m)=trim(snapshots(isn)%var(m))
    enddo

    ! Check file existence
    ! ====================
    call mpicheckfile(trim(dirDATA)//snapshotfile)

    ! Plane indices
    ! =============
    nx1sn= snapshots(isn)%tectype%nx1
    nx2sn= snapshots(isn)%tectype%nx2
    ny1sn= snapshots(isn)%tectype%ny1
    ny2sn= snapshots(isn)%tectype%ny2
    nz1sn= snapshots(isn)%tectype%nz1
    nz2sn= snapshots(isn)%tectype%nz2

    ! allocate dummy variable for datas
    ! ---------------------------------
    allocate(dummy(nx1sn:nx2sn,ny1sn:ny2sn,nz1sn:nz2sn,ndata))

    ! fill varlist (pointer on datas)
    ! -------------------------------
    do m=1,ndata
       varlist(m)%data=>dummy(:,:,:,m)
       varlist(m)%name= dataname(m)
    enddo

    ! Read data [read_tec in mod_tecplot.f90]
    ! =========
    !if ((iproc.eq.0).and.(isn.eq.1)) print *,'Reading snapshot... '//trim(dirDATA)//snapshotfile
    call read_tec(trim(dirDATA)//snapshotfile,snapshots(isn)%tectype)
    call MPI_BARRIER(MPI_COMM_WORLD,info)

    ! Fill data in temporary file dummy
    ! =================================
    ! List of possible variables: prs,uu,vv,ww,rho,Tmp,div,Mach,Gamma,Frhov,Grhow,udf

    ! Initialize counter for user-defined (udf) variables
    cpt=1

    do m=1,ndata
       select case (trim(snapshots(isn)%var(m)))
       case('prs') ! pressure
          if (allocated(prs)) then
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      prs(i,j,k)=dummy(i,j,k,m)
                   enddo
                enddo
             enddo
          endif
       case( 'uu') ! streamwise velocity
          if (allocated(uu)) then
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      uu(i,j,k)=dummy(i,j,k,m)
                   enddo
                enddo
             enddo
          endif
       case( 'vv') ! crossflow velocity
          if (allocated(vv)) then
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      vv(i,j,k)=dummy(i,j,k,m)
                   enddo
                enddo
             enddo
          endif
       case( 'ww') ! spanwise velocity
          if (allocated(ww)) then
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      ww(i,j,k)=dummy(i,j,k,m)
                   enddo
                enddo
             enddo
          endif
       case('rho') ! density
          if (allocated(rho)) then
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      rho(i,j,k)=dummy(i,j,k,m)
                   enddo
                enddo
             enddo
          endif
       case('Tmp') ! temperature
          if (allocated(Tmp)) then
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      Tmp(i,j,k)=dummy(i,j,k,m)
                   enddo
                enddo
             enddo
          endif
       case('div') ! divergence
          call mpistop('no array for div',0)
!!$          if (allocated(div)) then
!!$             do k=nz1sn,nz2sn
!!$                do j=ny1sn,ny2sn
!!$                   do i=nx1sn,nx2sn
!!$                      div(i,j,k)=dummy(i,j,k,m)
!!$                   enddo
!!$                enddo
!!$             enddo
!!$          endif
       case('Mach') ! local Mach number
          call mpistop('no array for Mach',0)
!!$          if (allocated(Mach_)) then
!!$             do k=nz1sn,nz2sn
!!$                do j=ny1sn,ny2sn
!!$                   do i=nx1sn,nx2sn
!!$                      Mach_(i,j,k)=dummy(i,j,k,m)
!!$                   enddo
!!$                enddo
!!$             enddo
!!$          endif
       case('Gamma') ! fundamental derivative of gas dynamics
          call mpistop('no array for Gamma',0)
!!$          if (allocated(Gamma_)) then
!!$             do k=nz1sn,nz2sn
!!$                do j=ny1sn,ny2sn
!!$                   do i=nx1sn,nx2sn
!!$                      Gamma_(i,j,k)=dummy(i,j,k,m)
!!$                   enddo
!!$                enddo
!!$             enddo
!!$          endif
       case('Frhou') ! Frhou
          if (allocated(Frhou)) then
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      Frhou(i,j,k)=dummy(i,j,k,m)
                   enddo
                enddo
             enddo
          endif
       case('Frhov') ! Frhov
          if (allocated(Frhov)) then
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      Frhov(i,j,k)=dummy(i,j,k,m)
                   enddo
                enddo
             enddo
          endif
       case('Frhow') ! Frhow
          if (allocated(Frhow)) then
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      Frhow(i,j,k)=dummy(i,j,k,m)
                   enddo
                enddo
             enddo
          endif
       case('Grhov') ! Grhov
          if (allocated(Grhov)) then
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      Grhov(i,j,k)=dummy(i,j,k,m)
                   enddo
                enddo
             enddo
          endif
       case('Grhow') ! Grhow
          if (allocated(Grhow)) then
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      Grhow(i,j,k)=dummy(i,j,k,m)
                   enddo
                enddo
             enddo
          endif
       case('Hrhow') ! Hrhow
          if (allocated(Hrhow)) then
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      Hrhow(i,j,k)=dummy(i,j,k,m)
                   enddo
                enddo
             enddo
          endif
       case('QQ') ! Q criterion
          if (allocated(QQ)) then
             do k=nz1sn,nz2sn
                do j=ny1sn,ny2sn
                   do i=nx1sn,nx2sn
                      QQ(i,j,k)=dummy(i,j,k,m)
                   enddo
                enddo
             enddo
          endif
       case('udf') ! user-defined variable
          call mpistop('no array for udf',0)
!!$          if (allocated(uvar)) then
!!$             do k=nz1sn,nz2sn
!!$                do j=ny1sn,ny2sn
!!$                   do i=nx1sn,nx2sn
!!$                      uvar(i,j,k,cpt)=dummy(i,j,k,m)
!!$                   enddo
!!$                enddo
!!$             enddo
!!$          endif
!!$          cpt=cpt+1
       case default
          call mpistop('variable '//trim(snapshots(isn)%var(m)) &
               //' not listed in write_snapshot [mod_io_snapshots.f90]',0)
       end select
    enddo

! TO BE CHANGED jva not implemented

!!$       ! fill data Not generic -> TO BE CHANGED  depend on the number of data !!
!!$       if (snapshots(ip)%plnorm.ne.ijva) then ! for non-jva snapshots first
!!$
!!$       else                                ! for jva snapshot
!!$          do k=nz1sn,nz2sn
!!$             do j=ny1sn,ny2sn
!!$                do i=nx1sn,nx2sn
!!$                   jj=snapshots(ip)%jvar(i)
!!$                   rho(i,jj,k)= dummy(i,j,k,1)
!!$                   div(i,jj,k)= dummy(i,j,k,2)
!!$                   prs(i,jj,k)= dummy(i,j,k,3)
!!$                   uu(i,jj,k) = dummy(i,j,k,4)
!!$                   vv(i,jj,k) = dummy(i,j,k,5)
!!$                   ww(i,jj,k) = dummy(i,j,k,6)
!!$                enddo
!!$             enddo
!!$          enddo
!!$       endif

    ! Free temporary memory
    ! =====================
    deallocate(varlist,dataname,dummy)

  end subroutine read_snapshot

  !==============================================================================================
  subroutine free_comm_io_snapshot(MPI_COMM_SNAPSHOT)
  !==============================================================================================
    !> author: Luca Sciacovelli
    !> date: April 2018 (modif XG June 2021 - AB May 2022)
    !> Free sub-communicators for handling I/O snapshots
  !==============================================================================================
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Input/Output arguments
    integer, intent(inout) :: MPI_COMM_SNAPSHOT
    ! -------------------------------------------------------------------------------------------

    ! Destroy snapshot sub-comminicators
    ! ==================================
    if (MPI_COMM_SNAPSHOT/=MPI_COMM_NULL) call MPI_COMM_FREE(MPI_COMM_SNAPSHOT,info)

  end subroutine free_comm_io_snapshot

  ! !==============================================================================================
  ! subroutine read_loc_max(jm_urms)
  ! !==============================================================================================
  !   !> author: XG
  !   !> date: April 2018
  !   !> read loc_max_urms.dat: index j for max(urms) / modal energy analysis
  ! !==============================================================================================
  !   use warnstop
  !   implicit none
  !   ! -------------------------------------------------------------------------------------------
  !   ! Input/Output arguments
  !   integer, dimension(nx), intent(out) :: jm_urms
  !   ! -------------------------------------------------------------------------------------------
  !   ! Local variables
  !   integer :: i
  !   logical :: iexist
  !   integer, dimension(ngx) :: jmg_urms
  !   ! -------------------------------------------------------------------------------------------

  !   ! Open loc_max_urms.dat
  !   ! =====================
  !   inquire(file='loc_max_urms.dat', exist=iexist)
  !   if (.not.iexist) then
  !      call mpistop('loc_max_urms.dat does not exist!', 0)
  !   endif

  !   ! Read
  !   ! ====
  !   open(31,file='loc_max_urms.dat')
  !   jm_urms=1
  !   do i=1,ngx
  !      read(31,*) jmg_urms(i)
  !   enddo
  !   close(31)

  !   ! MPI partitioning
  !   ! ================
  !   do i=1,nx
  !      jm_urms(i)=jmg_urms(i+coord(1)*nx)
  !   enddo

  ! end subroutine read_loc_max


end module mod_io_snapshots
