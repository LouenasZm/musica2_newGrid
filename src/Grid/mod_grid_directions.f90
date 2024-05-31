!===============================================================================
module mod_grid_directions
!===============================================================================
  !> Module to search swap and reverse directions between domains or blocks
!===============================================================================
  use mod_grid
  use mod_block
  use mod_mpi
  implicit none
  ! ----------------------------------------------------------------------------
  integer, private :: nbl ! current block
  ! ----------------------------------------------------------------------------
  ! Indicator at domain level:
  ! ==========================
  ! indicators for swap i-j directions
  logical, dimension(6) :: is_swapij,is_swapij2
  ! indicators for swap i-k directions
  logical, dimension(6) :: is_swapik,is_swapik2
  ! indicators for swap j-k directions
  logical, dimension(6) :: is_swapjk,is_swapjk2
  ! indicators for swap i-j-k directions
  logical, dimension(6) :: is_swapijk,is_swapijk2
  ! indicators for reverse directions
  logical, dimension(6,3) :: is_rev,is_rev2
  ! directions of neighboring faces
  integer, dimension(6) :: ndir
  ! ----------------------------------------------------------------------------
  ! Indicator at block level:
  ! =========================
  ! indicators for swap i-j block directions
  logical, dimension(6) :: is_swapij_bl,is_swapij2_bl
  ! indicators for swap i-k block directions
  logical, dimension(6) :: is_swapik_bl,is_swapik2_bl
  ! indicators for swap j-k block directions
  logical, dimension(6) :: is_swapjk_bl,is_swapjk2_bl
  ! indicators for swap i-j-k block directions
  logical, dimension(6) :: is_swapijk_bl,is_swapijk2_bl
  ! indicators for reverse block directions
  logical, dimension(6,3) :: is_rev_bl,is_rev2_bl
  ! directions of neighboring faces for blocks
  integer, dimension(6) :: ndir_bl
  ! ----------------------------------------------------------------------------

contains
  
  !===============================================================
  subroutine grid_directions(neighbors)
  !===============================================================
    !> Determine ij-swap or reverse directions between blocks
  !===============================================================
    use warnstop
    implicit none
    ! ------------------------------------------------------------
    integer, dimension(4), intent(in) :: neighbors    
    ! ------------------------------------------------------------
    integer :: n,nv,k,m,ip
    integer :: n_nb                  ! number of tested neighbors
    integer, dimension(4) :: list_nb ! list of tested neighbors
    integer, dimension(:,:), allocatable :: iproc_corner
    ! ------------------------------------------------------------
    ! recreate coord1,coord2 -> TO BE CHANGED
    integer, dimension(:), allocatable :: coord1,coord2
    
    ! Initializations
    ! ===============
    
    ! Initialize is_swapij (swap i-j directions)
    ! ------------------
    is_swapij=.false.
    is_swapij2=.false.
    is_swapij_bl=.false.
    is_swapij2_bl=.false.

    ! Initialize is_swapik (swap i-k directions)
    ! ------------------
    is_swapik=.false.
    is_swapik2=.false.
    is_swapik_bl=.false.
    is_swapik2_bl=.false.

    ! Initialize is_swapjk (swap j-k directions)
    ! ------------------
    is_swapjk=.false.
    is_swapjk2=.false.
    is_swapjk_bl=.false.
    is_swapjk2_bl=.false.
    
    ! Initialize is_rev (reverse direction)
    ! -----------------
    is_rev=.false.
    is_rev2=.false.
    is_rev_bl=.false.
    is_rev2_bl=.false.
    
    ! Initialize array of neighboring directions
    ! ------------------------------------------
    ndir(1)=2 ! by default 1/W exchange with 2/E
    ndir(2)=1 ! by default 2/E exchange with 1/W
    ndir(3)=4 ! by default 3/S exchange with 4/N
    ndir(4)=3 ! by default 4/N exchange with 3/S
    ndir(5)=6 ! by default 5/F exchange with 6/B
    ndir(6)=5 ! by default 6/B exchange with 5/F
    
    ! Initialize array of neighboring directions for blocks
    ! -----------------------------------------------------
    ndir_bl(1)=2 ! by default 1/W exchange with 2/E
    ndir_bl(2)=1 ! by default 2/E exchange with 1/W
    ndir_bl(3)=4 ! by default 3/S exchange with 4/N
    ndir_bl(4)=3 ! by default 4/N exchange with 3/S
    ndir_bl(5)=6 ! by default 5/F exchange with 6/B
    ndir_bl(6)=5 ! by default 6/B exchange with 5/F

    if (.not.is_curv) return

    ! Current block
    ! =============
    nbl=nob(iproc)
    
    ! Determine block corners
    ! =======================
    !  3----------------4
    !  |                |
    !  |    ^           |
    !  |    | j         |
    !  |    |           |
    !  |    ----->      |
    !  |     i          |
    !  |                |
    !  1----------------2
    
    ! ! Use global grid <~ To be suppressed
    ! ! ---------------
    ! bl(nbl)%xcorner(1)=xgc(1,1)
    ! bl(nbl)%xcorner(2)=xgc(ngx,1)
    ! bl(nbl)%xcorner(3)=xgc(1,ngy)
    ! bl(nbl)%xcorner(4)=xgc(ngx,ngy)

    ! bl(nbl)%ycorner(1)=ygc(1,1)
    ! bl(nbl)%ycorner(2)=ygc(ngx,1)
    ! bl(nbl)%ycorner(3)=ygc(1,ngy)
    ! bl(nbl)%ycorner(4)=ygc(ngx,ngy)

    ! call MPI_BARRIER(COMM_global,info)

    ! ! share block corner coordinate between procs
    ! do n=1,nbloc
    !    call MPI_BCAST(bl(n)%xcorner,4,MPI_DOUBLE_PRECISION,iproc_leader(n),COMM_global,info)
    !    call MPI_BCAST(bl(n)%ycorner,4,MPI_DOUBLE_PRECISION,iproc_leader(n),COMM_global,info)
    ! enddo
    ! ! ---------------

    ! Use local grid <~ new way
    ! --------------
    ! gather coord in the MPI_CART communicator
    allocate(coord1(0:nproc-1),coord2(0:nproc-1))
    call MPI_ALLGATHER(coord(1),1,MPI_INTEGER,coord1,1,MPI_INTEGER,COMM_global,info)
    call MPI_ALLGATHER(coord(2),1,MPI_INTEGER,coord2,1,MPI_INTEGER,COMM_global,info)

    ! determine who has corners
    allocate(iproc_corner(nbloc,4))

    do ip=0,nproc-1
       n=nob(ip)
       if ((coord1(ip)==0).and.(coord2(ip)==0))                         iproc_corner(n,1)=ip
       if ((coord1(ip)==bl(n)%ndomi-1).and.(coord2(ip)==0))             iproc_corner(n,2)=ip
       if ((coord1(ip)==0).and.(coord2(ip)==bl(n)%ndomj-1))             iproc_corner(n,3)=ip
       if ((coord1(ip)==bl(n)%ndomi-1).and.(coord2(ip)==bl(n)%ndomj-1)) iproc_corner(n,4)=ip
    enddo

    if ((coord(1)==0).and.(coord(2)==0)) then
       bl(nbl)%xcorner(1)=xc(1,1)
       bl(nbl)%ycorner(1)=yc(1,1)
    endif
    if ((coord(1)==ndomx-1).and.(coord(2)==0)) then
       bl(nbl)%xcorner(2)=xc(nx,1)
       bl(nbl)%ycorner(2)=yc(nx,1)
    endif
    if ((coord(1)==0).and.(coord(2)==ndomy-1)) then
       bl(nbl)%xcorner(3)=xc(1,ny)
       bl(nbl)%ycorner(3)=yc(1,ny)
    endif
    if ((coord(1)==ndomx-1).and.(coord(2)==ndomy-1)) then
       bl(nbl)%xcorner(4)=xc(nx,ny)
       bl(nbl)%ycorner(4)=yc(nx,ny)
    endif

    do n=1,nbloc
       do m=1,4
          call MPI_BCAST(bl(n)%xcorner(m),1,MPI_DOUBLE_PRECISION,iproc_corner(n,m),COMM_global,info)
          call MPI_BCAST(bl(n)%ycorner(m),1,MPI_DOUBLE_PRECISION,iproc_corner(n,m),COMM_global,info)
       enddo
    enddo
    ! --------------

    ! Current block
    ! =============
    n=nob(iproc)
    
100 format(1x,'iproc',i4,', list:',4i4,', n_nb ',i4)

    ! Search for swap and reverse directions for domains with lower rank
    ! ==================================================================
    
    ! Identify neighboring domains
    ! ----------------------------
    list_nb=0
    n_nb=0
    do k=1,4
       nv=bl(n)%BC(k)
       ! nv>0: connectivity with a domain, nv<n: lower rank,
       ! neighbor(k)==MPI_PROC_NULL means that the domain has no neighbor for the tested face
       ! because it is not yet determined => it belongs to another block (not intrablock)
       ! nv.ne.n: exclude comm with itself,
       if ((nv>0).and.(nv<n).and.(neighbors(k)==MPI_PROC_NULL).and.(nv.ne.n)) then
          list_nb(k)=bl(n)%BC(k)
          n_nb=n_nb+1
       endif
    enddo
    !print 100,iproc,list_nb(1:4),n_nb

    ! Determine indicators
    ! --------------------
    ! NOTA: ndir not used
    if (n_nb>0) call check_directions(list_nb,is_swapij,is_rev(1:4,1),is_rev(1:4,2),ndir)

    ! Search for swap and reverse directions for all domains
    ! ======================================================

    ! Identify neighboring domains
    ! ----------------------------
    list_nb=0
    n_nb=0
    do k=1,4
       nv=bl(n)%BC(k)
       ! nv>0: connectivity with a domain,
       ! neighbor(k)==MPI_PROC_NULL means that the domain has no neighbor for the tested face
       ! because it is not yet determined => it belongs to another block (not intrablock),
       ! nv.ne.n: exclude comm with itself
       if ((nv>0).and.(neighbors(k)==MPI_PROC_NULL).and.(nv.ne.n)) then
          list_nb(k)=bl(n)%BC(k)
          n_nb=n_nb+1
       endif
    enddo   
    !print 100,iproc,list_nb(1:4),n_nb

    ! Determine indicators
    ! --------------------
    ! NOTA: ndir2 is ndir
    if (n_nb>0) call check_directions(list_nb,is_swapij2,is_rev2(1:4,1),is_rev2(1:4,2),ndir)

    ! Search for swap and reverse directions for blocks with lower rank
    ! =================================================================

    ! Identify neighboring domains
    ! ----------------------------
    list_nb=0
    n_nb=0
    do k=1,4
       nv=bl(n)%BC(k)
       ! nv>0: connectivity with a block, nv<n: lower rank,
       ! nv.ne.n: exclude comm with itself
       if ((nv>0).and.(nv<n).and.(nv.ne.n)) then
          list_nb(k)=bl(n)%BC(k)
          n_nb=n_nb+1
       endif
    enddo    
    !print 100,iproc,list_nb(1:4),n_nb

    ! Determine indicators
    ! --------------------
    ! NOTA: ndir_bl not used; is_rev_bl(:,2) not used
    if (n_nb>0) call check_directions(list_nb,is_swapij_bl,is_rev_bl(1:4,1),is_rev_bl(1:4,2),ndir_bl)

    ! Search for swap and reverse directions for all blocks
    ! =====================================================

    ! Identify neighboring domains
    ! ----------------------------
    list_nb=0
    n_nb=0
    do k=1,4
       nv=bl(n)%BC(k)
       ! nv>0: connectivity with a block,
       ! nv.ne.n: exclude comm with itself
       if ((nv>0).and.(nv.ne.n)) then
          list_nb(k)=bl(n)%BC(k)
          n_nb=n_nb+1
       endif
    enddo  
    !print 100,iproc,list_nb(1:4),n_nb
      
    ! Determine indicators
    ! --------------------
    ! NOTA: ndir2_bl is ndir_bl
    if (n_nb>0) call check_directions(list_nb,is_swapij2_bl,is_rev2_bl(1:4,1),is_rev2_bl(1:4,2),ndir_bl)

    ! if (iproc==1) then
    ! print *,iproc,'is_swapij_bl',is_swapij_bl,'is_swapij2_bl',is_swapij2_bl
    ! print *,iproc,' is_rev_p',is_rev_bl(:,1),'is_rev2_p',is_rev2_bl(:,1)
    ! print *,iproc,' is_rev_n',is_rev_bl(:,2),'is_rev2_n',is_rev2_bl(:,2)
    ! endif
 
    ! call mpistop('pause check',0)

  end subroutine grid_directions
  
  !=======================================================================
  subroutine check_directions(list_nb, &
                              iswap_ij,irev_parallel,irev_normal,n_dir)
  !=======================================================================
    !> Subroutine to determine indicators swap and reverse
  !=======================================================================
    use warnstop
    implicit none
    ! --------------------------------------------------------------------
    ! INPUT: list of tested neighbors
    integer, dimension(4), intent(in) :: list_nb
    ! OUTPUT: indicator for swap between i and j directions
    logical, dimension(4), intent(inout)  :: iswap_ij
    ! OUTPUT: indicator for reverse direction (parallel and normal to face)
    logical, dimension(4), intent(inout)  :: irev_parallel,irev_normal
    ! OUTPUT: array of neighboring directions
    integer, dimension(6), intent(inout)  :: n_dir
    ! --------------------------------------------------------------------
    integer :: k  ! tested direction [1:imin;2:imax;3:jmin;4:jmax]
    integer :: nv ! number of neighboring block
    integer :: nb ! number of current block
    ! integer :: iA,jA,iB,jB,iC,jC          ! point indices <~ To be suppressed
    integer :: imA,imB                    ! point indices
    real(wp) :: xA,yA,xB,yB,xC,yC         ! points
    real(wp), dimension(4) :: distA,distB ! distances from corners
    real(wp) :: distC_imA,distA_C         ! distances in normal direction
    ! --------------------------------------------------------------------

    ! Current block
    ! =============
    nb=nob(iproc)
   
    ! Check directions
    ! ================
    do k=1,4 ! for each BC (2D here)

       if (list_nb(k)>0) then

          nv=list_nb(k)

          select case(k)

          ! ! Based on global grid <~ To be suppressed
          ! ! --------------------
          ! case (1)  ! Face imin/1/W
          !    iA=1   ! =============
          !    jA=1   !          B
          !           !       ^  | (iA,jA)=(1,1)
          !    iB=1   !     j |  |
          !    jB=ngy !       |  | (iB,jB)=(1,ngy)
          !           !  <-----  |
          !    iC=ngx !     i    | (iC,jC)=(ngx,1)
          !    jC=1   ! C--------A
            
          ! case (2)  ! Face imax/2/E
          !    iA=ngx ! =============
          !    jA=1   !          B
          !           !   ^      | (iA,jA)=(ngx,1)
          !    iB=ngx !   | j    |
          !    jB=ngy !   |      | (iB,jB)=(ngx,ngy)
          !           !   ---->  |
          !    iC=1   !     i    | (iC,jC)=(1,1)
          !    jC=1   ! C--------A
                                         
          ! case (3)  ! Face jmin/3/S
          !    iA=1   ! =============
          !    jA=1   !          B
          !           !       ^  | (iA,jA)=(1,1)
          !    iB=ngx !     i |  |
          !    jB=1   !       |  | (iB,jB)=(ngx,1)
          !           !  <-----  |
          !    iC=1   !     j    | (iC,jC)=(1,ngy)
          !    jC=ngy ! C--------A
                          
          ! case (4)  ! Face jmax/4/N
          !    iA=1   ! =============
          !    jA=ngy !          B
          !           !   ^      | (iA,jA)=(1,ngy)
          !    iB=ngx !   | i    |
          !    jB=ngy !   |      | (iB,jB)=(ngx,ngy)
          !           !   -----> |
          !    iC=1   !     j    | (iC,jC)=(1,1)
          !    jC=1   ! C--------A
             
          ! end select

          ! ! distance A-corner1
          ! distA(1)=(xgc(iA,jA)-bl(nv)%xcorner(1))**2 &
          !         +(ygc(iA,jA)-bl(nv)%ycorner(1))**2
          ! ! distance A-corner2
          ! distA(2)=(xgc(iA,jA)-bl(nv)%xcorner(2))**2 &
          !         +(ygc(iA,jA)-bl(nv)%ycorner(2))**2
          ! ! distance A-corner3
          ! distA(3)=(xgc(iA,jA)-bl(nv)%xcorner(3))**2 &
          !         +(ygc(iA,jA)-bl(nv)%ycorner(3))**2
          ! ! distance A-corner4
          ! distA(4)=(xgc(iA,jA)-bl(nv)%xcorner(4))**2 &
          !         +(ygc(iA,jA)-bl(nv)%ycorner(4))**2
          ! ! closest corner to A
          ! imA=minloc(distA,1)

          ! ! distance B-corner1
          ! distB(1)=(xgc(iB,jB)-bl(nv)%xcorner(1))**2 &
          !         +(ygc(iB,jB)-bl(nv)%ycorner(1))**2
          ! ! distance B-corner2
          ! distB(2)=(xgc(iB,jB)-bl(nv)%xcorner(2))**2 &
          !         +(ygc(iB,jB)-bl(nv)%ycorner(2))**2
          ! ! distance B-corner3
          ! distB(3)=(xgc(iB,jB)-bl(nv)%xcorner(3))**2 &
          !         +(ygc(iB,jB)-bl(nv)%ycorner(3))**2
          ! ! distance B-corner4
          ! distB(4)=(xgc(iB,jB)-bl(nv)%xcorner(4))**2 &
          !         +(ygc(iB,jB)-bl(nv)%ycorner(4))**2
          ! ! closest corner to B
          ! imB=minloc(distB,1)

          ! ! distance C-imA (closest corner to A)
          ! ! [used to check rev_normal]
          ! distC_imA=(xgc(iC,jC)-bl(nv)%xcorner(imA))**2 &
          !          +(ygc(iC,jC)-bl(nv)%ycorner(imA))**2

          ! ! distance A-C
          ! ! [used to check periodicity]
          ! distA_C=(xgc(iC,jC)-xgc(iA,jA))**2 &
          !        +(ygc(iC,jC)-ygc(iA,jA))**2
          ! ! --------------------


          ! Based on local grid <~ new way
          ! -------------------
          case (1)  ! Face imin/1/W
             ! iA=1   ! =============
             ! jA=1   !          B
             !        !       ^  | (iA,jA)=(1,1)
             ! iB=1   !     j |  |
             ! jB=ngy !       |  | (iB,jB)=(1,ngy)
             !        !  <-----  |
             ! iC=ngx !     i    | (iC,jC)=(ngx,1)
             ! jC=1   ! C--------A
             xA=bl(nb)%xcorner(1)
             yA=bl(nb)%ycorner(1)
             xB=bl(nb)%xcorner(3)
             yB=bl(nb)%ycorner(3)
             xC=bl(nb)%xcorner(2)
             yC=bl(nb)%ycorner(2)

          case (2)  ! Face imax/2/E
             ! iA=ngx ! =============
             ! jA=1   !          B
             !        !   ^      | (iA,jA)=(ngx,1)
             ! iB=ngx !   | j    |
             ! jB=ngy !   |      | (iB,jB)=(ngx,ngy)
             !        !   ---->  |
             ! iC=1   !     i    | (iC,jC)=(1,1)
             ! jC=1   ! C--------A
             xA=bl(nb)%xcorner(2)
             yA=bl(nb)%ycorner(2)
             xB=bl(nb)%xcorner(4)
             yB=bl(nb)%ycorner(4)
             xC=bl(nb)%xcorner(1)
             yC=bl(nb)%ycorner(1)

          case (3)  ! Face jmin/3/S
             ! iA=1   ! =============
             ! jA=1   !          B
             !        !       ^  | (iA,jA)=(1,1)
             ! iB=ngx !     i |  |
             ! jB=1   !       |  | (iB,jB)=(ngx,1)
             !        !  <-----  |
             ! iC=1   !     j    | (iC,jC)=(1,ngy)
             ! jC=ngy ! C--------A
             xA=bl(nb)%xcorner(1)
             yA=bl(nb)%ycorner(1)
             xB=bl(nb)%xcorner(2)
             yB=bl(nb)%ycorner(2)
             xC=bl(nb)%xcorner(3)
             yC=bl(nb)%ycorner(3)

          case (4)  ! Face jmax/4/N
             ! iA=1   ! =============
             ! jA=ngy !          B
             !        !   ^      | (iA,jA)=(1,ngy)
             ! iB=ngx !   | i    |
             ! jB=ngy !   |      | (iB,jB)=(ngx,ngy)
             !        !   -----> |
             ! iC=1   !     j    | (iC,jC)=(1,1)
             ! jC=1   ! C--------A
             xA=bl(nb)%xcorner(3)
             yA=bl(nb)%ycorner(3)
             xB=bl(nb)%xcorner(4)
             yB=bl(nb)%ycorner(4)
             xC=bl(nb)%xcorner(1)
             yC=bl(nb)%ycorner(1)

          end select

          ! distance A-corner1
          distA(1)=(xA-bl(nv)%xcorner(1))**2 &
                  +(yA-bl(nv)%ycorner(1))**2
          ! distance A-corner2
          distA(2)=(xA-bl(nv)%xcorner(2))**2 &
                  +(yA-bl(nv)%ycorner(2))**2
          ! distance A-corner3
          distA(3)=(xA-bl(nv)%xcorner(3))**2 &
                  +(yA-bl(nv)%ycorner(3))**2
          ! distance A-corner4
          distA(4)=(xA-bl(nv)%xcorner(4))**2 &
                  +(yA-bl(nv)%ycorner(4))**2
          ! closest corner to A
          imA=minloc(distA,1)

          ! distance B-corner1
          distB(1)=(xB-bl(nv)%xcorner(1))**2 &
                  +(yB-bl(nv)%ycorner(1))**2
          ! distance B-corner2
          distB(2)=(xB-bl(nv)%xcorner(2))**2 &
                  +(yB-bl(nv)%ycorner(2))**2
          ! distance B-corner3
          distB(3)=(xB-bl(nv)%xcorner(3))**2 &
                  +(yB-bl(nv)%ycorner(3))**2
          ! distance B-corner4
          distB(4)=(xB-bl(nv)%xcorner(4))**2 &
                  +(yB-bl(nv)%ycorner(4))**2
          ! closest corner to B
          imB=minloc(distB,1)

          ! distance C-imA (closest corner to A)
          ! [used to check rev_normal]
          distC_imA=(xC-bl(nv)%xcorner(imA))**2 &
                   +(yC-bl(nv)%ycorner(imA))**2

          ! distance A-C
          ! [used to check periodicity]
          distA_C=(xC-xA)**2 &
                 +(yC-yA)**2
          ! -------------------



          
          ! Case #1: ex pour imax (k=2)
          !              B   3----------------4
          !   ^          |   |                |  distance min between A and corners of neighboring domain
          !   | j        |   |    ^           |   index imA=1
          !   |          |   |    | j         |  distance min between B and corners of neighboring domain
          !   ----->     |   |    |           |   index imB=3
          !     i        |   |    ----->      |   is_swap=.false.
          !              |   |     i          |   is_rev(k,1)=.false. for parallel direction
          !              |   |                |   is_rev(k,2)=.false. for perpendicular direction
          !  C-----------A   1----------------2                       [dist(A,imA)<dist(C,imA)]
          
          if ((imA==1).and.(imB==3)) then
             ! if (distA(imA)>distA_C) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(k)=.true.
             ! else
             if ((distA(imA).le.distA_C).and.(.not.bl(nbl)%periodic(k))) then
                ! direction of neighbor is W (1)
                n_dir(k)=1
                ! if faces S or N => swap i-j
                if (k==3.or.k==4) iswap_ij(k)=.true.
                ! no reverse in parallel direction
                !
                ! use distance of A/C with imA to determine
                ! reverse in normal direction
                if ((k==1).or.(k==3)) then
                   if (distA(imA)<distC_imA) irev_normal(k)=.true.
                elseif ((k==2).or.(k==4)) then
                   if (distA(imA)>distC_imA) irev_normal(k)=.true.
                endif
             endif
          endif

          ! Case #2: ex pour imax
          !              B   1----------------2
          !   ^          |   |                |  distance min between A and corners of neighboring domain
          !   | j        |   |      i         |   index imA=3
          !   |          |   |    ----->      |  distance min between B and corners of neighboring domain
          !   ----->     |   |    |           |   index imB=1
          !     i        |   |    | j         |   is_swap=.false.
          !              |   |    v           |   is_rev(k,1)=.true.  for parallel direction
          !              |   |                |   is_rev(k,2)=.false. for perpendicular direction
          !  C-----------A   3----------------4                       [dist(C,imA)>dist(A,imA)]

          if ((imA==3).and.(imB==1)) then
             ! if (distA(imA)>distA_C) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(k)=.true.
             ! else
             if ((distA(imA).le.distA_C).and.(.not.bl(nbl)%periodic(k))) then
                ! direction of neighbor is W (1)
                n_dir(k)=1
                ! if faces S or N => swap i-j
                if (k==3.or.k==4) iswap_ij(k)=.true.
                ! reverse in parallel direction
                irev_parallel(k)=.true.
                ! use distance of A/C with imA to determine
                ! reverse in normal direction
                if ((k==1).or.(k==3)) then
                   if (distA(imA)<distC_imA) irev_normal(k)=.true.
                elseif ((k==2).or.(k==4)) then
                   if (distA(imA)>distC_imA) irev_normal(k)=.true.
                endif
             endif
          endif

          ! Case #3: ex pour imax
          !              B   2----------------4
          !   ^          |   |                |  distance min between A and corners of neighboring domain
          !   | j        |   |    ^           |   index imA=1
          !   |          |   |    | i         |  distance min between B and corners of neighboring domain
          !   ----->     |   |    |           |   index imB=2
          !     i        |   |    ----->      |   is_swap=.true.
          !              |   |     j          |   is_rev(k,1)=.false. for parallel direction
          !              |   |                |   is_rev(k,2)=.false. for perpendicular direction
          !  C-----------A   1----------------3                       [dist(C,imA)>dist(A,imA)]

          if ((imA==1).and.(imB==2)) then
             ! if (distA(imA)>distA_C) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(k)=.true.
             ! else
             if ((distA(imA).le.distA_C).and.(.not.bl(nbl)%periodic(k))) then
                ! direction of neighbor is S (3)
                n_dir(k)=3
                ! if faces E or W => swap i-j
                if (k==1.or.k==2) iswap_ij(k)=.true.
                ! no reverse in parallel direction
                !
                ! use distance of A/C with imA to determine
                ! reverse in normal direction
                if ((k==1).or.(k==3)) then
                   if (distA(imA)<distC_imA) irev_normal(k)=.true.
                elseif ((k==2).or.(k==4)) then
                   if (distA(imA)>distC_imA) irev_normal(k)=.true.
                endif
             endif
          endif

          ! Case #4: ex pour imax
          !              B   1----------------3
          !   ^          |   |                |  distance min between A and corners of neighboring domain
          !   | j        |   |      j         |   index imA=2
          !   |          |   |    ----->      |  distance min between B and corners of neighboring domain
          !   ----->     |   |    |           |   index imB=1
          !     i        |   |    | i         |   is_swap=.true.
          !              |   |    v           |   is_rev(k,1)=.true.  for parallel direction
          !              |   |                |   is_rev(k,2)=.false. for perpendicular direction
          !  C-----------A   2----------------4                       [dist(C,imA)>dist(A,imA)]

          if ((imA==2).and.(imB==1)) then
             ! if (distA(imA)>distA_C) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(k)=.true.
             ! else
             if ((distA(imA).le.distA_C).and.(.not.bl(nbl)%periodic(k))) then
                ! direction of neighbor is S (3)
                n_dir(k)=3
                ! if faces E or W => swap i-j
                if (k==1.or.k==2) iswap_ij(k)=.true.
                ! reverse in parallel direction
                irev_parallel(k)=.true.
                ! use distance of A/C with imA to determine
                ! reverse in normal direction
                if ((k==1).or.(k==3)) then
                   if (distA(imA)<distC_imA) irev_normal(k)=.true.
                elseif ((k==2).or.(k==4)) then
                   if (distA(imA)>distC_imA) irev_normal(k)=.true.
                endif
                !if (dist(1)<dist(2)) irev_normal(k)=.true.
             endif
          endif

          ! Case #5: ex pour imax
          !              B   4----------------3
          !   ^          |   |                |  distance min between A and corners of neighboring domain
          !   | j        |   |         ^      |   index imA=2
          !   |          |   |       j |      |  distance min between B and corners of neighboring domain
          !   ----->     |   |         |      |   index imB=4
          !     i        |   |    <-----      |   is_swap=.false.
          !              |   |     i          |   is_rev(k,1)=.false. for parallel direction
          !              |   |                |   is_rev(k,2)=.true.  for perpendicular direction
          !  C-----------A   2----------------1                       [dist(C,imA)<dist(A,imA)]

          if ((imA==2).and.(imB==4)) then
             ! if (distA(imA)>distA_C) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(k)=.true.
             ! else
             if ((distA(imA).le.distA_C).and.(.not.bl(nbl)%periodic(k))) then
                ! direction of neighbor is E (2)
                n_dir(k)=2
                ! if faces S or N => swap i-j
                if (k==3.or.k==4) iswap_ij(k)=.true.
                ! no reverse in parallel direction
                !
                ! use distance of A/C with imA to determine
                ! reverse in normal direction
                if ((k==1).or.(k==3)) then
                   if (distA(imA)>distC_imA) irev_normal(k)=.true.
                elseif ((k==2).or.(k==4)) then
                   if (distA(imA)<distC_imA) irev_normal(k)=.true.
                endif
             endif
          endif

          ! Case #6: ex pour imax
          !              B   2----------------1
          !   ^          |   |                |  distance min between A and corners of neighboring domain
          !   | j        |   |      i         |   index imA=4
          !   |          |   |    <-----      |  distance min between B and corners of neighboring domain
          !   ----->     |   |         |      |   index imB=2
          !     i        |   |       j |      |   is_swap=.false.
          !              |   |         v      |   is_rev(k,1)=.true.  for parallel direction
          !              |   |                |   is_rev(k,2)=.true.  for perpendicular direction
          !  C-----------A   4----------------3                       [dist(C,imA)<dist(A,imA)]

          if ((imA==4).and.(imB==2)) then
             ! if (distA(imA)>distA_C) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(k)=.true.
             ! else
             if ((distA(imA).le.distA_C).and.(.not.bl(nbl)%periodic(k))) then
                ! direction of neighbor is E (2)
                n_dir(k)=2
                ! if faces S or N => swap i-j
                if (k==3.or.k==4) iswap_ij(k)=.true.
                ! reverse in parallel direction
                irev_parallel(k)=.true.
                ! use distance of A/C with imA to determine
                ! reverse in normal direction
                if ((k==1).or.(k==3)) then
                   if (distA(imA)>distC_imA) irev_normal(k)=.true.
                elseif ((k==2).or.(k==4)) then
                   if (distA(imA)<distC_imA) irev_normal(k)=.true.
                endif
             endif
          endif

          ! Case #7: ex pour imax
          !              B   4----------------2
          !   ^          |   |                |  distance min between A and corners of neighboring domain
          !   | j        |   |         ^      |   index imA=3
          !   |          |   |       i |      |  distance min between B and corners of neighboring domain
          !   ----->     |   |         |      |   index imB=4
          !     i        |   |    <-----      |   is_swap=.true.
          !              |   |     j          |   is_rev(k,1)=.false. for parallel direction
          !              |   |                |   is_rev(k,2)=.true.  for perpendicular direction
          !  C-----------A   3----------------1                       [dist(C,imA)<dist(A,imA)]

          if ((imA==3).and.(imB==4)) then
             ! if (distA(imA)>distA_C) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(k)=.true.
             ! else
             if ((distA(imA).le.distA_C).and.(.not.bl(nbl)%periodic(k))) then
                ! direction of neighbor is N (4)
                n_dir(k)=4
                ! if faces E or W => swap i-j
                if (k==1.or.k==2) iswap_ij(k)=.true.
                ! no reverse in parallel direction
                !
                ! use distance of A/C with imA to determine
                ! reverse in normal direction
                if ((k==1).or.(k==3)) then
                   if (distA(imA)>distC_imA) irev_normal(k)=.true.
                elseif ((k==2).or.(k==4)) then
                   if (distA(imA)<distC_imA) irev_normal(k)=.true.
                endif
             endif
          endif

          ! Case #8: ex pour imax
          !              B   3----------------1
          !   ^          |   |                |  distance min between A and corners of neighboring domain
          !   | j        |   |      j         |   index imA=4
          !   |          |   |    <-----      |  distance min between B and corners of neighboring domain
          !   ----->     |   |         |      |   index imB=3
          !     i        |   |       i |      |   is_swap=.true.
          !              |   |         v      |   is_rev(k,1)=.true.  for parallel direction
          !              |   |                |   is_rev(k,2)=.true.  for perpendicular direction
          !  C-----------A   4----------------2                       [dist(C,imA)<dist(A,imA)]

          if ((imA==4).and.(imB==3)) then
             ! if (distA(imA)>distA_C) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(k)=.true.
             ! else
             if ((distA(imA).le.distA_C).and.(.not.bl(nbl)%periodic(k))) then
                ! direction of neighbor is N (4)
                n_dir(k)=4
                ! if faces E or W => swap i-j
                if (k==1.or.k==2) iswap_ij(k)=.true.
                ! reverse in parallel direction
                irev_parallel(k)=.true.
                ! use distance of A/C with imA to determine
                ! reverse in normal direction
                if ((k==1).or.(k==3)) then
                   if (distA(imA)>distC_imA) irev_normal(k)=.true.
                elseif ((k==2).or.(k==4)) then
                   if (distA(imA)<distC_imA) irev_normal(k)=.true.
                endif
             endif
          endif
       endif
    enddo

  end subroutine check_directions

end module mod_grid_directions
