!===============================================================================
module mod_grid_directions_3d
!===============================================================================
  !> Module to search swap and reverse directions between domains or blocks
!===============================================================================
  !use mod_grid
  !use mod_block
  !use mod_mpi
  use mod_grid_directions
  implicit none
  ! ----------------------------------------------------------------------------
  integer, private :: nbl ! current block
  ! ----------------------------------------------------------------------------
!!$  ! Indicator at domain level:
!!$  ! ==========================
!!$  ! indicators for swap i-j directions
!!$  logical, dimension(6) :: is_swapij,is_swapij2
!!$  ! indicators for swap i-k directions
!!$  logical, dimension(6) :: is_swapik,is_swapik2
!!$  ! indicators for swap j-k directions
!!$  logical, dimension(6) :: is_swapjk,is_swapjk2
!!$  ! indicators for swap i-j-k directions
!!$  logical, dimension(6) :: is_swapijk,is_swapijk2
!!$  ! indicators for reverse directions
!!$  logical, dimension(6,3) :: is_rev,is_rev2
!!$  ! directions of neighboring faces
!!$  integer, dimension(6) :: ndir
!!$  ! ----------------------------------------------------------------------------
!!$  ! Indicator at block level:
!!$  ! =========================
!!$  ! indicators for swap i-j block directions
!!$  logical, dimension(6) :: is_swapij_bl,is_swapij2_bl
!!$  ! indicators for swap i-k block directions
!!$  logical, dimension(6) :: is_swapik_bl,is_swapik2_bl
!!$  ! indicators for swap j-k block directions
!!$  logical, dimension(6) :: is_swapjk_bl,is_swapjk2_bl
!!$  ! indicators for swap i-j-k block directions
!!$  logical, dimension(6) :: is_swapijk_bl,is_swapijk2_bl
!!$  ! indicators for reverse block directions
!!$  logical, dimension(6,3) :: is_rev_bl,is_rev2_bl
!!$  ! directions of neighboring faces for blocks
!!$  integer, dimension(6) :: ndir_bl
  ! ----------------------------------------------------------------------------

contains
  
  !===============================================================
  subroutine grid_directions_3d(neighbors)
  !===============================================================
    !> Determine ij-swap or reverse directions between blocks
  !===============================================================
    use mod_constant
    use warnstop
    implicit none
    ! ------------------------------------------------------------
    integer, dimension(6), intent(in) :: neighbors    
    ! ------------------------------------------------------------
    integer :: n,m,nv,ip
    integer :: n_nb                  ! number of tested neighbors
    integer, dimension(6) :: list_nb ! list of tested neighbors
    integer, dimension(:,:), allocatable :: iproc_corner
    ! ------------------------------------------------------------
    ! /!\ cannot use mod_mpi_part, which includes mod_grid_directions
    ! recreate coord1,coord2,coord3 -> TO BE CHANGED
    integer, dimension(:), allocatable :: coord1,coord2,coord3
    
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
    
    ! Initialize is_swapijk (swap i-j-k directions)
    ! ------------------
    is_swapijk=.false.
    is_swapijk2=.false.
    is_swapijk_bl=.false.
    is_swapijk2_bl=.false.
    
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
    
    ! Current block
    ! =============
    nbl=nob(iproc)
    
    ! Determine block corners
    ! =======================
    !      6---------------8 
    !     /|               /
    !    / |              /|
    !   /  |             / |
    !  3----------------4  |
    !  |   | j ^  k     |  |
    !  |   |   | /      |  |
    !  |   |   |/       |  |
    !  |   |    ---> i  |  |
    !  |   5----------- |--7
    !  |  /             | /
    !  | /              |/
    !  1----------------2
 
    ! if (idepart==1) then ! use global grid

    !    bl(nbl)%xcorner(1)=xgc3(1,1,1)
    !    bl(nbl)%xcorner(2)=xgc3(ngx,1,1)
    !    bl(nbl)%xcorner(3)=xgc3(1,ngy,1)
    !    bl(nbl)%xcorner(4)=xgc3(ngx,ngy,1)
    !    bl(nbl)%xcorner(5)=xgc3(1,1,ngz)
    !    bl(nbl)%xcorner(6)=xgc3(1,ngy,ngz)
    !    bl(nbl)%xcorner(7)=xgc3(ngx,1,ngz)
    !    bl(nbl)%xcorner(8)=xgc3(ngx,ngy,ngz)

    !    bl(nbl)%ycorner(1)=ygc3(1,1,1)
    !    bl(nbl)%ycorner(2)=ygc3(ngx,1,1)
    !    bl(nbl)%ycorner(3)=ygc3(1,ngy,1)
    !    bl(nbl)%ycorner(4)=ygc3(ngx,ngy,1)
    !    bl(nbl)%ycorner(5)=ygc3(1,1,ngz)
    !    bl(nbl)%ycorner(6)=ygc3(1,ngy,ngz)
    !    bl(nbl)%ycorner(7)=ygc3(ngx,1,ngz)
    !    bl(nbl)%ycorner(8)=ygc3(ngx,ngy,ngz)

    !    bl(nbl)%zcorner(1)=zgc3(1,1,1)
    !    bl(nbl)%zcorner(2)=zgc3(ngx,1,1)
    !    bl(nbl)%zcorner(3)=zgc3(1,ngy,1)
    !    bl(nbl)%zcorner(4)=zgc3(ngx,ngy,1)
    !    bl(nbl)%zcorner(5)=zgc3(1,1,ngz)
    !    bl(nbl)%zcorner(6)=zgc3(1,ngy,ngz)
    !    bl(nbl)%zcorner(7)=zgc3(ngx,1,ngz)
    !    bl(nbl)%zcorner(8)=zgc3(ngx,ngy,ngz)

    !    ! share block corner coordinate between procs
    !    do n=1,nbloc
    !       call MPI_BCAST(bl(n)%xcorner,8,MPI_DOUBLE_PRECISION,iproc_leader(n),COMM_global,info)
    !       call MPI_BCAST(bl(n)%ycorner,8,MPI_DOUBLE_PRECISION,iproc_leader(n),COMM_global,info)
    !       call MPI_BCAST(bl(n)%zcorner,8,MPI_DOUBLE_PRECISION,iproc_leader(n),COMM_global,info)
    !    enddo

    ! else ! use local grid

       ! gather coord in the MPI_CART communicator
       allocate(coord1(0:nproc-1),coord2(0:nproc-1),coord3(0:nproc-1))
       call MPI_ALLGATHER(coord(1),1,MPI_INTEGER,coord1,1,MPI_INTEGER,COMM_global,info)
       call MPI_ALLGATHER(coord(2),1,MPI_INTEGER,coord2,1,MPI_INTEGER,COMM_global,info)
       call MPI_ALLGATHER(coord(3),1,MPI_INTEGER,coord3,1,MPI_INTEGER,COMM_global,info)

       ! determine who has corners
       allocate(iproc_corner(nbloc,8))

       do ip=0,nproc-1
          n=nob(ip)
          if ((coord1(ip)==0).and.(coord2(ip)==0).and.(coord3(ip)==0)) &
               iproc_corner(n,1)=ip
          if ((coord1(ip)==bl(n)%ndomi-1).and.(coord2(ip)==0).and.(coord3(ip)==0)) &
               iproc_corner(n,2)=ip
          if ((coord1(ip)==0).and.(coord2(ip)==bl(n)%ndomj-1).and.(coord3(ip)==0)) &
               iproc_corner(n,3)=ip
          if ((coord1(ip)==bl(n)%ndomi-1).and.(coord2(ip)==bl(n)%ndomj-1).and.(coord3(ip)==0)) &
               iproc_corner(n,4)=ip
          if ((coord1(ip)==0).and.(coord2(ip)==0).and.(coord3(ip)==bl(n)%ndomk-1)) &
               iproc_corner(n,5)=ip
          if ((coord1(ip)==0).and.(coord2(ip)==bl(n)%ndomj-1).and.(coord3(ip)==bl(n)%ndomk-1)) &
               iproc_corner(n,6)=ip
          if ((coord1(ip)==bl(n)%ndomi-1).and.(coord2(ip)==0).and.(coord3(ip)==bl(n)%ndomk-1)) &
               iproc_corner(n,7)=ip
          if ((coord1(ip)==bl(n)%ndomi-1).and.(coord2(ip)==bl(n)%ndomj-1).and.(coord3(ip)==bl(n)%ndomk-1)) &
               iproc_corner(n,8)=ip
       enddo

       if ((coord(1)==0).and.(coord(2)==0).and.(coord(3)==0)) then
          bl(nbl)%xcorner(1)=xc3(1,1,1)
          bl(nbl)%ycorner(1)=yc3(1,1,1)
          bl(nbl)%zcorner(1)=zc3(1,1,1)
       endif
       if ((coord(1)==ndomx-1).and.(coord(2)==0).and.(coord(3)==0)) then
          bl(nbl)%xcorner(2)=xc3(nx,1,1)
          bl(nbl)%ycorner(2)=yc3(nx,1,1)
          bl(nbl)%zcorner(2)=zc3(nx,1,1)
       endif
       if ((coord(1)==0).and.(coord(2)==ndomy-1).and.(coord(3)==0)) then
          bl(nbl)%xcorner(3)=xc3(1,ny,1)
          bl(nbl)%ycorner(3)=yc3(1,ny,1)
          bl(nbl)%zcorner(3)=zc3(1,ny,1)
       endif
       if ((coord(1)==ndomx-1).and.(coord(2)==ndomy-1).and.(coord(3)==0)) then
          bl(nbl)%xcorner(4)=xc3(nx,ny,1)
          bl(nbl)%ycorner(4)=yc3(nx,ny,1)
          bl(nbl)%zcorner(4)=zc3(nx,ny,1)
       endif
       if ((coord(1)==0).and.(coord(2)==0).and.(coord(3)==ndomz-1)) then
          bl(nbl)%xcorner(5)=xc3(1,1,nz)
          bl(nbl)%ycorner(5)=yc3(1,1,nz)
          bl(nbl)%zcorner(5)=zc3(1,1,nz)
       endif
       if ((coord(1)==0).and.(coord(2)==ndomy-1).and.(coord(3)==ndomz-1)) then
          bl(nbl)%xcorner(6)=xc3(1,ny,nz)
          bl(nbl)%ycorner(6)=yc3(1,ny,nz)
          bl(nbl)%zcorner(6)=zc3(1,ny,nz)
       endif
       if ((coord(1)==ndomx-1).and.(coord(2)==0).and.(coord(3)==ndomz-1)) then
          bl(nbl)%xcorner(7)=xc3(nx,1,nz)
          bl(nbl)%ycorner(7)=yc3(nx,1,nz)
          bl(nbl)%zcorner(7)=zc3(nx,1,nz)
       endif
       if ((coord(1)==ndomx-1).and.(coord(2)==ndomy-1).and.(coord(3)==ndomz-1)) then
          bl(nbl)%xcorner(8)=xc3(nx,ny,nz)
          bl(nbl)%ycorner(8)=yc3(nx,ny,nz)
          bl(nbl)%zcorner(8)=zc3(nx,ny,nz)
       endif

       do n=1,nbloc
          do m=1,8
             call MPI_BCAST(bl(n)%xcorner(m),1,MPI_DOUBLE_PRECISION,iproc_corner(n,m),COMM_global,info)
             call MPI_BCAST(bl(n)%ycorner(m),1,MPI_DOUBLE_PRECISION,iproc_corner(n,m),COMM_global,info)
             call MPI_BCAST(bl(n)%zcorner(m),1,MPI_DOUBLE_PRECISION,iproc_corner(n,m),COMM_global,info)
          enddo
       enddo

       !if (iproc==7) then
       !   do n=1,nbloc
       !      do m=1,8
       !         print *,n,m,iproc_corner(n,m)
       !         print *,bl(n)%xcorner(m),bl(n)%ycorner(m),bl(n)%zcorner(m)
       !      enddo
       !   enddo
       !endif

    ! endif

    call MPI_BARRIER(COMM_global,info)

    ! Search for swap and reverse directions for domains with lower rank
    ! ==================================================================
    
    ! Identify neighboring domains
    ! ----------------------------
    list_nb=0
    n_nb=0
    do n=1,6
       nv=bl(nbl)%BC(n)
       ! nv>0: connectivity with a domain, nv<n: lower rank,
       ! neighbor(n)==MPI_PROC_NULL means that the domain has no neighbor for the tested face
       ! because it is not yet determined => it belongs to another block (not intrablock)
       ! nv.ne.n: exclude comm with itself,
       if ((nv>0).and.(nv<nbl).and.(neighbors(n)==MPI_PROC_NULL).and.(nv.ne.nbl)) then
          list_nb(n)=bl(nbl)%BC(n)
          n_nb=n_nb+1
       endif
    enddo
    !if (iproc==6) print 100,iproc,list_nb(1:6),n_nb
100 format(1x,'iproc',i4,', list :',6i4,', n_nb ',i4)

    ! Determine indicators
    ! --------------------
    ! NOTA: ndir not used
    if (n_nb>0) call check_directions_3d(list_nb,is_swapij,is_swapik,is_swapjk, &
                                         is_swapijk,is_rev(:,2),is_rev(:,1),is_rev(:,3),ndir)

    ! Search for swap and reverse directions for all domains
    ! ======================================================

    ! Identify neighboring domains
    ! ----------------------------
    list_nb=0
    n_nb=0
    do n=1,6
       nv=bl(nbl)%BC(n)
       ! nv>0: connectivity with a domain,
       ! neighbor(n)==MPI_PROC_NULL means that the domain has no neighbor for the tested face
       ! because it is not yet determined => it belongs to another block (not intrablock),
       ! nv.ne.n: exclude comm with itself
       if ((nv>0).and.(neighbors(n)==MPI_PROC_NULL).and.(nv.ne.nbl)) then
          list_nb(n)=bl(nbl)%BC(n)
          n_nb=n_nb+1
       endif
    enddo   
    !if (iproc==2) print 101,iproc,list_nb(1:6),n_nb
101 format(1x,'iproc',i4,', list2:',6i4,', n_nb ',i4)

    ! Determine indicators
    ! --------------------
    ! NOTA: ndir2 is ndir
    if (n_nb>0) call check_directions_3d(list_nb,is_swapij2,is_swapik2,is_swapjk2, &
                                         is_swapijk2,is_rev2(:,2),is_rev2(:,1),is_rev2(:,3),ndir)
    
!!$    if (iproc==2) then
!!$    print *,iproc,'is_swapij ',is_swapij,'is_swapij2 ',is_swapij2 
!!$    print *,iproc,'is_swapik ',is_swapik,'is_swapik2 ',is_swapik2 
!!$    print *,iproc,'is_swapjk ',is_swapjk,'is_swapjk2 ',is_swapjk2 
!!$    print *,iproc,'is_swapijk',is_swapijk,'is_swapijk2',is_swapijk2 
!!$    print *,iproc,' is_rev_n ',is_rev(:,2),'is_rev2_n ',is_rev2(:,2)
!!$    print *,iproc,' is_rev_p1',is_rev(:,1),'is_rev2_p1',is_rev2(:,1)
!!$    print *,iproc,' is_rev_p2',is_rev(:,3),'is_rev2_p2',is_rev2(:,3)
!!$    endif
!!$    call mpistop('pause check in grid_directions_3d',0)
    
    ! Search for swap and reverse directions for blocks with lower rank
    ! =================================================================

    ! Identify neighboring domains
    ! ----------------------------
    list_nb=0
    n_nb=0
    do n=1,6
       nv=bl(nbl)%BC(n)
       ! nv>0: connectivity with a block, nv<n: lower rank,
       ! nv.ne.n: exclude comm with itself
       if ((nv>0).and.(nv<nbl).and.(nv.ne.nbl)) then
          list_nb(n)=bl(nbl)%BC(n)
          n_nb=n_nb+1
       endif
    enddo    
    !if ((iproc==3).or.(iproc==3)) print 102,iproc,list_nb(1:6),n_nb
102 format(1x,'iproc',i4,', list_bl :',6i4,', n_nb ',i4)

    ! Determine indicators
    ! --------------------
    ! NOTA: ndir_bl not used; is_rev_bl(:,2) not used
    if (n_nb>0) call check_directions_3d(list_nb,is_swapij_bl,is_swapik_bl,is_swapjk_bl, &
                                         is_swapijk_bl,is_rev_bl(:,2),is_rev_bl(:,1),is_rev_bl(:,3),ndir_bl)

    ! Search for swap and reverse directions for all blocks
    ! =====================================================

    ! Identify neighboring domains
    ! ----------------------------
    list_nb=0
    n_nb=0
    do n=1,6
       nv=bl(nbl)%BC(n)
       ! nv>0: connectivity with a block,
       ! nv.ne.n: exclude comm with itself
       if ((nv>0).and.(nv.ne.nbl)) then
          list_nb(n)=bl(nbl)%BC(n)
          n_nb=n_nb+1
       endif
    enddo  
    !if ((iproc==1).or.(iproc==1)) print 103,iproc,list_nb(1:6),n_nb
103 format(1x,'iproc',i4,', list2_bl:',6i4,', n_nb ',i4)
      
    ! Determine indicators
    ! --------------------
    ! NOTA: ndir2_bl is ndir_bl
    if (n_nb>0) call check_directions_3d(list_nb,is_swapij2_bl,is_swapik2_bl,is_swapjk2_bl, &
                                         is_swapijk2_bl,is_rev2_bl(:,2),is_rev2_bl(:,1),is_rev2_bl(:,3),ndir_bl)

!!$    if (iproc==6) then
!!$    print *,iproc,'is_swapij ',is_swapij,'is_swapij2 ',is_swapij2 
!!$    print *,iproc,'is_swapik ',is_swapik,'is_swapik2 ',is_swapik2 
!!$    print *,iproc,'is_swapjk ',is_swapjk,'is_swapjk2 ',is_swapjk2 
!!$    print *,iproc,'is_swapijk',is_swapijk,'is_swapijk2',is_swapijk2 
!!$    print *,iproc,' is_rev_n ',is_rev(:,2),'is_rev2_n ',is_rev2(:,2)
!!$    print *,iproc,' is_rev_p1',is_rev(:,1),'is_rev2_p1',is_rev2(:,1)
!!$    print *,iproc,' is_rev_p2',is_rev(:,3),'is_rev2_p2',is_rev2(:,3)
!!$    endif
!!$    call mpistop('pause check in grid_directions_3d',0)
    
!!$    if (iproc==5) then
!!$    print *,iproc,'is_swapij ',is_swapij,'is_swapij2 ',is_swapij2 
!!$    !print *,iproc,'is_swapik ',is_swapik,'is_swapik2 ',is_swapik2 
!!$    !print *,iproc,'is_swapjk ',is_swapjk,'is_swapjk2 ',is_swapjk2 
!!$    !print *,iproc,'is_swapijk',is_swapijk,'is_swapijk2',is_swapijk2 
!!$    print *,iproc,' is_rev_n ',is_rev(:,2),'is_rev2_n ',is_rev2(:,2)
!!$    print *,iproc,' is_rev_p1',is_rev(:,1),'is_rev2_p1',is_rev2(:,1)
!!$    print *,iproc,' is_rev_p2',is_rev(:,3),'is_rev2_p2',is_rev2(:,3)
!!$    endif
!!$    call mpistop('pause check in grid_directions_3d',0)
        
  end subroutine grid_directions_3d
  
  !=======================================================================
  subroutine check_directions_3d(list_nb,iswap_ij,iswap_ik,iswap_jk,iswap_ijk, &
                                 irev_normal,irev_parallel1,irev_parallel2,n_dir)
  !=======================================================================
    !> Subroutine to determine indicators swap and reverse
  !=======================================================================
    use mod_constant
    implicit none
    ! --------------------------------------------------------------------
    ! INPUT: list of tested neighbors
    integer, dimension(6) :: list_nb
    ! OUTPUT: indicator for swap between i and j directions
    logical, dimension(6) :: iswap_ij,iswap_ik,iswap_jk,iswap_ijk
    ! OUTPUT: indicator for reverse direction (parallel and normal to face)
    logical, dimension(6) :: irev_normal,irev_parallel1,irev_parallel2
    ! OUTPUT: array of neighboring directions
    integer, dimension(6) :: n_dir
    ! --------------------------------------------------------------------
    integer :: n  ! tested direction [1:imin;2:imax;3:jmin;4:jmax;5:kmin;6:kmax]
    integer :: nv ! number of neighboring block
    integer :: nb ! number of current block
    integer :: iA,jA,kA,iB,jB,kB ! point coord
    integer :: iC,jC,kC,iD,jD,kD ! point coord
    integer :: imA,imB,imC  ! point indices
    real(wp) :: xA,yA,zA,xB,yB,zB,xC,yC,zC,xD,yD,zD ! points
    real(wp), dimension(8) :: distA,distB,distC ! distances from corners
    real(wp) :: distD_imA,distA_D               ! distances in normal direction
    ! --------------------------------------------------------------------
    
    ! Current block
    ! =============
    nb=nob(iproc)
    
    ! Check directions
    ! ================
    do n=1,6 ! for each BC
    !do n=4,4 ! for each BC

       if (list_nb(n)>0) then

          nv=list_nb(n)

!           if (idepart==1) then ! use global grid

!              select case(n)

!              case (1)                ! Face imin (1)
!                 iA=1;  jA=1  ;kA=1   ! =============
!                                      !           B     (iA,jA,kA)=(1,1,1)
!                 iB=1;  jB=ngy;kB=1   !     j  k  |
!                                      !     | /   |   C (iB,jB,kB)=(1,ngy,1)
!                 iC=1;  jC=1;  kC=ngz ! i___|/    |  /
!                                      !           | /   (iC,jC,kC)=(1,1,ngz)
!                 iD=ngx;jD=1;  kD=1   !           |/
!                                      !  D--------A     (iD,jD,kD)=(ngx,1,1)

!              case (2)                ! Face imax (2)
!                 iA=ngx;jA=1  ;kA=1   ! =============
!                                      !           B     (iA,jA,kA)=(ngx,1,1)
!                 iB=ngx;jB=ngy;kB=1   !   j  k    |
!                                      !   | /     |   C (iB,jB,kB)=(ngx,ngy,1)
!                 iC=ngx;jC=1;  kC=ngz !   |/___i  |  /
!                                      !           | /   (iC,jC,kC)=(ngx,1,ngz)
!                 iD=1  ;jD=1;  kD=1   !           |/
!                                      !  D--------A     (iD,jD,kD)=(1,1,1)

!              case (3)                ! Face jmin (3)
!                 iA=1;  jA=1  ;kA=1   ! =============
!                                      !           B     (iA,jA,kA)=(1,1,1)
!                 iB=ngx;jB=1;  kB=1   !     i  k  |
!                                      !     | /   |   C (iB,jB,kB)=(ngx,1,1)
!                 iC=1;  jC=1;  kC=ngz ! j___|/    |  /
!                                      !           | /   (iC,jC,kC)=(1,1,ngz)
!                 iD=1  ;jD=ngy;kD=1   !           |/
!                                      !  D--------A     (iD,jD,kD)=(1,ngy,1)

!              case (4)                ! Face jmax (4)
!                 iA=1  ;jA=ngy;kA=1   ! =============
!                                      !           B     (iA,jA,kA)=(1,ngy,1)
!                 iB=ngx;jB=ngy;kB=1   !   i  k    |
!                                      !   | /     |   C (iB,jB,kB)=(ngx,ngy,1)
!                 iC=1  ;jC=ngy;kC=ngz !   |/___j  |  /
!                                      !           | /   (iC,jC,kC)=(1,1ngy,ngz)
!                 iD=1  ;jD=1;  kD=1   !           |/
!                                      !  D--------A     (iD,jD,kD)=(1,1,1)

!              case (5)                ! Face kmin (5)
!                 iA=1;  jA=1  ;kA=1   ! =============
!                                      !           B     (iA,jA,kA)=(1,1,1)
!                 iB=ngx;jB=1;  kB=1   !     i  j  |
!                                      !     | /   |   C (iB,jB,kB)=(1,ngy,1)
!                 iC=1;  jC=ngy;kC=1   ! k___|/    |  /
!                                      !           | /   (iC,jC,kC)=(1,1,ngz)
!                 iD=1  ;jD=1;  kD=ngz !           |/
!                                      !  D--------A     (iD,jD,kD)=(ngx,1,1)

!              case (6)                ! Face kmax (6)
!                 iA=1  ;jA=1  ;kA=ngz ! =============
!                                      !           B     (iA,jA,kA)=(ngx,1,1)
!                 iB=ngx;jB=1;  kB=ngz !   i  j    |
!                                      !   | /     |   C (iB,jB,kB)=(ngx,ngy,1)
!                 iC=1;  jC=ngy;kC=ngz !   |/___k  |  /
!                                      !           | /   (iC,jC,kC)=(ngx,1,ngz)
!                 iD=1  ;jD=1;  kD=1   !           |/
!                                      !  D--------A     (iD,jD,kD)=(1,1,1)
!              end select

!              ! distance A-corner1
!              distA(1)=(xgc3(iA,jA,kA)-bl(nv)%xcorner(1))**2 &
!                      +(ygc3(iA,jA,kA)-bl(nv)%ycorner(1))**2 &
!                      +(zgc3(iA,jA,kA)-bl(nv)%zcorner(1))**2
!              ! distance A-corner2
!              distA(2)=(xgc3(iA,jA,kA)-bl(nv)%xcorner(2))**2 &
!                      +(ygc3(iA,jA,kA)-bl(nv)%ycorner(2))**2 &
!                      +(zgc3(iA,jA,kA)-bl(nv)%zcorner(2))**2
!              ! distance A-corner3
!              distA(3)=(xgc3(iA,jA,kA)-bl(nv)%xcorner(3))**2 &
!                      +(ygc3(iA,jA,kA)-bl(nv)%ycorner(3))**2 &
!                      +(zgc3(iA,jA,kA)-bl(nv)%zcorner(3))**2
!              ! distance A-corner4
!              distA(4)=(xgc3(iA,jA,kA)-bl(nv)%xcorner(4))**2 &
!                      +(ygc3(iA,jA,kA)-bl(nv)%ycorner(4))**2 &
!                      +(zgc3(iA,jA,kA)-bl(nv)%zcorner(4))**2
!              ! distance A-corner5
!              distA(5)=(xgc3(iA,jA,kA)-bl(nv)%xcorner(5))**2 &
!                      +(ygc3(iA,jA,kA)-bl(nv)%ycorner(5))**2 &
!                      +(zgc3(iA,jA,kA)-bl(nv)%zcorner(5))**2
!              ! distance A-corner6
!              distA(6)=(xgc3(iA,jA,kA)-bl(nv)%xcorner(6))**2 &
!                      +(ygc3(iA,jA,kA)-bl(nv)%ycorner(6))**2 &
!                      +(zgc3(iA,jA,kA)-bl(nv)%zcorner(6))**2
!              ! distance A-corner7
!              distA(7)=(xgc3(iA,jA,kA)-bl(nv)%xcorner(7))**2 &
!                      +(ygc3(iA,jA,kA)-bl(nv)%ycorner(7))**2 &
!                      +(zgc3(iA,jA,kA)-bl(nv)%zcorner(7))**2
!              ! distance A-corner8
!              distA(8)=(xgc3(iA,jA,kA)-bl(nv)%xcorner(8))**2 &
!                      +(ygc3(iA,jA,kA)-bl(nv)%ycorner(8))**2 &
!                      +(zgc3(iA,jA,kA)-bl(nv)%zcorner(8))**2
!              ! closest corner to A
!              imA=minloc(distA,1)

!              ! distance B-corner1
!              distB(1)=(xgc3(iB,jB,kB)-bl(nv)%xcorner(1))**2 &
!                      +(ygc3(iB,jB,kB)-bl(nv)%ycorner(1))**2 &
!                      +(zgc3(iB,jB,kB)-bl(nv)%zcorner(1))**2
!              ! distance B-corner2
!              distB(2)=(xgc3(iB,jB,kB)-bl(nv)%xcorner(2))**2 &
!                      +(ygc3(iB,jB,kB)-bl(nv)%ycorner(2))**2 &
!                      +(zgc3(iB,jB,kB)-bl(nv)%zcorner(2))**2
!              ! distance B-corner3
!              distB(3)=(xgc3(iB,jB,kB)-bl(nv)%xcorner(3))**2 &
!                      +(ygc3(iB,jB,kB)-bl(nv)%ycorner(3))**2 &
!                      +(zgc3(iB,jB,kB)-bl(nv)%zcorner(3))**2
!              ! distance B-corner4
!              distB(4)=(xgc3(iB,jB,kB)-bl(nv)%xcorner(4))**2 &
!                      +(ygc3(iB,jB,kB)-bl(nv)%ycorner(4))**2 &
!                      +(zgc3(iB,jB,kB)-bl(nv)%zcorner(4))**2
!              ! distance B-corner5
!              distB(5)=(xgc3(iB,jB,kB)-bl(nv)%xcorner(5))**2 &
!                      +(ygc3(iB,jB,kB)-bl(nv)%ycorner(5))**2 &
!                      +(zgc3(iB,jB,kB)-bl(nv)%zcorner(5))**2
!              ! distance B-corner6
!              distB(6)=(xgc3(iB,jB,kB)-bl(nv)%xcorner(6))**2 &
!                      +(ygc3(iB,jB,kB)-bl(nv)%ycorner(6))**2 &
!                      +(zgc3(iB,jB,kB)-bl(nv)%zcorner(6))**2
!              ! distance B-corner7
!              distB(7)=(xgc3(iB,jB,kB)-bl(nv)%xcorner(7))**2 &
!                      +(ygc3(iB,jB,kB)-bl(nv)%ycorner(7))**2 &
!                      +(zgc3(iB,jB,kB)-bl(nv)%zcorner(7))**2
!              ! distance B-corner8
!              distB(8)=(xgc3(iB,jB,kB)-bl(nv)%xcorner(8))**2 &
!                      +(ygc3(iB,jB,kB)-bl(nv)%ycorner(8))**2 &
!                      +(zgc3(iB,jB,kB)-bl(nv)%zcorner(8))**2
!              ! closest corner to B
!              imB=minloc(distB,1)

!              ! distance C-corner1
!              distC(1)=(xgc3(iC,jC,kC)-bl(nv)%xcorner(1))**2 &
!                      +(ygc3(iC,jC,kC)-bl(nv)%ycorner(1))**2 &
!                      +(zgc3(iC,jC,kC)-bl(nv)%zcorner(1))**2
!              ! distance C-corner2
!              distC(2)=(xgc3(iC,jC,kC)-bl(nv)%xcorner(2))**2 &
!                      +(ygc3(iC,jC,kC)-bl(nv)%ycorner(2))**2 &
!                      +(zgc3(iC,jC,kC)-bl(nv)%zcorner(2))**2
!              ! distance C-corner3
!              distC(3)=(xgc3(iC,jC,kC)-bl(nv)%xcorner(3))**2 &
!                      +(ygc3(iC,jC,kC)-bl(nv)%ycorner(3))**2 &
!                      +(zgc3(iC,jC,kC)-bl(nv)%zcorner(3))**2
!              ! distance C-corner4
!              distC(4)=(xgc3(iC,jC,kC)-bl(nv)%xcorner(4))**2 &
!                      +(ygc3(iC,jC,kC)-bl(nv)%ycorner(4))**2 &
!                      +(zgc3(iC,jC,kC)-bl(nv)%zcorner(4))**2
!              ! distance C-corner5
!              distC(5)=(xgc3(iC,jC,kC)-bl(nv)%xcorner(5))**2 &
!                      +(ygc3(iC,jC,kC)-bl(nv)%ycorner(5))**2 &
!                      +(zgc3(iC,jC,kC)-bl(nv)%zcorner(5))**2
!              ! distance C-corner6
!              distC(6)=(xgc3(iC,jC,kC)-bl(nv)%xcorner(6))**2 &
!                      +(ygc3(iC,jC,kC)-bl(nv)%ycorner(6))**2 &
!                      +(zgc3(iC,jC,kC)-bl(nv)%zcorner(6))**2
!              ! distance C-corner7
!              distC(7)=(xgc3(iC,jC,kC)-bl(nv)%xcorner(7))**2 &
!                      +(ygc3(iC,jC,kC)-bl(nv)%ycorner(7))**2 &
!                      +(zgc3(iC,jC,kC)-bl(nv)%zcorner(7))**2
!              ! distance C-corner8
!              distC(8)=(xgc3(iC,jC,kC)-bl(nv)%xcorner(8))**2 &
!                      +(ygc3(iC,jC,kC)-bl(nv)%ycorner(8))**2 &
!                      +(zgc3(iC,jC,kC)-bl(nv)%zcorner(8))**2
!              ! closest corner to C
!              imC=minloc(distC,1)

! !!$             if ((iproc==2).and.(n==4)) then
! !!$                do i=1,8
! !!$                   print *,i,bl(nv)%xcorner(i),bl(nv)%ycorner(i),bl(nv)%zcorner(i)
! !!$                enddo
! !!$             endif

!              ! distance D-imA (closest corner to A)
!              ! [used to check rev_normal]
!              distD_imA=(xgc3(iD,jD,kD)-bl(nv)%xcorner(imA))**2 &
!                       +(ygc3(iD,jD,kD)-bl(nv)%ycorner(imA))**2 &
!                       +(zgc3(iD,jD,kD)-bl(nv)%zcorner(imA))**2

!              ! distance A-D
!              ! [used to check periodicity]
!              distA_D=(xgc3(iD,jD,kD)-xgc3(iA,jA,kA))**2 &
!                     +(ygc3(iD,jD,kD)-ygc3(iA,jA,kA))**2 &
!                     +(zgc3(iD,jD,kD)-zgc3(iA,jA,kA))**2

! !!$             if ((iproc==2).and.(n==4)) then
! !!$                print *,'neighbor',nv
! !!$                print *,'point A',xgc3(iA,jA,kA),ygc3(iA,jA,kA),zgc3(iA,jA,kA)
! !!$                print *,'point B',xgc3(iB,jB,kB),ygc3(iB,jB,kB),zgc3(iB,jB,kB)
! !!$                print *,'point C',xgc3(iC,jC,kC),ygc3(iC,jC,kC),zgc3(iC,jC,kC)
! !!$                print *,nv,distA(1),distA(2),distA(3),distA(4)
! !!$                print *,nv,distA(5),distA(6),distA(7),distA(8)
! !!$                print *,imA,imB,imC,distD_imA,distA_D
! !!$             endif

!           else ! use local grid

             select case(n)

             case (1)                 ! Face imin (1)
                !iA=1;  jA=1  ;kA=1   ! =============
                !                     !           B     (iA,jA,kA)=(1,1,1)  -> corner 1
                !iB=1;  jB=ny;kB=1    !     j  k  |     
                !                     !     | /   |   C (iB,jB,kB)=(1,ny,1) -> corner 3
                !iC=1;  jC=1;  kC=nz  ! i___|/    |  /  
                !                     !           | /   (iC,jC,kC)=(1,1,nz) -> corner 5
                !iD=nx;jD=1;  kD=1    !           |/     
                !                     !  D--------A     (iD,jD,kD)=(nx,1,1) -> corner 2
                xA=bl(nb)%xcorner(1)
                yA=bl(nb)%ycorner(1)
                zA=bl(nb)%zcorner(1)
                xB=bl(nb)%xcorner(3)
                yB=bl(nb)%ycorner(3)
                zB=bl(nb)%zcorner(3)
                xC=bl(nb)%xcorner(5)
                yC=bl(nb)%ycorner(5)
                zC=bl(nb)%zcorner(5)
                xD=bl(nb)%xcorner(2)
                yD=bl(nb)%ycorner(2)
                zD=bl(nb)%zcorner(2)

             case (2)                 ! Face imax (2)
                !iA=nx;jA=1  ;kA=1    ! =============
                !                     !           B     (iA,jA,kA)=(nx,1,1)  -> corner 2
                !iB=nx;jB=ny;kB=1     !   j  k    |     
                !                     !   | /     |   C (iB,jB,kB)=(nx,ny,1) -> corner 4
                !iC=nx;jC=1;  kC=nz   !   |/___i  |  /  
                !                     !           | /   (iC,jC,kC)=(nx,1,nz) -> corner 7
                !iD=1  ;jD=1;  kD=1   !           |/     
                !                     !  D--------A     (iD,jD,kD)=(1,1,1)   -> corner 1
                xA=bl(nb)%xcorner(2)
                yA=bl(nb)%ycorner(2)
                zA=bl(nb)%zcorner(2)
                xB=bl(nb)%xcorner(4)
                yB=bl(nb)%ycorner(4)
                zB=bl(nb)%zcorner(4)
                xC=bl(nb)%xcorner(7)
                yC=bl(nb)%ycorner(7)
                zC=bl(nb)%zcorner(7)
                xD=bl(nb)%xcorner(1)
                yD=bl(nb)%ycorner(1)
                zD=bl(nb)%zcorner(1)

             case (3)                 ! Face jmin (3)
                !iA=1;  jA=1  ;kA=1   ! =============
                !                     !           B     (iA,jA,kA)=(1,1,1)  -> corner 1
                !iB=nx;jB=1;  kB=1    !     i  k  |     
                !                     !     | /   |   C (iB,jB,kB)=(nx,1,1) -> corner 2
                !iC=1;  jC=1;  kC=nz  ! j___|/    |  /  
                !                     !           | /   (iC,jC,kC)=(1,1,nz) -> corner 5
                !iD=1  ;jD=ny;kD=1    !           |/     
                !                     !  D--------A     (iD,jD,kD)=(1,ny,1) -> corner 3
                xA=bl(nb)%xcorner(1)
                yA=bl(nb)%ycorner(1)
                zA=bl(nb)%zcorner(1)
                xB=bl(nb)%xcorner(2)
                yB=bl(nb)%ycorner(2)
                zB=bl(nb)%zcorner(2)
                xC=bl(nb)%xcorner(5)
                yC=bl(nb)%ycorner(5)
                zC=bl(nb)%zcorner(5)
                xD=bl(nb)%xcorner(3)
                yD=bl(nb)%ycorner(3)
                zD=bl(nb)%zcorner(3)

             case (4)                 ! Face jmax (4)
                !iA=1  ;jA=ny;kA=1    ! =============
                !                     !           B     (iA,jA,kA)=(1,ny,1)  -> corner 3
                !iB=nx;jB=ny;kB=1     !   i  k    |     
                !                     !   | /     |   C (iB,jB,kB)=(nx,ny,1) -> corner 4
                !iC=1  ;jC=ny;kC=nz   !   |/___j  |  /  
                !                     !           | /   (iC,jC,kC)=(1,ny,nz) -> corner 6
                !iD=1  ;jD=1;  kD=1   !           |/     
                !                     !  D--------A     (iD,jD,kD)=(1,1,1)   -> corner 1
                xA=bl(nb)%xcorner(3)
                yA=bl(nb)%ycorner(3)
                zA=bl(nb)%zcorner(3)
                xB=bl(nb)%xcorner(4)
                yB=bl(nb)%ycorner(4)
                zB=bl(nb)%zcorner(4)
                xC=bl(nb)%xcorner(6)
                yC=bl(nb)%ycorner(6)
                zC=bl(nb)%zcorner(6)
                xD=bl(nb)%xcorner(1)
                yD=bl(nb)%ycorner(1)
                zD=bl(nb)%zcorner(1)

             case (5)                 ! Face kmin (5)
                !iA=1;  jA=1  ;kA=1   ! =============
                !                     !           B     (iA,jA,kA)=(1,1,1)  -> corner 1
                !iB=nx;jB=1;  kB=1    !     i  j  |     
                !                     !     | /   |   C (iB,jB,kB)=(1,ny,1) -> corner 3
                !iC=1;  jC=ny;kC=1    ! k___|/    |  /  
                !                     !           | /   (iC,jC,kC)=(1,1,nz) -> corner 5
                !iD=1  ;jD=1;  kD=nz  !           |/     
                !                     !  D--------A     (iD,jD,kD)=(nx,1,1) -> corner 2
                xA=bl(nb)%xcorner(1)
                yA=bl(nb)%ycorner(1)
                zA=bl(nb)%zcorner(1)
                xB=bl(nb)%xcorner(3)
                yB=bl(nb)%ycorner(3)
                zB=bl(nb)%zcorner(3)
                xC=bl(nb)%xcorner(5)
                yC=bl(nb)%ycorner(5)
                zC=bl(nb)%zcorner(5)
                xD=bl(nb)%xcorner(2)
                yD=bl(nb)%ycorner(2)
                zD=bl(nb)%zcorner(2)

             case (6)                 ! Face kmax (6)
                !iA=1  ;jA=1  ;kA=nz  ! =============
                !                     !           B     (iA,jA,kA)=(nx,1,1)  -> corner 2
                !iB=nx;jB=1;  kB=nz   !   i  j    |     
                !                     !   | /     |   C (iB,jB,kB)=(nx,ny,1) -> corner 4
                !iC=1;  jC=ny;kC=nz   !   |/___k  |  /  
                !                     !           | /   (iC,jC,kC)=(nx,1,nz) -> corner 7
                !iD=1  ;jD=1;  kD=1   !           |/     
                !                     !  D--------A     (iD,jD,kD)=(1,1,1)   -> corner 1                        
                xA=bl(nb)%xcorner(2)
                yA=bl(nb)%ycorner(2)
                zA=bl(nb)%zcorner(2)
                xB=bl(nb)%xcorner(4)
                yB=bl(nb)%ycorner(4)
                zB=bl(nb)%zcorner(4)
                xC=bl(nb)%xcorner(7)
                yC=bl(nb)%ycorner(7)
                zC=bl(nb)%zcorner(7)
                xD=bl(nb)%xcorner(1)
                yD=bl(nb)%ycorner(1)
                zD=bl(nb)%zcorner(1)
             end select

             ! distance A-corner1
             distA(1)=(xA-bl(nv)%xcorner(1))**2 &
                     +(yA-bl(nv)%ycorner(1))**2 &
                     +(zA-bl(nv)%zcorner(1))**2
             ! distance A-corner2
             distA(2)=(xA-bl(nv)%xcorner(2))**2 &
                     +(yA-bl(nv)%ycorner(2))**2 &
                     +(zA-bl(nv)%zcorner(2))**2
             ! distance A-corner3
             distA(3)=(xA-bl(nv)%xcorner(3))**2 &
                     +(yA-bl(nv)%ycorner(3))**2 &
                     +(zA-bl(nv)%zcorner(3))**2
             ! distance A-corner4
             distA(4)=(xA-bl(nv)%xcorner(4))**2 &
                     +(yA-bl(nv)%ycorner(4))**2 &
                     +(zA-bl(nv)%zcorner(4))**2
             ! distance A-corner5
             distA(5)=(xA-bl(nv)%xcorner(5))**2 &
                     +(yA-bl(nv)%ycorner(5))**2 &
                     +(zA-bl(nv)%zcorner(5))**2
             ! distance A-corner6
             distA(6)=(xA-bl(nv)%xcorner(6))**2 &
                     +(yA-bl(nv)%ycorner(6))**2 &
                     +(zA-bl(nv)%zcorner(6))**2
             ! distance A-corner7
             distA(7)=(xA-bl(nv)%xcorner(7))**2 &
                     +(yA-bl(nv)%ycorner(7))**2 &
                     +(zA-bl(nv)%zcorner(7))**2
             ! distance A-corner8
             distA(8)=(xA-bl(nv)%xcorner(8))**2 &
                     +(yA-bl(nv)%ycorner(8))**2 &
                     +(zA-bl(nv)%zcorner(8))**2
             ! closest corner to A
             imA=minloc(distA,1)

             ! distance B-corner1
             distB(1)=(xB-bl(nv)%xcorner(1))**2 &
                     +(yB-bl(nv)%ycorner(1))**2 &
                     +(zB-bl(nv)%zcorner(1))**2
             ! distance B-corner2
             distB(2)=(xB-bl(nv)%xcorner(2))**2 &
                     +(yB-bl(nv)%ycorner(2))**2 &
                     +(zB-bl(nv)%zcorner(2))**2
             ! distance B-corner3
             distB(3)=(xB-bl(nv)%xcorner(3))**2 &
                     +(yB-bl(nv)%ycorner(3))**2 &
                     +(zB-bl(nv)%zcorner(3))**2
             ! distance B-corner4
             distB(4)=(xB-bl(nv)%xcorner(4))**2 &
                     +(yB-bl(nv)%ycorner(4))**2 &
                     +(zB-bl(nv)%zcorner(4))**2
             ! distance B-corner5
             distB(5)=(xB-bl(nv)%xcorner(5))**2 &
                     +(yB-bl(nv)%ycorner(5))**2 &
                     +(zB-bl(nv)%zcorner(5))**2
             ! distance B-corner6
             distB(6)=(xB-bl(nv)%xcorner(6))**2 &
                     +(yB-bl(nv)%ycorner(6))**2 &
                     +(zB-bl(nv)%zcorner(6))**2
             ! distance B-corner7
             distB(7)=(xB-bl(nv)%xcorner(7))**2 &
                     +(yB-bl(nv)%ycorner(7))**2 &
                     +(zB-bl(nv)%zcorner(7))**2
             ! distance B-corner8
             distB(8)=(xB-bl(nv)%xcorner(8))**2 &
                     +(yB-bl(nv)%ycorner(8))**2 &
                     +(zB-bl(nv)%zcorner(8))**2
             ! closest corner to B
             imB=minloc(distB,1)

             ! distance C-corner1
             distC(1)=(xC-bl(nv)%xcorner(1))**2 &
                     +(yC-bl(nv)%ycorner(1))**2 &
                     +(zC-bl(nv)%zcorner(1))**2
             ! distance C-corner2
             distC(2)=(xC-bl(nv)%xcorner(2))**2 &
                     +(yC-bl(nv)%ycorner(2))**2 &
                     +(zC-bl(nv)%zcorner(2))**2
             ! distance C-corner3
             distC(3)=(xC-bl(nv)%xcorner(3))**2 &
                     +(yC-bl(nv)%ycorner(3))**2 &
                     +(zC-bl(nv)%zcorner(3))**2
             ! distance C-corner4
             distC(4)=(xC-bl(nv)%xcorner(4))**2 &
                     +(yC-bl(nv)%ycorner(4))**2 &
                     +(zC-bl(nv)%zcorner(4))**2
             ! distance C-corner5
             distC(5)=(xC-bl(nv)%xcorner(5))**2 &
                     +(yC-bl(nv)%ycorner(5))**2 &
                     +(zC-bl(nv)%zcorner(5))**2
             ! distance C-corner6
             distC(6)=(xC-bl(nv)%xcorner(6))**2 &
                  +(yC-bl(nv)%ycorner(6))**2 &
                  +(zC-bl(nv)%zcorner(6))**2
             ! distance C-corner7
             distC(7)=(xC-bl(nv)%xcorner(7))**2 &
                     +(yC-bl(nv)%ycorner(7))**2 &
                     +(zC-bl(nv)%zcorner(7))**2
             ! distance C-corner8
             distC(8)=(xC-bl(nv)%xcorner(8))**2 &
                     +(yC-bl(nv)%ycorner(8))**2 &
                     +(zC-bl(nv)%zcorner(8))**2
             ! closest corner to C
             imC=minloc(distC,1)

!!$             if ((iproc==2).and.(n==4)) then
!!$                do i=1,8
!!$                   print *,i,bl(nv)%xcorner(i),bl(nv)%ycorner(i),bl(nv)%zcorner(i)
!!$                enddo
!!$             endif

             ! distance D-imA (closest corner to A)
             ! [used to check rev_normal]
             distD_imA=(xD-bl(nv)%xcorner(imA))**2 &
                      +(yD-bl(nv)%ycorner(imA))**2 &
                      +(zD-bl(nv)%zcorner(imA))**2

             ! distance A-D
             ! [used to check periodicity]
             distA_D=(xD-xA)**2 &
                    +(yD-yA)**2 &
                    +(zD-zA)**2

!!$             if ((iproc==2).and.(n==4)) then
!!$                print *,'neighbor',nv
!!$                print *,'point A',xA,yA,zA
!!$                print *,'point B',xB,yB,zB
!!$                print *,'point C',xC,yC,zC
!!$                print *,nv,distA(1),distA(2),distA(3),distA(4)
!!$                print *,nv,distA(5),distA(6),distA(7),distA(8)
!!$                print *,imA,imB,imC,distD_imA,distA_D
!!$             endif

          ! endif
       
          ! Case #1: ex pour imax (n=2)
          ! ========            6---------------8   distance min between A and corners of neighboring domain
          !                    / |              /|   index imA=1
          !            B      /  |             / |  distance min between B and corners of neighboring domain
          !            |     3----------------4  |   index imB=3
          !            |     |   | j ^  k     |  |  distance min between C and corners of neighboring domain
          !    j  k    |     |   |   | /      |  |   index imC=5
          !    | /     |   C |   |   |/       |  |  is_swapij=.false.
          !    |/___i  |  /  |   |    ---> i  |  |  is_swapik=.false.
          !            | /   |   5----------- |--7  is_swapjk=.false.
          !            |/    |  /             | /   is_rev(k,2)=.false. for normal dir.
          !   D--------A     | /              |/    is_rev(k,1)=.false. for parallel dir. #1
          !                  1----------------2     is_rev(k,3)=.false. for parallel dir. #2

          ! contact with imin face [direct orientation]
          ! -------------------------------------------            
          if ((imA==1).and.(imB==3).and.(imC==5)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 1
                n_dir(n)=1
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_ij(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ijk(n)=.true.
                ! no reverse in parallel direction
                !
                ! reverse in normal direction
                if ((n==1).or.(n==3).or.(n==5)) irev_normal(n)=.true.
             endif
          endif
          
          if ((imA==5).and.(imB==1).and.(imC==6)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 1
                n_dir(n)=1
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_jk(n)=.true.
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_ijk(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ik(n)=.true.
                ! reverse in parallel direction
                irev_parallel1(n)=.true. ! for parallel dir. #1
                ! reverse in normal direction
                if ((n==1).or.(n==3).or.(n==5)) irev_normal(n)=.true.
             endif
          endif
          
          if ((imA==6).and.(imB==5).and.(imC==3)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 1
                n_dir(n)=1
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_ij(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ijk(n)=.true.
                ! no reverse in parallel direction
                irev_parallel1(n)=.true. ! for parallel dir. #1
                irev_parallel2(n)=.true. ! for parallel dir. #2
                ! reverse in normal direction
                if ((n==1).or.(n==3).or.(n==5)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==3).and.(imB==6).and.(imC==1)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 1
                n_dir(n)=1
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_jk(n)=.true.
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_ijk(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ik(n)=.true.
                ! reverse in parallel direction
                irev_parallel2(n)=.true. ! for parallel dir. #2
                ! reverse in normal direction
                if ((n==1).or.(n==3).or.(n==5)) irev_normal(n)=.true.
             endif
          endif
          
          ! Case #2: ex pour imax (n=2)
          ! ========            5---------------7   distance min between A and corners of neighboring domain
          !                    / |              /|   index imA=3
          !            B      /  |             / |  distance min between B and corners of neighboring domain
          !            |     1----------------2  |   index imB=1
          !            |     |   |      k     |  |  distance min between C and corners of neighboring domain
          !    j  k    |     |   |     /      |  |   index imC=5
          !    | /     |   C |   |    /       |  |  is_swapij=.false.
          !    |/___i  |  /  |   |   | ---> i |  |  is_swapik=.false.
          !            | /   |   6---|--------|--8  is_swapjk=.false.
          !            |/    |  /    v j      | /   is_rev(k,2)=.false. for normal dir.
          !   D--------A     | /              |/    is_rev(k,1)=.true.  for parallel dir. #1
          !                  3----------------4     is_rev(k,3)=.false. for parallel dir. #2

          ! contact with imin face [indirect orientation]
          ! ---------------------------------------------
          if ((imA==3).and.(imB==1).and.(imC==6)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 1
                n_dir(n)=1
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_ij(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ijk(n)=.true.
                ! reverse in parallel direction
                irev_parallel1(n)=.true. ! for parallel dir. #1
                ! reverse in normal direction
                if ((n==1).or.(n==3).or.(n==5)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==1).and.(imB==5).and.(imC==3)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 1
                n_dir(n)=1
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_jk(n)=.true.              
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_ijk(n)=.true. 
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ik(n)=.true.
                ! no reverse in parallel direction
                !
                ! reverse in normal direction
                if ((n==1).or.(n==3).or.(n==5)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==5).and.(imB==6).and.(imC==1)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 1
                n_dir(n)=1
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_ij(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ijk(n)=.true.
                ! reverse in parallel direction
                irev_parallel2(n)=.true. ! for parallel dir. #2
                ! reverse in normal direction
                if ((n==1).or.(n==3).or.(n==5)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==6).and.(imB==3).and.(imC==5)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 1
                n_dir(n)=1
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_jk(n)=.true.              
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_ijk(n)=.true. 
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ik(n)=.true.
                ! reverse in parallel directions
                irev_parallel1(n)=.true. ! for parallel dir. #1
                irev_parallel2(n)=.true. ! for parallel dir. #2
                ! reverse in normal direction
                if ((n==1).or.(n==3).or.(n==5)) irev_normal(n)=.true.
             endif
          endif

          ! Case #3: ex pour imax (n=2)
          ! ========            5---------------6   distance min between A and corners of neighboring domain
          !                    / |              /|   index imA=2
          !            B      /  |             / |  distance min between B and corners of neighboring domain
          !            |     1----------------3  |   index imB=1
          !            |     |   |    / k     |  |  distance min between C and corners of neighboring domain
          !    j  k    |     |   |   /___ j   |  |    index imC=7
          !    | /     |   C |   |   |        |  |  is_swapij=.true.
          !    |/___i  |  /  |   |   v        |  |  is_swapik=.false.
          !            | /   |   7-- i -------|--8  is_swapjk=.false.
          !            |/    |  /             | /   is_rev(k,2)=.false. for normal dir.
          !   D--------A     | /              |/    is_rev(k,1)=.true.  for parallel dir. #1
          !                  2----------------4     is_rev(k,3)=.false. for parallel dir. #2

          ! contact with jmin face [direct orientation]
          ! -------------------------------------------
          if ((imA==2).and.(imB==1).and.(imC==7)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 3
                n_dir(n)=3
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_ij(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_jk(n)=.true.
                ! reverse in parallel direction
                irev_parallel1(n)=.true. ! for parallel dir. #1
                ! reverse in normal direction
                if ((n==1).or.(n==3).or.(n==5)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==1).and.(imB==5).and.(imC==2)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 3
                n_dir(n)=3
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_ijk(n)=.true.
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_ik(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ijk(n)=.true.
                ! no reverse in parallel direction
                !
                ! reverse in normal direction
                if ((n==1).or.(n==3).or.(n==5)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==5).and.(imB==7).and.(imC==1)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 3
                n_dir(n)=3
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_ij(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_jk(n)=.true.
                ! reverse in parallel direction
                irev_parallel2(n)=.true. ! for parallel dir. #2
                ! reverse in normal direction
                if ((n==1).or.(n==3).or.(n==5)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==7).and.(imB==2).and.(imC==5)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 3
                n_dir(n)=3
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_ijk(n)=.true.
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_ik(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ijk(n)=.true.
                ! no reverse in parallel direction
                irev_parallel1(n)=.true. ! for parallel dir. #1
                irev_parallel2(n)=.true. ! for parallel dir. #2
                ! reverse in normal direction
                if ((n==1).or.(n==3).or.(n==5)) irev_normal(n)=.true.
             endif
          endif

          ! Case #4: ex pour imax (n=2)
          ! ========            7---------------8   distance min between A and corners of neighboring domain
          !                    / |              /|   index imA=1
          !            B      /  |             / |  distance min between B and corners of neighboring domain
          !            |     2----------------4  |   index imB=2
          !            |     |   | i ^  k     |  |  distance min between C and corners of neighboring domain
          !    j  k    |     |   |   | /      |  |   index imC=5
          !    | /     |   C |   |   |/       |  |  is_swapij=.true.
          !    |/___i  |  /  |   |    ---> j  |  |  is_swapik=.false.
          !            | /   |   5----------- |--6  is_swapjk=.false.
          !            |/    |  /             | /   is_rev(k,2)=.false. for normal dir.
          !   D--------A     | /              |/    is_rev(k,1)=.false. for parallel dir. #1
          !                  1----------------3     is_rev(k,3)=.false. for parallel dir. #2

          ! contact with jmin face [indirect orientation]
          ! ---------------------------------------------
          if ((imA==1).and.(imB==2).and.(imC==5)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 3
                n_dir(n)=3
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_ij(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_jk(n)=.true.
                ! no reverse in parallel direction
                !
                ! reverse in normal direction
                if ((n==1).or.(n==3).or.(n==5)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==2).and.(imB==7).and.(imC==1)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 3
                n_dir(n)=3
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_ijk(n)=.true.
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_ik(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ijk(n)=.true.
                ! reverse in parallel direction
                irev_parallel2(n)=.true. ! for parallel dir. #2
                ! reverse in normal direction
                if ((n==1).or.(n==3).or.(n==5)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==7).and.(imB==5).and.(imC==2)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 3
                n_dir(n)=3
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_ij(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_jk(n)=.true.
                ! no reverse in parallel direction
                irev_parallel1(n)=.true. ! for parallel dir. #1
                irev_parallel2(n)=.true. ! for parallel dir. #2
                ! reverse in normal direction
                if ((n==1).or.(n==3).or.(n==5)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==5).and.(imB==1).and.(imC==7)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 3
                n_dir(n)=3
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_ijk(n)=.true.
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_ik(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ijk(n)=.true.
                ! reverse in parallel direction
                irev_parallel1(n)=.true. ! for parallel dir. #1
                ! reverse in normal direction
                if ((n==1).or.(n==3).or.(n==5)) irev_normal(n)=.true.
             endif
          endif

          ! Case #5: ex pour imax (n=2)
          ! ========            7---------------5   distance min between A and corners of neighboring domain
          !                    / |              /|   index imA=4
          !            B      /  |             / |  distance min between B and corners of neighboring domain
          !            |     2----------------1  |   index imB=2
          !            |     |   |       /k   |  |  distance min between C and corners of neighboring domain
          !    j  k    |     |   |  <---/     |  |   index imC=8
          !    | /     |   C |   |  i   |     |  |  is_swapij=.false.
          !    |/___i  |  /  |   |      |j    |  |  is_swapik=.false.
          !            | /   |   8----------- |--6  is_swapjk=.false.
          !            |/    |  /             | /   is_rev(k,2)=.true.  for normal dir.
          !   D--------A     | /              |/    is_rev(k,1)=.true.  for parallel dir. #1
          !                  4----------------3     is_rev(k,3)=.false. for parallel dir. #2

          ! contact with imax face [direct orientation]
          ! -------------------------------------------
          if ((imA==4).and.(imB==2).and.(imC==8)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 2
                n_dir(n)=2
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_ij(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ijk(n)=.true.
                ! reverse in parallel direction
                irev_parallel1(n)=.true. ! for parallel dir. #1
                ! reverse in normal direction
                if ((n==2).or.(n==4).or.(n==6)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==2).and.(imB==7).and.(imC==4)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 2
                n_dir(n)=2
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_jk(n)=.true.
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_ijk(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ik(n)=.true.
                ! reverse in parallel direction
                !
                ! reverse in normal direction
                if ((n==2).or.(n==4).or.(n==6)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==7).and.(imB==8).and.(imC==2)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 2
                n_dir(n)=2
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_ij(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ijk(n)=.true.
                ! reverse in parallel direction
                irev_parallel2(n)=.true. ! for parallel dir. #2
                ! reverse in normal direction
                if ((n==2).or.(n==4).or.(n==6)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==8).and.(imB==4).and.(imC==7)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 2
                n_dir(n)=2
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_jk(n)=.true.
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_ijk(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ik(n)=.true.
                ! reverse in parallel directions
                irev_parallel1(n)=.true. ! for parallel dir. #1
                irev_parallel2(n)=.true. ! for parallel dir. #2
                ! reverse in normal direction
                if ((n==2).or.(n==4).or.(n==6)) irev_normal(n)=.true.
             endif
          endif

          ! Case #6: ex pour imax (n=2)
          ! ========            8---------------6   distance min between A and corners of neighboring domain
          !                    / |              /|   index imA=2
          !            B      /  |             / |  distance min between B and corners of neighboring domain
          !            |     4----------------3  |   index imB=4
          !            |     |   |   j ^  k   |  |  distance min between C and corners of neighboring domain
          !    j  k    |     |   |     | /    |  |   index imC=7
          !    | /     |   C |   |     |/     |  |  is_swapij=.false.
          !    |/___i  |  /  |   | i<---      |  |  is_swapik=.false.
          !            | /   |   7----------- |--5  is_swapjk=.false.
          !            |/    |  /             | /   is_rev(k,2)=.true. for normal dir.
          !   D--------A     | /              |/    is_rev(k,1)=.false. for parallel dir. #1
          !                  2----------------1     is_rev(k,3)=.false. for parallel dir. #2

          ! contact with imax face [indirect orientation]
          ! ---------------------------------------------
          if ((imA==2).and.(imB==4).and.(imC==7)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 2
                n_dir(n)=2
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_ij(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ijk(n)=.true.
                ! reverse in parallel direction
                !
                ! reverse in normal direction
                if ((n==2).or.(n==4).or.(n==6)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==4).and.(imB==8).and.(imC==2)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 2
                n_dir(n)=2
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_jk(n)=.true.
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_ijk(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ik(n)=.true.
                ! reverse in parallel direction
                irev_parallel2(n)=.true. ! for parallel dir. #2
                ! reverse in normal direction
                if ((n==2).or.(n==4).or.(n==6)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==8).and.(imB==7).and.(imC==4)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 2
                n_dir(n)=2
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_ij(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ijk(n)=.true.
                ! reverse in parallel directions
                irev_parallel1(n)=.true. ! for parallel dir. #1
                irev_parallel2(n)=.true. ! for parallel dir. #2
                ! reverse in normal direction
                if ((n==2).or.(n==4).or.(n==6)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==7).and.(imB==2).and.(imC==8)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 2
                n_dir(n)=2
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_jk(n)=.true.
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_ijk(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ik(n)=.true.
                ! reverse in parallel direction
                irev_parallel1(n)=.true. ! for parallel dir. #1
                ! reverse in normal direction
                if ((n==2).or.(n==4).or.(n==6)) irev_normal(n)=.true.
             endif
          endif

          ! Case #7: ex pour imax (n=2)
          ! ========            8---------------7   distance min between A and corners of neighboring domain
          !                    / |              /|   index imA=3
          !            B      /  |             / |  distance min between B and corners of neighboring domain
          !            |     4----------------2  |   index imB=4
          !            |     |   |   i ^  k   |  |  distance min between C and corners of neighboring domain
          !    j  k    |     |   |     | /    |  |   index imC=6
          !    | /     |   C |   |     |/     |  |  is_swapij=.true.
          !    |/___i  |  /  |   | j<---      |  |  is_swapik=.false.
          !            | /   |   6----------- |--5  is_swapjk=.false.
          !            |/    |  /             | /   is_rev(k,2)=.true.  for normal dir.
          !   D--------A     | /              |/    is_rev(k,1)=.false. for parallel dir. #1
          !                  3----------------1     is_rev(k,3)=.false. for parallel dir. #2

          ! contact with jmax face [direct orientation]
          ! -------------------------------------------
          if ((imA==3).and.(imB==4).and.(imC==6)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 4
                n_dir(n)=4
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_ij(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_jk(n)=.true.
                ! no reverse in parallel direction
                !
                ! reverse in normal direction
                if ((n==2).or.(n==4).or.(n==6)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==4).and.(imB==8).and.(imC==3)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 4
                n_dir(n)=4
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_ijk(n)=.true.
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_ik(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ijk(n)=.true.
                ! reverse in parallel direction
                irev_parallel2(n)=.true. ! for parallel dir. #2
                ! reverse in normal direction
                if ((n==2).or.(n==4).or.(n==6)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==8).and.(imB==6).and.(imC==4)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 4
                n_dir(n)=4
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_ij(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_jk(n)=.true.
                ! reverse in parallel directions
                irev_parallel1(n)=.true. ! for parallel dir. #1
                irev_parallel2(n)=.true. ! for parallel dir. #2
                ! reverse in normal direction
                if ((n==2).or.(n==4).or.(n==6)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==6).and.(imB==3).and.(imC==8)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 4
                n_dir(n)=4
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_ijk(n)=.true.
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_ik(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ijk(n)=.true.
                ! reverse in parallel direction
                irev_parallel1(n)=.true. ! for parallel dir. #1
                ! reverse in normal direction
                if ((n==2).or.(n==4).or.(n==6)) irev_normal(n)=.true.
             endif
          endif

          ! Case #8: ex pour imax (n=2)
          ! ========            6---------------5   distance min between A and corners of neighboring domain
          !                    / |              /|   index imA=4
          !            B      /  |             / |  distance min between B and corners of neighboring domain
          !            |     3----------------1  |   index imB=8
          !            |     |   |      /k    |  |  distance min between C and corners of neighboring domain
          !    j  k    |     |   |     /      |  |   index imC=3
          !    | /     |   C |   |j<---|      |  |  is_swapij=.true.
          !    |/___i  |  /  |   |     |      |  |  is_swapik=.false.
          !            | /   |   8---- i -----|--7  is_swapjk=.false.
          !            |/    |  /             | /   is_rev(k,2)=.true.  for normal dir.
          !   D--------A     | /              |/    is_rev(k,1)=.true.  for parallel dir. #1
          !                  4----------------2     is_rev(k,3)=.false. for parallel dir. #2

          ! contact with jmax face [indirect orientation]
          ! ---------------------------------------------
          if ((imA==4).and.(imB==3).and.(imC==8)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 4
                n_dir(n)=4
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_ij(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_jk(n)=.true.
                ! reverse in parallel direction
                irev_parallel1(n)=.true. ! for parallel dir. #1
                ! reverse in normal direction
                if ((n==2).or.(n==4).or.(n==6)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==3).and.(imB==6).and.(imC==4)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 4
                n_dir(n)=4
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_ijk(n)=.true.
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_ik(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ijk(n)=.true.
                ! no reverse in parallel direction
                !
                ! reverse in normal direction
                if ((n==2).or.(n==4).or.(n==6)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==6).and.(imB==8).and.(imC==3)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 4
                n_dir(n)=4
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_ij(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_jk(n)=.true.
                ! reverse in parallel directions
                irev_parallel2(n)=.true. ! for parallel dir. #2
                ! reverse in normal direction
                if ((n==2).or.(n==4).or.(n==6)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==8).and.(imB==4).and.(imC==6)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 4
                n_dir(n)=4
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_ijk(n)=.true.
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_ik(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ijk(n)=.true.
                ! reverse in parallel directions
                irev_parallel1(n)=.true. ! for parallel dir. #1
                irev_parallel2(n)=.true. ! for parallel dir. #2
                ! reverse in normal direction
                if ((n==2).or.(n==4).or.(n==6)) irev_normal(n)=.true.
             endif
          endif

          ! Case #9: ex pour imax (n=2)
          ! ========            3---------------6   distance min between A and corners of neighboring domain
          !                    / |              /|   index imA=2
          !            B      /  |             / |  distance min between B and corners of neighboring domain
          !            |     4-------^--------8  |   index imB=4
          !            |     |   | j |        |  |  distance min between C and corners of neighboring domain
          !    j  k    |     |   |   |---> k  |  |   index imC=1
          !    | /     |   C |   |   /        |  |  is_swapij=.false.
          !    |/___i  |  /  |   |  /         |  |  is_swapik=.true.
          !            | /   |   1- i --------|--5  is_swapjk=.false.
          !            |/    |  /             | /   is_rev(k,2)=.false. for normal dir.
          !   D--------A     | /              |/    is_rev(k,1)=.false. for parallel dir. #1
          !                  2----------------7     is_rev(k,3)=.true.  for parallel dir. #2

          ! contact with kmin face [direct orientation]
          ! -------------------------------------------
          if ((imA==2).and.(imB==4).and.(imC==1)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 5
                n_dir(n)=5
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_ik(n)=.true.
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_ijk(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ij(n)=.true.
                ! reverse in parallel direction
                irev_parallel2(n)=.true. ! for parallel dir. #2
                ! reverse in normal direction
                if ((n==1).or.(n==3).or.(n==5)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==4).and.(imB==3).and.(imC==2)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 5
                n_dir(n)=5
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_ijk(n)=.true.
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_jk(n)=.true.
                ! reverse in parallel directions
                irev_parallel1(n)=.true. ! for parallel dir. #1
                irev_parallel2(n)=.true. ! for parallel dir. #2
                ! reverse in normal direction
                if ((n==1).or.(n==3).or.(n==5)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==3).and.(imB==1).and.(imC==4)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 5
                n_dir(n)=5
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_ik(n)=.true.
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_ijk(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ij(n)=.true.
                ! reverse in parallel direction
                irev_parallel1(n)=.true. ! for parallel dir. #1
                ! reverse in normal direction
                if ((n==1).or.(n==3).or.(n==5)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==1).and.(imB==2).and.(imC==3)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 5
                n_dir(n)=5
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_ijk(n)=.true.
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_jk(n)=.true.
                ! no reverse in parallel direction
                !
                ! reverse in normal direction
                if ((n==1).or.(n==3).or.(n==5)) irev_normal(n)=.true.
             endif
          endif

          ! Case #10: ex pour imax (n=2)
          ! ========            1---------------5   distance min between A and corners of neighboring domain
          !                    / |              /|   index imA=4
          !            B      /  |             / |  distance min between B and corners of neighboring domain
          !            |     2----------------7  |   index imB=2
          !            |     |   |     ---> k |  |  distance min between C and corners of neighboring domain
          !    j  k    |     |   |    /|      |  |   index imC=3
          !    | /     |   C |   | i / |      |  |  is_swapij=.false.
          !    |/___i  |  /  |   |      j     |  |  is_swapik=.true.
          !            | /   |   3----------- |--6  is_swapjk=.false.
          !            |/    |  /             | /   is_rev(k,2)=.false. for normal dir.
          !   D--------A     | /              |/    is_rev(k,1)=.true.  for parallel dir. #1
          !                  4----------------8     is_rev(k,3)=.true.  for parallel dir. #2

          ! contact with kmin face [indirect orientation]
          ! ---------------------------------------------
          if ((imA==4).and.(imB==2).and.(imC==3)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 5
                n_dir(n)=5
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_ik(n)=.true.
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_ijk(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ij(n)=.true.
                ! reverse in parallel directions
                irev_parallel1(n)=.true. ! for parallel dir. #1
                irev_parallel2(n)=.true. ! for parallel dir. #2
                ! reverse in normal direction
                if ((n==1).or.(n==3).or.(n==5)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==2).and.(imB==1).and.(imC==4)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 5
                n_dir(n)=5
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_ijk(n)=.true.
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_jk(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ij(n)=.true.
                ! reverse in parallel directions
                irev_parallel1(n)=.true. ! for parallel dir. #1
                ! reverse in normal direction
                if ((n==1).or.(n==3).or.(n==5)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==1).and.(imB==3).and.(imC==2)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 5
                n_dir(n)=5
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_ik(n)=.true.
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_ijk(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ij(n)=.true.
                ! no reverse in parallel direction
                !
                ! reverse in normal direction
                if ((n==1).or.(n==3).or.(n==5)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==3).and.(imB==4).and.(imC==1)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 5
                n_dir(n)=5
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_ijk(n)=.true.
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_jk(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ij(n)=.true.
                ! reverse in parallel direction
                irev_parallel2(n)=.true. ! for parallel dir. #2
                ! reverse in normal direction
                if ((n==1).or.(n==3).or.(n==5)) irev_normal(n)=.true.
             endif
          endif

          ! Case #11: ex pour imax (n=2)
          ! ========            7---------------2   distance min between A and corners of neighboring domain
          !                    / |              /|   index imA=6
          !            B      /  |             / |  distance min between B and corners of neighboring domain
          !            |     8----------------4  |   index imB=8
          !            |     |   |    i ^     |  |  distance min between C and corners of neighboring domain
          !    j  k    |     |   |      |     |  |   index imC=5
          !    | /     |   C |   |  k __|     |  |  is_swapij=.true.
          !    |/___i  |  /  |   |      /     |  |  is_swapik=.true. => is_swapijk=.true.
          !            | /   |   5-----/------|--1  is_swapjk=.true.
          !            |/    |  /     j       | /   is_rev(k,2)=.true.  for normal dir.
          !   D--------A     | /              |/    is_rev(k,1)=.false. for parallel dir. #1
          !                  6----------------3     is_rev(k,3)=.true.  for parallel dir. #2

          ! contact with kmax face [direct orientation]
          ! -------------------------------------------
          if ((imA==6).and.(imB==8).and.(imC==5)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 6
                n_dir(n)=6
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_ijk(n)=.true.
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_jk(n)=.true. 
                ! reverse in parallel direction
                irev_parallel2(n)=.true. ! for parallel dir. #2
                ! reverse in normal direction
                if ((n==2).or.(n==4).or.(n==6)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==8).and.(imB==7).and.(imC==6)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 6
                n_dir(n)=6
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_ik(n)=.true.
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_ijk(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ij(n)=.true.
                ! reverse in parallel directions
                irev_parallel1(n)=.true. ! for parallel dir. #1
                irev_parallel2(n)=.true. ! for parallel dir. #2
                ! reverse in normal direction
                if ((n==2).or.(n==4).or.(n==6)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==7).and.(imB==5).and.(imC==8)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 6
                n_dir(n)=6
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_ijk(n)=.true.
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_jk(n)=.true. 
                ! reverse in parallel direction
                irev_parallel1(n)=.true. ! for parallel dir. #1
                ! reverse in normal direction
                if ((n==2).or.(n==4).or.(n==6)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==5).and.(imB==6).and.(imC==7)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 6
                n_dir(n)=6
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_ik(n)=.true.
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_ijk(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ij(n)=.true.
                ! no reverse in parallel direction
                !
                ! reverse in normal direction
                if ((n==2).or.(n==4).or.(n==6)) irev_normal(n)=.true.
             endif
          endif

          ! Case #12: ex pour imax (n=2)
          ! ========            5---------------1   distance min between A and corners of neighboring domain
          !                    / |              /|   index imA=8
          !            B      /  |             / |  distance min between B and corners of neighboring domain
          !            |     6----------------3  |   index imB=6
          !            |     |   |  k____     |  |  distance min between C and corners of neighboring domain
          !    j  k    |     |   |      /|    |  |   index imC=7
          !    | /     |   C |   |     / |    |  |  is_swapij=.true.
          !    |/___i  |  /  |   |    j  i    |  |  is_swapik=.true. => is_swapijk=.true.
          !            | /   |   7----------- |--2  is_swapjk=.true.
          !            |/    |  /             | /   is_rev(k,2)=.true. for normal dir.
          !   D--------A     | /              |/    is_rev(k,1)=.true. for parallel dir. #1
          !                  8----------------4     is_rev(k,3)=.true. for parallel dir. #2

          ! contact with kmax face [indirect orientation]
          ! ---------------------------------------------
          if ((imA==8).and.(imB==6).and.(imC==7)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 6
                n_dir(n)=6
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_ijk(n)=.true.
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_jk(n)=.true. 
                ! reverse in parallel direction
                irev_parallel1(n)=.true. ! for parallel dir. #1
                irev_parallel2(n)=.true. ! for parallel dir. #2
                ! reverse in normal direction
                if ((n==2).or.(n==4).or.(n==6)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==6).and.(imB==5).and.(imC==8)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 6
                n_dir(n)=6                
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_ik(n)=.true.
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_ijk(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ij(n)=.true.
                ! reverse in parallel directions
                irev_parallel1(n)=.true. ! for parallel dir. #1
                ! reverse in normal direction
                if ((n==2).or.(n==4).or.(n==6)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==5).and.(imB==7).and.(imC==6)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 6
                n_dir(n)=6
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_ijk(n)=.true.
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_jk(n)=.true. 
                ! no reverse in parallel direction
                !
                ! reverse in normal direction
                if ((n==2).or.(n==4).or.(n==6)) irev_normal(n)=.true.
             endif
          endif

          if ((imA==7).and.(imB==8).and.(imC==5)) then
             ! if (distA(imA)>distA_D) then
             !    ! periodic condition: skip following (?)
             !    bl(nbl)%periodic(n)=.true.
             ! else
             if ((distA(imA).le.distA_D).and.(.not.bl(nbl)%periodic(n))) then
                ! direction of neighbor is 6
                n_dir(n)=6
                ! if faces 1 or 2
                if (n==1.or.n==2) iswap_ik(n)=.true.
                ! if faces 3 or 4
                if (n==3.or.n==4) iswap_ijk(n)=.true.
                ! if faces 5 or 6
                if (n==5.or.n==6) iswap_ij(n)=.true.
                ! reverse in parallel direction
                 irev_parallel2(n)=.true. ! for parallel dir. #2               
                ! reverse in normal direction
                if ((n==2).or.(n==4).or.(n==6)) irev_normal(n)=.true.
             endif
          endif

       endif

    enddo

  end subroutine check_directions_3d

end module mod_grid_directions_3d
