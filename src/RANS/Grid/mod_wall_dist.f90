!=================================================================================
module mod_wall_dist
!=================================================================================
  !> Module to compute distance to wall
!=================================================================================
  use mod_mpi
  use mod_flow
  implicit none
  !-------------------------------------------------------------------------------
  ! Wall distance computation
  ! =========================
  ! Local and shared wall coordinates for curvilinear extruded
  real(wp), dimension(:), allocatable :: xw,yw,zw,xwg,ywg,zwg

  ! Nb of elements to send and starting points inside buffer
  integer, dimension(:), allocatable :: counts,disp

  ! Array that indicates presence of wall for each proc, total nb of walls
  logical, dimension(:), allocatable :: walls
  integer                            :: nwall

  ! Arrays containing distances; d_ has wall distance from every wall at same loc,
  ! d is min wall distance for one wall loc, dg is min wall distance whatever the wall loc
  ! dk is used if wall is at k location and is_curv extruded
  ! dij is used if wall is at i/j and       "              "
  ! if is_curv3, d_ is used instead
  real(wp), dimension(:,:),     allocatable :: dk
  real(wp), dimension(:,:,:),   allocatable :: d,dij
  real(wp), dimension(:,:,:,:), allocatable :: dg,d_

  ! h_wn: grid-step in wall-normal direction to be used in IDDES
  real(wp), dimension(:,:),   allocatable :: h_wn
  ! ------------------------------------------------------------------------------

contains

  !==============================================================================
  subroutine wall_dist
  !==============================================================================
    !> Compute distance to closest wall
  !==============================================================================
    use mod_mpi
    use mod_grid
    use mod_flow
    use mod_constant
    use warnstop
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: ind,i,j,k,m,ind_loc
    logical, dimension(6) :: is_wall ! -> are there any walls located at i-j-kmin-max ?
    integer, dimension(6) :: locs_ijk,locs_minmax ! -> indices of is_bc_wall(locs_ijk,locs_minmax)
    integer, dimension(6) :: closest_wall
    character(len=6)  :: iproc_str
    character(len=30) :: wall_dir
    logical :: is_dir
    ! --------------------------------------------------------------------------
  
    if (iproc.eq.0) print*,'Computing distance to wall...'

    locs_ijk    = (/1,1,2,2,3,3/)
    locs_minmax = (/1,2,1,2,1,2/)
  
    ! Initialize with true
    is_wall=.true.
    allocate(walls(nproc))
    do m=1,6
       ! Check how may procs have walls at each location
       call MPI_ALLGATHER(is_bc_wall(locs_ijk(m),locs_minmax(m)),1,MPI_LOGICAL,walls,1,MPI_LOGICAL,COMM_global,info)
       nwall = count(walls,dim=1)
       if (nwall.eq.0) is_wall(m)=.false.
    enddo
    deallocate(walls)
  
    ind = count(is_wall,dim=1) ! -> nb of wall types (e.g =3 if 3 walls in same proc exist, wherever they are)
    allocate(dg(nx,ny,nz,ind))
    dg =0.0_wp
    ind=0
 
    ! Compute wall distance in every direction 
    call compute_wall_dist('imin',ind)
    call compute_wall_dist('imax',ind)
    call compute_wall_dist('jmin',ind)
    call compute_wall_dist('jmax',ind)
    call compute_wall_dist('kmin',ind)
    call compute_wall_dist('kmax',ind)

    ! Get distance from each point to the closest wall
    allocate(d(nx,ny,nz))
    do k=1,nz
       do j=1,ny
          do i=1,nx
             d(i,j,k) = minval(dg(i,j,k,:))
          enddo
       enddo
    enddo
    !! particular cases
    !! in case of inviscid wall in the first block
    !if ((nob(iproc).eq.1).and.is_slip_in) then
    !   do j=1,ny
    !      do i=1,nx!
    !         d(i,j,k) = sqrt((xc(i,j)-(-110.0_wp))**2+(yc(i,j)-1.0_wp)**2)
    !         d(i,j,k) = min( d(i,j,k),sqrt((xc(i,j)-(-110.0_wp))**2+(yc(i,j)-9.0_wp)**2) )
!   !          d(i,j,k) = sqrt((xc(i,j)-(-109.192226143_wp))**2+(yc(i,j)-1.0_wp)**2)
!   !          d(i,j,k) = min( d(i,j,k),sqrt((xc(i,j)-(-109.192226143_wp))**2+(yc(i,j)-9.0_wp)**2) )
    !      enddo
    !   enddo
    !endif
    ! in case the corner of the proc locates at wall without any wall BC (such as cavity edges)

    ! /!\ WORK ON THAT AND EXTEND TO 3D /!\ !    

    !if((d( 1, 1, 1).eq.0.0_wp).and.not(is_bc_wall(1,1)).and.not(is_bc_wall(2,1))) d( 1, 1, 1) = eps
    !if((d( 1,ny, 1).eq.0.0_wp).and.not(is_bc_wall(1,1)).and.not(is_bc_wall(2,2))) d( 1,ny, 1) = eps
    !if((d(nx, 1, 1).eq.0.0_wp).and.not(is_bc_wall(1,2)).and.not(is_bc_wall(2,1))) d(nx, 1, 1) = eps
    !if((d(nx,ny, 1).eq.0.0_wp).and.not(is_bc_wall(1,2)).and.not(is_bc_wall(2,2))) d(nx,ny, 1) = eps
    !if((d(nx,ny, 1).eq.0.0_wp).and.not(is_bc_wall(1,2)).and.not(is_bc_wall(2,2))) d(nx,ny, 1) = eps
    !if((d(nx,ny, 1).eq.0.0_wp).and.not(is_bc_wall(1,2)).and.not(is_bc_wall(2,2))) d(nx,ny, 1) = eps

    ! compute grid step in the wall-normal direction to be used in IDDES
    if (simulation_RANS.eq.'IDDES') then
      if (is_curv3) then
          call mpistop("IDDES not implemented in 3D curvilinear...",1)
      else if (is_curv) then
         closest_wall = 0
         i = 0
         do k=1,4
            if (is_wall(k)) then
               i = i + 1
               closest_wall(k) = i
            endif
         enddo

         allocate(h_wn(nx,ny))
         h_wn = 0.0_wp
         k=1
         do j=1,ny
            do i=1,nx
               ind_loc = minloc(dg(i,j,k,:), dim=1)
               if (closest_wall(1).eq.ind_loc) then
                  h_wn(i,j) = (xc(i+1,j)-xc(i,j))**2 + (yc(i+1,j)-yc(i,j))**2
               elseif (closest_wall(2).eq.ind_loc) then
                  h_wn(i,j) = (xc(i,j)-xc(i-1,j))**2 + (yc(i,j)-yc(i-1,j))**2
               elseif (closest_wall(3).eq.ind_loc) then
                  h_wn(i,j) = (xc(i,j+1)-xc(i,j))**2 + (yc(i,j+1)-yc(i,j))**2
               elseif (closest_wall(4).eq.ind_loc) then
                  h_wn(i,j) = (xc(i,j)-xc(i,j-1))**2 + (yc(i,j)-yc(i,j-1))**2
               endif
            enddo
         enddo
      else
          call mpistop("RANS not implemented in cartesian...",1)
      endif
    endif
 
    ! Write to file (very clumsy way) TEMP
    ! ====================================
    wall_dir = 'Wall_distance_check'
    write(iproc_str,'(I6.6)') iproc
    if (iproc.eq.0) then
       inquire(directory=wall_dir,exist=is_dir)
       if (.not.is_dir) call system('mkdir '//trim(wall_dir))
    endif
    call MPI_BARRIER(COMM_global,info)

    open(iproc,file=trim(wall_dir)//'/wall_dist_proc'//trim(iproc_str)//'.bin',status='replace',form='unformatted')
    write(iproc) nx,ny,nz
    if (is_curv) then
       write(iproc) d(1:nx,1:ny,1),xc(1:nx,1:ny),yc(1:nx,1:ny)
    elseif (is_curv3) then
       write(iproc) d(1:nx,1:ny,1:nz),xc3(1:nx,1:ny,1:nz),yc3(1:nx,1:ny,1:nz),zc3(1:nx,1:ny,1:nz)
    endif
    close(iproc)
   
    !open(iproc,file=trim(wall_dir)//'/wall_dist_proc'//trim(iproc_str)//'.x',status='replace',form='formatted')
    !if (is_curv) then 
    !   write(iproc,*) nx,ny,1
    !   write(iproc,*) ((xc(i,j),i=1,nx),j=1,ny)
    !   write(iproc,*) ((yc(i,j),i=1,nx),j=1,ny)
    !   write(iproc,*) (( d(i,j,1),i=1,nx),j=1,ny)
    !elseif (is_curv3) then
    !   write(iproc,*) nx,ny,nz
    !   write(iproc,*) (((xc3(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    !   write(iproc,*) (((yc3(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    !   write(iproc,*) (((zc3(i,j,k),i=1,nx),j=1,ny),k=1,nz)
    !   write(iproc,*) (((  d(i,j,k),i=1,nx),j=1,ny),k=1,nz)   
    !endif
 
    if (iproc.eq.0) print*,'Distance to wall OK'
    deallocate(dg)
  
  end subroutine wall_dist
   
  !==============================================================================
  subroutine compute_wall_dist(loc,ind)
  !==============================================================================
    !> Compute distance to closest wall when it is located at imin
  !==============================================================================
    use mod_mpi
    use mod_flow
    use warnstop
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k,m
    integer :: ind
    integer :: loc_ijk,loc_minmax
    integer :: nw
    integer :: ndxw,nfxw,ndyw,nfyw,ndzw,nfzw
    integer, dimension(:), allocatable :: nwg
    character(len=4)  :: loc
    ! --------------------------------------------------------------------------

    ! Index meaning:
    ! loc_ijk,loc_minmax -> see definition in subroutine wall_dist
    ! nw          -> nb of points at the wall for each proc, regardless of wall location
    ! ndxw,nfxw   -> start and finish for wall coordinates that will depend on where the
    !                wall is located (e.g for imin, ndxw=nfxw=1,ndyw=1,nfyw=ny,ndzw=1,nfzw=nz
    !                for jmax, ndxw=1,nfxw=nx,ndyw=nfyw=ny,ndzw=1,nfzw=nz)
    ! ind         -> used in case there are no walls at imin, imax, jmin, jmax, kmin or kmax
 

    if (loc.eq.'imin') then 
       loc_ijk    =1
       loc_minmax =1
       nw    =ny*nz
       ndxw  =1
       nfxw  =1
       ndyw  =1
       nfyw  =ny
       ndzw  =1
       nfzw  =nz
    
    elseif (loc.eq.'imax') then 
       loc_ijk    =1
       loc_minmax =2
       nw    =ny*nz
       ndxw  =nx
       nfxw  =nx
       ndyw  =1
       nfyw  =ny
       ndzw  =1
       nfzw  =nz
    
    elseif (loc.eq.'jmin') then 
       loc_ijk    =2
       loc_minmax =1
       nw    =nx*nz
       ndxw  =1
       nfxw  =nx
       ndyw  =1
       nfyw  =1
       ndzw  =1
       nfzw  =nz
  
    elseif (loc.eq.'jmax') then 
       loc_ijk    =2
       loc_minmax =2
       nw    =nx*nz
       ndxw  =1
       nfxw  =nx
       ndyw  =ny
       nfyw  =ny
       ndzw  =1
       nfzw  =nz
    
    elseif (loc.eq.'kmin') then 
       loc_ijk    =3
       loc_minmax =1
       nw    =nx*ny
       ndxw  =1
       nfxw  =nx
       ndyw  =1
       nfyw  =ny
       ndzw  =1
       nfzw  =1
    
    elseif (loc.eq.'kmax') then 
       loc_ijk    =3
       loc_minmax =2
       nw    =nx*ny
       ndxw  =1
       nfxw  =nx
       ndyw  =1
       nfyw  =ny
       ndzw  =nz
       nfzw  =nz
  
    endif

    ! Check how many procs have walls at location
    allocate(walls(nproc))
    call MPI_ALLGATHER(is_bc_wall(loc_ijk,loc_minmax),1,MPI_LOGICAL,walls,1,MPI_LOGICAL,COMM_global,info)
    nwall = count(walls,dim=1)

    ! If proc has at least 1 wall:
    if (nwall.gt.0) then
       ind=ind+1
       allocate(disp(0:nproc),counts(nproc))

       ! Case curvilinear extruded in 3rd direction
       ! ==========================================
       if (is_curv) then

          ! Allocate arrays and get wall coordinates
          ! ----------------------------------------
          if (is_bc_wall(loc_ijk,loc_minmax)) then
             ! Store wall coords
             if (loc_ijk.eq.1) then     ! imin-imax
                nw=ny
                allocate(xw(nw),yw(nw),zw(1))
                xw(:) = xc(ndxw,ndyw:nfyw)
                yw(:) = yc(ndxw,ndyw:nfyw)
                zw = 0.0_wp
             elseif (loc_ijk.eq.2) then ! jmin-jmax
                nw=nx
                allocate(xw(nw),yw(nw),zw(1))
                xw(:) = xc(ndxw:nfxw,ndyw)
                yw(:) = yc(ndxw:nfxw,ndyw)
                zw = 0.0_wp
             elseif (loc_ijk.eq.3) then ! kmin-kmax
                nw=0
                allocate(xw(nw),yw(nw),zw(1))
                zw = z(ndzw)
             endif

          else
             nw = 0
             allocate(xw(nw),yw(nw),zw(1))
             zw = 0.0_wp
          endif

          ! Every proc needs the counts and disp array
          ! ------------------------------------------
          allocate(nwg(nproc))
          call MPI_ALLGATHER(nw,1,MPI_INTEGER,nwg,1,MPI_INTEGER,COMM_global,info)
          disp = 0
          do i=1,nproc
             if (walls(i)) then
                ! Store displacement for receive buffer and size of message
                disp(i)   = nwg(i)+disp(i-1)
                counts(i) = nwg(i)
             else
                disp(i)   = 0+disp(i-1)
                counts(i) = 0
             endif
          enddo

          ! Give approriate size to total wall coordinate arrays: size = sum(nwg)
          allocate(xwg(sum(nwg)),ywg(sum(nwg)),zwg(nproc))

          ! Share wall coords
          ! -----------------
          call MPI_ALLGATHERV(xw,counts(iproc+1),MPI_DOUBLE_PRECISION,xwg,counts,disp,MPI_DOUBLE_PRECISION,COMM_global,info)
          call MPI_ALLGATHERV(yw,counts(iproc+1),MPI_DOUBLE_PRECISION,ywg,counts,disp,MPI_DOUBLE_PRECISION,COMM_global,info)
          call MPI_ALLGATHER(zw(1),1,MPI_DOUBLE_PRECISION,zwg,1,MPI_DOUBLE_PRECISION,COMM_global,info)
          
          ! Compute min distance to wall for every point
          ! --------------------------------------------
        
          if (loc_ijk.ne.3) then 
             allocate(dij(nx,ny,size(xwg)),d(nx,ny,nz))
             do i=1,nx
                do j=1,ny
                   do m=1,size(xwg)
                      dij(i,j,m) = sqrt((xc(i,j)-xwg(m))**2 + (yc(i,j)-ywg(m))**2)
                   enddo
                   d(i,j,:) = minval(dij(i,j,:))
                enddo
             enddo
          else
             allocate(dk(nz,size(zwg)),d(nx,ny,nz))
             ! store distance wrt to z-wall
             do k=1,nz
                do m=1,size(zwg)
                   dk(k,m) = abs(z(k)-zwg(m))
                enddo
                d(:,:,k) = minval(dk(k,:))
             enddo
          endif

       ! Case curvilinear in all 3 directions
       ! ====================================
       elseif (is_curv3) then

          ! Allocate arrays and get wall coordinates
          ! ----------------------------------------
          if (is_bc_wall(loc_ijk,loc_minmax)) then
             allocate(xw(nw),yw(nw),zw(nw))
             ! Store wall coords
             if (loc_ijk.eq.1) then      ! imin-imax
                xw(:) = pack(xc3(ndxw,ndyw:nfyw,ndzw:nfzw),.true.)
                yw(:) = pack(yc3(ndxw,ndyw:nfyw,ndzw:nfzw),.true.)
                zw(:) = pack(zc3(ndxw,ndyw:nfyw,ndzw:nfzw),.true.)
             elseif (loc_ijk.eq.2) then  ! jmin-jmax
                xw(:) = pack(xc3(ndxw:nfxw,ndyw,ndzw:nfzw),.true.)
                yw(:) = pack(yc3(ndxw:nfxw,ndyw,ndzw:nfzw),.true.)
                zw(:) = pack(zc3(ndxw:nfxw,ndyw,ndzw:nfzw),.true.)
             elseif (loc_ijk.eq.3) then  ! kmin-kmax
                xw(:) = pack(xc3(ndxw:nfxw,ndyw:nfyw,ndzw),.true.)
                yw(:) = pack(yc3(ndxw:nfxw,ndyw:nfyw,ndzw),.true.)
                zw(:) = pack(zc3(ndxw,ndyw:ndyw:nfyw,ndzw),.true.)
             endif
          else
             nw = 0
             allocate(xw(nw),yw(nw),zw(nw))
             xw = 0.0_wp
             yw = 0.0_wp
             zw = 0.0_wp
          endif
          
          ! Every proc needs the counts and disp array
          ! ------------------------------------------
          allocate(nwg(nproc))
          call MPI_ALLGATHER(nw,1,MPI_INTEGER,nwg,1,MPI_INTEGER,COMM_global,info)
          disp = 0
          do i=1,nproc
             if (walls(i)) then
                ! Store displacement for receive buffer and size of message
                disp(i)   = nwg(i)+disp(i-1)
                counts(i) = nwg(i)
             else
                disp(i)   = 0+disp(i-1)
                counts(i) = 0
             endif
          enddo

          ! Give approriate size to total wall coordinate arrays: size = sum(nwg)
          allocate(xwg(sum(nwg)),ywg(sum(nwg)),zwg(sum(nwg)))

          ! Share wall coords
          ! -----------------
          call MPI_ALLGATHERV(xw,counts(iproc+1),MPI_DOUBLE_PRECISION,xwg,counts,disp,MPI_DOUBLE_PRECISION,COMM_global,info)
          call MPI_ALLGATHERV(yw,counts(iproc+1),MPI_DOUBLE_PRECISION,ywg,counts,disp,MPI_DOUBLE_PRECISION,COMM_global,info)
          call MPI_ALLGATHERV(zw,counts(iproc+1),MPI_DOUBLE_PRECISION,zwg,counts,disp,MPI_DOUBLE_PRECISION,COMM_global,info)
          
          ! Compute min distance to wall for every point
          ! --------------------------------------------
          allocate(d_(nx,ny,nz,size(xwg)),d(nx,ny,nz))
        
          do i=1,nx
             do j=1,ny
                do k=1,nz
                   do m=1,size(xwg)
                      d_(i,j,k,m) = sqrt((xc3(i,j,k)-xwg(m))**2 + (yc3(i,j,k)-ywg(m))**2 + &
                                         (zc3(i,j,k)-zwg(m))**2)
                   enddo
                   d(i,j,k) = minval(d_(i,j,k,:))
                enddo
             enddo
          enddo
       
       endif

       dg(:,:,:,ind)=d(:,:,:)
  
       ! Clear memory space
       ! ==================
       deallocate(xwg,ywg,zwg,xw,yw,zw,counts,disp,d)
       if (is_curv3) then
          deallocate(d_)
       else
          if (loc_ijk.ne.3) then 
             deallocate(dij)
          else
             deallocate(dk)
          endif
       endif
 
    endif
    
    deallocate(walls)    
  
  end subroutine compute_wall_dist

end module mod_wall_dist
