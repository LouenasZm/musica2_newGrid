!==============================================================================
module mod_pp_extend
!==============================================================================
  !> Module to extend with 1 point the data
  !> Nota: not parallel -> to be run on a single proc per block
!==============================================================================
  ! use mod_pp_stats_var
  use mod_pp_var
  use warnstop
  use mod_io
  use mod_utils
  use mod_comm
  implicit none
  ! ---------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------

contains

  !============================================================================
  subroutine pp_extend_main
  !============================================================================
    !> main routine
  !============================================================================
    use mod_mpi
    use mod_constant
    use mod_block    ! for: block informations
    implicit none
    ! -------------------------------------------------------------------------
    integer :: ib_
    ! -------------------------------------------------------------------------

    ! Verification of 1 proc per block
    do ib_=1,nbloc
       if (bl(ib_)%nproc.gt.1) call mpistop('Only 1 proc per block can be used...',0)
    enddo

    ! Determination of block post-processed for each proc
    iblc_pp=0
    do ib_=1,nbloc_r
       if ((iproc.le.bl(ib_)%proc_max).and.(iproc.ge.bl(ib_)%proc_min)) then
          iblc_pp=nblr(ib_)
       endif
    enddo

    ! Type of extension
    select case(type_data)

    case(0) ! 0: point
       call mpistop('A point cannot be extended...',0)
    case(1) ! 1: line
       call mpistop('1 point extension not implemented yet for a line...',0)
    case(2) ! 2: plane
       call mpistop('1 point extension not implemented yet for a plane...',0)
    case(3) ! 3: volume
       call mpistop('1 point extension not implemented yet for a volume...',0)
    case(4) ! 4: stats

       call pp_extend_stats_2d
    case(5) ! 5: restart
       call pp_extend_restart
    case default
       call mpistop('bad choice of type_data for post-processing', 0)
    end select

  end subroutine pp_extend_main


  !============================================================================
  subroutine pp_extend_stats_2d
  !============================================================================
    !> 2D stats extension
  !============================================================================
    use mod_io_stats
    use mod_mpi
    use mod_constant
    implicit none
    ! -------------------------------------------------------------------------
    integer :: ivar,n,i,j,k,nvars
    integer :: strandID
    integer :: ni1,nj1,ni2,nj2
    character(len=5) :: zonename
    character(len=17) :: statsfile
    character(len=150), dimension(:), allocatable :: datanames
    real(wp), dimension(:,:,:), allocatable :: avg_t_ex
    real(wp), dimension(:,:,:,:), allocatable :: griddata
    real(wp), dimension(:,:), allocatable :: xc3m,yc3m
    ! -------------------------------------------------------------------------
    integer :: ierror,filetype
    integer :: fh      ! file handle
    real(kind=4) :: fzone,fhead
    character(len=8)  :: title
    character(len=50) :: bid
    real(wp), dimension(:), allocatable :: min_glob,max_glob
    integer, dimension(MPI_STATUS_SIZE) :: statut
    ! integer(kind=MPI_OFFSET_KIND) :: offset
    ! -------------------------------------------------------------------------------------------


    ! Initialization of extended communication
    ! ========================================
    call mpi_types_comm_ex(1)

    ! Read stats I/O file
    ! ===================
    call read_write_stats_xy(READ,iblc_pp)

    ! Creation of extended stats
    ! ==========================
    allocate(avg_t_ex(0:ngx+1,0:ngy+1,1:nstat))
    avg_t_ex(1:ngx,1:ngy,1:nstat) = avg_t(1:ngx,1:ngy,1:nstat)

    ! Communication of stats
    ! ======================
    do ivar=1,nstat
       call commun2d_ex(avg_t_ex(0:ngx+1,0:ngy+1,ivar))
    enddo

    ni1=1; ni2=ngx; nj1=1; nj2=ngy
    if (bl(nob(iproc))%BC(1)>0) ni1=0
    if (bl(nob(iproc))%BC(2)>0) ni2=ngx+1
    if (bl(nob(iproc))%BC(3)>0) nj1=0
    if (bl(nob(iproc))%BC(4)>0) nj2=ngy+1


    ! ===================
    ! Writting statistics
    ! ===================
    if (is_IOtec_write) then
       statsfile='stats_ex_bl'//trim(numchar(iblc_pp))//'.plt'
    else
       statsfile='stats_ex_bl'//trim(numchar(iblc_pp))//'.bin'
    endif

    if (iproc.eq.0) write(*,*) 'Writing stats... '//statsfile
    call MPI_FILE_OPEN(COMM_global,statsfile,MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL,fh,info)

    ! Tecplot header
    ! ==============
    if (is_IOtec_write) then
       ! Initialisation
       ! --------------
       ! ! offset initialisation
       ! offset=0_MPI_OFFSET_KIND
       ! datanames list definition
       nvars = nstat+3 ! var stats + X + Y + Z
       allocate(datanames(nvars))
       datanames(1) = "X"
       datanames(2) = "Y"
       datanames(3) = "Z"
       do n=1,nstat
          write(datanames(n+3),'(A3,I0)') 'var',n
       enddo
       ! delimiter
       fzone = 299.0
       fhead = 357.0
       title = 'mix'
       ! Tecplot needs to know minimum and maximum of each variable
       allocate(min_glob(nvars),max_glob(nvars))
       min_glob=0.0_wp
       max_glob=0.0_wp

       ! Grid
       allocate(griddata(ni1:ni2,nj1:nj2,1,3))
       if (is_curv3) then
          
          ! grid average
          allocate(xc3m(ni1:ni2,nj1:nj2),yc3m(ni1:ni2,nj1:nj2))
          xc3m=0.0_wp
          yc3m=0.0_wp
          do j=nj1,nj2
             do i=ni1,ni2
                do k=1,nz
                   xc3m(i,j)=xc3m(i,j)+xc3(i,j,k)
                   yc3m(i,j)=yc3m(i,j)+yc3(i,j,k)
                enddo              
             enddo
          enddo
          xc3m=xc3m/nz
          yc3m=yc3m/nz

          do j=nj1,nj2
             do i=ni1,ni2
                griddata(i,j,1,1) = xc3m(i,j)
                griddata(i,j,1,2) = yc3m(i,j)
                griddata(i,j,1,3) = 0.0_wp
             enddo
          enddo

          min_glob(1) = minval(xc3m(ni1:ni2,nj1:nj2))
          max_glob(1) = maxval(xc3m(ni1:ni2,nj1:nj2))
          min_glob(2) = minval(yc3m(ni1:ni2,nj1:nj2))
          max_glob(2) = maxval(yc3m(ni1:ni2,nj1:nj2))
          min_glob(3) = 0.0_wp
          max_glob(3) = 0.0_wp
       else
          do j=nj1,nj2
             do i=ni1,ni2
                griddata(i,j,1,1) = xgc(i,j)
                griddata(i,j,1,2) = ygc(i,j)
                griddata(i,j,1,3) = z(1)
             enddo
          enddo

          min_glob(1) = minval(xgc(ni1:ni2,nj1:nj2))
          max_glob(1) = maxval(xgc(ni1:ni2,nj1:nj2))
          min_glob(2) = minval(ygc(ni1:ni2,nj1:nj2))
          max_glob(2) = maxval(ygc(ni1:ni2,nj1:nj2))
          min_glob(3) = zg(1)
          max_glob(3) = zg(1)
       endif
       
       do n=1,nstat
          min_glob(n+3) = minval(avg_t_ex(ni1:ni2,nj1:nj2,n))
          max_glob(n+3) = maxval(avg_t_ex(ni1:ni2,nj1:nj2,n))
       enddo
       ! FileType (0:Full 1:Grid 2:Sol)
       filetype = 0
       ! zonename
       write(zonename,'(A5)') 'stats'
       ! strandID
       strandID = -2

       ! Writting header
       ! ---------------
       ! Magic Number
       call MPI_FILE_WRITE(fh,'#!TDV112',8,MPI_CHARACTER,statut,ierror)
       ! offset=offset+8
       ! Byte order of the reader
       call MPI_FILE_WRITE(fh,1,1,MPI_INTEGER,statut,ierror)
       ! offset=offset+4
       ! FileType (0:Full 1:Grid 2:Sol)
       call MPI_FILE_WRITE(fh,filetype,1,MPI_INTEGER,statut,ierror)
       ! offset=offset+4
       ! Title
       do i=1,len_trim(title)
          call MPI_FILE_WRITE(fh,ichar(title(i:i)),1,MPI_INTEGER,statut,ierror)
          ! offset=offset+4
       enddo
       ! Strings end with 0
       call MPI_FILE_WRITE(fh,0,1,MPI_INTEGER,statut,ierror)
       ! offset=offset+4
       ! Number of variables
       call MPI_FILE_WRITE(fh,nvars,1,MPI_INTEGER,statut,ierror)
       ! offset=offset+4
       ! Variable names
       do i=1,nvars
          bid=trim(datanames(i))
          ! Name of variable
          do j=1,len_trim(bid)
             call MPI_FILE_WRITE(fh,ichar(bid(j:j)),1,MPI_INTEGER,statut,ierror)
             ! offset=offset+4
          enddo
          ! Strings end with 0
          call MPI_FILE_WRITE(fh,0,1,MPI_INTEGER,statut,ierror)
          ! offset=offset+4
       enddo
       ! Zone Marker. Value = 299.0
       call MPI_FILE_WRITE(fh,fzone,1,MPI_REAL4,statut,ierror)
       ! offset=offset+4
       ! Zone name
       do i=1,len_trim(zonename)
          call MPI_FILE_WRITE(fh,ichar(zonename(i:i)),1,MPI_INTEGER,statut,ierror)
          ! offset=offset+4
       enddo
       ! Strings end with 0
       call MPI_FILE_WRITE(fh, 0,1,MPI_INTEGER,statut,ierror)
       ! offset=offset+4
       ! No parent Zone
       call MPI_FILE_WRITE(fh,-1,1,MPI_INTEGER,statut,ierror)
       ! offset=offset+4
       ! StrandID
       call MPI_FILE_WRITE(fh,field%strandID,1,MPI_INTEGER,statut,ierror)
       ! offset=offset+4
       ! Solution time
       call MPI_FILE_WRITE(fh,tstar,1,MPI_DOUBLE_PRECISION,statut,ierror)
       ! offset=offset+8
       ! Not used
       call MPI_FILE_WRITE(fh,-1,1,MPI_INTEGER,statut,ierror)
       ! offset=offset+4
       ! DataPacking (0:Block)
       call MPI_FILE_WRITE(fh, 0,1,MPI_INTEGER,statut,ierror)
       ! offset=offset+4
       ! Var location (0:Nodes)
       call MPI_FILE_WRITE(fh, 0,1,MPI_INTEGER,statut,ierror)
       ! offset=offset+4
       ! Are raw local 1-to-1 face neighbours supplied? (0: False)
       call MPI_FILE_WRITE(fh, 0,1,MPI_INTEGER,statut,ierror)
       ! offset=offset+4
       ! Number of miscellaneous user-defined face neighbor connections
       call MPI_FILE_WRITE(fh, 0,1,MPI_INTEGER,statut,ierror)
       ! offset=offset+4
       ! Imax
       call MPI_FILE_WRITE(fh,ni2-ni1+1,1,MPI_INTEGER,statut,ierror)
       ! offset=offset+4
       ! Jmax
       call MPI_FILE_WRITE(fh,nj2-nj1+1,1,MPI_INTEGER,statut,ierror)
       ! offset=offset+4
       ! Kmax
       call MPI_FILE_WRITE(fh,1,1,MPI_INTEGER,statut,ierror)
       ! offset=offset+4
       ! No more Auxiliary name/value pairs
       call MPI_FILE_WRITE(fh, 0,1,MPI_INTEGER,statut,ierror)
       ! offset=offset+4
       ! End of Header
       call MPI_FILE_WRITE(fh,fhead,1,MPI_REAL4,statut,ierror)
       ! offset=offset+4
       ! Data section
       call MPI_FILE_WRITE(fh,fzone,1,MPI_REAL4,statut,ierror)
       ! offset=offset+4
       ! Variable data format (2: double)
       do i=1,nvars
          call MPI_FILE_WRITE(fh,2,1,MPI_INTEGER,statut,ierror)
          ! offset=offset+4
       enddo
       ! Has passive variables (0: no)
       call MPI_FILE_WRITE(fh, 0,1,MPI_INTEGER,statut,ierror)
       ! offset=offset+4
       ! Has variable sharing (0: no)
       call MPI_FILE_WRITE(fh, 0,1,MPI_INTEGER,statut,ierror)
       ! offset=offset+4
       ! Share connectivity list with (-1: no)
       call MPI_FILE_WRITE(fh,-1,1,MPI_INTEGER,statut,ierror)
       ! offset=offset+4
       ! Min/Max of each variable
       do i=1,nvars
          call MPI_FILE_WRITE(fh,min_glob(i),1,MPI_DOUBLE_PRECISION,statut,ierror)
          call MPI_FILE_WRITE(fh,max_glob(i),1,MPI_DOUBLE_PRECISION,statut,ierror)
          ! offset=offset+2*8
       enddo

       deallocate(datanames)
       deallocate(min_glob,max_glob)
    endif

    ! Grid
    call MPI_FILE_WRITE(fh,griddata(:,:,:,1),size(griddata(:,:,:,1)),MPI_DOUBLE_PRECISION,statut,info)
    call MPI_FILE_WRITE(fh,griddata(:,:,:,2),size(griddata(:,:,:,2)),MPI_DOUBLE_PRECISION,statut,info)
    call MPI_FILE_WRITE(fh,griddata(:,:,:,3),size(griddata(:,:,:,3)),MPI_DOUBLE_PRECISION,statut,info)
    ! Statistics
    call MPI_FILE_WRITE(fh,avg_t_ex(ni1:ni2,nj1:nj2,1:nstat),size(avg_t_ex(ni1:ni2,nj1:nj2,1:nstat)),MPI_DOUBLE_PRECISION,statut,info)
    ! Closing file
    call MPI_FILE_CLOSE(fh,info)


  end subroutine pp_extend_stats_2d


  !============================================================================
  subroutine pp_extend_restart
  !============================================================================
    !> Restart extension
  !============================================================================
    use mod_io
    use mod_mpi
    use mod_interface
    ! use mod_constant
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i,j,k,nvars
    integer :: ni1,nj1,nk1,ni2,nj2,nk2
    character(len=20) :: restartfile
    ! -------------------------------------------------------------------------
    integer :: strandID
    character(len=5) :: zonename
    character(len=150), dimension(:), allocatable :: datanames
    real(wp), dimension(:,:,:,:), allocatable :: griddata
    integer :: ierror,filetype
    integer :: fh      ! file handle
    real(kind=4) :: fzone,fhead
    character(len=8)  :: title
    character(len=50) :: bid
    real(wp), dimension(:), allocatable :: min_glob,max_glob
    integer, dimension(MPI_STATUS_SIZE) :: statut
    ! -------------------------------------------------------------------------------------------

    ! Reading restart
    ! ===============
    allocate( rho(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(rhou(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(rhov(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(rhow(nx1:nx2,ny1:ny2,nz1:nz2))
    allocate(rhoe(nx1:nx2,ny1:ny2,nz1:nz2))
    if (is_RANS) allocate(nutil(nx1:nx2,ny1:ny2,nz1:nz2))
    ! if (trim(filestamp).eq.'0000_0000') then
       binfile= 'restart_bl'//trim(numchar(nob(iproc)))//filext_read
       TDfile = 'restartTD.bin'
    ! else
    !    binfile= 'restart_bl'//trim(numchar(nob(iproc)))//filext_read
    !    TDfile = 'restartTD'//filestamp//'.bin'
    ! endif

    call read_write_info(READ)
    call read_write_volume(binfile,TDfile,READ)

    ! Fill ghost cells
    ! ================
    call communication_(rho,rhou,rhov,rhow,rhoe)
    if (is_RANS) call communication_rans(nutil)

    ! Determination of borns
    ! ======================
    ni1=1; ni2=ngx; nj1=1; nj2=ngy; nk1=1; nk2=ngz
    if (bl(nob(iproc))%BC(1)>0) ni1=0
    if (bl(nob(iproc))%BC(2)>0) ni2=ngx+1
    if (bl(nob(iproc))%BC(3)>0) nj1=0
    if (bl(nob(iproc))%BC(4)>0) nj2=ngy+1
    if (bl(nob(iproc))%BC(5)>0) nk1=0
    if (bl(nob(iproc))%BC(6)>0) nk2=ngz+1

    ! ===================
    ! Writting restart_ex
    ! ===================
    if (is_IOtec_write) then
       restartfile='restart_ex_bl'//trim(numchar(iblc_pp))//'.plt'
    else
       restartfile='restart_ex_bl'//trim(numchar(iblc_pp))//'.bin'
    endif

    print *,'iproc',iproc,'iblc_pp',iblc_pp,restartfile

    !call mpistop('',0)

    if (iproc.eq.0) write(*,*) 'Writing restart... '//restartfile
    call MPI_FILE_OPEN(COMM_global,restartfile,MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL,fh,info)
    !call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(restartfile),MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL,fh,info)


!!$    open(194,file='test_bl'//trim(numchar(nob(iproc)))//'.bin',form='unformatted',status='unknown')
!!$    rewind(194)
!!$    write(194) ngx
!!$    write(194) ngy
!!$    write(194) ngz
!!$    close(194)
!!$
!!$    call mpistop('',0)


    ! Tecplot header
    ! ==============
    if (is_IOtec_write) then
       ! Initialisation
       ! --------------
       ! datanames list definition
       nvars = 8 ! X Y Z rho rhou rhov rhow rhoe
       if (is_RANS) nvars = 9 ! X Y Z rho rhou rhov rhow rhoe nutilde
       allocate(datanames(nvars))
       datanames(1) = "X"
       datanames(2) = "Y"
       datanames(3) = "Z"
       datanames(4) = "rho"
       datanames(5) = "rhou"
       datanames(6) = "rhov"
       datanames(7) = "rhow"
       datanames(8) = "rhoe"
       if (is_RANS) datanames(9) = "nutilde"
       ! delimiter
       fzone = 299.0
       fhead = 357.0
       title = 'mix'
       ! Tecplot needs to know minimum and maximum of each variable
       allocate(min_glob(nvars),max_glob(nvars))
       min_glob=0.0_wp
       max_glob=0.0_wp

       ! Grid
       allocate(griddata(ni1:ni2,nj1:nj2,nk1:nk2,3))
       if (is_curv3) then
          do k=nk1,nk2
             do j=nj1,nj2
                do i=ni1,ni2
                   griddata(i,j,k,1) = xgc3(i,j,k)
                   griddata(i,j,k,2) = ygc3(i,j,k)
                   griddata(i,j,k,3) = zgc3(i,j,k)
                enddo
             enddo
          enddo
       else
          do k=nk1,nk2
             do j=nj1,nj2
                do i=ni1,ni2
                   griddata(i,j,k,1) = xgc(i,j)
                   griddata(i,j,k,2) = ygc(i,j)
                   griddata(i,j,k,3) = z(k)
                enddo
             enddo
          enddo
       endif
    
       if (is_curv3) then
          min_glob(1) = minval(xgc3(ni1:ni2,nj1:nj2,nk1:nk2))
          max_glob(1) = maxval(xgc3(ni1:ni2,nj1:nj2,nk1:nk2))
          min_glob(2) = minval(ygc3(ni1:ni2,nj1:nj2,nk1:nk2))
          max_glob(2) = maxval(ygc3(ni1:ni2,nj1:nj2,nk1:nk2))
          min_glob(3) = minval(zgc3(ni1:ni2,nj1:nj2,nk1:nk2))
          max_glob(3) = maxval(zgc3(ni1:ni2,nj1:nj2,nk1:nk2))
       else
          min_glob(1) = minval(xgc(ni1:ni2,nj1:nj2)); max_glob(1) = maxval(xgc(ni1:ni2,nj1:nj2))
          min_glob(2) = minval(ygc(ni1:ni2,nj1:nj2)); max_glob(2) = maxval(ygc(ni1:ni2,nj1:nj2))
          min_glob(3) = minval(zg(nk1:nk2)); max_glob(3) = maxval(zg(nk1:nk2))
       endif
       min_glob(4) = minval(rho(ni1:ni2,nj1:nj2,nk1:nk2));  max_glob(4) = maxval(rho(ni1:ni2,nj1:nj2,nk1:nk2))
       min_glob(5) = minval(rhou(ni1:ni2,nj1:nj2,nk1:nk2)); max_glob(5) = maxval(rhou(ni1:ni2,nj1:nj2,nk1:nk2))
       min_glob(6) = minval(rhov(ni1:ni2,nj1:nj2,nk1:nk2)); max_glob(6) = maxval(rhov(ni1:ni2,nj1:nj2,nk1:nk2))
       min_glob(7) = minval(rhow(ni1:ni2,nj1:nj2,nk1:nk2)); max_glob(7) = maxval(rhow(ni1:ni2,nj1:nj2,nk1:nk2))
       min_glob(8) = minval(rhoe(ni1:ni2,nj1:nj2,nk1:nk2)); max_glob(8) = maxval(rhoe(ni1:ni2,nj1:nj2,nk1:nk2))
       if (is_RANS) then
          min_glob(9) = minval(nutil(ni1:ni2,nj1:nj2,nk1:nk2))
          max_glob(9) = maxval(nutil(ni1:ni2,nj1:nj2,nk1:nk2))
       endif

       ! FileType (0:Full 1:Grid 2:Sol)
       filetype = 0
       ! zonename
       write(zonename,'(A5)') 'stats'
       ! strandID
       strandID = -2

       ! Writting header
       ! ---------------
       ! Magic Number
       call MPI_FILE_WRITE(fh,'#!TDV112',8,MPI_CHARACTER,statut,ierror)
       ! Byte order of the reader
       call MPI_FILE_WRITE(fh,1,1,MPI_INTEGER,statut,ierror)
       ! FileType (0:Full 1:Grid 2:Sol)
       call MPI_FILE_WRITE(fh,filetype,1,MPI_INTEGER,statut,ierror)
       ! Title
       do i=1,len_trim(title)
          call MPI_FILE_WRITE(fh,ichar(title(i:i)),1,MPI_INTEGER,statut,ierror)
       enddo
       ! Strings end with 0
       call MPI_FILE_WRITE(fh,0,1,MPI_INTEGER,statut,ierror)
       ! Number of variables
       call MPI_FILE_WRITE(fh,nvars,1,MPI_INTEGER,statut,ierror)
       ! Variable names
       do i=1,nvars
          bid=trim(datanames(i))
          ! Name of variable
          do j=1,len_trim(bid)
             call MPI_FILE_WRITE(fh,ichar(bid(j:j)),1,MPI_INTEGER,statut,ierror)
          enddo
          ! Strings end with 0
          call MPI_FILE_WRITE(fh,0,1,MPI_INTEGER,statut,ierror)
       enddo
       ! Zone Marker. Value = 299.0
       call MPI_FILE_WRITE(fh,fzone,1,MPI_REAL4,statut,ierror)
       ! Zone name
       do i=1,len_trim(zonename)
          call MPI_FILE_WRITE(fh,ichar(zonename(i:i)),1,MPI_INTEGER,statut,ierror)
       enddo
       ! Strings end with 0
       call MPI_FILE_WRITE(fh, 0,1,MPI_INTEGER,statut,ierror)
       ! No parent Zone
       call MPI_FILE_WRITE(fh,-1,1,MPI_INTEGER,statut,ierror)
       ! StrandID
       call MPI_FILE_WRITE(fh,field%strandID,1,MPI_INTEGER,statut,ierror)
       ! Solution time
       call MPI_FILE_WRITE(fh,tstar,1,MPI_DOUBLE_PRECISION,statut,ierror)
       ! Not used
       call MPI_FILE_WRITE(fh,-1,1,MPI_INTEGER,statut,ierror)
       ! DataPacking (0:Block)
       call MPI_FILE_WRITE(fh, 0,1,MPI_INTEGER,statut,ierror)
       ! Var location (0:Nodes)
       call MPI_FILE_WRITE(fh, 0,1,MPI_INTEGER,statut,ierror)
       ! Are raw local 1-to-1 face neighbours supplied? (0: False)
       call MPI_FILE_WRITE(fh, 0,1,MPI_INTEGER,statut,ierror)
       ! Number of miscellaneous user-defined face neighbor connections
       call MPI_FILE_WRITE(fh, 0,1,MPI_INTEGER,statut,ierror)
       ! Imax
       call MPI_FILE_WRITE(fh,ni2-ni1+1,1,MPI_INTEGER,statut,ierror)
       ! Jmax
       call MPI_FILE_WRITE(fh,nj2-nj1+1,1,MPI_INTEGER,statut,ierror)
       ! Kmax
       call MPI_FILE_WRITE(fh,nk2-nk1+1,1,MPI_INTEGER,statut,ierror)
       ! No more Auxiliary name/value pairs
       call MPI_FILE_WRITE(fh, 0,1,MPI_INTEGER,statut,ierror)
       ! End of Header
       call MPI_FILE_WRITE(fh,fhead,1,MPI_REAL4,statut,ierror)
       ! Data section
       call MPI_FILE_WRITE(fh,fzone,1,MPI_REAL4,statut,ierror)
       ! Variable data format (2: double)
       do i=1,nvars
          call MPI_FILE_WRITE(fh,2,1,MPI_INTEGER,statut,ierror)
       enddo
       ! Has passive variables (0: no)
       call MPI_FILE_WRITE(fh, 0,1,MPI_INTEGER,statut,ierror)
       ! Has variable sharing (0: no)
       call MPI_FILE_WRITE(fh, 0,1,MPI_INTEGER,statut,ierror)
       ! Share connectivity list with (-1: no)
       call MPI_FILE_WRITE(fh,-1,1,MPI_INTEGER,statut,ierror)
       ! Min/Max of each variable
       do i=1,nvars
          call MPI_FILE_WRITE(fh,min_glob(i),1,MPI_DOUBLE_PRECISION,statut,ierror)
          call MPI_FILE_WRITE(fh,max_glob(i),1,MPI_DOUBLE_PRECISION,statut,ierror)
       enddo

       deallocate(datanames)
       deallocate(min_glob,max_glob)
    endif

    ! Grid
    call MPI_FILE_WRITE(fh,griddata(:,:,:,1),size(griddata(:,:,:,1)),MPI_DOUBLE_PRECISION,statut,info)
    call MPI_FILE_WRITE(fh,griddata(:,:,:,2),size(griddata(:,:,:,2)),MPI_DOUBLE_PRECISION,statut,info)
    call MPI_FILE_WRITE(fh,griddata(:,:,:,3),size(griddata(:,:,:,3)),MPI_DOUBLE_PRECISION,statut,info)
    call MPI_FILE_WRITE(fh, rho(ni1:ni2,nj1:nj2,nk1:nk2),size( rho(ni1:ni2,nj1:nj2,nk1:nk2)),MPI_DOUBLE_PRECISION,statut,info)
    call MPI_FILE_WRITE(fh,rhou(ni1:ni2,nj1:nj2,nk1:nk2),size(rhou(ni1:ni2,nj1:nj2,nk1:nk2)),MPI_DOUBLE_PRECISION,statut,info)
    call MPI_FILE_WRITE(fh,rhov(ni1:ni2,nj1:nj2,nk1:nk2),size(rhov(ni1:ni2,nj1:nj2,nk1:nk2)),MPI_DOUBLE_PRECISION,statut,info)
    call MPI_FILE_WRITE(fh,rhow(ni1:ni2,nj1:nj2,nk1:nk2),size(rhow(ni1:ni2,nj1:nj2,nk1:nk2)),MPI_DOUBLE_PRECISION,statut,info)
    call MPI_FILE_WRITE(fh,rhoe(ni1:ni2,nj1:nj2,nk1:nk2),size(rhoe(ni1:ni2,nj1:nj2,nk1:nk2)),MPI_DOUBLE_PRECISION,statut,info)
    if (is_RANS) call MPI_FILE_WRITE(fh,nutil(ni1:ni2,nj1:nj2,nk1:nk2),size(nutil(ni1:ni2,nj1:nj2,nk1:nk2)),MPI_DOUBLE_PRECISION,statut,info)
    ! Closing file
    call MPI_FILE_CLOSE(fh,info)


  end subroutine pp_extend_restart

end module mod_pp_extend
