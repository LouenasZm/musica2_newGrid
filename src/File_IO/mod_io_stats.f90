!==============================================================================================
module mod_io_stats
!==============================================================================================
  !> author: Luca Sciacovelli, XG
  !> date: 2018, 2021
  !> Module for Input/Output routines
!==============================================================================================
  use mod_io_snapshots
  implicit none
  ! -------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------

contains

  !==============================================================================================
  subroutine read_write_stats_xy(operation,ibl2read)
  !==============================================================================================
    !> author: Luca Sciacovelli, XG
    !> date: June 2018, 2021
    !> read/write stats data on xy plane (averaged in time and along z)
  !==============================================================================================
    use mod_utils
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Input/Output arguments
    integer, intent(in) :: operation,ibl2read
    ! -------------------------------------------------------------------------------------------
    ! Local variables
    character(len=17) :: stats1file,stats2file
    integer  :: n,nx1p,nx2p,ny1p,ny2p,nz1p,nz2p
    integer  :: filetype
    ! -------------------------------------------------------------------------------------------

    ! Only if stats are computed (iteration counter greater than ndeb)
    ! ================================================================
    if (ntotal.lt.ndeb) return

    ! Define plane for storing stats quantities
    ! =========================================
    nx1p=plane_stats%tectype%nx1
    nx2p=plane_stats%tectype%nx2
    ny1p=plane_stats%tectype%ny1
    ny2p=plane_stats%tectype%ny2
    nz1p=plane_stats%tectype%nz1
    nz2p=plane_stats%tectype%nz2

    filetype=2

    ! I/O file names [cut in 2 parts]
    ! ==============
    if (plane_stats%tectype%is_IOtec_read) then
       stats1file='stats1_bl'//trim(numchar(ibl2read))//'.plt'
       stats2file='stats2_bl'//trim(numchar(ibl2read))//'.plt'
    else
       stats1file='stats1_bl'//trim(numchar(ibl2read))//'.bin'
       stats2file='stats2_bl'//trim(numchar(ibl2read))//'.bin'
    endif

    select case (operation)

    case (WRITE)

       ! Take mean along span (if z-direction is homogeneous)
       ! ====================
       avg_s=0.0_wp
       call MPI_REDUCE(avg_t,avg_s,nx*ny*nstat,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMMXY,info)
       avg_s=avg_s/dble(ndomz)

       ! First slot of statistics
       ! ========================
       ! pack data
       ndata=23
       allocate(varlist(ndata),dummy(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p,ndata))

       dummy(1:nx,1:ny,1,1:ndata)=avg_s(1:nx,1:ny,1:ndata)

       do n=1,ndata
          write(varlist(n)%name,'(A3,I0)') 'var',n
          varlist(n)%data => dummy(:,:,:,n)
       enddo

       if (iproc.eq.0) write(*,*) 'Writing stats1... '// stats1file
       call write_tec(trim(stats1file),filetype,tstar,plane_stats%tectype)
 
       deallocate(varlist,dummy)

       ! Second slot of statistics
       ! =========================
       ndata=nstat-23
       allocate(varlist(ndata),dummy(nx,ny,1,ndata))

       dummy(1:nx,1:ny,1,1:ndata)=avg_s(1:nx,1:ny,24:nstat)

       do n=1,ndata
          write(varlist(n)%name,'(A3,I0)') 'var',n+23
          varlist(n)%data => dummy(:,:,:,n)
       enddo

       if (iproc.eq.0) write(*,*) 'Writing stats2... '// stats2file
       call write_tec(trim(stats2file),filetype,tstar,plane_stats%tectype )

       deallocate(varlist,dummy)

    case (READ)

       ! First slot of statistics
       ! ========================
       
       ! check file existence
       call mpicheckfile(trim(dirDATA)//trim(stats1file))
       ndata=23
       allocate(varlist(ndata),dummy(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p,ndata))
       dummy=0.0_wp

       do n=1,ndata
          varlist(n)%data => dummy(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p,n)
       enddo

       if (iproc.eq.0) write(*,*) 'Reading stats1... '// stats1file
       call read_tec(trim(dirDATA)//trim(stats1file),plane_stats%tectype)
 
       ! ! If on Big-Endian machine (Turing), swap bytes to write in Little-Endian
       ! ! if (ichar(transfer(1,'a')) == 0) then
       !    ! Taille du type de base MPI_DOUBLE_PRECISION
       !    call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, nbytes_dbl, ierror)
       !    do n=1,ndata
       !       call conv( dummy(:,:,:,n), nbytes_dbl, size(dummy(:,:,:,n)) )
       !    enddo
       ! ! endif

       avg_t(1:nx,1:ny,1:ndata)=dummy(1:nx,1:ny,1,1:ndata)

       deallocate(varlist,dummy)
       
       ! Second slot of statistics
       ! =========================
       ! check file existence
       call mpicheckfile(trim(dirDATA)//trim(stats2file))
       ndata=nstat-23
       allocate(varlist(ndata),dummy(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p,ndata))
       dummy=0.0_wp

       do n=1,ndata
          varlist(n)%data => dummy(:,:,:,n)
       enddo

       if (iproc.eq.0) write(*,*) 'Reading stats2... '// stats2file
       !call read_tec('LES_IRS/'//trim(stats2file),plane_stats%tectype)
       call read_tec(trim(dirDATA)//trim(stats2file),plane_stats%tectype)

       ! ! If on Big-Endian machine (Turing), swap bytes to write in Little-Endian
       ! ! if (ichar(transfer(1,'a')) == 0) then
       !    ! Taille du type de base MPI_DOUBLE_PRECISION
       !    call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, nbytes_dbl, ierror)
       !    do n=1,ndata
       !       call conv( dummy(:,:,:,n), nbytes_dbl, size(dummy(:,:,:,n)) )
       !    enddo
       ! ! endif

       avg_t(1:nx,1:ny,24:nstat) = dummy(1:nx,1:ny,1,1:ndata)

       call MPI_BCAST(avg_t,nx*ny*nstat,MPI_DOUBLE_PRECISION,0,COMMXY,info)

       deallocate(varlist,dummy)
       
    end select

    call MPI_BARRIER(COMM_global,info)

  end subroutine read_write_stats_xy

  !==============================================================================================
  subroutine read_write_stats_xyz(operation)
  !==============================================================================================
    !> author: XG
    !> date: April 2023
    !> read/write stats data for 3D case [is_curv3 for example, only time-averaged]
  !==============================================================================================
    use mod_io_snapshots
    use mod_utils
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Input/Output arguments
    integer, intent(in) :: operation
    ! -------------------------------------------------------------------------------------------
    ! Local variables
    character(len=17) :: statsfile
    integer  :: n,nx1v,nx2v,ny1v,ny2v,nz1v,nz2v
    integer  :: filetype
    ! -------------------------------------------------------------------------------------------

    ! Only if stats are computed (iteration counter greater than ndeb)
    ! ================================================================
    if (ntotal.lt.ndeb) return

    ! Define volume for storing stats quantities
    ! ==========================================
    nx1v=volume_stats%tectype%nx1
    nx2v=volume_stats%tectype%nx2
    ny1v=volume_stats%tectype%ny1
    ny2v=volume_stats%tectype%ny2
    nz1v=volume_stats%tectype%nz1
    nz2v=volume_stats%tectype%nz2

    filetype=2

    ! I/O file names
    ! ==============
    if (volume_stats%tectype%is_IOtec_read) then
       statsfile='stats_bl'//trim(numchar(nob(iproc)))//'.plt'
    else
       statsfile='stats_bl'//trim(numchar(nob(iproc)))//'.bin'
    endif

    select case (operation)

    case (WRITE)

       ! Write statistics
       ! ================
       ! pack data
       ndata=nstat
       allocate(varlist(ndata),dummy(nx1v:nx2v,ny1v:ny2v,nz1v:nz2v,ndata))

       dummy(1:nx,1:ny,1:nz,1:ndata)=avg_v(1:nx,1:ny,1:nz,1:ndata)

       do n=1,ndata
          write(varlist(n)%name,'(A3,I0)') 'var',n
          varlist(n)%data => dummy(:,:,:,n)
       enddo

       if (iproc.eq.0) write(*,*) 'Writing stats... '// statsfile
       call write_tec(trim(statsfile),filetype,tstar,volume_stats%tectype)
 
       deallocate(varlist,dummy)

    case (READ)

       ! Check file existence
       ! ====================
       call mpicheckfile(trim(dirDATA)//trim(statsfile))
       
       ! Read statistics
       ! ===============
       ndata=nstat
       allocate(varlist(ndata),dummy(nx1v:nx2v,ny1v:ny2v,nz1v:nz2v,ndata))
       dummy=0.0_wp

       do n=1,ndata
          varlist(n)%data => dummy(nx1v:nx2v,ny1v:ny2v,nz1v:nz2v,n)
       enddo

       if (iproc.eq.0) write(*,*) 'Reading stats... '// statsfile
       call read_tec(trim(dirDATA)//trim(statsfile),volume_stats%tectype)

       avg_v(1:nx,1:ny,1:nz,1:ndata)=dummy(1:nx,1:ny,1:nz,1:ndata)

       deallocate(varlist,dummy)
       
    end select

    call MPI_BARRIER(COMM_global,info)

  end subroutine read_write_stats_xyz

  !==============================================================================================
  subroutine read_write_stats_xy_xyz(operation,ibl2read)
  !==============================================================================================
    !> author: XG
    !> date: April 2023
    !> read/write stats data for both planes andd volumes
    !> 1. read/write stats on  xy plane  (averaged in time and along z)
    !> 2. read/write stats on  3D volume (averaged in time)
    !> 3. read/write stats on wall plane (averaged in time)
  !==============================================================================================
    use mod_io_snapshots
    use mod_utils
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Input/Output arguments
    integer, intent(in) :: operation,ibl2read
    ! -------------------------------------------------------------------------------------------
    ! Local variables
    character(len=17) :: stats1file,stats2file,statsfile,statswfile
    integer  :: n,nx1v,nx2v,ny1v,ny2v,nz1v,nz2v
    integer  :: nx1p,nx2p,ny1p,ny2p,nz1p,nz2p
    integer  :: nx1w,nx2w,ny1w,ny2w,nz1w,nz2w
    integer  :: filetype
    ! -------------------------------------------------------------------------------------------

    ! Only if stats are computed (iteration counter greater than ndeb)
    ! ================================================================
    if (ntotal.lt.ndeb) return

    ! Define plane/volume for storing stats quantities
    ! ================================================
    filetype=2

    ! 1. xy plane (averaged in time and along z)
    ! ------------------------------------------
    nx1p=plane_stats%tectype%nx1
    nx2p=plane_stats%tectype%nx2
    ny1p=plane_stats%tectype%ny1
    ny2p=plane_stats%tectype%ny2
    nz1p=plane_stats%tectype%nz1
    nz2p=plane_stats%tectype%nz2

    ! I/O file names [cut in 2 parts]
    if (plane_stats%tectype%is_IOtec_read) then
       stats1file='stats1_bl'//trim(numchar(ibl2read))//'.plt'
       stats2file='stats2_bl'//trim(numchar(ibl2read))//'.plt'
    else
       stats1file='stats1_bl'//trim(numchar(ibl2read))//'.bin'
       stats2file='stats2_bl'//trim(numchar(ibl2read))//'.bin'
    endif
    
    ! 2. xyz volume (averaged in time)
    ! --------------------------------
    nx1v=volume_stats%tectype%nx1
    nx2v=volume_stats%tectype%nx2
    ny1v=volume_stats%tectype%ny1
    ny2v=volume_stats%tectype%ny2
    nz1v=volume_stats%tectype%nz1
    nz2v=volume_stats%tectype%nz2

    ! I/O file names
    if (volume_stats%tectype%is_IOtec_read) then
       statsfile='stats_bl'//trim(numchar(ibl2read))//'.plt'
    else
       statsfile='stats_bl'//trim(numchar(ibl2read))//'.bin'
    endif

    ! 3. xz wall plane (averaged in time)
    ! -----------------------------------
    ! NOTA only wall at jmin implemented
    if (is_bc_wall(2,1)) then
       nx1w=wall_stats%tectype%nx1
       nx2w=wall_stats%tectype%nx2
       ny1w=wall_stats%tectype%ny1
       ny2w=wall_stats%tectype%ny2
       nz1w=wall_stats%tectype%nz1
       nz2w=wall_stats%tectype%nz2

       ! I/O file names [cut in 2 parts]
       if (wall_stats%tectype%is_IOtec_read) then
          statswfile='wstats_bl'//trim(numchar(ibl2read))//'.plt'
       else
          statswfile='wstats_bl'//trim(numchar(ibl2read))//'.bin'
       endif
    endif
    
    !======================
    select case (operation)
    !======================

    ! Write statistics
    ! ================
    case (WRITE)

       ! 1. xy plane (averaged in time and along z)
       ! ------------------------------------------
       
       ! take mean along spanwise MPI domains [homogeneous direction is z]
       avg_s=0.0_wp
       call MPI_REDUCE(avg_t,avg_s,nx*ny*nstat,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMMXY,info)
       avg_s=avg_s/dble(ndomz)

       ! pack data for first slot of statistics
       ndata=23
       allocate(varlist(ndata),dummy(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p,ndata))

       dummy(1:nx,1:ny,1,1:ndata)=avg_s(1:nx,1:ny,1:ndata)

       do n=1,ndata
          write(varlist(n)%name,'(A3,I0)') 'var',n
          varlist(n)%data => dummy(:,:,:,n)
       enddo

       if (iproc.eq.0) write(*,*) 'Writing stats1 ... '// stats1file
       call write_tec(trim(stats1file),filetype,tstar,plane_stats%tectype)
 
       ! free memory
       deallocate(varlist,dummy)

       ! pack data for second slot of statistics
       ndata=nstat-23
       allocate(varlist(ndata),dummy(nx,ny,1,ndata))

       dummy(1:nx,1:ny,1,1:ndata)=avg_s(1:nx,1:ny,24:nstat)

       do n=1,ndata
          write(varlist(n)%name,'(A3,I0)') 'var',n+23
          varlist(n)%data => dummy(:,:,:,n)
       enddo

       if (iproc.eq.0) write(*,*) 'Writing stats2 ... '// stats2file
       call write_tec(trim(stats2file),filetype,tstar,plane_stats%tectype )

       ! free memory
       deallocate(varlist,dummy)
       
       ! 2. xyz volume (averaged in time)
       ! --------------------------------
       
       ! pack data for volumes
       ndata=nstat_v
       allocate(varlist(ndata),dummy(nx1v:nx2v,ny1v:ny2v,nz1v:nz2v,ndata))

       dummy(1:nx,1:ny,1:nz,1:ndata)=avg_v(1:nx,1:ny,1:nz,1:ndata)

       do n=1,ndata
          write(varlist(n)%name,'(A3,I0)') 'var',n
          varlist(n)%data => dummy(:,:,:,n)
       enddo

       if (iproc.eq.0) write(*,*) 'Writing stats ... '// statsfile
       call write_tec(trim(statsfile),filetype,tstar,volume_stats%tectype)
 
       ! free memory
       deallocate(varlist,dummy)
       
       ! 3. xz wall plane (averaged in time)
       ! -----------------------------------
       ! NOTA only wall at jmin implemented
       if (iproc.eq.0) write(*,*) 'Writing wall stats ...'
       
       if (is_bc_wall(2,1)) then
          ! pack data for volumes
          ndata=nstat_w
          allocate(varlist(ndata),dummy(nx1w:nx2w,ny1w:ny2w,nz1w:nz2w,ndata))

          dummy(1:nx,1,1:nz,1:ndata)=avg_w(1:nx,1:nz,1:ndata)

          do n=1,ndata
             write(varlist(n)%name,'(A3,I0)') 'var',n
             varlist(n)%data => dummy(:,:,:,n)
          enddo

          call write_tec(trim(statswfile),filetype,tstar,wall_stats%tectype)
          
          ! free memory
          deallocate(varlist,dummy)
       endif
    
    ! Read statistics
    ! ===============
    case (READ)

       ! 1. xy plane (averaged in time and along z)
       ! ------------------------------------------
       
       ! check file existence
       call mpicheckfile(trim(dirDATA)//trim(stats1file))
       call mpicheckfile(trim(dirDATA)//trim(stats2file))

       ! read operation for first slot of statistics
       ndata=23
       allocate(varlist(ndata),dummy(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p,ndata))
       dummy=0.0_wp

       do n=1,ndata
          varlist(n)%data => dummy(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p,n)
       enddo

       if (iproc.eq.0) write(*,*) 'Reading stats1 ... '// stats1file
       call read_tec(trim(dirDATA)//trim(stats1file),plane_stats%tectype)
 
       avg_t(1:nx,1:ny,1:ndata)=dummy(1:nx,1:ny,1,1:ndata)

       deallocate(varlist,dummy)
       
       ! read operation for second slot of statistics
       ndata=nstat-23
       allocate(varlist(ndata),dummy(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p,ndata))
       dummy=0.0_wp

       do n=1,ndata
          varlist(n)%data => dummy(:,:,:,n)
       enddo

       if (iproc.eq.0) write(*,*) 'Reading stats2 ... '// stats2file
       call read_tec(trim(dirDATA)//trim(stats2file),plane_stats%tectype)

       avg_t(1:nx,1:ny,24:nstat) = dummy(1:nx,1:ny,1,1:ndata)

       ! send to all procs in spanwise direction
       call MPI_BCAST(avg_t,nx*ny*nstat,MPI_DOUBLE_PRECISION,0,COMMXY,info)

       ! free memory
       deallocate(varlist,dummy)
       
       ! 2. xyz volume (averaged in time)
       ! --------------------------------
       
       ! check file existence
       call mpicheckfile(trim(dirDATA)//trim(statsfile))
       
       ! read operation
       ndata=nstat_v
       allocate(varlist(ndata),dummy(nx1v:nx2v,ny1v:ny2v,nz1v:nz2v,ndata))
       dummy=0.0_wp

       do n=1,ndata
          varlist(n)%data => dummy(nx1v:nx2v,ny1v:ny2v,nz1v:nz2v,n)
       enddo

       if (iproc.eq.0) write(*,*) 'Reading stats ... '// statsfile
       call read_tec(trim(dirDATA)//trim(statsfile),volume_stats%tectype)

       avg_v(1:nx,1:ny,1:nz,1:ndata)=dummy(1:nx,1:ny,1:nz,1:ndata)

       ! free memory
       deallocate(varlist,dummy)
       
       ! 3. xz wall plane (averaged in time)
       ! -----------------------------------
       ! NOTA only wall at jmin implemented
       if (iproc.eq.0) write(*,*) 'Reading wall stats ...'
       
       if (is_bc_wall(2,1)) then
          ! check file existence
          call mpicheckfile(trim(dirDATA)//trim(statswfile))

          ! read operation
          ndata=nstat_w
          allocate(varlist(ndata),dummy(nx1w:nx2w,ny1w:ny2w,nz1w:nz2w,ndata))
          dummy=0.0_wp

          do n=1,ndata
             varlist(n)%data => dummy(nx1w:nx2w,ny1w:ny2w,nz1w:nz2w,n)
          enddo

          call read_tec(trim(dirDATA)//trim(statswfile),wall_stats%tectype)

          avg_w(1:nx,1:nz,1:ndata)=dummy(1:nx,1,1:nz,1:ndata)

          ! free memory
          deallocate(varlist,dummy)
       endif
       
    end select

    call MPI_BARRIER(COMM_global,info)

  end subroutine read_write_stats_xy_xyz

  !==============================================================================================
  subroutine read_write_stats_chan(operation)
  !==============================================================================================
    !> author: Luca Sciacovelli
    !> date: April 2018
    !> read/write stats data for channel flow case
    !==============================================================================================
    use mod_io ! <- for WRITE,READ
    !use mod_mpi      ! <- for BCAST [included in mod_io]
    !use mod_flow     ! <- for avg_tg,avg_t [included in mod_io]
    !use mod_time     ! <- for ntotal [included in mod_io]
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Input/Output arguments
    integer, intent(in) :: operation
    ! -------------------------------------------------------------------------------------------
    ! Local variables
    integer :: i,j
    real(wp) :: bid
    ! -------------------------------------------------------------------------------------------

    ! Only if stats are computed (iteration counter greater than ndeb)
    ! ================================================================
    if (ntotal.lt.ndeb) return

    select case (operation)

    case (WRITE)

       ! Global space- and time-averaged quantitites
       ! ===========================================
       avg_tg=0.0_wp
       do j=1,ny
          avg_tg(1,j+coord(2)*ny,:) = avg_t(1,j,:)
       enddo

       call MPI_ALLREDUCE(MPI_IN_PLACE,avg_tg,size(avg_tg),MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
       avg_tg = avg_tg/dble(ndomx*ndomz)

       ! Write stats quantitites [on proc 0]
       ! =======================
       if (iproc.eq.0) then
          open(123,file='stats.dat',form='formatted',status='replace')
          do j=1,ngy
             write(123,'(181(1pE19.12,1X))') yg(j)/hc, (avg_tg(1,j,i),i=1,nstat)
          enddo
          close(123)
       endif

    case (READ)

       ! Initialize global space- and time-averaged quantitites
       ! ======================================================
       avg_tg=0.0_wp
       
       ! Read stats quantitites [on proc 0]
       ! ======================
       ! /!\ for channel flow 181 stats [correct if necessary]
       if (iproc.eq.0) then
          open(123,file='stats.dat',form='formatted',status='old')
          do j=1,ngy
             read(123,'(181(1pE19.12,1X))') bid, (avg_tg(1,j,i),i=1,nstat)
          enddo
          close(123)
       endif

       ! Local space- and time-averaged quantitites
       ! ==========================================
       call MPI_BCAST(avg_tg,size(avg_tg),MPI_DOUBLE_PRECISION,0,COMM_global,info)
       avg_t=0.0_wp
       do j=1,ny
          avg_t(1,j,:) = avg_tg(1,j+coord(2)*ny,:)
       enddo

    end select

    call MPI_BARRIER(COMM_global,info)

  end subroutine read_write_stats_chan

end module mod_io_stats
