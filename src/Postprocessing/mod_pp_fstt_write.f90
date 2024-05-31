!==============================================================================
module mod_pp_fstt_write
!==============================================================================
  !> Module for writting outputs of fstt pp
!==============================================================================
  use warnstop
  use mod_constant
  use mod_flow
  use mod_mpi
  use mod_io_snapshots
  use mod_pp_var
  use mod_pp_fstt_comm
  implicit none
  ! ---------------------------------------------------------------------------
  type(snapshot_type) :: planes_stats_pp_fstt
  type(snapshot_type) :: volume_check_fstt
  type(snapshot_type), dimension(:), allocatable :: snap_pp_fstt
  ! -------------------------------------------------------------------------------------------

contains

!==============================================================================================
  subroutine read_write_stats_pp_fstt(operation)
  !==============================================================================================
    !> author: AB
    !> date: June 2022
    !> write discriminated stats data for pp fstt
  !==============================================================================================
    use mod_utils
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Input/Output arguments
    integer, intent(in) :: operation
    ! -------------------------------------------------------------------------------------------
    ! Local variables
    character(len=20) :: stats1file,stats2file,stats3file
    integer  :: n,nx1p,nx2p,ny1p,ny2p,nz1p,nz2p,nstat_pp
    integer  :: i,j,l
    integer  :: filetype
    ! -------------------------------------------------------------------------------------------

    nstat_pp = 10

    nx1p=planes_stats_pp_fstt%tectype%nx1
    nx2p=planes_stats_pp_fstt%tectype%nx2
    ny1p=planes_stats_pp_fstt%tectype%ny1
    ny2p=planes_stats_pp_fstt%tectype%ny2
    nz1p=planes_stats_pp_fstt%tectype%nz1
    nz2p=planes_stats_pp_fstt%tectype%nz2

    filetype=2

    ! I/O file names
    ! ==============
    if (planes_stats_pp_fstt%tectype%is_IOtec_read) then
       stats1file='stats_lam_bl'//trim(numchar(iblc_pp))//'.plt'
       stats2file='stats_turb_bl'//trim(numchar(iblc_pp))//'.plt'
       stats3file='stats_tot_bl'//trim(numchar(iblc_pp))//'.plt'
    else
       stats1file='stats_lam_bl'//trim(numchar(iblc_pp))//'.bin'
       stats2file='stats_turb_bl'//trim(numchar(iblc_pp))//'.bin'
       stats3file='stats_tot_bl'//trim(numchar(iblc_pp))//'.bin'
    endif

    select case (operation)

    case (WRITE)
       ! Take mean along span (z-direction)
       ! ====================
       ! ! if (iproc.eq.0) stats_lam(10,10:80,1) = 0.0_wp
       ! if (iproc.eq.0) print *,'stats_lam 0 before',stats_lam(10,59,1)
       ! ! if (iproc.eq.1) stats_lam(10,10:80,1) = 0.2_wp
       ! if (iproc.eq.1) print *,'stats_lam 1 before',stats_lam(10,59,1)

       ! Laminar statistics
       ! ------------------
       avg_s_lam=0.0_wp
       ! Ponderation of proc contribution
       do l=2,10
          do j=1,ny
             do i=1,nx
                stats_lam(i,j,l)  = stats_lam(i,j,l)*stats_lam(i,j,1)
             enddo
          enddo
       enddo
       call MPI_REDUCE(stats_lam(1:nx,1:ny,1:nstat_pp),avg_s_lam(1:nx,1:ny,1:nstat_pp),nx*ny*nstat_pp,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMMXY,info)
       ! if (iproc.eq.0) print *,'avg_s_lam 0 after sum',avg_s_lam(10,59,1)
       ! Normalisation by the sum of ponderation
       do l=2,10
          do j=1,ny
             do i=1,nx
                if (avg_s_lam(i,j,1).gt.0.0_wp) avg_s_lam(i,j,l)  = avg_s_lam(i,j,l)/avg_s_lam(i,j,1)
             enddo
          enddo
       enddo
       ! Mean of the intermittency
       avg_s_lam(:,:,1) = avg_s_lam(:,:,1)/dble(ndomz)
       ! if (iproc.eq.0) print *,'avg_s_lam 0 after mean',avg_s_lam(10,59,1)
       ! call mpistop('',0)
       ! ! Old
       ! call MPI_REDUCE(stats_lam,avg_s_lam,nx*ny*nstat_pp,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMMXY,info)
       ! avg_s_lam=avg_s_lam/dble(ndomz)

       ! Turbulent statistics
       ! --------------------
       avg_s_turb=0.0_wp
       ! Ponderation of proc contribution
       do l=2,10
          do j=1,ny
             do i=1,nx
                stats_turb(i,j,l)  = stats_turb(i,j,l)*stats_turb(i,j,1)
             enddo
          enddo
       enddo
       call MPI_REDUCE(stats_turb,avg_s_turb,nx*ny*nstat_pp,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMMXY,info)
       ! Normalisation by the sum of ponderation coeff.
       do l=2,10
          do j=1,ny
             do i=1,nx
                if (avg_s_turb(i,j,1).gt.0.0_wp) avg_s_turb(i,j,l)  = avg_s_turb(i,j,l)/avg_s_turb(i,j,1)
             enddo
          enddo
       enddo
       avg_s_turb(:,:,1) = avg_s_turb(:,:,1)/dble(ndomz)
       ! ! Old
       ! call MPI_REDUCE(stats_turb,avg_s_turb,nx*ny*nstat_pp,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMMXY,info)
       ! avg_s_turb=avg_s_turb/dble(ndomz)

       ! Global statistics (if tested)
       ! -----------------------------
       if (is_check_stats) then
          avg_s_tot=0.0_wp
          call MPI_REDUCE(stats_tot,avg_s_tot,nx*ny*nstat_pp,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMMXY,info)
          avg_s_tot=avg_s_tot/dble(ndomz)
       endif

       ! Laminar statistics
       ! ==================
       ndata=nstat_pp
       allocate(varlist(ndata),dummy(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p,ndata))

       ! dummy(1:nx,1:ny,1,1:ndata)=avg_s_lam(1:nx,1:ny,1:ndata)
       do n=1,ndata
          dummy(1:nx,1:ny,1,n)=avg_s_lam(1:nx,1:ny,n)
       enddo

       do n=1,ndata
          write(varlist(n)%name,'(A3,I0)') 'var',n
          varlist(n)%data => dummy(:,:,:,n)
       enddo

       if (iproc.eq.0) write(*,*) 'Writing laminar stats... '// stats1file
       call write_tec(trim(stats1file),filetype,tstar,planes_stats_pp_fstt%tectype)

       deallocate(varlist,dummy)
       call MPI_BARRIER(COMM_global,info)

       ! Turbulent statistics
       ! ====================
       ndata=nstat_pp
       allocate(varlist(ndata),dummy(nx,ny,1,ndata))

       dummy(1:nx,1:ny,1,1:ndata)=avg_s_turb(1:nx,1:ny,1:ndata)

       do n=1,ndata
          write(varlist(n)%name,'(A3,I0)') 'var',n
          varlist(n)%data => dummy(:,:,:,n)
       enddo

       if (iproc.eq.0) write(*,*) 'Writing turbulent stats... '// stats2file
       call write_tec(trim(stats2file),filetype,tstar,planes_stats_pp_fstt%tectype)

       deallocate(varlist,dummy)

       ! Global statistics (if tested)
       ! =============================
       if (is_check_stats) then
          ndata=nstat_pp
          allocate(varlist(ndata),dummy(nx,ny,1,ndata))

          dummy(1:nx,1:ny,1,1:ndata)=avg_s_tot(1:nx,1:ny,1:ndata)

          do n=1,ndata
             write(varlist(n)%name,'(A3,I0)') 'var',n
             varlist(n)%data => dummy(:,:,:,n)
          enddo

          if (iproc.eq.0) write(*,*) 'Writing global stats to check-up... '// stats3file
          call write_tec(trim(stats3file),filetype,tstar,planes_stats_pp_fstt%tectype)

          deallocate(varlist,dummy)
       endif

    case (READ)
       ! Check files existence
       ! =====================
       call mpicheckfile(stats1file)

       ndata=nstat_pp

       ! Laminar statistics
       ! ==================
       allocate(varlist(ndata),dummy(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p,ndata))
       dummy=0.0_wp

       do n=1,ndata
          write(varlist(n)%name,'(A3,I0)') 'var',n
          varlist(n)%data => dummy(:,:,:,n)
       enddo

       if (iproc.eq.0) write(*,*) 'Reading laminar stats... '// stats1file
       call read_tec(trim(stats1file),planes_stats_pp_fstt%tectype)

       stats_lam(1:nx,1:ny,1:ndata)=dummy(1:nx,1:ny,1,1:ndata)
       call MPI_BCAST(stats_lam,nx*ny*ndata,MPI_DOUBLE_PRECISION,0,COMMXY,info)
       deallocate(varlist,dummy)

       ! Turbulent statistics
       ! ====================
       allocate(varlist(ndata),dummy(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p,ndata))
       dummy=0.0_wp

       do n=1,ndata
          write(varlist(n)%name,'(A3,I0)') 'var',n
          varlist(n)%data => dummy(:,:,:,n)
       enddo

       if (iproc.eq.0) write(*,*) 'Reading turbulent stats... '// stats2file
       call read_tec(trim(stats2file),planes_stats_pp_fstt%tectype)

       stats_turb(1:nx,1:ny,1:ndata)=dummy(1:nx,1:ny,1,1:ndata)
       call MPI_BCAST(stats_turb,nx*ny*ndata,MPI_DOUBLE_PRECISION,0,COMMXY,info)
       deallocate(varlist,dummy)

       ! Global statistics
       ! =================
       if (is_check_stats) then
          allocate(varlist(ndata),dummy(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p,ndata))
          dummy=0.0_wp

          do n=1,ndata
             write(varlist(n)%name,'(A3,I0)') 'var',n
             varlist(n)%data => dummy(:,:,:,n)
          enddo

          if (iproc.eq.0) write(*,*) 'Reading global stats for check-up... '// stats3file
          call read_tec(trim(stats3file),planes_stats_pp_fstt%tectype)

          stats_tot(1:nx,1:ny,1:ndata)=dummy(1:nx,1:ny,1,1:ndata)
          call MPI_BCAST(stats_tot,nx*ny*ndata,MPI_DOUBLE_PRECISION,0,COMMXY,info)
          deallocate(varlist,dummy)
       endif

       ! if (iproc.eq.0) print *,'stats_lam 1 before',stats_lam(10,10:80,1)

    end select

    call MPI_BARRIER(COMM_global,info)

  end subroutine read_write_stats_pp_fstt

  !==============================================================================================
  subroutine init_volume_check(volume)
  !==============================================================================================
    !> author: AB
    !> date: May 2022
    !> Init volume to check
  !==============================================================================================
    use mod_utils
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Input/Output: plane type
    type(snapshot_type), intent(inout) :: volume
    ! -------------------------------------------------------------------------------------------
    logical :: iexist
    integer :: nbl
    integer :: i,volnprc
    integer :: GROUP_intrablock,GROUP_volume
    integer :: ivolume,iproc_has_volume
    integer, dimension(:), allocatable :: iproc_has_volume_glo,volume_ranks
    ! -------------------------------------------------------------------------------------------

    ! Verification of file existence
    ! ==============================
    if ((coord(1)==0).and.(coord(2)==0).and.(coord(3)==0)) then
       inquire( file='./check_volume_bl'//trim(numchar(iblc_pp))//'.bin', exist=iexist )
       if (iexist) then
          call system('rm check_volume_bl'//trim(numchar(iblc_pp))//'.bin')
       endif
    endif

    ! Number of current block
    ! =======================
    nbl=nob(iproc)

    ! Define number of points and number of processors
    ! ================================================
    ! Set the correct indexes
    ! -----------------------
    volume%tectype%ngx = volume%ind_i2 - volume%ind_i1 + 1
    volume%tectype%ngy = volume%ind_j2 - volume%ind_j1 + 1
    volume%tectype%ngz = volume%ind_k2 - volume%ind_k1 + 1

    ! Create volume communicators
    ! ===========================

    ! Initialize proc list for the volumes
    ! ------------------------------------
    iproc_has_volume=-1

    ! Determine procs belonging to a given volume
    ! -------------------------------------------
    if (((coord(1)*nil_interp.lt.volume%ind_i1.and.(coord(1)+1)*nil_interp.ge.volume%ind_i1) .or. &
         (coord(1)*nil_interp.lt.volume%ind_i2.and.(coord(1)+1)*nil_interp.ge.volume%ind_i2) .or. &
         (coord(1)*nil_interp.gt.volume%ind_i1.and. coord(1)   *nil_interp.lt.volume%ind_i2)).and. &
        ((coord(2)* nj_interp.lt.volume%ind_j1.and.(coord(2)+1)* nj_interp.ge.volume%ind_j1) .or. &
         (coord(2)* nj_interp.lt.volume%ind_j2.and.(coord(2)+1)* nj_interp.ge.volume%ind_j2) .or.&
         (coord(2)* nj_interp.gt.volume%ind_j1.and. coord(2)   * nj_interp.lt.volume%ind_j2)).and.&
        ((coord(3)* nk_interp.lt.volume%ind_k1.and.(coord(3)+1)* nk_interp.ge.volume%ind_k1) .or. &
         (coord(3)* nk_interp.lt.volume%ind_k1.and.(coord(3)+1)* nk_interp.ge.volume%ind_k1) .or. &
         (coord(3)* nk_interp.gt.volume%ind_k1.and. coord(3)   * nk_interp.lt.volume%ind_k2))) then
         iproc_has_volume = 1
    endif

    ! Communicate full list
    ! ---------------------
    allocate(iproc_has_volume_glo(bl(nbl)%nproc))
    iproc_has_volume_glo= -9999
    call MPI_ALLGATHER(iproc_has_volume,1,MPI_INTEGER,iproc_has_volume_glo, &
                       1,MPI_INTEGER,COMM_intrablock,info)

    ! Allocate volume list
    ! --------------------
    volnprc = 0
    do i=1,bl(nbl)%nproc
       if (iproc_has_volume_glo(i).gt.0) then
          volnprc=volnprc+1
       endif
    enddo
    allocate(volume_ranks(volnprc))

    ! Fill volume list
    ! ----------------
    ivolume= 1
    do i=1,bl(nbl)%nproc
       if (iproc_has_volume_glo(i).gt.0) then
          volume_ranks(ivolume)= i-1
          ivolume=ivolume+1
       endif
    enddo

    ! Create sub-communicator
    ! -----------------------
    ! Get world group
    call MPI_COMM_GROUP(COMM_intrablock,GROUP_intrablock,info)
    ! Create plane group
    call MPI_GROUP_INCL(GROUP_intrablock,volnprc,volume_ranks,GROUP_volume,info)
    ! Destroy world group
    call MPI_GROUP_FREE(GROUP_intrablock,info)
    ! Create communicator
    call MPI_COMM_CREATE(COMM_intrablock,GROUP_volume,volume%tectype%MPI_COMM,info)
    ! Destroy plane group
    call MPI_GROUP_FREE(GROUP_volume,info)

    ! free allocated memory
    deallocate(volume_ranks,iproc_has_volume_glo)

    ! Initialize plane local indices
    ! ==============================

    ! default initializations
    ! -----------------------
    volume%tectype%nx1 = 1;  volume%tectype%nx2 = nil_interp
    volume%tectype%ny1 = 1;  volume%tectype%ny2 = nj_interp
    volume%tectype%nz1 = 1;  volume%tectype%nz2 = nk_interp

    ! Initialize I/O type
    ! ===================
    if (iproc_has_volume.eq.1) then
       call mod_io_init(volume%tectype%ngx,volume%tectype%ngy,volume%tectype%ngz, &
                        volume%tectype%nx2-volume%tectype%nx1+1, &
                        volume%tectype%ny2-volume%tectype%ny1+1, &
                        volume%tectype%nz2-volume%tectype%nz1+1, &
                        volume%tectype%nx2-volume%tectype%nx1+1, &
                        volume%tectype%ny2-volume%tectype%ny1+1, &
                        volume%tectype%nz2-volume%tectype%nz1+1, &
                        ngh_pp,3,coord,is_IOtec_read,is_IOtec_write,volume%tectype)
    endif

    ! logical is_app (is append_file)
    volume%tectype%is_app=.true.
    ! define offset: number of points * number of variables * 8 bytes
    volume%tectype%disp=volume%tectype%ngx &
                       *volume%tectype%ngy &
                       *volume%tectype%ngz*8*volume%nvar
    ! inform about restart mode
    volume%tectype%restart=.false.

  end subroutine init_volume_check

  !==============================================================================================
  subroutine write_volume_check
  !==============================================================================================
    !> author: AB
    !> date: May 2022
    !> Write volume to check values on interpolated grid
  !==============================================================================================
    use mod_utils
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Local variables
    integer :: i,j,k,m,filetype,ndata,cpt
    integer :: nx1l,nx2l,ny1l,ny2l,nz1l,nz2l
    ! -------------------------------------------------------------------------------------------
    character(25) :: volumefile
    ! -------------------------------------------------------------------------------------------

    ! Determine volume name
    ! ====================
    volumefile='./check_volume_bl'//trim(numchar(iblc_pp))//'.bin'

    ! Prepare number and name of data
    ! ===============================
    ndata=volume_check_fstt%nvar
    allocate(varlist(ndata),dataname(ndata))
    do m=1,ndata
       dataname(m)=trim(volume_check_fstt%var(m))
    enddo

    ! volume indices
    ! ==============
    nx1l= volume_check_fstt%tectype%nx1
    nx2l= volume_check_fstt%tectype%nx2
    ny1l= volume_check_fstt%tectype%ny1
    ny2l= volume_check_fstt%tectype%ny2
    nz1l= volume_check_fstt%tectype%nz1
    nz2l= volume_check_fstt%tectype%nz2

    ! allocate dummy variable for datas
    ! ---------------------------------
    allocate(dummy(nx1l:nx2l,ny1l:ny2l,nz1l:nz2l,ndata))

    ! Fill data in temporary file dummy
    ! =================================
    ! List of possible variables: prs,uu,vv,ww,rho,Tmp,div,Mach,Gamma,Frhov,Grhow,udf

    ! Initialize counter for user-defined (udf) variables
    cpt=1

    do m=1,ndata
       select case (trim(volume_check_fstt%var(m)))
       case('check') ! extr_density
          do k=nz1l,nz2l
             do j=ny1l,ny2l
                do i=nx1l,nx2l
                   dummy(i,j,k,m) = extr_density(i,j,k)
                   ! dummy(i,j,k,m) = extr_density(i,j,k)
                enddo
             enddo
          enddo
       case default
          call mpistop('variable '//trim(volume_check_fstt%var(m)) &
               //' not listed in write_volume_check [mod_pp_fstt_write.f90]',0)
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
    if (iproc.eq.0) print *,'Writing volume... '//volumefile
    filetype=2 ! [0:full; 1:only grid; 2:only solution]
    call write_tec(volumefile,filetype,tstar,volume_check_fstt%tectype)

    ! Free temporary memory
    ! =====================
    deallocate(varlist,dataname,dummy)

  end subroutine write_volume_check


  !==============================================================================================
  subroutine write_grid_pp_fstt
  !==============================================================================================
    !> author: AB
    !> date: May 2022
    !> Write grid volume in tecplot format
  !==============================================================================================
    use mod_utils
    implicit none
    ! -------------------------------------------------------------------------------------------
    integer :: i,j,k,m,filetype
    real(wp) :: soltime
    logical :: iexist
    character(len=30) :: gridfile
    ! -------------------------------------------------------------------------------------------

    ! =========================================
    ! Write GRID tecplot format
    ! =========================================

    ! Name of the grid file
    ! =====================
    ! can be defined as optional input of the subroutine but never used at the present time
    gridfile='./grid_bl'//trim(numchar(iblc_pp))//'.plt'

    ! Prepare number and name of data: x,y,z coordinates (Att. CARTESIAN only)
    ! ===============================
    ndata=3
    allocate(dataname(ndata),varlist(ndata))
    dataname(1:ndata) = (/'X','Y','Z'/)

    ! Write volume grid
    ! =================
    ! check if grid file already exists
    ! ---------------------------------
    if (iproc==0) then
       inquire(file=gridfile,exist=iexist)
       if (iexist) then
          call system('rm '//trim(gridfile))
       endif
    endif

    ! if not: write GRID
    ! ------------------
    ! allocate and initialize dummy variable for datas
    allocate(dummy(nx,ny,nz,ndata))
    dummy = 0.0_wp

    ! fill data x,y,z in dummy array
    if (is_curv) then
       do i=1,nx
          do j=1,ny
             dummy(i,j,:,1)=xc(i,j)
             dummy(i,j,:,2)=yc(i,j)
          enddo
       enddo
    else
       do i=1,nx
          dummy(i,:,:,1)=x(i)
       enddo
       do j=1,ny
          dummy(:,j,:,2)=y(j)
       enddo
    endif
    do k=1,nz
       dummy(:,:,k,3)=z(k)
    enddo

    ! fill varlist (pointer on datas)
    do m=1,ndata
       varlist(m)%data=>dummy(:,:,:,m)
       varlist(m)%name= dataname(m)
    enddo

    ! filetype [0:full; 1:only grid; 2:only solution]
    filetype=1
    ! non-dimensional time
    soltime=0.0_wp

    ! write grid block in separate file
    if (iproc.eq.0) write(*,*) 'Writing gridfile ~>',trim(gridfile)
    call write_tec(gridfile,filetype,soltime,field(nob(iproc)))

    ! free temporary data
    deallocate(dataname,dummy,varlist)

  end subroutine write_grid_pp_fstt


  !==============================================================================================
  subroutine write_snapshot_pp_fstt(isn,split_num)
  !==============================================================================================
    !> author: AB
    !> date: May 2022
    !> Write volume # isn
  !==============================================================================================
    use mod_utils
    use mod_pp_vorticity
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Input: index of volume to be written
    integer, intent(in) :: isn,split_num
    ! -------------------------------------------------------------------------------------------
    ! Local variables
    integer :: i,j,k,m,filetype,ndata,cpt
    integer :: nx1l,nx2l,ny1l,ny2l,nz1l,nz2l
    ! -------------------------------------------------------------------------------------------
    character(len=30) :: split_name='./volume_001_sol_bl1_partx.bin'
    character(:), allocatable :: snapshotfile
    ! -------------------------------------------------------------------------------------------

    ! Determine snapshot name
    ! =======================
    if ((split_num.eq.0).or.(snap_pp_fstt(isn)%type.ne.3)) then
       if (snap_pp_fstt(isn)%stamp) then
          ! name of the snapshot with timestamp
          ! -----------------------------------
          cpt = ivol_2read
          if (.not.snap_pp_fstt(isn)%tectype%is_IOtec_write) then
             snapshotfile=trim(snap_pp_fstt(isn)%type_name)//'_'//snap_pp_fstt(isn)%tectype%zonename(6:8)//'_bl'//trim(numchar(iblc_pp))//'_pp_fstt_'//trim(numchar(cpt))//'.bin'
          else
             snapshotfile=trim(snap_pp_fstt(isn)%type_name)//'_'//snap_pp_fstt(isn)%tectype%zonename(6:8)//'_bl'//trim(numchar(iblc_pp))//'_pp_fstt_'//trim(numchar(cpt))//'.plt'
          endif
       else
          ! name of the snapshot without timestamp
          ! --------------------------------------
          if (.not.snap_pp_fstt(isn)%tectype%is_IOtec_write) then
             snapshotfile=trim(snap_pp_fstt(isn)%type_name)//'_'//snap_pp_fstt(isn)%tectype%zonename(6:8)//'_bl'//trim(numchar(iblc_pp))//'_pp_fstt.bin'
          else
             snapshotfile=trim(snap_pp_fstt(isn)%type_name)//'_'//snap_pp_fstt(isn)%tectype%zonename(6:8)//'_bl'//trim(numchar(iblc_pp))//'_pp_fstt.plt'
          endif
       endif
    else
       snapshotfile = trim(split_name)
       snapshotfile(26:26) = trim(numchar(split_num))
    endif

    ! Prepare number and name of data
    ! ===============================
    ndata=snap_pp_fstt(isn)%nvar
    allocate(varlist(ndata),dataname(ndata))
    do m=1,ndata
       dataname(m)=trim(snap_pp_fstt(isn)%var(m))
    enddo

    ! volume indices
    ! ==============
    nx1l= snap_pp_fstt(isn)%tectype%nx1
    nx2l= snap_pp_fstt(isn)%tectype%nx2
    ny1l= snap_pp_fstt(isn)%tectype%ny1
    ny2l= snap_pp_fstt(isn)%tectype%ny2
    nz1l= snap_pp_fstt(isn)%tectype%nz1
    nz2l= snap_pp_fstt(isn)%tectype%nz2

    ! allocate dummy variable for datas
    ! ---------------------------------
    allocate(dummy(nx1l:nx2l,ny1l:ny2l,nz1l:nz2l,ndata))

    ! Fill data in temporary file dummy
    ! =================================
    ! List of possible variables: uu,l2

    ! Initialize counter for user-defined (udf) variables
    cpt=1

    do m=1,ndata
       select case (trim(snap_pp_fstt(isn)%var(m)))
       case( 'uu') ! streamwise velocity
          do k=nz1l,nz2l
             do j=ny1l,ny2l
                do i=nx1l,nx2l
                   dummy(i,j,k,m)= uu(i,j,k)
                enddo
             enddo
          enddo
       case( 'vv') ! crossflow velocity
          do k=nz1l,nz2l
             do j=ny1l,ny2l
                do i=nx1l,nx2l
                   dummy(i,j,k,m)= vv(i,j,k)
                enddo
             enddo
          enddo
       case( 'ww') ! spanwise velocity
          do k=nz1l,nz2l
             do j=ny1l,ny2l
                do i=nx1l,nx2l
                   dummy(i,j,k,m)= ww(i,j,k)
                enddo
             enddo
          enddo
       case( 'up') ! streamwise velocity fluctuations
          do k=nz1l,nz2l
             do j=ny1l,ny2l
                do i=nx1l,nx2l
                   dummy(i,j,k,m)= uu(i,j,k) - stats_proc(i,j,2)
                enddo
             enddo
          enddo
       case( 'uut') ! tangential velocity fluctuations
          do k=nz1l,nz2l
             do j=ny1l,ny2l
                do i=nx1l,nx2l
                   ! dummy(i,j,k,m)= (uu(i,j,k)-stats_proc(i,j,2))*nyn_jmin(i) - (vv(i,j,k)-stats_proc(i,j,3))*nxn_jmin(i)
                   dummy(i,j,k,m) = uut(i,j,k)
                enddo
             enddo
          enddo
       case( 'uun') ! normal velocity fluctuations
          do k=nz1l,nz2l
             do j=ny1l,ny2l
                do i=nx1l,nx2l
                   ! dummy(i,j,k,m)= (uu(i,j,k)-stats_proc(i,j,2))*nxn_jmin(i) + (vv(i,j,k)-stats_proc(i,j,3))*nyn_jmin(i)
                   dummy(i,j,k,m) = uun(i,j,k)
                enddo
             enddo
          enddo
       case( 'vort') ! vorticity
          allocate(vort(nx,ny,nz))
          call compute_vorticity3d
          do k=nz1l,nz2l
             do j=ny1l,ny2l
                do i=nx1l,nx2l
                   dummy(i,j,k,m)= vort(i,j,k)
                enddo
             enddo
          enddo
          deallocate(vort)
       case( 'QQ') ! Q vorticity criterion
          allocate(QQ(nx,ny,nz))
          call compute_Q_criterion
          do k=nz1l,nz2l
             do j=ny1l,ny2l
                do i=nx1l,nx2l
                   dummy(i,j,k,m)= QQ(i,j,k)
                enddo
             enddo
          enddo
          deallocate(QQ)
       case( 'l2') ! lambda2 vorticity criterion
          allocate(lambda2(nx,ny,nz))
          call compute_l2_criterion
          do k=nz1l,nz2l
             do j=ny1l,ny2l
                do i=nx1l,nx2l
                   dummy(i,j,k,m)= lambda2(i,j,k)
                enddo
             enddo
          enddo
          deallocate(lambda2)
       case('ltbl') ! binary discrimination of laminar-turbulent region
          do k=nz1l,nz2l
             do j=ny1l,ny2l
                do i=nx1l,nx2l
                   dummy(i,j,k,m)=ltbl(i,j,k)
                enddo
             enddo
          enddo
       case('extr') ! extremum density
          do k=nz1l,nz2l
             do j=ny1l,ny2l
                do i=nx1l,nx2l
                   dummy(i,j,k,m)=extr_density2(i,j,k)
                enddo
             enddo
          enddo
       case default
          call mpistop('variable '//trim(snap_pp_fstt(isn)%var(m)) &
               //' not listed in write_snapshot_pp_fstt [mod_pp_fstt_write.f90]',0)
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
    call write_tec(snapshotfile,filetype,tstar,snap_pp_fstt(isn)%tectype)

    ! Free temporary memory
    ! =====================
    deallocate(varlist,dataname,dummy)

  end subroutine write_snapshot_pp_fstt

  !==============================================================================================
  subroutine check_snapshot_pp_fstt(isn)
  !==============================================================================================
    !> author: AB
    !> date: June 2022
    !> Check file existence snapshot_pp_fstt # isn
  !==============================================================================================
    use mod_utils
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Input: index of snapshot to be written
    integer, intent(in) :: isn
    ! -------------------------------------------------------------------------------------------
    integer :: iexist
    ! ! -------------------------------------------------------------------------------------------
    character(:), allocatable :: snapshotfile
    ! -------------------------------------------------------------------------------------------

    if (snap_pp_fstt(isn)%stamp) return

    ! Determine snapshot name
    ! =======================
    if (snap_pp_fstt(isn)%type.ne.3) then
       ! name of the snapshot without timestamp
       ! --------------------------------------
       if (.not.snap_pp_fstt(isn)%tectype%is_IOtec_write) then
          snapshotfile=trim(snap_pp_fstt(isn)%type_name)//'_'//snap_pp_fstt(isn)%tectype%zonename(6:8)//'_bl'//trim(numchar(iblc_pp))//'_pp_fstt.bin'
       else
          snapshotfile=trim(snap_pp_fstt(isn)%type_name)//'_'//snap_pp_fstt(isn)%tectype%zonename(6:8)//'_bl'//trim(numchar(iblc_pp))//'_pp_fstt.plt'
       endif
    endif

    ! Verification of file existence
    if (iproc==0) then
       inquire( file=trim(snapshotfile), exist=iexist )
       if (iexist) then
          call system('rm '//trim(snapshotfile))
       endif
    endif

  end subroutine check_snapshot_pp_fstt

end module mod_pp_fstt_write
