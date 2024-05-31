!==============================================================================================
module mod_tecplot
!==============================================================================================
  !> author: Luca Sciacovelli
  !> date: April 2018
  !> Module to read and write tecplot binary files
!==============================================================================================
  use mpi
  use precision
  use warnstop
  implicit none
  ! -------------------------------------------------------------------------------------------
  type teciotype
     integer              :: ngx,ngy,ngz             ! Full domain size
     integer              :: nx,ny,nz                ! Local domain size
     integer              :: nx1,nx2,ny1,ny2,nz1,nz2 ! Local bounds
     integer              :: MPI_COMM                ! MPI communicator
     integer, allocatable :: coord(:)                ! coord
     ! Variables for MPI_IO
     integer :: type_mat,type_mat_nogh,type_mat_view
     character(len=8) :: zonename='field'
     integer :: strandID
     logical :: is_IOtec_read,is_IOtec_write
     ! added by XG for appending existing file
     logical :: is_app  ! append file
     logical :: restart ! restarting mode
     integer :: fh      ! file handle
     integer(8) :: disp    ! displacement for 1 record
     integer(kind=MPI_OFFSET_KIND) :: offset ! file offset
  end type teciotype
  ! -------------------------------------------------------------------------------------------
  type listtype
     real(wp), dimension(:,:,:), pointer :: data
     character(len=150) :: name
  end type listtype
  ! -------------------------------------------------------------------------------------------
  type(listtype), dimension(:), allocatable :: varlist
  ! -------------------------------------------------------------------------------------------

contains

  !==============================================================================================
  subroutine mod_io_init(ngx,ngy,ngz,nx,ny,nz,nxr,nyr,nzr,ngh,ndim,coord,is_IOtec_read,is_IOtec_write,field)
  !==============================================================================================
    !> author: Luca Sciacovelli
    !> date: April 2018 (modif XG 2020)
    !> Initializations of MPI Input/Outputs
    !> ~> definitions of derived types for MPI_IO
  !==============================================================================================
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Input/output variables
    integer, intent(in)            :: ngx,ngy,ngz,nx,ny,nz,nxr,nyr,nzr,ngh,ndim
    integer, intent(in)            :: coord(:)
    logical, intent(in)            :: is_IOtec_read,is_IOtec_write
    type(teciotype), intent(inout) :: field
    ! -------------------------------------------------------------------------------------------
    ! Local variables
    integer :: ierror,nbytes_dbl
    integer, dimension(3) :: shape_wgh,shape_nogh,start_wgh,start_nogh
    integer, dimension(3) :: shape_glob,shape_local,start_local
    ! -------------------------------------------------------------------------------------------

    ! Get size of the type MPI_DOUBLE_PRECISION
    call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,nbytes_dbl,ierror)

    ! Creation of the derived type type_local that defines the array without ghost cells
    ! ==================================================================================
    
    ! Shape of the array with ghost cells
    ! -----------------------------------
    shape_wgh= (/ nx+2*ngh, ny+2*ngh, nz+2*ngh /)
    if (nz.eq.1) shape_wgh= (/ nx+2*ngh, ny+2*ngh, 1 /)

    ! Shape of the array without ghost cells
    ! --------------------------------------
    shape_nogh= (/ nx, ny, nz /)

    ! Starting coordinates for array without ghost cells (w.r.t. global array)
    ! --------------------------------------------------
    start_wgh= (/ ngh, ngh, ngh /)
    if (nz.eq.1) start_wgh = (/ ngh, ngh, 0 /)

    ! Starting coordinates for whole array (i.e. array=subarray)
    ! ------------------------------------
    start_nogh= (/ 0, 0, 0 /)

    ! Creation of derived type type_mat_nogh
    ! (This allows to remove ghost-cells from initial array)
    ! ------------------------------------------------------
    call MPI_TYPE_CREATE_SUBARRAY(ndim,shape_nogh,shape_nogh,start_nogh,&
         MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,field%type_mat_nogh,ierror)
    ! Commit of type_mat_nogh
    call MPI_TYPE_COMMIT(field%type_mat_nogh,ierror)

    ! Creation of derived type type_mat
    ! ---------------------------------
    call MPI_TYPE_CREATE_SUBARRAY(ndim,shape_wgh,shape_nogh,start_wgh,&
         MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,field%type_mat,ierror)
    ! Commit of type_mat
    call MPI_TYPE_COMMIT(field%type_mat,ierror)

    ! Creation of type type_mat_nogh_view to set the view on the file
    ! ===============================================================

    ! Shape of the global array
    ! -------------------------
    shape_glob= (/ ngx, ngy, ngz /)

    ! Shape of the local array
    ! ------------------------
    shape_local= (/ nx, ny, nz /)

    ! Starting coordinates for local array
    ! ------------------------------------
    start_local= (/ coord(1)*nxr, coord(2)*nyr, coord(3)*nzr /)
    if (nx==1) start_local(1) = 0
    if (ny==1) start_local(2) = 0
    if (nz==1) start_local(3) = 0
    
    ! Creation of derived type type_mat_view
    ! --------------------------------------
    call MPI_TYPE_CREATE_SUBARRAY(ndim,shape_glob,shape_local,start_local,&
         MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,field%type_mat_view,ierror)
    ! Commit of type_mat_view
    call MPI_TYPE_COMMIT(field%type_mat_view,ierror)

    ! Initialize field structure
    ! ==========================
    
    ! Initialize dimensions
    field%ngx=ngx
    field%ngy=ngy
    field%ngz=ngz
    field%nx =nx
    field%ny =ny
    field%nz =nz

    ! Initialize coord
    if (.not.allocated(field%coord)) allocate(field%coord(ndim))
    field%coord=coord

    ! Initialize StrandID [-2 for Static StrandID in tecplot binary]
    field%strandID=-2

    ! Initialize IO type (bin or tec)
    field%is_IOtec_read =is_IOtec_read
    field%is_IOtec_write=is_IOtec_write

    ! Initialize attributes to append files
    field%is_app=.false.
    field%disp=0
    field%restart=.false.
    field%offset=0_MPI_OFFSET_KIND
    
  end subroutine mod_io_init

  !==============================================================================================
  subroutine write_tec(filename,filetype,soltime,field)
  !==============================================================================================
    !> author: Luca Sciacovelli
    !> date: April 2018 (modified XG to append existing file)
    !> Write a tecplot/fortran binary file
  !==============================================================================================
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Input/Output arguments
    character(len=*), intent(in) :: filename
    integer, intent(in) :: filetype
    real(wp) :: soltime
    type(teciotype) , intent(inout) :: field
    ! -------------------------------------------------------------------------------------------
    ! Local variables
    integer :: i,j,nvars,ierror,fh,myid_local
    real(kind=4) :: fzone,fhead
    character(len=8)  :: title
    character(len=50) :: bid
    real(wp), dimension(:), allocatable :: min_glob,max_glob
    integer, dimension(MPI_STATUS_SIZE) :: statut
    integer(kind=MPI_OFFSET_KIND) :: offset
    logical :: iexist
    ! -------------------------------------------------------------------------------------------

    ! Find rank in MPI communicator
    ! =============================
    if (field%MPI_COMM.eq.MPI_COMM_NULL) return

    call MPI_COMM_RANK(field%MPI_COMM,myid_local,ierror)

    ! Open file for writing and set offset to append in an existing file
    ! ==================================================================
    if (field%is_app) then  
       inquire(file=trim(filename),exist=iexist)
       if (.not.(iexist)) then
          ! Open file for writing
          call MPI_FILE_OPEN(field%MPI_COMM,trim(filename), &
               MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierror)
          ! Test if open fails
          if (ierror/=MPI_SUCCESS) then
             call mpistop('MPI_IO: Error in opening file '//trim(filename),1)
          endif
          field%fh=fh
          ! Initialize offset
          offset=0_MPI_OFFSET_KIND
          field%restart=.false.
       else
          if (field%restart) then
             ! Open file for writing
             call MPI_FILE_OPEN(field%MPI_COMM,trim(filename), &
                  MPI_MODE_WRONLY+MPI_MODE_APPEND,MPI_INFO_NULL,fh,ierror)
             ! Test if open fails
             if (ierror /= MPI_SUCCESS) then
                call mpistop('MPI_IO: Error in opening file '//trim(filename),1)
             endif
             field%fh=fh
             ! Returns the current position of the individual file pointer
             call MPI_FILE_GET_POSITION(fh,offset,ierror)
             field%restart=.false.
             field%offset=offset
          else
             fh=field%fh
             offset=field%offset
          endif
       endif
    else
       ! Open file for writing
       call MPI_FILE_OPEN(field%MPI_COMM,trim(filename), &
            MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierror)

       ! Test if open fails
       if (ierror/=MPI_SUCCESS) then
          call mpistop('MPI_IO: Error in opening file '// trim(filename), 1)
       endif
       ! Initialize offset
       offset=0_MPI_OFFSET_KIND
    endif
    
!!$    ! Open file for writing
!!$    ! =====================
!!$    call MPI_FILE_OPEN(field%MPI_COMM,trim(filename), &
!!$         MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierror)
!!$
!!$    ! Test if open fails
!!$    if (ierror/= MPI_SUCCESS) then
!!$       call mpistop('MPI_IO: Error in opening file '//trim(filename),1)
!!$    endif
!!$    
!!$    ! Initialize offset
!!$    ! =================
!!$    offset=0_MPI_OFFSET_KIND

    ! Initialize nvars
    ! ================
    nvars = size(varlist)

    ! Read TECPLOT header [only for .tec files]
    ! ===================
    if (field%is_IOtec_write) then

       fzone = 299.0
       fhead = 357.0

       title = 'mix'

       ! Tecplot needs to know minimum and maximum of each variable
       ! ----------------------------------------------------------
       allocate(min_glob(nvars),max_glob(nvars))
       min_glob=0.0_wp
       max_glob=0.0_wp

       do i=1,nvars
          call MPI_REDUCE(minval(varlist(i)%data),min_glob(i),1,&
               MPI_DOUBLE_PRECISION,MPI_MIN,0,field%MPI_COMM,ierror)

          call MPI_REDUCE( maxval(varlist(i)%data),max_glob(i),1,&
               MPI_DOUBLE_PRECISION,MPI_MAX,0,field%MPI_COMM,ierror)
       enddo

       ! Start reading Tecplot header
       ! ----------------------------
       if (myid_local.eq.0) then
          ! Magic Number
          call MPI_FILE_WRITE(fh,'#!TDV112',8,MPI_CHARACTER,statut,ierror)
          offset=offset+8
          ! Byte order of the reader
          call MPI_FILE_WRITE(fh,1,1,MPI_INTEGER,statut,ierror)
          offset=offset+4
          ! FileType (0:Full 1:Grid 2:Sol)
          call MPI_FILE_WRITE(fh,filetype,1,MPI_INTEGER,statut,ierror)
          offset=offset+4
          ! Title
          do i=1,len_trim(title)
             call MPI_FILE_WRITE(fh,ichar(title(i:i)),1,MPI_INTEGER,statut,ierror)
             offset=offset+4
          enddo
          ! Strings end with 0
          call MPI_FILE_WRITE(fh,0,1,MPI_INTEGER,statut,ierror)
          offset=offset+4
          ! Number of variables
          call MPI_FILE_WRITE(fh,nvars,1,MPI_INTEGER,statut,ierror)
          offset=offset+4
          ! Variable names
          do i=1,nvars
             bid=trim(varlist(i)%name)
             ! Name of variable
             do j=1,len_trim(bid)
                call MPI_FILE_WRITE(fh,ichar(bid(j:j)),1,MPI_INTEGER,statut,ierror)
                offset=offset+4
             enddo
             ! Strings end with 0
             call MPI_FILE_WRITE(fh,0,1,MPI_INTEGER,statut,ierror)
             offset=offset+4
          enddo
          ! Zone Marker. Value = 299.0
          call MPI_FILE_WRITE(fh,fzone,1,MPI_REAL4,statut,ierror)
          offset=offset+4
          ! Zone name
          do i=1,len_trim(field%zonename)
             call MPI_FILE_WRITE(fh,ichar(field%zonename(i:i)),1,MPI_INTEGER,statut,ierror)
             offset=offset+4
          enddo
          ! Strings end with 0
          call MPI_FILE_WRITE(fh, 0,1,MPI_INTEGER,statut,ierror)
          offset=offset+4
          ! No parent Zone
          call MPI_FILE_WRITE(fh,-1,1,MPI_INTEGER,statut,ierror)
          offset=offset+4
          ! StrandID
          call MPI_FILE_WRITE(fh,field%strandID,1,MPI_INTEGER,statut,ierror)
          offset=offset+4
          ! Solution time
          call MPI_FILE_WRITE(fh,soltime,1,MPI_DOUBLE_PRECISION,statut,ierror)
          offset=offset+8
          ! Not used
          call MPI_FILE_WRITE(fh,-1,1,MPI_INTEGER,statut,ierror)
          offset=offset+4
          ! DataPacking (0:Block)
          call MPI_FILE_WRITE(fh, 0,1,MPI_INTEGER,statut,ierror)
          offset=offset+4
          ! Var location (0:Nodes)
          call MPI_FILE_WRITE(fh, 0,1,MPI_INTEGER,statut,ierror)
          offset=offset+4
          ! Are raw local 1-to-1 face neighbours supplied? (0: False)
          call MPI_FILE_WRITE(fh, 0,1,MPI_INTEGER,statut,ierror)
          offset=offset+4
          ! Number of miscellaneous user-defined face neighbor connections
          call MPI_FILE_WRITE(fh, 0,1,MPI_INTEGER,statut,ierror)
          offset=offset+4
          ! Imax
          call MPI_FILE_WRITE(fh,field%ngx,1,MPI_INTEGER,statut,ierror)
          offset=offset+4
          ! Jmax
          call MPI_FILE_WRITE(fh,field%ngy,1,MPI_INTEGER,statut,ierror)
          offset=offset+4
          ! Kmax
          call MPI_FILE_WRITE(fh,field%ngz,1,MPI_INTEGER,statut,ierror)
          offset=offset+4
          ! No more Auxiliary name/value pairs
          call MPI_FILE_WRITE(fh, 0,1,MPI_INTEGER,statut,ierror)
          offset=offset+4
          ! End of Header
          call MPI_FILE_WRITE(fh,fhead,1,MPI_REAL4,statut,ierror)
          offset=offset+4
          ! Data section
          call MPI_FILE_WRITE(fh,fzone,1,MPI_REAL4,statut,ierror)
          offset=offset+4
          ! Variable data format (2: double)
          do i=1,nvars
             call MPI_FILE_WRITE(fh,2,1,MPI_INTEGER,statut,ierror)
             offset=offset+4
          enddo
          ! Has passive variables (0: no)
          call MPI_FILE_WRITE(fh, 0,1,MPI_INTEGER,statut,ierror)
          offset=offset+4
          ! Has variable sharing (0: no)
          call MPI_FILE_WRITE(fh, 0,1,MPI_INTEGER,statut,ierror)
          offset=offset+4
          ! Share connectivity list with (-1: no)
          call MPI_FILE_WRITE(fh,-1,1,MPI_INTEGER,statut,ierror)
          offset=offset+4
          ! Min/Max of each variable
          do i=1,nvars
             call MPI_FILE_WRITE(fh,min_glob(i),1,MPI_DOUBLE_PRECISION,statut,ierror)
             call MPI_FILE_WRITE(fh,max_glob(i),1,MPI_DOUBLE_PRECISION,statut,ierror)
             offset=offset+2*8
          enddo
       endif
       deallocate(min_glob,max_glob)
       ! share offset value on all procs of the communicator
       ! ---------------------------------------------------
       call MPI_BCAST(offset,1,MPI_INTEGER,0,field%MPI_COMM,ierror)
    endif

    ! Setting the view on the file
    ! ============================
    call MPI_FILE_SET_VIEW(fh,offset,MPI_DOUBLE_PRECISION,field%type_mat_view, &
                           'native',MPI_INFO_NULL,ierror)

    ! Write data for all vars in list
    ! ===============================
    do i=1,nvars
       call MPI_FILE_WRITE_ALL(fh,varlist(i)%data,1,field%type_mat_nogh,statut,ierror)
       if (ierror /= MPI_SUCCESS) then
          call mpistop('MPI_IO: Error in writing file '//trim(filename),1)
       endif
    enddo

    ! Set offset for new writing or Close file
    ! ========================================
    if (field%is_app) then
       field%offset=field%offset+field%disp
       !!print *,'offset',field%offset
    else
       call MPI_FILE_CLOSE(fh,ierror) ! Close file
    endif

  end subroutine write_tec

  !==============================================================================================
  subroutine read_tec(filename,field)
  !==============================================================================================
    !> author: Luca Sciacovelli
    !> date: April 2018
    !> Read a tecplot/fortran binary file
  !==============================================================================================
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Input/Output arguments
    character(len=*), intent(in) :: filename
    type(teciotype) , intent(inout) :: field
     ! -------------------------------------------------------------------------------------------
    ! Local variables
    integer :: i,j,nvars,ierror,fh,dummy,myid_local
    real(kind=4) :: fzone,fhead
    integer, dimension(MPI_STATUS_SIZE) :: statut
    integer(kind=MPI_OFFSET_KIND) :: offset
    ! -------------------------------------------------------------------------------------------

    ! Find rank in MPI communicator
    ! =============================
    if (field%MPI_COMM.eq.MPI_COMM_NULL) return

    call MPI_COMM_RANK(field%MPI_COMM,myid_local,ierror)

    ! Initialize nvars
    ! ================
    nvars=size(varlist)

    ! Open file for reading and set offset to append in an existing file
    ! ==================================================================
    if (field%is_app) then  
       if (field%restart) then
          ! Open file for reading
          call MPI_FILE_OPEN(field%MPI_COMM,trim(filename), &
               MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierror)
          ! Test if open fails
          if (ierror /= MPI_SUCCESS) then
             call mpistop('MPI_IO: Error in opening file '//trim(filename),1)
          endif
          field%fh=fh
          ! Returns the current position of the individual file pointer
          ! field%offset equal 0 by default, =/= if user defined
          !call MPI_FILE_GET_POSITION(fh,offset,ierror)
          ! offset=0_MPI_OFFSET_KIND
          ! field%offset=offset
          offset=field%offset

          field%restart=.false.
       else
          fh=field%fh
          offset=field%offset
       endif
    else
       ! Open file for reading
       call MPI_FILE_OPEN(field%MPI_COMM,trim(filename),MPI_MODE_RDONLY, &
            MPI_INFO_NULL,fh,ierror)
       
       ! Test if open fails
       if (ierror /= MPI_SUCCESS) then
          call mpistop('MPI_IO: Error in opening file '// trim(filename),1)
       endif

       ! Initialize offset
       offset=0_MPI_OFFSET_KIND
    endif

!!$    ! Open file for reading
!!$    ! =====================
!!$    call MPI_FILE_OPEN(field%MPI_COMM,trim(filename),MPI_MODE_RDONLY, &
!!$                       MPI_INFO_NULL,fh,ierror)
!!$    if (ierror /= MPI_SUCCESS) then
!!$       call mpistop('MPI_IO: Error in opening file '// trim(filename),1)
!!$    endif
!!$    
!!$    ! Initialize offset
!!$    ! =================
!!$    offset=0_MPI_OFFSET_KIND

    ! Write TECPLOT header [only for .tec files]
    ! ====================
    if (field%is_IOtec_read) then

       if (myid_local.eq.0) then
          fzone=299.0
          fhead=357.0
          ! We have to read part of the header because we don't know the length of
          ! title,variable and zone names,and we need them to compute the offset
          dummy=999
          do while (dummy.ne.0)
             ! read useless stuff,up to number of variables
             call MPI_FILE_READ(fh,dummy,1,MPI_INTEGER,statut,ierror)
             offset=offset+4
          enddo
          ! Number of variables
          ! -------------------
          call MPI_FILE_READ(fh,nvars,1,MPI_INTEGER,statut,ierror)
          offset=offset+4
          ! Name of variables
          ! -----------------
          do i = 1,nvars
             varlist(i)%name = ''
             dummy = 999
             j=0
             do
                j=j+1
                call MPI_FILE_READ(fh,dummy,1,MPI_INTEGER,statut,ierror)
                offset=offset+4
                if (dummy.eq.0) exit
                varlist(i)%name(j:j) = char( dummy )
             enddo
          enddo
          ! Zone Marker. Value = 299.0
          ! -----------
          call MPI_FILE_READ(fh,fzone,1,MPI_REAL4,statut,ierror)
          offset=offset + 4
          ! Zone name
          ! ---------
          dummy=999
          do while (dummy.ne.0)
             call MPI_FILE_READ( fh,dummy,1,MPI_INTEGER,statut,ierror)
             offset=offset+4
          enddo
          ! The rest of the header is fixed [see write_tec if you need]
          ! -------------------------------
          offset=offset+(1+2*nvars)*8+(16+nvars)*4
       endif
       ! share offset value on all procs of the communicator
       ! ---------------------------------------------------
       call MPI_BCAST(offset,1,MPI_INTEGER,0,field%MPI_COMM,ierror)
    endif

    ! Setting the view on the file
    ! ============================
    call MPI_FILE_SET_VIEW(fh,offset,MPI_DOUBLE_PRECISION,field%type_mat_view,&
                           'native',MPI_INFO_NULL,ierror)

    ! Read data for all vars in list
    ! ==============================
    do i=1,nvars
       !if (myid_local.eq.0) write(*,*) 'Reading Variable ~> ',i!,varlist(i)%name
       call MPI_FILE_READ_ALL(fh,varlist(i)%data,1,field%type_mat_nogh,statut,ierror)
    enddo

    ! Set offset for new writing or Close file
    ! ========================================
    if (field%is_app) then
       field%offset=field%offset+field%disp
    else
       call MPI_FILE_CLOSE(fh,ierror) ! Close file
    endif

  end subroutine read_tec

end module mod_tecplot
