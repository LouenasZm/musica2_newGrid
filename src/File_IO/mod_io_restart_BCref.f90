!==============================================================================================
module mod_io_restart_BCref
!==============================================================================================
  !> author: AB
  !> date: June 2023
  !> Module to read/write restart_BCref.bin
  !> Adapted from mod_io_restartTD
!==============================================================================================
  use mod_mpi_part
  use mod_flow     ! <- for grid & solution
  use warnstop
  implicit none
  ! -------------------------------------------------------------------------------------------
  integer, parameter :: nvar=10
  integer(kind=MPI_OFFSET_KIND), dimension(3,2) :: offset_bcRef
  integer, dimension(3,2) :: type_mat_bc_nogh,type_mat_bc_view
  ! -------------------------------------------------------------------------------------------
  
contains

  !==============================================================================================
  subroutine init_io_restart_BCref
  !==============================================================================================
    !> author: AB
    !> date: June 2023
    !> Initialize the sub-communicator an I/O-type for restart_BCref outputs
    !> Adapted from init_io_restartTD
  !==============================================================================================
    use mod_block
    use mod_constant ! <- for is_init_2D_3D indicator
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Local variables
    integer :: nb,nbl,m1,m2,nz_,ngz_
    integer :: nbytes_dbl
    integer :: bsize_i,bsize_j,bsize_k
    integer :: offset
    integer, dimension(:), allocatable :: offset_bl
    integer, dimension(3) :: shape_nogh,start_nogh
    integer, dimension(3) :: shape_glob,shape_local,start_local
    ! -------------------------------------------------------------------------------------------

    ! if ((is_init_2D3D).and.(.not.is_BCref_init)) then
    !    nz_=1
    !    ngz_=1
    ! else
       nz_=nz
       ngz_=ngz
    ! endif

    ! Define offsets
    ! ==============

    ! Number of current block
    ! -----------------------
    nbl=nob(iproc)

    ! Buffer size for faces
    ! ---------------------
    
    ! Get size of the type MPI_DOUBLE_PRECISION
    call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,nbytes_dbl,info)
    ! i-direction (imin,imax)
    bsize_i=ngy*ngz_*nvar*nbytes_dbl
    ! j-direction (jmin,jmax)
    bsize_j=ngx*ngz_*nvar*nbytes_dbl
    ! k-direction (kmin,kmax)
    bsize_k=ngx*ngy*nvar*nbytes_dbl

    ! First pass to define global offset per block
    ! --------------------------------------------
    offset=0
    ! add all active face for current block
    if (BC_face(1,1)%is_mean_ref) offset=offset+bsize_i
    if (BC_face(1,2)%is_mean_ref) offset=offset+bsize_i
    if (BC_face(2,1)%is_mean_ref) offset=offset+bsize_j
    if (BC_face(2,2)%is_mean_ref) offset=offset+bsize_j
    if (BC_face(3,1)%is_mean_ref) offset=offset+bsize_k
    if (BC_face(3,2)%is_mean_ref) offset=offset+bsize_k
    ! share value
    allocate(offset_bl(0:nproc-1))
    offset_bl=0
    call MPI_ALLGATHER(offset,1,MPI_INTEGER,offset_bl,1,MPI_INTEGER,COMM_global,info)

    ! Define global offset for previous blocks
    ! ----------------------------------------
    ! sum for blocks with lower rank
    offset=0
    do nb=1,nbl-1
       offset=offset+offset_bl(iproc_leader(nb))
    enddo
    
    ! Define offsets for each faces
    ! -----------------------------
    offset_bcRef = 0_MPI_OFFSET_KIND+offset

    offset=0
    if (BC_face(1,1)%is_mean_ref) offset=bsize_i
    offset_bcRef(1,2) = offset_bcRef(1,2)+offset
    if (BC_face(1,2)%is_mean_ref) offset=offset+bsize_i
    offset_bcRef(2,1) = offset_bcRef(2,1)+offset
    if (BC_face(2,1)%is_mean_ref) offset=offset+bsize_j
    offset_bcRef(2,2) = offset_bcRef(2,2)+offset
    if (BC_face(2,2)%is_mean_ref) offset=offset+bsize_j
    offset_bcRef(3,1) = offset_bcRef(3,1)+offset
    if (BC_face(3,1)%is_mean_ref) offset=offset+bsize_k
    offset_bcRef(3,2) = offset_bcRef(3,2)+offset
 
    ! Creation of derived types for MPI-IO
    ! ====================================

    m1=1 ! faces imin and imax
    do m2=1,2

       ! Creation of the derived type type_local that defines the array without ghost cells
       ! ----------------------------------------------------------------------------------

       ! Shape of the array without ghost cells
       shape_nogh= (/ 1, ny, nz_ /)

       ! Starting coordinates for whole array (i.e. array=subarray)
       start_nogh= (/ 0, 0, 0 /)

       ! Creation of derived type type_mat_bc_nogh
       ! (This allows to remove ghost-cells from initial array)
       call MPI_TYPE_CREATE_SUBARRAY(3,shape_nogh,shape_nogh,start_nogh,&
            MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,type_mat_bc_nogh(m1,m2),info)
       ! Commit of type_mat_bc_nogh
       call MPI_TYPE_COMMIT(type_mat_bc_nogh(m1,m2),info)

       ! Creation of type type_mat_nogh_view to set the view on the file
       ! ---------------------------------------------------------------

       ! Shape of the global array
       shape_glob= (/ 1, ngy, ngz_ /)

       ! Shape of the local array
       shape_local= (/ 1, ny, nz_ /)

       ! Starting coordinates for local array
       start_local= (/ 0, coord(2)*ny, coord(3)*nz_ /)
       if (nx ==1) start_local(1) = 0
       if (ny ==1) start_local(2) = 0
       if (nz_==1) start_local(3) = 0

       ! Creation of derived type type_mat_bc_view
       call MPI_TYPE_CREATE_SUBARRAY(3,shape_glob,shape_local,start_local,&
            MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,type_mat_bc_view(m1,m2),info)
       ! Commit of type_mat_bc_view
       call MPI_TYPE_COMMIT(type_mat_bc_view(m1,m2),info)

    enddo

    m1=2 ! faces jmin and jmax
    do m2=1,2

       ! Creation of the derived type type_local that defines the array without ghost cells
       ! ----------------------------------------------------------------------------------

       ! Shape of the array without ghost cells
       shape_nogh= (/ nx, 1, nz_ /)

       ! Starting coordinates for whole array (i.e. array=subarray)
       start_nogh= (/ 0, 0, 0 /)

       ! Creation of derived type type_mat_bc_nogh
       ! (This allows to remove ghost-cells from initial array)
       call MPI_TYPE_CREATE_SUBARRAY(3,shape_nogh,shape_nogh,start_nogh,&
            MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,type_mat_bc_nogh(m1,m2),info)
       ! Commit of type_mat_bc_nogh
       call MPI_TYPE_COMMIT(type_mat_bc_nogh(m1,m2),info)

       ! Creation of type type_mat_nogh_view to set the view on the file
       ! ---------------------------------------------------------------

       ! Shape of the global array
       shape_glob= (/ ngx, 1, ngz_ /)

       ! Shape of the local array
       shape_local= (/ nx, 1, nz_ /)

       ! Starting coordinates for local array
       start_local= (/ coord(1)*nx, 0, coord(3)*nz_ /)
       if (nx ==1) start_local(1) = 0
       if (ny ==1) start_local(2) = 0
       if (nz_==1) start_local(3) = 0

       ! Creation of derived type type_mat_bc_view
       call MPI_TYPE_CREATE_SUBARRAY(3,shape_glob,shape_local,start_local,&
            MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,type_mat_bc_view(m1,m2),info)
       ! Commit of type_mat_bc_view
       call MPI_TYPE_COMMIT(type_mat_bc_view(m1,m2),info)

    enddo

   m1=3 ! faces kmin and kmax
   do m2=1,2
       ! Creation of the derived type type_local that defines the array without ghost cells
       ! ----------------------------------------------------------------------------------

       ! Shape of the array without ghost cells
       shape_nogh= (/ nx, ny, 1 /)

       ! Starting coordinates for whole array (i.e. array=subarray)
       start_nogh= (/ 0, 0, 0 /)

       ! Creation of derived type type_mat_bc_nogh
       ! (This allows to remove ghost-cells from initial array)
       call MPI_TYPE_CREATE_SUBARRAY(3,shape_nogh,shape_nogh,start_nogh,&
           MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,type_mat_bc_nogh(m1,m2),info)
       ! Commit of type_mat_bc_nogh
       call MPI_TYPE_COMMIT(type_mat_bc_nogh(m1,m2),info)

       ! Creation of type type_mat_nogh_view to set the view on the file
       ! ---------------------------------------------------------------

       ! Shape of the global array
       shape_glob= (/ ngx, ngy, 1 /)

       ! Shape of the local array
       shape_local= (/ nx, ny, 1 /)

       ! Starting coordinates for local array
       start_local= (/ coord(1)*nx, coord(2)*ny, 0 /)
       if (nx ==1) start_local(1) = 0
       if (ny ==1) start_local(2) = 0
       if (nz_==1) start_local(3) = 0

       ! Creation of derived type type_mat_bc_view
       call MPI_TYPE_CREATE_SUBARRAY(3,shape_glob,shape_local,start_local,&
           MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,type_mat_bc_view(m1,m2),info)
       ! Commit of type_mat_bc_view
       call MPI_TYPE_COMMIT(type_mat_bc_view(m1,m2),info)

    enddo

  end subroutine init_io_restart_BCref

  !==============================================================================================
  subroutine write_restart_BCref(filename)
  !==============================================================================================
    !> author: AB
    !> date: June 2023
    !> Write variables for BCref
    !> Adapted from write_restartTD
  !==============================================================================================
    implicit none
    ! -------------------------------------------------------------------------------------------
    character(len=*), intent(in) :: filename
    ! -------------------------------------------------------------------------------------------
    ! Local variables
    integer :: n,m1,m2,fh,n1_,n2_
    integer, dimension(MPI_STATUS_SIZE) :: statut
    ! -------------------------------------------------------------------------------------------

    if (iproc.eq.0) print *,'Writing restart_BCref.bin ... '

    ! Initialization of IO
    ! ====================
    call init_io_restart_BCref

    ! Open file for writing
    ! =====================
    call MPI_FILE_OPEN(COMM_global,trim(filename), &
         MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,fh,info)

    ! Write data
    ! ==========
    do m1=1,3
       do m2=1,2
          ! Setting the view on the file
          call MPI_FILE_SET_VIEW(fh,offset_bcRef(m1,m2),MPI_DOUBLE_PRECISION, &
               type_mat_bc_view(m1,m2),'native',MPI_INFO_NULL,info)

          ! Write each variable if BC is_mean_ref
          if (BC_face(m1,m2)%is_mean_ref) then
             ! Determination of array size
             if (m1.eq.1) then
               n1_=ny; n2_=nz
             else if (m1.eq.2) then
               n1_=nx; n2_=nz
             else
               n1_=nx; n2_=ny
             endif

             ! MPI_FILE_WRITE
             do n=1,nvar
                call MPI_FILE_WRITE(fh,BC_face(m1,m2)%Uref(1:n1_,1:n2_,n),1,&
                     type_mat_bc_nogh(m1,m2),statut,info)
             enddo
          endif
       enddo
    enddo

    ! Close file
    ! ==========
    call MPI_FILE_CLOSE(fh,info)

    ! Free IO
    ! =======
    call free_restart_BCref

  end subroutine write_restart_BCref

  !==============================================================================================
  subroutine read_restart_BCref(filename)
  !==============================================================================================
    !> author: AB
    !> date: June 2023
    !> Read variables for BCref
    !> Adapted from read_restartTD
  !==============================================================================================
    implicit none
    ! -------------------------------------------------------------------------------------------
    character(len=*), intent(in) :: filename
    ! -------------------------------------------------------------------------------------------
    ! Local variables
    integer :: n,m1,m2,fh,n1_,n2_
    integer, dimension(MPI_STATUS_SIZE) :: statut
    ! -------------------------------------------------------------------------------------------

    if (iproc.eq.0) print *,'Reading restart_BCref.bin ... '

    ! Initialization of IO
    ! ====================
    call init_io_restart_BCref

    ! Open file for reading
    ! =====================
    call MPI_FILE_OPEN(COMM_global,trim(filename), &
         MPI_MODE_RDONLY,MPI_INFO_NULL,fh,info)

    if (info/=MPI_SUCCESS) then
       call mpistop('MPI_IO: Error in opening restart_BCref.bin',1)
    endif

    ! Read data
    ! =========
    do m1=1,3
       do m2=1,2
          ! Setting the view on the file
          call MPI_FILE_SET_VIEW(fh,offset_bcRef(m1,m2),MPI_DOUBLE_PRECISION, &
               type_mat_bc_view(m1,m2),'native',MPI_INFO_NULL,info)
          ! Read each variable if BC is T&D
          if (BC_face(m1,m2)%is_mean_ref) then
             ! Determination of array size
             if (m1.eq.1) then
               n1_=ny; n2_=nz
             else if (m1.eq.2) then
               n1_=nx; n2_=nz
             else
               n1_=nx; n2_=ny
             endif
             do n=1,nvar
                call MPI_FILE_READ(fh,BC_face(m1,m2)%Uref(1:n1_,1:n2_,n),1,&
                     type_mat_bc_nogh(m1,m2),statut,info)
             enddo
          endif
       enddo
    enddo

    ! Close file
    ! ==========
    call MPI_FILE_CLOSE(fh,info)

    ! Free IO
    ! =======
    call free_restart_BCref

  end subroutine read_restart_BCref

  !==============================================================================================
  subroutine free_restart_BCref
  !==============================================================================================
    !> author: XG
    !> date: May 2022
    !> Free MPI-IO types for BCref
    !> Adapted from free_restartTD
  !==============================================================================================
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Local variables
    integer :: m1,m2
    ! -------------------------------------------------------------------------------------------

    ! Free types for MPI-IO of restart_BCref
    ! --------------------------------------
    do m1=1,3
       do m2=1,2
          call MPI_TYPE_FREE(type_mat_bc_nogh(m1,m2),info)
          call MPI_TYPE_FREE(type_mat_bc_view(m1,m2),info)
       enddo
    enddo

  end subroutine free_restart_BCref

end module mod_io_restart_BCref
