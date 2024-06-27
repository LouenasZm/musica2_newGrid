!==============================================================================================
module mod_io_restartR
!==============================================================================================
  !> author: CM
  !> date: January 2024
  !> Module to read/write restartR.bin (restart file for perturbative Riemann BCs)
  !> Adapted from mod_io_restartTD_1pt
!==============================================================================================
  use mod_mpi_part
  use mod_flow             ! <- for grid & solution
  use warnstop

  implicit none
  ! -------------------------------------------------------------------------------------------
  ! restartR TYPE
  ! -------------------------------------------------------------------------------------------
  type restartR_type
     ! local indices
     integer :: i1,i2,j1,j2,k1,k2
     ! variables for MPI_IO
     integer :: type_mat_nogh,type_mat_view
     integer(kind=MPI_OFFSET_KIND) :: offset
  end type restartR_type
  ! -------------------------------------------------------------------------------------------
  integer, parameter :: nvarR=3
  ! !> restartR is defined as a restartTD_type_1pt type for simplicity, eventhough
  !    no ghost cells are used so the subarray created IS the array we want to send
  type(restartR_type), dimension(3,2) :: restartR
  ! -------------------------------------------------------------------------------------------
  
contains

  !==============================================================================================
  subroutine init_io_restartR
  !==============================================================================================
    !> author: CM
    !> date: January 2024
    !> Initialize the sub-communicator an I/O-type for restartR outputs
  !==============================================================================================
    use mod_block
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Local variables
    integer :: nb,nbl,m1,m2
    integer :: nbytes_dbl
    integer :: bsize_i,bsize_j,bsize_k
    integer :: offset
    integer, dimension(:), allocatable :: offset_bl
    integer, dimension(3) :: shape_wgh,shape_nogh,start_wgh,start_nogh
    integer, dimension(3) :: shape_glob,shape_local,start_local
    ! -------------------------------------------------------------------------------------------
    ! Initialize local indices
    ! ========================
    do m1=1,3
       do m2=1,2
          restartR(m1,m2)%i1=1
          restartR(m1,m2)%i2=nx
          restartR(m1,m2)%j1=1
          restartR(m1,m2)%j2=ny
          restartR(m1,m2)%k1=1
          restartR(m1,m2)%k2=nz
       enddo
    enddo

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
    bsize_i=ngy*ngz*nvarR*nbytes_dbl
    ! j-direction (jmin,jmax)
    bsize_j=ngx*ngz*nvarR*nbytes_dbl
    ! k-direction (kmin,kmax)
    bsize_k=ngx*ngy*nvarR*nbytes_dbl
    ! First pass to define global offset per block
    ! --------------------------------------------
    offset=0
    ! add all active face for current block
    if (bl(nbl)%BC(1).eq.-41) offset=offset+bsize_i
    if (bl(nbl)%BC(2).eq.-41) offset=offset+bsize_i
    if (bl(nbl)%BC(3).eq.-41) offset=offset+bsize_j
    if (bl(nbl)%BC(4).eq.-41) offset=offset+bsize_j
    if (bl(nbl)%BC(5).eq.-41) offset=offset+bsize_k
    if (bl(nbl)%BC(6).eq.-41) offset=offset+bsize_k
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
    restartR%offset=0_MPI_OFFSET_KIND+offset

    offset=0
    if (bl(nbl)%BC(1).eq.-41) offset=bsize_i
    restartR(1,2)%offset=restartR(1,2)%offset+offset
    if (bl(nbl)%BC(2).eq.-41) offset=offset+bsize_i
    restartR(2,1)%offset=restartR(2,1)%offset+offset
    if (bl(nbl)%BC(3).eq.-41) offset=offset+bsize_j
    restartR(2,2)%offset=restartR(2,2)%offset+offset
    if (.not.is_2d) then
       if (bl(nbl)%BC(4).eq.-41) offset=offset+bsize_j
       restartR(3,1)%offset=restartR(3,1)%offset+offset
       if (bl(nbl)%BC(5).eq.-41) offset=offset+bsize_k
       restartR(3,2)%offset=restartR(3,2)%offset+offset
    endif
 
    ! Creation of derived types for MPI-IO
    ! ====================================

    m1=1 ! faces imin and imax
    do m2=1,2

       ! Creation of the derived type type_mat_nogh that defines the array to be sent
       ! ----------------------------------------------------------------------------

       ! Shape of the array
       shape_nogh= (/ 1, ny, nz /)

       ! Starting coordinates for whole array (i.e. array=subarray)
       start_nogh= (/ 0, 0, 0 /)

       ! Creation of derived type type_mat_nogh
       call MPI_TYPE_CREATE_SUBARRAY(3,shape_nogh,shape_nogh,start_nogh,&
            MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,restartR(m1,m2)%type_mat_nogh,info)
       ! Commit of type_mat_nogh
       call MPI_TYPE_COMMIT(restartR(m1,m2)%type_mat_nogh,info)

       ! Creation of type type_mat_nogh_view to set the view on the file
       ! ---------------------------------------------------------------

       ! Shape of the global array
       shape_glob= (/ 1, ngy, ngz /)

       ! Shape of the local array
       shape_local= (/ 1, ny, nz /)

       ! Starting coordinates for local array
       start_local= (/ 0, coord(2)*ny, coord(3)*nz /)
       if (nx==1) start_local(1) = 0
       if (ny==1) start_local(2) = 0
       if (nz==1) start_local(3) = 0

       ! Creation of derived type type_mat_view
       call MPI_TYPE_CREATE_SUBARRAY(3,shape_glob,shape_local,start_local,&
            MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,restartR(m1,m2)%type_mat_view,info)
       ! Commit of type_mat_view
       call MPI_TYPE_COMMIT(restartR(m1,m2)%type_mat_view,info)

    enddo

    m1=2 ! faces jmin and jmax
    do m2=1,2

       ! Creation of the derived type type_mat_nogh that defines the array to be sent
       ! ----------------------------------------------------------------------------

       ! Shape of the array
       shape_nogh= (/ nx, 1, nz /)

       ! Starting coordinates for whole array (i.e. array=subarray)
       start_nogh= (/ 0, 0, 0 /)

       ! Creation of derived type type_mat_nogh
       call MPI_TYPE_CREATE_SUBARRAY(3,shape_nogh,shape_nogh,start_nogh,&
            MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,restartR(m1,m2)%type_mat_nogh,info)
       ! Commit of type_mat_nogh
       call MPI_TYPE_COMMIT(restartR(m1,m2)%type_mat_nogh,info)

       ! Creation of type type_mat_nogh_view to set the view on the file
       ! ---------------------------------------------------------------

       ! Shape of the global array
       shape_glob= (/ ngx, 1, ngz /)

       ! Shape of the local array
       shape_local= (/ nx, 1, nz /)
        
       ! Starting coordinates for local array
       start_local= (/ coord(1)*nx, 0, coord(3)*nz /)
       if (nx==1) start_local(1) = 0
       if (ny==1) start_local(2) = 0
       if (nz==1) start_local(3) = 0

       ! Creation of derived type type_mat_view
       call MPI_TYPE_CREATE_SUBARRAY(3,shape_glob,shape_local,start_local,&
            MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,restartR(m1,m2)%type_mat_view,info)
       ! Commit of type_mat_view
       call MPI_TYPE_COMMIT(restartR(m1,m2)%type_mat_view,info)

    enddo

    if (.not.is_2d) then
        m1=3 ! faces kmin and kmax
        do m2=1,2
          ! Creation o f the derived type type_mat_nogh that defines the array to be sent
          ! ----------------------------------------------------------------------------

          ! Shape of the array
          shape_nogh= (/ nx, ny, 1 /)

          ! Starting coordinates for whole array (i.e. array=subarray)
          start_nogh= (/ 0, 0, 0 /)

          ! Creation of derived type type_mat_nogh
          call MPI_TYPE_CREATE_SUBARRAY(3,shape_nogh,shape_nogh,start_nogh,&
               MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,restartR(m1,m2)%type_mat_nogh,info)
        
          ! Commit of type_mat_nogh
          call MPI_TYPE_COMMIT(restartR(m1,m2)%type_mat_nogh,info)

          ! Creation of type type_mat_nogh_view to set the view on the file
          ! ---------------------------------------------------------------

          ! Shape of the global array
          shape_glob= (/ ngx, ngy, 1 /)

          ! Shape of the local array
          shape_local= (/ nx, ny, 1 /)

          ! Starting coordinates for local array
          start_local= (/ coord(1)*nx, coord(2)*ny, 0 /)

          if (nx==1) start_local(1) = 0
          if (ny==1) start_local(2) = 0
          if (nz==1) start_local(3) = 0

          ! Creation of derived type type_mat_view
          call MPI_TYPE_CREATE_SUBARRAY(3,shape_glob,shape_local,start_local,&
               MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,restartR(m1,m2)%type_mat_view,info)
          ! Commit of type_mat_view
          call MPI_TYPE_COMMIT(restartR(m1,m2)%type_mat_view,info)
       enddo
    endif
    
  end subroutine init_io_restartR

  !==============================================================================================
  subroutine write_restartR
  !==============================================================================================
    !> author: XG
    !> date: January 2022
    !> Write variables for restartR
  !==============================================================================================
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Local variables
    integer :: n,m1,m2,fh,dim
    integer :: i1,i2,j1,j2,k1,k2
    integer, dimension(MPI_STATUS_SIZE) :: statut
    ! -------------------------------------------------------------------------------------------
    
    if (iproc.eq.0) print *,'Writing restartR.bin ... '

    ! Initialization of IO
    ! ====================
    call init_io_restartR
    
    ! Open file for writing
    ! =====================
    call MPI_FILE_OPEN(COMM_global,'restartR.bin', &
         MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,fh,info)

    ! Write data
    ! ==========
    dim=2
    if (.not.is_2d) dim=3

    do m1=1,dim
       do m2=1,2
          ! Setting the view on the file

          call MPI_FILE_SET_VIEW(fh,restartR(m1,m2)%offset,MPI_DOUBLE_PRECISION, &
               restartR(m1,m2)%type_mat_view,'native',MPI_INFO_NULL,info)
          ! Write each variable if BC is Riemann perturbed
          if (BC_face(m1,m2)%sort==-41) then
             if (m1.eq.1) restartR(m1,m2)%i2=1
             if (m1.eq.2) restartR(m1,m2)%j2=1
             if (m1.eq.3) restartR(m1,m2)%k2=1
             i1=restartR(m1,m2)%i1
             i2=restartR(m1,m2)%i2
             j1=restartR(m1,m2)%j1
             j2=restartR(m1,m2)%j2
             k1=restartR(m1,m2)%k1
             k2=restartR(m1,m2)%k2
             do n=1,nvarR
                call MPI_FILE_WRITE(fh,BC_face(m1,m2)%U0R(i1:i2,j1:j2,k1:k2,n),1,&
                     restartR(m1,m2)%type_mat_nogh,statut,info)
             enddo
          endif
       enddo
    enddo
    
    ! Close file
    ! ==========
    call MPI_FILE_CLOSE(fh,info) 
     
  end subroutine write_restartR

  !==============================================================================================
  subroutine read_restartR
  !==============================================================================================
    !> author: XG
    !> date: January 2022
    !> Read variables for restartR
  !==============================================================================================
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Local variables
    integer :: n,m1,m2,fh,dim
    integer :: i1,i2,j1,j2,k1,k2
    integer, dimension(MPI_STATUS_SIZE) :: statut
    logical :: Rfile_exists
    ! -------------------------------------------------------------------------------------------
    
    if (iproc.eq.0) print *,'Reading restartR.bin ... '

    ! Initialization of IO
    ! ====================
    call init_io_restartR
    
    ! Check if file exists
    ! ====================
    ! if yes, read content
    ! if not, the time-averages U0R are initialized in init_U0R_bc_inlet_outlet (mod_bc_inlet_outlet)
    ! and the file is written at the end of the computation. 
    inquire(file="restartR.bin",exist=Rfile_exists)

    if (Rfile_exists) then

       ! Open file for reading
       ! =====================
       call MPI_FILE_OPEN(COMM_global,'restartR.bin', &
            MPI_MODE_RDONLY,MPI_INFO_NULL,fh,info)


       if (info/=MPI_SUCCESS) then
          call mpistop('MPI_IO: Error in opening restartR.bin',1)
       endif

       ! Read data
       ! =========
       dim=2
       if (.not.is_2d) dim=3

       do m1=1,dim
          do m2=1,2
             ! Setting the view on the file
             call MPI_FILE_SET_VIEW(fh,restartR(m1,m2)%offset,MPI_DOUBLE_PRECISION, &
                  restartR(m1,m2)%type_mat_view,'native',MPI_INFO_NULL,info)
             ! Read each variable if BC is Riemann perturbed
             if (BC_face(m1,m2)%sort==-41) then
                if (m1.eq.1) restartR(m1,m2)%i2=1
                if (m1.eq.2) restartR(m1,m2)%j2=1
                if (m1.eq.3) restartR(m1,m2)%k2=1
                i1=restartR(m1,m2)%i1
                i2=restartR(m1,m2)%i2
                j1=restartR(m1,m2)%j1
                j2=restartR(m1,m2)%j2
                k1=restartR(m1,m2)%k1
                k2=restartR(m1,m2)%k2
                do n=1,nvarR
                   call MPI_FILE_READ(fh,BC_face(m1,m2)%U0R(i1:i2,j1:j2,k1:k2,n),1,&
                        restartR(m1,m2)%type_mat_nogh,statut,info)
                enddo
             endif
          enddo
       enddo
    
       ! Close file
       ! ==========
       call MPI_FILE_CLOSE(fh,info) 

    endif
     
  end subroutine read_restartR

end module mod_io_restartR
