!==============================================================================================
module mod_io_restartTD
!==============================================================================================
  !> author: XG
  !> date: January 2022
  !> Module to read/write restartTD.bin (restart file for T&D' BC)
!==============================================================================================
  use mod_mpi_part
  use mod_flow     ! <- for grid & solution
  use warnstop

  implicit none
  ! -------------------------------------------------------------------------------------------
  ! restartTD TYPE
  ! -------------------------------------------------------------------------------------------
  type restartTD_type
     ! local indices
     integer :: i1,i2,j1,j2,k1,k2
     ! variables for MPI_IO
     integer :: type_mat,type_mat_nogh,type_mat_view
     integer(kind=MPI_OFFSET_KIND) :: offset
  end type restartTD_type
  ! -------------------------------------------------------------------------------------------
  integer, parameter :: nvar=6
  type(restartTD_type), dimension(3,2) :: restartTD,restartTD_o
  ! -------------------------------------------------------------------------------------------
  
contains

  !==============================================================================================
  subroutine init_io_restartTD(restartTD_,itype)
  !==============================================================================================
    !> author: XG
    !> date: January 2022
    !> Initialize the sub-communicator an I/O-type for restartTD outputs
  !==============================================================================================
    use mod_block
    ! use mod_constant ! <- for is_init_2D_3D indicator
    implicit none
    ! -------------------------------------------------------------------------------------------
    integer, intent(in) :: itype ! 0: nz ngz, 1: 1 1 (AB: don't use is_init_2D_3D /!\)
    type(restartTD_type), dimension(3,2) :: restartTD_
    ! -------------------------------------------------------------------------------------------
    ! Local variables
    integer :: nb,nbl,m1,m2,nz_,ngz_
    integer :: nbytes_dbl
    integer :: bsize_i,bsize_j,bsize_k
    integer :: offset
    integer :: npt_TD ! number of points for T&D BC (either ngh or 1)
    integer :: npt_ex ! nb points for T&D BC extended for 7-pt stencil
    integer, dimension(:), allocatable :: offset_bl
    integer, dimension(3) :: shape_wgh,shape_nogh,start_wgh,start_nogh
    integer, dimension(3) :: shape_glob,shape_local,start_local
    ! -------------------------------------------------------------------------------------------

    npt_TD=ngh
    npt_ex=npt_TD+3

    if (is_TamDong_1pt) then
       npt_TD=1
       npt_ex=npt_TD+4
    endif

    ! ! BEUG:
    ! ! If is_init_2D3D, already extruded with extrude_2D_field subroutine
    ! ! If ngz put to 1, beug when writting because writting 2D file
    ! ! -> init_io_restartTD is only called once at the beginning
    ! ! Beug not revealed before because:
    ! ! -> init_io_restartTD called in musica_main
    ! ! -> read_write_info was called in solver in the previous version
    ! ! -> Therefore, "is_init_2D3D" was always .false. when entering here
    ! ! -> AB recently modified this to know is_init_2D3D for grid generation
    if (itype.eq.1) then
    ! if (is_init_2D3D) then
       nz_=1
       ngz_=1
    else
       nz_=nz
       ngz_=ngz
    endif
       
    ! Initialize restartTD local indices
    ! ==================================
    
    ! default initializations
    ! -----------------------
    restartTD_%i1=1
    restartTD_%i2=nx
    restartTD_%j1=1
    restartTD_%j2=ny
    restartTD_%k1=1
    restartTD_%k2=nz_

    ! face(1,1), imin
    restartTD_(1,1)%i1=1
    restartTD_(1,1)%i2=ngh+3      
    ! face(1,2), imax
    restartTD_(1,2)%i1=-2
    restartTD_(1,2)%i2=ngh       
    ! face(2,1), jmin
    restartTD_(2,1)%j1=1
    restartTD_(2,1)%j2=ngh+3       
    ! face(2,2), jmax
    restartTD_(2,2)%j1=-2
    restartTD_(2,2)%j2=ngh
    ! face(3,1), kmin
    restartTD_(3,1)%k1=1
    restartTD_(3,1)%k2=ngh+3       
    ! face(3,2), kmax
    restartTD_(3,2)%k1=-2
    restartTD_(3,2)%k2=ngh
    
    if (is_TamDong_1pt) then
       restartTD_(1,1)%i2=npt_TD+4
       restartTD_(2,1)%j2=npt_TD+4
       restartTD_(3,1)%k2=npt_TD+4
       restartTD_(1,2)%i1=-3
       restartTD_(1,2)%i2=npt_TD
       restartTD_(2,2)%j1=-3
       restartTD_(2,2)%j2=npt_TD
       restartTD_(3,2)%k1=-3
       restartTD_(3,2)%k2=npt_TD
    endif


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
    bsize_i=npt_ex*ngy*ngz_*nvar*nbytes_dbl
    ! j-direction (jmin,jmax)
    bsize_j=ngx*npt_ex*ngz_*nvar*nbytes_dbl
    ! k-direction (kmin,kmax)
    bsize_k=ngx*ngy*npt_ex*nvar*nbytes_dbl

    ! First pass to define global offset per block
    ! --------------------------------------------
    offset=0
    ! add all active face for current block
    if ((bl(nbl)%BC(1).eq.-1).or.(bl(nbl)%BC(1).eq.-2)) offset=offset+bsize_i
    if ((bl(nbl)%BC(2).eq.-1).or.(bl(nbl)%BC(2).eq.-2)) offset=offset+bsize_i
    if ((bl(nbl)%BC(3).eq.-1).or.(bl(nbl)%BC(3).eq.-2)) offset=offset+bsize_j
    if ((bl(nbl)%BC(4).eq.-1).or.(bl(nbl)%BC(4).eq.-2)) offset=offset+bsize_j
    if ((bl(nbl)%BC(5).eq.-1).or.(bl(nbl)%BC(5).eq.-2)) offset=offset+bsize_k
    if ((bl(nbl)%BC(6).eq.-1).or.(bl(nbl)%BC(6).eq.-2)) offset=offset+bsize_k
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
    restartTD_%offset=0_MPI_OFFSET_KIND+offset

    offset=0
    if ((bl(nbl)%BC(1).eq.-1).or.(bl(nbl)%BC(1).eq.-2)) offset=bsize_i
    restartTD_(1,2)%offset=restartTD_(1,2)%offset+offset
    if ((bl(nbl)%BC(2).eq.-1).or.(bl(nbl)%BC(2).eq.-2)) offset=offset+bsize_i
    restartTD_(2,1)%offset=restartTD_(2,1)%offset+offset
    if ((bl(nbl)%BC(3).eq.-1).or.(bl(nbl)%BC(3).eq.-2)) offset=offset+bsize_j
    restartTD_(2,2)%offset=restartTD_(2,2)%offset+offset
    if (is_TamDong3D) then
       if ((bl(nbl)%BC(4).eq.-1).or.(bl(nbl)%BC(4).eq.-2)) offset=offset+bsize_j
       restartTD_(3,1)%offset=restartTD_(3,1)%offset+offset
       if ((bl(nbl)%BC(5).eq.-1).or.(bl(nbl)%BC(5).eq.-2)) offset=offset+bsize_k
       restartTD_(3,2)%offset=restartTD_(3,2)%offset+offset
    endif
 
    ! Creation of derived types for MPI-IO
    ! ====================================

    m1=1 ! faces imin and imax
    do m2=1,2

       ! Creation of the derived type type_local that defines the array without ghost cells
       ! ----------------------------------------------------------------------------------

       ! Shape of the array with ghost cells
       shape_wgh= (/ npt_ex, ny+2*ngh, nz_+2*ngh /)
       if (nz_.eq.1) shape_wgh= (/ npt_ex, ny+2*ngh, 1 /)

       ! Shape of the array without ghost cells
       shape_nogh= (/ npt_ex, ny, nz_ /)

       ! Starting coordinates for array without ghost cells (w.r.t. global array)
       start_wgh= (/ 0, ngh, ngh /)
       if (nz_.eq.1) start_wgh = (/ 0, ngh, 0 /)

       ! Starting coordinates for whole array (i.e. array=subarray)
       start_nogh= (/ 0, 0, 0 /)

       ! Creation of derived type type_mat_nogh
       ! (This allows to remove ghost-cells from initial array)
       call MPI_TYPE_CREATE_SUBARRAY(3,shape_nogh,shape_nogh,start_nogh,&
            MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,restartTD_(m1,m2)%type_mat_nogh,info)
       ! Commit of type_mat_nogh
       call MPI_TYPE_COMMIT(restartTD_(m1,m2)%type_mat_nogh,info)

       ! Creation of derived type type_mat
       call MPI_TYPE_CREATE_SUBARRAY(3,shape_wgh,shape_nogh,start_wgh,&
            MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,restartTD_(m1,m2)%type_mat,info)
       ! Commit of type_mat
       call MPI_TYPE_COMMIT(restartTD_(m1,m2)%type_mat,info)

       ! Creation of type type_mat_nogh_view to set the view on the file
       ! ---------------------------------------------------------------

       ! Shape of the global array
       shape_glob= (/ npt_ex, ngy, ngz_ /)

       ! Shape of the local array
       shape_local= (/ npt_ex, ny, nz_ /)

       ! Starting coordinates for local array
       start_local= (/ 0, coord(2)*ny, coord(3)*nz_ /)
       if (nx ==1) start_local(1) = 0
       if (ny ==1) start_local(2) = 0
       if (nz_==1) start_local(3) = 0

       ! Creation of derived type type_mat_view
       call MPI_TYPE_CREATE_SUBARRAY(3,shape_glob,shape_local,start_local,&
            MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,restartTD_(m1,m2)%type_mat_view,info)
       ! Commit of type_mat_view
       call MPI_TYPE_COMMIT(restartTD_(m1,m2)%type_mat_view,info)

    enddo

    m1=2 ! faces jmin and jmax
    do m2=1,2

       ! Creation of the derived type type_local that defines the array without ghost cells
       ! ----------------------------------------------------------------------------------

       ! Shape of the array with ghost cells
       shape_wgh= (/ nx+2*ngh, npt_ex, nz_+2*ngh /)
       if (nz_.eq.1) shape_wgh= (/ nx+2*ngh, npt_ex, 1 /)

       ! Shape of the array without ghost cells
       shape_nogh= (/ nx, npt_ex, nz_ /)

       ! Starting coordinates for array without ghost cells (w.r.t. global array)
       start_wgh= (/ ngh, 0, ngh /)
       if (nz_.eq.1) start_wgh = (/ ngh, 0, 0 /)

       ! Starting coordinates for whole array (i.e. array=subarray)
       start_nogh= (/ 0, 0, 0 /)

       ! Creation of derived type type_mat_nogh
       ! (This allows to remove ghost-cells from initial array)
       call MPI_TYPE_CREATE_SUBARRAY(3,shape_nogh,shape_nogh,start_nogh,&
            MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,restartTD_(m1,m2)%type_mat_nogh,info)
       ! Commit of type_mat_nogh
       call MPI_TYPE_COMMIT(restartTD_(m1,m2)%type_mat_nogh,info)

       ! Creation of derived type type_mat
       call MPI_TYPE_CREATE_SUBARRAY(3,shape_wgh,shape_nogh,start_wgh,&
            MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,restartTD_(m1,m2)%type_mat,info)
       ! Commit of type_mat
       call MPI_TYPE_COMMIT(restartTD_(m1,m2)%type_mat,info)

       ! Creation of type type_mat_nogh_view to set the view on the file
       ! ---------------------------------------------------------------

       ! Shape of the global array
       shape_glob= (/ ngx, npt_ex, ngz_ /)

       ! Shape of the local array
       shape_local= (/ nx, npt_ex, nz_ /)

       ! Starting coordinates for local array
       start_local= (/ coord(1)*nx, 0, coord(3)*nz_ /)
       if (nx ==1) start_local(1) = 0
       if (ny ==1) start_local(2) = 0
       if (nz_==1) start_local(3) = 0

       ! Creation of derived type type_mat_view
       call MPI_TYPE_CREATE_SUBARRAY(3,shape_glob,shape_local,start_local,&
            MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,restartTD_(m1,m2)%type_mat_view,info)
       ! Commit of type_mat_view
       call MPI_TYPE_COMMIT(restartTD_(m1,m2)%type_mat_view,info)

    enddo

    if (is_TamDong3D) then
       m1=3 ! faces kmin and kmax
       do m2=1,2

          ! Creation of the derived type type_local that defines the array without ghost cells
          ! ----------------------------------------------------------------------------------

          ! Shape of the array with ghost cells
          shape_wgh= (/ nx+2*ngh, ny+2*ngh, npt_ex /)

          ! Shape of the array without ghost cells
          shape_nogh= (/ nx, ny, npt_ex /)

          ! Starting coordinates for array without ghost cells (w.r.t. global array)
          start_wgh= (/ ngh, ngh, 0 /)

          ! Starting coordinates for whole array (i.e. array=subarray)
          start_nogh= (/ 0, 0, 0 /)

          ! Creation of derived type type_mat_nogh
          ! (This allows to remove ghost-cells from initial array)
          call MPI_TYPE_CREATE_SUBARRAY(3,shape_nogh,shape_nogh,start_nogh,&
               MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,restartTD_(m1,m2)%type_mat_nogh,info)
          ! Commit of type_mat_nogh
          call MPI_TYPE_COMMIT(restartTD_(m1,m2)%type_mat_nogh,info)

          ! Creation of derived type type_mat
          call MPI_TYPE_CREATE_SUBARRAY(3,shape_wgh,shape_nogh,start_wgh,&
               MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,restartTD_(m1,m2)%type_mat,info)
          ! Commit of type_mat
          call MPI_TYPE_COMMIT(restartTD_(m1,m2)%type_mat,info)

          ! Creation of type type_mat_nogh_view to set the view on the file
          ! ---------------------------------------------------------------

          ! Shape of the global array
          shape_glob= (/ ngx, ngy, npt_ex /)

          ! Shape of the local array
          shape_local= (/ nx, ny, npt_ex /)

          ! Starting coordinates for local array
          start_local= (/ coord(1)*nx, coord(2)*ny, 0 /)
          if (nx ==1) start_local(1) = 0
          if (ny ==1) start_local(2) = 0
          if (nz_==1) start_local(3) = 0

          ! Creation of derived type type_mat_view
          call MPI_TYPE_CREATE_SUBARRAY(3,shape_glob,shape_local,start_local,&
               MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,restartTD_(m1,m2)%type_mat_view,info)
          ! Commit of type_mat_view
          call MPI_TYPE_COMMIT(restartTD_(m1,m2)%type_mat_view,info)

       enddo
    endif
 
    if (iproc==0) print *,'init I/O T&D OK'

  end subroutine init_io_restartTD

  !==============================================================================================
  subroutine write_restartTD(filename)
  !==============================================================================================
    !> author: XG
    !> date: January 2022
    !> Write variables for restartTD
  !==============================================================================================
    implicit none
    ! -------------------------------------------------------------------------------------------
    character(len=*), intent(in) :: filename
    ! -------------------------------------------------------------------------------------------
    ! Local variables
    integer :: n,m1,m2,fh,dim
    integer :: i1,i2,j1,j2,k1,k2
    integer, dimension(MPI_STATUS_SIZE) :: statut
    ! -------------------------------------------------------------------------------------------
    
    if (iproc.eq.0) print *,'Writing restartTD.bin ... '
    
    ! Open file for writing
    ! =====================
    call MPI_FILE_OPEN(COMM_global,trim(filename), &
         MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,fh,info)

    ! Write data
    ! ==========
    dim=2
    if (is_TamDong3D) dim=3

    do m1=1,dim
       do m2=1,2
          ! Setting the view on the file
          call MPI_FILE_SET_VIEW(fh,restartTD(m1,m2)%offset,MPI_DOUBLE_PRECISION, &
               restartTD(m1,m2)%type_mat_view,'native',MPI_INFO_NULL,info)
          ! Write each variable if BC is T&D
          if (is_bc_TD(m1,m2)) then
             i1=restartTD(m1,m2)%i1
             i2=restartTD(m1,m2)%i2
             j1=restartTD(m1,m2)%j1
             j2=restartTD(m1,m2)%j2
             k1=restartTD(m1,m2)%k1
             k2=restartTD(m1,m2)%k2
             do n=1,nvar
                call MPI_FILE_WRITE(fh,BC_face(m1,m2)%U0(i1:i2,j1:j2,k1:k2,n),1,&
                     restartTD(m1,m2)%type_mat_nogh,statut,info)
             enddo
          endif
       enddo
    enddo
    
    ! Close file
    ! ==========
    call MPI_FILE_CLOSE(fh,info) 
     
  end subroutine write_restartTD

  !==============================================================================================
  subroutine read_restartTD(filename,restartTD_)
  !==============================================================================================
    !> author: XG
    !> date: January 2022
    !> Read variables for restartTD
  !==============================================================================================
    implicit none
    ! -------------------------------------------------------------------------------------------
    character(len=*), intent(in) :: filename
    type(restartTD_type), dimension(3,2) :: restartTD_
    ! -------------------------------------------------------------------------------------------
    ! Local variables
    integer :: n,m1,m2,fh,dim
    integer :: i1,i2,j1,j2,k1,k2
    integer, dimension(MPI_STATUS_SIZE) :: statut
    ! -------------------------------------------------------------------------------------------
    
    if (iproc.eq.0) print *,'Reading restartTD.bin ... '
    
    ! Open file for reading
    ! =====================
    call MPI_FILE_OPEN(COMM_global,trim(filename), &
         MPI_MODE_RDONLY,MPI_INFO_NULL,fh,info)
    
    if (info/=MPI_SUCCESS) then
       call mpistop('MPI_IO: Error in opening restartTD.bin',1)
    endif

    ! Read data
    ! =========
    dim=2
    if (is_TamDong3D) dim=3

    do m1=1,dim
       do m2=1,2
          ! Setting the view on the file
          call MPI_FILE_SET_VIEW(fh,restartTD_(m1,m2)%offset,MPI_DOUBLE_PRECISION, &
               restartTD_(m1,m2)%type_mat_view,'native',MPI_INFO_NULL,info)
          ! Read each variable if BC is T&D
          if (is_bc_TD(m1,m2)) then
             i1=restartTD_(m1,m2)%i1
             i2=restartTD_(m1,m2)%i2
             j1=restartTD_(m1,m2)%j1
             j2=restartTD_(m1,m2)%j2
             k1=restartTD_(m1,m2)%k1
             k2=restartTD_(m1,m2)%k2
             do n=1,nvar
                call MPI_FILE_READ(fh,BC_face(m1,m2)%U0(i1:i2,j1:j2,k1:k2,n),1,&
                     restartTD_(m1,m2)%type_mat_nogh,statut,info)
             enddo
          endif
       enddo
    enddo
  !  if (iproc==68) print*,'vdalfasda'
  !  if (iproc==68) print *,BC_face(1,1)%U0(1,1,1,5),'after read'
   ! call mpistop('',0)
    
    ! Close file
    ! ==========
    call MPI_FILE_CLOSE(fh,info) 
     
  end subroutine read_restartTD

  !==============================================================================================
  subroutine free_restartTD(restartTD_)
  !==============================================================================================
    !> author: XG
    !> date: May 2022
    !> Free MPI-IO types for restartTD
  !==============================================================================================
    implicit none
    ! -------------------------------------------------------------------------------------------
    type(restartTD_type), dimension(3,2) :: restartTD_
    ! -------------------------------------------------------------------------------------------
    ! Local variables
    integer :: m1,m2,dim
    ! -------------------------------------------------------------------------------------------

    ! select dimension
    ! ----------------
    dim=2
    if (is_TamDong3D) dim=3

    ! Free types for MPI-IO of 2D restartTD
    ! -------------------------------------
    do m1=1,dim
       do m2=1,2
          call MPI_TYPE_FREE(restartTD_(m1,m2)%type_mat_nogh,info)
          call MPI_TYPE_FREE(restartTD_(m1,m2)%type_mat,info)
          call MPI_TYPE_FREE(restartTD_(m1,m2)%type_mat_view,info)
       enddo
    enddo

  end subroutine free_restartTD

end module mod_io_restartTD
