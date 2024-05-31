!===============================================================================
module mod_flow_o
!===============================================================================
  !> Module to handle old flow field, which serves as DONOR for interpolation 
  !> Rq: put old grid and restart in a directory named "old"
!===============================================================================
  use mpi
  use mod_ngh
  implicit none
  ! ----------------------------------------------------------------------------
  ! Global donor grid (old grid)
  integer  :: ngx_o,ngy_o,ngz_o
  real(wp) :: deltaz_o
  real(wp), dimension(:), allocatable :: xg_o,yg_o,zg_o
  real(wp), dimension(:,:), allocatable :: xgc_o,ygc_o
  ! ----------------------------------------------------------------------------
  ! Local donor grid (old grid)
  integer  :: nx_o,ny_o,nz_o
  integer :: nx_or,ny_or,nz_or
  real(wp), dimension(:), allocatable :: x_o,y_o,z_o
  real(wp), dimension(:,:), allocatable :: xc_o,yc_o 
  ! ----------------------------------------------------------------------------
  ! Local donor field (old grid)
  real(wp), dimension(:,:,:), allocatable :: rho_o,rhou_o,rhov_o,rhow_o,rhoe_o
  ! ----------------------------------------------------------------------------
  ! MPI types for communications
  integer :: type_facex_o,type_facey_o,type_facez_o,type_linez_o
  ! type face (ngh layers of cells) for each direction [flux_euler & cie]
  integer :: type_faceW_o,type_faceE_o
  integer :: type_faceS_o,type_faceN_o
  integer :: type_faceF_o,type_faceB_o
  ! type face+edges (ngh layers of cells) for each direction [NOT USED ANYMORE]
  integer :: type_e_faceW,type_e_faceE
  integer :: type_e_faceS,type_e_faceN
  integer :: type_e_faceF,type_e_faceB
  ! ----------------------------------------------------------------------------
  !integer, dimension(MPI_STATUS_SIZE) :: status
  
  ! indices to start communications in the direction [p]arallel to the face
  integer :: ipW_o,ipE_o,ipS_o,ipN_o
  ! indices to start communications in the direction [n]ormal to the face
  integer :: inW_o,inE_o,inS_o,inN_o
  !--------for the cases of adjoint block interfaces
  ! indices to start send routines in the direction [n]ormal to the face
  integer :: inWs_o,inEs_o,inSs_o,inNs_o
  ! indices to start receive routines in the direction [n]ormal to the face
  integer :: inWr_o,inEr_o,inSr_o,inNr_o
  
  ! indices to start communications in the direction [p]arallel for face+edges
  integer :: ipW_e,ipE_e,ipS_e,ipN_e
  ! indices to start communications in the direction [n]ormal for face+edges
  integer :: inW_e,inE_e,inS_e,inN_e
  !-----------------------------------------------------------------------------
!!$  ! indices to start communications in the direction [p]arallel for face+edges
!!$  integer :: ipW_e,ipE_e,ipS_e,ipN_e
!!$  ! indices to start communications in the direction [n]ormal for face+edges
!!$  integer :: inW_e,inE_e,inS_e,inN_e
  ! Communications of unknown vector along faces
  ! ============================================
  ! 3D version
  integer, parameter :: size_3d=12*5 ! 2(send/recv)*6(neighbors)*5(variables)
  integer, dimension(size_3d) :: request_3d
  integer, dimension(MPI_STATUS_SIZE,size_3d) :: status_3d
  ! 2D version
  integer, parameter :: size_2d=8*4 ! 2(send/recv)*4(neighbors)*4(variables)
  integer, dimension(size_2d) :: request_2d
  integer, dimension(MPI_STATUS_SIZE,size_2d) :: status_2d

contains

  !=============================================================================
  subroutine read_donor_field(dir_old)
  !=============================================================================
    !> author: XG
    !> date: April 2021
    !> read donor field on old grid to be interpolated on new grid
  !=============================================================================
    use mod_mpi
    use mod_block
    use mod_io
    use mod_utils
    use warnstop
    implicit none
    ! --------------------------------------------------------------------------
    character(len=*) :: dir_old
    ! --------------------------------------------------------------------------
    integer :: i,j
    
    ! Read old grid
    ! =============
    !call read_grid_o
    call read_grid_ex_o(dir_old)

    ! MPI Partitioning
    ! ==================
    call mpi_part_o
    ! create MPI types for communications
    call mpi_types_comm_edges_o
    !call mpi_types_comm_o

    !print *,iproc,z_o
    
    !call mpistop('check z',0)
    
    
    ! Allocate old field
    ! ==================
    allocate( rho_o(1-ngh:nx_o+ngh,1-ngh:ny_o+ngh,1-ngh:nz_o+ngh))
    allocate(rhou_o(1-ngh:nx_o+ngh,1-ngh:ny_o+ngh,1-ngh:nz_o+ngh))
    allocate(rhov_o(1-ngh:nx_o+ngh,1-ngh:ny_o+ngh,1-ngh:nz_o+ngh))
    allocate(rhow_o(1-ngh:nx_o+ngh,1-ngh:ny_o+ngh,1-ngh:nz_o+ngh))
    allocate(rhoe_o(1-ngh:nx_o+ngh,1-ngh:ny_o+ngh,1-ngh:nz_o+ngh))
    
!!$    rho_o=rho_ref
!!$    rhou_o=0.0_wp
!!$    rhov_o=0.0_wp
!!$    rhow_o=0.0_wp
!!$    rhoe_o=p_ref/0.4
    
    ! Create MPI_IO types and field structure for old field [field_o]
    ! =====================================================
    allocate(field_o(nbloc))
    field_o(nob(iproc))%MPI_COMM = COMM_intrablock
    call mod_io_init(ngx_o,ngy_o,ngz_o,nx_o,ny_o,nz_o,nx_or,ny_or,nz_or, &
         ngh,3,coord,is_IOtec_read,is_IOtec_write,field_o(nob(iproc)))
    
    ! Read old field
    ! ==============
    binfile = trim(dir_old)//'/restart_bl'//trim(numchar(nob(iproc)))//'.bin'
    call read_vol_o(binfile, rho_o(1:nx_o,1:ny_o,1:nz_o) &
                           ,rhou_o(1:nx_o,1:ny_o,1:nz_o) &
                           ,rhov_o(1:nx_o,1:ny_o,1:nz_o) &
                           ,rhow_o(1:nx_o,1:ny_o,1:nz_o) &
                           ,rhoe_o(1:nx_o,1:ny_o,1:nz_o) )

    ! Communicate old field in ghost cells
    ! ====================================

    if (is_2d) then
       call commun2d_edges(rho_o,1)
       call commun2d_edges(rhou_o,2)
       call commun2d_edges(rhov_o,3)
       call commun2d_edges(rhoe_o,4)
       call MPI_WAITALL(size_2d,request_2d,status_2d,info)
       ! second time for edge order
       call commun2d_edges(rho_o,1)
       call commun2d_edges(rhou_o,2)
       call commun2d_edges(rhov_o,3)
       call commun2d_edges(rhoe_o,4)
       call MPI_WAITALL(size_2d,request_2d,status_2d,info)
    else
       call commun3d_edges(rho_o,1)
       call commun3d_edges(rhou_o,2)
       call commun3d_edges(rhov_o,3)
       call commun3d_edges(rhow_o,4)
       call commun3d_edges(rhoe_o,5)
       call MPI_WAITALL(size_3d,request_3d,status_3d,info)
       ! second time for edge order
       call commun3d_edges(rho_o,1)
       call commun3d_edges(rhou_o,2)
       call commun3d_edges(rhov_o,3)
       call commun3d_edges(rhow_o,4)
       call commun3d_edges(rhoe_o,5)
       call MPI_WAITALL(size_3d,request_3d,status_3d,info)
    endif
    
!!$    if (iproc.eq.7) then
!!$       open(194,file='sol_proc8_ex.bin',form='unformatted',status='unknown')
!!$       rewind(194)
!!$       write(194) nx_o+2*ngh
!!$       write(194) ny_o+2*ngh
!!$       write(194) ((xc_o(i,j),i=1-ngh,nx_o+ngh),j=1-ngh,ny_o+ngh)
!!$       write(194) ((yc_o(i,j),i=1-ngh,nx_o+ngh),j=1-ngh,ny_o+ngh)
!!$       write(194) ((rho_o(i,j,1),i=1-ngh,nx_o+ngh),j=1-ngh,ny_o+ngh)
!!$       write(194) ((rhou_o(i,j,1),i=1-ngh,nx_o+ngh),j=1-ngh,ny_o+ngh)
!!$       close(194)
!!$    endif
    
!!$    if (iproc.eq.8) then
!!$       open(194,file='sol_proc9_ex.bin',form='unformatted',status='unknown')
!!$       rewind(194)
!!$       write(194) nx_o+2*ngh
!!$       write(194) ny_o+2*ngh
!!$       write(194) ((xc_o(i,j),i=1-ngh,nx_o+ngh),j=1-ngh,ny_o+ngh)
!!$       write(194) ((yc_o(i,j),i=1-ngh,nx_o+ngh),j=1-ngh,ny_o+ngh)
!!$       write(194) ((rho_o(i,j,1),i=1-ngh,nx_o+ngh),j=1-ngh,ny_o+ngh)
!!$       write(194) ((rhou_o(i,j,1),i=1-ngh,nx_o+ngh),j=1-ngh,ny_o+ngh)
!!$       close(194)
!!$    endif

    !call mpistop('check',0)
    
!!$
!!$    if (is_2d) then
!!$       call commun2d_o(rho_o,1)
!!$       call commun2d_o(rhou_o,2)
!!$       call commun2d_o(rhov_o,3)
!!$       call commun2d_o(rhoe_o,4)
!!$       call MPI_WAITALL(size_2d,request_2d,status_2d,info)
!!$    else
!!$       call commun3d_o(rho_o,1)
!!$       call commun3d_o(rhou_o,2)
!!$       call commun3d_o(rhov_o,3)
!!$       call commun3d_o(rhow_o,4)
!!$       call commun3d_o(rhoe_o,5)
!!$       call MPI_WAITALL(size_3d,request_3d,status_3d,info)
!!$    endif

    !call mpistop('pause',0)

!!$    if (is_2d) then
!!$       ! communicate faces
!!$       call communic2d_o(rho_o)
!!$       call communic2d_o(rhou_o)
!!$       call communic2d_o(rhov_o)
!!$       call communic2d_o(rhow_o)
!!$       call communic2d_o(rhoe_o)
!!$    else
!!$       ! communicate faces
!!$       call communic3d_o(rho_o)
!!$       call communic3d_o(rhou_o)
!!$       call communic3d_o(rhov_o)
!!$       call communic3d_o(rhow_o)
!!$       call communic3d_o(rhoe_o)
!!$    endif
!!$    ! communicate edges along z
!!$    call communic_edgesz_o(rho_o)
!!$    call communic_edgesz_o(rhou_o)
!!$    call communic_edgesz_o(rhov_o)
!!$    call communic_edgesz_o(rhow_o)
!!$    call communic_edgesz_o(rhoe_o)
    
  end subroutine read_donor_field

  !=============================================================================
  subroutine read_grid_ex_o(dir_old)
  !=============================================================================
    !> Read donor grid
  !=============================================================================
    use mod_mpi
    use mod_grid ! <~ for is_curv
    use mod_utils
    use warnstop
    implicit none
    ! --------------------------------------------------------------------------
    character(len=*) :: dir_old
    ! --------------------------------------------------------------------------
    integer :: i,j,k
    integer :: ngx_ex,ngy_ex,ngz_ex
    real(wp) :: delta_z
    ! --------------------------------------------------------------------------

    ! Read old grid
    ! =============
    open(194,file=trim(dir_old)//'/grid_bl'//trim(numchar(nob(iproc)))//'.bin',form='unformatted',status='unknown')
    rewind(194)
    ! read dims for block
    read(194) ngx_o
    read(194) ngy_o
    read(194) ngz_o
    ngz_o=300
    if (iproc==iproc_leader(nob(iproc))) then
       print *,'Nx x Ny x Nz old:', ngx_o,ngy_o,ngz_o,' for block',nob(iproc)
       print *,'Nx x Ny x Nz new:', ngx,ngy,ngz,' for block',nob(iproc)
    endif
    ! read x,y for block
    if (is_curv) then
       allocate(xgc_o(1-ngh:ngx_o+ngh,1-ngh:ngy_o+ngh))
       allocate(ygc_o(1-ngh:ngx_o+ngh,1-ngh:ngy_o+ngh))
       xgc_o=0.0_wp
       ygc_o=0.0_wp
       read(194) ((xgc_o(i,j),i=1,ngx_o),j=1,ngy_o)
       read(194) ((ygc_o(i,j),i=1,ngx_o),j=1,ngy_o)
    else
       allocate(xg_o(1-ngh:ngx_o+ngh),yg_o(1-ngh:ngy_o+ngh))  
       read(194) (xg_o(i),i=1,ngx_o)
       read(194) (yg_o(j),j=1,ngy_o)
    endif
    ! read z for block
    allocate(zg_o(1-ngh:ngz_o+ngh))
    read(194) (zg_o(k),k=1,ngz_o)
    close(194)

    ! Compute deltaz in old grid
    ! ==========================
    deltaz_o=deltaz
    do k=2,ngz_o
       delta_z=zg_o(k)-zg_o(k-1)
       if ((iproc==0).and.(abs(delta_z-deltaz_o)>1.0e-16_wp)) &
            print *,'deltaz_o different from deltaz',delta_z,deltaz_o
       deltaz_o=delta_z
    enddo

         ! old grid
     deltaz_o=deltaz*2.0_wp
     zg_o(ngz_o/2)=0.0_wp
     do k=ngz_o/2+1,ngz_o
        zg_o(k)=zg_o(k-1)+deltaz_o
     enddo
     do k=ngz_o/2-1,1,-1
        zg_o(k)=zg_o(k+1)-deltaz_o
     enddo

     if ((iproc==0).and.(abs(delta_z-deltaz_o)>1.0e-16_wp)) &
          print *,'deltaz_o different from deltaz',delta_z,deltaz_o

    ! Read old grid_ex (with ghost points)
    ! ================
    open(194,file=trim(dir_old)//'/grid_bl'//trim(numchar(nob(iproc)))//'_ex.bin',form='unformatted',status='unknown')
    rewind(194)
    ! read dims for block
    read(194) ngx_ex
    read(194) ngy_ex
    read(194) ngz_ex
    if (iproc==iproc_leader(nob(iproc))) then
       if (ngx_ex.ne.ngx_o+2*ngh) then
          print *,'Dimension mismatch in grid and grid_ex: ngx_ex=',ngx_ex,' ngx_o+2*ngh=',ngx_o+2*ngh
          call mpistop('check please',0)
       endif
       if (ngy_ex.ne.ngy_o+2*ngh) then
          print *,'Dimension mismatch in grid and grid_ex: ngy_ex=',ngy_ex,' ngy_o+2*ngh=',ngy_o+2*ngh
          call mpistop('check please',0)
       endif
    endif
    ! read x,y for block
    if (is_curv) then
       !allocate(xgc_o(1-ngh:ngx_o+ngh,1-ngh:ngy_o+ngh))
       !allocate(ygc_o(1-ngh:ngx_o+ngh,1-ngh:ngy_o+ngh))
       xgc_o=0.0_wp
       ygc_o=0.0_wp
       read(194) ((xgc_o(i,j),i=1-ngh,ngx_o+ngh),j=1-ngh,ngy_o+ngh)
       read(194) ((ygc_o(i,j),i=1-ngh,ngx_o+ngh),j=1-ngh,ngy_o+ngh)
    else
       !allocate(xg_o(1-ngh:ngx_o+ngh),yg_o(1-ngh:ngy_o+ngh))  
       read(194) (xg_o(i),i=1-ngh,ngx_o+ngh)
       read(194) (yg_o(j),j=1-ngh,ngy_o+ngh)
    endif
    ! read z for block
    !allocate(zg_o(1-ngh:ngz_o+ngh))
    !read(194) (zg_o(k),k=1,ngz_o)

  end subroutine read_grid_ex_o

  !=============================================================================
  subroutine read_grid_o
  !=============================================================================
    !> Read donor grid
  !=============================================================================
    use mod_mpi
    use mod_grid ! <~ for is_curv
    use mod_utils
    implicit none
    ! --------------------------------------------------------------------------
    integer :: i,j,k
    ! --------------------------------------------------------------------------

    ! Read old grid
    ! =============
    open(194,file='old/grid_bl'//trim(numchar(nob(iproc)))//'.bin',form='unformatted',status='unknown')
    rewind(194)
    ! read dims for block
    read(194) ngx_o
    read(194) ngy_o
    read(194) ngz_o
    if (iproc==iproc_leader(nob(iproc))) print *,'Nx x Ny x Nz old:', ngx_o,ngy_o,ngz_o
    ! read x,y for block
    if (is_curv) then
       allocate(xgc_o(1-ngh:ngx_o+ngh,1-ngh:ngy_o+ngh))
       allocate(ygc_o(1-ngh:ngx_o+ngh,1-ngh:ngy_o+ngh))
       xgc_o=0.0_wp
       ygc_o=0.0_wp
       read(194) ((xgc_o(i,j),i=1,ngx_o),j=1,ngy_o)
       read(194) ((ygc_o(i,j),i=1,ngx_o),j=1,ngy_o)
    else
       allocate(xg_o(1-ngh:ngx_o+ngh),yg_o(1-ngh:ngy_o+ngh))  
       read(194) (xg_o(i),i=1,ngx_o)
       read(194) (yg_o(j),j=1,ngy_o)
    endif
    ! read z for block
    allocate(zg_o(1-ngh:ngz_o+ngh))
    read(194) (zg_o(k),k=1,ngz_o)

    ! Fill ghost points and apply periodicity for the grid
    ! ====================================================
    call grid_comm_o

  end subroutine read_grid_o

  !===============================================================================
  subroutine grid_comm_o
    !===============================================================================
    !> Communicate global grids in interblock communicator
    !> and then gather in intrablock communicator (version old field)
    !===============================================================================
    !!use mod_mpi_part
    use mod_block ! <~ for nbloc POUSSIF TO BE CHANGED
    use mod_grid  ! <~ needed for is_curv & is_2D TO BE CHANGED
                  ! <~ also for Lx,Ly,Lz Periodicity TO BE CHANGED
    use mod_grid_comm
    use mod_constant ! <~ needed for CHAN  TO BE CHANGED
    use mod_bc_periodicity ! <~ needed for Lxp,Lyp,Lzp  TO BE CHANGED
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    ! Definition of MPI types for grid communications
    ! ===============================================
    call grid_comm_type_cart
    if (is_curv) call grid_comm_type_curv(ngx_o,ngy_o)
    call MPI_BARRIER(COMM_global,info)

    ! Interblock communications of global grid
    ! ========================================
    if (iproc.eq.iproc_leader(nob(iproc))) then
       if (is_curv) then
          call grid_comm_interblock_curv(xgc_o,ngx_o,ngy_o)
          call grid_comm_interblock_curv(ygc_o,ngx_o,ngy_o)
       else
          call grid_comm_interblock_cart(xg_o,ngx_o,1)
          call grid_comm_interblock_cart(yg_o,ngy_o,2)
       endif

       if (.not.(is_2D)) then
          call grid_comm_interblock_cart(zg_o,ngz_o,3)
       endif

       ! correction for periodicity ????????????? TO BE CHANGED
       if (PHILL) then   
          if (iproc==iproc_leader(1)) then
             xgc_o(-4:0,:)= xgc_o(-4:0,:)-Lxp
          endif
          if (iproc==iproc_leader(nbloc)) then
             xgc_o(ngx_o+1:ngx_o+5,:)=xgc_o(ngx_o+1:ngx_o+5,:)+Lxp
          endif
       endif

!!$       ! correction for periodicity along x
!!$       if (periods(1)==.true.) then
!!$          xg_o(-4:0)= xg_o(-4:0)-Lx!-deltax
!!$          xg_o(ngx_o+1:ngx_o+5)= xg_o(ngx_o+1:ngx_o+5)+Lx!+deltax
!!$       endif
!!$       ! correction for periodicity along z
!!$       if (periods(2)==.true.) then
!!$          yg_o(-4:0)= yg_o(-4:0)-Ly!-deltay
!!$          yg_o(ngy_o+1:ngy_o+5)= yg_o(ngy_o+1:ngy_o+5)+Ly!+deltay
!!$       endif
!!$       ! correction for periodicity along z
!!$       if (periods(3)==.true.) then
!!$          zg_o(-4:0)= zg_o(-4:0)-Lz
!!$          zg_o(ngz_o+1:ngz_o+5)= zg_o(ngz_o+1:ngz_o+5)+Lz
!!$       endif

    endif

    call MPI_BARRIER(COMM_global,info)

    ! Intrablock communications of global grid
    ! ========================================
    if (is_curv) then
       call grid_comm_intrablock_curv(xgc_o,ngx_o,ngy_o)
       call grid_comm_intrablock_curv(ygc_o,ngx_o,ngy_o)
    else
       call grid_comm_intrablock_cart(xg_o,ngx_o)
       call grid_comm_intrablock_cart(yg_o,ngy_o)
    endif

    if (.not.(is_2D)) then
       call grid_comm_intrablock_cart(zg_o,ngz_o)
    endif

    ! Free MPI types for grid comm
    ! ============================
    call free_grid_comm_type_cart
    if (is_curv) call free_grid_comm_type_curv

  end subroutine grid_comm_o

  !=============================================================================
  subroutine mpi_part_o
  !=============================================================================
    !> MPI partitioning of donor grid (old grid)
  !=============================================================================
    use mod_mpi_part
    use mod_grid
    use warnstop
    implicit none
    ! --------------------------------------------------------------------------
    integer :: i,j,k
    ! --------------------------------------------------------------------------

    ! Dimension partitioning
    ! ======================
    if (coord(1)==ndomx-1) then
       nx_o=ngx_o-(ndomx-1)*(ngx_o/ndomx)
    else
       nx_o=ngx_o/ndomx
    endif
    if (coord(2)==ndomy-1) then
       ny_o=ngy_o-(ndomy-1)*(ngy_o/ndomy)
    else
       ny_o=ngy_o/ndomy
    endif
    if (coord(3)==ndomz-1) then
       nz_o=ngz_o-(ndomz-1)*(ngz_o/ndomz)
    else
       nz_o=ngz_o/ndomz
    endif

!!$    ! Regular dimension (for all procs except, potentially, last in a given direction)
!!$    ! =================
!!$    if (iproc==0) then
!!$       nx_or=nx_o
!!$       ny_or=ny_o
!!$       nz_or=nz_o
!!$    endif
!!$    ! Share regular dims (send to procs 1~>proc-1)
!!$    call MPI_BCAST(nx_or,1,MPI_INTEGER,0,MPI_COMM_WORLD,info)
!!$    call MPI_BCAST(ny_or,1,MPI_INTEGER,0,MPI_COMM_WORLD,info)
!!$    call MPI_BCAST(nz_or,1,MPI_INTEGER,0,MPI_COMM_WORLD,info)

    ! Regular dimension (for all procs except, potentially, last in a given direction)
    ! =================
    if (iproc==iproc_leader(nob(iproc))) then
       nx_or=nx_o
       ny_or=ny_o
       nz_or=nz_o
    endif
    ! Share regular dims (send to procs 1~>proc-1)
    call MPI_BCAST(nx_or,1,MPI_INTEGER,0,COMM_intrablock,info)
    call MPI_BCAST(ny_or,1,MPI_INTEGER,0,COMM_intrablock,info)
    call MPI_BCAST(nz_or,1,MPI_INTEGER,0,COMM_intrablock,info)

!!$    print *,iproc,nx_o,ny_o,nx_or,ny_or
!!$    call mpistop('h',0)
    
    ! Partitioning of the grid on MPI proc
    ! ====================================
    if (is_curv) then
       ! Partitioning of the curvilinear grid
       ! ------------------------------------
       allocate(xc_o(1-ngh:nx_o+ngh,1-ngh:ny_o+ngh))
       allocate(yc_o(1-ngh:nx_o+ngh,1-ngh:ny_o+ngh))
       do i=1-ngh,nx_o+ngh
          do j=1-ngh,ny_o+ngh
             xc_o(i,j)=xgc_o(i+coord(1)*nx_or,j+coord(2)*ny_or)
             yc_o(i,j)=ygc_o(i+coord(1)*nx_or,j+coord(2)*ny_or)
          enddo
       enddo
    else
       ! Partitioning of the Cartesian grid
       ! ----------------------------------
       allocate(x_o(1-ngh:nx_o+ngh),y_o(1-ngh:ny_o+ngh))  
       do i=1-ngh,nx_o+ngh
          x_o(i)=xg_o(i+coord(1)*nx_o)
       enddo

       do j=1-ngh,ny_o+ngh
          y_o(j)=yg_o(j+coord(2)*ny_o)
       enddo
    end if

    ! Only Cartesian third direction is accepted
    ! ------------------------------------------
    allocate(z_o(1-ngh:nz_o+ngh))
    if (is_2D) then
       z_o(1)=zg_o(1) ! useful ???????????????????????????
    else
       do k=1-ngh,nz_o+ngh
          z_o(k)=zg_o(k+coord(3)*nz_o)
       enddo
    endif

  end subroutine mpi_part_o

  !===============================================================
  subroutine mpi_types_comm_o
  !===============================================================
    !> Definition of MPI types for communications
  !===============================================================
    use mod_grid_directions
    use mod_mpi_part
    implicit none
    ! ------------------------------------------------------------
    integer :: nex,ney ! extended sizes
    integer :: sign_i,sign_j
    integer :: memjump ! jump in memory between block starts
    ! strides in memory between block starts
    integer(kind=MPI_ADDRESS_KIND) :: stride,stride2
    integer :: type_base,type_edge
    !integer :: sizeofreal ! size of MPI real numbers
    ! ------------------------------------------------------------

    ! Extended sizes (+ ghost cells)
    ! ------------------------------
    nex=nx_o+2*ngh
    ney=ny_o+2*ngh

    ! Determine 'sizeofreal'
    ! ---------------------
    !call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION,sizeofreal,info)
    
    ! stride between edges along x or y
    ! ---------------------------------
    stride2=nex*ney*sizeofreal

    ! MPI-type construction for imin face "along y" (W: west)
    ! -------------------------------------------------------
    ! for imin/W, parallel dir is j and normal dir is i
    inW_o=1
    inWs_o=1
    inwr_o=1-ngh
    if (neighbor(nW).ne.MPI_PROC_NULL) then
      if (nob(neighbor(nW))/=nob(iproc)) inWs_o=2
    endif
    ipW_o=1  
    sign_i=1
    sign_j=1
    ! if rev parallel -> +j => -j
    if (is_rev(1,1)) then
       sign_j=-1
       ipW_o=ny_o
    endif
    ! if rev  normal  -> +i => -i
    if (is_rev(1,2)) then
       sign_i=-1
       inW_o=ngh
       inWs_o=ngh+1
       inwr_o=0
    endif
    if (is_swapij(1)) then
       ! if swap invert role of i and j
       memjump=sign_j*nex  ! j
       stride =sign_i*sizeofreal ! i
       call MPI_TYPE_VECTOR(ny_o,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride,type_base,type_edge,info)
    else
       ! default
       memjump=sign_i ! i
       stride =sign_j*nex*sizeofreal ! j
       call MPI_TYPE_VECTOR(ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ny_o,1,stride,type_base,type_edge,info)
    endif
    call MPI_TYPE_CREATE_HVECTOR(nz_o,1,stride2,type_edge,type_faceW_o,info)
    call MPI_TYPE_COMMIT(type_faceW_o,info)

    ! MPI-type construction for imax face "along y" (E: east)
    ! -------------------------------------------------------
    ! for imax/E, parallel dir is j and normal dir is i
    inE_o=nx_o+1
    inEs_o=nx_o+1-ngh
    inEr_o=nx_o+1
    if (neighbor(nE).ne.MPI_PROC_NULL) then
      if (nob(neighbor(nE))/=nob(iproc)) inEs_o=nx_o-ngh
    endif
    ipE_o=1
    sign_i=1
    sign_j=1
    ! if rev parallel -> +j => -j
    if (is_rev(2,1)) then
       sign_j=-1
       ipE_o=ny_o
    endif
    ! if rev  normal  -> +i => -i
    if (is_rev(2,2)) then
       sign_i=-1
       inE_o=nx_o+ngh
       inEs_o=nx_o-1
       inEr_o=nx_o+ngh
    endif
    if (is_swapij(2)) then
       ! if swap invert role of i and j
       memjump=sign_j*nex  ! j
       stride =sign_i*sizeofreal ! i
       call MPI_TYPE_VECTOR(ny_o,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride,type_base,type_edge,info)
    else
       ! default
       memjump=sign_i ! i
       stride =sign_j*nex*sizeofreal ! j
       call MPI_TYPE_VECTOR(ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ny_o,1,stride,type_base,type_edge,info)
    endif
    call MPI_TYPE_CREATE_HVECTOR(nz_o,1,stride2,type_edge,type_faceE_o,info)
    call MPI_TYPE_COMMIT(type_faceE_o,info)

    ! MPI-type construction for jmin face "along x" (S: south)
    ! --------------------------------------------------------
    ! for jmin/S, parallel dir is i and normal dir is j
    ipS_o=1
    inS_o=1
    inSs_o=1
    inSr_o=1-ngh
    if (neighbor(nS).ne.MPI_PROC_NULL) then
      if (nob(neighbor(nS))/=nob(iproc)) inSs_o=2
    endif
    sign_i=1
    sign_j=1
    ! if rev parallel -> +i => -i
    if (is_rev(3,1)) then
       sign_i=-1
       ipS_o=nx_o
    endif
    ! if rev  normal  -> +j => -j
    if (is_rev(3,2)) then
       sign_j=-1
       inS_o=ngh
       inSs_o=ngh+1
       inSr_o=0
    endif
    if (is_swapij(3)) then
       ! if swap invert role of i and j
       memjump=sign_j*nex  ! j
       stride =sign_i*sizeofreal ! i
       call MPI_TYPE_VECTOR(ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(nx_o,1,stride,type_base,type_edge,info)
    else
       ! default
       memjump=sign_i ! i
       stride =sign_j*nex*sizeofreal ! j
       call MPI_TYPE_VECTOR(nx_o,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride,type_base,type_edge,info)
    endif
    call MPI_TYPE_CREATE_HVECTOR(nz_o,1,stride2,type_edge,type_faceS_o,info)
    call MPI_TYPE_COMMIT(type_faceS_o,info)

    ! MPI-type construction for jmax face "along x" (N: north)
    ! --------------------------------------------------------
    ! for jmax/N, parallel dir is i and normal dir is j
    ipN_o=1
    inN_o=ny_o+1
    inNs_o=ny_o+1-ngh
    inNr_o=ny_o+1
    if (neighbor(nN).ne.MPI_PROC_NULL) then
      if (nob(neighbor(nN))/=nob(iproc)) inNs_o=ny_o-ngh
    endif
    sign_i=1
    sign_j=1
    ! if rev parallel -> +i => -i
    if (is_rev(4,1)) then
       sign_i=-1
       ipN_o=nx_o
    endif
    ! if rev  normal  -> +j => -j
    if (is_rev(4,2)) then
       sign_j=-1
       inN_o=ny_o+ngh
       inNs_o=ny_o-1
       inNr_o=ny_o+ngh
    endif
    if (is_swapij(4)) then
       ! if swap invert role of i and j
       memjump=sign_j*nex  ! j
       stride =sign_i*sizeofreal ! i
       call MPI_TYPE_VECTOR(ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(nx_o,1,stride,type_base,type_edge,info)
    else
       ! default
       memjump=sign_i ! i
       stride =sign_j*nex*sizeofreal ! j
       call MPI_TYPE_VECTOR(nx_o,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride,type_base,type_edge,info)
    endif
    call MPI_TYPE_CREATE_HVECTOR(nz_o,1,stride2,type_edge,type_faceN_o,info)
    call MPI_TYPE_COMMIT(type_faceN_o,info)

    stride=nex*sizeofreal
    call MPI_TYPE_VECTOR(nx_o,1,1,MPI_DOUBLE_PRECISION,type_base,info)

    ! MPI-type construction for kmin face "along z" (F: forward)
    ! --------------------------------------------------------
    call MPI_TYPE_CREATE_HVECTOR( ny_o,1, stride,type_base,type_edge,info)
    call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride2,type_edge,type_faceF_o,info)
    call MPI_TYPE_COMMIT(type_faceF_o,info)
    
    ! MPI-type construction for kmax face "along z" (B: backward)
    ! -----------------------------------------------------------
    call MPI_TYPE_CREATE_HVECTOR( ny_o,1, stride,type_base,type_edge,info)
    call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride2,type_edge,type_faceB_o,info)
    call MPI_TYPE_COMMIT(type_faceB_o,info)

  end subroutine mpi_types_comm_o

  !===============================================================
  subroutine mpi_types_comm_edges_o
  !===============================================================
    !> Definition of MPI types for communications
  !===============================================================
    use mod_grid_directions
    use mod_mpi_part
    implicit none
    ! ------------------------------------------------------------
    integer :: nex,ney ! extended sizes
    integer :: sign_i,sign_j
    integer :: memjump ! jump in memory between block starts
    ! strides in memory between block starts
    integer(kind=MPI_ADDRESS_KIND) :: stride,stride2
    integer :: type_base,type_edge
    ! ------------------------------------------------------------

    ! MPI types for faces+edges (ngh ghost cells)
    ! =========================

    ! Extended sizes (+ ghost cells)
    ! ------------------------------
    nex=nx_o+2*ngh
    ney=ny_o+2*ngh

    ! stride between edges along x or y
    ! ---------------------------------
    stride2=nex*ney*sizeofreal

    !!!
    ! REMARK ! changes in ipW_e but no changes in inW_e
    !!!
    
    ! MPI-type construction for imin face "along y" (W: west)
    ! -------------------------------------------------------
    ! for imin/W, parallel dir is j and normal dir is i
    inW_e=1
    ipW_e=1-ngh
    sign_i=1
    sign_j=1
    ! if rev parallel -> +j => -j
    if (is_rev(1,1)) then
       sign_j=-1
       ipW_e=ny_o+ngh
    endif
    ! if rev  normal  -> +i => -i
    if (is_rev(1,2)) then
       sign_i=-1
       inW_e=ngh
    endif
    if (is_swapij(1)) then
       ! if swap invert role of i and j
       memjump=sign_j*nex  ! j
       stride =sign_i*sizeofreal ! i
       call MPI_TYPE_VECTOR(ny_o+2*ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride,type_base,type_edge,info)
    else
       ! default
       memjump=sign_i ! i
       stride =sign_j*nex*sizeofreal ! j
       call MPI_TYPE_VECTOR(ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ny_o+2*ngh,1,stride,type_base,type_edge,info)
    endif
    call MPI_TYPE_CREATE_HVECTOR(nz_o,1,stride2,type_edge,type_e_faceW,info)
    call MPI_TYPE_COMMIT(type_e_faceW,info)

    ! MPI-type construction for imax face "along y" (E: east)
    ! -------------------------------------------------------
    ! for imax/E, parallel dir is j and normal dir is i
    inE_e=nx_o+1
    ipE_e=1-ngh
    sign_i=1
    sign_j=1
    ! if rev parallel -> +j => -j
    if (is_rev(2,1)) then
       sign_j=-1
       ipE_e=ny_o+ngh
    endif
    ! if rev  normal  -> +i => -i
    if (is_rev(2,2)) then
       sign_i=-1
       inE_e=nx_o+ngh
    endif
    if (is_swapij(2)) then
       ! if swap invert role of i and j
       memjump=sign_j*nex  ! j
       stride =sign_i*sizeofreal ! i
       call MPI_TYPE_VECTOR(ny_o+2*ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride,type_base,type_edge,info)
    else
       ! default
       memjump=sign_i ! i
       stride =sign_j*nex*sizeofreal ! j
       call MPI_TYPE_VECTOR(ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ny_o+2*ngh,1,stride,type_base,type_edge,info)
    endif
    call MPI_TYPE_CREATE_HVECTOR(nz_o,1,stride2,type_edge,type_e_faceE,info)
    call MPI_TYPE_COMMIT(type_e_faceE,info)

    ! MPI-type construction for jmin face "along x" (S: south)
    ! --------------------------------------------------------
    ! for jmin/S, parallel dir is i and normal dir is j
    ipS_e=1-ngh
    inS_e=1
    sign_i=1
    sign_j=1
    ! if rev parallel -> +i => -i
    if (is_rev(3,1)) then
       sign_i=-1
       ipS_e=nx_o+ngh
    endif
    ! if rev  normal  -> +j => -j
    if (is_rev(3,2)) then
       sign_j=-1
       inS_e=ngh
    endif
    if (is_swapij(3)) then
       ! if swap invert role of i and j
       memjump=sign_j*nex  ! j
       stride =sign_i*sizeofreal ! i
       call MPI_TYPE_VECTOR(ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(nx_o+2*ngh,1,stride,type_base,type_edge,info)
    else
       ! default
       memjump=sign_i ! i
       stride =sign_j*nex*sizeofreal ! j
       call MPI_TYPE_VECTOR(nx_o+2*ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride,type_base,type_edge,info)
    endif
    call MPI_TYPE_CREATE_HVECTOR(nz_o,1,stride2,type_edge,type_e_faceS,info)
    call MPI_TYPE_COMMIT(type_e_faceS,info)

    ! MPI-type construction for jmax face "along x" (N: north)
    ! --------------------------------------------------------
    ! for jmax/N, parallel dir is i and normal dir is j
    ipN_e=1-ngh
    inN_e=ny_o+1
    sign_i=1
    sign_j=1
    ! if rev parallel -> +i => -i
    if (is_rev(4,1)) then
       sign_i=-1
       ipN_e=nx_o+ngh
    endif
    ! if rev  normal  -> +j => -j
    if (is_rev(4,2)) then
       sign_j=-1
       inN_e=ny_o+ngh
    endif
    if (is_swapij(4)) then
       ! if swap invert role of i and j
       memjump=sign_j*nex  ! j
       stride =sign_i*sizeofreal ! i
       call MPI_TYPE_VECTOR(ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(nx_o+2*ngh,1,stride,type_base,type_edge,info)
    else
       ! default
       memjump=sign_i ! i
       stride =sign_j*nex*sizeofreal ! j
       call MPI_TYPE_VECTOR(nx_o+2*ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride,type_base,type_edge,info)
    endif
    call MPI_TYPE_CREATE_HVECTOR(nz_o,1,stride2,type_edge,type_e_faceN,info)
    call MPI_TYPE_COMMIT(type_e_faceN,info)

    stride=nex*sizeofreal
    call MPI_TYPE_VECTOR(nx_o+2*ngh,1,1,MPI_DOUBLE_PRECISION,type_base,info)
    
    ! MPI-type construction for kmin face "along z" (F: forward)
    ! --------------------------------------------------------    
    call MPI_TYPE_CREATE_HVECTOR(ny_o+2*ngh,1,stride,type_base,type_edge,info)
    call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride2,type_edge,type_e_faceF,info)
    call MPI_TYPE_COMMIT(type_e_faceF,info)
    
    ! MPI-type construction for kmax face "along z" (B: backward)
    ! -----------------------------------------------------------
    call MPI_TYPE_CREATE_HVECTOR(ny_o+2*ngh,1,stride,type_base,type_edge,info)
    call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride2,type_edge,type_e_faceB,info)
    call MPI_TYPE_COMMIT(type_e_faceB,info)

  end subroutine mpi_types_comm_edges_o

  !===============================================================
  subroutine mpi_types_comm_o_old
  !===============================================================
    !> Define MPI types for communications (version old field)
  !===============================================================
    use mod_mpi
    implicit none
    ! ------------------------------------------------------------
    integer(kind=MPI_ADDRESS_KIND) :: pas1,pas2
    integer :: sizeofreal
    integer :: type_base1,type_base2,type_base3
    integer :: type_line1,type_line2,type_column
    ! ------------------------------------------------------------

    ! Determine 'sizeofreal'
    ! ======================
    call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION,sizeofreal,info)

    ! MPI types for faces (ngh ghost cells)
    ! ===================
    pas1 = (nx_o+2*ngh)*sizeofreal
    pas2 = (nx_o+2*ngh)*(ny_o+2*ngh)*sizeofreal
    
    ! Construction of base types
    call MPI_TYPE_VECTOR(nx_o,1,1,MPI_DOUBLE_PRECISION,type_base1,info)
    call MPI_TYPE_VECTOR(ngh,1,1,MPI_DOUBLE_PRECISION,type_base2,info)

    ! Construction of type "face along x" (north or south)
    call MPI_TYPE_CREATE_HVECTOR(ngh,1,pas1,type_base1,type_line1,info)
    call MPI_TYPE_CREATE_HVECTOR( nz_o,1,pas2,type_line1,type_facex_o,info)
    call MPI_TYPE_COMMIT(type_facex_o,info)

    ! Construction of type "face along y" (west or east)
    call MPI_TYPE_CREATE_HVECTOR(ny_o,1,pas1,type_base2,type_column,info)
    call MPI_TYPE_CREATE_HVECTOR(nz_o,1,pas2,type_column,type_facey_o,info)
    call MPI_TYPE_COMMIT(type_facey_o,info)

    ! Construction of type "face along z" (front or back)
    call MPI_TYPE_CREATE_HVECTOR( ny_o,1,pas1,type_base1,type_line2,info)
    call MPI_TYPE_CREATE_HVECTOR(ngh,1,pas2,type_line2,type_facez_o,info)
    call MPI_TYPE_COMMIT(type_facez_o,info)
  
    ! MPI types for z-edges 
    ! =====================
    call MPI_TYPE_CREATE_HVECTOR( ngh,1,pas1,type_base2,type_base3,info)
    call MPI_TYPE_CREATE_HVECTOR(nz_o,1,pas2,type_base3,type_linez_o,info)
    call MPI_TYPE_COMMIT(type_linez_o,info)

  end subroutine mpi_types_comm_o_old

  !===============================================================================
  subroutine communic2d_o(var)
  !===============================================================================
    !> NSEW communications using blocking SENDRECV (old field 2D)
  !===============================================================================
    use mod_mpi_part
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), dimension(1-ngh:nx_o+ngh,1-ngh:ny_o+ngh,1-ngh:nz_o+ngh), intent(inout) :: var
    ! ----------------------------------------------------------------------------

    ! Send to neighbor N and reception from neighbor S
    call MPI_SENDRECV(var(1,ny_o-ngh+1,1),1,type_facex_o,neighbor(nN),tag, &
                      var(1,    -ngh+1,1),1,type_facex_o,neighbor(nS),tag,COMM_global,status,info)

    ! Send to neighbor S and reception from neighbor N
    call MPI_SENDRECV(var(1,     1,1),1,type_facex_o,neighbor(nS),tag, &
                      var(1,ny_o+1,1),1,type_facex_o,neighbor(nN),tag,COMM_global,status,info)

    ! Send to neighbor E and reception from neighbor W
    call MPI_SENDRECV(var(nx_o-ngh+1,1,1),1,type_facey_o,neighbor(nE),tag, &
                      var(    -ngh+1,1,1),1,type_facey_o,neighbor(nW),tag,COMM_global,status,info)

    ! Send to neighbor W and reception from neighbor E
    call MPI_SENDRECV(var(     1,1,1),1,type_facey_o,neighbor(nW),tag, &
                      var(nx_o+1,1,1),1,type_facey_o,neighbor(nE),tag,COMM_global,status,info)

  end subroutine communic2d_o

  !===============================================================================
  subroutine communic3d_o(var)
  !===============================================================================
    !> NSEW communications using blocking SENDRECV (old field 3D)
  !===============================================================================
    use mod_mpi_part
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), dimension(1-ngh:nx_o+ngh,1-ngh:ny_o+ngh,1-ngh:nz_o+ngh), intent(inout) :: var
    ! ----------------------------------------------------------------------------

    ! Send to neighbor N and reception from neighbor S
    call MPI_SENDRECV(var(1,ny_o-ngh+1,1),1,type_facex_o,neighbor(nN),tag, &
                      var(1,    -ngh+1,1),1,type_facex_o,neighbor(nS),tag,COMM_global,status,info)

    ! Send to neighbor S and reception from neighbor N
    call MPI_SENDRECV(var(1,     1,1),1,type_facex_o,neighbor(nS),tag, &
                      var(1,ny_o+1,1),1,type_facex_o,neighbor(nN),tag,COMM_global,status,info)

    ! Send to neighbor E and reception from neighbor W
    call MPI_SENDRECV(var(nx_o-ngh+1,1,1),1,type_facey_o,neighbor(nE),tag, &
                      var(    -ngh+1,1,1),1,type_facey_o,neighbor(nW),tag,COMM_global,status,info)

    ! Send to neighbor W and reception from neighbor E
    call MPI_SENDRECV(var(     1,1,1),1,type_facey_o,neighbor(nW),tag, &
                      var(nx_o+1,1,1),1,type_facey_o,neighbor(nE),tag,COMM_global,status,info)

    ! Send to neighbor B and reception from neighbor F
    call MPI_SENDRECV(var(1,1,nz_o-ngh+1),1,type_facez_o,neighbor(nB),tag, &
                      var(1,1,    -ngh+1),1,type_facez_o,neighbor(nF),tag,COMM_global,status,info)

    ! Send to neighbor F and reception from neighbor B
    call MPI_SENDRECV(var(1,1,     1),1,type_facez_o,neighbor(nF),tag, &
                      var(1,1,nz_o+1),1,type_facez_o,neighbor(nB),tag,COMM_global,status,info)

  end subroutine communic3d_o

  !===============================================================================
  subroutine commun2d_o(var,i)
  !===============================================================================
    !> NSEW communications using non-blocking ISEND/IRECV (2D faces only)
  !===============================================================================
    use mod_mpi_part
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,k
    real(wp), dimension(1-ngh:nx_o+ngh,1-ngh:ny_o+ngh,1-ngh:nz_o+ngh), intent(inout) :: var
    ! ----------------------------------------------------------------------------

    ! Counter: 2(send/recv)*4(neighbors)=8*(variable number)
    k=8*(i-1)
    
    if (is_adjoint_block) then ! for adjoint block interfaces
    
      ! Send to neighbor W and reception from neighbor W
      call MPI_ISEND(var(inWs_o,ipW_o,1),1,type_faceW_o,neighbor(nW),tags(nW),COMM_global,request_2d(k+1),info)
      call MPI_IRECV(var(inWr_o,ipW_o,1),1,type_faceW_o,neighbor(nW),tagr(nW),COMM_global,request_2d(k+2),info)
      
      ! Send to neighbor E and reception from neighbor E
      call MPI_ISEND(var(inEs_o,ipE_o,1),1,type_faceE_o,neighbor(nE),tags(nE),COMM_global,request_2d(k+3),info)
      call MPI_IRECV(var(inEr_o,ipE_o,1),1,type_faceE_o,neighbor(nE),tagr(nE),COMM_global,request_2d(k+4),info)
      
      ! Send to neighbor S and reception from neighbor S
      call MPI_ISEND(var(ipS_o,inSs_o,1),1,type_faceS_o,neighbor(nS),tags(nS),COMM_global,request_2d(k+5),info)
      call MPI_IRECV(var(ipS_o,inSr_o,1),1,type_faceS_o,neighbor(nS),tagr(nS),COMM_global,request_2d(k+6),info)
      
      ! Send to neighbor N and reception from neighbor N
      call MPI_ISEND(var(ipN_o,inNs_o,1),1,type_faceN_o,neighbor(nN),tags(nN),COMM_global,request_2d(k+7),info)
      call MPI_IRECV(var(ipN_o,inNr_o,1),1,type_faceN_o,neighbor(nN),tagr(nN),COMM_global,request_2d(k+8),info)
    
    else

      ! Send to neighbor W and reception from neighbor W
      call MPI_ISEND(var(inW_o    ,ipW_o,1),1,type_faceW_o,neighbor(nW),tags(nW),COMM_global,request_2d(k+1),info)
      call MPI_IRECV(var(inW_o-ngh,ipW_o,1),1,type_faceW_o,neighbor(nW),tagr(nW),COMM_global,request_2d(k+2),info)
    
      ! Send to neighbor E and reception from neighbor E
      call MPI_ISEND(var(inE_o-ngh,ipE_o,1),1,type_faceE_o,neighbor(nE),tags(nE),COMM_global,request_2d(k+3),info)
      call MPI_IRECV(var(inE_o    ,ipE_o,1),1,type_faceE_o,neighbor(nE),tagr(nE),COMM_global,request_2d(k+4),info)
    
      ! Send to neighbor S and reception from neighbor S
      call MPI_ISEND(var(ipS_o,inS_o    ,1),1,type_faceS_o,neighbor(nS),tags(nS),COMM_global,request_2d(k+5),info)
      call MPI_IRECV(var(ipS_o,inS_o-ngh,1),1,type_faceS_o,neighbor(nS),tagr(nS),COMM_global,request_2d(k+6),info)
    
      ! Send to neighbor N and reception from neighbor N
      call MPI_ISEND(var(ipN_o,inN_o-ngh,1),1,type_faceN_o,neighbor(nN),tags(nN),COMM_global,request_2d(k+7),info)
      call MPI_IRECV(var(ipN_o,inN_o    ,1),1,type_faceN_o,neighbor(nN),tagr(nN),COMM_global,request_2d(k+8),info)
      
    endif  

  end subroutine commun2d_o
  
  !===============================================================================
  subroutine commun3d_o(var,i)
  !===============================================================================
    !> NSEWFB communications using non-blocking ISEND/IRECV (faces only)
  !===============================================================================
    use mod_mpi_part
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,k
    real(wp), dimension(1-ngh:nx_o+ngh,1-ngh:ny_o+ngh,1-ngh:nz_o+ngh), intent(inout) :: var
    ! ----------------------------------------------------------------------------

    ! Counter: 2(send/recv)*6(neighbors)=12*(variable number)
    k=12*(i-1)
    
    if (is_adjoint_block) then ! for adjoint block interfaces

      ! Send to neighbor W and reception from neighbor W
      call MPI_ISEND(var(inWs_o,ipW_o,1),1,type_faceW_o,neighbor(nW),tags(nW),COMM_global,request_3d(k+1),info)
      call MPI_IRECV(var(inWr_o,ipW_o,1),1,type_faceW_o,neighbor(nW),tagr(nW),COMM_global,request_3d(k+2),info)
      
      ! Send to neighbor E and reception from neighbor E
      call MPI_ISEND(var(inEs_o,ipE_o,1),1,type_faceE_o,neighbor(nE),tags(nE),COMM_global,request_3d(k+3),info)
      call MPI_IRECV(var(inEr_o,ipE_o,1),1,type_faceE_o,neighbor(nE),tagr(nE),COMM_global,request_3d(k+4),info)
      
      ! Send to neighbor S and reception from neighbor S
      call MPI_ISEND(var(ipS_o,inSs_o,1),1,type_faceS_o,neighbor(nS),tags(nS),COMM_global,request_3d(k+5),info)
      call MPI_IRECV(var(ipS_o,inSr_o,1),1,type_faceS_o,neighbor(nS),tagr(nS),COMM_global,request_3d(k+6),info)
      
      ! Send to neighbor N and reception from neighbor N
      call MPI_ISEND(var(ipN_o,inNs_o,1),1,type_faceN_o,neighbor(nN),tags(nN),COMM_global,request_3d(k+7),info)
      call MPI_IRECV(var(ipN_o,inNr_o,1),1,type_faceN_o,neighbor(nN),tagr(nN),COMM_global,request_3d(k+8),info)
    
    else
    
      ! Send to neighbor W and reception from neighbor W
      call MPI_ISEND(var(inW_o    ,ipW_o,1),1,type_faceW_o,neighbor(nW),tags(nW),COMM_global,request_3d(k+1),info)
      call MPI_IRECV(var(inW_o-ngh,ipW_o,1),1,type_faceW_o,neighbor(nW),tagr(nW),COMM_global,request_3d(k+2),info)
    
      ! Send to neighbor E and reception from neighbor E
      call MPI_ISEND(var(inE_o-ngh,ipE_o,1),1,type_faceE_o,neighbor(nE),tags(nE),COMM_global,request_3d(k+3),info)
      call MPI_IRECV(var(inE_o    ,ipE_o,1),1,type_faceE_o,neighbor(nE),tagr(nE),COMM_global,request_3d(k+4),info)
    
      ! Send to neighbor S and reception from neighbor S
      call MPI_ISEND(var(ipS_o,inS_o    ,1),1,type_faceS_o,neighbor(nS),tags(nS),COMM_global,request_3d(k+5),info)
      call MPI_IRECV(var(ipS_o,inS_o-ngh,1),1,type_faceS_o,neighbor(nS),tagr(nS),COMM_global,request_3d(k+6),info)
    
      ! Send to neighbor N and reception from neighbor N
      call MPI_ISEND(var(ipN_o,inN_o-ngh,1),1,type_faceN_o,neighbor(nN),tags(nN),COMM_global,request_3d(k+7),info)
      call MPI_IRECV(var(ipN_o,inN_o    ,1),1,type_faceN_o,neighbor(nN),tagr(nN),COMM_global,request_3d(k+8),info)
    
    endif

    ! Send to neighbor F and reception from neighbor F
    call MPI_ISEND(var(1,1,    1),1,type_faceF_o,neighbor(nF),tags(nF),COMM_global,request_3d(k+9),info)
    call MPI_IRECV(var(1,1,1-ngh),1,type_faceF_o,neighbor(nF),tagr(nF),COMM_global,request_3d(k+10),info)

    ! Send to neighbor B and reception from neighbor B
    call MPI_ISEND(var(1,1,nz_o-ngh+1),1,type_faceB_o,neighbor(nB),tags(nB),COMM_global,request_3d(k+11),info)
    call MPI_IRECV(var(1,1,    nz_o+1),1,type_faceB_o,neighbor(nB),tagr(nB),COMM_global,request_3d(k+12),info)

  end subroutine commun3d_o
  
  !===============================================================================
  subroutine commun3d_edges(var,i)
  !===============================================================================
    !> NSEWFB communications using non-blocking ISEND/IRECV (faces only)
  !===============================================================================
    use mod_mpi_part
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,k
    real(wp), dimension(1-ngh:nx_o+ngh,1-ngh:ny_o+ngh,1-ngh:nz_o+ngh), intent(inout) :: var
    ! ----------------------------------------------------------------------------

    ! Counter: 2(send/recv)*6(neighbors)=12*(variable number)
    k=12*(i-1)

    ! Send to neighbor W and reception from neighbor W
    call MPI_ISEND(var(inW_e    ,ipW_e,1),1,type_e_faceW,neighbor(nW),tags(nW),COMM_global,request_3d(k+1),info)
    call MPI_IRECV(var(inW_e-ngh,ipW_e,1),1,type_e_faceW,neighbor(nW),tagr(nW),COMM_global,request_3d(k+2),info)

    ! Send to neighbor E and reception from neighbor E
    call MPI_ISEND(var(inE_e-ngh,ipE_e,1),1,type_e_faceE,neighbor(nE),tags(nE),COMM_global,request_3d(k+3),info)
    call MPI_IRECV(var(inE_e    ,ipE_e,1),1,type_e_faceE,neighbor(nE),tagr(nE),COMM_global,request_3d(k+4),info)

    ! Send to neighbor S and reception from neighbor S
    call MPI_ISEND(var(ipS_e,inS_e    ,1),1,type_e_faceS,neighbor(nS),tags(nS),COMM_global,request_3d(k+5),info)
    call MPI_IRECV(var(ipS_e,inS_e-ngh,1),1,type_e_faceS,neighbor(nS),tagr(nS),COMM_global,request_3d(k+6),info)

    ! Send to neighbor N and reception from neighbor N
    call MPI_ISEND(var(ipN_e,inN_e-ngh,1),1,type_e_faceN,neighbor(nN),tags(nN),COMM_global,request_3d(k+7),info)
    call MPI_IRECV(var(ipN_e,inN_e    ,1),1,type_e_faceN,neighbor(nN),tagr(nN),COMM_global,request_3d(k+8),info)

    ! Send to neighbor F and reception from neighbor F
    call MPI_ISEND(var(1-ngh,1-ngh,    1),1,type_e_faceF,neighbor(nF),tags(nF),COMM_global,request_3d(k+ 9),info)
    call MPI_IRECV(var(1-ngh,1-ngh,1-ngh),1,type_e_faceF,neighbor(nF),tagr(nF),COMM_global,request_3d(k+10),info)

    ! Send to neighbor B and reception from neighbor B
    call MPI_ISEND(var(1-ngh,1-ngh,nz_o-ngh+1),1,type_e_faceB,neighbor(nB),tags(nB),COMM_global,request_3d(k+11),info)
    call MPI_IRECV(var(1-ngh,1-ngh,    nz_o+1),1,type_e_faceB,neighbor(nB),tagr(nB),COMM_global,request_3d(k+12),info)

  end subroutine commun3d_edges

  !===============================================================================
  subroutine commun2d_edges(var,i)
  !===============================================================================
    !> NSEW communications using non-blocking ISEND/IRECV (2D faces only)
  !===============================================================================
    use mod_mpi_part
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,k
    real(wp), dimension(1-ngh:nx_o+ngh,1-ngh:ny_o+ngh,1-ngh:nz_o+ngh), intent(inout) :: var
    ! ----------------------------------------------------------------------------

    ! Counter: 2(send/recv)*4(neighbors)=8*(variable number)
    k=8*(i-1)

    ! Send to neighbor W and reception from neighbor W
    call MPI_ISEND(var(inW_e    ,ipW_e,1),1,type_e_faceW,neighbor(nW),tags(nW),COMM_global,request_2d(k+1),info)
    call MPI_IRECV(var(inW_e-ngh,ipW_e,1),1,type_e_faceW,neighbor(nW),tagr(nW),COMM_global,request_2d(k+2),info)

    ! Send to neighbor E and reception from neighbor E
    call MPI_ISEND(var(inE_e-ngh,ipE_e,1),1,type_e_faceE,neighbor(nE),tags(nE),COMM_global,request_2d(k+3),info)
    call MPI_IRECV(var(inE_e    ,ipE_e,1),1,type_e_faceE,neighbor(nE),tagr(nE),COMM_global,request_2d(k+4),info)

    ! Send to neighbor S and reception from neighbor S
    call MPI_ISEND(var(ipS_e,inS_e    ,1),1,type_e_faceS,neighbor(nS),tags(nS),COMM_global,request_2d(k+5),info)
    call MPI_IRECV(var(ipS_e,inS_e-ngh,1),1,type_e_faceS,neighbor(nS),tagr(nS),COMM_global,request_2d(k+6),info)

    ! Send to neighbor N and reception from neighbor N
    call MPI_ISEND(var(ipN_e,inN_e-ngh,1),1,type_e_faceN,neighbor(nN),tags(nN),COMM_global,request_2d(k+7),info)
    call MPI_IRECV(var(ipN_e,inN_e    ,1),1,type_e_faceN,neighbor(nN),tagr(nN),COMM_global,request_2d(k+8),info)

  end subroutine commun2d_edges
  
  !===============================================================================
  subroutine communic_edgesz_o(var)
  !===============================================================================
    !> Communications of edgesz using blocking SENDRECV (old field)
  !===============================================================================
    use mod_mpi_part
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), dimension(1-ngh:nx_o+ngh,1-ngh:ny_o+ngh,1-ngh:nz_o+ngh), intent(inout) :: var
    ! ----------------------------------------------------------------------------

    ! Send to neighbor NE and reception from neighbor SW
    call MPI_SENDRECV(var(nx_o-ngh+1,ny_o-ngh+1,1),1,type_linez_o,neighbor2(nNE),tag, &
                      var(    -ngh+1,    -ngh+1,1),1,type_linez_o,neighbor2(nSW),tag,COMM_global,status,info)

    ! Send to neighbor SW and reception from neighbor NE
    call MPI_SENDRECV(var(     1,     1,1),1,type_linez_o,neighbor2(nSW),tag, &
                      var(nx_o+1,ny_o+1,1),1,type_linez_o,neighbor2(nNE),tag,COMM_global,status,info)

    ! Send to neighbor NW and reception from neighbor SE
    call MPI_SENDRECV(var(     1,ny_o-ngh+1,1),1,type_linez_o,neighbor2(nNW),tag, &
                      var(nx_o+1,    -ngh+1,1),1,type_linez_o,neighbor2(nSE),tag,COMM_global,status,info)

    ! Send to neighbor SE and reception from neighbor NW
    call MPI_SENDRECV(var(nx_o-ngh+1,     1,1),1,type_linez_o,neighbor2(nSE),tag, &
                      var(    -ngh+1,ny_o+1,1),1,type_linez_o,neighbor2(nNW),tag,COMM_global,status,info)

  end subroutine communic_edgesz_o

  !===============================================================
  subroutine free_field_o
  !===============================================================
    !> Free variables (allocated & MPI types) for old field
  !===============================================================
    use mod_mpi
    use mod_grid
    use mod_io
    implicit none
    ! ------------------------------------------------------------
    ! ------------------------------------------------------------
    
    ! Deallocate variable for old field
    ! =================================
    if (is_curv) then
       deallocate(xgc_o,ygc_o,xc_o,yc_o)
    else
       deallocate(xg_o,yg_o,x_o,y_o)
    endif
    deallocate(zg_o,z_o)
    deallocate(rho_o,rhou_o,rhov_o,rhow_o,rhoe_o)
   
    ! Deallocate I/O structure for old field (TO BE CHANGED - added sept2021)
    ! ======================================
    deallocate(field_o)
    
    ! Free MPI types for communications of old field
    ! ==============================================
    call MPI_TYPE_FREE(type_facex_o,info)
    call MPI_TYPE_FREE(type_facey_o,info)
    call MPI_TYPE_FREE(type_facez_o,info)
    call MPI_TYPE_FREE(type_linez_o,info)
        
  end subroutine free_field_o
    
end module mod_flow_o
