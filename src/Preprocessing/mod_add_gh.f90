!=================================================================================
module mod_add_gh
!=================================================================================
  !> Module to generate extended grid for cartesian and 2D curvilinear grid
  !> Created with some adapted NASA grid utilities routines
!=================================================================================
  use mod_grid
  implicit none
  ! ----------------------------------------------------------------------------
  ! neighbor grid
  integer :: COMM_r_grid ! communicator to read neighboring grids in MPI-IO
  integer, private :: ngxn,ngyn,ngzn
  real(wp), dimension(:,:), allocatable, private :: xgcn,ygcn
  ! ----------------------------------------------------------------------------

contains

  !===============================================================================
  subroutine add_ghost_cells_grid
  !===============================================================================
    !> add ghost cells at block junctions in cartesian & 2D curvilinear grid
  !===============================================================================
    use mod_block
    use mod_mpi_part
    use mod_grid_directions
    use mod_bc_periodicity
    use mod_constant ! for: dirGRID,nameGRID
    use mod_utils
    use mod_io
    use warnstop
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k
    integer :: ni1,nj1,ni2,nj2
    integer :: ngx_,ngy_,ngz_
    integer :: n,nv,mb
    integer :: nv1,nv2
    logical :: is_swap_ij
    logical :: is_rev_n,is_rev_p1
    logical, dimension(:,:), allocatable :: is_swapn_ij
    logical, dimension(:,:), allocatable :: is_revn_n,is_revn_p1
    real(wp) :: Lzp0,z_mid,dz0
    character(len=200) :: gridname
    ! ----------------------------------------------------------------------------

    ! Print message
    ! =============
    if ((iproc==0).and.(verbose)) then
       print *,repeat('=',70)
       print *,'Add ghost-cells at block junctions'
       print *,repeat('=',70)
    endif

    ! Get swap/reverse indicators for all blocks
    ! ==========================================
    allocate(is_swapn_ij(6,nbloc),is_revn_n(6,nbloc),is_revn_p1(6,nbloc))
    do k=1,6
       call MPI_ALLGATHER(is_swapij2_bl(k),1,MPI_LOGICAL,is_swapn_ij(k,:),1,MPI_LOGICAL,COMM_interblock,info)
       call MPI_ALLGATHER(is_rev2_bl(k,2),1,MPI_LOGICAL,is_revn_n(k,:),1,MPI_LOGICAL,COMM_interblock,info)
       call MPI_ALLGATHER(is_rev2_bl(k,1),1,MPI_LOGICAL,is_revn_p1(k,:),1,MPI_LOGICAL,COMM_interblock,info)
    enddo

    ! current block number
    ! -------------------
    n=nob(iproc)

    if (iproc==iproc_leader(n)) then
       ! Fill ghost cells for block faces
       ! ================================

       ! face imin
       ! ---------
       if (bl(n)%BC(1)>0) then
          ! neighbor
          nv=bl(n)%BC(1)
          ! read grid in neighboring block nv
          ! call read_neighbor(trim(gridname)//trim(numchar(nv))//'.x')
          call read_neighbor('grid_bl'//trim(numchar(nv))//'.bin',nv)
          ! fill ghost cells in face #1
          call fill_ghost_cells(1,1,ngy,is_swapij2_bl(1),is_rev2_bl(1,1),is_rev2_bl(1,2))
          ! free temporary array
          deallocate(xgcn,ygcn)
       endif

       ! face imax
       ! ---------
       if (bl(n)%BC(2)>0) then
          ! neighbor
          nv=bl(n)%BC(2)
          ! read grid in neighboring block nv
          ! call read_neighbor(trim(gridname)//trim(numchar(nv))//'.x')
          call read_neighbor('grid_bl'//trim(numchar(nv))//'.bin',nv)
          ! fill ghost cells in face #2
          call fill_ghost_cells(2,1,ngy,is_swapij2_bl(2),is_rev2_bl(2,1),is_rev2_bl(2,2))
          ! free temporary array
          deallocate(xgcn,ygcn)
       endif

       ! face jmin
       ! ---------
       if (bl(n)%BC(3)>0) then
          ! neighbor
          nv=bl(n)%BC(3)
          ! read grid in neighboring block nv
          ! call read_neighbor(trim(gridname)//trim(numchar(nv))//'.x')
          call read_neighbor('grid_bl'//trim(numchar(nv))//'.bin',nv)
          ! fill ghost cells in face #3
          call fill_ghost_cells(3,1,ngx,is_swapij2_bl(3),is_rev2_bl(3,1),is_rev2_bl(3,2))
          ! free temporary array
          deallocate(xgcn,ygcn)
       endif

       ! face jmax
       ! ---------
       if (bl(n)%BC(4)>0) then
          ! neighbor
          nv=bl(n)%BC(4)
          ! read grid in neighboring block nv
          ! call read_neighbor(trim(gridname)//trim(numchar(nv))//'.x')
          call read_neighbor('grid_bl'//trim(numchar(nv))//'.bin',nv)
          ! fill ghost cells in face #4
          call fill_ghost_cells(4,1,ngx,is_swapij2_bl(4),is_rev2_bl(4,1),is_rev2_bl(4,2))
          ! free temporary array
          deallocate(xgcn,ygcn)
       endif

       ! Fill ghost cells for block edges
       ! ================================
100    format(1x,' neighbors are different (',i3,' and ',i3,') ~> create artificial edge')

       ! edge imin-jmin
       ! --------------
       if ((bl(n)%BC(1)>0).and.(bl(n)%BC(3)>0)) then
          ! start from jmin
          nv=bl(n)%BC(3)
          nv2=bl(nv)%BC(1)
          ! start from imin
          nv=bl(n)%BC(1)
          nv1=bl(nv)%BC(3)
          if (nv1.ne.nv2) then
             ! construct artificial edge imin-jmin
             call cons_grid(1,1)
          else
             ! read grid in neighboring block nv
             ! call read_neighbor(trim(gridname)//trim(numchar(nv1))//'.x')
          call read_neighbor('grid_bl'//trim(numchar(nv))//'.bin',nv)
             ! determine swap/reverse for corner block
             is_swap_ij=is_swapij2_bl(1).or.is_swapn_ij(3,nv)
             is_rev_n  = is_rev2_bl(1,2).or.is_revn_p1(3,nv)
             is_rev_p1 = is_rev2_bl(1,1).or.is_revn_n(3,nv)
             ! fill ghost cells in edge
             call fill_ghost_cells(1,1-ngh,0,is_swap_ij,is_rev_p1,is_rev_n)
             ! free temporary array
             deallocate(xgcn,ygcn)
          endif
       endif

       ! edge imin-jmax
       ! --------------
       if ((bl(n)%BC(1)>0).and.(bl(n)%BC(4)>0)) then
          ! start from jmax
          nv=bl(n)%BC(4)
          nv2=bl(nv)%BC(1)
          ! start from imin
          nv=bl(n)%BC(1)
          nv1=bl(nv)%BC(4)
          if (nv1.ne.nv2) then
             !print 100,nv1,nv2
             ! construct artificial edge imin-jmax
             call cons_grid(1,2)
          else
             ! read grid in neighboring block nv
             ! call read_neighbor(trim(gridname)//trim(numchar(nv1))//'.x')
             call read_neighbor('grid_bl'//trim(numchar(nv))//'.bin',nv)
             ! determine swap/reverse for corner block
             is_swap_ij=is_swapij2_bl(1).or.is_swapn_ij(4,nv)
             is_rev_n  = is_rev2_bl(1,2).or.is_revn_p1(4,nv)
             is_rev_p1 = is_rev2_bl(1,1).or.is_revn_n(4,nv)
             ! fill ghost cells in edge
             call fill_ghost_cells(1,ngy+1,ngy+ngh,is_swap_ij,is_rev_p1,is_rev_n)
             ! free temporary array
             deallocate(xgcn,ygcn)
          endif
       endif

       ! edge imax-jmin
       ! --------------
       if ((bl(n)%BC(2)>0).and.(bl(n)%BC(3)>0)) then
          ! start from jmin
          nv=bl(n)%BC(3)
          nv2=bl(nv)%BC(2)
          ! start from imax
          nv=bl(n)%BC(2)
          nv1=bl(nv)%BC(3)
          if (nv1.ne.nv2) then
             !print 100,nv1,nv2
             ! construct artificial edge imax-jmin
             call cons_grid(2,1)
          else
             ! read grid in neighboring block nv
             ! call read_neighbor(trim(gridname)//trim(numchar(nv1))//'.x')
             call read_neighbor('grid_bl'//trim(numchar(nv))//'.bin',nv)
             ! determine swap/reverse for corner block
             is_swap_ij=is_swapij2_bl(2).or.is_swapn_ij(3,nv)
             is_rev_n  = is_rev2_bl(2,2).or.is_revn_p1(3,nv)
             is_rev_p1 = is_rev2_bl(2,1).or.is_revn_n(3,nv)
             ! fill ghost cells in edge
             call fill_ghost_cells(2,1-ngh,0,is_swap_ij,is_rev_p1,is_rev_n)
             ! free temporary array
             deallocate(xgcn,ygcn)
          endif
       endif

       ! edge imax-jmax
       ! --------------
       if ((bl(n)%BC(2)>0).and.(bl(n)%BC(4)>0)) then
          ! start from jmax
          nv=bl(n)%BC(4)
          nv2=bl(nv)%BC(2)
          ! start from imax
          nv=bl(n)%BC(2)
          nv1=bl(nv)%BC(4)
          if (nv1.ne.nv2) then
             !print 100,nv1,nv2
             ! construct artificial edge imax-jmax
             call cons_grid(2,2)
          else
             ! read grid in neighboring block nv
             ! call read_neighbor(trim(gridname)//trim(numchar(nv1))//'.x')
             call read_neighbor('grid_bl'//trim(numchar(nv))//'.bin',nv)
             ! determine swap/reverse for corner block
             is_swap_ij=is_swapij2_bl(2).or.is_swapn_ij(4,nv)
             is_rev_n  = is_rev2_bl(2,2).or.is_revn_p1(4,nv)
             is_rev_p1 = is_rev2_bl(2,1).or.is_revn_n(4,nv)
             ! fill ghost cells in edge
             call fill_ghost_cells(2,ngy+1,ngy+ngh,is_swap_ij,is_rev_p1,is_rev_n)
             ! free temporary array
             deallocate(xgcn,ygcn)
          endif
       endif

       ! Extended dimensions
       ! -------------------
       ni1=1
       ni2=ngx
       nj1=1
       nj2=ngy
       if (bl(n)%BC(1)>0) ni1=1-ngh
       if (bl(n)%BC(2)>0) ni2=ngx+ngh
       if (bl(n)%BC(3)>0) nj1=1-ngh
       if (bl(n)%BC(4)>0) nj2=ngy+ngh

       ! Periodicity
       ! -----------
       if (bl(n)%periodic(1)) then ! at imin
          if (verbose) print *,'translation periodicity imin for block ',n
          ! apply translation of first rows of ghost cells
          do j=nj1,nj2
             xgce(1-ngh:0,j)=xgce(ngx-ngh+1:ngx,j)+Lxp
             ygce(1-ngh:0,j)=ygce(ngx-ngh+1:ngx,j)+Lyp
          enddo
       endif
       if (bl(n)%periodic(2)) then ! at imax
          if (verbose) print *,'translation periodicity imax for block ',n
          ! apply translation of last rows of ghost cells
          do j=nj1,nj2
             xgce(ngx+1:ngx+ngh,j)=xgce(1:ngh,j)+Lxp
             ygce(ngx+1:ngx+ngh,j)=ygce(1:ngh,j)+Lyp
          enddo
       endif
       if (bl(n)%periodic(3)) then ! at jmin
          print *,'translation periodicity jmin for block ',n
          ! apply translation of first rows of ghost cells
          do i=ni1,ni2
             xgce(i,1-ngh:0)=xgce(i,1-ngh:0)+Lxp
             ygce(i,1-ngh:0)=ygce(i,1-ngh:0)+Lyp
          enddo
       endif
       if (bl(n)%periodic(4)) then ! at jmax
          if (verbose) print *,'translation periodicity jmax for block ',n
          ! apply translation of last rows of ghost cells
          do i=ni1,ni2
             xgce(i,ngy+1:ngy+ngh)=xgce(i,ngy+1:ngy+ngh)+Lxp
             ygce(i,ngy+1:ngy+ngh)=ygce(i,ngy+1:ngy+ngh)+Lyp
          enddo
       endif
       if (ngz.gt.1) then
          if (bl(n)%periodic(5)) then ! at kmin
             if (verbose) print *,'translation periodicity kmin for block ',n
             ! apply translation of first rows of ghost cells
             zmin=minval(zg(1:ngz),1)
             zmax=maxval(zg(1:ngz),1)
             Lzp0=zmax-zmin+abs(zg(2)-zg(1))
             zg(1-ngh:0)=zg(ngz-ngh+1:ngz)-Lzp0
          endif
          if (bl(n)%periodic(6)) then ! at kmax
             if (verbose) print *,'translation periodicity kmax for block ',n
             ! apply translation of last rows of ghost cells
             zmin=minval(zg(1:ngz),1)
             zmax=maxval(zg(1:ngz),1)
             if (Lzp.eq.0.0_wp) then
                Lzp0=zmax-zmin+abs(zg(2)-zg(1))
             else
                Lzp0=Lzp
             endif
             zg(ngz+1:ngz+ngh)=zg(1:ngh)+Lzp0
          endif
       endif

       ! free grid arrays
       ! ----------------
       deallocate(is_swapn_ij,is_revn_n,is_revn_p1)

       ! Write extended global grid
       ! --------------------------
       gridname='grid_bl'//trim(numchar(nob(iproc)))//'_ngh'//trim(numchar(ngh))//'.bin'
       ! init MPI-IO for grid write global grid
       ! & read extended grid
       ! & free MPI IO memory
       if (is_curv) then
          call read_write_grid_curv(gridname,ngh,WRITE_LEADER)
       else
          xg(1-ngh:ngx+ngh) = xgce(1-ngh:ngx+ngh,1)
          yg(1-ngh:ngy+ngh) = ygce(1,1-ngh:ngy+ngh)
          call read_write_grid_cart(gridname,ngh,WRITE_LEADER)
          deallocate(xg,yg)
       endif

       ! Deallocation
       ! ------------
       deallocate(xgce,ygce,zg)
    endif

  end subroutine add_ghost_cells_grid

  !===============================================================================
  subroutine read_neighbor(gridname,nv)
  !===============================================================================
    !> read grid in neighboring block nv
  !===============================================================================
    use mod_mpi
    use mod_block
    use mod_utils
    use warnstop
    implicit none
    ! ----------------------------------------------------------------------------
    integer, intent(in) :: nv
    character(len=*), intent(in) :: gridname
    ! Local variables IO
    integer :: type_mat_cartx,type_mat_viewx
    integer :: type_mat_carty,type_mat_viewy
    integer :: type_mat_curv,type_mat_view
    integer :: ierror,nbytes_dbl,fh
    integer(kind=MPI_OFFSET_KIND) :: offset ! file offset
    integer, dimension(1) :: shape_grx,shape_gry,start_grx
    integer, dimension(1) :: shape_globx,shape_localx,start_localx
    integer, dimension(1) :: shape_globy,shape_localy
    integer, dimension(2) :: shape_gr,start_gr
    integer, dimension(2) :: shape_glob,shape_local,start_local
    integer, dimension(MPI_STATUS_SIZE) :: statut
    ! Local variables read
    logical  :: iexist
    ! ----------------------------------------------------------------------------
    integer :: i,j,mb
    ! ----------------------------------------------------------------------------

    ! Get size of the type MPI_DOUBLE_PRECISION
    call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,nbytes_dbl,ierror)

    ! Size of neighboring block
    ngxn = bl(nv)%ni; ngyn = bl(nv)%nj

    ! Allocate grid arrays
    allocate(xgcn(ngxn,ngyn),ygcn(ngxn,ngyn))

    if (is_curv) then
       ! Creation of derived type type_mat_curv
       shape_gr= (/ ngxn, ngyn /); start_gr= (/ 0, 0 /)
       call MPI_TYPE_CREATE_SUBARRAY(2,shape_gr,shape_gr,start_gr,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,type_mat_curv,ierror)
       call MPI_TYPE_COMMIT(type_mat_curv,ierror)
       ! Creation of derived type type_mat_view
       shape_glob=(/ ngxn, ngyn /); shape_local=(/ ngxn, ngyn /); start_local=(/ 0, 0 /)
       call MPI_TYPE_CREATE_SUBARRAY(2,shape_glob,shape_local,start_local,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,type_mat_view,ierror)
       call MPI_TYPE_COMMIT(type_mat_view,ierror)
       ! Open file
       call MPI_FILE_OPEN(COMM_r_grid,trim(gridname),MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierror)
       if (ierror/=MPI_SUCCESS) call mpistop('MPI_IO: Error in read_neighbor in opening gridfile '// trim(gridname),1)
       ! Setting the view on the file
       call MPI_FILE_SET_VIEW(fh,0_MPI_OFFSET_KIND,MPI_DOUBLE_PRECISION,type_mat_view,'native',MPI_INFO_NULL,ierror)
       ! Reading xgcn
       call MPI_FILE_READ_ALL(fh,xgcn(1:ngxn,1:ngyn),1,type_mat_curv,statut,ierror)
       if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_neighbor xgcn in mod_add_gh for file '//trim(gridname),1)
       ! Reading ygcn
       call MPI_FILE_READ_ALL(fh,ygcn(1:ngxn,1:ngyn),1,type_mat_curv,statut,ierror)
       if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_neighbor ygcn in mod_add_gh for file '//trim(gridname),1)
       ! Close file
       call MPI_FILE_CLOSE(fh,ierror)
       ! Free types for MPI-IO of grid
       call MPI_TYPE_FREE(type_mat_curv,info); call MPI_TYPE_FREE(type_mat_view,info)
    else
       ! Creation of derived type type_mat_curv
       shape_grx= (/ ngxn /); shape_gry= (/ ngyn /); start_grx= (/ 0 /)
       print *,shape_grx,shape_gry,start_grx
       call MPI_TYPE_CREATE_SUBARRAY(1,shape_grx,shape_grx,start_grx,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,type_mat_cartx,ierror)
       call MPI_TYPE_COMMIT(type_mat_cartx,ierror)
       call MPI_TYPE_CREATE_SUBARRAY(1,shape_gry,shape_gry,start_grx,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,type_mat_carty,ierror)
       call MPI_TYPE_COMMIT(type_mat_carty,ierror)
       ! Creation of derived type type_mat_view
       shape_globx=(/ ngxn /); shape_localx=(/ ngxn /); start_localx=(/ 0 /)
       call MPI_TYPE_CREATE_SUBARRAY(1,shape_globx,shape_localx,start_localx,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,type_mat_viewx,ierror)
       call MPI_TYPE_COMMIT(type_mat_viewx,ierror)
       shape_globy=(/ ngyn /); shape_localy=(/ ngyn /)
       call MPI_TYPE_CREATE_SUBARRAY(1,shape_globy,shape_localy,start_localx,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,type_mat_viewy,ierror)
       call MPI_TYPE_COMMIT(type_mat_viewy,ierror)
       ! Open file
       call MPI_FILE_OPEN(COMM_r_grid,trim(gridname),MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierror)
       if (ierror/=MPI_SUCCESS) call mpistop('MPI_IO: Error in read_neighbor in opening gridfile '// trim(gridname),1)
       ! Reading xgcn
       call MPI_FILE_SET_VIEW(fh,0_MPI_OFFSET_KIND,MPI_DOUBLE_PRECISION,type_mat_viewx,'native',MPI_INFO_NULL,ierror)
       call MPI_FILE_READ_ALL(fh,xgcn(1:ngxn,1),1,type_mat_cartx,statut,ierror)
       if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_neighbor xgcn in mod_add_gh for file '//trim(gridname),1)
       do j=2,ngyn
          xgcn(1:ngxn,j) = xgcn(1:ngxn,1)
       enddo
       ! Reading ygcn
       offset = 8*ngxn ! Setting the offset
       call MPI_FILE_SET_VIEW(fh,offset,MPI_DOUBLE_PRECISION,type_mat_viewy,'native',MPI_INFO_NULL,ierror)
       call MPI_FILE_READ_ALL(fh,ygcn(1,1:ngyn),1,type_mat_carty,statut,ierror)
       if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_neighbor ygcn in mod_add_gh for file '//trim(gridname),1)
       do i=2,ngxn
          ygcn(i,1:ngyn) = ygcn(1,1:ngyn)
       enddo
       ! Close file
       call MPI_FILE_CLOSE(fh,ierror)
       ! Free types for MPI-IO of grid
       call MPI_TYPE_FREE(type_mat_cartx,info); call MPI_TYPE_FREE(type_mat_viewx,info)
       call MPI_TYPE_FREE(type_mat_carty,info); call MPI_TYPE_FREE(type_mat_viewy,info)
    endif

  end subroutine read_neighbor

  !===============================================================================
  subroutine cons_grid(ind2,ind3)
  !===============================================================================
    !> construct artificial grid for egdes of curvilinear or cartesian grids
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer,intent(in) :: ind2,ind3
    ! ----------------------------------------------------------------------------
    integer :: i,j,i1,j1,ii,jj
    real(wp) :: x1,x2,x3,y1,y2,y3
    ! ----------------------------------------------------------------------------

    ! Definition of edges
    ! ===================
    ! 1,1 : imin-jmin
    ! 1,2 : imin-jmax
    ! 2,1 : imax-jmin
    ! 2,2 : imax-jmax
        
    ! origin point
    ! ------------
    if (ind2==1) i1=1
    if (ind2==2) i1=ngx
    if (ind3==1) j1=1
    if (ind3==2) j1=ngy
    x1=xgce(i1,j1)
    y1=ygce(i1,j1)

    do i=1,ngh
       do j=1,ngh
          if (ind2==1) ii=i-ngh
          if (ind2==2) ii=ngx+i
          if (ind3==1) jj=j-ngh
          if (ind3==2) jj=ngy+j
          ! corner point #1
          x2=xgce(ii,j1)
          y2=ygce(ii,j1)
          ! corner point #2
          x3=xgce(i1,jj)
          y3=ygce(i1,jj)
          ! reconstruction
          xgce(ii,jj)=x2+x3-x1
          ygce(ii,jj)=y2+y3-y1
       enddo
    enddo

  end subroutine cons_grid


  !===============================================================================
  subroutine fill_ghost_cells(n_f,n_1,n_2,is_swap_ij,is_rev_p1,is_rev_n)
  !===============================================================================
    !> fill ghost cells from neighboring block
  !===============================================================================
    use mod_mpi
    implicit none
    ! ----------------------------------------------------------------------------
    integer, intent(in) :: n_f ! face number
    integer, intent(in) :: n_1,n_2 ! sizes of parallel direction
    logical, intent(in) :: is_swap_ij ! swap direction
    logical, intent(in) :: is_rev_p1,is_rev_n ! reverse directions
    ! ----------------------------------------------------------------------------
    integer :: i,j,i2,j2
    ! ----------------------------------------------------------------------------

    ! Indices n_1 & n_2 are used for parallel direction #1
    ! ====================================================
    ! for imin- imax-: n_1 n_2 is direction j
    ! for jmin- jmax-: n_1 n_2 is direction i

    ! Fill ghost cells for block faces
    ! ================================

    ! face imin
    ! ---------
    select case (n_f)

    case (1)
    if (is_swap_ij) then
          do j=n_1,n_2
             do i=1-ngh,0
                j2=ngyn+i
                if (n_1<0) then
                   i2=ngxn+j
                elseif (n_1>ngy) then
                   i2=j-ngy
                else
                   i2=j
                endif
                ! reverse normal
                if (is_rev_n) j2=1-i
                ! reverse parallel 1
                if (is_rev_p1) then
                   if (n_1<0) then
                      i2=1-j
                   elseif (n_1>ngy) then
                      !i2=2*ngxn-j+1
                      i2=ngxn-i2+1
                   else
                      i2=ngxn-j+1
                   endif
                endif
                xgce(i,j)=xgcn(i2,j2)
                ygce(i,j)=ygcn(i2,j2)
             enddo
          enddo
       else
          do j=n_1,n_2
             do i=1-ngh,0
                i2=ngxn+i
                if (n_1<0) then
                   j2=ngyn+j
                elseif (n_1>ngy) then
                   j2=j-ngy
                else
                   j2=j
                endif
                ! reverse normal
                if (is_rev_n) i2=1-i
                ! reverse parallel 1
                if (is_rev_p1) then
                   if (n_1<0) then
                      j2=1-j
                   elseif (n_1>ngy) then
                      !j2=2*ngyn-j+1
                      j2=ngyn-j2+1
                   else
                      j2=ngyn-j+1
                   endif
                endif
                xgce(i,j)=xgcn(i2,j2)
                ygce(i,j)=ygcn(i2,j2)
             enddo
          enddo
       endif

    ! face imax
    ! ---------
    case (2)
       if (is_swap_ij) then
          do j=n_1,n_2
             do i=ngx+1,ngx+ngh
                j2=i-ngx
                if (n_1<0) then
                   i2=ngxn+j
                elseif (n_1>ngy) then
                   i2=j-ngy
                else
                   i2=j
                endif
                ! reverse normal
                !if (is_rev_n) j2=2*ngyn-j+1
                if (is_rev_n) j2=ngyn-j2+1
                ! reverse parallel 1
                if (is_rev_p1) then
                   if (n_1<0) then
                      i2=1-j
                   elseif (n_1>ngy) then
                      !i2=2*ngxn-j+1
                      i2=ngxn-i2+1
                   else
                      i2=ngxn-j+1
                   endif
                endif
                xgce(i,j)=xgcn(i2,j2)
                ygce(i,j)=ygcn(i2,j2)
             enddo
          enddo
       else
          do j=n_1,n_2
             do i=ngx+1,ngx+ngh
                i2=i-ngx
                if (n_1<0) then
                   j2=ngyn+j
                elseif (n_1>ngy) then
                   j2=j-ngy
                else
                   j2=j
                endif
                ! reverse normal
                !if (is_rev_n) i2=2*ngxn-i+1
                if (is_rev_n) i2=ngxn-i2+1
                ! reverse parallel 1
                if (is_rev_p1) then
                   if (n_1<0) then
                      j2=1-j
                   elseif (n_1>ngy) then
                      !j2=2*ngyn-j+1
                      j2=ngyn-j2+1
                   else
                      j2=ngyn-j+1
                   endif
                endif
                xgce(i,j)=xgcn(i2,j2)
                ygce(i,j)=ygcn(i2,j2)
             enddo
          enddo
       endif

    ! face jmin
    ! ---------
    case (3)
       if (is_swap_ij) then
          do j=1-ngh,0
             do i=n_1,n_2
                if (n_1<0) then
                   j2=ngyn+i
                elseif (n_1>ngx) then
                   j2=i-ngx
                else
                   j2=i
                endif
                i2=ngxn+j
                ! reverse normal
                if (is_rev_n) i2=1-j
                ! reverse parallel 1
                if (is_rev_p1) then
                   if (n_1<0) then
                      j2=1-i
                   elseif (n_1>ngx) then
                      !j2=2*ngyn-i+1
                      j2=ngyn-j2+1
                   else
                      j2=ngyn-i+1
                   endif
                endif
                xgce(i,j)=xgcn(i2,j2)
                ygce(i,j)=ygcn(i2,j2)
             enddo
          enddo
       else
          do j=1-ngh,0
             do i=n_1,n_2
                if (n_1<0) then
                   i2=ngxn+i
                elseif (n_1>ngx) then
                   i2=i-ngx
                else
                   i2=i
                endif
                j2=ngyn+j
                ! reverse normal
                if (is_rev_n) j2=1-j
                ! reverse parallel 1
                if (is_rev_p1) then
                   if (n_1<0) then
                      i2=1-i
                   elseif (n_1>ngx) then
                      !i2=2*ngxn-i+1
                      i2=ngxn-i2+1
                   else
                      i2=ngxn-i+1
                   endif
                endif
                xgce(i,j)=xgcn(i2,j2)
                ygce(i,j)=ygcn(i2,j2)
             enddo
          enddo
       endif

    ! face jmax
    ! ---------
    case (4)
       if (is_swap_ij) then
          do j=ngy+1,ngy+ngh
             do i=n_1,n_2
                if (n_1<0) then
                   j2=ngyn+i
                elseif (n_1>ngx) then
                   j2=i-ngx
                else
                   j2=i
                endif
                i2=j-ngy
                ! reverse normal
                !if (is_rev_n) i2=2*ngxn-j+1
                if (is_rev_n) i2=ngxn-i2+1
                ! reverse parallel 1
                if (is_rev_p1) then
                   if (n_1<0) then
                      j2=1-i
                   elseif (n_1>ngx) then
                      !j2=2*ngyn-i+1
                      j2=ngyn-j2+1
                   else
                      j2=ngyn-i+1
                   endif
                endif
                xgce(i,j)=xgcn(i2,j2)
                ygce(i,j)=ygcn(i2,j2)
             enddo
          enddo
       else
          do j=ngy+1,ngy+ngh
             do i=n_1,n_2
                if (n_1<0) then
                   i2=ngxn+i
                elseif (n_1>ngx) then
                   i2=i-ngx
                else
                   i2=i
                endif
                j2=j-ngy
                ! reverse normal
                !if (is_rev_n) j2=2*ngyn-j+1
                if (is_rev_n) j2=ngyn-j2+1
                ! reverse parallel 1
                if (is_rev_p1) then
                   if (n_1<0) then
                      i2=1-i
                   elseif (n_1>ngx) then
                      !i2=2*ngxn-i+1
                      i2=ngxn-i2+1
                   else
                      i2=ngxn-i+1
                   endif
                endif
                xgce(i,j)=xgcn(i2,j2)
                ygce(i,j)=ygcn(i2,j2)
             enddo
          enddo
       endif

    end select

  end subroutine fill_ghost_cells

end module mod_add_gh
