!=================================================================================
module mod_add_gh3d
!=================================================================================
  !> Module to generate extended grid for 3D curvilinear grid
  !> Created with some adapted NASA grid utilities routines
!=================================================================================
  use mod_grid
  implicit none
  ! ----------------------------------------------------------------------------
  ! current grid
  real(wp), dimension(:,:,:), allocatable, private :: xgc31,ygc31,zgc31
  ! neighbor grid
  integer, private :: ngxn,ngyn,ngzn
  real(wp), dimension(:,:,:), allocatable, private :: xgc3n,ygc3n,zgc3n
  ! ----------------------------------------------------------------------------

contains
  
  !===============================================================================
  subroutine add_ghost_cells_grid3d
  !===============================================================================
    !> add ghost cells at block junctions
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
    integer :: ni,nj,nk,ni1,nj1,nk1,ni2,nj2,nk2
    integer :: n,nv,nv_,mb
    integer :: nv1,nv2,nv3,nv4,nv5,nv6
    logical :: is_swap_ij,is_swap_ik,is_swap_jk
    logical :: is_rev_n,is_rev_p1,is_rev_p2
    logical, dimension(:,:), allocatable :: is_swapn_ij,is_swapn_ik,is_swapn_jk
    logical, dimension(:,:), allocatable :: is_revn_n,is_revn_p1,is_revn_p2
    real(wp), dimension(:), allocatable :: y_p,z_p
    real(wp) :: Lzp0
    character(len=200) :: gridname
    ! ----------------------------------------------------------------------------

    ! Print message
    ! =============
    if ((iproc==0).and.(verbose)) then
       print *,repeat('=',70)
       print *,'Add ghost-cells at block junctions (modif oct 2023)'
       print *,repeat('=',70)
    endif

    ! Get swap/reverse indicators for all blocks
    ! ==========================================
    allocate(is_swapn_ij(6,nbloc),is_swapn_ik(6,nbloc),is_swapn_jk(6,nbloc))
    allocate(is_revn_n(6,nbloc),is_revn_p1(6,nbloc),is_revn_p2(6,nbloc))
    do k=1,6
       call MPI_ALLGATHER(is_swapij2_bl(k),1,MPI_LOGICAL,is_swapn_ij(k,:),1,MPI_LOGICAL,COMM_interblock,info)
       call MPI_ALLGATHER(is_swapik2_bl(k),1,MPI_LOGICAL,is_swapn_ik(k,:),1,MPI_LOGICAL,COMM_interblock,info)
       call MPI_ALLGATHER(is_swapjk2_bl(k),1,MPI_LOGICAL,is_swapn_jk(k,:),1,MPI_LOGICAL,COMM_interblock,info)
       call MPI_ALLGATHER(is_rev2_bl(k,2),1,MPI_LOGICAL,is_revn_n(k,:),1,MPI_LOGICAL,COMM_interblock,info)
       call MPI_ALLGATHER(is_rev2_bl(k,1),1,MPI_LOGICAL,is_revn_p1(k,:),1,MPI_LOGICAL,COMM_interblock,info)
       call MPI_ALLGATHER(is_rev2_bl(k,3),1,MPI_LOGICAL,is_revn_p2(k,:),1,MPI_LOGICAL,COMM_interblock,info)
    enddo
 
    ! current block number
    ! -------------------
    n=nob(iproc)
    
    !if (iproc==0) then
    if (iproc==iproc_leader(n)) then

       ! ! Interior points of current block  <~ Already in xgc3e !!
       ! ! ================================
       
       ! ! read grid without ghost cells
       ! ! -----------------------------
       ! if (is_adjoint_block) then
       !    gridname=trim(dirGRID)//'/'//trim(nameGRID)//'_bl'
       ! else
       !    gridname=trim(dirGRID)//'/'//trim(nameGRID)//'_mod_bl'
       ! endif
       
       ! open(50,file=trim(gridname)//trim(numchar(n))//'.x',form='formatted')
       ! rewind(50)
       ! read(50,*) mb
       ! ! block dimensions
       ! read(50,*) ngx,ngy,ngz
       ! allocate(xgc31(ngx,ngy,ngz),ygc31(ngx,ngy,ngz),zgc31(ngx,ngy,ngz))
       ! ! grid
       ! read(50,*) (((xgc31(i,j,k),i=1,ngx),j=1,ngy),k=1,ngz)
       ! read(50,*) (((ygc31(i,j,k),i=1,ngx),j=1,ngy),k=1,ngz)
       ! read(50,*) (((zgc31(i,j,k),i=1,ngx),j=1,ngy),k=1,ngz)
       ! close(50)
       
       ! ! extended grid with ghost cells
       ! ! ------------------------------
       ! allocate(xgc3e(1-ngh:ngx+ngh,1-ngh:ngy+ngh,1-ngh:ngz+ngh))
       ! allocate(ygc3e(1-ngh:ngx+ngh,1-ngh:ngy+ngh,1-ngh:ngz+ngh))
       ! allocate(zgc3e(1-ngh:ngx+ngh,1-ngh:ngy+ngh,1-ngh:ngz+ngh))
       ! xgc3e(1:ngx,1:ngy,1:ngz)=xgc31
       ! ygc3e(1:ngx,1:ngy,1:ngz)=ygc31
       ! zgc3e(1:ngx,1:ngy,1:ngz)=zgc31

       ! Fill ghost cells for block faces
       ! ================================

       ! face imin
       ! ---------
       if (bl(n)%BC(1)>0) then
          ! neighbor
          nv=bl(n)%BC(1)
          ! read grid in neighboring block nv
          call read_neighbor3d(trim(gridname)//trim(numchar(nv))//'.x')
          ! fill ghost cells in face #1
          call fill_ghost_cells3d(1,1,ngy,1,ngz, &
               is_swapij2_bl(1),is_swapik2_bl(1),is_swapjk2_bl(1), &
               is_rev2_bl(1,1),is_rev2_bl(1,2),is_rev2_bl(1,3))
          ! free temporary array
          deallocate(xgc3n,ygc3n,zgc3n)
       endif

       ! face imax
       ! ---------
       if (bl(n)%BC(2)>0) then
          ! neighbor
          nv=bl(n)%BC(2)
          ! read grid in neighboring block nv
          call read_neighbor3d(trim(gridname)//trim(numchar(nv))//'.x')
          ! fill ghost cells in face #2
          call fill_ghost_cells3d(2,1,ngy,1,ngz, &
               is_swapij2_bl(2),is_swapik2_bl(2),is_swapjk2_bl(2), &
               is_rev2_bl(2,1),is_rev2_bl(2,2),is_rev2_bl(2,3))
          ! free temporary array
          deallocate(xgc3n,ygc3n,zgc3n)
       endif

       ! face jmin
       ! ---------
       if (bl(n)%BC(3)>0) then
          ! neighbor
          nv=bl(n)%BC(3)
          ! read grid in neighboring block nv
          call read_neighbor3d(trim(gridname)//trim(numchar(nv))//'.x')
          ! fill ghost cells in face #3
          call fill_ghost_cells3d(3,1,ngx,1,ngz, &
               is_swapij2_bl(3),is_swapik2_bl(3),is_swapjk2_bl(3), &
               is_rev2_bl(3,1),is_rev2_bl(3,2),is_rev2_bl(3,3))
          ! free temporary array
          deallocate(xgc3n,ygc3n,zgc3n)
       endif

       ! face jmax
       ! ---------
       if (bl(n)%BC(4)>0) then
          ! neighbor
          nv=bl(n)%BC(4)
          ! read grid in neighboring block nv
          call read_neighbor3d(trim(gridname)//trim(numchar(nv))//'.x')
          ! fill ghost cells in face #4
          call fill_ghost_cells3d(4,1,ngx,1,ngz, &
               is_swapij2_bl(4),is_swapik2_bl(4),is_swapjk2_bl(4), &
               is_rev2_bl(4,1),is_rev2_bl(4,2),is_rev2_bl(4,3))
          ! free temporary array
          deallocate(xgc3n,ygc3n,zgc3n)
       endif

       ! face kmin
       ! ---------
       if (bl(n)%BC(5)>0) then
          ! neighbor
          nv=bl(n)%BC(5)
          ! read grid in neighboring block nv
          call read_neighbor3d(trim(gridname)//trim(numchar(nv))//'.x')
          ! fill ghost cells in face #5
          call fill_ghost_cells3d(5,1,ngx,1,ngy, &
               is_swapij2_bl(5),is_swapik2_bl(5),is_swapjk2_bl(5), &
               is_rev2_bl(5,1),is_rev2_bl(5,2),is_rev2_bl(5,3))
          ! free temporary array
          deallocate(xgc3n,ygc3n,zgc3n)
       endif

       ! face kmax
       ! ---------
       if (bl(n)%BC(6)>0) then
          ! neighbor
          nv=bl(n)%BC(6)
          ! read grid in neighboring block nv
          call read_neighbor3d(trim(gridname)//trim(numchar(nv))//'.x')
          ! fill ghost cells in face #6
          call fill_ghost_cells3d(6,1,ngx,1,ngy, &
               is_swapij2_bl(6),is_swapik2_bl(6),is_swapjk2_bl(6), &
               is_rev2_bl(6,1),is_rev2_bl(6,2),is_rev2_bl(6,3))    
          ! free temporary array
          deallocate(xgc3n,ygc3n,zgc3n)
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
             !print 100,nv1,nv2
             ! construct artificial edge imin-jmin
             call cons_grid3d(1,1,1,.true.)
          else
             ! read grid in neighboring block nv
             call read_neighbor3d(trim(gridname)//trim(numchar(nv1))//'.x')
             ! determine swap/reverse for corner block
             is_swap_ij=is_swapij2_bl(1).or.is_swapn_ij(3,nv)
             is_swap_ik=is_swapik2_bl(1).or.is_swapn_ik(3,nv)
             is_swap_jk=is_swapjk2_bl(1).or.is_swapn_jk(3,nv)
             is_rev_n  = is_rev2_bl(1,2).or.is_revn_p1(3,nv)
             is_rev_p1 = is_rev2_bl(1,1).or.is_revn_n(3,nv)
             is_rev_p2 = is_rev2_bl(1,3).or.is_revn_p2(3,nv)
             ! fill ghost cells in edge
             call fill_ghost_cells3d(1,1-ngh,0,1,ngz, &
                  is_swap_ij,is_swap_ik,is_swap_jk, &
                  is_rev_p1,is_rev_n,is_rev_p2)
             ! free temporary array
             deallocate(xgc3n,ygc3n,zgc3n)
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
             call cons_grid3d(1,1,2,.true.)
          else
             ! read grid in neighboring block nv
             call read_neighbor3d(trim(gridname)//trim(numchar(nv1))//'.x')
             ! determine swap/reverse for corner block
             is_swap_ij=is_swapij2_bl(1).or.is_swapn_ij(4,nv)
             is_swap_ik=is_swapik2_bl(1).or.is_swapn_ik(4,nv)
             is_swap_jk=is_swapjk2_bl(1).or.is_swapn_jk(4,nv)
             is_rev_n  = is_rev2_bl(1,2).or.is_revn_p1(4,nv)
             is_rev_p1 = is_rev2_bl(1,1).or.is_revn_n(4,nv)
             is_rev_p2 = is_rev2_bl(1,3).or.is_revn_p2(4,nv)
             ! fill ghost cells in edge
             call fill_ghost_cells3d(1,ngy+1,ngy+ngh,1,ngz, &
                  is_swap_ij,is_swap_ik,is_swap_jk, &
                  is_rev_p1,is_rev_n,is_rev_p2)
             ! free temporary array
             deallocate(xgc3n,ygc3n,zgc3n)
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
             call cons_grid3d(1,2,1,.true.)
          else
             ! read grid in neighboring block nv
             call read_neighbor3d(trim(gridname)//trim(numchar(nv1))//'.x')
             ! determine swap/reverse for corner block
             is_swap_ij=is_swapij2_bl(2).or.is_swapn_ij(3,nv)
             is_swap_ik=is_swapik2_bl(2).or.is_swapn_ik(3,nv)
             is_swap_jk=is_swapjk2_bl(2).or.is_swapn_jk(3,nv)
             is_rev_n  = is_rev2_bl(2,2).or.is_revn_p1(3,nv)
             is_rev_p1 = is_rev2_bl(2,1).or.is_revn_n(3,nv)
             is_rev_p2 = is_rev2_bl(2,3).or.is_revn_p2(3,nv)
             ! fill ghost cells in edge
             call fill_ghost_cells3d(2,1-ngh,0,1,ngz, &
                  is_swap_ij,is_swap_ik,is_swap_jk, &
                  is_rev_p1,is_rev_n,is_rev_p2)
             ! free temporary array
             deallocate(xgc3n,ygc3n,zgc3n)
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
             call cons_grid3d(1,2,2,.true.)
          else
             ! read grid in neighboring block nv
             call read_neighbor3d(trim(gridname)//trim(numchar(nv1))//'.x')
             ! determine swap/reverse for corner block
             is_swap_ij=is_swapij2_bl(2).or.is_swapn_ij(4,nv)
             is_swap_ik=is_swapik2_bl(2).or.is_swapn_ik(4,nv)
             is_swap_jk=is_swapjk2_bl(2).or.is_swapn_jk(4,nv)
             is_rev_n  = is_rev2_bl(2,2).or.is_revn_p1(4,nv)
             is_rev_p1 = is_rev2_bl(2,1).or.is_revn_n(4,nv)
             is_rev_p2 = is_rev2_bl(2,3).or.is_revn_p2(4,nv)
             ! fill ghost cells in edge
             call fill_ghost_cells3d(2,ngy+1,ngy+ngh,1,ngz, &
                  is_swap_ij,is_swap_ik,is_swap_jk, &
                  is_rev_p1,is_rev_n,is_rev_p2)
             ! free temporary array
             deallocate(xgc3n,ygc3n,zgc3n)
          endif
       endif

       ! edge jmin-kmin
       ! --------------
       if ((bl(n)%BC(3)>0).and.(bl(n)%BC(5)>0)) then
          ! start from kmin
          nv=bl(n)%BC(5)
          nv2=bl(nv)%BC(3)
          ! start from jmin
          nv=bl(n)%BC(3)
          nv1=bl(nv)%BC(5)
          if (nv1.ne.nv2) then
             !print 100,nv1,nv2
             ! construct artificial edge jmin-kmin
             call cons_grid3d(2,1,1,.true.)
          else
             ! read grid in neighboring block nv
             call read_neighbor3d(trim(gridname)//trim(numchar(nv1))//'.x')
             ! determine swap/reverse for corner block
             is_swap_ij=is_swapij2_bl(3).or.is_swapn_ij(5,nv)
             is_swap_ik=is_swapik2_bl(3).or.is_swapn_ik(5,nv)
             is_swap_jk=is_swapjk2_bl(3).or.is_swapn_jk(5,nv)
             is_rev_n  = is_rev2_bl(3,2).or.is_revn_p2(5,nv)
             is_rev_p1 = is_rev2_bl(3,1).or.is_revn_p1(5,nv)
             is_rev_p2 = is_rev2_bl(3,3).or.is_revn_n(5,nv)
             ! fill ghost cells in edge
             call fill_ghost_cells3d(3,1,ngx,1-ngh,0, &
                  is_swap_ij,is_swap_ik,is_swap_jk, &
                  is_rev_p1,is_rev_n,is_rev_p2)
             ! free temporary array
             deallocate(xgc3n,ygc3n,zgc3n)
          endif
       endif

       ! edge jmin-kmax
       ! --------------
       if ((bl(n)%BC(3)>0).and.(bl(n)%BC(6)>0)) then
          ! start from kmax
          nv=bl(n)%BC(6)
          nv2=bl(nv)%BC(3)
          ! start from jmin
          nv=bl(n)%BC(3)
          nv1=bl(nv)%BC(6)
          if (nv1.ne.nv2) then
             !print 100,nv1,nv2
             ! construct artificial edge jmin-kmax
             call cons_grid3d(2,1,2,.true.)
          else
             ! read grid in neighboring block nv
             call read_neighbor3d(trim(gridname)//trim(numchar(nv1))//'.x')
             ! determine swap/reverse for corner block
             is_swap_ij=is_swapij2_bl(3).or.is_swapn_ij(6,nv)
             is_swap_ik=is_swapik2_bl(3).or.is_swapn_ik(6,nv)
             is_swap_jk=is_swapjk2_bl(3).or.is_swapn_jk(6,nv)
             is_rev_n  = is_rev2_bl(3,2).or.is_revn_p2(6,nv)
             is_rev_p1 = is_rev2_bl(3,1).or.is_revn_p1(6,nv)
             is_rev_p2 = is_rev2_bl(3,3).or.is_revn_n(6,nv)
             ! fill ghost cells in edge
             call fill_ghost_cells3d(3,1,ngx,ngz+1,ngz+ngh, &
                  is_swap_ij,is_swap_ik,is_swap_jk, &
                  is_rev_p1,is_rev_n,is_rev_p2)
             ! free temporary array
             deallocate(xgc3n,ygc3n,zgc3n)
          endif
       endif

       ! edge jmax-kmin
       ! --------------
       if ((bl(n)%BC(4)>0).and.(bl(n)%BC(5)>0)) then
          ! start from kmin
          nv=bl(n)%BC(5)
          nv2=bl(nv)%BC(4)
          ! start from jmax
          nv=bl(n)%BC(4)
          nv1=bl(nv)%BC(5)
          if (nv1.ne.nv2) then
             !print 100,nv1,nv2
             ! construct artificial edge jmax-kmin
             call cons_grid3d(2,2,1,.true.)
          else
             ! read grid in neighboring block nv
             call read_neighbor3d(trim(gridname)//trim(numchar(nv1))//'.x')
             ! determine swap/reverse for corner block
             is_swap_ij=is_swapij2_bl(4).or.is_swapn_ij(5,nv)
             is_swap_ik=is_swapik2_bl(4).or.is_swapn_ik(5,nv)
             is_swap_jk=is_swapjk2_bl(4).or.is_swapn_jk(5,nv)
             is_rev_n  = is_rev2_bl(4,2).or.is_revn_p2(5,nv)
             is_rev_p1 = is_rev2_bl(4,1).or.is_revn_p1(5,nv)
             is_rev_p2 = is_rev2_bl(4,3).or.is_revn_n(5,nv)
             ! fill ghost cells in edge
             call fill_ghost_cells3d(4,1,ngx,1-ngh,0, &
                  is_swap_ij,is_swap_ik,is_swap_jk, &
                  is_rev_p1,is_rev_n,is_rev_p2)
             ! free temporary array
             deallocate(xgc3n,ygc3n,zgc3n)
          endif
       endif

       ! edge jmax-kmax
       ! --------------
       if ((bl(n)%BC(4)>0).and.(bl(n)%BC(6)>0)) then
          ! start from kmax
          nv=bl(n)%BC(6)
          nv2=bl(nv)%BC(4)
          ! start from jmax
          nv=bl(n)%BC(4)
          nv1=bl(nv)%BC(6)
          if (nv1.ne.nv2) then
             !print 100,nv1,nv2
             ! construct artificial edge jmax-kmax
             call cons_grid3d(2,2,2,.true.)
          else
             ! read grid in neighboring block nv
             call read_neighbor3d(trim(gridname)//trim(numchar(nv1))//'.x')
             ! determine swap/reverse for corner block
             is_swap_ij=is_swapij2_bl(4).or.is_swapn_ij(6,nv)
             is_swap_ik=is_swapik2_bl(4).or.is_swapn_ik(6,nv)
             is_swap_jk=is_swapjk2_bl(4).or.is_swapn_jk(6,nv)
             is_rev_n  = is_rev2_bl(4,2).or.is_revn_p2(6,nv)
             is_rev_p1 = is_rev2_bl(4,1).or.is_revn_p1(6,nv)
             is_rev_p2 = is_rev2_bl(4,3).or.is_revn_n(6,nv)
             ! fill ghost cells in edge
             call fill_ghost_cells3d(4,1,ngx,ngz+1,ngz+ngh, &
                  is_swap_ij,is_swap_ik,is_swap_jk, &
                  is_rev_p1,is_rev_n,is_rev_p2)
             ! free temporary array
             deallocate(xgc3n,ygc3n,zgc3n)
          endif
       endif

       ! edge kmin-imin
       ! --------------
       if ((bl(n)%BC(5)>0).and.(bl(n)%BC(1)>0)) then
          ! start from imin
          nv=bl(n)%BC(1)
          nv2=bl(nv)%BC(5)
          ! start from kmin
          nv=bl(n)%BC(5)
          nv1=bl(nv)%BC(1)
          if (nv1.ne.nv2) then
             !print 100,nv1,nv2
             ! construct artificial edge kmin-imin
             call cons_grid3d(3,1,1,.true.)
          else
             ! read grid in neighboring block nv
             call read_neighbor3d(trim(gridname)//trim(numchar(nv1))//'.x')
             ! determine swap/reverse for corner block
             is_swap_ij=is_swapij2_bl(5).or.is_swapn_ij(1,nv)
             is_swap_ik=is_swapik2_bl(5).or.is_swapn_ik(1,nv)
             is_swap_jk=is_swapjk2_bl(5).or.is_swapn_jk(1,nv)
             is_rev_n  = is_rev2_bl(5,2).or.is_revn_p2(1,nv)
             is_rev_p1 = is_rev2_bl(5,1).or.is_revn_n(1,nv)
             is_rev_p2 = is_rev2_bl(5,3).or.is_revn_p1(1,nv)
             ! fill ghost cells in edge
             call fill_ghost_cells3d(5,1-ngh,0,1,ngy, &
                  is_swap_ij,is_swap_ik,is_swap_jk, &
                  is_rev_p1,is_rev_n,is_rev_p2)
             ! free temporary array
             deallocate(xgc3n,ygc3n,zgc3n)
          endif
       endif

       ! edge kmin-imax
       ! --------------
       if ((bl(n)%BC(5)>0).and.(bl(n)%BC(2)>0)) then
          ! start from imax
          nv=bl(n)%BC(2)
          nv2=bl(nv)%BC(5)
          ! start from kmin
          nv=bl(n)%BC(5)
          nv1=bl(nv)%BC(2)
          if (nv1.ne.nv2) then
             !print 100,nv1,nv2
             ! construct artificial edge kmin-imax
             call cons_grid3d(3,1,2,.true.)
          else
             ! read grid in neighboring block nv
             call read_neighbor3d(trim(gridname)//trim(numchar(nv1))//'.x')
             ! determine swap/reverse for corner block
             is_swap_ij=is_swapij2_bl(5).or.is_swapn_ij(2,nv)
             is_swap_ik=is_swapik2_bl(5).or.is_swapn_ik(2,nv)
             is_swap_jk=is_swapjk2_bl(5).or.is_swapn_jk(2,nv)
             is_rev_n  = is_rev2_bl(5,2).or.is_revn_p2(2,nv)
             is_rev_p1 = is_rev2_bl(5,1).or.is_revn_n(2,nv)
             is_rev_p2 = is_rev2_bl(5,3).or.is_revn_p1(2,nv)
             ! fill ghost cells in edge
             call fill_ghost_cells3d(5,ngx+1,ngx+ngh,1,ngy, &
                  is_swap_ij,is_swap_ik,is_swap_jk, &
                  is_rev_p1,is_rev_n,is_rev_p2)
             ! free temporary array
             deallocate(xgc3n,ygc3n,zgc3n)
          endif
       endif

       ! edge kmax-imin
       ! --------------
       if ((bl(n)%BC(6)>0).and.(bl(n)%BC(1)>0)) then
          ! start from imin
          nv=bl(n)%BC(1)
          nv2=bl(nv)%BC(6)
          ! start from kmax
          nv=bl(n)%BC(6)
          nv1=bl(nv)%BC(1)
          if (nv1.ne.nv2) then
             !print 100,nv1,nv2
             ! construct artificial edge kmax-imin
             call cons_grid3d(3,2,1,.true.)
          else
             ! read grid in neighboring block nv
             call read_neighbor3d(trim(gridname)//trim(numchar(nv1))//'.x')
             ! determine swap/reverse for corner block
             is_swap_ij=is_swapij2_bl(6).or.is_swapn_ij(1,nv)
             is_swap_ik=is_swapik2_bl(6).or.is_swapn_ik(1,nv)
             is_swap_jk=is_swapjk2_bl(6).or.is_swapn_jk(1,nv)
             is_rev_n  = is_rev2_bl(6,2).or.is_revn_p2(1,nv)
             is_rev_p1 = is_rev2_bl(6,1).or.is_revn_n(1,nv)
             is_rev_p2 = is_rev2_bl(6,3).or.is_revn_p1(1,nv)
             ! fill ghost cells in edge
             call fill_ghost_cells3d(6,1-ngh,0,1,ngy, &
                  is_swap_ij,is_swap_ik,is_swap_jk, &
                  is_rev_p1,is_rev_n,is_rev_p2)
             ! free temporary array
             deallocate(xgc3n,ygc3n,zgc3n)
          endif
       endif

       ! edge kmax-imax
       ! --------------
       if ((bl(n)%BC(6)>0).and.(bl(n)%BC(2)>0)) then
          ! start from imax
          nv=bl(n)%BC(2)
          nv2=bl(nv)%BC(6)
          ! start from kmax
          nv=bl(n)%BC(6)
          nv1=bl(nv)%BC(2)
          if (nv1.ne.nv2) then
             !print 100,nv1,nv2
             ! construct artificial edge kmax-imax
             call cons_grid3d(3,2,2,.true.)
          else
             ! read grid in neighboring block nv
             call read_neighbor3d(trim(gridname)//trim(numchar(nv1))//'.x')
             ! determine swap/reverse for corner block
             is_swap_ij=is_swapij2_bl(6).or.is_swapn_ij(2,nv)
             is_swap_ik=is_swapik2_bl(6).or.is_swapn_ik(2,nv)
             is_swap_jk=is_swapjk2_bl(6).or.is_swapn_jk(2,nv)
             is_rev_n  = is_rev2_bl(6,2).or.is_revn_p2(2,nv)
             is_rev_p1 = is_rev2_bl(6,1).or.is_revn_n(2,nv)
             is_rev_p2 = is_rev2_bl(6,3).or.is_revn_p1(2,nv)
             ! fill ghost cells in edge          
             call fill_ghost_cells3d(6,ngx+1,ngx+ngh,1,ngy, &
                  is_swap_ij,is_swap_ik,is_swap_jk, &
                  is_rev_p1,is_rev_n,is_rev_p2)
             ! free temporary array
             deallocate(xgc3n,ygc3n,zgc3n)
          endif
       endif

       ! Fill ghost cells for block corners
       ! ==================================
101    format(1x,' neighbors are different (',i3,' and ',i3,') ~> create artificial corner')

       ! corner imin-jmin-kmin (1,1,1)
       ! ---------------------
       if ((bl(n)%BC(1)>0).and.(bl(n)%BC(3)>0).and.(bl(n)%BC(5)>0)) then
          !print *,'corner imin-jmin-kmin'
!!$          ! order kmin-jmin-imin
!!$          nv=bl(n)%BC(5)
!!$          nv_=bl(nv)%BC(3)
!!$          nv6=bl(nv_)%BC(1)
!!$          ! order kmin-imin-jmin
!!$          nv=bl(n)%BC(5)
!!$          nv_=bl(nv)%BC(1)
!!$          nv5=bl(nv_)%BC(3)
!!$          ! order jmin-kmin-imin
!!$          nv=bl(n)%BC(3)
!!$          nv_=bl(n)%BC(5)
!!$          nv4=bl(nv_)%BC(1)
!!$          ! order jmin-imin-kmin
!!$          nv=bl(n)%BC(3)
!!$          nv_=bl(n)%BC(1)
!!$          nv3=bl(nv_)%BC(5)
!!$          ! order imin-kmin-jmin
!!$          nv=bl(n)%BC(1)
!!$          nv_=bl(nv)%BC(5)
!!$          nv2=bl(nv_)%BC(3)
!!$          ! order imin-jmin-kmin
!!$          nv=bl(n)%BC(1)
!!$          nv_=bl(nv)%BC(3)
!!$          nv1=bl(nv_)%BC(5)
!!$          !print *,nv1,nv2,nv3,nv4,nv5,nv6
!!$          if ((nv1.ne.nv2).or.(nv1.ne.nv3).or.(nv1.ne.nv4).or.(nv1.ne.nv5).or.(nv1.ne.nv6)) then
             ! construct artificial corner imax-jmin-kmin
             call cons_grid3d(1,1,1,.false.)
          !endif
       endif

       ! corner imin-jmin-kmax (1,1,2)
       ! ---------------------
       if ((bl(n)%BC(1)>0).and.(bl(n)%BC(3)>0).and.(bl(n)%BC(6)>0)) then
          !print *,'corner imin-jmin-kmax'
!!$          ! order kmin-jmin-imin
!!$          nv=bl(n)%BC(6)
!!$          nv_=bl(nv)%BC(3)
!!$          nv6=bl(nv_)%BC(1)
!!$          ! order kmax-imin-jmin
!!$          nv=bl(n)%BC(6)
!!$          nv_=bl(nv)%BC(1)
!!$          nv5=bl(nv_)%BC(3)
!!$          ! order jmin-kmax-imin
!!$          nv=bl(n)%BC(3)
!!$          nv_=bl(n)%BC(6)
!!$          nv4=bl(nv_)%BC(1)
!!$          ! order jmin-imin-kmax
!!$          nv=bl(n)%BC(3)
!!$          nv_=bl(n)%BC(1)
!!$          nv3=bl(nv_)%BC(6)
!!$          ! order imin-kmax-jmin
!!$          nv=bl(n)%BC(1)
!!$          nv_=bl(nv)%BC(6)
!!$          nv2=bl(nv_)%BC(3)
!!$          ! order imin-jmin-kmax
!!$          nv=bl(n)%BC(1)
!!$          nv_=bl(nv)%BC(3)
!!$          nv1=bl(nv_)%BC(6)
!!$          !print *,nv1,nv2,nv3,nv4,nv5,nv6
!!$          if ((nv1.ne.nv2).or.(nv1.ne.nv3).or.(nv1.ne.nv4).or.(nv1.ne.nv5).or.(nv1.ne.nv6)) then
             ! construct artificial corner imax-jmin-kmin
             call cons_grid3d(1,1,2,.false.)
          !endif
       endif

       ! corner imin-jmax-kmin (1,2,1)
       ! ---------------------
       if ((bl(n)%BC(1)>0).and.(bl(n)%BC(4)>0).and.(bl(n)%BC(5)>0)) then
          !print *,'corner imin-jmax-kmin'
!!$          ! order kmin-jmax-imin
!!$          nv=bl(n)%BC(5)
!!$          nv_=bl(nv)%BC(4)
!!$          nv6=bl(nv_)%BC(1)
!!$          ! order kmin-imin-jmax
!!$          nv=bl(n)%BC(5)
!!$          nv_=bl(nv)%BC(1)
!!$          nv5=bl(nv_)%BC(4)
!!$          ! order jmax-kmin-imin
!!$          nv=bl(n)%BC(4)
!!$          nv_=bl(n)%BC(5)
!!$          nv4=bl(nv_)%BC(1)
!!$          ! order jmax-imin-kmin
!!$          nv=bl(n)%BC(4)
!!$          nv_=bl(n)%BC(1)
!!$          nv3=bl(nv_)%BC(5)
!!$          ! order imin-kmin-jmax
!!$          nv=bl(n)%BC(1)
!!$          nv_=bl(nv)%BC(5)
!!$          nv2=bl(nv_)%BC(4)
!!$          ! order imin-jmax-kmin
!!$          nv=bl(n)%BC(1)
!!$          nv_=bl(nv)%BC(4)
!!$          nv1=bl(nv_)%BC(5)
!!$          !print *,nv1,nv2,nv3,nv4,nv5,nv6
!!$          if ((nv1.ne.nv2).or.(nv1.ne.nv3).or.(nv1.ne.nv4).or.(nv1.ne.nv5).or.(nv1.ne.nv6)) then
             ! construct artificial corner imax-jmin-kmin
             call cons_grid3d(1,2,1,.false.)
          !endif
       endif

       ! corner imin-jmax-kmax (1,2,2)
       ! ---------------------
       if ((bl(n)%BC(1)>0).and.(bl(n)%BC(4)>0).and.(bl(n)%BC(6)>0)) then
          !print *,'corner imin-jmax-kmax'
!!$          ! order kmax-jmax-imin
!!$          nv=bl(n)%BC(6)
!!$          nv_=bl(nv)%BC(4)
!!$          nv6=bl(nv_)%BC(1)
!!$          ! order kmax-imin-jmax
!!$          nv=bl(n)%BC(6)
!!$          nv_=bl(nv)%BC(1)
!!$          nv5=bl(nv_)%BC(4)
!!$          ! order jmax-kmax-imin
!!$          nv=bl(n)%BC(4)
!!$          nv_=bl(n)%BC(6)
!!$          nv4=bl(nv_)%BC(1)
!!$          ! order jmax-imin-kmax
!!$          nv=bl(n)%BC(4)
!!$          nv_=bl(n)%BC(1)
!!$          nv3=bl(nv_)%BC(6)
!!$          ! order imin-kmax-jmax
!!$          nv=bl(n)%BC(1)
!!$          nv_=bl(nv)%BC(6)
!!$          nv2=bl(nv_)%BC(4)
!!$          ! order imin-jmax-kmax
!!$          nv=bl(n)%BC(1)
!!$          nv_=bl(nv)%BC(4)
!!$          nv1=bl(nv_)%BC(6)
!!$          !print *,nv1,nv2,nv3,nv4,nv5,nv6
!!$          if ((nv1.ne.nv2).or.(nv1.ne.nv3).or.(nv1.ne.nv4).or.(nv1.ne.nv5).or.(nv1.ne.nv6)) then
             ! construct artificial corner imax-jmin-kmin
             call cons_grid3d(1,2,2,.false.)
          !endif
       endif

       ! corner imax-jmin-kmin (2,1,1)
       ! ---------------------
       if ((bl(n)%BC(2)>0).and.(bl(n)%BC(3)>0).and.(bl(n)%BC(5)>0)) then
          !print *,'corner imax-jmin-kmin'
!!$          ! order kmin-jmin-imax
!!$          nv=bl(n)%BC(5)
!!$          nv_=bl(nv)%BC(3)
!!$          nv6=bl(nv_)%BC(2)
!!$          ! order kmin-imax-jmin
!!$          nv=bl(n)%BC(5)
!!$          nv_=bl(nv)%BC(2)
!!$          nv5=bl(nv_)%BC(3)
!!$          ! order jmin-kmin-imax
!!$          nv=bl(n)%BC(3)
!!$          nv_=bl(n)%BC(5)
!!$          nv4=bl(nv_)%BC(2)
!!$          ! order jmin-imax-kmin
!!$          nv=bl(n)%BC(3)
!!$          nv_=bl(n)%BC(2)
!!$          nv3=bl(nv_)%BC(5)
!!$          ! order imax-kmin-jmin
!!$          nv=bl(n)%BC(2)
!!$          nv_=bl(nv)%BC(5)
!!$          nv2=bl(nv_)%BC(3)
!!$          ! order imax-jmin-kmin
!!$          nv=bl(n)%BC(2)
!!$          print *,'nv',nv
!!$          nv_=bl(nv)%BC(3)
!!$          print *,'nv_',nv_
!!$          nv1=bl(nv_)%BC(5)
!!$          print *,'nv1',nv1
!!$          !print *,nv1,nv2,nv3,nv4,nv5,nv6
!!$          if ((nv1.ne.nv2).or.(nv1.ne.nv3).or.(nv1.ne.nv4).or.(nv1.ne.nv5).or.(nv1.ne.nv6)) then
             ! construct artificial corner imax-jmin-kmin
             call cons_grid3d(2,1,1,.false.)
          !endif
       endif

       ! corner imax-jmin-kmax (2,1,2)
       ! ---------------------
       if ((bl(n)%BC(2)>0).and.(bl(n)%BC(3)>0).and.(bl(n)%BC(6)>0)) then
          !print *,'corner imax-jmin-kmax'
!!$          ! order kmax-jmin-imax
!!$          nv=bl(n)%BC(6)
!!$          nv_=bl(nv)%BC(3)
!!$          nv6=bl(nv_)%BC(2)
!!$          ! order kmax-imax-jmin
!!$          nv=bl(n)%BC(6)
!!$          nv_=bl(nv)%BC(2)
!!$          nv5=bl(nv_)%BC(3)
!!$          ! order jmin-kmax-imax
!!$          nv=bl(n)%BC(3)
!!$          nv_=bl(n)%BC(6)
!!$          nv4=bl(nv_)%BC(2)
!!$          ! order jmin-imax-kmax
!!$          nv=bl(n)%BC(3)
!!$          nv_=bl(n)%BC(2)
!!$          nv3=bl(nv_)%BC(6)
!!$          ! order imax-kmax-jmin
!!$          nv=bl(n)%BC(2)
!!$          nv_=bl(nv)%BC(6)
!!$          nv2=bl(nv_)%BC(3)
!!$          ! order imax-jmin-kmax
!!$          nv=bl(n)%BC(2)
!!$          nv_=bl(nv)%BC(3)
!!$          nv1=bl(nv_)%BC(6)
!!$          !print *,nv1,nv2,nv3,nv4,nv5,nv6
!!$          if ((nv1.ne.nv2).or.(nv1.ne.nv3).or.(nv1.ne.nv4).or.(nv1.ne.nv5).or.(nv1.ne.nv6)) then
             ! construct artificial corner imax-jmin-kmin
             call cons_grid3d(2,1,2,.false.)
          !endif
       endif

       ! corner imax-jmax-kmin (2,2,1)
       ! ---------------------
       if ((bl(n)%BC(2)>0).and.(bl(n)%BC(4)>0).and.(bl(n)%BC(5)>0)) then
          !print *,'corner imax-jmax-kmin'
!!$          ! order kmin-jmax-imax
!!$          nv=bl(n)%BC(5)
!!$          nv_=bl(nv)%BC(4)
!!$          nv6=bl(nv_)%BC(2)
!!$          ! order kmin-imax-jmax
!!$          nv=bl(n)%BC(5)
!!$          nv_=bl(nv)%BC(2)
!!$          nv5=bl(nv_)%BC(4)
!!$          ! order jmax-kmin-imax
!!$          nv=bl(n)%BC(4)
!!$          nv_=bl(n)%BC(5)
!!$          nv4=bl(nv_)%BC(2)
!!$          ! order jmax-imax-kmin
!!$          nv=bl(n)%BC(4)
!!$          nv_=bl(n)%BC(2)
!!$          nv3=bl(nv_)%BC(5)
!!$          ! order imax-kmin-jmax
!!$          nv=bl(n)%BC(2)
!!$          nv_=bl(nv)%BC(5)
!!$          nv2=bl(nv_)%BC(4)
!!$          ! order imax-jmax-kmin
!!$          nv=bl(n)%BC(2)
!!$          nv_=bl(nv)%BC(4)
!!$          nv1=bl(nv_)%BC(5)
!!$          !print *,nv1,nv2,nv3,nv4,nv5,nv6
!!$          if ((nv1.ne.nv2).or.(nv1.ne.nv3).or.(nv1.ne.nv4).or.(nv1.ne.nv5).or.(nv1.ne.nv6)) then
             ! construct artificial corner imax-jmin-kmin
             call cons_grid3d(2,2,1,.false.)
          !endif
       endif

       ! corner imax-jmax-kmax (2,2,2)
       ! ---------------------
       if ((bl(n)%BC(2)>0).and.(bl(n)%BC(4)>0).and.(bl(n)%BC(6)>0)) then
!!$          !print *,'corner imax-jmax-kmax'
!!$          ! order kmax-jmax-imax
!!$          nv=bl(n)%BC(6)
!!$          nv_=bl(nv)%BC(4)
!!$          nv6=bl(nv_)%BC(2)
!!$          ! order kmax-imax-jmax
!!$          nv=bl(n)%BC(6)
!!$          nv_=bl(nv)%BC(2)
!!$          nv5=bl(nv_)%BC(4)
!!$          ! order jmax-kmax-imax
!!$          nv=bl(n)%BC(4)
!!$          nv_=bl(n)%BC(6)
!!$          nv4=bl(nv_)%BC(2)
!!$          ! order jmax-imax-kmax
!!$          nv=bl(n)%BC(4)
!!$          nv_=bl(n)%BC(2)
!!$          nv3=bl(nv_)%BC(6)
!!$          ! order imax-kmax-jmax
!!$          nv=bl(n)%BC(2)
!!$          nv_=bl(nv)%BC(6)
!!$          nv2=bl(nv_)%BC(4)
!!$          ! order imax-jmax-kmax
!!$          nv=bl(n)%BC(2)
!!$          nv_=bl(nv)%BC(4)
!!$          nv1=bl(nv_)%BC(6)
!!$          !print *,nv1,nv2,nv3,nv4,nv5,nv6
!!$          if ((nv1.ne.nv2).or.(nv1.ne.nv3).or.(nv1.ne.nv4).or.(nv1.ne.nv5).or.(nv1.ne.nv6)) then
             ! construct artificial corner imax-jmin-kmin
             call cons_grid3d(2,2,2,.false.)
          !endif
       endif

       ! Extended dimensions
       ! -------------------
       ni1=1
       ni2=ngx
       nj1=1
       nj2=ngy
       nk1=1
       nk2=ngz
       if (bl(n)%BC(1)>0) ni1=1-ngh
       if (bl(n)%BC(2)>0) ni2=ngx+ngh
       if (bl(n)%BC(3)>0) nj1=1-ngh
       if (bl(n)%BC(4)>0) nj2=ngy+ngh
       if (bl(n)%BC(5)>0) nk1=1-ngh
       if (bl(n)%BC(6)>0) nk2=ngz+ngh
       
       ! Apply angular periodicity
       ! -------------------------
       if (theta_period.ne.0.0_wp) then
          if (bl(n)%periodic(1)) then ! at imin
             print *,'angular periodicity imin for block ',n
             ! apply rotation of first rows of ghost cells
             allocate(y_p(1-ngh:0),z_p(1-ngh:0))
             do k=nk1,nk2
                do j=nj1,nj2
                   y_p=ygc3e(1-ngh:0,j,k)
                   z_p=zgc3e(1-ngh:0,j,k)
                   ygc3e(1-ngh:0,j,k)= costp*y_p-sintp(1)*z_p
                   zgc3e(1-ngh:0,j,k)= sintp(1)*y_p+costp*z_p
                enddo
             enddo
             deallocate(y_p,z_p)
          endif
          if (bl(n)%periodic(2)) then ! at imax
             print *,'angular periodicity imax for block ',n
             ! apply rotation of last rows of ghost cells
             allocate(y_p(ngx+1:ngx+ngh),z_p(ngx+1:ngx+ngh))
             do k=nk1,nk2
                do j=nj1,nj2
                   y_p=ygc3e(ngx+1:ngx+ngh,j,k)
                   z_p=zgc3e(ngx+1:ngx+ngh,j,k)
                   ygc3e(ngx+1:ngx+ngh,j,k)= costp*y_p-sintp(2)*z_p
                   zgc3e(ngx+1:ngx+ngh,j,k)= sintp(2)*y_p+costp*z_p
                enddo
             enddo
             deallocate(y_p,z_p)
          endif
          if (bl(n)%periodic(3)) then ! at jmin
             print *,'angular periodicity jmin for block ',n
             ! apply rotation of first rows of ghost cells
             allocate(y_p(1-ngh:0),z_p(1-ngh:0))
             do k=nk1,nk2
                do i=ni1,ni2
                   y_p=ygc3e(i,1-ngh:0,k)
                   z_p=zgc3e(i,1-ngh:0,k)
                   ygc3e(i,1-ngh:0,k)= costp*y_p-sintp(3)*z_p
                   zgc3e(i,1-ngh:0,k)= sintp(3)*y_p+costp*z_p
                enddo
             enddo
             deallocate(y_p,z_p)
          endif
          if (bl(n)%periodic(4)) then ! at jmax
             print *,'angular periodicity jmax for block ',n
             ! apply rotation of last rows of ghost cells
             allocate(y_p(ngy+1:ngy+ngh),z_p(ngy+1:ngy+ngh))
             do k=nk1,nk2
                do i=ni1,ni2
                   y_p=ygc3e(i,ngy+1:ngy+ngh,k)
                   z_p=zgc3e(i,ngy+1:ngy+ngh,k)
                   ygc3e(i,ngy+1:ngy+ngh,k)= costp*y_p-sintp(4)*z_p
                   zgc3e(i,ngy+1:ngy+ngh,k)= sintp(4)*y_p+costp*z_p
                enddo
             enddo
             deallocate(y_p,z_p)
          endif
          if (bl(n)%periodic(5)) then ! at kmin
             print *,'angular periodicity kmin for block ',n
             ! apply rotation of first rows of ghost cells
             allocate(y_p(1-ngh:0),z_p(1-ngh:0))
             do j=nj1,nj2
                do i=ni1,ni2
                   y_p=ygc3e(i,j,1-ngh:0)
                   z_p=zgc3e(i,j,1-ngh:0)
                   ygc3e(i,j,1-ngh:0)= costp*y_p-sintp(5)*z_p
                   zgc3e(i,j,1-ngh:0)= sintp(5)*y_p+costp*z_p
                enddo
             enddo
             deallocate(y_p,z_p)
          endif
          if (bl(n)%periodic(6)) then ! at kmax
             print *,'angular periodicity kmax for block ',n
             ! apply rotation of last rows of ghost cells
             allocate(y_p(ngz+1:ngz+ngh),z_p(ngz+1:ngz+ngh))
             do j=nj1,nj2
                do i=ni1,ni2
                   y_p=ygc3e(i,j,ngz+1:ngz+ngh)
                   z_p=zgc3e(i,j,ngz+1:ngz+ngh)
                   ygc3e(i,j,ngz+1:ngz+ngh)= costp*y_p-sintp(6)*z_p
                   zgc3e(i,j,ngz+1:ngz+ngh)= sintp(6)*y_p+costp*z_p
                enddo
             enddo
             deallocate(y_p,z_p)
          endif
       ! Apply translation periodicity
       ! ------------------------------
       else
          if (bl(n)%periodic(1)) then ! at imin
             if (verbose) print *,'translation periodicity imin for block ',n
             ! apply translation of first rows of ghost cells
             do k=nk1,nk2
                do j=nj1,nj2
                   xgc3e(1-ngh:0,j,k)=xgc3e(1-ngh:0,j,k)+Lxp
                   ygc3e(1-ngh:0,j,k)=ygc3e(1-ngh:0,j,k)+Lyp
                   zgc3e(1-ngh:0,j,k)=zgc3e(1-ngh:0,j,k)+Lzp
                enddo
             enddo
          endif
          if (bl(n)%periodic(2)) then ! at imax
             if (verbose) print *,'translation periodicity imax for block ',n
             ! apply translation of last rows of ghost cells
             do k=nk1,nk2
                do j=nj1,nj2
                   xgc3e(ngx+1:ngx+ngh,j,k)=xgc3e(ngx+1:ngx+ngh,j,k)+Lxp
                   ygc3e(ngx+1:ngx+ngh,j,k)=ygc3e(ngx+1:ngx+ngh,j,k)+Lyp
                   zgc3e(ngx+1:ngx+ngh,j,k)=zgc3e(ngx+1:ngx+ngh,j,k)+Lzp
                enddo
             enddo
          endif
          if (bl(n)%periodic(3)) then ! at jmin
             if (verbose) print *,'translation periodicity jmin for block ',n
             ! apply translation of first rows of ghost cells
             do k=nk1,nk2
                do i=ni1,ni2
                   xgc3e(i,1-ngh:0,k)=xgc3e(i,1-ngh:0,k)+Lxp
                   ygc3e(i,1-ngh:0,k)=ygc3e(i,1-ngh:0,k)+Lyp
                   zgc3e(i,1-ngh:0,k)=zgc3e(i,1-ngh:0,k)+Lzp
                enddo
             enddo
          endif
          if (bl(n)%periodic(4)) then ! at jmax
             if (verbose) print *,'translation periodicity jmax for block ',n
             ! apply translation of last rows of ghost cells
             do k=nk1,nk2
                do i=ni1,ni2
                   xgc3e(i,ngy+1:ngy+ngh,k)=xgc3e(i,ngy+1:ngy+ngh,k)+Lxp
                   ygc3e(i,ngy+1:ngy+ngh,k)=ygc3e(i,ngy+1:ngy+ngh,k)+Lyp
                   zgc3e(i,ngy+1:ngy+ngh,k)=zgc3e(i,ngy+1:ngy+ngh,k)+Lzp
                enddo
             enddo
          endif
          if (bl(n)%periodic(5)) then ! at kmin
             if (verbose) print *,'translation periodicity kmin for block ',n
             ! apply translation of first rows of ghost cells
             zmin=minval(minval(minval(zgc31(1:ngx,1:ngy,1:ngz),1),1),1)
             zmax=maxval(maxval(maxval(zgc31(1:ngx,1:ngy,1:ngz),1),1),1)      
             Lzp0=zmax-zmin+abs(zgc31(1,1,2)-zgc31(1,1,1))
             do j=nj1,nj2
                do i=ni1,ni2
                   !xgc3e(i,j,1-ngh:0)=xgc3e(i,j,1-ngh:0)+0.0_wp
                   !ygc3e(i,j,1-ngh:0)=ygc3e(i,j,1-ngh:0)+0.0_wp
                   zgc3e(i,j,1-ngh:0)=zgc3e(i,j,1-ngh:0)-Lzp0
                enddo
             enddo
          endif
          if (bl(n)%periodic(6)) then ! at kmax
             if (verbose) print *,'translation periodicity kmax for block ',n
             ! apply translation of last rows of ghost cells
             zmin=minval(minval(minval(zgc31(1:ngx,1:ngy,1:ngz),1),1),1)
             zmax=maxval(maxval(maxval(zgc31(1:ngx,1:ngy,1:ngz),1),1),1)      
             Lzp0=zmax-zmin+abs(zgc31(1,1,2)-zgc31(1,1,1))
             do j=nj1,nj2
                do i=ni1,ni2
                   !xgc3e(i,j,ngz+1:ngz+ngh)=xgc3e(i,j,ngz+1:ngz+ngh)+0.0_wp
                   !ygc3e(i,j,ngz+1:ngz+ngh)=ygc3e(i,j,ngz+1:ngz+ngh)+0.0_wp
                   zgc3e(i,j,ngz+1:ngz+ngh)=zgc3e(i,j,ngz+1:ngz+ngh)+Lzp0
                enddo
             enddo
          endif
       endif

       ! Reconstruct ghost cells of degenerate edges only for j directions
       ! -----------------------------------------------------------------
       !call correct_edges_i
       !if (theta_period==0.0_wp) call correct_edges_j
      
!!$       ! write grid after modifications for block n
!!$       ! -------------------------------------------
!!$       n=nob(iproc)
!!$       if (is_adjoint_block) then
!!$          gridname=trim(dirGRID)//'/'//trim(nameGRID)//'_ngh'//trim(numchar(ngh))//'_bl'
!!$       else
!!$          gridname=trim(dirGRID)//'/'//trim(nameGRID)//'_ngh'//trim(numchar(ngh))//'_mod_bl'
!!$       endif       
!!$       open(50,file=trim(gridname)//trim(numchar(n))//'.x',form='formatted')
!!$       write(50,*) 1
!!$       ni=ni2-ni1+1
!!$       nj=nj2-nj1+1
!!$       nk=nk2-nk1+1
!!$       write(50,*) ni,nj,nk
!!$       write(50,*) (((xgc3e(i,j,k),i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
!!$       write(50,*) (((ygc3e(i,j,k),i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
!!$       write(50,*) (((zgc3e(i,j,k),i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)

       ! free grid arrays
       ! ----------------
       ! deallocate(xgc31,ygc31,zgc31)
       deallocate(is_swapn_ij,is_swapn_ik,is_swapn_jk)
       deallocate(is_revn_n,is_revn_p1,is_revn_p2)

       xc3=xgc3e
       yc3=ygc3e
       zc3=zgc3e
       
        ! init MPI-IO for grid read/write tata
        call init_io_grid_ex3d

        ! gridname
        n=nob(iproc)
        if (verbose) print *,'block',n,iproc
        if (is_adjoint_block) then
           gridname=trim(dirGRID)//'/'//trim(nameGRID)// &
                '_ngh'//trim(numchar(ngh))//'_bl'//trim(numchar(n))//filext_write
        else
           gridname=trim(dirGRID)//'/'//trim(nameGRID)// &
                '_ngh'//trim(numchar(ngh))//'_mod_bl'//trim(numchar(n))//filext_write
        endif
        
        ! read extended grid
        call read_write_grid3d(gridname,WRITE)

        ! free memory
        call free_grid_ex3d
       deallocate(xgc3e,ygc3e,zgc3e)
    endif

    !call mpistop('test...', 0)

  end subroutine add_ghost_cells_grid3d
  
  !===============================================================================
  subroutine read_neighbor3d(gridname)
  !===============================================================================
    !> read grid in neighboring block nv
  !===============================================================================
    use mod_utils
    implicit none
    ! ----------------------------------------------------------------------------
    character(len=*), intent(in) :: gridname
    ! ----------------------------------------------------------------------------
    integer :: i,j,k,mb
    ! ----------------------------------------------------------------------------

    ! Open plot3D file
    ! ----------------
    open(50,file=gridname,form='formatted')
    rewind(50)
    read(50,*) mb
    
    ! Read dimensions
    ! ---------------
    read(50,*) ngxn,ngyn,ngzn

    ! Allocate grid arrays
    ! --------------------
    allocate(xgc3n(ngxn,ngyn,ngzn),ygc3n(ngxn,ngyn,ngzn),zgc3n(ngxn,ngyn,ngzn))

    ! Read neighbor grid
    ! ------------------
    read(50,*) (((xgc3n(i,j,k),i=1,ngxn),j=1,ngyn),k=1,ngzn)
    read(50,*) (((ygc3n(i,j,k),i=1,ngxn),j=1,ngyn),k=1,ngzn)
    read(50,*) (((zgc3n(i,j,k),i=1,ngxn),j=1,ngyn),k=1,ngzn)
    close(50)
    
  end subroutine read_neighbor3d

  !===============================================================================
  subroutine cons_grid3d(ind1,ind2,ind3,is_edge)
  !===============================================================================
    !> construct artificial grid for egdes or corners
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer,intent(in) :: ind1,ind2,ind3
    logical :: is_edge
    ! ----------------------------------------------------------------------------
    integer :: i,j,k,i1,j1,k1,ii,jj,kk
    integer :: ideb,ifin
    real(wp) :: x1,x2,x3,y1,y2,y3,z1,z2,z3
    ! ----------------------------------------------------------------------------

    ! Definition of edges
    ! ===================
    ! 1,1,1 : imin-jmin
    ! 1,1,2 : imin-jmax
    ! 1,2,1 : imax-jmin
    ! 1,2,2 : imax-jmax
    ! 2,1,1 : jmin-kmin
    ! 2,1,2 : jmin-kmax
    ! 2,2,1 : jmax-kmin
    ! 2,2,2 : jmax-kmax
    ! 3,1,1 : kmin-imin
    ! 3,1,2 : kmin-imax
    ! 3,2,1 : kmax-imin
    ! 3,2,2 : kmax-imax

    if (is_edge) then
       ideb=1
       if (ind1==1) ifin=ngz
       if (ind1==2) ifin=ngx
       if (ind1==3) ifin=ngy
    endif
    
    ! Definition of corners
    ! ===================
    ! 1,1,1 : imin-jmin-kmin
    ! 1,1,2 : imin-jmin-kmax
    ! 1,2,1 : imin-jmax-kmin
    ! 1,2,2 : imin-jmax-kmax
    ! 2,1,1 : imax-jmin-kmin
    ! 2,1,2 : imax-jmin-kmax
    ! 2,2,1 : imax-jmax-kmin
    ! 2,2,2 : imax-jmax-kmax

    if (.not.is_edge) then
       if (ind1==1) then
          if (ind3==1) then
             ideb=1-ngh
             ifin=0
          else
             ideb=ngz+1
             ifin=ngz+ngh
          endif
       else
          ideb=ngx+1
          ifin=ngx+ngh
       endif
    endif
        
    ! Edges ij / Corners imin--
    ! -------------------------
    if (ind1==1) then
       if (is_edge) then
          if (ind2==1) i1=1
          if (ind2==2) i1=ngx
          if (ind3==1) j1=1
          if (ind3==2) j1=ngy
       else
          if (ind1==1) i1=1
          if (ind1==2) i1=ngx
          if (ind2==1) j1=1
          if (ind2==2) j1=ngy
       endif
       ! construct edge or corner
       do k=ideb,ifin
                    ! origin point
          x1=xgc3e(i1,j1,k)
          y1=ygc3e(i1,j1,k)
          z1=zgc3e(i1,j1,k)

          do i=1,ngh
             do j=1,ngh
                if (is_edge) then
                   if (ind2==1) ii=i-ngh
                   if (ind2==2) ii=ngx+i
                   if (ind3==1) jj=j-ngh
                   if (ind3==2) jj=ngy+j
                else
                   if (ind1==1) ii=i-ngh
                   if (ind1==2) ii=ngx+i
                   if (ind2==1) jj=j-ngh
                   if (ind2==2) jj=ngy+j
                endif
                ! corner point #1
                x2=xgc3e(ii,j1,k)
                y2=ygc3e(ii,j1,k)
                z2=zgc3e(ii,j1,k)
                ! corner point #2
                x3=xgc3e(i1,jj,k)
                y3=ygc3e(i1,jj,k)
                z3=zgc3e(i1,jj,k)
                ! reconstruction
                xgc3e(ii,jj,k)=x2+x3-x1
                ygc3e(ii,jj,k)=y2+y3-y1
                zgc3e(ii,jj,k)=z2+z3-z1
             enddo
          enddo
       enddo
       
    ! Edges jk / Corners imax--
    ! -------------------------
    elseif (ind1==2) then
       if (ind2==1) j1=1
       if (ind2==2) j1=ngy
       if (ind3==1) k1=1
       if (ind3==2) k1=ngz
       ! construct edge or corner
       do i=ideb,ifin
          ! origin point
          x1=xgc3e(i,j1,k1)
          y1=ygc3e(i,j1,k1)
          z1=zgc3e(i,j1,k1)

          do j=1,ngh
             do k=1,ngh
                if (ind2==1) jj=j-ngh
                if (ind2==2) jj=ngy+j
                if (ind3==1) kk=k-ngh
                if (ind3==2) kk=ngz+k
                ! corner point #1
                x2=xgc3e(i,jj,k1)
                y2=ygc3e(i,jj,k1)
                z2=zgc3e(i,jj,k1)
                ! corner point #2
                x3=xgc3e(i,j1,kk)
                y3=ygc3e(i,j1,kk)
                z3=zgc3e(i,j1,kk)
                ! reconstruction
                xgc3e(i,jj,kk)=x2+x3-x1
                ygc3e(i,jj,kk)=y2+y3-y1
                zgc3e(i,jj,kk)=z2+z3-z1
                !if ((ind1==2).and.(ind2==2).and.(ind3==2).and.(.not.is_edge)) print *,i,jj,kk,xgc3e(i,jj,kk)
             enddo
          enddo
       enddo
       
    ! Edges ik
    ! --------
    elseif (ind1==3) then
       if (ind2==1) k1=1
       if (ind2==2) k1=ngz
       if (ind3==1) i1=1
       if (ind3==2) i1=ngx
       ! construct edge or corner
       do j=ideb,ifin
          ! origin point
          x1=xgc3e(i1,j,k1)
          y1=ygc3e(i1,j,k1)
          z1=zgc3e(i1,j,k1)

          do i=1,ngh
             do k=1,ngh
                if (ind2==1) kk=k-ngh
                if (ind2==2) kk=ngz+k
                if (ind3==1) ii=i-ngh
                if (ind3==2) ii=ngx+i
                ! corner point #1
                x2=xgc3e(ii,j,k1)
                y2=ygc3e(ii,j,k1)
                z2=zgc3e(ii,j,k1)
                ! corner point #2
                x3=xgc3e(i1,j,kk)
                y3=ygc3e(i1,j,kk)
                z3=zgc3e(i1,j,kk)
                ! reconstruction
                xgc3e(ii,j,kk)=x2+x3-x1
                ygc3e(ii,j,kk)=y2+y3-y1
                zgc3e(ii,j,kk)=z2+z3-z1
             enddo
          enddo
       enddo
    endif

  end subroutine cons_grid3d

!!$  !===============================================================================
!!$  subroutine cons_grid3d(ind1,ind2,ind3,is_edge)
!!$  !===============================================================================
!!$    !> construct artificial grid for egdes or corners
!!$  !===============================================================================
!!$    implicit none
!!$    ! ----------------------------------------------------------------------------
!!$    integer,intent(in) :: ind1,ind2,ind3
!!$    logical :: is_edge
!!$    ! ----------------------------------------------------------------------------
!!$    integer :: i,j,k,i1,j1,k1,ii,jj,kk
!!$    integer :: ideb,ifin
!!$    real(wp) :: x1,x2,x3,y1,y2,y3,z1,z2,z3
!!$    ! ----------------------------------------------------------------------------
!!$
!!$       print*,'bb',ngx,ngy,ngz,ngx1,ngy1,ngz1
!!$    ! Definition of edges
!!$    ! ===================
!!$    ! 1,1,1 : imin-jmin
!!$    ! 1,1,2 : imin-jmax
!!$    ! 1,2,1 : imax-jmin
!!$    ! 1,2,2 : imax-jmax
!!$    ! 2,1,1 : jmin-kmin
!!$    ! 2,1,2 : jmin-kmax
!!$    ! 2,2,1 : jmax-kmin
!!$    ! 2,2,2 : jmax-kmax
!!$    ! 3,1,1 : kmin-imin
!!$    ! 3,1,2 : kmin-imax
!!$    ! 3,2,1 : kmax-imin
!!$    ! 3,2,2 : kmax-imax
!!$
!!$    if (is_edge) then
!!$       ideb=1
!!$       if (ind1==1) ifin=ngz1
!!$       if (ind1==2) ifin=ngx1
!!$       if (ind1==3) ifin=ngy1
!!$    endif
!!$    
!!$    ! Definition of corners
!!$    ! ===================
!!$    ! 1,1,1 : imin-jmin-kmin
!!$    ! 1,1,2 : imin-jmin-kmax
!!$    ! 1,2,1 : imin-jmax-kmin
!!$    ! 1,2,2 : imin-jmax-kmax
!!$    ! 2,1,1 : imax-jmin-kmin
!!$    ! 2,1,2 : imax-jmin-kmax
!!$    ! 2,2,1 : imax-jmax-kmin
!!$    ! 2,2,2 : imax-jmax-kmax
!!$
!!$    if (.not.is_edge) then
!!$       if (ind1==1) then
!!$          if (ind3==1) then
!!$             ideb=-4
!!$             ifin=0
!!$          else
!!$             ideb=ngz1+1
!!$             ifin=ngz1+5
!!$          endif
!!$       else
!!$          ideb=ngx1+1
!!$          ifin=ngx1+5
!!$       endif
!!$    endif
!!$        
!!$    ! Edges ij / Corners imin--
!!$    ! -------------------------
!!$    if (ind1==1) then
!!$       if (is_edge) then
!!$          if (ind2==1) i1=1
!!$          if (ind2==2) i1=ngx1
!!$          if (ind3==1) j1=1
!!$          if (ind3==2) j1=ngy1
!!$       else
!!$          if (ind1==1) i1=1
!!$          if (ind1==2) i1=ngx1
!!$          if (ind2==1) j1=1
!!$          if (ind2==2) j1=ngy1
!!$       endif
!!$       ! construct edge or corner
!!$       do k=ideb,ifin
!!$                    ! origin point
!!$          x1=xgc3e(i1,j1,k)
!!$          y1=ygc3e(i1,j1,k)
!!$          z1=zgc3e(i1,j1,k)
!!$
!!$          do i=1,5
!!$             do j=1,5
!!$                if (is_edge) then
!!$                   if (ind2==1) ii=i-5
!!$                   if (ind2==2) ii=ngx1+i
!!$                   if (ind3==1) jj=j-5
!!$                   if (ind3==2) jj=ngy1+j
!!$                else
!!$                   if (ind1==1) ii=i-5
!!$                   if (ind1==2) ii=ngx1+i
!!$                   if (ind2==1) jj=j-5
!!$                   if (ind2==2) jj=ngy1+j
!!$                endif
!!$                ! corner point #1
!!$                x2=xgc3e(ii,j1,k)
!!$                y2=ygc3e(ii,j1,k)
!!$                z2=zgc3e(ii,j1,k)
!!$                ! corner point #2
!!$                x3=xgc3e(i1,jj,k)
!!$                y3=ygc3e(i1,jj,k)
!!$                z3=zgc3e(i1,jj,k)
!!$                ! reconstruction
!!$                xgc3e(ii,jj,k)=x2+x3-x1
!!$                ygc3e(ii,jj,k)=y2+y3-y1
!!$                zgc3e(ii,jj,k)=z2+z3-z1
!!$             enddo
!!$          enddo
!!$       enddo
!!$       
!!$    ! Edges jk / Corners imax--
!!$    ! -------------------------
!!$    elseif (ind1==2) then
!!$       if (ind2==1) j1=1
!!$       if (ind2==2) j1=ngy1
!!$       if (ind3==1) k1=1
!!$       if (ind3==2) k1=ngz1
!!$       ! construct edge or corner
!!$       do i=ideb,ifin
!!$          ! origin point
!!$          x1=xgc3e(i,j1,k1)
!!$          y1=ygc3e(i,j1,k1)
!!$          z1=zgc3e(i,j1,k1)
!!$
!!$          do j=1,5
!!$             do k=1,5
!!$                if (ind2==1) jj=j-5
!!$                if (ind2==2) jj=ngy1+j
!!$                if (ind3==1) kk=k-5
!!$                if (ind3==2) kk=ngz1+k
!!$                ! corner point #1
!!$                x2=xgc3e(i,jj,k1)
!!$                y2=ygc3e(i,jj,k1)
!!$                z2=zgc3e(i,jj,k1)
!!$                ! corner point #2
!!$                x3=xgc3e(i,j1,kk)
!!$                y3=ygc3e(i,j1,kk)
!!$                z3=zgc3e(i,j1,kk)
!!$                ! reconstruction
!!$                xgc3e(i,jj,kk)=x2+x3-x1
!!$                ygc3e(i,jj,kk)=y2+y3-y1
!!$                zgc3e(i,jj,kk)=z2+z3-z1
!!$                !if ((ind1==2).and.(ind2==2).and.(ind3==2).and.(.not.is_edge)) print *,i,jj,kk,xgc3e(i,jj,kk)
!!$             enddo
!!$          enddo
!!$       enddo
!!$       
!!$    ! Edges ik
!!$    ! --------
!!$    elseif (ind1==3) then
!!$       if (ind2==1) k1=1
!!$       if (ind2==2) k1=ngz1
!!$       if (ind3==1) i1=1
!!$       if (ind3==2) i1=ngx1
!!$       ! construct edge or corner
!!$       do j=ideb,ifin
!!$          ! origin point
!!$          x1=xgc3e(i1,j,k1)
!!$          y1=ygc3e(i1,j,k1)
!!$          z1=zgc3e(i1,j,k1)
!!$
!!$          do i=1,5
!!$             do k=1,5
!!$                if (ind2==1) kk=k-5
!!$                if (ind2==2) kk=ngz1+k
!!$                if (ind3==1) ii=i-5
!!$                if (ind3==2) ii=ngx1+i
!!$                ! corner point #1
!!$                x2=xgc3e(ii,j,k1)
!!$                y2=ygc3e(ii,j,k1)
!!$                z2=zgc3e(ii,j,k1)
!!$                ! corner point #2
!!$                x3=xgc3e(i1,j,kk)
!!$                y3=ygc3e(i1,j,kk)
!!$                z3=zgc3e(i1,j,kk)
!!$                ! reconstruction
!!$                xgc3e(ii,j,kk)=x2+x3-x1
!!$                ygc3e(ii,j,kk)=y2+y3-y1
!!$                zgc3e(ii,j,kk)=z2+z3-z1
!!$             enddo
!!$          enddo
!!$       enddo
!!$    endif
!!$
!!$  end subroutine cons_grid3d
  
  !===============================================================================
  subroutine fill_ghost_cells3d(n_f,n_1,n_2,n_3,n_4, &
       is_swap_ij,is_swap_ik,is_swap_jk,is_rev_p1,is_rev_n,is_rev_p2)
  !===============================================================================
    !> fill ghost cells from neighboring block
  !===============================================================================
    use mod_mpi
    implicit none
    ! ----------------------------------------------------------------------------
    integer, intent(in) :: n_f ! face number
    integer, intent(in) :: n_1,n_2,n_3,n_4 ! sizes of parallel directions
    logical, intent(in) :: is_swap_ij,is_swap_ik,is_swap_jk ! swap directions
    logical, intent(in) :: is_rev_p1,is_rev_n,is_rev_p2 ! reverse directions
    ! ----------------------------------------------------------------------------
    integer :: i,j,k,i2,j2,k2,lp1,lp2
    ! ----------------------------------------------------------------------------

    ! Indices n_1 & n_2 are used for parallel direction #1
    ! ====================================================
    ! for imin- imax-: n_1 n_2 is direction j
    ! for jmin- jmax-: n_1 n_2 is direction i
    ! for kmin- kmax-: n_1 n_2 is direction i

    ! Indices n_3 & n_4 are used for parallel direction #2
    ! ====================================================
    ! for imin- imax-: n_3 n_4 is direction k
    ! for jmin- jmax-: n_3 n_4 is direction k
    ! for kmin- kmax-: n_3 n_4 is direction j

    ! Number of points in parallel directions
    ! =======================================
    lp1=n_2-n_1+1
    lp2=n_4-n_3+1
   
    ! Fill ghost cells for block faces
    ! ================================

    ! face imin
    ! ---------
    select case (n_f)

    case (1)
    if (is_swap_ij) then
          do k=n_3,n_4
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
                   k2=k
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
                   ! reverse parallel 2
                   if (is_rev_p2) then
                      k2=ngzn-k+1
                   endif
                   xgc3e(i,j,k)=xgc3n(i2,j2,k2)
                   ygc3e(i,j,k)=ygc3n(i2,j2,k2)
                   zgc3e(i,j,k)=zgc3n(i2,j2,k2)
                enddo
             enddo
          enddo
       elseif (is_swap_ik) then
          do k=n_3,n_4
             do j=n_1,n_2
                do i=1-ngh,0
                   k2=ngzn+i
                   if (n_1<0) then
                      j2=ngyn+j
                   elseif (n_1>ngy) then
                      j2=j-ngy
                   else
                      j2=j
                   endif
                   i2=k
                   ! reverse normal
                   if (is_rev_n) k2=1-i
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
                   ! reverse parallel 2
                   if (is_rev_p2) then
                      i2=ngxn-k+1
                   endif
                   xgc3e(i,j,k)=xgc3n(i2,j2,k2)
                   ygc3e(i,j,k)=ygc3n(i2,j2,k2)
                   zgc3e(i,j,k)=zgc3n(i2,j2,k2)
                enddo
             enddo
          enddo
       elseif (is_swap_jk) then
          do k=n_3,n_4
             do j=n_1,n_2
                do i=1-ngh,0
                   i2=ngxn+i
                   if (n_1<0) then
                      k2=ngzn+j
                   elseif (n_1>ngy) then
                      k2=j-ngy
                   else
                      k2=j
                   endif
                   j2=k
                   ! reverse normal
                   if (is_rev_n) i2=1-i
                   ! reverse parallel 1
                   if (is_rev_p1) then
                      if (n_1<0) then
                         k2=1-j
                      elseif (n_1>ngy) then
                         !k2=2*ngzn-j+1
                         k2=ngzn-k2+1
                      else
                         k2=ngzn-j+1
                      endif
                   endif
                   ! reverse parallel 2
                   if (is_rev_p2) then
                      j2=ngyn-k+1
                   endif
                   xgc3e(i,j,k)=xgc3n(i2,j2,k2)
                   ygc3e(i,j,k)=ygc3n(i2,j2,k2)
                   zgc3e(i,j,k)=zgc3n(i2,j2,k2)
                enddo
             enddo
          enddo
       else
          do k=n_3,n_4
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
                   k2=k
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
                   ! reverse parallel 2
                   if (is_rev_p2) then
                      k2=ngzn-k+1
                   endif
                   xgc3e(i,j,k)=xgc3n(i2,j2,k2)
                   ygc3e(i,j,k)=ygc3n(i2,j2,k2)
                   zgc3e(i,j,k)=zgc3n(i2,j2,k2)
                enddo
             enddo
          enddo
       endif

    ! face imax
    ! ---------
    case (2)
       if (is_swap_ij) then
          do k=n_3,n_4
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
                   k2=k

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
                   ! reverse parallel 2
                   if (is_rev_p2) then
                      k2=ngzn-k+1
                   endif
                   
                   xgc3e(i,j,k)=xgc3n(i2,j2,k2)
                   ygc3e(i,j,k)=ygc3n(i2,j2,k2)
                   zgc3e(i,j,k)=zgc3n(i2,j2,k2)
                enddo
             enddo
          enddo
       elseif (is_swap_ik) then
          do k=n_3,n_4
             do j=n_1,n_2
                do i=ngx+1,ngx+ngh
                   k2=i-ngx
                   if (n_1<0) then
                      j2=ngyn+j
                   elseif (n_1>ngy) then
                      j2=j-ngy
                   else
                      j2=j
                   endif
                   i2=k
                   ! reverse normal
                   !if (is_rev_n) k2=2*ngzn-i+1
                   if (is_rev_n) k2=ngzn-k2+1
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
                   ! reverse parallel 2
                   if (is_rev_p2) then
                      i2=ngxn-k+1
                   endif
                   xgc3e(i,j,k)=xgc3n(i2,j2,k2)
                   ygc3e(i,j,k)=ygc3n(i2,j2,k2)
                   zgc3e(i,j,k)=zgc3n(i2,j2,k2)
                enddo
             enddo
          enddo
       elseif (is_swap_jk) then
          do k=n_3,n_4
             do j=n_1,n_2
                do i=ngx+1,ngx+ngh
                   i2=i-ngx
                   if (n_1<0) then
                      k2=ngzn+j
                   elseif (n_1>ngy) then
                      k2=j-ngy
                   else
                      k2=j
                   endif
                   j2=k
                   ! reverse normal
                   !if (is_rev_n) i2=2*ngxn-i+1
                   if (is_rev_n) i2=ngxn-i2+1
                   ! reverse parallel 1
                   if (is_rev_p1) then
                      if (n_1<0) then
                         k2=1-j
                      elseif (n_1>ngy) then
                         !k2=2*ngzn-j+1
                         k2=ngzn-k2+1
                      else
                         k2=ngzn-j+1
                      endif
                   endif
                   ! reverse parallel 2
                   if (is_rev_p2) then
                      j2=ngyn-k+1
                   endif
                   xgc3e(i,j,k)=xgc3n(i2,j2,k2)
                   ygc3e(i,j,k)=ygc3n(i2,j2,k2)
                   zgc3e(i,j,k)=zgc3n(i2,j2,k2)
                enddo
             enddo
          enddo
       else
          do k=n_3,n_4
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
                   k2=k
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
                   ! reverse parallel 2
                   if (is_rev_p2) then
                      k2=ngzn-k+1
                   endif
                   xgc3e(i,j,k)=xgc3n(i2,j2,k2)
                   ygc3e(i,j,k)=ygc3n(i2,j2,k2)
                   zgc3e(i,j,k)=zgc3n(i2,j2,k2)
                enddo
             enddo
          enddo
       endif

    ! face jmin
    ! ---------
    case (3)
       if (is_swap_ij) then
          do k=n_3,n_4
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
                   if (n_3<0) then
                      k2=ngzn+k
                   elseif (n_3>ngz) then
                      k2=k-ngz
                   else
                      k2=k
                   endif
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
                   ! reverse parallel 2
                   if (is_rev_p2) then
                      if (n_3<0) then
                         k2=1-k
                      elseif (n_3>ngz) then
                         !k2=2*ngzn-k+1
                         k2=ngzn-k2+1
                      else
                         k2=ngzn-k+1
                      endif
                   endif
                   xgc3e(i,j,k)=xgc3n(i2,j2,k2)
                   ygc3e(i,j,k)=ygc3n(i2,j2,k2)
                   zgc3e(i,j,k)=zgc3n(i2,j2,k2)
                enddo
             enddo
          enddo
       elseif (is_swap_ik) then
          do k=n_3,n_4
             do j=1-ngh,0
                do i=n_1,n_2
                   if (n_1<0) then
                      k2=ngzn+i
                   elseif (n_1>ngx) then
                      k2=i-ngx
                   else
                      k2=i
                   endif
                   j2=ngyn+j
                   if (n_3<0) then
                      i2=ngxn+k
                   elseif (n_3>ngz) then
                      i2=k-ngz
                   else
                      i2=k
                   endif
                   ! reverse normal
                   if (is_rev_n) j2=1-j
                   ! reverse parallel 1
                   if (is_rev_p1) then
                      if (n_1<0) then
                         k2=1-i
                      elseif (n_1>ngx) then
                         !k2=2*ngzn-i+1
                         k2=ngzn-k2+1
                      else
                         k2=ngzn-i+1
                      endif
                   endif
                   ! reverse parallel 2
                   if (is_rev_p2) then
                      if (n_3<0) then
                         i2=1-k
                      elseif (n_3>ngz) then
                         !i2=2*ngxn-k+1
                         i2=ngxn-i2+1
                      else
                         i2=ngxn-k+1
                      endif
                   endif
                   xgc3e(i,j,k)=xgc3n(i2,j2,k2)
                   ygc3e(i,j,k)=ygc3n(i2,j2,k2)
                   zgc3e(i,j,k)=zgc3n(i2,j2,k2)
                enddo
             enddo
          enddo
       elseif (is_swap_jk) then
          do k=n_3,n_4
             do j=1-ngh,0
                do i=n_1,n_2
                   if (n_1<0) then
                      i2=ngxn+i
                   elseif (n_1>ngx) then
                      i2=i-ngx
                   else
                      i2=i
                   endif
                   k2=ngzn+j
                   if (n_3<0) then
                      j2=ngyn+k
                   elseif (n_3>ngz) then
                      j2=k-ngz
                   else
                      j2=k
                   endif
                   ! reverse normal
                   if (is_rev_n) k2=1-j
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
                   ! reverse parallel 2
                   if (is_rev_p2) then
                      if (n_3<0) then
                         j2=1-k
                      elseif (n_3>ngz) then
                         !j2=2*ngyn-k+1
                         j2=ngyn-j2+1
                      else
                         j2=ngyn-k+1
                      endif
                   endif
                   xgc3e(i,j,k)=xgc3n(i2,j2,k2)
                   ygc3e(i,j,k)=ygc3n(i2,j2,k2)
                   zgc3e(i,j,k)=zgc3n(i2,j2,k2)
                enddo
             enddo
          enddo
       else
          !print *,'jmin',n_3,n_4,1-ngh,0,n_1,n_2,nob(iproc)
          do k=n_3,n_4
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
                   if (n_3<0) then
                      k2=ngzn+k
                   elseif (n_3>ngz) then
                      k2=k-ngz
                   else
                      k2=k
                   endif
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
                   ! reverse parallel 2
                   if (is_rev_p2) then
                      if (n_3<0) then
                         k2=1-k
                      elseif (n_3>ngz) then
                         !k2=2*ngzn-k+1
                         k2=ngzn-k2+1
                      else
                         k2=ngzn-k+1
                      endif
                   endif
                   xgc3e(i,j,k)=xgc3n(i2,j2,k2)
                   ygc3e(i,j,k)=ygc3n(i2,j2,k2)
                   zgc3e(i,j,k)=zgc3n(i2,j2,k2)
                enddo
             enddo
          enddo
       endif

    ! face jmax
    ! ---------
    case (4)
       if (is_swap_ij) then
          do k=n_3,n_4
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
                   if (n_3<0) then
                      k2=ngzn+k
                   elseif (n_3>ngz) then
                      k2=k-ngz
                   else
                      k2=k
                   endif
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
                   ! reverse parallel 2
                   if (is_rev_p2) then
                      if (n_3<0) then
                         k2=1-k
                      elseif (n_3>ngz) then
                         !k2=2*ngzn-k+1
                         k2=ngzn-k2+1
                      else
                         k2=ngzn-k+1
                      endif
                   endif
                   xgc3e(i,j,k)=xgc3n(i2,j2,k2)
                   ygc3e(i,j,k)=ygc3n(i2,j2,k2)
                   zgc3e(i,j,k)=zgc3n(i2,j2,k2)
                enddo
             enddo
          enddo
       elseif (is_swap_ik) then
          do k=n_3,n_4
             do j=ngy+1,ngy+ngh
                do i=n_1,n_2
                   if (n_1<0) then
                      k2=ngzn+i
                   elseif (n_1>ngx) then
                      k2=i-ngx
                   else
                      k2=i
                   endif
                   j2=j-ngy
                   if (n_3<0) then
                      i2=ngxn+k
                   elseif (n_3>ngz) then
                      i2=k-ngz
                   else
                      i2=k
                   endif
                   ! reverse normal
                   !if (is_rev_n) j2=2*ngyn-j+1
                   if (is_rev_n) j2=ngyn-j2+1
                   ! reverse parallel 1
                   if (is_rev_p1) then
                      if (n_1<0) then
                         k2=1-i
                      elseif (n_1>ngx) then
                         !k2=2*ngzn-i+1
                         k2=ngzn-k2+1
                      else
                         k2=ngzn-i+1
                      endif
                   endif
                   ! reverse parallel 2
                   if (is_rev_p2) then
                      if (n_3<0) then
                         i2=1-k
                      elseif (n_3>ngz) then
                         !i2=2*ngxn-k+1
                         i2=ngxn-i2+1
                      else
                         i2=ngxn-k+1
                      endif
                   endif
                   xgc3e(i,j,k)=xgc3n(i2,j2,k2)
                   ygc3e(i,j,k)=ygc3n(i2,j2,k2)
                   zgc3e(i,j,k)=zgc3n(i2,j2,k2)
                enddo
             enddo
          enddo
       elseif (is_swap_jk) then
          do k=n_3,n_4
             do j=ngy+1,ngy+ngh
                do i=n_1,n_2
                   if (n_1<0) then
                      i2=ngxn+i
                   elseif (n_1>ngx) then
                      i2=i-ngx
                   else
                      i2=i
                   endif
                   k2=j-ngy
                   if (n_3<0) then
                      j2=ngyn+k
                   elseif (n_3>ngz) then
                      j2=k-ngz
                   else
                      j2=k
                   endif
                   ! reverse normal
                   !if (is_rev_n) k2=2*ngzn-j+1
                   if (is_rev_n) k2=ngzn-k2+1
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
                   ! reverse parallel 2
                   if (is_rev_p2) then
                      if (n_3<0) then
                         j2=1-k
                      elseif (n_3>ngz) then
                         !j2=2*ngyn-k+1
                         j2=ngyn-j2+1
                      else
                         j2=ngyn-k+1
                      endif
                   endif
                   xgc3e(i,j,k)=xgc3n(i2,j2,k2)
                   ygc3e(i,j,k)=ygc3n(i2,j2,k2)
                   zgc3e(i,j,k)=zgc3n(i2,j2,k2)
                enddo
             enddo
          enddo
       else
          do k=n_3,n_4
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
                   if (n_3<0) then
                      k2=ngzn+k
                   elseif (n_3>ngz) then
                      k2=k-ngz
                   else
                      k2=k
                   endif
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
                   ! reverse parallel 2
                   if (is_rev_p2) then
                      if (n_3<0) then
                         k2=1-k
                      elseif (n_3>ngz) then
                         !k2=2*ngzn-k+1
                         k2=ngzn-k2+1
                      else
                         k2=ngzn-k+1
                      endif
                   endif
                   xgc3e(i,j,k)=xgc3n(i2,j2,k2)
                   ygc3e(i,j,k)=ygc3n(i2,j2,k2)
                   zgc3e(i,j,k)=zgc3n(i2,j2,k2)
                enddo
             enddo
          enddo
       endif

    ! face kmin
    ! ---------
    case (5)
       if (is_swap_ij) then
          do k=1-ngh,0
             do j=n_3,n_4
                do i=n_1,n_2
                   if (n_1<0) then
                      j2=ngyn+i
                   elseif (n_1>ngx) then
                      j2=i-ngx
                   else
                      j2=i
                   endif
                   i2=j
                   k2=ngzn+k
                   ! reverse normal
                   if (is_rev_n) k2=1-k
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
                   ! reverse parallel 2
                   if (is_rev_p2) then
                      i2=ngxn-j+1
                   endif
                   xgc3e(i,j,k)=xgc3n(i2,j2,k2)
                   ygc3e(i,j,k)=ygc3n(i2,j2,k2)
                   zgc3e(i,j,k)=zgc3n(i2,j2,k2)
                enddo
             enddo
          enddo
       elseif (is_swap_ik) then
          do k=1-ngh,0
             do j=n_3,n_4
                do i=n_1,n_2
                   if (n_1<0) then
                      k2=ngzn+i
                   elseif (n_1>ngx) then
                      k2=i-ngx
                   else
                      k2=i
                   endif
                   j2=j
                   i2=ngxn+k
                   ! reverse normal
                   if (is_rev_n) i2=1-k
                   ! reverse parallel 1
                   if (is_rev_p1) then
                      if (n_1<0) then
                         k2=1-i
                      elseif (n_1>ngx) then
                         !k2=2*ngzn-i+1
                         k2=ngzn-k2+1
                      else
                         k2=ngzn-i+1
                      endif
                   endif
                   ! reverse parallel 2
                   if (is_rev_p2) then
                      j2=ngyn-j+1
                   endif
                   xgc3e(i,j,k)=xgc3n(i2,j2,k2)
                   ygc3e(i,j,k)=ygc3n(i2,j2,k2)
                   zgc3e(i,j,k)=zgc3n(i2,j2,k2)
                enddo
             enddo
          enddo
       elseif (is_swap_jk) then
          do k=1-ngh,0
             do j=n_3,n_4
                do i=n_1,n_2
                   if (n_1<0) then
                      i2=ngxn+i
                   elseif (n_1>ngx) then
                      i2=i-ngx
                   else
                      i2=i
                   endif
                   k2=j
                   j2=ngyn+k
                   ! reverse normal
                   if (is_rev_n) j2=1-k
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
                   ! reverse parallel 2
                   if (is_rev_p2) then
                      k2=ngzn-j+1
                   endif
                   xgc3e(i,j,k)=xgc3n(i2,j2,k2)
                   ygc3e(i,j,k)=ygc3n(i2,j2,k2)
                   zgc3e(i,j,k)=zgc3n(i2,j2,k2)
                enddo
             enddo
          enddo
       else
          do k=1-ngh,0
             do j=n_3,n_4
                do i=n_1,n_2
                   if (n_1<0) then
                      i2=ngxn+i
                   elseif (n_1>ngx) then
                      i2=i-ngx
                   else
                      i2=i
                   endif
                   j2=j
                   k2=ngzn+k
                   ! reverse normal
                   if (is_rev_n) k2=1-k
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
                   ! reverse parallel 2
                   if (is_rev_p2) then
                      j2=ngyn-j+1
                   endif
                   xgc3e(i,j,k)=xgc3n(i2,j2,k2)
                   ygc3e(i,j,k)=ygc3n(i2,j2,k2)
                   zgc3e(i,j,k)=zgc3n(i2,j2,k2)
                enddo
             enddo
          enddo
       endif

    ! face kmax
    ! ---------
    case (6)
       if (is_swap_ij) then
          do k=ngz+1,ngz+ngh
             do j=n_3,n_4
                do i=n_1,n_2
                   if (n_1<0) then
                      j2=ngyn+i
                   elseif (n_1>ngx) then
                      j2=i-ngx
                   else
                      j2=i
                   endif
                   i2=j
                   k2=k-ngz
                   ! reverse normal
                   !if (is_rev_n) k2=2*ngzn-k+1
                   if (is_rev_n) k2=ngzn-k2+1
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
                   ! reverse parallel 2
                   if (is_rev_p2) then
                      i2=ngxn-j+1
                   endif
                   xgc3e(i,j,k)=xgc3n(i2,j2,k2)
                   ygc3e(i,j,k)=ygc3n(i2,j2,k2)
                   zgc3e(i,j,k)=zgc3n(i2,j2,k2)
                enddo
             enddo
          enddo
       elseif (is_swap_ik) then
          do k=ngz+1,ngz+ngh
             do j=n_3,n_4
                do i=n_1,n_2
                   if (n_1<0) then
                      k2=ngzn+i
                   elseif (n_1>ngx) then
                      k2=i-ngx
                   else
                      k2=i
                   endif
                   j2=j
                   i2=k-ngz
                   ! reverse normal
                   !if (is_rev_n) i2=2*ngxn-k+1
                   if (is_rev_n) i2=ngxn-i2+1
                   ! reverse parallel 1
                   if (is_rev_p1) then
                      if (n_1<0) then
                         k2=1-i
                      elseif (n_1>ngx) then
                         !k2=2*ngzn-i+1
                         k2=ngzn-k2+1
                      else
                         k2=ngzn-i+1
                      endif
                   endif
                   ! reverse parallel 2
                   if (is_rev_p2) then
                      j2=ngyn-j+1
                   endif
                   xgc3e(i,j,k)=xgc3n(i2,j2,k2)
                   ygc3e(i,j,k)=ygc3n(i2,j2,k2)
                   zgc3e(i,j,k)=zgc3n(i2,j2,k2)
                enddo
             enddo
          enddo
       elseif (is_swap_jk) then
          do k=ngz+1,ngz+ngh
             do j=n_3,n_4
                do i=n_1,n_2
                   if (n_1<0) then
                      i2=ngxn+i
                   elseif (n_1>ngx) then
                      i2=i-ngx
                   else
                      i2=i
                   endif
                   k2=j
                   j2=k-ngz
                   ! reverse normal
                   !if (is_rev_n) j2=2*ngyn-k+1
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
                   ! reverse parallel 2
                   if (is_rev_p2) then
                      k2=ngzn-j+1
                   endif
                   xgc3e(i,j,k)=xgc3n(i2,j2,k2)
                   ygc3e(i,j,k)=ygc3n(i2,j2,k2)
                   zgc3e(i,j,k)=zgc3n(i2,j2,k2)
                enddo
             enddo
          enddo
       else
          do k=ngz+1,ngz+ngh
             do j=n_3,n_4
                do i=n_1,n_2
                   if (n_1<0) then
                      i2=ngxn+i
                   elseif (n_1>ngx) then
                      i2=i-ngx
                   else
                      i2=i
                   endif
                   j2=j
                   k2=k-ngz
                   ! reverse normal
                   !if (is_rev_n) k2=2*ngzn-k+1
                   if (is_rev_n) k2=ngzn-k2+1
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
                   ! reverse parallel 2
                   if (is_rev_p2) then
                      j2=ngyn-j+1
                   endif
                   xgc3e(i,j,k)=xgc3n(i2,j2,k2)
                   ygc3e(i,j,k)=ygc3n(i2,j2,k2)
                   zgc3e(i,j,k)=zgc3n(i2,j2,k2)
                enddo
             enddo
          enddo
       endif

    end select
    
  end subroutine fill_ghost_cells3d

  !===============================================================================
  subroutine correct_edges_i
  !===============================================================================
    !> reconstruct ghost cells of degenerate edges for i direction
  !===============================================================================
    use mod_block
    use mod_mpi_part
    use mod_grid_directions
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k
    integer :: n,nv,nv_
    ! ----------------------------------------------------------------------------

    ! current block number
    ! -------------------
    n=nob(iproc)

    ! Reconstruct ghost cells of degenerate edges only for j directions
    ! ------------------------------------------------------------------
    
    !print *,'================================'
    !print *,'BC',bl(n)%BC

    ndzt=1
    nfzt=ngz
    if (bl(n)%BC(5)>0) ndzt=1-ngh
    if (bl(n)%BC(6)>0) nfzt=ngz+ngh

    ! imin-jmin
    ! ---------
    nv=bl(n)%BC(1) ! block imin
    if (is_swapij2_bl(1)) then
       if (is_rev2(1,1)) then
          nv_=bl(nv)%BC(2) ! block imax of imin
       else
          nv_=bl(nv)%BC(1) ! block imin of imin
       endif
    else
       nv_=bl(nv)%BC(3) ! block jmin of imin
    endif
    !print *,'imin-jmin',nv,nv_,bl(n)%BC(3),bl(n)%BC(4)
    !print *,'proc',coord(1),coord(2)
    if ((nv_>0).and.((nv_==bl(n)%BC(3)).or.(nv_==bl(n)%BC(4))).and.((coord(1)==0).and.(coord(2)==0))) then
       !print *,'do it'
       do k=ndzt,nfzt
          do j=1-ngh,0
             do i=1-ngh,0
                xgc3e(i,j,k)=xgc3e(j,1-i,k)
                ygc3e(i,j,k)=ygc3e(j,1-i,k)
                zgc3e(i,j,k)=zgc3e(j,1-i,k)
             enddo
          enddo
       enddo
    endif

    ! imin-jmax
    ! ---------
    nv=bl(n)%BC(1) ! block imin
    if (is_swapij2_bl(1)) then
       if (is_rev2(1,1)) then
          nv_=bl(nv)%BC(1) ! block imin of imin
       else
          nv_=bl(nv)%BC(2) ! block imax of imin
       endif
    else
       nv_=bl(nv)%BC(4) ! block jmax of imin
    endif
    !print *,'imin-jmax',nv,nv_,bl(n)%BC(3),bl(n)%BC(4)
    !print *,'proc',coord(1),coord(2)
    if ((nv_>0).and.((nv_==bl(n)%BC(3)).or.(nv_==bl(n)%BC(4))).and.((coord(1)==0).and.(coord(2)==ndomy-1))) then
       !print *,'do it'
       do k=ndzt,nfzt
          do j=1,5
             do i=1-ngh,0
                xgc3e(i,ngy+j,k)=xgc3e(1-j,ngy+i,k)
                ygc3e(i,ngy+j,k)=ygc3e(1-j,ngy+i,k)
                zgc3e(i,ngy+j,k)=zgc3e(1-j,ngy+i,k)
             enddo
          enddo
       enddo
    endif

    ! imax-jmin
    ! ---------
    nv=bl(n)%BC(2) ! block imax
    if (is_swapij2_bl(2)) then
       if (is_rev2(2,1)) then
          nv_=bl(nv)%BC(2) ! block imax of imax
       else
          nv_=bl(nv)%BC(1) ! block imin of imax
       endif
    else
       nv_=bl(nv)%BC(3) ! block jmin of imax
    endif
    !print *,'imax-jmin',nv,nv_,bl(n)%BC(3),bl(n)%BC(4)
    !print *,'proc',coord(1),coord(2)
    if ((nv_>0).and.((nv_==bl(n)%BC(3)).or.(nv_==bl(n)%BC(4))).and.((coord(1)==ndomx-1).and.(coord(2)==0))) then
       !print *,'do it'
       do k=ndzt,nfzt
          do j=1-ngh,0
             do i=1,5
                xgc3e(ngx+i,j,k)=xgc3e(ngx-j+1,i,k)
                ygc3e(ngx+i,j,k)=ygc3e(ngx-j+1,i,k)
                zgc3e(ngx+i,j,k)=zgc3e(ngx-j+1,i,k)
             enddo
          enddo
       enddo
    endif

    ! imax-jmax
    ! ---------
    nv=bl(n)%BC(2) ! block imax
    if (is_swapij2_bl(2)) then
       if (is_rev2(2,1)) then
          nv_=bl(nv)%BC(1) ! block imin of imax
       else
          nv_=bl(nv)%BC(2) ! block imax of imax
       endif
    else
       nv_=bl(nv)%BC(4) ! block jmax of imax
    endif
    !print *,'imax-jmax',nv,nv_,bl(n)%BC(3),bl(n)%BC(4)
    !print *,'proc',coord(1),coord(2)
    if ((nv_>0).and.((nv_==bl(n)%BC(3)).or.(nv_==bl(n)%BC(4))).and.((coord(1)==ndomx-1).and.(coord(2)==ndomy-1))) then
       !print *,'do it'
       do k=ndzt,nfzt
          do j=1,5
             do i=1,5
                xgc3e(ngx+i,ngy+j,k)=xgc3e(ngx+j,ngy-i+1,k)
                ygc3e(ngx+i,ngy+j,k)=ygc3e(ngx+j,ngy-i+1,k)
                zgc3e(ngx+i,ngy+j,k)=zgc3e(ngx+j,ngy-i+1,k)
             enddo
          enddo
       enddo
    endif

    ndxt=1
    nfxt=ngx
    if (bl(n)%BC(1)>0) ndxt=1-ngh
    if (bl(n)%BC(2)>0) nfxt=ngx+ngh

    ! jmin-kmin
    ! ---------
    nv=bl(n)%BC(3) ! block jmin
    if (is_swapjk2_bl(3)) then
       if (is_rev2(3,3)) then
          nv_=bl(nv)%BC(4) ! block jmax of jmin
       else
          nv_=bl(nv)%BC(3) ! block jmin of jmin
       endif
    else
       nv_=bl(nv)%BC(5) ! block kmin of jmin
    endif
    !print *,'jmin-kmin',nv,nv_,bl(n)%BC(5),bl(n)%BC(6)
    !print *,'proc',coord(2),coord(3)
    if ((nv_>0).and.((nv_==bl(n)%BC(5)).or.(nv_==bl(n)%BC(6))).and.((coord(2)==0).and.(coord(3)==0))) then
       !print *,'do it'
       do i=ndxt,nfxt
          do k=1-ngh,0
             do j=1-ngh,0
                xgc3e(i,j,k)=xgc3e(i,k,1-j)
                ygc3e(i,j,k)=ygc3e(i,k,1-j)
                zgc3e(i,j,k)=zgc3e(i,k,1-j)
             enddo
          enddo
       enddo
    endif

    ! jmin-kmax
    ! ---------
    nv=bl(n)%BC(3) ! block jmin
    if (is_swapjk2_bl(3)) then
       if (is_rev2(3,3)) then
          nv_=bl(nv)%BC(3) ! block jmin of jmin
       else
          nv_=bl(nv)%BC(4) ! block jmax of jmin
       endif
    else
       nv_=bl(nv)%BC(6) ! block kmax of jmin
    endif
    !print *,'jmin-kmax',nv,nv_,bl(n)%BC(5),bl(n)%BC(6)
    !print *,'proc',coord(2),coord(3)
    if ((nv_>0).and.((nv_==bl(n)%BC(5)).or.(nv_==bl(n)%BC(6))).and.((coord(2)==0).and.(coord(3)==ndomz-1))) then
       !print *,'do it'
       do i=ndxt,nfxt
          do k=1,5
             do j=1-ngh,0
                xgc3e(i,j,ngz+k)=xgc3e(i,1-k,ngz+j)
                ygc3e(i,j,ngz+k)=ygc3e(i,1-k,ngz+j)
                zgc3e(i,j,ngz+k)=zgc3e(i,1-k,ngz+j)
             enddo
          enddo
       enddo
    endif

    ! jmax-kmin
    ! ---------
    nv=bl(n)%BC(4) ! block jmax
    if (is_swapjk2_bl(4)) then
       if (is_rev2(4,3)) then
          nv_=bl(nv)%BC(4) ! block jmax of jmax
       else
          nv_=bl(nv)%BC(3) ! block jmin of jmax
       endif
    else
       nv_=bl(nv)%BC(5) ! block kmin of jmax
    endif
    !print *,'jmax-kmin',nv,nv_,bl(n)%BC(5),bl(n)%BC(6)
    !print *,'proc',coord(2),coord(3)
    if ((nv_>0).and.((nv_==bl(n)%BC(5)).or.(nv_==bl(n)%BC(6))).and.((coord(2)==ndomy-1).and.(coord(3)==0))) then
       !print *,'do it'
       do i=ndxt,nfxt
          do k=1-ngh,0
             do j=1,5
                xgc3e(i,ngy+j,k)=xgc3e(i,ngy-k+1,j)
                ygc3e(i,ngy+j,k)=ygc3e(i,ngy-k+1,j)
                zgc3e(i,ngy+j,k)=zgc3e(i,ngy-k+1,j)
             enddo
          enddo
       enddo
    endif

    ! jmax-kmax
    ! ---------
    nv=bl(n)%BC(4) ! block jmax
    if (is_swapjk2_bl(4)) then
       if (is_rev2(4,3)) then
          nv_=bl(nv)%BC(3) ! block jmin of jmax
       else
          nv_=bl(nv)%BC(4) ! block jmax of jmax
       endif
    else
       nv_=bl(nv)%BC(6) ! block kmax of jmax
    endif
    !print *,'jmax-kmax',nv,nv_,bl(n)%BC(5),bl(n)%BC(6)
    !print *,'proc',coord(2),coord(3)
    if ((nv_>0).and.((nv_==bl(n)%BC(5)).or.(nv_==bl(n)%BC(6))).and.((coord(2)==ndomy-1).and.(coord(3)==ndomz-1))) then
       !print *,'do it'
       do i=ndxt,nfxt
          do k=1,5
             do j=1,5
                xgc3e(i,ngy+j,ngz+k)=xgc3e(i,ngy+k,ngz-j+1)
                ygc3e(i,ngy+j,ngz+k)=ygc3e(i,ngy+k,ngz-j+1)
                zgc3e(i,ngy+j,ngz+k)=zgc3e(i,ngy+k,ngz-j+1)
             enddo
          enddo
       enddo
    endif

    ndyt=1
    nfyt=ngy
    if (bl(n)%BC(3)>0) ndyt=1-ngh
    if (bl(n)%BC(4)>0) nfyt=ngy+ngh

    ! imin-kmin
    ! ---------
    nv=bl(n)%BC(1) ! block imin
    if (is_swapik2_bl(1)) then
       if (is_rev2(1,3)) then
          nv_=bl(nv)%BC(2) ! block imax of imin
       else
          nv_=bl(nv)%BC(1) ! block imin of imin
       endif
    else
       nv_=bl(nv)%BC(5) ! block kmin of imin
    endif
    !print *,'imin-kmin',nv,nv_,bl(n)%BC(5),bl(n)%BC(6)
    !print *,'proc',coord(1),coord(3)
    if ((nv_>0).and.((nv_==bl(n)%BC(5)).or.(nv_==bl(n)%BC(6))).and.((coord(1)==0).and.(coord(3)==0))) then
       !print *,'do it'
       do j=ndyt,nfyt
          do k=1-ngh,0
             do i=1-ngh,0
                xgc3e(i,j,k)=xgc3e(k,j,1-i)
                ygc3e(i,j,k)=ygc3e(k,j,1-i)
                zgc3e(i,j,k)=zgc3e(k,j,1-i)
             enddo
          enddo
       enddo
    endif

    ! imin-kmax
    ! ---------
    nv=bl(n)%BC(1) ! block imin
    if (is_swapik2_bl(1)) then
       if (is_rev2(1,3)) then
          nv_=bl(nv)%BC(1) ! block imin of imin
       else
          nv_=bl(nv)%BC(2) ! block imax of imin
       endif
    else
       nv_=bl(nv)%BC(6) ! block kmax of imin
    endif
    !print *,'imin-kmax',nv,nv_,bl(n)%BC(5),bl(n)%BC(6)
    !print *,'proc',coord(1),coord(3)
    if ((nv_>0).and.((nv_==bl(n)%BC(5)).or.(nv_==bl(n)%BC(6))).and.((coord(1)==0).and.(coord(3)==ndomz-1))) then
       !print *,'do it'
       do j=ndyt,nfyt
          do k=1,5
             do i=1-ngh,0
                xgc3e(i,j,ngz+k)=xgc3e(1-k,j,ngz+i)
                ygc3e(i,j,ngz+k)=ygc3e(1-k,j,ngz+i)
                zgc3e(i,j,ngz+k)=zgc3e(1-k,j,ngz+i)
             enddo
          enddo
       enddo
    endif

    ! imax-kmin
    ! ---------
    nv=bl(n)%BC(2) ! block imax
    if (is_swapik2_bl(2)) then
       if (is_rev2(2,3)) then
          nv_=bl(nv)%BC(2) ! block imax of imax
       else
          nv_=bl(nv)%BC(1) ! block imin of imax
       endif
    else
       nv_=bl(nv)%BC(5) ! block kmin of imax
    endif
    !print *,'imax-kmin',nv,nv_,bl(n)%BC(5),bl(n)%BC(6)
    !print *,'proc',coord(1),coord(3)
    if ((nv_>0).and.((nv_==bl(n)%BC(5)).or.(nv_==bl(n)%BC(6))).and.((coord(1)==ndomx-1).and.(coord(3)==0))) then
       !print *,'do it'
       do j=ndyt,nfyt
          do k=1-ngh,0
             do i=1,5
                xgc3e(ngx+i,j,k)=xgc3e(ngx-k+1,j,i)
                ygc3e(ngx+i,j,k)=ygc3e(ngx-k+1,j,i)
                zgc3e(ngx+i,j,k)=zgc3e(ngx-k+1,j,i)
             enddo
          enddo
       enddo
    endif

    ! imax-kmax
    ! ---------
    nv=bl(n)%BC(2) ! block imax
    if (is_swapik2_bl(2)) then
       if (is_rev2(2,3)) then
          nv_=bl(nv)%BC(1) ! block imin of imax
       else
          nv_=bl(nv)%BC(2) ! block imax of imax
       endif      
    else
       nv_=bl(nv)%BC(6) ! block kmax of imax
    endif
    !print *,'jmax-kmax',nv,nv_,bl(n)%BC(5),bl(n)%BC(6)
    !print *,'proc',coord(1),coord(3)
    if ((nv_>0).and.((nv_==bl(n)%BC(5)).or.(nv_==bl(n)%BC(6))).and.((coord(1)==ndomx-1).and.(coord(3)==ndomz-1))) then
       !print *,'do it'
       do j=ndyt,nfyt
          do k=1,5
             do i=1,5
                xgc3e(ngx+i,j,ngz+k)=xgc3e(ngx+k,j,ngz-i+1)
                ygc3e(ngx+i,j,ngz+k)=ygc3e(ngx+k,j,ngz-i+1)
                zgc3e(ngx+i,j,ngz+k)=zgc3e(ngx+k,j,ngz-i+1)
             enddo
          enddo
       enddo
    endif

  end subroutine correct_edges_i
    
  !===============================================================================
  subroutine correct_edges_j
  !===============================================================================
    !> reconstruct ghost cells of degenerate edges for j direction
  !===============================================================================
    use mod_block
    use mod_mpi_part
    use mod_grid_directions
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k
    integer :: n,nv,nv_
    ! ----------------------------------------------------------------------------

    ! current block number
    ! -------------------
    n=nob(iproc)

    ! Reconstruct ghost cells of degenerate edges only for j directions
    ! ------------------------------------------------------------------
    
    !print *,'================================'
    !print *,'BC',bl(n)%BC

    ndzt=1
    nfzt=ngz
    if (bl(n)%BC(5)>0) ndzt=1-ngh
    if (bl(n)%BC(6)>0) nfzt=ngz+ngh

    ! imin-jmin
    ! ---------
    nv=bl(n)%BC(3) ! block jmin
    if (is_swapij2_bl(3)) then
       if (is_rev2(3,1)) then
          nv_=bl(nv)%BC(4) ! block jmax of jmin
       else
          nv_=bl(nv)%BC(3) ! block jmin of jmin
       endif
    else
       nv_=bl(nv)%BC(1) ! block imin of jmin
    endif
    !print *,'imin-jmin',nv,nv_,bl(n)%BC(1),bl(n)%BC(2)
    !print *,'proc',coord(1),coord(2)
    if ((nv_>0).and.((nv_==bl(n)%BC(1)).or.(nv_==bl(n)%BC(2))).and.((coord(1)==0).and.(coord(2)==0))) then
       !print *,'do it'
       do k=ndzt,nfzt
          do j=1-ngh,0
             do i=1-ngh,0
                xgc3e(i,j,k)=xgc3e(1-j,i,k)
                ygc3e(i,j,k)=ygc3e(1-j,i,k)
                zgc3e(i,j,k)=zgc3e(1-j,i,k)
             enddo
          enddo
       enddo
    endif

    ! imin-jmax
    ! ---------
    nv=bl(n)%BC(4) ! block jmax
    if (is_swapij2_bl(4)) then
       if (is_rev2(4,1)) then
          nv_=bl(nv)%BC(4) ! block jmax of jmax
       else
          nv_=bl(nv)%BC(3) ! block jmin of jmax
       endif
    else
       nv_=bl(nv)%BC(1) ! block imin of jmax
    endif
    !print *,'imin-jmax',nv,nv_,bl(n)%BC(1),bl(n)%BC(2)
    !print *,'proc',coord(1),coord(2)
    if ((nv_>0).and.((nv_==bl(n)%BC(1)).or.(nv_==bl(n)%BC(2))).and.((coord(1)==0).and.(coord(2)==ndomy-1))) then
       !print *,'do it'
       do k=ndzt,nfzt
          do j=1,5
             do i=1-ngh,0
                xgc3e(i,ngy+j,k)=xgc3e(j,ngy-i+1,k)
                ygc3e(i,ngy+j,k)=ygc3e(j,ngy-i+1,k)
                zgc3e(i,ngy+j,k)=zgc3e(j,ngy-i+1,k)
             enddo
          enddo
       enddo
    endif

    ! imax-jmin
    ! ---------
    nv=bl(n)%BC(3) ! block jmin
    if (is_swapij2_bl(3)) then
       if (is_rev2(3,1)) then
          nv_=bl(nv)%BC(3) ! block jmin of jmin
       else
          nv_=bl(nv)%BC(4) ! block jmax of jmin
       endif
    else
       nv_=bl(nv)%BC(2) ! block imax of jmin
    endif
    !print *,'imax-jmin',nv,nv_,bl(n)%BC(1),bl(n)%BC(2)
    !print *,'proc',coord(1),coord(2)
    if ((nv_>0).and.((nv_==bl(n)%BC(1)).or.(nv_==bl(n)%BC(2))).and.((coord(1)==ndomx-1).and.(coord(2)==0))) then
       !print *,'do it'
       do k=ndzt,nfzt
          do j=1-ngh,0
             do i=1,5
                xgc3e(ngx+i,j,k)=xgc3e(ngx+j,1-i,k)
                ygc3e(ngx+i,j,k)=ygc3e(ngx+j,1-i,k)
                zgc3e(ngx+i,j,k)=zgc3e(ngx+j,1-i,k)
             enddo
          enddo
       enddo
    endif

    ! imax-jmax
    ! ---------
    nv=bl(n)%BC(4) ! block jmax
    if (is_swapij2_bl(4)) then
       if (is_rev2(4,1)) then
          nv_=bl(nv)%BC(3) ! block jmin of jmax
       else
          nv_=bl(nv)%BC(4) ! block jmax of jmax
       endif
    else
       nv_=bl(nv)%BC(2) ! block imax of jmax
    endif
    !print *,'imax-jmax',nv,nv_,bl(n)%BC(1),bl(n)%BC(2)
    !print *,'proc',coord(1),coord(2)
    if ((nv_>0).and.((nv_==bl(n)%BC(1)).or.(nv_==bl(n)%BC(2))).and.((coord(1)==ndomx-1).and.(coord(2)==ndomy-1))) then
       !print *,'do it'
       do k=ndzt,nfzt
          do j=1,5
             do i=1,5
                xgc3e(ngx+i,ngy+j,k)=xgc3e(ngx-j+1,ngy+i,k)
                ygc3e(ngx+i,ngy+j,k)=ygc3e(ngx-j+1,ngy+i,k)
                zgc3e(ngx+i,ngy+j,k)=zgc3e(ngx-j+1,ngy+i,k)
             enddo
          enddo
       enddo
    endif

    ndxt=1
    nfxt=ngx
    if (bl(n)%BC(1)>0) ndxt=1-ngh
    if (bl(n)%BC(2)>0) nfxt=ngx+ngh

    ! jmin-kmin
    ! ---------
    nv=bl(n)%BC(5) ! block kmin
    if (is_swapjk2_bl(5)) then
       if (is_rev2(5,3)) then
          nv_=bl(nv)%BC(6) ! block kmax of kmin
       else
          nv_=bl(nv)%BC(5) ! block kmin of kmin
       endif
    else
       nv_=bl(nv)%BC(3) ! block jmin of kmin
    endif
    !print *,'jmin-kmin',nv,nv_,bl(n)%BC(3),bl(n)%BC(4)
    !print *,'proc',coord(2),coord(3)
    if ((nv_>0).and.((nv_==bl(n)%BC(3)).or.(nv_==bl(n)%BC(4))).and.((coord(2)==0).and.(coord(3)==0))) then
       !print *,'do it'
       do i=ndxt,nfxt
          do k=1-ngh,0
             do j=1-ngh,0
                xgc3e(i,j,k)=xgc3e(i,1-k,j)
                ygc3e(i,j,k)=ygc3e(i,1-k,j)
                zgc3e(i,j,k)=zgc3e(i,1-k,j)
             enddo
          enddo
       enddo
    endif

    ! jmin-kmax
    ! ---------
    nv=bl(n)%BC(6) ! block kmax
    if (is_swapjk2_bl(6)) then
       if (is_rev2(6,3)) then
          nv_=bl(nv)%BC(6) ! block kmax of kmax
       else
          nv_=bl(nv)%BC(5) ! block kmin of kmax
       endif
    else
       nv_=bl(nv)%BC(3) ! block jmin of kmax
    endif
    !print *,'jmin-kmax',nv,nv_,bl(n)%BC(3),bl(n)%BC(4)
    !print *,'proc',coord(2),coord(3)
    if ((nv_>0).and.((nv_==bl(n)%BC(3)).or.(nv_==bl(n)%BC(4))).and.((coord(2)==0).and.(coord(3)==ndomz-1))) then
       !print *,'do it'
       do i=ndxt,nfxt
          do k=1,5
             do j=1-ngh,0
                xgc3e(i,j,ngz+k)=xgc3e(i,k,ngz-j+1)
                ygc3e(i,j,ngz+k)=ygc3e(i,k,ngz-j+1)
                zgc3e(i,j,ngz+k)=zgc3e(i,k,ngz-j+1)
             enddo
          enddo
       enddo
    endif

    ! jmax-kmin
    ! ---------
    nv=bl(n)%BC(5) ! block kmin
    if (is_swapjk2_bl(5)) then
       if (is_rev2(5,3)) then
          nv_=bl(nv)%BC(5) ! block kmin of kmin
       else
          nv_=bl(nv)%BC(6) ! block kmax of kmin
       endif
    else
       nv_=bl(nv)%BC(4) ! block jmax of kmin
    endif
    !print *,'jmax-kmin',nv,nv_,bl(n)%BC(3),bl(n)%BC(4)
    !print *,'proc',coord(2),coord(3)
    if ((nv_>0).and.((nv_==bl(n)%BC(3)).or.(nv_==bl(n)%BC(4))).and.((coord(2)==ndomy-1).and.(coord(3)==0))) then
       !print *,'do it'
       do i=ndxt,nfxt
          do k=1-ngh,0
             do j=1,5
                xgc3e(i,ngy+j,k)=xgc3e(i,ngy+k,1-j)
                ygc3e(i,ngy+j,k)=ygc3e(i,ngy+k,1-j)
                zgc3e(i,ngy+j,k)=zgc3e(i,ngy+k,1-j)
             enddo
          enddo
       enddo
    endif

    ! jmax-kmax
    ! ---------
    nv=bl(n)%BC(6) ! block kmax
    if (is_swapjk2_bl(6)) then
       if (is_rev2(6,3)) then
          nv_=bl(nv)%BC(5) ! block kmin of kmax
       else
          nv_=bl(nv)%BC(6) ! block kmax of kmax
       endif
    else
       nv_=bl(nv)%BC(4) ! block jmax of kmax
    endif
    !print *,'jmax-kmax',nv,nv_,bl(n)%BC(3),bl(n)%BC(4)
    !print *,'proc',coord(2),coord(3)
    if ((nv_>0).and.((nv_==bl(n)%BC(3)).or.(nv_==bl(n)%BC(4))).and.((coord(2)==ndomy-1).and.(coord(3)==ndomz-1))) then
       !print *,'do it'
       do i=ndxt,nfxt
          do k=1,5
             do j=1,5
                xgc3e(i,ngy+j,ngz+k)=xgc3e(i,ngy-k+1,ngz+j)
                ygc3e(i,ngy+j,ngz+k)=ygc3e(i,ngy-k+1,ngz+j)
                zgc3e(i,ngy+j,ngz+k)=zgc3e(i,ngy-k+1,ngz+j)
             enddo
          enddo
       enddo
    endif

    ndyt=1
    nfyt=ngy
    if (bl(n)%BC(3)>0) ndyt=1-ngh
    if (bl(n)%BC(4)>0) nfyt=ngy+ngh

    ! imin-kmin
    ! ---------
    nv=bl(n)%BC(5) ! block kmin
    if (is_swapik2_bl(5)) then
       if (is_rev2(5,1)) then
          nv_=bl(nv)%BC(6) ! block kmax of kmin
       else
          nv_=bl(nv)%BC(5) ! block kmin of kmin
       endif
    else
       nv_=bl(nv)%BC(1) ! block imin of kmin
    endif
    !print *,'imin-kmin',nv,nv_,bl(n)%BC(1),bl(n)%BC(2)
    !print *,'proc',coord(1),coord(3)
    if ((nv_>0).and.((nv_==bl(n)%BC(1)).or.(nv_==bl(n)%BC(2))).and.((coord(1)==0).and.(coord(3)==0))) then
       !print *,'do it'
       do j=ndyt,nfyt
          do k=1-ngh,0
             do i=1-ngh,0
                xgc3e(i,j,k)=xgc3e(1-k,j,i)
                ygc3e(i,j,k)=ygc3e(1-k,j,i)
                zgc3e(i,j,k)=zgc3e(1-k,j,i)
             enddo
          enddo
       enddo
    endif

    ! imin-kmax
    ! ---------
    nv=bl(n)%BC(6) ! block kmax
    if (is_swapik2_bl(6)) then
       if (is_rev2(6,1)) then
          nv_=bl(nv)%BC(6) ! block kmax of kmax
       else
          nv_=bl(nv)%BC(5) ! block kmin of kmax
       endif
    else
       nv_=bl(nv)%BC(1) ! block imin of kmax
    endif
    !print *,'imin-kmax',nv,nv_,bl(n)%BC(1),bl(n)%BC(2)
    !print *,'proc',coord(1),coord(3)
    if ((nv_>0).and.((nv_==bl(n)%BC(1)).or.(nv_==bl(n)%BC(2))).and.((coord(1)==0).and.(coord(3)==ndomz-1))) then
       !print *,'do it'
       do j=ndyt,nfyt
          do k=1,5
             do i=1-ngh,0
                xgc3e(i,j,ngz+k)=xgc3e(k,j,ngz-i+1)
                ygc3e(i,j,ngz+k)=ygc3e(k,j,ngz-i+1)
                zgc3e(i,j,ngz+k)=zgc3e(k,j,ngz-i+1)
             enddo
          enddo
       enddo
    endif

    ! imax-kmin
    ! ---------
    nv=bl(n)%BC(5) ! block kmin
    if (is_swapik2_bl(5)) then
       if (is_rev2(5,1)) then
          nv_=bl(nv)%BC(5) ! block kmin of kmin
       else
          nv_=bl(nv)%BC(6) ! block kmax of kmin
       endif
    else
       nv_=bl(nv)%BC(2) ! block imax of kmin
    endif
    !print *,'imax-kmin',nv,nv_,bl(n)%BC(1),bl(n)%BC(2)
    !print *,'proc',coord(1),coord(3)
    if ((nv_>0).and.((nv_==bl(n)%BC(1)).or.(nv_==bl(n)%BC(2))).and.((coord(1)==ndomx-1).and.(coord(3)==0))) then
       !print *,'do it'
       do j=ndyt,nfyt
          do k=1-ngh,0
             do i=1,5
                xgc3e(ngx+i,j,k)=xgc3e(ngx+k,j,1-i)
                ygc3e(ngx+i,j,k)=ygc3e(ngx+k,j,1-i)
                zgc3e(ngx+i,j,k)=zgc3e(ngx+k,j,1-i)
             enddo
          enddo
       enddo
    endif

    ! imax-kmax
    ! ---------
    nv=bl(n)%BC(6) ! block kmax
    if (is_swapik2_bl(6)) then
       if (is_rev2(6,1)) then
          nv_=bl(nv)%BC(5) ! block kmin of kmax
       else
          nv_=bl(nv)%BC(6) ! block kmax of kmax
       endif       
    else
       nv_=bl(nv)%BC(2) ! block imax of kmax
    endif
    !print *,'jmax-kmax',nv,nv_,bl(n)%BC(1),bl(n)%BC(2)
    !print *,'proc',coord(1),coord(3)
    if ((nv_>0).and.((nv_==bl(n)%BC(1)).or.(nv_==bl(n)%BC(2))).and.((coord(1)==ndomx-1).and.(coord(3)==ndomz-1))) then
       !print *,'do it'
       do j=ndyt,nfyt
          do k=1,5
             do i=1,5
                xgc3e(ngx+i,j,ngz+k)=xgc3e(ngx-k+1,j,ngz+i)
                ygc3e(ngx+i,j,ngz+k)=ygc3e(ngx-k+1,j,ngz+i)
                zgc3e(ngx+i,j,ngz+k)=zgc3e(ngx-k+1,j,ngz+i)
             enddo
          enddo
       enddo
    endif

  end subroutine correct_edges_j
    
end module mod_add_gh3d
