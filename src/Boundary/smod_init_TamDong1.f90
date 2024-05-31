!===============================================================================
submodule (mod_init_TamDong) smod_init_TamDong1
!===============================================================================
  !> author: XG
  !> date: January 2024
  !> Initialization of Tam and Dong's boundary conditions
  !> routines: init_bc_TD, init_U0_TD, comm_U0_TD, mean0_U0_TD
!===============================================================================
  implicit none

contains

  !=============================================================================
  module subroutine init_bc_TD
  !=============================================================================
    !> Initialization of Tam & Dong's boundary conditions
  !=============================================================================
    use mod_RFM
    use mod_grid
    implicit none
    !---------------------------------------------------------------------------
    !integer :: i1,j1 ! radiation center
    !---------------------------------------------------------------------------
    integer :: i,j,jg,i_,j_
    !---------------------------------------------------------------------------

    if (iproc==0) write(*,*) repeat('=',70)
    if (iproc==0) print *,'Init Tam & Dong''s BC'

    ! Definition of radiation center (indices % global grid) TO BE CHANGED
    ! ------------------------------
!!$    i1=i0
!!$    j1=j0
!!$    !j1=1 ! for STBL case
!!$    if (is_curv) then
!!$       xcr=xgc(i1,j1)
!!$       ycr=ygc(i1,j1)
!!$    else
!!$       xcr=xg(i1)
!!$       ycr=yg(j1)
!!$    endif

    ! Defined in param.ini now
    ! xcr=xcr_ ! defined temporarily in grid
    ! ycr=ycr_ ! defined temporarily in grid
    ! zcr=zcr_ ! defined temporarily in grid
!!$    if (is_pulse) then
!!$       xcr=x_pulse
!!$       ycr=y_pulse
!!$    endif

    if (is_TamDong3D) then
       zcr=zcr_ ! defined temporarily in grid
       if (is_curv) zcr=zg(ngz/2)
       if (iproc==0) print *,'~> radiation center',xcr,ycr,zcr
    else
       if (iproc==0) print *,'~> radiation center',xcr,ycr
    endif

    ! Define local radiation center
    ! =============================
    ! BC imin (left in Cartesian)
    if (is_bc_TD(1,1)) then
       allocate(BC_face(1,1)%xrc(ny),BC_face(1,1)%yrc(ny))
       BC_face(1,1)%xrc=xcr
       BC_face(1,1)%yrc=ycr
       if (is_TamDong3D) then
          allocate(BC_face(1,1)%zrc(ny))
          BC_face(1,1)%zrc=zcr
       endif
    endif
    ! BC imax (right in Cartesian)
    if (is_bc_TD(1,2)) then
       allocate(BC_face(1,2)%xrc(ny),BC_face(1,2)%yrc(ny))
       BC_face(1,2)%xrc=xcr
       BC_face(1,2)%yrc=ycr
       if (is_TamDong3D) then
          allocate(BC_face(1,2)%zrc(ny))
          BC_face(1,2)%zrc=zcr
       endif
    endif
    ! BC jmin (bottom in Cartesian)
    if (is_bc_TD(2,1)) then
       allocate(BC_face(2,1)%xrc(nx),BC_face(2,1)%yrc(nx))
       BC_face(2,1)%xrc=xcr
       BC_face(2,1)%yrc=ycr
       if (is_TamDong3D) then
          allocate(BC_face(2,1)%zrc(nx))
          BC_face(2,1)%zrc=zcr
       endif
    endif
    ! BC jmax (top in Cartesian)
    if (is_bc_TD(2,2)) then
       allocate(BC_face(2,2)%xrc(nx),BC_face(2,2)%yrc(nx))
       BC_face(2,2)%xrc=xcr
       BC_face(2,2)%yrc=ycr
       if (is_TamDong3D) then
          allocate(BC_face(2,2)%zrc(nx))
          BC_face(2,2)%zrc=zcr
       endif
    endif
    if (is_TamDong3D) then
       ! BC kmin (front in Cartesian)
       if (is_bc_TD(3,1)) then
          allocate(BC_face(3,1)%xrc(nz),BC_face(3,1)%yrc(nz),BC_face(3,1)%zrc(nz))
          BC_face(3,1)%xrc=xcr
          BC_face(3,1)%yrc=ycr
          BC_face(3,1)%zrc=zcr
       endif
       ! BC kmax (back in Cartesian)
       if (is_bc_TD(3,2)) then
          allocate(BC_face(3,2)%xrc(nz),BC_face(3,2)%yrc(nz),BC_face(3,2)%zrc(nz))
          BC_face(3,2)%xrc=xcr
          BC_face(3,2)%yrc=ycr
          BC_face(3,2)%zrc=zcr
       endif
    endif

    ! ====================================
    ! Variable center of radiation for RFM
    ! ====================================
    ! To remain parallel to the boundary conditions
    ! ~> Specially important for the introduction of FST
    ! ~> Same realized for other boundary conditions to stay coherent
    ! /!\ If wall, at inlet/outlet, outflow angle needs to be parallel
    !     to the wall to avoid numerical issues, as observed in
    !     leading-edge case.
    if (is_RFM) then
       ! -----------------
       ! 3-D curv periodic ~> roughness, ...
       ! -----------------
       if (is_curv3) then
          if (is_bc_TD(1,1)) then
             do j=1,ny
                jg = j + coord(2)*ny
                if ((jg.le.ndy_rfm-5).and.(is_bc_wall(2,1))) then
                   BC_face(1,1)%xrc(j)=xc3(10,j,1)
                   BC_face(1,1)%yrc(j)=yc3(10,j,1)
                else if (jg.le.ndy_rfm-5) then
                   j_=min(ny,ndy_rfm-5)
                   BC_face(1,1)%xrc(j)=xc3(10,j_,1)
                   BC_face(1,1)%yrc(j)=yc3(10,j_,1)
                else if ((jg.ge.nfy_rfm+5).and.(is_bc_wall(2,2))) then
                   BC_face(1,1)%xrc(j)=xc3(10,j,1)
                   BC_face(1,1)%yrc(j)=yc3(10,j,1)
                else if (jg.ge.nfy_rfm+5) then
                   j_=max(nfy_rfm+5-coord(2)*ny,1)
                   BC_face(1,1)%xrc(j)=xc3(10,j_,1)
                   BC_face(1,1)%yrc(j)=yc3(10,j_,1)
                else
                   BC_face(1,1)%xrc(j)=xc3(10,j,1)
                   BC_face(1,1)%yrc(j)=yc3(10,j,1)
                endif
             enddo
          endif
          if (is_bc_TD(1,2)) then
             do j=1,ny
                jg = j + coord(2)*ny
                if ((jg.le.ndy_rfm-5).and.(is_bc_wall(2,1))) then
                   BC_face(1,2)%xrc(j)=xc3(nx-9,j,1)
                   BC_face(1,2)%yrc(j)=yc3(nx-9,j,1)
                else if (jg.le.ndy_rfm-5) then
                   j_=min(ny,ndy_rfm-5)
                   BC_face(1,2)%xrc(j)=xc3(nx-9,j_,1)
                   BC_face(1,2)%yrc(j)=yc3(nx-9,j_,1)
                else if ((jg.ge.nfy_rfm+5).and.(is_bc_wall(2,2))) then
                   BC_face(1,2)%xrc(j)=xc3(nx-9,j,1)
                   BC_face(1,2)%yrc(j)=yc3(nx-9,j,1)
                else if (jg.ge.nfy_rfm+5) then
                   j_=max(nfy_rfm+5-coord(2)*ny,1)
                   BC_face(1,2)%xrc(j)=xc3(nx-9,j_,1)
                   BC_face(1,2)%yrc(j)=yc3(nx-9,j_,1)
                else
                   BC_face(1,2)%xrc(j)=xc3(nx-9,j,1)
                   BC_face(1,2)%yrc(j)=yc3(nx-9,j,1)
                endif
             enddo
          endif
          if (is_bc_TD(2,1)) then
             j=1
             if (ndy_rfm.ne.1) j=min(ndy_rfm-5,ny)
             do i=1,nx
                if (is_bc_TD(1,1)) i_ = max(i,10)
                if (is_bc_TD(1,2)) i_ = min(i,nx-9)
                BC_face(2,1)%xrc(i)=xc3(i_,j,1)
                BC_face(2,1)%yrc(i)=yc3(i_,j,1)
             enddo
          endif
          if (is_bc_TD(2,2)) then
             j=ny
             if (nfy_rfm.ne.ngy) j=max(nfy_rfm+5-coord(2)*ny,1)
             do i=1,nx
                if (is_bc_TD(1,1)) i_ = max(i,10)
                if (is_bc_TD(1,2)) i_ = min(i,nx-9)
                BC_face(2,2)%xrc(i)=xc3(i_,j,1)
                BC_face(2,2)%yrc(i)=yc3(i_,j,1)
             enddo
          endif

       ! -------
       ! 2D curv
       ! -------
       else if (is_curv) then
          if (is_bc_TD(1,1)) then
             do j=1,ny
                jg = j + coord(2)*ny
                if ((jg.le.ndy_rfm-5).and.(is_bc_wall(2,1))) then
                   BC_face(1,1)%xrc(j)=xc(10,j)
                   BC_face(1,1)%yrc(j)=yc(10,j)
                else if (jg.le.ndy_rfm-5) then
                   j_=min(ny,ndy_rfm-5)
                   BC_face(1,1)%xrc(j)=xc(10,j_)
                   BC_face(1,1)%yrc(j)=yc(10,j_)
                else if ((jg.ge.nfy_rfm+5).and.(is_bc_wall(2,2))) then
                   BC_face(1,1)%xrc(j)=xc(10,j)
                   BC_face(1,1)%yrc(j)=yc(10,j)
                else if (jg.ge.nfy_rfm+5) then
                   j_=max(nfy_rfm+5-coord(2)*ny,1)
                   BC_face(1,1)%xrc(j)=xc(10,j_)
                   BC_face(1,1)%yrc(j)=yc(10,j_)
                else
                   BC_face(1,1)%xrc(j)=xc(10,j)
                   BC_face(1,1)%yrc(j)=yc(10,j)
                endif
             enddo
          endif
          if (is_bc_TD(1,2)) then
             do j=1,ny
                jg = j + coord(2)*ny
                if ((jg.le.ndy_rfm-5).and.(is_bc_wall(2,1))) then
                   BC_face(1,2)%xrc(j)=xc(nx-9,j)
                   BC_face(1,2)%yrc(j)=yc(nx-9,j)
                else if (jg.le.ndy_rfm-5) then
                   j_=min(ny,ndy_rfm-5)
                   BC_face(1,2)%xrc(j)=xc(nx-9,j_)
                   BC_face(1,2)%yrc(j)=yc(nx-9,j_)
                else if ((jg.ge.nfy_rfm+5).and.(is_bc_wall(2,2))) then
                   BC_face(1,2)%xrc(j)=xc(nx-9,j)
                   BC_face(1,2)%yrc(j)=yc(nx-9,j)
                else if (jg.ge.nfy_rfm+5) then
                   j_=max(nfy_rfm+5-coord(2)*ny,1)
                   BC_face(1,2)%xrc(j)=xc(nx-9,j_)
                   BC_face(1,2)%yrc(j)=yc(nx-9,j_)
                else
                   BC_face(1,2)%xrc(j)=xc(nx-9,j)
                   BC_face(1,2)%yrc(j)=yc(nx-9,j)
                endif
             enddo
          endif
          if (is_bc_TD(2,1)) then
             j=10
             if (ndy_rfm.ne.1) j=min(ndy_rfm-5,ny)
             do i=1,nx
                if (is_bc_TD(1,1)) i_ = max(i,10)
                if (is_bc_TD(1,2)) i_ = min(i,nx-9)
                BC_face(2,1)%xrc(i)=xc(i_,j)
                BC_face(2,1)%yrc(i)=yc(i_,j)
             enddo
          endif
          if (is_bc_TD(2,2)) then
             j=ny-10
             if (nfy_rfm.ne.ngy) j=max(nfy_rfm+5-coord(2)*ny,1)
             do i=1,nx
                if (is_bc_TD(1,1)) i_ = max(i,10)
                if (is_bc_TD(1,2)) i_ = min(i,nx-9)
                BC_face(2,2)%xrc(i)=xc(i_,j)
                BC_face(2,2)%yrc(i)=yc(i_,j)
             enddo
          endif

       ! ---------
       ! Cartesian grids
       ! ---------
       else
          if (is_bc_TD(1,1)) then
             BC_face(1,1)%xrc=x(10)
             do j=1,ny
                jg = j + coord(2)*ny
                if ((jg.le.ndy_rfm-5).and.(is_bc_wall(2,1))) then
                   BC_face(1,1)%yrc(j)=y(j)
                else if (jg.le.ndy_rfm-5) then
                   j_=min(ny,ndy_rfm-5)
                   BC_face(1,1)%yrc(j)=y(j_)
                else if ((jg.ge.nfy_rfm+5).and.(is_bc_wall(2,2))) then
                   BC_face(1,1)%yrc(j)=y(j)
                else if (jg.ge.nfy_rfm+5) then
                   j_=max(nfy_rfm+5-coord(2)*ny,1)
                   BC_face(1,1)%yrc(j)=y(j_)
                else
                   BC_face(1,1)%yrc(j)=y(j)
                endif
             enddo
          endif
          if (is_bc_TD(1,2)) then
             BC_face(1,2)%xrc=x(nx-9)
             do j=1,ny
                jg = j + coord(2)*ny
                if ((jg.le.ndy_rfm-5).and.(is_bc_wall(2,1))) then
                   BC_face(1,2)%yrc(j)=y(j)
                else if (jg.le.ndy_rfm-5) then
                   j_=min(ny,ndy_rfm-5)
                   BC_face(1,2)%yrc(j)=y(j_)
                else if ((jg.ge.nfy_rfm+5).and.(is_bc_wall(2,2))) then
                   BC_face(1,2)%yrc(j)=y(j)
                else if (jg.ge.nfy_rfm+5) then
                   j_=max(nfy_rfm+5-coord(2)*ny,1)
                   BC_face(1,2)%yrc(j)=y(j_)
                else
                   BC_face(1,2)%yrc(j)=y(j)
                endif
             enddo
          endif
          if (is_bc_TD(2,1)) then
             j=10
             if (ndy_rfm.ne.1) j=min(ndy_rfm-5,ny)
             BC_face(2,1)%yrc=y(j)
             do i=1,nx
                if (is_bc_TD(1,1)) i_ = max(i,10)
                if (is_bc_TD(1,2)) i_ = min(i,nx-9)
                BC_face(2,1)%xrc(i)=x(i_)
             enddo
          endif
          if (is_bc_TD(2,2)) then
             j=ny-10
             if (nfy_rfm.ne.ngy) j=max(nfy_rfm+5-coord(2)*ny,1)
             BC_face(2,2)%yrc=y(j)
             do i=1,nx
                if (is_bc_TD(1,1)) i_ = max(i,10)
                if (is_bc_TD(1,2)) i_ = min(i,nx-9)
                BC_face(2,2)%xrc(i)=x(i_)
             enddo
          endif

       endif
    endif

!!$    ! Variable center of radiation for 3-D curv periodic
!!$    ! ==================================================
!!$    if (iproc==0) then
!!$       if (is_bc_TD(1,1)) then
!!$          do k=1,nz
!!$             BC_face(1,1)%zrc(k)=zc3(1,1,k)
!!$          enddo
!!$       endif
!!$       if (is_bc_TD(1,2)) then
!!$          do k=1,nz
!!$             BC_face(1,2)%zrc(k)=zc3(1,1,k)
!!$          enddo
!!$       endif
!!$       if (is_bc_TD(2,1)) then
!!$          do k=1,nz
!!$             BC_face(2,1)%zrc(k)=zc3(1,1,k)
!!$          enddo
!!$       endif
!!$       if (is_bc_TD(2,2)) then
!!$          do k=1,nz
!!$             BC_face(2,2)%zrc(k)=zc3(1,1,k)
!!$          enddo
!!$       endif
!!$       if (is_bc_TD(3,1)) then
!!$          do k=1,nz
!!$             BC_face(3,1)%zrc(k)=zc3(1,1,k)
!!$          enddo
!!$       endif
!!$       if (is_bc_TD(3,2)) then
!!$          do k=1,nz
!!$             BC_face(3,2)%zrc(k)=zc3(1,1,k)
!!$          enddo
!!$       endif
!!$    endif
!!$    if (iproc==1) then
!!$       print *,nx,nz
!!$       if (is_bc_TD(1,1)) then
!!$          do i=1,nx
!!$             BC_face(1,1)%zrc(i)=zc3(i,1,1)
!!$          enddo
!!$       endif
!!$       if (is_bc_TD(1,2)) then
!!$          do i=1,nx
!!$             BC_face(1,2)%zrc(i)=zc3(i,1,1)
!!$          enddo
!!$       endif
!!$       if (is_bc_TD(2,1)) then
!!$          do i=1,nx
!!$             BC_face(2,1)%zrc(i)=zc3(i,1,1)
!!$          enddo
!!$       endif
!!$       if (is_bc_TD(2,2)) then
!!$          do i=1,nx
!!$             BC_face(2,2)%zrc(i)=zc3(i,1,1)
!!$          enddo
!!$       endif
!!$       if (is_bc_TD(3,1)) then
!!$          do i=1,nx
!!$             BC_face(3,1)%zrc(i)=zc3(i,1,1)
!!$          enddo
!!$       endif
!!$       if (is_bc_TD(3,2)) then
!!$          do i=1,nx
!!$             BC_face(3,2)%zrc(i)=zc3(i,1,1)
!!$          enddo
!!$       endif
!!$    endif

    ! Useful indices
    ! ==============
    nghp3=ngh+3
    nxmngh=nx-ngh
    nymngh=ny-ngh
    nzmngh=nz-ngh
    nxmnghp1=nx-ngh+1
    nymnghp1=ny-ngh+1
    nzmnghp1=nz-ngh+1

    ! For Tam & Dong on 1 point
    if (is_TamDong_1pt) then
       nghp3=5
       nxmngh=nx-2
       nymngh=ny-2
       nzmngh=nz-2
       nxmnghp1=nxmngh+1
       nymnghp1=nymngh+1
       nzmnghp1=nzmngh+1
    endif

    ! ! Initialize reference profiles <- Initialization based on initial field in mod_init_flow
    ! ! =============================
    ! if (is_read_ref) then
    !    call init_ref_TD
    ! endif
    ! ~> Done in mod_io_restart_BCref.f90 now

    
    ! Initialize time-averaged field
    ! ==============================
    if (idepart==FROM_SCRATCH) then
       call init_U0_TD
    elseif (idepart==FROM_FILE) then
       !call mean0_U0_TD
       call comm_U0_TD
    endif

    ! Initialize parameters for polar coordinates
    ! ===========================================
    if (is_TamDong3D) then
       call init_param_TD3d
    else if (is_TamDong_1pt) then
       call init_param_TD2d_1pt
    else
       call init_param_TD2d
    endif

    ! Free arrays for radiation center
    ! ================================
    ! BC imin
    if (is_bc_TD(1,1)) then
       deallocate(BC_face(1,1)%xrc,BC_face(1,1)%yrc)
       if (is_TamDong3D) deallocate(BC_face(1,1)%zrc)
    endif
    ! BC imax
    if (is_bc_TD(1,2)) then
       deallocate(BC_face(1,2)%xrc,BC_face(1,2)%yrc)
       if (is_TamDong3D) deallocate(BC_face(1,2)%zrc)
    endif
    ! BC jmin
    if (is_bc_TD(2,1)) then
       deallocate(BC_face(2,1)%xrc,BC_face(2,1)%yrc)
       if (is_TamDong3D) deallocate(BC_face(2,1)%zrc)
    endif
    ! BC jmax
    if (is_bc_TD(2,2)) then
       deallocate(BC_face(2,2)%xrc,BC_face(2,2)%yrc)
       if (is_TamDong3D) deallocate(BC_face(2,2)%zrc)
    endif
    if (is_TamDong3D) then
       ! BC kmin
       if (is_bc_TD(3,1)) then
          deallocate(BC_face(3,1)%xrc,BC_face(3,1)%yrc,BC_face(3,1)%zrc)
       endif
       ! BC kmax
       if (is_bc_TD(3,2)) then
          deallocate(BC_face(3,2)%xrc,BC_face(3,2)%yrc,BC_face(3,2)%zrc)
       endif
    endif
   
    if (iproc==0) print *,'~> init Tam & Dong OK'
    if (iproc==0) write(*,*) repeat('=',70)

  end subroutine init_bc_TD

  !=============================================================================
  module subroutine init_U0_TD
  !=============================================================================
    !> initialization of time-averaged field
  !=============================================================================
    implicit none
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    ! BC_face%U0 is allocated in Main/alloc.f90
    ! =========================================
    
    ! Initialize time-averaged primitive variables
    ! ============================================
    
    if (is_bc_TD(1,1)) then
       ! imin (left) boundary condition
       ! ------------------------------
       BC_face(1,1)%U0(:,:,:,1)=rho(1:nghp3,ny1:ny2,nz1:nz2)
       BC_face(1,1)%U0(:,:,:,2)= uu(1:nghp3,ny1:ny2,nz1:nz2)
       BC_face(1,1)%U0(:,:,:,3)= vv(1:nghp3,ny1:ny2,nz1:nz2)
       BC_face(1,1)%U0(:,:,:,4)= ww(1:nghp3,ny1:ny2,nz1:nz2)
       BC_face(1,1)%U0(:,:,:,5)=prs(1:nghp3,ny1:ny2,nz1:nz2)
       BC_face(1,1)%U0(:,:,:,6)= c_(1:nghp3,ny1:ny2,nz1:nz2)**2
    endif

    if (is_bc_TD(1,2)) then
       ! imax (right) boundary condition
       ! -------------------------------
       BC_face(1,2)%U0(:,:,:,1)=rho(nxmngh-2:nx,ny1:ny2,nz1:nz2)
       BC_face(1,2)%U0(:,:,:,2)= uu(nxmngh-2:nx,ny1:ny2,nz1:nz2)
       BC_face(1,2)%U0(:,:,:,3)= vv(nxmngh-2:nx,ny1:ny2,nz1:nz2)
       BC_face(1,2)%U0(:,:,:,4)= ww(nxmngh-2:nx,ny1:ny2,nz1:nz2)
       BC_face(1,2)%U0(:,:,:,5)=prs(nxmngh-2:nx,ny1:ny2,nz1:nz2)
       BC_face(1,2)%U0(:,:,:,6)= c_(nxmngh-2:nx,ny1:ny2,nz1:nz2)**2
    endif

    if (is_bc_TD(2,1)) then
       ! jmin (bottom) boundary condition
       ! --------------------------------
       BC_face(2,1)%U0(:,:,:,1)=rho(nx1:nx2,1:nghp3,nz1:nz2)
       BC_face(2,1)%U0(:,:,:,2)= uu(nx1:nx2,1:nghp3,nz1:nz2)
       BC_face(2,1)%U0(:,:,:,3)= vv(nx1:nx2,1:nghp3,nz1:nz2)
       BC_face(2,1)%U0(:,:,:,4)= ww(nx1:nx2,1:nghp3,nz1:nz2)
       BC_face(2,1)%U0(:,:,:,5)=prs(nx1:nx2,1:nghp3,nz1:nz2)
       BC_face(2,1)%U0(:,:,:,6)= c_(nx1:nx2,1:nghp3,nz1:nz2)**2
    endif

    if (is_bc_TD(2,2)) then
       ! jmax (top) boundary condition
       ! -----------------------------
       BC_face(2,2)%U0(:,:,:,1)=rho(nx1:nx2,nymngh-2:ny,nz1:nz2)
       BC_face(2,2)%U0(:,:,:,2)= uu(nx1:nx2,nymngh-2:ny,nz1:nz2)
       BC_face(2,2)%U0(:,:,:,3)= vv(nx1:nx2,nymngh-2:ny,nz1:nz2)
       BC_face(2,2)%U0(:,:,:,4)= ww(nx1:nx2,nymngh-2:ny,nz1:nz2)
       BC_face(2,2)%U0(:,:,:,5)=prs(nx1:nx2,nymngh-2:ny,nz1:nz2)
       BC_face(2,2)%U0(:,:,:,6)= c_(nx1:nx2,nymngh-2:ny,nz1:nz2)**2
    endif
    
    if (is_TamDong3D) then ! only for T&D 3D

       if (is_bc_TD(3,1)) then
          ! kmin (front) boundary condition
          ! -------------------------------
          BC_face(3,1)%U0(:,:,:,1)=rho(nx1:nx2,ny1:ny2,1:nghp3)
          BC_face(3,1)%U0(:,:,:,2)= uu(nx1:nx2,ny1:ny2,1:nghp3)
          BC_face(3,1)%U0(:,:,:,3)= vv(nx1:nx2,ny1:ny2,1:nghp3)
          BC_face(3,1)%U0(:,:,:,4)= ww(nx1:nx2,ny1:ny2,1:nghp3)
          BC_face(3,1)%U0(:,:,:,5)=prs(nx1:nx2,ny1:ny2,1:nghp3)
          BC_face(3,1)%U0(:,:,:,6)= c_(nx1:nx2,ny1:ny2,1:nghp3)**2
       endif

       if (is_bc_TD(3,2)) then
          ! kmax (back) boundary condition
          ! ------------------------------
          BC_face(3,2)%U0(:,:,:,1)=rho(nx1:nx2,ny1:ny2,nzmngh-2:nz)
          BC_face(3,2)%U0(:,:,:,2)= uu(nx1:nx2,ny1:ny2,nzmngh-2:nz)
          BC_face(3,2)%U0(:,:,:,3)= vv(nx1:nx2,ny1:ny2,nzmngh-2:nz)
          BC_face(3,2)%U0(:,:,:,4)= ww(nx1:nx2,ny1:ny2,nzmngh-2:nz)
          BC_face(3,2)%U0(:,:,:,5)=prs(nx1:nx2,ny1:ny2,nzmngh-2:nz)
          BC_face(3,2)%U0(:,:,:,6)= c_(nx1:nx2,ny1:ny2,nzmngh-2:nz)**2
       endif
    
    endif ! only for T&D 3D
    
  end subroutine init_U0_TD

  !=============================================================================
  module subroutine comm_U0_TD
  !=============================================================================
    !> initial communication of time-averaged field
  !=============================================================================
    use mod_interface
    implicit none
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
   
    if (is_bc_TD(1,1)) then
     rho_n(1:nghp3,:,:)=BC_face(1,1)%U0(:,:,:,1)
        uu(1:nghp3,:,:)=BC_face(1,1)%U0(:,:,:,2)
        vv(1:nghp3,:,:)=BC_face(1,1)%U0(:,:,:,3)
        ww(1:nghp3,:,:)=BC_face(1,1)%U0(:,:,:,4)
       prs(1:nghp3,:,:)=BC_face(1,1)%U0(:,:,:,5)
    endif
    if (is_bc_TD(1,2)) then
     rho_n(nxmngh-2:nx,:,:)=BC_face(1,2)%U0(:,:,:,1)
        uu(nxmngh-2:nx,:,:)=BC_face(1,2)%U0(:,:,:,2)
        vv(nxmngh-2:nx,:,:)=BC_face(1,2)%U0(:,:,:,3)
        ww(nxmngh-2:nx,:,:)=BC_face(1,2)%U0(:,:,:,4)
       prs(nxmngh-2:nx,:,:)=BC_face(1,2)%U0(:,:,:,5)
    endif
    if (is_bc_TD(2,1)) then
     rho_n(:,1:nghp3,:)=BC_face(2,1)%U0(:,:,:,1)
        uu(:,1:nghp3,:)=BC_face(2,1)%U0(:,:,:,2)
        vv(:,1:nghp3,:)=BC_face(2,1)%U0(:,:,:,3)
        ww(:,1:nghp3,:)=BC_face(2,1)%U0(:,:,:,4)
       prs(:,1:nghp3,:)=BC_face(2,1)%U0(:,:,:,5)
    endif
    if (is_bc_TD(2,2)) then
     rho_n(:,nymngh-2:ny,:)=BC_face(2,2)%U0(:,:,:,1)
        uu(:,nymngh-2:ny,:)=BC_face(2,2)%U0(:,:,:,2)
        vv(:,nymngh-2:ny,:)=BC_face(2,2)%U0(:,:,:,3)
        ww(:,nymngh-2:ny,:)=BC_face(2,2)%U0(:,:,:,4)
       prs(:,nymngh-2:ny,:)=BC_face(2,2)%U0(:,:,:,5)
    endif
    if (is_TamDong3D) then ! only for T&D 3D
       if (is_bc_TD(3,1)) then
        rho_n(:,:,1:nghp3)=BC_face(3,1)%U0(:,:,:,1)
           uu(:,:,1:nghp3)=BC_face(3,1)%U0(:,:,:,2)
           vv(:,:,1:nghp3)=BC_face(3,1)%U0(:,:,:,3)
           ww(:,:,1:nghp3)=BC_face(3,1)%U0(:,:,:,4)
          prs(:,:,1:nghp3)=BC_face(3,1)%U0(:,:,:,5)
       endif
       if (is_bc_TD(3,2)) then
        rho_n(:,:,nzmngh-2:nz)=BC_face(3,2)%U0(:,:,:,1)
           uu(:,:,nzmngh-2:nz)=BC_face(3,2)%U0(:,:,:,2)
           vv(:,:,nzmngh-2:nz)=BC_face(3,2)%U0(:,:,:,3)
           ww(:,:,nzmngh-2:nz)=BC_face(3,2)%U0(:,:,:,4)
          prs(:,:,nzmngh-2:nz)=BC_face(3,2)%U0(:,:,:,5)
       endif
    endif

    call communication_(rho_n,uu,vv,ww,prs)

    if (is_bc_TD(1,1)) then
       BC_face(1,1)%U0(:,:,:,1)=rho_n(1:nghp3,:,:)
       BC_face(1,1)%U0(:,:,:,2) =  uu(1:nghp3,:,:)
       BC_face(1,1)%U0(:,:,:,3) =  vv(1:nghp3,:,:)
       BC_face(1,1)%U0(:,:,:,4) =  ww(1:nghp3,:,:)
       BC_face(1,1)%U0(:,:,:,5) = prs(1:nghp3,:,:)
    endif
    if (is_bc_TD(1,2)) then
       BC_face(1,2)%U0(:,:,:,1)=rho_n(nxmngh-2:nx,:,:)
       BC_face(1,2)%U0(:,:,:,2) =  uu(nxmngh-2:nx,:,:)
       BC_face(1,2)%U0(:,:,:,3) =  vv(nxmngh-2:nx,:,:)
       BC_face(1,2)%U0(:,:,:,4) =  ww(nxmngh-2:nx,:,:)
       BC_face(1,2)%U0(:,:,:,5) = prs(nxmngh-2:nx,:,:)
    endif
    if (is_bc_TD(2,1)) then
       BC_face(2,1)%U0(:,:,:,1)=rho_n(:,1:nghp3,:)
       BC_face(2,1)%U0(:,:,:,2) =  uu(:,1:nghp3,:)
       BC_face(2,1)%U0(:,:,:,3) =  vv(:,1:nghp3,:)
       BC_face(2,1)%U0(:,:,:,4) =  ww(:,1:nghp3,:)
       BC_face(2,1)%U0(:,:,:,5) = prs(:,1:nghp3,:)
    endif
    if (is_bc_TD(2,2)) then
       BC_face(2,2)%U0(:,:,:,1)=rho_n(:,nymngh-2:ny,:)
       BC_face(2,2)%U0(:,:,:,2) =  uu(:,nymngh-2:ny,:)
       BC_face(2,2)%U0(:,:,:,3) =  vv(:,nymngh-2:ny,:)
       BC_face(2,2)%U0(:,:,:,4) =  ww(:,nymngh-2:ny,:)
       BC_face(2,2)%U0(:,:,:,5) = prs(:,nymngh-2:ny,:)
    endif
    if (is_TamDong3D) then ! only for T&D 3D
       if (is_bc_TD(3,1)) then
          BC_face(3,1)%U0(:,:,:,1)=rho_n(:,:,1:nghp3)
          BC_face(3,1)%U0(:,:,:,2) =  uu(:,:,1:nghp3)
          BC_face(3,1)%U0(:,:,:,3) =  vv(:,:,1:nghp3)
          BC_face(3,1)%U0(:,:,:,4) =  ww(:,:,1:nghp3)
          BC_face(3,1)%U0(:,:,:,5) = prs(:,:,1:nghp3)
       endif
       if (is_bc_TD(3,2)) then
          BC_face(3,2)%U0(:,:,:,1)=rho_n(:,:,nzmngh-2:nz)
          BC_face(3,2)%U0(:,:,:,2) =  uu(:,:,nzmngh-2:nz)
          BC_face(3,2)%U0(:,:,:,3) =  vv(:,:,nzmngh-2:nz)
          BC_face(3,2)%U0(:,:,:,4) =  ww(:,:,nzmngh-2:nz)
          BC_face(3,2)%U0(:,:,:,5) = prs(:,:,nzmngh-2:nz)
       endif
    endif
   
  end subroutine comm_U0_TD

  !=============================================================================
  module subroutine mean0_U0_TD
  !=============================================================================
    !> initial communication of time-averaged field
  !=============================================================================
    use mod_flow0
    implicit none
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    if (is_bc_TD(1,1)) then
       BC_face(1,1)%U0(:,:,:,1)=rho0(1:nghp3,:,:)
       BC_face(1,1)%U0(:,:,:,2) = u0(1:nghp3,:,:)
       BC_face(1,1)%U0(:,:,:,3) = v0(1:nghp3,:,:)
       BC_face(1,1)%U0(:,:,:,4) = w0(1:nghp3,:,:)
       BC_face(1,1)%U0(:,:,:,5) = p0(1:nghp3,:,:)
       BC_face(1,1)%U0(:,:,:,6)=gam*rg*T0(1:nghp3,:,:)
    endif
    if (is_bc_TD(1,2)) then
       BC_face(1,2)%U0(:,:,:,1)=rho0(nxmngh-2:nx,:,:)
       BC_face(1,2)%U0(:,:,:,2) = u0(nxmngh-2:nx,:,:)
       BC_face(1,2)%U0(:,:,:,3) = v0(nxmngh-2:nx,:,:)
       BC_face(1,2)%U0(:,:,:,4) = w0(nxmngh-2:nx,:,:)
       BC_face(1,2)%U0(:,:,:,5) = p0(nxmngh-2:nx,:,:)
       BC_face(1,2)%U0(:,:,:,6)=gam*rg*T0(nxmngh-2:nx,:,:)
    endif
    if (is_bc_TD(2,1)) then
       BC_face(2,1)%U0(:,:,:,1)=rho0(:,1:nghp3,:)
       BC_face(2,1)%U0(:,:,:,2) = u0(:,1:nghp3,:)
       BC_face(2,1)%U0(:,:,:,3) = v0(:,1:nghp3,:)
       BC_face(2,1)%U0(:,:,:,4) = w0(:,1:nghp3,:)
       BC_face(2,1)%U0(:,:,:,5) = p0(:,1:nghp3,:)
       BC_face(2,1)%U0(:,:,:,6)=gam*rg*T0(:,1:nghp3,:)
    endif
    if (is_bc_TD(2,2)) then
       BC_face(2,2)%U0(:,:,:,1)=rho0(:,nymngh-2:ny,:)
       BC_face(2,2)%U0(:,:,:,2) = u0(:,nymngh-2:ny,:)
       BC_face(2,2)%U0(:,:,:,3) = v0(:,nymngh-2:ny,:)
       BC_face(2,2)%U0(:,:,:,4) = w0(:,nymngh-2:ny,:)
       BC_face(2,2)%U0(:,:,:,5) = p0(:,nymngh-2:ny,:)
       BC_face(2,2)%U0(:,:,:,6)=gam*rg*T0(:,nymngh-2:ny,:)
    endif
   
    deallocate(rho0,p0,T0,u0,v0,w0)
    is_mean0=.false.

  end subroutine mean0_U0_TD
  
end submodule smod_init_TamDong1
