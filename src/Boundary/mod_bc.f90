!===============================================================
module mod_bc
!===============================================================
  !> Module for Boundary Conditions (BC)
!===============================================================
  use precision
  implicit none
  !---------------------------------------------------------------------------
  type border_face
     integer :: sort
     ! reference quantities imposed at boundary
     logical :: is_mean_ref
     ! number of IRS ghost points in MPI domain
     integer, dimension(6) :: ngh_irs
     ! CFL condition in interface for IRS scheme
     real(wp) :: cflmax(3)
     ! radiation center for Tam & Dong BC
     real(wp), dimension(:), pointer :: xrc,yrc,zrc
     ! 2D (periodic) polar coordinates for Tam & Dong BC
     real(wp), dimension(:,:), pointer :: ir,cosphi,sinphi
     ! 3D spherical coordinates for Tam & Dong BC
     real(wp), dimension(:,:,:), pointer :: i_r,cosp,sinp,cost,sint
     real(wp), dimension(:,:,:), pointer :: costcosp,costsinp,sintcosp,sintsinp
     ! Reference 2D profiles for boundary conditions
     ! [last index 1...10 corresponds to rho,u,v,w,p,c2,rhou,rhov,rhow,rhoe]
     real(wp), dimension(:,:,:), pointer :: Uref
     ! time-averaged primitive variables for Tam & Dong BC
     ! [last index 1...6 corresponds to rho0,u0,v0,w0,p0,c02]
     real(wp), dimension(:,:,:,:), pointer :: U0
  end type border_face
  !---------------------------------------------------------------------------
  type border_edge
     integer :: sort
     ! 2D(periodic) polar coordinates for Tam & Dong BC
     real(wp), dimension(:,:), pointer :: ir,cosphi,sinphi
     ! 3D spherical coordinates for Tam & Dong BC
     real(wp), dimension(:,:,:), pointer :: i_r,cosp,sinp,cost,sint
     real(wp), dimension(:,:,:), pointer :: costcosp,costsinp,sintcosp,sintsinp
  end type border_edge
  !---------------------------------------------------------------------------
  type border_corner
     integer :: sort
     ! 3D spherical coordinates for Tam & Dong BC
     real(wp), dimension(:,:,:), pointer :: i_r,cosp,sinp,cost,sint
     real(wp), dimension(:,:,:), pointer :: costcosp,costsinp,sintcosp,sintsinp
  end type border_corner
  !---------------------------------------------------------------------------
  ! Boundary conditions
  logical :: is_boundary(3,2),is_TamDong,is_TamDong3D,is_TamDong_1pt
  logical :: is_inlet_outlet
  logical :: is_bc_wall(3,2),is_bc_wall2(3,2)
  logical :: is_bc_TD(3,2),is_bc_1pt(3,2)
  logical :: is_bc_slip,is_slip(3,2) ! slip wall BC
  type(border_face), dimension(3,2) :: BC_face
  type(border_edge), dimension(3,2,2) :: BC_edge
  type(border_corner), dimension(2,2,2) :: BC_corner  
  !---------------------------------------------------------------------------
  real(wp) :: T_wall ! TO BE CHANGED ???? -> mod_constant ?
  !---------------------------------------------------------------------------
  logical :: is_BC_ref,is_read_ref,is_BCref_init ! TO BE CHANGED ??

contains

  !===============================================================================
  subroutine bc_define
  !===============================================================================
  !* Define Boundary Conditions
  ! ----------------------------
  ! Types of boundary conditions
  ! ----------------------------
  !    0 : wall BC
  !   -1 : Tam & Dong's radiation BC
  !  -11 : Tam & Dong's radiation BC on 1 point
  !   -2 : Tam & Dong's outflow BC
  !   -3 : non-reflecting characteristics BC
  !   -4 : inlet riemann conditions
  !   -5 : outlet riemann conditions
  !    1 : interior points
  !    2 : overlapping (interpolation)
  !*
  !===============================================================================
    use mod_mpi
    use mod_block
    use mod_grid
    use mod_constant
    use warnstop
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: nbl
    integer :: i,j,k
    integer :: nbl_W,nbl_E,nbl_S,nbl_N
    logical :: is_TamDong_proc,is_TamDong_1pt_proc,is_BCref_proc
    logical :: is_inlet_outlet_proc,is_bc_slip_proc
    logical :: is_errBC,is_errBC_proc
    ! ---------------------------------------------------------------------------

    ! Definition of faces
    ! ===================
    if (iproc==0) print *,'init BC faces'

    ! Initializations
    ! ---------------
    is_BCref_proc=.false.
    is_slip=.false.
    do i=1,3
       do j=1,2
          BC_face(i,j)%sort=1
          is_boundary(i,j)=.false.
          BC_face(i,j)%is_mean_ref=.false.
       enddo
    enddo

    ! Determine face BC
    ! -----------------
    nbl=nob(iproc)
    if (coord(1)==0) then
       if (bl(nbl)%BC(1)<=0) then
          BC_face(1,1)%sort=bl(nbl)%BC(1)
          is_boundary(1,1)=.true.
          if (bl(nbl)%flag(1)=='r') then
             is_BCref_proc=.true.
             BC_face(1,1)%is_mean_ref=.true.
          endif
          if (bl(nbl)%flag(1)=='s') is_slip(1,1)=.true.
       endif
    endif
    if (coord(1)==ndomx-1) then
       if (bl(nbl)%BC(2)<=0) then
          BC_face(1,2)%sort=bl(nbl)%BC(2)
          is_boundary(1,2)=.true.
          if (bl(nbl)%flag(2)=='r') then
             is_BCref_proc=.true.
             BC_face(1,2)%is_mean_ref=.true.
          endif
          if (bl(nbl)%flag(2)=='s') is_slip(1,2)=.true.
       endif
    endif
    if (coord(2)==0) then
       if (bl(nbl)%BC(3)<=0) then
          BC_face(2,1)%sort=bl(nbl)%BC(3)
          is_boundary(2,1)=.true.
          if (bl(nbl)%flag(3)=='r') then
             is_BCref_proc=.true.
             BC_face(2,1)%is_mean_ref=.true.
          endif
          if (bl(nbl)%flag(3)=='s') is_slip(2,1)=.true.
       endif
    endif
    if (coord(2)==ndomy-1) then
       if (bl(nbl)%BC(4)<=0) then
          BC_face(2,2)%sort=bl(nbl)%BC(4)
          is_boundary(2,2)=.true.
          if (bl(nbl)%flag(4)=='r') then
             is_BCref_proc=.true.
             BC_face(2,2)%is_mean_ref=.true.
          endif
          if (bl(nbl)%flag(4)=='s') is_slip(2,2)=.true.
       endif
    endif
    if (coord(3)==0) then
       if (bl(nbl)%BC(5)<=0) then
          BC_face(3,1)%sort=bl(nbl)%BC(5)
          is_boundary(3,1)=.true.
          if (bl(nbl)%flag(5)=='r') then
             is_BCref_proc=.true.
             BC_face(3,1)%is_mean_ref=.true.
          endif
          if (bl(nbl)%flag(5)=='s') is_slip(3,1)=.true.
       endif
    endif
    if (coord(3)==ndomz-1) then
       if (bl(nbl)%BC(6)<=0) then
          BC_face(3,2)%sort=bl(nbl)%BC(6)
          is_boundary(3,2)=.true.
          if (bl(nbl)%flag(6)=='r') then
             is_BCref_proc=.true.
             BC_face(3,2)%is_mean_ref=.true.
          endif
          if (bl(nbl)%flag(6)=='s') is_slip(3,2)=.true.
       endif
    endif

    is_bc_slip_proc= is_slip(1,1).or.is_slip(1,2) &
                 .or.is_slip(2,1).or.is_slip(2,2) &
                 .or.is_slip(3,1).or.is_slip(3,2)

    ! Set boolean indicator is_BC_ref & is_bc_slip
    ! --------------------------------------------
    call MPI_ALLREDUCE(is_BCref_proc,is_BC_ref,1,MPI_LOGICAL,MPI_LOR,COMM_global,info)
    call MPI_ALLREDUCE(is_bc_slip_proc,is_bc_slip,1,MPI_LOGICAL,MPI_LOR,COMM_global,info)

    ! Determine indicators is_bc_wall & is_bc_1pt
    ! -------------------------------------------
    is_bc_wall=.false.
    is_bc_1pt =.false.
    do i=1,3
       do j=1,2
          if (BC_face(i,j)%sort== 0) is_bc_wall(i,j)=.true.
          if (BC_face(i,j)%sort<=-3) is_bc_1pt(i,j) =.true.
          if (is_bc_wall(i,j)) is_bc_1pt(i,j) =.true.
       enddo
    enddo

    ! is_bc_wall2: advancement of wall boundary points
    if (is_wall2) then
       is_bc_wall2=is_bc_wall
    else
       is_bc_wall2=.false.
    endif
    
    ! Determine indicator is_inlet_outlet
    ! --------------------------------------
    is_inlet_outlet=.false.
    do i=1,3
       do j=1,2
          if ((BC_face(i,j)%sort==-4).or.(BC_face(i,j)%sort==-5)) is_inlet_outlet_proc=.true.
       enddo
    enddo
    call MPI_ALLREDUCE(is_inlet_outlet_proc,is_inlet_outlet,1,MPI_LOGICAL,MPI_LOR,COMM_global,info)

    ! Determine indicators is_TamDong
    ! -------------------------------
    ! true if there is one BC among block that use Tam & Dong's BC
    is_TamDong_proc=.false.
    is_TamDong_1pt_proc=.false.
    do i=1,3
       do j=1,2
          if ((BC_face(i,j)%sort==-1).or.(BC_face(i,j)%sort==-2).or.(BC_face(i,j)%sort==-11)) then
             is_TamDong_proc=.true.
             is_bc_TD(i,j)=.true.
          endif
          if (BC_face(i,j)%sort==-11) is_TamDong_1pt_proc=.true.
       enddo
    enddo
    call MPI_ALLREDUCE(is_TamDong_proc,is_TamDong,1,MPI_LOGICAL,MPI_LOR,COMM_global,info)
    call MPI_ALLREDUCE(is_TamDong_1pt_proc,is_TamDong_1pt,1,MPI_LOGICAL,MPI_LOR,COMM_global,info)
    is_TamDong3D=.false.
    !is_TamDong3D=.true.
    if (is_curv3) is_TamDong3D=.true.

    ! Compute mean field for Tam & Dong BC
    ! ====================================
    is_mean0=.false.
    !is_mean0=.true.
    ! old version [commented]
    !! if (is_TamDong) is_mean0=.true.
   
!!$    !print *,'iproc',iproc,'BC',bl(nbl)%BC(1),bl(nbl)%BC(2),bl(nbl)%BC(3),bl(nbl)%BC(4)
!!$    print *,'iproc',iproc,'BC_face',BC_face(1,1)%sort,BC_face(1,2)%sort,BC_face(2,1)%sort,BC_face(2,2)%sort
!!$    call mpistop('stop BC!', 0)

!!$    do i=1,3
!!$    !do i=1,2
!!$       do j=1,2
!!$          print *,'iproc',iproc,'BC_face',i,j,BC_face(i,j)%sort,is_bc_TD(i,j)
!!$       enddo
!!$    enddo
!!$    call mpistop('stop BC!', 0)

!!$    !do i=1,3
!!$    do i=1,2
!!$       do j=1,2
!!$          print *,'iproc',iproc,'BC_face',i,j,BC_face(i,j)%sort
!!$       enddo
!!$    enddo

    ! Definition of edges
    ! ===================
    ! ----------------------------------------------------
    ! 1,1,1 : imin-jmin / left-bottom
    ! 1,1,2 : imin-jmax / left-top
    ! 1,2,1 : imax-jmin / right-bottom
    ! 1,2,2 : imax-jmax / right-top
    ! 2,1,1 : jmin-kmin / bottom-front
    ! 2,1,2 : jmin-kmax / bottom-back
    ! 2,2,1 : jmax-kmin / top-front
    ! 2,2,2 : jmax-kmax / top-back
    ! 3,1,1 : kmin-imin / front-left
    ! 3,1,2 : kmin-imax / front-right
    ! 3,2,1 : kmax-imin / back-left
    ! 3,2,2 : kmax-imax / back-right
    ! ----------------------------------------------------
    if (iproc==0) print *,'init BC edges'

    ! Initializations
    ! ---------------
    do i=1,3
       do j=1,2
          do k=1,2
             BC_edge(i,j,k)%sort=1
          enddo
       enddo
    enddo

    ! Edges along z (face 1: left)
    ! ----------------------------
    ! bottom
    if ((is_boundary(1,1)).and.(is_boundary(2,1))) then
       BC_edge(1,1,1)%sort=min(BC_face(1,1)%sort,BC_face(2,1)%sort)
    endif
    ! top
    if ((is_boundary(1,1)).and.(is_boundary(2,2))) then
       BC_edge(1,1,2)%sort=min(BC_face(1,1)%sort,BC_face(2,2)%sort)
    endif

    ! Edges along z (face 2: right)
    ! -----------------------------
    ! bottom
    if ((is_boundary(1,2)).and.(is_boundary(2,1))) then
       BC_edge(1,2,1)%sort=min(BC_face(1,2)%sort,BC_face(2,1)%sort)
    endif
    ! top
    if ((is_boundary(1,2)).and.(is_boundary(2,2))) then
       BC_edge(1,2,2)%sort=min(BC_face(1,2)%sort,BC_face(2,2)%sort)
    endif

    ! Edges along x (face 1: bottom)
    ! ------------------------------
    ! front
    if ((is_boundary(2,1)).and.(is_boundary(3,1))) then
       BC_edge(2,1,1)%sort=min(BC_face(2,1)%sort,BC_face(3,1)%sort)
    endif
    ! back
    if ((is_boundary(2,1)).and.(is_boundary(3,2))) then
       BC_edge(2,1,2)%sort=min(BC_face(2,1)%sort,BC_face(3,2)%sort)
    endif

    ! Edges along x (face 2: top)
    ! ---------------------------
    ! front
    if ((is_boundary(2,2)).and.(is_boundary(3,1))) then
       BC_edge(2,2,1)%sort=min(BC_face(2,2)%sort,BC_face(3,1)%sort)
    endif
    ! back
    if ((is_boundary(2,2)).and.(is_boundary(3,2))) then
       BC_edge(2,2,2)%sort=min(BC_face(2,2)%sort,BC_face(3,2)%sort)
    endif

    ! Edges along y (face 1: front)
    ! -----------------------------
    ! left
    if ((is_boundary(3,1)).and.(is_boundary(1,1))) then
       BC_edge(3,1,1)%sort=min(BC_face(3,1)%sort,BC_face(1,1)%sort)
    endif
    ! right
    if ((is_boundary(3,1)).and.(is_boundary(1,2))) then
       BC_edge(3,1,2)%sort=min(BC_face(3,1)%sort,BC_face(1,2)%sort)
    endif

    ! Edges along y (face 2: back)
    ! ----------------------------
    ! left
    if ((is_boundary(3,2)).and.(is_boundary(1,1))) then
       BC_edge(3,2,1)%sort=min(BC_face(3,2)%sort,BC_face(1,1)%sort)
    endif
    ! right
    if ((is_boundary(3,2)).and.(is_boundary(1,2))) then
       BC_edge(3,2,2)%sort=min(BC_face(3,2)%sort,BC_face(1,2)%sort)
    endif

    ! Determination of edge imin-jmin
    ! -------------------------------
    if ((coord(1)==0).and.(coord(2)==0).and.(bl(nbl)%BC(1)>0).and.(bl(nbl)%BC(3)>0)) then
       nbl_W=bl(nbl)%BC(1)
       nbl_S=bl(nbl)%BC(3)
       if ((bl(nbl_W)%BC(3)==0).and.(bl(nbl_S)%BC(1)==0)) BC_edge(1,1,1)%sort=2
       if ((bl(nbl_W)%BC(3)==0).and.(bl(nbl_S)%BC(1)==0)) print *,iproc,'BC_edge(1,1,1)',BC_edge(1,1,1)%sort
    endif

    ! Determination of edge imax-jmin
    ! -------------------------------
    if ((coord(1)==ndomx-1).and.(coord(2)==0).and.(bl(nbl)%BC(2)>0).and.(bl(nbl)%BC(3)>0)) then
       nbl_E=bl(nbl)%BC(2)
       nbl_S=bl(nbl)%BC(3)
       if ((bl(nbl_E)%BC(3)==0).and.(bl(nbl_S)%BC(2)==0)) BC_edge(1,2,1)%sort=2
       if ((bl(nbl_E)%BC(3)==0).and.(bl(nbl_S)%BC(2)==0)) print *,iproc,'BC_edge(1,2,1)',BC_edge(1,2,1)%sort
    endif

    ! Determination of edge imin-jmax
    ! -------------------------------
    if ((coord(1)==0).and.(coord(2)==ndomy-1).and.(bl(nbl)%BC(1)>0).and.(bl(nbl)%BC(4)>0)) then
       nbl_W=bl(nbl)%BC(1)
       nbl_N=bl(nbl)%BC(4)
       if ((bl(nbl_W)%BC(4)==0).and.(bl(nbl_N)%BC(1)==0)) BC_edge(1,1,2)%sort=2
       if ((bl(nbl_W)%BC(4)==0).and.(bl(nbl_N)%BC(1)==0)) print *,iproc,'BC_edge(1,1,2)',BC_edge(1,1,2)%sort
    endif

    ! Determination of edge imax-jmax
    ! -------------------------------
    if ((coord(1)==ndomx-1).and.(coord(2)==ndomy-1).and.(bl(nbl)%BC(2)>0).and.(bl(nbl)%BC(4)>0)) then
       nbl_E=bl(nbl)%BC(2)
       nbl_N=bl(nbl)%BC(4)
       if ((bl(nbl_E)%BC(4)==0).and.(bl(nbl_N)%BC(2)==0)) BC_edge(1,2,2)%sort=2
       if ((bl(nbl_E)%BC(4)==0).and.(bl(nbl_N)%BC(2)==0)) print *,iproc,'BC_edge(1,2,2)',BC_edge(1,2,2)%sort
    endif

!!$    ! ********************************************************
!!$    ! ATTENTION TEMPORARY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    ! ********************************************************
!!$    if (TURB) then
!!$       if (nbloc==4) then
!!$          if (iproc==0) then
!!$             ! edge imin-jmin
!!$             BC_edge(1,1,1)%sort=2
!!$          endif
!!$
!!$          if (iproc==3) then
!!$             ! edge imin-jmax
!!$             BC_edge(1,1,2)%sort=2
!!$          endif
!!$       endif
!!$       if (nbloc==6) then
!!$          if (iproc==0) then
!!$             ! edge imin-jmin
!!$             BC_edge(1,1,1)%sort=2
!!$          endif
!!$
!!$          if (iproc==5) then
!!$             ! edge imin-jmax
!!$             BC_edge(1,1,2)%sort=2
!!$          endif
!!$
!!$          if (iproc==2) then
!!$             ! edge imax-jmin
!!$             BC_edge(1,2,1)%sort=2
!!$          endif
!!$
!!$          if (iproc==3) then
!!$             ! edge imax-jmax
!!$             BC_edge(1,2,2)%sort=2
!!$          endif
!!$       endif
!!$    endif
    ! ********************************************************

!!$    !do i=1,3
!!$    do i=1,1
!!$       do j=1,2
!!$          do k=1,2
!!$             print *,'iproc',iproc,'BC_edge',i,j,k,BC_edge(i,j,k)%sort
!!$          enddo
!!$       enddo
!!$    enddo
!!$    call mpistop('stop BC!', 0)

    ! Definition of corners
    ! ===================
    ! ----------------------------------------------------
    ! 1,1,1 : imin-jmin-kmin / left-bottom-front
    ! 1,1,2 : imin-jmin-kmax / left-bottom-back
    ! 1,2,1 : imin-jmax-kmin / left-top-front
    ! 1,2,2 : imin-jmax-kmax / left-top-back
    ! 2,1,1 : imax-jmin-kmin / right-bottom-front
    ! 2,1,2 : imax-jmin-kmax / right-bottom-back
    ! 2,2,1 : imax-jmax-kmin / right-top-front
    ! 2,2,2 : imax-jmax-kmax / right-top-back
    ! ----------------------------------------------------
    if (iproc==0) print *,'init BC corners'

    ! Initializations
    ! ---------------
    do i=1,2
       do j=1,2
          do k=1,2
             BC_corner(i,j,k)%sort=1
          enddo
       enddo
    enddo

    ! Corner left-bottom-front
    ! ------------------------
    if ((is_boundary(1,1)).and.(is_boundary(2,1)).and.(is_boundary(3,1))) then
       BC_corner(1,1,1)%sort=min(BC_face(1,1)%sort,BC_face(2,1)%sort,BC_face(3,1)%sort)
    endif
    ! Corner left-top-front
    ! ---------------------
    if ((is_boundary(1,1)).and.(is_boundary(2,2)).and.(is_boundary(3,1))) then
       BC_corner(1,2,1)%sort=min(BC_face(1,1)%sort,BC_face(2,2)%sort,BC_face(3,1)%sort)
    endif
    ! Corner right-bottom-front
    ! -------------------------
    if ((is_boundary(1,2)).and.(is_boundary(2,1)).and.(is_boundary(3,1))) then
       BC_corner(2,1,1)%sort=min(BC_face(1,2)%sort,BC_face(2,1)%sort,BC_face(3,1)%sort)
    endif
    ! Corner right-top-front
    ! ----------------------
    if ((is_boundary(1,2)).and.(is_boundary(2,2)).and.(is_boundary(3,1))) then
       BC_corner(2,2,1)%sort=min(BC_face(1,2)%sort,BC_face(2,2)%sort,BC_face(3,1)%sort)
    endif
    ! Corner left-bottom-back
    ! -----------------------
    if ((is_boundary(1,1)).and.(is_boundary(2,1)).and.(is_boundary(3,2))) then
       BC_corner(1,1,2)%sort=min(BC_face(1,1)%sort,BC_face(2,1)%sort,BC_face(3,2)%sort)
    endif
    ! Corner left-top-back
    ! --------------------
    if ((is_boundary(1,1)).and.(is_boundary(2,2)).and.(is_boundary(3,2))) then
       BC_corner(1,2,2)%sort=min(BC_face(1,1)%sort,BC_face(2,2)%sort,BC_face(3,2)%sort)
    endif
    ! Corner right-bottom-back
    ! ------------------------
    if ((is_boundary(1,2)).and.(is_boundary(2,1)).and.(is_boundary(3,2))) then
       BC_corner(2,1,2)%sort=min(BC_face(1,2)%sort,BC_face(2,1)%sort,BC_face(3,2)%sort)
    endif
    ! Corner right-top-back
    ! ---------------------
    if ((is_boundary(1,2)).and.(is_boundary(2,2)).and.(is_boundary(3,2))) then
       BC_corner(2,2,2)%sort=min(BC_face(1,2)%sort,BC_face(2,2)%sort,BC_face(3,2)%sort)
    endif

!!$    do i=1,2
!!$       do j=1,2
!!$          do k=1,2
!!$             print *,'iproc',iproc,'BC_corner',i,j,k,BC_corner(i,j,k)%sort
!!$          enddo
!!$       enddo
!!$    enddo
!!$    !!call mpistop('problem definition of corner # 1', 0)

    ! Protection and consistency check for supersonic flow
    ! ====================================================
    if (Mach.ge.1) then
       ! 1. supersonic flow is not compatible with Tam&Dong's BCs
       ! --------------------------------------------------------
       is_errBC_proc=.false.
       if ((BC_face(1,1)%sort==-1).or.(BC_face(2,1)%sort==-1).or. &
           (BC_face(2,1)%sort==-1).or.(BC_face(2,2)%sort==-1)) then
          is_errBC_proc=.true.
       endif
       call MPI_ALLREDUCE(is_errBC_proc,is_errBC,1,MPI_LOGICAL,MPI_LOR,COMM_global,info)
       if (is_errBC) &
            call mpistop('You cannot use Tam&Dong''s BCs for supersonic flow. Use characteristics (-3)',0)
       
!!$       ! 2. supersonic inlet uses BC_face(1,1)%Uref: need 'r' flag
!!$       ! ---------------------------------------------------------
!!$       is_errBC_proc=.false.
!!$       if ((BC_face(1,1)%sort==-3).and.(.not.BC_face(1,1)%is_mean_ref)) then
!!$          is_errBC_proc=.true.
!!$       endif
!!$       call MPI_ALLREDUCE(is_errBC_proc,is_errBC,1,MPI_LOGICAL,MPI_LOR,COMM_global,info)
!!$       if (is_errBC) &
!!$            call mpistop('You need to activate the use of reference values (''r'' flag in param_blocks for inlet/imin)',0)
    endif

  end subroutine bc_define

end module mod_bc

