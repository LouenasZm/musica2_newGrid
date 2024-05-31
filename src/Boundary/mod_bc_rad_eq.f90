!================================================================================
module mod_bc_rad_eq
!================================================================================
  !> Module to apply radial equilibrium approximation at exit BCs
!================================================================================
  use mod_flow ! for: flow variables
  implicit none
  !------------------------------------------------------------------------------
  abstract interface
     subroutine rea_typ(dp)
       use mod_mpi
       use mod_flow
       implicit none
       real(wp), dimension(ny,nz), intent(out) :: dp
     end subroutine rea_typ
  end interface
  !------------------------------------------------------------------------------
  procedure(rea_typ), pointer :: bc_rea_ox
  !------------------------------------------------------------------------------
  ! radial direction & tangential components
  integer :: ir1,ir2 ! starting/ending indices of radial direction
  integer :: na ! index of BC location
  real(wp), dimension(:,:), allocatable, private :: rr    ! radial dir.
  real(wp), dimension(:,:), allocatable, private :: t1,t2 ! tangential comp.
  !------------------------------------------------------------------------------
  ! MPI communicator for procs with exit BC
  integer, private :: COMM_exit
  integer, private :: iproc_exit ! nb of procs
  integer, private :: nproc_exit ! number of procs
  integer, dimension(:), allocatable, private :: rank_exit ! rank of procs
  !------------------------------------------------------------------------------
  ! Choice of reference value (pivot) for exit BC
  character(len=3) :: rea_ref
  integer, private :: iproc_send ! proc having min/max
  integer, private :: ind_m ! either k_min=1 or k_max=nz
  !------------------------------------------------------------------------------

contains

  !==============================================================================
  subroutine init_rea
  !==============================================================================
    !> Init procedure for radial equilibrium (axis is Ox)
    !> * check compatibility with param options
    !> * initialize communicator for exit BC
    !> * determine radial direction
    !> * assign procedures for REA
  !==============================================================================
    use warnstop
    use mod_mpi_part
    use mod_block
    use mod_constant ! for pi
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j,k,n,nc,nbl,ierr
    integer :: radial_dir,transv_dir,proc_in_BCexit
    integer, dimension(:), allocatable :: n_rad
    logical :: is_backpressure,is_backpressure_proc
    logical :: is_exit_imin,is_exit_imax
    logical :: is_exit_jmin,is_exit_jmax
    logical :: is_exit_kmin,is_exit_kmax
    logical :: bloc_has_wall,bloc_has_2walls
    logical, dimension(:), allocatable :: mk ! mask for ordering
    real(wp) :: rmax,r1,r2,t_c,t_s
    real(wp), dimension(:), allocatable :: r_max
    ! not used yet [lavoro in corso]
    !integer :: nbv,nbs
    !integer :: rad_dir_nbv
    !integer, dimension(:), allocatable :: rad_dir
    !real(wp) :: dwall
    !real(wp), dimension(:), allocatable :: dwall_bl
    ! ---------------------------------------------------------------------------
    integer :: COMM_allexit
    integer :: nproc_allexit ! number of procs
    integer :: cpt,cpt1,cpt2
    integer, dimension(:), allocatable :: nth_c,nth_s
    real(wp) :: min_val,max_val
    real(wp), dimension(:), allocatable :: th_c,th_s
    ! ---------------------------------------------------------------------------

    if (iproc==0) then
       print *,repeat('=',70)
       print *,'Init radial equilibrium'
       print *,repeat('=',70)
    endif

    ! Check compatibility with param options
    ! ======================================
    if (.not.is_curv3) call mpistop('REA only for 3D curvilinear solver (is_curv3)',0)

    ! Initialize MPI communicator
    ! ===========================

    ! number of block
    nbl=nob(iproc)
    
    ! init local boolean indicator for 'exit BC'
    is_exit_imin=.false.
    is_exit_imax=.false.
    is_exit_jmin=.false.
    is_exit_jmax=.false.
    is_exit_kmin=.false.
    is_exit_kmax=.false.
    
    ! init local boolean indicator 'has a backpressure outflow BC'
    is_backpressure_proc=.false.

    ! init radial direction (-1 -> not an exit BC)
    radial_dir=-1

    ! init indicator 'block has wall BC'
    bloc_has_wall=.false.
    bloc_has_2walls=.false.
    
    if (BC_face(1,1)%sort==-5) then
       ! exit BC is at imin
       is_backpressure_proc=.true.
       is_exit_imin=.true.
       
       ! should be the only one
       if ( (BC_face(1,2)%sort==-5).or.(BC_face(2,1)%sort==-5).or.&
            (BC_face(2,2)%sort==-5).or.(BC_face(3,1)%sort==-5).or.&
            (BC_face(3,2)%sort==-5) ) call mpistop('Outflow BC -5 on two faces: not implemented',0)

       ! determine radial direction
       ! --------------------------
       ! if wall condition is present in block, this is the radial direction
       radial_dir=0
       ! at jmin/jmax => radial direction is j
       !              => transverse direction is k
       if ((bl(nbl)%BC(3)==0).or.(bl(nbl)%BC(4)==0)) then
          transv_dir=3
          r1=sqrt(yc3(1, 1,1)**2+zc3(1, 1,1)**2)
          r2=sqrt(yc3(1,ny,1)**2+zc3(1,ny,1)**2)
          if (r2<r1) then
             radial_dir=3 ! -j
          else
             radial_dir=4 ! +j
          endif
       endif
       ! at kmin/kmax => radial direction is k
       !              => transverse direction is j
       if ((bl(nbl)%BC(5)==0).or.(bl(nbl)%BC(6)==0)) then
          transv_dir=2
          r1=sqrt(yc3(1,1, 1)**2+zc3(1,1, 1)**2)
          r2=sqrt(yc3(1,1,nz)**2+zc3(1,1,nz)**2)
          if (r2<r1) then
             radial_dir=5 ! -k
          else
             radial_dir=6 ! +k
          endif
       endif
       ! set indicator 'block has wall BC'
       if (radial_dir>0) bloc_has_wall=.true.
       ! if no wall condition (interior block), check connectivity later
    endif

    if (BC_face(1,2)%sort==-5) then
       ! exit BC is at imax
       is_backpressure_proc=.true.
       is_exit_imax=.true.
       
       ! should be the only one
       if ( (BC_face(1,1)%sort==-5).or.(BC_face(2,1)%sort==-5).or.&
            (BC_face(2,2)%sort==-5).or.(BC_face(3,1)%sort==-5).or.&
            (BC_face(3,2)%sort==-5) ) call mpistop('Outflow BC -5 on two faces: not implemented',0)

       ! determine radial direction
       ! --------------------------
       ! if wall condition is present, this is the radial direction
       radial_dir=0
       ! at jmin/jmax => radial direction is j
       !              => transverse direction is k
       if ((bl(nbl)%BC(3)==0).or.(bl(nbl)%BC(4)==0)) then
          transv_dir=3
          r1=sqrt(yc3(nx, 1,1)**2+zc3(nx, 1,1)**2)
          r2=sqrt(yc3(nx,ny,1)**2+zc3(nx,ny,1)**2)
          if (r2<r1) then
             radial_dir=3 ! -j
          else
             radial_dir=4 ! +j
          endif
       endif
       ! at kmin/kmax => radial direction is k
       !              => transverse direction is j
       if ((bl(nbl)%BC(5)==0).or.(bl(nbl)%BC(6)==0)) then
          transv_dir=2
          r1=sqrt(yc3(nx,1, 1)**2+zc3(nx,1, 1)**2)
          r2=sqrt(yc3(nx,1,nz)**2+zc3(nx,1,nz)**2)
          if (r2<r1) then
             radial_dir=5 ! -k
          else
             radial_dir=6 ! +k
          endif
       endif
       ! set indicator 'block has wall BC'
       if (radial_dir>0) bloc_has_wall=.true.
       ! if no wall condition (interior block), check connectivity later
    endif

    if (BC_face(2,1)%sort==-5) then
       ! exit BC is at jmin
       is_backpressure_proc=.true.
       is_exit_jmin=.true.
       
       ! should be the only one
       if ( (BC_face(1,1)%sort==-5).or.(BC_face(1,2)%sort==-5).or.&
            (BC_face(2,2)%sort==-5).or.(BC_face(3,1)%sort==-5).or.&
            (BC_face(3,2)%sort==-5) ) call mpistop('Outflow BC -5 on two faces: not implemented',0)

       ! determine radial direction
       ! --------------------------
       ! if wall condition is present, this is the radial direction
       radial_dir=0
       ! at imin/imax => radial direction is i
       !              => transverse direction is k
       if ((bl(nbl)%BC(1)==0).or.(bl(nbl)%BC(2)==0)) then
          transv_dir=3
          r1=sqrt(yc3( 1,1,1)**2+zc3( 1,1,1)**2)
          r2=sqrt(yc3(nx,1,1)**2+zc3(nx,1,1)**2)
          if (r2<r1) then
             radial_dir=1 ! -i
          else
             radial_dir=2 ! +i
          endif
       endif
       ! at kmin/kmax => radial direction is k
       !              => transverse direction is i
       if ((bl(nbl)%BC(5)==0).or.(bl(nbl)%BC(6)==0)) then
          transv_dir=1
          r1=sqrt(yc3(1,1, 1)**2+zc3(1,1, 1)**2)
          r2=sqrt(yc3(1,1,nz)**2+zc3(1,1,nz)**2)
          if (r2<r1) then
             radial_dir=5 ! -k
          else
             radial_dir=6 ! +k
          endif
       endif
       ! set indicator 'block has wall BC'
       if (radial_dir>0) bloc_has_wall=.true.
       ! if no wall condition (interior block), check connectivity later
    endif

    if (BC_face(2,2)%sort==-5) then
       ! exit BC is at jmax
       is_backpressure_proc=.true.
       is_exit_jmax=.true.
       
       ! should be the only one
       if ( (BC_face(1,1)%sort==-5).or.(BC_face(1,2)%sort==-5).or.&
            (BC_face(2,1)%sort==-5).or.(BC_face(3,1)%sort==-5).or.&
            (BC_face(3,2)%sort==-5) ) call mpistop('Outflow BC -5 on two faces: not implemented',0)

       ! determine radial direction
       ! --------------------------
       ! if wall condition is present, this is the radial direction
       radial_dir=0
       ! at imin/imax => radial direction is i
       !              => transverse direction is k
       if ((bl(nbl)%BC(1)==0).or.(bl(nbl)%BC(2)==0)) then
          transv_dir=3
          r1=sqrt(yc3( 1,ny,1)**2+zc3( 1,ny,1)**2)
          r2=sqrt(yc3(nx,ny,1)**2+zc3(nx,ny,1)**2)
          if (r2<r1) then
             radial_dir=1 ! -i
          else
             radial_dir=2 ! +i
          endif
       endif
       ! at kmin/kmax => radial direction is k
       !              => transverse direction is i
       if ((bl(nbl)%BC(5)==0).or.(bl(nbl)%BC(6)==0)) then
          transv_dir=1
          r1=sqrt(yc3(1,ny, 1)**2+zc3(1,ny, 1)**2)
          r2=sqrt(yc3(1,ny,nz)**2+zc3(1,ny,nz)**2)
          if (r2<r1) then
             radial_dir=5 ! -k
          else
             radial_dir=6 ! +k
          endif
       endif
       ! set indicator 'block has wall BC'
       if (radial_dir>0) bloc_has_wall=.true.
       ! if no wall condition (interior block), check connectivity later
    endif

    if (BC_face(3,1)%sort==-5) then
       ! exit BC is at kmin
       is_backpressure_proc=.true.
       is_exit_kmin=.true.
       
       ! should be the only one
       if ( (BC_face(1,1)%sort==-5).or.(BC_face(1,2)%sort==-5).or.&
            (BC_face(2,1)%sort==-5).or.(BC_face(2,2)%sort==-5).or.&
            (BC_face(3,2)%sort==-5) ) call mpistop('Outflow BC -5 on two faces: not implemented',0)

       ! determine radial direction
       ! --------------------------
       ! if wall condition is present, this is the radial direction
       radial_dir=0
       ! at imin/imax => radial direction is i
       !              => transverse direction is j
       if ((bl(nbl)%BC(1)==0).or.(bl(nbl)%BC(2)==0)) then
          transv_dir=2
          r1=sqrt(yc3( 1,1,1)**2+zc3( 1,1,1)**2)
          r2=sqrt(yc3(nx,1,1)**2+zc3(nx,1,1)**2)
          if (r2<r1) then
             radial_dir=1 ! -i
          else
             radial_dir=2 ! +i
          endif
       endif
       ! at jmin/jmax => radial direction is j
       !              => transverse direction is i
       if ((bl(nbl)%BC(3)==0).or.(bl(nbl)%BC(4)==0)) then
          transv_dir=1
          r1=sqrt(yc3(1, 1,1)**2+zc3(1, 1,1)**2)
          r2=sqrt(yc3(1,ny,1)**2+zc3(1,ny,1)**2)
          if (r2<r1) then
             radial_dir=3 ! -j
          else
             radial_dir=4 ! +j
          endif
       endif
       ! set indicator 'block has wall BC'
       if (radial_dir>0) bloc_has_wall=.true.
       ! if no wall condition (interior block), check connectivity later
    endif

    if (BC_face(3,2)%sort==-5) then
       ! exit BC is at kmax
       is_backpressure_proc=.true.
       is_exit_kmax=.true.
       
       ! should be the only one
       if ( (BC_face(1,1)%sort==-5).or.(BC_face(1,2)%sort==-5).or.&
            (BC_face(2,1)%sort==-5).or.(BC_face(2,2)%sort==-5).or.&
            (BC_face(3,1)%sort==-5) ) call mpistop('Outflow BC -5 on two faces: not implemented',0)

       ! determine radial direction
       ! --------------------------
       ! if wall condition is present, this is the radial direction
       radial_dir=0
       ! at imin/imax => radial direction is i
       !              => transverse direction is j
       if ((bl(nbl)%BC(1)==0).or.(bl(nbl)%BC(2)==0)) then
          transv_dir=2
          r1=sqrt(yc3( 1,1,1)**2+zc3( 1,1,1)**2)
          r2=sqrt(yc3(nx,1,1)**2+zc3(nx,1,1)**2)
          if (r2<r1) then
             radial_dir=1 ! -i
          else
             radial_dir=2 ! +i
          endif
       endif
       ! at jmin/jmax => radial direction is j
       !              => transverse direction is i
       if ((bl(nbl)%BC(3)==0).or.(bl(nbl)%BC(4)==0)) then
          transv_dir=1
          r1=sqrt(yc3(1, 1,1)**2+zc3(1, 1,1)**2)
          r2=sqrt(yc3(1,ny,1)**2+zc3(1,ny,1)**2)
          if (r2<r1) then
             radial_dir=3 ! -j
          else
             radial_dir=4 ! +j
          endif
       endif
       ! set indicator 'block has wall BC'
       if (radial_dir>0) bloc_has_wall=.true.
       ! if no wall condition (interior block), check connectivity later
    endif

    ! global boolean indicator 'has a backpressure outflow BC'
    ! --------------------------------------------------------
    call MPI_ALLREDUCE(is_backpressure_proc,is_backpressure,1,MPI_LOGICAL,MPI_LOR,COMM_global,info)

    if (.not.is_backpressure) call mpistop('REA only implemented for Back Pressure Outflow BC yet',0)

    ! Nota: number of COMM_exit is number of different theta
    ! ======================================================
    
    ! create temporary communicator for all procs having exit
    ! -------------------------------------------------------
    if (is_backpressure_proc) then
       proc_in_BCexit=1
    else
       proc_in_BCexit=0
    endif
    call MPI_COMM_SPLIT(COMM_global,proc_in_BCexit,iproc,COMM_allexit,info)
    ! number of procs in the allexit communicator
    call MPI_COMM_SIZE(COMM_allexit,nproc_allexit,info)

    if (is_backpressure_proc) then

       ! compute theta from cos and sin for corner point
       ! ------------------------------------------------
       if (is_exit_imin) then
          t_c=yc3(1,1,1)/sqrt(yc3(1,1,1)**2+zc3(1,1,1)**2)
          t_c=abs(acos(t_c)*180.0_wp/pi)
          t_s=zc3(1,1,1)/sqrt(yc3(1,1,1)**2+zc3(1,1,1)**2)
          t_s=abs(asin(t_s)*180.0_wp/pi)
       elseif (is_exit_imax) then
          t_c=yc3(nx,1,1)/sqrt(yc3(nx,1,1)**2+zc3(nx,1,1)**2)
          t_c=abs(acos(t_c)*180.0_wp/pi)
          t_s=zc3(nx,1,1)/sqrt(yc3(nx,1,1)**2+zc3(nx,1,1)**2)
          t_s=abs(asin(t_s)*180.0_wp/pi)
       elseif (is_exit_jmin) then
          t_c=yc3(1,1,1)/sqrt(yc3(1,1,1)**2+zc3(1,1,1)**2)
          t_c=abs(acos(t_c)*180.0_wp/pi)
          t_s=zc3(1,1,1)/sqrt(yc3(1,1,1)**2+zc3(1,1,1)**2)
          t_s=abs(asin(t_s)*180.0_wp/pi)
       elseif (is_exit_jmax) then
          t_c=yc3(1,ny,1)/sqrt(yc3(1,1,1)**2+zc3(1,ny,1)**2)
          t_c=abs(acos(t_c)*180.0_wp/pi)
          t_s=zc3(1,ny,1)/sqrt(yc3(1,1,1)**2+zc3(1,ny,1)**2)
          t_s=abs(asin(t_s)*180.0_wp/pi)
       elseif (is_exit_kmin) then
          t_c=yc3(1,1,1)/sqrt(yc3(1,1,1)**2+zc3(1,1,1)**2)
          t_c=abs(acos(t_c)*180.0_wp/pi)
          t_s=zc3(1,1,1)/sqrt(yc3(1,1,1)**2+zc3(1,1,1)**2)
          t_s=abs(asin(t_s)*180.0_wp/pi)
       elseif (is_exit_kmax) then
          t_c=yc3(1,1,nz)/sqrt(yc3(1,1,nz)**2+zc3(1,1,nz)**2)
          t_c=abs(acos(t_c)*180.0_wp/pi)
          t_s=zc3(1,1,nz)/sqrt(yc3(1,1,nz)**2+zc3(1,1,nz)**2)
          t_s=abs(asin(t_s)*180.0_wp/pi)
       endif

       ! share theta (from cos & sin) among procs of COMM_allexit
       ! --------------------------------------------------------
       allocate(th_c(1:nproc_allexit),th_s(1:nproc_allexit))
       call MPI_ALLGATHER(t_c,1,MPI_DOUBLE_PRECISION,th_c,1,MPI_DOUBLE_PRECISION,COMM_allexit,info)
       call MPI_ALLGATHER(t_s,1,MPI_DOUBLE_PRECISION,th_s,1,MPI_DOUBLE_PRECISION,COMM_allexit,info)

       ! count number of different theta
       ! -------------------------------
       ! compute nearest integer values
       allocate(nth_c(1:nproc_allexit),nth_s(1:nproc_allexit))
       do n=1,nproc_allexit
          nth_c(n)=nint(th_c(n))
          nth_s(n)=nint(th_s(n))
          !print *,'val',n,nth_c(n),nth_s(n)
       enddo

       ! count number of different theta (from cos)
       cpt1=0
       min_val=minval(nth_c)-1
       max_val=maxval(nth_c)
       do while (min_val<max_val)
          cpt1=cpt1+1
          min_val=minval(nth_c,mask=nth_c>min_val)
          th_c(cpt1)=min_val
       enddo
       !print *,cpt1,':',nint(th_c(1:cpt1))

       ! count number of different theta (from sin)
       cpt2=0
       min_val=minval(nth_s)-1
       max_val=maxval(nth_s)
       do while (min_val<max_val)
          cpt2=cpt2+1
          min_val=minval(nth_s,mask=nth_s>min_val)
          th_s(cpt2)=min_val
       enddo
       !print *,cpt2,':',nint(th_s(1:cpt2))

       ! global counter
       cpt=max(cpt1,cpt2)
       print *,iproc,'number',cpt

    endif
    
    ! create communicator of all processors with exit condition having same theta
    ! ---------------------------------------------------------------------------
    if (is_backpressure_proc) then
       !!proc_in_BCexit=coord(transv_dir)+1
       do n=1,cpt
          if ((nint(t_c)==nint(th_c(n))).and.(nint(t_s)==nint(th_s(n)))) then
             proc_in_BCexit=n
          endif
       enddo
    else
       proc_in_BCexit=0
    endif
    !if (is_backpressure_proc) print *,iproc,'my color is',proc_in_BCexit,nint(t_c),nint(t_s)

    call MPI_COMM_SPLIT(COMM_global,proc_in_BCexit,iproc,COMM_exit,info)

    ! number of procs in the exit communicator
    call MPI_COMM_SIZE(COMM_exit,nproc_exit,info)

    ! No proc in COMM_exit
    call MPI_COMM_RANK(COMM_exit,iproc_exit,info)
    
    !print *,'check init REA',is_backpressure
    !if (is_backpressure_proc) then
    !   print *,'proc',iproc,':',is_backpressure_proc,bloc_has_wall,transv_dir,radial_dir,proc_in_BCexit
    !   print *,'proc',iproc,':',is_backpressure_proc,proc_in_BCexit,nproc_exit,iproc_exit
    !endif
    !if (is_backpressure_proc) print *,iproc,proc_in_BCexit,nproc_exit,iproc_exit,transv_dir
        
    ! correct radial direction
    ! ------------------------    
!!$       ! share radial dir. with all procs in the exit communicator
!!$       allocate(rad_dir(nproc_exit))
!!$       call MPI_ALLGATHER(radial_dir,1,MPI_INTEGER,rad_dir,1,MPI_INTEGER,COMM_exit,info)
   

!!$    ! if not determined yet (because no wall condition in interior block), check connectivity
!!$    if (radial_dir==0) then
!!$       ! check neighboring block
!!$       do n=1,6
!!$          ! neighboring block
!!$          nbv=neighbor_bl(n)
!!$          ! its radial direction
!!$          rad_dir_nbv=rad_dir(iproc_leader(nbv))
!!$          ! if -1 neighbouring block is not an exit block
!!$          ! if  0 radial dir. of neighbouring block is undtermined
!!$          if ((rad_dir_nbv==-1).or.(rad_dir_nbv==0)) cycle
!!$          ! determine radial direction taking into account eventual swap
!!$          if (rad_dir_nbv==1) then
!!$             if (is_swapij2_bl(n)) then
!!$                radial_dir=2
!!$             elseif (is_swapik2_bl(n)) then
!!$                radial_dir=3
!!$             else
!!$                radial_dir=1
!!$             endif
!!$          endif
!!$          if (rad_dir_nbv==2) then
!!$             if (is_swapij2_bl(n)) then
!!$                radial_dir=1
!!$             elseif (is_swapjk2_bl(n)) then
!!$                radial_dir=3
!!$             else
!!$                radial_dir=2
!!$             endif
!!$          endif
!!$          if (rad_dir_nbv==3) then
!!$             if (is_swapik2_bl(n)) then
!!$                radial_dir=1
!!$             elseif (is_swapjk2_bl(n)) then
!!$                radial_dir=2
!!$             else
!!$                radial_dir=3
!!$             endif
!!$          endif
!!$       enddo
!!$    endif


!!$    if (proc_in_BCexit==1) print *,iproc,'radial_dir,is_bc_wall(2,:)',radial_dir,is_bc_wall(2,1),is_bc_wall(2,2)
!!$    call mpistop('nnnnnnnnn',0)
    
    ! For procs belonging to exit BC
    ! ==============================
    if (is_backpressure_proc) then
       
       ! 1/ compute radial direction & tangential components
       !    ------------------------------------------------

       ! radial direction is i
       ! ---------------------
       if ((radial_dir==1).or.(radial_dir==2)) then
          ! starting and ending indices for radial dir.
          ir1=1
          ir2=nx
          ! -i: add 1 pt before if imin is not a wall
          if ((radial_dir==1).and.(.not.is_bc_wall(1,1))) ir1=ir1-1
          ! +i: add 1 pt after if imax is not a wall
          if ((radial_dir==2).and.(.not.is_bc_wall(1,2))) ir2=ir2+1

          ! for jmin/jmax BC, k is azimuthal direction
          if (is_exit_jmin) na=1
          if (is_exit_jmax) na=ngy
          if ((is_exit_jmin).or.(is_exit_jmax)) then
             ! allocate
             allocate(rr(ir1:ir2,nz),t1(ir1:ir2,nz),t2(ir1:ir2,nz))
             ! use global coordinate for points at nr+1 
             do k=1,nz
                do i=ir1,ir2
                   r2=ygc3(i+coord(1)*nx,na,k+coord(3)*nz)**2 &
                     +zgc3(i+coord(1)*nx,na,k+coord(3)*nz)**2
                   rr(i,k)=sqrt(r2)
                   t1(i,k)= zgc3(i+coord(1)*nx,na,k+coord(3)*nz)/rr(i,k)
                   t2(i,k)=-ygc3(i+coord(1)*nx,na,k+coord(3)*nz)/rr(i,k)
                enddo
             enddo
          endif

          ! for kmin/kmax BC, j is azimuthal direction
          if (is_exit_kmin) na=1
          if (is_exit_kmax) na=ngz
          if ((is_exit_kmin).or.(is_exit_kmax)) then
             ! allocate
             allocate(rr(ir1:ir2,ny),t1(ir1:ir2,ny),t2(ir1:ir2,ny))
             ! use global coordinate for points at nr+1 
             do j=1,ny
                do i=ir1,ir2
                   r2=ygc3(i+coord(1)*nx,j+coord(2)*ny,na)**2 &
                     +zgc3(i+coord(1)*nx,j+coord(2)*ny,na)**2
                   rr(i,j)=sqrt(r2)
                   t1(i,j)= zgc3(i+coord(1)*nx,j+coord(2)*ny,na)/rr(i,j)
                   t2(i,j)=-ygc3(i+coord(1)*nx,j+coord(2)*ny,na)/rr(i,j)
                enddo
             enddo
          endif
       endif

       ! radial direction is j
       ! ---------------------
       if ((radial_dir==3).or.(radial_dir==4)) then
          ! starting and ending indices for radial dir.
          ir1=1
          ir2=ny
          ! -j: add 1 pt before if jmin is not a wall
          if ((radial_dir==3).and.(.not.is_bc_wall(2,1))) then
             if (proc_in_BCexit==1) print *,'radial_dir,is_bc_wall(2,1)',radial_dir,is_bc_wall(2,1)
             ir1=ir1-1
          endif
          ! +j: add 1 pt after if jmax is not a wall
          if ((radial_dir==4).and.(.not.is_bc_wall(2,2))) then
             if (proc_in_BCexit==1) print *,'radial_dir,is_bc_wall(2,2)',radial_dir,is_bc_wall(2,2)
             ir2=ir2+1
          endif

          ! for imin/imax BC, k is azimuthal direction
          if (is_exit_imin) na=1
          if (is_exit_imax) na=ngx
          if ((is_exit_imin).or.(is_exit_imax)) then
             ! allocate
             allocate(rr(ir1:ir2,nz),t1(ir1:ir2,nz),t2(ir1:ir2,nz))
             ! use global coordinate for points at nr+1 
             do k=1,nz
                do j=ir1,ir2
                   r2=ygc3(na,j+coord(2)*ny,k+coord(3)*nz)**2 &
                      +zgc3(na,j+coord(2)*ny,k+coord(3)*nz)**2
                   rr(j,k)=sqrt(r2)
                   t1(j,k)= zgc3(na,j+coord(2)*ny,k+coord(3)*nz)/rr(j,k)
                   t2(j,k)=-ygc3(na,j+coord(2)*ny,k+coord(3)*nz)/rr(j,k)
                enddo
             enddo
          endif

          ! for kmin/kmax BC, i is azimuthal direction
          if (is_exit_kmin) na=1
          if (is_exit_kmax) na=ngz
          if ((is_exit_kmin).or.(is_exit_kmax)) then
             ! allocate
             allocate(rr(ir1:ir2,nx),t1(ir1:ir2,nx),t2(ir1:ir2,nx))
             ! use global coordinate for points at nr+1 
             do i=1,nx
                do j=ir1,ir2
                   r2=ygc3(i+coord(1)*nx,j+coord(2)*ny,na)**2 &
                      +zgc3(i+coord(1)*nx,j+coord(2)*ny,na)**2
                   rr(j,i)=sqrt(r2)
                   t1(j,i)= zgc3(i+coord(1)*nx,j+coord(2)*ny,na)/rr(j,i)
                   t2(j,i)=-ygc3(i+coord(1)*nx,j+coord(2)*ny,na)/rr(j,i)
                enddo
             enddo
          endif
       endif

       ! radial direction is k
       ! ---------------------
       if ((radial_dir==5).or.(radial_dir==6)) then
          ! starting and ending indices for radial dir.
          ir1=1
          ir2=nz
          ! -k: add 1 pt before if jmin is not a wall
          if ((radial_dir==5).and.(.not.is_bc_wall(3,1))) ir1=ir1-1
          ! +k: add 1 pt after if jmax is not a wall
          if ((radial_dir==6).and.(.not.is_bc_wall(3,2))) ir2=ir2+1

          !print *,iproc,radial_dir,ir1,ir2,is_exit_imin,is_exit_imax
          
          ! for imin/imax BC, j is azimuthal direction
          if (is_exit_imin) na=1
          if (is_exit_imax) na=ngx
          if ((is_exit_imin).or.(is_exit_imax)) then
             ! allocate
             allocate(rr(ir1:ir2,ny),t1(ir1:ir2,ny),t2(ir1:ir2,ny))
             ! use global coordinate for points at nr+1 
             do j=1,ny
                do k=ir1,ir2
                   r2=ygc3(na,j+coord(2)*ny,k+coord(3)*nz)**2 &
                     +zgc3(na,j+coord(2)*ny,k+coord(3)*nz)**2
                   rr(k,j)=sqrt(r2)
                   t1(k,j)= zgc3(na,j+coord(2)*ny,k+coord(3)*nz)/rr(k,j)
                   t2(k,j)=-ygc3(na,j+coord(2)*ny,k+coord(3)*nz)/rr(k,j)
                enddo
             enddo
          endif

          ! for jmin/jmax BC, i is azimuthal direction
          if (is_exit_jmin) na=1
          if (is_exit_jmax) na=ngy
          if ((is_exit_jmin).or.(is_exit_jmax)) then
             ! allocate
             allocate(rr(ir1:ir2,nx),t1(ir1:ir2,nx),t2(ir1:ir2,nx))
             ! use global coordinate for points at nr+1 
             do i=1,nx
                do k=ir1,ir2
                   r2=ygc3(i+coord(1)*nx,na,k+coord(3)*nz)**2 &
                      +zgc3(i+coord(1)*nx,na,k+coord(3)*nz)**2
                   rr(k,i)=sqrt(r2)
                   t1(k,i)= zgc3(i+coord(1)*nx,na,k+coord(3)*nz)/rr(k,i)
                   t2(k,i)=-ygc3(i+coord(1)*nx,na,k+coord(3)*nz)/rr(k,i)
                enddo
             enddo
          endif

          ! for kmin/kmax BC, i is azimuthal direction
          if (is_exit_kmin) na=1
          if (is_exit_kmax) na=ngz
          if ((is_exit_kmin).or.(is_exit_kmax)) then
             ! allocate
             allocate(rr(ir1:ir2,nx),t1(ir1:ir2,nx),t2(ir1:ir2,nx))
             ! use global coordinate for points at nr+1 
             do i=1,nx
                do j=ir1,ir2
                   r2=ygc3(i+coord(1)*nx,j+coord(2)*ny,na)**2 &
                     +zgc3(i+coord(1)*nx,j+coord(2)*ny,na)**2
                   rr(j,i)=sqrt(r2)
                   t1(j,i)= zgc3(i+coord(1)*nx,j+coord(2)*ny,na)/rr(j,i)
                   t2(j,i)=-ygc3(i+coord(1)*nx,j+coord(2)*ny,na)/rr(j,i)
                enddo
             enddo
          endif
       endif

    !if (proc_in_BCexit==coord(transv_dir)+1) then
    
       ! 2/ Ascending re-ordering of procs in COMM_exit
       !    -------------------------------------------
       ! for each sub-communicator
       do nc=1,cpt
          if (proc_in_BCexit==nc) then
             
             ! determine maximum radius for each domain of exit BC
             select case (radial_dir)
             case (1)
                if (is_exit_jmin) rmax=sqrt(yc3( 1, 1, 1)**2+zc3( 1, 1, 1)**2)
                if (is_exit_jmax) rmax=sqrt(yc3( 1,ny, 1)**2+zc3( 1,ny, 1)**2)
                if (is_exit_kmin) rmax=sqrt(yc3( 1, 1, 1)**2+zc3( 1, 1, 1)**2)
                if (is_exit_kmax) rmax=sqrt(yc3( 1, 1,nz)**2+zc3( 1, 1,nz)**2)
             case (2)
                if (is_exit_jmin) rmax=sqrt(yc3(nx, 1, 1)**2+zc3(nx, 1, 1)**2)
                if (is_exit_jmax) rmax=sqrt(yc3(nx,ny, 1)**2+zc3(nx,ny, 1)**2)
                if (is_exit_kmin) rmax=sqrt(yc3(nx, 1, 1)**2+zc3(nx, 1, 1)**2)
                if (is_exit_kmax) rmax=sqrt(yc3(nx, 1,nz)**2+zc3(nx, 1,nz)**2)
             case (3)
                if (is_exit_imin) rmax=sqrt(yc3( 1, 1, 1)**2+zc3( 1, 1, 1)**2)
                if (is_exit_imax) rmax=sqrt(yc3(nx, 1, 1)**2+zc3(nx, 1, 1)**2)
                if (is_exit_kmin) rmax=sqrt(yc3( 1, 1, 1)**2+zc3( 1, 1, 1)**2)
                if (is_exit_kmax) rmax=sqrt(yc3( 1, 1,nz)**2+zc3( 1, 1,nz)**2)
             case (4)
                if (is_exit_imin) rmax=sqrt(yc3( 1,ny, 1)**2+zc3( 1,ny, 1)**2)
                if (is_exit_imax) rmax=sqrt(yc3(nx,ny, 1)**2+zc3(nx,ny, 1)**2)
                if (is_exit_kmin) rmax=sqrt(yc3( 1,ny, 1)**2+zc3( 1,ny, 1)**2)
                if (is_exit_kmax) rmax=sqrt(yc3( 1,ny,nz)**2+zc3( 1,ny,nz)**2)
             case (5)
                if (is_exit_imin) rmax=sqrt(yc3( 1, 1, 1)**2+zc3( 1, 1, 1)**2)
                if (is_exit_imax) rmax=sqrt(yc3(nx, 1, 1)**2+zc3(nx, 1, 1)**2)
                if (is_exit_jmin) rmax=sqrt(yc3( 1, 1, 1)**2+zc3( 1, 1, 1)**2)
                if (is_exit_jmax) rmax=sqrt(yc3( 1,ny, 1)**2+zc3( 1,ny, 1)**2)
             case (6)
                if (is_exit_imin) rmax=sqrt(yc3( 1, 1,nz)**2+zc3( 1, 1,nz)**2)
                if (is_exit_imax) rmax=sqrt(yc3(nx, 1,nz)**2+zc3(nx, 1,nz)**2)
                if (is_exit_jmin) rmax=sqrt(yc3( 1, 1,nz)**2+zc3( 1, 1,nz)**2)
                if (is_exit_jmax) rmax=sqrt(yc3( 1,ny,nz)**2+zc3( 1,ny,nz)**2)
             end select

             ! share rmax among procs of COMM_exit
             allocate(r_max(1:nproc_exit))
             call MPI_ALLGATHER(rmax,1,MPI_DOUBLE_PRECISION,r_max,1,MPI_DOUBLE_PRECISION,COMM_exit,info)

             ! create mask for sorting rmax
             allocate(mk(1:nproc_exit),rank_exit(0:nproc_exit-1))
             mk=.true.

             ! ordering rmax for exit procs
             do n=1,nproc_exit
                k=minloc(r_max,1,mk) ! index of minimum
                rank_exit(n-1)=k-1 ! numbering of procs starts at 0
                mk(k)=.false. ! suppress rank in mask
             enddo

             !print *,'proc',iproc,'rad dir:',rad_dir
             !print *,'bl',nbl,'proc',iproc,'r_max:',r_max
             !print *,'bl',nbl,'proc',iproc,iproc_exit,'rank:',rank_exit(iproc_exit)

          endif
       enddo
    endif
    
!!$    ! reordering rank in exit communicator
!!$    ! ------------------------
!!$
!!$    ! compute distance from axis Ox to wall BC
!!$    dwall=1.0e9_wp
!!$    if (bloc_has_wall) then
!!$       if (radial_dir==1) then
!!$          if (bl(nbl)%BC(1)==0) then
!!$             if (is_exit_jmin) dwall=xgc3(1,1,1)
!!$             if (is_exit_jmax) dwall=xgc3(1,ngy,1)
!!$             if (is_exit_kmin) dwall=xgc3(1,1,1)
!!$             if (is_exit_kmax) dwall=xgc3(1,1,ngz)
!!$          elseif  (bl(nbl)%BC(2)==0) then
!!$             if (is_exit_jmin) dwall=xgc3(ngx,1,1)
!!$             if (is_exit_jmax) dwall=xgc3(ngx,ngy,1)
!!$             if (is_exit_kmin) dwall=xgc3(ngx,1,1)
!!$             if (is_exit_kmax) dwall=xgc3(ngx,1,ngz)
!!$          else
!!$             print *,'pb in algorithm'
!!$          endif
!!$       endif
!!$       if (radial_dir==2) then
!!$          if (bl(nbl)%BC(3)==0) then
!!$             if (is_exit_imin) dwall=ygc3(1,1,1)
!!$             if (is_exit_imax) dwall=ygc3(ngx,1,1)
!!$             if (is_exit_kmin) dwall=ygc3(1,1,1)
!!$             if (is_exit_kmax) dwall=ygc3(1,1,ngz)
!!$          elseif  (bl(nbl)%BC(4)==0) then
!!$             if (is_exit_imin) dwall=ygc3(1,ngy,1)
!!$             if (is_exit_imax) dwall=ygc3(ngx,ngy,1)
!!$             if (is_exit_kmin) dwall=ygc3(1,ngy,1)
!!$             if (is_exit_kmax) dwall=ygc3(1,ngy,ngz)
!!$          else
!!$             print *,'pb in algorithm'
!!$          endif
!!$       endif
!!$       if (radial_dir==3) then
!!$          if (bl(nbl)%BC(5)==0) then
!!$             if (is_exit_imin) dwall=zgc3(1,1,1)
!!$             if (is_exit_imax) dwall=zgc3(ngx,1,1)
!!$             if (is_exit_jmin) dwall=zgc3(1,1,1)
!!$             if (is_exit_jmax) dwall=zgc3(1,ngy,1)
!!$          elseif  (bl(nbl)%BC(6)==0) then
!!$             if (is_exit_imin) dwall=zgc3(1,1,ngz)
!!$             if (is_exit_imax) dwall=zgc3(ngx,1,ngz)
!!$             if (is_exit_jmin) dwall=zgc3(1,1,ngz)
!!$             if (is_exit_jmax) dwall=zgc3(1,ngy,ngz)
!!$          else
!!$             print *,'pb in algorithm'
!!$          endif
!!$       endif
!!$    endif
!!$    ! determine starting block
!!$    ! [has a wall BC and closest to axis]
!!$    if (iproc_leader(nbl)) then
!!$       allocate(dwall_bl(nbloc))
!!$       call MPI_ALLGATHER(dwall,1,MPI_DOUBLE_PRECISION,dwall_bl,1,MPI_DOUBLE_PRECISION,COMM_interblock,info)
!!$       if (iproc==0) print *,'dwall_bl',dwall_bl
!!$       ! number of starting block
!!$       nbs=minloc(dwall_bl,1)
!!$    endif
!!$    ! broadcast the number of the starting block
!!$    call MPI_BCAST(nbs,1,MPI_INTEGER,iproc_leader(1),COMM_global,info)

!!$    ! renumbering is done by starting block
!!$    if (iproc_leader(nbs)) then
!!$       allocate(ip_exit(nproc_exit))
!!$
!!$       ! start from
!!$       select case (radial_dir)
!!$       case (1)
!!$          do n=1,ndomx
!!$             ip_exit(n)=n-1
!!$          enddo
!!$       case (2)
!!$          do n=1,ndomy
!!$             ip_exit(n)=n-1
!!$          enddo
!!$       case (3)
!!$          do n=1,ndomx
!!$             ip_exit(n)=n-1
!!$          enddo
!!$       end select
!!$    endif
   
    ! Assign procedures
    ! =================
!!$    if (is_backpressure_proc) then
!!$       select case (radial_dir)
!!$       case (1)
!!$          if ((is_exit_jmin).or.(is_exit_jmax)) bc_rea_ox=>bc_rea_ox_ip
!!$          if ((is_exit_kmin).or.(is_exit_kmax)) 
!!$       case (2)
!!$          if ((is_exit_jmin).or.(is_exit_jmax)) 
!!$          if ((is_exit_kmin).or.(is_exit_kmax)) 
!!$       case (3)
!!$          if ((is_exit_imin).or.(is_exit_imax)) 
!!$          if ((is_exit_kmin).or.(is_exit_kmax)) 
!!$       case (4)
!!$          if ((is_exit_imin).or.(is_exit_imax)) 
!!$          if ((is_exit_kmin).or.(is_exit_kmax)) 
!!$       case (5)
!!$          if ((is_exit_imin).or.(is_exit_imax)) 
!!$          if ((is_exit_jmin).or.(is_exit_jmax)) 
!!$       case (6)
!!$          if ((is_exit_imin).or.(is_exit_imax)) 
!!$          if ((is_exit_jmin).or.(is_exit_jmax)) 
!!$       end select
!!$    endif

    ierr=0
    if (radial_dir>0) then   
       if ((is_exit_imax).and.(radial_dir==4)) then
          print *,'exit imax & radial direction +j',iproc,'bl',nob(iproc)
          bc_rea_ox=>bc_rea_ox_jp
       elseif ((is_exit_imax).and.(radial_dir==6)) then
          print *,'exit imax & radial direction +k',iproc,'bl',nob(iproc)
          bc_rea_ox=>bc_rea_ox_kp
       else
          print *,'radial dir',radial_dir      
          if (is_exit_imin) then
             ierr=1
          elseif (is_exit_jmin) then
             ierr=2
          elseif (is_exit_jmax) then
             ierr=3
          elseif (is_exit_kmin) then
             ierr=4
          elseif (is_exit_kmax) then
             ierr=5
          else
             ierr=6
          endif
       endif
    endif

    select case (ierr)
    case (1)
       call mpistop('REA for imin exit not implemented yet',0)
    case (2)
       call mpistop('REA for jmin exit not implemented yet',0)
    case (3)
       call mpistop('REA for jmax exit not implemented yet',0)
    case (4)
       call mpistop('REA for kmin exit not implemented yet',0)
    case (5)
       call mpistop('REA for kmax exit not implemented yet',0)
    case (6)
       call mpistop('REA for this direction not implemented yet',0)
    end select

    ! Check that all blocks involved in BC_exit has the same dimension
    ! ----------------------------------------------------------------
    ! Nota: important to avoid blocking MPI_ALLREDUCE for global mean calculation in bc_rea_ox_*
    if (is_backpressure_proc) then
       if ((radial_dir==1).or.(radial_dir==2)) then !+/-i
          ! share dim in radial dir. with all procs in the exit communicator
          allocate(n_rad(nproc_exit))
          call MPI_ALLGATHER(nx,1,MPI_INTEGER,n_rad,1,MPI_INTEGER,COMM_exit,info)
          do i=2,nproc_exit
             if (n_rad(1).ne.n_rad(i)) ierr=7
          enddo 
       elseif ((radial_dir==3).or.(radial_dir==4)) then !+/-j
          ! share dim in radial dir. with all procs in the exit communicator
          allocate(n_rad(nproc_exit))
          call MPI_ALLGATHER(ny,1,MPI_INTEGER,n_rad,1,MPI_INTEGER,COMM_exit,info)
          do i=2,nproc_exit
             if (n_rad(1).ne.n_rad(i)) ierr=7
          enddo 
       elseif ((radial_dir==5).or.(radial_dir==6)) then !+/-k
          ! share dim in radial dir. with all procs in the exit communicator
          allocate(n_rad(nproc_exit))
          call MPI_ALLGATHER(nz,1,MPI_INTEGER,n_rad,1,MPI_INTEGER,COMM_exit,info)
          do i=2,nproc_exit
             if (n_rad(1).ne.n_rad(i)) ierr=7
          enddo 
       endif
    endif

    if (ierr==7) then
       print *,'problem in dimensions for radial direction',n_rad
       call mpistop('change partitioning for REA',0)
    endif
        
    ! Choice of reference value (pivot) for exit BC
    ! =============================================

    ! Choice --> TO BE CHANGED --> in param.ini
    !rea_ref='min'
    rea_ref='max'
    !rea_ref='moy'
    
    if (is_backpressure_proc) then
       if (rea_ref=='min') then
          ! proc having min/max
          do n=0,nproc_exit-1
             if (rank_exit(n)==0) iproc_send=n
          enddo          
          ind_m=1
       elseif (rea_ref=='max') then
          ! proc having min/max
          do n=0,nproc_exit-1
             if (rank_exit(n)==nproc_exit-1) iproc_send=n
          enddo          
          if (radial_dir==4) ind_m=ny
          if (radial_dir==6) ind_m=nz
       endif
       !print *,iproc,iproc_exit,iproc_send,r_max(iproc_exit),rank_exit(iproc_exit)
    endif

    !call mpistop('check in init_rea',0)
    
  end subroutine init_rea

  !==============================================================================
  subroutine bc_rea_ox_jp(dpr)
  !==============================================================================
    !> Compute pressure distribution consistent with radial equilibrium
    !> * apply Radial Equilibrium Approximation (REA): dp/dr=ro*(vth^2)/r
    !> * compute pressure correction to be added to reference pressure (p_exit)
    !> * face imax; rotational axis ox (along i); radial direction j
    !> * only for 3D curvilinear (is_curv3)
  !==============================================================================
    use mod_mpi
    use mod_constant
    implicit none
    ! ---------------------------------------------------------------------------
    ! pressure correction
    real(wp), dimension(ny,nz), intent(out) :: dpr 
    ! ---------------------------------------------------------------------------
    integer :: j,k,ip
    real(wp) :: vth,dprm
    real(wp), dimension(ir1:ir2) :: dpdr
    ! MPI comm
    real(wp), dimension(nz) :: dpr_max
    real(wp), dimension(nz,0:nproc_exit-1) :: dpr_max_ip
    ! ---------------------------------------------------------------------------

    ! Radial Equilibrium Approximation
    ! ================================    
    do k=1,nz

!!$       ! compute dp/dr=ro*(vth^2)/r
!!$       ! --------------------------
!!$       do j=1,nyp
!!$          rr2=yc3(nx,j,k)**2+zc3(nx,j,k)**2
!!$          if ((k==1).and.(j==ny+1)) print *,j,rr2,yc3(nx,j,k),zc3(nx,j,k)
!!$          rr(j)=sqrt(rr2)
!!$          ty= zc3(nx,j,k)/rr(j)
!!$          tz=-yc3(nx,j,k)/rr(j)
!!$          vth=(rhov_n(nx,j,k)*ty+rhow_n(nx,j,k)*tz)/rho_n(nx,j,k)
!!$          dpdr(j)=rho_n(nx,j,k)*vth**2/rr(j)
!!$       enddo

       ! compute dp/dr=ro*(vth^2)/r
       ! --------------------------
       do j=ir1,ir2
          vth=(rhov_n(nx,j,k)*t1(j,k)+rhow_n(nx,j,k)*t2(j,k))/rho_n(nx,j,k)
          dpdr(j)=rho_n(nx,j,k)*vth**2/rr(j,k)
       enddo

       ! integrate with trapezoidal rule
       ! -------------------------------
!!$       !dpr(1,k)=0.5*(dpdr(2)+dpdr(1))*(rr(2)-rr(1))
!!$       !dpr(1,k)=0.5*(dpdr(2)+dpdr(1))*(rr(2)-rr(1))
!!$       dpr(1,k)=0.0_wp
!!$       do j=2,ny-1
!!$          dpr(j+1,k)=dpr(j-1+1,k)
!!$          dpr(j+1,k)=dpr(j+1,k)+0.5*(dpdr(j+1)+dpdr(j))*(rr(j+1)-rr(j))
!!$       enddo
!!$       !dpr(ny,k)=dpr(ny-1,k)+dpdr(ny)*(rr(ny)-rr(ny-1))

       dpr(1,k)=0.0_wp
       do j=2,ny
          dpr(j,k)=dpr(j-1,k)
          dpr(j,k)=dpr(j,k)+0.5*(dpdr(j)+dpdr(j-1))*(rr(j,k)-rr(j-1,k))
       enddo

!!$       ! store integral for ny+1
!!$       ! -----------------------
!!$       if (.not.is_bc_wall(2,1)) then
!!$          !print *,'not a bord',iproc,iproc_exit
!!$          dpr_max(k)=dpr(ny,k)+0.5*(dpdr(ny+1)+dpdr(ny))*(rr(ny+1,k)-rr(ny,k))
!!$       endif
       
       ! store integral for ny
       ! ---------------------
       dpr_max(k)=dpr(ir2-1,k)+0.5*(dpdr(ir2)+dpdr(ir2-1))*(rr(ir2,k)-rr(ir2-1,k))

       !dpr_max(k)=dpr(ir2,k)

    enddo

    ! Communication of dpr_max for COMM_exit
    ! ======================================

    ! gather all dpr_max
    ! ------------------
    call MPI_ALLGATHER(dpr_max,nz,MPI_DOUBLE_PRECISION,dpr_max_ip,nz,MPI_DOUBLE_PRECISION,COMM_exit,info)

    !print *,'dpr_max',dpr_max(1)
    !print *,'dpr_max_ip',dpr_max_ip(1,0),dpr_max_ip(1,1)
    
    ! add dpr_max for procs with lower r_max
    ! --------------------------------------
    do k=1,nz
       do ip=0,nproc_exit-1
          !if (k==1) print *,ip,rank_exit(ip),iproc_exit,rank_exit(iproc_exit),iproc
          if (rank_exit(ip)<rank_exit(iproc_exit)) then
             !print *,'add',ip,dpr_max_ip(k,rank_exit(ip))
             dpr(:,k)=dpr(:,k)+dpr_max_ip(k,rank_exit(ip))
          endif
       enddo
    enddo

    !dpr=dpr!+p_ref
    !if (iproc==2) print *,dpr(1:10,1)

    ! Subtract mean value
    ! ===================
    do k=1,nz

       ! local mean
       ! ----------
       dprm=0.0_wp
       do j=1,ny
          dprm=dprm+dpr(j,k)
       enddo
       
       ! global mean
       ! -----------
       call MPI_ALLREDUCE(MPI_IN_PLACE,dprm,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_exit,info)
       dprm=dprm/dble(ny-1)/nproc_exit
       !if (k==1) print *,dprm
    
       ! subtract mean value
       ! -------------------
       !dpr(:,k)=dpr(:,k)*p_ref/dprm!-dprm
       dpr(:,k)=dpr(:,k)-dprm
       
    enddo

!!$    ! Check mean value
!!$    ! ================
!!$    do k=1,nz
!!$       ! local mean
!!$       ! ----------
!!$       dprm=0.0_wp
!!$       do j=1,ny
!!$          dprm=dprm+dpr(j,k)
!!$       enddo       
!!$       ! global mean
!!$       ! -----------
!!$       call MPI_ALLREDUCE(MPI_IN_PLACE,dprm,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_exit,info)
!!$       dprm=dprm/dble(ny-1)/nproc_exit
!!$       if (k==1) print *,dprm/p_ref          
!!$    enddo

    !dpr=dpr-(1.-0.999196914621331)*p_ref
    ! annulus
    !!dpr=dpr-(1.-0.9999)*p_ref
    !if (iproc==2) print *,'subtract mean'
    !if (iproc==2) print *,dpr(1:10,1)
    
  end subroutine bc_rea_ox_jp

  !==============================================================================
  subroutine bc_rea_ox_kp(dpr)
  !==============================================================================
    !> Compute pressure distribution consistent with radial equilibrium
    !> * apply Radial Equilibrium Approximation (REA): dp/dr=ro*(vth^2)/r
    !> * compute pressure correction to be added to reference pressure (p_exit)
    !> * face imax; rotational axis ox (along i); radial direction k
    !> * only for 3D curvilinear (is_curv3)
  !==============================================================================
    use mod_mpi
    use mod_constant
    implicit none
    ! ---------------------------------------------------------------------------
    ! pressure correction
    real(wp), dimension(ny,nz), intent(out) :: dpr 
    ! ---------------------------------------------------------------------------
    integer :: j,k,ip
    real(wp) :: vth,dprm
    real(wp), dimension(ir1:ir2) :: dpdr
    ! MPI comm
    real(wp), dimension(ny) :: dpr_max,dprmo
    real(wp), dimension(ny,0:nproc_exit-1) :: dpr_max_ip
    ! ---------------------------------------------------------------------------

    ! Radial Equilibrium Approximation
    ! ================================    
    do j=1,ny

       ! compute dp/dr=ro*(vth^2)/r
       ! --------------------------
       do k=ir1,ir2
          vth=(rhov_n(nx,j,k)*t1(k,j)+rhow_n(nx,j,k)*t2(k,j))/rho_n(nx,j,k)
          dpdr(k)=rho_n(nx,j,k)*vth**2/rr(k,j)
       enddo

       ! integrate with trapezoidal rule
       ! -------------------------------
       dpr(j,1)=0.0_wp
       do k=2,nz
          dpr(j,k)=dpr(j,k-1)
          dpr(j,k)=dpr(j,k)+0.5*(dpdr(k)+dpdr(k-1))*(rr(k,j)-rr(k-1,j))
       enddo
       
       ! store integral for ny
       ! ---------------------
       dpr_max(j)=dpr(j,ir2-1)+0.5*(dpdr(ir2)+dpdr(ir2-1))*(rr(ir2,j)-rr(ir2-1,j))

    enddo

    ! Communication of dpr_max for COMM_exit
    ! ======================================

    ! gather all dpr_max
    ! ------------------
    call MPI_ALLGATHER(dpr_max,ny,MPI_DOUBLE_PRECISION,dpr_max_ip,ny,MPI_DOUBLE_PRECISION,COMM_exit,info)

    ! add dpr_max for procs with lower r_max
    ! --------------------------------------
    do j=1,ny
       do ip=0,nproc_exit-1
          if (rank_exit(ip)<rank_exit(iproc_exit)) then
             dpr(j,:)=dpr(j,:)+dpr_max_ip(j,rank_exit(ip))
          endif
       enddo
    enddo

    if (rea_ref=='moy') then
       
       ! Subtract mean value
       ! ===================
       do j=1,ny
          ! local mean
          ! ----------
          dprm=0.0_wp
          do k=1,nz
             dprm=dprm+dpr(j,k)
          enddo

          ! global mean
          ! -----------
          call MPI_ALLREDUCE(MPI_IN_PLACE,dprm,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_exit,info)
          dprm=dprm/dble(nz-1)/nproc_exit

          ! subtract mean value
          ! -------------------
          dpr(j,:)=dpr(j,:)-dprm
       enddo
       
    else
       
       ! Impose min/max value
       ! ====================
       ! iproc_send is proc having either k_min/k_max
       do j=1,ny
          if (iproc_exit==iproc_send) dprm=dpr(j,ind_m)

          ! send to other procs in communicator the value at k_min/k_max
          ! ------------------------------------------------------------
          call MPI_BCAST(dprm,1,MPI_DOUBLE_PRECISION,iproc_send,COMM_exit,info)
          
          ! impose min/max value
          ! --------------------
          dpr(j,:)=dpr(j,:)-dprm
       enddo
       
    endif
    
    ! annulus
    !!!dpr=dpr+1.123136382328114e+03
    !!!dpr=dpr-(1.-0.9999)*p_ref
    
!!$    ! Subtract mean value
!!$    ! ===================
!!$    ! local mean
!!$    ! ----------
!!$    dprmo=0.0_wp
!!$    do j=1,ny
!!$       do k=1,nz
!!$          dprmo(j)=dprmo(j)+dpr(j,k)
!!$       enddo
!!$    enddo
!!$    
!!$    ! global mean
!!$    ! -----------
!!$    !print *,iproc,nproc_exit
!!$    do j=1,ny
!!$       call MPI_ALLREDUCE(MPI_IN_PLACE,dprmo(j),1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_exit,info)
!!$    enddo
!!$    dprmo=dprmo/dble(nz-1)/nproc_exit
!!$
!!$    ! subtract mean value
!!$    ! -------------------
!!$    dpr(j,:)=dpr(j,:)-dprmo(j)
    
  end subroutine bc_rea_ox_kp

  !==============================================================================
  subroutine bc_rea_analytical(dpr,iflow)
  !==============================================================================
    !> Compute analytical radial pressure equilibrium for annulus case
    !> * analytical solution for free vortex flow    (iflow=1)
    !> * analytical solution for solid body rotation (iflow=2)
    !> /!\ ONLY 1 proc in radial direction (not parallelized) & imax face
  !==============================================================================
    use mod_fluid    ! for gas constant rg
    use mod_constant ! for reference quantities
    implicit none
    ! ---------------------------------------------------------------------------
    ! type of flow
    integer, intent(in) :: iflow
    ! pressure correction
    real(wp), dimension(ny,nz), intent(out) :: dpr 
    ! ---------------------------------------------------------------------------
    integer :: j,k
    real(wp) :: R1,R2 ! annulus radii
    real(wp) :: fac,den1,den2,rp2,rp ! working variables
    real(wp) :: ks,alpha_1,Ei1,Ei2 ! free vortex flow
    ! ---------------------------------------------------------------------------
    
    ! Annulus dimensions
    ! ==================
    R1=0.2_wp
    R2=0.28_wp

    if (iflow==1) then
       ! Free vortex flow
       ! ================

       ! strength
       ! --------
       ks=20.11_wp
       fac=ks**2/(2.0_wp*rg*T_ref)

       ! compute exponential integrals (in src/Mathematics)
       ! -----------------------------
       call e1xb(fac/R1**2,Ei1)
       call e1xb(fac/R2**2,Ei2)

       ! compute integration constant
       ! ----------------------------
       den1=-fac*Ei1+R1**2*exp(-fac/R1**2)
       den2=-fac*Ei2+R2**2*exp(-fac/R2**2)
       alpha_1=p_ref*(R2**2-R1**2)/(den2-den1)

       ! compute radial pressure profile
       ! -------------------------------
       do k=1,nz
          do j=1,ny
             rp2=yc3(nx,j,k)**2+zc3(nx,j,k)**2
             rp=sqrt(rp2)
             dpr(j,k)=alpha_1*exp(-fac/rp**2)-p_ref
          enddo
       enddo

    endif

  end subroutine bc_rea_analytical

end module mod_bc_rad_eq
