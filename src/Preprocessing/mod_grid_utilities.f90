!=================================================================================
module mod_grid_utilities
!=================================================================================
  !> Module with some adapted NASA grid utilities routines
!=================================================================================
  use precision
  implicit none
  ! ----------------------------------------------------------------------------
  ! Local constants
  real(wp), parameter, private :: zero=0.0_wp,one=1.0_wp
  real(wp), parameter, private :: two=2.0_wp,three=3.0_wp
  real(wp), parameter, private :: third=one/3.0_wp
  ! ----------------------------------------------------------------------------
  ! Grid
  integer, private :: ngx1,ngy1,ngz1
  real(wp), dimension(:,:), allocatable, private :: xgc1,ygc1
  real(wp), dimension(:,:,:), allocatable, private :: xgc31,ygc31,zgc31
  ! ----------------------------------------------------------------------------

contains
  
  !===============================================================================
  subroutine modif_grid
  !===============================================================================
    !> modif grid with NASA utilities: suppress half-cells at block junctions
  !===============================================================================
    use mod_mpi
    use mod_block
    use mod_constant ! for: dirGRID,nameGRID
    use mod_utils
    use warnstop
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,n,n_start,n_end
    logical :: check_grid
    ! intermediate grid
    integer :: ngx2,ngy2
    real(wp), dimension(:,:), allocatable :: xgc2,ygc2
    character(len=200) :: gridname_
    ! ----------------------------------------------------------------------------

    ! Print message
    ! =============
    if (iproc==0) then
       print *,repeat('=',70)
       print *,'Suppress half-cells at block junctions'
       print *,'(adapted from NASA grid utilities)'
       print *,repeat('=',70)
    endif

    ! Only possible on 1 proc
    ! =======================
    if (nproc.ne.1) then
       if (iproc==0) print *,'you can only run the grid modification with nproc=1'
       call mpistop('Change nproc to 1 or simply, ./musica. Shutting down...', 0)
    endif

    ! Indicator to check grid
    ! =======================
    check_grid=.false.
    if (check_grid) print *,'.. check_grid activated'

    ! Loop for all blocks
    ! ===================
    do n=1,nbloc

       ! block dimensions
       ! ----------------
       ngx1=bl(n)%ni
       ngy1=bl(n)%nj
       ngz1=bl(n)%nk

       ! read grid before modifications for block n
       ! -------------------------------------------
       allocate(xgc1(ngx1,ngy1),ygc1(ngx1,ngy1))

       if (is_coarse_grid) then
          gridname_=trim(dirGRID)//'/'//trim(nameGRID)//'c_bl'//trim(numchar(n))
       else
          gridname_=trim(dirGRID)//'/'//trim(nameGRID)//'_bl'//trim(numchar(n))
       endif

       call read_grid(trim(gridname_)//'.x')

       ! write grid in fortran bin for check with Matlab/Python
       ! ------------------------------------------------------
       if (check_grid) then
          open(194,file=trim(gridname_)//'.bin',form='unformatted',status='unknown')
          rewind(194)
          write(194) ngx1
          write(194) ngy1
          write(194) ((xgc1(i,j),i=1,ngx1),j=1,ngy1)
          write(194) ((ygc1(i,j),i=1,ngx1),j=1,ngy1)
          close(194)
       endif

       ! Grid modifications are carried out with changen2d [from NASA utilities]
       ! =================================================

       ! syntax of changen2d [all routines written in fortran90 are in mod_grid_utilities]
       ! --------------------
       ! changen2d: interpolates X1 & Y1 along a curve between points i1 & i2
       ! in order to redistribute them as a different number of points iA:iB
       ! with the end points matching.
       !     METHOD   C*1(Uppercase) Type of fit to be used:
       !            'M' means Monotonic piecewise cubics;
       !            'B' means non-monotonic "Bessel"-type piecewise cubics (looser fit);
       !            'L' means piecewise Linear fit;
       !            'C' means Cyclic (periodic) end conditions: loose fit assumed.

       ! Modifications in direction I (called here x)
       ! ============================================

       ! 0. check if it is necessary, ie if there is a neighboring block at imin and/or imax
       !    ------------------------
       if ((bl(n)%BC(1)>0).or.(bl(n)%BC(2)>0)) then
       !if (((bl(n)%BC(1)>0).and.(bl(n)%BC(1)/=n)).or.((bl(n)%BC(2)>0).and.(bl(n)%BC(2)/=n))) then

          ! 1. double the number of points along I (x-direction)
          !    -----------------------------------
          ngx2=ngx1*2
          ! allocate refined grid
          allocate(xgc2(ngx2,ngy1),ygc2(ngx2,ngy1))
          ! interpolate new grid with methods B (the most accurate)
          do j=1,ngy1
             call changen2d(1,ngx1,xgc1(:,j),ygc1(:,j),1,ngx2,xgc2(:,j),ygc2(:,j),'B')
          enddo

          ! 2. perform backward interpolation deleting half-cells near block junctions
          !    -----------------------------------------------------------------------
          ! bounds for backward interpolation
          n_start=1
          n_end=ngx2
          ! if there is a neighboring block at imin, suppress a half-cell
          if (bl(n)%BC(1)>0) n_start=n_start+1
          ! if there is a neighboring block at imax, suppress a half-cell
          if (bl(n)%BC(2)>0) n_end=n_end-1
          ! interpolate refined grid on old dimensions
          do j=1,ngy1
             call changen2d(1,n_end-n_start+1,xgc2(n_start:n_end,j),ygc2(n_start:n_end,j),1,ngx1,xgc1(:,j),ygc1(:,j),'B')
          enddo
          ! free arrays
          deallocate(xgc2,ygc2)

       endif

       ! Modifications in direction J (called here y)
       ! ============================================

       ! 0. check if it is necessary, ie if there is a neighboring block at jmin and/or jmax
       !    ------------------------
       if ((bl(n)%BC(3)>0).or.(bl(n)%BC(4)>0)) then
       !if (((bl(n)%BC(3)>0).and.(bl(n)%BC(4)/=n)).or.((bl(n)%BC(3)>0).and.(bl(n)%BC(4)/=n))) then

          ! 1. double the number of points along I (x-direction)
          !    -----------------------------------
          ngy2=ngy1*2
          ! allocate refined grid
          allocate(xgc2(ngx1,ngy2),ygc2(ngx1,ngy2))
          ! interpolate new grid with methods B (the most accurate)
          do i=1,ngx1
             call changen2d(1,ngy1,xgc1(i,:),ygc1(i,:),1,ngy2,xgc2(i,:),ygc2(i,:),'B')
          enddo

          ! 2. perform backward interpolation deleting half-cells near block junctions
          !    -----------------------------------------------------------------------
          ! bounds for backward interpolation
          n_start=1
          n_end=ngy2
          ! if there is a neighboring block at jmin, suppress a half-cell
          if (bl(n)%BC(3)>0) n_start=n_start+1
          ! if there is a neighboring block at jmax, suppress a half-cell
          if (bl(n)%BC(4)>0) n_end=n_end-1
          ! interpolate refined grid on old dimensions
          do i=1,ngx1
             call changen2d(1,n_end-n_start+1,xgc2(i,n_start:n_end),ygc2(i,n_start:n_end),1,ngy1,xgc1(i,:),ygc1(i,:),'B')
          enddo
          ! free arrays
          deallocate(xgc2,ygc2)

       endif

       ! write modified grid in fortran bin for check with Matlab/Python
       ! ---------------------------------------------------------------
       if (check_grid) then
          open(194,file=trim(gridname_)//'_mod.bin',form='unformatted',status='unknown')
          rewind(194)
          write(194) ngx1
          write(194) ngy1
          write(194) ((xgc1(i,j),i=1,ngx1),j=1,ngy1)
          write(194) ((ygc1(i,j),i=1,ngx1),j=1,ngy1)
          close(194)
       endif

       ! write grid after modifications for block n
       ! -------------------------------------------
       gridname_=trim(dirGRID)//'/'//trim(nameGRID)//'_mod_bl'//trim(numchar(n))

       call write_grid(trim(gridname_)//'.x',ngx1,ngy1,ngz1)

       ! free grid arrays
       ! ----------------
       deallocate(xgc1,ygc1)

    enddo ! end loop blocks

  end subroutine modif_grid

  !===============================================================================
  subroutine coarse_grid
  !===============================================================================
    !> create coarse grid (divided by a factor of 2 in each direction)
  !===============================================================================
    use mod_mpi
    use mod_block
    use mod_constant ! for: dirGRID,nameGRID
    use mod_utils
    use warnstop
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,n
    ! intermediate grid
    integer :: ngx2,ngy2
    real(wp), dimension(:,:), allocatable :: xgc2,ygc2
    character(len=200) :: gridname_
    ! ----------------------------------------------------------------------------

    ! Print message
    ! =============
    if (iproc==0) then
       print *,repeat('=',70)
       print *,'Create coarse grid 2D (divided by a factor of 2 in each direction)'
       print *,'(adapted from NASA grid utilities)'
       print *,repeat('=',70)
    endif

    ! Only possible on 1 proc
    ! =======================
    if (nproc.ne.1) then
       if (iproc==0) print *,'you can only run the grid modification with nproc=1'
       call mpistop('Change nproc to 1 or simply, ./musica. Shutting down...', 0)
    endif

    ! Loop for all blocks
    ! ===================
    do n=1,nbloc

       ! block dimensions
       ! ----------------
       ngx1=bl(n)%ni
       ngy1=bl(n)%nj
       ngz1=bl(n)%nk

       ! read grid before modifications for block n
       ! -------------------------------------------
       allocate(xgc1(ngx1,ngy1),ygc1(ngx1,ngy1))

       gridname_=trim(dirGRID)//'/'//trim(nameGRID)//'_bl'//trim(numchar(n))

       call read_grid(trim(gridname_)//'.x')

       ! Grid modifications are carried out with changen2d [from NASA utilities]
       ! =================================================

       ! syntax of changen2d [all routines written in fortran90 are in mod_grid_utilities]
       ! --------------------
       ! changen2d: interpolates X1 & Y1 along a curve between points i1 & i2
       ! in order to redistribute them as a different number of points iA:iB
       ! with the end points matching.
       !     METHOD   C*1(Uppercase) Type of fit to be used:
       !            'M' means Monotonic piecewise cubics;
       !            'B' means non-monotonic "Bessel"-type piecewise cubics (looser fit);
       !            'L' means piecewise Linear fit;
       !            'C' means Cyclic (periodic) end conditions: loose fit assumed.

       ! Modifications in direction I (called here x)
       ! ============================================

       ! divide by 2 the number of points along I
       ! ----------------------------------------
       ngx2=ngx1/2

       ! allocate refined grid
       allocate(xgc2(ngx2,ngy1),ygc2(ngx2,ngy1))

       ! interpolate new grid with methods B (the most accurate)
       ! -----------------------------------
       do j=1,ngy1
          call changen2d(1,ngx1,xgc1(:,j),ygc1(:,j),1,ngx2,xgc2(:,j),ygc2(:,j),'B')
       enddo

       ! store new grid in original array
       deallocate(xgc1,ygc1)
       allocate(xgc1(ngx2,ngy1),ygc1(ngx2,ngy1))
       xgc1=xgc2
       ygc1=ygc2
       deallocate(xgc2,ygc2)

       ! Modifications in direction J (called here y)
       ! ============================================

       ! divide by 2 the number of points along J
       ! ----------------------------------------
       ngy2=ngy1/2

       ! allocate refined grid
       allocate(xgc2(ngx2,ngy2),ygc2(ngx2,ngy2))

       ! interpolate new grid with methods B (the most accurate)
       ! -----------------------------------
       do i=1,ngx2
          call changen2d(1,ngy1,xgc1(i,:),ygc1(i,:),1,ngy2,xgc2(i,:),ygc2(i,:),'B')
       enddo

       ! store new grid in original array
       deallocate(xgc1,ygc1)
       allocate(xgc1(ngx2,ngy2),ygc1(ngx2,ngy2))
       xgc1=xgc2
       ygc1=ygc2
       deallocate(xgc2,ygc2)

       ! write grid after modifications for block n
       ! -------------------------------------------
       gridname_=trim(dirGRID)//'/'//trim(nameGRID)//'c_bl'//trim(numchar(n))

       call write_grid(trim(gridname_)//'.x',ngx2,ngy2,ngz1)

       ! free grid arrays
       ! ----------------
       deallocate(xgc1,ygc1)

    enddo ! end loop blocks

  end subroutine coarse_grid

  !=============================================================================
  subroutine read_grid(gridname)
  !=============================================================================
    !> Read grids in PLOT3D format
  !=============================================================================
    implicit none
    ! --------------------------------------------------------------------------
    character(len=*), intent(in) :: gridname
    ! --------------------------------------------------------------------------
    ! Local variables
    integer :: i,j
    integer :: mb,ni,nj,nk
    ! --------------------------------------------------------------------------

    open(50,file=gridname,form='formatted')
    rewind(50)
    read(50,*) mb
    read(50,*) ni,nj,nk
    read(50,*) ((xgc1(i,j),i=1,ni),j=1,nj)
    read(50,*) ((ygc1(i,j),i=1,ni),j=1,nj)
    close(50)

  end subroutine read_grid

  !=============================================================================
  subroutine write_grid(gridname,ni,nj,nk)
  !=============================================================================
    !> Write grids in PLOT3D format
  !=============================================================================
    implicit none
    ! --------------------------------------------------------------------------
    character(len=*), intent(in) :: gridname
    ! store infos read in plot3d files
    integer, intent(in) :: ni,nj,nk
    ! --------------------------------------------------------------------------
    ! Local variables
    integer :: i,j
    ! --------------------------------------------------------------------------

    open(50,file=gridname,form='formatted')
    rewind(50)
    write(50,*) 1
    write(50,*) ni,nj,nk
    write(50,*) ((xgc1(i,j),i=1,ni),j=1,nj)
    write(50,*) ((ygc1(i,j),i=1,ni),j=1,nj)
    close(50)

  end subroutine write_grid

  !===============================================================================
  subroutine modif_grid3d
  !===============================================================================
    !> modif grid with NASA utilities: suppress half-cells at block junctions
  !===============================================================================
    use mod_mpi
    use mod_block
    use mod_constant ! for: dirGRID,nameGRID
    use mod_utils
    use warnstop
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k,n,n_start,n_end
    logical :: check_grid
    ! intermediate grid
    integer :: ngx2,ngy2,ngz2
    real(wp), dimension(:,:,:), allocatable :: xgc32,ygc32,zgc32
    character(len=200) :: gridname_
    ! ----------------------------------------------------------------------------

    ! Print message
    ! =============
    if (iproc==0) then
       print *,repeat('=',70)
       print *,'Suppress half-cells at block junctions'
       print *,'(adapted from NASA grid utilities)'
       print *,repeat('=',70)
    endif

    ! Only possible on 1 proc
    ! =======================
    if (nproc.ne.1) then
       if (iproc==0) print *,'you can only run the grid modification with nproc=1'
       call mpistop('Change nproc to 1 or simply, ./musica. Shutting down...', 0)
    endif

    ! Indicator to check grid
    ! =======================
    check_grid=.false.
    if (check_grid) print *,'.. check_grid activated'

    ! Loop for all blocks
    ! ===================
    do n=1,nbloc

       ! block dimensions
       ! ----------------
       ngx1=bl(n)%ni
       ngy1=bl(n)%nj
       ngz1=bl(n)%nk

       ! read grid before modifications for block n
       ! -------------------------------------------
       allocate(xgc31(ngx1,ngy1,ngz1),ygc31(ngx1,ngy1,ngz1),zgc31(ngx1,ngy1,ngz1))

       if (is_coarse_grid) then
          gridname_=trim(dirGRID)//'/'//trim(nameGRID)//'c_bl'//trim(numchar(n))
       else
          gridname_=trim(dirGRID)//'/'//trim(nameGRID)//'_bl'//trim(numchar(n))
       endif

       call read_grid3d(trim(gridname_)//'.x')

       ! write grid in fortran bin for check with Matlab/Python
       ! ------------------------------------------------------
       if (check_grid) then
          open(194,file=trim(gridname_)//'.bin',form='unformatted',status='unknown')
          rewind(194)
          write(194) ngx1
          write(194) ngy1
          write(194) ((xgc31(i,j,1),i=1,ngx1),j=1,ngy1)
          write(194) ((ygc31(i,j,1),i=1,ngx1),j=1,ngy1)
          close(194)
       endif

       ! Grid modifications are carried out with changen3d [from NASA utilities]
       ! =================================================

       ! syntax of changen3d [all routines written in fortran90 are in mod_grid_utilities]
       ! --------------------
       ! changen3d: interpolates X1, Y1 & Z1 along a curve between points i1 & i2
       ! in order to redistribute them as a different number of points iA:iB
       ! with the end points matching.
       !     METHOD   C*1(Uppercase) Type of fit to be used:
       !            'M' means Monotonic piecewise cubics;
       !            'B' means non-monotonic "Bessel"-type piecewise cubics (looser fit);
       !            'L' means piecewise Linear fit;
       !            'C' means Cyclic (periodic) end conditions: loose fit assumed.
       ! /!\ use UPPERCASE letters (see changen3d)

       ! Modifications in direction I (called here x)
       ! ============================================

       if ((bl(n)%BC(1)>0).or.(bl(n)%BC(2)>0)) then
          ! 0. check if it is necessary, ie if there is a neighboring block at imin and/or imax
          !    ------------------------

          ! 1. double the number of points along I
          !    -----------------------------------
          ngx2=ngx1*2

          ! allocate refined grid
          allocate(xgc32(ngx2,ngy1,ngz1),ygc32(ngx2,ngy1,ngz1),zgc32(ngx2,ngy1,ngz1))
          ! interpolate new grid with methods B (the most accurate)
          do k=1,ngz1
             do j=1,ngy1
                call changen3d(1,ngx1,xgc31(:,j,k),ygc31(:,j,k),zgc31(:,j,k), &
                               1,ngx2,xgc32(:,j,k),ygc32(:,j,k),zgc32(:,j,k),'B')
             enddo
          enddo

          ! 2. perform backward interpolation deleting half-cells near block junctions
          !    -----------------------------------------------------------------------
          ! bounds for backward interpolation
          n_start=1
          n_end=ngx2
          ! if there is a neighboring block at imin, suppress a half-cell
          if (bl(n)%BC(1)>0) n_start=n_start+1
          ! if there is a neighboring block at imax, suppress a half-cell
          if (bl(n)%BC(2)>0) n_end=n_end-1
          ! interpolate refined grid on old dimensions
          do k=1,ngz1
             do j=1,ngy1
                call changen3d(1,n_end-n_start+1, &
                     xgc32(n_start:n_end,j,k),ygc32(n_start:n_end,j,k),zgc32(n_start:n_end,j,k), &
                     1,ngx1,xgc31(:,j,k),ygc31(:,j,k),zgc31(:,j,k),'B')
             enddo
          enddo
          ! free arrays
          deallocate(xgc32,ygc32,zgc32)

       endif

       ! Modifications in direction J (called here y)
       ! ============================================

       if ((bl(n)%BC(3)>0).or.(bl(n)%BC(4)>0)) then
          ! 0. check if it is necessary, ie if there is a neighboring block at jmin and/or jmax
          !    ------------------------

          ! 1. double the number of points along J
          !    -----------------------------------
          ngy2=ngy1*2
          ! allocate refined grid
          allocate(xgc32(ngx1,ngy2,ngz1),ygc32(ngx1,ngy2,ngz1),zgc32(ngx1,ngy2,ngz1))
          ! interpolate new grid with methods B (the most accurate)
          do k=1,ngz1
             do i=1,ngx1
                call changen3d(1,ngy1,xgc31(i,:,k),ygc31(i,:,k),zgc31(i,:,k), &
                               1,ngy2,xgc32(i,:,k),ygc32(i,:,k),zgc32(i,:,k),'B')
             enddo
          enddo

          ! 2. perform backward interpolation deleting half-cells near block junctions
          !    -----------------------------------------------------------------------
          ! bounds for backward interpolation
          n_start=1
          n_end=ngy2
          ! if there is a neighboring block at jmin, suppress a half-cell
          if (bl(n)%BC(3)>0) n_start=n_start+1
          ! if there is a neighboring block at jmax, suppress a half-cell
          if (bl(n)%BC(4)>0) n_end=n_end-1
          ! interpolate refined grid on old dimensions
          do k=1,ngz1
             do i=1,ngx1
                call changen3d(1,n_end-n_start+1, &
                     xgc32(i,n_start:n_end,k),ygc32(i,n_start:n_end,k),zgc32(i,n_start:n_end,k), &
                     1,ngy1,xgc31(i,:,k),ygc31(i,:,k),zgc31(i,:,k),'B')
             enddo
          enddo
          ! free arrays
          deallocate(xgc32,ygc32,zgc32)

       endif

       ! Modifications in direction K (called here z)
       ! ============================================

       if ((bl(n)%BC(5)>0).or.(bl(n)%BC(6)>0)) then
          ! 0. check if it is necessary, ie if there is a neighboring block at kmin and/or kmax
          !    ------------------------

          ! 1. double the number of points along K
          !    -----------------------------------
          ngz2=ngz1*2
          ! allocate refined grid
          allocate(xgc32(ngx1,ngy1,ngz2),ygc32(ngx1,ngy1,ngz2),zgc32(ngx1,ngy1,ngz2))
          ! interpolate new grid with methods B (the most accurate)
          do j=1,ngy1
             do i=1,ngx1
                call changen3d(1,ngz1,xgc31(i,j,:),ygc31(i,j,:),zgc31(i,j,:), &
                               1,ngz2,xgc32(i,j,:),ygc32(i,j,:),zgc32(i,j,:),'B')
             enddo
          enddo

          ! 2. perform backward interpolation deleting half-cells near block junctions
          !    -----------------------------------------------------------------------
          ! bounds for backward interpolation
          n_start=1
          n_end=ngz2
          ! if there is a neighboring block at kmin, suppress a half-cell
          if (bl(n)%BC(5)>0) n_start=n_start+1
          ! if there is a neighboring block at kmax, suppress a half-cell
          if (bl(n)%BC(6)>0) n_end=n_end-1
          ! interpolate refined grid on old dimensions
          do j=1,ngy1
             do i=1,ngx1
                call changen3d(1,n_end-n_start+1, &
                     xgc32(i,j,n_start:n_end),ygc32(i,j,n_start:n_end),zgc32(i,j,n_start:n_end), &
                     1,ngz1,xgc31(i,j,:),ygc31(i,j,:),zgc31(i,j,:),'B')
             enddo
          enddo
          ! free arrays
          deallocate(xgc32,ygc32,zgc32)

       endif

       ! write modified grid in fortran bin for check with Matlab/Python
       ! ---------------------------------------------------------------
       if (check_grid) then
          open(194,file=trim(gridname_)//'_mod.bin',form='unformatted',status='unknown')
          rewind(194)
          write(194) ngx1
          write(194) ngy1
          write(194) ((xgc31(i,j,1),i=1,ngx1),j=1,ngy1)
          write(194) ((ygc31(i,j,1),i=1,ngx1),j=1,ngy1)
          close(194)
       endif

       ! write grid after modifications for block n
       ! -------------------------------------------
       gridname_=trim(dirGRID)//'/'//trim(nameGRID)//'_mod_bl'//trim(numchar(n))

       call write_grid3d(trim(gridname_)//'.x',ngx1,ngy1,ngz1)

       ! free grid arrays
       ! ----------------
       deallocate(xgc31,ygc31,zgc31)

    enddo ! end loop blocks

  end subroutine modif_grid3d

  !===============================================================================
  subroutine coarse_grid3d
  !===============================================================================
    !> create coarse grid (divided by a factor of 2 in each direction)
  !===============================================================================
    use mod_mpi
    use mod_block
    use mod_constant ! for: dirGRID,nameGRID
    use mod_utils
    use warnstop
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k,n
    ! intermediate grid
    integer :: ngx2,ngy2,ngz2
    real(wp), dimension(:,:,:), allocatable :: xgc32,ygc32,zgc32
    character(len=200) :: gridname_
    ! ----------------------------------------------------------------------------

    ! Print message
    ! =============
    if (iproc==0) then
       print *,repeat('=',70)
       print *,'Create coarse grid 3D (divided by a factor of 2 in each direction)'
       print *,'(adapted from NASA grid utilities)'
       print *,repeat('=',70)
    endif

    ! Only possible on 1 proc
    ! =======================
    if (nproc.ne.1) then
       if (iproc==0) print *,'you can only run the grid modification with nproc=1'
       call mpistop('Change nproc to 1 or simply, ./musica. Shutting down...', 0)
    endif

    ! Loop for all blocks
    ! ===================
    do n=1,nbloc

       ! block dimensions
       ! ----------------
       ngx1=bl(n)%ni
       ngy1=bl(n)%nj
       ngz1=bl(n)%nk

       ! read grid before modifications for block n
       ! -------------------------------------------
       allocate(xgc31(ngx1,ngy1,ngz1),ygc31(ngx1,ngy1,ngz1),zgc31(ngx1,ngy1,ngz1))

       gridname_=trim(dirGRID)//'/'//trim(nameGRID)//'_bl'//trim(numchar(n))

       call read_grid3d(trim(gridname_)//'.x')

       ! Grid modifications are carried out with changen3d [from NASA utilities]
       ! =================================================

       ! syntax of changen3d [all routines written in fortran90 are in mod_grid_utilities]
       ! --------------------
       ! changen3d: interpolates X1, Y1 & Z1 along a curve between points i1 & i2
       ! in order to redistribute them as a different number of points iA:iB
       ! with the end points matching.
       !     METHOD   C*1(Uppercase) Type of fit to be used:
       !            'M' means Monotonic piecewise cubics;
       !            'B' means non-monotonic "Bessel"-type piecewise cubics (looser fit);
       !            'L' means piecewise Linear fit;
       !            'C' means Cyclic (periodic) end conditions: loose fit assumed.
       ! /!\ use UPPERCASE letters (see changen3d)

       ! Modifications in direction I (called here x)
       ! ============================================

       ! divide by 2 the number of points along I
       ! ----------------------------------------
       ngx2=ngx1/2

       ! allocate refined grid
       allocate(xgc32(ngx2,ngy1,ngz1),ygc32(ngx2,ngy1,ngz1),zgc32(ngx2,ngy1,ngz1))

       ! interpolate new grid with methods B (the most accurate)
       ! -----------------------------------
       do k=1,ngz1
          do j=1,ngy1
             call changen3d(1,ngx1,xgc31(:,j,k),ygc31(:,j,k),zgc31(:,j,k), &
                            1,ngx2,xgc32(:,j,k),ygc32(:,j,k),zgc32(:,j,k),'B')
          enddo
       enddo

       ! store new grid in original array
       deallocate(xgc31,ygc31,zgc31)
       allocate(xgc31(ngx2,ngy1,ngz1),ygc31(ngx2,ngy1,ngz1),zgc31(ngx2,ngy1,ngz1))
       xgc31=xgc32
       ygc31=ygc32
       zgc31=zgc32
       deallocate(xgc32,ygc32,zgc32)

       ! Modifications in direction J (called here y)
       ! ============================================

       ! divide by 2 the number of points along J
       ! ----------------------------------------
       ngy2=ngy1/2

       ! allocate refined grid
       allocate(xgc32(ngx2,ngy2,ngz1),ygc32(ngx2,ngy2,ngz1),zgc32(ngx2,ngy2,ngz1))

       ! interpolate new grid with methods B (the most accurate)
       ! -----------------------------------
       do k=1,ngz1
          do i=1,ngx2
             call changen3d(1,ngy1,xgc31(i,:,k),ygc31(i,:,k),zgc31(i,:,k), &
                            1,ngy2,xgc32(i,:,k),ygc32(i,:,k),zgc32(i,:,k),'B')
          enddo
       enddo

       ! store new grid in original array
       deallocate(xgc31,ygc31,zgc31)
       allocate(xgc31(ngx2,ngy2,ngz1),ygc31(ngx2,ngy2,ngz1),zgc31(ngx2,ngy2,ngz1))
       xgc31=xgc32
       ygc31=ygc32
       zgc31=zgc32
       deallocate(xgc32,ygc32,zgc32)

       ! Modifications in direction K (called here z)
       ! ============================================

       ! divide by 2 the number of points along J
       ! ----------------------------------------
       ngz2=ngz1/2

       ! allocate refined grid
       allocate(xgc32(ngx2,ngy2,ngz2),ygc32(ngx2,ngy2,ngz2),zgc32(ngx2,ngy2,ngz2))

       ! interpolate new grid with methods B (the most accurate)
       ! -----------------------------------
       do j=1,ngy2
          do i=1,ngx2
             call changen3d(1,ngz1,xgc31(i,j,:),ygc31(i,j,:),zgc31(i,j,:), &
                            1,ngz2,xgc32(i,j,:),ygc32(i,j,:),zgc32(i,j,:),'B')
          enddo
       enddo

       ! store new grid in original array
       deallocate(xgc31,ygc31,zgc31)
       allocate(xgc31(ngx2,ngy2,ngz2),ygc31(ngx2,ngy2,ngz2),zgc31(ngx2,ngy2,ngz2))
       xgc31=xgc32
       ygc31=ygc32
       zgc31=zgc32
       deallocate(xgc32,ygc32,zgc32)

       ! write grid after modifications for block n
       ! -------------------------------------------
       gridname_=trim(dirGRID)//'/'//trim(nameGRID)//'c_bl'//trim(numchar(n))

       call write_grid3d(trim(gridname_)//'.x',ngx2,ngy2,ngz2)

       ! free grid arrays
       ! ----------------
       deallocate(xgc31,ygc31,zgc31)

    enddo ! end loop blocks

    !call modif_grid3d(trim(griddir)//'/',trim(gridname)//'c')

  end subroutine coarse_grid3d

  !=============================================================================
  subroutine read_grid3d(gridname)
  !=============================================================================
    !> Read grids in PLOT3D format
  !=============================================================================
    implicit none
    ! --------------------------------------------------------------------------
    character(len=*), intent(in) :: gridname
    ! --------------------------------------------------------------------------
    ! Local variables
    integer :: i,j,k
    integer :: mb,ni,nj,nk
    ! --------------------------------------------------------------------------

    open(50,file=gridname,form='formatted')
    rewind(50)
    read(50,*) mb
    read(50,*) ni,nj,nk
    read(50,*) (((xgc31(i,j,k),i=1,ni),j=1,nj),k=1,nk)
    read(50,*) (((ygc31(i,j,k),i=1,ni),j=1,nj),k=1,nk)
    read(50,*) (((zgc31(i,j,k),i=1,ni),j=1,nj),k=1,nk)
    close(50)

  end subroutine read_grid3d

  !=============================================================================
  subroutine write_grid3d(gridname,ni,nj,nk)
  !=============================================================================
    !> Write grids in PLOT3D format
  !=============================================================================
    implicit none
    ! --------------------------------------------------------------------------
    character(len=*), intent(in) :: gridname
    ! store infos read in plot3d files
    integer, intent(in) :: ni,nj,nk
    ! --------------------------------------------------------------------------
    ! Local variables
    integer :: i,j,k
    ! --------------------------------------------------------------------------

    open(50,file=gridname,form='formatted')
    rewind(50)
    write(50,*) 1
    write(50,*) ni,nj,nk
    write(50,*) (((xgc31(i,j,k),i=1,ni),j=1,nj),k=1,nk)
    write(50,*) (((ygc31(i,j,k),i=1,ni),j=1,nj),k=1,nk)
    write(50,*) (((zgc31(i,j,k),i=1,ni),j=1,nj),k=1,nk)
    close(50)

  end subroutine write_grid3d

  !===============================================================================
  subroutine add_sponge(griddir,gridname)
  !===============================================================================
    !> modif grid with NASA utilities: suppress half-cells at block junctions
  !===============================================================================
    use mod_mpi
    use mod_block
    use mod_constant ! for: dirGRID,nameGRID
    use mod_utils
    use warnstop
    implicit none
    ! ----------------------------------------------------------------------------
    character(len=*), intent(in) :: griddir,gridname
    ! ----------------------------------------------------------------------------
    integer :: i,j,n,i1,i2,j1,j2
    logical :: check_grid
    ! new grid
    real(wp), dimension(:,:), allocatable :: xgc2,ygc2
    ! relative arc-length increments for all lines in the i and j directions
    real(wp), dimension(:,:,:), allocatable :: s0
    ! storage arrays for edge perturbations
    real(wp), dimension(:,:,:), allocatable :: dedgei,dedgej
    ! stretching rate and bissection parameter
    real(wp) :: dx0,rx,rxa,rxb
    integer :: it,nb
    real(wp) :: error,val_end,val
    real(wp), parameter :: epsi=1.0e-12_wp ! tolerance for bissection
    ! ----------------------------------------------------------------------------

    ! Print message
    ! =============
    if (iproc==0) then
       print *,repeat('=',70)
       print *,'Add stretching zone for exit blocks'
       print *,repeat('=',70)
    endif

    ! Only possible on 1 proc
    ! =======================
    if (nproc.ne.1) then
       if (iproc==0) print *,'you can only run the grid modification with nproc=1'
       call mpistop('Change nproc to 1 or simply, ./musica. Shutting down...', 0)
    endif

    ! Indicator to check grid
    ! =======================
    check_grid=.true.
    if (check_grid) print *,'.. check_grid activated'
    
    ! Loop for all blocks
    ! ===================
    do n=1,nbloc
    if ((n==5).or.(n==6).or.(n==12)) then

       ! block dimensions
       ! ----------------
       ngx1=bl(n)%ni
       ngy1=bl(n)%nj
       ngz1=bl(n)%nk

       ! read grid before modifications for block n
       ! -------------------------------------------
       allocate(xgc1(ngx1,ngy1),ygc1(ngx1,ngy1))
       allocate(xgc2(ngx1,ngy1),ygc2(ngx1,ngy1))
       
       call read_grid(trim(griddir)//'/'//trim(gridname)//'_bl'//trim(numchar(n))//'.x')

       ! write grid in fortran bin for check with Matlab/Python
       ! ------------------------------------------------------
       if (check_grid) then
          open(194,file=trim(gridname)//'_bl'//trim(numchar(n))//'.bin',form='unformatted',status='unknown')
          rewind(194)
          write(194) ngx1
          write(194) ngy1
          write(194) ((xgc1(i,j),i=1,ngx1),j=1,ngy1)
          write(194) ((ygc1(i,j),i=1,ngx1),j=1,ngy1)
          close(194)
       endif

       ! Grid modifications are carried out with warp2d [from NASA utilities]
       ! ==============================================

       ! allocate working arrays
       ! -----------------------
       allocate(s0(1:ngx1,1:ngy1,2))
       allocate(dedgei(2,1:ngy1,2),dedgej(2,1:ngx1,2))
       
       ! perturb grid edges
       ! ------------------
       ! copy unperturbed grid
       xgc2=xgc1
       ygc2=ygc1
       ! modify location of exit
       val_end=23.0_wp
       print *,'new exit at',val_end

       if (n==5) then
          print *,'block 5'
          print *,'-------'
          
          ! active zone (part of the grid to be perturbed)
          i1=1
          i2=ngx1
          j1=ngy1-20
          j2=ngy1
          
          ! modify imin edge (at i1, bottom)
          ! --------------------------------
          dx0=xgc1(i1,j1)-xgc1(i1,j1-1)
          val=val_end-xgc2(i1,j1-1)
          nb=j2-j1+1

          ! bissection to determine rx
          rxa=1.0000000001
          rxb=1.5
          it=1
          error=1.0_wp
          do while ((error>epsi).and.(it<1000))
             rx=(rxa+rxb)/2.

             if (f_rx(rxa)*f_rx(rx)>0) then
                rxa=rx
                rxb=rxb
             else
                rxa=rxa
                rxb=rx
             endif
             error=abs(f_rx(rx))
             it=it+1
          enddo
          if (it>=1000) call mpistop('failed to converge bissection (it>1000).', 0)
          print *,'at i1 (bottom): rx',rx
          
          do j=j1,j2
             xgc2(i1,j)=xgc2(i1,j-1)+dx0
             dx0=rx*dx0
          enddo
          
          ! modify imax edge (at i2, top)
          ! -----------------------------
          dx0=xgc1(i1,j1)-xgc1(i1,j1-1)
          val=val_end-xgc2(i2,j1-1)
          nb=j2-j1+1

          ! bissection to determine rx
          rxa=1.0000000001
          rxb=1.5
          it=1
          error=1.0_wp
          do while ((error>epsi).and.(it<1000))
             rx=(rxa+rxb)/2.

             if (f_rx(rxa)*f_rx(rx)>0) then
                rxa=rx
                rxb=rxb
             else
                rxa=rxa
                rxb=rx
             endif
             error=abs(f_rx(rx))
             it=it+1
          enddo
          if (it>=1000) call mpistop('failed to converge bissection (it>1000).', 0)
          print *,'at i2 (top): rx',rx
          
          do j=j1,j2
             xgc2(i2,j)=xgc2(i2,j-1)+dx0
             dx0=rx*dx0
          enddo
                    
          ! modify jmax edge (at j2, right)
          ! -------------------------------
          xgc2(i1:i2,j2)=val_end

       elseif (n==6) then
          print *,'block 6'
          print *,'-------'

          ! active zone (part of the grid to be perturbed)
          i1=1
          i2=21
          j1=1
          j2=ngy1
          
          !dx0=xgc1(i2,j1)-xgc1(i2+1,j1)
          !rx=1.06
          !print *,xgc2(i2+1,j1),xgc2(i2+1,j1)+dx0*(1-rx**21)/(1-rx)
          !print *,xgc2(i2+1,j2),xgc2(i2+1,j2)+dx0*(1-rx**21)/(1-rx)

          ! modify jmin edge (at j1, bottom)
          ! --------------------------------
          dx0=xgc1(i2,j1)-xgc1(i2+1,j1)
          val=val_end-xgc2(i2+1,j1)
          nb=i2-i1+1

          ! bissection to determine rx
          rxa=1.0000000001
          rxb=1.5
          it=1
          error=1.0_wp
          do while ((error>epsi).and.(it<1000))
             rx=(rxa+rxb)/2.

             if (f_rx(rxa)*f_rx(rx)>0) then
                rxa=rx
                rxb=rxb
             else
                rxa=rxa
                rxb=rx
             endif
             error=abs(f_rx(rx))
             it=it+1
          enddo
          if (it>=1000) call mpistop('failed to converge bissection (it>1000).', 0)
          print *,'at j1 (bottom): rx',rx
         
          do i=i2,i1,-1
             xgc2(i,j1)=xgc2(i+1,j1)+dx0
             dx0=rx*dx0
          enddo

          ! modify jmax edge (at j2, top)
          ! -----------------------------
          dx0=xgc1(i2,j1)-xgc1(i2+1,j1)
          val=val_end-xgc2(i2+1,j2)
          nb=i2-i1+1

          ! bissection to determine rx
          rxa=1.0000000001
          rxb=1.5
          it=1
          error=1.0_wp
          do while ((error>epsi).and.(it<1000))
             rx=(rxa+rxb)/2.

             if (f_rx(rxa)*f_rx(rx)>0) then
                rxa=rx
                rxb=rxb
             else
                rxa=rxa
                rxb=rx
             endif
             error=abs(f_rx(rx))
             it=it+1
          enddo
          if (it>=1000) call mpistop('failed to converge bissection (it>1000).', 0)
          print *,'at j2 (top): rx',rx

          do i=i2,i1,-1
             xgc2(i,j2)=xgc2(i+1,j2)+dx0
             dx0=rx*dx0
          enddo

          ! modify imin edge (at i1, right)
          ! -------------------------------
          xgc2(i1,j1:j2)=val_end
          
       elseif (n==12) then
          print *,'block 12'
          print *,'--------'
          
          ! active zone (part of the grid to be perturbed)
          i1=ngx1-20
          i2=ngx1
          j1=1
          j2=ngy1

          !dx0=xgc1(i1,j1)-xgc1(i1-1,j1)
          !rx=1.06
          !print *,xgc2(i1-1,j1),xgc2(i1-1,j1)+dx0*(1-rx**21)/(1-rx)
          !print *,xgc2(i1-1,j2),xgc2(i1-1,j2)+dx0*(1-rx**21)/(1-rx)
          
          ! modify jmin edge (at j1, top)
          ! -----------------------------
          dx0=xgc1(i1,j1)-xgc1(i1-1,j1)
          val=val_end-xgc2(i1-1,j1)
          nb=i2-i1+1
          
          ! bissection to determine rx
          rxa=1.0000000001
          rxb=1.5
          it=1
          error=1.0_wp
          do while ((error>epsi).and.(it<1000))
             rx=(rxa+rxb)/2.

             if (f_rx(rxa)*f_rx(rx)>0) then
                rxa=rx
                rxb=rxb
             else
                rxa=rxa
                rxb=rx
             endif
             error=abs(f_rx(rx))
             it=it+1
          enddo
          if (it>=1000) call mpistop('failed to converge bissection (it>1000).', 0)
          print *,'at j1 (top): rx',rx
          
          do i=i1,i2
             xgc2(i,j1)=xgc2(i-1,j1)+dx0
             dx0=rx*dx0
          enddo

          ! modify jmax edge (at j2, bottom)
          ! --------------------------------
          dx0=xgc1(i1,j1)-xgc1(i1-1,j1)
          val=val_end-xgc2(i1-1,j2)
          nb=i2-i1+1
          
          ! bissection to determine rx
          rxa=1.0000000001
          rxb=1.5
          it=1
          error=1.0_wp
          do while ((error>epsi).and.(it<1000))
             rx=(rxa+rxb)/2.

             if (f_rx(rxa)*f_rx(rx)>0) then
                rxa=rx
                rxb=rxb
             else
                rxa=rxa
                rxb=rx
             endif
             error=abs(f_rx(rx))
             it=it+1
          enddo
          if (it>=1000) call mpistop('failed to converge bissection (it>1000).', 0)
          print *,'at j2 (bottom): rx',rx

          do i=i1,i2
             xgc2(i,j2)=xgc2(i-1,j2)+dx0
             dx0=rx*dx0
          enddo
           
          ! modify imax edge (at i2, right)
          ! -------------------------------
          xgc2(i2,j1:j2)=val_end
       endif

       ! define relative arc-lengths
       ! ---------------------------
       call paramxy(1,ngx1,1,ngy1,i1,i2,j1,j2,xgc1,ygc1,s0)

       ! perturb interior of a plane grid given new edges
       ! ------------------------------------------------
       call warp2d(1,ngx1,1,ngy1,i1,i2,j1,j2,xgc1,ygc1,s0,dedgei,dedgej,xgc2,ygc2)

       ! replace grid
       ! ------------
       xgc1=xgc2
       ygc1=ygc2
    
       ! write modified grid in fortran bin for check with Matlab/Python
       ! ---------------------------------------------------------------
       if (check_grid) then
          open(194,file=trim(gridname)//'_bl'//trim(numchar(n))//'_mod.bin',form='unformatted',status='unknown')
          rewind(194)
          write(194) ngx1
          write(194) ngy1
          write(194) ((xgc1(i,j),i=1,ngx1),j=1,ngy1)
          write(194) ((ygc1(i,j),i=1,ngx1),j=1,ngy1)
          close(194)
       endif

       ! write grid after modifications for block n
       ! -------------------------------------------
       call write_grid(trim(griddir)//'/'//trim(gridname)//'_mod_bl'//trim(numchar(n))//'.x',ngx1,ngy1,ngz1)

       ! free arrays
       ! -----------
       deallocate(xgc1,ygc1,xgc2,ygc2)
       deallocate(s0,dedgei,dedgej)

    endif
    enddo ! end loop blocks

  contains

    !=================================================================================
    ! Bissection method to determine the stretching rate
    !=================================================================================
    function f_rx(rx)
      implicit none
      real(wp), intent(in) :: rx
      real(wp) :: f_rx

      f_rx=val-dx0*(1-rx**nb)/(1-rx)

    end function f_rx

  end subroutine add_sponge

  !! **************************
  !! NASA ROUTINES STARTS HERE (-> converted in f90 by XG)
  !! **************************

  !===============================================================================
  subroutine changen2d(i1,i2,X1,Y1,iA,iB,X2,Y2,method)
  !===============================================================================
    !> Change # pts. on a 2-sp. curve;  same relative distrbn.
  !===============================================================================
    !  changen2d interpolates X1 & Y1 along a curve between points I1 & I2
    !  in order to redistribute them as a different number of points iA:iB
    !  with the end points matching.
    !
    !  12/22/97  DAS  3-space grid line utility changen wrapped around plscrv3d.
    !  10/23/02   "   2-space variant: changen2d.
    !
    !  Author:  David Saunders, ELORET/NASA Ames Research Ctr., Moffett Field, CA
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! I/O arguments
    integer :: i1,i2          ! I First & last data points
    real(wp) :: X1(i2),Y1(i2) ! I Data points
    integer :: iA,iB          ! I First & last redistributed points
    real(wp) :: X2(iB),Y2(iB) ! O Redistributed coordinates
    character :: method*1     ! I Type of fits to be used by LCSFIT, q.v.
    ! ----------------------------------------------------------------------------
    ! Local variables.
    integer :: i,iP,ndata,ninterp
    real(wp) :: t(i1:i2),teval(iA:iB),deriv(iA:iB),p,r,ri1,rip,total
    logical, parameter :: inew=.true.,norm=.false.
    ! ----------------------------------------------------------------------------

    ndata=i2-i1+1

    call chords2d(ndata,X1(i1),Y1(i1),norm,total,t)

    r=dble(i2-i1)/dble(iB-iA)
    ri1=dble(i1)

    do i=iA,iB-1
       rip= ri1 + r*dble(i-iA)
       ip= int(rip)
       p= rip - dble(ip)
       teval(i)= (one-p)*t(ip) + p*t(ip+1)
    end do

    ninterp=iB-iA

    call lcsfit(ndata,t(i1),X1(i1),inew,method,ninterp,teval,X2(iA),deriv)

    X2(iB)=X1(i2)

    call lcsfit(ndata,t(i1),Y1(i1),inew,method,ninterp,teval,Y2(iA),deriv)

    Y2(iB)=Y1(i2)

  end subroutine changen2d

  !===============================================================================
  subroutine chords2d(n,x,y,normaliz,total,chord)
  !===============================================================================
    !> Cumulative [relative] chord-lengths for 2-space geometric curve
  !===============================================================================
    !  DESCRIPTION:
    !        chordS2D computes the cumulative Euclidean distances for two or
    !     more points on a 2-space curve represented by two arrays, with an
    !     option to normalize all values by the total distance.  Thus chord(1)
    !     is returned as 0. and chord(n) may be either 1. or total depending
    !     on whether normaliz is true or false.
    !
    !        chordS2D was introduced in spite of the existing chord function
    !     because use of chord with the double precision DTNURBS library in a
    !     single precision application such as SMOOTH would clash with prior
    !     use of the single precision version of chord.  The functionality
    !     is a little different (multiple values per call, and an option to
    !     normalize), and since NURBS are intended for geometric data (i.e.,
    !     x and y expected to have similar units), the careful safeguarding
    !     of chord is eschewed.
    !
    !  ARGUMENTS:
    !     Name    Type/Dimension I/O/S  Description
    !     n         I            I      Number of points on the curve. n >= 2.
    !     x         R(n)         I      Array of abscissas.
    !     y         R(n)         I      Array of ordinates.
    !     normaliz  L            I      .TRUE. means normalize results
    !                                   to the interval [0, 1].
    !     total     R              O    Total chord length, returned
    !                                   because it is otherwise lost
    !                                   if results are normalized.
    !     chord     R(n)           O    Cumulative chord lengths as
    !                                   described above.
    !
    !  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
    !  HISTORY:
    !     13 Mar. 1992   DAS   Adapted from chord for use with DTNURBS.
    !     08 Dec. 1993    "    Added total argument.  
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! I/O arguments
    integer :: n
    real(wp) :: chord(n),total,x(n),y(n)
    logical :: normaliz
    ! ----------------------------------------------------------------------------
    ! Local variables.
    integer :: i
    real(wp) :: dinv
    ! ----------------------------------------------------------------------------

    chord(1)=zero

    do i=2,n
       chord(i)=chord(i-1)+sqrt((x(i)-x(i-1))**2+(y(i)-y(i-1))**2)
    enddo

    total=chord(n)

    if (normaliz) then
       dinv= one/total
       do i=2,n
          chord(i)=chord(i)*dinv
       enddo
       chord(n)=one
    endif

  end subroutine chords2d

  !===============================================================================
  subroutine changen3d(i1,i2,X1,Y1,Z1,iA,iB,X2,Y2,Z2,method)
  !===============================================================================
    !> Change # pts. on a 3-sp. curve;  same relative distrbn.
  !===============================================================================
    !  changen3d interpolates X1, Y1 & Z1 along a curve between points I1 & I2
    !  in order to redistribute them as a different number of points iA:iB
    !  with the end points matching.
    !
    ! use chords3d : Arc-length utility
    !     plscrv3d : 3-space local cubic spline utility
    !     upcase   : Character string utility ensuring uppercase
    !
    !  12/22/97  DAS  3-space grid line utility changen wrapped around plscrv3d.
    !  01/11/08   "   Added option to make the output distribution uniform.
    !                 See expanded use of 'method'.
    !
    !  Author:  David Saunders, ELORET/NASA Ames Research Ctr., Moffett Field, CA
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! I/O arguments
    integer :: i1,i2                 ! I First & last data points
    real(wp) :: X1(i2),Y1(i2),Z1(i2) ! I Data points
    integer :: iA,iB                 ! I First & last redistributed points
    real(wp) :: X2(iB),Y2(iB),Z2(iB) ! O Redistributed coordinates
    character :: method*1            ! I Type of fits to be used by LCSFIT, q.v.
                                     ! M (monotonic), B ("loose"), or L (linear);
                                     ! lower case m, b, or l indicates that the
                                     ! output line is (essentially) uniform
    ! ----------------------------------------------------------------------------
    ! Local variables.
    integer :: i,iP,ieval,ndata
    real(wp) :: t(i2),dt,p,r,ri1,rip,teval,total
    real(wp) :: derivs(3)
    logical :: inew,uniform
    logical, parameter :: norm=.false.,closed=.false.
    character :: mupper*1
    ! ----------------------------------------------------------------------------

    ! suppresses derivatives
    derivs=-999.0_wp

    mupper=method

    call upcase(mupper)

    uniform= (mupper.ne.method)

    ndata=i2-i1+1

    call chords3d(ndata,X1(i1),Y1(i1),Z1(i1),norm,total,t(i1))

    ieval=i1
    inew=.true.

    X2(iA)=X1(i1)
    Y2(iA)=Y1(i1)
    Z2(iA)=Z1(i1)

    if (uniform) then

       dt=total/dble(iB-iA)

       do i=iA+1,iB-1

          teval=dt*dble(i-iA)

          call plscrv3d(ndata,X1(i1),Y1(i1),Z1(i1),t(i1), &
               mupper,inew,closed,teval,ieval,X2(i),Y2(i),Z2(i),derivs)
          inew=.false.

       enddo

    else ! Preserve original form of spacing as much as possible

       r=dble(i2-i1)/dble(iB-iA)
       ri1=dble(i1)

       do i=iA+1,iB-1

          rip=ri1+r*dble(i-iA)
          ip=int(rip)
          p=rip-dble(ip)
          teval=(one-p)*t(ip)+P*t(ip+1)

          call plscrv3d(ndata,X1(i1),Y1(i1),Z1(i1),t(i1),&
                        mupper,inew,closed,teval,ieval,X2(i),Y2(i),Z2(i),derivs)
          inew=.false.

       enddo

    endif

    X2(iB)=X1(i2)
    Y2(iB)=Y1(i2)
    Z2(iB)=Z1(i2)

  end subroutine changen3d

  !===============================================================================
  subroutine chords3d(n,x,y,z,normaliz,total,chord)
  !===============================================================================
    !> Cumulative [relative] chord-lengths for 3-space geometric curve
  !===============================================================================
    !  DESCRIPTION:
    !        chordS3d computes the cumulative Euclidean distances for two or
    !     more points on a 3-space curve represented by three arrays, with an
    !     option to normalize all values by the TOTAL distance.  Thus chord(1)
    !     is returned as 0. and chord(N) may be either 1. or TOTAL depending
    !     on whether NORMALIZ is true or false.
    !
    !        chordS3d was introduced in spite of the existing chord3d function
    !     because use of chord3d with the double precision DTNURBS library in a
    !     single precision application such as SMOOTH would clash with prior
    !     use of the single precision version of chord3d.  The functionality
    !     is a little different (multiple values per call, and an option to
    !     normalize), and since NURBS are intended for geometric data (i.e.,
    !     X and Y expected to have similar units), the careful safeguarding
    !     of chord3d is eschewed.
    !
    !  ARGUMENTS:
    !     Name    Type/Dimension  I/O/S  Description
    !     n         I             I      Number of points on the curve. n >= 2.
    !     x         R(n)          I      Coordinate x of data points.
    !     y         R(n)          I      Coordinate x of data points.
    !     z         R(n)          I      Coordinate x of data points.
    !     normaliz  L             I      .true. means normalize results
    !                                    to the interval [0, 1].
    !     total     R               O    Total chord length, returned
    !                                    because it is otherwise lost
    !                                    if results are normalized.
    !     chord     R(n)            O    Cumulative chord lengths as
    !
    !  AUTHOR:  David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
    !  HISTORY:
    !     13 Mar. 1992   DAS   Analog of chords2d for use with DTNURBS.
    !     08 Dec. 1993    "    Added total argument.
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! I/O arguments
    integer :: n
    real(wp) :: chord(n),total,x(n),y(n),z(n)
    logical :: normaliz
    ! ----------------------------------------------------------------------------
    ! Local variables.
    integer :: i
    real(wp) :: dinv
    ! ----------------------------------------------------------------------------

    chord(1)=zero

    do i=2,n
       chord(i)=chord(i-1)+sqrt((x(i)-x(i-1))**2+(y(i)-y(i-1))**2+(z(i)-z(i-1))**2)
    enddo

    total=chord(n)

    if (normaliz) then
       dinv= one/total
       do i=2,n
          chord(i)=chord(i)*dinv
       enddo
       chord(n)=one
    endif

  end subroutine chords3d

  !===============================================================================
  subroutine plscrv3d (ndata,X,Y,Z,t,method,inew,closed, &
                       teval,ieval,Xeval,Yeval,Zeval,derivs)
  !===============================================================================
    !> Storage-efficient local cubic spline fit (3-space)
    !> (monotonic and piecewise linear options too)
  !===============================================================================
    ! Description and usage:
    !        plscrv3d interpolates X,Y,Z along a curve at a specified teval using
    !     piecewise linear or local cubic spline techniques.  It is an adaptation of
    !     the earlier PLSFIT3D and PLSCURVE intended to avoid internal arc length
    !     calculations by assuming that the cumulative chord lengths associated with
    !     the data points are supplied by the application. (PLSFIT3D avoided storing
    !     all the Ts as part of its space-efficient graphics-oriented approach, but
    !     in practical grid generation applications the Ts are needed anyway.)
    !
    !        plscrv3d restores the method options that were dropped by PLSCURVE,
    !     in anticipation of a requirement for piecewise linear results.  The option
    !     to return derivatives for intersection calculations influenced the choice
    !     between one and many target Ts: a single T wins.  Use of normalized Ts is
    !     recommended for intersection calculations, so T is known to be in [0, 1],
    !     but other applications may prefer not to normalize, and this is permitted.
    !     T(*) may also be in descending order.
    !
    !        Arithmetic for unwanted derivatives is avoided by calls with input
    !     derivs(1) = -999.
    !
    !        See Robert Kennelly's original 2-space PLSFIT for more on compact local
    !     spline methods.

    !     USE:
    !        bessel  : 1st deriv. (central 3-point formula)
    !        brodlie : 1st deriv. (central) adjusted for monotonicity
    !        butland : 1st deriv. (non-central)  "    "    "
    !        threept : 1st deriv. (non-central 3-point formula)
    !        interval: Interpolatory 1-D search utility
    !
    !     Author:  David Saunders, Sterling Software and ELORET/NASA Ames Research Ctr., Moffett Field, CA
    !     09/22/97  DAS  Initial adaptation of PLSFIT3D/PLSCURVE with more efficient
    !                    curve/surface intersections and wing surface grids in mind,
    !                    especially ndata = 2 and piecewise linear cases.
    !     11/17/97   "   Allowed for non-normalized input Ts.
    !     05/06/98   "   Minor Fortran 90 updates.
    !     07/25/06   "   Allowed for descending t(*) order.
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! I/O arguments
    ! -------------
    ! Length of X, Y, Z input data arrays.
    integer, intent(in) :: ndata
    ! Input curve coordinates. Successive data points must be distinct.
    ! If closed, then first and last points must agree (not checked here).
    real(wp), intent(in) :: X(ndata),Y(ndata),Z(ndata)
    ! Parametric variable values associated with the data points,
    ! usually normalized cumulative chord lengths, but they need
    ! not be normalized or increasing except for intersection calculations
    ! via INTSEC4 or -5.
    real(wp), intent(in) :: t(ndata)
    ! The type of parametric fit to be used:
    !    'M' means monotonic piecewise cubics;
    !    'B' means non-monotonic "Bessel"-type piecewise cubics (looser fit);
    !    'L' means piecewise linear (only way to guarantee no overshoots).
    character, intent(in) :: method*1
    ! If inew=T, the search for a bracket starts from scratch,
    ! otherwise locally-saved search and fit information will be assumed to be
    ! correct. If calling plscrv3d from within a loop,
    ! set inew=F after the first call.
    logical, intent(in) :: inew
    ! closed=T means periodic boundary conditions are to be used.
    ! The curve should wrap around smoothly on itself.
    ! The calling routine must ensure that the end points agree.
    logical, intent(in) :: closed
    ! Target value of t at which to evaluate X,Y,Z.
    real, intent(in) :: teval
    ! Input with the data index at which to start searching for teval;
    ! output with actual index
    integer, intent(inout) :: ieval
    ! Interpolated coordinates.
    real, intent(out) :: Xeval, Yeval, Zeval
    ! Partial derivatives of X, Y, Z with respect to t.
    ! Enter derivs(1)=-999. to suppress these,
    ! else  derivs(1)=0. (say) on the first call and something other than -999. thereafter.
    real(wp), intent(inout) :: derivs(3)
    ! ----------------------------------------------------------------------------
    ! Local variables
    integer    ind (-1:2),j,left,right
    logical :: deriv,linear,loose,memory,twopts
    real(wp) :: arrow,bx(0:1),by(0:1),bz(0:1),cx,cy,cz,dt
    real(wp) :: delx(-1:1),dely(-1:1),delz(-1:1),dx,dy,dz
    real(wp) :: h(-1:1),rh,tbyarrow,tleft,tright
    ! ----------------------------------------------------------------------------
    ! Storage
    save arrow,bx,by,bz,cx,cy,cz,dx,dy,dz,left,right, &
         tleft,tright,deriv,linear,loose,twopts
    ! ----------------------------------------------------------------------------

    ! The 2-data-point case is common enough for intersections that it
    ! is treated specially. Minor code duplication is OK.

    if (inew) then

       deriv = (derivs(1).ne.-999.)
       twopts= (ndata==2)

       if (twopts) then
          rh   = one/(t(2)-t(1))
          bx(0)= (X(2)-X(1))*rh
          by(0)= (Y(2)-Y(1))*rh
          bz(0)= (Z(2)-Z(1))*rh

          dt   = teval-t(1)
          Xeval= X(1)+dt*bx(0)
          Yeval= Y(1)+dt*by(0)
          Zeval= Z(1)+dt*bz(0)

          if (deriv) then
             derivs(1)= bx(0)
             derivs(2)= by(0)
             derivs(3)= bz(0)
          endif

          return
       else
          memory= .false.
          linear= (method== 'L')
          loose = (method== 'B') ! Else monotonic
          arrow = sign(one,t(2)-t(1))
       endif

    elseif (twopts) then

       dt   = teval-t(1)
       Xeval= X(1)+dt*bx(0)
       Yeval= Y(1)+dt*by(0)
       Zeval= Z(1)+dt*bz(0)

       if (deriv) then  ! INTSEC5 overwrites previous derivatives
          derivs(1)= bx(0)
          derivs(2)= by(0)
          derivs(3)= bz(0)
       endif

       return

    else ! inew=F and ndata>2

       ! We can save time when plscrv3d is being called from within a loop
       ! by setting memory if possible.  The out-of-range checking relies on
       ! the fact that right = left+1 after the last return.  Cater to the
       ! more likely case of teval in the previous (same) interior interval.

       tbyarrow= teval*arrow
       memory= (tbyarrow>=tleft*arrow).and.(tbyarrow<tright*arrow)

       if (.not.memory) then
          memory= ((left==1).and.(tbyarrow<tright*arrow)) &
              .or.((right==ndata).and.(tbyarrow>=tleft*arrow))
       endif

    endif

    if (.not.memory) then

       ! Interpolation search for bracketing interval:

       left= ieval

       call interval(ndata,t,teval,arrow,left)

       ieval = left
       tleft = t(left)
       right = left+1
       tright= t(right)

       ! -------------------------------------------
       ! |   1 <= left < right = left+1 <= ndata   |
       ! -------------------------------------------

       if (linear) then ! Piecewise linear, ndata > 2

          rh   = one/(tright-tleft)
          bx(0)= (X(right)-X(left))*rh
          by(0)= (Y(right)-Y(left))*rh
          bz(0)= (Z(right)-Z(left))*rh

       else ! Piecewise cubic

          ! Three cases are handled together:(1) both endpoints are
          ! interior to the data range, or the interval is at the
          ! beginning/end of the range with either(2) periodic,
          ! or(3) free end conditions.  For special cases(2) and(3),
          ! the initial calculations will be overridden below.

          ind(-1)= left-1
          ind( 0)= left
          ind(+1)= right
          ind(+2)= right+1

          ! Patch index array for periodic case(2). (Later, the free end
          ! case will be overridden again, following derivative calculation.)

          if (left==1) then
             ind(-1)= ndata-1 ! Left side wrap-around boundary condition
          else if (right==ndata) then
             ind(2)= 2        ! Right side
          endif

          ! Interval and derivative approximations:

          h(-1)= t(ind(-1)+1)-t(ind(-1))
          h( 0)= t(right)-t(left)
          h(+1)= t(ind(2))-t(ind(2)-1)

          do j=-1,1
             rh     = one/h(j)
             delx(j)= (X(ind(j+1))-X(ind(j)))*rh
             dely(j)= (Y(ind(j+1))-Y(ind(j)))*rh
             delz(j)= (Z(ind(j+1))-Z(ind(j)))*rh
          enddo

          ! Select the interpolation scheme, and compute [adjusted] first
          ! derivatives at both left- and right-hand endpoints of the interval:

          if (loose) then ! Bessel-use the central formula at 0,+1

             bx(0)= bessel(0,h,delx)
             bx(1)= bessel(1,h,delx)
             by(0)= bessel(0,h,dely)
             by(1)= bessel(1,h,dely)
             bz(0)= bessel(0,h,delz)
             bz(1)= bessel(1,h,delz)

          else ! Monotone-use the Brodlie modification of Butland's formula

             bx(0)= brodlie(0,h,delx)
             bx(1)= brodlie(1,h,delx)
             by(0)= brodlie(0,h,dely)
             by(1)= brodlie(1,h,dely)
             bz(0)= brodlie(0,h,delz)
             bz(1)= brodlie(1,h,delz)

          endif

          !  Patch initial/final derivatives if not periodic, case(3):

          if (.not.closed) then
             if (left==1) then
                if (loose) then
                   bx(0)= threept(0,h,delx)
                   by(0)= threept(0,h,dely)
                   bz(0)= threept(0,h,delz)
                else
                   bx(0)= butland(0,h,delx)
                   by(0)= butland(0,h,dely)
                   bz(0)= butland(0,h,delz)
                endif
             elseif (right==ndata) then
                if (loose) then
                   bx(1)= threept(1,h,delx)
                   by(1)= threept(1,h,dely)
                   bz(1)= threept(1,h,delz)
                else
                   bx(1)= butland(1,h,delx)
                   by(1)= butland(1,h,dely)
                   bz(1)= butland(1,h,delz)
                endif
             endif
          endif

          ! Compute the remaining cubic coefficients relative to the left-hand endpoint.

          rh= one/h(0)
          cx= (three*delx(0)-two*bx(0)-bx(1))*rh
          cy= (three*dely(0)-two*by(0)-by(1))*rh
          cz= (three*delz(0)-two*bz(0)-bz(1))*rh
          dx= (-two*delx(0)+bx(0)+bx(1))*rh**2
          dy= (-two*dely(0)+by(0)+by(1))*rh**2
          dz= (-two*delz(0)+bz(0)+bz(1))*rh**2

       endif

    endif

    ! Evaluate the identified polynomials:

    dt= teval-tleft

    if (linear) then

       Xeval= X(left)+dt*bx(0)
       Yeval= Y(left)+dt*by(0)
       Zeval= Z(left)+dt*bz(0)

       if (deriv) then
          derivs(1)= bx(0)
          derivs(2)= by(0)
          derivs(3)= bz(0)
       endif

    else ! Parametric cubics

       Xeval= X(left)+dt*(bx(0)+dt*(cx+dt*dx))
       Yeval= Y(left)+dt*(by(0)+dt*(cy+dt*dy))
       Zeval= Z(left)+dt*(bz(0)+dt*(cz+dt*dz))

       if (deriv) then
          derivs(1)= bx(0)+dt*(two*cx+dt*three*dx)
          derivs(2)= by(0)+dt*(two*cy+dt*three*dy)
          derivs(3)= bz(0)+dt*(two*cz+dt*three*dz)
       endif

    endif

  end subroutine plscrv3d

  !===============================================================================
  subroutine lcsfit(ndata,x,y,inew,method,neval,xeval,yeval,ypeval)
  !===============================================================================
    !> Storage-efficient local cubic spline fit (2-space)
    !> (monotonic and piecewise linear options too)
  !===============================================================================
    ! Description and usage:
    !    lcsfit is the non-parametric analog of PLSFIT (parametric).
    ! It is intended for spline applications which do not require the
    ! spline coefficients as output.  It is efficient for repeated
    ! calls with the same data, so repeated use with neval = 1 may be
    ! preferable to storing vectors of results.
    !
    !    lcsfit offers monotonic spline and piecewise linear options
    ! also.  And it returns an interpolated first derivative along
    ! with the function value.  (The second derivative is omitted
    ! because y" is not guaranteed to be continuous by local methods.)
    !
    !    See PLSFIT for more details on local methods.  As with most
    ! numerical methods, scaling of the data to the unit interval (and
    ! unscaling of the result) is recommended to avoid unnecessary
    ! effects attributable to the data units.  Utilities GETSCALE and
    ! USESCALE from the present authors are appropriate.  The data
    ! abscissas should be distinct and either ascending or descending.
    ! PROTECT is available to check this.  Extrapolation is permitted
    ! (mainly in case of round-off; it is normally inadvisable).
    !
    !    The CSFIT/CSEVAL or CSDVAL pair are probably preferable if
    ! efficiency is not an issue, since CSFIT gives y" continuity.
    !
    ! Arguments:
    !     Name    Type/Dimension  I/O/S  Description
    !     ndata   I               I      Length of x, y input data arrays.
    !     x,      R (ndata)       I      Input data coordinates.  The xs
    !     y                              must be distinct and monotonic,
    !                                    either ascending or descending.
    !                                    (No check here.) 
    !     inew     L               I      If control flag inew is .TRUE., the
    !                                    search for a bracket starts from
    !                                    scratch, otherwise locally-saved
    !                                    search and fit information will be
    !                                    assumed to be correct. If calling
    !                                    L!SFIT from within a loop, set
    !                                    inew = .FALSE. after the first call.
    !     method   C*1            I      (Uppercase) Type of fit to be used:
    !                                    'M' means Monotonic piecewise cubics;
    !                                    'B' means non-monotonic "Bessel"-type
    !                                        piecewise cubics (looser fit);
    !                                    'L' means piecewise Linear fit;
    !                                    'C' means Cyclic (periodic) end
    !                                        conditions: loose fit assumed.
    !     neval   I               I      Number of interpolations requested.
    !                                    neval >= 1.  One call per result
    !                                    (neval = 1) may save storage, and is
    !                                    not too inefficient as long as inew
    !                                    is set to .FALSE. after the first.
    !     xeval   R (neval)       I      Abscissa(s) to interpolate to.  These
    !                                    are normally in the data range, but
    !                                    extrapolation - probably due to
    !                                    round-off - is not prevented.
    !     yeval   R (neval)       O      Interpolated function value(s).
    !     ypeval  R (neval)       O      Interpolated 1st derivative value(s).
    !                                    Pass the same storage as for yeval
    !                                    if no derivatives are required.
    !
    ! Significant local variables:
    !     memory         Indicates that coefficients are correct for the
    !                    current point.
    !
    !     h, del         Delta x and forward difference derivative arrays.
    !
    !     B, C, D        Coefficients of cubic on the bracketing interval.
    !
    !     Procedures:
    !     interval  1-D "interpolation" search.
    !     bessel    First derivative (central 3-point formula).
    !     brodlie   First derivative (central), adjusted for monotonicity.
    !     butland   First derivative (non-central), adjusted for monotonicity.
    !     threept   First derivative (non-central 3-point formula).
    !
    ! Notes:
    !     Since many of the calculations must be repeated at both ends
    !     of an interval, the various finite difference quantities used
    !     are stored as arrays. The following "map" of a typical interior
    !     interval and its neighbors should help in understanding the
    !     notation.  The local array indices are all numbered relative
    !     to the left-hand end of the interval which brackets the point
    !     to be evaluated.
    !                                  left       right
    !
    !          Point         -1          0         +1          +2
    !
    !          Data           x -------- x -------- x --------- x
    !
    !          Interval      -1          0         +1
    !
    ! Author: Robert Kennelly, Sterling Software/NASA Ames (PLSFIT)
    !  27 Feb. 1987  R.A.Kennelly  Initial implementation of PLSFIT.
    !  23 Aug. 1989  D.A.Saunders  lcsfit adapted as non-parametric form,
    !                              for embedding in other utilities where
    !                              minimizing work-space is desirable.
    !  20 June 1991    "    "      threept (monotonic) renamed butland;
    !                              threept (pure 3-pt. formula) now used
    !                              for nonmonotonic end-point handling;
    !                              method='C' case belatedly added, as
    !                              needed by PLSINTRP for closed curves.
    !  23 July 1991    "    "      The tests for being in the same interval
    !                              as before were not allowing for the
    !                              descending-xs case.
    !  06 May  1998    "    "      Minor Fortran 90 updates.
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! I/O arguments
    integer, intent(in) :: ndata,neval
    real(wp), intent(in) :: x(ndata),y(ndata),xeval(neval)
    real(wp), intent(out) :: yeval(neval),ypeval(neval)
    logical, intent(in) :: inew
    character, intent(in) :: method*1
    ! ----------------------------------------------------------------------------
    ! Local variables
    integer :: ieval,j,k,left,right
    logical ::  icyclic,linear,memory,mono
    real(wp) :: arrow,B(0:1),C,dely(-1:1),D,dx,h(-1:1),xbyarrow,xe
    ! ----------------------------------------------------------------------------
    ! Storage
    save arrow,B,C,D,left,right
    ! ----------------------------------------------------------------------------

    mono   = method=='M'
    icyclic= method=='C'
    linear = method=='L'

    if (icyclic) then
       if (y(ndata)/=y(1)) stop 'lcsfit: End points must match.'
    endif

    ! Initialize search or avoid it if possible:
    if (inew) then
       memory=.false.
       arrow =sign(one,x(2)-x(1))
       left  =1
    endif

    ieval= 1
    xe= xeval(1)
    xbyarrow= xe*arrow

    if (.not.inew) then
       ! We can save a lot of time when lcsfit is being called from within
       ! a loop by setting memory if possible. The out-of-range checking
       ! relies on the fact that right = left + 1 after the last return.
       ! Cater to the more likely case of xe in the previous, interior
       ! interval.
       memory=xbyarrow>=x(left)*arrow .and. xbyarrow<x(right)*arrow

       if (.not.memory) then
          memory = left==1 .and. xbyarrow<x(right)*arrow .or. &
               right==ndata .and. xbyarrow>=x(left)*arrow
       endif
    endif

    if (memory) GO TO 70 ! Skip the bulk of the computation

    ! Loop over evaluation points requiring a inew search:
    ! ---------------------------------------------------

10  continue

    ! Interpolation search for bracketing interval:
    ! ---------------------------------------------
    call interval(ndata,x,xe,arrow,left)

    right=left+1

    !  -------------------------------------------
    ! |                                           |
    ! |   1 <= left < right = left + 1 <= ndata   |
    ! |                                           |
    !  -------------------------------------------

    ! Compute derivatives by finite-differences:
    ! ------------------------------------------

    if (ndata>2 .and. .not.linear) then

       ! Interval and derivative approximations:
       ! ---------------------------------------

       ! The following duplicates more code than PLSFIT's approach,
       ! but eliminates some indirection - no need to wrap-around here.
       ! Handle the end conditions first to minimize testing left, right.

       if (left==1) then

          h(0)=x(2)-x(1)
          dely(0)=(y(2)-y(1))/h(0)
          h(1)=x(3)-x(2)
          dely(1)=(y(3)-y(2))/h(1)

          if (icyclic) then ! Loose fit assumed
             h(-1)=x(ndata)-x(ndata-1)
             dely(-1)=(y(ndata)-y(ndata-1))/h(-1)
             B(0)=bessel(0,h,dely)
             B(1)=bessel(1,h,dely)
          else
             if (mono) then
                B(0)=butland(0,h,dely)
                B(1)=brodlie(1,h,dely)
             else
                B(0)=threept(0,h,dely)
                B(1)= bessel(1,h,dely)
             endif
          endif

       elseif (right==ndata) then

          h(-1)=x(left)-x(left-1)
          dely(-1)=(y(left)-y(left-1))/h(-1)
          h(0)=x(right)-x(left)
          dely(0)=(y(right)-y(left))/h(0)

          if (icyclic) then
             h(1)=x(2)-x(1)
             dely(1)=(y(2)-y(1))/h(1)
             B(0)=bessel(0,h,dely)
             B(1)=bessel(1,h,dely)
          else

             if (mono) then
                B(0)=brodlie(0,h,dely)
                B(1)=butland(1,h,dely)
             else
                B(0)= bessel(0,h,dely)
                B(1)=threept(1,h,dely)
             endif
          endif

       else
          k=left
          do j=-1,+1
             h(j)=x(k)-x(k-1)
             dely(j)=(y(k)-y(k-1))/h(j)
             k=k+1
          enddo

          ! Select interpolation scheme:
          ! ----------------------------

          ! Compute (possibly adjusted) first derivatives at both
          ! left- and right-hand endpoints of the interval.

          if (mono) then
             ! Monotone - use Brodlie modification of Butland's
             ! formula to adjust the derivatives at the knots.
             B(0)=brodlie(0,h,dely)
             B(1)=brodlie(1,h,dely)
          else ! method = 'B'
             ! Bessel - use central difference formula at the knots.
             B(0)=bessel(0,h,dely)
             B(1)=bessel(1,h,dely)
          endif

       endif

       ! Compute the remaining cubic coefficients.

       C=(three*dely(0)-two*B(0)-B(1))/h(0)
       D=(-two*dely(0)+B(0)+B(1))/h(0)**2

    else ! ndata=2 .or. method='L'

       ! Degenerate case (linear).
       ! -------------------------

       B(0)=(y(right)-y(left))/(x(right)-x(left))
       C=zero
       D=zero

    endif

    ! Evaluate the cubic (derivative first in case only yeval is reqd.):
    ! ------------------------------------------------------------------

70  continue ! Start of same-interval loop inside new-interval loop

    dx=xe-x(left)
    ypeval(ieval)=B(0)+dx*(two*C+dx*three*D)
    yeval(ieval)=y(left)+dx*(B(0)+dx*(C+dx*D))

    ! The next evaluation point may be in the same interval:
    ! ------------------------------------------------------

    if (ieval < neval) then ! Skips this if neval = 1

       ieval=ieval+1
       xe=xeval(ieval)
       xbyarrow=xe*arrow
       if (xbyarrow>=x(left)*arrow .and. &
            xbyarrow< x(right)*arrow) GO TO 70

       GO TO 10 ! Else much more work required.

    endif

  end subroutine lcsfit

  !===============================================================================
  subroutine interval(nx,x,xfind,arrow,left)
  !===============================================================================
    !> Interpolation search for interval containing a point
  !===============================================================================
    ! Description and usage:
    ! Written primarily for interval-based interpolations such as
    ! piecewise linear or cubic spline, interval performs a search to
    ! locate the best interval for evaluating the interpolant at a
    ! given point. The normal case returns the "left-hand" endpoint of
    ! the interval bracketing the point, but for the out-of-range cases
    ! below or above the range of the knots, the interval to be used is
    ! the first or last. The array of knots must be monotonic, either
    ! increasing or decreasing. Diagrammatically, left is returned as
    ! shown below for the normal case (no extrapolation):
    !
    !          x (1)  ...   x (left)   x (left+1)   ...      x (nx)
    !                                ^
    !                              xfind
    !
    !     And for extrapolation:
    !
    !                     x (left = 1)  ...   x (nx)
    !             ^
    !           xfind
    !
    !     or,
    !                x (1)  ...   x (left = nx-1)    x (nx)
    !                                                           ^
    !                                                         xfind
    !
    !     If the point to be bracketed (xfind) matches one of the knots, the
    !     index of that knot is returned as left, i.e., the condition for a
    !     bracket of an interior point is:
    !
    !        x (left) <= xfind < x (left+1)  if  arrow = +1.0,  or
    !        x (left) >= xfind > x (left+1)  if  arrow = -1.0.
    !
    !        This is a low-level routine with minimal error checking. The
    !     calling program is assumed to have verified the following:
    !
    !     (1)  nx >= 2
    !     (2)  x strictly monotonic
    !     (3)  arrow = +1.0 or -1.0
    !
    !     Subroutine PROTECT is available from the author for easily checking
    !     conditions (2) and (3). left is verified on input, but efficiency in
    !     loops will benefit from passing the best estimate available, usually
    !     just the result of the last call.
    !
    !        interval was originally written for use with CSEVAL and TABLE1.
    !     The interpolation search was adapted from ideas in Sedgewick's book
    !     referenced below.
    !
    ! Arguments:
    !     Name  Dimension  Type  I/O/S  Description
    !     nx                I    I      Number of points in array x; must
    !                                   be >= 2 (no check performed).
    !     x        nx       R    I      Array of points defining the set
    !                                   of intervals to be examined. Only
    !                                   the first nx-1 points are required.
    !     xfind             R    I      The point for which a bracketing
    !                                   interval is sought.
    !     arrow             R    I      Monotonicity indicator for input
    !                                   array x:
    !                                     -1.0  strictly decreasing
    !                                      0.0  NOT ALLOWED!
    !                                     +1.0  strictly increasing
    !                                   Supplied by the calling routine for
    !                                   reasons of speed (not checked).
    !     left              I    I/O    Input: guessed index of left-hand
    !                                   endpoint of the interval containing
    !                                   the specified point.
    !
    !                                   Output: index of the largest array
    !                                   value <= specified point (if arrow=+1.0).
    !                                   Special case for data out of range:
    !                                   return left endpoint of closest interval.
    !                                   Thus, left = 1 for xfind < x (2), and
    !                                   left = nx-1 for xfind >= x (nx-1).
    !                                   (If arrow=-1.0, reverse the inequalities.)
    ! Notes:
    !     In speed-critical applications, it might be a good idea to build
    !     this algorithm in-line since it is typically called many times
    !     from within a loop. Another potential speed-up is removal of the
    !     arrow multiplies, which restricts the method to increasing data.
    !     So far, the simplicity of separating out the messy search details
    !     and the generality of bi-directional searching have outweighed
    !     the modest speed penalty incurred.
    !
    ! Bibliography:
    !   Sedgewick, R.  Algorithms.  Reading: Addison-Wesley, 1983 (Chap. 14)
    !
    ! Author:  Robert Kennelly and David Saunders, Sterling Federal Systems
    ! 20 Oct. 1987    RAK    Interpolation search adapted (with mods.
    !                        for bidirectional search and some minor
    !                        repair) from CSEVAL (RAK) and TABLE1 (DAS).
    ! 08 Aug. 1988    DAS    Clarified descriptions of bracketing, where
    !                        the inequalities depend upon arrow.
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! I/O arguments
    integer :: left,nx
    real(wp) :: arrow,x(nx),xfind
    ! ----------------------------------------------------------------------------
    ! Local variables
    integer :: length, nxless1, right, trial
    real(wp) :: xbyarrow
    ! ----------------------------------------------------------------------------

    xbyarrow=xfind*arrow

    ! Simplify things by disposing of two important special cases so that
    ! x (left) and x (right) can really bracket xfind. As a by-product,
    ! this also takes care of the nx = 2, 3 cases.

    nxless1=nx-1

    if (xbyarrow.ge.x(nxless1)*arrow) then
       left= nxless1
       GO TO 990
    else if (xbyarrow.lt.x(2)*arrow) then
       left= 1
       GO TO 990
    endif

    !      -------------------------------
    !     |                               |
    !     |   x (2) <= xfind < x (nx-1)   |
    !     |            - or -             |
    !     |   x (2) > xfind >= x (nx-1)   |
    !     |                               |
    !     |   nx > 3                      |
    !     |                               |
    !      -------------------------------

    ! Adjust the pointers. We hope that the calling routine has provided
    ! a reasonable guess (since it's probably working on an ordered array
    ! of points to evaluate), but check anyway.

    left= min(max(2,left),nx-2)

    if (xbyarrow.ge.x(left)*arrow) then
       if (xbyarrow.lt.x(left+1)*arrow) then
          !  xfind is in the original guessed-at interval.
          GO TO 990
       else
          ! We'll look farther to the right. Evidently left was < nx-2.
          right=nxless1
          left =left+1
       endif
    else
       !        Look to the left of the guess. Evidently left was > 2.
       right=left
       left =2
    endif

    !      --------------------------------
    !     |                                |
    !     |   2 <= left < right <= nx-1    |
    !     |                                |
    !      --------------------------------

    ! The interval length must decrease each time through - terminate
    ! when the correct interval is found or when the interval length
    ! cannot be decreased.

10  continue

    length=right-left
    if (length.gt.1) then
       ! The trial value is a "linear" estimate of the left-hand endpoint
       ! of the interval bracketing the target xfind, with protection
       ! against round-off (which can affect convergence).
       trial=min(right-1,left+max(0,int(dble(length)*(xfind-x(left))/(x(right)-x(left)))))

       !    ----------------------------------------
       !   |                                        |
       !   |   2 <= left <= trial < right <= nx-1   |
       !   |                                        |
       !    ----------------------------------------

       !  Adjust pointers. Increase left or decrease right until done.
       if (xbyarrow.ge.x(trial+1)*arrow) then
          left =trial+1
       else if (xbyarrow.lt.x(trial)*arrow) then
          right=trial
       else
          ! We're done: xfind is in the interval [x (trial), x (trial+1)).
          left =trial
          GO TO 990
       endif
       GO TO 10

    endif

990 continue

  end subroutine interval

  !===============================================================================
  function bessel(j,h,del)
  !===============================================================================
    !> First derivative using central 3-point formula
  !===============================================================================
    ! Description and usage:
    ! Computes a first derivative approximation using the central
    ! 3-point formula.  The data must be in the form of arrays containing
    ! finite difference interval lengths and 2-point forward difference
    ! derivatives.  bessel is intended to be used by PLSFIT for determin-
    ! ing end conditions on an interval for (non-monotonic) interpolation
    ! by piecewise cubics.  See the PLSFIT header for more details.
    !
    ! Arguments:
    !     Name    Type/Dimension  I/O/S  Description
    !     j       I               I      Indicates at which end of the
    !                                    interval the derivative is to be
    !                                    estimated. j = 0 means left-hand
    !                                    side, j = 1 means right.
    !     h       R (-1:1)        I      Array of interval lengths. The 0th
    !                                    element is the length of the interval
    !                                    on which the cubic is to be deter-
    !                                    mined.
    !     del     R (-1:1)        I      Array of derivative estimates. The
    !                                    0th element is the forward difference
    !                                    derivative over the interval on which
    !                                    the cubic is to be determined.
    !     bessel  R                 O    The function value is the adjusted
    !                                    derivative.
    !
    ! Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
    ! 18 Feb. 1987    RAK    Initial design and coding.
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! I/O arguments
    integer :: j
    real(wp) :: h(-1:1),del(-1:1),bessel
    ! ----------------------------------------------------------------------------
    ! Local variables
    real(wp) :: weight
    ! ----------------------------------------------------------------------------

    ! Estimate first derivative on left (j=0) or right side (j=1) of an interval
    weight= h(j) / (h(j)+h(j-1))
    bessel= weight*del(j-1) + (one-weight)*del(j)

  end function bessel

  !===============================================================================
  function butland(j,h,del)
  !===============================================================================
    !> First derivative, non-central 3-point formula, adjusted
  !===============================================================================
    ! Description and usage:
    ! Computes a first derivative approximation for PLSFIT over an
    ! interval at a data boundary, using a modified forward or backward
    ! 3-point formula.  The data must be in the form of arrays containing
    ! finite difference interval lengths and 2-point forward difference
    ! derivatives, and the differencing direction is controlled by a flag.
    ! See PLSFIT for more details, or THREEPT for the pure 3-pt. formula.
    !
    ! The "shape preserving adjustments" are from PCHIP, a monotone
    ! piecewise cubic interpolation package by F. N. Fritsch.
    !
    ! Arguments:
    !     Name    Type/Dimension  I/O/S  Description
    !     j       I               I      Indicates at which end of the
    !                                    interval the derivative is to be
    !                                    estimated. J = 0 means left-hand
    !                                    side, J = 1 means right.
    !     h       R (-1:1)        I      Array of interval lengths. The 0th
    !                                    element is the length of the interval
    !                                    on which the cubic is to be deter-
    !                                    mined.
    !     del     R (-1:1)        I      Array of derivative estimates. The
    !                                    0th element is the forward difference
    !                                    derivative over the interval on which
    !                                    the cubic is to be determined.
    !     butland R                 O    The function value is the adjusted
    !                                    derivative.
    !  Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
    !  18 Feb. 1987    RAK    Initial design and coding, as THREEPT.
    !  20 June 1991    DAS    Monotonic form renamed BUTLAND; THREEPT
    !                         is now the pure 3-point formula.
    !  04 Dec. 2002     "     SIGN work-arounds for Intel compiler bug.
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! I/O arguments
    integer :: j
    real(wp) :: h(-1:1),del(-1:1),butland
    ! ----------------------------------------------------------------------------
    ! Local variables
    integer :: step
    real(wp) :: dmax,weight
    logical :: constrain
    ! ----------------------------------------------------------------------------

    ! Estimate first derivative on a left-hand boundary using a 3-point
    ! forward difference (step = +1), or with a backward difference for
    ! the right-hand boundary (step = -1).
    step = 1-j-j   ! j here is consistent with related modules.

    ! In {h,del} form, the derivative looks like a weighted average.

    !***  Avoid zero as the second argument of SIGN: Intel compiler misbehaves.

    if (del(0)==zero) then

       butland=zero

    else ! butland below cannot be zero

       weight=-h(0)/(h(0)+h(step))
       butland=weight*del(step) + (one-weight)*del(0)

       ! Shape-preserving adjustments.  Note that we try to avoid overflow
       ! by not multiplying quantities directly.

       if (sign (one,butland)/=sign(one,del(0))) then

          ! Defer to the estimate closest to the boundary.
          butland = zero

          !******  elseif (sign(one,del(0)).ne.sign(one,del(step))) then
       else

          if (del(step)==zero) then
             constrain=del(0)<zero
          else
             constrain=sign(one,del(0))/=sign(one,del(step))
          endif

          if (constrain) then
             ! If the monotonicity switches, may need to bound the estimate.
             dmax= three*del(0)
             if (abs(butland)>abs(dmax)) butland= dmax
          endif

       endif

    endif

  end function butland

  !===============================================================================
  function brodlie(j,h,del)
  !===============================================================================
    !> First derivative, adjusted for monotonicity
  !===============================================================================
    ! Description and usage:
    ! brodlie is intended to be used by PLSFIT for determining end
    ! conditions on an interval for monotonic interpolation by piecewise
    ! cubics. The data must be in the form of arrays containing finite
    ! difference interval lengths and 2-point forward difference deriva-
    ! tives. See the PLSFIT header for more details.
    !
    !        The method is due to Brodlie, Butland, Carlson, and Fritsch,
    !     as referenced in the PLSFIT header.
    !
    ! Arguments:
    !     Name    Type/Dimension  I/O/S  Description
    !     j       I               I      Indicates at which end of the
    !                                    interval the derivative is to be
    !                                    estimated. j = 0 means left-hand
    !                                    side, j = 1 means right.
    !     h       R (-1:1)        I      Array of interval lengths. The 0th
    !                                    element is the length of the interval
    !                                    on which the cubic is to be determined.
    !     del     R (-1:1)        I      Array of derivative estimates. The
    !                                    0th element is the forward difference
    !                                    derivative over the interval on which
    !                                    the cubic is to be determined.
    !     brodlie R                 O    The function value is the adjusted
    !                                    derivative.
    ! Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
    ! 18 Feb. 1987    RAK    Initial design and coding.
    ! 04 Dec. 2002    DAS    SIGN work-around for Intel compiler bug.
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! I/O arguments
    integer :: j
    real(wp) :: brodlie,h(-1:1),del(-1:1)
    ! ----------------------------------------------------------------------------
    ! Local variables
    real(wp) :: alpha, product
    ! ----------------------------------------------------------------------------

    ! Compare the algebraic signs of the two del's.  Have to test that
    ! at least one is positive to avoid a zero denominator (this fancy
    ! test permits one term to be zero, but the answer below is zero
    ! anyway in these cases).  The trick is to work around the SIGN
    ! function, which returns positive even if its 2nd argument is zero.

    !**** NO:  SIGN misbehaves on Intel systems when the 2nd argument is zero.

    product=del(j-1)*del(j)

    if (product==zero) then
       brodlie=zero
    elseif (sign (one, -del (j - 1)) .ne. sign (one, del (j))) then
       ! form "weighted harmonic mean" of the two finite-difference
       ! derivative approximations.  Note that we try to avoid overflow
       ! by not multiplying them together directly.
       alpha= third*(one + h(j)/(h(j-1)+h(j)))
       brodlie= product/(alpha*del(j)+(one-alpha)*del(j-1))
    else
       ! the signs differ, so make this point a local extremum.
       brodlie=zero
    endif

  end function brodlie

  !===============================================================================
  function threept(j,h,del)
  !===============================================================================
    !> First derivative, non-central 3-point formula
  !===============================================================================
    ! Description and usage:
    ! Computes a first derivative approximation for PLSFIT over an
    ! interval at a data boundary, using a forward or backward 3-point
    ! formula.  The data must be in the form of arrays containing finite
    ! difference interval lengths and 2-point forward difference deriva-
    ! tives, and the differencing direction is controlled by a flag. See
    ! PLSFIT for more details.
    !
    ! See module BUTLAND for a version with "shape-preserving" adjustments.
    !
    !     Arguments:
    !     ----------
    !
    !     Name    Type/Dimension  I/O/S  Description
    !     j       I               I      Indicates at which end of the
    !                                    interval the derivative is to be
    !                                    estimated. j = 0 means left-hand
    !                                    side, j = 1 means right. 
    !     h       R (-1:1)        I      Array of interval lengths. The 0th
    !                                    element is the length of the interval
    !                                    on which the cubic is to be deter-
    !                                    mined.
    !     del     R (-1:1)        I      Array of derivative estimates. The
    !                                    0th element is the forward difference
    !                                    derivative over the interval on which
    !                                    the cubic is to be determined.
    !     threept R                 O    The function value is the derivative.
    !
    ! Author:  Robert Kennelly, Sterling Federal Systems/NASA-Ames
    ! 18 Feb. 1987    RAK    Initial design and coding.
    ! 06 June 1991    DAS    Original threept renamed BUTLAND; threept
    !                        now gives unmodified 1-sided 3-pt. results.
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! I/O arguments
    integer :: j
    real(wp) :: h(-1:1),del(-1:1),threept
    ! ----------------------------------------------------------------------------
    ! Local variables
    integer :: step
    real(wp) :: weight
    ! ----------------------------------------------------------------------------

    ! Estimate first derivative on a left-hand boundary using a 3-point
    ! forward difference (step = +1), or with a backward difference for
    ! the right-hand boundary (step = -1).
    step = 1-j-j   ! j here is consistent with related modules.

    ! In {h,del} form, the derivative looks like a weighted average.

    weight= -h(0)/(h(0)+h(step))
    threept= weight*del(step) + (one-weight)*del(0)

  end function threept

  !===============================================================================
  subroutine warp2d(imin,imax,jmin,jmax,i1,i2,j1,j2,x0,y0,s0,dedgei,dedgej,x,y)
  !===============================================================================
    !> Perturb interior of a plane grid given new edges
  !===============================================================================
    ! Description:
    !       WARP2D perturbs the interior points of a 2-space grid block
    !    given the original grid and perturbed edges. All four edges are
    !    assumed to be perturbed (though some may of course be fixed).
    !    If all corners are found to be unperturbed, considerably less
    !    work is performed.  In general, an intermediate perturbation is
    !    required to account for the corner motion, using edges derived
    !    from the original edges, then a second perturbation accounts
    !    for differences between the intermediate edges and the specified
    !    new edges.
    !      The perturbed edges should be input as edges of the desired
    !    output grid.  The original relative arc-length increments in
    !    each index direction should also be input.  See PARAMxy for
    !    setting them up in preparation for multiple perturbations.
    !
    ! Author: David Saunders/James Reuther, NASA Ames, Mt. View, CA.
    ! ----------------------------------------------------------------------------
    implicit none
    ! ----------------------------------------------------------------------------
    ! I/O arguments:
    ! --------------
    !  I  Grid array dimensions.
    integer :: imin,imax,jmin,jmax
    !  I  Active area is (i1:i2,j1:j2).
    integer :: i1,i2,j1,j2
    !  I  Original grid coordinates.
    real(wp) :: x0(imin:imax,jmin:jmax),y0(imin:imax,jmin:jmax)
    !  I  Relative arc-length increments for all lines in the i and j directions
    ! (see paramxy). If from a grid larger than the active subgrid, S0 is transformed here.
    real(wp) :: s0(imin:imax,jmin:jmax,2)
    !  S  Storage for edge perturbations: dedgei (1:2,j1:j2,1)=dx,dy along the i=i1 edge, etc.
    real(wp) :: dedgei(2,jmin:jmax,2),dedgej(2,imin:imax,2)
    ! I/O Grid coordinates: input with the edges perturbed; output fully perturbed.
    real(wp) :: x(imin:imax,jmin:jmax),y(imin:imax,jmin:jmax)
    ! ----------------------------------------------------------------------------
    ! Local variables:
    ! ----------------
    integer :: i,j
    ! eps safeguards a divide by zero (presumably only if result is zero)
    real(wp), parameter :: eps=1.0e-8_wp,one=1.0_wp
    real(wp) :: deli,delj,wti1,wti2,wtj1,wtj2
    logical :: same
    ! ----------------------------------------------------------------------------

    ! Have any corners moved?  Set up the perturbations in available
    ! storage because they are needed anyway if the corners have moved.
    same=.true.
    i=i1

    do j=1,2
       dedgei(1,j1,j)=x(i,j1)-x0(i,j1)
       dedgei(2,j1,j)=y(i,j1)-y0(i,j1)
       dedgei(1,j2,j)=x(i,j2)-x0(i,j2)
       dedgei(2,j2,j)=y(i,j2)-y0(i,j2)

       same=(dedgei(1,j1,j).eq.0.).and.(dedgei(2,j1,j).eq.0.).and. &
            (dedgei(1,j2,j).eq.0.).and.(dedgei(2,j2,j).eq.0.).and.same
       i=i2
    enddo

    if (.not.same) then  ! Corners have moved

       ! Set up intermediate edges with the final corners but
       ! otherwise derived from the original edges.  Actually,
       ! all we need are the intermediate edge perturbations.

       ! i=i1 and i2 intermediate edges:

       do j=j1+1,j2-1
          wtj2=(s0(i1,j, 2)-s0(i1,j1,2))/(s0(i1,j2,2)-s0(i1,j1,2))
          wtj1=one-wtj2
          dedgei(1,j,1)=wtj1*dedgei(1,j1,1)+wtj2*dedgei(1,j2,1)
          dedgei(2,j,1)=wtj1*dedgei(2,j1,1)+wtj2*dedgei(2,j2,1)

          wtj2=(s0(i2,j, 2)-s0(i2,j1,2))/(s0(i2,j2,2)-s0(i2,j1,2))
          wtj1=one-wtj2
          dedgei(1,j,2)=wtj1*dedgei(1,j1,2)+wtj2*dedgei(1,j2,2)
          dedgei(2,j,2)=wtj1*dedgei(2,j1,2)+wtj2*dedgei(2,j2,2)
       enddo

       ! j=j1 and j2 intermediate edges:

       do i=i1+1,i2-1
          wti2=(s0(i, j1,1)-s0(i1,j1,1))/(s0(i2,j1,1)-s0(i1,j1,1))
          wti1=one-wti2
          dedgej(1,i,1)=wti1*dedgei(1,j1,1)+wti2*dedgei(1,j1,2)
          dedgej(2,i,1)=wti1*dedgei(2,j1,1)+wti2*dedgei(2,j1,2)

          wti2=(s0(i, j2,1)-s0(i1,j2,1))/(s0(i2,j2,1)-s0(i1,j2,1))
          wti1=one-wti2
          dedgej(1,i,2)=wti1*dedgei(1,j2,1)+wti2*dedgei(1,j2,2)
          dedgej(2,i,2)=wti1*dedgei(2,j2,1)+wti2*dedgei(2,j2,2)
       enddo


       ! Perturb the interior based on the original interior and
       ! compatible intermediate edges connecting the moved corners.
       ! The contributions from each index direction are NOT
       ! independent,so they are combined as a weighted average.

       do j=j1+1,j2-1
          do i=i1+1,i2-1
             wti2=(s0(i, j,1)-s0(i1,j,1))/(s0(i2,j,1)-s0(i1,j,1))
             wti1=one-wti2
             wtj2=(s0(i,j, 2)-s0(i,j1,2))/(s0(i,j2,2)-s0(i,j1,2))
             wtj1=one-wtj2

             deli=wti1*dedgei(1,j,1)+wti2*dedgei(1,j,2)
             delj=wtj1*dedgej(1,i,1)+wtj2*dedgej(1,i,2)

             x(i,j)=x0(i,j)+(abs(deli)*deli+abs(delj)*delj) &
                  / max(abs(deli)+abs(delj),eps)

             deli=wti1*dedgei(2,j,1)+wti2*dedgei(2,j,2)
             delj=wtj1*dedgej(2,i,1)+wtj2*dedgej(2,i,2)

             y(i,j)=y0(i,j)+(abs(deli)*deli+abs(delj)*delj) &
                  / max(abs(deli)+abs(delj),eps)
          enddo
       enddo

       ! Set up the second-phase edge perturbations from the
       ! intermediate edges to the final edges:

       ! i=i1 and i2 edges:

       do j=j1+1,j2-1
          dedgei(1,j,1)=x(i1,j)-x0(i1,j)-dedgei(1,j,1)
          dedgei(2,j,1)=y(i1,j)-y0(i1,j)-dedgei(2,j,1)
          dedgei(1,j,2)=x(i2,j)-x0(i2,j)-dedgei(1,j,2)
          dedgei(2,j,2)=y(i2,j)-y0(i2,j)-dedgei(2,j,2)
       enddo

       ! j=j1 and j2 edges:

       do i=i1+1,i2-1
          dedgej(1,i,1)=x(i,j1)-x0(i,j1)-dedgej(1,i,1)
          dedgej(2,i,1)=y(i,j1)-y0(i,j1)-dedgej(2,i,1)
          dedgej(1,i,2)=x(i,j2)-x0(i,j2)-dedgej(1,i,2)
          dedgej(2,i,2)=y(i,j2)-y0(i,j2)-dedgej(2,i,2)
       enddo

    else

       ! just set up the edge perturbations for a single phase.

       ! i=i1 and i2 edges:

       do j=j1+1,j2-1
          dedgei(1,j,1)=x(i1,j)-x0(i1,j)
          dedgei(2,j,1)=y(i1,j)-y0(i1,j)
          dedgei(1,j,2)=x(i2,j)-x0(i2,j)
          dedgei(2,j,2)=y(i2,j)-y0(i2,j)
       enddo

       ! j=j1 and j2 edges:

       do i=i1+1,i2-1
          dedgej(1,i,1)=x(i,j1)-x0(i,j1)
          dedgej(2,i,1)=y(i,j1)-y0(i,j1)
          dedgej(1,i,2)=x(i,j2)-x0(i,j2)
          dedgej(2,i,2)=y(i,j2)-y0(i,j2)
       enddo

       ! Transfer the original interior to avoid duplicating code below:

       do j=j1+1,j2-1
          do i=i1+1,i2-1
             x(i,j)=x0(i,j)
             y(i,j)=y0(i,j)
          enddo
       enddo

    endif

    ! Update the interior with contributions from all edges.  The two
    ! directions are independent,so just accumulate the contributions.

    do j=j1+1,j2-1
       do i=i1+1,i2-1
          wti2=(s0(i, j,1)-s0(i1,j,1))/(s0(i2,j,1)-s0(i1,j,1))
          wti1=one-wti2
          wtj2=(s0(i,j, 2)-s0(i,j1,2))/(s0(i,j2,2)-s0(i,j1,2))
          wtj1=one-wtj2

          deli=wti1*dedgei(1,j,1)+wti2*dedgei(1,j,2)
          delj=wtj1*dedgej(1,i,1)+wtj2*dedgej(1,i,2)

          x(i,j)=x(i,j)+deli+delj

          deli=wti1*dedgei(2,j,1)+wti2*dedgei(2,j,2)
          delj=wtj1*dedgej(2,i,1)+wtj2*dedgej(2,i,2)

          y(i,j)=y(i,j)+deli+delj
       enddo
    enddo

  end subroutine warp2d

  !===============================================================================
  subroutine paramxy(imin,imax,jmin,jmax,i1,i2,j1,j2,x,y,s)
  !===============================================================================
    !> Relative arc-lengths for all lines of a 2-space grid
  !===============================================================================
    ! Description:
    !    PARAMXY parameterizes a regular 2-space grid by setting up
    !    the normalized arc-length increments in both index directions.
    !
    ! Author: David Saunders, Sterling Software/NASA Ames, Mt. View, CA.
    ! ----------------------------------------------------------------------------
    implicit none
    ! ----------------------------------------------------------------------------
    ! I/O arguments:
    ! --------------
    !  I  Grid array dimensions.
    integer :: imin,imax,jmin,jmax
    !  I  Active area is (i1:i2,j1:j2).
    integer :: i1,i2,j1,j2
    !  I  Grid coordinates.
    real(wp), dimension(imin:imax,jmin:jmax) :: x,y
    !  O  Relative arc-length increments for lines in the i and j directions
    !     with S(i1,j,1)=S(i,j1,2)=0. and S(i2,j,1)=S(i,j2,2)=1.
    real(wp), dimension(imin:imax,jmin:jmax,2) :: s
    ! ----------------------------------------------------------------------------
    ! Local variables:
    ! ----------------
    integer :: i,j
    real(wp), dimension(imin:imax,jmin:jmax) :: deli,delj
    ! ----------------------------------------------------------------------------

    do i=i1+1,i2
       do j=j1,j2
          deli(i,j)=sqrt((x(i,j)-x(i-1,j))**2+(y(i,j)-y(i-1,j))**2)
       enddo
    enddo

    do i=i1,i2
       do j=j1+1,j2
          delj(i,j)=sqrt((x(i,j)-x(i,j-1))**2+(y(i,j)-y(i,j-1))**2)
       enddo
    enddo

    ! Zero the two low-end edges.
    do j=j1,j2
       s(i1,j,1)=0.0_wp
    enddo

    do i=i1,i2
       s(i,j1,2)=0.0_wp
    enddo

    ! set up the low-end edge lines because they are missed by the
    ! following loops over most of the interior:

    do i=i1+1,i2
       s(i,j1,1)=s(i-1,j1,1)+deli(i,j1)
    enddo

    do j=j1+1,j2
       s(i1,j,2)=s(i1,j-1,2)+delj(i1,j)
    enddo

    ! Traverse the grid just once for all lines except those within the low-end edges.

    do j=j1+1,j2
       do i=i1+1,i2
          s(i,j,1)=s(i-1,j,1)+deli(i,j)
          s(i,j,2)=s(i,j-1,2)+delj(i,j)
       enddo
    enddo

    ! Normalizing requires another pass through the grid:

    do j=j1,j2
       do i=i1,i2
          s(i,j,1)=s(i,j,1)/s(i2,j,1)
          s(i,j,2)=s(i,j,2)/s(i,j2,2)
       enddo
    enddo

    ! Finally,precise 1s for the two high-end edges:

    do j=j1,j2
       s(i2,j,1)=1.0_wp
    enddo

    do i=i1,i2
       s(i,j2,2)=1.0_wp
    enddo

  end subroutine paramxy

  !===============================================================================
  subroutine upcase(STRING)
  !===============================================================================
    !> changes all lower case letters in the given character STRING to upper case
  !===============================================================================
    ! Description:
    !    Each character in STRING is treated in turn. The intrinsic
    !    function INDEX effectively allows a table lookup, with
    !    the local strings LOW and UPP acting as two tables.
    !    This method avoids the use of CHAR and ICHAR, which appear
    !    to behave differently on ASCII and EBCDIC machines.
    !
    ! Author: Michael Saunders, Systems Optimization Lab., Stanford, 9/10/85
    !         (rewrite of version by Hooper/Kennelly, Informatics, 1983)
    ! ----------------------------------------------------------------------------
    implicit none
    ! ----------------------------------------------------------------------------
    ! I/O arguments:
    ! --------------
    ! character string possibly containing some lower-case letters on input;
    ! strictly upper-case letters on output with no change to any non-alphabetic characters.
    character(len=*), intent(inout) :: STRING(*)
    ! ----------------------------------------------------------------------------
    ! Local variables:
    ! ----------------
    integer :: i,j
    character :: c
    character(len=26) :: low,upp
    ! ----------------------------------------------------------------------------

    low='abcdefghijklmnopqrstuvwxyz'
    upp='ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    do j=1,len(STRING)
       !c=STRING(j:j)
       c=STRING(j)
       if ((c.ge.'a').and.(c.le.'z')) then
          i=index(low,c)
          if (i.gt.0) STRING(j:j)=upp(i:i)
          !if (i.gt.0) STRING(j)=upp(i)
       endif
    enddo

 end subroutine upcase

end module mod_grid_utilities
