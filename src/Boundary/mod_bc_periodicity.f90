!================================================================================
module mod_bc_periodicity
!================================================================================
  !> Module to manage periodicity Boundary Conditions
!================================================================================
  use mod_flow ! for: flow variables
  use mod_mpi  ! for: coord
  implicit none
  !------------------------------------------------------------------------------
  ! boolean indicator for periodicity
  logical :: is_periodic(6) ! periodic BC
  ! direction where periodicity is applied
  integer :: period_dir(6)
  ! translation vector for planar periodicity
  real(wp) :: Lxp,Lyp,Lzp
  ! angle for angular periodicity
  real(wp) :: theta_period
  ! sin/cos of periodicity angle
  real(wp) :: costp,sintp(6)
  !------------------------------------------------------------------------------

contains

  !==============================================================================
  subroutine init_periodicity
  !==============================================================================
    !> Init periodicity Boundary Conditions
    !> --> determine periodic BCs
    !> --> init angular periodicity (sin/cos + rotation orientation)
    !> [ called in mod_mpi_part ]
  !==============================================================================
    use mod_mpi      ! for: iproc, coord
    use mod_constant ! for: pi
    use mod_block    ! for: block informations
    use warnstop
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,nbl,nbv
    real(wp) :: d_o,d_rp,d_rm
    real(wp) :: x_r,y_r,z_r,x_n,y_n,z_n
    real(wp) :: x_rp,y_rp,z_rp,x_rm,y_rm,z_rm
    ! ---------------------------------------------------------------------------
    
    ! Current block
    ! =============
    nbl=nob(iproc)

    ! Determine periodic faces
    ! ========================
    ! initializations
    is_periodic=.false.
    bl(nbl)%periodic = .false.

    ! use flag set in param_blocks.ini for block attribute
    if (bl(nbl)%flag(1)=='p') then
       bl(nbl)%periodic(1) = .true.
    endif
    if (bl(nbl)%flag(2)=='p') then
       bl(nbl)%periodic(2) = .true.
    endif
    if (bl(nbl)%flag(3)=='p') then
       bl(nbl)%periodic(3) = .true.
    endif
    if (bl(nbl)%flag(4)=='p') then
       bl(nbl)%periodic(4) = .true.
    endif
    if (bl(nbl)%flag(5)=='p') then
       bl(nbl)%periodic(5) = .true.
    endif
    if (bl(nbl)%flag(6)=='p') then
       bl(nbl)%periodic(6) = .true.
    endif


    ! use flag set in param_blocks.ini for local periodicity
    if ((coord(1)==0).and.(bl(nbl)%flag(1)=='p')) then
       is_periodic(1)=.true.
       period_dir(1)=1
    endif
    if ((coord(1)==ndomx-1).and.(bl(nbl)%flag(2)=='p')) then
       is_periodic(2)=.true.
       period_dir(2)=1
    endif
    if ((coord(2)==0).and.(bl(nbl)%flag(3)=='p')) then
       is_periodic(3)=.true.
       period_dir(3)=2
    endif
    if ((coord(2)==ndomy-1).and.(bl(nbl)%flag(4)=='p')) then
       is_periodic(4)=.true.
       period_dir(4)=2
    endif
    if ((coord(3)==0).and.(bl(nbl)%flag(5)=='p')) then
       is_periodic(5)=.true.
       period_dir(5)=3
    endif
    if ((coord(3)==ndomz-1).and.(bl(nbl)%flag(6)=='p')) then
       is_periodic(6)=.true.
       period_dir(6)=3
    endif
     
    if (((Lxp.ne.0.0_wp).or.(Lyp.ne.0.0_wp).or.(Lzp.ne.0.0_wp)) &
         .and.(theta_period==0.0_wp).and.(is_curv3).and.(idepart==1)) then
    !     ============================
    ! --> init translation periodicity
    !     ============================
       
       ! determine sign of translation
       ! ==========================
       
       ! based on block corners
       ! ----------------------
       call share_block_corners

       ! imin face
       ! ---------
       if (is_periodic(1)) then

          ! reference point in block
          x_r=bl(nbl)%xcorner(1)
          y_r=bl(nbl)%ycorner(1)
          z_r=bl(nbl)%zcorner(1)
          ! corresponding point in neighboring block
          nbv=bl(nbl)%BC(1)      
          x_n=bl(nbv)%xcorner(2)
          y_n=bl(nbv)%ycorner(2)
          z_n=bl(nbv)%zcorner(2)

          ! original distance
          d_o=sqrt((x_r-x_n)**2+(y_r-y_n)**2+(z_r-z_n)**2)

          ! apply translation with plus sign
          x_rp=x_r+Lxp
          y_rp=y_r+Lyp
          z_rp=z_r+Lzp
          d_rp=sqrt((x_rp-x_n)**2+(y_rp-y_n)**2+(z_rp-z_n)**2)

          ! apply translation with minus sign
          x_rm=x_r-Lxp
          y_rm=y_r-Lyp
          z_rm=z_r-Lzp
          d_rm=sqrt((x_rm-x_n)**2+(y_rm-y_n)**2+(z_rm-z_n)**2)

          ! change sign of translation
          if ((d_rp/d_o<1).and.(d_rm/d_o>1)) then
             Lxp=-Lxp
             Lyp=-Lyp
             Lzp=-Lzp
          endif

          !print *,'imin: block',nbl,'neighbor',nbv
          !print *,'dist',d_rp/d_o,d_rm/d_o
       endif
       
       ! imax face
       ! ---------
       if (is_periodic(2)) then

          ! reference point in block
          x_r=bl(nbl)%xcorner(2)
          y_r=bl(nbl)%ycorner(2)
          z_r=bl(nbl)%zcorner(2)
          ! corresponding point in neighboring block
          nbv=bl(nbl)%BC(2)      
          x_n=bl(nbv)%xcorner(1)
          y_n=bl(nbv)%ycorner(1)
          z_n=bl(nbv)%zcorner(1)

          ! original distance
          d_o=sqrt((x_r-x_n)**2+(y_r-y_n)**2+(z_r-z_n)**2)

          ! apply translation with plus sign
          x_rp=x_r+Lxp
          y_rp=y_r+Lyp
          z_rp=z_r+Lzp
          d_rp=sqrt((x_rp-x_n)**2+(y_rp-y_n)**2+(z_rp-z_n)**2)

          ! apply translation with minus sign
          x_rm=x_r-Lxp
          y_rm=y_r-Lyp
          z_rm=z_r-Lzp
          d_rm=sqrt((x_rm-x_n)**2+(y_rm-y_n)**2+(z_rm-z_n)**2)

          ! change sign of translation
          if ((d_rp/d_o<1).and.(d_rm/d_o>1)) then
             Lxp=-Lxp
             Lyp=-Lyp
             Lzp=-Lzp
          endif

          !print *,'imax: block',nbl,'neighbor',nbv
          !print *,'dist',d_rp/d_o,d_rm/d_o
       endif

       ! jmin face
       ! ---------
       if (is_periodic(3)) then

          ! reference point in block
          x_r=bl(nbl)%xcorner(1)
          y_r=bl(nbl)%ycorner(1)
          z_r=bl(nbl)%zcorner(1)
          ! corresponding point in neighboring block
          nbv=bl(nbl)%BC(3)      
          x_n=bl(nbv)%xcorner(3)
          y_n=bl(nbv)%ycorner(3)
          z_n=bl(nbv)%zcorner(3)

          ! original distance
          d_o=sqrt((x_r-x_n)**2+(y_r-y_n)**2+(z_r-z_n)**2)

          ! apply translation with plus sign
          x_rp=x_r+Lxp
          y_rp=y_r+Lyp
          z_rp=z_r+Lzp
          d_rp=sqrt((x_rp-x_n)**2+(y_rp-y_n)**2+(z_rp-z_n)**2)

          ! apply translation with minus sign
          x_rm=x_r-Lxp
          y_rm=y_r-Lyp
          z_rm=z_r-Lzp
          d_rm=sqrt((x_rm-x_n)**2+(y_rm-y_n)**2+(z_rm-z_n)**2)

          ! change sign of translation
          if ((d_rp/d_o<1).and.(d_rm/d_o>1)) then
             Lxp=-Lxp
             Lyp=-Lyp
             Lzp=-Lzp
          endif

          !print *,'jmin: block',nbl,'neighbor',nbv
          !print *,'dist',d_rp/d_o,d_rm/d_o
       endif

       ! jmax face
       ! ---------
       if (is_periodic(4)) then

          ! reference point in block
          x_r=bl(nbl)%xcorner(3)
          y_r=bl(nbl)%ycorner(3)
          z_r=bl(nbl)%zcorner(3)
          ! corresponding point in neighboring block
          nbv=bl(nbl)%BC(4)      
          x_n=bl(nbv)%xcorner(1)
          y_n=bl(nbv)%ycorner(1)
          z_n=bl(nbv)%zcorner(1)

          ! original distance
          d_o=sqrt((x_r-x_n)**2+(y_r-y_n)**2+(z_r-z_n)**2)

          ! apply translation with plus sign
          x_rp=x_r+Lxp
          y_rp=y_r+Lyp
          z_rp=z_r+Lzp
          d_rp=sqrt((x_rp-x_n)**2+(y_rp-y_n)**2+(z_rp-z_n)**2)

          ! apply translation with minus sign
          x_rm=x_r-Lxp
          y_rm=y_r-Lyp
          z_rm=z_r-Lzp
          d_rm=sqrt((x_rm-x_n)**2+(y_rm-y_n)**2+(z_rm-z_n)**2)

          ! change sign of translation
          if ((d_rp/d_o<1).and.(d_rm/d_o>1)) then
             Lxp=-Lxp
             Lyp=-Lyp
             Lzp=-Lzp
          endif

          !print *,'jmax: block',nbl,'neighbor',nbv
          !print *,'dist',d_rp/d_o,d_rm/d_o
       endif
    
    endif
   
    if (theta_period==0.0_wp) return
    !     ========================
    ! --> init angular periodicity
    !     ========================
    
    ! sine and cosine of periodicity angle
    ! ====================================
    costp=cos(theta_period*pi/180.0_wp)
    sintp=sin(theta_period*pi/180.0_wp)
    
    ! determine sign of rotation
    ! ==========================
    
    ! based on block corners
    ! ----------------------
    call share_block_corners
    
    ! imin face
    ! ---------
    if (is_periodic(1)) then
       
       ! reference point in block
       x_r=bl(nbl)%xcorner(1)
       y_r=bl(nbl)%ycorner(1)
       z_r=bl(nbl)%zcorner(1)
       ! corresponding point in neighboring block
       nbv=bl(nbl)%BC(1)      
       x_n=bl(nbv)%xcorner(2)
       y_n=bl(nbv)%ycorner(2)
       z_n=bl(nbv)%zcorner(2)

       ! original distance
       d_o=sqrt((x_r-x_n)**2+(y_r-y_n)**2+(z_r-z_n)**2)
       
       ! apply rotation with plus sign
       y_rp= costp*y_r-sintp(1)*z_r
       z_rp= sintp(1)*y_r+costp*z_r
       d_rp=sqrt((x_r-x_n)**2+(y_rp-y_n)**2+(z_rp-z_n)**2)

       ! apply rotation with minus sign
       y_rm= costp*y_r+sintp(1)*z_r
       z_rm=-sintp(1)*y_r+costp*z_r
       d_rm=sqrt((x_r-x_n)**2+(y_rm-y_n)**2+(z_rm-z_n)**2)

       ! change sign of sin(theta_period) if negative rotation
       !if ((d_rp/d_o>1).and.(d_rm/d_o<1)) sintp(1)=-sintp(1)
       if ((d_rp/d_o<1).and.(d_rm/d_o>1)) then
          sintp(1)=-sintp(1)
          print *,'imin: change sign',nbl,iproc
       endif
       
       !print *,'imin: block',nbl,'neighbor',nbv
       !print *,'dist',d_rp/d_o,d_rm/d_o
    endif

    ! imax face
    ! ---------
    if (is_periodic(2)) then
       
       ! reference point in block
       x_r=bl(nbl)%xcorner(2)
       y_r=bl(nbl)%ycorner(2)
       z_r=bl(nbl)%zcorner(2)
       ! corresponding point in neighboring block
       nbv=bl(nbl)%BC(2)      
       x_n=bl(nbv)%xcorner(1)
       y_n=bl(nbv)%ycorner(1)
       z_n=bl(nbv)%zcorner(1)

       ! original distance
       d_o=sqrt((x_r-x_n)**2+(y_r-y_n)**2+(z_r-z_n)**2)
       
       ! apply rotation with plus sign
       y_rp= costp*y_r-sintp(2)*z_r
       z_rp= sintp(2)*y_r+costp*z_r
       d_rp=sqrt((x_r-x_n)**2+(y_rp-y_n)**2+(z_rp-z_n)**2)

       ! apply rotation with minus sign
       y_rm= costp*y_r+sintp(2)*z_r
       z_rm=-sintp(2)*y_r+costp*z_r
       d_rm=sqrt((x_r-x_n)**2+(y_rm-y_n)**2+(z_rm-z_n)**2)
       
       ! change sign of sin(theta_period) if negative rotation
       !if ((d_rp/d_o>1).and.(d_rm/d_o<1)) sintp(2)=-sintp(2)
       if ((d_rp/d_o<1).and.(d_rm/d_o>1)) then
          sintp(2)=-sintp(2)
          print *,'imax: change sign',nbl,iproc
       endif
       
       !print *,'imax: block',nbl,'neighbor',nbv
       !print *,'dist',d_rp/d_o,d_rm/d_o
    endif

    ! jmin face
    ! ---------
    if (is_periodic(3)) then
       
       ! reference point in block
       x_r=bl(nbl)%xcorner(1)
       y_r=bl(nbl)%ycorner(1)
       z_r=bl(nbl)%zcorner(1)
       ! corresponding point in neighboring block
       nbv=bl(nbl)%BC(3)      
       x_n=bl(nbv)%xcorner(3)
       y_n=bl(nbv)%ycorner(3)
       z_n=bl(nbv)%zcorner(3)

       ! original distance
       d_o=sqrt((x_r-x_n)**2+(y_r-y_n)**2+(z_r-z_n)**2)
       
       ! apply rotation with plus sign
       y_rp= costp*y_r-sintp(3)*z_r
       z_rp= sintp(3)*y_r+costp*z_r
       d_rp=sqrt((x_r-x_n)**2+(y_rp-y_n)**2+(z_rp-z_n)**2)

       ! apply rotation with minus sign
       y_rm= costp*y_r+sintp(3)*z_r
       z_rm=-sintp(3)*y_r+costp*z_r
       d_rm=sqrt((x_r-x_n)**2+(y_rm-y_n)**2+(z_rm-z_n)**2)
       
       ! change sign of sin(theta_period) if negative rotation
       !if ((d_rp/d_o>1).and.(d_rm/d_o<1)) sintp(3)=-sintp(3)
       if ((d_rp/d_o<1).and.(d_rm/d_o>1)) then
          sintp(3)=-sintp(3)
          print *,'jmin: change sign',nbl,iproc
       endif
       
       !print *,'jmin: block',nbl,'neighbor',nbv
       !print *,'dist',d_rp/d_o,d_rm/d_o
    endif

    ! jmax face
    ! ---------
    if (is_periodic(4)) then
       
       ! reference point in block
       x_r=bl(nbl)%xcorner(3)
       y_r=bl(nbl)%ycorner(3)
       z_r=bl(nbl)%zcorner(3)
       ! corresponding point in neighboring block
       nbv=bl(nbl)%BC(4)      
       x_n=bl(nbv)%xcorner(1)
       y_n=bl(nbv)%ycorner(1)
       z_n=bl(nbv)%zcorner(1)

       ! original distance
       d_o=sqrt((x_r-x_n)**2+(y_r-y_n)**2+(z_r-z_n)**2)
       
       ! apply rotation with plus sign
       y_rp= costp*y_r-sintp(4)*z_r
       z_rp= sintp(4)*y_r+costp*z_r
       d_rp=sqrt((x_r-x_n)**2+(y_rp-y_n)**2+(z_rp-z_n)**2)

       ! apply rotation with minus sign
       y_rm= costp*y_r+sintp(4)*z_r
       z_rm=-sintp(4)*y_r+costp*z_r
       d_rm=sqrt((x_r-x_n)**2+(y_rm-y_n)**2+(z_rm-z_n)**2)
       
       ! change sign of sin(theta_period) if negative rotation
       !if ((d_rp/d_o>1).and.(d_rm/d_o<1)) sintp(4)=-sintp(4)
       if ((d_rp/d_o<1).and.(d_rm/d_o>1)) then
          sintp(4)=-sintp(4)
          print *,'jmax: change sign',nbl,iproc
       endif
       
       !print *,'jmax: block',nbl,'neighbor',nbv
       !print *,'dist',d_rp/d_o,d_rm/d_o
    endif

    ! kmin face
    ! ---------
    if (is_periodic(5)) then
       
       ! reference point in block
       x_r=bl(nbl)%xcorner(1)
       y_r=bl(nbl)%ycorner(1)
       z_r=bl(nbl)%zcorner(1)
       ! corresponding point in neighboring block
       nbv=bl(nbl)%BC(5)      
       x_n=bl(nbv)%xcorner(5)
       y_n=bl(nbv)%ycorner(5)
       z_n=bl(nbv)%zcorner(5)

       ! original distance
       d_o=sqrt((x_r-x_n)**2+(y_r-y_n)**2+(z_r-z_n)**2)
       
       ! apply rotation with plus sign
       y_rp= costp*y_r-sintp(5)*z_r
       z_rp= sintp(5)*y_r+costp*z_r
       d_rp=sqrt((x_r-x_n)**2+(y_rp-y_n)**2+(z_rp-z_n)**2)

       ! apply rotation with minus sign
       y_rm= costp*y_r+sintp(5)*z_r
       z_rm=-sintp(5)*y_r+costp*z_r
       d_rm=sqrt((x_r-x_n)**2+(y_rm-y_n)**2+(z_rm-z_n)**2)
       
       ! change sign of sin(theta_period) if negative rotation
       !if ((d_rp/d_o>1).and.(d_rm/d_o<1)) sintp(5)=-sintp(5)
       if ((d_rp/d_o<1).and.(d_rm/d_o>1)) then
          sintp(5)=-sintp(5)
          print *,'kmin: change sign',nbl,iproc
       endif
       
       !print *,'kmin: block',nbl,'neighbor',nbv
       !print *,'dist',d_rp/d_o,d_rm/d_o
    endif

    ! kmax face
    ! ---------
    if (is_periodic(6)) then
       
       ! reference point in block
       x_r=bl(nbl)%xcorner(5)
       y_r=bl(nbl)%ycorner(5)
       z_r=bl(nbl)%zcorner(5)
       ! corresponding point in neighboring block
       nbv=bl(nbl)%BC(6)      
       x_n=bl(nbv)%xcorner(1)
       y_n=bl(nbv)%ycorner(1)
       z_n=bl(nbv)%zcorner(1)

       ! original distance
       d_o=sqrt((x_r-x_n)**2+(y_r-y_n)**2+(z_r-z_n)**2)
       
       ! apply rotation with plus sign
       y_rp= costp*y_r-sintp(6)*z_r
       z_rp= sintp(6)*y_r+costp*z_r
       d_rp=sqrt((x_r-x_n)**2+(y_rp-y_n)**2+(z_rp-z_n)**2)

       ! apply rotation with minus sign
       y_rm= costp*y_r+sintp(6)*z_r
       z_rm=-sintp(6)*y_r+costp*z_r
       d_rm=sqrt((x_r-x_n)**2+(y_rm-y_n)**2+(z_rm-z_n)**2)
       
       ! change sign of sin(theta_period) if negative rotation
       !if ((d_rp/d_o>1).and.(d_rm/d_o<1)) sintp(6)=-sintp(6)
       if ((d_rp/d_o<1).and.(d_rm/d_o>1)) then
          sintp(6)=-sintp(6)
          print *,'kmax: change sign',nbl,iproc
       endif
       
       !print *,'kmax: block',nbl,'neighbor',nbv
       !print *,'dist',d_rp/d_o,d_rm/d_o
    endif

    !call mpistop('test rot',0)

    
    ! determine periodicity direction
    ! ===============================
!!$    if (period_dir==0) &
!!$         call mpistop('Pb in defining angular periodicity direction. Please check.',0)

    if (iproc==0) print *,repeat('=',70)
    if (iproc==0) print 100,theta_period
    do i=1,6
       select case (period_dir(i))
       case(1)
          print 101,iproc,nbl
       case(2)
          print 102,iproc,nbl
       case(3)
          print 103,iproc,nbl
       end select
    enddo
    if (iproc==0) print *,repeat('=',70)
    
100 format(1x,'~> apply angular periodicity with angle ',f6.2,'°')
101 format(1x,'   in direction i for proc',i4,' of block',i4)
102 format(1x,'   in direction j for proc',i4,' of block',i4)
103 format(1x,'   in direction k for proc',i4,' of block',i4)
 
  end subroutine init_periodicity

  !==============================================================================
  subroutine share_block_corners
  !==============================================================================
    !> Assign and share block corners
    !> [ attribute of blocks ]
  !==============================================================================
    use mod_grid
    use mod_block
    use mod_mpi
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: n,nbl
    ! ---------------------------------------------------------------------------
    
    ! Current block
    ! =============
    nbl=nob(iproc)
    
    ! Determine block corners
    ! =======================
    !      6---------------8 
    !     /|               /
    !    / |              /|
    !   /  |             / |
    !  3----------------4  |
    !  |   | j ^  k     |  |
    !  |   |   | /      |  |
    !  |   |   |/       |  |
    !  |   |    ---> i  |  |
    !  |   5----------- |--7
    !  |  /             | /
    !  | /              |/
    !  1----------------2
    
    bl(nbl)%xcorner(1)=xgc3(1,1,1)
    bl(nbl)%xcorner(2)=xgc3(ngx,1,1)
    bl(nbl)%xcorner(3)=xgc3(1,ngy,1)
    bl(nbl)%xcorner(4)=xgc3(ngx,ngy,1)    
    bl(nbl)%xcorner(5)=xgc3(1,1,ngz)
    bl(nbl)%xcorner(6)=xgc3(1,ngy,ngz)
    bl(nbl)%xcorner(7)=xgc3(ngx,1,ngz)
    bl(nbl)%xcorner(8)=xgc3(ngx,ngy,ngz)

    bl(nbl)%ycorner(1)=ygc3(1,1,1)
    bl(nbl)%ycorner(2)=ygc3(ngx,1,1)
    bl(nbl)%ycorner(3)=ygc3(1,ngy,1)
    bl(nbl)%ycorner(4)=ygc3(ngx,ngy,1)
    bl(nbl)%ycorner(5)=ygc3(1,1,ngz)
    bl(nbl)%ycorner(6)=ygc3(1,ngy,ngz)
    bl(nbl)%ycorner(7)=ygc3(ngx,1,ngz)
    bl(nbl)%ycorner(8)=ygc3(ngx,ngy,ngz)

    bl(nbl)%zcorner(1)=zgc3(1,1,1)
    bl(nbl)%zcorner(2)=zgc3(ngx,1,1)
    bl(nbl)%zcorner(3)=zgc3(1,ngy,1)
    bl(nbl)%zcorner(4)=zgc3(ngx,ngy,1)
    bl(nbl)%zcorner(5)=zgc3(1,1,ngz)
    bl(nbl)%zcorner(6)=zgc3(1,ngy,ngz)
    bl(nbl)%zcorner(7)=zgc3(ngx,1,ngz)
    bl(nbl)%zcorner(8)=zgc3(ngx,ngy,ngz)

    call MPI_BARRIER(COMM_global,info)

    ! share block corner coordinate between procs
    do n=1,nbloc
       call MPI_BCAST(bl(n)%xcorner,8,MPI_DOUBLE_PRECISION,iproc_leader(n),COMM_global,info)
       call MPI_BCAST(bl(n)%ycorner,8,MPI_DOUBLE_PRECISION,iproc_leader(n),COMM_global,info)
       call MPI_BCAST(bl(n)%zcorner,8,MPI_DOUBLE_PRECISION,iproc_leader(n),COMM_global,info)
    enddo

  end subroutine share_block_corners
    
  !==============================================================================
  subroutine bc_angular_period
  !==============================================================================
    !> Main routine to apply "angular periodicity" Boundary Conditions
    !> --> correction of rhov & rhow after comm.
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    ! ---------------------------------------------------------------------------

    if (is_periodic(1)) call bc_angular_period_ox_imin

    if (is_periodic(2)) call bc_angular_period_ox_imax

    if (is_periodic(3)) call bc_angular_period_ox_jmin

    if (is_periodic(4)) call bc_angular_period_ox_jmax

    if (is_periodic(5)) call bc_angular_period_ox_kmin

    if (is_periodic(6)) call bc_angular_period_ox_kmax

  end subroutine bc_angular_period
  
  !==============================================================================
  subroutine bc_angular_period_n
  !==============================================================================
    !> Main routine to apply "angular periodicity" Boundary Conditions
    !> --> correction of rhov_n & rhow_n after comm.
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    ! ---------------------------------------------------------------------------

    if (is_periodic(1)) call bc_angular_period_ox_n_imin

    if (is_periodic(2)) call bc_angular_period_ox_n_imax

    if (is_periodic(3)) call bc_angular_period_ox_n_jmin

    if (is_periodic(4)) call bc_angular_period_ox_n_jmax

    if (is_periodic(5)) call bc_angular_period_ox_n_kmin

    if (is_periodic(6)) call bc_angular_period_ox_n_kmax

  end subroutine bc_angular_period_n
  
  !==============================================================================
  subroutine bc_angular_period_grad
  !==============================================================================
    !> Main routine to apply "angular periodicity" Boundary Conditions
    !> --> correction of velocity gradients after comm.
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    ! ---------------------------------------------------------------------------

    if (is_periodic(1)) call bc_angular_period_ox_g_imin

    if (is_periodic(2)) call bc_angular_period_ox_g_imax

    if (is_periodic(3)) call bc_angular_period_ox_g_jmin

    if (is_periodic(4)) call bc_angular_period_ox_g_jmax

    if (is_periodic(5)) call bc_angular_period_ox_g_kmin

    if (is_periodic(6)) call bc_angular_period_ox_g_kmax

  end subroutine bc_angular_period_grad
  
  !==============================================================================
  subroutine bc_angular_period_ox_imin
  !==============================================================================
    !> Apply "angular periodicity" Boundary Conditions at boundary imin
    !> based on a rotation of the coordinate system for vectorial quantities
    !> [rotational axis ox, angle theta_period (positive in clockwise direction)]
    !> --> correction of rhov & rhow after comm.
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: j,k
    real(wp), dimension(1-ngh:0) :: rov_m,row_m
    ! ---------------------------------------------------------------------------
    
    do k=1,nz
       do j=1,ny
          rov_m=rhov(1-ngh:0,j,k)
          row_m=rhow(1-ngh:0,j,k)
          rhov(1-ngh:0,j,k)= costp*rov_m-sintp(1)*row_m
          rhow(1-ngh:0,j,k)= sintp(1)*rov_m+costp*row_m
       enddo
    enddo

  end subroutine bc_angular_period_ox_imin

  !==============================================================================
  subroutine bc_angular_period_ox_n_imin
  !==============================================================================
    !> Apply "angular periodicity" Boundary Conditions at boundary imin
    !> based on a rotation of the coordinate system for vectorial quantities
    !> [rotational axis ox, angle theta_period (positive in clockwise direction)]
    !> --> correction of rhov_n & rhow_n after comm.
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: j,k
    real(wp), dimension(1-ngh:0) :: rov_m,row_m
    ! ---------------------------------------------------------------------------

    do k=1,nz
       do j=1,ny
          rov_m=rhov_n(1-ngh:0,j,k)
          row_m=rhow_n(1-ngh:0,j,k)
          rhov_n(1-ngh:0,j,k)= costp*rov_m-sintp(1)*row_m
          rhow_n(1-ngh:0,j,k)= sintp(1)*rov_m+costp*row_m
       enddo
    enddo

  end subroutine bc_angular_period_ox_n_imin

  !==============================================================================
  subroutine bc_angular_period_ox_g_imin
  !==============================================================================
    !> Apply "angular periodicity" Boundary Conditions at boundary imin
    !> based on a rotation of the coordinate system for vectorial quantities
    !> [rotational axis ox, angle theta_period (positive in clockwise direction)]
    !> --> correction of velocity gradients after comm.
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: j,k
    real(wp), dimension(1-ngh_v:0) :: dv_m,dw_m
    ! ---------------------------------------------------------------------------

    ! x-gradients for v,w
    do k=1,nz
       do j=1,ny
          dv_m=dvx(1-ngh_v:0,j,k)
          dw_m=dwx(1-ngh_v:0,j,k)
          dvx(1-ngh_v:0,j,k)= costp*dv_m-sintp(1)*dw_m
          dwx(1-ngh_v:0,j,k)= sintp(1)*dv_m+costp*dw_m
       enddo
    enddo
    ! y-gradients for v,w
    do k=1,nz
       do j=1,ny
          dv_m=dvy(1-ngh_v:0,j,k)
          dw_m=dwy(1-ngh_v:0,j,k)
          dvy(1-ngh_v:0,j,k)= costp*dv_m-sintp(1)*dw_m
          dwy(1-ngh_v:0,j,k)= sintp(1)*dv_m+costp*dw_m
       enddo
    enddo
    ! z-gradients for v,w
    do k=1,nz
       do j=1,ny
          dv_m=dvz(1-ngh_v:0,j,k)
          dw_m=dwz(1-ngh_v:0,j,k)
          dvz(1-ngh_v:0,j,k)= costp*dv_m-sintp(1)*dw_m
          dwz(1-ngh_v:0,j,k)= sintp(1)*dv_m+costp*dw_m
       enddo
    enddo

  end subroutine bc_angular_period_ox_g_imin

  !==============================================================================
  subroutine bc_angular_period_ox_imax
  !==============================================================================
    !> Apply "angular periodicity" Boundary Conditions at boundary imax
    !> based on a rotation of the coordinate system for vectorial quantities
    !> [rotational axis ox, angle theta_period (positive in clockwise direction)]
    !> --> correction of rhov & rhow after comm.
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: j,k
    real(wp), dimension(nx+1:nx+ngh) :: rov_p,row_p
    ! ---------------------------------------------------------------------------

    do k=1,nz
       do j=1,ny
          rov_p=rhov(nx+1:nx+ngh,j,k)
          row_p=rhow(nx+1:nx+ngh,j,k)
          rhov(nx+1:nx+ngh,j,k)= costp*rov_p-sintp(2)*row_p
          rhow(nx+1:nx+ngh,j,k)= sintp(2)*rov_p+costp*row_p
       enddo
    enddo
   
  end subroutine bc_angular_period_ox_imax

  !==============================================================================
  subroutine bc_angular_period_ox_n_imax
  !==============================================================================
    !> Apply "angular periodicity" Boundary Conditions at boundary imax
    !> based on a rotation of the coordinate system for vectorial quantities
    !> [rotational axis ox, angle theta_period (positive in clockwise direction)]
    !> --> correction of rhov_n & rhow_n after comm.
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: j,k
    real(wp), dimension(nx+1:nx+ngh) :: rov_p,row_p
    ! ---------------------------------------------------------------------------

    do k=1,nz
       do j=1,ny
          rov_p=rhov_n(nx+1:nx+ngh,j,k)
          row_p=rhow_n(nx+1:nx+ngh,j,k)
          rhov_n(nx+1:nx+ngh,j,k)= costp*rov_p-sintp(2)*row_p
          rhow_n(nx+1:nx+ngh,j,k)= sintp(2)*rov_p+costp*row_p
       enddo
    enddo
   
  end subroutine bc_angular_period_ox_n_imax

  !==============================================================================
  subroutine bc_angular_period_ox_g_imax
  !==============================================================================
    !> Apply "angular periodicity" Boundary Conditions at boundary imax
    !> based on a rotation of the coordinate system for vectorial quantities
    !> [rotational axis ox, angle theta_period (positive in clockwise direction)]
    !> --> correction of velocity gradients after comm.
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: j,k
    real(wp), dimension(nx+1:nx+ngh_v) :: dv_p,dw_p
    ! ---------------------------------------------------------------------------

    ! x-gradients for v,w
    do k=1,nz
       do j=1,ny
          dv_p=dvx(nx+1:nx+ngh_v,j,k)
          dw_p=dwx(nx+1:nx+ngh_v,j,k)
          dvx(nx+1:nx+ngh_v,j,k)= costp*dv_p-sintp(2)*dw_p
          dwx(nx+1:nx+ngh_v,j,k)= sintp(2)*dv_p+costp*dw_p
       enddo
    enddo
    ! y-gradients for v,w
    do k=1,nz
       do j=1,ny
          dv_p=dvy(nx+1:nx+ngh_v,j,k)
          dw_p=dwy(nx+1:nx+ngh_v,j,k)
          dvy(nx+1:nx+ngh_v,j,k)= costp*dv_p-sintp(2)*dw_p
          dwy(nx+1:nx+ngh_v,j,k)= sintp(2)*dv_p+costp*dw_p
       enddo
    enddo
    ! z-gradients for v,w
    do k=1,nz
       do j=1,ny
          dv_p=dvz(nx+1:nx+ngh_v,j,k)
          dw_p=dwz(nx+1:nx+ngh_v,j,k)
          dvz(nx+1:nx+ngh_v,j,k)= costp*dv_p-sintp(2)*dw_p
          dwz(nx+1:nx+ngh_v,j,k)= sintp(2)*dv_p+costp*dw_p
       enddo
    enddo
   
  end subroutine bc_angular_period_ox_g_imax

  !==============================================================================
  subroutine bc_angular_period_ox_jmin
  !==============================================================================
    !> Apply "angular periodicity" Boundary Conditions at boundary jmin
    !> based on a rotation of the coordinate system for vectorial quantities
    !> [rotational axis ox, angle theta_period (positive in clockwise direction)]
    !> --> correction of rhov & rhow after comm.
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,k
    real(wp), dimension(1-ngh:0) :: rov_m,row_m
    ! ---------------------------------------------------------------------------
    
    do k=1,nz
       do i=1,nx
          rov_m=rhov(i,1-ngh:0,k)
          row_m=rhow(i,1-ngh:0,k)
          rhov(i,1-ngh:0,k)= costp*rov_m-sintp(3)*row_m
          rhow(i,1-ngh:0,k)= sintp(3)*rov_m+costp*row_m
       enddo
    enddo

  end subroutine bc_angular_period_ox_jmin

  !==============================================================================
  subroutine bc_angular_period_ox_n_jmin
  !==============================================================================
    !> Apply "angular periodicity" Boundary Conditions at boundary jmin
    !> based on a rotation of the coordinate system for vectorial quantities
    !> [rotational axis ox, angle theta_period (positive in clockwise direction)]
    !> --> correction of rhov_n & rhow_n after comm.
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,k
    real(wp), dimension(1-ngh:0) :: rov_m,row_m
    ! ---------------------------------------------------------------------------

    do k=1,nz
       do i=1,nx
          rov_m=rhov_n(i,1-ngh:0,k)
          row_m=rhow_n(i,1-ngh:0,k)
          rhov_n(i,1-ngh:0,k)= costp*rov_m-sintp(3)*row_m
          rhow_n(i,1-ngh:0,k)= sintp(3)*rov_m+costp*row_m
       enddo
    enddo

  end subroutine bc_angular_period_ox_n_jmin

  !==============================================================================
  subroutine bc_angular_period_ox_g_jmin
  !==============================================================================
    !> Apply "angular periodicity" Boundary Conditions at boundary jmin
    !> based on a rotation of the coordinate system for vectorial quantities
    !> [rotational axis ox, angle theta_period (positive in clockwise direction)]
    !> --> correction of velocity gradients after comm.
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,k
    real(wp), dimension(1-ngh_v:0) :: dv_m,dw_m
    ! ---------------------------------------------------------------------------

    ! x-gradients for v,w
    do k=1,nz
       do i=1,nx
          dv_m=dvx(i,1-ngh_v:0,k)
          dw_m=dwx(i,1-ngh_v:0,k)
          dvx(i,1-ngh_v:0,k)= costp*dv_m-sintp(3)*dw_m
          dwx(i,1-ngh_v:0,k)= sintp(3)*dv_m+costp*dw_m
       enddo
    enddo
    ! y-gradients for v,w
    do k=1,nz
       do i=1,nx
          dv_m=dvy(i,1-ngh_v:0,k)
          dw_m=dwy(i,1-ngh_v:0,k)
          dvy(i,1-ngh_v:0,k)= costp*dv_m-sintp(3)*dw_m
          dwy(i,1-ngh_v:0,k)= sintp(3)*dv_m+costp*dw_m
       enddo
    enddo
    ! z-gradients for v,w
    do k=1,nz
       do i=1,nx
          dv_m=dvz(i,1-ngh_v:0,k)
          dw_m=dwz(i,1-ngh_v:0,k)
          dvz(i,1-ngh_v:0,k)= costp*dv_m-sintp(3)*dw_m
          dwz(i,1-ngh_v:0,k)= sintp(3)*dv_m+costp*dw_m
       enddo
    enddo

  end subroutine bc_angular_period_ox_g_jmin

  !==============================================================================
  subroutine bc_angular_period_ox_jmax
  !==============================================================================
    !> Apply "angular periodicity" Boundary Conditions at boundary jmax
    !> based on a rotation of the coordinate system for vectorial quantities
    !> [rotational axis ox, angle theta_period (positive in clockwise direction)]
    !> --> correction of rhov & rhow after comm.
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,k
    real(wp), dimension(ny+1:ny+ngh) :: rov_p,row_p
    ! ---------------------------------------------------------------------------

    do k=1,nz
       do i=1,nx
          rov_p=rhov(i,ny+1:ny+ngh,k)
          row_p=rhow(i,ny+1:ny+ngh,k)
          rhov(i,ny+1:ny+ngh,k)= costp*rov_p-sintp(4)*row_p
          rhow(i,ny+1:ny+ngh,k)= sintp(4)*rov_p+costp*row_p
       enddo
    enddo
   
  end subroutine bc_angular_period_ox_jmax

  !==============================================================================
  subroutine bc_angular_period_ox_n_jmax
  !==============================================================================
    !> Apply "angular periodicity" Boundary Conditions at boundary jmax
    !> based on a rotation of the coordinate system for vectorial quantities
    !> [rotational axis ox, angle theta_period (positive in clockwise direction)]
    !> --> correction of rhov_n & rhow_n after comm.
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,k
    real(wp), dimension(ny+1:ny+ngh) :: rov_p,row_p
    ! ---------------------------------------------------------------------------

    do k=1,nz
       do i=1,nx
          rov_p=rhov_n(i,ny+1:ny+ngh,k)
          row_p=rhow_n(i,ny+1:ny+ngh,k)
          rhov_n(i,ny+1:ny+ngh,k)= costp*rov_p-sintp(4)*row_p
          rhow_n(i,ny+1:ny+ngh,k)= sintp(4)*rov_p+costp*row_p
       enddo
    enddo
   
  end subroutine bc_angular_period_ox_n_jmax

  !==============================================================================
  subroutine bc_angular_period_ox_g_jmax
  !==============================================================================
    !> Apply "angular periodicity" Boundary Conditions at boundary jmax
    !> based on a rotation of the coordinate system for vectorial quantities
    !> [rotational axis ox, angle theta_period (positive in clockwise direction)]
    !> --> correction of velocity gradients after comm.
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,k
    real(wp), dimension(ny+1:ny+ngh_v) :: dv_p,dw_p
    ! ---------------------------------------------------------------------------

    ! x-gradients for v,w
    do k=1,nz
       do i=1,nx
          dv_p=dvx(i,ny+1:ny+ngh_v,k)
          dw_p=dwx(i,ny+1:ny+ngh_v,k)
          dvx(i,ny+1:ny+ngh_v,k)= costp*dv_p-sintp(4)*dw_p
          dwx(i,ny+1:ny+ngh_v,k)= sintp(4)*dv_p+costp*dw_p
       enddo
    enddo
    ! y-gradients for v,w
    do k=1,nz
       do i=1,nx
          dv_p=dvy(i,ny+1:ny+ngh_v,k)
          dw_p=dwy(i,ny+1:ny+ngh_v,k)
          dvy(i,ny+1:ny+ngh_v,k)= costp*dv_p-sintp(4)*dw_p
          dwy(i,ny+1:ny+ngh_v,k)= sintp(4)*dv_p+costp*dw_p
       enddo
    enddo
    ! z-gradients for v,w
    do k=1,nz
       do i=1,nx
          dv_p=dvz(i,ny+1:ny+ngh_v,k)
          dw_p=dwz(i,ny+1:ny+ngh_v,k)
          dvz(i,ny+1:ny+ngh_v,k)= costp*dv_p-sintp(4)*dw_p
          dwz(i,ny+1:ny+ngh_v,k)= sintp(4)*dv_p+costp*dw_p
       enddo
    enddo
   
  end subroutine bc_angular_period_ox_g_jmax

  !==============================================================================
  subroutine bc_angular_period_ox_kmin
  !==============================================================================
    !> Apply "angular periodicity" Boundary Conditions at boundary kmin
    !> based on a rotation of the coordinate system for vectorial quantities
    !> [rotational axis ox, angle theta_period (positive in clockwise direction)]
    !> --> correction of rhov & rhow after comm.
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j
    real(wp), dimension(1-ngh:0) :: rov_m,row_m
    ! ---------------------------------------------------------------------------

    do j=1,ny
       do i=1,nx
          rov_m=rhov(i,j,1-ngh:0)
          row_m=rhow(i,j,1-ngh:0)
          rhov(i,j,1-ngh:0)= costp*rov_m-sintp(5)*row_m
          rhow(i,j,1-ngh:0)= sintp(5)*rov_m+costp*row_m

       enddo
    enddo

  end subroutine bc_angular_period_ox_kmin

  !==============================================================================
  subroutine bc_angular_period_ox_n_kmin
  !==============================================================================
    !> Apply "angular periodicity" Boundary Conditions at boundary kmin
    !> based on a rotation of the coordinate system for vectorial quantities
    !> [rotational axis ox, angle theta_period (positive in clockwise direction)]
    !> --> correction of rhov_n & rhow_n after comm.
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j
    real(wp), dimension(1-ngh:0) :: rov_m,row_m
    ! ---------------------------------------------------------------------------

    do j=1,ny
       do i=1,nx
          rov_m=rhov_n(i,j,1-ngh:0)
          row_m=rhow_n(i,j,1-ngh:0)
          rhov_n(i,j,1-ngh:0)= costp*rov_m-sintp(5)*row_m
          rhow_n(i,j,1-ngh:0)= sintp(5)*rov_m+costp*row_m
       enddo
    enddo

  end subroutine bc_angular_period_ox_n_kmin

  !==============================================================================
  subroutine bc_angular_period_ox_g_kmin
  !==============================================================================
    !> Apply "angular periodicity" Boundary Conditions at boundary kmin
    !> based on a rotation of the coordinate system for vectorial quantities
    !> [rotational axis ox, angle theta_period (positive in clockwise direction)]
    !> --> correction of velocity gradients after comm.
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j
    real(wp), dimension(1-ngh_v:0) :: dv_m,dw_m
    ! ---------------------------------------------------------------------------

    ! x-gradients for v,w
    do j=1,ny
       do i=1,nx
          dv_m=dvx(i,j,1-ngh_v:0)
          dw_m=dwx(i,j,1-ngh_v:0)
          dvx(i,j,1-ngh_v:0)= costp*dv_m-sintp(5)*dw_m
          dwx(i,j,1-ngh_v:0)= sintp(5)*dv_m+costp*dw_m
       enddo
    enddo
    ! y-gradients for v,w
    do j=1,ny
       do i=1,nx
          dv_m=dvy(i,j,1-ngh_v:0)
          dw_m=dwy(i,j,1-ngh_v:0)
          dvy(i,j,1-ngh_v:0)= costp*dv_m-sintp(5)*dw_m
          dwy(i,j,1-ngh_v:0)= sintp(5)*dv_m+costp*dw_m
       enddo
    enddo
    ! z-gradients for v,w
    do j=1,ny
       do i=1,nx
          dv_m=dvz(i,j,1-ngh_v:0)
          dw_m=dwz(i,j,1-ngh_v:0)
          dvz(i,j,1-ngh_v:0)= costp*dv_m-sintp(5)*dw_m
          dwz(i,j,1-ngh_v:0)= sintp(5)*dv_m+costp*dw_m
       enddo
    enddo

  end subroutine bc_angular_period_ox_g_kmin

  !==============================================================================
  subroutine bc_angular_period_ox_kmax
  !==============================================================================
    !> Apply "angular periodicity" Boundary Conditions at boundary kmax
    !> based on a rotation of the coordinate system for vectorial quantities
    !> [rotational axis ox, angle theta_period (positive in clockwise direction)]
    !> --> correction of rhov & rhow after comm.
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j
    real(wp), dimension(nz+1:nz+ngh) :: rov_p,row_p
    ! ---------------------------------------------------------------------------
    
    do j=1,ny
       do i=1,nx
          rov_p=rhov(i,j,nz+1:nz+ngh)
          row_p=rhow(i,j,nz+1:nz+ngh)
          rhov(i,j,nz+1:nz+ngh)= costp*rov_p-sintp(6)*row_p
          rhow(i,j,nz+1:nz+ngh)= sintp(6)*rov_p+costp*row_p
       enddo
    enddo

  end subroutine bc_angular_period_ox_kmax

  !==============================================================================
  subroutine bc_angular_period_ox_n_kmax
  !==============================================================================
    !> Apply "angular periodicity" Boundary Conditions at boundary kmax
    !> based on a rotation of the coordinate system for vectorial quantities
    !> [rotational axis ox, angle theta_period (positive in clockwise direction)]
    !> --> correction of rhov_n & rhow_n after comm.
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j
    real(wp), dimension(nz+1:nz+ngh) :: rov_p,row_p
    ! ---------------------------------------------------------------------------

    do j=1,ny
       do i=1,nx
          rov_p=rhov_n(i,j,nz+1:nz+ngh)
          row_p=rhow_n(i,j,nz+1:nz+ngh)
          rhov_n(i,j,nz+1:nz+ngh)= costp*rov_p-sintp(6)*row_p
          rhow_n(i,j,nz+1:nz+ngh)= sintp(6)*rov_p+costp*row_p
       enddo
    enddo
   
  end subroutine bc_angular_period_ox_n_kmax

  !==============================================================================
  subroutine bc_angular_period_ox_g_kmax
  !==============================================================================
    !> Apply "angular periodicity" Boundary Conditions at boundary kmax
    !> based on a rotation of the coordinate system for vectorial quantities
    !> [rotational axis ox, angle theta_period (positive in clockwise direction)]
    !> --> correction of velocity gradients after comm.
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j
    real(wp), dimension(nz+1:nz+ngh_v) :: dv_p,dw_p
    ! ---------------------------------------------------------------------------

    ! x-gradients for v,w
    do j=1,ny
       do i=1,nx
          dv_p=dvx(i,j,nz+1:nz+ngh_v)
          dw_p=dwx(i,j,nz+1:nz+ngh_v)
          dvx(i,j,nz+1:nz+ngh_v)= costp*dv_p-sintp(6)*dw_p
          dwx(i,j,nz+1:nz+ngh_v)= sintp(6)*dv_p+costp*dw_p
       enddo
    enddo
    ! y-gradients for v,w
    do j=1,ny
       do i=1,nx
          dv_p=dvy(i,j,nz+1:nz+ngh_v)
          dw_p=dwy(i,j,nz+1:nz+ngh_v)
          dvy(i,j,nz+1:nz+ngh_v)= costp*dv_p-sintp(6)*dw_p
          dwy(i,j,nz+1:nz+ngh_v)= sintp(6)*dv_p+costp*dw_p
       enddo
    enddo
    ! z-gradients for v,w
    do j=1,ny
       do i=1,nx
          dv_p=dvz(i,j,nz+1:nz+ngh_v)
          dw_p=dwz(i,j,nz+1:nz+ngh_v)
          dvz(i,j,nz+1:nz+ngh_v)= costp*dv_p-sintp(6)*dw_p
          dwz(i,j,nz+1:nz+ngh_v)= sintp(6)*dv_p+costp*dw_p
       enddo
    enddo
   
  end subroutine bc_angular_period_ox_g_kmax

end module mod_bc_periodicity
