!==============================================================================================
module mod_io
!==============================================================================================
  !> author: Luca Sciacovelli
  !> date: April 2018
  !> Module for Input/Output routines
!==============================================================================================
  use mod_constant ! <- for is_mean0 , STBL in mod_io_planes for stats TO BE CHANGED
  use mod_mpi_part ! ???? mod_mpi sufficient ?? maybe for COMM_intrablock in planes
  use mod_flow     ! <- for grid & solution
  use mod_time     ! <- for tstar
  use mod_tecplot  ! <- module defining Tecplot/Bin I/O
  use mod_io_restartTD
  implicit none
  ! -------------------------------------------------------------------------------------------
  ! options for creating binary file
  integer, parameter :: WRITE=1,READ=2,WRITE_LEADER=3,READ_LEADER=4,READ_DEB=5
  ! -------------------------------------------------------------------------------------------
  ! Value of timestamp and type of operation to perform on binary file
  character(len=50) :: binfile='./restart0000_0000.plt'
  character(len=50) :: TDfile
  character(len=9) :: filestamp
  character(len=4) :: filext_read,filext_write
  logical :: is_IOtec_read,is_IOtec_write
  logical :: is_timestamp
  ! -------------------------------------------------------------------------------------------
  ! Variables for MPI writing
  integer :: ndata
  real(wp), dimension(:,:,:,:), allocatable, target :: dummy
  character(len=10), dimension(:), allocatable, target :: dataname
  ! -------------------------------------------------------------------------------------------
  ! Derived type with all informations for one block (dim is number of blocks)
  type(teciotype), dimension(:), allocatable  :: field
  ! -------------------------------------------------------------------------------------------
  ! Derived type with all informations for one block (dim is number of blocks)
  ! ~> this type is used in INTERP mode to store old field to be interpolated
  type(teciotype), dimension(:), allocatable  :: field_o
  ! -------------------------------------------------------------------------------------------
  ! Derived type with all informations for one block (dim is number of blocks)
  ! ~> this type is used to read/write extended grid
  type(teciotype), dimension(:), allocatable  :: grid_ex_
  ! -------------------------------------------------------------------------------------------
  integer :: COMM_rw_grid

contains

  !==============================================================================================
  subroutine read_write_volume(binfile,TDfile,operation,gridname)
  !==============================================================================================
    !> author: Luca Sciacovelli
    !> date: April 2018
    !> read/write binary restart volume
  !==============================================================================================
    use mod_utils
    use mod_flow0    ! <- for time-averaged solution
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Input/Output arguments
    character(len=*), intent(in) :: binfile,TDfile
    integer, intent(in) :: operation
    character(len=*), optional :: gridname !  TO BE CHANGED ??? -> never used at the present time
    ! -------------------------------------------------------------------------------------------
    ! Local variables
    integer  :: i,j,k,m,ndata,filetype
    logical  :: iexist
    real(wp) :: soltime
    character(len=30) :: gridfile
    ! -------------------------------------------------------------------------------------------
    ! for LITTLE_ENDIAN <-> BIG_ENDIAN conversion *** COMMENTED ***
    ! integer :: nbytes_dbl,ierror
    ! -------------------------------------------------------------------------------------------

    ! =========================================
    ! Write/Read GRID (only for tecplot format)
    ! =========================================

    ! Name of the grid file ! TO BE CHANGED ???
    ! =====================
    ! can be defined as optional input of the subroutine but never used at the present time
    if (present(gridname)) then
       gridfile=gridname
    else
       gridfile='./grid_bl'//trim(numchar(nob(iproc)))//'.plt'
    endif

    ! Prepare number and name of data: x,y,z coordinates (Att. CARTESIAN only)
    ! ===============================
    ndata=3
    allocate(dataname(ndata),varlist(ndata))
    dataname(1:ndata) = (/'X','Y','Z'/)

    ! Write volume grid : done only if tecplot format is selected
    ! ===================
    if ((operation==WRITE).and.(is_IOtec_write)) then

       ! check if grid file already exists
       ! ---------------------------------
       inquire(file=gridfile,exist=iexist)

       ! if not: write GRID
       ! ------------------
       if ((.not.iexist).or.(idepart==FROM_SCRATCH)) then

          ! allocate and initialize dummy variable for datas
          allocate(dummy(nx,ny,nz,ndata))
          dummy = 0.0_wp

          ! fill data x,y,z in dummy array
          if (is_curv3) then
             do k=1,nz
                do j=1,ny
                   do i=1,nx
                      dummy(i,j,k,1)=xc3(i,j,k)
                      dummy(i,j,k,2)=yc3(i,j,k)
                      dummy(i,j,k,3)=zc3(i,j,k)
                   enddo
                enddo
             enddo
          else
             if (is_curv) then
                do j=1,ny
                   do i=1,nx
                      dummy(i,j,:,1)=xc(i,j)
                      dummy(i,j,:,2)=yc(i,j)
                   enddo
                enddo
             else
                do i=1,nx
                   dummy(i,:,:,1)=x(i)
                enddo
                do j=1,ny
                   dummy(:,j,:,2)=y(j)
                enddo
             endif
             do k=1,nz
                dummy(:,:,k,3)=z(k)
             enddo
          endif

          ! fill varlist (pointer on datas)
          do m=1,ndata
             varlist(m)%data=>dummy(:,:,:,m)
             varlist(m)%name= dataname(m)
          enddo

          ! filetype [0:full; 1:only grid; 2:only solution]
          filetype=1
          ! non-dimensional time
          soltime=tstar

          ! write grid block in separate file
          if (iproc.eq.0) write(*,*) 'Writing gridfile ~>',trim(gridfile)
          call write_tec(gridfile,filetype,soltime,field(nob(iproc)))

          ! free temporary data
          deallocate(dummy)

       endif ! end if test iexist

    endif ! end write grid if tecplot format

    ! Read volume grid : done only if tecplot format is selected
    ! ==================
    if ((operation==READ).and.(is_IOtec_read)) then

       ! allocate and initialize dummy variable for datas
       allocate(dummy(nx,ny,nz,ndata))
       dummy = 0.0_wp

       ! fill varlist (pointer on datas)
       do m=1,ndata
          varlist(m)%data=>dummy(:,:,:,m)
          varlist(m)%name= dataname(m)
       enddo

       ! check file gridfile
       call mpicheckfile(gridfile)

       ! read grid block
       if (iproc.eq.0) write(*,*) 'Reading gridfile ~>', trim(gridfile)
       call read_tec(gridfile,field(nob(iproc)))

       ! fill data x,y,z from dummy array (Att. CARTESIAN only) TO BE CHANGED for curvilinear
       if (is_curv3) then
          xc3(1:nx,1:ny,1:nz)=dummy(1:nx,1:ny,1:nz,1)
          yc3(1:nx,1:ny,1:nz)=dummy(1:nx,1:ny,1:nz,2)
          zc3(1:nx,1:ny,1:nz)=dummy(1:nx,1:ny,1:nz,3)
       else
          if (is_curv) then
             xc(1:nx,1:ny)=dummy(1:nx,1:ny,1,1)
             yc(1:nx,1:ny)=dummy(1:nx,1:ny,1,2)
          else
             x(1:nx)=dummy(1:nx,1,1,1)
             y(1:ny)=dummy(1,1:ny,1,2)
          endif
          z(1:nz)=dummy(1,1,1:nz,3)
       endif

       ! free temporary data
       deallocate(dummy)

    endif ! end read grid if tecplot format

    ! free pointer varlist & dataname
    ! ===============================
    deallocate(varlist,dataname)

    ! ==================================================
    ! Read/Write instantaneous solution in restart files
    ! ==================================================

    ! Restart files store conservative variables in global volume
    ! ===========================================================
    if ((is_RANS).and.(ndeb_RANS.le.ntotal)) then
       ndata = 6
       allocate(varlist(ndata),dataname(ndata),dummy(nx,ny,nz,ndata))
       dummy=0.0_wp
       if (model_RANS.eq.'SA') then
          dataname(1:ndata)=(/'rro','rou','rov','row','roe','nutil'/)
       endif
    else
       ndata = 5
       !ndata = 7
       allocate(varlist(ndata),dataname(ndata),dummy(nx,ny,nz,ndata))
       dummy=0.0_wp
       dataname(1:ndata)=(/'rro','rou','rov','row','roe'/)
       !dataname(1:ndata)=(/'rro','rou','rov','row','roe','cfi','cfj'/)
    endif
    
    ! fill varlist (pointer on datas)
    ! ===============================
    do m=1,ndata
       varlist(m)%data=>dummy(:,:,:,m)
       varlist(m)%name= dataname(m)
    enddo

    ! Write volume instantaneous solution
    ! ===================================
    if (operation==WRITE) then

       ! fill data: instantaneous conservative variables
       dummy(:,:,:,1)= rho(1:nx,1:ny,1:nz)
       dummy(:,:,:,2)=rhou(1:nx,1:ny,1:nz)
       dummy(:,:,:,3)=rhov(1:nx,1:ny,1:nz)
       dummy(:,:,:,4)=rhow(1:nx,1:ny,1:nz)
       dummy(:,:,:,5)=rhoe(1:nx,1:ny,1:nz)
       !dummy(:,:,:,6)=cfl_i
       !dummy(:,:,:,7)=cfl_j
       if ((is_RANS).and.(ndeb_RANS.le.ntotal)) then
          if (model_RANS.eq.'SA') dummy(:,:,:,6)=nutil(1:nx,1:ny,1:nz)
       endif

       ! filetype [0:full; 1:only grid; 2:only solution]
       filetype=2
       ! non-dimensional time
       soltime=tstar

       ! write instantaneous solution block
       if (iproc.eq.0) write(*,*) 'Writing restart ~>',trim(binfile)
       call write_tec(binfile,filetype,soltime,field(nob(iproc)))

    endif ! end write instantaneous solution

    ! Read volume instantaneous solution
    ! ==================================
    if (operation==READ) then

       ! check file existence
       call mpicheckfile(binfile)

       ! read instantaneous solution block
       if (iproc.eq.0) write(*,*) 'Reading restart ~>',trim(binfile)
       call read_tec(binfile,field(nob(iproc)))

       ! If on Big-Endian machine (Turing), swap bytes to write in Little-Endian
       ! =======================================================================
       !if (ichar(transfer(1,'a'))==0) then
       !   ! Taille du type de base MPI_DOUBLE_PRECISION
       !   call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,nbytes_dbl,ierror)
       !   do m=1,ndata
       !      call conv(dummy(:,:,:,m),nbytes_dbl,size(dummy(:,:,:,m)))
       !   enddo
       !endif

       ! fill data: instantaneous conservative variables
        rho(1:nx,1:ny,1:nz)=dummy(:,:,:,1)
       rhou(1:nx,1:ny,1:nz)=dummy(:,:,:,2)
       rhov(1:nx,1:ny,1:nz)=dummy(:,:,:,3)
       rhow(1:nx,1:ny,1:nz)=dummy(:,:,:,4)
       rhoe(1:nx,1:ny,1:nz)=dummy(:,:,:,5)
       !cfl_i=dummy(:,:,:,6)
       !cfl_j=dummy(:,:,:,7)
       if ((is_RANS).and.(ndeb_RANS.le.ntotal)) then
          if (model_RANS.eq.'SA') nutil(1:nx,1:ny,1:nz)=dummy(:,:,:,6)
       endif

    endif ! end read instantaneous solution

    ! free pointer varlist, dataname, dummy
    ! =====================================
    deallocate(varlist,dataname,dummy)

    ! ===================================================
    ! Read/Write online mean fields (primitive variables)
    ! ===================================================
    ! ~~> useful when Tam & Dong's BC are applied
    if (is_mean0) then

       ! Mean0 file store primitive variables (T is useful for non-ideal gases)
       ! ====================================
       ndata = 6
       allocate(varlist(ndata),dataname(ndata),dummy(nx,ny,nz,ndata))
       dataname(1:ndata) = (/'u','v','w','p','r','T'/)

       ! fill varlist (pointer on datas)
       ! ===============================
       do m=1,ndata
          varlist(m)%data=>dummy(:,:,:,m)
          varlist(m)%name= dataname(m)
       enddo

       ! Write volume mean solution
       ! ==========================
       if (operation==WRITE) then

          ! fill data: mean primitive variables
          dummy(:,:,:,1)=  u0(1:nx,1:ny,1:nz)
          dummy(:,:,:,2)=  v0(1:nx,1:ny,1:nz)
          dummy(:,:,:,3)=  w0(1:nx,1:ny,1:nz)
          dummy(:,:,:,4)=  p0(1:nx,1:ny,1:nz)
          dummy(:,:,:,5)=rho0(1:nx,1:ny,1:nz)
          dummy(:,:,:,6)=  T0(1:nx,1:ny,1:nz)

          ! filetype [0:full; 1:only grid; 2:only solution]
          filetype=2
          ! non-dimensional time
          soltime=tstar

          ! write mean solution block
          if (iproc.eq.0) write(*,*) 'Writing averaged volumes'
          call write_tec('mean0_bl'//trim(numchar(nob(iproc)))//filext_write,filetype,soltime,field(nob(iproc)))

       endif ! end write mean solution

       ! Read volume mean solution
       ! =========================
       if (operation==READ) then

          ! check file existence
          call mpicheckfile('mean0_bl'//trim(numchar(nob(iproc)))//filext_read)

          ! read mean solution block
          if (iproc.eq.0) write(*,*) 'Reading averaged volumes'
          call read_tec('mean0_bl'//trim(numchar(nob(iproc)))//filext_read,field(nob(iproc)))

          ! fill data: mean primitive variables
            u0(1:nx,1:ny,1:nz)=dummy(:,:,:,1)
            v0(1:nx,1:ny,1:nz)=dummy(:,:,:,2)
            w0(1:nx,1:ny,1:nz)=dummy(:,:,:,3)
            p0(1:nx,1:ny,1:nz)=dummy(:,:,:,4)
          rho0(1:nx,1:ny,1:nz)=dummy(:,:,:,5)
            T0(1:nx,1:ny,1:nz)=dummy(:,:,:,6)

       endif ! end read mean solution

       ! free pointer varlist, dataname, dummy
       ! =====================================
       deallocate(varlist, dataname, dummy)
       
    endif ! end if is_mean0

    ! ===================================================
    ! Read/Write time-averaged variables for T&D BC
    ! ===================================================
    if (is_TamDong) then

       if (operation==WRITE) call write_restartTD(TDfile)
       
       if (operation==READ) call read_restartTD(TDfile,restartTD)
       
    endif ! end if is_TamDong
    
  end subroutine read_write_volume

  !==============================================================================================
  subroutine read_write_vol(binfile,operation,U_var)
  !==============================================================================================
    !> author: XG
    !> date: April 2022
    !> read/write volume for prescribed n_var variables
  !==============================================================================================
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Input/Output arguments
    character(len=*), intent(in) :: binfile
    integer, intent(in) :: operation
    real(wp), dimension(:,:,:,:), target :: U_var ! variable array
    ! -------------------------------------------------------------------------------------------
    ! Local variables
    integer  :: m,ndata,filetype
    real(wp) :: soltime
    ! -------------------------------------------------------------------------------------------

    ! Number & name of datas
    ! ======================
    ndata=size(U_var,4)
    allocate(varlist(ndata),dataname(ndata))
    dataname(1:ndata)='x'

    ! fill varlist (pointer on datas)
    ! ===============================
    do m=1,ndata
       varlist(m)%data=> U_var(:,:,:,m)
       varlist(m)%name= dataname(m)
    enddo

    ! Write volume for datas
    ! ======================
    if (operation==WRITE) then

       ! filetype [0:full; 1:only grid; 2:only solution]
       filetype=2
       ! non-dimensional time
       soltime=tstar

       ! write operation
       if (iproc.eq.0) write(*,*) 'Writing ', trim(binfile)
       call write_tec(binfile,filetype,soltime,field(nob(iproc)))

    endif

    ! Read volume for datas
    ! =====================
    if (operation==READ) then

       ! check file existence
       call mpicheckfile(binfile)

       ! read operation
       if (iproc.eq.0) write(*,*) 'Reading ',trim(binfile)
       call read_tec(binfile,field(nob(iproc)))

    endif

    ! free pointer varlist, dataname
    ! ==============================
    deallocate(varlist,dataname)

  end subroutine read_write_vol

!!$  !==============================================================================================
!!$  subroutine read_write_vol(binfile,operation)
!!$  !==============================================================================================
!!$    !> author: XG
!!$    !> date: April 2021
!!$    !> read/write volume [added for metrics commutations]
!!$  !==============================================================================================
!!$    implicit none
!!$    ! -------------------------------------------------------------------------------------------
!!$    ! Input/Output arguments
!!$    character(len=*), intent(in) :: binfile
!!$    integer, intent(in) :: operation
!!$    ! -------------------------------------------------------------------------------------------
!!$    ! Local variables
!!$    integer  :: m,ndata,filetype
!!$    real(wp) :: soltime
!!$    ! -------------------------------------------------------------------------------------------
!!$
!!$    ! Number & name of datas
!!$    ! ======================
!!$    ndata=4
!!$    allocate(varlist(ndata),dataname(ndata),dummy(nx,ny,nz,ndata))
!!$    dummy=0.0_wp
!!$    dataname(1:4) = (/'a','b','c','d'/)
!!$
!!$    ! fill varlist (pointer on datas)
!!$    ! ===============================
!!$    do m=1,ndata
!!$       varlist(m)%data=> dummy(:,:,:,m)
!!$       varlist(m)%name= dataname(m)
!!$    enddo
!!$
!!$    ! Write volume for metrics commutation
!!$    ! ====================================
!!$    if (operation==WRITE) then
!!$
!!$       ! fill data: double derivatives of metrics
!!$       dummy(:,:,:,1) = x_ksi_eta(1:nx,1:ny,1:nz)
!!$       dummy(:,:,:,2) = y_ksi_eta(1:nx,1:ny,1:nz)
!!$       dummy(:,:,:,3) = x_eta_ksi(1:nx,1:ny,1:nz)
!!$       dummy(:,:,:,4) = y_eta_ksi(1:nx,1:ny,1:nz)
!!$
!!$       ! filetype [0:full; 1:only grid; 2:only solution]
!!$       filetype=2
!!$       ! non-dimensional time
!!$       soltime=tstar
!!$
!!$       ! write double derivatives of metrics
!!$       if (iproc.eq.0) write(*,*) 'Writing restart ~>', trim(binfile)
!!$       call write_tec(binfile,filetype,soltime,field(nob(iproc)))
!!$
!!$    endif ! end write for metrics commutation
!!$
!!$    ! Read volume for metrics commutation
!!$    ! ===================================
!!$    if (operation==READ) then
!!$
!!$       ! check file existence
!!$       call mpicheckfile(binfile)
!!$
!!$       ! read double derivatives of metrics
!!$       if (iproc.eq.0) write(*,*) 'Reading restart ~>',trim(binfile)
!!$       call read_tec(binfile,field(nob(iproc)))
!!$
!!$       ! fill data: double derivatives of metrics
!!$       x_ksi_eta(1:nx,1:ny,1:nz) = dummy(:,:,:,1)
!!$       y_ksi_eta(1:nx,1:ny,1:nz) = dummy(:,:,:,2)
!!$       x_eta_ksi(1:nx,1:ny,1:nz) = dummy(:,:,:,3)
!!$       y_eta_ksi(1:nx,1:ny,1:nz) = dummy(:,:,:,4)
!!$
!!$    endif ! end read for metrics commutation
!!$
!!$    ! free pointer varlist, dataname, dummy
!!$    ! =====================================
!!$    deallocate(varlist,dataname,dummy)
!!$
!!$  end subroutine read_write_vol

  !==============================================================================================
  subroutine read_vol_o(binfile,rho_o,rhou_o,rhov_o,rhow_o,rhoe_o)
  !==============================================================================================
    !> author: XG
    !> date: April 2021
    !> read/write volume for old field [added for interp2d]
  !==============================================================================================
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Input/Output arguments
    character(len=*), intent(in) :: binfile
    real(wp), dimension(:,:,:) :: rho_o,rhou_o,rhov_o,rhow_o,rhoe_o
    ! -------------------------------------------------------------------------------------------
    ! Local variables
    integer  :: m,ndata
    ! -------------------------------------------------------------------------------------------

    ! Restart files store conservative variables in global volume
    ! ===========================================================
    ndata=5
    allocate(varlist(ndata),dataname(ndata),dummy(size(rho_o,1),size(rho_o,2),size(rho_o,3),ndata))
    dummy=0.0_wp
    dataname(1:5)=(/'rro','rou','rov','row','roe'/)

    ! fill varlist (pointer on datas)
    ! ===============================
    do m=1,ndata
       varlist(m)%data=>dummy(:,:,:,m)
       varlist(m)%name= dataname(m)
    enddo

    ! check file existence
    ! ====================
    call mpicheckfile(binfile)

    ! read instantaneous solution to be interpolated
    ! ==============================================
    if (iproc.eq.0) write(*,*) 'Reading old restart ~>', trim(binfile)
    call read_tec(binfile,field_o(nob(iproc)))

    ! fill data: instantaneous conservative variables
    ! ==========
     rho_o= dummy(:,:,:,1)
    rhou_o= dummy(:,:,:,2)
    rhov_o= dummy(:,:,:,3)
    rhow_o= dummy(:,:,:,4)
    rhoe_o= dummy(:,:,:,5)

    ! free pointer varlist, dataname, dummy
    ! =====================================
    deallocate(varlist,dataname,dummy)

  end subroutine read_vol_o

  !==============================================================================================
  subroutine read_vol_o_rans(binfile,rho_o,rhou_o,rhov_o,rhow_o,rhoe_o,var1_rans_o)
  !==============================================================================================
    !> author: XG
    !> date: April 2021
    !> read/write volume for old field [added for interp2d]
  !==============================================================================================
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Input/Output arguments
    character(len=*), intent(in) :: binfile
    real(wp), dimension(:,:,:) :: rho_o,rhou_o,rhov_o,rhow_o,rhoe_o,var1_rans_o
    ! -------------------------------------------------------------------------------------------
    ! Local variables
    integer  :: m,ndata
    ! -------------------------------------------------------------------------------------------

    ! Restart files store conservative variables in global volume
    ! ===========================================================
    if (model_RANS.eq.'SA') then
       ndata=6
       allocate(varlist(ndata),dataname(ndata),dummy(size(rho_o,1),size(rho_o,2),size(rho_o,3),ndata))
       dataname(1:6)=(/'rro','rou','rov','row','roe','nutil'/)
    endif
    dummy=0.0_wp

    ! fill varlist (pointer on datas)
    ! ===============================
    do m=1,ndata
       varlist(m)%data=>dummy(:,:,:,m)
       varlist(m)%name= dataname(m)
    enddo

    ! check file existence
    ! ====================
    call mpicheckfile(binfile)

    ! read instantaneous solution to be interpolated
    ! ==============================================
    if (iproc.eq.0) write(*,*) 'Reading old restart ~>', trim(binfile)
    call read_tec(binfile,field_o(nob(iproc)))

    ! fill data: instantaneous conservative variables
    ! ==========
     rho_o= dummy(:,:,:,1)
    rhou_o= dummy(:,:,:,2)
    rhov_o= dummy(:,:,:,3)
    rhow_o= dummy(:,:,:,4)
    rhoe_o= dummy(:,:,:,5)
    if (model_RANS.eq.'SA') var1_rans_o= dummy(:,:,:,6)

    ! free pointer varlist, dataname, dummy
    ! =====================================
    deallocate(varlist,dataname,dummy)

  end subroutine read_vol_o_rans

  !==============================================================================================
  subroutine read_mean_o(binfile,rho0_o,u0_o,v0_o,w0_o,p0_o,T0_o)
  !==============================================================================================
    !> author: XG
    !> date: April 2021
    !> read/write volume for old field [added for interp2d]
  !==============================================================================================
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Input/Output arguments
    character(len=*), intent(in) :: binfile
    real(wp), dimension(:,:,:) :: rho0_o,u0_o,v0_o,w0_o,p0_o,T0_o
    ! -------------------------------------------------------------------------------------------
    ! Local variables
    integer  :: m,ndata
    ! -------------------------------------------------------------------------------------------

    ! Restart files store conservative variables in global volume
    ! ===========================================================
    ndata=6
    allocate(varlist(ndata),dataname(ndata),dummy(size(rho0_o,1),size(rho0_o,2),size(rho0_o,3),ndata))
    dummy=0.0_wp
    dataname(1:ndata)=(/'u','v','w','p','r','T'/)

    ! fill varlist (pointer on datas)
    ! ===============================
    do m=1,ndata
       varlist(m)%data=>dummy(:,:,:,m)
       varlist(m)%name= dataname(m)
    enddo

    ! check file existence
    ! ====================
    call mpicheckfile(binfile)

    ! read instantaneous solution to be interpolated
    ! ==============================================
    if (iproc.eq.0) write(*,*) 'Reading old mean0 ~>', trim(binfile)
    call read_tec(binfile,field_o(nob(iproc)))

    ! fill data: instantaneous conservative variables
    ! ==========
      u0_o=dummy(:,:,:,1)
      v0_o=dummy(:,:,:,2)
      w0_o=dummy(:,:,:,3)
      p0_o=dummy(:,:,:,4)
    rho0_o=dummy(:,:,:,5)
      T0_o=dummy(:,:,:,6)

    ! free pointer varlist, dataname, dummy
    ! =====================================
    deallocate(varlist,dataname,dummy)

  end subroutine read_mean_o

    !==============================================================================================
  subroutine read_write_grid_cart(gridname,ngh_,operation)
  !==============================================================================================
    !> author: AB
    !> date: March 2024
    !> MPI-IO read/write binary grid for 2D/3D cartesian grid: xg, yg if 2D/ xg, yg, zg if 3D
  !==============================================================================================
    use mod_block
    use mod_utils
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Input/Output arguments
    character(len=*), intent(in) :: gridname
    integer, intent(in) :: operation,ngh_
    ! -------------------------------------------------------------------------------------------
    ! Local variables IO
    integer :: type_mat_cartx,type_mat_viewx
    integer :: type_mat_carty,type_mat_viewy
    integer :: type_mat_cartz,type_mat_viewz
    integer :: ierror,nbytes_dbl,fh
    integer(kind=MPI_OFFSET_KIND) :: offset ! file offset
    integer, dimension(1) :: shape_grx,shape_gry,shape_grz,start_gr
    integer, dimension(1) :: shape_globx,shape_localx,start_localx
    integer, dimension(1) :: shape_globy,shape_localy,start_localy
    integer, dimension(1) :: shape_globz,shape_localz,start_localz
    integer, dimension(MPI_STATUS_SIZE) :: statut
    ! Local variables read/write
    logical  :: iexist
    ! -------------------------------------------------------------------------------------------

    ! Get size of the type MPI_DOUBLE_PRECISION
    call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,nbytes_dbl,ierror)

    ! Creation of the derived type type_local that defines the array with/without ghost cells
    ! ======================================================================================
    if ((operation==WRITE_LEADER).or.(operation==READ_LEADER)) then
       ! Shape of x, y & z array with ghost cells
       shape_grx= (/ ngx+2*ngh_ /); shape_gry= (/ ngy+2*ngh_ /); shape_grz= (/ ngz+2*ngh_ /)
    else if ((operation==READ).or.(operation==WRITE)) then
       ! Shape of xc & yc array with ghost cells
       shape_grx= (/ nx+2*ngh_ /); shape_gry= (/ ny+2*ngh_ /); shape_grz= (/ nz+2*ngh_ /)
    else
       call mpistop("Unknown type of operation in read_write_grid_cart in mod_io.f90",1)
    endif
    ! Starting coordinates for array (directly taken into account in MPI_FILE_SET_VIEW)
    start_gr= (/ 0 /)

    ! Creation & commit of derived type type_mat_cart
    ! -----------------------------------------------
    ! type_mat_cartx
    call MPI_TYPE_CREATE_SUBARRAY(1,shape_grx,shape_grx,start_gr,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,type_mat_cartx,ierror)
    call MPI_TYPE_COMMIT(type_mat_cartx,ierror)
    ! type_mat_carty
    call MPI_TYPE_CREATE_SUBARRAY(1,shape_gry,shape_gry,start_gr,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,type_mat_carty,ierror)
    call MPI_TYPE_COMMIT(type_mat_carty,ierror)
    if (nz.gt.1) then
       ! Creation & commit of derived type type_mat_cartz
       ! ------------------------------------------------
       call MPI_TYPE_CREATE_SUBARRAY(1,shape_grz,shape_grz,start_gr,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,type_mat_cartz,ierror)
       call MPI_TYPE_COMMIT(type_mat_cartz,ierror)
    endif

    ! Creation of type type_mat_nogh_view to set the view on the file
    ! ===============================================================
    ! Shape of the global array for cartesian grids
    shape_globx=(/ ngx+2*ngh_ /);  shape_globy=(/ ngy+2*ngh_ /);  shape_globz=(/ ngz+2*ngh_ /)
    if ((operation==WRITE_LEADER).or.(operation==READ_LEADER)) then
       ! Shape of the local array for cartesian grids
       shape_localx=(/ ngx+2*ngh_ /); shape_localy=(/ ngy+2*ngh_ /);  shape_localz=(/ ngz+2*ngh_ /)
       ! Starting coordinates for local array
       start_localx=(/ 0 /); start_localy=(/ 0 /); start_localz=(/ 0 /)
    else
       ! Shape of the local array for cartesian grids
       shape_localx=(/ nx+2*ngh_ /); shape_localy=(/ ny+2*ngh_ /);  shape_localz=(/ nz+2*ngh_ /)
       ! Starting coordinates for local array
       start_localx=(/ coord(1)*nx /); start_localy=(/ coord(2)*ny /); start_localz=(/ coord(3)*nz /)
    endif

    ! Creation & commit of derived type type_mat_view
    ! -----------------------------------------------
    ! type_mat_viewx
    call MPI_TYPE_CREATE_SUBARRAY(1,shape_globx,shape_localx,start_localx,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,type_mat_viewx,ierror)
    call MPI_TYPE_COMMIT(type_mat_viewx,ierror)
    ! type_mat_viewy
    call MPI_TYPE_CREATE_SUBARRAY(1,shape_globy,shape_localy,start_localy,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,type_mat_viewy,ierror)
    call MPI_TYPE_COMMIT(type_mat_viewy,ierror)
    if (nz.gt.1) then
       ! Creation & commit of derived type type_mat_viewz
       ! ------------------------------------------------
       call MPI_TYPE_CREATE_SUBARRAY(1,shape_globz,shape_localz,start_localz,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,type_mat_viewz,ierror)
       call MPI_TYPE_COMMIT(type_mat_viewz,ierror)
    endif

    ! ===========================================
    ! Read/Write GRID (with ngh_ extended points)
    ! ===========================================
    ! Create temporary communicator
    ! -----------------------------
    if ((operation==WRITE_LEADER).or.(operation==READ_LEADER)) then
       call MPI_COMM_SPLIT(COMM_global,iproc,iproc,COMM_rw_grid,info)
    else
       COMM_rw_grid=COMM_intrablock
    endif

    ! Open file for reading/writing
    ! -----------------------------
    if ((operation==WRITE_LEADER).or.(operation==WRITE)) then
       call MPI_FILE_OPEN(COMM_rw_grid,trim(gridname),MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierror)
    else
       call MPI_FILE_OPEN(COMM_rw_grid,trim(gridname),MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierror)
    endif
    if (ierror/=MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_cart in opening gridfile '// trim(gridname),1)

    ! Reading/Writting x
    ! ------------------
    ! Setting the view on the file
    call MPI_FILE_SET_VIEW(fh,0_MPI_OFFSET_KIND,MPI_DOUBLE_PRECISION,type_mat_viewx,'native',MPI_INFO_NULL,ierror)
    if (operation==WRITE_LEADER) then
       ! Writting x
       call MPI_FILE_WRITE_ALL(fh,xg(1-ngh_:ngx+ngh_),1,type_mat_cartx,statut,ierror)
       if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_cart in writing gridfile '//trim(gridname),1)
    else if (operation==READ_LEADER) then
       ! Reading x
       call MPI_FILE_READ_ALL(fh,xg(1-ngh_:ngx+ngh_),1,type_mat_cartx,statut,ierror)
       if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_cart in reading gridfile '//trim(gridname),1)
    else if (operation==WRITE) then
       ! Writting x
       call MPI_FILE_WRITE_ALL(fh,x(1-ngh_:nx+ngh_),1,type_mat_cartx,statut,ierror)
       if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_cart in writing gridfile '//trim(gridname),1)
    else
       ! Reading x
       call MPI_FILE_READ_ALL(fh,x(1-ngh_:nx+ngh_),1,type_mat_cartx,statut,ierror)
       if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_cart in reading gridfile '//trim(gridname),1)
    endif

    ! Reading/Writting y
    ! ------------------
    ! Setting the offset
    offset = nbytes_dbl*(ngx+2*ngh_)
    ! Setting the view on the file
    call MPI_FILE_SET_VIEW(fh,offset,MPI_DOUBLE_PRECISION,type_mat_viewy,'native',MPI_INFO_NULL,ierror)
    if (operation==WRITE_LEADER) then
       ! Writting y
       call MPI_FILE_WRITE_ALL(fh,yg(1-ngh_:ngy+ngh_),1,type_mat_carty,statut,ierror)
       if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_cart in writing gridfile '//trim(gridname),1)
    else if (operation==READ_LEADER) then
       ! Reading y
       call MPI_FILE_READ_ALL(fh,yg(1-ngh_:ngy+ngh_),1,type_mat_carty,statut,ierror)
       if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_cart in reading gridfile '//trim(gridname),1)
    else if (operation==WRITE) then
       ! Writting y
       call MPI_FILE_WRITE_ALL(fh,y(1-ngh_:ny+ngh_),1,type_mat_carty,statut,ierror)
       if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_cart in writing gridfile '//trim(gridname),1)
    else
       ! Reading y
       call MPI_FILE_READ_ALL(fh,y(1-ngh_:ny+ngh_),1,type_mat_carty,statut,ierror)
       if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_cart in reading gridfile '//trim(gridname),1)
    endif

    ! Reading/Writting z
    ! ------------------
    if (nz.gt.1) then
       ! Setting the offset
       offset = nbytes_dbl*(ngx+2*ngh_) + nbytes_dbl*(ngy+2*ngh_)
       ! Setting the view on the file
       call MPI_FILE_SET_VIEW(fh,offset,MPI_DOUBLE_PRECISION,type_mat_viewz,'native',MPI_INFO_NULL,ierror)
       if (operation==WRITE_LEADER) then
          ! Writting z
          call MPI_FILE_WRITE_ALL(fh,zg(1-ngh_:ngz+ngh_),1,type_mat_cartz,statut,ierror)
          if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_cart in writing gridfile '//trim(gridname),1)
       else if (operation==READ_LEADER) then
          ! Reading z
          call MPI_FILE_READ_ALL(fh,zg(1-ngh_:ngz+ngh_),1,type_mat_cartz,statut,ierror)
          if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_cart in reading gridfile '//trim(gridname),1)
       else if (operation==WRITE) then
          ! Writting z
          call MPI_FILE_WRITE_ALL(fh,z(1-ngh_:nz+ngh_),1,type_mat_cartz,statut,ierror)
          if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_cart in writing gridfile '//trim(gridname),1)
       else
          ! Reading z
          call MPI_FILE_READ_ALL(fh,z(1-ngh_:nz+ngh_),1,type_mat_cartz,statut,ierror)
          if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_cart in reading gridfile '//trim(gridname),1)
       endif
    endif

    ! Close file
    ! ----------
    call MPI_FILE_CLOSE(fh,ierror) ! Close file

    ! Free types for MPI-IO of grid
    ! =============================
    call MPI_TYPE_FREE(type_mat_viewx,info)
    call MPI_TYPE_FREE(type_mat_viewy,info)
    call MPI_TYPE_FREE(type_mat_cartx,info); call MPI_TYPE_FREE(type_mat_carty,info)
    if (nz.gt.1) then
       call MPI_TYPE_FREE(type_mat_viewz,info); call MPI_TYPE_FREE(type_mat_cartz,info)
    endif

  end subroutine read_write_grid_cart

  !==============================================================================================
  subroutine read_write_grid_curv(gridname,ngh_,operation)
  !==============================================================================================
    !> author: AB
    !> date: February 2024
    !> MPI-IO read/write binary grid for 2D curvilinear grid: xgc, ygc if 2D/ xgc, ygc, zg if 3D
  !==============================================================================================
    use mod_block
    use mod_utils
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Input/Output arguments
    character(len=*), intent(in) :: gridname
    integer, intent(in) :: operation,ngh_
    ! -------------------------------------------------------------------------------------------
    ! Local variables IO
    integer :: type_mat_curv,type_mat_cart,type_mat_view1,type_mat_view2
    integer :: ierror,nbytes_dbl,fh
    integer(kind=MPI_OFFSET_KIND) :: offset ! file offset
    integer, dimension(1) :: shape_gr2,start_gr2
    integer, dimension(1) :: shape_glob2,shape_local2,start_local2
    integer, dimension(2) :: shape_gr1,start_gr1
    integer, dimension(2) :: shape_glob1,shape_local1,start_local1
    integer, dimension(MPI_STATUS_SIZE) :: statut
    ! Local variables read/write
    logical  :: iexist
    ! -------------------------------------------------------------------------------------------

    ! Get size of the type MPI_DOUBLE_PRECISION
    call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,nbytes_dbl,ierror)

    ! Creation of the derived type type_local that defines the array with/without ghost cells
    ! ======================================================================================
    if ((operation==WRITE_LEADER).or.(operation==READ_LEADER)) then
       ! Shape of xc & yc array with ghost cells
       shape_gr1= (/ ngx+2*ngh_, ngy+2*ngh_ /)
       ! Shape of z array with ghost cells
       shape_gr2= (/ ngz+2*ngh_ /)
    else if ((operation==READ).or.(operation==WRITE)) then
       ! Shape of xc & yc array with ghost cells
       shape_gr1= (/ nx+2*ngh_, ny+2*ngh_ /)
       ! Shape of z array with ghost cells
       shape_gr2= (/ nz+2*ngh_ /)
    else
       call mpistop("Unknown type of operation in read_write_grid_curv in mod_io.f90",1)
    endif
    ! Starting coordinates for array (directly taken into account in MPI_FILE_SET_VIEW)
    start_gr1= (/ 0, 0 /)
    start_gr2= (/ 0 /)

    ! Creation of derived type type_mat_curv
    ! --------------------------------------
    call MPI_TYPE_CREATE_SUBARRAY(2,shape_gr1,shape_gr1,start_gr1,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,type_mat_curv,ierror)
    ! Commit of type_mat_curv
    call MPI_TYPE_COMMIT(type_mat_curv,ierror)
    if (nz.gt.1) then
       ! Creation of derived type type_mat_cart
       ! --------------------------------------
       call MPI_TYPE_CREATE_SUBARRAY(1,shape_gr2,shape_gr2,start_gr2,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,type_mat_cart,ierror)
       ! Commit of type_mat_cart
       call MPI_TYPE_COMMIT(type_mat_cart,ierror)
    endif

    ! Creation of type type_mat_nogh_view to set the view on the file
    ! ===============================================================
    if ((operation==WRITE_LEADER).or.(operation==READ_LEADER)) then
       ! Shape of the global & local array for 2D curvilinear grid
       shape_glob1=(/ ngx+2*ngh_, ngy+2*ngh_ /); shape_local1=(/ ngx+2*ngh_, ngy+2*ngh_ /)
       ! Shape of the global & local array for spanwise cartesian grid
       shape_glob2=(/ ngz+2*ngh_ /); shape_local2=(/ ngz+2*ngh_ /)
       ! Starting coordinates for local array
       start_local1=(/ 0, 0 /); start_local2=(/ 0 /)
    else
       ! Shape of the global & local array for 2D curvilinear grid
       shape_glob1=(/ ngx+2*ngh_, ngy+2*ngh_ /); shape_local1=(/ nx+2*ngh_, ny+2*ngh_ /)
       ! Shape of the global & local array for spanwise cartesian grid
       shape_glob2=(/ ngz+2*ngh_ /); shape_local2=(/ nz+2*ngh_ /)
       ! Starting coordinates for local array
       start_local1=(/ coord(1)*nx, coord(2)*ny /); start_local2=(/ coord(3)*nz /)
    endif

    ! Creation of derived type type_mat_view1
    ! ---------------------------------------
    call MPI_TYPE_CREATE_SUBARRAY(2,shape_glob1,shape_local1,start_local1,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,type_mat_view1,ierror)
    ! Commit of type_mat_view1
    call MPI_TYPE_COMMIT(type_mat_view1,ierror)
    ! Creation of derived type type_mat_view2
    ! ---------------------------------------
    if (nz.gt.1) then
       call MPI_TYPE_CREATE_SUBARRAY(1,shape_glob2,shape_local2,start_local2,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,type_mat_view2,ierror)
       ! Commit of type_mat_view2
       call MPI_TYPE_COMMIT(type_mat_view2,ierror)
    endif

    ! ===========================================
    ! Read/Write GRID (with ngh_ extended points)
    ! ===========================================
    ! Create temporary communicator
    ! -----------------------------
    if ((operation.eq.WRITE_LEADER).or.(operation==READ_LEADER)) then
       call MPI_COMM_SPLIT(COMM_global,iproc,iproc,COMM_rw_grid,info)
    else
       COMM_rw_grid=COMM_intrablock
    endif



    ! Open file for reading/writing
    ! -----------------------------
    if ((operation==WRITE_LEADER).or.(operation==WRITE)) then
       call MPI_FILE_OPEN(COMM_rw_grid,trim(gridname),MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierror)
    else
       call MPI_FILE_OPEN(COMM_rw_grid,trim(gridname),MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierror)
    endif
    if (ierror/=MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_curv in opening gridfile '// trim(gridname),1)

    ! Reading/Writting xc & yc
    ! ------------------------
    ! Setting the view on the file
    call MPI_FILE_SET_VIEW(fh,0_MPI_OFFSET_KIND,MPI_DOUBLE_PRECISION,type_mat_view1,'native',MPI_INFO_NULL,ierror)
    if (operation==WRITE_LEADER) then
       ! Writting xc
       call MPI_FILE_WRITE_ALL(fh,xgce(1-ngh_:ngx+ngh_,1-ngh_:ngy+ngh_),1,type_mat_curv,statut,ierror)
       if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_curv in writing gridfile '//trim(gridname),1)
       ! Writting yc
       call MPI_FILE_WRITE_ALL(fh,ygce(1-ngh_:ngx+ngh_,1-ngh_:ngy+ngh_),1,type_mat_curv,statut,ierror)
       if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_curv in writing gridfile '//trim(gridname),1)
    else if (operation==READ_LEADER) then
       ! Reading xc
       call MPI_FILE_READ_ALL(fh,xgce(1-ngh_:ngx+ngh_,1-ngh_:ngy+ngh_),1,type_mat_curv,statut,ierror)
       if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_curv in reading gridfile '//trim(gridname),1)
       ! Reading yc
       call MPI_FILE_READ_ALL(fh,ygce(1-ngh_:ngx+ngh_,1-ngh_:ngy+ngh_),1,type_mat_curv,statut,ierror)
       if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_curv in reading gridfile '//trim(gridname),1)
    else if (operation==WRITE) then
       ! Writting xc
       call MPI_FILE_WRITE_ALL(fh,xc(1-ngh_:nx+ngh_,1-ngh_:ny+ngh_),1,type_mat_curv,statut,ierror)
       if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_curv in writing gridfile '//trim(gridname),1)
       ! Writting yc
       call MPI_FILE_WRITE_ALL(fh,yc(1-ngh_:nx+ngh_,1-ngh_:ny+ngh_),1,type_mat_curv,statut,ierror)
       if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_curv in writing gridfile '//trim(gridname),1)
    else
       ! Reading xc
       call MPI_FILE_READ_ALL(fh,xc(1-ngh_:nx+ngh_,1-ngh_:ny+ngh_),1,type_mat_curv,statut,ierror)
       if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_curv in reading gridfile '//trim(gridname),1)
       ! Reading yc
       call MPI_FILE_READ_ALL(fh,yc(1-ngh_:nx+ngh_,1-ngh_:ny+ngh_),1,type_mat_curv,statut,ierror)
       if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_curv in reading gridfile '//trim(gridname),1)
    endif

    ! Reading/Writting z
    ! ------------------
    if (nz.gt.1) then
       ! Setting the offset
       offset = nbytes_dbl*(ngx+2*ngh_)*(ngy+2*ngh_)*2
       ! Setting the view on the file
       call MPI_FILE_SET_VIEW(fh,offset,MPI_DOUBLE_PRECISION,type_mat_view2,'native',MPI_INFO_NULL,ierror)
       if (operation==WRITE_LEADER) then
          ! Writting z
          call MPI_FILE_WRITE_ALL(fh,zg(1-ngh_:ngz+ngh_),1,type_mat_cart,statut,ierror)
          if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_curv in writing gridfile '//trim(gridname),1)
       else if (operation==WRITE_LEADER) then
          ! Reading z
          call MPI_FILE_READ_ALL(fh,zg(1-ngh_:ngz+ngh_),1,type_mat_cart,statut,ierror)
          if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_curv in reading gridfile '//trim(gridname),1)
       else if (operation==WRITE) then
          ! Writting z
          call MPI_FILE_WRITE_ALL(fh,z(1-ngh_:nz+ngh_),1,type_mat_cart,statut,ierror)
          if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_curv in writing gridfile '//trim(gridname),1)
       else
          ! Reading z
          call MPI_FILE_READ_ALL(fh,z(1-ngh_:nz+ngh_),1,type_mat_cart,statut,ierror)
          if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_curv in reading gridfile '//trim(gridname),1)
       endif
    endif

    ! Close file
    ! ----------
    call MPI_FILE_CLOSE(fh,ierror) ! Close file

    ! Free types for MPI-IO of grid
    ! =============================
    call MPI_TYPE_FREE(type_mat_curv,info); call MPI_TYPE_FREE(type_mat_view1,info)
    if (nz.gt.1) then
       call MPI_TYPE_FREE(type_mat_cart,info); call MPI_TYPE_FREE(type_mat_view2,info)
    endif

  end subroutine read_write_grid_curv

  !==============================================================================================
  subroutine read_write_grid_curv3(gridname,ngh_,operation)
  !==============================================================================================
    !> author: AB
    !> date: March 2024
    !> MPI-IO read/write binary grid for 3D curvilinear grid: xgc3, ygc3, zgc3
  !==============================================================================================
    use mod_block
    use mod_constant ! <- for is_read_ex
    use mod_utils
    ! use mod_grid
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Input/Output arguments
    character(len=*), intent(in) :: gridname
    integer, intent(in) :: operation,ngh_
    ! -------------------------------------------------------------------------------------------
    ! Local variables IO
    integer :: type_mat_curv3,type_mat_view
    integer :: ierror,fh
    integer, dimension(3) :: shape_gr,start_gr
    integer, dimension(3) :: shape_glob,shape_local,start_local
    integer, dimension(MPI_STATUS_SIZE) :: statut
    ! Local variables read/write
    logical  :: iexist
    ! -------------------------------------------------------------------------------------------

    ! Creation of the derived type type_local that defines the array with/without ghost cells
    ! ======================================================================================
    if ((operation==WRITE_LEADER).or.(operation==READ_LEADER)) then
       ! Shape of xc3, yc3 & zc3 array with ghost cells
       shape_gr= (/ ngx+2*ngh_, ngy+2*ngh_, ngz+2*ngh_ /)
    else if ((operation==READ).or.(operation==WRITE)) then
       ! Shape of xc3, yc3 & zc3 array with ghost cells
       shape_gr= (/ nx+2*ngh_, ny+2*ngh_ , nz+2*ngh_ /)
    else
       call mpistop("Unknown type of operation in read_write_grid_curv3 in mod_io.f90",1)
    endif
    ! Starting coordinates for array (directly taken into account in MPI_FILE_SET_VIEW)
    start_gr= (/ 0, 0 , 0 /)

    ! Creation of derived type type_mat_curv3
    ! ---------------------------------------
    call MPI_TYPE_CREATE_SUBARRAY(3,shape_gr,shape_gr,start_gr,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,type_mat_curv3,ierror)
    ! Commit of type_mat_curv3
    call MPI_TYPE_COMMIT(type_mat_curv3,ierror)

    ! Creation of type type_mat_nogh_view to set the view on the file
    ! ===============================================================
    if ((operation==WRITE_LEADER).or.(operation==READ_LEADER)) then
       ! Shape of the global & local array for 2D curvilinear grid
       shape_glob=(/ ngx+2*ngh_, ngy+2*ngh_, ngz+2*ngh_ /)
       shape_local=(/ ngx+2*ngh_, ngy+2*ngh_, ngz+2*ngh_ /)
       ! Starting coordinates for local array
       start_local=(/ 0, 0, 0 /)
    else
       ! Shape of the global & local array for 2D curvilinear grid
       shape_glob=(/ ngx+2*ngh_, ngy+2*ngh_, ngz+2*ngh_ /)
       shape_local=(/ nx+2*ngh_, ny+2*ngh_, nz+2*ngh_ /)
       ! Starting coordinates for local array
       start_local=(/ coord(1)*nx, coord(2)*ny, coord(3)*nz /)
    endif

    ! Creation of derived type type_mat_view
    ! --------------------------------------
    call MPI_TYPE_CREATE_SUBARRAY(3,shape_glob,shape_local,start_local,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,type_mat_view,ierror)
    ! Commit of type_mat_view
    call MPI_TYPE_COMMIT(type_mat_view,ierror)

    ! ===========================================
    ! Read/Write GRID (with ngh_ extended points)
    ! ===========================================
    ! Create temporary communicator
    ! -----------------------------
    if ((operation==WRITE_LEADER).or.(operation==READ_LEADER)) then
       call MPI_COMM_SPLIT(COMM_global,iproc,iproc,COMM_rw_grid,info)
    else
       COMM_rw_grid=COMM_intrablock
    endif

    ! Open file for reading/writing
    ! -----------------------------
    if ((operation==WRITE_LEADER).or.(operation==WRITE)) then
       call MPI_FILE_OPEN(COMM_rw_grid,trim(gridname),MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierror)
    else
       call MPI_FILE_OPEN(COMM_rw_grid,trim(gridname),MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierror)
    endif
    if (ierror/=MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_curv3 in opening gridfile '// trim(gridname),1)

    ! Reading/Writting xc3, yc3 & zc3
    ! -------------------------------
    ! Setting the view on the file
    call MPI_FILE_SET_VIEW(fh,0_MPI_OFFSET_KIND,MPI_DOUBLE_PRECISION,type_mat_view,'native',MPI_INFO_NULL,ierror)
    if (operation==WRITE_LEADER) then
       ! Writting xc3
       call MPI_FILE_WRITE_ALL(fh,xgc3e(1-ngh_:ngx+ngh_,1-ngh_:ngy+ngh_,1-ngh_:ngz+ngh_),1,type_mat_curv3,statut,ierror)
       if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_curv3 in writing gridfile '//trim(gridname),1)
       ! Writting yc3
       call MPI_FILE_WRITE_ALL(fh,ygc3e(1-ngh_:ngx+ngh_,1-ngh_:ngy+ngh_,1-ngh_:ngz+ngh_),1,type_mat_curv3,statut,ierror)
       if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_curv3 in writing gridfile '//trim(gridname),1)
       ! Writting zc3
       call MPI_FILE_WRITE_ALL(fh,zgc3e(1-ngh_:ngx+ngh_,1-ngh_:ngy+ngh_,1-ngh_:ngz+ngh_),1,type_mat_curv3,statut,ierror)
       if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_curv3 in writing gridfile '//trim(gridname),1)
    else if (operation==READ_LEADER) then
       ! Reading xc3
       call MPI_FILE_READ_ALL(fh,xgc3e(1-ngh_:ngx+ngh_,1-ngh_:ngy+ngh_,1-ngh_:ngz+ngh_),1,type_mat_curv3,statut,ierror)
       if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_curv3 in reading gridfile '//trim(gridname),1)
       ! Reading yc3
       call MPI_FILE_READ_ALL(fh,ygc3e(1-ngh_:ngx+ngh_,1-ngh_:ngy+ngh_,1-ngh_:ngz+ngh_),1,type_mat_curv3,statut,ierror)
       if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_curv3 in reading gridfile '//trim(gridname),1)
       ! Reading zc3
       call MPI_FILE_READ_ALL(fh,zgc3e(1-ngh_:ngx+ngh_,1-ngh_:ngy+ngh_,1-ngh_:ngz+ngh_),1,type_mat_curv3,statut,ierror)
       if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_curv3 in reading gridfile '//trim(gridname),1)
    else if (operation==WRITE) then
       ! Writting xc3
       call MPI_FILE_WRITE_ALL(fh,xc3(1-ngh_:nx+ngh_,1-ngh_:ny+ngh_,1-ngh_:nz+ngh_),1,type_mat_curv3,statut,ierror)
       if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_curv3 in writing gridfile '//trim(gridname),1)
       ! Writting yc3
       call MPI_FILE_WRITE_ALL(fh,yc3(1-ngh_:nx+ngh_,1-ngh_:ny+ngh_,1-ngh_:nz+ngh_),1,type_mat_curv3,statut,ierror)
       if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_curv3 in writing gridfile '//trim(gridname),1)
       ! Writting zc3
       call MPI_FILE_WRITE_ALL(fh,zc3(1-ngh_:nx+ngh_,1-ngh_:ny+ngh_,1-ngh_:nz+ngh_),1,type_mat_curv3,statut,ierror)
       if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_curv3 in writing gridfile '//trim(gridname),1)
    else
       ! Reading xc3
       call MPI_FILE_READ_ALL(fh,xc3(1-ngh_:nx+ngh_,1-ngh_:ny+ngh_,1-ngh_:nz+ngh_),1,type_mat_curv3,statut,ierror)
       if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_curv3 in reading gridfile '//trim(gridname),1)
       ! Reading yc3
       call MPI_FILE_READ_ALL(fh,yc3(1-ngh_:nx+ngh_,1-ngh_:ny+ngh_,1-ngh_:nz+ngh_),1,type_mat_curv3,statut,ierror)
       if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_curv3 in reading gridfile '//trim(gridname),1)
       ! Reading zc3
       call MPI_FILE_READ_ALL(fh,zc3(1-ngh_:nx+ngh_,1-ngh_:ny+ngh_,1-ngh_:nz+ngh_),1,type_mat_curv3,statut,ierror)
       if (ierror /= MPI_SUCCESS) call mpistop('MPI_IO: Error in read_write_grid_curv3 in reading gridfile '//trim(gridname),1)
    endif

    ! Close file
    ! ----------
    call MPI_FILE_CLOSE(fh,ierror) ! Close file

    ! Free types for MPI-IO of grid
    ! =============================
    call MPI_TYPE_FREE(type_mat_curv3,info); call MPI_TYPE_FREE(type_mat_view,info)


  end subroutine read_write_grid_curv3


  !==============================================================================================
  subroutine init_io_grid_ex3d
  !==============================================================================================
    !> author:
    !> date: October 2023
    !> Initialize the sub-communicator an I/O-type for extended grid inputs/outputs
  !==============================================================================================
    use mod_block
    use mod_constant ! <- for is_read_ex
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Local variables
    integer :: ierror,nbytes_dbl
    integer, dimension(3) :: shape_gr,start_gr
    integer, dimension(3) :: shape_glob,shape_local,start_local
    ! -------------------------------------------------------------------------------------------

    allocate(grid_ex_(nbloc))

    grid_ex_(nob(iproc))%MPI_COMM = COMM_intrablock

    ! Get size of the type MPI_DOUBLE_PRECISION
    call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,nbytes_dbl,ierror)

    ! Creation of the derived type type_local that defines the array with/without ghost cells
    ! ======================================================================================
    if (is_read_ex) then
       ! Shape of the array with ghost cells
       shape_gr= (/ nx+2*ngh, ny+2*ngh, nz+2*ngh /)
    else
       ! Shape of the array without ghost cells
       shape_gr= (/ nx, ny, nz /)
    endif
    ! Starting coordinates for array without ghost cells (w.r.t. global array)
    start_gr= (/ 0, 0, 0 /)

    ! Creation of derived type type_mat_nogh
    ! --------------------------------------
    call MPI_TYPE_CREATE_SUBARRAY(3,shape_gr,shape_gr,start_gr,&
         MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,grid_ex_(nob(iproc))%type_mat_nogh,ierror)
    ! Commit of type_mat_nogh
    call MPI_TYPE_COMMIT(grid_ex_(nob(iproc))%type_mat_nogh,ierror)

    ! Creation of type type_mat_nogh_view to set the view on the file
    ! ===============================================================

    if (is_read_ex) then
       ! Shape of the global array
       shape_glob= (/ ngx+2*ngh, ngy+2*ngh, ngz+2*ngh /)
       ! Shape of the local array
       shape_local= (/ nx+2*ngh, ny+2*ngh, nz+2*ngh /)
    else
       ! Shape of the global array
       shape_glob= (/ ngx, ngy, ngz /)
       ! Shape of the local array
       shape_local= (/ nx, ny, nz /)
    endif

    ! Starting coordinates for local array
    ! ------------------------------------
    start_local= (/ coord(1)*nx, coord(2)*ny, coord(3)*nz /)
    if (nz==1) start_local(3) = 0

    ! Creation of derived type type_mat_view
    ! --------------------------------------
    call MPI_TYPE_CREATE_SUBARRAY(3,shape_glob,shape_local,start_local,&
         MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,grid_ex_(nob(iproc))%type_mat_view,ierror)
    ! Commit of type_mat_view
    call MPI_TYPE_COMMIT(grid_ex_(nob(iproc))%type_mat_view,ierror)

    ! Initialize grid_ex_ structure
    ! =============================

    ! Initialize dimensions
    if (is_read_ex) then
       grid_ex_(nob(iproc))%ngx=ngx+2*ngh
       grid_ex_(nob(iproc))%ngy=ngy+2*ngh
       grid_ex_(nob(iproc))%ngz=ngz+2*ngh
       grid_ex_(nob(iproc))%nx =nx +2*ngh
       grid_ex_(nob(iproc))%ny =ny +2*ngh
       grid_ex_(nob(iproc))%nz =nz +2*ngh
    else
       grid_ex_(nob(iproc))%ngx=ngx
       grid_ex_(nob(iproc))%ngy=ngy
       grid_ex_(nob(iproc))%ngz=ngz
       grid_ex_(nob(iproc))%nx =nx
       grid_ex_(nob(iproc))%ny =ny
       grid_ex_(nob(iproc))%nz =nz
    endif

    ! Initialize coord
    if (.not.allocated(grid_ex_(nob(iproc))%coord)) allocate(grid_ex_(nob(iproc))%coord(3))
    grid_ex_(nob(iproc))%coord=coord

    ! Initialize StrandID [-2 for Static StrandID in tecplot binary]
    grid_ex_(nob(iproc))%strandID=-2

    ! Initialize IO type (bin or tec)
    grid_ex_(nob(iproc))%is_IOtec_read =is_IOtec_read
    grid_ex_(nob(iproc))%is_IOtec_write=is_IOtec_write

    ! Initialize attributes to append files
    grid_ex_(nob(iproc))%is_app=.false.
    grid_ex_(nob(iproc))%disp=0
    grid_ex_(nob(iproc))%restart=.false.
    grid_ex_(nob(iproc))%offset=0_MPI_OFFSET_KIND

    if ((iproc==0).and.(verbose)) then
       if (is_read_ex) then
          print *,'init I/O GRID EXTENDED OK'
       else
          print *,'init I/O GRID OK'
       endif
    endif

  end subroutine init_io_grid_ex3d

  !==============================================================================================
  subroutine read_write_grid3d(gridname,operation)
  !==============================================================================================
    !> author: XG
    !> date: October 2023
    !> MPI-IO read/write binary grid
  !==============================================================================================
    use mod_utils
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Input/Output arguments
    character(len=*), intent(in) :: gridname
    integer, intent(in) :: operation
    ! -------------------------------------------------------------------------------------------
    ! Local variables
    integer  :: i,j,k,m,ndata,filetype
    integer :: nex,ney,nez,ngh_
    logical  :: iexist
    real(wp) :: soltime
    ! -------------------------------------------------------------------------------------------

    ! =========================================
    ! Write/Read GRID (with/without extended points)
    ! =========================================

    ! Prepare number and name of data: x,y,z coordinates
    ! ===============================
    ndata=3
    allocate(dataname(ndata),varlist(ndata))
    dataname(1:ndata) = (/'X','Y','Z'/)

    ! Write volume grid :
    ! ===================
    if (operation==WRITE) then

       ! check if grid file already exists
       ! ---------------------------------
       inquire(file=gridname,exist=iexist)

       ! if not: write GRID
       ! ------------------
       if ((.not.iexist).or.(idepart==FROM_SCRATCH)) then

          ! if extended grid with ghost points
          if (is_read_ex) then
             ngh_=ngh
          else
             ngh_=0
          endif

          ! extended sizes (+ ghost cells)
          nex=nx+2*ngh_
          ney=ny+2*ngh_
          nez=nz+2*ngh_

          ! allocate and initialize dummy variable for datas
          allocate(dummy(nex,ney,nez,ndata))
          dummy = 0.0_wp

          ! fill data x,y,z in dummy array
          if (is_curv3) then
             do k=1,nez
                do j=1,ney
                   do i=1,nex
                      dummy(i,j,k,1)=xc3(i-ngh_,j-ngh_,k-ngh_)
                      dummy(i,j,k,2)=yc3(i-ngh_,j-ngh_,k-ngh_)
                      dummy(i,j,k,3)=zc3(i-ngh_,j-ngh_,k-ngh_)
                   enddo
                enddo
             enddo
          else
             if (is_curv) then
                do j=1,ney
                   do i=1,nex
                      dummy(i,j,:,1)=xc(i-ngh_,j-ngh_)
                      dummy(i,j,:,2)=yc(i-ngh_,j-ngh_)
                   enddo
                enddo
             else
                do i=1,nex
                   dummy(i,:,:,1)=x(i-ngh_)
                enddo
                do j=1,ney
                   dummy(:,j,:,2)=y(j-ngh_)
                enddo
             endif
             do k=1,nez
                dummy(:,:,k,3)=z(k-ngh_)
             enddo
          endif

          ! fill varlist (pointer on datas)
          do m=1,ndata
             varlist(m)%data=>dummy(:,:,:,m)
             varlist(m)%name= dataname(m)
          enddo

          ! filetype [0:full; 1:only grid; 2:only solution]
          filetype=1
          ! non-dimensional time
          soltime=tstar

          ! write grid block in separate file
          if (iproc.eq.0) write(*,*) 'Writing grid file ~>',trim(gridname)
          call write_tec(gridname,filetype,soltime,grid_ex_(nob(iproc)))

          ! free temporary data
          deallocate(dummy)

       endif ! end if test iexist

    endif ! end write grid if tecplot format

    ! Read volume grid : done only if tecplot format is selected
    ! ==================
    if (operation==READ) then

       ! if extended grid with ghost points
       if (is_read_ex) then
          ngh_=ngh
       else
          ngh_=0
       endif

       ! extended sizes (+ ghost cells)
       nex=nx+2*ngh_
       ney=ny+2*ngh_
       nez=nz+2*ngh_

       ! allocate and initialize dummy variable for datas
       allocate(dummy(nex,ney,nez,ndata))
       dummy = 0.0_wp

       ! fill varlist (pointer on datas)
       do m=1,ndata
          varlist(m)%data=>dummy(:,:,:,m)
          varlist(m)%name= dataname(m)
       enddo

       ! check file gridfile
       call mpicheckfile(gridname)

       ! read grid block
       if (iproc.eq.0) write(*,*) 'Reading grid file ~>', trim(gridname)
       call read_tec(gridname,grid_ex_(nob(iproc)))

       ! fill data x,y,z from dummy array (Att. CARTESIAN only) TO BE CHANGED for curvilinear
       if (is_curv3) then
          xc3(1-ngh_:nx+ngh_,1-ngh_:ny+ngh_,1-ngh_:nz+ngh_)=dummy(1:nex,1:ney,1:nez,1)
          yc3(1-ngh_:nx+ngh_,1-ngh_:ny+ngh_,1-ngh_:nz+ngh_)=dummy(1:nex,1:ney,1:nez,2)
          zc3(1-ngh_:nx+ngh_,1-ngh_:ny+ngh_,1-ngh_:nz+ngh_)=dummy(1:nex,1:ney,1:nez,3)
       else
          if (is_curv) then
             xc(1-ngh_:nx+ngh_,1-ngh_:ny+ngh_)=dummy(1:nex,1:ney,1,1)
             yc(1-ngh_:nx+ngh_,1-ngh_:ny+ngh_)=dummy(1:nex,1:ney,1,2)
          else
             x(1-ngh_:nx+ngh_)=dummy(1:nex,1,1,1)
             y(1-ngh_:ny+ngh_)=dummy(1,1:ney,1,2)
          endif
          z(1-ngh_:nz+ngh_)=dummy(1,1,1:nez,3)
       endif

       ! free temporary data
       deallocate(dummy)

    endif ! end read grid if tecplot format

    ! free pointer varlist & dataname
    ! ===============================
    deallocate(varlist,dataname)

  end subroutine read_write_grid3d

  !==============================================================================================
  subroutine free_grid_ex3d
  !==============================================================================================
    !> author: XG
    !> date: October 2023
    !> Free MPI-IO types for grid
  !==============================================================================================
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! -------------------------------------------------------------------------------------------

    ! Free types for MPI-IO of grid
    ! -----------------------------
    call MPI_TYPE_FREE(grid_ex_(nob(iproc))%type_mat_nogh,info)
    call MPI_TYPE_FREE(grid_ex_(nob(iproc))%type_mat_view,info)

    ! Deallocate grid structure
    ! --------------------------
    deallocate(grid_ex_)

  end subroutine free_grid_ex3d

  !==============================================================================================
  subroutine mpicheckfile(filen)
  !==============================================================================================
    !> author: Luca Sciacovelli
    !> date: April 2018
    !> MPI subroutine to check existence of a file
  !==============================================================================================
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Inputs/Outputs
    character(len=*) :: filen
    ! -------------------------------------------------------------------------------------------
    ! Local variables
    logical :: filex
    integer :: myid_local
    ! -------------------------------------------------------------------------------------------

    ! Determine rank in general communicator
    ! --------------------------------------
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid_local,info)

    ! Check existence of file
    ! -----------------------
    inquire(file=trim(filen),exist=filex)

    if (.not.filex) then
       call mpistop('Cannot find '//trim(filen)//' Shutting down...', 0)
    endif

  end subroutine mpicheckfile

end module mod_io
