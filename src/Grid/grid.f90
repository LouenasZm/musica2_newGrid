!===============================================================================
subroutine grid_define
!===============================================================================
  !> author: AB & XG
  !> date: February 2024
  !> Definition of grid:
  !> * if FROM_FIELD (idepart=2) and same ghost points as before, reading
  !>   extended grid directly
  !> * if FROM_SCRATCH (idepart=1) or FROM_INTERP (idepart=3) or number of ghost
  !>   points has changed compared to the previous run, 2 possibilities:
  !>   -> if not user_defined (if user prescribed), grid directly read (file .x)
  !>   -> if user_defined, needs param_grid.ini file. Only relatively simple
  !>      geometry grid are handled.
!===============================================================================
  use mod_grid
  use mod_constant
  use mod_io
  use mod_mpi
  use warnstop
  use mod_utils
  implicit none
  ! ---------------------------------------------------------------------------
  logical  :: iexist
  integer :: i,j,k,mb
  integer :: ngx_,ngy_,ngz_
  integer :: ibl2read
  real(wp) :: z_mid,dz0
  character(len=200) :: gridname
  ! ---------------------------------------------------------------------------

  ! ------------------------
  ! Old way of defining grid <~ TO BE SUPPRESSED
  ! ------------------------
  if (is_grid_old) then
     call grid_define_old
     return
  endif

  ! ---------------
  ! Grid definition <~ New version
  ! ---------------
  ! Only the extended grid is considered now
  is_read_ex=.true.  ! obsolete ? Kept for compatibility

  ! Block number
  ! ------------
  ibl2read = nob(iproc)
  if (idepart==POST_PROCESSING) ibl2read = iblc_pp

  ! Allocate 3D curvilinear grid
  ! ----------------------------
  if (is_curv3) then
     ! Allocation of local grid
     allocate(xc3(nx1:nx2,ny1:ny2,nz1:nz2))
     allocate(yc3(nx1:nx2,ny1:ny2,nz1:nz2))
     allocate(zc3(nx1:nx2,ny1:ny2,nz1:nz2))
     xc3=0.0_wp; yc3=0.0_wp; zc3=0.0_wp

  ! Allocate 2D curvilinear grid
  ! ----------------------------
  else if (is_curv) then
     ! Allocation of local grid
     allocate(xc(nx1:nx2,ny1:ny2),yc(nx1:nx2,ny1:ny2),z(nz1:nz2))
     xc=0.0_wp ; yc=0.0_wp; z=0.0_wp
     ! for inviscid fluxes TO BE CHANGED why extended in ghost points? for IRS????
     allocate(idz(1-ngh:nz+ngh))
     ! for viscous fluxes (extended to ghost points for double derivatives)
     allocate(idz_v(1-ngh_v:nz+ngh_v))

  ! Allocate Cartesian grid
  ! -----------------------
  else
     ! Allocation of local grid
     allocate(x(nx1:nx2),y(ny1:ny2),z(nz1:nz2))
     x=0.0_wp; y=0.0_wp; z=0.0_wp
     ! for inviscid fluxes TO BE CHANGED why extended in ghost points? for IRS????
     allocate(idx(1-ngh:nx+ngh),idy(1-ngh:ny+ngh),idz(1-ngh:nz+ngh))
     ! for viscous fluxes (extended to ghost points for double derivatives)
     allocate(idx_v(1-ngh_v:nx+ngh_v),idy_v(1-ngh_v:ny+ngh_v),idz_v(1-ngh_v:nz+ngh_v))
  endif


  ! Check existence of file
  ! -----------------------
  ! Even if FROM_FIELD, may not exist if ngh changed w.r.t previous run
  if (is_curv3) then
     gridname='grid_bl'//trim(numchar(ibl2read))//'_ngh'//trim(numchar(ngh))//filext_read
  else
     gridname='grid_bl'//trim(numchar(ibl2read))//'_ngh'//trim(numchar(ngh))//'.bin'
  endif
  inquire(file=gridname,exist=iexist)

  ! Grid necessary for post-processing mode
  if ((idepart==POST_PROCESSING).and.(.not.iexist)) call mpistop("Grid "//gridname//" not found...",1)

  ! Generation of global grid without ghost point by leader processor if necessary
  ! ------------------------------------------------------------------------------
  if ((idepart.eq.FROM_SCRATCH).or.(.not.iexist).or.(is_init_2D3D)) then
     if (iproc.eq.iproc_leader(nob(iproc))) then
        ! Global grid of the block & Writting of grid in binary file in MPI_IO without ghost points
        ! -----------------------------------------------------------------------------------------
        ! If is_curv3, necessary to have local grid for init_periodicity & grid_directions_3d in mpi_connect
        if (is_curv3) then
           ! Allocate extended global grid
           allocate(xgc3e(1-ngh:ngx+ngh,1-ngh:ngy+ngh,1-ngh:ngz+ngh))
           allocate(ygc3e(1-ngh:ngx+ngh,1-ngh:ngy+ngh,1-ngh:ngz+ngh))
           allocate(zgc3e(1-ngh:ngx+ngh,1-ngh:ngy+ngh,1-ngh:ngz+ngh))

           ! Reading global grid in 3d
           ! -------------------------
           if (is_adjoint_block) then
              gridname=trim(dirGRID)//'/'//trim(nameGRID)//'_bl'
           else
              gridname=trim(dirGRID)//'/'//trim(nameGRID)//'_mod_bl'
           endif
           open(50,file=trim(gridname)//trim(numchar(nob(iproc)))//'.x',form='formatted')
           rewind(50)
           read(50,*) mb
           read(50,*) ngx_,ngy_,ngz_
           read(50,*) (((xgc3e(i,j,k),i=1,ngx_),j=1,ngy_),k=1,ngz_)
           read(50,*) (((ygc3e(i,j,k),i=1,ngx_),j=1,ngy_),k=1,ngz_)
           read(50,*) (((zgc3e(i,j,k),i=1,ngx_),j=1,ngy_),k=1,ngz_)
           close(50)

           ! Writting grid without ghost point
           ! ---------------------------------
           gridname='grid_bl'//trim(numchar(nob(iproc)))//'.bin'
           call read_write_grid_curv3(gridname,0,WRITE_LEADER)
        else
           ! Creation of spanwise grid direction
           ! -----------------------------------
           dz0=deltaz
           if (is_2D) then
              allocate(zg(1:1))
              zg=0.0_wp
           else
              allocate(zg(1-ngh:ngz+ngh))
              if (nstretchz.eq.0) then
                 zg(ngz/2)=0.0_wp
                 do k=ngz/2+1,ngz
                    zg(k)=zg(k-1)+dz0
                 enddo
                 do k=ngz/2-1,1,-1
                    zg(k)=zg(k+1)-dz0
                 enddo
              else
                 ! Initialization
                 zg(1)=0.0_wp
                 ! First grid points, if stretching doesn't begin straight at 1 (nrz1(1)>1)
                 if (nrz1(1)>1) then
                    do k=2,nrz1(1)
                       zg(k)=zg(k-1)+dz0
                    enddo
                 endif
                 ! Loop on number of stretching
                 do i=1,nstretchz
                    ! Stretching
                    do k=max(nrz1(i),2),nrz2(i)
                       zg(k)=zg(k-1)+dz0
                       dz0=rsz(i)*dz0
                    enddo
                    ! Uniform until next stretching or ngz
                    if (i.ne.nstretchz) then
                       do k=nrz2(i),nrz1(i+1)
                          zg(k)=zg(k-1)+dz0
                       enddo
                    else
                       do k=nrz2(i),ngz
                          zg(k)=zg(k-1)+dz0
                       enddo
                    endif
                 enddo
                 ! Grid centered on zg(ngz/2)
                 z_mid=zg(ngz/2)
                 do k=1,ngz
                    zg(k)=zg(k)-z_mid
                 enddo
              endif
           endif

           if (is_def_grid) then
              ! Generating global grid
              ! ----------------------
              call grid_generation
           else
              ! Reading global grid in 2d
              ! -------------------------
              allocate(xgce(1-ngh:ngx+ngh,1-ngh:ngy+ngh),ygce(1-ngh:ngx+ngh,1-ngh:ngy+ngh))
              if (is_adjoint_block) then
                 gridname=trim(dirGRID)//'/'//trim(nameGRID)//'_bl'
              else
                 gridname=trim(dirGRID)//'/'//trim(nameGRID)//'_mod_bl'
              endif
              open(50,file=trim(gridname)//trim(numchar(nob(iproc)))//'.x',form='formatted')
              rewind(50)
              read(50,*) mb
              read(50,*) ngx_,ngy_,ngz_
              read(50,*) ((xgce(i,j),i=1,ngx_),j=1,ngy_)
              read(50,*) ((ygce(i,j),i=1,ngx_),j=1,ngy_)
              close(50)

              if (.not.is_curv) then
                 allocate(xg(1-ngh:ngx+ngh),yg(1-ngh:ngy+ngh))
                 xg(1-ngh:ngx+ngh) = xgce(1-ngh:ngx+ngh,1)
                 yg(1-ngh:ngy+ngh) = ygce(1,1-ngh:ngy+ngh)
              endif
           endif

           ! Writting grid without ghost point
           ! ---------------------------------
           gridname='grid_bl'//trim(numchar(nob(iproc)))//'.bin'
           if (is_curv) then
              call read_write_grid_curv(gridname,0,WRITE_LEADER)
           else
              call read_write_grid_cart(gridname,0,WRITE_LEADER)
           endif

        endif

     ! if not proc_leader
     else
        ! Create temporary communicator
        call MPI_COMM_SPLIT(COMM_global,iproc,iproc,COMM_rw_grid,info)
     endif
  endif

  ! Waiting for everybody
  call MPI_BARRIER(COMM_global,info)

  ! Reading of local grid by every processor
  ! ----------------------------------------
  if ((iexist).and.(idepart.ne.FROM_SCRATCH).and.(.not.is_init_2D3D)) then
     ! If file with extended points already exist, directly read here
     ! --------------------------------------------------------------
     gridname='grid_bl'//trim(numchar(ibl2read))//'_ngh'//trim(numchar(ngh))//'.bin'
     ! init MPI-IO for grid write global grid & read extended grid & free MPI IO memory
     if (is_curv3) then
        call read_write_grid_curv3(gridname,ngh,READ)
     else if (is_curv) then
        call read_write_grid_curv(gridname,ngh,READ)
     else
        call read_write_grid_cart(gridname,ngh,READ)
     endif
  else
     ! Else, local grid without ghost point is read, to allow mpi_connect
     ! and call to add_ghost_cells_grid performed after call to mpi_connect
     ! --------------------------------------------------------------------
     gridname='grid_bl'//trim(numchar(nob(iproc)))//'.bin'
     if (is_curv3) then
        call read_write_grid_curv3(gridname,0,READ)
     else if (is_curv) then
        call read_write_grid_curv(gridname,0,READ)
     else
        call read_write_grid_cart(gridname,0,READ)
     endif
  endif

  ! Waiting for everybody
  call MPI_BARRIER(COMM_global,info)

end subroutine grid_define


!===============================================================================
subroutine grid_extend
!===============================================================================
  !> author: AB & XG
  !> date: February 2024
  !> Finalization of grid:
  !> * if FROM_SCRATCH (idepart=1) or FROM_INTERP (idepart=3) or number of ghost
  !>   points has changed compared to the previous run:
  !>   -> extension of grid with ghost point by leader processor (needs "call
  !>      mpi_connect" to have been done)
  !>   -> writting grid_bl*_ngh*.bin file by leader processor
  !>   -> every proc read extended grid
  !> * Scaling local grid by Lgrid
!===============================================================================
  use mod_constant
  use mod_grid
  use mod_mpi
  use mod_io
  use mod_utils
  use warnstop
  use mod_add_gh
  use mod_add_gh3d
  implicit none
  ! ---------------------------------------------------------------------------
  logical  :: iexist
  integer :: ibl2read
  character(len=200) :: gridname
  ! ---------------------------------------------------------------------------

  ! Block number
  ! ------------
  ibl2read = nob(iproc)
  if (idepart==POST_PROCESSING) ibl2read = iblc_pp

  ! Check existence of file with extended ghost points
  ! --------------------------------------------------
  ! Even if FROM_FIELD, may not exist if ngh changed w.r.t previous run
  if (is_curv3) then
     gridname='grid_bl'//trim(numchar(ibl2read))//'_ngh'//trim(numchar(ngh))//filext_read
  else
     gridname='grid_bl'//trim(numchar(ibl2read))//'_ngh'//trim(numchar(ngh))//'.bin'
  endif
  inquire(file=gridname,exist=iexist)

  ! Generation of global grid with ghost point by leader processor if necessary
  ! ---------------------------------------------------------------------------
  if ((idepart.eq.FROM_SCRATCH).or.(.not.iexist).or.(is_init_2D3D)) then
    ! Temporary communicator
    ! ----------------------
     call MPI_COMM_SPLIT(COMM_global,iproc,iproc,COMM_r_grid,info)

     ! if  proc_leader
     if (iproc.eq.iproc_leader(nob(iproc))) then
        ! Extension & Writting of grid in binary file in MPI_IO with ghost points
        ! -----------------------------------------------------------------------
        if (is_curv3) then
           call add_ghost_cells_grid3d
        else
           call add_ghost_cells_grid
        endif
        ! Waiting for everybody before rm grid_blx.bin file
        call MPI_BARRIER(COMM_global,info)
        ! Remove temporary grid_bl*.bin file generated in grid_define
        gridname='grid_bl'//trim(numchar(nob(iproc)))//'.bin'
        if (iproc.eq.iproc_leader(nob(iproc))) call system('rm '//trim(gridname))

     ! if not proc_leader
     else
        ! Create temporary communicator for WRITE_LEADER grid with ngh
        call MPI_COMM_SPLIT(COMM_global,iproc,iproc,COMM_rw_grid,info)
        ! Waiting for everybody before rm grid_blx.bin file
        call MPI_BARRIER(COMM_global,info)
     endif

     ! Reading of extended local grid by every processor
     ! -------------------------------------------------
     gridname='grid_bl'//trim(numchar(nob(iproc)))//'_ngh'//trim(numchar(ngh))//'.bin'
     ! init MPI-IO for grid write global grid & read extended grid & free MPI IO memory
     if (is_curv3) then
        call read_write_grid_curv3(gridname,ngh,READ)
     else if (is_curv) then
        call read_write_grid_curv(gridname,ngh,READ)
     else
        call read_write_grid_cart(gridname,ngh,READ)
     endif
  endif

  ! Scale dimension for grid
  ! ------------------------
  ! Lgrid checked in setupref
  if (is_curv3) then
     xc3=xc3*Lgrid
     yc3=yc3*Lgrid
     zc3=zc3*Lgrid
  else
     if (is_curv) then
        xc=xc*Lgrid
        yc=yc*Lgrid
     else
        x=x*Lgrid
        y=y*Lgrid
     endif
     z=z*Lgrid
  endif

  ! Computation of extrapolation coefficients for borders_wall
  ! ----------------------------------------------------------
  ! Done in mod_coeff_deriv by init_coeff_wall_cart, called in musica_main.f90

end subroutine grid_extend



!===============================================================================
subroutine grid_generation
!===============================================================================
  !> author: AB
  !> date: March 2024
  !> Generation in the solver of basic grid
  !> 1: Cartesian, with stretching
  !> 2: polar grid for cylinder [to be completed...]
!===============================================================================
  use mod_constant
  use mod_grid
  use mod_mpi
  use warnstop
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j,k
  real(wp) :: z_mid,dz0
  ! ---------------------------------------------------------------------------

  ! Initialization <~ Obsolete ?
  deltax=deltax0(nob(iproc)); deltay=deltay0(nob(iproc))

  ! Grid in direction i & j
  ! -----------------------
  select case(num_grid)

  case(1) ! 1: Cartesian grids
     call grid_cartesian
  case(2) ! 2: Polar grid for cylinder
     call grid_polar_cyl
  case default
     call mpistop('Bad choice of type of grid in param_grid.ini', 1)
  end select

end subroutine grid_generation




!===============================================================================
subroutine grid_cartesian
!===============================================================================
  !> author: AB
  !> date: March 2024
  !> Generation of cartesian grid, with stretching
!===============================================================================
  use mod_block
  use mod_constant
  use mod_grid
  use mod_mpi
  use warnstop
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j,k,nbl
  real(wp) :: dx0,dy0
  ! ---------------------------------------------------------------------------
  integer :: ngx_,ngy_,ngz_

  if ((is_curv).or.(is_curv3)) &
    call mpistop("Wrong type of grid & solver combination, solver must be cartesian for cartesian grid generation",1)

  ! Block number
  nbl=nob(iproc)

  ! Initialization
  ! --------------
  dx0=deltax0(nbl); dy0=deltay0(nbl)

  ! Generation of xg
  ! ----------------
  allocate(xg(1-ngh:ngx+ngh))
  if (nstretchx(nbl).eq.0) then
     xg(1)=x0(nbl)
     do i=2,ngx
        xg(i)=xg(i-1)+dx0
     enddo
  else
     ! Initialization
     xg(1)=x0(nbl)
     ! First grid points, if stretching doesn't begin straight at 1 (nrx1(1)>1)
     if (bl(nbl)%nrx1(1)>1) then
        do i=2,bl(nbl)%nrx1(1)
           xg(i)=xg(i-1)+dx0
        enddo
     endif
     ! Loop on number of stretching
     do k=1,nstretchx(nbl)
        ! Stretching
        do i=max(bl(nbl)%nrx1(k),2),bl(nbl)%nrx2(k)
           xg(i)=xg(i-1)+dx0
           dx0=bl(nbl)%rsx(k)*dx0
        enddo
        ! Uniform until next stretching or ngx
        if (k.ne.nstretchx(nbl)) then
           do i=bl(nbl)%nrx2(k),bl(nbl)%nrx1(k+1)
              xg(i)=xg(i-1)+dx0
           enddo
        else
           do i=bl(nbl)%nrx2(k),ngx
              xg(i)=xg(i-1)+dx0
           enddo
        endif
     enddo
  endif

  ! Generation of yg
  ! ----------------
  allocate(yg(1-ngh:ngy+ngh))
  if (nstretchy(nbl).eq.0) then
     yg(1)=y0(nbl)
     do i=2,ngy
        yg(i)=yg(i-1)+dy0
     enddo
  else
     ! Initialization
     yg(1)=y0(nbl)
     ! First grid points, if stretching doesn't begin straight at 1 (nry1(1)>1)
     if (bl(nbl)%nry1(1)>1) then
        do i=2,bl(nbl)%nry1(1)
           yg(i)=yg(i-1)+dy0
        enddo
     endif
     ! Loop on number of stretching
     do k=1,nstretchy(nbl)
        ! Stretching
        do i=max(bl(nbl)%nry1(k),2),bl(nbl)%nry2(k)
           yg(i)=yg(i-1)+dy0
           dy0=bl(nbl)%rsy(k)*dy0
        enddo
        ! Uniform until next stretching or ngy
        if (k.ne.nstretchy(nbl)) then
           do i=bl(nbl)%nry2(k),bl(nbl)%nry1(k+1)
              yg(i)=yg(i-1)+dy0
           enddo
        else
           do i=bl(nbl)%nry2(k),ngy
              yg(i)=yg(i-1)+dy0
           enddo
        endif
     enddo
  endif

  open(51,file=trim(dirDATA)//'vites_u.bin',form='unformatted',status='unknown')
  rewind(51)
  ! global dim
  ! ==========
  read(51) ngx_
  read(51) ngy_
  read(51) ngz_
  ! global grid
  ! ===========
  read(51) (xg(i),i=1,ngx)
  read(51) (yg(j),j=1,ngy)
  read(51) (zg(k),k=1,ngz)
  close(51)
  deltax=xg(2)-xg(1)
  deltaz=zg(2)-zg(1)

  ! Allocation
  allocate(xgce(1-ngh:ngx+ngh,1-ngh:ngy+ngh))
  allocate(ygce(1-ngh:ngx+ngh,1-ngh:ngy+ngh))
  xgce(1-ngh:ngx+ngh,1)=xg(1-ngh:ngx+ngh)
  ygce(1,1-ngh:ngy+ngh)=yg(1-ngh:ngy+ngh)

end subroutine grid_cartesian



!===============================================================================
subroutine grid_polar_cyl
!===============================================================================
  !> author: AB (adapted from grid_cyl)
  !> date: March 2024
  !> Generation of circular grid around cylinder
!===============================================================================
  use mod_block
  use mod_constant
  use mod_grid
  use mod_mpi_part
  use warnstop
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j,k,nbl
  real(wp) :: radius,dr,dteta
  real(wp), dimension(:), allocatable :: r ! radial grid
  ! ---------------------------------------------------------------------------

  ! Block number
  nbl=nob(iproc)

  if (nbloc.gt.1) call mpistop("grid_polar_cyl only working in monoblock for the moment",1)

  ! Cylinder radius =diameter(Lref)/2
  ! ---------------------------------
  radius=1/2.0_wp

  ! Allocation
  allocate(xgce(1-ngh:ngx+ngh,1-ngh:ngy+ngh))
  allocate(ygce(1-ngh:ngx+ngh,1-ngh:ngy+ngh))


  ! Radial distribution
  ! -------------------
  ! Stretching rate in radial direction
  ! <~ specified in param_grid in direction normal to wall
  if (bl(nbl)%BC(1).eq.0) then       ! if wall at imin
     allocate(r(ngx))
     ! Initialization
     r(1)=radius; dr=deltax0(nbl)
     ! If no stretching
     if (nstretchx(nbl).eq.0) then
        do i=2,ngx
           r(i)=r(i-1)+dr
        enddo
     ! If stretching
     else
        ! First grid points, if stretching doesn't begin straight at 1 (nrx1(1)>1)
        if (bl(nbl)%nrx1(1)>1) then
           do i=2,bl(nbl)%nrx1(1)
              r(i)=r(i-1)+dr
           enddo
        endif
        ! Loop on number of stretching
        do k=1,nstretchx(nbl)
           ! Stretching
           do i=max(bl(nbl)%nrx1(k),2),bl(nbl)%nrx2(k)
              r(i)=r(i-1)+dr
              dr=bl(nbl)%rsx(k)*dr
           enddo
           ! Uniform until next stretching or ngx
           if (k.ne.nstretchx(nbl)) then
              do i=bl(nbl)%nrx2(k),bl(nbl)%nrx1(k+1)
                 r(i)=r(i-1)+dr
              enddo
           else
              do i=bl(nbl)%nrx2(k),ngx
                 r(i)=r(i-1)+dr
              enddo
           endif
        enddo
     endif
  else if (bl(nbl)%BC(3).eq.0) then  ! if wall at jmin
     allocate(r(ngy))
     ! Initialization
     r(1)=radius; dr=deltay0(nbl)
     ! If no stretching
     if (nstretchy(nbl).eq.0) then
        do i=2,ngy
           r(i)=r(i-1)+dr
        enddo
     ! If stretching
     else
        ! First grid points, if stretching doesn't begin straight at 1 (nry1(1)>1)
        if (bl(nbl)%nry1(1)>1) then
           do i=2,bl(nbl)%nry1(1)
              r(i)=r(i-1)+dr
           enddo
        endif
        ! Loop on number of stretching
        do k=1,nstretchy(nbl)
           ! Stretching
           do i=max(bl(nbl)%nry1(k),2),bl(nbl)%nry2(k)
              r(i)=r(i-1)+dr
              dr=bl(nbl)%rsy(k)*dr
           enddo
           ! Uniform until next stretching or ngy
           if (k.ne.nstretchy(nbl)) then
              do i=bl(nbl)%nry2(k),bl(nbl)%nry1(k+1)
                 r(i)=r(i-1)+dr
              enddo
           else
              do i=bl(nbl)%nry2(k),ngy
                 r(i)=r(i-1)+dr
              enddo
           endif
        enddo
     endif
  endif

  ! Azimuthal step
  ! --------------
  if (is_adjoint_block) then
     dteta=2.*pi/dble(ngx-1)
  else
     dteta=2.*pi/dble(ngx)
  endif

  ! Polar grid
  ! ----------
  do i=1,ngx
     do j=1,ngy
        xgce(i,j)= r(j)*cos((i-1)*dteta+2*pi/2) + x0(nbl)
        ygce(i,j)= r(j)*sin((i-1)*dteta+2*pi/2) + y0(nbl)
     enddo
  enddo

end subroutine grid_polar_cyl
