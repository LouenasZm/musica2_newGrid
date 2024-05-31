!===============================================================================
subroutine grid_define_old
!===============================================================================
  !> Old way of creating grid
  !> Grid generation directly written in subroutines
!===============================================================================
  use mod_grid
  use mod_constant
  use mod_mpi
  use mod_io
  use warnstop
  use mod_utils
  use mod_pp_var
  use mod_bc_periodicity ! <~ needed for Lxp,Lyp,Lzp  TO BE CHANGED
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j,k
  integer :: ngx_,ngy_,ngz_,ibl2read
  real(wp) :: chord
  character(len=120) :: gridname
  ! ---------------------------------------------------------------------------

  ! allocate 3D curvilinear grid
  ! ============================
  if (is_curv3) then
     ! global grid
     if ((.not.is_read_ex).or.(idepart.eq.1)) then
        allocate(xgc3(1-ngh:ngx+ngh,1-ngh:ngy+ngh,1-ngh:ngz+ngh))
        allocate(ygc3(1-ngh:ngx+ngh,1-ngh:ngy+ngh,1-ngh:ngz+ngh))
        allocate(zgc3(1-ngh:ngx+ngh,1-ngh:ngy+ngh,1-ngh:ngz+ngh))
        xgc3=0.0_wp
        ygc3=0.0_wp
        zgc3=0.0_wp
     endif
     ! local grid
     allocate(xc3(nx1:nx2,ny1:ny2,nz1:nz2))
     allocate(yc3(nx1:nx2,ny1:ny2,nz1:nz2))
     allocate(zc3(nx1:nx2,ny1:ny2,nz1:nz2))
     xc3=0.0_wp; yc3=0.0_wp; zc3=0.0_wp

  ! allocate 2D curvilinear grid
  ! ============================
  else if (is_curv) then
     allocate(xgc(1-ngh:ngx+ngh,1-ngh:ngy+ngh),xc(nx1:nx2,ny1:ny2))
     allocate(ygc(1-ngh:ngx+ngh,1-ngh:ngy+ngh),yc(nx1:nx2,ny1:ny2))
     allocate(z(nz1:nz2))
     xgc=0.0_wp; xc=0.0_wp
     ygc=0.0_wp; yc=0.0_wp
     if (is_2D) then
        allocate(zg(1:1))
     else
        allocate(zg(1-ngh:ngz+ngh))
     endif
     zg=0.0_wp

  ! allocate Cartesian grid
  ! =======================
  else
     allocate(xg(1-ngh:ngx+ngh),yg(1-ngh:ngy+ngh))
     allocate(x(nx1:nx2),y(ny1:ny2),z(nz1:nz2))
     xg=0.0_wp; x=0.0_wp
     yg=0.0_wp; y=0.0_wp
     if (is_2D) then
        allocate(zg(1:1))
     else
        allocate(zg(1-ngh:ngz+ngh))
     endif
     z=0.0_wp; zg=0.0_wp
  endif

  ! allocate Cartesian metrics
  ! ==========================
  ! for inviscid fluxes TO BE CHANGED why extended in ghost points? for IRS????
  allocate(idx(1-ngh:nx+ngh),idy(1-ngh:ny+ngh),idz(1-ngh:nz+ngh))
  ! for viscous fluxes (extended to ghost points for double derivatives)
  allocate(idx_v(1-ngh_v:nx+ngh_v),idy_v(1-ngh_v:ny+ngh_v),idz_v(1-ngh_v:nz+ngh_v))

  ! init igrd
  igrd=1

  if (idepart==POST_PROCESSING) then
     ! Read global grid per block (fortran IEEE binary file - Att big or little endian)
     ! ==========================
     ibl2read=nblr(nob(iproc))
     if (is_curv3) then
        if ((is_read_ex).and.(idepart.ne.1)) then
           ! init MPI-IO for grid read/write tata
           call init_io_grid_ex3d
           ! gridname
           gridname=trim(dirGRID)//'/'//trim(nameGRID)// &
                '_ngh'//trim(numchar(ngh))//'_mod_bl'//trim(numchar(nob(iproc)))//filext_read
           ! read extended grid
           call read_write_grid3d(trim(gridname),READ)
           ! free memory
           call free_grid_ex3d
        else
           call grid_read_ex(trim(dirGRID)//'/'//trim(nameGRID)// &
                '_ngh'//trim(numchar(ngh))//'_mod_bl'//trim(numchar(nob(iproc)))//'.x')
        endif
        ! bypass of comm grid
        if (is_curv3) igrd=6
     else
        if (iproc==iproc_leader(nob(iproc))) print *,'read '//trim(dirDATA)//'grid_bl'//trim(numchar(nob(iproc)))//'.bin'
        open(194,file=trim(dirDATA)//'grid_bl'//trim(numchar(ibl2read))//'.bin', &
             form='unformatted',status='unknown')
        rewind(194)
        read(194) ngx_
        read(194) ngy_
        read(194) ngz_
        if (is_curv3) then
           read(194) (((xgc3(i,j,k),i=1,ngx),j=1,ngy),k=1,ngz)
           read(194) (((ygc3(i,j,k),i=1,ngx),j=1,ngy),k=1,ngz)
           read(194) (((zgc3(i,j,k),i=1,ngx),j=1,ngy),k=1,ngz)
        else
           if (is_curv) then
              read(194) ((xgc(i,j),i=1,ngx),j=1,ngy)
              read(194) ((ygc(i,j),i=1,ngx),j=1,ngy)
           else
              read(194) (xg(i),i=1,ngx)
              read(194) (yg(j),j=1,ngy)
           endif
           read(194) (zg(k),k=1,ngz)
        endif
        close(194)
     endif

     if ((TURB).and.(nbloc==9)) then
        chord=L_ref
        Lxp=0.0_wp*chord
        Lyp=0.709666666666667_wp*chord
        if (is_curv3) then
           if ((is_read_ex).and.(idepart.ne.1)) then
              xc3=xc3*chord
              yc3=yc3*chord
              zc3=zc3*chord
           else
              xgc3=xgc3*chord
              ygc3=ygc3*chord
              zgc3=zgc3*chord
           endif
        endif
     endif

     if (iproc==0) print *,'Global grid OK'

     return
  endif

  ! -----------
  ! define grid <~ old way
  ! -----------
  if ((TGV).or.(CHIT)) then
     call grid_periodic
  elseif (CHAN) then
     call grid_chan
  elseif (PHILL) then
     call grid_phill
  elseif (STBL) then
     call grid_stbl
     !call grid_flatplate
     !call grid_backstep
  elseif (CYL) then
     call grid_cyl
  elseif ((SRC).or.(CAV)) then ! fourre-tout
     call grid_src
  elseif (ACT) then
     call grid_act
  elseif (SHIT) then
     call grid_shit
  elseif (TURB) then
     call grid_turb
  elseif (LE) then
     call grid_le
  elseif (TE) then
     call grid_TE
  endif

end subroutine grid_define_old

!===============================================================================
subroutine grid_periodic
!===============================================================================
  !> predefined grid for THI, TGV cases ...
!===============================================================================
  use mod_mpi
  use mod_flow
  use mod_constant
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j,k
  real(wp) :: Lx,Ly,Lz
  ! ---------------------------------------------------------------------------

  ! Construction of global grid
  ! ---------------------------  
!!$  Lx = 2.0_wp*pi;   deltax= Lx/dble(ngx-1)
!!$  Ly = 2.0_wp*pi;   deltay= Ly/dble(ngy-1)
!!$  Lz = 2.0_wp*pi;   deltaz= Lz/dble(ngz-1)
  Lx=2.0_wp*pi*L_ref;   deltax=Lx/dble(ngx)
  Ly=2.0_wp*pi*L_ref;   deltay=Ly/dble(ngy)
  Lz=2.0_wp*pi*L_ref;   deltaz=Lz/dble(ngz)
 
  do i=1,ngx
     xg(i)=dble(i-1)*deltax
  enddo
  do j=1,ngy
     yg(j)=dble(j-1)*deltay
  enddo
  do k=1,ngz
     zg(k)=dble(k-1)*deltaz
  enddo

  ! curvilinear grid (for test)
  if (is_curv) then
     do i=1,ngx
        do j=1,ngy
           xgc(i,j)=xg(i)
           ygc(i,j)=yg(j)
        enddo
     enddo
  endif

end subroutine grid_periodic

!===============================================================================
subroutine grid_src
!===============================================================================
  !> predefined grids for src, cav, actuator, ...
!===============================================================================
  use mod_flow
  use mod_constant
  use mod_mpi
  use mod_block
  use mod_utils
  use mod_mpi_part
  use mod_bc_periodicity ! <~ needed for Lxp,Lyp,Lzp  TO BE CHANGED
  use warnstop
  implicit none
  ! ---------------------------------------------------------------------------
  integer  :: i,j,k
  real(wp) :: dx0,dy0,dz0,ry
  ! real(wp) :: x0,alp
  real(wp) :: Length,Depth,LsD
  integer  :: hxm,hym,it
  real(wp) :: rx1,rx3,rx4,rx5,ry1,ry2,ry3
  real(wp) :: rxa,rxb,rxc,rya,ryb,ryc,epsi,erreur
  real(wp) :: Lx,Ly,Lz
  ! ---------------------------------------------------------------------------
  ! sinusoidal grid
  real(wp) :: ax,ay,az,enex,eney
  real(wp) :: enxy,enxz,enyx,enyz,enzx,enzy
  ! ---------------------------------------------------------------------------

  ! Construction of global grid
  ! ---------------------------
  dx0 = deltax
  dy0 = deltay
  dz0 = deltaz

  ry = 1.0_wp

  ! Some predefined grids TO BE CHANGED (preliminary XG)
  ! =====================
  ! igrd=1: 1 block (cart or curv)
  ! igrd=2: 2 blocks (cart or curv) if non zero alp -> curv (2nd block is skewed)
  ! igrd=3: 4 blocks acoustics scatt in cavity (cart but can be run curv)
  ! igrd=4: 4 blocks low-Re cavity flow (cart but can be run in curv)
  ! igrd=5: actuator (14 blocks) (necessarily curv)
  ! igrd=6: 1 block curv: skewed or sinusoidal grid
  ! igrd=7: 4 blocks curv: test 4 blocks skewed
  ! igrd=8: 5 blocks curv: test 5 blocks multiple point
  ! igrd=9: 1 block curv3: sinusoidal grid
  ! igrd=10: 5 blocks curv3: 3D nozzle
  ! igrd=11: 3 blocks test IRS
  if ((SRC).and.(nbloc==1)) igrd=1
  if ((SRC).and.(nbloc==2)) igrd=2
!!$  if ((SRC).and.(nbloc==1).and.(is_curv)) igrd=6 ! (sinusoidal grid)
  if ((SRC).and.(nbloc==4)) igrd=7
  if ((SRC).and.(nbloc==5)) igrd=8
  if ((SRC).and.(nbloc==1).and.(is_curv3)) igrd=9
  if ((SRC).and.(nbloc==5).and.(is_curv3)) igrd=10
  if ((SRC).and.(nbloc==3)) igrd=11

  if (CAV) then
     if (iorder_visc==0) then
        igrd=3
     else
        igrd=4
     endif
  endif
  if (ACT) igrd=5
  ! for the source case, just read the grid file.
  if (is_curv3) then
     !if (SRC) igrd=9
     if (SRC) igrd=10
     !if (SRC) igrd=2
     !if (SRC) igrd=12
     !if (SRC) igrd=13
 endif
  if ((igrd==1).and.(nbloc.ne. 1)) call mpistop('grid#1 for nbloc=1 !', 0)
  if ((igrd==2).and.(nbloc.ne. 2)) call mpistop('grid#2 for nbloc=2 !', 0)
  if ((igrd==3).and.(nbloc.ne. 4)) call mpistop('grid#3 for nbloc=4 !', 0)
  if ((igrd==4).and.(nbloc.ne. 4)) call mpistop('grid#4 for nbloc=4 !', 0)
  !if ((igrd==5).and.(nbloc.ne.14)) call mpistop('grid#5 for nbloc=14 !', 0)
  if ((igrd==6).and.(nbloc.ne. 1)) call mpistop('grid#6 for nbloc=1 !', 0)
  if ((igrd==7).and.(nbloc.ne. 4)) call mpistop('grid#7 for nbloc=4 !', 0)
  if ((igrd==8).and.(nbloc.ne. 5)) call mpistop('grid#8 for nbloc=5 !', 0)

  !!igrd=1
  
  select case (igrd)
  !*********************************************************************   
  case (1) ! 
  !*********************************************************************
     ! along x
     xg(ngx/2) = 0.0_wp
     do i=ngx/2+1,ngx
        xg(i)=xg(i-1)+dx0
     enddo
     do i=ngx/2-1,1,-1
        xg(i)=xg(i+1)-dx0
     enddo
     ! along y
     yg(ngy/2) = 0.0_wp
     do j=ngy/2+1,ngy
        yg(j)=yg(j-1)+dy0
     enddo
     do j=ngy/2-1,1,-1
        yg(j)=yg(j+1)-dy0
     enddo
     
!!$     yg(1) = 0.0_wp
!!$     do j=2,ngy
!!$        yg(j)=yg(j-1)+dy0
!!$     enddo
     
!!$     ! along x
!!$     xg(1)=0.0_wp
!!$     do i=2,ngx
!!$        xg(i)=xg(i-1)+dx0
!!$     enddo
!!$     ! along y
!!$     yg(1)=0.0_wp
!!$     do j=2,ngy
!!$        yg(j)=yg(j-1)+dy0
!!$        dy0=ry*dy0
!!$     enddo
     
!!$     ! test-case shock-vortex interaction
!!$     ! ----------------------------------
!!$     ! along x
!!$     xg(1)=0.0_wp
!!$     do i=2,ngx
!!$        xg(i)=xg(i-1)+dx0
!!$     enddo
!!$     ! along y
!!$     yg(1)=0.0_wp
!!$     do j=2,ngy
!!$        yg(j)=yg(j-1)+dy0
!!$        dy0=ry*dy0
!!$     enddo

     !deltax=60.e-3
     !xg=xg*deltax
     !yg=yg*deltax
     
     ! Tam & Dong radiation center (default) TO BE CHANGED
     xcr_=xg(ngx/2)
     ycr_=yg(ngy/2)    
!!$     ycr_=yg(1)    
     ! convert to curvilinear grid
     if (is_curv) then
        do i=1,ngx
           do j=1,ngy
              xgc(i,j)=xg(i)
              ygc(i,j)=yg(j)
           enddo
        enddo
     endif

!!$     print *,'xg'
!!$     print *,xg
!!$     print *,'yg'
!!$     print *,yg
!!$     stop

     is_curv=.true.
     !call grid_read('grid_bl'//trim(numchar(nob(iproc)))//'.x')
     !call grid_read('bosse2_imax.x')
     !call grid_read('bosse2.x')
     !call grid_read('conv2.x')
     !call grid_read('conv2_imax.x')
     !call grid_read('conv4_imin.x')
     if (is_adjoint_block) then
        call grid_read('conv5_bl'//trim(numchar(nob(iproc)))//'.x')
     else
        !call grid_read('conv5_mod_bl'//trim(numchar(nob(iproc)))//'.x')
        call grid_read('Grid_LES_GO/ls59c_mod_bl'//trim(numchar(nob(iproc)))//'.x')
     endif
     !call grid_read('cdnoz2Df_jmin.x')
     !call grid_read('nozSU2.x')
     !call grid_read('nozSU2_jmin.x')
     deltax=1.e-6
     deltax=60.e-3

!!$     ! along x
!!$     xg(1)=xgc(1,1)
!!$     dx0=xgc(2,1)-xgc(1,1)
!!$     do i=2,ngx
!!$        xg(i)=xg(i-1)+dx0
!!$     enddo
!!$     ! along y
!!$     do j=1,ngy
!!$        yg(j)=ygc(1,j)
!!$     enddo

!!$     ! 2D
!!$     do i=1,ngx
!!$        xg(i)=xgc(i,1)
!!$     enddo
!!$     do j=1,ngy
!!$        yg(j)=ygc(1,j)
!!$     enddo

!!$     ! convert to curvilinear grid
!!$     if (is_curv) then
!!$        do i=1,ngx
!!$           do j=1,ngy
!!$              xgc(i,j)=xg(i)
!!$              ygc(i,j)=yg(j)
!!$           enddo
!!$        enddo
!!$     endif

     !xg=xg*deltax
     !yg=yg*deltax

     xgc=xgc*deltax
     ygc=ygc*deltax

     Lxp=0.0
     Lyp=0.7096666666666666*deltax

     ! Tam & Dong radiation center (default) TO BE CHANGED
     !xcr_=0.0_wp
     xcr_=xgc(ngx/2,1)
     ycr_=ygc(ngx/2,1)
     xcr_=-0.015
     ycr_= 0.005

!!$     is_curv=.false.
!!$     print *,is_curv, xcr_,ycr_

 !*********************************************************************   
  case (2) ! 2 blocks (cart or curv) if non zero alp -> curv (2nd block is skewed)
  !*********************************************************************
!!$     ! skew angle
!!$     alp=0.0_wp*pi/180.0_wp
!!$     if (alp.ne.0.0_wp) is_curv=.true.
!!$     if (nob(iproc)==1) then
!!$        xg(1)=0.0_wp
!!$        do i=2,ngx
!!$           xg(i)=xg(i-1)+dx0
!!$        enddo
!!$        yg(1)=0.0_wp
!!$        do j=2,ngy
!!$           yg(j)=yg(j-1)+dy0
!!$           dy0=ry*dy0
!!$        enddo
!!$        ! convert to curvilinear grid
!!$        if (is_curv) then
!!$           do i=1,ngx
!!$              do j=1,ngy
!!$                 xgc(i,j)=xg(i)
!!$                 ygc(i,j)=yg(j)
!!$              enddo
!!$           enddo
!!$        endif
!!$     elseif (nob(iproc)==2) then
!!$        if (is_curv) then
!!$           x0=120.0_wp
!!$           if (is_adjoint_block) x0=119.0_wp
!!$           xgc(1,1)=x0
!!$           ygc(1,1)=0.0_wp
!!$           do i=2,ngx
!!$              xgc(i,1)=xgc(i-1,1)+dx0
!!$              ygc(i,1)=(xgc(i,1)-xgc(1,1))*tan(alp)
!!$           enddo
!!$           do j=2,ngy
!!$              xgc(:,j)=xgc(:,j-1)
!!$              ygc(:,j)=ygc(:,j-1)+dy0
!!$           enddo
!!$        else
!!$           xg(1)=120.0_wp
!!$           if (is_adjoint_block) xg(1)=119.0_wp
!!$           do i=2,ngx
!!$              xg(i)=xg(i-1)+dx0
!!$           enddo
!!$           yg(1)=0.0_wp
!!$           do j=2,ngy
!!$              yg(j)=yg(j-1)+dy0
!!$              dy0=ry*dy0
!!$           enddo
!!$        endif
!!$     endif
!!$     ! Tam & Dong radiation center (default) TO BE CHANGED
!!$     xcr_=100.0_wp
!!$     ycr_=60.0_wp

     if (is_curv3) then
        if (is_adjoint_block) then
           call grid_read('tv3_bl'//trim(numchar(nob(iproc)))//'.x')
        else
           call grid_read('tv3_mod_bl'//trim(numchar(nob(iproc)))//'.x')
        endif
     else
        is_curv=.true.
        !call grid_read('tr_bl'//trim(numchar(nob(iproc)))//'.x')
        !call grid_read('tv_bl'//trim(numchar(nob(iproc)))//'.x')

        if (is_adjoint_block) then
           call grid_read32('tv3_bl'//trim(numchar(nob(iproc)))//'.x')
        else
           call grid_read32('tv3_mod_bl'//trim(numchar(nob(iproc)))//'.x')
        endif
     endif

     deltax=1.0_wp
     
     ! Tam & Dong radiation center (default) TO BE CHANGED
     xcr_=0.0_wp
     ycr_=0.0_wp
     zcr_=0.0_wp
  !*********************************************************************   
  case (3) ! 4 blocks acoustics scatt in cavity (cart but can be run curv)
  !*********************************************************************   
     dx0=deltax
     dy0=deltay
     if (nob(iproc)==1) then
        ! along x
        xg(1)=-60.0_wp
        do i=2,ngx
           xg(i)=xg(i-1)+dx0
        enddo
        ! along y
        yg(1)=0.0_wp
        do j=2,ngy
           yg(j)=yg(j-1)+dy0
           dy0=ry*dy0
        enddo
     elseif (nob(iproc)==2) then
        ! along x
        xg(1)=0.0_wp
        do i=2,ngx
           xg(i)=xg(i-1)+dx0
        enddo
        ! along y
        yg(1)=0.0_wp
        do j=2,ngy
           yg(j)=yg(j-1)+dy0
           dy0=ry*dy0
        enddo
     elseif (nob(iproc)==3) then
        ! along x
        xg(1)=120.0_wp
        do i=2,ngx
           xg(i)=xg(i-1)+dx0
        enddo
        ! along y
        yg(1)=0.0_wp
        do j=2,ngy
           yg(j)=yg(j-1)+dy0
           dy0=ry*dy0
        enddo
     elseif (nob(iproc)==4) then
        ! along x
        xg(1)=0.0_wp
        do i=2,ngx
           xg(i)=xg(i-1)+dx0
        enddo
        ! along y
        yg(1)=-60.0_wp
        do j=2,ngy
           yg(j)=yg(j-1)+dy0
           dy0=ry*dy0
        enddo
     endif
     ! convert to curvilinear grid
     if (is_curv) then
        do i=1,ngx
           do j=1,ngy
              xgc(i,j)=xg(i)
              ygc(i,j)=yg(j)
           enddo
        enddo
     endif
     ! Tam & Dong radiation center (default) TO BE CHANGED
     xcr_=119._wp
     ycr_=0.0_wp
  !*********************************************************************   
  case (4) ! 4 blocks low-Re cavity flow (cart but can be run in curv)
  !*********************************************************************   
     ! Cavity aspect ratio: Length/Depth
     LsD=2.0_wp
     ! Cavity depth --> ReD=1500
     Depth=1.927472504432046e-4_wp
     ! Cavity depth --> ReD=2000
     Depth=2.561992378660378e-4_wp
     ! Cavity length
     Length=LsD*Depth
     ! Stretching rates
     rx1=1.04_wp
     rx3=1.028314333_wp
     rx4=1.04_wp
     rx5=1.06_wp
     ry1=1.03202_wp
     ry2=1.038_wp
     ry3=1.04_wp

     if (nob(iproc)==1) then
        ! along x
        xg(ngx)=0.0_wp-dx0
        dx0=rx1*dx0
        do i=ngx-1,1,-1
           xg(i) = xg(i+1)-dx0 
           dx0=rx1*dx0
        enddo
        ! along y
        yg(1)=0.0_wp
        do j=2,ngy
           yg(j)=yg(j-1)+dy0
           dy0=ry2*dy0
        enddo
     elseif (nob(iproc)==2) then
        ! bissection to determine rx3
        hxm=ngx/2
        rxa=1.0000000001
        rxb=1.5
        it=1
        epsi=1e-12
        erreur=1.
        do while (erreur>epsi)
           rxc=(rxa+rxb)/2.

           if (f_rx(rxa)*f_rx(rxc)>0) then
              rxa=rxc
              rxb=rxb
           else
              rxa=rxa
              rxb=rxc
           endif
           erreur=abs(f_rx(rxc))
           it=it+1
        enddo
        !print *,rx3,rxc,it
        rx3=rxc
        ! along x
        xg(1)=0.0_wp
        do i=2,hxm
           xg(i)=xg(i-1)+dx0
           dx0=rx3*dx0
        enddo
        do i=hxm+1,ngx
           xg(i)=xg(i-1)+dx0
           dx0=dx0/rx3
        enddo
        ! along y
        yg(1)=0.0_wp
        do j=2,ngy
           yg(j)=yg(j-1)+dy0
           dy0=ry2*dy0
        enddo
     elseif (nob(iproc)==3) then
        ! along x
        xg(1)=Length+dx0
        dx0=rx4*dx0
        do i=2,ngx-15
           xg(i)=xg(i-1)+dx0
           dx0=rx4*dx0
        enddo
        do i=ngx-14,ngx
           xg(i)=xg(i-1)+dx0
           dx0=rx5*dx0
        enddo
        ! along y
        yg(1)=0.0_wp
        do j=2,ngy
           yg(j)=yg(j-1)+dy0
           dy0=ry2*dy0
        enddo
     elseif (nob(iproc)==4) then
        ! bissection to determine rx3
        hxm=ngx/2
        rxa=1.0000000001
        rxb=1.5
        it=1
        epsi=1e-12
        erreur=1.
        do while (erreur>epsi)
           rxc=(rxa+rxb)/2.

           if (f_rx(rxa)*f_rx(rxc)>0) then
              rxa=rxc
              rxb=rxb
           else
              rxa=rxa
              rxb=rxc
           endif
           erreur=abs(f_rx(rxc))
           it=it+1
        enddo
        !print *,rx3,rxc,it
        rx3=rxc
        ! along x
        xg(1)=0.0_wp
        do i=2,hxm
           xg(i)=xg(i-1)+dx0
           dx0=rx3*dx0
        enddo
        do i=hxm+1,ngx
           xg(i)=xg(i-1)+dx0
           dx0=dx0/rx3
        enddo
        ! bissection to determine ry1
        hym=ngy/2
        rya=1.0000000001
        ryb=1.5
        it=1
        epsi=1e-12
        erreur=1.
        do while (erreur>epsi)
           ryc=(rya+ryb)/2.

           if (f_ry(rya)*f_ry(ryc)>0) then
              rya=ryc
              ryb=ryb
           else
              rya=rya
              ryb=ryc
           endif
           erreur=abs(f_ry(ryc))
           it=it+1
        enddo
        print *,ry1,ryc,it
        ry1=ryc
        ! along y
        yg(ngy)=0.0_wp-dy0
        dy0=ry1*dy0
        do j=ngy-1,hym+1,-1
           yg(j)=yg(j+1)-dy0
           dy0=ry1*dy0
        enddo
        do j=hym,1,-1
           yg(j)=yg(j+1)-dy0
           dy0=dy0/ry1
        enddo
     endif
     ! convert to curvilinear grid
     if (is_curv) then
        do i=1,ngx
           do j=1,ngy
              xgc(i,j)=xg(i)
              ygc(i,j)=yg(j)
           enddo
        enddo
     endif
     ! Tam & Dong radiation center (default) TO BE CHANGED
     xcr_=Length
     ycr_=0.0_wp
  !*********************************************************************   
  case (5) ! actuator (14 blocks) (necessarily curv)
  !*********************************************************************
     is_curv=.true.
     call grid_read('grid_actu/actu_bl'//trim(numchar(nob(iproc)))//'.x')
     deltax=0.1
     ! Tam & Dong radiation center (default) TO BE CHANGED
     xcr_=0.0_wp
     ycr_=31.7_wp
  !*********************************************************************   
  case (6) ! skewed grid & sinusoidal grid
  !*********************************************************************
!!$     ! skewed grid
!!$     alp=20.0_wp*pi/180.0_wp
!!$     x0=0.0_wp
!!$     xgc(1,1)=x0
!!$     ygc(1,1)=0.0_wp
!!$     Lxp=120.0_wp
!!$     Lyp=120.0_wp
!!$     do i=2,ngx+5
!!$        xgc(i,1)=xgc(i-1,1)+dx0
!!$        ygc(i,1)=(xgc(i,1)-xgc(1,1))*tan(alp)
!!$     enddo
!!$     do i=0,-4,-1
!!$        xgc(i,1)=xgc(i+1,1)-dx0
!!$        ygc(i,1)=(xgc(i,1)-xgc(1,1))*tan(alp)
!!$     enddo
!!$     do j=2,ngy+5
!!$        xgc(:,j)=xgc(:,j-1)
!!$        ygc(:,j)=ygc(:,j-1)+dy0
!!$     enddo
!!$     do j=0,-4,-1
!!$        xgc(:,j)=xgc(:,j+1)
!!$        ygc(:,j)=ygc(:,j+1)-dy0
!!$     enddo
     ! sinusoidal grid
!!$     xmin=-60.0_wp
!!$     ymin=-60.0_wp
!!$     xmax= 60.0_wp
!!$     ymax= 60.0_wp
     xmin=-50.0_wp
     ymin=-50.0_wp
     xmax= 50.0_wp
     ymax= 50.0_wp
     Lx= xmax-xmin
     Lxp=Lx
     Ly= ymax-ymin
     Lyp=Ly
     ax= 3.0_wp
     ay= 3.0_wp
     enex= 6.0_wp  ! 3 sine periods
     eney= 6.0_wp  ! 3 sine periods
     deltax= Lx/(ngx-1)
     deltay= Ly/(ngy-1)
     dx0=deltax
     dy0=deltay

     do i=-4,ngx+5
        do j=-4,ngy+5
           xgc(i,j)= xmin+dx0*((i-1)+ax*sin(enex*pi*(j-1)*dy0/Ly))
           ygc(i,j)= ymin+dy0*((j-1)+ay*sin(eney*pi*(i-1)*dx0/Lx))
        enddo
     enddo
     ! Tam & Dong radiation center (default) TO BE CHANGED
     xcr_=xgc(ngx/2,ngy/2)
     ycr_=ygc(ngx/2,ngy/2)    
     !xcr_=xgc(ngx/2,1)
     !ycr_=ygc(ngx/2,1)    
  !*********************************************************************   
  case (7) ! test 4 blocks skewed
  !*********************************************************************
     is_curv=.true.
     call grid_read('t2_bl'//trim(numchar(nob(iproc)))//'.x')
     deltax=1.0_wp

     xgc=xgc*1.e-4
     ygc=ygc*1.e-4
     deltax=1.0e-4_wp
     
     ! Tam & Dong radiation center (default) TO BE CHANGED
     xcr_=0.0_wp
     ycr_=0.0_wp
  !*********************************************************************   
  case (8) ! test 5 blocks multiple point
  !*********************************************************************
     is_curv=.true.
     call grid_read('t3_bl'//trim(numchar(nob(iproc)))//'.x')
     deltax=1.0_wp

     !xgc=xgc*1.e-4
     !ygc=ygc*1.e-4
     !deltax=1.0e-4_wp
     
     ! Tam & Dong radiation center (default) TO BE CHANGED
     xcr_=0.0_wp
     ycr_=0.0_wp
  !*********************************************************************
  case (9) ! 1 block sinusoidal grid
  !*********************************************************************
     xmin=-4.0_wp
     ymin=-4.0_wp
     zmin=-4.0_wp
     Lx=8.0_wp
     Lxp=Lx
     Ly=8.0_wp
     Lyp=Ly
     Lz=8.0_wp
     Lzp=Lz
     ax=1.0_wp
     ay=1.0_wp
     az=1.0_wp
     ! 4 sine periods
     enxy=8.0_wp
     enxz=8.0_wp
     enyx=8.0_wp
     enyz=8.0_wp
     enzx=8.0_wp
     enzy=8.0_wp
     deltax= Lx/ngx
     deltay= Ly/ngy
     deltaz= Lz/ngz
     dx0=deltax
     dy0=deltay
     dz0=deltaz

     do k=-4,ngz+5
        do j=-4,ngy+5
           do i=-4,ngx+5
              xgc3(i,j,k)=xmin+dx0*(dble(i-1)+ax*sin(enxy*pi*dble(j-1)*dy0/Ly)*sin(enxz*pi*dble(k-1)*dz0/Lz))
              ygc3(i,j,k)=ymin+dy0*(dble(j-1)+ay*sin(enyx*pi*dble(i-1)*dx0/Lx)*sin(enyz*pi*dble(k-1)*dz0/Lz))
              zgc3(i,j,k)=zmin+dz0*(dble(k-1)+az*sin(enzx*pi*dble(i-1)*dx0/Lx)*sin(enzy*pi*dble(j-1)*dy0/Ly))
!!$              xgc3(i,j,k)=xmin+dx0*(dble(i-1)+ax*sin(enxy*pi*dble(j-1)*dy0/Ly))
!!$              ygc3(i,j,k)=ymin+dy0*(dble(j-1)+ay*sin(enyx*pi*dble(i-1)*dx0/Lx))
!!$              zgc3(i,j,k)=zmin+dz0*dble(k-1)
           enddo
        enddo
     enddo

     ! bypass of comm grid
     igrd=6

     ! Tam & Dong radiation center (default) TO BE CHANGED
     xcr_=0.0_wp
     ycr_=0.0_wp

  !*********************************************************************
  case (10) ! 5 blocks 3D nozzle
  !*********************************************************************
     if (is_adjoint_block) then
        call grid_read('cdnoz3D2_bl'//trim(numchar(nob(iproc)))//'.x')
     else
        !call grid_read('cdnoz3D2_mod_bl'//trim(numchar(nob(iproc)))//'.x')
        call grid_read_ex('cdnoz3D2_e_bl'//trim(numchar(nob(iproc)))//'.x')
        !call grid_read('cdnoz3D2_bl'//trim(numchar(nob(iproc)))//'.x')
     endif
     deltax=1.

     igrd=6

     ! Tam & Dong radiation center (default) TO BE CHANGED
     !xcr_=0.0_wp
     xcr_=xgc(ngx/2,1)
     ycr_=ygc(ngx/2,1)
  !*********************************************************************
  case (11) ! test 3 blocks multiple point for IRS
  !*********************************************************************
     is_curv=.true.
     if (is_adjoint_block) then
        call grid_read('tourb32_bl'//trim(numchar(nob(iproc)))//'.x')
     else
        call grid_read('tourb32_mod_bl'//trim(numchar(nob(iproc)))//'.x')
     endif

     deltax=1.0_wp

     ! Tam & Dong radiation center (default) TO BE CHANGED
     xcr_=0.0_wp
     ycr_=0.0_wp
  !*********************************************************************
  case (12) ! test 4 blocks curvilinear 3D
  !*********************************************************************
     if (is_adjoint_block) then
        call grid_read('test_bl'//trim(numchar(nob(iproc)))//'.x')
     else
        call grid_read('test_mod_bl'//trim(numchar(nob(iproc)))//'.x')
     endif

     deltax=1.0_wp

     ! Tam & Dong radiation center (default) TO BE CHANGED
     xcr_=0.0_wp
     ycr_=0.0_wp
     zcr_=0.0_wp
  !*********************************************************************
  case (13) ! test 4 blocks curvilinear 3D
  !*********************************************************************
     !call grid_read_ex('sphereBCc_e_bl'//trim(numchar(nob(iproc)))//'.x')
     call grid_read('sphereBCc_mod_bl'//trim(numchar(nob(iproc)))//'.x')
     !call grid_read('sphereBC2_mod_bl'//trim(numchar(nob(iproc)))//'.x')
     !call grid_read('BC_bl'//trim(numchar(nob(iproc)))//'.x')
     !igrd=6
     deltax=0.03_wp

     ! Tam & Dong radiation center (default) TO BE CHANGED
     xcr_=19.0_wp
     ycr_=0.0_wp
     zcr_=0.0_wp
  case default
     call mpistop('not defined!', 0)
  end select

  if (.not.is_curv3) then
     ! regular grid along z
     if (is_2D) then
        zg=0.0_wp
     else
        zg(ngz/2)=0.0_wp
        do k=ngz/2+1,ngz
           zg(k)=zg(k-1)+dz0
        enddo
        do k=ngz/2-1,1,-1
           zg(k)=zg(k+1)-dz0
        enddo
     endif
  endif

  ! ! ---------------------------------------------------------------------------
  ! ! Computation of extrapolation coefficients for borders_wall
  ! cextrp2 = (1.0_wp+ry)**2/(ry*(2.0_wp+ry))
  ! cextrp3 =  1.0_wp       /(ry*(2.0_wp+ry))
  ! if (iproc.eq.0) then
  !    write(*,*) 'Grid Stretching      :', ry
  !    write(*,*) 'c2/c3 for border_wall:', cextrp2, cextrp3
  ! endif
  ! ! ~> Done by init_coeff_wall_cart in mod_coeff_deriv.f90 now
  ! ---------------------------------------------------------------------------
  
contains

  !=================================================================================
  ! Bissection method to determine the stretching rate
  !=================================================================================
  function f_rx(rx)
    implicit none
    real(wp), intent(in) :: rx
    real(wp) :: f_rx

    f_rx=LsD-(deltax*(1-rx**(hxm-1))/(1-rx)/Depth &
         +deltax*rx**(hxm-1)*(1-(1/rx)**(ngx-hxm))/(1-1/rx)/Depth)

  end function f_rx
  !=================================================================================
  function f_ry(ry)
    implicit none
    real(wp), intent(in) :: ry
    real(wp) :: f_ry

    f_ry=1.-(deltay*(1.-ry**hym)/(1-ry)/Depth &
            +deltay*ry**hym*(1.-(1./ry)**hym)/(1-(1./ry))/Depth)

  end function f_ry
  !=================================================================================

end subroutine grid_src

!===============================================================================
subroutine grid_chan
!===============================================================================
  !> compute predefined grid for channel flow case
!===============================================================================
  use mpi
  use mod_flow
  use mod_constant
  use mod_mpi
  use warnstop
  implicit none
  ! ----------------------------------------------------------------------------
  integer :: i,j,k,hym,it
  real(wp) :: dx0,dy0,dz0,ry
  real(wp) :: erreur,rya,ryb,ryc,lplusm1
  real(wp), parameter :: threshold=1.e-12_wp
  ! ----------------------------------------------------------------------------
  real(wp) :: Lx,Ly,Lz

  ! Construction of global grid
  !----------------------------
  Lx = longx*pi*L_ref
  dx0 = Lx/dble(ngx)
  deltax = dx0

  Lz = longz*pi*L_ref
  dz0 = Lz/dble(ngz)
  deltaz = dz0
  if (is_2d) then
     deltaz=1.0_wp
     Lz=1.0_wp
  endif

  xmin = 0.0_wp
  do i=-4,ngx+5
     xg(i) = xmin + dble(i-1)*dx0
  enddo

  zmin = - Lz/2.0_wp
  do k=-4,ngz+5
     zg(k) = zmin + dble(k-1)*dz0
  enddo

  Ly = 2.0_wp*L_ref

  ! Construction of wall-normal grid
  ! --------------------------------
  if (is_2D) then
     dy0 = 2.*hc/dble(ngy)/2.
  else
     dy0 = dy0p*mu_ref/(rho_wall*utau)
  endif

  ! compute geometric stretching rate ry in wall normal direction using bissection
  ! ------------------------------------------------------------------------------
  ry  = 1.0_wp
  hym = ngy/2
  rya = 1.0000000001_wp
  ryb = 1.5_wp
  it = 1
  erreur = 1.0_wp

  do while ((erreur>threshold).and.(it<1000))
     ryc=(rya+ryb)/2.0_wp
     if (f_ry(rya)*f_ry(ryc)>0.0_wp) then
        rya=ryc
        ryb=ryb
     else
        rya=rya
        ryb=ryc
     endif
     erreur = abs(f_ry(ryc))
     it = it+1
  enddo

  if (it==1000) call mpistop('dichotomy for grid computation not converged!!!', 0)

  ry = ryc
  yg(1) = -L_ref
  ! Lower half
  do j=2,hym
     yg(j) = yg(j-1) + dy0
     dy0= ry*dy0
  enddo
  ! Upper half
  do j=hym+1,ngy
     yg(j)=abs(yg(ngy-j+1))
  enddo

  deltay=abs(yg(2)-yg(1))

  ! ! Computation of extrapolation coefficients for borders_wall
  ! cextrp2 = (1.0_wp+ry)**2/(ry*(2.0_wp+ry))
  ! cextrp3 =  1.0_wp       /(ry*(2.0_wp+ry))
  ! if (iproc.eq.0) then
  !    write(*,*) 'Grid Stretching, iter:', ry, it
  !    write(*,*) 'c2/c3 for border_wall:', cextrp2, cextrp3
  ! endif
  ! ! ~> Done by init_coeff_wall_cart in mod_coeff_deriv.f90 now

  ! curvilinear grid (for test)
  if (is_curv) then
     do i=-4,ngx+5
        do j=1,ngy
           xgc(i,j)=xg(i)
           ygc(i,j)=yg(j)
        enddo
     enddo
  endif
     
  ! ---------------------------------------------------------------------------
  if (iproc==0) then
     lplusm1 = utau/mu_ref*rho_wall
     write(*,*) 'dx+ :', (xg(2)-xg(1))*lplusm1
     write(*,*) 'dyw+:', (yg(2)-yg(1))*lplusm1
     write(*,*) 'dz+ :', (zg(2)-zg(1))*lplusm1
     write(*,*) 'dyc+:', (yg(hym+1)-yg(hym))*lplusm1
  endif

  xmin = xg(1); xmax = xg(ngx)
  ymin = yg(1); ymax = yg(ngy)
  zmin = zg(1); zmax = zg(ngz)
  
contains

  !=================================================================================
  ! Bissection method to determine the wall-normal stretching rate
  !=================================================================================
  function f_ry(ry)
    use mod_constant
    implicit none
    real(wp), intent(in) :: ry
    real(wp) :: f_ry

    f_ry=1.0_wp-(dy0*(1.0_wp-ry**(hym-1.0_wp))/(1.0_wp-ry)/Ly &
         +dy0*ry**(hym-1.0_wp)*(1.0_wp-(1.0_wp/ry)**(ngy-hym))/(1.0_wp-1.0_wp/ry)/Ly)

  end function f_ry
  !=================================================================================

end subroutine grid_chan

!===============================================================================
subroutine grid_phill
!===============================================================================
  !> read predefined grid for periodic hill flows benchmark
!===============================================================================
  use mpi
  use mod_block
  use mod_flow
  use mod_constant
  use mod_mpi
  use warnstop
  use mod_comm
  use mod_bc_periodicity ! <~ needed for Lxp,Lyp,Lzp  TO BE CHANGED
  implicit none
  ! ----------------------------------------------------------------------------
  integer :: i,j,k
  ! ----------------------------------------------------------------------------
  real(wp) :: Lx,Ly,Lz
 
  if (iproc==0) print *, '  .. read periodic hill grid'
  
  ! Read global grid
  ! ================
  if (nbloc==1) then
     if (ngx==64) then
        open(10,file='grid_phill/phill_64x33x32_FD.dat',form='formatted',status='unknown')
     elseif (ngx==128) then
        open(10,file='grid_phill/phill_128x64x64_FD.dat',form='formatted',status='unknown')
     elseif (ngx==256) then
        open(10,file='grid_phill/phill_256x128x128_FD.dat',form='formatted',status='unknown')
     elseif (ngx==512) then
        open(10,file='grid_phill/phill_512x256x256_FD.dat',form='formatted',status='unknown')
     else
        call mpistop('no predefined grid for this value of ngx', 0)
     endif
     rewind(10)
     do j=1,ngy
        do i=1,ngx
           read(10,*) xgc(i,j),ygc(i,j)
        enddo
     enddo
     do k=1,ngz
        read(10,*) zg(k)
     enddo
     close(10)

!!$     ! used to test periodic hill on two blocks (dev XG)
!!$     open(10,file='grid_phill/phill_256x128x128_bl1.dat',form='formatted',status='unknown')
!!$     rewind(10)
!!$     do j=1,ngy
!!$        do i=1,ngx/2
!!$           write(10,*) xgc(i,j),ygc(i,j)
!!$        enddo
!!$     enddo
!!$     do k=1,ngz
!!$        write(10,*) zg(k)
!!$     enddo
!!$     close(10)
!!$     open(10,file='grid_phill/phill_256x128x128_bl2.dat',form='formatted',status='unknown')
!!$     rewind(10)
!!$     do j=1,ngy
!!$        do i=ngx/2+1,ngx
!!$           write(10,*) xgc(i,j),ygc(i,j)
!!$        enddo
!!$     enddo
!!$     do k=1,ngz
!!$        write(10,*) zg(k)
!!$     enddo
!!$     close(10)
!!$     stop
  else
     if (nob(iproc)==1) then
        if (ngx==128/2) then
           open(10,file='grid_phill/phill_128x64x64_bl1.dat',form='formatted',status='unknown')
        elseif (ngx==256/2) then
           open(10,file='grid_phill/phill_256x128x128_bl1.dat',form='formatted',status='unknown')
        else
           call mpistop('no predefined grid for this value of ngx', 0)
        endif
        rewind(10)
        do j=1,ngy
           do i=1,ngx
              read(10,*) xgc(i,j),ygc(i,j)
           enddo
        enddo
        do k=1,ngz
           read(10,*) zg(k)
        enddo
        close(10)
     elseif (nob(iproc)==2) then
        if (ngx==128/2) then
           open(10,file='grid_phill/phill_128x64x64_bl2.dat',form='formatted',status='unknown')
        elseif (ngx==256/2) then
           open(10,file='grid_phill/phill_256x128x128_bl2.dat',form='formatted',status='unknown')
        else
           call mpistop('no predefined grid for this value of ngx', 0)
        endif
        rewind(10)
        do j=1,ngy
           do i=1,ngx
              read(10,*) xgc(i,j),ygc(i,j)
           enddo
        enddo
        do k=1,ngz
           read(10,*) zg(k)
        enddo
        close(10)
     endif
  endif
  
  ! Dimensionalize grid
  ! ===================
  xgc=L_ref*xgc
  ygc=L_ref*ygc
  zg =L_ref*zg

  ! Domain sizes
  ! ============
  ! along x
  !Lx=xgc(ngx,ngy)-xgc(1,ngy)+xgc(2,ngy)-xgc(1,ngy)
  Lx=9.0_wp*L_ref
  Lxp=9.0_wp*L_ref
  ! along y
  !Ly=ygc(1,ngy)-ygc(1,1)
  Ly=2.035*L_ref
  ! along z
  if (is_2d) then
     deltaz=1.0_wp
     Lz=1.0_wp
  else
     deltaz=zg(2)-zg(1)
     Lz=ngz*deltaz
  endif
  if (iproc==0) print *,'before',Lx,Ly,Lz

!!$  ! Grid periodicity along x
!!$  ! ========================
!!$  xgc(ngx+1:ngx+5,:)=xgc(1:5,:)+Lx
!!$  ygc(ngx+1:ngx+5,:)=ygc(1:5,:)
!!$
!!$  xgc(-4:0,:)=xgc(ngx-4:ngx,:)-Lx
!!$  ygc(-4:0,:)=ygc(ngx-4:ngx,:)

  ! Grid periodicity along z
  ! ========================
  if (.not.(is_2d)) then
     zg(ngz+1:ngz+5)=zg(1:5)+Lz
     zg(-4:0)=zg(ngz-4:ngz)-Lz
  endif

!!$  ! false Cartesian counterpart for tests
!!$  do i=-4,ngx+5
!!$     xg(i)=xgc(i,ny)
!!$  enddo
!!$  do j=1,ngy
!!$     yg(j)=ygc(1,j)
!!$  enddo

end subroutine grid_phill

!===============================================================================
subroutine grid_stbl
!===============================================================================
  !> compute predefined grid for boundary layer cases
!===============================================================================
  use mod_block
  use mod_flow
  use mod_constant
  use mod_mpi
  implicit none
  ! ----------------------------------------------------------------------------
  integer :: i,j,k
  ! integer :: nxsv_1,nxsv_2,nysv_1,nysv_2,nzsv_1,nzsv_2,iv,ib
  real(wp) :: dx0,dy0,dz0
  real(wp) :: ry!,rx
  ! real(wp) :: ry0,deltax1 ! for hypersonic DG
  ! ----------------------------------------------------------------------------
  ! old M09
  ! integer :: ngxr,ngyr,ngzr
  ! real, dimension(:), allocatable :: xgr,ygr,zgr
  ! real, dimension(:,:), allocatable :: xgrc,ygrc
  ! ----------------------------------------------------------------------------
  real(wp) :: Lx,Ly,Lz
  ! ----------------------------------------------------------------------------

  ! ! For freestream turbulence induced transition post-treatment
  ! ! ===========================================================
  ! if (idepart.eq.POST_PROCESSING) then
  !    if (type_pp.eq.7) then
  !    ! if (is_RFM) then
  !       if (iproc.eq.0) print *,"Reading subvolume grid from file..."

  !       ! Grid is read and surface corresponding to the sub-volume is taken
  !       ! =================================================================
  !       ! Reading of all grid
  !       open(193,file='grid_bl1.bin',form='unformatted',status='unknown')
  !       rewind(193)
  !       read(193) ngxr
  !       read(193) ngyr
  !       read(193) ngzr
  !       if (is_curv) then
  !          allocate(xgrc(ngxr,ngyr),ygrc(ngxr,ngyr))
  !          read(193) ((xgrc(i,j),i=1,ngxr),j=1,ngyr)
  !          read(193) ((ygrc(i,j),i=1,ngxr),j=1,ngyr)
  !       else
  !          allocate(xgr(ngxr),ygr(ngyr))
  !          read(193) (xgr(i),i=1,ngx)
  !          read(193) (ygr(j),j=1,ngy)
  !       endif
  !       allocate(zgr(ngzr))
  !       read(193) (zgr(k),k=1,ngzr)
  !       close(193)

  !       ! Selection of the sub-volume grid
  !       ib = iblc_pp; iv = nplr
  !       nxsv_1 = bl(ib)%volume(iv)%ind_i1; nxsv_2 = bl(ib)%volume(iv)%ind_i2
  !       nysv_1 = bl(ib)%volume(iv)%ind_j1; nysv_2 = bl(ib)%volume(iv)%ind_j2
  !       nzsv_1 = bl(ib)%volume(iv)%ind_k1; nzsv_2 = bl(ib)%volume(iv)%ind_k2

  !       if (is_curv) then
  !          xgc = xgrc(nxsv_1:nxsv_2,nysv_1:nysv_2)
  !          ygc = ygrc(nxsv_1:nxsv_2,nysv_1:nysv_2)
  !          deallocate(xgrc,ygrc)
  !       else
  !          xg(nxsv_1:nxsv_2) = xgr(nxsv_1:nxsv_2)
  !          yg(nysv_1:nysv_2) = ygr(nysv_1:nysv_2)
  !          deallocate(xgr,ygr)
  !       endif
  !       zg(nzsv_1:nzsv_2) = zgr(nzsv_1:nzsv_2)
  !       deallocate(zgr)
  !       return
  !    endif
  ! endif

  ! Stretching rate
  ! ===============
  !ry=1.0_wp
  ry=1.015_wp

!!$  ! Hypersonic DG case (x-direction)
!!$  ! ==================
!!$  if (is_leading_edge) then
!!$     deltax1=3.6e-7_wp
!!$
!!$     xg(6)=0.0_wp
!!$     dx0=deltax1/20.0_wp
!!$     do i=5,1,-1
!!$        xg(i)=xg(i+1)-dx0
!!$        dx0=1.2_wp*dx0
!!$     enddo
!!$     dx0=deltax1/20.0_wp
!!$     do i=7,ngx
!!$        !if (xg(i-1)<7.9e-4_wp) then
!!$        if (xg(i-1)<12.e-4_wp) then
!!$           if (dx0<deltax1) then
!!$              xg(i)=xg(i-1)+dx0
!!$              dx0=1.2_wp*dx0
!!$           else
!!$              xg(i)=xg(i-1)+dx0
!!$           endif
!!$        else
!!$           if (dx0>deltax) then
!!$              xg(i)=xg(i-1)+dx0
!!$              dx0=0.98_wp*dx0
!!$           else
!!$              xg(i)=xg(i-1)+deltax
!!$           endif
!!$        endif
!!$     enddo
!!$     !$! SUPERSONIC
!!$     !$!      xg(6) = 0.0_wp
!!$     !$!      dx0=deltax/20.0_wp
!!$     !$!      do i=5,1,-1
!!$     !$!         xg(i) = xg(i+1) - dx0
!!$     !$!         dx0 = 1.2_wp * dx0
!!$     !$!      enddo
!!$     !$!      dx0=deltax/20.0_wp
!!$     !$!      do i=7,ngx
!!$     !$!         if (dx0<deltax) then
!!$     !$!            xg(i) = xg(i-1) + dx0
!!$     !$!            dx0 = 1.2_wp * dx0
!!$     !$!         else
!!$     !$!            xg(i) = xg(i-1) + dx0
!!$     !$!         endif
!!$     !$!      enddo
!!$  else
!!$     dx0=deltax
!!$     xg(1)=0.0_wp
!!$     do i=2,ngx
!!$        xg(i)=xg(i-1) + dx0
!!$     enddo
!!$  endif
!!$
!!$  ! ---------------------------------------------------------------------------
!!$
!!$  ! ! do i=ngx-49,ngx-10
!!$  ! !    xg(i) = xg(i-1) + dx0
!!$  ! !    !dx0 = 1.03 * dx0
!!$  ! !    !write(20,*) i,xg(i),dx0
!!$  ! ! enddo
!!$
!!$  ! ! do i=ngx-9,ngx
!!$  ! !    xg(i) = xg(i-1) + dx0
!!$  ! !    !dx0 = 1.1 * dx0
!!$  ! !    !write(20,*) i,xg(i),dx0
!!$  ! ! enddo
!!$
!!$  ! ry = 1.015_wp
!!$  ! deltax1=4.8e-7
!!$
!!$  ! xg(6)=0.0_wp
!!$  ! dx0=deltax1/20.0_wp
!!$  ! do i=5,1,-1
!!$  !    xg(i)=xg(i+1)-dx0
!!$  !    dx0=1.2_wp*dx0
!!$  ! enddo
!!$  ! dx0=deltax1/20.0_wp
!!$  ! do i=7,ngx
!!$  !    if (xg(i-1)<2.3e-4_wp) then
!!$  !       if (dx0<deltax1) then
!!$  !          xg(i)=xg(i-1)+dx0
!!$  !          dx0=1.2_wp*dx0
!!$  !       else
!!$  !          xg(i)=xg(i-1)+dx0
!!$  !       endif
!!$  !    else
!!$  !       if (dx0>deltax) then
!!$  !          xg(i)=xg(i-1)+dx0
!!$  !          dx0=0.98_wp*dx0
!!$  !       else
!!$  !          xg(i)=xg(i-1)+deltax
!!$  !       endif
!!$  !    endif
!!$  ! enddo
!!$  ! ---------------------------------------------------------------------------
!!$
!!$  ! dy0=deltay
!!$  ! yg(1) = 0.0_wp
!!$  ! do j=2,ngy
!!$  ! yg(j) = yg(j-1) + dy0
!!$  ! dy0 = ry * dy0
!!$  ! enddo
!!$  
!!$  ! Hypersonic DG case (y-direction)
!!$  ! ==================
!!$  dy0=deltay
!!$  yg(1)=0.0_wp
!!$  do j=2,ngy
!!$     yg(j)=yg(j-1)+dy0
!!$     !dy0=ry*dy0
!!$     if (yg(j)<2.4e-5_wp) dy0=ry*dy0
!!$  enddo
!!$  !$! SUPERSONIC
!!$  !$!   dy0=deltay
!!$  !$!   yg(1) = 0.0_wp
!!$  !$!   do j=2,ngy-30
!!$  !$!      yg(j) = yg(j-1) + dy0
!!$  !$!      dy0 = ry * dy0
!!$  !$!   enddo
!!$  !$!
!!$  !$!   ry0 = ry
!!$  !$!   do j=ngy-29,ngy
!!$  !$!      ry = ry0 + (1.05_wp-ry0)*dble(j-ngy+30)/30.0_wp
!!$  !$!      dy0 = ry * dy0
!!$  !$!      yg(j) = yg(j-1) + dy0
!!$  !$!   enddo

  ! simple grid along x
  ! ===================
!!$  dx0=deltax/3
!!$  rx=1.02
!!$  xg(1)=0.0_wp
!!$  do i=2,ngx
!!$     xg(i)=xg(i-1)+dx0
!!$     if (dx0<deltax) dx0=dx0*rx
!!$  enddo

  ! With sponge zone
  ! ================
  dx0=deltax
  xg(1)=0.0_wp
  do i=2,ngx-50
     xg(i)=xg(i-1)+dx0
  enddo

  do i=ngx-49,ngx-10
     xg(i)=xg(i-1)+dx0
     dx0=1.02*dx0
  enddo

  do i=ngx-9,ngx
     xg(i)=xg(i-1)+dx0
     dx0=1.06*dx0
  enddo

  ! print *,"             Modifier grid pour STBL !    "
  ! dx0=deltax
  ! xg(1)=0.0_wp
  ! do i=2,ngx
  !    xg(i)=xg(i-1)+dx0
  ! enddo

!!$  rx=1.1
!!$  do i=5,1,-1
!!$     dx0=dx0*rx
!!$     xg(i)=xg(i+1)-dx0
!!$  enddo

  ! ! simple grid along y
  ! ! ===================
  ! dy0=deltay
  ! !ry=1.01_wp
  ! yg(1)=0.0_wp
  ! do j=2,ngy
  !    yg(j)=yg(j-1)+dy0
  !    dy0=ry*dy0
  ! enddo

  ! ! GRID DNS Novec M=0.9
  ! ! ====================
  ! dy0=deltay
  ! ry=1.015_wp
  ! yg(1)=0.0_wp
  ! do j=2,190
  !    yg(j)=yg(j-1)+dy0
  !    dy0=ry*dy0
  ! enddo
  ! do j=191,278
  !    yg(j)=yg(j-1)+dy0
  ! enddo
  ! ry=1.025_wp
  ! do j=279,ngy
  !    yg(j)=yg(j-1)+dy0
  !    dy0=ry*dy0
  ! enddo

!!$  ! GRID LES Novec M=0.9
!!$  ! ====================
!!$  dy0=deltay
!!$  ry=1.015_wp
!!$  yg(1)=0.0_wp
!!$  do j=2,210
!!$     yg(j)=yg(j-1)+dy0
!!$     dy0=ry*dy0
!!$  enddo
!!$  do j=211,278
!!$     yg(j)=yg(j-1)+dy0
!!$  enddo
!!$  ry=1.025_wp
!!$  do j=279,ngy
!!$     yg(j)=yg(j-1)+dy0
!!$     dy0=ry*dy0
!!$  enddo

!!$  dy0=deltay
!!$  ry=1.02_wp
!!$  yg(1)=0.0_wp
!!$  do j=2,163
!!$     yg(j)=yg(j-1)+dy0
!!$     dy0=ry*dy0
!!$  enddo
!!$  do j=164,180
!!$     yg(j)=yg(j-1)+dy0
!!$  enddo

  ! ! FST transition T3A (+ wall-model test)
  ! ! ==================
  ! dy0=deltay
  ! ry=1.02_wp
  ! yg(1)=0.0_wp
  ! if (ngy-28.ge.128) then
  !    do j=2,128
  !       yg(j)=yg(j-1)+dy0
  !       dy0=ry*dy0
  !    enddo
  !    do j=128,ngy-28
  !       yg(j)=yg(j-1)+dy0
  !    enddo
  !    ry=1.025_wp
  !    do j=ngy-27,ngy
  !       yg(j)=yg(j-1)+dy0
  !       dy0=ry*dy0
  !    enddo
  ! else if (ngy.ge.128) then
  !     do j=2,128
  !        yg(j)=yg(j-1)+dy0
  !        dy0=ry*dy0
  !     enddo
  !     do j=128,ngy
  !        yg(j)=yg(j-1)+dy0
  !     enddo
  ! else
  !     do j=2,ngy
  !        yg(j)=yg(j-1)+dy0
  !        dy0=ry*dy0
  !     enddo
  !  endif

  ! ! Wall-model
  ! ! ==========
  ! dy0=deltay
  ! ry=1.02_wp
  ! yg(1)=0.0_wp
  ! do j=2,22
  !    yg(j)=yg(j-1)+dy0
  !    dy0=ry*dy0
  ! enddo
  ! do j=22,ngy-9
  !    yg(j)=yg(j-1)+dy0
  ! enddo
  ! ry=1.025_wp
  ! do j=ngy-9,ngy
  !    yg(j)=yg(j-1)+dy0
  !    dy0=ry*dy0
  ! enddo

  !! LES Novec M0p9 modes obliques
  !! =============================
  !dy0=deltay
  !ry=1.015_wp
  !yg(1)=0.0_wp
  !do j=2,210
  !   yg(j)=yg(j-1)+dy0
  !   dy0=ry*dy0
  !enddo
  !do j=211,278
  !   yg(j)=yg(j-1)+dy0
  !enddo
  !ry=1.025_wp
  !do j=279,ngy
  !   yg(j)=yg(j-1)+dy0
  !   dy0=ry*dy0
  !enddo

  ! LES Novec M0p9 FSTT T3A old /!\
  ! =======================
  ! jmax=200
  !dy0=deltay
  !ry=1.015_wp
  !yg(1)=0.0_wp
  !do j=2,160
  !  yg(j)=yg(j-1)+dy0
  !  dy0=ry*dy0
  !enddo
  !do j=161,170
  !  yg(j)=yg(j-1)+dy0
  !enddo
  !ry=1.025_wp
  !do j=171,ngy
  !  yg(j)=yg(j-1)+dy0
  !  dy0=ry*dy0
  !enddo

  ! LES Novec M0p9 FSTT T3A
  ! =======================
  ! jmax=280
  dy0=deltay
  ry=1.015_wp
  yg(1)=0.0_wp
  do j=2,175
    yg(j)=yg(j-1)+dy0
    dy0=ry*dy0
  enddo
  do j=176,250
    yg(j)=yg(j-1)+dy0
  enddo
  ry=1.025_wp
  do j=251,ngy
    yg(j)=yg(j-1)+dy0
    dy0=ry*dy0
  enddo

  ! ! LES Novec M0p9 FSTT T3A high Tu
  ! ! ===============================
  ! ! jmax=320
  ! dy0=deltay
  ! ry=1.015_wp
  ! yg(1)=0.0_wp
  ! do j=2,175
  !   yg(j)=yg(j-1)+dy0
  !   dy0=ry*dy0
  ! enddo
  ! do j=176,290
  !   yg(j)=yg(j-1)+dy0
  ! enddo
  ! ry=1.025_wp
  ! do j=291,ngy
  !   yg(j)=yg(j-1)+dy0
  !   dy0=ry*dy0
  ! enddo

  ! ! LES Novec M0p9 FSTT 10xT3A
  ! ! ==========================
  ! ! jmax=480
  ! dy0=deltay
  ! ry=1.015_wp
  ! yg(1)=0.0_wp
  ! do j=2,230
  !   yg(j)=yg(j-1)+dy0
  !   dy0=ry*dy0
  ! enddo
  ! do j=231,450
  !   yg(j)=yg(j-1)+dy0
  ! enddo
  ! ry=1.025_wp
  ! do j=451,ngy
  !   yg(j)=yg(j-1)+dy0
  !   dy0=ry*dy0
  ! enddo

  ! ! LES Novec M0p9 FSTT 10xT3A DNS
  ! ! ==============================
  ! ! jmax=840
  ! dy0=deltay
  ! ry=1.015_wp
  ! yg(1)=0.0_wp
  ! do j=2,170
  !   yg(j)=yg(j-1)+dy0
  !   dy0=ry*dy0
  ! enddo
  ! do j=171,810
  !   yg(j)=yg(j-1)+dy0
  ! enddo
  !   ry=1.025_wp
  ! do j=811,ngy
  !   yg(j)=yg(j-1)+dy0
  !   dy0=ry*dy0
  ! enddo



  ! ! LES Novec M0p9 FSTT 10xT3A DNS
  ! ! ==============================
  ! ! jmax=840
  ! dy0=deltay
  ! ry=1.015_wp
  ! yg(1)=0.0_wp
  ! do j=2,170
  !   yg(j)=yg(j-1)+dy0
  !   dy0=ry*dy0
  ! enddo
  ! do j=171,810
  !   yg(j)=yg(j-1)+dy0
  ! enddo
  ! ry=1.025_wp
  ! do j=811,ngy
  !   yg(j)=yg(j-1)+dy0
  !   dy0=ry*dy0
  ! enddo

  ! ! Temp
  ! ! ====
  ! ! jmax=440
  ! dy0=deltay
  ! ry=1.015_wp
  ! yg(1)=0.0_wp
  ! do j=2,170
  !   yg(j)=yg(j-1)+dy0
  !   dy0=ry*dy0
  ! enddo
  ! do j=171,410
  !   yg(j)=yg(j-1)+dy0
  ! enddo
  ! ry=1.025_wp
  ! do j=411,ngy
  !   yg(j)=yg(j-1)+dy0
  !   dy0=ry*dy0
  ! enddo

!!$  ! old M09 case
!!$  ! ============
!!$  delta=yg(jdel)
!!$
!!$  ! read old grid to avoid problem
!!$  open(193,file='maillage2dm1.bin',form='unformatted',status='unknown')
!!$  rewind(193)
!!$  read(193) ngxr
!!$  read(193) ngyr
!!$  allocate(xgr(ngxr),ygr(ngyr))
!!$  read(193) (xgr(j),j=1,ngxr)
!!$  read(193) (ygr(j),j=1,ngyr)
!!$  close(193)
!!$  xg=xgr(1:ngx)
!!$  yg=ygr(1:ngy)
!!$  deallocate(xgr,ygr)
!!$  deltax=abs(xg(101)-xg(100))
!!$  deltay=abs(yg(2)-yg(1))
!!$
!!$  delta=yg(jdel)-yg(18)

!!$  ! Read POINTWISE grid
!!$  ! -------------------
!!$  call grid_read('fp_bl1.x')
  
  ! ! Computation of extrapolation coefficients for borders_wall
  ! ! ==========================================================
  ! cextrp2 = (1.0_wp+ry)**2/(ry*(2.0_wp+ry))
  ! cextrp3 =  1.0_wp       /(ry*(2.0_wp+ry))
  ! if (iproc.eq.0) then
  !    write(*,*) 'Grid Stretching      :', ry
  !    write(*,*) 'c2/c3 for border_wall:', cextrp2, cextrp3
  ! endif
  ! ! ~> Done by init_coeff_wall_cart in mod_coeff_deriv.f90 now

  ! simple grid along z
  ! ===================
  dz0=deltaz
  if (is_2D) then
     zg=0.0_wp
  else
     zg(ngz/2)=0.0_wp
     do k=ngz/2+1,ngz
        zg(k)=zg(k-1)+dz0
     enddo
     do k=ngz/2-1,1,-1
        zg(k)=zg(k+1)-dz0
     enddo
  endif

  ! Domain sizes
  ! ============
  Lx = xg(ngx)-xg(1)
  Ly = yg(ngy)-yg(1)
  Lz = zg(ngz)-zg(1)

  ! Tam & Dong radiation center (default) TO BE CHANGED
  xcr_=xg(ngx/2)
  xcr_=xg(10)
  ycr_=yg(1)    

end subroutine grid_stbl

!===============================================================================
subroutine grid_cyl
!===============================================================================
  !> compute predefined grid for cylinder cases
!===============================================================================
  use mod_block
  use mod_flow
  use mod_constant
  use mod_mpi_part
  use mod_utils
  use warnstop
  implicit none
  ! ----------------------------------------------------------------------------
  integer :: i,j,k
  real(wp) :: rr,dr,dteta,dx,dz,dphi
  real(wp) :: radius
  real(wp), dimension(:), allocatable :: r ! radial grid
  ! ----------------------------------------------------------------------------
  real(wp) :: Lx,Ly,Lz
  ! ----------------------------------------------------------------------------
  real(wp) :: rra,rrb,rrc,erreur
  integer :: it

  if (is_curv3) then
  !if ((is_curv3).or.(is_curv)) then

     if (is_adjoint_block) then
        call grid_read('sphere5_bl'//trim(numchar(nob(iproc)))//'.x')
     else
        !call grid_read('Grid5/sphere5_mod_bl'//trim(numchar(nob(iproc)))//'.x')
        !call grid_read_ex('Grid/sphere5_e_bl'//trim(numchar(nob(iproc)))//'.x')
        !call grid_read('Grid3/sphere3_mod_bl'//trim(numchar(nob(iproc)))//'.x')
        !call grid_read_ex('sphere3_e_bl'//trim(numchar(nob(iproc)))//'.x')
        !call grid_read('Grid_50x50x100/sphere2_mod_bl'//trim(numchar(nob(iproc)))//'.x')
        call grid_read_ex('sphere2_e_bl'//trim(numchar(nob(iproc)))//'.x')
        !call grid_read('sphere_mod_bl'//trim(numchar(nob(iproc)))//'.x')
     endif

     ! bypass of comm grid
     igrd=6

!!$     if (is_adjoint_block) then
!!$        call grid_read('cyl3_bl'//trim(numchar(nob(iproc)))//'.x')
!!$     else
!!$        call grid_read('cyl3_mod_bl'//trim(numchar(nob(iproc)))//'.x')
!!$     endif

     ! Dimensionalize grid
     ! -------------------
     xgc3=xgc3*L_ref
     ygc3=ygc3*L_ref
     zgc3=zgc3*L_ref
     deltax=0.1*L_ref

!!$     dz=L_ref/ngz
!!$     zgc3(:,:,ngz/2)=0.0_wp
!!$     do k=ngz/2+1,ngz
!!$        zgc3(:,:,k)=zgc3(:,:,k-1)+dz
!!$     enddo
!!$     do k=ngz/2-1,1,-1
!!$        zgc3(:,:,k)=zgc3(:,:,k+1)-dz
!!$     enddo

!!$     if (is_adjoint_block) then
!!$        call grid_read('cyl3_bl'//trim(numchar(nob(iproc)))//'.x')
!!$     else
!!$        call grid_read('cyl3_mod_bl'//trim(numchar(nob(iproc)))//'.x')
!!$     endif
!!$
!!$     ! Dimensionalize grid
!!$     ! -------------------
!!$     xgc=xgc*L_ref
!!$     ygc=ygc*L_ref
!!$     deltax=0.1*L_ref
!!$
!!$     dz=L_ref/ngz
!!$     if (is_2D) then
!!$        zg=0.0_wp
!!$     else
!!$        zg(ngz/2)=0.0_wp
!!$        do k=ngz/2+1,ngz
!!$           zg(k)=zg(k-1)+dz
!!$        enddo
!!$        do k=ngz/2-1,1,-1
!!$           zg(k)=zg(k+1)-dz
!!$        enddo
!!$     endif

!!$     ! Sphere radius =diameter(Lref)/2
!!$     ! -------------------------------
!!$     radius=L_ref/2.0_wp
!!$
!!$     ! Stretching rate in radial direction
!!$     ! -----------------------------------
!!$     rr=1.02_wp
!!$     dr=0.05*L_ref
!!$
!!$     ! Radial distribution
!!$     ! -------------------
!!$     allocate(r(ngy))
!!$     r(1)=radius
!!$     do j=2,ngy
!!$        r(j)=r(j-1)+dr
!!$        dr=dr*rr
!!$     enddo
!!$
!!$     !print *,r(ngy)/L_ref
!!$
!!$     ! Azimuthal step
!!$     ! --------------
!!$     dphi=2.0_wp*pi/dble(ngx)
!!$
!!$     ! Spherical step
!!$     ! --------------
!!$     dteta=pi/dble(ngz-1)
!!$
!!$     ! Polar grid
!!$     ! ----------
!!$     do k=1,ngz
!!$        do i=-4,ngx+5
!!$           do j=1,ngy
!!$              zgc3(i,j,k)= r(j)*sin((k-1)*dteta)*cos((i-1)*dphi)
!!$              ygc3(i,j,k)= r(j)*sin((k-1)*dteta)*sin((i-1)*dphi)
!!$              xgc3(i,j,k)=-r(j)*cos((k-1)*dteta)
!!$           enddo
!!$        enddo
!!$     enddo

!!$     dteta=2.0_wp*pi/dble(ngx)
!!$
!!$     do i=-4,ngx+5
!!$        do j=1,ngy
!!$           xgc3(i,j,:)= r(j)*cos((i-1)*dteta+2*pi/2)
!!$           ygc3(i,j,:)= r(j)*sin((i-1)*dteta+2*pi/2)
!!$        enddo
!!$     enddo
!!$
!!$     dz=L_ref/ngz
!!$     zgc3(:,:,ngz/2)=0.0_wp
!!$     do k=ngz/2+1,ngz+5
!!$        zgc3(:,:,k)=zgc3(:,:,k-1)+dz
!!$     enddo
!!$     do k=ngz/2-1,-4,-1
!!$        zgc3(:,:,k)=zgc3(:,:,k+1)-dz
!!$     enddo

!!$     do i=-4,ngx+5
!!$        do j=1,ngy
!!$           xgc(i,j)= r(j)*cos((i-1)*dteta+2*pi/2)
!!$           ygc(i,j)= r(j)*sin((i-1)*dteta+2*pi/2)
!!$        enddo
!!$     enddo
!!$
!!$     dz=L_ref/ngz
!!$     if (is_2D) then
!!$        zg=0.0_wp
!!$     else
!!$        zg(ngz/2)=0.0_wp
!!$        do k=ngz/2+1,ngz
!!$           zg(k)=zg(k-1)+dz
!!$        enddo
!!$        do k=ngz/2-1,1,-1
!!$           zg(k)=zg(k+1)-dz
!!$        enddo
!!$     endif

!!$     do i=-4,ngx+5
!!$        do j=1,ngy
!!$           zgc3(i,j,:)= r(j)*cos((i-1)*dteta+2*pi/2)
!!$           ygc3(i,j,:)= r(j)*sin((i-1)*dteta+2*pi/2)
!!$        enddo
!!$     enddo
!!$
!!$     dz=L_ref/ngz
!!$     xgc3(:,:,ngz/2)=0.0_wp
!!$     do k=ngz/2+1,ngz+5
!!$        xgc3(:,:,k)=xgc3(:,:,k-1)+dz
!!$     enddo
!!$     do k=ngz/2-1,-4,-1
!!$        xgc3(:,:,k)=xgc3(:,:,k+1)-dz
!!$     enddo

!!$     dteta=2.0_wp*pi/dble(ngz)
!!$
!!$     do k=-4,ngz+5
!!$        do j=1,ngy
!!$           zgc3(:,j,k)= r(j)*cos((k-1)*dteta+2*pi/2)
!!$           ygc3(:,j,k)= r(j)*sin((k-1)*dteta+2*pi/2)
!!$        enddo
!!$     enddo
!!$
!!$     dx=L_ref/ngx
!!$     xgc3(ngx/2,:,:)=0.0_wp
!!$     do i=ngx/2+1,ngx+5
!!$        xgc3(i,:,:)=xgc3(i-1,:,:)+dx
!!$     enddo
!!$     do i=ngx/2-1,-4,-1
!!$        xgc3(i,:,:)=xgc3(i+1,:,:)-dx
!!$     enddo

!!$     ! bypass of comm grid
!!$     igrd=6

  else

  if (nbloc==1) then
     ! Case 1 block: POLAR GRID
     ! ========================

     ! Cylinder radius =diameter(Lref)/2
     ! ---------------------------------
     radius=L_ref/2.0_wp ! *40.

     ! Stretching rate in radial direction
     ! -----------------------------------
     rr=1.025_wp
     rr=1.02_wp
     rr=1.01_wp
     dr=radius/20.0_wp!/0.5!/10.
     rr=1.01_wp
     dr=radius/20.0_wp!/0.5!/10.
     !rr=1.0_wp
     !dr=1.0_wp
     ! for timestep calculation
     deltay=dr

     ! compute geometric stretching rate ry in wall normal direction using bissection
     ! ------------------------------------------------------------------------------
     dr=0.0028*L_ref
     rr =1.0_wp
     rra=1.0000000001_wp
     rrb=1.5_wp
     it=1
     erreur=1.0_wp

     do while ((erreur>1.0e-12_wp).and.(it<1000))
        rrc=(rra+rrb)/2.0_wp
        if (f_rr(rra)*f_rr(rrc)>0.0_wp) then
           rra=rrc
           rrb=rrb
        else
           rra=rra
           rrb=rrc
        endif
        erreur = abs(f_rr(rrc))
        it = it+1
     enddo
     rr=rrc

     print *,'conv. dichotomy',it,rr

     ! Radial distribution
     ! -------------------
     allocate(r(ngy))
     r(1)=radius
     do j=2,ngy
        r(j)=r(j-1)+dr
        dr=dr*rr
     enddo

     print *,r(ngy),20.*L_ref

!!$     do j=2,(ngy-1)/2
!!$        r(j)=r(j-1)+dr
!!$        dr=dr*rr
!!$     enddo
!!$     do j=(ngy-1)/2+1,ngy-20
!!$        r(j)=r(j-1)+dr
!!$     enddo
!!$     do j=ngy-19,ngy
!!$        r(j)=r(j-1)+dr
!!$        dr=dr*1.06_wp
!!$     enddo

     ! Azimuthal step
     ! --------------
     if (is_adjoint_block) then
        dteta=2.*pi/dble(ngx-1)
     else
        dteta=2.*pi/dble(ngx)
     endif

     ! Polar grid
     ! ----------
     do i=-4,ngx+5
        do j=1,ngy
           xgc(i,j)= r(j)*cos((i-1)*dteta+2*pi/2)
           ygc(i,j)= r(j)*sin((i-1)*dteta+2*pi/2)
        enddo
     enddo

     ! bypass of comm grid
     igrd=6

     ! ! Computation of extrapolation coefficients for borders_wall
     ! ! ==========================================================
     ! cextrp2 = (1.0_wp+rr)**2/(rr*(2.0_wp+rr))
     ! cextrp3 =  1.0_wp       /(rr*(2.0_wp+rr))
     ! if (iproc.eq.0) then
     !    write(*,*) 'Grid Stretching      :', rr
     !    write(*,*) 'c2/c3 for border_wall:', cextrp2, cextrp3
     ! endif
     ! ! ~> Done by init_coeff_wall_cart in mod_coeff_deriv.f90 now

  elseif (nbloc==6) then

     ! Case 6 blockss: C-H GRID
     ! ========================
     is_curv=.true.

     ! Read POINTWISE grid
     ! -------------------
!!     call grid_read('cyl_bl'//trim(numchar(nob(iproc)))//'.x')
     call grid_read('Grid_cyl_n/cyl_mod_bl'//trim(numchar(nob(iproc)))//'.x')

!!     ! Dimensionalize grid
!!     ! -------------------
!!     radius=L_ref/2.0_wp
!!     xgc=xgc*radius
!!     ygc=ygc*radius
!!     deltax=0.1*radius
     
     ! Dimensionalize grid
     ! -------------------
     xgc=xgc*L_ref
     ygc=ygc*L_ref
     deltax=0.1*L_ref
     
  elseif (nbloc==12) then
     
     ! Case 12 blocks: H-O-H GRID
     ! ========================
     is_curv=.true.

     ! Read POINTWISE grid
     ! -------------------
!!$     call grid_read('Grids_for_python_x/cyl_mod_bl'//trim(numchar(nob(iproc)))//'.x')
     
!!$     if (is_adjoint_block) then
!!$        call grid_read('cyl_blocks_F_large_O/cyl_bl'//trim(numchar(nob(iproc)))//'.x')
!!$     else
!!$        call grid_read('cyl_blocks_F_large_O/cyl_mod_bl'//trim(numchar(nob(iproc)))//'.x')
!!$     endif
     ! call grid_read('cyl_blocks_F_large_O_sponge_zone/cyl_mod_bl'//trim(numchar(nob(iproc)))//'.x')
     
     if (is_adjoint_block) then
        call grid_read('Grid/cylt_bl'//trim(numchar(nob(iproc)))//'.x')
     else
        call grid_read('Grid/cylt2_mod_bl'//trim(numchar(nob(iproc)))//'.x')
     endif
!!$     if (is_adjoint_block) then
!!$        call grid_read('Grid/cyl_bl'//trim(numchar(nob(iproc)))//'.x')
!!$     else
!!$        call grid_read('Grid/cyl_mod_bl'//trim(numchar(nob(iproc)))//'.x')
!!$     endif

     ! Dimensionalize grid
     ! -------------------
     xgc=xgc*L_ref!/20.0_wp
     ygc=ygc*L_ref!/20.0_wp

!!$     print *,L_ref*1e3,ygc(300,1)
!!$     if (iproc==0) then
!!$     do i=1,nx
!!$        print *,i,sqrt(xgc(i,1)**2+ygc(i,1)**2)
!!$     enddo
!!$     endif
!!$     call mpistop('grid',0)
     
     deltax=0.1*L_ref
  endif

  ! simple grid along z
  ! ===================
  dz=deltaz
  if (is_2D) then
     zg=0.0_wp
  else
     zg(ngz/2)=0.0_wp
     do k=ngz/2+1,ngz
        zg(k)=zg(k-1)+dz
     enddo
     do k=ngz/2-1,1,-1
        zg(k)=zg(k+1)-dz
     enddo
  endif

  endif

  ! Domain sizes
  ! ============
!!$  Lx = 2*pi
!!$  Ly = r(ngy)-r(1)
!!$  Lz = zg(ngz)-zg(1)

  Lx = 1.
  Ly = 1.
  Lz = 1.

  ! Tam & Dong radiation center (default) TO BE CHANGED
  !xcr_=20.0_wp
  xcr_=0.0_wp
  ycr_=0.0_wp
  zcr_=0.0_wp

contains
  
  !=============================================================================
  function f_rr(rr)
    implicit none
    real(wp), intent(in) :: rr
    real(wp) :: f_rr

    f_rr=20.0_wp*L_ref-dr*(1.-rr**(ngy-1))/(1.-rr)

  end function f_rr
  !=============================================================================
  
end subroutine grid_cyl

!===============================================================================
subroutine grid_turb
!===============================================================================
  !> compute predefined grid for turbine cases
!===============================================================================
  use mod_block
  use mod_flow
  use mod_constant
  use mod_mpi_part
  use mod_bc_periodicity ! <~ needed for Lxp,Lyp,Lzp  TO BE CHANGED
  use mod_io
  use mod_utils
  use warnstop
  implicit none
  ! ----------------------------------------------------------------------------
  integer :: i,j,k
  integer :: nxread,nyread
  real(wp) :: dx0,dy0,dz0,ry,chord
  ! ----------------------------------------------------------------------------
  real(wp) :: Lx,Ly,Lz
  ! annular duct case
  real(wp) :: Length,R1,R2,dr,dteta
  real(wp), dimension(:), allocatable :: r ! radial grid
  ! ----------------------------------------------------------------------------
  character(len=120) :: gridname
  ! ----------------------------------------------------------------------------

  if (nbloc==1) then ! annular duct
      ! Temporary ~> TO REMOVE /!\
      ! call grid_read('Grid_DDES_4interp/grid1_bl'//trim(numchar(nob(iproc)))//'.x')

     ! free vortex flow
     !Length=64.0e-3_wp*5.0_wp
     !R1=Length/2.0_wp
     !R2=R1+Length/4.0_wp

     R1=0.2_wp
     R2=0.28_wp
     Length=0.32_wp

     ! Axial direction
     ! ---------------
     dx0=Length/dble(ngx-1)
     xgc3(1,:,:)=0.0_wp
     do i=2,ngx
        xgc3(i,:,:)=xgc3(i-1,:,:)+dx0
     enddo

     ! Radial direction
     ! ----------------
     dr=(R2-R1)/dble(ngy-1)
     ! for timestep calculation (?)
     deltay=dr

     ! Radial distribution
     ! -------------------
     allocate(r(ngy))
     r(1)=R1
     do j=2,ngy
        r(j)=r(j-1)+dr
     enddo

     ! Azimuthal step
     ! --------------
     if (is_adjoint_block) then
        dteta=2.*pi/dble(ngz-1)
     else
        if (theta_period==0.0_wp) then
           dteta=2.*pi/dble(ngz)
        else
           dteta=theta_period*pi/180.0_wp/dble(ngz)
        endif
     endif

     ! Polar grid
     ! ----------
     do k=-4,ngz+5
        do j=1,ngy
           ygc3(:,j,k)= r(j)*cos((k-1)*dteta+pi/4.0_wp)
           zgc3(:,j,k)= r(j)*sin((k-1)*dteta+pi/4.0_wp)
        enddo
     enddo

     ! bypass of comm grid
     igrd=6

     ! Tam & Dong radiation center (default) TO BE CHANGED
     xcr_=xgc3(ngx/2,1,1)
     ycr_=0.0_wp

  elseif (nbloc==2) then ! sector of annular duct / 2 blocks

     if (is_adjoint_block) then
        call grid_read('ann_sect_bl'//trim(numchar(nob(iproc)))//'.x')
     else
        !call grid_read('ann_sect_mod_bl'//trim(numchar(nob(iproc)))//'.x')
        !call grid_read_ex('ann_sect_mod_e_bl'//trim(numchar(nob(iproc)))//'.x')

        !call grid_read(trim(dirGRID)//'/'//trim(nameGRID)// &
        !               '_mod_bl'//trim(numchar(nob(iproc)))//'.x')

        call grid_read_ex(trim(dirGRID)//'/'//trim(nameGRID)// &
                          '_ngh'//trim(numchar(ngh))//'_mod_bl'//trim(numchar(nob(iproc)))//'.x')
     endif

     ! bypass of comm grid
     igrd=6

     ! Tam & Dong radiation center (default) TO BE CHANGED
     xcr_=xgc3(ngx/2,1,1)
     ycr_=0.0_wp

  elseif (nbloc==4) then

!!$     dx0 = deltax
!!$     dy0 = deltay/5.0_wp
!!$     ry=1.02_wp
!!$
!!$     ! Semi-infinite plane (half_plane)
!!$     ! ===================
!!$     ! Bloc #1
!!$     if (nob(iproc)==1) then
!!$        ! along x
!!$        xg(1)=0.0_wp
!!$        do i=2,ngx
!!$           xg(i)=xg(i-1)+dx0
!!$        enddo
!!$        ! along y
!!$        yg(1)=0.0_wp
!!$        do j=2,ngy
!!$           yg(j)=yg(j-1)+dy0
!!$           dy0=ry*dy0
!!$        enddo
!!$     elseif (nob(iproc)==2) then
!!$        ! along x
!!$        xg(ngx)=0.0_wp
!!$        do i=ngx-1,1,-1
!!$           xg(i)=xg(i+1)-dx0
!!$        enddo
!!$        ! along y
!!$        yg(1)=0.0_wp
!!$        do j=2,ngy
!!$           yg(j)=yg(j-1)+dy0
!!$           dy0=ry*dy0
!!$        enddo
!!$     elseif (nob(iproc)==3) then
!!$        ! along x
!!$        xg(ngx)=0.0_wp
!!$        do i=ngx-1,1,-1
!!$           xg(i)=xg(i+1)-dx0
!!$        enddo
!!$        ! along y
!!$        yg(ngy)=0.0_wp
!!$        do j=ngy-1,1,-1
!!$           yg(j)=yg(j+1)-dy0
!!$           dy0=ry*dy0
!!$        enddo
!!$     elseif (nob(iproc)==4) then
!!$        ! along x
!!$        xg(1)=0.0_wp
!!$        do i=2,ngx
!!$           xg(i)=xg(i-1)+dx0
!!$        enddo
!!$        ! along y
!!$        yg(ngy)=0.0_wp
!!$        do j=ngy-1,1,-1
!!$           yg(j)=yg(j+1)-dy0
!!$           dy0=ry*dy0
!!$        enddo
!!$     endif
!!$
!!$     ! convert to curvilinear grid
!!$     if (is_curv) then
!!$        do i=1,ngx
!!$           do j=1,ngy
!!$              xgc(i,j)=xg(i)
!!$              ygc(i,j)=yg(j)
!!$           enddo
!!$        enddo
!!$     endif

!!$     if (is_adjoint_block) then
!!$        call grid_read('grid/wedge_bl'//trim(numchar(nob(iproc)))//'.x')
!!$     else
!!$        call grid_read('grid/wedge_mod_bl'//trim(numchar(nob(iproc)))//'.x')
!!$     endif
!!$
!!$     ! Tam & Dong radiation center (default) TO BE CHANGED
!!$     xcr_=0.0_wp
!!$     ycr_=0.0_wp

     if (is_adjoint_block) then
        call grid_read('inlet_bl'//trim(numchar(nob(iproc)))//'.x')
     else
        call grid_read('inletc_mod_bl'//trim(numchar(nob(iproc)))//'.x')
     endif

     chord=60.e-3_wp

     xgc=xgc*chord
     ygc=ygc*chord

     Lxp=0.0_wp*chord
     Lyp=0.709666666666667_wp*chord

     ! Tam & Dong radiation center (default) TO BE CHANGED
     xcr_=0.0_wp
     ycr_=0.0_wp
     xcr_=-0.015
     ycr_= 0.005

  elseif (nbloc==6) then

     ! Thin plate
     ! ==========
     dx0 = deltax
     dy0 = deltay/1.0_wp
     ry=1.0_wp

!!$     if (nob(iproc)==1) then
!!$        ! along x
!!$        xg(1)=49.0_wp
!!$        do i=2,ngx
!!$           xg(i)=xg(i-1)+dx0
!!$        enddo
!!$        ! along y
!!$        yg(1)=0.0_wp
!!$        do j=2,ngy
!!$           yg(j)=yg(j-1)+dy0
!!$           dy0=ry*dy0
!!$        enddo
!!$     elseif (nob(iproc)==2) then
!!$        ! along x
!!$        xg(1)=0.0_wp
!!$        do i=2,ngx
!!$           xg(i)=xg(i-1)+dx0
!!$        enddo
!!$        ! along y
!!$        yg(1)=0.0_wp
!!$        do j=2,ngy
!!$           yg(j)=yg(j-1)+dy0
!!$           dy0=ry*dy0
!!$        enddo
!!$     elseif (nob(iproc)==3) then
!!$        ! along x
!!$        xg(ngx)=0.0_wp
!!$        do i=ngx-1,1,-1
!!$           xg(i)=xg(i+1)-dx0
!!$        enddo
!!$        ! along y
!!$        yg(1)=0.0_wp
!!$        do j=2,ngy
!!$           yg(j)=yg(j-1)+dy0
!!$           dy0=ry*dy0
!!$        enddo
!!$     elseif (nob(iproc)==4) then
!!$        ! along x
!!$        xg(ngx)=0.0_wp
!!$        do i=ngx-1,1,-1
!!$           xg(i)=xg(i+1)-dx0
!!$        enddo
!!$        ! along y
!!$        yg(ngy)=0.0_wp
!!$        do j=ngy-1,1,-1
!!$           yg(j)=yg(j+1)-dy0
!!$           dy0=ry*dy0
!!$        enddo
!!$     elseif (nob(iproc)==5) then
!!$        ! along x
!!$        xg(1)=0.0_wp
!!$        do i=2,ngx
!!$           xg(i)=xg(i-1)+dx0
!!$        enddo
!!$        ! along y
!!$        yg(ngy)=0.0_wp
!!$        do j=ngy-1,1,-1
!!$           yg(j)=yg(j+1)-dy0
!!$           dy0=ry*dy0
!!$        enddo
!!$     elseif (nob(iproc)==6) then
!!$        ! along x
!!$        xg(1)=49.0_wp
!!$        do i=2,ngx
!!$           xg(i)=xg(i-1)+dx0
!!$        enddo
!!$        ! along y
!!$        yg(ngy)=0.0_wp
!!$        do j=ngy-1,1,-1
!!$           yg(j)=yg(j+1)-dy0
!!$           dy0=ry*dy0
!!$        enddo
!!$     endif

!!$      if (nob(iproc)==1) then
!!$        ! along x
!!$        xg(1)=49.0_wp
!!$        do i=2,ngx+ngh
!!$           xg(i)=xg(i-1)+dx0
!!$        enddo
!!$        do i=0,1-ngh,-1
!!$           xg(i)=xg(i+1)-dx0
!!$        enddo
!!$        ! along y
!!$        yg(1)=0.0_wp
!!$        do j=2,ngy+ngh
!!$           yg(j)=yg(j-1)+dy0
!!$           dy0=ry*dy0
!!$        enddo
!!$        do j=0,1-ngh,-1
!!$           yg(j)=yg(j+1)-dy0
!!$           dy0=ry*dy0
!!$        enddo
!!$     elseif (nob(iproc)==2) then
!!$        ! along x
!!$        xg(1)=0.0_wp
!!$        do i=2,ngx+ngh
!!$           xg(i)=xg(i-1)+dx0
!!$        enddo
!!$        do i=0,1-ngh,-1
!!$           xg(i)=xg(i+1)-dx0
!!$        enddo
!!$        ! along y
!!$        yg(1)=0.0_wp
!!$        do j=2,ngy+ngh
!!$           yg(j)=yg(j-1)+dy0
!!$           dy0=ry*dy0
!!$        enddo
!!$        do j=0,1-ngh,-1
!!$           yg(j)=yg(j+1)-dy0
!!$           dy0=ry*dy0
!!$        enddo
!!$     elseif (nob(iproc)==3) then
!!$        ! along x
!!$        xg(ngx)=0.0_wp
!!$        do i=ngx-1,1-ngh,-1
!!$           xg(i)=xg(i+1)-dx0
!!$        enddo
!!$        do i=ngx+1,ngx+ngh
!!$           xg(i)=xg(i-1)+dx0
!!$        enddo
!!$        ! along y
!!$        yg(1)=0.0_wp
!!$        do j=2,ngy+ngh
!!$           yg(j)=yg(j-1)+dy0
!!$           dy0=ry*dy0
!!$        enddo
!!$        do j=0,1-ngh,-1
!!$           yg(j)=yg(j+1)-dy0
!!$           dy0=ry*dy0
!!$        enddo
!!$     elseif (nob(iproc)==4) then
!!$        ! along x
!!$        xg(ngx)=0.0_wp
!!$        do i=ngx-1,1-ngh,-1
!!$           xg(i)=xg(i+1)-dx0
!!$        enddo
!!$        do i=ngx+1,ngx+ngh
!!$           xg(i)=xg(i-1)+dx0
!!$        enddo
!!$        ! along y
!!$        yg(ngy)=0.0_wp
!!$        do j=ngy-1,1-ngh,-1
!!$           yg(j)=yg(j+1)-dy0
!!$           dy0=ry*dy0
!!$        enddo
!!$        do j=ngy+1,ngy+ngh
!!$           yg(j)=yg(j-1)+dy0
!!$           dy0=ry*dy0
!!$        enddo
!!$     elseif (nob(iproc)==5) then
!!$        ! along x
!!$        xg(1)=0.0_wp
!!$        do i=2,ngx+ngh
!!$           xg(i)=xg(i-1)+dx0
!!$        enddo
!!$        do i=0,1-ngh,-1
!!$           xg(i)=xg(i+1)-dx0
!!$        enddo
!!$        ! along y
!!$        yg(ngy)=0.0_wp
!!$        do j=ngy-1,1-ngh,-1
!!$           yg(j)=yg(j+1)-dy0
!!$           dy0=ry*dy0
!!$        enddo
!!$        do j=ngy+1,ngy+ngh
!!$           yg(j)=yg(j-1)+dy0
!!$           dy0=ry*dy0
!!$        enddo
!!$     elseif (nob(iproc)==6) then
!!$        ! along x
!!$        xg(1)=49.0_wp
!!$        do i=2,ngx+ngh
!!$           xg(i)=xg(i-1)+dx0
!!$        enddo
!!$        do i=0,1-ngh,-1
!!$           xg(i)=xg(i+1)-dx0
!!$        enddo
!!$       ! along y
!!$        yg(ngy)=0.0_wp
!!$        do j=ngy-1,1-ngh,-1
!!$           yg(j)=yg(j+1)-dy0
!!$           dy0=ry*dy0
!!$        enddo
!!$        do j=ngy+1,ngy+ngh
!!$           yg(j)=yg(j-1)+dy0
!!$           dy0=ry*dy0
!!$        enddo
!!$     endif
!!$     ! periodicity already filled
!!$     igrd=6
!!$
!!$    ! convert to curvilinear grid
!!$     if (is_curv) then
!!$        do i=1-ngh,ngx+ngh
!!$           do j=1-ngh,ngy+ngh
!!$              xgc(i,j)=xg(i)
!!$              ygc(i,j)=yg(j)
!!$           enddo
!!$        enddo
!!$     endif

     ! Ailette
     ! =======
     open(194,file='grid/grid_bl'//trim(numchar(nob(iproc)))//'_ex.bin',form='unformatted',status='unknown')
     rewind(194)
     read(194) nxread
     read(194) nyread
     read(194) ((xgc(i,j),i=-4,ngx+5),j=-4,ngy+5)
     read(194) ((ygc(i,j),i=-4,ngx+5),j=-4,ngy+5)
     close(194)
     igrd=6
     xgc=xgc/1.0_wp
     ygc=ygc/1.0_wp

     ! Tam & Dong radiation center (default) TO BE CHANGED
     xcr_=0.0_wp
     ycr_=0.0_wp

  elseif (nbloc==1) then

     ! Nozzle flow (TROVA)
     ! ===================
     dx0 = deltax
     dy0 = deltay/1.0_wp
     ry=1.0_wp

     ! LS59-Euler
     ! ==========
     if (is_adjoint_block) then
        call grid_read('Grid_Euler/turb_bl'//trim(numchar(nob(iproc)))//'.x')
     else
        call grid_read('noz_mod_bl'//trim(numchar(nob(iproc)))//'.x')
     endif

     ! Tam & Dong radiation center (default) TO BE CHANGED
     xcr_=0.0_wp
     ycr_=0.0_wp

  elseif (nbloc==9) then
  !TEST elseif (nbloc==2) then

     ! Thin plate
     ! ==========
     dx0 = deltax
     dy0 = deltay/1.0_wp
     ry=1.0_wp

     ! LS59-Euler
     ! ==========
     if ((is_read_ex).and.(idepart.ne.1)) then

        ! init MPI-IO for grid read/write tata
        call init_io_grid_ex3d

        ! gridname
        gridname=trim(dirGRID)//'/'//trim(nameGRID)// &
             '_ngh'//trim(numchar(ngh))//'_mod_bl'//trim(numchar(nob(iproc)))//filext_read

        ! read extended grid
        call read_write_grid3d(trim(gridname),READ)

        ! free memory
        call free_grid_ex3d

     else
        if (is_adjoint_block) then
           call grid_read('Grid_Euler/turb_bl'//trim(numchar(nob(iproc)))//'.x')
        else
           !call grid_read('Grid_Euler/turb_mod_bl'//trim(numchar(nob(iproc)))//'.x')
           !call grid_read('Grid_RANS/turb_mod_bl'//trim(numchar(nob(iproc)))//'.x')
           !call grid_read('Grid_LES/turb_mod_bl'//trim(numchar(nob(iproc)))//'.x')
           !call grid_read('Grid_LES_GO/ls59_mod_bl'//trim(numchar(nob(iproc)))//'.x')
           !call grid_read('Grid_LES_GO2/ls59_mod_bl'//trim(numchar(nob(iproc)))//'.x')
           !call grid_read('Grid_LES_blunt/ls59c_mod_bl'//trim(numchar(nob(iproc)))//'.x')
           !call grid_read('Grid_LES_rough/ls59c_mod_bl'//trim(numchar(nob(iproc)))//'.x')
           !call grid_read('Grid_RANS_d_WakeF/ls59_mod_bl'//trim(numchar(nob(iproc)))//'.x')
           !call grid_read('Grid_RANS_CLOWT_WakeF/ls59_mod_bl'//trim(numchar(nob(iproc)))//'.x')
           !call grid_read('Grid_LES_GO_trip/ls59t_mod_bl'//trim(numchar(nob(iproc)))//'.x')
           !call grid_read('Grid_RANS_GO2/ls59_mod_bl'//trim(numchar(nob(iproc)))//'.x')
           !call grid_read('Grid_RANS_GO/ls593_mod_bl'//trim(numchar(nob(iproc)))//'.x')
           !call grid_read_ex('Grid_RANS_GO/ls593_mod_e_bl'//trim(numchar(nob(iproc)))//'.x')
           !call grid_read_ex('Grid_RANS_GO/ls593_mod_ngh'//trim(numchar(ngh))//'_bl'//trim(numchar(nob(iproc)))//'.x')
           !call grid_read_ex('Grid_RANS_GO/ls59w_mod_ngh'//trim(numchar(ngh))//'_bl'//trim(numchar(nob(iproc)))//'.x')
           !call grid_read('Grid_LES_CLOWT_fine/ls59_mod_bl'//trim(numchar(nob(iproc)))//'.x')
           !call grid_read('Grid_LES_CLOWT_fine2/ls59_mod_bl'//trim(numchar(nob(iproc)))//'.x')

           !call grid_read(trim(dirGRID)//'/'//trim(nameGRID)// &
           !               '_mod_bl'//trim(numchar(nob(iproc)))//'.x')

           call grid_read_ex(trim(dirGRID)//'/'//trim(nameGRID)// &
                '_ngh'//trim(numchar(ngh))//'_mod_bl'//trim(numchar(nob(iproc)))//'.x')

           !call grid_read('Grid_RANS_WakeF/ls59_mod_bl'//trim(numchar(nob(iproc)))//'.x')

           !call grid_read('Grid_RANS/turb_mod_bl'//trim(numchar(nob(iproc)))//'.x')
           !call grid_read_ex('Grid_RANS/turb_mod_e_bl'//trim(numchar(nob(iproc)))//'.x')
        endif
     endif

     ! bypass of comm grid
     if (is_curv3) igrd=6

     !chord=32.6e-4_wp
     !chord=100.e-3_wp
     ! GO
     chord=60.e-3_wp
     ! CLOWT
     chord=30.5e-3_wp
     ! axial chord Baumgartner
     !chord=16.e-3_wp
     chord=L_ref
     !chord=1.0_wp
     if ((idepart.eq.1).and.(is_curv3)) chord=1.0_wp

     xgc=xgc*chord
     ygc=ygc*chord
     if (is_curv3) then
        if ((is_read_ex).and.(idepart.ne.1)) then
           xc3=xc3*chord
           yc3=yc3*chord
           zc3=zc3*chord
        else
           xgc3=xgc3*chord
           ygc3=ygc3*chord
           zgc3=zgc3*chord
        endif
     endif

     !!!! Baumgartner !!!
     !Lxp=0.0_wp*L_ref
     !!!Lyp=0.9162978573_wp*L_ref ! pitch at mid-span
     !!!Lyp=0.7853981634_wp*L_ref ! pitch at hub
     !Lyp=1.047197551_wp*L_ref ! pitch at shroud

     Lxp=0.0_wp*chord
     Lyp=0.709666666666667_wp*chord
     !Lyp=0.721311475409836_wp*chord
     !Lyp=0.6954_wp*chord
     !Lyp=0.8495_wp*chord

     !Lxp=0.0_wp
     !Lyp=0.0_wp
     if (iproc==0) print *,'pitch',Lyp

     ! Tam & Dong radiation center (default) TO BE CHANGED
     xcr_=0.0_wp
     ycr_=0.0_wp
     !xcr_=-0.015
     !ycr_= 0.005

  elseif ((nbloc==12).or.(nbloc==13)) then

     ! Thin plate
     ! ==========
     dx0 = deltax
     dy0 = deltay/1.0_wp
     ry=1.0_wp

     ! LS59-Euler
     ! ==========
     if (is_adjoint_block) then
        call grid_read('Grid_LES_HU2/ls59_bl'//trim(numchar(nob(iproc)))//'.x')
     else
        call grid_read('Grid_LES_HU2/ls59c_mod_bl'//trim(numchar(nob(iproc)))//'.x')

        !call grid_read(trim(dirGRID)//'/'//trim(nameGRID)// &
        !               '_mod_bl'//trim(numchar(nob(iproc)))//'.x')

        !call grid_read_ex(trim(dirGRID)//'/'//trim(nameGRID)// &
        !                  '_ngh'//trim(numchar(ngh))//'_mod_bl'//trim(numchar(nob(iproc)))//'.x')
     endif

     ! bypass of comm grid
     !igrd=6

     ! HSU
     chord=61.2e-3_wp
     !chord=1.0_wp
     !chord=L_ref

     xgc=xgc*chord
     ygc=ygc*chord
     if (is_curv3) then
        xgc3=xgc3*chord
        ygc3=ygc3*chord
        zgc3=zgc3*chord
     endif

     Lxp=0.0_wp*chord
     Lyp=0.709666666666667_wp*chord
     if (iproc==0) print *,'pitch',Lyp

     ! Tam & Dong radiation center (default) TO BE CHANGED
     xcr_=0.0_wp
     ycr_=0.0_wp

  elseif (nbloc==14) then

     ! Thin plate
     ! ==========
     dx0 = deltax
     dy0 = deltay/1.0_wp
     ry=1.0_wp

     ! LS89 3D vane passage
     ! =====================
     if (is_adjoint_block) then
        call grid_read('Grid_stator/stator_bl'//trim(numchar(nob(iproc)))//'.x')
     else
        call grid_read('Grid_stator/stator_mod_bl'//trim(numchar(nob(iproc)))//'.x')
        !call grid_read_ex('Grid_stator/stator_mod_e_bl'//trim(numchar(nob(iproc)))//'.x')
     endif

     ! bypass of comm grid
     !igrd=6

     ! already dimensional but in mm
     chord=41.16_wp ! mm
     xgc=xgc*1e-3_wp
     ygc=ygc*1e-3_wp

     Lxp=0.0_wp
     Lyp=0.0_wp

     ! Tam & Dong radiation center (default) TO BE CHANGED
     xcr_=0.0_wp
     ycr_=0.0_wp

  ! Idealized blade vane configuration
  ! ==================================
  elseif (nbloc==32) then
     ! RANS
     ! ====
     if (is_adjoint_block) then
        ! call grid_read('Grid_RANS_Passmann/idealized_blade_Passmann_bl'//trim(numchar(nob(iproc)))//'.x')
        call mpistop('Problem with adjoint block',0)
     else
        ! call grid_read('Grid_RANS/idealized_blade_mod_bl'//trim(numchar(nob(iproc)))//'.x')
        ! call grid_read('Grid_RANS_Passmann/idealized_blade_Passmann_mod_bl'//trim(numchar(nob(iproc)))//'.x')
        ! call grid_read('Grid_Passm_1/grid1_mod_bl'//trim(numchar(nob(iproc)))//'.x')
        ! call grid_read('Grid_Passm_c1/grid1_mod_bl'//trim(numchar(nob(iproc)))//'.x')
        ! call grid_read('Grid_Passm_c2/grid1_mod_bl'//trim(numchar(nob(iproc)))//'.x')
        ! call grid_read('Grid_RANS_DLES_vfinal/grid1_mod_bl'//trim(numchar(nob(iproc)))//'.x')
        ! call grid_read('Grid_RANS_ORC23_0deg/grid1_mod_bl'//trim(numchar(nob(iproc)))//'.x')
        ! call grid_read('Grid_RANS_ORC23_0deg_2/grid1_mod_bl'//trim(numchar(nob(iproc)))//'.x')
        call grid_read('Grid_RANS_ORC23_0deg_3/grid1_mod_bl'//trim(numchar(nob(iproc)))//'.x')
        ! call grid_read('Grid_RANS_Coarse_ORC23_0deg_3/grid1_mod_bl'//trim(numchar(nob(iproc)))//'.x')
        ! call grid_read('Grid_RANS_VeryCoarse_ORC23_0deg_3/grid1_mod_bl'//trim(numchar(nob(iproc)))//'.x')
        ! call grid_read('Grid_RANS_SlightlyCoarse_ORC23_0deg_3/grid1_mod_bl'//trim(numchar(nob(iproc)))//'.x')
        ! call grid_read('Grid_Passm_2/grid2_mod_bl'//trim(numchar(nob(iproc)))//'.x')
        ! call grid_read('Grid_Passm_3/grid3_mod_bl'//trim(numchar(nob(iproc)))//'.x')
        ! call grid_read('Grid_Passm_1_ex/grid1_mod_bl'//trim(numchar(nob(iproc)))//'.x')
        ! call grid_read('Grid_Passm_1_ex2/grid1_mod_bl'//trim(numchar(nob(iproc)))//'.x')
        ! call grid_read('Grid_Passm_debeuger/debeug_mod_bl'//trim(numchar(nob(iproc)))//'.x')
        ! call grid_read('Grid_Passm_debeuger2/debeug_mod_bl'//trim(numchar(nob(iproc)))//'.x')
        ! call grid_read('Grid_Passm_debeuger3/debeug_mod_bl'//trim(numchar(nob(iproc)))//'.x')
        ! call mpistop('Problem without adjoint block',0)
     endif

     chord=L_ref

     xgc=xgc*chord
     ygc=ygc*chord

     Lxp=0.0_wp
     Lyp=0.0_wp

     ! Tam & Dong radiation center (default) TO BE CHANGED
     xcr_=0.0_wp
     ycr_=0.0_wp
  endif

  ! simple grid along z
  ! ===================
  dz0=deltaz
  if (is_2D) then
     zg=0.0_wp
  else
     zg(ngz/2)=0.0_wp
     do k=ngz/2+1,ngz
        zg(k)=zg(k-1)+dz0
     enddo
     do k=ngz/2-1,1,-1
        zg(k)=zg(k+1)-dz0
     enddo
  endif

!!$  if (is_curv3) then
!!$     do i=1,ngx
!!$        do j=1,ngy
!!$           do k=1,ngz
!!$              zgc3(i,j,k)=zg(k)
!!$           enddo
!!$        enddo
!!$     enddo
!!$
!!$     open(50,file='Grid_RANS_GO/ls593_mod_bl'//trim(numchar(nob(iproc)))//'.x',form='formatted')
!!$     write(50,*) 1
!!$     write(50,*) ngx,ngy,ngz
!!$     write(50,*) (((xgc3(i,j,k),i=1,ngx),j=1,ngy),k=1,ngz)
!!$     write(50,*) (((ygc3(i,j,k),i=1,ngx),j=1,ngy),k=1,ngz)
!!$     write(50,*) (((zgc3(i,j,k),i=1,ngx),j=1,ngy),k=1,ngz)
!!$  endif

  ! Domain sizes
  ! ============
  Lx = 1.
  Ly = 1.
  Lz = 1.

end subroutine grid_turb

!===============================================================================
subroutine grid_act
!===============================================================================
  !> read the actuator grid from the pointwise files
!===============================================================================
  use mod_block
  use mod_flow
  use mod_constant
  use mod_mpi
  use mod_utils
  implicit none
  ! ----------------------------------------------------------------------------
  integer :: k
  real(wp) :: dz
  ! ----------------------------------------------------------------------------

     call grid_read('grid/grid_bl'//trim(numchar(nob(iproc)))//'.x')

     if (.not.is_2d) then
        ! extruding 2D grid along z, uniformly (periodic span)
        ! ====================================================
        ! Note that nozzle exit throat is 0.2 mm
        dz=0.015_wp ! in mm. It corresponds to Delta_z^+ ~ 15
        zg(ngz/2)=0.0_wp
        do k=ngz/2+1,ngz
           zg(k)=zg(k-1)+dz
        enddo
        do k=ngz/2-1,1,-1
           zg(k)=zg(k+1)-dz
        enddo
     endif

     ! rescale the dimensions (mm to m)
     xgc=xgc/1000.0_wp
     ygc=ygc/1000.0_wp
     zg =zg /1000.0_wp

     ! Tam & Dong radiation center (in meters)
     xcr_=0.0_wp
!      ycr_=0.03499258199729956_wp
     if(iproc==iproc_leader(19)) ycr_=ygc(nx/2,1)
     call MPI_BCAST(ycr_,1,MPI_DOUBLE_PRECISION,iproc_leader(19),COMM_global,info)

end subroutine grid_act

!===============================================================================
subroutine grid_flatplate
!===============================================================================
  !> read the actuator grid from the pointwise files
!===============================================================================
  use mod_block
  use mod_flow
  use mod_constant
  use mod_mpi
  use mod_utils
  implicit none
  ! ----------------------------------------------------------------------------
  integer :: k
  real(wp) :: dz
  ! ----------------------------------------------------------------------------

     call grid_read('grid/grid1_bl'//trim(numchar(nob(iproc)))//'.x')

     if (.not.is_2d) then
        ! extruding 2D grid along z, uniformly (periodic span)
        ! ====================================================
        dz=0.007843137_wp
        zg(ngz/2)=0.0_wp
        do k=ngz/2+1,ngz
           zg(k)=zg(k-1)+dz
        enddo
        do k=ngz/2-1,1,-1
           zg(k)=zg(k+1)-dz
        enddo
     endif

     ! Tam & Dong radiation center
     xcr_=1.0_wp
     ycr_=0.0_wp

end subroutine grid_flatplate

!===============================================================================
subroutine grid_backstep
!===============================================================================
  !> read the actuator grid from the pointwise files
!===============================================================================
  use mod_block
  use mod_flow
  use mod_constant
  use mod_mpi
  use mod_utils
  implicit none
  ! ----------------------------------------------------------------------------
  integer :: k
  real(wp) :: dz
  ! ----------------------------------------------------------------------------

     call grid_read('grid/grid_3levdn_bl'//trim(numchar(nob(iproc)))//'.x')

     if (.not.is_2d) then
        ! extruding 2D grid along z, uniformly (periodic span)
        ! ====================================================
        dz=0.08_wp
        zg(ngz/2)=0.0_wp
        do k=ngz/2+1,ngz
           zg(k)=zg(k-1)+dz
        enddo
        do k=ngz/2-1,1,-1
           zg(k)=zg(k+1)-dz
        enddo
     endif

     ! Tam & Dong radiation center
     xcr_=-110.0_wp
     ycr_=5.0_wp

end subroutine grid_backstep

!===============================================================================
subroutine grid_shit
!===============================================================================
  !> compute predefined grid for spatial HIT cases
!===============================================================================
  use warnstop
  use mod_flow
  use mod_constant
  use mod_block
  use mod_mpi_part
  use mod_mpi
  use mod_bc_periodicity ! <~ needed for Lxp,Lyp,Lzp  TO BE CHANGED
  use mod_utils
  implicit none
  ! ----------------------------------------------------------------------------
  integer :: i,j,k
  real(wp) :: dx0,dy0,dz0,dx0_old,chord
  ! ----------------------------------------------------------------------------

  if (nbloc==4) then

     if (is_adjoint_block) then
        call grid_read('inlet_bl'//trim(numchar(nob(iproc)))//'.x')
     else
        call grid_read('inletc_mod_bl'//trim(numchar(nob(iproc)))//'.x')
     endif

     chord=60.0e-3

     xgc=xgc*chord
     ygc=ygc*chord

     Lxp=0.0_wp*chord
     Lyp=0.709666666666667_wp*chord

     ! Tam & Dong radiation center (default) TO BE CHANGED
     xcr_=0.0_wp
     ycr_=0.0_wp
     xcr_=-0.015
     ycr_= 0.005

  else

     ! Vortex
     ! ======
     ! dx0=deltax
     ! xg(1)=0.0_wp
     ! do i=2,ngx
     !    xg(i)=xg(i-1)+dx0
     ! enddo

     ! dy0=deltay
     ! yg(ngy/2) = 0.0_wp
     ! do j=ngy/2+1,ngy
     !    yg(j)=yg(j-1)+dy0
     ! enddo
     ! do j=ngy/2-1,1,-1
     !    yg(j)=yg(j+1)-dy0
     ! enddo

     ! dz0=deltaz
     ! zg(ngz/2) = 0.0_wp
     ! do k=ngz/2+1,ngz
     !    zg(k)=zg(k-1)+dz0
     ! enddo
     ! do k=ngz/2-1,1,-1
     !    zg(k)=zg(k+1)-dz0
     ! enddo

     ! Spatial HIT
     ! ===========
     if (.not.is_curv) then
        dx0=deltax
        xg(1)=0.0_wp
        do i=2,ngx-50
           xg(i)=xg(i-1)+dx0
        enddo

        do i=ngx-49,ngx-20
           xg(i)=xg(i-1)+dx0
           dx0=1.02*dx0
        enddo

        do i=ngx-19,ngx
           xg(i)=xg(i-1)+dx0
           dx0=1.06*dx0
        enddo

        dy0=deltay
        yg(ngy/2) = 0.0_wp
        do j=ngy/2+1,ngy
           yg(j)=yg(j-1)+dy0
        enddo
        do j=ngy/2-1,1,-1
           yg(j)=yg(j+1)-dy0
        enddo

        xcr_ = xg(10)
        ycr_ = 0.0_wp

     else
        dx0=deltax
        xgc(1,:)=0.0_wp
        do j=1,ngy
           do i=2,ngx-50
              xgc(i,j)=xgc(i-1,j)+dx0
           enddo
        enddo

        dx0_old = dx0
        do j=1,ngy
           dx0 = dx0_old
           do i=ngx-49,ngx-20
              xgc(i,j)=xgc(i-1,j)+dx0
              dx0=1.02_wp*dx0
           enddo
        enddo
        dx0_old = dx0

        do j=1,ngy
           dx0 = dx0_old
           do i=ngx-19,ngx
              xgc(i,j)=xgc(i-1,j)+dx0
              dx0=1.06*dx0
           enddo
        enddo

        dy0=deltay
        ygc(:,ngy/2) = 0.0_wp
        do j=ngy/2+1,ngy
           do i=1,ngx
              ygc(i,j)=ygc(i,j-1)+dy0
           enddo
        enddo
        do j=ngy/2-1,1,-1
           do i=1,ngx
              ygc(i,j)=ygc(i,j+1)-dy0
           enddo
        enddo

        xcr_ = xgc(10,1)
        ycr_ = 0.0_wp
        dxmax = max(deltax,deltay)
     endif
  endif


  dz0=deltaz
  if (is_2D) then
     zg=0.0_wp
  else
     zg(ngz/2)=0.0_wp
     do k=ngz/2+1,ngz
        zg(k)=zg(k-1)+dz0
     enddo
     do k=ngz/2-1,1,-1
        zg(k)=zg(k+1)-dz0
     enddo
  endif

end subroutine grid_shit


!===============================================================================
subroutine grid_le
!===============================================================================
  !> Read grid for the vane configuration
!===============================================================================
  use mod_flow
  use mod_block
  use mod_constant
  use mod_mpi
  use mod_utils
  use mod_mpi_part
  use warnstop
  implicit none
  ! ----------------------------------------------------------------------------
  integer :: jg,k
  ! integer :: i,j,ibl2read,ngx_,ngy_,ngz_
  real(wp) :: dz,xmin_loc,xmax_loc,dxmax_loc
  ! real(wp), dimension(:,:), allocatable :: xgc_temp,ygc_temp
  ! ----------------------------------------------------------------------------

   ! is_curv=.true.

  ! if (idepart==POST_PROCESSING) then
  !    if (iblc_pp.ne.0) then
  !       ibl2read=iblc_pp
  !       if (type_pp.ne.7) then
  !          idem_pp = 1; iend_pp = bl(iblc_pp)%ni
  !          jdem_pp = 1; jend_pp = bl(iblc_pp)%nj
  !       endif
  !    else
  !       ibl2read=nob(iproc)
  !    endif

  !    ! open(194,file=trim(dirDATA)//'grid_bl'//trim(numchar(ibl2read))//"_proc"//trim(numchar(iproc))//'.bin', &
  !    !   form='unformatted',status='unknown')
  !    open(194,file=trim(dirDATA)//'grid_bl'//trim(numchar(ibl2read))//'.bin', &
  !      form='unformatted',status='unknown')
  !    rewind(194)
  !    read(194) ngx_
  !    read(194) ngy_
  !    read(194) ngz_
  !    allocate(xgc_temp(1:ngx_,1:ngy_))
  !    allocate(ygc_temp(1:ngx_,1:ngy_))
  !    read(194) ((xgc_temp(i,j),i=1,ngx_),j=1,ngy_)
  !    read(194) ((ygc_temp(i,j),i=1,ngx_),j=1,ngy_)
  !    read(194) (zg(k),k=1,ngz_)
  !    close(194)

  !    do i=idem_pp,iend_pp
  !       do j=jdem_pp,jend_pp
  !          xgc(i+1-idem_pp,j+1-jdem_pp) = xgc_temp(i,j)
  !          ygc(i+1-idem_pp,j+1-jdem_pp) = ygc_temp(i,j)
  !       enddo
  !    enddo
  !    deallocate(xgc_temp,ygc_temp)

  !    ! if (iproc.eq.7) print *,"xgc",xgc()

  !    return
  ! endif

  ! Read grid
  ! ---------
  if (is_adjoint_block) then
     call mpistop('is_adjoint_block is set to true !',0)
     ! call grid_read('Grid/FSTT_LE_bl'//trim(numchar(nob(iproc)))//'.x')
  else
     call grid_read('Grid/FSTT_LE_mod_bl'//trim(numchar(nob(iproc)))//'.x')
  endif

  ! Dimensionalize grid
  ! -------------------
  xgc=xgc*L_ref
  ygc=ygc*L_ref

  ! simple grid along z
  ! ===================
  dz=deltaz
  if (is_2D) then
     zg=0.0_wp
  else
     zg(ngz/2)=0.0_wp
     do k=ngz/2+1,ngz
        zg(k)=zg(k-1)+dz
     enddo
     do k=ngz/2-1,1,-1
        zg(k)=zg(k+1)-dz
     enddo
  endif

  ! Tam & Dong radiation center  TO BE CHANGED
  ! xcr_=0.0_wp
  xcr_ = xgc(10,ngy/2)
  ycr_=0.0_wp

  ! Definition of limit of domain
  xmin_loc = MINVAL(xgc)
  call MPI_Allreduce(xmin_loc,xmin,1,MPI_DOUBLE_PRECISION,MPI_MIN,COMM_global,info)
  xmax_loc = MAXVAL(xgc)
  call MPI_Allreduce(xmax_loc,xmax,1,MPI_DOUBLE_PRECISION,MPI_MIN,COMM_global,info)

  ! Determination of grid resolution at inlet plane
  dxmax_loc = -1e6;
  ! if inlet
  if (bl(nob(iproc))%BC(1).lt.0) then
     if (.not.is_2D) dxmax_loc = deltaz
     do jg=1,ngy
        dxmax_loc = max(abs(xgc(2,jg)-xgc(1,jg)),dxmax_loc)
     enddo
  endif
  call MPI_Allreduce(dxmax_loc,dxmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,COMM_global,info)

  if (iproc.eq.0) print *,"grid step maximal at inlet:",dxmax
end subroutine grid_le

!===============================================================================
subroutine grid_TE
!===============================================================================
  !> compute predefined grid for turbine cases
!===============================================================================
  use mod_block
  use mod_flow
  use mod_constant
  use mod_mpi_part
  use mod_bc_periodicity ! <~ needed for Lxp,Lyp,Lzp  TO BE CHANGED
  use mod_utils
  use warnstop
  implicit none
  ! ----------------------------------------------------------------------------
  integer :: k
  real(wp) :: dz0,chord
  ! ----------------------------------------------------------------------------
  real(wp) :: Lx,Ly,Lz
  ! ----------------------------------------------------------------------------

  ! TE configuration
  ! ================
  if (nbloc==8) then
     if (is_adjoint_block) then
        ! call grid_read('Grid_RANS_Passmann/idealized_blade_Passmann_bl'//trim(numchar(nob(iproc)))//'.x')
        call mpistop('Problem with adjoint block',0)
     else
        ! call grid_read('blade_0deg_TE_DDES/TE_DDES_mod_bl'//trim(numchar(nob(iproc)))//'.x')
        call grid_read('blade_0deg_TE_DDES_vCoarse/TE_DDES_mod_bl'//trim(numchar(nob(iproc)))//'.x')
        ! call mpistop('Problem without adjoint block',0)
     endif

     chord=L_ref

     xgc=xgc*chord
     ygc=ygc*chord

     Lxp=0.0_wp
     Lyp=0.0_wp

     ! Tam & Dong radiation center (default) TO BE CHANGED
     xcr_=0.0_wp
     ycr_=0.0_wp
  else
     call mpistop('Grid in TE configuration not managed...',0)
  endif

  ! simple grid along z
  ! ===================
  dz0=deltaz
  if (is_2D) then
     zg=0.0_wp
  else
     zg(ngz/2)=0.0_wp
     do k=ngz/2+1,ngz
        zg(k)=zg(k-1)+dz0
     enddo
     do k=ngz/2-1,1,-1
        zg(k)=zg(k+1)-dz0
     enddo
  endif

  ! Domain sizes
  ! ============
  Lx = 1.
  Ly = 1.
  Lz = 1.

end subroutine grid_TE
