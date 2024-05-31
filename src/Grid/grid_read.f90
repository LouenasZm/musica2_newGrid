!==============================================================================
subroutine grid_read(gridname)
!==============================================================================
  !> Read grids in PLOT3D format (ASCII)
!==============================================================================
  use mod_grid
  use mod_mpi
  implicit none
  ! ---------------------------------------------------------------------------
  character(len=*) :: gridname
  ! ---------------------------------------------------------------------------
  ! Local variables
  integer :: i,j,k,m
  integer :: ni,nj,nk
  ! ---------------------------------------------------------------------------

  open(50,file=gridname,form='formatted')
  rewind(50)

  if (is_curv3) then ! read a 3-D curvilinear grid file
     read(50,*) m
     read(50,*) ni,nj,nk
     !print *,'dim',ni,nj,nk,' for block',nob(iproc)
     read(50,*) (((xgc3(i,j,k),i=1,ni),j=1,nj),k=1,nk)
     read(50,*) (((ygc3(i,j,k),i=1,ni),j=1,nj),k=1,nk)
     read(50,*) (((zgc3(i,j,k),i=1,ni),j=1,nj),k=1,nk)
     
!!$     read(50,*) ((xgc(i,j),i=1,ni),j=1,nj)
!!$     read(50,*) ((ygc(i,j),i=1,ni),j=1,nj)
!!$     do k=1,ngz
!!$        xgc3(:,:,k)=xgc
!!$        ygc3(:,:,k)=ygc
!!$     enddo
  else
     if (is_curv) then ! read a 2-D curvilinear grid file
        read(50,*) m
        read(50,*) ni,nj,nk
        !print *,'dim',ni,nj,' for block',nob(iproc)
        read(50,*) ((xgc(i,j),i=1,ni),j=1,nj)
        read(50,*) ((ygc(i,j),i=1,ni),j=1,nj)
     else ! read a Cartesian grid file
        read(50,*) m
        read(50,*) ni,nj,nk
        !print *,'dim',ni,nj,' for block',nob(iproc)
        read(50,*) ((xg(i),i=1,ni),j=1,nj)
        read(50,*) ((yg(j),i=1,ni),j=1,nj)
     endif
  endif

  close(50)

end subroutine grid_read

!==============================================================================
subroutine grid_read_ex(gridname)
!==============================================================================
  !> Read grids in PLOT3D format (ASCII) - Extended grid
!==============================================================================
  use mod_grid
  use mod_block
  use mod_mpi
  implicit none
  ! ---------------------------------------------------------------------------
  character(len=*) :: gridname
  ! ---------------------------------------------------------------------------
  ! Local variables
  integer :: i,j,k,n
  integer :: ni,nj,nk
  integer :: ni1,nj1,nk1,ni2,nj2,nk2
  ! ---------------------------------------------------------------------------
!!$  real(wp), dimension(:,:,:), allocatable :: xgc3n,ygc3n,zgc3n

  open(50,file=gridname,form='formatted')
  rewind(50)

  if (is_curv3) then ! read a 3-D curvilinear grid file
     read(50,*) n
     read(50,*) ni,nj,nk

     n=nob(iproc)
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

     read(50,*) (((xgc3(i,j,k),i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
     read(50,*) (((ygc3(i,j,k),i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
     read(50,*) (((zgc3(i,j,k),i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)

!!$     n=nob(iproc)
!!$     ni1=1
!!$     ni2=ngx
!!$     nj1=1
!!$     nj2=ngy
!!$     nk1=1
!!$     nk2=90
!!$     if (bl(n)%BC(1)>0) ni1=1-ngh
!!$     if (bl(n)%BC(2)>0) ni2=ngx+ngh
!!$     if (bl(n)%BC(3)>0) nj1=1-ngh
!!$     if (bl(n)%BC(4)>0) nj2=ngy+ngh
!!$     if (bl(n)%BC(5)>0) nk1=1-ngh
!!$     if (bl(n)%BC(6)>0) nk2=ngz+ngh
!!$
!!$     allocate(xgc3n(1-ngh:ngx+ngh,1-ngh:ngy+ngh,1-ngh:90+ngh))
!!$     allocate(ygc3n(1-ngh:ngx+ngh,1-ngh:ngy+ngh,1-ngh:90+ngh))
!!$     allocate(zgc3n(1-ngh:ngx+ngh,1-ngh:ngy+ngh,1-ngh:90+ngh))
!!$
!!$     read(50,*) (((xgc3n(i,j,k),i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
!!$     read(50,*) (((ygc3n(i,j,k),i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
!!$     read(50,*) (((zgc3n(i,j,k),i=ni1,ni2),j=nj1,nj2),k=nk1,nk2)
!!$     
!!$     do k=1-ngh,ngz+ngh
!!$        do j=1-ngh,ngy+ngh
!!$           do i=1-ngh,ngx+ngh
!!$              xgc3(i,j,k)=xgc3n(i,j,k)
!!$              ygc3(i,j,k)=ygc3n(i,j,k)
!!$              zgc3(i,j,k)=zgc3n(i,j,k)
!!$           enddo
!!$        enddo
!!$     enddo
!!$
!!$     deallocate(xgc3n,ygc3n,zgc3n)

  endif

  close(50)
  
end subroutine grid_read_ex

!==============================================================================
subroutine grid_read32(gridname)
!==============================================================================
  !> Read grids in PLOT3D format (ASCII)
!==============================================================================
  use mod_grid
  use mod_mpi
  implicit none
  ! ---------------------------------------------------------------------------
  character(len=*) :: gridname
  ! ---------------------------------------------------------------------------
  ! Local variables
  integer :: i,j,k,m
  integer :: ni,nj,nk
  ! ---------------------------------------------------------------------------

  open(50,file=gridname,form='formatted')

  print *,'sals'
  allocate(xgc3(1-ngh:ngx+ngh,1-ngh:ngy+ngh,1-ngh:ngz+ngh))
  allocate(ygc3(1-ngh:ngx+ngh,1-ngh:ngy+ngh,1-ngh:ngz+ngh))
  allocate(zgc3(1-ngh:ngx+ngh,1-ngh:ngy+ngh,1-ngh:ngz+ngh))
  
  read(50,*) m
  read(50,*) ni,nj,nk
  !print *,'dim',ni,nj,nk,' for block',nob(iproc)
  read(50,*) (((xgc3(i,j,k),i=1,ni),j=1,nj),k=1,nk)
  read(50,*) (((ygc3(i,j,k),i=1,ni),j=1,nj),k=1,nk)
  read(50,*) (((zgc3(i,j,k),i=1,ni),j=1,nj),k=1,nk)

  xgc=xgc3(:,:,1)
  ygc=ygc3(:,:,1)
  zg=zgc3(1,1,:)

  deallocate(xgc3,ygc3,zgc3)
  
end subroutine grid_read32
