!===============================================================================
subroutine grid_finalize
!===============================================================================
  !> Finalize the grid elements
  !> 1. determine grid center (i0,j0,k0) per block
  !> 2. write grid
  !> 3. communicate global grid among procs
  !> 4. partitioning of the grid on local proc
  !> 5. compute grid metrics
  !> 6. communicate metrics
  !> 7. compute wall normals
!===============================================================================
  use mod_mpi_part
  use mod_grid
  use mod_constant
  use warnstop
  use mod_utils
  use mod_grid_directions ! TO BE CHANGED used to determine deltay -> must be suppressed
  implicit none
  ! ----------------------------------------------------------------------------
  integer :: i,j,k
  ! ----------------------------------------------------------------------------

  ! Determine grid center (i0,j0,k0) per domain
  ! ===========================================
  if (mod(nx,2)==0) then
     i0=nx/2
  else
     i0=(nx+1)/2
  endif
  if (mod(ny,2)==0) then
     j0=ny/2
  else
     j0=(ny+1)/2
  endif
  if (mod(nz,2)==0) then
     k0=nz/2
  else
     k0=(nz+1)/2
  endif
  
  ! Write global grid per block (fortran IEEE binary file - Att big or little endian)
  ! ===========================
  if (.not.is_read_ex) then

     if (iproc.eq.iproc_leader(nob(iproc))) then
        open(194,file='grid_bl'//trim(numchar(nob(iproc)))//'.bin',form='unformatted',status='unknown')
        rewind(194)
        write(194) ngx
        write(194) ngy
        write(194) ngz
        if (is_curv3) then
           write(194) (((xgc3(i,j,k),i=1,ngx),j=1,ngy),k=1,ngz)
           write(194) (((ygc3(i,j,k),i=1,ngx),j=1,ngy),k=1,ngz)
           write(194) (((zgc3(i,j,k),i=1,ngx),j=1,ngy),k=1,ngz)
        else
           if (is_curv) then
              write(194) ((xgc(i,j),i=1,ngx),j=1,ngy)
              write(194) ((ygc(i,j),i=1,ngx),j=1,ngy)
           else
              write(194) (xg(i),i=1,ngx)
              write(194) (yg(j),j=1,ngy)
           endif
           write(194) (zg(k),k=1,ngz)
        endif

        ! write swap/rev indicators to fill ghost points in post-pro
        ! (added at the end for compatibility with old files)
        do i=1,4
           if (is_swapij2_bl(i)) then
              write(194) 1
           else
              write(194) 0
           endif
           if (is_rev2_bl(i,1)) then
              write(194) 1
           else
              write(194) 0
           endif
           if (is_rev2_bl(i,2)) then
              write(194) 1
           else
              write(194) 0
           endif
        enddo
        close(194)
     endif

     if (iproc==0) print *,'Global grid OK'

     ! Determine global xmin,xmax,ymin,ymax,zmin,zmax
     ! ==============================================
     if (iproc.eq.iproc_leader(nob(iproc))) then

        if (is_curv3) then
           xmin=minval(minval(minval(xgc3(1:ngx,1:ngy,1:ngz),1),1),1)
           xmax=maxval(maxval(maxval(xgc3(1:ngx,1:ngy,1:ngz),1),1),1)
           ymin=minval(minval(minval(ygc3(1:ngx,1:ngy,1:ngz),1),1),1)
           ymax=maxval(maxval(maxval(ygc3(1:ngx,1:ngy,1:ngz),1),1),1)
           zmin=minval(minval(minval(zgc3(1:ngx,1:ngy,1:ngz),1),1),1)
           zmax=maxval(maxval(maxval(zgc3(1:ngx,1:ngy,1:ngz),1),1),1)
        else
           if (is_curv) then
              xmin=minval(minval(xgc(1:ngx,1:ngy),1),1)
              xmax=maxval(maxval(xgc(1:ngx,1:ngy),1),1)
              ymin=minval(minval(ygc(1:ngx,1:ngy),1),1)
              ymax=maxval(maxval(ygc(1:ngx,1:ngy),1),1)
           else
              xmin=minval(xg(1:ngx),1)
              xmax=maxval(xg(1:ngx),1)
              ymin=minval(yg(1:ngy),1)
              ymax=maxval(yg(1:ngy),1)
           endif
           zmin=minval(zg(1:ngz),1)
           zmax=maxval(zg(1:ngz),1)
        endif

        call MPI_ALLREDUCE(MPI_IN_PLACE,xmin,1,MPI_DOUBLE_PRECISION,MPI_MIN,COMM_interblock,info)
        call MPI_ALLREDUCE(MPI_IN_PLACE,xmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,COMM_interblock,info)
        call MPI_ALLREDUCE(MPI_IN_PLACE,ymin,1,MPI_DOUBLE_PRECISION,MPI_MIN,COMM_interblock,info)
        call MPI_ALLREDUCE(MPI_IN_PLACE,ymax,1,MPI_DOUBLE_PRECISION,MPI_MAX,COMM_interblock,info)
        call MPI_ALLREDUCE(MPI_IN_PLACE,zmin,1,MPI_DOUBLE_PRECISION,MPI_MIN,COMM_interblock,info)
        call MPI_ALLREDUCE(MPI_IN_PLACE,zmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,COMM_interblock,info)
     endif
  else
     ! Determine local xmin,xmax,ymin,ymax,zmin,zmax
     ! ==============================================
     if (is_curv3) then
        xmin=minval(minval(minval(xc3(1:nx,1:ny,1:nz),1),1),1)
        xmax=maxval(maxval(maxval(xc3(1:nx,1:ny,1:nz),1),1),1)
        ymin=minval(minval(minval(yc3(1:nx,1:ny,1:nz),1),1),1)
        ymax=maxval(maxval(maxval(yc3(1:nx,1:ny,1:nz),1),1),1)
        zmin=minval(minval(minval(zc3(1:nx,1:ny,1:nz),1),1),1)
        zmax=maxval(maxval(maxval(zc3(1:nx,1:ny,1:nz),1),1),1)
     else
        if (is_curv) then
           xmin=minval(minval(xc(1:nx,1:ny),1),1)
           xmax=maxval(maxval(xc(1:nx,1:ny),1),1)
           ymin=minval(minval(yc(1:nx,1:ny),1),1)
           ymax=maxval(maxval(yc(1:nx,1:ny),1),1)
        else
           xmin=minval(x(1:nx),1)
           xmax=maxval(x(1:nx),1)
           ymin=minval(y(1:ny),1)
           ymax=maxval(y(1:ny),1)
        endif
        zmin=minval(z(1:nz),1)
        zmax=maxval(z(1:nz),1)
     endif

     call MPI_ALLREDUCE(MPI_IN_PLACE,xmin,1,MPI_DOUBLE_PRECISION,MPI_MIN,COMM_global,info)
     call MPI_ALLREDUCE(MPI_IN_PLACE,xmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,COMM_global,info)
     call MPI_ALLREDUCE(MPI_IN_PLACE,ymin,1,MPI_DOUBLE_PRECISION,MPI_MIN,COMM_global,info)
     call MPI_ALLREDUCE(MPI_IN_PLACE,ymax,1,MPI_DOUBLE_PRECISION,MPI_MAX,COMM_global,info)
     call MPI_ALLREDUCE(MPI_IN_PLACE,zmin,1,MPI_DOUBLE_PRECISION,MPI_MIN,COMM_global,info)
     call MPI_ALLREDUCE(MPI_IN_PLACE,zmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,COMM_global,info)

  endif
     
   ! ! Determine global Lx,Ly,Lz <~ Necessary for some test cases (done in grid.f90 before, which is now grid_old.f90)
   ! ! =========================
   ! Lx=xmax-xmin
   ! Ly=ymax-ymin
   ! Lz=zmax-zmin
   ! if (iproc==0) print *,'after',Lx,Ly,Lz
   !if (iproc==0) print *,'after',Lx+deltax,Ly+deltay,Lz+deltaz
   ! if (iproc==0) print *,'after',Lx+xgc(2,ngy)-xgc(1,ngy),Ly,Lz

  ! Interblock and Intrablock communications for grids
  ! ==================================================
  if ((igrd.ne.6).and.(is_grid_old)) call grid_comm
  ! for SRC and igrd=6 (skewed and sinusoidal grids) : pb with double periodicity
  ! in general, igrd=6 is used to bypass grid_comm

  if (.not.is_read_ex) then
     if (iproc.eq.iproc_leader(nob(iproc))) then
        open(194,file='grid_bl'//trim(numchar(nob(iproc)))//'_ex.bin',form='unformatted',status='unknown')
        rewind(194)
        write(194) ngx+2*ngh
        write(194) ngy+2*ngh
        write(194) ngz+2*ngh
        if (is_curv3) then
           write(194) (((xgc3(i,j,k),i=1-ngh,ngx+ngh),j=1-ngh,ngy+ngh),k=1-ngh,ngz+ngh)
           write(194) (((ygc3(i,j,k),i=1-ngh,ngx+ngh),j=1-ngh,ngy+ngh),k=1-ngh,ngz+ngh)
           write(194) (((zgc3(i,j,k),i=1-ngh,ngx+ngh),j=1-ngh,ngy+ngh),k=1-ngh,ngz+ngh)
        else
           if (is_curv) then
              write(194) ((xgc(i,j),i=1-ngh,ngx+ngh),j=1-ngh,ngy+ngh)
              write(194) ((ygc(i,j),i=1-ngh,ngx+ngh),j=1-ngh,ngy+ngh)
           else
              write(194) (xg(i),i=1-ngh,ngx+ngh)
              write(194) (yg(j),j=1-ngh,ngy+ngh)
           endif
        endif
        close(194)
     endif

     if (iproc==0) print *,'Comm grid OK'
     !call mpistop('stop here',0)

     ! Partitioning of the grid
     ! ========================
     call grid_local
  endif

  if (iproc==0) print *,'Local grid OK'

  ! Compute metrics
  ! ===============
  call grid_metrics
  
  if (iproc==0) print *,'Metrics grid OK'

  ! Check metrics commutations
  ! ==========================
!!$  if (is_curv3) then
!!$     call grid_metrics3d_commutation
!!$  endif

!!$  if (is_curv) then
!!$     select case (stencil)
!!$        case(11)
!!$           call grid_metrics_commutation_11pts
!!$        case(5)
!!$           !call grid_metrics_commutation_5pts
!!$        case default
!!$           call mpistop('Only 11pts and 5pts metrics commutations yet.',0)
!!$     end select
!!$  endif
  
  if ((iproc==0).and.(is_curv)) print *,'check commutations OK'
  
  ! Enforce metrics
  if ((TGV).or.(CHIT)) then
     idx= 1.0_wp/deltax;   idx_v= 1.0_wp/deltax
     idy= 1.0_wp/deltay;   idy_v= 1.0_wp/deltay
     idz= 1.0_wp/deltaz;   idz_v= 1.0_wp/deltaz
  endif
  if ((CHAN).or.(PHILL)) then
     idz  = 1.0_wp/deltaz
     idz_v= 1.0_wp/deltaz
  endif

  ! Computation of wall normals and BC parameters
  ! =============================================
  if ((is_curv).or.(is_curv3)) then
     call grid_normals
     if (iproc==0) print *,'Grid normals OK'
  endif
  
  if (is_curv) then
     if ((.not.CYL).and.(.not.TURB).and.(.not.TE).and.(.not.T3C)) &
          deltay=abs(ygc(1,2)-ygc(1,1))
     if (is_swapij_bl(3).or.is_swapij_bl(4)) &
          ! deltay=abs(ygc(2,1)-ygc(1,1))  <~ Useful ?? What is it for ?
          deltay=abs(yc(2,1)-yc(1,1))
     if (ACT) deltay = 1d-6!4d-7
  else
     !deltay=1.0_wp/dyg(1)
     ! deltay=yg(2)-yg(1)  <~ Useful ?? What is it for ?
     deltay=y(2)-y(1)
  endif
  ! added by XG 04/03/22 for post-processing old channel flow databases
!!$  if (CHAN) then
!!$     if (coord(2)==0) deltay=1.0_wp/idy(1)
!!$  endif
  
end subroutine grid_finalize
 
!===============================================================================
subroutine grid_local
!===============================================================================
  !> Partitioning of the grid on MPI proc
!===============================================================================
  use mod_mpi
  use mod_flow
  use mod_constant
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j,k
  ! ---------------------------------------------------------------------------

  if (is_curv3) then

     ! Partitioning of the 3D curvilinear grid
     ! ---------------------------------------
     do k=1-ngh,nz+ngh
        do j=1-ngh,ny+ngh
           do i=1-ngh,nx+ngh
              xc3(i,j,k)=xgc3(i+coord(1)*nx,j+coord(2)*ny,k+coord(3)*nz)
              yc3(i,j,k)=ygc3(i+coord(1)*nx,j+coord(2)*ny,k+coord(3)*nz)
              zc3(i,j,k)=zgc3(i+coord(1)*nx,j+coord(2)*ny,k+coord(3)*nz)
           enddo
        enddo
     enddo

  else
     if (is_curv) then

        ! Partitioning of the curvilinear grid
        ! ------------------------------------
        do j=1-ngh,ny+ngh
           do i=1-ngh,nx+ngh
              xc(i,j)=xgc(i+coord(1)*nx,j+coord(2)*ny)
              yc(i,j)=ygc(i+coord(1)*nx,j+coord(2)*ny)
           enddo
        enddo

!!$        ! FOR COMPATIBILITY with is_curv3
!!$        do k=1-ngh,nz+ngh
!!$           do j=1-ngh,ny+ngh
!!$              do i=1-ngh,nx+ngh
!!$                 xc3(i,j,k)=xgc3(i+coord(1)*nx,j+coord(2)*ny,k+coord(3)*nz)
!!$                 yc3(i,j,k)=ygc3(i+coord(1)*nx,j+coord(2)*ny,k+coord(3)*nz)
!!$                 zc3(i,j,k)=zgc3(i+coord(1)*nx,j+coord(2)*ny,k+coord(3)*nz)
!!$              enddo
!!$           enddo
!!$        enddo

     else

        ! Partitioning of the Cartesian grid
        ! ----------------------------------
        do i=ndxt,nfxt
           x(i)=xg(i+coord(1)*nx)
        enddo

        do j=ndyt,nfyt
           y(j)=yg(j+coord(2)*ny)
        enddo

     end if

     ! Third direction is Cartesian
     ! ----------------------------
     if (is_2D) then
        z(1)=zg(1)
     else
        do k=ndzt,nfzt
           z(k)=zg(k+coord(3)*nz)
        enddo
     endif

  endif

end subroutine grid_local

!===============================================================================
subroutine grid_metrics
!===============================================================================
  !> Compute grid metrics
!===============================================================================
  use mod_grid
  use mod_bc
  use mod_constant ! <- for is_SBP
  use mod_grid_metrics_c3
  implicit none
  ! ---------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------
  
  if (is_curv3) then

     ! metrics with GCL
     ! ----------------
     call grid_metrics_3d

     call grid_metrics_gradients_3d
     
  else

     if (is_curv) then

        ! for inviscid fluxes (ngh ghost cells)
        ! -------------------------------------
        call grid_metrics_ksi
        call grid_metrics_eta

        ! for viscous fluxes (ngh_v ghost cells)
        ! --------------------------------------
        call grid_metrics_ksi_v
        call grid_metrics_eta_v

        ! curvilinear metrics for wall BC
        ! -------------------------------
        if (is_SBP) then
           if (BC_face(1,1)%sort==0) call grid_metrics_wall_imin_SBP4
           if (BC_face(1,2)%sort==0) call grid_metrics_wall_imax_SBP4
           if (BC_face(2,1)%sort==0) call grid_metrics_wall_jmin_SBP4
           if (BC_face(2,2)%sort==0) call grid_metrics_wall_jmax_SBP4
           if (BC_face(1,1)%sort<=-3) call grid_metrics_wall_imin_SBP4
           if (BC_face(1,2)%sort<=-3) call grid_metrics_wall_imax_SBP4
           if (BC_face(2,1)%sort<=-3) call grid_metrics_wall_jmin_SBP4
           if (BC_face(2,2)%sort<=-3) call grid_metrics_wall_jmax_SBP4
           if (BC_edge(1,1,1)%sort==2) call grid_metrics_wall_imin_jmin_SBP4
           if (BC_edge(1,1,2)%sort==2) call grid_metrics_wall_imin_jmax_SBP4
           if (BC_edge(1,2,1)%sort==2) call grid_metrics_wall_imax_jmin_SBP4
           if (BC_edge(1,2,2)%sort==2) call grid_metrics_wall_imax_jmax_SBP4
        else
           if (BC_face(1,1)%sort==0) call grid_metrics_wall_imin
           if (BC_face(1,2)%sort==0) call grid_metrics_wall_imax
           if (BC_face(2,1)%sort==0) call grid_metrics_wall_jmin
           if (BC_face(2,2)%sort==0) call grid_metrics_wall_jmax
           if (BC_face(1,1)%sort<=-3) call grid_metrics_wall_imin
           if (BC_face(1,2)%sort<=-3) call grid_metrics_wall_imax
           if (BC_face(2,1)%sort<=-3) call grid_metrics_wall_jmin
           if (BC_face(2,2)%sort<=-3) call grid_metrics_wall_jmax
           if (BC_edge(1,1,1)%sort==2) call grid_metrics_wall_imin_jmin
           if (BC_edge(1,1,2)%sort==2) call grid_metrics_wall_imin_jmax
           if (BC_edge(1,2,1)%sort==2) call grid_metrics_wall_imax_jmin
           if (BC_edge(1,2,2)%sort==2) call grid_metrics_wall_imax_jmax
        endif

        ! communications of metrics for ghost points
        ! ------------------------------------------
        call grid_comm_metrics_curv

        ! inverse Jacobian
        ! ----------------
        ! for inviscid fluxes (ngh ghost cells)
        call grid_metrics_ijacob
        ! for viscous fluxes (ngh_v ghost cells)
        call grid_metrics_ijacob_v

        ! norms of gradients metrics [for artificial viscosity] IRS TO BE CHANGED
        ! --------------------------
        call grid_metrics_gradients
     else
        ! Cartesian metrics along x
        ! -------------------------
        call grid_metrics_cart(idx(1:nx),x,nx,1,ngh)   ! for inviscid fluxes
        call grid_comm_metrics_cart(idx,nx,1,ngh)
        call grid_metrics_cart(idx_v(1:nx),x,nx,1,ngh_v) ! for viscous fluxes
        call grid_comm_metrics_cart(idx_v,nx,1,ngh_v)

        ! Cartesian metrics along y
        ! -------------------------
        call grid_metrics_cart(idy(1:ny),y,ny,2,ngh)   ! for inviscid fluxes
        call grid_comm_metrics_cart(idy,ny,2,ngh)
        call grid_metrics_cart(idy_v(1:ny),y,ny,2,ngh_v) ! for viscous fluxes
        call grid_comm_metrics_cart(idy_v,ny,2,ngh_v)
     endif

     ! Cartesian metrics along z
     ! -------------------------
     if (is_2d) then
        idz  =1.0_wp/deltaz
        idz_v=1.0_wp/deltaz
     else
        call grid_metrics_cart(idz(1:nz),z,nz,3,ngh)   ! for inviscid fluxes
        call grid_comm_metrics_cart(idz,nz,3,ngh)
        call grid_metrics_cart(idz_v(1:nz),z,nz,3,ngh_v) ! for viscous fluxes
        call grid_comm_metrics_cart(idz_v,nz,3,ngh_v)
     endif

     ! compute near-wall metrics for Cartesian directions
     ! --------------------------------------------------
     if (is_SBP) then
        call grid_metrics_wall_cart_SBP4
     else
        call grid_metrics_wall_cart
     endif

  endif

end subroutine grid_metrics

!===============================================================================
subroutine grid_pp
!===============================================================================
  !> Finalize the grid elements
  !> 1. determine grid center (i0,j0,k0) per block
  !> 2. write grid
  !> 3. communicate global grid among procs
  !> 4. partitioning of the grid on local proc
  !> 5. compute grid metrics
  !> 6. communicate metrics
  !> 7. compute wall normals
!===============================================================================
  use mod_mpi_part
  use mod_grid
  use mod_constant
  use warnstop
  use mod_utils
  use mod_grid_directions ! TO BE CHANGED used to determine deltay -> must be suppressed
  implicit none
  ! ----------------------------------------------------------------------------
  ! integer :: i,j,k
  ! integer :: ibl2read
  ! integer :: ngx_,ngy_,ngz_
  ! ! ----------------------------------------------------------------------------

  ! ! Read global grid per block (fortran IEEE binary file - Att big or little endian)
  ! ! ==========================
  ! if (iproc==iproc_leader(nob(iproc))) print *,'read '//trim(dirDATA)//'grid_bl'//trim(numchar(nob(iproc)))//'.bin'
  ! ! open(194,file=trim(dirDATA)//'grid_bl'//trim(numchar(nob(iproc)))//'.bin', &
  ! !      form='unformatted',status='unknown')
  ! ibl2read=nblr(nob(iproc))
  ! open(194,file=trim(dirDATA)//'grid_bl'//trim(numchar(ibl2read))//'.bin', &
  !      form='unformatted',status='unknown')
  ! rewind(194)
  ! read(194) ngx_
  ! read(194) ngy_
  ! read(194) ngz_
  ! if (is_curv3) then
  !    read(194) (((xgc3(i,j,k),i=1,ngx),j=1,ngy),k=1,ngz)
  !    read(194) (((ygc3(i,j,k),i=1,ngx),j=1,ngy),k=1,ngz)
  !    read(194) (((zgc3(i,j,k),i=1,ngx),j=1,ngy),k=1,ngz)
  ! else
  !    if (is_curv) then
  !       read(194) ((xgc(i,j),i=1,ngx),j=1,ngy)
  !       read(194) ((ygc(i,j),i=1,ngx),j=1,ngy)
  !    else
  !       read(194) (xg(i),i=1,ngx)
  !       read(194) (yg(j),j=1,ngy)
  !    endif
  !    read(194) (zg(k),k=1,ngz)
  ! endif
  ! close(194)

  ! if (iproc==0) print *,'Global grid OK'

  ! Determine global xmin,xmax,ymin,ymax,zmin,zmax
  ! ==============================================
  ! AB: useful for mod_init_flow, init_vel_cyl ?

  !is_read_ex=.false.
  
!!$  if (iproc.eq.iproc_leader(nob(iproc))) then
!!$
!!$     print *,is_read_ex
!!$     
!!$     if (.not.is_read_ex) then
!!$        if (is_curv3) then
!!$           xmin=minval(minval(minval(xgc3(1:ngx,1:ngy,1:ngz),1),1),1)
!!$           xmax=maxval(maxval(maxval(xgc3(1:ngx,1:ngy,1:ngz),1),1),1)
!!$           ymin=minval(minval(minval(ygc3(1:ngx,1:ngy,1:ngz),1),1),1)
!!$           ymax=maxval(maxval(maxval(ygc3(1:ngx,1:ngy,1:ngz),1),1),1)
!!$           zmin=minval(minval(minval(zgc3(1:ngx,1:ngy,1:ngz),1),1),1)
!!$           zmax=maxval(maxval(maxval(zgc3(1:ngx,1:ngy,1:ngz),1),1),1)
!!$        else
!!$           if (is_curv) then
!!$              xmin=minval(minval(xgc(1:ngx,1:ngy),1),1)
!!$              xmax=maxval(maxval(xgc(1:ngx,1:ngy),1),1)
!!$              ymin=minval(minval(ygc(1:ngx,1:ngy),1),1)
!!$              ymax=maxval(maxval(ygc(1:ngx,1:ngy),1),1)
!!$           else
!!$              xmin=minval(xg(1:ngx),1)
!!$              xmax=maxval(xg(1:ngx),1)
!!$              ymin=minval(yg(1:ngy),1)
!!$              ymax=maxval(yg(1:ngy),1)
!!$           endif
!!$           zmin=minval(zg(1:ngz),1)
!!$           zmax=maxval(zg(1:ngz),1)
!!$        endif
!!$
!!$        call MPI_ALLREDUCE(MPI_IN_PLACE,xmin,1,MPI_DOUBLE_PRECISION,MPI_MIN,COMM_interblock,info)
!!$        call MPI_ALLREDUCE(MPI_IN_PLACE,xmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,COMM_interblock,info)
!!$        call MPI_ALLREDUCE(MPI_IN_PLACE,ymin,1,MPI_DOUBLE_PRECISION,MPI_MIN,COMM_interblock,info)
!!$        call MPI_ALLREDUCE(MPI_IN_PLACE,ymax,1,MPI_DOUBLE_PRECISION,MPI_MAX,COMM_interblock,info)
!!$        call MPI_ALLREDUCE(MPI_IN_PLACE,zmin,1,MPI_DOUBLE_PRECISION,MPI_MIN,COMM_interblock,info)
!!$        call MPI_ALLREDUCE(MPI_IN_PLACE,zmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,COMM_interblock,info)
!!$     else
!!$        ! Determine local xmin,xmax,ymin,ymax,zmin,zmax
!!$        ! ==============================================
!!$        if (is_curv3) then
!!$           xmin=minval(minval(minval(xc3(1:nx,1:ny,1:nz),1),1),1)
!!$           xmax=maxval(maxval(maxval(xc3(1:nx,1:ny,1:nz),1),1),1)
!!$           ymin=minval(minval(minval(yc3(1:nx,1:ny,1:nz),1),1),1)
!!$           ymax=maxval(maxval(maxval(yc3(1:nx,1:ny,1:nz),1),1),1)
!!$           zmin=minval(minval(minval(zc3(1:nx,1:ny,1:nz),1),1),1)
!!$           zmax=maxval(maxval(maxval(zc3(1:nx,1:ny,1:nz),1),1),1)
!!$        else
!!$           if (is_curv) then
!!$              xmin=minval(minval(xc(1:nx,1:ny),1),1)
!!$              xmax=maxval(maxval(xc(1:nx,1:ny),1),1)
!!$              ymin=minval(minval(yc(1:nx,1:ny),1),1)
!!$              ymax=maxval(maxval(yc(1:nx,1:ny),1),1)
!!$           else
!!$              xmin=minval(x(1:nx),1)
!!$              xmax=maxval(x(1:nx),1)
!!$              ymin=minval(y(1:ny),1)
!!$              ymax=maxval(y(1:ny),1)
!!$           endif
!!$           zmin=minval(z(1:nz),1)
!!$           zmax=maxval(z(1:nz),1)
!!$        endif
!!$
!!$        call MPI_ALLREDUCE(MPI_IN_PLACE,xmin,1,MPI_DOUBLE_PRECISION,MPI_MIN,COMM_global,info)
!!$        call MPI_ALLREDUCE(MPI_IN_PLACE,xmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,COMM_global,info)
!!$        call MPI_ALLREDUCE(MPI_IN_PLACE,ymin,1,MPI_DOUBLE_PRECISION,MPI_MIN,COMM_global,info)
!!$        call MPI_ALLREDUCE(MPI_IN_PLACE,ymax,1,MPI_DOUBLE_PRECISION,MPI_MAX,COMM_global,info)
!!$        call MPI_ALLREDUCE(MPI_IN_PLACE,zmin,1,MPI_DOUBLE_PRECISION,MPI_MIN,COMM_global,info)
!!$        call MPI_ALLREDUCE(MPI_IN_PLACE,zmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,COMM_global,info)
!!$     endif
!!$  endif

  ! Interblock and Intrablock communications for grids
  ! ==================================================
  if (iproc==0) print *,'comm grid',is_curv3
  !if (igrd.ne.6) call grid_comm
  ! for SRC and igrd=6 (skewed and sinusoidal grids) : pb with double periodicity

  if (iproc==0) print *,'Comm grid OK'
  !call mpistop('stop here',0)

  ! Partitioning of the grid
  ! ========================
  if (.not.is_read_ex) call grid_local

  if (iproc==0) print *,'Local grid OK'

  ! Compute metrics
  ! ===============
  call grid_metrics

  if (iproc==0) print *,'Metrics grid OK'

  ! Check metrics commutations
  ! ==========================
  if (is_curv) call grid_metrics_commutation_11pts

  if ((iproc==0).and.(is_curv)) print *,'check commutations OK'

  ! Enforce metrics
  if ((TGV).or.(CHIT)) then
     idx= 1.0_wp/deltax;   idx_v= 1.0_wp/deltax
     idy= 1.0_wp/deltay;   idy_v= 1.0_wp/deltay
     idz= 1.0_wp/deltaz;   idz_v= 1.0_wp/deltaz
  endif
  if ((CHAN).or.(PHILL)) then
     idz  = 1.0_wp/deltaz
     idz_v= 1.0_wp/deltaz
  endif

  ! Computation of wall normals and BC parameters
  ! =============================================
  if (is_curv) call grid_normals
  
  if (.not.is_curv3) then
  if (is_curv) then
     if ((.not.CYL).and.(.not.TURB).and.(.not.TE)) &
          deltay=abs(ygc(1,2)-ygc(1,1))
     if (is_swapij_bl(3).or.is_swapij_bl(4)) &
          deltay=abs(ygc(2,1)-ygc(1,1))
     if (ACT) deltay = 1d-6!4d-7
  else
     !deltay=1.0_wp/dyg(1)
     !deltay=yg(2)-yg(1)
  endif
  endif
  ! added by XG 04/03/22 for post-processing old channel flow databases
  if (CHAN) then
     if (coord(2)==0) deltay=1.0_wp/idy(1)
  endif

end subroutine grid_pp

