!===============================================================================
subroutine grid_metrics3d_commutation
!===============================================================================
  !> Check 3-D curvilinear metrics conservation
!===============================================================================
  use mod_mpi   ! for nob
  use mod_utils ! for numchar
  use mod_grid_metrics_c3 ! for derivative_ksi,eta,phi
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j,k
  ! relationships Ri=0 to be satisfied
  real(wp), dimension(:,:,:), allocatable :: R1,R2,R3
  ! ---------------------------------------------------------------------------

  ! ==============================================
  ! R1=d(ksi_x)/dksi+d(eta_x)/deta+d(phi_x)/dphi=0
  ! R2=d(ksi_y)/dksi+d(eta_y)/deta+d(phi_y)/dphi=0
  ! R3=d(ksi_z)/dksi+d(eta_z)/deta+d(phi_z)/dphi=0
  ! ==============================================

  ! Metrics on stencil -ngh:+ngh for inviscid fluxes
  ! ================================================

  ! Initializations
  ! ---------------
  allocate(R1(nx1:nx2,ny1:ny2,nz1:nz2))
  allocate(R2(nx1:nx2,ny1:ny2,nz1:nz2))
  allocate(R3(nx1:nx2,ny1:ny2,nz1:nz2))
  R1=0.0_wp
  R2=0.0_wp
  R3=0.0_wp

  ! R1=d(ksi_x)/dksi, ...
  ! ---------------------
  call derivative_ksi(ksi_x,R1,ksi_y,R2,ksi_z,R3)

  ! R1=R1+d(eta_x)/deta, ...
  ! ------------------------
  call derivative_eta(eta_x,R1,eta_y,R2,eta_z,R3)

  ! R1=R1+d(phi_x)/dphi, ...
  ! ------------------------
  call derivative_phi(phi_x,R1,phi_y,R2,phi_z,R3)

  ! Write results for check
  ! =======================
  open(194,file='check_metrics_bl'//trim(numchar(nob(iproc)))//'.bin',form='unformatted',status='unknown')
  rewind(194)
  write(194) ngx
  write(194) ngy
  write(194) ngz
  write(194) (((R1(i,j,k),i=1,ngx),j=1,ngy),k=1,ngz)
  write(194) (((R2(i,j,k),i=1,ngx),j=1,ngy),k=1,ngz)
  write(194) (((R3(i,j,k),i=1,ngx),j=1,ngy),k=1,ngz)
  close(194)

  deallocate(R1,R2,R3)

  call mpistop('come back after check ...',0)

  ! Metrics on stencil -ngh_v:+ngh_v for viscous fluxes
  ! ===================================================

  ! Initializations
  ! ---------------
  allocate(R1(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
  allocate(R2(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
  allocate(R3(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
  R1=0.0_wp
  R2=0.0_wp
  R3=0.0_wp

  ! R1=d(ksi_x)/dksi, ...
  ! ---------------------
  call derivative_ksi_v(ksi_x_v,R1,ksi_y_v,R2,ksi_z_v,R3)

  ! R1=R1+d(eta_x)/deta, ...
  ! ------------------------
  call derivative_eta_v(eta_x_v,R1,eta_y_v,R2,eta_z_v,R3)

  ! R1=R1+d(phi_x)/dphi, ...
  ! ------------------------
  call derivative_phi_v(phi_x_v,R1,phi_y_v,R2,phi_z_v,R3)

  ! Write results for check
  ! =======================
  open(194,file='check_metrics_bl'//trim(numchar(nob(iproc)))//'.bin',form='unformatted',status='unknown')
  rewind(194)
  write(194) ngx
  write(194) ngy
  write(194) ngz
  write(194) (((R1(i,j,k),i=1,ngx),j=1,ngy),k=1,ngz)
  write(194) (((R2(i,j,k),i=1,ngx),j=1,ngy),k=1,ngz)
  write(194) (((R3(i,j,k),i=1,ngx),j=1,ngy),k=1,ngz)
  close(194)

  deallocate(R1,R2,R3)

  call mpistop('come back after check ...',0)

end subroutine grid_metrics3d_commutation

!===============================================================================
subroutine grid_metrics_commutation_11pts
!===============================================================================
  !> Check curvilinear metrics commutations
  !> / 11-point stencil metrics / DEV
!===============================================================================
  use mod_mpi
  use mod_grid
  use mod_coeff_deriv
  use mod_flow
  use mod_utils
  use warnstop
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j!,l
  ! ---------------------------------------------------------------------------

  ! Initializations
  ! ===============
  allocate(x_ksi_eta(nx1:nx2,ny1:ny2,1),y_ksi_eta(nx1:nx2,ny1:ny2,1))
  allocate(x_eta_ksi(nx1:nx2,ny1:ny2,1),y_eta_ksi(nx1:nx2,ny1:ny2,1))
  x_ksi_eta=0.0_wp
  y_ksi_eta=0.0_wp
  x_eta_ksi=0.0_wp
  y_eta_ksi=0.0_wp
  
  ! derivatives along ksi
  ! =====================
  if (BC_face(1,1)%sort==0) then
     i=2 ! wall layer points i=2
     do j=ndy_e,nfy_e
        x_eta_ksi(i,j,1)= 0.5_wp*(x_eta(i+1,j)-x_eta(i-1,j))
        y_eta_ksi(i,j,1)= 0.5_wp*(y_eta(i+1,j)-y_eta(i-1,j))
     enddo

     i=3 ! wall layer points i=3
     do j=ndy_e,nfy_e
        x_eta_ksi(i,j,1)= a5(1)*(x_eta(i+1,j)-x_eta(i-1,j)) &
                        + a5(2)*(x_eta(i+2,j)-x_eta(i-2,j))
        y_eta_ksi(i,j,1)= a5(1)*(y_eta(i+1,j)-y_eta(i-1,j)) &
                        + a5(2)*(y_eta(i+2,j)-y_eta(i-2,j))
     enddo

     i=4 ! wall layer points i=4
     do j=ndy_e,nfy_e
        x_eta_ksi(i,j,1)= a7(1) * ( x_eta(i+1,j)-x_eta(i-1,j) ) &
                        + a7(2) * ( x_eta(i+2,j)-x_eta(i-2,j) ) &
                        + a7(3) * ( x_eta(i+3,j)-x_eta(i-3,j) )
        y_eta_ksi(i,j,1)= a7(1) * ( y_eta(i+1,j)-y_eta(i-1,j) ) &
                        + a7(2) * ( y_eta(i+2,j)-y_eta(i-2,j) ) &
                        + a7(3) * ( y_eta(i+3,j)-y_eta(i-3,j) )
     enddo

     i=5 ! wall layer points i=5
     do j=ndy_e,nfy_e
        x_eta_ksi(i,j,1)= a9(1) * ( x_eta(i+1,j)-x_eta(i-1,j) ) &
                        + a9(2) * ( x_eta(i+2,j)-x_eta(i-2,j) ) &
                        + a9(3) * ( x_eta(i+3,j)-x_eta(i-3,j) ) &
                        + a9(4) * ( x_eta(i+4,j)-x_eta(i-4,j) )      
        y_eta_ksi(i,j,1)= a9(1) * ( y_eta(i+1,j)-y_eta(i-1,j) ) &
                        + a9(2) * ( y_eta(i+2,j)-y_eta(i-2,j) ) &
                        + a9(3) * ( y_eta(i+3,j)-y_eta(i-3,j) ) &
                        + a9(4) * ( y_eta(i+4,j)-y_eta(i-4,j) )      
     enddo
  endif

  do i=ndx,nfx
     do j=ndy_e,nfy_e
        x_eta_ksi(i,j,1)= a11(1)*( x_eta(i+1,j)-x_eta(i-1,j) ) &
                        + a11(2)*( x_eta(i+2,j)-x_eta(i-2,j) ) &
                        + a11(3)*( x_eta(i+3,j)-x_eta(i-3,j) ) &
                        + a11(4)*( x_eta(i+4,j)-x_eta(i-4,j) ) &
                        + a11(5)*( x_eta(i+5,j)-x_eta(i-5,j) )    
        y_eta_ksi(i,j,1)= a11(1)*( y_eta(i+1,j)-y_eta(i-1,j) ) &
                        + a11(2)*( y_eta(i+2,j)-y_eta(i-2,j) ) &
                        + a11(3)*( y_eta(i+3,j)-y_eta(i-3,j) ) &
                        + a11(4)*( y_eta(i+4,j)-y_eta(i-4,j) ) &
                        + a11(5)*( y_eta(i+5,j)-y_eta(i-5,j) )
     enddo
  enddo
  
  if (BC_face(1,2)%sort==0) then
     i=nx-4 ! wall layer points i=nx-4
     do j=ndy_e,nfy_e
        x_eta_ksi(i,j,1)= a9(1) * ( x_eta(i+1,j)-x_eta(i-1,j) ) &
                        + a9(2) * ( x_eta(i+2,j)-x_eta(i-2,j) ) &
                        + a9(3) * ( x_eta(i+3,j)-x_eta(i-3,j) ) &
                        + a9(4) * ( x_eta(i+4,j)-x_eta(i-4,j) )      
        y_eta_ksi(i,j,1)= a9(1) * ( y_eta(i+1,j)-y_eta(i-1,j) ) &
                        + a9(2) * ( y_eta(i+2,j)-y_eta(i-2,j) ) &
                        + a9(3) * ( y_eta(i+3,j)-y_eta(i-3,j) ) &
                        + a9(4) * ( y_eta(i+4,j)-y_eta(i-4,j) )      
     enddo

     i=nx-3 ! wall layer points i=nx-3
     do j=ndy_e,nfy_e
        x_eta_ksi(i,j,1)= a7(1) * ( x_eta(i+1,j)-x_eta(i-1,j) ) &
                        + a7(2) * ( x_eta(i+2,j)-x_eta(i-2,j) ) &
                        + a7(3) * ( x_eta(i+3,j)-x_eta(i-3,j) )
        y_eta_ksi(i,j,1)= a7(1) * ( y_eta(i+1,j)-y_eta(i-1,j) ) &
                        + a7(2) * ( y_eta(i+2,j)-y_eta(i-2,j) ) &
                        + a7(3) * ( y_eta(i+3,j)-y_eta(i-3,j) )
     enddo

     i=nx-2 ! wall layer points i=nx-2
     do j=ndy_e,nfy_e
        x_eta_ksi(i,j,1)= a5(1)*(x_eta(i+1,j)-x_eta(i-1,j)) &
                        + a5(2)*(x_eta(i+2,j)-x_eta(i-2,j))
        y_eta_ksi(i,j,1)= a5(1)*(y_eta(i+1,j)-y_eta(i-1,j)) &
                        + a5(2)*(y_eta(i+2,j)-y_eta(i-2,j))
     enddo

     i=nx-1 ! wall layer points i=nx-1
     do j=ndy_e,nfy_e
        x_eta_ksi(i,j,1)= 0.5_wp*(x_eta(i+1,j)-x_eta(i-1,j))
        y_eta_ksi(i,j,1)= 0.5_wp*(y_eta(i+1,j)-y_eta(i-1,j))
     enddo
  endif

  ! derivatives along eta
  ! =====================
  if (BC_face(2,1)%sort==0) then
     j=2 ! wall layer points j=2
     do i=ndx_e,nfx_e
        x_ksi_eta(i,j,1)= 0.5_wp*(x_ksi(i,j+1)-x_ksi(i,j-1))
        y_ksi_eta(i,j,1)= 0.5_wp*(y_ksi(i,j+1)-y_ksi(i,j-1))
     enddo

     j=3 ! wall layer points j=3
     do i=ndx_e,nfx_e
        x_ksi_eta(i,j,1)= a5(1)*(x_ksi(i,j+1)-x_ksi(i,j-1)) &
                        + a5(2)*(x_ksi(i,j+2)-x_ksi(i,j-2))
        y_ksi_eta(i,j,1)= a5(1)*(y_ksi(i,j+1)-y_ksi(i,j-1)) &
                        + a5(2)*(y_ksi(i,j+2)-y_ksi(i,j-2))
     enddo

     j=4 ! wall layer points j=4
     do i=ndx_e,nfx_e
        x_ksi_eta(i,j,1)= a7(1) * ( x_ksi(i,j+1)-x_ksi(i,j-1) ) &
                        + a7(2) * ( x_ksi(i,j+2)-x_ksi(i,j-2) ) &
                        + a7(3) * ( x_ksi(i,j+3)-x_ksi(i,j-3) )
        y_ksi_eta(i,j,1)= a7(1) * ( y_ksi(i,j+1)-y_ksi(i,j-1) ) &
                        + a7(2) * ( y_ksi(i,j+2)-y_ksi(i,j-2) ) &
                        + a7(3) * ( y_ksi(i,j+3)-y_ksi(i,j-3) )
     enddo

     j=5 ! wall layer points j=5
     do i=ndx_e,nfx_e
        x_ksi_eta(i,j,1)= a9(1) * ( x_ksi(i,j+1)-x_ksi(i,j-1) ) &
                        + a9(2) * ( x_ksi(i,j+2)-x_ksi(i,j-2) ) &
                        + a9(3) * ( x_ksi(i,j+3)-x_ksi(i,j-3) ) &
                        + a9(4) * ( x_ksi(i,j+4)-x_ksi(i,j-4) )      
        y_ksi_eta(i,j,1)= a9(1) * ( y_ksi(i,j+1)-y_ksi(i,j-1) ) &
                        + a9(2) * ( y_ksi(i,j+2)-y_ksi(i,j-2) ) &
                        + a9(3) * ( y_ksi(i,j+3)-y_ksi(i,j-3) ) &
                        + a9(4) * ( y_ksi(i,j+4)-y_ksi(i,j-4) )      
     enddo
  endif

  do i=ndx_e,nfx_e
     do j=ndy,nfy
        x_ksi_eta(i,j,1)= a11(1) * ( x_ksi(i,j+1)-x_ksi(i,j-1) ) &
                        + a11(2) * ( x_ksi(i,j+2)-x_ksi(i,j-2) ) &
                        + a11(3) * ( x_ksi(i,j+3)-x_ksi(i,j-3) ) &
                        + a11(4) * ( x_ksi(i,j+4)-x_ksi(i,j-4) ) &
                        + a11(5) * ( x_ksi(i,j+5)-x_ksi(i,j-5) )
        y_ksi_eta(i,j,1)= a11(1) * ( y_ksi(i,j+1)-y_ksi(i,j-1) ) &
                        + a11(2) * ( y_ksi(i,j+2)-y_ksi(i,j-2) ) &
                        + a11(3) * ( y_ksi(i,j+3)-y_ksi(i,j-3) ) &
                        + a11(4) * ( y_ksi(i,j+4)-y_ksi(i,j-4) ) &
                        + a11(5) * ( y_ksi(i,j+5)-y_ksi(i,j-5) )
      enddo
  enddo
     
  if (BC_face(2,2)%sort==0) then
     j=ny-4 ! wall layer points j=ny-4
     do i=ndx_e,nfx_e
        x_ksi_eta(i,j,1)= a9(1) * ( x_ksi(i,j+1)-x_ksi(i,j-1) ) &
                        + a9(2) * ( x_ksi(i,j+2)-x_ksi(i,j-2) ) &
                        + a9(3) * ( x_ksi(i,j+3)-x_ksi(i,j-3) ) &
                        + a9(4) * ( x_ksi(i,j+4)-x_ksi(i,j-4) )      
        y_ksi_eta(i,j,1)= a9(1) * ( y_ksi(i,j+1)-y_ksi(i,j-1) ) &
                        + a9(2) * ( y_ksi(i,j+2)-y_ksi(i,j-2) ) &
                        + a9(3) * ( y_ksi(i,j+3)-y_ksi(i,j-3) ) &
                        + a9(4) * ( y_ksi(i,j+4)-y_ksi(i,j-4) )      
     enddo

     j=ny-3 ! wall layer points j=ny-3
     do i=ndx_e,nfx_e
        x_ksi_eta(i,j,1)= a7(1) * ( x_ksi(i,j+1)-x_ksi(i,j-1) ) &
                        + a7(2) * ( x_ksi(i,j+2)-x_ksi(i,j-2) ) &
                        + a7(3) * ( x_ksi(i,j+3)-x_ksi(i,j-3) )
        y_ksi_eta(i,j,1)= a7(1) * ( y_ksi(i,j+1)-y_ksi(i,j-1) ) &
                        + a7(2) * ( y_ksi(i,j+2)-y_ksi(i,j-2) ) &
                        + a7(3) * ( y_ksi(i,j+3)-y_ksi(i,j-3) )
     enddo

     j=ny-2 ! wall layer points j=ny-2
     do i=ndx_e,nfx_e
        x_ksi_eta(i,j,1)= a5(1)*(x_ksi(i,j+1)-x_ksi(i,j-1)) &
                        + a5(2)*(x_ksi(i,j+2)-x_ksi(i,j-2))
        y_ksi_eta(i,j,1)= a5(1)*(y_ksi(i,j+1)-y_ksi(i,j-1)) &
                        + a5(2)*(y_ksi(i,j+2)-y_ksi(i,j-2))
     enddo

     j=ny-1 ! wall layer points j=ny-1
     do i=ndx_e,nfx_e
        x_ksi_eta(i,j,1)= 0.5_wp*(x_ksi(i,j+1)-x_ksi(i,j-1))
        y_ksi_eta(i,j,1)= 0.5_wp*(y_ksi(i,j+1)-y_ksi(i,j-1))
     enddo
  endif

!!$  if (iproc==2) then
!!$     j=2
!!$     do i=nx-8,nx
!!$     print *,i,j,x_ksi_eta(i,j,1)-x_eta_ksi(i,j,1)
!!$     enddo
!!$  endif

  !print *,'commutation x',maxval(x_ksi_eta-x_eta_ksi)
  !print *,'commutation y',maxval(y_ksi_eta-y_eta_ksi)


!!$  call mpistop('hhhh',0)
  
!!$  ! write
!!$  ! =====================
!!$  uvar(:,:,:,1)=x_ksi_eta(1:nx,1:ny,:)-x_eta_ksi(1:nx,1:ny,:)
!!$  uvar(:,:,:,2)=y_ksi_eta(1:nx,1:ny,:)-y_eta_ksi(1:nx,1:ny,:)
!!$  call write_plane(1)
!!$  !call read_write_vol('commut_bl'//trim(numchar(nob(iproc)))//filext_write,WRITE)
  
  ! free allocated arrays
  ! =====================
  deallocate(x_ksi_eta,y_ksi_eta)
  deallocate(x_eta_ksi,y_eta_ksi)

end subroutine grid_metrics_commutation_11pts

!!$!===============================================================================
!!$subroutine grid_metrics_commutation_5pts
!!$!===============================================================================
!!$  !> Check curvilinear metrics commutations
!!$  !> / 11-point stencil metrics / DEV
!!$!===============================================================================
!!$  use mod_mpi
!!$  use mod_grid
!!$  use mod_coeff_deriv
!!$  use mod_io_planes
!!$  use mod_utils
!!$  use warnstop
!!$  implicit none
!!$  ! ---------------------------------------------------------------------------
!!$  integer :: i,j,l
!!$  ! ---------------------------------------------------------------------------
!!$
!!$  ! Initializations
!!$  ! ===============
!!$  allocate(x_ksi_eta(nx1:nx2,ny1:ny2,1),y_ksi_eta(nx1:nx2,ny1:ny2,1))
!!$  allocate(x_eta_ksi(nx1:nx2,ny1:ny2,1),y_eta_ksi(nx1:nx2,ny1:ny2,1))
!!$  x_ksi_eta=0.0_wp
!!$  y_ksi_eta=0.0_wp
!!$  x_eta_ksi=0.0_wp
!!$  y_eta_ksi=0.0_wp
!!$
!!$  ! derivatives along ksi
!!$  ! =====================
!!$  if (BC_face(1,1)%sort==0) then
!!$     i=2 ! wall layer points i=2
!!$     do j=ndy_e,nfy_e
!!$        x_eta_ksi(i,j,1)= 0.5_wp*(x_eta(i+1,j)-x_eta(i-1,j))
!!$        y_eta_ksi(i,j,1)= 0.5_wp*(y_eta(i+1,j)-y_eta(i-1,j))
!!$     enddo
!!$  endif
!!$
!!$  do i=ndx,nfx
!!$     do j=ndy_e,nfy_e
!!$        x_eta_ksi(i,j,1)= a5(1)*(x_eta(i+1,j)-x_eta(i-1,j)) &
!!$                        + a5(2)*(x_eta(i+2,j)-x_eta(i-2,j))
!!$        y_eta_ksi(i,j,1)= a5(1)*(y_eta(i+1,j)-y_eta(i-1,j)) &
!!$                        + a5(2)*(y_eta(i+2,j)-y_eta(i-2,j))
!!$     enddo
!!$  enddo
!!$
!!$  if (BC_face(1,2)%sort==0) then
!!$     i=nx-1 ! wall layer points i=nx-1
!!$     do j=ndy_e,nfy_e
!!$        x_eta_ksi(i,j,1)= 0.5_wp*(x_eta(i+1,j)-x_eta(i-1,j))
!!$        y_eta_ksi(i,j,1)= 0.5_wp*(y_eta(i+1,j)-y_eta(i-1,j))
!!$     enddo
!!$  endif
!!$
!!$  ! derivatives along eta
!!$  ! =====================
!!$  if (BC_face(2,1)%sort==0) then
!!$     j=2 ! wall layer points j=2
!!$     do i=ndx_e,nfx_e
!!$        x_ksi_eta(i,j,1)= 0.5_wp*(x_ksi(i,j+1)-x_ksi(i,j-1))
!!$        y_ksi_eta(i,j,1)= 0.5_wp*(y_ksi(i,j+1)-y_ksi(i,j-1))
!!$     enddo
!!$  endif
!!$
!!$  do i=ndx_e,nfx_e
!!$     do j=ndy,nfy
!!$        x_ksi_eta(i,j,1)= a5(1)*(x_ksi(i,j+1)-x_ksi(i,j-1)) &
!!$                        + a5(2)*(x_ksi(i,j+2)-x_ksi(i,j-2))
!!$        y_ksi_eta(i,j,1)= a5(1)*(y_ksi(i,j+1)-y_ksi(i,j-1)) &
!!$                        + a5(2)*(y_ksi(i,j+2)-y_ksi(i,j-2))
!!$      enddo
!!$  enddo
!!$
!!$  if (BC_face(2,2)%sort==0) then
!!$     j=ny-1 ! wall layer points j=ny-1
!!$     do i=ndx_e,nfx_e
!!$        x_ksi_eta(i,j,1)= 0.5_wp*(x_ksi(i,j+1)-x_ksi(i,j-1))
!!$        y_ksi_eta(i,j,1)= 0.5_wp*(y_ksi(i,j+1)-y_ksi(i,j-1))
!!$     enddo
!!$  endif
!!$
!!$  ! free allocated arrays
!!$  ! =====================
!!$  deallocate(x_ksi_eta,y_ksi_eta)
!!$  deallocate(x_eta_ksi,y_eta_ksi)
!!$
!!$end subroutine grid_metrics_commutation_5pts
