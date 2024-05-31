!===============================================================================
subroutine init_interp
!===============================================================================
  !> author: XG
  !> date: April 2021
  !> Initialize flow field by interpolating from old grid
!===============================================================================
  use mod_mpi
  use mod_flow
  use mod_flow_o
  use mod_constant ! <~ for L_ref
  use mod_interp2d_c
  use warnstop
  implicit none
  ! ----------------------------------------------------------------------------
  integer :: nvar,n,k
  integer :: dw_z
  real(wp), dimension(:,:,:), allocatable :: var_iz
  ! ----------------------------------------------------------------------------
  ! options
  logical :: is_interp_z,is_duplicate_z,is_duplicate_only,is_midpoint
  ! ----------------------------------------------------------------------------

  is_duplicate_only=.false.
  
  ! Read old field [DONOR]
  ! ==============
  call read_donor_field('old')

  ! Non-dimensionalize grids before interpolations
  ! ==============================================
  xc_o=xc_o/L_ref
  yc_o=yc_o/L_ref
  z_o=z_o/L_ref

  xc=xc/L_ref
  yc=yc/L_ref
  z=z/L_ref

  ! Choose interpolation stencil size
  ! =================================
  dw=4
  dw_z=2 ! (for z direction if Lagrangian interp)

  ! Perform interpolations
  ! ======================  
  if (is_2D) then
     
     nvar=4 ! 4 conservative variables in 3D
     
     ! Perform 2D curvilinear interpolations
     ! -------------------------------------
     do n=1,nvar
        if (n==1) call interp2d_c(nx_o,ny_o,xc_o,yc_o,rho_o(:,:,1), &
                                  xc(1:nx,1:ny),yc(1:nx,1:ny),rho(1:nx,1:ny,1))
        if (n==2) call interp2d_c(nx_o,ny_o,xc_o,yc_o,rhou_o(:,:,1), &
                                  xc(1:nx,1:ny),yc(1:nx,1:ny),rhou(1:nx,1:ny,1))
        if (n==3) call interp2d_c(nx_o,ny_o,xc_o,yc_o,rhov_o(:,:,1), &
                                  xc(1:nx,1:ny),yc(1:nx,1:ny),rhov(1:nx,1:ny,1))
        if (n==4) call interp2d_c(nx_o,ny_o,xc_o,yc_o,rhoe_o(:,:,1), &
                                  xc(1:nx,1:ny),yc(1:nx,1:ny),rhoe(1:nx,1:ny,1))
     enddo
     rhow=0.0_wp
     
  else
     
     nvar=5 ! 5 conservative variables in 3D

     ! Treat z-direction first
     ! =======================
     
     ! options
     ! -------
     ! 1/ interpolate first in z-direction (mid-point or Lagrangian interp)
     is_interp_z=.true.
     
     if (is_interp_z.and.(abs(deltaz-deltaz_o)<1.0e-16_wp)) then
        if (iproc==0) print *,'deltaz same as deltaz_o, are you sure you want to interpolate in z?'
        call mpistop('Please check',0)
     endif

     if (ndomz>1) then
        call mpistop('ndomz should be 1 for interp_z mode',0)
     endif

     if (is_interp_z.and.(ngz==2*ngz_o).and.(abs(deltaz-0.5_wp*deltaz_o)<1.0e-16_wp)) then
        is_midpoint=.true.
     else
        is_midpoint=.false.
     endif
     
     ! 2/ duplicate z (if same deltaz)     
     is_duplicate_z=.false.

     if (is_duplicate_z.and.(abs(deltaz-deltaz_o)>1.0e-16_wp)) then
        if (iproc==0) print *,'deltaz different from deltaz_o, are you sure you want to duplicate in z?'
        call mpistop('Please check',0)
     endif
     if (is_duplicate_z.and.(nz==nz_o)) then
        if (iproc==0) print *,'nz=nz_o ! Are you sure you want to duplicate in z?'
        call mpistop('Please check',0)
     endif
             
     if (is_duplicate_only) then
           
        ! Duplicate first in z direction
        rho(1:nx,1:ny,1:nz_o)=rho_o(1:nx,1:ny,1:nz_o)
        rho(1:nx,1:ny,nz_o+1:2*nz_o)=rho_o(1:nx,1:ny,1:nz_o)
        rhou(1:nx,1:ny,1:nz_o)=rhou_o(1:nx,1:ny,1:nz_o)
        rhou(1:nx,1:ny,nz_o+1:2*nz_o)=rhou_o(1:nx,1:ny,1:nz_o)
        rhov(1:nx,1:ny,1:nz_o)=rhov_o(1:nx,1:ny,1:nz_o)
        rhov(1:nx,1:ny,nz_o+1:2*nz_o)=rhov_o(1:nx,1:ny,1:nz_o)
        rhow(1:nx,1:ny,1:nz_o)=rhow_o(1:nx,1:ny,1:nz_o)
        rhow(1:nx,1:ny,nz_o+1:2*nz_o)=rhow_o(1:nx,1:ny,1:nz_o)
        rhoe(1:nx,1:ny,1:nz_o)=rhoe_o(1:nx,1:ny,1:nz_o)
        rhoe(1:nx,1:ny,nz_o+1:2*nz_o)=rhoe_o(1:nx,1:ny,1:nz_o)

     else

        ! allocate intermediate field (for interpolation/duplication along z)
        allocate(var_iz(1-ngh:nx_o+ngh,1-ngh:ny_o+ngh,nz))

        do n=1,nvar

           if (is_interp_z) then

              ! Interpolate first along z direction
              ! -----------------------------------
              if (is_midpoint) then
                 ! Simple mid-point interpolation
                 if (n==1) call interp_z( rho_o,var_iz)
                 if (n==2) call interp_z(rhou_o,var_iz)
                 if (n==3) call interp_z(rhov_o,var_iz)
                 if (n==4) call interp_z(rhow_o,var_iz)
                 if (n==5) call interp_z(rhoe_o,var_iz)
              else
                 if (n==1) call interp_z_lagrange( rho_o,var_iz,dw_z)
                 if (n==2) call interp_z_lagrange(rhou_o,var_iz,dw_z)
                 if (n==3) call interp_z_lagrange(rhov_o,var_iz,dw_z)
                 if (n==4) call interp_z_lagrange(rhow_o,var_iz,dw_z)
                 if (n==5) call interp_z_lagrange(rhoe_o,var_iz,dw_z)
                 !call mpistop('Lagrangian interpolation in z not yet implemented',1)
              endif

           elseif (is_duplicate_z) then

              ! Duplicate first in z direction
              ! ------------------------------
              if ((ngz==2*ngz_o).and.(ndomz==1)) then           
                 if (n==1) then
                    var_iz(:,:,1:ngz_o)=rho_o(:,:,1:ngz_o)
                    var_iz(:,:,ngz_o+1:2*ngz_o)=rho_o(:,:,1:ngz_o)
                 endif
                 if (n==2) then
                    var_iz(:,:,1:ngz_o)=rhou_o(:,:,1:ngz_o)
                    var_iz(:,:,ngz_o+1:2*ngz_o)=rhou_o(:,:,1:ngz_o)
                 endif
                 if (n==3) then
                    var_iz(:,:,1:ngz_o)=rhov_o(:,:,1:ngz_o)
                    var_iz(:,:,ngz_o+1:2*ngz_o)=rhov_o(:,:,1:ngz_o)
                 endif
                 if (n==4) then
                    var_iz(:,:,1:ngz_o)=rhow_o(:,:,1:ngz_o)
                    var_iz(:,:,ngz_o+1:2*ngz_o)=rhow_o(:,:,1:ngz_o)
                 endif
                 if (n==5) then
                    var_iz(:,:,1:ngz_o)=rhoe_o(:,:,1:ngz_o)
                    var_iz(:,:,ngz_o+1:2*ngz_o)=rhoe_o(:,:,1:ngz_o)
                 endif
              else
                 if (ngz.ne.2*ngz_o) then
                    call mpistop('cannot duplicate, check dim',1)
                 else
                    call mpistop('ndomz should be 1 for duplicate mode',1)
                 endif
              endif

           else

              ! Simple copy in intermediate array
              ! ---------------------------------
              if (n==1) then
                 var_iz(:,:,1:nz_o)=rho_o(:,:,1:nz_o)
              endif
              if (n==2) then
                 var_iz(:,:,1:nz_o)=rhou_o(:,:,1:nz_o)
              endif
              if (n==3) then
                 var_iz(:,:,1:nz_o)=rhov_o(:,:,1:nz_o)
              endif
              if (n==4) then
                 var_iz(:,:,1:nz_o)=rhow_o(:,:,1:nz_o)
              endif
              if (n==5) then
                 var_iz(:,:,1:nz_o)=rhoe_o(:,:,1:nz_o)
              endif

           endif

           ! Perform 2D curvilinear interpolations
           ! -------------------------------------
           do k=1,nz
              !print *,'iproc',n,k,iproc
              if (n==1) call interp2d_c(nx_o,ny_o,xc_o,yc_o,var_iz(:,:,k), &
                   xc(1:nx,1:ny),yc(1:nx,1:ny),rho(1:nx,1:ny,k))
              if (n==2) call interp2d_c(nx_o,ny_o,xc_o,yc_o,var_iz(:,:,k), &
                   xc(1:nx,1:ny),yc(1:nx,1:ny),rhou(1:nx,1:ny,k))
              if (n==3) call interp2d_c(nx_o,ny_o,xc_o,yc_o,var_iz(:,:,k), &
                   xc(1:nx,1:ny),yc(1:nx,1:ny),rhov(1:nx,1:ny,k))
              if (n==4) call interp2d_c(nx_o,ny_o,xc_o,yc_o,var_iz(:,:,k), &
                   xc(1:nx,1:ny),yc(1:nx,1:ny),rhow(1:nx,1:ny,k))
              if (n==5) call interp2d_c(nx_o,ny_o,xc_o,yc_o,var_iz(:,:,k), &
                   xc(1:nx,1:ny),yc(1:nx,1:ny),rhoe(1:nx,1:ny,k))
           enddo

           if (iproc_leader(nob(iproc))) print *,'interp',n,iproc
        enddo
     
     endif
     
  endif
  
  ! Re-dimensionalize current grid
  ! ===============================
  xc=xc*L_ref
  yc=yc*L_ref
  z=z*L_ref

  ! Deallocate old field
  ! ====================
  !call free_field_o
  
end subroutine init_interp

!=============================================================================
subroutine interp_z(u1,u1i)
!=============================================================================
  !> Mid-point interpolations along z
  !> if nz2=2*nz1 (doubling regular grid along z)
!=============================================================================
  use mod_grid
  use mod_flow_o
  implicit none
  ! --------------------------------------------------------------------------
  real(wp), dimension(1-ngh:nx_o+ngh,1-ngh:ny_o+ngh,1-ngh:nz_o+ngh) :: u1
  real(wp), dimension(1-ngh:nx_o+ngh,1-ngh:ny_o+ngh,nz) :: u1i
  ! --------------------------------------------------------------------------
  integer :: i,j,k
  ! --------------------------------------------------------------------------

  do i=1-ngh,nx_o+ngh
     do j=1-ngh,ny_o+ngh

        do k=1,nz
           if (mod(k+1,2)==0) then 
              ! coincident planes
              !print *,k,(k+1)/2
              u1i(i,j,k)=u1(i,j,(k+1)/2)
           else
              !print *,k,k/2,k/2+1
              u1i(i,j,k)=0.5_wp*(u1(i,j,k/2)+u1(i,j,k/2+1))
           endif
        enddo

     enddo
  enddo

end subroutine interp_z

!=============================================================================
subroutine interp_z_lagrange(u1,u1i,dw_z)
!=============================================================================
  !> Mid-point interpolations along z
  !> if nz2=2*nz1 (doubling regular grid along z)
!=============================================================================
  use mod_grid
  use mod_flow_o
  use mod_mpi_part
  implicit none
  ! --------------------------------------------------------------------------
  integer, intent(in) :: dw_z
  real(wp), dimension(1-ngh:nx_o+ngh,1-ngh:ny_o+ngh,1-ngh:nz_o+ngh) :: u1
  real(wp), dimension(1-ngh:nx_o+ngh,1-ngh:ny_o+ngh,nz) :: u1i
  ! --------------------------------------------------------------------------
  integer :: i,j,k,kk,kp,ks,l,m
  real(wp) :: min_d,dist
  real(wp), dimension(dw_z) :: coeff_lagz
  ! --------------------------------------------------------------------------
 
  do i=1-ngh,nx_o+ngh
     do j=1-ngh,ny_o+ngh

        do k=1,nz

           ! closest point in the donor grid (#1)
           ! -------------------------------
           min_d=1.e9
           do kk=1,nz_o
              dist=abs(z_o(kk)-z(k))
              if (dist<min_d) then
                 kp=kk
                 min_d=dist
              endif
           enddo

           if (z_o(kp)<=z(k)) then
              ks=kp-dw_z/2+2-1
           else
              ks=kp-dw_z/2+1-1
           endif

           ! z-direction
           coeff_lagz=1.0_wp
           do l=1,dw_z
              do m=1,dw_z
                 if (m.ne.l) then
                    coeff_lagz(l)=coeff_lagz(l)*(z(k)-z_o(ks+(m-1)))/(z_o(ks+(l-1))-z_o(ks+(m-1)))
                 endif
              enddo
           enddo

           u1i(i,j,k)=0.0_wp
           do l=1,dw_z
              u1i(i,j,k)=u1i(i,j,k)+coeff_lagz(l)*u1(i,j,ks+(l-1))
           enddo
           
        enddo

     enddo
  enddo

end subroutine interp_z_lagrange
