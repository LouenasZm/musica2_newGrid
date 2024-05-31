!===============================================================================
subroutine timestep_new
!===============================================================================
  !> Computation of the global timestep
  ! new writing May 2022 by XG
!===============================================================================
  use mod_mpi
  use mod_time
  use mod_flow     ! <- for uu,vv,ww,c_
  use mod_constant ! <- for ..
  use warnstop
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j,k
  ! minimum mesh size
  integer :: ind_min(2),iz_min(1),ip_min(1),ind3_min(3)
  integer, dimension(:), allocatable :: ind_i,ind_j,ind_k
  real(wp) :: delta_min,d_min
  real(wp), dimension(nx,ny) :: d2D_min
  real(wp), dimension(nx,ny,nz) :: d3D_min
  real(wp), dimension(nz) :: dz_min
  real(wp), dimension(:), allocatable :: dp_min
  logical :: is_min_z
  ! maximum velocity
  real(wp) :: u_max,u_maxp,tf
  ! ---------------------------------------------------------------------------
  real(wp) :: deltat_visq
  
  if (iproc==0) then
     print *,repeat('=',70)
     print *,'Compute timestep'
  endif

  ! Determination of minimal mesh size
  ! ==================================

  ! In 2-D slices (third direction is always extruded for the moment)
  ! -------------
  ! initialize array with large values
  d3D_min=1.0e6_wp
  d2D_min=1.0e6_wp
  dz_min =1.0e6_wp

  if (is_curv3) then ! for full 3D curvilinear grids
     ! minimal mesh size along i
     do k=1,nz
        do j=1,ny
           do i=1,nx-1
              d3D_min(i,j,k)=min(d3D_min(i,j,k),sqrt((xc3(i+1,j,k)-xc3(i,j,k))**2 &
                                                    +(yc3(i+1,j,k)-yc3(i,j,k))**2 &
                                                    +(zc3(i+1,j,k)-zc3(i,j,k))**2))
           enddo
        enddo
     enddo
     ! minimal mesh size along j
     do k=1,nz
        do j=1,ny-1
           do i=1,nx
              d3D_min(i,j,k)=min(d3D_min(i,j,k),sqrt((xc3(i,j+1,k)-xc3(i,j,k))**2 &
                                                    +(yc3(i,j+1,k)-yc3(i,j,k))**2 &
                                                    +(zc3(i,j+1,k)-zc3(i,j,k))**2))
           enddo
        enddo
     enddo
     ! minimal mesh size along k
     do k=1,nz-1
        do j=1,ny
           do i=1,nx
              d3D_min(i,j,k)=min(d3D_min(i,j,k),sqrt((xc3(i,j,k+1)-xc3(i,j,k))**2 &
                                                    +(yc3(i,j,k+1)-yc3(i,j,k))**2 &
                                                    +(zc3(i,j,k+1)-zc3(i,j,k))**2))
           enddo
        enddo
     enddo
     ! minimal mesh size along first ij diagonal
     do k=1,nz
        do j=1,ny-1
           do i=1,nx-1
              d3D_min(i,j,k)=min(d3D_min(i,j,k),sqrt((xc3(i+1,j+1,k)-xc3(i,j,k))**2 &
                                                    +(yc3(i+1,j+1,k)-yc3(i,j,k))**2 &
                                                    +(zc3(i+1,j+1,k)-zc3(i,j,k))**2))
           enddo
        enddo
     enddo
     ! minimal mesh size along second ij diagonal
     do k=1,nz
        do j=1,ny-1
           do i=1,nx-1
              d3D_min(i,j,k)=min(d3D_min(i,j,k),sqrt((xc3(i,j+1,k)-xc3(i+1,j,k))**2 &
                                                    +(yc3(i,j+1,k)-yc3(i+1,j,k))**2 &
                                                    +(zc3(i,j+1,k)-zc3(i+1,j,k))**2))
           enddo
        enddo
     enddo
     ! minimal mesh size along first ik diagonal
     do k=1,nz
        do j=1,ny-1
           do i=1,nx-1
              d3D_min(i,j,k)=min(d3D_min(i,j,k),sqrt((xc3(i+1,j,k+1)-xc3(i,j,k))**2 &
                                                    +(yc3(i+1,j,k+1)-yc3(i,j,k))**2 &
                                                    +(zc3(i+1,j,k+1)-zc3(i,j,k))**2))
           enddo
        enddo
     enddo
     ! minimal mesh size along second ik diagonal
     do k=1,nz
        do j=1,ny-1
           do i=1,nx-1
              d3D_min(i,j,k)=min(d3D_min(i,j,k),sqrt((xc3(i,j,k+1)-xc3(i+1,j,k))**2 &
                                                    +(yc3(i,j,k+1)-yc3(i+1,j,k))**2 &
                                                    +(zc3(i,j,k+1)-zc3(i+1,j,k))**2))
           enddo
        enddo
     enddo
     ! minimal mesh size along first jk diagonal
     do k=1,nz
        do j=1,ny-1
           do i=1,nx-1
              d3D_min(i,j,k)=min(d3D_min(i,j,k),sqrt((xc3(i,j+1,k+1)-xc3(i,j,k))**2 &
                                                    +(yc3(i,j+1,k+1)-yc3(i,j,k))**2 &
                                                    +(zc3(i,j+1,k+1)-zc3(i,j,k))**2))
           enddo
        enddo
     enddo
     ! minimal mesh size along second jk diagonal
     do k=1,nz
        do j=1,ny-1
           do i=1,nx-1
              d3D_min(i,j,k)=min(d3D_min(i,j,k),sqrt((xc3(i,j+1,k)-xc3(i,j,k+1))**2 &
                                                    +(yc3(i,j+1,k)-yc3(i,j,k+1))**2 &
                                                    +(zc3(i,j+1,k)-zc3(i,j,k+1))**2))
           enddo
        enddo
     enddo

     ! Search smallest mesh size (and its location)
     ! --------------------------------------------
     ! location of min for a proc
     ind3_min=minloc(d3D_min)

     ! global minimum
     d_min=d3D_min(ind3_min(1),ind3_min(2),ind3_min(3))

  else
     if (is_curv) then ! for curvilinear grids
        ! minimal mesh size along i
        do j=1,ny
           do i=1,nx-1
              d2D_min(i,j)=min(d2D_min(i,j),sqrt((xc(i+1,j)-xc(i,j))**2+(yc(i+1,j)-yc(i,j))**2))
           enddo
        enddo
        ! minimal mesh size along j
        do j=1,ny-1
           do i=1,nx
              d2D_min(i,j)=min(d2D_min(i,j),sqrt((xc(i,j+1)-xc(i,j))**2+(yc(i,j+1)-yc(i,j))**2))
           enddo
        enddo
        ! minimal mesh size along first diagonal
        do j=1,ny-1
           do i=1,nx-1
              d2D_min(i,j)=min(d2D_min(i,j),sqrt((xc(i+1,j+1)-xc(i,j))**2+(yc(i+1,j+1)-yc(i,j))**2))
           enddo
        enddo
        ! minimal mesh size along second diagonal
        do j=1,ny-1
           do i=1,nx-1
              d2D_min(i,j)=min(d2D_min(i,j),sqrt((xc(i,j+1)-xc(i+1,j))**2+(yc(i,j+1)-yc(i+1,j))**2))
           enddo
        enddo
     else ! for Cartesian grids
        ! minimal mesh size along i
        do j=1,ny
           do i=1,nx-1
              d2D_min(i,j)=min(d2D_min(i,j),x(i+1)-x(i))
           enddo
        enddo
        ! minimal mesh size along j
        do j=1,ny-1
           do i=1,nx
              d2D_min(i,j)=min(d2D_min(i,j),y(j+1)-y(j))
           enddo
        enddo
     endif

     ! Along z (third direction)
     ! -------
     ! minimal mesh size along k
     if (.not.is_2D) then
        do k=1,nz-1
           dz_min(k)=z(k+1)-z(k)
        enddo
     endif

     ! Search smallest mesh size (and its location)
     ! --------------------------------------------
     ! location of min for a proc
     ind_min=minloc(d2D_min)
     iz_min=minloc(dz_min)

     ! minimum taking into account deltaz
     d_min=min(d2D_min(ind_min(1),ind_min(2)),dz_min(iz_min(1)))

     ! check if z has the limiting dimension
     is_min_z=.false.
     if (d_min==dz_min(iz_min(1))) is_min_z=.true.

  endif

  ! search min among MPI domains
  call MPI_ALLREDUCE(d_min,delta_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,COMM_global,info)

!!$  !!!!! TEMP
!!$  if (is_curv3) then
!!$     !d_min=0.002*L_ref
!!$     !d_min=0.005*L_ref
!!$     !d_min=0.012*L_ref
!!$     !d_min=0.1
!!$     !d_min=1.
!!$     ! annulus
!!$     d_min=1.0e-3_wp
!!$     ! Baumgartner
!!$     d_min=1.0e-3_wp/60.0_wp*L_ref
!!$     ! LS59
!!$     !d_min=1.0e-6_wp
!!$  endif

!!$  do k=1,nz
!!$     do j=1,ny
!!$        do i=1,nx
!!$           uvar(i,j,k,2)=d3D_min(i,j,k)
!!$           uvar(i,j,k,3)=(1./ijacob3(i,j,k))**(1./3.)
!!$        enddo
!!$     enddo
!!$  enddo

  ! Indicate to user where is the minimal mesh size
  ! -----------------------------------------------
  allocate(ind_i(0:nproc-1),ind_j(0:nproc-1),ind_k(0:nproc-1))
  allocate(dp_min(1:nproc))

  ! gather d_min on proc 0
  call MPI_GATHER(d_min,1,MPI_DOUBLE_PRECISION,dp_min,1,MPI_DOUBLE_PRECISION,0,COMM_global,info)
  ! gather indices on proc 0
  if (is_curv3) then ! for full 3D curvilinear grids
     call MPI_GATHER(ind3_min(1),1,MPI_INTEGER,ind_i,1,MPI_INTEGER,0,COMM_global,info)
     call MPI_GATHER(ind3_min(2),1,MPI_INTEGER,ind_j,1,MPI_INTEGER,0,COMM_global,info)
     call MPI_GATHER(ind3_min(3),1,MPI_INTEGER,ind_k,1,MPI_INTEGER,0,COMM_global,info)
  else
     call MPI_GATHER(ind_min(1),1,MPI_INTEGER,ind_i,1,MPI_INTEGER,0,COMM_global,info)
     call MPI_GATHER(ind_min(2),1,MPI_INTEGER,ind_j,1,MPI_INTEGER,0,COMM_global,info)
     call MPI_GATHER( iz_min   ,1,MPI_INTEGER,ind_k,1,MPI_INTEGER,0,COMM_global,info)
  endif

  if (iproc==0) then
     
     ! first proc which has min
     ip_min=minloc(dp_min)-1

     ! proc 0 print min at screen
     !!print *,'comp',dp_min(ip_min(1)+1),delta_min
     print 100,dp_min(ip_min(1)+1)
100  format(1x,'~> minimum mesh size= ',e18.11)
     print 101,ind_i(ip_min(1)),ind_j(ip_min(1)),ind_k(ip_min(1))
101  format(1x,'~> located at (i,j,k)=(',i5,',',i5,',',i5,')')
     print 102,nob(ip_min(1)),ip_min(1)
102  format(1x,'~> in block ',i4,' (proc ',i5,')')

  endif
  deallocate(ind_i,ind_j,ind_k,dp_min)

  ! Determination of max velocity norm
  ! ==================================
  u_maxp=0.0_wp
  do k=1,nz
     do j=1,ny
        do i=1,nx
           u_max=max(u_maxp,sqrt(uu(i,j,k)**2+vv(i,j,k)**2+ww(i,j,k)**2)+c_(i,j,k))
           u_maxp=u_max
        enddo
     enddo
  enddo

  ! search max among MPI domains
  call MPI_ALLREDUCE(MPI_IN_PLACE,u_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,COMM_global,info)

  ! Compute deltat
  ! ==============
  deltat= CFL*delta_min/u_max

  if (iorder_visc.ne.0) then
     if (iproc==0) print *,deltat,'CFL',CFL

     deltat_visq= 0.5*delta_min**2*rho_ref/mu_ref

     if (iproc==0) print *,deltat_visq,'Fourier'
  endif

!!$  deltat=1.444378773634171E-003

!!deltat=1.E-6

  ! Select between viscous/inviscid timestep condition
  
  ! Nondimensional timestep
  ! -----------------------
  dtstar = deltat/tscale

  ! Determine nmax from timemax
  ! -----------------------
  tf = timemax*tscale
  !print *,delta_min,deltat,tf,timemax,tscale!,tf/deltat

  nmax=min(nmax,int(tf/deltat)+1)

  if (int(dtstar*nmax).ge.10000) call mpistop("/!\ tstar by the end of the simulation will be greater than 9999, problem with timestamp.",0)

  if (mod(ntime,nprint).eq.0) then  !!!!!  FOR dtvar TO BE CHANGED
     if (iproc==0) then
        print *,'Dt [s], Dt* [-]:',deltat,dtstar
     endif
  endif
  
  if (iproc==0) print *,repeat('=',70)

  ! local timestep
  ! ==============
  if (is_dtlocal) then
     allocate(dt_local(nx,ny,nz))

     ! tests compatibility: .not.is_RFM  .not.is_eigenmode

     !!!! TO BE WRITTEN
  endif
  ! /!\ some options are not compatible with the use of local timestep
  ! where time-dependent source terms or boundary conditions are imposed
  ! -> is_source; is_eigenmode; suction and blowing; is_RFM in BC ......

  
end subroutine timestep_new

!===============================================================================
subroutine local_cfl
!===============================================================================
  !> Computation of local CFL
!===============================================================================
  use mod_time  ! <- for deltat
  use mod_flow  ! <- for uu,vv,ww,c_,metrics
  implicit none
  ! ----------------------------------------------------------------------------
  integer :: i,j,k
  ! ----------------------------------------------------------------------------

  ! Compute local CFL number per direction
  ! ======================================
  if (is_curv3) then

     do k=1,nz
        do j=1,ny
           do i=1,nx
              ! calculation of max characteristic speed
              cfl_i(i,j,k)=abs((uu(i,j,k)*ksi_x(i,j,k)+vv(i,j,k)*ksi_y(i,j,k)+ww(i,j,k)*ksi_z(i,j,k))*ijacob3(i,j,k)) &
                          +c_(i,j,k)*g3_ksi(i,j,k)
              ! compute cfl_ksi
              cfl_i(i,j,k)=cfl_i(i,j,k)*deltat

              ! calculation of max characteristic speed
              cfl_j(i,j,k)=abs((uu(i,j,k)*eta_x(i,j,k)+vv(i,j,k)*eta_y(i,j,k)+ww(i,j,k)*eta_z(i,j,k))*ijacob3(i,j,k)) &
                          +c_(i,j,k)*g3_eta(i,j,k)
              ! compute cfl_eta
              cfl_j(i,j,k)=cfl_j(i,j,k)*deltat

              ! calculation of max characteristic speed
              cfl_k(i,j,k)=abs((uu(i,j,k)*phi_x(i,j,k)+vv(i,j,k)*phi_y(i,j,k)+ww(i,j,k)*phi_z(i,j,k))*ijacob3(i,j,k)) &
                          +c_(i,j,k)*g3_phi(i,j,k)
              ! compute cfl_phi
              cfl_k(i,j,k)=cfl_k(i,j,k)*deltat
           enddo
        enddo
     enddo

  else

     ! Compute local CFL number along i-direction
     ! ==========================================
     if (is_curv) then
        do k=1,nz
           do j=1,ny
              do i=1,nx
                 ! calculation of max characteristic speed
                 cfl_i(i,j,k)=abs((uu(i,j,k)*y_eta(i,j)-vv(i,j,k)*x_eta(i,j))*ijacob(i,j)) &
                             +c_(i,j,k)*g_ksi(i,j)
                 ! compute cfl_ksi
                 cfl_i(i,j,k)=cfl_i(i,j,k)*deltat
              enddo
           enddo
        enddo
     else
        do k=1,nz
           do j=1,ny
              do i=1,nx
                 cfl_i(i,j,k)=(abs(uu(i,j,k))+c_(i,j,k))*deltat*idx(i)
              enddo
           enddo
        enddo
     endif

     ! Compute local CFL number along j-direction
     ! ==========================================
     if (is_curv) then
        do k=1,nz
           do j=1,ny
              do i=1,nx
                 ! calculation of max characteristic speed
                 cfl_j(i,j,k)=abs((vv(i,j,k)*x_ksi(i,j)-uu(i,j,k)*y_ksi(i,j))*ijacob(i,j)) &
                             +c_(i,j,k)*g_eta(i,j)
                 ! compute cfl_eta
                 cfl_j(i,j,k)=cfl_j(i,j,k)*deltat
              enddo
           enddo
        enddo
     else
        do k=1,nz
           do j=1,ny
              do i=1,nx
                 cfl_j(i,j,k)=(abs(vv(i,j,k))+c_(i,j,k))*deltat*idy(i)
              enddo
           enddo
        enddo
     endif

  endif

end subroutine local_cfl

!===============================================================================
subroutine mean_local_cfl
!===============================================================================
  !> Computation of mean local CFL
!===============================================================================
  use mod_time  ! <- for deltat
  use mod_flow  ! <- for uu,vv,ww,c_,metrics
  implicit none
  ! ----------------------------------------------------------------------------
  integer :: i,j,k
  real(wp) :: nn,inn
  ! ----------------------------------------------------------------------------
  
  nn = dble(ntotal)
  inn= 1.0_wp/nn

  ! Compute local CFL number
  ! ========================
  if (is_curv3) then

     do k=1,nz
        do j=1,ny
           do i=1,nx
              ! along i-direction
              cfl_i(i,j,k)=( (nn-1.0_wp)*cfl_i(i,j,k) + &
                   (abs((uu(i,j,k)*ksi_x(i,j,k)+vv(i,j,k)*ksi_y(i,j,k)+ww(i,j,k)*ksi_z(i,j,k))*ijacob3(i,j,k)) &
                   +c_(i,j,k)*g3_ksi(i,j,k))*deltat )*inn
              ! along j-direction
              cfl_j(i,j,k)=( (nn-1.0_wp)*cfl_j(i,j,k) + &
                   (abs((uu(i,j,k)*eta_x(i,j,k)+vv(i,j,k)*eta_y(i,j,k)+ww(i,j,k)*eta_z(i,j,k))*ijacob3(i,j,k)) &
                   +c_(i,j,k)*g3_eta(i,j,k))*deltat )*inn
              ! along k-direction
              cfl_k(i,j,k)=( (nn-1.0_wp)*cfl_k(i,j,k) + &
                   (abs((uu(i,j,k)*phi_x(i,j,k)+vv(i,j,k)*phi_y(i,j,k)+ww(i,j,k)*phi_z(i,j,k))*ijacob3(i,j,k)) &
                   +c_(i,j,k)*g3_phi(i,j,k))*deltat )*inn
           enddo
        enddo
     enddo

  else

     if (is_curv) then
        do k=1,nz
           do j=1,ny
              do i=1,nx
                 ! along i-direction
                 cfl_i(i,j,k)=( (nn-1.0_wp)*cfl_i(i,j,k) + &
                      (abs((uu(i,j,k)*y_eta(i,j)-vv(i,j,k)*x_eta(i,j))*ijacob(i,j)) &
                      +c_(i,j,k)*g_ksi(i,j))*deltat )*inn
                 ! along j-direction
                 cfl_j(i,j,k)=( (nn-1.0_wp)*cfl_j(i,j,k) + &
                      (abs((vv(i,j,k)*x_ksi(i,j)-uu(i,j,k)*y_ksi(i,j))*ijacob(i,j)) &
                      +c_(i,j,k)*g_eta(i,j))*deltat )*inn
              enddo
           enddo
        enddo
     else
        do k=1,nz
           do j=1,ny
              do i=1,nx
                 ! along i-direction
                 cfl_i(i,j,k)=( (nn-1.0_wp)*cfl_i(i,j,k) + &
                      (abs(uu(i,j,k))+c_(i,j,k))*deltat*idx(i) )*inn
                 ! along j-direction
                 cfl_j(i,j,k)=( (nn-1.0_wp)*cfl_j(i,j,k) + &
                      (abs(vv(i,j,k))+c_(i,j,k))*deltat*idy(j) )*inn
              enddo
           enddo
        enddo
     endif

  endif

end subroutine mean_local_cfl

!===============================================================================
subroutine local_timestep
!===============================================================================
  !> Computation of local timestep (with constant CFL and Fourier criteria)
  !> Option #1: min per direction
  !> Option #2: min of d_min and max of local characteristic speed at (i,j,k)
!===============================================================================
  use mod_time  ! <- for dt_local
  use mod_flow  ! <- for uu,vv,ww,c_,metrics
  use mod_mpi
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j,k
  real(wp) :: sr_i,sr_j,sr_k ! spectral radii per direction
  real(wp) :: dt_min ! minimum dt in computation
  real(wp) :: sigma ! control for local timestep (0 or 1)
  real(wp) :: alpha ! ratio of dt_local/dt_min
  ! ---------------------------------------------------------------------------

  ! ~~> Option #1: min per direction

  if (is_curv3) then

     ! Compute local timestep along i,j,k directions
     ! =============================================
     do k=1,nz
        do j=1,ny
           do i=1,nx
              ! calculation of max characteristic speed (along ksi)
              sr_i= abs((uu(i,j,k)*ksi_x(i,j,k)+vv(i,j,k)*ksi_y(i,j,k)+ww(i,j,k)*ksi_z(i,j,k))*ijacob3(i,j,k)) &
                  + c_(i,j,k)*g3_ksi(i,j,k)
              ! calculation of max characteristic speed (along eta)
              sr_j= abs((uu(i,j,k)*eta_x(i,j,k)+vv(i,j,k)*eta_y(i,j,k)+ww(i,j,k)*eta_z(i,j,k))*ijacob3(i,j,k)) &
                  + c_(i,j,k)*g3_eta(i,j,k)
              ! calculation of max characteristic speed (along eta)
              sr_k= abs((uu(i,j,k)*phi_x(i,j,k)+vv(i,j,k)*phi_y(i,j,k)+ww(i,j,k)*phi_z(i,j,k))*ijacob3(i,j,k)) &
                  + c_(i,j,k)*g3_phi(i,j,k)
              ! compute local timestep
              dt_local(i,j,k)=min(CFL/sr_i,CFL/sr_j,CFL/sr_k,100.*deltat)

              !if (dt_local(i,j,k)>100.*deltat) dt_local(i,j,k)=100.*deltat

              !uvar(i,j,k,1) = dt_local(i,j,k)/deltat


!!$              ! calculation of max speed
!!$              sr_i=sqrt((uu(i,j,k)*ksi_x(i,j,k)+vv(i,j,k)*ksi_y(i,j,k)+ww(i,j,k)*ksi_z(i,j,k))**2 &
!!$                       +(uu(i,j,k)*eta_x(i,j,k)+vv(i,j,k)*eta_y(i,j,k)+ww(i,j,k)*eta_z(i,j,k))**2 &
!!$                       +(uu(i,j,k)*phi_x(i,j,k)+vv(i,j,k)*phi_y(i,j,k)+ww(i,j,k)*phi_z(i,j,k))**2)*ijacob3(i,j,k) &
!!$                       +c_(i,j,k)*sqrt(g3_ksi(i,j,k)**2+g3_eta(i,j,k)**2+g3_phi(i,j,k)**2)
!!$              ! compute local timestep
!!$              dt_local(i,j,k)=CFL/sr_i
           enddo
        enddo
     enddo

!!$     ! Control how much dt_local is of dt_min
!!$     ! ======================================
!!$
!!$     ! Get min dt in whole computations
!!$     call MPI_ALLREDUCE(minval(dt_local),dt_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,COMM_global,info)
!!$
!!$     ! sigma = 1 ~> dt_local = dt_local
!!$     ! sigma = 0 ~> dt_local = dt_min
!!$     sigma = 1.0_wp
!!$     do k=1,nz
!!$        do j=1,ny
!!$           do i=1,nx
!!$              alpha = dt_local(i,j,k)/dt_min
!!$              dt_local(i,j,k) = dt_min*(1.0_wp-sigma) + sigma*alpha*dt_min
!!$              !uvar(i,j,k,7) = dt_local(i,j,k)
!!$           enddo
!!$        enddo
!!$     enddo

  else

     ! Compute local timestep along i & j directions
     ! =============================================
     if (is_curv) then
        do k=1,nz
           do j=1,ny
              do i=1,nx
                 ! calculation of max characteristic speed (along ksi)
                 sr_i= abs((uu(i,j,k)*y_eta(i,j)-vv(i,j,k)*x_eta(i,j))*ijacob(i,j)) &
                     + c_(i,j,k)*g_ksi(i,j)
                 ! calculation of max characteristic speed (along eta)
                 sr_j= abs((vv(i,j,k)*x_ksi(i,j)-uu(i,j,k)*y_ksi(i,j))*ijacob(i,j)) &
                     + c_(i,j,k)*g_eta(i,j)
                 ! compute local timestep
                 dt_local(i,j,k)=min(CFL/sr_i,CFL/sr_j)
              enddo
           enddo
        enddo
     else
        do k=1,nz
           do j=1,ny
              do i=1,nx
                 ! calculation of max characteristic speed (along x)
                 sr_i= (abs(uu(i,j,k))+c_(i,j,k))*idx(i)
                 ! calculation of max characteristic speed (along y)
                 sr_j= (abs(vv(i,j,k))+c_(i,j,k))*idy(j)
                 ! compute local timestep
                 dt_local(i,j,k)=min(CFL/sr_i,CFL/sr_j)
              enddo
           enddo
        enddo
     endif

     ! Control how much dt_local is of dt_min
     ! ======================================

     ! Get min dt in whole computations
     call MPI_ALLREDUCE(minval(dt_local),dt_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,COMM_global,info)

     ! sigma = 1 ~> dt_local = dt_local
     ! sigma = 0 ~> dt_local = dt_min
     sigma = 1.0_wp
     do k=1,nz
        do j=1,ny
           do i=1,nx
              alpha = dt_local(i,j,k)/dt_min
              dt_local(i,j,k) = dt_min*(1.0_wp-sigma) + sigma*alpha*dt_min
              !uvar(i,j,k,7) = dt_local(i,j,k)
           enddo
        enddo
     enddo

     !****************
     if (is_2D) return
     !****************

     ! Compute local timestep along k-direction
     ! ========================================
     do k=1,nz
        do j=1,ny
           do i=1,nx
              ! calculation of max characteristic speed (along z)
              sr_k= (abs(ww(i,j,k))+c_(i,j,k))*idz(k)
              ! compute local timestep
              dt_local(i,j,k)=min(dt_local(i,j,k),CFL/sr_k)
           enddo
        enddo
     enddo

  endif

end subroutine local_timestep

!===============================================================================
subroutine check_max_cfl
!===============================================================================
  !> Computation of maximum CFL
!===============================================================================
  use mod_mpi
  use mod_time  ! <- for deltat
  use mod_flow  ! <- for uu,vv,ww,c_,metrics
  use warnstop
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j,k
  integer :: ip_max(1)
  real(wp) :: val
  real(wp) :: cflmax_i,cflmax_j,cflmax_k
  real(wp), dimension(:), allocatable :: cfl_max
  ! ---------------------------------------------------------------------------
  
  ! Compute max CFL number along i-direction
  ! ========================================
  cflmax_i=0.0_wp
  if (is_curv3) then
     do k=1,nz
        do j=1,ny
           do i=1,nx
              ! calculation of max characteristic speed
              val=abs((uu(i,j,k)*ksi_x(i,j,k)+vv(i,j,k)*ksi_y(i,j,k)+ww(i,j,k)*ksi_z(i,j,k))*ijacob3(i,j,k)) &
                   +c_(i,j,k)*g3_ksi(i,j,k)
              ! compute cflmax_ksi
              cflmax_i=max(cflmax_i,val*deltat)
              !if (val*deltat>100.) then
              !   print *,iproc,i,j,k,ijacob3(i,j,k)
              !   print *,ksi_x(i,j,k),ksi_y(i,j,k),ksi_z(i,j,k)
              !   print *,abs((uu(i,j,k)*ksi_x(i,j,k)+vv(i,j,k)*ksi_y(i,j,k)+ww(i,j,k)*ksi_z(i,j,k))*ijacob3(i,j,k))
              !   print *,c_(i,j,k),g3_ksi(i,j,k)
              !   print *,abs((uu(i,j,k)*ksi_x(i,j,k)+vv(i,j,k)*ksi_y(i,j,k)+ww(i,j,k)*ksi_z(i,j,k))*ijacob3(i,j,k))*deltat
              !   print *,c_(i,j,k)*g3_ksi(i,j,k)*deltat
              !endif
           enddo
        enddo
     enddo
  elseif (is_curv) then
     do k=1,nz
        do j=1,ny
           do i=1,nx
              !print *,i,j,sqrt((y_eta(i,j)**2+x_eta(i,j)**2)*ijacob(i,j)**2)
              ! calculation of max characteristic speed
              val=abs((uu(i,j,k)*y_eta(i,j)-vv(i,j,k)*x_eta(i,j))*ijacob(i,j)) &
                   +c_(i,j,k)*g_ksi(i,j)
              ! compute cflmax_ksi
              cflmax_i=max(cflmax_i,val*deltat)
           enddo
        enddo
     enddo
  else
     do k=1,nz
        do j=1,ny
           do i=1,nx
              cflmax_i=max(cflmax_i,(abs(uu(i,j,k))+c_(i,j,k))*deltat*idx(i))
           enddo
        enddo
     enddo
  endif

  !call mpistop('pause in CFL i',0)

  ! Gather max value on proc 0
  ! --------------------------
  allocate(cfl_max(1:nproc))

  call MPI_GATHER(cflmax_i,1,MPI_DOUBLE_PRECISION,cfl_max,1,MPI_DOUBLE_PRECISION,0,COMM_global,info)
  
  if (iproc==0) then
     ! proc which has max
     ip_max=maxloc(cfl_max)-1
     ! proc 0 print max at screen
     print 100,cfl_max(ip_max+1),nob(ip_max),ip_max
100  format(2x,'maximum CFL in i-direction=',f12.6,' in block ',i4,' (proc ',i5,')')
  endif
 
  ! Print message at screen if cfmax greater than explicit limit
  ! ------------------------------------------------------------
  if ((.not.is_irs_i).and.(cflmax_i>cfl_limit)) &
       call mpiwarn('ATTENTION! max CFL in i-direction is greater than explicit limit of RK!! Check.',0)

  ! Compute max CFL number along j-direction
  ! ========================================
  cflmax_j=0.0_wp
  if (is_curv3) then
     do k=1,nz
        do j=1,ny
           do i=1,nx
              ! calculation of max characteristic speed
              val=abs((uu(i,j,k)*eta_x(i,j,k)+vv(i,j,k)*eta_y(i,j,k)+ww(i,j,k)*eta_z(i,j,k))*ijacob3(i,j,k)) &
                   +c_(i,j,k)*g3_eta(i,j,k)
              ! compute cflmax_eta
              cflmax_j=max(cflmax_j,val*deltat)
           enddo
        enddo
     enddo
  elseif (is_curv) then
     do k=1,nz
        do j=1,ny
           do i=1,nx
              ! calculation of max characteristic speed
              val=abs((vv(i,j,k)*x_ksi(i,j)-uu(i,j,k)*y_ksi(i,j))*ijacob(i,j)) &
                   +c_(i,j,k)*g_eta(i,j)
              ! compute cflmax_eta
              cflmax_j=max(cflmax_j,val*deltat)
           enddo
        enddo
     enddo
  else
     do k=1,nz
        do j=1,ny
           do i=1,nx
              cflmax_j=max(cflmax_j,(abs(vv(i,j,k))+c_(i,j,k))*deltat*idy(j))
           enddo
        enddo
     enddo
  endif

  ! Gather max value on proc 0
  ! --------------------------
  call MPI_GATHER(cflmax_j,1,MPI_DOUBLE_PRECISION,cfl_max,1,MPI_DOUBLE_PRECISION,0,COMM_global,info)
  
  if (iproc==0) then
     ! proc which has max
     ip_max=maxloc(cfl_max)-1
     ! proc 0 print max at screen
     print 101,cfl_max(ip_max+1),nob(ip_max),ip_max
101  format(2x,'maximum CFL in j-direction=',f12.6,' in block ',i4,' (proc ',i5,')')
  endif
 
  ! Print message at screen if cfmax greater than explicit limit
  ! ------------------------------------------------------------
  if ((.not.is_irs_j).and.(cflmax_i>cfl_limit)) &
       call mpiwarn('ATTENTION! max CFL in j-direction is greater than explicit limit of RK!! Check.',0)

  ! Compute max CFL number along k-direction
  ! ========================================
  cflmax_k=0.0_wp
  if (is_curv3) then
     do k=1,nz
        do j=1,ny
           do i=1,nx
              ! calculation of max characteristic speed
              val=abs((uu(i,j,k)*phi_x(i,j,k)+vv(i,j,k)*phi_y(i,j,k)+ww(i,j,k)*phi_z(i,j,k))*ijacob3(i,j,k)) &
                   +c_(i,j,k)*g3_phi(i,j,k)
              ! compute cflmax_eta
              cflmax_k=max(cflmax_k,val*deltat)
           enddo
        enddo
     enddo
  else
     do k=1,nz
        do j=1,ny
           do i=1,nx
              cflmax_k=max(cflmax_k,(abs(ww(i,j,k))+c_(i,j,k))*deltat*idz(k))
           enddo
        enddo
     enddo
  endif

  ! Gather max value on proc 0
  ! --------------------------
  call MPI_GATHER(cflmax_k,1,MPI_DOUBLE_PRECISION,cfl_max,1,MPI_DOUBLE_PRECISION,0,COMM_global,info)

  if (iproc==0) then
     ! proc which has max
     ip_max=maxloc(cfl_max)-1
     ! proc 0 print max at screen
     print 102,cfl_max(ip_max+1),nob(ip_max),ip_max
102  format(2x,'maximum CFL in k-direction=',f12.6,' in block ',i4,' (proc ',i5,')')
  endif

  ! Print message at screen if cfmax greater than explicit limit
  ! ------------------------------------------------------------
  if ((.not.is_irs_k).and.(cflmax_k>cfl_limit)) &
       call mpiwarn('ATTENTION! max CFL in k-direction is greater than explicit limit of RK!! Check.',0)

  deallocate(cfl_max)

end subroutine check_max_cfl

!===============================================================================
subroutine check_max_cfl_interface
!===============================================================================
  !> Computation of maximum CFL in interfaces
!===============================================================================
  use mod_bc
  use mod_time  ! <- for deltat
  use mod_flow  ! <- for uu,vv,ww,c_,metrics
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: i,j,k,ii,jj,ndim,cfldir
  integer :: i1,i2,j1,j2,k1,k2
  real(wp) :: val
  real(wp) :: cflmax_i,cflmax_j,cflmax_k
  ! ---------------------------------------------------------------------------

  ! Dimension
  ! =========
  ndim=3
  if (is_2D) ndim=2

  ! Compute max CFL for each interface (block or MPI)
  ! ==================================
  do ii=1,ndim
     do jj=1,2
        
        ! initialize
        ! ----------
        BC_face(ii,jj)%cflmax=0.0_wp
        
        if (BC_face(ii,jj)%sort>0) then

           ! Fill indices for each faces
           ! ---------------------------
           if (ii==1) then     ! i-direction
              j1=1
              j2=ny
              k1=1
              if (BC_face(3,1)%sort==-8) k1=k1+1
              k2=nz
              if (BC_face(3,2)%sort==-8) k2=k2-1
              if (jj==1) then     ! face imin
                 i1=1
                 i2=1
              elseif (jj==2) then ! face imax
                 i1=nx
                 i2=nx
              endif
           elseif (ii==2) then ! j-direction
              i1=1
              i2=nx
              k1=1
              if (BC_face(3,1)%sort==-8) k1=k1+1
              k2=nz
              if (BC_face(3,2)%sort==-8) k2=k2-1
              if (jj==1) then     ! face jmin
                 j1=1
                 j2=1
              elseif (jj==2) then ! face jmax
                 j1=ny
                 j2=ny
              endif
           elseif (ii==3) then ! k-direction
              i1=1
              i2=nx
              j1=1
              j2=ny
              if (jj==1) then     ! face kmin
                 k1=1
                 k2=1
              elseif (jj==2) then ! face kmax
                 k1=nz
                 k2=nz
              endif
           endif

           ! Compute max CFL number along i-direction
           ! ----------------------------------------
           cflmax_i=0.0_wp
           if (is_curv3) then
              do k=k1,k2
                 do j=j1,j2
                    do i=i1,i2
                       ! calculation of max characteristic speed
                       val=abs((uu(i,j,k)*ksi_x(i,j,k)+vv(i,j,k)*ksi_y(i,j,k)+ww(i,j,k)*ksi_z(i,j,k))*ijacob3(i,j,k)) &
                            +c_(i,j,k)*g3_ksi(i,j,k)
                       ! compute cflmax_ksi
                       cflmax_i=max(cflmax_i,val*deltat)
                    enddo
                 enddo
              enddo
           else
              if (is_curv) then
                 do k=k1,k2
                    do j=j1,j2
                       do i=i1,i2
                          ! calculation of max characteristic speed
                          val=abs((uu(i,j,k)*y_eta(i,j)-vv(i,j,k)*x_eta(i,j))*ijacob(i,j)) &
                               +c_(i,j,k)*g_ksi(i,j)
                          ! compute cflmax_ksi
                          cflmax_i=max(cflmax_i,val*deltat)
                       enddo
                    enddo
                 enddo
              else
                 do k=k1,k2
                    do j=j1,j2
                       do i=i1,i2
                          cflmax_i=max(cflmax_i,(abs(uu(i,j,k))+c_(i,j,k))*deltat*idx(i))
                       enddo
                    enddo
                 enddo
              endif
           endif
           BC_face(ii,jj)%cflmax(1)=cflmax_i

           ! Compute max CFL number along j-direction
           ! ----------------------------------------
           cflmax_j=0.0_wp
           if (is_curv3) then
              do k=k1,k2
                 do j=j1,j2
                    do i=i1,i2
                       ! calculation of max characteristic speed
                       val=abs((uu(i,j,k)*eta_x(i,j,k)+vv(i,j,k)*eta_y(i,j,k)+ww(i,j,k)*eta_z(i,j,k))*ijacob3(i,j,k)) &
                            +c_(i,j,k)*g3_eta(i,j,k)
                       ! compute cflmax_ksi
                       cflmax_j=max(cflmax_j,val*deltat)
                    enddo
                 enddo
              enddo
           else
              if (is_curv) then
                 do k=k1,k2
                    do j=j1,j2
                       do i=i1,i2
                          ! calculation of max characteristic speed
                          val=abs((vv(i,j,k)*x_ksi(i,j)-uu(i,j,k)*y_ksi(i,j))*ijacob(i,j)) &
                               +c_(i,j,k)*g_eta(i,j)
                          ! compute cflmax_eta
                          cflmax_j=max(cflmax_j,val*deltat)
                       enddo
                    enddo
                 enddo
              else
                 do k=k1,k2
                    do j=j1,j2
                       do i=i1,i2
                          cflmax_j=max(cflmax_j,(abs(vv(i,j,k))+c_(i,j,k))*deltat*idy(j))
                       enddo
                    enddo
                 enddo
              endif
           endif
           BC_face(ii,jj)%cflmax(2)=cflmax_j

           ! Compute max CFL number along k-direction
           ! ----------------------------------------
           cflmax_k=0.0_wp
           if (is_curv3) then
              do k=k1,k2
                 do j=j1,j2
                    do i=i1,i2
                       ! calculation of max characteristic speed
                       val=abs((uu(i,j,k)*phi_x(i,j,k)+vv(i,j,k)*phi_y(i,j,k)+ww(i,j,k)*phi_z(i,j,k))*ijacob3(i,j,k)) &
                            +c_(i,j,k)*g3_phi(i,j,k)
                       ! compute cflmax_ksi
                       cflmax_k=max(cflmax_k,val*deltat)
                    enddo
                 enddo
              enddo
           else
              if (.not.is_2D) then
                 do k=k1,k2
                    do j=j1,j2
                       do i=i1,i2
                          cflmax_k=max(cflmax_k,(abs(ww(i,j,k))+c_(i,j,k))*deltat*idz(k))
                       enddo
                    enddo
                 enddo
              endif
           endif
           BC_face(ii,jj)%cflmax(3)=cflmax_k

           ! Select max CFL among directions (and save direction of max)
           ! -----------------------------------------------------------
           cfldir=1
           if (cflmax_j>cflmax_i) cfldir=2
           if (cflmax_k>cflmax_j) cfldir=3

           ! Print info at screen
           ! --------------------
!!$           if (cfldir==1) then
!!$              print 101,ii,jj,BC_face(ii,jj)%cflmax(1)
!!$              if (is_irs_i) print 104,int(2*BC_face(ii,jj)%cflmax(1)+1)
!!$           elseif (cfldir==2) then
!!$              print 102,ii,jj,BC_face(ii,jj)%cflmax(2)
!!$              if (is_irs_j) print 104,int(2*BC_face(ii,jj)%cflmax(2)+1)
!!$           elseif (cfldir==3) then
!!$              print 103,ii,jj,BC_face(ii,jj)%cflmax(3)
!!$              if (is_irs_k) print 104,int(2*BC_face(ii,jj)%cflmax(3)+1)
!!$           endif
101  format(2x,'max CFL in face(',i1,',',i1,')=',f12.6,' in i-direction')
102  format(2x,'max CFL in face(',i1,',',i1,')=',f12.6,' in j-direction')
103  format(2x,'max CFL in face(',i1,',',i1,')=',f12.6,' in k-direction')
104  format(2x,'[IRS hint: interface overlap of 2CFL+1 should be',i3,' points]')
           
        endif
     enddo
  enddo

end subroutine check_max_cfl_interface
