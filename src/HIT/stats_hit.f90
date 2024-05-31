!===============================================================================
subroutine stats_hit
!===============================================================================
  !> Compute some stats for HIT
!===============================================================================
  use mod_mpi
  use mod_constant
  use mod_flow
  use mod_time
  use, intrinsic :: iso_c_binding ! C++ interface for FFTW routines
  implicit none
  ! ----------------------------------------------------------------------------
  include 'fftw3-mpi.f03'
  ! ----------------------------------------------------------------------------
  integer :: i,j,k,n_shell
  integer :: ii,jj,kk
  real(wp) :: fac,fac2
  real(wp) :: sctmp
  real(wp), dimension(:), allocatable :: e_spec,e_spec1

  real(wp) :: energy,eps_v,eta,etakmax,re_lambda,u_var,x_length
  real(wp) :: lambda,tau_e!,enstrophy
  real(wp) :: nu
  ! ----------------------------------------------------------------------------
  ! wavenumber vectors ("a" added for real arrays)
  integer :: kmax
  real(wp), dimension(:), allocatable :: akx,aky,akz
  real(wp) :: ak,ak2 ! wavenumber norm, squared
  ! ----------------------------------------------------------------------------
  ! FFTW variables
  ! ==============
  type(C_PTR) :: plan_f_x,plan_f_y,plan_f_z ! FFTW planes for forward transforms
  type(C_PTR) :: cus,cvs,cws ! pointers on velocity data
  complex(C_DOUBLE_COMPLEX), dimension(:,:,:), pointer :: us,vs,ws ! velocity data
  ! ----------------------------------------------------------------------------

  ! Compute velocity spectra
  ! ========================

  ! Init FFTW routines
  ! ------------------
  call x_fftw_init

  ! Compute FFT of velocity field
  ! -----------------------------
  us=uu(1:nx,1:ny,1:nz)
  vs=vv(1:nx,1:ny,1:nz)
  ws=ww(1:nx,1:ny,1:nz)
  call x_fftw3d_f(us) 
  call x_fftw3d_f(vs) 
  call x_fftw3d_f(ws) 
  
  ! max wavenumber
  ! --------------
  kmax=ngx/2
  
  ! Compute non-dimensional wavenumbers
  ! -----------------------------------
  allocate(akx(nx),aky(ny),akz(nz))
  do i=1,nx
     ii=i+coord(1)*nx-1
     if (ii.le.ngx/2) then
        akx(i)=dble(ii)
     else
        akx(i)=dble(ii-ngx)
     endif
  enddo
  do j=1,ny
     jj=j+coord(2)*ny-1
     if (jj.le.ngy/2) then
        aky(j)=dble(jj)
     else
        aky(j)=dble(jj-ngy)
     endif
  enddo
  do k=1,nz
     kk=k+coord(3)*nz-1
     if (kk.le.ngz/2) then
        akz(k)=dble(kk)
     else
        akz(k)=dble(kk-ngz)
     endif
  enddo

  ! Normalization factor (because the FFT is unnormalized)
  ! --------------------
  fac=0.5_wp/real(ngx*ngy*ngz)**2
  ! [0.5 because kinetic energy is 0.5*ui^2]

  ! Allocate energy spectra
  ! -------------------------
  allocate(e_spec(kmax),e_spec1(kmax))

  ! Assembling the total energy in each shell
  ! -----------------------------------------
  e_spec=0.0_wp
  do k=1,nz
     do j=1,ny
        do i=1,nx
           ak2=akx(i)**2+aky(j)**2+akz(k)**2
           n_shell=nint(sqrt(ak2))
           if ((n_shell.gt.0).and.(n_shell.le.kmax)) then
              fac2=fac*(abs(us(i,j,k))**2+abs(vs(i,j,k))**2+abs(ws(i,j,k))**2)
              e_spec(n_shell)=e_spec(n_shell)+fac2
           endif
        enddo
     enddo
  enddo
  call MPI_ALLREDUCE(MPI_IN_PLACE,e_spec,kmax,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
  e_spec=e_spec*L_ref

  ! Compute dissipation spectrum
  ! ----------------------------
  nu=mu_ref/rho_ref
  do k=1,kmax
     ak=dble(k)/L_ref
     e_spec1(k)=2.0_wp*nu*e_spec(k)*ak**2
  end do

  ! Write spectra [-> in es.dat opened in stats.f90]
  ! =============

  ! outputting the energy spectrum
  !open(900,file='es.dat',position='append')
  !write(900,"()")
  !write(900,"()")
  !write(900,"('# ITIME=',i7,' TIME=',e17.8)") ITIME, TIME
  if (iproc==0) then
     do k=1,kmax
        write(100,*) k/L_ref,e_spec(k),e_spec1(k)
     enddo
  endif
  
  ! Compute velocity stats
  ! ======================

  ! Total energy
  ! ------------
  energy=sum(e_spec(1:kmax))/L_ref

  ! Total dissipation
  ! -----------------
  nu=mu_ref/rho_ref
  do k=1,kmax
     ak=dble(k)/L_ref
     e_spec1(k)=2.0_wp*nu*e_spec(k)*ak**2
  end do
  eps_v=sum(e_spec1(1:kmax))/L_ref

  ! Kolmogorov scale
  ! ----------------
  eta=(nu**3/eps_v)**0.25
  etakmax=eta*dble(kmax)/L_ref

  ! Variance
  ! --------
  u_var=2.0_wp/3.0_wp*energy

  ! Integral length scale
  ! ---------------------
  sctmp=0.0_wp
  do k=1,kmax
     ak=dble(k)/L_ref
     sctmp=sctmp+e_spec(k)/ak
  enddo
  x_length=pi/2.0_wp*sctmp/u_var/L_ref

  ! Taylor microscale
  ! -----------------
  lambda=sqrt(15.0_wp*u_var*nu/eps_v)

  ! Taylor-Reynolds number
  ! ----------------------
  !re_lambda =u_var*sqrt(15.0_wp/eps_v*RE)
  re_lambda=sqrt(u_var)*lambda/nu

  ! Eddy turnover time
  ! ------------------
  tau_e=x_length/sqrt(u_var)

  ! Write velocity stats
  ! ====================

  ! outputting all this in the stat2 file
!!$  inquire(file='stat2.dat', exist=there, opened=there2)
!!$  if (.not.there) then
!!$     open(70,file='stat2.dat',form='formatted')
!!$     write(70,'(A)') '# 1.itime  2.time         3.int LS       4. lambda      5.R_lambda1    6.tau_e        7.etakmax'
!!$  end if
!!$  if(there.and..not.there2) then
!!$     open(70,file='stat2.dat',position='append')
!!$  end if

  if (iproc==0) then

     ! Print at screen
     ! ---------------
     print *,'   >> energy int(E(k),dk)',energy
     print *,'   >> dissipation 2*nu*int(k^2*E(k),dk)',eps_v,eps_v/(1e2)*(1e2)**3
     print *,'   >> Kolmogorov scale (cm)',eta*1e2
     print *,'   >> urms (cm/s)',sqrt(u_var)*1e2
     print *,'   >> integral scale (cm)',x_length*1e2
     print *,'   >> Taylor scale (cm)',lambda*1e2
     print *,'   >> Re_Taylor',re_lambda
     print *,'   >> Eddy turnover time',tau_e
     
     ! Write in file stat1.dat [opened in stats.f90]
     ! ---------------------------------------------
     write(199,'(10(f15.10,1x))') tstar,energy,eps_v,eta,u_var,x_length,lambda,re_lambda,tau_e,etakmax

  endif

  ! FFTW: Free planes and memory from malloc
  ! ========================================
  call x_fftw_exit
  
contains

  !===============================================================================
  subroutine x_fftw_init
  !===============================================================================
    !> Stats HIT: subroutine that initializes the auxiliary arrays for FFT
  !===============================================================================
    use mod_mpi_part
    implicit none
    ! ----------------------------------------------------------------------------
    integer(C_INTPTR_T) :: l_nx,l_ny,l_nz ! local dimensions
    integer(C_INTPTR_T) :: g_nx,g_ny,g_nz ! global dimensions
    integer(C_INTPTR_T) :: alloc_size
    ! ----------------------------------------------------------------------------

    ! Local and global dimensions in C format
    ! =======================================
    g_nx=ngx
    g_ny=ngy
    g_nz=ngz
    l_nx=nx
    l_ny=ny
    l_nz=nz

    ! init MPI for FFTW
    ! =================
    call fftw_mpi_init()

    ! Get local data size and allocate arrays
    ! =======================================
    alloc_size=l_nx*l_ny*l_nz
    cus=fftw_alloc_complex(alloc_size)
    call c_f_pointer(cus,us,[l_nx,l_ny,l_nz]) ! pointer on us
    cvs=fftw_alloc_complex(alloc_size)
    call c_f_pointer(cvs,vs,[l_nx,l_ny,l_nz]) ! pointer on vs
    cws=fftw_alloc_complex(alloc_size)
    call c_f_pointer(cws,ws,[l_nx,l_ny,l_nz]) ! pointer on ws

    ! Create MPI planes for in-place forward DFT
    ! ==========================================
    plan_f_x=fftw_mpi_plan_dft_1d(g_nx,us(:,1,1),us(:,1,1),COMMYZ,FFTW_FORWARD,FFTW_ESTIMATE)
    plan_f_y=fftw_mpi_plan_dft_1d(g_ny,us(1,:,1),us(1,:,1),COMMXZ,FFTW_FORWARD,FFTW_ESTIMATE)
    plan_f_z=fftw_mpi_plan_dft_1d(g_nz,us(1,1,:),us(1,1,:),COMMXY,FFTW_FORWARD,FFTW_ESTIMATE)

  end subroutine x_fftw_init

  !===============================================================================
  subroutine x_fftw3d_f(var)
  !===============================================================================
    !> Stats HIT: 3x1D (forward) Fourier transforms
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k
    complex(C_DOUBLE_COMPLEX), dimension(:,:,:), pointer :: var
    ! ----------------------------------------------------------------------------

    do k=1,nz
       do j=1,ny
          call fftw_mpi_execute_dft(plan_f_x,var(:,j,k),var(:,j,k))
       enddo
    enddo
    do k=1,nz
       do i=1,nx
          call fftw_mpi_execute_dft(plan_f_y,var(i,:,k),var(i,:,k))
       enddo
    enddo
    do i=1,nx
       do j=1,ny
          call fftw_mpi_execute_dft(plan_f_z,var(i,j,:),var(i,j,:))
       enddo
    enddo

  end subroutine x_fftw3d_f

  !===============================================================================
  subroutine x_fftw_exit
  !===============================================================================
    !> Stats HIT: subroutine to free FFTW planes and memory
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    ! Destroy FFTW planes
    ! ===================
    call fftw_destroy_plan(plan_f_x)
    call fftw_destroy_plan(plan_f_y)
    call fftw_destroy_plan(plan_f_z)

    ! Free memory from malloc
    ! =======================
    call fftw_free(cus)
    call fftw_free(cvs)
    call fftw_free(cws)

  end subroutine x_fftw_exit

end subroutine stats_hit
