!=================================================================================
module mod_init_hit
!=================================================================================
  !> Module to initialize HIT (Homogeneous Isotropic Turbulence)
!=================================================================================
  use mod_mpi
  use mod_constant ! for diffscale TO BE CHANGED ??
  use mod_flow
  use warnstop
  use, intrinsic :: iso_c_binding ! C++ interface for FFTW routines
  implicit none
  !-------------------------------------------------------------------------------
  include 'fftw3-mpi.f03'
  !-------------------------------------------------------------------------------
  ! HIT parameters [read in param_hit.ini]
  ! ==============
  ! Initialization methods
  ! ----------------------
  ! Choice of method to generate velocities
  character(len=2) :: init_vel
  ! Compute solenoidal pressure by solving Poisson problem
  logical :: is_poisson
  ! Compressible case: init non-zero divergence part of velocities
  real(wp) :: chi0 ! Compressibility ratio
  logical :: is_compr
  ! Energy spectrum model
  ! ---------------------
  character(len=2) :: spectr_model ! Spectrum model
  ! ['Ko':Kolmogorov inertial range; 'PP':Passot-Pouquet; 'Ex': exponential sp; 'VK'; von Karman]
  character(len=1) :: dissip_range ! Dissipation range ['N':none; 'P':Pao; 'S':Saffman]
  character(len=1) :: bottleneck   ! Bottleneck correction of Kang, Chester & Meneveau ['N':none; 'B': bottleneck]
  real(wp) :: ak0     ! Peak wavenumber [for Passot-Pouquet or exponential spectra]
  real(wp) :: Lf_int  ! Integral length scale [for von Karman spectrum]
  real(wp) :: eta_kol ! Kolmogorov microscale [for dissipation range or bottleneck correction]
  real(wp) :: eps_v, tau_e   ! Dissipation rate, eddy time turnover
  !-------------------------------------------------------------------------------
  ! FFTW variables
  ! ==============
  type(C_PTR), private :: plan_b_x,plan_b_y,plan_b_z ! FFTW planes for backward transforms
  type(C_PTR), private :: plan_f_x,plan_f_y,plan_f_z ! FFTW planes for forward transforms
  type(C_PTR), private :: cus,cvs,cws ! pointers on velocity data
  complex(C_DOUBLE_COMPLEX), dimension(:,:,:), pointer, private :: us,vs,ws ! velocity data
  type(C_PTR), private :: cus_x,cvs_x,cws_x ! pointers on x-derivatives
  type(C_PTR), private :: cus_y,cvs_y,cws_y ! pointers on y-derivatives
  type(C_PTR), private :: cus_z,cvs_z,cws_z ! pointers on z-derivatives
  complex(C_DOUBLE_COMPLEX), dimension(:,:,:), pointer, private :: us_x,vs_x,ws_x ! x-derivatives
  complex(C_DOUBLE_COMPLEX), dimension(:,:,:), pointer, private :: us_y,vs_y,ws_y ! y-derivatives
  complex(C_DOUBLE_COMPLEX), dimension(:,:,:), pointer, private :: us_z,vs_z,ws_z ! z-derivatives
  type(C_PTR) :: cps ! pointers on pressure & RHS
  complex(C_DOUBLE_COMPLEX), dimension(:,:,:), pointer, private :: ps ! pressure & RHS for Poisson problem
  !-------------------------------------------------------------------------------  
  ! constants
  real(wp), private :: two_pi
  complex(wp), parameter, private :: i_=(0.0_wp,1.0_wp)
  real(wp), private :: norm ! normalization factor
  !-------------------------------------------------------------------------------  
  ! wavenumber vectors ("a" added for real arrays)
  integer, private :: kmax
  real(wp), dimension(:), allocatable, private :: akx,aky,akz
  real(wp), private :: ak,ak2 ! wavenumber norm, squared
  !-------------------------------------------------------------------------------  
  ! random number generator
  real(wp), private :: rseed
  !-------------------------------------------------------------------------------  
  real(wp) :: Uref
  !-------------------------------------------------------------------------------  

contains
  
  !===============================================================================
  subroutine read_param_hit
  !===============================================================================
    !> Read file for HIT parameters
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    logical :: iexist
    ! ----------------------------------------------------------------------------

    if (iproc==0) then
       print *,repeat('=',70)
       print *,'HIT initialization'
       print *,repeat('=',70)
    endif

    ! constant
    two_pi=2.0_wp*pi

    ! New reference velocity TO BE CHANGED
    ! ======================
    Uref=u_ref/sqrt(2.0_wp)
           
    ! Read param_hit.ini
    ! ===================
    inquire(file='param_hit.ini', exist=iexist)
    if (.not.iexist) then
       call mpistop('Paramfile param_hit.ini does not exist!', 0)
    endif

    !=============================================================================
    open(30,file=trim('param_hit.ini'))
    !=============================================================================
    read(30,*)! ==================================================================
    read(30,*)! HOMOGENEOUS ISOTROPIC TURBULENCE (HIT) PARAMETERS
    read(30,*)! ==================================================================
    read(30,*)! Initialization methods
    read(30,*)! ==================================================================
    read(30,*)! Choice of method to generate velocities
    read(30,*)! 'C': version adapted from HIT3D code (Chumakov et al.): incompressible
    read(30,*)! 'PG': version adapted from Pirozzoli and Grasso: compressible or not
    read(30,*) init_vel
    read(30,*)! Compute solenoidal pressure by solving Poisson problem 
    read(30,*) is_poisson
    read(30,*)! Compressible case: init non-zero divergence part of velocities
    read(30,*)! Compressibility ratio [compressible if chi0.ne.0]
    read(30,*) chi0
    read(30,*)! ==================================================================
    read(30,*)! Energy spectrum model
    read(30,*)! ==================================================================
    read(30,*)! Spectrum model: ['Ko':Kolmogorov inertial range; 
    read(30,*)! 'PP':Passot-Pouquet; 'Ex': exponential sp; 'VK'; von Karman]
    read(30,*) spectr_model
    read(30,*)! Define viscous dissipation range ['N':none; 'P':Pao; 'S':Saffman]
    read(30,*) dissip_range
    read(30,*)! Bottleneck correction of Kang, Chester & Meneveau ['N':none; 'B': bottleneck]
    read(30,*) bottleneck
    read(30,*)! Peak Wavenumber: ak0 [for Passot-Pouquet or exponential spectra]
    read(30,*) ak0
    read(30,*)! Integral length scale: Lf_int [for von Karman spectrum]
    read(30,*) Lf_int
    ! non-dimensionalize integral length scale
    Lf_int=Lf_int/L_ref
    ! if (L_ref==1.0_wp) Lf_int=Lf_int/0.087318767977937 ! /Lref=5.08e-2*10.8/(2pi)
    read(30,*)! Kolmogorov microscale: eta_kol [for dissipation range or bottleneck correction] 
    read(30,*) eta_kol
    ! non-dimensionalize Kolmogorov's scale
    eta_kol=eta_kol/L_ref
    ! if (L_ref==1.0_wp) eta_kol=eta_kol/0.087318767977937 ! /Lref=5.08e-2*10.8/(2pi)
    close(30)
    !=============================================================================

    ! Check parameters
    ! ================
    if (trim(init_vel)=='C') then
       is_compr=.false.
       if (chi0.ne.0.0_wp) then
          if (iproc==0) print *,'WARNING: no compressibility initialization for this method of velocity generation'
          chi0=0.0_wp
       endif
    elseif (trim(init_vel)=='PG') then
       if (chi0.gt.1.e-3_wp) then
          is_compr=.true.
       else
          is_compr=.false.
          if (iproc==0) print *,'   only initialization of solenoidal velocity'
       endif
    else
       call mpistop('Wrong method to generate HIT velocities !', 0)
    endif
    
    ! Check compatibility of spectrum options
    ! =======================================
    if (spectr_model.ne.'VK') then
       if ((dissip_range.ne.'N').and.(iproc==0)) &
          print *,'WARNING: dissipation model is ignored for this choice of spectrum model'
       if ((bottleneck.ne.'N').and.(iproc==0)) &
          print *,'WARNING: bottleneck correction is ignored for this choice of spectrum model'
    endif

    ! Print spectrum model at screen
    ! ==============================
    if (spectr_model=='Ko') then
       if (iproc==0) print *,'   Energy spectrum model: Kolmogorov (inertial range)'
    elseif (spectr_model=='PP') then
       if (iproc==0) print *,'   Energy spectrum model: Passot-Pouquet'
    elseif (spectr_model=='Ex') then
       if (iproc==0) print *,'   Energy spectrum model: exponential'
    elseif (spectr_model=='VK') then       
       if (dissip_range=='N') then
          if (bottleneck=='N') then
             if (iproc==0) print *,'   Spectrum model: von Karman'
          elseif (bottleneck=='B') then
             if (iproc==0) print *,'   Spectrum model: von Karman + bottleneck'
          else
             call mpistop('Wrong bottleneck model !', 0)
          endif
       elseif (dissip_range=='P') then
          if (bottleneck=='N') then
             if (iproc==0) print *,'   Spectrum model: von Karman-Pao'
          elseif (bottleneck=='B') then
             if (iproc==0) print *,'   Spectrum model: von Karman-Pao + bottleneck'
          else
             call mpistop('Wrong bottleneck model !', 0)
          endif
       elseif (dissip_range=='S') then
          if (bottleneck=='N') then
             if (iproc==0) print *,'   Spectrum model: von Karman-Saffman'
          elseif (bottleneck=='B') then
             if (iproc==0) print *,'   Spectrum model: von Karman-Saffman + bottleneck'
          else
             call mpistop('Wrong bottleneck model !', 0)
          endif
       else
          call mpistop('Wrong dissipation range model !', 0)
       endif       
    else
       call mpistop('Wrong spectrum model !', 0)
    endif

  end subroutine read_param_hit

  !===============================================================================
  subroutine init_hit
  !===============================================================================
    !> Initialization for CHIT case (MPI version)
  !===============================================================================
    use mod_io        ! to write restart file
    use mod_interface ! for grad_velT & communications TO BE CHANGED
    use mod_deriv
    use mod_eos       ! for tcal_pro & ecal_tro
    use mod_time      ! for tscale in initialization part TO BE CHANGED
    implicit none
    ! ----------------------------------------------------------------------------
    ! loop counters
    integer :: i,j,k,ii,jj,kk
    ! constants
    logical :: iexist
    ! variables for normalization
    real(wp) :: akmax2
    real(wp) :: sfac
    real(wp) :: uu1,vv1,ww1,dudx2m
    real(wp) :: urms,vrms,wrms,prms,qq,max_rms
    ! ----------------------------------------------------------------------------
    real(wp) :: L_taylor,Atau,tau_eddy
    ! ----------------------------------------------------------------------------

    ! Parameters settings
    ! ===================
    call read_param_hit
    
    ! FFT initializations
    ! ===================
    
    ! normalization factor for FFT
    ! ----------------------------
    norm=1.0_wp/dble(ngx*ngy*ngz)

    ! max wavenumber
    ! --------------
    kmax=ngx/2
    if (iproc==0) write(6,'(''    kmax='',i5)') kmax

    ! init FFTW routines
    ! ------------------
    call x_fftw_init
    
    ! Compute non-dimensional wavenumbers
    ! ===================================
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

    ! Compute non-dimensional velocity components in Fourier space
    ! ============================================================
    call init_velocity_hit
    
    ! Check velocity stats
    ! ====================
    call stat_velocity

    ! Compute derivatives of velocity components in Fourier space
    ! ===========================================================
    do i=1,nx
       us_x(i,:,:)=i_*akx(i)*us(i,:,:)
       vs_x(i,:,:)=i_*akx(i)*vs(i,:,:)
       ws_x(i,:,:)=i_*akx(i)*ws(i,:,:)
    enddo
    do j=1,ny
       us_y(:,j,:)=i_*aky(j)*us(:,j,:)
       vs_y(:,j,:)=i_*aky(j)*vs(:,j,:)
       ws_y(:,j,:)=i_*aky(j)*ws(:,j,:)
    enddo
    do k=1,nz
       us_z(:,:,k)=i_*akz(k)*us(:,:,k)
       vs_z(:,:,k)=i_*akz(k)*vs(:,:,k)
       ws_z(:,:,k)=i_*akz(k)*ws(:,:,k)
    enddo

    ! Mode truncation to ensure isotropy (for all modes that are higher than kmax)
    ! ==================================
    akmax2=dble(kmax)**2
    do k=1,nz
       do j=1,ny
          do i=1,nx
             ak2=akx(i)**2+aky(j)**2+akz(k)**2
             if (ak2.gt.akmax2) then
                us(i,j,k)=0.0_wp
                vs(i,j,k)=0.0_wp
                ws(i,j,k)=0.0_wp
             endif
          enddo
       enddo
    enddo
   
    ! Compute inverse Fourier transforms of velocity and derivatives
    ! ==============================================================    
    call x_fftw3d_b(us) 
    call x_fftw3d_b(vs) 
    call x_fftw3d_b(ws)
    call x_fftw3d_b(us_x) 
    call x_fftw3d_b(vs_x) 
    call x_fftw3d_b(ws_x) 
    call x_fftw3d_b(us_y) 
    call x_fftw3d_b(vs_y) 
    call x_fftw3d_b(ws_y) 
    call x_fftw3d_b(us_z) 
    call x_fftw3d_b(vs_z) 
    call x_fftw3d_b(ws_z)
        
    ! Rescaling velocity components
    ! =============================

    ! Compute q=sqrt(<u'^2+v'^2+w'^2>), rms components and max fluctuations
    ! ---------------------------------------------------------------------
    urms=0.0_wp
    vrms=0.0_wp
    wrms=0.0_wp
    qq=0.0_wp
    max_rms=0.0_wp
    do k=1,nz
       do j=1,ny
          do i=1,nx
             uu1=real(us(i,j,k))
             vv1=real(vs(i,j,k))
             ww1=real(ws(i,j,k))
             sfac=uu1**2+vv1**2+ww1**2
             qq=qq+sfac
             max_rms=max(max_rms,sqrt(sfac))
             urms=urms+uu1**2
             vrms=vrms+vv1**2
             wrms=wrms+ww1**2
          enddo
       enddo
    enddo
    call MPI_ALLREDUCE(MPI_IN_PLACE,urms,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
    call MPI_ALLREDUCE(MPI_IN_PLACE,vrms,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
    call MPI_ALLREDUCE(MPI_IN_PLACE,wrms,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
    call MPI_ALLREDUCE(MPI_IN_PLACE, qq ,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
    call MPI_ALLREDUCE(MPI_IN_PLACE, max_rms,1,MPI_DOUBLE_PRECISION,MPI_MAX,COMM_global,info)
    urms=sqrt(urms*norm)
    vrms=sqrt(vrms*norm)
    wrms=sqrt(wrms*norm)
    qq=sqrt(qq*norm)
    
    ! Rescaling factor [u_ref=Mach_turbulent*c_ref=q=sqrt(3.)*u_rms=sqrt(2.*K)]
    ! ----------------
    if (trim(init_vel)=='C') then
       sfac=Uref    
    elseif (trim(init_vel)=='PG') then
       sfac = sqrt(1.0_wp-chi0)*u_ref/qq
    endif
    if (iproc==0) write(6,'(''    [rescaling factor:'',F8.2,F8.2,'']'')') sfac,Uref
    
    ! Rescale variables
    ! -----------------
    ! rms velocities
    qq=qq*sfac
    max_rms=max_rms*sfac
    urms=urms*sfac
    vrms=vrms*sfac
    wrms=wrms*sfac    
    ! velocity components
    us=real(us)*sfac
    vs=real(vs)*sfac
    ws=real(ws)*sfac
    ! velocity derivatives
    us_x=real(us_x)*sfac
    vs_x=real(vs_x)*sfac
    ws_x=real(ws_x)*sfac
    us_y=real(us_y)*sfac
    vs_y=real(vs_y)*sfac
    ws_y=real(ws_y)*sfac
    us_z=real(us_z)*sfac
    vs_z=real(vs_z)*sfac
    ws_z=real(ws_z)*sfac
    
    ! Compute dudx2m=<(du/dx)^2+(dv/dy)^2+(dw/dz)^2> in Fourier space
    ! ---------------------------------------------------------------
    dudx2m=0.0_wp
    do k=1,nz
       do j=1,ny
          do i=1,nx
             dudx2m=dudx2m+real(us_x(i,j,k))**2+real(vs_y(i,j,k))**2+real(ws_z(i,j,k))**2
          enddo
       enddo
    enddo
    call MPI_ALLREDUCE(MPI_IN_PLACE,dudx2m,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
    dudx2m=dudx2m*norm/L_ref**2
    
    ! Taylor's microscale
    ! -------------------
    L_taylor=qq/sqrt(dudx2m)
    
    ! Determine solenoidal pressure (solving the Poisson problem in dimensional space)
    ! =============================
    if (is_poisson) then

       call init_pressure_hit
       
       ! check rms pressure
       prms=0.0_wp
       do k=1,nz
          do j=1,ny
             do i=1,nx
                prms =prms+real(ps(i,j,k))**2
             enddo
          enddo
       enddo
       call MPI_ALLREDUCE(MPI_IN_PLACE,prms,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
       prms=sqrt(prms*norm)
       
    endif
 
    ! Check urms
    ! ----------    
    if (iproc==0) then
       write(6,'(''    u_rms='',F8.4)') urms
       write(6,'(''    v_rms='',F8.4)') vrms
       write(6,'(''    w_rms='',F8.4)') wrms
       if (is_poisson) &
       write(6,'(''    p_rms/u_rms^2='',F8.4)') prms/urms**2
       write(6,'(''    q=<ui^2>='',F8.4,''/ q[target]='',F8.4)') qq,u_ref
       write(6,'(''    M_t='',F10.6,''/ M_t[target]='',F10.6)') qq/c_ref,Mach
       write(6,'(''    max fluct (max_rms)='',F8.4)') max_rms
    endif  

    open(51,file='vites2.bin',form='unformatted',status='unknown')
    rewind(51)
    write(51) nx
    write(51) ny
    write(51) nz
    write(51) ((real(us(i,j,nz/2)),i=1,nx),j=1,ny)
    write(51) ((real(vs(i,j,nz/2)),i=1,nx),j=1,ny)
    if (is_poisson) then
       write(51) ((real(ps(i,j,nz/2)),i=1,nx),j=1,ny)
    else
       write(51) ((real(ws(i,j,nz/2)),i=1,nx),j=1,ny)
    endif
    write(51) ((real(us(i,ny/2,k)),i=1,nx),k=1,nz)
    write(51) ((real(vs(i,ny/2,k)),i=1,nx),k=1,nz)
    if (is_poisson) then
       write(51) ((real(ps(i,ny/2,k)),i=1,nx),k=1,nz)
    else
       write(51) ((real(ws(i,ny/2,k)),i=1,nx),k=1,nz)
    endif
    write(51) ((real(us(nx/2,j,k)),j=1,ny),k=1,nz)
    write(51) ((real(vs(nx/2,j,k)),j=1,ny),k=1,nz)
    if (is_poisson) then
       write(51) ((real(ps(nx/2,j,k)),j=1,ny),k=1,nz)
    else
       write(51) ((real(ws(nx/2,j,k)),j=1,ny),k=1,nz)
    endif
    close(51)

    ! Fill conservative variables and compute max turbulent velocity
    ! ==============================================================
    if (iproc==0) print *,'~> compute conservative variables'
    rho=rho_ref
    do k=1,nz
       do j=1,ny
          do i=1,nx
             uu1=real(us(i,j,k))
             vv1=real(vs(i,j,k))
             ww1=real(ws(i,j,k))
             if (is_poisson) then
                prs(i,j,k)=p_ref+real(ps(i,j,k))
             else
                prs(i,j,k)=p_ref
             endif
             ! momenta
             rhou(i,j,k)=rho_ref*uu1
             rhov(i,j,k)=rho_ref*vv1
             rhow(i,j,k)=rho_ref*ww1
             ! temperature
             Tmp(i,j,k) =tcalc_pro(prs(i,j,k),rho(i,j,k),T_ref)
             ! total energy
             rhoe(i,j,k)= rho(i,j,k)*(ecalc_tro(Tmp(i,j,k),rho(i,j,k)) &
                        + 0.5_wp*(uu1**2+vv1**2+ww1**2))
          enddo
       enddo
    enddo

    ! Evaluate Taylor's microscale from velocity gradients
    ! ====================================================
    if (iproc==0) print *,'~> compute Taylor''s microscale ..'

    ! Values computed in Fourier space
    ! --------------------------------
    if (iproc==0) then
       write(6,'(''    * in Fourier space *'')') 
       write(6,'(''    <(dui/dxi)^2>='',F8.2)') dudx2m
       write(6,'(''    computed Taylor''''s microscale:'',F8.4,'' cm'')') L_taylor*1.e2_wp
       write(6,'(''    Taylor''''s Reynolds number:'',F12.2)') rho_ref*L_taylor*u_ref/sqrt(3.0_wp)/muw_ref
    endif
       
    ! Compute dudx2m=<(du/dx)^2+(dv/dy)^2+(dw/dz)^2> in physical space
    ! ----------------------------------------------------------------
    ! fill ghost points
    call communication_(rho,rhou,rhov,rhow,rhoe)
    ! velocity components
    uu = rhou/rho
    vv = rhov/rho
    ww = rhow/rho
    ! compute velocity derivatives
    call deriv_x_11pts(uu,dux)
    call deriv_y_11pts(vv,dvy)
    call deriv_z_11pts(ww,dwz)
    ! compute dudx2m
    dudx2m=0.0_wp
    do k=1,nz
       do j=1,ny
          do i=1,nx
             dudx2m=dudx2m+dux(i,j,k)**2+dvy(i,j,k)**2+dwz(i,j,k)**2
          enddo
       enddo
    enddo
    call MPI_ALLREDUCE(MPI_IN_PLACE,dudx2m,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
    dudx2m=dudx2m*norm
    ! Taylor's microscale
    ! -------------------
    L_taylor=qq/sqrt(dudx2m)
    
    if (iproc==0) then
       write(6,'(''    * in physical space *'')') 
       write(6,'(''    <(dui/dxi)^2>='',F8.2)') dudx2m
       write(6,'(''    computed Taylor''''s microscale:'',F8.4,'' cm'')') L_taylor*1.e2_wp
       write(6,'(''    Taylor''''s Reynolds number:'',F12.2)') rho_ref*L_taylor*u_ref/sqrt(3.0_wp)/muw_ref
    endif
   
    ! Rescaling viscosity
    ! ===================
    if (L_ref==1.0_wp) then
       ! if L_ref is 2*pi, we need to rescale viscosity
       if (iproc==0) print *,'~> rescale viscosity to get the correct Reynolds number ..'
       ! The Reynolds number being based on Taylor's microscale:
       ! Re_ref=rho_ref*L_taylor*u_ref/sqrt(3.)/muw_ref
       !       =rho_ref*L_taylor*q/sqrt(3.)/muw_ref
       mu_ref = rho_ref*L_taylor*qq/sqrt(3.0_wp)/Re_ref
    else
       mu_ref = rho_ref*L_ref*qq/Re_ref
       mu_ref = rho_ref*L_ref*u_ref/Re_ref
    endif
    
    ! Rescaling factor for viscosity
    ! ------------------------------
    diffscale=mu_ref/muw_ref
    if (iproc==0) write(6,'(''    Rescaling factor for viscosity:'',F19.12)') diffscale

    ! Compute time scale
    ! ==================
    
    if (spectr_model=='PP') then
       ! for Passot-Pouquet spectrum:
       ! rescaling is applied to obtain the target Re_lambda
       ! and an analytical formula can be used to obtain the eddy turnover time
       ! ~> it is used as time scale for the simulation
       
       ! Compute eddy turnover time
       ! --------------------------
       qq=qq**2
       Atau = 64.0_wp*(0.5_wp*qq)/(3.0_wp*sqrt(two_pi)*ak0**5)
       tau_eddy = sqrt(32.0_wp/Atau)*(two_pi)**0.25*ak0**(-3.5_wp)
       if (iproc==0) write(6,'(''    tau_eddy:'',F19.12)') tau_eddy
       
       ! non-dimensional time scale
       ! --------------------------
       tscale = tau_eddy
       
    elseif (spectr_model=='VK') then
       ! -------------------------------
       !       CBC initialisation
       ! -------------------------------
       ! for CBC experiments, the time scale is based on non-dimensional time M/U_flow
       ! check dimensional quantities for CBC
       ! ------------------------------------
       call stat_velocity_dim
       ! ! non-dimensional time scale
       ! ! --------------------------
       ! tscale=5.08e-2_wp/10.0_wp ! M∕Uo
       ! ! if (L_ref==1.0_wp) tscale=tscale/0.087318767977937 ! rescale length

       ! -------------------------------
       !     Standard initialisation
       ! -------------------------------
       ! Compute eddy turnover time
       ! --------------------------
       ! tau_eddy = eps_v / qq**2         ! Paper Saffman 1963 - tau_eddy = epsi / u'**2
       tau_eddy = tau_e * L_ref / u_ref
       tscale = tau_eddy

       if (iproc==0) write(6,'(''     tscale:'',F19.12)') tscale
    else
       call mpistop('NOT IMPLEMENTED YET ...',0)
    endif
    
    if (iproc==0) then
       deltat=1.*deltay/(u_ref+c_ref)
       print *,'deltat',deltat

       print *,int(tscale*(98-42)/deltat)+1 
       print *,int(tscale*(171-42)/deltat)+1 
    endif
    !call mpistop('stop at the end of init HIT',0)  

    ! Define new u_ref for timestep calculation
    ! =========================================
    ! (based on the maximum fuctuations in the synthetic initial velocity field)
    u_ref=max_rms
    
    ! FFTW: Free planes and memory from malloc
    ! ========================================
    if (iproc==0) print *,'~> free FFTW'
    call x_fftw_exit

    ! Write information file for restart
    ! ==================================    
    if (iproc==0) then
       inquire(file='vn2dudx2.dat',exist=iexist)
       if (iexist) then
          open(69,file='vn2dudx2.dat',status='replace')
       else
          open(69,file='vn2dudx2.dat',status='new')
       endif
       write(69,*) qq
       write(69,*) dudx2m
       write(69,*) max_rms+c_ref
       !write(69,*) vort2_moy
       close(69)

       !write(*,*) 'Initial resolution:',dble(nx/2)*sqrt(muw_ref*diffscale/rho_ref)*vort2_moy**(-0.25_wp)
    endif

!!$    ! Write restart file with MPI-IO
!!$    ! ==============================
!!$    call read_write_volume(binfile,WRITE)
!!$    call mpistop('stop at the end of init HIT',0)

  end subroutine init_hit

  !===============================================================================
  subroutine init_velocity_hit
  !===============================================================================
    !> Initialization of velocity components in Fourier space
    !> for incompressible/compressible HIT case (MPI version)
    ! Two versions are available:
    ! ---------------------------
    ! 'C': version adapted from HIT3D code [spectral incompressible HIT code
    !      developed by Sergei Chumakov (Stanford U), Natalia Vladimirova (U of Chicago)
    !      and Misha Stepanov (U of Arizona) https://github.com/fbusabiaga/hit3d]
    !      * SG Chumakov, A priori study of subgrid-scale flux of a passive scalar
    !        in turbulence, Phys.Rev.E, 78 15563 (2008)
    !      * SG Chumakov, Scaling properties of subgrid-scale energy dissipation,
    !        Phys. Fluids, 19 058104 (2007)
    ! 'PG': version adapted from Pirozzoli and Grasso [courtesy of Francesco Grasso]
    !      * S Pirozzoli, F Grasso, Direct numerical simulations of isotropic compressible
    !        turbulence: influence of compressibility on dynamics and structures,
    !        Phys. Fluids, 16 4386-4407 (2004)
    !===============================================================================
    use ifport ! used for random number interface (intel fortran portability)
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------
    ! loop counters
    integer :: i,j,k,n  
    ! random number generator [version C]
    integer*8 :: seed1
    complex(wp), dimension(:,:,:,:), allocatable :: wrk
    ! variables for solenoidal field [version PG]
    real(wp) :: phi1,phi2,phi3,phi4,phi5,phi6 ! random phases
    real(wp) :: aa,bb,alphar,alphai,betar,betai
    complex(wp) :: alphat,betat
    real(wp) :: ak12 ! 2D wavenumber norm
    ! energy spectrum
    integer :: n_shell
    real(wp) :: fac,fac2
    real(wp), dimension(:), allocatable :: e_spec,e_spec1
    ! ----------------------------------------------------------------------------
    
    if (iproc==0) print *,'~> initialization of velocity components in Fourier space'
    
    if (trim(init_vel)=='C') then
       
       if (iproc==0) print *,'   [adapted from Chumakov et al.]'

       ! Initialize random number generator
       ! ==================================
       seed1=67748937.
       rseed=dble(seed1)
       fac=d_random(-rseed)
       do i=1,100
          fac=d_random(rseed)
       enddo
       ! advance the random number sequence for each proc
       do i=1,nx*ny*nz*6*iproc
          fac=d_random(rseed)
       enddo

       ! Compute velocity components in Fourier space
       ! ============================================

       ! Filling the working arrays wrk1...wrk6 with random numbers
       ! ----------------------------------------------------------
       allocate(wrk(nx,ny,nz,6))
       do n=1,6
          do k=1,nz
             do j=1,ny
                do i=1,nx+2 ! to compare with Chumakov routine TO BE CHANGED
                   if (i<=nx) then
                      wrk(i,j,k,n)=d_random(rseed)
                   else
                      fac=d_random(rseed)
                   endif
                enddo
             enddo
          enddo
       enddo

       ! Making three random arrays with Gaussian PDF out of the six arrays that we generated
       ! --------------------------------------------
       ! [reuse arrays for derivatives as temporary arrays for fft]
       us_x=sqrt(-2.0_wp*log(wrk(:,:,:,1)))*sin(two_pi*wrk(:,:,:,4))
       vs_x=sqrt(-2.0_wp*log(wrk(:,:,:,2)))*sin(two_pi*wrk(:,:,:,5))
       ws_x=sqrt(-2.0_wp*log(wrk(:,:,:,3)))*sin(two_pi*wrk(:,:,:,6))
       deallocate(wrk)

       ! Go to Fourier space
       ! -------------------
       call x_fftw3d_f(us_x)
       call x_fftw3d_f(vs_x)
       call x_fftw3d_f(ws_x)

       ! Making three arrays that have the incompressibility property
       ! ------------------------------------------------------------
       do k=1,nz
          do j=1,ny
             do i=1,nx
                us(i,j,k)=i_*aky(j)*ws_x(i,j,k)-i_*akz(k)*vs_x(i,j,k)
                vs(i,j,k)=i_*akz(k)*us_x(i,j,k)-i_*akx(i)*ws_x(i,j,k)
                ws(i,j,k)=i_*akx(i)*vs_x(i,j,k)-i_*aky(j)*us_x(i,j,k)
             enddo
          enddo
       enddo

    elseif (trim(init_vel)=='PG') then

       if (iproc==0) print *,'   [adapted from Pirozzoli and Grasso]'
       
       ! Initialize random number generator
       ! ==================================
       call srand(1)
       ! advance the random number sequence for each proc
       do i=1,nx*ny*nz*6*iproc
          phi1 = rand()
       enddo

       ! Construct velocity components in Fourier space
       ! ==============================================
       do k=1,nz
          do j=1,ny
             do i=1,nx
                ! wavenumber module squared
                ak2=akx(i)*akx(i)+aky(j)*aky(j)+akz(k)*akz(k)
                if ((ak2).lt.0.5_wp) cycle
                ! wavenumber module
                ak=sqrt(ak2)
                ! random generated phases
                phi1 = two_pi*rand()
                phi2 = two_pi*rand()
                phi3 = two_pi*rand()
                phi4 = two_pi*rand()
                phi5 = two_pi*rand()
                phi6 = two_pi*rand()
                ! solenoidal velocity field
                aa=cos(phi3)
                bb=sin(phi3)
                alphar=aa*cos(phi1)
                alphai=aa*sin(phi1)
                betar =bb*cos(phi2)
                betai =bb*sin(phi2)
                alphat=cmplx(alphar,alphai)
                betat =cmplx(betar,betai)
                ak12=sqrt(akx(i)*akx(i)+aky(j)*aky(j))
                if (ak12.le.1e-12_wp) then
                   us(i,j,k)= alphat
                   vs(i,j,k)= betat
                   ws(i,j,k)= (0.0_wp,0.0_wp)
                else
                   us(i,j,k)= (alphat*ak*aky(j)+betat*akx(i)*akz(k))/(ak*ak12)
                   vs(i,j,k)= (betat*aky(j)*akz(k)-alphat*ak*akx(i))/(ak*ak12)
                   ws(i,j,k)= -betat*ak12/ak
                endif
             enddo
          enddo
       enddo
    endif

    ! Making the spectrum to be what it should
    ! ========================================

    ! Normalization factor (because the FFT is unnormalized)
    ! --------------------
    fac=0.5_wp/real(ngx*ngy*ngz)**2
    ! [0.5 because kinetic energy is 0.5*ui^2]
 
    ! Allocate energy spectra
    ! -----------------------
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

    ! Compute desired energy spectrum
    ! -------------------------------    
    do k=1,kmax
       ak=dble(k)
       e_spec1(k)=spectr(ak)
    enddo
    ! normalize it so it has the unit total energy
    e_spec1=e_spec1/sum(e_spec1(1:kmax))
    
    do k=1,nz
       do j=1,ny
          do i=1,nx
             ak2=akx(i)**2+aky(j)**2+akz(k)**2
             n_shell=nint(sqrt(ak2))
             if ((n_shell.gt.0).and.(n_shell.le.kmax)) then
                if (e_spec(n_shell).gt.0.0) then
                   us(i,j,k)=us(i,j,k)*sqrt(e_spec1(n_shell)/e_spec(n_shell))
                   vs(i,j,k)=vs(i,j,k)*sqrt(e_spec1(n_shell)/e_spec(n_shell))
                   ws(i,j,k)=ws(i,j,k)*sqrt(e_spec1(n_shell)/e_spec(n_shell))
                else
                   us(i,j,k)=0.0_wp
                   vs(i,j,k)=0.0_wp
                   ws(i,j,k)=0.0_wp
                endif
             else
                us(i,j,k)=0.0_wp
                vs(i,j,k)=0.0_wp
                ws(i,j,k)=0.0_wp
             endif
          enddo
       enddo
    enddo

    if (iproc==0) print *,"~> solenoidal velocities generated in Fourier space"

  end subroutine init_velocity_hit
  
  !===============================================================================
  subroutine init_pressure_hit
  !===============================================================================
    !> Initialization of solenoidal pressure (MPI version)
    !> Solve Poisson problem with RHS from solenoidal velocities
    !
    ! subroutine adapted from HIT3D code [spectral incompressible HIT code
    ! developed by Sergei Chumakov (Stanford U), Natalia Vladimirova (U of Chicago)
    ! and Misha Stepanov (U of Arizona) https://github.com/fbusabiaga/hit3d]
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! loop counters
    integer :: i,j,k
    ! variables
    real(wp) :: akmax2
    ! ----------------------------------------------------------------------------
    
    if (iproc==0) print *,'~> solve Poisson problem for pressure ..'

    ! Compute RHS for Poisson problem (dimensional quantities)
    ! ===============================
    ! (Att. array ps is used for rhs)
    do k=1,nz
       do j=1,ny
          do i=1,nx
             ps(i,j,k)=-( us_x(i,j,k)**2+vs_y(i,j,k)**2+ws_z(i,j,k)**2 &
                        + 2.0_wp*us_y(i,j,k)*vs_x(i,j,k) &
                        + 2.0_wp*us_z(i,j,k)*ws_x(i,j,k) &
                        + 2.0_wp*ws_y(i,j,k)*vs_z(i,j,k) )*rho_ref
          enddo
       enddo
    enddo

    ! RHS: go to Fourier space
    ! ========================
    ! (Att. array ps is used for rhs)
    call x_fftw3d_f(ps) 

    ! Solve Poisson problem in Fourier space
    ! ======================================
    ! (Att. array ps is used for rhs)
    do k=1,nz
       do j=1,ny
          do i=1,nx
             ! inverse Laplace operator
             ak2=akx(i)**2+aky(j)**2+akz(k)**2
             if (ak2.eq.0.0_wp) ak2=9.e20_wp

             ! calculating pressure
             ps(i,j,k)=-ps(i,j,k)/ak2
          enddo
       enddo
    enddo

    ! Mode truncation to ensure isotropy (for all modes that are higher than kmax)
    ! ==================================
    akmax2=dble(kmax)**2
    do k=1,nz
       do j=1,ny
          do i=1,nx
             ak2=akx(i)**2+aky(j)**2+akz(k)**2
             if (ak2.gt.akmax2) ps(i,j,k)=0.0_wp
          enddo
       enddo
    enddo

    ! Pressure: go back from Fourier space
    ! ====================================
    call x_fftw3d_b(ps) 

  end subroutine init_pressure_hit

  !===============================================================================
  function spectr(ak)
  !===============================================================================
    !> Init HIT: Models of HIT energy spectrum
    ! spectrum models:
    ! ----------------
    ! 'Ko': plain Kolmogorov spectrum
    ! 'PP': Passot-Pouquet spectrum [need peak wavenumber ak0]
    ! 'Ex': exponential spectrum [need peak wavenumber ak0]
    ! 'VK': von Karman spectrum
    ! 'VK+P': von Karman-Pao spectrum
    ! 'VK+S': von Karman-Saffman spectrum
    ! 'VK+P+B': von Karman-Pao spectrum + bottleneck correction
    ! 'VK+S+B': von Karman-Saffman spectrum + bottleneck correction
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: ak ! INPUT: wavenumber
    real(wp) :: spectr         ! OUTPUT: spectrum
    ! ----------------------------------------------------------------------------
    ! local variables
    real(wp) :: fac,ratio,kkol,ke
    ! ----------------------------------------------------------------------------

    if (spectr_model=='Ko') then
       
       ! Plain Kolmogorov spectrum
       ! =========================
       spectr=ak**(-5.0_wp/3.0_wp)
       
    elseif (spectr_model=='PP') then
       
       ! Passot-Pouquet spectrum
       ! =======================
       ratio=ak/ak0
       spectr=ak**4*exp(-2.0_wp*ratio**2)
       
    elseif (spectr_model=='Ex') then
       
       ! Exponential spectrum
       ! ====================
       ratio=ak/ak0
       spectr = ratio**3/ak0*exp(-3.0_wp*ratio)
       
    elseif (spectr_model=='VK') then
       
       ! Von Karman spectrum
       ! ===================
       ! version 1
       !fac =two_pi*ratio
       !spectr = fac**4/(1.0_wp+fac**2)**3
       ! version 2
       ke=0.747_wp/Lf_int
       ! Nota: Uref=sqrt(3./2.)*urms and uf=urms/Uref -> uf=sqrt(2./3.)
       ! -> 1.453*uf**2=1.453*2./3.=0.9687
       ratio=ak/ke
       spectr=0.9687_wp/ke*ratio**4/exp(17.0_wp/6.0_wp*log(1.0_wp+ratio**2))

       ! Define viscous dissipation range
       ! ================================
       ! Nota: Kolmogorov's constant: cK=1.613       
       ! Kolmogorov's wavenumber
       kkol=1.0_wp/eta_kol
       
       ratio=ak/kkol       
       if (dissip_range=='P') then
          ! Pao model
          ! ---------
          spectr=spectr*exp(-1.5_wp*1.613_wp*(ratio**(4.0_wp/3.0_wp)))
       elseif (dissip_range=='S') then
          ! Saffman model
          ! -------------
          spectr=spectr*exp(-1.5_wp*1.613_wp*(ratio**2))          
       endif

       ! Add bottleneck correction [Kang, Chester & Meneveau, JFM]
       ! =========================
       if (bottleneck=='B') then
          fac=1.0_wp+0.522_wp*(0.5_wp+atan(10.0_wp*log10(ak*eta_kol)+12.58_wp)/pi)
          spectr=spectr*fac
       endif
    
    endif

  end function spectr

  !===============================================================================
  subroutine x_fftw_init
  !===============================================================================
    !> Init HIT: subroutine that initializes the auxiliary arrays for FFT
  !===============================================================================
    use mod_mpi_part
    implicit none
    ! ----------------------------------------------------------------------------
    integer(C_INTPTR_T) :: l_nx,l_ny,l_nz ! local dimensions
    integer(C_INTPTR_T) :: g_nx,g_ny,g_nz ! global dimensions
    integer(C_INTPTR_T) :: alloc_size
    ! ----------------------------------------------------------------------------

    if (iproc==0) print *,'~> initialize FFT arrays ..'

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
    ! x-derivatives
    cus_x=fftw_alloc_complex(alloc_size)
    call c_f_pointer(cus_x,us_x,[l_nx,l_ny,l_nz]) ! pointer on us_x
    cvs_x=fftw_alloc_complex(alloc_size)
    call c_f_pointer(cvs_x,vs_x,[l_nx,l_ny,l_nz]) ! pointer on vs_x
    cws_x=fftw_alloc_complex(alloc_size)
    call c_f_pointer(cws_x,ws_x,[l_nx,l_ny,l_nz]) ! pointer on ws_x
    ! y-derivatives
    cus_y=fftw_alloc_complex(alloc_size)
    call c_f_pointer(cus_y,us_y,[l_nx,l_ny,l_nz]) ! pointer on us_y
    cvs_y=fftw_alloc_complex(alloc_size)
    call c_f_pointer(cvs_y,vs_y,[l_nx,l_ny,l_nz]) ! pointer on vs_y
    cws_y=fftw_alloc_complex(alloc_size)
    call c_f_pointer(cws_y,ws_y,[l_nx,l_ny,l_nz]) ! pointer on ws_y
    ! z-derivatives
    cus_z=fftw_alloc_complex(alloc_size)
    call c_f_pointer(cus_z,us_z,[l_nx,l_ny,l_nz]) ! pointer on us_z
    cvs_z=fftw_alloc_complex(alloc_size)
    call c_f_pointer(cvs_z,vs_z,[l_nx,l_ny,l_nz]) ! pointer on vs_z
    cws_z=fftw_alloc_complex(alloc_size)
    call c_f_pointer(cws_z,ws_z,[l_nx,l_ny,l_nz]) ! pointer on ws_z
    ! pressure & RHS for Poisson problem
    if (is_poisson) then
       cps=fftw_alloc_complex(alloc_size)
       call c_f_pointer(cps,ps,[l_nx,l_ny,l_nz]) ! pointer on ps
    endif
    
    ! Create MPI planes for in-place backward DFT
    ! ===========================================
    plan_b_x=fftw_mpi_plan_dft_1d(g_nx,us(:,1,1),us(:,1,1),COMMYZ,FFTW_BACKWARD,FFTW_ESTIMATE)
    plan_b_y=fftw_mpi_plan_dft_1d(g_ny,us(1,:,1),us(1,:,1),COMMXZ,FFTW_BACKWARD,FFTW_ESTIMATE)
    plan_b_z=fftw_mpi_plan_dft_1d(g_nz,us(1,1,:),us(1,1,:),COMMXY,FFTW_BACKWARD,FFTW_ESTIMATE)
    
    ! Create MPI planes for in-place forward DFT
    ! ==========================================
    plan_f_x=fftw_mpi_plan_dft_1d(g_nx,us(:,1,1),us(:,1,1),COMMYZ,FFTW_FORWARD,FFTW_ESTIMATE)
    plan_f_y=fftw_mpi_plan_dft_1d(g_ny,us(1,:,1),us(1,:,1),COMMXZ,FFTW_FORWARD,FFTW_ESTIMATE)
    plan_f_z=fftw_mpi_plan_dft_1d(g_nz,us(1,1,:),us(1,1,:),COMMXY,FFTW_FORWARD,FFTW_ESTIMATE)

  end subroutine x_fftw_init
  
  !===============================================================================
  subroutine x_fftw3d_f(var)
  !===============================================================================
    !> Init HIT: 3x1D (forward) Fourier transforms
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
  subroutine x_fftw3d_b(var)
  !===============================================================================
    !> Init HIT: 3x1D inverse (backward) Fourier transforms + FFT normalization
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k
    complex(C_DOUBLE_COMPLEX), dimension(:,:,:), pointer :: var
    ! ----------------------------------------------------------------------------
    
    ! Perform successive backward FFT
    ! ===============================
    do k=1,nz
       do j=1,ny
          call fftw_mpi_execute_dft(plan_b_x,var(:,j,k),var(:,j,k))
       enddo
    enddo
    do k=1,nz
       do i=1,nx
          call fftw_mpi_execute_dft(plan_b_y,var(i,:,k),var(i,:,k))
       enddo
    enddo
    do i=1,nx
       do j=1,ny
          call fftw_mpi_execute_dft(plan_b_z,var(i,j,:),var(i,j,:))
       enddo
    enddo

    ! FFT normalization
    ! =================
    var=var*norm

  end subroutine x_fftw3d_b
  
  !===============================================================================
  subroutine x_fftw_exit
  !===============================================================================
    !> Init HIT: subroutine to free FFTW planes and memory
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------
    
    ! Destroy FFTW planes
    ! ===================
    call fftw_destroy_plan(plan_b_x)
    call fftw_destroy_plan(plan_b_y)
    call fftw_destroy_plan(plan_b_z)
    call fftw_destroy_plan(plan_f_x)
    call fftw_destroy_plan(plan_f_y)
    call fftw_destroy_plan(plan_f_z)
    
    ! Free memory from malloc
    ! =======================
    call fftw_free(cus)
    call fftw_free(cvs)
    call fftw_free(cws)
    call fftw_free(cus_x)
    call fftw_free(cvs_x)
    call fftw_free(cws_x)
    call fftw_free(cus_y)
    call fftw_free(cvs_y)
    call fftw_free(cws_y)
    call fftw_free(cus_z)
    call fftw_free(cvs_z)
    call fftw_free(cws_z)
    if (is_poisson) then
       call fftw_free(cps)
    endif

  end subroutine x_fftw_exit

  !===============================================================================
  function d_random(idum3)
  !===============================================================================
    !> Double precision random number generator
    !
    ! subroutine taken from HIT3D code [spectral incompressible HIT code
    ! developed by Sergei Chumakov (Stanford U), Natalia Vladimirova (U of Chicago)
    ! and Misha Stepanov (U of Arizona) https://github.com/fbusabiaga/hit3d]
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
    real(wp) :: d_random,am,eps,rnmx,idum3
    parameter (im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1, &
         ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,  &
         ir2=3791,ntab=32,ndiv=1+imm1/ntab,eps=1.2e-7,rnmx=1.0-eps)
    integer :: idum2,j,k,iv(ntab),iy
    ! ----------------------------------------------------------------------------

    save iv,iy,idum2
    data idum2/123456789/, iv/ntab*0/, iy/0/

    ! initialize with idum<0; Then keep idum>0 unchanged for the same sequence
    idum=nint(idum3/987)
    if (idum.le.0) then
       idum=max(-idum,1)
       idum2=idum
       do j=ntab+8,1,-1
          k=idum/iq1
          idum=ia1*(idum-k*iq1)-k*ir1
          if(idum.lt.0) idum=idum+im1
          if(j.le.ntab) iv(j)=idum
       enddo
    endif

    ! start here when not initializing
    k=idum/iq1
    idum=ia1*(idum-k*iq1)-k*ir1
    if (idum.lt.0) idum=idum+im1
    k=idum2/iq2
    idum2=ia2*(idum2-k*iq2)-k*ir2
    if (idum2.lt.0) idum2=idum2+im2
    j=1+iy/ndiv
    iy=iv(j)-idum2
    iv(j)=idum
    if(iy.lt.1) iy=iy+imm1
    d_random=min(am*iy,rnmx)

    return

  end function d_random

  !===============================================================================
  subroutine stat_velocity
  !===============================================================================
    !> Compute velocity statististics [non-dimensional version]
    !
    ! subroutine inspired by HIT3D code [spectral incompressible HIT code
    ! developed by Sergei Chumakov (Stanford U), Natalia Vladimirova (U of Chicago)
    ! and Misha Stepanov (U of Arizona) https://github.com/fbusabiaga/hit3d]
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k,n_shell
    real(wp) :: fac,fac2
    real(wp) :: sctmp
    real(wp), dimension(:), allocatable :: e_spec,e_spec1

    real(wp) :: energy,eps_v,eta,etakmax,re_lambda,uvar,x_length
    real(wp) :: lambda!,enstrophy
    ! real(wp) :: lambda,tau_e!,enstrophy
    real(wp) :: nu,Re
    ! ----------------------------------------------------------------------------

    ! getting the enstrophy
    !call get_gradient_statistics

    ! getting the energy spectrum e_spec to the main process
    !call get_e_spec
    
    ! Normalization factor (because the FFT is unnormalized)
    ! --------------------
    fac=0.5_wp/real(ngx*ngy*ngz)**2
    ! [0.5 because kinetic energy is 0.5*ui^2]
   
    ! Allocate energy spectra
    ! -------------------------
    allocate(e_spec(kmax),e_spec1(kmax))

    ! Assembling the energy spectrum in each shell
    ! --------------------------------------------
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

    ! Write spectrum
    ! --------------
    if (iproc==0) then
       open(unit=49,file='spectr.dat',status='unknown',form='formatted')
       rewind(49)
       do k=1,kmax
          write(49,*) k,e_spec(k)
       enddo
       close(49)
    endif

    ! Total energy
    ! ------------
    energy=sum(e_spec(1:kmax))
    if (iproc==0) print 100,energy,energy*uref**2

    ! Dissipation spectrum and total dissipation
    ! ------------------------------------------
    Re=rho_ref*Uref*L_ref/muw_ref
    nu=1.0_wp/Re
    do k=1,kmax
       e_spec1(k)=2.0_wp*nu*e_spec(k)*dble(k**2)
    end do
    eps_v=sum(e_spec1(1:kmax))
    if (iproc==0) print 101,eps_v,eps_v/(1e2*L_ref)*(1e2*Uref)**3

    ! Kolmogorov scale
    ! ----------------
    eta=(nu**3/eps_v)**0.25
    etakmax=eta*dble(kmax)
    if (iproc==0) print 102,eta,eta*L_ref*1e2

    ! Variance
    ! --------
    uvar=2.0_wp/3.0_wp*energy
    if (iproc==0) print 103,uvar,sqrt(uvar)*uref*1e2
    
    ! Integral length scale
    ! ---------------------
    sctmp=0.0_wp
    do k=1,kmax
       sctmp=sctmp+e_spec(k)/dble(k)
    enddo
    x_length=pi/2.0_wp*sctmp/uvar
    if (iproc==0) print 104,x_length*L_ref*1e2

    ! Taylor microscale
    ! -----------------
    lambda=sqrt(15.0_wp*uvar*nu/eps_v)
    if (iproc==0) print 105,lambda*L_ref*1e2

    ! Taylor-Reynolds number
    ! ----------------------
    !re_lambda =uvar*sqrt(15.0_wp/eps_v*RE)
    re_lambda=sqrt(uvar)*lambda/nu
    if (iproc==0) print 106,re_lambda

    ! Eddy turnover time
    ! ------------------
    tau_e=x_length/sqrt(uvar)
    if (iproc==0) print 107,tau_e

100 format(4x,'>> energy int(E(k),dk)',F8.4,' (redim ',F8.4,' cm^2/s^2)')
101 format(4x,'>> dissipation 2*nu*int(k^2*E(k),dk)',F8.4,' (redim ',F10.4,' cm^2/s^3)')  
103 format(4x,'>> u variance',F8.4,' -> urms',F12.4,' cm/s')
102 format(4x,'>> Kolmogorov''s microscale',F8.4,' (redim ',F8.4,' cm)')
104 format(4x,'>> integral length scale (redim ',F8.4,' cm)')
105 format(4x,'>> Taylor''s microscale (redim ',F8.4,' cm)')
106 format(4x,'>> Re_Taylor',F12.2)
107 format(4x,'>> Eddy turnover time',F8.4)

   end subroutine stat_velocity

  !===============================================================================
  subroutine stat_velocity_dim
  !===============================================================================
    !> Compute velocity statististics [dimensional version]
    !
    ! subroutine inspired by HIT3D code [spectral incompressible HIT code
    ! developed by Sergei Chumakov (Stanford U), Natalia Vladimirova (U of Chicago)
    ! and Misha Stepanov (U of Arizona) https://github.com/fbusabiaga/hit3d]
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k,n_shell
    real(wp) :: fac,fac2
    real(wp) :: sctmp
    real(wp), dimension(:), allocatable :: e_spec,e_spec1

    real(wp) :: energy,eps_v,eta,etakmax,re_lambda,uvar,x_length
    real(wp) :: lambda,tau_e!,enstrophy
    real(wp) :: nu
    ! ----------------------------------------------------------------------------
    
    call x_fftw3d_f(us) 
    call x_fftw3d_f(vs) 
    call x_fftw3d_f(ws) 
!!$    call stat_velocity_dim

    ! getting the enstrophy
    !call get_gradient_statistics

    ! getting the energy spectrum e_spec to the main process
    !call get_e_spec
    
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

    ! Write spectrum
    ! --------------
    if (iproc==0) then
       open(unit=49,file='spectr_dim.dat',status='unknown',form='formatted')
       rewind(49)
       do k=1,kmax
          write(49,*) k/L_ref,e_spec(k)
       enddo
       close(49)
    endif
 
    ! Total energy
    ! ------------
    energy=sum(e_spec(1:kmax))/L_ref
    if (iproc==0) print *,'   >> energy int(E(k),dk)',energy

    ! Dissipation spectrum and total dissipation
    ! ------------------------------------------
    nu=mu_ref/rho_ref
    do k=1,kmax
       ak=dble(k)/L_ref
       e_spec1(k)=2.0_wp*nu*e_spec(k)*ak**2
    end do
    eps_v=sum(e_spec1(1:kmax))/L_ref
    if (iproc==0) print *,'   >> dissipation 2*nu*int(k^2*E(k),dk)',eps_v,eps_v/(1e2)*(1e2)**3

    ! Kolmogorov scale
    ! ----------------
    eta=(nu**3/eps_v)**0.25
    etakmax=eta*dble(kmax)/L_ref
    if (iproc==0) print *,'   >> Kolmogorov scale (cm)',eta*1e2

    ! Variance
    ! --------
    uvar=2.0_wp/3.0_wp*energy
    if (iproc==0) print *,'   >> urms (cm/s)',sqrt(uvar)*1e2
    
    ! Integral length scale
    ! ---------------------
    sctmp=0.0_wp
    do k=1,kmax
       ak=dble(k)/L_ref
       sctmp=sctmp+e_spec(k)/ak
    enddo
    x_length=pi/2.0_wp*sctmp/uvar/L_ref
    if (iproc==0) print *,'   >> integral scale (cm)',x_length*1e2

    ! Taylor microscale
    ! -----------------
    lambda=sqrt(15.0_wp*uvar*nu/eps_v)
    if (iproc==0) print *,'   >> Taylor scale (cm)',lambda*1e2

    ! Taylor-Reynolds number
    ! ----------------------
    !re_lambda =uvar*sqrt(15.0_wp/eps_v*RE)
    re_lambda=sqrt(uvar)*lambda/nu
    if (iproc==0) print *,'   >> Re_Taylor',re_lambda

    ! Eddy turnover time
    ! ------------------
    tau_e=x_length/sqrt(uvar)
    if (iproc==0) print *,'   >> Eddy turnover time',tau_e

  end subroutine stat_velocity_dim

!!$  !===============================================================================
!!$  subroutine get_e_spec
!!$  !===============================================================================
!!$    !> Compute the energy spectrum
!!$  !===============================================================================
!!$    implicit none
!!$    ! ----------------------------------------------------------------------------
!!$    integer :: i,j,k,n_shell
!!$    real(wp) :: fac,fac2
!!$    ! ----------------------------------------------------------------------------
!!$    
!!$    ! normalization factor (because the FFT is unnormalized)
!!$    ! --------------------
!!$    fac=0.5_wp/real(nx*ny*nz)**2
!!$    ! [0.5 because kinetic energy is 0.5*ui^2]
!!$
!!$    e_spec=0.0_wp
!!$
!!$    ! assembling the total energy in each shell
!!$    ! -----------------------------------------
!!$    do k=1,nz
!!$       do j=1,ny
!!$          do i=1,nx
!!$             ak2=akx(i)**2+aky(j)**2+akz(k)**2
!!$             n_shell=nint(sqrt(ak2))
!!$             if ((n_shell.gt.0).and.(n_shell.le.kmax)) then
!!$                fac2=fac*(abs(us(i,j,k))**2+abs(vs(i,j,k))**2+abs(ws(i,j,k))**2)
!!$                e_spec(n_shell)=e_spec(n_shell)+fac2
!!$             endif
!!$          enddo
!!$       enddo
!!$    enddo
!!$
!!$  end subroutine get_e_spec

end module mod_init_hit
