!=================================================================================
module mod_eigenmode
!=================================================================================
  !> Read and enter eigenfunctions of unstable modes in the inlet plane
  ! Preliminary version XG 28/02/2020
  ! restricted to perfect gas and 5 points BC
  ! restricted to 2D (periodic) configuration  
!=================================================================================
  use mod_constant
  use mod_mpi
  use mod_flow
  use mod_eos
  implicit none
  !-------------------------------------------------------------------------------
  ! indicator to enter eigenmodes at inlet
  logical :: is_eigenmode
  !-------------------------------------------------------------------------------
  ! from stability code
  character(len=50) :: ncas ! Name of the eigenfunctions
  integer :: ndata          ! Number of points discretizing the eigenfunctions
  logical :: is_oblique     ! Oblique modes (+/- beta)
  real(wp) :: limy          ! upper limit of the definition of eigenfunctions
  logical :: isro           ! if eigenfunction for rho is available
  real(wp), dimension(:), allocatable :: Qr,Qi,yd
  !-------------------------------------------------------------------------------
  ! eigenmode type
  type eigenmode
     real(wp) :: om,kr,ki,beta,ampl
     real(wp), dimension(:), pointer :: ur,ui,dur,dui
     real(wp), dimension(:), pointer :: vr,vi,dvr,dvi
     real(wp), dimension(:), pointer :: wr,wi,dwr,dwi
     real(wp), dimension(:), pointer :: pr,pi,dpr,dpi
     real(wp), dimension(:), pointer :: rr,ri,drr,dri
     real(wp), dimension(:), pointer :: tr,ti,dtr,dti
  end type eigenmode
  !-------------------------------------------------------------------------------
  ! eigenmodes
  integer :: n_eig
  type(eigenmode), dimension(:), pointer :: eig
  !-------------------------------------------------------------------------------
  
contains

  !===============================================================================
  subroutine init_eigenmodes
  !===============================================================================
    !> Initialization of eigenfunctions
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: k!,j
    ! ----------------------------------------------------------------------------

    if (iproc==0) then
       print *,repeat('=',70)
       print *,'Init eigenmodes'
       print *,repeat('=',70)
    endif

    ! Read unstable mode parameters to define eigenmodes
    ! ==================================================
    call read_param_stab

    ! Dimensionalize eigenmodes wavenumbers and frequency
    ! ===================================================
    ! Nota : L_ref=h_ref for step case and L_ref=deltas_forc for STBL
    
    ! ---> dans setupref.f90
    !omeg_sb = omeg_sb/deltas_forc*u_ref
    !if (is_2d) then
    !   beta_sb = 0.0_wp
    !else
    !   beta_sb = beta_sb/deltas_forc
    !endif
    
    do k=1,n_eig
       eig(k)%om=eig(k)%om/L_ref*u_ref
       eig(k)%kr=eig(k)%kr/L_ref
       eig(k)%ki=eig(k)%ki/L_ref
       eig(k)%beta=eig(k)%beta/L_ref
    enddo

    ! Allocate memory for eigenfunctions
    ! ==================================
    do k=1,n_eig
       allocate(eig(k)%ur(ny),eig(k)%ui(ny),eig(k)%dur(ny),eig(k)%dui(ny))
       allocate(eig(k)%vr(ny),eig(k)%vi(ny),eig(k)%dvr(ny),eig(k)%dvi(ny))
       allocate(eig(k)%wr(ny),eig(k)%wi(ny),eig(k)%dwr(ny),eig(k)%dwi(ny))
       allocate(eig(k)%pr(ny),eig(k)%pi(ny),eig(k)%dpr(ny),eig(k)%dpi(ny))
       allocate(eig(k)%rr(ny),eig(k)%ri(ny),eig(k)%drr(ny),eig(k)%dri(ny))
       allocate(eig(k)%tr(ny),eig(k)%ti(ny),eig(k)%dtr(ny),eig(k)%dti(ny))
    enddo
        
    allocate(Qr(ndata),Qi(ndata),yd(ndata))
    ! isro=F -> classical LST solver: T profile is given and rho is deduced
    ! isro=T -> LST solver with 5 variables: rho profile directly available
    isro=.true.

    ! Eigenfunctions
    ! ==============
    do k=1,n_eig
       !if (k==1) ncas='_M0p9_w0p1_3d'
       !if (k==2) ncas='_M0p9_w0p1_3d'
       !if (k==3) ncas='_M0p9_w0p2_2d'
       call def_eigenfunction(eig(k))
    enddo

!!$    ! check
!!$    if (iproc==0) then
!!$       open(40,file='eig_interp2.bin',form='unformatted',status='unknown')
!!$       rewind(40)
!!$       write(40) ny
!!$       write(40) (y(j),j=1,ny)
!!$       write(40) (eig(1)%ur(j),j=1,ny)
!!$       write(40) (eig(1)%ui(j),j=1,ny)
!!$       write(40) (eig(1)%vr(j),j=1,ny)
!!$       write(40) (eig(1)%vi(j),j=1,ny)
!!$       write(40) (eig(1)%wr(j),j=1,ny)
!!$       write(40) (eig(1)%wi(j),j=1,ny)
!!$       write(40) (eig(1)%pr(j),j=1,ny)
!!$       write(40) (eig(1)%pi(j),j=1,ny)
!!$       write(40) (eig(1)%tr(j),j=1,ny)
!!$       write(40) (eig(1)%ti(j),j=1,ny)
!!$       write(40) (eig(1)%rr(j),j=1,ny)
!!$       write(40) (eig(1)%ri(j),j=1,ny)
!!$       close(40)
!!$    endif
!!$    call mpistop('check eigenmodes',0)
    
    ! free memory
    deallocate(Qr,Qi,yd)
    if (iproc==0) print *,repeat('=',70)

  end subroutine init_eigenmodes
  
  !===============================================================================
  subroutine read_param_stab
  !===============================================================================
    !> Read parameters in param_stab.ini
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    logical :: iexist
    integer :: k,km,neig_read
    real(wp) :: psi
    ! ----------------------------------------------------------------------------

    ! Read param_stab.ini
    ! ===================
    inquire(file=trim(dirDATA)//'param_stab.ini', exist=iexist)
    if (.not.iexist) then
       call mpistop('Paramfile param_stab.ini does not exist!', 0)
    endif
    
    !=============================================================================
    open(30,file=trim('param_stab.ini'))
    !=============================================================================
    read(30,*)! ==================================================================
    read(30,*)! UNSTABLE MODE PARAMETERS (determined with stability solver)
    read(30,*)! ==================================================================
    read(30,*)! 
    read(30,*)! Name of the eigenfunctions
    read(30,*) ncas
    ncas='_'//trim(ncas)
    read(30,*)! Number of unstable modes
    read(30,*) neig_read
    read(30,*)! Oblique modes (+/- beta)
    read(30,*) is_oblique
    if (is_oblique) then
       n_eig=neig_read+1
    else
       n_eig=neig_read
    endif
    read(30,*)! Number of points discretizing the eigenfunctions
    read(30,*) ndata
    read(30,*)! upper limit of the definition of eigenfunctions
    read(30,*) limy

    ! Print number and name of eigenmodes
    ! ===================================
    if (iproc==0) then
       write(6,'(''    Number of eigenmodes:'',i3)') n_eig
       print *,'~> eigenmode name: fpropre_*'//trim(ncas)//'.dat'
    endif

    ! Allocate eigmode type
    ! =====================
    allocate(eig(n_eig))
 
    ! Read eigenmode characteristics
    ! ==============================
    km=1
    do k=1,neig_read
       
       read(30,*)! ==================================================================
       read(30,*)! Parameters of mode
       read(30,*)! ==================================================================
       read(30,*)! angular frequency
       read(30,*) eig(km)%om
       read(30,*)! real part of streamwise wavenumber
       read(30,*) eig(km)%kr
       read(30,*)! imaginary part of streamwise wavenumber
       read(30,*) eig(km)%ki
       read(30,*)!~> define 'wave angle' OR 'spanwise wavenumber'
                 !   (the first non-zero is taken into account)
       read(30,*)! wave angle: psi [in degrees]
       read(30,*) psi
       if (psi.ne.0.0_wp) then
          ! compute spanwise wavenumber from wave angle
          ! -------------------------------------------
          eig(km)%beta=eig(km)%kr*tan(psi*pi/180.0_wp)
          read(30,*)! spanwise wavenumber is ignored
          read(30,*) 
       else
          read(30,*)! spanwise wavenumber: beta
          read(30,*) eig(km)%beta
       endif
       read(30,*)! amplitude of mode excitation 
       read(30,*) eig(km)%ampl
       ! duplicate mode with -beta if is_oblique
       ! ---------------------------------------
       if (eig(km)%beta.ne.0.0_wp) then
          if (is_oblique) then
             eig(km+1)%om  = eig(km)%om
             eig(km+1)%kr  = eig(km)%kr
             eig(km+1)%ki  = eig(km)%ki
             eig(km+1)%beta=-eig(km)%beta
             eig(km+1)%ampl= eig(km)%ampl
             km=km+2
          else
             km=km+1
          endif
       else
          km=km+1
       endif
       
    enddo
    close(30)
    !=============================================================================

    ! Print nondimensional wavenumbers and frequency
    ! ==============================================
    if (iproc==0) then
       do k=1,n_eig
          write(6,*) '   ------------------'
          write(6,20) k
          write(6,*) '   ------------------'
          write(6,21) eig(k)%om
          write(6,22) eig(k)%kr,eig(k)%ki
          write(6,23) eig(k)%beta
          write(6,24) eig(k)%ampl
          write(6,25) 2.*pi/eig(k)%kr
          write(6,26) 2.*pi/eig(k)%kr*L_ref/(x(2)-x(1))
       enddo
    endif
20  format(4x,'Mode',i3)
21  format(4x,'angular frequency     (wL/U):',f8.4)
22  format(4x,'streamwise wavenumber  (k L):',f8.5,'+i*',f8.5)
23  format(4x,'spanwise wavenumber (beta L):',f8.2)
24  format(4x,'amplitude of mode      (eps):',f8.5)
25  format(4x,'wavelength        (lambda/L):',f8.2)
26  format(4x,'wavelength       (lambda/dx):',f8.2)

  endsubroutine read_param_stab

  !===============================================================================
  subroutine def_eigenfunction(eg)
  !===============================================================================
    !> Define eigenfunctions (read, interpolate, derive)
    !===================================================
    ! we need eigenfunction profiles for primitive variables (ui,p,rho)
    ! isro=F -> classical LST solver: T profile is given and rho is deduced
    ! isro=T -> LST solver with 5 variables: rho profile directly available
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    type(eigenmode) :: eg
    ! ----------------------------------------------------------------------------
    integer :: j
    real(wp) :: p_infty ! pressure adim
    ! ----------------------------------------------------------------------------

    if (iproc==0) print *,'~> read and interpolate eigenfunctions'

    ! Velocity component u
    ! ====================
    allocate(eg%ur(ny),eg%ui(ny),eg%dur(ny),eg%dui(ny))
    ! read, interpolate, derive
    call read_eigenfunction('u',eg%ur,eg%ui,eg%dur,eg%dui)
    ! dimensionalize
    eg%ur =eg%ampl*U_ref*eg%ur
    eg%ui =eg%ampl*U_ref*eg%ui
    eg%dur=eg%ampl*U_ref*eg%dur
    eg%dui=eg%ampl*U_ref*eg%dui
    
    ! Velocity component v
    ! ====================
    allocate(eg%vr(ny),eg%vi(ny),eg%dvr(ny),eg%dvi(ny))
    ! read, interpolate, derive
    call read_eigenfunction('v',eg%vr,eg%vi,eg%dvr,eg%dvi)
    ! dimensionalize
    eg%vr =eg%ampl*U_ref*eg%vr
    eg%vi =eg%ampl*U_ref*eg%vi
    eg%dvr=eg%ampl*U_ref*eg%dvr
    eg%dvi=eg%ampl*U_ref*eg%dvi
    
    ! Velocity component w
    ! ====================
    allocate(eg%wr(ny),eg%wi(ny),eg%dwr(ny),eg%dwi(ny))
    if (is_2d) then
       eg%wr =0.0_wp
       eg%wi =0.0_wp
       eg%dwr=0.0_wp
       eg%dwi=0.0_wp
    else
       ! read, interpolate, derive
       call read_eigenfunction('w',eg%wr,eg%wi,eg%dwr,eg%dwi)
       ! dimensionalize
       eg%wr =eg%ampl*U_ref*eg%wr
       eg%wi =eg%ampl*U_ref*eg%wi
       eg%dwr=eg%ampl*U_ref*eg%dwr
       eg%dwi=eg%ampl*U_ref*eg%dwi
    endif
    
    ! Pressure p
    ! ==========
    allocate(eg%pr(ny),eg%pi(ny),eg%dpr(ny),eg%dpi(ny))
    ! read, interpolate, derive
    call read_eigenfunction('p',eg%pr,eg%pi,eg%dpr,eg%dpi)
    ! dimensionalize
    p_infty=rho_ref*U_ref**2
    eg%pr =eg%ampl*p_infty*eg%pr
    eg%pi =eg%ampl*p_infty*eg%pi
    eg%dpr=eg%ampl*p_infty*eg%dpr
    eg%dpi=eg%ampl*p_infty*eg%dpi

    ! Temperature T
    ! =============
    allocate(eg%Tr(ny),eg%Ti(ny),eg%dTr(ny),eg%dTi(ny))
    ! read, interpolate, derive
    call read_eigenfunction('T',eg%Tr,eg%Ti,eg%dTr,eg%dTi)
    ! dimensionalize
    eg%Tr =eg%ampl*T_ref*eg%Tr
    eg%Ti =eg%ampl*T_ref*eg%Ti
    eg%dTr=eg%ampl*T_ref*eg%dTr
    eg%dTi=eg%ampl*T_ref*eg%dTi

    ! Density rho
    ! ===========
    allocate(eg%rr(ny),eg%ri(ny),eg%drr(ny),eg%dri(ny))
    
    if (isro) then
       ! read, interpolate, derive
       call read_eigenfunction('rho',eg%rr,eg%ri,eg%drr,eg%dri)
       ! dimensionalize
       eg%rr =eg%ampl*rho_ref*eg%rr
       eg%ri =eg%ampl*rho_ref*eg%ri
       eg%drr=eg%ampl*rho_ref*eg%drr
       eg%dri=eg%ampl*rho_ref*eg%dri
    else
       ! reconstruct rho from T and p
       eg%rr=(p_ref+eg%pr)/(T_ref+eg%Tr)/rg-rho_ref
       eg%ri=(p_ref+eg%pi)/(T_ref+eg%Ti)/rg-rho_ref
       ! reconstruct derivatives
       eg%drr=eg%dpr-eg%dTr*(p_ref+eg%pr)/(T_ref+eg%Tr)
       eg%drr=eg%drr/(T_ref+eg%Tr)/rg
       eg%dri=eg%dpi-eg%dTi*(p_ref+eg%pi)/(T_ref+eg%Ti)
       eg%dri=eg%dri/(T_ref+eg%Ti)/rg
    endif
    
    if (iproc==0) print *,'~> read eigenfunctions OK'

  end subroutine def_eigenfunction

  !===============================================================================
  subroutine read_eigenfunction(varname,Qr_,Qi_,dQr_,dQi_)
  !===============================================================================
    !> Read, interpolate and derive an eigenfunction on global grid
  !===============================================================================
    use mod_utils
    use mod_interp1 ! <- module for 1D interpolation [Mathematics]
    implicit none
    ! ----------------------------------------------------------------------------
    character(len=*)  :: varname
    real(wp), dimension(ny), intent(out) :: Qr_,Qi_,dQr_,dQi_
    ! ----------------------------------------------------------------------------
    ! local var
    integer :: j,uid,ndy_,nfy_
    real(wp), dimension(ny+2) :: Qr_e,Qi_e,dQr_e,dQi_e
    ! ----------------------------------------------------------------------------

    ! Load eigenfunction from file
    ! ============================
    call get_free_unit(uid)
    if (iproc==0) print *,'   fpropre_'//trim(varname)//trim(ncas)//'.dat'
    open(unit=uid,file='fpropre_'//trim(varname)//trim(ncas)//'.dat',form='formatted',status='unknown')
    read(uid,*)
    read(uid,*)
    do j=1,ndata
       read(uid,*) yd(j),Qr(j),Qi(j)
    enddo
    close(uid)
    
    ! Interpolate data on extended grid
    ! ================
    call interp1(y(0:ny+1)/deltas_in,Qr_e,ny+2,yd,Qr,ndata,'linear')
    call interp1(y(0:ny+1)/deltas_in,Qi_e,ny+2,yd,Qi,ndata,'linear')

    ! modify data that are too far from the wall
    do j=1,ngy
       if (y(j)/deltas_in.ge.limy) then
          Qr_e(j)=0.0_wp
          Qi_e(j)=0.0_wp
       endif
    enddo

    ! Restrict to ny
    ! ==============
    Qr_=Qr_e(1:ny)
    Qi_=Qi_e(1:ny)
    
    ! Derivative of eigenfunction with respect to y (second order)
    ! =============================================
    ndy_=1
    nfy_=ny
    ! at jmin
    if (BC_face(2,1)%sort.le.0) then
       dQr_(1)=(Qr_e(2)-Qr_e(1))/(y(2)-y(1))
       dQi_(1)=(Qi_e(2)-Qi_e(1))/(y(2)-y(1))
       ndy_=2
    endif
    ! at jmax
    if (BC_face(2,2)%sort.le.0) then
       dQr_(ny)=0.0_wp
       dQi_(ny)=0.0_wp
       nfy_=ny-1
    endif
    ! interior
    do j=ndy_,nfy_
       dQr_(j)=(Qr_e(j+1)-Qr_e(j-1))/(y(j+1)-y(j-1))
       dQi_(j)=(Qi_e(j+1)-Qi_e(j-1))/(y(j+1)-y(j-1))
    enddo
    
  end subroutine read_eigenfunction

  !===============================================================================
  subroutine eig_disturb1(i,j,k,um,vm,wm,rm,tm)
  !===============================================================================
    !> Compute eigenmode disturbances (for bc_supersonic_inlet)
  !===============================================================================
    use mod_time ! for: time,deltat,ck
    implicit none
    ! ----------------------------------------------------------------------------
    integer, intent(in) :: i,j,k ! grid point coord.
    real(wp), intent(out) :: um,vm,wm,rm,tm ! disturbances
    ! ----------------------------------------------------------------------------
    ! local var
    integer :: nm ! mode number 
    real(wp) :: trk,phi,ampl ! time,phase,amplitude
    ! ----------------------------------------------------------------------------
    
    ! time at RK stage
    ! ================
    trk =time+ck(irk)*deltat
    
    ! initialize derivatives of disturbances
    ! ======================================
    um=0.0_wp
    vm=0.0_wp
    wm=0.0_wp
    rm=0.0_wp
    tm=0.0_wp
    
    ! compute for all eigenmodes at point (i,j,k)
    ! ===========================================
    do nm=1,n_eig

       ! argument and amplitude
       ! ----------------------
       phi =eig(nm)%kr*(x(i)-x(1))-eig(nm)%om*trk+eig(nm)%beta*z(k)
       ampl=exp(-eig(nm)%ki*(x(i)-x(1)))

       ! distub for velocity u
       ! ---------------------
       um = um + ampl*(eig(nm)%ur(j)*cos(phi)-eig(nm)%ui(j)*sin(phi))

       ! distub for velocity v
       ! ---------------------
       vm = vm + ampl*(eig(nm)%vr(j)*cos(phi)-eig(nm)%vi(j)*sin(phi))

       ! distub for velocity w
       ! ---------------------
       wm = wm + ampl*(eig(nm)%wr(j)*cos(phi)-eig(nm)%wi(j)*sin(phi))

       ! distub for density rho
       ! ----------------------
       rm = rm + ampl*(eig(nm)%rr(j)*cos(phi)-eig(nm)%ri(j)*sin(phi))

       ! distub for temperature T 
       ! ------------------------
       tm = tm + ampl*(eig(nm)%tr(j)*cos(phi)-eig(nm)%ti(j)*sin(phi))
    enddo

  end subroutine eig_disturb1
             
  !===============================================================================
  subroutine eig_disturb1_dt(i,j,k,dumt,dvmt,dwmt,drmt,dpmt)
  !===============================================================================
    !> Compute time derivative of eigenmode disturbances (for bc_supersonic_inlet)
  !===============================================================================
    use mod_time ! for: time,deltat,ck
    implicit none
    ! ----------------------------------------------------------------------------
    integer, intent(in) :: i,j,k ! grid point coord.
    real(wp), intent(out) :: dumt,dvmt,dwmt,drmt,dpmt ! time derivatives of disturbances
    ! ----------------------------------------------------------------------------
    ! local var
    integer :: nm ! mode number 
    real(wp) :: trk,phi,ampl ! time,phase,amplitude
    ! ----------------------------------------------------------------------------
    
    ! time at RK stage
    ! ================
    trk =time+ck(irk)*deltat

    ! initialize derivatives of disturbances
    ! ======================================
    dumt=0.0_wp
    dvmt=0.0_wp
    dwmt=0.0_wp
    drmt=0.0_wp
    dpmt=0.0_wp
    
    ! compute for all eigenmodes at point (i,j,k)
    ! ===========================================
    do nm=1,n_eig

       ! argument and amplitude
       ! ----------------------
       phi =eig(nm)%kr*(x(i)-x(1))-eig(nm)%om*trk+eig(nm)%beta*z(k)
       ampl=eig(nm)%om*exp(-eig(nm)%ki*(x(i)-x(1)))

       ! time derivative of distub for velocity u
       ! ----------------------------------------
       dumt = dumt + ampl*(eig(nm)%ur(j)*sin(phi)+eig(nm)%ui(j)*cos(phi))

       ! time derivative of distub for velocity v
       ! ----------------------------------------
       dvmt = dvmt + ampl*(eig(nm)%vr(j)*sin(phi)+eig(nm)%vi(j)*cos(phi))

       ! time derivative of distub for velocity w
       ! ----------------------------------------
       dwmt = dwmt + ampl*(eig(nm)%wr(j)*sin(phi)+eig(nm)%wi(j)*cos(phi))

       ! time derivative of distub for density rho
       ! -----------------------------------------
       drmt = drmt + ampl*(eig(nm)%rr(j)*sin(phi)+eig(nm)%ri(j)*cos(phi))

       ! time derivative of distub for pressure p
       ! ----------------------------------------
       dpmt = dpmt + ampl*(eig(nm)%pr(j)*sin(phi)+eig(nm)%pi(j)*cos(phi))
    enddo

  end subroutine eig_disturb1_dt
             
  !===============================================================================
  subroutine eig_disturb2_imin(i,j,k,pt_in,ut_in,vt_in,wt_in,rt_in,dp_in,du_in,dv_in,dw_in,dr_in)
  !===============================================================================
    !> Impose eigenmode disturbances at inlet (imin Tam & Dong's BC)
  !===============================================================================
    use mod_time ! for: time,deltat,ck
    implicit none
    ! ----------------------------------------------------------------------------
    integer, intent(in)   :: i,j,k
    ! eigenmode disturbances
    real(wp), intent(out) :: pt_in,ut_in,vt_in,wt_in,rt_in
    real(wp), intent(out) :: dp_in,du_in,dv_in,dw_in,dr_in
    ! ----------------------------------------------------------------------------
    ! local var
    integer :: nm
    real(wp) :: trk,phi,ampl
    real(wp) :: qm,dqmx,dqmy,dqmt
    real(wp), dimension(:,:), pointer :: ir,cosphi,sinphi
    ! ----------------------------------------------------------------------------
    
    ir=>BC_face(1,1)%ir
    cosphi=>BC_face(1,1)%cosphi
    sinphi=>BC_face(1,1)%sinphi

    ! time at RK stage
    ! ================
    trk =time+ck(irk)*deltat

    ! pressure
    qm  =0.0_wp
    dqmx=0.0_wp
    dqmy=0.0_wp
    dqmt=0.0_wp
    do nm=1,n_eig   
       phi =eig(nm)%kr*(x(i)-x(1))-eig(nm)%om*trk+eig(nm)%beta*z(k)
       ampl=exp(-eig(nm)%ki*(x(i)-x(1)))

       qm   = (            eig(nm)%pr(j)*cos(phi)-           eig(nm)%pi(j)*sin(phi))*ampl + qm
       dqmx = (-eig(nm)%kr*eig(nm)%pr(j)*sin(phi)-eig(nm)%kr*eig(nm)%pi(j)*cos(phi))*ampl &
            + (            eig(nm)%pr(j)*cos(phi)-           eig(nm)%pi(j)*sin(phi))*(-eig(nm)%ki*ampl) + dqmx
       dqmy = (           eig(nm)%dpr(j)*cos(phi)-          eig(nm)%dpi(j)*sin(phi))*ampl + dqmy
       dqmt = ( eig(nm)%om*eig(nm)%pr(j)*sin(phi)+eig(nm)%om*eig(nm)%pi(j)*cos(phi))*ampl + dqmt
    enddo
    dp_in=dqmx*cosphi(i,j)+dqmy*sinphi(i,j)+qm*ir(i,j)
    pt_in=dqmt
    !dp_in=0.0_wp
    !pt_in=0.0_wp
    
    ! velocity u
    qm  =0.0_wp
    dqmx=0.0_wp
    dqmy=0.0_wp
    dqmt=0.0_wp
    do nm=1,n_eig   
       phi =eig(nm)%kr*(x(i)-x(1))-eig(nm)%om*trk+eig(nm)%beta*z(k)
       ampl=exp(-eig(nm)%ki*(x(i)-x(1)))

       qm   = (            eig(nm)%ur(j)*cos(phi)-           eig(nm)%ui(j)*sin(phi))*ampl + qm
       dqmx = (-eig(nm)%kr*eig(nm)%ur(j)*sin(phi)-eig(nm)%kr*eig(nm)%ui(j)*cos(phi))*ampl &
            + (            eig(nm)%ur(j)*cos(phi)-           eig(nm)%ui(j)*sin(phi))*(-eig(nm)%ki*ampl) + dqmx
       dqmy = (           eig(nm)%dur(j)*cos(phi)-          eig(nm)%dui(j)*sin(phi))*ampl + dqmy
       dqmt = ( eig(nm)%om*eig(nm)%ur(j)*sin(phi)+eig(nm)%om*eig(nm)%ui(j)*cos(phi))*ampl + dqmt
    enddo
    du_in=dqmx*cosphi(i,j)+dqmy*sinphi(i,j)+qm*ir(i,j)
    ut_in=dqmt
    
    ! velocity v
    qm  =0.0_wp
    dqmx=0.0_wp
    dqmy=0.0_wp
    dqmt=0.0_wp
    do nm=1,n_eig   
       phi =eig(nm)%kr*(x(i)-x(1))-eig(nm)%om*trk+eig(nm)%beta*z(k)
       ampl=exp(-eig(nm)%ki*(x(i)-x(1)))

       qm   = (            eig(nm)%vr(j)*cos(phi)-           eig(nm)%vi(j)*sin(phi))*ampl + qm
       dqmx = (-eig(nm)%kr*eig(nm)%vr(j)*sin(phi)-eig(nm)%kr*eig(nm)%vi(j)*cos(phi))*ampl &
            + (            eig(nm)%vr(j)*cos(phi)-           eig(nm)%vi(j)*sin(phi))*(-eig(nm)%ki*ampl) + dqmx
       dqmy = (           eig(nm)%dvr(j)*cos(phi)-          eig(nm)%dvi(j)*sin(phi))*ampl + dqmy
       dqmt = ( eig(nm)%om*eig(nm)%vr(j)*sin(phi)+eig(nm)%om*eig(nm)%vi(j)*cos(phi))*ampl + dqmt
    enddo
    dv_in=dqmx*cosphi(i,j)+dqmy*sinphi(i,j)+qm*ir(i,j)
    vt_in=dqmt
    
    ! velocity w
    qm  =0.0_wp
    dqmx=0.0_wp
    dqmy=0.0_wp
    dqmt=0.0_wp
    do nm=1,n_eig   
       phi =eig(nm)%kr*(x(i)-x(1))-eig(nm)%om*trk+eig(nm)%beta*z(k)
       ampl=exp(-eig(nm)%ki*(x(i)-x(1)))

       qm   = (            eig(nm)%wr(j)*cos(phi)-           eig(nm)%wi(j)*sin(phi))*ampl + qm
       dqmx = (-eig(nm)%kr*eig(nm)%wr(j)*sin(phi)-eig(nm)%kr*eig(nm)%wi(j)*cos(phi))*ampl &
            + (            eig(nm)%wr(j)*cos(phi)-           eig(nm)%wi(j)*sin(phi))*(-eig(nm)%ki*ampl) + dqmx
       dqmy = (           eig(nm)%dwr(j)*cos(phi)-          eig(nm)%dwi(j)*sin(phi))*ampl + dqmy
       dqmt = ( eig(nm)%om*eig(nm)%wr(j)*sin(phi)+eig(nm)%om*eig(nm)%wi(j)*cos(phi))*ampl + dqmt
    enddo
    dw_in=dqmx*cosphi(i,j)+dqmy*sinphi(i,j)+qm*ir(i,j)
    wt_in=dqmt
    !dw_in=0.
    !wt_in=0.
    
    ! density
    qm  =0.0_wp
    dqmx=0.0_wp
    dqmy=0.0_wp
    dqmt=0.0_wp
    do nm=1,n_eig   
       phi =eig(nm)%kr*(x(i)-x(1))-eig(nm)%om*trk+eig(nm)%beta*z(k)
       ampl=exp(-eig(nm)%ki*(x(i)-x(1)))

       qm   = (            eig(nm)%rr(j)*cos(phi)-           eig(nm)%ri(j)*sin(phi))*ampl + qm
       dqmx = (-eig(nm)%kr*eig(nm)%rr(j)*sin(phi)-eig(nm)%kr*eig(nm)%ri(j)*cos(phi))*ampl &
            + (            eig(nm)%rr(j)*cos(phi)-           eig(nm)%ri(j)*sin(phi))*(-eig(nm)%ki*ampl) + dqmx
       dqmy = (           eig(nm)%drr(j)*cos(phi)-          eig(nm)%dri(j)*sin(phi))*ampl + dqmy
       dqmt = ( eig(nm)%om*eig(nm)%rr(j)*sin(phi)+eig(nm)%om*eig(nm)%ri(j)*cos(phi))*ampl + dqmt
    enddo
    dr_in=dqmx*cosphi(i,j)+dqmy*sinphi(i,j)+qm*ir(i,j)
    rt_in=dqmt
    !dr_in=0.
    !rt_in=0.

  end subroutine eig_disturb2_imin
             
  !===============================================================================
  subroutine eig_disturb2_imin_jmin(i,j,k,pt_in,ut_in,vt_in,wt_in,rt_in,dp_in,du_in,dv_in,dw_in,dr_in)
  !===============================================================================
    !> Impose eigenmode disturbances at inlet (bottom-left corner: imin_jmin)
  !===============================================================================
    use mod_time 
    implicit none
    ! ----------------------------------------------------------------------------
    integer, intent(in)   :: i,j,k
    ! eigenmode disturbances
    real(wp), intent(out) :: pt_in,ut_in,vt_in,wt_in,rt_in
    real(wp), intent(out) :: dp_in,du_in,dv_in,dw_in,dr_in
    ! ----------------------------------------------------------------------------
    ! local var
    integer :: nm
    real(wp) :: trk,phi,ampl
    real(wp) :: qm,dqmx,dqmy,dqmt
    real(wp), dimension(:,:), pointer :: ir,cosphi,sinphi
    ! ----------------------------------------------------------------------------
    
    ir=>BC_edge(1,1,1)%ir
    cosphi=>BC_edge(1,1,1)%cosphi
    sinphi=>BC_edge(1,1,1)%sinphi

    trk =time+ck(irk)*deltat

    ! pressure
    qm  =0.0_wp
    dqmx=0.0_wp
    dqmy=0.0_wp
    dqmt=0.0_wp
    do nm=1,n_eig   
       phi =eig(nm)%kr*(x(i)-x(1))-eig(nm)%om*trk+eig(nm)%beta*z(k)
       ampl=exp(-eig(nm)%ki*(x(i)-x(1)))

       qm   = (            eig(nm)%pr(j)*cos(phi)-           eig(nm)%pi(j)*sin(phi))*ampl + qm
       dqmx = (-eig(nm)%kr*eig(nm)%pr(j)*sin(phi)-eig(nm)%kr*eig(nm)%pi(j)*cos(phi))*ampl &
            + (            eig(nm)%pr(j)*cos(phi)-           eig(nm)%pi(j)*sin(phi))*(-eig(nm)%ki*ampl) + dqmx
       dqmy = (           eig(nm)%dpr(j)*cos(phi)-          eig(nm)%dpi(j)*sin(phi))*ampl + dqmy
       dqmt = ( eig(nm)%om*eig(nm)%pr(j)*sin(phi)+eig(nm)%om*eig(nm)%pi(j)*cos(phi))*ampl + dqmt
    enddo
    dp_in=dqmx*cosphi(i,j)+dqmy*sinphi(i,j)+qm*ir(i,j)
    pt_in=dqmt
    !dp_in=0.
    !pt_in=0.
    
    ! velocity u
    qm  =0.0_wp
    dqmx=0.0_wp
    dqmy=0.0_wp
    dqmt=0.0_wp
    do nm=1,n_eig   
       phi =eig(nm)%kr*(x(i)-x(1))-eig(nm)%om*trk+eig(nm)%beta*z(k)
       ampl=exp(-eig(nm)%ki*(x(i)-x(1)))

       qm   = (            eig(nm)%ur(j)*cos(phi)-           eig(nm)%ui(j)*sin(phi))*ampl + qm
       dqmx = (-eig(nm)%kr*eig(nm)%ur(j)*sin(phi)-eig(nm)%kr*eig(nm)%ui(j)*cos(phi))*ampl &
            + (            eig(nm)%ur(j)*cos(phi)-           eig(nm)%ui(j)*sin(phi))*(-eig(nm)%ki*ampl) + dqmx
       dqmy = (           eig(nm)%dur(j)*cos(phi)-          eig(nm)%dui(j)*sin(phi))*ampl + dqmy
       dqmt = ( eig(nm)%om*eig(nm)%ur(j)*sin(phi)+eig(nm)%om*eig(nm)%ui(j)*cos(phi))*ampl + dqmt
    enddo
    du_in=dqmx*cosphi(i,j)+dqmy*sinphi(i,j)+qm*ir(i,j)
    ut_in=dqmt
    
    ! velocity v
    qm  =0.0_wp
    dqmx=0.0_wp
    dqmy=0.0_wp
    dqmt=0.0_wp
    do nm=1,n_eig   
       phi =eig(nm)%kr*(x(i)-x(1))-eig(nm)%om*trk+eig(nm)%beta*z(k)
       ampl=exp(-eig(nm)%ki*(x(i)-x(1)))

       qm   = (            eig(nm)%vr(j)*cos(phi)-           eig(nm)%vi(j)*sin(phi))*ampl + qm
       dqmx = (-eig(nm)%kr*eig(nm)%vr(j)*sin(phi)-eig(nm)%kr*eig(nm)%vi(j)*cos(phi))*ampl &
            + (            eig(nm)%vr(j)*cos(phi)-           eig(nm)%vi(j)*sin(phi))*(-eig(nm)%ki*ampl) + dqmx
       dqmy = (           eig(nm)%dvr(j)*cos(phi)-          eig(nm)%dvi(j)*sin(phi))*ampl + dqmy
       dqmt = ( eig(nm)%om*eig(nm)%vr(j)*sin(phi)+eig(nm)%om*eig(nm)%vi(j)*cos(phi))*ampl + dqmt
    enddo
    dv_in=dqmx*cosphi(i,j)+dqmy*sinphi(i,j)+qm*ir(i,j)
    vt_in=dqmt
    
    ! velocity w
    qm  =0.0_wp
    dqmx=0.0_wp
    dqmy=0.0_wp
    dqmt=0.0_wp
    do nm=1,n_eig   
       phi =eig(nm)%kr*(x(i)-x(1))-eig(nm)%om*trk+eig(nm)%beta*z(k)
       ampl=exp(-eig(nm)%ki*(x(i)-x(1)))

       qm   = (            eig(nm)%wr(j)*cos(phi)-           eig(nm)%wi(j)*sin(phi))*ampl + qm
       dqmx = (-eig(nm)%kr*eig(nm)%wr(j)*sin(phi)-eig(nm)%kr*eig(nm)%wi(j)*cos(phi))*ampl &
            + (            eig(nm)%wr(j)*cos(phi)-           eig(nm)%wi(j)*sin(phi))*(-eig(nm)%ki*ampl) + dqmx
       dqmy = (           eig(nm)%dwr(j)*cos(phi)-          eig(nm)%dwi(j)*sin(phi))*ampl + dqmy
       dqmt = ( eig(nm)%om*eig(nm)%wr(j)*sin(phi)+eig(nm)%om*eig(nm)%wi(j)*cos(phi))*ampl + dqmt
    enddo
    dw_in=dqmx*cosphi(i,j)+dqmy*sinphi(i,j)+qm*ir(i,j)
    wt_in=dqmt
    !dw_in=0.
    !wt_in=0.
    
    ! density
    qm  =0.0_wp
    dqmx=0.0_wp
    dqmy=0.0_wp
    dqmt=0.0_wp
    do nm=1,n_eig   
       phi =eig(nm)%kr*(x(i)-x(1))-eig(nm)%om*trk+eig(nm)%beta*z(k)
       ampl=exp(-eig(nm)%ki*(x(i)-x(1)))

       qm   = (            eig(nm)%rr(j)*cos(phi)-           eig(nm)%ri(j)*sin(phi))*ampl + qm
       dqmx = (-eig(nm)%kr*eig(nm)%rr(j)*sin(phi)-eig(nm)%kr*eig(nm)%ri(j)*cos(phi))*ampl &
            + (            eig(nm)%rr(j)*cos(phi)-           eig(nm)%ri(j)*sin(phi))*(-eig(nm)%ki*ampl) + dqmx
       dqmy = (           eig(nm)%drr(j)*cos(phi)-          eig(nm)%dri(j)*sin(phi))*ampl + dqmy
       dqmt = ( eig(nm)%om*eig(nm)%rr(j)*sin(phi)+eig(nm)%om*eig(nm)%ri(j)*cos(phi))*ampl + dqmt
    enddo
    dr_in=dqmx*cosphi(i,j)+dqmy*sinphi(i,j)+qm*ir(i,j)
    rt_in=dqmt
    !dr_in=0.
    !rt_in=0.

  end subroutine eig_disturb2_imin_jmin
             
!!$  !===============================================================
!!$  ! Compute eigenmode disturbances (for Tam & Dong)
!!$  !===============================================================
!!$  subroutine eig_disturb_imin
!!$    use mod_time
!!$    implicit none
!!$    !---------------------------------------------------------------------------
!!$    integer   :: i,j,k
!!$    ! eigenmode disturbances
!!$    real(wp) :: pt_in,ut_in,vt_in,wt_in,rt_in
!!$    real(wp) :: dp_in,du_in,dv_in,dw_in,dr_in
!!$    ! local var
!!$    integer :: nm
!!$    real(wp) :: trk,phi,ampl
!!$    real(wp) :: qm,dqmx,dqmy,dqmt
!!$    real(wp), dimension(:,:), pointer :: ir,cosphi,sinphi,vg
!!$    !---------------------------------------------------------------------------
!!$    
!!$    ir=>BC_face(1,1)%ir
!!$    cosphi=>BC_face(1,1)%cosphi
!!$    sinphi=>BC_face(1,1)%sinphi
!!$    vg=>BC_face(1,1)%vg
!!$
!!$    trk =time+ck(irk)*deltat
!!$
!!$    do i=1,5
!!$       do j=ndy,nfy
!!$          do k=1,nz
!!$             ! pressure
!!$             qm  =0.0_wp
!!$             dqmx=0.0_wp
!!$             dqmy=0.0_wp
!!$             dqmt=0.0_wp
!!$             do nm=1,n_eig   
!!$                phi =eig(nm)%kr*(x(i)-x(1))-eig(nm)%om*trk+eig(nm)%beta*z(k)
!!$                ampl=exp(-eig(nm)%ki*(x(i)-x(1)))
!!$
!!$                qm   = ( eig(nm)%pr(j)*cos(phi)- eig(nm)%pi(j)*sin(phi))*ampl + qm
!!$                dqmx =-( eig(nm)%pr(j)*sin(phi)+ eig(nm)%pi(j)*cos(phi))*eig(nm)%kr*ampl &
!!$                      -( eig(nm)%pr(j)*cos(phi)- eig(nm)%pi(j)*sin(phi))*eig(nm)%ki*ampl + dqmx
!!$                dqmy = (eig(nm)%dpr(j)*cos(phi)-eig(nm)%dpi(j)*sin(phi))*ampl + dqmy
!!$                dqmt = ( eig(nm)%pr(j)*sin(phi)+ eig(nm)%pi(j)*cos(phi))*eig(nm)%om*ampl + dqmt
!!$             enddo
!!$             dp_in=dqmx*cosphi(i,j)+dqmy*sinphi(i,j)+qm*ir(i,j)
!!$             pt_in=dqmt
!!$
!!$             ! velocity u
!!$             qm  =0.0_wp
!!$             dqmx=0.0_wp
!!$             dqmy=0.0_wp
!!$             dqmt=0.0_wp
!!$             do nm=1,n_eig   
!!$                phi =eig(nm)%kr*(x(i)-x(1))-eig(nm)%om*trk+eig(nm)%beta*z(k)
!!$                ampl=exp(-eig(nm)%ki*(x(i)-x(1)))
!!$
!!$                qm   = ( eig(nm)%ur(j)*cos(phi)- eig(nm)%ui(j)*sin(phi))*ampl + qm
!!$                dqmx =-( eig(nm)%ur(j)*sin(phi)+ eig(nm)%ui(j)*cos(phi))*eig(nm)%kr*ampl &
!!$                      -( eig(nm)%ur(j)*cos(phi)- eig(nm)%ui(j)*sin(phi))*eig(nm)%ki*ampl + dqmx
!!$                dqmy = (eig(nm)%dur(j)*cos(phi)-eig(nm)%dui(j)*sin(phi))*ampl + dqmy
!!$                dqmt = ( eig(nm)%ur(j)*sin(phi)+ eig(nm)%ui(j)*cos(phi))*eig(nm)%om*ampl + dqmt
!!$             enddo
!!$             du_in=dqmx*cosphi(i,j)+dqmy*sinphi(i,j)+qm*ir(i,j)
!!$             ut_in=dqmt
!!$
!!$             ! velocity v
!!$             qm  =0.0_wp
!!$             dqmx=0.0_wp
!!$             dqmy=0.0_wp
!!$             dqmt=0.0_wp
!!$             do nm=1,n_eig   
!!$                phi =eig(nm)%kr*(x(i)-x(1))-eig(nm)%om*trk+eig(nm)%beta*z(k)
!!$                ampl=exp(-eig(nm)%ki*(x(i)-x(1)))
!!$
!!$                qm   = ( eig(nm)%vr(j)*cos(phi)- eig(nm)%vi(j)*sin(phi))*ampl + qm
!!$                dqmx =-( eig(nm)%vr(j)*sin(phi)+ eig(nm)%vi(j)*cos(phi))*eig(nm)%kr*ampl &
!!$                      -( eig(nm)%vr(j)*cos(phi)- eig(nm)%vi(j)*sin(phi))*eig(nm)%ki*ampl + dqmx
!!$                dqmy = (eig(nm)%dvr(j)*cos(phi)-eig(nm)%dvi(j)*sin(phi))*ampl + dqmy
!!$                dqmt = ( eig(nm)%vr(j)*sin(phi)+ eig(nm)%vi(j)*cos(phi))*eig(nm)%om*ampl + dqmt
!!$             enddo
!!$             dv_in=dqmx*cosphi(i,j)+dqmy*sinphi(i,j)+qm*ir(i,j)
!!$             vt_in=dqmt
!!$
!!$             ! velocity w
!!$             qm  =0.0_wp
!!$             dqmx=0.0_wp
!!$             dqmy=0.0_wp
!!$             dqmt=0.0_wp
!!$             do nm=1,n_eig   
!!$                phi =eig(nm)%kr*(x(i)-x(1))-eig(nm)%om*trk+eig(nm)%beta*z(k)
!!$                ampl=exp(-eig(nm)%ki*(x(i)-x(1)))
!!$
!!$                qm   = ( eig(nm)%wr(j)*cos(phi)- eig(nm)%wi(j)*sin(phi))*ampl + qm
!!$                dqmx =-( eig(nm)%wr(j)*sin(phi)+ eig(nm)%wi(j)*cos(phi))*eig(nm)%kr*ampl &
!!$                      -( eig(nm)%wr(j)*cos(phi)- eig(nm)%wi(j)*sin(phi))*eig(nm)%ki*ampl + dqmx
!!$                dqmy = (eig(nm)%dwr(j)*cos(phi)-eig(nm)%dwi(j)*sin(phi))*ampl + dqmy
!!$                dqmt = ( eig(nm)%wr(j)*sin(phi)+ eig(nm)%wi(j)*cos(phi))*eig(nm)%om*ampl + dqmt
!!$             enddo
!!$             dw_in=dqmx*cosphi(i,j)+dqmy*sinphi(i,j)+qm*ir(i,j)
!!$             wt_in=dqmt
!!$
!!$             ! density
!!$             qm  =0.0_wp
!!$             dqmx=0.0_wp
!!$             dqmy=0.0_wp
!!$             dqmt=0.0_wp
!!$             do nm=1,n_eig   
!!$                phi =eig(nm)%kr*(x(i)-x(1))-eig(nm)%om*trk+eig(nm)%beta*z(k)
!!$                ampl=exp(-eig(nm)%ki*(x(i)-x(1)))
!!$
!!$                qm   = ( eig(nm)%rr(j)*cos(phi)- eig(nm)%ri(j)*sin(phi))*ampl + qm
!!$                dqmx =-( eig(nm)%rr(j)*sin(phi)+ eig(nm)%ri(j)*cos(phi))*eig(nm)%kr*ampl &
!!$                      -( eig(nm)%rr(j)*cos(phi)- eig(nm)%ri(j)*sin(phi))*eig(nm)%ki*ampl + dqmx
!!$                dqmy = (eig(nm)%drr(j)*cos(phi)-eig(nm)%dri(j)*sin(phi))*ampl + dqmy
!!$                dqmt = ( eig(nm)%rr(j)*sin(phi)+ eig(nm)%ri(j)*cos(phi))*eig(nm)%om*ampl + dqmt
!!$             enddo
!!$             dr_in=dqmx*cosphi(i,j)+dqmy*sinphi(i,j)+qm*ir(i,j)
!!$             rt_in=dqmt
!!$
!!$!!             pt(i,j,k) = pt(i,j,k) - vg(i,j,k)*dp_in - pt_in
!!$!!             ut(i,j,k) = ut(i,j,k) - vg(i,j,k)*du_in - ut_in
!!$!!             vt(i,j,k) = vt(i,j,k) - vg(i,j,k)*dv_in - vt_in
!!$!!             wt(i,j,k) = wt(i,j,k) - vg(i,j,k)*dw_in - wt_in
!!$!!             rt(i,j,k) = rt(i,j,k) - vg(i,j,k)*dr_in - rt_in
!!$          enddo
!!$       enddo
!!$    enddo
!!$    
!!$  end subroutine eig_disturb_imin
             
end module mod_eigenmode
