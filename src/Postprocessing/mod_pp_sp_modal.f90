!==============================================================================
module mod_pp_sp_modal
!==============================================================================
  !> Module for modal analysis of transition
!==============================================================================
  use mod_mpi
  use mod_pp_var
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: nmod      ! number of modes
  integer :: n_period  ! number of time periods
  integer :: n_lambdaz ! number of spanwise wavelengthes
  real(wp) :: om0      ! fundamental frequency
  real(wp) :: beta0    ! fundamental wavenumber
  type mode
     integer :: om,beta
  end type mode
  type(mode), dimension(:), allocatable :: imode ! mode indices
  ! ---------------------------------------------------------------------------

contains

  !============================================================================
  subroutine comp_modes
  !============================================================================
    !> Compute selected modes
  !============================================================================
    use mod_constant  ! <- for twopi
    use mod_grid      ! <- for dimensions
    use mod_time      ! <- for deltat
    use mod_eigenmode ! <- for excitation mode characteristics
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i
    integer :: n_om,n_beta,cpt
    integer, dimension(:), allocatable :: iom,ibeta
    integer :: T_period
    real(wp) :: Lz
    ! -------------------------------------------------------------------------
    
    ! Read unstable mode parameters to define eigenmodes
    ! ==================================================
    call read_param_stab

    ! dimensionalize
    ! --------------
    do i=1,n_eig
       eig(i)%om=eig(i)%om/L_ref*u_ref
       eig(i)%kr=eig(i)%kr/L_ref
       eig(i)%ki=eig(i)%ki/L_ref
       eig(i)%beta=eig(i)%beta/L_ref
    enddo

    om0=eig(1)%om
    beta0=eig(1)%beta

    if (iproc==0) write(*,'(A)') repeat("-",80)

    ! Number of wavelengthes in the spanwise dimension
    ! ================================================
    Lz=zg(ngz)-zg(1)+deltaz ! ???
    Lz=zg(ngz)-zg(1) ! ???
    n_lambdaz=anint(beta0*Lz/twopi)
    if (iproc.eq.0) print *,' ~> n_lambdaz',beta0*Lz/twopi,n_lambdaz

    ! Number of time periods in the time dimension
    ! ============================================
    if (iproc.eq.0) print *,'  T/dt            =', twopi/om0/deltat
    T_period=anint(twopi/om0/deltat/10.0_wp)*10
    deltat=twopi/om0/dble(T_period)
    if (iproc.eq.0) print *,'it per period (T/dt)=', T_period
    if (iproc.eq.0) print *,'new Dt [s]:', deltat
    n_period=anint(dble(nmax)/dble(T_period))
    if (iproc.eq.0) print *,'nb of periods (nmax/T_period)',n_period,dble(nmax)/dble(T_period)
    ngt=n_period*T_period
    if (iproc.eq.0) print *,'ngt=',ngt,'sur nmax=',nmax

    if (iproc.eq.0) print *,'check even number of periods'
    if ((n_period.gt.2).and.(mod(n_period,2).ne.0)) then
       n_period=n_period-1
       if (iproc.eq.0) print *,'nb of periods (nmax/T_period)',n_period,nmax/T_period,dble(nmax)/dble(T_period)
       ngt=n_period*T_period
       if (iproc.eq.0) print *,'ngt=',ngt,'sur nmax=',nmax
    else
       if (iproc.eq.0) print *,'nothing to do'
    endif
    
    if (beta0.ne.0.0_wp) then
       ! mode (0,n_beta)
       n_beta=7
       allocate(ibeta(n_beta))
       ibeta(1)=1
       do i=1,n_beta-1
          ibeta(i+1)=i*n_lambdaz+1
       enddo
    else
       n_beta=1
       allocate(ibeta(n_beta))
       ibeta(1)=1
    endif
    if (iproc.eq.0) print *,'ibeta:',ibeta
    
    ! subharmonic only if n_period is even
    if (mod(n_period,2)==0) then
       n_om=10
       allocate(iom(n_om))
       do i=1,n_om
          iom(i)=i*n_period/2+1
       enddo
    else
       n_om=6
       allocate(iom(n_om))
       do i=1,n_om
          iom(i)=i*n_period+1
       enddo
    endif
    if (iproc.eq.0) print *,'iom:',iom
    
    nmod=0
    ! modes (0,ibeta)
    nmod=nmod+n_beta
    ! modes (iom,0)
    nmod=nmod+n_om
    ! modes (iom,1)
    nmod=nmod+n_om
    ! modes (iom,2)
    nmod=nmod+n_om
    ! modes (iom,3)
    nmod=nmod+n_om

    allocate(imode(nmod))
    cpt=0
    ! modes (0,ibeta)
    do i=1,n_beta
       cpt=cpt+1
       imode(cpt)%om=1
       imode(cpt)%beta=ibeta(i)
    enddo
    ! modes (iom,0)
    do i=1,n_om
       cpt=cpt+1
       imode(cpt)%om=iom(i)
       imode(cpt)%beta=ibeta(1)
    enddo
    ! modes (iom,1)
    do i=1,n_om
       cpt=cpt+1
       imode(cpt)%om=iom(i)
       imode(cpt)%beta=ibeta(2)
    enddo   
    ! modes (iom,2)
    do i=1,n_om
       cpt=cpt+1
       imode(cpt)%om=iom(i)
       imode(cpt)%beta=ibeta(3)
    enddo
    ! modes (iom,3)
    do i=1,n_om
       cpt=cpt+1
       imode(cpt)%om=iom(i)
       imode(cpt)%beta=ibeta(4)
    enddo

    if (iproc==0) print *,'total number of modes',nmod
    if (iproc==0) print *,'omega of modes:',imode%om   
    if (iproc==0) print *,'beta of modes:',imode%beta   
   
  end subroutine comp_modes

end module mod_pp_sp_modal
