!=================================================================================
module mod_analytical_sol
!=================================================================================
  !> Module to compute analytical solution and error
!=================================================================================
  use mod_io        ! <- write analytical solution & error
                    ! <- include mpi,flow,constant,time
  use mod_init_flow ! <- infos for vortex or pulse
  implicit none
  !-------------------------------------------------------------------------------
  complex(wp), parameter :: ii=(0,1)
  !-------------------------------------------------------------------------------

contains
  
  !==============================================================================
  subroutine sol_vortex
  !==============================================================================
    !> Analytical solution for 2-D vortex advection + error computation
  !==============================================================================
    use mod_vortex_model ! <- for vortex models
    use mod_utils        ! <- for numchar
    implicit none
    ! ---------------------------------------------------------------------------
    integer, parameter :: n_var=3 ! 3 variables: u,v velocities and pressure p
    integer :: i,j,k
    real(wp) :: u_advec,v_advec,x_vort,y_vort
    real(wp), dimension(n_var) :: err
    real(wp), dimension(nx,ny) :: u_vort,v_vort,p_vort
    real(wp), dimension(nx,ny,nz,n_var) :: U_an
    ! ---------------------------------------------------------------------------

    ! Advection velocity
    ! ==================
    u_advec=Uc1*u_ref
    v_advec=Uc2*u_ref

    ! New vortex position
    ! ===================
    x_vort=x_vortex+u_advec*time
    y_vort=y_vortex+v_advec*time
   
    ! Compute analytical solution
    ! ===========================
    ! vortex model
    call vortex_2d_model(type_vortex,x_vort,y_vort,u_vort,v_vort,p_vort)
    ! add mean flow
    do k=1,nz
       do j=1,ny
          do i=1,nx
             U_an(i,j,k,1)=u_advec + u_vort(i,j)
             U_an(i,j,k,2)=v_advec + v_vort(i,j)
             U_an(i,j,k,3)=  p_ref + p_vort(i,j)
          enddo
       enddo
    enddo

    ! Compute error (L2-norm)
    ! =======================

    ! local error (per proc)
    err=0.0_wp
    do k=1,nz
       do j=1,ny
          do i=1,nx
             err(1)=err(1)+( uu(i,j,k)-U_an(i,j,k,1))**2
             err(2)=err(2)+( vv(i,j,k)-U_an(i,j,k,2))**2
             err(3)=err(3)+(prs(i,j,k)-U_an(i,j,k,3))**2
          enddo
       enddo
    enddo
    
    ! global error
    call MPI_ALLREDUCE(MPI_IN_PLACE,err,3,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
    err=err/ngx/ngy/ngz

    ! Write error
    ! ===========
    
    ! Write error on proc 0
    ! ---------------------
    if (iproc==0) then
       !open(94,file='err.dat',position='append')
       open(94,file='error.dat')
       rewind(94)
       write(94,'(5(E14.6,1x))') deltax,deltat,err
       print *,'ERROR vortex on u,v,p:',err
    endif

    ! MPI-IO write of volume
    ! ----------------------
    call read_write_vol('analytical_sol_bl'//trim(numchar(nob(iproc)))//filext_write,WRITE,U_an)
    
  end subroutine sol_vortex
  
  !==============================================================================
  subroutine sol_pulse
  !==============================================================================
    !> Analytical solution for 2-D acoustic pulse + error computation
    !> with or w/o wall reflection (on a single plane wall)
    !> (-> AMOS library: double precision) 
  !==============================================================================
    use mod_time      ! <- for deltat
    use mod_utils     ! <- for numchar
    use mod_quadpack  ! <- for numerical quadrature (Gauss-Kronrod)
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,j
    ! pulse
    real(wp) :: ar ! pulse half-width
    real(wp) :: eta
    ! mirror pulse
    logical :: is_wall
    integer :: iw
    real(wp) :: x_mirror,y_mirror,eta_mirror,dist_wall
    real(wp), dimension(:), allocatable :: distw
    ! GK quadrature
    integer, parameter :: key=6
    integer :: neval,ier
    real(wp) :: abserr
    real(wp), parameter :: epsabs=0.0_wp
    real(wp), parameter :: epsrel=0.0001_wp
    real(wp) :: val ! result
    ! solution
    real(wp) :: err ! error
    real(wp), dimension(nx,ny,1,1) :: p_an ! analytical solution
    ! ---------------------------------------------------------------------------
    
    ! Determine if wall reflection is present
    ! =======================================
    is_wall=.false.
    do i=1,2
       do j=1,2
          is_wall=is_wall.or.is_bc_wall(i,j)
       enddo
    enddo
    
    ! Pulse characteristics (defined in param.ini & mod_init_flow)
    ! =====================

    ! Pulse Gaussian half-width
    ar=log(2.0_wp)/(b_pulse*deltax)**2

    ! Find mirror pulse location
    ! ==========================
    if (is_wall) then

       ! left wall
       ! ---------
       if (is_bc_wall(1,1)) then
          
          if (is_curv) then
             ! distances to wall points
             allocate(distw(ngy))
             do j=1,ngy
                distw(j)=sqrt((xgc(1,j)-x_pulse)**2+(ygc(1,j)-y_pulse)**2)
             enddo

             ! index of minimal distance to wall
             iw=minloc(distw,1)
             deallocate(distw)

             ! mirror pulse coordinates
             x_mirror=2.0_wp*xgc(1,iw)-x_pulse
             y_mirror=2.0_wp*ygc(1,iw)-y_pulse
          else
             ! distance to wall
             dist_wall=abs(x_pulse-xg(1))
             
             ! mirror pulse coordinates
             x_mirror=xg(1)-dist_wall
             y_mirror=y_pulse
          endif
          
       ! right wall
       ! ----------
       elseif (is_bc_wall(1,2)) then 
          if (is_curv) then
             ! distances to wall points
             allocate(distw(ngy))
             do j=1,ngy
                distw(j)=sqrt((xgc(ngx,j)-x_pulse)**2+(ygc(ngx,j)-y_pulse)**2)
             enddo

             ! index of minimal distance to wall
             iw=minloc(distw,1)
             deallocate(distw)

             ! mirror pulse coordinates
             x_mirror=2.0_wp*xgc(ngx,iw)-x_pulse
             y_mirror=2.0_wp*ygc(ngx,iw)-y_pulse

          else
             ! distance to wall
             dist_wall=abs(x_pulse-xg(ngx))

             ! mirror pulse coordinates
             x_mirror=xg(ngx)+dist_wall
             y_mirror=y_pulse
          endif

       ! bottom wall
       ! -----------
       elseif (is_bc_wall(2,1)) then
          if (is_curv) then
             ! distances to wall points
             allocate(distw(ngx))
             do i=1,ngx
                distw(i)=sqrt((xgc(i,1)-x_pulse)**2+(ygc(i,1)-y_pulse)**2)
             enddo
             
             ! index of minimal distance to wall
             iw=minloc(distw,1)
             deallocate(distw)

             ! mirror pulse coordinates
             x_mirror=2.0_wp*xgc(iw,1)-x_pulse
             y_mirror=2.0_wp*ygc(iw,1)-y_pulse

          else
             ! distance to wall
             dist_wall=abs(y_pulse-yg(1))
             
             ! mirror pulse coordinates
             x_mirror=x_pulse
             y_mirror=yg(1)-dist_wall
          endif
          
       ! top wall
       ! --------
       elseif (is_bc_wall(2,2)) then ! top wall
          if (is_curv) then
             ! distances to wall points
             allocate(distw(ngx))
             do i=1,ngx
                distw(i)=sqrt((xgc(i,ngy)-x_pulse)**2+(ygc(i,ngy)-y_pulse)**2)
             enddo

             ! index of minimal distance to wall
             iw=minloc(distw,1)
             deallocate(distw)

             ! mirror pulse coordinates
             x_mirror=2.0_wp*xgc(iw,ngy)-x_pulse
             y_mirror=2.0_wp*ygc(iw,ngy)-y_pulse
          else
             ! distance to wall
             dist_wall=abs(y_pulse-yg(ngy))

             ! mirror pulse coordinates
             x_mirror=x_pulse
             y_mirror=yg(ngy)+dist_wall
          endif
       endif
    endif
 
    ! Print infos at screen
    ! ---------------------
    if (iproc==0) then
       print *,'Pulse location x_pulse =',x_pulse
       print *,'               y_pulse =',y_pulse
       if (is_wall) then
          print *,'Pulse location x_mirror=',x_mirror
          print *,'               y_mirror=',y_mirror
       endif
       print *,'Pulse amplitude      A =',ampl_pulse
    endif
    
    ! Compute analytical solution
    ! ===========================

    if (iproc==0) print *,'  ~~> compute & write pulse solution at time ',time
    
    if (is_wall) then
       do i=1,nx
          do j=1,ny
             if (is_curv) then
                eta=sqrt((xc(i,j)-x_pulse-Uc1*u_ref*time)**2+(yc(i,j)-y_pulse-Uc2*u_ref*time)**2)
                eta_mirror=sqrt((xc(i,j)-x_mirror-Uc1*u_ref*time)**2 &
                               +(yc(i,j)-y_mirror-Uc2*u_ref*time)**2)
             else
                eta=sqrt((x(i)-x_pulse-Uc1*u_ref*time)**2+(y(j)-y_pulse-Uc2*u_ref*time)**2)
                eta_mirror=sqrt((x(i)-x_mirror-Uc1*u_ref*time)**2 &
                               +(y(j)-y_mirror-Uc2*u_ref*time)**2)
             endif
             ! call Gauss-Kronrod quadrature (precision and order[key] are defined as parameters)
             call qag(funcw,0.001_wp,10.0_wp,epsabs,epsrel,key,val,abserr,neval,ier)
             p_an(i,j,1,1)=ampl_pulse/(2.0_wp*ar)*val
          enddo
       enddo
    else    
       do i=1,nx
          do j=1,ny
             if (is_curv) then
                eta=sqrt((xc(i,j)-x_pulse-Uc1*u_ref*time)**2+(yc(i,j)-y_pulse-Uc2*u_ref*time)**2)
             else
                eta=sqrt((x(i)-x_pulse-Uc1*u_ref*time)**2+(y(j)-y_pulse-Uc2*u_ref*time)**2)
             endif
             ! call Gauss-Kronrod quadrature (precision and order[key] are defined as parameters)
             call qag(func,0.001_wp,10.0_wp,epsabs,epsrel,key,val,abserr,neval,ier)
             p_an(i,j,1,1)=ampl_pulse/(2.0_wp*ar)*val
          enddo
       enddo
    endif
    
    ! Compute error (L2-norm)
    ! =======================

    ! local error (per proc)
    err=0.0_wp
    do j=1,ny
       do i=1,nx
          err=err+(prs(i,j,1)-p_ref-p_an(i,j,1,1))**2
       enddo
    enddo
    ! global error
    call MPI_ALLREDUCE(MPI_IN_PLACE,err,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
    err=err/ngx/ngy

    ! Write error
    ! ===========
    
    ! Write error on proc 0
    ! ---------------------
    if (iproc==0) then
       !open(94,file='err.dat',position='append')
       open(94,file='error.dat')
       rewind(94)
       write(94,'(3(E14.7,1x))') deltax,deltat,err
       print *,'ERROR pulse on p:',err
    endif

    ! MPI-IO write of volume
    ! ----------------------
    call read_write_vol('analytical_sol_bl'//trim(numchar(nob(iproc)))//filext_write,WRITE,p_an)

  contains

    !============================================================================
    function func(x_)
    !============================================================================
      !> integrand function for pulse analytical solution
    !============================================================================
      !use warnstop
      implicit none
      ! -------------------------------------------------------------------------
      integer :: ifail,n_z
      real(wp) :: func,x_
      real(wp) :: J0r,J0i
      ! -------------------------------------------------------------------------

      ! Compute J_0(Z) (Bessel 1st kind)
      call ZBESJ(x_*eta,0.0_wp,0.0_wp,1,1,J0r,J0i,n_z,ifail)
      !if (ifail>0) call mpiwarn('ATTENTION! Pb convergence Bessel function',0)

      ! Compute integrand
      func=x_*exp(-x_**2/(4.0*ar))*cos(c_ref*time*x_)*J0r

    end function func
      
    !============================================================================
    function funcw(x_)
    !============================================================================
      !> integrand function for pulse analytical solution with wall
    !============================================================================
      !use warnstop
      implicit none
      ! -------------------------------------------------------------------------
      integer :: ifail,n_z
      real(wp) :: funcw,x_
      real(wp) :: J0r,J0m,J0i
      ! -------------------------------------------------------------------------

      ! Compute J_0(Z) (Bessel 1st kind)
      call ZBESJ(x_*eta,0.0_wp,0.0_wp,1,1,J0r,J0i,n_z,ifail)
      call ZBESJ(x_*eta_mirror,0.0_wp,0.0_wp,1,1,J0m,J0i,n_z,ifail)
      !if (ifail>0) call mpiwarn('ATTENTION! Pb convergence Bessel function',0)

      ! Compute integrand
      funcw=x_*exp(-x_**2/(4.0*ar))*cos(c_ref*time*x_)*(J0r+J0m)

    end function funcw
      
  end subroutine sol_pulse
  
!!$  !==============================================================================
!!$  subroutine sol_pulse_old
!!$  !==============================================================================
!!$    !> Analytical solution for 2-D acoustic pulse + error computation
!!$    !> (-> AMOS library: double precision) 
!!$  !==============================================================================
!!$    use mod_time      ! <- for deltat
!!$    use mod_utils     ! <- for numchar
!!$    use mod_integrate ! <- for numerical quadrature (Gauss-Legendre)
!!$    use mod_quadpack  ! <- for numerical quadrature (Gauss-Kronrod)
!!$    implicit none
!!$    ! ---------------------------------------------------------------------------
!!$    integer, parameter :: ndim=5000 ! number of quadrature points
!!$    integer :: i,j,k
!!$    real(wp) :: ar ! pulse half-width
!!$    real(wp) :: eta,arg
!!$    real(wp), dimension(ndim) :: xi,yi
!!$    real(wp) :: err ! error
!!$    real(wp), dimension(nx,ny,1,1) :: p_an ! analytical solution
!!$    ! ---------------------------------------------------------------------------
!!$    !real(wp) :: xl,xu,xm,xc,x_i,y_i
!!$    real(wp):: val ! result
!!$    !real(wp) :: BJ
!!$    integer, parameter :: key=6
!!$    integer :: neval,ier
!!$    real(wp) :: abserr
!!$    real(wp), parameter :: epsabs=0.0_wp
!!$    real(wp), parameter :: epsrel=0.0001_wp
!!$    
!!$    ! Pulse characteristics (defined in param.ini & mod_init_flow)
!!$    ! =====================
!!$
!!$    ! Pulse Gaussian half-width
!!$    ar=log(2.0_wp)/(b_pulse*deltax)**2
!!$
!!$    if (iproc==0) then
!!$       print *,'Pulse location x_pulse=',x_pulse
!!$       print *,'               y_pulse=',y_pulse
!!$       print *,'Pulse amplitude     A =',ampl_pulse
!!$    endif
!!$    
!!$    ! Compute analytical solution
!!$    ! ===========================
!!$
!!$!!    ! definition of integration variable
!!$!!    do k=1,ndim
!!$!!       xi(k)=0.001_wp*dble(k-1)
!!$!!    enddo
!!$!!
!!$!!    ! init Gauss-Legendre coefficients
!!$!!    call init_GL(128)
!!$!!    ! integration bounds
!!$!!    xl=0.001_wp
!!$!!    xu=10.0_wp
!!$!!    ! bounds transformation [l u]->[-1 1]
!!$!!    xm=0.5_wp*(xu-xl)
!!$!!    xc=0.5_wp*(xu+xl)
!!$
!!$    if (iproc==0) print *,'  ~~> compute & write pulse solution at time ',time
!!$
!!$!!    i=80
!!$!!    j=50
!!$!!    do i=1,nx
!!$!!       do j=1,ny
!!$!!    eta=sqrt((x(i)-x_pulse-u_ref*time)**2+(y(j)-y_pulse)**2)
!!$!!    print *,eta
!!$!!          do k=1,ndim
!!$!!             yi(k)=xi(k)*exp(-xi(k)**2/(4.0_wp*ar))    &
!!$!!                  *cos(c_ref*time*xi(k))*BJ(xi(k)*eta,0.0_wp)
!!$!!          enddo
!!$!!          call trapz(ndim,xi,yi,arg)
!!$!!          p_an(i,j,1,1)=ampl_pulse/(2.0_wp*ar)*arg
!!$!!       enddo
!!$!!    enddo
!!$!!          print *,p_an(i,j,1,1)
!!$          
!!$    do i=1,nx
!!$       do j=1,ny
!!$          if (is_curv) then
!!$             eta=sqrt((xc(i,j)-x_pulse-Uc1*u_ref*time)**2+(yc(i,j)-y_pulse-Uc2*u_ref*time)**2)
!!$          else
!!$             eta=sqrt((x(i)-x_pulse-Uc1*u_ref*time)**2+(y(j)-y_pulse-Uc2*u_ref*time)**2)
!!$          endif
!!$          call qag(func,0.001_wp,10.0_wp,epsabs,epsrel,key,val,abserr,neval,ier)
!!$          p_an(i,j,1,1)=ampl_pulse/(2.0_wp*ar)*val
!!$       enddo
!!$    enddo
!!$!!          print *,p_an(i,j,1,1)
!!$!!          print*,'sals'
!!$!!          call qag(func,0.001_wp,10.0_wp,epsabs,epsrel,key,val,abserr,neval,ier)
!!$!!          write(*,'(a,g17.12)') '  Estimated integral is         ', ampl_pulse/(2.*ar)*val
!!$!!          write(*,'(a,g14.6)') '  Estimated integral error =    ', abserr
!!$!!          write(*,'(a,i8)') '  Number of function evaluations, neval = ', neval
!!$!!          write(*,'(a,i8)') '  Error return code ier = ', ier
!!$!!          stop
!!$
!!$!!    do i=1,nx
!!$!!       do j=1,ny
!!$!!          eta=sqrt((x(i)-x_pulse-u_ref*time)**2+(y(j)-y_pulse)**2)
!!$!!          ! G-L quadrature
!!$!!          val=0.0_wp
!!$!!          do k=1,128
!!$!!             !x_i=xl+(xu-xl)*snx(k)
!!$!!             x_i=xm*xabsc(k)+xc
!!$!!             y_i=x_i*exp(-x_i**2/(4.0_wp*ar))    &
!!$!!                    *cos(c_ref*time*x_i)*BJ(x_i*eta,0.0_wp)
!!$!!             !val=val+snwt(k)*y_i
!!$!!             val=val+weig(k)*y_i
!!$!!          enddo
!!$!!          !val=val*(xu-xl)
!!$!!          val=val*xm
!!$!!          p_an(i,j,1,1)=ampl_pulse/(2.0_wp*ar)*val
!!$!!       enddo
!!$!!    enddo
!!$!!      
!!$    ! Compute error (L2-norm)
!!$    ! =======================
!!$
!!$    ! local error (per proc)
!!$    err=0.0_wp
!!$    do j=1,ny
!!$       do i=1,nx
!!$          err=err+(prs(i,j,1)-p_ref-p_an(i,j,1,1))**2
!!$       enddo
!!$    enddo
!!$    ! global error
!!$    call MPI_ALLREDUCE(MPI_IN_PLACE,err,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
!!$    err=err/ngx/ngy
!!$
!!$    ! Write error
!!$    ! ===========
!!$    
!!$    ! Write error on proc 0
!!$    ! ---------------------
!!$    if (iproc==0) then
!!$       !open(94,file='err.dat',position='append')
!!$       open(94,file='error.dat')
!!$       rewind(94)
!!$       write(94,'(3(E14.7,1x))') deltax,deltat,err
!!$       print *,'ERROR pulse on p:',err
!!$    endif
!!$
!!$    ! MPI-IO write of volume
!!$    ! ----------------------
!!$    call read_write_vol('analytical_sol_bl'//trim(numchar(nob(iproc)))//filext_write,WRITE,p_an)
!!$
!!$  contains
!!$
!!$    !============================================================================
!!$    function func(x)
!!$    !============================================================================
!!$      !> integrand function for pulse analytical solution
!!$    !============================================================================
!!$      !use warnstop
!!$      implicit none
!!$      ! -------------------------------------------------------------------------
!!$      integer :: ifail,n_z
!!$      real(wp) :: func,x
!!$      real(wp) :: J0r,J0i
!!$      ! -------------------------------------------------------------------------
!!$
!!$      ! Compute J_0(Z) (Bessel 1st kind)
!!$      call ZBESJ(x*eta,0.0_wp,0.0_wp,1,1,J0r,J0i,n_z,ifail)
!!$      !if (ifail>0) call mpiwarn('ATTENTION! Pb convergence Bessel function',0)
!!$
!!$      ! Compute integrand
!!$      func=x*exp(-x**2/(4.0*ar))*cos(c_ref*time*x)*J0r
!!$
!!$    end function func
!!$      
!!$  end subroutine sol_pulse_old
!!$
!!$  !==============================================================================
!!$  subroutine sol_pulse_wall_old
!!$  !==============================================================================
!!$    !> Analytical solution for 2-D acoustic pulse + error computation
!!$    !> (-> AMOS library: double precision) 
!!$  !==============================================================================
!!$    use mod_time      ! <- for deltat
!!$    use mod_utils     ! <- for numchar
!!$    use mod_integrate ! <- for numerical quadrature (Gauss-Legendre)
!!$    use mod_quadpack  ! <- for numerical quadrature (Gauss-Kronrod)
!!$    implicit none
!!$    ! ---------------------------------------------------------------------------
!!$    integer, parameter :: ndim=5000 ! number of quadrature points
!!$    integer :: i,j,k
!!$    real(wp) :: ar ! pulse half-width
!!$    real(wp) :: eta,arg
!!$    ! mirror source
!!$    real(wp) :: x_mirror,y_mirror,eta_mirror,dist_wall
!!$    real(wp), dimension(ndim) :: xi,yi
!!$    real(wp) :: err ! error
!!$    real(wp), dimension(nx,ny,1,1) :: p_an ! analytical solution
!!$    ! ---------------------------------------------------------------------------
!!$    !real(wp) :: xl,xu,xm,xc,x_i,y_i
!!$    real(wp):: val ! result
!!$    !real(wp) :: BJ
!!$    integer, parameter :: key=6
!!$    integer :: neval,ier
!!$    real(wp) :: abserr
!!$    real(wp), parameter :: epsabs=0.0_wp
!!$    real(wp), parameter :: epsrel=0.0001_wp
!!$    integer :: iw
!!$    real(wp), dimension(ngx) :: distw
!!$    
!!$    ! Pulse characteristics (defined in param.ini & mod_init_flow)
!!$    ! =====================
!!$
!!$    ! Pulse Gaussian half-width
!!$    ar=log(2.0_wp)/(b_pulse*deltax)**2
!!$
!!$    x_mirror=x_pulse
!!$    dist_wall=abs(y_pulse-yg(1))
!!$    y_mirror=yg(1)-dist_wall
!!$
!!$    ! distances to wall points
!!$    do i=1,ngx
!!$       distw(i)=sqrt((xgc(i,1)-x_pulse)**2+(ygc(i,1)-y_pulse)**2)
!!$    enddo
!!$    ! minimal distance to wall
!!$    iw=minloc(distw,1)
!!$    print *,'iw',iw
!!$    dist_wall=distw(iw)
!!$    print *,'dist_wall',dist_wall
!!$    x_mirror=2.0_wp*xgc(iw,1)-x_pulse
!!$    y_mirror=2.0_wp*ygc(iw,1)-y_pulse
!!$    
!!$    if (iproc==0) then
!!$       print *,'Pulse location x_pulse =',x_pulse
!!$       print *,'               y_pulse =',y_pulse
!!$       print *,'Pulse location x_mirror=',x_mirror
!!$       print *,'               y_mirror=',y_mirror
!!$       print *,'Pulse amplitude      A =',ampl_pulse
!!$    endif
!!$    
!!$    ! Compute analytical solution
!!$    ! ===========================
!!$
!!$!!    ! definition of integration variable
!!$!!    do k=1,ndim
!!$!!       xi(k)=0.001_wp*dble(k-1)
!!$!!    enddo
!!$!!
!!$!!    ! init Gauss-Legendre coefficients
!!$!!    call init_GL(128)
!!$!!    ! integration bounds
!!$!!    xl=0.001_wp
!!$!!    xu=10.0_wp
!!$!!    ! bounds transformation [l u]->[-1 1]
!!$!!    xm=0.5_wp*(xu-xl)
!!$!!    xc=0.5_wp*(xu+xl)
!!$
!!$    if (iproc==0) print *,'  ~~> compute & write pulse solution at time ',time
!!$
!!$!!    do i=1,nx
!!$!!       do j=1,ny
!!$!!          eta=sqrt((x(i)-x_pulse-u_ref*time)**2+(y(j)-y_pulse)**2)
!!$!!          eta_mirror=sqrt((x(i)-x_mirror-u_ref*time)**2+(y(j)-y_mirror)**2)
!!$!!          do k=1,ndim
!!$!!             yi(k)=xi(k)*exp(-xi(k)**2/(4.0_wp*ar))    &
!!$!!                  *cos(c_ref*time*xi(k))*(BJ(xi(k)*eta,0.0_wp)+BJ(xi(k)*eta_mirror,0.0_wp))
!!$!!          enddo
!!$!!          call trapz(ndim,xi,yi,arg)
!!$!!          p_an(i,j,1,1)=ampl_pulse/(2.0_wp*ar)*arg
!!$!!       enddo
!!$!!    enddo
!!$
!!$    do i=1,nx
!!$       do j=1,ny
!!$          if (is_curv) then
!!$             eta=sqrt((xc(i,j)-x_pulse-Uc1*u_ref*time)**2+(yc(i,j)-y_pulse-Uc2*u_ref*time)**2)
!!$             eta_mirror=sqrt((xc(i,j)-x_mirror-Uc1*u_ref*time)**2+(yc(i,j)-y_mirror-Uc2*u_ref*time)**2)
!!$          else
!!$             eta=sqrt((x(i)-x_pulse-Uc1*u_ref*time)**2+(y(j)-y_pulse-Uc2*u_ref*time)**2)
!!$             eta_mirror=sqrt((x(i)-x_mirror-Uc1*u_ref*time)**2+(y(j)-y_mirror-Uc2*u_ref*time)**2)
!!$          endif
!!$          call qag(func,0.001_wp,10.0_wp,epsabs,epsrel,key,val,abserr,neval,ier)
!!$          p_an(i,j,1,1)=ampl_pulse/(2.0_wp*ar)*val
!!$       enddo
!!$    enddo
!!$
!!$!!    do i=1,nx
!!$!!       do j=1,ny
!!$!!          eta=sqrt((x(i)-x_pulse-u_ref*time)**2+(y(j)-y_pulse)**2)
!!$!!          eta_mirror=sqrt((x(i)-x_mirror-u_ref*time)**2+(y(j)-y_mirror)**2)
!!$!!          ! G-L quadrature
!!$!!          val=0.0_wp
!!$!!          do k=1,128
!!$!!             !x_i=xl+(xu-xl)*snx(k)
!!$!!             x_i=xm*xabsc(k)+xc
!!$!!             y_i=x_i*exp(-x_i**2/(4.0_wp*ar))    &
!!$!!                    *cos(c_ref*time*x_i)*(BJ(x_i*eta,0.0_wp)+BJ(x_i*eta_mirror,0.0_wp))
!!$!!             !val=val+snwt(k)*y_i
!!$!!             val=val+weig(k)*y_i
!!$!!          enddo
!!$!!          !val=val*(xu-xl)
!!$!!          val=val*xm
!!$!!          p_an(i,j,1,1)=ampl_pulse/(2.0_wp*ar)*val
!!$!!       enddo
!!$!!    enddo
!!$      
!!$    ! Compute error (L2-norm)
!!$    ! =======================
!!$
!!$    ! local error (per proc)
!!$    err=0.0_wp
!!$    do j=1,ny
!!$       do i=1,nx
!!$          err=err+(prs(i,j,1)-p_ref-p_an(i,j,1,1))**2
!!$       enddo
!!$    enddo
!!$    ! global error
!!$    call MPI_ALLREDUCE(MPI_IN_PLACE,err,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
!!$    err=err/ngx/ngy
!!$
!!$    ! Write error
!!$    ! ===========
!!$    
!!$    ! Write error on proc 0
!!$    ! ---------------------
!!$    if (iproc==0) then
!!$       !open(94,file='err.dat',position='append')
!!$       open(94,file='error.dat')
!!$       rewind(94)
!!$       write(94,'(3(E14.7,1x))') deltax,deltat,err
!!$       print *,'ERROR pulse on p:',err
!!$    endif
!!$
!!$    ! MPI-IO write of volume
!!$    ! ----------------------
!!$    call read_write_vol('analytical_sol_bl'//trim(numchar(nob(iproc)))//filext_write,WRITE,p_an)
!!$
!!$  contains
!!$
!!$    !============================================================================
!!$    function func(x)
!!$    !============================================================================
!!$      !> integrand function for pulse analytical solution
!!$    !============================================================================
!!$      !use warnstop
!!$      implicit none
!!$      ! -------------------------------------------------------------------------
!!$      integer :: ifail,n_z
!!$      real(wp) :: func,x
!!$      real(wp) :: J0r,J0m,J0i
!!$      ! -------------------------------------------------------------------------
!!$
!!$      ! Compute J_0(Z) (Bessel 1st kind)
!!$      call ZBESJ(x*eta,0.0_wp,0.0_wp,1,1,J0r,J0i,n_z,ifail)
!!$      call ZBESJ(x*eta_mirror,0.0_wp,0.0_wp,1,1,J0m,J0i,n_z,ifail)
!!$      !if (ifail>0) call mpiwarn('ATTENTION! Pb convergence Bessel function',0)
!!$
!!$      ! Compute integrand
!!$      func=x*exp(-x**2/(4.0*ar))*cos(c_ref*time*x)*(J0r+J0m)
!!$
!!$    end function func
!!$      
!!$  end subroutine sol_pulse_wall_old

!!$  !==============================================================================
!!$  function BJ(Z,fnu)
!!$  !==============================================================================
!!$    !> Pulse solution: real function to compute Bessel J of order nu
!!$  !==============================================================================
!!$    use warnstop
!!$    implicit none
!!$    ! ---------------------------------------------------------------------------
!!$    real(wp), intent(in) :: Z,fnu ! argument & order
!!$    real(wp) :: BJ ! OUTPUT
!!$    ! ---------------------------------------------------------------------------
!!$    integer :: ifail,n_z
!!$    real(wp) :: cyr,cyi
!!$    ! ---------------------------------------------------------------------------
!!$
!!$    ! Compute J_nu(Z) (Bessel 1st kind)
!!$    ! ===============
!!$    call ZBESJ(Z,0.,fnu,1,1,cyr,cyi,n_z,ifail)
!!$    BJ=cyr
!!$
!!$    if (ifail>0) call mpiwarn('ATTENTION! Pb convergence Bessel function',0)
!!$
!!$  end function BJ

!!$  !==============================================================================
!!$  real(wp) function BJ(Z,fnu)
!!$  !==============================================================================
!!$    !> Pulse solution: real function to compute Bessel J of order nu
!!$  !==============================================================================
!!$    use warnstop
!!$    implicit none
!!$    ! ---------------------------------------------------------------------------
!!$    real(wp), intent(in) :: Z,fnu
!!$    ! ---------------------------------------------------------------------------
!!$    integer :: ifail,nz
!!$    real(wp) :: cyr,cyi
!!$    ! ---------------------------------------------------------------------------
!!$
!!$    ! Compute J_nu(Z) (Bessel 1st kind)
!!$    ! ===============
!!$    call ZBESJ(Z,0.,fnu,1,1,cyr,cyi,nz,ifail)
!!$    BJ=cyr
!!$
!!$    if (ifail>0) call mpiwarn('ATTENTION! Pb convergence Bessel function',0)
!!$    
!!$  end function BJ
!!$
!!$  !==============================================================================
!!$  complex(wp) function J_nu(Z,fnu)
!!$  !==============================================================================
!!$    !> Pulse solution: complex function to compute Bessel J of order nu
!!$  !==============================================================================
!!$    use warnstop
!!$    implicit none
!!$    ! ---------------------------------------------------------------------------
!!$    complex(wp), intent(in) :: Z
!!$    real(wp), intent(in) :: fnu
!!$    ! ---------------------------------------------------------------------------
!!$    integer, parameter :: n=1
!!$    integer :: m,ifail,n_Z
!!$    complex(wp), dimension(n) :: cy
!!$    !real(wp) :: Zr,Zi
!!$    !real(wp), dimension(n) :: cyr,cyi
!!$    ! ---------------------------------------------------------------------------
!!$
!!$    ! Compute J_0(Z) and J_1(Z) (Bessel 1st kind)
!!$    ! =========================
!!$    !fnu=0.00
!!$    !Zr=real(Z)
!!$    !Zi=aimag(Z)
!!$
!!$    call CBESJ(Z,fnu,1,n,cy,nz,ifail)
!!$    !print *,Z,cy
!!$    J_nu=cy(1)
!!$    !J_1=cy(2)
!!$    !call CBESJ(Zr,Zi,fnu,1,n,cyr,cyi,nz,ifail)
!!$    !J_0=cyr(1)+ii*cyi(1)
!!$    !J_1=cyr(2)+ii*cyi(2)
!!$
!!$    if (ifail>0) call mpiwarn('ATTENTION! Pb convergence Bessel function',0)
!!$    
!!$  end function J_nu
!!$
!!$  !==============================================================================
!!$  subroutine besselJ(Z,J_0,J_1)
!!$  !==============================================================================
!!$    !> Pulse solution: compute Bessel function (Bessel of first kind)
!!$  !==============================================================================
!!$    use warnstop
!!$    implicit none
!!$    ! ---------------------------------------------------------------------------
!!$    complex(wp), intent(in) :: Z
!!$    complex(wp), intent(out) :: J_0,J_1
!!$    ! ---------------------------------------------------------------------------
!!$    integer, parameter :: n=1
!!$    integer :: m,ifail,n_Z
!!$    real(wp) :: fnu,Zr,Zi
!!$    real(wp), dimension(n) :: cyr,cyi
!!$    complex(wp), dimension(n) :: cy
!!$    ! ---------------------------------------------------------------------------
!!$
!!$    ! Compute J_0(Z) and J_1(Z) (Bessel 1st kind)
!!$    ! =========================
!!$    fnu=0.0_wp
!!$    ifail=0
!!$    Zr=real(Z)
!!$    Zi=aimag(Z)
!!$
!!$    call CBESJ(Z,fnu,1,n,cy,n_Z,ifail)
!!$    !print *,Z,cy
!!$    J_0=cy(1)
!!$    J_1=cy(2)
!!$    !call ZBESJ(Zr,Zi,fnu,1,n,cyr,cyi,n_Z,ifail)
!!$    !J_0=cyr(1)+ii*cyi(1)
!!$    !J_1=cyr(2)+ii*cyi(2)
!!$
!!$    if (ifail>0) call mpiwarn('ATTENTION! Pb convergence Bessel function',0)
!!$
!!$  end subroutine besselJ
!!$
!!$  !==============================================================================
!!$  subroutine besselY(Z,Y_0,Y_1)
!!$  !==============================================================================
!!$    !> Pulse solution: compute Neumann function (Bessel of second kind)
!!$  !==============================================================================
!!$    use warnstop
!!$    implicit none
!!$    ! ---------------------------------------------------------------------------
!!$    complex(wp), intent(in) :: Z
!!$    complex(wp), intent(out) :: Y_0,Y_1
!!$    ! ---------------------------------------------------------------------------
!!$    integer, parameter :: n=2
!!$    integer :: m,ifail,n_Z
!!$    real(wp) :: fnu,Zr,Zi
!!$    real(wp), dimension(n) :: cyr,cyi
!!$    real(wp), dimension(n+2) :: cwrkr,cwrki
!!$    ! ---------------------------------------------------------------------------
!!$
!!$    ! Compute Y_0(Z) and Y_1(Z) (Bessel 2nd kind)
!!$    ! =========================
!!$    fnu=0.0_wp
!!$    ifail=0
!!$    Zr=real(Z)
!!$    Zi=aimag(Z)
!!$
!!$    call ZBESY(Zr,Zi,fnu,1,n,cyr,cyi,n_Z,cwrkr,cwrki,ifail)
!!$    Y_0=cyr(1)+ii*cyi(1)
!!$    Y_1=cyr(2)+ii*cyi(2)
!!$
!!$    if (ifail>0) call mpiwarn('ATTENTION! Pb convergence Neumann function',0)
!!$
!!$  end subroutine besselY
!!$
!!$  !==============================================================================
!!$  subroutine hankel(Z,H2_0,H2_1)
!!$  !==============================================================================
!!$    !> Pulse solution: compute Hankel function of second kind
!!$  !==============================================================================
!!$    use warnstop
!!$    implicit none
!!$    ! ---------------------------------------------------------------------------
!!$    complex(wp), intent(in) :: Z
!!$    complex(wp), intent(out) :: H2_0,H2_1
!!$    ! ---------------------------------------------------------------------------
!!$    integer, parameter :: n=2
!!$    integer :: m,ifail,n_Z
!!$    real(wp) :: fnu,Zr,Zi
!!$    complex(wp) :: val
!!$    real(wp), dimension(n) :: cyr,cyi
!!$    ! ---------------------------------------------------------------------------
!!$
!!$    ! Compute H2_0(Z) and H2_1(Z) (2nd kind)
!!$    ! ===========================
!!$    m=2
!!$    fnu=0.0_wp
!!$    ifail=0
!!$    Zr=real(Z)
!!$    Zi=aimag(Z)
!!$
!!$    if (abs(Z)<100.0_wp) then
!!$       call ZBESH(Zr,Zi,fnu,1,m,n,cyr,cyi,n_Z,ifail)
!!$       H2_0=cyr(1)+ii*cyi(1)
!!$       H2_1=cyr(2)+ii*cyi(2)
!!$    else
!!$       val=sqrt(2.0_wp/(pi*Z))*exp(-ii*(Z-pi/4.0_wp))
!!$       H2_0=val
!!$       H2_1=ii*val
!!$    endif
!!$    
!!$    if (ifail>0) call mpiwarn('ATTENTION! Pb convergence Hankel function',0)
!!$
!!$  end subroutine hankel

end module mod_analytical_sol
