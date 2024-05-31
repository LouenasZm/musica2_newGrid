!==============================================================================
module mod_pp_emd
!==============================================================================
  !> author: XG
  !> date: 2020
  !> Module for Empirical Mode Decomposition (EMD)
  !> adapted from C routines in Flandrin et al. EMD package
!==============================================================================
  use precision
  implicit none
  ! ---------------------------------------------------------------------------
  ! default values 
  ! ==============
  integer, parameter :: MAX_ITERATIONS=1000
  integer, parameter :: LIM_GMP=30000
  integer, parameter :: NBSYM=2
  real(wp), parameter :: DEFAULT_THRESHOLD=0.05_wp
  real(wp), parameter :: DEFAULT_TOLERANCE=0.05_wp
  ! ---------------------------------------------------------------------------
  ! derived types
  ! =============
  ! structure for stopping criteria
  type stop_t
     real(wp) :: threshold,tolerance
  end type stop_t
  ! structure used to store an IMF and the associated number of iterations
  type imf_t
     integer :: nb_iterations
     real(wp), dimension(:), allocatable :: next
  end type imf_t
  ! structure of the IMF list
  type imf_list_t
     type(imf_t) :: first
     type(imf_t) :: last
     integer :: m,n
  end type imf_list_t
  ! structure of extrema
  type extrema_t
     integer :: n_min,n_max
     real(wp), dimension(:), allocatable :: x_min,y_min,x_max,y_max
  end type extrema_t
  ! structure of envelope
  type envelope_t
     integer :: n
     real(wp), dimension(:), allocatable :: e_min,e_max,tmp1,tmp2
  end type envelope_t
  ! ---------------------------------------------------------------------------

contains

  !============================================================================
  subroutine emdx(y)
  !============================================================================
    !> main subroutine: call EMD decomposition
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    real(wp), dimension(0:) :: y
    ! -------------------------------------------------------------------------
    integer :: i,n,nb_imfs,max_imfs,iteration_counter
    logical :: stop_status,stop_EMD,allocated_x
    type(extrema_t) :: ex
    type(envelope_t) :: env
    type(stop_t) :: stop_params
    real(wp), dimension(:), allocatable :: x,z,m,a
    ! -------------------------------------------------------------------------

    ! Initializations
    ! ===============   
    ! size of inputs
    n=size(y,1)
    print *,n
    n=n-1

    ! define default parameters
    stop_params.threshold = DEFAULT_THRESHOLD
    stop_params.tolerance = DEFAULT_TOLERANCE
    !error_flag=.false.
    max_imfs=0

    ! define discretization (time or space)
    allocate(x(0:n))
    allocated_x=.true.
    do i=0,n
       x(i) = real(i+1)
    enddo

    ! allocations of memory
    call init_extr(n+2*NBSYM,ex)
    !X call init_imf_list(n,list)
    allocate(z(0:n),m(0:n),a(0:n))
    call init_local_mean(n+2*NBSYM,env)

    ! open files
    open(50,file='data_ex.bin',form='unformatted',status='unknown')
    rewind(50)
    write(50) n+1
    write(50) (x(i),i=0,n)
    write(50) (y(i),i=0,n)

    open(51,file='number_imf.bin',form='unformatted',status='unknown')
    rewind(51)

    open(52,file='data_imf.bin',form='unformatted',status='unknown')
    rewind(52)

    ! Main loop
    ! =========
    print *,'start main loop' 
    nb_imfs=0
    stop_EMD=.false.

    do while ( ((.not.(max_imfs)).or.(nb_imfs < max_imfs)) .and. (.not.(stop_EMD)) )

       ! initialisation
       do i=0,n
          z(i)=y(i)
       enddo
       do i=0,n
          m(i)=y(i)
       enddo
       iteration_counter=0

       call mean_and_amplitude(x,z,m,a,n,ex,env,stop_status)
       write(50) ex%n_min
       write(50) (ex%x_min(i),i=0,ex%n_min-1)
       write(50) (ex%y_min(i),i=0,ex%n_min-1)
       write(50) ex%n_max
       write(50) (ex%x_max(i),i=0,ex%n_max-1)
       write(50) (ex%y_max(i),i=0,ex%n_max-1) 
       write(50) (env%e_max(i),i=0,n)
       write(50) (env%e_min(i),i=0,n)
       write(50) (m(i),i=0,n)

       ! Sifting loop
       !print *,'iteration sifting loop',iteration_counter
       !print *,'stop_status',stop_status
       !print *,'stop_sifting',stop_sifting(m,a,ex,stop_params,n,iteration_counter)
       do while ((.not.(stop_status)).and.(.not.(stop_sifting(m,a,ex,stop_params,n,iteration_counter))))
          ! subtract the local mean
          do i=0,n
             z(i)=z(i)-m(i)
          enddo
          iteration_counter=iteration_counter+1

          call mean_and_amplitude(x,z,m,a,n,ex,env,stop_status)
       enddo ! end sifting loop

       print *,'fin',iteration_counter

       ! save current IMF into list if at least
       ! one sifting iteration has been performed
       if (iteration_counter>0) then
          !X add_imf(&list,z,iteration_counter)
          nb_imfs=nb_imfs+1
          print *,'imf number',nb_imfs
          write(52) (z(i),i=0,n)
          do i=0,n
             y(i)=y(i)-z(i)
          enddo
       else
          stop_EMD=.true.
       endif

    enddo ! end main loop

    ! save the residual into list
    !X add_imf(&list,y,0)

    ! output into a MATLAB array
    !!write_output(list,plhs)
    write(52) (z(i),i=0,n)

    write(51) nb_imfs+1
    ! close files
    close(50)
    close(51)
    close(52)
    !stop

    ! free allocated memory
    if (allocated_x) deallocate(x)
    deallocate(m,a,z)
    call free_local_mean(env)
    !X call free_imf_list(list)
    call free_extr(ex)

  end subroutine emdx

  !============================================================================
  subroutine emd_filter(y)
  !============================================================================
    !> main subroutine: call EMD filtering
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    real(wp), dimension(0:) :: y
    ! -------------------------------------------------------------------------
    integer :: i,n,nb_imfs,max_imfs,iteration_counter
    logical :: stop_status,stop_EMD
    type(extrema_t) :: ex
    type(envelope_t) :: env
    type(stop_t) :: stop_params
    real(wp), dimension(:), allocatable :: x,z,m,a,yf
    ! -------------------------------------------------------------------------

    ! Initializations
    ! ===============   
    ! size of inputs
    n=size(y,1)
    n=n-1

    ! define default parameters
    stop_params.threshold = DEFAULT_THRESHOLD
    stop_params.tolerance = DEFAULT_TOLERANCE
    !error_flag=.false.
    max_imfs=0

    ! define discretization (time or space)
    allocate(x(0:n))
    do i=0,n
       x(i) = real(i+1)
    enddo

    ! allocations of memory
    allocate(z(0:n),m(0:n),a(0:n),yf(0:n))
    yf=0.0_wp  
    call init_extr(n+2*NBSYM,ex)
    call init_local_mean(n+2*NBSYM,env)

    ! Main loop
    ! =========
    !print *,'start EMD filter ...' 
    nb_imfs=0
    stop_EMD=.false.

    do while ( ((.not.(max_imfs)).or.(nb_imfs < max_imfs)) .and. (.not.(stop_EMD)) )

       ! initialisation
       z=y
       m=y
       iteration_counter=0

       call mean_and_amplitude(x,z,m,a,n,ex,env,stop_status)

       ! Sifting loop
       do while ((.not.(stop_status)).and.(.not.(stop_sifting(m,a,ex,stop_params,n,iteration_counter))))
          ! subtract the local mean
          z=z-m
          iteration_counter=iteration_counter+1
          call mean_and_amplitude(x,z,m,a,n,ex,env,stop_status)
       enddo ! end sifting loop

       ! save current IMF into list if at least
       ! one sifting iteration has been performed
       if (iteration_counter>0) then
          nb_imfs=nb_imfs+1
          ! reconstruction
          yf=yf+z
          ! new signal (- imf)
          y=y-z
       else
          stop_EMD=.true.
       endif

    enddo ! end main loop

    ! free allocated memory
    deallocate(x,m,a,z)
    call free_local_mean(env)
    call free_extr(ex)

    ! filtered signal
    y=yf

  end subroutine emd_filter

  !============================================================================
  function emd_fabs(x)
  !============================================================================
    !> Function Absolute value
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    real(wp) :: x,emd_fabs
    ! -------------------------------------------------------------------------

    if (x<0) then
       emd_fabs=-x
    else
       emd_fabs=x
    endif

  end function emd_fabs

  !============================================================================
  function stop_sifting(m,a,ex,sp,n,counter)
  !============================================================================
    !> Stop test for the sifting loop
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: n,counter
    real(wp), dimension(0:) :: m,a
    type(extrema_t) :: ex
    type(stop_t) :: sp
    logical :: stop_sifting
    ! -------------------------------------------------------------------------
    integer :: i,count
    real(wp) :: tol,eps
    ! -------------------------------------------------------------------------

    tol = sp%tolerance*n
    eps = sp%threshold
    count = 0;

    if (counter >= MAX_ITERATIONS) then
       stop_sifting=.true.
       return
    endif
    do i=0,ex%n_min
       if (ex%y_min(i) > 0) then
          stop_sifting=.false.
          return
       endif
    enddo
    do i=0,ex%n_max
       if (ex%y_max(i) < 0) then
          stop_sifting=.false.
          return
       endif
    enddo
    do i=0,n
       if ((emd_fabs(m(i)) > eps*emd_fabs(a(i)))) then
          count=count+1
          if (count>tol) then
             stop_sifting=.false.
             return
          endif
       endif
    enddo
    stop_sifting=.true.

  end function stop_sifting

  !============================================================================
  subroutine init_extr(n,ex)
  !============================================================================
    !> Initialization of extrema structure
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: n
    type(extrema_t) :: ex
    ! -------------------------------------------------------------------------

    allocate(ex%x_min(0:n))
    allocate(ex%y_min(0:n))
    allocate(ex%x_max(0:n))
    allocate(ex%y_max(0:n))

  end subroutine init_extr

  !============================================================================
  subroutine free_extr(ex)
  !============================================================================
    !> Free allocated memory for extrema structure
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
     type(extrema_t) :: ex
   ! -------------------------------------------------------------------------

    deallocate(ex%x_min)
    deallocate(ex%y_min)
    deallocate(ex%x_max)
    deallocate(ex%y_max)

  end subroutine free_extr

  !============================================================================
  subroutine init_local_mean(n,env)
  !============================================================================
    !> Allocate memory for the envelopes and temporary data
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: n
    type(envelope_t) :: env
    ! -------------------------------------------------------------------------

    allocate(env%e_min(0:n))
    allocate(env%e_max(0:n))
    allocate(env%tmp1(0:n))
    allocate(env%tmp2(0:n))

  end subroutine init_local_mean

  !============================================================================
  subroutine free_local_mean(env)
  !============================================================================
    !> Free allocated memory for envelope structure
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    type(envelope_t) :: env
    ! -------------------------------------------------------------------------

    deallocate(env%e_min)
    deallocate(env%e_max)
    deallocate(env%tmp1)
    deallocate(env%tmp2)

  end subroutine free_local_mean

  !============================================================================
  subroutine extr(x,y,n,ex)
  !============================================================================
    !> Detection of local extrema
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: n
    type(extrema_t) :: ex
    real(wp), dimension(:) :: x,y
    ! -------------------------------------------------------------------------
    integer :: cour
    ! -------------------------------------------------------------------------

    ex%n_min=0
    ex%n_max=0

    do cour=2,n
       if ( (y(cour)<=y(cour-1)).and.(y(cour)<=y(cour+1)) ) then ! local minimum
          ex%x_min(ex%n_min+NBSYM)=x(cour)
          ex%y_min(ex%n_min+NBSYM)=y(cour)
          ex%n_min=ex%n_min+1
       endif
       if ( (y(cour)>=y(cour-1)).and.(y(cour)>=y(cour+1)) ) then ! local maximum
          ex%x_max(ex%n_max+NBSYM)=x(cour)
          ex%y_max(ex%n_max+NBSYM)=y(cour)
          ex%n_max=ex%n_max+1
       endif
    enddo

  end subroutine extr

  !============================================================================
  subroutine boundary_conditions(x,y,n,ex)
  !============================================================================
    !> Extrapolation of extrema to limit border effects
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    real(wp), dimension(0:) :: x,y
    integer :: n
    type(extrema_t) :: ex
    ! -------------------------------------------------------------------------
    integer :: cour,nbsym_
    ! -------------------------------------------------------------------------

    nbsym_= NBSYM
    ! reduce the number of symmetrized points if there is not enough extrema
    do while ( (ex%n_min<nbsym_+1).and.(ex%n_max<nbsym_+1) )
       nbsym_=nbsym_-1
    enddo

    if (nbsym_ < NBSYM) then
       do cour=0,ex%n_max
          ex%x_max(nbsym_+cour) = ex%x_max(NBSYM+cour)
          ex%y_max(nbsym_+cour) = ex%y_max(NBSYM+cour)
       enddo
       do cour=0,ex%n_min
          ex%x_min(nbsym_+cour) = ex%x_min(NBSYM+cour)
          ex%y_min(nbsym_+cour) = ex%y_min(NBSYM+cour)
       enddo
    endif

    ! select the symmetrized points and the axis of symmetry at the beginning of the signal
    ! =====================================================================================
    if (ex%x_max(nbsym_) < ex%x_min(nbsym_)) then ! first = max
       if (y(0) > ex%y_min(nbsym_)) then ! the edge is not a min
          if (2*ex%x_max(nbsym_)-ex%x_min(2*nbsym_-1) > x(0)) then ! symmetrized parts are too short
             do cour=0,nbsym_
                ex%x_max(cour) = 2*x(0)-ex%x_max(2*nbsym_-1-cour)
                ex%y_max(cour) = ex%y_max(2*nbsym_-1-cour)
                ex%x_min(cour) = 2*x(0)-ex%x_min(2*nbsym_-1-cour)
                ex%y_min(cour) = ex%y_min(2*nbsym_-1-cour)
             enddo
          else ! symmetrized parts are long enough
             do cour=0,nbsym_
                ex%x_max(cour) = 2*ex%x_max(nbsym_)-ex%x_max(2*nbsym_-cour)
                ex%y_max(cour) = ex%y_max(2*nbsym_-cour)
                ex%x_min(cour) = 2*ex%x_max(nbsym_)-ex%x_min(2*nbsym_-1-cour)
                ex%y_min(cour) = ex%y_min(2*nbsym_-1-cour)
             enddo
          endif
       else ! edge is a min % sym with respect to the edge
          do cour=0,nbsym_
             ex%x_max(cour) = 2*x(0)-ex%x_max(2*nbsym_-1-cour)
             ex%y_max(cour) = ex%y_max(2*nbsym_-1-cour)
          enddo
          do cour=0,nbsym_-1
             ex%x_min(cour) = 2*x(0)-ex%x_min(2*nbsym_-2-cour)
             ex%y_min(cour) = ex%y_min(2*nbsym_-2-cour)
          enddo
          ex%x_min(nbsym_-1) = x(0)
          ex%y_min(nbsym_-1) = y(0)
       endif
    else ! first = min 

       if (y(0) < ex%y_max(nbsym_)) then ! the edge is not a max
          if (2*ex%x_min(nbsym_)-ex%x_max(2*nbsym_-1) > x(0)) then ! symmetrized parts are too short
             do cour=0,nbsym_
                ex%x_max(cour) = 2*x(0)-ex%x_max(2*nbsym_-1-cour)
                ex%y_max(cour) = ex%y_max(2*nbsym_-1-cour)
                ex%x_min(cour) = 2*x(0)-ex%x_min(2*nbsym_-1-cour)
                ex%y_min(cour) = ex%y_min(2*nbsym_-1-cour)
             enddo
          else ! symmetrized parts are long enough
             do cour=0,nbsym_
                ex%x_max(cour) = 2*ex%x_min(nbsym_)-ex%x_max(2*nbsym_-1-cour)
                ex%y_max(cour) = ex%y_max(2*nbsym_-1-cour)
                ex%x_min(cour) = 2*ex%x_min(nbsym_)-ex%x_min(2*nbsym_-cour)
                ex%y_min(cour) = ex%y_min(2*nbsym_-cour)
             enddo
          endif
       else ! edge is a max % sym with respect to the edge
          do cour=0,nbsym_
             ex%x_min(cour) = 2*x(0)-ex%x_min(2*nbsym_-1-cour)
             ex%y_min(cour) = ex%y_min(2*nbsym_-1-cour)
          enddo
          do cour=0,nbsym_-1
             ex%x_max(cour) = 2*x(0)-ex%x_max(2*nbsym_-2-cour)
             ex%y_max(cour) = ex%y_max(2*nbsym_-2-cour)
          enddo
          ex%x_max(nbsym_-1) = x(0)
          ex%y_max(nbsym_-1) = y(0)
       endif
    endif

    ex%n_min= ex%n_min+nbsym_-1
    ex%n_max= ex%n_max+nbsym_-1

    ! select the symmetrized points and the axis of symmetry at the end of the signal
    ! ===============================================================================
    if (ex%x_max(ex%n_max) < ex%x_min(ex%n_min)) then ! last is a min
       if (y(n) < ex%y_max(ex%n_max)) then ! the edge is not a max
          if (2*ex%x_min(ex%n_min)-ex%x_max(ex%n_max-nbsym_+1) < x(n)) then ! symmetrized parts are too short
             do cour=0,nbsym_
                ex%x_max(ex%n_max+1+cour) = 2*x(n)-ex%x_max(ex%n_max-cour)
                ex%y_max(ex%n_max+1+cour) = ex%y_max(ex%n_max-cour)
                ex%x_min(ex%n_min+1+cour) = 2*x(n)-ex%x_min(ex%n_min-cour)
                ex%y_min(ex%n_min+1+cour) = ex%y_min(ex%n_min-cour)
             enddo
          else ! symmetrized parts are long enough
             do cour=0,nbsym_
                ex%x_max(ex%n_max+1+cour) = 2*ex%x_min(ex%n_min)-ex%x_max(ex%n_max-cour)
                ex%y_max(ex%n_max+1+cour) = ex%y_max(ex%n_max-cour)
                ex%x_min(ex%n_min+1+cour) = 2*ex%x_min(ex%n_min)-ex%x_min(ex%n_min-1-cour)
                ex%y_min(ex%n_min+1+cour) = ex%y_min(ex%n_min-1-cour)
             enddo
          endif
       else ! edge is a max % sym with respect to the edge
          do cour=0,nbsym_
             ex%x_min(ex%n_min+1+cour) = 2*x(n)-ex%x_min(ex%n_min-cour)
             ex%y_min(ex%n_min+1+cour) = ex%y_min(ex%n_min-cour)
          enddo
          do cour=0,nbsym_-1
             ex%x_max(ex%n_max+2+cour) = 2*x(n)-ex%x_max(ex%n_max-cour)
             ex%y_max(ex%n_max+2+cour) = ex%y_max(ex%n_max-cour)
          enddo
          ex%x_max(ex%n_max+1) = x(n)
          ex%y_max(ex%n_max+1) = y(n)
       endif
    else ! last is a max
       if (y(n) > ex%y_min(ex%n_min)) then ! the edge is not a min
          if (2*ex%x_max(ex%n_max)-ex%x_min(ex%n_min-nbsym_+1) < x(n)) then ! symmetrized parts are too short
             do cour=0,nbsym_
                ex%x_max(ex%n_max+1+cour) = 2*x(n)-ex%x_max(ex%n_max-cour)
                ex%y_max(ex%n_max+1+cour) = ex%y_max(ex%n_max-cour)
                ex%x_min(ex%n_min+1+cour) = 2*x(n)-ex%x_min(ex%n_min-cour)
                ex%y_min(ex%n_min+1+cour) = ex%y_min(ex%n_min-cour)
             enddo
          else ! symmetrized parts are long enough
             do cour=0,nbsym_
                ex%x_max(ex%n_max+1+cour) = 2*ex%x_max(ex%n_max)-ex%x_max(ex%n_max-1-cour)
                ex%y_max(ex%n_max+1+cour) = ex%y_max(ex%n_max-1-cour)
                ex%x_min(ex%n_min+1+cour) = 2*ex%x_max(ex%n_max)-ex%x_min(ex%n_min-cour)
                ex%y_min(ex%n_min+1+cour) = ex%y_min(ex%n_min-cour)
             enddo
          endif
       else ! edge is a min % sym with respect to the edge
          do cour=0,nbsym_
             ex%x_max(ex%n_max+1+cour) = 2*x(n)-ex%x_max(ex%n_max-cour)
             ex%y_max(ex%n_max+1+cour) = ex%y_max(ex%n_max-cour)
          enddo
          do cour=0,nbsym_-1
             ex%x_min(ex%n_min+2+cour) = 2*x(n)-ex%x_min(ex%n_min-cour)
             ex%y_min(ex%n_min+2+cour) = ex%y_min(ex%n_min-cour)
          enddo
          ex%x_min(ex%n_min+1) = x(n)
          ex%y_min(ex%n_min+1) = y(n)
       endif
    endif

    ex%n_min = ex%n_min + nbsym_ + 1
    ex%n_max = ex%n_max + nbsym_ + 1

  end subroutine boundary_conditions

  !============================================================================
  subroutine interpolation(y,xs,ys,n,x,nx,ys2,temp)
  !============================================================================
    !> Interpolates the sequence (xs,ys) using cubic spline
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: n,nx
    real(wp), dimension(0:) :: x,y,xs,ys,ys2,temp
    ! -------------------------------------------------------------------------
    integer :: i,j,jfin,cur,prev
    real(wp) :: p,sig,a,b,c,d,e,f,g,a0,a1,a2,a3,xc
    ! -------------------------------------------------------------------------

    ! Compute second derivatives at the knots
    ys2(0)=0.0_wp
    temp(0)=0.0_wp
    do i=1,n-1
       sig=(xs(i)-xs(i-1))/(xs(i+1)-xs(i-1))
       p=sig*ys2(i-1)+2.0_wp
       ys2(i)=(sig-1.0_wp)/p
       temp(i)=(ys(i+1)-ys(i))/(xs(i+1)-xs(i))-(ys(i)-ys(i-1))/(xs(i)-xs(i-1))
       temp(i)=(6.0_wp*temp(i)/(xs(i+1)-xs(i-1))-sig*temp(i-1))/p
    enddo
    ys2(n-1)=0.0_wp
    do j=n-2,0,-1
       ys2(j)=ys2(j)*ys2(j+1)+temp(j)
    enddo

    ! Compute the spline coefficients
    cur=0
    j=0
    jfin=n-1
    do while (xs(j+2)<x(0))
       j=j+1
    enddo
    do while (xs(jfin)>x(nx-1))
       jfin=jfin-1
    enddo
    do j=1,jfin
       ! Compute the coefficients of the polynomial between two knots
       a=xs(j)
       b=xs(j+1)
       c=b-a
       d=ys(j)
       e=ys(j+1)
       f=ys2(j)
       g=ys2(j+1)
       a0=(b*d-a*e+b**3*f/6.0_wp-a**3*g/6.0_wp)/c+c*(a*g-b*f)/6.0_wp
       a1=(e-d-b**2*f/2.0_wp+a**2*g/2.0_wp)/c+c*(f-g)/6.0_wp
       a2=(b*f-a*g)/(2.0_wp*c)
       a3=(g-f)/(6.0_wp*c)

       prev=cur
       do while ((cur<nx).and.((j==jfin).or.(x(cur)<xs(j+1))))
          cur=cur+1
       enddo

       ! Compute the value of the spline at the sampling times x(i)
       do i=prev,cur
          xc=x(i)
          y(i)=a0+a1*xc+a2*xc**2+a3*xc**3
       enddo
    enddo

  end subroutine interpolation

  !============================================================================
  subroutine mean_and_amplitude(x,z,m,a,n,ex,env,istop)
  !============================================================================
    !> Computes mean, envelopes and amplitude of the current IMF
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: n
    real(wp), dimension(0:) :: x,z,m,a
    type(extrema_t) :: ex
    type(envelope_t) :: env
    logical :: istop
    ! -------------------------------------------------------------------------
    integer :: i
    ! -------------------------------------------------------------------------

    ! detect maxima and minima
    call extr(x,z,n,ex)

    ! if not enough extrema % stop
    if (ex%n_min+ex%n_max<6) then
       istop=.true.
       return
    endif

    ! add extra points at the edges
    call boundary_conditions(x,z,n,ex)

    ! interpolation - upper envelope
    call interpolation(env%e_max,ex%x_max,ex%y_max,ex%n_max,x,n,env%tmp1,env%tmp2)
    ! interpolation - lower envelope
    call interpolation(env%e_min,ex%x_min,ex%y_min,ex%n_min,x,n,env%tmp1,env%tmp2)

    ! compute the mean
    do i=0,n
       m(i)=(env%e_max(i)+env%e_min(i))/2.0_wp
    enddo

    ! compute the amplitude
    do i=0,n
       a(i)=(env%e_max(i)-env%e_min(i))/2.0_wp
    enddo

    istop=.false.

  end subroutine mean_and_amplitude

end module mod_pp_emd
