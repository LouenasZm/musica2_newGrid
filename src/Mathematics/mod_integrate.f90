!=================================================================================
module mod_integrate
!=================================================================================
  !> Module for spatial integrations (quadrature)
!=================================================================================
  use precision
  implicit none
  !-------------------------------------------------------------------------------
  ! Gauss-Legendre quadrature coefficients
  ! --------------------------------------
  integer, private :: ngp ! # of Gauss Points
  real(wp), dimension(:), allocatable :: xabsc,weig ! abscissae & weights
  ! 16-point version
  real(wp), dimension(8) :: ti,ci ! # of Gauss points/2
  ! 32-point version
  real(wp), dimension(32) :: snwt,snx
  !-------------------------------------------------------------------------------

contains
  
!!$  !===============================================================================
!!$  subroutine init_integ(dim)
!!$  !===============================================================================
!!$    !> Initialize integration elements for trapezoidal rule
!!$  !===============================================================================
!!$    use mod_grid
!!$    implicit none
!!$    ! ----------------------------------------------------------------------------
!!$    integer :: dim ! direction for which the integration is applied
!!$    ! ----------------------------------------------------------------------------
!!$    integer :: i,j,k
!!$    real(wp), dimension(:), allocatable :: dygi,dzgi ! integration element
!!$    real(wp) :: ro
!!$    ! ----------------------------------------------------------------------------
!!$
!!$    ! Integration elements at inlet (to compute bulk quantities)
!!$    ! =============================
!!$    allocate(dygi(ngy))
!!$    if (is_curv) then
!!$       dygi(1)=(ygc(1,2)-ygc(1,1))*0.5_wp
!!$       do j=2,ngy-1
!!$          dygi(j)=(ygc(1,j+1)-ygc(1,j-1))*0.5_wp
!!$       enddo
!!$       dygi(ngy)=(ygc(1,ngy)-ygc(1,ngy-1))*0.5_wp
!!$    else
!!$       dygi(1)=(yg(2)-yg(1))*0.5_wp
!!$       do j=2,ngy-1
!!$          dygi(j)=(yg(j+1)-yg(j-1))*0.5_wp
!!$       enddo
!!$       dygi(ngy)=(yg(ngy)-yg(ngy-1))*0.5_wp
!!$    endif
!!$
!!$    allocate(dzgi(ngz))
!!$    dzgi=deltaz
!!$
!!$    ! Partitioning of integration elements
!!$    ! ====================================
!!$    if (allocated(dyi)) deallocate(dyi)
!!$    allocate(dyi(ny))
!!$    do j=1,ny
!!$       dyi(j)=dygi(j+coord(2)*ny)
!!$    enddo
!!$
!!$    if (allocated(dzi)) deallocate(dzi)
!!$    allocate(dzi(nz))
!!$    do k=1,nz
!!$       dzi(k)=dzgi(k+coord(3)*nz)
!!$    enddo
!!$
!!$  end subroutine init_integ

  !==============================================================================
  subroutine trapz(ni,x,f,xint)
  !==============================================================================
    !> Pulse solution: numerical quadrature (trapezoidal rule)
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer, intent(in) :: ni
    real(wp), intent(in) :: x(ni),f(ni) ! abscissae, function
    real(wp), intent(out) :: xint       ! result
    ! ---------------------------------------------------------------------------
    integer :: i
    ! ---------------------------------------------------------------------------

    ! initialization
    xint=0.0_wp

    ! Trapezoidal rule
    do i=1,ni-1
       xint=xint+(f(i)+f(i+1))*(x(i+1)-x(i))*0.5_wp
    enddo
    
  end subroutine trapz
  
  !===============================================================================
  subroutine init_GL(n_gl)
  !===============================================================================
    !> Initialization of Gauss-Legendre abscissas and weights
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer, intent(in) :: n_gl ! # of Gauss Points
    ! ----------------------------------------------------------------------------

    if (n_gl==16) then
       ngp=8
       ti=[0.0950125098, 0.2816035507, 0.4580167776, 0.6178762444, &  
           0.7554044083, 0.8656312023, 0.9445750230, 0.9894009349]
       ci=[0.1894506104, 0.1826034150, 0.1691565193, 0.1495959888, &
           0.1246289712, 0.0951585116, 0.0622535239, 0.0271524594]
    elseif (n_gl==32) then
       ngp=32
       snx=[.001368073,.0071942277,.0176188818,.0325469453,&
            .0518394344,.0753161897,.1027581033,.1339089403,.1684778666,&
            .2061421214,.2465500455,.2893243619,.3340656989,.3803563189,&
            .4277640192,.4758461672,.5241538328,.5722359808,.6196436811,&
            .6659343011,.7106756381,.7534499545,.7938578786,.8315221334,&
            .8660910597,.8972418967,.9246838103,.9481605656,.9674530547,&
            .9823811182,.9928057723,.9986319268]
       snwt=[.003509312,.0081372115,.0126959771,.0171370005,&
            .0214179378,.0254990642,.0293420289,.0329111118,.0361728923,&
            .0390969435,.0416559573,.0438260415,.0455869341,.0469221942,&
            .0478193546,.0482700387,.0482700387,.0478193546,.0469221942,&
            .0455869341,.0438260415,.0416559573,.0390969435,.0361728923,&
            .0329111118,.0293420289,.0254990642,.0214179378,.0171370005,&
            .0126959771,.0081372115,.0035093120]
    else
       ngp=n_gl
       allocate(xabsc(ngp),weig(ngp))
       ! compute coefficients & weights
       call gauleg
    endif
      
  end subroutine init_GL
  
  !===============================================================================
  subroutine gauleg
  !===============================================================================
    !> Calculation of Gauss-Legendre abscissas and weights for Gaussian quadrature
    !> integration of polynomial functions.
    ! ----------------------------------------------------------------------------
    ! For normalized lower and upper limits of integration -1.0 & 1.0, and given n,
    ! this routine calculates, arrays xabsc(1:n) and  weig(1:n) of length n,
    ! containing the abscissas and weights of the Gauss-Legendre n-point quadrature
    ! formula.  For detailed explanations finding weights & abscissas, see
    ! "Numerical Recipes in Fortran */ [adapté de gaussm3.f90]
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,m
    real(wp) :: p1,p2,p3,pp,z,z1
    real(wp), parameter  :: eps=3.0e-15_wp,M_PI=3.141592654_wp
    ! ----------------------------------------------------------------------------

    m=(ngp+1)/2
    ! Roots are symmetric in the interval - so only need to find half of them

    do i=1,m ! Loop over the desired roots

       z=cos(M_PI*(i-0.25_wp)/(ngp+0.5_wp))
       ! Starting with the above approximation to the ith root,
       ! we enter the main loop of refinement by Newton'S method
100    p1=1.0_wp
       p2=0.0_wp
       ! Loop up the recurrence relation to get the Legendre
       ! polynomial evaluated at z 
       do j=1,ngp
          p3=p2
          p2=p1
          p1=((2.0_wp*j-1.0_wp)*z*p2-(j-1.0_wp)*p3)/dble(j)
       enddo

       ! p1 is now the desired Legendre polynomial. We next compute pp,
       ! its derivative, by a standard relation involving also p2, the
       ! polynomial of one lower order.
       pp= ngp*(z*p1-p2)/(z*z-1.0_wp)
       z1= z
       z = z1-p1/pp ! Newton's Method

       if (abs(z-z1).gt.eps) goto 100

       xabsc(i)=-z                     ! Roots will be bewteen -1.0 & 1.0
       xabsc(ngp+1-i)=z                ! and symmetric about the origin
       weig(i)=2.0_wp/((1.-z*z)*pp*pp) ! Compute the weight and its
       weig(ngp+1-i)=weig(i)           ! symmetric counterpart

    enddo ! i loop

  end subroutine gauleg

  !===============================================================================
  subroutine integ_GL(xl,xu,f,val)
  !===============================================================================
    !> Numerical integration - Gauss-Legendre quadrature (general complex version)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: xl,xu ! integration bounds
    complex(wp), intent(out):: val ! result
    complex(wp) :: f ! function to integrate
    external :: f
    ! ----------------------------------------------------------------------------
    integer(wp) :: i
    real(wp) :: xm,xc
    complex(wp) :: x
    ! ----------------------------------------------------------------------------

    ! initialization
    val=0.0_wp

    ! bounds transformation [l u]->[-1 1]
    xm=0.5_wp*(xu-xl)
    xc=0.5_wp*(xu+xl)
    
    ! G-L quadrature
    do i=1,ngp
       x=xm*xabsc(i)+xc
       val=val+weig(i)*f(x)
    enddo
    val=val*xm

  end subroutine integ_GL

  !===============================================================================
  subroutine integ_GL16(xl,xu,f,val)
  !===============================================================================
    !> Numerical integration - Gauss-Legendre quadrature (16-point real version)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: xl,xu ! integration bounds
    real(wp), intent(out):: val ! result
    real(wp) :: f ! function to integrate
    external :: f
    ! ----------------------------------------------------------------------------
    integer(wp) :: i
    real(wp) :: xm,xc,x1,x2
    ! ----------------------------------------------------------------------------
        
    ! initialization
    val=0.0_wp

    ! bounds transformation [l u]->[-1 1]
    xm=0.5_wp*(xu-xl)
    xc=0.5_wp*(xu+xl)
    
    ! G-L quadrature
    do i=1,ngp
       x1=-xm*ti(i)+xc
       x2= xm*ti(i)+xc
       val=val+ci(i)*(f(x1)+f(x2))
    end do
    val=val*xm

  end subroutine integ_GL16

  !===============================================================================
  subroutine integ_GL32(xl,xu,f,val)
  !===============================================================================
    !> Numerical integration - Gauss-Legendre quadrature (32-point real version)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in) :: xl,xu ! integration bounds
    real(wp), intent(out):: val ! result
    real(wp) :: f ! function to integrate
    external :: f
    ! ----------------------------------------------------------------------------
    integer(wp) :: i
    real(wp) :: x
    ! ----------------------------------------------------------------------------

    ! initialization
    val=0.0_wp
    
    ! G-L quadrature
    do i=1,ngp
       x=xl+(xu-xl)*snx(i)
       val=val+snwt(i)*f(x)
    enddo
    val=val*(xu-xl)

  end subroutine integ_GL32

end module mod_integrate

