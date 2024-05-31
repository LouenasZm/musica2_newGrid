!==============================================================================
module mod_blasius
!==============================================================================
  !> author: XG
  !> date: February 2018
  !> Similarity solution for compressible laminar boundary layer
!==============================================================================
  use mod_mpi
  use mod_constant
  use mod_eos
  use mod_tranprop
  implicit none
  ! ---------------------------------------------------------------------------
  integer, parameter :: neta=500
  real(wp) :: G_w,G_M1,matto
  real(wp), dimension(:), allocatable :: eta,etai
  real(wp), dimension(:), allocatable :: ffeta,feta,u_bl,v_bl,T_bl,rho_bl
  real(wp), dimension(:), allocatable :: mu_bl,G_bl,M_bl,cp_bl,lambda_bl,Pr_bl
  ! ---------------------------------------------------------------------------

contains

  !==============================================================================
  subroutine compute_blasius
  !==============================================================================
    !> Similarity solution for compressible laminar boundary layer
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer  :: j,it
    real(wp) :: deta                            ! increment for eta
    real(wp) :: C,Ec                            ! physical constant
    real(wp), dimension(:,:), allocatable :: U  ! solution vector
    real(wp), dimension(5) :: k1,k2,k3,k4       ! R-K terms
    real(wp) :: TsTe,rho,mu,cp,lambda,Pra
    ! Parameters for quasi-Newton
    real(wp), parameter :: epsilon = 1.e-12_wp
    real(wp) :: error                           ! for stopping criterion
    real(wp) :: alpha0,alpha1                   ! first shooting parameter
    real(wp) :: beta0,beta1                     ! second shooting parameter
    real(wp) :: g1,g0,h1,h0,aa
    real(wp) :: dgda,dgdb,dhda,dhdb,D
    real(wp) :: da,db,dg,dh,dn2
    ! for post-proc
    real(wp) :: Tj,roj
    
    logical :: is_nan ! obsolete ??
    ! ---------------------------------------------------------------------------

    is_nan = .false.

    if (iproc.eq.0) then
       print *,repeat('=',70)
       print *,'Similarity solution for compressible laminar boundary layer'
       print *,repeat('=',70)
    endif

    s_ref = scalc_tro(T_ref,rho_ref)
    !P0 = p0calc(T_ref,p_ref,rho_ref,Mach)
    h_ref  = ecalc_pro(p_ref, rho_ref, T_ref) + p_ref/rho_ref ! External hentalpy
    Ec     = (Mach*c_ref)**2/h_ref                            ! Eckert parameter

    if (iproc.eq.0) then
       if (eos_type.eq.'pfg') then
          write(*,*) '    T_ref  :', T_ref
          write(*,*) '    p_ref  :', p_ref
          write(*,*) '    rho_ref:', rho_ref
       else
          write(*,*) '    T_ref  :', T_ref/Tc
          write(*,*) '    p_ref  :', p_ref/pc
          write(*,*) '    rho_ref:', rho_ref/roc
       endif
       write(*,*) '    h_ref  :', h_ref
       write(*,*) '    c_ref  :', c_ref
       write(*,*) '    Ec     :', Ec
       write(*,*) '    mu_ref :', mu_ref
       write(*,*) '    Mach   :', Mach
       write(*,*) '    G_ref  :', G_ref
       write(*,*) '    s_ref  :', s_ref
       write(*,*) '    Re /m  :', rho_ref*u_ref/mu_ref
    endif

    ! eta discretization
    ! ==================
    ! dimension

    deta = 9.0_wp/dble(neta-1)      ! step size in eta
    eta(1) = 0
    do j=2,neta
       eta(j) = eta(j-1) + deta
    enddo

    ! allocations
    ! ===========
    allocate( U(5,neta) )
    ! U(1,:)=C*f''(eta)
    ! U(2,:)=f'(eta)=u(eta)/U_e
    ! U(3,:)=f(eta)
    ! U(4,:)=C/Pr*g'(eta)
    ! U(5,:)=g(eta)=h/h_ref  =T/T_ref (<-PFG assumption)

    ! initialisation quasi-Newon
    ! ==========================
    error = 2.0_wp*epsilon

    ! first shooting parameter alpha (to impose u_infty/U_e=1)
    alpha1 = 0.4_wp
    alpha0 = alpha1*0.999_wp            ! <- small perturbation
    ! second shooting parameter beta (to impose T_infty/T_ref=1)
    aa = 1.0_wp + 0.5_wp*sqrt(Pr)*Ec    ! <- first guess with Crocco's relation  <- ?
    beta1 = aa
    beta0 = beta1*0.999_wp              ! <- small perturbation

    ! if (iproc.eq.0) then
    !    write(*,*) 'Shooting parameters:'
    !    write(*,*) 'alpha0:', alpha0,'alpha1:', alpha1
    !    write(*,*) ' beta0:', beta0 ,' beta1:', beta1
    ! endif

    g1 = 0.0_wp
    g0 = 0.0_wp
    h1 = 0.0_wp
    h0 = 0.0_wp
    dgda = (g1-g0)/(alpha1-alpha0)
    dgdb = (g1-g0)/(beta1-beta0)
    dhda = (h1-h0)/(alpha1-alpha0)
    dhdb = (h1-h0)/(beta1-beta0)
    
    ! Shooting method (quasi-Newton)
    ! ==============================
    it = 0
    ! first initialisation vector
    U(:,1) = 0.0_wp
    U(1,1) = alpha1 ! init shooting for f''
    U(5,1) =  beta1 ! init shooting fot g

    blloop: do while ((error>epsilon).and.(it<1000)) ! stopping criteria

       do j = 2, neta

          ! Runge-Kutta integration
          ! -----------------------
          ! T∕Te substep #1
          TsTe = tcalc_ph(p_ref, U(5,j-1)*h_ref, rho_ref, T_ref)/T_ref
          call thermo( TsTe*T_ref, rho, mu, C, cp, lambda, Pra )
          k1 = dU(U(:,j-1),C,Pra,Ec)


          ! T∕Te substep #2
          TsTe = tcalc_ph(p_ref, U(5,j-1)*h_ref, rho_ref, T_ref)/T_ref + deta*k1(5)/2.0_wp
          call thermo( TsTe*T_ref, rho, mu, C, cp, lambda, Pra )
          k2 = dU( U(:,j-1)+deta*k1/2.0_wp, C, Pra, Ec )

          ! T∕Te substep #3
          TsTe = tcalc_ph(p_ref, U(5,j-1)*h_ref, rho_ref, T_ref)/T_ref + deta*k2(5)/2.0_wp
          call thermo( TsTe*T_ref, rho, mu, C, cp, lambda, Pra )
          k3 = dU( U(:,j-1)+deta*k2/2.0_wp, C, Pra, Ec )


          ! T∕Te substep #4
          TsTe = tcalc_ph(p_ref, U(5,j-1)*h_ref, rho_ref, T_ref)/T_ref+deta*k3(5)
          call thermo( TsTe*T_ref, rho, mu, C, cp, lambda, Pra )
          k4 = dU( U(:,j-1)+deta*k3, C, Pra, Ec )

          ! update R-K
          U(:,j) = U(:,j-1) + deta/6.0_wp*(k1 + 2.0_wp*k2 + 2.0_wp*k3 + k4)
       enddo

       ! quasi-Newton
       ! ------------
       g1 = U(2,neta)
       h1 = U(5,neta)

       ! compute increment
       da = alpha1 - alpha0
       db =  beta1 -  beta0
       !if (is_1D)  db = 0.0_wp       ! <- 1D degenerescence
       dg = g1 - g0
       dh = h1 - h0

       ! update Jacobian with Broyden's method
       dn2 = da**2+db**2
       dgda = dgda + da/dn2*(dg - dgda*da - dgdb*db)
       dgdb = dgdb + da/dn2*(dh - dhda*da - dhdb*db)
       dhda = dhda + db/dn2*(dg - dgda*da - dgdb*db)
       dhdb = dhdb + db/dn2*(dh - dhda*da - dhdb*db)

       ! Argh!!! : delete cross-terms (quasi-quasi)
       dgdb = 0.0_wp
       dhda = 0.0_wp

       ! determinant
       D = dgda*dhdb - dgdb*dhda
       !write(*,*) dgda, dgdb, dhda, dhdb, 'D:', D

       ! store old values
       alpha0 = alpha1
       beta0 = beta1
       g0 = g1
       h0 = h1

       ! update shooting paramters
       ! 2D version (explicit inverse of 2x2 Jacobian matrix)
       alpha1= alpha1 - ((g1-1.0_wp)*dhdb-(h1-1.0_wp)*dgdb)/D
       beta1 = beta1 - (-(g1-1.0_wp)*dhda+(h1-1.0_wp)*dgda)/D
       ! compute error
       error = sqrt((g1-1.0_wp)**2+(h1-1.0_wp)**2)

       !if (isnan(error)) then
       !   is_nan=.true.
       !   exit blloop
       !endif

       ! new initialisation vector
       U(:,1) = 0.0_wp
       U(1,1) = alpha1
       U(5,1) = beta1

       ! iteration increment
       it = it+1
    enddo blloop

    if (iproc.eq.0)  write(*,'(4X,A,2X,i0,2X,A,3X,1pE19.12)') &
         'number of iterations:',it, 'error:', error

    ! ---------------------------------------------------------------------------
    if (it.eq.1000) is_nan = .true.

    if (is_nan) return
    ! ---------------------------------------------------------------------------

    ! streamwise velocity: u∕U_e
    u_bl = U(2,:)

    ! Compute primitive variables
    ! ===========================
    do j=1,neta
       ! temperature: T/T_ref
       T_bl(j) = tcalc_ph(p_ref,U(5,j)*h_ref,rho_ref,T_ref)/T_ref
       Tj   = T_bl(j)*T_ref

       ! density: rho/rho_ref
       rho_bl(j) = rocalc_pt(p_ref, Tj, rho_ref)/rho_ref
       roj = rho_bl(j)*rho_ref

       ! Fundamental derivative
       G_bl(j) = gcalc_tro(Tj, roj)

       ! Mach number
       M_bl(j) = u_bl(j)*Mach*c_ref / sqrt(c2calc_tro(Tj, roj) )

       ! Incompressible eta
       ! etai(j) = xint(j,eta(1:j),U(5,1:j))
       etai(j) = xint( j, eta(1:j), 1.0_wp/rho_bl(1:j) )

       ! dynamic viscosity: mu∕mu_ref
       mu_bl(j) = viscosity_law(Tj, roj)/mu_ref

       ! cp
       cp_bl(j) = cpcalc_tro(Tj,roj)

       ! lambda
       lambda_bl(j) = thconductivity(mu_bl(j)*mu_ref,Tj,roj)

       ! Prandtl
       Pr_bl(j) = cp_bl(j)*mu_bl(j)*mu_ref/lambda_bl(j)
    enddo
    etai = sqrt(2.0_wp)*etai

    ! integrant = U(5,:) - U(2,:)
    matto = sqrt(2.0_wp)*xint( neta, eta, 1.0_wp/rho_bl - U(2,:) )

    ! vertical velocity: sqrt(Re_x)*v/U_e
    v_bl = rho_bl/sqrt(2.0_wp)*( eta*U(2,:)-U(3,:) )

    ! f(eta)
    feta = U(3,:)

    ! f''(eta)
    ffeta = U(1,:)/C
    ! ---------------------------------------------------------------------------
    G_w = G_bl(1)
    gmloop: do j=1,neta
       if (M_bl(j).ge.1.0_wp) then
          G_M1 = G_bl(j)
          exit gmloop
       endif
    enddo gmloop
    ! ---------------------------------------------------------------------------
    if (iproc.eq.0) then
       print *,'    Compressible similarity solution'
       print *,'    --------------------------------'
       print *,'    h_neta/h_ref=', U(5,neta)
       print *,'    T_w         =', T_bl(1)*T_ref
       print *,'    mu_w/mu_ref =', mu_bl(1)
       print *,'    G_w         =', G_bl(1)
       print *,'    sqrt(Re_x)deltas/x=',matto
       print *,repeat('=', 70)

       open(13, file='sim_sol.dat', form='formatted', status='replace')
       write(13,'(A)') 'VARIABLES = "`h" "`h_i" "u" "v" "f(`h)" "T" "`r" "`m" "M" "`G" "`etaff" "Pr"'
       do j=1,neta
          write(13,'(11(f15.10,1x),f17.5)') eta(j),  etai(j),  u_bl(j),  v_bl(j) , feta(j)&
               ,T_bl(j),rho_bl(j),mu_bl(j),M_bl(j), G_bl(j), ffeta(j),Pr_bl(j)
       enddo
       close(13)
    endif

    deallocate(U)

  end subroutine compute_blasius

  !==============================================================================
  subroutine thermo(T,rho,mu,C,cp,lambda,Pr)
  !==============================================================================
    !> subroutine thermo for compressible Blasius
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), intent(in)  :: T
    real(wp), intent(out) :: rho, mu, C, cp, lambda, Pr
    ! ---------------------------------------------------------------------------

    rho =rocalc_pt(p_ref,T,rho_ref)

    mu = viscosity_law(T,rho)

    C = (rho*mu)/(rho_ref*mu_ref)

    cp = cpcalc_tro(T,rho)

    lambda = thconductivity(mu,T,rho)

    Pr = cp*mu/lambda

  end subroutine thermo

  !==============================================================================
  function dU(U,C,Pra,Ec)
  !==============================================================================
    !> function differential system for compressible Blasius
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(5) :: U,dU
    real(wp) :: C,Pra,Ec
    ! ---------------------------------------------------------------------------

    dU(1) = -U(1)*U(3)/C
    dU(2) = U(1)/C
    dU(3) = U(2)
    dU(4) = -Ec/C*U(1)**2-U(3)*U(4)*Pra/C
    dU(5) = U(4)*Pra/C
    
  end function dU

  !==============================================================================
  function xint(ni,x,f)
  !==============================================================================
    !> function integration with trapezoidal rule
  !==============================================================================  
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: ni,i
    real(wp) :: x(ni), f(ni), xint
    ! ---------------------------------------------------------------------------

    xint = 0.0_wp
    do i=1,ni-1
       xint = xint + (f(i+1) + f(i))*(x(i+1) - x(i))*0.5_wp
    enddo
    
  end function xint

end module mod_blasius
