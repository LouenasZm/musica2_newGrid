!================================================================================
module mod_wall_model
!================================================================================
  !> author: Aurelien Bienner
  !> date: June 2023
  !> Module to apply wall model -> algebraic or ODE
!================================================================================
  use mod_flow
  use mod_constant
  use mod_grid
  use mod_mpi
  use warnstop
  implicit none
  !------------------------------------------------------------------------------
  integer :: wm_ind ! wall-model static interface
  integer , parameter :: nwm = 50 ! Number of points for ODE
  real(wp), parameter :: res_wm   = 1.e-4_wp
  real(wp) :: dir_jmin,dir_jmax
  !------------------------------------------------------------------------------
  real(wp), dimension(0:nwm+1) :: uwm, twm, rwm, mutwm
  real(wp), dimension(0:nwm+1) :: ycwm
  ! Arrays for storing values beetween iterations (algebraic)
  real(wp), dimension(:,:), allocatable :: utau_low, utau_upp
  ! Arrays for storing values beetween iterations (ODE)
  real(wp), dimension(:,:,:), allocatable :: uwm_low, twm_low, rwm_low &
                                            , uwm_upp, twm_upp, rwm_upp

contains

  function fsuth(t11)
  !----------------------------------------------------------------------------
  ! Compute gas viscosity - Sutherland's Law
  !----------------------------------------------------------------------------
    implicit none
    real(wp) :: fsuth,t11
    real(wp), parameter :: asuth = 1.4579326545176255e-6_wp
    real(wp), parameter :: bsuth = 110.4_wp
    fsuth = asuth*sqrt(t11*t11*t11)/(t11+bsuth)
  end function fsuth

  !==============================================================================
  subroutine init_wm
  !==============================================================================
    !> Initialization of wall-model
  !==============================================================================
    use mod_bc
    use mod_mpi
    use mod_eos
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: j0,j1,i,j,k,ind,l,n
    real(wp) :: u_e,T_e,r_e,p_e,r_wall,hwm
    real(wp) :: dnk,res,loc,delta
    real(wp), parameter :: stch_wm = 2.0_wp
    real(wp), dimension(0:nwm+1) :: sf,sc
    ! ---------------------------------------------------------------------------

    if (iproc.eq.0) print *,"======================================================================"
    if (iproc.eq.0) print *,"Init wall-model"
    if (iproc.eq.0) print *,"======================================================================"

    ! Specified here for the moment
    wm_ind = 4

    ! Initialization of utau_jmin & utau_jmax <~ Only for chan for the moment
    if (is_bc_wall(2,1)) then
       allocate(utau_jmin(1:nx,1:nz))
       utau_jmin(:,:) = utau

       ! Direction grid
       j0 = 1; j1 = j0+wm_ind
       if (y(j1) - y(j0).gt.0) then
          dir_jmin = 1.0_wp
       else
          dir_jmin = -1.0_wp
       endif
     endif
    if (is_bc_wall(2,2)) then
       allocate(utau_jmax(1:nx,1:nz))
       utau_jmax(:,:) = utau

       ! Direction grid
       j0 = ny; j1 = j0-wm_ind
       if (y(j1) - y(j0).gt.0) then
          dir_jmax = 1.0_wp
       else
          dir_jmax = -1.0_wp
       endif
    endif



    ! Initialization for ODE
    ! ----------------------
    if (wm_model_type.eq."ODE") then
       if (iproc.eq.0) print *,"   wm_model_type: ODE"

       if ((is_bc_wall(2,1)).or.(is_bc_wall(2,2))) then
          ! Create phantom grid for wall model
          ! ----------------------------------
          hwm = (yg(1+wm_ind) - yg(1))
          ! hwm = (yg(1+wm_ind) - yg(1))*dir_jmin
          do k = nwm+1,0,-1
             dnk   = dble(nwm+1-k)/dble(nwm+1)
             sf(k) = 1.0_wp - tanh(stch_wm*dnk)/tanh(stch_wm)
          enddo
          do k = 1,nwm+1
             sc(k) = (sf(k)+sf(k-1))*0.5_wp
          enddo
          sc(0) = -sc(1)

          ! Adjust nodes at top so that the top ghost cell center
          ! matches the LES exchange location cell center
          ! ---------------------------------------------
          sc(nwm-1) = sc(nwm-2) + (sf(nwm+1)-sc(nwm-2))*ONE_THIRD
          sc(nwm  ) = sc(nwm-1) + (sf(nwm+1)-sc(nwm-2))*ONE_THIRD
          sc(nwm+1) = sf(nwm+1)

          ycwm(:) = hwm*sc(:) ! cell centers
          ! do k = 0,nwm
          !    yfwm(k) = 0.5_wp*(ycwm(k)+ycwm(k+1)) ! face centers 'j+1/2'
          ! enddo

          do j=2,wm_ind
             res = 10000_wp
             do n = 1,nwm
                loc = 0.5_wp*(ycwm(n)+ycwm(n-1))
                delta = abs(abs(yg(j)-yg(1))-loc)
                if (delta <  res ) then
                   res = delta
                   ind = n
                endif
             enddo
             ! indici(j-1) = ind
          enddo
       endif

       if (is_bc_wall(2,1)) then
          ! Allocate arrays for storing values
          ! ----------------------------------
          allocate(uwm_low(0:nwm+1,1:nx,1:nz) &
                  ,rwm_low(0:nwm+1,1:nx,1:nz) &
                  ,twm_low(0:nwm+1,1:nx,1:nz))

          ! Initialize ODE arrays
          ! ---------------------
          j1 = 1+wm_ind

          do k=1,nz
             do i=1,nx
                u_e = uu(i,j1,k)
                ! u_e = sqrt(uu(i,j1,k)**2 + ww(i,j1,k)**2)
                ! u_e = sqrt( (vel(i,j1,k,1) + Uscale)**2 + vel(i,j1,k,3)**2 )
                T_e = Tmp(i,j1,k)
                p_e = prs(i,j1,k)
                r_e = rho(i,j1,k)
                uwm_low(1:nwm,i,k) = u_e*sc(1:nwm)
                twm_low(1:nwm,i,k) = T_e*sc(1:nwm) + T_wall*(1.0_wp - sc(1:nwm))
                do l=1,nwm
                   rwm_low(l,i,k) = rocalc_pt(uwm_low(l,i,k),twm_low(l,i,k),r_e)
                enddo
                !---------BCs for Isothermal wall-------------------
                r_wall = rwm_low(1,i,k)*twm_low(1,i,k)/T_wall
                uwm_low(    0,i,k) = -uwm_low(1,i,k)           ! no-slip wall
                rwm_low(    0,i,k) = 2.0_wp*r_wall - rwm_low(1,i,k)
                twm_low(    0,i,k) = 2.0_wp*T_wall - twm_low(1,i,k)
                uwm_low(nwm+1,i,k) = u_e                        ! dirichlet
                rwm_low(nwm+1,i,k) = r_e
                twm_low(nwm+1,i,k) = T_e
             enddo
          enddo
       endif

       if (is_bc_wall(2,2)) then
          ! Allocate arrays for storing values
          ! ----------------------------------
          allocate(twm_upp(0:nwm+1,1:nx,1:nz) &
                  ,uwm_upp(0:nwm+1,1:nx,1:nz) &
                  ,rwm_upp(0:nwm+1,1:nx,1:nz) )

          ! Initialize ODE arrays
          ! ---------------------
          j1 = ny-wm_ind

          do k=1,nz
             do i=1,nx
                u_e = uu(i,j1,k)
                T_e = Tmp(i,j1,k)
                p_e = prs(i,j1,k)
                r_e = rho(i,j1,k)
                uwm_upp(1:nwm,i,k) = u_e*sc(1:nwm)
                twm_upp(1:nwm,i,k) = T_e*sc(1:nwm) + T_wall*(1.0_wp - sc(1:nwm))
                do l=1,nwm
                   rwm_upp(l,i,k) = rocalc_pt(uwm_upp(l,i,k),twm_upp(l,i,k),r_e)
                enddo
                !---------BCs for Isothermal wall-------------------
                r_wall = rwm_upp(1,i,k)*twm_upp(1,i,k)/T_wall
                uwm_upp(    0,i,k) = -uwm_upp(1,i,k)           ! no-slip wall
                rwm_upp(    0,i,k) = 2.0_wp*r_wall - rwm_upp(1,i,k)
                twm_upp(    0,i,k) = 2.0_wp*T_wall - twm_upp(1,i,k)
                uwm_upp(nwm+1,i,k) = u_e                        ! dirichlet
                rwm_upp(nwm+1,i,k) = r_e
                twm_upp(nwm+1,i,k) = T_e
             enddo
          enddo
       endif
    else
       if (iproc.eq.0) print *,"   wm_model_type: ALG"
    endif


  end subroutine init_wm

  !=====================================================================================================================================================================
  !                                         Algebraic wall model
  !=====================================================================================================================================================================

  !==============================================================================
  subroutine bc_wm_alg_jmin
  !==============================================================================
    !> Use algebraic model
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,k,j0,j1
    real(wp) :: tau_w,u_e,T_e,r_e,T_w,r_w,mu_w,hwm,utau_wm, q_w
    ! ---------------------------------------------------------------------------
    ! ! Test stochastic forcing
    ! real(wp) :: tau11,tau22,tau33,tau12,tau13,tau23,trace,mu
    ! real(wp) :: coeff_kappa,var,yp_
    ! real(wp) :: dupx,dupy,dupz,dvpx,dvpy,dvpz,dwpx,dwpy,dwpz
    ! real(wp), parameter :: sp_kappa = 0.41_wp
    ! ! random number generator
    ! real(wp) :: randomnb
    ! integer :: nseed
    ! integer, dimension(:), allocatable :: initseed

    j0 = 1; j1 = j0+wm_ind
    hwm = abs(y(j1) - y(j0))

    do k=1,nz
       do i=1,nx
          ! Give first guess for Utau
          utau_wm = utau_jmin(i,k)
          ! Fill edge values
          u_e = (uu(i,j1,k)**2 + ww(i,j1,k)**2)**0.5_wp
          ! u_e = uu(i,j1,k)
          T_e = Tmp(i,j1,k)
          r_e = rho(i,j1,k)
          ! Fill wall values
          T_w = Tmp(i,j0,k)
          r_w = rho(i,j0,k)
          mu_w = visc(i,j0,k)

          ! Only implemented with spalding for the moment
          call wm_alg_spalding(u_e,T_e,r_e,T_w,r_w,mu_w,hwm,utau_wm,q_w)

          ! Store new value
          utau_jmin(i,k) = utau_wm

          ! Calculation of tau_w
          tau_w = rho(i,j0,k)*utau_jmin(i,k)**2

          ! Impose wall-model values
          Frhou(i,j0,k) = 0.0_wp; Frhow(i,j0,k) = 0.0_wp
          Grhov(i,j0,k) = 0.0_wp; Grhow(i,j0,k) = 0.0_wp
          Hrhou(i,j0,k) = 0.0_wp; Hrhov(i,j0,k) = 0.0_wp; Hrhow(i,j0,k) = 0.0_wp; Hrhoe(i,j0,k) = 0.0_wp
          Frhov(i,j0,k) = -dir_jmin*tau_w
          Grhou(i,j0,k) = -dir_jmin*tau_w
          Frhoe(i,j0,k) = -dir_jmin*q_w
          Grhoe(i,j0,k) = -dir_jmin*q_w
       enddo
    enddo

 !    ! Stochastic forcing on j=2
 !    ! ==================
 !    ! j=2
 !    j=j1
 !    coeff_kappa = 1.0_wp/(15*sp_kappa)**0.5
 !    var = 0.0_wp

 !    ! Initialization of random number generator
 !    !------------------------------------------
 !    call random_seed(size=nseed)
 !    allocate(initseed(nseed))
 !    do k=1,iproc+1
 !       call random_number(randomnb)
 !    enddo
 !    do k=1,nseed
 !       call random_number(randomnb)
 !       initseed(k) = int(randomnb*10**8) * (-1)**k * iproc
 !    enddo
 !    call random_seed(put=initseed)
 !    deallocate(initseed)

 !    ! if ((iproc.eq.0).or.(iproc.eq.1)) then
 !    do k=1,nz
 !       do i=1,nx
 !          yp_ = abs(y(j)-y(j0))*utau_jmin(i,k)*rho(i,j0,k)/visc(i,j0,k)
 !          tau_w = rho(i,j0,k)*utau_jmin(i,k)**2
 !          tau_w = tau_w*coeff_kappa/(yp_)**0.5/visc(i,j0,k)
 !          call random_number(var)
 !          dupx = tau_w * (2)**0.5 * ERF(2*var-1) ! du'/dx
 !          ! print *,"dupx",iproc,dupx
 !          ! print *,"dux",dux(i,j,k)
 !          ! print *,"duy",iproc,duy(i,j,k)
 !          ! print *,"duz",duz(i,j,k)
 !          ! print *,"dvx",dvx(i,j,k)
 !          ! print *,"dvy",dvy(i,j,k)
 !          ! print *,"dvz",dvz(i,j,k)
 !          ! print *,"dwx",dwx(i,j,k)
 !          ! print *,"dwy",dwy(i,j,k)
 !          ! print *,"dwz",dwz(i,j,k)
 !          call random_number(var)
 !          dupy = tau_w * 2 * ERF(2*var-1) ! du'/dy
 !          call random_number(var)
 !          dupz = tau_w * 2 * ERF(2*var-1) ! du'/dz
 !          call random_number(var)
 !          dvpx = tau_w * 2 * ERF(2*var-1) ! dv'/dx
 !          call random_number(var)
 !          dvpy = tau_w * (2)**0.5 * ERF(2*var-1) ! dv'/dy
 !          call random_number(var)
 !          dvpz = tau_w * 2 * ERF(2*var-1) ! dv'/dz
 !          call random_number(var)
 !          dwpx = tau_w * 2 * ERF(2*var-1) ! dw'/dx
 !          call random_number(var)
 !          dwpy = tau_w * 2 * ERF(2*var-1) ! dw'/dy
 !          call random_number(var)
 !          dwpz = tau_w * (2)**0.5 * ERF(2*var-1) ! dw'/dz
 !          ! print *,"dup",iproc,dupx,dupy,dupz
 !          ! print *,"dvp",iproc,dvpx,dvpy,dvpz
 !          ! print *,"dwp",iproc,dwpx,dwpy,dwpz

 !          ! print *,"Frhou bef",Frhou(i,j,k)
 !          ! print *,"Frhov bef",Frhov(i,j,k)


 !          ! compute S_ij
 !          tau11 = dupx
 !          tau22 = dvpy
 !          tau33 = dwpz
 !          tau12 = 0.5_wp*(dupy + dvpx)
 !          tau13 = 0.5_wp*(dupz + dwpx)
 !          tau23 = 0.5_wp*(dvpz + dwpy)
 !          trace = ONE_THIRD*(tau11+tau22+tau33)

 !          ! compute -tau_ij
 !          mu =-2.0_wp*visc(i,j,k)
 !          tau11=mu*(tau11-trace)
 !          tau22=mu*(tau22-trace)
 !          tau33=mu*(tau33-trace)
 !          tau12=mu*tau12
 !          tau13=mu*tau13
 !          tau23=mu*tau23

 !          ! viscous fluxes along x
 !          Frhou(i,j,k) = Frhou(i,j,k) + tau11
 !          Frhov(i,j,k) = Frhov(i,j,k) + tau12
 !          Frhow(i,j,k) = Frhow(i,j,k) + tau13
 !          Frhoe(i,j,k) = Frhoe(i,j,k) + (uu(i,j,k)*tau11 + vv(i,j,k)*tau12 + ww(i,j,k)*tau13)

 !          ! viscous fluxes along y
 !          Grhou(i,j,k) = Grhou(i,j,k) + tau12
 !          Grhov(i,j,k) = Grhov(i,j,k) + tau22
 !          Grhow(i,j,k) = Grhow(i,j,k) + tau23
 !          Grhoe(i,j,k) = Grhoe(i,j,k) + (uu(i,j,k)*tau12 + vv(i,j,k)*tau22 + ww(i,j,k)*tau23)

 !          ! viscous fluxes along z
 !          Hrhou(i,j,k) = Hrhou(i,j,k) + tau13
 !          Hrhov(i,j,k) = Hrhov(i,j,k) + tau23
 !          Hrhow(i,j,k) = Hrhow(i,j,k) + tau33
 !          Hrhoe(i,j,k) = Hrhoe(i,j,k) + (uu(i,j,k)*tau13 + vv(i,j,k)*tau23 + ww(i,j,k)*tau33)

 !          ! print *,"Frhou af",Frhou(i,j,k)
 !          ! print *,"Frhov af",Frhov(i,j,k)

 !          ! call mpistop('',0)
 !       enddo
 !    enddo
 ! ! endif
 !          ! call mpistop('',0)


  end subroutine bc_wm_alg_jmin

  !==============================================================================
  subroutine bc_wm_alg_jmax
  !==============================================================================
    !> Use algebraic model
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer :: i,k,j0,j1
    real(wp) :: tau_w,u_e,T_e,r_e,T_w,r_w,mu_w,hwm,utau_wm, q_w
    ! ---------------------------------------------------------------------------
    ! ! Test stochastic forcing
    ! real(wp) :: tau11,tau22,tau33,tau12,tau13,tau23,trace,mu
    ! real(wp) :: coeff_kappa,var,yp_
    ! real(wp) :: dupx,dupy,dupz,dvpx,dvpy,dvpz,dwpx,dwpy,dwpz
    ! real(wp), parameter :: sp_kappa = 0.41_wp
    ! ! random number generator
    ! real(wp) :: randomnb
    ! integer :: nseed
    ! integer, dimension(:), allocatable :: initseed

    j0 = ny; j1 = j0-wm_ind
    hwm = abs(y(j1) - y(j0))

    do k=1,nz
       do i=1,nx
          ! Give first guess for Utau
          utau_wm = utau_jmax(i,k)
          ! Fill edge values
          u_e = (uu(i,j1,k)**2 + ww(i,j1,k)**2)**0.5_wp
          ! u_e = uu(i,j1,k)
          T_e = Tmp(i,j1,k)
          r_e = rho(i,j1,k)
          ! Fill wall values
          T_w = Tmp(i,j0,k)
          r_w = rho(i,j0,k)
          mu_w = visc(i,j0,k)

          ! Only implemented with spalding for the moment
          call wm_alg_spalding(u_e,T_e,r_e,T_w,r_w,mu_w,hwm,utau_wm,q_w)

          ! Store new value
          utau_jmax(i,k) = utau_wm

          ! Calculation of tau_w
          tau_w = rho(i,j0,k)*utau_jmax(i,k)**2

          ! Impose wall-model values
          Frhou(i,j0,k) = 0.0_wp; Frhow(i,j0,k) = 0.0_wp
          Grhov(i,j0,k) = 0.0_wp; Grhow(i,j0,k) = 0.0_wp
          Hrhou(i,j0,k) = 0.0_wp; Hrhov(i,j0,k) = 0.0_wp; Hrhow(i,j0,k) = 0.0_wp; Hrhoe(i,j0,k) = 0.0_wp
          Frhov(i,j0,k) = -dir_jmax*tau_w
          Grhou(i,j0,k) = -dir_jmax*tau_w
          Frhoe(i,j0,k) = -dir_jmax*q_w
          Grhoe(i,j0,k) = -dir_jmax*q_w
       enddo
    enddo

 !    ! Stochastic forcing on j=ny-1
 !    ! ==================
 !    ! j=ny-1
 !    j=j1
 !    coeff_kappa = 1.0_wp/(15*sp_kappa)**0.5
 !    var = 0.0_wp

 !    ! Initialization of random number generator
 !    !------------------------------------------
 !    call random_seed(size=nseed)
 !    allocate(initseed(nseed))
 !    do k=1,iproc+1
 !       call random_number(randomnb)
 !    enddo
 !    do k=1,nseed
 !       call random_number(randomnb)
 !       initseed(k) = int(randomnb*10**8) * (-1)**k * iproc
 !    enddo
 !    call random_seed(put=initseed)
 !    deallocate(initseed)

 !    ! if ((iproc.eq.0).or.(iproc.eq.1)) then
 !    do k=1,nz
 !       do i=1,nx
 !          yp_ = abs(y(j)-y(j0))*utau_jmax(i,k)*rho(i,j0,k)/visc(i,j0,k)
 !          tau_w = rho(i,j0,k)*utau_jmax(i,k)**2
 !          tau_w = tau_w*coeff_kappa/(yp_)**0.5/visc(i,j0,k)
 !          call random_number(var)
 !          dupx = tau_w * (2)**0.5 * ERF(2*var-1) ! du'/dx
 !          call random_number(var)
 !          dupy = tau_w * 2 * ERF(2*var-1) ! du'/dy
 !          call random_number(var)
 !          dupz = tau_w * 2 * ERF(2*var-1) ! du'/dz
 !          call random_number(var)
 !          dvpx = tau_w * 2 * ERF(2*var-1) ! dv'/dx
 !          call random_number(var)
 !          dvpy = tau_w * (2)**0.5 * ERF(2*var-1) ! dv'/dy
 !          call random_number(var)
 !          dvpz = tau_w * 2 * ERF(2*var-1) ! dv'/dz
 !          call random_number(var)
 !          dwpx = tau_w * 2 * ERF(2*var-1) ! dw'/dx
 !          call random_number(var)
 !          dwpy = tau_w * 2 * ERF(2*var-1) ! dw'/dy
 !          call random_number(var)
 !          dwpz = tau_w * (2)**0.5 * ERF(2*var-1) ! dw'/dz
 !          ! print *,"dup",iproc,dupx,dupy,dupz
 !          ! print *,"dvp",iproc,dvpx,dvpy,dvpz
 !          ! print *,"dwp",iproc,dwpx,dwpy,dwpz

 !          ! print *,"Frhou bef",Frhou(i,j,k)
 !          ! print *,"Frhov bef",Frhov(i,j,k)


 !          ! compute S_ij
 !          tau11 = dupx
 !          tau22 = dvpy
 !          tau33 = dwpz
 !          tau12 = 0.5_wp*(dupy + dvpx)
 !          tau13 = 0.5_wp*(dupz + dwpx)
 !          tau23 = 0.5_wp*(dvpz + dwpy)
 !          trace = ONE_THIRD*(tau11+tau22+tau33)

 !          ! compute -tau_ij
 !          mu =-2.0_wp*visc(i,j,k)
 !          tau11=mu*(tau11-trace)
 !          tau22=mu*(tau22-trace)
 !          tau33=mu*(tau33-trace)
 !          tau12=mu*tau12
 !          tau13=mu*tau13
 !          tau23=mu*tau23

 !          ! viscous fluxes along x
 !          Frhou(i,j,k) = Frhou(i,j,k) + tau11
 !          Frhov(i,j,k) = Frhov(i,j,k) + tau12
 !          Frhow(i,j,k) = Frhow(i,j,k) + tau13
 !          Frhoe(i,j,k) = Frhoe(i,j,k) + (uu(i,j,k)*tau11 + vv(i,j,k)*tau12 + ww(i,j,k)*tau13)

 !          ! viscous fluxes along y
 !          Grhou(i,j,k) = Grhou(i,j,k) + tau12
 !          Grhov(i,j,k) = Grhov(i,j,k) + tau22
 !          Grhow(i,j,k) = Grhow(i,j,k) + tau23
 !          Grhoe(i,j,k) = Grhoe(i,j,k) + (uu(i,j,k)*tau12 + vv(i,j,k)*tau22 + ww(i,j,k)*tau23)

 !          ! viscous fluxes along z
 !          Hrhou(i,j,k) = Hrhou(i,j,k) + tau13
 !          Hrhov(i,j,k) = Hrhov(i,j,k) + tau23
 !          Hrhow(i,j,k) = Hrhow(i,j,k) + tau33
 !          Hrhoe(i,j,k) = Hrhoe(i,j,k) + (uu(i,j,k)*tau13 + vv(i,j,k)*tau23 + ww(i,j,k)*tau33)

 !          ! print *,"Frhou af",Frhou(i,j,k)
 !          ! print *,"Frhov af",Frhov(i,j,k)

 !          ! call mpistop('',0)
 !       enddo
 !    enddo
 ! ! endif

  end subroutine bc_wm_alg_jmax

  !===============================================================================
  subroutine wm_alg_spalding(u_e,T_e,r_e,T_w,r_w,mu_w,hwm,utau_wm,q_w)
  !===============================================================================
    !> Apply Spalding's and Kader's laws
  !===============================================================================
    use mod_eos
    use mod_tranprop
    implicit none
    ! ----------------------------------------------------------------------------
    ! Input/Output arguments
    real(wp), intent(in) :: u_e,T_e,r_e,T_w,r_w,mu_w,hwm
    real(wp), intent(inout) :: utau_wm, q_w
    ! ----------------------------------------------------------------------------
    ! Local variables
    integer  :: n
    real(wp), parameter :: sp_kappa = 0.4_wp, sp_B = 5.5_wp
    real(wp) :: exp_mkappaB,kup,uplus,fun,der,u_test
    real(wp) :: gamma,beta,yplus,Tplus,cpav,inv_nu,Prdtl
    ! ----------------------------------------------------------------------------

    inv_nu = r_w/mu_w
    cpav = cpcalc_tro(T_w,r_w)
    ! cpav = comp(1)%cp0ig/mw(1)
    Prdtl = mu_w*cpav/thconductivity(mu_w,T_w,r_w)
    exp_mkappaB = exp(-sp_kappa*sp_B)

    ! Newton's algorithm
    ! ------------------
    newtonloop: do n = 1,100
       uplus = u_e/utau_wm
       kup = sp_kappa*uplus
       fun = uplus + exp_mkappaB*( exp(kup) - 1.0_wp - kup - 0.5_wp*kup**2 &
                                  - kup**3/6.0_wp ) - hwm*utau_wm*inv_nu
       der = -hwm*inv_nu -uplus/utau_wm - kup/utau_wm*exp_mkappaB &
                                 *( exp(kup) - 1.0_wp - kup - 0.5_wp*kup**2 )
       u_test = utau_wm - fun/der

       if (abs(u_test-utau_wm)/utau_wm.le.1.e-4_wp) then
          utau_wm = u_test
          exit newtonloop
       endif

       if (n.ge.100 ) then
          write(*,*) n, abs(u_test-utau_wm)/utau_wm, u_test
          call mpistop('Problem with convergence of newton loop in wm_alg_spalding...',0)
       endif
       utau_wm = u_test
    enddo newtonloop

    ! heat flux
    yplus = hwm*utau_wm*inv_nu
    ! Kader law
    gamma = 1e-2_wp*(Prdtl*yplus)**4/(1.0_wp+5.0_wp*Prdtl**3*yplus)
    beta  = (3.85_wp*Prdtl**ONE_THIRD - 1.3_wp)**2 + 2.12_wp*log(Prdtl)
    Tplus = Prdtl*yplus*exp(-gamma) + (2.12_wp*log(1.0_wp+yplus) + beta)*exp(-1.0_wp/gamma)

    q_w = (T_e - T_w)*r_w*cpav*utau_wm/Tplus

  end subroutine wm_alg_spalding

  !=================================================================================================================================================================
  !                                             Ordinary Differential Equation model
  !=================================================================================================================================================================

  !==============================================================================
  subroutine bc_wm_ODE_jmin
  !==============================================================================
    !> Wal-model ODE for jmin wall
  !==============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer  :: i,j0,j1,k
    real(wp) :: u_e,T_e,r_e,p_e,T_w,r_w,mu_w,tau_w,q_w
    ! ------------------------------------------------------------------------

    j0 = 1; j1 = j0+wm_ind

    do k=1,nz
       do i=1,nx
          ! Fill edge values
          u_e = uu(i,j1,k)
          T_e = Tmp(i,j1,k)
          r_e = rho(i,j1,k)
          p_e = prs(i,j1,k)
          ! Fill wall values
          T_w = Tmp(i,j0,k)
          r_w = rho(i,j0,k)
          mu_w = visc(i,j0,k)

          uwm(0:nwm+1) = uwm_low(0:nwm+1,i,k)
          twm(0:nwm+1) = twm_low(0:nwm+1,i,k)
          rwm(0:nwm+1) = rwm_low(0:nwm+1,i,k)

          call wm_ODE(u_e,p_e,T_e,r_e,T_w,r_w,mu_w,tau_w,q_w)

          ! Save new values
          uwm_low(0:nwm+1,i,k) = uwm(0:nwm+1)
          twm_low(0:nwm+1,i,k) = twm(0:nwm+1)
          rwm_low(0:nwm+1,i,k) = rwm(0:nwm+1)

          ! Store new value
          utau_jmin(i,k) = (tau_w/rho(i,j0,k))**0.5

          ! Impose wall-model values
          Frhou(i,j0,k) = 0.0_wp; Frhow(i,j0,k) = 0.0_wp
          Grhov(i,j0,k) = 0.0_wp; Grhow(i,j0,k) = 0.0_wp
          Hrhou(i,j0,k) = 0.0_wp; Hrhov(i,j0,k) = 0.0_wp; Hrhow(i,j0,k) = 0.0_wp; Hrhoe(i,j0,k) = 0.0_wp
          Frhov(i,j0,k) = -dir_jmin*tau_w
          Grhou(i,j0,k) = -dir_jmin*tau_w
          Frhoe(i,j0,k) = -dir_jmin*q_w
          Grhoe(i,j0,k) = -dir_jmin*q_w
       enddo
    enddo

  end subroutine bc_wm_ODE_jmin

  !==============================================================================
  subroutine bc_wm_ODE_jmax
  !==============================================================================
    !> Wal-model ODE for jmin wall
  !==============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer  :: i,j0,j1,k
    real(wp) :: u_e,T_e,r_e,p_e,T_w,r_w,mu_w,tau_w,q_w
    ! ------------------------------------------------------------------------

    j0 = ny; j1 = j0-wm_ind

    do k=1,nz
       do i=1,nx
          ! Fill edge values
          u_e = uu(i,j1,k)
          T_e = Tmp(i,j1,k)
          r_e = rho(i,j1,k)
          p_e = prs(i,j1,k)
          ! Fill wall values
          T_w = Tmp(i,j0,k)
          r_w = rho(i,j0,k)
          mu_w = visc(i,j0,k)

          uwm(0:nwm+1) = uwm_upp(0:nwm+1,i,k)
          twm(0:nwm+1) = twm_upp(0:nwm+1,i,k)
          rwm(0:nwm+1) = rwm_upp(0:nwm+1,i,k)

          call wm_ODE(u_e,p_e,T_e,r_e,T_w,r_w,mu_w,tau_w,q_w)

          ! Save new values
          uwm_upp(0:nwm+1,i,k) = uwm(0:nwm+1)
          twm_upp(0:nwm+1,i,k) = twm(0:nwm+1)
          rwm_upp(0:nwm+1,i,k) = rwm(0:nwm+1)

          ! Store new value
          utau_jmax(i,k) = (tau_w/rho(i,j0,k))**0.5

          ! Impose wall-model values
          Frhou(i,j0,k) = 0.0_wp; Frhow(i,j0,k) = 0.0_wp
          Grhov(i,j0,k) = 0.0_wp; Grhow(i,j0,k) = 0.0_wp
          Hrhou(i,j0,k) = 0.0_wp; Hrhov(i,j0,k) = 0.0_wp; Hrhow(i,j0,k) = 0.0_wp; Hrhoe(i,j0,k) = 0.0_wp
          Frhov(i,j0,k) = -dir_jmax*tau_w
          Grhou(i,j0,k) = -dir_jmax*tau_w
          Frhoe(i,j0,k) = -dir_jmax*q_w
          Grhoe(i,j0,k) = -dir_jmax*q_w
       enddo
    enddo

  end subroutine bc_wm_ODE_jmax


  !==============================================================================
  subroutine wm_ODE(u_e,p_e,T_e,r_e,T_w,r_w,mu_w,tau_w,q_w)
  !==============================================================================
    !> Solve the ODE
  !==============================================================================
    use mod_eos
    use mod_tranprop
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), intent(in)  :: u_e,p_e,T_e,r_e,T_w,r_w,mu_w
    real(wp), intent(out) :: tau_w,q_w
    ! ----------------------------------------------------------------------------
    integer  :: kwm,iter,l
    real(wp) :: rwall
    real(wp) :: tau_wm,rtau_w,vdamp,ystar
    real(wp) :: urms,trms
    real(wp) :: uwmf,twmf,rwmf,muwmf,mutwmf,dywm,duwm
    real(wp) :: Prdtl_inv, cpav
    ! Boundary conditions and free-stream values
    ! integer, parameter :: iwall   = 2
    real(wp), parameter :: prti    = 1.0d0/0.9d0
    ! Constants mixing length model
    real(wp), parameter :: kt_wm   = 0.41_wp
    real(wp), parameter :: apl1_wm = 17.0_wp
    !-------------------------------------------------------------------------
    ! Arrays for TDMA algorithm
    real(wp), dimension(nwm) :: a_low, b_dia, c_upp, d_rhs, xuwm, xtwm

    cpav = cpcalc_tro(T_w,r_w)
    Prdtl_inv = thconductivity(mu_w,T_w,r_w)/(mu_w*cpav)

    do iter=1,100
       ! Solve momentum equation
       ! -----------------------
       tau_wm = abs(mu_w*uwm(1)/ycwm(1))

       ! select case (iturb_wm)
       !    case(1) ! Mixing length model, JK
             do kwm = 1,nwm+1
                muwmf = fsuth(twm(kwm)) ! TO BE GENERALIZED
                rtau_w = sqrt(rwm(kwm)*tau_wm)
                ystar = ycwm(kwm)*rtau_w/muwmf       ! semi-local scaling
                vdamp = (1.0_wp-exp(-ystar/apl1_wm))**2
                mutwm(kwm) = kt_wm*ycwm(kwm)*rtau_w*vdamp
             enddo
             mutwm(0) = -mutwm(1)
       !    case default
       !       mutwm = 0
       ! end select
       !- Form matrix system
       d_rhs = 0.0_wp
       do kwm = 1,nwm
          twmf = 0.5_wp*( twm(kwm) + twm(kwm-1) )
          rwmf = 0.5_wp*( rwm(kwm) + rwm(kwm-1) )
          dywm   = ycwm(kwm) - ycwm(kwm-1)
          muwmf  = fsuth(twmf)
          mutwmf = 0.5_wp*( mutwm(kwm) + mutwm(kwm-1) )
          !- Lower
          a_low(kwm) = (muwmf + mutwmf)/dywm
          twmf  = 0.5_wp*( twm(kwm+1) + twm(kwm) )
          rwmf  = 0.5_wp*( rwm(kwm+1) + rwm(kwm) )
          dywm   = ycwm(kwm+1) - ycwm(kwm)
          muwmf  = fsuth(twmf)
          mutwmf = 0.5_wp*( mutwm(kwm+1) + mutwm(kwm) )
          !- Upper
          c_upp(kwm) = (muwmf + mutwmf)/dywm
          !- Diagonal
          b_dia(kwm) = -(a_low(kwm) + c_upp(kwm))
       enddo ! nwm

       !- Implicit BCs: no-slip wall
       b_dia(1) = b_dia(1) - a_low(1)
       d_rhs(1) = d_rhs(1) + 0.0_wp
       !- Implicit BCs: dirichlet top
       d_rhs(nwm) = d_rhs(nwm) - u_e*c_upp(nwm)
       !- Solve tridiagonal system for U
       call tdma(nwm,a_low,b_dia,c_upp,d_rhs,xuwm ) ! To be changed
       urms = sum((uwm(1:nwm) - xuwm(1:nwm))**2)
       urms = sqrt(urms/dble(nwm))/u_e
       ! Update - U
       uwm(1:nwm) = xuwm(1:nwm)
       ! Explicit BCs - U
       uwm(0)     = -uwm(1)     ! no-slip wall
       uwm(nwm+1) = u_e         ! dirichlet

       !----------------------------------------------------------------------
       ! Solve energy equation

       tau_wm = abs(mu_w*uwm(1)/ycwm(1))
       ! select case (iturb_wm)
          ! case(1) ! Mixing length model, JK
             do kwm = 1,nwm+1
                muwmf = fsuth(twm(kwm))
                rtau_w = sqrt(rwm(kwm)*tau_wm)
                ystar = ycwm(kwm)*rtau_w/muwmf       ! semi-local scaling
                vdamp = (1.0_wp-exp(-ystar/apl1_wm))**2
                mutwm(kwm) = kt_wm*ycwm(kwm)*rtau_w*vdamp
             enddo
             mutwm(0) = -mutwm(1)
       !    case default
       !       mutwm = 0
       ! end select
       !- Form matrix system
       d_rhs = 0.0_wp
       do kwm = 1,nwm
          uwmf = 0.5_wp*( uwm(kwm) + uwm(kwm-1) )
          twmf = 0.5_wp*( twm(kwm) + twm(kwm-1) )
          rwmf = 0.5_wp*( rwm(kwm) + rwm(kwm-1) )
          duwm   = uwm (kwm) - uwm (kwm-1)
          dywm   = ycwm(kwm) - ycwm(kwm-1)
          muwmf  = fsuth(twmf)
          mutwmf = 0.5_wp*( mutwm(kwm) + mutwm(kwm-1) )
          !- Lower
          a_low(kwm) = cpav*(muwmf*Prdtl_inv + mutwmf*prti)/dywm
          !- RHS
          d_rhs(kwm) = d_rhs(kwm) + (muwmf+mutwmf)*uwmf*duwm/dywm
          uwmf = 0.5_wp*( uwm(kwm+1) + uwm(kwm) )
          twmf = 0.5_wp*( twm(kwm+1) + twm(kwm) )
          rwmf = 0.5_wp*( rwm(kwm+1) + rwm(kwm) )
          duwm   = uwm (kwm+1) - uwm (kwm)
          dywm   = ycwm(kwm+1) - ycwm(kwm)
          muwmf  = fsuth(twmf)
          mutwmf = 0.5_wp*( mutwm(kwm+1) + mutwm(kwm) )
          !- Upper
          c_upp(kwm) = cpav*(muwmf*Prdtl_inv + mutwmf*prti)/dywm
          !- Diagonal
          b_dia(kwm) = -(a_low(kwm) + c_upp(kwm))
          !- RHS: viscous work
          d_rhs(kwm) = d_rhs(kwm) - (muwmf+mutwmf)*uwmf*duwm/dywm
       enddo ! nwm
       !- Implicit BCs: wall
       ! if (iwall==1) then
       !     ! Adiabatic
       !     b_dia(1) = b_dia(1) + a_low(1)
       !     d_rhs(1) = d_rhs(1) + 0.0_wp
       ! elseif (iwall==2 .or. iwall==0) then
          ! Isothermal
          b_dia(1) = b_dia(1) - a_low(1)
          d_rhs(1) = d_rhs(1) - 2.0_wp*T_w*a_low(1)
       ! endif
       !- Implicit BCs: dirichlet top
       d_rhs(nwm) = d_rhs(nwm) - T_e*c_upp(nwm)
       !- Solve tridiagonal system for T
       call tdma( nwm, a_low, b_dia, c_upp, d_rhs, xtwm )
       ! Get Trms
       trms = sum((twm(1:nwm) - xtwm(1:nwm))**2)
       trms = sqrt(trms/dble(nwm))/T_e
       ! Update - RHO & T
       twm(1:nwm) = xtwm(1:nwm)           ! T
       do l=1,nwm
          rwm(kwm) = rocalc_pt(uwm(kwm),twm(kwm),r_e)
       enddo
       ! rwm(1:nwm) = p_e/(Rair*twm(1:nwm)) ! Rho
       ! Explicit BCs - RHO & T
       ! if (iwall==1) then      ! adiabatic
          ! rwall  = rwm(1)
          ! T_w  = twm(1)
          ! mu_w = fsuth(T_w)
          ! rwm(0) = rwm(1)
          ! twm(0) = twm(1)
       ! elseif (iwall==2) then  ! isothermal
          rwall = rwm(1)*twm(1)/T_w
          rwm(0) = 2.0_wp*rwall - rwm(1)
          twm(0) = 2.0_wp*T_w - twm(1)
       ! endif
       rwm(nwm+1) = r_e        ! dirichlet
       twm(nwm+1) = T_e
       !- Check if urms/trms < 0.0001% free stream values
       if (urms<res_wm.and.trms<res_wm) then
          tau_w = mu_w*uwm(1)/ycwm(1)
          ! q_w   = max( 0.0_wp, lambdaref*(twm(1)-T_w)/ycwm(1) )
          q_w   = max(0.0_wp, u_e*tau_w + cpav*mutwm(nwm)*prti*(twm(nwm) - twm(nwm-1))/(ycwm(nwm) - ycwm(nwm-1)))
          ! do kwm=0,nwm+1
          !    write(44,*) ycwm(kwm), uwm(kwm)/u_e, Twm(kwm)/T_e
          ! enddo
          return
       endif

    enddo

    call mpistop('Problem in wall-model ODE resolution...', 1)

  end subroutine wm_ODE

subroutine tdma( n, a, b, c, d, x )
!-------------------------------------------------------------------------------
! Solves a tri-diagonal matrix system
! n - No. of columns/rows
! a - Lower diagonal
! b - Main  diagonal
! c - Upper diagonal
! d - Right hand side
! x - Solution vector
!-------------------------------------------------------------------------------
   implicit none
   !-------------------------------------------------------------------------
   ! Input/Output arguments
   integer , intent(in)  :: n
   real(wp), intent(in)  :: a(n), c(n), b(n), d(n)
   real(wp), intent(out) :: x(n)
   ! Local variables
   integer  :: i
   real(wp) :: m, cp(n), dp(n)
   !-------------------------------------------------------------------------
   cp(1) = c(1)/b(1)
   dp(1) = d(1)/b(1)
   !- Forward elimination
   do i = 2,n
     m = b(i)-cp(i-1)*a(i)
     cp(i) = c(i)/m
     dp(i) = (d(i)-dp(i-1)*a(i))/m
   enddo
   !- Back substitution
   x(n) = dp(n)
   do i = n-1,1,-1
     x(i) = dp(i)-cp(i)*x(i+1)
   enddo
end subroutine tdma

!=================================================================================================================================================================
!                                             Algebraic wall model supplemented with the Smagorinsky model at wall
!================================================================================================================================================================

 !==============================================================================
subroutine bc_wm_wale_jmin
    !==============================================================================
      !> Use algebraic model
    !==============================================================================
      implicit none
      ! ---------------------------------------------------------------------------
      integer :: i,k,j0,j1
      real(wp) :: tau_w,u_e,T_e,r_e,T_w,r_w,mu_w,hwm,utau_wm, q_w
      ! ---------------------------------------------------------------------------
      real(wp) :: visc_sgs, sum_gij, visco_sgs, deltac
      real(wp) :: S11_w, S22_w, S33_w, S12_w, S13_w, S23_w
      real(wp), parameter :: coeff = sqrt(10.6_wp)
      ! ---------------------------------------------------------------------------
      external :: filter_sij 
      
    deltac = (idx(1)*idy(1)*idz(1))**(-1./3.)
    
      j0 = 1; j1 = j0+wm_ind
      hwm = abs(y(j1) - y(j0))
  
      do k=1,nz
         do i=1,nx
            ! Give first guess for Utau
            utau_wm = utau_jmin(i,k)
            ! Fill edge values
            u_e = (uu(i,j1,k)**2 + ww(i,j1,k)**2)**0.5_wp
            ! u_e = uu(i,j1,k)
            T_e = Tmp(i,j1,k)
            r_e = rho(i,j1,k)
            ! Fill wall values
            T_w = Tmp(i,j0,k)
            r_w = rho(i,j0,k)
            mu_w = visc(i,j0,k)
  
            ! Only implemented with spalding for the moment
            call wm_alg_spalding(u_e,T_e,r_e,T_w,r_w,mu_w,hwm,utau_wm,q_w)
  
            ! Store new value
            utau_jmin(i,k) = utau_wm
  
            ! Calculation of tau_w
            tau_w = rho(i,j0,k)*utau_jmin(i,k)**2
  

            ! Impose wall-model values
            Frhou(i,j0,k) = 0.0_wp; Frhow(i,j0,k) = 0.0_wp
            Grhov(i,j0,k) = 0.0_wp; Grhow(i,j0,k) = 0.0_wp
            Hrhou(i,j0,k) = 0.0_wp; Hrhov(i,j0,k) = 0.0_wp; Hrhow(i,j0,k) = 0.0_wp; Hrhoe(i,j0,k) = 0.0_wp
            Frhov(i,j0,k) = -dir_jmin*tau_w
            Grhou(i,j0,k) = -dir_jmin*tau_w
            Frhoe(i,j0,k) = -dir_jmin*q_w
            Grhoe(i,j0,k) = -dir_jmin*q_w
         enddo
      enddo
    !==============================================================================
    !> Supplement with Smago at wall: 
    !==============================================================================
    do k=1,nz 
        do i=1,nx 
            ! Compute S_ij at wall: 
            call filtre_Sij
            S11_w   = S11f(i,j0,k)
            S22_w   = S22f(i,j0,k)
            S33_w   = S33f(i,j0,k)
            S12_w   = S12f(i,j0,k)
            S13_w   = S13f(i,j0,k)
            S23_w   = S23f(i,j0,k)
            ! Get nu_sgs: 
            visco_sgs       = rho(i,j0,k)*((Cs_SM*deltac)**2)*sqrt(2.*(S11_w**2+S22_w**2+S33_w**2 &
                                              + 2.*(S12_w**2 +S13_w**2+S23_w**2)))
            Frhov(i,j0,k) = Frhov(i,j0,k) - visc_sgs* 2.0_wp * S12_w
            Grhou(i,j0,k) = Grhou(i,j0,k) - visc_sgs* 2.0_wp * S12_w
        enddo
    enddo
    end subroutine bc_wm_wale_jmin
  
    !==============================================================================
    subroutine bc_wm_wale_jmax
    !==============================================================================
      !> Use algebraic model
    !==============================================================================
      implicit none
      ! ---------------------------------------------------------------------------
      integer :: i,k,j0,j1
      real(wp) :: tau_w,u_e,T_e,r_e,T_w,r_w,mu_w,hwm,utau_wm, q_w
      ! ---------------------------------------------------------------------------
      real(wp) :: visc_sgs, sum_gij, visco_sgs, deltac
      real(wp) :: S11_w, S22_w, S33_w, S12_w, S13_w, S23_w
      real(wp), parameter :: coeff = sqrt(10.6_wp)
      external :: filtre_Sij
      ! ! Test stochastic forcing
      ! real(wp) :: tau11,tau22,tau33,tau12,tau13,tau23,trace,mu
      ! real(wp) :: coeff_kappa,var,yp_
      ! real(wp) :: dupx,dupy,dupz,dvpx,dvpy,dvpz,dwpx,dwpy,dwpz
      ! real(wp), parameter :: sp_kappa = 0.41_wp
      ! ! random number generator
      ! real(wp) :: randomnb
      ! integer :: nseed
      ! integer, dimension(:), allocatable :: initseed
    deltac = (idx(1)*idy(1)*idz(1))**(-1./3.)
        
      j0 = ny; j1 = j0-wm_ind
      hwm = abs(y(j1) - y(j0))
  
      do k=1,nz
         do i=1,nx
            ! Give first guess for Utau
            utau_wm = utau_jmax(i,k)
            ! Fill edge values
            u_e = (uu(i,j1,k)**2 + ww(i,j1,k)**2)**0.5_wp
            ! u_e = uu(i,j1,k)
            T_e = Tmp(i,j1,k)
            r_e = rho(i,j1,k)
            ! Fill wall values
            T_w = Tmp(i,j0,k)
            r_w = rho(i,j0,k)
            mu_w = visc(i,j0,k)
  
            ! Only implemented with spalding for the moment
            call wm_alg_spalding(u_e,T_e,r_e,T_w,r_w,mu_w,hwm,utau_wm,q_w)
  
            ! Store new value
            utau_jmax(i,k) = utau_wm
  
            ! Calculation of tau_w
            tau_w = rho(i,j0,k)*utau_jmax(i,k)**2
  
            ! Impose wall-model values
            Frhou(i,j0,k) = 0.0_wp; Frhow(i,j0,k) = 0.0_wp
            Grhov(i,j0,k) = 0.0_wp; Grhow(i,j0,k) = 0.0_wp
            Hrhou(i,j0,k) = 0.0_wp; Hrhov(i,j0,k) = 0.0_wp; Hrhow(i,j0,k) = 0.0_wp; Hrhoe(i,j0,k) = 0.0_wp
            Frhov(i,j0,k) = -dir_jmax*tau_w
            Grhou(i,j0,k) = -dir_jmax*tau_w
            Frhoe(i,j0,k) = -dir_jmax*q_w
            Grhoe(i,j0,k) = -dir_jmax*q_w
         enddo
      enddo

    !==============================================================================
    !> Supplement with Smago at wall: 
    !==============================================================================
    do k=1,nz 
        do i=1,nx 
            ! Compute S_ij at wall: 
            call filtre_Sij
            S11_w   = S11f(i,j0,k)
            S22_w   = S22f(i,j0,k)
            S33_w   = S33f(i,j0,k)
            S12_w   = S12f(i,j0,k)
            S13_w   = S13f(i,j0,k)
            S23_w   = S23f(i,j0,k)
            ! Get nu_sgs: 
            visco_sgs   = rho(i,j0,k)*((Cs_SM*deltac)**2)*sqrt(2.*(S11_w**2+S22_w**2+S33_w**2 &
                                              + 2.*(S12_w**2 +S13_w**2+S23_w**2)))
            Frhov(i,j0,k) = Frhov(i,j0,k) - visc_sgs* 2.0_wp * S12_w
            Grhou(i,j0,k) = Grhou(i,j0,k) - visc_sgs* 2.0_wp * S12_w
        enddo
    enddo
    end subroutine bc_wm_wale_jmax
  

end module mod_wall_model
