!=================================================================================
module mod_vortex_model
!=================================================================================
  !> Module to define vortex models
!=================================================================================
  use mod_grid     ! <- spatial coordinates
  use mod_constant ! <- constant pi
  implicit none
  !-------------------------------------------------------------------------------
  ! vortex models
  integer, parameter :: Taylor   = 1
  integer, parameter :: Oseen    = 2
  integer, parameter :: Scully   = 3
  integer, parameter :: Batchelor= 4
  integer, parameter :: shock    = 5
  !-------------------------------------------------------------------------------

contains
  
  !==============================================================================
  subroutine vortex_2d_model(m,xv,yv,u_vortex,v_vortex,p_vortex)
  !==============================================================================
    !> Add vortex
    !==============================================================================
    use mod_eos ! <- constant gam
    use warnstop
    implicit none
    ! ---------------------------------------------------------------------------
    ! vortex type
    integer, intent(in) :: m
    ! vortex location (center)
    real(wp), intent(in) :: xv,yv
    ! vortex velocities & pressure
    real(wp), dimension(nx,ny), intent(out) :: u_vortex,v_vortex,p_vortex
    ! ---------------------------------------------------------------------------
    integer :: i,j
    real(wp) :: ar,ampl
    real(wp) :: beta,nu,r0,circ,r,utheta
    real(wp) :: Ap,Bp,Ei1,Ei2
    real(wp) :: uc,rc,c1,c2
    ! ---------------------------------------------------------------------------

    ! Define vortex
    ! =============
    select case (m)

    case (Taylor) ! Taylor's vortex
                  ! ---------------    

       ! Gaussian half-width
       !ar=log(2.0_wp)/(9.0_wp*deltax**2)
       !ar=log(2.0_wp)/9.0_wp
       ar=log(2.0_wp)/16.0_wp
       ! ar=log(2.0_wp)/(25.0_wp*deltax**2)
       ! ar=log(2.0_wp)/(10.0_wp*deltax)**2
       !ar=log(2.0_wp)/(2.0_wp*deltax)**2
       ! ar=log(2.0_wp)/(6.0_wp*deltax)**2
       
       ! amplitude
       ampl=5.0_wp
       ! ampl=1.0_wp
       
       if (is_curv3) then
          do j=1,ny
             do i=1,nx
                u_vortex(i,j)= ampl*(yc3(i,j,1)-yv)/deltax   &
                     * exp(-ar*((xc3(i,j,1)-xv)**2+(yc3(i,j,1)-yv)**2))
                v_vortex(i,j)=-ampl*(xc3(i,j,1)-xv)/deltax   &
                     * exp(-ar*((xc3(i,j,1)-xv)**2+(yc3(i,j,1)-yv)**2))
                p_vortex(i,j)=-rho_ref*ampl**2/(4.0_wp*ar*(deltax)**2) &
                     * exp(-2.0_wp*ar*((xc3(i,j,1)-xv)**2+(yc3(i,j,1)-yv)**2))
                enddo
             enddo
       else
          if (is_curv) then
             do j=1,ny
                do i=1,nx
                   u_vortex(i,j)= ampl*(yc(i,j)-yv)/deltax   &
                        * exp(-ar*((xc(i,j)-xv)**2+(yc(i,j)-yv)**2))
                   v_vortex(i,j)=-ampl*(xc(i,j)-xv)/deltax   &
                        * exp(-ar*((xc(i,j)-xv)**2+(yc(i,j)-yv)**2))
                   p_vortex(i,j)=-rho_ref*ampl**2/(4.0_wp*ar*(deltax)**2) &
                        * exp(-2.0_wp*ar*((xc(i,j)-xv)**2+(yc(i,j)-yv)**2))
                enddo
             enddo
          else
             do j=1,ny
                do i=1,nx
                   u_vortex(i,j)= ampl*(y(j)-yv)/deltax   &
                        * exp(-ar*((x(i)-xv)**2+(y(j)-yv)**2))
                   v_vortex(i,j)=-ampl*(x(i)-xv)/deltax   &
                        * exp(-ar*((x(i)-xv)**2+(y(j)-yv)**2))
                   p_vortex(i,j)=-rho_ref*ampl**2/(4.0_wp*ar*(deltax)**2) &
                        * exp(-2.0_wp*ar*((x(i)-xv)**2+(y(j)-yv)**2))
                enddo
             enddo
          endif
       endif

    case (Oseen) ! Lamb-Oseen's vortex
                 ! -------------------    

       ! parameters
       beta=0.7272_wp
       nu=16.0_wp ! *deltax
       r0=nu/4.0_wp
       ar=1.256431_wp
       circ=pi*c_ref*r0/4.0_wp/beta
       
       if (is_curv) then
          do j=1,ny
             do i=1,nx
                r=sqrt((xc(i,j)-xv)**2+(yc(i,j)-yv)**2)
                if (r==0.0_wp) r=1.0e-6
                utheta=circ/(2*pi*r)*(1.0_wp-exp(-ar*(r/r0)**2))
                u_vortex(i,j)=-utheta*(yc(i,j)-yv)/r
                v_vortex(i,j)= utheta*(xc(i,j)-xv)/r
                !p_vortex(i,j) = p_ref
                Ap=ar*(r/r0)**2
                call e1xb(2.0_wp*Ap,Ei1)
                call e1xb(Ap,Ei2)
                Bp=0.5-exp(-Ap)+0.5*exp(-2*Ap)-Ap*Ei1+Ap*Ei2
                p_vortex(i,j)=-p_ref*((gam-1)*circ**2/4/c_ref**2/pi**2/r**2*Bp)**(gam/(gam-1))                
             enddo
          enddo
       else
          do j=1,ny
             do i=1,nx
                r=sqrt((x(i)-xv)**2+(y(j)-yv)**2)
                if (r==0.0_wp) r=1.0e-6
                utheta=circ/(2*pi*r)*(1.0_wp-exp(-ar*(r/r0)**2))
                u_vortex(i,j)=-utheta*(y(j)-yv)/r
                v_vortex(i,j)= utheta*(x(i)-xv)/r
                !p_vortex(i,j) = p_ref
                Ap=ar*(r/r0)**2
                call e1xb(2.0_wp*Ap,Ei1)
                call e1xb(Ap,Ei2)
                Bp=0.5-exp(-Ap)+0.5*exp(-2*Ap)-Ap*Ei1+Ap*Ei2
                p_vortex(i,j)=-p_ref*((gam-1)*circ**2/4/c_ref**2/pi**2/r**2*Bp)**(gam/(gam-1))                
             enddo
          enddo
       endif
       
    case (Scully) ! Scully's vortex
                  ! ---------------
       
       call mpistop('Scully''s vortex not yet defined', 0)
       
    case (Batchelor) ! Batchelor's vortex
                     ! ------------------
       
       call mpistop('Batchelor''s vortex not yet defined', 0)
       
    case (shock) ! shock-vortex interaction
                 ! ------------------------

       ! parameters
       uc=0.25_wp 
       rc=0.075_wp 
       c1=uc/rc*exp(0.5_wp) 
       c2=0.5_wp/(rc**2)

       if (is_curv) then
          do j=1,ny
             do i=1,nx
                r=sqrt((xc(i,j)-xv)**2+(yc(i,j)-yv)**2)
                if (r==0.0_wp) r=1.0e-6_wp
                utheta=u_ref*c1*r*exp(-c2*r**2)
                u_vortex(i,j)=-utheta*(yc(i,j)-yv)/r
                v_vortex(i,j)= utheta*(xc(i,j)-xv)/r
                p_vortex(i,j)= 0.0_wp              
             enddo
          enddo
       else
          do j=1,ny
             do i=1,nx
                r=sqrt((x(i)-xv)**2+(y(j)-yv)**2)
                if (r==0.0_wp) r=1.0e-6_wp
                utheta=u_ref*c1*r*exp(-c2*r**2)
                u_vortex(i,j)=-utheta*(y(j)-yv)/r
                v_vortex(i,j)= utheta*(x(i)-xv)/r
                p_vortex(i,j)= 0.0_wp       
             enddo
          enddo
       endif
       
    case default
       
       call mpistop('vortex model not defined', 0)
       
    end select

  end subroutine vortex_2d_model

  !==============================================================================
  subroutine vortex_3d_model(m,xv,yv,u_vortex,v_vortex,p_vortex)
  !==============================================================================
    !> Add vortex
    !==============================================================================
    use mod_eos ! <- constant gam
    use warnstop
    implicit none
    ! ---------------------------------------------------------------------------
    ! vortex type
    integer, intent(in) :: m
    ! vortex location (center)
    real(wp), intent(in) :: xv,yv
    ! vortex velocities & pressure
    real(wp), dimension(nx,ny,nz), intent(out) :: u_vortex,v_vortex,p_vortex
    ! ---------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: ar,ampl
    ! ---------------------------------------------------------------------------

    ! Define vortex
    ! =============
    select case (m)

    case (Taylor) ! Taylor's vortex
                  ! ---------------    

       ! Gaussian half-width
       !ar=log(2.0_wp)/(9.0_wp*deltax**2)
       !ar=log(2.0_wp)/9.0_wp
       ar=log(2.0_wp)/16.0_wp
       ! ar=log(2.0_wp)/(25.0_wp*deltax**2)
       ! ar=log(2.0_wp)/(10.0_wp*deltax)**2
       !ar=log(2.0_wp)/(2.0_wp*deltax)**2
       ! ar=log(2.0_wp)/(6.0_wp*deltax)**2
       
       ar=log(2.0_wp)/1.5_wp

       ! amplitude
       ampl=5.0_wp
       ! ampl=1.0_wp
       
       if (is_curv3) then
          do k=1,nz
             do j=1,ny
                do i=1,nx
                   u_vortex(i,j,k)= ampl*(yc3(i,j,k)-yv)/deltax   &
                        * exp(-ar*((xc3(i,j,k)-xv)**2+(yc3(i,j,k)-yv)**2))
                   v_vortex(i,j,k)=-ampl*(xc3(i,j,k)-xv)/deltax   &
                        * exp(-ar*((xc3(i,j,k)-xv)**2+(yc3(i,j,k)-yv)**2))
                   p_vortex(i,j,k)=-rho_ref*ampl**2/(4.0_wp*ar*(deltax)**2) &
                        * exp(-2.0_wp*ar*((xc3(i,j,k)-xv)**2+(yc3(i,j,k)-yv)**2))
                enddo
             enddo
          enddo
       endif
       
    case default
       
       call mpistop('vortex model not defined', 0)
       
    end select

  end subroutine vortex_3d_model

end module mod_vortex_model
