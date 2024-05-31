!===========================================================================
module mod_turb_inlet_outlet
!===========================================================================
  !> authors: OY,CM
  !> date: dunno
  !> Turbomachine subsonic inlet and outlet boundary conditions
  !> Att. for perfect and PRS gas only
!===========================================================================
  use mod_mpi
  use mod_flow
  use mod_eos
  use mod_tranprop
  use mod_constant
  use mod_interface
  implicit none
  ! ---------------------------------------------------------------------------
  integer  :: i,j,k
  ! ---------------------------------------------------------------------------
  ! PFG
  ! ---------------------------------------------------------------------------
  real(wp) :: vel_1,vel_2,vel_3,Rp,Rm,atmp,btmp,ctmp,dtmp,ro1, &
              mach_1,c02,vc,rovc,vel_b,costheta,vel_int(3),a
  real(wp) :: sinf,cinf,Vinf(3),Vt(3),s,Vn_2,V_2(3),Vtmp,V_1(3)
  ! ---------------------------------------------------------------------------
  ! PRS
  ! ---------------------------------------------------------------------------
  ! inlet
  real(wp) :: nb(3),veltmp(3),theta ! direction and velocity of flow at inlet
  real(wp) :: nbc(3) ! boundary normal
  real(wp) :: p_tot,T_tot,rho_tot ! total thermos
  real(wp) :: ptent,Ttent,rotent  ! tentative variables for PRS Newton iterations
  real(wp) :: ptmp,Ttmp,rhotmp  ! dummy variables to loop over
  real(wp) :: pbc,Tbc,rhobc ! thermos at boundary
  real(wp) :: velbc,veli(3),ci ! velocities at boundary and inside domain
  real(wp) :: err,err1,err2,tol ! convergence criterions
  real(wp) :: h,Htmp,H_tot ! enthalpy
  real(wp) :: rhoi,pi_ ! neighbour thermos
  real(wp) :: dpdrho,dHtmpdrho,dveldrho ! derivatives w.r.t density for Newton algorithm
  real(wp) :: dpdT,dHtmpdT,dveldT ! derivatives w.r.t temperature for Newton algorithm
  integer  :: count_
  ! outlet
  real(wp) :: machloc,machbc,thetaloc
  ! ---------------------------------------------------------------------------

contains

  ! ====== !
  ! Inlets !
  ! ====== !

  ! PFG
  ! ===

!!$  subroutine bc_inlet1_jmin
!!$
!!$    ! set total pressure and temperature
!!$    p_tot = 2.0_wp*p_ref
!!$    T_tot = T_ref/0.833_wp
!!$  
!!$    ! direction of the velocity at the boundary (WARNING: THIS SHOULD BE COME FROM THE INPUT FILE)
!!$    nb = (/0.0_wp,1.0_wp,0.0_wp/)
!!$  
!!$    ! Update inlet points
!!$    ! ==================  
!!$    do k=1,nz
!!$       do i=ndx_c,nfx_c
!!$       
!!$          ! velocity vector at the interior nodes
!!$          veltmp = (/uu(i,2,k),vv(i,2,k),ww(i,2,k)/)
!!$  
!!$          ! compute total enthalpy from the interior nodes
!!$          ! total enthalpy is conserved at boundary due to adiabatic assumption
!!$          H_tot = prs(i,2,k)/rho_n(i,2,k)*(gam/gam1)+0.5_wp*dot_product(veltmp,veltmp)
!!$          
!!$          ! exrapolate R+ from the interior nodes
!!$          vel_2 = dot_product(veltmp,nb)
!!$          Rp = -vel_2-2.0_wp*c_(i,2,k)/gam1
!!$          
!!$          ! obtain speed of sound at the boundary
!!$          atmp = 1.0_wp+2.0_wp/gam1
!!$          btmp = 2.0_wp*Rp
!!$          ctmp = 0.5_wp*gam1*(Rp**2-2.0_wp*H_tot)
!!$          dtmp = sqrt(btmp**2-4.0_wp*atmp*ctmp+1d-20)
!!$          c_(i,1,k) = max(-btmp+dtmp,-btmp-dtmp)/(2.0_wp*atmp)
!!$          
!!$          ! update velocity at the boundary
!!$          vel_1 = -Rp-2.0_wp*c_(i,1,k)/gam1
!!$          
!!$          ! compute static pressure and static temperature from specified total values
!!$          prs(i,1,k) = p_tot*(1.0_wp+0.5_wp*gam1*(vel_1/c_(i,1,k))**2)**(-gam/gam1)
!!$          Tmp(i,1,k) = T_tot*(prs(i,1,k)/p_tot)**(gam1/gam)
!!$  
!!$          ! calc rho
!!$          ro1 = rocalc_pt(prs(i,1,k),Tmp(i,1,k),rho_n(i,2,k) )
!!$          rho_n(i,1,k)  = ro1     
!!$          
!!$          ! compute velocity and momentum variables
!!$          uu(i,1,k) = vel_1*nb(1)
!!$          vv(i,1,k) = vel_1*nb(2)
!!$          ww(i,1,k) = vel_1*nb(3)
!!$          rhou_n(i,1,k) = ro1*uu(i,1,k)
!!$          rhov_n(i,1,k) = ro1*vv(i,1,k)
!!$          rhow_n(i,1,k) = ro1*ww(i,1,k)
!!$          
!!$          ! compute the rest
!!$          visc(i,1,k)=viscosity_law(Tmp(i,1,k),ro1)*diffscale
!!$          cok(i,1,k) =thconductivity(visc(i,1,k),Tmp(i,1,k),ro1)
!!$  !         c_(i,1,k)  = sqrt(c2calc_tro(Tmp(i,1,k),ro1))
!!$          rhoe_n(i,1,k) = ro1*(ecalc_tro(Tmp(i,1,k),ro1) &
!!$            + 0.5_wp*(uu(i,1,k)**2+vv(i,1,k)**2+ww(i,1,k)**2))
!!$          
!!$          Krho(i,1,k) = 0.0_wp
!!$          Krhou(i,1,k) = 0.0_wp
!!$          Krhov(i,1,k) = 0.0_wp
!!$          Krhow(i,1,k) = 0.0_wp
!!$          Krhoe(i,1,k) = 0.0_wp
!!$  
!!$       enddo
!!$    enddo
!!$  
!!$  end subroutine bc_inlet1_jmin
  
  subroutine bc_inlet2_jmin
  
    ! set total pressure and temperature
    p_tot = 1.6_wp*p_ref
    T_tot = T_ref
    
    ! unit normal vector of the boundary
    nb = (/-1.0_wp,0.0_wp,0.0_wp/)
    
    ! cos(theta) where theta is the flow angle relative to the boundary
    costheta = 1.0_wp
    
    ! Update inlet points
    ! ==================
    do k=1,nz
       do i=ndx_c,nfx_c
               
          ! velocity vector at the interior node
          vel_int = (/rhou_n(i,2,k)/rho_n(i,2,k),rhov_n(i,2,k)/rho_n(i,2,k),rhow_n(i,2,k)/rho_n(i,2,k)/)
          
          ! R- from the interior node
          Rm = dot_product(vel_int,nb)-2.0_wp*c_(i,2,k)/gam1
          
          ! obtain speed of sound at the boundary
          c02 = c_(i,2,k)**2+0.5_wp*gam1*dot_product(vel_int,vel_int)
          c_(i,1,k) = -Rm*gam1/(gam1*costheta**2+2.0_wp)*(1.0_wp+costheta*sqrt((gam1*costheta**2+2.0_wp)*c02/(gam1*Rm**2)-0.5_wp*gam1))
          
          ! compute static pressure and static temperature from specified total values
          Tmp(i,1,k) = T_tot*c_(i,1,k)**2/c02
          prs(i,1,k) = p_tot*(Tmp(i,1,k)/T_tot)**(gam/gam1)
  !         if(i==25) print*, prs(i,1,k),Tmp(i,1,k)
          
          !calc rho
          ro1 = rocalc_pt(prs(i,1,k),Tmp(i,1,k),rho_n(i,2,k) )
          rho_n(i,1,k)  = ro1
  
          ! compute velocity at the boundary
          vel_b = sqrt(2.0_wp*cpfg*(T_tot-Tmp(i,1,k)+1d-10))
          
  !         if(i==25) print*, c_(i,1,k),vel_b,prs(i,1,k),Tmp(i,1,k),ro1
          
          ! compute velocity and momentum variables
          uu(i,1,k) = -vel_b*nb(1)
          vv(i,1,k) = -vel_b*nb(2)
          ww(i,1,k) = -vel_b*nb(3)
          rhou_n(i,1,k) = ro1*uu(i,1,k)
          rhov_n(i,1,k) = ro1*vv(i,1,k)
          rhow_n(i,1,k) = ro1*ww(i,1,k)
          
          ! compute the rest
          visc(i,1,k)=viscosity_law(Tmp(i,1,k),ro1)*diffscale
          cok(i,1,k) =thconductivity(visc(i,1,k),Tmp(i,1,k),ro1)
  !         c_(i,1,k)  = sqrt(c2calc_tro(Tmp(i,1,k),ro1))
          rhoe_n(i,1,k) = ro1*(ecalc_tro(Tmp(i,1,k),ro1) &
            + 0.5_wp*(uu(i,1,k)**2+vv(i,1,k)**2+ww(i,1,k)**2))
          
          Krho(i,1,k) = 0.0_wp
          Krhou(i,1,k) = 0.0_wp
          Krhov(i,1,k) = 0.0_wp
          Krhow(i,1,k) = 0.0_wp
          Krhoe(i,1,k) = 0.0_wp
          
       enddo
    enddo
  
  end subroutine bc_inlet2_jmin

!!$subroutine bc_inlet1_imin
!!$  use mod_mpi
!!$  use mod_flow
!!$  use mod_eos
!!$  use mod_tranprop
!!$  use mod_constant
!!$  use mod_interface
!!$  use mod_time
!!$  implicit none
!!$  ! ---------------------------------------------------------------------------
!!$  integer :: i,j,k
!!$  real(wp) :: T_tot,p_tot,nb(3),veltmp(3),vel_2,vel_3,H_tot,Rp,atmp,btmp,ctmp,dtmp,vel_1,ro1,mach_1,c02,vc,rovc
!!$  real(wp) :: rho_tot,sinf,cinf,Vinf(3),Vt(3),Rm,s,Vn_2,V_2(3),Vtmp,V_1(3)
!!$  ! ---------------------------------------------------------------------------
!!$
!!$  !p_in=0.87*pref
!!$  !M_in=(2./gam1)*((pref/p_in)**(gam1/gamma)-1.);
!!$  !M_in=sqrt(M_in)
!!$  !T_in=T_infty*(p_in/pref)**(gam1/gamma)
!!$  !rho_in=p_in/(rg*T_in)
!!$
!!$  ! set total pressure and temperature
!!$  !p_tot = 250000_wp!*p_ref
!!$  p_tot = p_ref*1.
!!$  T_tot = T_ref!/0.833_wp
!!$
!!$  ! direction of the velocity at the boundary (WARNING: THIS SHOULD BE COME FROM THE INPUT FILE)
!!$  !nb = (/1.0_wp,0.0_wp,0.0_wp/)
!!$  !nb = (/cos(36.0_wp/180.0_wp*pi),sin(36.0_wp/180.0_wp*pi),0.0_wp/)
!!$  nb = (/cos(30.0_wp/180.0_wp*pi),sin(30.0_wp/180.0_wp*pi),0.0_wp/)
!!$
!!$  ! Update inlet points
!!$  ! ==================  
!!$  do k=1,nz
!!$     do j=1,ny
!!$     
!!$        ! velocity vector at the interior nodes
!!$        veltmp = (/uu(2,j,k),vv(2,j,k),ww(2,j,k)/)
!!$
!!$        ! compute total enthalpy from the interior nodes
!!$        ! total enthalpy is conserved at boundary due to adiabatic assumption
!!$        H_tot = prs(2,j,k)/rho_n(2,j,k)*(gam/gam1)+0.5_wp*dot_product(veltmp,veltmp)
!!$        
!!$        ! exrapolate R+ from the interior nodes
!!$        vel_2 = dot_product(veltmp,nb)
!!$        Rp = -vel_2-2.0_wp*c_(2,j,k)/gam1
!!$        
!!$        ! obtain speed of sound at the boundary
!!$        atmp = 1.0_wp+2.0_wp/gam1
!!$        btmp = 2.0_wp*Rp
!!$        ctmp = 0.5_wp*gam1*(Rp**2-2.0_wp*H_tot)
!!$        dtmp = sqrt(btmp**2-4.0_wp*atmp*ctmp+1d-20)
!!$        c_(1,j,k) = max(-btmp+dtmp,-btmp-dtmp)/(2.0_wp*atmp)
!!$        
!!$        ! update velocity at the boundary
!!$        vel_1 = -Rp-2.0_wp*c_(1,j,k)/gam1
!!$        
!!$        ! compute static pressure and static temperature from specified total values
!!$        prs(1,j,k) = p_tot*(1.0_wp+0.5_wp*gam1*(vel_1/c_(1,j,k))**2)**(-gam/gam1)
!!$        Tmp(1,j,k) = T_tot*(prs(1,j,k)/p_tot)**(gam1/gam)
!!$
!!$        ! calc rho
!!$        ro1 = rocalc_pt(prs(1,j,k),Tmp(1,j,k),rho_n(2,j,k) )
!!$        rho_n(1,j,k)  = ro1     
!!$        
!!$        ! compute velocity and momentum variables
!!$        uu(1,j,k) = vel_1*nb(1)
!!$        vv(1,j,k) = vel_1*nb(2)
!!$        ww(1,j,k) = vel_1*nb(3)
!!$        rhou_n(1,j,k) = ro1*uu(1,j,k)
!!$        rhov_n(1,j,k) = ro1*vv(1,j,k)
!!$        rhow_n(1,j,k) = ro1*ww(1,j,k)
!!$        
!!$        ! compute the rests
!!$        !visc(1,j,k)=viscosity_law(Tmp(1,j,k),ro1)*diffscale
!!$        !cok(1,j,k) =thconductivity(visc(1,j,k),Tmp(1,j,k),ro1)
!!$!         c_(1,j,k)  = sqrt(c2calc_tro(Tmp(1,j,k),ro1))
!!$        rhoe_n(1,j,k) = ro1*(ecalc_tro(Tmp(1,j,k),ro1) &
!!$                      + 0.5_wp*(uu(1,j,k)**2+vv(1,j,k)**2+ww(1,j,k)**2))
!!$        
!!$        Krho(1,j,k) = 0.0_wp
!!$        Krhou(1,j,k) = 0.0_wp
!!$        Krhov(1,j,k) = 0.0_wp
!!$        Krhow(1,j,k) = 0.0_wp
!!$        Krhoe(1,j,k) = 0.0_wp
!!$
!!$     enddo
!!$  enddo
!!$
!!$  !if (mod(ntime,100)==0) print *,iproc,'p(1,15)/p_ref',prs(1,15,1)/p_ref
!!$  
!!$end subroutine bc_inlet1_imin

  subroutine bc_inlet1_imin
    use mod_time      
    real(wp) :: V2,vel_t,c2t

    ! set total pressure and temperature
!!    T_tot = 287.00_wp
!!    !p_ref = pcalc_tro(T_ref,rho_ref)
!!    !p_tot = p_ref*(1.0_wp+gam1/2.0_wp*Mach**2)**(gam/gam1)
!!    p_tot=129000.0_wp

    T_tot=T_ref
    p_tot=p_ref
    
    !p_tot = p_ref*(1.0_wp+gam1/2.0_wp*Mach**2)**(gam/gam1)
    !T_tot = T_ref*(1.0_wp+gam1/2.0_wp*Mach**2)
!!$    print *,p_tot,T_tot,p_ref,T_ref
!!$    stop

    ! boundary normal
    nb = (/1.0_wp,0.0_wp,0.0_wp/)
  
    ! direction of the velocity at the boundary
    theta = 30.0_wp
    !theta = 0.0_wp
    theta = theta*pi/180.0_wp
    nbc = (/cos(theta),sin(theta),0.0_wp/)
    !nb=nbc
    c2t=cos(theta)**2
  
    ! Update inlet points
    ! ===================
    do k=1,nz
       !do j=ndy_c,nfy_c
       do j=1,ny
       
!!$          ! velocity vector at the interior nodes
!!$          veltmp = (/uu(2,j,k),vv(2,j,k),ww(2,j,k)/)
!!$          V2=uu(2,j,k)*uu(2,j,k)+vv(2,j,k)*vv(2,j,k)+ww(2,j,k)*ww(2,j,k)
!!$  
!!$          ! compute total enthalpy from the interior nodes
!!$          ! total enthalpy is conserved at boundary due to adiabatic assumption
!!$          H_tot = prs(2,j,k)/rho_n(2,j,k)*(gam/gam1)+0.5_wp*V2
!!$          
!!$          ! exrapolate R+ from the interior nodes
!!$          !vel_2 = uu(2,j,k)*nb(1)+vv(2,j,k)*nb(2)+ww(2,j,k)*nb(3)
!!$          vel_2 = uu(2,j,k)*nb(1)+vv(2,j,k)*nb(2)+ww(2,j,k)*nb(3)
!!$          vel_t = sqrt(V2-vel_2**2)
!!$          Rp = vel_2-2.0_wp*c_(2,j,k)/gam1
!!$          
!!$          ! obtain speed of sound at the boundary
!!$          atmp = 1.0_wp+2.0_wp/gam1
!!$          btmp = 2.0_wp*Rp
!!$          ctmp = 0.5_wp*gam1*(Rp**2+vel_t**2-2.0_wp*H_tot)
!!$          dtmp = sqrt(btmp**2-4.0_wp*atmp*ctmp+1d-20)
!!$          c_(1,j,k) = max(-btmp+dtmp,-btmp-dtmp)/(2.0_wp*atmp)
!!$          
!!$          ! update velocity at the boundary
!!$          vel_1 =Rp+2.0_wp*c_(1,j,k)/gam1
!!$          !vel_1 = vel_1/cos(theta)
!!$          vel_1 = sqrt(vel_1**2+vel_t**2)
          
          ! velocity vector at the interior nodes
          veltmp = (/uu(2,j,k),vv(2,j,k),ww(2,j,k)/)
          V2=uu(2,j,k)*uu(2,j,k)+vv(2,j,k)*vv(2,j,k)+ww(2,j,k)*ww(2,j,k)
  
          ! compute total enthalpy from the interior nodes
          ! total enthalpy is conserved at boundary due to adiabatic assumption
          H_tot = prs(2,j,k)/rho_n(2,j,k)*(gam/gam1)+0.5_wp*V2
          
          ! exrapolate R+ from the interior nodes
          !vel_2 = uu(2,j,k)*nb(1)+vv(2,j,k)*nb(2)+ww(2,j,k)*nb(3)
          vel_2 = uu(2,j,k)*nb(1)+vv(2,j,k)*nb(2)+ww(2,j,k)*nb(3)
          !vel_t = sqrt(V2-vel_2**2)
          Rp = vel_2-2.0_wp*c_(2,j,k)/gam1
          
          ! obtain speed of sound at the boundary
!!$          atmp = c2t+2.0_wp/gam1
!!$          btmp = 2.0_wp*Rp
!!$          ctmp = 0.5_wp*gam1*(Rp**2-2.0_wp*c2t*H_tot)
!!$          dtmp = sqrt(btmp**2-4.0_wp*atmp*ctmp+1d-20)
!!$          c_(1,j,k) = max(-btmp+dtmp,-btmp-dtmp)/(2.0_wp*atmp)

          atmp=c2t+2.0_wp/gam1
          c_(1,j,k)=-Rp/atmp*(1+cos(theta)*sqrt(gam1*(atmp*H_tot/Rp**2-0.5_wp)))

          
          ! Old update
          ! ==========
!!$          ! update velocity at the boundary
!!$          vel_1 =Rp+2.0_wp*c_(1,j,k)/gam1
!!$          !vel_1 = sqrt(vel_1**2+vel_t**2)
!!$          vel_1 = vel_1/cos(theta)
!!$          
!!$          ! compute static pressure and static temperature from specified total values
!!$          prs(1,j,k) = p_tot*(1.0_wp+0.5_wp*gam1*(vel_1/c_(1,j,k))**2)**(-gam/gam1)
!!$          Tmp(1,j,k) = T_tot*(prs(1,j,k)/p_tot)**(gam1/gam)
!!$  
!!$          ! calc rho
!!$          ro1 = rocalc_pt(prs(1,j,k),Tmp(1,j,k),rho_n(2,j,k) )

          ! New update
          ! ==========
          Tmp(1,j,k)=T_tot*(c_(1,j,k)**2/gam1/H_tot)
          prs(1,j,k) = p_tot*(Tmp(1,j,k)/T_tot)**(gam/gam1)

          ! density
          ro1=prs(1,j,k)/rg/Tmp(1,j,k)

!!$          ! velocity modulus
!!$          vel_1=sqrt(2.0_wp*cpfg*(T_tot-Tmp(1,j,k)))          
          ! update velocity at the boundary
          vel_1 =Rp+2.0_wp*c_(1,j,k)/gam1
          !vel_1 = sqrt(vel_1**2+vel_t**2)
          vel_1 = vel_1/cos(theta)

          !if ((mod(ntotal,1000)==0).and.(irk==nrk).and.(j==25).and.(iproc<=3)) print *,'iproc',iproc,'M1',vel_1/c_(1,j,k)

          ! calc rho
          rho_n(1,j,k)  = ro1     
          
          ! compute velocity and momentum variables
          uu(1,j,k) = vel_1*nbc(1)
          vv(1,j,k) = vel_1*nbc(2)
          ww(1,j,k) = vel_1*nbc(3)
          rhou_n(1,j,k) = ro1*uu(1,j,k)
          rhov_n(1,j,k) = ro1*vv(1,j,k)
          rhow_n(1,j,k) = ro1*ww(1,j,k)
          
          ! compute the rest
          visc(1,j,k)=viscosity_law(Tmp(1,j,k),ro1)*diffscale
          cok(1,j,k) =thconductivity(visc(1,j,k),Tmp(1,j,k),ro1)
          rhoe_n(1,j,k) = ro1*(ecalc_tro(Tmp(1,j,k),ro1) &
            + 0.5_wp*(uu(1,j,k)**2+vv(1,j,k)**2+ww(1,j,k)**2))
          
          Krho(1,j,k) = 0.0_wp
          Krhou(1,j,k) = 0.0_wp
          Krhov(1,j,k) = 0.0_wp
          Krhow(1,j,k) = 0.0_wp
          Krhoe(1,j,k) = 0.0_wp
  
       enddo
    enddo
  
  end subroutine bc_inlet1_imin

  subroutine bc_inlet1_jmin

    ! set total pressure and temperature
    !T_tot=T_ref
    !p_tot=p_ref
    
    p_tot = p_ref*(1.0_wp+gam1/2.0_wp*Mach**2)**(gam/gam1)
    T_tot = T_ref*(1.0_wp+gam1/2.0_wp*Mach**2)

    ! boundary normal
    nb = (/1.0_wp,0.0_wp,0.0_wp/)
  
    ! direction of the velocity at the boundary
    !theta = 30.0_wp
    theta = 0.0_wp
    theta = theta*pi/180.0_wp
    nbc = (/cos(theta),sin(theta),0.0_wp/)
    nbc=nb
  
    ! Update inlet points
    ! ===================
    do k=1,nz
       do i=1,nx
       
          ! velocity vector at the interior nodes
          veltmp = (/uu(i,2,k),vv(i,2,k),ww(i,2,k)/)
  
          ! compute total enthalpy from the interior nodes
          ! total enthalpy is conserved at boundary due to adiabatic assumption
          H_tot = prs(i,2,k)/rho_n(i,2,k)*(gam/gam1)+0.5_wp*dot_product(veltmp,veltmp)
          
          ! exrapolate R+ from the interior nodes
          vel_2 = dot_product(veltmp,nb)
          Rp = -vel_2-2.0_wp*c_(i,2,k)/gam1
          
          ! obtain speed of sound at the boundary
          atmp = 1.0_wp+2.0_wp/gam1
          btmp = 2.0_wp*Rp
          ctmp = 0.5_wp*gam1*(Rp**2-2.0_wp*H_tot)
          dtmp = sqrt(btmp**2-4.0_wp*atmp*ctmp+1d-20)
          c_(i,1,k) = max(-btmp+dtmp,-btmp-dtmp)/(2.0_wp*atmp)
          
          ! update velocity at the boundary
          vel_1 = -Rp-2.0_wp*c_(i,1,k)/gam1
          
          ! compute static pressure and static temperature from specified total values
          prs(i,1,k) = p_tot*(1.0_wp+0.5_wp*gam1*(vel_1/c_(i,1,k))**2)**(-gam/gam1)
          Tmp(i,1,k) = T_tot*(prs(i,1,k)/p_tot)**(gam1/gam)
  
          ! calc rho
          ro1 = rocalc_pt(prs(i,1,k),Tmp(i,1,k),rho_n(i,2,k) )
          rho_n(i,1,k)  = ro1     
          
          ! compute velocity and momentum variables
          uu(i,1,k) = vel_1*nbc(1)
          vv(i,1,k) = vel_1*nbc(2)
          ww(i,1,k) = vel_1*nbc(3)
          rhou_n(i,1,k) = ro1*uu(i,1,k)
          rhov_n(i,1,k) = ro1*vv(i,1,k)
          rhow_n(i,1,k) = ro1*ww(i,1,k)
          
          ! compute the rest
          !visc(i,1,k)=viscosity_law(Tmp(i,1,k),ro1)*diffscale
          !cok(i,1,k) =thconductivity(visc(i,1,k),Tmp(i,1,k),ro1)
          rhoe_n(i,1,k) = ro1*(ecalc_tro(Tmp(i,1,k),ro1) &
            + 0.5_wp*(uu(i,1,k)**2+vv(i,1,k)**2+ww(i,1,k)**2))
          
          Krho(i,1,k) = 0.0_wp
          Krhou(i,1,k) = 0.0_wp
          Krhov(i,1,k) = 0.0_wp
          Krhow(i,1,k) = 0.0_wp
          Krhoe(i,1,k) = 0.0_wp
  
       enddo
    enddo
  
  end subroutine bc_inlet1_jmin

  ! PRS
  ! ===

  subroutine bc_inlet_prs_imin
  !==============================================================================
    !> boundary condition at imin for PRS eos
  !==============================================================================
  
    ! set total pressure and temperature
    !p_tot = 52200.0_wp ! Gori et al. discharge 2 experiment E2
    !T_tot = 502.88_wp
    p_tot = 129000.0_wp ! Gori et al. discharge 2 experiment D2
    T_tot = 504.33_wp
    rho_tot = rocalc_pt(p_tot,T_tot,rho_n(1,1,1))
    H_tot = ecalc_tro(T_tot,rho_tot) + p_tot/rho_tot
    s = scalc_tro(T_tot,rho_tot) ! s0 = sinf
  
    ! boundary normal
    nb = (/-1.0_wp,0.0_wp,0.0_wp/)
    
    ! direction of the velocity at the boundary 
    theta = 30.0_wp
    theta = theta*pi/180.0_wp
    nbc = (/-cos(theta),-sin(theta),0.0_wp/)
  
    ! Initialize
    tol = 1e-6_wp
    do k=1,nz
       do j=ndy_c,nfy_c
          err   = 1.0_wp
          err1  = 1.0_wp
          err2  = 1.0_wp
          count_= 0

          ! use neighbour point for test values
          rotent = rho_n(2,j,k)
          ptent  = prs(2,j,k)
          Ttent  = Tmp(2,j,k)

          ! boundary state
          rhobc = rho_n(2,j,k)
          pbc   = prs(2,j,k)
          Tbc   = Tmp(2,j,k)
  
          ! interior nodes
          veli = (/uu(2,j,k),vv(2,j,k),ww(2,j,k)/)
          ci = sqrt(c2calc_tro(Tmp(2,j,k),rho_n(2,j,k)))
          rhoi = rho_n(2,j,k)
          pi_ = prs(2,j,k)

          ! Newton iteration: rho(n+1) = rho(n) - (Htmp(rho(n)-H_tot))/dHtmp/drho
          ! Htmp(T,rho) = e(T,rho) + 0.5u**2(rho) + P(T,rho)/rho

          ! Newton iteration: T(n+1) = T(n) - (Htmp(T(n)-H_tot))/dHtmp/dT
          ! Htmp(T,rho) = e(T,rho) + 0.5u**2(rho) + P(T,rho)/rho
          !do while (err1.gt.tol.or.err2.gt.tol)
          do while (err2.gt.tol)

             ! temporary variables
             rhotmp = rhobc 
             Ttmp   = Tbc
             ptmp   = pbc

             ! d(P/rho)/drho = -rho*dP/dv - P/rho**2
             dpdrho = -rhotmp**2*dpdvcalc_tro(Ttmp,rhotmp)

             ! d(P/rho)/dT = dP/dT/rho
             dpdT = dpdicalc_tro(Ttmp,rhotmp)*cvcalc_tro(Ttmp,rhotmp)

             ! extrapolate R+ from the interior nodes, so dR=0
             ! dR = u_i-u_bc + c_i/rho_i*(rho_i-rho_bc)
             !velbc = dot_product(veli,nb) + ci/rhoi*(rhoi-rhotmp)
             ! alternative expression
             ! dR = u_i-u_bc + 1/(rho_i*c_i)(P_i-P_bc)
             velbc = dot_product(veli,nb) + 1.0_wp/rhoi/ci*(pi_-ptmp)

             ! velocity derivative w.r.t density
             dveldrho = -2.0_wp*ci/rhoi*(dot_product(veli,nb) + ci/rhoi*(rhoi-rhotmp))

             ! velocity derivative w.r.t temperature
             dveldT = -2.0_wp*velbc*dpdT/rhoi/ci

             ! total enthalpy derivative w.r.t density
             dHtmpdrho = dedrocalc_tro(Ttmp,rhotmp) + 0.5_wp*dveldrho + dpdrho/rhotmp-ptmp/rhotmp**2

             ! total enthalpy derivative w.r.t temperature
             !dHtmpdT = dedTcalc_tro(Ttmp,rhotmp) + dpdT/rhotmp
             dHtmpdT = dedTcalc_tro(Ttmp,rhotmp) + 0.5_wp*dveldT + dpdT/rhotmp

             ! total enthalpy
             Htmp = ecalc_tro(Ttmp,rhotmp) + 0.5_wp*velbc**2 + ptmp/rhotmp

             ! new thermos value
             !rhobc = rhotmp-(Htmp-H_tot)/dHtmpdrho
             !Tbc = tcalc_pro(ptmp,rhobc,Ttent)
             Tbc = Ttmp-(Htmp-H_tot)/dHtmpdT
             rhobc = rocalc_st(s,Tbc,rotent)
             pbc = pcalc_tro(Tbc,rhobc)

             ! check convergence 
             err1 = abs(rhobc-rhotmp)/rhotmp
             err2 = abs(Tbc-Ttmp)/Ttmp
             count_ = count_+1
             !print*,err1,err2,count_
          enddo
 
!          ! Loop over total quantities
!          do while (err.gt.tol)
!             ! initialise pressure
!             ptmp = pbc
!  
!             ! extrapolate R+ from the interior nodes, so dR=0
!             ! dR = u_i-u_bc + c_i/rho_i*(rho_i-rho_bc)
!             ci = sqrt(c2calc_tro(Tmp(2,j,k),rho_n(2,j,k)))
!             velbc = dot_product(veli,nb) + ci/rho_n(2,j,k)*(rho_n(2,j,k)-rhobc)
!     
!             ! compute static enthalpy
!             h = H_tot - 0.5_wp*velbc**2
!
!             ! update thermo variables
!             Tbc = tcalc_ph(ptmp,h,rotent,Ttent)
!             !rhobc = rocalc_pt(ptmp,Tbc,rotent)
!             !s = scalc_tro(Tbc,rhobc)
!             rhobc = rocalc_st(s,Tbc,rotent)
!             !Tbc = tcalc_sro(s,rhobc,Ttent)
!             pbc = pcalc_tro(Tbc,rhobc)
!         
!             ! check convergence 
!             err = abs(pbc-ptmp)/ptmp
!          enddo
  
          ! update thermo variables at boundary 
          rho_n(1,j,k) = rhobc
          prs(1,j,k)   = pbc 
          Tmp(1,j,k)   = Tbc 
          !rho_n(1,j,k) = rho_tot
          !prs(1,j,k)   = p_tot
          !Tmp(1,j,k)   = T_tot
  
          ! compute velocity and momentum variables
          uu(1,j,k) = velbc*nbc(1)
          vv(1,j,k) = velbc*nbc(2)
          ww(1,j,k) = velbc*nbc(3)
          rhou_n(1,j,k) = rho_n(1,j,k)*uu(1,j,k)
          rhov_n(1,j,k) = rho_n(1,j,k)*vv(1,j,k)
          rhow_n(1,j,k) = rho_n(1,j,k)*ww(1,j,k)
          
          ! compute the rest
          visc(1,j,k)= viscosity_law(Tmp(1,j,k),rho_n(1,j,k))
          cok(1,j,k) = thconductivity(visc(1,j,k),Tmp(1,j,k),rho_n(1,j,k))
          c_(1,j,k)  = sqrt(c2calc_tro(Tmp(1,j,k),rho_n(1,j,k)))
          rhoe_n(1,j,k) = rho_n(1,j,k)*(ecalc_tro(Tmp(1,j,k),rho_n(1,j,k)) &
            + 0.5_wp*(uu(1,j,k)**2+vv(1,j,k)**2+ww(1,j,k)**2))
          
          ! extrapolate to ghost points
          do a=1-ngh,1
             Tmp(a,j,k)    = Tmp(1,j,k)
             prs(a,j,k)    = prs(1,j,k)
             rho_n(a,j,k)  = rho_n(1,j,k)
             rhou_n(a,j,k) = rhou_n(1,j,k)
             rhov_n(a,j,k) = rhov_n(1,j,k)
             rhow_n(a,j,k) = rhow_n(1,j,k)
             rhoe_n(a,j,k) = rhoe_n(1,j,k)
          enddo
         
          ! set increments to 0
          Krho(1,j,k)  = 0.0_wp
          Krhou(1,j,k) = 0.0_wp
          Krhov(1,j,k) = 0.0_wp
          Krhow(1,j,k) = 0.0_wp
          Krhoe(1,j,k) = 0.0_wp

       enddo
    enddo
 
  end subroutine bc_inlet_prs_imin
  
  subroutine bc_inlet_prs_imax
  !==============================================================================
    !> boundary condition at imax for PRS eos
  !==============================================================================

    ! set total pressure and temperature
    p_tot = 250000_wp!*p_ref
    T_tot = T_ref/0.833_wp
  
    ! direction of the velocity at the boundary 
    nb = (/1.0_wp,0.0_wp,0.0_wp/)
  
    ! Initialize
    err = 1.0_wp
    tol = 1e-6_wp
    do k=1,nz
       do j=ndy_c,nfy_c
          ! Use neighbour point for test values
          rotent = rho_n(nx-1,j,k)
          ptent  = prs(nx-1,j,k)
          Ttent  = Tmp(nx-1,j,k)
          ! initialize with total values
          prs(nx,j,k) = p_tot
          Tmp(nx,j,k) = T_tot
          rho_tot = rocalc_pt(p_tot,T_tot,rotent)
          rho_n(nx,j,k) = rho_tot
       enddo
    enddo
    
    do k=1,nz
       do j=1,ny
          rhobc = rho_n(nx,j,k)
          pbc  = prs(nx,j,k)
          Tbc  = Tmp(nx,j,k)
  
          ! velocity vector at the interior nodes
          veli = (/uu(nx-1,j,k),vv(nx-1,j,k),ww(nx-1,j,k)/)
  
          ! compute total enthalpy without adiabatic condition, as H0_i = E0_i + P/rho
          ! H_tot = ecalc_pro(p_tot,rho_tot,T_tot) + prs(nx-1,j,k)/rho_n(nx-1,j,k)
  
          do while (err.gt.tol)
             
             ptmp = pbc
  
             ! exrapolate R+ from the interior nodes, so dR=0
             ! dR = u_i-u_bc + c_i/rho_i*(rho_i-rho_bc)
             ci = sqrt(c2calc_tro(Tmp(nx-1,j,k),rho_n(nx-1,j,k)))
             velbc = sqrt(dot_product(veli,veli)) + ci/rho_n(nx-1,j,k)*(rho_n(nx-1,j,k)-rhobc)
     
             ! compute total enthalpy, apply adiabaticity so H0_bc = H0_i
             H_tot = ecalc_pro(prs(nx-1,j,k),rho_n(nx-1,j,k),Ttent) + prs(nx-1,j,k)/rho_n(nx-1,j,k) &
                     + 0.5_wp*dot_product(veli,veli)
             h = H_tot - 0.5_wp*velbc**2
  
             ! update thermo variables
             Tbc  = tcalc_ph(pbc,h,rotent,Ttent)
             rhobc = rocalc_pt(pbc,Tbc,rotent)
             pbc  = pcalc_tro(Tbc,rhobc)
          
             err = (pbc-ptmp)/ptmp
          enddo
  
       ! update thermo variables at boundary 
       rho_n(nx,j,k) = rhobc
       prs(nx,j,k)   = pbc 
       Tmp(nx,j,k)   = Tbc 
  
       ! compute velocity and momentum variables
       uu(nx,j,k) = velbc*nb(1)
       vv(nx,j,k) = velbc*nb(2)
       ww(nx,j,k) = velbc*nb(3)
       rhou_n(nx,j,k) = rho_n(nx,j,k)*uu(nx,j,k)
       rhov_n(nx,j,k) = rho_n(nx,j,k)*vv(nx,j,k)
       rhow_n(nx,j,k) = rho_n(nx,j,k)*ww(nx,j,k)
       
       ! compute the rest
       visc(nx,j,k)= viscosity_law(Tmp(nx,j,k),rho_n(nx,j,k))
       cok(nx,j,k) = thconductivity(visc(nx,j,k),Tmp(nx,j,k),rho_n(nx,j,k))
       c_(nx,j,k)  = sqrt(c2calc_tro(Tmp(nx,j,k),rho_n(nx,j,k)))
       rhoe_n(nx,j,k) = rho_n(nx,j,k)*(ecalc_tro(Tmp(nx,j,k),rho_n(nx,j,k)) &
         + 0.5_wp*(uu(nx,j,k)**2+vv(nx,j,k)**2+ww(nx,j,k)**2))
       
       ! extrapolate to ghost points
       do a=nx,nx+ngh
          Tmp(a,j,k)    = Tmp(nx,j,k)
          prs(a,j,k)    = prs(nx,j,k)
          rho_n(a,j,k)  = rho_n(nx,j,k)
          rhou_n(a,j,k) = rhou_n(nx,j,k)
          rhov_n(a,j,k) = rhov_n(nx,j,k)
          rhow_n(a,j,k) = rhow_n(nx,j,k)
          rhoe_n(a,j,k) = rhoe_n(nx,j,k)
       enddo
      
       ! set increments to 0 
       Krho(nx,j,k)  = 0.0_wp
       Krhou(nx,j,k) = 0.0_wp
       Krhov(nx,j,k) = 0.0_wp
       Krhow(nx,j,k) = 0.0_wp
       Krhoe(nx,j,k) = 0.0_wp
  
       enddo
    enddo
  
  end subroutine bc_inlet_prs_imax
  
  
  subroutine bc_inlet_prs_jmin
  !==============================================================================
    !> boundary condition at jmin for PRS eos
  !==============================================================================

    ! set total pressure and temperature
    p_tot = 250000_wp!*p_ref
    T_tot = T_ref/0.833_wp
  
    ! direction of the velocity at the boundary 
    nb = (/1.0_wp,0.0_wp,0.0_wp/)
  
    ! Initialize
    err = 1.0_wp
    tol = 1e-6_wp
    do k=1,nz
       do i=ndx_c,nfx_c
          ! Use neighbour point for test values
          rotent = rho_n(i,2,k)
          ptent  = prs(i,2,k)
          Ttent  = Tmp(i,2,k)
          ! initialize with total values
          prs(i,1,k) = p_tot
          Tmp(i,1,k) = T_tot
          rho_tot = rocalc_pt(p_tot,T_tot,rotent)
          rho_n(i,1,k) = rho_tot
       enddo
    enddo
    
    do k=1,nz
       do i=ndx_c,nfx_c
          rhobc = rho_n(i,1,k)
          pbc  = prs(i,1,k)
          Tbc  = Tmp(i,1,k)
  
          ! velocity vector at the interior nodes
          veli = (/uu(i,2,k),vv(i,2,k),ww(i,2,k)/)
  
          ! compute total enthalpy without adiabatic condition, as H0_i = E0_i + P/rho
          ! H_tot = ecalc_pro(p_tot,rho_tot,T_tot) + prs(i,2,k)/rho_n(i,2,k)
  
          do while (err.gt.tol)
             
             ptmp = pbc
  
             ! exrapolate R+ from the interior nodes, so dR=0
             ! dR = u_i-u_bc + c_i/rho_i*(rho_i-rho_bc)
             ci = sqrt(c2calc_tro(Tmp(i,2,k),rho_n(i,2,k)))
             velbc = sqrt(dot_product(veli,veli)) + ci/rho_n(i,2,k)*(rho_n(i,2,k)-rhobc)
     
             ! compute total enthalpy, apply adiabaticity so H0_bc = H0_i
             H_tot = ecalc_pro(prs(i,2,k),rho_n(i,2,k),Ttent) + prs(i,2,k)/rho_n(i,2,k) &
                     + 0.5_wp*dot_product(veli,veli)
             h = H_tot - 0.5_wp*velbc**2
  
             ! update thermo variables
             Tbc  = tcalc_ph(pbc,h,rotent,Ttent)
             rhobc = rocalc_pt(pbc,Tbc,rotent)
             pbc  = pcalc_tro(Tbc,rhobc)
          
             err = (pbc-ptmp)/ptmp
          enddo
  
       ! update thermo variables at boundary 
       rho_n(i,1,k) = rhobc
       prs(i,1,k)   = pbc 
       Tmp(i,1,k)   = Tbc 
  
       ! compute velocity and momentum variables
       uu(i,1,k) = velbc*nb(1)
       vv(i,1,k) = velbc*nb(2)
       ww(i,1,k) = velbc*nb(3)
       rhou_n(i,1,k) = rho_n(i,1,k)*uu(i,1,k)
       rhov_n(i,1,k) = rho_n(i,1,k)*vv(i,1,k)
       rhow_n(i,1,k) = rho_n(i,1,k)*ww(i,1,k)
       
       ! compute the rest
       visc(i,1,k)= viscosity_law(Tmp(i,1,k),rho_n(i,1,k))
       cok(i,1,k) = thconductivity(visc(i,1,k),Tmp(i,1,k),rho_n(i,1,k))
       c_(i,1,k)  = sqrt(c2calc_tro(Tmp(i,1,k),rho_n(i,1,k)))
       rhoe_n(i,1,k) = rho_n(i,1,k)*(ecalc_tro(Tmp(i,1,k),rho_n(i,1,k)) &
         + 0.5_wp*(uu(i,1,k)**2+vv(i,1,k)**2+ww(i,1,k)**2))
       
       ! extrapolate to ghost points
       do a=1-ngh,1
          Tmp(i,a,k)    = Tmp(i,1,k)
          prs(i,a,k)    = prs(i,1,k)
          rho_n(i,a,k)  = rho_n(i,1,k)
          rhou_n(i,a,k) = rhou_n(i,1,k)
          rhov_n(i,a,k) = rhov_n(i,1,k)
          rhow_n(i,a,k) = rhow_n(i,1,k)
          rhoe_n(i,a,k) = rhoe_n(i,1,k)
       enddo
      
       ! set increments to 0
       Krho(i,1,k)  = 0.0_wp
       Krhou(i,1,k) = 0.0_wp
       Krhov(i,1,k) = 0.0_wp
       Krhow(i,1,k) = 0.0_wp
       Krhoe(i,1,k) = 0.0_wp
  
       enddo
    enddo
  
  end subroutine bc_inlet_prs_jmin
  
  
  subroutine bc_inlet_prs_jmax
  !==============================================================================
    !> boundary condition at jmax for PRS eos
  !==============================================================================

    ! set total pressure and temperature
    p_tot = 250000_wp!*p_ref
    T_tot = T_ref/0.833_wp
  
    ! direction of the velocity at the boundary 
    nb = (/1.0_wp,0.0_wp,0.0_wp/)
  
    ! Initialize
    err = 1.0_wp
    tol = 1e-6_wp
    do k=1,nz
       do i=ndx_c,nfx_c
          ! Use neighbour point for test values
          rotent = rho_n(i,ny-1,k)
          ptent  = prs(i,ny-1,k)
          Ttent  = Tmp(i,ny-1,k)
          ! initialize with total values
          prs(i,ny,k) = p_tot
          Tmp(i,ny,k) = T_tot
          rho_tot = rocalc_pt(p_tot,T_tot,rotent)
          rho_n(i,ny,k) = rho_tot
       enddo
    enddo
    
    do k=1,nz
       do i=ndx_c,nfx_c
          rhobc = rho_n(i,ny,k)
          pbc  = prs(i,ny,k)
          Tbc  = Tmp(i,ny,k)
  
          ! velocity vector at the interior nodes
          veli = (/uu(i,ny-1,k),vv(i,ny-1,k),ww(i,ny-1,k)/)
  
          ! compute total enthalpy without adiabatic condition, as H0_i = E0_i + P/rho
          ! H_tot = ecalc_pro(p_tot,rho_tot,T_tot) + prs(i,ny-1,k)/rho_n(i,ny-1,k)
  
          do while (err.gt.tol)
             
             ptmp = pbc
  
             ! exrapolate R+ from the interior nodes, so dR=0
             ! dR = u_i-u_bc + c_i/rho_i*(rho_i-rho_bc)
             ci = sqrt(c2calc_tro(Tmp(i,ny-1,k),rho_n(i,ny-1,k)))
             velbc = sqrt(dot_product(veli,veli)) + ci/rho_n(i,ny-1,k)*(rho_n(i,ny-1,k)-rhobc)
     
             ! compute total enthalpy, apply adiabaticity so H0_bc = H0_i
             H_tot = ecalc_pro(prs(i,ny-1,k),rho_n(i,ny-1,k),Ttent) + prs(i,ny-1,k)/rho_n(i,ny-1,k) &
                     + 0.5_wp*dot_product(veli,veli)
             h = H_tot - 0.5_wp*velbc**2
  
             ! update thermo variables
             Tbc  = tcalc_ph(pbc,h,rotent,Ttent)
             rhobc = rocalc_pt(pbc,Tbc,rotent)
             pbc  = pcalc_tro(Tbc,rhobc)
          
             err = (pbc-ptmp)/ptmp
          enddo
  
       ! update thermo variables at boundary 
       rho_n(i,ny,k) = rhobc
       prs(i,ny,k)   = pbc 
       Tmp(i,ny,k)   = Tbc 
  
       ! compute velocity and momentum variables
       uu(i,ny,k) = velbc*nb(1)
       vv(i,ny,k) = velbc*nb(2)
       ww(i,ny,k) = velbc*nb(3)
       rhou_n(i,ny,k) = rho_n(i,ny,k)*uu(i,ny,k)
       rhov_n(i,ny,k) = rho_n(i,ny,k)*vv(i,ny,k)
       rhow_n(i,ny,k) = rho_n(i,ny,k)*ww(i,ny,k)
       
       ! compute the rest
       visc(i,ny,k)= viscosity_law(Tmp(i,ny,k),rho_n(i,ny,k))
       cok(i,ny,k) = thconductivity(visc(i,ny,k),Tmp(i,ny,k),rho_n(i,ny,k))
       c_(i,ny,k)  = sqrt(c2calc_tro(Tmp(i,ny,k),rho_n(i,ny,k)))
       rhoe_n(i,ny,k) = rho_n(i,ny,k)*(ecalc_tro(Tmp(i,ny,k),rho_n(i,ny,k)) &
         + 0.5_wp*(uu(i,ny,k)**2+vv(i,ny,k)**2+ww(i,ny,k)**2))
       
       ! extrapolate to ghost points
       do a=ny,ny+ngh
          Tmp(i,a,k)    = Tmp(i,ny,k)
          prs(i,a,k)    = prs(i,ny,k)
          rho_n(i,a,k)  = rho_n(i,ny,k)
          rhou_n(i,a,k) = rhou_n(i,ny,k)
          rhov_n(i,a,k) = rhov_n(i,ny,k)
          rhow_n(i,a,k) = rhow_n(i,ny,k)
          rhoe_n(i,a,k) = rhoe_n(i,ny,k)
       enddo
      
       ! set increments to 0
       Krho(i,ny,k)  = 0.0_wp
       Krhou(i,ny,k) = 0.0_wp
       Krhov(i,ny,k) = 0.0_wp
       Krhow(i,ny,k) = 0.0_wp
       Krhoe(i,ny,k) = 0.0_wp
  
       enddo
    enddo
  
  end subroutine bc_inlet_prs_jmax
  
  ! ======= !
  ! Outlets !
  ! ======= !

  ! PFG
  ! === 
  
  subroutine bc_backpressure_imax
    use mod_time
    
    ! Update outlet points by applying back-pressure (ONLY FOR SUBSONIC OUTFLOW)
    ! ==================
    machbc= 0.932
    theta = -67.8_wp
    theta = theta*pi/180.0_wp
    !nb = (/cos(theta),sin(theta),0.0_wp/)
    nb = (/1.0_wp,0.0_wp,0.0_wp/)
    !p_tot=129000.0_wp
    p_tot=p_ref
    
    do k=1,nz
       !do j=ndy_c,nfy_c
       do j=1,ny
       
          vel_int = (/uu(nx-1,j,k),vv(nx-1,j,k),ww(nx-1,j,k)/)
          !machloc = sqrt(dot_product(vel_int,vel_int))/c_(nx-1,j,k) 
          machloc = dot_product(vel_int,nb)/c_(nx-1,j,k) ! normal to boundary

          !!!! Impose Mach number
          !! extrapolate temperature from the interior nodes
          !Tmp(nx,j,k) = Tmp(nx-1,j,k)

          !! extrapolate speed of sound
          !c_(nx,j,k) = sqrt(gam*Rg*Tmp(nx,j,k))

          !! new total pressure
          !p_tot = prs(nx-1,j,k)*(1.0_wp+gam1/2.0_wp*machloc**2)**(gam/gam1)
          !prs(nx,j,k) = p_tot/(1.0_wp+gam1/2.0_wp*machbc**2)**(gam/gam1)

          !! density
          !rho_n(nx,j,k) = rocalc_pt(prs(nx,j,k),Tmp(nx,j,k),rho_n(nx,j,k))

          !!! Impose back pressure 
          ! set pressure as back pressure (p_ref) when M_i<1
          ! if M_i>=1, pressure is extrapolated from the interior node
          if(machloc.lt.1.0_wp) then
             !prs(nx,j,k) = 10943.65_wp
             !prs(nx,j,k) = p_tot/1.82_wp
             !prs(nx,j,k) = p_ref*0.528281787717174_wp
             prs(nx,j,k) = p_ref*0.534460302143549_wp
          else
             prs(nx,j,k) = prs(nx-1,j,k)
          endif

          !if ((mod(ntotal,100)==0).and.(irk==nrk).and.(j==25).and.(iproc>=32)) print *,'iproc',iproc,'p_out/p_ref',prs(nx,j,1)/p_ref
       
          ! extrapolate entropy from the interior nodes
          s = scalc_tro(Tmp(nx-1,j,k),rho_n(nx-1,j,k))

          ! calc rho and Tmp
          ro1 = rocalc_ps(prs(nx,j,k),s,Tmp(nx-1,j,k))
          rho_n(nx,j,k) = ro1
          Tmp(nx,j,k) = tcalc_pro(prs(nx,j,k),rho_n(nx,j,k),Tmp(nx-1,j,k))
             
          ! extrapolate velocities from the interior nodes
          uu(nx,j,k) = uu(nx-1,j,k)
          vv(nx,j,k) = vv(nx-1,j,k)
          ww(nx,j,k) = ww(nx-1,j,k)
          rhou_n(nx,j,k) = rho_n(nx,j,k)*uu(nx,j,k)
          rhov_n(nx,j,k) = rho_n(nx,j,k)*vv(nx,j,k)
          rhow_n(nx,j,k) = rho_n(nx,j,k)*ww(nx,j,k)
          
          ! compute the rest
          !visc(nx,j,k)=viscosity_law(Tmp(nx,j,k),ro1)*diffscale
          !cok(nx,j,k) =thconductivity(visc(nx,j,k),Tmp(nx,j,k),ro1)
          !c_(nx,j,k)  = sqrt(c2calc_tro(Tmp(nx,j,k),ro1))
          rhoe_n(nx,j,k) = rho_n(nx,j,k)*(ecalc_tro(Tmp(nx,j,k),ro1) &
                         + 0.5_wp*(uu(nx,j,k)**2+vv(nx,j,k)**2+ww(nx,j,k)**2))

          ! extrapolate to ghost points
          !do a=nx,nx+ngh
          !   Tmp(a,j,k)    = Tmp(nx,j,k)
          !   prs(a,j,k)    = prs(nx,j,k)
          !   rho_n(a,j,k)  = rho_n(nx,j,k)
          !   rhou_n(a,j,k) = rhou_n(nx,j,k)
          !   rhov_n(a,j,k) = rhov_n(nx,j,k)
          !   rhow_n(a,j,k) = rhow_n(nx,j,k)
          !   rhoe_n(a,j,k) = rhoe_n(nx,j,k)
          !enddo
         
          ! set increments to 0 
          Krho(nx,j,k) = 0.0_wp
          Krhou(nx,j,k) = 0.0_wp
          Krhov(nx,j,k) = 0.0_wp
          Krhow(nx,j,k) = 0.0_wp
          Krhoe(nx,j,k) = 0.0_wp
          
       enddo
    enddo
  
!!$    ! Update outlet points by applying back-pressure (ONLY FOR SUBSONIC OUTFLOW)
!!$    ! ==================
!!$    nb = (/1.0_wp,0.0_wp,0.0_wp/)
!!$
!!$    do k=1,nz
!!$       do j=1,ny
!!$
!!$          ! if M_i>=1, pressure is extrapolated from the interior node
!!$          vel_int = (/uu(nx-1,j,k),vv(nx-1,j,k),ww(nx-1,j,k)/)
!!$          !machloc = sqrt(dot_product(vel_int,vel_int))/c_(nx-1,j,k)
!!$          machloc = dot_product(vel_int,nb)/c_(nx-1,j,k) ! normal to boundary
!!$                    
!!$          ! set pressure as back pressure (p_ref) when M_i<1
!!$          if (machloc.lt.1.0_wp) then
!!$             prs(nx,j,k) = p_ref
!!$             !prs(nx,j,k) = 0.92*p_ref
!!$             !prs(nx,j,k) = 0.528*p_ref
!!$          else
!!$             prs(nx,j,k) = prs(nx-1,j,k)
!!$          endif
!!$
!!$          ! extrapolate entropy from the interior nodes
!!$          s = scalc_tro(Tmp(nx-1,j,k),rho_n(nx-1,j,k))
!!$
!!$          ! calc rho and Tmp
!!$          ro1 = rocalc_ps(prs(nx,j,k),s,Tmp(nx-1,j,k))
!!$          rho_n(nx,j,k) = ro1
!!$          Tmp(nx,j,k) = tcalc_pro(prs(nx,j,k),rho_n(nx,j,k),Tmp(nx-1,j,k))
!!$             
!!$          ! extrapolate velocities from the interior nodes
!!$          uu(nx,j,k) = uu(nx-1,j,k)!2.0_wp*rhou_n(nx-1,j,k)/rho_n(nx-1,j,k)-rhou_n(nx-2,j,k)/rho_n(nx-2,j,k)
!!$          vv(nx,j,k) = vv(nx-1,j,k)!2.0_wp*rhov_n(nx-1,j,k)/rho_n(nx-1,j,k)-rhov_n(nx-2,j,k)/rho_n(nx-2,j,k)
!!$          ww(nx,j,k) = ww(nx-1,j,k)!2.0_wp*rhow_n(nx-1,j,k)/rho_n(nx-1,j,k)-rhow_n(nx-2,j,k)/rho_n(nx-2,j,k)
!!$          rhou_n(nx,j,k) = ro1*uu(nx,j,k)
!!$          rhov_n(nx,j,k) = ro1*vv(nx,j,k)
!!$          rhow_n(nx,j,k) = ro1*ww(nx,j,k)
!!$
!!$          ! compute the rest
!!$          !visc(nx,j,k)=viscosity_law(Tmp(nx,j,k),ro1)*diffscale
!!$          !cok(nx,j,k) =thconductivity(visc(nx,j,k),Tmp(nx,j,k),ro1)
!!$          !c_(nx,j,k)  = sqrt(c2calc_tro(Tmp(nx,j,k),ro1))
!!$          rhoe_n(nx,j,k) = ro1*(ecalc_tro(Tmp(nx,j,k),ro1) &
!!$               + 0.5_wp*(uu(nx,j,k)**2+vv(nx,j,k)**2+ww(nx,j,k)**2))
!!$          
!!$          Krho(nx,j,k) = 0.0_wp
!!$          Krhou(nx,j,k) = 0.0_wp
!!$          Krhov(nx,j,k) = 0.0_wp
!!$          Krhow(nx,j,k) = 0.0_wp
!!$          Krhoe(nx,j,k) = 0.0_wp
!!$
!!$       enddo
!!$    enddo

  end subroutine bc_backpressure_imax

  subroutine bc_backpressure_jmax
    
    ! Update outlet points by applying back-pressure (ONLY FOR SUBSONIC OUTFLOW)
    ! ==================
    nb = (/1.0_wp,0.0_wp,0.0_wp/)
    
    do k=1,nz
       do i=1,nx
       
          vel_int = (/uu(i,ny-1,k),vv(i,ny-1,k),ww(i,ny-1,k)/)
          !machloc = sqrt(dot_product(vel_int,vel_int))/c_(i,ny-1,k) 
          machloc = dot_product(vel_int,nb)/c_(i,ny-1,k) ! normal to boundary

          ! set pressure as back pressure (p_ref) when M_i<1
          ! if M_i>=1, pressure is extrapolated from the interior node
          if(machloc.lt.1.0_wp) then
            !prs(i,ny,k) = 10943.65_wp
            !prs(i,ny,k) = p_tot/1.82_wp
            prs(i,ny,k) = p_ref
          else
            prs(i,ny,k) = prs(i,ny-1,k)
          endif
       
          ! extrapolate entropy from the interior nodes
          s = scalc_tro(Tmp(i,ny-1,k),rho_n(i,ny-1,k))

          ! calc rho and Tmp
          ro1 = rocalc_ps(prs(i,ny,k),s,Tmp(i,ny-1,k))
          rho_n(i,ny,k) = ro1
          Tmp(i,ny,k) = tcalc_pro(prs(i,ny,k),rho_n(i,ny,k),Tmp(i,ny-1,k))
             
          ! extrapolate velocities from the interior nodes
          uu(i,ny,k) = uu(i,ny-1,k)
          vv(i,ny,k) = vv(i,ny-1,k)
          ww(i,ny,k) = ww(i,ny-1,k)
          rhou_n(i,ny,k) = rho_n(i,ny,k)*uu(i,ny,k)
          rhov_n(i,ny,k) = rho_n(i,ny,k)*vv(i,ny,k)
          rhow_n(i,ny,k) = rho_n(i,ny,k)*ww(i,ny,k)
          
          ! compute the rest
          !visc(i,ny,k)=viscosity_law(Tmp(i,ny,k),ro1)*diffscale
          !cok(i,ny,k) =thconductivity(visc(i,ny,k),Tmp(i,ny,k),ro1)
          !c_(i,ny,k)  = sqrt(c2calc_tro(Tmp(i,ny,k),ro1))
          rhoe_n(i,ny,k) = rho_n(i,ny,k)*(ecalc_tro(Tmp(i,ny,k),ro1) &
                         + 0.5_wp*(uu(i,ny,k)**2+vv(i,ny,k)**2+ww(i,ny,k)**2))

          ! extrapolate to ghost points
          !do a=nx,nx+ngh
          !   Tmp(a,j,k)    = Tmp(i,ny,k)
          !   prs(a,j,k)    = prs(i,ny,k)
          !   rho_n(a,j,k)  = rho_n(i,ny,k)
          !   rhou_n(a,j,k) = rhou_n(i,ny,k)
          !   rhov_n(a,j,k) = rhov_n(i,ny,k)
          !   rhow_n(a,j,k) = rhow_n(i,ny,k)
          !   rhoe_n(a,j,k) = rhoe_n(i,ny,k)
          !enddo
         
          ! set increments to 0 
          Krho(i,ny,k) = 0.0_wp
          Krhou(i,ny,k) = 0.0_wp
          Krhov(i,ny,k) = 0.0_wp
          Krhow(i,ny,k) = 0.0_wp
          Krhoe(i,ny,k) = 0.0_wp
          
       enddo
    enddo

  end subroutine bc_backpressure_jmax

!!$  subroutine bc_backpressure_jmax
!!$    
!!$    ! Update outlet points by applying back-pressure (ONLY FOR SUBSONIC OUTFLOW)
!!$    ! ==================
!!$    do k=1,nz
!!$       do i=ndx_c,nfx_c
!!$       
!!$          ! set pressure as back pressure (p_ref) when M_i<1
!!$          ! if M_i>=1, pressure is extrapolated from the interior node
!!$          vel_int = (/uu(i,ny-1,k),vv(i,ny-1,k),ww(i,ny-1,k)/)
!!$          machloc = sqrt(dot_product(vel_int,vel_int))/c_(i,ny-1,k)
!!$          if(machloc.lt.1.0_wp) then
!!$            prs(i,ny,k) = p_ref
!!$          else
!!$            prs(i,ny,k) = prs(i,ny-1,k)
!!$          endif
!!$       
!!$          ! extrapolate temperature from the interior nodes
!!$          Tmp(i,ny,k) = Tmp(i,ny-1,k)!2.0_wp*Tmp(i,ny-1,k)-Tmp(i,ny-2,k)
!!$          
!!$          ! calc rho
!!$          ro1 = rocalc_pt(prs(i,ny,k),Tmp(i,ny,k),rho_n(i,ny-1,k) )
!!$  
!!$          rho_n(i,ny,k) = ro1
!!$             
!!$          ! extrapolate velocities from the interior nodes
!!$          uu(i,ny,k) = uu(i,ny-1,k)
!!$          vv(i,ny,k) = vv(i,ny-1,k)
!!$          ww(i,ny,k) = ww(i,ny-1,k)
!!$          rhou_n(i,ny,k) = ro1*uu(i,ny,k)
!!$          rhov_n(i,ny,k) = ro1*vv(i,ny,k)
!!$          rhow_n(i,ny,k) = ro1*ww(i,ny,k)
!!$          
!!$          ! compute the rest
!!$          visc(i,ny,k)=viscosity_law(Tmp(i,ny,k),ro1)*diffscale
!!$          cok(i,ny,k) =thconductivity(visc(i,ny,k),Tmp(i,ny,k),ro1)
!!$          c_(i,ny,k)  = sqrt(c2calc_tro(Tmp(i,ny,k),ro1))
!!$          rhoe_n(i,ny,k) = rho_n(i,ny,k)*(ecalc_tro(Tmp(i,ny,k),ro1) &
!!$            + 0.5_wp*(uu(i,ny,k)**2+vv(i,ny,k)**2+ww(i,ny,k)**2))
!!$          
!!$
!!$          ! extrapolate to ghost points
!!$          do a=ny,ny+ngh
!!$             Tmp(i,a,k)    = Tmp(i,ny,k)
!!$             prs(i,a,k)    = prs(i,ny,k)
!!$             rho_n(i,a,k)  = rho_n(i,ny,k)
!!$             rhou_n(i,a,k) = rhou_n(i,ny,k)
!!$             rhov_n(i,a,k) = rhov_n(i,ny,k)
!!$             rhow_n(i,a,k) = rhow_n(i,ny,k)
!!$             rhoe_n(i,a,k) = rhoe_n(i,ny,k)
!!$          enddo
!!$         
!!$          ! set increments to 0
!!$          Krho(i,ny,k) = 0.0_wp
!!$          Krhou(i,ny,k) = 0.0_wp
!!$          Krhov(i,ny,k) = 0.0_wp
!!$          Krhow(i,ny,k) = 0.0_wp
!!$          Krhoe(i,ny,k) = 0.0_wp
!!$          
!!$       enddo
!!$    enddo
!!$  
!!$  end subroutine bc_backpressure_jmax



  ! PRS
  ! ===

  subroutine bc_backpressure_prs_imin
  !==============================================================================
    !> boundary condition at imin for PRS eos
  !==============================================================================
    
    ! Update outlet points by applying back-pressure 
    ! ==================
    do k=1,nz
       do j=ndy_c,nfy_c
  
          Ttent = Tmp(2,j,k)
       
          ! set pressure as back pressure (p_ref) when M_i<1
          ! if M_i>=1, pressure is extrapolated from the interior node
          veli = (/uu(2,j,k),vv(2,j,k),ww(2,j,k)/)
          machloc = sqrt(dot_product(veli,veli))/c_(2,j,k)
          if(machloc.lt.1.0_wp) then
            prs(1,j,k) = 100000_wp!p_ref
          else
            prs(1,j,k) = prs(2,j,k)
          endif
       
          ! extrapolate entropy from the interior nodes
          s = scalc_tro(Tmp(2,j,k),rho_n(2,j,k))
          
          ! compute density
          rho_n(1,j,k) = rocalc_ps(prs(1,j,k),s,Ttent)
          
          ! extrapolate velocities from the interior nodes
          uu(1,j,k) = uu(2,j,k)
          vv(1,j,k) = vv(2,j,k)
          ww(1,j,k) = ww(2,j,k)
          rhou_n(1,j,k) = rho_n(1,j,k)*uu(1,j,k)
          rhov_n(1,j,k) = rho_n(1,j,k)*vv(1,j,k)
          rhow_n(1,j,k) = rho_n(1,j,k)*ww(1,j,k)
          
          ! compute the rest
          visc(i,ny,k)= viscosity_law(Tmp(i,ny,k),rho_n(i,ny,k))*diffscale
          cok(i,ny,k) = thconductivity(visc(i,ny,k),Tmp(i,ny,k),rho_n(i,ny,k))
          c_(i,ny,k)  = sqrt(c2calc_tro(Tmp(i,ny,k),rho_n(i,ny,k)))
          rhoe_n(1,j,k) = rho_n(1,j,k)*(ecalc_tro(Tmp(1,j,k),rho_n(1,j,k)) &
            + 0.5_wp*(uu(1,j,k)**2+vv(1,j,k)**2+ww(1,j,k)**2))
          
          ! extrapolate to ghost points
          do a=1-ngh,1
             Tmp(a,j,k)    = Tmp(1,j,k)
             prs(a,j,k)    = prs(1,j,k)
             rho_n(a,j,k)  = rho_n(1,j,k)
             rhou_n(a,j,k) = rhou_n(1,j,k)
             rhov_n(a,j,k) = rhov_n(1,j,k)
             rhow_n(a,j,k) = rhow_n(1,j,k)
             rhoe_n(a,j,k) = rhoe_n(1,j,k)
          enddo
         
          ! set increments to 0
          Krho(i,ny,k)  = 0.0_wp
          Krhou(i,ny,k) = 0.0_wp
          Krhov(i,ny,k) = 0.0_wp
          Krhow(i,ny,k) = 0.0_wp
          Krhoe(i,ny,k) = 0.0_wp
          
       enddo
    enddo
  
  end subroutine bc_backpressure_prs_imin

  subroutine bc_backpressure_prs_imax
  !==============================================================================
    !> boundary condition at imax for PRS eos
  !==============================================================================
    
    ! Update outlet points by applying back-pressure 
    ! ==================
    p_tot = 129000.0_wp ! (inlet) Gori et al. discharge 2 experiment D2
    nb = (/1.0_wp,0.0_wp,0.0_wp/)
    do k=1,nz
       !do j=ndy_c,nfy_c
       do j=1,ny
  
          Ttent = Tmp(nx-1,j,k)
          rotent= rho_n(nx-1,j,k)

          ! set pressure as back pressure (p_ref) when M_i<1
          ! if M_i>=1, pressure is extrapolated from the interior node
          veli = (/uu(nx-1,j,k),vv(nx-1,j,k),ww(nx-1,j,k)/)
          machloc = abs(dot_product(veli,nb))/c_(nx-1,j,k) ! Normal to boundary
          !machloc = dot_product(veli,veli)/c_(nx-1,j,k) ! magnitude
          if(machloc.lt.1.0_wp) then
             !prs(nx,j,k) = 14500.0_wp ! Gori et al. discharge 2 experiment E2
             !prs(nx,j,k) = 38800.0_wp ! Gori et al. discharge 2 experiment D2
             prs(nx,j,k) =  p_tot/1.82 ! PR = 1.82 (Cinnella LS59)

             ! extrapolate entropy from the interior nodes
             s = scalc_tro(Tmp(nx-1,j,k),rho_n(nx-1,j,k))
              
             ! compute density and temperature
             rhobc = rocalc_ps(prs(nx,j,k),s,Ttent)
             rho_n(nx,j,k) = rhobc
             !Tmp(nx,j,k) = tcalc_pro(prs(nx,j,k),rho_n(nx,j,k),Ttent)
             Tmp(nx,j,k) = tcalc_sro(s,rho_n(nx,j,k),Ttent)
             !rho_n(nx,j,k) = rocalc_pt(prs(nx,j,k),Tmp(nx,j,k),rotent)

             ! exrapolate R- from the interior nodes, so dR=0
             ! dR = u_i-u_bc - c_i/rho_i*(rho_i-rho_bc)
             ci = sqrt(c2calc_tro(Tmp(nx-1,j,k),rho_n(nx-1,j,k)))
             velbc = dot_product(veli,nb) - ci/rho_n(nx-1,j,k)*(rho_n(nx-1,j,k)-rhobc)

             ! compute velocities
             uu(nx,j,k) = veli(1)+velbc*nb(1)-dot_product(veli,nb)*nb(1)
             vv(nx,j,k) = veli(2)+velbc*nb(2)-dot_product(veli,nb)*nb(2)
             ww(nx,j,k) = veli(3)+velbc*nb(3)-dot_product(veli,nb)*nb(3)
             !print*,sqrt(dot_product(veli,veli)),velbc
             !print*
          else
             prs(nx,j,k) = prs(nx-1,j,k)
             rho_n(nx,j,k) = rho_n(nx-1,j,k)
             Tmp(nx,j,k) = Tmp(nx-1,j,k)

             ! extrapolate velocities from the interior nodes
             uu(nx,j,k) = uu(nx-1,j,k)
             vv(nx,j,k) = vv(nx-1,j,k)
             ww(nx,j,k) = ww(nx-1,j,k)
          endif
          
          ! extrapolate velocities from the interior nodes
          uu(nx,j,k) = uu(nx-1,j,k)
          vv(nx,j,k) = vv(nx-1,j,k)
          ww(nx,j,k) = ww(nx-1,j,k)
          ! Compute conservative variables
          rhou_n(nx,j,k) = rho_n(nx,j,k)*uu(nx,j,k)
          rhov_n(nx,j,k) = rho_n(nx,j,k)*vv(nx,j,k)
          rhow_n(nx,j,k) = rho_n(nx,j,k)*ww(nx,j,k)
          
          ! compute the rest
          visc(nx,j,k)= viscosity_law(Tmp(nx,j,k),rho_n(nx,j,k))*diffscale
          cok(nx,j,k) = thconductivity(visc(nx,j,k),Tmp(nx,j,k),rho_n(nx,j,k))
          c_(nx,j,k)  = sqrt(c2calc_tro(Tmp(nx,j,k),rho_n(nx,j,k)))
          rhoe_n(nx,j,k) = rho_n(nx,j,k)*(ecalc_tro(Tmp(nx,j,k),rho_n(nx,j,k)) &
            + 0.5_wp*(uu(nx,j,k)**2+vv(nx,j,k)**2+ww(nx,j,k)**2))
                    
          ! extrapolate to ghost points
          do a=nx,nx+ngh
             Tmp(a,j,k)    = Tmp(nx,j,k)
             prs(a,j,k)    = prs(nx,j,k)
             rho_n(a,j,k)  = rho_n(nx,j,k)
             rhou_n(a,j,k) = rhou_n(nx,j,k)
             rhov_n(a,j,k) = rhov_n(nx,j,k)
             rhow_n(a,j,k) = rhow_n(nx,j,k)
             rhoe_n(a,j,k) = rhoe_n(nx,j,k)
          enddo
         
          ! set increments to 0 
          Krho(nx,j,k)  = 0.0_wp
          Krhou(nx,j,k) = 0.0_wp
          Krhov(nx,j,k) = 0.0_wp
          Krhow(nx,j,k) = 0.0_wp
          Krhoe(nx,j,k) = 0.0_wp
          
       enddo
    enddo
  
  end subroutine bc_backpressure_prs_imax

  subroutine bc_backpressure_prs_jmin
  !==============================================================================
    !> boundary condition at jmin for PRS eos
  !==============================================================================
    
    ! Update outlet points by applying back-pressure 
    ! ==================
    do k=1,nz
       do i=ndx_c,nfx_c
  
          Ttent = Tmp(i,2,k)
       
          ! set pressure as back pressure (p_ref) when M_i<1
          ! if M_i>=1, pressure is extrapolated from the interior node
          veli = (/uu(i,2,k),vv(i,2,k),ww(i,2,k)/)
          machloc = sqrt(dot_product(veli,veli))/c_(i,2,k)
          if(machloc.lt.1.0_wp) then
            prs(i,1,k) = 100000_wp!p_ref
          else
            prs(i,1,k) = prs(i,2,k)
          endif
       
          ! extrapolate entropy from the interior nodes
          s = scalc_tro(Tmp(i,2,k),rho_n(i,2,k))
          
          ! compute density
          rho_n(i,1,k) = rocalc_ps(prs(i,1,k),s,Ttent)
          
          ! extrapolate velocities from the interior nodes
          uu(i,1,k) = uu(i,2,k)
          vv(i,1,k) = vv(i,2,k)
          ww(i,1,k) = ww(i,2,k)
          rhou_n(i,1,k) = rho_n(i,1,k)*uu(i,1,k)
          rhov_n(i,1,k) = rho_n(i,1,k)*vv(i,1,k)
          rhow_n(i,1,k) = rho_n(i,1,k)*ww(i,1,k)
          
          ! compute the rest
          visc(i,1,k)= viscosity_law(Tmp(i,1,k),rho_n(i,1,k))*diffscale
          cok(i,1,k) = thconductivity(visc(i,1,k),Tmp(i,1,k),rho_n(i,1,k))
          c_(i,1,k)  = sqrt(c2calc_tro(Tmp(i,1,k),rho_n(i,1,k)))
          rhoe_n(i,1,k) = rho_n(i,1,k)*(ecalc_tro(Tmp(i,1,k),rho_n(i,1,k)) &
            + 0.5_wp*(uu(i,1,k)**2+vv(i,1,k)**2+ww(i,1,k)**2))
          
          ! extrapolate to ghost points
          do a=1-ngh,1
             Tmp(i,a,k)    = Tmp(i,1,k)
             prs(i,a,k)    = prs(i,1,k)
             rho_n(i,a,k)  = rho_n(i,1,k)
             rhou_n(i,a,k) = rhou_n(i,1,k)
             rhov_n(i,a,k) = rhov_n(i,1,k)
             rhow_n(i,a,k) = rhow_n(i,1,k)
             rhoe_n(i,a,k) = rhoe_n(i,1,k)
          enddo
         
          ! set increments to 0
          Krho(i,1,k)  = 0.0_wp
          Krhou(i,1,k) = 0.0_wp
          Krhov(i,1,k) = 0.0_wp
          Krhow(i,1,k) = 0.0_wp
          Krhoe(i,1,k) = 0.0_wp
          
       enddo
    enddo
  
  end subroutine bc_backpressure_prs_jmin

  subroutine bc_backpressure_prs_jmax
  !==============================================================================
    !> boundary condition at jmax for PRS eos
  !==============================================================================
    
    ! Update outlet points by applying back-pressure 
    ! ==================
    do k=1,nz
       do i=ndx_c,nfx_c
  
          Ttent = Tmp(i,ny-1,k)
       
          ! set pressure as back pressure (p_ref) when M_i<1
          ! if M_i>=1, pressure is extrapolated from the interior node
          veli = (/uu(i,ny-1,k),vv(i,ny-1,k),ww(i,ny-1,k)/)
          machloc = sqrt(dot_product(veli,veli))/c_(i,ny-1,k)
          if(machloc.lt.1.0_wp) then
            prs(i,ny,k) = 100000_wp!p_ref
          else
            prs(i,ny,k) = prs(i,ny-1,k)
          endif
       
          ! extrapolate entropy from the interior nodes
          s = scalc_tro(Tmp(i,ny-1,k),rho_n(i,ny-1,k))
          
          ! compute density
          rho_n(i,ny,k) = rocalc_ps(prs(i,ny,k),s,Ttent)
          
          ! extrapolate velocities from the interior nodes
          uu(i,ny,k) = uu(i,ny-1,k)
          vv(i,ny,k) = vv(i,ny-1,k)
          ww(i,ny,k) = ww(i,ny-1,k)
          rhou_n(i,ny,k) = rho_n(i,ny,k)*uu(i,ny,k)
          rhov_n(i,ny,k) = rho_n(i,ny,k)*vv(i,ny,k)
          rhow_n(i,ny,k) = rho_n(i,ny,k)*ww(i,ny,k)
          
          ! compute the rest
          visc(i,ny,k)= viscosity_law(Tmp(i,ny,k),rho_n(i,ny,k))*diffscale
          cok(i,ny,k) = thconductivity(visc(i,ny,k),Tmp(i,ny,k),rho_n(i,ny,k))
          c_(i,ny,k)  = sqrt(c2calc_tro(Tmp(i,ny,k),rho_n(i,ny,k)))
          rhoe_n(i,ny,k) = rho_n(i,ny,k)*(ecalc_tro(Tmp(i,ny,k),rho_n(i,ny,k)) &
            + 0.5_wp*(uu(i,ny,k)**2+vv(i,ny,k)**2+ww(i,ny,k)**2))
          
          ! extrapolate to ghost points
          do a=ny,ny+ngh
             Tmp(i,a,k)    = Tmp(i,ny,k)
             prs(i,a,k)    = prs(i,ny,k)
             rho_n(i,a,k)  = rho_n(i,ny,k)
             rhou_n(i,a,k) = rhou_n(i,ny,k)
             rhov_n(i,a,k) = rhov_n(i,ny,k)
             rhow_n(i,a,k) = rhow_n(i,ny,k)
             rhoe_n(i,a,k) = rhoe_n(i,ny,k)
          enddo
         
          ! set increments to 0
          Krho(i,ny,k)  = 0.0_wp
          Krhou(i,ny,k) = 0.0_wp
          Krhov(i,ny,k) = 0.0_wp
          Krhow(i,ny,k) = 0.0_wp
          Krhoe(i,ny,k) = 0.0_wp
          
       enddo
    enddo
  
  end subroutine bc_backpressure_prs_jmax

end module mod_turb_inlet_outlet
