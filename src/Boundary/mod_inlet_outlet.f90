!===========================================================================
module mod_inlet_outlet
!===========================================================================
  !> authors: CM
  !> date: 09/11/2022
  !> Subsonic inlet and oulet boundary conditions
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
  integer  :: i,j,k,a
  real(wp) :: Rp,Rm,s
  real(wp) :: nb(3),nbc(3),veltmp(3) ! direction, BC normal and velocity of flow at inlet
  real(wp) :: ptent,Ttent,rotent  ! tentative variables for PRS Newton iterations
  real(wp) :: pbc,Tbc,rhobc ! thermos at boundary
  real(wp) :: velbc,vel_ref(3),veli(3) ! velocities at boundary and inside domain
  real(wp) :: err,tol ! convergence criterions
  real(wp) :: h,H_tot ! enthalpy
  real(wp) :: ptmp ! pressure to loop over
  real(wp) :: ucbc,vcbc
  real(wp) :: uci,vci
  real(wp) :: uc_ref,vc_ref
  real(wp) :: theta ! inlet angle
  ! ---------------------------------------------------------------------------

contains

  ! ====== !
  ! Inlets !
  ! ====== !

  ! PFG
  ! ===

  subroutine bc_inlet_pfg_imin
  !==============================================================================
    !> Inlet boundary condition at imin for PFG eos
    !> Set free stream static quantities
  !==============================================================================

    ! Subscripts:
    ! ref: free stream values, prescribed in param.ini
    ! i  : interior of domain, so neighbour node
    ! bc : boundary value

    ! pressure from free stream
    p_ref = pcalc_tro(T_ref,rho_ref)
  
    ! direction of the velocity at the boundary
    theta = 30.0_wp
    theta = theta*pi/180.0_wp
    nb = (/cos(theta),sin(theta),0.0_wp/)
    nb = (/1.0_wp,0.0_wp,0.0_wp/)

    ! velocity vector at free stream (from inlet Mach)
    c_ref = sqrt(c2calc_tro(T_ref,rho_ref))
    vel_ref = c_ref*Mach*nb
  
    i=1
    ! Update inlet points
    ! ===================  
    do k=1,nz
       do j=ndy_e,nfy_e

          ! boundary normals
          nbc = (/-nxn_imin(j,k),-nyn_imin(j,k),nzn_imin(j,k)/)
          nbc = (/-1.0_wp,0.0_wp,0.0_wp/)
       
          ! velocity vector at the interior nodes
          veli = (/uu(i+1,j,k),vv(i+1,j,k),ww(i+1,j,k)/)

          ! Inlet
          if (dot_product(veli,nbc).le.0.0_wp) then
             ! extrapolate R+ from the interior nodes
             Rp = dot_product(veli,nbc)+2.0_wp*c_(i+1,j,k)/gam1

             ! extrapolate R- from the prescribed free stream
             Rm = dot_product(vel_ref,nbc)-2.0_wp*c_ref/gam1

             ! Compute velocity and speed of sound at boundary
             velbc = (Rp+Rm)/2.0_wp
             c_(i,j,k) = gam1*(Rp-Rm)/4.0_wp

             ! Velocity components
             uu(i,j,k) = vel_ref(1)+velbc*nbc(1)-dot_product(vel_ref,nbc)*nbc(1)
             vv(i,j,k) = vel_ref(2)+velbc*nbc(2)-dot_product(vel_ref,nbc)*nbc(2)
             ww(i,j,k) = vel_ref(3)+velbc*nbc(3)-dot_product(vel_ref,nbc)*nbc(3)

             ! Extrapolate entropy from free stream: sbc = s_ref
             s = scalc_tro(T_ref,rho_ref)
          ! Outlet
          else 
             ! exrapolate R- from the interior nodes
             Rm = dot_product(veli,nbc)-2.0_wp*c_(i+1,j,k)/gam1

             ! exrapolate R+ from the prescribed free stream
             Rp = dot_product(vel_ref,nbc)+2.0_wp*c_ref/gam1

             ! Compute velocity and speed of sound at boundary
             velbc = (Rp+Rm)/2.0_wp
             c_(i,j,k) = gam1*(Rp-Rm)/4.0_wp

             ! Velocity components
             uu(i,j,k) = veli(1)+velbc*nbc(1)-dot_product(veli,nbc)*nbc(1)
             vv(i,j,k) = veli(2)+velbc*nbc(2)-dot_product(veli,nbc)*nbc(2)
             ww(i,j,k) = veli(3)+velbc*nbc(3)-dot_product(veli,nbc)*nbc(3)

             ! Extrapolate entropy from interior: sbc = si
             s = scalc_tro(Tmp(i+1,j,k),rho_n(i+1,j,k))
          endif

          ! Compute density, pressure and temperature based on s and c
          rho_n(i,j,k) = (c_(i,j,k)**2/(gam*s))**(1.0_wp/gam1)
          prs(i,j,k) = rho_n(i,j,k)*c_(i,j,k)**2/gam
          Tmp(i,j,k) = tcalc_pro(prs(i,j,k),rho_n(i,j,k),Tmp(i+1,j,k))

          ! Conservative variables
          rhou_n(i,j,k) = rho_n(i,j,k)*uu(i,j,k)
          rhov_n(i,j,k) = rho_n(i,j,k)*vv(i,j,k)
          rhow_n(i,j,k) = rho_n(i,j,k)*ww(i,j,k)
          rhoe_n(i,j,k) = rho_n(i,j,k)*ecalc_tro(Tmp(i,j,k),rho_n(i,j,k))

          !if (is_bc_wall(2,1)) then
          !   uu(i,1,k)=0.0_wp
          !   vv(i,1,k)=0.0_wp
          !   rhou_n(i,1,k)=0.0_wp
          !   rhov_n(i,1,k)=0.0_wp
          !endif
          !if (is_bc_wall(2,2)) then
          !   uu(i,nx,k)=0.0_wp
          !   vv(i,nx,k)=0.0_wp
          !   rhou_n(i,nx,k)=0.0_wp
          !   rhov_n(i,nx,k)=0.0_wp
          !endif


          !rho_n(i+1,j,k) = rho_n(i,j,k)
          !rhou_n(i+1,j,k) = rho_n(i,j,k)*uu(i,j,k)
          !rhov_n(i+1,j,k) = rho_n(i,j,k)*vv(i,j,k)
          !rhow_n(i+1,j,k) = rho_n(i,j,k)*ww(i,j,k)
          !rhoe_n(i+1,j,k) = rho_n(i,j,k)*ecalc_tro(Tmp(i,j,k),rho_n(i,j,k))



          !rho_n(i,j,k) = rho_ref!rho_n(i,j-1,k)
          !prs(i,j,k) = p_ref!prs(i,j-1,k)
          !Tmp(i,j,k) = T_ref!Tmp(i,j-1,k)
          !rhou_n(i,j,k) = rho_ref*u_ref!rhou_n(i,j-1,k)
          !rhov_n(i,j,k) = 0!rhov_n(i,j-1,k)
          !rhoe_n(i,j,k) = rho_ref*cvfg*T_ref!rhoe_n(i,j-1,k)
          !uu(i,j,k) = u_ref 
          !vv(i,j,k) = 0     
          !ww(i,j,k) = 0  

          ! Neumann
          !rho_n(i,j,k) = rho_n(i+1,j,k)
          !prs(i,j,k) = prs(i+1,j,k)
          !Tmp(i,j,k) = Tmp(i+1,j,k)
          !rhou_n(i,j,k) = rhou_n(i+1,j,k)
          !rhov_n(i,j,k) = rhov_n(i+1,j,k)
          !rhoe_n(i,j,k) = rhoe_n(i+1,j,k)
          !uu(i,j,k) = uu(i+1,j,k) 
          !vv(i,j,k) = vv(i+1,j,k) 
          !ww(i,j,k) = ww(i+1,j,k) 



          ! Extrapolate to ghost points
          !do a=-5,-1
          !   rho_n(i+a,j,k) = rho_n(i,j,k)
          !   rhou_n(i+a,j,k) = rhou_n(i,j,k)
          !   rhov_n(i+a,j,k) = rhov_n(i,j,k)
          !   rhow_n(i+a,j,k) = rhow_n(i,j,k)
          !   rhoe_n(i+a,j,k) = rhoe_n(i,j,k)
          !enddo

          ! compute the rest
          visc(i,j,k) = viscosity_law(Tmp(i,j,k),rho_n(i,j,k))*diffscale
          cok(i,j,k)  = thconductivity(visc(i,j,k),Tmp(i,j,k),rho_n(i,j,k))
          
          Krho(i,j,k) = 0.0_wp
          Krhou(i,j,k) = 0.0_wp
          Krhov(i,j,k) = 0.0_wp
          Krhow(i,j,k) = 0.0_wp
          Krhoe(i,j,k) = 0.0_wp

          ! Krho(i+1,j,k) = 0.0_wp
          !Krhou(i+1,j,k) = 0.0_wp
          !Krhov(i+1,j,k) = 0.0_wp
          !Krhow(i+1,j,k) = 0.0_wp
          !Krhoe(i+1,j,k) = 0.0_wp
       enddo
    enddo
  
  end subroutine bc_inlet_pfg_imin

  subroutine bc_inlet_pfg_imax
  !==============================================================================
    !> Inlet boundary condition at imin for PFG eos
    !> Set free stream static quantities
  !==============================================================================

    ! Subscripts:
    ! ref: free stream values, prescribed in param.ini
    ! i  : interior of domain, so neighbour node
    ! bc : boundary value

    ! pressure from free stream
    p_ref = pcalc_tro(T_ref,rho_ref)
  
    ! direction of the velocity at the boundary (WARNING: THIS SHOULD COME FROM THE INPUT FILE)
    nb = (/1.0_wp,0.0_wp,0.0_wp/)

    ! velocity vector at free stream
    vel_ref = u_ref*nb
  
    i=nx
    ! Update inlet points
    ! ==================  
    do k=1,nz
       !do j=ndy_e,nfy_e
       do j=1-ngh,ny+ngh

          ! boundary normals
          nbc = (/nxn_imax(j,k),nyn_imax(j,k),nzn_imax(j,k)/)

          ! velocity vector at the interior nodes
          veli = (/uu(i-1,j,k),vv(i-1,j,k),ww(i-1,j,k)/)

          ! exrapolate R+ from the interior nodes
          Rp = dot_product(veli,nb)+2.0_wp*c_(i-1,j,k)/gam1

          ! exrapolate R- from the prescribed free stream
          c_ref = sqrt(c2calc_tro(T_ref,rho_ref))
          Rm = dot_product(vel_ref,nb)-2.0_wp*c_ref/gam1

          ! Compute velocity and speed of sound at boundary
          velbc = (Rp+Rm)/2.0_wp
          c_(nx,j,k) = gam1*(Rp-Rm)/4.0_wp

          ! Extrapolate entropy from free stream: sbc = s_ref
          s = scalc_tro(T_ref,rho_ref)

          ! Compute density, pressure and temperature based on s and c
          rho_n(i,j,k) = (c_(i,j,k)**2/(gam*s))**(1.0_wp/gam1)
          prs(i,j,k) = rho_n(i,j,k)*c_(i,j,k)**2/gam
          Tmp(i,j,k) = tcalc_pro(prs(i,j,k),rho_n(i,j,k),Tmp(i-1,j,k))

          ! Compute velocity and conservative variables
          !uu(i,j,k) = velbc*nb(1)
          !vv(i,j,k) = velbc*nb(2)
          !ww(i,j,k) = velbc*nb(3)
          uu(i,j,k) = vel_ref(1)+velbc*nbc(1)-dot_product(vel_ref,nbc)*nbc(1)
          vv(i,j,k) = vel_ref(2)+velbc*nbc(2)-dot_product(vel_ref,nbc)*nbc(2)
          ww(i,j,k) = vel_ref(3)+velbc*nbc(3)-dot_product(vel_ref,nbc)*nbc(3)
          rhou_n(i,j,k) = rho_n(i,j,k)*uu(i,j,k)
          rhov_n(i,j,k) = rho_n(i,j,k)*vv(i,j,k)
          rhow_n(i,j,k) = rho_n(i,j,k)*ww(i,j,k)
          rhoe_n(i,j,k) = rho_n(i,j,k)*ecalc_tro(Tmp(i,j,k),rho_n(i,j,k))



          rho_n(i,j,k) = rho_ref!rho_n(i,j-1,k)
          prs(i,j,k) = p_ref!prs(i,j-1,k)
          Tmp(i,j,k) = T_ref!Tmp(i,j-1,k)
          rhou_n(i,j,k) = rho_ref*u_ref!rhou_n(i,j-1,k)
          rhov_n(i,j,k) = 0!rhov_n(i,j-1,k)
          rhoe_n(i,j,k) = rho_ref*cvfg*T_ref!rhoe_n(i,j-1,k)
          uu(i,j,k) = u_ref 
          vv(i,j,k) = 0     
          ww(i,j,k) = 0    



          ! compute the rest
          visc(i,j,k) = viscosity_law(Tmp(i,j,k),rho_n(i,j,k))*diffscale
          cok(i,j,k)  = thconductivity(visc(i,j,k),Tmp(i,j,k),rho_n(i,j,k))
          
          Krho(i,j,k) = 0.0_wp
          Krhou(i,j,k) = 0.0_wp
          Krhov(i,j,k) = 0.0_wp
          Krhow(i,j,k) = 0.0_wp
          Krhoe(i,j,k) = 0.0_wp

       enddo
    enddo
  
  end subroutine bc_inlet_pfg_imax

  subroutine bc_inlet_pfg_jmin
  !==============================================================================
    !> Inlet boundary condition at imin for PFG eos
    !> Set free stream static quantities
  !==============================================================================

    ! Subscripts:
    ! ref: free stream values, prescribed in param.ini
    ! i  : interior of domain, so neighbour node
    ! bc : boundary value

    ! pressure from free stream
    p_ref = pcalc_tro(T_ref,rho_ref)
  
    ! direction of the velocity at the boundary (WARNING: THIS SHOULD COME FROM THE INPUT FILE)
    nb = (/1.0_wp,0.0_wp,0.0_wp/)

    ! velocity vector at free stream
    vel_ref = u_ref*nb

    j=1 
    ! Update inlet points
    ! ==================  
    do k=1,nz
       !do i=ndx_e,nfx_e
       do i=1-ngh,nx+ngh

          ! boundary normals
          nbc = (/nxn_jmin(i,k),nyn_jmin(i,k),nzn_jmin(i,k)/)
       
          ! velocity vector at the interior nodes
          veli = (/uu(i,j+1,k),vv(i,j+1,k),ww(i,j+1,k)/)

          ! exrapolate R+ from the interior nodes
          Rp = dot_product(veli,nb)+2.0_wp*c_(i,j+1,k)/gam1

          ! exrapolate R- from the prescribed free stream
          c_ref = sqrt(c2calc_tro(T_ref,rho_ref))
          Rm = dot_product(vel_ref,nb)-2.0_wp*c_ref/gam1

          ! Compute velocity and speed of sound at boundary
          velbc = (Rp+Rm)/2.0_wp
          c_(i,j,k) = gam1*(Rp-Rm)/4.0_wp

          ! Extrapolate entropy from free stream: sbc = s_ref
          s = scalc_tro(T_ref,rho_ref)

          ! Compute density, pressure and temperature based on s and c
          rho_n(i,j,k) = (c_(i,j,k)**2/(gam*s))**(1.0_wp/gam1)
          prs(i,j,k) = rho_n(i,j,k)*c_(i,j,k)**2/gam
          Tmp(i,j,k) = tcalc_pro(prs(i,j,k),rho_n(i,j,k),Tmp(i,j+1,k))

          ! Compute velocity and conservative variables
          !uu(i,j,k) = velbc*nb(1)
          !vv(i,j,k) = velbc*nb(2)
          !ww(i,j,k) = velbc*nb(3)
          uu(i,j,k) = vel_ref(1)+velbc*nbc(1)-dot_product(vel_ref,nbc)*nbc(1)
          vv(i,j,k) = vel_ref(2)+velbc*nbc(2)-dot_product(vel_ref,nbc)*nbc(2)
          ww(i,j,k) = vel_ref(3)+velbc*nbc(3)-dot_product(vel_ref,nbc)*nbc(3)
          rhou_n(i,j,k) = rho_n(i,j,k)*uu(i,j,k)
          rhov_n(i,j,k) = rho_n(i,j,k)*vv(i,j,k)
          rhow_n(i,j,k) = rho_n(i,j,k)*ww(i,j,k)
          rhoe_n(i,j,k) = rho_n(i,j,k)*ecalc_tro(Tmp(i,j,k),rho_n(i,j,k))



          rho_n(i,j,k) = rho_ref!rho_n(i,j-1,k)
          prs(i,j,k) = p_ref!prs(i,j-1,k)
          Tmp(i,j,k) = T_ref!Tmp(i,j-1,k)
          rhou_n(i,j,k) = rho_ref*u_ref!rhou_n(i,j-1,k)
          rhov_n(i,j,k) = 0!rhov_n(i,j-1,k)
          rhoe_n(i,j,k) = rho_ref*cvfg*T_ref!rhoe_n(i,j-1,k)
          uu(i,j,k) = u_ref 
          vv(i,j,k) = 0     
          ww(i,j,k) = 0    



          ! compute the rest
          visc(i,j,k) = viscosity_law(Tmp(i,j,k),rho_n(i,j,k))*diffscale
          cok(i,j,k)  = thconductivity(visc(i,j,k),Tmp(i,j,k),rho_n(i,j,k))
          
          Krho(i,j,k) = 0.0_wp
          Krhou(i,j,k) = 0.0_wp
          Krhov(i,j,k) = 0.0_wp
          Krhow(i,j,k) = 0.0_wp
          Krhoe(i,j,k) = 0.0_wp
  
       enddo
    enddo
  
  end subroutine bc_inlet_pfg_jmin

  subroutine bc_inlet_pfg_jmax
  !==============================================================================
    !> Inlet boundary condition at imin for PFG eos
    !> Set free stream static quantities
  !==============================================================================

    ! Subscripts:
    ! ref: free stream values, prescribed in param.ini
    ! i  : interior node
    ! bc : boundary value

    ! set total pressure and temperature
    ! pressure from free stream
    p_ref = pcalc_tro(T_ref,rho_ref)
  
    ! direction of the velocity at the boundary (WARNING: THIS SHOULD COME FROM THE INPUT FILE)
    nb = (/1.0_wp,0.0_wp,0.0_wp/)

    ! velocity vector at free stream
    vel_ref = u_ref*nb
  
    j=ny
    ! Update inlet points
    ! ==================  
    do k=1,nz
       !do i=ndx_e,nfx_e
       do i=1-ngh,nx+ngh

          ! boundary normals
          nbc = (/-nxn_jmax(i,k),nyn_jmax(i,k),nzn_jmax(i,k)/)

          ! velocity vector at the interior nodes
          veli = (/uu(i,j-1,k),vv(i,j-1,k),ww(i,j-1,k)/)

          ! exrapolate R+ from the interior nodes
          Rp = dot_product(veli,nbc)+2.0_wp*c_(i,j-1,k)/gam1
          !print*,dot_product(veli,nbc),i,iproc

          ! exrapolate R- from the prescribed free stream
          c_ref = sqrt(c2calc_tro(T_ref,rho_ref))
          Rm = dot_product(vel_ref,nbc)-2.0_wp*c_ref/gam1

          ! Compute velocity and speed of sound at boundary
          velbc = (Rp+Rm)/2.0_wp
          c_(i,j,k) = gam1*(Rp-Rm)/4.0_wp
          !print*,c_(i,j,k)

          ! Extrapolate entropy from free stream: sbc = s_ref
          s = scalc_tro(T_ref,rho_ref)

          ! Compute density, pressure and temperature based on s and c
          rho_n(i,j,k) = (c_(i,j,k)**2/(gam*s))**(1.0_wp/gam1)
          prs(i,j,k) = rho_n(i,j,k)*c_(i,j,k)**2/gam
          Tmp(i,j,k) = tcalc_pro(prs(i,j,k),rho_n(i,j,k),Tmp(i,j-1,k))

          ! Compute velocity and conservative variables
          !uu(i,j,k) = velbc*nb(1)
          !vv(i,j,k) = velbc*nb(2)
          !ww(i,j,k) = velbc*nb(3)
          uu(i,j,k) = vel_ref(1)+velbc*nbc(1)-dot_product(vel_ref,nbc)*nbc(1)
          vv(i,j,k) = vel_ref(2)+velbc*nbc(2)-dot_product(vel_ref,nbc)*nbc(2)
          ww(i,j,k) = vel_ref(3)+velbc*nbc(3)-dot_product(vel_ref,nbc)*nbc(3)
          !if (uu(i,j,k).ne.u_ref) then
          !   print*,i
          !endif
          print*,nbc
          rhou_n(i,j,k) = rho_n(i,j,k)*uu(i,j,k)
          rhov_n(i,j,k) = rho_n(i,j,k)*vv(i,j,k)
          rhow_n(i,j,k) = rho_n(i,j,k)*ww(i,j,k)
          rhoe_n(i,j,k) = rho_n(i,j,k)*ecalc_tro(Tmp(i,j,k),rho_n(i,j,k))

          !print*,veli(2),c_(i,j,k),rho_n(i,j,k)


          rho_n(i,j,k) = rho_ref!rho_n(i,j-1,k)
          prs(i,j,k) = p_ref!prs(i,j-1,k)
          Tmp(i,j,k) = T_ref!Tmp(i,j-1,k)
          rhou_n(i,j,k) = rho_ref*u_ref!rhou_n(i,j-1,k)
          rhov_n(i,j,k) = 0!rhov_n(i,j-1,k)
          rhoe_n(i,j,k) = rho_ref*cvfg*T_ref!rhoe_n(i,j-1,k)
          uu(i,j,k) = u_ref 
          vv(i,j,k) = 0     
          ww(i,j,k) = 0    



          ! compute the rest
          visc(i,ny,k) = viscosity_law(Tmp(i,ny,k),rho_n(i,ny,k))*diffscale
          cok(i,ny,k)  = thconductivity(visc(i,ny,k),Tmp(i,ny,k),rho_n(i,ny,k))
          
          Krho(i,ny,k) = 0.0_wp
          Krhou(i,ny,k) = 0.0_wp
          Krhov(i,ny,k) = 0.0_wp
          Krhow(i,ny,k) = 0.0_wp
          Krhoe(i,ny,k) = 0.0_wp

       enddo
    enddo
  
  end subroutine bc_inlet_pfg_jmax
 
 
!  ! PRS
!  ! ===
!
!  subroutine bc_inlet_prs_imin
!  !==============================================================================
!    !> Inlet boundary condition at imin for PRS eos
!    !> Set total pressure and temperature
!  !==============================================================================
!  
!    ! set total pressure and temperature
!    p_tot = 250000_wp!*p_ref
!    T_tot = T_ref/0.833_wp
!  
!    ! direction of the velocity at the boundary 
!    nb = (/1.0_wp,0.0_wp,0.0_wp/)
!  
!    ! Initialize
!    err = 1.0_wp
!    tol = 1e-6_wp
!    do k=1,nz
!       do j=1,ny
!          ! Use neighbour point for test values
!          rotent = rho_n(2,j,k)
!          ptent  = prs(2,j,k)
!          Ttent  = Tmp(2,j,k)
!          ! initialize with total values
!          prs(1,j,k) = p_tot
!          Tmp(1,j,k) = T_tot
!          rho_tot = rocalc_pt(p_tot,T_tot,rotent)
!          rho_n(1,j,k) = rho_tot
!       enddo
!    enddo
!    
!    do k=1,nz
!       do j=1,ny
!          rhobc = rho_n(1,j,k)
!          pbc  = prs(1,j,k)
!          Tbc  = Tmp(1,j,k)
!  
!          ! velocity vector at the interior nodes
!          veli = (/uu(2,j,k),vv(2,j,k),ww(2,j,k)/)
!  
!          ! compute total enthalpy without adiabatic condition, as H0_i = E0_i + P/rho
!          ! H_tot = ecalc_pro(p_tot,rho_tot,T_tot) + prs(2,j,k)/rho_n(2,j,k)
!  
!          do while (err.gt.tol)
!             
!             ptmp = pbc
!  
!             ! exrapolate R+ from the interior nodes, so dR=0
!             ! dR = u_i-u_bc + c_i/rho_i*(rho_i-rho_bc)
!             ci = sqrt(c2calc_tro(Tmp(2,j,k),rho_n(2,j,k)))
!             velbc = sqrt(dot_product(veli,veli)) + ci/rho_n(2,j,k)*(rho_n(2,j,k)-rhobc)
!     
!             ! compute total enthalpy, apply adiabaticity so H0_bc = H0_i
!             H_tot = ecalc_pro(prs(2,j,k),rho_n(2,j,k),Ttent) + prs(2,j,k)/rho_n(2,j,k) &
!                     + 0.5_wp*dot_product(veli,veli)
!             h = H_tot - 0.5_wp*velbc**2
!  
!             ! update thermo variables
!             Tbc  = tcalc_ph(pbc,h,rotent,Ttent)
!             rhobc = rocalc_pt(pbc,Tbc,rotent)
!             pbc  = pcalc_tro(Tbc,rhobc)
!          
!             err = (pbc-ptmp)/ptmp
!          enddo
!  
!       ! update thermo variables at boundary 
!       rho_n(1,j,k) = rhobc
!       prs(1,j,k)   = pbc 
!       Tmp(1,j,k)   = Tbc 
!  
!       ! compute velocity and momentum variables
!       uu(1,j,k) = velbc*nb(1)
!       vv(1,j,k) = velbc*nb(2)
!       ww(1,j,k) = velbc*nb(3)
!       rhou_n(1,j,k) = rho_n(1,j,k)*uu(1,j,k)
!       rhov_n(1,j,k) = rho_n(1,j,k)*vv(1,j,k)
!       rhow_n(1,j,k) = rho_n(1,j,k)*ww(1,j,k)
!       
!       ! compute the rest
!       visc(1,j,k)= viscosity_law(Tmp(1,j,k),rho_n(1,j,k))
!       cok(1,j,k) = thconductivity(visc(1,j,k),Tmp(1,j,k),rho_n(1,j,k))
!       c_(1,j,k)  = sqrt(c2calc_tro(Tmp(1,j,k),rho_n(1,j,k)))
!       rhoe_n(1,j,k) = rho_n(1,j,k)*ecalc_tro(Tmp(1,j,k),rho_n(1,j,k))
!       
!       Krho(1,j,k)  = 0.0_wp
!       Krhou(1,j,k) = 0.0_wp
!       Krhov(1,j,k) = 0.0_wp
!       Krhow(1,j,k) = 0.0_wp
!       Krhoe(1,j,k) = 0.0_wp
!  
!       enddo
!    enddo
!  
!  end subroutine bc_inlet_prs_imin
!  
!  subroutine bc_inlet_prs_imax
!  !==============================================================================
!    !> Inlet boundary condition at imax for PRS eos
!    !> Set total pressure and temperature
!  !==============================================================================
!
!    ! set total pressure and temperature
!    p_tot = 250000_wp!*p_ref
!    T_tot = T_ref/0.833_wp
!  
!    ! direction of the velocity at the boundary 
!    nb = (/1.0_wp,0.0_wp,0.0_wp/)
!  
!    ! Initialize
!    err = 1.0_wp
!    tol = 1e-6_wp
!    do k=1,nz
!       do j=1,ny
!          ! Use neighbour point for test values
!          rotent = rho_n(nx-1,j,k)
!          ptent  = prs(nx-1,j,k)
!          Ttent  = Tmp(nx-1,j,k)
!          ! initialize with total values
!          prs(nx,j,k) = p_tot
!          Tmp(nx,j,k) = T_tot
!          rho_tot = rocalc_pt(p_tot,T_tot,rotent)
!          rho_n(nx,j,k) = rho_tot
!       enddo
!    enddo
!    
!    do k=1,nz
!       do j=1,ny
!          rhobc = rho_n(nx,j,k)
!          pbc  = prs(nx,j,k)
!          Tbc  = Tmp(nx,j,k)
!  
!          ! velocity vector at the interior nodes
!          veli = (/uu(nx-1,j,k),vv(nx-1,j,k),ww(nx-1,j,k)/)
!  
!          ! compute total enthalpy without adiabatic condition, as H0_i = E0_i + P/rho
!          ! H_tot = ecalc_pro(p_tot,rho_tot,T_tot) + prs(nx-1,j,k)/rho_n(nx-1,j,k)
!  
!          do while (err.gt.tol)
!             
!             ptmp = pbc
!  
!             ! exrapolate R+ from the interior nodes, so dR=0
!             ! dR = u_i-u_bc + c_i/rho_i*(rho_i-rho_bc)
!             ci = sqrt(c2calc_tro(Tmp(nx-1,j,k),rho_n(nx-1,j,k)))
!             velbc = sqrt(dot_product(veli,veli)) + ci/rho_n(nx-1,j,k)*(rho_n(nx-1,j,k)-rhobc)
!     
!             ! compute total enthalpy, apply adiabaticity so H0_bc = H0_i
!             H_tot = ecalc_pro(prs(nx-1,j,k),rho_n(nx-1,j,k),Ttent) + prs(nx-1,j,k)/rho_n(nx-1,j,k) &
!                     + 0.5_wp*dot_product(veli,veli)
!             h = H_tot - 0.5_wp*velbc**2
!  
!             ! update thermo variables
!             Tbc  = tcalc_ph(pbc,h,rotent,Ttent)
!             rhobc = rocalc_pt(pbc,Tbc,rotent)
!             pbc  = pcalc_tro(Tbc,rhobc)
!          
!             err = (pbc-ptmp)/ptmp
!          enddo
!  
!       ! update thermo variables at boundary 
!       rho_n(nx,j,k) = rhobc
!       prs(nx,j,k)   = pbc 
!       Tmp(nx,j,k)   = Tbc 
!  
!       ! compute velocity and momentum variables
!       uu(nx,j,k) = velbc*nb(1)
!       vv(nx,j,k) = velbc*nb(2)
!       ww(nx,j,k) = velbc*nb(3)
!       rhou_n(nx,j,k) = rho_n(nx,j,k)*uu(nx,j,k)
!       rhov_n(nx,j,k) = rho_n(nx,j,k)*vv(nx,j,k)
!       rhow_n(nx,j,k) = rho_n(nx,j,k)*ww(nx,j,k)
!       
!       ! compute the rest
!       visc(nx,j,k)= viscosity_law(Tmp(nx,j,k),rho_n(nx,j,k))
!       cok(nx,j,k) = thconductivity(visc(nx,j,k),Tmp(nx,j,k),rho_n(nx,j,k))
!       c_(nx,j,k)  = sqrt(c2calc_tro(Tmp(nx,j,k),rho_n(nx,j,k)))
!       rhoe_n(nx,j,k) = rho_n(nx,j,k)*ecalc_tro(Tmp(nx,j,k),rho_n(nx,j,k))
!       
!       Krho(nx,j,k)  = 0.0_wp
!       Krhou(nx,j,k) = 0.0_wp
!       Krhov(nx,j,k) = 0.0_wp
!       Krhow(nx,j,k) = 0.0_wp
!       Krhoe(nx,j,k) = 0.0_wp
!  
!       enddo
!    enddo
!  
!  end subroutine bc_inlet_prs_imax
!  
!  
!  subroutine bc_inlet_prs_jmin
!  !==============================================================================
!    !> Inlet boundary condition at jmin for PRS eos
!    !> Set total pressure and temperature
!  !==============================================================================
!
!    ! set total pressure and temperature
!    p_tot = 250000_wp!*p_ref
!    T_tot = T_ref/0.833_wp
!  
!    ! direction of the velocity at the boundary 
!    nb = (/1.0_wp,0.0_wp,0.0_wp/)
!  
!    ! Initialize
!    err = 1.0_wp
!    tol = 1e-6_wp
!    do k=1,nz
!       do i=1,nx
!          ! Use neighbour point for test values
!          rotent = rho_n(i,2,k)
!          ptent  = prs(i,2,k)
!          Ttent  = Tmp(i,2,k)
!          ! initialize with total values
!          prs(i,1,k) = p_tot
!          Tmp(i,1,k) = T_tot
!          rho_tot = rocalc_pt(p_tot,T_tot,rotent)
!          rho_n(i,1,k) = rho_tot
!       enddo
!    enddo
!    
!    do k=1,nz
!       do i=1,nx
!          rhobc = rho_n(i,1,k)
!          pbc  = prs(i,1,k)
!          Tbc  = Tmp(i,1,k)
!  
!          ! velocity vector at the interior nodes
!          veli = (/uu(i,2,k),vv(i,2,k),ww(i,2,k)/)
!  
!          ! compute total enthalpy without adiabatic condition, as H0_i = E0_i + P/rho
!          ! H_tot = ecalc_pro(p_tot,rho_tot,T_tot) + prs(i,2,k)/rho_n(i,2,k)
!  
!          do while (err.gt.tol)
!             
!             ptmp = pbc
!  
!             ! exrapolate R+ from the interior nodes, so dR=0
!             ! dR = u_i-u_bc + c_i/rho_i*(rho_i-rho_bc)
!             ci = sqrt(c2calc_tro(Tmp(i,2,k),rho_n(i,2,k)))
!             velbc = sqrt(dot_product(veli,veli)) + ci/rho_n(i,2,k)*(rho_n(i,2,k)-rhobc)
!     
!             ! compute total enthalpy, apply adiabaticity so H0_bc = H0_i
!             H_tot = ecalc_pro(prs(i,2,k),rho_n(i,2,k),Ttent) + prs(i,2,k)/rho_n(i,2,k) &
!                     + 0.5_wp*dot_product(veli,veli)
!             h = H_tot - 0.5_wp*velbc**2
!  
!             ! update thermo variables
!             Tbc  = tcalc_ph(pbc,h,rotent,Ttent)
!             rhobc = rocalc_pt(pbc,Tbc,rotent)
!             pbc  = pcalc_tro(Tbc,rhobc)
!          
!             err = (pbc-ptmp)/ptmp
!          enddo
!  
!       ! update thermo variables at boundary 
!       rho_n(i,1,k) = rhobc
!       prs(i,1,k)   = pbc 
!       Tmp(i,1,k)   = Tbc 
!  
!       ! compute velocity and momentum variables
!       uu(i,1,k) = velbc*nb(1)
!       vv(i,1,k) = velbc*nb(2)
!       ww(i,1,k) = velbc*nb(3)
!       rhou_n(i,1,k) = rho_n(i,1,k)*uu(i,1,k)
!       rhov_n(i,1,k) = rho_n(i,1,k)*vv(i,1,k)
!       rhow_n(i,1,k) = rho_n(i,1,k)*ww(i,1,k)
!       
!       ! compute the rest
!       visc(i,1,k)= viscosity_law(Tmp(i,1,k),rho_n(i,1,k))
!       cok(i,1,k) = thconductivity(visc(i,1,k),Tmp(i,1,k),rho_n(i,1,k))
!       c_(i,1,k)  = sqrt(c2calc_tro(Tmp(i,1,k),rho_n(i,1,k)))
!       rhoe_n(i,1,k) = rho_n(i,1,k)*ecalc_tro(Tmp(i,1,k),rho_n(i,1,k))
!       
!       Krho(i,1,k)  = 0.0_wp
!       Krhou(i,1,k) = 0.0_wp
!       Krhov(i,1,k) = 0.0_wp
!       Krhow(i,1,k) = 0.0_wp
!       Krhoe(i,1,k) = 0.0_wp
!  
!       enddo
!    enddo
!  
!  end subroutine bc_inlet_prs_jmin
!  
!  
!  subroutine bc_inlet_prs_jmax
!  !==============================================================================
!    !> Inlet boundary condition at jmax for PRS eos
!    !> Set total pressure and temperature
!  !==============================================================================
!
!    ! set total pressure and temperature
!    p_tot = 250000_wp!*p_ref
!    T_tot = T_ref/0.833_wp
!  
!    ! direction of the velocity at the boundary 
!    nb = (/1.0_wp,0.0_wp,0.0_wp/)
!  
!    ! Initialize
!    err = 1.0_wp
!    tol = 1e-6_wp
!    do k=1,nz
!       do i=1,nx
!          ! Use neighbour point for test values
!          rotent = rho_n(i,ny-1,k)
!          ptent  = prs(i,ny-1,k)
!          Ttent  = Tmp(i,ny-1,k)
!          ! initialize with total values
!          prs(i,ny,k) = p_tot
!          Tmp(i,ny,k) = T_tot
!          rho_tot = rocalc_pt(p_tot,T_tot,rotent)
!          rho_n(i,ny,k) = rho_tot
!       enddo
!    enddo
!    
!    do k=1,nz
!       do i=1,nx
!          rhobc = rho_n(i,ny,k)
!          pbc  = prs(i,ny,k)
!          Tbc  = Tmp(i,ny,k)
!  
!          ! velocity vector at the interior nodes
!          veli = (/uu(i,ny-1,k),vv(i,ny-1,k),ww(i,ny-1,k)/)
!  
!          ! compute total enthalpy without adiabatic condition, as H0_i = E0_i + P/rho
!          ! H_tot = ecalc_pro(p_tot,rho_tot,T_tot) + prs(i,ny-1,k)/rho_n(i,ny-1,k)
!  
!          do while (err.gt.tol)
!             
!             ptmp = pbc
!  
!             ! exrapolate R+ from the interior nodes, so dR=0
!             ! dR = u_i-u_bc + c_i/rho_i*(rho_i-rho_bc)
!             ci = sqrt(c2calc_tro(Tmp(i,ny-1,k),rho_n(i,ny-1,k)))
!             velbc = sqrt(dot_product(veli,veli)) + ci/rho_n(i,ny-1,k)*(rho_n(i,ny-1,k)-rhobc)
!     
!             ! compute total enthalpy, apply adiabaticity so H0_bc = H0_i
!             H_tot = ecalc_pro(prs(i,ny-1,k),rho_n(i,ny-1,k),Ttent) + prs(i,ny-1,k)/rho_n(i,ny-1,k) &
!                     + 0.5_wp*dot_product(veli,veli)
!             h = H_tot - 0.5_wp*velbc**2
!  
!             ! update thermo variables
!             Tbc  = tcalc_ph(pbc,h,rotent,Ttent)
!             rhobc = rocalc_pt(pbc,Tbc,rotent)
!             pbc  = pcalc_tro(Tbc,rhobc)
!          
!             err = (pbc-ptmp)/ptmp
!          enddo
!  
!       ! update thermo variables at boundary 
!       rho_n(i,ny,k) = rhobc
!       prs(i,ny,k)   = pbc 
!       Tmp(i,ny,k)   = Tbc 
!  
!       ! compute velocity and momentum variables
!       uu(i,ny,k) = velbc*nb(1)
!       vv(i,ny,k) = velbc*nb(2)
!       ww(i,ny,k) = velbc*nb(3)
!       rhou_n(i,ny,k) = rho_n(i,ny,k)*uu(i,ny,k)
!       rhov_n(i,ny,k) = rho_n(i,ny,k)*vv(i,ny,k)
!       rhow_n(i,ny,k) = rho_n(i,ny,k)*ww(i,ny,k)
!       
!       ! compute the rest
!       visc(i,ny,k)= viscosity_law(Tmp(i,ny,k),rho_n(i,ny,k))
!       cok(i,ny,k) = thconductivity(visc(i,ny,k),Tmp(i,ny,k),rho_n(i,ny,k))
!       c_(i,ny,k)  = sqrt(c2calc_tro(Tmp(i,ny,k),rho_n(i,ny,k)))
!       rhoe_n(i,ny,k) = rho_n(i,ny,k)*ecalc_tro(Tmp(i,ny,k),rho_n(i,ny,k))
!       
!       Krho(i,ny,k)  = 0.0_wp
!       Krhou(i,ny,k) = 0.0_wp
!       Krhov(i,ny,k) = 0.0_wp
!       Krhow(i,ny,k) = 0.0_wp
!       Krhoe(i,ny,k) = 0.0_wp
!  
!       enddo
!    enddo
!  
!  end subroutine bc_inlet_prs_jmax
  
  ! ======= !
  ! Outlets !
  ! ======= !

  ! PFG
  ! === 
  
  subroutine bc_outlet_pfg_imin
  !==============================================================================
    !> Outlet boundary condition at imin for PFG eos
    !> Use interior node information
  !==============================================================================

    ! Subscripts:
    ! ref: free stream values, prescribed in param.ini
    ! i  : interior of domain, so neighbour node
    ! bc : boundary value

    ! pressure from free stream
    p_ref = pcalc_tro(T_ref,rho_ref)
  
    ! direction of the velocity at the boundary (WARNING: THIS SHOULD COME FROM THE INPUT FILE)
    nb = (/1.0_wp,0.0_wp,0.0_wp/)

    ! velocity vector at free stream
    vel_ref = u_ref*nb

    i=1
    ! Update outlet points
    ! ==================  
    do k=1,nz
       do j=1,ny

          ! boundary normals
          nbc = (/nxn_imin(j,k),nyn_imin(j,k),nzn_imin(j,k)/)
       
          ! velocity vector at the interior nodes
          veli = (/uu(i+i,j,k),vv(i+i,j,k),ww(i+i,j,k)/)

          ! exrapolate R- from the interior nodes
          Rm = dot_product(veli,nbc)-2.0_wp*c_(i+i,j,k)/gam1

          ! exrapolate R+ from the prescribed free stream
          c_ref = sqrt(c2calc_tro(T_ref,rho_ref))
          Rp = dot_product(vel_ref,nbc)+2.0_wp*c_ref/gam1

          ! Compute velocity and speed of sound at boundary
          velbc = (Rp+Rm)/2.0_wp
          c_(i,j,k) = gam1*(Rp-Rm)/4.0_wp

          ! Extrapolate entropy from interior: sbc = si
          s = scalc_tro(Tmp(2,j,k),rho_n(2,j,k))

          ! Compute density, pressure and temperature based on s and c
          rho_n(i,j,k) = (c_(i,j,k)**2/(gam*s))**(1.0_wp/gam1)
          prs(i,j,k) = rho_n(i,j,k)*c_(i,j,k)**2/gam
          Tmp(i,j,k) = tcalc_pro(prs(i,j,k),rho_n(i,j,k),Tmp(i+i,j,k))

          ! Compute velocity and conservative variables
          uu(i,j,k) = veli(1)+velbc*nbc(1)-dot_product(veli,nbc)*nbc(1)
          vv(i,j,k) = veli(2)+velbc*nbc(2)-dot_product(veli,nbc)*nbc(2)
          ww(i,j,k) = veli(3)+velbc*nbc(3)-dot_product(veli,nbc)*nbc(3)
          rhou_n(i,j,k) = rho_n(i,j,k)*uu(i,j,k)
          rhov_n(i,j,k) = rho_n(i,j,k)*vv(i,j,k)
          rhow_n(i,j,k) = rho_n(i,j,k)*ww(i,j,k)
          rhoe_n(i,j,k) = rho_n(i,j,k)*ecalc_tro(Tmp(i,j,k),rho_n(i,j,k))


          ! Dirichlet
          rho_n(i,j,k) = rho_ref!rho_n(i,j-1,k)
          prs(i,j,k) = p_ref!prs(i,j-1,k)
          Tmp(i,j,k) = T_ref!Tmp(i,j-1,k)
          rhou_n(i,j,k) = rho_ref*u_ref!rhou_n(i,j-1,k)
          rhov_n(i,j,k) = 0!rhov_n(i,j-1,k)
          rhoe_n(i,j,k) = rho_ref*cvfg*T_ref!rhoe_n(i,j-1,k)
          uu(i,j,k) = u_ref 
          vv(i,j,k) = 0     
          ww(i,j,k) = 0    

          ! Neumann
          rho_n(i,j,k) = rho_n(i+1,j,k)
          prs(i,j,k) = prs(i+1,j,k)
          Tmp(i,j,k) = Tmp(i+1,j,k)
          rhou_n(i,j,k) = rhou_n(i+1,j,k)
          rhov_n(i,j,k) = rhov_n(i+1,j,k)
          rhoe_n(i,j,k) = rhoe_n(i+1,j,k)
          uu(i,j,k) = uu(i+1,j,k) 
          vv(i,j,k) = vv(i+1,j,k) 
          ww(i,j,k) = ww(i+1,j,k) 



          ! compute the rest
          visc(i,j,k) = viscosity_law(Tmp(i,j,k),rho_n(i,j,k))*diffscale
          cok(i,j,k)  = thconductivity(visc(i,j,k),Tmp(i,j,k),rho_n(i,j,k))
          
          Krho(i,j,k) = 0.0_wp
          Krhou(i,j,k) = 0.0_wp
          Krhov(i,j,k) = 0.0_wp
          Krhow(i,j,k) = 0.0_wp
          Krhoe(i,j,k) = 0.0_wp

       enddo
    enddo
  
  end subroutine bc_outlet_pfg_imin

  subroutine bc_outlet_pfg_imax
  !==============================================================================
    !> Outlet boundary condition at imax for PFG eos
    !> Use interior node information
  !==============================================================================

    ! Subscripts:
    ! ref: free stream values, prescribed in param.ini
    ! i  : interior of domain, so neighbour node
    ! bc : boundary value

    ! pressure from free stream
    p_ref = pcalc_tro(T_ref,rho_ref)

    ! direction of the velocity at the boundary (WARNING: THIS SHOULD COME FROM THE INPUT FILE)
    nb = (/1.0_wp,0.0_wp,0.0_wp/)

    ! velocity vector at free stream
    vel_ref = u_ref*nb

    i=nx 
    ! Update outlet points
    ! ==================  
    do k=1,nz
       do j=ndy_e,nfy_e

          ! boundary normals
          nbc = (/nxn_imax(j,k),nyn_imax(j,k),nzn_imax(j,k)/)
          nbc = (/1.0_wp,0.0_wp,0.0_wp/)
       
          ! velocity vector at the interior nodes
          veli = (/uu(i-1,j,k),vv(i-1,j,k),ww(i-1,j,k)/)

          ! exrapolate R- from the interior nodes
          Rm = dot_product(veli,nbc)-2.0_wp*c_(i-1,j,k)/gam1

          ! exrapolate R+ from the prescribed free stream
          c_ref = sqrt(c2calc_tro(T_ref,rho_ref))
          Rp = dot_product(vel_ref,nbc)+2.0_wp*c_ref/gam1

          ! Compute velocity and speed of sound at boundary
          velbc = (Rp+Rm)/2.0_wp
          c_(i,j,k) = gam1*(Rp-Rm)/4.0_wp
          !print*,c_(i,j,k),velbc,dot_product(veli,nb)

          ! Extrapolate entropy from interior: sbc = si
          s = scalc_tro(Tmp(i-1,j,k),rho_n(i-1,j,k))
          !s = scalc_tro(T_ref,rho_ref)

          ! Compute density, pressure and temperature based on s and c
          rho_n(i,j,k) = (c_(i,j,k)**2/(gam*s))**(1.0_wp/gam1)
          prs(i,j,k) = rho_n(i,j,k)*c_(i,j,k)**2/gam
          Tmp(i,j,k) = tcalc_pro(prs(i,j,k),rho_n(i,j,k),Tmp(i-1,j,k))
          !print*,rho_n(i,j,k),Tmp(i,j,k),prs(i,j,k)

          ! Compute velocity and conservative variables
          uu(i,j,k) = veli(1)+velbc*nbc(1)-dot_product(veli,nbc)*nbc(1)
          vv(i,j,k) = veli(2)+velbc*nbc(2)-dot_product(veli,nbc)*nbc(2)
          ww(i,j,k) = veli(3)+velbc*nbc(3)-dot_product(veli,nbc)*nbc(3)
          !uu(i,j,k) = velbc*nb(1)
          !vv(i,j,k) = velbc*nb(2)
          !ww(i,j,k) = velbc*nb(3)
          rhou_n(i,j,k) = rho_n(i,j,k)*uu(i,j,k)
          rhov_n(i,j,k) = rho_n(i,j,k)*vv(i,j,k)
          rhow_n(i,j,k) = rho_n(i,j,k)*ww(i,j,k)
          rhoe_n(i,j,k) = rho_n(i,j,k)*ecalc_tro(Tmp(i,j,k),rho_n(i,j,k))



          !rho_n(i,j,k) = rho_ref!rho_n(i,ny-1,k)
          !prs(i,j,k) = p_ref!prs(i,ny-1,k)
          !Tmp(i,j,k) = T_ref!Tmp(i,ny-1,k)
          !c_(i,j,k) = c_ref!Tmp(i,ny-1,k)
          !rhou_n(i,j,k) = rho_ref*u_ref!rhou_n(i,ny-1,k)
          !rhov_n(i,j,k) = 0!rhov_n(i,ny-1,k)
          !rhoe_n(i,j,k) = rho_ref*ecalc_tro(T_ref,rho_ref)!rhoe_n(i,ny-1,k)
          !uu(i,j,k) = u_ref
          !vv(i,j,k) = 0
          !ww(i,j,k) = 0
    
          ! Neumann
          rho_n(i,j,k) = rho_n(i-1,j,k)
          prs(i,j,k) = prs(i-1,j,k)
          Tmp(i,j,k) = Tmp(i-1,j,k)
          c_(i,j,k) = c_(i-1,j,k)
          rhou_n(i,j,k) = rhou_n(i-1,j,k)
          rhov_n(i,j,k) = rhov_n(i-1,j,k)
          rhoe_n(i,j,k) = rhoe_n(i-1,j,k)
          uu(i,j,k) = uu(i-1,j,k) 
          vv(i,j,k) = vv(i-1,j,k) 
          ww(i,j,k) = ww(i-1,j,k) 

          ! compute the rest
          visc(i,j,k) = viscosity_law(Tmp(i,j,k),rho_n(i,j,k))*diffscale
          cok(i,j,k)  = thconductivity(visc(i,j,k),Tmp(i,j,k),rho_n(i,j,k))
          
          Krho(i,j,k) = 0.0_wp
          Krhou(i,j,k) = 0.0_wp
          Krhov(i,j,k) = 0.0_wp
          Krhow(i,j,k) = 0.0_wp
          Krhoe(i,j,k) = 0.0_wp

       enddo
    enddo
  
  end subroutine bc_outlet_pfg_imax

  subroutine bc_outlet_pfg_jmin
  !==============================================================================
    !> Outlet boundary condition at jmin for PFG eos
    !> Use interior node information
  !==============================================================================

    ! Subscripts:
    ! ref: free stream values, prescribed in param.ini
    ! i  : interior of domain, so neighbour node
    ! bc : boundary value

    ! pressure from free stream
    p_ref = pcalc_tro(T_ref,rho_ref)

    ! direction of the velocity at the boundary (WARNING: THIS SHOULD COME FROM THE INPUT FILE)
    nb = (/1.0_wp,0.0_wp,0.0_wp/)

    ! velocity vector at free stream
    vel_ref = u_ref*nb

    j=1 
    ! Update outlet points
    ! ==================  
    do k=1,nz
       do i=ndx_e,nfx_e
       !do i=1-ngh,nx+ngh

          ! boundary normals
          nbc = (/nxn_jmin(i,k),-nyn_jmin(i,k),nzn_jmin(i,k)/)
       
          ! velocity vector at the interior nodes
          veli = (/uu(i,j+1,k),vv(i,j+1,k),ww(i,j+1,k)/)

          ! exrapolate R- from the interior nodes
          Rm = dot_product(veli,nbc)-2.0_wp*c_(i,j+1,k)/gam1

          ! exrapolate R+ from the prescribed free stream
          c_ref = sqrt(c2calc_tro(T_ref,rho_ref))
          Rp = dot_product(vel_ref,nbc)+2.0_wp*c_ref/gam1

          ! Compute velocity and speed of sound at boundary
          velbc = (Rp+Rm)/2.0_wp
          c_(i,j,k) = gam1*(Rp-Rm)/4.0_wp

          ! Extrapolate entropy from interior: sbc = si
          s = scalc_tro(Tmp(i,j+1,k),rho_n(i,j+1,k))

          ! Compute density, pressure and temperature based on s and c
          rho_n(nx,j,k) = (c_(i,j,k)**2/(gam*s))**(1.0_wp/gam1)
          prs(nx,j,k) = rho_n(i,j,k)*c_(i,j,k)**2/gam
          Tmp(i,j,k) = tcalc_pro(prs(i,j,k),rho_n(i,j,k),Tmp(i,j+1,k))

          ! Compute velocity and conservative variables
          uu(i,j,k) = veli(1)+velbc*nbc(1)-dot_product(veli,nbc)*nbc(1)
          vv(i,j,k) = veli(2)+velbc*nbc(2)-dot_product(veli,nbc)*nbc(2)
          ww(i,j,k) = veli(3)+velbc*nbc(3)-dot_product(veli,nbc)*nbc(3)
          rhou_n(i,j,k) = rho_n(i,j,k)*uu(i,j,k)
          rhov_n(i,j,k) = rho_n(i,j,k)*vv(i,j,k)
          rhow_n(i,j,k) = rho_n(i,j,k)*ww(i,j,k)
          rhoe_n(i,j,k) = rho_n(i,j,k)*ecalc_tro(Tmp(i,j,k),rho_n(i,j,k))



          !rho_n(i,j,k) = rho_ref!rho_n(i,ny-1,k)
          !prs(i,j,k) = p_ref!prs(i,ny-1,k)
          !Tmp(i,j,k) = T_ref!Tmp(i,ny-1,k)
          !c_(i,j,k) = c_ref!Tmp(i,ny-1,k)
          !rhou_n(i,j,k) = rho_ref*u_ref!rhou_n(i,ny-1,k)
          !rhov_n(i,j,k) = 0!rhov_n(i,ny-1,k)
          !rhoe_n(i,j,k) = rho_ref*ecalc_tro(T_ref,rho_ref)!rhoe_n(i,ny-1,k)
          !uu(i,j,k) = u_ref
          !vv(i,j,k) = 0
          !ww(i,j,k) = 0

          ! Neumann
          rho_n(i,j,k) = rho_n(i,j+1,k)
          prs(i,j,k) = prs(i,j+1,k)
          Tmp(i,j,k) = Tmp(i,j+1,k)
          rhou_n(i,j,k) = rhou_n(i,j+1,k)
          rhov_n(i,j,k) = rhov_n(i,j+1,k)
          rhoe_n(i,j,k) = rhoe_n(i,j+1,k)
          uu(i,j,k) = uu(i,j+1,k) 
          vv(i,j,k) = vv(i,j+1,k) 
          ww(i,j,k) = ww(i,j+1,k) 

          ! compute the rest
          visc(i,j,k) = viscosity_law(Tmp(i,j,k),rho_n(i,j,k))*diffscale
          cok(i,j,k)  = thconductivity(visc(i,j,k),Tmp(i,j,k),rho_n(i,j,k))
          
          Krho(i,j,k) = 0.0_wp
          Krhou(i,j,k) = 0.0_wp
          Krhov(i,j,k) = 0.0_wp
          Krhow(i,j,k) = 0.0_wp
          Krhoe(i,j,k) = 0.0_wp
  
       enddo
    enddo
  
  end subroutine bc_outlet_pfg_jmin

  subroutine bc_outlet_pfg_jmax
  !==============================================================================
    !> Outlet boundary condition at jmax for PFG eos
    !> Use interior node information
  !==============================================================================

    ! Subscripts:
    ! ref: free stream values, prescribed in param.ini
    ! i  : interior of domain, so neighbour node
    ! bc : boundary value

    ! pressure from free stream
    p_ref = pcalc_tro(T_ref,rho_ref)

    ! direction of the velocity at the boundary (WARNING: THIS SHOULD COME FROM THE INPUT FILE)
    nb = (/1.0_wp,0.0_wp,0.0_wp/)

    ! velocity vector at free stream
    vel_ref = u_ref*nb
 
    j=nx 
    ! Update outlet points
    ! ==================  
    do k=1,nz
       do i=ndx_e,nfx_e

          ! boundary normals
          nbc = (/nxn_jmax(i,k),nyn_jmax(i,k),nzn_jmax(i,k)/)
       
          ! velocity vector at the interior nodes
          veli = (/uu(i,j-1,k),vv(i,j-1,k),ww(i,j-1,k)/)

          ! exrapolate R- from the interior nodes
          !!!!!! nbc changed to nb !!!!!!
          Rm = dot_product(veli,nbc)-2.0_wp*c_(i,j-1,k)/gam1

          ! exrapolate R+ from the prescribed free stream
          c_ref = sqrt(c2calc_tro(T_ref,rho_ref))
          Rp = dot_product(vel_ref,nbc)+2.0_wp*c_ref/gam1

          ! Compute velocity and speed of sound at boundary
          velbc = (Rp+Rm)/2.0_wp
          c_(i,j,k) = gam1*(Rp-Rm)/4.0_wp
          !print*,c_(i,j,k),velbc,dot_product(veli,nb)

          ! Extrapolate entropy from interior: sbc = si
          s = scalc_tro(Tmp(i,j-1,k),rho_n(i,j-1,k))

          ! Compute density, pressure and temperature based on s and c
          rho_n(i,j,k) = (c_(i,j,k)**2/(gam*s))**(1.0_wp/gam1)
          prs(i,j,k) = rho_n(i,j,k)*c_(i,j,k)**2/gam
          Tmp(i,j,k) = tcalc_pro(prs(i,j,k),rho_n(i,j,k),Tmp(i,j-1,k))
          !print*,rho_n(i,j,k),Tmp(i,j,k),prs(i,j,k)

          ! Compute velocity and conservative variables
          uu(i,j,k) = veli(1)+velbc*nbc(1)-dot_product(veli,nbc)*nbc(1)
          vv(i,j,k) = veli(2)+velbc*nbc(2)-dot_product(veli,nbc)*nbc(2)
          ww(i,j,k) = veli(3)+velbc*nbc(3)-dot_product(veli,nbc)*nbc(3)
          !uu(i,j,k) = velbc*nb(1)
          !vv(i,j,k) = velbc*nb(2)
          !ww(i,j,k) = velbc*nb(3)
          rhou_n(i,j,k) = rho_n(i,j,k)*uu(i,j,k)
          rhov_n(i,j,k) = rho_n(i,j,k)*vv(i,j,k)
          rhow_n(i,j,k) = rho_n(i,j,k)*ww(i,j,k)
          rhoe_n(i,j,k) = rho_n(i,j,k)*ecalc_tro(Tmp(i,j,k),rho_n(i,j,k))



          !rho_n(i,j,k) = rho_ref!rho_n(i,ny-1,k)
          !prs(i,j,k) = p_ref!prs(i,ny-1,k)
          !Tmp(i,j,k) = T_ref!Tmp(i,ny-1,k)
          !c_(i,j,k) = c_ref!Tmp(i,ny-1,k)
          !rhou_n(i,j,k) = rho_ref*u_ref!rhou_n(i,ny-1,k)
          !rhov_n(i,j,k) = 0!rhov_n(i,ny-1,k)
          !rhoe_n(i,j,k) = rho_ref*ecalc_tro(T_ref,rho_ref)!rhoe_n(i,ny-1,k)
          !uu(i,j,k) = u_ref
          !vv(i,j,k) = 0
          !ww(i,j,k) = 0

          ! Neumann
          rho_n(i,j,k) = rho_n(i,j-1,k)
          prs(i,j,k) = prs(i,j-1,k)
          Tmp(i,j,k) = Tmp(i,j-1,k)
          rhou_n(i,j,k) = rhou_n(i,j-1,k)
          rhov_n(i,j,k) = rhov_n(i,j-1,k)
          rhoe_n(i,j,k) = rhoe_n(i,j-1,k)
          uu(i,j,k) = uu(i,j-1,k) 
          vv(i,j,k) = vv(i,j-1,k) 
          ww(i,j,k) = ww(i,j-1,k) 

          ! compute the rest
          visc(i,j,k) = viscosity_law(Tmp(i,j,k),rho_n(i,j,k))*diffscale
          cok(i,j,k)  = thconductivity(visc(i,j,k),Tmp(i,j,k),rho_n(i,j,k))
          
          Krho(i,j,k) = 0.0_wp
          Krhou(i,j,k) = 0.0_wp
          Krhov(i,j,k) = 0.0_wp
          Krhow(i,j,k) = 0.0_wp
          Krhoe(i,j,k) = 0.0_wp

       enddo
    enddo
  
  end subroutine bc_outlet_pfg_jmax

end module mod_inlet_outlet
