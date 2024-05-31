!===============================================================================
module mod_turb_model_length_scale
!================================================================================
  !> author: OY
  !> date: April 2023
  !> Module to compute turbulence model length scale
  !> This length scale determines the simulation frameworks.
  !> The options are: RANS, DES97 (the initial version of Detached-Eddy Simulation),
  !                   DDES (Delayed DES), IDDES (Improved DDES),
  !                   DDES-SLA (DDES with Shear-Layer-Adapted (SLA) subgrid length scale),
  !                   IDDES-SLA (IDDES with SLA)
  !> Note that all DES options are based on S-A model equation.
  !> ATTENTION: In all cases except IDDES, the boundary layer meshing can be done using
  !             a stretching ratio of 1.2 as in the classical RANS meshing.
  !             In IDDES, use a stretching ratio of 1.14 at most, in WMLES regions!!
  ! NOTE-1: I propose to change the name of wall distance from "d" to "d_wall".
  ! NOTE-2: for now, d is a 2D array. so make changes accordingly after 3D version
  ! NOTE-3: The S-A model constants' values can be given in mod_rans.
  ! NOTE-4: IN SLA, the position vectors of vertices, r_n, can be put in a separate subroutine.
!================================================================================
  use mod_flow
  use mod_rans
  use mod_wall_dist
  implicit none
  !------------------------------------------------------------------------------
  real(wp):: C_des = 0.65_wp ! DES constant
  real(wp):: eps = 1e-20_wp ! to avoid zero division
  ! real(wp), dimension(:,:,:), allocatable :: lengthscale,delta_max
  !------------------------------------------------------------------------------

contains

  !==============================================================================
  subroutine RANS_length_scale
  !==============================================================================
    !> Apply RANS approach to obtain turbulent flow field
    !> Model all eddy scales
    !> l_RANS = d_wall where d_wall is the nearest wall distance
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer:: i,j,k
    ! ---------------------------------------------------------------------------

    do k=1,nz
       do j=1,ny
          do i=1,nx
             lengthscale(i,j,k) = d(i,j,k)
          enddo
       enddo
    enddo

  end subroutine RANS_length_scale

  !==============================================================================
  subroutine DES97_length_scale
  !==============================================================================
    !> Apply DES97 approach to obtain turbulent flow field
    !> Model eddies inside the ATTACHED boundary layer (RANS mode), and resolve the rest (LES mode)
    !> l_DES97 = min(d_wall,C_des*delta_max) where delta_max is the max cell dimension
    !> ATTENTION: Generate b/l grid considering delta_max >> b/l thickness
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer:: i,j,k
    ! ---------------------------------------------------------------------------

    do k=ndz_s_r,nfz_s_r
       do j=ndy_s_r,nfy_s_r
          do i=ndx_s_r,nfx_s_r
             lengthscale(i,j,k) = min(d(i,j,k),C_des*delta_max(i,j,k))
          enddo
       enddo
    enddo

  end subroutine DES97_length_scale

  !==============================================================================
  subroutine DDES_length_scale
  !==============================================================================
    !> Apply DDES approach to obtain turbulent flow field
    !> Model eddies inside the ATTACHED boundary layer (RANS mode), and resolve the rest (LES mode)
    !> This ensures to keep RANS mode inside b/l.
    !> l_DDES = d_wall-f_d*max(0,d_wall-Psi*C_des*delta_max)
    !> f_d is a delaying function, ensuring to keep RANS mode in b/l.
    !> Psi is a low Re correction term, preventing an excessive reduction of eddy viscosity in low Re regions.
    !> Note that the trip term, ft2, is ignored.
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer:: i,j,k
    real(wp):: r_d,f_d,Psi,nu,fw,cnu13
    ! ---------------------------------------------------------------------------

    ! Model constants & parameters
    ! ============================
    sig  = 2.0_wp/3.0_wp
    kap  = 0.41_wp
    cb1  = 0.1355_wp
    cb2  = 0.622_wp
    cnu1 = 7.1_wp
    cnu13=cnu1**3
    cw1  = cb1/kap**2+(1.0_wp+cb2)/sig
    fw   = 0.424_wp

    do k=ndz_s_r,nfz_s_r
       do j=ndy_s_r,nfy_s_r
          do i=ndx_s_r,nfx_s_r
             nu = visc(i,j,k)/rho_n(i,j,k)
             khi = nutil(i,j,k)/nu
             fnu1 = khi**3/(khi**3+cnu13)
             fnu2 = 1.0_wp-khi/(1.0_wp+khi*fnu1)
             r_d = (nutil(i,j,k)*fnu1+nu) / ( (kap*d(i,j,k))**2* &
                   max( sqrt(dux(i,j,k)*dux(i,j,k) + duy(i,j,k)*duy(i,j,k) + duz(i,j,k)*duz(i,j,k) + &
                             dvx(i,j,k)*dvx(i,j,k) + dvy(i,j,k)*dvy(i,j,k) + dvz(i,j,k)*dvz(i,j,k) + &
                             dwx(i,j,k)*dwx(i,j,k) + dwy(i,j,k)*dwy(i,j,k) + dwz(i,j,k)*dwz(i,j,k)),eps ) )
             Psi = sqrt(min(100.0_wp,(1.0_wp-cb1*fnu2/(cw1*kap**2*fw))/fnu1))
             f_d = 1.0_wp - tanh((8.0_wp*r_d)**3)
             lengthscale(i,j,k) = d(i,j,k)-f_d*max(0.0_wp,d(i,j,k)-Psi*C_des*delta_max(i,j,k))
          enddo
       enddo
    enddo

  end subroutine DDES_length_scale

  !==============================================================================
  subroutine DDES_SLA_length_scale
  !==============================================================================
    !> Apply DDES-SLA approach to obtain turbulent flow field
    !> Model eddies inside the ATTACHED boundary layer (RANS mode), and resolve the rest (LES mode)
    !> The SLA subgrid length scale (delta_sla) releases K-H instability waves, unlike delta_max.
    !> delta_sla accelerates the transition from RANS to LES mode in case of shear layer separation.
    !> After obtaining delta_sla (although not simple), simply replace delta_max with delta_sla.
    !> l_DDES-SLA = d_wall-f_d*max(0,d_wall-Psi*C_des*delta_sla)
  !==============================================================================
    use mod_mpi_part
    use mod_grid
    implicit none
    ! ---------------------------------------------------------------------------
    integer:: i,j,k,m,n
    real(wp):: r_d,f_d,Psi,nu,fw,delta_sla,cnu13,sqrt6
    ! real(wp):: w1,w2,w3,e11,e12,e13,e22,e23,e33                ! 1st option
    real(wp):: vort(3),strain(3,3),strain2(3,3),str_vort(3),str_vort_vort(3), &
               F_kh,VTM_averaged,fkh_maxmin_a2a1,n_vort(3),vort_magn, &
               r_n(8,3),oneoversqrt3,Inm(3),rnm(3)                   ! 2nd option
    real(wp):: F_kh_min=0.1_wp,F_kh_max=1.0_wp,a1=0.15_wp,a2=0.3_wp
    real(wp), dimension(:,:,:), allocatable :: VTM
    ! ---------------------------------------------------------------------------

    ! Model constants & parameters
    ! ============================
    sig  = 2.0_wp/3.0_wp
    kap  = 0.41_wp
    cb1  = 0.1355_wp
    cb2  = 0.622_wp
    cnu1 = 7.1_wp
    cnu13=cnu1**3
    cw1  = cb1/kap**2+(1.0_wp+cb2)/sig
    fw   = 0.424_wp
    sqrt6 = sqrt(6.0_wp)
    fkh_maxmin_a2a1 = (F_kh_max-F_kh_min)/(a2-a1)
    oneoversqrt3 = 1.0_wp/sqrt(3.0_wp)

    ! compute "Vorticity Tilting Measure" (VTM) and store it
    ! VTM is to detect initial regions of K-H instability
    allocate( VTM(0:nx+1,0:ny+1,0:nz+1) )
    do k=1,nz
       do j=1,ny
          do i=1,nx
!              ! FIRST OPTION (check which option is faster)
!              ! vorticity vector components
!              w1 = dwy(i,j,k)-dvz(i,j,k)
!              w2 = duz(i,j,k)-dwx(i,j,k)
!              w3 = dvx(i,j,k)-duy(i,j,k)
!
!              ! strain rate tensor components (only 6 independent ones)
!              e11 = dux(i,j,k)
!              e12 = 0.5_wp*(duy(i,j,k)+dvx(i,j,k))
!              e13 = 0.5_wp*(duz(i,j,k)+dwx(i,j,k))
!              e22 = dvy(i,j,k)
!              e23 = 0.5_wp*(dvz(i,j,k)+dwy(i,j,k))
!              e33 = dwz(i,j,k)
!
!              VTM(i,j,k) = sqrt6*sqrt((w2*(e11*w1 + e12*w2 + e13*w3) - w1*(e12*w1 + e22*w2 + e23*w3))**2 &
!                                    + (w3*(e11*w1 + e12*w2 + e13*w3) - w1*(e13*w1 + e23*w2 + e33*w3))**2 &
!                                    + (w3*(e12*w1 + e22*w2 + e23*w3) - w2*(e13*w1 + e23*w2 + e33*w3))**2) &
!                         / ((w1**2 + w2**2 + w3**2)*sqrt(3.0_wp(e11**2 + e22**2 + e33**2 + &
!                           2.0_wp*(e12**2 + e13**2 + e23**2)) - (e11 + e22 + e33)**2))

             ! SECOND OPTION (check which option is faster)
             ! vorticity vector
             vort(1) = dwy(i,j,k)-dvz(i,j,k)
             vort(2) = duz(i,j,k)-dwx(i,j,k)
             vort(3) = dvx(i,j,k)-duy(i,j,k)

             ! strain rate tensor
             strain(1,1) = dux(i,j,k)
             strain(1,2) = 0.5_wp*(duy(i,j,k)+dvx(i,j,k))
             strain(1,3) = 0.5_wp*(duz(i,j,k)+dwx(i,j,k))
             strain(2,1) = strain(1,2)
             strain(2,2) = dvy(i,j,k)
             strain(2,3) = 0.5_wp*(dvz(i,j,k)+dwy(i,j,k))
             strain(3,1) = strain(1,3)
             strain(3,2) = strain(2,3)
             strain(3,3) = dwz(i,j,k)

             ! some tensor and vector computations to be used in VTM
             str_vort = matmul(strain,vort)
             str_vort_vort(1) = str_vort(2)*vort(3)-str_vort(3)*vort(2)
             str_vort_vort(2) = str_vort(3)*vort(1)-str_vort(1)*vort(3)
             str_vort_vort(3) = str_vort(1)*vort(2)-str_vort(2)*vort(1)
             strain2 = matmul(strain,strain) ! square of a strain tensor, is it ????

             VTM(i,j,k) = sqrt6*sqrt(dot_product(str_vort_vort,str_vort_vort)) / &
                          max( dot_product(vort,vort)*sqrt(3.0_wp*(strain2(1,1)+strain2(2,2)+strain2(3,3))- &
                          (strain(1,1)+strain(2,2)+strain(3,3))**2),eps )
          enddo
       enddo
    enddo
    ! extrapolate VTM for one ghost node (required for averaging in the next step)
    VTM(0,:,:) = VTM(1,:,:);  VTM(nx+1,:,:) = VTM(nx,:,:)
    VTM(:,0,:) = VTM(:,1,:);  VTM(:,ny+1,:) = VTM(:,ny,:)
    VTM(:,:,0) = VTM(:,:,1);  VTM(:,:,nz+1) = VTM(:,:,nz)

    ! compute the length scale
    do k=ndz_s_r,nfz_s_r
       do j=ndy_s_r,nfy_s_r
          do i=ndx_s_r,nfx_s_r
             ! functions required for DDES
             nu = visc(i,j,k)/rho_n(i,j,k)
             khi = nutil(i,j,k)/nu
             fnu1 = khi**3/(khi**3+cnu13)
             fnu2 = 1.0_wp-khi/(1.0_wp+khi*fnu1)

             r_d = (nutil(i,j,k)*fnu1+nu) / ( (kap*d(i,j,k))**2* &
                   max( sqrt(dux(i,j,k)*dux(i,j,k) + duy(i,j,k)*duy(i,j,k) + duz(i,j,k)*duz(i,j,k) + &
                             dvx(i,j,k)*dvx(i,j,k) + dvy(i,j,k)*dvy(i,j,k) + dvz(i,j,k)*dvz(i,j,k) + &
                             dwx(i,j,k)*dwx(i,j,k) + dwy(i,j,k)*dwy(i,j,k) + dwz(i,j,k)*dwz(i,j,k)),eps ) )
             Psi = sqrt(min(100.0_wp,(1.0_wp-cb1*fnu2/(cw1*kap**2*fw))/fnu1))
             f_d = 1.0_wp - tanh((8.0_wp*r_d)**3)

             ! function to reduce subgrid length scale in K-H regions (like ILES treatment)
             if (f_d < (0.99_wp)) then ! recover DDES approach in case of attached b/l
                F_kh = 1.0_wp
             else
                VTM_averaged = (VTM(i,j,k)+VTM(i+1,j,k)+VTM(i-1,j,k)+VTM(i,j+1,k)+VTM(i,j-1,k)+VTM(i,j,k+1)+VTM(i,j,k-1))/7.0_wp
                ! modification for inviscid regions
                VTM_averaged = VTM_averaged * max(1.0_wp,0.2_wp*nu/(fnu1*max((nutil(i,j,k)-Nutref),1e-6_wp*Nutref)))
                F_kh = max(F_kh_min,min(F_kh_max,F_kh_min+fkh_maxmin_a2a1*(VTM_averaged-a1)))
             endif

             ! This approach uses vorticity dependent subgrid length scale, instead of delta_max
             ! vorticity vector and its magnitude
             vort(1) = dwy(i,j,k)-dvz(i,j,k)
             vort(2) = duz(i,j,k)-dwx(i,j,k)
             vort(3) = dvx(i,j,k)-duy(i,j,k)
             vort_magn = sqrt(dot_product(vort,vort))
             ! unit vector aligned with the vorticity vector
             if (vort_magn < eps) then
                n_vort = oneoversqrt3 ! to make delta_sla=delta_max in irrotational regions
             else
                n_vort = vort/vort_magn
             endif
             ! position vectors of vertices of a hexahedral cell
             if (is_curv) then
                r_n(1,:) = (/xc(i  ,j  ),yc(i  ,j  ),z(k  )/)
                r_n(2,:) = (/xc(i+1,j  ),yc(i+1,j  ),z(k  )/)
                r_n(3,:) = (/xc(i  ,j+1),yc(i  ,j+1),z(k  )/)
                r_n(4,:) = (/xc(i+1,j+1),yc(i+1,j+1),z(k  )/)
                r_n(5,:) = (/xc(i  ,j  ),yc(i  ,j  ),z(k+1)/)
                r_n(6,:) = (/xc(i+1,j  ),yc(i+1,j  ),z(k+1)/)
                r_n(7,:) = (/xc(i  ,j+1),yc(i  ,j+1),z(k+1)/)
                r_n(8,:) = (/xc(i+1,j+1),yc(i+1,j+1),z(k+1)/)
!              elseif (is_curv3) then
!                 r_n(1,:) = (/xc3(i  ,j  ,k  ),yc3(i  ,j  ,k  ),zc3(i  ,j  ,k  )/)
!                 r_n(2,:) = (/xc3(i+1,j  ,k  ),yc3(i+1,j  ,k  ),zc3(i+1,j  ,k  )/)
!                 r_n(3,:) = (/xc3(i  ,j+1,k  ),yc3(i  ,j+1,k  ),zc3(i  ,j+1,k  )/)
!                 r_n(4,:) = (/xc3(i+1,j+1,k  ),yc3(i+1,j+1,k  ),zc3(i+1,j+1,k  )/)
!                 r_n(5,:) = (/xc3(i  ,j  ,k+1),yc3(i  ,j  ,k+1),zc3(i  ,j  ,k+1)/)
!                 r_n(6,:) = (/xc3(i+1,j  ,k+1),yc3(i+1,j  ,k+1),zc3(i+1,j  ,k+1)/)
!                 r_n(7,:) = (/xc3(i  ,j+1,k+1),yc3(i  ,j+1,k+1),zc3(i  ,j+1,k+1)/)
!                 r_n(8,:) = (/xc3(i+1,j+1,k+1),yc3(i+1,j+1,k+1),zc3(i+1,j+1,k+1)/)
             else ! Cartesian
                r_n(1,:) = (/x(i  ),y(j  ),z(k  )/)
                r_n(2,:) = (/x(i+1),y(j  ),z(k  )/)
                r_n(3,:) = (/x(i  ),y(j+1),z(k  )/)
                r_n(4,:) = (/x(i+1),y(j+1),z(k  )/)
                r_n(5,:) = (/x(i  ),y(j  ),z(k+1)/)
                r_n(6,:) = (/x(i+1),y(j  ),z(k+1)/)
                r_n(7,:) = (/x(i  ),y(j+1),z(k+1)/)
                r_n(8,:) = (/x(i+1),y(j+1),z(k+1)/)
             endif
             ! compute the vorticity dependent subgrid length scale
             delta_sla = eps
             do m=1,7
                do n=m+1,8
                   rnm(:) = r_n(n,:)-r_n(m,:)
                   Inm(1) = n_vort(2)*rnm(3)-n_vort(3)*rnm(2)
                   Inm(2) = n_vort(3)*rnm(1)-n_vort(1)*rnm(3)
                   Inm(3) = n_vort(1)*rnm(2)-n_vort(2)*rnm(1)
                   delta_sla = max(delta_sla,dot_product(Inm,Inm))
                enddo
             enddo

             ! compute delta_sla
             delta_sla = sqrt(delta_sla)*oneoversqrt3*F_kh

             ! finally obtain the turb. model length scale
             lengthscale(i,j,k) = d(i,j,k)-f_d*max(0.0_wp,d(i,j,k)-Psi*C_des*delta_sla)
          enddo
       enddo
    enddo

    deallocate(VTM)

  end subroutine DDES_SLA_length_scale

  !==============================================================================
  subroutine IDDES_length_scale
  !==============================================================================
    !> Apply IDDES approach to obtain turbulent flow field
    !> Model eddies inside the ATTACHED boundary layer (RANS mode), and resolve the rest (LES mode)
    !> a combination of DDES and WMLES
    !> It behaves as WMLES starting from loglayer if there is an unsteady inflow turb content (as in reattachment).
    !> Otherwise, it equals to DDES.
    !> l_IDDES = f_d_tilde*(1+f_e)*d_wall+(1-f_d_tilde)*C_des*delta_iddes
    !> f_d_tilde is a function to make switching between DDES and WMLES.
    !> f_e is an elevating func to prevent excessive reduction of modeled Re stress in grey-area (around loglayer).
    !> delta_iddes provides a variation along the wall-normal direction as in the eddy viscosity levels.
    !> delta_iddes gives a reduction of subgrid length, thereby destabilizing the flow under potential instabilities.
  !==============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    integer:: i,j,k
    real(wp):: nu,fw,Psi,cnu13
    real(wp):: c_l,c_t,r_dl,r_dt,alpha,f_e,f_e1,f_e2,f_t,f_l,f_b,f_dt,f_d_til,delta_iddes,C_w
    ! ---------------------------------------------------------------------------

    ! Model constants & parameters
    ! ============================
    sig  = 2.0_wp/3.0_wp
    kap  = 0.41_wp
    cb1  = 0.1355_wp
    cb2  = 0.622_wp
    cnu1 = 7.1_wp
    cnu13=cnu1**3
    cw1  = cb1/kap**2+(1.0_wp+cb2)/sig
    fw   = 0.424_wp
    c_l = 3.55_wp
    c_t = 1.63_wp
    C_w = 0.15_wp

    do k=ndz_s_r,nfz_s_r
       do j=ndy_s_r,nfy_s_r
          do i=ndx_s_r,nfx_s_r
             nu = visc(i,j,k)/rho_n(i,j,k)
             khi = nutil(i,j,k)/nu
             fnu1 = khi**3/(khi**3+cnu13)
             fnu2 = 1.0_wp-khi/(1.0_wp+khi*fnu1)
             Psi = sqrt(min(100.0_wp,(1.0_wp-cb1*fnu2/(cw1*kap**2*fw))/fnu1))

             ! compute blending and elevating functions
             r_dl = nu / ( (kap*d(i,j,k))**2* &
                    max( sqrt(dux(i,j,k)*dux(i,j,k) + duy(i,j,k)*duy(i,j,k) + duz(i,j,k)*duz(i,j,k) + &
                              dvx(i,j,k)*dvx(i,j,k) + dvy(i,j,k)*dvy(i,j,k) + dvz(i,j,k)*dvz(i,j,k) + &
                              dwx(i,j,k)*dwx(i,j,k) + dwy(i,j,k)*dwy(i,j,k) + dwz(i,j,k)*dwz(i,j,k)),eps ) )
             r_dt = (r_dl/nu)*nutil(i,j,k)*fnu1
             f_dt = 1.0_wp - tanh((8.0_wp*r_dt)**3)
             alpha = 0.25_wp-d(i,j,k)/delta_max(i,j,k)
             f_e1 = 2.0_wp*exp(-11.09_wp*(max(0.0_wp,alpha))**2-9.0_wp*(min(0.0_wp,alpha))**2)
             f_l = tanh((c_l**2*r_dl)**8) !ATTENTION this should be tanh((cl**2*r_dl)**10), but 10 causes floating overflow
             f_t = tanh((c_t**2*r_dt)**3)
             f_e2 = 1.0_wp - max(f_l,f_t)
             f_e = max(f_e1-1.0_wp,0.0_wp)*Psi*f_e2
             f_b = min(2.0_wp*exp(-9.0_wp*alpha**2),1.0_wp)
             f_d_til = max(1.0_wp-f_dt,f_b)

             ! compute delta_iddes
             delta_iddes = min(max(C_w*d(i,j,k),C_w*delta_max(i,j,k),h_wn(i,j)),delta_max(i,j,k))

             ! finally obtain the turb. model length scale
             lengthscale(i,j,k) = f_d_til*(1.0_wp+f_e)*d(i,j,k)+(1.0_wp-f_d_til)*C_des*delta_iddes
          enddo
       enddo
    enddo

  end subroutine IDDES_length_scale

  !==============================================================================
  subroutine maximum_cell_dimension
  !==============================================================================
    !> Compute local max cell dimension of each grid node
  !==============================================================================
    use mod_mpi_part
    use mod_grid
    implicit none
    ! ---------------------------------------------------------------------------
    integer:: i,j,k
    real(wp):: dip1,dim1,djp1,djm1,dipm1,djpm1,dkpm1
    ! real(wp):: dkp1,dkm1
    ! ---------------------------------------------------------------------------

    allocate(delta_max(1:nx,1:ny,1:nz))

    if (is_curv) then
       if (is_2D) then
          do j=ndy_d_r,nfy_d_r
             do i=ndx_d_r,nfx_d_r
                dip1 = sqrt( (xc(i+1,j)-xc(i,j))**2 + (yc(i+1,j)-yc(i,j))**2 )
                dim1 = sqrt( (xc(i-1,j)-xc(i,j))**2 + (yc(i-1,j)-yc(i,j))**2 )
                djp1 = sqrt( (xc(i,j+1)-xc(i,j))**2 + (yc(i,j+1)-yc(i,j))**2 )
                djm1 = sqrt( (xc(i,j-1)-xc(i,j))**2 + (yc(i,j-1)-yc(i,j))**2 )
                delta_max(i,j,1) = 0.5_wp*max(dip1+dim1,djp1+djm1)
             enddo
          enddo
       else ! 3D
          do k=ndz_d_r,nfz_d_r
             do j=ndy_d_r,nfy_d_r
                do i=ndx_d_r,nfx_d_r
                   dip1 = sqrt( (xc(i+1,j)-xc(i,j))**2 + (yc(i+1,j)-yc(i,j))**2 )
                   dim1 = sqrt( (xc(i-1,j)-xc(i,j))**2 + (yc(i-1,j)-yc(i,j))**2 )
                   djp1 = sqrt( (xc(i,j+1)-xc(i,j))**2 + (yc(i,j+1)-yc(i,j))**2 )
                   djm1 = sqrt( (xc(i,j-1)-xc(i,j))**2 + (yc(i,j-1)-yc(i,j))**2 )
                   dkpm1 = z(k+1)-z(k-1)
                   delta_max(i,j,k) = 0.5_wp*max(dip1+dim1,djp1+djm1,dkpm1)
                enddo
             enddo
          enddo
       endif

!     elseif (is_curv3) then
!
!        do k=1,nz
!           do j=1,ny
!              do i=1,nx
!                 dip1 = sqrt( (xc3(i+1,j,k)-xc3(i,j,k))**2 + (yc3(i+1,j,k)-yc3(i,j,k))**2 + (zc3(i+1,j,k)-yc3(i,j,k))**2 )
!                 dim1 = sqrt( (xc3(i-1,j,k)-xc3(i,j,k))**2 + (yc3(i-1,j,k)-yc3(i,j,k))**2 + (zc3(i-1,j,k)-yc3(i,j,k))**2 )
!                 djp1 = sqrt( (xc3(i,j+1,k)-xc3(i,j,k))**2 + (yc3(i,j+1,k)-yc3(i,j,k))**2 + (zc3(i,j+1,k)-yc3(i,j,k))**2 )
!                 djm1 = sqrt( (xc3(i,j-1,k)-xc3(i,j,k))**2 + (yc3(i,j-1,k)-yc3(i,j,k))**2 + (zc3(i,j-1,k)-yc3(i,j,k))**2 )
!                 dkp1 = sqrt( (xc3(i,j,k+1)-xc3(i,j,k))**2 + (yc3(i,j,k+1)-yc3(i,j,k))**2 + (zc3(i,j,k+1)-yc3(i,j,k))**2 )
!                 dkm1 = sqrt( (xc3(i,j,k-1)-xc3(i,j,k))**2 + (yc3(i,j,k-1)-yc3(i,j,k))**2 + (zc3(i,j,k-1)-yc3(i,j,k))**2 )
!                 delta_max(i,j,k) = 0.5_wp*max(dip1+dim1,djp1+djm1,dkp1+dkm1)
!              enddo
!           enddo
!        enddo

    else ! Cartesian

       if (is_2D) then
          do j=ndy_d_r,nfy_d_r
             do i=ndx_d_r,nfx_d_r
                dipm1 = x(i+1)-x(i-1)
                djpm1 = y(j+1)-y(j-1)
                delta_max(i,j,1) = 0.5_wp*max(dipm1,djpm1)
             enddo
          enddo
       else ! 3D
          do k=ndz_d_r,nfz_d_r
             do j=ndy_d_r,nfy_d_r
                do i=ndx_d_r,nfx_d_r
                   dipm1 = x(i+1)-x(i-1)
                   djpm1 = y(j+1)-y(j-1)
                   dkpm1 = z(k+1)-z(k-1)
                   delta_max(i,j,k) = 0.5_wp*max(dipm1,djpm1,dkpm1)
                enddo
             enddo
          enddo
       endif

    endif

    ! For boundaries, delta_max put at the interior value
    if (is_boundary(1,1)) delta_max(1 ,ndy_d_r:nfy_d_r,ndz_d_r:nfz_d_r) = delta_max(2   ,ndy_d_r:nfy_d_r,ndz_d_r:nfz_d_r)
    if (is_boundary(1,2)) delta_max(nx,ndy_d_r:nfy_d_r,ndz_d_r:nfz_d_r) = delta_max(nx-1,ndy_d_r:nfy_d_r,ndz_d_r:nfz_d_r)
    if (is_boundary(2,1)) delta_max(ndx_d_r:nfx_d_r,1 ,ndz_d_r:nfz_d_r) = delta_max(ndx_d_r:nfx_d_r,2   ,ndz_d_r:nfz_d_r)
    if (is_boundary(2,2)) delta_max(ndx_d_r:nfx_d_r,ny,ndz_d_r:nfz_d_r) = delta_max(ndx_d_r:nfx_d_r,ny-1,ndz_d_r:nfz_d_r)
    if (.not.is_2D) then
       if (is_boundary(3,1)) delta_max(ndx_d_r:nfx_d_r,ndy_d_r:nfy_d_r,1 ) = delta_max(ndx_d_r:nfx_d_r,ndy_d_r:nfy_d_r,2   )
       if (is_boundary(3,2)) delta_max(ndx_d_r:nfx_d_r,ndy_d_r:nfy_d_r,nz) = delta_max(ndx_d_r:nfx_d_r,ndy_d_r:nfy_d_r,nz-1)
    endif

  end subroutine maximum_cell_dimension

end module mod_turb_model_length_scale
