!=================================================================================
submodule (mod_primitives) smod_primitives_mah
!=================================================================================
  !> author: XG
  !> date: January 2024
  !> Computation of primitive variables and viscosity for Martin & Hou EoS
!=================================================================================

contains

  !================================================================================
  module subroutine primitives_mah
  !================================================================================
    !> Computation of primitive variables from conservative ones
    !> - fully inlined version for Martin-Hou EoS -
    !> [applied on rho_n,rhou_n,... for inviscid fluxes]
    !> Nota: check saturation curve not activated
  !================================================================================
    use mod_ineos_mah ! for: constants of MAH model
    use warnstop      ! for: non-convergent Newton
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k,l,icont !,i_sat
    real(wp) :: rro,v,rroe
    real(wp) :: vbm,aaa1,cc1,cc2,T,T1,ekTr,fn,der_fn,err
    ! ----------------------------------------------------------------------------

    ! Compute primitive variables from conservative ones
    ! ==================================================
    icont= 0

    ! for all points + ghost cells
    do k=ndzt,nfzt
       do j=ndyt,nfyt
          !!!dir$ simd
          do i=ndxt,nfxt

             ! Density
             ! -------
             rro=rho_n(i,j,k)

             ! Specific volume
             ! ---------------
             v=1.0_wp/rro

             ! Velocity component
             ! ------------------
             uu(i,j,k)=rhou_n(i,j,k)*v
             vv(i,j,k)=rhov_n(i,j,k)*v
             ww(i,j,k)=rhow_n(i,j,k)*v

             ! total rhoe -> rhoe - kinetic energy
             ! -----------------------------------
             rroe=rhoe_n(i,j,k)-0.5_wp*rro*(uu(i,j,k)**2+vv(i,j,k)**2+ww(i,j,k)**2)

             ! Temperature [inlining instead of calling tcalc_roero]
             ! -----------------------------------------------------
             !! Tmp(i,j,k)=tcalc_roero(bid,rho_n(i,j,k),Tmp(i,j,k))

             vbm = v - bmah

             cc1 = (c2mh/vbm + c3mh/(2.0_wp*vbm**2) + c5mh/(4.0_wp*vbm**4))*kmh/tc
             cc2 = (c2mh/vbm + c3mh/(2.0_wp*vbm**2) + c5mh/(4.0_wp*vbm**4))
             aaa1= (a2mh/vbm + a3mh/(2.0_wp*vbm**2) + a4mh/(3.0_wp*vbm**3)) - rroe*v

             ! initial guess
             T=Tmp(i,j,k)

             ! Newton's loop
             tloop: do l=1,100

!!$                ! check that T does not correspond to a point in the saturation curve
!!$                if (T.lt.tc) then
!!$                   i_sat=nint((v-v_sat(1))/deltav_sat)
!!$                   if (i_sat.lt.1) then
!!$                      write(78,'(4(i0,1X),3(g0,1X))') ntime,i+coord(1)*nx,j+coord(2)*ny,k+coord(3)*nz &
!!$                           ,T/tc,roc/ro,roe*v
!!$                      ! correct with an average of neighboring points
!!$                      ro(i,j,k) =0.5_wp*( ro(i-1,j,k)+ ro(i+1,j,k))
!!$                      rou(i,j,k)=0.5_wp*(rou(i-1,j,k)+rou(i+1,j,k))
!!$                      rov(i,j,k)=0.5_wp*(rov(i-1,j,k)+rov(i+1,j,k))
!!$                      row(i,j,k)=0.5_wp*(row(i-1,j,k)+row(i+1,j,k))
!!$                      roe(i,j,k)=0.5_wp*(roe(i-1,j,k)+roe(i+1,j,k))
!!$                      Tmp(i,j,k)=0.5_wp*(Tmp(i-1,j,k)+Tmp(i+1,j,k))
!!$                      T = Tmp(i,j,k)
!!$                      exit tloop
!!$                   elseif (i_sat.gt.440) then
!!$                      write(79,'(4(i0,1X),3(g0,1X))') ntime,i+coord(1)*nx,j+coord(2)*ny,k+coord(3)*nz &
!!$                           ,T/tc,roc/ro,roe*v
!!$                   else
!!$                      if (T.lt.T_sat(i_sat)) then
!!$                         write(77,'(4(i0,1X),4(g0,1X))') ntime,i+coord(1)*nx,j+coord(2)*ny,k+coord(3)*nz &
!!$                              ,T/tc,T_sat(i_sat)/tc,roc/ro,roe*v
!!$                         T = T_sat(i_sat) + 0.05_wp ! correct T value
!!$                         Tmp(i,j,k)= T
!!$                         roe(i,j,k)=ro*(ecalc_tro(T,ro)+0.5_wp*(uu(i,j,k)**2+vv(i,j,k)**2+ww(i,j,k)**2))
!!$                         icont= 1
!!$                         exit tloop
!!$                      endif
!!$                   endif
!!$                endif

                ekTr= exp(-kmh*T/tc)
                ! function
                fn= aaa1 + cc1*T*ekTr + cc2*ekTr + cvinf*tc/(nexp+1.0_wp)*(abs(T/tc))**(nexp+1.0_wp)
                ! its derivative
                der_fn= (cc1 - kmh/tc*cc2)*ekTr - kmh/Tc*cc1*T*ekTr + cvinf*(abs(T/tc))**nexp

                ! Newton's update
                T1= T - fn/der_fn

                err= abs(T1-T)/T
                if (err.le.tol) then
                   Tmp(i,j,k)= T1
                   T= T1
                   icont= 1
                   exit tloop
                endif
                T= T1
             enddo tloop

             ! Pressure [inlining instead of calling pcalc_Tro]
             ! ------------------------------------------------
             !! prs(i,j,k) = pcalc_tro(Tmp(i,j,k),rho_n(i,j,k))

             ekTr= exp(-kmh*T/tc)

             prs(i,j,k)=                        rg*T/vbm    &
                  + (a2mh + b2mh*T + c2mh*ekTr)/vbm**2 &
                  + (a3mh + b3mh*T + c3mh*ekTr)/vbm**3 &
                  + (a4mh                     )/vbm**4 &
                  + (       b5mh*T + c5mh*ekTr)/vbm**5
          enddo
       enddo
    enddo

    if (icont==0) then
       call mpistop('Error in primitives mah:',1)
    endif

  end subroutine primitives_mah

  !================================================================================
  module subroutine primitives_visc_mah
  !================================================================================
    !> 1/ computation of primitive variables from conservative ones
    !> 2/ computation of thermo-physical properties (visc. & therm. cond.)
    !> 3/ computation of sound speed
    !> - fully inlined version for Martin-Hou EoS -
    !> [applied on rho,rhou,... for viscous fluxes]
    !> Nota: check saturation curve not activated
    !> ATTENTION IL FAUT ENLEVER LES TESTS !!
  !================================================================================
    use mod_constant  ! for: diffscale
    use mod_tranprop  ! for: visc. & therm. cond. parameters
    use mod_ineos_mah ! for: constants of MAH model
    use warnstop      ! for: debug mode & non-convergent Newton
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k,l,icont !,i_sat
    real(wp) :: rro,v,rroe
    real(wp) :: vbm,aaa1,cc1,cc2,T,T1,ekTr,fn,der_fn,err
    real(wp) :: cv,cvres,dpdv,dpdT
    real(wp) :: eta,t_st,omegav,Yv,alfav,betav,zetav,psiv
    real(wp) :: G1,G2,eta_0,eta_k,eta_p,H1,H2,lambda_0,lambda_k,lambda_p
    ! ----------------------------------------------------------------------------

    ! Compute primitive variables, viscosity, thermal cond. and sound speed
    ! =====================================================================
    icont=0

    ! for all points + ghost cells
    do k=ndzt,nfzt
       do j=ndyt,nfyt
          do i=ndxt,nfxt

             ! Density
             ! -------
             rro=rho(i,j,k)

             ! Specific volume
             ! ---------------
             v=1.0_wp/rro

             ! Velocity component
             ! ------------------
             uu(i,j,k)=rhou(i,j,k)*v
             vv(i,j,k)=rhov(i,j,k)*v
             ww(i,j,k)=rhow(i,j,k)*v

             ! total rhoe -> rhoe - kinetic energy
             ! -----------------------------------
             rroe=rhoe(i,j,k)-0.5_wp*rro*(uu(i,j,k)**2+vv(i,j,k)**2+ww(i,j,k)**2)

             ! Temperature [inlining instead of calling tcalc_roero]
             ! -----------------------------------------------------
             !! Tmp(i,j,k)=tcalc_roero(bid,ro(i,j,k),Tmp(i,j,k))

             vbm = v - bmah

             cc1 = (c2mh/vbm + c3mh/(2.0_wp*vbm**2) + c5mh/(4.0_wp*vbm**4))*kmh/tc
             cc2 = (c2mh/vbm + c3mh/(2.0_wp*vbm**2) + c5mh/(4.0_wp*vbm**4))
             aaa1= (a2mh/vbm + a3mh/(2.0_wp*vbm**2) + a4mh/(3.0_wp*vbm**3)) - rroe*v

             ! initial guess
             T=Tmp(i,j,k)

             ! Newton's loop
             tloop: do l=1,100

!!$                ! check that T does not correspond to a point in the saturation curve
!!$                if (T.lt.tc) then
!!$                   i_sat=nint((v-v_sat(1))/deltav_sat)
!!$                   if (i_sat.lt.1) then
!!$                      write(78,'(4(i0,1X),3(g0,1X))') ntime,i+coord(1)*nx,j+coord(2)*ny,k+coord(3)*nz &
!!$                           ,T/tc,roc/ro,roe*v
!!$                      ! correct with an average of neighboring points
!!$                      ro(i,j,k) =0.5_wp*( ro(i-1,j,k)+ ro(i+1,j,k))
!!$                      rou(i,j,k)=0.5_wp*(rou(i-1,j,k)+rou(i+1,j,k))
!!$                      rov(i,j,k)=0.5_wp*(rov(i-1,j,k)+rov(i+1,j,k))
!!$                      row(i,j,k)=0.5_wp*(row(i-1,j,k)+row(i+1,j,k))
!!$                      roe(i,j,k)=0.5_wp*(roe(i-1,j,k)+roe(i+1,j,k))
!!$                      Tmp(i,j,k)=0.5_wp*(Tmp(i-1,j,k)+Tmp(i+1,j,k))
!!$                      T = Tmp(i,j,k)
!!$                      exit tloop
!!$                   elseif (i_sat.gt.440) then
!!$                      write(79,'(4(i0,1X),3(g0,1X))') ntime,i+coord(1)*nx,j+coord(2)*ny,k+coord(3)*nz &
!!$                           ,T/tc,roc/ro,roe*v
!!$                   else
!!$                      if (T.lt.T_sat(i_sat)) then
!!$                         write(77,'(4(i0,1X),4(g0,1X))') ntime,i+coord(1)*nx,j+coord(2)*ny,k+coord(3)*nz &
!!$                              ,T/tc,T_sat(i_sat)/tc,roc/ro,roe*v
!!$                         T = T_sat(i_sat) + 0.05_wp ! correct T value
!!$                         Tmp(i,j,k)= T
!!$                         roe(i,j,k)=ro*(ecalc_tro(T,ro)+0.5_wp*(uu(i,j,k)**2+vv(i,j,k)**2+ww(i,j,k)**2))
!!$                         icont= 1
!!$                         exit tloop
!!$                      endif
!!$                   endif
!!$                endif

                ekTr= exp(-kmh*T/tc)
                ! function
                fn= aaa1 + cc1*T*ekTr + cc2*ekTr + cvinf*tc/(nexp+1.0_wp)*(abs(T/tc))**(nexp+1.0_wp)
                ! its derivative
                der_fn= (cc1 - kmh/tc*cc2)*ekTr - kmh/Tc*cc1*T*ekTr + cvinf*(abs(T/tc))**nexp

                ! Newton's update
                T1= T - fn/der_fn

                err= abs(T1-T)/T
                if (err.le.tol) then
                   Tmp(i,j,k)= T1
                   T= T1
                   icont= 1
                   exit tloop
                endif
                T= T1
             enddo tloop

             ! Pressure [inlining instead of calling pcalc_Tro]
             ! ------------------------------------------------
             !! prs(i,j,k) = pcalc_tro(Tmp(i,j,k),ro(i,j,k))

             ekTr= exp(-kmh*T/tc)

             prs(i,j,k)=                        rg*T/vbm    &
                  + (a2mh + b2mh*T + c2mh*ekTr)/vbm**2 &
                  + (a3mh + b3mh*T + c3mh*ekTr)/vbm**3 &
                  + (a4mh                     )/vbm**4 &
                  + (       b5mh*T + c5mh*ekTr)/vbm**5

             ! Sound speed [inlining instead of calling c2calc_Tro]
             ! ----------------------------------------------------
             ! [used to compute spectral radius, stats, characteristics]
             !! c_(i,j,k) = sqrt(c2calc_tro(Tmp(i,j,k),ro(i,j,k)))

             cvres= c2mh/vbm + c3mh/(2.0_wp*vbm**2) + c5mh/(4.0_wp*vbm**4)
             cv= cvinf*abs(T/Tc)**nexp - T*kmh**2/Tc**2*ekTr*cvres
             dpdv= -                             rg*T/vbm**2 &
                  - 2.0_wp*(a2mh + b2mh*T + c2mh*ekTr)/vbm**3 &
                  - 3.0_wp*(a3mh + b3mh*T + c3mh*ekTr)/vbm**4 &
                  - 4.0_wp*(a4mh                     )/vbm**5 &
                  - 5.0_wp*(       b5mh*T + c5mh*ekTr)/vbm**6

             dpdT=                        rg/vbm    &
                  + (b2mh - c2mh*kmh/Tc*ekTr)/vbm**2 &
                  + (b3mh - c3mh*kmh/Tc*ekTr)/vbm**3 &
                  + (b5mh - c5mh*kmh/Tc*ekTr)/vbm**5

             c_(i,j,k)= sqrt(v**2*(T/cv*dpdT**2 - dpdv))

             ! Viscosity from Chung-Lee-Starling's law [inlining instead of viscosity_law]
             ! ---------------------------------------------------------------------------
             !! visc(i,j,k) = viscosity_law(Tmp(i,j,k),ro(i,j,k))

             t_st= 1.2593_wp*T/tc
             omegav= avisc/t_st**(bvisc) + cvisc/exp(dvisc*t_st)   &
                  + evisc/exp(fvisc*t_st) &
                  + gvisc*t_st**(bvisc)*sin(svisc*t_st**(wvisc)-hvisc)

             eta_0= k_chung*sqrt(T)/omegav

             ! Yv = rro*vmolc/6.0_wp
             Yv= rro/(6.0_wp*roc)

             G1= (1.0_wp - 0.5_wp*Yv)/(1.0_wp-Yv)**3
             G2= (AA1(1)/Yv*(1.0_wp-exp(-AA1(4)*Yv)) + AA1(2)*G1*exp(AA1(5)*Yv) + AA1(3)*G1) / &
                  (AA1(1)*AA1(4) + AA1(2) + AA1(3))
             ! Dilute-Gas component
             eta_k= eta_0 * (1.0_wp/G2 + AA1(6)*Yv)

             ! Dense-Gas component
             eta_p= (36.344e-6_wp * sqrt(pmol*tc)/(vmolc)**(2.0_wp/3.0_wp)) * &
                  AA1(7)*Yv**2*G2*exp(AA1(8) + AA1(9)/t_st + AA1(10)/t_st**2)

             eta= 0.1_wp*(eta_k+eta_p)

             ! Result is scaled (*diffscale) to respect Reynolds analogy
             visc(i,j,k) = eta*diffscale

             ! Thermal conductivity from Chung-Lee-Starling's law [inlining instead of thconductivity]
             ! ---------------------------------------------------------------------------------------
             !! cok(i,j,k) = thconductivity(visc(i,j,k),Tmp(i,j,k),ro(i,j,k))

             alfav= cv/rg - 1.5_wp
             betav= 0.7862_wp - 0.7109_wp*om + 1.3168_wp*om**2
             zetav= 2.0_wp + 10.5_wp*(T/tc)**2
             psiv = 1.0_wp + alfav*(0.215_wp + 0.28288_wp*alfav - 1.061_wp*betav &
                  + 0.26665_wp*zetav)/(0.6366_wp+betav*zetav+1.061_wp*alfav*betav)

             lambda_0= 7.452_wp*eta_0/pmol*psiv

             H1= (1.0_wp - 0.5_wp*Yv)/(1.0_wp-Yv)**3
             H2= (BB1(1)/Yv*(1.0_wp-exp(-BB1(4)*Yv)) + BB1(2)*G1*exp(BB1(5)*Yv) + BB1(3)*G1) / &
                  (BB1(1)*BB1(4) + BB1(2) + BB1(3))

             ! Dilute-Gas component
             lambda_k= lambda_0 * (1.0_wp/H2 + BB1(6)*Yv)

             ! Dense-Gas component
             lambda_p= (3.039e-4_wp * sqrt(T/pmol)/(vmolc)**(2.0_wp/3.0_wp)) * BB1(7)*Yv**2*H2

             ! Result is scaled (*diffscale) to respect Reynolds analogy
             cok(i,j,k)= (lambda_k+lambda_p)*418.68_wp*diffscale

          enddo
       enddo
    enddo

    if (icont.eq.0) then
       call mpistop('Error in primitives mah:',1)
    endif

  end subroutine primitives_visc_mah

end submodule smod_primitives_mah
