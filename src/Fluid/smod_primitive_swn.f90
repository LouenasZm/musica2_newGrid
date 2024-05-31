!=================================================================================
submodule (mod_primitives) smod_primitives_swn
!=================================================================================
  !> author: XG
  !> date: January 2024
  !> Computation of primitive variables and viscosity for Span-Wagner EoS
  !> for non-polar compounds
!=================================================================================

contains

  !================================================================================
  module subroutine primitives_swn
  !================================================================================
    !> Computation of primitive variables from conservative ones
    !> - partially inlined version for non-polar Span-Wagner EoS -
    !> [applied on rho_n,rhou_n,... for inviscid fluxes]
    !> /!\ NOT FINALIZED and NOT TESTED
    !> Nota: viscosity law and thermal conductivity law NOT INLINED
  !================================================================================
    use mod_ineos_swn ! for: constants of non-polar SW model
    use warnstop      ! for: non-convergent Newton
    implicit none
    ! ----------------------------------------------------------------------------
    integer  :: i,j,k,l,icont
    real(wp) :: rro,v,rroe
    real(wp) :: ea,tau,del,dedT,func,ecalc
    real(wp) :: Ttent,T,T1,err
    real(wp) :: d1p0_tau1,d2p0_tau2,d1pr_tau1,d2pr_tau2,d1pr_del1
    ! ----------------------------------------------------------------------------

    ! ****************************
    ! /!\ Routine not updated /!\
    ! ****************************
     call mpistop('Routine primitives_swn not updated ! No calculation of the speed sound',0)

    ! Compute primitive variables from conservative ones
    ! ==================================================
    icont=0

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

             ! initial guess
             Ttent = Tmp(i,j,k)          
             T  = Ttent*Tci
             ea = rroe*v/(zc*pcvc)
             rro= rro*roci
             del= rro

             ! Newton's loop
             tloop: do l=1,1000
                tau= 1.0_wp/T

                d1p0_tau1= (eta1-1.0_wp)*(1.0_wp/tau   -1.0_wp) &
                     + eta2/2.0_wp*tc   *(1.0_wp/tau**2-1.0_wp) &
                     + eta3/3.0_wp*tc**2*(1.0_wp/tau**3-1.0_wp) &
                     + eta4/4.0_wp*tc**3*(1.0_wp/tau**4-1.0_wp)

                d2p0_tau2=-(eta1-1.0_wp)/tau**2 + eta2*tc*(-1.0_wp/tau**3) &
                     + eta3*tc**2*(-1.0_wp/tau**4) &
                     + eta4*tc**3*(-1.0_wp/tau**5)

                d1pr_tau1= n01*0.25_wp *del   *tau**(-0.75_wp)               &
                         + n02*1.125_wp*del   *tau**( 0.125_wp)              &
                         + n03*1.5_wp  *del   *tau**( 0.5_wp)                &
                         + n04*1.375_wp*del**2*tau**( 0.375_wp)              &
                         + n05*0.25_wp *del**3*tau**(-0.75_wp)               &
                         + n06*0.875_wp*del**7*tau**(-0.125_wp)              &
                         + n07*0.625_wp*del**2*tau**(-0.375_wp)*exp(-del   ) &
                         + n08*1.75_wp *del**5*tau**( 0.75_wp) *exp(-del   ) &
                         + n09*3.625_wp*del   *tau**( 2.625_wp)*exp(-del**2) &
                         + n10*3.625_wp*del**4*tau**( 2.625_wp)*exp(-del**2) &
                         + n11*14.5_wp *del**3*tau**( 13.5_wp) *exp(-del**3) &
                         + n12*12.0_wp *del**4*tau**( 11.0_wp) *exp(-del**3)

                d2pr_tau2=- n01*0.75_wp *0.25_wp *del   *tau**(-1.75_wp)               &
                          + n02*0.125_wp*1.125_wp*del   *tau**(-0.875_wp)              &
                          + n03*0.5_wp  *1.5_wp  *del   *tau**(-0.5_wp)                &
                          + n04*0.375_wp*1.375_wp*del**2*tau**(-0.625_wp)              &
                          - n05*0.75_wp *0.25_wp *del**3*tau**(-1.75_wp)               &
                          - n06*0.125_wp*0.875_wp*del**7*tau**(-1.125_wp)              &
                          - n07*0.375_wp*0.625_wp*del**2*tau**(-1.375_wp)*exp(-del   ) &
                          + n08*0.75_wp *1.75_wp *del**5*tau**(-0.25_wp) *exp(-del   ) &
                          + n09*2.625_wp*3.625_wp*del   *tau**( 1.625_wp)*exp(-del**2) &
                          + n10*2.625_wp*3.625_wp*del**4*tau**( 1.625_wp)*exp(-del**2) &
                          + n11*13.5_wp *14.5_wp *del**3*tau**( 12.5_wp) *exp(-del**3) &
                          + n12*11.0_wp *12.0_wp *del**4*tau**( 10.0_wp) *exp(-del**3)


                ecalc= 100.0_wp + zci*(d1pr_tau1+d1p0_tau1)
                dedT =          - zci*(d2pr_tau2+d2p0_tau2)/T**2

                ! function
                func= ecalc - ea

                ! Newton's update
                T1= T - func/dedT

                err= abs(T1-T)/T
                if (err.lt.tol) then
                   Tmp(i,j,k)= T1*tc
                   T= T1*tc
                   icont= 1
                   exit tloop
                endif
                T= T1
             enddo tloop

             ! Pressure [inlining instead of calling pcalc_Tro]
             ! ------------------------------------------------
             !! prs(i,j,k) = pcalc_tro(Tmp(i,j,k),rho_n(i,j,k))

             T= Tmp(i,j,k)*Tci
             tau= 1.0_wp/T

             d1pr_del1= n01              *tau**(0.25_wp)                                      &
                      + n02              *tau**(1.125_wp)                                     &
                      + n03              *tau**(1.5_wp)                                       &
                      + n04*2.0_wp*del   *tau**(1.375_wp)                                     &
                      + n05*3.0_wp*del**2*tau**(0.25_wp)                                      &
                      + n06*7.0_wp*del**6*tau**(0.875_wp)                                     &
                      + n07       *del   *tau**(0.625_wp)*exp(-del)   *(2.0_wp-       del   ) &
                      + n08       *del**4*tau**(1.75_wp) *exp(-del)   *(5.0_wp-       del   ) &
                      + n09              *tau**(3.625_wp)*exp(-del**2)*(1.0_wp-2.0_wp*del**2) &
                      + n10       *del**3*tau**(3.625_wp)*exp(-del**2)*(4.0_wp-2.0_wp*del**2) &
                      + n11*3.0_wp*del**2*tau**(14.5_wp) *exp(-del**3)*(1.0_wp-       del**3) &
                      + n12       *del**3*tau**(12.0_wp) *exp(-del**3)*(4.0_wp-3.0_wp*del**3)

             prs(i,j,k)= (del*zci/tau*(1.0_wp+del*d1pr_del1))*pc
          enddo
       enddo
    enddo

    if (icont==0) then
       call mpistop('Error in primitives swn:',1)
    endif

  end subroutine primitives_swn

  !================================================================================
  module subroutine primitives_visc_swn
  !================================================================================
    !> 1/ computation of primitive variables from conservative ones
    !> 2/ computation of thermo-physical properties (visc. & therm. cond.)
    !> 3/ computation of sound speed
    !> - partially inlined version for non-polar Span-Wagner EoS -
    !> [applied on rho,rhou,... for viscous fluxes]
    !> /!\ NOT FINALIZED (speed of sound) and NOT TESTED
  !================================================================================
    use mod_constant  ! for: diffscale
    use mod_tranprop  ! for: visc. & therm. cond. parameters
    use mod_eos_swn   ! for: cvcalc_tro   
    use mod_ineos_swn ! for: constants of non-polar SW model
    use warnstop      ! for: debug mode & non-convergent Newton
    implicit none
    ! ----------------------------------------------------------------------------
    integer  :: i,j,k,l,icont
    real(wp) :: rro,v,rroe
    real(wp) :: ea,tau,del,dedT,func,ecalc
    real(wp) :: Ttent,T,T1,err
    real(wp) :: eta,t_st,omegav,Yv,alfav,betav,zetav,psiv
    real(wp) :: G1,G2,eta_0,eta_k,eta_p,H1,H2,lambda_0,lambda_k,lambda_p
    real(wp) :: d1p0_tau1,d2p0_tau2,d1pr_tau1,d2pr_tau2,d1pr_del1
    ! ----------------------------------------------------------------------------

    ! ****************************
    ! /!\ Routine not updated /!\
    ! ****************************
    call mpistop('Routine primitives_swn not updated ! No calculation of the speed sound',0)

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

             ! initial guess
             Ttent = Tmp(i,j,k)          
             T  = Ttent*Tci
             ea = rroe/rro/(zc*pcvc)
             rro= rro*roci
             del= rro

             ! Newton's loop
             tloop: do l=1,1000
                tau= 1.0_wp/T

                d1p0_tau1= (eta1-1.0_wp)*(1.0_wp/tau   -1.0_wp) &
                          + eta2/2.0_wp*tc   *(1.0_wp/tau**2-1.0_wp) &
                          + eta3/3.0_wp*tc**2*(1.0_wp/tau**3-1.0_wp) &
                          + eta4/4.0_wp*tc**3*(1.0_wp/tau**4-1.0_wp)

                d2p0_tau2=-(eta1-1.0_wp)/tau**2 + eta2*tc*(-1.0_wp/tau**3) &
                          + eta3*tc**2*(-1.0_wp/tau**4) &
                          + eta4*tc**3*(-1.0_wp/tau**5)

                d1pr_tau1= n01*0.25_wp *del   *tau**(-0.75_wp)               &
                         + n02*1.125_wp*del   *tau**( 0.125_wp)              &
                         + n03*1.5_wp  *del   *tau**( 0.5_wp)                &
                         + n04*1.375_wp*del**2*tau**( 0.375_wp)              &
                         + n05*0.25_wp *del**3*tau**(-0.75_wp)               &
                         + n06*0.875_wp*del**7*tau**(-0.125_wp)              &
                         + n07*0.625_wp*del**2*tau**(-0.375_wp)*exp(-del   ) &
                         + n08*1.75_wp *del**5*tau**( 0.75_wp) *exp(-del   ) &
                         + n09*3.625_wp*del   *tau**( 2.625_wp)*exp(-del**2) &
                         + n10*3.625_wp*del**4*tau**( 2.625_wp)*exp(-del**2) &
                         + n11*14.5_wp *del**3*tau**( 13.5_wp) *exp(-del**3) &
                         + n12*12.0_wp *del**4*tau**( 11.0_wp) *exp(-del**3)

                d2pr_tau2=- n01*0.75_wp *0.25_wp *del   *tau**(-1.75_wp)               &
                          + n02*0.125_wp*1.125_wp*del   *tau**(-0.875_wp)              &
                          + n03*0.5_wp  *1.5_wp  *del   *tau**(-0.5_wp)                &
                          + n04*0.375_wp*1.375_wp*del**2*tau**(-0.625_wp)              &
                          - n05*0.75_wp *0.25_wp *del**3*tau**(-1.75_wp)               &
                          - n06*0.125_wp*0.875_wp*del**7*tau**(-1.125_wp)              &
                          - n07*0.375_wp*0.625_wp*del**2*tau**(-1.375_wp)*exp(-del   ) &
                          + n08*0.75_wp *1.75_wp *del**5*tau**(-0.25_wp) *exp(-del   ) &
                          + n09*2.625_wp*3.625_wp*del   *tau**( 1.625_wp)*exp(-del**2) &
                          + n10*2.625_wp*3.625_wp*del**4*tau**( 1.625_wp)*exp(-del**2) &
                          + n11*13.5_wp *14.5_wp *del**3*tau**( 12.5_wp) *exp(-del**3) &
                          + n12*11.0_wp *12.0_wp *del**4*tau**( 10.0_wp) *exp(-del**3)


                ecalc= 100.0_wp + zci*(d1pr_tau1+d1p0_tau1)
                dedT =          - zci*(d2pr_tau2+d2p0_tau2)/T**2

                ! function
                func= ecalc - ea

                ! Newton's update
                T1= T - func/dedT

                err= abs(T1-T)/T
                if (err.lt.tol) then
                   Tmp(i,j,k)= T1*tc
                   T= T1*tc
                   icont= 1
                   exit tloop
                endif
                T= T1
             enddo tloop

             ! Pressure [inlining instead of calling pcalc_Tro]
             ! ------------------------------------------------
             !! prs(i,j,k) = pcalc_tro(Tmp(i,j,k),ro(i,j,k))

             T= Tmp(i,j,k)*Tci
             tau= 1.0_wp/T

             d1pr_del1= n01              *tau**(0.25_wp)                                      &
                      + n02              *tau**(1.125_wp)                                     &
                      + n03              *tau**(1.5_wp)                                       &
                      + n04*2.0_wp*del   *tau**(1.375_wp)                                     &
                      + n05*3.0_wp*del**2*tau**(0.25_wp)                                      &
                      + n06*7.0_wp*del**6*tau**(0.875_wp)                                     &
                      + n07       *del   *tau**(0.625_wp)*exp(-del)   *(2.0_wp-       del   ) &
                      + n08       *del**4*tau**(1.75_wp) *exp(-del)   *(5.0_wp-       del   ) &
                      + n09              *tau**(3.625_wp)*exp(-del**2)*(1.0_wp-2.0_wp*del**2) &
                      + n10       *del**3*tau**(3.625_wp)*exp(-del**2)*(4.0_wp-2.0_wp*del**2) &
                      + n11*3.0_wp*del**2*tau**(14.5_wp) *exp(-del**3)*(1.0_wp-       del**3) &
                      + n12       *del**3*tau**(12.0_wp) *exp(-del**3)*(4.0_wp-3.0_wp*del**3)

             prs(i,j,k)= (del*zci/tau*(1.0_wp+del*d1pr_del1))*pc

             ! Sound speed [inlining instead of calling c2calc_Tro]
             ! ----------------------------------------------------
             ! [used to compute spectral radius, stats, characteristics]
             !! c2(i,j,k) = c2calc_tro(Tmp(i,j,k),ro(i,j,k))

             ! Viscosity from Chung-Lee-Starling's law [inlining instead of viscosity_law]
             ! ---------------------------------------------------------------------------
             !! visc(i,j,k) = viscosity_law(Tmp(i,j,k),ro(i,j,k))

             T= Tmp(i,j,k)
             rro= rho(i,j,k)

             t_st = 1.2593_wp*T/tc
             omegav = avisc/t_st**(bvisc) + cvisc/exp(dvisc*t_st)   &
                  + evisc/exp(fvisc*t_st) &
                  + gvisc*t_st**(bvisc)*sin(svisc*t_st**(wvisc)-hvisc)

             eta_0= k_chung*sqrt(T)/omegav

             ! Yv = ro*vmolc/6.0_wp
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

             alfav= cvcalc_tro_swn(T,rro)/rg - 1.5_wp
             betav= 0.7862_wp - 0.7109_wp*om + 1.3168_wp*om**2
             zetav= 2.0_wp + 10.5_wp*(T/tc)**2
             psiv = 1.0_wp + alfav*(0.215_wp + 0.28288_wp*alfav - 1.061_wp*betav       &
                  + 0.26665_wp*zetav)/(0.6366_wp+betav*zetav+1.061_wp*alfav*betav)
             !
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

    if (icont==0) then
       call mpistop('Error in primitives swn:',1)
    endif

  end subroutine primitives_visc_swn

end submodule smod_primitives_swn
