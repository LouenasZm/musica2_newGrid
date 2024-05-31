!=================================================================================
submodule (mod_primitives) smod_primitives_swp
!=================================================================================
  !> author: XG
  !> date: January 2024
  !> Computation of primitive variables and viscosity for Span-Wagner EoS
  !> for polar compounds
!=================================================================================

contains

  !================================================================================
  module subroutine primitives_swp
  !================================================================================
    !> Computation of primitive variables from conservative ones
    !> - partially inlined version for polar Span-Wagner EoS -
    !> [applied on rho_n,rhou_n,... for inviscid fluxes]
    !> /!\ NOT FINALIZED and NOT TESTED
    !> Nota: viscosity law and thermal conductivity law NOT INLINED
  !================================================================================
    use mod_ineos_swn ! for: constants of polar SW model
    use warnstop      ! for: non-convergent Newton
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k,l,icont
    real(wp) :: rro,v,rroe
    real(wp) :: ea,tau,del,dedT,func,ecalc
    real(wp) :: Ttent,T,T1,err
    ! real(wp) :: eta,lambda,t_st,omegav,Yv,alfav,betav,zetav,psiv
    ! real(wp) :: G1,G2,eta_0,eta_k,eta_p,H1,H2,lambda_0,lambda_k,lambda_p
    real(wp) :: d1p0_tau1,d2p0_tau2,d1pr_tau1,d2pr_tau2,d1pr_del1
    ! ----------------------------------------------------------------------------

    ! ****************************
    ! /!\ Routine not updated /!\
    ! ****************************
    call mpistop('Routine primitives_swp not updated ! No calculation of the speed sound',0)

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
             !! Tmp(i,j,k)=tcalc_roero(bid,ro(i,j,k),Tmp(i,j,k))

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

                d1pr_tau1= n01*0.25_wp *del   *tau**(-0.75_wp )              &
                         + n02*1.25_wp *del   *tau**( 0.25_wp )              &
                         + n03*1.5_wp  *del   *tau**( 0.5_wp  )              &
                         + n04*0.25_wp *del**3*tau**(-0.75_wp )              &
                         + n05*0.875_wp*del**7*tau**(-0.125_wp)              &
                         + n06*2.375_wp*del   *tau**( 1.375_wp)*exp(-del)    &
                         + n07*2.0_wp  *del**2*tau             *exp(-del)    &
                         + n08*2.125_wp*del**5*tau**( 1.125_wp)*exp(-del)    &
                         + n09*3.5_wp  *del   *tau**( 2.5_wp)  *exp(-del**2) &
                         + n10*6.5_wp  *del   *tau**( 5.5_wp)  *exp(-del**2) &
                         + n11*4.75_wp *del**4*tau**( 3.75_wp) *exp(-del**2) &
                         + n12*12.5_wp *del**2*tau**( 11.5_wp) *exp(-del**3)

                d2pr_tau2= n01*0.25_wp *(-0.75_wp) *del   *tau**(-1.75_wp)               &
                         + n02*1.25_wp *( 0.25_wp) *del   *tau**(-0.75_wp)               &
                         + n03*1.5_wp  *( 0.5_wp)  *del   *tau**(-0.5_wp)                &
                         + n04*0.25_wp *(-0.75_wp) *del**3*tau**(-1.75_wp)               &
                         + n05*0.875_wp*(-0.125_wp)*del**7*tau**(-1.125_wp)              &
                         + n06*2.375_wp*( 1.375_wp)*del   *tau**( 0.375_wp)*exp(-del)    &
                         + n07*2.0_wp              *del**2                 *exp(-del)    &
                         + n08*2.125_wp*( 1.125_wp)*del**5*tau**( 0.125_wp)*exp(-del)    &
                         + n09*3.5_wp  *( 2.5_wp)  *del   *tau**( 1.5_wp)  *exp(-del**2) &
                         + n10*6.5_wp  *( 5.5_wp)  *del   *tau**( 4.5_wp)  *exp(-del**2) &
                         + n11*4.75_wp *( 3.75_wp) *del**4*tau**( 2.75_wp) *exp(-del**2) &
                         + n12*12.5_wp *( 11.5_wp) *del**2*tau**( 10.5_wp) *exp(-del**3)


                ecalc= 100.0_wp + zci*(d1pr_tau1+d1p0_tau1)
                dedT =          - zci*(d2pr_tau2+d2p0_tau2)/T**2

                ! function
                func= ecalc - ea

                ! Newton's update
                T1= T - func/dedT

                err = abs(T1-T)/T
                if (err.lt.tol) then
                   Tmp(i,j,k)= T1*tc
                   T= T1*tc
                   icont= 1
                   exit tloop
                endif
                T = T1
             enddo tloop

             ! Pressure [inlining instead of calling pcalc_Tro]
             ! ------------------------------------------------
             !! prs(i,j,k) = pcalc_tro(Tmp(i,j,k),ro(i,j,k))

             T= Tmp(i,j,k)*Tci
             tau= 1.0_wp/T

             d1pr_del1= n01              *tau**(0.25_wp )                                     &
                      + n02              *tau**(1.25_wp )                                     &
                      + n03              *tau**(1.5_wp  )                                     &
                      + n04*3.0_wp*del**2*tau**(0.25_wp )                                     &
                      + n05*7.0_wp*del**6*tau**(0.875_wp)                                     &
                      + n06              *tau**(2.375_wp)*exp(-del   )*(1.0_wp-     del   )   &
                      + n07       *del   *tau**(2       )*exp(-del   )*(2.0_wp-     del   )   &
                      + n08       *del**4*tau**(2.125_wp)*exp(-del   )*(5.0_wp-     del   )   &
                      + n09              *tau**(3.5_wp  )*exp(-del**2)*(1.0_wp-2.0_wp*del**2) &
                      + n10              *tau**(6.5_wp  )*exp(-del**2)*(1.0_wp-2.0_wp*del**2) &
                      + n11       *del**3*tau**(4.75_wp )*exp(-del**2)*(4.0_wp-2.0_wp*del**2) &
                      + n12       *del   *tau**(12.5_wp )*exp(-del**3)*(2.0_wp-3.0_wp*del**3)

             prs(i,j,k)= (del*zci/tau*(1.0_wp+del*d1pr_del1))*pc
          enddo
       enddo
    enddo

    if (icont.eq.0) then
       call mpistop('Error in primitives swp:',1)
    endif

  end subroutine primitives_swp

  !================================================================================
  module subroutine primitives_visc_swp
  !================================================================================
    !> 1/ computation of primitive variables from conservative ones
    !> 2/ computation of thermo-physical properties (visc. & therm. cond.)
    !> 3/ computation of sound speed
    !> - partially inlined version for polar Span-Wagner EoS -
    !> [applied on rho,rhou,... for viscous fluxes]
    !> /!\ NOT FINALIZED (speed of sound) and NOT TESTED
    !> Nota: viscosity law and thermal conductivity law NOT INLINED
  !================================================================================
    use mod_constant  ! for: diffscale
    use mod_tranprop  ! for: visc. & therm. cond. parameters
    use mod_ineos_swp ! for: constants of polar SW model
    use warnstop      ! for: debug mode & non-convergent Newton
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k,l,icont
    real(wp) :: rro,v,rroe
    real(wp) :: ea,tau,del,dedT,func,ecalc
    real(wp) :: Ttent,T,T1,err
    !real(wp) :: eta,lambda,t_st,omegav,Yv,alfav,betav,zetav,psiv
    !real(wp) :: G1,G2,eta_0,eta_k,eta_p,H1,H2,lambda_0,lambda_k,lambda_p
    real(wp) :: d1p0_tau1,d2p0_tau2,d1pr_tau1,d2pr_tau2,d1pr_del1
    ! ----------------------------------------------------------------------------

    ! ****************************
    ! /!\ Routine not updated /!\
    ! ****************************
    call mpistop('Routine primitives_swp not updated ! No calculation of the speed sound',0)

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

                d1pr_tau1= n01*0.25_wp *del   *tau**(-0.75_wp )              &
                         + n02*1.25_wp *del   *tau**( 0.25_wp )              &
                         + n03*1.5_wp  *del   *tau**( 0.5_wp  )              &
                         + n04*0.25_wp *del**3*tau**(-0.75_wp )              &
                         + n05*0.875_wp*del**7*tau**(-0.125_wp)              &
                         + n06*2.375_wp*del   *tau**( 1.375_wp)*exp(-del)    &
                         + n07*2.0_wp  *del**2*tau             *exp(-del)    &
                         + n08*2.125_wp*del**5*tau**( 1.125_wp)*exp(-del)    &
                         + n09*3.5_wp  *del   *tau**( 2.5_wp)  *exp(-del**2) &
                         + n10*6.5_wp  *del   *tau**( 5.5_wp)  *exp(-del**2) &
                         + n11*4.75_wp *del**4*tau**( 3.75_wp) *exp(-del**2) &
                         + n12*12.5_wp *del**2*tau**( 11.5_wp) *exp(-del**3)

                d2pr_tau2= n01*0.25_wp *(-0.75_wp) *del   *tau**(-1.75_wp)               &
                         + n02*1.25_wp *( 0.25_wp) *del   *tau**(-0.75_wp)               &
                         + n03*1.5_wp  *( 0.5_wp)  *del   *tau**(-0.5_wp)                &
                         + n04*0.25_wp *(-0.75_wp) *del**3*tau**(-1.75_wp)               &
                         + n05*0.875_wp*(-0.125_wp)*del**7*tau**(-1.125_wp)              &
                         + n06*2.375_wp*( 1.375_wp)*del   *tau**( 0.375_wp)*exp(-del)    &
                         + n07*2.0_wp              *del**2                 *exp(-del)    &
                         + n08*2.125_wp*( 1.125_wp)*del**5*tau**( 0.125_wp)*exp(-del)    &
                         + n09*3.5_wp  *( 2.5_wp)  *del   *tau**( 1.5_wp)  *exp(-del**2) &
                         + n10*6.5_wp  *( 5.5_wp)  *del   *tau**( 4.5_wp)  *exp(-del**2) &
                         + n11*4.75_wp *( 3.75_wp) *del**4*tau**( 2.75_wp) *exp(-del**2) &
                         + n12*12.5_wp *( 11.5_wp) *del**2*tau**( 10.5_wp) *exp(-del**3)


                ecalc= 100.0_wp + zci*(d1pr_tau1+d1p0_tau1)
                dedT =          - zci*(d2pr_tau2+d2p0_tau2)/T**2

                ! function
                func= ecalc - ea

                ! Newton's update
                T1= T - func/dedT

                err = abs(T1-T)/T
                if (err.lt.tol) then
                   Tmp(i,j,k)= T1*tc
                   T= T1*tc
                   icont= 1
                   exit tloop
                endif
                T = T1
             enddo tloop

             ! Pressure [inlining instead of calling pcalc_Tro]
             ! ------------------------------------------------
             !! prs(i,j,k) = pcalc_tro(Tmp(i,j,k),ro(i,j,k))

             T= Tmp(i,j,k)*Tci
             tau= 1.0_wp/T

             d1pr_del1= n01              *tau**(0.25_wp )                                     &
                      + n02              *tau**(1.25_wp )                                     &
                      + n03              *tau**(1.5_wp  )                                     &
                      + n04*3.0_wp*del**2*tau**(0.25_wp )                                     &
                      + n05*7.0_wp*del**6*tau**(0.875_wp)                                     &
                      + n06              *tau**(2.375_wp)*exp(-del   )*(1.0_wp-     del   )   &
                      + n07       *del   *tau**(2       )*exp(-del   )*(2.0_wp-     del   )   &
                      + n08       *del**4*tau**(2.125_wp)*exp(-del   )*(5.0_wp-     del   )   &
                      + n09              *tau**(3.5_wp  )*exp(-del**2)*(1.0_wp-2.0_wp*del**2) &
                      + n10              *tau**(6.5_wp  )*exp(-del**2)*(1.0_wp-2.0_wp*del**2) &
                      + n11       *del**3*tau**(4.75_wp )*exp(-del**2)*(4.0_wp-2.0_wp*del**2) &
                      + n12       *del   *tau**(12.5_wp )*exp(-del**3)*(2.0_wp-3.0_wp*del**3)

             prs(i,j,k)= (del*zci/tau*(1.0_wp+del*d1pr_del1))*pc

             ! Sound speed [inlining instead of calling c2calc_Tro]
             ! ----------------------------------------------------
             ! [used to compute spectral radius, stats, characteristics]
             !! c2(i,j,k) = c2calc_tro(Tmp(i,j,k),ro(i,j,k))

             ! Viscosity from Chung-Lee-Starling's law [inlining instead of viscosity_law]
             ! ---------------------------------------------------------------------------
             ! Result is scaled (*diffscale) to respect Reynolds analogy
             visc(i,j,k)= viscosity_law(Tmp(i,j,k),rho(i,j,k))*diffscale

             ! Thermal conductivity from Chung-Lee-Starling's law [inlining instead of thconductivity]
             ! ---------------------------------------------------------------------------------------
             ! Result is scaled (*diffscale) to respect Reynolds analogy
             cok(i,j,k)= thconductivity(visc(i,j,k),Tmp(i,j,k),rho(i,j,k))*diffscale          
          enddo
       enddo
    enddo

    if (icont==0) then
       call mpistop('Error in primitives swp:',1)
    endif

  end subroutine primitives_visc_swp

end submodule smod_primitives_swp
