!=================================================================================
submodule (mod_primitives) smod_primitives_prs
!=================================================================================
  !> author: XG
  !> date: January 2024
  !> Computation of primitive variables and viscosity for Peng-Robinson-Stryjek-Vera EoS
!=================================================================================

contains

  !================================================================================
  module subroutine primitives_prs
  !================================================================================
    !> Computation of primitive variables from conservative ones
    !> - fully inlined version for Peng-Robinson EoS -
    !> [applied on rho_n,rhou_n,... for inviscid fluxes]
  !================================================================================
    use mod_ineos_prs ! for: constants of PRSV model
    use warnstop      ! for: non-convergent Newton
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k,l,icont
    real(wp) :: rro,v,rroe
    real(wp) :: T,T1,ekTr,fn,der_fn,err
    real(wp) :: alp,dalpdT,d2alpdT2,cc,ct1,ct2,iv1,iv2,cv
    ! ----------------------------------------------------------------------------

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
             
             cc=log(abs(2.0_wp*v+2.0_wp*bpr+bpr*sqr8)/abs(2.0_wp*v+2.0_wp*bpr-bpr*sqr8))/(bpr*sqr8)
             ! initial guess
             T=Tmp(i,j,k)
             ! Newton's loop
             tloop: do l=1,100
                ct1=1.0_wp+Kpr*(1.0_wp-sqrt(T/Tc))
                ct2=Kpr/sqrt(T*Tc)

                alp=ct1**2
                dalpdT=-ct1*ct2
                d2alpdT2=0.5_wp*ct2*(ct2+ct1/T)

                ! function & derivative
                fn= cvinf/(nexp+1.0_wp)*(abs(T/Tc))**(nexp+1.0_wp)*Tc &
                     - apr*(alp-T*dalpdT)*cc - rroe/rro

                der_fn= cvinf*abs(T/Tc)**nexp + apr*d2alpdT2*T*cc

                ! update solution
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
             ct1=1.0_wp+Kpr*(1.0_wp-sqrt(T/Tc))
             ct2=Kpr/sqrt(T*Tc)

             alp=ct1**2
             dalpdT=-ct1*ct2
             d2alpdT2=0.5_wp*ct2*(ct2+ct1/T)

             iv1=1.0_wp/(v-bpr)
             iv2=1.0_wp/(v**2+2.0_wp*v*bpr-bpr**2)

             prs(i,j,k)=rg*T*iv1 - apr*alp*iv2
          enddo
       enddo
    enddo

    if (icont==0) then
       call mpistop('Error in primitives prs',1)
    endif

  end subroutine primitives_prs

  !================================================================================
  module subroutine primitives_visc_prs
  !================================================================================
    !> 1/ computation of primitive variables from conservative ones
    !> 2/ computation of thermo-physical properties (visc. & therm. cond.)
    !> 3/ computation of sound speed
    !> - fully inlined version for Peng-Robinson EoS -
    !> [applied on rho,rhou,... for viscous fluxes]
    !> ATTENTION IL FAUT ENLEVER LES TESTS !!
  !================================================================================
    use mod_constant  ! for: diffscale
    use mod_tranprop  ! for: visc. & therm. cond. parameters
    use mod_ineos_prs ! for: constants of PRSV model
    use mod_mpi       ! for: iproc (print in debug mode)
    use warnstop      ! for: debug mode & non-convergent Newton
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k,l,icont
    real(wp) :: rro,v,rroe
    real(wp) :: T,T1,ekTr,fn,der_fn,err
    real(wp) :: alp,dalpdT,d2alpdT2,cc,ct1,ct2,iv1,iv2,cv
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
             cc=log(abs(2.0_wp*v+2.0_wp*bpr+bpr*sqr8)/abs(2.0_wp*v+2.0_wp*bpr-bpr*sqr8))/(bpr*sqr8)

             ! initial guess
             T=Tmp(i,j,k)

             ! Newton's loop
             tloop: do l=1,100

                ! DEBUG mode: test of negative temperatures
                ! ----------              
                if (T<0.0_wp) print *,T,i+coord(1)*nx,j+coord(2)*ny,k+coord(3)*nz,'proc',iproc,'bl',nob(iproc),'it',l

                ct1=1.0_wp+Kpr*(1.0_wp-sqrt(T/Tc))
                ct2=Kpr/sqrt(T*Tc)

                alp=ct1**2
                dalpdT=-ct1*ct2
                d2alpdT2=0.5_wp*ct2*(ct2+ct1/T)

                ! function & derivative
                fn= cvinf/(nexp+1.0_wp)*(abs(T/Tc))**(nexp+1.0_wp)*Tc &
                     - apr*(alp-T*dalpdT)*cc - rroe/rro

                der_fn= cvinf*abs(T/Tc)**nexp + apr*d2alpdT2*T*cc

                ! update solution
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
             ct1=1.0_wp+Kpr*(1.0_wp-sqrt(T/Tc))
             ct2=Kpr/sqrt(T*Tc)

             alp=ct1**2
             dalpdT=-ct1*ct2
             d2alpdT2=0.5_wp*ct2*(ct2+ct1/T)

             iv1=1.0_wp/(v-bpr)
             iv2=1.0_wp/(v**2+2.0_wp*v*bpr-bpr**2)

             prs(i,j,k)=rg*T*iv1 - apr*alp*iv2

             ! Sound speed [inlining instead of calling c2calc_Tro]
             ! ----------------------------------------------------
             ! [used to compute spectral radius, stats, characteristics]
             !! c_(i,j,k) = c2calc_tro(Tmp(i,j,k),ro(i,j,k))

             cv= cvinf*abs(T/Tc)**nexp

             c_(i,j,k)= sqrt((-v**2)*(-rg*T*iv1**2 &
                  + 2.0_wp*apr*alp*(v+bpr)*iv2**2 &
                  - T/cv*(rg*iv1-dalpdT*apr*iv2)**2))

             cv= cv+ apr*d2alpdT2*T*cc

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

    if (icont==0) then
       call mpistop('Error in primitives prs:',1)
    endif

  end subroutine primitives_visc_prs

end submodule smod_primitives_prs
