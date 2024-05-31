!=================================================================================
submodule (mod_primitives) smod_primitives_vdw
!=================================================================================
  !> author: XG
  !> date: January 2024
  !> Computation of primitive variables and viscosity for van der Waals EoS
!=================================================================================

contains

  !================================================================================
  module subroutine primitives_vdw
  !================================================================================
    !> Computation of primitive variables from conservative ones
    !> - fully inlined version for van der Waals EoS -
    !> [applied on rho_n,rhou_n,... for inviscid fluxes]
  !================================================================================
    use mod_ineos_vdw ! for: constants of VdW model
    implicit none
    ! ----------------------------------------------------------------------------
    integer  :: i,j,k
    real(wp) :: rro,v,rroe
    ! ----------------------------------------------------------------------------

    ! Compute primitive variables from conservative ones
    ! ==================================================

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
             Tmp(i,j,k)= gam1/rg*(rroe*v+rro*avw)

             ! Pressure [inlining instead of calling pcalc_Tro]
             ! ------------------------------------------------
             !! prs(i,j,k) = pcalc_tro(Tmp(i,j,k),rho_n(i,j,k))
             prs(i,j,k)= rro*rg*Tmp(i,j,k)/(1.0_wp-rro*bvw)-rro**2*avw

          enddo
       enddo
    enddo

  end subroutine primitives_vdw

  !================================================================================
  module subroutine primitives_visc_vdw
  !================================================================================
    !> 1/ computation of primitive variables from conservative ones
    !> 2/ computation of thermo-physical properties (visc. & therm. cond.)
    !> 3/ computation of sound speed
    !> - fully inlined version for van der Waals EoS -
    !> [applied on rho,rhou,... for viscous fluxes]
    !> ATTENTION IL FAUT ENLEVER LES TESTS !!
  !================================================================================
    use mod_constant  ! for: diffscale
    use mod_tranprop  ! for: visc. & therm. cond. parameters
    use mod_ineos_vdw ! for: constants of VdW model
    implicit none
    ! ----------------------------------------------------------------------------
    integer  :: i,j,k
    real(wp) :: rro,v,rroe
    ! ----------------------------------------------------------------------------

    ! Compute primitive variables, viscosity, thermal cond. and sound speed
    ! =====================================================================

    ! for all points + ghost cells
    do k=ndzt,nfzt
       do j=ndyt,nfyt
          !!!dir$ simd
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
             Tmp(i,j,k)= gam1/rg*(rroe/rro+rro*avw)

             ! Pressure [inlining instead of calling pcalc_Tro]
             ! ------------------------------------------------
             !! prs(i,j,k) = pcalc_tro(Tmp(i,j,k),ro(i,j,k))
             prs(i,j,k)= rro*rg*Tmp(i,j,k)/(1.0_wp-rro*bvw)-rro**2*avw

             ! Sound speed [inlining instead of calling c2calc_Tro]
             ! ----------------------------------------------------
             ! [used to compute spectral radius, stats, characteristics]
             !! c2(i,j,k) = c2calc_tro(Tmp(i,j,k),ro(i,j,k))
             c_(i,j,k)= sqrt(gam*rg*Tmp(i,j,k)/(1.0_wp-rro*bvw)**2-2.0_wp*avw*rro)

             ! Viscosity from power law [inlining instead of calling viscosity_law]
             ! ---------------------------------------------------------------------------
             !! visc(i,j,k) = viscosity_law(Tmp(i,j,k),ro(i,j,k))

             ! Result is scaled (*diffscale) to respect Reynolds analogy
             visc(i,j,k) = var3*Tmp(i,j,k)**nexp*diffscale

             ! Thermal conductivity from Chung-Lee-Starling's law [inlining instead of thconductivity]
             ! ---------------------------------------------------------------------------------------
             !! cok(i,j,k) = thconductivity(visc(i,j,k),Tmp(i,j,k),ro(i,j,k))

             ! Result is scaled (*diffscale) to respect Reynolds analogy
             cok(i,j,k) = var2*visc(i,j,k)*diffscale
          enddo
       enddo
    enddo

  end subroutine primitives_visc_vdw

end submodule smod_primitives_vdw
