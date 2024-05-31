!=================================================================================
submodule (mod_primitives) smod_primitives_pfg
!=================================================================================
  !> author: XG
  !> date: January 2024
  !> Computation of primitive variables and viscosity for perfect gas EoS
!=================================================================================

contains

  !================================================================================
  module subroutine primitives_pfg
  !================================================================================
    !> Computation of primitive variables from conservative ones
    !> - version for perfect gas EoS -
    !> [applied on rho_n,rhou_n,... for inviscid fluxes]
  !================================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: rro,v
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
             if (rro==0.) print *,i,j,k

             ! Specific volume
             ! ---------------
             v=1.0_wp/rro

             ! Velocity component
             ! ------------------
             uu(i,j,k)=rhou_n(i,j,k)*v
             vv(i,j,k)=rhov_n(i,j,k)*v
             ww(i,j,k)=rhow_n(i,j,k)*v

             ! Pressure
             ! --------
             prs(i,j,k)= gam1*(rhoe_n(i,j,k)-0.5_wp*rro*(uu(i,j,k)**2+vv(i,j,k)**2+ww(i,j,k)**2))

             ! Temperature !!! ??? only for wall condition ????
             ! -----------
             Tmp(i,j,k)= prs(i,j,k)*v/rg
          enddo
       enddo
    enddo

  end subroutine primitives_pfg

  !================================================================================
  module subroutine primitives_visc_pfg
  !================================================================================
    !> 1/ computation of primitive variables from conservative ones
    !> 2/ computation of thermo-physical properties (visc. & therm. cond.)
    !> 3/ computation of sound speed
    !> - version for perfect gas EoS -
    !> [applied on rho,rhou,... for viscous fluxes]
    !> ATTENTION IL FAUT ENLEVER LES TESTS !!
  !================================================================================
    use mod_constant ! for: diffscale
    use mod_tranprop ! for: visc. & therm. cond. parameters
    use mod_mpi      ! for: iproc (print in debug mode)
    use warnstop     ! for: debug mode
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: rro,v
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

             ! Pressure
             ! --------
             prs(i,j,k)= gam1*(rhoe(i,j,k)-0.5_wp*rro*(uu(i,j,k)**2+vv(i,j,k)**2+ww(i,j,k)**2))

             ! Temperature
             ! -----------
             Tmp(i,j,k)= prs(i,j,k)*v/rg

             ! DEBUG mode: test of negative temperatures
             ! ----------
             if (Tmp(i,j,k)<0.0_wp) then
                print *,Tmp(i,j,k),i+coord(1)*nx,j+coord(2)*ny,k+coord(3)*nz,'proc',iproc,'block',nob(iproc)
                call mpistop('negative temp',1)
             endif

             ! Sound speed [used to compute spectral radius, stats, characteristics]
             ! -----------
             c_(i,j,k)= sqrt(gam*rg*Tmp(i,j,k))

             ! Viscosity from Sutherland's law
             ! -------------------------------
             ! Result is scaled (*diffscale) to respect Reynolds analogy
             visc(i,j,k)= c0*sqrt(Tmp(i,j,k)**3)/(Tmp(i,j,k)+S0)*diffscale

             ! Thermal conductivity from constant Prandtl assumption
             ! -----------------------------------------------------
             cok(i,j,k)= var2*visc(i,j,k)

          enddo
       enddo
    enddo

  end subroutine primitives_visc_pfg

end submodule smod_primitives_pfg
