!===============================================================================
subroutine alloc_tab
!===============================================================================
  !> Dynamic allocation of flow arrays (and stats)
  !> NOT FINALIZED (to be cleaned)
!===============================================================================
  use mod_flow
  use mod_flow0
  use mod_constant
  use mod_mpi
  use mod_time
  use mod_io
  use mod_rans
  use mod_block    ! for nbloc
  ! use mod_turb_model_length_scale
  implicit none
  ! ----------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------

  if (iproc==0) write(*,*) 'Matrix allocation..'

  ! ---------------------------------------------------------------------------
  if (is_curv3) then ! 3D curvilinear geometry

     if (TURB) then
        ! 2D stats (plane) averaged along z
        nstat=177
        allocate(avg_s(nx,ny,nstat))
        allocate(avg_t(nx,ny,nstat))

        ! 3D stats (volume) + Initialize
        nstat_v=18
        allocate(avg_v(nx,ny,nz,nstat_v))
        avg_v=0.0_wp

        ! 2D stats on walls (plane) + Initialize
        ! /!\ ONLY implemented for walls at jmin
        if (is_bc_wall(2,1)) then
           nstat_w=18
           allocate(avg_w(nx,nz,nstat_w))
           avg_w=0.0_wp
        endif
     endif

     if (CYL) then ! in fact for SPHERE
        ! -> statistics in volume
        nstat=18
        allocate(avg_v(nx,ny,nz,nstat))
        ! Initialize
        avg_v=0.0_wp
     endif

!!$     if (CYL.or.TURB) then
!!$        ! 3D curvilinear geometry -> statistics in volume
!!$        nstat=18
!!$        allocate(avg_v(nx,ny,nz,nstat))
!!$        ! Initialize
!!$        avg_v=0.0_wp
!!$     endif
  else
     ! Stats
     if (CHIT.or.TGV) then
        nstat=200
        allocate( avg_s(1,1,nstat))
        allocate( avg_t(1,1,nstat))
        allocate(avg_tg(1,1,nstat))
     elseif (CHAN) then
        if (is_RANS) then
           nstat = 38
        else
           nstat = 180
        endif
        allocate(avg_s(1,ny,nstat))
        allocate(avg_t(1,ny,nstat))
        allocate(avg_tg(1,ngy,nstat))
     elseif (STBL.or.SRC.or.LE.or.T3C) then
        nstat=167
        allocate(avg_s(nx,ny,nstat))
        allocate(avg_t(nx,ny,nstat))
     elseif (TURB.or.TE) then
        nstat=177
        allocate( avg_s(nx,ny,nstat))
        allocate( avg_t(nx,ny,nstat))
     elseif (CYL) then
        nstat=167
        allocate(avg_s(nx,ny,nstat))
        allocate(avg_t(nx,ny,nstat))
     elseif (SHIT) then
        nstat=30
        allocate(avg_s(nx,ny,nstat))
        allocate(avg_t(nx,ny,nstat))
     elseif (ACT) then
        nstat=98
        allocate(avg_s(nx,ny,nstat))
        allocate(avg_t(nx,ny,nstat))
     endif

     if (CHIT.or.TGV.or.CHAN.or.STBL.or.SRC.or.CYL.or.SHIT.or.ACT.or.TURB.or.LE) then
        avg_s=0.0_wp
        avg_t=0.0_wp
     endif
  endif

  if (idepart.eq.POST_PROCESSING) then
     if (.not.is_curv3) deallocate(avg_s)
!!$  allocate( rho(nx1:nx2,ny1:ny2,nz1:nz2))
!!$  allocate( Tmp(nx1:nx2,ny1:ny2,nz1:nz2))
!!$  allocate( prs(nx1:nx2,ny1:ny2,nz1:nz2))
!!$  allocate(  uu(nx1:nx2,ny1:ny2,nz1:nz2))
!!$  allocate(  vv(nx1:nx2,ny1:ny2,nz1:nz2))
!!$  allocate(  ww(nx1:nx2,ny1:ny2,nz1:nz2))
     return
  endif   

  ! ---------------------------------------------------------------------------
  allocate( rho(nx1:nx2,ny1:ny2,nz1:nz2))
  allocate(rhou(nx1:nx2,ny1:ny2,nz1:nz2))
  allocate(rhov(nx1:nx2,ny1:ny2,nz1:nz2))
  allocate(rhow(nx1:nx2,ny1:ny2,nz1:nz2))
  allocate(rhoe(nx1:nx2,ny1:ny2,nz1:nz2))
  ! ---------------------------------------------------------------------------
  allocate( Tmp(nx1:nx2,ny1:ny2,nz1:nz2))
  allocate( prs(nx1:nx2,ny1:ny2,nz1:nz2))
  allocate( cok(nx1:nx2,ny1:ny2,nz1:nz2))
  allocate(visc(nx1:nx2,ny1:ny2,nz1:nz2))
  allocate(  c_(nx1:nx2,ny1:ny2,nz1:nz2))
  allocate(  uu(nx1:nx2,ny1:ny2,nz1:nz2))
  allocate(  vv(nx1:nx2,ny1:ny2,nz1:nz2))
  allocate(  ww(nx1:nx2,ny1:ny2,nz1:nz2))
  ! Velocity and Temperatyre gradients
  allocate(dux(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
  allocate(dvx(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
  allocate(dwx(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
  allocate(duy(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
  allocate(dvy(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
  allocate(dwy(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
  allocate(duz(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
  allocate(dvz(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
  allocate(dwz(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
  allocate(dTx(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
  allocate(dTy(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
  allocate(dTz(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
  ! initialize derivatives for corners
  ! ----------------------------------
  ! [used in calculation of Ducros' sensor in mod_artvisc_shock]
  dux=0.0_wp
  duy=0.0_wp
  duz=0.0_wp
  dvx=0.0_wp
  dvy=0.0_wp
  dvz=0.0_wp
  dwx=0.0_wp
  dwy=0.0_wp
  dwz=0.0_wp

  ! local CFL number
  ! ----------------
  allocate(cfl_i(nx,ny,nz),cfl_j(nx,ny,nz))

  ! ---------------------------------------------------------------------------
  ! if (coord(1)==0) allocate(uturb(ny,nz1:nz2),vturb(ny,nz1:nz2),wturb(ny,nz1:nz2))
  ! ---------------------------------------------------------------------------

  if (idepart.eq.POST_PROCESSING) return
  ! ---------------------------------------------------------------------------
  allocate( rho_n(nx1:nx2,ny1:ny2,nz1:nz2))
  allocate(rhoe_n(nx1:nx2,ny1:ny2,nz1:nz2))
  allocate(rhou_n(nx1:nx2,ny1:ny2,nz1:nz2))
  allocate(rhov_n(nx1:nx2,ny1:ny2,nz1:nz2))
  allocate(rhow_n(nx1:nx2,ny1:ny2,nz1:nz2))

!!$  allocate( rhon(nx1:nx2,ny1:ny2,nz1:nz2))
!!$  allocate(rhoen(nx1:nx2,ny1:ny2,nz1:nz2))
!!$  allocate(rhoun(nx1:nx2,ny1:ny2,nz1:nz2))
!!$  allocate(rhovn(nx1:nx2,ny1:ny2,nz1:nz2))
!!$  allocate(rhown(nx1:nx2,ny1:ny2,nz1:nz2))

  ! Increments arrays
  ! =================
  ! Nota: allocated in mod_init_irs for IRS (use of extended ghost points)
  if (.not.is_irs) then
     allocate( Krho(nx1:nx2,ny1:ny2,nz1:nz2))
     allocate(Krhou(nx1:nx2,ny1:ny2,nz1:nz2))
     allocate(Krhov(nx1:nx2,ny1:ny2,nz1:nz2))
     allocate(Krhow(nx1:nx2,ny1:ny2,nz1:nz2))
     allocate(Krhoe(nx1:nx2,ny1:ny2,nz1:nz2))
     ! added by Camille ?? why ??
     Krho=0.0_wp
     Krhou=0.0_wp
     Krhov=0.0_wp
     Krhow=0.0_wp
     Krhoe=0.0_wp
  endif
  
  allocate( Frho(nx1:nx2,ny1:ny2,nz1:nz2))
  allocate(Frhou(nx1:nx2,ny1:ny2,nz1:nz2))
  allocate(Frhov(nx1:nx2,ny1:ny2,nz1:nz2))
  allocate(Frhow(nx1:nx2,ny1:ny2,nz1:nz2))
  allocate(Frhoe(nx1:nx2,ny1:ny2,nz1:nz2))
  allocate( Grho(nx1:nx2,ny1:ny2,nz1:nz2))
  allocate(Grhou(nx1:nx2,ny1:ny2,nz1:nz2))
  allocate(Grhov(nx1:nx2,ny1:ny2,nz1:nz2))
  allocate(Grhow(nx1:nx2,ny1:ny2,nz1:nz2))
  allocate(Grhoe(nx1:nx2,ny1:ny2,nz1:nz2))
  allocate( Hrho(nx1:nx2,ny1:ny2,nz1:nz2))
  allocate(Hrhou(nx1:nx2,ny1:ny2,nz1:nz2))
  allocate(Hrhov(nx1:nx2,ny1:ny2,nz1:nz2))
  allocate(Hrhow(nx1:nx2,ny1:ny2,nz1:nz2))
  allocate(Hrhoe(nx1:nx2,ny1:ny2,nz1:nz2))

  ! Modif pour RANS
  if (is_RANS) then
     if (model_RANS.eq.'SA') then
        allocate(  nutil(nx1:nx2,ny1:ny2,nz1:nz2))
        allocate(    mut(nx1:nx2,ny1:ny2,nz1:nz2))
        allocate(    nut(nx1:nx2,ny1:ny2,nz1:nz2))
        allocate(nutil_n(nx1:nx2,ny1:ny2,nz1:nz2))
        allocate(dnutilx(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
        allocate(dnutily(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
        allocate(dnutilz(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
        allocate( Fnutil(nx1:nx2,ny1:ny2,nz1:nz2))
        allocate( Gnutil(nx1:nx2,ny1:ny2,nz1:nz2))
        allocate( Hnutil(nx1:nx2,ny1:ny2,nz1:nz2))
        allocate( Knutil(nx1:nx2,ny1:ny2,nz1:nz2))
        allocate(  Sterm(nx1:nx2,ny1:ny2,nz1:nz2))
        allocate(lengthscale(nx,ny,nz))
     endif
  endif


  ! Modif pour Smago
  if (is_sgs_model) then
     allocate(S11(nx1:nx2,ny1:ny2,nz1:nz2))
     allocate(S12(nx1:nx2,ny1:ny2,nz1:nz2))
     allocate(S13(nx1:nx2,ny1:ny2,nz1:nz2))
     allocate(S23(nx1:nx2,ny1:ny2,nz1:nz2))
     allocate(S22(nx1:nx2,ny1:ny2,nz1:nz2))
     allocate(S33(nx1:nx2,ny1:ny2,nz1:nz2))
     allocate(S11f(nx1:nx2,ny1:ny2,nz1:nz2))
     allocate(S12f(nx1:nx2,ny1:ny2,nz1:nz2))
     allocate(S13f(nx1:nx2,ny1:ny2,nz1:nz2))
     allocate(S23f(nx1:nx2,ny1:ny2,nz1:nz2))
     allocate(S22f(nx1:nx2,ny1:ny2,nz1:nz2))
     allocate(S33f(nx1:nx2,ny1:ny2,nz1:nz2))
  end if

  ! for Tam & Dong's BC
  ! ===================

  ! old version [kept temporarily, eg for filter_fluct,...]
  if (is_mean0) then
     allocate(rho0(nx1:nx2,ny1:ny2,nz1:nz2))
     allocate(  p0(nx1:nx2,ny1:ny2,nz1:nz2))
     allocate(  T0(nx1:nx2,ny1:ny2,nz1:nz2))
     allocate(  u0(nx1:nx2,ny1:ny2,nz1:nz2))
     allocate(  v0(nx1:nx2,ny1:ny2,nz1:nz2))
     allocate(  w0(nx1:nx2,ny1:ny2,nz1:nz2))
  end if

  ! new version
  if (is_TamDong) then
     if (is_TamDong_1pt) then
        if (is_bc_TD(1,1)) allocate(BC_face(1,1)%U0( 1:5,ny1:ny2,nz1:nz2,6))
        if (is_bc_TD(1,2)) allocate(BC_face(1,2)%U0(-3:1,ny1:ny2,nz1:nz2,6))
        if (is_bc_TD(2,1)) allocate(BC_face(2,1)%U0(nx1:nx2, 1:5,nz1:nz2,6))
        if (is_bc_TD(2,2)) allocate(BC_face(2,2)%U0(nx1:nx2,-3:1,nz1:nz2,6))
        if (is_TamDong3D) then
           if (is_bc_TD(3,1)) allocate(BC_face(3,1)%U0(nx1:nx2,ny1:ny2, 1:5,6))
           if (is_bc_TD(3,2)) allocate(BC_face(3,2)%U0(nx1:nx2,ny1:ny2,-3:1,6))
        endif
     else
        if (is_bc_TD(1,1)) allocate(BC_face(1,1)%U0(1:ngh+3,ny1:ny2,nz1:nz2,6))
        if (is_bc_TD(1,2)) allocate(BC_face(1,2)%U0( -2:ngh,ny1:ny2,nz1:nz2,6))
        if (is_bc_TD(2,1)) allocate(BC_face(2,1)%U0(nx1:nx2,1:ngh+3,nz1:nz2,6))
        if (is_bc_TD(2,2)) allocate(BC_face(2,2)%U0(nx1:nx2, -2:ngh,nz1:nz2,6))
        if (is_TamDong3D) then
           if (is_bc_TD(3,1)) allocate(BC_face(3,1)%U0(nx1:nx2,ny1:ny2,1:ngh+3,6))
           if (is_bc_TD(3,2)) allocate(BC_face(3,2)%U0(nx1:nx2,ny1:ny2, -2:ngh,6))
        endif
     endif
  endif

    ! for perturbed Riemann BC
  ! ========================

  if (BC_face(1,1)%sort==-41) allocate(BC_face(1,1)%U0R(1,1:ny,1:nz,3))
  if (BC_face(1,2)%sort==-41) allocate(BC_face(1,2)%U0R(1,1:ny,1:nz,3))
  if (BC_face(2,1)%sort==-41) allocate(BC_face(2,1)%U0R(1:nx,1,1:nz,3))
  if (BC_face(2,2)%sort==-41) allocate(BC_face(2,2)%U0R(1:nx,1,1:nz,3))
  if (.not.is_2d) then
     if (BC_face(3,1)%sort==-41) allocate(BC_face(3,1)%U0R(1:nx,1:ny,1,3))
     if (BC_face(3,2)%sort==-41) allocate(BC_face(3,2)%U0R(1:nx,1:ny,1,3))
  endif

end subroutine alloc_tab

!===============================================================================
subroutine lib_tab
!===============================================================================
  !> Free dynamic allocation of flow arrays (and stats) and MPI-types
  !> NOT UP-TO-DATE
!===============================================================================
  use mod_flow
  use mod_constant
  use mod_mpi_part
  use mod_mpi_types_two_sided
  use mod_rans
  use mod_turb_model_length_scale
  implicit none
  ! ----------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------

!!$  deallocate(x,y,z,dx,dy,dz)
!!$  
!!$  call mpistop('', 0)
!!$  
!!$  deallocate(x,y,z,dx,dy,dz,dx2,dy2,dz2)
!!$
!!$  call mpistop('', 0)
  
  deallocate(rho,rhou,rhov,rhow,rhoe)
  deallocate(rho_n,rhou_n,rhov_n,rhow_n,rhoe_n)
  if (is_sgs_model) then
     deallocate(S11,S12,S13,S22,S23,S33)
     deallocate(S11f,S12f,S13f,S22f,S23f,S33f)
  end if
  ! deallocate(prs,uu,vv,ww,Tmp)
  deallocate(prs,Tmp)

  deallocate(Krho,Krhou,Krhov,Krhow,Krhoe)
  deallocate(Frho,Frhou,Frhov,Frhow,Frhoe)
  deallocate(Grho,Grhou,Grhov,Grhow,Grhoe)
  deallocate(Hrho,Hrhou,Hrhov,Hrhow,Hrhoe)

  if (is_curv3) then
     deallocate(xc,yc,zc3)
  else
     if (is_curv) then
        ! deallocate(xgc,ygc,xc,yc) ! to be suppressed
        deallocate(xc,yc)
     else
        ! deallocate(xg,yg,x,y) ! to be suppressed
        deallocate(x,y)
     endif
     deallocate(z)
     ! deallocate(z,zg) ! to be suppressed

     ! Modif pour RANS
     if (is_RANS) then
        if (model_RANS.eq.'SA') then
           deallocate(nutil,mut,nut,nutil_n)
           deallocate(dnutilx,dnutily,dnutilz)
           deallocate(Fnutil,Gnutil,Hnutil,Knutil)
           deallocate(lengthscale)
        endif
     endif
  endif



  ! Free MPI types
  ! ==============
  !call free_mpi_types_comm
  !call free_mpi_types_IO

end subroutine lib_tab
