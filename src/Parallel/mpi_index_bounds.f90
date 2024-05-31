!===============================================================
subroutine mpi_index_bounds
!===============================================================
  !> Index bounds for arrays in MPI
  !
  !  ngh : number of ghost cells for Eulerian fluxes
  ! ngh_v: number of ghost cells for viscous fluxes
!===============================================================
  use mod_grid
  use mod_bc
  use mod_constant ! <- for is_wall2
  use mod_mpi
  implicit none
  ! ------------------------------------------------------------
  ! ------------------------------------------------------------

  ! Index bounds for inviscid fluxes
  ! ================================

  ! Definitions of nd_e and nf_e
  ! ----------------------------
  ! This index corresponds to the points where Eulerian terms
  ! are advanced [in flux_euler_XX], by default between 1 and nx
  ! (for interior points).
  ! -> not advanced in BCs
  !    (1 cell layer for BCs coded on 1 point [is_bc_1pt] and
  !     5 cell layers for BCs coded on 5 points)
  ndx_e=1
  if (BC_face(1,1)%sort==0) then
     if (.not.is_wall2) ndx_e=ndx_e+1
  elseif ((BC_face(1,1)%sort<=-3)) then
     ndx_e=ndx_e+1
  elseif ((BC_face(1,1)%sort==-1).or.(BC_face(1,1)%sort==-2)) then
     ndx_e=ndx_e+ngh
  endif
  nfx_e=nx
  if (BC_face(1,2)%sort==0) then
     if (.not.is_wall2) nfx_e=nfx_e-1
  elseif ((BC_face(1,2)%sort<=-3)) then
     nfx_e=nfx_e-1
  elseif ((BC_face(1,2)%sort==-1).or.(BC_face(1,2)%sort==-2)) then
     nfx_e=nfx_e-ngh
  endif
  ndy_e=1
  if (BC_face(2,1)%sort==0) then
     if (.not.is_wall2) ndy_e=ndy_e+1
  elseif ((BC_face(2,1)%sort<=-3)) then
     ndy_e=ndy_e+1
  elseif ((BC_face(2,1)%sort==-1).or.(BC_face(2,1)%sort==-2)) then
     ndy_e=ndy_e+ngh
  endif
  nfy_e=ny
  if (BC_face(2,2)%sort==0) then
     if (.not.is_wall2) nfy_e=nfy_e-1
  elseif ((BC_face(2,2)%sort<=-3)) then
     nfy_e=nfy_e-1
  elseif ((BC_face(2,2)%sort==-1).or.(BC_face(2,2)%sort==-2)) then
     nfy_e=nfy_e-ngh
  endif
  ndz_e=1
  if (BC_face(3,1)%sort==0) then
     if (.not.is_wall2) ndz_e=ndz_e+1
  elseif ((BC_face(3,1)%sort<=-3)) then
     ndz_e=ndz_e+1
  elseif ((BC_face(3,1)%sort==-1).or.(BC_face(3,1)%sort==-2)) then
     ndz_e=ndz_e+ngh
  endif
  nfz_e=nz
  if (BC_face(3,2)%sort==0) then
     if (.not.is_wall2) nfz_e=nfz_e-1
  elseif ((BC_face(3,2)%sort<=-3)) then
     nfz_e=nfz_e-1
  elseif ((BC_face(3,2)%sort==-1).or.(BC_face(3,2)%sort==-2)) then
     nfz_e=nfz_e-ngh
  endif

  ! Definitions of nd and nf
  ! ------------------------
  ! This index corresponds to the points where Eulerian terms
  ! are advanced [in flux_euler_XX] for interior points
  ! for which boundary points are already advanced.
  ndx=1
  nfx=nx
  if (is_boundary(1,1)) ndx= ndx+ngh
  if (is_boundary(1,2)) nfx= nfx-ngh
  ndy=1
  nfy=ny
  if (is_boundary(2,1)) ndy= ndy+ngh
  if (is_boundary(2,2)) nfy= nfy-ngh
  ndz=1
  nfz=nz
  if (is_boundary(3,1)) ndz= ndz+ngh
  if (is_boundary(3,2)) nfz= nfz-ngh

  ! Definitions of ndt and nft
  ! --------------------------
  ! This index is the index of interior points
  ! augmented by the number of ghost cells (-/+ ngh).
  ndxt=ndx-ngh
  nfxt=nfx+ngh
  ndyt=ndy-ngh
  nfyt=nfy+ngh
  ndzt=ndz-ngh
  nfzt=nfz+ngh

  ! Index bounds for viscous fluxes
  ! ===============================
  ! Definitions of nd_v and nf_v
  ! ----------------------------
  ! This index corresponds to the points where viscous terms
  ! are advanced [in flux_visc_XX], by default between 1 and nx
  ! (for interior points and wall BCs).
  ! -> not advanced in non-reflecting BCs
  !    (1 cell layer for BCs coded on 1 point [is_bc_1pt] and
  !     ngh cell layers for BCs coded on ngh points, ie T&D)
  ndx_v=1
  if ((is_bc_1pt(1,1)).and.(.not.is_bc_wall(1,1))) then
     ndx_v= ndx_v+1
  elseif ((BC_face(1,1)%sort==-1).or.(BC_face(1,1)%sort==-2)) then
     ndx_v= ndx_v+ngh
  endif
  nfx_v=nx
  if ((is_bc_1pt(1,2)).and.(.not.is_bc_wall(1,2))) then
     nfx_v= nfx_v-1
  elseif ((BC_face(1,2)%sort==-1).or.(BC_face(1,2)%sort==-2)) then
     nfx_v= nfx_v-ngh
  endif
  ndy_v=1
  if ((is_bc_1pt(2,1)).and.(.not.is_bc_wall(2,1))) then
     ndy_v= ndy_v+1
  elseif ((BC_face(2,1)%sort==-1).or.(BC_face(2,1)%sort==-2)) then
     ndy_v= ndy_v+ngh
  endif
  nfy_v=ny
  if ((is_bc_1pt(2,2)).and.(.not.is_bc_wall(2,2))) then
     nfy_v= nfy_v-1
  elseif ((BC_face(2,2)%sort==-1).or.(BC_face(2,2)%sort==-2)) then
     nfy_v= nfy_v-ngh
  endif
  ndz_v=1
  if ((is_bc_1pt(3,1)).and.(.not.is_bc_wall(3,1))) then
     ndz_v= ndz_v+1
  elseif ((BC_face(3,1)%sort==-1).or.(BC_face(3,1)%sort==-2)) then
     ndz_v= ndz_v+ngh
  endif
  nfz_v=nz
  if ((is_bc_1pt(3,2)).and.(.not.is_bc_wall(3,2))) then
     nfz_v= nfz_v-1
  elseif ((BC_face(3,2)%sort==-1).or.(BC_face(3,2)%sort==-2)) then
     nfz_v= nfz_v-ngh
  endif

  ! Definitions of nd_vi and nf_vi
  ! ------------------------------
  ! This index corresponds to the points where viscous terms
  ! are advanced [in flux_visc_XX] for interior points
  ! for which boundary points are already advanced.
  ! -> for bc_wall, it is nd_v+ngh_v or nf_v-ngh_v
  ! -> for bc_1pt , it is nd_v+ngh_v-1 or nf_v-ngh_v+1
  ! -> for bc on ngh pts, nd_v or nf_v unchanged
  ndx_vi=ndx_v
  nfx_vi=nfx_v
  if (is_bc_wall(1,1)) ndx_vi= ndx_vi+ngh_v
  if (is_bc_wall(1,2)) nfx_vi= nfx_vi-ngh_v
  if ((is_bc_1pt(1,1)).and.(.not.is_bc_wall(1,1))) ndx_vi= ndx_vi+ngh_v-1
  if ((is_bc_1pt(1,2)).and.(.not.is_bc_wall(1,2))) nfx_vi= nfx_vi-ngh_v+1
  ndy_vi=ndy_v
  nfy_vi=nfy_v
  if (is_bc_wall(2,1)) ndy_vi= ndy_vi+ngh_v
  if (is_bc_wall(2,2)) nfy_vi= nfy_vi-ngh_v
  if ((is_bc_1pt(2,1)).and.(.not.is_bc_wall(2,1))) ndy_vi= ndy_vi+ngh_v-1
  if ((is_bc_1pt(2,2)).and.(.not.is_bc_wall(2,2))) nfy_vi= nfy_vi-ngh_v+1
  ndz_vi=ndz_v
  nfz_vi=nfz_v
  if (is_bc_wall(3,1)) ndz_vi= ndz_vi+ngh_v
  if (is_bc_wall(3,2)) nfz_vi= nfz_vi-ngh_v
  if ((is_bc_1pt(3,1)).and.(.not.is_bc_wall(3,1))) ndz_vi= ndz_vi+ngh_v-1
  if ((is_bc_1pt(3,2)).and.(.not.is_bc_wall(3,2))) nfz_vi= nfz_vi-ngh_v+1

  ! Definitions of ndt_v and nft_v
  ! -------------------------------
  ! This index is the index of interior points
  ! augmented by the number of ghost cells (-/+ ngh_v).
  ndxt_v=ndx_vi-ngh_v
  nfxt_v=nfx_vi+ngh_v
  ndyt_v=ndy_vi-ngh_v
  nfyt_v=nfy_vi+ngh_v
  ndzt_v=ndz_vi-ngh_v
  nfzt_v=nfz_vi+ngh_v

  ! Definitions of ndt_v1 and nft_v1 & ndt_v2 and nft_v2
  ! ----------------------------------------------------
  ! [indices for double derivatives in grad_velT]

  ndx_v1= 1-ngh_v
  nfx_v1=nx+ngh_v
  if (is_boundary(1,1)) ndx_v1=ndx_v1+ngh_v
  if (is_boundary(1,2)) nfx_v1=nfx_v1-ngh_v
  ndy_v1= 1-ngh_v
  nfy_v1=ny+ngh_v
  if (is_boundary(2,1)) ndy_v1=ndy_v1+ngh_v
  if (is_boundary(2,2)) nfy_v1=nfy_v1-ngh_v
  ndz_v1= 1-ngh_v
  nfz_v1=nz+ngh_v
  if (is_boundary(3,1)) ndz_v1=ndz_v1+ngh_v
  if (is_boundary(3,2)) nfz_v1=nfz_v1-ngh_v

  ndx_v2= 1-ngh_v
  nfx_v2=nx+ngh_v
  if (is_boundary(1,1)) ndx_v2=ndx_v2+2*ngh_v
  if (is_boundary(1,2)) nfx_v2=nfx_v2-2*ngh_v
  ndy_v2= 1-ngh_v
  nfy_v2=ny+ngh_v
  if (is_boundary(2,1)) ndy_v2=ndy_v2+2*ngh_v
  if (is_boundary(2,2)) nfy_v2=nfy_v2-2*ngh_v
  ndz_v2= 1-ngh_v
  nfz_v2=nz+ngh_v
  if (is_boundary(3,1)) ndz_v2=ndz_v2+2*ngh_v
  if (is_boundary(3,2)) nfz_v2=nfz_v2-2*ngh_v

  ndx_v3= 1
  nfx_v3=nx
  if (is_boundary(1,1)) ndx_v3=ndx_v3+ngh_v
  if (is_boundary(1,2)) nfx_v3=nfx_v3-ngh_v
  ndy_v3= 1
  nfy_v3=ny
  if (is_boundary(2,1)) ndy_v3=ndy_v3+ngh_v
  if (is_boundary(2,2)) nfy_v3=nfy_v3-ngh_v
  ndz_v3= 1
  nfz_v3=nz
  if (is_boundary(3,1)) ndz_v3=ndz_v3+ngh_v
  if (is_boundary(3,2)) nfz_v3=nfz_v3-ngh_v
  
  ! Index bounds for numerical dissipation
  ! ======================================
  
  ! Definitions of nd_d and nf_d
  ! ----------------------------
  ! This index corresponds to the points where numerical dissipation
  ! is applied [in mod_artvisc], by default between 1 and nx.
  ! -> exclude first and last points at boundaries
  ! -> possibly exclude 3 and 5 pts stencils (here coeff are set to zero)
  ndx_d= 1
  nfx_d=nx 
  if (is_boundary(1,1)) ndx_d=ndx_d+1
  if (is_boundary(1,2)) nfx_d=nfx_d-1
  ndy_d= 1
  nfy_d=ny
  if (is_boundary(2,1)) ndy_d=ndy_d+1
  if (is_boundary(2,2)) nfy_d=nfy_d-1
  ndz_d= 1
  nfz_d=nz
  if (is_boundary(3,1)) ndz_d=ndz_d+1
  if (is_boundary(3,2)) nfz_d=nfz_d-1

  ! Definitions of nd_d1 and nf_d1
  ! ----------------------------
  ! extended index for sensor calculation (shock-capturing)
  ndx_d1=ndx_d-1
  nfx_d1=nfx_d+1
  ndy_d1=ndy_d-1
  nfy_d1=nfy_d+1
  ndz_d1=ndz_d-1
  nfz_d1=nfz_d+1
  
  ! Definitions of nd_di and nf_di
  ! ----------------------------
  ! This index corresponds to the interior points where numerical
  ! dissipation is applied [in mod_artvisc]
  ! -> start at zero because evaluation at faces (n+1 faces for n points/nodes)
  ! -> +/- ngh in boundaries
  ndx_di= 0
  nfx_di=nx 
  if (is_boundary(1,1)) ndx_di=ndx_di+ngh
  if (is_boundary(1,2)) nfx_di=nfx_di-ngh
  ndy_di= 0
  nfy_di=ny
  if (is_boundary(2,1)) ndy_di=ndy_di+ngh
  if (is_boundary(2,2)) nfy_di=nfy_di-ngh
  ndz_di= 0
  nfz_di=nz
  if (is_boundary(3,1)) ndz_di=ndz_di+ngh
  if (is_boundary(3,2)) nfz_di=nfz_di-ngh

  ! Index bounds for characteristic conditions
  ! AND Riemann based
  ! ==========================================
  ! Definitions of nd_c and nf_c
  ! ----------------------------
  ! This index corresponds to the points where characteristic
  ! conditions are applied [in mod_charac]
  ! -> not advanced in wall BCs if characteristic BC on 1 point
  !    advanced if characteristic BC on 5 points (use of 1,ny)
  ! -> not advanced in corners in the j-direction (already
  !    advanced in the i-direction)
  ! -> advanced in ghost points with neighboor block
  ! CORRECTION: not advanced in ghost points
  ndx_c=1
  nfx_c=nx
  ndy_c=1
  nfy_c=ny
  if (BC_face(1,1)%sort.lt.0) then
     if (is_bc_1pt(1,1)) then
        ndx_c = ndx_c+1
        if (is_bc_wall(2,1)) then
           ndy_c = ndy_c+1
        else
           ndy_c = ndy_c!-ngh
        endif
        if (is_bc_wall(2,2)) then
           nfy_c = nfy_c-1
        else
           nfy_c = nfy_c!+ngh
        endif
     else
        ndx_c = ndx_c+ngh
     endif
  endif
  if (BC_face(1,2)%sort.lt.0) then
     if (is_bc_1pt(1,2)) then
        nfx_c = nfx_c-1
        if (is_bc_wall(2,1)) then
           ndy_c = ndy_c+1
        else
           ndy_c = ndy_c!-ngh
        endif
        if (is_bc_wall(2,2)) then
           nfy_c = nfy_c-1
        else
           nfy_c = nfy_c!+ngh
        endif
     else
        nfx_c = nfx_c-ngh
     endif
  endif
  if (BC_face(2,1)%sort.lt.0) then
     if (is_bc_1pt(2,1)) then
        ndy_c = ndy_c+1
        if (is_bc_wall(1,1)) then
           ndx_c = ndx_c+1
        else
           ndx_c = ndx_c!-ngh
        endif
        if (is_bc_wall(1,2)) then
           nfx_c = nfx_c-1
        else
           nfx_c = nfx_c!+ngh
        endif
     else
        ndy_c = ndy_c+ngh
     endif
  endif
  if (BC_face(2,2)%sort.lt.0) then
     if (is_bc_1pt(2,2)) then
        nfy_c = nfy_c-1
        if (is_bc_wall(1,1)) then
           ndx_c = ndx_c+1
        else
           ndx_c = ndx_c!-ngh
        endif
        if (is_bc_wall(1,2)) then
           nfx_c = nfx_c-1
        else
           nfx_c = nfx_c!+ngh
        endif
     else
        nfy_c = nfy_c-ngh
     endif
  endif


  ! Indices for TamDong on 1 point
  ! ==============================
  ! Definitions of ndx_td1m1 and nfx_td1p1
  ! --------------------------------------
  ! This index corresponds to the points where TamDong
  ! conditions are applied [in mod_TamDong2d_1pt]
  ! -> not advanced in corners if wall
  ! -> treated separetely in bc_TD2d_1pt_imin
  !    and bc_TD2d_1pt_imax if TamDong
  ndx_td1=1
  nfx_td1=nx
  if (BC_face(1,1)%sort.le.0) ndx_td1 = ndx_td1+2
  if (BC_face(1,2)%sort.le.0) nfx_td1 = nfx_td1-2
  ndy_td1=1
  nfy_td1=ny
  if (BC_face(2,1)%sort.le.0) ndy_td1 = ndy_td1+2
  if (BC_face(2,2)%sort.le.0) nfy_td1 = nfy_td1-2
  ndx_td1m1 = max(ndx_td1-1,1);nfx_td1p1 = min(nfx_td1+1,nx)
  ndy_td1m1 = max(ndy_td1-1,1);nfy_td1p1 = min(nfy_td1+1,ny)
  if (.not.is_bc_wall(2,1)) ndy_td1m1=1
  if (.not.is_bc_wall(2,2)) nfy_td1p1=ny


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! for non-regression purpose [in mod_artvisc/artvisc_o9o]
  ndx_f= 1
  nfx_f=nx 
  if (is_boundary(1,1)) ndx_f=ndx_f+ngh
  if (is_boundary(1,2)) nfx_f=nfx_f-ngh
  ndy_f= 1
  nfy_f=ny
  if (is_boundary(2,1)) ndy_f=ndy_f+ngh
  if (is_boundary(2,2)) nfy_f=nfy_f-ngh
  ndz_f= 1
  nfz_f=nz
  if (is_boundary(3,1)) ndz_f=ndz_f+ngh
  if (is_boundary(3,2)) nfz_f=nfz_f-ngh

  ! <- in MOD_ART_VISC *** TO BE CHANGED Nota used for g_ksi (grid_metrics_curv) & in mod_irs_(v2,v3)
  ! kept for compatibility with old routine [in mod_artvisc/artvisc_o9o]
  ndx1 =    0;   ndy1 =    0;   ndz1 =    0
  nfx1 = nx+1;   nfy1 = ny+1;   nfz1 = nz+1
  if (is_boundary(1,1)) ndx1 = ndx1+1
  if (is_boundary(1,2)) nfx1 = nfx1-1
  if (is_boundary(2,1)) ndy1 = ndy1+1
  if (is_boundary(2,2)) nfy1 = nfy1-1
  if (is_boundary(3,1)) ndz1 = ndz1+1
  if (is_boundary(3,2)) nfz1 = nfz1-1

  ! <- in MOD_CHARAC & in MOD_SPONGE  *** TO BE CHANGED
  ! charac not yet rewritten
  ndx3 =  1;   ndy3 =  1;   ndz3 = 1
  nfx3 = nx;   nfy3 = ny;   nfz3 = nz
  if (is_boundary(1,1)) ndx3 = ndx3+1
  if (is_boundary(1,2)) nfx3 = nfx3-1
  if (is_boundary(2,1)) ndy3 = ndy3+1
  if (is_boundary(2,2)) nfy3 = nfy3-1
  if (is_boundary(3,1)) ndz3 = ndz3+1
  if (is_boundary(3,2)) nfz3 = nfz3-1

  ! Modify all z indices for 2D runs
  ! ================================
  if (is_2d) then
     ! index bounds for inviscid fluxes
     ndz=1;  ndzt=1;  ndz_e=1
     nfz=1;  nfzt=1;  nfz_e=1

     ! index bounds for viscous fluxes
     ndz_v=1;  ndz_vi=1;  ndz_v1=1;  ndz_v2=1;  ndzt_v=1
     nfz_v=1;  nfz_vi=1;  nfz_v1=1;  nfz_v2=1;  nfzt_v=1
     
     ! index bounds for numerical dissipation
     ndz_d=1;  ndz_d1=1;  ndz_di=1
     nfz_d=1;  nfz_d1=1;  nfz_di=1
  endif

  ! ! Definitions of ndw and nfw
  ! ! --------------------------
  ! ! This index is used in mod_bc_wall_c
  ! !!! Modification CM mars 2023 !!!
  ! ! if the same proc has two walls on opposite frontiers, it does nd[]w+1 2 times
  ! ! ends up with e.g. for i in [1,nx]: i in [3,nx-2] instead of [2,nx-1]
  ! ndyw=1
  ! nfyw=ny
  ! if (is_bc_wall(1,1).and.BC_face(2,1)%sort<0) then
  !    if(is_bc_1pt(2,1)) ndyw=ndyw+1
  !    if(.not.is_bc_1pt(2,1)) ndyw=ndyw+ngh
  ! endif
  ! if (is_bc_wall(1,1).and.BC_face(2,2)%sort<0) then
  !    if(is_bc_1pt(2,2)) nfyw=nfyw-1
  !    if(.not.is_bc_1pt(2,2)) nfyw=nfyw-ngh
  ! endif
  ! if (is_bc_wall(1,2).and.BC_face(2,1)%sort<0.and.ndyw.eq.1) then
  !    if(is_bc_1pt(2,1)) ndyw=ndyw+1
  !    if(.not.is_bc_1pt(2,1)) ndyw=ndyw+ngh
  ! endif
  ! if (is_bc_wall(1,2).and.BC_face(2,2)%sort<0.and.nfyw.eq.ny) then
  !    if(is_bc_1pt(2,2)) nfyw=nfyw-1
  !    if(.not.is_bc_1pt(2,2)) nfyw=nfyw-ngh
  ! endif

  ! ndxw=1
  ! nfxw=nx
  ! if (is_bc_wall(2,1).and.BC_face(1,1)%sort<0) then
  !    if(is_bc_1pt(1,1)) ndxw=ndxw+1
  !    if(.not.is_bc_1pt(1,1)) ndxw=ndxw+ngh
  ! endif
  ! if (is_bc_wall(2,1).and.BC_face(1,2)%sort<0) then
  !    if(is_bc_1pt(1,2)) nfxw=nfxw-1
  !    if(.not.is_bc_1pt(1,2)) nfxw=nfxw-ngh
  ! endif
  ! if (is_bc_wall(2,2).and.BC_face(1,1)%sort<0.and.ndxw.eq.1) then
  !    if(is_bc_1pt(1,1)) ndxw=ndxw+1
  !    if(.not.is_bc_1pt(1,1)) ndxw=ndxw+ngh
  ! endif
  ! if (is_bc_wall(2,2).and.BC_face(1,2)%sort<0.and.nfxw.eq.nx) then
  !    if(is_bc_1pt(1,2)) nfxw=nfxw-1
  !    if(.not.is_bc_1pt(1,2)) nfxw=nfxw-ngh
  ! endif

  ! ! Other version where wall points are imposed
  ! ! even if the edge is shared with an inlet/outlet but
  ! ! only in the case of is_bc_1pt=.true.
  ! ! TEMPORARY
  ! ndyw=1
  ! nfyw=ny
  ! if (is_bc_wall(1,1).and.BC_face(2,1)%sort<0) then
  !    if(.not.is_bc_1pt(2,1)) ndyw=ndyw+ngh
  ! endif
  ! if (is_bc_wall(1,1).and.BC_face(2,2)%sort<0) then
  !    if(.not.is_bc_1pt(2,2)) nfyw=nfyw-ngh
  ! endif
  ! if (is_bc_wall(1,2).and.BC_face(2,1)%sort<0.and.ndyw.eq.1) then
  !    if(.not.is_bc_1pt(2,1)) ndyw=ndyw+ngh
  ! endif
  ! if (is_bc_wall(1,2).and.BC_face(2,2)%sort<0.and.nfyw.eq.ny) then
  !    if(.not.is_bc_1pt(2,2)) nfyw=nfyw-ngh
  ! endif

  ! ndxw=1
  ! nfxw=nx
  ! if (is_bc_wall(2,1).and.BC_face(1,1)%sort<0) then
  !    if(.not.is_bc_1pt(1,1)) ndxw=ndxw+ngh
  ! endif
  ! if (is_bc_wall(2,1).and.BC_face(1,2)%sort<0) then
  !    if(.not.is_bc_1pt(1,2)) nfxw=nfxw-ngh
  ! endif
  ! if (is_bc_wall(2,2).and.BC_face(1,1)%sort<0.and.ndxw.eq.1) then
  !    if(.not.is_bc_1pt(1,1)) ndxw=ndxw+ngh
  ! endif
  ! if (is_bc_wall(2,2).and.BC_face(1,2)%sort<0.and.nfxw.eq.nx) then
  !    if(.not.is_bc_1pt(1,2)) nfxw=nfxw-ngh
  ! endif

  ! Version AB - correction ?
  ! Wall points advanced everywhere
  ! If BC min (max), first (last) point not advanced in loop but aside
  ! Index ndw and nfw indicates loop bounds

!!$  ! comment XG: for jmin/jmax
!!$  ndxw=1
!!$  nfxw=nx
!!$  if (is_bc_wall(2,1).and.BC_face(1,1)%sort<0) ndxw=2
!!$  if (is_bc_wall(2,1).and.BC_face(1,2)%sort<0) nfxw=nx-1
!!$  if (is_bc_wall(2,2).and.BC_face(1,1)%sort<0) ndxw=2
!!$  if (is_bc_wall(2,2).and.BC_face(1,2)%sort<0) nfxw=nx-1
!!$
!!$  ! comment XG: for imin/imax
!!$  ndyw=1
!!$  nfyw=ny
!!$  if (is_bc_wall(1,1).and.BC_face(2,1)%sort<0) ndyw=2
!!$  if (is_bc_wall(1,1).and.BC_face(2,2)%sort<0) nfyw=ny-1
!!$  if (is_bc_wall(1,2).and.BC_face(2,1)%sort<0) ndyw=2
!!$  if (is_bc_wall(1,2).and.BC_face(2,2)%sort<0) nfyw=ny-1

  ! New version by XG
  ! Wall points advanced everywhere

  ! Indices for wall boundary conditions
  ! ====================================
  ! This index is used to compute corners between wall BCs
  ! and other BCs (wall, free, etc) with coordinate transform (curvilinear)
  ! -> the problem is the computation of pressure and temperature gradients
  !    in parallel directions (1 direction in 2D curv and 2 directions in 3D curv)
  ! -> for example at imin wall BC, for the parallel direction j (y)
  !    * the gradients use centered second-order between ndy_imin and nfy_imin
  !    * the gradients use non-centered first-order at j=1 and/or j=ny (if they are corners)

  ! wall BCs at i=cste
  ! ------------------
  ! at imin: parallel direction j (ndy/nfy)
  ndy_imin=1
  nfy_imin=ny
  if (is_bc_wall(1,1).and.BC_face(2,1)%sort.le.0) ndy_imin=2
  if (is_bc_wall(1,1).and.BC_face(2,2)%sort.le.0) nfy_imin=ny-1
  if (is_curv3) then
     ! at imin: parallel direction k (ndz/nfz) - only in 3D
     ndz_imin=1
     nfz_imin=nz
     if (is_bc_wall(1,1).and.BC_face(3,1)%sort.le.0) ndz_imin=2
     if (is_bc_wall(1,1).and.BC_face(3,2)%sort.le.0) nfz_imin=nz-1
  endif
  ! at imax: parallel direction j (ndy/nfy)
  ndy_imax=1
  nfy_imax=ny
  if (is_bc_wall(1,2).and.BC_face(2,1)%sort.le.0) ndy_imax=2
  if (is_bc_wall(1,2).and.BC_face(2,2)%sort.le.0) nfy_imax=ny-1
  if (is_curv3) then
     ! at imax: parallel direction k (ndz/nfz) - only in 3D
     ndz_imax=1
     nfz_imax=nz
     if (is_bc_wall(1,2).and.BC_face(3,1)%sort.le.0) ndz_imax=2
     if (is_bc_wall(1,2).and.BC_face(3,2)%sort.le.0) nfz_imax=nz-1
  endif

  ! wall BCs at j=cste
  ! ------------------
  ! at jmin: parallel direction i (ndx/nfx)
  ndx_jmin=1
  nfx_jmin=nx
  if (is_bc_wall(2,1).and.BC_face(1,1)%sort.le.0) ndx_jmin=2
  if (is_bc_wall(2,1).and.BC_face(1,2)%sort.le.0) nfx_jmin=nx-1
  if (is_curv3) then
     ! at jmin: parallel direction k (ndz/nfz) - only in 3D
     ndz_jmin=1
     nfz_jmin=nz
     if (is_bc_wall(2,1).and.BC_face(3,1)%sort.le.0) ndz_jmin=2
     if (is_bc_wall(2,1).and.BC_face(3,2)%sort.le.0) nfz_jmin=nz-1
  endif
  ! at jmax: parallel direction i (ndx/nfx)
  ndx_jmax=1
  nfx_jmax=nx
  if (is_bc_wall(2,2).and.BC_face(1,1)%sort.le.0) ndx_jmax=2
  if (is_bc_wall(2,2).and.BC_face(1,2)%sort.le.0) nfx_jmax=nx-1
  if (is_curv3) then
     ! at jmax: parallel direction k (ndz/nfz) - only in 3D
     ndz_jmax=1
     nfz_jmax=nz
     if (is_bc_wall(2,2).and.BC_face(3,1)%sort.le.0) ndz_jmax=2
     if (is_bc_wall(2,2).and.BC_face(3,2)%sort.le.0) nfz_jmax=nz-1
  endif

  ! wall BCs at k=cste
  ! ------------------
  if (is_curv3) then ! - only in 3D
     ! at kmin: parallel direction i (ndx/nfx)
     ndx_kmin=1
     nfx_kmin=nx
     if (is_bc_wall(3,1).and.BC_face(1,1)%sort.le.0) ndx_kmin=2
     if (is_bc_wall(3,1).and.BC_face(1,2)%sort.le.0) nfx_kmin=nx-1
     ! at kmin: parallel direction j (ndy/nfy)
     ndy_kmin=1
     nfy_kmin=ny
     if (is_bc_wall(3,1).and.BC_face(2,1)%sort.le.0) ndy_kmin=2
     if (is_bc_wall(3,1).and.BC_face(2,2)%sort.le.0) nfy_kmin=ny-1
     ! at kmax: parallel direction i (ndx/nfx)
     ndx_kmax=1
     nfx_kmax=nx
     if (is_bc_wall(3,2).and.BC_face(1,1)%sort.le.0) ndx_kmax=2
     if (is_bc_wall(3,2).and.BC_face(1,2)%sort.le.0) nfx_kmax=nx-1
     ! at kmax: parallel direction j (ndy/nfy)
     ndy_kmax=1
     nfy_kmax=ny
     if (is_bc_wall(3,2).and.BC_face(2,1)%sort.le.0) ndy_kmax=2
     if (is_bc_wall(3,2).and.BC_face(2,2)%sort.le.0) nfy_kmax=ny-1
  endif



  ! RANS indices
  ! ============

  if (is_RANS) then


     ! Index bounds for inviscid fluxes
     ! ================================

     ! Definitions of nd_e and nf_e
     ! ----------------------------
     ! This index corresponds to the points where Eulerian terms
     ! are advanced [in flux_euler_XX], by default between 1 and nx
     ! (for interior points).
     ! -> not advanced in BCs
     !    (1 cell layer for BCs coded on 1 point [is_bc_1pt] and
     !     5 cell layers for BCs coded on 5 points)
     ndx_e_r=1
     if (is_bc_1pt(1,1)) then
        ndx_e_r=ndx_e_r+1
     elseif ((BC_face(1,1)%sort==-1).or.(BC_face(1,1)%sort==-2)) then
        ndx_e_r=ndx_e_r!+ngh_r
     endif
     nfx_e_r=nx
     if (is_bc_1pt(1,2)) then
        nfx_e_r=nfx_e_r-1
     elseif ((BC_face(1,2)%sort==-1).or.(BC_face(1,2)%sort==-2)) then
        nfx_e_r=nfx_e_r!-ngh_r
     endif
     ndy_e_r=1
     if (is_bc_1pt(2,1)) then
        ndy_e_r=ndy_e_r+1
     elseif ((BC_face(2,1)%sort==-1).or.(BC_face(2,1)%sort==-2)) then
        ndy_e_r=ndy_e_r!+ngh_r
     endif
     nfy_e_r=ny
     if (is_bc_1pt(2,2)) then
        nfy_e_r=nfy_e_r-1
     elseif ((BC_face(2,2)%sort==-1).or.(BC_face(2,2)%sort==-2)) then
        nfy_e_r=nfy_e_r!-ngh_r
     endif
     ndz_e_r=1
     if (is_bc_1pt(3,1)) then
        ndz_e_r=ndz_e_r+1
     elseif ((BC_face(3,1)%sort==-1).or.(BC_face(3,1)%sort==-2)) then
        ndz_e_r=ndz_e!+ngh_r
     endif
     nfz_e_r=nz
     if (is_bc_1pt(3,2)) then
        nfz_e_r=nfz_e_r-1
     elseif ((BC_face(3,2)%sort==-1).or.(BC_face(3,2)%sort==-2)) then
        nfz_e_r=nfz_e_r!-ngh_r
     endif


     ! Definitions of nd and nf
     ! ------------------------
     ! This index corresponds to the points where Eulerian terms
     ! are advanced [in flux_euler_XX] for interior points
     ! for which boundary points are already advanced.
     ndx_r=1
     nfx_r=nx
     if (is_boundary(1,1)) ndx_r= ndx_r+ngh_r
     if (is_boundary(1,2)) nfx_r= nfx_r-ngh_r
     ndy_r=1
     nfy_r=ny
     if (is_boundary(2,1)) ndy_r= ndy_r+ngh_r
     if (is_boundary(2,2)) nfy_r= nfy_r-ngh_r
     ndz_r=1
     nfz_r=nz
     if (is_boundary(3,1)) ndz_r= ndz_r+ngh_r
     if (is_boundary(3,2)) nfz_r= nfz_r-ngh_r

     ! Definitions of ndt and nft
     ! --------------------------
     ! This index is the index of interior points
     ! augmented by the number of ghost cells (-/+ ngh).
     ndxt_r=ndx_r-ngh_r
     nfxt_r=nfx_r+ngh_r
     ndyt_r=ndy_r-ngh_r
     nfyt_r=nfy_r+ngh_r
     ndzt_r=ndz_r-ngh_r
     nfzt_r=nfz_r+ngh_r

     ! Index bounds for viscous fluxes
     ! ===============================
     ! Definitions of nd_v and nf_v
     ! ----------------------------
     ! This index corresponds to the points where viscous terms
     ! are advanced [in flux_visc_XX], by default between 1 and nx
     ! (for interior points and wall BCs).
     ! -> not advanced in non-reflecting BCs
     !    (1 cell layer for BCs coded on 1 point [is_bc_1pt] and
     !     5 cell layers for BCs coded on 5 points)
     ndx_v_r=1
     if ((is_bc_1pt(1,1)).and.(.not.is_bc_wall(1,1))) then
        ndx_v_r= ndx_v_r+1
     elseif ((BC_face(1,1)%sort==-1).or.(BC_face(1,1)%sort==-2).or.(BC_face(1,1)%sort==-4)) then
        ndx_v_r= ndx_v_r!+ngh_r
     endif
     nfx_v_r=nx
     if ((is_bc_1pt(1,2)).and.(.not.is_bc_wall(1,2))) then
        nfx_v_r= nfx_v_r-1
     elseif ((BC_face(1,2)%sort==-1).or.(BC_face(1,2)%sort==-2).or.(BC_face(1,2)%sort==-4)) then
        nfx_v_r= nfx_v_r!-ngh_r
     endif
     ndy_v_r=1
     if ((is_bc_1pt(2,1)).and.(.not.is_bc_wall(2,1))) then
        ndy_v_r= ndy_v_r+1
     elseif ((BC_face(2,1)%sort==-1).or.(BC_face(2,1)%sort==-2).or.(BC_face(2,1)%sort==-4)) then
        ndy_v_r= ndy_v_r!+ngh_r
     endif
     nfy_v_r=ny
     if ((is_bc_1pt(2,2)).and.(.not.is_bc_wall(2,2))) then
        nfy_v_r= nfy_v_r-1
     elseif ((BC_face(2,2)%sort==-1).or.(BC_face(2,2)%sort==-2).or.(BC_face(2,2)%sort==-4)) then
        nfy_v_r= nfy_v_r!-ngh_r
     endif
     ndz_v_r=1
     if ((is_bc_1pt(3,1)).and.(.not.is_bc_wall(3,1))) then
        ndz_v_r= ndz_v_r+1
     elseif ((BC_face(3,1)%sort==-1).or.(BC_face(3,1)%sort==-2).or.(BC_face(3,1)%sort==-4)) then
        ndz_v_r= ndz_v_r!+ngh_r
     endif
     nfz_v_r=nz
     if ((is_bc_1pt(3,2)).and.(.not.is_bc_wall(3,2))) then
        nfz_v_r= nfz_v_r-1
     elseif ((BC_face(3,2)%sort==-1).or.(BC_face(3,2)%sort==-2).or.(BC_face(3,2)%sort==-4)) then
        nfz_v_r= nfz_v_r!-ngh_r
     endif

     ! Definitions of nd_vi and nf_vi
     ! ------------------------------
     ! This index corresponds to the points where viscous terms
     ! are advanced [in flux_visc_XX] for interior points
     ! for which boundary points are already advanced.
     ! -> for bc_wall, it is nd_v+ngh_v or nf_v-ngh_v
     ! -> for bc_1pt , it is nd_v+ngh_v-1 or nf_v-ngh_v+1
     ! -> for bc on 5 pts, nd_v or nf_v unchanged
     ndx_vi_r=ndx_v_r
     nfx_vi_r=nfx_v_r
     if (is_bc_wall(1,1)) ndx_vi_r= ndx_vi_r+ngh_v_r
     if (is_bc_wall(1,2)) nfx_vi_r= nfx_vi_r-ngh_v_r
     if ((is_bc_1pt(1,1)).and.(.not.is_bc_wall(1,1)))  ndx_vi_r= ndx_vi_r+ngh_v_r-1
     if ((is_bc_1pt(1,2)).and.(.not.is_bc_wall(1,2)))  nfx_vi_r= nfx_vi_r-ngh_v_r+1
     ndy_vi_r=ndy_v_r
     nfy_vi_r=nfy_v_r
     if (is_bc_wall(2,1)) ndy_vi_r= ndy_vi_r+ngh_v_r
     if (is_bc_wall(2,2)) nfy_vi_r= nfy_vi_r-ngh_v_r
     if ((is_bc_1pt(2,1)).and.(.not.is_bc_wall(2,1)))  ndy_vi_r= ndy_vi_r+ngh_v_r-1
     if ((is_bc_1pt(2,2)).and.(.not.is_bc_wall(2,2)))  nfy_vi_r= nfy_vi_r-ngh_v_r+1
     ndz_vi_r=ndz_v_r
     nfz_vi_r=nfz_v_r
     if (is_bc_wall(3,1)) ndz_vi_r= ndz_vi_r+ngh_v_r
     if (is_bc_wall(3,2)) nfz_vi_r= nfz_vi_r-ngh_v_r
     if ((is_bc_1pt(3,1)).and.(.not.is_bc_wall(3,1)))  ndz_vi_r= ndz_vi_r+ngh_v_r-1
     if ((is_bc_1pt(3,2)).and.(.not.is_bc_wall(3,2)))  nfz_vi_r= nfz_vi_r-ngh_v_r+1

     ! Definitions of ndt_v and nft_v
     ! -------------------------------
     ! This index is the index of interior points
     ! augmented by the number of ghost cells (-/+ ngh_v).
     ndxt_v_r=ndx_vi_r-ngh_v_r
     nfxt_v_r=nfx_vi_r+ngh_v_r
     ndyt_v_r=ndy_vi_r-ngh_v_r
     nfyt_v_r=nfy_vi_r+ngh_v_r
     ndzt_v_r=ndz_vi_r-ngh_v_r
     nfzt_v_r=nfz_vi_r+ngh_v_r

     ! Definitions of ndt_v1 and nft_v1 & ndt_v2 and nft_v2
     ! ----------------------------------------------------
     ! [indices for double derivatives in grad_velT]

     ndx_v1_r= 1-ngh_v_r
     nfx_v1_r=nx+ngh_v_r
     if (is_boundary(1,1)) ndx_v1_r=ndx_v1_r+ngh_v_r
     if (is_boundary(1,2)) nfx_v1_r=nfx_v1_r-ngh_v_r
     ndy_v1_r= 1-ngh_v_r
     nfy_v1_r=ny+ngh_v_r
     if (is_boundary(2,1)) ndy_v1_r=ndy_v1_r+ngh_v_r
     if (is_boundary(2,2)) nfy_v1_r=nfy_v1_r-ngh_v_r
     ndz_v1_r= 1-ngh_v_r
     nfz_v1_r=nz+ngh_v_r
     if (is_boundary(3,1)) ndz_v1_r=ndz_v1_r+ngh_v_r
     if (is_boundary(3,2)) nfz_v1_r=nfz_v1_r-ngh_v_r

     ndx_v2_r= 1-ngh_v_r
     nfx_v2_r=nx+ngh_v_r
     if (is_boundary(1,1)) ndx_v2_r=ndx_v2_r+2*ngh_v_r
     if (is_boundary(1,2)) nfx_v2_r=nfx_v2_r-2*ngh_v_r
     ndy_v2_r= 1-ngh_v_r
     nfy_v2_r=ny+ngh_v_r
     if (is_boundary(2,1)) ndy_v2_r=ndy_v2_r+2*ngh_v_r
     if (is_boundary(2,2)) nfy_v2_r=nfy_v2_r-2*ngh_v_r
     ndz_v2_r= 1-ngh_v_r
     nfz_v2_r=nz+ngh_v_r
     if (is_boundary(3,1)) ndz_v2_r=ndz_v2_r+2*ngh_v_r
     if (is_boundary(3,2)) nfz_v2_r=nfz_v2_r-2*ngh_v_r

     ndx_v3_r= 1
     nfx_v3_r=nx
     if (is_boundary(1,1)) ndx_v3_r=ndx_v3_r+ngh_v_r
     if (is_boundary(1,2)) nfx_v3_r=nfx_v3_r-ngh_v_r
     ndy_v3_r= 1
     nfy_v3_r=ny
     if (is_boundary(2,1)) ndy_v3_r=ndy_v3_r+ngh_v_r
     if (is_boundary(2,2)) nfy_v3_r=nfy_v3_r-ngh_v_r
     ndz_v3_r= 1
     nfz_v3_r=nz
     if (is_boundary(3,1)) ndz_v3_r=ndz_v3_r+ngh_v_r
     if (is_boundary(3,2)) nfz_v3_r=nfz_v3_r-ngh_v_r

     ! Index bounds for numerical dissipation
     ! ======================================

     ! Definitions of nd_d and nf_d
     ! ----------------------------
     ! This index corresponds to the points where numerical dissipation
     ! is applied [in mod_artvisc], by default between 1 and nx.
     ! -> exclude first and last points at boundaries
     ! -> possibly exclude 3 and 5 pts stencils (here coeff are set to zero)
     ndx_d_r= 1
     nfx_d_r=nx
     if (is_boundary(1,1)) ndx_d_r=ndx_d_r+1
     if (is_boundary(1,2)) nfx_d_r=nfx_d_r-1
     ndy_d_r= 1
     nfy_d_r=ny
     if (is_boundary(2,1)) ndy_d_r=ndy_d_r+1
     if (is_boundary(2,2)) nfy_d_r=nfy_d_r-1
     ndz_d_r= 1
     nfz_d_r=nz
     if (is_boundary(3,1)) ndz_d_r=ndz_d_r+1
     if (is_boundary(3,2)) nfz_d_r=nfz_d_r-1

     ! Definitions of nd_d1 and nf_d1
     ! ----------------------------
     ! extended index for sensor calculation (shock-capturing)
     ndx_d1_r=ndx_d_r-1
     nfx_d1_r=nfx_d_r+1
     ndy_d1_r=ndy_d_r-1
     nfy_d1_r=nfy_d_r+1
     ndz_d1_r=ndz_d_r-1
     nfz_d1_r=nfz_d_r+1

     ! Definitions of nd_di_r and nf_di_r
     ! ----------------------------
     ! This index corresponds to the interior points where numerical
     ! dissipation is applied [in mod_artvisc] in the case of RANS where the stencil
     ! is the same for turbulent variable Euler and viscous fluxes, but smaller than
     ! the one for the mean field
     ! -> start at zero because evaluation at faces (n+1 faces for n points/nodes)
     ! -> +/- ngh_v in boundaries
     ndx_di_r= 0
     nfx_di_r=nx
     if (is_boundary(1,1)) ndx_di_r=ndx_di_r+ngh_r
     if (is_boundary(1,2)) nfx_di_r=nfx_di_r-ngh_r
     ndy_di_r= 0
     nfy_di_r=ny
     if (is_boundary(2,1)) ndy_di_r=ndy_di_r+ngh_r
     if (is_boundary(2,2)) nfy_di_r=nfy_di_r-ngh_r
     ndz_di_r= 0
     nfz_di_r=nz
     if (is_boundary(3,1)) ndz_di_r=ndz_di_r+ngh_r
     if (is_boundary(3,2)) nfz_di_r=nfz_di_r-ngh_r

     ! Index bounds for characteristic conditions
     ! ==========================================
     ! Definitions of nd_c and nf_c
     ! ----------------------------
     ! This index corresponds to the points where characteristic
     ! conditions are applied [in mod_charac]
     ! -> not advanced in wall BCs
     ! -> not advanced in corners in the j-direction (already
     !    advanced in the i-direction)
!     ndx_c=1
!     nfx_c=nx
!     if (BC_face(1,1)%sort.le.0) ndx_c = ndx_c+1
!     if (BC_face(1,2)%sort.le.0) nfx_c = nfx_c-1
!     ndy_c=1
!     nfy_c=ny
!     if (BC_face(2,1)%sort.eq.0) ndy_c = ndy_c+1
!     if (BC_face(2,2)%sort.eq.0) nfy_c = nfy_c-1
!     ndz_c=1
!     nfz_c=nz


     ! Modify all z indices for 2D runs
     ! ================================
     if (is_2d) then
        ! index bounds for inviscid fluxes
        ndz_r=1;  ndzt_r=1;  ndz_e_r=1
        nfz_r=1;  nfzt_r=1;  nfz_e_r=1

        ! index bounds for viscous fluxes
        ndz_v_r=1;  ndz_vi_r=1;  ndz_v1_r=1;  ndz_v2_r=1;  ndzt_v_r=1
        nfz_v_r=1;  nfz_vi_r=1;  nfz_v1_r=1;  nfz_v2_r=1;  nfzt_v_r=1

        ! index bounds for numerical dissipation
        ndz_d_r=1;  ndz_d1_r=1;  ndz_di_r=1
        nfz_d_r=1;  nfz_d1_r=1;  nfz_di_r=1
     endif

     ! Index bounds for source terms
     ! =============================
     ! apply everywhere in interior except at wall
     ndx_s_r=1
     nfx_s_r=nx
     ndy_s_r=1
     nfy_s_r=ny
     ndz_s_r=1
     nfz_s_r=nz
     if (is_bc_wall(1,1)) ndx_s_r=ndx_s_r+1
     if (is_bc_wall(1,2)) nfx_s_r=nfx_s_r-1
     if (is_bc_wall(2,1)) ndy_s_r=ndy_s_r+1
     if (is_bc_wall(2,2)) nfy_s_r=nfy_s_r-1


  endif

end subroutine mpi_index_bounds
