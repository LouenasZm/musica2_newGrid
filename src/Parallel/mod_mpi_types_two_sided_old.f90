!==================================================================================
module mod_mpi_types_two_sided_old
!==================================================================================
  !> Module to create MPI types ** OLD VERSION **
!==================================================================================
  use mod_mpi_part
  implicit none
  ! -------------------------------------------------------------------------------
  ! type face (ngh layers of cells) for each direction [flux_euler & cie]
  integer :: type_faceW,type_faceE
  integer :: type_faceS,type_faceN
  integer :: type_faceF,type_faceB
  integer :: type_mW,type_mE,type_mS,type_mN ! for metrics
  ! -------------------------------------------------------------------------------
  ! type face (ngh_v layers of cells) for each direction [flux_visc & cie]
  integer :: type_faceW_v,type_faceE_v
  integer :: type_faceS_v,type_faceN_v
  integer :: type_faceF_v,type_faceB_v
  integer :: type_mW_v,type_mE_v,type_mS_v,type_mN_v ! for metrics
  ! -------------------------------------------------------------------------------
  ! type face (nghirs layers of cells) for each direction [increments for IRS]
  integer :: type_faceWinc,type_faceEinc
  integer :: type_faceSinc,type_faceNinc
  integer :: type_faceFinc,type_faceBinc
!!$  ! -------------------------------------------------------------------------------
!!$  ! type face+edges (ngh layers of cells) for each direction [NOT USED ANYMORE]
!!$  integer :: type_e_faceW,type_e_faceE
!!$  integer :: type_e_faceS,type_e_faceN
!!$  integer :: type_e_faceF,type_e_faceB
!!$  ! -------------------------------------------------------------------------------
  
contains

  !===============================================================
  subroutine mpi_types_comm_old
  !===============================================================
    !> Definition of MPI types for communications ** OLD VERSION **
  !===============================================================
    use mod_grid_directions
    implicit none
    ! ------------------------------------------------------------
    integer :: nex,ney ! extended sizes
    integer :: sign_i,sign_j
    integer :: memjump ! jump in memory between block starts
    ! strides in memory between block starts
    integer(kind=MPI_ADDRESS_KIND) :: stride,stride2
    integer :: type_base,type_edge
    ! ------------------------------------------------------------

    ! Extended sizes (+ ghost cells)
    ! ------------------------------
    nex=nx+2*ngh
    ney=ny+2*ngh

    ! stride between edges along x or y
    ! ---------------------------------
    stride2=nex*ney*sizeofreal

    ! MPI-type construction for imin face "along y" (W: west)
    ! -------------------------------------------------------
    ! for imin/W, parallel dir is j and normal dir is i
    inW=1
    inWs=1
    inwr=1-ngh
    if (neighbor(nW).ne.MPI_PROC_NULL) then
       if (nob(neighbor(nW))/=nob(iproc)) inWs=2
       if ((nob(neighbor(nW))==nob(iproc)).and.(coord(1)==0)) inWs=2
    endif
    ipW=1  
    sign_i=1
    sign_j=1
    ! if rev parallel -> +j => -j
    if (is_rev(1,1)) then
       sign_j=-1
       ipW=ny
    endif
    ! if rev  normal  -> +i => -i
    if (is_rev(1,2)) then
       sign_i=-1
       inW=ngh
       inWs=ngh+1
       inwr=0
!        if(nob(neighbor(nW))/=nob(iproc)) inWs=ngh+1
    endif
    if (is_swapij(1)) then
       ! if swap invert role of i and j
       memjump=sign_j*nex  ! j
       stride =sign_i*sizeofreal ! i
       call MPI_TYPE_VECTOR(ny,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride,type_base,type_mW,info)
    else
       ! default
       memjump=sign_i ! i
       stride =sign_j*nex*sizeofreal ! j
       call MPI_TYPE_VECTOR(ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ny,1,stride,type_base,type_mW,info)
    endif
    call MPI_TYPE_COMMIT(type_mW,info)
    call MPI_TYPE_CREATE_HVECTOR(nz,1,stride2,type_mW,type_faceW,info)
    call MPI_TYPE_COMMIT(type_faceW,info)

    ! MPI-type construction for imax face "along y" (E: east)
    ! -------------------------------------------------------
    ! for imax/E, parallel dir is j and normal dir is i
    inE=nx+1
    inEs=nx+1-ngh
    inEr=nx+1
    if (neighbor(nE).ne.MPI_PROC_NULL) then
       if (nob(neighbor(nE))/=nob(iproc)) inEs=nx-ngh
       if ((nob(neighbor(nE))==nob(iproc)).and.(coord(1)==ndomx-1)) inEs=nx-ngh
    endif
    ipE=1
    sign_i=1
    sign_j=1
    ! if rev parallel -> +j => -j
    if (is_rev(2,1)) then
       sign_j=-1
       ipE=ny
    endif
    ! if rev  normal  -> +i => -i
    if (is_rev(2,2)) then
       sign_i=-1
       inE=nx+ngh
       inEs=nx-1
       inEr=nx+ngh
!        if(nob(neighbor(nE))/=nob(iproc)) inEs=nx-1
    endif
    if (is_swapij(2)) then
       ! if swap invert role of i and j
       memjump=sign_j*nex  ! j
       stride =sign_i*sizeofreal ! i
       call MPI_TYPE_VECTOR(ny,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride,type_base,type_mE,info)
    else
       ! default
       memjump=sign_i ! i
       stride =sign_j*nex*sizeofreal ! j
       call MPI_TYPE_VECTOR(ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ny,1,stride,type_base,type_mE,info)
    endif
    call MPI_TYPE_COMMIT(type_mE,info)
    call MPI_TYPE_CREATE_HVECTOR(nz,1,stride2,type_mE,type_faceE,info)
    call MPI_TYPE_COMMIT(type_faceE,info)

    ! MPI-type construction for jmin face "along x" (S: south)
    ! --------------------------------------------------------
    ! for jmin/S, parallel dir is i and normal dir is j
    ipS=1
    inS=1
    inSs=1
    inSr=1-ngh
    if (neighbor(nS).ne.MPI_PROC_NULL) then
       if (nob(neighbor(nS))/=nob(iproc)) inSs=2
       if ((nob(neighbor(nS))==nob(iproc)).and.(coord(2)==0)) inSs=2
    endif
    sign_i=1
    sign_j=1
    ! if rev parallel -> +i => -i
    if (is_rev(3,1)) then
       sign_i=-1
       ipS=nx
    endif
    ! if rev  normal  -> +j => -j
    if (is_rev(3,2)) then
       sign_j=-1
       inS=ngh
       inSs=ngh+1
       inSr=0
!        if(nob(neighbor(nS))/=nob(iproc)) inSs=ngh+1
    endif
    if (is_swapij(3)) then
       ! if swap invert role of i and j
       memjump=sign_j*nex  ! j
       stride =sign_i*sizeofreal ! i
       call MPI_TYPE_VECTOR(ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(nx,1,stride,type_base,type_mS,info)
    else
       ! default
       memjump=sign_i ! i
       stride =sign_j*nex*sizeofreal ! j
       call MPI_TYPE_VECTOR(nx,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride,type_base,type_mS,info)
    endif
    call MPI_TYPE_COMMIT(type_mS,info)
    call MPI_TYPE_CREATE_HVECTOR(nz,1,stride2,type_mS,type_faceS,info)
    call MPI_TYPE_COMMIT(type_faceS,info)

    ! MPI-type construction for jmax face "along x" (N: north)
    ! --------------------------------------------------------
    ! for jmax/N, parallel dir is i and normal dir is j
    ipN=1
    inN=ny+1
    inNs=ny+1-ngh
    inNr=ny+1
    if (neighbor(nN).ne.MPI_PROC_NULL) then
       if (nob(neighbor(nN))/=nob(iproc)) inNs=ny-ngh
       if ((nob(neighbor(nN))==nob(iproc)).and.(coord(2)==ndomy-1)) inNs=ny-ngh
    endif
    sign_i=1
    sign_j=1
    ! if rev parallel -> +i => -i
    if (is_rev(4,1)) then
       sign_i=-1
       ipN=nx
    endif
    ! if rev  normal  -> +j => -j
    if (is_rev(4,2)) then
       sign_j=-1
       inN=ny+ngh
       inNs=ny-1
       inNr=ny+ngh
!        if(nob(neighbor(nN))/=nob(iproc)) inNs=ny-1
    endif
    if (is_swapij(4)) then
       ! if swap invert role of i and j
       memjump=sign_j*nex  ! j
       stride =sign_i*sizeofreal ! i
       call MPI_TYPE_VECTOR(ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(nx,1,stride,type_base,type_mN,info)
    else
       ! default
       memjump=sign_i ! i
       stride =sign_j*nex*sizeofreal ! j
       call MPI_TYPE_VECTOR(nx,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride,type_base,type_mN,info)
    endif
    call MPI_TYPE_COMMIT(type_mN,info)
    call MPI_TYPE_CREATE_HVECTOR(nz,1,stride2,type_mN,type_faceN,info)
    call MPI_TYPE_COMMIT(type_faceN,info)

    stride=nex*sizeofreal
    call MPI_TYPE_VECTOR(nx,1,1,MPI_DOUBLE_PRECISION,type_base,info)

    ! MPI-type construction for kmin face "along z" (F: forward)
    ! --------------------------------------------------------
    call MPI_TYPE_CREATE_HVECTOR( ny,1, stride,type_base,type_edge,info)
    call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride2,type_edge,type_faceF,info)
    call MPI_TYPE_COMMIT(type_faceF,info)
    
    ! MPI-type construction for kmax face "along z" (B: backward)
    ! -----------------------------------------------------------
    call MPI_TYPE_CREATE_HVECTOR( ny,1, stride,type_base,type_edge,info)
    call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride2,type_edge,type_faceB,info)
    call MPI_TYPE_COMMIT(type_faceB,info)

  end subroutine mpi_types_comm_old

  !===============================================================
  subroutine mpi_types_comm_v_old
  !===============================================================
    !> Definition of MPI types for communications
    !> of viscous terms ** OLD VERSION **
  !===============================================================
    use mod_grid_directions
    implicit none
    ! ------------------------------------------------------------
    integer :: nex,ney ! extended sizes
    integer :: sign_i,sign_j
    integer :: memjump ! jump in memory between block starts
    ! strides in memory between block starts
    integer(kind=MPI_ADDRESS_KIND) :: stride,stride2
    integer :: type_base,type_edge
    ! ------------------------------------------------------------

    ! MPI types for faces (ngh_v ghost cells)
    ! ===================
    
    ! Extended sizes (+ ghost cells)
    ! ------------------------------
    nex=nx+2*ngh_v
    ney=ny+2*ngh_v

    ! stride between edges along x or y
    ! ---------------------------------
    stride2=nex*ney*sizeofreal

    ! MPI-type construction for imin face "along y" (W: west)
    ! -------------------------------------------------------
    ! for imin/W, parallel dir is j and normal dir is i
    inW_v=1
    inWs_v=1
    inwr_v=1-ngh_v
    if(neighbor(nW).ne.MPI_PROC_NULL) then
      if(nob(neighbor(nW))/=nob(iproc)) inWs_v=2
    endif
    ipW=1  
    sign_i=1
    sign_j=1
    ! if rev parallel -> +j => -j
    if (is_rev(1,1)) then
       sign_j=-1
       ipW=ny
    endif
    ! if rev  normal  -> +i => -i
    if (is_rev(1,2)) then
       sign_i=-1
       inW_v=ngh_v
       inWs_v=ngh_v+1
       inwr_v=0
!        if(nob(neighbor(nW))/=nob(iproc)) inWs_v=ngh_v+1
    endif
    if (is_swapij(1)) then
       ! if swap invert role of i and j
       memjump=sign_j*nex  ! j
       stride =sign_i*sizeofreal ! i
       call MPI_TYPE_VECTOR(ny,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ngh_v,1,stride,type_base,type_mW_v,info)
    else
       ! default
       memjump=sign_i ! i
       stride =sign_j*nex*sizeofreal ! j
       call MPI_TYPE_VECTOR(ngh_v,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ny,1,stride,type_base,type_mW_v,info)
    endif
    call MPI_TYPE_COMMIT(type_mW_v,info)
    call MPI_TYPE_CREATE_HVECTOR(nz,1,stride2,type_mW_v,type_faceW_v,info)
    call MPI_TYPE_COMMIT(type_faceW_v,info)

    ! MPI-type construction for imax face "along y" (E: east)
    ! -------------------------------------------------------
    ! for imax/E, parallel dir is j and normal dir is i
    inE_v=nx+1
    inEs_v=nx+1-ngh_v
    inEr_v=nx+1
    if(neighbor(nE).ne.MPI_PROC_NULL) then
      if(nob(neighbor(nE))/=nob(iproc)) inEs_v=nx-ngh_v
    endif
    ipE=1
    sign_i=1
    sign_j=1
    ! if rev parallel -> +j => -j
    if (is_rev(2,1)) then
       sign_j=-1
       ipE=ny
    endif
    ! if rev  normal  -> +i => -i
    if (is_rev(2,2)) then
       sign_i=-1
       inE_v=nx+ngh_v
       inEs_v=nx-1
       inEr_v=nx+ngh_v
!        if(nob(neighbor(nE))/=nob(iproc)) inEs_v=nx-1
    endif
    if (is_swapij(2)) then
       ! if swap invert role of i and j
       memjump=sign_j*nex  ! j
       stride =sign_i*sizeofreal ! i
       call MPI_TYPE_VECTOR(ny,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ngh_v,1,stride,type_base,type_mE_v,info)
    else
       ! default
       memjump=sign_i ! i
       stride =sign_j*nex*sizeofreal ! j
       call MPI_TYPE_VECTOR(ngh_v,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ny,1,stride,type_base,type_mE_v,info)
    endif
    call MPI_TYPE_COMMIT(type_mE_v,info)
    call MPI_TYPE_CREATE_HVECTOR(nz,1,stride2,type_mE_v,type_faceE_v,info)
    call MPI_TYPE_COMMIT(type_faceE_v,info)

    ! MPI-type construction for jmin face "along x" (S: south)
    ! --------------------------------------------------------
    ! for jmin/S, parallel dir is i and normal dir is j
    ipS=1
    inS_v=1
    inSs_v=1
    inSr_v=1-ngh_v
    if(neighbor(nS).ne.MPI_PROC_NULL) then
      if(nob(neighbor(nS))/=nob(iproc)) inSs_v=2
    endif
    sign_i=1
    sign_j=1
    ! if rev parallel -> +i => -i
    if (is_rev(3,1)) then
       sign_i=-1
       ipS=nx
    endif
    ! if rev  normal  -> +j => -j
    if (is_rev(3,2)) then
       sign_j=-1
       inS_v=ngh_v
       inSs_v=ngh_v+1
       inSr_v=0
!        if(nob(neighbor(nS))/=nob(iproc)) inSs_v=ngh_v+1
    endif
    if (is_swapij(3)) then
       ! if swap invert role of i and j
       memjump=sign_j*nex  ! j
       stride =sign_i*sizeofreal ! i
       call MPI_TYPE_VECTOR(ngh_v,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(nx,1,stride,type_base,type_mS_v,info)
    else
       ! default
       memjump=sign_i ! i
       stride =sign_j*nex*sizeofreal ! j
       call MPI_TYPE_VECTOR(nx,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ngh_v,1,stride,type_base,type_mS_v,info)
    endif
    call MPI_TYPE_COMMIT(type_mS_v,info)
    call MPI_TYPE_CREATE_HVECTOR(nz,1,stride2,type_mS_v,type_faceS_v,info)
    call MPI_TYPE_COMMIT(type_faceS_v,info)

    ! MPI-type construction for jmax face "along x" (N: north)
    ! --------------------------------------------------------
    ! for jmax/N, parallel dir is i and normal dir is j
    ipN=1
    inN_v=ny+1
    inNs_v=ny+1-ngh_v
    inNr_v=ny+1
    if(neighbor(nN).ne.MPI_PROC_NULL) then
      if(nob(neighbor(nN))/=nob(iproc)) inNs_v=ny-ngh_v
    endif
    sign_i=1
    sign_j=1
    ! if rev parallel -> +i => -i
    if (is_rev(4,1)) then
       sign_i=-1
       ipN=nx
    endif
    ! if rev  normal  -> +j => -j
    if (is_rev(4,2)) then
       sign_j=-1
       inN_v=ny+ngh_v
       inNs_v=ny-1
       inNr_v=ny+ngh_v
!        if(nob(neighbor(nN))/=nob(iproc)) inNs_v=ny-1
    endif
    if (is_swapij(4)) then
       ! if swap invert role of i and j
       memjump=sign_j*nex  ! j
       stride =sign_i*sizeofreal ! i
       call MPI_TYPE_VECTOR(ngh_v,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(nx,1,stride,type_base,type_mN_v,info)
    else
       ! default
       memjump=sign_i ! i
       stride =sign_j*nex*sizeofreal ! j
       call MPI_TYPE_VECTOR(nx,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ngh_v,1,stride,type_base,type_mN_v,info)
    endif
    call MPI_TYPE_COMMIT(type_mN_v,info)
    call MPI_TYPE_CREATE_HVECTOR(nz,1,stride2,type_mN_v,type_faceN_v,info)
    call MPI_TYPE_COMMIT(type_faceN_v,info)

    stride=nex*sizeofreal
    call MPI_TYPE_VECTOR(nx,1,1,MPI_DOUBLE_PRECISION,type_base,info)

    ! MPI-type construction for kmin face "along z" (F: forward)
    ! --------------------------------------------------------
    call MPI_TYPE_CREATE_HVECTOR(   ny,1, stride,type_base,type_edge,info)
    call MPI_TYPE_CREATE_HVECTOR(ngh_v,1,stride2,type_edge,type_faceF_v,info)
    call MPI_TYPE_COMMIT(type_faceF_v,info)
    
    ! MPI-type construction for kmax face "along z" (B: backward)
    ! -----------------------------------------------------------
    call MPI_TYPE_CREATE_HVECTOR(   ny,1, stride,type_base,type_edge,info)
    call MPI_TYPE_CREATE_HVECTOR(ngh_v,1,stride2,type_edge,type_faceB_v,info)
    call MPI_TYPE_COMMIT(type_faceB_v,info)

  end subroutine mpi_types_comm_v_old

  !===============================================================
  subroutine mpi_types_comm_inc_old
  !===============================================================
    !> Definition of MPI types for communications
    !> of increments (nghirs ghost cells) ** OLD VERSION **
  !===============================================================
    use mod_grid_directions
    implicit none
    ! ------------------------------------------------------------
    integer :: nex,ney ! extended sizes
    integer :: sign_i,sign_j
    integer :: memjump ! jump in memory between block starts
    ! strides in memory between block starts
    integer(kind=MPI_ADDRESS_KIND) :: stride,stride2
    integer :: type_base,type_edge
    ! ------------------------------------------------------------
   
    ! Extended sizes (+ ghost cells)
    ! ------------------------------
    nex=nx+2*nghirs
    ney=ny+2*nghirs

    ! stride between edges along x or y
    ! ---------------------------------
    stride2=nex*ney*sizeofreal

    ! MPI-type construction for imin face "along y" (W: west)
    ! -------------------------------------------------------
    ! for imin/W, parallel dir is j and normal dir is i
    inWi=1
    inWsi=1
    inWri=1-nghirs
    if (neighbor(nW).ne.MPI_PROC_NULL) then
       if (nob(neighbor(nW))/=nob(iproc)) inWsi=2
    endif
    ipWi=1  
    sign_i=1
    sign_j=1
    ! if rev parallel -> +j => -j
    if (is_rev(1,1)) then
       sign_j=-1
       ipWi=ny
    endif
    ! if rev  normal  -> +i => -i
    if (is_rev(1,2)) then
       sign_i=-1
       inWi=nghirs
       inWsi=nghirs+1
       inWri=0
    endif
    if (is_swapij(1)) then
       ! if swap invert role of i and j
       memjump=sign_j*nex  ! j
       stride =sign_i*sizeofreal ! i
       call MPI_TYPE_VECTOR(ny,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(nghirs,1,stride,type_base,type_edge,info)
    else
       ! default
       memjump=sign_i ! i
       stride =sign_j*nex*sizeofreal ! j
       call MPI_TYPE_VECTOR(nghirs,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ny,1,stride,type_base,type_edge,info)
    endif
    call MPI_TYPE_CREATE_HVECTOR(nz,1,stride2,type_edge,type_faceWinc,info)
    call MPI_TYPE_COMMIT(type_faceWinc,info)

    ! MPI-type construction for imax face "along y" (E: east)
    ! -------------------------------------------------------
    ! for imax/E, parallel dir is j and normal dir is i
    inEi=nx+1
    inEsi=nx+1-nghirs
    inEri=nx+1
    if (neighbor(nE).ne.MPI_PROC_NULL) then
       if (nob(neighbor(nE))/=nob(iproc)) inEsi=nx-nghirs
    endif
    ipEi=1
    sign_i=1
    sign_j=1
    ! if rev parallel -> +j => -j
    if (is_rev(2,1)) then
       sign_j=-1
       ipEi=ny
    endif
    ! if rev  normal  -> +i => -i
    if (is_rev(2,2)) then
       sign_i=-1
       inEi=nx+nghirs
       inEsi=nx-1
       inEri=nx+nghirs
    endif
    if (is_swapij(2)) then
       ! if swap invert role of i and j
       memjump=sign_j*nex  ! j
       stride =sign_i*sizeofreal ! i
       call MPI_TYPE_VECTOR(ny,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(nghirs,1,stride,type_base,type_edge,info)
    else
       ! default
       memjump=sign_i ! i
       stride =sign_j*nex*sizeofreal ! j
       call MPI_TYPE_VECTOR(nghirs,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ny,1,stride,type_base,type_edge,info)
    endif
    call MPI_TYPE_CREATE_HVECTOR(nz,1,stride2,type_edge,type_faceEinc,info)
    call MPI_TYPE_COMMIT(type_faceEinc,info)

    ! MPI-type construction for jmin face "along x" (S: south)
    ! --------------------------------------------------------
    ! for jmin/S, parallel dir is i and normal dir is j
    ipSi=1
    inSi=1
    inSsi=1
    inSri=1-nghirs
    if (neighbor(nS).ne.MPI_PROC_NULL) then
       if (nob(neighbor(nS))/=nob(iproc)) inSsi=2
    endif
    sign_i=1
    sign_j=1
    ! if rev parallel -> +i => -i
    if (is_rev(3,1)) then
       sign_i=-1
       ipSi=nx
    endif
    ! if rev  normal  -> +j => -j
    if (is_rev(3,2)) then
       sign_j=-1
       inSi=nghirs
       inSsi=nghirs+1
       inSri=0
    endif
    if (is_swapij(3)) then
       ! if swap invert role of i and j
       memjump=sign_j*nex  ! j
       stride =sign_i*sizeofreal ! i
       call MPI_TYPE_VECTOR(nghirs,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(nx,1,stride,type_base,type_edge,info)
    else
       ! default
       memjump=sign_i ! i
       stride =sign_j*nex*sizeofreal ! j
       call MPI_TYPE_VECTOR(nx,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(nghirs,1,stride,type_base,type_edge,info)
    endif
    call MPI_TYPE_CREATE_HVECTOR(nz,1,stride2,type_edge,type_faceSinc,info)
    call MPI_TYPE_COMMIT(type_faceSinc,info)

    ! MPI-type construction for jmax face "along x" (N: north)
    ! --------------------------------------------------------
    ! for jmax/N, parallel dir is i and normal dir is j
    ipNi=1
    inNi=ny+1
    inNsi=ny+1-nghirs
    inNri=ny+1
    if (neighbor(nN).ne.MPI_PROC_NULL) then
       if (nob(neighbor(nN))/=nob(iproc)) inNsi=ny-nghirs
    endif
    sign_i=1
    sign_j=1
    ! if rev parallel -> +i => -i
    if (is_rev(4,1)) then
       sign_i=-1
       ipNi=nx
    endif
    ! if rev  normal  -> +j => -j
    if (is_rev(4,2)) then
       sign_j=-1
       inNi=ny+nghirs
       inNsi=ny-1
       inNri=ny+nghirs
    endif
    if (is_swapij(4)) then
       ! if swap invert role of i and j
       memjump=sign_j*nex  ! j
       stride =sign_i*sizeofreal ! i
       call MPI_TYPE_VECTOR(nghirs,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(nx,1,stride,type_base,type_edge,info)
    else
       ! default
       memjump=sign_i ! i
       stride =sign_j*nex*sizeofreal ! j
       call MPI_TYPE_VECTOR(nx,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(nghirs,1,stride,type_base,type_edge,info)
    endif
    call MPI_TYPE_CREATE_HVECTOR(nz,1,stride2,type_edge,type_faceNinc,info)
    call MPI_TYPE_COMMIT(type_faceNinc,info)

    stride=nex*sizeofreal
    call MPI_TYPE_VECTOR(nx,1,1,MPI_DOUBLE_PRECISION,type_base,info)

    ! MPI-type construction for kmin face "along z" (F: forward)
    ! --------------------------------------------------------
    call MPI_TYPE_CREATE_HVECTOR( ny,1, stride,type_base,type_edge,info)
    call MPI_TYPE_CREATE_HVECTOR(nghirs,1,stride2,type_edge,type_faceFinc,info)
    call MPI_TYPE_COMMIT(type_faceFinc,info)
    
    ! MPI-type construction for kmax face "along z" (B: backward)
    ! -----------------------------------------------------------
    call MPI_TYPE_CREATE_HVECTOR( ny,1, stride,type_base,type_edge,info)
    call MPI_TYPE_CREATE_HVECTOR(nghirs,1,stride2,type_edge,type_faceBinc,info)
    call MPI_TYPE_COMMIT(type_faceBinc,info)

  end subroutine mpi_types_comm_inc_old

!!$  !===============================================================
!!$  subroutine mpi_types_comm_edges_old
!!$  !===============================================================
!!$    !> Definition of MPI types for communications ** OLD VERSION **
!!$  !===============================================================
!!$    use mod_grid
!!$    use mod_grid_directions
!!$    implicit none
!!$    ! ------------------------------------------------------------
!!$    integer :: nex,ney ! extended sizes
!!$    integer :: sign_i,sign_j
!!$    integer :: memjump ! jump in memory between block starts
!!$    ! strides in memory between block starts
!!$    integer(kind=MPI_ADDRESS_KIND) :: stride,stride2
!!$    integer :: type_base,type_edge
!!$    ! ------------------------------------------------------------
!!$
!!$    ! MPI types for faces+edges (ngh ghost cells)
!!$    ! =========================
!!$
!!$    ! Extended sizes (+ ghost cells)
!!$    ! ------------------------------
!!$    nex=nx+2*ngh
!!$    ney=ny+2*ngh
!!$
!!$    ! stride between edges along x or y
!!$    ! ---------------------------------
!!$    stride2=nex*ney*sizeofreal
!!$
!!$    !!!
!!$    ! REMARK ! changes in ipW_e but no changes in inW_e
!!$    !!!
!!$    
!!$    ! MPI-type construction for imin face "along y" (W: west)
!!$    ! -------------------------------------------------------
!!$    ! for imin/W, parallel dir is j and normal dir is i
!!$    inW_e=1
!!$    ipW_e=1-ngh 
!!$    sign_i=1
!!$    sign_j=1
!!$    ! if rev parallel -> +j => -j
!!$    if (is_rev(1,1)) then
!!$       sign_j=-1
!!$       ipW_e=ny+ngh
!!$    endif
!!$    ! if rev  normal  -> +i => -i
!!$    if (is_rev(1,2)) then
!!$       sign_i=-1
!!$       inW_e=ngh
!!$    endif
!!$    if (is_swapij(1)) then
!!$       ! if swap invert role of i and j
!!$       memjump=sign_j*nex  ! j
!!$       stride =sign_i*sizeofreal ! i
!!$       call MPI_TYPE_VECTOR(ny+2*ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
!!$       call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride,type_base,type_edge,info)
!!$    else
!!$       ! default
!!$       memjump=sign_i ! i
!!$       stride =sign_j*nex*sizeofreal ! j
!!$       call MPI_TYPE_VECTOR(ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
!!$       call MPI_TYPE_CREATE_HVECTOR(ny+2*ngh,1,stride,type_base,type_edge,info)
!!$    endif
!!$    call MPI_TYPE_CREATE_HVECTOR(nz,1,stride2,type_edge,type_e_faceW,info)
!!$    call MPI_TYPE_COMMIT(type_e_faceW,info)
!!$
!!$    ! MPI-type construction for imax face "along y" (E: east)
!!$    ! -------------------------------------------------------
!!$    ! for imax/E, parallel dir is j and normal dir is i
!!$    inE_e=nx+1
!!$    ipE_e=1-ngh
!!$    sign_i=1
!!$    sign_j=1
!!$    ! if rev parallel -> +j => -j
!!$    if (is_rev(2,1)) then
!!$       sign_j=-1
!!$       ipE_e=ny+ngh
!!$    endif
!!$    ! if rev  normal  -> +i => -i
!!$    if (is_rev(2,2)) then
!!$       sign_i=-1
!!$       inE_e=nx+ngh
!!$    endif
!!$    if (is_swapij(2)) then
!!$       ! if swap invert role of i and j
!!$       memjump=sign_j*nex  ! j
!!$       stride =sign_i*sizeofreal ! i
!!$       call MPI_TYPE_VECTOR(ny+2*ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
!!$       call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride,type_base,type_edge,info)
!!$    else
!!$       ! default
!!$       memjump=sign_i ! i
!!$       stride =sign_j*nex*sizeofreal ! j
!!$       call MPI_TYPE_VECTOR(ngh+2*ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
!!$       call MPI_TYPE_CREATE_HVECTOR(ny,1,stride,type_base,type_edge,info)
!!$    endif
!!$    call MPI_TYPE_CREATE_HVECTOR(nz,1,stride2,type_edge,type_e_faceE,info)
!!$    call MPI_TYPE_COMMIT(type_e_faceE,info)
!!$
!!$    ! MPI-type construction for jmin face "along x" (S: south)
!!$    ! --------------------------------------------------------
!!$    ! for jmin/S, parallel dir is i and normal dir is j
!!$    ipS_e=1-ngh
!!$    inS_e=1
!!$    sign_i=1
!!$    sign_j=1
!!$    ! if rev parallel -> +i => -i
!!$    if (is_rev(3,1)) then
!!$       sign_i=-1
!!$       ipS_e=nx+ngh
!!$    endif
!!$    ! if rev  normal  -> +j => -j
!!$    if (is_rev(3,2)) then
!!$       sign_j=-1
!!$       inS_e=ngh
!!$    endif
!!$    if (is_swapij(3)) then
!!$       ! if swap invert role of i and j
!!$       memjump=sign_j*nex  ! j
!!$       stride =sign_i*sizeofreal ! i
!!$       call MPI_TYPE_VECTOR(ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
!!$       call MPI_TYPE_CREATE_HVECTOR(nx+2*ngh,1,stride,type_base,type_edge,info)
!!$    else
!!$       ! default
!!$       memjump=sign_i ! i
!!$       stride =sign_j*nex*sizeofreal ! j
!!$       call MPI_TYPE_VECTOR(nx+2*ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
!!$       call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride,type_base,type_edge,info)
!!$    endif
!!$    call MPI_TYPE_CREATE_HVECTOR(nz,1,stride2,type_edge,type_e_faceS,info)
!!$    call MPI_TYPE_COMMIT(type_e_faceS,info)
!!$
!!$    ! MPI-type construction for jmax face "along x" (N: north)
!!$    ! --------------------------------------------------------
!!$    ! for jmax/N, parallel dir is i and normal dir is j
!!$    ipN_e=1-ngh
!!$    inN_e=ny+1
!!$    sign_i=1
!!$    sign_j=1
!!$    ! if rev parallel -> +i => -i
!!$    if (is_rev(4,1)) then
!!$       sign_i=-1
!!$       ipN_e=nx+ngh
!!$    endif
!!$    ! if rev  normal  -> +j => -j
!!$    if (is_rev(4,2)) then
!!$       sign_j=-1
!!$       inN_e=ny+ngh
!!$    endif
!!$    if (is_swapij(4)) then
!!$       ! if swap invert role of i and j
!!$       memjump=sign_j*nex  ! j
!!$       stride =sign_i*sizeofreal ! i
!!$       call MPI_TYPE_VECTOR(ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
!!$       call MPI_TYPE_CREATE_HVECTOR(nx+2*ngh,1,stride,type_base,type_edge,info)
!!$    else
!!$       ! default
!!$       memjump=sign_i ! i
!!$       stride =sign_j*nex*sizeofreal ! j
!!$       call MPI_TYPE_VECTOR(nx+2*ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
!!$       call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride,type_base,type_edge,info)
!!$    endif
!!$    call MPI_TYPE_CREATE_HVECTOR(nz,1,stride2,type_edge,type_e_faceN,info)
!!$    call MPI_TYPE_COMMIT(type_e_faceN,info)
!!$
!!$    stride=nex*sizeofreal
!!$    call MPI_TYPE_VECTOR(nx+2*ngh,1,1,MPI_DOUBLE_PRECISION,type_base,info)
!!$    
!!$    ! MPI-type construction for kmin face "along z" (F: forward)
!!$    ! --------------------------------------------------------    
!!$    call MPI_TYPE_CREATE_HVECTOR(ny+2*ngh,1,stride,type_base,type_edge,info)
!!$    call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride2,type_edge,type_e_faceF,info)
!!$    call MPI_TYPE_COMMIT(type_e_faceF,info)
!!$    
!!$    ! MPI-type construction for kmax face "along z" (B: backward)
!!$    ! -----------------------------------------------------------
!!$    call MPI_TYPE_CREATE_HVECTOR(ny+2*ngh,1,stride,type_base,type_edge,info)
!!$    call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride2,type_edge,type_e_faceB,info)
!!$    call MPI_TYPE_COMMIT(type_e_faceB,info)
!!$
!!$  end subroutine mpi_types_comm_edges_old
  
  !===============================================================
  subroutine free_mpi_types_comm_old
  !===============================================================
    !> Free MPI types for communications ** OLD VERSION **
  !===============================================================
    implicit none
    ! ------------------------------------------------------------
    ! ------------------------------------------------------------
    
    ! MPI types for communications of conservative variables
    ! [inviscid fluxes]
    call MPI_TYPE_FREE(type_faceW,info)
    call MPI_TYPE_FREE(type_faceE,info)
    call MPI_TYPE_FREE(type_faceS,info)
    call MPI_TYPE_FREE(type_faceN,info)
    call MPI_TYPE_FREE(type_faceF,info)
    call MPI_TYPE_FREE(type_faceB,info)
    
    ! MPI types for communications of velocity derivatives
    ! [viscous fluxes]
    call MPI_TYPE_FREE(type_faceW_v,info)
    call MPI_TYPE_FREE(type_faceE_v,info)
    call MPI_TYPE_FREE(type_faceS_v,info)
    call MPI_TYPE_FREE(type_faceN_v,info)
    call MPI_TYPE_FREE(type_faceF_v,info)
    call MPI_TYPE_FREE(type_faceB_v,info)
        
  end subroutine free_mpi_types_comm_old
  
!!$  !===============================================================
!!$  subroutine free_mpi_types_comm_edges_old
!!$  !===============================================================
!!$    !> Free MPI types for communications ** OLD VERSION **
!!$  !===============================================================
!!$    implicit none
!!$    ! ------------------------------------------------------------
!!$    ! ------------------------------------------------------------
!!$    
!!$    ! MPI types for communications of faces+edges
!!$    call MPI_TYPE_FREE(type_e_faceW,info)
!!$    call MPI_TYPE_FREE(type_e_faceE,info)
!!$    call MPI_TYPE_FREE(type_e_faceS,info)
!!$    call MPI_TYPE_FREE(type_e_faceN,info)
!!$    call MPI_TYPE_FREE(type_e_faceF,info)
!!$    call MPI_TYPE_FREE(type_e_faceB,info)
!!$            
!!$  end subroutine free_mpi_types_comm_edges_old
  
  !===============================================================
  subroutine free_mpi_types_metrics_old
  !===============================================================
    !> Free MPI types for communications ** OLD VERSION **
  !===============================================================
    implicit none
    ! ------------------------------------------------------------
    ! ------------------------------------------------------------

    ! MPI types for communications of 2D curvilinear metrics
    ! [inviscid fluxes]
    call MPI_TYPE_FREE(type_mW,info)
    call MPI_TYPE_FREE(type_mE,info)
    call MPI_TYPE_FREE(type_mS,info)
    call MPI_TYPE_FREE(type_mN,info)
    
    ! MPI types for communications of 2D curvilinear metrics
    ! [viscous fluxes]
    call MPI_TYPE_FREE(type_mW_v,info)
    call MPI_TYPE_FREE(type_mE_v,info)
    call MPI_TYPE_FREE(type_mS_v,info)
    call MPI_TYPE_FREE(type_mN_v,info)
        
  end subroutine free_mpi_types_metrics_old

end module mod_mpi_types_two_sided_old
