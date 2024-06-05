!===============================================================
! Programme
! AUTHOR : Xavier Gloerfelt
!          Laboratoire DynFluid - Arts et Metiers
!===============================================================
program MUSICAA
  use mod_mpi
  use mod_block ! <- nbloc for field ?
  use mod_flow
  use mod_eos
  use mod_constant
  use mod_comm
  use mod_time
  use mod_interface
  use mod_io
  use mod_io_snapshots
  use mod_tecplot
  use mod_tranprop
  use mod_filtering
  use mod_artvisc
  use mod_coeff_deriv
  use mod_init_flow ! <- for is_pulse, x_pulse TO BE CHANGED
  use mod_routines  ! <- for is_src            TO BE CHANGED
  use mod_saturation_curve
  use mod_grid_utilities
  ! use mod_add_gh3d  ! Not necessary anymore
  use mod_pp_main
  use warnstop
  use mod_comm1     ! for: communication1 (one-sided RMA)
  use mod_mpi_types_two_sided
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: ndim ! <- was in mod_grid but used only for init_io
  ! ---------------------------------------------------------------------------

  ! First initialisations
  ! =====================
  ! MPI initialisation
  call init_mpi
  ! initialize CPU time counters
  call init_time
  ! print welcome message
  if (iproc.eq.0) call welcome

  ! Read of parameter files
  !========================
  call read_param('param.ini')
  if (idepart.ne.POST_PROCESSING) dirDATA='./'
  call read_param_blocks(trim(dirDATA)//'param_blocks.ini')

  ! Grid pre-processing
  ! =====================
  if (idepart==PRE_PROCESSING) then
     if ((is_add_sponge.or.is_coarse_grid.or.is_half_cell.or.is_satur_curve).and.&
         (nproc.ne.1)) call mpistop("Pre-processing mode selected needs to be run on 1 processor",0)

     ! Saturation curve mode
     ! ---------------------
     if (is_satur_curve) call mpistop("Saturation curve mode needs to be properly implemented back",0)

     ! Linear stability solver
     ! -----------------------
     if (is_LST) call mpistop("Linear Stability solver mode needs to be properly implemented back",0)

     ! Streching of grid at the exit
     ! -----------------------------
     if (is_add_sponge) then
        !call add_sponge('Grids_for_python_x','cyl_mod')
        call mpistop("Subroutine add_sponge needs to be properly implemented back",0)
     endif

     ! Creation of coarse grid
     ! -----------------------
     if (is_coarse_grid) then
        if (is_curv3) then
           call coarse_grid3d
        else
           call coarse_grid
        endif
        call mpistop('Generation of coarse grid completed... Update param_blocks.ini',0)
     endif

     ! Half-cell suppression (if not is_adjoint_blocks)
     ! ---------------------
     ! If is_coarse_grid, half-cell suppression is realized on coarse grid
     if (is_half_cell) then
        if (is_curv3) then
           call modif_grid3d
        else
           ! print *,'dirGRID:',trim(dirGRID),' nameGRID:',trim(nameGRID)
           call modif_grid
        endif
        call mpistop('Generation of half-cell grid completed...',0)
     endif


     !call compute_clapeyron

     call mpistop('End of pre-processing',0)
  endif

  ! Post-processing
  ! ===============
  if (idepart==POST_PROCESSING) call read_param_pp

  ! Init equation of state & transport properties
  ! =============================================
  call init_eos
  call init_viscosity

!!$  call pre_thermo
!!$  call mpistop('',0)

  ! MPI partitioning: set dimensions
  ! ================
  call mpi_set_dim
  call MPI_BARRIER(COMM_global,info)
  
  ! Setup of reference quantities
  ! =============================
  if (CHIT.or.CHAN.or.PHILL.or.STBL.or.CYL.or.SHIT.or.TURB.or.SRC.or.LE.or.TE.or.T3C) call setupref

  ! MPI partitioning: intrablock communicator
  ! ================
  call mpi_intrablock

  ! is_init_2D3D must be determined for grid_define
  if (idepart.eq.FROM_FILE) call read_write_info(READ_DEB)

  ! Define grid
  ! ===========
  call grid_define

  ! MPI partitioning: connectivity
  ! ================
  call mpi_connect
  call MPI_BARRIER(COMM_global,info)

  ! Add ghost cells <~ TO BE CHANGED
  ! ===============
  !!call add_ghost_cells_grid3d!('Grid5/','sphere5c')
  !!call add_ghost_cells_grid3d!('Grid5/','sphere5')
  !!call add_ghost_cells_grid3d!('Grid3/','sphere3')
  !!call add_ghost_cells_grid3d!('Grid_50x50x100/','sphere2')
  !!call add_ghost_cells_grid3d!('./','cdnoz3D2')
  !!call add_ghost_cells_grid3d!('./','ann_sect_mod')
  !!call add_ghost_cells_grid3d!('Grid_stator/','stator_mod')
  !!call add_ghost_cells_grid3d!('Grid_RANS/','turb_mod')
  !!call add_ghost_cells_grid3d!('Grid_RANS_GO/','ls593_mod')
  ! if ((is_curv3).and.(idepart==1)) call add_ghost_cells_grid3d
  !call mpistop('',0)

  ! Scaling of grid
  ! ===============
  if (.not.is_grid_old) call grid_extend

  ! Definition of MPI types for communications
  ! ==========================================
  !call mpi_types_comm_old
  !call mpi_types_comm_v_old
  call mpi_types_comm(ngh,'r')
  call mpi_types_comm(ngh_v,'v')
  call mpi_types_comm_ex(ngh)
  if (is_comm_onesided) call mpi_types_comm1

  call MPI_BARRIER(COMM_global,info)

  ! Define boundary conditions
  ! ==========================
  call bc_define
  
  ! Define MPI index bounds
  ! =======================
  call mpi_index_bounds

!!$  if (iproc==0) then
!!$     print *,'ndx_v ', ndx_v,nfx_v,ndy_v,nfy_v,ndz_v,nfz_v
!!$     print *,'ndx_vi',ndx_vi,nfx_vi,ndy_vi,nfy_vi,ndz_vi,nfz_vi
!!$     print *,'ndxt_v',ndxt_v,nfxt_v,ndyt_v,nfyt_v,ndzt_v,nfzt_v
!!$     print *,'ndx_v1',ndx_v1,nfx_v1,ndy_v1,nfy_v1,ndz_v1,nfz_v1
!!$     print *,'ndx_v2',ndx_v2,nfx_v2,ndy_v2,nfy_v2,ndz_v2,nfz_v2
!!$     print *,'ndx_v3',ndx_v3,nfx_v3,ndy_v3,nfy_v3,ndz_v3,nfz_v3
!!$  endif
  
!!$  if ((ndx_vi.ne.ndx_v3).or.(nfx_vi.ne.nfx_v3).or.(ndy_vi.ne.ndy_v3).or.(nfy_vi.ne.nfy_v3)) then
!!$     print *,'iproc',iproc,ndx_vi,nfx_vi,ndy_vi,nfy_vi
!!$     print *,'iproc',iproc,ndx_v3,nfx_v3,ndy_v3,nfy_v3
!!$  endif
!!$  
!!$  call mpistop('check bounds!', 0)

  ! Initializations of scheme coefficients
  ! ======================================
  ! Finite-differences coefficients
  ! -------------------------------
  call init_coeff_deriv(stencil,is_DRP)

  ! Summation-By-Parts (SBP) schemes
  ! --------------------------------
  call init_coeff_SBP(6)
  ! TEST SBP
  !call test_coeff_SBP(6)
  !call mpistop('stop SBP!', 0)

  ! Extrapolation coefficients for bc_wall in cartesian
  ! ---------------------------------------------------
  if ((.not.is_curv3).and.(.not.is_curv)) call init_coeff_wall_cart

  ! Runge-Kutta coefficients
  ! ------------------------
  call init_RK

  ! Init IO
  ! =======
  allocate(field(nbloc))
  field(nob(iproc))%MPI_COMM = COMM_intrablock
  ndim=3 ! not useful in modules !!!!!!!!!!!!!
  call mod_io_init(ngx,ngy,ngz,nx,ny,nz,nx,ny,nz,ngh,ndim,coord &
       ,is_IOtec_read,is_IOtec_write,field(nob(iproc)))
  if (iproc.eq.0) print *,"Freq. vol., plane, lines & points:",freq_volume,freq_plane,freq_line,freq_point
  call init_io_snapshots
  if (is_TamDong) call init_io_restartTD(restartTD,0)
  
  ! Mesh construction
  ! =================
  ! if (idepart==POST_PROCESSING) then
  !    call grid_pp
  ! else
     call grid_finalize
  ! endif

  ! if ((STBL.or.SHIT).and.(.not.is_curv).and.(nbloc.eq.1)) tscale = (xg(ngx))/u_ref
  if ((STBL.or.SHIT).and.(.not.is_curv).and.(nbloc.eq.1)) tscale = abs(xmax-xmin)/u_ref

  ! Setup of reference quantities
  ! =============================
  if (.not.(CHIT.or.CHAN.or.PHILL.or.CYL.or.TURB)) call setupref

  ! Assignment of procedure pointers
  ! ================================
  call assign_procedures

  ! Allocation of flow arrays
  ! =========================
  call alloc_tab

  ! Define windows on allocated arrays
  ! ==================================
  if (is_comm_onesided) call mpi_win_comm1
  
  call MPI_BARRIER(COMM_global,info)

  ! Call solver routines
  ! ====================
  call solver

  ! Closing sequence
  ! ================
  if (iproc.eq.0) write(*,*) 'End of computation..'

  t_end = MPI_WTIME()
  call MPI_BARRIER(COMM_global, info)
  if (iproc.eq.0) write(*,'(A,i4,A,g0)') 'proc', iproc, ' t_end-t_start ', t_end-t_start

  if (is_comm_onesided) call mpi_close_win_comm1
  !if ((is_irs).and.(is_comm_onesided)) call mpi_close_win_comm1_inc

  ! /!\ TO BE UPDATED
  call lib_tab
  call MPI_COMM_FREE(COMM_global, info)

  call mpistop('', 0)

end program MUSICAA
