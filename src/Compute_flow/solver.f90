!===============================================================================
subroutine solver
!===============================================================================
  !> call solver
!===============================================================================
  use mod_time
  use mod_eos
  use mod_init_flow
  use mod_init_hit
  use mod_forcing_bulk
  use mod_eigenmode
  use mod_sponge
  use mod_interface
  use mod_bc_apply
  use mod_init_TamDong
  use mod_bc_inlet_outlet
  use mod_io
  use mod_io_restart_BCref
  use mod_io_snapshots
  use mod_io_stats
  use mod_utils
  use mod_filtering
  use mod_filtering_shock
  use mod_artvisc
  use mod_artvisc_shock
  use mod_filtering_inc
  use mod_artvisc_inc
  use mod_artvisc_shock_inc
  use mod_init_2D_3D
  use mod_comm ! for communication_2d_edges not yet proc pointer
  use mod_flow0
  use mod_pp_main
  use mod_analytical_sol
  use mod_init_irs
  use mod_rfm ! for init_vel_rfm
  use mod_rans
  use mod_wall_dist
  use mod_artvisc_rans
  use mod_artvisc_shock_rans
  use mod_turb_model_length_scale
  use mod_wall_model

  implicit none
  ! ----------------------------------------------------------------------------
  logical :: iexist
  integer :: isn
  !character(len=60) :: dirDATA
  real(wp) :: Tround
  ! ----------------------------------------------------------------------------
  !integer :: k

  select case (idepart)

  case (FROM_SCRATCH)
     
     if (iproc.eq.0) then
        inquire( file='time.ini', exist=iexist )
        if (iexist) then
           call system('mv time.ini time_bak.ini')
        endif
     endif
     ntotal = 0
     ntotal_old = 0
     time  = 0.0_wp
     tstar = 0.0_wp
     
     if (CHIT) then
        if (trim(filestamp).eq.'0000_0000') then
           binfile = 'restart_bl'//trim(numchar(nob(iproc)))//filext_read
           TDfile = 'restartTD.bin'
        else
           binfile = 'restart'//filestamp//'_bl'//trim(numchar(nob(iproc)))//filext_write
           TDfile = 'restartTD'//filestamp//'.bin'
        endif
        call read_write_info(WRITE)
        call init_thermo
        call init_hit
     else
        call init_flow
     endif
     call stats_init

  case (FROM_FILE)
     
     call init_thermo

     if (trim(filestamp).eq.'0000_0000') then
        binfile= 'restart_bl'//trim(numchar(nob(iproc)))//filext_read
        TDfile = 'restartTD.bin'
     else
        binfile= 'restart'//filestamp//'_bl'//trim(numchar(nob(iproc)))//filext_read
        TDfile = 'restartTD'//filestamp//'.bin'
     endif

     call read_write_info(READ)

     if (is_init_2D3D) then
        call extrude_2D_field
     else
        call read_write_volume(binfile,TDfile,READ)
     endif
     
     ! Read stats quantities if ntotal> ndeb
     ! =====================================
     if (is_curv3) then
        if (TURB) call read_write_stats_xy_xyz(READ,nob(iproc))
        ! for starting new stats in LS59 rough [old are 2D only]
        !if (TURB) call read_write_stats_xy(READ,nob(iproc))
        if (CYL) call read_write_stats_xyz(READ)
        !if (CYL.or.TURB) call read_write_stats_xyz(READ)
     else
        if (STBL.or.CYL.or.SHIT.or.ACT.or.TURB.or.LE.or.TE) call read_write_stats_xy(READ,nob(iproc))
        if (CHAN) call read_write_stats_chan(READ)
     endif

     ! Init stats.dat (HIT stats or forces) TO BE CHANGED
     ! ====================================
     if (CHIT.or.CYL.or.ACT.or.TURB) call stats_init

  case (FROM_INTERP)
     
     call init_thermo
     
     if (trim(filestamp).eq.'0000_0000') then
        binfile = 'restart_bl'//trim(numchar(nob(iproc)))//filext_read
        TDfile = 'restartTD.bin'
     else
        binfile = 'restart'//filestamp//'_bl'//trim(numchar(nob(iproc)))//filext_read
        TDfile = 'restartTD'//filestamp//'.bin'
     endif

     call init_interp

     call read_write_volume('restart_bl'//trim(numchar(nob(iproc)))//filext_write, &
                            'restartTD.bin',WRITE)
     
     ntotal = 0
     ntotal_old = 0
     time  = 0.0_wp
     tstar = 0.0_wp
     
     call read_write_info(WRITE)

     call primitives_visc!(rho,rhou,rhov,rhow,rhoe)

     call mpistop('TEST INTERP (lavoro in corso)',0)

     call stats_init
    
  case (POST_PROCESSING)

     call read_write_info(READ)

     call pp_main

     call mpistop('end Postprocessing', 0)

  case default
     call mpistop('error in depart choice!', 0)
  end select

  call MPI_BARRIER(COMM_global,info)

  ! Initialization BC reference values
  ! ==================================
  if (is_BC_ref) call init_BC_ref

  call MPI_BARRIER(COMM_global,info)

  ! Initialization Random Fourier Modes
  ! ===================================
  if (is_RFM) call init_RFM

  ! Initialization of wall-model
  ! ============================
  if (is_wall_model) call init_wm

  ! Fill ghost cells
  ! ================
  call communication_(rho,rhou,rhov,rhow,rhoe)
  ! apply angular periodicity if necessary
  call bc_angular_periodicity
  if (is_mean0) call communication_(rho0,u0,v0,w0,p0)
  if (is_RANS) call communication_rans(nutil)

  call MPI_BARRIER(COMM_global,info)

  ! Compute primitives variables (with ghost cells)
  ! ============================
  call primitives_visc!(rho,rhou,rhov,rhow,rhoe)
  
  ! Init Tam & Dong's BC
  ! ====================
  if (is_TamDong) call init_bc_TD
  if (is_inlet_outlet) call init_bc_inlet_outlet

  ! Init conservative variables arrays
  ! ==================================
  rho_n  = rho
  rhou_n = rhou
  rhov_n = rhov
  rhow_n = rhow
  rhoe_n = rhoe
  
  ! Time step computation (old value read in info.ini)
  ! =====================
  if ((ntotal<ndeb).and.(CFL>0)) then
     !call timestep
     call timestep_new
  else
     if ((is_dtlocal).and.(.not.allocated(dt_local))) allocate(dt_local(nx,ny,nz))
  endif

  !call mpistop('check',0)

  ! IRS initialization
  ! ==================
  if (is_irs) call init_irs

  ! Check max CFL
  ! =============
  call check_max_cfl
  !!call check_max_cfl_interface
  !!call mpistop('on arrete',0)
  
!!$  if (coord(1)==0) then
!!$     print *,is_BC_ref,BC_face(1,1)%is_mean_ref
!!$     do k=1,nz
!!$        do j=1,ny
!!$           rho(1,j,k)=BC_face(1,1)%Uref(j,k,1)
!!$           uu(1,j,k)=BC_face(1,1)%Uref(j,k,2)
!!$           vv(1,j,k)=BC_face(1,1)%Uref(j,k,3)
!!$           ww(1,j,k)=BC_face(1,1)%Uref(j,k,4)
!!$           rhou(1,j,k)=rho(1,j,k)*uu(1,j,k)
!!$           rhov(1,j,k)=rho(1,j,k)*vv(1,j,k)
!!$           rhow(1,j,k)=rho(1,j,k)*ww(1,j,k)
!!$           prs(1,j,k)=BC_face(1,1)%Uref(j,k,5)
!!$           rhoe(1,j,k)= prs(1,j,k)/gam1+0.5_wp*rho(1,j,k)*(uu(1,j,k)**2+vv(1,j,k)**2+ww(1,j,k)**2)
!!$        enddo
!!$     enddo
!!$  endif

  ! Write initial field OR NOT (TO BE CHANGED)
  ! ===================
  if (idepart.eq.FROM_SCRATCH) then
     ! call read_write_volume('restart'//filestamp//'_bl'//trim(numchar(nob(iproc)))//filext_write, &
     !                       'restartTD'//filestamp//'.bin',WRITE)
!!$     call read_write_volume('restart_bl'//trim(numchar(nob(iproc)))//filext_write, &
!!$                            'restartTD.bin',WRITE)
     call read_write_info(WRITE)
     if (CHIT) call stats_hit
  endif

!!$  ! check for angular periodicity
!!$  open(194,file='sol_bl'//trim(numchar(nob(iproc)))//'_ex.bin',form='unformatted',status='unknown')
!!$  rewind(194)
!!$  write(194) ngy
!!$  write(194) ngz+2*ngh
!!$  write(194) ((ygc3(10,j,k),j=1,ngy),k=1-ngh,ngz+ngh)
!!$  write(194) ((zgc3(10,j,k),j=1,ngy),k=1-ngh,ngz+ngh)
!!$  write(194) ((rhov(10,j,k),j=1,ngy),k=1-ngh,ngz+ngh)
!!$  write(194) ((rhow(10,j,k),j=1,ngy),k=1-ngh,ngz+ngh)
!!$  close(194)

!!$  call mpistop('after init in solver',0)
    
  ! Initializations for forcing mass flow rate
  ! ==========================================
  if (is_forcing_bulk) call init_forcing_bulk
  
  ! Initializations for entering eigenmodes
  ! =======================================
  if (is_eigenmode) then
     if (coord(1)==0) call init_eigenmodes
  endif

  if (STBL.and.is_eigenmode) then
     if (coord(1)==0) then
        if (iproc.eq.0) write(6,*) '  T/dt            =', 2.0_wp*pi/eig(1)%om/deltat
        Tround=anint(2.0_wp*pi/eig(1)%om/deltat/10.0_wp)*10.0_wp
        deltat=2.0_wp*pi/eig(1)%om/Tround
        dtstar = deltat/tscale
        if (iproc.eq.0) write(6,*) 'round(T/dt)       =', Tround
        if (iproc.eq.0) write(6,*) 'new Dt [s], Dt* [-]:', deltat, dtstar
     endif
     call MPI_BCAST(deltat,1,MPI_DOUBLE_PRECISION,0,COMM_global,info)
     call MPI_BCAST(dtstar,1,MPI_DOUBLE_PRECISION,0,COMM_global,info)
  endif
  
  ! Initializations of numerical dissipation coefficients
  ! =====================================================
  if (is_dissip_in_increments) then
     if (is_shock) then ! Artificial viscosity (DNC-Jameson)
        call init_artvisc_shock_inc(dissip_coeff,dissip_shock)
     else  
        if (is_SF) then ! Selective filtering
           call init_filtering_inc(stencil,is_DRP,dissip_coeff)
        else ! Artificial viscosity (DNC)
           call init_artvisc_inc(dissip_coeff)
        endif
     endif
  else
     if (is_shock) then ! Artificial viscosity (DNC-Jameson)
        if (is_SF) then ! Selective filtering
           call init_filter_shock(dissip_coeff,dissip_shock)
        else ! Artificial viscosity (DNC)
           call init_artvisc_shock(dissip_coeff,dissip_shock)
        endif
        if (is_RANS) call init_artvisc_shock_rans(dissip_coeff,dissip_shock)
     else  
        if (is_SF) then ! Selective filtering
           call init_filtering(stencil,is_DRP,dissip_coeff)
        else ! Artificial viscosity (DNC)
           call init_artvisc(dissip_coeff)
        endif
        if (is_RANS) call init_artvisc_rans(dissip_coeff)
     endif
  endif

  ! Initializations for sponge zones
  ! ================================
  if (is_curv3) then
     call init_sponge_3d
  else
     call init_sponge
  endif

  ! Initializations for suction and blowing
  ! =======================================
  !!!if (STBL.and.is_forc_sb) call init_forcing

  ! Info TO BE CHANGED
  ! ====
  if (CHAN.or.PHILL) then
     if (iproc.eq.0) write(*,*) 'Iterations/Turn :', (xmax-xmin)/u_ref/deltat

  elseif (TGV) then
     ndeb = max(nmax/200, 1)
     if (iproc.eq.0) write(*,*) 'TOTAL ITERATIONS:', nmax
     
  elseif (CHIT) then
     ndeb = max(nmax/freq_stats, 1)
     freq_stats=ndeb
     if (iproc.eq.0) write(*,*) 'TOTAL ITERATIONS:', nmax

  elseif (STBL.or.SRC.or.CAV.or.ACT.or.SHIT.or.LE.or.TE) then
     if (iproc.eq.0) write(*,*) 'TOTAL ITERATIONS:', nmax
     if (iproc.eq.0) write(*,*) 'Iterations/Turn :', (xmax-xmin)/c_ref/deltat
     
  elseif ((CYL).or.(TURB)) then
     !call write_surf_norm
     if (iproc.eq.0) write(*,*) 'TOTAL ITERATIONS:', nmax
     if (u_ref.ne.0.0_wp) then
        if (iproc.eq.0) write(*,*) 'Iterations/Turn :', 2.0_wp*(xmax-xmin)/u_ref/deltat
     endif
  else
     call mpistop('Timestep not implemented for this flowtype!!', 0)
  endif

  ! RANS modeling !! TEMP !!
  ! ------------------------
  if (is_RANS) then
     call wall_dist
     call init_rans

     ! compute and store max cell dimensions in case of DES97, DDES, and IDDES approaches
     ! in case of moving mesh in the future, call it in "runge_kutta_rans" subroutine
     if ((simulation_RANS.eq.'DES97').or.(simulation_RANS.eq.'DDES').or.(simulation_RANS.eq.'IDDES')) &
     call maximum_cell_dimension

     ! Init conservative variables arrays
     ! ==================================
     nutil_n = nutil
  endif

  is_RANS=.false.

  ! ---------------------------------------------------------------------------
  call MPI_BARRIER(COMM_global,info)
  if (iproc==0) print *,'Init ok. Starting Time Loop..'

  ! Snapshots of the field without any iteration
  if (nmax.eq.0) then
     do isn=1,nsnapshots
        call write_snapshot(isn)
     end do
  endif

  ! ===========================================================================
  ! ===========================================================================
  !                             MAIN LOOP
  ! ===========================================================================
  ! ===========================================================================
  looptime: do ntime=1, nmax

     ! recompute timestep if is_dtvar
     ! ------------------------------
     if (is_dtvar) call timestep
     if (is_dtlocal) call local_timestep

     ! start measure of elapsed CPU time
     ! ---------------------------------
     time_before=MPI_WTIME()

     ! wait for RANS to kick in
     ! ------------------------
     if (ntotal.ge.ndeb_RANS) is_RANS=.true.

     ! increment iteration counter
     ! ---------------------------
     ntotal=ntotal+1

     ! call time integration
     ! ---------------------
     if (ntime==1) call start_runge_kutta ! start first RK step
     if (is_RANS.and.ntotal.eq.ndeb_RANS) call start_runge_kutta_rans
     call runge_kutta

     ! increment physical and nondimensional time
     ! -------------------------------------------
     time  = time + deltat
     tstar = time/tscale

     ! check forcing
     ! -------------
     if ((is_forcing_bulk).and.(.not.(is_2d))) then
        if (mod(ntime,nprint).eq.0) call check_forcing
        ! if (mod(ntime,10).eq.0) call check_forcing
     endif


     !if ((is_forcing_bulk).and.(mod(ntime,10).eq.0)) call check_forcing
!!     if (mod(ntime,nprint).eq.0) then  !!!!!  TO BE CHANGED
!!        call check_max_cfl
!!        call check_max_cfl_interface
!!     endif

     ! handle outputs
     ! --------------
     !!call MPI_BARRIER(COMM_global,info)

     !call mean_local_cfl
     !call local_cfl

     printstep = .false.
     if (mod(ntotal,freq_stats).eq.0) printstep(1)=.true.
     if (mod(ntotal,freq_plane).eq.0) printstep(2)=.true.
     if (mod(ntotal,freq_field).eq.0) printstep(3)=.true.

     if (printstep(1)) call stats
     if (is_curv.and.CYL.and.printstep(1)) call compute_cl_cd
     if (is_curv3.and.CYL.and.printstep(1)) call compute_cl_cd3

     ! Writting, according to the frequency, of snapshot outputs
     call calcfilestamp(tstar,filestamp)
     do isn=1,nsnapshots
        if (mod(ntime,snapshots(isn)%freq).eq.0) call write_snapshot(isn)
     end do

     if (printstep(3)) then
        if (is_curv3) then
           if (TURB) call read_write_stats_xy_xyz(WRITE,nob(iproc))
           if (CYL) call read_write_stats_xyz(WRITE)
           !if (CYL.or.TURB) call read_write_stats_xyz(WRITE)
        else
           if (STBL.or.CYL.or.SHIT.or.ACT.or.TURB.or.LE.or.TE) call read_write_stats_xy(WRITE,nob(iproc))
           if (CHAN) call read_write_stats_chan(WRITE)
        endif
        call calcfilestamp(tstar,filestamp)
        binfile = 'restart'//filestamp//'_bl'//trim(numchar(nob(iproc)))//filext_write
        TDfile = 'restartTD'//filestamp//'.bin'
        call read_write_info(WRITE)
        call read_write_volume(binfile,TDfile,WRITE)
     endif

     if ((CHIT).and.((ntime==int(tscale*(98-42)/deltat)+1).or. &
          (ntime==int(tscale*(171-42)/deltat)+1))) call stats_hit

     ! if ((CHIT).and.(ntime==int(10*0.0003918567081900893/deltat)+1)) call stats_hit

     ! if ((CYL).and.(mod(ntotal,freq_plane).eq.0)) call compute_cl_cd

     ! end measure of elapsed CPU time
     ! ---------------------------------
     time_after = MPI_WTIME()
     
     cputime = time_after - time_before
     cputot = cputot + cputime
     cpurun = cpurun + cputime

     ! write informations at screen
     ! ----------------------------
     if (mod(ntime,nprint).eq.0) then
        call calcfilestamp( tstar, filestamp )
        if (iproc.eq.0) then
           !print *, ntotal_old,nmax,ntotal,timemax,tstar,dtstar
           write(*,'(a)') repeat('=',80)
!!$           write(*,'(1x,a,3(i0, a),10x,a,f15.5)') 'Iteration: ', ntime,   &
!!$                ' ( ', ntotal, '/', min( ntotal_old + nmax, ntotal+int((timemax-tstar)/dtstar)), ' )' &
!!$                , 'Tstar:', tstar
           write(*,'(1x,a,3(i0, a),10x,a,f15.5)') 'Iteration: ', ntime,   &
                ' ( ', ntotal, '/', ntotal_old + nmax, ' )' &
                , 'Tstar:', tstar
!!$           write(*,'(1x,a,2(i0, a),10x,a,f15.5)') 'Iteration: ', ntime,   &
!!$                ' ( ', ntotal, ' )' &
!!$                , 'Tstar:', tstar
           write(*,20) cputime, int(cpurun), nproc*cputot/3600.0_wp
        endif
     endif
20   format(10x,'CPU time/step:',f12.5,' s, Run time: ',i0,' s, Total CPU:',f10.3,' h')

     ! handle exit conditions
     ! ----------------------
     ! ~> maximum simulation time reached
     if (tstar.gt.timemax) then
        if (iproc.eq.0) write(*,*) 'Max Tstar reached = ', timemax
        exit looptime
     endif
     ! ~> maximum computational time reached
     if (cpurun.gt.cpumax) then
        if (iproc.eq.0) write(*,*) 'Max cpu time reached = ', cpumax
        exit looptime
     endif
     ! ~> exit on external input
     if (mod(ntime,1000).eq.0) then
        if (iproc.eq.0) call MPI_FILE_DELETE('./stop', MPI_INFO_NULL, icheck)
        call MPI_BCAST(icheck,1,MPI_INTEGER,0,COMM_global,info)
        if (icheck.eq.0) then
           if (iproc.eq.0) write(*,*) 'Stop file found. Exiting timeloop..'
           call MPI_BARRIER(COMM_global,info)
           exit looptime
        endif
     endif

  enddo looptime
  ! ===========================================================================
  ! ===========================================================================

  if (iproc.eq.0) write(*,*) 'End of Time Loop'

!!$  print *,'Tmp-wall',Tmp(:,1,1)

!!$  i=1
!!$  open(55,file='profil_M0p9_R500.dat')
!!$  do j=1,ny
!!$     write(55,'(5f16.8)') yg(j),u0(i,j,1),v0(i,j,1),p0(i,j,1),rho0(i,j,1)
!!$  enddo
  
!!$  open(194,file='sol_bl'//trim(numchar(nob(iproc)))//'.bin',form='unformatted',status='unknown')
!!$  rewind(194)
!!$  write(194) ((prs(i,j,1),i=-4,ngx+5),j=-4,ngy+5)
!!$  close(194)
  
  ! write output files at last step only if not printstep=.true.
  ! ------------------------------------------------------------
  if (.not.printstep(1)) call stats
  !if (.not.is_curv3.and.CYL.and.printstep(1)) call compute_cl_cd
!!$  if ((CYL.or.TURB).and.printstep(1)) call compute_cl_cd
  !if (.not.printstep(2)) then
  !   call calcfilestamp(tstar,filestamp)
  !   do ip=1,npl
  !      call write_plane(ip)
  !   enddo
  !endif
  if (.not.printstep(3)) then
     if (is_curv3) then
        if (TURB) call read_write_stats_xy_xyz(WRITE,nob(iproc))
        if (CYL) call read_write_stats_xyz(WRITE)
        !if (CYL.or.TURB) call read_write_stats_xyz(WRITE)
     else
        if (STBL.or.CYL.or.SHIT.or.ACT.or.LE.or.TURB.or.TE) call read_write_stats_xy(WRITE,nob(iproc))
        if (CHAN) call read_write_stats_chan(WRITE)
     endif
     call calcfilestamp(tstar,filestamp)
     call read_write_info(WRITE)
     call read_write_volume('restart_bl'//trim(numchar(nob(iproc)))//filext_write, &
                            'restartTD.bin',WRITE)
  else
     if (iproc.eq.iproc_leader(nob(iproc))) then
        call system('mv '//binfile//'  restart_bl'//trim(numchar(nob(iproc)))//filext_write)
        call system('mv '//TDfile//'  restartTD.bin')
     endif
  endif

  if ((is_BC_ref).and.(is_BCref_init)) call write_restart_BCref("restart_BCref.bin")

  ! compute and write analytical solution and error
  ! -----------------------------------------------
  if (SRC.and.is_vortex.and.(.not.is_curv3)) call sol_vortex
  if (SRC.and.is_pulse.and.(.not.is_curv3)) call sol_pulse
  
end subroutine solver
