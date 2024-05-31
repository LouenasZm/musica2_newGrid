!==============================================================================
module mod_pp_main
!==============================================================================
  !>  Main post-processing module
!==============================================================================
  use warnstop
  use mod_constant
  use mod_pp_var
  implicit none
  ! ---------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------

contains

  !============================================================================
  subroutine pp_main
  !============================================================================
    !> author: XG 
    !> date: March 2022
    !> Main subroutine for post-processing
    !============================================================================
    use mod_utils        ! <- for numchar
    use mod_io           ! <- for read_write_volume
    use mod_interface    ! <- for communication
    use mod_eos          ! <- for primitives_visc
    use mod_init_flow    ! <- for init_thermo
    use mod_pp_sp_main
    use mod_pp_stats_main
    use mod_pp_corr
    use mod_pp_fstt_main
    use mod_pp_vorticity 
    use mod_pp_extend
    use mod_pp_convert
    implicit none
    ! -------------------------------------------------------------------------
    ! -------------------------------------------------------------------------

    ! Selection of corresponding snapshot number
    ! ------------------------------------------
    if (type_data.eq.0) then
       nsr = ipt_2_isn(nplr)
    else if (type_data.eq.1) then
       nsr = ili_2_isn(nplr)
    else if (type_data.eq.2) then
       nsr = ipl_2_isn(nplr)
    else if (type_data.eq.3) then
       nsr = ivl_2_isn(nplr)
    endif

    select case (type_pp)
    case(1) ! 1: spectra
       call pp_spectra
    case(2) ! 2: budgets
       !call mpistop('not defined yet', 0)
       call pp_stats_main
    case(3) ! 3: vorticity
       !call init_thermo
       allocate( rho(nx1:nx2,ny1:ny2,nz1:nz2))
       allocate(rhou(nx1:nx2,ny1:ny2,nz1:nz2))
       allocate(rhov(nx1:nx2,ny1:ny2,nz1:nz2))
       allocate(rhow(nx1:nx2,ny1:ny2,nz1:nz2))
       allocate(rhoe(nx1:nx2,ny1:ny2,nz1:nz2))
       is_mean0=.false.
       is_TamDong=.false. ! avoid reading restartTD.bin
       call read_write_volume(trim(dirDATA)//'restart_bl'//trim(numchar(nob(iproc)))//filext_read, &
            'restartTD.bin',READ)
       call communication_(rho,rhou,rhov,rhow,rhoe)
       ! apply angular periodicity if necessary
       call bc_angular_periodicity

       !call pp_clim
       if (iproc==0) print *,'Set primitive variables ...'
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

       !call mpistop('here',0)
       call init_thermo
       call primitives_visc!(rho,rhou,rhov,rhow,rhoe)
       deallocate(cok,visc,c_)

       ! if (iproc==0) print *,'pp_vorticity ...'
       call calc_vort_crit

       ! ww fluctuations
       !call write_ww('ww_bl'//trim(numchar(nob(iproc)))//filext_write)

    case(4) ! 4: correlation
       call pp_correl
    case(5) ! 5: PDF
       call mpistop('not defined yet', 0)
    case(6) ! 6: modal energy
       call mpistop('not defined yet', 0)
    case(7) ! 7: FST induced transition
       call pp_fstt
    case(8) ! 8: Extension of 1 point
       call pp_extend_main
    case(9) ! 9: Conversion binary <-> tecplot
       call pp_convert
    case default
       call mpistop('bad choice: post-processing not defined', 0)
    end select
    
    call MPI_BARRIER(COMM_global,info)

    ! Closing sequence
    ! ================
    if (iproc==0) write(*,'(A)') repeat("-",80)
    if (iproc==0) write(*,*) 'End of computation..'

    t_end = MPI_WTIME()
    call MPI_BARRIER(COMM_global,info)
    if (iproc.eq.0) write(*,'(A,i4,A,g0)') 'proc',iproc,' t_end-t_start ',t_end-t_start

!!$    !!call lib_tab
!!$    call MPI_COMM_FREE(COMM_global, info)

    call mpistop('', 0)

  end subroutine pp_main

  !============================================================================
  subroutine read_param_pp
  !============================================================================
    !> Determine post-processing parameters
    !============================================================================
    use mod_pp_fstt_main
    implicit none
    ! -------------------------------------------------------------------------
    ! -------------------------------------------------------------------------

    ! Directories (default)
    ! ===========
    !dirDATA='./'
    dirRESU='./'

    ! Read post-processing parameters
    ! ===============================
    open(30,file='param_pp.ini')
    
    read(30,*) ! =====================================================================================================
    read(30,*) ! SELECT POST-PROCESSING
    read(30,*) ! ======================
    read(30,*) ! Possible postprocessing functions:
    read(30,*) ! 1: spectra
    read(30,*) ! 2: budgets
    read(30,*) ! 3: vorticity
    read(30,*) ! 4: correlation
    read(30,*) ! 5: PDF
    read(30,*) ! 6: modal energy
    read(30,*) ! 7: FSTT
    read(30,*) ! 8: extend restart,stats,plane
    read(30,*) !
    read(30,*) ! 1/ you need to provide regular input files: param.ini, param_block.ini, [param_lines.ini,] info.ini
    read(30,*) ! 2/ set number of procs per direction in param_block.ini
    read(30,*) ! 3/ for planes post-processing: informations for plane definition
    read(30,*) !                                ~> in param_block.ini (nb, index, vars, ...)
    read(30,*) !                                ~> in param.ini (freq_planes, is_timestamp, ...)
    read(30,*) ! =====================================================================================================
    read(30,*) ! Select the type of postprocessing: type_pp
    read(30,*) type_pp
    read(30,*) ! Select the type of data: type_data=[0: points, 1: lines; 2: planes; 3: volumes; 4: stats; 5: restart]
    read(30,*) type_data
    read(30,*) ! Number of blocks read: nbloc_r
    read(30,*) nbloc_r
    read(30,*) ! Block number(s)
    allocate(nblr(nbloc_r))
    read(30,*) nblr
    ! restrict blocks to read blocks
    ! ------------------------------
    iblc_pp = nblr(1) ! Correct definition if post-processing only on 1 block ~> something different is made in init_pp_fstt, to be generalized if necessary
    if ((nbloc_r>1).and.(type_pp.lt.7).and.(type_pp.ne.3)) call mpistop('Except for mod_pp_fstt and vorticity, multiblock not available for post-processing',0)
    if (nbloc_r<nbloc) call restrict_blocks


!!$    print *,nbloc
!!$    do n=1,nbloc
!!$    print *,bl(n).ni,bl(n).nj,bl(n).nk
!!$    enddo
!!$    call mpistop('check',0)

    read(30,*) ! volume/plane/line/point number: nplr
    read(30,*) nplr
    read(30,*) ! Directory where data are stored: dirDATA
    read(30,*) !dirDATA
    read(30,*) ! Directory where results are stored: dirRESU
    read(30,*) dirRESU
    read(30,*) ! Name to append to output files
    read(30,*) name_output
    
    read(30,*) ! =====================================================================================================
    read(30,*) ! PARAMETERS for SPECTRA
    read(30,*) ! ======================
    read(30,*) ! number of samples: set by nmax in param.ini
    read(30,*) !
    read(30,*) ! list of possible directions: t, x,y,z [i,j,k in curvilinear]
    read(30,*) ! list of possible variables: prs,uu,vv,ww,rho,Tmp,Frhov,Grhow (,udf)
    read(30,*) !
    read(30,*) ! windowing:
    read(30,*) ! 1/ type_win='T': Tukey window; param_win in [0.,1.] (0.:constant window ~~> 1.:Hann window)
    read(30,*) !    for  Hann window, set type_win='T' & param_win=1.
    read(30,*) !    for rect. window, set type_win='T' & param_win=0.
    read(30,*) ! 2/ type_win='K': Kaiser-Bessel window; param_win=alpha (narrow for large alpha) 
    read(30,*) !    recommanded value for alpha: 3. (3.*pi)
    read(30,*) !    require AMOS library for Bessel functions
    read(30,*) !
    read(30,*) ! =====================================================================================================

    read(30,*) ! dimension: dim_spec=[1: 1D-spectra; 2: 2D-spectra; 3: 3D-spectra]
    read(30,*) sp%dim
    
    allocate(sp%d(sp%dim))
    sp%is_overlap=.false.
    read(30,*) ! directions for which PSD is applied: dir_spec(1:dim_spectra) 
    read(30,*) sp%d%name
    read(30,*) ! number of directions for averaging spectra: ndir_averag 
    read(30,*) ndir_av
    allocate(dir_av(ndir_av),i_av(ndir_av)) 
    read(30,*) ! directions for averaging spectra: dir_averag(1:ndir_averag) 
    read(30,*) dir_av
    read(30,*) ! number variables: nvar
    read(30,*) sp%nvar
    read(30,*) ! variables: varname(1:nvar) [see list of possible variables]
    read(30,*) sp%varname(1:sp%nvar)
    read(30,*) ! size of PSD blocks: lbloc(1:dim) [set to 0 if the full length is taken]
    read(30,*) sp%d%lbloc
    read(30,*) ! overlapping (Welch method): is_overlap
    read(30,*) sp%is_overlap
    read(30,*) ! define overlapping length: loverlap
    read(30,*) sp%loverlap
    read(30,*) ! windowing type: type_win(1:dim_spec)
    read(30,*) sp%d%type_win
    read(30,*) ! parameter for windowing: param_win(1:dim_spec)
    read(30,*) sp%d%param_win
    read(30,*) ! is_onesided: T: one-sided spectra / F: two-sided spectra
    read(30,*) sp%is_onesided
    read(30,*) ! use of Capon spectral estimator: is_capon (T: Capon / F: Fourier)
    read(30,*) sp%is_capon
    read(30,*) ! Capon filter length (Capon order): ncapon
    read(30,*) sp%ncapon
    read(30,*) ! modal analysis: output of selected modes
    read(30,*) is_modal
    read(30,*) ! number of slices for each direction [only for 2D or 3D spectra]
    read(30,*) n_sl
    if (sp%dim==1) n_sl=0
    read(30,*) ! slice indices for direction 1
    if (n_sl(1).ne.0) then
       allocate(i1_sl(n_sl(1)))
       read(30,*) i1_sl
    else
       read(30,*)
    endif
    read(30,*) ! slice indices for direction 2
    if (n_sl(2).ne.0) then
       allocate(i2_sl(n_sl(2)))
       read(30,*) i2_sl
    else
       read(30,*)
    endif
    read(30,*) ! slice indices for direction 3
    if (n_sl(3).ne.0) then
       allocate(i3_sl(n_sl(3)))
       read(30,*) i3_sl
    else
       read(30,*)
    endif
    if (n_sl(1)+n_sl(2)+n_sl(3)>0) then
       is_slice=.true.
    else
       is_slice=.false.
    endif
    read(30,*) ! initial time for samples: time_ini [useful if is_timestamp in param.ini]
    read(30,*) time_ini
    
    if (type_pp.eq.7) call init_pp_fstt
    ! if (type_pp.eq.7) call mpistop('mod_pp_fstt = lavoro in corso...',0)

  end subroutine read_param_pp

end module mod_pp_main
