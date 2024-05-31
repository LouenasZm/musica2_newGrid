!==============================================================================
module mod_pp_fstt_main
!==============================================================================
  !> Post-processing module for freestream turbulence induced transition (main)
!==============================================================================
  use warnstop
  use mod_pp_var
  use mod_mpi
  use mod_io_snapshots
  use mod_utils       ! <- for numchar
  use mod_pp_fstt_interpol
  use mod_pp_fstt_interpol_c
  use mod_pp_fstt_comm
  use mod_pp_fstt_ltbl
  use mod_pp_fstt_streaks
  use mod_pp_fstt_write
  implicit none
  ! ---------------------------------------------------------------------------
  integer :: dir_pp
  integer :: ialloc ! allocation cut in 2 parts
  integer(8) :: nout_prev
  logical :: is_calc_uut_uun
  integer(8), dimension(:), allocatable :: nout_array
  ! ---------------------------------------------------------------------------

contains

  !============================================================================
  subroutine pp_fstt
  !============================================================================
    !> author: AB
    !> date: May 2022
    !> Main subroutine for free-stream turbulence induced transition pp
    !============================================================================
    use mod_time        ! <- for deltat
    use mod_interface
    use mod_comm
    implicit none
    ! -------------------------------------------------------------------------
    logical :: is_compute_grad
    logical :: split_vol,is_3dOutputs
    integer :: i,j,k,l,i2,j2,isn,split_num
    integer :: ipt,ili,ipl,ivl
    integer(8) :: size_lim
    real(wp) :: n2_eta
    ! -------------------------------------------------------------------------


    split_num = 0
    split_vol = .false.
    ! split_vol = .true.
    size_lim = 3848290697000 ! Maximal file size on IRENE store (4 To, taken here equal to 3.5 To)

    is_check_stats = .true.

    is_3dOutputs = .false.
    ! is_3dOutputs = .true.

    ! Limitation for the detection of extremums to 2*delta_99
    is_lim_fst = .false.

    is_calc_uut_uun = .false.

    if (is_3dOutputs) then
       is_pp_discr=.false.
       is_pp_streaks=.false.
       is_check_stats=.false.
    endif
    if ((is_pp_discr).and.(is_curv)) is_calc_uut_uun = .true.
    ! if (is_pp_streaks) is_check_stats=.false.
    if ((is_curv).and.(is_pp_streaks))  call mpistop('streaks detection not implemented yet in curvilinear !',0)
    if ((is_curv).and.(is_lim_fst))  call mpistop('is_lim_fst not implemented yet in curvilinear !',0)

    if (iproc==0) print *,repeat('=',70)
    if (iproc==0) print *,"Post-processing for FST induced transition"
    if (iproc==0) print *,"------------------------------------------"
    if (iproc==0) print *,""

    if ((is_3dOutputs).and.(iproc==0))       print *,'PP mode ~> 3D outputs'
    if ((split_vol).and.(iproc==0))          print *,'PP mode ~> Splitting of saved volumes'
    if ((is_pp_streaks).and.(iproc==0))      print *,'PP mode ~> streaks detection'
    if ((is_pp_streaks).and.(is_pp_discr).and.(iproc==0)) print *,'        ~> laminar-turbulent discrimination'
    if ((is_pp_streaks).and.(.not.is_pp_discr)) call mpistop('Detection of streaks with turbulent laminar discrimination read from a file not implemented yet !',0)
    if ((.not.is_pp_streaks).and.(is_pp_discr).and.(iproc==0)) print *,'PP mode ~> laminar-turbulent discrimination'
    if ((is_check_stats).and.(iproc==0))     print *,'        ~> Global statistics check-up activated'
    if (iproc==0) print *,""

    ! Attribution of correct interpolated dimensions
    ! ==============================================
    nit_interp = nit_interp_g(iblc_pp+1-nblr(1)); nil_interp = nil_interp_g(iblc_pp+1-nblr(1))
    nj_interp = nj_interp_g(iblc_pp+1-nblr(1)); nkt_interp = nkt_interp_g(iblc_pp+1-nblr(1))
    nk_interp = nk_interp_g(iblc_pp+1-nblr(1))

    ! Allocation of the necessary arrays #1
    ! ==================================
    ialloc = 1
    call alloc_pp_fstt

    ! Initialisation of interpolation for FSTT
    ! ========================================
    if (iproc==0) print *,"Initialisation of interpolation..."
    if (is_curv) then
       call init_fstt_interpol_c
    else
       call init_fstt_interpol
    endif

    ! Initialisation of communication for interpolated var
    ! ====================================================
    if (iproc==0) print *,""
    if (iproc==0) print *,"Initialisation of specific communication for fstt..."
    call mpi_types_comm_fstt

    ! Time step computation
    ! =====================
    ! call timestep

    if (.not.is_curv) then
       ! Calculation of Bl values
       ! ========================
       u_lim = 0.95_wp*U_ref

       ! Calculation of U_max at each streamwise position
       do i=1,nx
          ! U0e(i) = maxval(stats_full(i,:,2))
          U0e(i) = U_ref
       enddo

       ! Calculation of 99% thickness BL
       do i=1,nx
          j=1
          do while ((stats_full(i+nx*coord(1),j,2).lt.0.99_wp*U0e(i)).and.(j.lt.ngy))
             j=j+1
          enddo
          d99(i) = yg(j)
       enddo

       ! Calculation of j_lim for ltbl
       do i=1,nx
          j=2
          do while ((y_interp(j).lt.1.5*d99(i)).and.(j.lt.nj_interp))
             j=j+1
          enddo
          j_1p5d99(i) = j
       enddo
    endif

    ! Filled with value from entire stats1 file
    ! -----------------------------------------
    do l=1,23
       do j=1,ny
          j2 = j + bl_glob(nob(iproc))%snapshot(nsr)%ind_j1 - 1 + coord(2)*ny
          do i=1,nx
             i2 = i + bl_glob(nob(iproc))%snapshot(nsr)%ind_i1 - 1 + coord(1)*nx
             stats_proc(i,j,l) = stats_full(i2,j2,l)
          enddo
       enddo
    enddo

    ! Allocation of the necessary arrays #2
    ! ==================================
    ialloc = 2
    call alloc_pp_fstt

    ! Interpolation of stats
    ! ======================
    if (is_curv) then
       call fstt_interpol_3d_c(1)
    else
       call fstt_interpol_3d(1)
    endif

    ! Initial settings for sub-volume reading
    ! =======================================
    ! Setting initial offset for sub-volume
    snapshots(nsr)%tectype%offset = 0_MPI_OFFSET_KIND
    ! define offset for 1 instance: number of points * number of variables * 8 bytes
    snapshots(nsr)%tectype%disp=snapshots(nsr)%tectype%ngx ! in 2 steps to correctly have disp (> max(integer(4)) sometimes)
    snapshots(nsr)%tectype%disp=snapshots(nsr)%tectype%disp &
                                *snapshots(nsr)%tectype%ngy &
                                *snapshots(nsr)%tectype%ngz*8*snapshots(nsr)%nvar
    ! Not a check volume ~> forced to .false. (but should already be .false. anyway)
    snapshots(nsr)%is_check = .false.

    ! Initialisation of the kernel coefficients
    ! =========================================
    if (is_curv) then
       call init_kernel_c
    else
       call init_kernel
    endif

    ! Initialisation of the streaks detection subroutine
    ! ==================================================
    if (is_pp_streaks) call init_streaks_detection

    ! Discriminated stats
    ! ===================
    ! Prepare a xy-plane
    ! ------------------
    ! planes_stats_pp_fstt%normal=3
    ! planes_stats_pp_fstt%index=1
    planes_stats_pp_fstt%ind_i1=1; planes_stats_pp_fstt%ind_i2=ngx
    planes_stats_pp_fstt%ind_j1=1; planes_stats_pp_fstt%ind_j2=ngy
    planes_stats_pp_fstt%ind_k1=1; planes_stats_pp_fstt%ind_k2=1
    planes_stats_pp_fstt%freq=0

    ! plane variables (number & names)
    planes_stats_pp_fstt%nvar=1
    allocate(planes_stats_pp_fstt%var(planes_stats_pp_fstt%nvar))
    planes_stats_pp_fstt%var='mean'
    ! plane timestamp
    planes_stats_pp_fstt%stamp=is_timestamp

    write(planes_stats_pp_fstt%tectype%zonename,'(A5)') 'stats'
    planes_stats_pp_fstt%tectype%strandID = -2

    ! Not a check plane ~> .false.
    planes_stats_pp_fstt%is_check = .false.

    call init_snapshot(planes_stats_pp_fstt)

    ! To check the interpolated field
    ! ===============================
    if (is_check_fstt_vl) then
       ! Prepare a volume check
       ! ----------------------
       ! Fixed here
       volume_check_fstt%nvar = 1
       allocate(volume_check_fstt%var(volume_check_fstt%nvar))
       volume_check_fstt%var = 'check'

       ! Volume check over interpolated mesh
       volume_check_fstt%ind_i1 = 1
       volume_check_fstt%ind_i2 = nit_interp
       volume_check_fstt%ind_j1 = 1
       volume_check_fstt%ind_j2 = nj_interp
       volume_check_fstt%ind_k1 = 1
       volume_check_fstt%ind_k2 = nkt_interp

       ! volume timestamp
       volume_check_fstt%stamp=is_timestamp

       volume_check_fstt%tectype%strandID = -30

       ! Not a check volume ~> .false.
       volume_check_fstt%is_check = .false.

       ! Appel init_volume_check car maillage interpolation
       call init_volume_check(volume_check_fstt)

       ! inform about restart mode
       volume_check_fstt%tectype%restart=.false.
    endif

    ! Snapshot outputs
    ! ================
    is_compute_grad = .false.
    if (is_snap_pp_fstt) then
       do isn=1,nsnap_pp_fstt
          ! Prepare volume ouputs
          ! ---------------------
          snap_pp_fstt(isn)%tectype%is_IOtec_write = .false.
          ! Snapshot timestamp: put to true to write tecplot file
          snap_pp_fstt(isn)%tectype%strandID = isn
          write(snap_pp_fstt(isn)%tectype%zonename,'(A5,I3.3)') 'plane',isn

          ! Not a check snapshot ~> .false.
          snap_pp_fstt(isn)%is_check = .false.

          ! Initialisation snapshot
          call init_snapshot(snap_pp_fstt(isn))

          ! Snapshot timestamp: put to true to write tecplot file
          if (split_vol) then
             split_num = 1
             snap_pp_fstt(isn)%stamp = .false.
          else if (snap_pp_fstt(isn)%type.eq.3) then
             split_num = 0
             snap_pp_fstt(isn)%stamp = .true.
          else if ((snap_pp_fstt(isn)%type.eq.2).and.(is_3dOutputs)) then
             split_num = 0
             snap_pp_fstt(isn)%stamp = .true.
          else
             split_num = 0
             snap_pp_fstt(isn)%stamp = .false.
          endif
       enddo

       ! Save link between snapshot num. and vol,... + attribution of new numbers
       ! ------------------------------------------------------------------------
       allocate(ipt_2_isn_pp(1:npoints_pp))
       allocate(ili_2_isn_pp(1:nlines_pp))
       allocate(ipl_2_isn_pp(1:nplanes_pp))
       allocate(ivl_2_isn_pp(1:nvolumes_pp))
       ipt=0;ili=0;ipl=0;ivl=0
       do isn=1,nsnap_pp_fstt
          if (snap_pp_fstt(isn)%type.eq.0) then
             ipt = ipt + 1
             ipt_2_isn_pp(ipt) = isn
             write(snap_pp_fstt(isn)%tectype%zonename,'(A5,I3.3)') 'point',ipt
          else if (snap_pp_fstt(isn)%type.eq.1) then
             ili = ili + 1
             ili_2_isn_pp(ili) = isn
             write(snap_pp_fstt(isn)%tectype%zonename,'(A5,I3.3)') 'line',ili
          else if (snap_pp_fstt(isn)%type.eq.2) then
             ipl = ipl + 1
             ipl_2_isn_pp(ipl) = isn
             write(snap_pp_fstt(isn)%tectype%zonename,'(A5,I3.3)') 'plane',ipl
          else if (snap_pp_fstt(isn)%type.eq.3) then
             ivl = ivl + 1
             ivl_2_isn_pp(ivl) = isn
             write(snap_pp_fstt(isn)%tectype%zonename,'(A5,I3.3)') 'volume',ivl
          endif
       enddo

       do isn=1,nsnap_pp_fstt
          ! Part added to allow appending in the same file
          ! ==============================================
          if (.not.snap_pp_fstt(isn)%stamp) then
             ! logical is_app (is append_file)
             snap_pp_fstt(isn)%tectype%is_app=.true.
             ! define offset: number of points * number of variables * 8 bytes
             snap_pp_fstt(isn)%tectype%disp=snap_pp_fstt(isn)%tectype%ngx ! in 2 steps to correctly have disp (> max(integer(4)) sometimes)
             snap_pp_fstt(isn)%tectype%disp=snap_pp_fstt(isn)%tectype%disp &
                                           *snap_pp_fstt(isn)%tectype%ngy &
                                           *snap_pp_fstt(isn)%tectype%ngz*8*snap_pp_fstt(isn)%nvar

             ! inform about restart mode
             if (is_restart_pp) then
                ! To allow appending to the file
                snap_pp_fstt(isn)%tectype%restart=.true.
             else
                ! If not is_restart_pp, suppression of snapshot
                snap_pp_fstt(isn)%tectype%restart=.false.
                ! To remove if file exists already
                if (.not.snap_pp_fstt(isn)%stamp) call check_snapshot_pp_fstt(isn)
             endif
          endif

          ! Determine if necessity to derivate velocities (for vorticity)
          do l=1,snap_pp_fstt(isn)%nvar
             if ((trim(snap_pp_fstt(isn)%var(l)).eq.'l2').or.&
                 (trim(snap_pp_fstt(isn)%var(l)).eq.'QQ').or.&
                 (trim(snap_pp_fstt(isn)%var(l)).eq.'vort')) is_compute_grad=.true.
          enddo

          ! Determine wall normals if necessary (for tangential and normal velocities)
          do l=1,snap_pp_fstt(isn)%nvar
             if ((trim(snap_pp_fstt(isn)%var(l)).eq.'uut').or.&
                 (trim(snap_pp_fstt(isn)%var(l)).eq.'uun')) then
                if (.not.is_curv) call mpistop('uut equivalent to uu in cartesian ~> put var to "uu"',0)
                if (bl(nob(iproc))%BC(3).ne.0) call mpistop('Tangential velocity for wall not in jmin is not implemented yet !',iproc)
                ! if not done in grid_normals, done here
                ! --------------------------------------
                if (.not.allocated(nxn_jmin)) then
                   allocate(nxn_jmin(1:nx,1:nz),nyn_jmin(1:nx,1:nz))
                   do k=1,nz
                      do i=1,nx
                         ! J²*||grad(eta)||²
                         n2_eta=eta_x(i,1,k)**2+eta_y(i,1,k)**2+eta_z(i,1,k)**2
                         ! BC normal components: J*grad(eta)/J*||grad(eta)||
                         nxn_jmin(i,k) = eta_x(i,1,k)/sqrt(n2_eta)
                         nyn_jmin(i,k) = eta_y(i,1,k)/sqrt(n2_eta)
                      enddo
                   enddo
                endif
                is_calc_uut_uun = .true.
             endif
          enddo
       enddo
    endif

    ! if vorticity calculated, allocation of derivative arrays
    if (is_compute_grad) then
       allocate(dux(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
       allocate(dvx(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
       allocate(dwx(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
       allocate(duy(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
       allocate(dvy(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
       allocate(dwy(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
       allocate(duz(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
       allocate(dvz(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
       allocate(dwz(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v))
    endif

    ! Allocation if necessary
    ! -----------------------
    if (is_calc_uut_uun) then
       allocate(uut(nx1:nx2,ny1:ny2,nz1:nz2))
       allocate(uun(nx1:nx2,ny1:ny2,nz1:nz2))
    endif

    ! ************************************
    ! Loop over number of outputs selected
    ! ************************************
    ! Total number of outputs for sub-volume
    noutputs = nmax/freq_volume
    ! Output number where to begin PP
    nbeg_pp = nit_beg/freq_volume + 1
    ! Output number where to end PP
    nend_pp = nit_end/freq_volume

    ! Check to see if loop bounds coherent
    if (nbeg_pp.gt.nend_pp) call mpistop('Post-processing loop bounds in param_pp_fstt.ini incoherent',0)
    ! Check to see if enough outputs for the loop bounds selected
    if (nend_pp.gt.noutputs) call mpistop('Bounds of Post-processing loop greater than number of total subvolume'//&
                                          ' outputs indicated in param.ini',0)
    ! Check to see if restart at minimal born => incoherent
    if ((is_restart_pp).and.(nbeg_pp.eq.1)) call mpistop('Restart post-processing is true while minimal iteration'//&
                                                         ' indicated is 1...',0)
    ! If checks ok, offset for reading subvolume modified to correspond to output number nbeg_pp
    snapshots(nsr)%tectype%offset = snapshots(nsr)%tectype%disp*(nbeg_pp-1)

    ! Generation of the list of volumes to read
    if (is_3dOutputs) then
       if (iproc==0) print *,"Reading vols2read.bin..."
       open(30,file='vols2read.dat',form="formatted")
       rewind(30)
       read(30,*) noutputs
       allocate(nout_array(noutputs))
       do k=1,noutputs
          read(30,*) nout_array(k)
       enddo
       close(30)
    else
       ! Number of outputs selected for post-processing
       noutputs = nend_pp - nbeg_pp + 1
       allocate(nout_array(noutputs))
       ! loop over noutputs
       do k=1,noutputs
          nout_array(k) = k + nbeg_pp - 1
       enddo
    endif

    ! snapshot frequency output
    ! -------------------------
    if (is_snap_pp_fstt) then
       do isn=1,nsnap_pp_fstt
          if (snap_pp_fstt(isn)%freq.eq.0) then
             snap_pp_fstt(isn)%freq = 1
          else if (snap_pp_fstt(isn)%freq.lt.0) then
             snap_pp_fstt(isn)%freq = - noutputs/snap_pp_fstt(isn)%freq
             if (snap_pp_fstt(isn)%freq.eq.0) snap_pp_fstt(isn)%freq=noutputs
          endif
       enddo
    endif

    if (is_restart_pp) then
       ! Reading of the discriminated stats
       ! ----------------------------------
       call read_write_stats_pp_fstt(READ)

       ! Initialisation of the cpt for averaging statistics
       ! --------------------------------------------------
       do j=1,ny
          do i=1,nx
             stats_cpt(1,i,j) =  stats_lam(i,j,1)*dble(nbeg_pp-1)
             stats_cpt(2,i,j) = stats_turb(i,j,1)*dble(nbeg_pp-1)
          enddo
       enddo
    endif


    if (iproc==0) print *,""
    if (iproc==0) write(*,*) 'Init ok. Starting post-processing loop...'
    if (iproc==0)    print *,"----------------------------------------"
    if (iproc==0) print *,""

    nout_prev = nout_array(1) - 1
    do nout=1,noutputs
    ! do nout=1,61,15
       ! start measure of elapsed CPU time
       ! ---------------------------------
       time_before=MPI_WTIME()

       ! Volume to be read
       ! -----------------
       ivol_2read = nout_array(nout)

       ! Reading of sub-volume
       ! ---------------------
       ! Addition of an offset
       snapshots(nsr)%tectype%offset = snapshots(nsr)%tectype%offset + (ivol_2read-(nout_prev+1))*snapshots(nsr)%tectype%disp

       ! MPI-IO read
       call read_snapshot(nplr,iblc_pp,dirDATA)

       ! Temporary <- to comment if necessary or split_vol=.false.
       ! Writting sub-volume in multiple files
       ! call write_subvol_mult
       if (split_vol) then
          if (snap_pp_fstt(1)%tectype%offset.ge.size_lim) then
             split_num = split_num + 1
             snap_pp_fstt(1)%tectype%offset = 0_MPI_OFFSET_KIND
          endif
          call write_snapshot_pp_fstt(1,split_num)
          if ((iproc==0).and.(mod(nout,10).eq.0)) print *,nout,noutputs
          cycle
       endif

       if (is_curv) then
          ! Calculation of normal and tangential velocity
          if (is_calc_uut_uun) call calc_uut_uun

          if (is_pp_discr) then
             ! Communication for interpolation
             call communication_3d(uu,vv,ww,uut,uun)

             ! Interpolation on rougher mesh
             ! -----------------------------
             call fstt_interpol_3d_c(2)

             ! Calculation of fluctuating field
             ! --------------------------------
             call fstt_fluct_field

             ! Communication of values on gh
             ! -----------------------------
             call communication_fstt3(uut_interp,uun_interp,ww_interp)
          else if (is_compute_grad) then
             ! Communication for vorticity
             call communication_3d(uu,vv,ww,uut,uun)
          endif
       else
          if ((is_pp_streaks).or.(is_pp_discr)) then
             ! Communication for interpolation
             call communication_3d(uu,vv,ww,uu,vv)

             ! Interpolation on rougher mesh
             ! -----------------------------
             call fstt_interpol_3d(2)

             ! Calculation of BL edge
             ! ----------------------
             if ((is_pp_streaks).or.(is_lim_fst)) call fstt_bl_edge

             ! Calculation of fluctuating field
             ! --------------------------------
             call fstt_fluct_field

             ! Communication of values on gh
             ! -----------------------------
             call communication_fstt3(uu_interp,vv_interp,ww_interp)

          endif

          if ((is_compute_grad).or.(is_pp_streaks).or.&
              (is_pp_discr))        call communication_3d(uu,vv,ww,uu_fluct,uu)
       endif


       if (is_pp_discr) then
          ! Determination of laminar and turbulent regions -> stored in ltbl
          ! ----------------------------------------------
          ! Computation of extremums density distribution
          call fstt_lam_turb

          ! Interpolation of density over the computationnal grid
          if (is_curv) then
             call fstt_interpol_3d_inv_c
          else
             call fstt_interpol_3d_inv
          endif

          ! Application of thresold (Otsu, 1979)
          call fstt_thresold

          ! Calculation of discriminated statistics
          ! ---------------------------------------
          call fstt_stats_ltbl

          ! Volume to check on interpolated grid
          ! ------------------------------------
          if (is_check_fstt_vl) call write_volume_check

       endif

       ! Snapshot outputs
       ! ----------------
       if (is_snap_pp_fstt) then
          if (is_compute_grad) call grad_vel
          do isn=1,nsnap_pp_fstt
             if (mod(nout,snap_pp_fstt(isn)%freq).eq.0) call write_snapshot_pp_fstt(isn,split_num)
          enddo
       endif

       ! Deactivated for volumes output
       if (is_pp_streaks) then
          ! Local streaks detection
          ! -----------------------
          call local_streaks_detection

          ! Communication of streaks to master proc, merge & writting
          ! ---------------------------------------------------------
          call merging_streaks
       endif

       ! Save number of last volume read
       nout_prev = ivol_2read

       ! end measure of elapsed CPU time
       ! -------------------------------
       time_after = MPI_WTIME()

       cputime = time_after - time_before
       cputot = cputot + cputime
       cpurun = cpurun + cputime

       if ((iproc==0).and.(mod(nout,1).eq.0)) write(*,20) nout,noutputs,nout+nbeg_pp-1,nend_pp,int(cpurun)
20     format(10x,'~> Iterations:   ',i0,' / ',i0,' (Tot: ',i0,' / ',i0,'), Run time: ',i0,' s')
       if ((iproc==0).and.(mod(nout,1).eq.0)) print *,""

    enddo

    if (iproc.eq.0) print *,''
    if (iproc.eq.0) write(*,*) 'End of Post-processing Loop'
    if (iproc.eq.0)    print *,"---------------------------"
    if (iproc.eq.0) print *,""

    ! Writting of the discriminated stats
    ! --------------------------------
    call read_write_stats_pp_fstt(WRITE)

    call mpistop("End of post-processing for FSTT novec.",0)

  end subroutine pp_fstt

  !============================================================================
  subroutine init_pp_fstt
  !============================================================================
    !> Set particular dimensions of domain for PP FST induced transition
    !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: isn,ib,ny_v,ib_,count_vol
    logical :: iexist
    ! -------------------------------------------------------------------------
    integer, dimension(MPI_STATUS_SIZE) :: statut
    integer :: fh

    if (iproc==0) print *,'FSTT post-processing: dimensions set to subvolume'
    if (iproc==0) print *,'--------------------'

    ! is_curv forced
    if (LE) is_curv=.true.

    ! Read parameter file for PP FST
    ! ------------------------------
    call read_param_pp_fstt

    ! Save of the global block values
    ! ===============================
    allocate(bl_glob(1:nbloc_r))
    do ib=1,nbloc_r
       bl_glob(ib)%ni    = bl(ib)%ni
       bl_glob(ib)%nj    = bl(ib)%nj
       bl_glob(ib)%nk    = bl(ib)%nk
       bl_glob(ib)%ndomi = bl(ib)%ndomi
       bl_glob(ib)%ndomj = bl(ib)%ndomj
       bl_glob(ib)%ndomk = bl(ib)%ndomk
       bl_glob(ib)%nproc = bl_glob(ib)%ndomi*bl_glob(ib)%ndomj*bl_glob(ib)%ndomk

       ! Save snapshots
       bl_glob(ib)%nsnapshot = bl(ib)%nsnapshot
       allocate(bl_glob(ib)%snapshot(bl_glob(ib)%nsnapshot))
       do isn=1,bl_glob(ib)%nsnapshot
          bl_glob(ib)%snapshot(isn)%ind_i1 = bl(ib)%snapshot(isn)%ind_i1
          bl_glob(ib)%snapshot(isn)%ind_i2 = bl(ib)%snapshot(isn)%ind_i2
          bl_glob(ib)%snapshot(isn)%ind_j1 = bl(ib)%snapshot(isn)%ind_j1
          bl_glob(ib)%snapshot(isn)%ind_j2 = bl(ib)%snapshot(isn)%ind_j2
          bl_glob(ib)%snapshot(isn)%ind_k1 = bl(ib)%snapshot(isn)%ind_k1
          bl_glob(ib)%snapshot(isn)%ind_k2 = bl(ib)%snapshot(isn)%ind_k2
          bl_glob(ib)%snapshot(isn)%freq   = bl(ib)%snapshot(isn)%freq
          bl_glob(ib)%snapshot(isn)%nvar   = bl(ib)%snapshot(isn)%nvar
          bl_glob(ib)%snapshot(isn)%var(1:bl(ib)%snapshot(isn)%nvar) = bl(ib)%snapshot(isn)%var(1:bl(ib)%snapshot(isn)%nvar)
       enddo
    enddo

    ! Identification of the correct nsr <~ Done later at the beginning of mod_pp_main.f90
    ! =================================
    nsr=0; count_vol=0
    loopnsr: do isn=1,bl_glob(1)%nsnapshot
       nsr = nsr+1
       if ((bl(1)%snapshot(isn)%ind_i1.ne.bl(1)%snapshot(isn)%ind_i2).and.&
           (bl(1)%snapshot(isn)%ind_j1.ne.bl(1)%snapshot(isn)%ind_j2).and.&
           (bl(1)%snapshot(isn)%ind_k1.ne.bl(1)%snapshot(isn)%ind_k2)) then
          count_vol = count_vol+1
          if (count_vol.eq.nplr) exit loopnsr
       endif
    enddo loopnsr

    ! Modification of subvolume indices considered for pp
    ! ===================================================
    do ib=1,nbloc_r
       bl(ib)%snapshot(nsr)%ind_i1 = 1
       bl(ib)%snapshot(nsr)%ind_i2 = bl_glob(ib)%snapshot(nsr)%ind_i2 - bl_glob(ib)%snapshot(nsr)%ind_i1 + 1
       bl(ib)%snapshot(nsr)%ind_j1 = 1
       bl(ib)%snapshot(nsr)%ind_j2 = bl_glob(ib)%snapshot(nsr)%ind_j2 - bl_glob(ib)%snapshot(nsr)%ind_j1 + 1
       bl(ib)%snapshot(nsr)%ind_k1 = 1
       bl(ib)%snapshot(nsr)%ind_k2 = bl_glob(ib)%snapshot(nsr)%ind_k2 - bl_glob(ib)%snapshot(nsr)%ind_k1 + 1
    enddo


    do ib=1,nbloc_r
       ! Modification of block dimensions to subvolume considered for pp
       ! ===============================================================
       ! Verification of sub-volume existence
       ! ------------------------------------
       if (bl(ib)%nsnapshot.lt.nsr) call mpistop("  ~> Problem in FSTT PP: selected volume number doesn't exist",0)

       ! Verification of sub-volume size
       ! -------------------------------
       if ((bl(ib)%snapshot(nsr)%ind_i1).ge.(bl(ib)%snapshot(nsr)%ind_i2)) call mpistop('  ~> Problem for sub-volume definition: I1 greater or equal to I2',0)
       if ((bl(ib)%snapshot(nsr)%ind_j1).ge.(bl(ib)%snapshot(nsr)%ind_j2)) call mpistop('  ~> Problem for sub-volume definition: J1 greater or equal to J2',0)
       if ((bl(ib)%snapshot(nsr)%ind_k1).ge.(bl(ib)%snapshot(nsr)%ind_k2)) call mpistop('  ~> Problem for sub-volume definition: K1 greater or equal to K2',0)

       ! Verification if only 1 proc in j direction
       ! ------------------------------------------
       if (bl(ib)%ndomj.gt.1) call mpistop("  ~> Problem in FSTT PP: only 1 proc authorized in j-direction",0)

       ! Modification of block dimensions to sub-volume considered for pp
       ! ----------------------------------------------------------------
       bl(ib)%ni = bl(ib)%snapshot(nsr)%ind_i2 - bl(ib)%snapshot(nsr)%ind_i1 + 1
       bl(ib)%nj = bl(ib)%snapshot(nsr)%ind_j2 - bl(ib)%snapshot(nsr)%ind_j1 + 1
       bl(ib)%nk = bl(ib)%snapshot(nsr)%ind_k2 - bl(ib)%snapshot(nsr)%ind_k1 + 1
       ny_v = bl(ib)%nj/bl(ib)%ndomj

       ! Temporary: ny_interp fixed to be a divisor number of ny_v
       ndivy = ny_v/nj_interp_g(ib)
       nj_interp_g(ib) = ny_v/ndivy

       if (is_curv) then
          ! Temporary: nx_interp fixed to be a divisor number of bl(ib)%ni
          ndivx = bl(ib)%ni/nit_interp_g(ib)
          nit_interp_g(ib) = bl(ib)%ni/ndivx
          if (mod(nit_interp_g(ib),bl(ib)%ndomi).ne.0) call mpistop('WARNING: not possible (nit_interp must be divisible by ndomi)', 0)
          nil_interp_g(ib) = nit_interp_g(ib)/bl(ib)%ndomi
       endif
    enddo

    ! Find block corresponding to proc = ib
    ! =====================================

    iblc_pp=0
    do ib_=1,nbloc_r
       if ((iproc.le.bl(ib_)%proc_max).and.(iproc.ge.bl(ib_)%proc_min)) then
          iblc_pp=nblr(ib_)
          ib = ib_
       endif
    enddo
    if (iblc_pp.eq.0) call mpistop('Problem with attribution of iblc_pp to proc !',0)

    ! Creation of stats file to right size
    ! ====================================
    ! Verification of the file existence
    !-----------------------------------
    inquire(file='stats1_bl'//trim(numchar(iblc_pp))//'.bin',exist=iexist)
    if (.not.iexist) then
       call mpistop('stats1_bl'//trim(numchar(iblc_pp))//'.bin file not found !', 0)
    endif

    ! Reading all stats
    ! -----------------
    ! Temporary allocation for stats data
    allocate(stats_full(1:bl_glob(ib)%ni,1:bl_glob(ib)%nj,1:23))

    ! Reading of stats file
    call MPI_FILE_OPEN(COMM_global,'stats1_bl'//trim(numchar(iblc_pp))//'.bin', &
      MPI_MODE_RDONLY, MPI_INFO_NULL,fh,info)

    ! Reading var
    call MPI_FILE_READ(fh,stats_full(:,:,:),size(stats_full(:,:,:)),MPI_DOUBLE_PRECISION,statut,info)

    ! Closing file
    call MPI_FILE_CLOSE(fh,info)

    ! Beginning and end of post-processing zone
    idem_pp = bl_glob(ib)%snapshot(nsr)%ind_i1; iend_pp = bl_glob(ib)%snapshot(nsr)%ind_i2
    jdem_pp = bl_glob(ib)%snapshot(nsr)%ind_j1; jend_pp = bl_glob(ib)%snapshot(nsr)%ind_j2
    ! idem_pp = bl(ib)%snapshot(nsr)%ind_i1; iend_pp = bl(ib)%snapshot(nsr)%ind_i2
    ! jdem_pp = bl(ib)%snapshot(nsr)%ind_j1; jend_pp = bl(ib)%snapshot(nsr)%ind_j2

    ! Modification of BCs & sponge zone ~> put in hard for the moment
    ! ---------------------------------
    if (nbloc_r.gt.1) then
       bl(1)%BC(1)=-1; bl(1)%BC(2)=2;  bl(1)%BC(3)=0; bl(1)%BC(4)=-1; bl(1)%BC(5)=1; bl(1)%BC(6)=1 ! Block 7
       bl(2)%BC(1)=1;  bl(2)%BC(2)=-1; bl(2)%BC(3)=0; bl(2)%BC(4)=-1; bl(2)%BC(5)=2; bl(2)%BC(6)=2 ! Block 8

       bl(1)%is_sponge = "F"; bl(2)%is_sponge = "F"
    endif

    ! Snapshots ouputs attribution
    nsnap_pp_fstt = blc_snap_pp(ib)%nsn_pp
    if (nsnap_pp_fstt.gt.0) then
       is_snap_pp_fstt = .true.
       allocate(snap_pp_fstt(nsnap_pp_fstt))
       do isn=1,nsnap_pp_fstt
          allocate(snap_pp_fstt(isn)%var(10))
          snap_pp_fstt(isn)%ind_i1 = blc_snap_pp(ib)%snap(isn)%ind_i1
          snap_pp_fstt(isn)%ind_i2 = blc_snap_pp(ib)%snap(isn)%ind_i2
          snap_pp_fstt(isn)%ind_j1 = blc_snap_pp(ib)%snap(isn)%ind_j1
          snap_pp_fstt(isn)%ind_j2 = blc_snap_pp(ib)%snap(isn)%ind_j2
          snap_pp_fstt(isn)%ind_k1 = blc_snap_pp(ib)%snap(isn)%ind_k1
          snap_pp_fstt(isn)%ind_k2 = blc_snap_pp(ib)%snap(isn)%ind_k2
          snap_pp_fstt(isn)%freq   = blc_snap_pp(ib)%snap(isn)%freq
          snap_pp_fstt(isn)%nvar   = blc_snap_pp(ib)%snap(isn)%nvar
          snap_pp_fstt(isn)%var(1:snap_pp_fstt(isn)%nvar)= blc_snap_pp(ib)%snap(isn)%var(1:snap_pp_fstt(isn)%nvar)
       enddo
    endif
    deallocate(blc_snap_pp)

  end subroutine init_pp_fstt


  !============================================================================
  subroutine read_param_pp_fstt
  !============================================================================
    !> Determine post-processing parameters for FST
    !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: isn,ibb
    logical :: iexist
    ! -------------------------------------------------------------------------

    ! Read post-processing parameters for FST
    ! =======================================
    inquire(file='param_pp_fstt.ini',exist=iexist)
    if (.not.iexist) then
       call mpistop("Paramfile for FST PP doesn't exist !", 0)
    endif
    open(30,file='param_pp_fstt.ini')
    rewind(30)

    read(30,*) !=============================================================
    read(30,*) !=============================================================
    read(30,*) ! MUSICA2 : fill post-processing options for fstt
    read(30,*) !=============================================================
    read(30,*) !=============================================================
    read(30,*) ! Time iterations to begin and end post-processing (all if 0 0)
    read(30,*) nit_beg, nit_end
    if (nit_beg.le.0) nit_beg = 1
    if (nit_end.le.0) nit_end = nmax
    read(30,*) ! Restart post-processing (is_restart_pp)
    read(30,*) is_restart_pp
    read(30,*) !=============================================================
    read(30,*) ! Parameters for laminar-turbulent discrimination
    read(30,*) !=============================================================
    read(30,*) ! is_pp_discr (if F, read from volume outputs)
    read(30,*) is_pp_discr
    read(30,*) ! Number of points for interpolation
    allocate(nit_interp_g(nbloc_r))
    allocate(nil_interp_g(nbloc_r))
    allocate(nj_interp_g(nbloc_r))
    allocate(nkt_interp_g(nbloc_r))
    allocate(nk_interp_g(nbloc_r))
    do ib=1,nbloc_r
       read(30,*) nit_interp_g(ib), nj_interp_g(ib), nkt_interp_g(ib)
       if (mod(nit_interp_g(ib),bl(ib)%ndomi).ne.0) call mpistop('WARNING: not possible (nit_interp must be divisible by ndomi)', 0)
       nil_interp_g(ib) = nit_interp_g(ib)/bl(ib)%ndomi
       if (mod(nkt_interp_g(ib),bl(ib)%ndomk).ne.0) call mpistop('WARNING: not possible (nkt_interp must be divisible by ndomk)', 0)
       nk_interp_g(ib) = nkt_interp_g(ib)/bl(ib)%ndomk
    enddo
    read(30,*) ! Kernel: smoothing parameter (~ 1-2)
    read(30,*) h_kern
    read(30,*) ! Levels for the thresold
    read(30,*) lvl_th
    read(30,*) !=============================================================
    read(30,*) ! Parameters for streaks detections
    read(30,*) !=============================================================
    read(30,*) ! is_pp_streaks (to detect streaks in laminar boundary layer)
    read(30,*) is_pp_streaks
    read(30,*) ! Streaks direction: 1=i, 2=j
    read(30,*) dir_pp
    read(30,*) ! Definition of distances used in PP for FSTT
    read(30,*) ! Order advice: dist_max_extr < dist_max_xy < dist_max < l_min_str
    read(30,*) ! dist_max_extr: max. dist. btwn 2 extremums pts in y-z plane (m)
    read(30,*) dist_max_extr
    read(30,*) ! dist_max_yz: max. dist. btwn 2 consecutives streaks pts in y-z plane (m)
    read(30,*) dist_max_yz
    read(30,*) ! dist_max_x: max. streamwise dist. btwn 2 consecutives streaks tip (m)
    read(30,*) dist_max_x
    read(30,*) ! l_min_str: min. streaks length (m)
    read(30,*) l_min_str
    read(30,*) !=============================================================
    read(30,*) ! Supplementary outputs
    read(30,*) !=============================================================
    read(30,*) ! Volume check: T (True) / F (False)
    read(30,*) is_check_fstt_vl
    allocate(blc_snap_pp(nbloc_r))
    do ibb=1,nbloc_r
       read(30,*) !=============================================================
       read(30,*) ! Define output snapshots of block #iblc_pp
       read(30,*) !=============================================================
       read(30,*) blc_snap_pp(ibb)%nsn_pp
       read(30,*) !----|----|----|----|----|----|------|------|--------------
       read(30,*) ! I1 | I2 | J1 | J2 | K1 | K2 | freq | nvar | list var
       read(30,*) !    |    |    |    |    |    |      |   n  | name ...
       read(30,*) !----|----|----|----|----|----|------|------|--------------
       if (blc_snap_pp(ibb)%nsn_pp.gt.0) then
          ! allocate snapshot type
          allocate(blc_snap_pp(ibb)%snap(blc_snap_pp(ibb)%nsn_pp))
          ! read infos
          do isn=1,blc_snap_pp(ibb)%nsn_pp
             read(30,*) blc_snap_pp(ibb)%snap(isn)%ind_i1, &
                        blc_snap_pp(ibb)%snap(isn)%ind_i2, &
                        blc_snap_pp(ibb)%snap(isn)%ind_j1, &
                        blc_snap_pp(ibb)%snap(isn)%ind_j2, &
                        blc_snap_pp(ibb)%snap(isn)%ind_k1, &
                        blc_snap_pp(ibb)%snap(isn)%ind_k2, &
                        blc_snap_pp(ibb)%snap(isn)%freq,   &
                        blc_snap_pp(ibb)%snap(isn)%nvar,   &
                        blc_snap_pp(ibb)%snap(isn)%var(1:blc_snap_pp(ibb)%snap(isn)%nvar)
          enddo
       endif
    enddo

    close(30)

  end subroutine read_param_pp_fstt

  !============================================================================
  subroutine alloc_pp_fstt
  !============================================================================
    !> Allocation for FSTT
    !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: nx1v,nx2v,ny1v,ny2v,nz1v,nz2v,iv,nv
    ! -------------------------------------------------------------------------

    iv = nplr

    if (ialloc.eq.1) then
       ! Stats for each proc
       allocate(stats_proc(1:nx,1:ny,23))
       ! Stats usefull
       allocate(U0e(nx))
       allocate(d99(nx))
       allocate(j_1p5d99(nx))
       if (is_curv) then
          allocate( bk_ic(1:nil_interp,2:nj_interp,2))
          allocate( bk_jc(1:nil_interp,2:nj_interp,2))
       else
          allocate(born_kern_j(2:nj_interp,2))
       endif
    else if (ialloc.eq.2) then
       ! Deallocation of entire stats1 data
       ! ----------------------------------
       deallocate(stats_full)

       ! Allocation of stats interpolated
       ! --------------------------------
       allocate(stats_interp(1-ngh_pp:nil_interp+ngh_pp,1:nj_interp,23))

       ! allocate variables that are read
       ! --------------------------------
       ! volume indices
       nx1v= snapshots(nsr)%tectype%nx1
       nx2v= snapshots(nsr)%tectype%nx2
       ny1v= snapshots(nsr)%tectype%ny1
       ny2v= snapshots(nsr)%tectype%ny2
       nz1v= snapshots(nsr)%tectype%nz1
       nz2v= snapshots(nsr)%tectype%nz2
       ! allocate snapshots(nsr)%var
       do nv=1,snapshots(nsr)%nvar
          select case(snapshots(nsr)%var(nv))
          case('rho')
             allocate(rho(nx1v:nx2v,ny1v:ny2v,nz1v:nz2v))
          case('Tmp')
             allocate(Tmp(nx1v:nx2v,ny1v:ny2v,nz1v:nz2v))
          case('prs')
             allocate(prs(nx1v:nx2v,ny1v:ny2v,nz1v:nz2v))
          case('uu')
             allocate(uu(nx1:nx2,ny1:ny2,nz1:nz2))
          case('vv')
             allocate(vv(nx1:nx2,ny1:ny2,nz1:nz2))
          case('ww')
             allocate(ww(nx1:nx2,ny1:ny2,nz1:nz2))
          case('Frhov')
             allocate(Frhov(nx1v:nx2v,ny1v:ny2v,nz1v:nz2v))
          case('Grhow')
             allocate(Grhow(nx1v:nx2v,ny1v:ny2v,nz1v:nz2v))
          end select
       enddo

       ! allocate var for fluctuating u
       if (is_pp_streaks) allocate(uu_fluct(nx1:nx2,ny1:ny2,nz1:nz2))

       ! allocate interpolated var
       if (is_curv) then
          allocate(uut_interp(1-ngh_pp:nil_interp+ngh_pp,1:nj_interp,1-nkernel_k:nk_interp+nkernel_k))
          allocate(uun_interp(1-ngh_pp:nil_interp+ngh_pp,1:nj_interp,1-nkernel_k:nk_interp+nkernel_k))
       else
          allocate(uu_interp(1-ngh_pp:nil_interp+ngh_pp,1:nj_interp,1-nkernel_k:nk_interp+nkernel_k))
          allocate(vv_interp(1-ngh_pp:nil_interp+ngh_pp,1:nj_interp,1-nkernel_k:nk_interp+nkernel_k))
       endif
       allocate(ww_interp(1-ngh_pp:nil_interp+ngh_pp,1:nj_interp,1-nkernel_k:nk_interp+nkernel_k))

       ! allocate y coordinate of the BL edge
       ! allocate(j_edge(1:nil_interp,0:nk_interp+1)) ! ICI MON GARS
       allocate(j_edge(1:nil_interp,ndzt_pp:nfzt_pp))
       allocate(j_edge2(1:nx,1:nz))

       ! allocate laminar-turbulent discriminator arrays
       allocate(extr_density2(1:nx,1:ny,1:nz)) ! Extremums density distribution for computationnal grid
       allocate(extr_density(1-ngh_pp:nil_interp+ngh_pp,1:nj_interp,1-nkernel_k:nk_interp+nkernel_k)) ! for interpolated grid
       allocate(ltbl(1:nx,1:ny,1:nz)) ! 0 -> laminar, 1 -> turbulent
       if (is_curv) then
          allocate(coeff_kc(1:nil_interp,2:nj_interp,-ngh_pp:ngh_pp,2:nj_interp,-nkernel_k:nkernel_k))
          allocate(coeff_kc2(1:nx,2:ny,2:ny))
          allocate(bk_jc2(1:nx,2:ny,2))
       else
          allocate(coeff_kernel(2:nj_interp,-ngh_pp:ngh_pp,2:nj_interp,-nkernel_k:nkernel_k))
          allocate(coeff_kernel2(2:ny,2:ny))
          allocate(born_kern_j2(2:ny,2))
       endif
       allocate(stats_lam(nx,ny,10))
       allocate(stats_turb(nx,ny,10))
       allocate(avg_s_lam(nx,ny,10))
       allocate(avg_s_turb(nx,ny,10))
       if (is_check_stats) allocate(stats_tot(nx,ny,10))
       if (is_check_stats) allocate(avg_s_tot(nx,ny,10))
       allocate(stats_cpt(2,nx,ny))
       stats_cpt = 0.0_wp

    endif

  end subroutine alloc_pp_fstt

end module mod_pp_fstt_main
