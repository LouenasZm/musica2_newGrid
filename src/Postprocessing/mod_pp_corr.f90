!==============================================================================
module mod_pp_corr
!==============================================================================
  !> Post-processing module for spectral evaluation (main)
!==============================================================================
  use warnstop
  use mod_pp_var
  use mod_pp_sp_write
  implicit none
  ! ---------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------

contains

  !============================================================================
  subroutine pp_correl
  !============================================================================
    !> author: AB (directly adapted from pp_sp -> XG)
    !> date: February 2023
    !> Main subroutine for correlation
    !> Only autocorrelation in time and in spanwise direction implemented
    !============================================================================
    use mod_io_snapshots   ! <- for snapshots attribute
    use mod_pp_dfft        ! <- for init_dim_DFFT
    use mod_pp_mpi         ! <- for mpi_types_sp TO BE CHANGED
    use mod_pp_sp_main     ! <- for spectra_param
    use mod_pp_sp_read     ! <- for read plane info
    implicit none
    ! -------------------------------------------------------------------------
    integer :: nbl,nv
    integer :: nx1p,nx2p,ny1p,ny2p,nz1p,nz2p
    ! -------------------------------------------------------------------------

    if (iproc==0) write(*,'(A)') repeat("-",80)

    ! Timestep ~> in info.ini

    ! Define spectra parameters
    ! =========================
    ! Parameters checked here
    if (sp%dim.ne.1) call mpistop('Autocorrelation only implemented in 1D ~> modify sp%dim...',0)
    if ((sp%d(1)%name.ne.'z').and.(sp%d(1)%name.ne.'t')) call mpistop('Autocorrelation only implemented for z and t ~> modify sp%name...',0)
    sp%is_capon = .false.
    !
    call spectra_param
    if (iproc==0) write(*,'(A)') repeat("-",80)

    ! Definition of MPI types for communications
    ! ==========================================
    call mpi_types_sp
    call MPI_BARRIER(COMM_global,info)

    ! Initialize dFFTpack routines
    ! ============================
    call init_dim_DFFT

    ! Allocations
    ! ===========

    ! allocate only variables that are read
    ! -------------------------------------
    ! plane indices
    nx1p= snapshots(nsr)%tectype%nx1
    nx2p= snapshots(nsr)%tectype%nx2
    ny1p= snapshots(nsr)%tectype%ny1
    ny2p= snapshots(nsr)%tectype%ny2
    nz1p= snapshots(nsr)%tectype%nz1
    nz2p= snapshots(nsr)%tectype%nz2
    ! allocate sp%var
    do nv=1,sp%nvar
       select case(sp%varname(nv))
       case('rho')
          allocate(rho(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p))
       case('rhou')
          if (ANY(snapshots(nsr)%var.eq."rhou")) then
             allocate(rhou(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p))
          else if ((ANY(snapshots(nsr)%var.eq."rho")).and.(ANY(snapshots(nsr)%var.eq."uu"))) then
             allocate(uu(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p))
             allocate(rho(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p))
          endif
       case('Tmp')
          allocate(Tmp(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p))
       case('prs')
          allocate(prs(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p))
       case('uu')
          allocate(uu(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p))
       case('vv')
          allocate(vv(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p))
       case('ww')
          allocate(ww(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p))
       case('Frhov')
          allocate(Frhov(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p))
       case('Grhow')
          allocate(Grhow(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p))
       case('uut')    ! u tangential to the wall
          if (.not.is_curv) call mpistop('uut equivalent to uu in cartesian ~> put var to "uu"',0)
          if (bl(1)%BC(3).ne.0) call mpistop('Correl. on tangential velocity for wall not in jmin is not implemented yet !',iproc)
          if (iproc.eq.0) print *,"Correlation on tangential velocity..."
          allocate(uu(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p))
          allocate(vv(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p))
       case('uun')    ! u normal to the wall
          if (.not.is_curv) call mpistop('uun equivalent to uu in cartesian ~> put var to "uu"',0)
          if (bl(1)%BC(3).ne.0) call mpistop('Correl. on normal velocity for wall not in jmin is not implemented yet !',iproc)
          if (iproc.eq.0) print *,"Correlation on normal velocity..."
          allocate(uu(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p))
          allocate(vv(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p))
       end select
    enddo

    ! allocate real array var_r [to read plane]
    ! -------------------------
    allocate(var_r(nl1,nl2,nl3))

    ! allocate complex array var_c [main working array]
    ! ----------------------------
    allocate(var_c(nl1,nl2,nl3))

    ! allocate varm_r [main array]
    ! ---------------
    ! -> block averaging of var_r for autospectra [only allocated if nbloc>1]
    if (sp%nbloc>1) then
       allocate(varm_r(nl1,nl2,nl3))
       varm_r=0.0_wp
    endif

    ! allocate var0_r [store half-block already read]
    ! ---------------
    !if (sp%is_overlap) allocate(var0_r(nl1,nl2,nl3))
    if (sp%is_overlap) allocate(var0_r(nl1,nl2,nl3-sp%loverlap))

    ! **********************
    ! Loop over time blocks  (Welsh's periodogram method)
    ! **********************
    do nbl=1,sp%nbloc
       if (iproc==0) write(*,'(A)') repeat("-",80)

       ! Read planes using MPI-IO routine
       ! ================================
       call read_snapshot_part(nbl)

       ! Initialize complex array
       ! ========================
       var_c=var_r

       ! Compute autocorrelation
       ! =======================
       call compute_autocorr

       ! Averaging
       ! =========
       if (sp%nbloc>1) then
          ! average over blocks
          varm_r = varm_r + real(var_c(:,:,:))
       endif
    enddo
    ! *************************
    ! End loop over time blocks

    ! Average over blocks, normalization and write ouputs
    ! ===================================================
    ! free memory
    if (sp%is_overlap) deallocate(var0_r)
    deallocate(rms2_1var)

    ! autospectra
    if ((sp%dim==1).and.(sp%nbloc==1)) then
       allocate(varm_r(nl1,nl2,nl3))
       varm_r = real(var_c(:,:,:))
    endif

    ! free memory
    deallocate(var_r)

    ! average over blocks
    varm_r=varm_r/sp%nbloc

    ! writing
    call write_corr

    ! Closing sequence
    ! ================

    ! free memory
    deallocate(varm_r)
    deallocate(var_c)

    ! delete local MPI types
    call MPI_TYPE_FREE(type_p,info)
    call MPI_TYPE_FREE(type_pg,info)

    call MPI_BARRIER(COMM_global,info)

  end subroutine pp_correl

  !============================================================================
  subroutine compute_autocorr
  !============================================================================
    !> Compute autospectra (1-D)
    !>  * use real arrays *
  !============================================================================
    use mod_pp_dfft
    use mod_pp_transpose
    use mod_pp_sp_wind
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i1,i2,n
    ! -------------------------------------------------------------------------

    if (iproc==0) then
       print *,'=========================================================='
       print *,' autocorrelation evaluation'
       print *,'=========================================================='
       print *,'~> autocorrelation in '//sp%d(1)%name//'-direction ...'
    endif

    select case (sp%d(1)%name)
    case('z')
       if (ndomz>1) then ! (parallel implementation)

          ! Transposition
          ! -------------
          if (iproc==0) print *,'   - zt-transpose'
          call transpose_zt_c

          ! ! Windowing
          ! ! ---------
          ! call windowing_c(var_ct,1,3)

          ! ZFFT in z (third direction after transposition)
          ! ---------
          if (iproc==0) print *,'   - ZFFT'
          call ZFFT(var_ct,3)
          ! normalization of FFT
          ! --------------------
          var_ct=var_ct/dble(sp%d(1)%lbloc)
          ! compute module
          ! --------------
          var_ct=var_ct*conjg(var_ct)

          ! INVZFFT in z
          ! ---------
          if (iproc==0) print *,'   - INVZFFT'
          call INVZFFT(var_ct,3)

          ! Reverse transposition
          ! ---------------------
          if (iproc==0) print *,'   - zt-untranspose'
          call untranspose_zt_c

       else ! (serial implementation)
          ! ! Windowing
          ! ! ---------
          ! call windowing_c(var_c,1,2)

          ! ZFFT in z
          ! ---------
          if (iproc==0) print *,'   - ZFFT'
          call ZFFT(var_c,sp%d(1)%i)
          ! normalization of FFT
          ! --------------------
          var_c=var_c/dble(sp%d(1)%lbloc)
          ! compute module
          ! --------------
          var_c=var_c*conjg(var_c)

          ! INVZFFT in z
          ! ---------
          if (iproc==0) print *,'   - INVZFFT'
          call INVZFFT(var_c,sp%d(1)%i)

       endif

       ! Normalisation by rms -> To be generalize
       ! --------------------
       do i1=1,nl1
          do n=1,nl3
             var_c(i1,:,n) = var_c(i1,:,n)/rms2_1var(i1,n)
          enddo
       enddo

    case('t')
       ! ! Windowing
       ! ! ---------
       ! call windowing_c(var_c,1,3)

       ! ZFFT in t
       ! ---------
       if (iproc==0) print *,'   - ZFFT'
       call ZFFT(var_c,3)
       ! normalization of FFT
       ! --------------------
       var_c=var_c/dble(sp%d(1)%lbloc)
       ! compute module
       ! --------------
       var_c=var_c*conjg(var_c)

       ! INVZFFT in z
       ! ---------
       if (iproc==0) print *,'   - INVZFFT'
       call INVZFFT(var_c,sp%d(1)%i)

       ! Normalisation by rms -> To be generalize
       ! --------------------
       do i1=1,nl1
          do i2=1,nl2
             var_c(i1,i2,:) = var_c(i1,i2,:)/rms2_1var(i1,i2)
          enddo
       enddo

    end select

  end subroutine compute_autocorr

  !============================================================================
  subroutine write_corr
  !============================================================================
    !> Write correlation
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    ! integer :: i1,i2
    integer :: n
    ! ! file units
    ! integer :: uid,uid_sl(3),uid_m
    ! ! averaging
    ! real(wp), dimension(:,:), allocatable :: var_av
    ! ! slices
    ! real(wp), dimension(:,:,:), allocatable :: E_sl1,E_sl2,E_sl3
    ! -------------------------------------------------------------------------

    ! Open binary file and write header
    ! =================================
    if (iproc==0) call write_corr_header

    ! Allocate averaging array
    ! ========================
    if ((iproc==0).and.(ndir_av>0)) then
       ! allocate
       ! --------
       select case(i_av(1))
       case(1)
          allocate(var_av(ng2,ng3))
       case(2)
          allocate(var_av(ng1,ng3))
       case(3)
          allocate(var_av(ng1,ng2))
       end select

       ! initialize
       ! ----------
       var_av=0.0_wp
    endif

    ! MPI reconstruction & write
    ! ==========================
    do n=1,nl3
       call write_sp_varg(varm_r(:,:,n),n)
    enddo

    ! Proc 0 writes output files
    ! ==========================
    if (iproc==0) then

       if (ndir_av>0) then
          ! End of averaging
          ! ----------------
          var_av=var_av/n_av

          ! Write averaged spectrum
          ! -----------------------
          write(uid) var_av
          deallocate(var_av)
       endif

       ! Close file unit
       ! ---------------
       close(uid)

       ! ! Write slices
       ! ! ------------
       ! if (is_slice) then
       !    ! direction 1
       !    if (n_sl(1)>0) then
       !       do i1=1,n_sl(1)
       !          write(uid_sl(1)) ((E_sl1(i2,n,i1),i2=1,ng2),n=1,ng3)
       !       enddo
       !       close(uid_sl(1))
       !    endif

       !    ! direction 2
       !    if (n_sl(2)>0) then
       !       do i2=1,n_sl(2)
       !          write(uid_sl(2)) ((E_sl2(i1,n,i2),i1=1,ng1),n=1,ng3)
       !       enddo
       !       close(uid_sl(2))
       !    endif

       !    ! direction 3
       !    if (n_sl(3)>0) then
       !       do n=1,n_sl(3)
       !          write(uid_sl(3)) ((E_sl3(i1,i2,n),i1=1,ng1),i2=1,ng2)
       !       enddo
       !       close(uid_sl(3))
       !    endif
       ! endif

    endif

  end subroutine write_corr

  !============================================================================
  subroutine write_corr_header
  !============================================================================
    !> Write header for correlation output files
  !============================================================================
    use mod_grid  ! <- for dimensions
    use mod_io_snapshots
    use mod_utils ! <- for get_free_unit
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i,j,k,nd,cpt
    ! integer :: nlx,nly,nsl
    integer :: nlz,nlt
    real(wp), dimension(:), allocatable :: zg_l,dt_l
    ! character(len=1) :: icap,isl
    character(len=2) :: itype
    character(len=2), dimension(3) :: ispd
    character(len=6) :: isp
    character(len=3) :: inpl
    character(len=10) :: ivar
    character(len=70) :: namefile
    ! character(len=70) :: namefile_sl
    ! -------------------------------------------------------------------------

    ! Counter for uid number
    ! ======================
    cpt=0

    ! Define output file name
    ! =======================

    ! plane number -> character
    ! ------------
    write(inpl,FMT='(I3.3)') nplr

    ! PSD directions
    ! --------------
    ispd=''
    do nd=1,sp%dim
       ispd(nd)=sp%d(nd)%name
    enddo
    isp=trim(ispd(1))//trim(ispd(2))//trim(ispd(3))

    ! variable name
    ! -------------
    select case (trim(sp%varname(1)))
    case('prs') ! pressure
       ivar='pp'
    case( 'uu') ! streamwise velocity
       ivar='uu'
    case( 'vv') ! crossflow velocity
       ivar='vv'
    case( 'ww') ! spanwise velocity
       ivar='ww'
    case('rho') ! density
       ivar='rhorho'
    case('rhou') ! density*uu
       ivar='rhou'
    case('Tmp') ! temperature
       ivar='TT'
    case('Frhov') ! Frhov
       ivar='tauwtauw'
    case('Grhow') ! Grhow
       ivar='GrhowGrhow'
    case('udf') ! user-defined variable
       ivar='udfudf'
    case('uut') ! tangential velocity
       ivar='uut'
    case('uun') ! tangential velocity
       ivar='uun'
    end select

    ! type of snapshot
    ! ----------------
    if (snapshots(nsr)%type.eq.1) then
       itype='_l'
    else if (snapshots(nsr)%type.eq.2) then
       itype='_p'
    endif

    ! file name
    ! ---------
    namefile=trim(dirRESU)//'R'//trim(ivar)//'_'//trim(isp)// &
             trim(itype)//trim(inpl)//'_bl'//trim(numchar(iblc_pp))//'.bin'

    print *,'~> writing correlation file '//trim(namefile)

    ! Open binary file for full spectrum
    ! ==================================
    call get_free_unit(uid)
    cpt=cpt+1
    open(uid,file=trim(namefile),form='unformatted',status='unknown')
    rewind(uid)

    ! write info for non-homogeneous direction
    ! ----------------------------------------
    ! if one direction is inhomogeneous
    if (i_in>0) then
       if (snapshots(nsr)%type.eq.1) then
          select case(name_in)
          case('x')
             write(uid) 1
             write(uid) xg(snapshots(nsr)%ind_i1)
          case('y')
             write(uid) 1
             write(uid) yg(snapshots(nsr)%ind_j1)
          case('z')
             write(uid) 1
             write(uid) zg(snapshots(nsr)%ind_k1)
          case('t')
             write(uid) ngt
             write(uid) dt_spec
          end select
       else if (snapshots(nsr)%type.eq.2) then
          print *,'   - inhomogeneous direction ',name_in
          select case(name_in)
          case('x')
             write(uid) ngx
             write(uid) (xg(i),i=1,ngx)
          case('y')
             write(uid) ngy
             write(uid) (yg(j),j=1,ngy)
          case('z')
             write(uid) ngz
             write(uid) (zg(k),k=1,ngz)
          case('t')
             write(uid) ngt
             write(uid) dt_spec
          end select
       endif
    endif

    ! write info for spectrum directions: only positive frequencies/wavenumbers
    ! ----------------------------------
    if (iproc==0) print *,'   - correlation direction ',sp%d(1)%name
    select case(sp%d(1)%name)
    case('z')
       nlz=ngz/2
       write(uid) nlz
       allocate(zg_l(nlz))
       do k=1,nlz
          zg_l(k)=dble(k-1)*deltaz;
       enddo
       write(uid) (zg_l(k),k=1,nlz)
       deallocate(zg_l)
    case('t')
       nlt=nl3/2
       write(uid) nlt
       allocate(dt_l(nlt))
       do k=1,nlt
          dt_l(k)=dble(k-1)*dt_spec;
       enddo
       write(uid) (dt_l(k),k=1,nlt)
       deallocate(dt_l)
    end select

    if (ndir_av>0) print *,'   - averaging direction ',dir_av(1)

    ! ! Open binary file for slice outputs
    ! ! ==================================
    ! if (is_slice) then

    !    do nsl=1,3 ! * loop on slices *

    !       if (n_sl(nsl)>0) then

    !          ! slice number
    !          ! ------------
    !          write(isl,FMT='(I1.1)') nsl

    !          ! file name
    !          ! ---------
    !          if (i_in>0) then
    !             namefile_sl=trim(dirRESU)//'E'//trim(ivar)//'_'//trim(isp)//'_slices'//trim(isl) &
    !                         //'_'//name_in//'_p'//trim(inpl)//'_'//trim(icap)//'.bin'
    !          else
    !             namefile_sl=trim(dirRESU)//'E'//trim(ivar)//'_'//trim(isp)//'_slices'//trim(isl) &
    !                         //'_p'//trim(inpl)//'_'//trim(icap)//'.bin'
    !          endif

    !          ! file unit & open binary
    !          ! -----------------------
    !          uid_sl(nsl)=uid+cpt
    !          cpt=cpt+1
    !          open(uid_sl(nsl),file=trim(namefile_sl),form='unformatted',status='unknown')
    !          rewind(uid_sl(nsl))

    !          print *,'~> writing spectrum slices '//trim(namefile_sl)

    !          ! write info for spectrum directions
    !          ! ----------------------------------
    !          do nd=1,sp%dim
    !             select case(sp%d(nd)%name)
    !             case('x')
    !                write(uid_sl(nsl)) ngx
    !                write(uid_sl(nsl)) (kx(k),k=-nlx+1,nlx)
    !             case('y')
    !                write(uid_sl(nsl)) ngy
    !                write(uid_sl(nsl)) (ky(k),k=-nly+1,nly)
    !             case('z')
    !                write(uid_sl(nsl)) ngz
    !                write(uid_sl(nsl)) (kz(k),k=-nlz+1,nlz)
    !             case('t')
    !                write(uid_sl(nsl)) nl3
    !                write(uid_sl(nsl)) (ff(k),k=-nlt+1,nlt)
    !             end select
    !          enddo

    !          ! allocate slice array
    !          ! --------------------
    !          select case(nsl)
    !          case(1)
    !             allocate(E_sl1(ng2,ng3,n_sl(nsl)))
    !          case(2)
    !             allocate(E_sl2(ng1,ng3,n_sl(nsl)))
    !          case(3)
    !             allocate(E_sl3(ng1,ng2,n_sl(nsl)))
    !          end select

    !       endif

    !    enddo ! * end loop on slices *

    ! endif

    if (allocated(zg_l)) deallocate(zg_l)
    if (allocated(dt_l)) deallocate(dt_l)

  end subroutine write_corr_header

end module mod_pp_corr
