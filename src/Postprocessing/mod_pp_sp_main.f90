!==============================================================================
module mod_pp_sp_main
!==============================================================================
  !> Post-processing module for spectral evaluation (main)
!==============================================================================
  use warnstop
  use mod_pp_var
  implicit none
  ! ---------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------

contains

  !============================================================================
  subroutine pp_spectra
  !============================================================================
    !> author: XG 
    !> date: March 2022
    !> Main subroutine for spectrum analysis
    !============================================================================
    use mod_time         ! <- for deltat,dtstar
    use mod_io_snapshots ! <- for snapshots attribute
    use mod_eigenmode    ! <- for to recompute time step
    use mod_pp_dfft      ! <- for init_dim_DFFT
    use mod_pp_mpi       ! <- for mpi_types_sp TO BE CHANGED
    use mod_pp_sp_read   ! <- for read plane info
    use mod_pp_sp_eval   ! <- for PSD estimation (1D, 2D or 3D)
    use mod_pp_sp_modal  ! <- for modal analysis
    use mod_pp_sp_write  ! <- for write spectra
    !use mod_pp_kw_spectrum
    implicit none
    ! -------------------------------------------------------------------------
    integer :: nbl,nv
    integer :: nx1p,nx2p,ny1p,ny2p,nz1p,nz2p
    real(wp) :: Tround
    ! -------------------------------------------------------------------------
    integer :: i,j,k,ngx_,ngy_,ngz_

    ! Read parameters of old simulations
    ! ==================================
    !call read_param_old
    if (iproc==0) write(*,'(A)') repeat("-",80)

    ! Time step computation
    ! =====================
    !!call timestep

    if (CHAN) then
       !!deltat = 1.581980679150934E-009
       !!dtstar = deltat/tscale
       !!dtstar = 9.890879081397594E-005
       deltat = 4.827069251767853E-009
       dtstar = 5.543239924738873E-005
       !!if (iproc.eq.0) write(*,*) 'dyg(1)',dyg(1)
       if (iproc.eq.0) write(6,*) 'dyg(1)',idy(1)
       if (iproc.eq.0) write(6,*) 'Dt [s], Dt* [-]:',deltat,dtstar

       open(51,file=trim(dirDATA)//'vites_u.bin',form='unformatted',status='unknown')
       rewind(51)
       ! global dim
       ! ==========
       read(51) ngx_
       read(51) ngy_
       read(51) ngz_
       ! global grid
       ! ===========
       allocate(xg(ngx),yg(ngy),zg(ngz))
       read(51) (xg(i),i=1,ngx)
       read(51) (yg(j),j=1,ngy)
       read(51) (zg(k),k=1,ngz)
       close(51)

    endif

    if (STBL.and.is_eigenmode) then
       if (coord(1)==0) then
          ! read unstable mode parameters to define eigenmodes
          ! --------------------------------------------------
          call read_param_stab
          ! dimensionalize
          ! --------------
          do nv=1,n_eig
             eig(nv)%om=eig(nv)%om/L_ref*u_ref
             eig(nv)%kr=eig(nv)%kr/L_ref
             eig(nv)%ki=eig(nv)%ki/L_ref
             eig(nv)%beta=eig(nv)%beta/L_ref
          enddo
          ! adjust time step
          ! ----------------
          if (iproc.eq.0) write(6,*) '  T/dt            =', twopi/eig(1)%om/deltat
          Tround=anint(twopi/eig(1)%om/deltat/10.0_wp)*10.0_wp
          deltat=twopi/eig(1)%om/Tround
          dtstar = deltat/tscale
          if (iproc.eq.0) write(6,*) 'round(T/dt)       =', Tround
          if (iproc.eq.0) write(6,*) 'new Dt [s], Dt* [-]:', deltat, dtstar
       endif
       call MPI_BCAST(deltat,1,MPI_DOUBLE_PRECISION,0,COMM_global,info)
       call MPI_BCAST(dtstar,1,MPI_DOUBLE_PRECISION,0,COMM_global,info)
    endif
    if (iproc==0) write(*,'(A)') repeat("-",80)

    ! Determine modes for modal analysis
    ! ==================================
    if (is_modal) then
       call comp_modes
       if (iproc==0) write(*,'(A)') repeat("-",80)
    endif
    
    ! Define spectra parameters
    ! =========================
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
       case('Tmp')  
          allocate(Tmp(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p))
       case('prs')  
          allocate(prs(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p))
       case('rhou')
          if (ANY(snapshots(nsr)%var.eq."rhou")) then
             allocate(rhou(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p))
          else if ((ANY(snapshots(nsr)%var.eq."rho")).and.(ANY(snapshots(nsr)%var.eq."uu"))) then
             allocate(uu(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p))
             allocate(rho(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p))
          endif
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
          if (bl(1)%BC(3).ne.0) call mpistop('Spectra on tangential velocity for wall not in jmin is not implemented yet !',iproc)
          if (iproc.eq.0) print *,"FFT on tangential velocity..."
          ! if (.not.allocated(uu)) allocate(uu(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p))
          ! if (.not.allocated(vv)) allocate(vv(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p))
          allocate(uu(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p))
          allocate(vv(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p))
       case('uun')    ! u tangential to the wall
          if (.not.is_curv) call mpistop('uun equivalent to uu in cartesian ~> put var to "uu"',0)
          if (bl(1)%BC(3).ne.0) call mpistop('Spectra on normal velocity for wall not in jmin is not implemented yet !',iproc)
          if (iproc.eq.0) print *,"FFT on normal velocity..."
          ! if (.not.allocated(uu)) allocate(uu(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p))
          ! if (.not.allocated(vv)) allocate(vv(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p))
          allocate(uu(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p))
          allocate(vv(nx1p:nx2p,ny1p:ny2p,nz1p:nz2p))
       end select
    enddo

    ! allocate real array var_r [main array for autospectra / used for rms for multiD-spectra]
    ! -------------------------
    allocate(var_r(nl1,nl2,nl3))
    
    ! allocate complex array var_c for multiD-spectra
    ! -----------------------------------------------
    if (sp%dim>1) allocate(var_c(nl1,nl2,nl3))
    
    ! allocate varm_r [main array for autospectra / store abs for multiD-spectra]
    ! ---------------
    ! -> block averaging of var_r for autospectra [only allocated if nbloc>1]
    ! -> block averaging of module for multiD-spectra
    if ((sp%dim>1).or.(sp%nbloc>1)) then
       allocate(varm_r(nl1,nl2,nl3))
       varm_r=0.0_wp
    endif

    ! allocate var0_r [store half-block already read]
    ! ---------------
    !if (sp%is_overlap) allocate(var0_r(nl1,nl2,nl3))
    if (sp%is_overlap) allocate(var0_r(nl1,nl2,nl3-sp%loverlap))
   
    ! allocate rms squared, skewness and kurtosis arrays
    ! --------------------------------------------------
    ! if one direction is inhomogeneous
    if (i_in>0) then
       allocate(rms2_in(n_in),skew_in(n_in),kurt_in(n_in))
       ! averaging over blocks
       if (sp%nbloc>1) then
          allocate(rms2m_in(n_in),skewm_in(n_in),kurtm_in(n_in))
          rms2m_in=0.0_wp
          skewm_in=0.0_wp
          kurtm_in=0.0_wp
       endif
    endif

    if (sp%nbloc==0) call mpistop('no block: check parameters !',0)

    ! **********************
    ! Loop over time blocks  (Welsh's periodogram method)
    ! **********************
    do nbl=1,sp%nbloc
       if (iproc==0) write(*,'(A)') repeat("-",80)

       ! Read planes using MPI-IO routine
       ! ================================
       call read_snapshot_part(nbl)

       if (sp%dim==1) then
          ! Compute autospectra
          ! ===================
          call compute_autospectrum
          if (sp%nbloc>1) then
             ! average over blocks
             varm_r=varm_r+var_r
          endif
       else 
          ! Compute multi-D spectra
          ! =======================
          call compute_multiD_spectrum
          ! module average over blocks
          if (sp%is_capon) then
             varm_r=varm_r+abs(var_c)
          else
             varm_r=varm_r+var_c*conjg(var_c)
          endif
       endif
       
       ! average over blocks
       ! -------------------
       if ((i_in>0).and.(sp%nbloc>1)) then
          ! if one direction is inhomogeneous
          rms2m_in=rms2m_in+rms2_in
          skewm_in=skewm_in+skew_in
          kurtm_in=kurtm_in+kurt_in
       else
          rms2m=rms2m+rms2
       endif
       
    enddo
    ! *************************
    ! End loop over time blocks
    ! *************************

    ! Average over blocks, normalization and write ouputs
    ! ===================================================
    ! free memory
    if (sp%is_overlap) deallocate(var0_r)
    if (sp%dim==1) deallocate(rms2_1var)
    
    ! autospectra 
    if ((sp%dim==1).and.(sp%nbloc==1)) then      
       allocate(varm_r(nl1,nl2,nl3))
       varm_r=var_r
    endif
      
    ! if one direction is inhomogeneous
    if (i_in>0) then
       if (sp%nbloc==1) then
          allocate(rms2m_in(n_in),skewm_in(n_in),kurtm_in(n_in))
          rms2m_in=rms2_in
          skewm_in=skew_in
          kurtm_in=kurt_in
       else
          ! average over blocks
          rms2m_in=rms2m_in/sp%nbloc
          skewm_in=skewm_in/sp%nbloc
          kurtm_in=kurtm_in/sp%nbloc
       endif
    else
       ! average over blocks
       rms2m=rms2m/sp%nbloc
    endif
    ! free memory
    if (i_in>0) deallocate(rms2_in,skew_in,kurt_in)
    deallocate(var_r)
   
    ! average over blocks
    varm_r=varm_r/sp%nbloc

    ! normalization
    call spectrum_normalization

    ! writing
    call write_spectrum

    ! Closing sequence
    ! ================
    
    ! free memory
    if (i_in>0) deallocate(rms2m_in,skewm_in,kurtm_in)
    deallocate(varm_r)
    if (sp%dim>1) deallocate(var_c)
    
    ! delete local MPI types
    call MPI_TYPE_FREE(type_p,info)
    call MPI_TYPE_FREE(type_pg,info)
    if (i_in>0) then
       call MPI_TYPE_FREE(type_in,info)
       call MPI_TYPE_FREE(typeg_in,info)
    endif
    
    call MPI_BARRIER(COMM_global,info)
    
  end subroutine pp_spectra

  ! !============================================================================
  ! subroutine read_param_pp_new
  ! !============================================================================
  !   !> Determine post-processing parameters
  ! !============================================================================
  !   use mod_block
  !   implicit none
  !   ! -------------------------------------------------------------------------
  !   integer :: lbloc
  !   logical :: ispwall=.false.    ! k-w spectra wall pressure
  !   logical :: ispwall_autosp     ! autospectra wall pressure
  !   ! -------------------------------------------------------------------------

  !   ! Directories
  !   ! ===========
  !   dirRESU='./'

  !   ! Read post-processing parameters
  !   ! ===============================
  !   open(30,file='param_pp.ini')
  !   read(30,*) ! Post-processing parameters
  !   read(30,*) ! ==========================
  !   read(30,*) ! number of selected classes: nclass
  !   read(30,*) nclass
  !   read(30,*) ! size of temporal blocks: lbloc
  !   read(30,*) lbloc
  !   read(30,*) ! overlapping (Welch's method): is_overlap / loverlap
  !   read(30,*) sp%is_overlap, sp%loverlap
  !   read(30,*) ! Capon's method: is_capon
  !   read(30,*) sp%is_capon
  !   read(30,*) ! Capon filter length (Capon's order)
  !   read(30,*) sp%ncapon
  !   read(30,*) ! number of MPI domains per direction: ndomx,ndomy,ndomz
  !   !read(30,*) ndomx
  !   !read(30,*) ndomy
  !   !read(30,*) ndomz
  !   read(30,*) bl(1)%ndomi
  !   read(30,*) bl(1)%ndomj
  !   read(30,*) bl(1)%ndomk
  !   bl(1)%nproc=bl(1)%ndomi*bl(1)%ndomj*bl(1)%ndomk
  !   read(30,*) ! plane number: nplr
  !   read(30,*) nplr
  !   read(30,*) ! directory: dirDATA
  !   read(30,*) dirDATA
  !   read(30,*) ! initial time for samples
  !   read(30,*) time_ini
  !   read(30,*) ! Choice of actions
  !   read(30,*) ! -----------------
  !   read(30,*) ! logical ispwall
  !   read(30,*) ispwall
  !   read(30,*) ! logical ispwall_autosp / compute autospectra
  !   read(30,*) ispwall_autosp
  !   close(30)

  !   if (ispwall_autosp) sp%dim=1
  !   if (ispwall) sp%dim=3
    
  !   allocate(sp%d(sp%dim))
    
  !   sp%d(1)%lbloc=lbloc
  !   sp%d(1)%type_win='T'
  !   sp%d(1)%param_win=1.0_wp
  !   sp%is_onesided=.true.

  ! end subroutine read_param_pp_new

  !============================================================================
  subroutine spectra_param
  !============================================================================
    !> Determine & print spectra parameters
  !============================================================================
    use mod_io_snapshots
    implicit none
    ! -------------------------------------------------------------------------
    integer :: nd,nd_t,ntt
    integer :: nm1,nm2,nm3,nav
    integer(8) :: ng_tot
    real(wp) :: dom,dkx,dky,dkz ! frequency and wavenumber resolutions  
    ! -------------------------------------------------------------------------

    if ((snapshots(nsr)%type.ne.2).and.(snapshots(nsr)%type.ne.1)) then
       call mpistop('Spectra post-processing only available for planes & lines', 0)
    endif

    
    ! Time dimension
    ! ==============
    ! sampling time step
    dt_spec=deltat*snapshots(nsr)%freq
    ! number of classes
    nclass=1
    ! number of samples
    ntt=nmax/snapshots(nsr)%freq
    ! total number of samples
    ngt=nclass*ntt


    ! set lbloc to full dimension if 0 in param_pp.ini
    ! ================================================
    do nd=1,sp%dim
       if (sp%d(nd)%lbloc==0) then
          if (sp%d(nd)%name=='x') sp%d(nd)%lbloc=ngx
          if (sp%d(nd)%name=='y') sp%d(nd)%lbloc=ngy
          if (sp%d(nd)%name=='z') sp%d(nd)%lbloc=ngz
          if (sp%d(nd)%name=='t') sp%d(nd)%lbloc=ngt
       endif
    enddo
    
    ! Check parameters
    ! ================

    ! 1/ Check number of procs #1: only 1 proc in direction not in plane
    ! ----------------------------
    if (snapshots(nsr)%type.eq.2) then
       if ((snapshots(nsr)%normal==1).and.(ndomx.ne.1)) &
            call mpistop('post-pro of yz-plane: ndomx must be 1',0)
       if ((snapshots(nsr)%normal==2).and.(ndomy.ne.1)) &
            call mpistop('post-pro of xz-plane: ndomy must be 1',0)
       if ((snapshots(nsr)%normal==3).and.(ndomz.ne.1)) &
            call mpistop('post-pro of xy-plane: ndomz must be 1',0)
    else if (snapshots(nsr)%type.eq.1) then
       if ((snapshots(nsr)%dir==1).and.((ndomy.ne.1).or.(ndomz.ne.1))) &
            call mpistop('post-pro of x-line: ndomy and ndomz must be 1',0)
       if ((snapshots(nsr)%dir==2).and.((ndomx.ne.1).or.(ndomz.ne.1))) &
            call mpistop('post-pro of y-line: ndomx and ndomz must be 1',0)
       if ((snapshots(nsr)%dir==3).and.((ndomx.ne.1).or.(ndomy.ne.1))) &
            call mpistop('post-pro of z-line: ndomx and ndomy must be 1',0)
    endif


    ! 2/ Modify eventually parameters for homogeneous directions
    ! ----------------------------------------------------------
    do nd=1,sp%dim
       if ( ((sp%d(nd)%name=='x').and.(periods(1))).or. &
            ((sp%d(nd)%name=='y').and.(periods(2))).or. &
            ((sp%d(nd)%name=='z').and.(periods(3))) ) then
          ! for homogeneous directions -> change to no window
          ! -------------------------------------------------
          if (sp%d(nd)%type_win=='K') then
             if (iproc==0) then
                print *,'%%%%%%%%%%%%%%%%%%%%%%%%%%'
                print *,sp%d(nd)%name//'-dimension is homogeneous'
                print *,'%%%%%%%%%%%%%%%%%%%%%%%%%%'
                print *,'--> windowing parameters set to NO WINDOW'
             endif
             sp%d(nd)%type_win='T'
             sp%d(nd)%param_win=0.0_wp
          endif
          if ((sp%d(nd)%type_win=='T').and.(sp%d(nd)%param_win.ne.0.0_wp)) then
             if (iproc==0) then
                print *,'%%%%%%%%%%%%%%%%%%%%%%%%%%'
                print *,sp%d(nd)%name//'-dimension is homogeneous'
                print *,'%%%%%%%%%%%%%%%%%%%%%%%%%%'
                print *,'--> windowing parameters set to NO WINDOW'
             endif
             sp%d(nd)%param_win=0.0_wp
          endif
!!$          ! for homogeneous directions -> not Capon a priori
!!$          ! ------------------------------------------------
!!$          if ((sp%is_capon).and.(nd==sp%dim)) then
!!$             if (iproc==0) then
!!$                print *,'%%%%%%%%%%%%%%%%%%%%%%%%%%'
!!$                print *,sp%d(nd)%name//'-dimension is homogeneous'
!!$                print *,'%%%%%%%%%%%%%%%%%%%%%%%%%%'
!!$                print *,'--> use FOURIER instead of CAPON for that direction'
!!$             endif
!!$             sp%is_capon=.false.
!!$          endif
          ! for homogeneous directions -> direction should not be cut into segments
          ! -----------------------------------------------------------------------
          if ( ((sp%d(nd)%name=='x').and.(sp%d(nd)%lbloc.ne.ngx)).or. &
               ((sp%d(nd)%name=='y').and.(sp%d(nd)%lbloc.ne.ngy)).or. &
               ((sp%d(nd)%name=='z').and.(sp%d(nd)%lbloc.ne.ngz)) ) then
             if (iproc==0) then
                print *,'%%%%%%%%%%%%%%%%%%%%%%%%%%'
                print *,sp%d(nd)%name//'-dimension is homogeneous'
                print *,'%%%%%%%%%%%%%%%%%%%%%%%%%%'
                print *,'--> direction should not be cut into segments'
                print *,'--> new block length is ng'//sp%d(nd)%name
             endif
              if (sp%d(nd)%name=='x') sp%d(nd)%lbloc=ngx
              if (sp%d(nd)%name=='y') sp%d(nd)%lbloc=ngy
              if (sp%d(nd)%name=='z') sp%d(nd)%lbloc=ngz
          endif
       endif
    enddo

    ! 3/ Welch's method only implemented for time / Determine number of blocks
    ! ------------------------------------------------------------------------
    ! check is_overlap
    if (sp%is_overlap) then
       if ((sp%dim==1).and.(sp%d(1)%name.ne.'t')) call mpistop('overlap Welsh method only for time dimension. Please change !',0)
       if (sp%dim==2) then
          if ((sp%d(1)%name.ne.'t').and.(sp%d(2)%name.ne.'t')) call mpistop('overlap Welsh method only for time dimension. Please change !',0)
       endif
    endif

    ! check and determine number of blocks
    do nd=1,sp%dim
       ! number of blocks
       select case (sp%d(nd)%name)
       case('x')          
          sp%nbloc=ngx/sp%d(nd)%lbloc
          if (sp%nbloc.ne.1) call mpistop('Welsh method only for t-dim. Change block length for x-dim!',0)
       case('y')          
          sp%nbloc=ngy/sp%d(nd)%lbloc
          if (sp%nbloc.ne.1) call mpistop('Welsh method only for t-dim. Change block length for y-dim!',0)
       case('z')          
          sp%nbloc=ngz/sp%d(nd)%lbloc
          if (sp%nbloc.ne.1) call mpistop('Welsh method only for t-dim. Change block length for z-dim!',0)
       case('t')          
          sp%nbloc=ngt/sp%d(nd)%lbloc
          ! time segment overlap (for Welch's method)
          if (sp%is_overlap) sp%nbloc=(ngt-sp%d(nd)%lbloc)/sp%loverlap+1
          ! if not Welch for time, check dimension
          if ((sp%nbloc==1).and.(sp%d(nd)%lbloc.gt.ngt)) &
               call mpistop('Problem in time block length: greater than number of samples!',0)   
       end select
    enddo

    ! 4/ init window factors
    ! ----------------------
    do nd=1,sp%dim
       sp%d(nd)%Cw=1.0_wp
    enddo

    ! Determine spectral resolution in each directions
    ! ================================================

    ! initialization at zero (if used and not defined -> divide by zero)
    ! ----------------------
    dom=0.0_wp
    dkx=0.0_wp
    dky=0.0_wp
    dkz=0.0_wp
    
    deltax=x(2)-x(1)
    deltaz=z(2)-z(1)
    
    do nd=1,sp%dim
       select case (sp%d(nd)%name)
       case('x')          
          ! wavenumber increment
          ! ---------------------
          !dkx=2.*pi/((ngx-1)*deltax)
          print *,sp%d(nd)%lbloc,'deltax',deltax
          dkx=twopi/(sp%d(nd)%lbloc*deltax)
       case('y')          
          ! wavenumber increment
          ! ---------------------
          dky=twopi/(sp%d(nd)%lbloc*deltay)         
       case('z')                    
          ! wavenumber increment
          ! ---------------------
          !dkz=2.*pi/((ngz-1)*deltaz)
          dkz=twopi/(sp%d(nd)%lbloc*deltaz)
       case('t')                       
          ! frequency resolution
          ! --------------------
          ! angular freq.
          dom=twopi/(sp%d(nd)%lbloc*dt_spec)
       end select
    enddo

    ! Determine dimension per direction
    ! =================================
    allocate(coordx1(0:nproc-1),coordx2(0:nproc-1))

    ! initialize n_av
    ! ---------------
    n_av=1 ! default value if no averaging
    ndom_sp=1 ! default value



    ! Determine dimension per direction for YZ,XZ or XY planes
    ! ========================================================
    if (snapshots(nsr)%type.eq.2) then

       select case(snapshots(nsr)%normal)

       ! plane YZ (with normal 1 along X)
       ! ================================
       case (1)
          ! PSD directions
          ! --------------
          do nd=1,sp%dim
             select case(sp%d(nd)%name)
             case ('y')
                sp%d(nd)%i=1
             case ('z')
                sp%d(nd)%i=2
             case ('t')
                sp%d(nd)%i=3
                nd_t=nd
             case default
                call mpistop('bad choice for spectra direction: yz-plane -> choose y,z or t', 0)
             end select
          enddo

          ! MPI cutting of spectral direction (for 1D spectra)
          ! ---------------------------------
          select case(sp%d(1)%name)
          case ('y')
             ndom_sp=ndomy
             i_sp=2
          case ('z')
             ndom_sp=ndomz
             i_sp=3
          end select

          ! coordxi in MPI_CART
          ! -------------------
          coordx1=coordy
          coordx2=coordz

          ! ndom
          ! ----
          ndom1=ndomy
          ndom2=ndomz

          ! dim per direction
          ! -----------------
          nl1=ny
          nl2=nz
          if (sp%nbloc>1) then
             nl3=sp%d(nd_t)%lbloc
          else
             nl3=ngt
          endif
          ng1=ngy
          ng2=ngz
          ng3=nl3

          ! find dim for averaging and inhomogeneous direction
          ! --------------------------------------------------
          nm1=1
          nm2=1
          nm3=1
          i_in=0
          nav=0
          do nd=1,sp%dim
             if (sp%d(nd)%name=='y') then
                nm1=nl1
                nav=nav+1
             endif
          enddo
          do nd=1,ndir_av
             if (dir_av(nd)=='y') then
                i_av(nd)=1
                n_av=ngy
                nm1=nl1
                nav=nav+1
             endif
          enddo
          if (nm1==1) then
             i_in=1
             n_in=nl1
             name_in='y'
             ng_in=ngy
          endif
          do nd=1,sp%dim
             if (sp%d(nd)%name=='z') then
                nm2=nl2
                nav=nav+1
             endif
          enddo
          do nd=1,ndir_av
             if (dir_av(nd)=='z') then
                i_av(nd)=2
                n_av=ngz
                nm2=nl2
                nav=nav+1
             endif
          enddo
          if (nm2==1) then
             i_in=2
             n_in=nl2
             name_in='z'
             ng_in=ngz
          endif
          do nd=1,sp%dim
             if (sp%d(nd)%name=='t') then
                nm3=nl3
                nav=nav+1
             endif
          enddo
          do nd=1,ndir_av
             if (dir_av(nd)=='t') then
                i_av(nd)=3
                n_av=ngt
                nm3=nl3
                nav=nav+1
             endif
          enddo
          if (nm3==1) then
             i_in=3
             n_in=nl3
             name_in='t'
             ng_in=ngt
          endif

       ! plane XZ (with normal 2 along Y)
       ! ================================
       case (2)
          ! PSD directions
          ! --------------
          do nd=1,sp%dim
             select case(sp%d(nd)%name)
             case ('x')
                sp%d(nd)%i=1
             case ('z')
                sp%d(nd)%i=2
             case ('t')
                sp%d(nd)%i=3
                nd_t=nd
             case default
                call mpistop('bad choice for spectra direction: xz-plane -> choose x,z or t', 0)
             end select
          enddo

          ! MPI cutting of spectral direction (for 1D spectra)
          ! ---------------------------------
          select case(sp%d(1)%name)
          case ('x')
             ndom_sp=ndomx
             i_sp=1
          case ('z')
             ndom_sp=ndomz
             i_sp=3
          end select

          ! coordxi in MPI_CART
          ! -------------------
          coordx1=coordx
          coordx2=coordz

          ! ndom
          ! ----
          ndom1=ndomx
          ndom2=ndomz

          ! dim per direction
          ! -----------------
          nl1=nx
          nl2=nz
          if (sp%nbloc>1) then
             nl3=sp%d(nd_t)%lbloc
          else
             nl3=ngt
          endif
          ng1=ngx
          ng2=ngz
          ng3=nl3

          ! find dim for averaging and inhomogeneous direction
          ! --------------------------------------------------
          nm1=1
          nm2=1
          nm3=1
          i_in=0
          nav=0
          do nd=1,sp%dim
             if (sp%d(nd)%name=='x') then
                nm1=nl1
                nav=nav+1
             endif
          enddo
          do nd=1,ndir_av
             if (dir_av(nd)=='x') then
                i_av(nd)=1
                n_av=ngx
                nm1=nl1
                nav=nav+1
             endif
          enddo
          if (nm1==1) then
             i_in=1
             n_in=nl1
             name_in='x'
             ng_in=ngx
          endif
          do nd=1,sp%dim
             if (sp%d(nd)%name=='z') then
                nm2=nl2
                nav=nav+1
             endif
          enddo
          do nd=1,ndir_av
             if (dir_av(nd)=='z') then
                i_av(nd)=2
                n_av=ngz
                nm2=nl2
                nav=nav+1
             endif
          enddo
          if (nm2==1) then
             i_in=2
             n_in=nl2
             name_in='z'
             ng_in=ngz
          endif
          do nd=1,sp%dim
             if (sp%d(nd)%name=='t') then
                nm3=nl3
                nav=nav+1
             endif
          enddo
          do nd=1,ndir_av
             if (dir_av(nd)=='t') then
                i_av(nd)=3
                n_av=ngt
                nm3=nl3
                nav=nav+1
             endif
          enddo
          if (nm3==1) then
             i_in=3
             n_in=nl3
             name_in='t'
             ng_in=ngt
          endif

       ! plane XY (with normal 3 along Z)
       ! ================================
       case (3)
          ! PSD directions
          ! --------------
          do nd=1,sp%dim
             select case(sp%d(nd)%name)
             case ('x')
                sp%d(nd)%i=1
                if (nd==1) ndom_sp=ndomx
             case ('y')
                sp%d(nd)%i=2
                if (nd==1) ndom_sp=ndomy
             case ('t')
                sp%d(nd)%i=3
                nd_t=nd
             case default
                call mpistop('bad choice for spectra direction: yz-plane -> choose x,y or t', 0)
             end select
          enddo

          ! MPI cutting of spectral direction (for 1D spectra)
          ! ---------------------------------
          select case(sp%d(1)%name)
          case ('x')
             ndom_sp=ndomx
             i_sp=1
          case ('y')
             ndom_sp=ndomy
             i_sp=2
          end select

          ! coordxi in MPI_CART
          ! -------------------
          coordx1=coordx
          coordx2=coordy

          ! ndom
          ! ----
          ndom1=ndomx
          ndom2=ndomy

          ! dim per direction
          ! -----------------
          nl1=nx
          nl2=ny
          if (sp%nbloc>1) then
             nl3=sp%d(nd_t)%lbloc
          else
             nl3=ngt
          endif
          ng1=ngx
          ng2=ngy
          ng3=nl3

          ! find dim for averaging and inhomogeneous direction
          ! --------------------------------------------------
          nm1=1
          nm2=1
          nm3=1
          i_in=0
          nav=0
          do nd=1,sp%dim
             if (sp%d(nd)%name=='x') then
                nm1=nl1
                nav=nav+1
             endif
          enddo
          do nd=1,ndir_av
             if (dir_av(nd)=='x') then
                i_av(nd)=1
                n_av=ngx
                nm1=nl1
                nav=nav+1
             endif
          enddo
          if (nm1==1) then
             i_in=1
             n_in=nl1
             name_in='x'
             ng_in=ngx
          endif
          do nd=1,sp%dim
             if (sp%d(nd)%name=='y') then
                nm2=nl2
                nav=nav+1
             endif
          enddo
          do nd=1,ndir_av
             if (dir_av(nd)=='y') then
                i_av(nd)=2
                n_av=ngy
                nm2=nl2
                nav=nav+1
             endif
          enddo
          if (nm2==1) then
             i_in=2
             n_in=nl2
             name_in='y'
             ng_in=ngy
          endif
          do nd=1,sp%dim
             if (sp%d(nd)%name=='t') then
                nm3=nl3
                nav=nav+1
             endif
          enddo
          do nd=1,ndir_av
             if (dir_av(nd)=='t') then
                i_av(nd)=3
                n_av=ngt
                nm3=nl3
                nav=nav+1
             endif
          enddo
          if (nm3==1) then
             i_in=3
             n_in=nl3
             name_in='t'
             ng_in=ngt
          endif

       end select

    ! Determine dimension per direction for X,Y or Z lines
    ! ====================================================
    else if (snapshots(nsr)%type.eq.1) then
       select case(snapshots(nsr)%dir)

       ! line X (normal to plane YZ)
       ! ===========================
       case (1)
          ! PSD directions
          ! --------------
          do nd=1,sp%dim
             select case(sp%d(nd)%name)
             case ('x')
                sp%d(nd)%i=1
             case ('t')
                sp%d(nd)%i=3
                nd_t=nd
             case default
                call mpistop('bad choice for spectra direction: x-line -> choose x or t', 0)
             end select
          enddo

          ! MPI cutting of spectral direction (for 1D spectra)
          ! ---------------------------------
          select case(sp%d(1)%name)
          case ('x')
             ndom_sp=ndomx
             i_sp=1
          end select

          ! coordxi in MPI_CART
          ! -------------------
          coordx1=coordx
          coordx2=coordz

          ! ndom
          ! ----
          ndom1=ndomx
          ndom2=1

          ! dim per direction
          ! -----------------
          nl1=nx
          nl2=1
          if (sp%nbloc>1) then
             nl3=sp%d(nd_t)%lbloc
          else
             nl3=ngt
          endif
          ng1=ngx
          ng2=1
          ng3=nl3

          ! find dim for averaging and inhomogeneous direction
          ! --------------------------------------------------
          nm1=1
          nm2=1
          nm3=1
          i_in=0
          nav=0
          do nd=1,sp%dim
             if (sp%d(nd)%name=='x') then
                nm1=nl1
                nav=nav+1
             endif
          enddo
          do nd=1,ndir_av
             if (dir_av(nd)=='x') then
                i_av(nd)=1
                n_av=ngx
                nm1=nl1
                nav=nav+1
             endif
          enddo
          ! Inhomogeneous direction superficialy taken a y
          i_in=2
          n_in=nl2
          name_in='y'
          ng_in=ngy
          !
          if (nm1==1) then
             i_in=1
             n_in=nl1
             name_in='x'
             ng_in=ngx
          endif
          do nd=1,sp%dim
             if (sp%d(nd)%name=='t') then
                nm3=nl3
                nav=nav+1
             endif
          enddo
          do nd=1,ndir_av
             if (dir_av(nd)=='t') then
                i_av(nd)=3
                n_av=ngt
                nm3=nl3
                nav=nav+1
             endif
          enddo
          if (nm3==1) then
             i_in=3
             n_in=nl3
             name_in='t'
             ng_in=ngt
          endif

       ! line Y (normal to plane XZ)
       ! ===========================
       case (2)
          ! PSD directions
          ! --------------
          do nd=1,sp%dim
             select case(sp%d(nd)%name)
             case ('y')
                sp%d(nd)%i=1
             case ('t')
                sp%d(nd)%i=3
                nd_t=nd
             case default
                call mpistop('bad choice for spectra direction: y-line -> choose y or t', 0)
             end select
          enddo

          ! MPI cutting of spectral direction (for 1D spectra)
          ! ---------------------------------
          select case(sp%d(1)%name)
          case ('y')
             ndom_sp=ndomy
             i_sp=2
          end select

          ! coordxi in MPI_CART
          ! -------------------
          coordx1=coordy
          coordx2=coordz

          ! ndom
          ! ----
          ndom1=ndomy
          ndom2=1

          ! dim per direction
          ! -----------------
          nl1=ny
          nl2=1
          if (sp%nbloc>1) then
             nl3=sp%d(nd_t)%lbloc
          else
             nl3=ngt
          endif
          ng1=ngy
          ng2=1
          ng3=nl3

          ! find dim for averaging and inhomogeneous direction
          ! --------------------------------------------------
          nm1=1
          nm2=1
          nm3=1
          i_in=0
          nav=0
          do nd=1,sp%dim
             if (sp%d(nd)%name=='y') then
                nm1=nl1
                nav=nav+1
             endif
          enddo
          do nd=1,ndir_av
             if (dir_av(nd)=='y') then
                i_av(nd)=1
                n_av=ngy
                nm1=nl1
                nav=nav+1
             endif
          enddo
          if (nm1==1) then
             i_in=1
             n_in=nl1
             name_in='y'
             ng_in=ngy
          endif
          ! Inhomogeneous direction superficialy taken a z
          i_in=2
          n_in=1
          name_in='z'
          ng_in=1
          !
          do nd=1,sp%dim
             if (sp%d(nd)%name=='t') then
                nm3=nl3
                nav=nav+1
             endif
          enddo
          do nd=1,ndir_av
             if (dir_av(nd)=='t') then
                i_av(nd)=3
                n_av=ngt
                nm3=nl3
                nav=nav+1
             endif
          enddo
          if (nm3==1) then
             i_in=3
             n_in=nl3
             name_in='t'
             ng_in=ngt
          endif


       ! line Z (normal to plane XY)
       ! ===========================
       case (3)
          ! PSD directions
          ! --------------
          do nd=1,sp%dim
             select case(sp%d(nd)%name)
             case ('z')
                sp%d(nd)%i=2
             case ('t')
                sp%d(nd)%i=3
                nd_t=nd
             case default
                call mpistop('bad choice for spectra direction: z-line -> choose z or t', 0)
             end select
          enddo

          ! MPI cutting of spectral direction (for 1D spectra)
          ! ---------------------------------
          select case(sp%d(1)%name)
          case ('z')
             ndom_sp=ndomz
             i_sp=3
          end select

          ! coordxi in MPI_CART
          ! -------------------
          coordx1=coordy
          coordx2=coordz

          ! ndom
          ! ----
          ndom1=1
          ndom2=ndomz

          ! dim per direction
          ! -----------------
          nl1=1
          nl2=nz
          if (sp%nbloc>1) then
             nl3=sp%d(nd_t)%lbloc
          else
             nl3=ngt
          endif
          ng1=1
          ng2=ngz
          ng3=nl3

          ! find dim for averaging and inhomogeneous direction
          ! --------------------------------------------------
          nm1=1
          nm2=1
          nm3=1
          i_in=0
          nav=0
          ! Inhomogeneous direction superficialy taken a y
          i_in=1
          n_in=1
          name_in='y'
          ng_in=1
          !
          do nd=1,sp%dim
             if (sp%d(nd)%name=='z') then
                nm2=nl2
                nav=nav+1
             endif
          enddo
          do nd=1,ndir_av
             if (dir_av(nd)=='z') then
                i_av(nd)=2
                n_av=ngz
                nm2=nl2
                nav=nav+1
             endif
          enddo
          if (nm2==1) then
             i_in=2
             n_in=nl2
             name_in='z'
             ng_in=ngz
          endif
          do nd=1,sp%dim
             if (sp%d(nd)%name=='t') then
                nm3=nl3
                nav=nav+1
             endif
          enddo
          do nd=1,ndir_av
             if (dir_av(nd)=='t') then
                i_av(nd)=3
                n_av=ngt
                nm3=nl3
                nav=nav+1
             endif
          enddo
          if (nm3==1) then
             i_in=3
             n_in=nl3
             name_in='t'
             ng_in=ngt
          endif
       end select

    endif


    ! Number of averaging direction can not be greater than 3
    ! -------------------------------------------------------
    if (nav>3) call mpistop('problem in pp_spectra: too many averaging !', 0)

    ! Check number of procs #2: division for transposition
    ! -------------------------
    do nd=1,sp%dim
       if ((sp%d(nd)%name=='x').and.(ndomx>1)) then
          ! Cutting of time dimension by ndomx (for transposition x <-> t)
          nt_txt=nl3/ndomx
          if (ndomx*nt_txt.ne.nl3) &
               call mpistop('transposition xt not possible: nmax/freq_plane is not divisible by ndomx',0)
       endif
       if ((sp%d(nd)%name=='z').and.(ndomz>1)) then
          ! Cutting of time dimension by ndomz (for transposition z <-> t)
          nt_tzt=nl3/ndomz
          if (ndomz*nt_tzt.ne.nl3) &
               call mpistop('transposition zt not possible: nmax/freq_plane is not divisible by ndomz',0)
       endif
    enddo

    ! Total number of samples (nb of points * nb of timesteps)
    ! ========================================================
    ng_tot=nl1*nl2*nl3*nproc

    ! Number of samples for averaging
    ! ===============================
    if (n_in>0) then
       ng_av=ng_tot/ng_in
    else
       ng_av=ng_tot ! default for 3D spectra
    endif

    ! Spectral resolution which includes one/two-sided factor
    ! =======================================================
    sp%dk=1.0_wp
    do nd=1,sp%dim
       if (sp%d(nd)%name=='x') sp%dk=sp%dk*dkx
       if (sp%d(nd)%name=='y') sp%dk=sp%dk*dky
       if (sp%d(nd)%name=='z') sp%dk=sp%dk*dkz
       if (sp%d(nd)%name=='t') sp%dk=sp%dk*dom
    enddo
    if (sp%is_onesided) sp%dk=sp%dk*0.5_wp    
    
    ! Determine communicator for PSD/averaging directions
    ! ===================================================
    select case(i_in)
       
    ! inhomogeneous direction is 1 (y for yz-plane/x for xz-plane or xy-plane)
    ! ----------------------------
    case (1)
       if (sp%dim==1) then
          if ((sp%d(1)%name=='y').or.(dir_av(1)=='y')) then            
             COMM_in=COMMXZ ! (MPI comm along y-direction)
          endif
          if ((sp%d(1)%name=='z').or.(dir_av(1)=='z')) then
             COMM_in=COMMXY ! (MPI comm along z-direction)
          endif
       elseif (sp%dim==2) then
          if ((sp%d(1)%name=='y').or.(sp%d(2)%name=='y')) then            
             COMM_in=COMMXZ ! (MPI comm along y-direction)
          endif
          if ((sp%d(1)%name=='z').or.(sp%d(2)%name=='z')) then
             COMM_in=COMMXY ! (MPI comm along z-direction)
          endif
       endif

    ! inhomogeneous direction is 2 (z for yz-plane or xz-plane/y for xy-plane)
    ! ----------------------------
    case (2) 
       if (sp%dim==1) then
          if ((sp%d(1)%name=='x').or.(dir_av(1)=='x')) then
             COMM_in=COMMYZ ! (MPI comm along x-direction)
          endif
          if ((sp%d(1)%name=='y').or.(dir_av(1)=='y')) then
             COMM_in=COMMXZ ! (MPI comm along y-direction)
          endif
       elseif (sp%dim==2) then
          if ((sp%d(1)%name=='x').or.(sp%d(2)%name=='x')) then
             COMM_in=COMMYZ ! (MPI comm along x-direction)
          endif
          if ((sp%d(1)%name=='y').or.(sp%d(2)%name=='y')) then
             COMM_in=COMMXZ ! (MPI comm along y-direction)
          endif
       endif

    ! inhomogeneous direction is 3 (t)
    ! ----------------------------
    case (3) 
       if (sp%dim==1) then
          if ((sp%d(1)%name=='x').or.(dir_av(1)=='x')) then
             COMM_in=COMMYZ ! (MPI comm along x-direction)
          endif
          if ((sp%d(1)%name=='y').or.(dir_av(1)=='y')) then
             COMM_in=COMMXZ ! (MPI comm along y-direction)
          endif
          if ((sp%d(1)%name=='z').or.(dir_av(1)=='z')) then
             COMM_in=COMMXY ! (MPI comm along z-direction)
          endif
       elseif (sp%dim==2) then
          if ((sp%d(1)%name=='x').or.(sp%d(2)%name=='x')) then
             COMM_in=COMMYZ ! (MPI comm along x-direction)
          endif
          if ((sp%d(1)%name=='y').or.(sp%d(2)%name=='y')) then
             COMM_in=COMMXZ ! (MPI comm along y-direction)
          endif
          if ((sp%d(1)%name=='z').or.(sp%d(2)%name=='z')) then
             COMM_in=COMMXY ! (MPI comm along z-direction)
          endif
       endif
       
    end select

    ! Print at screen spectra parameters
    ! ==================================
    if (iproc==0) then
       write(6,*) 'Spectra parameters:'
       write(6,*) '-------------------'
       if (snapshots(nsr)%type.eq.1) write(6,fmt='(1x,''~> line number:'',i4)') nplr
       if (snapshots(nsr)%type.eq.2) write(6,fmt='(1x,''~> plane number:'',i4)') nplr
       write(6,*) '   in directory: ',dirDATA
       if (snapshots(nsr)%normal==1) write(6,fmt='(1x,''   yz-plane with '',i7,'' samples'')') ngt
       if (snapshots(nsr)%normal==2) write(6,fmt='(1x,''   xz-plane with '',i7,'' samples'')') ngt
       if (snapshots(nsr)%normal==3) write(6,fmt='(1x,''   xy-plane with '',i7,'' samples'')') ngt
       if (sp%dim==1) print *,'~> autospectra in ',sp%d(1)%name,'-dimension for variable(s) ',sp%varname(1:sp%nvar)
       if (sp%dim >1) print *,'~> multidimensional spectra in ',sp%d%name,'-dimensions for variable(s) ',sp%varname(1:sp%nvar)
       if (sp%is_onesided) then
          print *,'   (one-sided spectrum)'
       else
          print *,'   (two-sided spectrum)'
       endif
       do nd=1,sp%dim
          write(6,fmt='(2x,''* direction #'',i1,'': '',a)') nd,sp%d(nd)%name
          ! block length
          write(6,fmt='(3x,''-> block length:'',i6)') sp%d(nd)%lbloc
          ! windowing informations
          if (sp%d(nd)%type_win=='T') then
             if (sp%d(nd)%param_win==0.0_wp) then
                print *,'  -> rectangular windowing'
             elseif (sp%d(nd)%param_win==1.0_wp) then
                print *,'  -> Hann windowing'
             else
                write(6,fmt='(3x,''-> Tukey windowing with coeff: '',f3.1)') sp%d(nd)%param_win
             endif
          elseif (sp%d(nd)%type_win=='K') then
             write(6,fmt='(3x,''-> Kaiser-Bessel windowing with coeff: '',f3.1)') sp%d(nd)%param_win
          endif
          if ((nd==sp%dim).and.(sp%is_capon)) write(6,fmt='(3x,''-> Capon MVSE with filter length: '',i4)') sp%ncapon
          ! info on resolution
          select case (sp%d(nd)%name)
          case('x')
             write(6,fmt='(3x,''-> wavenumber resolution: dkx='',f10.3,''/m'')') dkx             
          case('y')          
             write(6,fmt='(3x,''-> wavenumber resolution: dky='',f10.3,''/m'')') dky
          case('z')                    
             write(6,fmt='(3x,''-> wavenumber resolution: dkz='',f10.3,''/m'')') dkz
          case('t')
             if (sp%nbloc>1) then
                print *,'  -> apply Welch''s method'
                if (sp%is_overlap) then
                   write(6,fmt='(3x,''-> overlapping:'',l2,''-> of length:'',i6)') sp%is_overlap,sp%loverlap
                   write(6,fmt='(3x,''-> number of blocks after overlap:'',i3)') sp%nbloc
                else
                   write(6,fmt='(3x,''-> overlapping:'',l2)') sp%is_overlap
                   write(6,fmt='(3x,''-> number of blocks:'',i3)') sp%nbloc
                endif
             endif
             write(6,fmt='(3x,''-> frequency resolution: df='',f10.3,'' Hz'')') dom/twopi
          end select         
       enddo
    endif

  end subroutine spectra_param

end module mod_pp_sp_main
