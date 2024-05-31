!==============================================================================
module mod_pp_sp_read_old
!==============================================================================
  !> Module to read param file and planes for post-processing (spectra)
  !> * old format for compatibility with old databases *
!==============================================================================
  use warnstop
  use mod_flow
  use mod_time
  use mod_constant
  use mod_pp_mpi
  use mod_pp_var
  implicit none
  ! ---------------------------------------------------------------------------
  integer, private :: nit
  integer, private :: ngx_,ngz_ ! <- for old reli_plane
  ! boolean indicators
  logical, private :: isxy=.false.,isxz=.false.
  logical, private :: isvites=.false.,ispress=.false.
  logical, private :: ismodal_nrj=.false.,ismoyz=.false.
  logical, private :: ispwall=.false.    ! k-w spectra wall pressure
  logical, private :: ispwall_tr=.false. ! k-w spectra wall pressure / per class
  logical, private :: ispwall_autosp     ! autospectra wall pressure
  logical, private :: isreli=.false.     ! relecture/ecriture sans calcul
  logical, private :: decoup_orig        ! use old MPI partitioning
  character(len=30), private :: cas      ! character strings
  ! ---------------------------------------------------------------------------

contains

  !============================================================================
  subroutine read_param_old
  !============================================================================
    !> Read parameters of old simulations
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: icase
    integer :: hx1,hx2,hy1 ! indices of triggering step
    integer :: nprocx_old,nprocy_old,nprocz_old ! MPI partitioning of the old sim.
    integer :: nclass_tot,ndeb_class ! class indices
    !integer :: nplan_xz
    real(wp) :: hm ! height of the triggering step
    real(wp) :: Mach,Uo ! Mach and flow speed
    logical :: icurv=.false.  
    ! -------------------------------------------------------------------------
    
    ! Configurations
    ! ==============
    icase=1
    !******
    if (icase==1) then     ! M=0.5 ZPG
       cas='M05_ZPG'
    elseif (icase==2) then ! M=0.7 ZPG
       cas='M07_ZPG'
    elseif (icase==3) then ! M=0.9 ZPG
       cas='M09_ZPG'
    elseif (icase==4) then ! M=0.3 ZPG
       cas='M03_ZPG'
    elseif (icase==5) then ! M=0.5 APGw
       cas='M05_APGw'
    elseif (icase==6) then ! M=0.5 APGs
       cas='M05_APGs'
    elseif (icase==7) then ! M=0.5 FPGw
       cas='M05_FPGw'
    elseif (icase==8) then ! M=0.5 FPGs
       cas='M05_FPGs'
    elseif (icase==9) then ! M=0.5 FPGs
       cas='M07L'
    endif
    ! indicateur maillage cartesien ou curviligne
    if ((icase==1).or.(icase==2).or.(icase==3).or.(icase==4).or.(icase==9)) then
       icurv=.false.
    elseif ((icase==5).or.(icase==6).or.(icase==7).or.(icase==8)) then
       icurv=.true.
    endif

    ! Read parameters
    ! ===============

    open(30,file='param_'//trim(cas)//'.ini')
    read(30,*) !'parameters of the old simulation'
    read(30,*) !'--------------------------------'
    read(30,*) !'Mach number'
    read(30,*) Mach
    read(30,*) !'free-stream velocity: Uo'
    read(30,*) Uo
    read(30,*) !'old global mesh: ngx,ngy,ngz'
    read(30,*) ngx,ngy,ngz
    read(30,*) !'old MPI partitioning: nprocx,nprocy,nprocz'
    read(30,*) nprocx_old,nprocy_old,nprocz_old
    read(30,*) !'grid sizes: deltax,deltay,deltaz'
    read(30,*) deltax,deltay,deltaz
    read(30,*) !'small triggering step: hx1,hx2,hy1'
    read(30,*) hx1,hx2,hy1
    read(30,*) !'height of the step: hm'
    read(30,*) hm
    read(30,*) !'characteristics of simulation classes:'
    read(30,*) !'-> time step: deltat'
    read(30,*) deltat
    read(30,*) ! -> number of iterations per class: nit
    read(30,*) nit
    read(30,*) ! -> output frequency for planes: freq_plane
    read(30,*) freq_plane
    read(30,*) ! -> number of recorded classes: nclass
    read(30,*) nclass_tot
    read(30,*) ! -> index of the first class: ndeb_class
    read(30,*) ndeb_class
    close(30)

    ! Print parameters of old simulations
    ! ===================================
    if (iproc==0) then
       write(6,*) 'Informations about the old simulation'
       write(6,*) '--------------------------------'
       write(6,fmt='(1x,''Mach number='',F5.2,'', Uo='',F10.4)') Mach,Uo
       write(6,fmt='(1x,''old global mesh: '',i5,i5,i5)') ngx,ngy,ngz
       write(6,fmt='(1x,''old MPI partitioning:'',i4,i4,i4)') nprocx_old,nprocy_old,nprocz_old
       write(6,fmt='(1x,''small triggering step: hx1='',i3,'', hx2='',i3,'', hy1='',i3)') hx1,hx2,hy1
       write(6,fmt='(1x,''            of height:'',E14.7)') hm    
       write(6,fmt='(1x,''-> grid sizes:, Dx='',E14.7,'', Dy='',E14.7,'', Dz='',E14.7)') deltax,deltay,deltaz
       write(6,fmt='(1x,''-> time step:'',E14.7)') deltat
       write(6,fmt='(1x,''-> number of iterations per class:'',i6)') nit
       write(6,fmt='(1x,''-> output frequency for planes:'',i4)') freq_plane
       write(6,fmt='(1x,''-> number of recorded classes:'',i4)') nclass_tot
       write(6,fmt='(1x,''-> index of the first class:'',i4)') ndeb_class
    endif

    ! Streamwise dimension after plane rewriting
    ! ==========================================
    if (icase==9) then
       ngx=90*55
       !ngx=90*8
       !ngz=400
    else
       ngx=ngx/2 
    endif

  end subroutine read_param_old

  !============================================================================
  subroutine read_param_pp_old
  !============================================================================
    !> Determine post-processing parameters
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: lbloc
    integer :: ndebx,ndeby,ndebz
    ! -------------------------------------------------------------------------

    ! Directories
    ! ===========
    dirDATA='./'
    dirRESU='./'
    dirDATA = '/run/media/xavier/anterak27/plan_xz_pression_parietale/'
    !dirDATA = '../../CLT_MACH/DATA_M07/'

    ! Read post-processing parameters
    ! ===============================
    open(30,file='param_pp_'//trim(cas)//'.ini')
    read(30,*) ! Post-processing parameters
    read(30,*) ! ===============================================
    read(30,*) ! number of selected classes: nclass
    read(30,*) nclass
    read(30,*) ! size of temporal blocks: lbloc
    read(30,*) lbloc
    read(30,*) ! overlapping (Welch's method): is_overlap
    read(30,*) sp%is_overlap, sp%loverlap
    read(30,*) ! Capon's method: is_capon
    read(30,*) sp%is_capon
    read(30,*) ! Capon filter length (Capon's order)
    read(30,*) sp%ncapon
    read(30,*) ! number of MPI domains per direction: ndomx,ndomy,ndomz
    read(30,*) ndomx
    read(30,*) ndomy
    read(30,*) ndomz
    read(30,*) ! index of the first block in each directions: ndebx,ndeby,ndebz (ex: ndebz=1; 1st bottom block)
    read(30,*) ndebx
    read(30,*) ndeby
    read(30,*) ndebz
    read(30,*) ! Choice of actions
    read(30,*) ! -----------------
    read(30,*) ! logical decoup_orig
    read(30,*) decoup_orig
    read(30,*) ! logical isvites
    read(30,*) isvites
    read(30,*) ! logical ispress
    read(30,*) ispress
    read(30,*) ! logical isxy
    read(30,*) isxy
    read(30,*) ! logical isxz
    read(30,*) isxz
    read(30,*) ! logical ispwall
    read(30,*) ispwall
    read(30,*) ! logical ispwall_tr / per class
    read(30,*) ispwall_tr
    read(30,*) ! logical ispwall_autosp / compute autospectra
    read(30,*) ispwall_autosp
    read(30,*) ! logical isreli
    read(30,*) isreli
    read(30,*) ! modal energy (for transition): ismodal_nrj 
    read(30,*) ismodal_nrj
    close(30)

    if (isreli) decoup_orig=.true.

    if (ispwall_autosp) sp%dim=1
    if (ispwall) sp%dim=3
    if (ispwall_tr) call mpistop('ispwall_tr NOT YET INCLUDED',0)
    
    allocate(sp%d(sp%dim))
    
    sp%d(1)%lbloc=lbloc
    sp%d(1)%type_win='T'
    sp%d(1)%param_win=1.0_wp
    sp%is_onesided=.true.

  end subroutine read_param_pp_old

  !============================================================================
  subroutine read_plane_part_old(nbl)
  !============================================================================
    !> Read plane for 1 time block and MPI partitioning
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer, intent(in) :: nbl
    ! -------------------------------------------------------------------------
    integer :: i,k,n,nn,n1,n2
    integer :: ind
    real(wp), dimension(:,:), allocatable :: varxz
    ! -------------------------------------------------------------------------

    ! print screen
    if (iproc==0) print 103,nbl

    if (nbl==1) then
       !open(51,file=trim(dirDATA)//'pp_'//trim(cas)//'_lf.bin',form='unformatted',status='unknown')
       open(51,file=trim(dirDATA)//'pp_'//trim(cas)//'.bin',form='unformatted',status='unknown')
       rewind(51)
       read(51) ngx_
       read(51) ngz_
       read(51) ind
       read(51)

       allocate(xg(ngx_),yg(ngy),zg(ngz_))
       read(51) (xg(i),i=1,ngx_)
       read(51) (zg(k),k=1,ngz_)

       if (iproc==0) print*,ngx_,ngz_,ind

!!$       allocate(var_r(nx,nz,nl3))
!!$       if ((ispwall).or.(ispwall_tr)) allocate(var_c(nx,nz,nl3))
!!$       ! array for transposition #1: z <-> t
!!$       if ((ispwall).or.(ispwall_tr)) allocate(var_ct1(nx,nt,ngz))
!!$       ! array for transposition #2: x <-> z
!!$       if (ispwall) allocate(var_ct2(nz_t,nt,ngx),Swk1k3(nz_t,nt,ngx))
!!$       ! autospectra & moments
!!$       if (ispwall_autosp) allocate(varm_r(nx,nz,nl3))
!!$       if (ispwall_autosp) allocate(var_rms(nx),var_skew(nx),var_kurt(nx))
!!$       if (ispwall_autosp) allocate(varm_rms(nx),varm_skew(nx),varm_kurt(nx))
    endif

    allocate(varxz(ngx_,ngz_))

    if (sp%is_overlap) then
       if (nbl==1) then
          n1=1
          n2=nl3
       else
          var_r(:,:,1:nl3-sp%loverlap)=var_r(:,:,sp%loverlap+1:nl3)
          nn=(nbl-1)*sp%loverlap
          n1=nn+1+nl3-sp%loverlap;
          n2=nn+nl3;
       endif
    else
       n1=(nbl-1)*nl3+1
       n2=nbl*nl3
    endif
    if (iproc==0) print*,n1,n2

    do n=n1,n2
       if ((sp%is_overlap).and.(nbl>1)) then
          nn=n-(nbl-1)*sp%loverlap
       else
          nn=n-(nbl-1)*nl3
       endif
       if ((iproc==0).and.(mod(n,100)==0)) print *,'   it:',n,nn
       read(51) ((varxz(i,k),i=1,ngx_),k=1,ngz_)
       call MPI_BARRIER(COMM_global,info)

       do k=1,nz
          do i=1,nx
             var_r(i,k,nn)=varxz(i+coord(1)*nx,k+coord(3)*nz)
          enddo
       enddo
       call MPI_BARRIER(COMM_global,info)
    enddo

    if ((ispwall).or.(ispwall_tr)) var_c=var_r
    
    deallocate(varxz)

    ! print screen
    if (iproc==0) print 104,nbl
103 format(1x,'reading block: ',i3)
104 format(1x,'block ',i3,' read')

  end subroutine read_plane_part_old
  
end module mod_pp_sp_read_old
