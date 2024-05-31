!===============================================================================
module mod_block
!===============================================================================
  !> Module for Block partitioning
!===============================================================================
  use precision
  implicit none
  ! ----------------------------------------------------------------------------
  ! Type for "Block snapshots" attributes
  ! ----------------------------------------------------------------------------
  type snapshot_bloc
     integer :: ind_i1, ind_i2          ! volume position: i indexes
     integer :: ind_j1, ind_j2          ! volume position: j indexes
     integer :: ind_k1, ind_k2          ! volume position: k indexes
     integer :: nvar                    ! number of variables
     integer :: freq                    ! writing frequency
     character(6), dimension(20) :: var ! name of variables
  end type snapshot_bloc
  ! ----------------------------------------------------------------------------
  !> Type for Block attributes
  ! ----------------------------------------------------------------------------
  type bloc
     ! grid sizes in i,j,k-directions
     integer :: ni,nj,nk
     ! number of procs in i,j,k-directions
     integer :: ndomi,ndomj,ndomk
     ! total number of procs for the block
     integer :: nproc
     ! global rank of first and last procs in a block
     integer :: proc_min,proc_max ! used to define interblock neighbors
     ! types of Boundary Conditions for each faces of the block
     ! 1 : imin / i=1  / left   / W(est)
     ! 2 : imax / i=nx / right  / E(ast)
     ! 3 : jmin / j=1  / bottom / S(outh)
     ! 4 : jmax / j=ny / top    / N(orth)
     ! 5 : kmin / k=1  / front  / F(ront)
     ! 6 : kmax / k=nz / back   / B(ack)
     integer, dimension(6) :: BC
     character, dimension(6) :: flag
     logical, dimension(6) :: periodic
     ! corner coordinates
     real(wp), dimension(8) :: xcorner,ycorner,zcorner
     ! definition of sponge zone
     logical :: is_sponge
     integer :: is1,is2,js1,js2,ks1,ks2,d_is,d_js,d_ks
     ! definition of outputs for each block
     integer :: nsnapshot
     type(snapshot_bloc), dimension(:), allocatable :: snapshot
     ! stretching param. for def. of grid in the solver (is_def_grid)
     integer, dimension(:), allocatable :: nrx1,nrx2,nry1,nry2
     real(wp), dimension(:), allocatable :: rsx,rsy
  end type bloc
  ! ----------------------------------------------------------------------------
  !> Type for Block PP attributes
  ! ----------------------------------------------------------------------------
  type bloc_snap_pp
     type(snapshot_bloc), dimension(:), allocatable :: snap
     integer :: nsn_pp
  end type bloc_snap_pp
  ! ----------------------------------------------------------------------------
  ! number of blocks
  integer :: nbloc
  ! block attributes
  type(bloc), dimension(:), allocatable, target :: bl,blr
  type(bloc_snap_pp), dimension(:), allocatable, target :: blc_snap_pp
  ! ----------------------------------------------------------------------------
  !> Block restriction for postprocessing
  ! ----------------------------------------------------------------------------
  integer :: nbloc_r ! number of blocks which are read
  integer, dimension(:), allocatable :: nblr ! number of the read blocks
  ! ----------------------------------------------------------------------------

contains
  
  !============================================================================
  subroutine read_param_blocks(paramfile)
  !============================================================================
    !> Reading block partitioning in param_blocks.ini
    !> TO BE CHANGED: protect reading from wrong formating
    !  (cf readline in read_write_info)
    !============================================================================
    use mod_ngh
    use warnstop
    use mod_mpi
    use mod_constant
    implicit none
    ! --------------------------------------------------------------------------
    character(len=*), intent(in) :: paramfile
    ! --------------------------------------------------------------------------
    logical :: iexist,ierr
    integer :: n,i,isn
    character(len=4) :: eol
    ! --------------------------------------------------------------------------

    ! Check if param_blocks.ini exists
    ! ================================
    inquire(file=trim(paramfile),exist=iexist)
    if (.not.iexist) then
       call mpistop('Paramfile for blocks does not exist!', 0)
    endif

    ! Read param_blocks.ini
    ! =====================
    open(30,file=trim(paramfile))

    ! Documentation part
    ! ------------------

    read(30,*)!=============================================================
    read(30,*)!=============================================================
    read(30,*)! MUSICA2 : fill Block definitions
    read(30,*)!=============================================================
    read(30,*)!=============================================================
    read(30,*)! 1/ Specify number of block
    read(30,*)! 2/ Paste and copy the 'Block #x' section as many times
    read(30,*)!    as the number of blocks
    read(30,*)! 3/ Fill parameters for each block
    read(30,*)!    --> number of points and number of processors
    read(30,*)!    --> boundary conditions and connectivity
    read(30,*)!    --> sponge zone characteristics (if present)
    read(30,*)!    --> define output snapshots (if present)
    read(30,*)!        ~> list as many snapshot as necessary by paste and
    read(30,*)!           copy of the corresponding line
    read(30,*)!        ~> can be: point, line, plane or volume
    read(30,*)!        ~> can be over full block or subdomain
    read(30,*)!=============================================================
    read(30,*)! BOUNDARY CONDITIONS & CONNECTIVITY:
    read(30,*)! --- first line ---------------------------------------------
    read(30,*)! n (>0): number of the neighbouring block [connectivity]
    read(30,*)!   0   : wall boundary condition
    read(30,*)!  -1   : non-reflecting condition of Tam & Dong
    read(30,*)!  -2   : outflow condition of Tam & Dong
    read(30,*)!  -3   : non-reflecting characteristic condition
    read(30,*)!  -4   : turb. inflow BC (imposed p_tot,T_tot,velocity dir.)
    read(30,*)!  -5   : back-pressure outflow BC
    read(30,*)! --- second line --------------------------------------------
    read(30,*)!   p   : periodic boundary condition [only for n>0]
    read(30,*)!   s   : slip wall [only for n=0]
    read(30,*)!   r   : impose reference quantities in BC [only for n<0]
    read(30,*)!   -   : no flag
    read(30,*)!=============================================================
    read(30,*)! SPONGE ZONE DEFINITION:
    read(30,*)! 1/ logical is_sponge (T/F): T if block contains sponge zone
    read(30,*)! 2/ bound indices in grid (before MPI partitioning)
    read(30,*)!    [is1:is2]x[js1:js2]x[ks1:ks2]
    read(30,*)!        in I-direction: from is1 to is2
    read(30,*)!        in J-direction: from js1 to js2
    read(30,*)!        in K-direction: from ks1 to ks2
    read(30,*)! 3/ number of points on which a progressive damping is applied
    read(30,*)!        in I-direction: d_is
    read(30,*)!        in J-direction: d_js
    read(30,*)!        in K-direction: d_ks
    read(30,*)!=============================================================
    read(30,*)! SNAPSHOT DEFINITION:
    read(30,*)! positions (index in block): I1 to I2, J1 to J2, K1 to K2
    read(30,*)! ---------------------------
    read(30,*)!    If point: I1=I2 and J1=J2 and K1=K2
    read(30,*)!    If line: I1=I2/J1=J2/K1.ne.K2 or I1=I2/J1.ne.J2/K1=K2 or ...
    read(30,*)!    If plane: I1=I2/J1.ne.J2/K1.ne.K2 or ...
    read(30,*)!    If volume: I1.ne.I2 and J1.ne.J2 and K1.ne.K2
    read(30,*)! nfreq: frequency for writting snapshots
    read(30,*)! ------
    read(30,*)!    if nfreq>0, output each nfreq iterations
    read(30,*)!    if 0, put to freq_point if point or freq_plane if plane or ...
    read(30,*)!    if nfreq<0, output each ntotal/(-nfreq) iterations
    read(30,*)! nvar: number of variables (must be < 20)
    read(30,*)! -----
    read(30,*)! list of possible variables:
    read(30,*)! ---------------------------
    read(30,*)!    prs,uu,vv,ww,rho,Tmp,div,Mach,Gamma,Frhov,Grhow,udf
    read(30,*)!    [if not in the list, should be added in mod_io_snapshots.f90]
    read(30,*)!    [udf if for user-defined variables stored in uvar(:,:,:,n)]
    read(30,*)!
    read(30,*)!=============================================================
    read(30,*)!=============================================================

    ! Number of blocks
    ! ----------------
    read(30,*)! nbloc: Number of Blocks
    read(30,*) nbloc
    allocate(bl(nbloc))

    ! Block definitions
    ! -----------------
    do n=1,nbloc
       if (n==1) &
       read(30,*)!==========================================================
       read(30,*)! Block # n
       read(30,*)!==========================================================
       read(30,*)! dimensions + MPI Partitioning
       read(30,*) bl(n)%ni, bl(n)%ndomi
       read(30,*) bl(n)%nj, bl(n)%ndomj
       read(30,*) bl(n)%nk, bl(n)%ndomk
       
       ! Numbers of procs per block
       ! --------------------------
       bl(n)%nproc=bl(n)%ndomi*bl(n)%ndomj*bl(n)%ndomk
       
       read(30,*)!----------------------------------------------------------
       read(30,*)! Boundary conditions & connectivity
       read(30,*)! imin,imax,jmin,jmax,kmin,kmax
       read(30,*) bl(n)%BC
       read(30,*) bl(n)%flag
       
       read(30,*)!----------------------------------------------------------
       read(30,*)! Definition of Sponge zone
       read(30,*) bl(n)%is_sponge, &
                  bl(n)%is1,bl(n)%is2,bl(n)%js1,bl(n)%js2,bl(n)%ks1,bl(n)%ks2, &
                  bl(n)%d_is, bl(n)%d_js, bl(n)%d_ks
       
       read(30,*)!----------------------------------------------------------
       read(30,*)! Define output snapshots:
       read(30,*) bl(n)%nsnapshot
       read(30,*)!----|----|----|----|----|----|------|------|--------------
       read(30,*)! I1 | I2 | J1 | J2 | K1 | K2 | freq | nvar | list var
       read(30,*)!    |    |    |    |    |    |      |   n  | name ...
       read(30,*)!----|----|----|----|----|----|------|------|--------------
       if (bl(n)%nsnapshot.gt.0) then
          ! allocate plane type for block
          allocate(bl(n)%snapshot(bl(n)%nsnapshot))
          ! read infos
          do isn=1,bl(n)%nsnapshot
             read(30,*) bl(n)%snapshot(isn)%ind_i1, &
                        bl(n)%snapshot(isn)%ind_i2, &
                        bl(n)%snapshot(isn)%ind_j1, &
                        bl(n)%snapshot(isn)%ind_j2, &
                        bl(n)%snapshot(isn)%ind_k1, &
                        bl(n)%snapshot(isn)%ind_k2, &
                        bl(n)%snapshot(isn)%freq,   &
                        bl(n)%snapshot(isn)%nvar,   &
                        bl(n)%snapshot(isn)%var(1:bl(n)%snapshot(isn)%nvar)
          enddo
       endif

       ! Finish reading block
       ! --------------------
       eol='' ! read separation line for next block
       do while (eol.ne.'!===')
          read(30,*) eol
       enddo
 
    enddo ! end loop block
    close(30)

    ! Few consistency checks
    ! ======================

    ! check consistency of flag indicators
    ierr=.false.

    if (iproc==0) then
       do n=1,nbloc
          do i=1,6
             if ((bl(n)%flag(i)=='p').and.(bl(n)%BC(i)<=0)) then
                print 99,i,n
                print *,'~> periodicity can be enforced only for interior BC [n>0]'
                ierr=.true.
             endif
             if ((bl(n)%flag(i)=='s').and.(bl(n)%BC(i).ne.0)) then
                print 99,i,n
                print *,'~> slip wall can be enforced only for wall BC [n=0]'
                ierr=.true.
             endif
             if ((bl(n)%flag(i)=='r').and.(bl(n)%BC(i)>=0)) then
                print 99,i,n
                print *,'~> reference quantities can be imposed only for exterior BC [n<0]'
                ierr=.true.
             endif
          enddo
       enddo
    endif
99  format(1x,'Problem in flag definition for BC #',i1,' in block',i4)

    if (ierr) call mpistop('Please check param_blocks.ini',0)
    
    ! all blocks should have the same dimension along k (except in 3D curvilinear)
    if (.not.is_curv3) then
       do n=2,nbloc
          if (bl(n)%nk.ne.bl(1)%nk) &
               call mpistop('all blocks should have the same dimension in z !',0)
       enddo
    endif

    ! Determine proc_min and proc_max
    ! -------------------------------
    bl(1)%proc_min=0
    bl(1)%proc_max=bl(1)%nproc-1
    do n=2,nbloc
       bl(n)%proc_min=bl(n-1)%proc_max+1
       bl(n)%proc_max=bl(n)%proc_min+bl(n)%nproc-1
    enddo

    ! If grid is generated in the solver
    if (is_def_grid) call read_param_grid

  end subroutine read_param_blocks

  !============================================================================
  subroutine restrict_blocks
  !============================================================================
    !> Reading block partitioning in param_blocks.ini
    !> TO BE CHANGED: protect reading from wrong formating
    !  (cf readline in read_write_info)
    !============================================================================
    use mod_ngh
    use warnstop
    implicit none
    ! --------------------------------------------------------------------------
    integer :: n,nb,ibc,nt,isn
    ! --------------------------------------------------------------------------

    ! Create temporary structure for read blocks
    ! ==========================================
    ! nbloc_r & nblr are read n param_pp.ini
    allocate(blr(nbloc_r))

    ! Set attributes to read blocks
    ! =============================
    nb=1
    do while (nb<=nbloc_r)
       do n=1,nbloc
          if (n==nblr(nb)) then

             ! dimensions
             ! ----------
             blr(nb)%ni=bl(n)%ni
             blr(nb)%nj=bl(n)%nj
             blr(nb)%nk=bl(n)%nk

             blr(nb)%ndomi=bl(n)%ndomi
             blr(nb)%ndomj=bl(n)%ndomj
             blr(nb)%ndomk=bl(n)%ndomk

             blr(nb)%nproc=bl(n)%nproc
             if (nb.eq.1) then
                blr(nb)%proc_min = 0
                blr(nb)%proc_max = blr(nb)%nproc-1
             else
                blr(nb)%proc_min = blr(nb-1)%proc_max + 1
                blr(nb)%proc_max = blr(nb)%proc_min + blr(nb)%nproc - 1
             endif

             ! boundary conditions
             ! -------------------
             blr(nb)%BC=bl(n)%BC

             ! Modification of BCs & sponge zone
             ! ---------------------------------
             ! Kept identical if wall, or periodic, or neighbor of a post-processed block
             ! Put to -1 if not
             do ibc=1,6
                if (blr(nb)%BC(ibc).eq.0) then ! wall
                   blr(nb)%BC(ibc)=0
                else if (blr(nb)%BC(ibc).eq.n) then ! periodic
                   blr(nb)%BC(ibc) = nb
                else
                   isn=0
                   loopnt: do nt=1,nbloc_r
                      if (blr(nb)%BC(ibc).eq.nblr(nt)) then
                         blr(nb)%BC(ibc) = nt ! neighbor of another post-processed block
                         isn=1
                         exit loopnt
                      endif
                   enddo loopnt
                   if (isn.eq.0) blr(nb)%BC(ibc) = -1 ! if not, put to -1
                endif
             enddo
             blr(nb)%is_sponge = .false.

             ! snapshot characteristics
             ! ------------------------
             blr(nb)%nsnapshot=bl(n)%nsnapshot
             allocate(blr(nb)%snapshot(bl(n)%nsnapshot))
             ! read infos
             do isn=1,bl(n)%nsnapshot
                blr(nb)%snapshot(isn)%ind_i1=bl(n)%snapshot(isn)%ind_i1
                blr(nb)%snapshot(isn)%ind_i2=bl(n)%snapshot(isn)%ind_i2
                blr(nb)%snapshot(isn)%ind_j1=bl(n)%snapshot(isn)%ind_j1
                blr(nb)%snapshot(isn)%ind_j2=bl(n)%snapshot(isn)%ind_j2
                blr(nb)%snapshot(isn)%ind_k1=bl(n)%snapshot(isn)%ind_k1
                blr(nb)%snapshot(isn)%ind_k2=bl(n)%snapshot(isn)%ind_k2
                blr(nb)%snapshot(isn)%freq  =bl(n)%snapshot(isn)%freq
                blr(nb)%snapshot(isn)%nvar  =bl(n)%snapshot(isn)%nvar
                blr(nb)%snapshot(isn)%var   =bl(n)%snapshot(isn)%var
             enddo

             ! increment read block number
             ! ---------------------------
             nb=nb+1

             if (nb>nbloc_r) exit
          endif
       enddo
    enddo

    ! Delete bl structure and replace by blr
    ! ======================================
    deallocate(bl)
    nbloc=nbloc_r
    allocate(bl(nbloc))

    ! copy structure
    do n=1,nbloc

       ! dimensions
       ! ----------
       bl(n)%ni=blr(n)%ni
       bl(n)%nj=blr(n)%nj
       bl(n)%nk=blr(n)%nk

       bl(n)%ndomi=blr(n)%ndomi
       bl(n)%ndomj=blr(n)%ndomj
       bl(n)%ndomk=blr(n)%ndomk

       bl(n)%nproc=blr(n)%nproc
       bl(n)%proc_min=blr(n)%proc_min
       bl(n)%proc_max=blr(n)%proc_max

       ! boundary conditions
       ! -------------------
       bl(n)%BC=blr(n)%BC

       ! snapshot characteristics
       ! ------------------------
       bl(n)%nsnapshot=blr(n)%nsnapshot
       allocate(bl(n)%snapshot(bl(n)%nsnapshot))
       ! read infos
       do isn=1,bl(n)%nsnapshot
          bl(n)%snapshot(isn)%ind_i1=blr(n)%snapshot(isn)%ind_i1
          bl(n)%snapshot(isn)%ind_i2=blr(n)%snapshot(isn)%ind_i2
          bl(n)%snapshot(isn)%ind_j1=blr(n)%snapshot(isn)%ind_j1
          bl(n)%snapshot(isn)%ind_j2=blr(n)%snapshot(isn)%ind_j2
          bl(n)%snapshot(isn)%ind_k1=blr(n)%snapshot(isn)%ind_k1
          bl(n)%snapshot(isn)%ind_k2=blr(n)%snapshot(isn)%ind_k2
          bl(n)%snapshot(isn)%freq  =blr(n)%snapshot(isn)%freq
          bl(n)%snapshot(isn)%nvar  =blr(n)%snapshot(isn)%nvar
          bl(n)%snapshot(isn)%var   =blr(n)%snapshot(isn)%var
       enddo
    enddo

    deallocate(blr)

  end subroutine restrict_blocks


  !===============================================================================
  subroutine read_param_grid
  !===============================================================================
    !> author: AB
    !> date: March 2024
    !> Reading grid parameters if grid is generated in the solver (is_def_grid)
  !===============================================================================
    use mod_grid
    use mod_utils
    use warnstop
    implicit none
    ! ---------------------------------------------------------------------------
    logical :: iexist
    integer :: n,i
    ! ---------------------------------------------------------------------------
    real(wp), dimension(:,:), allocatable :: temp_array

    ! Read param_grid.ini file
    ! ========================
    inquire(file="param_grid.ini", exist=iexist)
    if (.not.iexist) then
       call mpistop('param_grid file for simple grid does not exist!', 0)
    endif

    ! Allocation
    allocate(deltax0(nbloc)); allocate(deltay0(nbloc))
    allocate(x0(nbloc)); allocate(y0(nbloc))
    allocate(nstretchx(nbloc)); allocate(nstretchy(nbloc))

    ! ---------------
    ! TO BE COMPLETED
    ! ---------------
    open(30,file="param_grid.ini")
    rewind(30)
    read(30,*) !===============================================================================================
    read(30,*) !===============================================================================================
    read(30,*) ! MUSICA2 : fill advanced grid parameters for simple geometry
    read(30,*) !===============================================================================================
    read(30,*) !===============================================================================================
    read(30,*) ! ***** WORK IN PROGRESS *****  TO BE DEVELOPED
    read(30,*) ! ----------------------------------------------------------------------------------------------
    read(30,*) ! Grid settings: generate basic meshes in the solver
    read(30,*) ! ----------------------------------------------------------------------------------------------
    read(30,*) !   -> Cartesian grid, with different strechings in the different directions
    read(30,*) !   -> Circular grid
    read(30,*) !   -> Sinusoidal grid
    read(30,*) !   -> And others to be implemented...
    read(30,*) !===============================================================================================
    read(30,*) !===============================================================================================
    read(30,*) !                                         Grid settings
    read(30,*) !===============================================================================================
    read(30,*) ! Generic type of grid: nb_grid ***** WORK IN PROGRESS *****
    read(30,*) ! 1: cartesian
    read(30,*) ! 2: polar grid for cylinder [to be completed...]
    read(30,*) num_grid
    read(30,*) ! ----------------------------------------------------------------------------------------------
    read(30,*) ! General grid informations
    read(30,*) ! ----------------------------------------------------------------------------------------------
    read(30,*) ! Grid sizes at origin & origin coordinates: deltax, deltay, x0, y0
    read(30,*) !  --> Paste and copy this line as many times as the number of blocks
    do n=1,nbloc
       read(30,*) deltax0(n), deltay0(n), x0(n), y0(n) ! # block n
    enddo
    read(30,*) ! Spanwise extrusion origin: index k
    read(30,*) ! k_extr
    read(30,*) ! ----------------------------------------------------------------------------------------------
    read(30,*) ! Stretching parameters:  ***** WORK IN PROGRESS *****
    read(30,*) !  --> Paste and copy the 'Block #x' line as many times as the number of blocks
    read(30,*) !  --> 1 line for direction i & 1 line for direction j
    read(30,*) !  --> 1 line is:
    read(30,*) !      ~> number of different strechings for this direction (nstretch)
    read(30,*) !      ~> for 1 streching, index where to begin (nr1) and where to stop (nr2), along with
    read(30,*) !         streching coefficient (r)
    read(30,*) !      ~> [nr1, nr2, r] replicate as the number of strechings (nstretch)
    read(30,*) ! ----------------------------------------------------------------------------------------------
    read(30,*) ! Block #1:
    read(30,*) ! Stretching in direction i: nstretchx,  nri1,nri2,ri (X nstretchx)
    do n=1,nbloc
       read(30,*) nstretchx(n)
       ! If stretching specified in the x direction
       if (nstretchx(n).gt..0) then
          rewind(30)
          do i=1,39+n+nbloc !number of lines before this one
             read(30,*) !
          enddo
          ! Allocation of temporary array
          allocate(temp_array(3,nstretchx(n)))
          ! Read line
          read(30,*) nstretchx(n), temp_array
          ! Rearrange array with new ones
          allocate(bl(n)%nrx1(nstretchx(n)),bl(n)%nrx2(nstretchx(n)),bl(n)%rsx(nstretchx(n)))
          do i=1,nstretchx(n)
             bl(n)%nrx1(i)=temp_array(1,i); bl(n)%nrx2(i)=temp_array(2,i); bl(n)%rsx(i)=temp_array(3,i)
          enddo
          deallocate(temp_array)
          ! Check
          do i=1,nstretchx(n)
             if (bl(n)%nrx2(i).le.bl(n)%nrx1(i)) call mpistop("Stretching number "//trim(numchar(i))//" of block "//trim(numchar(n))//" in x: nrx2 must be greater than nrx1",0)
             if (bl(n)%rsx(i).le.0.0_wp) call mpistop("Stretching number "//trim(numchar(i))//" of block "//trim(numchar(n))//" in x: stretching ratio rsx must be > 0",0)
             if (bl(n)%rsx(i).gt.2.0_wp) call mpistop("Stretching number "//trim(numchar(i))//" of block "//trim(numchar(n))//" in x: stretching ratio rsx seems relatively elevated...",0)
             if (i.gt.1) then
                if (bl(n)%nrx1(i).lt.bl(n)%nrx2(i-1)) call mpistop("Beginning index of streching num."//trim(numchar(i))//" ("//trim(numchar(bl(n)%nrx1(i)))//&
                                                       ")  must not be lower than end index of streching num."//trim(numchar(i-1))//" ("//trim(numchar(bl(n)%nrx2(i-1)))//")",0)
             endif
          enddo
       endif
    enddo
    read(30,*) ! Stretching in direction j: nstretchy,  nrj1,nrj2,rj (X nstretchy)
    do n=1,nbloc
       read(30,*) nstretchy(n)
       ! If stretching specified in the y direction
       if (nstretchy(n).gt..0) then
          rewind(30)
          do i=1,40+2*nbloc+n !number of lines before this one
             read(30,*) !
          enddo
          ! Allocation of temporary array
          allocate(temp_array(3,nstretchy(n)))
          ! Read line
          read(30,*) nstretchy(n), temp_array
          ! Rearrange array with new ones
          allocate(bl(n)%nry1(nstretchy(n)),bl(n)%nry2(nstretchy(n)),bl(n)%rsy(nstretchy(n)))
          do i=1,nstretchy(n)
             bl(n)%nry1(i)=temp_array(1,i); bl(n)%nry2(i)=temp_array(2,i); bl(n)%rsy(i)=temp_array(3,i)
          enddo
          deallocate(temp_array)
          ! Check
          do i=1,nstretchy(n)
             if (bl(n)%nry2(i).le.bl(n)%nry1(i)) call mpistop("Stretching number "//trim(numchar(i))//" of block "//trim(numchar(n))//" in x: nrx2 must be greater than nry1",0)
             if (bl(n)%rsy(i).le.0.0_wp) call mpistop("Stretching number "//trim(numchar(i))//" of block "//trim(numchar(n))//" in x: stretching ratio rsy must be > 0",0)
             if (bl(n)%rsy(i).gt.2.0_wp) call mpistop("Stretching number "//trim(numchar(i))//" of block "//trim(numchar(n))//" in x: stretching ratio rsy seems relatively elevated...",0)
             if (i.gt.1) then
                if (bl(n)%nry1(i).lt.bl(n)%nry2(i-1)) call mpistop("Beginning index of streching num."//trim(numchar(i))//" ("//trim(numchar(bl(n)%nry1(i)))//&
                                                       ")  must not be lower than end index of streching num."//trim(numchar(i-1))//" ("//trim(numchar(bl(n)%nry2(i-1)))//")",0)
             endif
          enddo
       endif
    enddo
    close(30)

  end subroutine read_param_grid

end module mod_block
