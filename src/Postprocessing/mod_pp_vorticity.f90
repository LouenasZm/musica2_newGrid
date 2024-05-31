!==============================================================================
module mod_pp_vorticity
!==============================================================================
  !> Post-processing module for vorticity computation
!==============================================================================
  use mod_mpi
  use mod_flow
  use mod_constant
  use mod_io
  use warnstop
  implicit none
  ! ----------------------------------------------------------------------------
  real(wp) :: ss11,ss12,ss13,ss22,ss23,ss33
  real(wp) :: oo12,oo13,oo23
  real(wp), dimension(:,:,:), allocatable :: vort,QQ,lambda2
  ! ----------------------------------------------------------------------------

contains

  !===============================================================================
  subroutine calc_vort_crit
  !===============================================================================
    !> Compute vorticity and vortex identification criteria
  !===============================================================================
    use mod_interface
    use mod_utils
    implicit none
    ! ----------------------------------------------------------------------------

    ! Compute velocity gradients
    ! ==========================
    if (iproc==0) print *,'Compute velocity gradients ...'

    call grad_vel

    ! Compute vorticity
    ! =================
    allocate(vort(nx,ny,nz))
    call compute_vorticity3d

    ! Compute vorticity criteria
    ! ==========================
    if (iproc==0) print *,'Compute vorticity criteria ...'

    ! allocations
    ! -----------
    allocate(QQ(nx,ny,nz),lambda2(nx,ny,nz))

    ! Compute Q-criterion
    ! ===================
    call compute_Q_criterion

    !QQ=sqrt(QQ)*2.

    ! Compute lambda2-criterion
    ! =========================
    if (iproc==0) print *,'Diagonalization'
    call compute_l2_criterion

    ! Write vorticity file
    ! ====================
    if (iproc==0) print *,'write vorticity binary files ...'

    call write_vorticity('vort_bl'//trim(numchar(nob(iproc)))//filext_write)

    ! Free memory
    ! ===========
    deallocate(vort,QQ,lambda2)

  end subroutine calc_vort_crit


  !===============================================================================
  subroutine compute_vorticity3d
  !===============================================================================
    !> Compute 3D vorticity norm
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------

    vort=(dvx(1:nx,1:ny,1:nz)-duy(1:nx,1:ny,1:nz))**2
    vort=(dwx(1:nx,1:ny,1:nz)-duz(1:nx,1:ny,1:nz))**2+vort
    vort=(dwy(1:nx,1:ny,1:nz)-dvz(1:nx,1:ny,1:nz))**2+vort
    vort=sqrt(vort)

  end subroutine compute_vorticity3d

  !===============================================================================
  subroutine compute_Q_criterion
  !===============================================================================
    !> Compute Q-criterion
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k
    ! ----------------------------------------------------------------------------

    ! Q criterion:  Q=(Oij*Oij - Sij*Sij)/2.
    ! ------------
    do i=1,nx
       do j=1,ny
          do k=1,nz

             ss11 = dux(i,j,k)
             ss22 = dvy(i,j,k)
             ss33 = dwz(i,j,k)
             ss12 = 0.5_wp*(duy(i,j,k)+dvx(i,j,k))
             ss13 = 0.5_wp*(duz(i,j,k)+dwx(i,j,k))
             ss23 = 0.5_wp*(dvz(i,j,k)+dwy(i,j,k))

             oo12 = 0.5_wp*(duy(i,j,k)-dvx(i,j,k))
             oo13 = 0.5_wp*(duz(i,j,k)-dwx(i,j,k))
             oo23 = 0.5_wp*(dvz(i,j,k)-dwy(i,j,k))

             QQ(i,j,k)=oo12**2+oo13**2+oo23**2-(0.5_wp*(ss11**2+ss22**2+ss33**2)+ss12**2+ss13**2+ss23**2)
             !QQ(i,j,k)=oo12**2+oo13**2+oo23**2

          enddo
       enddo
    enddo

  end subroutine compute_Q_criterion

  !===============================================================================
  subroutine compute_l2_criterion
  !===============================================================================
    !> Compute lambda2 criterion
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k,l,m
    ! variables for diagonalization
    real, dimension(3,3) :: O,S,SO
    integer :: ldvl,ldvr
    real(wp), dimension(3) :: wi,wr
    real(wp), dimension(3,3) :: vr
    real(wp), dimension(1,3)  :: vl
    integer :: lwork
    real(wp), dimension(:), allocatable :: work
    real(wp) :: key
    ! ----------------------------------------------------------------------------

    ! initialization for diagonalization subroutine
    ! ---------------------------------------------
    ldvl=1
    ldvr=3
    lwork=5*3
    allocate(work(lwork))

    do i=1,nx
       if (mod(i,100)==0) print *,'index',i
       do j=1,ny
          do k=1,nz

             ss11 = dux(i,j,k)
             ss22 = dvy(i,j,k)
             ss33 = dwz(i,j,k)
             ss12 = 0.5_wp*(duy(i,j,k)+dvx(i,j,k))
             ss13 = 0.5_wp*(duz(i,j,k)+dwx(i,j,k))
             ss23 = 0.5_wp*(dvz(i,j,k)+dwy(i,j,k))

             oo12 = 0.5_wp*(duy(i,j,k)-dvx(i,j,k))
             oo13 = 0.5_wp*(duz(i,j,k)-dwx(i,j,k))
             oo23 = 0.5_wp*(dvz(i,j,k)-dwy(i,j,k))

             S(1,1)=ss11
             S(1,2)=ss12
             S(1,3)=ss13

             S(2,1)=ss12
             S(2,2)=ss22
             S(2,3)=ss23

             S(3,1)=ss13
             S(3,2)=ss23
             S(3,3)=ss33

             O(1,1)=0.
             O(1,2)=oo12
             O(1,3)=oo13

             O(2,1)=-oo12
             O(2,2)=0.
             O(2,3)=oo23

             O(3,1)=-oo13
             O(3,2)=-oo23
             O(3,3)=0.

             SO=matmul(S,S)+matmul(O,O)

             !print *,SO(1,1),SO(1,2),SO(1,3)
             !print *,SO(2,1),SO(2,2),SO(2,3)
             !print *,SO(3,1),SO(3,2),SO(3,3)

             call dgeev('n','v',3,SO,3,wr,wi,vl,ldvl,vr,ldvr,work,lwork,info)

             ! sort eigenvalues
             ! -----------------------
             !print *,'sort eigenvalues'

             do m=2,3
                key=wr(m)
                l=m-1
                !if (l==0) print *,'zero l',i,j,k,iproc
!!$                do while (l>0)
!!$                   if (l==0) print *,'zero l',i,j,k,iproc
!!$                   do while (wr(l)<key)
!!$                      wr(l+1)=wr(l)
!!$                      l=l-1
!!$                   enddo
!!$                enddo
                do while (wr(l)<key)
                   wr(l+1)=wr(l)
                   l=l-1
                   if (l==0) exit
                enddo
                wr(l+1)=key
             enddo

             !print *,'eigenvalues:'
             !print *,'------------'
             !do i=1,3
             !   print *,'lambda',i,': ',wr(i)
             !enddo

             lambda2(i,j,k)=wr(2)

          enddo
       enddo
    enddo

    deallocate(work)

  end subroutine compute_l2_criterion


  !==============================================================================================
  subroutine write_vorticity(namefile)
  !==============================================================================================
    !> write volume for vorticity
  !==============================================================================================
    use mod_utils
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Input/Output arguments
    character(len=*), intent(in) :: namefile
    ! -------------------------------------------------------------------------------------------
    ! Local variables
    integer  :: m,ndata,filetype
    real(wp) :: soltime
    integer :: i,j,k
    logical  :: iexist
    character(len=30) :: gridfile
    ! -------------------------------------------------------------------------------------------

    ! Write tecplot grid
    ! ==================

    ! Name of the grid file ! TO BE CHANGED ???
    ! =====================
    ! can be defined as optional input of the subroutine but never used at the present time
    gridfile='./grid_bl'//trim(numchar(nob(iproc)))//'.plt'
    
    ! check if grid file already exists
    ! ---------------------------------
    inquire(file=gridfile,exist=iexist)

    if (is_IOtec_write) then
       ! if not: write GRID
       ! ------------------
       if (.not.iexist) then

          ! prepare number and name of data: x,y,z coordinates
          ndata=3
          allocate(dataname(ndata),varlist(ndata))
          dataname(1:ndata) = (/'X','Y','Z'/)

          ! allocate and initialize dummy variable for datas
          allocate(dummy(nx,ny,nz,ndata))
          dummy = 0.0_wp

          ! fill data x,y,z in dummy array (Att. CARTESIAN only) TO BE CHANGED for curvilinear
          if (is_curv) then
             do i=1,nx
                do j=1,ny
                   dummy(i,j,:,1)=xc(i,j)
                   dummy(i,j,:,2)=yc(i,j)
                enddo
             enddo
          else
             do i=1,nx
                dummy(i,:,:,1)=x(i)
             enddo
             do j=1,ny
                dummy(:,j,:,2)=y(j)
             enddo
          endif
          do k=1,nz
             dummy(:,:,k,3)=z(k)
          enddo

          ! fill varlist (pointer on datas)
          do m=1,ndata
             varlist(m)%data=>dummy(:,:,:,m)
             varlist(m)%name= dataname(m)
          enddo

          ! filetype [0:full; 1:only grid; 2:only solution]
          filetype=1
          ! non-dimensional time
          soltime=tstar

          ! write grid block in separate file
          if (iproc.eq.0) write(*,*) 'Writing gridfile ~>',trim(gridfile)
          call write_tec(gridfile,filetype,soltime,field(nob(iproc)))

          ! free temporary data
          deallocate(dummy,varlist,dataname)

       endif ! end if test iexist
    endif

    ! Number & name of datas
    ! ======================
    ndata=4
    allocate(varlist(ndata),dataname(ndata),dummy(nx,ny,nz,ndata))
    dummy=0.0_wp
    dataname(1:ndata) = (/'vort','Qcri','lamb','uadi'/)

    ! fill varlist (pointer on datas)
    ! ===============================
    do m=1,ndata
       varlist(m)%data=> dummy(:,:,:,m)
       varlist(m)%name= dataname(m)
    enddo

    ! Write volume for metrics commutation
    ! ====================================

    ! fill data: double derivatives of metrics
    dummy(:,:,:,1) = vort(1:nx,1:ny,1:nz)
    dummy(:,:,:,2) = QQ(1:nx,1:ny,1:nz)
    dummy(:,:,:,3) = lambda2(1:nx,1:ny,1:nz)
    dummy(:,:,:,4) = uu(1:nx,1:ny,1:nz)/u_ref

    ! filetype [0:full; 1:only grid; 2:only solution]
    filetype=2
    ! non-dimensional time
    soltime=tstar

    ! write double derivatives of metrics
    if (iproc.eq.0) write(*,*) 'Writing vorticity ~>', trim(namefile)
    call write_tec(namefile,filetype,soltime,field(nob(iproc)))

    ! free pointer varlist, dataname, dummy
    ! =====================================
    deallocate(varlist,dataname,dummy)

  end subroutine write_vorticity

  !==============================================================================================
  subroutine write_ww(namefile)
  !==============================================================================================
    !> write volume for vorticity
  !==============================================================================================
    use mod_utils
    implicit none
    ! -------------------------------------------------------------------------------------------
    ! Input/Output arguments
    character(len=*), intent(in) :: namefile
    ! -------------------------------------------------------------------------------------------
    ! Local variables
    integer  :: m,ndata,filetype
    real(wp) :: soltime
    ! -------------------------------------------------------------------------------------------

    ! ww fluctuations
    ! ---------------
    if (iproc==0) print *,'ww fluctuations ...'
    ! Number & name of datas
    ! ======================
    ndata=1
    allocate(varlist(ndata),dataname(ndata),dummy(nx,ny,nz,ndata))
    dummy=0.0_wp
    dataname(1:ndata) = (/''/)

    ! fill varlist (pointer on datas)
    ! ===============================
    do m=1,ndata
       varlist(m)%data=> dummy(:,:,:,m)
       varlist(m)%name= dataname(m)
    enddo

    ! Write volume for metrics commutation
    ! ====================================

    ! fill data: double derivatives of metrics
    dummy(:,:,:,1) = ww(1:nx,1:ny,1:nz)

    ! filetype [0:full; 1:only grid; 2:only solution]
    filetype=2
    ! non-dimensional time
    soltime=tstar

    ! write double derivatives of metrics
    if (iproc.eq.0) write(*,*) 'Writing ww ~>', trim(namefile)
    call write_tec(namefile,filetype,soltime,field(nob(iproc)))

    ! free pointer varlist, dataname, dummy
    ! =====================================
    deallocate(varlist,dataname,dummy)

  end subroutine write_ww


end module mod_pp_vorticity
