!==============================================================================
module mod_pp_sp_write
!==============================================================================
  !> Module to compute autospectrum
!==============================================================================
  use mod_pp_mpi
  use mod_pp_var
  use mod_pp_sp_modal ! <- for modal analysis
  implicit none
  ! ---------------------------------------------------------------------------
  ! file units
  integer :: uid,uid_sl(3),uid_m
  ! averaging
  real(wp), dimension(:,:), allocatable :: var_av
  ! slices
  real(wp), dimension(:,:,:), allocatable :: E_sl1,E_sl2,E_sl3
  real(wp), dimension(:,:), allocatable :: E_mod
  ! ---------------------------------------------------------------------------

contains

  !============================================================================
  subroutine write_spectrum
  !============================================================================
    !> Write frequency-wavenumber spectra
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i1,i2,n
    ! -------------------------------------------------------------------------

    ! Open binary file and write header
    ! =================================
    if (iproc==0) call write_sp_header

    ! Reconstruction and writing of variance, skewness and kurtosis
    ! =============================================================
    ! if one direction is inhomogeneous
    if (i_in>0) then
       if (iproc==0) print *,'~> writing moments along inhomogeneous direction'
       call write_var_in(sqrt(rms2m_in))
       call write_var_in(skewm_in)
       call write_var_in(kurtm_in)
    endif
    
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
       
       ! Write slices
       ! ------------
       if (is_slice) then
          ! direction 1
          if (n_sl(1)>0) then
             do i1=1,n_sl(1)
                write(uid_sl(1)) ((E_sl1(i2,n,i1),i2=1,ng2),n=1,ng3)
             enddo
             close(uid_sl(1))
          endif

          ! direction 2
          if (n_sl(2)>0) then
             do i2=1,n_sl(2)
                write(uid_sl(2)) ((E_sl2(i1,n,i2),i1=1,ng1),n=1,ng3)
             enddo
             close(uid_sl(2))
          endif

          ! direction 3
          if (n_sl(3)>0) then
             do n=1,n_sl(3)
                write(uid_sl(3)) ((E_sl3(i1,i2,n),i1=1,ng1),i2=1,ng2)
             enddo
             close(uid_sl(3))
          endif
       endif

       ! Write modes
       ! -----------
       if (is_modal) then
          do n=1,nmod
             write(uid_m) (E_mod(i1,n),i1=1,ng1)
          enddo
          close(uid_m)
       endif
       
    endif
    
  end subroutine write_spectrum

  !============================================================================
  subroutine write_sp_header
  !============================================================================
    !> Write header for spectra output files
  !============================================================================
    use mod_grid  ! <- for dimensions
    use mod_utils ! <- for get_free_unit
    use mod_io_snapshots
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i,j,k,nd,cpt,nsl
    integer :: nlx,nly,nlz,nlt
    real(wp), dimension(:), allocatable :: kx,ky,kz,ff
    character(len=1) :: icap,isl
    character(len=2) :: itype
    character(len=2), dimension(3) :: ispd
    character(len=6) :: isp
    character(len=3) :: inpl
    character(len=10) :: ivar
    character(len=70) :: namefile,namefile_sl,namefile_m
    ! -------------------------------------------------------------------------

    ! Counter for uid number
    ! ======================
    cpt=0
    
    ! Define output file name
    ! =======================
    
    ! plane number -> character
    ! ------------
    write(inpl,FMT='(I3.3)') nplr
    
    ! is_capon -> character
    ! --------
    if (sp%is_capon) then
       icap='C'
    else
       icap='F'
    endif

    ! PSD directions
    ! --------------
    ispd=''
    do nd=1,sp%dim
       if (sp%d(nd)%name=='t') then
          ispd(nd)='f'
       else
          ispd(nd)='k'//sp%d(nd)%name
       endif
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
    case('uun') ! normal velocity
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
    if ((i_in>0).and.(snapshots(nsr)%type.eq.2)) then
       namefile=trim(dirRESU)//'E'//trim(ivar)//'_'//trim(isp)//'_'//name_in// &
                trim(itype)//trim(inpl)//'_bl'//trim(numchar(iblc_pp))//'_'//trim(icap)//'.bin'
    else
       namefile=trim(dirRESU)//'E'//trim(ivar)//'_'//trim(isp)// &
                trim(itype)//trim(inpl)//'_bl'//trim(numchar(iblc_pp))//'_'//trim(icap)//'.bin'
    endif

    print *,'~> writing spectrum file '//trim(namefile)

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
	     ! Temporary
             write(uid) yc(snapshots(nsr)%ind_i1,snapshots(nsr)%ind_j1)
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

    if (sp%dim==1) then
       
       ! write info for spectrum directions: only positive frequencies/wavenumbers
       ! ----------------------------------
       if (iproc==0) print *,'   - spectrum direction ',sp%d(1)%name
       select case(sp%d(1)%name)
       case('x')
          nlx=ngx/2
          write(uid) nlx
          allocate(kx(nlx))
          do k=1,nlx
             kx(k)=twopi*dble(k)/deltax/dble(ngx);
          enddo
          write(uid) (kx(k),k=1,nlx)
          deallocate(kx)
       case('y')
          nly=ngy/2
          write(uid) nly
          allocate(ky(nly))
          do k=1,nly
             ky(k)=twopi*dble(k)/deltay/dble(ngy);
          enddo
          write(uid) (ky(k),k=1,nly)
          deallocate(ky)
       case('z')
          nlz=ngz/2
          write(uid) nlz
          allocate(kz(nlz))
          do k=1,nlz
             kz(k)=twopi*dble(k)/deltaz/dble(ngz);
          enddo
          write(uid) (kz(k),k=1,nlz)
          deallocate(kz)
       case('t')
          nlt=nl3/2
          write(uid) nlt
          allocate(ff(nlt))
          do k=1,nlt
             ff(k)=dble(k)/dt_spec/dble(nl3);
          enddo
          write(uid) (ff(k),k=1,nlt)
          deallocate(ff)
       end select

       if (ndir_av>0) print *,'   - averaging direction ',dir_av(1)

    else

       ! write info for spectrum directions
       ! ----------------------------------
       do nd=1,sp%dim
          print *,'   - spectrum direction ',sp%d(nd)%name
          select case(sp%d(nd)%name)
          case('x')
             nlx=ngx/2
             write(uid) ngx
             allocate(kx(-nlx+1:nlx))
             do k=-nlx+1,nlx
                kx(k)=twopi*dble(k)/deltax/dble(ngx);
             enddo
             write(uid) (kx(k),k=-nlx+1,nlx)
          case('y')
             nly=ngy/2
             write(uid) ngy
             allocate(ky(-nly+1:nly))
             do k=-nly+1,nly
                ky(k)=twopi*dble(k)/deltay/dble(ngy);
             enddo
             write(uid) (ky(k),k=-nly+1,nly)
          case('z')
             nlz=ngz/2
             write(uid) ngz
             allocate(kz(-nlz+1:nlz))
             do k=-nlz+1,nlz
                kz(k)=twopi*dble(k)/deltaz/dble(ngz);
             enddo
             write(uid) (kz(k),k=-nlz+1,nlz)
          case('t')
             nlt=nl3/2
             write(uid) nl3
             allocate(ff(-nlt+1:nlt))
             do k=-nlt+1,nlt
                ff(k)=dble(k)/dt_spec/dble(nl3);
             enddo
             write(uid) (ff(k),k=-nlt+1,nlt)
          end select
       enddo

       if (ndir_av>0) print *,'   - averaging direction ',dir_av(1)

       ! Open binary file for slice outputs
       ! ==================================
       if (is_slice) then

          do nsl=1,3 ! * loop on slices *

             if (n_sl(nsl)>0) then

                ! slice number
                ! ------------
                write(isl,FMT='(I1.1)') nsl

                ! file name
                ! ---------
                if (i_in>0) then
                   namefile_sl=trim(dirRESU)//'E'//trim(ivar)//'_'//trim(isp)//'_slices'//trim(isl) &
                               //'_'//name_in//'_p'//trim(inpl)//'_'//trim(icap)//'.bin'
                else
                   namefile_sl=trim(dirRESU)//'E'//trim(ivar)//'_'//trim(isp)//'_slices'//trim(isl) &
                               //'_p'//trim(inpl)//'_'//trim(icap)//'.bin'
                endif

                ! file unit & open binary
                ! -----------------------
                uid_sl(nsl)=uid+cpt
                cpt=cpt+1
                open(uid_sl(nsl),file=trim(namefile_sl),form='unformatted',status='unknown')
                rewind(uid_sl(nsl))

                print *,'~> writing spectrum slices '//trim(namefile_sl)
                
                ! write info for spectrum directions
                ! ----------------------------------
                do nd=1,sp%dim
                   select case(sp%d(nd)%name)
                   case('x')
                      write(uid_sl(nsl)) ngx
                      write(uid_sl(nsl)) (kx(k),k=-nlx+1,nlx)
                   case('y')
                      write(uid_sl(nsl)) ngy
                      write(uid_sl(nsl)) (ky(k),k=-nly+1,nly)
                   case('z')
                      write(uid_sl(nsl)) ngz
                      write(uid_sl(nsl)) (kz(k),k=-nlz+1,nlz)
                   case('t')
                      write(uid_sl(nsl)) nl3
                      write(uid_sl(nsl)) (ff(k),k=-nlt+1,nlt)
                   end select
                enddo

                ! allocate slice array
                ! --------------------
                select case(nsl)
                case(1)
                   allocate(E_sl1(ng2,ng3,n_sl(nsl)))
                case(2)
                   allocate(E_sl2(ng1,ng3,n_sl(nsl)))
                case(3)
                   allocate(E_sl3(ng1,ng2,n_sl(nsl)))
                end select

             endif

          enddo ! * end loop on slices *
    
       endif

       ! Open binary file for modes [modal analysis]
       ! ===========================
       if (is_modal) then
          
          if (i_in>0) then
             namefile_m=trim(dirRESU)//'E'//trim(ivar)//'_'//trim(isp)//'_modes_'//name_in// &
                  '_p'//trim(inpl)//'_'//trim(icap)//'.bin'
          else
             namefile_m=trim(dirRESU)//'E'//trim(ivar)//'_'//trim(isp)//'_modes'// &
                  '_p'//trim(inpl)//'_'//trim(icap)//'.bin'
          endif
       
          uid_m=uid+cpt
          open(uid_m,file=trim(namefile_m),form='unformatted',status='unknown')
          rewind(uid_m)
          
          print *,'~> writing modal analysis '//trim(namefile_m)

          ! write info for non-homogeneous direction
          ! ----------------------------------------
          ! if one direction is inhomogeneous
          print *,'   - mode evolution in direction ',name_in
          select case(name_in)
          case('x')
             write(uid_m) ngx
             write(uid_m) (xg(i),i=1,ngx)
          case('y')
             write(uid_m) ngy
             write(uid_m) (yg(j),j=1,ngy)
          case('z')
             write(uid_m) ngz
             write(uid_m) (zg(k),k=1,ngz)
          case('t')
             write(uid_m) ngt
             write(uid_m) dt_spec
          end select

          ! write info for modes
          ! --------------------
          write(uid_m) nmod
          ! write om & beta as (0,0), (0,1), etc...
          write(uid_m) ((imode(i)%om-1)/n_period,i=1,nmod)
          write(uid_m) ((imode(i)%beta-1)/n_lambdaz,i=1,nmod)
          ! modify indices! (0,0) at (nlz,nlt)
          imode%om=nlt+imode%om-1
          imode%beta=nlz+imode%beta-1

          ! allocate mode array
          ! -------------------
          allocate(E_mod(ng_in,nmod))

       endif
       
    endif

    if (allocated(kx)) deallocate(kx)
    if (allocated(ky)) deallocate(ky)
    if (allocated(kz)) deallocate(kz)
    if (allocated(ff)) deallocate(ff)

  end subroutine write_sp_header

  !============================================================================
  subroutine write_sp_varg(tab,n)
  !============================================================================
    !> MPI recontruction on proc 0 and write
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer, intent(in) :: n ! time/freq index
    ! var to be reconstructed and written
    real(wp), dimension(nl1,nl2), intent(in) :: tab  
    ! -------------------------------------------------------------------------
    integer :: ip,i1,i2,m
    real(wp), dimension(ng1,ng2) :: varg
    ! -------------------------------------------------------------------------
    
    ! Reconstruction and writing
    ! ==========================
    if (iproc.ne.0) then
       
       ! everyone (except 0) send to proc 0
       ! ----------------------------------
       call MPI_SEND(tab,1,type_p,0,tag,COMM_global,info)
       
    else
       ! reconstruction on the global mesh
       ! ---------------------------------
       varg(1:nl1,1:nl2)=tab

       do ip=1,nproc-1
          call MPI_RECV(varg(coordx1(ip),coordx2(ip)),1,type_pg,ip,tag,COMM_global,status,info)
       enddo
       
       if (ndir_av>0) then
          ! averaging on proc 0
          ! -------------------
          select case(i_av(1))
          case(1)
             do i1=1,ng1
                do i2=1,ng2
                   var_av(i2,n)=var_av(i2,n)+varg(i1,i2)
                enddo
             enddo
          case(2)
             do i1=1,ng1
                do i2=1,ng2
                   var_av(i1,n)=var_av(i1,n)+varg(i1,i2)
                enddo
             enddo
          case(3)
             var_av=var_av+varg
          end select
       else
          ! write binary file
          ! -----------------
          write(uid) ((varg(i1,i2),i1=1,ng1),i2=1,ng2)
       endif
       
       ! store slices in direction 1
       ! ---------------------------
       if (n_sl(1)>0) then
          do m=1,n_sl(1)
             E_sl1(:,n,m)=varg(i1_sl(m),:)
          enddo
       endif

       ! store slices in direction 2
       ! ---------------------------
       if (n_sl(2)>0) then
          do m=1,n_sl(2)
             E_sl2(:,n,m)=varg(:,i2_sl(m))
          enddo
       endif

       ! store slices in direction 3
       ! ---------------------------
       if (n_sl(3)>0) then
          do m=1,n_sl(3)
             if (n==i3_sl(m)) E_sl3(:,:,m)=varg
          enddo
       endif

       ! store modes (modal analysis)
       ! ----------------------------
       if (is_modal) then
          do m=1,nmod
             if (n==imode(m)%om) E_mod(:,m)=varg(:,imode(m)%beta)
          enddo
       endif

    endif

  end subroutine write_sp_varg

  !============================================================================
  subroutine write_var_in(var)
  !============================================================================
    !> MPI reconstruction in inhomogeneous direction & write 
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i,ip
    real(wp), dimension(n_in) :: var
    real(wp), dimension(ng_in) :: varg
    ! -------------------------------------------------------------------------

    ! MPI reconstruction & write on proc 0
    ! ============================
    
    if (iproc.ne.0) then
       if (is_in(iproc)) &
            call MPI_SEND(var,1,type_in,0,tag,COMM_global,info)
    else
       ! reconstruction on the whole domain
       if (is_in(iproc)) varg(1:n_in)=var
       do ip=1,nproc-1
          if (is_in(ip)) &
               call MPI_RECV(varg(coordin(ip)),1,typeg_in,ip,tag,COMM_global,status,info)
       enddo

       ! Write moments
       ! =============
       write(uid) (varg(i),i=1,ng_in)
      
    endif

  end subroutine write_var_in
     
end module mod_pp_sp_write
