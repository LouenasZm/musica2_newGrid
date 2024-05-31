!==============================================================================
module mod_pp_sp_read
!==============================================================================
  !> Module to read plane for post-processing (spectra)
!==============================================================================
  use mod_flow
  use mod_pp_mpi
  use mod_pp_var
  implicit none
  ! ---------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------

contains

  !============================================================================
  subroutine read_snapshot_part(nbl)
  !============================================================================
    !> Read plane for 1 time block and MPI partitioning
  !============================================================================
    use mod_time
    use mod_io_snapshots
    use mod_pp_sp_mean
    use mod_pp_emd ! for detrending
    use warnstop
    implicit none
    ! -------------------------------------------------------------------------
    integer, intent(in) :: nbl
    ! -------------------------------------------------------------------------
    integer :: i1,i2,i,k
    integer :: n,nn,n1,n2
    real(wp) :: time_rec,tstar_rec,n2_eta
    ! -------------------------------------------------------------------------
    ! from time.ini: 0034_9179 =  1000000 179819.5484027917 .2736598238224434E-01
    !time_ini=0.2736598238224434E-01 ! R=180 M=0.1
    ! from time.ini: 0076_4569 =   900000 161394.4001420772 .6657890969725914E-02
    !time_ini=0.6657890969725914E-02 ! R=180 M=0.3
    ! from time.ini: 0158_1577 =  1000000 180726.7783547812 .2529628623518157E-02
    !time_ini=0.2529628623518157E-02 ! R=180 M=0.7
    ! from time.ini: 0080_0948 =  1250000 252517.0771675049 .3030538347403331E-02
    !time_ini=0.3030538347403331E-02 ! R=360 M=0.7
    ! from time.ini: 0080_0948 =  1250000 252517.0771675049 .3030538347403331E-02
    !time_ini=0.3030538347403331E-02 ! R=360 M=0.7
    ! from time.ini: 0017_1108 =   900000 305982.1162062156 .8212071366630540E-02
    !time_ini=0.8212071366630540E-02 ! R=720L M=0.3
    ! from time.ini: 0035_1898 =   800000 185535.5538490151 .7249136876193624E-02
    !time_ini=0.7249136876193624E-02 ! R=360 M=0.3 irene/chan2
    ! from time.ini: 0039_8655 =   950000 345073.5371079335 .8212336081027335E-02
    !time_ini=0.8212336081027335E-02 ! R=360G3 M=0.3 irene/chan1
    ! from time.ini: 0019_9937 =   450000 344197.0191025645 .4118734918041008E-02
    ! from time.ini: 0019_9940 =   450010 344208.8135702535 .4118794178260440E-02
    !time_ini=0.4118734918041008E-02 ! R=360L M=0.3 irene/chan0

    ! Print screen
    ! ============
    if (iproc==0) print 103,nbl
103 format(1x,'reading block: ',i3)


    ! Define reading bounds n1 & n2 (with eventual overlapping)
    ! =============================
    if (sp%is_overlap) then
       if (nbl==1) then
          n1=1
          n2=nl3
       else
          var_r(:,:,1:nl3-sp%loverlap)=var0_r
          nn=(nbl-1)*sp%loverlap
          n1=nn+1+nl3-sp%loverlap;
          n2=nn+nl3;
       endif
    else
       n1=(nbl-1)*nl3+1
       n2=nbl*nl3
    endif
    if (iproc==0) print*,n1,n2

    if ((sp%varname(1).eq."uut").or.(sp%varname(1).eq."uun")) then
       ! if not done in grid_normals, done here
       ! --------------------------------------
       if (.not.allocated(nxn_jmin)) then
          allocate(nxn_jmin(1:nx,1:nz),nyn_jmin(1:nx,1:nz))
          do k=1,nz
             do i=1,nx
                n2_eta=y_ksi(i,1)**2+x_ksi(i,1)**2
                nxn_jmin(i,k) = -y_ksi(i,1)/sqrt(n2_eta)
                nyn_jmin(i,k) =  x_ksi(i,1)/sqrt(n2_eta)
             enddo
          enddo
       endif
    endif

    ! Read plane between samples n1 & n2
    ! ==================================
    do n=n1,n2

       ! time index (after overlapping)
       ! ----------
       if ((sp%is_overlap).and.(nbl>1)) then
          nn=n-(nbl-1)*sp%loverlap
       else
          nn=n-(nbl-1)*nl3
       endif

       ! determine filestamp [not useful if is_timestamp=.false.]
       ! -------------------
       time_rec=time_ini+n*dt_spec
       tstar_rec=time_rec/tscale
       call calcfilestamp(tstar_rec,filestamp)
       ! correct filestamp
       !if (n==20450) filestamp='0178_3846'
       !if ((iproc==0).and.(mod(n,1)==0)) print *,'   it:',n,nn,'filestamp:',filestamp
       ! if ((iproc==0).and.(mod(n,100)==0)) print *,'   it:',n,nn,'filestamp:',filestamp
       if ((iproc==0).and.(mod(n,n2/10)==0)) print *,'   it:',n,nn,'filestamp:',filestamp

       ! MPI-IO read
       ! -----------
       call read_snapshot(nsr,iblc_pp,dirDATA)
       
       ! fill variable in var_r
       ! ----------------------

       ! reading planes
       ! --------------
       if (snapshots(nsr)%type.eq.2) then
          select case(snapshots(nsr)%normal)

          case (1) ! plane YZ (with normal 1 along X)
             select case (trim(sp%varname(1)))
             case('prs') ! pressure
                do i2=1,nl2
                   do i1=1,nl1
                      var_r(i1,i2,nn)=prs(snapshots(nsr)%index,i1,i2)
                   enddo
                enddo
             case( 'uu') ! streamwise velocity
                do i2=1,nl2
                   do i1=1,nl1
                      var_r(i1,i2,nn)=uu(snapshots(nsr)%index,i1,i2)
                   enddo
                enddo
             case( 'vv') ! crossflow velocity
                do i2=1,nl2
                   do i1=1,nl1
                      var_r(i1,i2,nn)=vv(snapshots(nsr)%index,i1,i2)
                   enddo
                enddo
             case( 'ww') ! spanwise velocity
                do i2=1,nl2
                   do i1=1,nl1
                      var_r(i1,i2,nn)=ww(snapshots(nsr)%index,i1,i2)
                   enddo
                enddo
             case( 'uut') ! tangential velocity
                do i2=1,nl2
                   do i1=1,nl1
                      var_r(i1,i2,nn)=uu(snapshots(nsr)%index,i1,i2)*nyn_jmin(snapshots(nsr)%index,i2) - vv(snapshots(nsr)%index,i1,i2)*nxn_jmin(snapshots(nsr)%index,i2)
                   enddo
                enddo
             case( 'uun') ! normal velocity
                do i2=1,nl2
                   do i1=1,nl1
                      var_r(i1,i2,nn)=uu(snapshots(nsr)%index,i1,i2)*nxn_jmin(snapshots(nsr)%index,i2) + vv(snapshots(nsr)%index,i1,i2)*nyn_jmin(snapshots(nsr)%index,i2)
                   enddo
                enddo
             case('rho') ! density
                call mpistop('NOT IMPLEMENTED YET',0)
             case('rhou') ! density*u
                if (allocated(rhou)) then
                   do i2=1,nl2
                      do i1=1,nl1
                         var_r(i1,i2,nn)=rhou(snapshots(nsr)%index,i1,i2)
                      enddo
                   enddo
                else if (allocated(rho).and.allocated(uu)) then
                   do i2=1,nl2
                      do i1=1,nl1
                         var_r(i1,i2,nn)=rho(snapshots(nsr)%index,i1,i2)*uu(snapshots(nsr)%index,i1,i2)
                      enddo
                   enddo
                else
                   call mpistop('rhou cannot be constructed from var in snapshot being PP',0)
                endif
             case('Tmp') ! temperature
                call mpistop('NOT IMPLEMENTED YET',0)
             case('Frhov') ! Frhov
                call mpistop('NOT IMPLEMENTED YET',0)
             case('Grhow') ! Grhow
                call mpistop('NOT IMPLEMENTED YET',0)
             case('udf') ! user-defined variable
                call mpistop('NOT IMPLEMENTED YET',0)
             case default
                call mpistop('variable '//trim(sp%varname(1)) &
                     //' not listed in read_snapshot_part [mod_pp_sp_read.f90]',0)
             end select

          case (2) ! plane XZ (with normal 2 along Y)
             select case (trim(sp%varname(1)))
             case('prs') ! pressure
                do i2=1,nl2
                   do i1=1,nl1
                      var_r(i1,i2,nn)=prs(i1,snapshots(nsr)%index,i2)
                   enddo
                enddo
             case( 'uu') ! streamwise velocity
                do i2=1,nl2
                   do i1=1,nl1
                      var_r(i1,i2,nn)=uu(i1,snapshots(nsr)%index,i2)
                   enddo
                enddo
             case( 'vv') ! crossflow velocity
                do i2=1,nl2
                   do i1=1,nl1
                      var_r(i1,i2,nn)=vv(i1,snapshots(nsr)%index,i2)
                   enddo
                enddo
             case( 'ww') ! spanwise velocity
                do i2=1,nl2
                   do i1=1,nl1
                      var_r(i1,i2,nn)=ww(i1,snapshots(nsr)%index,i2)
                   enddo
                enddo
             case( 'uut') ! tangential velocity
                do i2=1,nl2
                   do i1=1,nl1
                      var_r(i1,i2,nn)=uu(i1,snapshots(nsr)%index,i2)*nyn_jmin(i1,i2) - vv(i1,snapshots(nsr)%index,i2)*nxn_jmin(i1,i2)
                   enddo
                enddo
             case( 'uun') ! normal velocity
                do i2=1,nl2
                   do i1=1,nl1
                      var_r(i1,i2,nn)=uu(i1,snapshots(nsr)%index,i2)*nxn_jmin(i1,i2) + vv(i1,snapshots(nsr)%index,i2)*nyn_jmin(i1,i2)
                   enddo
                enddo
             case('rho') ! density
                call mpistop('NOT IMPLEMENTED YET',0)
             case('Tmp') ! temperature
                call mpistop('NOT IMPLEMENTED YET',0)
             case('Frhov') ! Frhov
                call mpistop('NOT IMPLEMENTED YET',0)
             case('Grhow') ! Grhow
                call mpistop('NOT IMPLEMENTED YET',0)
             case('udf') ! user-defined variable
                call mpistop('NOT IMPLEMENTED YET',0)
             case default
                call mpistop('variable '//trim(sp%varname(1)) &
                     //' not listed in read_snapshot_part [mod_pp_sp_read.f90]',0)
             end select

          case (3) ! plane XY (with normal 3 along Z)
             select case (trim(sp%varname(1)))
             case('prs') ! pressure
                do i2=1,nl2
                   do i1=1,nl1
                      var_r(i1,i2,nn)=prs(i1,i2,snapshots(nsr)%index)
                   enddo
                enddo
             case( 'uu') ! streamwise velocity
                do i2=1,nl2
                   do i1=1,nl1
                      var_r(i1,i2,nn)=uu(i1,i2,snapshots(nsr)%index)
                   enddo
                enddo
             case( 'vv') ! crossflow velocity
                do i2=1,nl2
                   do i1=1,nl1
                      var_r(i1,i2,nn)=vv(i1,i2,snapshots(nsr)%index)
                   enddo
                enddo
             case( 'ww') ! spanwise velocity
                do i2=1,nl2
                   do i1=1,nl1
                      var_r(i1,i2,nn)=ww(i1,i2,snapshots(nsr)%index)
                   enddo
                enddo
             case( 'uut') ! tangential velocity
                do i2=1,nl2
                   do i1=1,nl1
                      var_r(i1,i2,nn)=uu(i1,i2,snapshots(nsr)%index)*nyn_jmin(i1,snapshots(nsr)%index) - vv(i1,i2,snapshots(nsr)%index)*nxn_jmin(i1,snapshots(nsr)%index)
                   enddo
                enddo
             case( 'uun') ! normal velocity
                do i2=1,nl2
                   do i1=1,nl1
                      var_r(i1,i2,nn)=uu(i1,i2,snapshots(nsr)%index)*nxn_jmin(i1,snapshots(nsr)%index) + vv(i1,i2,snapshots(nsr)%index)*nyn_jmin(i1,snapshots(nsr)%index)
                   enddo
                enddo
             case('rho') ! density
                call mpistop('NOT IMPLEMENTED YET',0)
             case('Tmp') ! temperature
                call mpistop('NOT IMPLEMENTED YET',0)
             case('Frhov') ! Frhov
                call mpistop('NOT IMPLEMENTED YET',0)
             case('Grhow') ! Grhow
                call mpistop('NOT IMPLEMENTED YET',0)
             case('udf') ! user-defined variable
                call mpistop('NOT IMPLEMENTED YET',0)
             case default
                call mpistop('variable '//trim(sp%varname(1)) &
                     //' not listed in read_snapshot_part [mod_pp_sp_read.f90]',0)
             end select

          case default

             call mpistop('error in snapshot number selection in param_pp.ini', 0)

          end select

       ! reading lines
       ! -------------
       else if (snapshots(nsr)%type.eq.1) then
          select case(snapshots(nsr)%dir)

          case (1) ! line X (normal to plane YZ)
             select case (trim(sp%varname(1)))
             case('prs') ! pressure
                do i1=1,nl1
                   var_r(i1,1,nn)=prs(i1,snapshots(nsr)%tectype%ny1,snapshots(nsr)%tectype%nz1)
                enddo
             case( 'uu') ! streamwise velocity
                do i1=1,nl1
                   var_r(i1,1,nn)=uu(i1,snapshots(nsr)%tectype%ny1,snapshots(nsr)%tectype%nz1)
                enddo
             case( 'vv') ! crossflow velocity
                do i1=1,nl1
                   var_r(i1,1,nn)=vv(i1,snapshots(nsr)%tectype%ny1,snapshots(nsr)%tectype%nz1)
                enddo
             case( 'ww') ! spanwise velocity
                do i1=1,nl1
                   var_r(i1,1,nn)=ww(i1,snapshots(nsr)%tectype%ny1,snapshots(nsr)%tectype%nz1)
                enddo
             case( 'uut') ! tangential velocity
                do i1=1,nl1
                   var_r(i1,1,nn)=uu(i1,snapshots(nsr)%tectype%ny1,snapshots(nsr)%tectype%nz1)*nyn_jmin(i1,snapshots(nsr)%tectype%nz1) - vv(i1,snapshots(nsr)%tectype%ny1,snapshots(nsr)%tectype%nz1)*nxn_jmin(i1,snapshots(nsr)%tectype%nz1)
                enddo
             case( 'uun') ! normal velocity
                do i1=1,nl1
                   var_r(i1,1,nn)=uu(i1,snapshots(nsr)%tectype%ny1,snapshots(nsr)%tectype%nz1)*nxn_jmin(i1,snapshots(nsr)%tectype%nz1) + vv(i1,snapshots(nsr)%tectype%ny1,snapshots(nsr)%tectype%nz1)*nyn_jmin(i1,snapshots(nsr)%tectype%nz1)
                enddo
             case('rho') ! density
                call mpistop('NOT IMPLEMENTED YET',0)
             case('Tmp') ! temperature
                call mpistop('NOT IMPLEMENTED YET',0)
             case('Frhov') ! Frhov
                call mpistop('NOT IMPLEMENTED YET',0)
             case('Grhow') ! Grhow
                call mpistop('NOT IMPLEMENTED YET',0)
             case('udf') ! user-defined variable
                call mpistop('NOT IMPLEMENTED YET',0)
             case default
                call mpistop('variable '//trim(sp%varname(1)) &
                     //' not listed in read_snapshot_part [mod_pp_sp_read.f90]',0)
             end select

          case (2) ! line Y (normal to plane XZ)
             select case (trim(sp%varname(1)))
             case('prs') ! pressure
                do i1=1,nl1
                   var_r(i1,1,nn)=prs(snapshots(nsr)%tectype%nx1,i1,snapshots(nsr)%tectype%nz1)
                enddo
             case( 'uu') ! streamwise velocity
                do i1=1,nl1
                   var_r(i1,1,nn)=uu(snapshots(nsr)%tectype%nx1,i1,snapshots(nsr)%tectype%nz1)
                enddo
             case( 'vv') ! crossflow velocity
                do i1=1,nl1
                   var_r(i1,1,nn)=vv(snapshots(nsr)%tectype%nx1,i1,snapshots(nsr)%tectype%nz1)
                enddo
             case( 'ww') ! spanwise velocity
                do i1=1,nl1
                   var_r(i1,1,nn)=ww(snapshots(nsr)%tectype%nx1,i1,snapshots(nsr)%tectype%nz1)
                enddo
             case( 'uut') ! tangential velocity
                do i1=1,nl1
                   var_r(i1,1,nn)=uu(snapshots(nsr)%tectype%nx1,i1,snapshots(nsr)%tectype%nz1)*nyn_jmin(snapshots(nsr)%tectype%nx1,snapshots(nsr)%tectype%nz1) - vv(snapshots(nsr)%tectype%nx1,i1,snapshots(nsr)%tectype%nz1)*nxn_jmin(snapshots(nsr)%tectype%nx1,snapshots(nsr)%tectype%nz1)
                enddo
             case( 'uun') ! normal velocity
                do i1=1,nl1
                   var_r(i1,1,nn)=uu(snapshots(nsr)%tectype%nx1,i1,snapshots(nsr)%tectype%nz1)*nxn_jmin(snapshots(nsr)%tectype%nx1,snapshots(nsr)%tectype%nz1) + vv(snapshots(nsr)%tectype%nx1,i1,snapshots(nsr)%tectype%nz1)*nyn_jmin(snapshots(nsr)%tectype%nx1,snapshots(nsr)%tectype%nz1)
                enddo
             case('rho') ! density
                call mpistop('NOT IMPLEMENTED YET',0)
             case('Tmp') ! temperature
                call mpistop('NOT IMPLEMENTED YET',0)
             case('Frhov') ! Frhov
                call mpistop('NOT IMPLEMENTED YET',0)
             case('Grhow') ! Grhow
                call mpistop('NOT IMPLEMENTED YET',0)
             case('udf') ! user-defined variable
                call mpistop('NOT IMPLEMENTED YET',0)
             case default
                call mpistop('variable '//trim(sp%varname(1)) &
                     //' not listed in read_snapshot_part [mod_pp_sp_read.f90]',0)
             end select

          case (3) ! line Z (normal to plane XY)
             select case (trim(sp%varname(1)))
             case('prs') ! pressure
                do i2=1,nl2
                   var_r(1,i2,nn)=prs(snapshots(nsr)%tectype%nx1,snapshots(nsr)%tectype%ny1,i2)
                enddo
             case( 'uu') ! streamwise velocity
                do i2=1,nl2
                   var_r(1,i2,nn)=uu(snapshots(nsr)%tectype%nx1,snapshots(nsr)%tectype%ny1,i2)
                enddo
             case( 'vv') ! crossflow velocity
                do i2=1,nl2
                   var_r(1,i2,nn)=vv(snapshots(nsr)%tectype%nx1,snapshots(nsr)%tectype%ny1,i2)
                enddo
             case( 'ww') ! spanwise velocity
                do i2=1,nl2
                   var_r(1,i2,nn)=ww(snapshots(nsr)%tectype%nx1,snapshots(nsr)%tectype%ny1,i2)
                enddo
             case( 'uut') ! tangential velocity
                do i2=1,nl2
                   var_r(1,i2,nn)=uu(snapshots(nsr)%tectype%nx1,snapshots(nsr)%tectype%ny1,i2)*nyn_jmin(snapshots(nsr)%tectype%nx1,i2) - vv(snapshots(nsr)%tectype%nx1,snapshots(nsr)%tectype%ny1,i2)*nxn_jmin(snapshots(nsr)%tectype%nx1,i2)
                enddo
             case( 'uun') ! normal velocity
                do i2=1,nl2
                   var_r(1,i2,nn)=uu(snapshots(nsr)%tectype%nx1,snapshots(nsr)%tectype%ny1,i2)*nxn_jmin(snapshots(nsr)%tectype%nx1,i2) + vv(snapshots(nsr)%tectype%nx1,snapshots(nsr)%tectype%ny1,i2)*nyn_jmin(snapshots(nsr)%tectype%nx1,i2)
                enddo
             case('rho') ! density
                call mpistop('NOT IMPLEMENTED YET',0)
             case('Tmp') ! temperature
                call mpistop('NOT IMPLEMENTED YET',0)
             case('Frhov') ! Frhov
                call mpistop('NOT IMPLEMENTED YET',0)
             case('Grhow') ! Grhow
                call mpistop('NOT IMPLEMENTED YET',0)
             case('udf') ! user-defined variable
                call mpistop('NOT IMPLEMENTED YET',0)
             case default
                call mpistop('variable '//trim(sp%varname(1)) &
                     //' not listed in read_snapshot_part [mod_pp_sp_read.f90]',0)
             end select

          case default

             call mpistop('error in snapshot number selection in param_pp.ini', 0)

          end select

       endif
       
       call MPI_BARRIER(COMM_global,info)
       
    enddo ! * end loop on samples *

    ! Store previous field (for overlapping)
    ! ====================    
    if (sp%is_overlap) var0_r=var_r(:,:,sp%loverlap+1:nl3)

    ! Print screen
    ! ============
    if (iproc==0) print 104,nbl
104 format(1x,'block ',i3,' read')

    ! Detrending the signal: using EMD decomposition 
    ! ======================
    !if (iproc==0) print *,'detrend with EMD ...'
    !do i=1,nx
    !   if ((iproc==0).and.(mod(i,10)==0)) print *,'i',i
    !   do k=1,nz
    !      call emd_filter(p1(i,k,:))
    !   enddo
    !enddo

    ! Centering the signal: subtraction of average + compute moments
    ! =====================
    if (sp%dim==1) call subtract_mean_1var
    if (sp%dim==2) call subtract_mean_2var
    if (sp%dim==3) call subtract_mean_3var
 
    ! Initialize complex array
    ! ========================
    if (sp%dim>1) var_c=var_r

  end subroutine read_snapshot_part

end module mod_pp_sp_read
