!=================================================================
! Read/Write info (files info.ini and time.ini)
! Author: L. Sciacovelli (lucasciacovelli@gmail.com)
!=================================================================
subroutine read_write_info(operation)
  use mod_io
  use mod_mpi
  use mod_block
  use mod_flow
  use mod_constant
  use mod_time
  use mod_forcing_bulk
  implicit none
  ! --------------------------------------------------------------
  ! Input/Output arguments
  integer, intent(in) :: operation
  ! --------------------------------------------------------------
  ! Local variables
  integer :: n,i
  logical :: iexist
  integer :: uni_io
  integer, dimension(nbloc) :: ni_in,nj_in,nk_in
  character(len=120) :: line
  character(len=4) :: n_bl
  ! --------------------------------------------------------------
  
  if (operation.eq.WRITE) then
     
     if (iproc.eq.0) then
        uni_io = 99

        ! Write file info.ini
        ! ===================
        open(unit=uni_io, file=trim(dirDATA)//'info.ini',status='replace',form='formatted')
        write(uni_io,'(A,1X,I5,1X,L1)')   'nbloc & is_curv   =', nbloc, is_curv
        do n=1,nbloc
           write(n_bl,'(I4.4)') n
           write(uni_io,'(A,3(1X,I5))')   'NX NY NZ bl'//n_bl//'   =', bl(n)%ni,bl(n)%nj,bl(n)%nk
        enddo
        write(uni_io,'(A,2(1X,1E19.12))') 'Etot0 mgtot0      =', ektot0, mgtot0
        write(uni_io,'(A,3(1X,1E19.12))') 'xmin ymin zmin    =', xmin, ymin, zmin
        write(uni_io,'(A,3(1X,1E19.12))') 'xmax ymax zmax    =', xmax, ymax, zmax
        write(uni_io,'(A,2(1X,1E19.12))') 'Mref Reref        =', Mach, Re_ref
        write(uni_io,'(A,2(1X,1E19.12))') 'Mupref Muref      =', muw_ref, mu_ref
        write(uni_io,'(A,3(1X,1E19.12))') 'Roref Pref Tref   =', rho_ref, p_ref, T_ref
        write(uni_io,'(A,3(1X,1E19.12))') 'Uref cref Tscale  =', u_ref, c_ref, tscale
        write(uni_io,'(A,1X,1E19.12)')    'time step deltat  =', deltat
        if (is_forcing_bulk) &
        write(uni_io,'(A,1(1X,1E19.12))') 'Bulk forcing      =', forc_rhou
        close(unit=uni_io)

        ! Write file time.ini
        ! ===================
        inquire( file=trim(dirDATA)//'time.ini', exist=iexist )
        if (iexist) then
           open( unit=uni_io, file=trim(dirDATA)//'time.ini', status='old', position='append', &
                form='formatted' )
        else
           open( unit=uni_io, file=trim(dirDATA)//'time.ini', status='replace', form='formatted' )
           write(uni_io,'(A,5X,A,6X,A,15X,A)') 'Timestamp', 'Ntotal', 'cputot', 'time'
        endif
        write(uni_io,'(2(A,1X),I8,2(1X,g0))') filestamp, '=', ntotal, cputot, time
        close(uni_io)

     endif

  elseif ((operation.eq.READ).or.(operation.eq.READ_DEB)) then

     if (iproc==0) print *,'From info.ini & time.ini:'
     if (iproc==0) print *,'-------------------------'
     
     ! Read file info.ini
     ! ==================
     inquire( file=trim(dirDATA)//'info.ini', exist=iexist )
     if (.not.iexist) then
        call mpistop('file info.ini does not exist.',0)
     else
        ! open file
        uni_io = 99+iproc
        open( unit=uni_io, file=trim(dirDATA)//'info.ini',status='old',form='formatted')
        do n=1,nbloc
           write(n_bl,'(I4.4)') n
           call readline(uni_io,'NX NY NZ bl'//n_bl,line)
           read(line,*) ni_in(n),nj_in(n),nk_in(n)
        enddo
        if (TGV.or.CHIT) then
           call readline(uni_io,'Etot0 mgtot0',line)
           read(line,*) ektot0, mgtot0
           call readline(uni_io,'Uref cref Tscale',line)
           read(line,*) u_ref, c_ref, tscale
        endif
        ! read old time step
        call readline(uni_io,'time step deltat',line)
        read(line,*) deltat
        if (iproc==0) print *,'deltat read in info.ini: deltat=',deltat
        if (is_forcing_bulk) then
           call readline(uni_io,'Bulk forcing',line)
           read(line,*) forc_rhou
        endif
        close( unit=uni_io )
     endif
     ! check size for consistency
     if ((ngx/=ni_in(nob(iproc))).or.(ngy/=nj_in(nob(iproc))).or.(ngz/=nk_in(nob(iproc)))) then
        if (iproc.eq.0) then
           write(*,*) ngx, ni_in(nob(iproc))
           write(*,*) ngy, nj_in(nob(iproc))
           write(*,*) ngz, nk_in(nob(iproc))
        endif
        if ((ngx==ni_in(nob(iproc))).and.(ngy==nj_in(nob(iproc))).and.(nk_in(nob(iproc))==1)) then
           if (iproc==0) print *,'2D->3D extrusion of restart files ...'
           is_init_2D3D=.true.
        else
           call mpiwarn('ATTENTION! SIZE OF RESTART FILE NOT CONSISTENT!! Check.',0)
        endif
     endif

     if (operation.eq.READ_DEB) return ! Early call to read_write_info
     
     ! Read file time.ini
     ! ==================
     inquire( file=trim(dirDATA)//'time.ini', exist=iexist )

     if (.not.iexist) then
        call mpistop('file time.ini does not exist.', 0)
     else
        ! open file
        uni_io = 99+iproc
        if (trim(filestamp).eq.'0000_0000') then
           ! Filestamp not provided in input_card, we read the last line of time.ini
           open(unit=uni_io,file=trim(dirDATA)//'time.ini',status='old',form='formatted',position='append')
           backspace(uni_io)
           read(uni_io,'(A)') line
           filestamp = line(1:9)
           ! Don't use readline function -> issues with rewind
           do i = 1, len(line)
              if (line(i:i).eq.'=') then
                 line = adjustl(trim(line(i+1:)))
              endif
           enddo
        else
           ! Filestamp from input_card, we do not need to search the last filestamp
           open(unit=uni_io,file=trim(dirDATA)//'time.ini',status='old',form='formatted')
           call readline(uni_io, filestamp, line)
        endif
        read(line,*) ntotal, cputot, time
        ntotal_old = ntotal
        ! update tstar
        tstar = time / tscale
        close( unit=uni_io )
     endif
  endif

end subroutine read_write_info

!===============================================================
!  1) Read file From unit *, skipping commented lines
!  2) Search for values identified by the string 'key'
!  3) Return the line without the string key
!===============================================================
subroutine readline(uni,key,line)
  use mpi
  use warnstop
  implicit none
  ! ------------------------------------------------------------
  ! Input/Output arguments
  integer, intent(in) :: uni
  character(len=*), intent(in)  :: key
  character(len=*), intent(out) :: line
  ! Local variables
  integer :: eof,ikey,ierror,i,j,myid_local
  ! ------------------------------------------------------------
  !
  ! whoami
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid_local,ierror)

  ikey = len(key)
  eof = 0
  rewind(uni)
  do while (eof.eq.0)
     read(uni, '(A)',iostat=eof) line
     ! Handling errors
     if (eof.gt.0) then
        call mpistop('Check input file. Something went wrong.', 0)
        ! Handling end of file
     elseif (eof.lt.0) then
        if (myid_local.eq.0) then
           write(*,*) 'Key "', key, '" not found in input/ini file! Please add it!'
        endif
        call mpistop('End of file reached in input/ini file.', 0)
        ! Handling line
     else
        ! Check if line is not a comment
        if (line(1:1).ne.'!') then
           ! Check if line is the key of interest
           if (line(1:ikey).eq.key) then
              !if (myid_local.eq.0) write(*,*) trim(line)
              ! Search for "=" delimiter and remove left part
              do i = 1, len(line)
                 if (line(i:i).eq.'=') then
                    line = adjustl(trim(line(i+1:)))
                    ! Check if comments are present on line and remove them
                    do j = 1, len(line)
                       if (line(j:j).eq.'!') then
                          line = line(1:j-1)
                          return
                       endif
                    enddo
                    ! if no comments, return
                    return
                 endif
              enddo
              rewind(uni)
              return
           endif
        endif
     endif
  enddo

end subroutine readline

!===============================================================
! Subroutine for calculating filestamp, sequential string.
! format 0000.0000 ~> 0000_0000
!
! AUTHOR: L. Sciacovelli (lucasciacovelli@gmail.com)
!===============================================================
subroutine calcfilestamp(time,filestamp)
  use precision
  implicit none
  ! ------------------------------------------------------------
  real(wp), intent(in) :: time
  character(len=*), intent(out) :: filestamp
  ! ------------------------------------------------------------

  ! write( *,* ) time, int(100.*(time-int(time)))
  !!write( filestamp, '(I5.5,A,I3.3)' ) int(time), '_', int(1000._wp*(time-int(time)))
  write( filestamp, '(I4.4,A,I4.4)' ) int(time), '_', int(10000._wp*(time-int(time)))

end subroutine calcfilestamp
