!===============================================================
module warnstop
!===============================================================
  !> author: Luca Sciacovelli
  !> date: 2018
  !> make serial or MPI stop or warning
  ! 
  ! This module is used to make serialstop, mpistop, serialwarn
  ! and mpiwarn independent of other modules.
!===============================================================
  implicit none
  ! ------------------------------------------------------------
  integer, private :: nwarn=0
  integer, private :: maxwarning=100
  ! ------------------------------------------------------------

contains

  !===============================================================
  subroutine mpistop(err_msg,flag)
  !===============================================================
    !> author: Luca Sciacovelli
    !> date: 2018
    !> this subroutine stop the mpi processes
    ! 
    ! err_msg is a string that explain what went wrong.
    ! no error message and flag == 0 --> normal shut down
    ! error message and flag == 0 --> controlled shut down
    ! flag /= 0 --> catastrophic shut down
  !===============================================================
    use mpi
    implicit none
    ! ------------------------------------------------------------
    ! Input variables
    character(len=*), intent(in) :: err_msg
    integer, intent(in) :: flag
    ! ------------------------------------------------------------
    ! Local variables
    character(len=10) :: s
    character(len=3) :: sflag
    character(len=4) :: smyid
    integer :: myid,ierror,masterid
    integer :: nwarntot
    ! ------------------------------------------------------------

    masterid=0

    ! whoami
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierror)

    nwarntot=0
    ! reduce nwarn
    if (flag==0) then
       call MPI_REDUCE(nwarn,nwarntot,1,mpi_integer,mpi_sum,masterid &
                      ,MPI_COMM_WORLD,ierror )
    endif

    if ((myid == 0).and.(nwarntot>0)) then
       write(s,'(I8)') nwarntot
       if (nwarntot==1) then
          call mpiwarn('There was '//trim( s )//' warning.',0)
       else
          call mpiwarn('There were '//trim( s )//' warnings.',0)
       endif
    endif

    if ((myid==0).and.(flag==0)) then
       write(*,'(A)') repeat("-",80)
       if ( len( trim(err_msg) ) > 1 ) then
          write(*,*) 'ERROR!    '//trim( err_msg )
       else
          write(*,*) "You got lucky, or at least program didn't crash... time for a drink!"
       endif
       write(*,'(A)') repeat( "-", 80 )
    endif

    if (flag/=0) then
       write(sflag,'(I3)') flag
       write(smyid,'(I4)') myid
       write(*,*) 'ERROR!    '//trim( err_msg )//' (flag ='//sflag//', myid ='//smyid//')'
       call MPI_ABORT(MPI_COMM_WORLD,flag,ierror)
    else
       call MPI_FINALIZE(ierror)
    endif

    stop

  end subroutine mpistop

  !===============================================================
  subroutine mpiwarn(err_msg,idflag)
  !===============================================================
    !> author: Luca Sciacovelli
    !> date: 2018
    !> this subroutine issues warnings in a uniform fashion
    ! 
    ! err_msg is a string that explain what went wrong.
    ! if idflag >= 0 only idflag process prints, otherwise all processes print
    ! nwarn is a global counter for the warnings.
    ! In oreder to avoid huge log files when there is a recurrent warning on large
    ! simulations, warnings are printed only if the number of warnings is less then
    ! maxwarning. The counter will continue to keep track of the warnings.
  !===============================================================
    use mpi
    implicit none
    ! ------------------------------------------------------------
    ! Inputs/Outputs
    character(len=*), intent(in) :: err_msg
    integer, intent(in) :: idflag
    ! ------------------------------------------------------------
    ! Local variables
    character(len=5) :: idtag
    integer :: myid,ierror
    ! ------------------------------------------------------------

    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierror)

    if (idflag<0) then
       ! print all
       if (nwarn<maxwarning) then
          write(idtag,'(I4)') myid
          write(*,'(1X,A)') 'WARNING!  '//trim(err_msg)//' (process #'//trim(idtag)//')'
       elseif(nwarn==maxwarning) then
          write(idtag,'(I4)') myid
          write(*,'(1X,A)') 'WARNING!  Max number of warnings reached. Going to silent mode.' &
               //' (process #'//trim( idtag )//')'
       endif
       nwarn = nwarn + 1
    else
       if (myid==0) then
          if (nwarn<maxwarning) then
             write(idtag,'(I4)') idflag
             write(*,'(1X,A)') 'WARNING!  '//trim(err_msg)//' (process #'//trim(idtag)//')'
          elseif(nwarn==maxwarning) then
             write(idtag,'(I4)') idflag
             write(*,'(1X,A)') 'WARNING!  Max number of warnings reached. Going to silent mode.' &
                             //' (process #'//trim(idtag)//')'
          endif
          nwarn = nwarn + 1
       endif
    endif

  end subroutine mpiwarn

end module warnstop
