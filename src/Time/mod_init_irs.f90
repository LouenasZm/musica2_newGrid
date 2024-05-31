!=================================================================================
module mod_init_irs
!=================================================================================
  !> author: XG,AB
  !> date: 2021-2023
  !> Module for Implicit Residual Smoothing, parallelized with ghost cells
  !> (variable ngh per direction only for one-sided comms)
!=================================================================================
  use precision
  use warnstop
  implicit none
  ! ------------------------------------------------------------------------------
  ! ------------------------------------------------------------------------------

contains

  !===============================================================================
  subroutine init_irs
  !===============================================================================
    !> Initialisation of Implicit Residual Smoothing
  !===============================================================================
    use mod_flow          ! <- for ngh_irs,dimensions,increments
    use mod_time          ! <- for IRS parameters
    use mod_mpi_part      ! <- for iproc
    use mod_mpi_types_one_sided
    use mod_mpi_types_two_sided
    use mod_comm1
    use mod_grid_directions
    use mod_constant ! <- is_RANS
    use warnstop
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k,max_nghirs
    ! ----------------------------------------------------------------------------

    if (iproc==0) then
       print *,repeat('=',70)
       print *,'Initialization of IRS...'
    endif

    if (iproc==0) then
       print *,'*******************************************************************************'
       print *,'current version (split per direction with only one-sided comm)'
       print *,'TEMPORARY two-sided but not optimized -> ngh_max everywhere'
       print *,'*******************************************************************************'
    endif
    
    ! Safety checks
    ! =============
    if (is_2d) is_irs_k=.false.

    ! Parallelisation with ghost points only implemented for IRS2 or IRS4
    if (type_para.ne.1) then
       call mpistop('Only ghost points parallelisation for IRS2-IRS4. Shutting down...',0)
    endif
    if ((iirs.ne.2).and.(iirs.ne.4)) then
       call mpistop('Only IRS2-IRS4 implemented yet. Shutting down...',0)
    endif

    ! Assign procedure pointers
    ! =========================
    call assign_irs
       
    ! Parallelisation with ghost points
    ! =================================
    if (iproc==0) print *,' Parallelisation of IRS with ghost points'

    ! Determination of ghost points for each face
    ! -------------------------------------------

    ! initialization
    ngh_irs=ngh
    nghirs=ngh

    ! compute max CFL in each faces
    call check_max_cfl_interface

    ! fill ngh_irs(i) corresponding to 2*CFL+1
    do i=1,dim
       do j=1,2
          k=2*(i-1)+j
          ngh_irs(k)=max(ngh,1*int(2*BC_face(i,j)%cflmax(i)+1)) ! 2*CFL_normal+1
          !print *,'iproc',iproc,'ngh_irs',k,':',ngh_irs(k)
       enddo
    enddo

    !!ngh_irs=nghirs

    ! Maximum number of ghost cells for IRS
    ! -------------------------------------
    max_nghirs=maxval(ngh_irs,1)
    call MPI_ALLREDUCE(MPI_IN_PLACE,max_nghirs,1,MPI_INTEGER,MPI_MAX,COMM_global,info)
    if (iproc==0) print *,' Maximum number of ghost points used for IRS : ', max_nghirs

    call MPI_BARRIER(COMM_global,info)

    nghirs=max_nghirs
    ngh_irs=nghirs

    call MPI_BARRIER(COMM_global,info)

    ! Communicate numbers of ghost points for each face's neighbor
    ! ------------------------------------------------------------

    ! initialization of ngh_irs in neighboring domain
    do i=1,3
       do j=1,2
          BC_face(i,j)%ngh_irs=[ngh,ngh,ngh,ngh,ngh,ngh]
       enddo
    enddo

    ! send ngh_irs to existing neighbors
    do i=1,dim
       do j=1,2
          k=2*(i-1)+j
          if (neighbor(k)>=0) call MPI_SEND(ngh_irs,6,MPI_INTEGER,neighbor(k),tag,COMM_global,info)
       enddo
    enddo

    ! receive ngh_irs from existing neighbors
    do i=1,dim
       do j=1,2
          k=2*(i-1)+j
          if (neighbor(k)>=0) call MPI_RECV(BC_face(i,j)%ngh_irs,6,MPI_INTEGER,neighbor(k),tag,COMM_global,status,info)
       enddo
    enddo

    ! Set local grid sizes extended to ghost points
    ! ---------------------------------------------
    nx1_irs= 1-ngh_irs(1)
    nx2_irs=nx+ngh_irs(2)
    ny1_irs= 1-ngh_irs(3)
    ny2_irs=ny+ngh_irs(4)
    nz1_irs= 1-ngh_irs(5)
    nz2_irs=nz+ngh_irs(6)

    if (is_2d) then
       nz1_irs=1
       nz2_irs=1
    endif

    ! Allocate increments arrays with extended ghost points
    ! -----------------------------------------------------
    allocate( Krho(nx1_irs:nx2_irs,ny1_irs:ny2_irs,nz1_irs:nz2_irs))
    allocate(Krhou(nx1_irs:nx2_irs,ny1_irs:ny2_irs,nz1_irs:nz2_irs))
    allocate(Krhov(nx1_irs:nx2_irs,ny1_irs:ny2_irs,nz1_irs:nz2_irs))
    allocate(Krhow(nx1_irs:nx2_irs,ny1_irs:ny2_irs,nz1_irs:nz2_irs))
    allocate(Krhoe(nx1_irs:nx2_irs,ny1_irs:ny2_irs,nz1_irs:nz2_irs))
    ! initialize
    Krho=0.0_wp
    Krhou=0.0_wp
    Krhov=0.0_wp
    Krhow=0.0_wp
    Krhoe=0.0_wp
    ! RANS
    if (is_RANS) then
       deallocate(Knutil)
       allocate(Knutil(nx1_irs:nx2_irs,ny1_irs:ny2_irs,nz1_irs:nz2_irs))
       Knutil=0.0_wp
    endif

    ! Allocate CFL arrays with extended ghost points
    ! ----------------------------------------------
    allocate(cfl_l(nx1_irs:nx2_irs,ny1_irs:ny2_irs,nz1_irs:nz2_irs))
    cfl_l=0.0_wp

    ! Definition of MPI types for communications
    ! ------------------------------------------
    !call mpi_types_comm_inc_old
    call mpi_types_comm(nghirs,'i')
    !call mpi_types_comm1_inc

    ! Define windows on allocated arrays (for one-sided comms)
    ! ----------------------------------
    !call mpi_win_comm1_inc
    !call mpi_win_comm1_cfl
 
  end subroutine init_irs

  !===============================================================================
  subroutine assign_irs
  !===============================================================================
    !> 2nd-order Implicit Residual Smoothing (IRS4) * per direction *
    !> Parallel inversion of pentadiagonal matrix using ghost points
  !===============================================================================
    use mod_ngh        ! <- for is_curv
    use mod_time       ! <- for is_irs
    use mod_interface  ! <- procedure pointer interfaces
    use mod_utils      ! <- procedure voidp
    use mod_irs_d      ! <- for irs2/4_ngh_i/j/k
    use mod_irs_d_rans ! <- for irs2/4_ngh_i/j/k RANS
    use mod_comm       ! <- for communication_inc
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------
    
    ! Assign directions
    ! =================
    if (is_irs_i) then
       if (iirs==2) apply_irs_1 => irs2_ngh_i
       if (iirs==4) apply_irs_1 => irs4_ngh_i
    else
       apply_irs_1 => voidp
    endif

    if (is_irs_j) then
       if (iirs==2) apply_irs_2 => irs2_ngh_j
       if (iirs==4) apply_irs_2 => irs4_ngh_j
    else
       apply_irs_2 => voidp
    endif

    if (is_irs_k) then
       if (iirs==2) apply_irs_3 => irs2_ngh_k
       if (iirs==4) apply_irs_3 => irs4_ngh_k
    else
       apply_irs_3 => voidp
    endif

    ! Main IRS routine
    ! ================
    irs => irs_routine
  
    ! communication routine for increments
    ! ====================================
    if (is_2D) then
       communication_inc => communication_inc2d
    else
       communication_inc => communication_inc3d
    endif
  
    ! RANS irs pointers
    ! =================
    if (is_RANS) then
       if (is_irs_i) then
          if (iirs==2) apply_irs_1_rans => irs2_ngh_i_rans
          if (iirs==4) apply_irs_1_rans => irs4_ngh_i_rans
       else
          apply_irs_1_rans => voidp
       endif
       if (is_irs_j) then
          if (iirs==2) apply_irs_2_rans => irs2_ngh_j_rans
          if (iirs==4) apply_irs_2_rans => irs4_ngh_j_rans
       else
          apply_irs_2_rans => voidp
       endif
       if (is_irs_k) then
          if (iirs==2) apply_irs_3_rans => irs2_ngh_k_rans
          if (iirs==4) apply_irs_3_rans => irs4_ngh_k_rans
       else
          apply_irs_3_rans => voidp
       endif

       ! Main IRS routine
       ! ================
       irs_rans => irs_routine_rans
  
       ! communication routine for increments
       ! ====================================
       if (is_2D) then
          communication_inc_rans => communication_inc2d_rans
       else
          communication_inc_rans => communication_inc3d_rans
       endif
    endif

  end subroutine assign_irs

  !===============================================================================
  subroutine irs_routine
  !===============================================================================
    !> Application of the Implicit Residual Smoothing (IRS) ** TO BE CHANGED
  !===============================================================================
    use mod_irs_d ! <- for irs_ngh_d
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    call update_var_in_dw

    call irs_ngh_d

    call update_var_of_dw

  end subroutine irs_routine

  !===============================================================================
  subroutine irs_routine_rans
  !===============================================================================
    !> Application of the Implicit Residual Smoothing (IRS) ** TO BE CHANGED
  !===============================================================================
    use mod_irs_d_rans ! <- for irs_ngh_d_rans
    use mod_interface  ! <- for update_var_rans
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    call irs_ngh_d_rans

    call update_var_rans

  end subroutine irs_routine_rans

end module mod_init_irs
