!===============================================================================
module mod_mpi
!===============================================================================
  !> Module for main MPI variables
!===============================================================================
  use mpi
  implicit none
  ! ----------------------------------------------------------------------------
  !> 1. Global communicator
  integer :: COMM_global ! name of communicator
  integer :: nproc       ! number of procs
  integer :: iproc       ! proc ID in global comm
  integer :: info        ! MPI variable
  integer :: icheck      ! MPI variable used in solver.f90
  integer :: ndomx,ndomy,ndomz ! number of MPI dom in each direction of a block
  integer, dimension(3) :: coord ! MPI cart coordinates within a block
  ! ----------------------------------------------------------------------------
  !> 2. MPI constants for Blocks
  integer, dimension(:), allocatable :: nob          ! block ID for a proc
  integer, dimension(:), allocatable :: iproc_leader ! proc ID of leader in one block
  ! ----------------------------------------------------------------------------

contains

  !===============================================================
  subroutine init_mpi
  !===============================================================
    !> Initialization of MPI
  !===============================================================
    implicit none
    ! ------------------------------------------------------------
    ! ------------------------------------------------------------

    ! Initialize MPI
    ! ==============
    call MPI_INIT(info)
    
    ! Rank of processes in the generic communicator MPI_COMM_WORLD
    ! ============================================================
    call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,info)
    
    ! Total number of processes (from mpirun -n nproc ...)
    ! =========================
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,info)

    ! Rename global communicator
    ! ==========================
    call MPI_COMM_DUP(MPI_COMM_WORLD,COMM_global,info)    

  end subroutine init_mpi

end module mod_mpi
