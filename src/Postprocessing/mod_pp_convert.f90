!==============================================================================
module mod_pp_convert
!==============================================================================
  !> Module to convert format: binary <-> tecplot
  !> Nota: bin/tec is set by input/output formats in param.ini
!==============================================================================
  implicit none
  ! ---------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------

contains

  !============================================================================
  subroutine pp_convert
  !============================================================================
    !> main routine
  !============================================================================
    use mod_io_snapshots
    implicit none
    ! -------------------------------------------------------------------------

    ! Read plane between samples n1 & n2
    ! ==================================
    allocate(QQ(nx,ny,nz))
    allocate(ww(nx1:nx2,ny1:ny2,nz1:nz2))
    
    iblc_pp=nob(iproc)

    print *,iproc,nsr,iblc_pp

    ! MPI-IO read
    ! -----------
    call read_snapshot(nsr,iblc_pp,dirDATA)    

    call write_grid_snapshot(nsr)    
    !call write_snapshot(nsr)    

    !enddo
    
  end subroutine pp_convert

end module mod_pp_convert
