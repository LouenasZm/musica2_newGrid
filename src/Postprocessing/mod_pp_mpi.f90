!==============================================================================
module mod_pp_mpi
!==============================================================================
  !> Module for specific MPI parallelization of post-processing
!==============================================================================
  use mod_mpi_part
  implicit none
  ! ---------------------------------------------------------------------------
  ! MPI types for plane reconstruction
  integer :: type_p,type_pg
  integer :: type_in,typeg_in ! for global moments
  logical, dimension(:), allocatable :: is_in
  integer, dimension(:), allocatable :: coordin
  ! ---------------------------------------------------------------------------

contains

  !============================================================================
  subroutine mpi_types_sp
  !============================================================================
    !> Definition of MPI types for communications (frequency-wavenumber spectra)
  !============================================================================
    use warnstop
    use mod_grid
    use mod_io_snapshots
    use mod_pp_var
    implicit none
    ! -------------------------------------------------------------------------
    integer :: ip,type_base
    integer(kind=MPI_ADDRESS_KIND) :: pas
    ! -------------------------------------------------------------------------

    if (iproc==0) print *,'MPI addendum ...'
    
    ! MPI types for spectrum MPI reconstruction
    ! =========================================

    ! MPI types for planes
    ! ====================
    if (snapshots(nsr)%type.eq.2) then
       select case(snapshots(nsr)%normal)

       case (1) ! plane YZ (with normal 1 along X)

          call MPI_TYPE_VECTOR(ny,1,1,MPI_DOUBLE_PRECISION,type_base,info)

          ! Construction of type "yz-portion on global mesh"
          pas=ngy*sizeofreal
          call MPI_TYPE_CREATE_HVECTOR(nz,1,pas,type_base,type_pg,info)
          call MPI_TYPE_COMMIT(type_pg,info)

          ! Construction of type "yz-portion on local mesh"
          pas=ny*sizeofreal
          call MPI_TYPE_CREATE_HVECTOR(nz,1,pas,type_base,type_p,info)
          call MPI_TYPE_COMMIT(type_p,info)

       case (2) ! plane XZ (with normal 2 along Y)

          call MPI_TYPE_VECTOR(nx,1,1,MPI_DOUBLE_PRECISION,type_base,info)

          ! Construction of type "xz-portion on global mesh"
          pas=ngx*sizeofreal
          call MPI_TYPE_CREATE_HVECTOR(nz,1,pas,type_base,type_pg,info)
          call MPI_TYPE_COMMIT(type_pg,info)

          ! Construction of type "xz-portion on local mesh"
          pas=nx*sizeofreal
          call MPI_TYPE_CREATE_HVECTOR(nz,1,pas,type_base,type_p,info)
          call MPI_TYPE_COMMIT(type_p,info)

       case (3) ! plane XY (with normal 3 along Z)

          call MPI_TYPE_VECTOR(nx,1,1,MPI_DOUBLE_PRECISION,type_base,info)

          ! Construction of type "xy-portion on global mesh"
          pas=ngx*sizeofreal
          call MPI_TYPE_CREATE_HVECTOR(ny,1,pas,type_base,type_pg,info)
          call MPI_TYPE_COMMIT(type_pg,info)

          ! Construction of type "xy-portion on local mesh"
          pas=nx*sizeofreal
          call MPI_TYPE_CREATE_HVECTOR(ny,1,pas,type_base,type_p,info)
          call MPI_TYPE_COMMIT(type_p,info)

       case default

        call mpistop('error in snapshot number selection in param_pp.ini', 0)

       end select

    ! MPI types for lines
    ! ===================
    ! <~ To be changed, directly adapted from planes by AB
    else if (snapshots(nsr)%type.eq.1) then

       select case(snapshots(nsr)%dir)

       case (1) ! line X (normal to plane YZ)

          call MPI_TYPE_VECTOR(nx,1,1,MPI_DOUBLE_PRECISION,type_base,info)

          ! Construction of type "x-portion on global mesh"
          pas=ngx*sizeofreal
          call MPI_TYPE_CREATE_HVECTOR(1,1,pas,type_base,type_pg,info)
          call MPI_TYPE_COMMIT(type_pg,info)

          ! Construction of type "x-portion on local mesh"
          pas=nx*sizeofreal
          call MPI_TYPE_CREATE_HVECTOR(1,1,pas,type_base,type_p,info)
          call MPI_TYPE_COMMIT(type_p,info)

       case (2) ! line Y (normal to plane XZ)

          call MPI_TYPE_VECTOR(ny,1,1,MPI_DOUBLE_PRECISION,type_base,info)

          ! Construction of type "y-portion on global mesh"
          pas=ngy*sizeofreal
          call MPI_TYPE_CREATE_HVECTOR(1,1,pas,type_base,type_pg,info)
          call MPI_TYPE_COMMIT(type_pg,info)

          ! Construction of type "y-portion on local mesh"
          pas=ny*sizeofreal
          call MPI_TYPE_CREATE_HVECTOR(1,1,pas,type_base,type_p,info)
          call MPI_TYPE_COMMIT(type_p,info)

       case (3) ! line Z (normal to plane XY)

          call MPI_TYPE_VECTOR(1,1,1,MPI_DOUBLE_PRECISION,type_base,info)

          ! Construction of type "yz-portion on global mesh"
          pas=1*sizeofreal
          call MPI_TYPE_CREATE_HVECTOR(nz,1,pas,type_base,type_pg,info)
          call MPI_TYPE_COMMIT(type_pg,info)

          ! Construction of type "yz-portion on local mesh"
          pas=1*sizeofreal
          call MPI_TYPE_CREATE_HVECTOR(nz,1,pas,type_base,type_p,info)
          call MPI_TYPE_COMMIT(type_p,info)

       case default

        call mpistop('error in snapshot number selection in param_pp.ini', 0)

       end select

    endif
    
    ! MPI types for MPI reconstruction in inhomogeneous direction (if present)
    ! ===========================================================
    if (i_in>0) then

       ! Determination of procs that communicate
       ! ---------------------------------------
       allocate(is_in(0:nproc-1),coordin(0:nproc-1))

       is_in=.false.
       select case(name_in) ! inhomogeneous direction
       case('x')
          if ((coord(2)==0).and.(coord(3)==0)) is_in(iproc)=.true.
          coordin=coordx
       case('y')
          if ((coord(1)==0).and.(coord(3)==0)) is_in(iproc)=.true.
          coordin=coordy
       case('z')
          if ((coord(1)==0).and.(coord(2)==0)) is_in(iproc)=.true.
          coordin=coordz
       case('t')
          call mpistop('inhomogeneous time direction NOT IMPLEMENTED',0)       
       end select
       do ip=0,nproc-1
          call MPI_BCAST(is_in(ip),1,MPI_INTEGER,ip,COMM_global,info)
       enddo

       ! Create MPI types for communications
       ! -----------------------------------
       call MPI_TYPE_VECTOR(n_in,1,1,MPI_DOUBLE_PRECISION,type_base,info)
       ! create type for local grid
       call MPI_TYPE_HVECTOR(1,1,n_in*sizeofreal,type_base,type_in,info)
       call MPI_TYPE_COMMIT(type_in,info)
       ! create type for global grid
       call MPI_TYPE_HVECTOR(1,1,ng_in*sizeofreal,type_base,typeg_in,info)
       call MPI_TYPE_COMMIT(typeg_in,info)
       
    endif
    
  end subroutine mpi_types_sp

end module mod_pp_mpi
