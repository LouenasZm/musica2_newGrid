!==================================================================================
module mod_mpi_types_one_sided
!==================================================================================
  !> Module to create MPI types
!==================================================================================
  use mod_mpi_part
  use mod_grid_directions
  implicit none
  ! -------------------------------------------------------------------------------
  ! MPI types for conservative variables (with ngh ghost points)
  ! ------------------------------------
  ! types for sending from origin
  integer :: type_facei,type_facej,type_facek
  ! types & displacements in target
  integer :: type_face_target1,type_face_target2
  integer :: type_face_target3,type_face_target4
  integer(kind=MPI_ADDRESS_KIND) :: disp_target(6)
  ! -------------------------------------------------------------------------------
  ! MPI types for increments (with ngh_irs(1:6) ghost points)
  ! ------------------------
  ! types for sending from origin
  !integer :: type_inci,type_incj,type_inck
  integer :: type_inc1,type_inc2,type_inc3
  integer :: type_inc4,type_inc5,type_inc6
  ! types & displacements in target
  integer :: type_inc_target1,type_inc_target2
  integer :: type_inc_target3,type_inc_target4
  integer(kind=MPI_ADDRESS_KIND) :: disp_inc_target(6)
  ! -------------------------------------------------------------------------------
contains

  !================================================================================
  subroutine mpi_types_comm1
  !================================================================================
    !> Create an MPI types for one-sided communications
    !> (curvilinear version with swap and reverse directions)
    ! Nota:
    ! -----
    ! face imin : numbered 1
    ! face imax : numbered 2
    ! face jmin : numbered 3
    ! face jmax : numbered 4
    ! face kmin : numbered 5
    ! face kmax : numbered 6
  !================================================================================
    implicit none
    ! -----------------------------------------------------------------------------  
    integer :: i
    integer :: nex,ney ! extended sizes
    integer :: sign_i,sign_j
    integer :: memjump ! jump in memory between block starts
    ! strides to construct MPI types
    integer(kind=MPI_ADDRESS_KIND) :: stride,stride2
    ! intermediate MPI types (not commited)
    integer :: type_base,type_edge
    ! target displacements
    integer, dimension(6) :: displ
    ! -----------------------------------------------------------------------------  

    ! Regular MPI types for one-sided comm
    ! ====================================
    ! used for sending data in directions i,j
    ! used for sending & receiving in direction k

    ! Extended sizes (+ ghost cells)
    ! ------------------------------
    nex=nx+2*ngh
    ney=ny+2*ngh

    ! stride in i- or j-directions
    ! ----------------------------
    stride =nex*sizeofreal
    stride2=nex*ney*sizeofreal
    
    ! MPI-type construction for imin/imax faces (normal to i-direction)
    ! -----------------------------------------
    call MPI_TYPE_VECTOR(ngh,1,1,MPI_DOUBLE_PRECISION,type_base,info)
    call MPI_TYPE_CREATE_HVECTOR(ny,1,stride, type_base,type_edge,info)
    call MPI_TYPE_CREATE_HVECTOR(nz,1,stride2,type_edge,type_facei,info)
    call MPI_TYPE_COMMIT(type_facei,info)

    ! MPI-type construction for jmin/jmax faces (normal to j-direction)
    ! -----------------------------------------
    call MPI_TYPE_VECTOR(nx,1,1,MPI_DOUBLE_PRECISION,type_base,info)
    call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride, type_base,type_edge,info)
    call MPI_TYPE_CREATE_HVECTOR( nz,1,stride2,type_edge,type_facej,info)
    call MPI_TYPE_COMMIT(type_facej,info)

    ! MPI-type construction for kmin/kmax faces (normal to k-direction)
    ! -----------------------------------
    call MPI_TYPE_VECTOR(nx,1,1,MPI_DOUBLE_PRECISION,type_base,info)
    call MPI_TYPE_CREATE_HVECTOR( ny,1,stride, type_base,type_edge,info)
    call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride2,type_edge,type_facek,info)
    call MPI_TYPE_COMMIT(type_facek,info)
    
    ! Define sending indices
    ! ======================
    in1s=1
    in2s=nx-ngh+1
    in3s=1
    in4s=ny-ngh+1
    if (is_adjoint_block) then ! for adjoint block interfaces
       if (neighbor(1).ne.MPI_PROC_NULL) then
          if (nob(neighbor(1))/=nob(iproc)) in1s=in1s+1
       endif
       if (neighbor(2).ne.MPI_PROC_NULL) then
          if (nob(neighbor(2))/=nob(iproc)) in2s=in2s-1
       endif
       if (neighbor(3).ne.MPI_PROC_NULL) then
          if (nob(neighbor(3))/=nob(iproc)) in3s=in3s+1
       endif
       if (neighbor(4).ne.MPI_PROC_NULL) then
          if (nob(neighbor(4))/=nob(iproc)) in4s=in4s-1
       endif
    endif

    ! Compute target displacements
    ! ============================
    call mpi_disp_comm1(displ)
    disp_target=displ
    
    ! Define target MPI types
    ! =======================
    ! Hint: since we assume for the moment that number of points in one direction is proportional
    !       to the number of procs in that direction (nx,ny) is the same in all subdomains
    !       For instance, in case of swap, nx_neighbor=ny and ny_neighbor=nx
    !       This particular property is used to avoid getting informations by SEND/RECV about
    !       the dimensions of the neighboring subdomain, which is the target of PUT.
    ! Note the use of is_swap2 and is_rev2 (bilateral versions of the indicators determined
    ! in mod_grid_directions because origin MPI types are not modified.

    do i=1,4
       ! Target #i
       ! --------- 
       ! origin is imin (i=1), imax (i=2), jmin (i=3) or jmax (i=4) and
       ! receiver depends on swap and reverse on i,j directions
       !          ~> i=1: the regular target would be imax without swap & reverse
       !          ~> i=2: the regular target would be imin without swap & reverse
       !          ~> i=3: the regular target would be jmax without swap & reverse
       !          ~> i=4: the regular target would be jmin without swap & reverse
       sign_i=1
       sign_j=1
       ! for i=1,2 faces along i: origin parallel dir. is j and normal dir. is i
       ! -----------------------------------------------------------------------
       if (i==1.or.i==2) then
          ! if rev parallel -> +j => -j
          if (is_rev2(i,1)) sign_j=-1
          ! if rev  normal  -> +i => -i
          if (is_rev2(i,2)) sign_i=-1
          ! for i=3,4 faces along j: origin parallel dir. is i and normal dir. is j
          ! -----------------------------------------------------------------------
       elseif (i==3.or.i==4) then
          ! if rev parallel -> +i => -i
          if (is_rev2(i,1)) sign_i=-1
          ! if rev  normal  -> +j => -j
          if (is_rev2(i,2)) sign_j=-1
       endif

       ! Define memory jumps
       ! -------------------
       if (is_swapij2(i)) then
          ! if swap invert role of i and j (nx_neighbor=ny & ny_neighbor=nx)
          ! memjump=sign_i*nex ! (nex_neighbor=ney)
          memjump=sign_i*ney   ! (sign_j_neighbor=sign_i)     
          stride =sign_j*sizeofreal ! i (sign_i_neighbor=sign_j)
          ! extended sizes (+ ghost cells) relative to receiving domain (target=neighbor)
          nex=ny+2*ngh ! extended i-dim (with ghost cells). Note that nx_neighbor=ny
          ney=nx+2*ngh ! extended j-dim (with ghost cells). Note that ny_neighbor=nx
       else
          ! default
          memjump=sign_i ! i
          stride =sign_j*nex*sizeofreal ! j
          ! extended sizes (+ ghost cells) relative to receiving domain (target=neighbor)
          nex=nx+2*ngh ! extended i-dim (with ghost cells).
          ney=ny+2*ngh ! extended j-dim (with ghost cells).
       endif

       ! stride between edges in i- or j-directions
       ! ------------------------------------------
       stride2=nex*ney*sizeofreal

       ! MPI-type construction for target
       ! --------------------------------
       select case (i)
       case (1)
          call MPI_TYPE_VECTOR(ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
          call MPI_TYPE_CREATE_HVECTOR(ny,1,stride,type_base,type_edge,info)
          ! extend in 3D and commit
          call MPI_TYPE_CREATE_HVECTOR(nz,1,stride2,type_edge,type_face_target1,info)
          call MPI_TYPE_COMMIT(type_face_target1,info)
       case (2)
          call MPI_TYPE_VECTOR(ngh,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
          call MPI_TYPE_CREATE_HVECTOR(ny,1,stride,type_base,type_edge,info)
          ! extend in 3D and commit
          call MPI_TYPE_CREATE_HVECTOR(nz,1,stride2,type_edge,type_face_target2,info)
          call MPI_TYPE_COMMIT(type_face_target2,info)
       case (3)
          call MPI_TYPE_VECTOR(nx,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
          call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride,type_base,type_edge,info)
          ! extend in 3D and commit
          call MPI_TYPE_CREATE_HVECTOR(nz,1,stride2,type_edge,type_face_target3,info)
          call MPI_TYPE_COMMIT(type_face_target3,info)
       case (4)
          call MPI_TYPE_VECTOR(nx,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
          call MPI_TYPE_CREATE_HVECTOR(ngh,1,stride,type_base,type_edge,info)
          ! extend in 3D and commit
          call MPI_TYPE_CREATE_HVECTOR(nz,1,stride2,type_edge,type_face_target4,info)
          call MPI_TYPE_COMMIT(type_face_target4,info)
       end select
    enddo

  end subroutine mpi_types_comm1

  !================================================================================
  subroutine mpi_disp_comm1(displ)
  !================================================================================
    !> Compute MPI displacements in target for one-sided communications
    !> for a given number of ghost points ngh
    !> (curvilinear version with swap and reverse directions)
    !
    ! Nota:
    ! -----
    ! face imin : numbered 1 -> displ(1)
    ! face imax : numbered 2 -> displ(2)
    ! face jmin : numbered 3 -> displ(3)
    ! face jmax : numbered 4 -> displ(4)
    ! face kmin : numbered 5 -> displ(5)
    ! face kmax : numbered 6 -> displ(6) 
  !================================================================================
    implicit none
    ! -----------------------------------------------------------------------------
    integer, dimension(6), intent(out) :: displ
    ! -----------------------------------------------------------------------------  
    integer :: nex,ney ! extended sizes
    integer, dimension(6) :: dispb ! base displacements
    ! -----------------------------------------------------------------------------  

    ! Base displacements (i.e. without swap & reverse)
    ! ==================
    ! give locations of the initial point where data are put
    
    ! Extended sizes (+ ghost cells)
    ! ------------------------------
    nex=nx+2*ngh
    ney=ny+2*ngh

    ! reception in imin target (origin from imax without swap,rev)
    ! ------------------------
    dispb(1)=nex*ngh ! skip lower ghost points (ngh extended lines)
    
    ! reception in imax target (origin from imax without swap,rev)
    ! ------------------------
    dispb(2)=nex*ngh+ngh+nx ! skip lower ghost points & traverse ngh+nx
    
    ! reception in jmin target (origin from jmax without swap,rev)
    ! ------------------------
    dispb(3)=ngh ! skip ngh ghost points (left)
    
    ! reception in jmax target (origin from jmin without swap,rev)
    ! ------------------------
    dispb(4)=nex*(ngh+ny)+ngh ! traverse all domain (ngh at bottom + ny full lines) + ngh
    
    ! reception in kmin target (origin from kmax)
    ! ------------------------
    dispb(5)=nex*ngh+ngh ! skip all bottom layer of ghost points + position after ngh points
    
    ! reception in kmax target (origin from kmin)
    ! ------------------------
    dispb(6)=nex*ngh+nex*ney*(nz+ngh)+ngh ! traverse volume (nex*ney*(nz+ngh)) + layer of ghost points + ngh
    
    ! Define target displacements
    ! ===========================
    ! Hint: since we assume for the moment that number of points in one direction is proportional
    !       to the number of procs in that direction (nx,ny) is the same in all subdomains
    !       For instance, in case of swap, nx_neighbor=ny and ny_neighbor=nx
    !       This particular property is used to avoid getting informations by SEND/RECV about
    !       the dimensions of the neighboring subdomain, which is the target of PUT.
    ! Note the use of is_swap2 and is_rev2 (bilateral versions of the indicators determined
    ! in mod_grid_directions) because origin MPI types are not modified.

    ! Target #1: origin is imin and receiver depends on swap and reverse on i,j directions
    ! ---------
   
    ! Determine displacement in target depending on its direction
    ! the direction is given by ndir(i), filled in mod_grid_directions    
    if (is_swapij2(1)) then
       if (ndir(1)==3) displ(1)=dispb(3) ! swap (neighbor dim not used)
       if (ndir(1)==4) displ(1)=ngh+ney*(nx+ngh) ! swap (neighbor dim: nx_neighbor=ny & ny_neighbor=nx)
                                                                          ! skip nex_v*(ny_v+nghbottom)
       ! in case reverse parallel (i), receive target from right
       if (is_rev2(1,1)) displ(1)=displ(1)+ny-1
       ! in case reverse normal (j), receive target from top
       if (is_rev2(1,2)) displ(1)=displ(1)+(ngh-1)*ney ! add ngh-1 full lines
    else
       if (ndir(1)==1) displ(1)=dispb(1) ! reverse i
       if (ndir(1)==2) displ(1)=dispb(2) ! regular value
       ! in case reverse parallel (j), receive target from top
       if (is_rev2(1,1)) displ(1)=displ(1)+nex*(ny-1)
       ! in case reverse normal (i), receive target from right
       if (is_rev2(1,2)) displ(1)=displ(1)+ngh-1
    endif
 
    ! Target #2: origin is imax and receiver depends on swap and reverse on i,j directions
    ! ---------

    ! Determine displacement in target depending on its direction
    ! the direction is given by ndir(i), filled in mod_grid_directions
    if (is_swapij2(2)) then
       if (ndir(2)==3) displ(2)=dispb(3) ! swap (neighbor dim not used)
       if (ndir(2)==4) displ(2)=ngh+ney*(nx+ngh) ! swap (neighbor dim: nx_neighbor=ny & ny_neighbor=nx)
       ! in case reverse parallel (i), receive target from right
       if (is_rev2(2,1)) displ(2)=displ(2)+ny-1
       ! in case reverse normal (j), receive target from top
       if (is_rev2(2,2)) displ(2)=displ(2)+(ngh-1)*ney ! add ngh-1 full lines
    else
       if (ndir(2)==1) displ(2)=dispb(1) ! regular value
       if (ndir(2)==2) displ(2)=dispb(2) ! reverse i
       ! in case reverse parallel (j), receive target from top
       if (is_rev2(2,1)) displ(2)=displ(2)+nex*(ny-1)
       ! in case reverse normal (i), receive target from right
       if (is_rev2(2,2)) displ(2)=displ(2)+ngh-1
    endif
 
    ! Target #3: origin is jmin and receiver depends on swap and reverse on i,j directions
    ! ---------
    
    ! Determine displacement in target depending on its direction
    ! the direction is given by ndir(i), filled in mod_grid_directions    
    if (is_swapij2(3)) then
       if (ndir(3)==1) displ(3)=ney*ngh ! swap (neighbor dim: nx_neighbor=ny)
       if (ndir(3)==2) displ(3)=ney*ngh+ngh+ny ! swap (neighbor dim: nx_neighbor=ny)
       !  in case reverse parallel (i), receive target from right
       if (is_rev2(3,1)) displ(3)=displ(3)+ney*(nx-1)
       !  in case reverse normal (j), receive target from right
       if (is_rev2(3,2)) displ(3)=displ(3)+ngh-1
    else    
       if (ndir(3)==3) displ(3)=dispb(3) ! inverse
       if (ndir(3)==4) displ(3)=dispb(4) ! regular value
       ! in case reverse normal (j), receive target from top
       if (is_rev2(3,2)) displ(3)=displ(3)+nex*(ngh-1)
       ! in case reverse parallel (i), receive target from right
       if (is_rev2(3,1)) displ(3)=displ(3)+nx-1
    endif

    ! Target #4: origin is jmax and receiver depends on swap and reverse on i,j directions
    ! ---------
    
    ! Determine displacement in target depending on its direction
    ! the direction is given by ndir(i), filled in mod_grid_directions
    if (is_swapij2(4)) then
       if (ndir(4)==1) displ(4)=ney*ngh ! swap (neighbor dim: nx_neighbor=ny)
       if (ndir(4)==2) displ(4)=ney*ngh+ngh+ny ! swap (neighbor dim: nx_neighbor=ny)
       !  in case reverse parallel (i), receive target from right
       if (is_rev2(4,1)) displ(4)=displ(4)+ney*(nx-1)
       !  in case reverse normal (j), receive target from right
       if (is_rev2(4,2)) displ(4)=displ(4)+ngh-1
    else
       if (ndir(4)==3) displ(4)=dispb(3) ! regular value
       if (ndir(4)==4) displ(4)=dispb(4) ! inverse
       ! in case reverse normal (j), receive target from top
       if (is_rev2(4,2)) displ(4)=displ(4)+nex*(ngh-1)
       !if (is_rev2(4,2)) displ(4)=nex*(ney-1)+ngh
       ! in case reverse parallel (i), receive target from right
       if (is_rev2(4,1)) displ(4)=displ(4)+nx-1
    endif

    ! in 3D, skip ghost layers in front face (k=-4,0)
    ! -----------------------------------------------
    if (.not.is_2D) displ(1:4)=displ(1:4)+ngh*nex*ney

    ! Target #5: origin is kmin and receiver (target) is kmax
    ! ---------
    ! displacement
    displ(5)=dispb(6)
    
    ! Target #6: origin is kmax and receiver (target) is kmin
    ! ---------
    ! displacement
    displ(6)=dispb(5)
 
  end subroutine mpi_disp_comm1

  !================================================================================
  subroutine mpi_types_comm1_inc
  !================================================================================
    !> Create MPI types for one-sided communications of increments
    !> (use ngh_irs(1:6) ghost points for interfaces)
    !> (curvilinear version with swap and reverse directions)
    !
    ! Nota:
    ! -----
    ! face imin : numbered 1
    ! face imax : numbered 2
    ! face jmin : numbered 3
    ! face jmax : numbered 4
    ! face kmin : numbered 5
    ! face kmax : numbered 6
  !================================================================================
    use mod_bc
    implicit none
    ! -----------------------------------------------------------------------------
    integer :: i,j,k,dim
    integer :: nex,ney ! extended sizes
    ! strides to construct MPI types
    integer(kind=MPI_ADDRESS_KIND) :: stride,stride2
    ! intermediate MPI types (not commited)
    integer :: type_base,type_edge
    ! -----------------------------------------------------------------------------  

    ! Extended sizes (+ ghost cells)
    ! ------------------------------
    nex=nx+ngh_irs(1)+ngh_irs(2)
    ney=ny+ngh_irs(3)+ngh_irs(4)

    ! Regular MPI types for one-sided comm (relative to sending MPI domain)
    ! ====================================
    ! used for sending data in directions i,j
    ! used for sending & receiving in direction k
    
    ! stride in i- or j-directions
    ! ----------------------------
    stride =nex*sizeofreal
    stride2=nex*ney*sizeofreal
    
    ! MPI-type construction for imin faces (normal to i-direction)
    ! ------------------------------------
    call MPI_TYPE_VECTOR(ngh_irs(1),1,1,MPI_DOUBLE_PRECISION,type_base,info)
    call MPI_TYPE_CREATE_HVECTOR(ny,1,stride,type_base,type_edge,info)
    call MPI_TYPE_CREATE_HVECTOR(nz,1,stride2,type_edge,type_inc1,info)
    call MPI_TYPE_COMMIT(type_inc1,info)

    ! MPI-type construction for imax faces (normal to i-direction)
    ! ------------------------------------
    call MPI_TYPE_VECTOR(ngh_irs(2),1,1,MPI_DOUBLE_PRECISION,type_base,info)
    call MPI_TYPE_CREATE_HVECTOR(ny,1,stride,type_base,type_edge,info)
    call MPI_TYPE_CREATE_HVECTOR(nz,1,stride2,type_edge,type_inc2,info)
    call MPI_TYPE_COMMIT(type_inc2,info)

    ! MPI-type construction for jmin faces (normal to j-direction)
    ! ------------------------------------
    call MPI_TYPE_VECTOR(nx,1,1,MPI_DOUBLE_PRECISION,type_base,info)
    call MPI_TYPE_CREATE_HVECTOR(ngh_irs(3),1,stride,type_base,type_edge,info)
    call MPI_TYPE_CREATE_HVECTOR(nz,1,stride2,type_edge,type_inc3,info)
    call MPI_TYPE_COMMIT(type_inc3,info)

    ! MPI-type construction for jmax faces (normal to j-direction)
    ! ------------------------------------
    call MPI_TYPE_VECTOR(nx,1,1,MPI_DOUBLE_PRECISION,type_base,info)
    call MPI_TYPE_CREATE_HVECTOR(ngh_irs(4),1,stride,type_base,type_edge,info)
    call MPI_TYPE_CREATE_HVECTOR(nz,1,stride2,type_edge,type_inc4,info)
    call MPI_TYPE_COMMIT(type_inc4,info)

    ! MPI-type construction for kmin faces (normal to k-direction)
    ! ------------------------------------
    call MPI_TYPE_VECTOR(nx,1,1,MPI_DOUBLE_PRECISION,type_base,info)
    call MPI_TYPE_CREATE_HVECTOR(ny,1,stride,type_base,type_edge,info)
    call MPI_TYPE_CREATE_HVECTOR(ngh_irs(5),1,stride2,type_edge,type_inc5,info)
    call MPI_TYPE_COMMIT(type_inc5,info)
    
    ! MPI-type construction for kmax faces (normal to k-direction)
    ! ------------------------------------
    call MPI_TYPE_VECTOR(nx,1,1,MPI_DOUBLE_PRECISION,type_base,info)
    call MPI_TYPE_CREATE_HVECTOR(ny,1,stride,type_base,type_edge,info)
    call MPI_TYPE_CREATE_HVECTOR(ngh_irs(6),1,stride2,type_edge,type_inc6,info)
    call MPI_TYPE_COMMIT(type_inc6,info)

    ! Define sending indices
    ! ======================
    ii1s=1
    ii2s=nx-ngh_irs(2)+1
    ii3s=1
    ii4s=ny-ngh_irs(4)+1
    ii5s=1
    ii6s=nz-ngh_irs(6)+1
    if (is_adjoint_block) then ! for adjoint block interfaces
       if (neighbor(1).ne.MPI_PROC_NULL) then
          if (nob(neighbor(1))/=nob(iproc)) ii1s=ii1s+1
       endif
       if (neighbor(2).ne.MPI_PROC_NULL) then
          if (nob(neighbor(2))/=nob(iproc)) ii2s=ii2s-1
       endif
       if (neighbor(3).ne.MPI_PROC_NULL) then
          if (nob(neighbor(3))/=nob(iproc)) ii3s=ii3s+1
       endif
       if (neighbor(4).ne.MPI_PROC_NULL) then
          if (nob(neighbor(4))/=nob(iproc)) ii4s=ii4s-1
       endif
    endif

    ! Compute target displacement & create MPI type for comm1
    ! =======================================================
    if (is_2D) then
       dim=2
    else
       dim=3
    endif
    
    do i=1,dim
       do j=1,2
          k=2*(i-1)+j
          call create_target(k,BC_face(i,j)%ngh_irs)
       enddo
    enddo
    
  end subroutine mpi_types_comm1_inc
  
  !================================================================================
  subroutine create_target(i,ngh_)
  !================================================================================
    !> Compute MPI displacement & type for targets in one-sided communications
    !> for a given number of ghost points ngh_ in neighboring block
    !> (curvilinear version with swap and reverse directions)
    !
    ! Nota:
    ! -----
    ! face imin : numbered 1 -> target1
    ! face imax : numbered 2 -> target2
    ! face jmin : numbered 3 -> target3
    ! face jmax : numbered 4 -> target4
    ! face kmin : numbered 5 -> target5
    ! face kmax : numbered 6 -> target6
  !================================================================================
    implicit none
    ! -----------------------------------------------------------------------------
    integer, intent(in) :: i ! target number
    integer, dimension(6), intent(in) :: ngh_ ! number of ghost cells in neighbor
    ! -----------------------------------------------------------------------------  
    integer :: nex,ney ! extended sizes
    integer :: displ
    integer :: sign_i,sign_j
    integer :: memjump ! jump in memory between block starts
    ! strides to construct MPI types
    integer(kind=MPI_ADDRESS_KIND) :: stride,stride2
    ! intermediate MPI types (not commited)
    integer :: type_base,type_edge
    ! -----------------------------------------------------------------------------     
    
    ! Hint: since we assume for the moment that number of points in one direction is proportional
    !       to the number of procs in that direction (nx,ny) is the same in all subdomains
    !       For instance, in case of swap, nx_neighbor=ny and ny_neighbor=nx
    !       This particular property is used to avoid getting informations by SEND/RECV about
    !       the dimensions of the neighboring subdomain, which is the target of PUT.
    ! Note the use of is_swap2 and is_rev2 (bilateral versions of the indicators determined
    ! in mod_grid_directions because origin MPI types are not modified.
    
    ! Determine displacement in target depending on its direction
    ! ===========================================================
    ! the direction is given by ndir(i), filled in mod_grid_directions
    
    ! for i=1,2 faces along i: origin parallel dir. is j and normal dir. is i
    ! -----------------------------------------------------------------------
    if (i==1.or.i==2) then
       if (is_swapij2(i)) then
          nex=ny+ngh_(1)+ngh_(2) ! extended i-dim (with ghost cells). Note that nx_neighbor=ny
          ney=nx+ngh_(3)+ngh_(4) ! extended j-dim (with ghost cells). Note that ny_neighbor=nx

          if (ndir(i)==3) displ=ngh_(1) ! swap (neighbor dim not used)
          if (ndir(i)==4) displ=ngh_(1)+nex*(nx+ngh_(3)) ! swap (neighbor dim: nx_neighbor=ny & ny_neighbor=nx)
          ! skip nex_v*(ny_v+ngh_bottom)
          ! in case reverse parallel, receive target from right
          if (is_rev2(i,1)) displ=displ+ny-1 ! (nx_neighbor=ny)
          ! in case reverse normal, receive target from top
          if (is_rev2(i,2).and.i==1) displ=displ+(ngh_(3)-1)*nex ! add ngh_-1 full lines
          if (is_rev2(i,2).and.i==2) displ=displ+(ngh_(4)-1)*nex ! add ngh_-1 full lines
       else
          nex=nx+ngh_(1)+ngh_(2) ! extended i-dim (with ghost cells)
          ney=ny+ngh_(3)+ngh_(4) ! extended j-dim (with ghost cells)

          if (ndir(i)==1) displ=nex*ngh_(3) ! skip lower ghost points (ngh_ extended lines)
          if (ndir(i)==2) displ=nex*ngh_(3)+ngh_(1)+nx ! skip lower ghost points & traverse ngh_+nx
          ! in case reverse parallel, receive target from top
          if (is_rev2(i,1)) displ=displ+nex*(ny-1)
          ! in case reverse normal, receive target from right
          if (is_rev2(i,2).and.i==1) displ=displ+ngh_(1)-1
          if (is_rev2(i,2).and.i==2) displ=displ+ngh_(2)-1
       endif

    ! for i=3,4 faces along j: origin parallel dir. is i and normal dir. is j
    ! -----------------------------------------------------------------------
    elseif (i==3.or.i==4) then
       if (is_swapij2(i)) then
          nex=ny+ngh_(1)+ngh_(2) ! extended i-dim (with ghost cells). Note that nx_neighbor=ny
          ney=nx+ngh_(3)+ngh_(4) ! extended j-dim (with ghost cells). Note that ny_neighbor=nx

          if (ndir(i)==1) displ=nex*ngh_(3) ! swap (neighbor dim: nx_neighbor=ny)
          if (ndir(i)==2) displ=nex*ngh_(3)+ngh_(1)+ny ! swap (neighbor dim: nx_neighbor=ny)
          !  in case reverse parallel (i), receive target from right
          if (is_rev2(i,1)) displ=displ+nex*(nx-1)
          !  in case reverse normal (j), receive target from right
          if (is_rev2(i,2).and.i==3) displ=displ+ngh_(1)-1
          if (is_rev2(i,2).and.i==4) displ=displ+ngh_(2)-1
       else    
          nex=nx+ngh_(1)+ngh_(2) ! extended i-dim (with ghost cells)
          ney=ny+ngh_(3)+ngh_(4) ! extended j-dim (with ghost cells)

          if (ndir(i)==3) displ=ngh_(1)
          if (ndir(i)==4) displ=nex*(ngh_(3)+ny)+ngh_(1)
          ! in case reverse parallel (i), receive target from right
          if (is_rev2(i,1)) displ=displ+nx-1
          ! in case reverse normal (j), receive target from top
          if (is_rev2(i,2).and.i==3) displ=displ+nex*(ngh_(3)-1)
          if (is_rev2(i,2).and.i==4) displ=displ+nex*(ngh_(4)-1)
       endif
    endif
    
    ! in 3D, skip ghost layers in front face (k=-4,0)
    ! -----------------------------------------------    
    if (.not.is_2D) displ=displ+ngh_(5)*nex*ney

    ! for i=5,6 faces along k
    ! -----------------------
    if (i==5) then
       ! Target #5: origin is kmin and receiver (target) is kmax
       nex=nx+ngh_(1)+ngh_(2) ! extended i-dim (with ghost cells)
       ney=ny+ngh_(3)+ngh_(4) ! extended j-dim (with ghost cells)
       ! traverse volume (nex*ney*(nz+ngh_(5))) + layer of ghost points + ngh_(1)
       displ=nex*ney*(nz+ngh_(5))+nex*ngh_(3)+ngh_(1) 
    elseif (i==6) then
       ! Target #6: origin is kmax and receiver (target) is kmin
       nex=nx+ngh_(1)+ngh_(2) ! extended i-dim (with ghost cells)
       ! skip all bottom layer of ghost points + position after ngh_(1) points
       displ=nex*ngh_(3)+ngh_(1)
    endif
    
    ! Displacement in target
    ! ----------------------
    select case (i)
    case (1)
       disp_inc_target(1)=displ
    case (2)
       disp_inc_target(2)=displ
    case (3)
       disp_inc_target(3)=displ
    case (4)
       disp_inc_target(4)=displ
    case (5)
       disp_inc_target(5)=displ
    case (6)
       disp_inc_target(6)=displ
    end select

    ! no particular receiving types for kmin/kmax
    ! -------------------------------------------
    if (i==5.or.i==6) return

    ! Define target MPI types
    ! =======================

    ! Target #i
    ! --------- 
    ! origin is imin (i=1), imax (i=2), jmin (i=3) or jmax (i=4) and
    ! receiver depends on swap and reverse on i,j directions
    !          ~> i=1: the regular target would be imax without swap & reverse
    !          ~> i=2: the regular target would be imin without swap & reverse
    !          ~> i=3: the regular target would be jmax without swap & reverse
    !          ~> i=4: the regular target would be jmin without swap & reverse    
    sign_i=1
    sign_j=1
    ! for i=1,2 faces along i: origin parallel dir. is j and normal dir. is i
    ! -----------------------------------------------------------------------
    if (i==1.or.i==2) then
       ! if rev parallel -> +j => -j
       if (is_rev2(i,1)) sign_j=-1
       ! if rev  normal  -> +i => -i
       if (is_rev2(i,2)) sign_i=-1
    ! for i=3,4 faces along j: origin parallel dir. is i and normal dir. is j
    ! -----------------------------------------------------------------------
    elseif (i==3.or.i==4) then
       ! if rev parallel -> +i => -i
       if (is_rev2(i,1)) sign_i=-1
       ! if rev  normal  -> +j => -j
       if (is_rev2(i,2)) sign_j=-1
    endif
    
    ! Define memory jumps
    ! -------------------
    if (is_swapij2(i)) then
       ! if swap invert role of i and j (nx_neighbor=ny & ny_neighbor=nx)
       ! memjump=sign_i*nex  ! j (sign_j_neighbor=sign_i)
       !         with nex=nx_neighbor+ngh_(1)+ngh_(2)=ny+ngh_(1)+ngh_(2)
       memjump=sign_i*(ny+ngh_(1)+ngh_(2))  ! j (sign_j_neighbor=sign_i)     
       stride =sign_j*sizeofreal ! i (sign_i_neighbor=sign_j)
       ! extended sizes (+ ghost cells) relative to receiving domain (target=neighbor)
       nex=ny+ngh_(1)+ngh_(2) ! extended i-dim (with ghost cells). Note that nx_neighbor=ny
       ney=nx+ngh_(3)+ngh_(4) ! extended j-dim (with ghost cells). Note that ny_neighbor=nx
    else
       ! default
       memjump=sign_i ! i
       stride =sign_j*nex*sizeofreal ! j
       ! extended sizes (+ ghost cells) relative to receiving domain (target=neighbor)
       nex=nx+ngh_(1)+ngh_(2) ! extended i-dim (with ghost cells).
       ney=ny+ngh_(3)+ngh_(4) ! extended j-dim (with ghost cells).
    endif
    
    ! stride between edges in i- or j-directions
    ! ------------------------------------------
    stride2=nex*ney*sizeofreal
    
    ! MPI-type construction for target
    ! --------------------------------
    select case (i)
    case (1)
       call MPI_TYPE_VECTOR(ngh_(ndir(i)),1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ny,1,stride,type_base,type_edge,info)
       ! extend in 3D and commit
       call MPI_TYPE_CREATE_HVECTOR(nz,1,stride2,type_edge,type_inc_target1,info)
       call MPI_TYPE_COMMIT(type_inc_target1,info)
    case (2)
       call MPI_TYPE_VECTOR(ngh_(ndir(i)),1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ny,1,stride,type_base,type_edge,info)
       ! extend in 3D and commit
       call MPI_TYPE_CREATE_HVECTOR(nz,1,stride2,type_edge,type_inc_target2,info)
       call MPI_TYPE_COMMIT(type_inc_target2,info)
    case (3)
       call MPI_TYPE_VECTOR(nx,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ngh_(ndir(i)),1,stride,type_base,type_edge,info)
       ! extend in 3D and commit
       call MPI_TYPE_CREATE_HVECTOR(nz,1,stride2,type_edge,type_inc_target3,info)
       call MPI_TYPE_COMMIT(type_inc_target3,info)
    case (4)
       call MPI_TYPE_VECTOR(nx,1,memjump,MPI_DOUBLE_PRECISION,type_base,info)
       call MPI_TYPE_CREATE_HVECTOR(ngh_(ndir(i)),1,stride,type_base,type_edge,info)
       ! extend in 3D and commit
       call MPI_TYPE_CREATE_HVECTOR(nz,1,stride2,type_edge,type_inc_target4,info)
       call MPI_TYPE_COMMIT(type_inc_target4,info)
    end select

  end subroutine create_target
  
end module mod_mpi_types_one_sided
