!=================================================================================
module mod_comm1
!=================================================================================
  !> Module for one-sided RMA MPI communications
!=================================================================================
  use mod_mpi_types_one_sided
  use mod_flow
  implicit none
  !-------------------------------------------------------------------------------
  ! Communications one-sided RMA MPI
  ! ============================================
  integer :: win_Krho,win_Krhou,win_Krhov,win_Krhow,win_Krhoe
  integer :: win_rho,win_rhou,win_rhov,win_rhow,win_rhoe
  integer :: win_cfl
  integer :: group_Krho,group_Krhou,group_Krhov,group_Krhow,group_Krhoe
  !-------------------------------------------------------------------------------

contains

  !===============================================================================
  !===============================================================================
  ! I - Communications of faces for conservative variables
  !===============================================================================
  !===============================================================================
  
  !===============================================================================
  subroutine mpi_win_comm1
  !===============================================================================
    !> Window & Displacements for one-sided RMA communications
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: nex,ney,nez ! extended sizes
    integer(KIND=MPI_ADDRESS_KIND) :: win_buffer_size
    ! ----------------------------------------------------------------------------

    ! Extended sizes (+ ghost cells)
    ! ------------------------------
    nex=nx+2*ngh
    ney=ny+2*ngh
    if (is_2D) then
       nez=1
    else
       nez=nz+2*ngh
    endif
       
    ! Create window for increments
    ! ----------------------------
    win_buffer_size=nex*ney*nez*sizeofreal
    call MPI_WIN_CREATE(rho_n ,win_buffer_size,sizeofreal,MPI_INFO_NULL,COMM_global,win_rho ,info)
    call MPI_WIN_CREATE(rhou_n,win_buffer_size,sizeofreal,MPI_INFO_NULL,COMM_global,win_rhou,info)
    call MPI_WIN_CREATE(rhov_n,win_buffer_size,sizeofreal,MPI_INFO_NULL,COMM_global,win_rhov,info)
    if (.not.is_2D) call MPI_WIN_CREATE(rhow_n,win_buffer_size,sizeofreal,MPI_INFO_NULL,COMM_global,win_rhow,info)
    call MPI_WIN_CREATE(rhoe_n,win_buffer_size,sizeofreal,MPI_INFO_NULL,COMM_global,win_rhoe,info)

!!$    call MPI_WIN_ALLOCATE(win_buffer_size,sizeofreal,MPI_INFO_NULL,COMM_global,win_rho ,rho_n, info)
!!$    call MPI_WIN_ALLOCATE(win_buffer_size,sizeofreal,MPI_INFO_NULL,COMM_global,win_rhou,rhou_n,info)
!!$    call MPI_WIN_ALLOCATE(win_buffer_size,sizeofreal,MPI_INFO_NULL,COMM_global,win_rhov,rhov_n,info)
!!$    if (is_2D) call MPI_WIN_ALLOCATE(win_buffer_size,sizeofreal,MPI_INFO_NULL,COMM_global,win_rhow,rhow_n,info)
!!$    call MPI_WIN_ALLOCATE(win_buffer_size,sizeofreal,MPI_INFO_NULL,COMM_global,win_rhoe,rhoe_n,info)

  end subroutine mpi_win_comm1
  
  !===============================================================================
  subroutine mpi_close_win_comm1
  !===============================================================================
    !> CLose open MPI memory windows used for one-sided communications
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------
    
    ! close open memory windows
    call MPI_WIN_FREE(win_rho,info)
    call MPI_WIN_FREE(win_rhou,info)
    call MPI_WIN_FREE(win_rhov,info)
    if (.not.is_2D) call MPI_WIN_FREE(win_rhow,info)
    call MPI_WIN_FREE(win_rhoe,info)

  end subroutine mpi_close_win_comm1
  
  !===============================================================================
  subroutine communication1
  !===============================================================================
    !> Call one-sided RMA communications for conservative variables (rho_n,...)
    !> [1/ set window; 2/ comm; 3/ check]
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------
    
!    ! start check
    call MPI_WIN_FENCE(0,win_rho, info)
    call MPI_WIN_FENCE(0,win_rhou,info)
    call MPI_WIN_FENCE(0,win_rhov,info)
    if (.not.is_2D) call MPI_WIN_FENCE(0,win_rhow,info)
    call MPI_WIN_FENCE(0,win_rhoe,info)

    ! start commnunications
    call commun1(rho_n, win_rho)
    call commun1(rhou_n,win_rhou)
    call commun1(rhov_n,win_rhov)
    if (.not.is_2D) call commun1(rhow_n,win_rhow)
    call commun1(rhoe_n,win_rhoe)
    
    ! check end of commnunication
    call MPI_WIN_FENCE(0,win_rho, info)
    call MPI_WIN_FENCE(0,win_rhou,info)
    call MPI_WIN_FENCE(0,win_rhov,info)
    if (.not.is_2D) call MPI_WIN_FENCE(0,win_rhow,info)
    call MPI_WIN_FENCE(0,win_rhoe,info)

  end subroutine communication1
 
  !===============================================================================
  subroutine commun1(var,win)
  !===============================================================================
    !> communications of 3D variable (faces) using one-sided RMA PUT
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer, intent(in) :: win
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(inout) :: var
    ! ----------------------------------------------------------------------------
    
    ! put data from imin (dir.1) into ghost cells of target
    call MPI_PUT(var(in1s,1,1),1,type_facei,neighbor(1), &
                disp_target(1),1,type_face_target1,win,info)
 
    ! put data from imax (dir.2) into ghost cells of target
    call MPI_PUT(var(in2s,1,1),1,type_facei,neighbor(2), &
                disp_target(2),1,type_face_target2,win,info)
    
    ! put data from jmin (dir.3) into ghost cells of target
    call MPI_PUT(var(1,in3s,1),1,type_facej,neighbor(3), &
                disp_target(3),1,type_face_target3,win,info)

    ! put data from jmax (dir.4) into ghost cells of target
    call MPI_PUT(var(1,in4s,1),1,type_facej,neighbor(4), &
                disp_target(4),1,type_face_target4,win,info)

    if (.not.is_2D) then
       ! put data from kmin (dir.5) into ghost cells of target
       call MPI_PUT(var(1,1,1),1,type_facek,neighbor(5), &
                disp_target(5),1,type_facek,win,info)

        ! put data from kmax (dir.6) into ghost cells of target
       call MPI_PUT(var(1,1,nz-ngh+1),1,type_facek,neighbor(6), &
                       disp_target(6),1,type_facek,win,info)
    endif
        
  end subroutine commun1
  
  !===============================================================================
  subroutine window_comm1
  !===============================================================================
    !> Set/Check target windows (MPI_WIN_FENCE)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    ! start/check
    call MPI_WIN_FENCE(0,win_rho, info)
    call MPI_WIN_FENCE(0,win_rhou,info)
    call MPI_WIN_FENCE(0,win_rhov,info)
    if (.not.is_2D) call MPI_WIN_FENCE(0,win_rhow,info)
    call MPI_WIN_FENCE(0,win_rhoe,info)

  end subroutine window_comm1
 
  !===============================================================================
  subroutine communication1_i
  !===============================================================================
    !> Call one-sided RMA communications in i-direction only
    !> [1/ set window; 2/ comm; 3/ check]
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    ! start check
    call MPI_WIN_FENCE(0,win_rho, info)
    call MPI_WIN_FENCE(0,win_rhou,info)
    call MPI_WIN_FENCE(0,win_rhov,info)
    if (.not.is_2D) call MPI_WIN_FENCE(0,win_rhow,info)
    call MPI_WIN_FENCE(0,win_rhoe,info)

    ! start commnunications
    call commun1_i(rho_n, win_rho)
    call commun1_i(rhou_n,win_rhou)
    call commun1_i(rhov_n,win_rhov)
    if (.not.is_2D) call commun1_i(rhow_n,win_rhow)
    call commun1_i(rhoe_n,win_rhoe)

    ! check end of commnunication
    call MPI_WIN_FENCE(0,win_rho, info)
    call MPI_WIN_FENCE(0,win_rhou,info)
    call MPI_WIN_FENCE(0,win_rhov,info)
    if (.not.is_2D) call MPI_WIN_FENCE(0,win_rhow,info)
    call MPI_WIN_FENCE(0,win_rhoe,info)

  end subroutine communication1_i
 
  !===============================================================================
  subroutine communication_1_i
  !===============================================================================
    !> Call one-sided RMA communications in i-direction only (comm. only)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    ! perform communications
    call commun1_i(rho_n, win_rho)
    call commun1_i(rhou_n,win_rhou)
    call commun1_i(rhov_n,win_rhov)
    if (.not.is_2D) call commun1_i(rhow_n,win_rhow)
    call commun1_i(rhoe_n,win_rhoe)

  end subroutine communication_1_i
 
  !===============================================================================
  subroutine commun1_i(var,win)
  !===============================================================================
    !> imin/imax communications using one-sided RMA PUT
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer, intent(in) :: win
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(inout) :: var
    ! ----------------------------------------------------------------------------
        
    ! put data from imin (dir.1) into ghost cells of target
    call MPI_PUT(var(in1s,1,1),1,type_facei,neighbor(1), &
                disp_target(1),1,type_face_target1,win,info)
    
    ! put data from imax (dir.2) into ghost cells of target
    call MPI_PUT(var(in2s,1,1),1,type_facei,neighbor(2), &
                disp_target(2),1,type_face_target2,win,info)
    
  end subroutine commun1_i
  
  !===============================================================================
  subroutine communication1_j
  !===============================================================================
    !> Call one-sided RMA communications in j-direction only
    !> [1/ set window; 2/ comm; 3/ check]
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    ! start check
    call MPI_WIN_FENCE(0,win_rho, info)
    call MPI_WIN_FENCE(0,win_rhou,info)
    call MPI_WIN_FENCE(0,win_rhov,info)
    if (.not.is_2D) call MPI_WIN_FENCE(0,win_rhow,info)
    call MPI_WIN_FENCE(0,win_rhoe,info)

    ! start commnunications
    call commun1_j(rho_n, win_rho)
    call commun1_j(rhou_n,win_rhou)
    call commun1_j(rhov_n,win_rhov)
    if (.not.is_2D) call commun1_j(rhow_n,win_rhow)
    call commun1_j(rhoe_n,win_rhoe)

    ! check end of commnunication
    call MPI_WIN_FENCE(0,win_rho, info)
    call MPI_WIN_FENCE(0,win_rhou,info)
    call MPI_WIN_FENCE(0,win_rhov,info)
    if (.not.is_2D) call MPI_WIN_FENCE(0,win_rhow,info)
    call MPI_WIN_FENCE(0,win_rhoe,info)

  end subroutine communication1_j
  
  !===============================================================================
  subroutine communication_1_j
  !===============================================================================
    !> Call one-sided RMA communications in j-direction only (comm. only)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    ! perform communications
    call commun1_j(rho_n, win_rho)
    call commun1_j(rhou_n,win_rhou)
    call commun1_j(rhov_n,win_rhov)
    if (.not.is_2D) call commun1_j(rhow_n,win_rhow)
    call commun1_j(rhoe_n,win_rhoe)

  end subroutine communication_1_j
 
  !===============================================================================
  subroutine commun1_j(var,win)
  !===============================================================================
    !> jmin/jmax communications using one-sided RMA PUT
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer, intent(in) :: win
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(inout) :: var
    ! ----------------------------------------------------------------------------
        
    ! put data from jmin (dir.3) into ghost cells of target
    call MPI_PUT(var(1,in3s,1),1,type_facej,neighbor(3), &
                disp_target(3),1,type_face_target3,win,info)

    ! put data from jmax (dir.4) into ghost cells of target
    call MPI_PUT(var(1,in4s,1),1,type_facej,neighbor(4), &
                disp_target(4),1,type_face_target4,win,info)

  end subroutine commun1_j
  
  !===============================================================================
  subroutine communication1_k
  !===============================================================================
    !> Call one-sided RMA communications in k-direction only
    !> [1/ set window; 2/ comm; 3/ check]
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    ! start check
    call MPI_WIN_FENCE(0,win_rho, info)
    call MPI_WIN_FENCE(0,win_rhou,info)
    call MPI_WIN_FENCE(0,win_rhov,info)
    if (.not.is_2D) call MPI_WIN_FENCE(0,win_rhow,info)
    call MPI_WIN_FENCE(0,win_rhoe,info)

    ! start commnunications
    call commun1_k(rho_n, win_rho)
    call commun1_k(rhou_n,win_rhou)
    call commun1_k(rhov_n,win_rhov)
    if (.not.is_2D) call commun1_k(rhow_n,win_rhow)
    call commun1_k(rhoe_n,win_rhoe)

    ! check end of commnunication
    call MPI_WIN_FENCE(0,win_rho, info)
    call MPI_WIN_FENCE(0,win_rhou,info)
    call MPI_WIN_FENCE(0,win_rhov,info)
    if (.not.is_2D) call MPI_WIN_FENCE(0,win_rhow,info)
    call MPI_WIN_FENCE(0,win_rhoe,info)

  end subroutine communication1_k
 
  !===============================================================================
  subroutine communication_1_k
  !===============================================================================
    !> Call one-sided RMA communications in k-direction only (comm. only)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    ! perform communications
    call commun1_k(rho_n, win_rho)
    call commun1_k(rhou_n,win_rhou)
    call commun1_k(rhov_n,win_rhov)
    if (.not.is_2D) call commun1_k(rhow_n,win_rhow)
    call commun1_k(rhoe_n,win_rhoe)

  end subroutine communication_1_k
 
  !===============================================================================
  subroutine commun1_k(var,win)
  !===============================================================================
    !> kmin/kmax communications using one-sided RMA PUT
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer, intent(in) :: win
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(inout) :: var
    ! ----------------------------------------------------------------------------
        
    ! put data from kmin (dir.5) into ghost cells of target
    call MPI_PUT(var(1,1,1),1,type_facek,neighbor(5), &
             disp_target(5),1,type_facek,win,info)

    ! put data from kmax (dir.6) into ghost cells of target
    call MPI_PUT(var(1,1,nz-ngh+1),1,type_facek,neighbor(6), &
                    disp_target(6),1,type_facek,win,info)
       
  end subroutine commun1_k
  
  !===============================================================================
  !===============================================================================
  ! II - Communications of faces for increments (extended to ngh_irs ghost points)
  !===============================================================================
  !===============================================================================
  
  !===============================================================================
  subroutine mpi_win_comm1_inc
  !===============================================================================
    !> Window for one-sided RMA communications of increments (IRS method)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: nex,ney,nez ! extended sizes
    integer(KIND=MPI_ADDRESS_KIND) :: win_buffer_size
    integer :: infos
    integer :: group,ningroup,cpt,i
    integer, dimension(:), allocatable :: ranks
    ! ----------------------------------------------------------------------------

    ! Extended sizes (+ ghost cells)
    ! ------------------------------
    nex=nx+ngh_irs(1)+ngh_irs(2)
    ney=ny+ngh_irs(3)+ngh_irs(4)
    if (is_2D) then
       nez=1
    else
       nez=nz+ngh_irs(5)+ngh_irs(6)
    endif

    ! Create window for increments
    ! ----------------------------
    win_buffer_size=nex*ney*nez*sizeofreal
!!$    call MPI_WIN_CREATE(Krho ,win_buffer_size,sizeofreal,MPI_INFO_NULL,COMM_global,win_Krho ,info)
!!$    call MPI_WIN_CREATE(Krhou,win_buffer_size,sizeofreal,MPI_INFO_NULL,COMM_global,win_Krhou,info)
!!$    call MPI_WIN_CREATE(Krhov,win_buffer_size,sizeofreal,MPI_INFO_NULL,COMM_global,win_Krhov,info)
!!$    if (.not.is_2D) call MPI_WIN_CREATE(Krhow,win_buffer_size,sizeofreal,MPI_INFO_NULL,COMM_global,win_Krhow,info)
!!$    call MPI_WIN_CREATE(Krhoe,win_buffer_size,sizeofreal,MPI_INFO_NULL,COMM_global,win_Krhoe,info)

    call MPI_INFO_CREATE(infos,info)
    call MPI_INFO_SET(infos,'no_locks','true',info)
    call MPI_WIN_CREATE(Krho ,win_buffer_size,sizeofreal,infos,COMM_global,win_Krho ,info)
    call MPI_WIN_CREATE(Krhou,win_buffer_size,sizeofreal,infos,COMM_global,win_Krhou,info)
    call MPI_WIN_CREATE(Krhov,win_buffer_size,sizeofreal,infos,COMM_global,win_Krhov,info)
    if (.not.is_2D) call MPI_WIN_CREATE(Krhow,win_buffer_size,sizeofreal,infos,COMM_global,win_Krhow,info)
    call MPI_WIN_CREATE(Krhoe,win_buffer_size,sizeofreal,infos,COMM_global,win_Krhoe,info)
    call MPI_INFO_FREE(infos,info)
    
    ! Create groups
    ! ----------------------------
    ! Krho
    call MPI_WIN_GET_GROUP(win_Krho,group,info)
    ningroup=0
    do i=1,4
       if (neighbor(i)/=MPI_PROC_NULL) ningroup=ningroup+1
    enddo
    allocate(ranks(ningroup))
    cpt=1
    do i=1,4
       if (neighbor(i)/=MPI_PROC_NULL) then
          ranks(cpt)=neighbor(i)
          cpt=cpt+1
       endif
    enddo
    call MPI_GROUP_INCL(group,ningroup,ranks,group_Krho,info)
    call MPI_GROUP_FREE(group,info)
    deallocate(ranks)
    ! Krhou
    call MPI_WIN_GET_GROUP(win_Krhou,group,info)
    ningroup=0
    do i=1,4
       if (neighbor(i)/=MPI_PROC_NULL) ningroup=ningroup+1
    enddo
    allocate(ranks(ningroup))
    cpt=1
    do i=1,4
       if (neighbor(i)/=MPI_PROC_NULL) then
          ranks(cpt)=neighbor(i)
          cpt=cpt+1
       endif
    enddo
    call MPI_GROUP_INCL(group,ningroup,ranks,group_Krhou,info)
    call MPI_GROUP_FREE(group,info)
    deallocate(ranks)
    ! Krhov
    call MPI_WIN_GET_GROUP(win_Krhov,group,info)
    ningroup=0
    do i=1,4
       if (neighbor(i)/=MPI_PROC_NULL) ningroup=ningroup+1
    enddo
    allocate(ranks(ningroup))
    cpt=1
    do i=1,4
       if (neighbor(i)/=MPI_PROC_NULL) then
          ranks(cpt)=neighbor(i)
          cpt=cpt+1
       endif
    enddo
    call MPI_GROUP_INCL(group,ningroup,ranks,group_Krhov,info)
    call MPI_GROUP_FREE(group,info)
    deallocate(ranks)
    ! Krhow
    if (.not.is_2D) then
       call MPI_WIN_GET_GROUP(win_Krhow,group,info)
       ningroup=0
       do i=1,4
          if (neighbor(i)/=MPI_PROC_NULL) ningroup=ningroup+1
       enddo
       allocate(ranks(ningroup))
       cpt=1
       do i=1,4
          if (neighbor(i)/=MPI_PROC_NULL) then
             ranks(cpt)=neighbor(i)
             cpt=cpt+1
          endif
       enddo
       call MPI_GROUP_INCL(group,ningroup,ranks,group_Krhow,info)
       call MPI_GROUP_FREE(group,info)
       deallocate(ranks)
    endif
    ! Krhoe
    call MPI_WIN_GET_GROUP(win_Krhoe,group,info)
    ningroup=0
    do i=1,4
       if (neighbor(i)/=MPI_PROC_NULL) ningroup=ningroup+1
    enddo
    allocate(ranks(ningroup))
    cpt=1
    do i=1,4
       if (neighbor(i)/=MPI_PROC_NULL) then
          ranks(cpt)=neighbor(i)
          cpt=cpt+1
       endif
    enddo
    call MPI_GROUP_INCL(group,ningroup,ranks,group_Krhoe,info)
    call MPI_GROUP_FREE(group,info)
    deallocate(ranks)
    
  end subroutine mpi_win_comm1_inc
 
  !===============================================================================
  subroutine mpi_close_win_comm1_inc
  !===============================================================================
    !> CLose open MPI memory windows used for one-sided communications
    !> of increments (IRS method)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------
    
    ! close open memory windows
    call MPI_WIN_FREE(win_Krho,info)
    call MPI_WIN_FREE(win_Krhou,info)
    call MPI_WIN_FREE(win_Krhov,info)
    if (.not.is_2D) call MPI_WIN_FREE(win_Krhow,info)
    call MPI_WIN_FREE(win_Krhoe,info)

  end subroutine mpi_close_win_comm1_inc
 
  !===============================================================================
  subroutine communication1_inc
  !===============================================================================
    !> Call one-sided RMA communications for increments (Krho,...)
    !> [1/ set window; 2/ comm; 3/ check]
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------
    
    ! start check
    call MPI_WIN_FENCE(0,win_Krho, info)
    call MPI_WIN_FENCE(0,win_Krhou,info)
    call MPI_WIN_FENCE(0,win_Krhov,info)
    if (.not.is_2D) call MPI_WIN_FENCE(0,win_Krhow,info)
    call MPI_WIN_FENCE(0,win_Krhoe,info)

    ! start commnunications
    call commun1_inc(Krho,win_Krho)
    call commun1_inc(Krhou,win_Krhou)
    call commun1_inc(Krhov,win_Krhov)
    if (.not.is_2D) call commun1_inc(Krhow,win_Krhow)
    call commun1_inc(Krhoe,win_Krhoe)

    ! check end of commnunication
    call MPI_WIN_FENCE(0,win_Krho, info)
    call MPI_WIN_FENCE(0,win_Krhou,info)
    call MPI_WIN_FENCE(0,win_Krhov,info)
    if (.not.is_2D) call MPI_WIN_FENCE(0,win_Krhow,info)
    call MPI_WIN_FENCE(0,win_Krhoe,info)

  end subroutine communication1_inc
 
  !===============================================================================
  subroutine communication_1_inc
  !===============================================================================
    !> Call one-sided RMA communications for increments (Krho,...)
    !> [1/ set window; 2/ comm; 3/ check]
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------
    
    ! start commnunications
    call commun1_inc(Krho,win_Krho)
    call commun1_inc(Krhou,win_Krhou)
    call commun1_inc(Krhov,win_Krhov)
    if (.not.is_2D) call commun1_inc(Krhow,win_Krhow)
    call commun1_inc(Krhoe,win_Krhoe)

  end subroutine communication_1_inc
 
  !===============================================================================
  subroutine commun1_inc(var,win)
  !===============================================================================
    !> communications of 3D variable (increments) using one-sided RMA PUT
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer, intent(in) :: win
    real(wp), dimension(nx1_irs:nx2_irs,ny1_irs:ny2_irs,nz1_irs:nz2_irs), intent(inout) :: var
    ! ----------------------------------------------------------------------------
    
    ! put data from imin (dir.1) into ghost cells of target
    call MPI_PUT(var(ii1s,1,1),1,type_inc1,neighbor(1), &
            disp_inc_target(1),1,type_inc_target1,win,info)
    
    ! put data from imax (dir.2) into ghost cells of target
    call MPI_PUT(var(ii2s,1,1),1,type_inc2,neighbor(2), &
            disp_inc_target(2),1,type_inc_target2,win,info)
    
    ! put data from jmin (dir.3) into ghost cells of target
    call MPI_PUT(var(1,ii3s,1),1,type_inc3,neighbor(3), &
            disp_inc_target(3),1,type_inc_target3,win,info)

    ! put data from jmax (dir.4) into ghost cells of target
    call MPI_PUT(var(1,ii4s,1),1,type_inc4,neighbor(4), &
            disp_inc_target(4),1,type_inc_target4,win,info)

    if (.not.is_2D) then
       ! put data from kmin (dir.5) into ghost cells of target
       call MPI_PUT(var(1,1,ii5s),1,type_inc5,neighbor(5), &
               disp_inc_target(5),1,type_inc5,win,info)

        ! put data from kmax (dir.6) into ghost cells of target
       call MPI_PUT(var(1,1,ii6s),1,type_inc6,neighbor(6), &
               disp_inc_target(6),1,type_inc6,win,info)
    endif
    
  end subroutine commun1_inc
  
  !===============================================================================
  subroutine window_comm1_inc(K_rho,K_rhou,K_rhov,K_rhow,K_rhoe)
  !===============================================================================
    !> Set/Check target windows for increments (MPI_WIN_FENCE)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp) :: K_rho(*),K_rhou(*),K_rhov(*),K_rhow(*),K_rhoe(*)
    ! ----------------------------------------------------------------------------

    ! start/check
    call MPI_WIN_FENCE(0,win_Krho, info)
    call MPI_WIN_FENCE(0,win_Krhou,info)
    call MPI_WIN_FENCE(0,win_Krhov,info)
    if (.not.is_2D) call MPI_WIN_FENCE(0,win_Krhow,info)
    call MPI_WIN_FENCE(0,win_Krhoe,info)

  end subroutine window_comm1_inc
 
  !===============================================================================
  subroutine communication1_inc_ij
  !===============================================================================
    !> Call one-sided RMA communications in i-direction only
    !> [1/ set window; 2/ comm; 3/ check]
    !> for increment (ngh_irs ghost points)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    ! start check
    call MPI_WIN_FENCE(0,win_Krho, info)
    call MPI_WIN_FENCE(0,win_Krhou,info)
    call MPI_WIN_FENCE(0,win_Krhov,info)
    if (.not.is_2D) call MPI_WIN_FENCE(0,win_Krhow,info)
    call MPI_WIN_FENCE(0,win_Krhoe,info)

    ! start commnunications
    call commun1_inc_ij(Krho,win_Krho)
    call commun1_inc_ij(Krhou,win_Krhou)
    call commun1_inc_ij(Krhov,win_Krhov)
    if (.not.is_2D) call commun1_inc_i(Krhow,win_Krhow)
    call commun1_inc_ij(Krhoe,win_Krhoe)

    ! check end of commnunication
    call MPI_WIN_FENCE(0,win_Krho, info)
    call MPI_WIN_FENCE(0,win_Krhou,info)
    call MPI_WIN_FENCE(0,win_Krhov,info)
    if (.not.is_2D) call MPI_WIN_FENCE(0,win_Krhow,info)
    call MPI_WIN_FENCE(0,win_Krhoe,info)

!!$    ! 3 times SLOWER !! / 2times with MPI_MODE_NOCHECK / POST,START,COMPLETE,WAIT PSCW 1.5 times SLOWER
!!$    ! start check
!!$    call MPI_WIN_POST (group_Krho,MPI_MODE_NOSTORE+MPI_MODE_NOPUT+MPI_MODE_NOCHECK,win_Krho,info)
!!$    call MPI_WIN_POST (group_Krhou,MPI_MODE_NOSTORE+MPI_MODE_NOPUT+MPI_MODE_NOCHECK,win_Krhou,info)
!!$    call MPI_WIN_POST (group_Krhov,MPI_MODE_NOSTORE+MPI_MODE_NOPUT+MPI_MODE_NOCHECK,win_Krhov,info)
!!$    call MPI_WIN_POST (group_Krhoe,MPI_MODE_NOSTORE+MPI_MODE_NOPUT+MPI_MODE_NOCHECK,win_Krhoe,info)
!!$    call MPI_WIN_START(group_Krho,MPI_MODE_NOCHECK,win_Krho,info)
!!$    call MPI_WIN_START(group_Krhou,MPI_MODE_NOCHECK,win_Krhou,info)
!!$    call MPI_WIN_START(group_Krhov,MPI_MODE_NOCHECK,win_Krhov,info)
!!$    call MPI_WIN_START(group_Krhoe,MPI_MODE_NOCHECK,win_Krhoe,info)
!!$    
!!$    if (.not.is_2D) then
!!$       call MPI_WIN_POST (group_Krhow,MPI_MODE_NOSTORE+MPI_MODE_NOPUT+MPI_MODE_NOCHECK,win_Krhow,info)
!!$       call MPI_WIN_START(group_Krhow,MPI_MODE_NOCHECK,win_Krhow,info)
!!$    endif
!!$
!!$    ! start commnunications
!!$    call commun1_inc_ij(Krho,win_Krho)
!!$    call commun1_inc_ij(Krhou,win_Krhou)
!!$    call commun1_inc_ij(Krhov,win_Krhov)
!!$    if (.not.is_2D) call commun1_inc_i(Krhow,win_Krhow)
!!$    call commun1_inc_ij(Krhoe,win_Krhoe)
!!$
!!$    ! check end of commnunication
!!$    call MPI_WIN_COMPLETE(win_Krho,info)
!!$    call MPI_WIN_COMPLETE(win_Krhou,info)
!!$    call MPI_WIN_COMPLETE(win_Krhov,info)
!!$    call MPI_WIN_COMPLETE(win_Krhoe,info)
!!$    
!!$    call MPI_WIN_WAIT(win_Krho,info)
!!$    call MPI_WIN_WAIT(win_Krhou,info)
!!$    call MPI_WIN_WAIT(win_Krhov,info)
!!$    call MPI_WIN_WAIT(win_Krhoe,info)
!!$    
!!$    if (.not.is_2D) then
!!$       call MPI_WIN_COMPLETE(win_Krhow,info)
!!$       call MPI_WIN_WAIT(win_Krhow,info)
!!$    endif
    
!!$    ! 5 times SLOWER !!
!!$    ! start check
!!$    call MPI_WIN_POST (group_Krho,0,win_Krho,info)
!!$    call MPI_WIN_START(group_Krho,0,win_Krho,info)
!!$    call commun1_inc_ij(Krho,win_Krho)
!!$    call MPI_WIN_COMPLETE(win_Krho,info)
!!$    call MPI_WIN_WAIT(win_Krho,info)
!!$    
!!$    call MPI_WIN_POST (group_Krhou,0,win_Krhou,info)
!!$    call MPI_WIN_START(group_Krhou,0,win_Krhou,info)
!!$    call commun1_inc_ij(Krhou,win_Krhou)
!!$    call MPI_WIN_COMPLETE(win_Krhou,info)
!!$    call MPI_WIN_WAIT(win_Krhou,info)
!!$    
!!$    call MPI_WIN_POST (group_Krhov,0,win_Krhov,info)
!!$    call MPI_WIN_START(group_Krhov,0,win_Krhov,info)
!!$    call commun1_inc_ij(Krhov,win_Krhov)
!!$    call MPI_WIN_COMPLETE(win_Krhov,info)
!!$    call MPI_WIN_WAIT(win_Krhov,info)
!!$    
!!$    call MPI_WIN_POST (group_Krhoe,0,win_Krhoe,info)
!!$    call MPI_WIN_START(group_Krhoe,0,win_Krhoe,info)
!!$    call commun1_inc_ij(Krhoe,win_Krhoe)
!!$    call MPI_WIN_COMPLETE(win_Krhoe,info)
!!$    call MPI_WIN_WAIT(win_Krhoe,info)
!!$   
!!$    if (.not.is_2D) then
!!$       call MPI_WIN_POST (group_Krhow,0,win_Krhow,info)
!!$       call MPI_WIN_START(group_Krhow,0,win_Krhow,info)
!!$       call commun1_inc_i(Krhow,win_Krhow)
!!$       call MPI_WIN_COMPLETE(win_Krhow,info)
!!$       call MPI_WIN_WAIT(win_Krhow,info)
!!$    endif
    
!!$    ! --> NO GAIN with this version
!!$    ! start check
!!$    call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE,win_Krho, info)
!!$    call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE,win_Krhou,info)
!!$    call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE,win_Krhov,info)
!!$    if (.not.is_2D) call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE,win_Krhow,info)
!!$    call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE,win_Krhoe,info)
!!$
!!$    ! start commnunications
!!$    call commun1_inc_ij(Krho,win_Krho)
!!$    call commun1_inc_ij(Krhou,win_Krhou)
!!$    call commun1_inc_ij(Krhov,win_Krhov)
!!$    if (.not.is_2D) call commun1_inc_i(Krhow,win_Krhow)
!!$    call commun1_inc_ij(Krhoe,win_Krhoe)
!!$
!!$    ! check end of commnunication
!!$    call MPI_WIN_FENCE(MPI_MODE_NOSTORE+MPI_MODE_NOPUT+MPI_MODE_NOSUCCEED,win_Krho, info)
!!$    call MPI_WIN_FENCE(MPI_MODE_NOSTORE+MPI_MODE_NOPUT+MPI_MODE_NOSUCCEED,win_Krhou,info)
!!$    call MPI_WIN_FENCE(MPI_MODE_NOSTORE+MPI_MODE_NOPUT+MPI_MODE_NOSUCCEED,win_Krhov,info)
!!$    if (.not.is_2D) call MPI_WIN_FENCE(MPI_MODE_NOSTORE+MPI_MODE_NOPUT+MPI_MODE_NOSUCCEED,win_Krhow,info)
!!$    call MPI_WIN_FENCE(MPI_MODE_NOSTORE+MPI_MODE_NOPUT+MPI_MODE_NOSUCCEED,win_Krhoe,info)


!!$    ! --> this version is four times SLOWER !!!!
!!$    ! start/comm/fullfill
!!$    call MPI_WIN_FENCE(0,win_Krho, info)
!!$    call commun1_inc_ij(Krho,win_Krho)
!!$    call MPI_WIN_FENCE(0,win_Krho, info)
!!$    
!!$    call MPI_WIN_FENCE(0,win_Krhou,info)
!!$    call commun1_inc_ij(Krhou,win_Krhou)
!!$    call MPI_WIN_FENCE(0,win_Krhou,info)
!!$    
!!$    call MPI_WIN_FENCE(0,win_Krhov,info)
!!$    call commun1_inc_ij(Krhov,win_Krhov)
!!$    call MPI_WIN_FENCE(0,win_Krhov,info)
!!$
!!$    call MPI_WIN_FENCE(0,win_Krhoe,info)
!!$    call commun1_inc_ij(Krhoe,win_Krhoe)
!!$    call MPI_WIN_FENCE(0,win_Krhoe,info)
!!$    
!!$    if (.not.is_2D) then
!!$       call MPI_WIN_FENCE(0,win_Krhow,info)
!!$       call commun1_inc_i(Krhow,win_Krhow)
!!$       call MPI_WIN_FENCE(0,win_Krhow,info)
!!$    endif

  end subroutine communication1_inc_ij
 
  !===============================================================================
  subroutine communication_1_inc_ij
  !===============================================================================
    !> Call one-sided RMA communications in i-direction only (comm. only)
    !> for increment (ngh_irs ghost points)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    ! perform communications
    call commun1_inc_ij(Krho,win_Krho)
    call commun1_inc_ij(Krhou,win_Krhou)
    call commun1_inc_ij(Krhov,win_Krhov)
    if (.not.is_2D) call commun1_inc_i(Krhow,win_Krhow)
    call commun1_inc_i(Krhoe,win_Krhoe)

  end subroutine communication_1_inc_ij
 
  !===============================================================================
  subroutine commun1_inc_ij(var,win)
  !===============================================================================
    !> imin/imax communications using one-sided RMA PUT
    !> for increment (ngh_irs ghost points)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer, intent(in) :: win
    real(wp), dimension(nx1_irs:nx2_irs,ny1_irs:ny2_irs,nz1_irs:nz2_irs), intent(inout) :: var
    ! ----------------------------------------------------------------------------
        
    ! put data from imin (dir.1) into ghost cells of target
    call MPI_PUT(var(ii1s,1,1),1,type_inc1,neighbor(1), &
            disp_inc_target(1),1,type_inc_target1,win,info)
    
    ! put data from imax (dir.2) into ghost cells of target
    call MPI_PUT(var(ii2s,1,1),1,type_inc2,neighbor(2), &
            disp_inc_target(2),1,type_inc_target2,win,info)
    
    ! put data from jmin (dir.3) into ghost cells of target
    call MPI_PUT(var(1,ii3s,1),1,type_inc3,neighbor(3), &
            disp_inc_target(3),1,type_inc_target3,win,info)

    ! put data from jmax (dir.4) into ghost cells of target
    call MPI_PUT(var(1,ii4s,1),1,type_inc4,neighbor(4), &
            disp_inc_target(4),1,type_inc_target4,win,info)

  end subroutine commun1_inc_ij
  
  !===============================================================================
  subroutine communication1_inc_i
  !===============================================================================
    !> Call one-sided RMA communications in i-direction only
    !> [1/ set window; 2/ comm; 3/ check]
    !> for increment (ngh_irs ghost points)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    ! start check
    call MPI_WIN_FENCE(0,win_Krho, info)
    call MPI_WIN_FENCE(0,win_Krhou,info)
    call MPI_WIN_FENCE(0,win_Krhov,info)
    if (.not.is_2D) call MPI_WIN_FENCE(0,win_Krhow,info)
    call MPI_WIN_FENCE(0,win_Krhoe,info)

    ! start commnunications
    call commun1_inc_i(Krho,win_Krho)
    call commun1_inc_i(Krhou,win_Krhou)
    call commun1_inc_i(Krhov,win_Krhov)
    if (.not.is_2D) call commun1_inc_i(Krhow,win_Krhow)
    call commun1_inc_i(Krhoe,win_Krhoe)

    ! check end of commnunication
    call MPI_WIN_FENCE(0,win_Krho, info)
    call MPI_WIN_FENCE(0,win_Krhou,info)
    call MPI_WIN_FENCE(0,win_Krhov,info)
    if (.not.is_2D) call MPI_WIN_FENCE(0,win_Krhow,info)
    call MPI_WIN_FENCE(0,win_Krhoe,info)

  end subroutine communication1_inc_i
 
  !===============================================================================
  subroutine communication_1_inc_i
  !===============================================================================
    !> Call one-sided RMA communications in i-direction only (comm. only)
    !> for increment (ngh_irs ghost points)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    ! perform communications
    call commun1_inc_i(Krho,win_Krho)
    call commun1_inc_i(Krhou,win_Krhou)
    call commun1_inc_i(Krhov,win_Krhov)
    if (.not.is_2D) call commun1_inc_i(Krhow,win_Krhow)
    call commun1_inc_i(Krhoe,win_Krhoe)

  end subroutine communication_1_inc_i
 
  !===============================================================================
  subroutine commun1_inc_i(var,win)
  !===============================================================================
    !> imin/imax communications using one-sided RMA PUT
    !> for increment (ngh_irs ghost points)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer, intent(in) :: win
    real(wp), dimension(nx1_irs:nx2_irs,ny1_irs:ny2_irs,nz1_irs:nz2_irs), intent(inout) :: var
    ! ----------------------------------------------------------------------------
        
    ! put data from imin (dir.1) into ghost cells of target
    call MPI_PUT(var(ii1s,1,1),1,type_inc1,neighbor(1), &
            disp_inc_target(1),1,type_inc_target1,win,info)
    
    ! put data from imax (dir.2) into ghost cells of target
    call MPI_PUT(var(ii2s,1,1),1,type_inc2,neighbor(2), &
            disp_inc_target(2),1,type_inc_target2,win,info)
    
  end subroutine commun1_inc_i
  
  !===============================================================================
  subroutine communication1_inc_j
  !===============================================================================
    !> Call one-sided RMA communications in j-direction only
    !> [1/ set window; 2/ comm; 3/ check]
    !> for increment (ngh_irs ghost points)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    ! start check
    call MPI_WIN_FENCE(0,win_Krho, info)
    call MPI_WIN_FENCE(0,win_Krhou,info)
    call MPI_WIN_FENCE(0,win_Krhov,info)
    if (.not.is_2D) call MPI_WIN_FENCE(0,win_Krhow,info)
    call MPI_WIN_FENCE(0,win_Krhoe,info)

    ! start commnunications
    call commun1_inc_j(Krho,win_Krho)
    call commun1_inc_j(Krhou,win_Krhou)
    call commun1_inc_j(Krhov,win_Krhov)
    if (.not.is_2D) call commun1_inc_j(Krhow,win_Krhow)
    call commun1_inc_j(Krhoe,win_Krhoe)

    ! check end of commnunication
    call MPI_WIN_FENCE(0,win_Krho, info)
    call MPI_WIN_FENCE(0,win_Krhou,info)
    call MPI_WIN_FENCE(0,win_Krhov,info)
    if (.not.is_2D) call MPI_WIN_FENCE(0,win_Krhow,info)
    call MPI_WIN_FENCE(0,win_Krhoe,info)

  end subroutine communication1_inc_j
 
  !===============================================================================
  subroutine communication_1_inc_j
  !===============================================================================
    !> Call one-sided RMA communications in j-direction only (comm. only)
    !> for increment (ngh_irs ghost points)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    ! perform communications
    call commun1_inc_j(Krho,win_Krho)
    call commun1_inc_j(Krhou,win_Krhou)
    call commun1_inc_j(Krhov,win_Krhov)
    if (.not.is_2D) call commun1_inc_j(Krhow,win_Krhow)
    call commun1_inc_j(Krhoe,win_Krhoe)

  end subroutine communication_1_inc_j
 
  !===============================================================================
  subroutine commun1_inc_j(var,win)
  !===============================================================================
    !> jmin/jmax communications using one-sided RMA PUT
    !> for increment (ngh_irs ghost points)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer, intent(in) :: win
    real(wp), dimension(nx1_irs:nx2_irs,ny1_irs:ny2_irs,nz1_irs:nz2_irs), intent(inout) :: var
    ! ----------------------------------------------------------------------------
        
    ! put data from jmin (dir.3) into ghost cells of target
    call MPI_PUT(var(1,ii3s,1),1,type_inc3,neighbor(3), &
            disp_inc_target(3),1,type_inc_target3,win,info)

    ! put data from jmax (dir.4) into ghost cells of target
    call MPI_PUT(var(1,ii4s,1),1,type_inc4,neighbor(4), &
            disp_inc_target(4),1,type_inc_target4,win,info)

  end subroutine commun1_inc_j
  
  !===============================================================================
  subroutine communication1_inc_k
  !===============================================================================
    !> Call one-sided RMA communications in k-direction only
    !> [1/ set window; 2/ comm; 3/ check]
    !> for increment (ngh_irs ghost points)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    ! start check
    call MPI_WIN_FENCE(0,win_Krho, info)
    call MPI_WIN_FENCE(0,win_Krhou,info)
    call MPI_WIN_FENCE(0,win_Krhov,info)
    if (.not.is_2D) call MPI_WIN_FENCE(0,win_Krhow,info)
    call MPI_WIN_FENCE(0,win_Krhoe,info)

    ! start commnunications
    call commun1_inc_k(Krho,win_Krho)
    call commun1_inc_k(Krhou,win_Krhou)
    call commun1_inc_k(Krhov,win_Krhov)
    if (.not.is_2D) call commun1_inc_k(Krhow,win_Krhow)
    call commun1_inc_k(Krhoe,win_Krhoe)

    ! check end of commnunication
    call MPI_WIN_FENCE(0,win_Krho, info)
    call MPI_WIN_FENCE(0,win_Krhou,info)
    call MPI_WIN_FENCE(0,win_Krhov,info)
    if (.not.is_2D) call MPI_WIN_FENCE(0,win_Krhow,info)
    call MPI_WIN_FENCE(0,win_Krhoe,info)

  end subroutine communication1_inc_k
 
  !===============================================================================
  subroutine communication_1_inc_k
  !===============================================================================
    !> Call one-sided RMA communications in k-direction only (comm. only)
    !> for increment (ngh_irs ghost points)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    ! perform communications
    call commun1_inc_k(Krho,win_Krho)
    call commun1_inc_k(Krhou,win_Krhou)
    call commun1_inc_k(Krhov,win_Krhov)
    if (.not.is_2D) call commun1_inc_k(Krhow,win_Krhow)
    call commun1_inc_k(Krhoe,win_Krhoe)

  end subroutine communication_1_inc_k
 
  !===============================================================================
  subroutine commun1_inc_k(var,win)
  !===============================================================================
    !> kmin/kmax communications using one-sided RMA PUT
    !> for increment (ngh_irs ghost points)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer, intent(in) :: win
    real(wp), dimension(nx1_irs:nx2_irs,ny1_irs:ny2_irs,nz1_irs:nz2_irs), intent(inout) :: var
    ! ----------------------------------------------------------------------------
        
    ! put data from kmin (dir.5) into ghost cells of target
    call MPI_PUT(var(1,1,ii5s),1,type_inc5,neighbor(5), &
            disp_inc_target(5),1,type_inc5,win,info)

    ! put data from kmax (dir.6) into ghost cells of target
    call MPI_PUT(var(1,1,ii6s),1,type_inc6,neighbor(6), &
            disp_inc_target(6),1,type_inc6,win,info)
       
  end subroutine commun1_inc_k
 
  !===============================================================================
  !===============================================================================
  ! III - Communications of faces for CFL (extended to ngh_irs ghost points)
  !===============================================================================
  !===============================================================================
  
  !===============================================================================
  subroutine mpi_win_comm1_cfl
  !===============================================================================
    !> Window for one-sided RMA communications of increments (IRS method)
  !===============================================================================
    use mod_time
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: nex,ney,nez ! extended sizes
    integer(KIND=MPI_ADDRESS_KIND) :: win_buffer_size
    integer :: infos
    ! ----------------------------------------------------------------------------

    ! Extended sizes (+ ghost cells)
    ! ------------------------------
    nex=nx+ngh_irs(1)+ngh_irs(2)
    ney=ny+ngh_irs(3)+ngh_irs(4)
    if (is_2D) then
       nez=1
    else
       nez=nz+ngh_irs(5)+ngh_irs(6)
    endif

    ! Create window for increments
    ! ----------------------------
    win_buffer_size=nex*ney*nez*sizeofreal
    call MPI_INFO_CREATE(infos,info)
    call MPI_INFO_SET(infos,'no_locks','true',info)
    !call MPI_WIN_CREATE(cfl_l ,win_buffer_size,sizeofreal,MPI_INFO_NULL,COMM_global,win_cfl ,info)
    call MPI_WIN_CREATE(cfl_l ,win_buffer_size,sizeofreal,infos,COMM_global,win_cfl ,info)
    call MPI_INFO_FREE(infos,info)

  end subroutine mpi_win_comm1_cfl
 
  !===============================================================================
  subroutine mpi_close_win_comm1_cfl
  !===============================================================================
    !> CLose open MPI memory windows used for one-sided communications
    !> of increments (IRS method)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------
    
    ! close open memory windows
    call MPI_WIN_FREE(win_cfl,info)

  end subroutine mpi_close_win_comm1_cfl
 
  !===============================================================================
  subroutine window_comm1_cfl
  !===============================================================================
    !> Set/Check target windows for increments (MPI_WIN_FENCE)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    ! start/check
    call MPI_WIN_FENCE(0,win_cfl,info)

  end subroutine window_comm1_cfl
 
  !===============================================================================
  subroutine communication1_cfl_ij
  !===============================================================================
    !> Call one-sided RMA communications in i-direction only
    !> [1/ set window; 2/ comm; 3/ check]
    !> for CFL (ngh_irs ghost points)
  !===============================================================================
    use mod_time
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    ! start check
    call MPI_WIN_FENCE(0,win_cfl,info)

    ! start commnunications
    call commun1_inc_ij(cfl_l,win_cfl)

    ! check end of commnunication
    call MPI_WIN_FENCE(0,win_cfl,info)

  end subroutine communication1_cfl_ij
 
  !===============================================================================
  subroutine communication1_cfl_i
  !===============================================================================
    !> Call one-sided RMA communications in i-direction only
    !> [1/ set window; 2/ comm; 3/ check]
    !> for CFL (ngh_irs ghost points)
  !===============================================================================
    use mod_time
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    ! start check
    call MPI_WIN_FENCE(0,win_cfl,info)

    ! start commnunications
    call commun1_inc_i(cfl_l,win_cfl)

    ! check end of commnunication
    call MPI_WIN_FENCE(0,win_cfl, info)

  end subroutine communication1_cfl_i
 
  !===============================================================================
  subroutine communication_1_cfl_i
  !===============================================================================
    !> Call one-sided RMA communications in i-direction only (comm. only)
    !> for CFL (ngh_irs ghost points)
  !===============================================================================
    use mod_time
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    ! perform communications
    call commun1_inc_i(cfl_l,win_cfl)

  end subroutine communication_1_cfl_i
   
  !===============================================================================
  subroutine communication1_cfl_j
  !===============================================================================
    !> Call one-sided RMA communications in j-direction only
    !> [1/ set window; 2/ comm; 3/ check]
    !> for CFL (ngh_irs ghost points)
  !===============================================================================
    use mod_time
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    ! start check
    call MPI_WIN_FENCE(0,win_cfl,info)

    ! start commnunications
    call commun1_inc_j(cfl_l,win_cfl)

    ! check end of commnunication
    call MPI_WIN_FENCE(0,win_cfl,info)

  end subroutine communication1_cfl_j
 
  !===============================================================================
  subroutine communication_1_cfl_j
  !===============================================================================
    !> Call one-sided RMA communications in j-direction only (comm. only)
    !> for CFL (ngh_irs ghost points)
  !===============================================================================
    use mod_time
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    ! perform communications
    call commun1_inc_j(cfl_l,win_cfl)

  end subroutine communication_1_cfl_j
 
  !===============================================================================
  subroutine communication1_cfl_k
  !===============================================================================
    !> Call one-sided RMA communications in k-direction only
    !> [1/ set window; 2/ comm; 3/ check]
    !> for CFL (ngh_irs ghost points)
  !===============================================================================
    use mod_time
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    ! start check
    call MPI_WIN_FENCE(0,win_cfl,info)

    ! start commnunications
    call commun1_inc_k(cfl_l,win_cfl)

    ! check end of commnunication
    call MPI_WIN_FENCE(0,win_cfl,info)

  end subroutine communication1_cfl_k
 
  !===============================================================================
  subroutine communication_1_cfl_k
  !===============================================================================
    !> Call one-sided RMA communications in k-direction only (comm. only)
    !> for CFL (ngh_irs ghost points)
  !===============================================================================
    use mod_time
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    ! perform communications
    call commun1_inc_k(cfl_l,win_cfl)

  end subroutine communication_1_cfl_k
 
end module mod_comm1
  
