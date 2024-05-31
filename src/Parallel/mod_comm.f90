!=================================================================================
module mod_comm
!=================================================================================
  !> Module for non-blocking ISEND/IRECV MPI communications
  !> * "two-sided communications" *
!=================================================================================
  use mod_mpi_types_two_sided
  use mod_flow
  implicit none
  !-------------------------------------------------------------------------------
  ! Communications of unknown vector along faces
  ! ============================================
  ! 2D version
  integer, parameter :: size_2d=8*4 ! 2(send/recv)*4(neighbors)*4(variables)
  integer, dimension(size_2d) :: request_2d
  integer, dimension(MPI_STATUS_SIZE,size_2d) :: status_2d
  ! 3D version
  integer, parameter :: size_3d=12*5 ! 2(send/recv)*6(neighbors)*5(variables)
  integer, dimension(size_3d) :: request_3d
  integer, dimension(MPI_STATUS_SIZE,size_3d) :: status_3d
  !-------------------------------------------------------------------------------
  ! Communications of velocity/temperature derivatives along faces
  ! ==============================================================
  ! 2D version
  integer, parameter :: size_v_2d=8*6 ! 2(send/recv)*4(neighbors)*6(variables)
  integer, dimension(size_v_2d) :: request_v_2d
  integer, dimension(MPI_STATUS_SIZE,size_v_2d) :: status_v_2d
  ! 3D version
  integer, parameter :: size_v=12*12 ! 2(send/recv)*6(neighbors)*12(variables)
  integer, dimension(size_v) :: request_v
  integer, dimension(MPI_STATUS_SIZE,size_v) :: status_v
  !-------------------------------------------------------------------------------
  ! Communications of increments along faces
  ! ========================================
  ! 2D version
  integer, parameter :: size_2di=8*5 ! 2(send/recv)*4(neighbors)*5(variables)
  integer, dimension(size_2di) :: request_2di
  integer, dimension(MPI_STATUS_SIZE,size_2di) :: status_2di
  ! 3D version
  integer, parameter :: size_3di=12*6 ! 2(send/recv)*4(neighbors)*6(variables)
  integer, dimension(size_3di) :: request_3di
  integer, dimension(MPI_STATUS_SIZE,size_3di) :: status_3di
  ! k-direction
  integer, parameter :: size_ik=4*6 ! 2(send/recv)*4(neighbors)*6(variables)
  integer, dimension(size_ik) :: request_ik
  integer, dimension(MPI_STATUS_SIZE,size_ik) :: status_ik
  !-------------------------------------------------------------------------------
  !!!!!!!!!!!!!!!!!!!! Modif for RANS !!!!!!!!!!!!!!!!!!!
  !-------------------------------------------------------------------------------
  ! Communications of one variable along faces
  ! ==========================================
  ! 2D version
  integer, parameter :: size2d=8 ! 2(send/recv)*4(neighbors)*1(variable)
  integer, dimension(size2d) :: request2d
  integer, dimension(MPI_STATUS_SIZE,size2d) :: status2d
  ! 3D version
  integer, parameter :: size3d=12 ! 2(send/recv)*6(neighbors)*1(variable)
  integer, dimension(size3d) :: request3d
  integer, dimension(MPI_STATUS_SIZE,size3d) :: status3d

  ! Communications of turbulent variable derivatives along faces
  ! ============================================================
  ! 2D version comm
  integer, parameter :: size_2d_rans=8*2 ! 2(send/recv)*4(neighbors)*2(variables)
  integer, dimension(size_2d_rans) :: request_2d_rans
  integer, dimension(MPI_STATUS_SIZE,size_2d_rans) :: status_2d_rans
  ! 2D version increments
  integer, parameter :: size_2di_rans=8*1 ! 2(send/recv)*4(neighbors)*1(variables)
  integer, dimension(size_2di_rans) :: request_2di_rans
  integer, dimension(MPI_STATUS_SIZE,size_2di_rans) :: status_2di_rans
  ! 3D version comm
  integer, parameter :: size_3d_rans=12*3 ! 2(send/recv)*6(neighbors)*3(variables)
  integer, dimension(size_3d_rans) :: request_3d_rans
  integer, dimension(MPI_STATUS_SIZE,size_3d_rans) :: status_3d_rans
  ! 3D version increments
  integer, parameter :: size_3di_rans=12*1 ! 2(send/recv)*6(neighbors)*1(variables)
  integer, dimension(size_3di_rans) :: request_3di_rans
  integer, dimension(MPI_STATUS_SIZE,size_3di_rans) :: status_3di_rans

contains

  !!!!!!!!!!!!! Generic communications of unknown variables along faces !!!!!!!!!!

  !===============================================================================
  subroutine communication_2d(var1,var2,var3,var4,var5)
  !===============================================================================
    !> Call non-blocking communications in 2D (faces)
    !> * kept for compatibility *
    ! used in:
    !  * solver.f90 for initial comm of rho,rhou,... & mean0 (compatibility)
    !  * mod_init_TamDong.f90 for comm_U0_TD (restart mode)
    !  * mod_irs_*_v0-2.f90 for increments (compatibility with old versions)
    !  * flux_visc_5pts_SM.f90 for Sij (LES model) TO BE CHANGED
    !  * mod_init_hit.f90 for conservative variables rho,rhou,...
    !  * mod_pp_main.f90 for conservative variables rho,rhou,...
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(inout) :: var1,var2,var3,var4,var5
    ! ----------------------------------------------------------------------------

    call commun2d(var1,1)
    call commun2d(var2,2)
    call commun2d(var3,3)
    call commun2d(var5,4)

    call MPI_WAITALL(size_2d,request_2d,status_2d,info)

  end subroutine communication_2d

  !===============================================================================
  subroutine start_communication_2d(var1,var2,var3,var4,var5)
  !===============================================================================
    !> Start non-blocking communications in 2D (faces)
    !> NOT USED [for comm/computation overlap]
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(inout) :: var1,var2,var3,var4,var5
    ! ----------------------------------------------------------------------------

    call commun2d(var1,1)
    call commun2d(var2,2)
    call commun2d(var3,3)
    call commun2d(var5,4)

  end subroutine start_communication_2d

  !===============================================================================
  subroutine end_communication_2d
  !===============================================================================
    !> Close non-blocking communications in 2D (faces)
    !> NOT USED [for comm/computation overlap]
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    call MPI_WAITALL(size_2d,request_2d,status_2d,info)

  end subroutine end_communication_2d

  !===============================================================================
  subroutine communication_3d(var1,var2,var3,var4,var5)
  !===============================================================================
    !> Call non-blocking communications in 3D (faces)
    !> * kept for compatibility *
    ! used in:
    !  * solver.f90 for initial comm of rho,rhou,... & mean0 (compatibility)
    !  * mod_init_TamDong.f90 for comm_U0_TD (restart mode)
    !  * mod_irs_*_v0-2.f90 for increments (compatibility with old versions)
    !  * flux_visc_5pts_SM.f90 for Sij (LES model) TO BE CHANGED
    !  * mod_init_hit.f90 for conservative variables rho,rhou,...
    !  * mod_pp_main.f90 for conservative variables rho,rhou,...
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(inout) :: var1,var2,var3,var4,var5
    ! ----------------------------------------------------------------------------

    call commun3d(var1,1)
    call commun3d(var2,2)
    call commun3d(var3,3)
    call commun3d(var4,4)
    call commun3d(var5,5)

    call MPI_WAITALL(size_3d,request_3d,status_3d,info)

  end subroutine communication_3d

  !===============================================================================
  subroutine start_communication_3d(var1,var2,var3,var4,var5)
  !===============================================================================
    !> Start non-blocking communications in 3D (faces)
    !> NOT USED [for comm/computation overlap]
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(inout) :: var1,var2,var3,var4,var5
    ! ----------------------------------------------------------------------------

    call commun3d(var1,1)
    call commun3d(var2,2)
    call commun3d(var3,3)
    call commun3d(var4,4)
    call commun3d(var5,5)

  end subroutine start_communication_3d

  !===============================================================================
  subroutine end_communication_3d
  !===============================================================================
    !> Close non-blocking communications in 3D (faces)
    !> NOT USED [for comm/computation overlap]
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    call MPI_WAITALL(size_3d,request_3d,status_3d,info)

  end subroutine end_communication_3d

  !!!!!!!!!!!!!!!!! Communications of unknown vector along faces !!!!!!!!!!!!!!!!!

  !===============================================================================
  subroutine communication2d
  !===============================================================================
    !> Call non-blocking communications in 2D (faces)
    !> new version: June 2022, only for rho_n,rhou_n,...
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    call commun2d(rho_n,1)
    call commun2d(rhou_n,2)
    call commun2d(rhov_n,3)
    call commun2d(rhoe_n,4)

    call MPI_WAITALL(size_2d,request_2d,status_2d,info)

  end subroutine communication2d

  !===============================================================================
  subroutine commun2d(var,i)
  !===============================================================================
    !> NSEW communications using non-blocking ISEND/IRECV (2D faces)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,k,n
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(inout) :: var
    ! ----------------------------------------------------------------------------

    ! Counter: 2(send/recv)*6(neighbors)=12*(variable number)
    k=4*dim*(i-1)

    do n=1,2*dim
       call MPI_ISEND(var(iis(n),ijs(n),iks(n)),1,type_face(n), &
                      neighbor(n),tags(n),COMM_global,request_2d(k+2*n-1),info)
       call MPI_IRECV(var(iir(n),ijr(n),ikr(n)),1,type_face(n), &
                      neighbor(n),tagr(n),COMM_global,request_2d(k+2*n),info)
    enddo

  end subroutine commun2d

  !===============================================================================
  subroutine communication3d
  !===============================================================================
    !> Call non-blocking communications in 3D (faces)
    !> new version: June 2022, only for rho_n,rhou_n,...
    !  (written for easy switch with one-sided communications)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    call commun3d(rho_n,1)
    call commun3d(rhou_n,2)
    call commun3d(rhov_n,3)
    call commun3d(rhow_n,4)
    call commun3d(rhoe_n,5)

    call MPI_WAITALL(size_3d,request_3d,status_3d,info)

  end subroutine communication3d

  !===============================================================================
  subroutine commun3d(var,i)
  !===============================================================================
    !> Communications using non-blocking ISEND/IRECV (3D faces)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,k,n
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(inout) :: var
    ! ----------------------------------------------------------------------------

    ! Counter: 2(send/recv)*6(neighbors)=12*(variable number)
    k=4*dim*(i-1)

    do n=1,2*dim
       call MPI_ISEND(var(iis(n),ijs(n),iks(n)),1,type_face(n), &
                      neighbor(n),tags(n),COMM_global,request_3d(k+2*n-1),info)
       call MPI_IRECV(var(iir(n),ijr(n),ikr(n)),1,type_face(n), &
                      neighbor(n),tagr(n),COMM_global,request_3d(k+2*n),info)
    enddo

  end subroutine commun3d


  !!!!!!!!! Communications of velocity/temperature derivatives along faces !!!!!!!

  !===============================================================================
  subroutine communication_2dv
  !===============================================================================
    !> Call non-blocking communications for velocity/temperature derivatives
    !> in 2D (faces)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    call commun_2dv(dTx,1)
    call commun_2dv(dTy,2)
    call commun_2dv(dux,3)
    call commun_2dv(duy,4)
    call commun_2dv(dvx,5)
    call commun_2dv(dvy,6)

    call MPI_WAITALL(size_v_2d,request_v_2d,status_v_2d,info)

  end subroutine communication_2dv

  !===============================================================================
  subroutine start_communication_2dv
  !===============================================================================
    !> Start non-blocking communications for velocity/temperature derivatives
    !> in 2D (faces) * NOT USED (overlap comm/calc) *
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    ! Communicate faces for velocity derivatives
    ! ==========================================
    call commun_2dv(dTx,1)
    call commun_2dv(dTy,2)
    call commun_2dv(dux,3)
    call commun_2dv(duy,4)
    call commun_2dv(dvx,5)
    call commun_2dv(dvy,6)

  end subroutine start_communication_2dv

  !===============================================================================
  subroutine end_communication_2dv
  !===============================================================================
    !> End non-blocking communications for velocity/temperature derivatives
    !> in 2D (faces) * NOT USED (overlap comm/calc) *
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    call MPI_WAITALL(size_v_2d,request_v_2d,status_v_2d,info)

  end subroutine end_communication_2dv

  !===============================================================================
  subroutine commun_2dv(var,i)
  !===============================================================================
    !> Communications using non-blocking ISEND/IRECV (2D faces)
    !> for velocity/temperature derivatives
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,k,n
    real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v) :: var
    ! ----------------------------------------------------------------------------

    ! Counter: 2(send/recv)*4(neighbors)=8*(variable number)
    k=4*dim*(i-1)

    do n=1,2*dim
       call MPI_ISEND(var(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                      neighbor(n),tags(n),COMM_global,request_v_2d(k+2*n-1),info)
       call MPI_IRECV(var(iir_v(n),ijr_v(n),ikr_v(n)),1,type_face_v(n), &
                      neighbor(n),tagr(n),COMM_global,request_v_2d(k+2*n),info)
    enddo

  end subroutine commun_2dv

  !===============================================================================
  subroutine communication_3dv
  !===============================================================================
    !> Call non-blocking communications for velocity/temperature derivatives
    !> in 3D (faces)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    ! Communicate faces NSEWFB for velocity derivatives
    ! =================================================
    call commun_3dv(dTx,1)
    call commun_3dv(dTy,2)
    call commun_3dv(dTz,3)
    call commun_3dv(dux,4)
    call commun_3dv(duy,5)
    call commun_3dv(duz,6)
    call commun_3dv(dvx,7)
    call commun_3dv(dvy,8)
    call commun_3dv(dvz,9)
    call commun_3dv(dwx,10)
    call commun_3dv(dwy,11)
    call commun_3dv(dwz,12)

    call MPI_WAITALL(size_v,request_v,status_v,info)

  end subroutine communication_3dv

  !===============================================================================
  subroutine start_communication_3dv
  !===============================================================================
    !> Start non-blocking communications for velocity/temperature derivatives
    !> in 3D (faces) * NOT USED (overlap comm/calc) *
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    ! Communicate faces for velocity derivatives
    ! ==========================================
    call commun_3dv(dTx,1)
    call commun_3dv(dTy,2)
    call commun_3dv(dTz,3)
    call commun_3dv(dux,4)
    call commun_3dv(duy,5)
    call commun_3dv(duz,6)
    call commun_3dv(dvx,7)
    call commun_3dv(dvy,8)
    call commun_3dv(dvz,9)
    call commun_3dv(dwx,10)
    call commun_3dv(dwy,11)
    call commun_3dv(dwz,12)

  end subroutine start_communication_3dv

  !===============================================================================
  subroutine end_communication_3dv
  !===============================================================================
    !> End non-blocking communications for velocity/temperature derivatives
    !> in 3D (faces) * NOT USED (overlap comm/calc) *
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------

    call MPI_WAITALL(size_v,request_v,status_v,info)

  end subroutine end_communication_3dv

  !===============================================================================
  subroutine commun_3dv(var,i)
  !===============================================================================
    !> Communications using non-blocking ISEND/IRECV (3D faces)
    !> for velocity/temperature derivatives
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,k,n
    real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v) :: var
    ! ----------------------------------------------------------------------------

    ! Counter: 2(send/recv)*6(neighbors)=12*(variable number)
    k=4*dim*(i-1)

    do n=1,2*dim
       call MPI_ISEND(var(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                      neighbor(n),tags(n),COMM_global,request_v(k+2*n-1),info)
       call MPI_IRECV(var(iir_v(n),ijr_v(n),ikr_v(n)),1,type_face_v(n), &
                      neighbor(n),tagr(n),COMM_global,request_v(k+2*n),info)
    enddo
      
  end subroutine commun_3dv


  !!!!!!!!!!!!!!!!!!! Communications of increments along faces !!!!!!!!!!!!!!!!!!!

  !===============================================================================
  subroutine communication_inc2d(var1,var2,var3,var4,var5,var6)
  !===============================================================================
    !> Call non-blocking communications for 2D increments (faces)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), dimension(nx1_irs:nx2_irs,ny1_irs:ny2_irs,nz1_irs:nz2_irs)&
            , intent(inout) :: var1,var2,var3,var4,var5,var6
    ! ----------------------------------------------------------------------------

    call commun_inc2d(var1,1)
    call commun_inc2d(var2,2)
    call commun_inc2d(var3,3)
    call commun_inc2d(var5,4)
    call commun_inc2d(var6,5)

    call MPI_WAITALL(size_2di,request_2di,status_2di,info)

  end subroutine communication_inc2d

  !===============================================================================
  subroutine commun_inc2d(var,i)
  !===============================================================================
    !> Communications using non-blocking ISEND/IRECV (2D increments)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,k,n
    real(wp), dimension(nx1_irs:nx2_irs,ny1_irs:ny2_irs,nz1_irs:nz2_irs), intent(inout) :: var
    ! ----------------------------------------------------------------------------

    ! Counter: 2(send/recv)*6(neighbors)=12*(variable number)
    k=8*(i-1)

    do n=1,4
       call MPI_ISEND(var(iis_i(n),ijs_i(n),iks_i(n)),1,type_face_i(n), &
                      neighbor(n),tags(n),COMM_global,request_2di(k+2*n-1),info)
       call MPI_IRECV(var(iir_i(n),ijr_i(n),ikr_i(n)),1,type_face_i(n), &
                      neighbor(n),tagr(n),COMM_global,request_2di(k+2*n),info)
    enddo

  end subroutine commun_inc2d

  !===============================================================================
  subroutine communication_inc3d(var1,var2,var3,var4,var5,var6)
  !===============================================================================
    !> Call non-blocking communications for 3D increments (faces)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), dimension(nx1_irs:nx2_irs,ny1_irs:ny2_irs,nz1_irs:nz2_irs)&
            , intent(inout) :: var1,var2,var3,var4,var5,var6
    ! ----------------------------------------------------------------------------

    call commun_inc3d(var1,1)
    call commun_inc3d(var2,2)
    call commun_inc3d(var3,3)
    call commun_inc3d(var4,4)
    call commun_inc3d(var5,5)
    call commun_inc3d(var6,6)

    call MPI_WAITALL(size_3di,request_3di,status_3di,info)

  end subroutine communication_inc3d

  !===============================================================================
  subroutine commun_inc3d(var,i)
  !===============================================================================
    !> Communications using non-blocking ISEND/IRECV (3D increments)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,k,n
    real(wp), dimension(nx1_irs:nx2_irs,ny1_irs:ny2_irs,nz1_irs:nz2_irs), intent(inout) :: var
    ! ----------------------------------------------------------------------------
    
    ! Counter: 2(send/recv)*6(neighbors)=12*(variable number)
    k=12*(i-1)
    
    do n=1,6
       call MPI_ISEND(var(iis_i(n),ijs_i(n),iks_i(n)),1,type_face_i(n), &
                      neighbor(n),tags(n),COMM_global,request_3di(k+2*n-1),info)
       call MPI_IRECV(var(iir_i(n),ijr_i(n),ikr_i(n)),1,type_face_i(n), &
                      neighbor(n),tagr(n),COMM_global,request_3di(k+2*n),info)
    enddo

  end subroutine commun_inc3d

  !===============================================================================
  subroutine communication_inc_k(var1,var2,var3,var4,var5,var6)
  !===============================================================================
    !> Call non-blocking communications for 3D increments in k-direction
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), dimension(nx1_irs:nx2_irs,ny1_irs:ny2_irs,nz1_irs:nz2_irs)&
            , intent(inout) :: var1,var2,var3,var4,var5,var6
    ! ----------------------------------------------------------------------------

    call commun_inc_k(var1,1)
    call commun_inc_k(var2,2)
    call commun_inc_k(var3,3)
    call commun_inc_k(var4,4)
    call commun_inc_k(var5,5)
    call commun_inc_k(var6,6)

    call MPI_WAITALL(size_ik,request_ik,status_ik,info)

  end subroutine communication_inc_k

  !===============================================================================
  subroutine commun_inc_k(var,i)
  !===============================================================================
    !> Communications using non-blocking ISEND/IRECV (3D increments in k-dir)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,k,n
    real(wp), dimension(nx1_irs:nx2_irs,ny1_irs:ny2_irs,nz1_irs:nz2_irs), intent(inout) :: var
    ! ----------------------------------------------------------------------------

    ! Counter: 2(send/recv)*6(neighbors)=12*(variable number)
    k=4*(i-1)

    do n=5,6
       call MPI_ISEND(var(iis_i(n),ijs_i(n),iks_i(n)),1,type_face_i(n), &
                      neighbor(n),tags(n),COMM_global,request_ik(k+2*n-1),info)
       call MPI_IRECV(var(iir_i(n),ijr_i(n),ikr_i(n)),1,type_face_i(n), &
                      neighbor(n),tagr(n),COMM_global,request_ik(k+2*n),info)
    enddo

  end subroutine commun_inc_k

  !!!!!!!!!!!!!!!!!! Communications of RANS variables along faces !!!!!!!!!!!!!!!!

  !===============================================================================
  subroutine communic2d(var)
  !===============================================================================
    !> Communications using non-blocking ISEND/IRECV (2D faces)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: n
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(inout) :: var
    ! ----------------------------------------------------------------------------

    do n=1,2*dim
       call MPI_ISEND(var(iis(n),ijs(n),iks(n)),1,type_face(n), &
                      neighbor(n),tags(n),COMM_global,request2d(2*n-1),info)
       call MPI_IRECV(var(iir(n),ijr(n),ikr(n)),1,type_face(n), &
                      neighbor(n),tagr(n),COMM_global,request2d(2*n),info)
    enddo

    call MPI_WAITALL(size2d,request2d,status2d,info)

  end subroutine communic2d

  !===============================================================================
  subroutine communication_2d_grad_rans(dvarx,dvary,dvarz)
  !===============================================================================
    !> Call non-blocking communications in 2D (faces)
    !> Modif CM 28/04/23: communication_2d_grad_rans takes 3 arguments but only uses
    !  two, because pointer "communication_grad_rans" is the same as for 3D, and
    !  requires an abstract interface. So had to keep the nb of arguments the same.
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v), intent(in) :: dvarx,dvary,dvarz
    ! ----------------------------------------------------------------------------

    ! Communicate faces NSEW for scalar derivatives
    ! =============================================
    call commun_2d_rans(dvarx,1)
    call commun_2d_rans(dvary,2)

    call MPI_WAITALL(size_2d_rans,request_2d_rans,status_2d_rans,info)

  end subroutine communication_2d_grad_rans

  !===============================================================================
  subroutine commun_2d_rans(var,i)
  !===============================================================================
    !> NSEW communications using non-blocking ISEND/IRECV (2D faces)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,k,n
    real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v) :: var
    ! ----------------------------------------------------------------------------

    ! Counter: 2(send/recv)*4(neighbors)=8*(variable number)
    k=4*dim*(i-1)

    do n=1,2*dim
       call MPI_ISEND(var(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                      neighbor(n),tags(n),COMM_global,request_2d_rans(k+2*n-1),info)
       call MPI_IRECV(var(iir_v(n),ijr_v(n),ikr_v(n)),1,type_face_v(n), &
                      neighbor(n),tagr(n),COMM_global,request_2d_rans(k+2*n),info)
    enddo

  end subroutine commun_2d_rans

  !===============================================================================
  subroutine communication_inc2d_rans(var)
  !===============================================================================
    !> Call non-blocking communications in 2D/3D (faces)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), dimension(nx1_irs:nx2_irs,ny1_irs:ny2_irs,nz1_irs:nz2_irs)&
            , intent(inout) :: var
    ! ----------------------------------------------------------------------------

    call commun_inc2d_rans(var,1)

    call MPI_WAITALL(size_2di_rans,request_2di_rans,status_2di_rans,info)

  end subroutine communication_inc2d_rans

  !===============================================================================
  subroutine commun_inc2d_rans(var,i)
  !===============================================================================
    !> NSEW communications using non-blocking ISEND/IRECV (2D faces)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,k,n
    real(wp), dimension(nx1_irs:nx2_irs,ny1_irs:ny2_irs,nz1_irs:nz2_irs), intent(inout) :: var
    ! ----------------------------------------------------------------------------

    ! Counter: 2(send/recv)*6(neighbors)=12*(variable number)
    k=8*(i-1)

    do n=1,4
       call MPI_ISEND(var(iis_i(n),ijs_i(n),iks_i(n)),1,type_face_i(n), &
                      neighbor(n),tags(n),COMM_global,request_2di_rans(k+2*n-1),info)
       call MPI_IRECV(var(iir_i(n),ijr_i(n),ikr_i(n)),1,type_face_i(n), &
                      neighbor(n),tagr(n),COMM_global,request_2di_rans(k+2*n),info)
    enddo

  end subroutine commun_inc2d_rans

  !===============================================================================
  subroutine communic3d(var)
  !===============================================================================
    !> Communications using non-blocking ISEND/IRECV (2D faces)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: n
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(inout) :: var
    ! ----------------------------------------------------------------------------

    do n=1,2*dim
       call MPI_ISEND(var(iis(n),ijs(n),iks(n)),1,type_face(n), &
                      neighbor(n),tags(n),COMM_global,request3d(2*n-1),info)
       call MPI_IRECV(var(iir(n),ijr(n),ikr(n)),1,type_face(n), &
                      neighbor(n),tagr(n),COMM_global,request3d(2*n),info)
    enddo

    call MPI_WAITALL(size3d,request3d,status3d,info)

  end subroutine communic3d

  !===============================================================================
  subroutine communication_3d_grad_rans(dvarx,dvary,dvarz)
  !===============================================================================
    !> Call non-blocking communications in 3D (faces)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v), intent(in) :: dvarx,dvary,dvarz
    ! ----------------------------------------------------------------------------

    ! Communicate faces NSEW for scalar derivatives
    ! =============================================
    call commun_3d_rans(dvarx,1)
    call commun_3d_rans(dvary,2)
    call commun_3d_rans(dvarz,3)

    call MPI_WAITALL(size_3d_rans,request_3d_rans,status_3d_rans,info)

  end subroutine communication_3d_grad_rans

  !===============================================================================
  subroutine commun_3d_rans(var,i)
  !===============================================================================
    !> NSEWFB communications using non-blocking ISEND/IRECV (faces + edges)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,k,n
    real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v) :: var
    ! ----------------------------------------------------------------------------

    ! Counter: 2(send/recv)*6(neighbors)=12*(variable number)
    k=4*dim*(i-1)

    do n=1,2*dim
       call MPI_ISEND(var(iis_v(n),ijs_v(n),iks_v(n)),1,type_face_v(n), &
                      neighbor(n),tags(n),COMM_global,request_3d_rans(k+2*n-1),info)
       call MPI_IRECV(var(iir_v(n),ijr_v(n),ikr_v(n)),1,type_face_v(n), &
                      neighbor(n),tagr(n),COMM_global,request_3d_rans(k+2*n),info)
    enddo

  end subroutine commun_3d_rans

  !===============================================================================
  subroutine communication_inc3d_rans(var)
  !===============================================================================
    !> Call non-blocking communications in 3D (faces)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    real(wp), dimension(nx1_irs:nx2_irs,ny1_irs:ny2_irs,nz1_irs:nz2_irs)&
            , intent(inout) :: var
    ! ----------------------------------------------------------------------------

    call commun_inc3d_rans(var,1)

    call MPI_WAITALL(size_3di_rans,request_3di_rans,status_3di_rans,info)

  end subroutine communication_inc3d_rans

  !===============================================================================
  subroutine commun_inc3d_rans(var,i)
  !===============================================================================
    !> NSEW communications using non-blocking ISEND/IRECV (2D faces)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,k,n
    real(wp), dimension(nx1_irs:nx2_irs,ny1_irs:ny2_irs,nz1_irs:nz2_irs), intent(inout) :: var
    ! ----------------------------------------------------------------------------

    ! Counter: 2(send/recv)*6(neighbors)=12*(variable number)
    k=12*(i-1)

    do n=1,6
       call MPI_ISEND(var(iis_i(n),ijs_i(n),iks_i(n)),1,type_face_i(n), &
                      neighbor(n),tags(n),COMM_global,request_3di_rans(k+2*n-1),info)
       call MPI_IRECV(var(iir_i(n),ijr_i(n),ikr_i(n)),1,type_face_i(n), &
                      neighbor(n),tagr(n),COMM_global,request_3di_rans(k+2*n),info)
    enddo

  end subroutine commun_inc3d_rans


  !===============================================================================
  subroutine commun2d_ex(var)
  !===============================================================================
    !> NSEW communications using non-blocking ISEND/IRECV (2D faces)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: n
    real(wp), dimension(0:ngx+1,0:ngy+1), intent(inout) :: var
    ! ----------------------------------------------------------------------------

    ! Counter: 2(send/recv)*4(neighbors)=8*(variable number)

    do n=1,4
       ! call MPI_ISEND(var(iis_e(n),ijs_e(n),iks_e(n)),1,type_face_e2d(n), &
       !                neighbor(n),tags(n),COMM_global,request2d(2*n-1),info)
       ! call MPI_IRECV(var(iir_e(n),ijr_e(n),ikr_e(n)),1,type_face_e2d(n), &
       !                neighbor(n),tagr(n),COMM_global,request2d(2*n),info)
       call MPI_ISEND(var(iis_e(n),ijs_e(n)),1,type_face_e2d(n), &
                      neighbor(n),tags(n),COMM_global,request2d(2*n-1),info)
       call MPI_IRECV(var(iir_e(n),ijr_e(n)),1,type_face_e2d(n), &
                      neighbor(n),tagr(n),COMM_global,request2d(2*n),info)
    enddo

    call MPI_WAITALL(size2d,request2d,status2d,info)

  end subroutine commun2d_ex

end module mod_comm
  
