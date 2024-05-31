!==============================================================================
module mod_pp_transpose
!==============================================================================
  !> Module for array transpositions (used for FFT parallelization)
!==============================================================================
  use mod_constant
  use mod_flow
  use mod_mpi_part ! not in mod_pp_mpi ???
  use mod_pp_mpi
  use mod_pp_var
  use mod_io_snapshots
  implicit none
  ! ---------------------------------------------------------------------------
  integer, private :: ip,ip2,sizeofcplx
  integer(kind=MPI_ADDRESS_KIND), private :: pas
  integer, private :: type_base1,type_base2,type_tr
  ! ---------------------------------------------------------------------------

contains

  !============================================================================
  subroutine transpose_zt_r
  !============================================================================
    !> Transposition z <-> t - for arrays of real numbers
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: k,n
    ! -------------------------------------------------------------------------

    ! Allocate transposed array
    ! =========================
    allocate(var_rt(nl1,nt_tzt,ngz))

    ! Definition of MPI types for communications in transposition step
    ! ==========================================
    call MPI_TYPE_VECTOR(nl1,1,1,MPI_DOUBLE_PRECISION,type_base1,info)
    pas=nl1*sizeofreal
    call MPI_TYPE_CREATE_HVECTOR(nt_tzt,1,pas,type_base1,type_base2,info)
    pas=nl1*nt_tzt*sizeofreal
    call MPI_TYPE_CREATE_HVECTOR(nz,1,pas,type_base2,type_tr,info)
    call MPI_TYPE_COMMIT(type_tr,info)
    call MPI_BARRIER(COMM_global,info)

    ! Local transposition
    ! ===================
    do ip=0,ndomz-1
       do k=1,nz
          do n=ip*nt_tzt+1,(ip+1)*nt_tzt
             var_rt(:,n-ip*nt_tzt,k+ip*nz)=var_r(:,k,n)
          enddo
       enddo
    enddo

    ! Global transposition (MPI communications)
    ! ====================
    do ip=0,nproc-1
       if (coordx1(ip)==coordx1(iproc)) then
          if (snapshots(nsr)%type.eq.1) then
             ip2=ip
          else
             ip2=ip-coordx1(ip)/nl1*ndomz
          endif
          call MPI_SENDRECV_REPLACE(var_rt(1,1,ip2*nz+1),1,type_tr,ip,tag,ip,tag, &
                                    COMM_global,status,info)
       endif
    enddo
    call MPI_BARRIER(COMM_global,info)

  end subroutine transpose_zt_r

  !============================================================================
  subroutine untranspose_zt_r
  !============================================================================
    !> Reverse transposition z <-> t - for arrays of real numbers
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: k,n
    ! -------------------------------------------------------------------------

    ! Definition of MPI types for communications in transposition step
    ! ==========================================
    call MPI_TYPE_VECTOR(nl1,1,1,MPI_DOUBLE_PRECISION,type_base1,info)
    pas=nl1*sizeofreal
    call MPI_TYPE_CREATE_HVECTOR(nz,1,pas,type_base1,type_base2,info)
    pas=nl1*nz*sizeofreal
    call MPI_TYPE_CREATE_HVECTOR(nt_tzt,1,pas,type_base2,type_tr,info)
    call MPI_TYPE_COMMIT(type_tr,info)
    call MPI_BARRIER(COMM_global,info)

    ! Local transposition
    ! ===================
    do ip=0,ndomz-1
       do n=1,nt_tzt
          do k=ip*nz+1,(ip+1)*nz
             var_r(:,k-ip*nz,n+ip*nt_tzt)=var_rt(:,n,k)
          enddo
       enddo
    enddo

    ! Global transposition (MPI communications)
    ! ====================
    do ip=0,nproc-1
       if (coordx1(ip)==coordx1(iproc)) then
          if (snapshots(nsr)%type.eq.1) then
             ip2=ip
          else
             ip2=ip-coordx1(ip)/nl1*ndomz
          endif
          call MPI_SENDRECV_REPLACE(var_r(1,1,ip2*nt_tzt+1),1,type_tr,ip,tag,ip,tag, &
                                    COMM_global,status,info)
       endif
    enddo
    call MPI_BARRIER(COMM_global,info)

    ! Free transposed array
    ! =====================
    deallocate(var_rt)

  end subroutine untranspose_zt_r

  !============================================================================
  subroutine transpose_zt_c
  !============================================================================
    !> Transposition z <-> t - for arrays of complex numbers
  !============================================================================
    use warnstop
    implicit none
    ! -------------------------------------------------------------------------
    integer :: k,n
    ! -------------------------------------------------------------------------

    ! if (iproc.eq.0) then

    ! Allocate transposed array
    ! =========================
    allocate(var_ct(nl1,nt_tzt,ngz))

    ! Definition of MPI types for communications in transposition step
    ! ==========================================
    call MPI_TYPE_EXTENT(MPI_DOUBLE_COMPLEX,sizeofcplx,info)
    call MPI_TYPE_VECTOR(nl1,1,1,MPI_DOUBLE_COMPLEX,type_base1,info)
    pas=nl1*sizeofcplx
    call MPI_TYPE_CREATE_HVECTOR(nt_tzt,1,pas,type_base1,type_base2,info)
    pas=nl1*nt_tzt*sizeofcplx
    call MPI_TYPE_CREATE_HVECTOR(nz,1,pas,type_base2,type_tr,info)
    call MPI_TYPE_COMMIT(type_tr,info)
    call MPI_BARRIER(COMM_global,info)

    ! print *,"nt_tzt",nt_tzt
    ! print *,"nz",nz
    ! ! print *,"var_c(1,1,n)",var_c(1,1,:)

    ! Local transposition
    ! ===================
    do ip=0,ndomz-1
       do k=1,nz
          do n=ip*nt_tzt+1,(ip+1)*nt_tzt
             var_ct(:,n-ip*nt_tzt,k+ip*nz)=var_c(:,k,n)
          enddo
       enddo
    enddo

    ! ip=0
    ! print *,"ip",ip
    ! print *,"ip*nt_tzt+1,(ip+1)*nt_tzt",ip*nt_tzt+1,(ip+1)*nt_tzt
    ! ! print *,"var_ct(1,:)",var_ct(1,:)


    ! Global transposition (MPI communications)
    ! ====================
    do ip=0,nproc-1
       if (coordx1(ip)==coordx1(iproc)) then
          if (snapshots(nsr)%type.eq.1) then
             ip2=ip
          else
             ip2=ip-coordx1(ip)/nl1*ndomz
          endif
          call MPI_SENDRECV_REPLACE(var_ct(1,1,ip2*nz+1),1,type_tr,ip,tag,ip,tag, &
                                    COMM_global,status,info)
       endif
    enddo
    call MPI_BARRIER(COMM_global,info)

    ! endif

    ! call mpistop('',0)

  end subroutine transpose_zt_c

  !============================================================================
  subroutine untranspose_zt_c
  !============================================================================
    !> Reverse transposition z <-> t - for arrays of complex numbers
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: k,n
    ! -------------------------------------------------------------------------

    ! Definition of MPI types for communications in transposition step
    ! ==========================================
    call MPI_TYPE_EXTENT(MPI_DOUBLE_COMPLEX,sizeofcplx,info)
    call MPI_TYPE_VECTOR(nl1,1,1,MPI_DOUBLE_COMPLEX,type_base1,info)
    pas=nl1*sizeofcplx
    call MPI_TYPE_CREATE_HVECTOR(nz,1,pas,type_base1,type_base2,info)
    pas=nl1*nz*sizeofcplx
    call MPI_TYPE_CREATE_HVECTOR(nt_tzt,1,pas,type_base2,type_tr,info)
    call MPI_TYPE_COMMIT(type_tr,info)
    call MPI_BARRIER(COMM_global,info)

    ! Local transposition
    ! ===================
    do ip=0,ndomz-1
       do n=1,nt_tzt
          do k=ip*nz+1,(ip+1)*nz
             var_c(:,k-ip*nz,n+ip*nt_tzt)=var_ct(:,n,k)
          enddo
       enddo
    enddo

    ! Global transposition (MPI communications)
    ! ====================
    do ip=0,nproc-1
       if (coordx1(ip)==coordx1(iproc)) then
          if (snapshots(nsr)%type.eq.1) then
             ip2=ip
          else
             ip2=ip-coordx1(ip)/nl1*ndomz
          endif
          call MPI_SENDRECV_REPLACE(var_c(1,1,ip2*nt_tzt+1),1,type_tr,ip,tag,ip,tag, &
                                    COMM_global,status,info)
       endif
    enddo
    call MPI_BARRIER(COMM_global,info)

    ! Free transposed array
    ! =====================
    deallocate(var_ct)

  end subroutine untranspose_zt_c

  !============================================================================
  subroutine transpose_xt_r
  !============================================================================
    !> Transposition x <-> t - for arrays of real numbers
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i,n
    ! -------------------------------------------------------------------------

    ! Allocate transposed array
    ! =========================
    allocate(var_rt(nt_txt,nl2,ngx))

    ! Definition of MPI types for communications in transposition step
    ! ==========================================
    call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION,sizeofreal,info)
    call MPI_TYPE_VECTOR(nt_txt,1,1,MPI_DOUBLE_PRECISION,type_base1,info)
    pas=nt_txt*sizeofreal
    call MPI_TYPE_CREATE_HVECTOR(nl2,1,pas,type_base1,type_base2,info)
    pas=nt_txt*nl2*sizeofreal
    call MPI_TYPE_CREATE_HVECTOR(nx,1,pas,type_base2,type_tr,info)
    call MPI_TYPE_COMMIT(type_tr,info)
    call MPI_BARRIER(COMM_global,info)

    ! Local transposition
    ! ===================
    do ip=0,ndomx-1
       do i=1,nx
          do n=ip*nt_txt+1,(ip+1)*nt_txt
             var_rt(n-ip*nt_txt,:,i+ip*nx)=var_r(i,:,n)
          enddo
       enddo
    enddo

    ! Global transposition (MPI communications)
    ! ====================
    do ip=0,nproc-1
       if (coordx2(ip)==coordx2(iproc)) then
          if (snapshots(nsr)%type.eq.1) then
             ip2=0
          else
             ip2=coordx2(ip)/nl2
          endif
          ! ip2=coordx2(ip)/nl2
          call MPI_SENDRECV_REPLACE(var_rt(1,1,ip2*nx+1),1,type_tr,ip,tag,ip,tag, &
                                    COMM_global,status,info)
       endif
    enddo
    call MPI_BARRIER(COMM_global,info)

  end subroutine transpose_xt_r

  !============================================================================
  subroutine untranspose_xt_r
  !============================================================================
    !> Reverse transposition x <-> t - for arrays of real numbers
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i,n
    ! -------------------------------------------------------------------------

    ! Definition of MPI types for communications in transposition step
    ! ==========================================
    call MPI_TYPE_VECTOR(nx,1,1,MPI_DOUBLE_PRECISION,type_base1,info)
    pas=nx*sizeofreal
    call MPI_TYPE_CREATE_HVECTOR(nl2,1,pas,type_base1,type_base2,info)
    pas=nx*nl2*sizeofreal
    call MPI_TYPE_CREATE_HVECTOR(nt_txt,1,pas,type_base2,type_tr,info)
    call MPI_TYPE_COMMIT(type_tr,info)
    call MPI_BARRIER(COMM_global,info)

    ! Local transposition
    ! ===================
    do ip=0,ndomx-1
       do n=1,nt_txt
          do i=ip*nx+1,(ip+1)*nx
             var_r(i-ip*nx,:,n+ip*nt_txt)=var_rt(n,:,i)
          enddo
       enddo
    enddo

    ! Global transposition (MPI communications)
    ! ====================
    do ip=0,nproc-1
       if (coordx2(ip)==coordx2(iproc)) then
          if (snapshots(nsr)%type.eq.1) then
             ip2=0
          else
             ip2=coordx2(ip)/nl2
          endif
          ! ip2=coordx2(ip)/nl2
          call MPI_SENDRECV_REPLACE(var_r(1,1,ip2*nt_txt+1),1,type_tr,ip,tag,ip,tag, &
                                    COMM_global,status,info)
       endif
    enddo
    call MPI_BARRIER(COMM_global,info)

    ! Free transposed array
    ! =====================
    deallocate(var_rt)

  end subroutine untranspose_xt_r

  !============================================================================
  subroutine transpose_xt_c
  !============================================================================
    !> Transposition x <-> t - for arrays of complex numbers
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i,n
    ! -------------------------------------------------------------------------

    ! Allocate transposed array
    ! =========================
    allocate(var_ct(nt_txt,nl2,ngx))

    ! Definition of MPI types for communications in transposition step
    ! ==========================================
    call MPI_TYPE_EXTENT(MPI_DOUBLE_COMPLEX,sizeofcplx,info)
    call MPI_TYPE_VECTOR(nt_txt,1,1,MPI_DOUBLE_COMPLEX,type_base1,info)
    pas=nt_txt*sizeofcplx
    call MPI_TYPE_CREATE_HVECTOR(nl2,1,pas,type_base1,type_base2,info)
    pas=nt_txt*nl2*sizeofcplx
    call MPI_TYPE_CREATE_HVECTOR(nx,1,pas,type_base2,type_tr,info)
    call MPI_TYPE_COMMIT(type_tr,info)
    call MPI_BARRIER(COMM_global,info)

    ! Local transposition
    ! ===================
    do ip=0,ndomx-1
       do i=1,nx
          do n=ip*nt_txt+1,(ip+1)*nt_txt
             var_ct(n-ip*nt_txt,:,i+ip*nx)=var_c(i,:,n)
          enddo
       enddo
    enddo

    ! Global transposition (MPI communications)
    ! ====================
    do ip=0,nproc-1
       if (coordx2(ip)==coordx2(iproc)) then
          if (snapshots(nsr)%type.eq.1) then
             ip2=0
          else
             ip2=coordx2(ip)/nl2
          endif
          ! ip2=coordx2(ip)/nl2
          call MPI_SENDRECV_REPLACE(var_ct(1,1,ip2*nx+1),1,type_tr,ip,tag,ip,tag, &
                                    COMM_global,status,info)
       endif
    enddo
    call MPI_BARRIER(COMM_global,info)

  end subroutine transpose_xt_c

  !============================================================================
  subroutine untranspose_xt_c
  !============================================================================
    !> Reverse transposition x <-> t - for arrays of complex numbers
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i,n
    ! -------------------------------------------------------------------------

    ! Definition of MPI types for communications in transposition step
    ! ==========================================
    call MPI_TYPE_EXTENT(MPI_DOUBLE_COMPLEX,sizeofcplx,info)
    call MPI_TYPE_VECTOR(nx,1,1,MPI_DOUBLE_COMPLEX,type_base1,info)
    pas=nx*sizeofcplx
    call MPI_TYPE_CREATE_HVECTOR(nl2,1,pas,type_base1,type_base2,info)
    pas=nx*nl2*sizeofcplx
    call MPI_TYPE_CREATE_HVECTOR(nt_txt,1,pas,type_base2,type_tr,info)
    call MPI_TYPE_COMMIT(type_tr,info)
    call MPI_BARRIER(COMM_global,info)

    ! Local transposition
    ! ===================
    do ip=0,ndomx-1
       do n=1,nt_txt
          do i=ip*nx+1,(ip+1)*nx
             var_c(i-ip*nx,:,n+ip*nt_txt)=var_ct(n,:,i)
          enddo
       enddo
    enddo

    ! Global transposition (MPI communications)
    ! ====================
    do ip=0,nproc-1
       if (coordx2(ip)==coordx2(iproc)) then
          if (snapshots(nsr)%type.eq.1) then
             ip2=0
          else
             ip2=coordx2(ip)/nl2
          endif
          ! ip2=coordx2(ip)/nl2
          call MPI_SENDRECV_REPLACE(var_c(1,1,ip2*nt_txt+1),1,type_tr,ip,tag,ip,tag, &
                                    COMM_global,status,info)
       endif
    enddo
    call MPI_BARRIER(COMM_global,info)

    ! Free transposed array
    ! =====================
    deallocate(var_ct)

  end subroutine untranspose_xt_c

end module mod_pp_transpose
