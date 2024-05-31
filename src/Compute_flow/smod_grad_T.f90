!=================================================================================
submodule (mod_gradient) smod_grad_T
!=================================================================================
  !> author: XG
  !> date: January 2024
  !> Cartesian version - routines to compute gradients of T 
!=================================================================================

contains

  !===============================================================================
  module subroutine grad_T_5pts
  !===============================================================================
    !> Compute temperature derivatives
    !> 5-point stencil - Cartesian version -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k
    ! ----------------------------------------------------------------------------

    ! Derivatives of T along x
    ! ========================
    if (is_boundary(1,1)) then
       i=1
       do k=1,nz
          do j=1,ny
             dTx(i,j,k)= ( a04(1)*Tmp(i  ,j,k)+a04(2)*Tmp(i+1,j,k) &
                          +a04(3)*Tmp(i+2,j,k)+a04(4)*Tmp(i+3,j,k) &
                          +a04(5)*Tmp(i+4,j,k) )*idx_v(i)
          enddo
       enddo
       ! enforce adiabaticity
       if ((BC_face(1,1)%sort==0).and.(is_adiab)) then
          do k=1,nz
             do j=1,ny
                dTx(1,j,k)=0.0_wp
             enddo
          enddo
       endif

       i=2
       do k=1,nz
          do j=1,ny
             dTx(i,j,k)= ( a13(1)*Tmp(i-1,j,k)+a13(2)*Tmp(i  ,j,k) &
                          +a13(3)*Tmp(i+1,j,k)+a13(4)*Tmp(i+2,j,k) &
                          +a13(5)*Tmp(i+3,j,k) )*idx_v(i)
          enddo
       enddo
    endif

    do k=1,nz
       do j=1,ny
          do i=ndx_v3,nfx_v3
             dTx(i,j,k)= ( a5(1)*(Tmp(i+1,j,k)-Tmp(i-1,j,k)) &
                         + a5(2)*(Tmp(i+2,j,k)-Tmp(i-2,j,k)) )*idx_v(i)
          enddo
       enddo
    enddo

    if (is_boundary(1,2)) then
       i=nx-1
       do k=1,nz
          do j=1,ny
             dTx(i,j,k)= ( a31(5)*Tmp(i-3,j,k)+a31(4)*Tmp(i-2,j,k) &
                          +a31(3)*Tmp(i-1,j,k)+a31(2)*Tmp(i  ,j,k) &
                          +a31(1)*Tmp(i+1,j,k) )*idx_v(i)
          enddo
       enddo

       i=nx
       do k=1,nz
          do j=1,ny
             dTx(i,j,k)= ( a40(5)*Tmp(i-4,j,k)+a40(4)*Tmp(i-3,j,k) &
                          +a40(3)*Tmp(i-2,j,k)+a40(2)*Tmp(i-1,j,k) &
                          +a40(1)*Tmp(i  ,j,k) )*idx_v(i)
          enddo
       enddo

       ! enforce adiabaticity
       if ((BC_face(1,2)%sort==0).and.(is_adiab)) then       
          do k=1,nz
             do j=1,ny
                dTx(nx,j,k)=0.0_wp
             enddo
          enddo
       endif
    endif

    ! Derivatives of T along y
    ! ========================
    if (is_boundary(2,1)) then
       j=1
       do k=1,nz
          do i=1,nx
             dTy(i,j,k)= ( a04(1)*Tmp(i,j  ,k)+a04(2)*Tmp(i,j+1,k) &
                          +a04(3)*Tmp(i,j+2,k)+a04(4)*Tmp(i,j+3,k) &
                          +a04(5)*Tmp(i,j+4,k) )*idy_v(j)
          enddo
       enddo

       ! enforce adiabaticity
       if ((BC_face(2,1)%sort==0).and.(is_adiab)) then
          do k=1,nz
             do i=1,nx
                dTy(i,1,k)=0.0_wp
             enddo
          enddo
       endif

       j=2
       do k=1,nz
          do i=1,nx
             dTy(i,j,k)= ( a13(1)*Tmp(i,j-1,k)+a13(2)*Tmp(i,j  ,k) &
                          +a13(3)*Tmp(i,j+1,k)+a13(4)*Tmp(i,j+2,k) &
                          +a13(5)*Tmp(i,j+3,k) )*idy_v(j)
          enddo
       enddo
    endif

    do k=1,nz
       do j=ndy_v3,nfy_v3
          do i=1,nx
             dTy(i,j,k)= ( a5(1)*(Tmp(i,j+1,k)-Tmp(i,j-1,k)) &
                         + a5(2)*(Tmp(i,j+2,k)-Tmp(i,j-2,k)) )*idy_v(j)
          enddo
       enddo
    enddo

    if (is_boundary(2,2)) then
       j=ny-1
       do k=1,nz
          do i=1,nx
             dTy(i,j,k)= ( a31(5)*Tmp(i,j-3,k)+a31(4)*Tmp(i,j-2,k) &
                          +a31(3)*Tmp(i,j-1,k)+a31(2)*Tmp(i,j  ,k) &
                          +a31(1)*Tmp(i,j+1,k) )*idy_v(j)
          enddo
       enddo

       j=ny
       do k=1,nz
          do i=1,nx
             dTy(i,j,k)= ( a40(5)*Tmp(i,j-4,k)+a40(4)*Tmp(i,j-3,k) &
                          +a40(3)*Tmp(i,j-2,k)+a40(2)*Tmp(i,j-1,k) &
                          +a40(1)*Tmp(i,j  ,k) )*idy_v(j)
          enddo
       enddo
       ! enforce adiabaticity
       if ((BC_face(2,2)%sort==0).and.(is_adiab)) then
          do k=1,nz
             do i=1,nx
                dTy(i,ny,k)=0.0_wp
             enddo
          enddo
       endif
    endif

    !****************
    if (is_2D) return
    !****************

    ! Derivatives of T along z
    ! ========================
    if (is_boundary(3,1)) then
       k=1
       do j=1,ny
          do i=1,nx
             dTz(i,j,k)= ( a04(1)*Tmp(i,j,k  )+a04(2)*Tmp(i,j,k+1) &
                          +a04(3)*Tmp(i,j,k+2)+a04(4)*Tmp(i,j,k+3) &
                          +a04(5)*Tmp(i,j,k+4) )*idz_v(k)
          enddo
       enddo
       ! enforce adiabaticity
       if ((BC_face(3,1)%sort==0).and.(is_adiab)) then
          do j=1,ny
             do i=1,nx
                dTz(i,j,1)=0.0_wp
             enddo
          enddo
       endif

       k=2
       do j=1,ny
          do i=1,nx
             dTz(i,j,k)= ( a13(1)*Tmp(i,j,k-1)+a13(2)*Tmp(i,j,k  ) &
                          +a13(3)*Tmp(i,j,k+1)+a13(4)*Tmp(i,j,k+2) &
                          +a13(5)*Tmp(i,j,k+3) )*idz_v(k)
          enddo
       enddo
    endif

    do k=ndz_v3,nfz_v3
       do j=1,ny
          do i=1,nx
             dTz(i,j,k)= ( a5(1)*(Tmp(i,j,k+1)-Tmp(i,j,k-1)) &
                         + a5(2)*(Tmp(i,j,k+2)-Tmp(i,j,k-2)) )*idz_v(k)
          enddo
       enddo
    enddo

    if (is_boundary(3,2)) then
       k=nz-1
       do j=1,ny
          do i=1,nx
             dTz(i,j,k)= ( a31(5)*Tmp(i,j,k-3)+a31(4)*Tmp(i,j,k-2) &
                          +a31(3)*Tmp(i,j,k-1)+a31(2)*Tmp(i,j,k  ) &
                          +a31(1)*Tmp(i,j,k+1) )*idz_v(k)
          enddo
       enddo

       k=nz
       do j=1,ny
          do i=1,nx
             dTz(i,j,k)= ( a40(5)*Tmp(i,j,k-4)+a40(4)*Tmp(i,j,k-3) &
                          +a40(3)*Tmp(i,j,k-2)+a40(2)*Tmp(i,j,k-1) &
                          +a40(1)*Tmp(i,j,k  ) )*idz_v(k)
          enddo
       enddo
       ! enforce adiabaticity
       if ((BC_face(3,2)%sort==0).and.(is_adiab)) then
          do j=1,ny
             do i=1,nx
                dTz(i,j,nz)=0.0_wp
             enddo
          enddo
       endif
    endif

  end subroutine grad_T_5pts

  !===============================================================================
  module subroutine grad_T_3pts
  !===============================================================================
    !> Compute temperature derivatives
    !> 3-point stencil - Cartesian version -
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k
    ! ----------------------------------------------------------------------------

    ! Derivatives of T along x
    ! ========================
    if (is_boundary(1,1)) then
       i=1
       do k=1,nz
          do j=1,ny
             dTx(i,j,k)= ( a02(1)*Tmp(i  ,j,k)+a02(2)*Tmp(i+1,j,k) &
                          +a02(3)*Tmp(i+2,j,k))*idx_v(i)
          enddo
       enddo

       ! enforce adiabaticity
       if ((BC_face(1,1)%sort==0).and.(is_adiab)) then
          do k=1,nz
             do j=1,ny
                dTx(1,j,k)=0.0_wp
             enddo
          enddo
       endif
    endif

    do k=1,nz
       do j=1,ny
          do i=ndx_v3,nfx_v3
             dTx(i,j,k) = 0.5_wp*(Tmp(i+1,j,k)-Tmp(i-1,j,k))*idx_v(i)
          enddo
       enddo
    enddo

    if (is_boundary(1,2)) then
       i=nx
       do k=1,nz
          do j=1,ny
             dTx(i,j,k)= ( a20(3)*Tmp(i-2,j,k)+a20(2)*Tmp(i-1,j,k) &
                          +a20(1)*Tmp(i  ,j,k) )*idx_v(i)
          enddo
       enddo

       ! enforce adiabaticity
       if ((BC_face(1,2)%sort==0).and.(is_adiab)) then       
          do k=1,nz
             do j=1,ny
                dTx(nx,j,k)=0.0_wp
             enddo
          enddo
       endif
    endif

    ! Derivatives of T along y
    ! ========================
    if (is_boundary(2,1)) then
       j=1
       do k=1,nz
          do i=1,nx
             dTy(i,j,k)= ( a02(1)*Tmp(i,j  ,k)+a02(2)*Tmp(i,j+1,k) &
                          +a02(3)*Tmp(i,j+2,k))*idy_v(j)
          enddo
       enddo

       ! enforce adiabaticity
       if ((BC_face(2,1)%sort==0).and.(is_adiab)) then
          do k=1,nz
             do i=1,nx
                dTy(i,1,k) = 0.0_wp
             enddo
          enddo
       endif
    endif

    do k=1,nz
       do j=ndy_v3,nfy_v3
          do i=1,nx
             dTy(i,j,k) = 0.5_wp*(Tmp(i,j+1,k)-Tmp(i,j-1,k))*idy_v(j)
          enddo
       enddo
    enddo

    if (is_boundary(2,2)) then
       j=ny
       do k=1,nz
          do i=1,nx
             dTy(i,j,k)= ( a20(3)*Tmp(i,j-2,k)+a20(2)*Tmp(i,j-1,k) &
                          +a20(1)*Tmp(i,j  ,k) )*idy_v(j)
          enddo
       enddo

       ! enforce adiabaticity
       if ((BC_face(2,2)%sort==0).and.(is_adiab)) then
          do k=1,nz
             do i=1,nx
                dTy(i,ny,k) = 0.0_wp
             enddo
          enddo
       endif
    endif

    !****************
    if (is_2D) return
    !****************

    ! Derivatives of T along z
    ! ========================
    if (is_boundary(3,1)) then
       k=1
       do j=1,ny
          do i=1,nx
             dTz(i,j,k)= ( a02(1)*Tmp(i,j,k  )+a02(2)*Tmp(i,j,k+1) &
                          +a02(3)*Tmp(i,j,k+2))*idz_v(k)
          enddo
       enddo

       ! enforce adiabaticity
       if ((BC_face(3,1)%sort==0).and.(is_adiab)) then
          do j=1,ny
             do i=1,nx
                dTz(i,j,1)=0.0_wp
             enddo
          enddo
       endif
    endif

    do k=ndz_v3,nfz_v3
       do j=1,ny
          do i=1,nx
             dTz(i,j,k) = 0.5_wp*(Tmp(i,j,k+1)-Tmp(i,j,k-1))*idz_v(k)
          enddo
       enddo
    enddo

    if (is_boundary(3,2)) then
       k=nz
       do j=1,ny
          do i=1,nx
             dTz(i,j,k)= ( a20(3)*Tmp(i,j,k-2)+a20(2)*Tmp(i,j,k-1) &
                          +a20(1)*Tmp(i,j,k  ) )*idz_v(k)
          enddo
       enddo

       ! enforce adiabaticity
       if ((BC_face(3,2)%sort==0).and.(is_adiab)) then
          do j=1,ny
             do i=1,nx
                dTz(i,j,nz) = 0.0_wp
             enddo
          enddo
       endif
    endif

  end subroutine grad_T_3pts

end submodule smod_grad_T
