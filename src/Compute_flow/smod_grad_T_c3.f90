!=================================================================================
submodule (mod_gradient) smod_grad_T_c3
!=================================================================================
  !> author: XG
  !> date: January 2024
  !> Full 3D curvilinear version - routines to compute gradients of T 
!=================================================================================

contains

  !===============================================================================
  module subroutine grad_T_5pts_c3
  !===============================================================================
    !> Compute temperature derivatives
    !> 5-point stencil - 3D curvilinear version
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: dvardx,dvardy,dvardz
    real(wp), dimension(ny1_v:ny2_v,nz1_v:nz2_v) :: dTeta_imin,dTeta_imax
    real(wp), dimension(ny1_v:ny2_v,nz1_v:nz2_v) :: dTphi_imin,dTphi_imax
    real(wp), dimension(nx1_v:nx2_v,nz1_v:nz2_v) :: dTksi_jmin,dTksi_jmax
    real(wp), dimension(nx1_v:nx2_v,nz1_v:nz2_v) :: dTphi_jmin,dTphi_jmax
    real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v) :: dTksi_kmin,dTksi_kmax
    real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v) :: dTeta_kmin,dTeta_kmax
    ! ---------------------------------------------------------------------------

    ! Derivatives of T along ksi (Caution: dTx is dT/dksi)
    ! ==========================
    if (is_boundary(1,1)) then
       i=1
       do k=1,nz
          do j=1,ny
             dTx(i,j,k)= a04(1)*Tmp(i  ,j,k)+a04(2)*Tmp(i+1,j,k) &
                        +a04(3)*Tmp(i+2,j,k)+a04(4)*Tmp(i+3,j,k) &
                        +a04(5)*Tmp(i+4,j,k)
          enddo
       enddo

       i=2
       do k=1,nz
          do j=1,ny
             dTx(i,j,k)= a13(1)*Tmp(i-1,j,k)+a13(2)*Tmp(i  ,j,k) &
                        +a13(3)*Tmp(i+1,j,k)+a13(4)*Tmp(i+2,j,k) &
                        +a13(5)*Tmp(i+3,j,k)
          enddo
       enddo
    endif

    do k=1,nz
       do j=1,ny
          do i=ndx_v3,nfx_v3
             dTx(i,j,k) = a5(1)*(Tmp(i+1,j,k)-Tmp(i-1,j,k)) &
                        + a5(2)*(Tmp(i+2,j,k)-Tmp(i-2,j,k))
          enddo
       enddo
    enddo

    if (is_boundary(1,2)) then
       i=nx-1
       do k=1,nz
          do j=1,ny
             dTx(i,j,k)= a31(5)*Tmp(i-3,j,k)+a31(4)*Tmp(i-2,j,k) &
                        +a31(3)*Tmp(i-1,j,k)+a31(2)*Tmp(i  ,j,k) &
                        +a31(1)*Tmp(i+1,j,k)
          enddo
       enddo

       i=nx
       do k=1,nz
          do j=1,ny
             dTx(i,j,k)= a40(5)*Tmp(i-4,j,k)+a40(4)*Tmp(i-3,j,k) &
                        +a40(3)*Tmp(i-2,j,k)+a40(2)*Tmp(i-1,j,k) &
                        +a40(1)*Tmp(i  ,j,k)
          enddo
       enddo
    endif

    ! store ksi-derivatives on boundaries to enforce adiabaticity
    dTksi_jmin=dTx(:,1,:)
    dTksi_jmax=dTx(:,ny,:)
    dTksi_kmin=dTx(:,:,1)
    dTksi_kmax=dTx(:,:,nz)

    ! Derivatives of T along eta (Caution: dTy is dT/deta)
    ! ==========================
    if (is_boundary(2,1)) then
       j=1
       do k=1,nz
          do i=1,nx
             dTy(i,j,k)= a04(1)*Tmp(i,j  ,k)+a04(2)*Tmp(i,j+1,k) &
                        +a04(3)*Tmp(i,j+2,k)+a04(4)*Tmp(i,j+3,k) &
                        +a04(5)*Tmp(i,j+4,k)
          enddo
       enddo

       j=2
       do k=1,nz
          do i=1,nx
             dTy(i,j,k)= a13(1)*Tmp(i,j-1,k)+a13(2)*Tmp(i,j  ,k) &
                        +a13(3)*Tmp(i,j+1,k)+a13(4)*Tmp(i,j+2,k) &
                        +a13(5)*Tmp(i,j+3,k)
          enddo
       enddo
    endif

    do k=1,nz
       do j=ndy_v3,nfy_v3
          do i=1,nx
             dTy(i,j,k) = a5(1)*(Tmp(i,j+1,k)-Tmp(i,j-1,k)) &
                        + a5(2)*(Tmp(i,j+2,k)-Tmp(i,j-2,k))
          enddo
       enddo
    enddo

    if (is_boundary(2,2)) then
       j=ny-1
       do k=1,nz
          do i=1,nx
             dTy(i,j,k)= a31(5)*Tmp(i,j-3,k)+a31(4)*Tmp(i,j-2,k) &
                        +a31(3)*Tmp(i,j-1,k)+a31(2)*Tmp(i,j  ,k) &
                        +a31(1)*Tmp(i,j+1,k)
          enddo
       enddo

       j=ny
       do k=1,nz
          do i=1,nx
             dTy(i,j,k)= a40(5)*Tmp(i,j-4,k)+a40(4)*Tmp(i,j-3,k) &
                        +a40(3)*Tmp(i,j-2,k)+a40(2)*Tmp(i,j-1,k) &
                        +a40(1)*Tmp(i,j  ,k)
          enddo
       enddo
    endif

    ! store eta-derivatives on boundaries to enforce adiabaticity
    dTeta_imin=dTy(1,:,:)
    dTeta_imax=dTy(nx,:,:)
    dTeta_kmin=dTy(:,:,1)
    dTeta_kmax=dTy(:,:,nz)

    ! Derivatives of T along phi (Caution: dTz is dT/deta)
    ! ==========================
    if (is_boundary(3,1)) then
       k=1
       do j=1,ny
          do i=1,nx
             dTz(i,j,k)= a04(1)*Tmp(i,j,k  )+a04(2)*Tmp(i,j,k+1) &
                        +a04(3)*Tmp(i,j,k+2)+a04(4)*Tmp(i,j,k+3) &
                        +a04(5)*Tmp(i,j,k+4)
          enddo
       enddo

       k=2
       do j=1,ny
          do i=1,nx
             dTz(i,j,k)= a13(1)*Tmp(i,j,k-1)+a13(2)*Tmp(i,j,k  ) &
                        +a13(3)*Tmp(i,j,k+1)+a13(4)*Tmp(i,j,k+2) &
                        +a13(5)*Tmp(i,j,k+3)
          enddo
       enddo
    endif

    do k=ndz_v3,nfz_v3
       do j=1,ny
          do i=1,nx
             dTz(i,j,k) = a5(1)*(Tmp(i,j,k+1)-Tmp(i,j,k-1)) &
                        + a5(2)*(Tmp(i,j,k+2)-Tmp(i,j,k-2))
          enddo
       enddo
    enddo

    if (is_boundary(3,2)) then
       k=nz-1
       do j=1,ny
          do i=1,nx
             dTz(i,j,k)= a31(5)*Tmp(i,j,k-3)+a31(4)*Tmp(i,j,k-2) &
                        +a31(3)*Tmp(i,j,k-1)+a31(2)*Tmp(i,j,k  ) &
                        +a31(1)*Tmp(i,j,k+1)
          enddo
       enddo

       k=nz
       do j=1,ny
          do i=1,nx
             dTz(i,j,k)= a40(5)*Tmp(i,j,k-4)+a40(4)*Tmp(i,j,k-3) &
                        +a40(3)*Tmp(i,j,k-2)+a40(2)*Tmp(i,j,k-1) &
                        +a40(1)*Tmp(i,j,k  )
          enddo
       enddo
    endif

    ! store phi-derivatives on boundaries to enforce adiabaticity
    dTphi_imin=dTz(1,:,:)
    dTphi_imax=dTz(nx,:,:)
    dTphi_jmin=dTz(:,1,:)
    dTphi_jmax=dTz(:,ny,:)

    ! Compute derivatives of T in Cartesian coordinates
    ! =================================================
    do k=1,nz
       do j=1,ny
          do i=1,nx
             dvardx=dTx(i,j,k)*ksi_x_v(i,j,k)+dTy(i,j,k)*eta_x_v(i,j,k)+dTz(i,j,k)*phi_x_v(i,j,k)
             dvardy=dTx(i,j,k)*ksi_y_v(i,j,k)+dTy(i,j,k)*eta_y_v(i,j,k)+dTz(i,j,k)*phi_y_v(i,j,k)
             dvardz=dTx(i,j,k)*ksi_z_v(i,j,k)+dTy(i,j,k)*eta_z_v(i,j,k)+dTz(i,j,k)*phi_z_v(i,j,k)
             ! dT/dx
             dTx(i,j,k)=dvardx*ijacob3_v(i,j,k)
             ! dT/dy
             dTy(i,j,k)=dvardy*ijacob3_v(i,j,k)
             ! dT/dz
             dTz(i,j,k)=dvardz*ijacob3_v(i,j,k)
          enddo
       enddo
    enddo

    ! enforce adiabaticity condition at walls
    ! ---------------------------------------
    ! BC imin
    if ((BC_face(1,1)%sort==0).and.(is_adiab)) then
       do k=1,nz
          do j=1,ny
             ! tangent flux along x
             dTx(1,j,k)=dTeta_imin(j,k)*txn_eta_imin(j,k)+dTphi_imin(j,k)*txn_phi_imin(j,k)
             ! tangent flux along y
             dTy(1,j,k)=dTeta_imin(j,k)*tyn_eta_imin(j,k)+dTphi_imin(j,k)*tyn_phi_imin(j,k)
             ! tangent flux along y
             dTz(1,j,k)=dTeta_imin(j,k)*tzn_eta_imin(j,k)+dTphi_imin(j,k)*tzn_phi_imin(j,k)
          enddo
       enddo
    endif

    ! BC imax
    if ((BC_face(1,2)%sort==0).and.(is_adiab)) then
       do k=1,nz
          do j=1,ny
             ! tangent flux along x
             dTx(nx,j,k)=dTeta_imax(j,k)*txn_eta_imax(j,k)+dTphi_imax(j,k)*txn_phi_imax(j,k)
             ! tangent flux along y
             dTy(nx,j,k)=dTeta_imax(j,k)*tyn_eta_imax(j,k)+dTphi_imax(j,k)*tyn_phi_imax(j,k)
             ! tangent flux along y
             dTz(nx,j,k)=dTeta_imax(j,k)*tzn_eta_imax(j,k)+dTphi_imax(j,k)*tzn_phi_imax(j,k)
          enddo
       enddo
    endif

    ! BC jmin
    if ((BC_face(2,1)%sort==0).and.(is_adiab)) then
       do k=1,nz
          do i=1,nx
             ! tangent flux along x
             dTx(i,1,k)=dTksi_jmin(i,k)*txn_ksi_jmin(i,k)+dTphi_jmin(i,k)*txn_phi_jmin(i,k)
             ! tangent flux along y
             dTy(i,1,k)=dTksi_jmin(i,k)*tyn_ksi_jmin(i,k)+dTphi_jmin(i,k)*tyn_phi_jmin(i,k)
             ! tangent flux along z
             dTz(i,1,k)=dTksi_jmin(i,k)*tzn_ksi_jmin(i,k)+dTphi_jmin(i,k)*tzn_phi_jmin(i,k)
          enddo
       enddo
    endif

    ! BC jmax
    if ((BC_face(2,2)%sort==0).and.(is_adiab)) then
       do k=1,nz
          do i=1,nx
             ! tangent flux along x
             dTx(i,ny,k)=dTksi_jmax(i,k)*txn_ksi_jmax(i,k)+dTphi_jmax(i,k)*txn_phi_jmax(i,k)
             ! tangent flux along y
             dTy(i,ny,k)=dTksi_jmax(i,k)*tyn_ksi_jmax(i,k)+dTphi_jmax(i,k)*tyn_phi_jmax(i,k)
             ! tangent flux along z
             dTz(i,ny,k)=dTksi_jmax(i,k)*tzn_ksi_jmax(i,k)+dTphi_jmax(i,k)*tzn_phi_jmax(i,k)
          enddo
       enddo
    endif

    ! BC kmin
    if ((BC_face(3,1)%sort==0).and.(is_adiab)) then
       do j=1,ny
          do i=1,nx
             ! tangent flux along x
             dTx(i,j,1)=dTksi_kmin(i,j)*txn_ksi_kmin(i,j)+dTeta_kmin(i,j)*txn_eta_kmin(i,j)
             ! tangent flux along y
             dTy(i,j,1)=dTksi_kmin(i,j)*tyn_ksi_kmin(i,j)+dTeta_kmin(i,j)*tyn_eta_kmin(i,j)
             ! tangent flux along z
             dTz(i,j,1)=dTksi_kmin(i,j)*tzn_ksi_kmin(i,j)+dTeta_kmin(i,j)*tzn_eta_kmin(i,j)
          enddo
       enddo
    endif

    ! BC kmax
    if ((BC_face(3,2)%sort==0).and.(is_adiab)) then
       do j=1,ny
          do i=1,nx
             ! tangent flux along x
             dTx(i,j,nz)=dTksi_kmax(i,j)*txn_ksi_kmax(i,j)+dTeta_kmax(i,j)*txn_eta_kmax(i,j)
             ! tangent flux along y
             dTy(i,j,nz)=dTksi_kmax(i,j)*tyn_ksi_kmax(i,j)+dTeta_kmax(i,j)*tyn_eta_kmax(i,j)
             ! tangent flux along z
             dTz(i,j,nz)=dTksi_kmax(i,j)*tzn_ksi_kmax(i,j)+dTeta_kmax(i,j)*tzn_eta_kmax(i,j)
          enddo
       enddo
    endif

  end subroutine grad_T_5pts_c3

  !===============================================================================
  module subroutine grad_T_3pts_c3
  !===============================================================================
    !> Compute temperature derivatives
    !> 3-point stencil - 3D curvilinear version
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: dvardx,dvardy,dvardz
    real(wp), dimension(ny1_v:ny2_v,nz1_v:nz2_v) :: dTeta_imin,dTeta_imax
    real(wp), dimension(ny1_v:ny2_v,nz1_v:nz2_v) :: dTphi_imin,dTphi_imax
    real(wp), dimension(nx1_v:nx2_v,nz1_v:nz2_v) :: dTksi_jmin,dTksi_jmax
    real(wp), dimension(nx1_v:nx2_v,nz1_v:nz2_v) :: dTphi_jmin,dTphi_jmax
    real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v) :: dTksi_kmin,dTksi_kmax
    real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v) :: dTeta_kmin,dTeta_kmax
    ! ---------------------------------------------------------------------------

    ! Derivatives of T along ksi (Caution: dTx is dT/dksi)
    ! ==========================
    if (is_boundary(1,1)) then
       i=1
       do k=1,nz
          do j=1,ny
             dTx(i,j,k)= a02(1)*Tmp(i  ,j,k)+a02(2)*Tmp(i+1,j,k) &
                        +a02(3)*Tmp(i+2,j,k)
          enddo
       enddo
    endif

    do k=1,nz
       do j=1,ny
          do i=ndx_v3,nfx_v3
             dTx(i,j,k) = 0.5_wp*(Tmp(i+1,j,k)-Tmp(i-1,j,k))
          enddo
       enddo
    enddo

    if (is_boundary(1,2)) then
       i=nx
       do k=1,nz
          do j=1,ny
             dTx(i,j,k)= a20(3)*Tmp(i-2,j,k)+a20(2)*Tmp(i-1,j,k) &
                        +a20(1)*Tmp(i  ,j,k)
          enddo
       enddo
    endif

    ! store ksi-derivatives on boundaries to enforce adiabaticity
    dTksi_jmin=dTx(:,1,:)
    dTksi_jmax=dTx(:,ny,:)
    dTksi_kmin=dTx(:,:,1)
    dTksi_kmax=dTx(:,:,nz)

    ! Derivatives of T along eta (Caution: dTy is dT/deta)
    ! ==========================
    if (is_boundary(2,1)) then
       j=1
       do k=1,nz
          do i=1,nx
             dTy(i,j,k)= a02(1)*Tmp(i,j  ,k)+a02(2)*Tmp(i,j+1,k) &
                        +a02(3)*Tmp(i,j+2,k)
          enddo
       enddo
    endif

    do k=1,nz
       do j=ndy_v3,nfy_v3
          do i=1,nx
             dTy(i,j,k) = 0.5_wp*(Tmp(i,j+1,k)-Tmp(i,j-1,k))
          enddo
       enddo
    enddo

    if (is_boundary(2,2)) then
       j=ny
       do k=1,nz
          do i=1,nx
             dTy(i,j,k)= a20(3)*Tmp(i,j-2,k)+a20(2)*Tmp(i,j-1,k) &
                        +a20(1)*Tmp(i,j  ,k)
          enddo
       enddo
    endif

    ! store eta-derivatives on boundaries to enforce adiabaticity
    dTeta_imin=dTy(1,:,:)
    dTeta_imax=dTy(nx,:,:)
    dTeta_kmin=dTy(:,:,1)
    dTeta_kmax=dTy(:,:,nz)

    ! Derivatives of T along phi (Caution: dTz is dT/deta)
    ! ==========================
    if (is_boundary(3,1)) then
       k=1
       do j=1,ny
          do i=1,nx
             dTz(i,j,k)= a02(1)*Tmp(i,j,k  )+a02(2)*Tmp(i,j,k+1) &
                        +a02(3)*Tmp(i,j,k+2)
          enddo
       enddo
    endif

    do k=ndz_v3,nfz_v3
       do j=1,ny
          do i=1,nx
             dTz(i,j,k) = 0.5_wp*(Tmp(i,j,k+1)-Tmp(i,j,k-1))
          enddo
       enddo
    enddo

    if (is_boundary(3,2)) then
       k=nz
       do j=1,ny
          do i=1,nx
             dTz(i,j,k)= a20(3)*Tmp(i,j,k-2)+a20(2)*Tmp(i,j,k-1) &
                        +a20(1)*Tmp(i,j,k  )
          enddo
       enddo
    endif

    ! store phi-derivatives on boundaries to enforce adiabaticity
    dTphi_imin=dTz(1,:,:)
    dTphi_imax=dTz(nx,:,:)
    dTphi_jmin=dTz(:,1,:)
    dTphi_jmax=dTz(:,ny,:)

    ! Compute derivatives of T in Cartesian coordinates
    ! =================================================
    do k=1,nz
       do j=1,ny
          do i=1,nx
             dvardx=dTx(i,j,k)*ksi_x_v(i,j,k)+dTy(i,j,k)*eta_x_v(i,j,k)+dTz(i,j,k)*phi_x_v(i,j,k)
             dvardy=dTx(i,j,k)*ksi_y_v(i,j,k)+dTy(i,j,k)*eta_y_v(i,j,k)+dTz(i,j,k)*phi_y_v(i,j,k)
             dvardz=dTx(i,j,k)*ksi_z_v(i,j,k)+dTy(i,j,k)*eta_z_v(i,j,k)+dTz(i,j,k)*phi_z_v(i,j,k)
             ! dT/dx
             dTx(i,j,k)=dvardx*ijacob3_v(i,j,k)
             ! dT/dy
             dTy(i,j,k)=dvardy*ijacob3_v(i,j,k)
             ! dT/dz
             dTz(i,j,k)=dvardz*ijacob3_v(i,j,k)
          enddo
       enddo
    enddo

    ! enforce adiabaticity condition at walls
    ! ---------------------------------------
    ! BC imin
    if ((BC_face(1,1)%sort==0).and.(is_adiab)) then
       do k=1,nz
          do j=1,ny
             ! tangent flux along x
             dTx(1,j,k)=dTeta_imin(j,k)*txn_eta_imin(j,k)+dTphi_imin(j,k)*txn_phi_imin(j,k)
             ! tangent flux along y
             dTy(1,j,k)=dTeta_imin(j,k)*tyn_eta_imin(j,k)+dTphi_imin(j,k)*tyn_phi_imin(j,k)
             ! tangent flux along y
             dTz(1,j,k)=dTeta_imin(j,k)*tzn_eta_imin(j,k)+dTphi_imin(j,k)*tzn_phi_imin(j,k)
          enddo
       enddo
    endif

    ! BC imax
    if ((BC_face(1,2)%sort==0).and.(is_adiab)) then
       do k=1,nz
          do j=1,ny
             ! tangent flux along x
             dTx(nx,j,k)=dTeta_imax(j,k)*txn_eta_imax(j,k)+dTphi_imax(j,k)*txn_phi_imax(j,k)
             ! tangent flux along y
             dTy(nx,j,k)=dTeta_imax(j,k)*tyn_eta_imax(j,k)+dTphi_imax(j,k)*tyn_phi_imax(j,k)
             ! tangent flux along y
             dTz(nx,j,k)=dTeta_imax(j,k)*tzn_eta_imax(j,k)+dTphi_imax(j,k)*tzn_phi_imax(j,k)
          enddo
       enddo
    endif

    ! BC jmin
    if ((BC_face(2,1)%sort==0).and.(is_adiab)) then
       do k=1,nz
          do i=1,nx
             ! tangent flux along x
             dTx(i,1,k)=dTksi_jmin(i,k)*txn_ksi_jmin(i,k)+dTphi_jmin(i,k)*txn_phi_jmin(i,k)
             ! tangent flux along y
             dTy(i,1,k)=dTksi_jmin(i,k)*tyn_ksi_jmin(i,k)+dTphi_jmin(i,k)*tyn_phi_jmin(i,k)
             ! tangent flux along z
             dTz(i,1,k)=dTksi_jmin(i,k)*tzn_ksi_jmin(i,k)+dTphi_jmin(i,k)*tzn_phi_jmin(i,k)
          enddo
       enddo
    endif

    ! BC jmax
    if ((BC_face(2,2)%sort==0).and.(is_adiab)) then
       do k=1,nz
          do i=1,nx
             ! tangent flux along x
             dTx(i,ny,k)=dTksi_jmax(i,k)*txn_ksi_jmax(i,k)+dTphi_jmax(i,k)*txn_phi_jmax(i,k)
             ! tangent flux along y
             dTy(i,ny,k)=dTksi_jmax(i,k)*tyn_ksi_jmax(i,k)+dTphi_jmax(i,k)*tyn_phi_jmax(i,k)
             ! tangent flux along z
             dTz(i,ny,k)=dTksi_jmax(i,k)*tzn_ksi_jmax(i,k)+dTphi_jmax(i,k)*tzn_phi_jmax(i,k)
          enddo
       enddo
    endif

    ! BC kmin
    if ((BC_face(3,1)%sort==0).and.(is_adiab)) then
       do j=1,ny
          do i=1,nx
             ! tangent flux along x
             dTx(i,j,1)=dTksi_kmin(i,j)*txn_ksi_kmin(i,j)+dTeta_kmin(i,j)*txn_eta_kmin(i,j)
             ! tangent flux along y
             dTy(i,j,1)=dTksi_kmin(i,j)*tyn_ksi_kmin(i,j)+dTeta_kmin(i,j)*tyn_eta_kmin(i,j)
             ! tangent flux along z
             dTz(i,j,1)=dTksi_kmin(i,j)*tzn_ksi_kmin(i,j)+dTeta_kmin(i,j)*tzn_eta_kmin(i,j)
          enddo
       enddo
    endif

    ! BC kmax
    if ((BC_face(3,2)%sort==0).and.(is_adiab)) then
       do j=1,ny
          do i=1,nx
             ! tangent flux along x
             dTx(i,j,nz)=dTksi_kmax(i,j)*txn_ksi_kmax(i,j)+dTeta_kmax(i,j)*txn_eta_kmax(i,j)
             ! tangent flux along y
             dTy(i,j,nz)=dTksi_kmax(i,j)*tyn_ksi_kmax(i,j)+dTeta_kmax(i,j)*tyn_eta_kmax(i,j)
             ! tangent flux along z
             dTz(i,j,nz)=dTksi_kmax(i,j)*tzn_ksi_kmax(i,j)+dTeta_kmax(i,j)*tzn_eta_kmax(i,j)
          enddo
       enddo
    endif

  end subroutine grad_T_3pts_c3

end submodule smod_grad_T_c3
