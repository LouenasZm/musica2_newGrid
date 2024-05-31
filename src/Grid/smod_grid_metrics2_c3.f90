!=================================================================================
submodule (mod_grid_metrics_c3) smod_grid_metrics2_c3
!=================================================================================
  !> Module to compute full 3D curvilinear metrics
  !> with Geometric Conservation Law (GCL)
!=================================================================================

contains

  !===============================================================================
  module subroutine correct_deriv_sign(dksi,deta,dphi)
  !===============================================================================
    !> Correct metrics derivatives due to swap and reverse in neighbors
    !===============================================================================
    use mod_grid_directions
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(inout) :: dksi,deta,dphi
    ! ---------------------------------------------------------------------------
    integer :: sgn_ksi,sgn_eta,sgn_phi
    ! ---------------------------------------------------------------------------

    ! Correct values in interface 1 (imin)
    ! ====================================
    sgn_ksi=1
    if (is_rev2(1,2)) sgn_ksi=-1
    dksi(1-ngh:0,:,:)=sgn_ksi*dksi(1-ngh:0,:,:)

    ! Correct values in interface 2 (imax)
    ! ====================================
    sgn_ksi=1
    if (is_rev2(2,2)) sgn_ksi=-1
    dksi(nx+1:nx+ngh,:,:)=sgn_ksi*dksi(nx+1:nx+ngh,:,:)
    
    ! Correct values in interface 3 (jmin)
    ! ====================================
    sgn_eta=1
    if (is_rev2(3,2)) sgn_eta=-1
    deta(:,1-ngh:0,:)=sgn_eta*deta(:,1-ngh:0,:)

    ! Correct values in interface 4 (jmax)
    ! ====================================
    sgn_eta=1
    if (is_rev2(4,2)) sgn_eta=-1
    deta(:,ny+1:ny+ngh,:)=sgn_eta*deta(:,ny+1:ny+ngh,:)

    ! Correct values in interface 5 (kmin)
    ! ====================================
    sgn_phi=1
    if (is_rev2(5,2)) sgn_phi=-1
    dphi(:,:,1-ngh:0)=sgn_phi*dphi(:,:,1-ngh:0)

    ! Correct values in interface 6 (kmax)
    ! ====================================
    sgn_phi=1
    if (is_rev2(6,2)) sgn_phi=-1
    dphi(:,:,nz+1:nz+ngh)=sgn_phi*dphi(:,:,nz+1:nz+ngh)

  end subroutine correct_deriv_sign
  
  !===============================================================================
  module subroutine correct_deriv_sign_v(dksi,deta,dphi)
  !===============================================================================
    !> Correct metrics derivatives due to swap and reverse in neighbors
    !===============================================================================
    use mod_grid_directions
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v), intent(inout) :: dksi,deta,dphi
    ! ---------------------------------------------------------------------------
    integer :: sgn_ksi,sgn_eta,sgn_phi
    ! ---------------------------------------------------------------------------

    ! Correct values in interface 1 (imin)
    ! ====================================
    sgn_ksi=1
    if (is_rev2(1,2)) sgn_ksi=-1
    dksi(-1:0,:,:)=sgn_ksi*dksi(-1:0,:,:)

    ! Correct values in interface 2 (imax)
    ! ====================================
    sgn_ksi=1
    if (is_rev2(2,2)) sgn_ksi=-1
    dksi(nx+1:nx+2,:,:)=sgn_ksi*dksi(nx+1:nx+2,:,:)
    
    ! Correct values in interface 3 (jmin)
    ! ====================================
    sgn_eta=1
    if (is_rev2(3,2)) sgn_eta=-1
    deta(:,-1:0,:)=sgn_eta*deta(:,-1:0,:)

    ! Correct values in interface 4 (jmax)
    ! ====================================
    sgn_eta=1
    if (is_rev2(4,2)) sgn_eta=-1
    deta(:,ny+1:ny+2,:)=sgn_eta*deta(:,ny+1:ny+2,:)

    ! Correct values in interface 5 (kmin)
    ! ====================================
    sgn_phi=1
    if (is_rev2(5,2)) sgn_phi=-1
    dphi(:,:,-1:0)=sgn_phi*dphi(:,:,-1:0)

    ! Correct values in interface 6 (kmax)
    ! ====================================
    sgn_phi=1
    if (is_rev2(6,2)) sgn_phi=-1
    dphi(:,:,nz+1:nz+2)=sgn_phi*dphi(:,:,nz+1:nz+2)

  end subroutine correct_deriv_sign_v

  !===============================================================================
  module subroutine correct_dderiv_sign(dksi1,deta1,dphi1,dksi2,deta2,dphi2)
  !===============================================================================
    !> Correct metrics double derivatives due to swap and reverse in neighbors
    !===============================================================================
    use mod_grid_directions
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(inout) :: dksi1,deta1,dphi1
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(inout) :: dksi2,deta2,dphi2
    ! ---------------------------------------------------------------------------
    integer :: sgn_ksi,sgn_eta,sgn_phi
    ! ---------------------------------------------------------------------------

    ! Correct values in interface 1 (imin)
    ! ====================================
    sgn_ksi=1
    sgn_eta=1
    sgn_phi=1
    ! correct signs for reverse direction
    if (is_rev2(1,2)) sgn_ksi=-1
    if (is_rev2(1,1)) sgn_eta=-1
    if (is_rev2(1,3)) sgn_phi=-1
    deta1(1-ngh:0,:,:)=sgn_phi*sgn_ksi*deta1(1-ngh:0,:,:)
    dphi2(1-ngh:0,:,:)=sgn_ksi*sgn_eta*dphi2(1-ngh:0,:,:)

    ! Correct values in interface 2 (imax)
    ! ====================================
    sgn_ksi=1
    sgn_eta=1
    sgn_phi=1
    ! correct signs for reverse direction
    if (is_rev2(2,2)) sgn_ksi=-1
    if (is_rev2(2,1)) sgn_eta=-1
    if (is_rev2(2,3)) sgn_phi=-1
    deta1(nx+1:nx+ngh,:,:)=sgn_phi*sgn_ksi*deta1(nx+1:nx+ngh,:,:)
    dphi2(nx+1:nx+ngh,:,:)=sgn_ksi*sgn_eta*dphi2(nx+1:nx+ngh,:,:)

    ! Correct values in interface 3 (jmin)
    ! ====================================
    sgn_ksi=1
    sgn_eta=1
    sgn_phi=1
    ! correct signs for reverse direction
    if (is_rev2(3,1)) sgn_ksi=-1
    if (is_rev2(3,2)) sgn_eta=-1
    if (is_rev2(3,3)) sgn_phi=-1
    dphi1(:,1-ngh:0,:)=sgn_ksi*sgn_eta*dphi1(:,1-ngh:0,:)
    dksi2(:,1-ngh:0,:)=sgn_eta*sgn_phi*dksi2(:,1-ngh:0,:)

    ! Correct values in interface 4 (jmax)
    ! ====================================
    sgn_ksi=1
    sgn_eta=1
    sgn_phi=1
    ! correct signs for reverse direction
    if (is_rev2(4,1)) sgn_ksi=-1
    if (is_rev2(4,2)) sgn_eta=-1
    if (is_rev2(4,3)) sgn_phi=-1
    dphi1(:,ny+1:ny+ngh,:)=sgn_ksi*sgn_eta*dphi1(:,ny+1:ny+ngh,:)
    dksi2(:,ny+1:ny+ngh,:)=sgn_eta*sgn_phi*dksi2(:,ny+1:ny+ngh,:)

    ! Correct values in interface 5 (kmin)
    ! ====================================
    sgn_ksi=1
    sgn_eta=1
    sgn_phi=1
    ! correct signs for reverse direction
    if (is_rev2(5,1)) sgn_ksi=-1
    if (is_rev2(5,3)) sgn_eta=-1
    if (is_rev2(5,2)) sgn_phi=-1
    dksi1(:,:,1-ngh:0)=sgn_eta*sgn_phi*dksi1(:,:,1-ngh:0)
    deta2(:,:,1-ngh:0)=sgn_phi*sgn_ksi*deta2(:,:,1-ngh:0)
    
    ! Correct values in interface 6 (kmax)
    ! ====================================
    sgn_ksi=1
    sgn_eta=1
    sgn_phi=1
    ! correct signs for reverse direction
    if (is_rev2(6,1)) sgn_ksi=-1
    if (is_rev2(6,3)) sgn_eta=-1
    if (is_rev2(6,2)) sgn_phi=-1
    dksi1(:,:,nz+1:nz+ngh)=sgn_eta*sgn_phi*dksi1(:,:,nz+1:nz+ngh)
    deta2(:,:,nz+1:nz+ngh)=sgn_phi*sgn_ksi*deta2(:,:,nz+1:nz+ngh)

  end subroutine correct_dderiv_sign
 
  !===============================================================================
  module subroutine correct_dderiv_sign_v(dksi1,deta1,dphi1,dksi2,deta2,dphi2)
  !===============================================================================
    !> Correct metrics double derivatives due to swap and reverse in neighbors
    !===============================================================================
    use mod_grid_directions
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v), intent(inout) :: dksi1,deta1,dphi1
    real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v), intent(inout) :: dksi2,deta2,dphi2
    ! ---------------------------------------------------------------------------
    integer :: sgn_ksi,sgn_eta,sgn_phi
    ! ---------------------------------------------------------------------------

    ! Correct values in interface 1 (imin)
    ! ====================================
    sgn_ksi=1
    sgn_eta=1
    sgn_phi=1
    ! correct signs for reverse direction
    if (is_rev2(1,2)) sgn_ksi=-1
    if (is_rev2(1,1)) sgn_eta=-1
    if (is_rev2(1,3)) sgn_phi=-1
    deta1(-1:0,:,:)=sgn_phi*sgn_ksi*deta1(-1:0,:,:)
    dphi2(-1:0,:,:)=sgn_ksi*sgn_eta*dphi2(-1:0,:,:)

    ! Correct values in interface 2 (imax)
    ! ====================================
    sgn_ksi=1
    sgn_eta=1
    sgn_phi=1
    ! correct signs for reverse direction
    if (is_rev2(2,2)) sgn_ksi=-1
    if (is_rev2(2,1)) sgn_eta=-1
    if (is_rev2(2,3)) sgn_phi=-1
    deta1(nx+1:nx+2,:,:)=sgn_phi*sgn_ksi*deta1(nx+1:nx+2,:,:)
    dphi2(nx+1:nx+2,:,:)=sgn_ksi*sgn_eta*dphi2(nx+1:nx+2,:,:)

    ! Correct values in interface 3 (jmin)
    ! ====================================
    sgn_ksi=1
    sgn_eta=1
    sgn_phi=1
    ! correct signs for reverse direction
    if (is_rev2(3,1)) sgn_ksi=-1
    if (is_rev2(3,2)) sgn_eta=-1
    if (is_rev2(3,3)) sgn_phi=-1
    dphi1(:,-1:0,:)=sgn_ksi*sgn_eta*dphi1(:,-1:0,:)
    dksi2(:,-1:0,:)=sgn_eta*sgn_phi*dksi2(:,-1:0,:)

    ! Correct values in interface 4 (jmax)
    ! ====================================
    sgn_ksi=1
    sgn_eta=1
    sgn_phi=1
    ! correct signs for reverse direction
    if (is_rev2(4,1)) sgn_ksi=-1
    if (is_rev2(4,2)) sgn_eta=-1
    if (is_rev2(4,3)) sgn_phi=-1
    dphi1(:,ny+1:ny+2,:)=sgn_ksi*sgn_eta*dphi1(:,ny+1:ny+2,:)
    dksi2(:,ny+1:ny+2,:)=sgn_eta*sgn_phi*dksi2(:,ny+1:ny+2,:)

    ! Correct values in interface 5 (kmin)
    ! ====================================
    sgn_ksi=1
    sgn_eta=1
    sgn_phi=1
    ! correct signs for reverse direction
    if (is_rev2(5,1)) sgn_ksi=-1
    if (is_rev2(5,3)) sgn_eta=-1
    if (is_rev2(5,2)) sgn_phi=-1
    dksi1(:,:,-1:0)=sgn_eta*sgn_phi*dksi1(:,:,-1:0)
    deta2(:,:,-1:0)=sgn_phi*sgn_ksi*deta2(:,:,-1:0)
    
    ! Correct values in interface 6 (kmax)
    ! ====================================
    sgn_ksi=1
    sgn_eta=1
    sgn_phi=1
    ! correct signs for reverse direction
    if (is_rev2(6,1)) sgn_ksi=-1
    if (is_rev2(6,3)) sgn_eta=-1
    if (is_rev2(6,2)) sgn_phi=-1
    dksi1(:,:,nz+1:nz+2)=sgn_eta*sgn_phi*dksi1(:,:,nz+1:nz+2)
    deta2(:,:,nz+1:nz+2)=sgn_phi*sgn_ksi*deta2(:,:,nz+1:nz+2)

  end subroutine correct_dderiv_sign_v
 
  !===============================================================================
  module subroutine correct_edges_i
  !===============================================================================
    !> reconstruct ghost cells of degenerate edges for i direction
  !===============================================================================
    use mod_block
    use mod_mpi_part
    use mod_grid_directions
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k
    integer :: n,nv,nv_
    ! ----------------------------------------------------------------------------

    ! current block number
    ! -------------------
    n=nob(iproc)

    ! Reconstruct ghost cells of degenerate edges only for j directions
    ! ------------------------------------------------------------------

    !print *,'================================'
    !print *,'BC',bl(n)%BC

    ! imin-jmin
    ! ---------
    nv=bl(n)%BC(1) ! block imin
    if (nv>0) then
       if (is_swapij2_bl(1)) then
          if (is_rev2(1,1)) then
             nv_=bl(nv)%BC(2) ! block imax of imin
          else
             nv_=bl(nv)%BC(1) ! block imin of imin
          endif
       else
          nv_=bl(nv)%BC(3) ! block jmin of imin
       endif
    else
       nv_=nv
    endif
    !print *,'imin-jmin',nv,nv_,bl(n)%BC(3),bl(n)%BC(4)
    !print *,'proc',coord(1),coord(2)
    if ((nv_>0).and.((nv_==bl(n)%BC(3)).or.(nv_==bl(n)%BC(4))).and.((coord(1)==0).and.(coord(2)==0))) then
       !print *,'do it'
       do k=ndzt,nfzt
          do j=1-ngh,0
             do i=1-ngh,0
                xc3(i,j,k)=xc3(j,1-i,k)
                yc3(i,j,k)=yc3(j,1-i,k)
                zc3(i,j,k)=zc3(j,1-i,k)
             enddo
          enddo
       enddo
    endif

    ! imin-jmax
    ! ---------
    nv=bl(n)%BC(1) ! block imin
    if (nv>0) then
       if (is_swapij2_bl(1)) then
          if (is_rev2(1,1)) then
             nv_=bl(nv)%BC(1) ! block imin of imin
          else
             nv_=bl(nv)%BC(2) ! block imax of imin
          endif
       else
          nv_=bl(nv)%BC(4) ! block jmax of imin
       endif
    else
       nv_=nv
    endif
    !print *,'imin-jmax',nv,nv_,bl(n)%BC(3),bl(n)%BC(4)
    !print *,'proc',coord(1),coord(2)
    if ((nv_>0).and.((nv_==bl(n)%BC(3)).or.(nv_==bl(n)%BC(4))).and.((coord(1)==0).and.(coord(2)==ndomy-1))) then
       !print *,'do it'
       do k=ndzt,nfzt
          do j=1,ngh
             do i=1-ngh,0
                xc3(i,ny+j,k)=xc3(1-j,ny+i,k)
                yc3(i,ny+j,k)=yc3(1-j,ny+i,k)
                zc3(i,ny+j,k)=zc3(1-j,ny+i,k)
             enddo
          enddo
       enddo
    endif

    ! imax-jmin
    ! ---------
    nv=bl(n)%BC(2) ! block imax
    if (nv>0) then
       if (is_swapij2_bl(2)) then
          if (is_rev2(2,1)) then
             nv_=bl(nv)%BC(2) ! block imax of imax
          else
             nv_=bl(nv)%BC(1) ! block imin of imax
          endif
       else
          nv_=bl(nv)%BC(3) ! block jmin of imax
       endif
    else
       nv_=nv
    endif
    !print *,'imax-jmin',nv,nv_,bl(n)%BC(3),bl(n)%BC(4)
    !print *,'proc',coord(1),coord(2)
    if ((nv_>0).and.((nv_==bl(n)%BC(3)).or.(nv_==bl(n)%BC(4))).and.((coord(1)==ndomx-1).and.(coord(2)==0))) then
       !print *,'do it'
       do k=ndzt,nfzt
          do j=1-ngh,0
             do i=1,ngh
                xc3(nx+i,j,k)=xc3(nx-j+1,i,k)
                yc3(nx+i,j,k)=yc3(nx-j+1,i,k)
                zc3(nx+i,j,k)=zc3(nx-j+1,i,k)
             enddo
          enddo
       enddo
    endif

    ! imax-jmax
    ! ---------
    nv=bl(n)%BC(2) ! block imax
    if (nv>0) then
       if (is_swapij2_bl(2)) then
          if (is_rev2(2,1)) then
             nv_=bl(nv)%BC(1) ! block imin of imax
          else
             nv_=bl(nv)%BC(2) ! block imax of imax
          endif
       else
          nv_=bl(nv)%BC(4) ! block jmax of imax
       endif
    else
       nv_=nv
    endif
    !print *,'imax-jmax',nv,nv_,bl(n)%BC(3),bl(n)%BC(4)
    !print *,'proc',coord(1),coord(2)
    if ((nv_>0).and.((nv_==bl(n)%BC(3)).or.(nv_==bl(n)%BC(4))).and.((coord(1)==ndomx-1).and.(coord(2)==ndomy-1))) then
       !print *,'do it'
       do k=ndzt,nfzt
          do j=1,ngh
             do i=1,ngh
                xc3(nx+i,ny+j,k)=xc3(nx+j,ny-i+1,k)
                yc3(nx+i,ny+j,k)=yc3(nx+j,ny-i+1,k)
                zc3(nx+i,ny+j,k)=zc3(nx+j,ny-i+1,k)
             enddo
          enddo
       enddo
    endif

    ! jmin-kmin
    ! ---------
    nv=bl(n)%BC(3) ! block jmin
    if (nv>0) then
       if (is_swapjk2_bl(3)) then
          if (is_rev2(3,3)) then
             nv_=bl(nv)%BC(4) ! block jmax of jmin
          else
             nv_=bl(nv)%BC(3) ! block jmin of jmin
          endif
       else
          nv_=bl(nv)%BC(5) ! block kmin of jmin
       endif
    else
       nv_=nv
    endif
    !print *,'jmin-kmin',nv,nv_,bl(n)%BC(5),bl(n)%BC(6)
    !print *,'proc',coord(2),coord(3)
    if ((nv_>0).and.((nv_==bl(n)%BC(5)).or.(nv_==bl(n)%BC(6))).and.((coord(2)==0).and.(coord(3)==0))) then
       !print *,'do it'
       do i=ndxt,nfxt
          do k=1-ngh,0
             do j=1-ngh,0
                xc3(i,j,k)=xc3(i,k,1-j)
                yc3(i,j,k)=yc3(i,k,1-j)
                zc3(i,j,k)=zc3(i,k,1-j)
             enddo
          enddo
       enddo
    endif

    ! jmin-kmax
    ! ---------
    nv=bl(n)%BC(3) ! block jmin
    if (nv>0) then
       if (is_swapjk2_bl(3)) then
          if (is_rev2(3,3)) then
             nv_=bl(nv)%BC(3) ! block jmin of jmin
          else
             nv_=bl(nv)%BC(4) ! block jmax of jmin
          endif
       else
          nv_=bl(nv)%BC(6) ! block kmax of jmin
       endif
    else
       nv_=nv
    endif
    !print *,'jmin-kmax',nv,nv_,bl(n)%BC(5),bl(n)%BC(6)
    !print *,'proc',coord(2),coord(3)
    if ((nv_>0).and.((nv_==bl(n)%BC(5)).or.(nv_==bl(n)%BC(6))).and.((coord(2)==0).and.(coord(3)==ndomz-1))) then
       !print *,'do it'
       do i=ndxt,nfxt
          do k=1,ngh
             do j=1-ngh,0
                xc3(i,j,nz+k)=xc3(i,1-k,nz+j)
                yc3(i,j,nz+k)=yc3(i,1-k,nz+j)
                zc3(i,j,nz+k)=zc3(i,1-k,nz+j)
             enddo
          enddo
       enddo
    endif

    ! jmax-kmin
    ! ---------
    nv=bl(n)%BC(4) ! block jmax
    if (nv>0) then
       if (is_swapjk2_bl(4)) then
          if (is_rev2(4,3)) then
             nv_=bl(nv)%BC(4) ! block jmax of jmax
          else
             nv_=bl(nv)%BC(3) ! block jmin of jmax
          endif
       else
          nv_=bl(nv)%BC(5) ! block kmin of jmax
       endif
    else
       nv_=nv
    endif
    !print *,'jmax-kmin',nv,nv_,bl(n)%BC(5),bl(n)%BC(6)
    !print *,'proc',coord(2),coord(3)
    if ((nv_>0).and.((nv_==bl(n)%BC(5)).or.(nv_==bl(n)%BC(6))).and.((coord(2)==ndomy-1).and.(coord(3)==0))) then
       !print *,'do it'
       do i=ndxt,nfxt
          do k=1-ngh,0
             do j=1,ngh
                xc3(i,ny+j,k)=xc3(i,ny-k+1,j)
                yc3(i,ny+j,k)=yc3(i,ny-k+1,j)
                zc3(i,ny+j,k)=zc3(i,ny-k+1,j)
             enddo
          enddo
       enddo
    endif

    ! jmax-kmax
    ! ---------
    nv=bl(n)%BC(4) ! block jmax
    if (nv>0) then
       if (is_swapjk2_bl(4)) then
          if (is_rev2(4,3)) then
             nv_=bl(nv)%BC(3) ! block jmin of jmax
          else
             nv_=bl(nv)%BC(4) ! block jmax of jmax
          endif
       else
          nv_=bl(nv)%BC(6) ! block kmax of jmax
       endif
    else
       nv_=nv
    endif
    !print *,'jmax-kmax',nv,nv_,bl(n)%BC(5),bl(n)%BC(6)
    !print *,'proc',coord(2),coord(3)
    if ((nv_>0).and.((nv_==bl(n)%BC(5)).or.(nv_==bl(n)%BC(6))).and.((coord(2)==ndomy-1).and.(coord(3)==ndomz-1))) then
       !print *,'do it'
       do i=ndxt,nfxt
          do k=1,ngh
             do j=1,ngh
                xc3(i,ny+j,nz+k)=xc3(i,ny+k,nz-j+1)
                yc3(i,ny+j,nz+k)=yc3(i,ny+k,nz-j+1)
                zc3(i,ny+j,nz+k)=zc3(i,ny+k,nz-j+1)
             enddo
          enddo
       enddo
    endif

    ! imin-kmin
    ! ---------
    nv=bl(n)%BC(1) ! block imin
    if (nv>0) then
       if (is_swapik2_bl(1)) then
          if (is_rev2(1,3)) then
             nv_=bl(nv)%BC(2) ! block imax of imin
          else
             nv_=bl(nv)%BC(1) ! block imin of imin
          endif
       else
          nv_=bl(nv)%BC(5) ! block kmin of imin
       endif
    else
       nv_=nv
    endif
    !print *,'imin-kmin',nv,nv_,bl(n)%BC(5),bl(n)%BC(6)
    !print *,'proc',coord(1),coord(3)
    if ((nv_>0).and.((nv_==bl(n)%BC(5)).or.(nv_==bl(n)%BC(6))).and.((coord(1)==0).and.(coord(3)==0))) then
       !print *,'do it'
       do j=ndyt,nfyt
          do k=1-ngh,0
             do i=1-ngh,0
                xc3(i,j,k)=xc3(k,j,1-i)
                yc3(i,j,k)=yc3(k,j,1-i)
                zc3(i,j,k)=zc3(k,j,1-i)
             enddo
          enddo
       enddo
    endif

    ! imin-kmax
    ! ---------
    nv=bl(n)%BC(1) ! block imin
    if (nv>0) then
       if (is_swapik2_bl(1)) then
          if (is_rev2(1,3)) then
             nv_=bl(nv)%BC(1) ! block imin of imin
          else
             nv_=bl(nv)%BC(2) ! block imax of imin
          endif
       else
          nv_=bl(nv)%BC(6) ! block kmax of imin
       endif
    else
       nv_=nv
    endif
    !print *,'imin-kmax',nv,nv_,bl(n)%BC(5),bl(n)%BC(6)
    !print *,'proc',coord(1),coord(3)
    if ((nv_>0).and.((nv_==bl(n)%BC(5)).or.(nv_==bl(n)%BC(6))).and.((coord(1)==0).and.(coord(3)==ndomz-1))) then
       !print *,'do it'
       do j=ndyt,nfyt
          do k=1,ngh
             do i=1-ngh,0
                xc3(i,j,nz+k)=xc3(1-k,j,nz+i)
                yc3(i,j,nz+k)=yc3(1-k,j,nz+i)
                zc3(i,j,nz+k)=zc3(1-k,j,nz+i)
             enddo
          enddo
       enddo
    endif

    ! imax-kmin
    ! ---------
    nv=bl(n)%BC(2) ! block imax
    if (nv>0) then
       if (is_swapik2_bl(2)) then
          if (is_rev2(2,3)) then
             nv_=bl(nv)%BC(2) ! block imax of imax
          else
             nv_=bl(nv)%BC(1) ! block imin of imax
          endif
       else
          nv_=bl(nv)%BC(5) ! block kmin of imax
       endif
    else
       nv_=nv
    endif
    !print *,'imax-kmin',nv,nv_,bl(n)%BC(5),bl(n)%BC(6)
    !print *,'proc',coord(1),coord(3)
    if ((nv_>0).and.((nv_==bl(n)%BC(5)).or.(nv_==bl(n)%BC(6))).and.((coord(1)==ndomx-1).and.(coord(3)==0))) then
       !print *,'do it'
       do j=ndyt,nfyt
          do k=1-ngh,0
             do i=1,ngh
                xc3(nx+i,j,k)=xc3(nx-k+1,j,i)
                yc3(nx+i,j,k)=yc3(nx-k+1,j,i)
                zc3(nx+i,j,k)=zc3(nx-k+1,j,i)
             enddo
          enddo
       enddo
    endif

    ! imax-kmax
    ! ---------
    nv=bl(n)%BC(2) ! block imax
    if (nv>0) then
       if (is_swapik2_bl(2)) then
          if (is_rev2(2,3)) then
             nv_=bl(nv)%BC(1) ! block imin of imax
          else
             nv_=bl(nv)%BC(2) ! block imax of imax
          endif
       else
          nv_=bl(nv)%BC(6) ! block kmax of imax
       endif
    else
       nv_=nv
    endif
    !print *,'jmax-kmax',nv,nv_,bl(n)%BC(5),bl(n)%BC(6)
    !print *,'proc',coord(1),coord(3)
    if ((nv_>0).and.((nv_==bl(n)%BC(5)).or.(nv_==bl(n)%BC(6))).and.((coord(1)==ndomx-1).and.(coord(3)==ndomz-1))) then
       !print *,'do it'
       do j=ndyt,nfyt
          do k=1,ngh
             do i=1,ngh
                xc3(nx+i,j,nz+k)=xc3(nx+k,j,nz-i+1)
                yc3(nx+i,j,nz+k)=yc3(nx+k,j,nz-i+1)
                zc3(nx+i,j,nz+k)=zc3(nx+k,j,nz-i+1)
             enddo
          enddo
       enddo
    endif

  end subroutine correct_edges_i

  !===============================================================================
  module subroutine correct_edges_j
  !===============================================================================
    !> reconstruct ghost cells of degenerate edges for j direction
  !===============================================================================
    use mod_block
    use mod_mpi_part
    use mod_grid_directions
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k
    integer :: n,nv,nv_
    ! ----------------------------------------------------------------------------

    ! current block number
    ! -------------------
    n=nob(iproc)

    ! Reconstruct ghost cells of degenerate edges only for j directions
    ! ------------------------------------------------------------------

    !print *,'================================'
    !print *,'BC',bl(n)%BC

    ndzt=1
    nfzt=nz
    if (bl(n)%BC(5)>0) ndzt=1-ngh
    if (bl(n)%BC(6)>0) nfzt=nz+ngh

    ! imin-jmin
    ! ---------
    nv=bl(n)%BC(3) ! block jmin
    if (nv>0) then
       if (is_swapij2_bl(3)) then
          if (is_rev2(3,1)) then
             nv_=bl(nv)%BC(4) ! block jmax of jmin
          else
             nv_=bl(nv)%BC(3) ! block jmin of jmin
          endif
       else
          nv_=bl(nv)%BC(1) ! block imin of jmin
       endif
    else
       nv_=nv
    endif
    !print *,'imin-jmin',nv,nv_,bl(n)%BC(1),bl(n)%BC(2)
    !print *,'proc',coord(1),coord(2)
    if ((nv_>0).and.((nv_==bl(n)%BC(1)).or.(nv_==bl(n)%BC(2))).and.((coord(1)==0).and.(coord(2)==0))) then
       !print *,'do it'
       do k=ndzt,nfzt
          do j=1-ngh,0
             do i=1-ngh,0
                xc3(i,j,k)=xc3(1-j,i,k)
                yc3(i,j,k)=yc3(1-j,i,k)
                zc3(i,j,k)=zc3(1-j,i,k)
             enddo
          enddo
       enddo
    endif

    ! imin-jmax
    ! ---------
    nv=bl(n)%BC(4) ! block jmax
    if (nv>0) then
       if (is_swapij2_bl(4)) then
          if (is_rev2(4,1)) then
             nv_=bl(nv)%BC(4) ! block jmax of jmax
          else
             nv_=bl(nv)%BC(3) ! block jmin of jmax
          endif
       else
          nv_=bl(nv)%BC(1) ! block imin of jmax
       endif
    else
       nv_=nv
    endif
    !print *,'imin-jmax',nv,nv_,bl(n)%BC(1),bl(n)%BC(2)
    !print *,'proc',coord(1),coord(2)
    if ((nv_>0).and.((nv_==bl(n)%BC(1)).or.(nv_==bl(n)%BC(2))).and.((coord(1)==0).and.(coord(2)==ndomy-1))) then
       !print *,'do it'
       do k=ndzt,nfzt
          do j=1,ngh
             do i=1-ngh,0
                xc3(i,ny+j,k)=xc3(j,ny-i+1,k)
                yc3(i,ny+j,k)=yc3(j,ny-i+1,k)
                zc3(i,ny+j,k)=zc3(j,ny-i+1,k)
             enddo
          enddo
       enddo
    endif

    ! imax-jmin
    ! ---------
    nv=bl(n)%BC(3) ! block jmin
    if (nv>0) then
       if (is_swapij2_bl(3)) then
          if (is_rev2(3,1)) then
             nv_=bl(nv)%BC(3) ! block jmin of jmin
          else
             nv_=bl(nv)%BC(4) ! block jmax of jmin
          endif
       else
          nv_=bl(nv)%BC(2) ! block imax of jmin
       endif
    else
       nv_=nv
    endif
    !print *,'imax-jmin',nv,nv_,bl(n)%BC(1),bl(n)%BC(2)
    !print *,'proc',coord(1),coord(2)
    if ((nv_>0).and.((nv_==bl(n)%BC(1)).or.(nv_==bl(n)%BC(2))).and.((coord(1)==ndomx-1).and.(coord(2)==0))) then
       !print *,'do it'
       do k=ndzt,nfzt
          do j=1-ngh,0
             do i=1,ngh
                xc3(nx+i,j,k)=xc3(nx+j,1-i,k)
                yc3(nx+i,j,k)=yc3(nx+j,1-i,k)
                zc3(nx+i,j,k)=zc3(nx+j,1-i,k)
             enddo
          enddo
       enddo
    endif

    ! imax-jmax
    ! ---------
    nv=bl(n)%BC(4) ! block jmax
    if (nv>0) then
       if (is_swapij2_bl(4)) then
          if (is_rev2(4,1)) then
             nv_=bl(nv)%BC(3) ! block jmin of jmax
          else
             nv_=bl(nv)%BC(4) ! block jmax of jmax
          endif
       else
          nv_=bl(nv)%BC(2) ! block imax of jmax
       endif
    else
       nv_=nv
    endif
    !print *,'imax-jmax',nv,nv_,bl(n)%BC(1),bl(n)%BC(2)
    !print *,'proc',coord(1),coord(2)
    if ((nv_>0).and.((nv_==bl(n)%BC(1)).or.(nv_==bl(n)%BC(2))).and.((coord(1)==ndomx-1).and.(coord(2)==ndomy-1))) then
       !print *,'do it'
       do k=ndzt,nfzt
          do j=1,ngh
             do i=1,ngh
                xc3(nx+i,ny+j,k)=xc3(nx-j+1,ny+i,k)
                yc3(nx+i,ny+j,k)=yc3(nx-j+1,ny+i,k)
                zc3(nx+i,ny+j,k)=zc3(nx-j+1,ny+i,k)
             enddo
          enddo
       enddo
    endif

    ndxt=1
    nfxt=nx
    if (bl(n)%BC(1)>0) ndxt=1-ngh
    if (bl(n)%BC(2)>0) nfxt=nx+ngh

    ! jmin-kmin
    ! ---------
    nv=bl(n)%BC(5) ! block kmin
    if (nv>0) then
       if (is_swapjk2_bl(5)) then
          if (is_rev2(5,3)) then
             nv_=bl(nv)%BC(6) ! block kmax of kmin
          else
             nv_=bl(nv)%BC(5) ! block kmin of kmin
          endif
       else
          nv_=bl(nv)%BC(3) ! block jmin of kmin
       endif
    else
       nv_=nv
    endif
    !print *,'jmin-kmin',nv,nv_,bl(n)%BC(3),bl(n)%BC(4)
    !print *,'proc',coord(2),coord(3)
    if ((nv_>0).and.((nv_==bl(n)%BC(3)).or.(nv_==bl(n)%BC(4))).and.((coord(2)==0).and.(coord(3)==0))) then
       !print *,'do it'
       do i=ndxt,nfxt
          do k=1-ngh,0
             do j=1-ngh,0
                xc3(i,j,k)=xc3(i,1-k,j)
                yc3(i,j,k)=yc3(i,1-k,j)
                zc3(i,j,k)=zc3(i,1-k,j)
             enddo
          enddo
       enddo
    endif

    ! jmin-kmax
    ! ---------
    nv=bl(n)%BC(6) ! block kmax
    if (nv>0) then
       if (is_swapjk2_bl(6)) then
          if (is_rev2(6,3)) then
             nv_=bl(nv)%BC(6) ! block kmax of kmax
          else
             nv_=bl(nv)%BC(5) ! block kmin of kmax
          endif
       else
          nv_=bl(nv)%BC(3) ! block jmin of kmax
       endif
       !print *,'jmin-kmax',nv,nv_,bl(n)%BC(3),bl(n)%BC(4)
       !print *,'proc',coord(2),coord(3)
    else
       nv_=nv
    endif
    if ((nv_>0).and.((nv_==bl(n)%BC(3)).or.(nv_==bl(n)%BC(4))).and.((coord(2)==0).and.(coord(3)==ndomz-1))) then
       !print *,'do it'
       do i=ndxt,nfxt
          do k=1,ngh
             do j=1-ngh,0
                xc3(i,j,nz+k)=xc3(i,k,nz-j+1)
                yc3(i,j,nz+k)=yc3(i,k,nz-j+1)
                zc3(i,j,nz+k)=zc3(i,k,nz-j+1)
             enddo
          enddo
       enddo
    endif

    ! jmax-kmin
    ! ---------
    nv=bl(n)%BC(5) ! block kmin
    if (nv>0) then
       if (is_swapjk2_bl(5)) then
          if (is_rev2(5,3)) then
             nv_=bl(nv)%BC(5) ! block kmin of kmin
          else
             nv_=bl(nv)%BC(6) ! block kmax of kmin
          endif
       else
          nv_=bl(nv)%BC(4) ! block jmax of kmin
       endif
    else
       nv_=nv
    endif
    !print *,'jmax-kmin',nv,nv_,bl(n)%BC(3),bl(n)%BC(4)
    !print *,'proc',coord(2),coord(3)
    if ((nv_>0).and.((nv_==bl(n)%BC(3)).or.(nv_==bl(n)%BC(4))).and.((coord(2)==ndomy-1).and.(coord(3)==0))) then
       !print *,'do it'
       do i=ndxt,nfxt
          do k=1-ngh,0
             do j=1,ngh
                xc3(i,ny+j,k)=xc3(i,ny+k,1-j)
                yc3(i,ny+j,k)=yc3(i,ny+k,1-j)
                zc3(i,ny+j,k)=zc3(i,ny+k,1-j)
             enddo
          enddo
       enddo
    endif

    ! jmax-kmax
    ! ---------
    nv=bl(n)%BC(6) ! block kmax
    if (nv>0) then
       if (is_swapjk2_bl(6)) then
          if (is_rev2(6,3)) then
             nv_=bl(nv)%BC(5) ! block kmin of kmax
          else
             nv_=bl(nv)%BC(6) ! block kmax of kmax
          endif
       else
          nv_=bl(nv)%BC(4) ! block jmax of kmax
       endif
    else
       nv_=nv
    endif
    !print *,'jmax-kmax',nv,nv_,bl(n)%BC(3),bl(n)%BC(4)
    !print *,'proc',coord(2),coord(3)
    if ((nv_>0).and.((nv_==bl(n)%BC(3)).or.(nv_==bl(n)%BC(4))).and.((coord(2)==ndomy-1).and.(coord(3)==ndomz-1))) then
       !print *,'do it'
       do i=ndxt,nfxt
          do k=1,ngh
             do j=1,ngh
                xc3(i,ny+j,nz+k)=xc3(i,ny-k+1,nz+j)
                yc3(i,ny+j,nz+k)=yc3(i,ny-k+1,nz+j)
                zc3(i,ny+j,nz+k)=zc3(i,ny-k+1,nz+j)
             enddo
          enddo
       enddo
    endif

    ndyt=1
    nfyt=ny
    if (bl(n)%BC(3)>0) ndyt=1-ngh
    if (bl(n)%BC(4)>0) nfyt=ny+ngh

    ! imin-kmin
    ! ---------
    nv=bl(n)%BC(5) ! block kmin
    if (nv>0) then
       if (is_swapik2_bl(5)) then
          if (is_rev2(5,1)) then
             nv_=bl(nv)%BC(6) ! block kmax of kmin
          else
             nv_=bl(nv)%BC(5) ! block kmin of kmin
          endif
       else
          nv_=bl(nv)%BC(1) ! block imin of kmin
       endif
    else
       nv_=nv
    endif
    !print *,'imin-kmin',nv,nv_,bl(n)%BC(1),bl(n)%BC(2)
    !print *,'proc',coord(1),coord(3)
    if ((nv_>0).and.((nv_==bl(n)%BC(1)).or.(nv_==bl(n)%BC(2))).and.((coord(1)==0).and.(coord(3)==0))) then
       !print *,'do it'
       do j=ndyt,nfyt
          do k=1-ngh,0
             do i=1-ngh,0
                xc3(i,j,k)=xc3(1-k,j,i)
                yc3(i,j,k)=yc3(1-k,j,i)
                zc3(i,j,k)=zc3(1-k,j,i)
             enddo
          enddo
       enddo
    endif

    ! imin-kmax
    ! ---------
    nv=bl(n)%BC(6) ! block kmax
    if (nv>0) then
       if (is_swapik2_bl(6)) then
          if (is_rev2(6,1)) then
             nv_=bl(nv)%BC(6) ! block kmax of kmax
          else
             nv_=bl(nv)%BC(5) ! block kmin of kmax
          endif
       else
          nv_=bl(nv)%BC(1) ! block imin of kmax
       endif
    else
       nv_=nv
    endif
    !print *,'imin-kmax',nv,nv_,bl(n)%BC(1),bl(n)%BC(2)
    !print *,'proc',coord(1),coord(3)
    if ((nv_>0).and.((nv_==bl(n)%BC(1)).or.(nv_==bl(n)%BC(2))).and.((coord(1)==0).and.(coord(3)==ndomz-1))) then
       !print *,'do it'
       do j=ndyt,nfyt
          do k=1,ngh
             do i=1-ngh,0
                xc3(i,j,nz+k)=xc3(k,j,nz-i+1)
                yc3(i,j,nz+k)=yc3(k,j,nz-i+1)
                zc3(i,j,nz+k)=zc3(k,j,nz-i+1)
             enddo
          enddo
       enddo
    endif

    ! imax-kmin
    ! ---------
    nv=bl(n)%BC(5) ! block kmin
    if (nv>0) then
       if (is_swapik2_bl(5)) then
          if (is_rev2(5,1)) then
             nv_=bl(nv)%BC(5) ! block kmin of kmin
          else
             nv_=bl(nv)%BC(6) ! block kmax of kmin
          endif
       else
          nv_=bl(nv)%BC(2) ! block imax of kmin
       endif
    else
       nv_=nv
    endif
    !print *,'imax-kmin',nv,nv_,bl(n)%BC(1),bl(n)%BC(2)
    !print *,'proc',coord(1),coord(3)
    if ((nv_>0).and.((nv_==bl(n)%BC(1)).or.(nv_==bl(n)%BC(2))).and.((coord(1)==ndomx-1).and.(coord(3)==0))) then
       !print *,'do it'
       do j=ndyt,nfyt
          do k=1-ngh,0
             do i=1,ngh
                xc3(nx+i,j,k)=xc3(nx+k,j,1-i)
                yc3(nx+i,j,k)=yc3(nx+k,j,1-i)
                zc3(nx+i,j,k)=zc3(nx+k,j,1-i)
             enddo
          enddo
       enddo
    endif

    ! imax-kmax
    ! ---------
    nv=bl(n)%BC(6) ! block kmax
    if (nv>0) then
       if (is_swapik2_bl(6)) then
          if (is_rev2(6,1)) then
             nv_=bl(nv)%BC(5) ! block kmin of kmax
          else
             nv_=bl(nv)%BC(6) ! block kmax of kmax
          endif
       else
          nv_=bl(nv)%BC(2) ! block imax of kmax
       endif
    else
       nv_=nv
    endif
    !print *,'jmax-kmax',nv,nv_,bl(n)%BC(1),bl(n)%BC(2)
    !print *,'proc',coord(1),coord(3)
    if ((nv_>0).and.((nv_==bl(n)%BC(1)).or.(nv_==bl(n)%BC(2))).and.((coord(1)==ndomx-1).and.(coord(3)==ndomz-1))) then
       !print *,'do it'
       do j=ndyt,nfyt
          do k=1,ngh
             do i=1,ngh
                xc3(nx+i,j,nz+k)=xc3(nx-k+1,j,nz+i)
                yc3(nx+i,j,nz+k)=yc3(nx-k+1,j,nz+i)
                zc3(nx+i,j,nz+k)=zc3(nx-k+1,j,nz+i)
             enddo
          enddo
       enddo
    endif

  end subroutine correct_edges_j
  
end submodule smod_grid_metrics2_c3
