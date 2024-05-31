!=================================================================================
submodule (mod_grid_metrics_c3) smod_grid_metrics_deriv_phi_c3
!=================================================================================
  !> Module to compute full 3D curvilinear metrics
  !> with Geometric Conservation Law (GCL)
!=================================================================================

contains

  !===============================================================================
  module subroutine derivative_phi(var1,dvar1,var2,dvar2,var3,dvar3)
  !===============================================================================
    !> Compute derivatives along phi for metrics (full 3D version)
    !> - stencil -ngh:+ngh for inviscid fluxes + BCs -
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(in) :: var1,var2,var3
    real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(inout) :: dvar1,dvar2,dvar3
    ! ---------------------------------------------------------------------------
    integer :: i,j,k,l
    real(wp), dimension(-ngh:ngh) :: a_ ! coefficients of interior FD scheme
    ! ---------------------------------------------------------------------------

    ! Metrics for BC at kmin
    ! ======================
    
    ! For wall & charac, inlet and outlet BCs
    ! ---------------------------------------
    if ((BC_face(3,1)%sort==0).or.(BC_face(3,1)%sort<=-3)) then

       ! SBP boundary schemes
       ! --------------------
       if (is_SBP) then
          if (ngh<4) then
             call mpistop('3D SBP metrics need yet at least 9pt-stencil !', 0)
          else
             ! point #1: stencil [o x x x]
             !k=1
             do j=ndyt,nfyt
                do i=ndxt,nfxt
                   dvar1(i,j,1)= as4p0(1)*var1(i,j,1)+as4p0(2)*var1(i,j,2) &
                               + as4p0(3)*var1(i,j,3)+as4p0(4)*var1(i,j,4) &
                               + dvar1(i,j,1)
                   dvar2(i,j,1)= as4p0(1)*var2(i,j,1)+as4p0(2)*var2(i,j,2) &
                               + as4p0(3)*var2(i,j,3)+as4p0(4)*var2(i,j,4) &
                               + dvar2(i,j,1)
                   dvar3(i,j,1)= as4p0(1)*var3(i,j,1)+as4p0(2)*var3(i,j,2) &
                               + as4p0(3)*var3(i,j,3)+as4p0(4)*var3(i,j,4) &
                               + dvar3(i,j,1)
                enddo
             enddo
             ! point #2: stencil [x o x x]
             !k=2
             do j=ndyt,nfyt
                do i=ndxt,nfxt
                   dvar1(i,j,2)= as4p1(1)*var1(i,j,1)+as4p1(2)*var1(i,j,2) &
                               + as4p1(3)*var1(i,j,3)+as4p1(4)*var1(i,j,4) &
                               + dvar1(i,j,2)
                   dvar2(i,j,2)= as4p1(1)*var2(i,j,1)+as4p1(2)*var2(i,j,2) &
                               + as4p1(3)*var2(i,j,3)+as4p1(4)*var2(i,j,4) &
                               + dvar2(i,j,2)
                   dvar3(i,j,2)= as4p1(1)*var3(i,j,1)+as4p1(2)*var3(i,j,2) &
                               + as4p1(3)*var3(i,j,3)+as4p1(4)*var3(i,j,4) &
                               + dvar3(i,j,2)
                enddo
             enddo
             ! point #3: stencil [x x o x x]
             !k=3
             do j=ndyt,nfyt
                do i=ndxt,nfxt
                   dvar1(i,j,3)= as4p2(1)*var1(i,j,1)+as4p2(2)*var1(i,j,2) &
                               + as4p2(3)*var1(i,j,3)+as4p2(4)*var1(i,j,4) &
                               + as4p2(5)*var1(i,j,5) + dvar1(i,j,3)
                   dvar2(i,j,3)= as4p2(1)*var2(i,j,1)+as4p2(2)*var2(i,j,2) &
                               + as4p2(3)*var2(i,j,3)+as4p2(4)*var2(i,j,4) &
                               + as4p2(5)*var2(i,j,5) + dvar2(i,j,3)
                   dvar3(i,j,3)= as4p2(1)*var3(i,j,1)+as4p2(2)*var3(i,j,2) &
                               + as4p2(3)*var3(i,j,3)+as4p2(4)*var3(i,j,4) &
                               + as4p2(5)*var3(i,j,5) + dvar3(i,j,3)
                enddo
             enddo
             ! point #4: stencil [x x x o x x]
             !k=4
             do j=ndyt,nfyt
                do i=ndxt,nfxt
                   dvar1(i,j,4)= as4p3(1)*var1(i,j,1)+as4p3(2)*var1(i,j,2) &
                               + as4p3(3)*var1(i,j,3)+as4p3(4)*var1(i,j,4) &
                               + as4p3(5)*var1(i,j,5)+as4p3(6)*var1(i,j,6) &
                               + dvar1(i,j,4)
                   dvar2(i,j,4)= as4p3(1)*var2(i,j,1)+as4p3(2)*var2(i,j,2) &
                               + as4p3(3)*var2(i,j,3)+as4p3(4)*var2(i,j,4) &
                               + as4p3(5)*var2(i,j,5)+as4p3(6)*var2(i,j,6) &
                               + dvar2(i,j,4)
                   dvar3(i,j,4)= as4p3(1)*var3(i,j,1)+as4p3(2)*var3(i,j,2) &
                               + as4p3(3)*var3(i,j,3)+as4p3(4)*var3(i,j,4) &
                               + as4p3(5)*var3(i,j,5)+as4p3(6)*var3(i,j,6) &
                               + dvar3(i,j,4)
                enddo
             enddo
          endif
          if (ngh>=5) then
             ! eighth-order (9-pt stencil) for point k=5
             k=5
             do j=ndyt,nfyt
                do i=ndxt,nfxt
                   dvar1(i,j,k)= a9(1)*(var1(i,j,k+1)-var1(i,j,k-1)) &
                               + a9(2)*(var1(i,j,k+2)-var1(i,j,k-2)) &
                               + a9(3)*(var1(i,j,k+3)-var1(i,j,k-3)) &
                               + a9(4)*(var1(i,j,k+4)-var1(i,j,k-4)) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= a9(1)*(var2(i,j,k+1)-var2(i,j,k-1)) &
                               + a9(2)*(var2(i,j,k+2)-var2(i,j,k-2)) &
                               + a9(3)*(var2(i,j,k+3)-var2(i,j,k-3)) &
                               + a9(4)*(var2(i,j,k+4)-var2(i,j,k-4)) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= a9(1)*(var3(i,j,k+1)-var3(i,j,k-1)) &
                               + a9(2)*(var3(i,j,k+2)-var3(i,j,k-2)) &
                               + a9(3)*(var3(i,j,k+3)-var3(i,j,k-3)) &
                               + a9(4)*(var3(i,j,k+4)-var3(i,j,k-4)) &
                               + dvar3(i,j,k)
                enddo
             enddo
          endif

       ! Reduced order stencil
       ! ---------------------
       else       
          if (ngh>=1) then
             ! first-order (2-pt stencil) for point k=1
             k=1
             do j=ndyt,nfyt
                do i=ndxt,nfxt
                   !do l=0,1
                   !   dvar1(i,j,k)=dvar1(i,j,k)+a01(1+l)*var1(i,j,k+l)
                   !   dvar2(i,j,k)=dvar2(i,j,k)+a01(1+l)*var2(i,j,k+l)
                   !   dvar3(i,j,k)=dvar3(i,j,k)+a01(1+l)*var3(i,j,k+l)
                   !enddo
                   do l=0,2
                      dvar1(i,j,k)=dvar1(i,j,k)+a02(1+l)*var1(i,j,k+l)
                      dvar2(i,j,k)=dvar2(i,j,k)+a02(1+l)*var2(i,j,k+l)
                      dvar3(i,j,k)=dvar3(i,j,k)+a02(1+l)*var3(i,j,k+l)
                   enddo
                enddo
             enddo
          endif
          if (ngh>=2) then
             ! second-order (3-pt stencil) for point k=2
             k=2
             do j=ndyt,nfyt
                do i=ndxt,nfxt
                   dvar1(i,j,k)= a3(1)*(var1(i,j,k+1)-var1(i,j,k-1)) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= a3(1)*(var2(i,j,k+1)-var2(i,j,k-1)) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= a3(1)*(var3(i,j,k+1)-var3(i,j,k-1)) &
                               + dvar3(i,j,k)
                enddo
             enddo
          endif
          if (ngh>=3) then
             ! fourth-order (5-pt stencil) for point k=3
             k=3
             do j=ndyt,nfyt
                do i=ndxt,nfxt
                   dvar1(i,j,k)= a5(1)*(var1(i,j,k+1)-var1(i,j,k-1)) &
                               + a5(2)*(var1(i,j,k+2)-var1(i,j,k-2)) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= a5(1)*(var2(i,j,k+1)-var2(i,j,k-1)) &
                               + a5(2)*(var2(i,j,k+2)-var2(i,j,k-2)) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= a5(1)*(var3(i,j,k+1)-var3(i,j,k-1)) &
                               + a5(2)*(var3(i,j,k+2)-var3(i,j,k-2)) &
                               + dvar3(i,j,k)
                enddo
             enddo
          endif
          if (ngh>=4) then
             ! sixth-order (7-pt stencil) for point k=4
             k=4
             do j=ndyt,nfyt
                do i=ndxt,nfxt
                   dvar1(i,j,k)= a7(1)*(var1(i,j,k+1)-var1(i,j,k-1)) &
                               + a7(2)*(var1(i,j,k+2)-var1(i,j,k-2)) &
                               + a7(3)*(var1(i,j,k+3)-var1(i,j,k-3)) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= a7(1)*(var2(i,j,k+1)-var2(i,j,k-1)) &
                               + a7(2)*(var2(i,j,k+2)-var2(i,j,k-2)) &
                               + a7(3)*(var2(i,j,k+3)-var2(i,j,k-3)) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= a7(1)*(var3(i,j,k+1)-var3(i,j,k-1)) &
                               + a7(2)*(var3(i,j,k+2)-var3(i,j,k-2)) &
                               + a7(3)*(var3(i,j,k+3)-var3(i,j,k-3)) &
                               + dvar3(i,j,k)
                enddo
             enddo
          endif
          if (ngh==5) then
             ! eighth-order (9-pt stencil) for point k=5
             k=5
             do j=ndyt,nfyt
                do i=ndxt,nfxt
                   dvar1(i,j,k)= a9(1)*(var1(i,j,k+1)-var1(i,j,k-1)) &
                               + a9(2)*(var1(i,j,k+2)-var1(i,j,k-2)) &
                               + a9(3)*(var1(i,j,k+3)-var1(i,j,k-3)) &
                               + a9(4)*(var1(i,j,k+4)-var1(i,j,k-4)) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= a9(1)*(var2(i,j,k+1)-var2(i,j,k-1)) &
                               + a9(2)*(var2(i,j,k+2)-var2(i,j,k-2)) &
                               + a9(3)*(var2(i,j,k+3)-var2(i,j,k-3)) &
                               + a9(4)*(var2(i,j,k+4)-var2(i,j,k-4)) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= a9(1)*(var3(i,j,k+1)-var3(i,j,k-1)) &
                               + a9(2)*(var3(i,j,k+2)-var3(i,j,k-2)) &
                               + a9(3)*(var3(i,j,k+3)-var3(i,j,k-3)) &
                               + a9(4)*(var3(i,j,k+4)-var3(i,j,k-4)) &
                               + dvar3(i,j,k)
                enddo
             enddo
          endif
       endif
    endif

    ! 7-point stencil for Tam&Dong BC
    ! -------------------------------
    if ((BC_face(3,1)%sort==-1).or.(BC_face(3,1)%sort==-2)) then
       if (ngh<3) then
          call mpistop('Tam&Dong BCs need at least 7pt-stencil !', 0)
       else
          k=1
          do j=ndyt,nfyt
             do i=ndxt,nfxt
                do l=0,6
                   dvar1(i,j,k)=dvar1(i,j,k)+a06(1+l)*var1(i,j,k+l)
                   dvar2(i,j,k)=dvar2(i,j,k)+a06(1+l)*var2(i,j,k+l)
                   dvar3(i,j,k)=dvar3(i,j,k)+a06(1+l)*var3(i,j,k+l)
                enddo
             enddo
          enddo
          k=2
          do j=ndyt,nfyt
             do i=ndxt,nfxt
                do l=-1,5
                   dvar1(i,j,k)=dvar1(i,j,k)+a15(2+l)*var1(i,j,k+l)
                   dvar2(i,j,k)=dvar2(i,j,k)+a15(2+l)*var2(i,j,k+l)
                   dvar3(i,j,k)=dvar3(i,j,k)+a15(2+l)*var3(i,j,k+l)
                enddo
             enddo
          enddo
          k=3
          do j=ndyt,nfyt
             do i=ndxt,nfxt
                do l=-2,4
                   dvar1(i,j,k)=dvar1(i,j,k)+a24(3+l)*var1(i,j,k+l)
                   dvar2(i,j,k)=dvar2(i,j,k)+a24(3+l)*var2(i,j,k+l)
                   dvar3(i,j,k)=dvar3(i,j,k)+a24(3+l)*var3(i,j,k+l)
                enddo
             enddo
          enddo
       endif
       if (ngh>=4) then
          ! sixth-order (7-pt stencil) for point k=4
          k=4
          do j=ndyt,nfyt
             do i=ndxt,nfxt
                dvar1(i,j,k)= a7(1)*(var1(i,j,k+1)-var1(i,j,k-1)) &
                            + a7(2)*(var1(i,j,k+2)-var1(i,j,k-2)) &
                            + a7(3)*(var1(i,j,k+3)-var1(i,j,k-3)) &
                            + dvar1(i,j,k)
                dvar2(i,j,k)= a7(1)*(var2(i,j,k+1)-var2(i,j,k-1)) &
                            + a7(2)*(var2(i,j,k+2)-var2(i,j,k-2)) &
                            + a7(3)*(var2(i,j,k+3)-var2(i,j,k-3)) &
                            + dvar2(i,j,k)
                dvar3(i,j,k)= a7(1)*(var3(i,j,k+1)-var3(i,j,k-1)) &
                            + a7(2)*(var3(i,j,k+2)-var3(i,j,k-2)) &
                            + a7(3)*(var3(i,j,k+3)-var3(i,j,k-3)) &
                            + dvar3(i,j,k)
             enddo
          enddo
       endif
       if (ngh==5) then
          ! sixth-order (7-pt stencil) for point k=5
          k=5
          do j=ndyt,nfyt
             do i=ndxt,nfxt
                dvar1(i,j,k)= a7(1)*(var1(i,j,k+1)-var1(i,j,k-1)) &
                            + a7(2)*(var1(i,j,k+2)-var1(i,j,k-2)) &
                            + a7(3)*(var1(i,j,k+3)-var1(i,j,k-3)) &
                            + dvar1(i,j,k)
                dvar2(i,j,k)= a7(1)*(var2(i,j,k+1)-var2(i,j,k-1)) &
                            + a7(2)*(var2(i,j,k+2)-var2(i,j,k-2)) &
                            + a7(3)*(var2(i,j,k+3)-var2(i,j,k-3)) &
                            + dvar2(i,j,k)
                dvar3(i,j,k)= a7(1)*(var3(i,j,k+1)-var3(i,j,k-1)) &
                            + a7(2)*(var3(i,j,k+2)-var3(i,j,k-2)) &
                            + a7(3)*(var3(i,j,k+3)-var3(i,j,k-3)) &
                            + dvar3(i,j,k)
             enddo
          enddo
       endif
    endif
          
    ! Metrics for BC at kmax
    ! ======================
    
    ! For wall & charac, inlet and outlet BCs
    ! ---------------------------------------
    if ((BC_face(3,2)%sort==0).or.(BC_face(3,2)%sort<=-3)) then
       
       ! SBP boundary schemes
       ! --------------------
       if (is_SBP) then
          if (ngh<4) then
             call mpistop('3D SBP metrics need yet at least 9pt-stencil !', 0)
          else
             ! point #1: stencil [o x x x]
             k=nz
             do j=ndyt,nfyt
                do i=ndxt,nfxt
                   dvar1(i,j,k)= as4m0(1)*var1(i,j,nz  )+as4m0(2)*var1(i,j,nz-1) &
                               + as4m0(3)*var1(i,j,nz-2)+as4m0(4)*var1(i,j,nz-3) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= as4m0(1)*var2(i,j,nz  )+as4m0(2)*var2(i,j,nz-1) &
                               + as4m0(3)*var2(i,j,nz-2)+as4m0(4)*var2(i,j,nz-3) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= as4m0(1)*var3(i,j,nz  )+as4m0(2)*var3(i,j,nz-1) &
                               + as4m0(3)*var3(i,j,nz-2)+as4m0(4)*var3(i,j,nz-3) &
                               + dvar3(i,j,k)
                enddo
             enddo
             ! point #2: stencil [x o x x]
             k=nz-1
             do j=ndyt,nfyt
                do i=ndxt,nfxt
                   dvar1(i,j,k)= as4m1(1)*var1(i,j,nz  )+as4m1(2)*var1(i,j,nz-1) &
                               + as4m1(3)*var1(i,j,nz-2)+as4m1(4)*var1(i,j,nz-3) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= as4m1(1)*var2(i,j,nz  )+as4m1(2)*var2(i,j,nz-1) &
                               + as4m1(3)*var2(i,j,nz-2)+as4m1(4)*var2(i,j,nz-3) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= as4m1(1)*var3(i,j,nz  )+as4m1(2)*var3(i,j,nz-1) &
                               + as4m1(3)*var3(i,j,nz-2)+as4m1(4)*var3(i,j,nz-3) &
                               + dvar3(i,j,k)
                enddo
             enddo
             ! point #3: stencil [x x o x x]
             k=nz-2
             do j=ndyt,nfyt
                do i=ndxt,nfxt
                   dvar1(i,j,k)= as4m2(1)*var1(i,j,nz  )+as4m2(2)*var1(i,j,nz-1) &
                               + as4m2(3)*var1(i,j,nz-2)+as4m2(4)*var1(i,j,nz-3) &
                               + as4m2(5)*var1(i,j,nz-4) + dvar1(i,j,k)
                   dvar2(i,j,k)= as4m2(1)*var2(i,j,nz  )+as4m2(2)*var2(i,j,nz-1) &
                               + as4m2(3)*var2(i,j,nz-2)+as4m2(4)*var2(i,j,nz-3) &
                               + as4m2(5)*var2(i,j,nz-4) + dvar2(i,j,k)
                   dvar3(i,j,k)= as4m2(1)*var3(i,j,nz  )+as4m2(2)*var3(i,j,nz-1) &
                               + as4m2(3)*var3(i,j,nz-2)+as4m2(4)*var3(i,j,nz-3) &
                               + as4m2(5)*var3(i,j,nz-4) + dvar3(i,j,k)
                enddo
             enddo
             ! point #4: stencil [x x x o x x]
             k=nz-3
             do j=ndyt,nfyt
                do i=ndxt,nfxt
                   dvar1(i,j,k)= as4m3(1)*var1(i,j,nz  )+as4m3(2)*var1(i,j,nz-1) &
                               + as4m3(3)*var1(i,j,nz-2)+as4m3(4)*var1(i,j,nz-3) &
                               + as4m3(5)*var1(i,j,nz-4)+as4m3(6)*var1(i,j,nz-5) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= as4m3(1)*var2(i,j,nz  )+as4m3(2)*var2(i,j,nz-1) &
                               + as4m3(3)*var2(i,j,nz-2)+as4m3(4)*var2(i,j,nz-3) &
                               + as4m3(5)*var2(i,j,nz-4)+as4m3(6)*var2(i,j,nz-5) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= as4m3(1)*var3(i,j,nz  )+as4m3(2)*var3(i,j,nz-1) &
                               + as4m3(3)*var3(i,j,nz-2)+as4m3(4)*var3(i,j,nz-3) &
                               + as4m3(5)*var3(i,j,nz-4)+as4m3(6)*var3(i,j,nz-5) &
                               + dvar3(i,j,k)
                enddo
             enddo
          endif
          if (ngh>=5) then
             ! eighth-order (9-pt stencil) for point k=nz-4
             k=nz-4
             do j=ndyt,nfyt
                do i=ndxt,nfxt
                   dvar1(i,j,k)= a9(1)*(var1(i,j,k+1)-var1(i,j,k-1)) &
                               + a9(2)*(var1(i,j,k+2)-var1(i,j,k-2)) &
                               + a9(3)*(var1(i,j,k+3)-var1(i,j,k-3)) &
                               + a9(4)*(var1(i,j,k+4)-var1(i,j,k-4)) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= a9(1)*(var2(i,j,k+1)-var2(i,j,k-1)) &
                               + a9(2)*(var2(i,j,k+2)-var2(i,j,k-2)) &
                               + a9(3)*(var2(i,j,k+3)-var2(i,j,k-3)) &
                               + a9(4)*(var2(i,j,k+4)-var2(i,j,k-4)) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= a9(1)*(var3(i,j,k+1)-var3(i,j,k-1)) &
                               + a9(2)*(var3(i,j,k+2)-var3(i,j,k-2)) &
                               + a9(3)*(var3(i,j,k+3)-var3(i,j,k-3)) &
                               + a9(4)*(var3(i,j,k+4)-var3(i,j,k-4)) &
                               + dvar3(i,j,k)
                enddo
             enddo
          endif
          
       ! Reduced order stencil
       ! ---------------------
       else       
          if (ngh>=1) then
             ! first-order (2-pt stencil) for point k=nz
             k=nz
             do j=ndyt,nfyt
                do i=ndxt,nfxt
                   !do l=-1,0
                   !   dvar1(i,j,k)=dvar1(i,j,k)-a01(1-l)*var1(i,j,k+l)
                   !   dvar2(i,j,k)=dvar2(i,j,k)-a01(1-l)*var2(i,j,k+l)
                   !   dvar3(i,j,k)=dvar3(i,j,k)-a01(1-l)*var3(i,j,k+l)
                   !enddo
                   do l=-2,0
                      dvar1(i,j,k)=dvar1(i,j,k)-a02(1-l)*var1(i,j,k+l)
                      dvar2(i,j,k)=dvar2(i,j,k)-a02(1-l)*var2(i,j,k+l)
                      dvar3(i,j,k)=dvar3(i,j,k)-a02(1-l)*var3(i,j,k+l)
                   enddo
                enddo
             enddo
          endif
          if (ngh>=2) then
             ! second-order (3-pt stencil) for point k=nz-1
             k=nz-1
             do j=ndyt,nfyt
                do i=ndxt,nfxt
                   dvar1(i,j,k)= a3(1)*(var1(i,j,k+1)-var1(i,j,k-1)) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= a3(1)*(var2(i,j,k+1)-var2(i,j,k-1)) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= a3(1)*(var3(i,j,k+1)-var3(i,j,k-1)) &
                               + dvar3(i,j,k)
                enddo
             enddo
          endif
          if (ngh>=3) then
             ! fourth-order (5-pt stencil) for point k=nz-2
             k=nz-2
             do j=ndyt,nfyt
                do i=ndxt,nfxt
                   dvar1(i,j,k)= a5(1)*(var1(i,j,k+1)-var1(i,j,k-1)) &
                               + a5(2)*(var1(i,j,k+2)-var1(i,j,k-2)) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= a5(1)*(var2(i,j,k+1)-var2(i,j,k-1)) &
                               + a5(2)*(var2(i,j,k+2)-var2(i,j,k-2)) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= a5(1)*(var3(i,j,k+1)-var3(i,j,k-1)) &
                               + a5(2)*(var3(i,j,k+2)-var3(i,j,k-2)) &
                               + dvar3(i,j,k)
                enddo
             enddo
          endif
          if (ngh>=4) then
             ! sixth-order (7-pt stencil) for point k=nz-3
             k=nz-3
             do j=ndyt,nfyt
                do i=ndxt,nfxt
                   dvar1(i,j,k)= a7(1)*(var1(i,j,k+1)-var1(i,j,k-1)) &
                               + a7(2)*(var1(i,j,k+2)-var1(i,j,k-2)) &
                               + a7(3)*(var1(i,j,k+3)-var1(i,j,k-3)) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= a7(1)*(var2(i,j,k+1)-var2(i,j,k-1)) &
                               + a7(2)*(var2(i,j,k+2)-var2(i,j,k-2)) &
                               + a7(3)*(var2(i,j,k+3)-var2(i,j,k-3)) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= a7(1)*(var3(i,j,k+1)-var3(i,j,k-1)) &
                               + a7(2)*(var3(i,j,k+2)-var3(i,j,k-2)) &
                               + a7(3)*(var3(i,j,k+3)-var3(i,j,k-3)) &
                               + dvar3(i,j,k)
                enddo
             enddo
          endif
          if (ngh==5) then
             ! eighth-order (9-pt stencil) for point k=nz-4
             k=nz-4
             do j=ndyt,nfyt
                do i=ndxt,nfxt
                   dvar1(i,j,k)= a9(1)*(var1(i,j,k+1)-var1(i,j,k-1)) &
                               + a9(2)*(var1(i,j,k+2)-var1(i,j,k-2)) &
                               + a9(3)*(var1(i,j,k+3)-var1(i,j,k-3)) &
                               + a9(4)*(var1(i,j,k+4)-var1(i,j,k-4)) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= a9(1)*(var2(i,j,k+1)-var2(i,j,k-1)) &
                               + a9(2)*(var2(i,j,k+2)-var2(i,j,k-2)) &
                               + a9(3)*(var2(i,j,k+3)-var2(i,j,k-3)) &
                               + a9(4)*(var2(i,j,k+4)-var2(i,j,k-4)) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= a9(1)*(var3(i,j,k+1)-var3(i,j,k-1)) &
                               + a9(2)*(var3(i,j,k+2)-var3(i,j,k-2)) &
                               + a9(3)*(var3(i,j,k+3)-var3(i,j,k-3)) &
                               + a9(4)*(var3(i,j,k+4)-var3(i,j,k-4)) &
                               + dvar3(i,j,k)
                enddo
             enddo
          endif
       endif
    endif

    ! 7-point stencil for Tam&Dong BC
    ! -------------------------------
    if ((BC_face(3,2)%sort==-1).or.(BC_face(3,2)%sort==-2)) then
       if (ngh<3) then
          call mpistop('Tam&Dong BCs need at least 7pt-stencil !', 0)
       else
          k=nz
          do j=ndyt,nfyt
             do i=ndxt,nfxt
                do l=-6,0
                   dvar1(i,j,k)=dvar1(i,j,k)-a06(1-l)*var1(i,j,k+l)
                   dvar2(i,j,k)=dvar2(i,j,k)-a06(1-l)*var2(i,j,k+l)
                   dvar3(i,j,k)=dvar3(i,j,k)-a06(1-l)*var3(i,j,k+l)
                enddo
             enddo
          enddo
          k=nz-1
          do j=ndyt,nfyt
             do i=ndxt,nfxt
                do l=-5,1
                   dvar1(i,j,k)=dvar1(i,j,k)-a15(2-l)*var1(i,j,k+l)
                   dvar2(i,j,k)=dvar2(i,j,k)-a15(2-l)*var2(i,j,k+l)
                   dvar3(i,j,k)=dvar3(i,j,k)-a15(2-l)*var3(i,j,k+l)
                enddo
             enddo
          enddo
          k=nz-2
          do j=ndyt,nfyt
             do i=ndxt,nfxt
                do l=-4,2
                   dvar1(i,j,k)=dvar1(i,j,k)-a24(3-l)*var1(i,j,k+l)
                   dvar2(i,j,k)=dvar2(i,j,k)-a24(3-l)*var2(i,j,k+l)
                   dvar3(i,j,k)=dvar3(i,j,k)-a24(3-l)*var3(i,j,k+l)
                enddo
             enddo
          enddo
       endif
       if (ngh>=4) then
          ! sixth-order (7-pt stencil) for point k=nz-3
          k=nz-3
          do j=ndyt,nfyt
             do i=ndxt,nfxt
                dvar1(i,j,k)= a7(1)*(var1(i,j,k+1)-var1(i,j,k-1)) &
                            + a7(2)*(var1(i,j,k+2)-var1(i,j,k-2)) &
                            + a7(3)*(var1(i,j,k+3)-var1(i,j,k-3)) &
                            + dvar1(i,j,k)
                dvar2(i,j,k)= a7(1)*(var2(i,j,k+1)-var2(i,j,k-1)) &
                            + a7(2)*(var2(i,j,k+2)-var2(i,j,k-2)) &
                            + a7(3)*(var2(i,j,k+3)-var2(i,j,k-3)) &
                            + dvar2(i,j,k)
                dvar3(i,j,k)= a7(1)*(var3(i,j,k+1)-var3(i,j,k-1)) &
                            + a7(2)*(var3(i,j,k+2)-var3(i,j,k-2)) &
                            + a7(3)*(var3(i,j,k+3)-var3(i,j,k-3)) &
                            + dvar3(i,j,k)
             enddo
          enddo
       endif
       if (ngh==5) then
          ! sixth-order (7-pt stencil) for point k=nz-4
          k=nz-4
          do j=ndyt,nfyt
             do i=ndxt,nfxt
                dvar1(i,j,k)= a7(1)*(var1(i,j,k+1)-var1(i,j,k-1)) &
                            + a7(2)*(var1(i,j,k+2)-var1(i,j,k-2)) &
                            + a7(3)*(var1(i,j,k+3)-var1(i,j,k-3)) &
                            + dvar1(i,j,k)
                dvar2(i,j,k)= a7(1)*(var2(i,j,k+1)-var2(i,j,k-1)) &
                            + a7(2)*(var2(i,j,k+2)-var2(i,j,k-2)) &
                            + a7(3)*(var2(i,j,k+3)-var2(i,j,k-3)) &
                            + dvar2(i,j,k)
                dvar3(i,j,k)= a7(1)*(var3(i,j,k+1)-var3(i,j,k-1)) &
                            + a7(2)*(var3(i,j,k+2)-var3(i,j,k-2)) &
                            + a7(3)*(var3(i,j,k+3)-var3(i,j,k-3)) &
                            + dvar3(i,j,k)
             enddo
          enddo
       endif
    endif

    ! Metrics at interior points
    ! ==========================
    select case (ngh)
    case (5)
       a_=a11
    case (4)
       a_=a9
    case (3)
       a_=a7
    case (2)
       a_=a5
    case (1)
       a_=a3        
    case default
       call mpistop('max ngh is 5 (11pts-stencil)!', 0)
    end select

    do k=ndz,nfz
       do j=ndyt,nfyt
          do i=ndxt,nfxt
             do l=-ngh,ngh
                dvar1(i,j,k)= dvar1(i,j,k) + a_(l)*var1(i,j,k+l)
                dvar2(i,j,k)= dvar2(i,j,k) + a_(l)*var2(i,j,k+l)
                dvar3(i,j,k)= dvar3(i,j,k) + a_(l)*var3(i,j,k+l)
             enddo
          enddo
       enddo
    enddo

  end subroutine derivative_phi

  !===============================================================================
  module subroutine derivative_phi_v(var1,dvar1,var2,dvar2,var3,dvar3)
  !===============================================================================
    !> Compute derivatives along phi for metrics (full 3D version)
    !> - stencil -ngh_v:+ngh_v for viscous fluxes -
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v), intent(in) :: var1,var2,var3
    real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v), intent(inout) :: dvar1,dvar2,dvar3
    ! ---------------------------------------------------------------------------
    integer :: i,j,k,l
    integer :: ndzv,nfzv
    real(wp), dimension(-ngh_v:ngh_v) :: a_ ! coefficients of interior FD scheme
    ! ---------------------------------------------------------------------------

    ! Metrics for BC at kmin
    ! ======================
    ndzv=1
    
    if (BC_face(3,1)%sort<=0) then
       select case (ngh_v)
       case (2)
          ndzv=3
          k=1
          do j=ndy_v1,nfy_v1
             do i=ndx_v1,nfx_v1
                do l=0,4
                   dvar1(i,j,k)=dvar1(i,j,k)+a04(1+l)*var1(i,j,k+l)
                   dvar2(i,j,k)=dvar2(i,j,k)+a04(1+l)*var2(i,j,k+l)
                   dvar3(i,j,k)=dvar3(i,j,k)+a04(1+l)*var3(i,j,k+l)
                enddo
             enddo
          enddo
          k=2
          do j=ndy_v1,nfy_v1
             do i=ndx_v1,nfx_v1
                do l=-1,3
                   dvar1(i,j,k)=dvar1(i,j,k)+a13(2+l)*var1(i,j,k+l)
                   dvar2(i,j,k)=dvar2(i,j,k)+a13(2+l)*var2(i,j,k+l)
                   dvar3(i,j,k)=dvar3(i,j,k)+a13(2+l)*var3(i,j,k+l)
                enddo
             enddo
          enddo
       case (1)
          ndzv=2
          k=1
          do j=ndy_v1,nfy_v1
             do i=ndx_v1,nfx_v1
                do l=0,1
                   dvar1(i,j,k)=dvar1(i,j,k)+a01(1+l)*var1(i,j,k+l)
                   dvar2(i,j,k)=dvar2(i,j,k)+a01(1+l)*var2(i,j,k+l)
                   dvar3(i,j,k)=dvar3(i,j,k)+a01(1+l)*var3(i,j,k+l)
                enddo
             enddo
          enddo
       case default
          call mpistop('max ngh_v is 2 (5pts-stencil)!', 0)
       end select
    endif

    ! Metrics for BC at kmax
    ! ======================
    nfzv=nz
    
    if (BC_face(3,2)%sort<=0) then
       select case (ngh_v)
       case (2)
          nfzv=nz-2
          k=nz-1
          do j=ndy_v1,nfy_v1
             do i=ndx_v1,nfx_v1
                do l=-3,1
                   dvar1(i,j,k)=dvar1(i,j,k)-a13(2-l)*var1(i,j,k+l)
                   dvar2(i,j,k)=dvar2(i,j,k)-a13(2-l)*var2(i,j,k+l)
                   dvar3(i,j,k)=dvar3(i,j,k)-a13(2-l)*var3(i,j,k+l)
                enddo
             enddo
          enddo
          k=nz
          do j=ndy_v1,nfy_v1
             do i=ndx_v1,nfx_v1
                do l=-4,0
                   dvar1(i,j,k)=dvar1(i,j,k)-a04(1-l)*var1(i,j,k+l)
                   dvar2(i,j,k)=dvar2(i,j,k)-a04(1-l)*var2(i,j,k+l)
                   dvar3(i,j,k)=dvar3(i,j,k)-a04(1-l)*var3(i,j,k+l)
                enddo
             enddo
          enddo
       case (1)
          nfzv=nz-1
          k=nz
          do j=ndy_v1,nfy_v1
             do i=ndx_v1,nfx_v1
                do l=-1,0
                   dvar1(i,j,k)=dvar1(i,j,k)-a01(1-l)*var1(i,j,k+l)
                   dvar2(i,j,k)=dvar2(i,j,k)-a01(1-l)*var2(i,j,k+l)
                   dvar3(i,j,k)=dvar3(i,j,k)-a01(1-l)*var3(i,j,k+l)
                enddo
             enddo
          enddo
       case default
          call mpistop('max ngh_v is 2 (5pts-stencil)!', 0)
       end select
    endif

    ! Metrics at interior points
    ! ==========================
    select case (ngh_v)
    case (2)
       a_=a5
    case (1)
       a_=a3        
    case default
       call mpistop('max ngh_v is 2 (5pts-stencil)!', 0)
    end select
    
    do k=ndzv,nfzv
       do j=ndy_v1,nfy_v1
          do i=ndx_v1,nfx_v1
             do l=-ngh_v,ngh_v
                dvar1(i,j,k)= dvar1(i,j,k) + a_(l)*var1(i,j,k+l)
                dvar2(i,j,k)= dvar2(i,j,k) + a_(l)*var2(i,j,k+l)
                dvar3(i,j,k)= dvar3(i,j,k) + a_(l)*var3(i,j,k+l)
             enddo
          enddo
       enddo
    enddo   

  end subroutine derivative_phi_v
  
end submodule smod_grid_metrics_deriv_phi_c3
