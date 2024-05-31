!=================================================================================
submodule (mod_grid_metrics_c3) smod_grid_metrics_deriv_ksi_c3
!=================================================================================
  !> Module to compute full 3D curvilinear metrics
  !> with Geometric Conservation Law (GCL)
!=================================================================================

contains

  !===============================================================================
  module subroutine derivative_ksi(var1,dvar1,var2,dvar2,var3,dvar3)
  !===============================================================================
    !> Compute derivatives along ksi for metrics (full 3D version)
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
   
    ! Metrics for BC at imin
    ! ======================
    
    ! For wall & charac, inlet and outlet BCs
    ! ---------------------------------------
    if ((BC_face(1,1)%sort==0).or.(BC_face(1,1)%sort<=-3)) then
       
       ! SBP boundary schemes
       ! --------------------
       if (is_SBP) then
          if (ngh<4) then
             call mpistop('3D SBP metrics need yet at least 9pt-stencil !', 0)
          else
             ! Point #1: stencil [o x x x]
             !i=1
             do k=ndzt,nfzt
                do j=ndyt,nfyt
                   dvar1(1,j,k)= as4p0(1)*var1(1,j,k)+as4p0(2)*var1(2,j,k) &
                               + as4p0(3)*var1(3,j,k)+as4p0(4)*var1(4,j,k) &
                               + dvar1(1,j,k)
                   dvar2(1,j,k)= as4p0(1)*var2(1,j,k)+as4p0(2)*var2(2,j,k) &
                               + as4p0(3)*var2(3,j,k)+as4p0(4)*var2(4,j,k) &
                               + dvar2(1,j,k)
                   dvar3(1,j,k)= as4p0(1)*var3(1,j,k)+as4p0(2)*var3(2,j,k) &
                               + as4p0(3)*var3(3,j,k)+as4p0(4)*var3(4,j,k) &
                               + dvar3(1,j,k)
                enddo
             enddo
             ! point #2: stencil [x o x x]
             !i=2
             do k=ndzt,nfzt
                do j=ndyt,nfyt
                   dvar1(2,j,k)= as4p1(1)*var1(1,j,k)+as4p1(2)*var1(2,j,k) &
                               + as4p1(3)*var1(3,j,k)+as4p1(4)*var1(4,j,k) &
                               + dvar1(2,j,k)
                   dvar2(2,j,k)= as4p1(1)*var2(1,j,k)+as4p1(2)*var2(2,j,k) &
                               + as4p1(3)*var2(3,j,k)+as4p1(4)*var2(4,j,k) &
                               + dvar2(2,j,k)
                   dvar3(2,j,k)= as4p1(1)*var3(1,j,k)+as4p1(2)*var3(2,j,k) &
                               + as4p1(3)*var3(3,j,k)+as4p1(4)*var3(4,j,k) &
                               + dvar3(2,j,k)
                enddo
             enddo
             ! point #3: stencil [x x o x x]
             !i=3
             do k=ndzt,nfzt
                do j=ndyt,nfyt
                   dvar1(3,j,k)= as4p2(1)*var1(1,j,k)+as4p2(2)*var1(2,j,k) &
                               + as4p2(3)*var1(3,j,k)+as4p2(4)*var1(4,j,k) &
                               + as4p2(5)*var1(5,j,k) + dvar1(3,j,k)
                   dvar2(3,j,k)= as4p2(1)*var2(1,j,k)+as4p2(2)*var2(2,j,k) &
                               + as4p2(3)*var2(3,j,k)+as4p2(4)*var2(4,j,k) &
                               + as4p2(5)*var2(5,j,k) + dvar2(3,j,k)
                   dvar3(3,j,k)= as4p2(1)*var3(1,j,k)+as4p2(2)*var3(2,j,k) &
                               + as4p2(3)*var3(3,j,k)+as4p2(4)*var3(4,j,k) &
                               + as4p2(5)*var3(5,j,k) + dvar3(3,j,k)
                enddo
             enddo
             ! point #4: stencil [x x x o x x]
             !i=4
             do k=ndzt,nfzt
                do j=ndyt,nfyt
                   dvar1(4,j,k)= as4p3(1)*var1(1,j,k)+as4p3(2)*var1(2,j,k) &
                               + as4p3(3)*var1(3,j,k)+as4p3(4)*var1(4,j,k) &
                               + as4p3(5)*var1(5,j,k)+as4p3(6)*var1(6,j,k) &
                               + dvar1(4,j,k)
                   dvar2(4,j,k)= as4p3(1)*var2(1,j,k)+as4p3(2)*var2(2,j,k) &
                               + as4p3(3)*var2(3,j,k)+as4p3(4)*var2(4,j,k) &
                               + as4p3(5)*var2(5,j,k)+as4p3(6)*var2(6,j,k) &
                               + dvar2(4,j,k)
                   dvar3(4,j,k)= as4p3(1)*var3(1,j,k)+as4p3(2)*var3(2,j,k) &
                               + as4p3(3)*var3(3,j,k)+as4p3(4)*var3(4,j,k) &
                               + as4p3(5)*var3(5,j,k)+as4p3(6)*var3(6,j,k) &
                               + dvar3(4,j,k)
                enddo
             enddo             
          endif
          if (ngh>=5) then
             ! eighth-order (9-pt stencil) for point i=5
             i=5
             do k=ndzt,nfzt
                do j=ndyt,nfyt
                   dvar1(i,j,k)= a9(1)*(var1(i+1,j,k)-var1(i-1,j,k)) &
                               + a9(2)*(var1(i+2,j,k)-var1(i-2,j,k)) &
                               + a9(3)*(var1(i+3,j,k)-var1(i-3,j,k)) &
                               + a9(4)*(var1(i+4,j,k)-var1(i-4,j,k)) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= a9(1)*(var2(i+1,j,k)-var2(i-1,j,k)) &
                               + a9(2)*(var2(i+2,j,k)-var2(i-2,j,k)) &
                               + a9(3)*(var2(i+3,j,k)-var2(i-3,j,k)) &
                               + a9(4)*(var2(i+4,j,k)-var2(i-4,j,k)) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= a9(1)*(var3(i+1,j,k)-var3(i-1,j,k)) &
                               + a9(2)*(var3(i+2,j,k)-var3(i-2,j,k)) &
                               + a9(3)*(var3(i+3,j,k)-var3(i-3,j,k)) &
                               + a9(4)*(var3(i+4,j,k)-var3(i-4,j,k)) &
                               + dvar3(i,j,k)
                enddo
             enddo
          endif
          
       ! Reduced order stencil
       ! ---------------------
       else       
          if (ngh>=1) then
             ! first-order (2-pt stencil) for point i=1
             i=1
             do k=ndzt,nfzt
                do j=ndyt,nfyt
                   !do l=0,1
                   !   dvar1(i,j,k)=dvar1(i,j,k)+a01(1+l)*var1(i+l,j,k)
                   !   dvar2(i,j,k)=dvar2(i,j,k)+a01(1+l)*var2(i+l,j,k)
                   !   dvar3(i,j,k)=dvar3(i,j,k)+a01(1+l)*var3(i+l,j,k)
                   !enddo
                   do l=0,2
                      dvar1(i,j,k)=dvar1(i,j,k)+a02(1+l)*var1(i+l,j,k)
                      dvar2(i,j,k)=dvar2(i,j,k)+a02(1+l)*var2(i+l,j,k)
                      dvar3(i,j,k)=dvar3(i,j,k)+a02(1+l)*var3(i+l,j,k)
                   enddo
                enddo
             enddo
          endif
          if (ngh>=2) then
             ! second-order (3-pt stencil) for point i=2
             i=2
             do k=ndzt,nfzt
                do j=ndyt,nfyt
                   dvar1(i,j,k)= a3(1)*(var1(i+1,j,k)-var1(i-1,j,k)) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= a3(1)*(var2(i+1,j,k)-var2(i-1,j,k)) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= a3(1)*(var3(i+1,j,k)-var3(i-1,j,k)) &
                               + dvar3(i,j,k)
                enddo
             enddo
          endif
          if (ngh>=3) then
             ! fourth-order (5-pt stencil) for point i=3
             i=3
             do k=ndzt,nfzt
                do j=ndyt,nfyt
                   dvar1(i,j,k)= a5(1)*(var1(i+1,j,k)-var1(i-1,j,k)) &
                               + a5(2)*(var1(i+2,j,k)-var1(i-2,j,k)) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= a5(1)*(var2(i+1,j,k)-var2(i-1,j,k)) &
                               + a5(2)*(var2(i+2,j,k)-var2(i-2,j,k)) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= a5(1)*(var3(i+1,j,k)-var3(i-1,j,k)) &
                               + a5(2)*(var3(i+2,j,k)-var3(i-2,j,k)) &
                               + dvar3(i,j,k)
                enddo
             enddo
          endif
          if (ngh>=4) then
             ! sixth-order (7-pt stencil) for point i=4
             i=4
             do k=ndzt,nfzt
                do j=ndyt,nfyt
                   dvar1(i,j,k)= a7(1)*(var1(i+1,j,k)-var1(i-1,j,k)) &
                               + a7(2)*(var1(i+2,j,k)-var1(i-2,j,k)) &
                               + a7(3)*(var1(i+3,j,k)-var1(i-3,j,k)) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= a7(1)*(var2(i+1,j,k)-var2(i-1,j,k)) &
                               + a7(2)*(var2(i+2,j,k)-var2(i-2,j,k)) &
                               + a7(3)*(var2(i+3,j,k)-var2(i-3,j,k)) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= a7(1)*(var3(i+1,j,k)-var3(i-1,j,k)) &
                               + a7(2)*(var3(i+2,j,k)-var3(i-2,j,k)) &
                               + a7(3)*(var3(i+3,j,k)-var3(i-3,j,k)) &
                               + dvar3(i,j,k)
                enddo
             enddo
          endif
          if (ngh==5) then
             ! eighth-order (9-pt stencil) for point i=5
             i=5
             do k=ndzt,nfzt
                do j=ndyt,nfyt
                   dvar1(i,j,k)= a9(1)*(var1(i+1,j,k)-var1(i-1,j,k)) &
                               + a9(2)*(var1(i+2,j,k)-var1(i-2,j,k)) &
                               + a9(3)*(var1(i+3,j,k)-var1(i-3,j,k)) &
                               + a9(4)*(var1(i+4,j,k)-var1(i-4,j,k)) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= a9(1)*(var2(i+1,j,k)-var2(i-1,j,k)) &
                               + a9(2)*(var2(i+2,j,k)-var2(i-2,j,k)) &
                               + a9(3)*(var2(i+3,j,k)-var2(i-3,j,k)) &
                               + a9(4)*(var2(i+4,j,k)-var2(i-4,j,k)) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= a9(1)*(var3(i+1,j,k)-var3(i-1,j,k)) &
                               + a9(2)*(var3(i+2,j,k)-var3(i-2,j,k)) &
                               + a9(3)*(var3(i+3,j,k)-var3(i-3,j,k)) &
                               + a9(4)*(var3(i+4,j,k)-var3(i-4,j,k)) &
                               + dvar3(i,j,k)
                enddo
             enddo
          endif
       endif
    endif
    
    ! 7-point stencil for Tam&Dong BC
    ! -------------------------------
    if ((BC_face(1,1)%sort==-1).or.(BC_face(1,1)%sort==-2)) then
       if (ngh<3) then
          call mpistop('Tam&Dong BCs need at least 7pt-stencil !', 0)
       else
          i=1
          do k=ndzt,nfzt
             do j=ndyt,nfyt
                do l=0,6
                   dvar1(i,j,k)=dvar1(i,j,k)+a06(1+l)*var1(i+l,j,k)
                   dvar2(i,j,k)=dvar2(i,j,k)+a06(1+l)*var2(i+l,j,k)
                   dvar3(i,j,k)=dvar3(i,j,k)+a06(1+l)*var3(i+l,j,k)
                enddo
             enddo
          enddo
          i=2
          do k=ndzt,nfzt
             do j=ndyt,nfyt
                do l=-1,5
                   dvar1(i,j,k)=dvar1(i,j,k)+a15(2+l)*var1(i+l,j,k)
                   dvar2(i,j,k)=dvar2(i,j,k)+a15(2+l)*var2(i+l,j,k)
                   dvar3(i,j,k)=dvar3(i,j,k)+a15(2+l)*var3(i+l,j,k)
                enddo
             enddo
          enddo
          i=3
          do k=ndzt,nfzt
             do j=ndyt,nfyt
                do l=-2,4
                   dvar1(i,j,k)=dvar1(i,j,k)+a24(3+l)*var1(i+l,j,k)
                   dvar2(i,j,k)=dvar2(i,j,k)+a24(3+l)*var2(i+l,j,k)
                   dvar3(i,j,k)=dvar3(i,j,k)+a24(3+l)*var3(i+l,j,k)
                enddo
             enddo
          enddo
       endif
       if (ngh>=4) then
          ! sixth-order (7-pt stencil) for point i=4
          i=4
          do k=ndzt,nfzt
             do j=ndyt,nfyt
                dvar1(i,j,k)= a7(1)*(var1(i+1,j,k)-var1(i-1,j,k)) &
                            + a7(2)*(var1(i+2,j,k)-var1(i-2,j,k)) &
                            + a7(3)*(var1(i+3,j,k)-var1(i-3,j,k)) &
                            + dvar1(i,j,k)
                dvar2(i,j,k)= a7(1)*(var2(i+1,j,k)-var2(i-1,j,k)) &
                            + a7(2)*(var2(i+2,j,k)-var2(i-2,j,k)) &
                            + a7(3)*(var2(i+3,j,k)-var2(i-3,j,k)) &
                            + dvar2(i,j,k)
                dvar3(i,j,k)= a7(1)*(var3(i+1,j,k)-var3(i-1,j,k)) &
                            + a7(2)*(var3(i+2,j,k)-var3(i-2,j,k)) &
                            + a7(3)*(var3(i+3,j,k)-var3(i-3,j,k)) &
                            + dvar3(i,j,k)
             enddo
          enddo
       endif
       if (ngh==5) then
          ! sixth-order (7-pt stencil) for point i=5
          i=5
          do k=ndzt,nfzt
             do j=ndyt,nfyt
                dvar1(i,j,k)= a7(1)*(var1(i+1,j,k)-var1(i-1,j,k)) &
                            + a7(2)*(var1(i+2,j,k)-var1(i-2,j,k)) &
                            + a7(3)*(var1(i+3,j,k)-var1(i-3,j,k)) &
                            + dvar1(i,j,k)
                dvar2(i,j,k)= a7(1)*(var2(i+1,j,k)-var2(i-1,j,k)) &
                            + a7(2)*(var2(i+2,j,k)-var2(i-2,j,k)) &
                            + a7(3)*(var2(i+3,j,k)-var2(i-3,j,k)) &
                            + dvar2(i,j,k)
                dvar3(i,j,k)= a7(1)*(var3(i+1,j,k)-var3(i-1,j,k)) &
                            + a7(2)*(var3(i+2,j,k)-var3(i-2,j,k)) &
                            + a7(3)*(var3(i+3,j,k)-var3(i-3,j,k)) &
                            + dvar3(i,j,k)
             enddo
          enddo
       endif
    endif

    ! Metrics for BC at imax
    ! ======================
    
    ! For wall & charac, inlet and outlet BCs
    ! ---------------------------------------
    if ((BC_face(1,2)%sort==0).or.(BC_face(1,2)%sort<=-3)) then
       
       ! SBP boundary schemes
       ! --------------------
       if (is_SBP) then
          if (ngh<4) then
             call mpistop('3D SBP metrics need yet at least 9pt-stencil !', 0)
          else
             ! point #1: stencil [o x x x]
             i=nx
             do k=ndzt,nfzt
                do j=ndyt,nfyt
                   dvar1(i,j,k)= as4m0(1)*var1(nx  ,j,k)+as4m0(2)*var1(nx-1,j,k) &
                               + as4m0(3)*var1(nx-2,j,k)+as4m0(4)*var1(nx-3,j,k) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= as4m0(1)*var2(nx  ,j,k)+as4m0(2)*var2(nx-1,j,k) &
                               + as4m0(3)*var2(nx-2,j,k)+as4m0(4)*var2(nx-3,j,k) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= as4m0(1)*var3(nx  ,j,k)+as4m0(2)*var3(nx-1,j,k) &
                               + as4m0(3)*var3(nx-2,j,k)+as4m0(4)*var3(nx-3,j,k) &
                               + dvar3(i,j,k)
                enddo
             enddo
             ! point #2: stencil [x o x x]
             i=nx-1
             do k=ndzt,nfzt
                do j=ndyt,nfyt
                   dvar1(i,j,k)= as4m1(1)*var1(nx  ,j,k)+as4m1(2)*var1(nx-1,j,k) &
                               + as4m1(3)*var1(nx-2,j,k)+as4m1(4)*var1(nx-3,j,k) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= as4m1(1)*var2(nx  ,j,k)+as4m1(2)*var2(nx-1,j,k) &
                               + as4m1(3)*var2(nx-2,j,k)+as4m1(4)*var2(nx-3,j,k) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= as4m1(1)*var3(nx  ,j,k)+as4m1(2)*var3(nx-1,j,k) &
                               + as4m1(3)*var3(nx-2,j,k)+as4m1(4)*var3(nx-3,j,k) &
                               + dvar3(i,j,k)
                enddo
             enddo
             ! point #3: stencil [x x o x x]
             i=nx-2
             do k=ndzt,nfzt
                do j=ndyt,nfyt
                   dvar1(i,j,k)= as4m2(1)*var1(nx  ,j,k)+as4m2(2)*var1(nx-1,j,k) &
                               + as4m2(3)*var1(nx-2,j,k)+as4m2(4)*var1(nx-3,j,k) &
                               + as4m2(5)*var1(nx-4,j,k) + dvar1(i,j,k)
                   dvar2(i,j,k)= as4m2(1)*var2(nx  ,j,k)+as4m2(2)*var2(nx-1,j,k) &
                               + as4m2(3)*var2(nx-2,j,k)+as4m2(4)*var2(nx-3,j,k) &
                               + as4m2(5)*var2(nx-4,j,k) + dvar2(i,j,k)
                   dvar3(i,j,k)= as4m2(1)*var3(nx  ,j,k)+as4m2(2)*var3(nx-1,j,k) &
                               + as4m2(3)*var3(nx-2,j,k)+as4m2(4)*var3(nx-3,j,k) &
                               + as4m2(5)*var3(nx-4,j,k) + dvar3(i,j,k)
                enddo
             enddo
             ! point #4: stencil [x x x o x x]
             i=nx-3
             do k=ndzt,nfzt
                do j=ndyt,nfyt
                   dvar1(i,j,k)= as4m3(1)*var1(nx  ,j,k)+as4m3(2)*var1(nx-1,j,k) &
                               + as4m3(3)*var1(nx-2,j,k)+as4m3(4)*var1(nx-3,j,k) &
                               + as4m3(5)*var1(nx-4,j,k)+as4m3(6)*var1(nx-5,j,k) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= as4m3(1)*var2(nx  ,j,k)+as4m3(2)*var2(nx-1,j,k) &
                               + as4m3(3)*var2(nx-2,j,k)+as4m3(4)*var2(nx-3,j,k) &
                               + as4m3(5)*var2(nx-4,j,k)+as4m3(6)*var2(nx-5,j,k) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= as4m3(1)*var3(nx  ,j,k)+as4m3(2)*var3(nx-1,j,k) &
                               + as4m3(3)*var3(nx-2,j,k)+as4m3(4)*var3(nx-3,j,k) &
                               + as4m3(5)*var3(nx-4,j,k)+as4m3(6)*var3(nx-5,j,k) &
                               + dvar3(i,j,k)
                enddo
             enddo
          endif
          if (ngh>=5) then
             ! eighth-order (9-pt stencil) for point i=nx-4
             i=nx-4
             do k=ndzt,nfzt
                do j=ndyt,nfyt
                   dvar1(i,j,k)= a9(1)*(var1(i+1,j,k)-var1(i-1,j,k)) &
                               + a9(2)*(var1(i+2,j,k)-var1(i-2,j,k)) &
                               + a9(3)*(var1(i+3,j,k)-var1(i-3,j,k)) &
                               + a9(4)*(var1(i+4,j,k)-var1(i-4,j,k)) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= a9(1)*(var2(i+1,j,k)-var2(i-1,j,k)) &
                               + a9(2)*(var2(i+2,j,k)-var2(i-2,j,k)) &
                               + a9(3)*(var2(i+3,j,k)-var2(i-3,j,k)) &
                               + a9(4)*(var2(i+4,j,k)-var2(i-4,j,k)) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= a9(1)*(var3(i+1,j,k)-var3(i-1,j,k)) &
                               + a9(2)*(var3(i+2,j,k)-var3(i-2,j,k)) &
                               + a9(3)*(var3(i+3,j,k)-var3(i-3,j,k)) &
                               + a9(4)*(var3(i+4,j,k)-var3(i-4,j,k)) &
                               + dvar3(i,j,k)
                enddo
             enddo
          endif
          
       ! Reduced order stencil
       ! ---------------------
       else       
          if (ngh>=1) then
             ! first-order (2-pt stencil) for point i=nx
             i=nx
             do k=ndzt,nfzt
                do j=ndyt,nfyt
                   !do l=-1,0
                   !   dvar1(i,j,k)=dvar1(i,j,k)-a01(1-l)*var1(i+l,j,k)
                   !   dvar2(i,j,k)=dvar2(i,j,k)-a01(1-l)*var2(i+l,j,k)
                   !   dvar3(i,j,k)=dvar3(i,j,k)-a01(1-l)*var3(i+l,j,k)
                   !enddo
                   do l=-2,0
                      dvar1(i,j,k)=dvar1(i,j,k)-a02(1-l)*var1(i+l,j,k)
                      dvar2(i,j,k)=dvar2(i,j,k)-a02(1-l)*var2(i+l,j,k)
                      dvar3(i,j,k)=dvar3(i,j,k)-a02(1-l)*var3(i+l,j,k)
                   enddo
                enddo
             enddo
          endif
          if (ngh>=2) then
             ! second-order (3-pt stencil) for point i=nx-1
             i=nx-1
             do k=ndzt,nfzt
                do j=ndyt,nfyt
                   dvar1(i,j,k)= a3(1)*(var1(i+1,j,k)-var1(i-1,j,k)) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= a3(1)*(var2(i+1,j,k)-var2(i-1,j,k)) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= a3(1)*(var3(i+1,j,k)-var3(i-1,j,k)) &
                               + dvar3(i,j,k)
                enddo
             enddo
          endif
          if (ngh>=3) then
             ! fourth-order (5-pt stencil) for point i=nx-2
             i=nx-2
             do k=ndzt,nfzt
                do j=ndyt,nfyt
                   dvar1(i,j,k)= a5(1)*(var1(i+1,j,k)-var1(i-1,j,k)) &
                               + a5(2)*(var1(i+2,j,k)-var1(i-2,j,k)) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= a5(1)*(var2(i+1,j,k)-var2(i-1,j,k)) &
                               + a5(2)*(var2(i+2,j,k)-var2(i-2,j,k)) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= a5(1)*(var3(i+1,j,k)-var3(i-1,j,k)) &
                               + a5(2)*(var3(i+2,j,k)-var3(i-2,j,k)) &
                               + dvar3(i,j,k)
                enddo
             enddo
          endif
          if (ngh>=4) then
             ! sixth-order (7-pt stencil) for point i=nx-3
             i=nx-3
             do k=ndzt,nfzt
                do j=ndyt,nfyt
                   dvar1(i,j,k)= a7(1)*(var1(i+1,j,k)-var1(i-1,j,k)) &
                               + a7(2)*(var1(i+2,j,k)-var1(i-2,j,k)) &
                               + a7(3)*(var1(i+3,j,k)-var1(i-3,j,k)) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= a7(1)*(var2(i+1,j,k)-var2(i-1,j,k)) &
                               + a7(2)*(var2(i+2,j,k)-var2(i-2,j,k)) &
                               + a7(3)*(var2(i+3,j,k)-var2(i-3,j,k)) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= a7(1)*(var3(i+1,j,k)-var3(i-1,j,k)) &
                               + a7(2)*(var3(i+2,j,k)-var3(i-2,j,k)) &
                               + a7(3)*(var3(i+3,j,k)-var3(i-3,j,k)) &
                               + dvar3(i,j,k)
                enddo
             enddo
          endif
          if (ngh==5) then
             ! eighth-order (9-pt stencil) for point i=nx-4
             i=nx-4
             do k=ndzt,nfzt
                do j=ndyt,nfyt
                   dvar1(i,j,k)= a9(1)*(var1(i+1,j,k)-var1(i-1,j,k)) &
                               + a9(2)*(var1(i+2,j,k)-var1(i-2,j,k)) &
                               + a9(3)*(var1(i+3,j,k)-var1(i-3,j,k)) &
                               + a9(4)*(var1(i+4,j,k)-var1(i-4,j,k)) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= a9(1)*(var2(i+1,j,k)-var2(i-1,j,k)) &
                               + a9(2)*(var2(i+2,j,k)-var2(i-2,j,k)) &
                               + a9(3)*(var2(i+3,j,k)-var2(i-3,j,k)) &
                               + a9(4)*(var2(i+4,j,k)-var2(i-4,j,k)) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= a9(1)*(var3(i+1,j,k)-var3(i-1,j,k)) &
                               + a9(2)*(var3(i+2,j,k)-var3(i-2,j,k)) &
                               + a9(3)*(var3(i+3,j,k)-var3(i-3,j,k)) &
                               + a9(4)*(var3(i+4,j,k)-var3(i-4,j,k)) &
                               + dvar3(i,j,k)
                enddo
             enddo
          endif
       endif
    endif
    
    ! 7-point stencil for Tam&Dong BC
    ! -------------------------------
    if ((BC_face(1,2)%sort==-1).or.(BC_face(1,2)%sort==-2)) then
       if (ngh<3) then
          call mpistop('Tam&Dong BCs need at least 7pt-stencil !', 0)
       else
          i=nx
          do k=ndzt,nfzt
             do j=ndyt,nfyt
                do l=-6,0
                   dvar1(i,j,k)=dvar1(i,j,k)-a06(1-l)*var1(i+l,j,k)
                   dvar2(i,j,k)=dvar2(i,j,k)-a06(1-l)*var2(i+l,j,k)
                   dvar3(i,j,k)=dvar3(i,j,k)-a06(1-l)*var3(i+l,j,k)
                enddo
             enddo
          enddo
          i=nx-1
          do k=ndzt,nfzt
             do j=ndyt,nfyt
                do l=-5,1
                   dvar1(i,j,k)=dvar1(i,j,k)-a15(2-l)*var1(i+l,j,k)
                   dvar2(i,j,k)=dvar2(i,j,k)-a15(2-l)*var2(i+l,j,k)
                   dvar3(i,j,k)=dvar3(i,j,k)-a15(2-l)*var3(i+l,j,k)
                enddo
             enddo
          enddo
          i=nx-2
          do k=ndzt,nfzt
             do j=ndyt,nfyt
                do l=-4,2
                   dvar1(i,j,k)=dvar1(i,j,k)-a24(3-l)*var1(i+l,j,k)
                   dvar2(i,j,k)=dvar2(i,j,k)-a24(3-l)*var2(i+l,j,k)
                   dvar3(i,j,k)=dvar3(i,j,k)-a24(3-l)*var3(i+l,j,k)
                enddo
             enddo
          enddo
       endif
       if (ngh>=4) then
          ! sixth-order (7-pt stencil) for point i=nx-3
          i=nx-3
          do k=ndzt,nfzt
             do j=ndyt,nfyt
                dvar1(i,j,k)= a7(1)*(var1(i+1,j,k)-var1(i-1,j,k)) &
                            + a7(2)*(var1(i+2,j,k)-var1(i-2,j,k)) &
                            + a7(3)*(var1(i+3,j,k)-var1(i-3,j,k)) &
                            + dvar1(i,j,k)
                dvar2(i,j,k)= a7(1)*(var2(i+1,j,k)-var2(i-1,j,k)) &
                            + a7(2)*(var2(i+2,j,k)-var2(i-2,j,k)) &
                            + a7(3)*(var2(i+3,j,k)-var2(i-3,j,k)) &
                            + dvar2(i,j,k)
                dvar3(i,j,k)= a7(1)*(var3(i+1,j,k)-var3(i-1,j,k)) &
                            + a7(2)*(var3(i+2,j,k)-var3(i-2,j,k)) &
                            + a7(3)*(var3(i+3,j,k)-var3(i-3,j,k)) &
                            + dvar3(i,j,k)
             enddo
          enddo
       endif
       if (ngh==5) then
          ! sixth-order (7-pt stencil) for point i=nx-4
          i=nx-4
          do k=ndzt,nfzt
             do j=ndyt,nfyt
                dvar1(i,j,k)= a7(1)*(var1(i+1,j,k)-var1(i-1,j,k)) &
                            + a7(2)*(var1(i+2,j,k)-var1(i-2,j,k)) &
                            + a7(3)*(var1(i+3,j,k)-var1(i-3,j,k)) &
                            + dvar1(i,j,k)
                dvar2(i,j,k)= a7(1)*(var2(i+1,j,k)-var2(i-1,j,k)) &
                            + a7(2)*(var2(i+2,j,k)-var2(i-2,j,k)) &
                            + a7(3)*(var2(i+3,j,k)-var2(i-3,j,k)) &
                            + dvar2(i,j,k)
                dvar3(i,j,k)= a7(1)*(var3(i+1,j,k)-var3(i-1,j,k)) &
                            + a7(2)*(var3(i+2,j,k)-var3(i-2,j,k)) &
                            + a7(3)*(var3(i+3,j,k)-var3(i-3,j,k)) &
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
       call mpistop('max ngh is 5 (11pt-stencil)!', 0)
    end select

    do k=ndzt,nfzt
       do j=ndyt,nfyt
          do i=ndx,nfx
             do l=-ngh,ngh
                dvar1(i,j,k)= dvar1(i,j,k) + a_(l)*var1(i+l,j,k)
                dvar2(i,j,k)= dvar2(i,j,k) + a_(l)*var2(i+l,j,k)
                dvar3(i,j,k)= dvar3(i,j,k) + a_(l)*var3(i+l,j,k)
             enddo
          enddo
       enddo
    enddo

  end subroutine derivative_ksi
  
  !===============================================================================
  module subroutine derivative_ksi_v(var1,dvar1,var2,dvar2,var3,dvar3)
  !===============================================================================
    !> Compute derivatives along ksi for metrics (full 3D version)
    !> - stencil -ngh_v:+ngh_v for viscous fluxes -
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v), intent(in) :: var1,var2,var3
    real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v), intent(inout) :: dvar1,dvar2,dvar3
    ! ---------------------------------------------------------------------------
    integer :: i,j,k,l
    integer :: ndxv,nfxv
    real(wp), dimension(-ngh_v:ngh_v) :: a_ ! coefficients of interior FD scheme
    ! ---------------------------------------------------------------------------
   
    ! Metrics for BC at imin
    ! ======================
    ndxv=1

    if (BC_face(1,1)%sort<=0) then     
       select case (ngh_v)
       case (2)
          ndxv=3
          i=1
          do k=ndz_v1,nfz_v1
             do j=ndy_v1,nfy_v1
                do l=0,4
                   dvar1(i,j,k)=dvar1(i,j,k)+a04(1+l)*var1(i+l,j,k)
                   dvar2(i,j,k)=dvar2(i,j,k)+a04(1+l)*var2(i+l,j,k)
                   dvar3(i,j,k)=dvar3(i,j,k)+a04(1+l)*var3(i+l,j,k)
                enddo
             enddo
          enddo
          i=2
          do k=ndz_v1,nfz_v1
             do j=ndy_v1,nfy_v1
                do l=-1,3
                   dvar1(i,j,k)=dvar1(i,j,k)+a13(2+l)*var1(i+l,j,k)
                   dvar2(i,j,k)=dvar2(i,j,k)+a13(2+l)*var2(i+l,j,k)
                   dvar3(i,j,k)=dvar3(i,j,k)+a13(2+l)*var3(i+l,j,k)
                enddo
             enddo
          enddo
       case (1)
          ndxv=2
          i=1
          do k=ndz_v1,nfz_v1
             do j=ndy_v1,nfy_v1
                do l=0,1
                   dvar1(i,j,k)=dvar1(i,j,k)+a01(1+l)*var1(i+l,j,k)
                   dvar2(i,j,k)=dvar2(i,j,k)+a01(1+l)*var2(i+l,j,k)
                   dvar3(i,j,k)=dvar3(i,j,k)+a01(1+l)*var3(i+l,j,k)
                enddo
             enddo
          enddo
       case default
          call mpistop('max ngh_v is 2 (5pts-stencil)!', 0)
       end select
    endif

    ! Metrics for BC at imax
    ! ======================
    nfxv=nx
    
    if (BC_face(1,2)%sort<=0) then
       select case (ngh_v)
       case (2)
          nfxv=nx-2
          i=nx-1
          do k=ndz_v1,nfz_v1
             do j=ndy_v1,nfy_v1
                do l=-3,1
                   dvar1(i,j,k)=dvar1(i,j,k)-a13(2-l)*var1(i+l,j,k)
                   dvar2(i,j,k)=dvar2(i,j,k)-a13(2-l)*var2(i+l,j,k)
                   dvar3(i,j,k)=dvar3(i,j,k)-a13(2-l)*var3(i+l,j,k)
                enddo
             enddo
          enddo
          i=nx
          do k=ndz_v1,nfz_v1
             do j=ndy_v1,nfy_v1
                do l=-4,0
                   dvar1(i,j,k)=dvar1(i,j,k)-a04(1-l)*var1(i+l,j,k)
                   dvar2(i,j,k)=dvar2(i,j,k)-a04(1-l)*var2(i+l,j,k)
                   dvar3(i,j,k)=dvar3(i,j,k)-a04(1-l)*var3(i+l,j,k)
                enddo
             enddo
          enddo
       case (1)
          nfxv=nx-1
          i=nx
          do k=ndz_v1,nfz_v1
             do j=ndy_v1,nfy_v1
                do l=-1,0
                   dvar1(i,j,k)=dvar1(i,j,k)-a01(1-l)*var1(i+l,j,k)
                   dvar2(i,j,k)=dvar2(i,j,k)-a01(1-l)*var2(i+l,j,k)
                   dvar3(i,j,k)=dvar3(i,j,k)-a01(1-l)*var3(i+l,j,k)
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

    do k=ndz_v1,nfz_v1
       do j=ndy_v1,nfy_v1
          do i=ndxv,nfxv
             do l=-ngh_v,ngh_v
                dvar1(i,j,k)= dvar1(i,j,k) + a_(l)*var1(i+l,j,k)
                dvar2(i,j,k)= dvar2(i,j,k) + a_(l)*var2(i+l,j,k)
                dvar3(i,j,k)= dvar3(i,j,k) + a_(l)*var3(i+l,j,k)
             enddo
          enddo
       enddo
    enddo

  end subroutine derivative_ksi_v
  
end submodule smod_grid_metrics_deriv_ksi_c3
