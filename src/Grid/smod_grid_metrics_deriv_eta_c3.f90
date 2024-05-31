!=================================================================================
submodule (mod_grid_metrics_c3) smod_grid_metrics_deriv_eta_c3
!=================================================================================
  !> Module to compute full 3D curvilinear metrics
  !> with Geometric Conservation Law (GCL)
!=================================================================================

contains

  !===============================================================================
  module subroutine derivative_eta(var1,dvar1,var2,dvar2,var3,dvar3)
  !===============================================================================
    !> Compute derivatives along eta for metrics (full 3D version)
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

    ! Metrics for BC at jmin
    ! ======================
    
    ! For wall & charac, inlet and outlet BCs
    ! ---------------------------------------
    if ((BC_face(2,1)%sort==0).or.(BC_face(2,1)%sort<=-3)) then

       ! SBP boundary schemes
       ! --------------------
       if (is_SBP) then
          if (ngh<4) then
             call mpistop('3D SBP metrics need yet at least 9pt-stencil !', 0)
          else
             ! point #1: stencil [o x x x]
             !j=1
             do k=ndzt,nfzt
                do i=ndxt,nfxt
                   dvar1(i,1,k)= as4p0(1)*var1(i,1,k)+as4p0(2)*var1(i,2,k) &
                               + as4p0(3)*var1(i,3,k)+as4p0(4)*var1(i,4,k) &
                               + dvar1(i,1,k)
                   dvar2(i,1,k)= as4p0(1)*var2(i,1,k)+as4p0(2)*var2(i,2,k) &
                               + as4p0(3)*var2(i,3,k)+as4p0(4)*var2(i,4,k) &
                               + dvar2(i,1,k)
                   dvar3(i,1,k)= as4p0(1)*var3(i,1,k)+as4p0(2)*var3(i,2,k) &
                               + as4p0(3)*var3(i,3,k)+as4p0(4)*var3(i,4,k) &
                               + dvar3(i,1,k)
                enddo
             enddo
             ! point #2: stencil [x o x x]
             !j=2
             do k=ndzt,nfzt
                do i=ndxt,nfxt
                   dvar1(i,2,k)= as4p1(1)*var1(i,1,k)+as4p1(2)*var1(i,2,k) &
                               + as4p1(3)*var1(i,3,k)+as4p1(4)*var1(i,4,k) &
                               + dvar1(i,2,k)
                   dvar2(i,2,k)= as4p1(1)*var2(i,1,k)+as4p1(2)*var2(i,2,k) &
                               + as4p1(3)*var2(i,3,k)+as4p1(4)*var2(i,4,k) &
                               + dvar2(i,2,k)
                   dvar3(i,2,k)= as4p1(1)*var3(i,1,k)+as4p1(2)*var3(i,2,k) &
                               + as4p1(3)*var3(i,3,k)+as4p1(4)*var3(i,4,k) &
                               + dvar3(i,2,k)
                enddo
             enddo
             ! point #3: stencil [x x o x x]
             !j=3
             do k=ndzt,nfzt
                do i=ndxt,nfxt
                   dvar1(i,3,k)= as4p2(1)*var1(i,1,k)+as4p2(2)*var1(i,2,k) &
                               + as4p2(3)*var1(i,3,k)+as4p2(4)*var1(i,4,k) &
                               + as4p2(5)*var1(i,5,k) + dvar1(i,3,k)
                   dvar2(i,3,k)= as4p2(1)*var2(i,1,k)+as4p2(2)*var2(i,2,k) &
                               + as4p2(3)*var2(i,3,k)+as4p2(4)*var2(i,4,k) &
                               + as4p2(5)*var2(i,5,k) + dvar2(i,3,k)
                   dvar3(i,3,k)= as4p2(1)*var3(i,1,k)+as4p2(2)*var3(i,2,k) &
                               + as4p2(3)*var3(i,3,k)+as4p2(4)*var3(i,4,k) &
                               + as4p2(5)*var3(i,5,k) + dvar3(i,3,k)
                enddo
             enddo
             ! point #4: stencil [x x x o x x]
             !j=4
             do k=ndzt,nfzt
                do i=ndxt,nfxt
                   dvar1(i,4,k)= as4p3(1)*var1(i,1,k)+as4p3(2)*var1(i,2,k) &
                               + as4p3(3)*var1(i,3,k)+as4p3(4)*var1(i,4,k) &
                               + as4p3(5)*var1(i,5,k)+as4p3(6)*var1(i,6,k) &
                               + dvar1(i,4,k)
                   dvar2(i,4,k)= as4p3(1)*var2(i,1,k)+as4p3(2)*var2(i,2,k) &
                               + as4p3(3)*var2(i,3,k)+as4p3(4)*var2(i,4,k) &
                               + as4p3(5)*var2(i,5,k)+as4p3(6)*var2(i,6,k) &
                               + dvar2(i,4,k)
                   dvar3(i,4,k)= as4p3(1)*var3(i,1,k)+as4p3(2)*var3(i,2,k) &
                               + as4p3(3)*var3(i,3,k)+as4p3(4)*var3(i,4,k) &
                               + as4p3(5)*var3(i,5,k)+as4p3(6)*var3(i,6,k) &
                               + dvar3(i,4,k)
                enddo
             enddo
          endif
          if (ngh>=5) then
             ! eighth-order (9-pt stencil) for point j=5
             j=5
             do k=ndzt,nfzt
                do i=ndxt,nfxt
                   dvar1(i,j,k)= a9(1)*(var1(i,j+1,k)-var1(i,j-1,k)) &
                               + a9(2)*(var1(i,j+2,k)-var1(i,j-2,k)) &
                               + a9(3)*(var1(i,j+3,k)-var1(i,j-3,k)) &
                               + a9(4)*(var1(i,j+4,k)-var1(i,j-4,k)) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= a9(1)*(var2(i,j+1,k)-var2(i,j-1,k)) &
                               + a9(2)*(var2(i,j+2,k)-var2(i,j-2,k)) &
                               + a9(3)*(var2(i,j+3,k)-var2(i,j-3,k)) &
                               + a9(4)*(var2(i,j+4,k)-var2(i,j-4,k)) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= a9(1)*(var3(i,j+1,k)-var3(i,j-1,k)) &
                               + a9(2)*(var3(i,j+2,k)-var3(i,j-2,k)) &
                               + a9(3)*(var3(i,j+3,k)-var3(i,j-3,k)) &
                               + a9(4)*(var3(i,j+4,k)-var3(i,j-4,k)) &
                               + dvar3(i,j,k)
                enddo
             enddo
          endif

       ! Reduced order stencil
       ! ---------------------
       else       
          if (ngh>=1) then
             ! first-order (2-pt stencil) for point j=1
             j=1
             do k=ndzt,nfzt
                do i=ndxt,nfxt
                   !do l=0,1
                   !   dvar1(i,j,k)=dvar1(i,j,k)+a01(1+l)*var1(i,j+l,k)
                   !   dvar2(i,j,k)=dvar2(i,j,k)+a01(1+l)*var2(i,j+l,k)
                   !   dvar3(i,j,k)=dvar3(i,j,k)+a01(1+l)*var3(i,j+l,k)
                   !enddo
                   do l=0,2
                      dvar1(i,j,k)=dvar1(i,j,k)+a02(1+l)*var1(i,j+l,k)
                      dvar2(i,j,k)=dvar2(i,j,k)+a02(1+l)*var2(i,j+l,k)
                      dvar3(i,j,k)=dvar3(i,j,k)+a02(1+l)*var3(i,j+l,k)
                   enddo
                enddo
             enddo
          endif
          if (ngh>=2) then
             ! second-order (3-pt stencil) for point j=2
             j=2
             do k=ndzt,nfzt
                do i=ndxt,nfxt
                   dvar1(i,j,k)= a3(1)*(var1(i,j+1,k)-var1(i,j-1,k)) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= a3(1)*(var2(i,j+1,k)-var2(i,j-1,k)) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= a3(1)*(var3(i,j+1,k)-var3(i,j-1,k)) &
                               + dvar3(i,j,k)
                enddo
             enddo
          endif
          if (ngh>=3) then
             ! fourth-order (5-pt stencil) for point j=3
             j=3
             do k=ndzt,nfzt
                do i=ndxt,nfxt
                   dvar1(i,j,k)= a5(1)*(var1(i,j+1,k)-var1(i,j-1,k)) &
                               + a5(2)*(var1(i,j+2,k)-var1(i,j-2,k)) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= a5(1)*(var2(i,j+1,k)-var2(i,j-1,k)) &
                               + a5(2)*(var2(i,j+2,k)-var2(i,j-2,k)) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= a5(1)*(var3(i,j+1,k)-var3(i,j-1,k)) &
                               + a5(2)*(var3(i,j+2,k)-var3(i,j-2,k)) &
                               + dvar3(i,j,k)
                enddo
             enddo
          endif
          if (ngh>=4) then
             ! sixth-order (7-pt stencil) for point j=4
             j=4
             do k=ndzt,nfzt
                do i=ndxt,nfxt
                   dvar1(i,j,k)= a7(1)*(var1(i,j+1,k)-var1(i,j-1,k)) &
                               + a7(2)*(var1(i,j+2,k)-var1(i,j-2,k)) &
                               + a7(3)*(var1(i,j+3,k)-var1(i,j-3,k)) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= a7(1)*(var2(i,j+1,k)-var2(i,j-1,k)) &
                               + a7(2)*(var2(i,j+2,k)-var2(i,j-2,k)) &
                               + a7(3)*(var2(i,j+3,k)-var2(i,j-3,k)) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= a7(1)*(var3(i,j+1,k)-var3(i,j-1,k)) &
                               + a7(2)*(var3(i,j+2,k)-var3(i,j-2,k)) &
                               + a7(3)*(var3(i,j+3,k)-var3(i,j-3,k)) &
                               + dvar3(i,j,k)
                enddo
             enddo
          endif
          if (ngh==5) then
             ! eighth-order (9-pt stencil) for point j=5
             j=5
             do k=ndzt,nfzt
                do i=ndxt,nfxt
                   dvar1(i,j,k)= a9(1)*(var1(i,j+1,k)-var1(i,j-1,k)) &
                               + a9(2)*(var1(i,j+2,k)-var1(i,j-2,k)) &
                               + a9(3)*(var1(i,j+3,k)-var1(i,j-3,k)) &
                               + a9(4)*(var1(i,j+4,k)-var1(i,j-4,k)) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= a9(1)*(var2(i,j+1,k)-var2(i,j-1,k)) &
                               + a9(2)*(var2(i,j+2,k)-var2(i,j-2,k)) &
                               + a9(3)*(var2(i,j+3,k)-var2(i,j-3,k)) &
                               + a9(4)*(var2(i,j+4,k)-var2(i,j-4,k)) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= a9(1)*(var3(i,j+1,k)-var3(i,j-1,k)) &
                               + a9(2)*(var3(i,j+2,k)-var3(i,j-2,k)) &
                               + a9(3)*(var3(i,j+3,k)-var3(i,j-3,k)) &
                               + a9(4)*(var3(i,j+4,k)-var3(i,j-4,k)) &
                               + dvar3(i,j,k)
                enddo
             enddo
          endif
       endif
    endif

    ! 7-point stencil for Tam&Dong BC
    ! -------------------------------
    if ((BC_face(2,1)%sort==-1).or.(BC_face(2,1)%sort==-2)) then
       if (ngh<3) then
          call mpistop('Tam&Dong BCs need at least 7pt-stencil !', 0)
       else
          j=1
          do k=ndzt,nfzt
             do i=ndxt,nfxt
                do l=0,6
                   dvar1(i,j,k)=dvar1(i,j,k)+a06(1+l)*var1(i,j+l,k)
                   dvar2(i,j,k)=dvar2(i,j,k)+a06(1+l)*var2(i,j+l,k)
                   dvar3(i,j,k)=dvar3(i,j,k)+a06(1+l)*var3(i,j+l,k)
                enddo
             enddo
          enddo
          j=2
          do k=ndzt,nfzt
             do i=ndxt,nfxt
                do l=-1,5
                   dvar1(i,j,k)=dvar1(i,j,k)+a15(2+l)*var1(i,j+l,k)
                   dvar2(i,j,k)=dvar2(i,j,k)+a15(2+l)*var2(i,j+l,k)
                   dvar3(i,j,k)=dvar3(i,j,k)+a15(2+l)*var3(i,j+l,k)
                enddo
             enddo
          enddo
          j=3
          do k=ndzt,nfzt
             do i=ndxt,nfxt
                do l=-2,4
                   dvar1(i,j,k)=dvar1(i,j,k)+a24(3+l)*var1(i,j+l,k)
                   dvar2(i,j,k)=dvar2(i,j,k)+a24(3+l)*var2(i,j+l,k)
                   dvar3(i,j,k)=dvar3(i,j,k)+a24(3+l)*var3(i,j+l,k)
                enddo
             enddo
          enddo
       endif
       if (ngh>=4) then
          ! sixth-order (7-pt stencil) for point j=4
          j=4
          do k=ndzt,nfzt
             do i=ndxt,nfxt
                dvar1(i,j,k)= a7(1)*(var1(i,j+1,k)-var1(i,j-1,k)) &
                            + a7(2)*(var1(i,j+2,k)-var1(i,j-2,k)) &
                            + a7(3)*(var1(i,j+3,k)-var1(i,j-3,k)) &
                            + dvar1(i,j,k)
                dvar2(i,j,k)= a7(1)*(var2(i,j+1,k)-var2(i,j-1,k)) &
                            + a7(2)*(var2(i,j+2,k)-var2(i,j-2,k)) &
                            + a7(3)*(var2(i,j+3,k)-var2(i,j-3,k)) &
                            + dvar2(i,j,k)
                dvar3(i,j,k)= a7(1)*(var3(i,j+1,k)-var3(i,j-1,k)) &
                            + a7(2)*(var3(i,j+2,k)-var3(i,j-2,k)) &
                            + a7(3)*(var3(i,j+3,k)-var3(i,j-3,k)) &
                            + dvar3(i,j,k)
             enddo
          enddo
       endif
       if (ngh==5) then
          ! sixth-order (7-pt stencil) for point j=5
          j=5
          do k=ndzt,nfzt
             do i=ndxt,nfxt
                dvar1(i,j,k)= a7(1)*(var1(i,j+1,k)-var1(i,j-1,k)) &
                            + a7(2)*(var1(i,j+2,k)-var1(i,j-2,k)) &
                            + a7(3)*(var1(i,j+3,k)-var1(i,j-3,k)) &
                            + dvar1(i,j,k)
                dvar2(i,j,k)= a7(1)*(var2(i,j+1,k)-var2(i,j-1,k)) &
                            + a7(2)*(var2(i,j+2,k)-var2(i,j-2,k)) &
                            + a7(3)*(var2(i,j+3,k)-var2(i,j-3,k)) &
                            + dvar2(i,j,k)
                dvar3(i,j,k)= a7(1)*(var3(i,j+1,k)-var3(i,j-1,k)) &
                            + a7(2)*(var3(i,j+2,k)-var3(i,j-2,k)) &
                            + a7(3)*(var3(i,j+3,k)-var3(i,j-3,k)) &
                            + dvar3(i,j,k)
             enddo
          enddo
       endif
    endif

    ! Metrics for BC at jmax
    ! ======================
    
    ! For wall & charac, inlet and outlet BCs
    ! ---------------------------------------
    if ((BC_face(2,2)%sort==0).or.(BC_face(2,2)%sort<=-3)) then
       
       ! SBP boundary schemes
       ! --------------------
       if (is_SBP) then
         if (ngh<4) then
             call mpistop('3D SBP metrics need yet at least 9pt-stencil !', 0)
          else
             ! point #1: stencil [o x x x]
             j=ny
             do k=ndzt,nfzt
                do i=ndxt,nfxt
                   dvar1(i,j,k)= as4m0(1)*var1(i,ny  ,k)+as4m0(2)*var1(i,ny-1,k) &
                               + as4m0(3)*var1(i,ny-2,k)+as4m0(4)*var1(i,ny-3,k) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= as4m0(1)*var2(i,ny  ,k)+as4m0(2)*var2(i,ny-1,k) &
                               + as4m0(3)*var2(i,ny-2,k)+as4m0(4)*var2(i,ny-3,k) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= as4m0(1)*var3(i,ny  ,k)+as4m0(2)*var3(i,ny-1,k) &
                               + as4m0(3)*var3(i,ny-2,k)+as4m0(4)*var3(i,ny-3,k) &
                               + dvar3(i,j,k)
                enddo
             enddo
             ! point #2: stencil [x o x x]
             j=ny-1
             do k=ndzt,nfzt
                do i=ndxt,nfxt
                   dvar1(i,j,k)= as4m1(1)*var1(i,ny  ,k)+as4m1(2)*var1(i,ny-1,k) &
                               + as4m1(3)*var1(i,ny-2,k)+as4m1(4)*var1(i,ny-3,k) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= as4m1(1)*var2(i,ny  ,k)+as4m1(2)*var2(i,ny-1,k) &
                               + as4m1(3)*var2(i,ny-2,k)+as4m1(4)*var2(i,ny-3,k) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= as4m1(1)*var3(i,ny  ,k)+as4m1(2)*var3(i,ny-1,k) &
                               + as4m1(3)*var3(i,ny-2,k)+as4m1(4)*var3(i,ny-3,k) &
                               + dvar3(i,j,k)
                enddo
             enddo
             ! point #3: stencil [x x o x x]
             j=ny-2
             do k=ndzt,nfzt
                do i=ndxt,nfxt
                   dvar1(i,j,k)= as4m2(1)*var1(i,ny  ,k)+as4m2(2)*var1(i,ny-1,k) &
                               + as4m2(3)*var1(i,ny-2,k)+as4m2(4)*var1(i,ny-3,k) &
                               + as4m2(5)*var1(i,ny-4,k) + dvar1(i,j,k)
                   dvar2(i,j,k)= as4m2(1)*var2(i,ny  ,k)+as4m2(2)*var2(i,ny-1,k) &
                               + as4m2(3)*var2(i,ny-2,k)+as4m2(4)*var2(i,ny-3,k) &
                               + as4m2(5)*var2(i,ny-4,k) + dvar2(i,j,k)
                   dvar3(i,j,k)= as4m2(1)*var3(i,ny  ,k)+as4m2(2)*var3(i,ny-1,k) &
                               + as4m2(3)*var3(i,ny-2,k)+as4m2(4)*var3(i,ny-3,k) &
                               + as4m2(5)*var3(i,ny-4,k) + dvar3(i,j,k)
                enddo
             enddo
             ! point #4: stencil [x x x o x x]
             j=ny-3
             do k=ndzt,nfzt
                do i=ndxt,nfxt
                   dvar1(i,j,k)= as4m3(1)*var1(i,ny  ,k)+as4m3(2)*var1(i,ny-1,k) &
                               + as4m3(3)*var1(i,ny-2,k)+as4m3(4)*var1(i,ny-3,k) &
                               + as4m3(5)*var1(i,ny-4,k)+as4m3(6)*var1(i,ny-5,k) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= as4m3(1)*var2(i,ny  ,k)+as4m3(2)*var2(i,ny-1,k) &
                               + as4m3(3)*var2(i,ny-2,k)+as4m3(4)*var2(i,ny-3,k) &
                               + as4m3(5)*var2(i,ny-4,k)+as4m3(6)*var2(i,ny-5,k) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= as4m3(1)*var3(i,ny  ,k)+as4m3(2)*var3(i,ny-1,k) &
                               + as4m3(3)*var3(i,ny-2,k)+as4m3(4)*var3(i,ny-3,k) &
                               + as4m3(5)*var3(i,ny-4,k)+as4m3(6)*var3(i,ny-5,k) &
                               + dvar3(i,j,k)
                enddo
             enddo
          endif
          if (ngh>=5) then
             ! eighth-order (9-pt stencil) for point j=ny-4
             j=ny-4
             do k=ndzt,nfzt
                do i=ndxt,nfxt
                   dvar1(i,j,k)= a9(1)*(var1(i,j+1,k)-var1(i,j-1,k)) &
                               + a9(2)*(var1(i,j+2,k)-var1(i,j-2,k)) &
                               + a9(3)*(var1(i,j+3,k)-var1(i,j-3,k)) &
                               + a9(4)*(var1(i,j+4,k)-var1(i,j-4,k)) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= a9(1)*(var2(i,j+1,k)-var2(i,j-1,k)) &
                               + a9(2)*(var2(i,j+2,k)-var2(i,j-2,k)) &
                               + a9(3)*(var2(i,j+3,k)-var2(i,j-3,k)) &
                               + a9(4)*(var2(i,j+4,k)-var2(i,j-4,k)) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= a9(1)*(var3(i,j+1,k)-var3(i,j-1,k)) &
                               + a9(2)*(var3(i,j+2,k)-var3(i,j-2,k)) &
                               + a9(3)*(var3(i,j+3,k)-var3(i,j-3,k)) &
                               + a9(4)*(var3(i,j+4,k)-var3(i,j-4,k)) &
                               + dvar3(i,j,k)
                enddo
             enddo
          endif

       ! Reduced order stencil
       ! ---------------------
       else       
          if (ngh>=1) then
             ! first-order (2-pt stencil) for point j=ny
             j=ny
             do k=ndzt,nfzt
                do i=ndxt,nfxt
                   !do l=-1,0
                   !   dvar1(i,j,k)=dvar1(i,j,k)-a01(1-l)*var1(i,j+l,k)
                   !   dvar2(i,j,k)=dvar2(i,j,k)-a01(1-l)*var2(i,j+l,k)
                   !   dvar3(i,j,k)=dvar3(i,j,k)-a01(1-l)*var3(i,j+l,k)
                   !enddo
                   do l=-2,0
                      dvar1(i,j,k)=dvar1(i,j,k)-a02(1-l)*var1(i,j+l,k)
                      dvar2(i,j,k)=dvar2(i,j,k)-a02(1-l)*var2(i,j+l,k)
                      dvar3(i,j,k)=dvar3(i,j,k)-a02(1-l)*var3(i,j+l,k)
                   enddo
                enddo
             enddo
          endif
          if (ngh>=2) then
             ! second-order (3-pt stencil) for point j=ny-1
             j=ny-1
             do k=ndzt,nfzt
                do i=ndxt,nfxt
                   dvar1(i,j,k)= a3(1)*(var1(i,j+1,k)-var1(i,j-1,k)) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= a3(1)*(var2(i,j+1,k)-var2(i,j-1,k)) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= a3(1)*(var3(i,j+1,k)-var3(i,j-1,k)) &
                               + dvar3(i,j,k)
                enddo
             enddo
          endif
          if (ngh>=3) then
             ! fourth-order (5-pt stencil) for point j=ny-2
             j=ny-2
             do k=ndzt,nfzt
                do i=ndxt,nfxt
                   dvar1(i,j,k)= a5(1)*(var1(i,j+1,k)-var1(i,j-1,k)) &
                               + a5(2)*(var1(i,j+2,k)-var1(i,j-2,k)) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= a5(1)*(var2(i,j+1,k)-var2(i,j-1,k)) &
                               + a5(2)*(var2(i,j+2,k)-var2(i,j-2,k)) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= a5(1)*(var3(i,j+1,k)-var3(i,j-1,k)) &
                               + a5(2)*(var3(i,j+2,k)-var3(i,j-2,k)) &
                               + dvar3(i,j,k)
                enddo
             enddo
          endif
          if (ngh>=4) then
             ! sixth-order (7-pt stencil) for point j=ny-3
             j=ny-3
             do k=ndzt,nfzt
                do i=ndxt,nfxt
                   dvar1(i,j,k)= a7(1)*(var1(i,j+1,k)-var1(i,j-1,k)) &
                               + a7(2)*(var1(i,j+2,k)-var1(i,j-2,k)) &
                               + a7(3)*(var1(i,j+3,k)-var1(i,j-3,k)) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= a7(1)*(var2(i,j+1,k)-var2(i,j-1,k)) &
                               + a7(2)*(var2(i,j+2,k)-var2(i,j-2,k)) &
                               + a7(3)*(var2(i,j+3,k)-var2(i,j-3,k)) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= a7(1)*(var3(i,j+1,k)-var3(i,j-1,k)) &
                               + a7(2)*(var3(i,j+2,k)-var3(i,j-2,k)) &
                               + a7(3)*(var3(i,j+3,k)-var3(i,j-3,k)) &
                               + dvar3(i,j,k)
                enddo
             enddo
          endif
          if (ngh==5) then
             ! eighth-order (9-pt stencil) for point j=ny-4
             j=ny-4
             do k=ndzt,nfzt
                do i=ndxt,nfxt
                   dvar1(i,j,k)= a9(1)*(var1(i,j+1,k)-var1(i,j-1,k)) &
                               + a9(2)*(var1(i,j+2,k)-var1(i,j-2,k)) &
                               + a9(3)*(var1(i,j+3,k)-var1(i,j-3,k)) &
                               + a9(4)*(var1(i,j+4,k)-var1(i,j-4,k)) &
                               + dvar1(i,j,k)
                   dvar2(i,j,k)= a9(1)*(var2(i,j+1,k)-var2(i,j-1,k)) &
                               + a9(2)*(var2(i,j+2,k)-var2(i,j-2,k)) &
                               + a9(3)*(var2(i,j+3,k)-var2(i,j-3,k)) &
                               + a9(4)*(var2(i,j+4,k)-var2(i,j-4,k)) &
                               + dvar2(i,j,k)
                   dvar3(i,j,k)= a9(1)*(var3(i,j+1,k)-var3(i,j-1,k)) &
                               + a9(2)*(var3(i,j+2,k)-var3(i,j-2,k)) &
                               + a9(3)*(var3(i,j+3,k)-var3(i,j-3,k)) &
                               + a9(4)*(var3(i,j+4,k)-var3(i,j-4,k)) &
                               + dvar3(i,j,k)
                enddo
             enddo
          endif
       endif
    endif

    ! 7-point stencil for Tam&Dong BC
    ! -------------------------------
    if ((BC_face(2,2)%sort==-1).or.(BC_face(2,2)%sort==-2)) then
       if (ngh<3) then
          call mpistop('Tam&Dong BCs need at least 7pt-stencil !', 0)
       else
          j=ny
          do k=ndzt,nfzt
             do i=ndxt,nfxt
                do l=-6,0
                   dvar1(i,j,k)=dvar1(i,j,k)-a06(1-l)*var1(i,j+l,k)
                   dvar2(i,j,k)=dvar2(i,j,k)-a06(1-l)*var2(i,j+l,k)
                   dvar3(i,j,k)=dvar3(i,j,k)-a06(1-l)*var3(i,j+l,k)
                enddo
             enddo
          enddo
          j=ny-1
          do k=ndzt,nfzt
             do i=ndxt,nfxt
                do l=-5,1
                   dvar1(i,j,k)=dvar1(i,j,k)-a15(2-l)*var1(i,j+l,k)
                   dvar2(i,j,k)=dvar2(i,j,k)-a15(2-l)*var2(i,j+l,k)
                   dvar3(i,j,k)=dvar3(i,j,k)-a15(2-l)*var3(i,j+l,k)
                enddo
             enddo
          enddo
          j=ny-2
          do k=ndzt,nfzt
             do i=ndxt,nfxt
                do l=-4,2
                   dvar1(i,j,k)=dvar1(i,j,k)-a24(3-l)*var1(i,j+l,k)
                   dvar2(i,j,k)=dvar2(i,j,k)-a24(3-l)*var2(i,j+l,k)
                   dvar3(i,j,k)=dvar3(i,j,k)-a24(3-l)*var3(i,j+l,k)
                enddo
             enddo
          enddo
       endif
       if (ngh>=4) then
          ! sixth-order (7-pt stencil) for point j=ny-3
          j=ny-3
          do k=ndzt,nfzt
             do i=ndxt,nfxt
                dvar1(i,j,k)= a7(1)*(var1(i,j+1,k)-var1(i,j-1,k)) &
                            + a7(2)*(var1(i,j+2,k)-var1(i,j-2,k)) &
                            + a7(3)*(var1(i,j+3,k)-var1(i,j-3,k)) &
                            + dvar1(i,j,k)
                dvar2(i,j,k)= a7(1)*(var2(i,j+1,k)-var2(i,j-1,k)) &
                            + a7(2)*(var2(i,j+2,k)-var2(i,j-2,k)) &
                            + a7(3)*(var2(i,j+3,k)-var2(i,j-3,k)) &
                            + dvar2(i,j,k)
                dvar3(i,j,k)= a7(1)*(var3(i,j+1,k)-var3(i,j-1,k)) &
                            + a7(2)*(var3(i,j+2,k)-var3(i,j-2,k)) &
                            + a7(3)*(var3(i,j+3,k)-var3(i,j-3,k)) &
                            + dvar3(i,j,k)
             enddo
          enddo
       endif
       if (ngh==5) then
          ! sixth-order (7-pt stencil) for point j=ny-4
          j=ny-4
          do k=ndzt,nfzt
             do i=ndxt,nfxt
                dvar1(i,j,k)= a7(1)*(var1(i,j+1,k)-var1(i,j-1,k)) &
                            + a7(2)*(var1(i,j+2,k)-var1(i,j-2,k)) &
                            + a7(3)*(var1(i,j+3,k)-var1(i,j-3,k)) &
                            + dvar1(i,j,k)
                dvar2(i,j,k)= a7(1)*(var2(i,j+1,k)-var2(i,j-1,k)) &
                            + a7(2)*(var2(i,j+2,k)-var2(i,j-2,k)) &
                            + a7(3)*(var2(i,j+3,k)-var2(i,j-3,k)) &
                            + dvar2(i,j,k)
                dvar3(i,j,k)= a7(1)*(var3(i,j+1,k)-var3(i,j-1,k)) &
                            + a7(2)*(var3(i,j+2,k)-var3(i,j-2,k)) &
                            + a7(3)*(var3(i,j+3,k)-var3(i,j-3,k)) &
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
       do j=ndy,nfy
          do i=ndxt,nfxt
             do l=-ngh,ngh
                dvar1(i,j,k)= dvar1(i,j,k) + a_(l)*var1(i,j+l,k)
                dvar2(i,j,k)= dvar2(i,j,k) + a_(l)*var2(i,j+l,k)
                dvar3(i,j,k)= dvar3(i,j,k) + a_(l)*var3(i,j+l,k)
             enddo
          enddo
       enddo
    enddo

  end subroutine derivative_eta
  
  !===============================================================================
  module subroutine derivative_eta_v(var1,dvar1,var2,dvar2,var3,dvar3)
  !===============================================================================
    !> Compute derivatives along eta for metrics (full 3D version)
    !> - stencil -ngh_v:+ngh_v for viscous fluxes -
  !===============================================================================
    implicit none
    ! ---------------------------------------------------------------------------
    real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v), intent(in) :: var1,var2,var3
    real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v), intent(inout) :: dvar1,dvar2,dvar3
    ! ---------------------------------------------------------------------------
    integer :: i,j,k,l
    integer :: ndyv,nfyv
    real(wp), dimension(-ngh_v:ngh_v) :: a_ ! coefficients of interior FD scheme
    ! ---------------------------------------------------------------------------

    ! Metrics for BC at jmin
    ! ======================
    ndyv=1
    
    if (BC_face(2,1)%sort<=0) then
       select case (ngh_v)
       case (2)
          ndyv=3
          j=1
          do k=ndz_v1,nfz_v1
             do i=ndx_v1,nfx_v1
                do l=0,4
                   dvar1(i,j,k)=dvar1(i,j,k)+a04(1+l)*var1(i,j+l,k)
                   dvar2(i,j,k)=dvar2(i,j,k)+a04(1+l)*var2(i,j+l,k)
                   dvar3(i,j,k)=dvar3(i,j,k)+a04(1+l)*var3(i,j+l,k)
                enddo
             enddo
          enddo
          j=2
          do k=ndz_v1,nfz_v1
             do i=ndx_v1,nfx_v1
                do l=-1,3
                   dvar1(i,j,k)=dvar1(i,j,k)+a13(2+l)*var1(i,j+l,k)
                   dvar2(i,j,k)=dvar2(i,j,k)+a13(2+l)*var2(i,j+l,k)
                   dvar3(i,j,k)=dvar3(i,j,k)+a13(2+l)*var3(i,j+l,k)
                enddo
             enddo
          enddo
       case (1)
          ndyv=2
          j=1
          do k=ndz_v1,nfz_v1
             do i=ndx_v1,nfx_v1
                do l=0,1
                   dvar1(i,j,k)=dvar1(i,j,k)+a01(1+l)*var1(i,j+l,k)
                   dvar2(i,j,k)=dvar2(i,j,k)+a01(1+l)*var2(i,j+l,k)
                   dvar3(i,j,k)=dvar3(i,j,k)+a01(1+l)*var3(i,j+l,k)
                enddo
             enddo
          enddo
       case default
          call mpistop('max ngh_v is 2 (5pts-stencil)!', 0)
       end select
    endif

    ! Metrics for BC at jmax
    ! ======================
    nfyv=ny
    
    if (BC_face(2,2)%sort<=0) then
       select case (ngh_v)
       case (2)
          nfyv=ny-2
          j=ny-1
          do k=ndz_v1,nfz_v1
             do i=ndx_v1,nfx_v1
                do l=-3,1
                   dvar1(i,j,k)=dvar1(i,j,k)-a13(2-l)*var1(i,j+l,k)
                   dvar2(i,j,k)=dvar2(i,j,k)-a13(2-l)*var2(i,j+l,k)
                   dvar3(i,j,k)=dvar3(i,j,k)-a13(2-l)*var3(i,j+l,k)
                enddo
             enddo
          enddo
          j=ny
          do k=ndz_v1,nfz_v1
             do i=ndx_v1,nfx_v1
                do l=-4,0
                   dvar1(i,j,k)=dvar1(i,j,k)-a04(1-l)*var1(i,j+l,k)
                   dvar2(i,j,k)=dvar2(i,j,k)-a04(1-l)*var2(i,j+l,k)
                   dvar3(i,j,k)=dvar3(i,j,k)-a04(1-l)*var3(i,j+l,k)
                enddo
             enddo
          enddo
       case (1)
          nfyv=ny-1
          j=ny
          do k=ndz_v1,nfz_v1
             do i=ndx_v1,nfx_v1
                do l=-1,0
                   dvar1(i,j,k)=dvar1(i,j,k)-a01(1-l)*var1(i,j+l,k)
                   dvar2(i,j,k)=dvar2(i,j,k)-a01(1-l)*var2(i,j+l,k)
                   dvar3(i,j,k)=dvar3(i,j,k)-a01(1-l)*var3(i,j+l,k)
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
       do j=ndyv,nfyv
          do i=ndx_v1,nfx_v1
             do l=-ngh_v,ngh_v
                dvar1(i,j,k)= dvar1(i,j,k) + a_(l)*var1(i,j+l,k)
                dvar2(i,j,k)= dvar2(i,j,k) + a_(l)*var2(i,j+l,k)
                dvar3(i,j,k)= dvar3(i,j,k) + a_(l)*var3(i,j+l,k)
             enddo
          enddo
       enddo
    enddo   

  end subroutine derivative_eta_v

end submodule smod_grid_metrics_deriv_eta_c3
