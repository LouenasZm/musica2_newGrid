!=================================================================================
module mod_grid_metrics_c3
!=================================================================================
  !> Module to compute full 3D curvilinear metrics
  !> with Geometric Conservation Law (GCL)
!=================================================================================
  use mod_coeff_deriv
  use mod_grid
  use mod_bc
  use mod_constant ! <- for is_SBP ** TO BE CHANGED in mod_coeff_deriv ??
  use warnstop
  implicit none
  !-------------------------------------------------------------------------------
  ! intermediate arrays for 3-D curvilinear metrics
  ! [Rq: used (and deallocated) in grid_normals.f90]
  real(wp), dimension(:,:,:), allocatable :: xdksi,ydksi,zdksi
  real(wp), dimension(:,:,:), allocatable :: xdeta,ydeta,zdeta
  real(wp), dimension(:,:,:), allocatable :: xdphi,ydphi,zdphi
  !-------------------------------------------------------------------------------

  interface
     !===============================================================================
     module subroutine grid_metrics_3d
     end subroutine grid_metrics_3d
     !===============================================================================
     module subroutine correct_deriv_sign(dksi,deta,dphi)
       use mod_grid ! for nx1,nx2,ny1,ny2,nz1,nz2
       real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(inout) :: dksi,deta,dphi
     end subroutine correct_deriv_sign
     !===============================================================================
     module subroutine correct_deriv_sign_v(dksi,deta,dphi)
       use mod_grid ! for nx1_v,nx2_v,ny1_v,ny2_v,nz1_v,nz2_v
       real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v), intent(inout) :: dksi,deta,dphi
     end subroutine correct_deriv_sign_v
     !===============================================================================
     module subroutine correct_dderiv_sign(dksi1,deta1,dphi1,dksi2,deta2,dphi2)
       use mod_grid ! for nx1,nx2,ny1,ny2,nz1,nz2
       real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(inout) :: dksi1,deta1,dphi1
       real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(inout) :: dksi2,deta2,dphi2
     end subroutine correct_dderiv_sign
     !===============================================================================
     module subroutine correct_dderiv_sign_v(dksi1,deta1,dphi1,dksi2,deta2,dphi2)
       use mod_grid ! for nx1_v,nx2_v,ny1_v,ny2_v,nz1_v,nz2_v
       real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v), intent(inout) :: dksi1,deta1,dphi1
       real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v), intent(inout) :: dksi2,deta2,dphi2
     end subroutine correct_dderiv_sign_v
     !===============================================================================
     module subroutine correct_edges_i
     end subroutine correct_edges_i
     !===============================================================================
     module subroutine correct_edges_j
     end subroutine correct_edges_j
     !===============================================================================
     module subroutine grid_metrics_ijacob_3d
     end subroutine grid_metrics_ijacob_3d
     !===============================================================================
     module subroutine grid_metrics_gradients_3d
     end subroutine grid_metrics_gradients_3d
     !===============================================================================
     module subroutine grid_metrics_ijacob_3d_v
     end subroutine grid_metrics_ijacob_3d_v
     !===============================================================================
     module subroutine grid_metrics_ijacob_3d2_v
     end subroutine grid_metrics_ijacob_3d2_v
     !===============================================================================
     module subroutine derivative_ksi(var1,dvar1,var2,dvar2,var3,dvar3)
       real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(in) :: var1,var2,var3
       real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(inout) :: dvar1,dvar2,dvar3
     end subroutine derivative_ksi
     !===============================================================================
     module subroutine derivative_eta(var1,dvar1,var2,dvar2,var3,dvar3)
       real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(in) :: var1,var2,var3
       real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(inout) :: dvar1,dvar2,dvar3
     end subroutine derivative_eta
     !===============================================================================
     module subroutine derivative_phi(var1,dvar1,var2,dvar2,var3,dvar3)
       real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(in) :: var1,var2,var3
       real(wp), dimension(nx1:nx2,ny1:ny2,nz1:nz2), intent(inout) :: dvar1,dvar2,dvar3
     end subroutine derivative_phi
     !===============================================================================
     module subroutine derivative_ksi_v(var1,dvar1,var2,dvar2,var3,dvar3)
       real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v), intent(in) :: var1,var2,var3
       real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v), intent(inout) :: dvar1,dvar2,dvar3
     end subroutine derivative_ksi_v
     !===============================================================================
     module subroutine derivative_eta_v(var1,dvar1,var2,dvar2,var3,dvar3)
       real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v), intent(in) :: var1,var2,var3
       real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v), intent(inout) :: dvar1,dvar2,dvar3
     end subroutine derivative_eta_v
     !===============================================================================
     module subroutine derivative_phi_v(var1,dvar1,var2,dvar2,var3,dvar3)
       real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v), intent(in) :: var1,var2,var3
       real(wp), dimension(nx1_v:nx2_v,ny1_v:ny2_v,nz1_v:nz2_v), intent(inout) :: dvar1,dvar2,dvar3
     end subroutine derivative_phi_v
     !===============================================================================
  end interface

end module mod_grid_metrics_c3
