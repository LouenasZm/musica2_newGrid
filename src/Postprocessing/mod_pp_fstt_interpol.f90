!==============================================================================
module mod_pp_fstt_interpol
!==============================================================================
  !> Module for interpolation for FSTT PP
!==============================================================================
  use mod_constant
  use mod_flow
  use mod_pp_mpi
  use mod_pp_var
  use mod_interp1
  implicit none
  ! ---------------------------------------------------------------------------

contains

  !============================================================================
  subroutine init_fstt_interpol
  !============================================================================
    !> Initialisation of interpolation for FSTT PP
    !============================================================================
    use mod_utils
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: Lx_g,Lz_g ! Length of domain per proc for interpolation
    real(wp) :: step,step3
    real(wp), dimension(:), allocatable :: x_interp_full,z_interp_full
    ! -------------------------------------------------------------------------

    if (iproc==0) write(*,18) nit_interp,nj_interp,nkt_interp
18  format(5x,"Number of points for interpolation: ",i0," x ",i0," x ",i0,"  (i x j x k)")
    if (iproc==0) write(*,19) nil_interp,nj_interp,nk_interp
19  format(5x,"Number of points per proc:          ",i0," x ",i0," x ",i0,"  (i x j x k)")

    ! Length of the domain per proc
    Lz_g = zg(ngz) - zg(1)

    ! Creation of the full interpolated grid in x direction
    allocate(x_interp_full(nit_interp))
    ! if ((bl(1)%is_sponge).and.(bl_glob%ni-50.lt.bl(1)%ni)) then
    !    i = ngx - (bl_glob%ni-50)
    !    Lx_g = xg(ngx-i) - xg(1)
    !    ratio = 1.0_wp*nit_interp/ngx
    !    k = (ngx-i)*ratio
    !    step = Lx_g/(k-0.5_wp)
    !    if (iproc.eq.0) print *,'step',step
    !    call mpistop('',0)
    ! else
    !    Lx_g = xg(ngx) - xg(1)
    !    step = Lx_g/(nit_interp-1)
    !    x_interp_full(1) = xg(1)
    !    do i=2,nit_interp
    !       x_interp_full(i) = x_interp_full(i-1) + step
    !    enddo
    ! endif


    Lx_g = xg(ngx) - xg(1)
    step = Lx_g/(nit_interp-1)
    x_interp_full(1) = xg(1)
    do i=2,nit_interp
       x_interp_full(i) = x_interp_full(i-1) + step
    enddo

    allocate(z_interp_full(nkt_interp))
    step = Lz_g/(nkt_interp-1)
    z_interp_full(1) = zg(1)
    do k=2,nkt_interp
       z_interp_full(k) = z_interp_full(k-1) + step
    enddo

    ! Determination of ghost points and index bounds
    ! ==============================================
    ! Calculated based on the kernel
    prct_kernel = 0.025_wp
    ngh_pp=0
    step = Lx_g/(nit_interp-1)
    if (iproc.eq.0) print *,"h_kern*step ~>",h_kern*step
    do while (exp(-0.5_wp*((ngh_pp*step)**2)**0.5/(h_kern*step))>prct_kernel)
       ngh_pp = ngh_pp+1
    enddo
    do j=2,nj_interp
       born_kern_j(j,1) = j
       do while ((exp(-0.5_wp*((y(ndivy*born_kern_j(j,1))-y(ndivy*j))**2)**0.5/(h_kern*step))>prct_kernel).and.(born_kern_j(j,1).gt.2))
          born_kern_j(j,1) = born_kern_j(j,1) - 1
       enddo
       born_kern_j(j,2) = j
       do while ((exp(-0.5_wp*((y(ndivy*born_kern_j(j,2))-y(ndivy*j))**2)**0.5/(h_kern*step))>prct_kernel).and.(born_kern_j(j,2).lt.nj_interp))
          born_kern_j(j,2) = born_kern_j(j,2) + 1
       enddo
    enddo
    nkernel_k=0
    step3 = Lz_g/(nkt_interp-1)
    do while (exp(-0.5_wp*((nkernel_k*step3)**2)**0.5/(h_kern*step))>prct_kernel)
       nkernel_k = nkernel_k+1
    enddo
    if (iproc==0) write(*,20) ngh_pp,nkernel_k
20  format(5x,"Ghost points for pp: ngh_pp = ",i0,", nkernel_k = ",i0)
    if (iproc==0) write(*,21) (2*ngh_pp+1),born_kern_j(2,2)-born_kern_j(2,1)+1,(2*nkernel_k+1)
21  format(5x,"Number of points for kernel:        ",i0," x ",i0," x ",i0,"  (i x j x k)")

    ! index bounds extended
    ndxt_pp = 1-ngh_pp; nfxt_pp = nil_interp+ngh_pp
    if (is_boundary(1,1)) ndxt_pp=1
    if (is_boundary(1,2)) nfxt_pp=nil_interp
    ndyt_pp = 1; nfyt_pp = nj_interp
    ndzt_pp = 1-nkernel_k; nfzt_pp = nk_interp+nkernel_k
    if (is_boundary(3,1)) ndzt_pp=1
    if (is_boundary(3,2)) nfzt_pp=nk_interp
    ! index bounds interior, without boundary points
    ndx_pp = 1; nfx_pp = nil_interp
    if (is_boundary(1,1)) ndx_pp=2
    if (is_boundary(1,2)) nfx_pp=nil_interp-1
    ndy_pp = 2; nfy_pp = nj_interp-1
    ndz_pp = 1; nfz_pp = nk_interp
    ! index bounds interior, minus ngh_pp at boundary
    ndxtpngh_pp = 1; nfxtmngh_pp = nil_interp
    ! if (is_boundary(1,1)) ndxtpngh_pp = 1 + ngh_pp
    if (is_boundary(1,1)) ndxtpngh_pp = min(1 + ngh_pp,nil_interp)
    if (is_boundary(1,2)) nfxtmngh_pp = max(nil_interp - ngh_pp,1)


    ! Creation of the grid for interpolation
    ! ======================================

    ! Allocation of grid for interpolation
    allocate(x_interp(0:nil_interp+1),y_interp(1:nj_interp),z_interp(0:nk_interp+1))

    ! Creation of local grid for interpolation
    ! i direction
    do i=1,nil_interp
       x_interp(i) = x_interp_full(i + coord(1)*nil_interp)
    enddo
    ! Artificial extension for inverse interpolation
    step = x_interp(2) - x_interp(1)
    x_interp(0) = x_interp(1) - step
    x_interp(nil_interp+1) = x_interp(nil_interp) + step

    ! j direction
    y_interp(1) = y(1)
    do j=2,nj_interp
       y_interp(j) = y(ndivy*j)
    enddo

    ! k direction
    do k=1,nk_interp
       z_interp(k) = z_interp_full(k + coord(3)*nk_interp)
    enddo
    ! Artificial extension for inverse interpolation
    step = z_interp(2) - z_interp(1)
    z_interp(0) = z_interp(1) - step
    z_interp(nk_interp+1) = z_interp(nk_interp) + step

    ! Writting of the interpolated grid
    ! =================================
    if (iproc.eq.0) then
       ! Writting
       open(194,file='grid_interp_bl'//trim(numchar(iblc_pp))//'.bin',form='unformatted',status='unknown')
       rewind(194)
       write(194) nit_interp
       write(194) nj_interp
       write(194) nkt_interp
       write(194) (x_interp_full(i),i=1,nit_interp)
       write(194) (y_interp(j),j=1,nj_interp)
       write(194) (z_interp_full(k),k=1,nkt_interp)
    endif

    deallocate(x_interp_full,z_interp_full)
    ! call mpistop('',0)

  end subroutine init_fstt_interpol

  !============================================================================
  subroutine fstt_interpol_3d(id_case)
  !============================================================================
    !> Trilinear interpolation for FSTT PP
    !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer, intent(in)  :: id_case
    ! -------------------------------------------------------------------------
    integer :: i,j,k,i0,j0,k0,i1,j1,k1,m
    real(wp) :: xd,yd,zd,cc,cc0,cc1,cc00,cc01,cc10,cc11
    ! -------------------------------------------------------------------------

    loopz: do k=1,nk_interp
      loopz1: do k1=1,nz
         if (z(k1).gt.z_interp(k)) exit loopz1
      enddo loopz1

      if (k1.gt.nz) then
         k1 = nz
         k0 = nz
         zd = 0.0_wp
      else
         k0 = k1-1
         zd = (z_interp(k)-z(k0))/(z(k1)-z(k0))
      endif
      if (zd.lt.0_wp  ) zd = 0.0_wp
      if (zd.gt.1.0_wp) zd = 1.0_wp

      do j=1,nj_interp
         loopy1: do j1=1,ny
            if (y(j1).gt.y_interp(j)) exit loopy1
         enddo loopy1

         if (j1.gt.ny) then
            j1 = ny
            j0 = ny
            yd = 0.0_wp
         else
            j0 = j1-1
            yd = (y_interp(j)-y(j0))/(y(j1)-y(j0))
         endif
         if (yd.lt.0_wp  ) yd = 0.0_wp
         if (yd.gt.1.0_wp) yd = 1.0_wp

         do i=1,nil_interp
            loopx1: do i1=ndx,nfx
               if (x(i1).gt.x_interp(i)) exit loopx1
            enddo loopx1

            if (i1.gt.nfx) then
               i1 = nfx
               i0 = nfx
               xd = 0.0_wp
            else if (i1.eq.ndx) then
               i0 = ndx
               xd = 0.0_wp
            else
               i0 = i1-1
               xd = (x_interp(i)-x(i0))/(x(i1)-x(i0))
            endif
            if (xd.lt.0_wp  ) xd = 0.0_wp
            if (xd.gt.1.0_wp) xd = 1.0_wp

            select case(id_case)
            case(1)
               ! Interpolation of stats
               ! ======================
               do m=1,23

                  cc0 = stats_proc(i0,j0,m)*(1.0_wp-xd) + stats_proc(i1,j0,m)*xd
                  cc1 = stats_proc(i0,j1,m)*(1.0_wp-xd) + stats_proc(i1,j1,m)*xd
                  !
                  cc = cc0*(1.0_wp-yd) + cc1*yd

                  stats_interp(i,j,m) = cc
               enddo
            case(2)
               ! Interpolation of u
               ! ==================
               cc00 = uu(i0,j0,k0)*(1.0_wp-xd) + uu(i1,j0,k0)*xd
               cc01 = uu(i0,j0,k1)*(1.0_wp-xd) + uu(i1,j0,k1)*xd
               cc10 = uu(i0,j1,k0)*(1.0_wp-xd) + uu(i1,j1,k0)*xd
               cc11 = uu(i0,j1,k1)*(1.0_wp-xd) + uu(i1,j1,k1)*xd
               !
               cc0 = cc00*(1.0_wp-yd) + cc10*yd
               cc1 = cc01*(1.0_wp-yd) + cc11*yd
               !
               cc  =  cc0*(1.0_wp-zd) +  cc1*zd

               uu_interp(i,j,k) = cc

               ! Interpolation of v
               ! ==================
               cc00 = vv(i0,j0,k0)*(1.0_wp-xd) + vv(i1,j0,k0)*xd
               cc01 = vv(i0,j0,k1)*(1.0_wp-xd) + vv(i1,j0,k1)*xd
               cc10 = vv(i0,j1,k0)*(1.0_wp-xd) + vv(i1,j1,k0)*xd
               cc11 = vv(i0,j1,k1)*(1.0_wp-xd) + vv(i1,j1,k1)*xd
               !
               cc0 = cc00*(1.0_wp-yd) + cc10*yd
               cc1 = cc01*(1.0_wp-yd) + cc11*yd
               !
               cc  =  cc0*(1.0_wp-zd) +  cc1*zd

               vv_interp(i,j,k) = cc

               ! Interpolation of w
               ! ==================
               cc00 = ww(i0,j0,k0)*(1.0_wp-xd) + ww(i1,j0,k0)*xd
               cc01 = ww(i0,j0,k1)*(1.0_wp-xd) + ww(i1,j0,k1)*xd
               cc10 = ww(i0,j1,k0)*(1.0_wp-xd) + ww(i1,j1,k0)*xd
               cc11 = ww(i0,j1,k1)*(1.0_wp-xd) + ww(i1,j1,k1)*xd
               !
               cc0 = cc00*(1.0_wp-yd) + cc10*yd
               cc1 = cc01*(1.0_wp-yd) + cc11*yd
               !
               cc  =  cc0*(1.0_wp-zd) +  cc1*zd

               ww_interp(i,j,k) = cc
            case default
               call mpistop('bad choice for interpolation', 0)
            end select
         enddo
      enddo
      if (id_case.eq.1) exit loopz
   enddo loopz

  end subroutine fstt_interpol_3d

  !============================================================================
  subroutine fstt_interpol_3d_inv
  !============================================================================
    !> Trilinear interpolation for FSTT PP
    !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i,j,k,i0,j0,k0,i1,j1,k1
    real(wp) :: xd,yd,zd,cc,cc0,cc1,cc00,cc01,cc10,cc11
    ! -------------------------------------------------------------------------

    do k=1,nz
      loopz1: do k1=1,nk_interp
         if (z_interp(k1).gt.z(k)) exit loopz1
      enddo loopz1

      if (k1.gt.nk_interp) then
         k1 = nk_interp+1
         k0 = nk_interp
         zd = (z(k)-z_interp(k0))/(z_interp(k1)-z_interp(k0))
      else
         k0 = k1-1
         zd = (z(k)-z_interp(k0))/(z_interp(k1)-z_interp(k0))
      endif
      if (zd.lt.0.0_wp  ) zd = 0.0_wp
      if (zd.gt.1.0_wp) zd = 1.0_wp

      do j=1,ny
         loopy1: do j1=1,nj_interp
            if (y_interp(j1).gt.y(j)) exit loopy1
         enddo loopy1

         if (j1.gt.nj_interp) then
            j1 = nj_interp
            j0 = nj_interp
            yd = 0.0_wp
         else
            j0 = j1-1
            yd = (y(j)-y_interp(j0))/(y_interp(j1)-y_interp(j0))
         endif
         if (yd.lt.0_wp  ) yd = 0.0_wp
         if (yd.gt.1.0_wp) yd = 1.0_wp

         do i=1,nx
            loopx1: do i1=1,nil_interp
               if (x_interp(i1).gt.x(i)) exit loopx1
            enddo loopx1

            if (i1.gt.nil_interp) then
               i1 = nil_interp+1
               i0 = nil_interp
               xd = (x(i)-x_interp(i0))/(x_interp(i1)-x_interp(i0))
             else
               i0 = i1-1
               xd = (x(i)-x_interp(i0))/(x_interp(i1)-x_interp(i0))
            endif
            if (xd.lt.0_wp  ) xd = 0.0_wp
            if (xd.gt.1.0_wp) xd = 1.0_wp

            ! Interpolation of extr_density
            ! =============================
            cc00 = extr_density(i0,j0,k0)*(1.0_wp-xd) + extr_density(i1,j0,k0)*xd
            cc01 = extr_density(i0,j0,k1)*(1.0_wp-xd) + extr_density(i1,j0,k1)*xd
            cc10 = extr_density(i0,j1,k0)*(1.0_wp-xd) + extr_density(i1,j1,k0)*xd
            cc11 = extr_density(i0,j1,k1)*(1.0_wp-xd) + extr_density(i1,j1,k1)*xd
            !
            cc0 = cc00*(1.0_wp-yd) + cc10*yd
            cc1 = cc01*(1.0_wp-yd) + cc11*yd
            !
            cc  =  cc0*(1.0_wp-zd) +  cc1*zd

            extr_density2(i,j,k) = cc
         enddo
      enddo
   enddo

  end subroutine fstt_interpol_3d_inv

  !============================================================================
  subroutine fstt_interpol_xy
  !============================================================================
    !> Interpolation for FSTT PP for xy plane
    !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    ! -------------------------------------------------------------------------

  end subroutine fstt_interpol_xy

  end module mod_pp_fstt_interpol
