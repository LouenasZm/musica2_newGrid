!===============================================================================
submodule (mod_init_TamDong) smod_init_TamDong2
!===============================================================================
  !> author: XG
  !> date: January 2024
  !> Initialization of Tam and Dong's boundary conditions
  !> routines: init_param_TD2d, init_param_TD2d_1pt, free_bc_TD2d
!===============================================================================
  implicit none

contains

  !=============================================================================
  module subroutine init_param_TD2d
  !=============================================================================
    !> Preprocessing for Tam & Dong's boundary conditions
    !> allocations and compute geometrical parameters (Att. Radiation center def)
  !=============================================================================
    implicit none
    !---------------------------------------------------------------------------
    integer :: i,j,l,m
    real(wp) :: xr,yr,ir
    !---------------------------------------------------------------------------

    ! To be changed
    !!if ((.not.is_vortex).and.(type_vortex.eq.-1)) is_mean0=.false.
    is_mean0=.false.

    ! imin (left) boundary condition
    ! ==============================
    if (is_bc_TD(1,1)) then
       
       ! polar coordinates
       ! -----------------
       allocate(BC_face(1,1)%ir(ngh,ndy:nfy))
       allocate(BC_face(1,1)%cosphi(ngh,ndy:nfy))
       allocate(BC_face(1,1)%sinphi(ngh,ndy:nfy))
       
       do i=1,ngh
          do j=ndy,nfy
             if (is_curv) then
                xr=xc(i,j)-BC_face(1,1)%xrc(j)
                yr=yc(i,j)-BC_face(1,1)%yrc(j)
             else
                xr=x(i)-BC_face(1,1)%xrc(j)
                yr=y(j)-BC_face(1,1)%yrc(j)
             endif
             if (SHIT) yr=0.0_wp
             ir=sqrt(xr**2+yr**2)
             if (ir==0.0_wp)  ir=1.0_wp
             ir=1.0_wp/ir
             BC_face(1,1)%ir(i,j)=ir
             BC_face(1,1)%cosphi(i,j)=xr*ir
             BC_face(1,1)%sinphi(i,j)=yr*ir
          enddo
       enddo
       !if (is_2D) BC_face(1,1)%ir=0.5_wp*BC_face(1,1)%ir
    endif

    ! imax (right) boundary condition
    ! ===============================
    if (is_bc_TD(1,2)) then

       ! polar coordinates
       ! -----------------
       allocate(BC_face(1,2)%ir(ngh,ndy:nfy))
       allocate(BC_face(1,2)%cosphi(ngh,ndy:nfy))
       allocate(BC_face(1,2)%sinphi(ngh,ndy:nfy))

       do i=nx-ngh+1,nx
          l=i-nx+ngh
          do j=ndy,nfy
             if (is_curv) then
                xr=xc(i,j)-BC_face(1,2)%xrc(j)
                yr=yc(i,j)-BC_face(1,2)%yrc(j)
             else
                xr=x(i)-BC_face(1,2)%xrc(j)
                yr=y(j)-BC_face(1,2)%yrc(j)
             endif
             if (SHIT) yr=0.0_wp
             ir=sqrt(xr**2+yr**2)
             if (ir==0.0_wp)  ir=1.0_wp
             ir=1.0_wp/ir
             BC_face(1,2)%ir(l,j)=ir
             BC_face(1,2)%cosphi(l,j)=xr*ir
             BC_face(1,2)%sinphi(l,j)=yr*ir
          enddo
       enddo
       !if (is_2D) BC_face(1,2)%ir=0.5_wp*BC_face(1,2)%ir
    endif

    ! jmin (bottom) boundary condition
    ! ================================
    if (is_bc_TD(2,1)) then

       ! polar coordinates
       ! -----------------
       allocate(BC_face(2,1)%ir(ndx:nfx,ngh))
       allocate(BC_face(2,1)%cosphi(ndx:nfx,ngh))
       allocate(BC_face(2,1)%sinphi(ndx:nfx,ngh))

       do j=1,ngh
          do i=ndx,nfx
             if (is_curv) then
                xr=xc(i,j)-BC_face(2,1)%xrc(i)
                yr=yc(i,j)-BC_face(2,1)%yrc(i)
             else
                xr=x(i)-BC_face(2,1)%xrc(i)
                yr=y(j)-BC_face(2,1)%yrc(i)
             endif
             ir=sqrt(xr**2+yr**2)
             if (ir==0.0_wp)  ir=1.0_wp
             ir=1.0_wp/ir
             BC_face(2,1)%ir(i,j)=ir
             BC_face(2,1)%cosphi(i,j)=xr*ir
             BC_face(2,1)%sinphi(i,j)=yr*ir
          enddo
       enddo
       !if (is_2D) BC_face(2,1)%ir=0.5_wp*BC_face(2,1)%ir
    endif

    ! jmax (top) boundary condition
    ! =============================
    if (is_bc_TD(2,2)) then

       ! polar coordinates
       ! -----------------
       allocate(BC_face(2,2)%ir(ndx:nfx,ngh))
       allocate(BC_face(2,2)%cosphi(ndx:nfx,ngh))
       allocate(BC_face(2,2)%sinphi(ndx:nfx,ngh))

       do j=ny-ngh+1,ny
          l=j-ny+ngh
          do i=ndx,nfx
             if (is_curv) then
                xr=xc(i,j)-BC_face(2,2)%xrc(i)
                yr=yc(i,j)-BC_face(2,2)%yrc(i)
             else
                xr=x(i)-BC_face(2,2)%xrc(i)
                yr=y(j)-BC_face(2,2)%yrc(i)
             endif
             ir=sqrt(xr**2+yr**2)
             if (ir==0.0_wp)  ir=1.0_wp
             ir=1.0_wp/ir
             BC_face(2,2)%ir(i,l)=ir
             BC_face(2,2)%cosphi(i,l)=xr*ir
             BC_face(2,2)%sinphi(i,l)=yr*ir
          enddo
       enddo
       !if (is_2D) BC_face(2,2)%ir=0.5_wp*BC_face(2,2)%ir
    endif

    ! imin-jmin (left-bottom) boundary condition
    ! ==========================================
    if ((BC_edge(1,1,1)%sort==-1).or.(BC_edge(1,1,1)%sort==-2)) then

       ! polar coordinates
       ! -----------------
       allocate(BC_edge(1,1,1)%ir(ngh,ngh))
       allocate(BC_edge(1,1,1)%cosphi(ngh,ngh))
       allocate(BC_edge(1,1,1)%sinphi(ngh,ngh))

       do i=1,ngh
          do j=1,ngh
             if (is_curv) then
                xr=xc(i,j)-BC_face(1,1)%xrc(j)
                yr=yc(i,j)-BC_face(1,1)%yrc(j)
             else
                xr=x(i)-BC_face(1,1)%xrc(j)
                yr=y(j)-BC_face(1,1)%yrc(j)
             endif
             ir=sqrt(xr**2+yr**2)
             if (ir==0.0_wp)  ir=1.0_wp
             ir=1.0_wp/ir
             BC_edge(1,1,1)%ir(i,j)=ir
             BC_edge(1,1,1)%cosphi(i,j)=xr*ir
             BC_edge(1,1,1)%sinphi(i,j)=yr*ir
          enddo
       enddo
       !if (is_2D) BC_edge(1,1,1)%ir=0.5_wp*BC_edge(1,1,1)%ir
    endif

    ! imin-jmax (left-top) boundary condition
    ! =======================================
    if ((BC_edge(1,1,2)%sort==-1).or.(BC_edge(1,1,2)%sort==-2)) then

       ! polar coordinates
       ! -----------------
       allocate(BC_edge(1,1,2)%ir(ngh,ngh))
       allocate(BC_edge(1,1,2)%cosphi(ngh,ngh))
       allocate(BC_edge(1,1,2)%sinphi(ngh,ngh))

       do i=1,ngh
          do j=ny-ngh+1,ny
             l=j-ny+ngh
             if (is_curv) then
                xr=xc(i,j)-BC_face(1,1)%xrc(j)
                yr=yc(i,j)-BC_face(1,1)%yrc(j)
             else
                xr=x(i)-BC_face(1,1)%xrc(j)
                yr=y(j)-BC_face(1,1)%yrc(j)
             endif
             ir=sqrt(xr**2+yr**2)
             if (ir==0.0_wp)  ir=1.0_wp
             ir=1.0_wp/ir
             BC_edge(1,1,2)%ir(i,l)=ir
             BC_edge(1,1,2)%cosphi(i,l)=xr*ir
             BC_edge(1,1,2)%sinphi(i,l)=yr*ir
          enddo
       enddo
       !if (is_2D) BC_edge(1,1,2)%ir=0.5_wp*BC_edge(1,1,2)%ir
    endif

    ! imax-jmin (right-bottom) boundary condition
    ! ===========================================
    if ((BC_edge(1,2,1)%sort==-1).or.(BC_edge(1,2,1)%sort==-2)) then

       ! polar coordinates
       ! -----------------
       allocate(BC_edge(1,2,1)%ir(ngh,ngh))
       allocate(BC_edge(1,2,1)%cosphi(ngh,ngh))
       allocate(BC_edge(1,2,1)%sinphi(ngh,ngh))

       do i=nx-ngh+1,nx
          l=i-nx+ngh
          do j=1,ngh
             if (is_curv) then
                xr=xc(i,j)-BC_face(1,2)%xrc(j)
                yr=yc(i,j)-BC_face(1,2)%yrc(j)
             else
                xr=x(i)-BC_face(1,2)%xrc(j)
                yr=y(j)-BC_face(1,2)%yrc(j)
             endif
             ir=sqrt(xr**2+yr**2)
             if (ir==0.0_wp)  ir=1.0_wp
             ir=1.0_wp/ir
             BC_edge(1,2,1)%ir(l,j)=ir
             BC_edge(1,2,1)%cosphi(l,j)=xr*ir
             BC_edge(1,2,1)%sinphi(l,j)=yr*ir
          enddo
       enddo
       !if (is_2D) BC_edge(1,2,1)%ir=0.5_wp*BC_edge(1,2,1)%ir
    endif

    ! imax-jmax (right-top) boundary condition
    ! ========================================
    if ((BC_edge(1,2,2)%sort==-1).or.(BC_edge(1,2,2)%sort==-2)) then

       ! polar coordinates
       ! -----------------
       allocate(BC_edge(1,2,2)%ir(ngh,ngh))
       allocate(BC_edge(1,2,2)%cosphi(ngh,ngh))
       allocate(BC_edge(1,2,2)%sinphi(ngh,ngh))

       do i=nx-ngh+1,nx
          l=i-nx+ngh
          do j=ny-ngh+1,ny
             m=j-ny+ngh
             if (is_curv) then
                xr=xc(i,j)-BC_face(1,2)%xrc(j)
                yr=yc(i,j)-BC_face(1,2)%yrc(j)
             else
                xr=x(i)-BC_face(1,2)%xrc(j)
                yr=y(j)-BC_face(1,2)%yrc(j)
             endif
             ir=sqrt(xr**2+yr**2)
             if (ir==0.0_wp)  ir=1.0_wp
             ir=1.0_wp/ir
             BC_edge(1,2,2)%ir(l,m)=ir
             BC_edge(1,2,2)%cosphi(l,m)=xr*ir
             BC_edge(1,2,2)%sinphi(l,m)=yr*ir
          enddo
       enddo
       !if (is_2D) BC_edge(1,2,2)%ir=0.5_wp*BC_edge(1,2,2)%ir
    endif

  end subroutine init_param_TD2d

  !=============================================================================
  module subroutine init_param_TD2d_1pt
  !=============================================================================
    !> Preprocessing for Tam & Dong's boundary conditions
    !> allocations and compute geometrical parameters (Att. Radiation center def)
  !=============================================================================
    implicit none
    !---------------------------------------------------------------------------
    integer :: i,j,l,m
    real(wp) :: xr,yr,ir
    !---------------------------------------------------------------------------

    ! To be changed
    !!if ((.not.is_vortex).and.(type_vortex.eq.-1)) is_mean0=.false.
    is_mean0=.false.

    ! imin (left) boundary condition
    ! ==============================
    if (is_bc_TD(1,1)) then
       ! Index of left boundary
       ! ----------------------
       i=1

       ! polar coordinates
       ! -----------------
       allocate(BC_face(1,1)%ir(1,ndy_td1:nfy_td1))
       allocate(BC_face(1,1)%cosphi(1,ndy_td1:nfy_td1))
       allocate(BC_face(1,1)%sinphi(1,ndy_td1:nfy_td1))

       do j=ndy_td1,nfy_td1
          if (is_curv) then
             xr=xc(i,j)-BC_face(1,1)%xrc(j)
             yr=yc(i,j)-BC_face(1,1)%yrc(j)
          else
             xr=x(i)-BC_face(1,1)%xrc(j)
             yr=y(j)-BC_face(1,1)%yrc(j)
          endif
          ! if (SHIT) yr=0.0_wp
          ir=sqrt(xr**2+yr**2)
          if (ir==0.0_wp)  ir=1.0_wp
          ir=1.0_wp/ir
          BC_face(1,1)%ir(i,j)=ir
          BC_face(1,1)%cosphi(i,j)=xr*ir
          BC_face(1,1)%sinphi(i,j)=yr*ir
       enddo
    endif

    ! imax (right) boundary condition
    ! ===============================
    if (is_bc_TD(1,2)) then
       ! Index of right boundary
       ! -----------------------
       i=nx;l=1

       ! polar coordinates
       ! -----------------
       allocate(BC_face(1,2)%ir(1,ndy_td1:nfy_td1))
       allocate(BC_face(1,2)%cosphi(1,ndy_td1:nfy_td1))
       allocate(BC_face(1,2)%sinphi(1,ndy_td1:nfy_td1))

       do j=ndy_td1,nfy_td1
          if (is_curv) then
             xr=xc(i,j)-BC_face(1,2)%xrc(j)
             yr=yc(i,j)-BC_face(1,2)%yrc(j)
          else
             xr=x(i)-BC_face(1,2)%xrc(j)
             yr=y(j)-BC_face(1,2)%yrc(j)
          endif
          ! if (SHIT) yr=0.0_wp
          ir=sqrt(xr**2+yr**2)
          if (ir==0.0_wp)  ir=1.0_wp
          ir=1.0_wp/ir
          BC_face(1,2)%ir(l,j)=ir
          BC_face(1,2)%cosphi(l,j)=xr*ir
          BC_face(1,2)%sinphi(l,j)=yr*ir
       enddo
    endif

    ! jmin (bottom) boundary condition
    ! ================================
    if (is_bc_TD(2,1)) then
       ! Index of bottom boundary
       ! ------------------------
       j=1

       ! polar coordinates
       ! -----------------
       allocate(BC_face(2,1)%ir(ndx_td1:nfx_td1,1))
       allocate(BC_face(2,1)%cosphi(ndx_td1:nfx_td1,1))
       allocate(BC_face(2,1)%sinphi(ndx_td1:nfx_td1,1))

       do i=ndx_td1,nfx_td1
          if (is_curv) then
             xr=xc(i,j)-BC_face(2,1)%xrc(i)
             yr=yc(i,j)-BC_face(2,1)%yrc(i)
          else
             xr=x(i)-BC_face(2,1)%xrc(i)
             yr=y(j)-BC_face(2,1)%yrc(i)
          endif
          ir=sqrt(xr**2+yr**2)
          if (ir==0.0_wp)  ir=1.0_wp
          ir=1.0_wp/ir
          BC_face(2,1)%ir(i,j)=ir
          BC_face(2,1)%cosphi(i,j)=xr*ir
          BC_face(2,1)%sinphi(i,j)=yr*ir
       enddo
    endif

    ! jmax (top) boundary condition
    ! =============================
    if (is_bc_TD(2,2)) then
       ! Index of top boundary
       ! ---------------------
       j=ny;l=1

       ! polar coordinates
       ! -----------------
       allocate(BC_face(2,2)%ir(ndx_td1:nfx_td1,1))
       allocate(BC_face(2,2)%cosphi(ndx_td1:nfx_td1,1))
       allocate(BC_face(2,2)%sinphi(ndx_td1:nfx_td1,1))

       do i=ndx_td1,nfx_td1
          if (is_curv) then
             xr=xc(i,j)-BC_face(2,2)%xrc(i)
             yr=yc(i,j)-BC_face(2,2)%yrc(i)
          else
             xr=x(i)-BC_face(2,2)%xrc(i)
             yr=y(j)-BC_face(2,2)%yrc(i)
          endif
          ir=sqrt(xr**2+yr**2)
          if (ir==0.0_wp)  ir=1.0_wp
          ir=1.0_wp/ir
          BC_face(2,2)%ir(i,l)=ir
          BC_face(2,2)%cosphi(i,l)=xr*ir
          BC_face(2,2)%sinphi(i,l)=yr*ir
       enddo
    endif

     ! imin-jmin (left-bottom) boundary condition
    ! ==========================================
    if (BC_edge(1,1,1)%sort==-11) then
       ! polar coordinates
       ! -----------------
       allocate(BC_edge(1,1,1)%ir(2,2))
       allocate(BC_edge(1,1,1)%cosphi(2,2))
       allocate(BC_edge(1,1,1)%sinphi(2,2))

       do i=1,2
          do j=1,2
             if (BC_face(1,1)%sort==-11) then
                if (is_curv) then
                   xr=xc(i,j)-BC_face(1,1)%xrc(j)
                   yr=yc(i,j)-BC_face(1,1)%yrc(j)
                else
                   xr=x(i)-BC_face(1,1)%xrc(j)
                   yr=y(j)-BC_face(1,1)%yrc(j)
                endif
             else
                if (is_curv) then
                   xr=xc(i,j)-BC_face(2,1)%xrc(i)
                   yr=yc(i,j)-BC_face(2,1)%yrc(i)
                else
                   xr=x(i)-BC_face(2,1)%xrc(i)
                   yr=y(j)-BC_face(2,1)%yrc(i)
                endif
             endif
             ir=sqrt(xr**2+yr**2)
             if (ir==0.0_wp)  ir=1.0_wp
             ir=1.0_wp/ir
             BC_edge(1,1,1)%ir(i,j)=ir
             BC_edge(1,1,1)%cosphi(i,j)=xr*ir
             BC_edge(1,1,1)%sinphi(i,j)=yr*ir
          enddo
       enddo
    endif

    ! imin-jmax (left-top) boundary condition
    ! =======================================
    if (BC_edge(1,1,2)%sort==-11) then
       ! polar coordinates
       ! -----------------
       allocate(BC_edge(1,1,2)%ir(2,2))
       allocate(BC_edge(1,1,2)%cosphi(2,2))
       allocate(BC_edge(1,1,2)%sinphi(2,2))

       do i=1,2
          do j=ny-1,ny
             l=j-ny+2
             if (BC_face(1,1)%sort==-11) then
                if (is_curv) then
                   xr=xc(i,j)-BC_face(1,1)%xrc(j)
                   yr=yc(i,j)-BC_face(1,1)%yrc(j)
                else
                   xr=x(i)-BC_face(1,1)%xrc(j)
                   yr=y(j)-BC_face(1,1)%yrc(j)
                endif
             else
                if (is_curv) then
                   xr=xc(i,j)-BC_face(2,2)%xrc(i)
                   yr=yc(i,j)-BC_face(2,2)%yrc(i)
                else
                   xr=x(i)-BC_face(2,2)%xrc(i)
                   yr=y(j)-BC_face(2,2)%yrc(i)
                endif
             endif
             ir=sqrt(xr**2+yr**2)
             if (ir==0.0_wp)  ir=1.0_wp
             ir=1.0_wp/ir
             BC_edge(1,1,2)%ir(i,l)=ir
             BC_edge(1,1,2)%cosphi(i,l)=xr*ir
             BC_edge(1,1,2)%sinphi(i,l)=yr*ir
          enddo
       enddo
    endif

    ! imax-jmin (right-bottom) boundary condition
    ! ===========================================
    if (BC_edge(1,2,1)%sort==-11) then
       ! polar coordinates
       ! -----------------
       allocate(BC_edge(1,2,1)%ir(2,2))
       allocate(BC_edge(1,2,1)%cosphi(2,2))
       allocate(BC_edge(1,2,1)%sinphi(2,2))

       do i=nx-1,nx
          l=i-nx+2
          do j=1,2
             if (BC_face(1,2)%sort==-11) then
                if (is_curv) then
                   xr=xc(i,j)-BC_face(1,2)%xrc(j)
                   yr=yc(i,j)-BC_face(1,2)%yrc(j)
                else
                   xr=x(i)-BC_face(1,2)%xrc(j)
                   yr=y(j)-BC_face(1,2)%yrc(j)
                endif
             else
                ! call mpistop('',0)
                if (is_curv) then
                   xr=xc(i,j)-BC_face(2,1)%xrc(i)
                   yr=yc(i,j)-BC_face(2,1)%yrc(i)
                else
                   xr=x(i)-BC_face(2,1)%xrc(i)
                   yr=y(j)-BC_face(2,1)%yrc(i)
                endif
             endif
             ir=sqrt(xr**2+yr**2)
             if (ir==0.0_wp)  ir=1.0_wp
             ir=1.0_wp/ir
             BC_edge(1,2,1)%ir(l,j)=ir
             BC_edge(1,2,1)%cosphi(l,j)=xr*ir
             BC_edge(1,2,1)%sinphi(l,j)=yr*ir
          enddo
       enddo
    endif

    ! imax-jmax (right-top) boundary condition
    ! ========================================
    if (BC_edge(1,2,2)%sort==-11) then
       i=nx;j=ny;l=1;m=1

       ! polar coordinates
       ! -----------------
       allocate(BC_edge(1,2,2)%ir(2,2))
       allocate(BC_edge(1,2,2)%cosphi(2,2))
       allocate(BC_edge(1,2,2)%sinphi(2,2))

       do i=nx-1,nx
          l=i-nx+2
          do j=ny-1,ny
             m=j-ny+2
             if (BC_face(1,2)%sort==-11) then
                if (is_curv) then
                   xr=xc(i,j)-BC_face(1,2)%xrc(j)
                   yr=yc(i,j)-BC_face(1,2)%yrc(j)
                else
                   xr=x(i)-BC_face(1,2)%xrc(j)
                   yr=y(j)-BC_face(1,2)%yrc(j)
                endif
             else
                if (is_curv) then
                   xr=xc(i,j)-BC_face(2,2)%xrc(i)
                   yr=yc(i,j)-BC_face(2,2)%yrc(i)
                else
                   xr=x(i)-BC_face(2,2)%xrc(i)
                   yr=y(j)-BC_face(2,2)%yrc(i)
                endif
             endif
             ir=sqrt(xr**2+yr**2)
             if (ir==0.0_wp)  ir=1.0_wp
             ir=1.0_wp/ir
             BC_edge(1,2,2)%ir(l,m)=ir
             BC_edge(1,2,2)%cosphi(l,m)=xr*ir
             BC_edge(1,2,2)%sinphi(l,m)=yr*ir
          enddo
       enddo
    endif

  end subroutine init_param_TD2d_1pt

  !===============================================================================
  module subroutine free_bc_TD2d
  !===============================================================================
    !> Closing of Tam & Dong's BC: deallocations of memory
  !===============================================================================
    implicit none
    !-----------------------------------------------------------------------------
    integer :: l,m
    !-----------------------------------------------------------------------------

    ! Free BC arrays for faces
    ! ========================
    do l=1,2
       do m=1,2
          if (is_bc_TD(l,m)) then
             deallocate(BC_face(l,m)%ir,BC_face(l,m)%U0)
             deallocate(BC_face(l,m)%cosphi,BC_face(l,m)%sinphi)
          endif
       enddo
    enddo

    ! Free BC arrays for edges
    ! ========================
    do l=1,2
       do m=1,2
          if (BC_edge(1,l,m)%sort==-1) then
             deallocate(BC_edge(1,l,m)%ir,BC_edge(1,l,m)%cosphi,BC_edge(1,l,m)%sinphi)
          endif
       enddo
    enddo

  end subroutine free_bc_TD2d
   
end submodule smod_init_TamDong2
