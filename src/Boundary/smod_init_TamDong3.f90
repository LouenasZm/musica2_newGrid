!===============================================================================
submodule (mod_init_TamDong) smod_init_TamDong3
!===============================================================================
  !> author: XG
  !> date: January 2024
  !> Initialization of Tam and Dong's boundary conditions
  !> routines: init_param_TD3d & free_bc_TD3d
!===============================================================================
  implicit none

contains

  !=============================================================================
  module subroutine init_param_TD3d
  !=============================================================================
    !> Preprocessing for Tam & Dong's boundary conditions
    !> allocations and compute geometrical parameters (Att. Radiation center def)
  !=============================================================================
    implicit none
    !---------------------------------------------------------------------------
    integer :: i,j,k,l,m,n
    real(wp) :: xr,yr,zr,ir,rp
    !---------------------------------------------------------------------------

    ! To be changed
    !!if ((.not.is_vortex).and.(type_vortex.eq.-1)) is_mean0=.false.
    is_mean0=.false.

    ! imin (left) boundary condition
    ! ==============================
    if (is_bc_TD(1,1)) then
       
       ! spherical coordinates
       ! ---------------------
       allocate(BC_face(1,1)%i_r(ngh,ndy:nfy,ndz:nfz))
       allocate(BC_face(1,1)%cosp(ngh,ndy:nfy,ndz:nfz))
       allocate(BC_face(1,1)%sinp(ngh,ndy:nfy,ndz:nfz))
       allocate(BC_face(1,1)%cost(ngh,ndy:nfy,ndz:nfz))
       allocate(BC_face(1,1)%sint(ngh,ndy:nfy,ndz:nfz))
       allocate(BC_face(1,1)%costcosp(ngh,ndy:nfy,ndz:nfz))
       allocate(BC_face(1,1)%costsinp(ngh,ndy:nfy,ndz:nfz))
       allocate(BC_face(1,1)%sintcosp(ngh,ndy:nfy,ndz:nfz))
       allocate(BC_face(1,1)%sintsinp(ngh,ndy:nfy,ndz:nfz))
       
       do i=1,ngh
          do j=ndy,nfy
             do k=ndz,nfz
                if (is_curv3) then
                   xr=xc3(i,j,k)-BC_face(1,1)%xrc(j)
                   yr=yc3(i,j,k)-BC_face(1,1)%yrc(j)
                   zr=zc3(i,j,k)-BC_face(1,1)%zrc(j)
                else
                   if (is_curv) then
                      xr=xc(i,j)-BC_face(1,1)%xrc(j)
                      yr=yc(i,j)-BC_face(1,1)%yrc(j)
                   else
                      xr=x(i)-BC_face(1,1)%xrc(j)
                      yr=y(j)-BC_face(1,1)%yrc(j)
                   endif
                   zr=z(k)-BC_face(1,1)%zrc(j)
                endif
                ir=sqrt(xr**2+yr**2+zr**2)
                rp=sqrt(xr**2+yr**2)
                if (ir==0.0_wp)  ir=1.0_wp
                if (rp==0.0_wp)  rp=1.0_wp
                ir=1.0_wp/ir
                BC_face(1,1)%i_r(i,j,k)=ir
                BC_face(1,1)%cosp(i,j,k)=xr/rp
                BC_face(1,1)%sinp(i,j,k)=yr/rp
                BC_face(1,1)%cost(i,j,k)=zr*ir
                BC_face(1,1)%sint(i,j,k)=rp*ir
             enddo
          enddo
       enddo
       BC_face(1,1)%sintsinp=BC_face(1,1)%sint*BC_face(1,1)%sinp
       BC_face(1,1)%sintcosp=BC_face(1,1)%sint*BC_face(1,1)%cosp
       BC_face(1,1)%costsinp=BC_face(1,1)%cost*BC_face(1,1)%sinp
       BC_face(1,1)%costcosp=BC_face(1,1)%cost*BC_face(1,1)%cosp
    endif

    ! imax (right) boundary condition
    ! ===============================
    if (is_bc_TD(1,2)) then

       ! spherical coordinates
       ! ---------------------
       allocate(BC_face(1,2)%i_r(ngh,ndy:nfy,ndz:nfz))
       allocate(BC_face(1,2)%cosp(ngh,ndy:nfy,ndz:nfz))
       allocate(BC_face(1,2)%sinp(ngh,ndy:nfy,ndz:nfz))
       allocate(BC_face(1,2)%cost(ngh,ndy:nfy,ndz:nfz))
       allocate(BC_face(1,2)%sint(ngh,ndy:nfy,ndz:nfz))
       allocate(BC_face(1,2)%costcosp(ngh,ndy:nfy,ndz:nfz))
       allocate(BC_face(1,2)%costsinp(ngh,ndy:nfy,ndz:nfz))
       allocate(BC_face(1,2)%sintcosp(ngh,ndy:nfy,ndz:nfz))
       allocate(BC_face(1,2)%sintsinp(ngh,ndy:nfy,ndz:nfz))
       
       do i=nx-ngh+1,nx
          l=i-nx+ngh
          do j=ndy,nfy
             do k=ndz,nfz
                if (is_curv3) then
                   xr=xc3(i,j,k)-BC_face(1,2)%xrc(j)
                   yr=yc3(i,j,k)-BC_face(1,2)%yrc(j)
                   zr=zc3(i,j,k)-BC_face(1,2)%zrc(j)
                else
                   if (is_curv) then
                      xr=xc(i,j)-BC_face(1,2)%xrc(j)
                      yr=yc(i,j)-BC_face(1,2)%yrc(j)
                   else
                      xr=x(i)-BC_face(1,2)%xrc(j)
                      yr=y(j)-BC_face(1,2)%yrc(j)
                   endif
                   zr=z(k)-BC_face(1,2)%zrc(j)
                endif
                ir=sqrt(xr**2+yr**2+zr**2)
                rp=sqrt(xr**2+yr**2)
                if (ir==0.0_wp)  ir=1.0_wp
                if (rp==0.0_wp)  rp=1.0_wp
                ir=1.0_wp/ir
                BC_face(1,2)%i_r(l,j,k)=ir
                BC_face(1,2)%cosp(l,j,k)=xr/rp
                BC_face(1,2)%sinp(l,j,k)=yr/rp
                BC_face(1,2)%cost(l,j,k)=zr*ir
                BC_face(1,2)%sint(l,j,k)=rp*ir
             enddo
          enddo
       enddo
       BC_face(1,2)%sintsinp=BC_face(1,2)%sint*BC_face(1,2)%sinp
       BC_face(1,2)%sintcosp=BC_face(1,2)%sint*BC_face(1,2)%cosp
       BC_face(1,2)%costsinp=BC_face(1,2)%cost*BC_face(1,2)%sinp
       BC_face(1,2)%costcosp=BC_face(1,2)%cost*BC_face(1,2)%cosp
    endif

    ! jmin (bottom) boundary condition
    ! ================================
    if (is_bc_TD(2,1)) then

       ! spherical coordinates
       ! ---------------------
       allocate(BC_face(2,1)%i_r(ndx:nfx,ngh,ndz:nfz))
       allocate(BC_face(2,1)%cosp(ndx:nfx,ngh,ndz:nfz))
       allocate(BC_face(2,1)%sinp(ndx:nfx,ngh,ndz:nfz))
       allocate(BC_face(2,1)%cost(ndx:nfx,ngh,ndz:nfz))
       allocate(BC_face(2,1)%sint(ndx:nfx,ngh,ndz:nfz))
       allocate(BC_face(2,1)%costcosp(ndx:nfx,ngh,ndz:nfz))
       allocate(BC_face(2,1)%costsinp(ndx:nfx,ngh,ndz:nfz))
       allocate(BC_face(2,1)%sintcosp(ndx:nfx,ngh,ndz:nfz))
       allocate(BC_face(2,1)%sintsinp(ndx:nfx,ngh,ndz:nfz))
       
       do j=1,ngh
          do i=ndx,nfx
             do k=ndz,nfz
                if (is_curv3) then
                   xr=xc3(i,j,k)-BC_face(2,1)%xrc(i)
                   yr=yc3(i,j,k)-BC_face(2,1)%yrc(i)
                   zr=zc3(i,j,k)-BC_face(2,1)%zrc(i)
                else
                   if (is_curv) then
                      xr=xc(i,j)-BC_face(2,1)%xrc(i)
                      yr=yc(i,j)-BC_face(2,1)%yrc(i)
                   else
                      xr=x(i)-BC_face(2,1)%xrc(i)
                      yr=y(j)-BC_face(2,1)%yrc(i)
                   endif
                   zr=z(k)-BC_face(2,1)%zrc(i)
                endif
                ir=sqrt(xr**2+yr**2+zr**2)
                rp=sqrt(xr**2+yr**2)
                if (ir==0.0_wp)  ir=1.0_wp
                if (rp==0.0_wp)  rp=1.0_wp
                ir=1.0_wp/ir
                BC_face(2,1)%i_r(i,j,k)=ir
                BC_face(2,1)%cosp(i,j,k)=xr/rp
                BC_face(2,1)%sinp(i,j,k)=yr/rp
                BC_face(2,1)%cost(i,j,k)=zr*ir
                BC_face(2,1)%sint(i,j,k)=rp*ir
             enddo
          enddo
       enddo
       BC_face(2,1)%sintsinp=BC_face(2,1)%sint*BC_face(2,1)%sinp
       BC_face(2,1)%sintcosp=BC_face(2,1)%sint*BC_face(2,1)%cosp
       BC_face(2,1)%costsinp=BC_face(2,1)%cost*BC_face(2,1)%sinp
       BC_face(2,1)%costcosp=BC_face(2,1)%cost*BC_face(2,1)%cosp
    endif

    ! jmax (top) boundary condition
    ! =============================
    if (is_bc_TD(2,2)) then

       ! spherical coordinates
       ! ---------------------
       allocate(BC_face(2,2)%i_r(ndx:nfx,ngh,ndz:nfz))
       allocate(BC_face(2,2)%cosp(ndx:nfx,ngh,ndz:nfz))
       allocate(BC_face(2,2)%sinp(ndx:nfx,ngh,ndz:nfz))
       allocate(BC_face(2,2)%cost(ndx:nfx,ngh,ndz:nfz))
       allocate(BC_face(2,2)%sint(ndx:nfx,ngh,ndz:nfz))
       allocate(BC_face(2,2)%costcosp(ndx:nfx,ngh,ndz:nfz))
       allocate(BC_face(2,2)%costsinp(ndx:nfx,ngh,ndz:nfz))
       allocate(BC_face(2,2)%sintcosp(ndx:nfx,ngh,ndz:nfz))
       allocate(BC_face(2,2)%sintsinp(ndx:nfx,ngh,ndz:nfz))
       
       do j=ny-ngh+1,ny
          l=j-ny+ngh
          do i=ndx,nfx
             do k=ndz,nfz
                if (is_curv3) then
                   xr=xc3(i,j,k)-BC_face(2,2)%xrc(i)
                   yr=yc3(i,j,k)-BC_face(2,2)%yrc(i)
                   zr=zc3(i,j,k)-BC_face(2,2)%zrc(i)
                else
                   if (is_curv) then
                      xr=xc(i,j)-BC_face(2,2)%xrc(i)
                      yr=yc(i,j)-BC_face(2,2)%yrc(i)
                   else
                      xr=x(i)-BC_face(2,2)%xrc(i)
                      yr=y(j)-BC_face(2,2)%yrc(i)
                   endif
                   zr=z(k)-BC_face(2,2)%zrc(i)
                endif
                ir=sqrt(xr**2+yr**2+zr**2)
                rp=sqrt(xr**2+yr**2)
                if (ir==0.0_wp)  ir=1.0_wp
                if (rp==0.0_wp)  rp=1.0_wp
                ir=1.0_wp/ir
                BC_face(2,2)%i_r(i,l,k)=ir
                BC_face(2,2)%cosp(i,l,k)=xr/rp
                BC_face(2,2)%sinp(i,l,k)=yr/rp
                BC_face(2,2)%cost(i,l,k)=zr*ir
                BC_face(2,2)%sint(i,l,k)=rp*ir
             enddo
          enddo
       enddo
       BC_face(2,2)%sintsinp=BC_face(2,2)%sint*BC_face(2,2)%sinp
       BC_face(2,2)%sintcosp=BC_face(2,2)%sint*BC_face(2,2)%cosp
       BC_face(2,2)%costsinp=BC_face(2,2)%cost*BC_face(2,2)%sinp
       BC_face(2,2)%costcosp=BC_face(2,2)%cost*BC_face(2,2)%cosp
    endif
    
    ! kmin (front) boundary condition
    ! ===============================
    if (is_bc_TD(3,1)) then

       ! spherical coordinates
       ! ---------------------
       allocate(BC_face(3,1)%i_r(ndx:nfx,ndy:nfy,ngh))
       allocate(BC_face(3,1)%cosp(ndx:nfx,ndy:nfy,ngh))
       allocate(BC_face(3,1)%sinp(ndx:nfx,ndy:nfy,ngh))
       allocate(BC_face(3,1)%cost(ndx:nfx,ndy:nfy,ngh))
       allocate(BC_face(3,1)%sint(ndx:nfx,ndy:nfy,ngh))
       allocate(BC_face(3,1)%costcosp(ndx:nfx,ndy:nfy,ngh))
       allocate(BC_face(3,1)%costsinp(ndx:nfx,ndy:nfy,ngh))
       allocate(BC_face(3,1)%sintcosp(ndx:nfx,ndy:nfy,ngh))
       allocate(BC_face(3,1)%sintsinp(ndx:nfx,ndy:nfy,ngh))
       
       do k=1,ngh
          do i=ndx,nfx
             do j=ndy,nfy
                if (is_curv3) then
                   xr=xc3(i,j,k)-BC_face(3,1)%xrc(k)
                   yr=yc3(i,j,k)-BC_face(3,1)%yrc(k)
                   zr=zc3(i,j,k)-BC_face(3,1)%zrc(k)
                else
                   if (is_curv) then
                      xr=xc(i,j)-BC_face(3,1)%xrc(k)
                      yr=yc(i,j)-BC_face(3,1)%yrc(k)
                   else
                      xr=x(i)-BC_face(3,1)%xrc(k)
                      yr=y(j)-BC_face(3,1)%yrc(k)
                   endif
                   zr=z(k)-BC_face(3,1)%zrc(k)
                endif
                ir=sqrt(xr**2+yr**2+zr**2)
                rp=sqrt(xr**2+yr**2)
                if (ir==0.0_wp)  ir=1.0_wp
                if (rp==0.0_wp)  rp=1.0_wp
                ir=1.0_wp/ir
                BC_face(3,1)%i_r(i,j,k)=ir
                BC_face(3,1)%cosp(i,j,k)=xr/rp
                BC_face(3,1)%sinp(i,j,k)=yr/rp
                BC_face(3,1)%cost(i,j,k)=zr*ir
                BC_face(3,1)%sint(i,j,k)=rp*ir
             enddo
          enddo
       enddo
       BC_face(3,1)%sintsinp=BC_face(3,1)%sint*BC_face(3,1)%sinp
       BC_face(3,1)%sintcosp=BC_face(3,1)%sint*BC_face(3,1)%cosp
       BC_face(3,1)%costsinp=BC_face(3,1)%cost*BC_face(3,1)%sinp
       BC_face(3,1)%costcosp=BC_face(3,1)%cost*BC_face(3,1)%cosp
    endif

    ! kmax (back) boundary condition
    ! ==============================
    if (is_bc_TD(3,2)) then

       ! spherical coordinates
       ! ---------------------
       allocate(BC_face(3,2)%i_r(ndx:nfx,ndy:nfy,ngh))
       allocate(BC_face(3,2)%cosp(ndx:nfx,ndy:nfy,ngh))
       allocate(BC_face(3,2)%sinp(ndx:nfx,ndy:nfy,ngh))
       allocate(BC_face(3,2)%cost(ndx:nfx,ndy:nfy,ngh))
       allocate(BC_face(3,2)%sint(ndx:nfx,ndy:nfy,ngh))
       allocate(BC_face(3,2)%costcosp(ndx:nfx,ndy:nfy,ngh))
       allocate(BC_face(3,2)%costsinp(ndx:nfx,ndy:nfy,ngh))
       allocate(BC_face(3,2)%sintcosp(ndx:nfx,ndy:nfy,ngh))
       allocate(BC_face(3,2)%sintsinp(ndx:nfx,ndy:nfy,ngh))
       
       do k=nz-ngh+1,nz
          l=k-nz+ngh
          do i=ndx,nfx
             do j=ndy,nfy
                if (is_curv3) then
                   xr=xc3(i,j,k)-BC_face(3,2)%xrc(k)
                   yr=yc3(i,j,k)-BC_face(3,2)%yrc(k)
                   zr=zc3(i,j,k)-BC_face(3,2)%zrc(k)
                else
                   if (is_curv) then
                      xr=xc(i,j)-BC_face(3,2)%xrc(k)
                      yr=yc(i,j)-BC_face(3,2)%yrc(k)
                   else
                      xr=x(i)-BC_face(3,2)%xrc(k)
                      yr=y(j)-BC_face(3,2)%yrc(k)
                   endif
                   zr=z(k)-BC_face(3,2)%zrc(k)
                endif
                ir=sqrt(xr**2+yr**2+zr**2)
                rp=sqrt(xr**2+yr**2)
                if (ir==0.0_wp)  ir=1.0_wp
                if (rp==0.0_wp)  rp=1.0_wp
                ir=1.0_wp/ir
                BC_face(3,2)%i_r(i,j,l)=ir
                BC_face(3,2)%cosp(i,j,l)=xr/rp
                BC_face(3,2)%sinp(i,j,l)=yr/rp
                BC_face(3,2)%cost(i,j,l)=zr*ir
                BC_face(3,2)%sint(i,j,l)=rp*ir
             enddo
          enddo
       enddo
       BC_face(3,2)%sintsinp=BC_face(3,2)%sint*BC_face(3,2)%sinp
       BC_face(3,2)%sintcosp=BC_face(3,2)%sint*BC_face(3,2)%cosp
       BC_face(3,2)%costsinp=BC_face(3,2)%cost*BC_face(3,2)%sinp
       BC_face(3,2)%costcosp=BC_face(3,2)%cost*BC_face(3,2)%cosp
    endif
    
    ! imin-jmin (left-bottom) edge boundary condition
    ! ===============================================
    if ((BC_edge(1,1,1)%sort==-1).or.(BC_edge(1,1,1)%sort==-2)) then

       ! spherical coordinates
       ! ---------------------
       allocate(BC_edge(1,1,1)%i_r(ngh,ngh,ndz:nfz))
       allocate(BC_edge(1,1,1)%cosp(ngh,ngh,ndz:nfz))
       allocate(BC_edge(1,1,1)%sinp(ngh,ngh,ndz:nfz))
       allocate(BC_edge(1,1,1)%cost(ngh,ngh,ndz:nfz))
       allocate(BC_edge(1,1,1)%sint(ngh,ngh,ndz:nfz))
       allocate(BC_edge(1,1,1)%costcosp(ngh,ngh,ndz:nfz))
       allocate(BC_edge(1,1,1)%costsinp(ngh,ngh,ndz:nfz))
       allocate(BC_edge(1,1,1)%sintcosp(ngh,ngh,ndz:nfz))
       allocate(BC_edge(1,1,1)%sintsinp(ngh,ngh,ndz:nfz))
       
       do i=1,ngh
          do j=1,ngh
             do k=ndz,nfz
                if (is_curv3) then
                   xr=xc3(i,j,k)-BC_face(1,1)%xrc(j)
                   yr=yc3(i,j,k)-BC_face(1,1)%yrc(j)
                   zr=zc3(i,j,k)-BC_face(1,1)%zrc(j)
                else
                   if (is_curv) then
                      xr=xc(i,j)-BC_face(1,1)%xrc(j)
                      yr=yc(i,j)-BC_face(1,1)%yrc(j)
                   else
                      xr=x(i)-BC_face(1,1)%xrc(j)
                      yr=y(j)-BC_face(1,1)%yrc(j)
                   endif
                   zr=z(k)-BC_face(1,1)%zrc(j)
                endif
                ir=sqrt(xr**2+yr**2+zr**2)
                rp=sqrt(xr**2+yr**2)
                if (ir==0.0_wp)  ir=1.0_wp
                if (rp==0.0_wp)  rp=1.0_wp
                ir=1.0_wp/ir
                BC_edge(1,1,1)%i_r(i,j,k)=ir
                BC_edge(1,1,1)%cosp(i,j,k)=xr/rp
                BC_edge(1,1,1)%sinp(i,j,k)=yr/rp
                BC_edge(1,1,1)%cost(i,j,k)=zr*ir
                BC_edge(1,1,1)%sint(i,j,k)=rp*ir
             enddo
          enddo
       enddo
       BC_edge(1,1,1)%sintsinp=BC_edge(1,1,1)%sint*BC_edge(1,1,1)%sinp
       BC_edge(1,1,1)%sintcosp=BC_edge(1,1,1)%sint*BC_edge(1,1,1)%cosp
       BC_edge(1,1,1)%costsinp=BC_edge(1,1,1)%cost*BC_edge(1,1,1)%sinp
       BC_edge(1,1,1)%costcosp=BC_edge(1,1,1)%cost*BC_edge(1,1,1)%cosp
    endif

    ! imin-jmax (left-top) edge boundary condition
    ! ============================================
    if ((BC_edge(1,1,2)%sort==-1).or.(BC_edge(1,1,2)%sort==-2)) then

       ! spherical coordinates
       ! ---------------------
       allocate(BC_edge(1,1,2)%i_r(ngh,ngh,ndz:nfz))
       allocate(BC_edge(1,1,2)%cosp(ngh,ngh,ndz:nfz))
       allocate(BC_edge(1,1,2)%sinp(ngh,ngh,ndz:nfz))
       allocate(BC_edge(1,1,2)%cost(ngh,ngh,ndz:nfz))
       allocate(BC_edge(1,1,2)%sint(ngh,ngh,ndz:nfz))
       allocate(BC_edge(1,1,2)%costcosp(ngh,ngh,ndz:nfz))
       allocate(BC_edge(1,1,2)%costsinp(ngh,ngh,ndz:nfz))
       allocate(BC_edge(1,1,2)%sintcosp(ngh,ngh,ndz:nfz))
       allocate(BC_edge(1,1,2)%sintsinp(ngh,ngh,ndz:nfz))
       
       do i=1,ngh
          do j=ny-ngh+1,ny
             l=j-ny+ngh
             do k=ndz,nfz
                if (is_curv3) then
                   xr=xc3(i,j,k)-BC_face(1,1)%xrc(j)
                   yr=yc3(i,j,k)-BC_face(1,1)%yrc(j)
                   zr=zc3(i,j,k)-BC_face(1,1)%zrc(j)
                else
                   if (is_curv) then
                      xr=xc(i,j)-BC_face(1,1)%xrc(j)
                      yr=yc(i,j)-BC_face(1,1)%yrc(j)
                   else
                      xr=x(i)-BC_face(1,1)%xrc(j)
                      yr=y(j)-BC_face(1,1)%yrc(j)
                   endif
                   zr=z(k)-BC_face(1,1)%zrc(j)
                endif
                ir=sqrt(xr**2+yr**2+zr**2)
                rp=sqrt(xr**2+yr**2)
                if (ir==0.0_wp)  ir=1.0_wp
                if (rp==0.0_wp)  rp=1.0_wp
                ir=1.0_wp/ir
                BC_edge(1,1,2)%i_r(i,l,k)=ir
                BC_edge(1,1,2)%cosp(i,l,k)=xr/rp
                BC_edge(1,1,2)%sinp(i,l,k)=yr/rp
                BC_edge(1,1,2)%cost(i,l,k)=zr*ir
                BC_edge(1,1,2)%sint(i,l,k)=rp*ir
             enddo
          enddo
       enddo
       BC_edge(1,1,2)%sintsinp=BC_edge(1,1,2)%sint*BC_edge(1,1,2)%sinp
       BC_edge(1,1,2)%sintcosp=BC_edge(1,1,2)%sint*BC_edge(1,1,2)%cosp
       BC_edge(1,1,2)%costsinp=BC_edge(1,1,2)%cost*BC_edge(1,1,2)%sinp
       BC_edge(1,1,2)%costcosp=BC_edge(1,1,2)%cost*BC_edge(1,1,2)%cosp
    endif

    ! imax-jmin (right-bottom) edge boundary condition
    ! ================================================
    if ((BC_edge(1,2,1)%sort==-1).or.(BC_edge(1,2,1)%sort==-2)) then

       ! spherical coordinates
       ! ---------------------
       allocate(BC_edge(1,2,1)%i_r(ngh,ngh,ndz:nfz))
       allocate(BC_edge(1,2,1)%cosp(ngh,ngh,ndz:nfz))
       allocate(BC_edge(1,2,1)%sinp(ngh,ngh,ndz:nfz))
       allocate(BC_edge(1,2,1)%cost(ngh,ngh,ndz:nfz))
       allocate(BC_edge(1,2,1)%sint(ngh,ngh,ndz:nfz))
       allocate(BC_edge(1,2,1)%costcosp(ngh,ngh,ndz:nfz))
       allocate(BC_edge(1,2,1)%costsinp(ngh,ngh,ndz:nfz))
       allocate(BC_edge(1,2,1)%sintcosp(ngh,ngh,ndz:nfz))
       allocate(BC_edge(1,2,1)%sintsinp(ngh,ngh,ndz:nfz))
       
       do i=nx-ngh+1,nx
          l=i-nx+ngh
          do j=1,ngh
             do k=ndz,nfz
                if (is_curv3) then
                   xr=xc3(i,j,k)-BC_face(1,2)%xrc(j)
                   yr=yc3(i,j,k)-BC_face(1,2)%yrc(j)
                   zr=zc3(i,j,k)-BC_face(1,2)%zrc(j)
                else
                   if (is_curv) then
                      xr=xc(i,j)-BC_face(1,2)%xrc(j)
                      yr=yc(i,j)-BC_face(1,2)%yrc(j)
                   else
                      xr=x(i)-BC_face(1,2)%xrc(j)
                      yr=y(j)-BC_face(1,2)%yrc(j)
                   endif
                   zr=z(k)-BC_face(1,2)%zrc(j)
                endif
                ir=sqrt(xr**2+yr**2+zr**2)
                rp=sqrt(xr**2+yr**2)
                if (ir==0.0_wp)  ir=1.0_wp
                if (rp==0.0_wp)  rp=1.0_wp
                ir=1.0_wp/ir
                BC_edge(1,2,1)%i_r(l,j,k)=ir
                BC_edge(1,2,1)%cosp(l,j,k)=xr/rp
                BC_edge(1,2,1)%sinp(l,j,k)=yr/rp
                BC_edge(1,2,1)%cost(l,j,k)=zr*ir
                BC_edge(1,2,1)%sint(l,j,k)=rp*ir
             enddo
          enddo
       enddo
       BC_edge(1,2,1)%sintsinp=BC_edge(1,2,1)%sint*BC_edge(1,2,1)%sinp
       BC_edge(1,2,1)%sintcosp=BC_edge(1,2,1)%sint*BC_edge(1,2,1)%cosp
       BC_edge(1,2,1)%costsinp=BC_edge(1,2,1)%cost*BC_edge(1,2,1)%sinp
       BC_edge(1,2,1)%costcosp=BC_edge(1,2,1)%cost*BC_edge(1,2,1)%cosp
    endif

    ! imax-jmax (right-top) edge boundary condition
    ! =============================================
    if ((BC_edge(1,2,2)%sort==-1).or.(BC_edge(1,2,2)%sort==-2)) then

       ! spherical coordinates
       ! ---------------------
       allocate(BC_edge(1,2,2)%i_r(ngh,ngh,ndz:nfz))
       allocate(BC_edge(1,2,2)%cosp(ngh,ngh,ndz:nfz))
       allocate(BC_edge(1,2,2)%sinp(ngh,ngh,ndz:nfz))
       allocate(BC_edge(1,2,2)%cost(ngh,ngh,ndz:nfz))
       allocate(BC_edge(1,2,2)%sint(ngh,ngh,ndz:nfz))
       allocate(BC_edge(1,2,2)%costcosp(ngh,ngh,ndz:nfz))
       allocate(BC_edge(1,2,2)%costsinp(ngh,ngh,ndz:nfz))
       allocate(BC_edge(1,2,2)%sintcosp(ngh,ngh,ndz:nfz))
       allocate(BC_edge(1,2,2)%sintsinp(ngh,ngh,ndz:nfz))
       
       do i=nx-ngh+1,nx
          l=i-nx+ngh
          do j=ny-ngh+1,ny
             m=j-ny+ngh
             do k=ndz,nfz
                if (is_curv3) then
                   xr=xc3(i,j,k)-BC_face(1,2)%xrc(j)
                   yr=yc3(i,j,k)-BC_face(1,2)%yrc(j)
                   zr=zc3(i,j,k)-BC_face(1,2)%zrc(j)
                else
                   if (is_curv) then
                      xr=xc(i,j)-BC_face(1,2)%xrc(j)
                      yr=yc(i,j)-BC_face(1,2)%yrc(j)
                   else
                      xr=x(i)-BC_face(1,2)%xrc(j)
                      yr=y(j)-BC_face(1,2)%yrc(j)
                   endif
                   zr=z(k)-BC_face(1,2)%zrc(j)
                endif
                ir=sqrt(xr**2+yr**2+zr**2)
                rp=sqrt(xr**2+yr**2)
                if (ir==0.0_wp)  ir=1.0_wp
                if (rp==0.0_wp)  rp=1.0_wp
                ir=1.0_wp/ir
                BC_edge(1,2,2)%i_r(l,m,k)=ir
                BC_edge(1,2,2)%cosp(l,m,k)=xr/rp
                BC_edge(1,2,2)%sinp(l,m,k)=yr/rp
                BC_edge(1,2,2)%cost(l,m,k)=zr*ir
                BC_edge(1,2,2)%sint(l,m,k)=rp*ir
             enddo
          enddo
       enddo
       BC_edge(1,2,2)%sintsinp=BC_edge(1,2,2)%sint*BC_edge(1,2,2)%sinp
       BC_edge(1,2,2)%sintcosp=BC_edge(1,2,2)%sint*BC_edge(1,2,2)%cosp
       BC_edge(1,2,2)%costsinp=BC_edge(1,2,2)%cost*BC_edge(1,2,2)%sinp
       BC_edge(1,2,2)%costcosp=BC_edge(1,2,2)%cost*BC_edge(1,2,2)%cosp
    endif

    ! jmin-kmin (bottom-front) edge boundary condition
    ! ================================================
    if ((BC_edge(2,1,1)%sort==-1).or.(BC_edge(2,1,1)%sort==-2)) then

       ! spherical coordinates
       ! ---------------------
       allocate(BC_edge(2,1,1)%i_r(ndx:nfx,ngh,ngh))
       allocate(BC_edge(2,1,1)%cosp(ndx:nfx,ngh,ngh))
       allocate(BC_edge(2,1,1)%sinp(ndx:nfx,ngh,ngh))
       allocate(BC_edge(2,1,1)%cost(ndx:nfx,ngh,ngh))
       allocate(BC_edge(2,1,1)%sint(ndx:nfx,ngh,ngh))
       allocate(BC_edge(2,1,1)%costcosp(ndx:nfx,ngh,ngh))
       allocate(BC_edge(2,1,1)%costsinp(ndx:nfx,ngh,ngh))
       allocate(BC_edge(2,1,1)%sintcosp(ndx:nfx,ngh,ngh))
       allocate(BC_edge(2,1,1)%sintsinp(ndx:nfx,ngh,ngh))
       
       do j=1,ngh
          do k=1,ngh
             do i=ndx,nfx
                if (is_curv3) then
                   xr=xc3(i,j,k)-BC_face(2,1)%xrc(i)
                   yr=yc3(i,j,k)-BC_face(2,1)%yrc(i)
                   zr=zc3(i,j,k)-BC_face(2,1)%zrc(i)
                else
                   if (is_curv) then
                      xr=xc(i,j)-BC_face(2,1)%xrc(i)
                      yr=yc(i,j)-BC_face(2,1)%yrc(i)
                   else
                      xr=x(i)-BC_face(2,1)%xrc(i)
                      yr=y(j)-BC_face(2,1)%yrc(i)
                   endif
                   zr=z(k)-BC_face(2,1)%zrc(i)
                endif
                ir=sqrt(xr**2+yr**2+zr**2)
                rp=sqrt(xr**2+yr**2)
                if (ir==0.0_wp)  ir=1.0_wp
                if (rp==0.0_wp)  rp=1.0_wp
                ir=1.0_wp/ir
                BC_edge(2,1,1)%i_r(i,j,k)=ir
                BC_edge(2,1,1)%cosp(i,j,k)=xr/rp
                BC_edge(2,1,1)%sinp(i,j,k)=yr/rp
                BC_edge(2,1,1)%cost(i,j,k)=zr*ir
                BC_edge(2,1,1)%sint(i,j,k)=rp*ir
             enddo
          enddo
       enddo
       BC_edge(2,1,1)%sintsinp=BC_edge(2,1,1)%sint*BC_edge(2,1,1)%sinp
       BC_edge(2,1,1)%sintcosp=BC_edge(2,1,1)%sint*BC_edge(2,1,1)%cosp
       BC_edge(2,1,1)%costsinp=BC_edge(2,1,1)%cost*BC_edge(2,1,1)%sinp
       BC_edge(2,1,1)%costcosp=BC_edge(2,1,1)%cost*BC_edge(2,1,1)%cosp
    endif

    ! jmin-kmax (bottom-back) edge boundary condition
    ! ===============================================
    if ((BC_edge(2,1,2)%sort==-1).or.(BC_edge(2,1,2)%sort==-2)) then

       ! spherical coordinates
       ! ---------------------
       allocate(BC_edge(2,1,2)%i_r(ndx:nfx,ngh,ngh))
       allocate(BC_edge(2,1,2)%cosp(ndx:nfx,ngh,ngh))
       allocate(BC_edge(2,1,2)%sinp(ndx:nfx,ngh,ngh))
       allocate(BC_edge(2,1,2)%cost(ndx:nfx,ngh,ngh))
       allocate(BC_edge(2,1,2)%sint(ndx:nfx,ngh,ngh))
       allocate(BC_edge(2,1,2)%costcosp(ndx:nfx,ngh,ngh))
       allocate(BC_edge(2,1,2)%costsinp(ndx:nfx,ngh,ngh))
       allocate(BC_edge(2,1,2)%sintcosp(ndx:nfx,ngh,ngh))
       allocate(BC_edge(2,1,2)%sintsinp(ndx:nfx,ngh,ngh))
       
       do j=1,ngh
          do k=nz-ngh+1,nz
             l=k-nz+ngh
             do i=ndx,nfx
                if (is_curv3) then
                   xr=xc3(i,j,k)-BC_face(2,1)%xrc(i)
                   yr=yc3(i,j,k)-BC_face(2,1)%yrc(i)
                   zr=zc3(i,j,k)-BC_face(2,1)%zrc(i)
                else
                   if (is_curv) then
                      xr=xc(i,j)-BC_face(2,1)%xrc(i)
                      yr=yc(i,j)-BC_face(2,1)%yrc(i)
                   else
                      xr=x(i)-BC_face(2,1)%xrc(i)
                      yr=y(j)-BC_face(2,1)%yrc(i)
                   endif
                   zr=z(k)-BC_face(2,1)%zrc(i)
                endif
                ir=sqrt(xr**2+yr**2+zr**2)
                rp=sqrt(xr**2+yr**2)
                if (ir==0.0_wp)  ir=1.0_wp
                if (rp==0.0_wp)  rp=1.0_wp
                ir=1.0_wp/ir
                BC_edge(2,1,2)%i_r(i,j,l)=ir
                BC_edge(2,1,2)%cosp(i,j,l)=xr/rp
                BC_edge(2,1,2)%sinp(i,j,l)=yr/rp
                BC_edge(2,1,2)%cost(i,j,l)=zr*ir
                BC_edge(2,1,2)%sint(i,j,l)=rp*ir
             enddo
          enddo
       enddo
       BC_edge(2,1,2)%sintsinp=BC_edge(2,1,2)%sint*BC_edge(2,1,2)%sinp
       BC_edge(2,1,2)%sintcosp=BC_edge(2,1,2)%sint*BC_edge(2,1,2)%cosp
       BC_edge(2,1,2)%costsinp=BC_edge(2,1,2)%cost*BC_edge(2,1,2)%sinp
       BC_edge(2,1,2)%costcosp=BC_edge(2,1,2)%cost*BC_edge(2,1,2)%cosp
    endif

    ! jmax-kmin (top-front) edge boundary condition
    ! =============================================
    if ((BC_edge(2,2,1)%sort==-1).or.(BC_edge(2,2,1)%sort==-2)) then

       ! spherical coordinates
       ! ---------------------
       allocate(BC_edge(2,2,1)%i_r(ndx:nfx,ngh,ngh))
       allocate(BC_edge(2,2,1)%cosp(ndx:nfx,ngh,ngh))
       allocate(BC_edge(2,2,1)%sinp(ndx:nfx,ngh,ngh))
       allocate(BC_edge(2,2,1)%cost(ndx:nfx,ngh,ngh))
       allocate(BC_edge(2,2,1)%sint(ndx:nfx,ngh,ngh))
       allocate(BC_edge(2,2,1)%costcosp(ndx:nfx,ngh,ngh))
       allocate(BC_edge(2,2,1)%costsinp(ndx:nfx,ngh,ngh))
       allocate(BC_edge(2,2,1)%sintcosp(ndx:nfx,ngh,ngh))
       allocate(BC_edge(2,2,1)%sintsinp(ndx:nfx,ngh,ngh))
       
       do j=ny-ngh+1,ny
          l=j-ny+ngh
          do k=1,ngh
             do i=ndx,nfx
                if (is_curv3) then
                   xr=xc3(i,j,k)-BC_face(2,2)%xrc(i)
                   yr=yc3(i,j,k)-BC_face(2,2)%yrc(i)
                   zr=zc3(i,j,k)-BC_face(2,2)%zrc(i)
                else
                   if (is_curv) then
                      xr=xc(i,j)-BC_face(2,2)%xrc(i)
                      yr=yc(i,j)-BC_face(2,2)%yrc(i)
                   else
                      xr=x(i)-BC_face(2,2)%xrc(i)
                      yr=y(j)-BC_face(2,2)%yrc(i)
                   endif
                   zr=z(k)-BC_face(2,2)%zrc(i)
                endif
                ir=sqrt(xr**2+yr**2+zr**2)
                rp=sqrt(xr**2+yr**2)
                if (ir==0.0_wp)  ir=1.0_wp
                if (rp==0.0_wp)  rp=1.0_wp
                ir=1.0_wp/ir
                BC_edge(2,2,1)%i_r(i,l,k)=ir
                BC_edge(2,2,1)%cosp(i,l,k)=xr/rp
                BC_edge(2,2,1)%sinp(i,l,k)=yr/rp
                BC_edge(2,2,1)%cost(i,l,k)=zr*ir
                BC_edge(2,2,1)%sint(i,l,k)=rp*ir
             enddo
          enddo
       enddo
       BC_edge(2,2,1)%sintsinp=BC_edge(2,2,1)%sint*BC_edge(2,2,1)%sinp
       BC_edge(2,2,1)%sintcosp=BC_edge(2,2,1)%sint*BC_edge(2,2,1)%cosp
       BC_edge(2,2,1)%costsinp=BC_edge(2,2,1)%cost*BC_edge(2,2,1)%sinp
       BC_edge(2,2,1)%costcosp=BC_edge(2,2,1)%cost*BC_edge(2,2,1)%cosp
    endif

    ! jmax-kmax (top-back) edge boundary condition
    ! ============================================
    if ((BC_edge(2,2,2)%sort==-1).or.(BC_edge(2,2,2)%sort==-2)) then

       ! spherical coordinates
       ! ---------------------
       allocate(BC_edge(2,2,2)%i_r(ndx:nfx,ngh,ngh))
       allocate(BC_edge(2,2,2)%cosp(ndx:nfx,ngh,ngh))
       allocate(BC_edge(2,2,2)%sinp(ndx:nfx,ngh,ngh))
       allocate(BC_edge(2,2,2)%cost(ndx:nfx,ngh,ngh))
       allocate(BC_edge(2,2,2)%sint(ndx:nfx,ngh,ngh))
       allocate(BC_edge(2,2,2)%costcosp(ndx:nfx,ngh,ngh))
       allocate(BC_edge(2,2,2)%costsinp(ndx:nfx,ngh,ngh))
       allocate(BC_edge(2,2,2)%sintcosp(ndx:nfx,ngh,ngh))
       allocate(BC_edge(2,2,2)%sintsinp(ndx:nfx,ngh,ngh))
       
       do j=ny-ngh+1,ny
          l=j-ny+ngh
          do k=nz-ngh+1,nz
             m=k-nz+ngh
             do i=ndx,nfx
                if (is_curv3) then
                   xr=xc3(i,j,k)-BC_face(2,2)%xrc(i)
                   yr=yc3(i,j,k)-BC_face(2,2)%yrc(i)
                   zr=zc3(i,j,k)-BC_face(2,2)%zrc(i)
                else
                   if (is_curv) then
                      xr=xc(i,j)-BC_face(2,2)%xrc(i)
                      yr=yc(i,j)-BC_face(2,2)%yrc(i)
                   else
                      xr=x(i)-BC_face(2,2)%xrc(i)
                      yr=y(j)-BC_face(2,2)%yrc(i)
                   endif
                   zr=z(k)-BC_face(2,2)%zrc(j)
                endif
                ir=sqrt(xr**2+yr**2+zr**2)
                rp=sqrt(xr**2+yr**2)
                if (ir==0.0_wp)  ir=1.0_wp
                if (rp==0.0_wp)  rp=1.0_wp
                ir=1.0_wp/ir
                BC_edge(2,2,2)%i_r(i,l,m)=ir
                BC_edge(2,2,2)%cosp(i,l,m)=xr/rp
                BC_edge(2,2,2)%sinp(i,l,m)=yr/rp
                BC_edge(2,2,2)%cost(i,l,m)=zr*ir
                BC_edge(2,2,2)%sint(i,l,m)=rp*ir
             enddo
          enddo
       enddo
       BC_edge(2,2,2)%sintsinp=BC_edge(2,2,2)%sint*BC_edge(2,2,2)%sinp
       BC_edge(2,2,2)%sintcosp=BC_edge(2,2,2)%sint*BC_edge(2,2,2)%cosp
       BC_edge(2,2,2)%costsinp=BC_edge(2,2,2)%cost*BC_edge(2,2,2)%sinp
       BC_edge(2,2,2)%costcosp=BC_edge(2,2,2)%cost*BC_edge(2,2,2)%cosp
    endif

    ! kmin-imin (front-left) edge boundary condition
    ! ==============================================
    if ((BC_edge(3,1,1)%sort==-1).or.(BC_edge(3,1,1)%sort==-2)) then

       ! spherical coordinates
       ! ---------------------
       allocate(BC_edge(3,1,1)%i_r(ngh,ndy:nfy,ngh))
       allocate(BC_edge(3,1,1)%cosp(ngh,ndy:nfy,ngh))
       allocate(BC_edge(3,1,1)%sinp(ngh,ndy:nfy,ngh))
       allocate(BC_edge(3,1,1)%cost(ngh,ndy:nfy,ngh))
       allocate(BC_edge(3,1,1)%sint(ngh,ndy:nfy,ngh))
       allocate(BC_edge(3,1,1)%costcosp(ngh,ndy:nfy,ngh))
       allocate(BC_edge(3,1,1)%costsinp(ngh,ndy:nfy,ngh))
       allocate(BC_edge(3,1,1)%sintcosp(ngh,ndy:nfy,ngh))
       allocate(BC_edge(3,1,1)%sintsinp(ngh,ndy:nfy,ngh))
       
       do i=1,ngh
          do k=1,ngh
             do j=ndy,nfy
                if (is_curv3) then
                   xr=xc3(i,j,k)-BC_face(1,1)%xrc(j)
                   yr=yc3(i,j,k)-BC_face(1,1)%yrc(j)
                   zr=zc3(i,j,k)-BC_face(1,1)%zrc(j)
                else
                   if (is_curv) then
                      xr=xc(i,j)-BC_face(1,1)%xrc(j)
                      yr=yc(i,j)-BC_face(1,1)%yrc(j)
                   else
                      xr=x(i)-BC_face(1,1)%xrc(j)
                      yr=y(j)-BC_face(1,1)%yrc(j)
                   endif
                   zr=z(k)-BC_face(1,1)%zrc(j)
                endif
                ir=sqrt(xr**2+yr**2+zr**2)
                rp=sqrt(xr**2+yr**2)
                if (ir==0.0_wp)  ir=1.0_wp
                if (rp==0.0_wp)  rp=1.0_wp
                ir=1.0_wp/ir
                BC_edge(3,1,1)%i_r(i,j,k)=ir
                BC_edge(3,1,1)%cosp(i,j,k)=xr/rp
                BC_edge(3,1,1)%sinp(i,j,k)=yr/rp
                BC_edge(3,1,1)%cost(i,j,k)=zr*ir
                BC_edge(3,1,1)%sint(i,j,k)=rp*ir
             enddo
          enddo
       enddo
       BC_edge(3,1,1)%sintsinp=BC_edge(3,1,1)%sint*BC_edge(3,1,1)%sinp
       BC_edge(3,1,1)%sintcosp=BC_edge(3,1,1)%sint*BC_edge(3,1,1)%cosp
       BC_edge(3,1,1)%costsinp=BC_edge(3,1,1)%cost*BC_edge(3,1,1)%sinp
       BC_edge(3,1,1)%costcosp=BC_edge(3,1,1)%cost*BC_edge(3,1,1)%cosp
    endif

    ! kmin-imax (front-right) edge boundary condition
    ! ===============================================
    if ((BC_edge(3,1,2)%sort==-1).or.(BC_edge(3,1,2)%sort==-2)) then

       ! spherical coordinates
       ! ---------------------
       allocate(BC_edge(3,1,2)%i_r(ngh,ndy:nfy,ngh))
       allocate(BC_edge(3,1,2)%cosp(ngh,ndy:nfy,ngh))
       allocate(BC_edge(3,1,2)%sinp(ngh,ndy:nfy,ngh))
       allocate(BC_edge(3,1,2)%cost(ngh,ndy:nfy,ngh))
       allocate(BC_edge(3,1,2)%sint(ngh,ndy:nfy,ngh))
       allocate(BC_edge(3,1,2)%costcosp(ngh,ndy:nfy,ngh))
       allocate(BC_edge(3,1,2)%costsinp(ngh,ndy:nfy,ngh))
       allocate(BC_edge(3,1,2)%sintcosp(ngh,ndy:nfy,ngh))
       allocate(BC_edge(3,1,2)%sintsinp(ngh,ndy:nfy,ngh))
       
       do i=nx-ngh+1,nx
          l=i-nx+ngh
          do k=1,ngh
             do j=ndy,nfy
                if (is_curv3) then
                   xr=xc3(i,j,k)-BC_face(1,2)%xrc(j)
                   yr=yc3(i,j,k)-BC_face(1,2)%yrc(j)
                   zr=zc3(i,j,k)-BC_face(1,2)%zrc(j)
                else
                   if (is_curv) then
                      xr=xc(i,j)-BC_face(1,2)%xrc(j)
                      yr=yc(i,j)-BC_face(1,2)%yrc(j)
                   else
                      xr=x(i)-BC_face(1,2)%xrc(j)
                      yr=y(j)-BC_face(1,2)%yrc(j)
                   endif
                   zr=z(k)-BC_face(1,2)%zrc(j)
                endif
                ir=sqrt(xr**2+yr**2+zr**2)
                rp=sqrt(xr**2+yr**2)
                if (ir==0.0_wp)  ir=1.0_wp
                if (rp==0.0_wp)  rp=1.0_wp
                ir=1.0_wp/ir
                BC_edge(3,1,2)%i_r(l,j,k)=ir
                BC_edge(3,1,2)%cosp(l,j,k)=xr/rp
                BC_edge(3,1,2)%sinp(l,j,k)=yr/rp
                BC_edge(3,1,2)%cost(l,j,k)=zr*ir
                BC_edge(3,1,2)%sint(l,j,k)=rp*ir
             enddo
          enddo
       enddo
       BC_edge(3,1,2)%sintsinp=BC_edge(3,1,2)%sint*BC_edge(3,1,2)%sinp
       BC_edge(3,1,2)%sintcosp=BC_edge(3,1,2)%sint*BC_edge(3,1,2)%cosp
       BC_edge(3,1,2)%costsinp=BC_edge(3,1,2)%cost*BC_edge(3,1,2)%sinp
       BC_edge(3,1,2)%costcosp=BC_edge(3,1,2)%cost*BC_edge(3,1,2)%cosp
    endif

    ! kmax-imin (back-left) edge boundary condition
    ! =============================================
    if ((BC_edge(3,2,1)%sort==-1).or.(BC_edge(3,2,1)%sort==-2)) then

       ! spherical coordinates
       ! ---------------------
       allocate(BC_edge(3,2,1)%i_r(ngh,ndy:nfy,ngh))
       allocate(BC_edge(3,2,1)%cosp(ngh,ndy:nfy,ngh))
       allocate(BC_edge(3,2,1)%sinp(ngh,ndy:nfy,ngh))
       allocate(BC_edge(3,2,1)%cost(ngh,ndy:nfy,ngh))
       allocate(BC_edge(3,2,1)%sint(ngh,ndy:nfy,ngh))
       allocate(BC_edge(3,2,1)%costcosp(ngh,ndy:nfy,ngh))
       allocate(BC_edge(3,2,1)%costsinp(ngh,ndy:nfy,ngh))
       allocate(BC_edge(3,2,1)%sintcosp(ngh,ndy:nfy,ngh))
       allocate(BC_edge(3,2,1)%sintsinp(ngh,ndy:nfy,ngh))
       
       do i=1,ngh
          do k=nz-ngh+1,nz
             l=k-nz+ngh
             do j=ndy,nfy
                if (is_curv3) then
                   xr=xc3(i,j,k)-BC_face(1,1)%xrc(j)
                   yr=yc3(i,j,k)-BC_face(1,1)%yrc(j)
                   zr=zc3(i,j,k)-BC_face(1,1)%zrc(j)
                else
                   if (is_curv) then
                      xr=xc(i,j)-BC_face(1,1)%xrc(j)
                      yr=yc(i,j)-BC_face(1,1)%yrc(j)
                   else
                      xr=x(i)-BC_face(1,1)%xrc(j)
                      yr=y(j)-BC_face(1,1)%yrc(j)
                   endif
                   zr=z(k)-BC_face(1,1)%zrc(j)
                endif
                ir=sqrt(xr**2+yr**2+zr**2)
                rp=sqrt(xr**2+yr**2)
                if (ir==0.0_wp)  ir=1.0_wp
                if (rp==0.0_wp)  rp=1.0_wp
                ir=1.0_wp/ir
                BC_edge(3,2,1)%i_r(i,j,l)=ir
                BC_edge(3,2,1)%cosp(i,j,l)=xr/rp
                BC_edge(3,2,1)%sinp(i,j,l)=yr/rp
                BC_edge(3,2,1)%cost(i,j,l)=zr*ir
                BC_edge(3,2,1)%sint(i,j,l)=rp*ir
             enddo
          enddo
       enddo
       BC_edge(3,2,1)%sintsinp=BC_edge(3,2,1)%sint*BC_edge(3,2,1)%sinp
       BC_edge(3,2,1)%sintcosp=BC_edge(3,2,1)%sint*BC_edge(3,2,1)%cosp
       BC_edge(3,2,1)%costsinp=BC_edge(3,2,1)%cost*BC_edge(3,2,1)%sinp
       BC_edge(3,2,1)%costcosp=BC_edge(3,2,1)%cost*BC_edge(3,2,1)%cosp
    endif

    ! kmax-imax (back-right) edge boundary condition
    ! ==============================================
    if ((BC_edge(3,2,2)%sort==-1).or.(BC_edge(3,2,2)%sort==-2)) then

       ! spherical coordinates
       ! ---------------------
       allocate(BC_edge(3,2,2)%i_r(ngh,ndy:nfy,ngh))
       allocate(BC_edge(3,2,2)%cosp(ngh,ndy:nfy,ngh))
       allocate(BC_edge(3,2,2)%sinp(ngh,ndy:nfy,ngh))
       allocate(BC_edge(3,2,2)%cost(ngh,ndy:nfy,ngh))
       allocate(BC_edge(3,2,2)%sint(ngh,ndy:nfy,ngh))
       allocate(BC_edge(3,2,2)%costcosp(ngh,ndy:nfy,ngh))
       allocate(BC_edge(3,2,2)%costsinp(ngh,ndy:nfy,ngh))
       allocate(BC_edge(3,2,2)%sintcosp(ngh,ndy:nfy,ngh))
       allocate(BC_edge(3,2,2)%sintsinp(ngh,ndy:nfy,ngh))
       
       do i=nx-ngh+1,nx
          l=i-nx+ngh
          do k=nz-ngh+1,nz
             m=k-nz+ngh
             do j=ndy,nfy
                if (is_curv3) then
                   xr=xc3(i,j,k)-BC_face(1,2)%xrc(j)
                   yr=yc3(i,j,k)-BC_face(1,2)%yrc(j)
                   zr=zc3(i,j,k)-BC_face(1,2)%zrc(j)
                else
                   if (is_curv) then
                      xr=xc(i,j)-BC_face(1,2)%xrc(j)
                      yr=yc(i,j)-BC_face(1,2)%yrc(j)
                   else
                      xr=x(i)-BC_face(1,2)%xrc(j)
                      yr=y(j)-BC_face(1,2)%yrc(j)
                   endif
                   zr=z(k)-BC_face(1,2)%zrc(j)
                endif
                ir=sqrt(xr**2+yr**2+zr**2)
                rp=sqrt(xr**2+yr**2)
                if (ir==0.0_wp)  ir=1.0_wp
                if (rp==0.0_wp)  rp=1.0_wp
                ir=1.0_wp/ir
                BC_edge(3,2,2)%i_r(l,j,m)=ir
                BC_edge(3,2,2)%cosp(l,j,m)=xr/rp
                BC_edge(3,2,2)%sinp(l,j,m)=yr/rp
                BC_edge(3,2,2)%cost(l,j,m)=zr*ir
                BC_edge(3,2,2)%sint(l,j,m)=rp*ir
             enddo
          enddo
       enddo
       BC_edge(3,2,2)%sintsinp=BC_edge(3,2,2)%sint*BC_edge(3,2,2)%sinp
       BC_edge(3,2,2)%sintcosp=BC_edge(3,2,2)%sint*BC_edge(3,2,2)%cosp
       BC_edge(3,2,2)%costsinp=BC_edge(3,2,2)%cost*BC_edge(3,2,2)%sinp
       BC_edge(3,2,2)%costcosp=BC_edge(3,2,2)%cost*BC_edge(3,2,2)%cosp
    endif

    ! imin-jmin-kmin (left-bottom-front) corner boundary condition
    ! ============================================================
    if ((BC_corner(1,1,1)%sort==-1).or.(BC_corner(1,1,1)%sort==-2)) then

       ! spherical coordinates
       ! ---------------------
       allocate(BC_corner(1,1,1)%i_r(ngh,ngh,ngh))
       allocate(BC_corner(1,1,1)%cosp(ngh,ngh,ngh))
       allocate(BC_corner(1,1,1)%sinp(ngh,ngh,ngh))
       allocate(BC_corner(1,1,1)%cost(ngh,ngh,ngh))
       allocate(BC_corner(1,1,1)%sint(ngh,ngh,ngh))
       allocate(BC_corner(1,1,1)%costcosp(ngh,ngh,ngh))
       allocate(BC_corner(1,1,1)%costsinp(ngh,ngh,ngh))
       allocate(BC_corner(1,1,1)%sintcosp(ngh,ngh,ngh))
       allocate(BC_corner(1,1,1)%sintsinp(ngh,ngh,ngh))
       
       do i=1,ngh
          do j=1,ngh
             do k=1,ngh
                if (is_curv3) then
                   xr=xc3(i,j,k)-BC_face(1,1)%xrc(j)
                   yr=yc3(i,j,k)-BC_face(1,1)%yrc(j)
                   zr=zc3(i,j,k)-BC_face(1,1)%zrc(j)
                else
                   if (is_curv) then
                      xr=xc(i,j)-BC_face(1,1)%xrc(j)
                      yr=yc(i,j)-BC_face(1,1)%yrc(j)
                   else
                      xr=x(i)-BC_face(1,1)%xrc(j)
                      yr=y(j)-BC_face(1,1)%yrc(j)
                   endif
                   zr=z(k)-BC_face(1,1)%zrc(j)
                endif
                ir=sqrt(xr**2+yr**2+zr**2)
                rp=sqrt(xr**2+yr**2)
                if (ir==0.0_wp)  ir=1.0_wp
                if (rp==0.0_wp)  rp=1.0_wp
                ir=1.0_wp/ir
                BC_corner(1,1,1)%i_r(i,j,k)=ir
                BC_corner(1,1,1)%cosp(i,j,k)=xr/rp
                BC_corner(1,1,1)%sinp(i,j,k)=yr/rp
                BC_corner(1,1,1)%cost(i,j,k)=zr*ir
                BC_corner(1,1,1)%sint(i,j,k)=rp*ir
             enddo
          enddo
       enddo
       BC_corner(1,1,1)%sintsinp=BC_corner(1,1,1)%sint*BC_corner(1,1,1)%sinp
       BC_corner(1,1,1)%sintcosp=BC_corner(1,1,1)%sint*BC_corner(1,1,1)%cosp
       BC_corner(1,1,1)%costsinp=BC_corner(1,1,1)%cost*BC_corner(1,1,1)%sinp
       BC_corner(1,1,1)%costcosp=BC_corner(1,1,1)%cost*BC_corner(1,1,1)%cosp
    endif

    ! imin-jmin-kmax (left-bottom-back) corner boundary condition
    ! ===========================================================
    if ((BC_corner(1,1,2)%sort==-1).or.(BC_corner(1,1,2)%sort==-2)) then

       ! spherical coordinates
       ! ---------------------
       allocate(BC_corner(1,1,2)%i_r(ngh,ngh,ngh))
       allocate(BC_corner(1,1,2)%cosp(ngh,ngh,ngh))
       allocate(BC_corner(1,1,2)%sinp(ngh,ngh,ngh))
       allocate(BC_corner(1,1,2)%cost(ngh,ngh,ngh))
       allocate(BC_corner(1,1,2)%sint(ngh,ngh,ngh))
       allocate(BC_corner(1,1,2)%costcosp(ngh,ngh,ngh))
       allocate(BC_corner(1,1,2)%costsinp(ngh,ngh,ngh))
       allocate(BC_corner(1,1,2)%sintcosp(ngh,ngh,ngh))
       allocate(BC_corner(1,1,2)%sintsinp(ngh,ngh,ngh))
       
       do i=1,ngh
          do j=1,ngh
             do k=nz-ngh+1,nz
                l=k-nz+ngh
                if (is_curv3) then
                   xr=xc3(i,j,k)-BC_face(1,1)%xrc(j)
                   yr=yc3(i,j,k)-BC_face(1,1)%yrc(j)
                   zr=zc3(i,j,k)-BC_face(1,1)%zrc(j)
                else
                   if (is_curv) then
                      xr=xc(i,j)-BC_face(1,1)%xrc(j)
                      yr=yc(i,j)-BC_face(1,1)%yrc(j)
                   else
                      xr=x(i)-BC_face(1,1)%xrc(j)
                      yr=y(j)-BC_face(1,1)%yrc(j)
                   endif
                   zr=z(k)-BC_face(1,1)%zrc(j)
                endif
                ir=sqrt(xr**2+yr**2+zr**2)
                rp=sqrt(xr**2+yr**2)
                if (ir==0.0_wp)  ir=1.0_wp
                if (rp==0.0_wp)  rp=1.0_wp
                ir=1.0_wp/ir
                BC_corner(1,1,2)%i_r(i,j,l)=ir
                BC_corner(1,1,2)%cosp(i,j,l)=xr/rp
                BC_corner(1,1,2)%sinp(i,j,l)=yr/rp
                BC_corner(1,1,2)%cost(i,j,l)=zr*ir
                BC_corner(1,1,2)%sint(i,j,l)=rp*ir
             enddo
          enddo
       enddo
       BC_corner(1,1,2)%sintsinp=BC_corner(1,1,2)%sint*BC_corner(1,1,2)%sinp
       BC_corner(1,1,2)%sintcosp=BC_corner(1,1,2)%sint*BC_corner(1,1,2)%cosp
       BC_corner(1,1,2)%costsinp=BC_corner(1,1,2)%cost*BC_corner(1,1,2)%sinp
       BC_corner(1,1,2)%costcosp=BC_corner(1,1,2)%cost*BC_corner(1,1,2)%cosp
    endif

    ! imin-jmax-kmin (left-top-front) corner boundary condition
    ! =========================================================
    if ((BC_corner(1,2,1)%sort==-1).or.(BC_corner(1,2,1)%sort==-2)) then

       ! spherical coordinates
       ! ---------------------
       allocate(BC_corner(1,2,1)%i_r(ngh,ngh,ngh))
       allocate(BC_corner(1,2,1)%cosp(ngh,ngh,ngh))
       allocate(BC_corner(1,2,1)%sinp(ngh,ngh,ngh))
       allocate(BC_corner(1,2,1)%cost(ngh,ngh,ngh))
       allocate(BC_corner(1,2,1)%sint(ngh,ngh,ngh))
       allocate(BC_corner(1,2,1)%costcosp(ngh,ngh,ngh))
       allocate(BC_corner(1,2,1)%costsinp(ngh,ngh,ngh))
       allocate(BC_corner(1,2,1)%sintcosp(ngh,ngh,ngh))
       allocate(BC_corner(1,2,1)%sintsinp(ngh,ngh,ngh))
       
       do i=1,ngh
          do j=ny-ngh+1,ny
             l=j-ny+ngh
             do k=1,ngh
                if (is_curv3) then
                   xr=xc3(i,j,k)-BC_face(1,1)%xrc(j)
                   yr=yc3(i,j,k)-BC_face(1,1)%yrc(j)
                   zr=zc3(i,j,k)-BC_face(1,1)%zrc(j)
                else
                   if (is_curv) then
                      xr=xc(i,j)-BC_face(1,1)%xrc(j)
                      yr=yc(i,j)-BC_face(1,1)%yrc(j)
                   else
                      xr=x(i)-BC_face(1,1)%xrc(j)
                      yr=y(j)-BC_face(1,1)%yrc(j)
                   endif
                   zr=z(k)-BC_face(1,1)%zrc(j)
                endif
                ir=sqrt(xr**2+yr**2+zr**2)
                rp=sqrt(xr**2+yr**2)
                if (rp==0.0_wp)  rp=1.0_wp
                if (ir==0.0_wp)  ir=1.0_wp
                ir=1.0_wp/ir
                BC_corner(1,2,1)%i_r(i,l,k)=ir
                BC_corner(1,2,1)%cosp(i,l,k)=xr/rp
                BC_corner(1,2,1)%sinp(i,l,k)=yr/rp
                BC_corner(1,2,1)%cost(i,l,k)=zr*ir
                BC_corner(1,2,1)%sint(i,l,k)=rp*ir
             enddo
          enddo
       enddo
       BC_corner(1,2,1)%sintsinp=BC_corner(1,2,1)%sint*BC_corner(1,2,1)%sinp
       BC_corner(1,2,1)%sintcosp=BC_corner(1,2,1)%sint*BC_corner(1,2,1)%cosp
       BC_corner(1,2,1)%costsinp=BC_corner(1,2,1)%cost*BC_corner(1,2,1)%sinp
       BC_corner(1,2,1)%costcosp=BC_corner(1,2,1)%cost*BC_corner(1,2,1)%cosp
    endif

    ! imin-jmax-kmax (left-top-back) corner boundary condition
    ! ========================================================
    if ((BC_corner(1,2,2)%sort==-1).or.(BC_corner(1,2,2)%sort==-2)) then

       ! spherical coordinates
       ! ---------------------
       allocate(BC_corner(1,2,2)%i_r(ngh,ngh,ngh))
       allocate(BC_corner(1,2,2)%cosp(ngh,ngh,ngh))
       allocate(BC_corner(1,2,2)%sinp(ngh,ngh,ngh))
       allocate(BC_corner(1,2,2)%cost(ngh,ngh,ngh))
       allocate(BC_corner(1,2,2)%sint(ngh,ngh,ngh))
       allocate(BC_corner(1,2,2)%costcosp(ngh,ngh,ngh))
       allocate(BC_corner(1,2,2)%costsinp(ngh,ngh,ngh))
       allocate(BC_corner(1,2,2)%sintcosp(ngh,ngh,ngh))
       allocate(BC_corner(1,2,2)%sintsinp(ngh,ngh,ngh))
       
       do i=1,ngh
          do j=ny-ngh+1,ny
             l=j-ny+ngh
             do k=nz-ngh+1,nz
                m=k-nz+ngh
                if (is_curv3) then
                   xr=xc3(i,j,k)-BC_face(1,1)%xrc(j)
                   yr=yc3(i,j,k)-BC_face(1,1)%yrc(j)
                   zr=zc3(i,j,k)-BC_face(1,1)%zrc(j)
                else
                   if (is_curv) then
                      xr=xc(i,j)-BC_face(1,1)%xrc(j)
                      yr=yc(i,j)-BC_face(1,1)%yrc(j)
                   else
                      xr=x(i)-BC_face(1,1)%xrc(j)
                      yr=y(j)-BC_face(1,1)%yrc(j)
                   endif
                   zr=z(k)-BC_face(1,1)%zrc(j)
                endif
                ir=sqrt(xr**2+yr**2+zr**2)
                rp=sqrt(xr**2+yr**2)
                if (ir==0.0_wp)  ir=1.0_wp
                if (rp==0.0_wp)  rp=1.0_wp
                ir=1.0_wp/ir
                BC_corner(1,2,2)%i_r(i,l,m)=ir
                BC_corner(1,2,2)%cosp(i,l,m)=xr/rp
                BC_corner(1,2,2)%sinp(i,l,m)=yr/rp
                BC_corner(1,2,2)%cost(i,l,m)=zr*ir
                BC_corner(1,2,2)%sint(i,l,m)=rp*ir
             enddo
          enddo
       enddo
       BC_corner(1,2,2)%sintsinp=BC_corner(1,2,2)%sint*BC_corner(1,2,2)%sinp
       BC_corner(1,2,2)%sintcosp=BC_corner(1,2,2)%sint*BC_corner(1,2,2)%cosp
       BC_corner(1,2,2)%costsinp=BC_corner(1,2,2)%cost*BC_corner(1,2,2)%sinp
       BC_corner(1,2,2)%costcosp=BC_corner(1,2,2)%cost*BC_corner(1,2,2)%cosp
    endif

    ! imax-jmin-kmin (right-bottom-front) corner boundary condition
    ! =============================================================
    if ((BC_corner(2,1,1)%sort==-1).or.(BC_corner(2,1,1)%sort==-2)) then

       ! spherical coordinates
       ! ---------------------
       allocate(BC_corner(2,1,1)%i_r(ngh,ngh,ngh))
       allocate(BC_corner(2,1,1)%cosp(ngh,ngh,ngh))
       allocate(BC_corner(2,1,1)%sinp(ngh,ngh,ngh))
       allocate(BC_corner(2,1,1)%cost(ngh,ngh,ngh))
       allocate(BC_corner(2,1,1)%sint(ngh,ngh,ngh))
       allocate(BC_corner(2,1,1)%costcosp(ngh,ngh,ngh))
       allocate(BC_corner(2,1,1)%costsinp(ngh,ngh,ngh))
       allocate(BC_corner(2,1,1)%sintcosp(ngh,ngh,ngh))
       allocate(BC_corner(2,1,1)%sintsinp(ngh,ngh,ngh))
       
       do i=nx-ngh+1,nx
          l=i-nx+ngh
          do j=1,ngh
             do k=1,ngh
                if (is_curv3) then
                   xr=xc3(i,j,k)-BC_face(1,2)%xrc(j)
                   yr=yc3(i,j,k)-BC_face(1,2)%yrc(j)
                   zr=zc3(i,j,k)-BC_face(1,2)%zrc(j)
                else
                   if (is_curv) then
                      xr=xc(i,j)-BC_face(1,2)%xrc(j)
                      yr=yc(i,j)-BC_face(1,2)%yrc(j)
                   else
                      xr=x(i)-BC_face(1,2)%xrc(j)
                      yr=y(j)-BC_face(1,2)%yrc(j)
                   endif
                   zr=z(k)-BC_face(1,2)%zrc(j)
                endif
                ir=sqrt(xr**2+yr**2+zr**2)
                rp=sqrt(xr**2+yr**2)
                if (ir==0.0_wp)  ir=1.0_wp
                if (rp==0.0_wp)  rp=1.0_wp
                ir=1.0_wp/ir
                BC_corner(2,1,1)%i_r(l,j,k)=ir
                BC_corner(2,1,1)%cosp(l,j,k)=xr/rp
                BC_corner(2,1,1)%sinp(l,j,k)=yr/rp
                BC_corner(2,1,1)%cost(l,j,k)=zr*ir
                BC_corner(2,1,1)%sint(l,j,k)=rp*ir
             enddo
          enddo
       enddo
       BC_corner(2,1,1)%sintsinp=BC_corner(2,1,1)%sint*BC_corner(2,1,1)%sinp
       BC_corner(2,1,1)%sintcosp=BC_corner(2,1,1)%sint*BC_corner(2,1,1)%cosp
       BC_corner(2,1,1)%costsinp=BC_corner(2,1,1)%cost*BC_corner(2,1,1)%sinp
       BC_corner(2,1,1)%costcosp=BC_corner(2,1,1)%cost*BC_corner(2,1,1)%cosp
    endif

    ! imax-jmin-kmax (right-bottom-back) corner boundary condition
    ! ============================================================
    if ((BC_corner(2,1,2)%sort==-1).or.(BC_corner(2,1,2)%sort==-2)) then

       ! spherical coordinates
       ! ---------------------
       allocate(BC_corner(2,1,2)%i_r(ngh,ngh,ngh))
       allocate(BC_corner(2,1,2)%cosp(ngh,ngh,ngh))
       allocate(BC_corner(2,1,2)%sinp(ngh,ngh,ngh))
       allocate(BC_corner(2,1,2)%cost(ngh,ngh,ngh))
       allocate(BC_corner(2,1,2)%sint(ngh,ngh,ngh))
       allocate(BC_corner(2,1,2)%costcosp(ngh,ngh,ngh))
       allocate(BC_corner(2,1,2)%costsinp(ngh,ngh,ngh))
       allocate(BC_corner(2,1,2)%sintcosp(ngh,ngh,ngh))
       allocate(BC_corner(2,1,2)%sintsinp(ngh,ngh,ngh))
       
       do i=nx-ngh+1,nx
          l=i-nx+ngh
          do j=1,ngh
             do k=nz-ngh+1,nz
                m=k-nz+ngh
                if (is_curv3) then
                   xr=xc3(i,j,k)-BC_face(1,2)%xrc(j)
                   yr=yc3(i,j,k)-BC_face(1,2)%yrc(j)
                   zr=zc3(i,j,k)-BC_face(1,2)%zrc(j)
                else
                   if (is_curv) then
                      xr=xc(i,j)-BC_face(1,2)%xrc(j)
                      yr=yc(i,j)-BC_face(1,2)%yrc(j)
                   else
                      xr=x(i)-BC_face(1,2)%xrc(j)
                      yr=y(j)-BC_face(1,2)%yrc(j)
                   endif
                   zr=z(k)-BC_face(1,2)%zrc(j)
                endif
                ir=sqrt(xr**2+yr**2+zr**2)
                rp=sqrt(xr**2+yr**2)
                if (ir==0.0_wp)  ir=1.0_wp
                if (rp==0.0_wp)  rp=1.0_wp
                ir=1.0_wp/ir
                BC_corner(2,1,2)%i_r(l,j,m)=ir
                BC_corner(2,1,2)%cosp(l,j,m)=xr/rp
                BC_corner(2,1,2)%sinp(l,j,m)=yr/rp
                BC_corner(2,1,2)%cost(l,j,m)=zr*ir
                BC_corner(2,1,2)%sint(l,j,m)=rp*ir
             enddo
          enddo
       enddo
       BC_corner(2,1,2)%sintsinp=BC_corner(2,1,2)%sint*BC_corner(2,1,2)%sinp
       BC_corner(2,1,2)%sintcosp=BC_corner(2,1,2)%sint*BC_corner(2,1,2)%cosp
       BC_corner(2,1,2)%costsinp=BC_corner(2,1,2)%cost*BC_corner(2,1,2)%sinp
       BC_corner(2,1,2)%costcosp=BC_corner(2,1,2)%cost*BC_corner(2,1,2)%cosp
    endif

    ! imax-jmax-kmin (right-top-front) corner boundary condition
    ! ============================================================
    if ((BC_corner(2,2,1)%sort==-1).or.(BC_corner(2,2,1)%sort==-2)) then

       ! spherical coordinates
       ! ---------------------
       allocate(BC_corner(2,2,1)%i_r(ngh,ngh,ngh))
       allocate(BC_corner(2,2,1)%cosp(ngh,ngh,ngh))
       allocate(BC_corner(2,2,1)%sinp(ngh,ngh,ngh))
       allocate(BC_corner(2,2,1)%cost(ngh,ngh,ngh))
       allocate(BC_corner(2,2,1)%sint(ngh,ngh,ngh))
       allocate(BC_corner(2,2,1)%costcosp(ngh,ngh,ngh))
       allocate(BC_corner(2,2,1)%costsinp(ngh,ngh,ngh))
       allocate(BC_corner(2,2,1)%sintcosp(ngh,ngh,ngh))
       allocate(BC_corner(2,2,1)%sintsinp(ngh,ngh,ngh))
       
       do i=nx-ngh+1,nx
          l=i-nx+ngh
          do j=ny-ngh+1,ny
             m=j-ny+ngh
             do k=1,ngh
                if (is_curv3) then
                   xr=xc3(i,j,k)-BC_face(1,2)%xrc(j)
                   yr=yc3(i,j,k)-BC_face(1,2)%yrc(j)
                   zr=zc3(i,j,k)-BC_face(1,2)%zrc(j)
                else
                   if (is_curv) then
                      xr=xc(i,j)-BC_face(1,2)%xrc(j)
                      yr=yc(i,j)-BC_face(1,2)%yrc(j)
                   else
                      xr=x(i)-BC_face(1,2)%xrc(j)
                      yr=y(j)-BC_face(1,2)%yrc(j)
                   endif
                   zr=z(k)-BC_face(1,2)%zrc(j)
                endif
                ir=sqrt(xr**2+yr**2+zr**2)
                rp=sqrt(xr**2+yr**2)
                if (ir==0.0_wp)  ir=1.0_wp
                if (rp==0.0_wp)  rp=1.0_wp
                ir=1.0_wp/ir
                BC_corner(2,2,1)%i_r(l,m,k)=ir
                BC_corner(2,2,1)%cosp(l,m,k)=xr/rp
                BC_corner(2,2,1)%sinp(l,m,k)=yr/rp
                BC_corner(2,2,1)%cost(l,m,k)=zr*ir
                BC_corner(2,2,1)%sint(l,m,k)=rp*ir
             enddo
          enddo
       enddo
       BC_corner(2,2,1)%sintsinp=BC_corner(2,2,1)%sint*BC_corner(2,2,1)%sinp
       BC_corner(2,2,1)%sintcosp=BC_corner(2,2,1)%sint*BC_corner(2,2,1)%cosp
       BC_corner(2,2,1)%costsinp=BC_corner(2,2,1)%cost*BC_corner(2,2,1)%sinp
       BC_corner(2,2,1)%costcosp=BC_corner(2,2,1)%cost*BC_corner(2,2,1)%cosp
    endif

    ! imax-jmax-kmax (right-top-back) corner boundary condition
    ! ============================================================
    if ((BC_corner(2,2,2)%sort==-1).or.(BC_corner(2,2,2)%sort==-2)) then

       ! spherical coordinates
       ! ---------------------
       allocate(BC_corner(2,2,2)%i_r(ngh,ngh,ngh))
       allocate(BC_corner(2,2,2)%cosp(ngh,ngh,ngh))
       allocate(BC_corner(2,2,2)%sinp(ngh,ngh,ngh))
       allocate(BC_corner(2,2,2)%cost(ngh,ngh,ngh))
       allocate(BC_corner(2,2,2)%sint(ngh,ngh,ngh))
       allocate(BC_corner(2,2,2)%costcosp(ngh,ngh,ngh))
       allocate(BC_corner(2,2,2)%costsinp(ngh,ngh,ngh))
       allocate(BC_corner(2,2,2)%sintcosp(ngh,ngh,ngh))
       allocate(BC_corner(2,2,2)%sintsinp(ngh,ngh,ngh))
       
       do i=nx-ngh+1,nx
          l=i-nx+ngh
          do j=ny-ngh+1,ny
             m=j-ny+ngh
             do k=nz-ngh+1,nz
                n=k-nz+ngh
                if (is_curv3) then
                   xr=xc3(i,j,k)-BC_face(1,2)%xrc(j)
                   yr=yc3(i,j,k)-BC_face(1,2)%yrc(j)
                   zr=zc3(i,j,k)-BC_face(1,2)%zrc(j)
                else
                   if (is_curv) then
                      xr=xc(i,j)-BC_face(1,2)%xrc(j)
                      yr=yc(i,j)-BC_face(1,2)%yrc(j)
                   else
                      xr=x(i)-BC_face(1,2)%xrc(j)
                      yr=y(j)-BC_face(1,2)%yrc(j)
                   endif
                   zr=z(k)-BC_face(1,2)%zrc(j)
                endif
                ir=sqrt(xr**2+yr**2+zr**2)
                rp=sqrt(xr**2+yr**2)
                if (ir==0.0_wp)  ir=1.0_wp
                if (rp==0.0_wp)  rp=1.0_wp
                ir=1.0_wp/ir
                BC_corner(2,2,2)%i_r(l,m,n)=ir
                BC_corner(2,2,2)%cosp(l,m,n)=xr/rp
                BC_corner(2,2,2)%sinp(l,m,n)=yr/rp
                BC_corner(2,2,2)%cost(l,m,n)=zr*ir
                BC_corner(2,2,2)%sint(l,m,n)=rp*ir
             enddo
          enddo
       enddo
       BC_corner(2,2,2)%sintsinp=BC_corner(2,2,2)%sint*BC_corner(2,2,2)%sinp
       BC_corner(2,2,2)%sintcosp=BC_corner(2,2,2)%sint*BC_corner(2,2,2)%cosp
       BC_corner(2,2,2)%costsinp=BC_corner(2,2,2)%cost*BC_corner(2,2,2)%sinp
       BC_corner(2,2,2)%costcosp=BC_corner(2,2,2)%cost*BC_corner(2,2,2)%cosp
    endif
    
  end subroutine init_param_TD3d
  
  !===============================================================================
  module subroutine free_bc_TD3d
  !===============================================================================
    !> Closing of Tam & Dong's BC: deallocations of memory
  !===============================================================================
    implicit none
    !-----------------------------------------------------------------------------
    integer :: l,m,n
    !-----------------------------------------------------------------------------

    ! Free BC arrays for faces
    ! ========================
    do l=1,3
       do m=1,2
          if (is_bc_TD(l,m)) then
             deallocate(BC_face(l,m)%i_r,BC_face(l,m)%U0)
             deallocate(BC_face(l,m)%cosp,BC_face(l,m)%sinp)
             deallocate(BC_face(l,m)%cost,BC_face(l,m)%sint)
             deallocate(BC_face(l,m)%costcosp,BC_face(l,m)%costsinp)
             deallocate(BC_face(l,m)%sintcosp,BC_face(l,m)%sintsinp)
          endif
       enddo
    enddo

    ! Free BC arrays for edges
    ! ========================
    do l=1,3
       do m=1,2
          do n=1,2
             if (BC_edge(l,m,n)%sort==-1) then
                deallocate(BC_edge(l,m,n)%i_r)
                deallocate(BC_edge(l,m,n)%cosp,BC_edge(l,m,n)%sinp)
                deallocate(BC_edge(l,m,n)%cost,BC_edge(l,m,n)%sint)
                deallocate(BC_edge(l,m,n)%costcosp,BC_edge(l,m,n)%costsinp)
                deallocate(BC_edge(l,m,n)%sintcosp,BC_edge(l,m,n)%sintsinp)
             endif
          enddo
       enddo
    enddo

    ! Free BC arrays for corners
    ! ========================
    do l=1,2
       do m=1,2
          do n=1,2
             if (BC_corner(l,m,n)%sort==-1) then
                deallocate(BC_corner(l,m,n)%i_r)
                deallocate(BC_corner(l,m,n)%cosp,BC_corner(l,m,n)%sinp)
                deallocate(BC_corner(l,m,n)%cost,BC_corner(l,m,n)%sint)
                deallocate(BC_corner(l,m,n)%costcosp,BC_corner(l,m,n)%costsinp)
                deallocate(BC_corner(l,m,n)%sintcosp,BC_corner(l,m,n)%sintsinp)
             endif
          enddo
       enddo
    enddo

  end subroutine free_bc_TD3d
  
end submodule smod_init_TamDong3
