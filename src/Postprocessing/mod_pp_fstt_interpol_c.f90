!==============================================================================
module mod_pp_fstt_interpol_c
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
  subroutine init_fstt_interpol_c
  !============================================================================
    !> Initialisation of interpolation for FSTT PP
    !============================================================================
    use mod_utils
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i,j,k,ngh_pp_temp,ncount,ngx_,ngy_,ngz_,ib,ideb
    real(wp) :: Lz_g ! Length in spanwise direction
    real(wp) :: step
    real(wp), dimension(:), allocatable :: z_interp_full
    ! -------------------------------------------------------------------------

    if ((coord(1)==0).and.(coord(3)==0)) write(*,18) iblc_pp,nit_interp,nj_interp,nkt_interp
18  format(5x,"Number of points for interpolation in block ",i0,":   ",i0," x ",i0," x ",i0,"  (i x j x k)")
    if ((coord(1)==0).and.(coord(3)==0)) write(*,19) iblc_pp,nil_interp,nj_interp,nk_interp
19  format(5x,"Number of points for interp per proc in block ",i0,": ",i0," x ",i0," x ",i0,"  (i x j x k)")

    idem_f_pp = bl_glob(nob(iproc))%snapshot(nsr)%ind_i1-1         ! Index imin of the full grid concerned by pp
    nx_full=0; nx_f_deb=0
    do ib=1,nbloc
       if (ib.eq.iblc_pp-nblr(1)+1) nx_f_deb = nx_full
       nx_full = nx_full + bl_glob(ib)%ni
    enddo
    if (nx_f_deb.eq.0) nx_f_deb=idem_f_pp
    iend_f_pp = nx_full - (bl_glob(nbloc_r)%ni - bl_glob(nob(iproc))%snapshot(nsr)%ind_i2)  ! Index imax

    ! Reading all grid
    ! ----------------
    allocate(xgc_f(1:nx_full,1:bl_glob(1)%nj))
    allocate(ygc_f(1:nx_full,1:bl_glob(1)%nj))
    ncount = 0
    do ib=nblr(1),nblr(nbloc_r)
       ! open(194,file=trim(dirDATA)//'grid_bl'//trim(numchar(ib))//"_proc"//trim(numchar(iproc))//'.bin', &
       !    form='unformatted',status='unknown')
       open(194,file=trim(dirDATA)//'grid_bl'//trim(numchar(ib))//'.bin', &
          form='unformatted',status='unknown')
       rewind(194)
       read(194) ngx_
       read(194) ngy_
       read(194) ngz_
       read(194) ((xgc_f(i+ncount,j),i=1,ngx_),j=1,ngy_)
       read(194) ((ygc_f(i+ncount,j),i=1,ngx_),j=1,ngy_)
       close(194)
       ncount = ncount + ngx_
    enddo

    ! Length of the spanwise domain
    ! -----------------------------
    Lz_g = zg(ngz) - zg(1)
    allocate(z_interp_full(nkt_interp))
    step = Lz_g/(nkt_interp-1)
    z_interp_full(1) = zg(1)
    do k=2,nkt_interp
       z_interp_full(k) = z_interp_full(k-1) + step
    enddo
    if (iproc.eq.0) print *,"h_kern*step ~>",h_kern*step

    ! Determination of ghost points and index bounds
    ! ==============================================
    ! Calculated based on the kernel
    prct_kernel = 0.025_wp
    ngh_pp_temp=0
    ! borns in i direction
    ideb = nx_f_deb + coord(1)*nx ! initial position in the full file of coordinates
    do j=2,nj_interp
       do i=1,nil_interp
          bk_ic(i,j,1) = i
          ncount = 0
          do while ((exp(-0.5_wp*((xgc_f(ndivx*bk_ic(i,j,1)+ideb,ndivy*j)-xgc_f(ndivx*i+ideb,ndivy*j))**2 + &
                     (ygc_f(ndivx*bk_ic(i,j,1)+ideb,ndivy*j)-ygc_f(ndivx*i+ideb,ndivy*j))**2)**0.5/(h_kern*step))>prct_kernel) &
                     .and.(ndivx*bk_ic(i,j,1)+ideb.gt.idem_f_pp))
             bk_ic(i,j,1) = bk_ic(i,j,1) - 1
             ncount = ncount+1
          enddo
          ! if (ndivx*bk_ic(i,j,1)+ideb.eq.idem_f_pp) print *,"min",iproc,i,j,bk_ic(i,j,1),nx_f_deb,ideb
          bk_ic(i,j,2) = i
          ncount = 0
          do while ((exp(-0.5_wp*((xgc_f(ndivx*bk_ic(i,j,2)+ideb,ndivy*j)-xgc_f(ndivx*i+ideb,ndivy*j))**2 + &
                     (ygc_f(ndivx*bk_ic(i,j,2)+ideb,ndivy*j)-ygc_f(ndivx*i+ideb,ndivy*j))**2)**0.5/(h_kern*step))>prct_kernel) &
                     .and.(ndivx*bk_ic(i,j,2)+ideb.lt.iend_f_pp))
             bk_ic(i,j,2) = bk_ic(i,j,2) + 1
             ncount = ncount+1
          enddo
          ! if (ndivx*bk_ic(i,j,2)+ideb.eq.iend_f_pp) print *,"max",iproc,ideb,coord(1),ndivx*bk_ic(i,j,2)+ideb,iend_f_pp
       enddo
       ngh_pp_temp = max(ngh_pp_temp,abs(bk_ic(1,j,1)-1),bk_ic(1,j,2)-nil_interp)
    enddo

    ! Communication of ghost points
    call MPI_ALLREDUCE(ngh_pp_temp,ngh_pp,1,MPI_INT,MPI_MAX,COMM_global,info)

    ! borns in j direction
    do i=1,nil_interp
       do j=2,nj_interp
          bk_jc(i,j,1) = j
          do while ((exp(-0.5_wp*((xgc_f(ndivx*i+ideb,ndivy*bk_jc(i,j,1))-xgc_f(ndivx*i+ideb,ndivy*j))**2 + &
                     (ygc_f(ndivx*i+ideb,ndivy*bk_jc(i,j,1))-ygc_f(ndivx*i+ideb,ndivy*j))**2)**0.5/(h_kern*step))>prct_kernel) &
                     .and.(bk_jc(i,j,1).gt.2))
             bk_jc(i,j,1) = bk_jc(i,j,1) - 1
          enddo
          bk_jc(i,j,2) = j
          do while ((exp(-0.5_wp*((xgc_f(ndivx*i+ideb,ndivy*bk_jc(i,j,2))-xgc_f(ndivx*i+ideb,ndivy*j))**2 + &
                     (ygc_f(ndivx*i+ideb,ndivy*bk_jc(i,j,2))-ygc_f(ndivx*i+ideb,ndivy*j))**2)**0.5/(h_kern*step))>prct_kernel) &
                     .and.(bk_jc(i,j,2).lt.nj_interp))
             bk_jc(i,j,2) = bk_jc(i,j,2) + 1
          enddo
       enddo
    enddo

    ! borns in k direction
    nkernel_k=0
    do while (exp(-0.5_wp*((nkernel_k*step)**2)**0.5/(h_kern*step))>prct_kernel)
       nkernel_k = nkernel_k+1
    enddo

    if (iproc==0) write(*,20) ngh_pp,nkernel_k
20  format(5x,"Ghost points for pp: ngh_pp = ",i0,", nkernel_k = ",i0)
    if (iproc==0) write(*,21) (2*ngh_pp+1),bk_jc(1,2,2)-bk_jc(1,2,1)+1,(2*nkernel_k+1)
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

    ! ! Creation of the grid for interpolation
    ! ! ======================================

    ! Allocation of grid for interpolation
    allocate(xc_interp(0:nil_interp+1,1:nj_interp),&
             yc_interp(0:nil_interp+1,1:nj_interp),&
             z_interp(0:nk_interp+1))

    ! Creation of local grid for interpolation
    ! ----------------------------------------
    ! i direction
    ! -----------
    ideb = nx_f_deb + coord(1)*nx ! initial position in the full file of coordinates
    do i=ndx_pp-1,nfx_pp+1
       do j=2,nj_interp
          xc_interp(i,j) = xgc_f(ndivx*i+ideb,ndivy*j)
       enddo
    enddo
    ! j=1
    do i=ndx_pp-1,nfx_pp+1
       xc_interp(i,1) = xgc_f(ndivx*i+ideb,1)
    enddo
    ! Artificial extension for inverse interpolation
    if (is_boundary(1,1)) then
       do j=2,nj_interp
          step = xc_interp(2,j) - xc_interp(1,j)
          xc_interp(0,j) = xc_interp(1,j) - step
       enddo
       ! j=1
       step = xc_interp(2,1) - xc_interp(1,1)
       xc_interp(0,1) = xc_interp(1,1) - step
    endif
    if (is_boundary(1,2)) then
       do j=2,nj_interp
          step = xc_interp(nil_interp,j) - xc_interp(nil_interp-1,j)
          xc_interp(nil_interp+1,j) = xc_interp(nil_interp,j) + step
       enddo
       ! j=1
       step = xc_interp(nil_interp,1) - xc_interp(nil_interp-1,1)
       xc_interp(nil_interp+1,1) = xc_interp(nil_interp,1) + step
    endif

    ! j direction
    ! -----------
    do i=ndx_pp-1,nfx_pp+1
       do j=2,nj_interp
          yc_interp(i,j) = ygc_f(ndivx*i+ideb,ndivy*j)
       enddo
    enddo
    ! j=1
    do i=ndx_pp-1,nfx_pp+1
       yc_interp(i,1) = ygc_f(ndivx*i+ideb,1)
    enddo
    ! Artificial extension for inverse interpolation
    if (is_boundary(1,1)) then
       do j=2,nj_interp
          step = yc_interp(2,j) - yc_interp(1,j)
          yc_interp(0,j) = yc_interp(1,j) - step
       enddo
       ! j=1
       step = yc_interp(2,1) - yc_interp(1,1)
       yc_interp(0,1) = yc_interp(1,1) - step
    endif
    if (is_boundary(1,2)) then
       do j=2,nj_interp
          step = yc_interp(nil_interp,j) - yc_interp(nil_interp-1,j)
          yc_interp(nil_interp+1,j) = yc_interp(nil_interp,j) + step
       enddo
       ! j=1
       step = yc_interp(nil_interp,1) - yc_interp(nil_interp-1,1)
       yc_interp(nil_interp+1,1) = yc_interp(nil_interp,1) + step
    endif

    ! k direction
    ! -----------
    do k=1,nk_interp
       z_interp(k) = z_interp_full(k + coord(3)*nk_interp)
    enddo
    ! Artificial extension for inverse interpolation
    step = z_interp(2) - z_interp(1)
    z_interp(0) = z_interp(1) - step
    z_interp(nk_interp+1) = z_interp(nk_interp) + step


    ! Writting of the interpolated grid
    ! =================================
    if ((coord(1)==0).and.(coord(3)==0)) then
       ! Writting
       open(194,file='grid_interp_bl'//trim(numchar(iblc_pp))//'.bin',form='unformatted',status='unknown')
       rewind(194)
       write(194) nit_interp
       write(194) nj_interp
       write(194) nkt_interp
       write(194) ((xgc_f(nx_f_deb+i*ndivx,j*ndivy),i=1,nit_interp),j=1,nj_interp)
       write(194) ((ygc_f(nx_f_deb+i*ndivx,j*ndivy),i=1,nit_interp),j=1,nj_interp)
       write(194) (z_interp_full(k),k=1,nkt_interp)
    endif

    deallocate(z_interp_full)

  end subroutine init_fstt_interpol_c

  !============================================================================
  subroutine fstt_interpol_3d_c(id_case)
  !============================================================================
    !> Interpolation (not really...) for FSTT LE
    !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer, intent(in)  :: id_case
    ! -------------------------------------------------------------------------
    integer :: i,j,k,k0,k1,m
    real(wp) :: zd
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

      select case(id_case)
      case(1)
         ! "Interpolation" of stats
         ! ========================
         do i=1,nil_interp
            do j=2,nj_interp
               do m=1,23
                  stats_interp(i,j,m) = stats_proc(i*ndivx,j*ndivy,m)
               enddo
            enddo
         enddo
         ! j=1
         do i=1,nil_interp
            do m=1,23
               stats_interp(i,1,m) = stats_proc(i*ndivx,1,m)
            enddo
         enddo
      case(2)
         ! "Interpolation" of u tangential
         ! ===============================
         do i=1,nil_interp
            do j=2,nj_interp
               uut_interp(i,j,k) = uut(i*ndivx,j*ndivy,k0)*(1.0_wp-zd) + uut(i*ndivx,j*ndivy,k1)*zd
            enddo
         enddo
         ! j=1
         do i=1,nil_interp
            uut_interp(i,1,k) = uut(i*ndivx,1,k0)*(1.0_wp-zd) + uut(i*ndivx,1,k1)*zd
         enddo

         ! "Interpolation" of u normal
         ! ===========================
         do i=1,nil_interp
            do j=2,nj_interp
               uun_interp(i,j,k) = uun(i*ndivx,j*ndivy,k0)*(1.0_wp-zd) + uun(i*ndivx,j*ndivy,k1)*zd
            enddo
         enddo
         ! j=1
         do i=1,nil_interp
            uun_interp(i,1,k) = uun(i*ndivx,1,k0)*(1.0_wp-zd) + uun(i*ndivx,1,k1)*zd
         enddo

         ! "Interpolation" of w
         ! ====================
         do i=1,nil_interp
            do j=2,nj_interp
               ww_interp(i,j,k) = ww(i*ndivx,j*ndivy,k0)*(1.0_wp-zd) + ww(i*ndivx,j*ndivy,k1)*zd
            enddo
         enddo
         ! j=1
         do i=1,nil_interp
            ww_interp(i,1,k) = ww(i*ndivx,1,k0)*(1.0_wp-zd) + ww(i*ndivx,1,k1)*zd
         enddo

      case default
         call mpistop('bad choice for interpolation', 0)
      end select

      if (id_case.eq.1) exit loopz
   enddo loopz


  end subroutine fstt_interpol_3d_c

  !============================================================================
  subroutine fstt_interpol_3d_inv_c
  !============================================================================
    !> Interpolation (not really...) for FSTT LE
    !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i,j,k,i0,j0,k0,i1,j1,k1,i_,j_
    real(wp) :: zd,cc,cc00,cc01,cc10,cc11,d00,d01,d10,d11,dtot
    ! -------------------------------------------------------------------------

    ! if (iproc.eq.0) then

      k=3
      i=17
      j=13
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
      if (zd.lt.0.0_wp) zd = 0.0_wp
      if (zd.gt.1.0_wp) zd = 1.0_wp

      do j=1,ny
         if (j.eq.1) then
            j0 = j; j1 = j
         elseif (REAL(j/REAL(ndivy)).eq.REAL(j/ndivy)) then
            j0 = j/ndivy; j1 = j/ndivy
         else
            loopj_: do j_=1,ndivy
               if (REAL((j-j_)/REAL(ndivy)).eq.REAL((j-j_)/ndivy)) exit loopj_
            enddo loopj_
            j0 = ((j-j_)/ndivy); j1 = j0+1
         endif

         if (j1.gt.nj_interp) then
            j1 = nj_interp
            j0 = nj_interp
            ! yd = 0.0_wp
         ! else
         !    j0 = j1-1
         !    yd = (y(j)-y_interp(j0))/(y_interp(j1)-y_interp(j0))
         endif
         ! if (yd.lt.0_wp  ) yd = 0.0_wp
         ! if (yd.gt.1.0_wp) yd = 1.0_wp

         do i=1,nx
            ! loopx1: do i1=1,nil_interp
            !    if (x_interp(i1).gt.x(i)) exit loopx1
            ! enddo loopx1

            if (REAL(i/REAL(ndivx)).eq.REAL(i/ndivx)) then
               i0 = i/ndivx; i1 = i/ndivx
            else
               loopi_: do i_=1,ndivx
                  if (REAL((i-i_)/REAL(ndivx)).eq.REAL((i-i_)/ndivx)) exit loopi_
               enddo loopi_
               i0 = ((i-i_)/ndivx); i1 = i0+1
            endif

            if (i1.gt.nil_interp) then
               i1 = nil_interp+1
               i0 = nil_interp
             !   xd = (x(i)-x_interp(i0))/(x_interp(i1)-x_interp(i0))
             ! else
             !   i0 = i1-1
             !   xd = (x(i)-x_interp(i0))/(x_interp(i1)-x_interp(i0))
            endif
            ! if (xd.lt.0_wp  ) xd = 0.0_wp
            ! if (xd.gt.1.0_wp) xd = 1.0_wp

            ! Interpolation of extr_density
            ! =============================
            d00 = ((xc(i,j)-xc_interp(i0,j0))**2 + (yc(i,j)-yc_interp(i0,j0))**2)**0.5
            d01 = ((xc(i,j)-xc_interp(i0,j1))**2 + (yc(i,j)-yc_interp(i0,j1))**2)**0.5
            d10 = ((xc(i,j)-xc_interp(i1,j0))**2 + (yc(i,j)-yc_interp(i1,j0))**2)**0.5
            d11 = ((xc(i,j)-xc_interp(i1,j1))**2 + (yc(i,j)-yc_interp(i1,j1))**2)**0.5
            dtot = d00 + d01 + d10 + d11
            if (dtot.eq.0.0_wp) then
               dtot=4.0_wp; d00=1.0_wp; d01=1.0_wp; d10=1.0_wp; d11=1.0_wp
            endif
            d00=d00/dtot; d01=d01/dtot; d10=d10/dtot; d11=d11/dtot

            ! Ponderated values
            cc00 = extr_density(i0,j0,k0)*(1.0_wp-zd) + extr_density(i0,j0,k1)*zd
            cc01 = extr_density(i0,j1,k0)*(1.0_wp-zd) + extr_density(i0,j1,k1)*zd
            cc10 = extr_density(i1,j0,k0)*(1.0_wp-zd) + extr_density(i1,j0,k1)*zd
            cc11 = extr_density(i1,j1,k0)*(1.0_wp-zd) + extr_density(i1,j1,k1)*zd
            !
            cc  =  cc00*d00 + cc01*d01 + cc10*d10  + cc11*d11

            extr_density2(i,j,k) = cc
         enddo
      enddo
   enddo

  end subroutine fstt_interpol_3d_inv_c

end module mod_pp_fstt_interpol_c
