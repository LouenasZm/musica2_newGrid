!==============================================================================
module mod_pp_fstt_ltbl
!==============================================================================
  !> Module for discrimination between laminar and turbulent regions
!==============================================================================
  use warnstop
  use mod_constant
  use mod_flow
  use mod_mpi
  use mod_pp_var
  use mod_pp_fstt_comm
  implicit none
  ! ---------------------------------------------------------------------------

contains

  !============================================================================
  subroutine init_kernel
  !============================================================================
    !> Initialization of kernel density coefficients
    !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i,j,k,j1
    real(wp) :: step,step3
    ! -------------------------------------------------------------------------

    ! For interpolated grid
    step  = x_interp(2) - x_interp(1)
    step3 = z_interp(2) - z_interp(1)

    do j1=2,nj_interp
       do k=-nkernel_k,nkernel_k
          do j=born_kern_j(j1,1),born_kern_j(j1,2)
             do i=-ngh_pp,ngh_pp
                coeff_kernel(j1,i,j,k) = exp(-0.5_wp*((i*step)**2 + (y_interp(j)-y_interp(j1))**2 + (k*step3)**2)**0.5/(h_kern*step))
             enddo
          enddo
       enddo
    enddo

    ! For computationnal grid
    step  = x(2) - x(1)
    do j=2,ny
       born_kern_j2(j,1) = j
       do while ((exp(-0.5_wp*((y(born_kern_j2(j,1))-y(j))**2)**0.5/(h_kern*step))>prct_kernel).and.(born_kern_j2(j,1).gt.2))
          born_kern_j2(j,1) = born_kern_j2(j,1) - 1
       enddo
       born_kern_j2(j,2) = j
       do while ((exp(-0.5_wp*((y(born_kern_j2(j,2))-y(j))**2)**0.5/(h_kern*step))>prct_kernel).and.(born_kern_j2(j,2).lt.ny))
          born_kern_j2(j,2) = born_kern_j2(j,2) + 1
       enddo
    enddo

    do j1=2,ny
       do j=born_kern_j2(j1,1),born_kern_j2(j1,2)
          coeff_kernel2(j1,j) = exp(-0.5_wp*((y(j)-y(j1))**2)**0.5/(h_kern*step))
       enddo
    enddo

  end subroutine init_kernel


  !============================================================================
  subroutine init_kernel_c
  !============================================================================
    !> Initialization of kernel density coefficients
    !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i,j,k,i1,j1,ideb
    real(wp) :: step
    ! -------------------------------------------------------------------------

    ! For interpolated grid
    step = z_interp(2) - z_interp(1)
    ideb = nx_f_deb + coord(1)*nx ! initial position in the full file of coordinates
    do i1=1,nil_interp
       do j1=2,nj_interp
          do k=-nkernel_k,nkernel_k
             do j=bk_jc(i1,j1,1),bk_jc(i1,j1,2)
                do i=bk_ic(i1,j1,1)-i1,bk_ic(i1,j1,2)-i1
                   coeff_kc(i1,j1,i,j,k) = exp(-0.5_wp*((xgc_f(ndivx*(i+i1)+ideb,ndivy*j)-xgc_f(ndivx*i1+ideb,ndivy*j1))**2 + (k*step)**2 + &
                                                        (ygc_f(ndivx*(i+i1)+ideb,ndivy*j)-ygc_f(ndivx*i1+ideb,ndivy*j1))**2)**0.5/(h_kern*step))

                enddo
             enddo
          enddo
       enddo
    enddo

    ! For computationnal grid
    step  = z(2) - z(1)
    do i=1,nx
       do j=2,ny
          bk_jc2(i,j,1) = j
          do while ((exp(-0.5_wp*((y(bk_jc2(i,j,1))-y(j))**2)**0.5/(h_kern*step))>prct_kernel).and.(bk_jc2(i,j,1).gt.2))
             bk_jc2(i,j,1) = bk_jc2(i,j,1) - 1
          enddo
          bk_jc2(i,j,2) = j
          do while ((exp(-0.5_wp*((y(bk_jc2(i,j,2))-y(j))**2)**0.5/(h_kern*step))>prct_kernel).and.(bk_jc2(i,j,2).lt.ny))
             bk_jc2(i,j,2) = bk_jc2(i,j,2) + 1
          enddo
       enddo
    enddo

    do i1=1,nx
       do j1=2,ny
          do j=bk_jc2(i1,j1,1),bk_jc2(i1,j1,2)
             coeff_kc2(i1,j1,j) = exp(-0.5_wp*((yc(i1,j)-yc(i1,j1))**2)**0.5/(h_kern*step))
          enddo
       enddo
    enddo

    deallocate(xgc_f,ygc_f)

  end subroutine init_kernel_c

  !============================================================================
  subroutine fstt_bl_edge
  !============================================================================
    !> Calculation of the BL edge altitude
    !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i,j,k
    ! -------------------------------------------------------------------------

    ! Initialisation
    j_edge=0; j_edge2=0

    ! Edge indices for the interpolated grid
    ! --------------------------------------
    do k=1,nk_interp
       do i=1,nil_interp
          j=1
          do while ((uu_interp(i,j,k).lt.u_lim).and.(j.lt.nj_interp))
             j=j+1
          enddo
          j_edge(i,k) = j
       enddo
    enddo

    ! Values filled in the ghost points
    do i=1,nil_interp
       j_edge(i,1-nkernel_k:0) = j_edge(i,nk_interp+1-nkernel_k:nk_interp)
       j_edge(i,nk_interp+1:nk_interp+nkernel_k) = j_edge(i,1:nkernel_k)
    enddo

    ! Edge indices for the base grid
    ! ------------------------------
    do k=1,nz
       do i=1,nx
          j=1
          do while ((uu(i,j,k).lt.u_lim).and.(j.lt.ny))
             j=j+1
          enddo
          j_edge2(i,k) = j
       enddo
    enddo

  end subroutine fstt_bl_edge

  !============================================================================
  subroutine fstt_fluct_field
  !============================================================================
    !> Calculation of fluctuating fields
    !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i,j,k
    ! -------------------------------------------------------------------------

    do k=1,nk_interp
       do j=1,nj_interp
          do i=1,nil_interp
             ww_interp(i,j,k) = ww_interp(i,j,k) - stats_interp(i,j,4)
          enddo
       enddo
    enddo

    if (.not.is_curv) then
       do k=1,nk_interp
          do j=1,nj_interp
             do i=1,nil_interp
                ! uu_interp(i,j,k) = uu_interp(i,j,k) - stats_interp(i,j,2)
                vv_interp(i,j,k) = vv_interp(i,j,k) - stats_interp(i,j,3)
             enddo
          enddo
       enddo
    endif

    if (is_pp_streaks) then
       do k=1,nz
          do j=1,ny
             do i=1,nx
                uu_fluct(i,j,k) = uu(i,j,k) - stats_proc(i,j,2)
             enddo
          enddo
       enddo
    endif

  end subroutine fstt_fluct_field

  !============================================================================
  subroutine calc_uut_uun
  !============================================================================
    !> Calculation of tangential and normal velocities fluctuations
    !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i,j,k
    ! -------------------------------------------------------------------------

    ! tangential velocity fluctuations
    do k=1,nz
       do j=1,ny
          do i=1,nx
             uut(i,j,k) = (uu(i,j,k)-stats_proc(i,j,2))*nyn_jmin(i,k) - (vv(i,j,k)-stats_proc(i,j,3))*nxn_jmin(i,k)
          enddo
       enddo
    enddo

    ! normal velocity fluctuations
    do k=1,nz
       do j=1,ny
          do i=1,nx
             uun(i,j,k) = (uu(i,j,k)-stats_proc(i,j,2))*nxn_jmin(i,k) + (vv(i,j,k)-stats_proc(i,j,3))*nyn_jmin(i,k)
          enddo
       enddo
    enddo

  end subroutine calc_uut_uun




  !============================================================================
  subroutine fstt_lam_turb
  !============================================================================
    !> Discrimination between laminar and turbulent regions
    !>   -> Computation of extremums density distribution over the
    !>      interpolated grid
    !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i,j,k,i1,j1,k1,nextr
    real(wp) :: temp,temp2,temp3,temp4
    real(wp), dimension(1-ngh_pp:nil_interp+ngh_pp,1:nj_interp,1-nkernel_k:nk_interp+nkernel_k) :: extr
    ! -------------------------------------------------------------------------

    nextr = 0
    extr_density = 0.0_wp
    extr = 0.0_wp

    ! Determination of the extremums in the field
    ! -------------------------------------------
    if (is_curv) then
       do k=1,nk_interp
          do j=2,nj_interp
             do i=ndx_pp,nfx_pp
                temp  = maxval(uun_interp(i-1:i+1,j:j,k-1:k+1))
                temp2 = minval(uun_interp(i-1:i+1,j:j,k-1:k+1))
                temp3 = maxval(ww_interp(i-1:i+1,j:j,k-1:k+1))
                temp4 = minval(ww_interp(i-1:i+1,j:j,k-1:k+1))
                if (uun_interp(i,j,k).eq.temp) then
                   extr(i,j,k) = abs(temp)
                   ! extr(i,j,k) = 1.0_wp
                else if (uun_interp(i,j,k).eq.temp2) then
                   extr(i,j,k) = abs(temp2)
                   ! extr(i,j,k) = 1.0_wp
                endif
                if (ww_interp(i,j,k).eq.temp3) then
                   extr(i,j,k) = extr(i,j,k) + abs(temp3)
                else if (ww_interp(i,j,k).eq.temp4) then
                   extr(i,j,k) = extr(i,j,k) + abs(temp4)
                endif
             enddo
          enddo
       enddo
    else
       ! Extremums positions stored in extr
       if (.not.is_lim_fst) then
          do k=1,nk_interp
             do j=2,nj_interp
                do i=ndx_pp,nfx_pp
                   temp  = maxval(vv_interp(i-1:i+1,j:j,k-1:k+1))
                   temp2 = minval(vv_interp(i-1:i+1,j:j,k-1:k+1))
                   temp3 = maxval(ww_interp(i-1:i+1,j:j,k-1:k+1))
                   temp4 = minval(ww_interp(i-1:i+1,j:j,k-1:k+1))
                   if (vv_interp(i,j,k).eq.temp) then
                      extr(i,j,k) = abs(temp)
                   else if (vv_interp(i,j,k).eq.temp2) then
                      extr(i,j,k) = abs(temp2)
                   endif
                   if (ww_interp(i,j,k).eq.temp3) then
                      extr(i,j,k) = extr(i,j,k) + abs(temp3)
                   else if (ww_interp(i,j,k).eq.temp4) then
                      extr(i,j,k) = extr(i,j,k) + abs(temp4)
                   endif
                enddo
             enddo
          enddo
       else
          do k=1,nk_interp
             do i=ndx_pp,nfx_pp
                ! do j=2,j_1p5d99(i)
                do j=2,j_edge(i,k)
                   temp  = maxval(vv_interp(i-1:i+1,j:j,k-1:k+1))
                   temp2 = minval(vv_interp(i-1:i+1,j:j,k-1:k+1))
                   temp3 = maxval(ww_interp(i-1:i+1,j:j,k-1:k+1))
                   temp4 = minval(ww_interp(i-1:i+1,j:j,k-1:k+1))
                   if (vv_interp(i,j,k).eq.temp) then
                      extr(i,j,k) = abs(temp)
                   else if (vv_interp(i,j,k).eq.temp2) then
                      extr(i,j,k) = abs(temp2)
                   endif
                   if (ww_interp(i,j,k).eq.temp3) then
                      extr(i,j,k) = extr(i,j,k) + abs(temp3)
                   else if (ww_interp(i,j,k).eq.temp4) then
                      extr(i,j,k) = extr(i,j,k) + abs(temp4)
                   endif
                enddo
             enddo
          enddo
       endif
    endif

    ! Communication of extremums
    ! --------------------------
    call communication_fstt1(extr)


    ! Application of the kernel over the extremums
    ! --------------------------------------------
    if (is_curv) then
       do k=1,nk_interp
          do j=2,nj_interp
             do i=ndxtpngh_pp,nfxtmngh_pp
                do k1=-nkernel_k,nkernel_k
                do j1=bk_jc(i,j,1),bk_jc(i,j,2)
                do i1=bk_ic(i,j,1)-i,bk_ic(i,j,2)-i
                extr_density(i,j,k) = extr_density(i,j,k) + extr(i+i1,j1,k+k1)*coeff_kc(i,j,i1,j1,k1)
                enddo
                enddo
                enddo
             enddo
          enddo
       enddo
    else
       do k=1,nk_interp
          do j=2,nj_interp
             do i=ndxtpngh_pp,nfxtmngh_pp
                do k1=-nkernel_k,nkernel_k
                do j1=born_kern_j(j,1),born_kern_j(j,2)
                do i1=-ngh_pp,ngh_pp
                extr_density(i,j,k) = extr_density(i,j,k) + extr(i+i1,j1,k+k1)*coeff_kernel(j,i1,j1,k1)
                enddo
                enddo
                enddo
             enddo
          enddo
       enddo
    endif

    ! For imin and imax
    if (is_boundary(1,1)) then
       do k=1,nk_interp
          do j=2,nj_interp
             do i=ndxt_pp,ndxtpngh_pp-1
                extr_density(i,j,k) = extr_density(ndxtpngh_pp,j,k)
             enddo
          enddo
       enddo
    endif
    if (is_boundary(1,2)) then
       do k=1,nk_interp
          do j=2,nj_interp
             do i=nfxtmngh_pp+1,nfxt_pp
                extr_density(i,j,k) = extr_density(nfxtmngh_pp,j,k)
             enddo
          enddo
       enddo
    endif

    ! For jmin
    do k=1,nk_interp
       do i=1,nil_interp
          extr_density(i,1,k) = extr_density(i,2,k)
       enddo
    enddo


    ! Communication for k=0 and k=nk_interp+1 and i=0 and i=nil_interp+1
    call communication_fstt1(extr_density)

  end subroutine fstt_lam_turb


  !============================================================================
  subroutine fstt_thresold
  !============================================================================
    !> Discrimination between laminar and turbulent regions
    !>   -> Application of thresold (Otsu, 1979)
    !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i,j,k,i1,j1,N_pop,lvl_,lvl_min,lvl_max
    real(wp) :: temp,max_D_loc
    integer, dimension(1) :: lvl_opt
    real(wp), dimension(lvl_th) :: p_lvl_loc,p_lvl
    real(wp), dimension(lvl_th-1) :: omega0,omega1,mu0,mu1,sigma0,sigma1,etat
    real(wp), dimension(2:ny) :: thres_2d,thres_2d_smooth
    real(wp), dimension(1:nx,2:ny) :: thres_2d_smc
    ! -------------------------------------------------------------------------

    ! Initialisation of array
    ltbl = 0.0_wp

    ! Application of thresold (Otsu 1979)
    ! Application in 2D
    ! N_pop = ngx*ngz
    N_pop = 0
    do i=1,nbloc
       N_pop = N_pop + bl(i)%ni*bl(i)%nk
    enddo

    loopj: do j=2,ny
       ! Initialisation
       p_lvl_loc = 0
       mu0 = 0.0_wp; mu1 = 0.0_wp
       sigma0 = 0.0_wp; sigma1 = 0.0_wp
       etat = 0.0_wp

       ! Max value for normalisation
       max_D_loc = maxval(extr_density2(:,j,:))
       call MPI_ALLREDUCE(max_D_loc,max_D,1,MPI_DOUBLE_PRECISION,MPI_MAX,COMM_global,info)
       max_D = max_D*1.00000001_wp


       ! if (iproc.eq.0) print *,"max_D",max_D

       if (max_D.eq.0.0_wp) cycle loopj

       ! Count of population for each discrete step
       do k=1,nz
          do i=1,nx
          ! do i=ndxtpngh_pp,nfxtmngh_pp
             i1 = int(extr_density2(i,j,k)*lvl_th/max_D) + 1
             p_lvl_loc(i1) = p_lvl_loc(i1) + 1
          enddo
       enddo
       call MPI_ALLREDUCE(p_lvl_loc,p_lvl,lvl_th,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)


       ! if (iproc.eq.0) print *,"p_lvl_loc",p_lvl_loc
       ! if (iproc.eq.0) print *,"p_lvl",p_lvl

       ! Proportion of population at each discrete step
       do lvl_=1,lvl_th
          p_lvl(lvl_) = p_lvl(lvl_)/N_pop
       enddo

       ! if (iproc.eq.0) print *,"p_lvl prop",p_lvl
       ! if (iproc.eq.0) print *,"SUM(p_lvl prop)",sum(p_lvl)

       ! Calculation of classes occurence, mean levels, variance, ...
       do lvl_=1,lvl_th-1
          omega0(lvl_) = sum(p_lvl(1:lvl_))
          omega1(lvl_) = sum(p_lvl(lvl_+1:lvl_th))
          ! if (iproc.eq.0) print *,"lvl",lvl_,omega0(lvl_),omega1(lvl_)
       enddo

       ! Restriction for optimal level: 0 < omega0 < 1
       ! Minimum level
       lvl_ = 1
       do while (omega0(lvl_).le.0.0_wp)
          lvl_ = lvl_ + 1
       enddo
       lvl_min = lvl_
       ! Maximum level
       lvl_ = lvl_th-1
       do while (omega1(lvl_).le.0.0_wp)
          lvl_ = lvl_ - 1
       enddo
       lvl_max = lvl_

       ! if (iproc.eq.0) print *," lvl_min lvl_max",lvl_min,lvl_max

       do lvl_=lvl_min,lvl_max
          do i=1,lvl_
             mu0(lvl_) = mu0(lvl_) + i*p_lvl(i)
          enddo
          mu0(lvl_) = mu0(lvl_)/omega0(lvl_)

          do k=lvl_+1,lvl_th
             mu1(lvl_) = mu1(lvl_) + k*p_lvl(k)
          enddo
          mu1(lvl_) = mu1(lvl_)/omega1(lvl_)
       enddo

       do lvl_=lvl_min,lvl_max
          do i=1,lvl_
             sigma0(lvl_) = sigma0(lvl_) + (i-mu0(lvl_))**2 * p_lvl(i)
          enddo
          sigma0(lvl_) = sigma0(lvl_)/omega0(lvl_)

          do k=lvl_+1,lvl_th
             sigma1(lvl_) = sigma1(lvl_) + (k-mu1(lvl_))**2 * p_lvl(k)
          enddo
          sigma1(lvl_) = sigma1(lvl_)/omega1(lvl_)
       enddo

       ! If all population is localized at extremities
       if (p_lvl(1)+p_lvl(lvl_th).eq.1.0_wp) then
          lvl_opt = lvl_th/2
       else
          do lvl_=lvl_min,lvl_max
             etat(lvl_) = (omega0(lvl_)*omega1(lvl_)*(mu1(lvl_) - mu0(lvl_))**2)/(omega0(lvl_)*sigma0(lvl_) + omega1(lvl_)*sigma1(lvl_))
          enddo
          ! Determination of thresold via maximisation of the inter-class variance
          lvl_opt = maxloc(etat(lvl_min:lvl_max))

          thres_2d(j) = max_D*lvl_opt(1)/lvl_th
       endif


       ! if (iproc.eq.0) print *,"thres_2d",j,thres_2d(j)

    enddo loopj


    if (is_curv) then
       thres_2d_smc = 0.0_wp
       do i=1,nx
          do j=2,ny
            ! Weighted average
             temp = 0
             do j1=bk_jc2(i,j,1),bk_jc2(i,j,2)
                thres_2d_smc(i,j) = thres_2d_smc(i,j) + coeff_kc2(i,j,j1)*thres_2d(j1)
                temp = temp + coeff_kc2(i,j,j1)
             enddo
             thres_2d_smc(i,j) = thres_2d_smc(i,j)/temp

             ! Discrimination between laminar and turbulent BL based on thresold
             do k=1,nz
                if (extr_density2(i,j,k).lt.thres_2d_smc(i,j)) then
                   ltbl(i,j,k) = 0.0_wp
                else
                   ltbl(i,j,k) = 1.0_wp
                endif
             enddo
          enddo
       enddo
    else
       thres_2d_smooth = 0.0_wp
       do j=2,ny
          ! Weighted average
          temp = 0
          do j1=born_kern_j2(j,1),born_kern_j2(j,2)
             thres_2d_smooth(j) = thres_2d_smooth(j) + coeff_kernel2(j,j1)*thres_2d(j1)
             temp = temp + coeff_kernel2(j,j1)
          enddo
          thres_2d_smooth(j) = thres_2d_smooth(j)/temp

          ! Discrimination between laminar and turbulent BL based on thresold
          do k=1,nz
             do i=1,nx
                if (extr_density2(i,j,k).lt.thres_2d_smooth(j)) then
                   ltbl(i,j,k) = 0.0_wp
                else
                   ltbl(i,j,k) = 1.0_wp
                endif
             enddo
          enddo
       enddo
    endif

    ! Wall values set to point above
    ! ------------------------------
    do k=1,nz
       do i=1,nx
          ltbl(i,1,k) = ltbl(i,2,k)
       enddo
    enddo

  end subroutine fstt_thresold



  ! !==============================================================================
  ! subroutine filtre_lam_turb(j)
  ! !==============================================================================
  !   !> Filtering of ltbl
  !   !> - Cartesian version -
  !   !==============================================================================
  !   use mod_coeff_deriv
  !   use mod_interface
  !   use mod_flow
  !   implicit none
  !   ! -------------------------------------------------------------------------
  !   integer, intent(in)  :: j
  !   integer :: i,k
  !   real :: d11t(0:5)
  !   real, dimension(nil_interp,nk_interp) :: ltbl_f
  !   ! -------------------------------------------------------------------------



  !   ! Initialisation
  !   ! ---------------

  !   ltbl_f = 0.

  !   ! ! Coupure en pi/3
  !   ! ! ---------------
  !   ! d11t(0)= 2./3.
  !   ! d11t(1)=-0.26775782
  !   ! d11t(2)=-0.12016956
  !   ! d11t(3)= 0.
  !   ! d11t(4)= 0.03683622
  !   ! d11t(5)= 0.01775782

  !   ! ! Standard filter order 10 (11pts)
  !   ! ! ------------------------
  !   ! d11t(0)= 252.0_wp/1024.0_wp
  !   ! d11t(1)=-210.0_wp/1024.0_wp
  !   ! d11t(2)= 120.0_wp/1024.0_wp
  !   ! d11t(3)=- 45.0_wp/1024.0_wp
  !   ! d11t(4)=  10.0_wp/1024.0_wp
  !   ! d11t(5)=-  1.0_wp/1024.0_wp

  !   ! Filtrage en x sur 11 points
  !   ! ---------------------------
  !   do k=1,nk_interp
  !      do i=1,nil_interp
  !         ltbl_f(i,k)= ( d11t(0) * ltbl(i,j,k)        &
  !                      + d11t(1) * ( ltbl(i+1,j,k)+ltbl(i-1,j,k) ) &
  !                      + d11t(2) * ( ltbl(i+2,j,k)+ltbl(i-2,j,k) ) &
  !                      + d11t(3) * ( ltbl(i+3,j,k)+ltbl(i-3,j,k) ) &
  !                      + d11t(4) * ( ltbl(i+4,j,k)+ltbl(i-4,j,k) ) &
  !                      + d11t(5) * ( ltbl(i+5,j,k)+ltbl(i-5,j,k) ) )
  !      enddo
  !   enddo


  !   ! Filtrage en z sur 11 points
  !   ! ---------------------------
  !   do k=1,nk_interp
  !      do i=1,nil_interp
  !         ltbl_f(i,k) =  ltbl_f(i,k) + ( d11t(0) * ltbl(i,j,k)  &
  !                        + d11t(1) * ( ltbl(i,j,k+1)+ltbl(i,j,k-1) ) &
  !                        + d11t(2) * ( ltbl(i,j,k+2)+ltbl(i,j,k-2) )     &
  !                        + d11t(3) * ( ltbl(i,j,k+3)+ltbl(i,j,k-3) )     &
  !                        + d11t(4) * ( ltbl(i,j,k+4)+ltbl(i,j,k-4) )     &
  !                        + d11t(5) * ( ltbl(i,j,k+5)+ltbl(i,j,k-5) ) )
  !      enddo
  !   enddo

  !   ! Moyenne pondérée
  !   ! ----------------
  !   d11t(0)=  1./6.*0.5
  !   d11t(1)= 1./12.*0.5
  !   d11t(2)= 1./12.*0.5
  !   d11t(3)= 1./12.*0.5
  !   d11t(4)= 1./12.*0.5
  !   d11t(5)= 1./12.*0.5



  !   ! Filtrage des ltbl
  !   ! ----------------

  !   do k=1,nk_interp
  !      do i=1,nil_interp
  !         ! if (iproc==27) print *,'ltbl',ltbl(i,j,k)
  !         ! if (iproc==27) print *,'ltbl_f',ltbl_f(i,k)
  !         ! if (iproc==27) print *,'2*ltbl_f',2*ltbl_f(i,k)
  !         ! ltbl(i,j,k) = ltbl(i,j,k) - int(2*ltbl_f(i,k))
  !         ! if (iproc==27) print *,'ltbl before',ltbl(i,j,k)
  !         ! if (iproc==27) print *,'ltbl_f',ltbl_f(i,k)
  !         ltbl(i,j,k) = 0.5_wp*ltbl_f(i,k)
  !         ! ltbl(i,j,k) = ltbl(i,j,k) - ltbl_f(i,k)
  !         ! if (iproc==27) print *,'ltbl after',ltbl(i,j,k)
  !         ! if (iproc==27) print *,'la',ltbl(i,j,k)
  !         ! if (iproc==27) print *,'la',int(2*(ltbl(i,j,k) + ltbl_f(i,k)))/2
  !         ! ltbl(i,j,k) = ltbl(i,j,k) + ltbl_f(i,k)
  !      enddo
  !   enddo

  !   ! call mpistop('',0)

  ! end subroutine filtre_lam_turb

  !============================================================================
  subroutine fstt_stats_ltbl
  !============================================================================
    !> Calculation of stats in laminar and turbulent regions separated
    !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer  :: i,j,k,l
    real(wp) :: nn1,inn1
    ! -------------------------------------------------------------------------

    nn1 = dble(nout + nbeg_pp - 1)
    inn1= 1.0_wp/nn1

    ! Spatial local and global averages at this iteration
    avg_s_lam  = 0.0_wp; avg_s_turb  = 0.0_wp

    do k=1,nz
       do j=1,ny
          do i=1,nx
             if (ltbl(i,j,k).lt.0.5_wp) then
                avg_s_lam(i,j, 1)= avg_s_lam(i,j, 1) + 1
                avg_s_lam(i,j, 2)= avg_s_lam(i,j, 2) + uu(i,j,k)
                avg_s_lam(i,j, 3)= avg_s_lam(i,j, 3) + vv(i,j,k)
                avg_s_lam(i,j, 4)= avg_s_lam(i,j, 4) + ww(i,j,k)
                avg_s_lam(i,j, 5)= avg_s_lam(i,j, 5) + uu(i,j,k)**2
                avg_s_lam(i,j, 6)= avg_s_lam(i,j, 6) + vv(i,j,k)**2
                avg_s_lam(i,j, 7)= avg_s_lam(i,j, 7) + ww(i,j,k)**2
                avg_s_lam(i,j, 8)= avg_s_lam(i,j, 8) + uu(i,j,k)*vv(i,j,k)
                avg_s_lam(i,j, 9)= avg_s_lam(i,j, 9) + uu(i,j,k)*ww(i,j,k)
                avg_s_lam(i,j,10)= avg_s_lam(i,j,10) + vv(i,j,k)*ww(i,j,k)

                ! avg_s_lam(i,j,)= avg_s_lam(i,j,) + rho_i
                ! avg_s_lam(i,j,)= avg_s_lam(i,j,) + p_i
                ! avg_s_lam(i,j,)= avg_s_lam(i,j,) + T_i
                ! avg_s_lam(i,j,)= avg_s_lam(i,j,) + rhou(i,j,k)
                ! avg_s_lam(i,j,)= avg_s_lam(i,j,) + rhov(i,j,k)
                ! avg_s_lam(i,j,)= avg_s_lam(i,j,) + rhow(i,j,k)
                ! avg_s_lam(i,j,)= avg_s_lam(i,j,) + rhoe(i,j,k)
                ! avg_s_lam(i,j,)= avg_s_lam(i,j,) + rho_i**2
                ! avg_s_lam(i,j,)= avg_s_lam(i,j,) + v_i*T_i
                ! avg_s_lam(i,j,)= avg_s_lam(i,j,) + p_i**2
                ! avg_s_lam(i,j,)= avg_s_lam(i,j,) + T_i**2
                ! avg_s_lam(i,j,)= avg_s_lam(i,j,) + mu
                ! avg_s_lam(i,j,)= avg_s_lam(i,j,) + divloc
                ! avg_s_lam(i,j,)= avg_s_lam(i,j,) + divloc**2
             else
                avg_s_turb(i,j, 1)= avg_s_turb(i,j, 1) + 1
                avg_s_turb(i,j, 2)= avg_s_turb(i,j, 2) + uu(i,j,k)
                avg_s_turb(i,j, 3)= avg_s_turb(i,j, 3) + vv(i,j,k)
                avg_s_turb(i,j, 4)= avg_s_turb(i,j, 4) + ww(i,j,k)
                avg_s_turb(i,j, 5)= avg_s_turb(i,j, 5) + uu(i,j,k)**2
                avg_s_turb(i,j, 6)= avg_s_turb(i,j, 6) + vv(i,j,k)**2
                avg_s_turb(i,j, 7)= avg_s_turb(i,j, 7) + ww(i,j,k)**2
                avg_s_turb(i,j, 8)= avg_s_turb(i,j, 8) + uu(i,j,k)*vv(i,j,k)
                avg_s_turb(i,j, 9)= avg_s_turb(i,j, 9) + uu(i,j,k)*ww(i,j,k)
                avg_s_turb(i,j,10)= avg_s_turb(i,j,10) + vv(i,j,k)*ww(i,j,k)
             endif
          enddo
       enddo
    enddo

    ! Spanwise averaging
    do l=2,10
       do j=1,ny
          do i=1,nx
             if (avg_s_lam(i,j,1).gt.0) avg_s_lam(i,j,l)  = avg_s_lam(i,j,l) / avg_s_lam(i,j,1)
             if (avg_s_turb(i,j,1).gt.0) avg_s_turb(i,j,l) = avg_s_turb(i,j,l)/avg_s_turb(i,j,1)
          enddo
       enddo
    enddo

    ! Spanwise averaging of the intermittency
    avg_s_lam(:,:,1)  = avg_s_lam(:,:,1)/dble(nz)
    avg_s_turb(:,:,1) = avg_s_turb(:,:,1)/dble(nz)

    ! Count contribution for this iteration at each position
    do j=1,ny
       do i=1,nx
          stats_cpt(1,i,j) = stats_cpt(1,i,j) + avg_s_lam(i,j,1)
          stats_cpt(2,i,j) = stats_cpt(2,i,j) + avg_s_turb(i,j,1)
       enddo
    enddo

    ! Temporal averaging
    do l=2,10
       do j=1,ny
          do i=1,nx
             if (stats_cpt(1,i,j).gt.0.0_wp) stats_lam(i,j,l)  = ((stats_cpt(1,i,j) -avg_s_lam(i,j,1)) *stats_lam(i,j,l)  + avg_s_lam(i,j,1)*avg_s_lam(i,j,l))/stats_cpt(1,i,j)
             if (stats_cpt(2,i,j).gt.0.0_wp) stats_turb(i,j,l)  = ((stats_cpt(2,i,j)-avg_s_turb(i,j,1))*stats_turb(i,j,l) + avg_s_turb(i,j,1)*avg_s_turb(i,j,l))/stats_cpt(2,i,j)
          enddo
       enddo
    enddo


    ! Temporal averaging for the intermittency
    l=1
    do j=1,ny
       do i=1,nx
          stats_lam(i,j,l)  = ((nn1-1.0_wp)*stats_lam(i,j,l)  + avg_s_lam(i,j,l))*inn1
          stats_turb(i,j,l) = ((nn1-1.0_wp)*stats_turb(i,j,l) + avg_s_turb(i,j,l))*inn1
       enddo
    enddo

    if (is_check_stats) then
       avg_s_tot  = 0.0_wp

       do k=1,nz
          do j=1,ny
             do i=1,nx
                ! ! Temporary
                avg_s_tot(i,j, 1)= avg_s_tot(i,j, 1) + 1
                avg_s_tot(i,j, 2)= avg_s_tot(i,j, 2) + uu(i,j,k)
                avg_s_tot(i,j, 3)= avg_s_tot(i,j, 3) + vv(i,j,k)
                avg_s_tot(i,j, 4)= avg_s_tot(i,j, 4) + ww(i,j,k)
                avg_s_tot(i,j, 5)= avg_s_tot(i,j, 5) + uu(i,j,k)**2
                avg_s_tot(i,j, 6)= avg_s_tot(i,j, 6) + vv(i,j,k)**2
                avg_s_tot(i,j, 7)= avg_s_tot(i,j, 7) + ww(i,j,k)**2
                avg_s_tot(i,j, 8)= avg_s_tot(i,j, 8) + uu(i,j,k)*vv(i,j,k)
                avg_s_tot(i,j, 9)= avg_s_tot(i,j, 9) + uu(i,j,k)*ww(i,j,k)
                avg_s_tot(i,j,10)= avg_s_tot(i,j,10) + vv(i,j,k)*ww(i,j,k)
             enddo
          enddo
       enddo

       avg_s_tot = avg_s_tot/dble(nz)

       stats_tot = ((nn1-1.0_wp)*stats_tot + avg_s_tot)*inn1

    endif

  end subroutine fstt_stats_ltbl

end module mod_pp_fstt_ltbl
