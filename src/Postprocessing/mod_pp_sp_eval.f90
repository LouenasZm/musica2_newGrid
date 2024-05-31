!==============================================================================
module mod_pp_sp_eval
!==============================================================================
  !> Module for spectrum evaluation and writing
!==============================================================================
  use mod_pp_mpi
  use mod_pp_var
  implicit none
  ! ---------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------

contains

  !============================================================================
  subroutine compute_autospectrum
  !============================================================================
    !> Compute autospectra (1-D)
    !>  * use real arrays *
  !============================================================================
    use mod_pp_dfft
    use mod_pp_capon
    use mod_pp_transpose
    use mod_pp_sp_wind
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i1,i2,n
    real(wp) :: eps ! MVSE regularization
    ! -------------------------------------------------------------------------
    
    ! small positive number for Tikhonov's regularization
    ! to prevent ill-conditioning of autocorrelation matrix
    eps=1.0e-12_wp 

    if (iproc==0) then
       print *,'=========================================================='
       print *,' autospectrum evaluation'
       print *,'=========================================================='
       print *,'~> spectrum in '//sp%d(1)%name//'-direction ...'
    endif

    select case (sp%d(1)%name)
    case('x')          
       if (ndomx>1) then ! (parallel implementation)
          ! Transposition
          ! -------------
          if (iproc==0) print *,'   - xt-transpose'
          call transpose_xt_r

          ! Windowing
          ! ---------
          call windowing_r(var_rt,1,3)

          if (sp%is_capon) then
             ! Capon's transform in x
             ! ----------------------
             if (iproc==0) print *,'   - Capon''s estimator (MVSE)'
             do i2=1,nl2
                if ((iproc==0).and.(mod(i2,10)==0)) print *,'     ind',i2
                do n=1,nt_txt
                   call fast_mvse_r(var_rt(n,i2,:),sp%d(1)%lbloc,sp%ncapon,eps)
                enddo
             enddo
             ! normalization of FFT
             ! --------------------
             var_rt=var_rt/sp%d(1)%lbloc 
          else
             ! DFFT in x (third direction after transposition)
             ! ---------
             if (iproc==0) print *,'   - DFFT'
             call DFFT(var_rt,3)
             ! normalization of FFT
             ! --------------------
             var_rt=var_rt/sp%d(1)%lbloc 
             ! compute module and phase
             ! ------------------------
             call sep_mod_phase(var_rt,3)
             var_rt=var_rt**2
          endif

          ! Reverse transposition
          ! ---------------------
          if (iproc==0) print *,'   - xt-untranspose'
          call untranspose_xt_r

       else ! (serial implementation)
          ! Windowing
          ! ---------
          call windowing_r(var_r,1,1)

          if (sp%is_capon) then
             ! Capon's transform in x
             ! ----------------------
             if (iproc==0) print *,'   - Capon''s estimator (MVSE)'
             do i2=1,nl2
                if ((iproc==0).and.(mod(i2,10)==0)) print *,'     ind',i2
                do n=1,ngt
                   call fast_mvse_r(var_r(:,i2,n),sp%d(1)%lbloc,sp%ncapon,eps)
                enddo
             enddo
             ! normalization of FFT
             ! --------------------
             var_r=var_r/sp%d(1)%lbloc 
          else
             ! DFFT in x
             ! ---------
             if (iproc==0) print *,'   - DFFT'
             call DFFT(var_r,sp%d(1)%i)
             ! normalization of FFT
             ! --------------------
             var_r=var_r/sp%d(1)%lbloc 
             ! compute module and phase
             ! ------------------------
             call sep_mod_phase(var_r,sp%d(1)%i)
             var_r=var_r**2
          endif

       endif

    case('z')          
       if (ndomz>1) then ! (parallel implementation)
          ! Transposition
          ! -------------
          if (iproc==0) print *,'   - zt-transpose'
          call transpose_zt_r

          ! Windowing
          ! ---------
          call windowing_r(var_rt,1,3)

          if (sp%is_capon) then
             ! Capon's transform in z
             ! ----------------------
             if (iproc==0) print *,'   - Capon''s estimator (MVSE)'
             do i1=1,nl1
                if ((iproc==0).and.(mod(i1,10)==0)) print *,'     ind',i1
                do n=1,nt_tzt
                   call fast_mvse_r(var_rt(i1,n,:),sp%d(1)%lbloc,sp%ncapon,eps)
                enddo
             enddo
             ! normalization of FFT
             ! --------------------
             var_rt=var_rt/sp%d(1)%lbloc 
          else
             ! DFFT in z (third direction after transposition)
             ! ---------
             if (iproc==0) print *,'   - DFFT'
             call DFFT(var_rt,3)
             ! normalization of FFT
             ! --------------------
             var_rt=var_rt/sp%d(1)%lbloc 
             ! compute module and phase
             ! ------------------------
             call sep_mod_phase(var_rt,3)
             var_rt=var_rt**2
          endif

          ! Reverse transposition
          ! ---------------------
          if (iproc==0) print *,'   - zt-untranspose'
          call untranspose_zt_r

       else ! (serial implementation)
          ! Windowing
          ! ---------
          call windowing_r(var_r,1,2)

          if (sp%is_capon) then
             ! Capon's transform in z
             ! ----------------------
             if (iproc==0) print *,'   - Capon''s estimator (MVSE)'
             do i1=1,nl1
                if ((iproc==0).and.(mod(i1,10)==0)) print *,'     ind',i1
                do n=1,ngt
                   call fast_mvse_r(var_r(i1,:,n),sp%d(1)%lbloc,sp%ncapon,eps)
                enddo
             enddo
             ! normalization of FFT
             ! --------------------
             var_r=var_r/sp%d(1)%lbloc 
          else
             ! DFFT in z
             ! ---------
             if (iproc==0) print *,'   - DFFT'
             call DFFT(var_r,sp%d(1)%i)
             ! normalization of FFT
             ! --------------------
             var_r=var_r/sp%d(1)%lbloc 
             ! compute module and phase
             ! ------------------------
             call sep_mod_phase(var_r,sp%d(1)%i)
             var_r=var_r**2
          endif

       endif

    case('t')          
       ! Windowing
       ! ---------
       call windowing_r(var_r,1,3)
       
       if (sp%is_capon) then
          ! Capon's transform in t
          ! ----------------------
          if (iproc==0) print *,'   - Capon''s estimator (MVSE)'
          do i2=1,nl2
             if ((iproc==0).and.(mod(i2,10)==0)) print *,'     ind',i2
             do i1=1,nl1
                call fast_mvse_r(var_r(i1,i2,:),sp%d(1)%lbloc,sp%ncapon,eps)
             enddo
          enddo
          ! normalization of FFT
          ! --------------------
          var_r=var_r/sp%d(1)%lbloc 
       else
          ! DFFT in t
          ! ---------
          if (iproc==0) print *,'   - DFFT'
          call DFFT(var_r,3)
          ! normalization of FFT
          ! --------------------
          var_r=var_r/sp%d(1)%lbloc 
          ! compute module and phase
          ! ------------------------
          call sep_mod_phase(var_r,sp%d(1)%i)
          var_r=var_r**2
       endif

    end select
    
  end subroutine compute_autospectrum

  !============================================================================
  subroutine compute_multiD_spectrum
  !============================================================================
    !> Compute multidimensional spectra (2-D or 3-D)
    !>  * use complex arrays *
  !============================================================================
    use mod_pp_dfft
    use mod_pp_capon
    use mod_pp_transpose
    use mod_pp_sp_wind
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i1,i2,n,nd    
    real(wp) :: eps ! MVSE regularization
    ! -------------------------------------------------------------------------
    
    ! small positive number for Tikhonov's regularization
    ! to prevent ill-conditioning of autocorrelation matrix
    eps=1.0e-12_wp 

    if (iproc==0) then
       if (sp%dim==2) then
          print *,'=========================================================='
          print *,'two-dimensional spectrum evaluation'
          print *,'=========================================================='
       elseif (sp%dim==3) then
          print *,'=========================================================='
          print *,'frequency-wavenumber spectrum evaluation'
          print *,'=========================================================='
       endif
    endif

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! loop on spectrum directions
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%
    do nd=1,sp%dim
       
       if (iproc==0) print *,'~> spectrum in '//sp%d(nd)%name//'-direction ...'
       
       select case (sp%d(nd)%name)
       case('x')          
          if (ndomx>1) then ! (parallel implementation)
             ! Transposition
             ! -------------
             if (iproc==0) print *,'   - xt-transpose'
             call transpose_xt_c
             
             ! Windowing
             ! ---------
             call windowing_c(var_ct,nd,3)

             if ((nd==sp%dim).and.(sp%is_capon)) then
                ! Capon's transform in x
                ! ----------------------
                if (iproc==0) print *,'   - Capon''s estimator (MVSE)'
                do i2=1,nl2
                   if ((iproc==0).and.(mod(i2,10)==0)) print *,'     ind',i2
                   do n=1,nt_txt
                      call fast_mvse_c(var_ct(n,i2,:),sp%d(nd)%lbloc,sp%ncapon,eps)
                   enddo
                enddo
             else
                ! DFFT in x (third direction after transposition)
                ! ---------
                if (iproc==0) print *,'   - ZFFT+shift'
                call ZFFT_shift(var_ct,3)
                var_ct=var_ct/sp%d(nd)%lbloc ! normalization of FFT
             endif

             ! Reverse transposition
             ! ---------------------
             if (iproc==0) print *,'   - xt-untranspose'
             call untranspose_xt_c

          else ! (serial implementation)
             ! Windowing
             ! ---------
             call windowing_c(var_c,nd,nd)

             if ((nd==sp%dim).and.(sp%is_capon)) then
                ! Capon's transform in x
                ! ----------------------
                if (iproc==0) print *,'   - Capon''s estimator (MVSE)'
                do i2=1,nl2
                   if ((iproc==0).and.(mod(i2,10)==0)) print *,'     ind',i2
                   do n=1,ngt
                      call fast_mvse_c(var_c(:,i2,n),sp%d(nd)%lbloc,sp%ncapon,eps)
                   enddo
                enddo
            else
                ! DFFT in x
                ! ---------
                if (iproc==0) print *,'   - ZFFT+shift'
                call ZFFT_shift(var_c,sp%d(nd)%i)
                var_c=var_c/sp%d(nd)%lbloc ! normalization of FFT
             endif

          endif

       case('z')          
          if (ndomz>1) then ! (parallel implementation)
             ! Transposition
             ! -------------
             if (iproc==0) print *,'   - zt-transpose'
             call transpose_zt_c
             
             ! Windowing
             ! ---------
             call windowing_c(var_ct,1,3)

             if ((nd==sp%dim).and.(sp%is_capon)) then
                ! Capon's transform in z
                ! ----------------------
                if (iproc==0) print *,'   - Capon''s estimator (MVSE)'
                do i1=1,nl1
                   if ((iproc==0).and.(mod(i1,10)==0)) print *,'     ind',i1
                   do n=1,nt_tzt
                      call fast_mvse_c(var_ct(i1,n,:),sp%d(nd)%lbloc,sp%ncapon,eps)
                   enddo
                enddo
             else
                ! DFFT in z (third direction after transposition)
                ! ---------
                if (iproc==0) print *,'   - ZFFT+shift'
                call ZFFT_shift(var_ct,3)
                var_ct=var_ct/sp%d(nd)%lbloc ! normalization of FFT
             endif

             ! Reverse transposition
             ! ---------------------
             if (iproc==0) print *,'   - zt-untranspose'
             call untranspose_zt_c

          else ! (serial implementation)
             ! Windowing
             ! ---------
             call windowing_c(var_c,nd,nd)
             
             if ((nd==sp%dim).and.(sp%is_capon)) then
                ! Capon's transform in z
                ! ----------------------
                if (iproc==0) print *,'   - Capon''s estimator (MVSE)'
                do i1=1,nl1
                   if ((iproc==0).and.(mod(i1,10)==0)) print *,'     ind',i1
                   do n=1,ngt
                      call fast_mvse_c(var_c(i1,:,n),sp%d(nd)%lbloc,sp%ncapon,eps)
                   enddo
                enddo
             else
                ! DFFT in z
                ! ---------
                if (iproc==0) print *,'   - ZFFT+shift'
                call ZFFT_shift(var_c,sp%d(nd)%i)
                var_c=var_c/sp%d(nd)%lbloc ! normalization of FFT
             endif

          endif

       case('t')          
          ! Windowing
          ! ---------
          call windowing_c(var_c,nd,nd)

          if ((nd==sp%dim).and.(sp%is_capon)) then
             ! Capon's transform in t
             ! ----------------------
             if (iproc==0) print *,'   - Capon''s estimator (MVSE)'
             do i2=1,nl2
                if ((iproc==0).and.(mod(i2,10)==0)) print *,'     ind',i2
                do i1=1,nl1
                   call fast_mvse_c(var_c(i1,i2,:),sp%d(nd)%lbloc,sp%ncapon,eps)
                enddo
             enddo
          else
             ! DFFT in t
             ! ---------
             if (iproc==0) print *,'   - ZFFT+shift'
             call ZFFT_shift(var_c,3)
             var_c=var_c/sp%d(nd)%lbloc ! normalization of FFT
          endif

       end select
       
    enddo
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! end loop on spectrum directions
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  end subroutine compute_multiD_spectrum
   
  !============================================================================
  subroutine spectrum_normalization
  !============================================================================
    !> Normalization of frequency-wavenumber spectra
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i1,i2,n,nd
    real(wp) :: rms2_new,fac,facm
    real(wp), dimension(n_in) :: rms2_nin
    ! -------------------------------------------------------------------------

    if (iproc==0) print *,'~> spectrum normalization'

    ! Apply correcting factor for windowing in each directions
    ! ========================================================
    do nd=1,sp%dim
       varm_r=varm_r/sp%d(nd)%Cw
    enddo
    
    ! Integration of spectra in each directions (except non-homogeneous one)
    ! =========================================   
    ! initialize
    facm=0.0_wp
    rms2_nin=0.0_wp

    ! if one direction is inhomogeneous
    if (i_in>0) then
       
       ! compute rms squared
       ! -------------------
       select case (i_in)
       case(1)
          ! compute local rms squared per proc
          ! ----------------------------------
                    
          ! for 1D spectrum, integration on half because real-valued FFT
          if (sp%dim==1) then
             if (sp%d(1)%i==2) then
                if (ndom_sp==1) then
                   do i1=1,nl1
                      do i2=1,nl2/2
                         do n=1,nl3
                            rms2_nin(i1)=rms2_nin(i1)+varm_r(i1,i2,n)
                         enddo
                      enddo
                   enddo
                else
                   ! in parallel, module is stored on first half of procs
                   if (coord(i_sp)<ndom_sp/2) then
                      do i1=1,nl1
                         do i2=1,nl2
                            do n=1,nl3
                               rms2_nin(i1)=rms2_nin(i1)+varm_r(i1,i2,n)
                            enddo
                         enddo
                      enddo
                   endif
                endif
             else
                ! for time dimension
                do i1=1,nl1
                   do i2=1,nl2
                      do n=1,nl3/2
                         rms2_nin(i1)=rms2_nin(i1)+varm_r(i1,i2,n)
                      enddo
                   enddo
                enddo
             endif
          else ! 2D spectra 
             do i1=1,nl1
                do i2=1,nl2
                   do n=1,nl3
                      rms2_nin(i1)=rms2_nin(i1)+varm_r(i1,i2,n)
                   enddo
                enddo
             enddo
          endif

          ! compute global rms squared: MPI reduce
          ! --------------------------------------
          call MPI_ALLREDUCE(MPI_IN_PLACE,rms2_nin,n_in,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_in,info)
          rms2_nin=rms2_nin/n_av
             
          do i1=1,nl1
             ! avoid division by zero
             ! (for velocity components that are null at the wall)
             if (rms2_nin(i1)==0.0_wp) rms2_nin(i1)=1.0_wp
             
             ! apply renormalization
             ! ---------------------
             ! factor
             fac=rms2m_in(i1)/rms2_nin(i1)
             facm=facm+fac
             ! renormalization
             varm_r(i1,:,:)=varm_r(i1,:,:)*fac
          enddo
          
       case(2)
          ! compute local rms squared per proc
          ! ----------------------------------
          
          ! for 1D spectrum, integration on half because real-valued FFT
          if (sp%dim==1) then
             if (sp%d(1)%i==1) then
                if (ndom_sp==1) then
                   do i2=1,nl2
                      do i1=1,nl1/2
                         do n=1,nl3
                            rms2_nin(i2)=rms2_nin(i2)+varm_r(i1,i2,n)
                         enddo
                      enddo
                   enddo
                else
                   ! in parallel, module is stored on first half of procs
                   if (coord(i_sp)<ndom_sp/2) then
                      do i2=1,nl2
                         do i1=1,nl1
                            do n=1,nl3
                               rms2_nin(i2)=rms2_nin(i2)+varm_r(i1,i2,n)
                            enddo
                         enddo
                      enddo
                   endif
                endif
             else
                ! for time dimension
                do i2=1,nl2
                   do i1=1,nl1
                      do n=1,nl3/2
                         rms2_nin(i2)=rms2_nin(i2)+varm_r(i1,i2,n)
                      enddo
                   enddo
                enddo
             endif
          else ! 2D spectra      
             do i2=1,nl2
                do n=1,nl3
                   do i1=1,nl1
                      rms2_nin(i2)=rms2_nin(i2)+varm_r(i1,i2,n)
                   enddo
                enddo
             enddo
          endif

          ! compute global rms squared: MPI reduce
          ! --------------------------------------
          call MPI_ALLREDUCE(MPI_IN_PLACE,rms2_nin,n_in,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_in,info)
          rms2_nin=rms2_nin/n_av

          do i2=1,nl2
             ! avoid division by zero
             ! (for velocity components that are null at the wall)
             if (rms2_nin(i2)==0.0_wp) rms2_nin(i2)=1.0_wp
             
             ! apply renormalization
             ! ---------------------
             ! factor
             fac=rms2m_in(i2)/rms2_nin(i2)
             facm=facm+fac
             ! renormalization
             varm_r(:,i2,:)=varm_r(:,i2,:)*fac
          enddo
          
       case(3)
          ! compute local rms squared per proc
          ! ----------------------------------
          
          ! for 1D spectrum, integration on half because real-valued FFT
          if (sp%dim==1) then
             if (sp%d(1)%i==1) then
                if (ndom_sp==1) then
                   do n=1,nl3
                      do i1=1,nl1/2
                         do i2=1,nl2
                            rms2_nin(n)=rms2_nin(n)+varm_r(i1,i2,n)
                         enddo
                      enddo
                   enddo
                else
                   ! in parallel, module is stored on first half of procs
                   if (coord(i_sp)<ndom_sp/2) then
                      do n=1,nl3
                         do i1=1,nl1
                            do i2=1,nl2
                               rms2_nin(n)=rms2_nin(n)+varm_r(i1,i2,n)
                            enddo
                         enddo
                      enddo
                   endif
                endif
             else ! sp%d(1)%i=2
                if (ndom_sp==1) then
                   do n=1,nl3
                      do i1=1,nl1
                         do i2=1,nl2/2
                            rms2_nin(n)=rms2_nin(n)+varm_r(i1,i2,n)
                         enddo
                      enddo
                   enddo
                else
                   ! in parallel, module is stored on first half of procs
                   if (coord(i_sp)<ndom_sp/2) then
                      do n=1,nl3
                         do i1=1,nl1
                            do i2=1,nl2
                               rms2_nin(n)=rms2_nin(n)+varm_r(i1,i2,n)
                            enddo
                         enddo
                      enddo
                   endif
                endif
             endif
          else ! 2D spectra
             do n=1,nl3
                do i1=1,nl1
                   do i2=1,nl2
                      rms2_nin(n)=rms2_nin(n)+varm_r(i1,i2,n)
                   enddo
                enddo
             enddo
          endif
             
          ! compute global rms squared: MPI reduce
          ! --------------------------------------
          call MPI_ALLREDUCE(MPI_IN_PLACE,rms2_nin,n_in,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_in,info)
          rms2_nin=rms2_nin/n_av
             
          do n=1,nl3
             ! avoid division by zero
             ! (for velocity components that are null at the wall)
             if (rms2_nin(n)==0.0_wp) rms2_nin(n)=1.0_wp
             
             ! apply renormalization
             ! ---------------------
             ! factor
             fac=rms2m_in(n)/rms2_nin(n)
             facm=facm+fac
             ! renormalization
             varm_r(:,:,n)=varm_r(:,:,n)*fac
          enddo
       end select

       ! averaged factor
       ! ---------------
       facm=facm/n_in
    
    else

       ! local rms squared
       ! -----------------
       rms2_new=0.0_wp

       do n=1,nl3
          do i1=1,nl1
             do i2=1,nl2
                rms2_new=rms2_new+varm_r(i1,i2,n)
             enddo
          enddo
       enddo

       ! global rms squared
       ! ------------------
       call MPI_ALLREDUCE(MPI_IN_PLACE,rms2_new,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
       
       ! apply normalization
       ! -------------------
       ! factor
       fac=rms2m/rms2_new
       facm=fac
       ! renormalization
       varm_r=varm_r*fac
       
    endif
    
    ! Normalization: divide by resolution (include one/two-sided factor)
    ! ===================================
    varm_r=varm_r/sp%dk

    ! Print normalization averaged factor
    ! ===================================
    if (iproc==0) print *,'   - normalization factor:',facm

  end subroutine spectrum_normalization

end module mod_pp_sp_eval
