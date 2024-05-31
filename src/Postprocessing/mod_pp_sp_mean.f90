!==============================================================================
module mod_pp_sp_mean
!==============================================================================
  !> Module to compute autospectrum
!==============================================================================
  use mod_mpi_part ! <- for iproc, REDUCE operations, COMMXY/XZ/YZ for 1var
  use mod_pp_var   ! <- 3D arrays and dimensions
  implicit none
  ! ---------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------

contains

  !============================================================================
  subroutine subtract_mean_3var
  !============================================================================
    !> subtract mean + compute rms in PSD/averaging directions
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i1,i2,n
    real(wp) :: mean
    ! -------------------------------------------------------------------------
    
    ! Subtraction of average
    ! ======================
    if (iproc==0) print *,'~> subtract global average ...'

    ! local mean
    ! ----------
    mean=0.0_wp
    do n=1,nl3
       do i2=1,nl2
          do i1=1,nl1
             mean=mean+var_r(i1,i2,n)
          enddo
       enddo
    enddo

    ! global mean
    ! -----------
    call MPI_ALLREDUCE(MPI_IN_PLACE,mean,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
    mean=mean/nl1/nl2/nl3/nproc

    ! subtract the mean value
    ! -----------------------
    do n=1,nl3
       do i2=1,nl2
          do i1=1,nl1
             var_r(i1,i2,n)=var_r(i1,i2,n)-mean
          enddo
       enddo
    enddo
       
    ! Compute rms squared (used for normalization)
    ! ===================
    
    ! local rms
    ! ---------
    rms2=0.0_wp
    do n=1,nl3
       do i2=1,nl2
          do i1=1,nl1
             rms2=rms2+var_r(i1,i2,n)**2
          enddo
       enddo
    enddo
    
    ! global rms
    ! ---------
    call MPI_ALLREDUCE(MPI_IN_PLACE,rms2,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info) 
    rms2=rms2/ng_av
    
    if (iproc==0) print *,'~> rms value:',sqrt(rms2)
    
  end subroutine subtract_mean_3var
  
  !============================================================================
  subroutine subtract_mean_2var
  !============================================================================
    !> subtract mean in PSD direction + averaging direction
  !============================================================================
    use mod_grid
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i1,i2,n
    real(wp) :: mean
    ! -------------------------------------------------------------------------
    
    ! Subtraction of average
    ! ======================
    if (iproc==0) then
       if (sp%dim==1) print *,'~> subtract '//sp%d(1)%name//dir_av(1)//'-average ...'
       if (sp%dim==2) print *,'~> subtract '//sp%d(1)%name//sp%d(2)%name//'-average ...'
    endif
    
    select case(i_in)
       
    ! inhomogeneous direction is 1 (y for yz-plane/x for xz-plane or xy-plane)
    ! ----------------------------
    case (1)
       do i1=1,nl1
          ! compute local mean
          mean=0.0_wp
          do n=1,nl3
             do i2=1,nl2
                mean=mean+var_r(i1,i2,n)
             enddo
          enddo
          
          ! global mean: MPI comm
          call MPI_ALLREDUCE(MPI_IN_PLACE,mean,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_in,info)
          
          ! divide by number of samples
          mean=mean/ng_av
          
          ! substract global mean
          do n=1,nl3
             do i2=1,nl2
                var_r(i1,i2,n)=var_r(i1,i2,n)-mean
             enddo
          enddo
       enddo

    ! inhomogeneous direction is 2 (z for yz-plane or xz-plane/y for xy-plane)
    ! ----------------------------
    case (2) 
       do i2=1,nl2
          ! compute local mean
          mean=0.0_wp
          do n=1,nl3
             do i1=1,nl1
                mean=mean+var_r(i1,i2,n)
             enddo
          enddo

          ! global mean: MPI comm
          call MPI_ALLREDUCE(MPI_IN_PLACE,mean,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_in,info)
          
          ! divide by number of samples
          mean=mean/ng_av
          
          ! substract global mean
          do n=1,nl3
             do i1=1,nl1
                var_r(i1,i2,n)=var_r(i1,i2,n)-mean
             enddo
          enddo
       enddo

    ! inhomogeneous direction is 3 (t)
    ! ----------------------------
    case (3) 
       do n=1,nl3
          ! compute local mean
          mean=0.0_wp
          do i1=1,nl1
             do i2=1,nl2
                mean=mean+var_r(i1,i2,n)
             enddo
          enddo

          ! global mean: MPI comm
          call MPI_ALLREDUCE(MPI_IN_PLACE,mean,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_in,info)
          
          ! divide by number of samples
          mean=mean/ng_av
          
          ! substract global mean
          do i1=1,nl1
             do i2=1,nl2
                var_r(i1,i2,n)=var_r(i1,i2,n)-mean
             enddo
          enddo
       enddo

    end select
       
    ! Compute rms squared (used for normalization)
    ! ===================
    rms2_in=0.0_wp
    skew_in=0.0_wp
    kurt_in=0.0_wp
    
    select case(i_in)
       
    ! inhomogeneous direction is 1 (y for yz-plane/x for xz-plane or xy-plane)
    ! ----------------------------
    case (1)
       do i1=1,nl1
          ! compute local rms
          do n=1,nl3
             do i2=1,nl2
                rms2_in(i1)=rms2_in(i1)+var_r(i1,i2,n)**2
                skew_in(i1)=skew_in(i1)+var_r(i1,i2,n)**3
                kurt_in(i1)=kurt_in(i1)+var_r(i1,i2,n)**4
             enddo
          enddo
          
          ! global rms: MPI comm
          call MPI_ALLREDUCE(MPI_IN_PLACE,rms2_in(i1),1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_in,info)
          call MPI_ALLREDUCE(MPI_IN_PLACE,skew_in(i1),1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_in,info)
          call MPI_ALLREDUCE(MPI_IN_PLACE,kurt_in(i1),1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_in,info)
          
          ! avoid division by zero
          ! (for velocity components that are null at the wall)
          if (rms2_in(i1)==0.0_wp) rms2_in(i1)=1.0_wp
       enddo

    ! inhomogeneous direction is 2 (z for yz-plane or xz-plane/y for xy-plane)
    ! ----------------------------
    case (2) 
       do i2=1,nl2
          ! compute local rms
          mean=0.0_wp
          do n=1,nl3
             do i1=1,nl1
                rms2_in(i2)=rms2_in(i2)+var_r(i1,i2,n)**2
                skew_in(i2)=skew_in(i2)+var_r(i1,i2,n)**3
                kurt_in(i2)=kurt_in(i2)+var_r(i1,i2,n)**4
             enddo
          enddo

          ! global rms: MPI comm
          call MPI_ALLREDUCE(MPI_IN_PLACE,rms2_in(i2),1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_in,info)
          call MPI_ALLREDUCE(MPI_IN_PLACE,skew_in(i2),1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_in,info)
          call MPI_ALLREDUCE(MPI_IN_PLACE,kurt_in(i2),1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_in,info)
          
          ! avoid division by zero
          ! (for velocity components that are null at the wall)
          if (rms2_in(i2)==0.0_wp) rms2_in(i2)=1.0_wp         
       enddo

    ! inhomogeneous direction is 3 (t)
    ! ----------------------------
    case (3) 
       do n=1,nl3
          ! compute local rms
          mean=0.0_wp
          do i1=1,nl1
             do i2=1,nl2
                rms2_in(n)=rms2_in(n)+var_r(i1,i2,n)**2
                skew_in(n)=skew_in(n)+var_r(i1,i2,n)**3
                kurt_in(n)=kurt_in(n)+var_r(i1,i2,n)**4
             enddo
          enddo

          ! global rms: MPI comm
          call MPI_ALLREDUCE(MPI_IN_PLACE,rms2_in(n),1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_in,info)
          call MPI_ALLREDUCE(MPI_IN_PLACE,skew_in(n),1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_in,info)
          call MPI_ALLREDUCE(MPI_IN_PLACE,kurt_in(n),1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_in,info)
          
          ! avoid division by zero
          ! (for velocity components that are null at the wall)
          if (rms2_in(n)==0.0_wp) rms2_in(n)=1.0_wp
       enddo
     
   end select
   
   ! divide by number of samples
   rms2_in=rms2_in/ng_av
   skew_in=skew_in/ng_av
   kurt_in=kurt_in/ng_av

   ! skewness
   skew_in=skew_in/rms2_in**(3.0_wp/2.0_wp) ! skewness <p'^3>/<p'^2>^3/2

   ! kurtosis (or flatness)
   kurt_in=kurt_in/rms2_in**2 ! kurtosis <p'^4>/<p'^2>^2

  end subroutine subtract_mean_2var
    
  !============================================================================
  subroutine subtract_mean_1var
  !============================================================================
    !> subtract mean in PSD direction
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i1,i2,n
    real(wp) :: mean
    ! -------------------------------------------------------------------------
    
    ! Subtraction of average
    ! ======================
    if (iproc==0) print *,'~> subtract '//sp%d(1)%name//'-average ...'

    select case(sp%d(1)%i)
       
    ! Averaging direction is 1 (y for yz-plane/x for xz-plane or xy-plane)
    ! ------------------------
    case (1) 
       do n=1,nl3
          do i2=1,nl2
             ! compute local mean
             mean=0.0_wp
             do i1=1,nl1
                mean=mean+var_r(i1,i2,n)
             enddo

             ! compute global mean
             if (sp%d(1)%name=='x') then
                ! (MPI comm along x-direction)
                call MPI_ALLREDUCE(MPI_IN_PLACE,mean,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMMYZ,info)
                mean=mean/nl1/ndomx
             endif
             if (sp%d(1)%name=='y') then
                ! (MPI comm along y-direction)
                call MPI_ALLREDUCE(MPI_IN_PLACE,mean,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXZ,info)
                mean=mean/nl1/ndomy
             endif

             ! substract global mean
             do i1=1,nl1
                var_r(i1,i2,n)=var_r(i1,i2,n)-mean
             enddo
          enddo
       enddo

    ! Averaging direction is 2 (z for yz-plane or xz-plane/y for xy-plane)
    ! ------------------------
    case (2) 
       do n=1,nl3
          do i1=1,nl1
             ! compute local mean
             mean=0.0_wp
             do i2=1,nl2
                mean=mean+var_r(i1,i2,n)
             enddo

             ! compute global mean
             if (sp%d(1)%name=='y') then
                ! (MPI comm along y-direction)
                call MPI_ALLREDUCE(MPI_IN_PLACE,mean,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXZ,info)
                mean=mean/nl2/ndomy
             endif
             if (sp%d(1)%name=='z') then
                ! (MPI comm along z-direction)
                call MPI_ALLREDUCE(MPI_IN_PLACE,mean,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXY,info)
                mean=mean/nl2/ndomz
             endif
             
             ! substract global mean
             do i2=1,nl2
                var_r(i1,i2,n)=var_r(i1,i2,n)-mean
             enddo
          enddo
       enddo

    ! Averaging direction is 3 (t)
    ! ------------------------
    case (3) 
       do i1=1,nl1
          do i2=1,nl2
             ! compute global mean
             mean=0.0_wp
             do n=1,nl3
                mean=mean+var_r(i1,i2,n)
             enddo
             mean=mean/nl3

             ! substract global mean
             do n=1,nl3
                var_r(i1,i2,n)=var_r(i1,i2,n)-mean
             enddo
          enddo
       enddo

    end select

    ! Compute rms squared (used for normalization)
    ! ===================

    select case(sp%d(1)%i)

    ! Averaging direction is 1 (y for yz-plane/x for xz-plane or xy-plane)
    ! ------------------------
    case (1)
       if (.not.allocated(rms2_1var)) allocate(rms2_1var(nl2,nl3))
       rms2_1var=0.0_wp
       do n=1,nl3
          do i2=1,nl2
             ! compute local rms
             do i1=1,nl1
                rms2_1var(i2,n) = rms2_1var(i2,n) + var_r(i1,i2,n)**2
             enddo

             ! compute global rms
             if (sp%d(1)%name=='x') then
                ! (MPI comm along x-direction)
                call MPI_ALLREDUCE(MPI_IN_PLACE,rms2_1var(i2,n),1,MPI_DOUBLE_PRECISION,MPI_SUM,COMMYZ,info)
                rms2_1var(i2,n) = rms2_1var(i2,n)/nl1/ndomx
             endif
             if (sp%d(1)%name=='y') then
                ! (MPI comm along y-direction)
                call MPI_ALLREDUCE(MPI_IN_PLACE,rms2_1var(i2,n),1,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXZ,info)
                rms2_1var(i2,n) = rms2_1var(i2,n)/nl1/ndomy
             endif

             ! avoid division by zero
             ! (for velocity components that are null at the wall)
             if (rms2_1var(i2,n)==0.0_wp) rms2_1var(i2,n) = 1.0_wp
          enddo
       enddo

    ! Averaging direction is 2 (z for yz-plane or xz-plane/y for xy-plane)
    ! ------------------------
    case (2)
       if (.not.allocated(rms2_1var)) allocate(rms2_1var(nl1,nl3))
       rms2_1var=0.0_wp
       do n=1,nl3
          do i1=1,nl1
             ! compute local rms
             do i2=1,nl2
                rms2_1var(i1,n) = rms2_1var(i1,n) + var_r(i1,i2,n)**2
             enddo

             ! compute global rms
             if (sp%d(1)%name=='y') then
                ! (MPI comm along y-direction)
                call MPI_ALLREDUCE(MPI_IN_PLACE,rms2_1var(i1,n),1,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXZ,info)
                rms2_1var(i1,n) = rms2_1var(i1,n)/nl2/ndomy
             endif
             if (sp%d(1)%name=='z') then
                ! (MPI comm along z-direction)
                call MPI_ALLREDUCE(MPI_IN_PLACE,rms2_1var(i1,n),1,MPI_DOUBLE_PRECISION,MPI_SUM,COMMXY,info)
                rms2_1var(i1,n) = rms2_1var(i1,n)/nl2/ndomz
             endif

             ! avoid division by zero
             ! (for velocity components that are null at the wall)
             if (rms2_1var(i1,n)==0.0_wp) rms2_1var(i1,n) = 1.0_wp
          enddo
       enddo

    ! Averaging direction is 3 (t)
    ! ------------------------
    case (3)
       if (.not.allocated(rms2_1var)) allocate(rms2_1var(nl1,nl2))
       rms2_1var=0.0_wp
       do i1=1,nl1
          do i2=1,nl2
             ! compute global rms
             do n=1,nl3
                rms2_1var(i1,i2) = rms2_1var(i1,i2) + var_r(i1,i2,n)**2
             enddo
             rms2_1var(i1,i2) = rms2_1var(i1,i2)/nl3

             ! avoid division by zero
             ! (for velocity components that are null at the wall)
             if (rms2_1var(i1,i2)==0.0_wp) rms2_1var(i1,i2) = 1.0_wp
          enddo
       enddo


    end select


    ! Compute rms squared (used for compatibility)
    !  ~> write in files for spectra but not for correl
    ! ===================
    if (type_pp.ne.4) then

       rms2_in=0.0_wp
       skew_in=0.0_wp
       kurt_in=0.0_wp

       select case(i_in)

       ! inhomogeneous direction is 1 (y for yz-plane/x for xz-plane or xy-plane)
       ! ----------------------------
       case (1)
          do i1=1,nl1
             ! compute local rms
             do n=1,nl3
                do i2=1,nl2
                   rms2_in(i1)=rms2_in(i1)+var_r(i1,i2,n)**2
                   skew_in(i1)=skew_in(i1)+var_r(i1,i2,n)**3
                   kurt_in(i1)=kurt_in(i1)+var_r(i1,i2,n)**4
                enddo
             enddo

             ! global rms: MPI comm
             call MPI_ALLREDUCE(MPI_IN_PLACE,rms2_in(i1),1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_in,info)
             call MPI_ALLREDUCE(MPI_IN_PLACE,skew_in(i1),1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_in,info)
             call MPI_ALLREDUCE(MPI_IN_PLACE,kurt_in(i1),1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_in,info)

             ! avoid division by zero
             ! (for velocity components that are null at the wall)
             if (rms2_in(i1)==0.0_wp) rms2_in(i1)=1.0_wp
          enddo

       ! inhomogeneous direction is 2 (z for yz-plane or xz-plane/y for xy-plane)
       ! ----------------------------
       case (2)
          do i2=1,nl2
             ! compute local rms
             do n=1,nl3
                do i1=1,nl1
                   rms2_in(i2)=rms2_in(i2)+var_r(i1,i2,n)**2
                   skew_in(i2)=skew_in(i2)+var_r(i1,i2,n)**3
                   kurt_in(i2)=kurt_in(i2)+var_r(i1,i2,n)**4
                enddo
             enddo

             ! global rms: MPI comm
             call MPI_ALLREDUCE(MPI_IN_PLACE,rms2_in(i2),1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_in,info)
             call MPI_ALLREDUCE(MPI_IN_PLACE,skew_in(i2),1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_in,info)
             call MPI_ALLREDUCE(MPI_IN_PLACE,kurt_in(i2),1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_in,info)

             ! avoid division by zero
             ! (for velocity components that are null at the wall)
             if (rms2_in(i2)==0.0_wp) rms2_in(i2)=1.0_wp
          enddo

       ! inhomogeneous direction is 3 (t)
       ! ----------------------------
       case (3)
          do n=1,nl3
             ! compute local rms
             do i1=1,nl1
                do i2=1,nl2
                   rms2_in(n)=rms2_in(n)+var_r(i1,i2,n)**2
                   skew_in(n)=skew_in(n)+var_r(i1,i2,n)**3
                   kurt_in(n)=kurt_in(n)+var_r(i1,i2,n)**4
                enddo
             enddo

             ! global rms: MPI comm
             call MPI_ALLREDUCE(MPI_IN_PLACE,rms2_in(n),1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_in,info)
             call MPI_ALLREDUCE(MPI_IN_PLACE,skew_in(n),1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_in,info)
             call MPI_ALLREDUCE(MPI_IN_PLACE,kurt_in(n),1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_in,info)

             ! avoid division by zero
             ! (for velocity components that are null at the wall)
             if (rms2_in(n)==0.0_wp) rms2_in(n)=1.0_wp
          enddo

       end select

       ! divide by number of samples
       rms2_in=rms2_in/ng_av
       skew_in=skew_in/ng_av
       kurt_in=kurt_in/ng_av

       ! skewness
       skew_in=skew_in/rms2_in**(3.0_wp/2.0_wp) ! skewness <p'^3>/<p'^2>^3/2

       ! kurtosis (or flatness)
       kurt_in=kurt_in/rms2_in**2 ! kurtosis <p'^4>/<p'^2>^2

    endif
       
  end subroutine subtract_mean_1var
    
end module mod_pp_sp_mean
