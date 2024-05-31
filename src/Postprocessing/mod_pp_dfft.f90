!==============================================================================
module mod_pp_dfft
!==============================================================================
  !> Fourier's transform routines based on DFFTPACK
!==============================================================================
  use precision
  implicit none
  ! ---------------------------------------------------------------------------
  integer, dimension(3,3) :: np
  ! ---------------------------------------------------------------------------

contains

  !============================================================================
  subroutine init_dim_DFFT
  !============================================================================
    !> Initialization of switch per direction
    !> for FFT of three-dimensional arrays
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    ! -------------------------------------------------------------------------

    ! if dim==1
    np(1,1)=1
    np(2,1)=2
    np(3,1)=3
    ! if dim==2
    np(1,2)=2
    np(2,2)=1
    np(3,2)=3
    ! if dim==3
    np(1,3)=3
    np(2,3)=1
    np(3,3)=2

  end subroutine init_dim_DFFT

  !============================================================================
  subroutine DFFT(PHI,dim)
  !============================================================================
    !> Direct FFT of 3-D array with respect to its direction 'dim'
    !> - Real version -
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer, intent(in) :: dim
    real(wp), dimension(:,:,:) :: PHI
    ! -------------------------------------------------------------------------
    integer :: i,j,l,n,m
    real(wp), dimension(:), allocatable :: tab
    real(wp), dimension(:), allocatable :: wsave ! work array of FFTPACK
    ! -------------------------------------------------------------------------

    ! Array dimensions
    ! ================
    l=size(PHI,np(1,dim))
    n=size(PHI,np(2,dim))
    m=size(PHI,np(3,dim))

    ! Init working array
    ! ==================
    allocate(tab(l),wsave(2*l+15))
    call dffti(l,wsave)

    ! Compute FFT
    ! ===========
    if (dim==1) then

       do i=1,n
          do j=1,m
             tab=PHI(:,i,j)
             call dfftf(l,tab,wsave)
             PHI(:,i,j)=tab       
          enddo
       enddo

    elseif (dim==2) then

       do i=1,n
          do j=1,m
             tab=PHI(i,:,j)
             call dfftf(l,tab,wsave)
             PHI(i,:,j)=tab       
          enddo
       enddo

    elseif (dim==3) then

       do i=1,n
          do j=1,m
             tab=PHI(i,j,:)
             call dfftf(l,tab,wsave)
             PHI(i,j,:)=tab       
          enddo
       enddo

    endif

    deallocate(tab,wsave)

  end subroutine DFFT

  !============================================================================
  subroutine INVDFFT(PHI,dim)
  !============================================================================
    !> Inverse FFT of 3-D array with respect to its direction 'dim'
    !> - Real version -
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer, intent(in) :: dim
    real(wp), dimension(:,:,:) :: PHI
    ! -------------------------------------------------------------------------
    integer :: i,j,l,n,m
    real(wp), dimension(:), allocatable :: tab
    real(wp), dimension(:), allocatable :: wsave ! work array of FFTPACK
    ! -------------------------------------------------------------------------

    ! Array dimensions
    ! ================
    l=size(PHI,np(1,dim))
    n=size(PHI,np(2,dim))
    m=size(PHI,np(3,dim))

    ! Init working array
    ! ==================
    allocate(tab(l),wsave(2*l+15))
    call dffti(l,wsave)

    ! Compute FFT
    ! ===========
    if (dim==1) then

       do i=1,n
          do j=1,m
             tab=PHI(:,i,j)
             call dfftb(l,tab,wsave)
             PHI(:,i,j)=tab       
          enddo
       enddo

    elseif (dim==2) then

       do i=1,n
          do j=1,m
             tab=PHI(i,:,j)
             call dfftb(l,tab,wsave)
             PHI(i,:,j)=tab       
          enddo
       enddo

    elseif (dim==3) then

       do i=1,n
          do j=1,m
             tab=PHI(i,j,:)
             call dfftb(l,tab,wsave)
             PHI(i,j,:)=tab       
          enddo
       enddo

    endif

    deallocate(tab,wsave)

  end subroutine INVDFFT

  !============================================================================
  subroutine ZFFT(PHI,dim)
  !============================================================================
    !> Direct FFT of 3-D array with respect to its direction 'dim'
    !> - Complex version -
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer, intent(in) :: dim
    complex(wp), dimension(:,:,:) :: PHI
    ! -------------------------------------------------------------------------
    integer :: i,j,l,n,m
    complex(wp), dimension(:), allocatable :: tab
    real(wp), dimension(:), allocatable :: wsave ! work array of FFTPACK
    ! -------------------------------------------------------------------------

    ! Array dimensions
    ! ================
    l=size(PHI,np(1,dim))
    n=size(PHI,np(2,dim))
    m=size(PHI,np(3,dim))

    ! Init working array
    ! ==================
    allocate(tab(l),wsave(4*l+15))
    call zffti(l,wsave)

    ! Compute FFT
    ! ===========
    if (dim==1) then

       do i=1,n
          do j=1,m
             tab=PHI(:,i,j)
             call zfftf(l,tab,wsave)
             PHI(:,i,j)=tab       
          enddo
       enddo

    elseif (dim==2) then

       do i=1,n
          do j=1,m
             tab=PHI(i,:,j)
             call zfftf(l,tab,wsave)
             PHI(i,:,j)=tab       
          enddo
       enddo

    elseif (dim==3) then

       do i=1,n
          do j=1,m
             tab=PHI(i,j,:)
             call zfftf(l,tab,wsave)
             PHI(i,j,:)=tab       
          enddo
       enddo

    endif

    deallocate(tab,wsave)

  end subroutine ZFFT

  !============================================================================
  subroutine INVZFFT(PHI,dim)
  !============================================================================
    !> Inverse FFT of 3-D array with respect to its direction 'dim'
    !> - Complex version -
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer, intent(in) :: dim
    complex(wp), dimension(:,:,:) :: PHI
    ! -------------------------------------------------------------------------
    integer :: i,j,l,n,m
    complex(wp), dimension(:), allocatable :: tab
    real(wp), dimension(:), allocatable :: wsave ! work array of FFTPACK
    ! -------------------------------------------------------------------------

    ! Array dimensions
    ! ================
    l=size(PHI,np(1,dim))
    n=size(PHI,np(2,dim))
    m=size(PHI,np(3,dim))

    ! Init working array
    ! ==================
    allocate(tab(l),wsave(4*l+15))
    call zffti(l,wsave)

    ! Compute FFT
    ! ===========
    if (dim==1) then

       do i=1,n
          do j=1,m
             tab=PHI(:,i,j)
             call zfftb(l,tab,wsave)
             PHI(:,i,j)=tab       
          enddo
       enddo

    elseif (dim==2) then

       do i=1,n
          do j=1,m
             tab=PHI(i,:,j)
             call zfftb(l,tab,wsave)
             PHI(i,:,j)=tab       
          enddo
       enddo

    elseif (dim==3) then

       do i=1,n
          do j=1,m
             tab=PHI(i,j,:)
             call zfftb(l,tab,wsave)
             PHI(i,j,:)=tab       
          enddo
       enddo

    endif

    deallocate(tab,wsave)

  end subroutine INVZFFT

  !============================================================================
  subroutine ZFFT_shift(PHI,dim)
  !============================================================================
    !> Direct FFT of 3-D array with respect to its direction 'dim'
    !> + Shift zero-frequency component to center of spectrum
    !> - Complex version -
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer, intent(in) :: dim
    complex(wp), dimension(:,:,:) :: PHI
    ! -------------------------------------------------------------------------
    integer :: i,j,n,m,l,q,r
    integer, dimension(:), allocatable :: idx
    complex(wp), dimension(:), allocatable :: tab
    complex(wp), dimension(size(PHI,1),size(PHI,2),size(PHI,3)) :: PHIb
    real(wp), dimension(:), allocatable :: wsave ! work array of FFTPACK
    ! -------------------------------------------------------------------------

    ! Array dimensions
    ! ================
    l=size(PHI,np(1,dim))
    n=size(PHI,np(2,dim))
    m=size(PHI,np(3,dim))

    ! Init working array
    ! ==================
    allocate(tab(l),idx(l),wsave(4*l+15))
    call zffti(l,wsave)

    ! Shift index
    ! ===========
    q = ceiling(real(l)/2.)
    r=mod(l,2)

    do i=1,q-1
       idx(i)=1+q-r+i
    enddo
    do i=q,l
       idx(i)=i-q+1
    enddo

    ! Compute FFT
    ! ===========
    if (dim==1) then

       do i=1,n
          do j=1,m
             tab=PHI(:,i,j)
             call zfftf(l,tab,wsave)
             PHIb(:,i,j)=tab       
          enddo
       enddo

       ! shift zero-frequency component to center of spectrum
       ! ----------------------------------------------------
       do i=1,l
          PHI(i,:,:)=PHIb(idx(i),:,:)
       enddo

    elseif (dim==2) then

       do i=1,n
          do j=1,m
             tab=PHI(i,:,j)
             call zfftf(l,tab,wsave)
             PHIb(i,:,j)=tab       
          enddo
       enddo

       ! shift zero-frequency component to center of spectrum
       ! ----------------------------------------------------
       do i=1,l
          PHI(:,i,:)=PHIb(:,idx(i),:)
       enddo

    elseif (dim==3) then

       do i=1,n
          do j=1,m
             tab=PHI(i,j,:)
             call zfftf(l,tab,wsave)
             PHIb(i,j,:)=tab       
          enddo
       enddo

       ! shift zero-frequency component to center of spectrum
       ! ----------------------------------------------------
       do i=1,l
          PHI(:,:,i)=PHIb(:,:,idx(i))
       enddo

    endif

    deallocate(tab,idx,wsave)

  end subroutine ZFFT_shift

  !============================================================================
  subroutine DFFT1(PHI)
  !============================================================================
    !> Direct FFT of 1-D array
    !> - Real version -
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    real(wp), dimension(:) :: PHI
    ! -------------------------------------------------------------------------
    integer :: l
    real(wp), dimension(:), allocatable :: wsave ! work array of FFTPACK
    ! -------------------------------------------------------------------------

    ! Array size
    ! ==========
    l=size(PHI)

    ! Init working array
    ! ==================
    allocate(wsave(2*l+15))
    call dffti(l,wsave)

    ! Compute FFT
    ! ===========
    call dfftf(l,PHI,wsave)   
    deallocate(wsave)

  end subroutine DFFT1

  !============================================================================
  subroutine ZFFT1(PHI)
  !============================================================================
    !> Direct FFT of 1-D array
    !> - Complex version -
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    complex(wp), dimension(:) :: PHI
    ! -------------------------------------------------------------------------
    integer :: l
    real(wp), dimension(:), allocatable :: wsave ! work array of FFTPACK
    ! -------------------------------------------------------------------------

    ! Array size
    ! ==========
    l=size(PHI)

    ! Init working array
    ! ==================
    allocate(wsave(4*l+15))
    call zffti(l,wsave)

    ! Compute FFT
    ! ===========
    call zfftf(l,PHI,wsave)   
    deallocate(wsave)

  end subroutine ZFFT1

  !============================================================================
  subroutine sep_mod_phase(var,dim)
  !============================================================================
    !> Rearrange array from dfft in module and phase
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer, intent(in) :: dim
    real(wp), dimension(:,:,:) :: var
    ! -------------------------------------------------------------------------
    integer :: i,j,k,l,n,m,q
    real(wp), dimension(:,:,:), allocatable :: re_var,im_var
    ! -------------------------------------------------------------------------

    ! Array dimensions
    ! ================
    l=size(var,1)
    n=size(var,2)
    m=size(var,3)

    if (dim==1) then
       
       q=l/2
       allocate(re_var(q,n,m),im_var(q,n,m))
       
       ! Separate real and imaginary parts 
       ! =================================
       do j=1,n
          do k=1,m
             re_var(1,j,k)=var(1,j,k)
             im_var(1,j,k)=0.0_wp
             do i=2,q
                re_var(i,j,k)=var(2*i-2,j,k)
                im_var(i,j,k)=var(2*i-1,j,k)
             enddo
          enddo
       enddo

       ! Compute module and phase
       ! ========================
       var(1:q,:,:)=sqrt(re_var(1:q,:,:)**2+im_var(1:q,:,:)**2)
       var(q+1:l,:,:)=atan2(im_var(1:q,:,:),re_var(1:q,:,:))
       
    elseif (dim==2) then
       
       q=n/2
       allocate(re_var(l,q,m),im_var(l,q,m))
       
       ! Separate real and imaginary parts 
       ! =================================
       do i=1,l
          do k=1,m
             re_var(i,1,k)=var(i,1,k)
             im_var(i,1,k)=0.0_wp
             do j=2,q
                re_var(i,j,k)=var(i,2*j-2,k)
                im_var(i,j,k)=var(i,2*j-1,k)
             enddo
          enddo
       enddo

       ! Compute module and phase
       ! ========================
       var(:,1:q,:)=sqrt(re_var(:,1:q,:)**2+im_var(:,1:q,:)**2)
       var(:,q+1:n,:)=atan2(im_var(:,1:q,:),re_var(:,1:q,:))
       
    elseif (dim==3) then
       
       q=m/2
       allocate(re_var(l,n,q),im_var(l,n,q))
       
       ! Separate real and imaginary parts 
       ! =================================
       do i=1,l
          do j=1,n
             re_var(i,j,1)=var(i,j,1)
             im_var(i,j,1)=0.0_wp
             do k=2,q
                re_var(i,j,k)=var(i,j,2*k-2)
                im_var(i,j,k)=var(i,j,2*k-1)
             enddo
          enddo
       enddo

       ! Compute module and phase
       ! ========================
       var(:,:,1:q)=sqrt(re_var(:,:,1:q)**2+im_var(:,:,1:q)**2)
       var(:,:,q+1:m)=atan2(im_var(:,:,1:q),re_var(:,:,1:q))
    endif
    
    deallocate(re_var,im_var)

  end subroutine sep_mod_phase

end module mod_pp_dfft
