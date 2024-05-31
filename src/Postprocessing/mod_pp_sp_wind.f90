!==============================================================================
module mod_pp_sp_wind
!==============================================================================
  !> Module to compute autospectrum
!==============================================================================
  implicit none
  ! ---------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------

contains

  !============================================================================
  subroutine windowing_r(var,dir,dimw)
  !============================================================================
    !> Windowing and correction factor
  !============================================================================
    use mod_mpi    ! <- for iproc
    use mod_pp_var ! <- for sp structure
    implicit none
    ! -------------------------------------------------------------------------
    integer, intent(in) :: dir,dimw
    real(wp), dimension(:,:,:) :: var
    ! -------------------------------------------------------------------------
    integer :: n
    real(wp), allocatable, dimension(:) :: win
    ! -------------------------------------------------------------------------
    
    if (iproc==0) print *,'   - windowing'
    allocate(win(sp%d(dir)%lbloc))
    
    ! Windowing
    ! =========
    if (sp%d(dir)%type_win=='T') then
       ! Tukey window
       win = tukeywin(sp%d(dir)%lbloc,sp%d(dir)%param_win)
    elseif (sp%d(dir)%type_win=='K') then
       ! Kaiser-Bessel window
       win=kaiserbwin(sp%d(dir)%lbloc,sp%d(dir)%param_win)
    endif
    
    ! Correcting factor for windowing
    ! ===============================
    sp%d(dir)%Cw=window_factor(win,sp%d(dir)%lbloc)
    if (iproc==0) print *,'     (window factor:',sp%d(dir)%Cw,')'

    ! Apply windowing
    ! ===============
    select case (dimw)
    case(1)
       do n=1,sp%d(dir)%lbloc
          var(n,:,:)=var(n,:,:)*win(n)
       enddo
    case(2)
       do n=1,sp%d(dir)%lbloc
          var(:,n,:)=var(:,n,:)*win(n)
       enddo
    case(3)
       do n=1,sp%d(dir)%lbloc
          var(:,:,n)=var(:,:,n)*win(n)
       enddo
    end select
    
    deallocate(win)
   
  end subroutine windowing_r
  
  !============================================================================
  subroutine windowing_c(var,dir,dimw)
  !============================================================================
    !> Windowing and correction factor
  !============================================================================
    use mod_mpi    ! <- for iproc
    use mod_pp_var ! <- for sp structure
    implicit none
    ! -------------------------------------------------------------------------
    integer, intent(in) :: dir,dimw
    complex(wp), dimension(:,:,:) :: var
    ! -------------------------------------------------------------------------
    integer :: n
    real(wp), allocatable, dimension(:) :: win
    ! -------------------------------------------------------------------------
    
    if (iproc==0) print *,'   - windowing'
    allocate(win(sp%d(dir)%lbloc))
    
    ! Windowing
    ! =========
    if (sp%d(dir)%type_win=='T') then
       ! Tukey window
       win = tukeywin(sp%d(dir)%lbloc,sp%d(dir)%param_win)
    elseif (sp%d(dir)%type_win=='K') then
       ! Kaiser-Bessel window
       win=kaiserbwin(sp%d(dir)%lbloc,sp%d(dir)%param_win)
    endif
    
    ! Correcting factor for windowing
    ! ===============================
    sp%d(dir)%Cw=window_factor(win,sp%d(dir)%lbloc)
    if (iproc==0) print *,'     (window factor:',sp%d(dir)%Cw,')'

    ! Apply windowing
    ! ===============
    select case (dimw)
    case(1)
       do n=1,sp%d(dir)%lbloc
          var(n,:,:)=var(n,:,:)*win(n)
       enddo
    case(2)
       do n=1,sp%d(dir)%lbloc
          var(:,n,:)=var(:,n,:)*win(n)
       enddo
    case(3)
       do n=1,sp%d(dir)%lbloc
          var(:,:,n)=var(:,:,n)*win(n)
       enddo
    end select
    
    deallocate(win)
   
  end subroutine windowing_c
  
  !============================================================================
  function tukeywin(nt,r)
  !============================================================================
    !> Definition of Tukey's window (tapered cosine window)
    !> -> nt: signal length
    !> -> r: coeff (r=1-> idem Hann window ~~> r=0 constant window)
  !============================================================================
    use mod_constant ! <- for pi and precision
    implicit none
    ! -------------------------------------------------------------------------
    integer, intent(in) :: nt
    real(wp), intent(in) :: r
    real(wp), dimension(nt) :: tukeywin
    ! -------------------------------------------------------------------------
    integer :: n,tl,th
    real(wp) :: pat,per
    real(wp), dimension(nt) :: t
    ! -------------------------------------------------------------------------

    ! r=0.   !--> rectangular window
    if (r==0.0_wp) then
       tukeywin=1.0_wp
    else
    ! r=1.   !--> Hann's window
    ! r=0.05 !--> 95% constant
       pat=1.0_wp/dble(nt-1)
       do n=1,nt
          t(n)=dble(n-1)*pat
       enddo

       per = r/2.0_wp
       tl = floor(per*(nt-1))+1
       th = nt-tl+1

       do n=1,tl
          tukeywin(n)=((1.0_wp+cos(pi/per*(t(n)-per)))/2.0_wp)
       enddo
       do n=tl+1,th-1
          tukeywin(n)=1.0_wp
       enddo
       do n=th,nt
          tukeywin(n)=((1.0_wp+cos(pi/per*(t(n)-1+per)))/2.0_wp)
       enddo
    endif

  end function tukeywin

  !============================================================================
  function kaiserbwin(nt,alpha)
  !============================================================================
    !> Definition of Kaiser-Bessel's window
    !> ** use AMOS library for Bessel functions **
  !============================================================================
    use mod_constant ! <- for pi and precision
    implicit none
    ! -------------------------------------------------------------------------
    integer, intent(in) :: nt
    real(wp), intent(in) :: alpha
    real(wp), dimension(nt) :: kaiserbwin
    ! -------------------------------------------------------------------------
    integer :: n,ifail,nzero
    real(wp) :: fnu,Zr,Zi,cyr,cyi
    complex(wp) :: num,denom
    complex(wp), parameter :: ii=(0.0_wp,1.0_wp)
    ! -------------------------------------------------------------------------

    fnu=0.0_wp
    ifail=0
    !alpha=3.0_wp

    Zr=pi*alpha
    Zi=0.0_wp
    call ZBESI(Zr,Zi,fnu,1,1,cyr,cyi,nzero,ifail)
    denom=cyr+ii*cyi

    do n=1,nt
       Zr=pi*alpha*sqrt(1.0_wp-(2.0_wp*real(n)/real(nt)-1)**2)
       Zi=0.0_wp
       call ZBESI(Zr,Zi,fnu,1,1,cyr,cyi,nzero,ifail)
       num=cyr+ii*cyi
       kaiserbwin(n)=num/denom
    enddo

  end function kaiserbwin

  !============================================================================
  function window_factor(w,nt)
  !============================================================================
    !> Compute window factor
  !============================================================================
    use precision
    implicit none
    ! -------------------------------------------------------------------------
    integer, intent(in) :: nt
    real(wp), dimension(nt), intent(in) :: w
    real(wp) :: window_factor
    ! -------------------------------------------------------------------------
    integer :: n
    real(wp) :: sum
    ! -------------------------------------------------------------------------

    sum=0.0_wp
    do n=1,nt
       sum=sum+w(n)**2
    enddo

    window_factor=sum/nt

  end function window_factor

end module mod_pp_sp_wind
