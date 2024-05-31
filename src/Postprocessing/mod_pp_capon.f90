!==============================================================================
module mod_pp_capon
!==============================================================================
  !> Module for Capon's spectral estimation
  !> also called MSVE: Modified Variance Spectral Estimation
!==============================================================================
  use precision
  implicit none

contains

  !============================================================================
  subroutine mvse_c(x,nt,ip,mode,iflag)
  !============================================================================
    !> MVSE : modified variance spectral estimation
    !> Basic algorithm [Cholesky fact. from LAPACK] - Complex version
  !============================================================================
    use mod_constant
    implicit none
    ! -------------------------------------------------------------------------
    integer, intent(in) :: ip,nt,mode
    integer, intent(out) :: iflag
    complex(wp), dimension(nt), intent(inout) :: x
    ! -------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: f
    complex(wp) :: sum
    complex(wp), dimension(ip) :: e
    complex(wp), dimension(ip,ip) :: rinv,rmat
    ! -------------------------------------------------------------------------

    ! Compute autocorrelation matrix
    ! ==============================
    call corrmat(x,nt,ip,rmat)

    ! Inversion of autocorrelation matrix
    ! ===================================
    ! ~> using Cholesky factorization from LAPACK
    rinv=0.0_wp
    do i=1,ip
       rinv(i,i)=1.0_wp
    enddo
    call ZPOSV('L',ip,ip,rmat,ip,rinv,ip,iflag)

    ! Computation of Hermitian form
    ! =============================
    do k=1,nt
       f=-0.5_wp+dble(k-1)/dble(nt)
       do j=1,ip
          e(j)=exp(cmplx(0.0_wp,twopi*(j-1)*f))
       enddo
       sum=0.0_wp
       do i=1,ip
          do j=1,ip
             sum=sum+conjg(e(i))*rinv(i,j)*e(j)
          enddo
       enddo

       ! compute PSD values (with normalization)
       ! ------------------
       if (mode==0) x(k)=real(1./sum)
       if (mode==1) x(k)=real(ip/sum)

    enddo

  end subroutine mvse_c

  !============================================================================
  subroutine corrmat(x,nt,ip,rmat)
  !============================================================================
    !> Determination of the autocorrelation matrix
    !> Modified covariance method - Complex version
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer, intent(in) :: nt,ip
    complex(wp), dimension(nt), intent(in) :: x
    complex(wp), dimension(ip,ip), intent(out) :: rmat
    ! -------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: fact
    ! -------------------------------------------------------------------------

    ! Compute autocorrelation matrix estimate using modified covariance method
    ! ========================================================================
    fact=2.0_wp*dble(nt-ip)
    
    do j=1,ip
       do i=1,j
          rmat(i,j)=0.0_wp
          do k=ip+1,nt
             rmat(i,j)=rmat(i,j)+conjg(x(k-i))*x(k-j)+x(k-ip+i)*conjg(x(k-ip+j))
          enddo
          rmat(i,j)=rmat(i,j)/fact
       enddo
    enddo
    
    ! Impose Hermitian symmetry
    ! =========================
    do i=1,ip
       do j=1,i
          rmat(i,j)=conjg(rmat(j,i))
       enddo
    enddo

  end subroutine corrmat
  
  !============================================================================
  subroutine fast_mvse_c(x,nt,ip,eps)
  !============================================================================
    !> MVSE : modified variance spectral estimation
    !> Fast algorithm of Musicus - Complex version
  !============================================================================
    use mod_pp_dfft
    implicit none
    ! -------------------------------------------------------------------------
    integer, intent(in) :: ip,nt
    real(wp), intent(in) :: eps
    complex(wp), dimension(nt), intent(inout) :: x
    ! -------------------------------------------------------------------------
    integer :: N_p,i,k,n,m
    complex(wp) :: var,sum
    complex(wp), dimension(ip) :: r,a
    complex(wp), dimension(2*ip-1) :: rn,mu
    complex(wp), dimension(nt) :: mul
    complex(wp), dimension(ip,ip) :: aa
    ! -------------------------------------------------------------------------

    ! Compute samples of unbiased autocorrelation estimator
    ! =====================================================
    call acorr(x,nt,ip,r)

    ! Construct autocorrelation sequence
    ! ==================================
    N_p=ip-1
    do k=-N_p,-1
       rn(k+ip)=conjg(r(-k+1))
    enddo
    do k=0,N_p
       rn(k+ip)=r(k+1)
    enddo

    ! Tikhonov's regularization
    ! =========================
    rn(ip)=rn(ip)*(1.+eps)

    ! Levinson-Durbin algorithm to compute partial correlation
    ! ========================================================
    
    ! initialization of Levinson recursion
    ! ------------------------------------
    aa(1,1)=1.0_wp
    var=rn(ip)

    ! Levinson-Durbin recursion
    ! -------------------------
    do n=1,N_p
       sum=0.0_wp
       do i=0,n-1
          sum=sum+aa(i+1,n)*rn(i-n+ip)
       enddo
       aa(n+1,n+1)=-sum/var
       aa(1,n+1)=1.0_wp
       do i=1,n-1
          aa(i+1,n+1)=aa(i+1,n)+aa(n+1,n+1)*conjg(aa(n-i+1,n))
       enddo
       var=(1.0_wp-abs(aa(n+1,n+1))**2)*var
    enddo
    
    ! store final solution of Yule-Walker equations
    ! ---------------------------------------------
    do i=0,N_p
       a(i+1)=aa(i+1,ip)
    enddo

    ! Compute mu (Musicus method)
    ! ==========
    do k=0,N_p 
       mu(k+ip)=0.0_wp
       do m=0,N_p-k
          mu(k+ip)=mu(k+ip)+(N_p-k-2*m+1)*a(m+1)*conjg(a(m+k+1))
       enddo
       mu(k+ip)=mu(k+ip)/var
    enddo
    do k=-N_p,-1 
       mu(k+ip)=conjg(mu(-k+ip))
    enddo

    do k=-N_p,N_p
       mu(k+ip)=mu(k+ip)*(-1)**k
    enddo
    
    ! zero-padding
    ! ------------
    mul=0.0_wp
    mul(1:2*ip-1)=mu
    
    ! Fourier transform
    ! -----------------
    call ZFFT1(mul)

    ! Compute PSD values
    ! ==================
    x=real(ip)/mul

  end subroutine fast_mvse_c

  !============================================================================
  subroutine acorr(x,n,ip,r)
  !============================================================================
    !> Determination of covariance matrix
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer, intent(in) :: n,ip
    complex(wp), dimension(n), intent(in) :: x
    complex(wp), dimension(ip), intent(out) :: r
    ! -------------------------------------------------------------------------
    integer :: j,k,nk,mode
    complex(wp) :: sum
    ! -------------------------------------------------------------------------

    mode=1 ! <- unbiased version

    do k=0,ip-1

       nk=n-k
       sum=0.0_wp
       do j=1,nk
          sum=sum+conjg(x(j))*x(j+k)
       enddo

       if (mode==0) r(k+1)=sum/dble(nk)
       if (mode==1) r(k+1)=sum/dble(n)

    enddo

  end subroutine acorr
  
  !============================================================================
  subroutine fast_mvse_r(x,nt,ip,eps)
  !============================================================================
    !> MVSE : modified variance spectral estimation
    !> Fast algorithm of Musicus - Real version 
  !============================================================================
    use mod_pp_dfft
    implicit none
    ! -------------------------------------------------------------------------
    integer, intent(in) :: ip,nt
    real(wp), intent(in) :: eps
    real(wp), dimension(nt), intent(inout) :: x
    ! -------------------------------------------------------------------------
    integer :: N_p,i,k,n,m,nk
    real(wp) :: var,sum
    real(wp), dimension(ip) :: a
    real(wp), dimension(2*ip-1) :: rn,mu
    real(wp), dimension(nt) :: mul
    real(wp), dimension(ip,ip) :: aa
    real(wp), dimension(2*nt+15) :: wsave
    ! -------------------------------------------------------------------------

    ! Compute samples of unbiased autocorrelation estimator
    ! =====================================================
    do k=0,ip-1
       nk=nt-k
       sum=0.0_wp
       do i=1,nk
          sum=sum+x(i)*x(i+k)
       enddo
       rn(k+ip)=sum/dble(nt)
    enddo
    ! impose symmetry
    N_p=ip-1
    do k=-N_p,-1
       rn(k+ip)=rn(ip-k)
    enddo

    ! Tikhonov's regularization
    ! =========================
    rn(ip)=rn(ip)*(1.+eps)

    ! Levinson-Durbin algorithm to compute partial correlation
    ! ========================================================
    
    ! initialization of Levinson recursion
    ! ------------------------------------
    aa(1,1)=1.0_wp
    var=rn(ip)

    if (var==0) var=1

    ! Levinson-Durbin recursion
    ! -------------------------
    do n=1,N_p
       sum=0.0_wp
       do i=0,n-1
          sum=sum+aa(i+1,n)*rn(i-n+ip)
       enddo
       aa(n+1,n+1)=-sum/var
       aa(1,n+1)=1.0_wp
       do i=1,n-1
          aa(i+1,n+1)=aa(i+1,n)+aa(n+1,n+1)*aa(n-i+1,n)
       enddo
       var=(1.0_wp-abs(aa(n+1,n+1))**2)*var
    enddo
    
    ! store final solution of Yule-Walker equations
    ! ---------------------------------------------
    do i=0,N_p
       a(i+1)=aa(i+1,ip)
    enddo

    ! Compute mu (Musicus method)
    ! ==========
    do k=0,N_p 
       mu(k+ip)=0.0_wp
       do m=0,N_p-k
          mu(k+ip)=mu(k+ip)+(N_p-k-2*m+1)*a(m+1)*a(m+k+1)
       enddo
       mu(k+ip)=mu(k+ip)/var
    enddo
    do k=-N_p,-1 
       mu(k+ip)=mu(ip-k)
    enddo

    do k=-N_p,N_p
       mu(k+ip)=mu(k+ip)*(-1.0_wp)**k
    enddo

    ! zero-padding
    ! ------------
    mul=0.0_wp
    mul(1:2*ip-1)=mu

    ! Fourier transform
    ! -----------------
    call dffti(nt,wsave)
    call dfftf(nt,mul,wsave)   

    ! Compute PSD values
    ! ==================
    !check
    do k=2,nt/2
       !if (mul(2*k-2)==0) print *,'k,2*k-2',k,2*k-2
       !if (mul(2*k-1)==0) print *,'k,2*k-1',k,2*k-1
       if ((mul(2*k-2)==0).and.(mul(2*k-1)==0)) then
          print *,'k,2*k-2,2*k-1',k,2*k-2,2*k-1
          mul(2*k-2)=1.0_wp
          mul(2*k-1)=1.0_wp
       endif
    enddo
    
    x(nt/2)=dble(ip)/mul(1)
    do k=2,nt/2
       nk=nt/2-k+1
       x(nk)=dble(ip)/sqrt(mul(2*k-2)**2+mul(2*k-1)**2) 
    enddo

  end subroutine fast_mvse_r
  !===============================================================

end module mod_pp_capon
