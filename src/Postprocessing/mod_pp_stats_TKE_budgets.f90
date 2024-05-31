!==============================================================================
module mod_pp_stats_TKE_budgets
!==============================================================================
  !> Module to compute TKE budgets from stats files
!==============================================================================
  use mod_grid
  use mod_pp_stats_var
  implicit none
  ! ---------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------

contains

  !============================================================================
  subroutine pp_TKE_budgets_chan
  !============================================================================
    !> Compute 1-D TKE budgets for channel flow (homogeneous in x and z)
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    ! -------------------------------------------------------------------------

    ! not included yet
    
  end subroutine pp_TKE_budgets_chan

  !============================================================================
  subroutine pp_TKE_budgets_xy
  !============================================================================
    !> Compute 2-D TKE budgets for xy-stats (homogeneous in z)
  !============================================================================
    !
    ! VARIABLES (from stats) USED for 2D BUDGETS:
    ! -------------------------------------------
    ! rhom,um,vm,wm,pm,rhoum,rhovm,rhowm
    ! tau11m,tau12m,tau13m,tau22m,tau23m,tau33m
    ! rhouum,rhovvm,rhowwm,rhouvm,rhouwm,rhovwm
    ! rhouuum,rhouuv,rhouvvm,rhovvvm,rhowwum,rhowwv
    ! tau11um,tau12um,tau12vm,tau22vm,tau13wm,tau23w
    ! tau11duxm,tau12duym,tau13duzm
    ! tau12dvxm,tau22dvym,tau23dvzm
    ! tau13dwxm,tau23dwym,tau33dwzm
    ! pduxm,pdvym,pdwzm,upm,vpm
    !
    ! VARIABLES (calculated) USED for 2D BUDGETS:
    ! -------------------------------------------
    ! duxm,duym,duzm,dvxm,dvym,dvzm,dvxm,dvym,dvzm
    ! umfavre,vmfavre,wmfavre
    ! uffavre,vffavre,wffavre
  !============================================================================
    use mod_deriv2d      ! <- for derivatives
    use mod_utils        ! <- for numchar
    use warnstop
    use mod_constant
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i,j,n,ii,ime
    real(wp), dimension(ngx,ngy) :: wrk1,wrk2
    real(wp), dimension(ngx) :: normb
    real(wp), dimension(ngy) :: dy0
    real(wp), dimension(ngx,ngy,8) :: budget_uu,budget_vv,budget_ww
    ! averaging on a slice
    real(wp), dimension(ngy,8) :: budget_uum,budget_vvm,budget_wwm
    ! -------------------------------------------------------------------------

    ! ! for DNS
    ! ! =======
    ! ! average on a slice
    ! ime=100
    ! ! loc for Re_theta=4060
    ! !ii=6943
    ! ii=7185
    ! ! loc for Re_theta=3270
    ! !ii=5720

!!$    ! for LES
!!$    ! =======
!!$    ! average on a slice
!!$    ime=10
!!$    ! loc for Re_theta=3270
!!$    ii=2945

    ! ! For LES flat plates FSTT
    ! ime=10
    ! ii=10
    ! do while ((Re_theta(ii).lt.2000).and.(ii.lt.ngx-50))
    !    ii=ii+1
    ! enddo
    ! print *,'  ~> Budget realized at Re_theta',Re_theta(ii)

    ! For LE configuration
    ime=10
    ii=ngx/2

    ! Balance of <ρ u''u''>/2
    ! =======================
    
    ! Mean flow convection C11 = -d/dxj(1/2<ρ.u''u''>.[uj])
    ! -----------------------------------------------------
    ! <ρ.u''u''>
    wrk1=rhouum-rhoum*rhoum/rhom
    ! d/dx(1/2<ρ.u''u''>.[u])
    call deriv2_x_11pts(wrk1*umfavre,wrk2)
    ! update C11
    budget_uu(:,:,1)=-0.5_wp*wrk2
    ! d/dy(1/2<ρ.u''u''>.[v])
    call deriv2_y_11pts(wrk1*vmfavre,wrk2)
    ! update C11
    budget_uu(:,:,1)=budget_uu(:,:,1)-0.5_wp*wrk2

    ! Production rate P11 = -<ρ.u''uj''>.d[u]/dxj
    ! -------------------------------------------
    ! d[u]/dx
    call deriv2_x_11pts(umfavre,wrk1)
    ! <ρ.u''u''>
    wrk2=rhouum-rhoum*rhoum/rhom
    ! update P11
    budget_uu(:,:,2)=-wrk2*wrk1
    ! d[u]/dy
    call deriv2_y_11pts(umfavre,wrk1)
    ! <ρ.u''v''>
    wrk2=rhouvm-rhoum*rhovm/rhom
    ! update P11
    budget_uu(:,:,2)=budget_uu(:,:,2)-wrk2*wrk1

    ! Turbulent diffusion/transport TD11 = -1/2*d/dxj(<ρ.u''u''uj''>)
    ! -----------------------------------------------------
    ! <ρ.u''u''u''>
    wrk1=rhouuum-rhouum*umfavre                &
                -2.0_wp*rhouum*umfavre         &
                +2.0_wp*rhoum*umfavre*umfavre  &
                +rhoum*umfavre**2              &
                -rhom*umfavre**2*umfavre
    ! d/dx
    call deriv2_x_11pts(wrk1,wrk2)
    ! update TD11
    budget_uu(:,:,3)=-0.5_wp*wrk2
    ! <ρ.u''u''v''>
    wrk1=rhouuvm-rhouum*vmfavre                &
                -2.0_wp*rhouvm*umfavre         &
                +2.0_wp*rhoum*umfavre*vmfavre  &
                +rhovm*umfavre**2              &
                -rhom*umfavre**2*vmfavre
    ! d/dy
    call deriv2_y_11pts(wrk1,wrk2)
    ! update TD11
    budget_uu(:,:,3)=budget_uu(:,:,3)-0.5_wp*wrk2

    ! Viscous/molecular diffusion VD11 = d/dxj(<tau1j'u'>)
    ! ------------------------------------------
    ! <tau11'u'>
    wrk1=tau11um-tau11m*um
    ! d/dx
    call deriv2_x_11pts(wrk1,wrk2)
    ! update VD11
    budget_uu(:,:,4)=wrk2
    ! <tau12'u'>
    wrk1=tau12um-tau12m*um
    ! d/dy
    call deriv2_y_11pts(wrk1,wrk2)
    ! update VD11
    budget_uu(:,:,4)=budget_uu(:,:,4)+wrk2

    ! Mass flux variation M11 = <u''>.(-d/dx(<p> + d/dxj(<tau1j>))
    ! ------------------------------------------------------------
    ! d<p>/dx
    call deriv2_x_11pts(pm,wrk1)
    ! add
    wrk2=wrk1
    ! d<tau11>/dx
    call deriv2_x_11pts(tau11m,wrk1)
    ! add
    wrk2=wrk2+wrk1
    ! d<tau12>/dy
    call deriv2_y_11pts(tau12m,wrk1)
    ! add
    wrk2=wrk2+wrk1
    ! update M11
    budget_uu(:,:,5)=uffavre*wrk2

    ! Pressure-dilation strain PS11 = <p'.du'/dx>
    ! --------------------------------------------
    budget_uu(:,:,6)=pduxm-pm*duxm

    ! Pressure diffusion PD11 = -d/dx<p'.u'>
    ! --------------------------------------
    ! <p'.u'>
    wrk1=upm-um*pm
    ! d/dx
    call deriv2_x_11pts(wrk1,wrk2)
    budget_uu(:,:,7)=-wrk2

    ! Viscous dissipation DS11 = <tau1j'.d/dxj(u')>
    ! ----------------------------------------------
    wrk2=tau11duxm+tau12duym+tau13duzm-tau11m*duxm-tau12m*duym-tau13m*duzm
    budget_uu(:,:,8)=-wrk2
    dissip = 0.0_wp
    dissip(:,:) = dissip(:,:)-wrk2

    ! Add total to dissipation (different definition compress/incompress)
    do i=1,ngx
       do j=1,ngy
          budget_uu(i,j,8)=budget_uu(i,j,8)-sum(budget_uu(i,j,1:8))
       enddo
    enddo

    ! Name
    ! Schlatter: bud_4060_dns_uu.prof
    ! Jimenez: Re_theta.4000.bal.uu
    
    ! normalization
    normb=rhowall*u_tau**4/nuwall
    if (LE) normb=1.0_wp
    do j=1,ngy
       do n=1,8
          budget_uu(:,j,n)=budget_uu(:,j,n)/normb
       enddo
    enddo

    ! average
    budget_uum=0.0_wp
    do i=ii-ime,ii+ime
       do n=1,8
          budget_uum(:,n)=budget_uum(:,n)+budget_uu(i,:,n)
       enddo
    enddo
    budget_uum=budget_uum/(2*ime+1)

    ! open(194,file='budget_Re_'//trim(numchar(int(Re_theta(ii))))//'_uu.bin',form='unformatted',status='unknown')
    ! open(194,file='budget_Re_4060_dns_uu.bin',form='unformatted',status='unknown')
    !open(194,file='budget_Re_3270_dns_uu.bin',form='unformatted',status='unknown')
    !open(194,file='budget_Re_3270_LES_IRS_uu.bin',form='unformatted',status='unknown')
    !open(194,file='budget_Re_3270_LESn_uu.bin',form='unformatted',status='unknown')
    open(194,file='budget_i'//trim(numchar(ii))//'_bl'//trim(numchar(iblc_pp))//'_uu.bin',form='unformatted',status='unknown')
    rewind(194)
    write(194) ngy
    write(194) ((budget_uum(j,n),j=1,ngy),n=1,8)
    close(194)

    ! Balance of <ρ.v''v''>/2
    ! =======================
    
    ! Mean flow convection C22 = -d/dxj(1/2<ρ.v''v''>.[uj])
    ! -----------------------------------------------------
    ! <ρ.v''v''>
    wrk1=rhovvm-rhovm*rhovm/rhom
    ! d/dx(1/2<ρ.v''v''>.[u])
    call deriv2_x_11pts(wrk1*umfavre,wrk2)
    ! update C22
    budget_vv(:,:,1)=-0.5_wp*wrk2
    ! d/dy(1/2<ρ.v''v''>.[v])
    call deriv2_y_11pts(wrk1*vmfavre,wrk2)
    ! update C22
    budget_vv(:,:,1)=budget_vv(:,:,1)-0.5_wp*wrk2

    ! Production rate P22 = -<ρ.v''uj''>.d[v]/dxj
    ! -------------------------------------------
    ! d[v]/dx
    call deriv2_x_11pts(vmfavre,wrk1)
    ! <ρ.v''u''>
    wrk2=rhouvm-rhovm*rhoum/rhom
    ! update P22
    budget_vv(:,:,2)=-wrk2*wrk1
    ! d[v]/dy
    call deriv2_y_11pts(vmfavre,wrk1)
    ! <ρ.v''v''>
    wrk2=rhovvm-rhovm*rhovm/rhom
    ! update P22
    budget_vv(:,:,2)=budget_vv(:,:,2)-wrk2*wrk1

    ! Turbulent diffusion TD22 = -1/2*d/dxj(<ρ.u''u''uj''>)
    ! -----------------------------------------------------
    ! <ρ.v''v''u''>
    wrk1=rhouvvm-rhovvm*umfavre                &
                -2.0_wp*rhouvm*vmfavre         &
                +2.0_wp*rhovm*vmfavre*umfavre  &
                +rhoum*vmfavre**2              &
                -rhom*vmfavre**2*umfavre
    ! d/dx
    call deriv2_x_11pts(wrk1,wrk2)
    ! update TD22
    budget_vv(:,:,3)=-0.5_wp*wrk2 
    ! <ρ.v''v''v''>
    !!wrk1=rhovvvm-3.0_wp*rhovvm*vmfavre+3.0_wp*rhovm*vmfavre**2-rhom*vmfavre**3
    wrk1=rhovvvm-rhovvm*vmfavre                &
                -2.0_wp*rhovvm*vmfavre         &
                +2.0_wp*rhovm*vmfavre*vmfavre  &
                +rhovm*vmfavre**2              &
                -rhom*vmfavre**2*vmfavre
    ! d/dy
    call deriv2_y_11pts(wrk1,wrk2)
    ! update TD22
    budget_vv(:,:,3)=budget_vv(:,:,3)-0.5_wp*wrk2

    ! Viscous diffusion VD22 = d/dxj(<tau2j'v'>)
    ! ------------------------------------------
    ! <tau21'v'>
    wrk1=tau12vm-tau12m*vm
    ! d/dx
    call deriv2_x_11pts(wrk1,wrk2)
    ! update VD22
    budget_vv(:,:,4)=wrk2
    ! <tau22'v'>
    wrk1=tau22vm-tau22m*vm
    ! d/dy
    call deriv2_y_11pts(wrk1,wrk2)
    ! update VD22
    budget_vv(:,:,4)=budget_vv(:,:,4)+wrk2

    ! Mass flux variation M22 = <v''>.(-d/dy(<p> + d/dxj(<tau1j>))
    ! ------------------------------------------------------------
    ! d<p>/dy
    call deriv2_y_11pts(pm,wrk1)
    ! add
    wrk2=wrk1
    ! d<tau21>/dx
    call deriv2_x_11pts(tau12m,wrk1)
    ! add
    wrk2=wrk2+wrk1
    ! d<tau22>/dy
    call deriv2_y_11pts(tau22m,wrk1)
    ! add
    wrk2=wrk2+wrk1
    ! update M22
    budget_vv(:,:,5)=vffavre*wrk2

    ! Pressure-dilation strain PS22 = <p'.dv'/dy>
    ! --------------------------------------------
    budget_vv(:,:,6)=pdvym-pm*dvym

    ! Pressure diffusion PD22 = -d/dy<p'.v'>
    ! --------------------------------------
    ! <p'.v'>
    wrk1=vpm-vm*pm
    ! d/dy
    call deriv2_y_11pts(wrk1,wrk2)
    budget_vv(:,:,7)=-wrk2

    ! Viscous dissipation DS22 = <tau2j'.d/dxj(v')>
    ! ---------------------------------------------
    wrk2=tau12dvxm+tau22dvym+tau23dvzm-tau12m*dvxm-tau22m*dvym-tau23m*dvzm
    budget_vv(:,:,8)=-wrk2
    dissip(:,:) = dissip(:,:)-wrk2

    !----------------------------------------------------------------------------
    ! (!!TEST!!)
    ! Add total to dissipation (different definition compress/incompress)
    do i=1,ngx
       do j=1,ngy
          budget_vv(i,j,8)=budget_vv(i,j,8)-sum(budget_vv(i,j,1:8))
       enddo
    enddo

    ! normalization
    do j=1,ngy
       do n=1,8
          budget_vv(:,j,n)=budget_vv(:,j,n)/normb
       enddo
    enddo

    ! average
    budget_vvm=0.0_wp
    do i=ii-ime,ii+ime
       do n=1,8
          budget_vvm(:,n)=budget_vvm(:,n)+budget_vv(i,:,n)
       enddo
    enddo
    budget_vvm=budget_vvm/(2*ime+1)

    ! open(194,file='budget_Re_'//trim(numchar(int(Re_theta(ii))))//'_vv.bin',form='unformatted',status='unknown')
    ! open(194,file='budget_Re_4060_dns_vv.bin',form='unformatted',status='unknown')
    !open(194,file='budget_Re_3270_dns_vv.bin',form='unformatted',status='unknown')
    !open(194,file='budget_Re_3270_LES_IRS_vv.bin',form='unformatted',status='unknown')
    !open(194,file='budget_Re_3270_LESn_vv.bin',form='unformatted',status='unknown')
    open(194,file='budget_i'//trim(numchar(ii))//'_bl'//trim(numchar(iblc_pp))//'_vv.bin',form='unformatted',status='unknown')
    rewind(194)
    write(194) ngy
    write(194) ((budget_vvm(j,n),j=1,ngy),n=1,8)
    close(194)

    ! Balance of <ρ.w''w''>/2
    ! =======================
    
    ! Mean flow convection C33 = -d/dxj(1/2<ρ.w''w''>.[uj])
    ! -----------------------------------------------------
    ! <ρ.w''w''>
    wrk1=rhowwm-rhowm*rhowm/rhom
    ! d/dx(1/2<ρ.w''w''>.[u])
    call deriv2_x_11pts(wrk1*umfavre,wrk2)
    ! update C33
    budget_ww(:,:,1)=-0.5_wp*wrk2
    ! d/dy(1/2<ρ.w''w''>.[v])
    call deriv2_y_11pts(wrk1*vmfavre,wrk2)
    ! update C33
    budget_ww(:,:,1)=budget_ww(:,:,1)-0.5_wp*wrk2

    ! Production rate P33 = -<ρ.w''uj''>.d[w]/dxj
    ! -------------------------------------------
    ! d[w]/dx
    call deriv2_x_11pts(wmfavre,wrk1)
    ! <ρ.w''u''>
    wrk2=rhouwm-rhowm*rhoum/rhom
    ! update P33
    budget_ww(:,:,2)=-wrk2*wrk1
    ! d[w]/dy
    call deriv2_y_11pts(wmfavre,wrk1)
    ! <ρ.w''v''>
    wrk2=rhovwm-rhowm*rhovm/rhom
    ! update P33
    budget_ww(:,:,2)=budget_ww(:,:,2)-wrk2*wrk1

    ! Turbulent diffusion TD33 = -1/2*d/dxj(<ρ.w''w''uj''>)
    ! -----------------------------------------------------
    ! <ρ.w''w''u''>
    wrk1=rhowwum-rhowwm*umfavre                &
                -2.0_wp*rhouwm*wmfavre         &
                +2.0_wp*rhowm*wmfavre*umfavre  &
                +rhoum*wmfavre**2              &
                -rhom*wmfavre**2*umfavre
    ! d/dx
    call deriv2_x_11pts(wrk1,wrk2)
    ! update TD33
    budget_ww(:,:,3)=-0.5_wp*wrk2
    ! <ρ.w''w''v''>
    wrk1=rhowwvm-rhowwm*vmfavre                &
                -2.0_wp*rhovwm*wmfavre         &
                +2.0_wp*rhowm*wmfavre*vmfavre  &
                +rhovm*wmfavre**2              &
                -rhom*wmfavre**2*vmfavre
    ! d/dy
    call deriv2_y_11pts(wrk1,wrk2)
    ! update TD33
    budget_ww(:,:,3)=budget_ww(:,:,3)-0.5_wp*wrk2

    ! Viscous diffusion VD33 = d/dxj(<tau3j'w'>)
    ! ------------------------------------------
    ! <tau31'w'>
    wrk1=tau13wm-tau13m*wm
    ! d/dx
    call deriv2_x_11pts(wrk1,wrk2)
    ! update VD33
    budget_ww(:,:,4)=wrk2
    ! <tau32'w'>
    wrk1=tau23wm-tau23m*wm
    ! d/dy
    call deriv2_y_11pts(wrk1,wrk2)
    ! update VD33
    budget_ww(:,:,4)=budget_ww(:,:,4)+wrk2

    ! Mass flux variation M33 = <w''>.(-d/dz(<p> + d/dxj(<tau3j>))
    ! ------------------------------------------------------------
    ! d<tau31>/dx
    call deriv2_x_11pts(tau13m,wrk1)
    ! add
    wrk2=wrk1
    ! d<tau32>/dy
    call deriv2_y_11pts(tau23m,wrk1)
    ! add
    wrk2=wrk2+wrk1
    ! update M11
    budget_ww(:,:,5)=wffavre*wrk2

    ! Pressure-dilation strain PS33 = <p'.dw'/dz>
    ! --------------------------------------------
    budget_ww(:,:,6)=pdwzm-pm*dwzm

    ! Pressure diffusion PD33 = -d/dz<p'.w'>=0 (homogeneous in z)
    ! --------------------------------------
    budget_ww(:,:,7)=0.0_wp

    ! Viscous dissipation DS33 = <tau3j'.d/dxj(w')>
    ! ---------------------------------------------
    wrk2=tau13dwxm+tau23dwym+tau33dwzm-tau13m*dwxm-tau23m*dwym-tau33m*dwzm
    budget_ww(:,:,8)=-wrk2
    dissip(:,:) = dissip(:,:)-wrk2

    !----------------------------------------------------------------------------
    ! (!!TEST!!)
    ! Add total to dissipation (different definition compress/incompress)
    do i=1,ngx
       do j=1,ngy
          budget_ww(i,j,8)=budget_ww(i,j,8)-sum(budget_ww(i,j,1:8))
       enddo
    enddo

    ! normalization
    do j=1,ngy
       do n=1,8
          budget_ww(:,j,n)=budget_ww(:,j,n)/normb
       enddo
    enddo

    ! average
    budget_wwm=0.0_wp
    do i=ii-ime,ii+ime
       do n=1,8
          budget_wwm(:,n)=budget_wwm(:,n)+budget_ww(i,:,n)
       enddo
    enddo
    budget_wwm=budget_wwm/(2*ime+1)

    ! open(194,file='budget_Re_'//trim(numchar(int(Re_theta(ii))))//'_ww.bin',form='unformatted',status='unknown')
    ! open(194,file='budget_Re_4060_dns_ww.bin',form='unformatted',status='unknown')
    !open(194,file='budget_Re_3270_dns_ww.bin',form='unformatted',status='unknown')
    !open(194,file='budget_Re_3270_LES_IRS_ww.bin',form='unformatted',status='unknown')
    !open(194,file='budget_Re_3270_LESn_ww.bin',form='unformatted',status='unknown')
    open(194,file='budget_i'//trim(numchar(ii))//'_bl'//trim(numchar(iblc_pp))//'_ww.bin',form='unformatted',status='unknown')
    rewind(194)
    write(194) ngy
    write(194) ((budget_wwm(j,n),j=1,ngy),n=1,8)
    close(194)
      
    ! Averaged dissipation and lenghtscales
    ! -------------------------------------

    ! dissip already defined before
    do j=1,ngy
       do i=1,ngx
          ! dissip(i,j) = abs(budget_uu(i,j,8) + budget_vv(i,j,8) + budget_ww(i,j,8)) * normb(i)
          !dissip(i,j)=abs(budget_uu(i,j,8))*normb
          dissip(i,j)=abs(dissip(i,j))
       enddo
    enddo

    open(194,file='dissip_bl'//trim(numchar(iblc_pp))//'.bin',form='unformatted',status='unknown')
    rewind(194)
    write(194) ((dissip(i,j),i=1,ngx),j=1,ngy)
    close(194)

    do j=1,ngy
       do i=1,ngx
          if (dissip(i,j).eq.0.0_wp) print *,"Dissip nul",i,j
          if (dissip(i,j).eq.0.0_wp) dissip(i,j)=1.0_wp ! problem of calculation of dissip in corners of T&D
       enddo
    enddo

    do j=1,ngy
       do i=1,ngx
          kolmog(i,j) =((mum(i,j)/rhom(i,j))**3*rhom(i,j)/dissip(i,j))**0.25_wp
       enddo
    enddo

    do j=1,ngy-1
       dy0(j)=yg(j+1)-yg(j)
    enddo
    dy0(ngy)=dy0(ngy-1)

    open(194,file='kolmog_bl'//trim(numchar(iblc_pp))//'.bin',form='unformatted',status='unknown')
    rewind(194)
    write(194) ((kolmog(i,j),i=1,ngx),j=1,ngy)
    write(194) ((deltax/kolmog(i,j),i=1,ngx),j=1,ngy)
    write(194) ((dy0(j)/kolmog(i,j),i=1,ngx),j=1,ngy)
    write(194) ((deltaz/kolmog(i,j),i=1,ngx),j=1,ngy)
    close(194)

    ! First definition: lambda = sqrt(u_rms^2 / dudx_rms^2)
    l_tayl(:,1)=1.0_wp !Pb à la paroi sinon
    do j=2,ngy
       do i=1,ngx
          ! l_tayl(i,j) = ((u2m(i,j)-um(i,j)**2)/(dux2m(i,j) - (rhoduxm(i,j)/rhom(i,j))**2))**0.5_wp
          l_tayl(i,j) = (((u2m(i,j)-um(i,j)**2 + v2m(i,j)-vm(i,j)**2 + w2m(i,j)-wm(i,j)**2)/3.0)/(dux2m(i,j) - (rhoduxm(i,j)/rhom(i,j))**2))**0.5_wp
       enddo
    enddo

    open(194,file='taylor_bl'//trim(numchar(iblc_pp))//'.bin',form='unformatted',status='unknown')
    rewind(194)
    write(194) ((l_tayl(i,j),i=1,ngx),j=1,ngy)
    write(194) ((deltax/l_tayl(i,j),i=1,ngx),j=1,ngy)
    write(194) ((dy0(j)/l_tayl(i,j),i=1,ngx),j=1,ngy)
    write(194) ((deltaz/l_tayl(i,j),i=1,ngx),j=1,ngy)
    close(194)

    ! Second definition: lambda = sqrt(15*nu*u_rms^2 / epsilon)
    l_tayl(:,1)=1.0_wp !Pb à la paroi sinon
    do j=2,ngy
       do i=1,ngx
          ! l_tayl(i,j) = ((15.0_wp*mum(i,j)*(u2m(i,j)-um(i,j)**2))/dissip(i,j))**0.5_wp
          l_tayl(i,j) = ((15.0_wp*mum(i,j)*((u2m(i,j)-um(i,j)**2 + v2m(i,j)-vm(i,j)**2 + w2m(i,j)-wm(i,j)**2)/3.0))/dissip(i,j))**0.5_wp
       enddo
    enddo

    open(194,file='taylor2_bl'//trim(numchar(iblc_pp))//'.bin',form='unformatted',status='unknown')
    rewind(194)
    write(194) ((l_tayl(i,j),i=1,ngx),j=1,ngy)
    write(194) ((deltax/l_tayl(i,j),i=1,ngx),j=1,ngy)
    write(194) ((dy0(j)/l_tayl(i,j),i=1,ngx),j=1,ngy)
    write(194) ((deltaz/l_tayl(i,j),i=1,ngx),j=1,ngy)
    close(194)

    ! Definition of integral length scale: L_f = 0.9*u_rms^3/epsilon
    ! /!\ Provide poor results (?)
    do j=1,ngy
       do i=1,ngx
          ! L_f_scale(i,j) = 0.9_wp*((u2m(i,j)-um(i,j)**2)**1.5_wp)/dissip(i,j)
          L_f_scale(i,j) = 0.9_wp*(((u2m(i,j)-um(i,j)**2 + v2m(i,j)-vm(i,j)**2 + w2m(i,j)-wm(i,j)**2)/3.0)**1.5_wp)/dissip(i,j)
       enddo
    enddo
    if (is_bc_wall(2,1)) L_f_scale(:,1)=1.0_wp !Pb à la paroi sinon

    open(194,file='L_f_bl'//trim(numchar(iblc_pp))//'.bin',form='unformatted',status='unknown')
    rewind(194)
    write(194) ((L_f_scale(i,j),i=1,ngx),j=1,ngy)
    write(194) ((deltax/L_f_scale(i,j),i=1,ngx),j=1,ngy)
    write(194) ((dy0(j)/L_f_scale(i,j),i=1,ngx),j=1,ngy)
    write(194) ((deltaz/L_f_scale(i,j),i=1,ngx),j=1,ngy)
    close(194)

    ! Definition of Re_FST = L_f*u'_rms/nu
    do j=1,ngy
       do i=1,ngx
          Re_fst(i,j) = L_f_scale(i,j)*((u2m(i,j)-um(i,j)**2 + v2m(i,j)-vm(i,j)**2 + w2m(i,j)-wm(i,j)**2)/3.0)**0.5*rhom(i,j)/mum(i,j)
       enddo
    enddo

    open(194,file='Re_fst_bl'//trim(numchar(iblc_pp))//'.bin',form='unformatted',status='unknown')
    rewind(194)
    write(194) ((Re_fst(i,j),i=1,ngx),j=1,ngy)
    ! Definition of Re_Lf = L_f*U/nu
    do j=1,ngy
       do i=1,ngx
          Re_fst(i,j) = L_f_scale(i,j)*(um(i,j)**2 + vm(i,j)**2 + wm(i,j)**2)**0.5*rhom(i,j)/mum(i,j)
       enddo
    enddo
    write(194) ((Re_fst(i,j),i=1,ngx),j=1,ngy)
    close(194)


    
  end subroutine pp_TKE_budgets_xy

end module mod_pp_stats_TKE_budgets
