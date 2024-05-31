!==============================================================================
module mod_pp_stats_skin_friction
!==============================================================================
  !> Module to compute mean skin friction decomposition
  !> Nota: not parallel -> to be run on a single proc
  !>       only for Cartesian coordinate
!==============================================================================
  use mod_grid
  use mod_constant
  use mod_deriv2d
  use mod_pp_var
  use mod_pp_stats_var
  implicit none
  ! ---------------------------------------------------------------------------
  real(wp), dimension(:,:), allocatable :: Cfi   ! Cf decomposition integrated
  real(wp), dimension(:,:), allocatable :: Cfi_i ! Cf decomposition inner layer
  real(wp), dimension(:,:), allocatable :: Cfi_l ! Cf decomposition  log  layer
  real(wp), dimension(:,:), allocatable :: Cfi_o ! Cf decomposition outer layer
  ! ---------------------------------------------------------------------------

contains

  !============================================================================
  subroutine pp_FIK_chan
  !============================================================================
    !> FIK decomposition of mean skin friction (Fukagata, Iwamoto, Kasagi, PoF 2002)
    !> compressible version for channel flow (Gomez, Flutet, Sagaut, PRE 2009)
  !============================================================================
    !
    ! VARIABLES (from stats):
    ! -----------------------
    ! rhouvm rhoum rhovm rhom mum tau12m
    !
    ! VARIABLES (calculated):
    ! -----------------------
    ! duym ygw c_f rhobulk ubulk hc muwall Rebulkb
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: j
    real(wp) :: dyloc
    real(wp), dimension(:), allocatable :: wrk
    ! -------------------------------------------------------------------------

    ! Memory allocations
    ! ------------------
    ! decomposition before integration: size 1 x ngy
    allocate(Cfi(1,4))
    ! wrk array
    allocate(wrk(ngy))
    
    ! Fukagata-Iwamoto-Kasagi identity
    ! ================================
    ! -> Cf_laminar constant
    ! -> Cf(:,:,1): Cf_turbulent
    ! -> Cf(:,:,2): Cf_compressibility
    ! -> Cf(:,:,3): Cf_ct (compressible-turbulent interaction)
    ! -> Cf(:,:,4): Cf reconstruction (sum)

    ! Useful quantities
    ! -----------------   
    wrk = (rhouvm(1,:)-rhoum(1,:)*rhovm(1,:)/rhom(1,:))/(rhobulk*ubulk**2)

    ! integration (nondimensional)
    Cfi=0.0_wp
    do j=2,ngy/2+1
       dyloc= ygw(1,j)-ygw(1,j-1)
       Cfi(1,1)=Cfi(1,1)-(1.0_wp-ygw(1,j))*dyloc*wrk(j)
       Cfi(1,2)=Cfi(1,2)+(1.0_wp-ygw(1,j))*dyloc*(mum(1,j)/muwall(1)-1.0_wp)*duym(1,j)*hc/ubulk
       Cfi(1,3)=Cfi(1,3)+(1.0_wp-ygw(1,j))*dyloc*(tau12m(1,j)-mum(1,j)*duym(1,j))*hc/ubulk/muwall(1)
    enddo

    ! Normalization
    ! -------------
    Cfi(1,1)=Cfi(1,1)*6.0_wp
    Cfi(1,2)=Cfi(1,2)*6.0_wp/Rebulkb
    Cfi(1,3)=Cfi(1,3)*6.0_wp/Rebulkb

    ! Summation for check of Cf reconstruction
    ! ----------------------------------------
    Cfi(1,4)=6.0_wp/Rebulkb+Cfi(1,1)+Cfi(1,2)+Cfi(1,3) !+ 6.0_wp/bid

    ! Write ASCII file
    ! -----------------
    open(69, file='FIK_Cf.dat', status='replace', form='formatted')
    write(69,*) 'Clam', 6.0_wp/Rebulkb,6.0_wp/Rebulkb/Cfi(1,4)*100.0_wp ! Cf_laminar
    write(69,*) 'Ctur', Cfi(1,1), Cfi(1,1)/Cfi(1,4)*100.0_wp    ! Cf_turbulent
    write(69,*) 'Ccom', Cfi(1,2), Cfi(1,2)/Cfi(1,4)*100.0_wp    ! Cf_compressibility
    write(69,*) 'C_ct', Cfi(1,3), Cfi(1,3)/Cfi(1,4)*100.0_wp    ! Cf_ct (compressible-turbulent interaction)
    write(69,*) 'Ctot', Cfi(1,4),(Cfi(1,4)/c_f-1.0_wp)*100.0_wp ! Cf reconstruction
    close(69)

    ! free memory
    ! -----------
    deallocate(Cfi,wrk)
    
  end subroutine pp_FIK_chan

  !============================================================================
  subroutine pp_FIK_stbl
  !============================================================================
    !> FIK decomposition of mean skin friction (Fukagata, Iwamoto, Kasagi, PoF 2002)
    !> compressible version for boundary layer (Gomez, Flutet, Sagaut, PRE 2009)
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    ! -------------------------------------------------------------------------
    
  end subroutine pp_FIK_stbl

  !============================================================================
  subroutine pp_RD_stbl
  !============================================================================
    !> RD decomposition of mean skin friction (Renard, Deck, JFM 2016)
    !> compressible version for boundary layer (Fan, Li, Pirozzoli, PoF 2019)
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i,j,id
    integer, dimension(:), allocatable :: jb
    real(wp) :: u_delt,rho_delt
    real(wp) :: arg1,arg2,ra
    real(wp), dimension(:), allocatable :: fact
    real(wp), dimension(:,:), allocatable :: wrk1,wrk2
    ! derivative of Favre velocity components
    real(wp), dimension(:,:), allocatable :: duxfavre,duyfavre
    ! -------------------------------------------------------------------------

    ! Memory allocations
    ! ------------------
    ! Cf decomposition
    allocate(Cfi(ngx,4))
    allocate(Cfi_i(ngx,3),Cfi_l(ngx,3),Cfi_o(ngx,3))
    ! integration bound
    allocate(jb(ngx),fact(ngx))
    ! wrk arrays
    allocate(wrk1(ngx,ngy),wrk2(ngx,ngy))
    allocate(duxfavre(ngx,ngy),duyfavre(ngx,ngy))
    
    ! Renard-Deck decomposition
    ! =========================
    ! -> Cfi(:,1): Cf_v ('viscous' term)
    !           -> direct effect of viscous dissipation, transforming mechanical work into heat
    ! -> Cfi(:,2): Cf_t ('turbulent' term)
    !           -> power spent for turbulence kinetic energy production
    ! -> Cfi(:,3): Cf_g ('growth' term)
    !           -> accounts for spatial growth of the flow
    ! -> Cfi(:,4): Cf reconstruction (sum)

    ! Normalization factor
    ! --------------------
    id=0
    do i=1,ngx
       jb(i)=j99(i)+id
       !rho_delt=rho_ref
       !u_delt=u_ref
       rho_delt=rhom(i,jb(i))
       u_delt=um(i,jb(i))
       fact(i)=2.0_wp/rho_delt/u_delt**3
    enddo
    
    ! Useful quantities
    ! -----------------   
    call deriv2_x_11pts(umfavre,duxfavre)
    call deriv2_y_11pts(umfavre,duyfavre)

    ! <tau_xx> - <rho><u''u''>
    wrk1 = tau11m-rhouum+rhoum*rhoum/rhom
    ! d/dx(<tau_xx>-<rho><u''u''>)
    call deriv2_x_11pts(wrk1,wrk2) 

    ! <rho><u''v''>
    wrk1 = rhouvm-rhoum*rhovm/rhom
    
    ! Integration along y (with trapezoidal rule)
    ! -------------------
    Cfi=0.0_wp
    do i=1,ngx

       do j=1,jb(i)
          arg1=tau12m(i,j)*duyfavre(i,j)
          arg2=tau12m(i,j+1)*duyfavre(i,j+1)
          Cfi(i,1)=Cfi(i,1)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
          arg1=-wrk1(i,j)*duyfavre(i,j)
          arg2=-wrk1(i,j+1)*duyfavre(i,j+1)
          Cfi(i,2)=Cfi(i,2)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
          arg1=(umfavre(i,j)-u_ref)*( rhom(i,j)*umfavre(i,j)*duxfavre(i,j) &
                                     +rhom(i,j)*vmfavre(i,j)*duyfavre(i,j)-wrk2(i,j) )
          arg2=(umfavre(i,j+1)-u_ref)*( rhom(i,j+1)*umfavre(i,j+1)*duxfavre(i,j+1) &
                                       +rhom(i,j+1)*vmfavre(i,j+1)*duyfavre(i,j+1)-wrk2(i,j+1) )
          Cfi(i,3)=Cfi(i,3)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo
!!$       ra=(0.99_wp*u_ref-um(i,j99(i)-1))/(um(i,j99(i))-um(i,j99(i)-1))
!!$       
!!$       arg1=tau12m(i,j99(i)-1)*duyfavre(i,j99(i)-1)
!!$       arg2=tau12m(i,j99(i))*duyfavre(i,j99(i))
!!$       arg2=arg1+ra*(arg2-arg1)
!!$       Cfi(i,1)=Cfi(i,1)+0.5_wp*(arg1+arg2)*(delt99_(i)-yg(j99(i)-1))
!!$       
!!$       arg1=-wrk1(i,j99(i)-1)*duyfavre(i,j99(i)-1)
!!$       arg2=-wrk1(i,j99(i))*duyfavre(i,j99(i))
!!$       arg2=arg1+ra*(arg2-arg1)       
!!$       Cfi(i,2)=Cfi(i,2)+0.5_wp*(arg1+arg2)*(delt99_(i)-yg(j99(i)-1))
!!$       
!!$       arg1=(umfavre(i,j99(i)-1)-u_ref)*( rhom(i,j99(i)-1)*umfavre(i,j99(i)-1)*duxfavre(i,j99(i)-1) &
!!$                                         +rhom(i,j99(i)-1)*vmfavre(i,j99(i)-1)*duyfavre(i,j99(i)-1)-wrk2(i,j99(i)-1) )
!!$       arg2=(umfavre(i,j99(i))-u_ref)*( rhom(i,j99(i))*umfavre(i,j99(i))*duxfavre(i,j99(i)) &
!!$                                       +rhom(i,j99(i))*vmfavre(i,j99(i))*duyfavre(i,j99(i))-wrk2(i,j99(i)) )
!!$       arg2=arg1+ra*(arg2-arg1)      
!!$       Cfi(i,3)=Cfi(i,3)+0.5_wp*(arg1+arg2)*(delt99_(i)-yg(j99(i)-1))
       
       Cfi(i,1)=Cfi(i,1)*fact(i)
       Cfi(i,2)=Cfi(i,2)*fact(i)
       Cfi(i,3)=Cfi(i,3)*fact(i)
    enddo

    ! Integration along y for inner layer
    ! -----------------------------------
    Cfi_i=0.0_wp
    do i=1,ngx
       do j=1,j_inn(i)
          arg1=tau12m(i,j)*duyfavre(i,j)
          arg2=tau12m(i,j+1)*duyfavre(i,j+1)
          Cfi_i(i,1)=Cfi_i(i,1)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
          arg1=-wrk1(i,j)*duyfavre(i,j)
          arg2=-wrk1(i,j+1)*duyfavre(i,j+1)
          Cfi_i(i,2)=Cfi_i(i,2)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
          arg1=(umfavre(i,j)-u_ref)*( rhom(i,j)*umfavre(i,j)*duxfavre(i,j) &
                                     +rhom(i,j)*vmfavre(i,j)*duyfavre(i,j)-wrk2(i,j) )
          arg2=(umfavre(i,j+1)-u_ref)*( rhom(i,j+1)*umfavre(i,j+1)*duxfavre(i,j+1) &
                                       +rhom(i,j+1)*vmfavre(i,j+1)*duyfavre(i,j+1)-wrk2(i,j+1) )
          Cfi_i(i,3)=Cfi_i(i,3)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo
!!$       ra=(30*nuwall(i)/u_tau(i)-yg(j_inn(i)-1))/(yg(j_inn(i))-yg(j_inn(i)-1))
!!$             
!!$       arg1=tau12m(i,j_inn(i)-1)*duyfavre(i,j_inn(i)-1)
!!$       arg2=tau12m(i,j_inn(i))*duyfavre(i,j_inn(i))
!!$       arg2=arg1+ra*(arg2-arg1)
!!$       Cfi_i(i,1)=Cfi_i(i,1)+0.5_wp*(arg1+arg2)*(30*nuwall(i)/u_tau(i)-yg(j_inn(i)-1))
!!$       
!!$       arg1=-wrk1(i,j_inn(i)-1)*duyfavre(i,j_inn(i)-1)
!!$       arg2=-wrk1(i,j_inn(i))*duyfavre(i,j_inn(i))
!!$       arg2=arg1+ra*(arg2-arg1)       
!!$       Cfi_i(i,2)=Cfi_i(i,2)+0.5_wp*(arg1+arg2)*(30*nuwall(i)/u_tau(i)-yg(j_inn(i)-1))
!!$       
!!$       arg1=(umfavre(i,j_inn(i)-1)-u_ref)*( rhom(i,j_inn(i)-1)*umfavre(i,j_inn(i)-1)*duxfavre(i,j_inn(i)-1) &
!!$                                           +rhom(i,j_inn(i)-1)*vmfavre(i,j_inn(i)-1)*duyfavre(i,j_inn(i)-1)-wrk2(i,j_inn(i)-1) )
!!$       arg2=(umfavre(i,j_inn(i))-u_ref)*( rhom(i,j_inn(i))*umfavre(i,j_inn(i))*duxfavre(i,j_inn(i)) &
!!$                                         +rhom(i,j_inn(i))*vmfavre(i,j_inn(i))*duyfavre(i,j_inn(i))-wrk2(i,j_inn(i)) )
!!$       arg2=arg1+ra*(arg2-arg1)      
!!$       Cfi_i(i,3)=Cfi_i(i,3)+0.5_wp*(arg1+arg2)*(30*nuwall(i)/u_tau(i)-yg(j_inn(i)-1))
       
       Cfi_i(i,1)=Cfi_i(i,1)*fact(i)
       Cfi_i(i,2)=Cfi_i(i,2)*fact(i)
       Cfi_i(i,3)=Cfi_i(i,3)*fact(i)
    enddo

    ! Integration along y for log layer
    ! ---------------------------------
    Cfi_l=0.0_wp
    do i=1,ngx
       do j=j_inn(i)+1,j_log(i)
          arg1=tau12m(i,j)*duyfavre(i,j)
          arg2=tau12m(i,j+1)*duyfavre(i,j+1)
          Cfi_l(i,1)=Cfi_l(i,1)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
          arg1=-wrk1(i,j)*duyfavre(i,j)
          arg2=-wrk1(i,j+1)*duyfavre(i,j+1)
          Cfi_l(i,2)=Cfi_l(i,2)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
          arg1=(umfavre(i,j)-u_ref)*( rhom(i,j)*umfavre(i,j)*duxfavre(i,j) &
                                     +rhom(i,j)*vmfavre(i,j)*duyfavre(i,j)-wrk2(i,j) )
          arg2=(umfavre(i,j+1)-u_ref)*( rhom(i,j+1)*umfavre(i,j+1)*duxfavre(i,j+1) &
                                       +rhom(i,j+1)*vmfavre(i,j+1)*duyfavre(i,j+1)-wrk2(i,j+1) )
          Cfi_l(i,3)=Cfi_l(i,3)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo
!!$       ra=(30*nuwall(i)/u_tau(i)-yg(j_inn(i)-1))/(yg(j_inn(i))-yg(j_inn(i)-1))
!!$             
!!$       arg1=tau12m(i,j_inn(i)-1)*duyfavre(i,j_inn(i)-1)
!!$       arg2=tau12m(i,j_inn(i))*duyfavre(i,j_inn(i))
!!$       arg2=arg1+ra*(arg2-arg1)
!!$       Cfi_l(i,1)=Cfi_l(i,1)+0.5_wp*(arg1+arg2)*(30*nuwall(i)/u_tau(i)-yg(j_inn(i)-1))
!!$       
!!$       arg1=-wrk1(i,j_inn(i)-1)*duyfavre(i,j_inn(i)-1)
!!$       arg2=-wrk1(i,j_inn(i))*duyfavre(i,j_inn(i))
!!$       arg2=arg1+ra*(arg2-arg1)       
!!$       Cfi_l(i,2)=Cfi_l(i,2)+0.5_wp*(arg1+arg2)*(30*nuwall(i)/u_tau(i)-yg(j_inn(i)-1))
!!$       
!!$       arg1=(umfavre(i,j_inn(i)-1)-u_ref)*( rhom(i,j_inn(i)-1)*umfavre(i,j_inn(i)-1)*duxfavre(i,j_inn(i)-1) &
!!$                                           +rhom(i,j_inn(i)-1)*vmfavre(i,j_inn(i)-1)*duyfavre(i,j_inn(i)-1)-wrk2(i,j_inn(i)-1) )
!!$       arg2=(umfavre(i,j_inn(i))-u_ref)*( rhom(i,j_inn(i))*umfavre(i,j_inn(i))*duxfavre(i,j_inn(i)) &
!!$                                         +rhom(i,j_inn(i))*vmfavre(i,j_inn(i))*duyfavre(i,j_inn(i))-wrk2(i,j_inn(i)) )
!!$       arg2=arg1+ra*(arg2-arg1)      
!!$       Cfi_l(i,3)=Cfi_l(i,3)+0.5_wp*(arg1+arg2)*(30*nuwall(i)/u_tau(i)-yg(j_inn(i)-1))
!!$       
!!$       ra=(0.3*delt99_(i)-yg(j_log(i)-1))/(yg(j_log(i))-yg(j_log(i)-1))
!!$             
!!$       arg1=tau12m(i,j_log(i)-1)*duyfavre(i,j_log(i)-1)
!!$       arg2=tau12m(i,j_log(i))*duyfavre(i,j_log(i))
!!$       arg2=arg1+ra*(arg2-arg1)
!!$       Cfi_l(i,1)=Cfi_l(i,1)+0.5_wp*(arg1+arg2)*(0.3*delt99_(i)-yg(j_log(i)-1))
!!$       
!!$       arg1=-wrk1(i,j_log(i)-1)*duyfavre(i,j_log(i)-1)
!!$       arg2=-wrk1(i,j_log(i))*duyfavre(i,j_log(i))
!!$       arg2=arg1+ra*(arg2-arg1)       
!!$       Cfi_l(i,2)=Cfi_l(i,2)+0.5_wp*(arg1+arg2)*(0.3*delt99_(i)-yg(j_log(i)-1))
!!$       
!!$       arg1=(umfavre(i,j_log(i)-1)-u_ref)*( rhom(i,j_log(i)-1)*umfavre(i,j_log(i)-1)*duxfavre(i,j_log(i)-1) &
!!$                                           +rhom(i,j_log(i)-1)*vmfavre(i,j_log(i)-1)*duyfavre(i,j_log(i)-1)-wrk2(i,j_log(i)-1) )
!!$       arg2=(umfavre(i,j_log(i))-u_ref)*( rhom(i,j_log(i))*umfavre(i,j_log(i))*duxfavre(i,j_log(i)) &
!!$                                         +rhom(i,j_log(i))*vmfavre(i,j_log(i))*duyfavre(i,j_log(i))-wrk2(i,j_log(i)) )
!!$       arg2=arg1+ra*(arg2-arg1)      
!!$       Cfi_l(i,3)=Cfi_l(i,3)+0.5_wp*(arg1+arg2)*(0.3*delt99_(i)-yg(j_log(i)-1))
       
       Cfi_l(i,1)=Cfi_l(i,1)*fact(i)
       Cfi_l(i,2)=Cfi_l(i,2)*fact(i)
       Cfi_l(i,3)=Cfi_l(i,3)*fact(i)
    enddo

    ! Integration along y for outer layer
    ! -----------------------------------
    Cfi_o=0.0_wp
    do i=1,ngx
       do j=j_log(i)+1,jb(i)
          arg1=tau12m(i,j)*duyfavre(i,j)
          arg2=tau12m(i,j+1)*duyfavre(i,j+1)
          Cfi_o(i,1)=Cfi_o(i,1)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
          arg1=-wrk1(i,j)*duyfavre(i,j)
          arg2=-wrk1(i,j+1)*duyfavre(i,j+1)
          Cfi_o(i,2)=Cfi_o(i,2)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
          arg1=(umfavre(i,j)-u_ref)*( rhom(i,j)*umfavre(i,j)*duxfavre(i,j) &
                                     +rhom(i,j)*vmfavre(i,j)*duyfavre(i,j)-wrk2(i,j) )
          arg2=(umfavre(i,j+1)-u_ref)*( rhom(i,j+1)*umfavre(i,j+1)*duxfavre(i,j+1) &
                                       +rhom(i,j+1)*vmfavre(i,j+1)*duyfavre(i,j+1)-wrk2(i,j+1) )
          Cfi_o(i,3)=Cfi_o(i,3)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo
!!$       ra=(0.3*delt99_(i)-yg(j_log(i)-1))/(yg(j_log(i))-yg(j_log(i)-1))
!!$             
!!$       arg1=tau12m(i,j_log(i)-1)*duyfavre(i,j_log(i)-1)
!!$       arg2=tau12m(i,j_log(i))*duyfavre(i,j_log(i))
!!$       arg2=arg1+ra*(arg2-arg1)
!!$       Cfi_o(i,1)=Cfi_o(i,1)+0.5_wp*(arg1+arg2)*(0.3*delt99_(i)-yg(j_log(i)-1))
!!$       
!!$       arg1=-wrk1(i,j_log(i)-1)*duyfavre(i,j_log(i)-1)
!!$       arg2=-wrk1(i,j_log(i))*duyfavre(i,j_log(i))
!!$       arg2=arg1+ra*(arg2-arg1)       
!!$       Cfi_o(i,2)=Cfi_o(i,2)+0.5_wp*(arg1+arg2)*(0.3*delt99_(i)-yg(j_log(i)-1))
!!$       
!!$       arg1=(umfavre(i,j_log(i)-1)-u_ref)*( rhom(i,j_log(i)-1)*umfavre(i,j_log(i)-1)*duxfavre(i,j_log(i)-1) &
!!$                                           +rhom(i,j_log(i)-1)*vmfavre(i,j_log(i)-1)*duyfavre(i,j_log(i)-1)-wrk2(i,j_log(i)-1) )
!!$       arg2=(umfavre(i,j_log(i))-u_ref)*( rhom(i,j_log(i))*umfavre(i,j_log(i))*duxfavre(i,j_log(i)) &
!!$                                         +rhom(i,j_log(i))*vmfavre(i,j_log(i))*duyfavre(i,j_log(i))-wrk2(i,j_log(i)) )
!!$       arg2=arg1+ra*(arg2-arg1)      
!!$       Cfi_o(i,3)=Cfi_o(i,3)+0.5_wp*(arg1+arg2)*(0.3*delt99_(i)-yg(j_log(i)-1))
!!$
!!$       !ra=(0.99_wp*u_ref-um(i,j99(i)-1))/(um(i,j99(i))-um(i,j99(i)-1))
!!$       ra=(delt99_(i)-yg(j99(i)-1))/(yg(j99(i))-yg(j99(i)-1))
!!$      
!!$       arg1=tau12m(i,j99(i)-1)*duyfavre(i,j99(i)-1)
!!$       arg2=tau12m(i,j99(i))*duyfavre(i,j99(i))
!!$       arg2=arg1+ra*(arg2-arg1)
!!$       Cfi_o(i,1)=Cfi_o(i,1)+0.5_wp*(arg1+arg2)*(delt99_(i)-yg(j99(i)-1))
!!$       
!!$       arg1=-wrk1(i,j99(i)-1)*duyfavre(i,j99(i)-1)
!!$       arg2=-wrk1(i,j99(i))*duyfavre(i,j99(i))
!!$       arg2=arg1+ra*(arg2-arg1)       
!!$       Cfi_o(i,2)=Cfi_o(i,2)+0.5_wp*(arg1+arg2)*(delt99_(i)-yg(j99(i)-1))
!!$       
!!$       arg1=(umfavre(i,j99(i)-1)-u_ref)*( rhom(i,j99(i)-1)*umfavre(i,j99(i)-1)*duxfavre(i,j99(i)-1) &
!!$                                         +rhom(i,j99(i)-1)*vmfavre(i,j99(i)-1)*duyfavre(i,j99(i)-1)-wrk2(i,j99(i)-1) )
!!$       arg2=(umfavre(i,j99(i))-u_ref)*( rhom(i,j99(i))*umfavre(i,j99(i))*duxfavre(i,j99(i)) &
!!$                                       +rhom(i,j99(i))*vmfavre(i,j99(i))*duyfavre(i,j99(i))-wrk2(i,j99(i)) )
!!$       arg2=arg1+ra*(arg2-arg1)      
!!$       Cfi_o(i,3)=Cfi_o(i,3)+0.5_wp*(arg1+arg2)*(delt99_(i)-yg(j99(i)-1))
       
       Cfi_o(i,1)=Cfi_o(i,1)*fact(i)
       Cfi_o(i,2)=Cfi_o(i,2)*fact(i)
       Cfi_o(i,3)=Cfi_o(i,3)*fact(i)
    enddo

    ! Summation for check of Cf reconstruction
    ! ----------------------------------------
    Cfi(:,4)=Cfi(:,1)+Cfi(:,2)+Cfi(:,3)
    
    ! Write binary file
    ! -----------------
    open(69,file=trim(dirRESU)//'RD_Cf'//trim(name_output)//'.bin',form='unformatted',status='unknown')
    rewind(69)
    write(69) ngx
    write(69) (Cfi(i,1),i=1,ngx) ! Cf_v
    write(69) (Cfi(i,2),i=1,ngx) ! Cf_t
    write(69) (Cfi(i,3),i=1,ngx) ! Cf_g
    write(69) (Cfi(i,4),i=1,ngx) ! Cf reconstruction
    write(69) (c_f(i),i=1,ngx) ! Cf computed
    close(69)
    
    open(69,file=trim(dirRESU)//'RD_Cf'//trim(name_output)//'_layers.bin',form='unformatted',status='unknown')
    rewind(69)
    write(69) ngx
    write(69) (Cfi_i(i,1),i=1,ngx) ! Cf_v
    write(69) (Cfi_i(i,2),i=1,ngx) ! Cf_t
    write(69) (Cfi_i(i,3),i=1,ngx) ! Cf_g
    write(69) (Cfi_l(i,1),i=1,ngx) ! Cf_v
    write(69) (Cfi_l(i,2),i=1,ngx) ! Cf_t
    write(69) (Cfi_l(i,3),i=1,ngx) ! Cf_g
    write(69) (Cfi_o(i,1),i=1,ngx) ! Cf_v
    write(69) (Cfi_o(i,2),i=1,ngx) ! Cf_t
    write(69) (Cfi_o(i,3),i=1,ngx) ! Cf_g
    close(69)

    ! free memory
    ! -----------
    deallocate(Cfi,Cfi_i,Cfi_l,Cfi_o)
    deallocate(wrk1,wrk2,jb,fact)
    deallocate(duxfavre,duyfavre)
    
  end subroutine pp_RD_stbl

  !============================================================================
  subroutine pp_FIK2_stbl
  !============================================================================
    !> twofold-FIK decomposition of mean skin friction (Wenzel et al, JFM 2022)
    !> (see also Xu, Wang, Chen JFM 2022)
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i,j,id
    integer, dimension(:), allocatable :: jb
    real(wp) :: yb,u_delt,rho_delt
    real(wp) :: arg1,arg2
    real(wp), dimension(:), allocatable :: fact
    real(wp), dimension(:,:), allocatable :: wrk1,wrk2
    ! derivative of Favre velocity components
    real(wp), dimension(:,:), allocatable :: duxfavre,duyfavre
    ! -------------------------------------------------------------------------
    
    ! Memory allocations
    ! ------------------
    ! Cf decomposition
    allocate(Cfi(ngx,9))
    Cfi=0.0_wp
    allocate(Cfi_i(ngx,8),Cfi_l(ngx,8),Cfi_o(ngx,8))
    Cfi_i=0.0_wp
    Cfi_l=0.0_wp
    Cfi_o=0.0_wp
    ! integration bound
    allocate(jb(ngx),fact(ngx))
    ! wrk arrays
    allocate(wrk1(ngx,ngy),wrk2(ngx,ngy))
    allocate(duxfavre(ngx,ngy),duyfavre(ngx,ngy))
    
    ! twofold-FIK decomposition [Notations of Xu, Wang, Chen JFM 2022]
    ! =========================
    ! -> Cfi(:,1): Cf_B  (mean boundary-layer term)
    ! -> Cfi(:,2): Cf_V  (viscous boundary-layer term)
    ! -> Cfi(:,3): Cf_T  (turbulent-convection term)
    ! -> Cfi(:,4): Cf_M  (vertical mean-convection term)
    ! -> Cfi(:,5): Cf_D1 (streamwise mean-convection term)
    ! -> Cfi(:,6): Cf_D2 (spatial-development term #1)
    ! -> Cfi(:,7): Cf_D3 (spatial-development term #2)
    ! -> Cfi(:,8): Cf_D4 (spatial-development term #3)
    ! -> Cfi(:,9): Cf reconstruction (sum)

    ! Normalization factor
    ! --------------------
    id=0
    do i=1,ngx
       jb(i)=j99(i)+id
       yb=yg(jb(i))
       !rho_delt=rho_ref
       !u_delt=u_ref
       rho_delt=rhom(i,jb(i))
       u_delt=um(i,jb(i))
       fact(i)=2.0_wp/rho_delt/u_delt**2/yb
    enddo
 
    ! Useful quantities
    ! -----------------   
    call deriv2_x_11pts(umfavre,duxfavre)
    call deriv2_y_11pts(umfavre,duyfavre)
 
    ! Integration for boundary-layer terms
    ! ------------------------------------
    wrk1=mum*(duym+dvxm)
    wrk2=tau12m-wrk1
    
    do i=1,ngx
       do j=1,jb(i)
          arg1=wrk1(i,j)
          arg2=wrk1(i,j+1)
          Cfi(i,1)=Cfi(i,1)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
          arg1=wrk2(i,j)
          arg2=wrk2(i,j+1)
          Cfi(i,2)=Cfi(i,2)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo
       Cfi(i,1)=Cfi(i,1)*fact(i)
       Cfi(i,2)=Cfi(i,2)*fact(i)
    enddo
   
    ! inner layer
    do i=1,ngx
       do j=1,j_inn(i)
          arg1=wrk1(i,j)
          arg2=wrk1(i,j+1)
          Cfi_i(i,1)=Cfi_i(i,1)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
          arg1=wrk2(i,j)
          arg2=wrk2(i,j+1)
          Cfi_i(i,2)=Cfi_i(i,2)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo
       Cfi_i(i,1)=Cfi_i(i,1)*fact(i)
       Cfi_i(i,2)=Cfi_i(i,2)*fact(i)
    enddo
   
    ! log layer
    do i=1,ngx
       do j=j_inn(i)+1,j_log(i)
          arg1=wrk1(i,j)
          arg2=wrk1(i,j+1)
          Cfi_l(i,1)=Cfi_l(i,1)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
          arg1=wrk2(i,j)
          arg2=wrk2(i,j+1)
          Cfi_l(i,2)=Cfi_l(i,2)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo
       Cfi_l(i,1)=Cfi_l(i,1)*fact(i)
       Cfi_l(i,2)=Cfi_l(i,2)*fact(i)
    enddo
   
    ! outer layer
    do i=1,ngx
       do j=j_log(i)+1,jb(i)
          arg1=wrk1(i,j)
          arg2=wrk1(i,j+1)
          Cfi_o(i,1)=Cfi_o(i,1)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
          arg1=wrk2(i,j)
          arg2=wrk2(i,j+1)
          Cfi_o(i,2)=Cfi_o(i,2)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo
       Cfi_o(i,1)=Cfi_o(i,1)*fact(i)
       Cfi_o(i,2)=Cfi_o(i,2)*fact(i)
    enddo
   
    ! Integration for turbulent-convection term
    ! -----------------------------------------
    ! <rho> <u''v''>
    wrk1=rhouvm-rhoum*rhovm/rhom
    
    do i=1,ngx
       do j=1,jb(i)
          arg1=wrk1(i,j)
          arg2=wrk1(i,j+1)
          Cfi(i,3)=Cfi(i,3)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo       
       Cfi(i,3)=-Cfi(i,3)*fact(i)
    enddo

    ! inner layer
    do i=1,ngx
       do j=1,j_inn(i)
          arg1=wrk1(i,j)
          arg2=wrk1(i,j+1)
          Cfi_i(i,3)=Cfi_i(i,3)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo       
       Cfi_i(i,3)=-Cfi_i(i,3)*fact(i)
    enddo

    ! log layer
    do i=1,ngx
       do j=j_inn(i)+1,j_log(i)
          arg1=wrk1(i,j)
          arg2=wrk1(i,j+1)
          Cfi_l(i,3)=Cfi_l(i,3)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo       
       Cfi_l(i,3)=-Cfi_l(i,3)*fact(i)
    enddo

    ! outer layer
    do i=1,ngx
       do j=j_log(i)+1,jb(i)
          arg1=wrk1(i,j)
          arg2=wrk1(i,j+1)
          Cfi_o(i,3)=Cfi_o(i,3)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo       
       Cfi_o(i,3)=-Cfi_o(i,3)*fact(i)
    enddo

    ! Integration for vertical mean-convection term
    ! ---------------------------------------------
    ! <rho> [v] d[u]/dy
    wrk1=rhom*vmfavre*duyfavre
    
    do i=1,ngx
       do j=1,jb(i)
          arg1=(yg(jb(i))-yg(j))  *wrk1(i,j)
          arg2=(yg(jb(i))-yg(j+1))*wrk1(i,j+1)
          Cfi(i,4)=Cfi(i,4)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo       
       Cfi(i,4)=-Cfi(i,4)*fact(i)
    enddo

    ! inner layer
    do i=1,ngx
       do j=1,j_inn(i)
          arg1=(yg(jb(i))-yg(j))  *wrk1(i,j)
          arg2=(yg(jb(i))-yg(j+1))*wrk1(i,j+1)
          Cfi_i(i,4)=Cfi_i(i,4)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo       
       Cfi_i(i,4)=-Cfi_i(i,4)*fact(i)
    enddo

    ! log layer
    do i=1,ngx
       do j=j_inn(i)+1,j_log(i)
          arg1=(yg(jb(i))-yg(j))  *wrk1(i,j)
          arg2=(yg(jb(i))-yg(j+1))*wrk1(i,j+1)
          Cfi_l(i,4)=Cfi_l(i,4)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo       
       Cfi_l(i,4)=-Cfi_l(i,4)*fact(i)
    enddo

    ! outer layer
    do i=1,ngx
       do j=j_log(i)+1,jb(i)
          arg1=(yg(jb(i))-yg(j))  *wrk1(i,j)
          arg2=(yg(jb(i))-yg(j+1))*wrk1(i,j+1)
          Cfi_o(i,4)=Cfi_o(i,4)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo       
       Cfi_o(i,4)=-Cfi_o(i,4)*fact(i)
    enddo

    ! Integration for streamwise mean-convection term
    ! -----------------------------------------------
    ! <rho> [u] d[u]/dx
    wrk1=rhom*umfavre*duxfavre
    
    do i=1,ngx
       do j=1,jb(i)
          arg1=(yg(jb(i))-yg(j))  *wrk1(i,j)
          arg2=(yg(jb(i))-yg(j+1))*wrk1(i,j+1)
          Cfi(i,5)=Cfi(i,5)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo       
       Cfi(i,5)=-Cfi(i,5)*fact(i)
    enddo

    ! inner layer
    do i=1,ngx
       do j=1,j_inn(i)
          arg1=(yg(jb(i))-yg(j))  *wrk1(i,j)
          arg2=(yg(jb(i))-yg(j+1))*wrk1(i,j+1)
          Cfi_i(i,5)=Cfi_i(i,5)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo       
       Cfi_i(i,5)=-Cfi_i(i,5)*fact(i)
    enddo

    ! log layer
    do i=1,ngx
       do j=j_inn(i)+1,j_log(i)
          arg1=(yg(jb(i))-yg(j))  *wrk1(i,j)
          arg2=(yg(jb(i))-yg(j+1))*wrk1(i,j+1)
          Cfi_l(i,5)=Cfi_l(i,5)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo       
       Cfi_l(i,5)=-Cfi_l(i,5)*fact(i)
    enddo

    ! outer layer
    do i=1,ngx
       do j=j_log(i)+1,jb(i)
          arg1=(yg(jb(i))-yg(j))  *wrk1(i,j)
          arg2=(yg(jb(i))-yg(j+1))*wrk1(i,j+1)
          Cfi_o(i,5)=Cfi_o(i,5)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo       
       Cfi_o(i,5)=-Cfi_o(i,5)*fact(i)
    enddo

    ! Integration for spatial-development term #1
    ! -------------------------------------------
    ! d[p]/dx
    call deriv2_x_11pts(pm,wrk1)
    
    do i=1,ngx
       do j=1,jb(i)
          arg1=(yg(jb(i))-yg(j))  *wrk1(i,j)
          arg2=(yg(jb(i))-yg(j+1))*wrk1(i,j+1)
          Cfi(i,6)=Cfi(i,6)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo       
       Cfi(i,6)=-Cfi(i,6)*fact(i)
    enddo

    ! inner layer
    do i=1,ngx
       do j=1,j_inn(i)
          arg1=(yg(jb(i))-yg(j))  *wrk1(i,j)
          arg2=(yg(jb(i))-yg(j+1))*wrk1(i,j+1)
          Cfi_i(i,6)=Cfi_i(i,6)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo       
       Cfi_i(i,6)=-Cfi_i(i,6)*fact(i)
    enddo

    ! log layer
    do i=1,ngx
       do j=j_inn(i)+1,j_log(i)
          arg1=(yg(jb(i))-yg(j))  *wrk1(i,j)
          arg2=(yg(jb(i))-yg(j+1))*wrk1(i,j+1)
          Cfi_l(i,6)=Cfi_l(i,6)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo       
       Cfi_l(i,6)=-Cfi_l(i,6)*fact(i)
    enddo

    ! outer layer
    do i=1,ngx
       do j=j_log(i)+1,jb(i)
          arg1=(yg(jb(i))-yg(j))  *wrk1(i,j)
          arg2=(yg(jb(i))-yg(j+1))*wrk1(i,j+1)
          Cfi_o(i,6)=Cfi_o(i,6)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo       
       Cfi_o(i,6)=-Cfi_o(i,6)*fact(i)
    enddo

    ! Integration for spatial-development term #2
    ! -------------------------------------------
    ! d[tau_xx]/dx
    call deriv2_x_11pts(tau11m,wrk1)
    
    do i=1,ngx
       do j=1,jb(i)
          arg1=(yg(jb(i))-yg(j))  *wrk1(i,j)
          arg2=(yg(jb(i))-yg(j+1))*wrk1(i,j+1)
          Cfi(i,7)=Cfi(i,7)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo       
       Cfi(i,7)=Cfi(i,7)*fact(i)
    enddo

    ! inner layer
    do i=1,ngx
       do j=1,j_inn(i)
          arg1=(yg(jb(i))-yg(j))  *wrk1(i,j)
          arg2=(yg(jb(i))-yg(j+1))*wrk1(i,j+1)
          Cfi_i(i,7)=Cfi_i(i,7)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo       
       Cfi_i(i,7)=Cfi_i(i,7)*fact(i)
    enddo

    ! log layer
    do i=1,ngx
       do j=j_inn(i)+1,j_log(i)
          arg1=(yg(jb(i))-yg(j))  *wrk1(i,j)
          arg2=(yg(jb(i))-yg(j+1))*wrk1(i,j+1)
          Cfi_l(i,7)=Cfi_l(i,7)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo       
       Cfi_l(i,7)=Cfi_l(i,7)*fact(i)
    enddo

    ! outer layer
    do i=1,ngx
       do j=j_log(i)+1,jb(i)
          arg1=(yg(jb(i))-yg(j))  *wrk1(i,j)
          arg2=(yg(jb(i))-yg(j+1))*wrk1(i,j+1)
          Cfi_o(i,7)=Cfi_o(i,7)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo       
       Cfi_o(i,7)=Cfi_o(i,7)*fact(i)
    enddo

    ! Integration for spatial-development term #3
    ! -------------------------------------------
    ! <rho><u''u''>
    wrk1=rhouum-rhoum*rhoum/rhom
    ! d/dx(<rho><u''u''>)
    call deriv2_x_11pts(wrk1,wrk2)
    
    do i=1,ngx
       do j=1,jb(i)
          arg1=(yg(jb(i))-yg(j))  *wrk2(i,j)
          arg2=(yg(jb(i))-yg(j+1))*wrk2(i,j+1)
          Cfi(i,8)=Cfi(i,8)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo       
       Cfi(i,8)=-Cfi(i,8)*fact(i)
    enddo
 
    ! inner layer
    do i=1,ngx
       do j=1,j_inn(i)
          arg1=(yg(jb(i))-yg(j))  *wrk2(i,j)
          arg2=(yg(jb(i))-yg(j+1))*wrk2(i,j+1)
          Cfi_i(i,8)=Cfi_i(i,8)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo       
       Cfi_i(i,8)=-Cfi_i(i,8)*fact(i)
    enddo
 
    ! log layer
    do i=1,ngx
       do j=j_inn(i)+1,j_log(i)
          arg1=(yg(jb(i))-yg(j))  *wrk2(i,j)
          arg2=(yg(jb(i))-yg(j+1))*wrk2(i,j+1)
          Cfi_l(i,8)=Cfi_l(i,8)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo       
       Cfi_l(i,8)=-Cfi_l(i,8)*fact(i)
    enddo
 
    ! outer layer
    do i=1,ngx
       do j=j_log(i)+1,jb(i)
          arg1=(yg(jb(i))-yg(j))  *wrk2(i,j)
          arg2=(yg(jb(i))-yg(j+1))*wrk2(i,j+1)
          Cfi_o(i,8)=Cfi_o(i,8)+0.5_wp*(arg1+arg2)*(yg(j+1)-yg(j))
       enddo       
       Cfi_o(i,8)=-Cfi_o(i,8)*fact(i)
    enddo
 
    ! Summation for check of Cf reconstruction
    ! ----------------------------------------
    do j=1,8
       Cfi(:,9)=Cfi(:,9)+Cfi(:,j)
    enddo
    
    ! Write binary file
    ! -----------------
    open(69,file=trim(dirRESU)//'FIK2_Cf'//trim(name_output)//'.bin',form='unformatted',status='unknown')
    rewind(69)
    write(69) ngx
    do j=1,9
       write(69) (Cfi(i,j),i=1,ngx)
    enddo
    write(69) (c_f(i),i=1,ngx) ! Cf computed
    write(69) (Re_theta(i),i=1,ngx) ! Cf computed
    write(69) (Re_ds(i),i=1,ngx) ! Cf computed
    write(69) (Re_d99(i),i=1,ngx) ! Cf computed
    close(69)

    open(69,file=trim(dirRESU)//'FIK2_Cf'//trim(name_output)//'_layers.bin',form='unformatted',status='unknown')
    rewind(69)
    write(69) ngx
    do j=1,8
       write(69) (Cfi_i(i,j),i=1,ngx)
    enddo
    do j=1,8
       write(69) (Cfi_l(i,j),i=1,ngx)
    enddo
    do j=1,8
       write(69) (Cfi_o(i,j),i=1,ngx)
    enddo
    close(69)

    ! free memory
    ! -----------
    deallocate(Cfi,Cfi_i,Cfi_l,Cfi_o)
    deallocate(wrk1,wrk2,jb,fact)
    deallocate(duxfavre,duyfavre)
    
  end subroutine pp_FIK2_stbl

end module mod_pp_stats_skin_friction
