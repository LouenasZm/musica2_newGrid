!==============================================================================
module mod_pp_stats_read
!==============================================================================
  !> Module to read and assign stats I/O
!==============================================================================
  use mod_grid
  use mod_io
  use mod_io_stats
  use mod_pp_stats_var
  implicit none
  ! ---------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------

contains

  !============================================================================
  subroutine pp_stats_read_chan
  !============================================================================
    !> Read stats.dat and assign variables for channel flow
  !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    logical :: iexist
    ! -------------------------------------------------------------------------
    
!!    ubulk=0.0_wp !!!!!!!!????????????????????????????????????????????????????????????????????????????

!!$    ! Allocation
!!$    ! ----------
!!$    allocate( dumdy(ngy), dvmdy(ngy), dTmdy(ngy), dpmdy(ngy) &
!!$         ,dmumdy(ngy),dromdy(ngy), shear(ngy),   ygw(ngy) )
!!$    allocate(xgp(ngx), ygp(ngy), zgp(ngz))
!!$    allocate(xgs(ngy), ygs(ngy), zgs(ngy))

    ! Check if restart file exists, if not disable analysis of 3D fields
    ! ------------------------------------------------------------------
    !if (idepart.eq.POST_PROCESSING) then
    !  call read_start_mpi('     ')
    !else
!!    i3df=.false.
    !endif
 
    ! Read stats I/O file
    ! ===================
    call read_write_stats_chan(READ)

    ! Check if stats file exists, if yes allocate and fill tables
    ! ===========================================================
    inquire(file='stats.dat',exist=iexist)
    
    if (iexist) then
       
       ! Allocate variables for analysis
       ! ===============================
       allocate(rhom     (1,ngy))
       allocate(um       (1,ngy))
       allocate(vm       (1,ngy))
       allocate(wm       (1,ngy))
       allocate(pm       (1,ngy))
       allocate(Tm       (1,ngy))
       allocate(em       (1,ngy))
       allocate(hm       (1,ngy))
       allocate(cm       (1,ngy))
       allocate(sm       (1,ngy))
       allocate(Mm       (1,ngy))
       allocate(ktm      (1,ngy))
       allocate(Gm       (1,ngy))
       allocate(mum      (1,ngy))
       allocate(lam      (1,ngy))
       allocate(cpm      (1,ngy))
       allocate(cvm      (1,ngy))
       allocate(prm      (1,ngy))
       allocate(eckm     (1,ngy))
       allocate(divm     (1,ngy))
       allocate(rhoduxm  (1,ngy))
       allocate(rhoduym  (1,ngy))
       allocate(rhoduzm  (1,ngy))
       allocate(rhodvxm  (1,ngy))
       allocate(rhodvym  (1,ngy))
       allocate(rhodvzm  (1,ngy))
       allocate(rhodwxm  (1,ngy))
       allocate(rhodwym  (1,ngy))
       allocate(rhodwzm  (1,ngy))
       allocate(pdivm    (1,ngy))
       allocate(rhodivm  (1,ngy))
       allocate(vrtxm    (1,ngy))
       allocate(vrtym    (1,ngy))
       allocate(vrtzm    (1,ngy))
       allocate(Ttotm    (1,ngy))
       allocate(ducrm    (1,ngy))
       allocate(rhoum    (1,ngy))
       allocate(rhovm    (1,ngy))
       allocate(rhowm    (1,ngy))
       allocate(rhoem    (1,ngy))
       allocate(rhoTm    (1,ngy))
       allocate(rho2m    (1,ngy))
       allocate(u2m      (1,ngy))
       allocate(v2m      (1,ngy))
       allocate(w2m      (1,ngy))
       allocate(uvm      (1,ngy))
       allocate(uwm      (1,ngy))
       allocate(vwm      (1,ngy))
       allocate(p2m      (1,ngy))
       allocate(T2m      (1,ngy))
       allocate(e2m      (1,ngy))
       allocate(h2m      (1,ngy))
       allocate(ccm      (1,ngy))
       allocate(s2m      (1,ngy))
       allocate(M2m      (1,ngy))
       allocate(G2m      (1,ngy))
       allocate(mu2m     (1,ngy))
       allocate(la2m     (1,ngy))
       allocate(cv2m     (1,ngy))
       allocate(cp2m     (1,ngy))
       allocate(pr2m     (1,ngy))
       allocate(eck2m    (1,ngy))
       allocate(upm      (1,ngy))
       allocate(vpm      (1,ngy))
       allocate(uTm      (1,ngy))
       allocate(vTm      (1,ngy))
       allocate(usm      (1,ngy))
       allocate(vsm      (1,ngy))
       allocate(rhopm    (1,ngy))
       allocate(rhoTm    (1,ngy))
       allocate(rhohm    (1,ngy))
       allocate(pTm      (1,ngy))
       allocate(spm      (1,ngy))
       allocate(sTm      (1,ngy))
       allocate(rhosm    (1,ngy))
       allocate(Grhom    (1,ngy))
       allocate(Gpm      (1,ngy))
       allocate(Gsm      (1,ngy))
       allocate(GTm      (1,ngy))
       allocate(Gum      (1,ngy))
       allocate(Gvm      (1,ngy))
       allocate(pduxm    (1,ngy))
       allocate(pdvym    (1,ngy))
       allocate(pdwzm    (1,ngy))
       allocate(pduym    (1,ngy))
       allocate(pdvxm    (1,ngy))
       allocate(div2m    (1,ngy))
       allocate(rhodiv2m (1,ngy))
       allocate(dux2m    (1,ngy))
       allocate(duy2m    (1,ngy))
       allocate(duz2m    (1,ngy))
       allocate(dvx2m    (1,ngy))
       allocate(dvy2m    (1,ngy))
       allocate(dvz2m    (1,ngy))
       allocate(dwx2m    (1,ngy))
       allocate(dwy2m    (1,ngy))
       allocate(dwz2m    (1,ngy))
       allocate(vrtx2m   (1,ngy))
       allocate(vrty2m   (1,ngy))
       allocate(vrtz2m   (1,ngy))
       allocate(rhovrtxm (1,ngy))
       allocate(rhovrtym (1,ngy))
       allocate(rhovrtzm (1,ngy))
       allocate(Ttot2m   (1,ngy))
       allocate(ducr2m   (1,ngy))
       allocate(rhouum   (1,ngy))
       allocate(rhovvm   (1,ngy))
       allocate(rhowwm   (1,ngy))
       allocate(rhoTTm   (1,ngy))
       allocate(rhovrtx2m(1,ngy))
       allocate(rhovrty2m(1,ngy))
       allocate(rhovrtz2m(1,ngy))
       allocate(rhouvm   (1,ngy))
       allocate(rhovwm   (1,ngy))
       allocate(rhovTm   (1,ngy))
       allocate(rhouuvm  (1,ngy))
       allocate(rhovvvm  (1,ngy))
       allocate(rhowwvm  (1,ngy))
       allocate(rhouvvm  (1,ngy))
       allocate(rhodux2m (1,ngy))
       allocate(rhodvy2m (1,ngy))
       allocate(rhodwz2m (1,ngy))
       allocate(rhoduydvxm(1,ngy))
       allocate(rhoduzdwxm(1,ngy))
       allocate(rhodvzdwym(1,ngy))
       allocate(u3m      (1,ngy))
       allocate(p3m      (1,ngy))
       allocate(u4m      (1,ngy))
       allocate(p4m      (1,ngy))
       allocate(tau11m   (1,ngy))
       allocate(tau12m   (1,ngy))
       allocate(tau13m   (1,ngy))
       allocate(tau22m   (1,ngy))
       allocate(tau23m   (1,ngy))
       allocate(tau33m   (1,ngy))
       allocate(tau12um  (1,ngy))
       allocate(tau12vm  (1,ngy))
       allocate(tau22um  (1,ngy))
       allocate(tau22vm  (1,ngy))
       allocate(tau23wm  (1,ngy))
       allocate(tau11duxm(1,ngy))
       allocate(tau11dvxm(1,ngy))
       allocate(tau12duxm(1,ngy))
       allocate(tau12duym(1,ngy))
       allocate(tau12dvxm(1,ngy))
       allocate(tau12dvym(1,ngy))
       allocate(tau13duzm(1,ngy))
       allocate(tau13dvzm(1,ngy))
       allocate(tau13dwxm(1,ngy))
       allocate(tau22duym(1,ngy))
       allocate(tau22dvym(1,ngy))
       allocate(tau23duzm(1,ngy))
       allocate(tau23dvzm(1,ngy))
       allocate(tau23dwym(1,ngy))
       allocate(tau33dwzm(1,ngy))
       allocate(ladTxm   (1,ngy))
       allocate(ladTym   (1,ngy))
       allocate(ladTzm   (1,ngy))
       allocate(hum      (1,ngy))
       allocate(hvm      (1,ngy))
       allocate(hwm      (1,ngy))
       allocate(rhohum   (1,ngy))
       allocate(rhohvm   (1,ngy))
       allocate(rhohwm   (1,ngy))
       allocate(rhouuum  (1,ngy))
       !!allocate(rhovvvm  (1,ngy)) already defined
       allocate(rhowwwm  (1,ngy))
       allocate(rhofm    (1,ngy))
       ! allocate(Drhom   (1,ngy))
       ! allocate(Drhoum  (1,ngy))
       ! allocate(Drhovm  (1,ngy))
       ! allocate(Drhowm  (1,ngy))
       ! allocate(Drhoem  (1,ngy))
       ! allocate(uDrom   (1,ngy))
       ! allocate(vDrom   (1,ngy))
       ! allocate(wDrom   (1,ngy))
       ! allocate(uDroum  (1,ngy))
       ! allocate(vDrovm  (1,ngy))
       ! allocate(wDrowm  (1,ngy))

       allocate(rhoukuk  (1,ngy))
       allocate(kolmog   (1,ngy))
       allocate(dissip   (1,ngy))

       ! Assign variables for analysis
       ! =============================
       rhom      (1,:) = avg_tg(1,:,  1)
       um        (1,:) = avg_tg(1,:,  2)
       vm        (1,:) = avg_tg(1,:,  3)
       wm        (1,:) = avg_tg(1,:,  4)
       pm        (1,:) = avg_tg(1,:,  5)
       Tm        (1,:) = avg_tg(1,:,  6)
       em        (1,:) = avg_tg(1,:,  7)
       hm        (1,:) = avg_tg(1,:,  8)
       cm        (1,:) = avg_tg(1,:,  9)
       sm        (1,:) = avg_tg(1,:, 10)
       Mm        (1,:) = avg_tg(1,:, 11)
       ktm       (1,:) = avg_tg(1,:, 12)
       Gm        (1,:) = avg_tg(1,:, 13)
       mum       (1,:) = avg_tg(1,:, 14)
       lam       (1,:) = avg_tg(1,:, 15)
       cpm       (1,:) = avg_tg(1,:, 16)
       cvm       (1,:) = avg_tg(1,:, 17)
       prm       (1,:) = avg_tg(1,:, 18)
       eckm      (1,:) = avg_tg(1,:, 19)
       divm      (1,:) = avg_tg(1,:, 20)
       rhoduxm   (1,:) = avg_tg(1,:, 21)
       rhoduym   (1,:) = avg_tg(1,:, 22)
       rhoduzm   (1,:) = avg_tg(1,:, 23)
       rhodvxm   (1,:) = avg_tg(1,:, 24)
       rhodvym   (1,:) = avg_tg(1,:, 25)
       rhodvzm   (1,:) = avg_tg(1,:, 26)
       rhodwxm   (1,:) = avg_tg(1,:, 27)
       rhodwym   (1,:) = avg_tg(1,:, 28)
       rhodwzm   (1,:) = avg_tg(1,:, 29)
       pdivm     (1,:) = avg_tg(1,:, 30)
       rhodivm   (1,:) = avg_tg(1,:, 31)
       vrtxm     (1,:) = avg_tg(1,:, 32)
       vrtym     (1,:) = avg_tg(1,:, 33)
       vrtzm     (1,:) = avg_tg(1,:, 34)
       Ttotm     (1,:) = avg_tg(1,:, 35)
       ducrm     (1,:) = avg_tg(1,:, 36)
       rhoum     (1,:) = avg_tg(1,:, 37)
       rhovm     (1,:) = avg_tg(1,:, 38)
       rhowm     (1,:) = avg_tg(1,:, 39)
       rhoem     (1,:) = avg_tg(1,:, 40)
       rhoTm     (1,:) = avg_tg(1,:, 41)
       rho2m     (1,:) = avg_tg(1,:, 42)
       u2m       (1,:) = avg_tg(1,:, 43)
       v2m       (1,:) = avg_tg(1,:, 44)
       w2m       (1,:) = avg_tg(1,:, 45)
       uvm       (1,:) = avg_tg(1,:, 46)
       uwm       (1,:) = avg_tg(1,:, 47)
       vwm       (1,:) = avg_tg(1,:, 48)
       vTm       (1,:) = avg_tg(1,:, 49)
       p2m       (1,:) = avg_tg(1,:, 50)
       T2m       (1,:) = avg_tg(1,:, 51)
       e2m       (1,:) = avg_tg(1,:, 52)
       h2m       (1,:) = avg_tg(1,:, 53)
       ccm       (1,:) = avg_tg(1,:, 54)
       s2m       (1,:) = avg_tg(1,:, 55)
       M2m       (1,:) = avg_tg(1,:, 56)
       G2m       (1,:) = avg_tg(1,:, 57)
       mu2m      (1,:) = avg_tg(1,:, 58)
       la2m      (1,:) = avg_tg(1,:, 59)
       cv2m      (1,:) = avg_tg(1,:, 60)
       cp2m      (1,:) = avg_tg(1,:, 61)
       pr2m      (1,:) = avg_tg(1,:, 62)
       eck2m     (1,:) = avg_tg(1,:, 63)
       upm       (1,:) = avg_tg(1,:, 64)
       vpm       (1,:) = avg_tg(1,:, 65)
       uTm       (1,:) = avg_tg(1,:, 66)
       vTm       (1,:) = avg_tg(1,:, 67)
       usm       (1,:) = avg_tg(1,:, 68)
       vsm       (1,:) = avg_tg(1,:, 69)
       rhopm     (1,:) = avg_tg(1,:, 70)
       rhoTm     (1,:) = avg_tg(1,:, 71)
       rhohm     (1,:) = avg_tg(1,:, 72)
       pTm       (1,:) = avg_tg(1,:, 73)
       spm       (1,:) = avg_tg(1,:, 74)
       sTm       (1,:) = avg_tg(1,:, 75)
       rhosm     (1,:) = avg_tg(1,:, 76)
       Grhom     (1,:) = avg_tg(1,:, 77)
       Gpm       (1,:) = avg_tg(1,:, 78)
       Gsm       (1,:) = avg_tg(1,:, 79)
       GTm       (1,:) = avg_tg(1,:, 80)
       Gum       (1,:) = avg_tg(1,:, 81)
       Gvm       (1,:) = avg_tg(1,:, 82)
       pduxm     (1,:) = avg_tg(1,:, 83)
       pdvym     (1,:) = avg_tg(1,:, 84)
       pdwzm     (1,:) = avg_tg(1,:, 85)
       pduym     (1,:) = avg_tg(1,:, 86)
       pdvxm     (1,:) = avg_tg(1,:, 87)
       div2m     (1,:) = avg_tg(1,:, 88)
       rhodiv2m  (1,:) = avg_tg(1,:, 89)
       dux2m     (1,:) = avg_tg(1,:, 90)
       duy2m     (1,:) = avg_tg(1,:, 91)
       duz2m     (1,:) = avg_tg(1,:, 92)
       dvx2m     (1,:) = avg_tg(1,:, 93)
       dvy2m     (1,:) = avg_tg(1,:, 94)
       dvz2m     (1,:) = avg_tg(1,:, 95)
       dwx2m     (1,:) = avg_tg(1,:, 96)
       dwy2m     (1,:) = avg_tg(1,:, 97)
       dwz2m     (1,:) = avg_tg(1,:, 98)
       vrtx2m    (1,:) = avg_tg(1,:, 99)
       vrty2m    (1,:) = avg_tg(1,:,100)
       vrtz2m    (1,:) = avg_tg(1,:,101)
       rhovrtxm  (1,:) = avg_tg(1,:,102)
       rhovrtym  (1,:) = avg_tg(1,:,103)
       rhovrtzm  (1,:) = avg_tg(1,:,104)
       Ttot2m    (1,:) = avg_tg(1,:,105)
       ducr2m    (1,:) = avg_tg(1,:,106)
       rhouum    (1,:) = avg_tg(1,:,107)
       rhovvm    (1,:) = avg_tg(1,:,108)
       rhowwm    (1,:) = avg_tg(1,:,109)
       rhoTTm    (1,:) = avg_tg(1,:,110)
       rhovrtx2m (1,:) = avg_tg(1,:,111)
       rhovrty2m (1,:) = avg_tg(1,:,112)
       rhovrtz2m (1,:) = avg_tg(1,:,113)
       rhouvm    (1,:) = avg_tg(1,:,114)
       rhovwm    (1,:) = avg_tg(1,:,115)
       rhovTm    (1,:) = avg_tg(1,:,116)
       rhouuvm   (1,:) = avg_tg(1,:,117)
       rhovvvm   (1,:) = avg_tg(1,:,118)
       rhowwvm   (1,:) = avg_tg(1,:,119)
       rhouvvm   (1,:) = avg_tg(1,:,120)
       rhodux2m  (1,:) = avg_tg(1,:,121)
       rhodvy2m  (1,:) = avg_tg(1,:,122)
       rhodwz2m  (1,:) = avg_tg(1,:,123)
       rhoduydvxm(1,:) = avg_tg(1,:,124)
       rhoduzdwxm(1,:) = avg_tg(1,:,125)
       rhodvzdwym(1,:) = avg_tg(1,:,126)
       u3m       (1,:) = avg_tg(1,:,127)
       p3m       (1,:) = avg_tg(1,:,128)
       u4m       (1,:) = avg_tg(1,:,129)
       p4m       (1,:) = avg_tg(1,:,130)
       tau11m    (1,:) = -avg_tg(1,:,131)
       tau12m    (1,:) = -avg_tg(1,:,132)
       tau13m    (1,:) = -avg_tg(1,:,133)
       tau22m    (1,:) = -avg_tg(1,:,134)
       tau23m    (1,:) = -avg_tg(1,:,135)
       tau33m    (1,:) = -avg_tg(1,:,136)
       tau12um   (1,:) = -avg_tg(1,:,137)
       tau12vm   (1,:) = -avg_tg(1,:,138)
       tau22um   (1,:) = -avg_tg(1,:,139)
       tau22vm   (1,:) = -avg_tg(1,:,140)
       tau23wm   (1,:) = -avg_tg(1,:,141)
       tau11duxm (1,:) = -avg_tg(1,:,142)
       tau11dvxm (1,:) = -avg_tg(1,:,143)
       tau12duxm (1,:) = -avg_tg(1,:,144)
       tau12duym (1,:) = -avg_tg(1,:,145)
       tau12dvxm (1,:) = -avg_tg(1,:,146)
       tau12dvym (1,:) = -avg_tg(1,:,147)
       tau13duzm (1,:) = -avg_tg(1,:,148)
       tau13dvzm (1,:) = -avg_tg(1,:,149)
       tau13dwxm (1,:) = -avg_tg(1,:,150)
       tau22duym (1,:) = -avg_tg(1,:,151)
       tau22dvym (1,:) = -avg_tg(1,:,152)
       tau23duzm (1,:) = -avg_tg(1,:,153)
       tau23dvzm (1,:) = -avg_tg(1,:,154)
       tau23dwym (1,:) = -avg_tg(1,:,155)
       tau33dwzm (1,:) = -avg_tg(1,:,156)
       ladTxm    (1,:) = avg_tg(1,:,157)
       ladTym    (1,:) = avg_tg(1,:,158)
       ladTzm    (1,:) = avg_tg(1,:,159)
       hum       (1,:) = avg_tg(1,:,160)
       hvm       (1,:) = avg_tg(1,:,161)
       hwm       (1,:) = avg_tg(1,:,162)
       rhohum    (1,:) = avg_tg(1,:,163)
       rhohvm    (1,:) = avg_tg(1,:,164)
       rhohwm    (1,:) = avg_tg(1,:,165)
       rhouuum   (1,:) = avg_tg(1,:,166)
       rhovvvm   (1,:) = avg_tg(1,:,167)
       rhowwwm   (1,:) = avg_tg(1,:,168)
       rhofm     (1,:) = avg_tg(1,:,169)
       !Drhom     (1,:) = avg_tg(1,:,170)
       !Drhoum    (1,:) = avg_tg(1,:,171)
       !Drhovm    (1,:) = avg_tg(1,:,172)
       !Drhowm    (1,:) = avg_tg(1,:,173)
       !Drhoem    (1,:) = avg_tg(1,:,174)
       !uDrom     (1,:) = avg_tg(1,:,175)
       !vDrom     (1,:) = avg_tg(1,:,176)
       !wDrom     (1,:) = avg_tg(1,:,177)
       !uDroum    (1,:) = avg_tg(1,:,178)
       !vDrovm    (1,:) = avg_tg(1,:,179)
       !wDrowm    (1,:) = avg_tg(1,:,180)
    endif

    deallocate(avg_t)
     
!!$    ! Normalization w.r.t half-channel height
!!$    ! =======================================
!!$    hc=yg(ngy)
!!$    xg=xg/hc
!!$    yg=yg/hc
!!$    zg=zg/hc
!!$    Lx=xg(ngx)-xg(1)
!!$    Ly=yg(ngy)-yg(1)
!!$    Lz=zg(ngz)-zg(1)
!!$
!!$    if (iproc==0) then
!!$       print *,'dim',hc,Lx,Ly,Lz
!!$       print *,'i3df',i3df,'itec',itec
!!$    endif
  
  end subroutine pp_stats_read_chan

  !============================================================================
  subroutine pp_stats_read_xy
    !============================================================================
    !> Read stats.bin and assign variables for xy field
    !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    ! -------------------------------------------------------------------------

    ! Read stats I/O file
    ! ===================
    call read_write_stats_xy(READ,iblc_pp)

    ! Allocate variables for analysis
    ! ===============================
    allocate(rhom     (ngx,ngy))
    allocate(um       (ngx,ngy))
    allocate(vm       (ngx,ngy))
    allocate(wm       (ngx,ngy))
    allocate(pm       (ngx,ngy))
    allocate(Tm       (ngx,ngy))
    allocate(em       (ngx,ngy))
    allocate(hm       (ngx,ngy))
    allocate(cm       (ngx,ngy))
    allocate(sm       (ngx,ngy))
    allocate(Mm       (ngx,ngy))
    allocate(ktm      (ngx,ngy))
    allocate(Gm       (ngx,ngy))
    allocate(mum      (ngx,ngy))
    allocate(lam      (ngx,ngy))
    allocate(cpm      (ngx,ngy))
    allocate(cvm      (ngx,ngy))
    allocate(prm      (ngx,ngy))
    allocate(eckm     (ngx,ngy))
    allocate(divm     (ngx,ngy))
    allocate(rhoduxm  (ngx,ngy))
    allocate(rhoduym  (ngx,ngy))
    allocate(rhoduzm  (ngx,ngy))
    allocate(rhodvxm  (ngx,ngy))
    allocate(rhodvym  (ngx,ngy))
    allocate(rhodvzm  (ngx,ngy))
    allocate(rhodwxm  (ngx,ngy))
    allocate(rhodwym  (ngx,ngy))
    allocate(rhodwzm  (ngx,ngy))
    allocate(pdivm    (ngx,ngy))
    allocate(rhodivm  (ngx,ngy))
    allocate(vrtxm    (ngx,ngy))
    allocate(vrtym    (ngx,ngy))
    allocate(vrtzm    (ngx,ngy))
    allocate(rhoum    (ngx,ngy))
    allocate(rhovm    (ngx,ngy))
    allocate(rhowm    (ngx,ngy))
    allocate(rhoem    (ngx,ngy))
    allocate(rhoTm    (ngx,ngy))
    allocate(rho2m    (ngx,ngy))
    allocate(u2m      (ngx,ngy))
    allocate(v2m      (ngx,ngy))
    allocate(w2m      (ngx,ngy))
    allocate(uvm      (ngx,ngy))
    allocate(uwm      (ngx,ngy))
    allocate(vwm      (ngx,ngy))
    allocate(p2m      (ngx,ngy))
    allocate(T2m      (ngx,ngy))
    allocate(e2m      (ngx,ngy))
    allocate(h2m      (ngx,ngy))
    allocate(ccm      (ngx,ngy))
    allocate(s2m      (ngx,ngy))
    allocate(M2m      (ngx,ngy))
    allocate(G2m      (ngx,ngy))
    allocate(mu2m     (ngx,ngy))
    allocate(la2m     (ngx,ngy))
    allocate(cv2m     (ngx,ngy))
    allocate(cp2m     (ngx,ngy))
    allocate(pr2m     (ngx,ngy))
    allocate(eck2m    (ngx,ngy))
    allocate(upm      (ngx,ngy))
    allocate(vpm      (ngx,ngy))
    allocate(uTm      (ngx,ngy))
    allocate(vTm      (ngx,ngy))
    allocate(usm      (ngx,ngy))
    allocate(vsm      (ngx,ngy))
    allocate(rhopm    (ngx,ngy))
    allocate(rhohm    (ngx,ngy))
    allocate(pTm      (ngx,ngy))
    allocate(spm      (ngx,ngy))
    allocate(sTm      (ngx,ngy))
    allocate(rhosm    (ngx,ngy))
    allocate(Grhom    (ngx,ngy))
    allocate(Gpm      (ngx,ngy))
    allocate(Gsm      (ngx,ngy))
    allocate(GTm      (ngx,ngy))
    allocate(Gum      (ngx,ngy))
    allocate(Gvm      (ngx,ngy))
    allocate(pduxm    (ngx,ngy))
    allocate(pdvym    (ngx,ngy))
    allocate(pdwzm    (ngx,ngy))
    allocate(pduym    (ngx,ngy))
    allocate(pdvxm    (ngx,ngy))
    allocate(div2m    (ngx,ngy))
    allocate(rhodiv2m (ngx,ngy))
    allocate(dux2m    (ngx,ngy))
    allocate(duy2m    (ngx,ngy))
    allocate(duz2m    (ngx,ngy))
    allocate(dvx2m    (ngx,ngy))
    allocate(dvy2m    (ngx,ngy))
    allocate(dvz2m    (ngx,ngy))
    allocate(dwx2m    (ngx,ngy))
    allocate(dwy2m    (ngx,ngy))
    allocate(dwz2m    (ngx,ngy))
    allocate(vrtx2m   (ngx,ngy))
    allocate(vrty2m   (ngx,ngy))
    allocate(vrtz2m   (ngx,ngy))
    allocate(rhovrtxm (ngx,ngy))
    allocate(rhovrtym (ngx,ngy))
    allocate(rhovrtzm (ngx,ngy))
    allocate(rhouum   (ngx,ngy))
    allocate(rhovvm   (ngx,ngy))
    allocate(rhowwm   (ngx,ngy))
    allocate(rhoTTm   (ngx,ngy))
    allocate(rhovrtx2m(ngx,ngy))
    allocate(rhovrty2m(ngx,ngy))
    allocate(rhovrtz2m(ngx,ngy))
    allocate(rhouvm   (ngx,ngy))
    allocate(rhouwm   (ngx,ngy))
    allocate(rhovwm   (ngx,ngy))
    allocate(rhovTm   (ngx,ngy))
    allocate(rhouuvm  (ngx,ngy))
    allocate(rhovvvm  (ngx,ngy))
    allocate(rhowwum  (ngx,ngy))
    allocate(rhowwvm  (ngx,ngy))
    allocate(rhouvvm  (ngx,ngy))
    allocate(rhodux2m (ngx,ngy))
    allocate(rhodvy2m (ngx,ngy))
    allocate(rhodwz2m (ngx,ngy))
    allocate(rhoduydvxm(ngx,ngy))
    allocate(rhoduzdwxm(ngx,ngy))
    allocate(rhodvzdwym(ngx,ngy))
    allocate(u3m      (ngx,ngy))
    allocate(p3m      (ngx,ngy))
    allocate(u4m      (ngx,ngy))
    allocate(p4m      (ngx,ngy))
    allocate(tau11m   (ngx,ngy))
    allocate(tau12m   (ngx,ngy))
    allocate(tau13m   (ngx,ngy))
    allocate(tau22m   (ngx,ngy))
    allocate(tau23m   (ngx,ngy))
    allocate(tau33m   (ngx,ngy))
    allocate(tau11um  (ngx,ngy))
    allocate(tau12um  (ngx,ngy))
    allocate(tau12vm  (ngx,ngy))
    allocate(tau22um  (ngx,ngy))
    allocate(tau22vm  (ngx,ngy))
    allocate(tau13wm  (ngx,ngy))
    allocate(tau23wm  (ngx,ngy))
    allocate(tau11duxm(ngx,ngy))
    allocate(tau11dvxm(ngx,ngy))
    allocate(tau12duxm(ngx,ngy))
    allocate(tau12duym(ngx,ngy))
    allocate(tau12dvxm(ngx,ngy))
    allocate(tau12dvym(ngx,ngy))
    allocate(tau13duzm(ngx,ngy))
    allocate(tau13dvzm(ngx,ngy))
    allocate(tau13dwxm(ngx,ngy))
    allocate(tau22duym(ngx,ngy))
    allocate(tau22dvym(ngx,ngy))
    allocate(tau23duzm(ngx,ngy))
    allocate(tau23dvzm(ngx,ngy))
    allocate(tau23dwym(ngx,ngy))
    allocate(tau33dwzm(ngx,ngy))
    allocate(ladTxm   (ngx,ngy))
    allocate(ladTym   (ngx,ngy))
    allocate(ladTzm   (ngx,ngy))
    allocate(hum      (ngx,ngy))
    allocate(hvm      (ngx,ngy))
    allocate(hwm      (ngx,ngy))
    allocate(rhohum   (ngx,ngy))
    allocate(rhohvm   (ngx,ngy))
    allocate(rhohwm   (ngx,ngy))
    allocate(rhouuum  (ngx,ngy))
    !allocate(rhovvvm  (ngx,ngy)) already defined
    allocate(rhowwwm  (ngx,ngy))
    ! allocate(Drhom   (ngx,ngy))
    ! allocate(Drhoum  (ngx,ngy))
    ! allocate(Drhovm  (ngx,ngy))
    ! allocate(Drhowm  (ngx,ngy))
    ! allocate(Drhoem  (ngx,ngy))
    ! allocate(uDrom   (ngx,ngy))
    ! allocate(vDrom   (ngx,ngy))
    ! allocate(wDrom   (ngx,ngy))
    ! allocate(uDroum  (ngx,ngy))
    ! allocate(vDrovm  (ngx,ngy))
    ! allocate(wDrowm  (ngx,ngy))

    allocate(rhoukuk  (ngx,ngy))
    allocate(kolmog   (ngx,ngy))
    allocate(dissip   (ngx,ngy))

    ! Assign variables for analysis
    ! =============================
    rhom      (:,:) = avg_t(:,:,  1)     
    um        (:,:) = avg_t(:,:,  2)
    vm        (:,:) = avg_t(:,:,  3)
    wm        (:,:) = avg_t(:,:,  4)
    pm        (:,:) = avg_t(:,:,  5)
    Tm        (:,:) = avg_t(:,:,  6)
    rhoum     (:,:) = avg_t(:,:,  7)
    rhovm     (:,:) = avg_t(:,:,  8)
    rhowm     (:,:) = avg_t(:,:,  9)
    rhoem     (:,:) = avg_t(:,:, 10)
    rho2m     (:,:) = avg_t(:,:, 11)
    u2m       (:,:) = avg_t(:,:, 12)
    v2m       (:,:) = avg_t(:,:, 13)
    w2m       (:,:) = avg_t(:,:, 14)
    uvm       (:,:) = avg_t(:,:, 15)
    uwm       (:,:) = avg_t(:,:, 16)
    vwm       (:,:) = avg_t(:,:, 17)
    vTm       (:,:) = avg_t(:,:, 18)
    p2m       (:,:) = avg_t(:,:, 19)
    T2m       (:,:) = avg_t(:,:, 20)
    mum       (:,:) = avg_t(:,:, 21)
    divm      (:,:) = avg_t(:,:, 22)
    div2m     (:,:) = avg_t(:,:, 23)

    em        (:,:) = avg_t(:,:, 24)
    hm        (:,:) = avg_t(:,:, 25)
    cm        (:,:) = avg_t(:,:, 26)
    sm        (:,:) = avg_t(:,:, 27)
    Mm        (:,:) = avg_t(:,:, 28)
    ktm       (:,:) = avg_t(:,:, 29)
    Gm        (:,:) = avg_t(:,:, 30)
    lam       (:,:) = avg_t(:,:, 31)
    cpm       (:,:) = avg_t(:,:, 32)
    cvm       (:,:) = avg_t(:,:, 33)
    prm       (:,:) = avg_t(:,:, 34)
    eckm      (:,:) = avg_t(:,:, 35)
    rhoduxm   (:,:) = avg_t(:,:, 36)
    rhoduym   (:,:) = avg_t(:,:, 37)
    rhoduzm   (:,:) = avg_t(:,:, 38)
    rhodvxm   (:,:) = avg_t(:,:, 39)
    rhodvym   (:,:) = avg_t(:,:, 40)
    rhodvzm   (:,:) = avg_t(:,:, 41)
    rhodwxm   (:,:) = avg_t(:,:, 42)
    rhodwym   (:,:) = avg_t(:,:, 43)
    rhodwzm   (:,:) = avg_t(:,:, 44)
    pdivm     (:,:) = avg_t(:,:, 45)
    rhodivm   (:,:) = avg_t(:,:, 46)
    vrtxm     (:,:) = avg_t(:,:, 47)
    vrtym     (:,:) = avg_t(:,:, 48)
    vrtzm     (:,:) = avg_t(:,:, 49)
    rhoTm     (:,:) = avg_t(:,:, 50)
    uTm       (:,:) = avg_t(:,:, 51)
    vTm       (:,:) = avg_t(:,:, 52)
    e2m       (:,:) = avg_t(:,:, 53)
    h2m       (:,:) = avg_t(:,:, 54)
    ccm       (:,:) = avg_t(:,:, 55)
    s2m       (:,:) = avg_t(:,:, 56)
    M2m       (:,:) = avg_t(:,:, 57)
    G2m       (:,:) = avg_t(:,:, 58)
    mu2m      (:,:) = avg_t(:,:, 59)
    la2m      (:,:) = avg_t(:,:, 60)
    cv2m      (:,:) = avg_t(:,:, 61)
    cp2m      (:,:) = avg_t(:,:, 62)
    pr2m      (:,:) = avg_t(:,:, 63)
    eck2m     (:,:) = avg_t(:,:, 64)
    upm       (:,:) = avg_t(:,:, 65)
    vpm       (:,:) = avg_t(:,:, 66)
    usm       (:,:) = avg_t(:,:, 67)
    vsm       (:,:) = avg_t(:,:, 68)
    rhopm     (:,:) = avg_t(:,:, 69)
    rhohm     (:,:) = avg_t(:,:, 70)     
    pTm       (:,:) = avg_t(:,:, 71)
    spm       (:,:) = avg_t(:,:, 72)
    sTm       (:,:) = avg_t(:,:, 73)
    rhosm     (:,:) = avg_t(:,:, 74)
    Grhom     (:,:) = avg_t(:,:, 75)
    Gpm       (:,:) = avg_t(:,:, 76)
    Gsm       (:,:) = avg_t(:,:, 77)
    GTm       (:,:) = avg_t(:,:, 78)
    Gum       (:,:) = avg_t(:,:, 79)
    Gvm       (:,:) = avg_t(:,:, 80)
    pduxm     (:,:) = avg_t(:,:, 81)
    pdvym     (:,:) = avg_t(:,:, 82)
    pdwzm     (:,:) = avg_t(:,:, 83)
    pduym     (:,:) = avg_t(:,:, 84)
    pdvxm     (:,:) = avg_t(:,:, 85)
    rhodiv2m  (:,:) = avg_t(:,:, 86)
    dux2m     (:,:) = avg_t(:,:, 87)
    duy2m     (:,:) = avg_t(:,:, 88)
    duz2m     (:,:) = avg_t(:,:, 89)
    dvx2m     (:,:) = avg_t(:,:, 90)
    dvy2m     (:,:) = avg_t(:,:, 91)
    dvz2m     (:,:) = avg_t(:,:, 92)
    dwx2m     (:,:) = avg_t(:,:, 93)
    dwy2m     (:,:) = avg_t(:,:, 94)
    dwz2m     (:,:) = avg_t(:,:, 95)
    vrtx2m    (:,:) = avg_t(:,:, 96)
    vrty2m    (:,:) = avg_t(:,:, 97)
    vrtz2m    (:,:) = avg_t(:,:, 98)     
    rhovrtxm  (:,:) = avg_t(:,:, 99)
    rhovrtym  (:,:) = avg_t(:,:,100)
    rhovrtzm  (:,:) = avg_t(:,:,101)
    rhouum    (:,:) = avg_t(:,:,102)
    rhovvm    (:,:) = avg_t(:,:,103)
    rhowwm    (:,:) = avg_t(:,:,104)
    rhoTTm    (:,:) = avg_t(:,:,105)
    rhovrtx2m (:,:) = avg_t(:,:,106)
    rhovrty2m (:,:) = avg_t(:,:,107)
    rhovrtz2m (:,:) = avg_t(:,:,108)
    rhouvm    (:,:) = avg_t(:,:,109)
    rhouwm    (:,:) = avg_t(:,:,110)
    rhovwm    (:,:) = avg_t(:,:,111)
    rhovTm    (:,:) = avg_t(:,:,112)     
    rhouuvm   (:,:) = avg_t(:,:,113)
    rhovvvm   (:,:) = avg_t(:,:,114)
    rhowwvm   (:,:) = avg_t(:,:,115)
    rhouvvm   (:,:) = avg_t(:,:,116)   
    rhodux2m  (:,:) = avg_t(:,:,117)
    rhodvy2m  (:,:) = avg_t(:,:,118)
    rhodwz2m  (:,:) = avg_t(:,:,119)
    rhoduydvxm(:,:) = avg_t(:,:,120)
    rhoduzdwxm(:,:) = avg_t(:,:,121)
    rhodvzdwym(:,:) = avg_t(:,:,122)
    u3m       (:,:) = avg_t(:,:,123)
    p3m       (:,:) = avg_t(:,:,124)
    u4m       (:,:) = avg_t(:,:,125)
    p4m       (:,:) = avg_t(:,:,126)
    tau11m    (:,:) =-avg_t(:,:,127)
    tau12m    (:,:) =-avg_t(:,:,128)
    tau13m    (:,:) =-avg_t(:,:,129)
    tau22m    (:,:) =-avg_t(:,:,130)
    tau23m    (:,:) =-avg_t(:,:,131)
    tau33m    (:,:) =-avg_t(:,:,132)
    tau12um   (:,:) =-avg_t(:,:,133)
    tau11um   (:,:) =-avg_t(:,:,134)
    tau12vm   (:,:) =-avg_t(:,:,135)
    tau13wm   (:,:) =-avg_t(:,:,136)
    tau22um   (:,:) =-avg_t(:,:,137)
    tau22vm   (:,:) =-avg_t(:,:,138)
    tau23wm   (:,:) =-avg_t(:,:,139)     
    tau11duxm (:,:) =-avg_t(:,:,140)
    tau11dvxm (:,:) =-avg_t(:,:,141)
    tau12duxm (:,:) =-avg_t(:,:,142)
    tau12duym (:,:) =-avg_t(:,:,143)
    tau12dvxm (:,:) =-avg_t(:,:,144)
    tau12dvym (:,:) =-avg_t(:,:,145)
    tau13duzm (:,:) =-avg_t(:,:,146)
    tau13dvzm (:,:) =-avg_t(:,:,147)
    tau13dwxm (:,:) =-avg_t(:,:,148)
    tau22duym (:,:) =-avg_t(:,:,149)
    tau22dvym (:,:) =-avg_t(:,:,150)
    tau23duzm (:,:) =-avg_t(:,:,151)
    tau23dvzm (:,:) =-avg_t(:,:,152)
    tau23dwym (:,:) =-avg_t(:,:,153)
    tau33dwzm (:,:) =-avg_t(:,:,154)
    ladTxm    (:,:) = avg_t(:,:,155)
    ladTym    (:,:) = avg_t(:,:,156)
    ladTzm    (:,:) = avg_t(:,:,157)
    hum       (:,:) = avg_t(:,:,158)
    hvm       (:,:) = avg_t(:,:,159)
    hwm       (:,:) = avg_t(:,:,160)
    rhohum    (:,:) = avg_t(:,:,161)
    rhohvm    (:,:) = avg_t(:,:,162)
    rhohwm    (:,:) = avg_t(:,:,163)
    rhouuum   (:,:) = avg_t(:,:,164)
    rhovvvm   (:,:) = avg_t(:,:,165)
    rhowwwm   (:,:) = avg_t(:,:,166)
    rhowwum   (:,:) = avg_t(:,:,167)

    deallocate(avg_t)

  end subroutine pp_stats_read_xy

end module mod_pp_stats_read
