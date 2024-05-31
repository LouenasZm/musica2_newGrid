!===============================================================================
subroutine stats_init
!===============================================================================
  !> Compute stats for TGV
!===============================================================================
  use mod_mpi
  use mod_constant
  implicit none
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  
  if (iproc.eq.0) then
     if ((.not.CHIT).and.(.not.CYL).and.(.not.ACT).and.(.not.TURB)) then
        open(98,file='stats.dat',form='formatted',status='replace') ! TO BE CHANGED (-> stats_init)
        close(98)
     endif

     if (CHIT) then
        open(199,file='stat1.dat',form='formatted',status='replace')
        rewind(199)
        open(100,file='es.dat',form='formatted',status='replace')
        rewind(100)
     elseif ((CYL).or.(TURB)) then
        if (idepart==1) then
           open(199,file='coeff.dat',form='formatted',status='replace')
           rewind(199)
        else
           open(199,file='coeff.dat',form='formatted',position='append')
        endif
     endif
     
  endif

end subroutine stats_init
  
!===============================================================================
subroutine stats_tgv
!===============================================================================
  !> Compute stats for TGV
!===============================================================================
  use mod_mpi
  use mod_constant
  use mod_flow
  use mod_time
  implicit none
  !-----------------------------------------------------------------------------
  ! Local variables
  integer  :: i,j,k
  logical  :: iexist
  real(wp) :: ekloc ,mgloc ,wrloc,ektot ,mgtot ,wrtot
  !-----------------------------------------------------------------------------

  ! Vorticity: compute positive azimuthal vorticity and enstrophy
  mgloc = 0.0_wp;   mgtot = 0.0_wp
  ekloc = 0.0_wp;   ektot = 0.0_wp
  wrloc = 0.0_wp;   wrtot = 0.0_wp
  do k=1,nz
     do j=1,ny
        do i=1,nx
           mgloc = mgloc + rho(i,j,k)

           ekloc = ekloc + 0.5_wp*rho(i,j,k)*( uu(i,j,k)**2 &
                + vv(i,j,k)**2 &
                + ww(i,j,k)**2 )

           wrloc = wrloc + 0.5_wp*rho(i,j,k)*( ( dwy(i,j,k) - dvz(i,j,k) )**2 &
                + ( duz(i,j,k) - dwx(i,j,k) )**2 &
                + ( dvx(i,j,k) - duy(i,j,k) )**2 )
        enddo
     enddo
  enddo

  !print *,mgloc,'for',iproc
  
  call MPI_REDUCE(mgloc,mgtot,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMM_global,info)
  call MPI_REDUCE(ekloc,ektot,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMM_global,info)
  call MPI_REDUCE(wrloc,wrtot,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMM_global,info)

  ! initial values
  if (abs(time).lt.tiny.or.ntime.eq.1) then
     mgtot0 = mgtot
     ektot0 = ektot
     wrtot0 = wrtot
  endif
  !
  ! ---------------------------------------------------------------------------
  ! Write to file and close
  ! -----------------------
  if(iproc.eq.0) then
     mgtot = mgtot/mgtot0
     ektot = ektot/(mgtot0*u_ref**2)
     wrtot = wrtot/(mgtot0/tscale**2)*2.0_wp/Re_ref

     inquire( file='stats.dat', exist=iexist )
     if (.not.iexist) then
        ! Create and write headers
        ! ------------------------
        open( unit=uni, file='stats.dat', status='replace' )
        write(uni,150) Re_ref, Mach, ngx, ngy, ngz
        write(uni,'(a)') 'VARIABLES = time Ek eps Mtot'
     else
        open( unit=uni, file='stats.dat', status='old', position='append' )
     endif

     ! Write data
     write(uni,'(4(1x, e19.10))') tstar, ektot, wrtot, mgtot
     close(uni)
  endif

150 format('TITLE="RE0=',f8.1,' MACH0=',f5.3,' grid:',i0,'x',i0,'x',i0,'"')

end subroutine stats_tgv

!===============================================================================
subroutine stats_chit
!===============================================================================
  !> Compute stats for CHIT
!===============================================================================
  use mpi
  use mod_constant
  use mod_mpi
  use mod_eos
  use mod_flow
  !use mod_artvisc
  use mod_time
  use mod_tranprop
  implicit none
  ! ----------------------------------------------------------------------------
  integer  :: i,j,k
  real(wp) :: cc2,qq,g,s_i,ttot_i,divloc,cpp,cvv,p_i,T_i,rho_i
  real(wp) :: u_i,v_i,w_i,b1,b2,b3,mu,la,prr,eck,e_i,h_i,vrtloc
  real(wp) :: bid,ming,maxg
  ! ----------------------------------------------------------------------------

  ! Initialization
  ! --------------
  avg_s  = 0.0_wp

  ! Collecting statistics
  do k=1,nz
     do j=1,ny
        do i=1,nx

           rho_i= rho(i,j,k)
           u_i = uu(i,j,k)
           v_i = vv(i,j,k)
           w_i = ww(i,j,k)
           T_i = Tmp(i,j,k)
           p_i = prs(i,j,k)
           cc2 = c_(i,j,k)**2

           g  = gcalc_tro(T_i,rho_i)
           cvv= cvcalc_tro(T_i,rho_i)
           cpp= cpcalc_tro(T_i,rho_i)
           s_i= scalc_tro(T_i,rho_i)

           qq = u_i**2 + v_i**2 + w_i**2
           ttot_i = T_i+0.5_wp*qq/cpp

           e_i = rhoe_n(i,j,k)/rho_i - 0.5_wp*qq
           h_i = e_i + p_i/rho_i

           b1 = dwy(i,j,k) - dvz(i,j,k)
           b2 = duz(i,j,k) - dwx(i,j,k)
           b3 = dvx(i,j,k) - duy(i,j,k)

           divloc = dux(i,j,k) + dvy(i,j,k) + dwz(i,j,k)
           vrtloc = sqrt(b1**2 + b2**2 + b3**2)

           mu = visc(i,j,k)
           la = cok(i,j,k)

           prr = cpp*mu/la
           eck= qq/(cpp*T_i)

           avg_s(1,1,  1)= avg_s(1,1,  1) + rho_i
           avg_s(1,1,  2)= avg_s(1,1,  2) + u_i
           avg_s(1,1,  3)= avg_s(1,1,  3) + v_i
           avg_s(1,1,  4)= avg_s(1,1,  4) + w_i
           avg_s(1,1,  5)= avg_s(1,1,  5) + p_i
           avg_s(1,1,  6)= avg_s(1,1,  6) + T_i
           avg_s(1,1,  7)= avg_s(1,1,  7) + e_i
           avg_s(1,1,  8)= avg_s(1,1,  8) + h_i
           avg_s(1,1,  9)= avg_s(1,1,  9) + sqrt(cc2)
           avg_s(1,1, 10)= avg_s(1,1, 10) + s_i
           avg_s(1,1, 11)= avg_s(1,1, 11) + sqrt(qq/cc2)
           avg_s(1,1, 12)= avg_s(1,1, 12) + 0.5*qq
           avg_s(1,1, 13)= avg_s(1,1, 13) + g
           avg_s(1,1, 14)= avg_s(1,1, 14) + mu
           avg_s(1,1, 15)= avg_s(1,1, 15) + la
           avg_s(1,1, 16)= avg_s(1,1, 16) + cpp
           avg_s(1,1, 17)= avg_s(1,1, 17) + cvv
           avg_s(1,1, 18)= avg_s(1,1, 18) + prr
           avg_s(1,1, 19)= avg_s(1,1, 19) + eck
           avg_s(1,1, 20)= avg_s(1,1, 20) + divloc
           avg_s(1,1, 21)= avg_s(1,1, 21) + rho_i*dux(i,j,k)
           avg_s(1,1, 22)= avg_s(1,1, 22) + rho_i*duy(i,j,k)
           avg_s(1,1, 23)= avg_s(1,1, 23) + rho_i*duz(i,j,k)
           avg_s(1,1, 24)= avg_s(1,1, 24) + rho_i*dvx(i,j,k)
           avg_s(1,1, 25)= avg_s(1,1, 25) + rho_i*dvy(i,j,k)
           avg_s(1,1, 26)= avg_s(1,1, 26) + rho_i*dvz(i,j,k)
           avg_s(1,1, 27)= avg_s(1,1, 27) + rho_i*dwx(i,j,k)
           avg_s(1,1, 28)= avg_s(1,1, 28) + rho_i*dwy(i,j,k)
           avg_s(1,1, 29)= avg_s(1,1, 29) + rho_i*dwz(i,j,k)
           avg_s(1,1, 30)= avg_s(1,1, 30) + p_i*divloc
           avg_s(1,1, 31)= avg_s(1,1, 31) + rho_i*divloc
           avg_s(1,1, 32)= avg_s(1,1, 32) + b1
           avg_s(1,1, 33)= avg_s(1,1, 33) + b2
           avg_s(1,1, 34)= avg_s(1,1, 34) + b3
           avg_s(1,1, 35)= avg_s(1,1, 35) + ttot_i
           avg_s(1,1, 36)= avg_s(1,1, 36) !+ ducros(i,j,k)
           avg_s(1,1, 37)= avg_s(1,1, 37) + rhou(i,j,k)
           avg_s(1,1, 38)= avg_s(1,1, 38) + rhov(i,j,k)
           avg_s(1,1, 39)= avg_s(1,1, 39) + rhow(i,j,k)
           avg_s(1,1, 40)= avg_s(1,1, 40) + rhoe(i,j,k)
           avg_s(1,1, 41)= avg_s(1,1, 41) + rho_i*T_i
           avg_s(1,1, 42)= avg_s(1,1, 42) + rho_i**2
           avg_s(1,1, 43)= avg_s(1,1, 43) + u_i**2
           avg_s(1,1, 44)= avg_s(1,1, 44) + v_i**2
           avg_s(1,1, 45)= avg_s(1,1, 45) + w_i**2
           avg_s(1,1, 46)= avg_s(1,1, 46) + u_i*v_i
           avg_s(1,1, 47)= avg_s(1,1, 47) + u_i*w_i
           avg_s(1,1, 48)= avg_s(1,1, 48) + v_i*w_i
           avg_s(1,1, 49)= avg_s(1,1, 49) + v_i*T_i
           avg_s(1,1, 50)= avg_s(1,1, 50) + p_i**2
           avg_s(1,1, 51)= avg_s(1,1, 51) + T_i**2
           avg_s(1,1, 52)= avg_s(1,1, 52) + e_i**2
           avg_s(1,1, 53)= avg_s(1,1, 53) + h_i**2
           avg_s(1,1, 54)= avg_s(1,1, 54) + cc2
           avg_s(1,1, 55)= avg_s(1,1, 55) + s_i**2
           avg_s(1,1, 56)= avg_s(1,1, 56) + qq/cc2
           avg_s(1,1, 57)= avg_s(1,1, 57) + g**2
           avg_s(1,1, 58)= avg_s(1,1, 58) + mu**2
           avg_s(1,1, 59)= avg_s(1,1, 59) + la**2
           avg_s(1,1, 60)= avg_s(1,1, 60) + cvv**2
           avg_s(1,1, 61)= avg_s(1,1, 61) + cpp**2
           avg_s(1,1, 62)= avg_s(1,1, 62) + prr**2
           avg_s(1,1, 63)= avg_s(1,1, 63) + eck**2
           avg_s(1,1, 64)= avg_s(1,1, 64) + p_i*u_i
           avg_s(1,1, 65)= avg_s(1,1, 65) + p_i*v_i
           avg_s(1,1, 66)= avg_s(1,1, 66) + T_i*u_i
           avg_s(1,1, 67)= avg_s(1,1, 67) + T_i*v_i
           avg_s(1,1, 68)= avg_s(1,1, 68) + s_i*u_i
           avg_s(1,1, 69)= avg_s(1,1, 69) + s_i*v_i
           avg_s(1,1, 70)= avg_s(1,1, 70) + p_i*rho_i
           avg_s(1,1, 71)= avg_s(1,1, 71) + T_i*rho_i
           avg_s(1,1, 72)= avg_s(1,1, 72) + h_i*rho_i
           avg_s(1,1, 73)= avg_s(1,1, 73) + T_i*p_i
           avg_s(1,1, 74)= avg_s(1,1, 74) + p_i*s_i
           avg_s(1,1, 75)= avg_s(1,1, 75) + T_i*s_i
           avg_s(1,1, 76)= avg_s(1,1, 76) + rho_i*s_i
           avg_s(1,1, 77)= avg_s(1,1, 77) + g*rho_i
           avg_s(1,1, 78)= avg_s(1,1, 78) + g*p_i
           avg_s(1,1, 79)= avg_s(1,1, 79) + g*s_i
           avg_s(1,1, 80)= avg_s(1,1, 80) + g*T_i
           avg_s(1,1, 81)= avg_s(1,1, 81) + g*u_i
           avg_s(1,1, 82)= avg_s(1,1, 82) + g*v_i
           avg_s(1,1, 83)= avg_s(1,1, 83) + p_i*dux(i,j,k)
           avg_s(1,1, 84)= avg_s(1,1, 84) + p_i*dvy(i,j,k)
           avg_s(1,1, 85)= avg_s(1,1, 85) + p_i*dwz(i,j,k)
           avg_s(1,1, 86)= avg_s(1,1, 86) + p_i*duy(i,j,k)
           avg_s(1,1, 87)= avg_s(1,1, 87) + p_i*dvx(i,j,k)
           avg_s(1,1, 88)= avg_s(1,1, 88) + divloc**2
           avg_s(1,1, 89)= avg_s(1,1, 89) + rho_i*divloc**2
           avg_s(1,1, 90)= avg_s(1,1, 90) + dux(i,j,k)**2
           avg_s(1,1, 91)= avg_s(1,1, 91) + duy(i,j,k)**2
           avg_s(1,1, 92)= avg_s(1,1, 92) + duz(i,j,k)**2
           avg_s(1,1, 93)= avg_s(1,1, 93) + dvx(i,j,k)**2
           avg_s(1,1, 94)= avg_s(1,1, 94) + dvy(i,j,k)**2
           avg_s(1,1, 95)= avg_s(1,1, 95) + dvz(i,j,k)**2
           avg_s(1,1, 96)= avg_s(1,1, 96) + dwx(i,j,k)**2
           avg_s(1,1, 97)= avg_s(1,1, 97) + dwy(i,j,k)**2
           avg_s(1,1, 98)= avg_s(1,1, 98) + dwz(i,j,k)**2
           avg_s(1,1, 99)= avg_s(1,1, 99) + b1**2
           avg_s(1,1,100)= avg_s(1,1,100) + b2**2
           avg_s(1,1,101)= avg_s(1,1,101) + b3**2
           avg_s(1,1,102)= avg_s(1,1,102) + rho_i*b1
           avg_s(1,1,103)= avg_s(1,1,103) + rho_i*b2
           avg_s(1,1,104)= avg_s(1,1,104) + rho_i*b3
           avg_s(1,1,105)= avg_s(1,1,105) + ttot_i**2
           avg_s(1,1,106)= avg_s(1,1,106) !+ ducros(i,j,k)**2
           avg_s(1,1,107)= avg_s(1,1,107) + rho_i*u_i**2
           avg_s(1,1,108)= avg_s(1,1,108) + rho_i*v_i**2
           avg_s(1,1,109)= avg_s(1,1,109) + rho_i*w_i**2
           avg_s(1,1,110)= avg_s(1,1,110) + rho_i*T_i**2
           avg_s(1,1,111)= avg_s(1,1,111) + rho_i*b1**2
           avg_s(1,1,112)= avg_s(1,1,112) + rho_i*b2**2
           avg_s(1,1,113)= avg_s(1,1,113) + rho_i*b3**2
           avg_s(1,1,114)= avg_s(1,1,114) + rho_i*u_i*v_i
           avg_s(1,1,115)= avg_s(1,1,115) + rho_i*v_i*w_i
           avg_s(1,1,116)= avg_s(1,1,116) + rho_i*v_i*T_i
           avg_s(1,1,117)= avg_s(1,1,117) + rho_i*u_i**2*v_i
           avg_s(1,1,118)= avg_s(1,1,118) + rho_i*v_i**2*v_i
           avg_s(1,1,119)= avg_s(1,1,119) + rho_i*w_i**2*v_i
           avg_s(1,1,120)= avg_s(1,1,120) + rho_i*v_i**2*u_i
           avg_s(1,1,121)= avg_s(1,1,121) + rho_i*dux(i,j,k)**2
           avg_s(1,1,122)= avg_s(1,1,122) + rho_i*dvy(i,j,k)**2
           avg_s(1,1,123)= avg_s(1,1,123) + rho_i*dwz(i,j,k)**2
           avg_s(1,1,124)= avg_s(1,1,124) + rho_i*duy(i,j,k)*dvx(i,j,k)
           avg_s(1,1,125)= avg_s(1,1,125) + rho_i*duz(i,j,k)*dwx(i,j,k)
           avg_s(1,1,126)= avg_s(1,1,126) + rho_i*dvz(i,j,k)*dwy(i,j,k)
           avg_s(1,1,127)= avg_s(1,1,127) + u_i**3
           avg_s(1,1,128)= avg_s(1,1,128) + p_i**3
           avg_s(1,1,129)= avg_s(1,1,129) + u_i**4
           avg_s(1,1,130)= avg_s(1,1,130) + p_i**4
           avg_s(1,1,131)= avg_s(1,1,131) + Frhou(i,j,k)
           avg_s(1,1,132)= avg_s(1,1,132) + Frhov(i,j,k)
           avg_s(1,1,133)= avg_s(1,1,133) + Frhow(i,j,k)
           avg_s(1,1,134)= avg_s(1,1,134) + Grhov(i,j,k)
           avg_s(1,1,135)= avg_s(1,1,135) + Grhow(i,j,k)
           avg_s(1,1,136)= avg_s(1,1,136) + Hrhow(i,j,k)
           avg_s(1,1,137)= avg_s(1,1,137) + Frhov(i,j,k)*u_i
           avg_s(1,1,138)= avg_s(1,1,138) + Frhov(i,j,k)*v_i
           avg_s(1,1,139)= avg_s(1,1,139) + Grhov(i,j,k)*u_i
           avg_s(1,1,140)= avg_s(1,1,140) + Grhov(i,j,k)*v_i
           avg_s(1,1,141)= avg_s(1,1,141) + Grhow(i,j,k)*w_i
           avg_s(1,1,142)= avg_s(1,1,142) + Frhou(i,j,k)*dux(i,j,k)
           avg_s(1,1,143)= avg_s(1,1,143) + Frhou(i,j,k)*dvx(i,j,k)
           avg_s(1,1,144)= avg_s(1,1,144) + Frhov(i,j,k)*dux(i,j,k)
           avg_s(1,1,145)= avg_s(1,1,145) + Frhov(i,j,k)*duy(i,j,k)
           avg_s(1,1,146)= avg_s(1,1,146) + Frhov(i,j,k)*dvx(i,j,k)
           avg_s(1,1,147)= avg_s(1,1,147) + Frhov(i,j,k)*dvy(i,j,k)
           avg_s(1,1,148)= avg_s(1,1,148) + Frhow(i,j,k)*duz(i,j,k)
           avg_s(1,1,149)= avg_s(1,1,149) + Frhow(i,j,k)*dvz(i,j,k)
           avg_s(1,1,150)= avg_s(1,1,150) + Frhow(i,j,k)*dwx(i,j,k)
           avg_s(1,1,151)= avg_s(1,1,151) + Grhov(i,j,k)*duy(i,j,k)
           avg_s(1,1,152)= avg_s(1,1,152) + Grhov(i,j,k)*dvy(i,j,k)
           avg_s(1,1,153)= avg_s(1,1,153) + Grhow(i,j,k)*duz(i,j,k)
           avg_s(1,1,154)= avg_s(1,1,154) + Grhow(i,j,k)*dvz(i,j,k)
           avg_s(1,1,155)= avg_s(1,1,155) + Grhow(i,j,k)*dwy(i,j,k)
           avg_s(1,1,156)= avg_s(1,1,156) + Hrhow(i,j,k)*dwz(i,j,k)
           avg_s(1,1,157)= avg_s(1,1,157) + la*dTx(i,j,k)
           avg_s(1,1,158)= avg_s(1,1,158) + la*dTy(i,j,k)
           avg_s(1,1,159)= avg_s(1,1,159) + la*dTz(i,j,k)
           avg_s(1,1,160)= avg_s(1,1,160) + h_i*u_i
           avg_s(1,1,161)= avg_s(1,1,161) + h_i*v_i
           avg_s(1,1,162)= avg_s(1,1,162) + h_i*w_i
           avg_s(1,1,163)= avg_s(1,1,163) + rho_i*h_i*u_i
           avg_s(1,1,164)= avg_s(1,1,164) + rho_i*h_i*v_i
           avg_s(1,1,165)= avg_s(1,1,165) + rho_i*h_i*w_i
           avg_s(1,1,166)= avg_s(1,1,166) + rho_i*u_i**3
           avg_s(1,1,167)= avg_s(1,1,167) + rho_i*v_i**3
           avg_s(1,1,168)= avg_s(1,1,168) + rho_i*w_i**3
           avg_s(1,1,169)= avg_s(1,1,169) + 0.!!divsens(i,j,k)
           avg_s(1,1,170)= avg_s(1,1,170) + dux(i,j,k)**3
           avg_s(1,1,171)= avg_s(1,1,171) + duy(i,j,k)**3
           avg_s(1,1,172)= avg_s(1,1,172) + duz(i,j,k)**3
           avg_s(1,1,173)= avg_s(1,1,173) + dvx(i,j,k)**3
           avg_s(1,1,174)= avg_s(1,1,174) + dvy(i,j,k)**3
           avg_s(1,1,175)= avg_s(1,1,175) + dvz(i,j,k)**3
           avg_s(1,1,176)= avg_s(1,1,176) + dwx(i,j,k)**3
           avg_s(1,1,177)= avg_s(1,1,177) + dwy(i,j,k)**3
           avg_s(1,1,178)= avg_s(1,1,178) + dwz(i,j,k)**3
           avg_s(1,1,179)= avg_s(1,1,179) + dux(i,j,k)**4
           avg_s(1,1,180)= avg_s(1,1,180) + duy(i,j,k)**4
           avg_s(1,1,181)= avg_s(1,1,181) + duz(i,j,k)**4
           avg_s(1,1,182)= avg_s(1,1,182) + dvx(i,j,k)**4
           avg_s(1,1,183)= avg_s(1,1,183) + dvy(i,j,k)**4
           avg_s(1,1,184)= avg_s(1,1,184) + dvz(i,j,k)**4
           avg_s(1,1,185)= avg_s(1,1,185) + dwx(i,j,k)**4
           avg_s(1,1,186)= avg_s(1,1,186) + dwy(i,j,k)**4
           avg_s(1,1,187)= avg_s(1,1,187) + dwz(i,j,k)**4
           avg_s(1,1,188)= avg_s(1,1,188) + mu/rho_i
           avg_s(1,1,189)= avg_s(1,1,189) + 4.0_wp/3.0_wp*mu*divloc**2
           avg_s(1,1,190)= avg_s(1,1,190) + mu*vrtloc*vrtloc
        enddo
     enddo
  enddo

  call MPI_ALLREDUCE(MPI_IN_PLACE,avg_s,200,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)

  avg_s = avg_s/dble(ngx*ngy*ngz)

  ! Find minimum and maximum density value
  ! --------------------------------------
  bid = minval(rho(1:nx,1:ny,1:nz))
  call MPI_REDUCE(bid,ming,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,COMM_global,info)
  bid = maxval(rho(1:nx,1:ny,1:nz))
  call MPI_REDUCE(bid,maxg,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,COMM_global,info)

  ! Write statistics to file
  ! ------------------------
  if (iproc.eq.0) then
     avg_s(1,1,191) = maxg/ming
     open(123,file='stats.dat',form='formatted',status='unknown',position='append')
     write(123,'(201(1pE19.12,1X))') tstar, (avg_s(1,1,i), i=1,size(avg_s))
     close(123)
  endif

end subroutine stats_chit

!===============================================================================
subroutine stats_chan
!===============================================================================
  !> Compute stats for channel flow
!===============================================================================
  use mod_mpi
  use mod_eos
  use mod_flow
  use mod_time
  use mod_constant
  use mod_rans
  !use mod_artvisc
  implicit none
  ! ----------------------------------------------------------------------------
  integer  :: i,j,k
  real(wp) :: nn1,inn1,cc2,qq,cnu13
  real(wp) :: g,s_i,ttot_i,divloc,cp,cv,p_i,T_i,rho_i,nu
  real(wp) :: u_i,v_i,w_i,b1,b2,b3,mu,la,prr,eck,e_i,h_i
  ! ----------------------------------------------------------------------------

  ! Collect statistics
  ! ------------------
  if (ntotal.lt.ndeb) return
  nn1 = dble(int((ntotal-ndeb)/freq_stats)) + 1
  inn1= 1.0_wp/nn1

  ! Spatial local and global averages at this iteration
  avg_s  = 0.0_wp

  j0=coord(2)*ny

  cnu1=7.1_wp
  cnu13=cnu1**3

  ! tau_ij = 2mu*s_ij - 2/3*mu*s_kk*delta_ij
  do k=1,nz
     do j=1,ny
        do i=1,nx
           rho_i= rho(i,j,k)
           u_i = uu(i,j,k)
           v_i = vv(i,j,k)
           w_i = ww(i,j,k)
           T_i = Tmp(i,j,k)
           p_i = prs(i,j,k)

           g = gcalc_tro(T_i,rho_i)
           cc2= c2calc_tro(T_i,rho_i)
           cv= cvcalc_tro(T_i,rho_i)
           cp= cpcalc_tro(T_i,rho_i)
           s_i = scalc_tro(T_i,rho_i)

           qq = u_i**2 + v_i**2 + w_i**2
           ttot_i = T_i+0.5_wp*qq/cp

           e_i = rhoe_n(i,j,k)/rho_i - 0.5_wp*qq
           h_i = e_i + p_i/rho_i

           b1 = dwy(i,j,k) - dvz(i,j,k)
           b2 = duz(i,j,k) - dwx(i,j,k)
           b3 = dvx(i,j,k) - duy(i,j,k)

           divloc = dux(i,j,k) + dvy(i,j,k) + dwz(i,j,k)

           mu = visc(i,j,k)
           la = cok(i,j,k)

           prr = cp*mu/la
           eck= qq/(cp*T_i)

           ! RANS statistics
           if (is_RANS) then
              avg_s(1,j,  1)= avg_s(1,j,  1) + rho_i
              avg_s(1,j,  2)= avg_s(1,j,  2) + u_i
              avg_s(1,j,  3)= avg_s(1,j,  3) + v_i
              avg_s(1,j,  4)= avg_s(1,j,  4) + p_i
              avg_s(1,j,  5)= avg_s(1,j,  5) + T_i
              avg_s(1,j,  6)= avg_s(1,j,  6) + e_i
              avg_s(1,j,  7)= avg_s(1,j,  7) + h_i
              avg_s(1,j,  8)= avg_s(1,j,  8) + sqrt(cc2)
              avg_s(1,j,  9)= avg_s(1,j,  9) + s_i
              avg_s(1,j, 10)= avg_s(1,j, 10) + sqrt(qq/cc2)
              avg_s(1,j, 11)= avg_s(1,j, 11) + 0.5*qq
              avg_s(1,j, 12)= avg_s(1,j, 12) + g
              avg_s(1,j, 13)= avg_s(1,j, 13) + mu
              avg_s(1,j, 14)= avg_s(1,j, 14) + la
              avg_s(1,j, 15)= avg_s(1,j, 15) + cp
              avg_s(1,j, 16)= avg_s(1,j, 16) + cv
              avg_s(1,j, 17)= avg_s(1,j, 17) + prr
              avg_s(1,j, 18)= avg_s(1,j, 18) + divloc
              avg_s(1,j, 19)= avg_s(1,j, 19) + mu*dux(i,j,k)
              avg_s(1,j, 20)= avg_s(1,j, 20) + mu*duy(i,j,k)
              avg_s(1,j, 21)= avg_s(1,j, 21) + mu*dvx(i,j,k)
              avg_s(1,j, 22)= avg_s(1,j, 22) + mu*dvy(i,j,k)
              avg_s(1,j, 23)= avg_s(1,j, 23) + dux(i,j,k)
              avg_s(1,j, 24)= avg_s(1,j, 24) + duy(i,j,k)
              avg_s(1,j, 25)= avg_s(1,j, 25) + dvx(i,j,k)
              avg_s(1,j, 26)= avg_s(1,j, 26) + dvy(i,j,k)
              avg_s(1,j, 27)= avg_s(1,j, 27) + dux(i,j,k)**2
              avg_s(1,j, 28)= avg_s(1,j, 28) + duy(i,j,k)**2
              avg_s(1,j, 29)= avg_s(1,j, 29) + dvx(i,j,k)**2
              avg_s(1,j, 30)= avg_s(1,j, 30) + dvy(i,j,k)**2
              avg_s(1,j, 31)= avg_s(1,j, 31) + u_i**2
              avg_s(1,j, 32)= avg_s(1,j, 32) + v_i**2
              avg_s(1,j, 33)= avg_s(1,j, 33) + u_i*v_i
              nu  = visc(i,j,k)/rho(i,j,k)
              khi = nutil(i,j,k)/nu
              fnu1 = khi**3/(khi**3+cnu13)
              mut(i,j,k) = nutil(i,j,k)*fnu1*rho_i
              avg_s(1,j, 34)= avg_s(1,j, 34) + mut(i,j,k)
              avg_s(1,j, 35)= avg_s(1,j, 35) + mut(i,j,k)*dux(i,j,k)
              avg_s(1,j, 36)= avg_s(1,j, 36) + mut(i,j,k)*duy(i,j,k)
              avg_s(1,j, 37)= avg_s(1,j, 37) + mut(i,j,k)*dvx(i,j,k)
              avg_s(1,j, 38)= avg_s(1,j, 38) + mut(i,j,k)*dvy(i,j,k)
           else
              avg_s(1,j,  1)= avg_s(1,j,  1) + rho_i
              avg_s(1,j,  2)= avg_s(1,j,  2) + u_i
              avg_s(1,j,  3)= avg_s(1,j,  3) + v_i
              avg_s(1,j,  4)= avg_s(1,j,  4) + w_i
              avg_s(1,j,  5)= avg_s(1,j,  5) + p_i
              avg_s(1,j,  6)= avg_s(1,j,  6) + T_i
              avg_s(1,j,  7)= avg_s(1,j,  7) + e_i
              avg_s(1,j,  8)= avg_s(1,j,  8) + h_i
              avg_s(1,j,  9)= avg_s(1,j,  9) + sqrt(cc2)
              avg_s(1,j, 10)= avg_s(1,j, 10) + s_i
              avg_s(1,j, 11)= avg_s(1,j, 11) + sqrt(qq/cc2)
              avg_s(1,j, 12)= avg_s(1,j, 12) + 0.5*qq
              avg_s(1,j, 13)= avg_s(1,j, 13) + g
              avg_s(1,j, 14)= avg_s(1,j, 14) + mu
              avg_s(1,j, 15)= avg_s(1,j, 15) + la
              avg_s(1,j, 16)= avg_s(1,j, 16) + cp
              avg_s(1,j, 17)= avg_s(1,j, 17) + cv
              avg_s(1,j, 18)= avg_s(1,j, 18) + prr
              avg_s(1,j, 19)= avg_s(1,j, 19) + eck
              avg_s(1,j, 20)= avg_s(1,j, 20) + divloc
              avg_s(1,j, 21)= avg_s(1,j, 21) + rho_i*dux(i,j,k)
              avg_s(1,j, 22)= avg_s(1,j, 22) + rho_i*duy(i,j,k)
              avg_s(1,j, 23)= avg_s(1,j, 23) + rho_i*duz(i,j,k)
              avg_s(1,j, 24)= avg_s(1,j, 24) + rho_i*dvx(i,j,k)
              avg_s(1,j, 25)= avg_s(1,j, 25) + rho_i*dvy(i,j,k)
              avg_s(1,j, 26)= avg_s(1,j, 26) + rho_i*dvz(i,j,k)
              avg_s(1,j, 27)= avg_s(1,j, 27) + rho_i*dwx(i,j,k)
              avg_s(1,j, 28)= avg_s(1,j, 28) + rho_i*dwy(i,j,k)
              avg_s(1,j, 29)= avg_s(1,j, 29) + rho_i*dwz(i,j,k)
              avg_s(1,j, 30)= avg_s(1,j, 30) + p_i*divloc
              avg_s(1,j, 31)= avg_s(1,j, 31) + rho_i*divloc
              avg_s(1,j, 32)= avg_s(1,j, 32) + b1
              avg_s(1,j, 33)= avg_s(1,j, 33) + b2
              avg_s(1,j, 34)= avg_s(1,j, 34) + b3
              avg_s(1,j, 35)= avg_s(1,j, 35) + ttot_i
              avg_s(1,j, 36)= avg_s(1,j, 36) + 0.!ducros(i,j,k)
              avg_s(1,j, 37)= avg_s(1,j, 37) + rhou(i,j,k)
              avg_s(1,j, 38)= avg_s(1,j, 38) + rhov(i,j,k)
              avg_s(1,j, 39)= avg_s(1,j, 39) + rhow(i,j,k)
              avg_s(1,j, 40)= avg_s(1,j, 40) + rhoe(i,j,k)
              avg_s(1,j, 41)= avg_s(1,j, 41) + rho_i*T_i
              avg_s(1,j, 42)= avg_s(1,j, 42) + rho_i**2
              avg_s(1,j, 43)= avg_s(1,j, 43) + u_i**2
              avg_s(1,j, 44)= avg_s(1,j, 44) + v_i**2
              avg_s(1,j, 45)= avg_s(1,j, 45) + w_i**2
              avg_s(1,j, 46)= avg_s(1,j, 46) + u_i*v_i
              avg_s(1,j, 47)= avg_s(1,j, 47) + u_i*w_i
              avg_s(1,j, 48)= avg_s(1,j, 48) + v_i*w_i
              avg_s(1,j, 49)= avg_s(1,j, 49) + v_i*T_i
              avg_s(1,j, 50)= avg_s(1,j, 50) + p_i**2
              avg_s(1,j, 51)= avg_s(1,j, 51) + T_i**2
              avg_s(1,j, 52)= avg_s(1,j, 52) + e_i**2
              avg_s(1,j, 53)= avg_s(1,j, 53) + h_i**2
              avg_s(1,j, 54)= avg_s(1,j, 54) + cc2
              avg_s(1,j, 55)= avg_s(1,j, 55) + s_i**2
              avg_s(1,j, 56)= avg_s(1,j, 56) + qq/cc2
              avg_s(1,j, 57)= avg_s(1,j, 57) + g**2
              avg_s(1,j, 58)= avg_s(1,j, 58) + mu**2
              avg_s(1,j, 59)= avg_s(1,j, 59) + la**2
              avg_s(1,j, 60)= avg_s(1,j, 60) + cv**2
              avg_s(1,j, 61)= avg_s(1,j, 61) + cp**2
              avg_s(1,j, 62)= avg_s(1,j, 62) + prr**2
              avg_s(1,j, 63)= avg_s(1,j, 63) + eck**2
              avg_s(1,j, 64)= avg_s(1,j, 64) + p_i*u_i
              avg_s(1,j, 65)= avg_s(1,j, 65) + p_i*v_i
              avg_s(1,j, 66)= avg_s(1,j, 66) + T_i*u_i
              avg_s(1,j, 67)= avg_s(1,j, 67) + T_i*v_i
              avg_s(1,j, 68)= avg_s(1,j, 68) + s_i*u_i
              avg_s(1,j, 69)= avg_s(1,j, 69) + s_i*v_i
              avg_s(1,j, 70)= avg_s(1,j, 70) + p_i*rho_i
              avg_s(1,j, 71)= avg_s(1,j, 71) + T_i*rho_i
              avg_s(1,j, 72)= avg_s(1,j, 72) + h_i*rho_i
              avg_s(1,j, 73)= avg_s(1,j, 73) + T_i*p_i
              avg_s(1,j, 74)= avg_s(1,j, 74) + p_i*s_i
              avg_s(1,j, 75)= avg_s(1,j, 75) + T_i*s_i
              avg_s(1,j, 76)= avg_s(1,j, 76) + rho_i*s_i
              avg_s(1,j, 77)= avg_s(1,j, 77) + g*rho_i
              avg_s(1,j, 78)= avg_s(1,j, 78) + g*p_i
              avg_s(1,j, 79)= avg_s(1,j, 79) + g*s_i
              avg_s(1,j, 80)= avg_s(1,j, 80) + g*T_i
              avg_s(1,j, 81)= avg_s(1,j, 81) + g*u_i
              avg_s(1,j, 82)= avg_s(1,j, 82) + g*v_i
              avg_s(1,j, 83)= avg_s(1,j, 83) + p_i*dux(i,j,k)
              avg_s(1,j, 84)= avg_s(1,j, 84) + p_i*dvy(i,j,k)
              avg_s(1,j, 85)= avg_s(1,j, 85) + p_i*dwz(i,j,k)
              avg_s(1,j, 86)= avg_s(1,j, 86) + p_i*duy(i,j,k)
              avg_s(1,j, 87)= avg_s(1,j, 87) + p_i*dvx(i,j,k)
              avg_s(1,j, 88)= avg_s(1,j, 88) + divloc**2
              avg_s(1,j, 89)= avg_s(1,j, 89) + rho_i*divloc**2
              avg_s(1,j, 90)= avg_s(1,j, 90) + dux(i,j,k)**2
              avg_s(1,j, 91)= avg_s(1,j, 91) + duy(i,j,k)**2
              avg_s(1,j, 92)= avg_s(1,j, 92) + duz(i,j,k)**2
              avg_s(1,j, 93)= avg_s(1,j, 93) + dvx(i,j,k)**2
              avg_s(1,j, 94)= avg_s(1,j, 94) + dvy(i,j,k)**2
              avg_s(1,j, 95)= avg_s(1,j, 95) + dvz(i,j,k)**2
              avg_s(1,j, 96)= avg_s(1,j, 96) + dwx(i,j,k)**2
              avg_s(1,j, 97)= avg_s(1,j, 97) + dwy(i,j,k)**2
              avg_s(1,j, 98)= avg_s(1,j, 98) + dwz(i,j,k)**2
              avg_s(1,j, 99)= avg_s(1,j, 99) + b1**2
              avg_s(1,j,100)= avg_s(1,j,100) + b2**2
              avg_s(1,j,101)= avg_s(1,j,101) + b3**2
              avg_s(1,j,102)= avg_s(1,j,102) + rho_i*b1
              avg_s(1,j,103)= avg_s(1,j,103) + rho_i*b2
              avg_s(1,j,104)= avg_s(1,j,104) + rho_i*b3
              avg_s(1,j,105)= avg_s(1,j,105) + ttot_i**2
              avg_s(1,j,106)= avg_s(1,j,106) + 0.!ducros(i,j,k)**2
              avg_s(1,j,107)= avg_s(1,j,107) + rho_i*u_i**2
              avg_s(1,j,108)= avg_s(1,j,108) + rho_i*v_i**2
              avg_s(1,j,109)= avg_s(1,j,109) + rho_i*w_i**2
              avg_s(1,j,110)= avg_s(1,j,110) + rho_i*T_i**2
              avg_s(1,j,111)= avg_s(1,j,111) + rho_i*b1**2
              avg_s(1,j,112)= avg_s(1,j,112) + rho_i*b2**2
              avg_s(1,j,113)= avg_s(1,j,113) + rho_i*b3**2
              avg_s(1,j,114)= avg_s(1,j,114) + rho_i*u_i*v_i
              avg_s(1,j,115)= avg_s(1,j,115) + rho_i*v_i*w_i
              avg_s(1,j,116)= avg_s(1,j,116) + rho_i*v_i*T_i
              avg_s(1,j,117)= avg_s(1,j,117) + rho_i*u_i**2*v_i
              avg_s(1,j,118)= avg_s(1,j,118) + rho_i*v_i**2*v_i
              avg_s(1,j,119)= avg_s(1,j,119) + rho_i*w_i**2*v_i
              avg_s(1,j,120)= avg_s(1,j,120) + rho_i*v_i**2*u_i
              avg_s(1,j,121)= avg_s(1,j,121) + rho_i*dux(i,j,k)**2
              avg_s(1,j,122)= avg_s(1,j,122) + rho_i*dvy(i,j,k)**2
              avg_s(1,j,123)= avg_s(1,j,123) + rho_i*dwz(i,j,k)**2
              avg_s(1,j,124)= avg_s(1,j,124) + rho_i*duy(i,j,k)*dvx(i,j,k)
              avg_s(1,j,125)= avg_s(1,j,125) + rho_i*duz(i,j,k)*dwx(i,j,k)
              avg_s(1,j,126)= avg_s(1,j,126) + rho_i*dvz(i,j,k)*dwy(i,j,k)
              avg_s(1,j,127)= avg_s(1,j,127) + u_i**3
              avg_s(1,j,128)= avg_s(1,j,128) + p_i**3
              avg_s(1,j,129)= avg_s(1,j,129) + u_i**4
              avg_s(1,j,130)= avg_s(1,j,130) + p_i**4
              avg_s(1,j,131)= avg_s(1,j,131) + Frhou(i,j,k)
              avg_s(1,j,132)= avg_s(1,j,132) + Frhov(i,j,k)
              avg_s(1,j,133)= avg_s(1,j,133) + Frhow(i,j,k)
              avg_s(1,j,134)= avg_s(1,j,134) + Grhov(i,j,k)
              avg_s(1,j,135)= avg_s(1,j,135) + Grhow(i,j,k)
              avg_s(1,j,136)= avg_s(1,j,136) + Hrhow(i,j,k)
              avg_s(1,j,137)= avg_s(1,j,137) + Frhov(i,j,k)*u_i
              avg_s(1,j,138)= avg_s(1,j,138) + Frhov(i,j,k)*v_i
              avg_s(1,j,139)= avg_s(1,j,139) + Grhov(i,j,k)*u_i
              avg_s(1,j,140)= avg_s(1,j,140) + Grhov(i,j,k)*v_i
              avg_s(1,j,141)= avg_s(1,j,141) + Grhow(i,j,k)*w_i
              avg_s(1,j,142)= avg_s(1,j,142) + Frhou(i,j,k)*dux(i,j,k)
              avg_s(1,j,143)= avg_s(1,j,143) + Frhou(i,j,k)*dvx(i,j,k)
              avg_s(1,j,144)= avg_s(1,j,144) + Frhov(i,j,k)*dux(i,j,k)
              avg_s(1,j,145)= avg_s(1,j,145) + Frhov(i,j,k)*duy(i,j,k)
              avg_s(1,j,146)= avg_s(1,j,146) + Frhov(i,j,k)*dvx(i,j,k)
              avg_s(1,j,147)= avg_s(1,j,147) + Frhov(i,j,k)*dvy(i,j,k)
              avg_s(1,j,148)= avg_s(1,j,148) + Frhow(i,j,k)*duz(i,j,k)
              avg_s(1,j,149)= avg_s(1,j,149) + Frhow(i,j,k)*dvz(i,j,k)
              avg_s(1,j,150)= avg_s(1,j,150) + Frhow(i,j,k)*dwx(i,j,k)
              avg_s(1,j,151)= avg_s(1,j,151) + Grhov(i,j,k)*duy(i,j,k)
              avg_s(1,j,152)= avg_s(1,j,152) + Grhov(i,j,k)*dvy(i,j,k)
              avg_s(1,j,153)= avg_s(1,j,153) + Grhow(i,j,k)*duz(i,j,k)
              avg_s(1,j,154)= avg_s(1,j,154) + Grhow(i,j,k)*dvz(i,j,k)
              avg_s(1,j,155)= avg_s(1,j,155) + Grhow(i,j,k)*dwy(i,j,k)
              avg_s(1,j,156)= avg_s(1,j,156) + Hrhow(i,j,k)*dwz(i,j,k)
              avg_s(1,j,157)= avg_s(1,j,157) + la*dTx(i,j,k)
              avg_s(1,j,158)= avg_s(1,j,158) + la*dTy(i,j,k)
              avg_s(1,j,159)= avg_s(1,j,159) + la*dTz(i,j,k)
              avg_s(1,j,160)= avg_s(1,j,160) + h_i*u_i
              avg_s(1,j,161)= avg_s(1,j,161) + h_i*v_i
              avg_s(1,j,162)= avg_s(1,j,162) + h_i*w_i
              avg_s(1,j,163)= avg_s(1,j,163) + rho_i*h_i*u_i
              avg_s(1,j,164)= avg_s(1,j,164) + rho_i*h_i*v_i
              avg_s(1,j,165)= avg_s(1,j,165) + rho_i*h_i*w_i
              avg_s(1,j,166)= avg_s(1,j,166) + rho_i*u_i**3
              avg_s(1,j,167)= avg_s(1,j,167) + rho_i*v_i**3
              avg_s(1,j,168)= avg_s(1,j,168) + rho_i*w_i**3
              ! avg_s(1,j,169)= avg_s(1,j,169) + rho_i*forc_rhou
              ! avg_s(1,j,170)= avg_s(1,j,170) + Drhotot(i,j,k)
              ! avg_s(1,j,171)= avg_s(1,j,171) + Drhoutot(i,j,k)
              ! avg_s(1,j,172)= avg_s(1,j,172) + Drhovtot(i,j,k)
              ! avg_s(1,j,173)= avg_s(1,j,173) + Drhowtot(i,j,k)
              ! avg_s(1,j,174)= avg_s(1,j,174) + Drhoetot(i,j,k)
              ! avg_s(1,j,175)= avg_s(1,j,175) + u_i*Drhotot(i,j,k)
              ! avg_s(1,j,176)= avg_s(1,j,176) + v_i*Drhotot(i,j,k)
              ! avg_s(1,j,177)= avg_s(1,j,177) + w_i*Drhotot(i,j,k)
              ! avg_s(1,j,178)= avg_s(1,j,178) + u_i*Drhoutot(i,j,k)
              ! avg_s(1,j,179)= avg_s(1,j,179) + v_i*Drhovtot(i,j,k)
              ! avg_s(1,j,180)= avg_s(1,j,180) + w_i*Drhowtot(i,j,k)
           endif
        enddo
     enddo
  enddo

  avg_s = avg_s/dble(nx*nz)

  ! Updating time averaging
  ! -----------------------
  avg_t = ((nn1-1.0_wp)*avg_t + avg_s)*inn1

end subroutine stats_chan

!===============================================================================
subroutine stats_stbl
!===============================================================================
  !> Compute stats for boundary layers
!===============================================================================
  use mod_mpi
  use mod_eos
  use mod_flow
  use mod_time
  use mod_constant
  !use mod_artvisc
  implicit none
  ! ----------------------------------------------------------------------------
  integer  :: i,j,k
  real(wp) :: nn1,inn1,cc2,qq
  real(wp) :: g,s_i,divloc,cp,cv,p_i,T_i,rho_i
  real(wp) :: u_i,v_i,w_i,b1,b2,b3,mu,la,prr,eck,e_i,h_i
  ! ----------------------------------------------------------------------------

  ! Collect statistics
  ! ------------------
  if (ntotal.lt.ndeb) return
  nn1 = dble(int((ntotal-ndeb)/freq_stats)) + 1
  inn1= 1.0_wp/nn1

  ! Spatial local and global averages at this iteration
  avg_s  = 0.0_wp
 
  do k=1,nz
     do j=1,ny
        do i=1,nx
           rho_i= rho(i,j,k)
           u_i = uu(i,j,k)
           v_i = vv(i,j,k)
           w_i = ww(i,j,k)
           T_i = Tmp(i,j,k)
           p_i = prs(i,j,k)

           g = gcalc_tro(T_i,rho_i)
           cc2= c2calc_tro(T_i,rho_i)
           cv= cvcalc_tro(T_i,rho_i)
           cp= cpcalc_tro(T_i,rho_i)
           s_i = scalc_tro(T_i,rho_i)

           qq = u_i**2 + v_i**2 + w_i**2

           e_i = rhoe_n(i,j,k)/rho_i - 0.5_wp*qq
           h_i = e_i + p_i/rho_i

           b1 = dwy(i,j,k) - dvz(i,j,k)
           b2 = duz(i,j,k) - dwx(i,j,k)
           b3 = dvx(i,j,k) - duy(i,j,k)
           divloc = dux(i,j,k) + dvy(i,j,k) + dwz(i,j,k)

           mu = visc(i,j,k)
           la = cok(i,j,k)

           prr = cp*mu/la
           eck= qq/(cp*T_i)

           avg_s(i,j,  1)= avg_s(i,j,  1) + rho_i
           avg_s(i,j,  2)= avg_s(i,j,  2) + u_i
           avg_s(i,j,  3)= avg_s(i,j,  3) + v_i
           avg_s(i,j,  4)= avg_s(i,j,  4) + w_i
           avg_s(i,j,  5)= avg_s(i,j,  5) + p_i
           avg_s(i,j,  6)= avg_s(i,j,  6) + T_i
           avg_s(i,j,  7)= avg_s(i,j,  7) + rhou(i,j,k)
           avg_s(i,j,  8)= avg_s(i,j,  8) + rhov(i,j,k)
           avg_s(i,j,  9)= avg_s(i,j,  9) + rhow(i,j,k)
           avg_s(i,j, 10)= avg_s(i,j, 10) + rhoe(i,j,k)
           avg_s(i,j, 11)= avg_s(i,j, 11) + rho_i**2
           avg_s(i,j, 12)= avg_s(i,j, 12) + u_i**2
           avg_s(i,j, 13)= avg_s(i,j, 13) + v_i**2
           avg_s(i,j, 14)= avg_s(i,j, 14) + w_i**2
           avg_s(i,j, 15)= avg_s(i,j, 15) + u_i*v_i
           avg_s(i,j, 16)= avg_s(i,j, 16) + u_i*w_i
           avg_s(i,j, 17)= avg_s(i,j, 17) + v_i*w_i
           avg_s(i,j, 18)= avg_s(i,j, 18) + v_i*T_i
           avg_s(i,j, 19)= avg_s(i,j, 19) + p_i**2
           avg_s(i,j, 20)= avg_s(i,j, 20) + T_i**2
           avg_s(i,j, 21)= avg_s(i,j, 21) + mu
           avg_s(i,j, 22)= avg_s(i,j, 22) + divloc
           avg_s(i,j, 23)= avg_s(i,j, 23) + divloc**2

           avg_s(i,j, 24)= avg_s(i,j, 24) + e_i
           avg_s(i,j, 25)= avg_s(i,j, 25) + h_i
           avg_s(i,j, 26)= avg_s(i,j, 26) + sqrt(cc2)
           avg_s(i,j, 27)= avg_s(i,j, 27) + s_i
           avg_s(i,j, 28)= avg_s(i,j, 28) + sqrt(qq/cc2)
           avg_s(i,j, 29)= avg_s(i,j, 29) + 0.5*qq
           avg_s(i,j, 30)= avg_s(i,j, 30) + g
           avg_s(i,j, 31)= avg_s(i,j, 31) + la
           avg_s(i,j, 32)= avg_s(i,j, 32) + cp
           avg_s(i,j, 33)= avg_s(i,j, 33) + cv
           avg_s(i,j, 34)= avg_s(i,j, 34) + prr
           avg_s(i,j, 35)= avg_s(i,j, 35) + eck
           avg_s(i,j, 36)= avg_s(i,j, 36) + rho_i*dux(i,j,k)
           avg_s(i,j, 37)= avg_s(i,j, 37) + rho_i*duy(i,j,k)
           avg_s(i,j, 38)= avg_s(i,j, 38) + rho_i*duz(i,j,k)
           avg_s(i,j, 39)= avg_s(i,j, 39) + rho_i*dvx(i,j,k)
           avg_s(i,j, 40)= avg_s(i,j, 40) + rho_i*dvy(i,j,k)
           avg_s(i,j, 41)= avg_s(i,j, 41) + rho_i*dvz(i,j,k)
           avg_s(i,j, 42)= avg_s(i,j, 42) + rho_i*dwx(i,j,k)
           avg_s(i,j, 43)= avg_s(i,j, 43) + rho_i*dwy(i,j,k)
           avg_s(i,j, 44)= avg_s(i,j, 44) + rho_i*dwz(i,j,k)
           avg_s(i,j, 45)= avg_s(i,j, 45) + p_i*divloc
           avg_s(i,j, 46)= avg_s(i,j, 46) + rho_i*divloc
           avg_s(i,j, 47)= avg_s(i,j, 47) + b1
           avg_s(i,j, 48)= avg_s(i,j, 48) + b2
           avg_s(i,j, 49)= avg_s(i,j, 49) + b3
           avg_s(i,j, 50)= avg_s(i,j, 50) + rho_i*T_i
           avg_s(i,j, 51)= avg_s(i,j, 51) + u_i*T_i
           avg_s(i,j, 52)= avg_s(i,j, 52) + v_i*T_i
           avg_s(i,j, 53)= avg_s(i,j, 53) + e_i**2
           avg_s(i,j, 54)= avg_s(i,j, 54) + h_i**2
           avg_s(i,j, 55)= avg_s(i,j, 55) + cc2
           avg_s(i,j, 56)= avg_s(i,j, 56) + s_i**2
           avg_s(i,j, 57)= avg_s(i,j, 57) + qq/cc2
           avg_s(i,j, 58)= avg_s(i,j, 58) + g**2
           avg_s(i,j, 59)= avg_s(i,j, 59) + mu**2
           avg_s(i,j, 60)= avg_s(i,j, 60) + la**2
           avg_s(i,j, 61)= avg_s(i,j, 61) + cv**2
           avg_s(i,j, 62)= avg_s(i,j, 62) + cp**2
           avg_s(i,j, 63)= avg_s(i,j, 63) + prr**2
           avg_s(i,j, 64)= avg_s(i,j, 64) + eck**2
           avg_s(i,j, 65)= avg_s(i,j, 65) + p_i*u_i
           avg_s(i,j, 66)= avg_s(i,j, 66) + p_i*v_i
           avg_s(i,j, 67)= avg_s(i,j, 67) + s_i*u_i
           avg_s(i,j, 68)= avg_s(i,j, 68) + s_i*v_i
           avg_s(i,j, 69)= avg_s(i,j, 69) + p_i*rho_i
           avg_s(i,j, 70)= avg_s(i,j, 70) + h_i*rho_i
           avg_s(i,j, 71)= avg_s(i,j, 71) + T_i*p_i
           avg_s(i,j, 72)= avg_s(i,j, 72) + p_i*s_i
           avg_s(i,j, 73)= avg_s(i,j, 73) + T_i*s_i
           avg_s(i,j, 74)= avg_s(i,j, 74) + rho_i*s_i
           avg_s(i,j, 75)= avg_s(i,j, 75) + g*rho_i
           avg_s(i,j, 76)= avg_s(i,j, 76) + g*p_i
           avg_s(i,j, 77)= avg_s(i,j, 77) + g*s_i
           avg_s(i,j, 78)= avg_s(i,j, 78) + g*T_i
           avg_s(i,j, 79)= avg_s(i,j, 79) + g*u_i
           avg_s(i,j, 80)= avg_s(i,j, 80) + g*v_i
           avg_s(i,j, 81)= avg_s(i,j, 81) + p_i*dux(i,j,k)
           avg_s(i,j, 82)= avg_s(i,j, 82) + p_i*dvy(i,j,k)
           avg_s(i,j, 83)= avg_s(i,j, 83) + p_i*dwz(i,j,k)
           avg_s(i,j, 84)= avg_s(i,j, 84) + p_i*duy(i,j,k)
           avg_s(i,j, 85)= avg_s(i,j, 85) + p_i*dvx(i,j,k)
           avg_s(i,j, 86)= avg_s(i,j, 86) + rho_i*divloc**2
           avg_s(i,j, 87)= avg_s(i,j, 87) + dux(i,j,k)**2
           avg_s(i,j, 88)= avg_s(i,j, 88) + duy(i,j,k)**2
           avg_s(i,j, 89)= avg_s(i,j, 89) + duz(i,j,k)**2
           avg_s(i,j, 90)= avg_s(i,j, 90) + dvx(i,j,k)**2
           avg_s(i,j, 91)= avg_s(i,j, 91) + dvy(i,j,k)**2
           avg_s(i,j, 92)= avg_s(i,j, 92) + dvz(i,j,k)**2
           avg_s(i,j, 93)= avg_s(i,j, 93) + dwx(i,j,k)**2
           avg_s(i,j, 94)= avg_s(i,j, 94) + dwy(i,j,k)**2
           avg_s(i,j, 95)= avg_s(i,j, 95) + dwz(i,j,k)**2
           avg_s(i,j, 96)= avg_s(i,j, 96) + b1**2
           avg_s(i,j, 97)= avg_s(i,j, 97) + b2**2
           avg_s(i,j, 98)= avg_s(i,j, 98) + b3**2
           avg_s(i,j, 99)= avg_s(i,j, 99) + rho_i*b1
           avg_s(i,j,100)= avg_s(i,j,100) + rho_i*b2
           avg_s(i,j,101)= avg_s(i,j,101) + rho_i*b3
           avg_s(i,j,102)= avg_s(i,j,102) + rho_i*u_i**2
           avg_s(i,j,103)= avg_s(i,j,103) + rho_i*v_i**2
           avg_s(i,j,104)= avg_s(i,j,104) + rho_i*w_i**2
           avg_s(i,j,105)= avg_s(i,j,105) + rho_i*T_i**2
           avg_s(i,j,106)= avg_s(i,j,106) + rho_i*b1**2
           avg_s(i,j,107)= avg_s(i,j,107) + rho_i*b2**2
           avg_s(i,j,108)= avg_s(i,j,108) + rho_i*b3**2
           avg_s(i,j,109)= avg_s(i,j,109) + rho_i*u_i*v_i
           avg_s(i,j,110)= avg_s(i,j,110) + rho_i*u_i*w_i
           avg_s(i,j,111)= avg_s(i,j,111) + rho_i*v_i*w_i
           avg_s(i,j,112)= avg_s(i,j,112) + rho_i*v_i*T_i
           avg_s(i,j,113)= avg_s(i,j,113) + rho_i*u_i**2*v_i
           avg_s(i,j,114)= avg_s(i,j,114) + rho_i*v_i**2*v_i
           avg_s(i,j,115)= avg_s(i,j,115) + rho_i*w_i**2*v_i
           avg_s(i,j,116)= avg_s(i,j,116) + rho_i*v_i**2*u_i
           avg_s(i,j,117)= avg_s(i,j,117) + rho_i*dux(i,j,k)**2
           avg_s(i,j,118)= avg_s(i,j,118) + rho_i*dvy(i,j,k)**2
           avg_s(i,j,119)= avg_s(i,j,119) + rho_i*dwz(i,j,k)**2
           avg_s(i,j,120)= avg_s(i,j,120) + rho_i*duy(i,j,k)*dvx(i,j,k)
           avg_s(i,j,121)= avg_s(i,j,121) + rho_i*duz(i,j,k)*dwx(i,j,k)
           avg_s(i,j,122)= avg_s(i,j,122) + rho_i*dvz(i,j,k)*dwy(i,j,k)
           avg_s(i,j,123)= avg_s(i,j,123) + u_i**3
           avg_s(i,j,124)= avg_s(i,j,124) + p_i**3
           avg_s(i,j,125)= avg_s(i,j,125) + u_i**4
           avg_s(i,j,126)= avg_s(i,j,126) + p_i**4
           avg_s(i,j,127)= avg_s(i,j,127) + Frhou(i,j,k)
           avg_s(i,j,128)= avg_s(i,j,128) + Frhov(i,j,k)
           avg_s(i,j,129)= avg_s(i,j,129) + Frhow(i,j,k)
           avg_s(i,j,130)= avg_s(i,j,130) + Grhov(i,j,k)
           avg_s(i,j,131)= avg_s(i,j,131) + Grhow(i,j,k)
           avg_s(i,j,132)= avg_s(i,j,132) + Hrhow(i,j,k)
           avg_s(i,j,133)= avg_s(i,j,133) + Frhov(i,j,k)*u_i
           avg_s(i,j,134)= avg_s(i,j,134) + Frhou(i,j,k)*u_i
           avg_s(i,j,135)= avg_s(i,j,135) + Frhov(i,j,k)*v_i
           avg_s(i,j,136)= avg_s(i,j,136) + Frhow(i,j,k)*w_i
           avg_s(i,j,137)= avg_s(i,j,137) + Grhov(i,j,k)*u_i
           avg_s(i,j,138)= avg_s(i,j,138) + Grhov(i,j,k)*v_i
           avg_s(i,j,139)= avg_s(i,j,139) + Grhow(i,j,k)*w_i
           avg_s(i,j,140)= avg_s(i,j,140) + Frhou(i,j,k)*dux(i,j,k)
           avg_s(i,j,141)= avg_s(i,j,141) + Frhou(i,j,k)*dvx(i,j,k)
           avg_s(i,j,142)= avg_s(i,j,142) + Frhov(i,j,k)*dux(i,j,k)
           avg_s(i,j,143)= avg_s(i,j,143) + Frhov(i,j,k)*duy(i,j,k)
           avg_s(i,j,144)= avg_s(i,j,144) + Frhov(i,j,k)*dvx(i,j,k)
           avg_s(i,j,145)= avg_s(i,j,145) + Frhov(i,j,k)*dvy(i,j,k)
           avg_s(i,j,146)= avg_s(i,j,146) + Frhow(i,j,k)*duz(i,j,k)
           avg_s(i,j,147)= avg_s(i,j,147) + Frhow(i,j,k)*dvz(i,j,k)
           avg_s(i,j,148)= avg_s(i,j,148) + Frhow(i,j,k)*dwx(i,j,k)
           avg_s(i,j,149)= avg_s(i,j,149) + Grhov(i,j,k)*duy(i,j,k)
           avg_s(i,j,150)= avg_s(i,j,150) + Grhov(i,j,k)*dvy(i,j,k)
           avg_s(i,j,151)= avg_s(i,j,151) + Grhow(i,j,k)*duz(i,j,k)
           avg_s(i,j,152)= avg_s(i,j,152) + Grhow(i,j,k)*dvz(i,j,k)
           avg_s(i,j,153)= avg_s(i,j,153) + Grhow(i,j,k)*dwy(i,j,k)
           avg_s(i,j,154)= avg_s(i,j,154) + Hrhow(i,j,k)*dwz(i,j,k)
           avg_s(i,j,155)= avg_s(i,j,155) + la*dTx(i,j,k)
           avg_s(i,j,156)= avg_s(i,j,156) + la*dTy(i,j,k)
           avg_s(i,j,157)= avg_s(i,j,157) + la*dTz(i,j,k)
           avg_s(i,j,158)= avg_s(i,j,158) + h_i*u_i
           avg_s(i,j,159)= avg_s(i,j,159) + h_i*v_i
           avg_s(i,j,160)= avg_s(i,j,160) + h_i*w_i
           avg_s(i,j,161)= avg_s(i,j,161) + rho_i*h_i*u_i
           avg_s(i,j,162)= avg_s(i,j,162) + rho_i*h_i*v_i
           avg_s(i,j,163)= avg_s(i,j,163) + rho_i*h_i*w_i
           avg_s(i,j,164)= avg_s(i,j,164) + rho_i*u_i**3
           avg_s(i,j,165)= avg_s(i,j,165) + rho_i*v_i**3
           avg_s(i,j,166)= avg_s(i,j,166) + rho_i*w_i**3
           avg_s(i,j,167)= avg_s(i,j,167) + rho_i*w_i**2*u_i
        enddo
     enddo
  enddo

  avg_s = avg_s/dble(nz)

  avg_t = ((nn1-1.0_wp)*avg_t + avg_s)*inn1

  ! if (ntotal.lt.ndeb2) return

  ! nn1 = dble(ntotal+1-ndeb2)
  ! nn1 = dble(ntotal+1-ndeb)
  ! inn1= 1.0_wp/nn1

  !u0  = ((nn1-1.0_wp)*u0 + uu) *inn1
  !v0  = ((nn1-1.0_wp)*v0 + vv) *inn1
  !w0  = ((nn1-1.0_wp)*w0 + ww) *inn1

end subroutine stats_stbl

!===============================================================================
subroutine stats_cyl
!===============================================================================
  !> Compute stats for cylinder flow
!===============================================================================
  use mod_mpi
  use mod_eos
  use mod_flow
  use mod_time
  use mod_constant
  use mod_rans
  use mod_deriv_c
  use mod_deriv_c3
  implicit none
  ! ----------------------------------------------------------------------------
  integer  :: i,j,k
  real(wp) :: nn1,inn1,cc2,qq,cnu13
  real(wp) :: g,s_i,divloc,cp,cv,p_i,T_i,rho_i
  real(wp) :: u_i,v_i,w_i,b1,b2,b3_,mu,la,prr,eck,e_i,h_i,nu!,c_fact
  real(wp), dimension(:,:,:), allocatable :: drhox,drhoy
  ! ----------------------------------------------------------------------------

  ! Collect statistics
  ! ------------------
  ! ! Original version
  ! if (ntotal.lt.ndeb) return
  ! nn1 = dble(ntotal+1-ndeb)/dble(freq_stats)
  ! inn1= 1.0_wp/nn1

  ! ! Camille version
  ! nn1 = dble(ntotal+1-ndeb)
  ! nn1 = nn1+dble(ntotal-ndeb)/dble(freq_stats)*dble(1.0_wp-freq_stats)
  ! inn1= 1.0_wp/nn1

  ! ! Incorrect correction proposed by AB
  ! if (ntotal.le.ndeb) return
  ! nn1 = dble(ntotal-ndeb)/dble(freq_stats)
  ! inn1= 1.0_wp/nn1

  ! New correct
  if (ntotal.lt.ndeb) return
  nn1 = dble(int((ntotal-ndeb)/freq_stats)) + 1
  inn1= 1.0_wp/nn1

  ! Spatial local and global averages at this iteration
  avg_s  = 0.0_wp

  cnu1=7.1_wp
  cnu13=cnu1**3

  ! Calculation of rho gradient
  allocate(drhox(1:nx,1:ny,1:nz))
  allocate(drhoy(1:nx,1:ny,1:nz))
  if (is_curv3) then
     call deriv_x_5pts_c3(rho(nx1:nx2,ny1:ny2,nz1:nz2),drhox(1:nx,1:ny,1:nz),1,nx,1,ny,1,nz)
     call deriv_y_5pts_c3(rho(nx1:nx2,ny1:ny2,nz1:nz2),drhoy(1:nx,1:ny,1:nz),1,nx,1,ny,1,nz)
  elseif (is_curv) then
     call deriv_x_5pts_c(rho(nx1:nx2,ny1:ny2,nz1:nz2),drhox(1:nx,1:ny,1:nz),1,nx,1,ny,1,nz)
     call deriv_y_5pts_c(rho(nx1:nx2,ny1:ny2,nz1:nz2),drhoy(1:nx,1:ny,1:nz),1,nx,1,ny,1,nz)
  else
     call mpistop('Stats cyl not implemented in cartesian',0)
  endif

  do k=1,nz
     do j=1,ny
        do i=1,nx
           rho_i= rho(i,j,k)

           if (is_RANS) then
              nu  = visc(i,j,k)/rho(i,j,k)
              khi = nutil(i,j,k)/nu
              fnu1 = khi**3/(khi**3+cnu13)
              mut(i,j,k) = nutil(i,j,k)*fnu1*rho_i
           endif

           u_i = uu(i,j,k)
           v_i = vv(i,j,k)
           w_i = ww(i,j,k)
           T_i = Tmp(i,j,k)
           p_i = prs(i,j,k)

           g = gcalc_tro(T_i,rho_i)
           !c_fact = p_i/T_i/rho_i
           cc2= c2calc_tro(T_i,rho_i)
           cv= cvcalc_tro(T_i,rho_i)
           cp= cpcalc_tro(T_i,rho_i)
           s_i = scalc_tro(T_i,rho_i)

           qq = u_i**2 + v_i**2 + w_i**2

           e_i = rhoe_n(i,j,k)/rho_i - 0.5_wp*qq
           h_i = e_i + p_i/rho_i

           b1 = dwy(i,j,k) - dvz(i,j,k)
           b2 = duz(i,j,k) - dwx(i,j,k)
           b3_= dvx(i,j,k) - duy(i,j,k)
           divloc = dux(i,j,k) + dvy(i,j,k) + dwz(i,j,k)

           mu = visc(i,j,k)
           la = cok(i,j,k)

           prr = cp*mu/la
           eck= qq/(cp*T_i)

           avg_s(i,j,  1)= avg_s(i,j,  1) + rho_i
           avg_s(i,j,  2)= avg_s(i,j,  2) + u_i
           avg_s(i,j,  3)= avg_s(i,j,  3) + v_i
           avg_s(i,j,  4)= avg_s(i,j,  4) + w_i
           avg_s(i,j,  5)= avg_s(i,j,  5) + p_i
           avg_s(i,j,  6)= avg_s(i,j,  6) + T_i
           avg_s(i,j,  7)= avg_s(i,j,  7) + rhou(i,j,k)
           avg_s(i,j,  8)= avg_s(i,j,  8) + rhov(i,j,k)
           avg_s(i,j,  9)= avg_s(i,j,  9) + rhow(i,j,k)
           avg_s(i,j, 10)= avg_s(i,j, 10) + rhoe(i,j,k)
           avg_s(i,j, 11)= avg_s(i,j, 11) + rho_i**2
           avg_s(i,j, 12)= avg_s(i,j, 12) + u_i**2
           avg_s(i,j, 13)= avg_s(i,j, 13) + v_i**2
           avg_s(i,j, 14)= avg_s(i,j, 14) + w_i**2
           avg_s(i,j, 15)= avg_s(i,j, 15) + u_i*v_i
           avg_s(i,j, 16)= avg_s(i,j, 16) + u_i*w_i
           avg_s(i,j, 17)= avg_s(i,j, 17) + v_i*w_i
           avg_s(i,j, 18)= avg_s(i,j, 18) + v_i*T_i
           avg_s(i,j, 19)= avg_s(i,j, 19) + p_i**2
           avg_s(i,j, 20)= avg_s(i,j, 20) + T_i**2
           avg_s(i,j, 21)= avg_s(i,j, 21) + mu
           avg_s(i,j, 22)= avg_s(i,j, 22) + divloc
           if (is_RANS) then
              avg_s(i,j, 23)= avg_s(i,j, 23) + mut(i,j,k)
           else
              avg_s(i,j, 23)= avg_s(i,j, 23) + divloc**2
           endif

           avg_s(i,j, 24)= avg_s(i,j, 24) + e_i
           avg_s(i,j, 25)= avg_s(i,j, 25) + h_i
           avg_s(i,j, 26)= avg_s(i,j, 26) + sqrt(cc2)
           avg_s(i,j, 27)= avg_s(i,j, 27) + s_i
           avg_s(i,j, 28)= avg_s(i,j, 28) + sqrt(qq/cc2)
           avg_s(i,j, 29)= avg_s(i,j, 29) + 0.5*qq
           avg_s(i,j, 30)= avg_s(i,j, 30) + g
           avg_s(i,j, 31)= avg_s(i,j, 31) + la
           avg_s(i,j, 32)= avg_s(i,j, 32) + cp
           avg_s(i,j, 33)= avg_s(i,j, 33) + cv
           avg_s(i,j, 34)= avg_s(i,j, 34) + prr
           avg_s(i,j, 35)= avg_s(i,j, 35) + eck
           avg_s(i,j, 36)= avg_s(i,j, 36) + rho_i*dux(i,j,k)
           avg_s(i,j, 37)= avg_s(i,j, 37) + rho_i*duy(i,j,k)
           avg_s(i,j, 38)= avg_s(i,j, 38) + rho_i*duz(i,j,k)
           avg_s(i,j, 39)= avg_s(i,j, 39) + rho_i*dvx(i,j,k)
           avg_s(i,j, 40)= avg_s(i,j, 40) + rho_i*dvy(i,j,k)
           avg_s(i,j, 41)= avg_s(i,j, 41) + rho_i*dvz(i,j,k)
           avg_s(i,j, 42)= avg_s(i,j, 42) + rho_i*dwx(i,j,k)
           avg_s(i,j, 43)= avg_s(i,j, 43) + rho_i*dwy(i,j,k)
           avg_s(i,j, 44)= avg_s(i,j, 44) + rho_i*dwz(i,j,k)
           avg_s(i,j, 45)= avg_s(i,j, 45) + p_i*divloc
           avg_s(i,j, 46)= avg_s(i,j, 46) + rho_i*divloc
           avg_s(i,j, 47)= avg_s(i,j, 47) + b1
           avg_s(i,j, 48)= avg_s(i,j, 48) + b2
           avg_s(i,j, 49)= avg_s(i,j, 49) + b3_
           avg_s(i,j, 50)= avg_s(i,j, 50) + rho_i*T_i
           avg_s(i,j, 51)= avg_s(i,j, 51) + u_i*T_i
           avg_s(i,j, 52)= avg_s(i,j, 52) + v_i*T_i
           avg_s(i,j, 53)= avg_s(i,j, 53) + e_i**2
           avg_s(i,j, 54)= avg_s(i,j, 54) + h_i**2
           avg_s(i,j, 55)= avg_s(i,j, 55) + cc2
           avg_s(i,j, 56)= avg_s(i,j, 56) + s_i**2
           avg_s(i,j, 57)= avg_s(i,j, 57) + qq/cc2
           avg_s(i,j, 58)= avg_s(i,j, 58) + g**2
           avg_s(i,j, 59)= avg_s(i,j, 59) + mu**2
           avg_s(i,j, 60)= avg_s(i,j, 60) + la**2
           avg_s(i,j, 61)= avg_s(i,j, 61) + cv**2
           avg_s(i,j, 62)= avg_s(i,j, 62) + cp**2
           avg_s(i,j, 63)= avg_s(i,j, 63) + prr**2
           avg_s(i,j, 64)= avg_s(i,j, 64) + eck**2
           avg_s(i,j, 65)= avg_s(i,j, 65) + p_i*u_i
           avg_s(i,j, 66)= avg_s(i,j, 66) + p_i*v_i
           avg_s(i,j, 67)= avg_s(i,j, 67) + s_i*u_i
           avg_s(i,j, 68)= avg_s(i,j, 68) + s_i*v_i
           avg_s(i,j, 69)= avg_s(i,j, 69) + p_i*rho_i
           avg_s(i,j, 70)= avg_s(i,j, 70) + h_i*rho_i
           avg_s(i,j, 71)= avg_s(i,j, 71) + T_i*p_i
           avg_s(i,j, 72)= avg_s(i,j, 72) + p_i*s_i
           avg_s(i,j, 73)= avg_s(i,j, 73) + T_i*s_i
           avg_s(i,j, 74)= avg_s(i,j, 74) + rho_i*s_i
           avg_s(i,j, 75)= avg_s(i,j, 75) + g*rho_i
           avg_s(i,j, 76)= avg_s(i,j, 76) + g*p_i
           ! avg_s(i,j, 77)= avg_s(i,j, 77) + g*s_i
           ! avg_s(i,j, 78)= avg_s(i,j, 78) + g*T_i
           ! avg_s(i,j, 79)= avg_s(i,j, 79) + g*u_i
           ! avg_s(i,j, 80)= avg_s(i,j, 80) + g*v_i
           avg_s(i,j, 77)= avg_s(i,j, 77) + drhox(i,j,k)
           avg_s(i,j, 78)= avg_s(i,j, 78) + drhoy(i,j,k)
           avg_s(i,j, 79)= avg_s(i,j, 79) + drhox(i,j,k)**2
           avg_s(i,j, 80)= avg_s(i,j, 80) + drhoy(i,j,k)**2
           avg_s(i,j, 81)= avg_s(i,j, 81) + p_i*dux(i,j,k)
           avg_s(i,j, 82)= avg_s(i,j, 82) + p_i*dvy(i,j,k)
           avg_s(i,j, 83)= avg_s(i,j, 83) + p_i*dwz(i,j,k)
           avg_s(i,j, 84)= avg_s(i,j, 84) + p_i*duy(i,j,k)
           avg_s(i,j, 85)= avg_s(i,j, 85) + p_i*dvx(i,j,k)
           avg_s(i,j, 86)= avg_s(i,j, 86) + rho_i*divloc**2
           avg_s(i,j, 87)= avg_s(i,j, 87) + dux(i,j,k)**2
           avg_s(i,j, 88)= avg_s(i,j, 88) + duy(i,j,k)**2
           avg_s(i,j, 89)= avg_s(i,j, 89) + duz(i,j,k)**2
           avg_s(i,j, 90)= avg_s(i,j, 90) + dvx(i,j,k)**2
           avg_s(i,j, 91)= avg_s(i,j, 91) + dvy(i,j,k)**2
           avg_s(i,j, 92)= avg_s(i,j, 92) + dvz(i,j,k)**2
           avg_s(i,j, 93)= avg_s(i,j, 93) + dwx(i,j,k)**2
           avg_s(i,j, 94)= avg_s(i,j, 94) + dwy(i,j,k)**2
           avg_s(i,j, 95)= avg_s(i,j, 95) + dwz(i,j,k)**2
           avg_s(i,j, 96)= avg_s(i,j, 96) + b1**2
           avg_s(i,j, 97)= avg_s(i,j, 97) + b2**2
           avg_s(i,j, 98)= avg_s(i,j, 98) + b3_**2
           avg_s(i,j, 99)= avg_s(i,j, 99) + rho_i*b1
           avg_s(i,j,100)= avg_s(i,j,100) + rho_i*b2
           avg_s(i,j,101)= avg_s(i,j,101) + rho_i*b3_
           avg_s(i,j,102)= avg_s(i,j,102) + rho_i*u_i**2
           avg_s(i,j,103)= avg_s(i,j,103) + rho_i*v_i**2
           avg_s(i,j,104)= avg_s(i,j,104) + rho_i*w_i**2
           avg_s(i,j,105)= avg_s(i,j,105) + rho_i*T_i**2
           avg_s(i,j,106)= avg_s(i,j,106) + rho_i*b1**2
           avg_s(i,j,107)= avg_s(i,j,107) + rho_i*b2**2
           avg_s(i,j,108)= avg_s(i,j,108) + rho_i*b3_**2
           avg_s(i,j,109)= avg_s(i,j,109) + rho_i*u_i*v_i
           avg_s(i,j,110)= avg_s(i,j,110) + rho_i*u_i*w_i
           avg_s(i,j,111)= avg_s(i,j,111) + rho_i*v_i*w_i
           avg_s(i,j,112)= avg_s(i,j,112) + rho_i*v_i*T_i
           avg_s(i,j,113)= avg_s(i,j,113) + rho_i*u_i**2*v_i
           avg_s(i,j,114)= avg_s(i,j,114) + rho_i*v_i**2*v_i
           avg_s(i,j,115)= avg_s(i,j,115) + rho_i*w_i**2*v_i
           avg_s(i,j,116)= avg_s(i,j,116) + rho_i*v_i**2*u_i
           avg_s(i,j,117)= avg_s(i,j,117) + rho_i*dux(i,j,k)**2
           avg_s(i,j,118)= avg_s(i,j,118) + rho_i*dvy(i,j,k)**2
           avg_s(i,j,119)= avg_s(i,j,119) + rho_i*dwz(i,j,k)**2
           avg_s(i,j,120)= avg_s(i,j,120) + rho_i*duy(i,j,k)*dvx(i,j,k)
           avg_s(i,j,121)= avg_s(i,j,121) + rho_i*duz(i,j,k)*dwx(i,j,k)
           avg_s(i,j,122)= avg_s(i,j,122) + rho_i*dvz(i,j,k)*dwy(i,j,k)
           avg_s(i,j,123)= avg_s(i,j,123) + u_i**3
           avg_s(i,j,124)= avg_s(i,j,124) + p_i**3
           avg_s(i,j,125)= avg_s(i,j,125) + u_i**4
           avg_s(i,j,126)= avg_s(i,j,126) + p_i**4
           avg_s(i,j,127)= avg_s(i,j,127) + Frhou(i,j,k)
           avg_s(i,j,128)= avg_s(i,j,128) + Frhov(i,j,k)
           avg_s(i,j,129)= avg_s(i,j,129) + Frhow(i,j,k)
           avg_s(i,j,130)= avg_s(i,j,130) + Grhov(i,j,k)
           avg_s(i,j,131)= avg_s(i,j,131) + Grhow(i,j,k)
           avg_s(i,j,132)= avg_s(i,j,132) + Hrhow(i,j,k)
           avg_s(i,j,133)= avg_s(i,j,133) + Frhov(i,j,k)*u_i
           avg_s(i,j,134)= avg_s(i,j,134) + Frhou(i,j,k)*u_i
           avg_s(i,j,135)= avg_s(i,j,135) + Frhov(i,j,k)*v_i
           avg_s(i,j,136)= avg_s(i,j,136) + Frhow(i,j,k)*w_i
           avg_s(i,j,137)= avg_s(i,j,137) + Grhov(i,j,k)*u_i
           avg_s(i,j,138)= avg_s(i,j,138) + Grhov(i,j,k)*v_i
           avg_s(i,j,139)= avg_s(i,j,139) + Grhow(i,j,k)*w_i
           avg_s(i,j,140)= avg_s(i,j,140) + Frhou(i,j,k)*dux(i,j,k)
           avg_s(i,j,141)= avg_s(i,j,141) + Frhou(i,j,k)*dvx(i,j,k)
           avg_s(i,j,142)= avg_s(i,j,142) + Frhov(i,j,k)*dux(i,j,k)
           avg_s(i,j,143)= avg_s(i,j,143) + Frhov(i,j,k)*duy(i,j,k)
           avg_s(i,j,144)= avg_s(i,j,144) + Frhov(i,j,k)*dvx(i,j,k)
           avg_s(i,j,145)= avg_s(i,j,145) + Frhov(i,j,k)*dvy(i,j,k)
           avg_s(i,j,146)= avg_s(i,j,146) + Frhow(i,j,k)*duz(i,j,k)
           avg_s(i,j,147)= avg_s(i,j,147) + Frhow(i,j,k)*dvz(i,j,k)
           avg_s(i,j,148)= avg_s(i,j,148) + Frhow(i,j,k)*dwx(i,j,k)
           avg_s(i,j,149)= avg_s(i,j,149) + Grhov(i,j,k)*duy(i,j,k)
           avg_s(i,j,150)= avg_s(i,j,150) + Grhov(i,j,k)*dvy(i,j,k)
           avg_s(i,j,151)= avg_s(i,j,151) + Grhow(i,j,k)*duz(i,j,k)
           avg_s(i,j,152)= avg_s(i,j,152) + Grhow(i,j,k)*dvz(i,j,k)
           avg_s(i,j,153)= avg_s(i,j,153) + Grhow(i,j,k)*dwy(i,j,k)
           avg_s(i,j,154)= avg_s(i,j,154) + Hrhow(i,j,k)*dwz(i,j,k)
           avg_s(i,j,155)= avg_s(i,j,155) + la*dTx(i,j,k)
           avg_s(i,j,156)= avg_s(i,j,156) + la*dTy(i,j,k)
           avg_s(i,j,157)= avg_s(i,j,157) + la*dTz(i,j,k)
           avg_s(i,j,158)= avg_s(i,j,158) + h_i*u_i
           avg_s(i,j,159)= avg_s(i,j,159) + h_i*v_i
           avg_s(i,j,160)= avg_s(i,j,160) + h_i*w_i
           avg_s(i,j,161)= avg_s(i,j,161) + rho_i*h_i*u_i
           avg_s(i,j,162)= avg_s(i,j,162) + rho_i*h_i*v_i
           avg_s(i,j,163)= avg_s(i,j,163) + rho_i*h_i*w_i
           avg_s(i,j,164)= avg_s(i,j,164) + rho_i*u_i**3
           avg_s(i,j,165)= avg_s(i,j,165) + rho_i*v_i**3
           avg_s(i,j,166)= avg_s(i,j,166) + rho_i*w_i**3
           avg_s(i,j,167)= avg_s(i,j,167) + rho_i*w_i**2*u_i
           !avg_s(i,j,168)= avg_s(i,j,168) + dux(i,j,k)
           !avg_s(i,j,169)= avg_s(i,j,169) + duy(i,j,k)
           !avg_s(i,j,170)= avg_s(i,j,170) + duz(i,j,k)
           !avg_s(i,j,171)= avg_s(i,j,171) + dvx(i,j,k)
           !avg_s(i,j,172)= avg_s(i,j,172) + dvy(i,j,k)
           !avg_s(i,j,173)= avg_s(i,j,173) + dvz(i,j,k)
           !avg_s(i,j,174)= avg_s(i,j,174) + dwx(i,j,k)
           !avg_s(i,j,175)= avg_s(i,j,175) + dwy(i,j,k)
           !avg_s(i,j,176)= avg_s(i,j,176) + dwz(i,j,k)
        enddo
     enddo
  enddo

  avg_s = avg_s/dble(nz)

  avg_t = ((nn1-1.0_wp)*avg_t + avg_s)*inn1

end subroutine stats_cyl

!===============================================================================
subroutine stats_shit
!===============================================================================
  !> Compute stats for boundary layers
!===============================================================================
  use mod_mpi
  use mod_eos
  use mod_flow
  use mod_time
  use mod_constant
  !use mod_artvisc
  implicit none
  ! ----------------------------------------------------------------------------
  integer  :: i,j,k
  real(wp) :: nn1,inn1,cc2,qq
  real(wp) :: g,s_i,divloc,cp,cv,p_i,T_i,rho_i
  real(wp) :: u_i,v_i,w_i,b1,b2,b3,mu,la,prr,eck,e_i,h_i
  ! ----------------------------------------------------------------------------

  ! Collect statistics
  ! ------------------
  if (ntotal.le.ndeb) return
  if ((mod(ntime,nprint).eq.0).and.(iproc.eq.0)) print *,"~> Averaging for Spatial HIT"
  nn1 = dble(ntotal-ndeb)/dble(freq_stats)
  inn1= 1.0_wp/nn1

  ! Spatial local and global averages at this iteration
  avg_s  = 0.0_wp

  do k=1,nz
     do j=1,ny
        do i=1,nx
           rho_i= rho(i,j,k)
           u_i = uu(i,j,k)-U_ref
           v_i = vv(i,j,k)
           w_i = ww(i,j,k)
           T_i = Tmp(i,j,k)
           p_i = prs(i,j,k)

           g = gcalc_tro(T_i,rho_i)
           cc2= c2calc_tro(T_i,rho_i)
           cv= cvcalc_tro(T_i,rho_i)
           cp= cpcalc_tro(T_i,rho_i)
           s_i = scalc_tro(T_i,rho_i)

           qq = u_i**2 + v_i**2 + w_i**2

           e_i = rhoe_n(i,j,k)/rho_i - 0.5_wp*qq
           h_i = e_i + p_i/rho_i

           b1 = dwy(i,j,k) - dvz(i,j,k)
           b2 = duz(i,j,k) - dwx(i,j,k)
           b3 = dvx(i,j,k) - duy(i,j,k)
           divloc = dux(i,j,k) + dvy(i,j,k) + dwz(i,j,k)

           mu = visc(i,j,k)
           la = cok(i,j,k)

           prr = cp*mu/la
           eck= qq/(cp*T_i)

           avg_s(i,j,  1)= avg_s(i,j,  1) + rho_i
           avg_s(i,j,  2)= avg_s(i,j,  2) + u_i
           avg_s(i,j,  3)= avg_s(i,j,  3) + v_i
           avg_s(i,j,  4)= avg_s(i,j,  4) + w_i
           avg_s(i,j,  5)= avg_s(i,j,  5) + p_i
           avg_s(i,j,  6)= avg_s(i,j,  6) + T_i
           avg_s(i,j,  7)= avg_s(i,j,  7) + rhou(i,j,k)
           avg_s(i,j,  8)= avg_s(i,j,  8) + rhov(i,j,k)
           avg_s(i,j,  9)= avg_s(i,j,  9) + rhow(i,j,k)
           avg_s(i,j, 10)= avg_s(i,j, 10) + rhoe(i,j,k)
           avg_s(i,j, 11)= avg_s(i,j, 11) + rho_i**2
           avg_s(i,j, 12)= avg_s(i,j, 12) + u_i**2
           avg_s(i,j, 13)= avg_s(i,j, 13) + v_i**2
           avg_s(i,j, 14)= avg_s(i,j, 14) + w_i**2
           avg_s(i,j, 15)= avg_s(i,j, 15) + u_i*v_i
           avg_s(i,j, 16)= avg_s(i,j, 16) + u_i*w_i
           avg_s(i,j, 17)= avg_s(i,j, 17) + v_i*w_i
           avg_s(i,j, 18)= avg_s(i,j, 18) + v_i*T_i
           avg_s(i,j, 19)= avg_s(i,j, 19) + p_i**2
           avg_s(i,j, 20)= avg_s(i,j, 20) + T_i**2
           avg_s(i,j, 21)= avg_s(i,j, 21) + mu
           avg_s(i,j, 22)= avg_s(i,j, 22) + divloc
           avg_s(i,j, 23)= avg_s(i,j, 23) + divloc**2

           avg_s(i,j, 24)= avg_s(i,j, 24) + e_i
           avg_s(i,j, 25)= avg_s(i,j, 25) + h_i
           avg_s(i,j, 26)= avg_s(i,j, 26) + sqrt(cc2)
           avg_s(i,j, 27)= avg_s(i,j, 27) + s_i
           avg_s(i,j, 28)= avg_s(i,j, 28) + sqrt(qq/cc2)
           avg_s(i,j, 29)= avg_s(i,j, 29) + 0.5*qq
           avg_s(i,j, 30)= avg_s(i,j, 30) + sqrt(qq)/U_ref
        enddo
     enddo
  enddo

  avg_s = avg_s/dble(nz)

  avg_t = ((nn1-1.0_wp)*avg_t + avg_s)*inn1

  ! if (ntotal.lt.ndeb2) return

  ! nn1 = dble(ntotal+1-ndeb2)
  ! nn1 = dble(ntotal+1-ndeb)
  ! inn1= 1.0_wp/nn1

  !u0  = ((nn1-1.0_wp)*u0 + uu) *inn1
  !v0  = ((nn1-1.0_wp)*v0 + vv) *inn1
  !w0  = ((nn1-1.0_wp)*w0 + ww) *inn1

end subroutine stats_shit

!===============================================================================
subroutine stats_act
!===============================================================================
  !> Compute stats for the actuator problem
!===============================================================================
  use mod_mpi
  use mod_eos
  use mod_flow
  use mod_time
  use mod_constant
  !use mod_artvisc
  implicit none
  ! ----------------------------------------------------------------------------
  integer  :: i,j,k
  real(wp) :: nn1,inn1,cc2,qq
  real(wp) :: s_i,divloc,p_i,T_i,rho_i
  real(wp) :: u_i,v_i,w_i,b1,b2,b3,mu,e_i,h_i
  ! ----------------------------------------------------------------------------

  ! Collect statistics
  ! ------------------
  if (ntotal.lt.ndeb) return
  nn1 = dble(int((ntotal-ndeb)/freq_stats)) + 1
  inn1= 1.0_wp/nn1

  ! Spatial local and global averages at this iteration
  avg_s  = 0.0_wp

  do k=1,nz
     do j=1,ny
        do i=1,nx
           rho_i= rho(i,j,k)
           u_i = uu(i,j,k)
           v_i = vv(i,j,k)
           w_i = ww(i,j,k)
           T_i = Tmp(i,j,k)
           p_i = prs(i,j,k)

!            g = gcalc_tro(T_i,rho_i)
           cc2= c2calc_tro(T_i,rho_i)
!            cv= cvcalc_tro(T_i,rho_i)
!            cp= cpcalc_tro(T_i,rho_i)
           s_i = scalc_tro(T_i,rho_i)

           qq = u_i**2 + v_i**2 + w_i**2

           e_i = rhoe_n(i,j,k)/rho_i - 0.5_wp*qq
           h_i = e_i + p_i/rho_i

           b1 = dwy(i,j,k) - dvz(i,j,k)
           b2 = duz(i,j,k) - dwx(i,j,k)
           b3 = dvx(i,j,k) - duy(i,j,k)
           divloc = dux(i,j,k) + dvy(i,j,k) + dwz(i,j,k)

           mu = visc(i,j,k)
!            la = cok(i,j,k)

!            prr = cp*mu/la
!            eck= qq/(cp*T_i)

           avg_s(i,j,  1)= avg_s(i,j,  1) + rho_i
           avg_s(i,j,  2)= avg_s(i,j,  2) + u_i
           avg_s(i,j,  3)= avg_s(i,j,  3) + v_i
           avg_s(i,j,  4)= avg_s(i,j,  4) + w_i
           avg_s(i,j,  5)= avg_s(i,j,  5) + p_i
           avg_s(i,j,  6)= avg_s(i,j,  6) + T_i
           avg_s(i,j,  7)= avg_s(i,j,  7) + rhou(i,j,k)
           avg_s(i,j,  8)= avg_s(i,j,  8) + rhov(i,j,k)
           avg_s(i,j,  9)= avg_s(i,j,  9) + rhow(i,j,k)
           avg_s(i,j, 10)= avg_s(i,j, 10) + rhoe(i,j,k)
           avg_s(i,j, 11)= avg_s(i,j, 11) + rho_i**2
           avg_s(i,j, 12)= avg_s(i,j, 12) + u_i**2
           avg_s(i,j, 13)= avg_s(i,j, 13) + v_i**2
           avg_s(i,j, 14)= avg_s(i,j, 14) + w_i**2
           avg_s(i,j, 15)= avg_s(i,j, 15) + u_i*v_i
           avg_s(i,j, 16)= avg_s(i,j, 16) + u_i*w_i
           avg_s(i,j, 17)= avg_s(i,j, 17) + v_i*w_i
           avg_s(i,j, 18)= avg_s(i,j, 18) + v_i*T_i
           avg_s(i,j, 19)= avg_s(i,j, 19) + p_i**2
           avg_s(i,j, 20)= avg_s(i,j, 20) + T_i**2
           avg_s(i,j, 21)= avg_s(i,j, 21) + mu
           avg_s(i,j, 22)= avg_s(i,j, 22) + divloc
           avg_s(i,j, 23)= avg_s(i,j, 23) + divloc**2
           avg_s(i,j, 24)= avg_s(i,j, 24) + e_i
           avg_s(i,j, 25)= avg_s(i,j, 25) + h_i
           avg_s(i,j, 26)= avg_s(i,j, 26) + sqrt(cc2)
           avg_s(i,j, 27)= avg_s(i,j, 27) + s_i
           avg_s(i,j, 28)= avg_s(i,j, 28) + sqrt(qq/cc2)
           avg_s(i,j, 29)= avg_s(i,j, 29) + 0.5*qq
           avg_s(i,j, 30)= avg_s(i,j, 30) + rho_i*dux(i,j,k)
           avg_s(i,j, 31)= avg_s(i,j, 31) + rho_i*duy(i,j,k)
           avg_s(i,j, 32)= avg_s(i,j, 32) + rho_i*duz(i,j,k)
           avg_s(i,j, 33)= avg_s(i,j, 33) + rho_i*dvx(i,j,k)
           avg_s(i,j, 34)= avg_s(i,j, 34) + rho_i*dvy(i,j,k)
           avg_s(i,j, 35)= avg_s(i,j, 35) + rho_i*dvz(i,j,k)
           avg_s(i,j, 36)= avg_s(i,j, 36) + rho_i*dwx(i,j,k)
           avg_s(i,j, 37)= avg_s(i,j, 37) + rho_i*dwy(i,j,k)
           avg_s(i,j, 38)= avg_s(i,j, 38) + rho_i*dwz(i,j,k)
           avg_s(i,j, 39)= avg_s(i,j, 39) + p_i*divloc
           avg_s(i,j, 40)= avg_s(i,j, 40) + rho_i*divloc
           avg_s(i,j, 41)= avg_s(i,j, 41) + b1
           avg_s(i,j, 42)= avg_s(i,j, 42) + b2
           avg_s(i,j, 43)= avg_s(i,j, 43) + b3
           avg_s(i,j, 44)= avg_s(i,j, 44) + rho_i*T_i
           avg_s(i,j, 45)= avg_s(i,j, 45) + u_i*T_i
           avg_s(i,j, 46)= avg_s(i,j, 46) + v_i*T_i
           avg_s(i,j, 47)= avg_s(i,j, 47) + e_i**2
           avg_s(i,j, 48)= avg_s(i,j, 48) + h_i**2
           avg_s(i,j, 49)= avg_s(i,j, 49) + cc2
           avg_s(i,j, 50)= avg_s(i,j, 50) + s_i**2
           avg_s(i,j, 51)= avg_s(i,j, 51) + qq/cc2
           avg_s(i,j, 52)= avg_s(i,j, 52) + mu**2
           avg_s(i,j, 53)= avg_s(i,j, 53) + p_i*u_i
           avg_s(i,j, 54)= avg_s(i,j, 54) + p_i*v_i
           avg_s(i,j, 55)= avg_s(i,j, 55) + p_i*rho_i
           avg_s(i,j, 56)= avg_s(i,j, 56) + T_i*p_i
           avg_s(i,j, 57)= avg_s(i,j, 57) + p_i*dux(i,j,k)
           avg_s(i,j, 58)= avg_s(i,j, 58) + p_i*dvy(i,j,k)
           avg_s(i,j, 59)= avg_s(i,j, 59) + p_i*dwz(i,j,k)
           avg_s(i,j, 60)= avg_s(i,j, 60) + p_i*duy(i,j,k)
           avg_s(i,j, 61)= avg_s(i,j, 61) + p_i*dvx(i,j,k)
           avg_s(i,j, 62)= avg_s(i,j, 62) + rho_i*divloc**2
           avg_s(i,j, 63)= avg_s(i,j, 63) + dux(i,j,k)**2
           avg_s(i,j, 64)= avg_s(i,j, 64) + duy(i,j,k)**2
           avg_s(i,j, 65)= avg_s(i,j, 65) + duz(i,j,k)**2
           avg_s(i,j, 66)= avg_s(i,j, 66) + dvx(i,j,k)**2
           avg_s(i,j, 67)= avg_s(i,j, 67) + dvy(i,j,k)**2
           avg_s(i,j, 68)= avg_s(i,j, 68) + dvz(i,j,k)**2
           avg_s(i,j, 69)= avg_s(i,j, 69) + dwx(i,j,k)**2
           avg_s(i,j, 70)= avg_s(i,j, 70) + dwy(i,j,k)**2
           avg_s(i,j, 71)= avg_s(i,j, 71) + dwz(i,j,k)**2
           avg_s(i,j, 72)= avg_s(i,j, 72) + b1**2
           avg_s(i,j, 73)= avg_s(i,j, 73) + b2**2
           avg_s(i,j, 74)= avg_s(i,j, 74) + b3**2
           avg_s(i,j, 75)= avg_s(i,j, 75) + rho_i*b1
           avg_s(i,j, 76)= avg_s(i,j, 76) + rho_i*b2
           avg_s(i,j, 77)= avg_s(i,j, 77) + rho_i*b3
           avg_s(i,j, 78)= avg_s(i,j, 78) + rho_i*u_i**2
           avg_s(i,j, 79)= avg_s(i,j, 79) + rho_i*v_i**2
           avg_s(i,j, 80)= avg_s(i,j, 80) + rho_i*w_i**2
           avg_s(i,j, 81)= avg_s(i,j, 81) + rho_i*u_i*v_i
           avg_s(i,j, 82)= avg_s(i,j, 82) + rho_i*u_i*w_i
           avg_s(i,j, 83)= avg_s(i,j, 83) + rho_i*v_i*w_i
           avg_s(i,j, 84)= avg_s(i,j, 84) + p_i**3
           avg_s(i,j, 85)= avg_s(i,j, 85) + p_i**4
           avg_s(i,j, 86)= avg_s(i,j, 86) + Frhou(i,j,k)
           avg_s(i,j, 87)= avg_s(i,j, 87) + Frhov(i,j,k)
           avg_s(i,j, 88)= avg_s(i,j, 88) + Frhow(i,j,k)
           avg_s(i,j, 89)= avg_s(i,j, 89) + Grhov(i,j,k)
           avg_s(i,j, 90)= avg_s(i,j, 90) + Grhow(i,j,k)
           avg_s(i,j, 91)= avg_s(i,j, 91) + Hrhow(i,j,k)
           avg_s(i,j, 92)= avg_s(i,j, 92) + Frhov(i,j,k)*u_i
           avg_s(i,j, 93)= avg_s(i,j, 93) + Frhou(i,j,k)*u_i
           avg_s(i,j, 94)= avg_s(i,j, 94) + Frhov(i,j,k)*v_i
           avg_s(i,j, 95)= avg_s(i,j, 95) + Frhow(i,j,k)*w_i
           avg_s(i,j, 96)= avg_s(i,j, 96) + Grhov(i,j,k)*u_i
           avg_s(i,j, 97)= avg_s(i,j, 97) + Grhov(i,j,k)*v_i
           avg_s(i,j, 98)= avg_s(i,j, 98) + Grhow(i,j,k)*w_i
        enddo
     enddo
  enddo

  avg_s = avg_s/dble(nz)

  avg_t = ((nn1-1.0_wp)*avg_t + avg_s)*inn1

end subroutine stats_act

!===============================================================================
subroutine stats_sphere
!===============================================================================
  !> Compute stats for flow past a sphere
!===============================================================================
  use mod_mpi
  use mod_flow
  use mod_time
  use mod_constant
  implicit none
  ! ----------------------------------------------------------------------------
  integer  :: i,j,k
  real(wp) :: nn1,nm1,inn1
  ! real(wp) :: cc2,qq,g,s_i,divloc,cp,cv,mu,la,prr,eck,e_i,h_i,nu!,c_fact
  real(wp) :: p_i,T_i,rho_i
  real(wp) :: u_i,v_i,w_i,b1,b2,b3
  ! ----------------------------------------------------------------------------

  ! Collect statistics
  ! ------------------
  if (ntotal.lt.ndeb) return
  nn1 = dble(int((ntotal-ndeb)/freq_stats)) + 1
  inn1= 1.0_wp/nn1
  nm1=nn1-1.0_wp

  do k=1,nz
     do j=1,ny
        do i=1,nx
           rho_i= rho(i,j,k)
           u_i = uu(i,j,k)
           v_i = vv(i,j,k)
           w_i = ww(i,j,k)
           T_i = Tmp(i,j,k)
           p_i = prs(i,j,k)

           b1 = dwy(i,j,k) - dvz(i,j,k)
           b2 = duz(i,j,k) - dwx(i,j,k)
           b3 = dvx(i,j,k) - duy(i,j,k)

           avg_v(i,j,k, 1)= (nm1*avg_v(i,j,k, 1) + rho_i)*inn1
           avg_v(i,j,k, 2)= (nm1*avg_v(i,j,k, 2) + u_i)*inn1
           avg_v(i,j,k, 3)= (nm1*avg_v(i,j,k, 3) + v_i)*inn1
           avg_v(i,j,k, 4)= (nm1*avg_v(i,j,k, 4) + w_i)*inn1
           avg_v(i,j,k, 5)= (nm1*avg_v(i,j,k, 5) + p_i)*inn1
           avg_v(i,j,k, 6)= (nm1*avg_v(i,j,k, 6) + T_i)*inn1
           avg_v(i,j,k, 7)= (nm1*avg_v(i,j,k, 7) + rho_i**2)*inn1
           avg_v(i,j,k, 8)= (nm1*avg_v(i,j,k, 8) + u_i**2)*inn1
           avg_v(i,j,k, 9)= (nm1*avg_v(i,j,k, 9) + v_i**2)*inn1
           avg_v(i,j,k,10)= (nm1*avg_v(i,j,k,10) + w_i**2)*inn1
           avg_v(i,j,k,11)= (nm1*avg_v(i,j,k,11) + u_i*v_i)*inn1
           avg_v(i,j,k,12)= (nm1*avg_v(i,j,k,12) + u_i*w_i)*inn1
           avg_v(i,j,k,13)= (nm1*avg_v(i,j,k,13) + v_i*w_i)*inn1
           avg_v(i,j,k,14)= (nm1*avg_v(i,j,k,14) + p_i**2)*inn1
           avg_v(i,j,k,15)= (nm1*avg_v(i,j,k,15) + T_i**2)*inn1
           avg_v(i,j,k,16)= (nm1*avg_v(i,j,k,16) + b1)*inn1
           avg_v(i,j,k,17)= (nm1*avg_v(i,j,k,17) + b2)*inn1
           avg_v(i,j,k,18)= (nm1*avg_v(i,j,k,18) + b3)*inn1
        enddo
     enddo
  enddo

end subroutine stats_sphere

!===============================================================================
subroutine stats_turb
!===============================================================================
  !> Compute stats for turbine flow
!===============================================================================
  use mod_mpi
  use mod_eos
  use mod_flow
  use mod_time
  use mod_constant
  use mod_rans
  implicit none
  ! ----------------------------------------------------------------------------
  integer  :: i,j,k
  real(wp) :: nn1,inn1,cc2,qq,cnu13
  real(wp) :: g,s_i,divloc,cp,cv,p_i,T_i,rho_i
  real(wp) :: u_i,v_i,w_i,b1,b2,b3,mu,la,prr,eck,e_i,h_i,nu!,c_fact
  real(wp) :: h0,e0,s0,T0,p0,ro0,M,g_eq
  ! ----------------------------------------------------------------------------

  ! Collect statistics
  ! ------------------
  if (ntotal.lt.ndeb) return
  nn1 = dble(int((ntotal-ndeb)/freq_stats)) + 1
  inn1= 1.0_wp/nn1

  ! Spatial local and global averages at this iteration
  avg_s  = 0.0_wp

  cnu1=7.1_wp
  cnu13=cnu1**3

  do k=1,nz
     do j=1,ny
        do i=1,nx
           ! if ((i.ne.20).and.(j.ne.1)) return

           rho_i = rho(i,j,k)

           u_i = uu(i,j,k)
           v_i = vv(i,j,k)
           w_i = ww(i,j,k)
           T_i = Tmp(i,j,k)
           p_i = prs(i,j,k)

           g = gcalc_tro(T_i,rho_i)
           !c_fact = p_i/T_i/rho_i
           cc2 = c2calc_tro(T_i,rho_i)
           cv = cvcalc_tro(T_i,rho_i)
           cp = cpcalc_tro(T_i,rho_i)
           s_i = scalc_tro(T_i,rho_i)

           qq = u_i**2 + v_i**2 + w_i**2

           e_i = rhoe_n(i,j,k)/rho_i - 0.5_wp*qq
           h_i = e_i + p_i/rho_i

           b1 = dwy(i,j,k) - dvz(i,j,k)
           b2 = duz(i,j,k) - dwx(i,j,k)
           b3 = dvx(i,j,k) - duy(i,j,k)
           divloc = dux(i,j,k) + dvy(i,j,k) + dwz(i,j,k)

           mu = visc(i,j,k)
           la = cok(i,j,k)

           prr = cp*mu/la
           eck = qq/(cp*T_i)

           M=sqrt(qq/cc2)

           ! total quantities
           s0 = scalc_tro(T_i,rho_i) ! useless, to be suppressed
           h0 = h_i + qq/2
           ! First guess using gamma equivalent
           g_eq = cp/cv
           ! T0 = T_i*(1 + (g_eq-1)/2 * M**2)                     ! Not working
           ! ro0 = rho_i*(1 + (g_eq-1)/2 * M**2)**(1/(g_eq-1))    ! Not working
           p0 = p_ref; e0 = 0.0_wp
           T0 = T_ref; ro0 = rho_ref
           !call stagnation_calc(h0,s0,e0,T0,p0,ro0)

           avg_s(i,j,  1)= avg_s(i,j,  1) + rho_i
           avg_s(i,j,  2)= avg_s(i,j,  2) + u_i
           avg_s(i,j,  3)= avg_s(i,j,  3) + v_i
           avg_s(i,j,  4)= avg_s(i,j,  4) + w_i
           avg_s(i,j,  5)= avg_s(i,j,  5) + p_i
           avg_s(i,j,  6)= avg_s(i,j,  6) + T_i
           avg_s(i,j,  7)= avg_s(i,j,  7) + rhou(i,j,k)
           avg_s(i,j,  8)= avg_s(i,j,  8) + rhov(i,j,k)
           avg_s(i,j,  9)= avg_s(i,j,  9) + rhow(i,j,k)
           avg_s(i,j, 10)= avg_s(i,j, 10) + rhoe(i,j,k)
           avg_s(i,j, 11)= avg_s(i,j, 11) + rho_i**2
           avg_s(i,j, 12)= avg_s(i,j, 12) + u_i**2
           avg_s(i,j, 13)= avg_s(i,j, 13) + v_i**2
           avg_s(i,j, 14)= avg_s(i,j, 14) + w_i**2
           avg_s(i,j, 15)= avg_s(i,j, 15) + u_i*v_i
           avg_s(i,j, 16)= avg_s(i,j, 16) + u_i*w_i
           avg_s(i,j, 17)= avg_s(i,j, 17) + v_i*w_i
           avg_s(i,j, 18)= avg_s(i,j, 18) + v_i*T_i
           avg_s(i,j, 19)= avg_s(i,j, 19) + p_i**2
           avg_s(i,j, 20)= avg_s(i,j, 20) + T_i**2
           avg_s(i,j, 21)= avg_s(i,j, 21) + mu
           avg_s(i,j, 22)= avg_s(i,j, 22) + divloc
           avg_s(i,j, 23)= avg_s(i,j, 23) + divloc**2

           avg_s(i,j, 24)= avg_s(i,j, 24) + e_i
           avg_s(i,j, 25)= avg_s(i,j, 25) + h_i
           avg_s(i,j, 26)= avg_s(i,j, 26) + sqrt(cc2)
           avg_s(i,j, 27)= avg_s(i,j, 27) + s_i
           avg_s(i,j, 28)= avg_s(i,j, 28) + M
           avg_s(i,j, 29)= avg_s(i,j, 29) + 0.5*qq
           avg_s(i,j, 30)= avg_s(i,j, 30) + g
           avg_s(i,j, 31)= avg_s(i,j, 31) + la
           avg_s(i,j, 32)= avg_s(i,j, 32) + cp
           avg_s(i,j, 33)= avg_s(i,j, 33) + cv
           avg_s(i,j, 34)= avg_s(i,j, 34) + prr
           avg_s(i,j, 35)= avg_s(i,j, 35) + eck
           avg_s(i,j, 36)= avg_s(i,j, 36) + rho_i*dux(i,j,k)
           avg_s(i,j, 37)= avg_s(i,j, 37) + rho_i*duy(i,j,k)
           avg_s(i,j, 38)= avg_s(i,j, 38) + rho_i*duz(i,j,k)
           avg_s(i,j, 39)= avg_s(i,j, 39) + rho_i*dvx(i,j,k)
           avg_s(i,j, 40)= avg_s(i,j, 40) + rho_i*dvy(i,j,k)
           avg_s(i,j, 41)= avg_s(i,j, 41) + rho_i*dvz(i,j,k)
           avg_s(i,j, 42)= avg_s(i,j, 42) + rho_i*dwx(i,j,k)
           avg_s(i,j, 43)= avg_s(i,j, 43) + rho_i*dwy(i,j,k)
           avg_s(i,j, 44)= avg_s(i,j, 44) + rho_i*dwz(i,j,k)
           avg_s(i,j, 45)= avg_s(i,j, 45) + p_i*divloc
           avg_s(i,j, 46)= avg_s(i,j, 46) + rho_i*divloc
           avg_s(i,j, 47)= avg_s(i,j, 47) + b1
           avg_s(i,j, 48)= avg_s(i,j, 48) + b2
           avg_s(i,j, 49)= avg_s(i,j, 49) + b3
           avg_s(i,j, 50)= avg_s(i,j, 50) + rho_i*T_i
           avg_s(i,j, 51)= avg_s(i,j, 51) + u_i*T_i
           avg_s(i,j, 52)= avg_s(i,j, 52) + v_i*T_i
           avg_s(i,j, 53)= avg_s(i,j, 53) + e_i**2
           avg_s(i,j, 54)= avg_s(i,j, 54) + h_i**2
           avg_s(i,j, 55)= avg_s(i,j, 55) + cc2
           avg_s(i,j, 56)= avg_s(i,j, 56) + s_i**2
           avg_s(i,j, 57)= avg_s(i,j, 57) + qq/cc2
           avg_s(i,j, 58)= avg_s(i,j, 58) + g**2
           avg_s(i,j, 59)= avg_s(i,j, 59) + mu**2
           avg_s(i,j, 60)= avg_s(i,j, 60) + la**2
           avg_s(i,j, 61)= avg_s(i,j, 61) + cv**2
           avg_s(i,j, 62)= avg_s(i,j, 62) + cp**2
           avg_s(i,j, 63)= avg_s(i,j, 63) + prr**2
           avg_s(i,j, 64)= avg_s(i,j, 64) + eck**2
           avg_s(i,j, 65)= avg_s(i,j, 65) + p_i*u_i
           avg_s(i,j, 66)= avg_s(i,j, 66) + p_i*v_i
           avg_s(i,j, 67)= avg_s(i,j, 67) + s_i*u_i
           avg_s(i,j, 68)= avg_s(i,j, 68) + s_i*v_i
           avg_s(i,j, 69)= avg_s(i,j, 69) + p_i*rho_i
           avg_s(i,j, 70)= avg_s(i,j, 70) + h_i*rho_i
           avg_s(i,j, 71)= avg_s(i,j, 71) + T_i*p_i
           avg_s(i,j, 72)= avg_s(i,j, 72) + p_i*s_i
           avg_s(i,j, 73)= avg_s(i,j, 73) + T_i*s_i
           avg_s(i,j, 74)= avg_s(i,j, 74) + rho_i*s_i
           avg_s(i,j, 75)= avg_s(i,j, 75) + g*rho_i
           avg_s(i,j, 76)= avg_s(i,j, 76) + g*p_i
           avg_s(i,j, 77)= avg_s(i,j, 77) + g*s_i
           avg_s(i,j, 78)= avg_s(i,j, 78) + g*T_i
           avg_s(i,j, 79)= avg_s(i,j, 79) + g*u_i
           avg_s(i,j, 80)= avg_s(i,j, 80) + g*v_i
           avg_s(i,j, 81)= avg_s(i,j, 81) + p_i*dux(i,j,k)
           avg_s(i,j, 82)= avg_s(i,j, 82) + p_i*dvy(i,j,k)
           avg_s(i,j, 83)= avg_s(i,j, 83) + p_i*dwz(i,j,k)
           avg_s(i,j, 84)= avg_s(i,j, 84) + p_i*duy(i,j,k)
           avg_s(i,j, 85)= avg_s(i,j, 85) + p_i*dvx(i,j,k)
           avg_s(i,j, 86)= avg_s(i,j, 86) + rho_i*divloc**2
           avg_s(i,j, 87)= avg_s(i,j, 87) + dux(i,j,k)**2
           avg_s(i,j, 88)= avg_s(i,j, 88) + duy(i,j,k)**2
           avg_s(i,j, 89)= avg_s(i,j, 89) + duz(i,j,k)**2
           avg_s(i,j, 90)= avg_s(i,j, 90) + dvx(i,j,k)**2
           avg_s(i,j, 91)= avg_s(i,j, 91) + dvy(i,j,k)**2
           avg_s(i,j, 92)= avg_s(i,j, 92) + dvz(i,j,k)**2
           avg_s(i,j, 93)= avg_s(i,j, 93) + dwx(i,j,k)**2
           avg_s(i,j, 94)= avg_s(i,j, 94) + dwy(i,j,k)**2
           avg_s(i,j, 95)= avg_s(i,j, 95) + dwz(i,j,k)**2
           avg_s(i,j, 96)= avg_s(i,j, 96) + b1**2
           avg_s(i,j, 97)= avg_s(i,j, 97) + b2**2
           avg_s(i,j, 98)= avg_s(i,j, 98) + b3**2
           avg_s(i,j, 99)= avg_s(i,j, 99) + rho_i*b1
           avg_s(i,j,100)= avg_s(i,j,100) + rho_i*b2
           avg_s(i,j,101)= avg_s(i,j,101) + rho_i*b3
           avg_s(i,j,102)= avg_s(i,j,102) + rho_i*u_i**2
           avg_s(i,j,103)= avg_s(i,j,103) + rho_i*v_i**2
           avg_s(i,j,104)= avg_s(i,j,104) + rho_i*w_i**2
           avg_s(i,j,105)= avg_s(i,j,105) + rho_i*T_i**2
           avg_s(i,j,106)= avg_s(i,j,106) + rho_i*b1**2
           avg_s(i,j,107)= avg_s(i,j,107) + rho_i*b2**2
           avg_s(i,j,108)= avg_s(i,j,108) + rho_i*b3**2
           avg_s(i,j,109)= avg_s(i,j,109) + rho_i*u_i*v_i
           avg_s(i,j,110)= avg_s(i,j,110) + rho_i*u_i*w_i
           avg_s(i,j,111)= avg_s(i,j,111) + rho_i*v_i*w_i
           avg_s(i,j,112)= avg_s(i,j,112) + rho_i*v_i*T_i
           avg_s(i,j,113)= avg_s(i,j,113) + rho_i*u_i**2*v_i
           avg_s(i,j,114)= avg_s(i,j,114) + rho_i*v_i**2*v_i
           avg_s(i,j,115)= avg_s(i,j,115) + rho_i*w_i**2*v_i
           avg_s(i,j,116)= avg_s(i,j,116) + rho_i*v_i**2*u_i
           avg_s(i,j,117)= avg_s(i,j,117) + rho_i*dux(i,j,k)**2
           avg_s(i,j,118)= avg_s(i,j,118) + rho_i*dvy(i,j,k)**2
           avg_s(i,j,119)= avg_s(i,j,119) + rho_i*dwz(i,j,k)**2
           avg_s(i,j,120)= avg_s(i,j,120) + rho_i*duy(i,j,k)*dvx(i,j,k)
           avg_s(i,j,121)= avg_s(i,j,121) + rho_i*duz(i,j,k)*dwx(i,j,k)
           avg_s(i,j,122)= avg_s(i,j,122) + rho_i*dvz(i,j,k)*dwy(i,j,k)
           avg_s(i,j,123)= avg_s(i,j,123) + u_i**3
           avg_s(i,j,124)= avg_s(i,j,124) + p_i**3
           avg_s(i,j,125)= avg_s(i,j,125) + u_i**4
           avg_s(i,j,126)= avg_s(i,j,126) + p_i**4
           avg_s(i,j,127)= avg_s(i,j,127) + Frhou(i,j,k)
           avg_s(i,j,128)= avg_s(i,j,128) + Frhov(i,j,k)
           avg_s(i,j,129)= avg_s(i,j,129) + Frhow(i,j,k)
           avg_s(i,j,130)= avg_s(i,j,130) + Grhov(i,j,k)
           avg_s(i,j,131)= avg_s(i,j,131) + Grhow(i,j,k)
           avg_s(i,j,132)= avg_s(i,j,132) + Hrhow(i,j,k)
           avg_s(i,j,133)= avg_s(i,j,133) + Frhov(i,j,k)*u_i
           avg_s(i,j,134)= avg_s(i,j,134) + Frhou(i,j,k)*u_i
           avg_s(i,j,135)= avg_s(i,j,135) + Frhov(i,j,k)*v_i
           avg_s(i,j,136)= avg_s(i,j,136) + Frhow(i,j,k)*w_i
           avg_s(i,j,137)= avg_s(i,j,137) + Grhov(i,j,k)*u_i
           avg_s(i,j,138)= avg_s(i,j,138) + Grhov(i,j,k)*v_i
           avg_s(i,j,139)= avg_s(i,j,139) + Grhow(i,j,k)*w_i
           avg_s(i,j,140)= avg_s(i,j,140) + Frhou(i,j,k)*dux(i,j,k)
           avg_s(i,j,141)= avg_s(i,j,141) + Frhou(i,j,k)*dvx(i,j,k)
           avg_s(i,j,142)= avg_s(i,j,142) + Frhov(i,j,k)*dux(i,j,k)
           avg_s(i,j,143)= avg_s(i,j,143) + Frhov(i,j,k)*duy(i,j,k)
           avg_s(i,j,144)= avg_s(i,j,144) + Frhov(i,j,k)*dvx(i,j,k)
           avg_s(i,j,145)= avg_s(i,j,145) + Frhov(i,j,k)*dvy(i,j,k)
           avg_s(i,j,146)= avg_s(i,j,146) + Frhow(i,j,k)*duz(i,j,k)
           avg_s(i,j,147)= avg_s(i,j,147) + Frhow(i,j,k)*dvz(i,j,k)
           avg_s(i,j,148)= avg_s(i,j,148) + Frhow(i,j,k)*dwx(i,j,k)
           avg_s(i,j,149)= avg_s(i,j,149) + Grhov(i,j,k)*duy(i,j,k)
           avg_s(i,j,150)= avg_s(i,j,150) + Grhov(i,j,k)*dvy(i,j,k)
           avg_s(i,j,151)= avg_s(i,j,151) + Grhow(i,j,k)*duz(i,j,k)
           avg_s(i,j,152)= avg_s(i,j,152) + Grhow(i,j,k)*dvz(i,j,k)
           avg_s(i,j,153)= avg_s(i,j,153) + Grhow(i,j,k)*dwy(i,j,k)
           avg_s(i,j,154)= avg_s(i,j,154) + Hrhow(i,j,k)*dwz(i,j,k)
           avg_s(i,j,155)= avg_s(i,j,155) + la*dTx(i,j,k)
           avg_s(i,j,156)= avg_s(i,j,156) + la*dTy(i,j,k)
           avg_s(i,j,157)= avg_s(i,j,157) + la*dTz(i,j,k)
           avg_s(i,j,158)= avg_s(i,j,158) + h_i*u_i
           avg_s(i,j,159)= avg_s(i,j,159) + h_i*v_i
           avg_s(i,j,160)= avg_s(i,j,160) + h_i*w_i
           avg_s(i,j,161)= avg_s(i,j,161) + rho_i*h_i*u_i
           avg_s(i,j,162)= avg_s(i,j,162) + rho_i*h_i*v_i
           avg_s(i,j,163)= avg_s(i,j,163) + rho_i*h_i*w_i
           avg_s(i,j,164)= avg_s(i,j,164) + rho_i*s_i*u_i
           avg_s(i,j,165)= avg_s(i,j,165) + rho_i*s_i*v_i
           avg_s(i,j,166)= avg_s(i,j,166) + rho_i*s_i*w_i
           avg_s(i,j,167)= avg_s(i,j,167) + rho_i*u_i**3
           avg_s(i,j,168)= avg_s(i,j,168) + rho_i*v_i**3
           avg_s(i,j,169)= avg_s(i,j,169) + rho_i*w_i**3
           avg_s(i,j,170)= avg_s(i,j,170) + rho_i*w_i**2*u_i

           ! Total quantities
           avg_s(i,j,171)= avg_s(i,j,171) + h0
           avg_s(i,j,172)= avg_s(i,j,172) + e0
           avg_s(i,j,173)= avg_s(i,j,173) + s0
           avg_s(i,j,174)= avg_s(i,j,174) + T0
           avg_s(i,j,175)= avg_s(i,j,175) + p0
           avg_s(i,j,176)= avg_s(i,j,176) + ro0

           if (is_RANS) then
              nu  = visc(i,j,k)/rho(i,j,k)
              khi = nutil(i,j,k)/nu
              fnu1 = khi**3/(khi**3+cnu13)
              mut(i,j,k) = nutil(i,j,k)*fnu1*rho_i
              avg_s(i,j,177)= avg_s(i,j,177) + mut(i,j,k)

              ! avg_s(i,j,)= avg_s(i,j,) + mu*dux(i,j,k)
              ! avg_s(i,j,)= avg_s(i,j,) + mu*duy(i,j,k)
              ! avg_s(i,j,)= avg_s(i,j,) + mu*dvx(i,j,k)
              ! avg_s(i,j,)= avg_s(i,j,) + mu*dvy(i,j,k)
              ! avg_s(i,j,)= avg_s(i,j,) + mut(i,j,k)
              ! avg_s(i,j,)= avg_s(i,j,) + mut(i,j,k)*dux(i,j,k)
              ! avg_s(i,j,)= avg_s(i,j,) + mut(i,j,k)*duy(i,j,k)
              ! avg_s(i,j,)= avg_s(i,j,) + mut(i,j,k)*dvx(i,j,k)
              ! avg_s(i,j,)= avg_s(i,j,) + mut(i,j,k)*dvy(i,j,k)
           else
              avg_s(i,j,177)= avg_s(i,j,177) + 0.0_wp
           endif

        enddo
     enddo
  enddo

  avg_s = avg_s/dble(nz)

  avg_t = ((nn1-1.0_wp)*avg_t + avg_s)*inn1

end subroutine stats_turb


!===============================================================================
subroutine stats_turb3
!===============================================================================
  !> Compute stats for turbine flow
!===============================================================================
  use mod_mpi
  use mod_eos
  use mod_flow
  use mod_time
  use mod_constant
  use mod_rans
  implicit none
  ! ----------------------------------------------------------------------------
  integer  :: i,j,k
  real(wp) :: nn1,nm1,inn1,inm1,cc2,qq,cnu13
  real(wp) :: g,s_i,divloc,cp,cv,p_i,T_i,rho_i
  real(wp) :: u_i,v_i,w_i,b1,b2,b3,mu,la,prr,eck,e_i,h_i,nu!,c_fact
  real(wp) :: h0,e0,s0,T0,p0,ro0,M,g_eq
  real(wp) :: taup1,taup2,taup3
  real(wp) :: tauv1,tauv2,tauv3
  ! ----------------------------------------------------------------------------
  integer :: ndeb_reprise

  ! Collect statistics
  ! ------------------
  if (ntotal.lt.ndeb) return
  nn1 = dble(int((ntotal-ndeb)/freq_stats)) + 1
  inn1= 1.0_wp/nn1

  ! calcul LS59 rough
  ndeb_reprise=550001
  nm1=dble(ntotal+1-ndeb_reprise)/dble(freq_stats)
  inm1= 1.0_wp/nm1
  nm1=nm1-1.0_wp

  ! Spatial local and global averages at this iteration
  avg_s  = 0.0_wp

  cnu1=7.1_wp
  cnu13=cnu1**3

  do k=1,nz
     do j=1,ny
        do i=1,nx
           ! if ((i.ne.20).and.(j.ne.1)) return

           rho_i = rho(i,j,k)

           u_i = uu(i,j,k)
           v_i = vv(i,j,k)
           w_i = ww(i,j,k)
           T_i = Tmp(i,j,k)
           p_i = prs(i,j,k)

           g = gcalc_tro(T_i,rho_i)
           !c_fact = p_i/T_i/rho_i
           cc2 = c2calc_tro(T_i,rho_i)
           cv = cvcalc_tro(T_i,rho_i)
           cp = cpcalc_tro(T_i,rho_i)
           s_i = scalc_tro(T_i,rho_i)

           qq = u_i**2 + v_i**2 + w_i**2

           e_i = rhoe_n(i,j,k)/rho_i - 0.5_wp*qq
           h_i = e_i + p_i/rho_i

           b1 = dwy(i,j,k) - dvz(i,j,k)
           b2 = duz(i,j,k) - dwx(i,j,k)
           b3 = dvx(i,j,k) - duy(i,j,k)
           divloc = dux(i,j,k) + dvy(i,j,k) + dwz(i,j,k)

           mu = visc(i,j,k)
           la = cok(i,j,k)

           prr = cp*mu/la
           eck = qq/(cp*T_i)

           M=sqrt(qq/cc2)

           ! total quantities
           s0 = scalc_tro(T_i,rho_i) ! useless, to be suppressed
           h0 = h_i + qq/2
           ! First guess using gamma equivalent
           g_eq = cp/cv
           ! T0 = T_i*(1 + (g_eq-1)/2 * M**2)                     ! Not working
           ! ro0 = rho_i*(1 + (g_eq-1)/2 * M**2)**(1/(g_eq-1))    ! Not working
           p0 = p_ref; e0 = 0.0_wp
           T0 = T_ref; ro0 = rho_ref
           !call stagnation_calc(h0,s0,e0,T0,p0,ro0)

           ! wall shear stress
           ! /!\ ONLY implemented for walls at jmin
           if ((is_bc_wall(2,1)).and.(j==1)) then

              ! contribution of pressure to the shear stress
              taup1= prs(i,1,k)*nxn_jmin(i,k)
              taup2= prs(i,1,k)*nyn_jmin(i,k)
              taup3= prs(i,1,k)*nzn_jmin(i,k)

              ! contribution of viscous stresses to the shear stress
              ! tau.nx=tau11*eta_x+tau12*eta_y+tau13*eta_z -> Grhou in visc.
              ! tau.ny=tau12*eta_x+tau22*eta_y+tau23*eta_z -> Grhov in visc.
              ! tau.nz=tau13*eta_x+tau23*eta_y+tau33*eta_z -> Grhow in visc.
              tauv1= Grhou(i,1,k)*ineta_jmin(i,k)
              tauv2= Grhov(i,1,k)*ineta_jmin(i,k)
              tauv3= Grhow(i,1,k)*ineta_jmin(i,k)

              avg_w(i,k,1)= (nm1*avg_w(i,k,1) + taup1)*inm1
              avg_w(i,k,2)= (nm1*avg_w(i,k,2) + taup2)*inm1
              avg_w(i,k,3)= (nm1*avg_w(i,k,3) + taup3)*inm1
              avg_w(i,k,4)= (nm1*avg_w(i,k,4) + tauv1)*inm1
              avg_w(i,k,5)= (nm1*avg_w(i,k,5) + tauv2)*inm1
              avg_w(i,k,6)= (nm1*avg_w(i,k,6) + tauv3)*inm1

           endif

           s0 = scalc_tro(T_i,rho_i) ! useless, to be suppressed

           avg_v(i,j,k, 1)= (nm1*avg_v(i,j,k, 1) + rho_i)*inm1
           avg_v(i,j,k, 2)= (nm1*avg_v(i,j,k, 2) + u_i)*inm1
           avg_v(i,j,k, 3)= (nm1*avg_v(i,j,k, 3) + v_i)*inm1
           avg_v(i,j,k, 4)= (nm1*avg_v(i,j,k, 4) + w_i)*inm1
           avg_v(i,j,k, 5)= (nm1*avg_v(i,j,k, 5) + p_i)*inm1
           avg_v(i,j,k, 6)= (nm1*avg_v(i,j,k, 6) + rho_i**2)*inm1
           avg_v(i,j,k, 7)= (nm1*avg_v(i,j,k, 7) + u_i**2)*inm1
           avg_v(i,j,k, 8)= (nm1*avg_v(i,j,k, 8) + v_i**2)*inm1
           avg_v(i,j,k, 9)= (nm1*avg_v(i,j,k, 9) + w_i**2)*inm1
           avg_v(i,j,k,10)= (nm1*avg_v(i,j,k,10) + u_i*v_i)*inm1
           avg_v(i,j,k,11)= (nm1*avg_v(i,j,k,11) + u_i*w_i)*inm1
           avg_v(i,j,k,12)= (nm1*avg_v(i,j,k,12) + v_i*w_i)*inm1
           avg_v(i,j,k,13)= (nm1*avg_v(i,j,k,13) + p_i**2)*inm1

           avg_s(i,j,  1)= avg_s(i,j,  1) + rho_i
           avg_s(i,j,  2)= avg_s(i,j,  2) + u_i
           avg_s(i,j,  3)= avg_s(i,j,  3) + v_i
           avg_s(i,j,  4)= avg_s(i,j,  4) + w_i
           avg_s(i,j,  5)= avg_s(i,j,  5) + p_i
           avg_s(i,j,  6)= avg_s(i,j,  6) + T_i
           avg_s(i,j,  7)= avg_s(i,j,  7) + rhou(i,j,k)
           avg_s(i,j,  8)= avg_s(i,j,  8) + rhov(i,j,k)
           avg_s(i,j,  9)= avg_s(i,j,  9) + rhow(i,j,k)
           avg_s(i,j, 10)= avg_s(i,j, 10) + rhoe(i,j,k)
           avg_s(i,j, 11)= avg_s(i,j, 11) + rho_i**2
           avg_s(i,j, 12)= avg_s(i,j, 12) + u_i**2
           avg_s(i,j, 13)= avg_s(i,j, 13) + v_i**2
           avg_s(i,j, 14)= avg_s(i,j, 14) + w_i**2
           avg_s(i,j, 15)= avg_s(i,j, 15) + u_i*v_i
           avg_s(i,j, 16)= avg_s(i,j, 16) + u_i*w_i
           avg_s(i,j, 17)= avg_s(i,j, 17) + v_i*w_i
           avg_s(i,j, 18)= avg_s(i,j, 18) + v_i*T_i
           avg_s(i,j, 19)= avg_s(i,j, 19) + p_i**2
           avg_s(i,j, 20)= avg_s(i,j, 20) + T_i**2
           avg_s(i,j, 21)= avg_s(i,j, 21) + mu
           avg_s(i,j, 22)= avg_s(i,j, 22) + divloc
           avg_s(i,j, 23)= avg_s(i,j, 23) + divloc**2

           avg_s(i,j, 24)= avg_s(i,j, 24) + e_i
           avg_s(i,j, 25)= avg_s(i,j, 25) + h_i
           avg_s(i,j, 26)= avg_s(i,j, 26) + sqrt(cc2)
           avg_s(i,j, 27)= avg_s(i,j, 27) + s_i
           avg_s(i,j, 28)= avg_s(i,j, 28) + M
           avg_s(i,j, 29)= avg_s(i,j, 29) + 0.5*qq
           avg_s(i,j, 30)= avg_s(i,j, 30) + g
           avg_s(i,j, 31)= avg_s(i,j, 31) + la
           avg_s(i,j, 32)= avg_s(i,j, 32) + cp
           avg_s(i,j, 33)= avg_s(i,j, 33) + cv
           avg_s(i,j, 34)= avg_s(i,j, 34) + prr
           avg_s(i,j, 35)= avg_s(i,j, 35) + eck
           avg_s(i,j, 36)= avg_s(i,j, 36) + rho_i*dux(i,j,k)
           avg_s(i,j, 37)= avg_s(i,j, 37) + rho_i*duy(i,j,k)
           avg_s(i,j, 38)= avg_s(i,j, 38) + rho_i*duz(i,j,k)
           avg_s(i,j, 39)= avg_s(i,j, 39) + rho_i*dvx(i,j,k)
           avg_s(i,j, 40)= avg_s(i,j, 40) + rho_i*dvy(i,j,k)
           avg_s(i,j, 41)= avg_s(i,j, 41) + rho_i*dvz(i,j,k)
           avg_s(i,j, 42)= avg_s(i,j, 42) + rho_i*dwx(i,j,k)
           avg_s(i,j, 43)= avg_s(i,j, 43) + rho_i*dwy(i,j,k)
           avg_s(i,j, 44)= avg_s(i,j, 44) + rho_i*dwz(i,j,k)
           avg_s(i,j, 45)= avg_s(i,j, 45) + p_i*divloc
           avg_s(i,j, 46)= avg_s(i,j, 46) + rho_i*divloc
           avg_s(i,j, 47)= avg_s(i,j, 47) + b1
           avg_s(i,j, 48)= avg_s(i,j, 48) + b2
           avg_s(i,j, 49)= avg_s(i,j, 49) + b3
           avg_s(i,j, 50)= avg_s(i,j, 50) + rho_i*T_i
           avg_s(i,j, 51)= avg_s(i,j, 51) + u_i*T_i
           avg_s(i,j, 52)= avg_s(i,j, 52) + v_i*T_i
           avg_s(i,j, 53)= avg_s(i,j, 53) + e_i**2
           avg_s(i,j, 54)= avg_s(i,j, 54) + h_i**2
           avg_s(i,j, 55)= avg_s(i,j, 55) + cc2
           avg_s(i,j, 56)= avg_s(i,j, 56) + s_i**2
           avg_s(i,j, 57)= avg_s(i,j, 57) + qq/cc2
           avg_s(i,j, 58)= avg_s(i,j, 58) + g**2
           avg_s(i,j, 59)= avg_s(i,j, 59) + mu**2
           avg_s(i,j, 60)= avg_s(i,j, 60) + la**2
           avg_s(i,j, 61)= avg_s(i,j, 61) + cv**2
           avg_s(i,j, 62)= avg_s(i,j, 62) + cp**2
           avg_s(i,j, 63)= avg_s(i,j, 63) + prr**2
           avg_s(i,j, 64)= avg_s(i,j, 64) + eck**2
           avg_s(i,j, 65)= avg_s(i,j, 65) + p_i*u_i
           avg_s(i,j, 66)= avg_s(i,j, 66) + p_i*v_i
           avg_s(i,j, 67)= avg_s(i,j, 67) + s_i*u_i
           avg_s(i,j, 68)= avg_s(i,j, 68) + s_i*v_i
           avg_s(i,j, 69)= avg_s(i,j, 69) + p_i*rho_i
           avg_s(i,j, 70)= avg_s(i,j, 70) + h_i*rho_i
           avg_s(i,j, 71)= avg_s(i,j, 71) + T_i*p_i
           avg_s(i,j, 72)= avg_s(i,j, 72) + p_i*s_i
           avg_s(i,j, 73)= avg_s(i,j, 73) + T_i*s_i
           avg_s(i,j, 74)= avg_s(i,j, 74) + rho_i*s_i
           avg_s(i,j, 75)= avg_s(i,j, 75) + g*rho_i
           avg_s(i,j, 76)= avg_s(i,j, 76) + g*p_i
           avg_s(i,j, 77)= avg_s(i,j, 77) + g*s_i
           avg_s(i,j, 78)= avg_s(i,j, 78) + g*T_i
           avg_s(i,j, 79)= avg_s(i,j, 79) + g*u_i
           avg_s(i,j, 80)= avg_s(i,j, 80) + g*v_i
           avg_s(i,j, 81)= avg_s(i,j, 81) + p_i*dux(i,j,k)
           avg_s(i,j, 82)= avg_s(i,j, 82) + p_i*dvy(i,j,k)
           avg_s(i,j, 83)= avg_s(i,j, 83) + p_i*dwz(i,j,k)
           avg_s(i,j, 84)= avg_s(i,j, 84) + p_i*duy(i,j,k)
           avg_s(i,j, 85)= avg_s(i,j, 85) + p_i*dvx(i,j,k)
           avg_s(i,j, 86)= avg_s(i,j, 86) + rho_i*divloc**2
           avg_s(i,j, 87)= avg_s(i,j, 87) + dux(i,j,k)**2
           avg_s(i,j, 88)= avg_s(i,j, 88) + duy(i,j,k)**2
           avg_s(i,j, 89)= avg_s(i,j, 89) + duz(i,j,k)**2
           avg_s(i,j, 90)= avg_s(i,j, 90) + dvx(i,j,k)**2
           avg_s(i,j, 91)= avg_s(i,j, 91) + dvy(i,j,k)**2
           avg_s(i,j, 92)= avg_s(i,j, 92) + dvz(i,j,k)**2
           avg_s(i,j, 93)= avg_s(i,j, 93) + dwx(i,j,k)**2
           avg_s(i,j, 94)= avg_s(i,j, 94) + dwy(i,j,k)**2
           avg_s(i,j, 95)= avg_s(i,j, 95) + dwz(i,j,k)**2
           avg_s(i,j, 96)= avg_s(i,j, 96) + b1**2
           avg_s(i,j, 97)= avg_s(i,j, 97) + b2**2
           avg_s(i,j, 98)= avg_s(i,j, 98) + b3**2
           avg_s(i,j, 99)= avg_s(i,j, 99) + rho_i*b1
           avg_s(i,j,100)= avg_s(i,j,100) + rho_i*b2
           avg_s(i,j,101)= avg_s(i,j,101) + rho_i*b3
           avg_s(i,j,102)= avg_s(i,j,102) + rho_i*u_i**2
           avg_s(i,j,103)= avg_s(i,j,103) + rho_i*v_i**2
           avg_s(i,j,104)= avg_s(i,j,104) + rho_i*w_i**2
           avg_s(i,j,105)= avg_s(i,j,105) + rho_i*T_i**2
           avg_s(i,j,106)= avg_s(i,j,106) + rho_i*b1**2
           avg_s(i,j,107)= avg_s(i,j,107) + rho_i*b2**2
           avg_s(i,j,108)= avg_s(i,j,108) + rho_i*b3**2
           avg_s(i,j,109)= avg_s(i,j,109) + rho_i*u_i*v_i
           avg_s(i,j,110)= avg_s(i,j,110) + rho_i*u_i*w_i
           avg_s(i,j,111)= avg_s(i,j,111) + rho_i*v_i*w_i
           avg_s(i,j,112)= avg_s(i,j,112) + rho_i*v_i*T_i
           avg_s(i,j,113)= avg_s(i,j,113) + rho_i*u_i**2*v_i
           avg_s(i,j,114)= avg_s(i,j,114) + rho_i*v_i**2*v_i
           avg_s(i,j,115)= avg_s(i,j,115) + rho_i*w_i**2*v_i
           avg_s(i,j,116)= avg_s(i,j,116) + rho_i*v_i**2*u_i
           avg_s(i,j,117)= avg_s(i,j,117) + rho_i*dux(i,j,k)**2
           avg_s(i,j,118)= avg_s(i,j,118) + rho_i*dvy(i,j,k)**2
           avg_s(i,j,119)= avg_s(i,j,119) + rho_i*dwz(i,j,k)**2
           avg_s(i,j,120)= avg_s(i,j,120) + rho_i*duy(i,j,k)*dvx(i,j,k)
           avg_s(i,j,121)= avg_s(i,j,121) + rho_i*duz(i,j,k)*dwx(i,j,k)
           avg_s(i,j,122)= avg_s(i,j,122) + rho_i*dvz(i,j,k)*dwy(i,j,k)
           avg_s(i,j,123)= avg_s(i,j,123) + u_i**3
           avg_s(i,j,124)= avg_s(i,j,124) + p_i**3
           avg_s(i,j,125)= avg_s(i,j,125) + u_i**4
           avg_s(i,j,126)= avg_s(i,j,126) + p_i**4
           avg_s(i,j,127)= avg_s(i,j,127) + Frhou(i,j,k)
           avg_s(i,j,128)= avg_s(i,j,128) + Frhov(i,j,k)
           avg_s(i,j,129)= avg_s(i,j,129) + Frhow(i,j,k)
           avg_s(i,j,130)= avg_s(i,j,130) + Grhov(i,j,k)
           avg_s(i,j,131)= avg_s(i,j,131) + Grhow(i,j,k)
           avg_s(i,j,132)= avg_s(i,j,132) + Hrhow(i,j,k)
           avg_s(i,j,133)= avg_s(i,j,133) + Frhov(i,j,k)*u_i
           avg_s(i,j,134)= avg_s(i,j,134) + Frhou(i,j,k)*u_i
           avg_s(i,j,135)= avg_s(i,j,135) + Frhov(i,j,k)*v_i
           avg_s(i,j,136)= avg_s(i,j,136) + Frhow(i,j,k)*w_i
           avg_s(i,j,137)= avg_s(i,j,137) + Grhov(i,j,k)*u_i
           avg_s(i,j,138)= avg_s(i,j,138) + Grhov(i,j,k)*v_i
           avg_s(i,j,139)= avg_s(i,j,139) + Grhow(i,j,k)*w_i
           avg_s(i,j,140)= avg_s(i,j,140) + Frhou(i,j,k)*dux(i,j,k)
           avg_s(i,j,141)= avg_s(i,j,141) + Frhou(i,j,k)*dvx(i,j,k)
           avg_s(i,j,142)= avg_s(i,j,142) + Frhov(i,j,k)*dux(i,j,k)
           avg_s(i,j,143)= avg_s(i,j,143) + Frhov(i,j,k)*duy(i,j,k)
           avg_s(i,j,144)= avg_s(i,j,144) + Frhov(i,j,k)*dvx(i,j,k)
           avg_s(i,j,145)= avg_s(i,j,145) + Frhov(i,j,k)*dvy(i,j,k)
           avg_s(i,j,146)= avg_s(i,j,146) + Frhow(i,j,k)*duz(i,j,k)
           avg_s(i,j,147)= avg_s(i,j,147) + Frhow(i,j,k)*dvz(i,j,k)
           avg_s(i,j,148)= avg_s(i,j,148) + Frhow(i,j,k)*dwx(i,j,k)
           avg_s(i,j,149)= avg_s(i,j,149) + Grhov(i,j,k)*duy(i,j,k)
           avg_s(i,j,150)= avg_s(i,j,150) + Grhov(i,j,k)*dvy(i,j,k)
           avg_s(i,j,151)= avg_s(i,j,151) + Grhow(i,j,k)*duz(i,j,k)
           avg_s(i,j,152)= avg_s(i,j,152) + Grhow(i,j,k)*dvz(i,j,k)
           avg_s(i,j,153)= avg_s(i,j,153) + Grhow(i,j,k)*dwy(i,j,k)
           avg_s(i,j,154)= avg_s(i,j,154) + Hrhow(i,j,k)*dwz(i,j,k)
           avg_s(i,j,155)= avg_s(i,j,155) + la*dTx(i,j,k)
           avg_s(i,j,156)= avg_s(i,j,156) + la*dTy(i,j,k)
           avg_s(i,j,157)= avg_s(i,j,157) + la*dTz(i,j,k)
           avg_s(i,j,158)= avg_s(i,j,158) + h_i*u_i
           avg_s(i,j,159)= avg_s(i,j,159) + h_i*v_i
           avg_s(i,j,160)= avg_s(i,j,160) + h_i*w_i
           avg_s(i,j,161)= avg_s(i,j,161) + rho_i*h_i*u_i
           avg_s(i,j,162)= avg_s(i,j,162) + rho_i*h_i*v_i
           avg_s(i,j,163)= avg_s(i,j,163) + rho_i*h_i*w_i
           avg_s(i,j,164)= avg_s(i,j,164) + rho_i*s_i*u_i
           avg_s(i,j,165)= avg_s(i,j,165) + rho_i*s_i*v_i
           avg_s(i,j,166)= avg_s(i,j,166) + rho_i*s_i*w_i
           avg_s(i,j,167)= avg_s(i,j,167) + rho_i*u_i**3
           avg_s(i,j,168)= avg_s(i,j,168) + rho_i*v_i**3
           avg_s(i,j,169)= avg_s(i,j,169) + rho_i*w_i**3
           avg_s(i,j,170)= avg_s(i,j,170) + rho_i*w_i**2*u_i

           ! Total quantities
           avg_s(i,j,171)= avg_s(i,j,171) + h0
           avg_s(i,j,172)= avg_s(i,j,172) + e0
           avg_s(i,j,173)= avg_s(i,j,173) + s0
           avg_s(i,j,174)= avg_s(i,j,174) + T0
           avg_s(i,j,175)= avg_s(i,j,175) + p0
           avg_s(i,j,176)= avg_s(i,j,176) + ro0

           if (is_RANS) then
              nu  = visc(i,j,k)/rho(i,j,k)
              khi = nutil(i,j,k)/nu
              fnu1 = khi**3/(khi**3+cnu13)
              mut(i,j,k) = nutil(i,j,k)*fnu1*rho_i
              avg_s(i,j,177)= avg_s(i,j,177) + mut(i,j,k)

              ! avg_s(i,j,)= avg_s(i,j,) + mu*dux(i,j,k)
              ! avg_s(i,j,)= avg_s(i,j,) + mu*duy(i,j,k)
              ! avg_s(i,j,)= avg_s(i,j,) + mu*dvx(i,j,k)
              ! avg_s(i,j,)= avg_s(i,j,) + mu*dvy(i,j,k)
              ! avg_s(i,j,)= avg_s(i,j,) + mut(i,j,k)
              ! avg_s(i,j,)= avg_s(i,j,) + mut(i,j,k)*dux(i,j,k)
              ! avg_s(i,j,)= avg_s(i,j,) + mut(i,j,k)*duy(i,j,k)
              ! avg_s(i,j,)= avg_s(i,j,) + mut(i,j,k)*dvx(i,j,k)
              ! avg_s(i,j,)= avg_s(i,j,) + mut(i,j,k)*dvy(i,j,k)
           else
              avg_s(i,j,177)= avg_s(i,j,177) + 0.0_wp
           endif

        enddo
     enddo
  enddo

  avg_s = avg_s/dble(nz)

  avg_t = ((nn1-1.0_wp)*avg_t + avg_s)*inn1

end subroutine stats_turb3
