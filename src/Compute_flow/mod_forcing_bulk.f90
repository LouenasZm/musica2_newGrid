!=================================================================================
module mod_forcing_bulk
!=================================================================================
  !> Module to define forcing term to enforce mass flow rate
!=================================================================================
  use precision
  use mod_mpi
  use mod_constant
  use mod_flow
  use mod_time
  implicit none
  !-------------------------------------------------------------------------------
  logical  :: is_forcing_bulk
  ! Bulk quantities
  real(wp) :: rhou_bulk  
  ! Forcing terms
  real(wp) :: forc_rhou,forc_rho
  real(wp) :: rhob_TRG,rhoub_TRG,rhob_SCH,rhoub_SCH
  real(wp) :: Q_n,Q_i
  real(wp), private :: alpha_f,beta_f
  !-------------------------------------------------------------------------------
  ! Inlet surface (section)
  real(wp), private :: Surf_in
  ! Flow volume (for periodic hill flow 04/03/21)
  real(wp), private :: Vol
  ! integration elements for inlet
  real(wp), dimension(:), allocatable, private :: dyi,dzi 
  ! dimensions of inlet section
  real(wp), private :: Lx,Ly,Lz
  
contains
  
  !===============================================================================
  subroutine init_forcing_bulk
  !===============================================================================
    !> Initializations for forcing mass flow rate
    !> for channel flow and periodic hill flow
  !===============================================================================
    use mod_bc_periodicity ! <~ needed for Lxp  TO BE CHANGED
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: j,k
    real(wp), dimension(:), allocatable :: dygi,dzgi ! integration elements on global grid
    ! ----------------------------------------------------------------------------
    integer :: i
real(wp) rou
    ! Integration elements at inlet (to compute bulk quantities)
    ! =============================
    allocate(dygi(ngy))
    if (is_curv) then
       dygi(1)=(ygc(1,2)-ygc(1,1))*0.5_wp
       do j=2,ngy-1
          dygi(j)=(ygc(1,j+1)-ygc(1,j-1))*0.5_wp
       enddo
       dygi(ngy)=(ygc(1,ngy)-ygc(1,ngy-1))*0.5_wp
    else
       dygi(1)=(yg(2)-yg(1))*0.5_wp
       do j=2,ngy-1
          dygi(j)=(yg(j+1)-yg(j-1))*0.5_wp
       enddo
       dygi(ngy)=(yg(ngy)-yg(ngy-1))*0.5_wp
    endif

    allocate(dzgi(ngz))
    dzgi=deltaz

    ! Partitioning of integration elements
    ! ====================================
    if (allocated(dyi)) deallocate(dyi)
    allocate(dyi(ny))
    do j=1,ny
       dyi(j)=dygi(j+coord(2)*ny)
    enddo

    if (allocated(dzi)) deallocate(dzi)
    allocate(dzi(nz))
    do k=1,nz
       dzi(k)=dzgi(k+coord(3)*nz)
    enddo

    deallocate(dygi,dzgi)
   
    ! Dimensions of inlet section
    ! ===========================
    if (iproc==0) then
       Lx=Lxp
       if (is_curv) then
          Ly=ygc(1,ngy)-ygc(1,1)
       else
          Ly=yg(ngy)-yg(1)
       endif
       if (is_2D) then
          Lz=1.0_wp
       else   
          Lz=zg(ngz)-zg(1)+deltaz
          !if (periods(3)==.true.) then
          !   Lz=zg(ngz)-zg(1)+zg(2)-zg(1) ! <- +deltaz
          !else
          !   Lz=zg(ngz)-zg(1)
          !endif
       endif
    endif
    call MPI_BCAST(Lx,1,MPI_DOUBLE_PRECISION,0,COMM_global,info)
    call MPI_BCAST(Ly,1,MPI_DOUBLE_PRECISION,0,COMM_global,info)
    call MPI_BCAST(Lz,1,MPI_DOUBLE_PRECISION,0,COMM_global,info)
    if ((iproc==0).and.(is_curv)) then    
       print *,'Lx forcing',Lx,9.*L_ref
       print *,'Ly forcing',Ly,2.035*L_ref
       print *,'Lz forcing',Lz,ngz*deltaz
    endif
    
    ! Fluid volume
    ! ============
    if (PHILL) then
       Vol=25.403009038579210*Lz*hc**2 ! specific to periodic hill flow case
    else
       Vol=Lx*Ly*Lz
    end if
    
    ! Inlet section
    ! =============
    Surf_in=Ly*Lz

    if (is_2D) then
       rhou_bulk=rho_ref*u_ref/2.
       rhob_TRG =rho_ref      ! This is the target bulk density
       rhoub_TRG=rhou_bulk/2. ! This is the target bulk momentum
    else
       rhou_bulk=rho_ref*u_ref
       rhob_TRG =rho_ref     ! This is the target bulk density
       rhoub_TRG=rhou_bulk   ! This is the target bulk momentum
    endif
 
    ! Init forcing
    ! ============
    ! initial forcing term
    forc_rhou=0.0_wp
    ! initial mass flow rate
    Q_i=rhou_bulk*Surf_in
    if (iproc.eq.0) print *,Q_i,rhou_bulk
    
    ! Computation of <rhou>_bulk
    ! ==========================
    rou = 0.0_wp
    do k=1,nz
       do j=1,ny
          do i=1,nx
             rou=rou+rhou_n(i,j,k)*dyi(j)*dzi(k)
          enddo
       enddo
    enddo
    call MPI_ALLREDUCE(MPI_IN_PLACE,rou,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
    rhou_bulk=rou/(Ly*ngx*Lz)

    
   if (iproc.eq.0) print *,'int',rhou_bulk



!stop
    
    Q_n=Q_i
    ! coefficients for updating the forcing term
    alpha_f= 2.0_wp/deltat/Surf_in
    beta_f =-0.2_wp/deltat/Surf_in
 
    ! Need mean temporal field for periodic hill
    ! ==========================================
    !if (is_curv) is_mean0=.true.
    if (PHILL) is_mean0=.true.
    
    ! Initialization of the forcing term ! TO BE CHANGED
    ! ==================================
    forc_rhou=utau**2/hc
    if (is_2D) forc_rhou=-1.5*mu_ref*u_ref/hc**2

    if (iproc.eq.0) then
       open(98, file='convqp.dat', form='formatted', status='replace')
       close(98)
       write(*,*) 'Initial forcing:', forc_rhou
       write(*,*) 'Initialisation OK'
    endif
    
  end subroutine init_forcing_bulk

  !===============================================================================
  subroutine forcing_rhou
  !===============================================================================
    !> Compute forcing term for plane channel flow
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,k
    real(wp) :: Q_np1,gn,Q_in
    ! compute dudy
    real(wp) :: Uav3_jmin,Uav2_jmin,Uavg3_jmin,Uavg2_jmin
    real(wp) :: Uav3_jmax,Uav2_jmax,Uavg3_jmax,Uavg2_jmax
    real(wp) :: dudy,dudy_jmin,dudy_jmax,idyw
    real(wp) :: dudy_meang_jmin,dudy_meang_jmax,dudy_mean_jmin,dudy_mean_jmax
    ! ----------------------------------------------------------------------------

    if (is_wall_model) then
       ! Compute dudy on lower wall (jmin or j=1)
       ! ==========================
       dudy_mean_jmin=0.0_wp
       if (is_bc_wall(2,1)) then
          do i=1,nx
             do k=1,nz
                dudy_mean_jmin = dudy_mean_jmin + rho(i,1,k)/visc(i,1,k)*utau_jmin(i,k)**2
             enddo
          enddo
       endif

       ! compute global (spatial) mean
       call MPI_ALLREDUCE(dudy_mean_jmin,dudy_meang_jmin,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)

       ! compute global dudy_jmin
       dudy_jmin = dudy_meang_jmin/ngx/ngz

       call MPI_BARRIER(COMM_global,info)


       ! Compute dudy on upper wall (jmax or j=ny)
       ! ==========================
       dudy_mean_jmax=0.0_wp
       if (is_bc_wall(2,2)) then
          do i=1,nx
             do k=1,nz
                dudy_mean_jmax = dudy_mean_jmax + rho(i,ny,k)/visc(i,ny,k)*utau_jmax(i,k)**2
             enddo
          enddo
       endif

       ! compute global (spatial) mean
       call MPI_ALLREDUCE(dudy_mean_jmax,dudy_meang_jmax,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)

       ! compute global dudy_jmin
       dudy_jmax = dudy_meang_jmax/ngx/ngz

       ! Average dudy on upper and lower walls
       ! =====================================
       dudy = (dudy_jmin+dudy_jmax)/2

    else
       ! Compute dudy on lower wall (jmin or j=1)
       ! ==========================
       ! compute local (spatial) mean for j=2 and 3
       Uav2_jmin=0.0_wp
       Uav3_jmin=0.0_wp
       if (coord(2)==0) then
          do i=1,nx
             do k=1,nz
                Uav2_jmin=Uav2_jmin+uu(i,2,k)
             enddo
          enddo
          do i=1,nx
             do k=1,nz
                Uav3_jmin=Uav3_jmin+uu(i,3,k)
             enddo
          enddo
       endif
       ! compute global (spatial) mean for j=2 and 3
       call MPI_ALLREDUCE(Uav2_jmin,Uavg2_jmin,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
       call MPI_ALLREDUCE(Uav3_jmin,Uavg3_jmin,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
       idyw=1./abs(yg(2)-yg(1))/ngx/ngz
       if (is_curv) idyw=1./abs(ygc(1,2)-ygc(1,1))/ngx/ngz
       ! compute global dudy_jmin
       dudy_jmin=(2.0_wp*Uavg2_jmin-0.5_wp*Uavg3_jmin)*idyw

       !print *,'forcing'
       !!print *,uu(:,2,10)
       !!print *,uu(:,3,10)
       !print *,Uavg2_jmin,Uavg3_jmin,idyw

       call MPI_BARRIER(COMM_global,info)

       ! Compute dudy on upper wall (jmax or j=ny)
       ! ==========================
       ! compute local (spatial) mean for j=ny-1 and ny-2
       Uav2_jmax=0.0_wp
       Uav3_jmax=0.0_wp
       if (coord(2)==ndomy-1) then
          do i=1,nx
             do k=1,nz
                Uav2_jmax=Uav2_jmax+uu(i,ny-1,k)
             enddo
          enddo
          do i=1,nx
             do k=1,nz
                Uav3_jmax=Uav3_jmax+uu(i,ny-2,k)
             enddo
          enddo
       endif
       ! compute global (spatial) mean for j=2 and 3
       call MPI_ALLREDUCE (Uav2_jmax,Uavg2_jmax,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
       call MPI_ALLREDUCE (Uav3_jmax,Uavg3_jmax,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
       idyw=1./abs(yg(ngy)-yg(ngy-1))/ngx/ngz
       if (is_curv) idyw=1./abs(ygc(1,ngy)-ygc(1,ngy-1))/ngx/ngz
       ! compute global dudy_jmax
       dudy_jmax=(-2.0_wp*Uavg2_jmax+0.5_wp*Uavg3_jmax)*idyw

       ! Average dudy on upper and lower walls
       ! =====================================
       dudy=dudy_jmin-dudy_jmax
       !if (iproc==0) print *,'dudy',dudy,dudy_jmin,dudy_jmax
    endif


    ! Update mass flow rate at new iteration
    ! ======================================
    gn=Ly*Lz*forc_rhou + Lz*mu_ref*dudy

    Q_np1=Q_n-deltat*gn

    ! Update forcing term (target is Q_i)
    ! ===================
    forc_rhou=forc_rhou+alpha_f*(Q_np1-Q_i)+beta_f*(Q_n-Q_i)

    ! Mass-flow-rate target
    ! =====================
    Q_in=rhou_bulk*Ly*Lz ! initial value Q_i to be corrected by actual mass flow rate
    
    ! Check at screen
    ! ===============
    if ((iproc==0).and.(mod(ntime,nprint)==0).and.(irk==nrk)) then
    !if ((iproc==0).and.(mod(ntime,100)==0).and.(irk==nrk)) then
       write(6,*) 'Q/Qi=',Q_n,Q_n/Q_i,forc_rhou
       write(6,*) 'Q/Qinit=',rhou_bulk*Ly*Lz,Q_n/(rhou_bulk*Ly*Lz),Q_n/Q_in
    endif

    ! Apply relaxation term for the new mass flow rate
    ! ================================================
    Q_n=Q_np1+0.001*(Q_in-Q_i)

  end subroutine forcing_rhou

  !===============================================================================
  subroutine forcing_rhou_c
  !===============================================================================
    !> Compute forcing term for periodic hill flow
  !===============================================================================
    use mod_flow0 ! TO BE CHANGED not T&D for periodic hill !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: Q_np1,gn,Q_in
    real(wp) :: Q_inp ! local value on proc
    ! compute friction force on the upper and lower walls
    real(wp) :: fx_jmin,fxg_jmin,fx_jmax,fxg_jmax
    ! ----------------------------------------------------------------------------

    fx_jmin=0.
    fx_jmax=0.

    ! Compute the contribution of stress tensor on walls
    ! ==================================================
    ! (only upper and lower wall for periodic hill flow)

    ! lower wall (jmin or j=1)
    ! ----------------------
    if (BC_face(2,1)%sort==0) then 
       do k=1,nz
          do i=1,nx
             !fx_jmin = fx_jmin-Frhou(i,1,k)*nxndl_jmin(i)-Frhov(i,1,k)*nyndl_jmin(i)
             fx_jmin = fx_jmin-Grhou(i,1,k)*dl_jmin(i)
          enddo
       enddo
    endif
    ! upper wall (jmax or j=ny)
    ! ----------------------
    if (BC_face(2,2)%sort==0) then
       do k=1,nz
          do i=1,nx
             !fx_jmax = fx_jmax-Frhou(i,ny,k)*nxndl_jmax(i)-Frhov(i,ny,k)*nyndl_jmax(i)
             fx_jmax = fx_jmax-Grhou(i,ny,k)*dl_jmax(i)
          enddo
       enddo
    endif
    ! summation for all processes
    ! ---------------------------
    call MPI_ALLREDUCE(fx_jmin,fxg_jmin,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
    call MPI_ALLREDUCE(fx_jmax,fxg_jmax,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)

    fxg_jmin=fxg_jmin*deltaz
    fxg_jmax=fxg_jmax*deltaz

    ! Update mass flow rate at new iteration
    ! ======================================
    gn=(Vol*forc_rhou+fxg_jmin-fxg_jmax)/Lx

    Q_np1=Q_n-deltat*gn

    ! Update forcing term
    ! ===================
    forc_rhou=forc_rhou+alpha_f*(Q_np1-Q_i)+beta_f*(Q_n-Q_i)

    ! Check at screen
    ! ===============
    if ((iproc==0).and.(mod(ntime,nprint)==0).and.(irk==nrk)) then
       write(6,*) 'Q/Qi=',Q_n,Q_n/Q_i,forc_rhou
    endif

    ! Mass-flow-rate target
    ! =====================
    !Q_in=rhou_bulk*Ly*Lz ! initial value to be corrected by actual mass flow rate

    ! Integrate bulk in inlet plane (Att. use of temporal mean u0 ????????????????????)
    ! ============================= 
    Q_inp=0.
    if ((coord(1)==0).and.(nob(iproc)==1)) then ! only first block for multiblock tests
       do j=1,ny
          do k=1,nz
             Q_inp=Q_inp+rho_ref*u0(1,j,k)*dyi(j)*deltaz
          enddo
       enddo
    endif
    ! actual mass flow rate
    call MPI_ALLREDUCE(Q_inp,Q_in,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)

    ! Apply relaxation term for the new mass flow rate
    ! ================================================
    Q_n=Q_np1+0.001*(Q_in-Q_i)

  end subroutine forcing_rhou_c

  !===============================================================================
  subroutine forcing_rho
  !===============================================================================
    !> Normalize bulk density
    !> (useful to avoid a drift of density in supersonic cases)
  !===============================================================================
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: ro
    ! ----------------------------------------------------------------------------

    ! Computation of <rho> bulk
    ! =========================
    ! volume integration
    ro=0.0_wp
    do k=1,nz
       do j=1,ny
          do i=1,nx
             ro=ro+rho_n(i,j,k)*dyi(j)*dzi(k)
          enddo
       enddo
    enddo
    ! summation for all processes
    call MPI_ALLREDUCE(MPI_IN_PLACE,ro,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
    ! actual rho_bulk
    rhob_SCH  = ro/(Ly*ngx*Lz)

    ! Weighting term to correct rho
    ! =============================
    forc_rho = rhob_TRG/rhob_SCH

    ! Update conservative variables with density correction
    ! =====================================================
    rho_n  = rho_n *forc_rho
    rhou_n = rhou_n*forc_rho
    rhov_n = rhov_n*forc_rho
    rhow_n = rhow_n*forc_rho
    rhoe_n = rhoe_n*forc_rho

  end subroutine forcing_rho

  !===============================================================================
  subroutine check_forcing
  !===============================================================================
    !> Check of the forcing term
  !===============================================================================
    use mod_eos
    implicit none
    ! ----------------------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: rhoavg_w,muavg_w,mududy_w,c_w,rou
    real(wp) :: utau_mean,utau_meang
    real(wp), dimension(3) :: moyloc,moytot
    ! ----------------------------------------------------------------------------
    
    ! Computation of averaged friction Reynolds number
    ! ================================================
    moyloc = 0.0_wp
    moytot = 0.0_wp
    
    ! Computation of wall quantities
    ! ==============================
    if (coord(2)==0) then
       j=1
       do k=1,nz
          do i=1,nx
             moyloc(1)= moyloc(1)+visc(i,j,k)*duy(i,j,k)
             moyloc(2)= moyloc(2)+rho_n(i,j,k)
             moyloc(3)= moyloc(3)+visc(i,j,k)
          enddo
       enddo
    endif
    if (coord(2)==ndomy-1) then
       j=ny
       do k=1,nz
          do i=1,nx
             moyloc(1)= moyloc(1)-visc(i,j,k)*duy(i,j,k)
             moyloc(2)= moyloc(2)+rho_n(i,j,k)
             moyloc(3)= moyloc(3)+visc(i,j,k)
          enddo
       enddo
    endif

    call MPI_ALLREDUCE(moyloc,moytot,3,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)

    mududy_w= 0.5_wp*moytot(1)/(ngx*ngz)
    rhoavg_w= 0.5_wp*moytot(2)/(ngx*ngz)
    muavg_w = 0.5_wp*moytot(3)/(ngx*ngz)

    ! Computation of <rhou>_bulk
    ! ==========================
    rou = 0.0_wp
    do k=1,nz
       do j=1,ny
          do i=1,nx
             rou=rou+rhou_n(i,j,k)*dyi(j)*dzi(k)
          enddo
       enddo
    enddo
    call MPI_ALLREDUCE(MPI_IN_PLACE,rou,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
    rhoub_SCH=rou/(Ly*ngx*Lz)
    rhou_bulk=rhoub_SCH

    if (is_wall_model) then
       utau_mean = 0.0_wp; utau_meang=0.0_wp
       if (is_bc_wall(2,1)) utau_mean = utau_mean + SUM(utau_jmin(1:nx,1:nz))
       if (is_bc_wall(2,2)) utau_mean = utau_mean + SUM(utau_jmax(1:nx,1:nz))
       call MPI_ALLREDUCE(utau_mean,utau_meang,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
       utau_meang = utau_meang/2/ngx/ngz
    endif

    if (iproc.eq.0) then

       c_w = sqrt(c2calc_tro(T_wall,rhoavg_w))

       if (is_wall_model) then
          utau = utau_meang
       else
          utau = sqrt(mududy_w/rhoavg_w)
       endif

       open(unit=98, file='convqp.dat', form='formatted', status='unknown' &
            , position='append' )

       write(98,'(I7,1X,9(g0,2X))')  ntotal, tstar & ! 1,2
            ,  utau*hc*rhoavg_w/muavg_w   & ! 3
            ,  rhob_SCH/rhob_TRG          & ! 4
            , rhoavg_w / rhob_TRG         & ! 5
            , rhoub_SCH/rhoub_TRG         & ! 6
            , Mach*c_w*rho_bulk/rhoub_TRG & ! 7
            , rhoub_SCH/rhob_TRG/c_w      & ! 8
            , -forc_rhou                  & ! 9
            , rhoub_SCH*hc/muavg_w          ! 10
       close(98)
    endif

    !if (isnan(sum(rho))) call mpistop('NaN in density!! Abort', 1)

  end subroutine check_forcing

end module mod_forcing_bulk
