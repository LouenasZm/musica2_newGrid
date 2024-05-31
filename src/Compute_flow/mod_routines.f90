!=================================================================
module mod_routines
!=================================================================
  !> Subroutines: residuals, source, 
!=================================================================
  use precision
  implicit none
  ! --------------------------------------------------------------
  ! Source parameters
  logical :: is_src
  real(wp) :: ampl_src,omeg_src,b_src
  real(wp) :: x_src,y_src,z_src
  ! --------------------------------------------------------------

contains
  
  !===============================================================
  subroutine residuals
  !===============================================================
    !> Compute residuals (for conservative variables)
  !===============================================================
    use mod_mpi
    use mod_flow
    use mod_constant
    use mod_time
    use mod_routines_rans
    implicit none
    ! ------------------------------------------------------------
    integer :: i,j,k
    integer :: nb_pts
    real(wp) :: R_rho,R_rhou,R_rhov,R_rhoe!,R_nutil
    real(wp) :: Res_rho,Res_rhou,Res_rhov,Res_rhoe!,Res_nutil
    ! ------------------------------------------------------------

    ! init residuals
    ! ==============
    R_rho =0.0_wp
    R_rhou=0.0_wp
    R_rhov=0.0_wp
    R_rhoe=0.0_wp

    ! compute residuals (L2 norm)
    ! =================
    do k=ndz,nfz
       do j=ndy,nfy
          do i=ndx,nfx
             R_rho =R_rho +( rho(i,j,k)- rho_n(i,j,k))**2
             R_rhou=R_rhou+(rhou(i,j,k)-rhou_n(i,j,k))**2
             R_rhov=R_rhov+(rhov(i,j,k)-rhov_n(i,j,k))**2
             R_rhoe=R_rhoe+(rhoe(i,j,k)-rhoe_n(i,j,k))**2
          enddo
       enddo
    enddo

    nb_pts=(nfx-ndx+1)*(nfy-ndy+1)*(nfz-ndz+1)

    R_rho =sqrt(R_rho)/dble(nb_pts)
    R_rhou=sqrt(R_rhou)/dble(nb_pts)
    R_rhov=sqrt(R_rhov)/dble(nb_pts)
    R_rhoe=sqrt(R_rhoe)/dble(nb_pts)

    R_rho =max(R_rho ,1.e-20_wp)
    R_rhou=max(R_rhou,1.e-20_wp)
    R_rhov=max(R_rhov,1.e-20_wp)
    R_rhoe=max(R_rhoe,1.e-20_wp)

    ! MPI communications
    ! ==================
    call MPI_REDUCE(R_rho, Res_rho, 1,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMM_global,info)
    call MPI_REDUCE(R_rhou,Res_rhou,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMM_global,info)
    call MPI_REDUCE(R_rhov,Res_rhov,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMM_global,info)
    call MPI_REDUCE(R_rhoe,Res_rhoe,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMM_global,info)

    ! print residuals
    ! ===============
    if (iproc==0) then
       ! normalization
       Res_rho =Res_rho /dble(nproc)/rho_ref
       Res_rhou=Res_rhou/dble(nproc)/(rho_ref*u_ref)
       Res_rhov=Res_rhov/dble(nproc)/(rho_ref*u_ref)
       Res_rhoe=Res_rhoe/dble(nproc)/(rho_ref*e_ref)
       if (.not.is_RANS) then
         ! write at screen
          write(6,101) log10(Res_rho),log10(Res_rhou),log10(Res_rhov),log10(Res_rhoe)
       101 format(1x,'~> residuals: log10(drho)=',f8.4,' log10(drhou)=',f8.4,' log10(drhov)=',f8.4,&
       ' log10(drhoe)=',f8.4)
       endif

       ! write in file
       write(31) ntime
       write(31) log10(Res_rho)
       write(31) log10(Res_rhou)
       write(31) log10(Res_rhov)
       write(31) log10(Res_rhoe)

       if (is_RANS) then
          write(6,100) log10(Res_rho),log10(Res_rhou),log10(Res_rhov),log10(Res_rhoe),log10(Res_nutil)
          100 format(1x,'~> residuals: log10(drho)=',f8.4,' log10(drhou)=',f8.4,' log10(drhov)=',f8.4,&
          ' log10(drhoe)=',f8.4,' log10(dnutilde)=',f8.4)
          write(31) log10(Res_nutil)
       endif
    endif

  end subroutine residuals

  !===============================================================
  subroutine source
  !===============================================================
    !> Acoustic source
  !===============================================================
    use mod_mpi
    use mod_flow
    use mod_constant
    use mod_time
    use mod_eos
    implicit none
    ! ------------------------------------------------------------
    integer :: i,j,k
    real(wp) :: S,ar,trk,alpha,sign_x
    ! ------------------------------------------------------------

    ! parameters
    ar=-log(2.0_wp)/(b_src*deltax)**2
    z_src=0.0_wp

    ! time advancement
    alpha=crk(irk)*deltat
    trk=time+ck(irk)*deltat

    if (is_curv3) then
       do k=ndzt,nfzt
          do j=ndyt,nfyt
             do i=ndxt,nfxt
                S= omeg_src*ampl_src*sin(omeg_src*trk) &
                 *exp(ar*((xc3(i,j,k)-x_src)**2+(yc3(i,j,k)-y_src)**2+(zc3(i,j,k)-z_src)**2))
                prs(i,j,k)= prs(i,j,k)+ alpha*S
                rho_n(i,j,k)= rho_n(i,j,k)+ alpha*S/c_ref**2
             enddo
          enddo
       enddo
    else
       if (is_curv) then
          sign_x= 1.0_wp
          do k=ndzt,nfzt
             do j=ndyt,nfyt
                do i=ndxt,nfxt
                   S= omeg_src*ampl_src*sin(omeg_src*trk) &
                      *exp(ar*((xc(i,j)-sign_x*x_src-0.*u_ref*trk)**2+(yc(i,j)-y_src)**2+(z(k)-z_src)**2))
                   !if (((xc(i,j)-x_src)<1.e-3).and.((yc(i,j)-y_src)<1.e-3)) print *,iproc,S
                   prs(i,j,k)= prs(i,j,k)+ alpha*S
                   rho_n(i,j,k)= rho_n(i,j,k)+ alpha*S/c_ref**2
                enddo
             enddo
          enddo
       else
          do k=ndzt,nfzt
             do j=ndyt,nfyt
                do i=ndxt,nfxt
                   S= omeg_src*ampl_src*sin(omeg_src*trk) &
                      *exp(ar*((x(i)-x_src-0.*u_ref*trk)**2+(y(j)-y_src)**2+(z(k)-z_src)**2))
                   prs(i,j,k)= prs(i,j,k)+ alpha*S
                   rho_n(i,j,k)= rho_n(i,j,k)+ alpha*S/c_ref**2
                enddo
             enddo
          enddo
       endif
    endif

  end subroutine source

end module mod_routines

