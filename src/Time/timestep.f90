!===============================================================================
subroutine timestep
!===============================================================================
  !> computation of time step  TO BE CHANGED
!===============================================================================
   use mod_mpi
   use mod_time
   use mod_flow     ! <- for uu,vv,ww
   use mod_constant ! <- for ref values
   use mod_tranprop ! <- for visc TO BE CHANGED visc is in mod_flow
   use mod_eos      ! <- for c2calc_tro
   implicit none
   ! ---------------------------------------------------------------------------
   integer  :: i,j,k
   real(wp) :: deltamin,tf,deltatproc,csound,testc,testd
   real(wp) :: delta_min,deltm
   ! ---------------------------------------------------------------------------

   if (is_dtvar) then

      deltatproc = 1.e10_wp

      do k=1,nz
         do j=1,ny
            do i=1,nx
               deltamin = 1.0_wp/max(idx(i),idy(j),idz(k))

               csound=sqrt(c2calc_tro(Tmp(i,j,k),rho(i,j,k)))

               ! Convective condition
               ! testc = CFL/( ( abs(uu(i,j,k)) + csound )*idx(i) &
               !             + ( abs(vv(i,j,k)) + csound )*idy(j) &
               !             + ( abs(ww(i,j,k)) + csound )*idz(k) )

               testc = CFL/max( ( abs(uu(i,j,k)) + csound )*idx(i) &
                              , ( abs(vv(i,j,k)) + csound )*idy(j) &
                              , ( abs(ww(i,j,k)) + csound )*idz(k) )

               ! Viscous condition
               testd = CFL*deltamin**2*rho(i,j,k)/visc(i,j,k)

               deltatproc = min(deltatproc, testc, testd)

            enddo
         enddo
      enddo

      call MPI_ALLREDUCE(deltatproc,deltat,1,MPI_DOUBLE_PRECISION,MPI_MIN,COMM_global,info)

   else
      
      ! First mesh size (to compute timestep)
      ! ===============
      
      if (PHILL) then
         ! Definition de delta min
         deltm=1.

         do i=1,ngx-1
            do j=1,ngy
               deltm=min(deltm,sqrt((xgc(i+1,j)-xgc(i,j))**2+(ygc(i+1,j)-ygc(i,j))**2))
            end do
         end do

         do i=1,ngx
            do j=1,ngy-1
               deltm=min(deltm,sqrt((xgc(i,j+1)-xgc(i,j))**2+(ygc(i,j+1)-ygc(i,j))**2))
            end do
         end do

         do i=1,ngx-1
            do j=1,ngy-1
               deltm=min(deltm,sqrt((xgc(i+1,j+1)-xgc(i,j))**2+(ygc(i+1,j+1)-ygc(i,j))**2))
            end do
         end do

         ! Search smallest delta
         delta_min=1.
         delta_min=min(delta_min,deltm)
      else
         delta_min=deltay
      endif

      deltatproc = CFL*delta_min/(u_ref+c_ref)

      call MPI_ALLREDUCE(deltatproc,deltat,1,MPI_DOUBLE_PRECISION,MPI_MIN,COMM_global,info)
         
      ! if (CHIT) then
      !    !deltat = 0.15_wp*delta_min/(u_ref+c_ref)!u_ref
      !    deltat = 1.0_wp*delta_min/(u_ref+c_ref)!u_ref
      ! endif
      
      !deltat = 1.096115504379462E-008
   endif
!!$deltat =2.293205278826566E-008
   ! Select between viscous/inviscid time-step condition
   dtstar = deltat/tscale

   tf = timemax*tscale
   !print *,deltat,tf,timemax,tscale
   !print *,nmax,tf/deltat

   nmax = min(nmax, int(tf/deltat)+1)

   if (mod(ntime,nprint).eq.0) then
      if (iproc.eq.0) then
         write(*,*) 'Dt [s], Dt* [-]:',deltat,dtstar
      endif
   endif

end subroutine timestep
