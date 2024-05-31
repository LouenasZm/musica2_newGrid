!===============================================================================
subroutine update_var
!===============================================================================
  !> Update conservative variables at next RK step + forcing/source terms
!===============================================================================
  use mod_flow         ! for: conservative variables
  use mod_time         ! for: irk,crk,deltat
  use mod_forcing_bulk ! for: CHAN,PHILL cases
!!$  use mod_mpi_part     ! for: averaging in adjoint block arrangement
  implicit none
  ! ----------------------------------------------------------------------------
  integer  :: i,j,k
  real(wp) :: alpha
!!$  integer :: ifrst,ilast,idrct,jfrst,jlast,jdrct
  ! ----------------------------------------------------------------------------

  ! Update variables
  ! ================
  if (is_dtlocal) then
     do k=1,nz
        do j=1,ny
           do i=1,nx
              rho_n(i,j,k)  = rho(i,j,k)  - dt_local(i,j,k)*crk(irk)*Krho(i,j,k)
              rhou_n(i,j,k) = rhou(i,j,k) - dt_local(i,j,k)*crk(irk)*Krhou(i,j,k)
              rhov_n(i,j,k) = rhov(i,j,k) - dt_local(i,j,k)*crk(irk)*Krhov(i,j,k)
              rhow_n(i,j,k) = rhow(i,j,k) - dt_local(i,j,k)*crk(irk)*Krhow(i,j,k)
              rhoe_n(i,j,k) = rhoe(i,j,k) - dt_local(i,j,k)*crk(irk)*Krhoe(i,j,k)
           enddo
        enddo
     enddo
  else
     alpha=deltat*crk(irk)

     do k=1,nz
        do j=1,ny
           do i=1,nx
              rho_n(i,j,k)  = rho(i,j,k)  - alpha*Krho(i,j,k)
              rhou_n(i,j,k) = rhou(i,j,k) - alpha*Krhou(i,j,k)
              rhov_n(i,j,k) = rhov(i,j,k) - alpha*Krhov(i,j,k)
              rhow_n(i,j,k) = rhow(i,j,k) - alpha*Krhow(i,j,k)
              rhoe_n(i,j,k) = rhoe(i,j,k) - alpha*Krhoe(i,j,k)
           enddo
        enddo
     enddo
  endif

!!$  ! averaging the flow values
!!$  if (is_curv.and.is_adjoint_block) then
!!$
!!$     !print *,'bloc',nob(iproc),'no_dble_crnr',no_dble_crnr
!!$       
!!$     ! corner averaging
!!$     do k=1,no_dble_crnr
!!$        if (dble_crnr(k)%dc_exist) then
!!$           i = dble_crnr(k)%dc_index(1)
!!$           j = dble_crnr(k)%dc_index(2)
!!$           !do i=1,nx
!!$           call MPI_ALLREDUCE(MPI_IN_PLACE, rho_n(i,j,1),1,MPI_DOUBLE_PRECISION,MPI_SUM,dc_COMM(k,6),info)
!!$           call MPI_ALLREDUCE(MPI_IN_PLACE,rhou_n(i,j,1),1,MPI_DOUBLE_PRECISION,MPI_SUM,dc_COMM(k,6),info)
!!$           call MPI_ALLREDUCE(MPI_IN_PLACE,rhov_n(i,j,1),1,MPI_DOUBLE_PRECISION,MPI_SUM,dc_COMM(k,6),info)
!!$           call MPI_ALLREDUCE(MPI_IN_PLACE,rhow_n(i,j,1),1,MPI_DOUBLE_PRECISION,MPI_SUM,dc_COMM(k,6),info)
!!$           call MPI_ALLREDUCE(MPI_IN_PLACE,rhoe_n(i,j,1),1,MPI_DOUBLE_PRECISION,MPI_SUM,dc_COMM(k,6),info)
!!$           rho_n(i,j,1) =  rho_n(i,j,1)/2.0_wp
!!$           rhou_n(i,j,1) = rhou_n(i,j,1)/2.0_wp
!!$           rhov_n(i,j,1) = rhov_n(i,j,1)/2.0_wp
!!$           rhow_n(i,j,1) = rhow_n(i,j,1)/2.0_wp
!!$           rhoe_n(i,j,1) = rhoe_n(i,j,1)/2.0_wp
!!$           !enddo
!!$        endif
!!$     enddo
!!$
!!$  endif
  
!!$  ! averaging the flow values at the 5-point double corners (in case of adjoint block case)
!!$  !if (is_curv.and.is_adjoint_block.and..false.) then
!!$  if (is_curv.and.is_adjoint_block) then

!!$     ! face averaging
!!$     do i=1,no_dble_crnr
!!$        do j=1,5
!!$           if(dble_crnr(i)%face_exist(j)) then
!!$              ifrst = dble_crnr(i)%face_index(j,1)
!!$              ilast = dble_crnr(i)%face_index(j,2)
!!$              idrct = dble_crnr(i)%face_index(j,3)
!!$              jfrst = dble_crnr(i)%face_index(j,4)
!!$              jlast = dble_crnr(i)%face_index(j,5)
!!$              jdrct = dble_crnr(i)%face_index(j,6)
!!$              call MPI_ALLREDUCE(MPI_IN_PLACE, rho_n(ifrst:ilast:idrct,jfrst:jlast:jdrct,1),5,MPI_DOUBLE_PRECISION,MPI_SUM,dc_COMM(i,j),info)
!!$              call MPI_ALLREDUCE(MPI_IN_PLACE,rhou_n(ifrst:ilast:idrct,jfrst:jlast:jdrct,1),5,MPI_DOUBLE_PRECISION,MPI_SUM,dc_COMM(i,j),info)
!!$              call MPI_ALLREDUCE(MPI_IN_PLACE,rhov_n(ifrst:ilast:idrct,jfrst:jlast:jdrct,1),5,MPI_DOUBLE_PRECISION,MPI_SUM,dc_COMM(i,j),info)
!!$              call MPI_ALLREDUCE(MPI_IN_PLACE,rhow_n(ifrst:ilast:idrct,jfrst:jlast:jdrct,1),5,MPI_DOUBLE_PRECISION,MPI_SUM,dc_COMM(i,j),info)
!!$              call MPI_ALLREDUCE(MPI_IN_PLACE,rhoe_n(ifrst:ilast:idrct,jfrst:jlast:jdrct,1),5,MPI_DOUBLE_PRECISION,MPI_SUM,dc_COMM(i,j),info)
!!$              rho_n(ifrst:ilast:idrct,jfrst:jlast:jdrct,1) =  rho_n(ifrst:ilast:idrct,jfrst:jlast:jdrct,1)/2.0_wp
!!$              rhou_n(ifrst:ilast:idrct,jfrst:jlast:jdrct,1) = rhou_n(ifrst:ilast:idrct,jfrst:jlast:jdrct,1)/2.0_wp
!!$              rhov_n(ifrst:ilast:idrct,jfrst:jlast:jdrct,1) = rhov_n(ifrst:ilast:idrct,jfrst:jlast:jdrct,1)/2.0_wp
!!$              rhow_n(ifrst:ilast:idrct,jfrst:jlast:jdrct,1) = rhow_n(ifrst:ilast:idrct,jfrst:jlast:jdrct,1)/2.0_wp
!!$              rhoe_n(ifrst:ilast:idrct,jfrst:jlast:jdrct,1) = rhoe_n(ifrst:ilast:idrct,jfrst:jlast:jdrct,1)/2.0_wp
!!$           endif
!!$        enddo
!!$     enddo
!!$
!!$     ! corner averaging
!!$     do k=1,no_dble_crnr
!!$        if (dble_crnr(k)%dc_exist) then
!!$           i = dble_crnr(k)%dc_index(1)
!!$           j = dble_crnr(k)%dc_index(2)
!!$           call MPI_ALLREDUCE(MPI_IN_PLACE, rho_n(i,j,1),1,MPI_DOUBLE_PRECISION,MPI_SUM,dc_COMM(k,6),info)
!!$           call MPI_ALLREDUCE(MPI_IN_PLACE,rhou_n(i,j,1),1,MPI_DOUBLE_PRECISION,MPI_SUM,dc_COMM(k,6),info)
!!$           call MPI_ALLREDUCE(MPI_IN_PLACE,rhov_n(i,j,1),1,MPI_DOUBLE_PRECISION,MPI_SUM,dc_COMM(k,6),info)
!!$           call MPI_ALLREDUCE(MPI_IN_PLACE,rhow_n(i,j,1),1,MPI_DOUBLE_PRECISION,MPI_SUM,dc_COMM(k,6),info)
!!$           call MPI_ALLREDUCE(MPI_IN_PLACE,rhoe_n(i,j,1),1,MPI_DOUBLE_PRECISION,MPI_SUM,dc_COMM(k,6),info)
!!$           rho_n(i,j,1) =  rho_n(i,j,1)/5.0_wp
!!$           rhou_n(i,j,1) = rhou_n(i,j,1)/5.0_wp
!!$           rhov_n(i,j,1) = rhov_n(i,j,1)/5.0_wp
!!$           rhow_n(i,j,1) = rhow_n(i,j,1)/5.0_wp
!!$           rhoe_n(i,j,1) = rhoe_n(i,j,1)/5.0_wp
!!$        endif
!!$     enddo
!!$
!!$  endif

  if (is_2D) rhow_n=0.0_wp

  ! Enforce mass-flow rate
  ! ======================
  if (is_forcing_bulk) then
     ! compute forcing term
     if (irk==nrk) then
        if (PHILL) then
           call forcing_rhou_c
        else
           call forcing_rhou
        endif
     endif
     ! update variables
     alpha=-alpha*forc_rhou !!!!! a modifier pour pas de temps local ????????????????????????????????
     do k=1,nz
        do j=1,ny
           do i=1,nx
              rhou_n(i,j,k) = rhou_n(i,j,k) + alpha
              rhoe_n(i,j,k) = rhoe_n(i,j,k) + alpha*uu(i,j,k)
           enddo
        enddo
     enddo
  endif
  !!call forcing_rho

end subroutine update_var

!===============================================================================
subroutine update_var2N
!===============================================================================
  !> Update conservative variables at next RK step + forcing/source terms
  !> 2N implementation of RK46 of Berland et al.
!===============================================================================
  use mod_flow         ! for: conservative variables
  use mod_time         ! for: irk,crk,deltat
  use mod_forcing_bulk ! for: CHAN,PHILL cases
  implicit none
  ! ----------------------------------------------------------------------------
  integer  :: i,j,k
  real(wp) :: alpha
  ! ----------------------------------------------------------------------------

  ! Update variables
  ! ================
  alpha=deltat*crk(irk)

  if (iproc==0.and.ntotal==1) print*,irk,Ark(irk),Brk(irk)
  
  do k=1,nz
     do j=1,ny
        do i=1,nx
           rhon(i,j,k) = Ark(irk)*rhon(i,j,k) -deltat*Krho(i,j,k)
           rhoun(i,j,k)= Ark(irk)*rhoun(i,j,k)-deltat*Krhou(i,j,k)
           rhovn(i,j,k)= Ark(irk)*rhovn(i,j,k)-deltat*Krhov(i,j,k)
           rhown(i,j,k)= Ark(irk)*rhown(i,j,k)-deltat*Krhow(i,j,k)
           rhoen(i,j,k)= Ark(irk)*rhoen(i,j,k)-deltat*Krhoe(i,j,k)

           rho_n(i,j,k) = rho_n(i,j,k)  + Brk(irk)*rhon(i,j,k)
           rhou_n(i,j,k)= rhou_n(i,j,k) + Brk(irk)*rhoun(i,j,k)
           rhov_n(i,j,k)= rhov_n(i,j,k) + Brk(irk)*rhovn(i,j,k)
           rhow_n(i,j,k)= rhow_n(i,j,k) + Brk(irk)*rhown(i,j,k)
           rhoe_n(i,j,k)= rhoe_n(i,j,k) + Brk(irk)*rhoen(i,j,k)
        enddo
     enddo
  enddo
               
  if (is_2D) rhow_n=0.0_wp

  ! Enforce mass-flow rate
  ! ======================
  if (is_forcing_bulk) then
     ! compute forcing term
     if (irk==nrk) then
        if (PHILL) then
           call forcing_rhou_c
        else
           call forcing_rhou
        endif
     endif
     ! update variables
     alpha=-alpha*forc_rhou
     do k=1,nz
        do j=1,ny
           do i=1,nx
              rhou_n(i,j,k) = rhou_n(i,j,k) + alpha
              rhoe_n(i,j,k) = rhoe_n(i,j,k) + alpha*uu(i,j,k)
           enddo
        enddo
     enddo
  endif

end subroutine update_var2N

!===============================================================================
subroutine update_1
!===============================================================================
  !> Update conservative variables at all but last RK step (OBSOLETE)
!===============================================================================
  use mod_flow         ! for: conservative variables
  use mod_time         ! for: irk,crk,deltat
  implicit none
  ! ----------------------------------------------------------------------------
  integer  :: i,j,k
  real(wp) :: alpha
  ! ----------------------------------------------------------------------------

  ! Update variables
  ! ================
  alpha=deltat*crk(irk)

  do k=1,nz
     do j=1,ny
        do i=1,nx
           rho_n(i,j,k)  = rho(i,j,k)  - alpha*Krho(i,j,k)
           rhou_n(i,j,k) = rhou(i,j,k) - alpha*Krhou(i,j,k)
           rhov_n(i,j,k) = rhov(i,j,k) - alpha*Krhov(i,j,k)
           rhow_n(i,j,k) = rhow(i,j,k) - alpha*Krhow(i,j,k)
           rhoe_n(i,j,k) = rhoe(i,j,k) - alpha*Krhoe(i,j,k)
        enddo
     enddo
  enddo

  ! secure initialization of increments TO BE CHANGED ????
  Krho  = 0.0_wp
  Krhou = 0.0_wp
  Krhov = 0.0_wp
  Krhow = 0.0_wp
  Krhoe = 0.0_wp
  
end subroutine update_1

!===============================================================================
subroutine update_2
!===============================================================================
  !> Update conservative variables at last RK step (OBSOLETE)
!===============================================================================
  use mod_flow         ! for: conservative variables
  use mod_time         ! for: irk,crk,deltat
  use mod_forcing_bulk ! for: CHAN,PHILL cases
  implicit none
  ! ----------------------------------------------------------------------------
  integer  :: i,j,k
  real(wp) :: alpha
  ! ----------------------------------------------------------------------------

  ! Update variables
  ! ================
  alpha=deltat*crk(irk)

  do k=1,nz
     do j=1,ny
        do i=1,nx
           rho_n(i,j,k)  = rho_n(i,j,k)  - alpha*Krho(i,j,k)
           rhou_n(i,j,k) = rhou_n(i,j,k) - alpha*Krhou(i,j,k)
           rhov_n(i,j,k) = rhov_n(i,j,k) - alpha*Krhov(i,j,k)
           rhow_n(i,j,k) = rhow_n(i,j,k) - alpha*Krhow(i,j,k)
           rhoe_n(i,j,k) = rhoe_n(i,j,k) - alpha*Krhoe(i,j,k)
        enddo
     enddo
  enddo

  if (is_2D) rhow_n=0.0_wp

  ! Enforce mass-flow rate
  ! ======================
  if (is_forcing_bulk) then
     ! compute forcing term
     if (is_curv) then
        call forcing_rhou_c
     else
        call forcing_rhou
     endif
     ! update variables
     alpha=-alpha*forc_rhou
     do k=1,nz
        do j=1,ny
           do i=1,nx
              rhou_n(i,j,k) = rhou_n(i,j,k) + alpha
              rhoe_n(i,j,k) = rhoe_n(i,j,k) + alpha*uu(i,j,k)
           enddo
        enddo
     enddo
  endif

end subroutine update_2

!===============================================================================
subroutine update_var_in_dw
!===============================================================================
  !> Update increments to have dwi before application of IRS operator
  !>  + addition of forcing/source terms
!===============================================================================
  use mod_flow         ! for: conservative variables
  use mod_time         ! for: irk,crk,deltat
  use mod_forcing_bulk ! for: CHAN,PHILL cases
  implicit none
  ! ----------------------------------------------------------------------------
  integer  :: i,j,k
  real(wp) :: alpha
  ! ----------------------------------------------------------------------------

  ! Enforce mass-flow rate
  ! ======================
  if (is_forcing_bulk) then
     ! compute forcing term
     if (irk==nrk) then
        if (PHILL) then
           call forcing_rhou_c
        else
           call forcing_rhou
        endif
     endif
     ! update increment
     alpha=deltat*crk(irk)*forc_rhou !!!!! a modifier pour pas de temps local ????????????????????????????????
     do k=1,nz
        do j=1,ny
           do i=1,nx
              Krhou(i,j,k) = Krhou(i,j,k) + alpha
              Krhoe(i,j,k) = Krhoe(i,j,k) + alpha*uu(i,j,k)
           enddo
        enddo
     enddo
  endif

end subroutine update_var_in_dw

!===============================================================================
subroutine update_var_of_dw
!===============================================================================
  !> Update conservative variables after application of IRS operator
!===============================================================================
  use mod_flow         ! for: conservative variables
  use mod_time         ! for: irk,crk,deltat
  implicit none
  ! ----------------------------------------------------------------------------
  integer  :: i,j,k
  real(wp) :: alpha
  ! ----------------------------------------------------------------------------

  ! Update variables
  ! ================
  if (is_dtlocal) then
     do k=1,nz
        do j=1,ny
           do i=1,nx
              rho_n(i,j,k)  = rho(i,j,k)  -dt_local(i,j,k)*crk(irk)*Krho(i,j,k)
              rhou_n(i,j,k) = rhou(i,j,k) -dt_local(i,j,k)*crk(irk)*Krhou(i,j,k)
              rhov_n(i,j,k) = rhov(i,j,k) -dt_local(i,j,k)*crk(irk)*Krhov(i,j,k)
              rhow_n(i,j,k) = rhow(i,j,k) -dt_local(i,j,k)*crk(irk)*Krhow(i,j,k)
              rhoe_n(i,j,k) = rhoe(i,j,k) -dt_local(i,j,k)*crk(irk)*Krhoe(i,j,k)
           enddo
        enddo
     enddo
  else
     alpha=-deltat*crk(irk)

     do k=1,nz
        do j=1,ny
           do i=1,nx
              rho_n(i,j,k)  = rho(i,j,k)  + alpha*Krho(i,j,k)
              rhou_n(i,j,k) = rhou(i,j,k) + alpha*Krhou(i,j,k)
              rhov_n(i,j,k) = rhov(i,j,k) + alpha*Krhov(i,j,k)
              rhow_n(i,j,k) = rhow(i,j,k) + alpha*Krhow(i,j,k)
              rhoe_n(i,j,k) = rhoe(i,j,k) + alpha*Krhoe(i,j,k)
           enddo
        enddo
     enddo
  endif

  if (is_2D) rhow_n=0.0_wp

end subroutine update_var_of_dw

!===============================================================================
subroutine update_var_of_dw2N
!===============================================================================
  !> Update conservative variables after application of IRS operator
!===============================================================================
  use mod_flow         ! for: conservative variables
  use mod_time         ! for: irk,crk,deltat
  use mod_mpi
  implicit none
  ! ----------------------------------------------------------------------------
  integer  :: i,j,k
  ! ----------------------------------------------------------------------------

  if (iproc==0.and.ntotal==1) print*,irk,Ark(irk),Brk(irk)

  ! Update variables
  ! ================
  do k=1,nz
     do j=1,ny
        do i=1,nx
           rhon(i,j,k) = Ark(irk)*rhon(i,j,k) -deltat*Krho(i,j,k)
           rhoun(i,j,k)= Ark(irk)*rhoun(i,j,k)-deltat*Krhou(i,j,k)
           rhovn(i,j,k)= Ark(irk)*rhovn(i,j,k)-deltat*Krhov(i,j,k)
           rhown(i,j,k)= Ark(irk)*rhown(i,j,k)-deltat*Krhow(i,j,k)
           rhoen(i,j,k)= Ark(irk)*rhoen(i,j,k)-deltat*Krhoe(i,j,k)

           rho_n(i,j,k) = rho_n(i,j,k)  + Brk(irk)*rhon(i,j,k)
           rhou_n(i,j,k)= rhou_n(i,j,k) + Brk(irk)*rhoun(i,j,k)
           rhov_n(i,j,k)= rhov_n(i,j,k) + Brk(irk)*rhovn(i,j,k)
           rhow_n(i,j,k)= rhow_n(i,j,k) + Brk(irk)*rhown(i,j,k)
           rhoe_n(i,j,k)= rhoe_n(i,j,k) + Brk(irk)*rhoen(i,j,k)
        enddo
     enddo
  enddo

  if (is_2D) rhow_n=0.0_wp

end subroutine update_var_of_dw2N

!===============================================================================
subroutine update_var_in_dw_
!===============================================================================
  !> Update increments to have dwi before application of IRS operator
  !>  + addition of forcing/source terms
!===============================================================================
  use mod_flow         ! for: conservative variables
  use mod_time         ! for: irk,crk,deltat
  use mod_forcing_bulk ! for: CHAN,PHILL cases
  implicit none
  ! ----------------------------------------------------------------------------
  integer  :: i,j,k
  real(wp) :: alpha
  ! ----------------------------------------------------------------------------

  ! Update variables
  ! ================
  if (is_dtlocal) then
     do k=1,nz
        do j=1,ny
           do i=1,nx
              Krho(i,j,k)  = -dt_local(i,j,k)*crk(irk)*Krho(i,j,k)
              Krhou(i,j,k) = -dt_local(i,j,k)*crk(irk)*Krhou(i,j,k)
              Krhov(i,j,k) = -dt_local(i,j,k)*crk(irk)*Krhov(i,j,k)
              Krhow(i,j,k) = -dt_local(i,j,k)*crk(irk)*Krhow(i,j,k)
              Krhoe(i,j,k) = -dt_local(i,j,k)*crk(irk)*Krhoe(i,j,k)
           enddo
        enddo
     enddo
  else
     alpha=-deltat*crk(irk)

     do k=1,nz
        do j=1,ny
           do i=1,nx
              Krho(i,j,k)  = alpha*Krho(i,j,k)
              Krhou(i,j,k) = alpha*Krhou(i,j,k)
              Krhov(i,j,k) = alpha*Krhov(i,j,k)
              Krhow(i,j,k) = alpha*Krhow(i,j,k)
              Krhoe(i,j,k) = alpha*Krhoe(i,j,k)
           enddo
        enddo
     enddo
  endif

  if (is_2D) Krhow=0.0_wp

  ! Enforce mass-flow rate
  ! ======================
  if (is_forcing_bulk) then
     ! compute forcing term
     if (irk==nrk) then
        if (PHILL) then
           call forcing_rhou_c
        else
           call forcing_rhou
        endif
     endif
     ! update increment
     alpha=-alpha*forc_rhou !!!!! a modifier pour pas de temps local ????????????????????????????????
     do k=1,nz
        do j=1,ny
           do i=1,nx
              Krhou(i,j,k) = Krhou(i,j,k) + alpha
              Krhoe(i,j,k) = Krhoe(i,j,k) + alpha*uu(i,j,k)
           enddo
        enddo
     enddo
  endif

end subroutine update_var_in_dw_

!===============================================================================
subroutine update_var_of_dw_
!===============================================================================
  !> Update conservative variables after application of IRS operator
!===============================================================================
  use mod_flow         ! for: conservative variables
  implicit none
  ! ----------------------------------------------------------------------------
  integer  :: i,j,k
  ! ----------------------------------------------------------------------------

  ! Update variables
  ! ================

  do k=1,nz
     do j=1,ny
        do i=1,nx
           rho_n(i,j,k)  = rho(i,j,k)  + Krho(i,j,k)
           rhou_n(i,j,k) = rhou(i,j,k) + Krhou(i,j,k)
           rhov_n(i,j,k) = rhov(i,j,k) + Krhov(i,j,k)
           rhow_n(i,j,k) = rhow(i,j,k) + Krhow(i,j,k)
           rhoe_n(i,j,k) = rhoe(i,j,k) + Krhoe(i,j,k)
        enddo
     enddo
  enddo

  if (is_2D) rhow_n=0.0_wp

end subroutine update_var_of_dw_

!===============================================================================
subroutine update_varn_filter
!===============================================================================
  !> Update conservative variables after application of numerical dissipation
  !> * version to be used with filtering *
!===============================================================================
  use mod_flow         ! for: conservative variables
  use mod_time         ! for: irk,crk,deltat
  implicit none
  ! ----------------------------------------------------------------------------
  integer  :: i,j,k
  ! ----------------------------------------------------------------------------

  ! Update variables
  ! ================

  do k=1,nz
     do j=1,ny
        do i=1,nx
           rho_n(i,j,k)  = rho_n(i,j,k)  + Krho(i,j,k)
           rhou_n(i,j,k) = rhou_n(i,j,k) + Krhou(i,j,k)
           rhov_n(i,j,k) = rhov_n(i,j,k) + Krhov(i,j,k)
           rhow_n(i,j,k) = rhow_n(i,j,k) + Krhow(i,j,k)
           rhoe_n(i,j,k) = rhoe_n(i,j,k) + Krhoe(i,j,k)
        enddo
     enddo
  enddo

end subroutine update_varn_filter

!===============================================================================
subroutine update_varn_artvisc
!===============================================================================
  !> Update conservative variables after application of numerical dissipation
  !> * version to be used with artificial viscosity *
!===============================================================================
  use mod_flow         ! for: conservative variables
  use mod_time         ! for: irk,crk,deltat
  implicit none
  ! ----------------------------------------------------------------------------
  integer  :: i,j,k
  ! ----------------------------------------------------------------------------

  ! Update variables
  ! ================
  if (is_dtlocal) then
     do k=1,nz
        do j=1,ny
           do i=1,nx
              rho_n(i,j,k)  = rho_n(i,j,k)  + dt_local(i,j,k)*Krho(i,j,k)
              rhou_n(i,j,k) = rhou_n(i,j,k) + dt_local(i,j,k)*Krhou(i,j,k)
              rhov_n(i,j,k) = rhov_n(i,j,k) + dt_local(i,j,k)*Krhov(i,j,k)
              rhow_n(i,j,k) = rhow_n(i,j,k) + dt_local(i,j,k)*Krhow(i,j,k)
              rhoe_n(i,j,k) = rhoe_n(i,j,k) + dt_local(i,j,k)*Krhoe(i,j,k)
           enddo
        enddo
     enddo
  else
     do k=1,nz
        do j=1,ny
           do i=1,nx
              rho_n(i,j,k)  = rho_n(i,j,k)  + deltat*Krho(i,j,k)
              rhou_n(i,j,k) = rhou_n(i,j,k) + deltat*Krhou(i,j,k)
              rhov_n(i,j,k) = rhov_n(i,j,k) + deltat*Krhov(i,j,k)
              rhow_n(i,j,k) = rhow_n(i,j,k) + deltat*Krhow(i,j,k)
              rhoe_n(i,j,k) = rhoe_n(i,j,k) + deltat*Krhoe(i,j,k)
           enddo
        enddo
     enddo
  endif

end subroutine update_varn_artvisc
