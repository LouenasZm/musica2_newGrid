!=================================================================
module mod_routines_rans
!=================================================================
  !> Subroutines: residuals, source, 
!=================================================================
  use precision
  implicit none
  ! --------------------------------------------------------------
  real(wp) :: Res_nutil
  ! --------------------------------------------------------------

contains
 
  !===============================================================
  subroutine residuals_sa
  !===============================================================
    !> Compute residuals (for conservative variables)
  !===============================================================
    use mod_mpi
    use mod_flow
    use mod_constant
    use mod_time
    implicit none
    ! ------------------------------------------------------------
    integer :: i,j
    integer :: nb_pts
    real(wp), dimension(:,:), allocatable :: R_nutil_dim
    real(wp), dimension(:), allocatable :: maxR_nutil_dim
    integer , dimension(:), allocatable ::  maxR_i_dim,maxR_j_dim,maxR_ij,maxR_ind
    real(wp) :: maxR_nutil_g,R_nutil!,Res_nutil
    integer  :: maxR_i_g,maxR_j_g
    ! ------------------------------------------------------------

    !******************!
    ! Spalart-Allmaras !
    !******************!

    ! find max residual in computational domain and report its value, proc and index
    allocate(R_nutil_dim(ndx:nfx,ndy:nfy),maxR_nutil_dim(1:nproc))
    allocate(maxR_ij(1:2),maxR_i_dim(1:nproc),maxR_j_dim(1:nproc))
    allocate(maxR_ind(1))

    do i=ndx,nfx
       do j=ndy,nfy
          R_nutil_dim(i,j) = log10(max(abs(nutil(i,j,1)-nutil_n(i,j,1))/mu_ref*rho_ref,1.e-20_wp))
       enddo
    enddo

    maxR_ij=maxloc(R_nutil_dim)
    
    call MPI_GATHER(maxval(R_nutil_dim),1,MPI_DOUBLE_PRECISION,maxR_nutil_dim,1,MPI_DOUBLE_PRECISION,0,COMM_global,info)
    call MPI_GATHER(maxR_ij(1),1,MPI_INTEGER,maxR_i_dim,1,MPI_INTEGER,0,COMM_global,info)
    call MPI_GATHER(maxR_ij(2),1,MPI_INTEGER,maxR_j_dim,1,MPI_INTEGER,0,COMM_global,info)

    ! write max residual to screen
    if (iproc.eq.0) then
       maxR_nutil_g = maxval(maxR_nutil_dim)
       maxR_ind     = maxloc(maxR_nutil_dim)

       maxR_i_g = maxR_i_dim(maxR_ind(1))
       maxR_j_g = maxR_j_dim(maxR_ind(1))

       write(6,102) maxR_nutil_g,maxR_i_g,maxR_j_g,maxR_ind
    102 format(1x,'~> MAX residual log10(nutilde)=',f8.4,' i=',i3,' j=',i3,' proc=',i2)
    endif

    deallocate(R_nutil_dim,maxR_nutil_dim,maxR_ij,maxR_i_dim,maxR_j_dim,maxR_ind)

    R_nutil = 0.0_wp

    do i=ndx,nfx
       do j=ndy,nfy
          R_nutil = R_nutil +(( nutil(i,j,1)- nutil_n(i,j,1))/mu_ref*rho_ref)**2
       enddo
    enddo

    nb_pts=(nfx-ndx+1)*(nfy-ndy+1)

    R_nutil = sqrt(R_nutil)/dble(nb_pts)

    R_nutil = max(R_nutil ,1.e-20_wp)

    ! MPI communications
    ! ==================
    call MPI_REDUCE(R_nutil, Res_nutil, 1,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMM_global,info)

    ! print residuals
    ! ===============
    if (iproc==0) then
       Res_nutil =Res_nutil/dble(nproc)
    endif


  end subroutine residuals_sa

end module mod_routines_rans
