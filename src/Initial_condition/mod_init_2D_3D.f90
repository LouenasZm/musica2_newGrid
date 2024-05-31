!===============================================================================
module mod_init_2D_3D
!===============================================================================
  !> Module to restart 3D case from 2D one (simple extrusion)
  ! limited to periodic flow along z
!===============================================================================
  implicit none
  ! ----------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------

contains

  !=============================================================================
  subroutine extrude_2D_field
  !=============================================================================
    !> author: XG
    !> date: September 2021
    !> read 2D field and perform 2D->3D extrusion
  !=============================================================================
    use mod_ngh
    use mod_mpi
    use mod_block
    use mod_io
    use mod_utils
    use mod_flow0
    use mod_io_restartTD
    use mod_rans
    !use warnstop
    implicit none
    ! --------------------------------------------------------------------------
    integer :: k,n,m1,m2,dimTD
    integer :: i1,i2,j1,j2
    real(wp), dimension(:,:,:), allocatable :: rho_o,rhou_o,rhov_o,rhow_o,rhoe_o,nutil_o
    real(wp), dimension(:,:,:), allocatable :: rho0_o,u0_o,v0_o,w0_o,p0_o,T0_o
    ! --------------------------------------------------------------------------
    
    ! Allocate old field
    ! ==================
    allocate( rho_o(1-ngh:nx+ngh,1-ngh:ny+ngh,1))
    allocate(rhou_o(1-ngh:nx+ngh,1-ngh:ny+ngh,1))
    allocate(rhov_o(1-ngh:nx+ngh,1-ngh:ny+ngh,1))
    allocate(rhow_o(1-ngh:nx+ngh,1-ngh:ny+ngh,1))
    allocate(rhoe_o(1-ngh:nx+ngh,1-ngh:ny+ngh,1))
    if (is_RANS.and.model_RANS.eq.'SA') allocate(nutil_o(1-ngh_r:nx+ngh_r,1-ngh_r:ny+ngh_r,1))
    
    ! Create MPI_IO types and field structure for old field [field_o]
    ! =====================================================
    allocate(field_o(nbloc))
    field_o(nob(iproc))%MPI_COMM = COMM_intrablock
    call mod_io_init(ngx,ngy,1,nx,ny,1,nx,ny,1, &
         ngh,3,coord,is_IOtec_read,is_IOtec_write,field_o(nob(iproc)))

    ! Read old field
    ! ==============
    binfile='restart_bl'//trim(numchar(nob(iproc)))//'.bin'
    if (is_RANS.and.model_RANS.eq.'SA') then
       call read_vol_o_rans(binfile, rho_o(1:nx,1:ny,1:1) &
                                   ,rhou_o(1:nx,1:ny,1:1) &
                                   ,rhov_o(1:nx,1:ny,1:1) &
                                   ,rhow_o(1:nx,1:ny,1:1) &
                                   ,rhoe_o(1:nx,1:ny,1:1) &
                                   ,nutil_o(1:nx,1:ny,1:1))
    else
       call read_vol_o(binfile, rho_o(1:nx,1:ny,1:1) &
                              ,rhou_o(1:nx,1:ny,1:1) &
                              ,rhov_o(1:nx,1:ny,1:1) &
                              ,rhow_o(1:nx,1:ny,1:1) &
                              ,rhoe_o(1:nx,1:ny,1:1) )
    endif

    ! Extrude 2D field in 3D (assuming that the field is uniform in the spanwise direction)
    ! ======================
    do k=1,nz
        rho(1:nx,1:ny,k)= rho_o(1:nx,1:ny,1)
       rhou(1:nx,1:ny,k)=rhou_o(1:nx,1:ny,1)
       rhov(1:nx,1:ny,k)=rhov_o(1:nx,1:ny,1)
       rhow(1:nx,1:ny,k)=rhow_o(1:nx,1:ny,1)
       rhoe(1:nx,1:ny,k)=rhoe_o(1:nx,1:ny,1)
    enddo
    if (is_RANS.and.model_RANS.eq.'SA') then
       do k=1,nz
          nutil(1:nx,1:ny,k)=nutil_o(1:nx,1:ny,1)
       enddo
    endif

    ! Deallocate variable for old field
    ! =================================
    deallocate(rho_o,rhou_o,rhov_o,rhow_o,rhoe_o)
    if (is_RANS.and.model_RANS.eq.'SA') deallocate(nutil_o)

    ! Read and extrude old mean field
    ! ===============================
    if (is_mean0) then
       allocate(rho0_o(1-ngh:nx+ngh,1-ngh:ny+ngh,1))
       allocate(  u0_o(1-ngh:nx+ngh,1-ngh:ny+ngh,1))
       allocate(  v0_o(1-ngh:nx+ngh,1-ngh:ny+ngh,1))
       allocate(  w0_o(1-ngh:nx+ngh,1-ngh:ny+ngh,1))
       allocate(  p0_o(1-ngh:nx+ngh,1-ngh:ny+ngh,1))
       allocate(  T0_o(1-ngh:nx+ngh,1-ngh:ny+ngh,1))

       ! Read old mean field
       binfile='mean0_bl'//trim(numchar(nob(iproc)))//'.bin'
       call read_mean_o(binfile,rho0_o(1:nx,1:ny,1:1) &
                                 ,u0_o(1:nx,1:ny,1:1) &
                                 ,v0_o(1:nx,1:ny,1:1) &
                                 ,w0_o(1:nx,1:ny,1:1) &
                                 ,p0_o(1:nx,1:ny,1:1) &
                                 ,T0_o(1:nx,1:ny,1:1) )
       
       ! Extrude mean 2D field in 3D
       do k=1,nz
          rho0(1:nx,1:ny,k)=rho0_o(1:nx,1:ny,1)
            u0(1:nx,1:ny,k) = u0_o(1:nx,1:ny,1)
            v0(1:nx,1:ny,k) = v0_o(1:nx,1:ny,1)
            w0(1:nx,1:ny,k) = w0_o(1:nx,1:ny,1)
            p0(1:nx,1:ny,k) = p0_o(1:nx,1:ny,1)
            T0(1:nx,1:ny,k) = T0_o(1:nx,1:ny,1)
       enddo
       
       ! Deallocate mean variable for old field
       deallocate(rho0_o,u0_o,v0_o,w0_o,p0_o,T0_o)
       
     endif
 
    ! Deallocate I/O strcture for old field
    ! =====================================
    deallocate(field_o)

    if (iproc==0) print *,'2D->3D extrusion of restart files OK'
    
    ! 2D-3D extrusion of restartTD for Tam & Dong's BCs
    ! =================================================
    if (is_TamDong) then

       ! Initialize restartTD type for 2D field
       ! --------------------------------------
       call init_io_restartTD(restartTD_o,1)

       ! Read 2-D restartTD.bin file
       ! ---------------------------
       call read_restartTD('restartTD.bin',restartTD_o)
       
       ! Extrude time-averaged primitive variables
       ! -----------------------------------------
       dimTD=2
       if (is_TamDong3D) dimTD=3

       do m1=1,dimTD
          do m2=1,2
             ! Extrude each variable if BC is T&D
             if (is_bc_TD(m1,m2)) then
                i1=restartTD(m1,m2)%i1
                i2=restartTD(m1,m2)%i2
                j1=restartTD(m1,m2)%j1
                j2=restartTD(m1,m2)%j2
                do n=1,nvar
                   do k=1,nz
                      BC_face(m1,m2)%U0(i1:i2,j1:j2,k,n)=BC_face(m1,m2)%U0(i1:i2,j1:j2,1,n)
                   enddo
                enddo
             endif
          enddo
       enddo

       ! Free types for MPI-IO of 2D restartTD
       ! -------------------------------------
       call free_restartTD(restartTD_o)
       
       if (iproc==0) print *,'2D->3D extrusion of restartTD file OK'
    endif

  end subroutine extrude_2D_field
  
end module mod_init_2D_3D
