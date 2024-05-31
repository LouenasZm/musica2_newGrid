!===============================================================================
subroutine compute_cl_cd3
!===============================================================================
  !> Compute lift and drag coefficients in 3D curvilinear
!===============================================================================
  use mod_mpi_part
  use mod_mpi_types_two_sided ! for indices iis,ijs,iks
  use mod_constant
  use mod_flow
  use mod_time
  implicit none
  ! ----------------------------------------------------------------------------
  integer :: i,j,k,i1,j1,k1
  real :: cl,cd,cz,clv,cdv,czv,fac
  ! ----------------------------------------------------------------------------

  ! Initialize lift and drag (for all procs)
  ! ========================
  cl=0.0_wp
  cd=0.0_wp
  cz=0.0_wp
  clv=0.0_wp
  cdv=0.0_wp
  czv=0.0_wp
  
  if (is_adjoint_block) then ! adjoint-block approach
                             ! ======================

     ! To avoid counting two times interface points, we can use iis, ijs & iks
     ! -----------------------------------------------------------------------
     ! Rq: indices declared and filled in mod_mpi_type_two_sided
     !     it takes the value 2 at block interfaces
     ! along i: iis(1) is the sending index for MPI comm with adjoint blocks
     ! along j: ijs(3) is the sending index for MPI comm with adjoint blocks
     ! along k: iks(5) is the sending index for MPI comm with adjoint blocks

     ! For walls at imin index
     ! -----------------------
     if (is_bc_wall(1,1)) then

        ! starting indices
        j1=ijs(3)
        k1=iks(5)
        
        ! contribution of pressure to the force
        do k=k1,nz
           do j=j1,ny
              cd = cd + prs(1,j,k)*nxnds_imin(j,k)
              cl = cl + prs(1,j,k)*nynds_imin(j,k)
              cz = cz + prs(1,j,k)*nznds_imin(j,k)
           enddo
        enddo

        ! contribution of viscous stresses to the force
        ! tau.nx=tau11*ksi_x+tau12*ksi_y+tau13*ksi_z -> Frhou in visc.
        ! tau.ny=tau12*ksi_x+tau22*ksi_y+tau23*ksi_z -> Frhov in visc.
        ! tau.nz=tau13*ksi_x+tau23*ksi_y+tau33*ksi_z -> Frhow in visc.
        do k=k1,nz
           do j=j1,ny
              cdv = cdv + Frhou(1,j,k)*sgn_imin
              clv = clv + Frhov(1,j,k)*sgn_imin
              czv = czv + Frhow(1,j,k)*sgn_imin
           enddo
        enddo


     endif

     ! For walls at imax index
     ! -----------------------
     if (is_bc_wall(1,2)) then

        ! starting indices
        j1=ijs(3)
        k1=iks(5)

        ! contribution of pressure to the force
        do k=k1,nz
           do j=j1,ny
              cd = cd + prs(nx,j,k)*nxnds_imax(j,k)
              cl = cl + prs(nx,j,k)*nynds_imax(j,k)
              cz = cz + prs(nx,j,k)*nznds_imax(j,k)
           enddo
        enddo

        ! contribution of viscous stresses to the force
        ! tau.nx=tau11*ksi_x+tau12*ksi_y+tau13*ksi_z -> Frhou in visc.
        ! tau.ny=tau12*ksi_x+tau22*ksi_y+tau23*ksi_z -> Frhov in visc.
        ! tau.nz=tau13*ksi_x+tau23*ksi_y+tau33*ksi_z -> Frhow in visc.
        do k=k1,nz
           do j=j1,ny
              cdv = cdv + Frhou(nx,j,k)*sgn_imax
              clv = clv + Frhov(nx,j,k)*sgn_imax
              czv = czv + Frhow(nx,j,k)*sgn_imax
           enddo
        enddo

     endif

     ! For walls at jmin index
     ! -----------------------
     if (is_bc_wall(2,1)) then

        ! starting indices
        i1=iis(1)
        k1=iks(5)

        ! contribution of pressure to the force
        do k=k1,nz
           do i=i1,nx
              cd = cd + prs(i,1,k)*nxnds_jmin(i,k)
              cl = cl + prs(i,1,k)*nynds_jmin(i,k)
              cz = cz + prs(i,1,k)*nznds_jmin(i,k)
           enddo
        enddo

        ! contribution of viscous stresses to the force
        ! tau.nx=tau11*eta_x+tau12*eta_y+tau13*eta_z -> Grhou in visc.
        ! tau.ny=tau12*eta_x+tau22*eta_y+tau23*eta_z -> Grhov in visc.
        ! tau.nz=tau13*eta_x+tau23*eta_y+tau33*eta_z -> Grhow in visc.
        do k=k1,nz
           do i=i1,nx
              cdv = cdv + Grhou(i,1,k)*sgn_jmin
              clv = clv + Grhov(i,1,k)*sgn_jmin
              czv = czv + Grhow(i,1,k)*sgn_jmin
           enddo
        enddo

     endif

     ! For walls at jmax index
     ! -----------------------
     if (is_bc_wall(2,2)) then

        ! starting indices
        i1=iis(1)
        k1=iks(5)

        ! contribution of pressure to the force
        do k=k1,nz
           do i=i1,nx
              cd = cd + prs(i,ny,k)*nxnds_jmax(i,k)
              cl = cl + prs(i,ny,k)*nynds_jmax(i,k)
              cz = cz + prs(i,ny,k)*nznds_jmax(i,k)
           enddo
        enddo

        ! contribution of viscous stresses to the force
        ! tau.nx=tau11*eta_x+tau12*eta_y+tau13*eta_z -> Grhou in visc.
        ! tau.ny=tau12*eta_x+tau22*eta_y+tau23*eta_z -> Grhov in visc.
        ! tau.nz=tau13*eta_x+tau23*eta_y+tau33*eta_z -> Grhow in visc.
        do k=k1,nz
           do i=i1,nx
              cdv = cdv + Grhou(i,ny,k)*sgn_jmax
              clv = clv + Grhov(i,ny,k)*sgn_jmax
              czv = czv + Grhow(i,ny,k)*sgn_jmax
           enddo
        enddo

     endif

     ! For walls at kmin index
     ! -----------------------
     if (is_bc_wall(3,1)) then

        ! starting indices
        i1=iis(1)
        j1=ijs(3)

        ! contribution of pressure to the force
        do i=i1,nx
           do j=j1,ny
              cd = cd + prs(i,j,1)*nxnds_kmin(i,j)
              cl = cl + prs(i,j,1)*nynds_kmin(i,j)
              cz = cz + prs(i,j,1)*nznds_kmin(i,j)
           enddo
        enddo

        ! contribution of viscous stresses to the force
        ! tau.nx=tau11*phi_x+tau12*phi_y+tau13*phi_z -> Hrhou in visc.
        ! tau.ny=tau12*phi_x+tau22*phi_y+tau23*phi_z -> Hrhov in visc.
        ! tau.nz=tau13*phi_x+tau23*phi_y+tau33*phi_z -> Hrhow in visc.
        do i=i1,nx
           do j=j1,ny
              cdv = cdv + Hrhou(i,j,1)*sgn_kmin
              clv = clv + Hrhov(i,j,1)*sgn_kmin
              czv = czv + Hrhow(i,j,1)*sgn_kmin
           enddo
        enddo

     endif

     ! For walls at kmax index
     ! -----------------------
     if (is_bc_wall(3,2)) then

        ! starting indices
        i1=iis(1)
        j1=ijs(3)

        ! contribution of pressure to the force
        do i=i1,nx
           do j=j1,ny
              cd = cd + prs(i,j,nz)*nxnds_kmax(i,j)
              cl = cl + prs(i,j,nz)*nynds_kmax(i,j)
              cz = cz + prs(i,j,nz)*nznds_kmax(i,j)
           enddo
        enddo

        ! contribution of viscous stresses to the force
        ! tau.nx=tau11*phi_x+tau12*phi_y+tau13*phi_z -> Hrhou in visc.
        ! tau.ny=tau12*phi_x+tau22*phi_y+tau23*phi_z -> Hrhov in visc.
        ! tau.nz=tau13*phi_x+tau23*phi_y+tau33*phi_z -> Hrhow in visc.
        do i=i1,nx
           do j=j1,ny
              cdv = cdv + Hrhou(i,j,nz)*sgn_kmax
              clv = clv + Hrhov(i,j,nz)*sgn_kmax
              czv = czv + Hrhow(i,j,nz)*sgn_kmax
           enddo
        enddo

     endif

  else ! half-cell approach (regular)
       ! ============================

     ! For walls at imin index
     ! -----------------------
     if (is_bc_wall(1,1)) then

        ! contribution of pressure to the force
        do k=1,nz
           do j=1,ny
              cd = cd + prs(1,j,k)*nxnds_imin(j,k)
              cl = cl + prs(1,j,k)*nynds_imin(j,k)
              cz = cz + prs(1,j,k)*nznds_imin(j,k)
           enddo
        enddo

        ! contribution of viscous stresses to the force
        ! tau.nx=tau11*ksi_x+tau12*ksi_y+tau13*ksi_z -> Frhou in visc.
        ! tau.ny=tau12*ksi_x+tau22*ksi_y+tau23*ksi_z -> Frhov in visc.
        ! tau.nz=tau13*ksi_x+tau23*ksi_y+tau33*ksi_z -> Frhow in visc.
        do k=1,nz
           do j=1,ny
              cdv = cdv + Frhou(1,j,k)*sgn_imin
              clv = clv + Frhov(1,j,k)*sgn_imin
              czv = czv + Frhow(1,j,k)*sgn_imin
           enddo
        enddo

     endif

     ! For walls at imax index
     ! -----------------------
     if (is_bc_wall(1,2)) then

        ! contribution of pressure to the force
        do k=1,nz
           do j=1,ny
              cd = cd + prs(nx,j,k)*nxnds_imax(j,k)
              cl = cl + prs(nx,j,k)*nynds_imax(j,k)
              cz = cz + prs(nx,j,k)*nznds_imax(j,k)
           enddo
        enddo

        ! contribution of viscous stresses to the force
        ! tau.nx=tau11*ksi_x+tau12*ksi_y+tau13*ksi_z -> Frhou in visc.
        ! tau.ny=tau12*ksi_x+tau22*ksi_y+tau23*ksi_z -> Frhov in visc.
        ! tau.nz=tau13*ksi_x+tau23*ksi_y+tau33*ksi_z -> Frhow in visc.
        do k=1,nz
           do j=1,ny
              cdv = cdv + Frhou(nx,j,k)*sgn_imax
              clv = clv + Frhov(nx,j,k)*sgn_imax
              czv = czv + Frhow(nx,j,k)*sgn_imax
           enddo
        enddo

     endif

     ! For walls at jmin index
     ! -----------------------
     if (is_bc_wall(2,1)) then

        ! contribution of pressure to the force
        do k=1,nz
           do i=1,nx
              cd = cd + prs(i,1,k)*nxnds_jmin(i,k)
              cl = cl + prs(i,1,k)*nynds_jmin(i,k)
              cz = cz + prs(i,1,k)*nznds_jmin(i,k)
           enddo
        enddo

        ! contribution of viscous stresses to the force
        ! tau.nx=tau11*eta_x+tau12*eta_y+tau13*eta_z -> Grhou in visc.
        ! tau.ny=tau12*eta_x+tau22*eta_y+tau23*eta_z -> Grhov in visc.
        ! tau.nz=tau13*eta_x+tau23*eta_y+tau33*eta_z -> Grhow in visc.
        do k=1,nz
           do i=1,nx
              cdv = cdv + Grhou(i,1,k)*sgn_jmin
              clv = clv + Grhov(i,1,k)*sgn_jmin
              czv = czv + Grhow(i,1,k)*sgn_jmin
           enddo
        enddo

     endif

     ! For walls at jmax index
     ! -----------------------
     if (is_bc_wall(2,2)) then

        ! contribution of pressure to the force
        do k=1,nz
           do i=1,nx
              cd = cd + prs(i,ny,k)*nxnds_jmax(i,k)
              cl = cl + prs(i,ny,k)*nynds_jmax(i,k)
              cz = cz + prs(i,ny,k)*nznds_jmax(i,k)
           enddo
        enddo

        ! contribution of viscous stresses to the force
        ! tau.nx=tau11*eta_x+tau12*eta_y+tau13*eta_z -> Grhou in visc.
        ! tau.ny=tau12*eta_x+tau22*eta_y+tau23*eta_z -> Grhov in visc.
        ! tau.nz=tau13*eta_x+tau23*eta_y+tau33*eta_z -> Grhow in visc.
        do k=1,nz
           do i=1,nx
              cdv = cdv + Grhou(i,ny,k)*sgn_jmax
              clv = clv + Grhov(i,ny,k)*sgn_jmax
              czv = czv + Grhow(i,ny,k)*sgn_jmax
           enddo
        enddo

     endif

     ! For walls at kmin index
     ! -----------------------
     if (is_bc_wall(3,1)) then

        ! contribution of pressure to the force
        do i=1,nx
           do j=1,ny
              cd = cd + prs(i,j,1)*nxnds_kmin(i,j)
              cl = cl + prs(i,j,1)*nynds_kmin(i,j)
              cz = cz + prs(i,j,1)*nznds_kmin(i,j)
           enddo
        enddo

        ! contribution of viscous stresses to the force
        ! tau.nx=tau11*phi_x+tau12*phi_y+tau13*phi_z -> Hrhou in visc.
        ! tau.ny=tau12*phi_x+tau22*phi_y+tau23*phi_z -> Hrhov in visc.
        ! tau.nz=tau13*phi_x+tau23*phi_y+tau33*phi_z -> Hrhow in visc.
        do i=1,nx
           do j=1,ny
              cdv = cdv + Hrhou(i,j,1)*sgn_kmin
              clv = clv + Hrhov(i,j,1)*sgn_kmin
              czv = czv + Hrhow(i,j,1)*sgn_kmin
           enddo
        enddo

     endif

     ! For walls at kmax index
     ! -----------------------
     if (is_bc_wall(3,2)) then

        ! contribution of pressure to the force
        do i=1,nx
           do j=1,ny
              cd = cd + prs(i,j,nz)*nxnds_kmax(i,j)
              cl = cl + prs(i,j,nz)*nynds_kmax(i,j)
              cz = cz + prs(i,j,nz)*nznds_kmax(i,j)
           enddo
        enddo

        ! contribution of viscous stresses to the force
        ! tau.nx=tau11*phi_x+tau12*phi_y+tau13*phi_z -> Hrhou in visc.
        ! tau.ny=tau12*phi_x+tau22*phi_y+tau23*phi_z -> Hrhov in visc.
        ! tau.nz=tau13*phi_x+tau23*phi_y+tau33*phi_z -> Hrhow in visc.
        do i=1,nx
           do j=1,ny
              cdv = cdv + Hrhou(i,j,nz)*sgn_kmax
              clv = clv + Hrhov(i,j,nz)*sgn_kmax
              czv = czv + Hrhow(i,j,nz)*sgn_kmax
           enddo
        enddo

     endif

  endif

  ! Sum over all walls and processes
  ! ================================
  call MPI_ALLREDUCE(MPI_IN_PLACE,cl,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
  call MPI_ALLREDUCE(MPI_IN_PLACE,cd,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
  call MPI_ALLREDUCE(MPI_IN_PLACE,clv,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
  call MPI_ALLREDUCE(MPI_IN_PLACE,cdv,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)

  ! Normalize
  ! =========
  ! adim for sphere case
  fac=2.0_wp/(rho_ref*u_ref**2*L_ref*L_ref*pi/4.0_wp)

  cd=cd*fac
  cl=cl*fac
  cz=cz*fac
  cdv=cdv*fac
  clv=clv*fac
  czv=czv*fac

  ! Write  lift and drag [-> in coeff.dat opened in stats.f90]
  ! ====================
  if (iproc==0) write(199,'(i12,1x,8(f15.10,1x))') ntotal,time,tstar,cl,cd,cz,clv,cdv,czv

end subroutine compute_cl_cd3

!===============================================================================
subroutine compute_cl_cd
!===============================================================================
  !> Compute lift and drag coefficients in 3D curvilinear
!===============================================================================
  use mod_mpi_part
  use mod_mpi_types_two_sided ! for indices iis,ijs,iks
  use mod_constant
  use mod_flow
  use mod_time
  implicit none
  ! ----------------------------------------------------------------------------
  integer :: i,j,k,i1,j1,k1
  real :: cl,cd,clv,cdv,fac
  ! ----------------------------------------------------------------------------

  ! Initialize lift and drag (for all procs)
  ! ========================
  cl=0.0_wp
  cd=0.0_wp
  clv=0.0_wp
  cdv=0.0_wp

  if (is_adjoint_block) then ! adjoint-block approach
                             ! ======================

     ! To avoid counting two times interface points, we can use iis, ijs & iks
     ! -----------------------------------------------------------------------
     ! Rq: indices declared and filled in mod_mpi_type_two_sided
     !     it takes the value 2 at block interfaces
     ! along i: iis(1) is the sending index for MPI comm with adjoint blocks
     ! along j: ijs(3) is the sending index for MPI comm with adjoint blocks
     ! along k: iks(5) is the sending index for MPI comm with adjoint blocks

     ! For walls at imin index
     ! -----------------------
     if (is_bc_wall(1,1)) then
        ! starting indices
        j1=ijs(3)
        k1=iks(5)

        ! contribution of pressure to the force
        do k=k1,nz
           do j=j1,ny
              cd = cd + prs(1,j,k)*nxndl_imin(j)
              cl = cl + prs(1,j,k)*nyndl_imin(j)
           enddo
        enddo

        ! contribution of viscous stresses to the force
        ! tau.nx=tau11*ksi_x+tau12*ksi_y -> Frhou in visc.
        ! tau.ny=tau12*ksi_x+tau22*ksi_y -> Frhov in visc.
        do k=k1,nz
           do j=j1,ny
              cdv = cdv + Frhou(1,j,k)*sgn_imin
              clv = clv + Frhov(1,j,k)*sgn_imin
           enddo
        enddo
     endif

     ! For walls at imax index
     ! -----------------------
     if (is_bc_wall(1,2)) then
        ! starting indices
        j1=ijs(3)
        k1=iks(5)

        ! contribution of pressure to the force
        do k=k1,nz
           do j=j1,ny
              cd = cd + prs(nx,j,k)*nxndl_imax(j)
              cl = cl + prs(nx,j,k)*nyndl_imax(j)
           enddo
        enddo

        ! contribution of viscous stresses to the force
        ! tau.nx=tau11*ksi_x+tau12*ksi_y -> Frhou in visc.
        ! tau.ny=tau12*ksi_x+tau22*ksi_y -> Frhov in visc.
        do k=k1,nz
           do j=j1,ny
              cdv = cdv + Frhou(nx,j,k)*sgn_imax
              clv = clv + Frhov(nx,j,k)*sgn_imax
           enddo
        enddo
     endif

     ! For walls at jmin index
     ! -----------------------
     if (is_bc_wall(2,1)) then
        ! starting indices
        i1=iis(1)
        k1=iks(5)

        ! contribution of pressure to the force
        do k=k1,nz
           do i=i1,nx
              cd = cd + prs(i,1,k)*nxndl_jmin(i)
              cl = cl + prs(i,1,k)*nyndl_jmin(i)
           enddo
        enddo

        ! contribution of viscous stresses to the force
        ! Frhou=tau11; Frhov=tau12
        ! Grhou=tau12; Grhov=tau22
        do k=k1,nz
           do i=i1,nx
              cdv = cdv -Frhou(i,1,k)*nxndl_jmin(i)-Frhov(i,1,k)*nyndl_jmin(i)
              clv = clv -Grhou(i,1,k)*nxndl_jmin(i)-Grhov(i,1,k)*nyndl_jmin(i)
           enddo
        enddo
     endif

     ! For walls at jmax index
     ! -----------------------
     if (is_bc_wall(2,2)) then
        ! starting indices
        i1=iis(1)
        k1=iks(5)

        ! contribution of pressure to the force
        do k=k1,nz
           do i=i1,nx
              cd = cd + prs(i,ny,k)*nxndl_jmax(i)
              cl = cl + prs(i,ny,k)*nyndl_jmax(i)
           enddo
        enddo

        ! contribution of viscous stresses to the force
        ! tau.nx=tau11*eta_x+tau12*eta_y -> Grhou in visc.
        ! tau.ny=tau12*eta_x+tau22*eta_y -> Grhov in visc.
        do k=k1,nz
           do i=i1,nx
              cdv = cdv + Grhou(i,ny,k)*sgn_jmax
              clv = clv + Grhov(i,ny,k)*sgn_jmax
           enddo
        enddo
     endif

  else ! half-cell approach (regular)
       ! ============================

     ! For walls at imin index
     ! -----------------------
     if (is_bc_wall(1,1)) then
        ! contribution of pressure to the force
        do k=1,nz
           do j=1,ny
              cd = cd + prs(1,j,k)*nxndl_imin(j)
              cl = cl + prs(1,j,k)*nyndl_imin(j)
           enddo
        enddo

        ! contribution of viscous stresses to the force
        ! tau.nx=tau11*ksi_x+tau12*ksi_y -> Frhou in visc.
        ! tau.ny=tau12*ksi_x+tau22*ksi_y -> Frhov in visc.
        do k=1,nz
           do j=1,ny
              cdv = cdv + Frhou(1,j,k)*sgn_imin
              clv = clv + Frhov(1,j,k)*sgn_imin
           enddo
        enddo
     endif

     ! For walls at imax index
     ! -----------------------
     if (is_bc_wall(1,2)) then
        ! contribution of pressure to the force
        do k=1,nz
           do j=1,ny
              cd = cd + prs(nx,j,k)*nxndl_imax(j)
              cl = cl + prs(nx,j,k)*nyndl_imax(j)
           enddo
        enddo

        ! contribution of viscous stresses to the force
        ! tau.nx=tau11*ksi_x+tau12*ksi_y -> Frhou in visc.
        ! tau.ny=tau12*ksi_x+tau22*ksi_y -> Frhov in visc.
        do k=1,nz
           do j=1,ny
              cdv = cdv + Frhou(nx,j,k)*sgn_imax
              clv = clv + Frhov(nx,j,k)*sgn_imax
           enddo
        enddo
     endif

     ! For walls at jmin index
     ! -----------------------
     if (is_bc_wall(2,1)) then
        ! contribution of pressure to the force
        do k=1,nz
           do i=1,nx
              cd = cd + prs(i,1,k)*nxndl_jmin(i)
              cl = cl + prs(i,1,k)*nyndl_jmin(i)
           enddo
        enddo

        ! contribution of viscous stresses to the force
        ! tau.nx=tau11*eta_x+tau12*eta_y -> Grhou in visc.
        ! tau.ny=tau12*eta_x+tau22*eta_y -> Grhov in visc.
        do k=1,nz
           do i=1,nx
              cdv = cdv + Grhou(i,1,k)*sgn_jmin
              clv = clv + Grhov(i,1,k)*sgn_jmin
           enddo
        enddo
     endif

     ! For walls at jmax index
     ! -----------------------
     if (is_bc_wall(2,2)) then
        ! contribution of pressure to the force
        do k=1,nz
           do i=1,nx
              cd = cd + prs(i,ny,k)*nxndl_jmax(i)
              cl = cl + prs(i,ny,k)*nyndl_jmax(i)
           enddo
        enddo

        ! contribution of viscous stresses to the force
        ! tau.nx=tau11*eta_x+tau12*eta_y -> Grhou in visc.
        ! tau.ny=tau12*eta_x+tau22*eta_y -> Grhov in visc.
        do k=1,nz
           do i=1,nx
              cdv = cdv + Grhou(i,ny,k)*sgn_jmax
              clv = clv + Grhov(i,ny,k)*sgn_jmax
           enddo
        enddo
     endif

  endif

  ! Sum over processes
  ! ==================
  call MPI_ALLREDUCE(MPI_IN_PLACE,cl,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)  
  call MPI_ALLREDUCE(MPI_IN_PLACE,cd,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)  
  call MPI_ALLREDUCE(MPI_IN_PLACE,clv,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)
  call MPI_ALLREDUCE(MPI_IN_PLACE,cdv,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_global,info)

  ! Normalize
  ! =========

  ! adim for cylinder case
  fac=2.0_wp/(rho_ref*u_ref**2*L_ref*ngz)

  cd=cd*fac
  cl=cl*fac
  cdv=cdv*fac
  clv=clv*fac
  
  ! Write  lift and drag [-> in coeff.dat opened in stats.f90]
  ! ====================
  if (iproc==0) write(199,'(i12,1x,6(f15.10,1x))') ntotal,time,tstar,cl,cd,clv,cdv
 
end subroutine compute_cl_cd

!===============================================================================
subroutine write_surf_norm
!===============================================================================
  !> Write normal vector to surface to file
!===============================================================================
  use mod_mpi
  use mod_mpi_part
  use mod_constant
  use mod_flow
  use mod_time
  implicit none
  ! ----------------------------------------------------------------------------
  integer :: i,j,k,n,ntot,norm_file=200,planes,every,yn
  character (1) :: n_str
  real, dimension(:), allocatable :: norm,norm_g
  logical, dimension(:), allocatable :: walls
  ! ----------------------------------------------------------------------------

  n = 6
  planes = ngz/nz
  allocate(norm(n*nx),norm_g(n*nproc*nx))
  norm = 0.0_wp
  norm_g = 0.0_wp

  ! Check how many procs have walls at jmin /!\ modify for walls at all locations /!\
  allocate(walls(nproc))
  call MPI_ALLGATHER(is_bc_wall(2,1),1,MPI_LOGICAL,walls,1,MPI_LOGICAL,COMM_global,info)

  ! Write norms to arrays
  if (is_bc_wall(2,1)) then
    k = 0
    do i=1,nx
       norm(i+k)   = nxndl_jmin(i)
       norm(i+k+1) = nyndl_jmin(i)
       norm(i+k+2) = nxn_jmin(i,1)
       norm(i+k+3) = nyn_jmin(i,1)
       norm(i+k+4) = dl_jmin(i)
       k = k+(n-1)
    enddo
  endif

  ! Gather and make proc 0 write to file
  call MPI_GATHER(norm,nx*n,MPI_DOUBLE_PRECISION,norm_g,nx*n,MPI_DOUBLE_PRECISION,0,COMM_global,info)
  if (iproc.eq.0) then
     print*,'Writing surface normals'
     write(n_str,'(i1)') n
     open(norm_file,file='norm_surf.dat',form='formatted',status='replace')
     k = 0
     every = 0
     ntot = nproc*nx
     do i=1,ntot
        yn = (i-1)/nx+1 ! yes/no -> index of table "walls" to skip procs without one
        if (walls(yn)) then
           write(norm_file,'('//trim(n_str)//'(f15.10,1x))') (/ (norm_g(i+j+k+every*nx*n),j=0,(n-1)) /)
        endif
        k = k+(n-1)
        if (planes.gt.1.and.mod(i,(planes-1)*nx).eq.0) every = every+1
     enddo
  endif

end subroutine write_surf_norm

!===============================================================================
subroutine compute_nusselt
!===============================================================================
  !> Compute spatially averaged Nusselt
!===============================================================================
  use mod_mpi
  use mod_constant
  use mod_flow
  use mod_time
  ! ----------------------------------------------------------------------------
  implicit none
  integer :: i,k
  real :: nu,nu_g
  ! ----------------------------------------------------------------------------

  nu = 0.0_wp
  nu_g = 0.0_wp

  if (is_bc_wall(2,1)) then
    do i=1,nx
       do k=1,nz
          nu = nu - (dTx(i,1,k)*nxndl_jmin(i)+dTy(i,1,k)*nyndl_jmin(i))/pi/(T_ref-Tmp(i,1,k))/nz
       enddo
    enddo
  endif

  call MPI_REDUCE(nu,nu_g,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMM_global,info)
  if (iproc==0) write(201,'(i12,1x,3(f20.10,1x))') ntotal,time,tstar,nu_g
  if (iproc==0.and.ntotal==1) print*,'L_ref =',L_ref

end subroutine compute_nusselt

!===============================================================================
subroutine compute_vorticity2d(k)
!===============================================================================
  !> Compute 2D vorticity for check
!===============================================================================
  use mod_flow
  implicit none
  ! ----------------------------------------------------------------------------
  integer, intent(in) :: k
  ! ----------------------------------------------------------------------------
  integer :: i,j
  ! ----------------------------------------------------------------------------
    
  ! Compute 2-D vorticity=dv/dx-du/dy
  ! ================================
  ! use same derivatives calculated for viscous fluxes
  uvar(:,:,:,1)=0.
  
  do j=1,ny
     do i=1,nx
        uvar(i,j,k,1)=dvx(i,j,k)-duy(i,j,k) 
     enddo
  enddo

end subroutine compute_vorticity2d

