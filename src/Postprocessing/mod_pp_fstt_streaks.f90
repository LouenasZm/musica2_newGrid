!==============================================================================
module mod_pp_fstt_streaks
!==============================================================================
  !> Module for streaks detection in laminar region
!==============================================================================
  use warnstop
  use mod_constant
  use mod_flow
  use mod_mpi
  use mod_pp_var
  use mod_pp_fstt_comm
  implicit none
  ! -------------------------------------------------------------------------------------------
  ! STREAK TYPE defining its attributes
  ! -------------------------------------------------------------------------------------------
  type streak_gath_type  ! for preliminary creation of streaks
     integer, dimension(:), allocatable :: active     ! active (1) or not (0) for all streaks
     integer, dimension(:), allocatable :: length     ! length <~ number of pts for all streaks
     real(wp), dimension(:,:), allocatable :: data    ! positions in x=(:,1), y=(:,2), z=(:,3)
                                                      ! u' in (:,4)
  end type streak_gath_type
  type final_streak_type
     real(wp), dimension(:,:), allocatable :: data    ! positions in x=(:,1), y=(:,2), z=(:,3)
                                                      ! u' in (:,4)
  end type final_streak_type
  ! ---------------------------------------------------------------------------
  ! type(streak_type), dimension(:), allocatable :: streaks
  ! type(streak_type), dimension(:), allocatable :: streaks_tot
  type(streak_gath_type) :: streaks
  type(streak_gath_type), dimension(:), allocatable :: streaks_gather
  type(final_streak_type), dimension(:), allocatable :: streaks_merge,streaks_final
  ! -------------------------------------------------------------------------

contains

  !============================================================================
  subroutine init_streaks_detection
  !============================================================================
    !> Initialisation for streaks detection
    !============================================================================
    use mod_utils
    implicit none
    ! -------------------------------------------------------------------------
    logical :: iexist
    ! -------------------------------------------------------------------------

    ! East neighbors
    allocate(neighbor_E_fstt(0:bl(1)%nproc-1))
    call MPI_ALLGATHER(neighbor(nE),1,MPI_INTEGER,& ! send
                        neighbor_E_fstt,1,MPI_INTEGER,& ! receive
                        COMM_global,info)
    if (iproc.ne.0) deallocate(neighbor_E_fstt)

    ! If streaks file exist, suppressed if not is_restart_pp
    if ((iproc==0).and.(.not.is_restart_pp)) then
       inquire(file='streaks_bl'//trim(numchar(nob(iproc)))//'.bin',exist=iexist)
       if (iexist) then
          call system('rm streaks_bl'//trim(numchar(nob(iproc)))//'.bin')
       endif
    endif

  end subroutine init_streaks_detection

  !============================================================================
  subroutine local_streaks_detection
  !============================================================================
    !> Main subroutine for streaks detection
    !============================================================================
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i,j,k,iprime,l,lp,nstreaks_p2_loc,nstreaks_p2_tot,nstreaks_p,nstreaks_v,ind_link,length_tot,ibeg
    real(wp) :: ref_min,ref_max,dist,dist_test,z_decal,Lz_fstt,dist_min_perio
    real(wp) :: y_fstt_j,y_fstt_j2,z_fstt_k,z_fstt_k2
    logical :: is_kept
    real(wp), dimension(4,4*ngz) :: temp_streaks_p2
    real(wp), dimension(4,4*ngz) :: temp_streaks_p
    real(wp), dimension(2*nx*ngz,-1:nx,1:4) :: temp_streaks_v
    ! -------------------------------------------------------------------------
    integer, dimension(MPI_STATUS_SIZE) :: status_
    ! -------------------------------------------------------------------------

    ! Initialisation
    ! --------------
    temp_streaks_v = 0.0_wp
    nstreaks_v = 0

        do i=1,nx
       ! Initialisation
       ! --------------
       temp_streaks_p = 0
       nstreaks_p = 0

       ! Indice of association put to 0
       ! ------------------------------
       do l=1,nstreaks_v
          temp_streaks_v(l,-1,3) = 0
       enddo

       ref_min = -0.005_wp*U_ref
       ref_max =  0.005_wp*U_ref

       ! Detection of extremums
       ! ----------------------
       do k=1,nz
          do j=2,j_edge2(i,k)-1
             if (ltbl(i,j,k).gt.0.5_wp) then
                continue
             else if (uu_fluct(i,j,k).eq.max(maxval(uu_fluct(i,j-1:j+1,k-1:k+1)),ref_max)) then
                nstreaks_p = nstreaks_p + 1
                temp_streaks_p(1,nstreaks_p) = y(j)
                temp_streaks_p(2,nstreaks_p) = z(k)
                temp_streaks_p(3,nstreaks_p) = 1.0_wp
                temp_streaks_p(4,nstreaks_p) = uu_fluct(i,j,k)

             else if (uu_fluct(i,j,k).eq.min(minval(uu_fluct(i,j-1:j+1,k-1:k+1)),ref_min)) then
                nstreaks_p = nstreaks_p + 1
                temp_streaks_p(1,nstreaks_p) = y(j)
                temp_streaks_p(2,nstreaks_p) = z(k)
                temp_streaks_p(3,nstreaks_p) = -1.0_wp
                temp_streaks_p(4,nstreaks_p) = uu_fluct(i,j,k)
             endif
          enddo
       enddo

       ! Elimination of extremums too close from each other
       ! --------------------------------------------------
       temp_streaks_p2 = 0
       temp_streaks_p2(:,1) = temp_streaks_p(:,1)
       nstreaks_p2_loc = min(1,nstreaks_p)
       do l=1,nstreaks_p
          y_fstt_j=temp_streaks_p(1,l);z_fstt_k=temp_streaks_p(2,l)
          is_kept = .true.
          do lp=1,nstreaks_p2_loc
             if (temp_streaks_p(3,l).eq.temp_streaks_p2(3,lp)) then
                y_fstt_j2=temp_streaks_p2(1,lp);z_fstt_k2=temp_streaks_p2(2,lp)
                dist = ((y_fstt_j - y_fstt_j2)**2 +&
                        (z_fstt_k - z_fstt_k2)**2)**0.5
                if (dist.lt.dist_max_extr) then
                   if (abs(temp_streaks_p2(4,lp)).ge.abs(temp_streaks_p(4,l))) then
                      is_kept = .false.
                      temp_streaks_p2(:,lp) = temp_streaks_p(:,l)
                   else
                      is_kept = .false.
                   endif
                endif
             endif
          enddo
          if (is_kept) then
             nstreaks_p2_loc = nstreaks_p2_loc + 1
             temp_streaks_p2(:,nstreaks_p2_loc) = temp_streaks_p(:,l)
          endif
       enddo

       if (coord(3).eq.0) then
          nstreaks_p2_tot = nstreaks_p2_loc
          ! Communication of extremums to master coord3
          do k=1,bl(1)%ndomk-1
             call MPI_RECV(nstreaks_p2_loc,1,MPI_INTEGER,&
                           k + coord(1)*bl(1)%ndomk,tag,COMM_global,status_,info)
             call MPI_RECV(temp_streaks_p2(:,nstreaks_p2_tot+1:nstreaks_p2_tot+nstreaks_p2_loc),nstreaks_p2_loc*4,&
                           MPI_DOUBLE_PRECISION,k+coord(1)*bl(1)%ndomk,tag,COMM_global,status_,info)
             nstreaks_p2_tot = nstreaks_p2_tot + nstreaks_p2_loc
          enddo
       else
          call MPI_SEND(nstreaks_p2_loc,1,MPI_INTEGER,&
                        coord(1)*bl(1)%ndomk,tag,COMM_global,info)
          call MPI_SEND(temp_streaks_p2(:,1:nstreaks_p2_loc),nstreaks_p2_loc*4,&
                        MPI_DOUBLE_PRECISION,coord(1)*bl(1)%ndomk,tag,COMM_global,info)
       endif

       if (coord(3).eq.0) then
          ! Elimination of extremums too close from each other (rest from data received)
          ! --------------------------------------------------
          temp_streaks_p = 0
          temp_streaks_p(:,1) = temp_streaks_p2(:,1)
          nstreaks_p = min(1,nstreaks_p2_tot)
          do l=1,nstreaks_p2_tot
             y_fstt_j=temp_streaks_p2(1,l);z_fstt_k=temp_streaks_p2(2,l)
             is_kept = .true.
             do lp=1,nstreaks_p
                if (temp_streaks_p2(3,l).eq.temp_streaks_p(3,lp)) then
                   y_fstt_j2=temp_streaks_p(1,lp);z_fstt_k2=temp_streaks_p(2,lp)
                   dist = ((y_fstt_j - y_fstt_j2)**2 +&
                           (z_fstt_k - z_fstt_k2)**2)**0.5
                   if (dist.lt.dist_max_extr) then
                      if (abs(temp_streaks_p(4,lp)).ge.abs(temp_streaks_p2(4,l))) then
                         is_kept = .false.
                         temp_streaks_p(:,lp) = temp_streaks_p2(:,l)
                      else
                         is_kept = .false.
                      endif
                   endif
                endif
             enddo
             if (is_kept) then
                nstreaks_p = nstreaks_p + 1
                temp_streaks_p(:,nstreaks_p) = temp_streaks_p2(:,l)
             endif
          enddo


          ! Association of extremums to streaks
          ! -----------------------------------
          Lz_fstt = zg(ngz) - zg(1)
          dist_min_perio = Lz_fstt + (z(2) - z(1)) - dist_max_extr

          ! case i=1: initialisation of streaks
          if (i.eq.1) then
             nstreaks_v = nstreaks_p
             do l=1,nstreaks_v
                ! Creation of each streaks
                temp_streaks_v(l,-1,1)= 1 ! imin
                temp_streaks_v(l,-1,2)= 0 ! distance between between 2 last points of streaks (y-z plane)
                temp_streaks_v(l,-1,3)= 1 ! 1 if streaks already associated, 0 else
                temp_streaks_v(l,0,1) = 1 ! length (number of points)
                temp_streaks_v(l,0,2) = 1 ! active (1) or not (0)
                temp_streaks_v(l,0,3) = temp_streaks_p(3,l) ! type: positive (1) or negative (-1)
                temp_streaks_v(l,1,1) = x(i) ! position (imin,x(i))
                temp_streaks_v(l,1,2) = temp_streaks_p(1,l) ! position (imin,y(j))
                temp_streaks_v(l,1,3) = temp_streaks_p(2,l) ! position (imin,z(k))
                temp_streaks_v(l,1,4) = temp_streaks_p(4,l) ! u'
             enddo
          ! case i>1: association of current streaks with previous found
          else
             ! Do loop over all the extremums existing to calculate the distance between new and old streaks
             do lp=1,nstreaks_p
                y_fstt_j=temp_streaks_p(1,lp);z_fstt_k=temp_streaks_p(2,lp)
                dist = 1e6
                ind_link = -999
                do l=1,nstreaks_v
                   ! If streaks still active and same type (neg. or pos.)
                   if ((temp_streaks_v(l,0,2).eq.1).and.(temp_streaks_v(l,0,3).eq.temp_streaks_p(3,lp))) then
                      iprime = temp_streaks_v(l,-1,1) + temp_streaks_v(l,0,1) - 1
                      ! If streaks separated by a distance superior to Lz-dist_extr+dz, periodicity considered
                      if (((temp_streaks_v(l,iprime,3) - z_fstt_k)**2)**0.5.gt.dist_min_perio) then
                         if (temp_streaks_v(l,iprime,3).le.zg(ngz/2)) then
                            z_decal = temp_streaks_v(l,iprime,3) + Lz_fstt
                         else
                            z_decal = temp_streaks_v(l,iprime,3) - Lz_fstt
                         endif
                         dist_test = ((temp_streaks_v(l,iprime,2) - y_fstt_j)**2 + &
                                      (z_decal - z_fstt_k)**2)**0.5
                      else
                         dist_test = ((temp_streaks_v(l,iprime,2) - y_fstt_j)**2 + &
                                      (temp_streaks_v(l,iprime,3) - z_fstt_k)**2)**0.5
                      endif
                      ! If old streak closer to new streak than the one already found, replaced
                      if (dist.gt.dist_test) then
                         dist = dist_test
                         ind_link = l
                      endif
                   endif
                enddo

                ! To avoid associating extremum to streaks if too much distance
                ! Creation of a new streak if it is the case
                if (dist.gt.dist_max_yz) then
                   nstreaks_v = nstreaks_v + 1
                   l = nstreaks_v
                   ! Initialisation of streak
                   temp_streaks_v(l,-1,1)= i ! imin
                   temp_streaks_v(l,-1,2)= 0 ! distance between between 2 last points of streaks (y-z plane)
                   temp_streaks_v(l,-1,3)= 1 ! 1 if streaks already associated, 0 else
                   temp_streaks_v(l,0,1) = 1 ! length (number of points)
                   temp_streaks_v(l,0,2) = 1 ! active (1) or not (0)
                   temp_streaks_v(l,0,3) = temp_streaks_p(3,lp) ! type: positive (1) or negative (-1)
                   iprime = temp_streaks_v(l,-1,1) + temp_streaks_v(l,0,1) - 1
                   temp_streaks_v(l,iprime,1) = x(i) ! position (imin,x(i))
                   temp_streaks_v(l,iprime,2) = y_fstt_j ! position (imin,y(j))
                   temp_streaks_v(l,iprime,3) = z_fstt_k ! position (imin,z(k))
                   temp_streaks_v(l,iprime,4) = temp_streaks_p(4,lp) ! u'
                ! If close enough to a streak not associated yet, linkage between extremum and streak
                else if (temp_streaks_v(ind_link,-1,3).eq.0) then
                   l = ind_link
                   temp_streaks_v(l,-1,2)= dist ! distance between between 2 last points of streaks (y-z plane)
                   temp_streaks_v(l,-1,3)= 1 ! 1 if streaks already associated, 0 else
                   temp_streaks_v(l,0,1) = temp_streaks_v(l,0,1) + 1 ! length (number of points)
                   iprime = temp_streaks_v(l,-1,1) + temp_streaks_v(l,0,1) - 1
                   temp_streaks_v(l,iprime,1) = x(i) ! position (i,x(i))
                   temp_streaks_v(l,iprime,2) = y_fstt_j ! position (i,y(j))
                   temp_streaks_v(l,iprime,3) = z_fstt_k ! position (i,z(k))
                   temp_streaks_v(l,iprime,4) = temp_streaks_p(4,lp) ! u'
                ! If streak already associated with an extremum
                else
                   nstreaks_v = nstreaks_v + 1
                   temp_streaks_v(nstreaks_v,-1,1)= i ! imin
                   temp_streaks_v(nstreaks_v,-1,2)= 0 ! distance between between 2 last points of streaks (y-z plane)
                   temp_streaks_v(nstreaks_v,-1,3)= 1 ! 1 if streaks already associated, 0 else
                   temp_streaks_v(nstreaks_v,0,1) = 1 ! length (number of points)
                   temp_streaks_v(nstreaks_v,0,2) = 1 ! active (1) or not (0)
                   l = ind_link
                   ! 1st situation: this extremum is closer to the streak than the previous one
                   if (dist.lt.temp_streaks_v(l,-1,2)) then
                      iprime = temp_streaks_v(l,-1,1) + temp_streaks_v(l,0,1) - 1
                      ! New streak associated with the previous extremum
                      temp_streaks_v(nstreaks_v,0,3) = temp_streaks_v(l,0,3) ! type: positive (1) or negative (-1)
                      temp_streaks_v(nstreaks_v,iprime,1) = temp_streaks_v(l,iprime,1) ! position (i,x(i))
                      temp_streaks_v(nstreaks_v,iprime,2) = temp_streaks_v(l,iprime,2) ! position (i,y(j))
                      temp_streaks_v(nstreaks_v,iprime,3) = temp_streaks_v(l,iprime,3) ! position (i,z(k))
                      temp_streaks_v(nstreaks_v,iprime,4) = temp_streaks_v(l,iprime,4) ! u'
                      ! Last values for streak replaced by actual extremum
                      temp_streaks_v(l,-1,2)= dist ! distance between between 2 last points of streaks (y-z plane)
                      temp_streaks_v(l,-1,3)= 1 ! 1 if streaks already associated, 0 else
                      temp_streaks_v(l,0,3) = temp_streaks_p(3,lp) ! type: positive (1) or negative (-1)
                      temp_streaks_v(l,iprime,1) = x(i) ! position (i,x(i))
                      temp_streaks_v(l,iprime,2) = y_fstt_j ! position (i,y(j))
                      temp_streaks_v(l,iprime,3) = z_fstt_k ! position (i,z(k))
                      temp_streaks_v(l,iprime,4) = temp_streaks_p(4,lp) ! u'
                   ! 2nd situation: the previous extremum is closer to the streak than this one
                   else
                      iprime = temp_streaks_v(nstreaks_v,-1,1) + temp_streaks_v(nstreaks_v,0,1) - 1
                      ! New streak associated with the new extremum
                      temp_streaks_v(nstreaks_v,0,3) = temp_streaks_p(3,lp) ! type: positive (1) or negative (-1)
                      temp_streaks_v(nstreaks_v,iprime,1) = x(i) ! position (i,x(i))
                      temp_streaks_v(nstreaks_v,iprime,2) = y_fstt_j ! position (i,y(j))
                      temp_streaks_v(nstreaks_v,iprime,3) = z_fstt_k ! position (i,z(k))
                      temp_streaks_v(nstreaks_v,iprime,4) = temp_streaks_p(4,lp) ! u'
                   endif
                endif
             enddo
          endif

          ! Deactivation of streaks with no extremum associated
          do l=1,nstreaks_v
             if (temp_streaks_v(l,-1,3).eq.0) temp_streaks_v(l,0,2) = 0
          enddo

       endif
    enddo

    ! Initialisation
    nstreaks_loc = 0

    if (coord(3).eq.0) then
       ! Elimination of streaks with length inferior to l_min_str
       ! --------------------------------------------------------
       ! Count of the number of streaks long enough
       dist = x(2)-x(1) ! dx


       do l=1,nstreaks_v
          iprime = temp_streaks_v(l,-1,1)
          if (((temp_streaks_v(l,0,1)-1)*dist.ge.l_min_str).or.&
              (temp_streaks_v(l,iprime,1).eq.x(1)).or.&
              (temp_streaks_v(l,0,2).eq.1))        nstreaks_loc = nstreaks_loc +1
       enddo

       ! Allocation of streaks length array
       allocate(streaks%length(1:nstreaks_loc))
       allocate(streaks%active(1:nstreaks_loc))
       ! Filled with each streak length
       lp = 0
       do l=1,nstreaks_v
          iprime = temp_streaks_v(l,-1,1)
          if (((temp_streaks_v(l,0,1)-1)*dist.ge.l_min_str).or.&
              (temp_streaks_v(l,iprime,1).eq.x(1)).or.&
              (temp_streaks_v(l,0,2).eq.1)) then
             lp = lp + 1
             streaks%length(lp) = temp_streaks_v(l,0,1)
          endif
       enddo

       ! Calculation of the total cumulated length
       length_tot = sum(streaks%length)
       ! Allocation of streak%data to concatenate all streaks data
       allocate(streaks%data(length_tot,4))
       ! Filled with streaks values
       lp = 0
       length_tot = 0
       do l=1,nstreaks_v
          is_kept = .false.
          iprime = temp_streaks_v(l,-1,1)

          if (((temp_streaks_v(l,0,1)-1)*dist.ge.l_min_str).or.&
              (temp_streaks_v(l,iprime,1).eq.x(1)).or.&
              (temp_streaks_v(l,0,2).eq.1))   is_kept = .true.

          if (is_kept) then
             lp = lp + 1
             ! Save if streaks kept: _bc long enough (0)
             ! _bc point only on boundary W (1)
             ! _bc point only on boundary E (2)
             ! _bc point on boundary W and boundary E (3)
             if ((temp_streaks_v(l,iprime,1).eq.x(1)).and.(temp_streaks_v(l,0,2).eq.1)) then
                streaks%active(lp) = 3
             else if (temp_streaks_v(l,iprime,1).eq.x(1)) then
                streaks%active(lp) = 1
             else if (temp_streaks_v(l,0,2).eq.1) then
                streaks%active(lp) = 2
             else
                streaks%active(lp) = 0
             endif

             ! Save of the streak data
             ibeg = length_tot

             do i=1,streaks%length(lp)
                iprime = temp_streaks_v(l,-1,1) + i - 1
                streaks%data(i+ibeg,1) = temp_streaks_v(l,iprime,1)
                streaks%data(i+ibeg,2) = temp_streaks_v(l,iprime,2)
                streaks%data(i+ibeg,3) = temp_streaks_v(l,iprime,3)
                streaks%data(i+ibeg,4) = temp_streaks_v(l,iprime,4)
                length_tot = length_tot + 1
             enddo

          endif
       enddo
    endif

  end subroutine local_streaks_detection

  !============================================================================
  subroutine merging_streaks
  !============================================================================
    !> Communication of streaks to master proc & merge
    !============================================================================
    use mod_utils
    implicit none
    ! -------------------------------------------------------------------------
    integer :: i,j,ip,ip_c1,iprime,iprime2,l,lp,lp_save,ls,info,ibeg,iend,npoints
    integer :: length_tot,length_tot2,ind_link,nstreaks_tot
    real(wp) :: dist,dist_test,dist_test2,dist_tot,dist_min_perio
    real(wp) :: x_beg,x_end,y_beg,y_end,z_beg,z_end,Lz_fstt,z_decal
    logical :: has_to_test,switch
    integer, dimension(:), allocatable :: nstrk_per_proc,l_str_loc,activ_loc,table_link_f,assoc_f
    integer, dimension(:), allocatable :: to_be_tested,streaks_m_points,streaks_f_points
    real(wp), dimension(:,:), allocatable :: assoc
    ! -------------------------------------------------------------------------
    integer, dimension(MPI_STATUS_SIZE) :: status_
    ! -------------------------------------------------------------------------
    type streaks_link
       integer, dimension(:), allocatable :: link ! link of each streak with (proc,num. streak proc)
    end type streaks_link
    type(streaks_link), dimension(0:bl(1)%ndomi-1) :: table_link,streak_tip
    ! -------------------------------------------------------------------------

    ! Communication of number of streaks for each proc
    ! ------------------------------------------------
    allocate(nstrk_per_proc(0:bl(1)%nproc-1))
    if (iproc==0) then
       allocate(streaks_gather(0:bl(1)%ndomi-1))

       ! Receive nstreaks of each proc ~> stored in nstrk_per_proc
       call MPI_GATHER(nstreaks_loc,1,MPI_INTEGER,& ! send
                        nstrk_per_proc,1,MPI_INTEGER,& ! receive
                        0,COMM_global,info)

       ! Calculation of the total number of streaks
       nstreaks_tot = sum(nstrk_per_proc)
    else
       call MPI_GATHER(nstreaks_loc,1,MPI_INTEGER,& ! send
                        nstrk_per_proc,1,MPI_INTEGER,& ! receive
                        0,COMM_global,info)
       deallocate(nstrk_per_proc)
    endif


    ! Creation of an array with the length of each streaks
    ! ----------------------------------------------------
    if (coord(3).eq.0) then
       if (iproc.eq.0) then
          ! Data directly filled in streaks_gather for master
          allocate(streaks_gather(0)%length(nstrk_per_proc(0)))
          allocate(streaks_gather(0)%active(nstrk_per_proc(0)))
          do l=1,nstreaks_loc
             streaks_gather(0)%length(l) = streaks%length(l)
             streaks_gather(0)%active(l) = streaks%active(l)
          enddo
       else
          allocate(l_str_loc(nstreaks_loc))
          allocate(activ_loc(nstreaks_loc))
          do l=1,nstreaks_loc
             l_str_loc(l) = streaks%length(l)
             activ_loc(l) = streaks%active(l)
          enddo
       endif



       ! Communication of streaks data of each proc to master
       ! ----------------------------------------------------
       if (iproc.eq.0) then
          ! Communication of length and activity of each streak to master
          do i=1,bl(1)%ndomi-1
             allocate(streaks_gather(i)%length(nstrk_per_proc(i*bl(1)%ndomk)))
             call MPI_RECV(streaks_gather(i)%length,nstrk_per_proc(i*bl(1)%ndomk),MPI_INTEGER,&
                           i*bl(1)%ndomk,tag,COMM_global,status_,info)
             allocate(streaks_gather(i)%active(nstrk_per_proc(i*bl(1)%ndomk)))
             call MPI_RECV(streaks_gather(i)%active,nstrk_per_proc(i*bl(1)%ndomk),MPI_INTEGER,&
                           i*bl(1)%ndomk,tag,COMM_global,status_,info)
          enddo

          ! Data directly filled for master
          length_tot = sum(streaks_gather(0)%length)
          allocate(streaks_gather(0)%data(length_tot,4))
          streaks_gather(0)%data = streaks%data

          ! Communication of concatenated streaks data per proc
          do i=1,bl(1)%ndomi-1
             length_tot = sum(streaks_gather(i)%length)
             allocate(streaks_gather(i)%data(length_tot,4))
             call MPI_RECV(streaks_gather(i)%data,4*length_tot,MPI_DOUBLE_PRECISION,&
                           i*bl(1)%ndomk,tag,COMM_global,status_,info)
          enddo
       else
          ! Communication of length of each streak to master
          call MPI_SEND(l_str_loc,nstreaks_loc,MPI_INTEGER,&
                        0,tag,COMM_global,info)
          ! Communication of activity of each streak to master
          call MPI_SEND(activ_loc,nstreaks_loc,MPI_INTEGER,&
                        0,tag,COMM_global,info)

          ! Communication of concatenated streaks data per proc to master
          length_tot = sum(l_str_loc)
           call MPI_SEND(streaks%data,4*length_tot,MPI_DOUBLE_PRECISION,&
                        0,tag,COMM_global,info)

           deallocate(l_str_loc,activ_loc)
       endif

       ! Deallocation of streaks and length per streak
       deallocate(streaks%length,streaks%active,streaks%data)

       if (iproc==0) then
          Lz_fstt = zg(ngz) - zg(1)
          dist_min_perio = Lz_fstt + (z(2) - z(1)) - dist_max_extr

          ! 1/ Association between the streaks received by master
          ! -----------------------------------------------------
          loopproc: do i=0,bl(1)%ndomi-1
             ! West neighbor
             ip = neighbor_E_fstt(i*bl(1)%ndomk)
             ip_c1 = ip/bl(1)%ndomk

             ! allocation of the table of streak linkage accros processors
             allocate(table_link(i)%link(nstrk_per_proc(i*bl(1)%ndomk)))

             ! If coord(1) == 0 then all streaks are the tip
             if (i*bl(1)%ndomk.eq.0) then
                allocate(streak_tip(i)%link(nstrk_per_proc(i*bl(1)%ndomk)))
                streak_tip(i)%link(:) = 1
             endif

             ! if last proc in streamwise direction, no association to be made
             if (neighbor_E_fstt(i*bl(1)%ndomk).lt.0) then
                table_link(i)%link(:) = -1
                cycle loopproc
             endif

             ! To save if the neigboring streak is the beginning of the global streak
             allocate(streak_tip(ip_c1)%link(nstrk_per_proc(ip)))
             ! By default, yes (1)
             streak_tip(ip_c1)%link(:) = 1

             has_to_test = .true.
             ! To save if streak associated or not
             allocate(to_be_tested(nstrk_per_proc(i*bl(1)%ndomk)))
             ! Initially put to 1 = need to be tested
             to_be_tested(:) = 1

             ! To save the associations made
             allocate(assoc(nstrk_per_proc(ip),2))
             assoc(:,1) = -1
             assoc(:,2) = 1e6

             ! Determination of linkage until all association has been made
             do while (has_to_test)

             ! loop over all the streaks of the current proc
             ! to determine its linkage
             length_tot = 0
             loop1: do l=1,nstrk_per_proc(i*bl(1)%ndomk)
                ! Determination of the position in data of the current streak
                length_tot = length_tot + streaks_gather(i)%length(l)
                iprime = length_tot

                ! Verification if association has been made for this streak already
                if (to_be_tested(l).eq.0) cycle loop1

                ! Verification if the current streak can be associated or not
                if ((streaks_gather(i)%active(l).lt.2)) then
                   table_link(i)%link(l) = -1
                   to_be_tested(l) = 0
                   cycle loop1
                endif

                ! Loop over the following proc to determine the linkage
                dist = 1e6
                ind_link = -999
                length_tot2 = 0
                loop2: do lp=1,nstrk_per_proc(ip)
                   ! Determination of the position in data of the streak
                   length_tot2 = length_tot2 + streaks_gather(ip_c1)%length(lp)
                   iprime2 = length_tot2 - streaks_gather(ip_c1)%length(lp) + 1

                   ! If streak is not on west boundary, no association
                   if ((streaks_gather(ip_c1)%active(lp).ne.1).and.(streaks_gather(ip_c1)%active(lp).ne.3)) cycle loop2

                   ! If type of streak is different, no association
                   if ((streaks_gather(ip_c1)%data(iprime2,4) + streaks_gather(i)%data(iprime,4) .ne.&
                        abs(streaks_gather(ip_c1)%data(iprime2,4)) + abs(streaks_gather(i)%data(iprime,4))) .and. &
                       (streaks_gather(ip_c1)%data(iprime2,4) + streaks_gather(i)%data(iprime,4) .ne.&
                       -abs(streaks_gather(ip_c1)%data(iprime2,4)) - abs(streaks_gather(i)%data(iprime,4)))) cycle loop2

                   ! If streaks separated by a distance superior to Lz-dist_extr+dz, periodicity considered
                   if (((streaks_gather(ip_c1)%data(iprime2,3) - streaks_gather(i)%data(iprime,3))**2)**0.5.gt.dist_min_perio) then
                      if (streaks_gather(i)%data(iprime,3).le.zg(ngz/2)) then
                         z_decal = streaks_gather(i)%data(iprime,3) + Lz_fstt
                      else
                         z_decal = streaks_gather(i)%data(iprime,3) - Lz_fstt
                      endif
                      ! Calculation of the distance
                      dist_test = ((streaks_gather(ip_c1)%data(iprime2,2) - streaks_gather(i)%data(iprime,2))**2 + &
                                   (streaks_gather(ip_c1)%data(iprime2,3) - z_decal)**2)**0.5
                   else
                      ! Calculation of the distance
                      dist_test = ((streaks_gather(ip_c1)%data(iprime2,2) - streaks_gather(i)%data(iprime,2))**2 + &
                                   (streaks_gather(ip_c1)%data(iprime2,3) - streaks_gather(i)%data(iprime,3))**2)**0.5
                   endif
                   ! If old streak closer to new streak than the one already found, replaced
                   if (dist.gt.dist_test) then
                      dist = dist_test
                      ind_link = lp
                   endif

                enddo loop2

                ! To avoid associating streaks if too much distance
                if (dist.gt.dist_max_yz) then
                   to_be_tested(l) = 0
                   table_link(i)%link(l) = -1
                   cycle loop1
                ! If close enough to a streak not associated yet, linkage between current and neighbor streaks
                else if (assoc(ind_link,1).eq.-1) then
                   lp = ind_link
                   to_be_tested(l) = 0
                   table_link(i)%link(l) = lp
                   assoc(lp,1) = l
                   assoc(lp,2) = dist
                   streak_tip(ip_c1)%link(lp) = -1
                ! If neighbor streak already associated with a streak
                else
                   lp = ind_link
                   ! 1st situation: this streak is closer to the neighbor streak than the previous one
                   if (dist.lt.assoc(lp,2)) then
                      ! Indicate that the other streak needs to be try again for association
                      to_be_tested(assoc(lp,1)) = 1

                      ! Linkage between current and neighbor streaks
                      to_be_tested(l) = 0
                      table_link(i)%link(l) = lp
                      assoc(lp,1) = l
                      assoc(lp,2) = dist
                   ! 2nd situation: the previous streak is closer to the neighbor streak than this one
                   else
                      table_link(i)%link(l) = -1
                      to_be_tested(l) = 0
                   endif
                endif

             enddo loop1
             has_to_test = .false.
             do l=1,nstrk_per_proc(i*bl(1)%ndomk)
                if (to_be_tested(l).eq.1) has_to_test = .true.
             enddo
             enddo

             deallocate(to_be_tested,assoc)
          enddo loopproc

          ! 2/ Creation of the final streaks inventory
          ! ------------------------------------------
          ! Final streaks organisation before writting
          ! 2 steps performed here:
          !      a) Count of all the streaks long enough
          !         Elimination of streaks too short (< l_min_str)
          !    _ b) Association of full streak thanks to table_link
          !         and streak_tip

          ! a) Count of streaks long enough
          ! -------------------------------
          ! Initialisation of counter
          nstreaks_merge = 0
          ! Array to save the number of points for each full streak
          allocate(streaks_m_points(nstreaks_tot))
          ! Array to save the streaks retained
          allocate(assoc(0:bl(1)%ndomi-1,maxval(nstrk_per_proc)))
          assoc = 0
          do i=0,bl(1)%ndomi-1
             loop3: do l=1,nstrk_per_proc(i*bl(1)%ndomk)
                ! If not tip of streak, next
                if (streak_tip(i)%link(l).eq.-1) cycle loop3

                ! Initialisation of streak points counter
                npoints = 0

                has_to_test = .true.
                lp = l; ip = i*bl(1)%ndomk
                ip_c1 = ip/bl(1)%ndomk
                iend = sum(streaks_gather(i)%length(:l))
                ibeg = iend - streaks_gather(i)%length(l) + 1

                ! Determination of the length of the streak: xmin
                dist_tot = streaks_gather(ip_c1)%data(ibeg,1)
                ! Linkage until end of the streak
                loop4: do while (has_to_test)
                   ! Addition of points of current streak
                   npoints = npoints + streaks_gather(ip_c1)%length(lp)

                   ! Looking at following streak, or stopping if doesn't exist
                   if (table_link(ip_c1)%link(lp).eq.-1) then
                      has_to_test = .false.
                      cycle loop4
                   endif

                   ! Indices of the following streak
                   lp = table_link(ip_c1)%link(lp)
                   ip = neighbor_E_fstt(ip)
                   ip_c1 = ip/bl(1)%ndomk

                   ! Position of the following streak
                   iend = sum(streaks_gather(ip_c1)%length(:lp))
                   ibeg = iend - streaks_gather(ip_c1)%length(lp) + 1
                enddo loop4

                ! Determination of the length of the streak: abs(xmax - xmin)
                dist_tot = abs(streaks_gather(ip_c1)%data(iend,1)-dist_tot)

                ! Elimination of streaks with length inferior to l_min_str
                if (dist_tot.ge.l_min_str) then
                     nstreaks_merge = nstreaks_merge + 1
                     streaks_m_points(nstreaks_merge) = npoints
                     assoc(i,l) = 1
                endif
             enddo loop3
          enddo

          ! Allocation of streaks_merge
          allocate(streaks_merge(1:nstreaks_merge))

          ! b) Creation of full streak instance
          ! -----------------------------------
          ! Counter for streaks
          ls = 0
          do i=0,bl(1)%ndomi-1
             loop5: do l=1,nstrk_per_proc(i*bl(1)%ndomk)
                ! If not tip of streak, next
                if (assoc(i,l).eq.0) cycle loop5
                ! if (streak_tip(i)%link(l).eq.-1) cycle loop5

                ! Initialisation of streak points counter
                npoints = 0

                ! 1 added to counter
                ls = ls + 1

                ! Allocation to save streak data
                allocate(streaks_merge(ls)%data(1:streaks_m_points(ls),1:4))

                ! Initialisations
                has_to_test = .true.
                lp = l; ip = i*bl(1)%ndomk
                ip_c1 = ip/bl(1)%ndomk
                iend = sum(streaks_gather(i)%length(:l))
                ibeg = iend - streaks_gather(i)%length(l) + 1

                ! Linkage until end of the streak
                loop6: do while (has_to_test)
                   ! Addition of points of current streak
                   npoints = npoints + streaks_gather(ip_c1)%length(lp)

                   ! Data of current local streak saved
                   do j=1,4
                      iprime = npoints - streaks_gather(ip_c1)%length(lp) + 1
                      iprime2 = npoints
                      streaks_merge(ls)%data(iprime:iprime2,j) = streaks_gather(ip_c1)%data(ibeg:iend,j)
                   enddo

                   ! Looking at following streak, or stopping if doesn't exist
                   if (table_link(ip_c1)%link(lp).eq.-1) then
                      has_to_test = .false.
                      cycle loop6
                   endif

                   ! Indices of the following streak
                   lp = table_link(ip_c1)%link(lp)
                   ip = neighbor_E_fstt(ip)
                   ip_c1 = ip/bl(1)%ndomk

                   ! Position of the following streak
                   iend = sum(streaks_gather(ip_c1)%length(:lp))
                   ibeg = iend - streaks_gather(ip_c1)%length(lp) + 1

                enddo loop6
             enddo loop5
          enddo

          ! Deallocation
          do ip=0,bl(1)%ndomi-1
             deallocate(table_link(ip)%link)
             deallocate(streak_tip(ip)%link)
             deallocate(streaks_gather(ip)%length,streaks_gather(ip)%data)
          enddo
          deallocate(streaks_gather)
          deallocate(nstrk_per_proc,assoc)

          ! 3/ Association of streaks tip close to each other
          ! -------------------------------------------------
          ! Modification of the definition of dist_min_perio
          dist_min_perio = Lz_fstt - dist_max_yz

          ! a) Count of streaks that can be associated
          ! ------------------------------------------
          allocate(table_link_f(nstreaks_merge))
          allocate(assoc_f(nstreaks_merge))
          table_link_f = 0
          assoc_f = 0
          nstreaks_final = nstreaks_merge
          ! loop over the merged streaks
          loop7: do l=1,nstreaks_merge
             ! End of the streaks
             iend = streaks_m_points(l)
             x_end = streaks_merge(l)%data(iend,1)
             y_end = streaks_merge(l)%data(iend,2)
             z_end = streaks_merge(l)%data(iend,3)

             ! Initialisation
             dist = 1e6
             ind_link = -999

             ! Calculation of distance with each streaks tip
             ! Save of the closest
             loop8: do lp=1,nstreaks_merge
                ! Doesn't compare streak with itself
                if (l.eq.lp) cycle loop8

                ! Doesn't test if it is not the same type
                if ((streaks_merge(l)%data(iend,4)+streaks_merge(lp)%data(1,4).ne.&
                     abs(streaks_merge(l)%data(iend,4))+abs(streaks_merge(lp)%data(1,4))).and.&
                    (streaks_merge(l)%data(iend,4)+streaks_merge(lp)%data(1,4).ne.&
                    -abs(streaks_merge(l)%data(iend,4))-abs(streaks_merge(lp)%data(1,4)))) cycle loop8

                ! Calculation of 2 distances
                ! dist_test: distance btwn streaks in y-z plane
                ! dist_test2: streamwise distance btwn streaks tip
                ibeg = 1
                ! 1st way to do it: find the beginning indice for the other streak
                do while (streaks_merge(lp)%data(ibeg,1).le.streaks_merge(l)%data(streaks_m_points(l),1))
                   ibeg = ibeg + 1
                   if (ibeg.gt.streaks_m_points(lp)) cycle loop8
                enddo


                ! Definition of the coordinates to consider
                x_beg = streaks_merge(lp)%data(1,1)
                y_beg = streaks_merge(lp)%data(ibeg,2)
                z_beg = streaks_merge(lp)%data(ibeg,3)

                ! Calculation of dist_test2
                dist_test2 = ((x_end - x_beg)**2)**0.5
                ! If too much distance between streaks tip, cycle
                if (dist_test2.gt.dist_max_x) cycle loop8

                ! If streaks separated by a distance superior to Lz-dist_extr+dz, periodicity considered
                if (((z_beg - z_end)**2)**0.5.gt.dist_min_perio) then
                   if (z_beg.le.zg(ngz/2)) then
                      z_beg = z_beg + Lz_fstt + z(2)-z(1)
                   else
                      z_beg = z_beg - Lz_fstt - (z(2)-z(1))
                   endif
                endif

                ! Calculation of dist_test
                dist_test = ((y_end - y_beg)**2 + &
                             (z_end - z_beg)**2)**0.5
                if (dist.gt.dist_test) then
                   dist = dist_test
                   ind_link = lp
                endif
             enddo loop8


             ! If distance between tips too important, cycle
             if (dist.gt.dist_max_yz) then
                table_link_f(l) = -1
                cycle loop7
             ! Save of linkage between streaks, if no linkage made yet
             else if (assoc_f(ind_link).eq.0) then
                nstreaks_final = nstreaks_final - 1
                table_link_f(l) = ind_link
                assoc_f(ind_link) = l
             ! Save of linkage between streaks
             ! If already a linkage, the longest streak is retained
             else
                ! if previous streak longer
                if (streaks_m_points(assoc_f(ind_link)).gt.streaks_m_points(l)) then
                   ! current streak
                   table_link_f(l) = -1
                else
                   ! old streak
                   table_link_f(assoc_f(ind_link)) = -1
                   ! current streak
                   table_link_f(l) = ind_link
                   assoc_f(ind_link) = l
                endif
             endif
          enddo loop7

          ! Allocation of streaks_final
          allocate(streaks_final(1:nstreaks_final))
          allocate(streaks_f_points(nstreaks_final))
          streaks_f_points = 0

          ! b) Realisation of the final association
          ! ---------------------------------------
          ! Counter for streaks
          ls = 0
          loop9: do l=1,nstreaks_merge
             ! if associated, mean it is not the tip
             if (assoc_f(l).ne.0) cycle loop9

             ! 1 added to counter
             ls = ls + 1

             ! (i) Calculation of the length of the streak
             has_to_test = .true.
             lp = l
             lp_save = lp
             do while (has_to_test)
                ! Determination where to begin to add the next streak
                ! Comparison between x_end of the final streak and x_beg of the current streak
                ! (x_beg must be > to x_end)
                ibeg = 1
                if (lp.ne.l) then
                   do while (streaks_merge(lp)%data(ibeg,1).le.streaks_merge(lp_save)%data(streaks_m_points(lp_save),1))
                      ibeg = ibeg + 1
                   enddo
                endif
                streaks_f_points(ls) = streaks_f_points(ls) + (streaks_m_points(lp) - ibeg + 1)

                ! Next streak to associate
                lp_save = lp
                lp = table_link_f(lp)

                ! Stop when no association found after
                if (lp.eq.-1) has_to_test = .false.
             enddo

             ! Allocation to save streak data
             allocate(streaks_final(ls)%data(1:streaks_f_points(ls),1:4))

             ! (ii) Realisation of the final list
             has_to_test = .true.
             lp = l
             lp_save = lp
             npoints = 0
             do while (has_to_test)
                ! Determination where to begin to add the next streak
                ! Comparison between x_end of the final streak and x_beg of the current streak
                ! (x_beg must be > to x_end)
                ibeg = 1
                iend = streaks_m_points(lp)
                if (lp.ne.l) then
                   do while (streaks_merge(lp)%data(ibeg,1).le.streaks_merge(lp_save)%data(streaks_m_points(lp_save),1))
                      ibeg = ibeg + 1
                   enddo
                endif
                npoints = npoints + (streaks_m_points(lp) - ibeg + 1)

                ! Data of current local streak saved
                do j=1,4
                   iprime = npoints - (streaks_m_points(lp) - ibeg)
                   iprime2 = npoints
                   streaks_final(ls)%data(iprime:iprime2,j) = streaks_merge(lp)%data(ibeg:iend,j)
                enddo

                ! Next streak to associate
                lp_save = lp
                lp = table_link_f(lp)

                ! Stop when no association found after
                if (lp.eq.-1) has_to_test = .false.
             enddo
          enddo loop9

          ! Deallocation
          deallocate(streaks_merge,streaks_m_points,assoc_f,table_link_f)

          ! c) Enforce periodicity for the streak
          ! -------------------------------------
          do ls=1,nstreaks_final
             switch = .false.
             z_decal = streaks_final(ls)%data(1,3)
             if (streaks_final(ls)%data(1,3).ge.zg(ngz/2)) then
                do i=2,streaks_f_points(ls)
                   if (abs(streaks_final(ls)%data(i,3)-z_decal).gt.Lz_fstt/2) then
                      if (switch) then
                         switch = .false.
                      else
                         switch = .true.
                      endif
                   endif
                   z_decal = streaks_final(ls)%data(i,3)
                   if (switch) streaks_final(ls)%data(i,3) = streaks_final(ls)%data(i,3) + Lz_fstt
                enddo
             else
                do i=2,streaks_f_points(ls)
                   if (abs(streaks_final(ls)%data(i,3)-z_decal).gt.Lz_fstt/2) then
                      if (switch) then
                         switch = .false.
                      else
                         switch = .true.
                      endif
                   endif
                   z_decal = streaks_final(ls)%data(i,3)
                   if (switch) streaks_final(ls)%data(i,3) = streaks_final(ls)%data(i,3) - Lz_fstt
                enddo
             endif
          enddo

          ! 4/ Writting of streak file
          ! --------------------------
          open(180,file='streaks_bl'//trim(numchar(nob(iproc)))//'.bin',form='unformatted',status='unknown',position="append")
          write(180) nstreaks_final
          ! Loop over all streaks
          do l=1,nstreaks_final
             write(180) streaks_f_points(l)
             write(180) ((streaks_final(l)%data(i,j),i=1,streaks_f_points(l)),j=1,4)
          enddo
          close(180)

          ! Deallocation
          deallocate(streaks_final,streaks_f_points)

       endif
    endif

    ! call mpistop('',0)

  end subroutine merging_streaks

end module mod_pp_fstt_streaks
