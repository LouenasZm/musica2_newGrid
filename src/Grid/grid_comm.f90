!===============================================================================
subroutine grid_comm
!===============================================================================
  !> Communicate global grids in interblock communicator
  !> and then gather in intrablock communicator
!===============================================================================
  use mod_grid
  use mod_grid_comm
  use mod_constant ! <~ to know L_ref TO BE CHANGED
  use warnstop
  use mod_mpi_part
  use mod_bc_periodicity ! <~ needed for Lxp,Lyp,Lzp  TO BE CHANGED
  implicit none
  ! ----------------------------------------------------------------------------
  integer :: i,j
  ! ----------------------------------------------------------------------------
  
  ! Definition of MPI types for grid communications
  ! ===============================================
  if (is_curv3) then
     !call grid_comm_type_curv3(ngx,ngy,ngz)
     call grid_comm_type_curv_2(ngx,ngy,ngz)
  else
     call grid_comm_type_cart
     !if (is_curv) call grid_comm_type_curv(ngx,ngy)
     if (is_curv) call grid_comm_type_curv_2(ngx,ngy,ngz)
     call MPI_BARRIER(COMM_global,info)
  endif

  ! Interblock communications of global grid
  ! ========================================
  nbl=nob(iproc)

  if (iproc.eq.iproc_leader(nbl)) then
          
     if (is_curv3) then
!!$        call grid_comm_interblock_curv3(xgc3,ngx,ngy,ngz)
!!$        call grid_comm_interblock_curv3(ygc3,ngx,ngy,ngz)
!!$        call grid_comm_interblock_curv3(zgc3,ngx,ngy,ngz)
!!$        ! second pass for corners
!!$        call grid_comm_interblock_curv3(xgc3,ngx,ngy,ngz)
!!$        call grid_comm_interblock_curv3(ygc3,ngx,ngy,ngz)
!!$        call grid_comm_interblock_curv3(zgc3,ngx,ngy,ngz)

        call grid_comm_interblock_curv3_2(xgc3,ngx,ngy,ngz)
        call grid_comm_interblock_curv3_2(ygc3,ngx,ngy,ngz)
        call grid_comm_interblock_curv3_2(zgc3,ngx,ngy,ngz)
        call grid_comm_interblock_curv3_2(xgc3,ngx,ngy,ngz)
        call grid_comm_interblock_curv3_2(ygc3,ngx,ngy,ngz)
        call grid_comm_interblock_curv3_2(zgc3,ngx,ngy,ngz)
     else
        if (is_curv) then
           !call grid_comm_interblock_curv(xgc,ngx,ngy)
           !call grid_comm_interblock_curv(ygc,ngx,ngy)
           call grid_comm_interblock_curv_2(xgc,ngx,ngy)
           call grid_comm_interblock_curv_2(ygc,ngx,ngy)
        else
           call grid_comm_interblock_cart(xg,ngx,1)
           call grid_comm_interblock_cart(yg,ngy,2)
        endif

        if (.not.(is_2D)) then
           call grid_comm_interblock_cart(zg,ngz,3)
        endif
     endif

!if (iproc==0) print *,xgc3(:,25,1)
!if (iproc==0) print *,'a envoyer'
!if (iproc==1) print *,xgc3(25,:,0)
!if (iproc==1) print *,'-------'
!if (iproc==1) print *,xgc3(25,:,51)

     !call mpistop('rrr',0)

     ! BC imin
     if (bl(nbl)%BC(1)>0) then
        if ((Lxp.ne.0.0_wp).and.(abs(xgc(0,1)-xgc(1,1)).gt.0.9_wp*Lxp)) then
           print *,iproc,nob(iproc),(xgc(0,1)-xgc(1,1)),'pb imin'
           if ((xgc(0,1)-xgc(1,1)).gt.0.0_wp) then
              xgc(1-ngh:0,:)=xgc(1-ngh:0,:)-Lxp
           else
              xgc(1-ngh:0,:)=xgc(1-ngh:0,:)+Lxp
           endif
        endif
        if ((Lyp.ne.0.0_wp).and.(abs(ygc(0,1)-ygc(1,1)).gt.0.9_wp*Lyp)) then
           print *,iproc,nob(iproc),(ygc(0,1)-ygc(1,1)),'pb imin'
           if ((ygc(0,1)-ygc(1,1)).gt.0.0_wp) then
              ygc(1-ngh:0,:)=ygc(1-ngh:0,:)-Lyp
           else
              ygc(1-ngh:0,:)=ygc(1-ngh:0,:)+Lyp
           endif
        endif
     endif
     ! BC imax
     if (bl(nbl)%BC(2)>0) then
        if ((Lxp.ne.0.0_wp).and.(abs(xgc(ngx+1,1)-xgc(ngx,1)).gt.0.9_wp*Lxp)) then
           print *,iproc,nob(iproc),(xgc(ngx+1,1)-xgc(ngx,1)),'pb imax'
           if ((xgc(ngx+1,1)-xgc(ngx,1)).gt.0.0_wp) then
              xgc(ngx+1:ngx+ngh,:)=xgc(ngx+1:ngx+ngh,:)-Lxp
           else
              xgc(ngx+1:ngx+ngh,:)=xgc(ngx+1:ngx+ngh,:)+Lxp
           endif
        endif
        if ((Lyp.ne.0.0_wp).and.(abs(ygc(ngx+1,1)-ygc(ngx,1)).gt.0.9_wp*Lyp)) then
           print *,iproc,nob(iproc),(ygc(ngx+1,1)-ygc(ngx,1)),'pb imax'
           if ((ygc(ngx+1,1)-ygc(ngx,1)).gt.0.0_wp) then
              ygc(ngx+1:ngx+ngh,:)=ygc(ngx+1:ngx+ngh,:)-Lyp
           else
              ygc(ngx+1:ngx+ngh,:)=ygc(ngx+1:ngx+ngh,:)+Lyp
           endif
        endif
     endif
     ! BC jmin
     if (bl(nbl)%BC(3)>0) then
        if ((Lyp.ne.0.0_wp).and.(abs(ygc(1,0)-ygc(1,1)).gt.0.9_wp*Lyp)) then
           print *,iproc,nob(iproc),(ygc(1,0)-ygc(1,1)),'pb jmin'
           if ((ygc(1,0)-ygc(1,1)).gt.0.0_wp) then
              ygc(:,1-ngh:0)=ygc(:,1-ngh:0)-Lyp
           else
              ygc(:,1-ngh:0)=ygc(:,1-ngh:0)+Lyp
           endif
        endif
     endif
     ! BC jmax
     if (bl(nbl)%BC(4)>0) then
        if ((Lyp.ne.0.0_wp).and.(abs(ygc(1,ngy+1)-ygc(1,ngy)).gt.0.9_wp*Lyp)) then
           print *,iproc,nob(iproc),(ygc(1,ngy+1)-ygc(1,ngy)),'pb jmax'
           if ((ygc(1,ngy+1)-ygc(1,ngy)).gt.0.0_wp) then
              ygc(:,ngy+1:ngy+ngh)=ygc(:,ngy+1:ngy+ngh)-Lyp
           else
              ygc(:,ngy+1:ngy+ngh)=ygc(:,ngy+1:ngy+ngh)+Lyp
           endif
        endif
     endif

     if (.not.is_curv3) then
     !if (iproc==4) then
     if ((.not.CYL).and.(.not.ACT).and.(.not.TURB).and.(.not.TE)) then
        ! correction for periodicity
        ! --------------------------
        ! BC imin
        if (bl(nbl)%BC(1)>0) then
           if (is_curv3) then
              if (iproc==1) print *,xgc3(0,ngy,1),xgc3(1,ngy,1),xgc3(2,ngy,1)
              if (iproc==1) print *,zgc3(0,ngy,1),zgc3(1,ngy,1),zgc3(2,ngy,1)
              if ((xgc3(0,ngy,1)>xgc3(1,ngy,1)).and.(xgc3(2,ngy,1)>xgc3(1,ngy,1))) then
                 print *,iproc,'period (1,1)'
                 if (Lxp==0.0_wp) Lxp=xmax-xmin+abs(xgc3(2,ngy,1)-xgc3(1,ngy,1))
                 print *,iproc,'(1,1)',Lxp,Lxp/L_ref
                 xgc3(1-ngh:0,:,:)=xgc3(1-ngh:0,:,:)-Lxp
              endif
              if ((xgc3(0,ngy,1)<xgc3(1,ngy,1)).and.(xgc3(2,ngy,1)<xgc3(1,ngy,1))) then
                 print *,iproc,'period (1,1)'
                 if (Lxp==0.0_wp) Lxp=xmax-xmin+abs(xgc3(2,ngy,1)-xgc3(1,ngy,1))
                 print *,iproc,'(1,1)',Lxp,Lxp/L_ref
                 xgc3(1-ngh:0,:,:)=xgc3(1-ngh:0,:,:)+Lxp
              endif

!!$              if ((zgc3(0,ngy,1)>zgc3(1,ngy,1)).and.(zgc3(2,ngy,1)>zgc3(1,ngy,1))) then
!!$                 print *,iproc,'period (1,1)'
!!$                 if (Lzp==0.0_wp) Lzp=zmax-zmin+abs(zgc3(2,ngy,1)-zgc3(1,ngy,1))
!!$                 print *,iproc,'(1,1)',Lzp,Lzp/L_ref
!!$                 zgc3(1-ngh:0,:,:)=zgc3(1-ngh:0,:,:)-Lzp
!!$              endif
!!$              if ((zgc3(0,ngy,1)<zgc3(1,ngy,1)).and.(zgc3(2,ngy,1)<zgc3(1,ngy,1))) then
!!$                 print *,iproc,'period (1,1)'
!!$                 if (Lzp==0.0_wp) Lzp=zmax-zmin+abs(zgc3(2,ngy,1)-zgc3(1,ngy,1))
!!$                 print *,iproc,'(1,1)',Lzp,Lzp/L_ref
!!$                 zgc3(1-ngh:0,:,:)=zgc3(1-ngh:0,:,:)+Lzp
!!$              endif
           else
              if (is_curv) then
                 if ((xgc(0,ngy)>xgc(1,ngy)).and.(xgc(2,ngy)>xgc(1,ngy))) then
                    print *,iproc,'period (1,1)'
                    if (Lxp==0.0_wp) Lxp=xmax-xmin+abs(xgc(2,ngy)-xgc(1,ngy))
                    print *,iproc,'(1,1)',Lxp,Lxp/L_ref
                    xgc(1-ngh:0,:)=xgc(1-ngh:0,:)-Lxp
                 endif
                 if ((xgc(0,ngy)<xgc(1,ngy)).and.(xgc(2,ngy)<xgc(1,ngy))) then
                    print *,iproc,'period (1,1)'
                    if (Lxp==0.0_wp) Lxp=xmax-xmin+abs(xgc(2,ngy)-xgc(1,ngy))
                    print *,iproc,'(1,1)',Lxp,Lxp/L_ref
                    xgc(1-ngh:0,:)=xgc(1-ngh:0,:)+Lxp
                 endif
              else
                 if ((xg(0)>xg(1)).and.(xg(2)>xg(1))) then
                    print *,iproc,'period (1,1) cart'
                    if (Lxp==0.0_wp) Lxp=xmax-xmin+abs(xg(2)-xg(1))
                    print *,iproc,'(1,1)',Lxp,Lxp/L_ref
                    xg(1-ngh:0)=xg(1-ngh:0)-Lxp
                 endif
                 if ((xg(0)<xg(1)).and.(xg(2)<xg(1))) then
                    print *,iproc,'period (1,1) cart'
                    if (Lxp==0.0_wp) Lxp=xmax-xmin+abs(xg(2)-xg(1))
                    print *,iproc,'(1,1)',Lxp,Lxp/L_ref
                    xg(1-ngh:0)=xg(1-ngh:0)+Lxp
                 endif
              endif
           endif
        endif
        ! BC imax
        if (bl(nbl)%BC(2)>0) then
           if (is_curv3) then
              if (iproc==1) print *,xgc3(ngx+1,ngy,1),xgc3(ngx,ngy,1),xgc3(ngx-1,ngy,1)
              if (iproc==1) print *,zgc3(ngx+1,ngy,1),zgc3(ngx,ngy,1),zgc3(ngx-1,ngy,1)
              if ((xgc3(ngx+1,ngy,1)<xgc3(ngx,ngy,1)).and.(xgc3(ngx-1,ngy,1)<xgc3(ngx,ngy,1))) then
                 print *,iproc,'period (1,2)'
                 if (Lxp==0.0_wp) Lxp=xmax-xmin+abs(xgc3(ngx,ngy,1)-xgc3(ngx-1,ngy,1))
                 print *,iproc,'(1,2)',Lxp,Lxp/L_ref
                 xgc3(ngx+1:ngx+ngh,:,:)=xgc3(ngx+1:ngx+ngh,:,:)+Lxp
              endif
              if ((xgc3(ngx+1,ngy,1)>xgc3(ngx,ngy,1)).and.(xgc3(ngx-1,ngy,1)>xgc3(ngx,ngy,1))) then
                 print *,iproc,'period (1,2)'
                 if (Lxp==0.0_wp) Lxp=xmax-xmin+abs(xgc3(ngx,ngy,1)-xgc3(ngx-1,ngy,1))
                 print *,iproc,'(1,2)',Lxp,Lxp/L_ref
                 xgc3(ngx+1:ngx+ngh,:,:)=xgc3(ngx+1:ngx+ngh,:,:)-Lxp
              endif

!!$              if ((zgc3(ngx+1,ngy,1)<zgc3(ngx,ngy,1)).and.(zgc3(ngx-1,ngy,1)<zgc3(ngx,ngy,1))) then
!!$                 print *,iproc,'period (1,2)'
!!$                 if (Lzp==0.0_wp) Lzp=zmax-zmin+abs(zgc3(ngx,ngy,1)-zgc3(ngx-1,ngy,1))
!!$                 print *,iproc,'(1,2)',Lzp,Lzp/L_ref
!!$                 zgc3(ngx+1:ngx+ngh,:,:)=zgc3(ngx+1:ngx+ngh,:,:)+Lzp
!!$              endif
!!$              if ((zgc3(ngx+1,ngy,1)>zgc3(ngx,ngy,1)).and.(zgc3(ngx-1,ngy,1)>zgc3(ngx,ngy,1))) then
!!$                 print *,iproc,'period (1,2)'
!!$                 if (Lzp==0.0_wp) Lzp=zmax-zmin+abs(zgc3(ngx,ngy,1)-zgc3(ngx-1,ngy,1))
!!$                 print *,iproc,'(1,2)',Lzp,Lzp/L_ref
!!$                 zgc3(ngx+1:ngx+ngh,:,:)=zgc3(ngx+1:ngx+ngh,:,:)-Lzp
!!$              endif
           else
              if (is_curv) then
                 if ((xgc(ngx+1,ngy)<xgc(ngx,ngy)).and.(xgc(ngx-1,ngy)<xgc(ngx,ngy))) then
                    print *,iproc,'period (1,2)'
                    if (Lxp==0.0_wp) Lxp=xmax-xmin+abs(xgc(ngx,ngy)-xgc(ngx-1,ngy))
                    print *,iproc,'(1,2)',Lxp,Lxp/L_ref
                    xgc(ngx+1:ngx+ngh,:)=xgc(ngx+1:ngx+ngh,:)+Lxp
                 endif
                 if ((xgc(ngx+1,ngy)>xgc(ngx,ngy)).and.(xgc(ngx-1,ngy)>xgc(ngx,ngy))) then
                    print *,iproc,'period (1,2)'
                    if (Lxp==0.0_wp) Lxp=xmax-xmin+abs(xgc(ngx,ngy)-xgc(ngx-1,ngy))
                    print *,iproc,'(1,2)',Lxp,Lxp/L_ref
                    xgc(ngx+1:ngx+ngh,:)=xgc(ngx+1:ngx+ngh,:)-Lxp
                 endif
              else
                 if ((xg(ngx+1)<xg(ngx)).and.(xg(ngx-1)<xg(ngx))) then
                    print *,iproc,'period (1,2) cart'
                    if (Lxp==0.0_wp) Lxp=xmax-xmin+abs(xg(ngx)-xg(ngx-1))
                    print *,iproc,'(1,2)',Lxp,Lxp/L_ref
                    xg(ngx+1:ngx+ngh)=xg(ngx+1:ngx+ngh)+Lxp
                 endif
                 if ((xg(ngx+1)>xg(ngx)).and.(xg(ngx-1)>xg(ngx))) then
                    print *,iproc,'period (1,2) cart'
                    if (Lxp==0.0_wp) Lxp=xmax-xmin+abs(xg(ngx)-xg(ngx-1))
                    print *,iproc,'(1,2)',Lxp,Lxp/L_ref
                    xg(ngx+1:ngx+ngh)=xg(ngx+1:ngx+ngh)-Lxp
                 endif
              endif
           endif
        endif
        ! BC jmin
        if (bl(nbl)%BC(3)>0) then
           if (is_curv3) then
              if ((ygc3(ngx,0,1)>ygc3(ngx,1,1)).and.(ygc3(ngx,2,1)>ygc3(ngx,1,1))) then
                 print *,iproc,'period (2,1)'
                 if (Lyp==0.0_wp) Lyp=ymax-ymin+abs(ygc3(ngx,2,1)-ygc3(ngx,1,1))
                 print *,iproc,'(2,1)',Lyp,Lyp/L_ref
                 ygc3(:,1-ngh:0,:)=ygc3(:,1-ngh:0,:)-Lyp
              endif
              if ((ygc3(ngx,0,1)<ygc3(ngx,1,1)).and.(ygc3(ngx,2,1)<ygc3(ngx,1,1))) then
                 print *,iproc,'period (2,1)'
                 if (Lyp==0.0_wp) Lyp=ymax-ymin+abs(ygc3(ngx,2,1)-ygc3(ngx,1,1))
                 print *,iproc,'(2,1)',Lyp,Lyp/L_ref
                 ygc3(:,1-ngh:0,:)=ygc3(:,1-ngh:0,:)+Lyp
              endif
!!$
!!$              if ((zgc3(ngx,0,1)>zgc3(ngx,1,1)).and.(zgc3(ngx,2,1)>zgc3(ngx,1,1))) then
!!$                 print *,iproc,'period (2,1)'
!!$                 if (Lzp==0.0_wp) Lzp=zmax-zmin+abs(zgc3(ngx,2,1)-zgc3(ngx,1,1))
!!$                 print *,iproc,'(2,1)',Lzp,Lzp/L_ref
!!$                 zgc3(:,1-ngh:0,:)=zgc3(:,1-ngh:0,:)-Lzp
!!$              endif
!!$              if ((zgc3(ngx,0,1)<zgc3(ngx,1,1)).and.(zgc3(ngx,2,1)<zgc3(ngx,1,1))) then
!!$                 print *,iproc,'period (2,1)'
!!$                 if (Lzp==0.0_wp) Lzp=zmax-zmin+abs(zgc3(ngx,2,1)-zgc3(ngx,1,1))
!!$                 print *,iproc,'(2,1)',Lzp,Lzp/L_ref
!!$                 zgc3(:,1-ngh:0,:)=zgc3(:,1-ngh:0,:)+Lzp
!!$              endif
           else
              if (is_curv) then
                 if ((ygc(ngx,0)>ygc(ngx,1)).and.(ygc(ngx,2)>ygc(ngx,1))) then
                    print *,iproc,'period (2,1)'
                    if (Lyp==0.0_wp) Lyp=ymax-ymin+abs(ygc(ngx,2)-ygc(ngx,1))
                    print *,iproc,'(2,1)',Lyp,Lyp/L_ref
                    ygc(:,1-ngh:0)=ygc(:,1-ngh:0)-Lyp
                 endif
                 if ((ygc(ngx,0)<ygc(ngx,1)).and.(ygc(ngx,2)<ygc(ngx,1))) then
                    print *,iproc,'period (2,1)'
                    if (Lyp==0.0_wp) Lyp=ymax-ymin+abs(ygc(ngx,2)-ygc(ngx,1))
                    print *,iproc,'(2,1)',Lyp,Lyp/L_ref
                    ygc(:,1-ngh:0)=ygc(:,1-ngh:0)+Lyp
                 endif
              else
                 if ((yg(0)>yg(1)).and.(yg(2)>yg(1))) then
                    print *,iproc,'period (2,1) cart'
                    if (Lyp==0.0_wp) Lyp=ymax-ymin+abs(yg(2)-yg(1))
                    print *,iproc,'(2,1)',Lyp,Lyp/L_ref
                    yg(1-ngh:0)=yg(1-ngh:0)-Lyp
                 endif
                 if ((yg(0)<yg(1)).and.(yg(2)<yg(1))) then
                    print *,iproc,'period (2,1) cart'
                    if (Lyp==0.0_wp) Lyp=ymax-ymin+abs(yg(2)-yg(1))
                    print *,iproc,'(2,1)',Lyp,Lyp/L_ref
                    yg(1-ngh:0)=yg(1-ngh:0)+Lyp
                 endif
              endif
           endif
        endif
        ! BC jmax
        if (bl(nbl)%BC(4)>0) then
           if (is_curv3) then
              if ((ygc3(ngx,ngy+1,1)<ygc3(ngx,ngy,1)).and.(ygc3(ngx,ngy-1,1)<ygc3(ngx,ngy,1))) then
                 print *,iproc,'period (2,2)'
                 if (Lyp==0.0_wp) Lyp=ymax-ymin+abs(ygc3(ngx,ngy,1)-ygc3(ngx,ngy-1,1))
                 print *,iproc,'(2,2)',Lyp,Lyp/L_ref
                 ygc3(:,ngy+1:ngy+ngh,:)=ygc3(:,ngy+1:ngy+ngh,:)+Lyp
              endif
              if ((ygc3(ngx,ngy+1,1)>ygc3(ngx,ngy,1)).and.(ygc3(ngx,ngy-1,1)>ygc3(ngx,ngy,1))) then
                 print *,iproc,'period (2,2)'
                 if (Lyp==0.0_wp) Lyp=ymax-ymin+abs(ygc3(ngx,ngy,1)-ygc3(ngx,ngy-1,1))
                 print *,iproc,'(2,2)',Lyp,Lyp/L_ref
                 ygc3(:,ngy+1:ngy+ngh,:)=ygc3(:,ngy+1:ngy+ngh,:)-Lyp
              endif
!!$
!!$              if ((zgc3(ngx,ngy+1,1)<zgc3(ngx,ngy,1)).and.(zgc3(ngx,ngy-1,1)<zgc3(ngx,ngy,1))) then
!!$                 print *,iproc,'period (2,2)'
!!$                 if (Lzp==0.0_wp) Lzp=zmax-zmin+abs(zgc3(ngx,ngy,1)-zgc3(ngx,ngy-1,1))
!!$                 print *,iproc,'(2,2)',Lzp,Lzp/L_ref
!!$                 zgc3(:,ngy+1:ngy+ngh,:)=zgc3(:,ngy+1:ngy+ngh,:)+Lzp
!!$              endif
!!$              if ((zgc3(ngx,ngy+1,1)>zgc3(ngx,ngy,1)).and.(zgc3(ngx,ngy-1,1)>zgc3(ngx,ngy,1))) then
!!$                 print *,iproc,'period (2,2)'
!!$                 if (Lzp==0.0_wp) Lzp=zmax-zmin+abs(zgc3(ngx,ngy,1)-zgc3(ngx,ngy-1,1))
!!$                 print *,iproc,'(2,2)',Lzp,Lzp/L_ref
!!$                 zgc3(:,ngy+1:ngy+ngh,:)=zgc3(:,ngy+1:ngy+ngh,:)-Lzp
!!$              endif
           else
              if (is_curv) then
                 if ((ygc(ngx,ngy+1)<ygc(ngx,ngy)).and.(ygc(ngx,ngy-1)<ygc(ngx,ngy))) then
                    print *,iproc,'period (2,2)'
                    if (Lyp==0.0_wp) Lyp=ymax-ymin+abs(ygc(ngx,ngy)-ygc(ngx,ngy-1))
                    print *,iproc,'(2,2)',Lyp,Lyp/L_ref
                    ygc(:,ngy+1:ngy+ngh)=ygc(:,ngy+1:ngy+ngh)+Lyp
                 endif
                 if ((ygc(ngx,ngy+1)>ygc(ngx,ngy)).and.(ygc(ngx,ngy-1)>ygc(ngx,ngy))) then
                    print *,iproc,'period (2,2)'
                    if (Lyp==0.0_wp) Lyp=ymax-ymin+abs(ygc(ngx,ngy)-ygc(ngx,ngy-1))
                    print *,iproc,'(2,2)',Lyp,Lyp/L_ref
                    ygc(:,ngy+1:ngy+ngh)=ygc(:,ngy+1:ngy+ngh)-Lyp
                 endif
              else
                 if ((yg(ngy+1)<yg(ngy)).and.(yg(ngy-1)<yg(ngy))) then
                    print *,iproc,'period (2,2) cart'
                    if (Lyp==0.0_wp) Lyp=ymax-ymin+abs(yg(ngy)-yg(ngy-1))
                    print *,iproc,'(2,2)',Lyp,Lyp/L_ref
                    yg(ngy+1:ngy+ngh)=yg(ngy+1:ngy+ngh)+Lyp
                 endif
                 if ((yg(ngy+1)>yg(ngy)).and.(yg(ngy-1)>yg(ngy))) then
                    print *,iproc,'period (2,2) cart'
                    if (Lyp==0.0_wp) Lyp=ymax-ymin+abs(yg(ngy)-yg(ngy-1))
                    print *,iproc,'(2,2)',Lyp,Lyp/L_ref
                    yg(ngy+1:ngy+ngh)=yg(ngy+1:ngy+ngh)-Lyp
                 endif
              endif
           endif
        endif
     endif
     endif

!!$     call mpistop('rr',0)

     ! BC kmin
     !if ((.not.is_2D).and.(.not.is_curv3)) then
     if (.not.is_2D) then     
        if (bl(nbl)%BC(5)>0) then
           if (is_curv3) then
              if ((zgc3(1,1,0)>zgc3(1,1,1)).and.(zgc3(1,1,2)>zgc3(1,1,1))) then
                 print *,iproc,'period (3,1)'
                 if (Lzp==0.0_wp) Lzp=zmax-zmin+abs(zgc3(1,1,2)-zgc3(1,1,1))
                 print *,iproc,'(3,1)',Lzp,Lzp/L_ref
                 zgc3(:,:,1-ngh:0)=zgc3(:,:,1-ngh:0)-Lzp
              endif
              if ((zgc3(1,1,0)<zgc3(1,1,1)).and.(zgc3(1,1,2)<zgc3(1,1,1))) then
                 print *,iproc,'period (3,1)'
                 if (Lzp==0.0_wp) Lzp=zmax-zmin+abs(zgc3(1,1,2)-zgc3(1,1,1))
                 print *,iproc,'(3,1)',Lzp,Lzp/L_ref
                 zgc3(:,:,1-ngh:0)=zgc3(:,:,1-ngh:0)+Lzp
              endif
           else
              if ((zg(0)>zg(1)).and.(zg(2)>zg(1))) then
                 print *,iproc,'period (3,1)'
                 if (Lzp==0.0_wp) Lzp=zmax-zmin+abs(zg(2)-zg(1))
                 print *,iproc,'(3,1)',Lzp,Lzp/L_ref
                 zg(1-ngh:0)=zg(1-ngh:0)-Lzp
              endif
              if ((zg(0)<zg(1)).and.(zg(2)<zg(1))) then
                 print *,iproc,'period (3,1)'
                 if (Lzp==0.0_wp) Lzp=zmax-zmin+abs(zg(2)-zg(1))
                 print *,iproc,'(3,1)',Lzp,Lzp/L_ref
                 zg(1-ngh:0)=zg(1-ngh:0)+Lzp
              endif
           endif
        endif
        ! BC kmax
        if (bl(nbl)%BC(6)>0) then     
           if (is_curv3) then
              if ((zgc3(1,1,ngz+1)<zgc3(1,1,ngz)).and.(zgc3(1,1,ngz-1)<zgc3(1,1,ngz))) then
                 print *,iproc,'period (3,2)'
                 if (Lzp==0.0_wp) Lzp=zmax-zmin+abs(zgc3(1,1,ngz)-zgc3(1,1,ngz-1))
                 print *,iproc,'(3,2)',Lzp,Lzp/L_ref
                 zgc3(:,:,ngz+1:ngz+ngh)=zgc3(:,:,ngz+1:ngz+ngh)+Lzp
              endif
              if ((zgc3(1,1,ngz+1)>zgc3(1,1,ngz)).and.(zgc3(1,1,ngz-1)>zgc3(1,1,ngz))) then
                 print *,iproc,'period (3,2)'
                 if (Lzp==0.0_wp) Lzp=zmax-zmin+abs(zgc3(1,1,ngz)-zgc3(1,1,ngz-1))
                 print *,iproc,'(3,2)',Lzp,Lzp/L_ref
                 zgc3(:,:,ngz+1:ngz+ngh)=zgc3(:,:,ngz+1:ngz+ngh)-Lzp
              endif
           else
              if ((zg(ngz+1)<zg(ngz)).and.(zg(ngz-1)<zg(ngz))) then
                 print *,iproc,'period (3,2)'
                 if (Lzp==0.0_wp) Lzp=zmax-zmin+abs(zg(ngz)-zg(ngz-1))
                 print *,iproc,'(3,2)',Lzp,Lzp/L_ref
                 zg(ngz+1:ngz+ngh)=zg(ngz+1:ngz+ngh)+Lzp
              endif
              if ((zg(ngz+1)>zg(ngz)).and.(zg(ngz-1)>zg(ngz))) then
                 print *,iproc,'period (3,2)'
                 if (Lzp==0.0_wp) Lzp=zmax-zmin+abs(zg(ngz)-zg(ngz-1))
                 print *,iproc,'(3,2)',Lzp,Lzp/L_ref
                 zg(ngz+1:ngz+ngh)=zg(ngz+1:ngz+ngh)-Lzp
              endif
           endif
        endif
     endif

  endif
  
  call MPI_BARRIER(COMM_global,info)

  ! Regeneration of fictitious corners !! ATTENTION TEMPORAIRE
  ! ==================================
  if (is_curv) then
     if (iproc.eq.iproc_leader(nob(iproc))) then
        if (dble_corner(nNW)) then
           ! NW
           do j=ngy+1,ngy+ngh
              do i=0,1-ngh,-1
                 xgc(i,j)=-xgc(i+1,j-1)+xgc(i+1,j)+xgc(i,j-1)
                 ygc(i,j)=-ygc(i+1,j-1)+ygc(i+1,j)+ygc(i,j-1)
              enddo
           enddo
        endif
        if (dble_corner(nNE)) then
           ! NE
           do j=ngy+1,ngy+ngh
              do i=ngx+1,ngx+ngh
                 xgc(i,j)=-xgc(i-1,j-1)+xgc(i-1,j)+xgc(i,j-1)
                 ygc(i,j)=-ygc(i-1,j-1)+ygc(i-1,j)+ygc(i,j-1)
              enddo
           enddo
        endif
        if (dble_corner(nSE)) then
           ! SE
           do j=0,1-ngh,-1
              do i=ngx+1,ngx+ngh
                 xgc(i,j)=-xgc(i-1,j+1)+xgc(i-1,j)+xgc(i,j+1)
                 ygc(i,j)=-ygc(i-1,j+1)+ygc(i-1,j)+ygc(i,j+1)
              enddo
           enddo
        endif
        if (dble_corner(nSW)) then
           ! SW
           do j=0,1-ngh,-1
              do i=0,1-ngh,-1
                 xgc(i,j)=-xgc(i+1,j+1)+xgc(i+1,j)+xgc(i,j+1)
                 ygc(i,j)=-ygc(i+1,j+1)+ygc(i+1,j)+ygc(i,j+1)
              enddo
           enddo
        endif
     endif
  endif

  ! Intrablock communications of global grid
  ! ========================================
  if (is_curv3) then
     call grid_comm_intrablock_curv3(xgc3,ngx,ngy,ngz)
     call grid_comm_intrablock_curv3(ygc3,ngx,ngy,ngz)
     call grid_comm_intrablock_curv3(zgc3,ngx,ngy,ngz)
  else
     if (is_curv) then
        call grid_comm_intrablock_curv(xgc,ngx,ngy)
        call grid_comm_intrablock_curv(ygc,ngx,ngy)
     else
        call grid_comm_intrablock_cart(xg,ngx)
        call grid_comm_intrablock_cart(yg,ngy)
     endif

     if (.not.(is_2D)) then
        call grid_comm_intrablock_cart(zg,ngz)
     endif
  endif

  call MPI_BARRIER(COMM_global,info)

  ! Free MPI types for grid comm
  ! ============================
  if (is_curv3) then
     !call free_grid_comm_type_curv3
  else
     call free_grid_comm_type_cart
     !if (is_curv) call free_grid_comm_type_curv
  endif

end subroutine grid_comm
