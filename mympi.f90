module mympi
!
! load balancing mpi routines -- these are rewrites by K.E. Schmidt
! of the mpi routines written by Michael A. Lee and I. Lomonosov for the
! parallel version of the Schmidt and Lee electronic structure GFMC code
! 
   implicit none
   include 'mpif.h'
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: i8=selected_int_kind(15)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, save :: irank,iproc,npart
   integer(kind=i4), private, allocatable, save :: numfig(:)
   integer(kind=i4), private, save, allocatable :: nsend(:),nrecv(:),nmoves(:)
   integer(kind=i4), private, save, allocatable :: nwmax(:),nwmin(:)
   integer(kind=i4), private, save :: ncheck,nbal
   real(kind=r8), private, save :: dobal
   real(kind=r8), private, save :: t0,t1,tbal
   integer(kind=i4), private, save :: totalnw
   logical, private, save :: ivmc

interface bcast ! broadcast from process 0
   module procedure bcasti1,bcasti1d,bcasti2d,bcasti3d
   module procedure bcastr1,bcastr1d,bcastr2d,bcastr3d
   module procedure bcastl1,bcastl1d
   module procedure bcastc1d
end interface bcast

interface addall ! return sum to process 0
   module procedure addalli1,addalli1d
   module procedure addallr1,addallr1d,addallr2d
   module procedure addallc1,addallc1d,addallc2d
end interface addall

interface gather ! gather to process 0
   module procedure gatheri1,gatheri1d
   module procedure gatherr1,gatherr1d
end interface gather

contains
   subroutine init0 ! call this before anything else
   integer(kind=i4) :: ierror
   call mpi_init(ierror)
   call mpi_comm_rank(mpi_comm_world,irank,ierror)
   call mpi_comm_size(mpi_comm_world,iproc,ierror)
   allocate(numfig(0:iproc-1))
   if (mpi_integer8.eq.0) then
      write (6,'(''mpi_integer8 not defined'')')
      call abort
   endif
   allocate(nsend(0:iproc-1),nrecv(0:iproc-1),nmoves(0:iproc-1))
   allocate(nwmax(0:iproc-1),nwmin(0:iproc-1))
   nsend=0
   nrecv=0
   nmoves=0
   nwmax=0
   nwmin=0
   ncheck=0
   nbal=0
   totalnw=0
   call cpu_time(t0)
   end subroutine init0

   subroutine init1(npartin,dobalin,ivmcin) !call after everyone has these
   integer(kind=i4) :: npartin
   real(kind=r8) :: dobalin
   logical :: ivmcin
   npart=npartin
   dobal=dobalin
   ivmc=ivmcin
   end subroutine init1

   subroutine done ! wrapper for finalize routine
   integer(kind=i4) :: ierror
   call mpi_finalize(ierror)
   end subroutine done

   subroutine bcasti1(i)
   integer(kind=i4) :: i,ierror
   call mpi_bcast(i,1,mpi_integer,0,mpi_comm_world,ierror)
   return
   end subroutine bcasti1

   subroutine bcasti1d(i)
   integer(kind=i4) :: i(:),ierror
   call mpi_bcast(i,size(i),mpi_integer,0,mpi_comm_world,ierror)
   return
   end subroutine bcasti1d

   subroutine bcasti2d(i)
   integer(kind=i4) :: i(:,:),ierror
   call mpi_bcast(i,size(i),mpi_integer,0,mpi_comm_world,ierror)
   return
   end subroutine bcasti2d

   subroutine bcasti3d(i)
   integer(kind=i4) :: i(:,:,:),ierror
   call mpi_bcast(i,size(i),mpi_integer,0,mpi_comm_world,ierror)
   return
   end subroutine bcasti3d

   subroutine bcastr1d(r)
   integer(kind=i4) :: ierror
   real(kind=r8) :: r(:)
   call mpi_bcast(r,size(r),mpi_double_precision,0,mpi_comm_world,ierror)
   return
   end subroutine bcastr1d

   subroutine bcastr2d(r)
   integer(kind=i4) :: ierror
   real(kind=r8) :: r(:,:)
   call mpi_bcast(r,size(r),mpi_double_precision,0,mpi_comm_world,ierror)
   return
   end subroutine bcastr2d

   subroutine bcastr3d(r)
   integer(kind=i4) :: ierror
   real(kind=r8) :: r(:,:,:)
   call mpi_bcast(r,size(r),mpi_double_precision,0,mpi_comm_world,ierror)
   return
   end subroutine bcastr3d

   subroutine bcastr1(r)
   integer(kind=i4) :: ierror
   real(kind=r8) :: r
   call mpi_bcast(r,1,mpi_double_precision,0,mpi_comm_world,ierror)
   return
   end subroutine bcastr1

   subroutine bcastl1(l)
   integer(kind=i4) :: ierror
   logical :: l
   call mpi_bcast(l,1,mpi_logical,0,mpi_comm_world,ierror)
   return
   end subroutine bcastl1

   subroutine bcastl1d(l)
   integer(kind=i4) :: ierror
   logical :: l(:)
   call mpi_bcast(l,size(l),mpi_logical,0,mpi_comm_world,ierror)
   return
   end subroutine bcastl1d

   subroutine bcastc1d(c)
   integer(kind=i4) :: ierror
   character(len=*) :: c
   call mpi_bcast(c,len(c),mpi_character,0,mpi_comm_world,ierror)
   return
   end subroutine bcastc1d

   function myrank() ! which process am I?
   integer(kind=i4) :: myrank
   myrank=irank
   end function myrank

   function nproc() ! How many of use are there anyway?
   integer(kind=i4) :: nproc
   nproc=iproc
   end function nproc

   subroutine barrier ! wrapper for mpi_barrier
   integer(kind=i4) :: ierror
   call mpi_barrier(mpi_comm_world,ierror)
   end subroutine barrier

   subroutine addalli1(i,isum)
   integer(kind=i4) :: ierror,i,isum
   call mpi_reduce(i,isum,1,mpi_integer,mpi_sum,0,mpi_comm_world,ierror)
   return
   end subroutine addalli1

   subroutine addalli1d(i,isum)
   integer(kind=i4) :: ierror,i(:),isum(:)
   call mpi_reduce(i,isum,size(i),mpi_integer,mpi_sum,0,mpi_comm_world,ierror)
   return
   end subroutine addalli1d

   subroutine addallr1(r,rsum)
   integer(kind=i4) :: ierror
   real(kind=r8) :: r,rsum
   call mpi_reduce(r,rsum,1,mpi_double_precision,mpi_sum,0, &
      mpi_comm_world,ierror)
   return
   end subroutine addallr1

   subroutine addallr1d(r,rsum)
   real(kind=r8) :: r(:),rsum(:)
   integer(kind=i4) :: ierror
   call mpi_reduce(r,rsum,size(r),mpi_double_precision,mpi_sum,0, &
      mpi_comm_world,ierror)
   return
   end subroutine addallr1d

   subroutine addallr2d(r,rsum)
   real(kind=r8) :: r(:,:),rsum(:,:)
   integer(kind=i4) :: ierror
   call mpi_reduce(r,rsum,size(r),mpi_double_precision,mpi_sum,0, &
      mpi_comm_world,ierror)
   return
   end subroutine addallr2d

   subroutine addallc1(c,csum)
   integer(kind=i4) :: ierror
   complex(kind=r8) :: c,csum
   call mpi_reduce(c,csum,1,mpi_double_complex,mpi_sum,0, &
      mpi_comm_world,ierror)
   return
   end subroutine addallc1

   subroutine addallc1d(c,csum)
   complex(kind=r8) :: c(:),csum(:)
   integer(kind=i4) :: ierror
   call mpi_reduce(c,csum,size(c),mpi_double_complex,mpi_sum,0, &
      mpi_comm_world,ierror)
   return
   end subroutine addallc1d

   subroutine addallc2d(c,csum)
   complex(kind=r8) :: c(:,:),csum(:,:)
   integer(kind=i4) :: ierror
   call mpi_reduce(c,csum,size(c),mpi_double_complex,mpi_sum,0, &
      mpi_comm_world,ierror)
   return
   end subroutine addallc2d

   subroutine gatheri1(i,igather)
   integer(kind=i4) :: i,igather(:),ierror
   call mpi_gather(i,1,mpi_integer,igather,1,mpi_integer,0, &
      mpi_comm_world,ierror)
   return
   end subroutine gatheri1

   subroutine gatheri1d(i,igather)
   integer(kind=i4) :: i(:),igather(:,:),ierror
   call mpi_gather(i,size(i),mpi_integer,igather,size(i),mpi_integer,0, &
      mpi_comm_world,ierror)
   return
   end subroutine gatheri1d

   subroutine gatherr1(r,rgather)
   real(kind=r8) :: r,rgather(:)
   integer(kind=i4) :: ierror
   call mpi_gather(r,1,mpi_double_precision,rgather,1, &
      mpi_double_precision,0,mpi_comm_world,ierror)
   return
   end subroutine gatherr1

   subroutine gatherr1d(r,rgather)
   real(kind=r8) :: r(:),rgather(:,:)
   integer(kind=i4) :: ierror
   call mpi_gather(r,size(r),mpi_double_precision,rgather,size(r) &
      ,mpi_double_precision,0,mpi_comm_world,ierror)
   return
   end subroutine gatherr1d

   subroutine abort
   integer(kind=i4) :: ierror
   call mpi_abort(mpi_comm_world,ierror)
   return
   end subroutine abort

   subroutine movewalkers(istack,nwalk,ifrom,ito,id)
   use stack
!
! routine to move nwalk walkers on istack from process ifrom to
! process ito. id is an arbitrary integer identifier
!
   integer(kind=i4) :: nwalk,ifrom,ito,id,i,istack
   integer(kind=i4) :: istatus(mpi_status_size),ierror
   logical :: empty
   if (ifrom.eq.ito) return
   if (irank.eq.ito) then
      do i=1,nwalk
         call mpi_recv(w2%psil,1,mpi_double_precision,ifrom,id, &
            mpi_comm_world,istatus,ierror)
         call mpi_recv(w2%d2psi,1,mpi_double_precision,ifrom,id, &
            mpi_comm_world,istatus,ierror)
         call mpi_recv(w2%v,1,mpi_double_precision,ifrom,id, &
            mpi_comm_world,istatus,ierror)
         call mpi_recv(w2%weight,1,mpi_double_precision,ifrom,id, &
            mpi_comm_world,istatus,ierror)
         call mpi_recv(w2%psigl,1,mpi_double_precision,ifrom,id, &
            mpi_comm_world,istatus,ierror)
         call mpi_recv(w2%d2psig,1,mpi_double_precision,ifrom,id, &
            mpi_comm_world,istatus,ierror)
         call mpi_recv(w2%psil0,1,mpi_double_precision,ifrom,id, &
            mpi_comm_world,istatus,ierror)
         call mpi_recv(w2%psigl0,1,mpi_double_precision,ifrom,id, &
            mpi_comm_world,istatus,ierror)
         call mpi_recv(w2%vext,1,mpi_double_precision,ifrom,id, &
            mpi_comm_world,istatus,ierror)
         call mpi_recv(w2%x,size(w2%x),mpi_double_precision,ifrom,id, &
            mpi_comm_world,istatus,ierror)
         call mpi_recv(w2%dpsi,size(w2%dpsi),mpi_double_precision,ifrom,id, &
            mpi_comm_world,istatus,ierror)
         call mpi_recv(w2%dpsig,size(w2%dpsig),mpi_double_precision,ifrom,id, &
            mpi_comm_world,istatus,ierror)
         call mpi_recv(w2%is,1,mpi_integer,ifrom,id, &
            mpi_comm_world,istatus,ierror)
         call mpi_recv(w2%isg,1,mpi_integer,ifrom,id, &
            mpi_comm_world,istatus,ierror)
         call mpi_recv(w2%irn,1,mpi_integer8,ifrom,id, &
            mpi_comm_world,istatus,ierror)
         call mpi_recv(w2%x0,size(w2%x0),mpi_double_precision,ifrom,id, &
            mpi_comm_world,istatus,ierror)
         call mpi_recv(w2%wt0,1,mpi_double_precision,ifrom,id, &
            mpi_comm_world,istatus,ierror)
         call push(istack,w2)
      enddo
   else if (irank.eq.ifrom) then
      do i=1,nwalk
         call pop(istack,w2,empty)
         if (empty) then
            write (6,'(''stack empty in movewalker'')')
            call abort
         endif
         call mpi_send(w2%psil,1,mpi_double_precision,ito,id, &
            mpi_comm_world,ierror)
         call mpi_send(w2%d2psi,1,mpi_double_precision,ito,id, &
            mpi_comm_world,ierror)
         call mpi_send(w2%v,1,mpi_double_precision,ito,id, &
            mpi_comm_world,ierror)
         call mpi_send(w2%weight,1,mpi_double_precision,ito,id, &
            mpi_comm_world,ierror)
         call mpi_send(w2%psigl,1,mpi_double_precision,ito,id, &
            mpi_comm_world,ierror)
         call mpi_send(w2%d2psig,1,mpi_double_precision,ito,id, &
            mpi_comm_world,ierror)
         call mpi_send(w2%psil0,1,mpi_double_precision,ito,id, &
            mpi_comm_world,ierror)
         call mpi_send(w2%psigl0,1,mpi_double_precision,ito,id, &
            mpi_comm_world,ierror)
         call mpi_send(w2%vext,1,mpi_double_precision,ito,id, &
            mpi_comm_world,ierror)
         call mpi_send(w2%x,size(w2%x),mpi_double_precision,ito,id, &
            mpi_comm_world,ierror)
         call mpi_send(w2%dpsi,size(w2%dpsi),mpi_double_precision,ito,id, &
            mpi_comm_world,ierror)
         call mpi_send(w2%dpsig,size(w2%dpsig),mpi_double_precision,ito,id, &
            mpi_comm_world,ierror)
         call mpi_send(w2%is,1,mpi_integer,ito,id, &
            mpi_comm_world,ierror)
         call mpi_send(w2%isg,1,mpi_integer,ito,id, &
            mpi_comm_world,ierror)
         call mpi_send(w2%irn,1,mpi_integer8,ito,id, &
            mpi_comm_world,ierror)
         call mpi_send(w2%x0,size(w2%x0),mpi_double_precision,ito,id, &
            mpi_comm_world,ierror)
         call mpi_send(w2%wt0,1,mpi_double_precision,ito,id, &
            mpi_comm_world,ierror)
      enddo
   endif
   return
   end subroutine movewalkers

   subroutine load(istack)
   use stack
   integer(kind=i4) :: istack
   integer(kind=i4) :: n1,numfig(0:iproc-1),nt,ierror,nave,isum,k
   integer(kind=i4) :: nr(0:iproc-1),ns(0:iproc-1),npos,nneg,isend,irec,ntrans
   n1=numstack(istack)
   call mpi_gather(n1,1,mpi_integer,numfig,1,mpi_integer,0, &
      mpi_comm_world,ierror)
   call bcast(numfig) !everyone knows how many walkers the others have
   nt=sum(numfig)
   nave=(nt+iproc-1)/iproc
   numfig=numfig-nave
   isum=-sum(numfig)
   numfig(0:isum-1)=numfig(0:isum-1)+1 ! numfig = the number of excess walkers
   npos=0
   nneg=0
   do k=0,iproc-1
      if (numfig(k).gt.0) then
         ns(npos)=k
         npos=npos+1
      endif
      if (numfig(k).lt.0) then
         nr(nneg)=k
         nneg=nneg+1
      endif
   enddo
   isend=0
   irec=0
   do k=1,iproc+1
      if (isend.eq.npos) then
         if (irec.eq.nneg) then
            return
         else
            write (6,'(1x,''error in load'',4i5)') irec,nneg,isend,npos
         endif
      endif
      ntrans=min(numfig(ns(isend)),-numfig(nr(irec)))
      call movewalkers(istack,ntrans,ns(isend),nr(irec),k)
      nsend(ns(isend))=nsend(ns(isend))+1
      nrecv(nr(irec))=nrecv(nr(irec))+1
      numfig(ns(isend))=numfig(ns(isend))-ntrans
      numfig(nr(irec))=numfig(nr(irec))+ntrans
      if (numfig(ns(isend)).eq.0) isend=isend+1
      if (numfig(nr(irec)).eq.0) irec=irec+1
   enddo
   end subroutine load

   subroutine checkpop(istack)
   use stack
   integer(kind=i4) :: istack
   integer(kind=i4) :: n1,numfig(0:iproc-1),nt,ierror,nave
   integer(kind=i4) :: nmax,nmin
   real(kind=r8) :: frac
   real(kind=r8) :: tbal0,tbal1
   call cpu_time(tbal0)
   if (ivmc) then
      call cpu_time(tbal1)
      tbal=tbal+tbal1-tbal0
      return
   endif
   ncheck=ncheck+1
   n1=numstack(istack)
   totalnw=totalnw+n1
   if (n1.gt.nwmax(myrank())) nwmax(myrank())=n1
   if (n1.lt.nwmin(myrank()).or.nwmin(myrank()).eq.0) nwmin(myrank())=n1
   call mpi_gather(n1,1,mpi_integer,numfig,1,mpi_integer,0, &
      mpi_comm_world,ierror)
   call bcast(numfig) !everyone knows how many walkers the others have
   nmax=maxval(numfig)
   nmin=minval(numfig)
   nt=sum(numfig)
   nave=(nt+iproc-1)/iproc
   frac=real(nmax-nmin)/real(nave)
   if (.not.ivmc.and.frac.gt.dobal/100.0_r8) then
      call load(istack)
      nbal=nbal+1
   endif
   call cpu_time(tbal1)
   tbal=tbal+tbal1-tbal0
   return
   end subroutine checkpop

   subroutine printmpilog
   integer(kind=i4) :: i
   integer(kind=i4), allocatable :: nwmaxtot(:),nwmintot(:)
   integer(kind=i4) :: nchecktot,nbaltot
   integer(kind=i4) :: nwtot
   real(kind=r8) :: timetotal,timebal,time
   real(kind=r8) :: fact
   allocate(nwmaxtot(0:iproc-1),nwmintot(0:iproc-1))
   call addall(nwmax,nwmaxtot)
   call addall(nwmin,nwmintot)
   call addall(ncheck,nchecktot)
   call addall(nbal,nbaltot)
   call addall(totalnw,nwtot)
   call addall(tbal,timebal)
   call cpu_time(t1)
   time=t1-t0
   call addall(time,timetotal)
   if (myrank().eq.0) then
      open(unit=10,file='mpi.log',status='unknown')
      write(10,'(''total number of cpu = '',i6)') iproc
      write(10,'(''total time of execution (sec) = '',f10.3)') timetotal/nproc()
      write(10,'(''total time for rebalance (sec) = '',f10.3)') timebal/nproc()
      write(10,'(2/)')
      write(10,'(5x,''cpu #'',7x,''sent'',5x,''received'',2x,''maxnw'',7x,''minnw'',/)')
      do i=0,iproc-1
         write(10,'(5i10)') i,nsend(i),nrecv(i),nwmaxtot(i),nwmintot(i)
      enddo
      write(10,'(/,5x,''total'',4i10)') sum(nsend),sum(nrecv),sum(nwmaxtot),sum(nwmintot)
      write(10,'(2/,''checkpop calls  '',i10)') nchecktot/iproc
      write(10,'(''rebalance calls '',i10)') nbaltot/iproc
      fact=real(sum(nsend))/real(nwtot)
      write(10,'(/,''number of walkers '',i10)') nwtot
      write(10,'(/,''number of walkers moved '',f10.3,'' %'')') fact*100.0_r8
      close(10)
   endif
   end subroutine printmpilog

   subroutine zerotimebal
   tbal=0.0_r8
   return
   end subroutine zerotimebal

end module mympi
