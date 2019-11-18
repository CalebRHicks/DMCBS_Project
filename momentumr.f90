module pdistribution
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4) ,private, save :: ntab,ntot,nblk
   real(kind=r8), private, save :: el,h,dk,wtblk,wttot,pi
   real(kind=r8), private, allocatable, save :: rhor(:)
   real(kind=r8), private, allocatable, save :: enkrnow(:),enkr(:),cs(:)
   real(kind=r8), private, allocatable, save :: enkrbad(:),enkrbad2(:)

contains
   subroutine setupenk(ndim,ntabin,elin)
!  use fft3d
   use mympi
   integer(kind=i4) :: ntabin,ndim
   real(kind=r8) :: elin
   integer(kind=i4) :: i
   if (ndim.ne.3) call abort
   ntab=ntabin
   el=elin
   h=0.5_r8*el/ntab
   allocate(rhor(0:ntab),enkrnow(0:ntab),enkr(0:ntab),cs(0:ntab))
   allocate(enkrbad(0:ntab),enkrbad2(0:ntab))
!  call fft3dinit(16*ntab,0.5_r8*el)
   pi=4.0_r8*atan(1.0_r8)
   dk=2.0_r8*pi/el
   end subroutine setupenk

   subroutine addenk(w)
   use stack
   use random
   use wavefunction
   type (walker) w
   real(kind=r8) a(3),dx(3),xnew(3),psirat,enkrtmp(0:ntab),wt
   integer(kind=i4) n,i,ir,ira
   wt=w%weight*exp(w%psil-w%psigl)
   call hpsi(w)
   n=size(w%x)/3
   rhor=0.0_r8
   call setrn(w%irn)
   a=gaussian(3)
   a=a/sqrt(sum(a**2))
   rhor(0)=n
   do ir=-ntab,ntab
      ira=abs(ir)
      dx=ir*h*a
      do i=1,n
         xnew=w%x(:,i)+dx(:)
         call move1(i,xnew,w,psirat)
         rhor(ira)=rhor(ira)+psirat
      enddo
   enddo
   enkrnow=enkrnow+0.5_r8*rhor*wt
   wtblk=wtblk+wt
   end subroutine addenk

   subroutine zeroenk
   enkr=0.0_r8
   enkrnow=0.0_r8
   enkrbad=0.0_r8
   enkrbad2=0.0_r8
   wttot=0.0_r8
   wtblk=0.0_r8
   nblk=0
   end subroutine zeroenk

   subroutine updateenk
   use mympi
   real(kind=r8) :: anorm,val
   integer(kind=i4) :: i
   real(kind=r8) :: enktnow(0:ntab),wttblk
   call addall(enkrnow,enktnow)
   call addall(wtblk,wttblk)
   if (myrank().eq.0) then
      wtblk=wttblk
      enkrnow=enktnow
      anorm=1.0_r8/wtblk
      do i=0,ntab
         val=anorm*enkrnow(i)
         enkrbad(i)=enkrbad(i)+val
         enkrbad2(i)=enkrbad2(i)+val*val
      enddo
      enkr=enkr+enkrnow
      wttot=wttot+wtblk
      nblk=nblk+1
   endif
   wtblk=0.0_r8
   enkrnow=0.0_r8
   end subroutine updateenk

   subroutine writeenk(filenamer,filenamek)
   use mympi
!  use fft3d
   character(len=*) :: filenamer,filenamek
   integer(kind=i4) :: ir
   real(kind=r8) :: av,av2,err,r,ak
   real(kind=r8) :: enk(0:16*ntab),enkerr(0:16*ntab),enr(0:ntab),enrerr(0:ntab)
   real(kind=r8) :: enrx(0:16*ntab),enrerrx(0:16*ntab)
   if (myrank().ne.0) return
   rewind 78
   do ir=0,ntab
      r=ir*h
      enr(ir)=enkr(ir)/(wttot)
      av=enkrbad(ir)/(max(1,nblk))
      av2=enkrbad2(ir)/(max(1,nblk))
      enrerr(ir)=sqrt(abs(av2-av*av)/max(1,nblk-1))
   enddo
   enrx=0.0_r8
   enrerrx=0.0_r8
   enrx(0:ntab)=enr
   enrerrx(0:ntab)=enrerr
   open(unit=78,file=filenamer,form='formatted')
   do ir=0,ntab
      r=ir*h
      write (78,'(6e15.7)') r,enr(ir),enrerr(ir)
   enddo
   close(78)
!  open(unit=78,file=filenamek,form='formatted')
!  enk=xfmr2k0(enrx)
!  enkerr=xfmr2k0(enrerrx)
!  do ir=0,16*ntab
!     ak=ir*dk/16
!     write (78,'(6e15.7)') ak,enk(ir),enkerr(ir)
!  enddo
!  close(78)
   end subroutine writeenk

end module pdistribution
