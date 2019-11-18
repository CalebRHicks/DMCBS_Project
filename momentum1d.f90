!buggy
module pdistribution
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4) ,private, save :: ntab,ntot,nblk
   real(kind=r8), private, save :: el,h,dk,wtblk,wttot
   real(kind=r8), private, allocatable, save :: rhoz(:)
   real(kind=r8), private, allocatable, save :: enkznow(:),enkz(:),cs(:)
   real(kind=r8), private, allocatable, save :: enkzbad(:),enkzbad2(:)

contains
   subroutine setupenk(ndim,ntabin,elin)
   use trans1d
   use mympi
   integer(kind=i4) :: ntabin,ndim
   real(kind=r8) :: elin
   integer(kind=i4) :: i
   real(kind=r8) :: pi
   if (ndim.ne.3) call abort
   ntab=ntabin
   el=elin
   h=el/ntab
   allocate(rhoz(0:ntab),enkznow(0:ntab),enkz(0:ntab),cs(0:ntab))
   allocate(enkzbad(0:ntab),enkzbad2(0:ntab))
   call tinit(ntab,el)
   pi=4.0_r8*atan(1.0_r8)
   dk=2.0_r8*pi/el
   end subroutine setupenk

   subroutine addenk(w)
   use stack
   use random
   use trans1d
   type (walker) w
   real(kind=r8) a(3),dx(3),xnew(3),psirat,enkztmp(0:ntab)
   integer(kind=i4) n,i,iz,iza
   n=size(w%x)/3
   rhoz=0.0_r8
   a=randn(3)-0.5_r8
   a(1)=a(1)*el
   a(2)=a(2)*el
   a(3)=a(3)*h
   do iz=-ntab,ntab
      iza=abs(iz)
!fix this
      dx(1)=iz*h+a(3)
      dx(2)=iz*h+a(3)
      dx(3)=iz*h+a(3)
      do i=1,n
         xnew=w%x(:,i)+dx(:)
         call move1(i,xnew,w,psirat)
         rhoz(iza)=rhoz(iza)+psirat
      enddo
   enddo
   enkztmp=rtok(rhoz,.true.)
   do iz=0,ntab
      enkznow(iz)=enkznow(iz)+cos(iz*dk*a(3))*enkztmp(iz)
   enddo
   end subroutine addenk

   subroutine move1(i,xnew,w,psirat)
   use stack
   use wavefunction
   type (walker) w
   real(kind=r8) :: xnew(3),psirat
   integer(kind=i4) :: i
   w1=w
   w%x(:,i)=xnew
   call hpsi(w)
   psirat=w%is*w1%is*exp(w%psil-w1%psil)
   w=w1
   end subroutine move1

   subroutine zeroenk
   enkz=0.0_r8
   enkznow=0.0_r8
   enkzbad=0.0_r8
   enkzbad2=0.0_r8
   wttot=0.0_r8
   wtblk=0.0_r8
   nblk=0
   end subroutine zeroenk

   subroutine updateenk
   use mympi
   real(kind=r8) :: anorm,val
   integer(kind=i4) :: i
   real(kind=r8) :: enktnow(0:ntab),wttblk
   call addall(enkznow,enktnow)
   call addall(wtblk,wttblk)
   if (myrank().eq.0) then
      wtblk=wttblk
      enkznow=enktnow
      anorm=1.0_r8/wtblk
      do i=0,ntab
         val=anorm*enkznow(i)
         enkzbad(i)=enkzbad(i)+val
         enkzbad2(i)=enkzbad2(i)+val*val
      enddo
      enkz=enkz+enkznow
      wttot=wttot+wtblk
      nblk=nblk+1
   endif
   wtblk=0.0_r8
   enkznow=0.0_r8
   end subroutine updateenk

   subroutine writeenk(filename)
   use mympi
   character(len=*) :: filename
   integer(kind=i4) :: iz
   real(kind=r8) :: av,av2,err
   if (myrank().ne.0) return
   open(unit=78,file=filename,form='formatted')
   rewind 78
   do iz=0,ntab
      enkz(iz)=enkz(iz)/(max(1,nblk))
      av=enkzbad(iz)/(max(1,nblk))
      av2=enkzbad2(iz)/(max(1,nblk))
      err=sqrt(abs(av2-av*av)/max(1,nblk-1))
      write (78,'(i10,3e15.7)') iz,iz*dk,enkz(iz),err
      write (78,*) av,av2
   enddo
   close(78)
   end subroutine writeenk

end module pdistribution
