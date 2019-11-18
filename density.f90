module density
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, save :: npart,ntab
   integer(kind=i4), private, save :: nblk
   integer(kind=i4), private, save, allocatable :: nbin(:)
   real(kind=r8), private, save :: dr,scale
   real(kind=r8), private, save, allocatable :: rnow(:,:),rr(:,:),wttot(:)
   real(kind=r8), private, save, allocatable :: rbad(:,:),rbad2(:,:),wtblk(:)
   real(kind=r8), private, save, allocatable :: mass(:)
contains
   subroutine setuprho(ntabin,nbinin,rmax,massin)
   integer(kind=i4) :: ntabin
   integer(kind=i4) :: nbinin(:)
   real(kind=r8) :: rmax,massin(:)
   ntab=ntabin
   allocate(nbin(2))
   allocate(rnow(2,0:ntab),rbad(2,0:ntab),rbad2(2,0:ntab))
   allocate(rr(2,0:ntab),wttot(2),wtblk(2))
   nbin=nbinin
   npart=sum(nbin)
   call zerorho
   dr=rmax/ntab
   scale=1.0_r8/dr
   allocate(mass(npart))
   mass(1:nbin(1))=massin(1)
   if (nbin(1).ne.npart) mass(nbin(1)+1:npart)=massin(2)
   end subroutine setuprho

   subroutine zerorho
   wttot=0.0_r8
   wtblk=0.0_r8
   nblk=0
   rr=0.0_r8
   rnow=0.0_r8
   rbad=0.0_r8
   rbad2=0.0_r8
   end subroutine zerorho

   subroutine addrho(x,weight)
   real(kind=r8) :: x(3,npart),weight,r,xcm(3),massi
   integer(kind=i4) :: i,itab,index
   massi=1.0_r8/sum(mass(:))
   xcm=0.0_r8
!  xcm(1)=(sum(x(1,:)*mass(:)))*massi
!  xcm(2)=(sum(x(2,:)*mass(:)))*massi
!  xcm(3)=(sum(x(3,:)*mass(:)))*massi
   do i=1,npart
      r=sqrt(sum((x(:,i)-xcm(:))**2))
      index=scale*r
      index=min(index,ntab)
      if (i.le.nbin(1)) then
         itab=1 ! up particle
      else
         itab=2 ! down particle
      endif
      rnow(itab,index)=rnow(itab,index)+weight
      wtblk(itab)=wtblk(itab)+weight
   enddo
   end subroutine addrho

   subroutine updaterho
   use mympi
   real(kind=r8) :: anorm(2),val(2)
   integer(kind=i4) :: i
   real(kind=r8) :: rtnow(2,0:ntab),wttblk(2)
   call addall(rnow,rtnow)
   call addall(wtblk,wttblk)
   if (myrank().eq.0) then
      wtblk=wttblk
      rnow=rtnow
      anorm=1.0_r8/wtblk
      do i=0,ntab
         val=anorm*rnow(:,i)
         rbad(:,i)=rbad(:,i)+val
         rbad2(:,i)=rbad2(:,i)+val*val
      enddo
      rr=rr+rnow
      wttot=wttot+wtblk
      nblk=nblk+1
   endif
   wtblk=0.0_r8
   rnow=0.0_r8
   end subroutine updaterho

   subroutine writerho(filename)
   use mympi
   character(len=*) :: filename
   real(kind=r8) :: av(2),av2(2),err(2),r0,r1,r2,vol,pi
   integer(kind=i4) :: i
   if (myrank().ne.0) return
   pi=4.0_r8*atan(1.0_r8)
!  open (unit=78,file=filename,position='append')
   open (unit=78,file=filename)
   do i=0,ntab-1
      r0=(i+0.5_r8)*dr
      r1=i*dr
      r2=r1+dr
! this works in 3d only!!!
      vol=4.0_r8*pi*(r2**3-r1**3)/3.0_r8
      rr(:,i)=rr(:,i)/(wttot*vol)
      av=rbad(:,i)/max(1,nblk)
      av2=rbad2(:,i)/max(1,nblk)
      err=sqrt(abs(av2-av*av)/max(1,nblk-1))/vol
      write (78,'(7e15.7)') r0,rr(1,i),err(1),rr(2,i),err(2)
   enddo
   write(78,*) ''
   write(78,*) ''
   close(78)
   end subroutine writerho
end module density
