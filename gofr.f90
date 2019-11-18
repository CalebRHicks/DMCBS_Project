module gofr
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, save :: npart,nlikep,nunlikep,ntab,ndim
   integer(kind=i4), private, save, allocatable :: nbin(:)
   real(kind=r8), private, save :: el,dr,scale,wsblk,totalw
   real(kind=r8), private, save, allocatable :: gnow(:,:),g(:,:)
   real(kind=r8), private, save, allocatable :: gerr(:,:),ge(:,:)
   real(kind=r8), private, save, allocatable :: dnow(:,:),d(:,:)
   real(kind=r8), private, save, allocatable :: derr(:,:),de(:,:)
   real(kind=r8), private, save :: kf
contains
   subroutine setupgofr(ndin,ntabin,nbinin,elin,kfin)
   integer(kind=i4) :: ndin,ntabin,nbinin(:)
   real(kind=r8) :: elin,kfin
   integer(kind=i4) :: i
   ndim=ndin
   ntab=ntabin
   allocate(nbin(2))
   nbin=nbinin
   el=elin
   allocate(gnow(3,0:ntab),gerr(3,0:ntab))
   allocate(g(3,0:ntab),ge(3,0:ntab))
   allocate(dnow(6,-ntab:ntab),derr(6,-ntab:ntab))
   allocate(d(6,-ntab:ntab),de(6,-ntab:ntab))
   nlikep=0
   npart=0
   do i=1,2
      npart=npart+nbin(i)
      nlikep=nlikep+(nbin(i)*(nbin(i)-1))/2
   enddo
   nunlikep=(npart*(npart-1))/2-nlikep
   call zerogofr
   dr=0.5_r8*el/ntab
   scale=1.0_r8/dr
   kf=kfin
   end subroutine setupgofr

   subroutine zerogofr
   g=0.0_r8
   gnow=0.0_r8
   gerr=0.0_r8
   ge=0.0_r8
   wsblk=0.0_r8
   d=0.0_r8
   dnow=0.0_r8
   derr=0.0_r8
   de=0.0_r8
   end subroutine zerogofr

   subroutine addgofr(x,weight)
   real(kind=r8) :: x(:,:),weight
   integer(kind=i4) :: i,j,itab,index,is
   real(kind=r8) :: dx(ndim),r
   real(kind=r8) :: xcm(3)
   do i=1,npart-1
      do j=i+1,npart
         itab=1   ! up-down
         if (i.le.nbin(1).and.j.le.nbin(1)) itab=2  ! up-up
         if (i.gt.nbin(1).and.j.gt.nbin(1)) itab=3  ! down-down
         dx(:)=x(:,i)-x(:,j)
         dx=dx-el*nint(dx/el)
         r=sqrt(sum(dx**2))
         if (r.lt.0.5_r8*el) then
            index=scale*r
         else
            index=ntab
         endif
         gnow(itab,index)=gnow(itab,index)+2.0_r8*weight
         gerr(itab,index)=gerr(itab,index)+4.0_r8*weight
      enddo
   enddo
   wsblk=wsblk+weight
!  do i=1,3
!     xcm(i)=sum(x(i,:))/npart
!  enddo
   xcm=0.0_r8
   do i=1,npart
      dx=x(:,i)-xcm
      dx=dx-el*nint(dx/el)
      do is=1,3
         if (i.le.nbin(1)) then
            itab=is    ! up particles
         else
            itab=is+3  ! down particles
         endif
         index=scale*dx(is)
         if (index.eq.0) then
            dnow(itab,index)=dnow(itab,index)+0.5_r8*weight
            derr(itab,index)=derr(itab,index)+0.25_r8*weight
         else
            dnow(itab,index)=dnow(itab,index)+weight
            derr(itab,index)=derr(itab,index)+weight
         endif
      enddo
   enddo
   end subroutine addgofr

   subroutine updategofr
   use mympi
   real(kind=r8) :: norm,wtsblk
   real(kind=r8) :: gtnow(3,0:ntab),errtnow(3,0:ntab),err(3,0:ntab)
   real(kind=r8) :: dtnow(6,-ntab:ntab),errdtnow(6,-ntab:ntab),errd(6,-ntab:ntab)
   call addall(gnow,gtnow)
   call addall(gerr,errtnow)
   call addall(wsblk,wtsblk)
   call addall(dnow,dtnow)
   call addall(derr,errdtnow)
   if (myrank().eq.0) then
      wsblk=wtsblk
      totalw=wsblk
      norm=1.0_r8/wsblk
      gnow=gtnow
      gerr=errtnow
      gnow=gnow*norm
      g=gnow
      gerr=gerr*norm
      err=sqrt(abs(gerr-gnow**2)/wsblk)
      ge=err
      dnow=dtnow
      derr=errdtnow
      dnow=dnow*norm
      d=dnow
      derr=derr*norm
      errd=sqrt(abs(derr-dnow**2)/wsblk)
      de=errd
   endif
   wsblk=0.0_r8
   gnow=0.0_r8
   gerr=0.0_r8
   dnow=0.0_r8
   derr=0.0_r8
   end subroutine updategofr

   subroutine writegofr(filename,tau)
   use mympi
   character(len=*) :: filename
   real(kind=r8) :: pi,r0,r1,facl,r
   real(kind=r8) :: vol,rho,tau
   integer(kind=i4) :: i
   if (myrank().ne.0) return
   pi=4.0_r8*atan(1.0_r8)
   open (unit=77,file=filename,position='append')
   write(77,'(''# r*kF, g(r), err'')')
   write(77,'(''# total weight = '',e15.7)') totalw
   write(77,'(''# tau = '',e15.7)') tau
   rho=(npart/el**ndim)
   do i=0,ntab-1
      r0=i*dr
      r1=(i+1)*dr
      r=(i+0.5_r8)*dr
      if (ndim.eq.2) vol=pi*(r1**2-r0**2)
      if (ndim.eq.3) vol=4.0_r8*pi*(r1**3-r0**3)/3.0_r8
      facl=1.0_r8/(rho*npart*vol)
      write (77,'(7e15.7)') r,g(1,i)*facl,ge(1,i)*facl,g(2,i)*facl,ge(2,i)*facl &
               ,g(3,i)*facl,ge(3,i)*facl
   enddo
   write(77,*) ''
   write(77,*) ''
   close(77)
   open (unit=77,file='rho.out',position='append')
   write(77,'(''# r, xrhoup, err, xrhodn, err, yrhoup, err, yrhodn, err, zrhoup, err, zrhodn, err'')')
   write(77,'(''# total weight = '',e15.7)') totalw
   write(77,'(''# tau = '',e15.7)') tau
   rho=(npart/el**ndim)
   do i=-ntab,ntab
      r=i*dr
      facl=1.0_r8/npart
      d(:,i)=d(:,i)*facl
      de(:,i)=de(:,i)*facl
      write (77,'(13e15.7)') r,d(1,i),de(1,i),d(4,i),de(4,i),d(2,i),de(2,i),d(5,i),de(5,i),d(3,i),de(3,i),d(6,i),de(6,i)
   enddo
   write(77,*) ''
   write(77,*) ''
   close(77)
   end subroutine writegofr
end module gofr
