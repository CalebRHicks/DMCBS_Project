module pdistribution
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4) ,private, save :: ntab,nsamp
   real(kind=r8), private, save :: el,dk,pi,h
   real(kind=r8), private, save :: wsblk,totalw
   real(kind=r8), private, allocatable, save :: enknow(:),enk(:)
   real(kind=r8), private, allocatable, save :: enkerr(:),enke(:)
   integer(kind=i4), private, allocatable, save :: indk(:),nk(:),ik(:,:)
   real(kind=r8), private, save :: kf
contains
   subroutine setupenk(ndim,ntabin,elin,kfin)
   use mympi
   integer(kind=i4) :: ntabin,ndim
   real(kind=r8) :: elin,kfin
   integer(kind=i4) :: ix,iy,iz,nnk,itmp,i,j,itab,ntot,iktmp(3)
   integer(kind=i4), allocatable :: ik2(:)
   if (ndim.ne.3) call abort
   ntab=ntabin
!change?
   nsamp=10
   el=elin
   pi=4.0_r8*atan(1.0_r8)
   dk=2.0_r8*pi/el
   h=0.5*el/nsamp
   nnk=sqrt(real(ntab,r8))+1
   ntot=0
   do ix=0,nnk
      do iy=0,nnk
         do iz=0,nnk
            if (ix**2+iy**2+iz**2.gt.ntab) cycle
            ntot=ntot+1
         enddo
      enddo
   enddo
   allocate(ik(3,ntot),ik2(ntot),indk(0:ntab),nk(0:ntab))
   ntot=0
   do ix=0,nnk
      do iy=0,nnk
         do iz=0,nnk
            if (ix**2+iy**2+iz**2.gt.ntab) cycle
            ntot=ntot+1
            ik2(ntot)=ix**2+iy**2+iz**2
            ik(1,ntot)=ix
            ik(2,ntot)=iy
            ik(3,ntot)=iz
         enddo
      enddo
   enddo
! stupid bubble sort -- fix if it takes too long
   do i=2,ntot
      do j=1,i-1
         if (ik2(i).lt.ik2(j)) then
            itmp=ik2(i)
            ik2(i)=ik2(j)
            ik2(j)=itmp
            iktmp=ik(:,i)
            ik(:,i)=ik(:,j)
            ik(:,j)=iktmp
         endif
      enddo
   enddo
   nk=0
   itab=0
   indk(0)=0
   nk(0)=1
   do i=2,ntot
      if (ik2(i).ne.ik2(i-1)) then
         itab=itab+1
         indk(itab)=i-1
      endif
      nk(itab)=nk(itab)+1
   enddo
   ntab=itab
   deallocate(ik2)
   allocate(enk(0:ntab),enknow(0:ntab),enkerr(0:ntab),enke(0:ntab))
   kf=kfin
   end subroutine setupenk

   subroutine zeroenk
   enk=0.0_r8
   enknow=0.0_r8
   enkerr=0.0_r8
   enke=0.0_r8
   wsblk=0.0_r8
   end subroutine zeroenk

   subroutine addenk(w,wt)
   use stack
   use random
   use wavefunction
   type (walker) w
   real(kind=r8) :: dx(3),xnew(3),psirat,enktmp(0:ntab),wt
!  real(kind=r8) :: psir1,psir2,r2,rn(1),a(3),amax,amag
!  real(kind=r8) :: akr
   integer(kind=i4) n,i,ii,k,in,j
   call setrn(w%irn)
   call hpsi(w)
   n=size(w%x)/3
   enktmp=0.0_r8
   do i=1,n
      do ii=1,nsamp
         dx=(randn(3)-0.5_r8)*el
         xnew=w%x(:,i)+dx(:)
         call move1(i,xnew,w,psirat)
         do k=0,ntab
            in=indk(k)
            do j=1,nk(k)
               enktmp(k)=enktmp(k)+psirat*cos(dk*ik(1,in+j)*dx(1))* &
                  cos(dk*ik(2,in+j)*dx(2))*cos(dk*ik(3,in+j)*dx(3))
!              akr=dk*sum(ik(:,in+j)*dx(:))
!              enktmp(k)=enktmp(k)+psirat*cos(akr)
            enddo
         enddo
      enddo
   enddo
   enknow(:)=enknow(:)+enktmp/nk*wt
   enkerr(:)=enkerr(:)+(enktmp/nk)**2*wt
   wsblk=wsblk+wt*nsamp
!  enktmp=0.0_r8
!  do i=1,n
!     a=gaussian(3)
!     amax=maxval(abs(a))
!     a=h*a/amax
!     amag=sqrt(sum(a**2))
!     rn=randn(1)
!     do ii=1,nsamp
!        dx=(ii-1+rn(1))*a
!        r2=sum(dx**2)*amag
!        xnew=w%x(:,i)+dx(:)
!        call move1(i,xnew,w,psir1)
!        xnew=w%x(:,i)-dx(:)
!        call move1(i,xnew,w,psir2)
!        do k=0,ntab
!           in=indk(k)
!           do j=1,nk(k)
!              enktmp(k)=enktmp(k)+ &
!                 r2*(psir1+psir2)*cos(dk*ik(1,in+j)*dx(1))* &
!                 cos(dk*ik(2,in+j)*dx(2))*cos(dk*ik(3,in+j)*dx(3))
!           enddo
!        enddo
!     enddo
!  enddo
!  enknow(2,:)=enknow(2,:)+enktmp/nk
!  wtblk(2)=wtblk(2)+w%weight*el**3/(2.0_r8*pi)
   end subroutine addenk
 
   subroutine updateenk
   use mympi
   real(kind=r8) :: norm
   real(kind=r8) :: enktnow(0:ntab),errtnow(0:ntab),err(0:ntab),wtsblk
   call addall(enknow,enktnow)
   call addall(enkerr,errtnow)
   call addall(wsblk,wtsblk)
   if (myrank().eq.0) then
      wsblk=wtsblk
      totalw=wsblk
      norm=1.0_r8/wsblk
      enknow=enktnow
      enkerr=errtnow
      enknow=enknow*norm
      enk=enknow
      enkerr=enkerr*norm
      err=sqrt(abs(enkerr-enknow**2)/wsblk)
      enke=err
   endif
   wsblk=0.0_r8
   enknow=0.0_r8
   enkerr=0.0_r8
   end subroutine updateenk

   subroutine writeenk(filenamek,tau)
   use mympi
   character(len=*) :: filenamek
   integer(kind=i4) :: ir,in
!  real(kind=r8) :: av,av2,enkx,enkxerr
   real(kind=r8) :: tau
   if (myrank().ne.0) return
   open(unit=78,file=filenamek,position='append')
   write(78,'(''# k/kF, n(k), err'')')
   write(78,'(''# total weight = '',e15.7)') totalw
   write(78,'(''# tau = '',e15.7)') tau
   do ir=0,ntab
         in=indk(ir)+1
         write (78,'(6e15.7)') sqrt(dk**2*sum(ik(:,in)**2)),enk(ir),enke(ir)
   enddo
   write(78,*) ''
   write(78,*) ''
   close(78)
!  open(unit=78,file=filenamek2,form='formatted')
!  open(unit=78,file=filenamek2,position='append')
!  do ir=0,ntab
!        enkx=enk(2,ir)/(wttot(2))
!        av=enkbad(2,ir)/(max(1,nblk))
!        av2=enkbad2(2,ir)/(max(1,nblk))
!        enkxerr=sqrt(abs(av2-av*av)/max(1,nblk-1))
!        in=indk(ir)+1
!        write (78,'(6e15.7)') sqrt(dk**2*sum(ik(:,in)**2)),enkx,enkxerr
!  enddo
!  write(78,*) ''
!  write(78,*) ''
!  close(78)
   end subroutine writeenk

end module pdistribution



