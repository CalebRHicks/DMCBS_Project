module sofq
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, save :: npart,ndim,nbin(2)
   integer(kind=i4), private, save :: ntabq,nqtot
   integer(kind=i4), private, save, allocatable :: numsh(:)
   real(kind=r8), private, save :: wsblk,totalw
   real(kind=r8), private, save, allocatable :: snow(:,:),s(:,:),serr(:,:),se(:,:)
   real(kind=r8), private, save, pointer :: ak(:,:)
   complex(kind=r8), private, parameter :: ci = (0.0_r8,1.0_r8)
   complex(kind=r8), private, parameter :: czero = (0.0_r8,0.0_r8)
   real(kind=r8), private, save :: kf
contains
   subroutine setupsofq(ndin,nbinin,el,ntabqin,kfin)
   use kshell
   integer(kind=i4) :: ndin,nbinin(:),ntabqin,nk
   real(kind=r8) :: el,kfin
   real(kind=r8), pointer :: ak2(:)
   ndim=ndin
   nbin=nbinin
   npart=sum(nbin)
   ntabq=ntabqin
   kf=kfin
   allocate(numsh(ntabq))
   call setupk(el,ntabq,numsh,ak,ak2,nk,ndim)
   nqtot=sum(numsh)
   allocate(snow(3,nqtot),s(3,nqtot),serr(3,nqtot),se(3,nqtot))
   call zerosofq
   end subroutine setupsofq

   subroutine zerosofq
   s=0.0_r8
   snow=0.0_r8
   serr=0.0_r8
   se=0.0_r8
   wsblk=0.0_r8
   end subroutine zerosofq

   subroutine addsofq(x,weight)
   real(kind=r8) :: x(3,npart),weight
   integer(kind=i4) :: i,iq
   real(kind=r8) :: arg
   complex(kind=r8) :: rhokup,rhokdn
   do iq=1,nqtot
      rhokup=czero
      rhokdn=czero
      do i=1,nbin(1)
         arg=sum(ak(:,iq)*x(:,i))
         rhokup=rhokup+exp(-ci*arg)
      enddo
      do i=1+nbin(1),nbin(1)+nbin(2)
         arg=sum(ak(:,iq)*x(:,i))
         rhokdn=rhokdn+exp(-ci*arg)
      enddo
! S_upup,S_dndn,S_updn,S_dnup
      snow(1,iq)=snow(1,iq)+rhokup*conjg(rhokup)*weight
      serr(1,iq)=serr(1,iq)+(rhokup*conjg(rhokup))**2*weight
      snow(2,iq)=snow(2,iq)+rhokdn*conjg(rhokdn)*weight
      serr(2,iq)=serr(2,iq)+(rhokdn*conjg(rhokdn))**2*weight
      snow(3,iq)=snow(3,iq)+rhokup*conjg(rhokdn)*weight
      serr(3,iq)=serr(3,iq)+(rhokup*conjg(rhokdn))**2*weight
   enddo
   wsblk=wsblk+weight
   end subroutine addsofq

   subroutine updatesofq
   use mympi
   real(kind=r8) :: norm
   real(kind=r8) :: stnow(3,nqtot),errtnow(3,nqtot),err(3,nqtot),wtsblk
   call addall(snow,stnow)
   call addall(serr,errtnow)
   call addall(wsblk,wtsblk)
   if (myrank().eq.0) then
      wsblk=wtsblk
      totalw=wsblk
      norm=1.0_r8/wsblk
      snow=stnow
      serr=errtnow
      snow=snow*norm
      s=snow
      serr=serr*norm
      err=sqrt(abs(serr-snow**2)/wsblk)
      se=err
   endif
   wsblk=0.0_r8
   snow=0.0_r8
   serr=0.0_r8
   end subroutine updatesofq

   subroutine writesofq(filename,tau)
   use mympi
   character(len=*) :: filename
   real(kind=r8) :: ave(3),err(3)
   real(kind=r8) :: q,tau
   integer(kind=i4) :: iq,j,k
   if (myrank().ne.0) return
   open(unit=79,file=filename,position='append')
   write(79,'(''# q/kF, s(q), err'')')
   write(79,'(''# total weight = '',e15.7)') totalw
   write(79,'(''# tau = '',e15.7)') tau
   k=0
   do iq=1,ntabq
      ave=0.0_r8
      err=0.0_r8
      do j=1,numsh(iq)
         k=k+1
         ave=ave+s(:,k)
         err=err+se(:,k)
      enddo
      ave=ave/numsh(iq)/npart
      err=err/numsh(iq)/npart
      q=sqrt(sum(ak(:,k)**2))
      write (79,'(7e15.7)') q/kf,ave(1),err(1),ave(2),err(2),ave(3),err(3)
   enddo
   write(79,*) ' '
   write(79,*) ' '
   close(79)
   end subroutine writesofq
end module sofq
