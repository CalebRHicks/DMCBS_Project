module euclidean
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, save :: npart,ndim,nbin(2)
   integer(kind=i4), private, save :: nblk
   integer(kind=i4), private, save :: ntabq,nqtot
   integer(kind=i4), private, save, allocatable :: numsh(:)
   real(kind=r8), private, save :: wsblk,totalw
   complex(kind=r8), private, save, allocatable :: snow(:),s(:),serr(:),se(:)
   complex(kind=r8), private, save, allocatable :: spnow(:),sp(:),sperr(:),spe(:)
   complex(kind=r8), private, save, allocatable :: esnow(:),es(:),eserr(:),ese(:)
   complex(kind=r8), private, save, allocatable :: sinow(:),si(:),sierr(:),sie(:)
   complex(kind=r8), private, save, allocatable :: espnow(:),esp(:),esperr(:),espe(:)
   real(kind=r8), private, save, pointer :: ak(:,:)
   complex(kind=r8), private, parameter :: ci = (0.0_r8,1.0_r8)
   complex(kind=r8), private, parameter :: czero = (0.0_r8,0.0_r8)
   real(kind=r8), private, save :: el,kf
   real(kind=r8), private, save :: qmin
   integer(kind=i4), private, save :: iqmin
   real(kind=r8), private, save :: hbar,etrial
contains
   subroutine setupqeuc(ndin,nbinin,elin,ntabqin,kfin,hbin,etin)
   use kshell
   use mympi
   integer(kind=i4) :: ndin,nbinin(:),ntabqin,nk,i
   real(kind=r8) :: elin,kfin,q,hbin,etin
   real(kind=r8), pointer :: ak2(:)
!  real(kind=r8), pointer :: aktmp(:,:)
   ndim=ndin
   nbin=nbinin
   npart=sum(nbin)
   ntabq=ntabqin
   allocate(numsh(ntabq))
   el=elin
!  call setupk(el,ntabq,numsh,aktmp,ak2,nk,ndim)
   call setupk(el,ntabq,numsh,ak,ak2,nk,ndim)
   nqtot=sum(numsh)
!  allocate(ak(3,2*nqtot-1))
!  ak(:,1)=aktmp(:,1)
!  k=2
!  do i=2,nqtot
!     ak(:,k)=aktmp(:,i)
!     k=k+1
!     ak(:,k)=-aktmp(:,i)
!     k=k+1
!  enddo
!  numsh(2:ntabq)=2*numsh(2:ntabq)
!  nqtot=2*nqtot-1
   allocate(snow(nqtot),s(nqtot),serr(nqtot),se(nqtot))
   allocate(spnow(nqtot),sp(nqtot),sperr(nqtot),spe(nqtot))
   allocate(esnow(nqtot),es(nqtot),eserr(nqtot),ese(nqtot))
   allocate(sinow(nqtot),si(nqtot),sierr(nqtot),sie(nqtot))
   allocate(espnow(nqtot),esp(nqtot),esperr(nqtot),espe(nqtot))
   kf=kfin
   qmin=5.48_r8*kf  ! minimum value of |q|
!  qmin=0.0_r8*kf  ! minimum value of |q|
   call zeroqeuc
   iqmin=0
   do i=1,nqtot
      q=(sum(ak(:,i)**2))**0.5_r8
      if (q.lt.qmin) iqmin=i
   enddo
   iqmin=iqmin+1
   if (myrank().eq.0) then
      write(6,'(''minimum q/kF in euclidean response ='',t40,f10.5)') (sum(ak(:,iqmin)**2))**0.5_r8/kf
      write(6,'(''maximum q/kF in euclidean response ='',t40,f10.5)') (sum(ak(:,nqtot)**2))**0.5_r8/kf
      write(6,'(''total number of q vectors ='',t40,i10)') nqtot-iqmin
   endif
   hbar=hbin
   etrial=etin
   end subroutine setupqeuc

   subroutine zeroqeuc
   snow=0.0_r8
   serr=0.0_r8
   spnow=0.0_r8
   sperr=0.0_r8
   esnow=0.0_r8
   eserr=0.0_r8
   sinow=0.0_r8
   sierr=0.0_r8
   espnow=0.0_r8
   esperr=0.0_r8
   wsblk=0.0_r8
   end subroutine zeroqeuc

   subroutine addqeuc(w)
   use stack
   type (walker) :: w
   integer(kind=i4) :: iq,i
   real(kind=r8) :: arg,dx(ndim),wt,psi,ene
   complex(kind=r8) :: ex,ex0
   complex(kind=r8) :: rhok,rhok0,rhoksp0,rhoksp,sq,rhoen,dot,rhokin,dotsp,rhospen
   integer(kind=i4) :: sigz
   wt=w%weight*exp(w%psil-w%psigl)*w%is*w%wt0
   psi=exp(w%psil)
   ene=-hbar*w%d2psi+w%v-etrial
   do iq=iqmin,nqtot
      rhok=0.0_r8
      rhok0=0.0_r8
      rhoksp=0.0_r8
      rhoksp0=0.0_r8
      rhokin=0.0_r8
      dot=0.0_r8
      dotsp=0.0_r8
      do i=1,npart
         if (i.le.nbin(1)) then
            sigz=1
         else
            sigz=-1
         endif
         dx=w%x(:,i)
         dx=dx-el*nint(dx/el)
         arg=sum(ak(:,iq)*dx(:))
         ex=exp(ci*arg)
         rhok=rhok+ex
         rhoksp=rhoksp+sigz*ex
         dx=w%x0(:,i)
         dx=dx-el*nint(dx/el)
         arg=sum(ak(:,iq)*dx(:))
         ex0=exp(ci*arg)
         rhok0=rhok0+ex0
         rhoksp0=rhoksp0+sigz*ex0
!        dot=dot+dot_product(ak(:,iq),w%dpsi(:,i))
         dot=dot+dot_product(ak(:,iq),w%dpsi(:,i))*ex
         dx=w%x(:,i)-w%x0(:,i)
         dx=dx-el*nint(dx/el)
         arg=sum(ak(:,iq)*dx(:))
         ex=exp(ci*arg)
         rhokin=rhokin+ex
      enddo
      sq=dconjg(rhok0)*rhok
      snow(iq)=snow(iq)+wt*sq
      serr(iq)=serr(iq)+wt*sq**2
      sq=dconjg(rhoksp0)*rhoksp
      spnow(iq)=spnow(iq)+wt*sq
      sperr(iq)=sperr(iq)+wt*sq**2
      dotsp=-2.0_r8*hbar*ci*dot*rhoksp
!     dot=-2.0_r8*hbar*ci*dot*rhok
      dot=-2.0_r8*hbar*ci*dot
      rhospen=ene*rhoksp+dotsp+hbar*sum(ak(:,iq)**2)*rhoksp
      rhoen=ene*rhok+dot+hbar*sum(ak(:,iq)**2)*rhok
      sq=dconjg(rhok0)*rhoen
      esnow(iq)=esnow(iq)+wt*sq
      eserr(iq)=eserr(iq)+wt*sq**2
      sq=dconjg(rhoksp0)*rhospen
      espnow(iq)=espnow(iq)+wt*sq
      esperr(iq)=esperr(iq)+wt*sq**2
      sq=rhokin
      sinow(iq)=sinow(iq)+wt*sq
      sierr(iq)=sierr(iq)+wt*sq**2
   enddo
   wsblk=wsblk+wt
   end subroutine addqeuc

   subroutine updateqeuc
   use mympi
   complex(kind=r8) :: stnow(nqtot),errtnow(nqtot),err(nqtot)
   complex(kind=r8) :: sptnow(nqtot),errsptnow(nqtot),errsp(nqtot)
   complex(kind=r8) :: estnow(nqtot),eerrtnow(nqtot),eerr(nqtot)
   complex(kind=r8) :: stinow(nqtot),errtinow(nqtot),erri(nqtot)
   complex(kind=r8) :: esptnow(nqtot),eerrsptnow(nqtot),eerrsp(nqtot)
   real(kind=r8) :: wtsblk,norm
   call addall(snow,stnow)
   call addall(serr,errtnow)
   call addall(spnow,sptnow)
   call addall(sperr,errsptnow)
   call addall(esnow,estnow)
   call addall(eserr,eerrtnow)
   call addall(sinow,stinow)
   call addall(sierr,errtinow)
   call addall(espnow,esptnow)
   call addall(esperr,eerrsptnow)
   call addall(wsblk,wtsblk)
   if (myrank().eq.0) then
      wsblk=wtsblk
      totalw=wsblk
      norm=1.0_r8/wsblk
      snow=stnow
      serr=errtnow
      snow(:)=snow(:)*norm
      s(:)=snow(:)
      serr(:)=serr(:)*norm
      err(:)=sqrt(abs(serr(:)-snow(:)**2)/wsblk)
      se(:)=err
      spnow=sptnow
      sperr=errsptnow
      spnow(:)=spnow(:)*norm
      sp(:)=spnow(:)
      sperr(:)=sperr(:)*norm
      errsp(:)=sqrt(abs(sperr(:)-spnow(:)**2)/wsblk)
      spe(:)=errsp
      esnow=estnow
      eserr=eerrtnow
      esnow(:)=esnow(:)*norm
      es(:)=esnow(:)
      eserr(:)=eserr(:)*norm
      eerr(:)=sqrt(abs(eserr(:)-esnow(:)**2)/wsblk)
      ese(:)=eerr
      sinow=stinow
      sierr=errtinow
      sinow(:)=sinow(:)*norm
      si(:)=sinow(:)
      sierr(:)=sierr(:)*norm
      erri(:)=sqrt(abs(sierr(:)-sinow(:)**2)/wsblk)
      sie(:)=erri
      espnow=esptnow
      esperr=eerrsptnow
      espnow(:)=espnow(:)*norm
      esp(:)=espnow(:)
      esperr(:)=esperr(:)*norm
      eerrsp(:)=sqrt(abs(esperr(:)-espnow(:)**2)/wsblk)
      espe(:)=eerrsp
   endif
   wsblk=0.0_r8
   snow=0.0_r8
   serr=0.0_r8
   spnow=0.0_r8
   sperr=0.0_r8
   esnow=0.0_r8
   eserr=0.0_r8
   sinow=0.0_r8
   sierr=0.0_r8
   espnow=0.0_r8
   esperr=0.0_r8
   end subroutine updateqeuc

   subroutine writeqeuc(filename,filename2,filename3,tau)
   use mympi
   character(len=*) :: filename,filename2,filename3
   integer(kind=i4) :: iq,k
   real(kind=r8) :: tau,qval
   if (myrank().ne.0) return
   open(unit=78,file=filename,position='append')
   open(unit=79,file=filename2,position='append')
   open(unit=80,file=filename3,position='append')
   open(unit=81,file='euclideanall_inc',position='append')
   open(unit=82,file='euclideanall_sptau',position='append')
   write(78,'(''# tau, |q|/kf, q/kf, s(tau,q), err'')')
   write(78,'(''total weight = '',e15.7)') totalw
   write(79,'(''# tau, |q|/kf, q/kf, sp(tau,q), err'')')
   write(79,'(''total weight = '',e15.7)') totalw
   write(80,'(''# tau, |q/kf|, q/kf, ds/dtau(tau,q), err'')')
   write(80,'(''total weight = '',e15.7)') totalw
   write(81,'(''# tau, |q/kf|, q/kf, s(tau,q), err'')')
   write(81,'(''total weight = '',e15.7)') totalw
   write(82,'(''# tau, |q/kf|, q/kf, dsp/dtau(tau,q), err'')')
   write(82,'(''total weight = '',e15.7)') totalw
   k=0
   do iq=1,nqtot
      qval=sqrt(sum(ak(:,iq)**2))
      if (iq.ge.iqmin) then 
         write(78,'(7e15.7)') tau,qval/kf,ak(:,iq)/kf,real(s(iq)),real(se(iq))
         write(79,'(7e15.7)') tau,qval/kf,ak(:,iq)/kf,real(sp(iq)),real(spe(iq))
         write(80,'(7e15.7)') tau,qval/kf,ak(:,iq)/kf,real(es(iq)),real(ese(iq))
         write(81,'(7e15.7)') tau,qval/kf,ak(:,iq)/kf,real(si(iq)),real(sie(iq))
         write(82,'(7e15.7)') tau,qval/kf,ak(:,iq)/kf,real(esp(iq)),real(espe(iq))
      endif
   enddo
   write(78,*) ' '
   write(78,*) ' '
   close(78)
   write(79,*) ' '
   write(79,*) ' '
   close(79)
   write(80,*) ' '
   write(80,*) ' '
   close(80)
   write(81,*) ' '
   write(81,*) ' '
   close(81)
   write(82,*) ' '
   write(82,*) ' '
   close(82)
   call zeroqeuc
   end subroutine writeqeuc


   subroutine writeqeuc2(filename,filename2,filename3,tau)
   use mympi
   character(len=*) :: filename,filename2,filename3
   integer(kind=i4) :: iq,k,j
   real(kind=r8) :: tau,qval
   real(kind=r8) :: ave,err
   if (myrank().ne.0) return
   open(unit=78,file=filename,position='append')
   write(78,'(''# tau, q, s(tau,q), err'')')
   write(78,'(''total weight = '',e15.7)') totalw
   k=0
   do iq=1,ntabq
      ave=0.0_r8
      err=0.0_r8
      do j=1,numsh(iq)
         k=k+1
         ave=ave+real(s(k))
         err=err+real(se(k))
      enddo
      ave=ave/numsh(iq)
      err=err/numsh(iq)
      qval=sqrt(sum(ak(:,k)**2))
      if (k.ge.iqmin) write(78,'(4e15.7)') tau,qval/kf,ave,err
   enddo
   write(78,*) ' '
   write(78,*) ' '
   close(78)
   open(unit=78,file=filename2,position='append')
   write(78,'(''# tau, q, sp(tau,q), err'')')
   write(78,'(''total weight = '',e15.7)') totalw
   k=0
   do iq=1,ntabq
      ave=0.0_r8
      err=0.0_r8
      do j=1,numsh(iq)
         k=k+1
         ave=ave+real(sp(k))
         err=err+real(spe(k))
      enddo
      ave=ave/numsh(iq)
      err=err/numsh(iq)
      qval=sqrt(sum(ak(:,k)**2))
      if (k.ge.iqmin) write(78,'(4e15.7)') tau,qval/kf,ave,err
   enddo
   write(78,*) ' '
   write(78,*) ' '
   close(78)
   open(unit=78,file=filename3,position='append')
   write(78,'(''# tau, q, ds/dtau(tau,q), err'')')
   write(78,'(''total weight = '',e15.7)') totalw
   k=0
   do iq=1,ntabq
      ave=0.0_r8
      err=0.0_r8
      do j=1,numsh(iq)
         k=k+1
         ave=ave+real(es(k))
         err=err+real(ese(k))
      enddo
      ave=ave/numsh(iq)
      err=err/numsh(iq)
      qval=sqrt(sum(ak(:,k)**2))
      if (k.ge.iqmin) write(78,'(4e15.7)') tau,qval/kf,ave,err
   enddo
   write(78,*) ' '
   write(78,*) ' '
   close(78)
   open(unit=78,file='euclidean_inc',position='append')
   write(78,'(''# tau, q, s(tau,q), err'')')
   write(78,'(''total weight = '',e15.7)') totalw
   k=0
   do iq=1,ntabq
      ave=0.0_r8
      err=0.0_r8
      do j=1,numsh(iq)
         k=k+1
         ave=ave+real(si(k))
         err=err+real(sie(k))
      enddo
      ave=ave/numsh(iq)
      err=err/numsh(iq)
      qval=sqrt(sum(ak(:,k)**2))
      if (k.ge.iqmin) write(78,'(4e15.7)') tau,qval/kf,ave,err
   enddo
   write(78,*) ' '
   write(78,*) ' '
   close(78)
   open(unit=78,file='euclidean_sptau',position='append')
   write(78,'(''# tau, q, dsp/dtau(tau,q), err'')')
   write(78,'(''total weight = '',e15.7)') totalw
   k=0
   do iq=1,ntabq
      ave=0.0_r8
      err=0.0_r8
      do j=1,numsh(iq)
         k=k+1
         ave=ave+real(esp(k))
         err=err+real(espe(k))
      enddo
      ave=ave/numsh(iq)
      err=err/numsh(iq)
      qval=sqrt(sum(ak(:,k)**2))
      if (k.ge.iqmin) write(78,'(4e15.7)') tau,qval/kf,ave,err
   enddo
   write(78,*) ' '
   write(78,*) ' '
   close(78)
   call zeroqeuc
   end subroutine writeqeuc2

end module euclidean
