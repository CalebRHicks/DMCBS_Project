module optimizer
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   real(kind=r8), private, save :: delta0,eps
   real(kind=r8), private, save, allocatable :: params(:)
   real(kind=r8), private, save, allocatable :: pmat(:,:)
   real(kind=r8), private, save, allocatable :: ovec(:),ovecw(:)
   real(kind=r8), private, save :: esum,esum2
   integer(kind=i4), private, save :: nparm,npdet,npjas,istep,npart
   logical, private, allocatable, save :: iopt(:)
   real(kind=r8), pointer, save :: pjas(:),pdet(:)
contains
   subroutine setupopt(npartin,optiter)
   use wavefunction
   use potential
   use mympi
   use jaslong
   integer(kind=i4) :: npartin,optiter,np,i
   npart=npartin
   allocate(iopt(2))
   if (myrank().eq.0) then
      read (5,*) optiter ! iterations for the optimization
      read (5,*) iopt(1) ! optimize Jastrow
      read (5,*) iopt(2) ! optimize nodes
      read (5,*) delta0
      read (5,*) eps
   endif
   call bcast(optiter)
   call bcast(iopt)
   call bcast(delta0)
   call bcast(eps)
   nparm=0
   npdet=0
   npjas=0
   if (iopt(1)) then
      call getjasparam(npjas,pjas)
   endif
   if (iopt(2)) then
      call getdetparam(npdet,pdet)
   endif
   nparm=npjas+npdet
   allocate(params(nparm))
   np=1
   if (iopt(1)) then
      params(1:npjas)=pjas
      np=npjas+1
   endif
   if (iopt(2)) then
      params(np:np+npdet-1)=pdet
   endif
   if (myrank().eq.0) then
      write (6,'(''iterations for the optimization'',t40,i10)') optiter
      write (6,'(''optimize Jastrow'',t40,l10)') iopt(1)
      write (6,'(''optimize nodes'',t40,l10)') iopt(2)
      write (6,'(''delta ='',t40,f10.5)') delta0
      write (6,'(''epsilon ='',t40,f10.5)') eps
      write (6,'(''number of parameters ='',t40,i10)') nparm
      do i=1,nparm
         write(6,'(f25.15)') params(i)
      enddo
   endif
   allocate(pmat(nparm,nparm),ovec(nparm),ovecw(nparm))
   istep=0
   ovec=0.0_r8
   ovecw=0.0_r8
   esum=0.0_r8
   esum2=0.0_r8
   pmat=0.0_r8
   end subroutine setupopt

   subroutine updateparams
   use betafun
   use wavefunction
   use potential
   use jaslong
   integer(kind=i4) :: np
   np=1
   if (iopt(1)) then
      pjas=params(1:npjas)
      call setjasparam(pjas)
      np=npjas+1
   endif
   if (iopt(2)) then
      pdet=params(np:np+npdet-1)
      call setdetparam(pdet)
   endif
   end subroutine updateparams

   subroutine calder(w,ene0)
   use stack
   use wavefunction
   use potential
   use jaslong
   type (walker) w
   real(kind=r8) :: dpsi(npdet),djas(npjas),derpsi(nparm)
   real(kind=r8) :: ene0
   integer(kind=i4) :: i,j
   djas=0.0_r8
   dpsi=0.0_r8
   if (iopt(1)) call getjasder(w%x,djas)
   if (iopt(2)) call getderpsi(w%x,dpsi)
   derpsi(1:npjas)=djas
   derpsi(npjas+1:nparm)=dpsi
   ovec=ovec+derpsi
   esum=esum+ene0
   esum2=esum2+ene0**2
   ovecw=ovecw+derpsi*ene0
   do i=1,nparm
      do j=1,nparm
         pmat(i,j)=pmat(i,j)+derpsi(i)*derpsi(j)
      enddo
   enddo
   istep=istep+1
   return
   end subroutine calder

   subroutine updatepar(ene)
   use matrixmod
   use mympi
   integer(kind=i4) :: i,j,is,isteptot
   real(kind=r8) :: det,fvec(nparm),dpar(nparm),newpars(nparm),ene,delta,eerr
   real(kind=r8) :: otot(nparm),owtot(nparm),ptot(nparm,nparm),etot,etot2
   call addall(ovec,otot)
   call addall(ovecw,owtot)
   call addall(pmat,ptot)
   call addall(esum,etot)
   call addall(esum2,etot2)
   call addall(istep,isteptot)
   if (myrank().eq.0) then
      ovec=otot/isteptot
      ovecw=owtot/isteptot
      pmat=ptot/isteptot
      esum=etot/isteptot
      esum2=etot2/isteptot
      eerr=sqrt((abs(esum**2-esum2))/isteptot)
      fvec=-2.0_r8*ovecw+2.0_r8*ovec*esum
      do i=1,nparm
         do j=1,nparm
            pmat(i,j)=pmat(i,j)-ovec(i)*ovec(j)
         enddo
         pmat(i,i)=pmat(i,i)+eps
      enddo
      call matinv(pmat,det,is,nparm)
      do i=1,nparm
         dpar(i)=sum(pmat(i,:)*fvec(:))
      enddo
!     delta=0.5_r8/abs(delta0-esum-sum(ovec*dpar))
      delta=delta0
      newpars=delta*dpar
   endif
   ene=esum
   call bcast(newpars)
   ovec=0.0_r8
   ovecw=0.0_r8
   pmat=0.0_r8
   esum=0.0_r8
   esum2=0.0_r8
   istep=0
   params=params+newpars
   if (myrank().eq.0) then
      write (6,*) 'eave = ',ene,eerr
      write (6,'(''new parameters: '',50000f30.20)') (params(i),i=1,nparm)
      do i=1,nparm
         write (6,'(f30.20)') params(i)
      enddo
   endif
   return
   end subroutine

   subroutine chkderp(w,dp)
   use stack
   use wavefunction
   use potential
   use jaslong
   type (walker) w
   real(kind=r8) :: dpsi(nparm),dpsi1(npjas),dpsi2(npdet),dpsit(nparm)
   integer(kind=i4) :: i
   real(kind=r8) :: dp,psi1,psi2
   integer(kind=i4) :: np
   np=1
   if (iopt(1)) then
      do i=1,npjas
         pjas(i)=pjas(i)+dp
         call setjasparam(pjas)
         call hpsi(w)
         psi1=w%psil
         pjas(i)=pjas(i)-2.0_r8*dp
         call setjasparam(pjas)
         call hpsi(w)
         psi2=w%psil
         pjas(i)=pjas(i)+dp
         call setjasparam(pjas)
         call hpsi(w)
         dpsi(np)=(psi1-psi2)/(2.0_r8*dp)
         np=np+1
      enddo
      call getjasder(w%x,dpsi1)
   endif
   if (iopt(2)) then
      do i=1,npdet
         pdet(i)=pdet(i)+dp
         call setdetparam(pdet)
         call hpsi(w)
         psi1=w%psil
         pdet(i)=pdet(i)-2.0_r8*dp
         call setdetparam(pdet)
         call hpsi(w)
         psi2=w%psil
         pdet(i)=pdet(i)+dp
         call setdetparam(pdet)
         call hpsi(w)
         dpsi(np)=(psi1-psi2)/(2.0_r8*dp)
         np=np+1
      enddo
      call getderpsi(w%x,dpsi2)
   endif
   dpsit(1:npjas)=dpsi1
   dpsit(npjas+1:nparm)=dpsi2
   call hpsi(w)
   dpsi1=0.0_r8
   do i=1,nparm
      write (6,'(''i analytic numerical err'',i5,2e14.6,f25.20)') i,dpsit(i),dpsi(i),dpsi(i)-dpsit(i)
   enddo
   end subroutine chkderp
end module optimizer
