module step
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   real(kind=r8), private, save :: hbar,dt,driftm,etrial,el,nfac,siga,sigb,wtcut
   real(kind=r8), private, save, allocatable :: mass(:),sigma(:),dcut(:)
   integer(kind=i4), private, save :: idmc,npart,ndim,nbin(2)
   logical, private, save :: cutdrift,fn
   integer(kind=i4), private, save :: wtisto(500)
contains
   subroutine setstep(idmcin,ndin,npartin,nbinin,hbarin,dtin,etrialin,elin,massin,wtcutin,dcutin,eunit,fnin)
   integer(kind=i4) :: idmcin,npartin,nbinin(:),ndin
   real(kind=r8) :: hbarin,dtin,etrialin,elin,sig,massin(:),wtcutin,dcutin,eunit
   logical :: fnin
   idmc=idmcin
   npart=npartin
   nbin=nbinin
   ndim=ndin
! here I need hbar**2/2 without the mass
   hbar=hbarin
   dt=dtin
   el=elin
   etrial=etrialin
   nfac=1.0_r8/eunit
   if (el.ne.0.0_r8) nfac=nfac/npart
   allocate(mass(npart),sigma(npart),dcut(npart))
   mass(1:nbin(1))=massin(1)
   mass(nbin(1)+1:nbin(1)+nbin(2))=massin(2)
   sig=sqrt(2.0_r8*hbar*dt)
   sigma=sig/sqrt(mass)
   wtcut=wtcutin
   if (dcutin.ge.0.0_r8) then
      dcut=sigma*dcutin
      cutdrift=.true.
   else
      cutdrift=.false.
   endif
   wtisto=0
   fn=fnin
   return
   end subroutine setstep

   subroutine step1(istin,istout,opcomp)
   use stack
   use random
   use wavefunction
   use estimator
   integer(kind=i4) :: istin,istout,i
   real(kind=r8), dimension (ndim,npart) :: gauss,dx,drifto,driftn
   real(kind=r8), dimension(ndim) :: dx1,dx2
   real(kind=r8) :: eo,en,prob,rn(1),ejf,wt,dot,ego,egn
   real(kind=r8) :: rat1,rat2,wt1,wt2,wtg
   logical :: empty,opcomp
   type(walker) :: wtmp
   if (idmc.eq.10) allocate(wtmp%x(ndim,npart),wtmp%dpsi(ndim,npart),wtmp%dpsig(ndim,npart),wtmp%x0(ndim,npart))
   do while (.true.)
      call pop(istin,w1,empty)
      if (empty) exit
      if (el.ne.0.0_r8) w1%x=w1%x-el*nint(w1%x/el)
      call setrn(w1%irn)
      if (idmc.ne.10) then
         eo=-hbar*w1%d2psi+w1%v
         gauss=reshape(gaussian(ndim*npart),(/ndim,npart/))
         do i=1,npart
            gauss(:,i)=sigma(i)*gauss(:,i)
            if (cutdrift) then
               drifto(:,i)=min(dcut(i),max(-dcut(i),w1%dpsig(:,i)*sigma(i)**2*sqrt(mass(i))))
            else
               drifto(:,i)=w1%dpsig(:,i)*sigma(i)**2*sqrt(mass(i))
            endif
         enddo
         dx=drifto+gauss
         w2%x=w1%x+dx
         if (el.ne.0.0_r8) w2%x=w2%x-el*nint(w2%x/el)
         call hpsi(w2)
         en=-hbar*w2%d2psi+w2%v
         select case (idmc)
           case(-1:0)
               w2%weight=w1%weight
               prob=2.0_r8*(w2%psil-w1%psil)
               do i=1,npart
                  if (cutdrift) then
                     driftn(:,i)=min(dcut(i),max(-dcut(i),w2%dpsig(:,i)*sigma(i)**2*sqrt(mass(i))))
                  else
                     driftn(:,i)=w2%dpsig(:,i)*sigma(i)**2*sqrt(mass(i))
                  endif
                  dx1=drifto(:,i)+driftn(:,i)
                  dx2=drifto(:,i)-driftn(:,i)-2.0_r8*dx(:,i)
                  prob=prob+dot_product(dx1,dx2)/(2.0_r8*sigma(i)**2)
               enddo
               prob=exp(min(0.0_r8,prob))
               rn=randn(1)
               if (rn(1).lt.prob) then
                  dot=0.0_r8
                  do i=1,npart
                     dx1=w2%dpsi(:,i)
                     dot=dot+dot_product(dx1,dx1)/mass(i)
                  enddo
                  ejf=-0.5_r8*hbar*(w2%d2psi-dot)+w2%v
                  call addval(1,en*nfac,w2%weight)
                  call addval(2,(en-w2%v)*nfac,w2%weight)
                  call addval(3,w2%v*nfac,w2%weight)
                  call addval(4,w2%vext*nfac,w2%weight)
                  call addval(5,ejf*nfac,w2%weight)
                  call addval(6,1.0_r8,1.0_r8)
                  call savern(w2%irn)
                  w2%x0=w1%x0
!                 w2%wt0=w1%wt0
                  w2%wt0=w2%is  ! I don't want to have negative weights in VMC
                  call push(istout,w2)
               else
                  dot=0.0_r8
                  do i=1,npart
                     dx1=w1%dpsi(:,i)
                     dot=dot+dot_product(dx1,dx1)/mass(i)
                  enddo
                  ejf=-0.5_r8*hbar*(w1%d2psi-dot)+w1%v
                  call addval(1,eo*nfac,w1%weight)
                  call addval(2,(eo-w1%v)*nfac,w1%weight)
                  call addval(3,w1%v*nfac,w1%weight)
                  call addval(4,w1%vext*nfac,w1%weight)
                  call addval(5,ejf*nfac,w1%weight)
                  call addval(6,0.0_r8,1.0_r8)
                  call savern(w1%irn)
                  call push(istout,w1)
               endif
            case(-2,1)
               ego=-hbar*w1%d2psig+w1%v
               egn=-hbar*w2%d2psig+w2%v
               wt=exp(-(0.5_r8*(ego+egn)-etrial)*dt)
               w2%weight=min(wtcut,wt)
               if (w2%is.ne.w1%is) then
                  if (fn) then
                     w2%weight=0.0_r8
                     wt=0.0_r8
                  endif
                  call addval(5,1.0_r8,1.0_r8)
               else
                  call addval(5,0.0_r8,1.0_r8)
               endif
               wtg=wt*exp(w2%psil-w2%psigl)*w2%is*w1%wt0
               if (opcomp) then
                  call addval(1,en*nfac,wtg)
                  call addval(2,(en-w2%v)*nfac,wtg)
                  call addval(3,w2%v*nfac,wtg)
                  call addval(4,w2%vext*nfac,wtg)
                  call addval(6,0.5_r8*(ego+egn)*nfac,1.0_r8)
                  call addval(7,exp(w2%psil-w2%psigl),1.0_r8)
                  call addval(8,w2%is*exp(w2%psil),1.0_r8)
                  call addval(9,exp(w2%psigl),1.0_r8)
                  call addval(10,1.0_r8*w2%is,1.0_r8)
                  call addval(11,wt,1.0_r8)
                  call addval(12,wtg,1.0_r8)
                  call addval(13,sign(1.0_r8,wtg),1.0_r8)
                  call addval(14,sign(1.0_r8,wtg),wtg)
               endif
               call savern(w2%irn)
               w2%x0=w1%x0
               w2%wt0=w1%wt0
               call push(istout,w2)
            case(2)
               prob=2.0_r8*(w2%psigl-w1%psigl)
               do i=1,npart
                  if (cutdrift) then
                     driftn(:,i)=min(dcut(i),max(-dcut(i),w2%dpsig(:,i)*sigma(i)**2*sqrt(mass(i))))
                  else
                     driftn(:,i)=w2%dpsig(:,i)*sigma(i)**2*sqrt(mass(i))
                  endif
                  dx1=drifto(:,i)+driftn(:,i)
                  dx2=drifto(:,i)-driftn(:,i)-2.0_r8*dx(:,i)
                  prob=prob+dot_product(dx1,dx2)/(2.0_r8*sigma(i)**2)
               enddo
               prob=exp(min(0.0_r8,prob))
               rn=randn(1)
               if (rn(1).lt.prob) then
                  wt=exp(-(0.5_r8*(eo+en)-etrial)*dt)
                  w2%weight=min(wtcut,wt)
                  if (w2%is.ne.w1%is) then
                     if (fn) then
                        w2%weight=0.0_r8
                        wt=0.0_r8
                     endif
                     call addval(4,1.0_r8,1.0_r8)
                  else 
                     call addval(4,0.0_r8,1.0_r8)
                  endif
                  wtg=wt*exp(w2%psil-w2%psigl)*w2%is*w1%wt0
                  if (opcomp) then
                     call addval(1,en*nfac,wtg)
                     call addval(2,(en-w2%v)*nfac,wtg)
                     call addval(3,w2%v*nfac,wtg)
                     call addval(6,exp(w2%psil-w2%psigl),1.0_r8)
                     call addval(7,w2%is*exp(w2%psil),1.0_r8)
                     call addval(8,exp(w2%psigl),1.0_r8)
                     call addval(9,1.0_r8*w2%is,1.0_r8)
                     call addval(10,wt,1.0_r8)
                     call addval(11,wtg,1.0_r8)
                     call addval(12,sign(1.0_r8,wtg),1.0_r8)
                     call addval(13,sign(1.0_r8,wtg),wtg)
                  endif
                  call savern(w2%irn)
                  w2%x0=w1%x0
                  w2%wt0=w1%wt0
                  call addval(5,1.0_r8,1.0_r8)
                  call push(istout,w2)
               else
                  wt=exp(-(eo-etrial)*dt)
                  wtg=wt*exp(w1%psil-w1%psigl)*w1%is*w1%wt0
                  w1%weight=min(wtcut,wt*w1%weight)
                  if (opcomp) then
                     call addval(1,eo*nfac,wtg)
                     call addval(2,(eo-w1%v)*nfac,wtg)
                     call addval(3,w1%v*nfac,wtg)
                     call addval(6,exp(w1%psil-w1%psigl),1.0_r8)
                     call addval(7,w1%is*exp(w1%psil),1.0_r8)
                     call addval(8,exp(w1%psigl),1.0_r8)
                     call addval(9,1.0_r8*w1%is,1.0_r8)
                     call addval(10,wt,1.0_r8)
                     call addval(11,wtg,1.0_r8)
                     call addval(12,sign(1.0_r8,wtg),1.0_r8)
                     call addval(13,sign(1.0_r8,wtg),wtg)
                  endif
                  call savern(w1%irn)
                  call addval(5,0.0_r8,1.0_r8)
                  call push(istout,w1)
               endif
            case default
               write (6,'(''Illegal idmc value -- step'',i10)') idmc
               stop
         end select
      else  ! do direct sampling of GF with importance sampling
         gauss=reshape(gaussian(ndim*npart),(/ndim,npart/))
         do i=1,npart
            gauss(:,i)=sigma(i)*gauss(:,i)
         enddo
         dx=gauss
         w2%x=w1%x+dx
         call hpsi(w2)
         rat1=exp(w2%psigl-w1%psigl)
         wt1=exp(-w2%v*dt)*rat1
         wtmp=w2
         w2%x=w1%x-dx
         call hpsi(w2)
         rat2=exp(w2%psigl-w1%psigl)
         wt2=exp(-w2%v*dt)*rat2
         wt=wt1+wt2
         wt1=wt1/wt
         rn=randn(1)
         if (rn(1).lt.wt1.or.w2%is.ne.w1%is) w2=wtmp ! use first new walker
         if (w2%is.ne.w1%is) then
            if (fn) wt=0.0_r8
            call addval(5,1.0_r8,1.0_r8)
         else
            call addval(5,0.0_r8,1.0_r8)
         endif
         wt=0.5_r8*wt*exp(etrial*dt)
         w2%weight=wt
         en=-hbar*w2%d2psi+w2%v
         wtg=wt*exp(w2%psil-w2%psigl)*w2%is*w1%wt0
         if (opcomp) then
            call addval(1,en*nfac,wtg)
            call addval(2,(en-w2%v)*nfac,wtg)
            call addval(3,w2%v*nfac,wtg)
            call addval(4,w2%vext*nfac,wtg)
            call addval(6,0.5_r8*(ego+egn)*nfac,1.0_r8)
         endif
         call addval(7,exp(w2%psil-w2%psigl),1.0_r8)
         call addval(8,w2%is*exp(w2%psil),1.0_r8)
         call addval(9,exp(w2%psigl),1.0_r8)
         call addval(10,1.0_r8*w2%is,1.0_r8)
         call addval(11,wt,1.0_r8)
         call addval(12,wtg,1.0_r8)
         call addval(13,sign(1.0_r8,wtg),1.0_r8)
         call addval(14,sign(1.0_r8,wtg),wtg)
         call savern(w2%irn)
         w2%x0=w1%x0
         w2%wt0=w1%wt0
         call push(istout,w2)
      endif
   enddo
   if (idmc.eq.10) deallocate(wtmp%x,wtmp%dpsi,wtmp%dpsig)
   istout=istin
   istin=3-istin
   return
   end subroutine step1

   subroutine branch(istin,istout,factor)
   use stack
   use random
   use random2
   integer(kind=i4) :: istin,istout
   real(kind=r8) :: wt,rn(1),dummy,factor
   integer(kind=i4) :: i,iwt
   logical :: empty
   do while (.true.) 
      call pop(istin,w1,empty)
      if (empty) then
         istout=istin
         istin=3-istin
         return
      endif
      wt=w1%weight*factor
      w1%weight=sign(1.0_r8,w1%weight)
      call ran1(rn(1),w1%irn)
      iwt=wt+rn(1)
      do i=1,iwt
         call ran2(dummy,w1%irn)
         call push(istout,w1)
      enddo
      if (iwt.lt.500.and.iwt.ge.0) wtisto(iwt+1)=wtisto(iwt+1)+1
   enddo
   write (6,'(''You should never reach this point!'')')
   stop
   return
   end subroutine branch

   subroutine writewtisto
   use mympi
   integer(kind=i4) :: wttot(500),i
   if (idmc.lt.1) return
   call addall(wtisto,wttot)
   if (myrank().eq.0) then
      write (6,*) ''
      write (6,*) 'walkers multiplicity'
      do i=1,20
         write (6,*) i-1,wttot(i)
      enddo
      do i=21,500
         if (wttot(i).ne.0) write (6,*) i-1,wttot(i)
      enddo
      write (6,*) ''
      write (6,*) ''
   endif
   end subroutine writewtisto
end module step
