module slaterbcs
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, save :: npart,ndim
   integer(kind=i4), private, save :: npair,nup,ndn
   integer(kind=i4), private, save, allocatable :: nbin(:)
   real(kind=r8), private, save :: el
   integer(kind=i4), private, save :: nexup,nexdn,norb,nk
   real(kind=r8), private, save, allocatable :: akup(:,:),akdn(:,:)
   real(kind=r8), private, save, allocatable :: ak2up(:),ak2dn(:)
   real(kind=r8), private, save, pointer :: ak(:,:),ak2(:)
   logical, private, save, allocatable :: sinup(:),sindn(:)
   real(kind=r8), private, save, allocatable :: mass(:),vorb(:)
   integer(kind=i4), private, save, allocatable :: numsh(:),spin(:)
   integer(kind=i4), private, save :: npdet
   real(kind=r8), private, save, allocatable :: smatisave(:,:)
   logical, private, save :: optsw=.true. 
   logical, private, save :: optbeta=.false.
   real(kind=r8), private, save :: kbes
contains
   subroutine setupbcs(rho,npartin,nbinin,massin,ndin,ntab)
   use kshell
   use mympi
   use betafun
   use coshsolution
   integer(kind=i4) :: npartin,nbinin(:),ndin,ntab
   integer(kind=i4), allocatable :: ikup(:,:),ikdn(:,:)
   real(kind=r8) :: rho,massin(:)
   integer(kind=i4) :: i,ic,j
   real(kind=r8) :: pi
   npart=npartin
   ndim=ndin
   allocate(nbin(2))
   nbin=nbinin
   allocate(mass(npart))
   mass(1:nbin(1))=massin(1)
   mass(nbin(1)+1:npart)=massin(2)
   allocate(spin(npart))
! assign spin 1 to up particles and -1 to down ones
   spin(1:nbin(1))=1
   spin(nbin(1)+1:nbin(1)+nbin(2))=-1
   if (myrank().eq.0) then
      read (5,*) rho ! number density
      read (5,*) nexup
      allocate(ikup(ndim,nexup))
      do j=1,nexup
         read (5,*) ikup(:,j)
      enddo
      read (5,*) nexdn
      allocate(ikdn(ndim,nexdn))
      do j=1,nexdn
         read (5,*) ikdn(:,j)
      enddo
      read (5,*) norb  !number of orbitals for s-wave
      allocate(vorb(norb))
      do i=1,norb
         read (5,*) vorb(i) !amplitude
      enddo
      write (6,'(''rho ='',t40,f10.5)') rho ! number density
      el=(npart/rho)**(1.0_r8/ndim)
      write (6,'(''L/2 ='',t40,f10.5)') el*0.5_r8
      write (6,'(''unpaired up states ='',t40,i10)') nexup
      if (ndim.eq.3) write (6,'(''unpaired k values ='',(t35,3i5))') (ikup(:,j),j=1,nexup)
      if (ndim.eq.2) write (6,'(''unpaired k values ='',(t40,2i5))') (ikup(:,j),j=1,nexup)
      write (6,'(''unpaired down states ='',t40,i10)') nexdn
      if (ndim.eq.3) write (6,'(''unpaired k values ='',(t35,3i5))') (ikdn(:,j),j=1,nexdn)
      if (ndim.eq.2) write (6,'(''unpaired k values ='',(t40,2i5))') (ikdn(:,j),j=1,nexdn)
      write (6,'(''s-wave pairing parameters:'')')
      write (6,'(''number of functions: '',t40,i10)') norb
      write (6,'(''amplitudes: '',t40,f10.5)') (vorb(i),i=1,norb)
   endif
   call bcast(rho)
   call bcast(el)
   call bcast(nexup)
   call bcast(nexdn)
   call bcast(norb)
   if (myrank().ne.0) allocate(vorb(norb))
   call bcast(vorb)
   pi=4.0_r8*atan(1.0_r8)
! initialize single particle states assigned in input
   allocate(akup(ndim,nexup),akdn(ndim,nexdn))
   allocate(sinup(nexup),sindn(nexdn))
   allocate(ak2up(nexup),ak2dn(nexdn))
   if (myrank().eq.0) then
      do i=1,nexup
         akup(:,i)=ikup(:,i)*2.0_r8*pi/el
         ak2up(i)=sum(akup(:,i)**2)
         sinup(i)=.false.
         do ic=1,ndim
            if (ikup(ic,i).ne.0) then
               sinup(i)=(ikup(ic,i).lt.0) !sin if last nonzero entry is negative
            endif
         enddo
      enddo
      do i=1,nexdn
         akdn(:,i)=ikdn(:,i)*2.0_r8*pi/el
         ak2dn(i)=sum(akdn(:,i)**2)
         sindn(i)=.false.
         do ic=1,ndim
            if (ikdn(ic,i).ne.0) then
               sindn(i)=(ikdn(ic,i).lt.0) !sin if last nonzero entry is negative
            endif
         enddo
      enddo
   endif
   call bcast(akup)
   call bcast(ak2up)
   call bcast(sinup)
   call bcast(akdn)
   call bcast(ak2dn)
   call bcast(sindn)
! initialize s-wave pairing function
   allocate(numsh(norb))
   call setupk(el,norb,numsh,ak,ak2,nk,ndim)
   npair=nbin(1)+nbin(2)-nexup-nexdn
   if (mod(npair,2).ne.0) then
      write (6,'(''inconsistent pairing number'',3i5)') &
         nbin(1)+nbin(2),nexup,nexdn
      stop
   endif
   npair=npair/2
   nup=npair+nexup
   ndn=npair+nexdn
   allocate(smatisave(npair+nexup+nexdn,npair+nexup+nexdn))
!  call setbeta(ndim,ntab,el)
   call initcoshfun(ntab,0.5_r8*el,kbes)
   end subroutine setupbcs

   subroutine bcs(x,phil,is,dphi,d2phi)
   use matrixmod
   real(kind=r8), dimension(ndim,npart) :: x,dphi
   real(kind=r8) :: d2phi
   real(kind=r8), dimension(npair+nexup+nexdn,npair+nexup+nexdn) :: smati,d2ph
   real(kind=r8), dimension(ndim,npair+nexup+nexdn,npair+nexup+nexdn) :: dph
   real(kind=r8) :: arg,c,s,phil
   integer(kind=i4) i,j,ii,jj,nd,ic
   integer(kind=i4) :: is
   phil=0.0_r8
   is=1
   nd=npair+nexup+nexdn
   smati=0.0_r8
   dph=0.0_r8
   d2ph=0.0_r8
   call allpair(x,smati(1:nup,1:ndn),dph(:,1:nup,1:ndn),d2ph(1:nup,1:ndn))
   do i=1,npair+nexup
      do j=1,nexup
         jj=j+npair+nexdn
         arg=sum(akup(:,j)*x(:,i))
         c=cos(arg)
         s=sin(arg)
         if (sinup(j)) then
            smati(i,jj)=s
            dph(:,i,jj)=akup(:,j)*c
            d2ph(i,jj)=-ak2up(j)*s
         else
            smati(i,jj)=c
            dph(:,i,jj)=-akup(:,j)*s
            d2ph(i,jj)=-ak2up(j)*c
         endif
      enddo
   enddo
! use -gradient for downs to be consistent with bcs
   do i=1,npair+nexdn
      ii=i+npair+nexup
      do j=1,nexdn
         jj=j+npair+nexup
         arg=sum(akdn(:,j)*x(:,ii))
         c=cos(arg)
         s=sin(arg)
         if (sindn(j)) then
            smati(jj,i)=s
            dph(:,jj,i)=-akdn(:,j)*c
            d2ph(jj,i)=-ak2dn(j)*s
         else
            smati(jj,i)=c
            dph(:,jj,i)=akdn(:,j)*s
            d2ph(jj,i)=-ak2dn(j)*c
         endif
      enddo
   enddo
   call matinv(smati,phil,is,nd)
   dphi=0.0_r8
   d2phi=0.0_r8
   if (phil.lt.-1.e29_r8) return
   do i=1,npair+nexup
      do ic=1,ndim
         dphi(ic,i)=sum(smati(:,i)*dph(ic,i,:))/sqrt(mass(i))
      enddo
      d2phi=d2phi+sum(smati(:,i)*d2ph(i,:))/mass(i)
   enddo
   do i=1,npair+nexdn
      j=npair+nexup+i
      do ic=1,ndim
         dphi(ic,j)=-sum(smati(i,:)*dph(ic,:,i))/sqrt(mass(j))
      enddo
      d2phi=d2phi+sum(smati(i,:)*d2ph(:,i))/mass(j)
   enddo
   smatisave=smati
   end subroutine bcs

   subroutine allpair(x,smati,dph,d2ph)
!  use betafun
   use coshsolution
   real(kind=r8) :: x(:,:),smati(:,:),dph(:,:,:),d2ph(:,:)
   real(kind=r8) :: uporbs(nup,2*nk-1),dnorbs(2*nk-1,ndn)
   real(kind=r8) :: duporbs(ndim,nup,2*nk-1),d2uporbs(nup,2*nk-1)
   real(kind=r8) :: c,s,arg,v,dx(ndim),r
   integer(kind=i4) :: i,k,ksh,j,kk
!  real(kind=r8) :: f,df,d2f
   real(kind=r8) :: f,df(3),d2f,dfr,d2fr
!
! smati, dph, and d2ph are assumed here to be zeroed before entry.
!
! This routine calculates the pairing part of the matrix as matrix
! multiplies. It uses fewer calculations of sines and cosines,
! but more memory.
!
! The matrix uporbs(i,n) = Phi_n(r_i) with r_i an up particle,
! and downorbs(n,j) = Phi_n(r_j) with r_j a down particle. The derivatives
! are calculated with respect to the up particles, so for polarized systems
! it is a more efficient to take up to be the species with the smaller number
! of particles.
!
   if (norb.ne.0) then
      uporbs(:,1)=vorb(1)
      duporbs(:,:,1)=0.0_r8
      d2uporbs(:,1)=0.0_r8
      dnorbs(1,:)=1.0_r8
      k=1
      kk=0
      do ksh=2,norb
         do j=1,numsh(ksh)
            k=k+1
            kk=kk+2
            v=vorb(ksh)*2.0_r8
            do i=1,nup
               arg=sum(ak(:,k)*x(:,i))
               c=cos(arg)
               s=sin(arg)
               uporbs(i,kk)=v*c
               duporbs(:,i,kk)=-ak(:,k)*v*s
               d2uporbs(i,kk)=-ak2(k)*v*c
               uporbs(i,kk+1)=v*s
               duporbs(:,i,kk+1)=ak(:,k)*v*c
               d2uporbs(i,kk+1)=-ak2(k)*v*s
            enddo
            do i=1,ndn
               arg=sum(ak(:,k)*x(:,i+nup))
               dnorbs(kk,i)=cos(arg)
               dnorbs(kk+1,i)=sin(arg)
            enddo
         enddo
      enddo
      smati=matmul(uporbs,dnorbs)
      dph=reshape(matmul(reshape(duporbs(:,:,:),(/ndim*nup,2*nk-1/)),dnorbs) &
         ,shape(dph))
      d2ph=matmul(d2uporbs,dnorbs)
   endif
   do i=1,nup
      do j=1,ndn
         dx=x(:,i)-x(:,j+nup)
         dx=dx-el*nint(dx/el)
         call getcoshfun(dx,f,df,d2f)
         smati(i,j)=smati(i,j)+f
         dph(:,i,j)=dph(:,i,j)+df
         d2ph(i,j)=d2ph(i,j)+d2f
      enddo
   enddo
   end subroutine allpair

   subroutine bcsone(x,dr,i,phiratio)
   use betafun
   real(kind=r8) :: x(:,:),dr(:),phiratio,dx(ndim),arg,c,s,v,r,f,dummy
   real(kind=r8), dimension(npair+nexup+nexdn) :: ph
   integer(kind=i4) :: i,j,jj,k,ish,kk
! x contains all spacial coordinates, dr is xold(i)-xnew(i)
! I must have storaged smati
   x(:,i)=x(:,i)-dr(:)
   if (i.le.nbin(1)) then
      do j=1,npair+nexdn
         jj=j+npair+nexup
         dx(:)=x(:,i)-x(:,jj)
         dx=dx-el*nint(dx/el)
         ph(j)=vorb(1)
         k=1
         do ish=2,norb
            do kk=1,numsh(ish)
               k=k+1
               arg=sum(ak(:,k)*dx(:))
               c=cos(arg)
               s=sin(arg)
               v=vorb(ish)*2.0_r8
               ph(j)=ph(j)+v*c
            enddo
         enddo
         r=sqrt(sum(dx**2))
         call radialfun(r,f,dummy,dummy)
         ph(j)=ph(j)+f
      enddo
      jj=npair+nexdn
      do j=1,nexup
         jj=jj+1
         arg=sum(akup(:,j)*x(:,i))
         c=cos(arg)
         s=sin(arg)
         if (sinup(j)) then
            ph(jj)=s
         else
            ph(jj)=c
         endif
      enddo
   else
      do j=1,npair+nexup
         dx(:)=x(:,i)-x(:,j)
         dx=dx-el*nint(dx/el)
         ph(j)=vorb(1)
         k=1
         do ish=2,norb
            do kk=1,numsh(ish)
               k=k+1
               arg=sum(ak(:,k)*dx(:))
               c=cos(arg)
               s=sin(arg)
               v=vorb(ish)*2.0_r8
               ph(j)=ph(j)+v*c
            enddo
         enddo
         r=sqrt(sum(dx**2))
         call radialfun(r,f,dummy,dummy)
         ph(j)=ph(j)+f
      enddo
      jj=npair+nexup
      do j=1,nexdn
         jj=jj+1
         arg=sum(akdn(:,j)*x(:,i))
         c=cos(arg)
         s=sin(arg)
         if (sindn(j)) then
            ph(jj)=s
         else
            ph(jj)=c
         endif
      enddo
   endif
   if (i.le.nbin(1)) then
      phiratio=sum(smatisave(:,i)*ph)
   else
      phiratio=sum(smatisave(i-nbin(1),:)*ph)
   endif
   x(:,i)=x(:,i)+dr(:)
   end subroutine bcsone

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the following are subrotuines needed for the optimization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine setbcspar(params)
   use coshsolution
   use betafun
   real(kind=r8) :: params(:),bbeta,dbeta
   if (optsw.and..not.optbeta) then
      kbes=params(1)
      call setcoshfun(kbes)
      return
   else if (optsw.and.optbeta) then
      vorb=params(1:npdet-2)
      bbeta=params(npdet-1)
      dbeta=params(npdet)
      call initbeta(bbeta,dbeta)
      return
   else if (.not.optsw.and.optbeta) then
      bbeta=params(1)
      dbeta=params(2)
      call initbeta(bbeta,dbeta)
      return
   endif
   end subroutine setbcspar

   subroutine getbcspar(nparam,params)
   use betafun
   integer(kind=i4) :: nparam
   real(kind=r8), pointer :: params(:)
   real(kind=r8) :: bbeta,dbeta
   if (optsw.and..not.optbeta) then
      nparam=1
      allocate(params(nparam))
      params=kbes
      npdet=nparam
   else if (optsw.and.optbeta) then
      nparam=norb+2
      allocate(params(nparam))
      params(1:nparam-2)=vorb
      call getbetadata(bbeta,dbeta)
      params(nparam-1)=bbeta
      params(nparam)=dbeta
      npdet=nparam
   else if (.not.optsw.and.optbeta) then
      nparam=2
      allocate(params(nparam))
      call getbetadata(bbeta,dbeta)
      params(1)=bbeta
      params(2)=dbeta
      npdet=nparam
   endif
   end subroutine getbcspar

   subroutine getderbcs(x,dpsi)
   use matrixmod
   use coshsolution
   real(kind=r8) :: x(:,:),dpsi(:)
   real(kind=r8) :: philp,philm,d2phi,delta
   real(kind=r8), dimension(ndim,npair+nexup+nexdn,npair+nexup+nexdn) :: dphi
   integer(kind=i4) :: is
   dpsi=0.0_r8
   delta=1.0e-5_r8
   kbes=kbes+delta
   call setcoshfun(kbes)
   call bcs(x,philp,is,dphi,d2phi)
   kbes=kbes-2.0_r8*delta
   call setcoshfun(kbes)
   call bcs(x,philm,is,dphi,d2phi)
   dpsi(1)=(philp-philm)/(2.0_r8*delta)
   end subroutine getderbcs
end module slaterbcs
