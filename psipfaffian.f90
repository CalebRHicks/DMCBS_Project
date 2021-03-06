module pfaffian
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, save :: npart,ntabb,ndim
   integer(kind=i4), private, save, allocatable :: nbin(:)
   real(kind=r8), private, save :: el
   integer(kind=i4), private, save :: nexup,nexdn,norbsw,norbpw,norbpwup,norbpwdn,nksw
   real(kind=r8), private, save, allocatable :: akup(:,:),akdn(:,:)
   real(kind=r8), private, save, allocatable :: ak2up(:),ak2dn(:)
   real(kind=r8), private, save, pointer :: aksw(:,:),ak2sw(:),akpw(:,:),ak2pw(:)
   real(kind=r8), private, save, allocatable :: akpwup(:,:),ak2pwup(:),akpwdn(:,:),ak2pwdn(:)
   logical, private, save, allocatable :: sinup(:),sindn(:)
   real(kind=r8), private, save, allocatable :: mass(:),vorbsw(:),vorbpw(:),vorbpwup(:),vorbpwdn(:)
   integer(kind=i4), private, save, allocatable :: numshsw(:),numshpw(:),spin(:)
   integer(kind=i4), private, save, allocatable :: numshpwup(:),numshpwdn(:)
   integer(kind=i4), private, save :: npdet
   logical, private, save, allocatable :: iopt(:)
contains
   subroutine setuppfaf(rho,npartin,nbinin,massin,ndin)
   use kshell
   use mympi
   integer(kind=i4) :: npartin,nbinin(:),ndin
   integer(kind=i4), allocatable :: ikup(:,:),ikdn(:,:)
   real(kind=r8) :: rho,massin(:),pi
   integer(kind=i4) :: i,ic,nk,j
   real(kind=r8), pointer :: aktmp(:,:),ak2tmp(:)
   integer(kind=i4), allocatable :: numshtmp(:)
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
      read (5,*) norbsw  !number of orbitals for s-wave
      allocate(vorbsw(norbsw))
      do i=1,norbsw
         read (5,*) vorbsw(i) !amplitude
      enddo
      read (5,*) norbpw  !number of orbitals for up-down p-wave
      allocate(vorbpw(norbpw))
      do i=1,norbpw
         read (5,*) vorbpw(i) !amplitude
      enddo
      read (5,*) norbpwup  !number of orbitals up-up for p-wave
      allocate(vorbpwup(norbpwup)) ! amplitudes for up-up p-wave with no m=0 component
      do i=1,norbpwup
         read (5,*) vorbpwup(i) ! state and aplitude
      enddo
      read (5,*) norbpwdn  !number of orbitals dn-dn for p-wave
      allocate(vorbpwdn(norbpwdn)) ! amplitudes for dn-dn p-wave with no m=0 component
      do i=1,norbpwdn
         read (5,*) vorbpwdn(i) ! state and aplitude
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
      write (6,'(''number of functions: '',t40,i10)') norbsw
      write (6,'(''amplitudes: '',t40,f10.5)') (vorbsw(i),i=1,norbsw)
      write (6,'(''p-wave pairing parameters:'')')
      write (6,'(''number of functions for up-down: '',t40,i10)') norbpw
      write (6,'(''amplitude: '',t40,f10.5)') ((vorbpw(i)),i=1,norbpw)
      write (6,'(''number of functions for up-up: '',t40,i10)') norbpwup
      write (6,'(''amplitude: '',t40,f10.5)') ((vorbpwup(i)),i=1,norbpwup) 
      write (6,'(''number of functions for down-down: '',t40,i10)') norbpwdn
      write (6,'(''amplitude: '',t40,f10.5)') ((vorbpwdn(i)),i=1,norbpwdn)
   endif
   call bcast(rho)
   call bcast(el)
   call bcast(nexup)
   call bcast(nexdn)
   if (myrank().ne.0) allocate(ikup(ndim,nexup),ikdn(ndim,nexdn))
   call bcast(ikup)
   call bcast(ikdn)
   call bcast(norbsw)
   call bcast(norbpw)
   call bcast(norbpwup)
   call bcast(norbpwdn)
   if (myrank().ne.0) allocate(vorbsw(norbsw),vorbpw(norbpw),vorbpwup(norbpwup),vorbpwdn(norbpwdn))
   call bcast(vorbsw)
   call bcast(vorbpw)
   call bcast(vorbpwup)
   call bcast(vorbpwdn)
   pi=4.0_r8*atan(1.0_r8)
! initialize single particle states assigned in input
   allocate(akup(ndim,nexup),akdn(ndim,nexdn))
   allocate(sinup(nexup),sindn(nexdn))
   allocate(ak2up(nexup),ak2dn(nexdn))
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
! initialize s-wave pairing function
   allocate(numshsw(norbsw))
   call setupk(el,norbsw,numshsw,aksw,ak2sw,nksw,ndim)
! initialize p-wave pairing functions
! up-down
   allocate(numshpw(norbpw))
   allocate(numshtmp(norbpw+1))
   call setupk(el,norbpw+1,numshtmp,aktmp,ak2tmp,nk,ndim)
   numshpw(1:norbpw)=numshtmp(2:norbpw+1)
   allocate(akpw(ndim,nk-1),ak2pw(nk-1))
   akpw(:,1:nk-1)=aktmp(:,2:nk)
   ak2pw(1:nk-1)=ak2tmp(2:nk)
   deallocate(aktmp,ak2tmp,numshtmp)
! up-up
   allocate(numshpwup(norbpwup))
   allocate(numshtmp(norbpwup+1))
   call setupk(el,norbpwup+1,numshtmp,aktmp,ak2tmp,nk,ndim)
   numshpwup(1:norbpwup)=numshtmp(2:norbpwup+1)
   allocate(akpwup(ndim,nk-1),ak2pwup(nk-1))
   akpwup(:,1:nk-1)=aktmp(:,2:nk)
   ak2pwup(1:nk-1)=ak2tmp(2:nk)
   deallocate(aktmp,ak2tmp,numshtmp)
! down-down
   allocate(numshpwdn(norbpwdn))
   allocate(numshtmp(norbpwdn+1))
   call setupk(el,norbpwdn+1,numshtmp,aktmp,ak2tmp,nk,ndim)
   numshpwdn(1:norbpwdn)=numshtmp(2:norbpwdn+1)
   allocate(akpwdn(ndim,nk-1),ak2pwdn(nk-1))
   akpwdn(:,1:nk-1)=aktmp(:,2:nk)
   ak2pwdn(1:nk-1)=ak2tmp(2:nk)
   deallocate(aktmp,ak2tmp,numshtmp)
!!!!!!!!!!! UGLY !!!!!
   allocate(iopt(6))
   iopt=.false.
   if (norbsw.ne.0) iopt(3)=.true.
   if (norbpw.ne.0) iopt(4)=.true.
   if (norbpwup.ne.0) iopt(5)=.true.
   if (norbpwdn.ne.0) iopt(6)=.true.
   end subroutine setuppfaf

   subroutine pfaff(x,phil,is,dphi,d2phi)
   use matrixmod
   use pfaffianmod
   real(kind=r8), dimension(ndim,npart) :: x,dphi
   real(kind=r8) :: d2phi
   real(kind=r8), dimension(npart+nexup+nexdn,npart+nexup+nexdn) :: smati,d2ph,pfaftmp
   real(kind=r8), dimension(ndim,npart+nexup+nexdn,npart+nexup+nexdn) :: dph
   real(kind=r8), dimension(ndim) :: dx
   real(kind=r8) :: arg,c,s,phil,psi
   integer(kind=i4) i,j,jj,is,nd,ic
   nd=npart+nexup+nexdn
   smati=0.0_r8
   dph=0.0_r8
   d2ph=0.0_r8
   dphi=0.0_r8
   do i=1,npart-1
      do j=i+1,npart
         if (i.ne.j) then
            dx(:)=x(:,i)-x(:,j)
            dx=dx-el*nint(dx/el)
            call pairfn(dx,smati(i,j),dph(:,i,j),d2ph(i,j),spin(i)+spin(j))
            smati(j,i)=-smati(i,j)
            dph(:,j,i)=dph(:,i,j)
            d2ph(j,i)=-d2ph(i,j)
         endif
      enddo
   enddo
   do i=1,nbin(1)
      jj=npart
      do j=1,nexup
         jj=jj+1
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
         smati(jj,i)=-smati(i,jj)
         dph(:,jj,i)=dph(:,i,jj)
         d2ph(jj,i)=-d2ph(i,jj)
      enddo
   enddo
   do i=nbin(1)+1,nbin(1)+nbin(2)
      jj=npart+nexup
      do j=1,nexdn
         jj=jj+1
         arg=sum(akdn(:,j)*x(:,i))
         c=cos(arg)
         s=sin(arg)
         if (sindn(j)) then
            smati(i,jj)=s
            dph(:,i,jj)=akdn(:,j)*c
            d2ph(i,jj)=-ak2dn(j)*s
         else
            smati(i,jj)=c
            dph(:,i,jj)=-akdn(:,j)*s
            d2ph(i,jj)=-ak2dn(j)*c
         endif
         smati(jj,i)=-smati(i,jj)
         dph(:,jj,i)=dph(:,i,jj)
         d2ph(jj,i)=-d2ph(i,jj)
      enddo
   enddo
   smati(npart+1:nd,npart+1:nd)=0.0_r8
   dph(:,npart+1:nd,npart+1:nd)=0.0_r8
   d2ph(npart+1:nd,npart+1:nd)=0.0_r8
   pfaftmp=smati
   call pfaf(pfaftmp,nd,psi)
   call matinv(smati,phil,is,nd)
   phil=log(abs(psi))
   is=sign(1.0_r8,psi)
   d2phi=0.0_r8
   do i=1,npart
      do ic=1,ndim
         dphi(ic,i)=sum(smati(:,i)*dph(ic,i,:))/sqrt(mass(i))
      enddo
      d2phi=d2phi+sum(smati(:,i)*d2ph(i,:))/mass(i)
   enddo
   end subroutine pfaff

   subroutine pairfn(x,pf,dpf,d2pf,ispin)
   use betafun
   real(kind=r8), dimension(:) :: x,dpf
   real(kind=r8), dimension(ndim) :: dx
   real(kind=r8) :: pf,d2pf,r,c,s,arg,v
   integer(kind=i4) :: i,j,ish,ispin
   real(kind=r8) :: f,df,d2f
   real(kind=r8) :: bet,ex
   pf=0.0_r8
   dpf=0.0_r8
   d2pf=0.0_r8
   if (ispin.eq.0) then
      if (norbsw.ne.0) then
         pf=vorbsw(1)
         i=1
         do ish=2,norbsw
            do j=1,numshsw(ish)
               i=i+1
               arg=sum(aksw(:,i)*x(:))
               c=cos(arg)
               s=sin(arg)
               v=vorbsw(ish)*2.0_r8
               pf=pf+v*c
               dpf(:)=dpf(:)-v*s*aksw(:,i)
               d2pf=d2pf-v*c*ak2sw(i)
            enddo
         enddo
      endif
      dx=x
      r=sqrt(sum(dx**2))
      call radialfun(r,f,df,d2f)
      pf=pf+f
      dpf=dpf+df*dx/r
      d2pf=d2pf+d2f
      if (norbpw.ne.0) then
         i=0
         do ish=1,norbpw
            do j=1,numshpw(ish)
               i=i+1
               arg=sum(akpw(:,i)*x(:))
               c=cos(arg)
               s=sin(arg)
               v=vorbpw(ish)
               pf=pf+v*s
               dpf(:)=dpf(:)+v*c*akpw(:,i)
               d2pf=d2pf-v*s*ak2pw(i)
            enddo
         enddo
      endif
   endif
   if (ispin.eq.2) then  ! up-up p-wave pairing
      i=0
      do ish=1,norbpwup
         do j=1,numshpwup(ish)
            i=i+1
            arg=sum(akpwup(:,i)*x(:))
            c=cos(arg)
            s=sin(arg)
            v=vorbpwup(ish)
            pf=pf+v*s
            dpf(:)=dpf(:)+v*c*akpwup(:,i)
            d2pf=d2pf-v*s*ak2pwup(i)
         enddo
      enddo
! set beta in the Gaussian to have zero derivative at the edge of the box
      bet=0.5_r8*sqrt(2.0_r8/el)
      ex=exp(-bet*x(1)**2)
      pf=pf+x(1)*ex
      dpf(1)=dpf(1)+(1.0_r8-2.0_r8*bet*x(1)**2)*ex
      d2pf=d2pf+bet*x(1)*(4.0_r8*bet*x(1)**2-6.0_r8)*ex
   endif
   if (ispin.eq.-2) then
      i=0
      do ish=1,norbpwdn   ! down-down p-wave pairing
         do j=1,numshpwdn(ish)
            i=i+1
            arg=sum(akpwdn(:,i)*x(:))
            c=cos(arg)
            s=sin(arg)
            v=vorbpwdn(ish)
            pf=pf+v*s
            dpf(:)=dpf(:)+v*c*akpwdn(:,i)
            d2pf=d2pf-v*s*ak2pwdn(i)
         enddo
      enddo
      bet=0.5_r8*sqrt(2.0_r8/el)
      ex=exp(-bet*x(2)**2)
      pf=pf+x(2)*ex
      dpf(2)=dpf(2)+(1.0_r8-2.0_r8*bet*x(2)**2)*ex
      d2pf=d2pf+bet*x(2)*(4.0_r8*bet*x(2)**2-6.0_r8)*ex
   endif
   end subroutine pairfn

   subroutine setderpfaf(nparmdet,ioptin)
   integer(kind=i4) :: nparmdet
   logical :: ioptin(:)
   npdet=nparmdet
   allocate(iopt(size(ioptin)))
   iopt=ioptin
   end subroutine setderpfaf

   subroutine getderpfaf(x,dpsi)
   use matrixmod
   real(kind=r8) :: x(:,:),dpsi(:)
   real(kind=r8) :: dphi(npart+nexup+nexdn,npart+nexup+nexdn,npdet)
   real(kind=r8), dimension(npart+nexup+nexdn,npart+nexup+nexdn) :: smati,d2ph
   real(kind=r8), dimension(ndim,npart+nexup+nexdn,npart+nexup+nexdn) :: dph
   real(kind=r8), dimension(ndim) :: dx
   real(kind=r8) :: arg,phil
   integer(kind=i4) i,j,jj,is,nd
   nd=npart+nexup+nexdn
   smati=0.0_r8
   do i=1,npart-1
      do j=i+1,npart
         if (i.ne.j) then
            dx(:)=x(:,i)-x(:,j)
            dx=dx-el*nint(dx/el)
            call pairfn(dx,smati(i,j),dph(:,i,j),d2ph(i,j),spin(i)+spin(j))
            smati(j,i)=-smati(i,j)
         endif
      enddo
   enddo
   do i=1,nbin(1)
      jj=npart
      do j=1,nexup
         jj=jj+1
         arg=sum(akup(:,j)*x(:,i))
         if (sinup(j)) then
            smati(i,jj)=sin(arg)
         else
            smati(i,jj)=cos(arg)
         endif
         smati(jj,i)=-smati(i,jj)
      enddo
   enddo
   do i=nbin(1)+1,nbin(1)+nbin(2)
      jj=npart+nexup
      do j=1,nexdn
         jj=jj+1
         arg=sum(akdn(:,j)*x(:,i))
         if (sindn(j)) then
            smati(i,jj)=sin(arg)
         else
            smati(i,jj)=cos(arg)
         endif
         smati(jj,i)=-smati(i,jj)
      enddo
   enddo
   smati(npart+1:nd,npart+1:nd)=0.0_r8
   call matinv(smati,phil,is,nd)
   dphi=0.0_r8
   do i=1,npart-1
      do j=i+1,npart
         dx(:)=x(:,i)-x(:,j)
         dx=dx-el*nint(dx/el)
         call derpairfn(dx,dphi(i,j,:),spin(i)+spin(j))
      enddo
   enddo
   dpsi=0.0_r8
   do j=1,npdet
      do i=1,npart+nexup+nexdn
         dpsi(j)=dpsi(j)+sum(smati(:,i)*dphi(i,:,j))
      enddo
   enddo
   end subroutine getderpfaf

   subroutine derpairfn(x,dphi,ispin)
   real(kind=r8) :: x(:),dphi(:)
   real(kind=r8) :: c,arg
   integer(kind=i4) :: i,j,ish,np,ispin
   dphi=0.0_r8
   np=1
   if (iopt(3)) then
      if (ispin.eq.0) dphi(np)=1.0_r8
      np=np+1
      i=1
      do ish=2,norbsw
         if (ispin.eq.0) then
            do j=1,numshsw(ish)
               i=i+1
               arg=sum(aksw(:,i)*x(:))
               c=cos(arg)
               dphi(np)=dphi(np)+2.0_r8*c
            enddo
         endif
         np=np+1
      enddo
   endif
   if (iopt(4)) then
      i=0
      do ish=1,norbpw
         if (ispin.eq.0) then
            do j=1,numshpw(ish)
               i=i+1
               arg=sum(akpw(:,i)*x(:))
               dphi(np)=dphi(np)+sin(arg)
            enddo
         endif
         np=np+1
      enddo
   endif
   if (iopt(5)) then
      i=0
      do ish=1,norbpwup
         if (ispin.eq.2) then
            do j=1,numshpwup(ish)
               i=i+1
               arg=sum(akpwup(:,i)*x(:))
               dphi(np)=dphi(np)+sin(arg)
            enddo
         endif
         np=np+1
      enddo
   endif
   if (iopt(6)) then
      i=0
      do ish=1,norbpwdn
         if (ispin.eq.-2) then
            do j=1,numshpwdn(ish)
               i=i+1
               arg=sum(akpwdn(:,i)*x(:))
               dphi(np)=dphi(np)+sin(arg)
            enddo
         endif
         np=np+1
      enddo
   endif
   end subroutine derpairfn

   subroutine setvorbpfaf(params)
   real(kind=r8) :: params(:)
   vorbsw=params(1:norbsw)
   vorbpw=params(norbsw+1:norbsw+norbpw)
   vorbpwup=params(norbsw+norbpw+1:norbsw+norbpw+norbpwup)
   vorbpwdn=params(norbsw+norbpw+norbpwup+1:norbsw+norbpw+norbpwup+norbpwdn)
   end subroutine setvorbpfaf

   subroutine getpfafampl(nparam,params)
   integer(kind=i4) :: nparam
   real(kind=r8), pointer :: params(:)
   nparam=norbsw+norbpw+norbpwup+norbpwdn
   allocate(params(nparam))
   params(1:norbsw)=vorbsw
   params(norbsw+1:norbsw+norbpw)=vorbpw
   params(norbsw+norbpw+1:norbsw+norbpw+norbpwup)=vorbpwup
   params(norbsw+norbpw+norbpwup+1:norbsw+norbpw+norbpwup+norbpwdn)=vorbpwdn
   npdet=nparam
   end subroutine getpfafampl
end module pfaffian
