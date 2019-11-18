module slaterbcsloff
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, save :: npart,ntabb,ndim
   integer(kind=i4), private, save :: npair,nup,ndn
   integer(kind=i4), private, save, allocatable :: nbin(:)
   real(kind=r8), private, save :: el
   integer(kind=i4), private, save :: nexup,nexdn,norb,nk
   real(kind=r8), private, save, allocatable :: akup(:,:),akdn(:,:)
   real(kind=r8), private, save, allocatable :: ak2up(:),ak2dn(:)
   real(kind=r8), private, save, allocatable :: ak(:,:,:),ak2(:,:)
   logical, private, save, allocatable :: sinup(:),sindn(:)
   real(kind=r8), private, save, allocatable :: mass(:),vorb(:)
   integer(kind=i4), private, save, allocatable :: numsh(:),spin(:)
   integer(kind=i4), private, save :: npdet
   logical, private, save, allocatable :: iopt(:)
   real(kind=r8), private, save :: q(3),q2,fact
   integer(kind=i4), private, save, allocatable :: numshloff(:)
contains
   subroutine setupbcsloff(rho,npartin,nbinin,massin,ndin)
   use kshell
   use mympi
   integer(kind=i4) :: npartin,nbinin(:),ndin
   integer(kind=i4), allocatable :: ikup(:,:),ikdn(:,:)
   real(kind=r8) :: rho,massin(:)
   integer(kind=i4) :: i,ic,j,k
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
      open(unit=10,file='loff.dat',status='old')
      read (10,*) norb
      allocate(numshloff(norb))
      allocate(vorb(norb))
      do i=1,norb
         read (10,*) vorb(i),numshloff(i)
      enddo
      nk=sum(numshloff)
      allocate(ak(2,3,nk))
      do i=1,nk
         read (10,*) ak(1,:,i),ak(2,:,i)
      enddo
      close(unit=10)
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
      write (6,'(''number of functions: '',t40,i10)') nk
      k=0
      do i=1,norb
         write (6,'(''degeneracy and amplitude: '',t25,i10,f20.15)') numshloff(i),vorb(i)
         do j=1,numshloff(i)
            k=k+1
            write (6,'(''state: '',t40,3i5,t60,3i5)') int(ak(1,:,k)),int(ak(2,:,k))
         enddo
      enddo
   endif
   call bcast(rho)
   call bcast(el)
   call bcast(nexup)
   call bcast(nexdn)
   call bcast(norb)
   call bcast(nk)
   if (myrank().ne.0) then
      allocate(vorb(norb))
      allocate(numshloff(norb))
      allocate(ak(2,3,nk))
   endif
   call bcast(vorb)
   call bcast(numshloff)
   call bcast(ak)
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
   fact=2.0_r8*pi/el
   ak=ak*fact
   allocate(ak2(2,nk))
   do i=1,nk
      ak2(1,i)=sum(ak(1,:,i)**2)
      ak2(2,i)=sum(ak(2,:,i)**2)
   enddo
   npair=nbin(1)+nbin(2)-nexup-nexdn
   if (mod(npair,2).ne.0) then
      write (6,'(''inconsistent pairing number'',3i5)') &
         nbin(1)+nbin(2),nexup,nexdn
      stop
   endif
   npair=npair/2
   nup=npair+nexup
   ndn=npair+nexdn
   end subroutine setupbcsloff

   subroutine bcsloff(x,phil,is,dphi,d2phi)
   use matrixmod
   real(kind=r8), dimension(ndim,npart) :: x,dphi
   real(kind=r8) :: d2phi
   real(kind=r8), dimension(npair+nexup+nexdn,npair+nexup+nexdn) :: smati
   real(kind=r8), dimension(2,npair+nexup+nexdn,npair+nexup+nexdn) :: d2ph
   real(kind=r8), dimension(2,ndim,npair+nexup+nexdn,npair+nexup+nexdn) :: dph
   real(kind=r8) :: arg,c,s,phil
   integer(kind=i4) i,j,ii,jj,nd,ic
   integer(kind=i4) :: is
   phil=0.0_r8
   is=1
   nd=npair+nexup+nexdn
   smati=0.0_r8
   dph=0.0_r8
   d2ph=0.0_r8
   do i=1,npair+nexup
      do j=1,npair+nexdn
         jj=j+npair+nexup
         call pairfn(x(:,i),x(:,jj),smati(i,j),dph(:,:,i,j),d2ph(:,i,j))
      enddo
   enddo
   do i=1,npair+nexup
      do j=1,nexup
         jj=j+npair+nexdn
         arg=sum(akup(:,j)*x(:,i))
         c=cos(arg)
         s=sin(arg)
         if (sinup(j)) then
            smati(i,jj)=s
            dph(1,:,i,jj)=akup(:,j)*c
            d2ph(1,i,jj)=-ak2up(j)*s
         else
            smati(i,jj)=c
            dph(1,:,i,jj)=-akup(:,j)*s
            d2ph(1,i,jj)=-ak2up(j)*c
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
            dph(2,:,jj,i)=akdn(:,j)*c
            d2ph(2,jj,i)=-ak2dn(j)*s
         else
            smati(jj,i)=c
            dph(2,:,jj,i)=-akdn(:,j)*s
            d2ph(2,jj,i)=-ak2dn(j)*c
         endif
      enddo
   enddo
   call matinv(smati,phil,is,nd)
   dphi=0.0_r8
   d2phi=0.0_r8
   if (phil.lt.-1.e29_r8) return
   do i=1,npair+nexup
      do ic=1,ndim
         dphi(ic,i)=sum(smati(:,i)*dph(1,ic,i,:))/sqrt(mass(i))
      enddo
      d2phi=d2phi+sum(smati(:,i)*d2ph(1,i,:))/mass(i)
   enddo
   do i=1,npair+nexdn
      j=npair+nexup+i
      do ic=1,ndim
         dphi(ic,j)=sum(smati(i,:)*dph(2,ic,:,i))/sqrt(mass(j))
      enddo
      d2phi=d2phi+sum(smati(i,:)*d2ph(2,:,i))/mass(j)
   enddo
   end subroutine bcsloff

   subroutine pairfn(x1,x2,pf,dpf,d2pf)
   use betafun
   real(kind=r8), dimension(3) :: x1,x2,dx
   real(kind=r8) :: dpf(2,3),dpfb(3)
   real(kind=r8) :: pf,d2pf(2),r,c,s,arg,v,pfb,d2pfb
   integer(kind=i4) :: i,j,k
   real(kind=r8) :: f,df,d2f
   pf=0.0_r8
   dpf=0.0_r8
   d2pf=0.0_r8
   pfb=0.0_r8
   dpfb=0.0_r8
   d2pfb=0.0_r8
   if (norb.eq.0) return
   k=0
   do i=1,norb
      do j=1,numshloff(i)
         k=k+1
         arg=sum(ak(1,:,k)*x1(:)+ak(2,:,k)*x2(:))
         c=cos(arg)
         s=sin(arg)
         v=vorb(i)*2.0_r8
         pf=pf+v*c
         dpf(1,:)=dpf(1,:)-v*s*ak(1,:,k)
         dpf(2,:)=dpf(2,:)-v*s*ak(2,:,k)
         d2pf(1)=d2pf(1)-v*c*ak2(1,k)
         d2pf(2)=d2pf(2)-v*c*ak2(2,k)
      enddo
   enddo
! cm has q**2=1 
   arg=0.5_r8*fact*(x1(1)+x2(1))
   c=2.0_r8*cos(arg)
   s=2.0_r8*sin(arg)
   pfb=pfb+c
   dpfb(1)=dpfb(1)-0.5_r8*fact*s
   d2pfb=d2pfb-0.25_r8*fact**2*c

   arg=0.5_r8*fact*(x1(2)+x2(2))
   c=2.0_r8*cos(arg)
   s=2.0_r8*sin(arg)
   pfb=pfb+c
   dpfb(2)=dpfb(2)-0.5_r8*fact*s
   d2pfb=d2pfb-0.25_r8*fact**2*c

   arg=0.5_r8*fact*(x1(3)+x2(3))
   c=2.0_r8*cos(arg)
   s=2.0_r8*sin(arg)
   pfb=pfb+c
   dpfb(3)=dpfb(3)-0.5_r8*fact*s
   d2pfb=d2pfb-0.25_r8*fact**2*c
!pfb=1.0_r8

   dx=x1-x2
   dx=dx-el*nint(dx/el)
   r=sqrt(sum(dx**2))
   call radialfun(r,f,df,d2f)
   pf=pf+pfb*f
   dpf(1,:)=dpf(1,:)+pfb*df*dx/r+dpfb(:)*f
   dpf(2,:)=dpf(2,:)-pfb*df*dx/r+dpfb(:)*f
   d2pf(1)=d2pf(1)+pfb*d2f+d2pfb*f+2.0_r8*df*sum(dx(:)*dpfb(:))/r
   d2pf(2)=d2pf(2)+pfb*d2f+d2pfb*f-2.0_r8*df*sum(dx(:)*dpfb(:))/r
   return
   end subroutine pairfn

   subroutine setderbcsloff(nparmdet,ioptin)
   integer(kind=i4) :: nparmdet
   logical :: ioptin(:)
   npdet=nparmdet
   allocate(iopt(size(ioptin)))
   iopt=ioptin
   end subroutine setderbcsloff

   subroutine getderbcsloff(x,dpsi)
   use matrixmod
   use betafun
   real(kind=r8) :: x(:,:),dpsi(:)
   real(kind=r8) :: arg,phil,x1(3),x2(3)
   real(kind=r8), dimension(npair+nexup+nexdn,npair+nexup+nexdn,npdet) :: dphi
   real(kind=r8), dimension(npair+nexup+nexdn,npair+nexup+nexdn) :: smati
   real(kind=r8), dimension(2,npair+nexup+nexdn,npair+nexup+nexdn) :: d2ph
   real(kind=r8), dimension(2,ndim,npair+nexup+nexdn,npair+nexup+nexdn) :: dph
   integer(kind=i4) :: i,j,jj,ii,nd,is,np,ish,k,kk
   dpsi=0.0_r8
   nd=npair+nexup+nexdn
   smati=0.0_r8
   dph=0.0_r8
   d2ph=0.0_r8
   do i=1,npair+nexup
      do j=1,npair+nexdn
         jj=j+npair+nexup
         call pairfn(x(:,i),x(:,jj),smati(i,j),dph(:,:,i,j),d2ph(:,i,j))
      enddo
   enddo
   do i=1,npair+nexup
      do j=1,nexup
         jj=j+npair+nexdn
         arg=sum(akup(:,j)*x(:,i))
         if (sinup(j)) then
            smati(i,jj)=sin(arg)
         else
            smati(i,jj)=cos(arg)
         endif
      enddo
   enddo
   do i=1,npair+nexdn
      ii=i+npair+nexup
      do j=1,nexdn
         jj=j+npair+nexup
         arg=sum(akdn(:,j)*x(:,ii))
         if (sindn(j)) then
            smati(jj,i)=sin(arg)
         else
            smati(jj,i)=cos(arg)
         endif
      enddo
   enddo
   call matinv(smati,phil,is,nd)
   dphi=0.0_r8
   do i=1,npair+nexup
      do j=1,npair+nexdn
         jj=j+npair+nexup
         x1=x(:,i)
         x2=x(:,jj)
         np=1
         if (iopt(3)) then
            k=0
            do ish=1,norb
               do kk=1,numshloff(ish)
                  k=k+1
                  arg=sum(ak(1,:,k)*x1(:)+ak(2,:,k)*x2(:))
                  dphi(i,j,np)=dphi(i,j,np)+2.0_r8*cos(arg)
               enddo
               np=np+1
            enddo
         endif
      enddo
   enddo
   do j=1,npdet
      do i=1,npair+nexup+nexdn
         dpsi(j)=dpsi(j)+sum(smati(:,i)*dphi(i,:,j))
      enddo
   enddo
   end subroutine getderbcsloff

   subroutine setvorbbcsloff(v0)
   real(kind=r8) :: v0(:)
   vorb=v0
   end subroutine setvorbbcsloff

   subroutine getbcsamplloff(nsw,vsw)
   integer(kind=i4) :: nsw
   real(kind=r8), pointer :: vsw(:)
   nsw=norb
   allocate(vsw(nsw))
   vsw=vorb
   end subroutine getbcsamplloff




end module slaterbcsloff
