module slater
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, save :: npart,ntabb,ndim
   integer(kind=i4), private, save :: npair,nup,ndn
   integer(kind=i4), private, save, allocatable :: nbin(:)
   real(kind=r8), private, save :: el
   real(kind=r8), private, save, allocatable :: akup(:,:),akdn(:,:)
   real(kind=r8), private, save, allocatable :: ak2up(:),ak2dn(:)
   real(kind=r8), private, save, allocatable :: mass(:)
   logical, private, save, allocatable :: sinup(:),sindn(:)
   logical, private, save :: usebf
   integer(kind=i4), private, save :: npdet
!  real(kind=r8), private, save :: a2
contains
   subroutine setupslat(rho,npartin,nbinin,massin,ndin,usebfin,ntab)
   use mympi
   use backflow
   integer(kind=i4) :: npartin,nbinin(:),nexup,nexdn,ndin,ntab
   logical :: usebfin
   integer(kind=i4),allocatable :: ikup(:,:),ikdn(:,:)
   real(kind=r8) :: rho,massin(:)
   integer(kind=i4) :: i,ic,j
   real(kind=r8) :: pi
   npart=npartin
   ndim=ndin
   usebf=usebfin
   allocate(nbin(2))
   nbin=nbinin
   allocate(mass(npart))
   mass(1:nbin(1))=massin(1)
   if (nbin(1).ne.npart) mass(nbin(1)+1:npart)=massin(2)
   if (myrank().eq.0) then
      read (5,*) rho ! number density
!     read (5,*) a2
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
      write (6,'(''rho ='',t30,f20.10)') rho ! number density
      el=(npart/rho)**(1.0_r8/ndim)
      write (6,'(''L/2 ='',t40,f10.5)') el*0.5_r8
!     write (6,'(''a constant in the importance function ='',(t40,e12.5))') a2
      write (6,'(''unpaired up states ='',t40,i10)') nexup
      if (ndim.eq.3) write (6,'(''unpaired k values ='',(t35,3i5))') (ikup(:,j),j=1,nexup)
      if (ndim.eq.2) write (6,'(''unpaired k values ='',(t40,2i5))') (ikup(:,j),j=1,nexup)
      write (6,'(''unpaired down states ='',t40,i10)') nexdn
      if (ndim.eq.3) write (6,'(''unpaired k values ='',(t35,3i5))') (ikdn(:,j),j=1,nexdn)
      if (ndim.eq.2) write (6,'(''unpaired k values ='',(t40,2i5))') (ikdn(:,j),j=1,nexdn)
   endif
   call bcast(rho)
!  call bcast(a2)
   call bcast(el)
   if (usebf) call setbf(ntab,el)
   call bcast(nexup)
   call bcast(nexdn)
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
   end subroutine setupslat

   subroutine slat(x,phil,is,dphi,d2phi)
   real(kind=r8) :: x(:,:),phil,dphi(:,:),d2phi
   real(kind=r8) :: detup,ddetup(ndim,npart),det2up
   real(kind=r8) :: detdn,ddetdn(ndim,npart),det2dn
   real(kind=r8) :: xbig(ndim,npart),del(npart,ndim,npart,ndim),del2(ndim,npart)
   integer(kind=i4) :: is,isup,isdn
   call addbf(x,xbig,del,del2)
   call slatmat(xbig(:,1:nbin(1)),0,nbin(1),detup,isup,ddetup,det2up,sinup  &
       ,akup,ak2up,del(1:nbin(1),:,:,:),del2(:,1:nbin(1)))
   call slatmat(xbig(:,1+nbin(1):npart),nbin(1),nbin(2),detdn,isdn,ddetdn,det2dn,sindn  &
       ,akdn,ak2dn,del(nbin(1)+1:npart,:,:,:),del2(:,nbin(1)+1:npart))
   phil=detup+detdn
   dphi=ddetup+ddetdn
   d2phi=det2up+det2dn
   if (usebf) d2phi=d2phi+2.0_r8*sum(ddetup(:,:)*ddetdn(:,:))
   is=isup*isdn
   end subroutine slat

   subroutine slatmat(x,np0,npar,det,is,ddet,d2det,sinp,ak,ak2,del,del2)
   use matrixmod
   use backflow
   integer(kind=i4) :: np0,npar,i,j,is,ic,jc
   real(kind=r8) :: x(:,:),d2det,akr,c,s
   real(kind=r8) :: det,ddet(:,:)
   real(kind=r8), dimension(npar,npar) :: smati
   real(kind=r8), dimension(ndim,npar,npar) :: dph
   real(kind=r8) :: ak(:,:),ak2(:)
   integer(kind=i4) :: k,kc
   real(kind=r8) :: del(:,:,:,:),dmat(ndim,npar,npar)
   real(kind=r8) :: del2(:,:),d2mat(ndim,ndim,npar)
   real(kind=r8), dimension(ndim,ndim,npar,npar) :: d2ph
   logical :: sinp(:)
   smati=0.0_r8
   dph=0.0_r8
   d2ph=0.0_r8
   if (npar.eq.0) then
      det=0.0_r8
      ddet=0.0_r8
      d2det=0.0_r8
      is=1
      return
   endif
   do i=1,npar
      do j=1,npar
         akr=sum(x(:,i)*ak(:,j))
         c=cos(akr)
         s=sin(akr)
         if (sinp(j)) then
            smati(i,j)=s
            dph(:,i,j)=ak(:,j)*c
            if (usebf) then
               do ic=1,ndim
                  do jc=1,ndim
                     d2ph(ic,jc,i,j)=-ak(ic,j)*ak(jc,j)*s
                  enddo
               enddo
             endif
         else
            smati(i,j)=c
            dph(:,i,j)=-ak(:,j)*s
            if (usebf) then
               do ic=1,ndim
                  do jc=1,ndim
                     d2ph(ic,jc,i,j)=-ak(ic,j)*ak(jc,j)*c
                  enddo
               enddo
            endif
         endif
      enddo
   enddo
   call matinv(smati,det,is,npar)
   ddet=0.0_r8
   d2det=0.0_r8
   if (usebf) then
      dmat=0.0_r8
      d2mat=0.0_r8
      do i=1,npar
         do ic=1,ndim
            do jc=1,ndim
               d2mat(ic,jc,i)=sum(smati(:,i)*d2ph(ic,jc,i,:))
            enddo
         enddo
         do j=1,npar
            do ic=1,ndim
               dmat(ic,i,j)=sum(smati(:,i)*dph(ic,j,:))
            enddo
         enddo
      enddo
      do i=1,npart
         do ic=1,ndim
            do j=1,npar
               do jc=1,ndim
                  ddet(ic,i)=ddet(ic,i)+dmat(jc,j,j)*del(j,jc,i,ic)/sqrt(mass(i))
               enddo
            enddo
         enddo
      enddo
      do i=1,npar
         d2det=d2det+sum(dmat(:,i,i)*del2(:,i))/mass(i+np0)
         do ic=1,ndim
            do jc=1,ndim
               d2det=d2det+d2mat(ic,jc,i)*sum(del(i,ic,:,:)*del(i,jc,:,:))/mass(i+np0)
            enddo
         enddo
      enddo
      do j=1,npar
         do k=1,npar
            if (j.ne.k) then
               do jc=1,ndim
                  do kc=1,ndim
                     d2det=d2det+ &
                        (dmat(jc,j,j)*dmat(kc,k,k)-dmat(jc,k,j)*dmat(kc,j,k)) &
                        *sum(del(k,kc,:,:)*del(j,jc,:,:))/mass(j+np0)
                  enddo
               enddo
            endif
         enddo
      enddo
   else
      do i=1,npar
         do ic=1,ndim
            ddet(ic,i+np0)=sum(smati(:,i)*dph(ic,i,:))/sqrt(mass(i+np0))
         enddo
      enddo
      d2det=-sum(ak2(:)/mass(np0+1:np0+npar))
   endif
   end subroutine slatmat

   subroutine addbf(x,xbig,del,del2)
   use backflow
   real(kind=r8) :: x(:,:),xbig(:,:),del(:,:,:,:),del2(:,:)
   real(kind=r8) :: r,dx(ndim),bf,dbf,d2bf
   integer(kind=i4) :: i,j,ic,jc
   xbig=x
   del=0.0_r8
   del2=0.0_r8
   if (.not.usebf) return
   do i=1,npart-1
      do j=i+1,npart
         dx=x(:,i)-x(:,j)
         dx=dx-el*nint(dx/el)
         r=sqrt(sum(dx**2))
         call getbf(r,bf,dbf,d2bf)
         xbig(:,i)=xbig(:,i)+bf*dx
         xbig(:,j)=xbig(:,j)-bf*dx
         do ic=1,ndim
            del(j,ic,i,ic)=del(j,ic,i,ic)-bf
            do jc=1,ndim
               del(j,jc,i,ic)=del(j,jc,i,ic)-dbf*dx(ic)*dx(jc)/r
               del(i,ic,j,jc)=del(j,jc,i,ic)
            enddo
         enddo
! check that the following coefficient is correct for any system
         del2(:,i)=del2(:,i)+2.0_r8*(d2bf+(ndim+1.0_r8)*dbf/r)*dx(:)
         del2(:,j)=del2(:,j)-2.0_r8*(d2bf+(ndim+1.0_r8)*dbf/r)*dx(:)
      enddo
   enddo
   do i=1,npart
      do j=1,npart
         if (j.ne.i) then
            del(i,:,i,:)=del(i,:,i,:)-del(j,:,i,:)
         endif
      enddo
   enddo
   do i=1,npart
      do ic=1,ndim
         del(i,ic,i,ic)=del(i,ic,i,ic)+1.0_r8
      enddo
   enddo
   end subroutine addbf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the following are subrotuines needed by the optimization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine setslatpar(params)
   use backflow
   real(kind=r8) :: params(:),pbf(4)
   if (usebf) then
      pbf=params
      call initbf(pbf)
   endif
   end subroutine setslatpar

   subroutine getslatpar(nparam,params)
   use backflow
   integer(kind=i4) :: nparam
   real(kind=r8), pointer :: params(:)
   real(kind=r8) :: pbf(4)
   if (usebf) then
      nparam=4
      allocate(params(nparam))
      call getbfpar(pbf)
      params=pbf
      npdet=nparam
   else
      nparam=0
   endif
   end subroutine getslatpar

   subroutine getderslat(x,dpsi)
   use matrixmod
   use backflow
   real(kind=r8) :: x(:,:),dpsi(:)
   real(kind=r8) :: psi1,psi2
   real(kind=r8), parameter :: dp=0.0001_r8
   real(kind=r8) :: pbf(4)
   real(kind=r8) :: phil,dum1(3,npart),dum2
   integer(kind=i4) :: dum,i
   dpsi=0.0_r8
   if (.not.usebf) return
   call getbfpar(pbf)
   do i=1,npdet
      pbf(i)=pbf(i)+dp
      call initbf(pbf)
      call slat(x,phil,dum,dum1,dum2)
      psi1=phil
      pbf(i)=pbf(i)-2.0_r8*dp
      call initbf(pbf)
      call slat(x,phil,dum,dum1,dum2)
      psi2=phil
      pbf(i)=pbf(i)+dp
      call initbf(pbf)
      dpsi(i)=(psi1-psi2)/(2.0_r8*dp)
   enddo
   call slat(x,phil,dum,dum1,dum2)
   end subroutine getderslat


end module slater
