module bcstrap
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, save :: npart,ntab,ndim
   integer(kind=i4), private, save :: npair
   integer(kind=i4), private, save, allocatable :: nbin(:)
   integer(kind=i4), private, save :: nexup,nexdn,norb
   integer(kind=i4), private, save, pointer :: numsh(:),nquant(:,:)
   real(kind=r8), private, save, allocatable :: mass(:),vorb(:)
   integer(kind=i4), private, save, allocatable :: spin(:)
   integer(kind=i4), private, save :: npdet
   logical, private, save, allocatable :: iopt(:)
   real(kind=r8), private, save, allocatable :: orb(:,:,:),dorb(:,:,:),d2orb(:,:,:)
   real(kind=r8), private, save :: omega,b,bex
   integer(kind=i4), private, save, allocatable :: nkup(:,:),nkdn(:,:)
   integer(kind=i4), private, save :: ntabf
   real(kind=r8), private, save :: rangef,scalef
   real(kind=r8), private, save, allocatable :: ftab(:),dftab(:),d2ftab(:)
   integer(kind=i4), private, save :: radpair,lmax
   integer(kind=i4), allocatable :: nsh(:),lsh(:)
contains
   subroutine setupbcstrap(rho,npartin,nbinin,massin,ndin,ntabin)
   use mympi
   use betafun
   integer(kind=i4) :: npartin,nbinin(:),ndin,ntabin
   real(kind=r8) :: rho,massin(:),r
   integer(kind=i4) :: i,j
   rho=0.0_r8
   npart=npartin
   ndim=ndin
   ntab=ntabin
   allocate(nbin(2))
   nbin=nbinin
   if (massin(1).ne.1.0_r8.or.massin(2).ne.1.0_r8) then
      if (myrank().eq.0) write (6,'(''bcstrap works with just mass=1 by now'')')
      stop
   endif
   allocate(mass(npart))
   mass(1:nbin(1))=massin(1)
   mass(nbin(1)+1:npart)=massin(2)
   allocate(spin(npart))
! assign spin 1 to up particles and -1 to down ones
   spin(1:nbin(1))=1
   spin(nbin(1)+1:nbin(1)+nbin(2))=-1
   if (myrank().eq.0) then
      read (5,*) omega  ! external well frequency
      read (5,*) nexup
      allocate(nkup(2,nexup))
      do j=1,nexup
         read (5,*) nkup(:,j)
      enddo
      read (5,*) nexdn
      allocate(nkdn(2,nexdn))
      do j=1,nexdn
         read (5,*) nkdn(:,j)
      enddo
      read (5,*) bex ! parameter of one body unpaired states
      read (5,*) norb  !number of orbitals for s-wave
      allocate(vorb(norb),nsh(norb),lsh(norb))
      do i=1,norb
         read (5,*) nsh(i),lsh(i),vorb(i) !amplitude
      enddo
      read (5,*) radpair ! 0=no, 1=betafun, 2=fort.80
      read (5,*) b ! parameter of one body confinement
      if (radpair.eq.1) read (5,*) rangef ! range of betafun
      write (6,'(''omega external potential ='',t40,f10.5)') omega
      write (6,'(''unpaired up states ='',t40,i10)') nexup
      write (6,'(''quantum numbers ='',(t35,3i5))') (nkup(:,j),j=1,nexup)
      write (6,'(''unpaired down states ='',t40,i10)') nexdn
      write (6,'(''quantum numbers ='',(t35,3i5))') (nkdn(:,j),j=1,nexdn)
      write (6,'(''bex of unpaired states ='',t40,f10.5)') bex
      write (6,'(''s-wave pairing parameters:'')')
      write (6,'(''number of functions: '',t40,i10)') norb
      do i=1,norb
         write (6,'(''(n,l) and amplitude '',t20,2i5,f20.10)') nsh(i),lsh(i),vorb(i)
      enddo
      if (radpair.eq.0) write (6,'(''no radial function in pairing'')')
      if (radpair.eq.1) write (6,'(''using beta function with range =''t40,f10.5)') rangef
      if (radpair.eq.2) write (6,'(''using fort.80 in pairing'')')
      write (6,'(''b parameter of Jastrow ='',f10.5)') b
   endif
   call bcast(omega)
   call bcast(nexup)
   call bcast(nexdn)
   if (myrank().ne.0) allocate(nkup(2,nexup),nkdn(2,nexdn))
   call bcast(nkup)
   call bcast(nkdn)
   call bcast(bex)
   call bcast(norb)
   if (myrank().ne.0) allocate(vorb(norb),nsh(norb),lsh(norb))
   call bcast(vorb)
   call bcast(nsh)
   call bcast(lsh)
   call bcast(radpair)
   call bcast(b)
   if (radpair.eq.2) then
      if (myrank().eq.0) then
         rewind 80
         read (80,*) ntabf
         ntabf=ntabf-1
         allocate(ftab(0:ntabf),dftab(0:ntabf),d2ftab(0:ntabf))
         do i=1,ntabf
            read (80,('(4e25.15)')) r,ftab(i),dftab(i),d2ftab(i)
         enddo
         rangef=r
         scalef=real(ntabf)/rangef
      endif
      call bcast(ntabf)
      call bcast(rangef)
      call bcast(scalef)
      if (myrank().ne.0) allocate(ftab(0:ntabf),dftab(0:ntabf),d2ftab(0:ntabf))
      call bcast(ftab)
      call bcast(dftab)
      call bcast(d2ftab)
   else if (radpair.eq.1) then
      call bcast(rangef)
      call setbeta(ndim,ntab,2.0_r8*rangef)
   endif
   lmax=maxval(lsh)
   npair=nbin(1)+nbin(2)-nexup-nexdn
   if (mod(npair,2).ne.0) then
      write (6,'(''inconsistent pairing number'',3i5)') &
         nbin(1)+nbin(2),nexup,nexdn
      stop
   endif
   npair=npair/2
   end subroutine setupbcstrap
 
   subroutine radialfun2(r,f,df,d2f)
   real(kind=r8) :: r,f,df,d2f
   real(kind=r8) :: ddr,a1,a2,a3,a4
   integer(kind=i4) :: index
   if (r.gt.rangef) then
      f=0.0_r8
      df=0.0_r8
      d2f=0.0_r8
      return
   endif
   ddr=scalef*r
   index=ddr
   index=max(1,min(index,ntabf-2))
   ddr=ddr-index
   a1=-ddr*(ddr-1.0_r8)*(ddr-2.0_r8)/6.0_r8
   a2=(ddr+1.0_r8)*(ddr-1.0_r8)*(ddr-2.0_r8)/2.0_r8
   a3=-(ddr+1.0_r8)*ddr*(ddr-2.0_r8)/2.0_r8
   a4=(ddr+1.0_r8)*ddr*(ddr-1.0_r8)/6.0_r8
   f=a1*ftab(index-1)+a2*ftab(index)+a3*ftab(index+1)+a4*ftab(index+2)
   df=a1*dftab(index-1)+a2*dftab(index)+a3*dftab(index+1)+a4*dftab(index+2)
   d2f=a1*d2ftab(index-1)+a2*d2ftab(index)+a3*d2ftab(index+1)+a4*d2ftab(index+2)
   end subroutine radialfun2

   subroutine pairfn(x1,x2,pf,dpf,d2pf)
   use betafun
   use ylmmod
   real(kind=r8) :: x1(:),x2(:),pf,dpf(:,:),d2pf(:)
   complex(kind=r8) :: f1,df1(3),d2f1,f2,df2(3),d2f2
   complex(kind=r8) :: ylm1(0:(lmax+1)*(lmax+1)-1)
   complex(kind=r8) :: ylm2(0:(lmax+1)*(lmax+1)-1)
   complex(kind=r8) :: dylm1(0:2,0:(lmax+1)*(lmax+1)-1)
   complex(kind=r8) :: dylm2(0:2,0:(lmax+1)*(lmax+1)-1)
   integer(kind=i4) :: ish,i,m,lm1,lm2
   real(kind=r8) :: dx(ndim),r,f,df,d2f,r2a,r2b,rr1,rr2
   real(kind=r8) :: alp,cg,ex,fact,v
   real(kind=r8) :: lag,dlag,d2lag
   real(kind=r8) :: p1,dp1,d2p1
   real(kind=r8) :: p2,dp2,d2p2
   alp=0.5_r8*omega
   pf=0.0_r8
   dpf=0.0_r8
   d2pf=0.0_r8
   i=0
! compute all the spherical harmonics and store them
   call ylmcal(x1,lmax,ylm1,dylm1)
   call ylmcal(x2,lmax,ylm2,dylm2)
   r2a=sum(x1**2)
   rr1=sqrt(r2a)
   r2b=sum(x2**2)
   rr2=sqrt(r2b)
   do ish=1,norb
! get generalized Laguerre pol. L_n^(l+1/2)(2*mu*r**2)
! where n and l are quantum numbers, mu=m*omega/(2*hbar)
! this works only for isotropic harmonic oscillator.
! I assume hbar=m=1, fix it!
      call laguerre(nsh(ish),lsh(ish)+0.5_r8,2.0_r8*alp*r2a,lag)
      call laguerre(nsh(ish)-1,lsh(ish)+0.5_r8+1.0_r8,2.0_r8*alp*r2a,dlag)
      dlag=-dlag
      call laguerre(nsh(ish)-2,lsh(ish)+0.5_r8+2,2.0_r8*alp*r2a,d2lag)
      ex=exp(-alp*r2a)
      p1=lag*ex
      dp1=2.0_r8*alp*rr1*ex*(2.0_r8*dlag-lag)
      d2p1=2.0_r8*alp*ex*((2.0_r8*alp*r2a-1.0_r8)*lag+(2.0_r8-8.0_r8*alp*r2a)*dlag+8.0_r8*alp*r2a*d2lag)
! do the same for the second particle particle
      call laguerre(nsh(ish),lsh(ish)+0.5_r8,2.0_r8*alp*r2b,lag)
      call laguerre(nsh(ish)-1,lsh(ish)+0.5_r8+1.0_r8,2.0_r8*alp*r2b,dlag)
      dlag=-dlag
      call laguerre(nsh(ish)-2,lsh(ish)+0.5_r8+2,2.0_r8*alp*r2b,d2lag)
      ex=exp(-alp*r2b)
      p2=lag*ex
      dp2=2.0_r8*alp*rr2*ex*(2.0_r8*dlag-lag)
      d2p2=2.0_r8*alp*ex*((2.0_r8*alp*r2b-1.0_r8)*lag+(2.0_r8-8.0_r8*alp*r2b)*dlag+8.0_r8*alp*r2b*d2lag)
      fact=1.0_r8/sqrt(2.0_r8*lsh(ish)+1.0_r8)
      do m=-lsh(ish),lsh(ish)
         cg=(-1)**(lsh(ish)+m)*fact
         lm1=lsh(ish)*(lsh(ish)+1)+m
         lm2=lsh(ish)*(lsh(ish)+1)-m
         f1=p1*ylm1(lm1)
         df1(:)=dp1*x1(:)/rr1*ylm1(lm1)+p1*dylm1(0:2,lm1)
         d2f1=(d2p1+2.0_r8*dp1/rr1)*ylm1(lm1)+2.0_r8*dot_product(dp1*x1(:)/rr1,dylm1(0:2,lm1))
         f2=p2*ylm2(lm2)
         df2(:)=dp2*x2(:)/rr2*ylm2(lm2)+p2*dylm2(0:2,lm2)
         d2f2=(d2p2+2.0_r8*dp2/rr2)*ylm2(lm2)+2.0_r8*dot_product(dp2*x2(:)/rr2,dylm2(0:2,lm2))
         v=vorb(ish)*cg
         pf=pf+v*f1*f2
         dpf(1,:)=dpf(1,:)+v*df1(:)*f2
         dpf(2,:)=dpf(2,:)+v*f1*df2(:)
         d2pf(1)=d2pf(1)+v*d2f1*f2
         d2pf(2)=d2pf(2)+v*f1*d2f2
      enddo
   enddo
   if (radpair.eq.0) return
   dx=x1-x2
   r=sqrt(sum(dx**2))
   if (radpair.eq.1) then
      call radialfun(r,f,df,d2f)
      d2pf=d2pf+d2f
   else 
      call radialfun2(r,f,df,d2f)
      d2pf=d2pf+d2f+(ndim-1.0_r8)*df/r
   endif
   pf=pf+f
   dpf(1,:)=dpf(1,:)+df*dx/r
   dpf(2,:)=dpf(2,:)-df*dx/r
   end subroutine pairfn

   subroutine psibcstrap(x,phil,is,dphi,d2phi)
   use matrixmod
   real(kind=r8), dimension(ndim,npart) :: x,dphi
   real(kind=r8) :: d2phi
   real(kind=r8) :: smati(npair+nexup+nexdn,npair+nexup+nexdn)
   real(kind=r8) :: dph(2,ndim,npair+nexup+nexdn,npair+nexup+nexdn)
   real(kind=r8) :: d2ph(2,npair+nexup+nexdn,npair+nexup+nexdn)
   real(kind=r8) :: phil
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
         call getunpaired(nkup(:,j),x(:,i),smati(i,jj),dph(1,:,i,jj),d2ph(1,i,jj))
      enddo
   enddo
   do i=1,npair+nexdn
      ii=i+npair+nexup
      do j=1,nexdn
         jj=j+npair+nexup
         call getunpaired(nkdn(:,j),x(:,ii),smati(jj,i),dph(2,:,jj,i),d2ph(2,jj,i))
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
   end subroutine psibcstrap

   subroutine onebody(x,g,dg,d2g,vext)
! this subroutine works with hbar=1, fix it
   real(kind=r8) :: x(:,:),g,dg(:,:),d2g,vext,r2
   integer(kind=i4) :: i
   g=0.0_r8
   dg=0.0_r8
   d2g=0.0_r8
   do i=1,npart
      r2=sum(x(:,i)**2)
      g=g-r2
      dg(:,i)=-2.0_r8*b*x(:,i)/sqrt(mass(i))
      d2g=d2g-(6.0_r8-4.0_r8*b*r2)*b/mass(i)
   enddo
   g=b*g
   vext=0.0_r8
   do i=1,npart
      vext=vext+0.5_r8*mass(i)*sum((omega*x(:,i))**2)
   enddo
   end subroutine onebody

   subroutine getunpaired(nquant,x,phi,dphi,d2phi)
   integer(kind=i4) :: nquant(:),n,l
   real(kind=r8) :: x(:),phi,dphi(:),d2phi
   real(kind=r8) :: ylm,dylm(3)
   real(kind=r8) :: lag,dlag,d2lag
   real(kind=r8) :: p,dp,d2p,ex
   real(kind=r8) :: alp,x2,y2,z2,r,r2
! nquant is a vector with (n,l,m)
! here I suppose to have only one unpaired state
! then I pickup just one value of m, I don't care what
! this is ugly but quickly, fix it!
   l=nquant(2)
   n=nquant(1)
   select case (l)
   case (0) ! l=0
      ylm=1.0_r8
      dylm=0.0_r8
   case (1) ! l=1
      ylm=x(1)
      dylm(1)=1.0_r8
   case (2) ! l=2
      ylm=2.0_r8*x(3)**2-x(1)**2-x(2)**2
      dylm(1)=-2.0_r8*x(1)
      dylm(2)=-2.0_r8*x(2)
      dylm(3)=4.0_r8*x(3)
   case (3) ! l=3
      ylm=x(1)*x(2)*x(3)
      dylm(1)=x(2)*x(3)
      dylm(2)=x(1)*x(3)
      dylm(3)=x(1)*x(2)
   case (4) ! l=4
      x2=x(1)*x(1)
      y2=x(2)*x(2)
      z2=x(3)*x(3)
      ylm=8.0_r8*z2*z2-24.0_r8*x2*z2-24.0_r8*y2*z2+3.0_r8*x2*x2+6.0_r8*x2*y2+3.0_r8*y2*y2
      dylm(1)=-60.0_r8*x(1)*z2+12.0_r8*(x2+y2+z2)*x(1)
      dylm(2)=-60.0_r8*x(2)*z2+12.0_r8*(x2+y2+z2)*x(2)
      dylm(3)=80.0_r8*z2*x(3)-48.0_r8*x(3)*(x2+y2+z2)
   case default
      write (6,*) 'I can use up to l=4!'
      stop
   end select
   r2=sum(x**2)
   r=sqrt(r2)
   alp=0.5_r8*omega*bex
! get generalized Laguerre pol. L_n^(l+1/2)(2*mu*r**2)
! where n and l are quantum numbers, mu=m*omega/(2*hbar)
! this works only for isotropic harmonic oscillator.
! I assume hbar=m=1, fix it!
   call laguerre(n,l+0.5_r8,2.0_r8*alp*r2,lag)
   call laguerre(n-1,l+0.5_r8+1.0_r8,2.0_r8*alp*r2,dlag)
   dlag=-dlag
   call laguerre(n-2,l+0.5_r8+2,2.0_r8*alp*r2,d2lag)
   ex=exp(-alp*r2)
   p=lag*ex
   dp=2.0_r8*alp*r*ex*(2.0_r8*dlag-lag)
   d2p=2.0_r8*alp*ex*((2.0_r8*alp*r2-1.0_r8)*lag+(2.0_r8-8.0_r8*alp*r2)*dlag+8.0_r8*alp*r2*d2lag)
   phi=p*ylm
   dphi(:)=dp*x(:)/r*ylm+p*dylm(:)
   d2phi=(d2p+2.0_r8*dp/r)*ylm+2.0_r8*dot_product(dp*x(:)/r,dylm)
   end subroutine getunpaired

   subroutine laguerre(n,a,x,lag)
   integer(kind=i4) :: n,i
   real(kind=r8) :: a,x,lag,l0,l1
   if (n.lt.0) then
      lag=0.0_r8
      return
   else if (n.eq.0) then
      lag=1.0_r8
      return
   endif
   l1=0.0_r8
   lag=1.0_r8
   do i=1,n
      l0=l1
      l1=lag
      lag=((2*i-1.0_r8+a-x)*l1-(i-1.0_r8+a)*l0)/i
   enddo
   return
   end subroutine laguerre


   subroutine setderbcstrap(nparmdet,ioptin)
   integer(kind=i4) :: nparmdet
   logical :: ioptin(:)
   npdet=nparmdet
   allocate(iopt(size(ioptin)))
   iopt=ioptin
   end subroutine setderbcstrap

   subroutine getderbcstrap(x,dpsi)
   use matrixmod
   use ylmmod
   real(kind=r8) :: x(:,:),dpsi(:)
   real(kind=r8) :: phil
   real(kind=r8), dimension(npair+nexup+nexdn,npair+nexup+nexdn,npdet) :: dphi
   real(kind=r8), dimension(npair+nexup+nexdn,npair+nexup+nexdn) :: smati
   real(kind=r8), dimension(2,npair+nexup+nexdn,npair+nexup+nexdn) :: d2ph
   real(kind=r8), dimension(2,ndim,npair+nexup+nexdn,npair+nexup+nexdn) :: dph
   complex(kind=r8) :: f1,f2
   real(kind=r8) :: r2,rr1,rr2,alp,fact,cg,lag,p1,p2
   integer(kind=i4) :: i,j,jj,ii,nd,is,np,ish,m,lm1,lm2
   complex(kind=r8) :: ylm1(0:(lmax+1)*(lmax+1)-1)
   complex(kind=r8) :: ylm2(0:(lmax+1)*(lmax+1)-1)
   complex(kind=r8) :: dylm1(0:2,0:(lmax+1)*(lmax+1)-1)
   complex(kind=r8) :: dylm2(0:2,0:(lmax+1)*(lmax+1)-1)
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
         call getunpaired(nkup(:,j),x(:,i),smati(i,jj),dph(1,:,i,jj),d2ph(1,i,jj))
      enddo
   enddo
   do i=1,npair+nexdn
      ii=i+npair+nexup
      do j=1,nexdn
         jj=j+npair+nexup
         call getunpaired(nkdn(:,j),x(:,ii),smati(jj,i),dph(2,:,jj,i),d2ph(2,jj,i))
      enddo
   enddo
   call matinv(smati,phil,is,nd)
   dphi=0.0_r8
   alp=0.5_r8*omega
   do i=1,npair+nexup
      do j=1,npair+nexdn
         jj=j+npair+nexup
         np=1
         if (iopt(3)) then
            call ylmcal(x(:,i),lmax,ylm1,dylm1)
            call ylmcal(x(:,jj),lmax,ylm2,dylm2)
            do ish=1,norb
               r2=sum(x(:,i)**2)
               rr1=sqrt(r2)
               call laguerre(nsh(ish),lsh(ish)+0.5_r8,2.0_r8*alp*r2,lag)
               p1=lag*exp(-alp*r2)
               r2=sum(x(:,jj)**2)
               rr2=sqrt(r2)
               call laguerre(nsh(ish),lsh(ish)+0.5_r8,2.0_r8*alp*r2,lag)
               p2=lag*exp(-alp*r2)
               fact=1.0_r8/sqrt(2.0_r8*lsh(ish)+1.0_r8)
               do m=-lsh(ish),lsh(ish)
                  cg=(-1)**(lsh(ish)+m)*fact
                  lm1=lsh(ish)*(lsh(ish)+1)+m
                  lm2=lsh(ish)*(lsh(ish)+1)-m
                  f1=p1*ylm1(lm1)
                  f2=p2*ylm2(lm2)
                  dphi(i,j,np)=dphi(i,j,np)+cg*f1*f2
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
   end subroutine getderbcstrap

   subroutine setvorbbcstrap(v0)
   real(kind=r8) :: v0(:)
   vorb=v0
   end subroutine setvorbbcstrap

   subroutine getbcstrapampl(nsw,vsw)
   integer(kind=i4) :: nsw
   real(kind=r8), pointer :: vsw(:)
   nsw=norb
   allocate(vsw(nsw))
   vsw=vorb
   end subroutine getbcstrapampl
end module bcstrap
