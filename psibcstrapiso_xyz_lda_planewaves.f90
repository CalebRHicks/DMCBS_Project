module bcstrap
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
   logical, private, save, allocatable :: iopt(:)
   real(kind=r8), private, save, allocatable :: smatisave(:,:)
   real(kind=r8), private, save :: omega,b
   real(kind=r8), private, save :: afact,bfact,hbar,csi,efg,rmax,omegal
   real(kind=r8), private, save :: p1,p2,p3
contains
   subroutine setupbcstrap(rho,npartin,nbinin,massin,ndin,ntabin)
   use kshell
   use mympi
   use betafun
   use matrixmod
   integer(kind=i4) :: npartin,nbinin(:),ndin
   integer(kind=i4), allocatable :: ikup(:,:),ikdn(:,:)
   real(kind=r8) :: rho,massin(:)
   integer(kind=i4) :: i,ic,j,ntabin,is
   real(kind=r8) :: pi,el2
   real(kind=r8) :: vext,dvext,d2vext,mat(3,3),vec(3),dum
   rho=0.0_r8
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
      read (5,*) el2
      read (5,*) omega  ! external well frequency
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
      allocate(vorb(norb+2))
      do i=1,norb
         read (5,*) vorb(i) !amplitude
      enddo
      read (5,*) vorb(norb+1)
      read (5,*) vorb(norb+2)
      read (5,*) b ! parameter of one body confinement
      el=2.0_r8*el2
      write (6,'(''L/2 ='',t40,f10.5)') el*0.5_r8
      write (6,'(''omega external potential ='',t40,f10.5)') omega
      write (6,'(''unpaired up states ='',t40,i10)') nexup
      if (ndim.eq.3) write (6,'(''unpaired k values ='',(t35,3i5))') (ikup(:,j),j=1,nexup)
      if (ndim.eq.2) write (6,'(''unpaired k values ='',(t40,2i5))') (ikup(:,j),j=1,nexup)
      write (6,'(''unpaired down states ='',t40,i10)') nexdn
      if (ndim.eq.3) write (6,'(''unpaired k values ='',(t35,3i5))') (ikdn(:,j),j=1,nexdn)
      if (ndim.eq.2) write (6,'(''unpaired k values ='',(t40,2i5))') (ikdn(:,j),j=1,nexdn)
      write (6,'(''s-wave pairing parameters:'')')
      write (6,'(''number of functions: '',t40,i10)') norb
      write (6,'(''amplitudes: '',t40,f10.5)') (vorb(i),i=1,norb)
      write (6,'(''gamma1 = '',t40,f10.5)') vorb(norb+1)
      write (6,'(''gamma2 = '',t40,f10.5)') vorb(norb+2)
      write (6,'(''b parameter of Jastrow ='',f10.5)') b
   endif
   call bcast(el)
   call bcast(omega)
   call bcast(nexup)
   call bcast(nexdn)
   call bcast(norb)
   if (myrank().ne.0) allocate(vorb(norb+2))
   call bcast(vorb)
   call bcast(b)
   call setbeta(ndim,ntabin,el)
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
   hbar=0.5_r8
   csi=0.40_r8  ! FIX THIS
   efg=0.3_r8
   afact=1.0_r8/(hbar*csi)
   bfact=csi*efg*5.0_r8/3.0_r8
   omegal=omega
   rmax=sqrt(2.0_r8*bfact)/omegal
   if (myrank().eq.0) then
      write(6,'(''LDA csi ='',t40,f10.5)') csi
      write(6,'(''LDA E_fg ='',t40,f10.5)') efg
      write(6,'(''LDA gamma1 ='',t40,f10.5)') vorb(norb+1)
      write(6,'(''LDA gamma2 ='',t40,f10.5)') vorb(norb+2)
      write(6,'(''LDA rmax ='',t40,f10.5)') rmax
   endif
   rmax=0.67_r8*rmax
   vext=0.5_r8*(omegal*rmax)**2
   dvext=omegal**2*rmax
   d2vext=omegal**2
   vec(1)=sqrt(afact*(bfact-vext))
   vec(2)=-0.5_r8/vec(1)*afact*dvext
   vec(3)=-afact*d2vext/(2.0_r8*vec(1))-(afact*dvext)**2/(4.0_r8*vec(1)**3)
   mat(1,1)=1.0_r8/rmax
   mat(1,2)=1.0_r8/rmax**2
   mat(1,3)=1.0_r8/rmax**3
   mat(2,1)=-1.0_r8/rmax**2
   mat(2,2)=-2.0_r8/rmax**3
   mat(2,3)=-3.0_r8/rmax**4
   mat(3,1)=2.0_r8/rmax**3
   mat(3,2)=6.0_r8/rmax**4
   mat(3,3)=12.0_r8/rmax**5
   call rmatinv(mat,dum,is,3)
   vec=matmul(mat,vec)
   p1=vec(1)
   p2=vec(2)
   p3=vec(3)
   end subroutine setupbcstrap

   subroutine psibcstrap(x,phil,is,dphi,d2phi)
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
            dph(2,:,jj,i)=-akdn(:,j)*c
            d2ph(2,jj,i)=-ak2dn(j)*s
         else
            smati(jj,i)=c
            dph(2,:,jj,i)=akdn(:,j)*s
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
smatisave=smati
   end subroutine psibcstrap

   subroutine pairfn(x1,x2,pf,dpf,d2pf)
   use betafun
   real(kind=r8), dimension(:) :: x1,x2
   real(kind=r8) :: dpf(2,3)
   real(kind=r8), dimension(ndim) :: dx
   real(kind=r8) :: pf,d2pf(2),r,c,s,arg,v
   integer(kind=i4) :: i,j,ish
   real(kind=r8) :: f,df,d2f,gradf(2,3),delf(2)
   real(kind=r8) :: xcm(3),rcm,gradk(3),delk
   real(kind=r8) :: vext,dvext,d2vext,kf,dkf,d2kf
   real(kind=r8) :: g1,dg1,d2g1,gam1,gradg1(3),delg1
   real(kind=r8) :: g2,dg2,d2g2,gam2,gradg2(3),delg2
   real(kind=r8) :: g3,dg3,d2g3,gam3,gradg3(3),delg3
   real(kind=r8) :: phi,gradphi(2,3),delphi(2)
!integer(kind=i4) :: iii,iiii
!x1=0.0_r8
!x2=0.0_r8
!do iii=0,10
   xcm=0.5_r8*(x1+x2)
   rcm=sqrt(sum(xcm**2))
!rcm=iii*0.2_r8
!do iiii=0,100
!x1(1)=iiii*0.05
   if (rcm.lt.rmax) then
      vext=0.5_r8*(omegal*rcm)**2
      dvext=omegal**2*rcm
      d2vext=omegal**2
      kf=sqrt(afact*(bfact-vext))
      dkf=-0.5_r8/kf*afact*dvext
      d2kf=-afact*d2vext/(2.0_r8*kf)-(afact*dvext)**2/(4.0_r8*kf**3)
   else
      kf=p1/rcm+p2/rcm**2+p3/rcm**3
      dkf=-p1/rcm**2-2.0_r8*p2/rcm**3-3.0_r8*p3/rcm**4
      d2kf=2.0_r8*p1/rcm**3+6.0_r8*p2/rcm**4+12.0_r8*p3/rcm**5
   endif
   gradk=0.5_r8*xcm/rcm*dkf
   delk=0.25_r8*(d2kf+2.0_r8/rcm*dkf)
   phi=0.0_r8
   gradphi=0.0_r8
   delphi=0.0_r8
   dx=x1-x2
   if (norb.eq.0) return
   phi=vorb(1)
   i=1
   do ish=2,norb
      do j=1,numsh(ish)
         i=i+1
         arg=sum(ak(:,i)*dx(:))
         c=cos(kf*arg)
         s=sin(kf*arg)
         v=vorb(ish)*2.0_r8
         phi=phi+v*c
         gradphi(1,:)=gradphi(1,:)-v*s*(kf*ak(:,i)+gradk*arg)
         gradphi(2,:)=gradphi(2,:)-v*s*(-kf*ak(:,i)+gradk*arg)
         delphi(1)=delphi(1)-v*c*sum((gradk(:)*arg+kf*ak(:,i))**2)  &
                -v*s*(delk*arg+2.0_r8*sum(gradk*ak(:,i)))
         delphi(2)=delphi(2)-v*c*sum((gradk(:)*arg-kf*ak(:,i))**2)  &
                -v*s*(delk*arg-2.0_r8*sum(gradk*ak(:,i)))
      enddo
   enddo
   r=sqrt(sum(dx**2))
!  kf=1.0_r8 !!!
!  gradk=0.0_r8 !!!
!  delk=0.0_r8 !!!
   call radialfun(kf*r,f,df,d2f)
   if (kf.eq.0.0_r8) kf=1e-10_r8
   d2f=d2f-2.0_r8*df/(kf*r)
   gradf(1,:)=df*(kf*dx/r+r*gradk)
   gradf(2,:)=df*(-kf*dx/r+r*gradk)
   delf(1)=d2f*sum((dx*kf/r+r*gradk)**2)+df*(2.0_r8/r*kf+2.0_r8*sum(dx*gradk)/r+r*delk)
   delf(2)=d2f*sum((-dx*kf/r+r*gradk)**2)+df*(2.0_r8/r*kf-2.0_r8*sum(dx*gradk)/r+r*delk)

   gam1=vorb(norb+1)
   g1=exp(-gam1*rcm**2)
   dg1=-2.0_r8*rcm*gam1*g1
   gradg1=0.5_r8*dg1*xcm/rcm
   d2g1=-2.0_r8*gam1*g1+4.0_r8*(rcm*gam1)**2*g1
   delg1=0.25_r8*(d2g1+2.0_r8/rcm*dg1)
   gam2=vorb(norb+2)
   g2=1.0_r8*exp(-gam2*rcm**2)
   dg2=-2.0_r8*rcm*gam2*g2
   gradg2=0.5_r8*dg2*xcm/rcm
   d2g2=-2.0_r8*gam2*g2+4.0_r8*(rcm*gam2)**2*g2
   delg2=0.25_r8*(d2g2+2.0_r8/rcm*dg2)
   gam3=0.5_r8 ! I assume external trap has omega=1.0
   g3=exp(-gam3*rcm**2)
   dg3=-2.0_r8*rcm*gam3*g3
   gradg3=0.5_r8*dg3*xcm/rcm
   d2g3=-2.0_r8*gam3*g3+4.0_r8*(rcm*gam3)**2*g3
   delg3=0.25_r8*(d2g3+2.0_r8/rcm*dg3)

   pf=phi*g1*g2+f*g3*(1.0_r8-g2)
   dpf(1,:)=gradphi(1,:)*g1*g2+phi*gradg1*g2+phi*g1*gradg2   &
           +gradf(1,:)*g3*(1.0_r8-g2)+f*gradg3*(1.0_r8-g2)-f*g3*gradg2
   dpf(2,:)=gradphi(2,:)*g1*g2+phi*gradg1*g2+phi*g1*gradg2   &
           +gradf(2,:)*g3*(1.0_r8-g2)+f*gradg3*(1.0_r8-g2)-f*g3*gradg2
   d2pf(1)=delphi(1)*g1*g2+sum(gradphi(1,:)*gradg1(:))*g2+g1*sum(gradphi(1,:)*gradg2(:))  &
          +sum(gradphi(1,:)*gradg1(:))*g2+phi*delg1*g2+phi*sum(gradg1(:)*gradg2(:))    &
          +sum(gradphi(1,:)*gradg2(:))*g1+phi*sum(gradg1(:)*gradg2(:))+phi*g1*delg2    &
          +delf(1)*g3*(1.0_r8-g2)+sum(gradf(1,:)*gradg3(:))*(1.0_r8-g2)-g3*sum(gradf(1,:)*gradg2(:)) &
          +sum(gradf(1,:)*gradg3(:))*(1.0_r8-g2)+f*delg3*(1.0_r8-g2)-f*sum(gradg3(:)*gradg2(:))  &
          -g3*sum(gradf(1,:)*gradg2(:))-f*sum(gradg3(:)*gradg2(:))-f*g3*delg2
   d2pf(2)=delphi(2)*g1*g2+sum(gradphi(2,:)*gradg1(:))*g2+g1*sum(gradphi(2,:)*gradg2(:))  &
          +sum(gradphi(2,:)*gradg1(:))*g2+phi*delg1*g2+phi*sum(gradg1(:)*gradg2(:))    &
          +sum(gradphi(2,:)*gradg2(:))*g1+phi*sum(gradg1(:)*gradg2(:))+phi*g1*delg2    &
          +delf(2)*g3*(1.0_r8-g2)+sum(gradf(2,:)*gradg3(:))*(1.0_r8-g2)-g3*sum(gradf(2,:)*gradg2(:)) &
          +sum(gradf(2,:)*gradg3(:))*(1.0_r8-g2)+f*delg3*(1.0_r8-g2)-f*sum(gradg3(:)*gradg2(:))  &
          -g3*sum(gradf(2,:)*gradg2(:))-f*sum(gradg3(:)*gradg2(:))-f*g3*delg2

!  pf=phi+f
!  dpf(1,:)=gradphi(1,:)+gradf(1,:)
!  dpf(2,:)=gradphi(2,:)+gradf(2,:)
!  d2pf(1)=delphi(1)+delf(1)
!  d2pf(2)=delphi(2)+delf(2)


!write(99,*) rcm,r,pf
!enddo
!write(99,*) ' '
!write(99,*) ' '
!enddo
!stop
   end subroutine pairfn



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

   subroutine setderbcstrap(nparmdet,ioptin)
   integer(kind=i4) :: nparmdet
   logical :: ioptin(:)
   npdet=nparmdet
   allocate(iopt(size(ioptin)))
   iopt=ioptin
   end subroutine setderbcstrap

   subroutine getderbcstrap(x,dpsi)
   use matrixmod
   use betafun
   real(kind=r8) :: x(:,:),dpsi(:)
   real(kind=r8) :: dx(ndim),arg,phil
   real(kind=r8), dimension(npair+nexup+nexdn,npair+nexup+nexdn,npdet) :: dphi
   real(kind=r8), dimension(npair+nexup+nexdn,npair+nexup+nexdn) :: smati,d2ph
   real(kind=r8), dimension(ndim,npair+nexup+nexdn,npair+nexup+nexdn) :: dph
   real(kind=r8) :: r,dbet(2)
   integer(kind=i4) :: i,j,jj,ii,k,kk,nd,is,np,ish
   real(kind=r8) :: xcm(3),rcm,kf,vext,g1,g2,g3,gam,pairf
   real(kind=r8) :: f,dum1,dum2
   dpsi=0.0_r8
   nd=npair+nexup+nexdn
   smati=0.0_r8
   dph=0.0_r8
   d2ph=0.0_r8
   do i=1,npair+nexup
      do j=1,npair+nexdn
         jj=j+npair+nexup
         dx(:)=x(:,i)-x(:,jj)
         call pairfn(x(:,i),x(:,jj),smati(i,j),dph(:,i,j),d2ph(i,j))
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
         dx(:)=x(:,i)-x(:,jj)
         r=sqrt(sum(dx**2))
         xcm=0.5_r8*(x(:,i)+x(:,jj))
         rcm=sqrt(sum(xcm**2))
         if (rcm.lt.rmax) then
            vext=0.5_r8*(omegal*rcm)**2
            kf=sqrt(afact*(bfact-vext))
         else
            kf=p1/rcm+p2/rcm**2+p3/rcm**3
         endif
         gam=vorb(norb+1)
         g1=exp(-gam*rcm**2)
         gam=vorb(norb+2)
         g2=1.0_r8*exp(-gam*rcm**2)
         gam=0.5_r8 ! I assume external trap has omega=1.0
         g3=exp(-gam*rcm**2)
         np=1
         if (iopt(1).or.iopt(2)) then
            call derbeta(r,dbet,iopt)
            if (iopt(1)) then
               dphi(i,j,np)=dbet(1)
               np=np+1
            endif
            if (iopt(2)) then
               dphi(i,j,np)=dbet(2)
               np=np+1
            endif
         endif
         if (iopt(3)) then
            dphi(i,j,np)=1.0_r8*g1*g2
            pairf=vorb(1)
            np=np+1
            k=1
            do ish=2,norb
               do kk=1,numsh(ish)
                  k=k+1
                  arg=sum(ak(:,k)*dx(:))
                  dphi(i,j,np)=dphi(i,j,np)+2.0_r8*cos(kf*arg)*g1*g2
                  pairf=pairf+2.0_r8*vorb(ish)*cos(kf*arg)
               enddo
               np=np+1
            enddo
            call radialfun(kf*r,f,dum1,dum2)
            dphi(i,j,np)=-pairf*rcm**2*g1*g2
            np=np+1
            dphi(i,j,np)=-pairf*rcm**2*g1*g2+f*g3*rcm**2*g2
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
   nsw=norb+2
   allocate(vsw(nsw))
   vsw=vorb
   end subroutine getbcstrapampl
end module bcstrap
