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
   real(kind=r8), private, save :: omega,b,bex
   integer(kind=i4), private, save, allocatable :: nkup(:,:),nkdn(:,:)
   integer(kind=i4), private, save :: ntabf
   real(kind=r8), private, save :: rangef,scalef
   real(kind=r8), private, save, allocatable :: ftab(:),dftab(:),d2ftab(:)
   integer(kind=i4), private, save :: radpair
   real(kind=r8), private, save :: afact,bfact,hbar,csi,efg,rmax,omegal
   real(kind=r8), private, save :: p1,p2,p3
contains
   subroutine setupbcstrap(rho,npartin,nbinin,massin,ndin,ntabin)
   use mympi
   use betafun
   use matrixmod
   integer(kind=i4) :: npartin,nbinin(:),ndin,ntabin
   real(kind=r8) :: rho,massin(:),r
   real(kind=r8), pointer :: ene(:)
   integer(kind=i4) :: i,j,k,is
   real(kind=r8) :: vext,dvext,d2vext,mat(3,3),vec(3),dum
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
      allocate(vorb(norb+3))
      do i=1,norb
         read (5,*) vorb(i) !amplitude
      enddo
      read (5,*) vorb(norb+1)
      read (5,*) vorb(norb+2)
      read (5,*) vorb(norb+3)
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
   if (myrank().ne.0) allocate(vorb(norb+3))
   call bcast(vorb)
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
   call sethobasis(omega,norb,numsh,nquant,ene)
   if (myrank().eq.0) then
      k=0
      do i=1,norb
         write (6,'(''shell number, energy and amplitude'',t30,i10,f10.5,f20.10)') i,ene(k+1),vorb(i)
         do j=1,numsh(i)
            k=k+1
            write (6,'(''state '',t35,3i5)') nquant(:,k)
         enddo
      enddo
      write (6,'(''gamma1 = '',t40,f10.5)') vorb(norb+1)
      write (6,'(''gamma2 = '',t40,f10.5)') vorb(norb+2)
      write (6,'(''multiply omega by this in orbitals = '',t40,f10.5)') vorb(norb+3)
   endif
   npair=nbin(1)+nbin(2)-nexup-nexdn
   if (mod(npair,2).ne.0) then
      write (6,'(''inconsistent pairing number'',3i5)') &
         nbin(1)+nbin(2),nexup,nexdn
      stop
   endif
   npair=npair/2
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
 
   subroutine sethobasis(om,nsh,degen,ik,ene)
   real(kind=r8) :: om
   integer(kind=i4) :: nsh
   integer(kind=i4), pointer :: ik(:,:),degen(:)
   integer(kind=i4), allocatable :: iktmp(:,:)
   integer(kind=i4) :: nmax
   integer(kind=i4) :: nx,ny,nz,i,nk,ntmp(3),j,ish
   real(kind=r8) :: tmp
   real(kind=r8), pointer :: ene(:)
   nmax=nsh/2
   if (nsh.eq.1) nmax=nsh
   if (nsh.gt.15) nmax=nsh/5
   allocate(iktmp(3,(nmax+1)**3))
   iktmp=0
   nk=0
   do nx=0,nmax
      do ny=0,nmax
         do nz=0,nmax
            nk=nk+1
            iktmp(1,nk)=nx
            iktmp(2,nk)=ny
            iktmp(3,nk)=nz
         enddo
      enddo
   enddo
   allocate(ene(nk))
   do i=1,nk-1
      do j=i+1,nk
         ene(i)=(0.5_r8+iktmp(1,i))*om+(0.5_r8+iktmp(2,i))*om+(0.5_r8+iktmp(3,i))*om
         ene(j)=(0.5_r8+iktmp(1,j))*om+(0.5_r8+iktmp(2,j))*om+(0.5_r8+iktmp(3,j))*om
         if (ene(i).gt.ene(j)) then
            ntmp=iktmp(:,i)
            iktmp(:,i)=iktmp(:,j)
            iktmp(:,j)=ntmp
            tmp=ene(i)
            ene(i)=ene(j)
            ene(j)=tmp
         endif
      enddo
   enddo
   do i=1,nk-1
      do j=i+1,nk
         if (ene(i).eq.ene(j).and.sum(iktmp(:,i)**2).gt.sum(iktmp(:,j)**2)) then
            ntmp=iktmp(:,i)
            iktmp(:,i)=iktmp(:,j)
            iktmp(:,j)=ntmp
            tmp=ene(i)
            ene(i)=ene(j)
            ene(j)=tmp
         endif
      enddo
   enddo
   allocate(degen(nsh))
   degen=1
   ish=1
   do i=2,nk
      if (ene(i).eq.ene(i-1).and.sum(iktmp(:,i)**2).eq.sum(iktmp(:,i-1)**2)) then
         degen(ish)=degen(ish)+1
      else
         degen(ish)=degen(ish)
         ish=ish+1
      endif
      if (ish.gt.nsh) exit
   enddo
   if (nk.le.sum(degen)) write (6,*) 'error in hobasis, you should try to increase nmax'
   allocate(ik(3,sum(degen)))
   ik(:,1:sum(degen))=iktmp(:,1:sum(degen))
   deallocate(iktmp)
   end subroutine sethobasis

   subroutine hoorb(n,om,x,phi,dphi,d2phi)
   integer(kind=i4) :: n
   real(kind=r8) :: x,phi,dphi,d2phi
   real(kind=r8) :: arg,om,h,dh,d2h,ex
   arg=sqrt(om)*x
   call hermite(n,arg,h,dh,d2h)
   dh=dh*sqrt(om)
   d2h=d2h*om
   ex=exp(-0.5_r8*om*x*x)
   phi=h*ex
   dphi=(dh-x*om*h)*ex
   d2phi=(d2h-om*(2.0_r8*x*dh+h-om*x*x*h))*ex
   end subroutine hoorb

   subroutine hermite(n,x,h,dh,d2h)
   real(kind=r8) :: x,h,dh,d2h
   integer(kind=i4) :: n,k
   real(kind=r8) :: dn,hm1,hm2,dnn,ypm
   h=1.0_r8
   dh=0.0_r8
   d2h=0.0_r8
   if(n.eq.0) return
   h=x*2.0_r8
   dh=2.0_r8
   d2h=0.0_r8
   if(n.eq.1) return
   hm2=1.0_r8
   do k=2,n
      hm1=h
      h=x*2.0_r8*hm1-2.0_r8*(k-1)*hm2
      ypm=hm2
      hm2=hm1
   enddo
   dn=n*2.0_r8
   dnn=2.0_r8*(n-1)
   dh=dn*hm2
   d2h=dn*dnn*ypm
   end subroutine hermite

   subroutine getorb(n,om,x,phi,dphi,d2phi)
   integer(kind=i4) :: n(:),ic
   real(kind=r8) :: x(:),phi,dphi(:),d2phi
   real(kind=r8) :: f(3),df(3),d2f(3),om
   do ic=1,3
      call hoorb(n(ic),om,x(ic),f(ic),df(ic),d2f(ic))
   enddo
   phi=f(1)*f(2)*f(3)
   dphi(1)=df(1)*f(2)*f(3)
   dphi(2)=f(1)*df(2)*f(3)
   dphi(3)=f(1)*f(2)*df(3)
   d2phi=d2f(1)*f(2)*f(3)+f(1)*d2f(2)*f(3)+f(1)*f(2)*d2f(3)
   end subroutine getorb


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
   real(kind=r8) :: x1(:),x2(:),pf,dpf(:,:),d2pf(:)
   real(kind=r8) :: f1,df1(3),d2f1,f2,df2(3),d2f2
   integer(kind=i4) :: ish,i,j
   real(kind=r8) :: dx(ndim),r,f,df,d2f
   real(kind=r8) :: gradf(2,3),delf(2)
   real(kind=r8) :: xcm(3),rcm,gradk(3),delk
   real(kind=r8) :: vext,dvext,d2vext,kf,dkf,d2kf
   real(kind=r8) :: g1,dg1,d2g1,gam1,gradg1(3),delg1
   real(kind=r8) :: g2,dg2,d2g2,gam2,gradg2(3),delg2
   real(kind=r8) :: g3,dg3,d2g3,gam3,gradg3(3),delg3
   real(kind=r8) :: phi,gradphi(2,3),delphi(2)
   real(kind=r8) :: om
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
   i=0
   om=omega*vorb(norb+3)
   do ish=1,norb
      do j=1,numsh(ish)
         i=i+1
         call getorb(nquant(:,i),om,x1,f1,df1,d2f1)
         call getorb(nquant(:,i),om,x2,f2,df2,d2f2)
         phi=phi+vorb(ish)*f1*f2
         gradphi(1,:)=gradphi(1,:)+vorb(ish)*df1(:)*f2
         gradphi(2,:)=gradphi(2,:)+vorb(ish)*f1*df2(:)
         delphi(1)=delphi(1)+vorb(ish)*d2f1*f2
         delphi(2)=delphi(2)+vorb(ish)*f1*d2f2
      enddo
   enddo
   if (radpair.eq.0) return
   dx=x1-x2
   r=sqrt(sum(dx**2))
   if (radpair.eq.1) then
      call radialfun(kf*r,f,df,d2f)
   else 
      call radialfun2(r,f,df,d2f)
      d2pf=d2pf+d2f+(ndim-1.0_r8)*df/r
   endif
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
!write(99,*) rcm,r,pf
!enddo
!write(99,*) ' '
!write(99,*) ' '
!enddo
!stop
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
   use betafun
   real(kind=r8) :: x(:,:),dpsi(:)
   real(kind=r8) :: phil
   real(kind=r8), dimension(npair+nexup+nexdn,npair+nexup+nexdn,npdet) :: dphi
   real(kind=r8), dimension(npair+nexup+nexdn,npair+nexup+nexdn) :: smati
   real(kind=r8), dimension(2,npair+nexup+nexdn,npair+nexup+nexdn) :: d2ph
   real(kind=r8), dimension(2,ndim,npair+nexup+nexdn,npair+nexup+nexdn) :: dph
   real(kind=r8) :: f1,df1(3),d2f1,f2,df2(3),d2f2,dx(3),r
   integer(kind=i4) :: i,j,jj,ii,k,kk,nd,is,np,ish
   real(kind=r8) :: xcm(3),rcm,kf,vext,g1,g2,g3,gam,pairf
   real(kind=r8) :: f,dum1,dum2,om
real(kind=r8) :: pf1,pf2
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
         k=0
         pairf=0.0_r8
         pf1=0.0_r8
         pf2=0.0_r8
         if (iopt(3)) then
            om=omega*vorb(norb+3)
            do ish=1,norb
               do kk=1,numsh(ish)
                  k=k+1
                  call getorb(nquant(:,k),om,x(:,i),f1,df1,d2f1)
                  call getorb(nquant(:,k),om,x(:,jj),f2,df2,d2f2)
                  dphi(i,j,np)=dphi(i,j,np)+f1*f2
                  pairf=pairf+vorb(ish)*f1*f2
                  call getorb(nquant(:,k),omega*(vorb(norb+3)+0.01),x(:,i),f1,df1,d2f1)
                  call getorb(nquant(:,k),omega*(vorb(norb+3)+0.01),x(:,jj),f2,df2,d2f2)
                  pf1=pf1+vorb(ish)*f1*f2
                  call getorb(nquant(:,k),omega*(vorb(norb+3)-0.01),x(:,i),f1,df1,d2f1)
                  call getorb(nquant(:,k),omega*(vorb(norb+3)-0.01),x(:,jj),f2,df2,d2f2)
                  pf2=pf2+vorb(ish)*f1*f2
               enddo
               dphi(i,j,np)=dphi(i,j,np)*g1*g2
               np=np+1
            enddo
            call radialfun(kf*r,f,dum1,dum2)
            dphi(i,j,np)=-pairf*rcm**2*g1*g2
            np=np+1
            dphi(i,j,np)=-pairf*rcm**2*g1*g2+f*g3*rcm**2*g2
            np=np+1
            dphi(i,j,np)=(pf1-pf2)*g1*g2/0.02_r8
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
   nsw=norb+3
   allocate(vsw(nsw))
   vsw=vorb
   end subroutine getbcstrapampl
end module bcstrap
