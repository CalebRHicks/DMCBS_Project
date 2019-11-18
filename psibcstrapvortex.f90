module bcstrap
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, save :: npart,ntab,ndim
   integer(kind=i4), private, save :: npair
   integer(kind=i4), private, save, allocatable :: nbin(:)
   real(kind=r8), private, save :: range,scale
   integer(kind=i4), private, save :: nexup,nexdn,norb
   integer(kind=i4), private, save, pointer :: numsh(:),nquant(:,:)
   real(kind=r8), private, save, allocatable :: mass(:),vorb(:)
   integer(kind=i4), private, save, allocatable :: spin(:)
   integer(kind=i4), private, save :: npdet
   logical, private, save, allocatable :: iopt(:)
   real(kind=r8), private, save, allocatable :: orb(:,:,:),dorb(:,:,:),d2orb(:,:,:)
   real(kind=r8), private, save :: omega(3),b
   integer(kind=i4), private, save, allocatable :: nkup(:,:),nkdn(:,:)
   integer(kind=i4), private, save :: ntabf
   real(kind=r8), private, save :: rangef,scalef
   real(kind=r8), private, save, allocatable :: ftab(:),dftab(:),d2ftab(:)
   integer(kind=i4), private, save :: radpair
contains
   subroutine setupbcstrap(rho,npartin,nbinin,massin,ndin,ntabin)
   use mympi
   use betafun
   integer(kind=i4) :: npartin,nbinin(:),ndin,ntabin,nmax
   real(kind=r8) :: rho,massin(:),r,dr
   real(kind=r8), pointer :: ene(:)
   integer(kind=i4) :: i,ic,j,k
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
      read (5,*) range
      read (5,*) omega  ! external well frequency
      read (5,*) nexup
      allocate(nkup(ndim,nexup))
      do j=1,nexup
         read (5,*) nkup(:,j)
      enddo
      read (5,*) nexdn
      allocate(nkdn(ndim,nexdn))
      do j=1,nexdn
         read (5,*) nkdn(:,j)
      enddo
      read (5,*) norb  !number of orbitals for s-wave
      allocate(vorb(norb))
      do i=1,norb
         read (5,*) vorb(i) !amplitude
      enddo
      read (5,*) radpair ! 0=no, 1=betafun, 2=fort.80
      read (5,*) b ! parameter of one body confinement
      if (radpair.eq.1) read (5,*) rangef ! range of betafun
      write (6,'(''range of orbitals ='',t40,f10.5)') range
      write (6,'(''omega_x external potential ='',t40,f10.5)') omega(1)
      write (6,'(''omega_y external potential ='',t40,f10.5)') omega(2)
      write (6,'(''omega_z external potential ='',t40,f10.5)') omega(3)
      write (6,'(''unpaired up states ='',t40,i10)') nexup
      write (6,'(''quantum numbers ='',(t35,3i5))') (nkup(:,j),j=1,nexup)
      write (6,'(''unpaired down states ='',t40,i10)') nexdn
      write (6,'(''quantum numbers ='',(t35,3i5))') (nkdn(:,j),j=1,nexdn)
      write (6,'(''s-wave pairing parameters:'')')
      write (6,'(''number of functions: '',t40,i10)') norb
      write (6,'(''amplitudes: '',t40,f10.5)') (vorb(i),i=1,norb)
      if (radpair.eq.0) write (6,'(''no radial function in pairing'')')
      if (radpair.eq.1) write (6,'(''using beta function with range =''t40,f10.5)') rangef
      if (radpair.eq.2) write (6,'(''using fort.80 in pairing'')')
      write (6,'(''b parameter of Jastrow ='',f10.5)') b
   endif
   call bcast(range)
   call bcast(omega)
   call bcast(nexup)
   call bcast(nexdn)
   if (myrank().ne.0) allocate(nkup(ndim,nexup),nkdn(ndim,nexdn))
   call bcast(nkup)
   call bcast(nkdn)
   call bcast(norb)
   if (myrank().ne.0) allocate(vorb(norb))
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
   scale=ntab/range
   dr=1.0_r8/scale
   call sethobasis(omega,norb,numsh,nquant,ene)
   if (myrank().eq.0) then
      k=0
      do i=1,norb
         write (6,'(''shell number and energy '',t30,i10,f10.5)') i,ene(k+1)
         do j=1,numsh(i)
            k=k+1
            write (6,'(''state '',t35,3i5)') nquant(:,k)
         enddo
      enddo
   endif
   nmax=max(maxval(nkup),maxval(nkdn))
   nmax=max(maxval(nquant),nmax)
   allocate(orb(3,0:nmax,-ntab:ntab),dorb(3,0:nmax,-ntab:ntab),d2orb(3,0:nmax,-ntab:ntab))
! orbitals have an index for the coordinate as the frequency may be different
! NOTE: the harmonic oscillator orbitals are correct only if hbar=1 and mass=1
! as it is not yet included in this code!!! fix it!
! I must also put an additional index for different masses!
   do j=0,nmax
      do i=-ntab,ntab
         r=i*dr
         do ic=1,3
            call hoorb(j,omega(ic),r,orb(ic,j,i),dorb(ic,j,i),d2orb(ic,j,i))
         enddo
!        if (myrank().eq.0) write(60+j,'(f10.5,9e15.7)') r, &
!                   orb(1,j,i),dorb(1,j,i),d2orb(1,j,i), &
!                   orb(2,j,i),dorb(2,j,i),d2orb(2,j,i), &
!                   orb(3,j,i),dorb(3,j,i),d2orb(3,j,i)
      enddo
   enddo
   if (myrank().eq.0) write (6,'(''maximum value of orbitals at range ='',t35,f15.10)') &
         maxval(orb(:,:,ntab))
      npair=nbin(1)+nbin(2)-nexup-nexdn
   if (mod(npair,2).ne.0) then
      write (6,'(''inconsistent pairing number'',3i5)') &
         nbin(1)+nbin(2),nexup,nexdn
      stop
   endif
   npair=npair/2
   end subroutine setupbcstrap
 
   subroutine sethobasis(om,nsh,degen,ik,ene)
   real(kind=r8) :: om(:)
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
         ene(i)=(0.5_r8+iktmp(1,i))*om(1)+(0.5_r8+iktmp(2,i))*om(2)+(0.5_r8+iktmp(3,i))*om(3)
         ene(j)=(0.5_r8+iktmp(1,j))*om(1)+(0.5_r8+iktmp(2,j))*om(2)+(0.5_r8+iktmp(3,j))*om(3)
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

   subroutine getorb(n,x,phi,dphi,d2phi)
   integer(kind=i4) :: n(:),index,ic
   real(kind=r8) :: x(:),phi,dphi(:),d2phi
   real(kind=r8) :: f(3),df(3),d2f(3)
   real(kind=r8) :: dr,c1,c2,c3,c4
   do ic=1,3
      if (x(ic).gt.range) then
         f(ic)=0.0_r8
         df(ic)=0.0_r8
         d2f(ic)=0.0_r8
      else
         dr=scale*x(ic)
         index=dr
         index=min(ntab-2,max(-ntab+1,index))
         dr=dr-index
         c1=-dr*(dr-1.0_r8)*(dr-2.0_r8)/6.0_r8
         c2=(dr+1.0_r8)*(dr-1.0_r8)*(dr-2.0_r8)/2.0_r8
         c3=-(dr+1.0_r8)*dr*(dr-2.0_r8)/2.0_r8
         c4=(dr+1.0_r8)*dr*(dr-1.0_r8)/6.0_r8
         f(ic)=c1*orb(ic,n(ic),index-1)+c2*orb(ic,n(ic),index) &
              +c3*orb(ic,n(ic),index+1)+c4*orb(ic,n(ic),index+2)
         df(ic)=c1*dorb(ic,n(ic),index-1)+c2*dorb(ic,n(ic),index) &
               +c3*dorb(ic,n(ic),index+1)+c4*dorb(ic,n(ic),index+2)
         d2f(ic)=c1*d2orb(ic,n(ic),index-1)+c2*d2orb(ic,n(ic),index) &
                +c3*d2orb(ic,n(ic),index+1)+c4*d2orb(ic,n(ic),index+2)
      endif
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
   real(kind=r8) :: x1(:),x2(:)
   complex(kind=r8) :: pf,dpf(:,:),d2pf(:)
   real(kind=r8) :: f1,df1(3),d2f1,f2,df2(3),d2f2
   integer(kind=i4) :: ish,i,j
   real(kind=r8) :: dx(ndim),r,f,df,d2f
   complex(kind=r8) :: ex,dex(2,3),d2ex(2)
   pf=0.0_r8
   dpf=0.0_r8
   d2pf=0.0_r8
   i=0
   do ish=1,norb
      do j=1,numsh(ish)
         i=i+1
         call getorb(nquant(:,i),x1,f1,df1,d2f1)
         call getorb(nquant(:,i),x2,f2,df2,d2f2)
         pf=pf+vorb(ish)*f1*f2
         dpf(1,:)=dpf(1,:)+vorb(ish)*df1(:)*f2
         dpf(2,:)=dpf(2,:)+vorb(ish)*f1*df2(:)
         d2pf(1)=d2pf(1)+vorb(ish)*d2f1*f2
         d2pf(2)=d2pf(2)+vorb(ish)*f1*d2f2
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
! add pair-excitation
   call excitation(x1,x2,ex,dex,d2ex)
   d2pf(1)=d2pf(1)*ex+2.0_r8*sum(dpf(1,:)*dex(1,:))+pf*d2ex(1)
   d2pf(2)=d2pf(2)*ex+2.0_r8*sum(dpf(2,:)*dex(2,:))+pf*d2ex(2)
   dpf(1,:)=dpf(1,:)*ex+pf*dex(1,:)
   dpf(2,:)=dpf(2,:)*ex+pf*dex(2,:)
   pf=pf*ex
   end subroutine pairfn

   subroutine excitation(x1,x2,ex,dex,d2ex)
   real(kind=r8) :: x1(:),x2(:),rat,mod,xx,yy,xx2,yy2
   complex(kind=r8) :: ex,dex(:,:),d2ex(:)
   complex(kind=r8), parameter :: ci=(0.0_r8,1.0_r8)
! here I assume the two particles have the same mass
   xx=x1(1)+x2(1)
   yy=x1(2)+x2(2)
!  xx2=xx*xx
!  yy2=yy*yy
!  rat=yy/xx
!  mod=sqrt(1.0_r8+rat**2)
!  ex=cmplx(1.0_r8,rat,r8)/mod
!  dex(1,1)=-ci*yy/(xx2*mod)+(1.0_r8+ci*yy/xx)*yy2/(mod*mod*mod*xx2*xx)
!  dex(1,2)=ci/(xx*mod)-(1.0_r8+ci*yy/xx)*yy/(mod*mod*mod*xx2)
!  dex(2,1)=dex(1,1)
!  dex(2,2)=dex(1,2)
!  dex(:,3)=0.0_r8
!  d2ex(1)=-(xx+ci*yy)/((xx2+yy2)*sqrt((xx2+yy2)/xx2)*xx)
!  d2ex(2)=d2ex(1)
ex=xx+ci*yy
dex(1,1)=1.0_r8
dex(1,2)=ci
dex(2,1)=1.0_r8
dex(2,2)=ci
d2ex=0.0_r8
   end subroutine excitation

   subroutine psibcstrap(x,phil,is,dphi,d2phi)
   use matrixmod
   real(kind=r8), dimension(ndim,npart) :: x,dphi
   real(kind=r8) :: d2phi
   complex(kind=r8) :: smati(npair+nexup+nexdn,npair+nexup+nexdn)
   complex(kind=r8) :: dph(2,ndim,npair+nexup+nexdn,npair+nexup+nexdn)
   complex(kind=r8) :: d2ph(2,npair+nexup+nexdn,npair+nexup+nexdn)
   complex(kind=r8) :: det,ddet(3,npart),d2det,phi
   real(kind=r8) :: corb,dcorb(3),d2corb
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
         call getorb(nkup(:,j),x(:,i),corb,dcorb,d2corb)
         smati(i,jj)=corb
         dph(1,:,i,jj)=dcorb
         d2ph(1,i,jj)=d2corb
      enddo
   enddo
   do i=1,npair+nexdn
      ii=i+npair+nexup
      do j=1,nexdn
         jj=j+npair+nexup
         call getorb(nkdn(:,j),x(:,ii),corb,dcorb,d2corb)
         smati(jj,i)=corb
         dph(2,:,jj,i)=dcorb
         d2ph(2,jj,i)=d2corb
      enddo
   enddo
   call cmatinv(smati,det,nd)
   ddet=0.0_r8
   d2det=0.0_r8
   do i=1,npair+nexup
      do ic=1,ndim
         ddet(ic,i)=sum(smati(:,i)*dph(1,ic,i,:))/sqrt(mass(i))
      enddo
      d2det=d2det+sum(smati(:,i)*d2ph(1,i,:))/mass(i)
   enddo
   do i=1,npair+nexdn
      j=npair+nexup+i
      do ic=1,ndim
         ddet(ic,j)=sum(smati(i,:)*dph(2,ic,:,i))/sqrt(mass(j))
      enddo
      d2det=d2det+sum(smati(i,:)*d2ph(2,:,i))/mass(j)
   enddo
   phi=exp(det)
   dphi=real(ddet*phi)/real(phi)
   d2phi=real(d2det*phi)/real(phi)
   phil=log(abs(real(phi)))
   is=sign(1.0_r8,real(phi))
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
      vext=vext+0.5_r8*mass(i)*sum((omega(:)*x(:,i))**2)
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
   real(kind=r8) :: x(:,:),dpsi(:)
   real(kind=r8) :: phil
   real(kind=r8), dimension(npair+nexup+nexdn,npair+nexup+nexdn,npdet) :: dphi
   real(kind=r8), dimension(npair+nexup+nexdn,npair+nexup+nexdn) :: smati
   real(kind=r8), dimension(2,npair+nexup+nexdn,npair+nexup+nexdn) :: d2ph
   real(kind=r8), dimension(2,ndim,npair+nexup+nexdn,npair+nexup+nexdn) :: dph
   real(kind=r8) :: f1,df1(3),d2f1,f2,df2(3),d2f2
   integer(kind=i4) :: i,j,jj,ii,k,kk,nd,is,np,ish
   dpsi=0.0_r8
   nd=npair+nexup+nexdn
   smati=0.0_r8
   dph=0.0_r8
   d2ph=0.0_r8
   do i=1,npair+nexup
      do j=1,npair+nexdn
         jj=j+npair+nexup
!        call pairfn(x(:,i),x(:,jj),smati(i,j),dph(:,:,i,j),d2ph(:,i,j))
      enddo
   enddo
   do i=1,npair+nexup
      do j=1,nexup
         jj=j+npair+nexdn
!        call getorb(nkup(:,j),x(:,i),smati(i,jj),dph(1,:,i,jj),d2ph(1,i,jj))
      enddo
   enddo
   do i=1,npair+nexdn
      ii=i+npair+nexup
      do j=1,nexdn
         jj=j+npair+nexup
!        call getorb(nkdn(:,j),x(:,ii),smati(jj,i),dph(2,:,jj,i),d2ph(2,jj,i))
      enddo
   enddo
   call matinv(smati,phil,is,nd)
   dphi=0.0_r8
   do i=1,npair+nexup
      do j=1,npair+nexdn
         jj=j+npair+nexup
         np=1
         k=0
         if (iopt(3)) then
            do ish=1,norb
               do kk=1,numsh(ish)
                  k=k+1
                  call getorb(nquant(:,k),x(:,i),f1,df1,d2f1)
                  call getorb(nquant(:,k),x(:,jj),f2,df2,d2f2)
                  dphi(i,j,np)=dphi(i,j,np)+f1*f2
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
