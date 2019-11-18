module wavefunction
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, save :: npart,ndim
   logical, private, save :: wfbose,wfslat,wfbcs,wfpfaf,wfatom,wfbcstrap,wfloff
   integer(kind=i4), private, save :: nexup,nexdn
   integer(kind=i4), private, save :: nsw,npw,npwup,npwdn
   integer(kind=i4), allocatable :: ikup(:,:),ikdn(:,:)
   real(kind=r8), private, save, allocatable :: vsw(:),vpw(:),vpwup(:),vpwdn(:)
   real(kind=r8), private, save :: a2
contains
   subroutine setpsi(ndin,nbin,el,mass,ntab)
   use mympi
   use bose
   use slater
   use slaterbcs
   use slaterbcsloff
   use pfaffian
   use atom
   use bcstrap
   use betafun
   integer(kind=i4) :: ndin,nbin(:),ntab
   real(kind=r8) :: el,mass(:),rho
   character(len=8) :: nodes
   rho=0.0_r8
   ndim=ndin
   npart=sum(nbin)
   wfbose=.false.
   wfslat=.false.
   wfbcs=.false.
   wfpfaf=.false.
   wfatom=.false.
   wfbcstrap=.false.
   wfloff=.false.
   if (myrank().eq.0) read (5,'(a8)') nodes
   if (myrank().eq.0) read (5,*) a2
   call bcast(nodes)
   call bcast(a2)
   if (myrank().eq.0) write (6,'(''!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!'')')
   if (nodes.eq.'bose') then
      if (myrank().eq.0) write(6,'(''using symmetric w.f.'')')
      call setupbose(rho)
      wfbose=.true.
   else if (nodes.eq.'slater') then
      if (myrank().eq.0) write(6,'(''using Slater wavefunction'')')
      call setupslat(rho,npart,nbin,mass,ndim,.false.,ntab)
      wfslat=.true.
   else if (nodes.eq.'slaterbf') then
      if (myrank().eq.0) write(6,'(''using Slater wavefunction with backflow'')')
      call setupslat(rho,npart,nbin,mass,ndim,.true.,ntab)
      wfslat=.true.
   else if (nodes.eq.'bcs') then
      if (myrank().eq.0) write(6,'(''using SlaterBCS wavefunction'')')
      call setupbcs(rho,npart,nbin,mass,ndim,ntab)
      wfbcs=.true.
   else if (nodes.eq.'pfaffian') then
      if (myrank().eq.0) write(6,'(''using Pfaffian wavefunction'')')
      call setuppfaf(rho,npart,nbin,mass,ndim)
      el=(npart/rho)**(1.0_r8/ndim)
      call setbeta(ndim,ntab,el)
      wfpfaf=.true.
   else if (nodes.eq.'atom') then
      if (myrank().eq.0) write(6,'(''using atomic wavefunction'')')
      call setuppsiatom(rho,npart,nbin,mass,ndim)
      wfatom=.true.
   else if (nodes.eq.'bcstrap') then
      if (myrank().eq.0) write(6,'(''using SlaterBCStrap wavefunction'')')
      call setupbcstrap(rho,npart,nbin,mass,ndim,ntab)
      wfbcstrap=.true.
   else if (nodes.eq.'loff') then
      if (myrank().eq.0) write(6,'(''using SlaterBCSloff wavefunction'')')
      call setupbcsloff(rho,npart,nbin,mass,ndim)
      el=(npart/rho)**(1.0_r8/ndim)
      call setbeta(ndim,ntab,el)
      wfloff=.true.
   else
      if (myrank().eq.0) then
         write(6,'(''Available wavefunctions:'')')
         write(6,'(''bose - Bose gas'')')
         write(6,'(''slater - Fermi gas'')')
         write(6,'(''slaterbf - Fermi gas with backflow'')')
         write(6,'(''bcs - BCS Fermi gas'')')
         write(6,'(''bcsloff - BCSloff Fermi gas'')')
         write(6,'(''pfaffian - Pfaffian Fermi gas'')')
         write(6,'(''atom - cluster of Fermions'')')
         write(6,'(''bcstrap - trapped BCS Fermions'')')
      endif
   endif
   if (myrank().eq.0) write (6,'(''a constant in the importance function ='',(t40,e12.5))') a2
   if (rho.ne.0.0_r8) then
      el=(npart/rho)**(1.0_r8/ndim)
   else
      el=0.0_r8
   endif
   end subroutine setpsi

   subroutine getpsi(x,phil,is,dphi,d2phi)
   use slater
   use slaterbcs
   use slaterbcsloff
   use pfaffian
   use atom
   use bcstrap
   real(kind=r8) :: x(:,:),phil,dphi(:,:),d2phi
   integer(kind=i4) :: is
   if (wfslat) then
      call slat(x,phil,is,dphi,d2phi)
      return
   else if (wfbcs) then
      call bcs(x,phil,is,dphi,d2phi)
      return
   else if (wfpfaf) then
      call pfaff(x,phil,is,dphi,d2phi)
      return
   else if (wfbose) then
      phil=0.0_r8
      is=1
      dphi=0.0_r8
      d2phi=0.0_r8
   else if (wfatom) then
      call psiatom(x,phil,is,dphi,d2phi)
      return
   else if (wfbcstrap) then
      call psibcstrap(x,phil,is,dphi,d2phi)
      return
   else if (wfloff) then
      call bcsloff(x,phil,is,dphi,d2phi)
      return
   endif
   end subroutine getpsi

   subroutine getonebody(x,u1,du1,d2u1,v1)
   use bcstrap
   use slaterbcs
   real(kind=r8) :: x(:,:),u1,du1(:,:),d2u1,v1
   if (wfbcstrap) then
      call onebody(x,u1,du1,d2u1,v1)
   else if (wfbcs) then
      call bcsonebody(x,u1,du1,d2u1,v1)
   else
      u1=0.0_r8
      du1=0.0_r8
      d2u1=0.0_r8
      v1=0.0_r8
   endif
   end subroutine getonebody

   subroutine getphilong(x,u,du,d2u)
   use jaslong
   real(kind=r8) :: x(:,:),u,du(:,:),d2u
   if (wfbcs) then
      call getjaslong(x,u,du,d2u)
   else 
      u=0.0_r8
      du=0.0_r8
      d2u=0.0_r8
   endif
   end subroutine getphilong
      
   subroutine hpsi(w)
   use stack
   use potential
   type (walker) w
   real(kind=r8), dimension(ndim,npart) :: dphi,du,du1,dulong
   real(kind=r8) :: philog,d2phi,u,d2u,v1,u1,d2u1
   real(kind=r8) :: ulong,d2ulong
   real(kind=r8), dimension(ndim*npart) :: d1,d2,d3
   integer(kind=i4) :: isj
   real(kind=r8) :: psi,psig
   call getphi(w%x,u,isj,du,d2u,w%v)
   call getonebody(w%x,u1,du1,d2u1,v1)
   call getphilong(w%x,ulong,dulong,d2ulong)
   call getpsi(w%x,philog,w%is,dphi,d2phi)
   w%vext=v1
   w%v=w%v+v1
   w%is=isj*w%is
! there is also a new variable isg (unused so far) in the stack
   u=u+u1+ulong
   d1=reshape(du,(/ndim*npart/))
   d2=reshape(du1,(/ndim*npart/))
   d3=reshape(dulong,(/ndim*npart/))
   d2u=d2u+d2u1+d2ulong+2.0_r8*(dot_product(d1,d2)+dot_product(d1,d3)+dot_product(d2,d3))
   du=du+du1+dulong
   w%psil=philog+u
   w%dpsi=dphi+du
   d1=reshape(du,(/ndim*npart/))
   d2=reshape(dphi,(/ndim*npart/))
   w%d2psi=d2phi+d2u+2.0_r8*dot_product(d1,d2)
   if (a2.eq.0.0_r8) then
      w%dpsig=w%dpsi
      w%d2psig=w%d2psi
      w%psigl=w%psil
   else
      psi=exp(w%psil)
      psig=sqrt(psi**2+a2)
      w%dpsig=psi**2*w%dpsi/psig**2
      w%d2psig=(sum(w%dpsi**2)+w%d2psi)*psi**2/psig**2-psi**4*sum(w%dpsi**2)/psig**4
      w%psigl=log(psig)
   endif
   end subroutine hpsi

   subroutine getderpsi(x,dpsi)
   use slater
   use slaterbcs
   use pfaffian
   use bcstrap
   use slaterbcsloff
   real(kind=r8) :: x(:,:),dpsi(:)
   if (wfslat) then
      call getderslat(x,dpsi)
      return
   else if (wfbcs) then
      call getderbcs(x,dpsi)
      return
   else if (wfpfaf) then
      call getderpfaf(x,dpsi)
      return
   else if (wfbcstrap) then
      call getderbcstrap(x,dpsi)
      return
   else if (wfbose) then
      dpsi=0.0_r8
      return
   else if (wfloff) then
      call getderbcsloff(x,dpsi)
      return
   endif
   end subroutine getderpsi

   subroutine setdetparam(par)
   use slater
   use slaterbcs
   use pfaffian
   use bcstrap
   use slaterbcsloff
   real(kind=r8) :: par(:)
!real(kind=r8) :: v0(1),v1(1),v2(1),v3(1)
   if (wfslat) then
      call setslatpar(par)
      return
   else if (wfbcs) then
      call setbcspar(par)
      return
   else if (wfpfaf) then
!     call setvorbpfaf(v0,v1,v2,v3)
      call setvorbpfaf(par)
      return
   else if (wfbcstrap) then
      call setvorbbcstrap(par)
      return
   else if (wfbose) then
      return
   else if (wfloff) then
!     call setvorbbcsloff(v0)
      return
   endif
   end subroutine setdetparam

   subroutine getdetparam(npar,par)
   use slater
   use slaterbcs
   use pfaffian
   use bcstrap
   use slaterbcsloff
   integer(kind=i4) :: npar
   real(kind=r8), pointer :: par(:)
!integer(kind=i4) :: nsw,npw,npwup,npwdn
!real(kind=r8), pointer :: vsw(:),vpw(:),vpwup(:),vpwdn(:)
   if (wfslat) then
      call getslatpar(npar,par)
      return
   else if (wfbcs) then
      call getbcspar(npar,par)
      return
   else if (wfpfaf) then
!     call getpfafampl(nsw,npw,npwup,npwdn,vsw,vpw,vpwup,vpwdn)
      call getpfafampl(npar,par)
      return
   else if (wfbcstrap) then
      call getbcstrapampl(npar,par)
      return
   else if (wfloff) then
!     call getbcsamplloff(nsw,vsw)
      return
   endif
   end subroutine getdetparam

   subroutine chkder(w,dx,error,massin,nbin)
   use stack
   type (walker) w
   real(kind=r8) :: dx,error,massin(:)
   real(kind=r8) :: dnum,d2num,d2n
   real(kind=r8) :: mass(npart)
   integer(kind=i4) :: i,ic,nbin(:)
   mass(1:nbin(1))=massin(1)
   mass(nbin(1)+1:npart)=massin(2)
   d2num=0.0_r8
   do i=1,npart
      do ic=1,ndim
         w%x(ic,i)=w%x(ic,i)+dx
         call hpsi(w)
         dnum=w%psil
         d2n=w%dpsi(ic,i)*exp(w%psil)
         w%x(ic,i)=w%x(ic,i)-2.0_r8*dx
         call hpsi(w)
         dnum=(dnum-w%psil)/(2.0_r8*dx)
         d2n=(d2n-w%dpsi(ic,i)*exp(w%psil))/(2.0_r8*dx)
         w%x(ic,i)=w%x(ic,i)+dx
         call hpsi(w)  
         dnum=dnum/sqrt(mass(i))
         d2num=d2num+d2n/sqrt(mass(i))/exp(w%psil)
         if (abs(w%dpsi(ic,i)-dnum).gt.error) then
            write (6,'(1x,''dpsi ic i analytic numerical '',2i5,1p,2e14.6)') &
               ic,i,w%dpsi(ic,i),dnum
         endif
      enddo
   enddo
   if (abs(d2num-w%d2psi).gt.error) then
      write (6,'(1x,''d2psi analytic numerical '',1p,2e14.6)') w%d2psi,d2num
   endif
   if (a2.eq.0.0_r8) return
   d2num=0.0_r8
   do i=1,npart
      do ic=1,ndim
         w%x(ic,i)=w%x(ic,i)+dx
         call hpsi(w)
         dnum=w%psigl
         d2n=w%dpsig(ic,i)*exp(w%psigl)
         w%x(ic,i)=w%x(ic,i)-2.0_r8*dx
         call hpsi(w)
         dnum=(dnum-w%psigl)/(2.0_r8*dx)
         d2n=(d2n-w%dpsig(ic,i)*exp(w%psigl))/(2.0_r8*dx)
         w%x(ic,i)=w%x(ic,i)+dx
         call hpsi(w)  
         dnum=dnum/sqrt(mass(i))
         d2num=d2num+d2n/sqrt(mass(i))/exp(w%psigl)
         if (abs(w%dpsig(ic,i)-dnum).gt.error) then
            write (6,'(1x,''dpsig ic i analytic numerical '',2i5,1p,2e14.6)') &
               ic,i,w%dpsig(ic,i),dnum
         endif
      enddo
   enddo
   if (abs(d2num-w%d2psig).gt.error) then
      write (6,'(1x,''d2psig analytic numerical '',1p,2e14.6)') w%d2psig,d2num
   endif
!write (6,*) ' calling chkone'
!call chkone(w,1.0_r8)
   end subroutine chkder

   subroutine chkone(w,dx)
   use stack
   use slaterbcs
   use potential
   type (walker) w
   real(kind=r8) :: dx,psio,psin,rat1,rat2,dr(3)
   integer(kind=i4) :: i,ic,is,isn
   real(kind=r8) :: dummy(3,npart),dummy1
   if (.not.wfbcs) then
      write (6,'(''chkone only works for slaterbcs !'')')
      return
   endif
   write (6,*) 'check bcsone'
   do i=1,npart
      do ic=1,ndim
         w%x(ic,i)=w%x(ic,i)-dx
         call bcs(w%x,psin,isn,dummy,dummy1)
         w%x(ic,i)=w%x(ic,i)+dx
         call bcs(w%x,psio,is,dummy,dummy1)
         rat1=psin-psio
         dr=0.0_r8
         dr(ic)=dx
         call bcsone(w%x,dr,i,rat2)
         write (6,'(1x,'' log(psin/psio): ic i bcs bcs1 '',2i5,1p,2e14.6)') ic,i,exp(rat1)*is*isn,rat2
      enddo
   enddo
   write (6,*) 'check phione'
   do i=1,npart
      do ic=1,ndim
         w%x(ic,i)=w%x(ic,i)-dx
         call getphi(w%x,psin,is,dummy,dummy1,dummy1)
         w%x(ic,i)=w%x(ic,i)+dx
         call getphi(w%x,psio,is,dummy,dummy1,dummy1)
         rat1=psin-psio
         dr=0.0_r8
         dr(ic)=dx
         call getphione(w%x,dr,i,rat2)
         write (6,'(1x,'' log(phin/phio): ic i bcs bcs1 '',2i5,1p,2e14.6)') ic,i,rat1,log(rat2)
      enddo
   enddo
   end subroutine chkone
   
   subroutine move1(i,xnew,w,psirat)
   use stack
   use slaterbcs
   use potential
   type (walker) w
   integer(kind=i4) :: i
   real(kind=r8) :: xnew(3),psirat,drr(3),rat1,rat2
   drr=w%x(:,i)-xnew(:)
   call getphione(w%x,drr,i,rat1)
   call bcsone(w%x,drr,i,rat2)
   psirat=rat1*rat2
   end subroutine move1

end module wavefunction
