!buggy
module pdistribution
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4) ,private, save :: ntab,ntot,nblk
   real(kind=r8), private, save :: el,h,dk,wtblk,wttot
   real(kind=r8), private, allocatable, save :: enkcub(:,:,:),rho(:,:,:)
   real(kind=r8), private, allocatable, save :: enknow(:),enk(:),cs(:)
   real(kind=r8), private, allocatable, save :: enkbad(:),enkbad2(:)

contains
   subroutine setupenk(ndim,ntabin,elin)
   use trans3d
   use mympi
   integer(kind=i4) :: ntabin,ndim
   real(kind=r8) :: elin
   integer(kind=i4) :: i
   real(kind=r8) :: dk,pi
   if (ndim.ne.3) call abort
   ntab=ntabin
   ntot=((ntab+3)*(ntab+2)*(ntab+1))/6
   el=elin
   h=el/ntab
   allocate(enkcub(0:ntab,0:ntab,0:ntab),rho(0:ntab,0:ntab,0:ntab),cs(0:ntab))
   allocate(enknow(ntot),enk(ntot),enkbad(ntot),enkbad2(ntot))
   call tinit(ntab,ntab,ntab,(/el,el,el/))
   pi=4.0_r8*atan(1.0_r8)
   dk=2.0_r8*pi/el
   end subroutine setupenk

   subroutine addenk(w)
   use stack
   use random
   use trans3d
   type (walker) w
   real(kind=r8) a(3),dx(3,8),xnew(3),psirat
   integer(kind=i4) n,i,j,ix,iy,iz,imax,imin,imid,k
   n=size(w%x)/3
   rho=0.0_r8
   a=0.5_r8*h*randn(3)
   do ix=0,ntab
      dx(1,1:4)=ix*h+a(1)
      dx(1,5:8)=-ix*h+a(1)
      do iy=0,ntab
         dx(2,1:2)=iy*h+a(2)
         dx(2,3:4)=-iy*h+a(2)
         dx(2,5:6)=iy*h+a(2)
         dx(2,7:8)=-iy*h+a(2)
         do iz=0,ntab
            dx(3,1:7:2)=iz*h+a(3)
            dx(3,2:8:2)=-iz*h+a(3)
            do i=1,n
               do j=1,8
                  xnew=w%x(:,i)+dx(:,j)
                  call move1(i,xnew,w,psirat)
                  rho(ix,iy,iz)=rho(ix,iy,iz)+psirat
               enddo
            enddo
         enddo
      enddo
   enddo
   do ix=0,ntab
      cs(ix)=cos(ix*dk)
   enddo
   enkcub=rtok(rho,(/.true.,.true.,.true./))
   do iz=0,ntab
      do iy=0,ntab
         do ix=0,ntab
            imax=max(ix,iy,iz)
            imin=min(ix,iy,iz)
            imid=ix+iy+iz-imax-imin
            k=(imax*(imax+1)*(imax+2))/6+(imid*(imid+1))/2+imin+1
            enknow(k)=enknow(k)+cs(ix)*cs(iy)*cs(iz)*enkcub(ix,iy,iz)
         enddo
      enddo
   enddo
   end subroutine addenk

!   subroutine move1(i,xnew,w,psirat)
!   use stack
!   use wavefunction
!   type (walker) w
!   real(kind=r8) :: xnew(3),psirat
!   integer(kind=i4) :: i
!   w1=w
!   w%x(:,i)=xnew
!   call hpsi(w)
!   psirat=w%is*w1%is*exp(w%psil-w1%psil)
!   w=w1
!   end subroutine move1

   subroutine zeroenk
   enk=0.0_r8
   enknow=0.0_r8
   enkbad=0.0_r8
   enkbad2=0.0_r8
   wttot=0.0_r8
   wtblk=0.0_r8
   nblk=0
   end subroutine zeroenk

   subroutine updateenk
   use mympi
   real(kind=r8) :: anorm,val
   integer(kind=i4) :: i
   real(kind=r8) :: enktnow(ntot),wttblk
   call addall(enknow,enktnow)
   call addall(wtblk,wttblk)
   if (myrank().eq.0) then
      wtblk=wttblk
      enknow=enktnow
      anorm=1.0_r8/wtblk
      do i=0,ntot
         val=anorm*enknow(i)
         enkbad(i)=enkbad(i)+val
         enkbad2(i)=enkbad2(i)+val*val
      enddo
      enk=enk+enknow
      wttot=wttot+wtblk
      nblk=nblk+1
   endif
   wtblk=0.0_r8
   enknow=0.0_r8
   end subroutine updateenk

   subroutine writeenk(filename)
   use mympi
   character(len=*) :: filename
   integer(kind=i4) :: k,imax,imid,imin,ifacx,ifacy,ifacz
   real(kind=r8) :: av,av2,err
   if (myrank().ne.0) return
   open(unit=78,file=filename,form='formatted')
   rewind 78
   k=0
   do imax=0,ntab
      ifacx=1+min(1,imax)
      do imid=0,imax
         ifacy=1+min(1,imid)
         do imin=0,imid
            ifacz=1+min(1,imin)
            k=k+1
            enk(k)=enk(k)/(ifacx*ifacy*ifacz*max(1,nblk))
            av=enkbad(k)/(ifacx*ifacy*ifacz*max(1,nblk))
            av2=enkbad2(k)/(ifacx*ifacy*ifacz*max(1,nblk))
            err=sqrt(abs(av2-av*av)/max(1,nblk-1))
            write (78,'(3i5,2e15.7)') imax,imid,imin,enk(k),err
         enddo
      enddo
   enddo
   close(78)
   end subroutine writeenk

end module pdistribution
