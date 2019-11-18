module potential
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   real(kind=r8), private, parameter :: a=544850.4_r8
   real(kind=r8), private, parameter :: alp=13.353384_r8
   real(kind=r8), private, parameter :: d=1.241314_r8
   real(kind=r8), private, parameter :: c6=1.3732412_r8
   real(kind=r8), private, parameter :: c8=0.4253785_r8
   real(kind=r8), private, parameter :: c10=0.1781_r8
   real(kind=r8), private, parameter :: sigma=1.0_r8
   real(kind=r8), private, parameter :: rm=2.9673_r8
   real(kind=r8), private, parameter :: eps=10.8_r8
   integer(kind=i4), private, save :: npart,ntab,ndim
   real(kind=r8), private, save :: el,range,scale,el2
   integer(kind=i4), private, save :: mmcmillan
   real(kind=r8), private, save :: bmcmillan,rfit,alpha,u0
   real(kind=r8), private, save :: vtail
   real(kind=r8), private, allocatable, save :: vtab(:)
   real(kind=r8), private, save, allocatable :: utab(:),dutab(:),d2utab(:)
contains
   subroutine setphi(ndin,ntabin,nbinin,elin,massin,hbin)
   use mympi
   real(kind=r8) :: elin,massin(2),hbin,dummy(2),dummy1
   real(kind=r8) :: r,dr,rho
   integer(kind=i4) :: ntabin,nbinin(:),i,ndin
   if (ndin.ne.3) then
      if (myrank().eq.0) write(6,'(''This code is for 3D Lennard-Jones, stopping'')')
      stop
   endif
   dummy=massin
   dummy1=hbin
   ndim=3
   ntab=ntabin
   el=elin
   npart=sum(nbinin)
   if (myrank().eq.0) then
      read(5,*) bmcmillan
      read(5,*) mmcmillan
      write(6,'(''Aziz HFDHE2 potential'')')
      write(6,'(''McMillan Jastrow parameters:'')')
      write(6,'(''b ='',t40,f10.5)') bmcmillan
      write(6,'(''m ='',t45,i5)') mmcmillan
   endif
   call bcast(bmcmillan)
   call bcast(mmcmillan)
   el2=0.5_r8*el
   range=0.5_r8*el
   dr=range/ntab
   scale=1.0_r8/dr
   allocate(vtab(0:ntab))
   allocate(utab(0:ntab),dutab(0:ntab),d2utab(0:ntab))
   call initjas
   do i=0,ntab
      r=i*dr
      vtab(i)=vpot(r)
      call jastro(r,utab(i),dutab(i),d2utab(i))
   enddo
   if (myrank().eq.0) then
      do i=1,ntab
!        if (myrank().eq.0) write(78,'(f10.5,3e15.7)') &
!           i*dr*0.5_r8,utab(i),dutab(i),d2utab(i)
      enddo
   endif
   rho=npart/el**3
   call vpottail(vtail,npart,rho)
   if (myrank().eq.0) write(6,'(''tail of potential ='',t35,f15.10)') vtail
   end subroutine setphi

   subroutine initjas
   real(kind=r8) :: em1,em3,ufit
   em1=mmcmillan+1
   em3=mmcmillan+3
   rfit=em1*el2/em3
   ufit=0.5d0*(bmcmillan/rfit)**mmcmillan
   alpha=-mmcmillan*ufit/(3.0d0*rfit*(rfit-el2)**2)
   u0=alpha*(rfit-el2)**3-ufit
   end subroutine initjas


   subroutine jastro(r,u,up,upp)
   real(kind=r8) :: r,u,up,upp
   if (r.gt.el2) then
      u=0.0d0
      up=0.0d0
      upp=0.0d0
   else if (r.lt.rfit) then
      u=-0.5_r8*(bmcmillan/r)**mmcmillan
      up=-mmcmillan*u/r
      upp=mmcmillan*(mmcmillan+1)*u/r**2
      u=u-u0
   else
      u=-alpha*(r-el2)**3
      up=-alpha*3.0d0*(r-el2)**2
      upp=-alpha*6.0d0*(r-el2)
   endif
   end subroutine jastro

   subroutine getphi(x,u,is,du,d2u,v)
   real(kind=r8), dimension(ndim,npart) :: x,du
   real(kind=r8) :: v,u,d2u,uj,duj,d2uj
   integer(kind=i4) :: is,i,j,index
   real(kind=r8), dimension(ndim) :: dx
   real(kind=r8) :: r,dr,c1,c2,c3,c4
   u=0.0_r8
   du=0.0_r8
   d2u=0.0_r8
   v=0.0_r8
   do i=1,npart-1
      do j=i+1,npart
         dx(:)=x(:,i)-x(:,j)
         dx=dx-el*nint(dx/el)
         r=sqrt(sum(dx**2))
         if (r.lt.range) then
            dr=scale*r
            index=dr
            index=max(1,min(index,ntab-2))
            dr=dr-index
            c1=-dr*(dr-1.0_r8)*(dr-2.0_r8)/6.0_r8
            c2=(dr+1.0_r8)*(dr-1.0_r8)*(dr-2.0_r8)/2.0_r8
            c3=-(dr+1.0_r8)*dr*(dr-2.0_r8)/2.0_r8
            c4=(dr+1.0_r8)*dr*(dr-1.0_r8)/6.0_r8
            v=v+c1*vtab(index-1)+c2*vtab(index) &
              +c3*vtab(index+1)+c4*vtab(index+2)
!           uj=c1*utab(index-1)+c2*utab(index) &
!             +c3*utab(index+1)+c4*utab(index+2)
!           duj=c1*dutab(index-1)+c2*dutab(index) &
!              +c3*dutab(index+1)+c4*dutab(index+2)
!           d2uj=c1*d2utab(index-1)+c2*d2utab(index) &
!               +c3*d2utab(index+1)+c4*d2utab(index+2)
            call jastro(r,uj,duj,d2uj)
            u=u+uj
            du(:,i)=du(:,i)+dx(:)*duj/r
            du(:,j)=du(:,j)-dx(:)*duj/r
            d2u=d2u+d2uj+(ndim-1.0_r8)*duj/r
         endif
      enddo
   enddo
   is=1
   d2u=2.0_r8*d2u+sum(du*du)
   v=v+vtail
   end subroutine getphi

   function vpot(r)
   real(kind=r8) :: vpot,r
   real(kind=r8) :: y,v1,x2,f,dy,v2
   real(kind=r8), parameter :: small=1e-10_r8
   if (r.eq.0.0_r8) r=small
   y=rm/r
   v1=a*eps*exp(-alp/y)
   x2=y*y
   f=eps
   dy=d*y
   v2=(c6+(c8+c10*x2)*x2)*x2**3
   if(dy.gt.1.0_r8) f=f*exp(-(dy-1.0_r8)**2)
   vpot=v1-f*v2
   end function vpot

   subroutine vpottail(tail,npart,rho)
   real(kind=r8) :: tail,rho
   integer(kind=i4) :: npart
   real(kind=r8) :: pi,x2,vint,al2,alporm,el2,yfit
   el2=0.5_r8*(npart/rho)**(1.0_r8/3.0_r8)
   pi=4.0_r8*atan(1.0_r8)
   yfit=rm/(sigma*el2)
   x2=yfit**2
   vint=(c6/3.0_r8+(c8/5.0_r8+c10*x2/7.0_r8)*x2)*x2*yfit
   tail=-4.0_r8*eps*pi*vint*(rm/sigma)**3
   al2=el2*sigma
   alporm=alp/rm
   tail=tail+4.0_r8*pi*eps*a*(((al2+2.0_r8/alporm)*al2+2.0_r8/alporm**2)/alporm)*exp(-alporm*al2)
   tail=tail*0.5_r8*rho*npart
   end subroutine vpottail

   subroutine getphione(x,drr,i,phiratio)
   real(kind=r8) :: x(:,:),drr(:),phiratio
   integer(kind=i4) :: i
! fake, subroutine to be implemented
   end subroutine getphione

   subroutine getjasparam(npjas,pjas)
   integer(kind=i4) :: npjas
   real(kind=r8), pointer :: pjas(:)
   npjas=1  ! number of Jastrow parameters
   allocate(pjas(npjas))
   if (npjas.ne.0) pjas=0.0_r8
   pjas(1)=bmcmillan
   end subroutine getjasparam

   subroutine setjasparam(pjas)
   real(kind=r8) :: pjas(:)
   bmcmillan=pjas(1)
   call initjas
   end subroutine setjasparam

   subroutine getjasder(x,dpsi)
   real(kind=r8) :: x(:,:),dpsi(:)
   real(kind=r8) :: psi1,psi2
   real(kind=r8) :: u,dum2(ndim,npart),dum3,dum4
   real(kind=r8), parameter :: dp=0.0001_r8
   integer(kind=i4) :: dum1
   dpsi=0.0_r8
   bmcmillan=bmcmillan+dp
   call initjas
   call getphi(x,u,dum1,dum2,dum3,dum4)
   psi1=u
   bmcmillan=bmcmillan-2.0_r8*dp
   call initjas
   call getphi(x,u,dum1,dum2,dum3,dum4)
   psi2=u
   bmcmillan=bmcmillan+dp
   call initjas
   dpsi(1)=(psi1-psi2)/(2.0_r8*dp)
   end subroutine getjasder

end module potential
