module potential
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   real(kind=r8), private, parameter :: epsilon = 10.22_r8
   real(kind=r8), private, parameter :: sigma = 2.556_r8
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
   real(kind=r8) :: r,dr,em1,em3,ufit,rho
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
      write(6,'(''Lennard-Jones potential for Helium'')')
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
   em1=mmcmillan+1
   em3=mmcmillan+3
   rfit=em1*el2/em3
   ufit=0.5d0*(bmcmillan/rfit)**mmcmillan
   alpha=-mmcmillan*ufit/(3.0d0*rfit*(rfit-el2)**2)
   u0=alpha*(rfit-el2)**3-ufit
   do i=0,ntab
      r=i*dr
      vtab(i)=vpot(r)
      call jastro(r,utab(i),dutab(i),d2utab(i))
   enddo
   if (myrank().eq.0) then
      do i=1,ntab
         if (myrank().eq.0) write(78,'(f10.5,3e15.7)') &
            i*dr*0.5_r8,utab(i),dutab(i),d2utab(i)
      enddo
   endif
   rho=npart/el**3
   call vpottail(vtail,npart,rho)
   if (myrank().eq.0) write(6,'(''tail of potential ='',t35,f15.10)') vtail
   end subroutine setphi

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
            uj=c1*utab(index-1)+c2*utab(index) &
              +c3*utab(index+1)+c4*utab(index+2)
            duj=c1*dutab(index-1)+c2*dutab(index) &
               +c3*dutab(index+1)+c4*dutab(index+2)
            d2uj=c1*d2utab(index-1)+c2*d2utab(index) &
                +c3*d2utab(index+1)+c4*d2utab(index+2)
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
   real(kind=r8) :: vpot,r,r6
   real(kind=r8), parameter :: small=1e-10_r8
   if (r.eq.0.0_r8) r=small
   r6=(sigma/r)**6
   vpot=4.0_r8*epsilon*r6*(r6-1.0_r8)
   end function vpot

   subroutine vpottail(tail,npart,rho)
   real(kind=r8) :: tail,rho,pi
   integer(kind=i4) :: npart
   el2=0.5_r8*(npart/rho)**(1.0_r8/3.0_r8)
   pi=4.0_r8*atan(1.0_r8)
   tail=8.0_r8*epsilon*pi*rho*npart*sigma**3* &
      (1.0_r8/(9.0_r8*(el2/sigma)**9)-1.0_r8/(3.0_r8*(el2/sigma)**3))
   end subroutine vpottail
end module potential
