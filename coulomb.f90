module potential
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, save :: npart,ntab,nk,ndim
   real(kind=r8), private, save :: el,range,scale,el2
   real(kind=r8), private, save :: constene,selfene,a,b,f
   real(kind=r8), private, allocatable, save :: vrtab(:),vktab(:)
   real(kind=r8), private, pointer, save :: ak(:,:)
   real(kind=r8), private, save, allocatable :: utab(:),dutab(:),d2utab(:)
contains
   subroutine setphi(ndin,ntabin,nbinin,elin,massin,hbin)
   use mympi
   use kshell
   real(kind=r8) :: elin,massin(2),hbin,dummy(2)
   real(kind=r8) :: r,pi,dr,um,fpv,ak2,rs
   real(kind=r8), pointer :: dummy2(:)
   integer(kind=i4) :: ntabin,nbinin(:),i,ndin,nsh
   integer(kind=i4), allocatable :: dummy1(:)
   dummy=massin  ! this is very stupd, but avoid a warning message
   if (ndin.ne.2) then
      if (myrank().eq.0) write(6,'(''This code is for 2D egas, stopping'')')
      stop
   endif
   ndim=2
   ntab=ntabin
   el=elin
   npart=sum(nbinin)
   if (myrank().eq.0) then
      read(5,*) nsh
      read(5,*) a
      read(5,*) b
      read(5,*) f 
      write(6,'(''number of shells in Ewald sums'',t40,i10)') nsh
      write(6,'(''Jastrow parameters:'')')
      write(6,'(''a ='',t40,f10.5)') a
      write(6,'(''b ='',t40,f10.5)') b
      write(6,'(''f ='',t40,f10.5)') f
   endif
   call bcast(nsh)
   call bcast(a)
   call bcast(b)
   call bcast(f)
   allocate(dummy1(nsh))
   el2=0.5_r8*el
   range=0.5_r8*el
   pi=4.0_r8*atan(1.0_r8)
   dr=range/ntab
   scale=1.0_r8/dr
   call setupk(el,nsh,dummy1,ak,dummy2,nk,ndim)
   allocate(vktab(nk),vrtab(0:ntab))
   um=el/5.0205_r8
   rs=1.0_r8/hbin**0.5_r8
   fpv=4.0_r8*pi/(el**2*rs)
   constene=-2.0_r8*um*sqrt(pi)/el**2*npart**2/rs
   selfene=-2.0_r8/(um*sqrt(pi))*npart/rs
   do i=2,nk
      ak2=sqrt(sum(ak(:,i)**2))
      vktab(i-1)=erfc(um*ak2*0.5_r8)/ak2*fpv
   enddo
   do i=0,ntab
      r=i*dr
      vrtab(i)=2.0_r8*derfc(r/um)/rs
   enddo
   allocate(utab(ntab),dutab(ntab),d2utab(ntab))
   call jastro(utab,dutab,d2utab)
!  if (myrank().eq.0) then
!     do i=1,ntab
!        if (myrank().eq.0) write(78,'(f10.5,3e15.7)') &
!           i*dr,utab(i),dutab(i),d2utab(i)
!     enddo
!  endif
   end subroutine setphi


   subroutine jastro(u,up,upp)
   use matrixmod
   real(kind=r8) :: u(:),up(:),upp(:)
   real(kind=r8) :: r
   real(kind=r8) :: rbar,r2,r3,mat(2,2),vec(2),dum,der
   real(kind=r8) :: rbmax,rbmin,shift
   integer(kind=i4) :: i,is,it
   real(kind=r8), parameter :: convergence=1e-8_r8
   real(kind=r8), parameter :: dr=1e-4_r8
   rbmax=el2
   rbmin=0.0_r8
   rbar=0.5_r8*(rbmax+rbmin)
   it=0
   do while (.true.)
      r3=rbar**3
      r2=rbar**2
      r=rbar
      mat(1,1)=r3-3.0_r8*el2*r2+3.0_r8*el2**2*r
      mat(1,2)=1.0_r8
      mat(2,1)=3.0_r8*r2-6.0_r8*el2*r+3.0_r8*el2**2
      mat(2,2)=0.0_r8
      call rmatinv(mat,dum,is,2)
      call getjas(rbar,vec(1),vec(2),dum)
      vec=matmul(mat,vec)
      der=6.0_r8*vec(1)*r-6.0_r8*vec(1)*el2
      if ((der-dum).lt.0.0_r8) then
         rbmax=rbar
      else
         rbmin=rbar
      endif
      if (abs(der-dum).lt.convergence) exit
      rbar=0.5_r8*(rbmax+rbmin)
      it=it+1
      if (it.gt.1000.or.rbar.lt.0.1_r8*el2) then
         write (6,'(''no convergence in jastrop, stopping execution'')')
         stop
      endif
   enddo
   do i=1,ntab
      r=i*el2/ntab
      if (r.le.rbar) then
         call getjas(r,u(i),up(i),upp(i))
         call getjas(r,u(i),up(i),upp(i))
         call getjas(r,u(i),up(i),upp(i))
      else
         u(i)=vec(1)*r*(r**2-3.0_r8*el2*r+3.0_r8*el2**2)+vec(2)
         up(i)=3.0_r8*vec(1)*(r**2-2.0_r8*el2*r+el2**2)
         upp(i)=3.0_r8*vec(1)*(2.0_r8*r-2.0_r8*el2)
     endif
     shift=vec(1)*el2*(el2**2-3.0_r8*el2**2+3.0_r8*el2**2)+vec(2)
     u(i)=u(i)-shift
   enddo
   end subroutine jastro

   subroutine getjas(r,u,up,upp)
   real(kind=r8) :: r,u,up,upp
   real(kind=r8) :: ex1,ex2,fac
   if (r.gt.el2) then
      u=0.0_r8
      up=0.0_r8
      upp=0.0_r8
   else
      ex1=exp(-b*r)
      ex2=exp(-r/f)
      fac=1.0_r8-ex2
      u=-a*ex1*fac/r
      up=a*b*ex1*fac/r+a*ex1*fac/r**2-a*ex1*ex2/(r*f)
      upp=-a*b**2*ex1*fac/r-2.0_r8*a*b*ex1*fac/r**2 &
         +2.0_r8*a*b*ex1*ex2/(r*f)-2.0_r8*a*ex1*fac/r**3 &
         +2.0_r8*a*ex1*ex2/(r**2*f)+a*ex1*ex2/(r*f**2)
   endif
   end subroutine getjas

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
   v=vpot(x)
   end subroutine getphi

   function vpot(x)
   real(kind=r8) :: x(:,:),vpot,r,dr,dd,rhok(2*nk),kr
   integer(kind=i4) :: idx,ik,ik2,ik1,i,j
   real(kind=r8) :: vr,vk,vrsum,vksum,dx(2)
! Real space potential
   vrsum=0.0_r8
   vksum=0.0_r8
   do i=2,npart
      do j=1,i-1
         dx=x(:,i)-x(:,j)
         dx=dx-el*nint(dx/el)
         r=sqrt(sum(dx(:)**2))
         dr=r*scale
         idx=dr
         dd=dr-idx
         if (idx.le.ntab-1) then
            vr=(vrtab(idx)*(1.0_r8-dd)+vrtab(idx+1)*dd)/r
            vrsum=vrsum+vr
         endif
      enddo
   enddo
   rhok=0.0_r8
   do ik=1,nk-1
      ik2=2*ik
      ik1=ik2-1
      do i=1,npart
         kr=ak(1,ik+1)*x(1,i)+ak(2,ik+1)*x(2,i)
         rhok(ik1)=rhok(ik1)+cos(kr)
         rhok(ik2)=rhok(ik2)+sin(kr)
      enddo
   enddo
! k-space potential
   do ik=1,nk-1
      ik2=ik*2
      ik1=ik2-1
      vk=vktab(ik)*(rhok(ik1)*rhok(ik1)+rhok(ik2)*rhok(ik2))
      vksum=vksum+vk
   enddo
   vpot=vrsum+vksum+constene+selfene
   return
   end function vpot

end module potential
