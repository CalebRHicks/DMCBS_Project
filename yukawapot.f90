module potential
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, save :: npart,ntab,nk,ndim
   real(kind=r8), private, save :: el,range,scale,el2
   real(kind=r8), private, save :: v0,lambda,a,b,vtail,gau
   real(kind=r8), private, allocatable, save :: vtab(:)
   real(kind=r8), private, save, allocatable :: utab(:),dutab(:),d2utab(:)
   real(kind=r8), private, save, allocatable :: xsite(:,:)
contains
   subroutine setphi(ndin,ntabin,nbinin,elin,massin,hbin)
   use mympi
   use lattice
   real(kind=r8) :: elin,massin(2),hbin
   real(kind=r8) :: r,dr,pi,r0
   integer(kind=i4) :: ntabin,nbinin(:),i,ndin,z
   real(kind=r8), parameter :: alpha=0.007297352568_r8
   if (ndin.ne.3) then
      if (myrank().eq.0) write(6,'(''This code is for 3D yukawa, stopping'')')
      stop
   endif
   ndim=3
   ntab=ntabin
   el=elin
   npart=sum(nbinin)
   allocate(xsite(3,npart))
   xsite=(bestcubic(npart)-0.5_r8)*el
   xsite=xsite-el*nint(xsite/el)
   pi=4.0_r8*atan(1.0_r8)
   r0=el*(3.0_r8/(4.0_r8*pi*npart))**(1.0_r8/3.0_r8)
   if (myrank().eq.0) then
      read(5,*) z 
      read(5,*) lambda
      read(5,*) a
      read(5,*) b
      read(5,*) gau
      v0=z**2*alpha
      write(6,'(''Yukawa potential'')')
      write(6,'(''v0 ='',t40,f10.5)') v0
      write(6,'(''lambda ='',t40,f10.5)') lambda
      write(6,'(''Jastrow parameters:'')')
      write(6,'(''a ='',t40,f10.5)') a
      write(6,'(''b ='',t40,f10.5)') b
      write(6,'(''gaussian ='',t40,f10.5)') gau
   endif
   call bcast(v0)
   call bcast(lambda)
   call bcast(a)
   call bcast(b)
   el2=0.5_r8*el
   range=0.5_r8*el
   dr=range/ntab
   scale=1.0_r8/dr
   allocate(vtab(0:ntab))
   do i=0,ntab
      r=i*dr
      vtab(i)=vpot(r)
   enddo
   allocate(utab(ntab),dutab(ntab),d2utab(ntab))
   call jastro(utab,dutab,d2utab)
   if (myrank().eq.0) then
      do i=1,ntab
         if (myrank().eq.0) write(78,'(f10.5,3e15.7)') &
            i*dr*0.5_r8,utab(i),dutab(i),d2utab(i)
      enddo
   endif
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
   real(kind=r8) :: ex1
   if (r.gt.el2) then
      u=0.0_r8
      up=0.0_r8
      upp=0.0_r8
   else
      ex1=exp(-b*r)
      u=-a*ex1/r
      up=a*b*ex1/r+a*ex1/r**2
      upp=-a*b**2*ex1/r-2.0_r8*a*b*ex1/r**2-2.0_r8*a*ex1/r**3
   endif
   end subroutine getjas

   subroutine getphi(x,u,is,du,d2u,v)
   real(kind=r8), dimension(ndim,npart) :: x,du
   real(kind=r8) :: v,u,d2u,uj,duj,d2uj
   integer(kind=i4) :: is,i,j,index
   real(kind=r8), dimension(ndim) :: dx
   real(kind=r8) :: r,dr,c1,c2,c3,c4
   real(kind=r8) :: gauss,dgauss(ndim,npart),d2gauss
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
   d2u=d2u+sum(du*du)
   v=v+vtail
   gauss=0.0_r8
   dgauss=0.0_r8
   d2gauss=0.0_r8
   do i=1,npart
      dx=x(:,i)-xsite(:,i)
      dx=dx-nint(dx/el)*el
      gauss=gauss-sum(dx*dx)
      dgauss(:,i)=-2.0_r8*gau*dx
   enddo
   gauss=gau*gauss
   d2gauss=-6.0_r8*npart*gau**2
   u=u+gauss
   d2u=d2u+d2gauss+2.0_r8*sum(du*dgauss)
   du=du+dgauss
   end subroutine getphi

   function vpot(r)
   real(kind=r8) :: vpot,r,expo
   real(kind=r8), parameter :: small=1e-6_r8
   r=max(r,small)
   expo=exp(-r/lambda)
   vpot=v0*expo/r
   end function vpot

   subroutine vpottail(tail,npart,rho)
   real(kind=r8) :: tail,rho,pi
   integer(kind=i4) :: npart
   el2=0.5_r8*(npart/rho)**(1.0_r8/3.0_r8)
   pi=4.0_r8*atan(1.0_r8)
!  tail=4.0_r8*pi*eps*lambda*(lambda+el2)*exp(-el2/lambda)
!  tail=tail*0.5_r8*rho*npart
tail=0.0_r8
   end subroutine vpottail
end module potential
