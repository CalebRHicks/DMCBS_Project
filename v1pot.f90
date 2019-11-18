module potential
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, save, allocatable :: nbin(:),spin(:)
   integer(kind=i4), private, save :: npart,ndim
   real(kind=r8), private, save :: el,eli
   real(kind=r8), private, save :: range,a,b,cjas
   real(kind=r8), private, save, allocatable :: mass(:)
contains
   subroutine setphi(ndin,ntabin,nbinin,elin,massin,hbin)
   use mympi
   real(kind=r8) :: elin,massin(2),hbin
   real(kind=r8) :: dummy,r,f,df,d2f,ex
real(kind=r8) :: norm,d0
   integer(kind=i4) :: ntabin,nbinin(:),ndin,i
   dummy=ntabin
   dummy=hbin
   ndim=ndin
   allocate(nbin(2))
   nbin=nbinin
   el=elin
   if (el.ne.0.0_r8) then
      eli=1.0_r8/el
   else
      eli=0.0_r8
   endif
   npart=sum(nbin)
   if (myrank().eq.0) then
      read (5,*) cjas !jastrow strength
      read (5,*) a
      read (5,*) b
      write (6,'(''toy-model potential'')')
      write (6,'(''Gaussian Jastrow'')')
      write (6,'(''Jastrow strength ='',t40,f10.5)') cjas
      write (6,'(''Jastrow a parameter ='',t40,f10.5)') a
      write (6,'(''Jastrow b parameter ='',t40,f10.5)') b
   endif
   call bcast(cjas)
   call bcast(a)
   call bcast(b)
   allocate(spin(npart),mass(npart))
! assign spin 1 to up particles and -1 to down ones
   spin(1:nbin(1))=1
   spin(nbin(1)+1:nbin(1)+nbin(2))=-1
   mass(1:nbin(1))=massin(1)
   mass(nbin(1)+1:npart)=massin(2)
   range=0.5_r8*el
b=-40.0_r8
d0=0.5_r8
a=-2.0_r8*d0**2*b
norm=1.0_r8-a*d0**2-b*d0**4
range=d0
   do i=0,1000
      r=i*range/1000
      ex=exp(-b*r**2)
!     f=1.0_r8-a*ex
!     df=2.0_r8*a*b*r*ex
!     d2f=-2.0_r8*a*b*ex*(2.0_r8*b*r**2-1.0_r8)
f=-(1.0_r8-a*r**2-b*r**4)+norm
df=2.0_r8*(a+2.0_r8*b*r**2)*r
d2f=2.0_r8*(a+6.0_r8*b*r**2)
      write (91,'(5e25.15)') r,f,df,d2f,2993.0_r8*exp(-6.1_r8*r**2)-338.0_r8*exp(-1.17_r8*r**2)
   enddo
   end subroutine setphi

   subroutine getphi(x,u,is,du,d2u,v)
   real(kind=r8), dimension(ndim,npart) :: x,du
   real(kind=r8) :: v,u,d2u,uj,duj,d2uj
   integer(kind=i4) :: is,i,j
   real(kind=r8), dimension(ndim) :: dx
   real(kind=r8) :: r,mu,ex,f,df,d2f
real(kind=r8) :: b,d0,a,norm
   u=0.0_r8
   du=0.0_r8
   d2u=0.0_r8
   v=0.0_r8
b=-40.0_r8
d0=0.5_r8
a=-2.0_r8*d0**2*b
norm=1.0_r8-a*d0**2-b*d0**4
   do i=1,npart-1
      do j=i+1,npart
         dx(:)=x(:,i)-x(:,j)
         dx=dx-el*nint(dx*eli)
         r=sqrt(sum(dx**2))
         if (r.lt.range) then
            if (cjas.ne.0.0_r8) then
               ex=exp(-b*r**2)
!              f=1.0_r8-a*ex
!              df=2.0_r8*a*b*r*ex
!              d2f=2.0_r8*a*b*ex*(1.0_r8-2.0_r8*b*r**2)
     f=-(1.0_r8-a*r**2-b*r**4)+norm
     df=2.0_r8*(a+2.0_r8*b*r**2)*r
     d2f=2.0_r8*(a+6.0_r8*b*r**2)

               uj=log(f)
               duj=df/f
               d2uj=d2f/f-duj**2
               u=u+uj
               du(:,i)=du(:,i)+dx(:)*duj/(r*sqrt(mass(i)))
               du(:,j)=du(:,j)-dx(:)*duj/(r*sqrt(mass(j)))
               mu=mass(i)*mass(j)/(mass(i)+mass(j))
               d2u=d2u+(d2uj+(ndim-1.0_r8)*duj/r)/mu
            endif
!           v=v+2993.0_r8*exp(-6.1_r8*r**2)-338.0_r8*exp(-1.17_r8*r**2)
         endif
if (r.lt.0.5_r8*el) v=v+2993.0_r8*exp(-6.1_r8*r**2)-338.0_r8*exp(-1.17_r8*r**2)
      enddo
   enddo
   is=1
   d2u=d2u+sum(du*du)
   end subroutine getphi
end module potential
