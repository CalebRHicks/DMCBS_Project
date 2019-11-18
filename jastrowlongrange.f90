module jaslong
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, save :: numq,npart,nk
   integer(kind=i4), private, save, allocatable :: numsh(:)
   real(kind=r8), private, save, pointer :: ak(:,:),ak2(:)
   real(kind=r8), private, save :: beta,gamma
contains
   subroutine setjaslong(npartin,el,numqin,betain,gammain)
   use kshell
   integer(kind=i4) :: npartin,numqin
   real(kind=r8) :: el,betain,gammain
   npart=npartin
   numq=numqin
   beta=betain
   gamma=gammain
   allocate(numsh(numq))
   call setupk(el,numq,numsh,ak,ak2,nk,3)
   end subroutine setjaslong

   subroutine getjaslong(x,u,du,d2u)
   real(kind=r8) :: x(:,:),u,du(:,:),d2u
   real(kind=r8) :: dx(3),qmod,fact,arg,c,s
   integer(kind=i4) :: i,j,iq
   u=0.0_r8
   du=0.0_r8
   d2u=0.0_r8
   do i=1,npart-1
      do j=i+1,npart
         dx(:)=x(:,i)-x(:,j)
         do iq=2,nk  ! do not include q=0,0,0 state
            qmod=sqrt(ak2(iq))
            fact=gamma*exp(-beta*qmod)/qmod
            arg=sum(ak(:,iq)*dx(:))
            c=2.0_r8*cos(arg)
            s=2.0_r8*sin(arg)
            u=u+fact*c
            du(:,i)=du(:,i)-fact*ak(:,iq)*s
            du(:,j)=du(:,j)+fact*ak(:,iq)*s
            d2u=d2u-2.0_r8*fact*ak2(iq)*c
         enddo
      enddo
    enddo
    d2u=d2u+sum(du*du)
    end subroutine getjaslong

   subroutine getjasparam(npjas,pjas)
   integer(kind=i4) :: npjas
   real(kind=r8), pointer :: pjas(:)
   npjas=2  ! number of Jastrow parameters
   allocate(pjas(npjas))
   pjas(1)=beta
   pjas(2)=gamma
   end subroutine getjasparam

   subroutine setjasparam(pjas)
   real(kind=r8) :: pjas(:)
   beta=pjas(1)
   gamma=pjas(2)
   end subroutine setjasparam

   subroutine getjasder(x,dpsi)
   real(kind=r8) :: x(:,:),dpsi(:)
   real(kind=r8) :: u,du(3,npart),d2u,psi1,psi2
   real(kind=r8), parameter :: dp=0.0001_r8
   dpsi=0.0_r8
   beta=beta+dp
   call getjaslong(x,u,du,d2u)
   psi1=u
   beta=beta-2.0_r8*dp
   call getjaslong(x,u,du,d2u)
   psi2=u
   beta=beta+dp
   call getjaslong(x,u,du,d2u)
   dpsi(1)=(psi1-psi2)/(2.0_r8*dp)
   dpsi(2)=u/gamma
   end subroutine getjasder
end module jaslong
