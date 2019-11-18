module potential
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, save, allocatable :: nbin(:),spin(:)
   integer(kind=i4), private, save :: npart,ntab,ndim
   real(kind=r8), private, save :: el,eli
   real(kind=r8), private, save :: r0,w,a,djas,cjas
contains
   subroutine setphi(ndin,ntabin,nbinin,elin,massin,hbin)
   use mympi
   integer(kind=i4) :: ndin,ntabin,nbinin(:)
   real(kind=r8) :: elin,massin(:),hbin
   real(kind=r8) :: f1,f2
   real(kind=r8), parameter :: small=1e-8
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
   djas=0.5_r8*el
   if (myrank().eq.0) then
      read (5,*) r0     ! range of the potential
      read (5,*) cjas   ! opposite spin jastrow strength
      write (6,'(''!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!'')')
      write (6,'(''barrier potential parameter:'')')
      write (6,'(''r0 = ''t40,f10.5)') r0
      write (6,'(''Jastrow parameters:'')')
      write (6,'(''healing distance ='',t40,f10.5)') djas
   endif
   call bcast(r0)
   call bcast(cjas)
   w=small
   f1=w-1.0_r8/djas*tan(w*(djas-r0))
   do while(.true.)
      w=w+1.0e-8_r8
      f2=w-1.0_r8/djas*tan(w*(djas-r0))
      if (f1*f2.lt.0.0_r8) exit
      f1=f2
   enddo
   a=djas/sin(w*(djas-r0))
   allocate(spin(npart))
   spin(1:nbin(1))=1
   spin(nbin(1)+1:nbin(1)+nbin(2))=-1
   end subroutine setphi

   function vpot(r)
   real(kind=r8) :: vpot,r
   real(kind=r8), parameter :: huge=1e30_r8
   if (r.le.r0) then
      vpot=huge
      return
   endif
   vpot=0.0_r8
   return
   end function vpot

   subroutine getphi(x,u,is,du,d2u,v)
   real(kind=r8), dimension(ndim,npart) :: x,du
   real(kind=r8) :: v,u,d2u,uj,duj,d2uj
   integer(kind=i4) :: is,i,j
   real(kind=r8), dimension(ndim) :: dx
   real(kind=r8) :: r,mu,s,c
   real(kind=r8), parameter :: verynegative=-1e37_r8
   u=0.0_r8
   du=0.0_r8
   d2u=0.0_r8
   v=0.0_r8
   do i=1,npart-1
      do j=i+1,npart
         if (spin(i).ne.spin(j)) then
            dx(:)=x(:,i)-x(:,j)
            dx=dx-el*nint(dx*eli)
            r=sqrt(sum(dx**2))
            if (r.le.r0) then
               uj=verynegative
               duj=0.0_r8
               d2uj=0.0_r8
            else 
               if (r.lt.djas) then
                  s=sin(w*(r-r0))
                  c=cos(w*(r-r0))
                  uj=a*s/r
                  duj=a*(w*c/r-s/r**2)
                  d2uj=a*(-w**2*s/r-2.0_r8*w*c/r**2+2.0_r8*s/r**3)
                  duj=duj/uj
                  d2uj=d2uj/uj-duj**2
                  uj=log(uj)
               else 
                  uj=0.0_r8
                  duj=0.0_r8
                  d2uj=0.0_r8
               endif
            endif
            uj=cjas*uj
            duj=cjas*duj
            d2uj=cjas*d2uj
            u=u+uj
            du(:,i)=du(:,i)+dx(:)*duj/r
            du(:,j)=du(:,j)-dx(:)*duj/r
            d2u=d2u+d2uj+(ndim-1.0_r8)*duj/r
            v=v+vpot(r)
         endif
      enddo
   enddo
   is=1
   d2u=2.0_r8*d2u+sum(du*du)
   end subroutine getphi

   subroutine getphione(x,drr,i,phiratio)
   real(kind=r8) :: x(:,:),drr(:),phiratio
   integer(kind=i4) :: i
! fake, subroutine to be implemented
   end subroutine getphione
end module potential
