module potential
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, save, allocatable :: nbin(:),spin(:)
   integer(kind=i4), private, save :: npart,ntab,ndim
   real(kind=r8), private, save :: hbar,v0,amu,el,eli
   real(kind=r8), private, save, allocatable :: mass(:)
   real(kind=r8), private, save :: a,w
contains
   subroutine setphi(ndin,ntabin,nbinin,elin,massin,hbin)
   use rungekutta
   use mympi
   real(kind=r8) :: elin,massin(2),hbin,djas,f1,f2
   integer(kind=i4) :: ntabin,nbinin(:),ndin
   ndim=ndin
   ntab=ntabin
   allocate(nbin(2))
   nbin=nbinin
   el=elin
   if (el.ne.0.0_r8) then
      eli=1.0_r8/el
   else
      eli=0.0_r8
   endif
   npart=sum(nbin)
   hbar=hbin
   if (myrank().eq.0) then
      read (5,*) v0
      read (5,*) amu  !potential width parameter
      write (6,'(''v3 potential parameters:'')')
      write (6,'(''V0 of the interaction ='',t40,f10.5)') v0
      write (6,'(''Rc of the interaction ='',t40,f10.5)') amu
   endif
   call bcast(v0)
   call bcast(amu)
   allocate(spin(npart),mass(npart))
! assign spin 1 to up particles and -1 to down ones
   spin(1:nbin(1))=1
   spin(nbin(1)+1:nbin(1)+nbin(2))=-1
   mass(1:nbin(1))=massin(1)
   mass(nbin(1)+1:npart)=massin(2)
   djas=0.5_r8*el
   w=1e-8
   f1=w-1.0_r8/djas*tan(w*(djas-amu))
   do while(.true.)
      w=w+1.0e-8_r8
      f2=w-1.0_r8/djas*tan(w*(djas-amu))
      if (f1*f2.lt.0.0_r8) exit
      f1=f2
   enddo
   a=djas/sin(w*(djas-amu))
   end subroutine setphi

   function vpot(r)
   real(kind=r8) :: vpot,r
   real(kind=r8), parameter :: huge=1e30_r8
   if (r.lt.amu) then
      vpot=huge
   else
      vpot=0.0_r8
   endif
   return
   end function vpot

   subroutine getphi(x,u,is,du,d2u,v)
   real(kind=r8), dimension(ndim,npart) :: x,du
   real(kind=r8) :: v,u,d2u
   integer(kind=i4) :: is,i,j,k
   real(kind=r8), dimension(ndim) :: dxij,dxjk,dxik
   real(kind=r8) :: r,s,c,uj,duj,d2uj
   real(kind=r8), parameter :: verynegative=-1e37_r8
   u=0.0_r8
   du=0.0_r8
   d2u=0.0_r8
   v=0.0_r8
   do i=1,npart
      do j=1,npart
         do k=1,npart
            if (i.ne.j.and.i.ne.k.and.j.ne.k) then
               r=0.0_r8
               dxij(:)=x(:,i)-x(:,j)
               dxij=dxij-el*nint(dxij*eli)
               dxik(:)=x(:,i)-x(:,k)
               dxik=dxik-el*nint(dxik*eli)
               dxjk(:)=x(:,j)-x(:,k)
               dxjk=dxjk-el*nint(dxjk*eli)
               r=sqrt(sum(dxij**2+dxik**2+dxjk**2)/3.0_r8)
               v=v+vpot(r)
               if (r.lt.amu) then
                  uj=verynegative
                  duj=0.0_r8
                  d2uj=0.0_r8
               else
                  if (r.lt.0.5_r8*el) then
                     s=sin(w*(r-amu))
                     c=cos(w*(r-amu))
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
               u=u+uj
               du(:,i)=du(:,i)+duj*(dxij+dxik)/r
               d2u=d2u+d2uj/3.0_r8+2.0_r8*duj/r
            endif
         enddo
      enddo
   enddo
   is=1
   end subroutine getphi
end module potential
