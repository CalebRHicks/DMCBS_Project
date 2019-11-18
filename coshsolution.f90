module coshsolution
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   real(kind=r8), private, save :: kbes,amu,el2
   real(kind=r8), private, save, allocatable :: phi(:),dphi(:),d2phi(:)
   integer(kind=i4), private, save :: ntab,nsh
contains 
   subroutine initcoshfun(ntabin,el2in,kbesout)
   use mympi
   real(kind=r8) :: el2in,kbesout
   integer(kind=i4) :: ntabin
   if (myrank().eq.0) then
      read (5,*) kbes
      read (5,*) amu
      read (5,*) nsh
      write (6,'(''kbes in coshfun ='',t40,f10.5)') kbes
      write (6,'(''amu in coshfun ='',t40,f10.5)') amu
      write (6,'(''number of shell in Ewald sum ='',t40,i10)') nsh
   endif
   call bcast(kbes)
   call bcast(amu)
   call bcast(nsh)
   kbesout=kbes
   ntab=ntabin
   el2=el2in
   allocate(phi(0:ntab),dphi(0:ntab),d2phi(0:ntab))
   call setcoshfun(kbes)
   end subroutine initcoshfun

   subroutine setcoshfun(kbesin)
   use ewald
   real(kind=r8) :: kbesin
   real(kind=r8) :: r,dr,arg,s,c,th
   integer(kind=i4) :: i
   kbes=kbesin
   dr=el2/ntab
   phi(0)=amu*(kbes**2+1.0_r8)/kbes
   dphi(0)=0.0_r8
   d2phi(0)=-amu**3*(kbes**4+3.0_r8*kbes**2+2.0_r8)/(3.0_r8*kbes)
   do i=1,ntab
      r=i*dr
      arg=amu*kbes*r
      s=sin(arg)
      c=cos(arg)
      th=tanh(amu*r)
      phi(i)=(kbes*s+c*th)/(kbes*r)
      dphi(i)=(kbes**2*amu*c-amu*kbes*s*th+amu*c*(1.0_r8-th**2))/(kbes*r) &
        -(kbes*s+c*th)/(kbes*r**2)
      d2phi(i)=(-kbes**3*amu**2*s-kbes**2*amu**2*c*th-2.0_r8*amu**2*kbes*s*(1.0_r8-th**2) &
         -2.0_r8*amu**2*c*th*(1.0_r8-th**2))/(kbes*r) &
         -2.0_r8*(kbes**2*amu*c-amu*kbes*s*th+amu*c*(1.0_r8-th**2))/(kbes*r**2) &
         +2.0_r8*(kbes*s+c*th)/(kbes*r**3)
   enddo
   dphi=amu*dphi/phi(0)
   d2phi=amu*d2phi/phi(0)
   phi=amu*phi/phi(0)
   do while(abs(dphi(ntab)).gt.0.025_r8)
      phi=phi/1.1_r8
      dphi=dphi/1.1_r8
      d2phi=d2phi/1.1_r8
   enddo
!  do i=0,ntab
!     write(69,('(f10.5,3f28.10)')) i*dr,phi(i),dphi(i),d2phi(i)
!  enddo
   call setewald(nsh,ntab,el2,phi,dphi,d2phi)
!  call checkder(phi,dphi,d2phi,dr)
   end subroutine setcoshfun

   subroutine getcoshfun(x,ph,dph,d2ph)
   use ewald
   real(kind=r8) :: x(3),dph(3)
   real(kind=r8) :: ph,d2ph
   call phiewald(x,ph,dph,d2ph)
   end subroutine getcoshfun

   subroutine checkder(f,df,d2f,dr)
   real(kind=r8) :: f(:),df(:),d2f(:),dr
   integer(kind=i4) :: i
   do i=2,ntab-1
      write(10,*) i*dr,df(i),(f(i+1)-f(i-1))/(2.0_r8*dr)
      write(11,*) i*dr,d2f(i),(f(i-1)-2.0_r8*f(i)+f(i+1))/dr**2
   enddo
   end subroutine
end module coshsolution
