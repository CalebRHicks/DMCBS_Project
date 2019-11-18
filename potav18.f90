module potential
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, save, allocatable :: nbin(:),spin(:)
   integer(kind=i4), private, save :: npart,ntab,ndim
   real(kind=r8), private, save :: hbar,el,range,scale,eli
   real(kind=r8), private, save, allocatable :: utab(:),dutab(:),d2utab(:)
   real(kind=r8), private, save, allocatable :: vtab(:)
   real(kind=r8), private, save, allocatable :: mass(:)
   real(kind=r8), private, save :: cjas,djas
   logical :: iop,ieq
contains
   subroutine setphi(ndin,ntabin,nbinin,elin,massin,hbin)
   use rungekutta
   use mympi
   real(kind=r8) :: elin,massin(2),hbin
   real(kind=r8) :: dr,r,el2
   real(kind=r8), allocatable, dimension(:) :: f,df,d2f
   integer(kind=i4) :: ntabin,nbinin(:),i,ntab0,ndin
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
   hbar=2.0_r8*hbar !here I need hbar**2/m or hbar**2/2mu
   if (myrank().eq.0) then
      if (el.eq.0.0_r8) read (5,*) range
      read (5,*) cjas    !opposite spin jastrow strength
      read (5,*) djas    !jastrow healing distance
      write (6,'(''s-wave part of AV18 potential'')')
! wavefunction parameters:
      write (6,'(''Jastrow parameters:'')')
      if (el.eq.0.0_r8) write (6,'(''range ='',t40,f10.5)') range
      write (6,'(''opposite spin strength ='',t40,f10.5)') cjas   
      write (6,'(''healing distance ='',t40,f10.5)') djas
   endif
   if (el.eq.0.0_r8) call bcast(range)
   call bcast(cjas)
   call bcast(djas)
   allocate(spin(npart),mass(npart))
! assign spin 1 to up particles and -1 to down ones
   spin(1:nbin(1))=1
   spin(nbin(1)+1:nbin(1)+nbin(2))=-1
   mass(1:nbin(1))=massin(1)
   mass(nbin(1)+1:npart)=massin(2)
   el2=0.5_r8*el
   if (el.ne.0.0_r8) range=0.5_r8*el
   dr=range/ntab
   scale=1.0_r8/dr
   allocate(vtab(0:ntab))
   do i=0,ntab
      r=i*dr
      vtab(i)=vpot(r)
   enddo
   ntab0=ntab*djas/el2
   if (el.eq.0.0_r8) ntab0=ntab*djas/range
   allocate(utab(0:ntab),dutab(0:ntab),d2utab(0:ntab))
   allocate(f(0:ntab0),df(0:ntab0),d2f(0:ntab0))
! use an effective interaction for Jastrow
   if (cjas.ne.0.0_r8) then
      if (ndim.ne.3) then
         if (myrank().eq.0) write(6,'(''numerical jastrow working only in 3D'')')
         stop
      endif
      call orbital(f,df,d2f,hbar,ntab0,dr,vpot,0)
      utab(0:ntab0)=f(0:ntab0)/f(ntab0)
      dutab(0:ntab0)=df(0:ntab0)/f(ntab0)
      d2utab(0:ntab0)=d2f(0:ntab0)/f(ntab0)
   endif
   if (ntab0.lt.ntab) then
      utab(ntab0+1:ntab)=1.0_r8
      dutab(ntab0+1:ntab)=0.0_r8
      d2utab(ntab0+1:ntab)=0.0_r8
   endif
   if (cjas.ne.0.0_r8) then
      dutab(:)=dutab(:)/utab(:)
      d2utab(:)=d2utab(:)/utab(:)-dutab(:)**2
      utab(:)=log(utab(:))
   endif
   utab(:)=cjas*utab(:)
   dutab(:)=cjas*dutab(:)
   d2utab(:)=cjas*d2utab(:)
   if (myrank().eq.0) then
      do i=0,ntab
         if (myrank().eq.0) write(78,'(f10.5,5e15.7)') &
            i*dr,utab(i),dutab(i),d2utab(i),vtab(i)
      enddo
   endif
   end subroutine setphi

   function vpot(r)
   real(kind=r8) :: vpot,r,vv(2,2)
   call av18pw(1,1,1,2,1,-1,-1,r,vv)
   vpot=vv(1,1)
   return
   end function vpot

   subroutine getphi(x,u,is,du,d2u,v)
   real(kind=r8), dimension(ndim,npart) :: x,du
   real(kind=r8) :: v,u,d2u,vpot,uj,duj,d2uj
   integer(kind=i4) :: is,i,j,index
   real(kind=r8), dimension(ndim) :: dx
   real(kind=r8) :: r,dr,c1,c2,c3,c4,mu
   real(kind=r8), parameter :: tiny=1.0e-14_r8
   u=0.0_r8
   du=0.0_r8
   d2u=0.0_r8
   v=0.0_r8
   do i=1,npart-1
      do j=i+1,npart
         if (spin(i).ne.spin(j)) then
            dx(:)=x(:,i)-x(:,j)
            dx=dx-el*nint(dx*eli)
            r=max(tiny,sqrt(sum(dx**2)))
            if (r.lt.range) then
               dr=scale*r
               index=dr
               index=max(1,min(index,ntab-2))
               dr=dr-index
               c1=-dr*(dr-1.0_r8)*(dr-2.0_r8)/6.0_r8
               c2=(dr+1.0_r8)*(dr-1.0_r8)*(dr-2.0_r8)/2.0_r8
               c3=-(dr+1.0_r8)*dr*(dr-2.0_r8)/2.0_r8
               c4=(dr+1.0_r8)*dr*(dr-1.0_r8)/6.0_r8
               vpot=c1*vtab(index-1)+c2*vtab(index) &
                  +c3*vtab(index+1)+c4*vtab(index+2)
               uj=c1*utab(index-1)+c2*utab(index) &
                 +c3*utab(index+1)+c4*utab(index+2)
               duj=c1*dutab(index-1)+c2*dutab(index) &
                  +c3*dutab(index+1)+c4*dutab(index+2)
               d2uj=c1*d2utab(index-1)+c2*d2utab(index) &
                   +c3*d2utab(index+1)+c4*d2utab(index+2)
               v=v+vpot
               u=u+uj
               du(:,i)=du(:,i)+dx(:)*duj/(r*sqrt(mass(i)))
               du(:,j)=du(:,j)-dx(:)*duj/(r*sqrt(mass(j)))
               mu=mass(i)*mass(j)/(mass(i)+mass(j))
               d2u=d2u+(d2uj+(ndim-1.0_r8)*duj/r)/mu
            endif
         endif
      enddo
   enddo
   is=1
   d2u=d2u+sum(du*du)
   end subroutine getphi
end module potential
