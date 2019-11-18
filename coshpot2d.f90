module potential
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, save, allocatable :: nbin(:),spin(:)
   integer(kind=i4), private, save :: npart,ntab,ndim
   real(kind=r8), private, save :: hbar,v0,amu,el,range,scale,eli
   real(kind=r8), private, save, allocatable :: utab(:,:),dutab(:,:),d2utab(:,:)
   real(kind=r8), private, save, allocatable :: vtab(:,:)
   real(kind=r8), private, save, allocatable :: mass(:)
   real(kind=r8), private, save :: cjas(2),djas,v0op,v0eq
   logical :: iop,ieq
contains
   subroutine setphi(ndin,ntabin,nbinin,elin,massin,hbin)
   use rungekutta
   use mympi
   real(kind=r8) :: elin,massin(2),hbin
   real(kind=r8) :: dr,r,v0tmp,el2
   real(kind=r8), allocatable, dimension(:) :: f,df,d2f
   integer(kind=i4) :: ntabin,nbinin(:),i,ntab0,ic,ndin
   integer(kind=i4) :: jfac
   integer(kind=i4), parameter :: jntab=10000000
   real(kind=r8) :: jdr
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
      read (5,*) v0
      read (5,*) amu  !potential width parameter
      read (5,*) cjas(1) !opposite spin jastrow strength
      read (5,*) cjas(2) !equal spin jastrow strength
      read (5,*) djas    !jastrow healing distance
      read (5,*) v0op    !effective potential for opposite spin Jastrow
      read (5,*) v0eq    !effective potential for equal spin  Jastrow
      if (v0op.eq.0.0_r8) v0op=v0
      if (v0eq.eq.0.0_r8) v0eq=v0
      write (6,'(''!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!'')')
      write (6,'(''s-wave cosh potential parameters:'')')
!     write (6,'(''s-wave 2gauss potential parameters:'')')
!     write (6,'(''s-wave Morse potential parameters:'')')
      write (6,'(''V0 of the interaction ='',t40,f10.5)') v0
      write (6,'(''amu of the interaction ='',t40,f10.5)') amu
! wavefunction parameters:
      write (6,'(''!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!'')')
      write (6,'(''Jastrow parameters:'')')
      write (6,'(''opposite spin strength ='',t40,f10.5)') cjas(1)
      write (6,'(''equal spin strength ='',t40,f10.5)') cjas(2)
      write (6,'(''healing distance ='',t40,f10.5)') djas
      write (6,'(''opposite spin interaction strength='',t40,f10.5)') v0op
      write (6,'(''equal spin interaction strength='',t40,f10.5)') v0eq
   endif
   call bcast(v0)
   call bcast(amu)
   call bcast(cjas)
   call bcast(djas)
   call bcast(v0op)
   call bcast(v0eq)
   allocate(spin(npart),mass(npart))
! assign spin 1 to up particles and -1 to down ones
   spin(1:nbin(1))=1
   spin(nbin(1)+1:nbin(1)+nbin(2))=-1
   mass(1:nbin(1))=massin(1)
   mass(nbin(1)+1:npart)=massin(2)
   el2=0.5_r8*el
   if (el.ne.0.0_r8) then
      range=0.5_r8*el
   else
      range=djas
   endif
   dr=range/ntab
   scale=1.0_r8/dr
   iop=.false.
   ieq=.false.
   allocate(vtab(2,0:ntab))
   do i=0,ntab
      r=i*dr
      vtab(1,i)=vpot(r)
      vtab(2,i)=0.0_r8
!     vtab(2,i)=vpot(r)
!     vtab(1,i)=0.0_r8
      iop=.true.
   enddo
   ntab0=ntab*djas/el2
   if (el.eq.0.0_r8) ntab0=ntab*djas/range
   allocate(utab(2,0:ntab),dutab(2,0:ntab),d2utab(2,0:ntab))
   allocate(f(0:jntab),df(0:jntab),d2f(0:jntab))
! use an effective interaction for Jastrow
   v0tmp=v0
   if (cjas(1).ne.0.0_r8) then
      if (ndim.eq.1) then
         if (myrank().eq.0) write(6,'(''numerical jastrow working only in 3D'')')
         stop
      endif
      v0=v0op
      jdr=djas/jntab
      call orbital(f,df,d2f,hbar,jntab,jdr,vpot,0,ndim)
      do i=0,ntab0
         jfac=nint(i*real(jntab)/ntab0)
         utab(1,i)=f(jfac)/f(nint(ntab0*real(jntab)/ntab0))
         dutab(1,i)=df(jfac)/f(nint(ntab0*real(jntab)/ntab0))
         d2utab(1,i)=d2f(jfac)/f(nint(ntab0*real(jntab)/ntab0))
      end do
      if (myrank().eq.0) write(6,*) '---------------------------------------'
      if (myrank().eq.0) write(6,*) 'rungekutta 2D ntab     =',jntab
      if (myrank().eq.0) write(6,*) 'renormalization factor =',real(jntab)/ntab0
      if (myrank().eq.0) write(6,*) '---------------------------------------'
   endif
   if (cjas(2).ne.0.0_r8) then
      if (ndim.eq.1) then
         if (myrank().eq.0) write(6,'(''numerical jastrow working only in 3D'')')
         stop
      endif
      v0=v0eq
      call orbital(f,df,d2f,hbar,ntab0,dr,vpot,0,ndim)
      utab(2,0:ntab0)=f(0:ntab0)/f(ntab0)
      dutab(2,0:ntab0)=df(0:ntab0)/f(ntab0)
      d2utab(2,0:ntab0)=d2f(0:ntab0)/f(ntab0)
   endif
   v0=v0tmp
   if (ntab0.lt.ntab) then
      utab(:,ntab0+1:ntab)=1.0_r8
      dutab(:,ntab0+1:ntab)=0.0_r8
      d2utab(:,ntab0+1:ntab)=0.0_r8
   endif
   do ic=1,2
      if (cjas(ic).ne.0.0_r8) then
         dutab(ic,:)=dutab(ic,:)/utab(ic,:)
         d2utab(ic,:)=d2utab(ic,:)/utab(ic,:)-dutab(ic,:)**2
         utab(ic,:)=log(utab(ic,:))
      endif
      utab(ic,:)=cjas(ic)*utab(ic,:)
      dutab(ic,:)=cjas(ic)*dutab(ic,:)
      d2utab(ic,:)=cjas(ic)*d2utab(ic,:)
   enddo
   if (cjas(1).ne.0.0_r8) iop=.true.
   if (cjas(2).ne.0.0_r8) ieq=.true.
!  if (myrank().eq.0) then
!     do i=0,ntab
!        if (mod(i,10000).eq.0) write(78,'(f10.5,8e15.7)') &
!           i*dr,utab(1,i),dutab(1,i),d2utab(1,i),vtab(1,i), &
!           utab(2,i),dutab(2,i),d2utab(2,i),vtab(2,i)
!     enddo
!  endif
   end subroutine setphi

   function vpot(r)
   real(kind=r8) :: vpot,r
!   vpot=0.0_r8
!   if (r.le.0.0025_r8) vpot=-hbar*v0*10.0_r8**4
   vpot=-2.0_r8*v0*amu**2*hbar/cosh(amu*r)**2
!  vpot=-2.0_r8*v0*amu**2*hbar*exp(-0.5_r8*(amu*r)**2)
!  vpot=2.0_r8*amu**2*v0*hbar*(4.0_r8*exp(-(amu*r)**2)-exp(-(0.5_r8*amu*r)**2))
!  vpot=-2.0_r8*amu**2*v0*hbar*(exp(-amu*r)-2.0_r8*exp(-2.0_r8*amu*r))
   return
   end function vpot

   subroutine getphi(x,u,is,du,d2u,v)
   real(kind=r8), dimension(ndim,npart) :: x,du
   real(kind=r8) :: v,u,d2u,vpot,uj,duj,d2uj
   integer(kind=i4) :: is,i,j,index,itab
   real(kind=r8), dimension(ndim) :: dx
   real(kind=r8) :: r,dr,c1,c2,c3,c4,mu
   real(kind=r8), parameter :: tiny=1.0e-14_r8
   u=0.0_r8
   du=0.0_r8
   d2u=0.0_r8
   v=0.0_r8
   do i=1,npart-1
      do j=i+1,npart
         if (ieq.and.spin(i).eq.spin(j).or.iop.and.spin(i).ne.spin(j)) then
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
               if (spin(i).ne.spin(j)) then
                  itab=1
               else
                  itab=2
               endif
               vpot=c1*vtab(itab,index-1)+c2*vtab(itab,index) &
                  +c3*vtab(itab,index+1)+c4*vtab(itab,index+2)
               uj=c1*utab(itab,index-1)+c2*utab(itab,index) &
                 +c3*utab(itab,index+1)+c4*utab(itab,index+2)
               duj=c1*dutab(itab,index-1)+c2*dutab(itab,index) &
                  +c3*dutab(itab,index+1)+c4*dutab(itab,index+2)
               d2uj=c1*d2utab(itab,index-1)+c2*d2utab(itab,index) &
                   +c3*d2utab(itab,index+1)+c4*d2utab(itab,index+2)
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

   subroutine getphione(x,drr,i,phiratio)
   real(kind=r8) :: x(:,:),drr(:),phiratio
   real(kind=r8) :: u,uj
   integer(kind=i4) :: i,j,index,itab
   real(kind=r8), dimension(ndim) :: dx
   real(kind=r8) :: r,c1,c2,c3,c4,dr
   real(kind=r8), parameter :: tiny=1.0e-14_r8
   u=0.0_r8
   do j=1,npart
      if (j.ne.i) then
         if (ieq.and.spin(i).eq.spin(j).or.iop.and.spin(i).ne.spin(j)) then
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
               if (spin(i).ne.spin(j)) then
                  itab=1
               else
                  itab=2
               endif
               uj=c1*utab(itab,index-1)+c2*utab(itab,index) &
                 +c3*utab(itab,index+1)+c4*utab(itab,index+2)
               u=u-uj
            endif
            dx(:)=x(:,i)-x(:,j)-drr(:)
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
               if (spin(i).ne.spin(j)) then
                  itab=1
               else
                  itab=2
               endif
               uj=c1*utab(itab,index-1)+c2*utab(itab,index) &
                 +c3*utab(itab,index+1)+c4*utab(itab,index+2)
               u=u+uj
            endif
         endif
      endif
   enddo
   phiratio=exp(u)
   end subroutine getphione

!  subroutine getjasparam(npjas,pjas)
!  integer(kind=i4) :: npjas
!  real(kind=r8), pointer :: pjas(:)
!  npjas=1  ! number of Jastrow parameters
!  allocate(pjas(npjas))
!  if (npjas.ne.0) pjas=0.0_r8
!  pjas(1)=p1 ...
!  end subroutine getjasparam

!  subroutine setjasparam(pjas)
!  real(kind=r8) :: pjas(:)
! p1=pjas(1) ...
!  end subroutine setjasparam

!  subroutine getjasder(x,dpsi)
!  real(kind=r8) :: x(:,:),dpsi(:)
!  dpsi=0.0_r8
!  end subroutine getjasder

end module potential
