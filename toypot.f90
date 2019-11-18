module potential
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, save, allocatable :: nbin(:),spin(:)
   integer(kind=i4), private, save :: npart,ntab,ndim
   real(kind=r8), private, save :: hbar,v0,el,range,scale,eli
   real(kind=r8), private, save, allocatable :: utab(:),dutab(:),d2utab(:)
   real(kind=r8), private, save, allocatable :: vtab(:)
   real(kind=r8), private, save, allocatable :: mass(:)
   real(kind=r8), private, save :: cjas,djas
contains
   subroutine setphi(ndin,ntabin,nbinin,elin,massin,hbin)
   use rungekutta
   use mympi
   real(kind=r8) :: elin,massin(2),hbin
   real(kind=r8) :: dr,r,quench
   integer(kind=i4) :: ntabin,nbinin(:),i,ndin
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
      read (5,*) cjas !jastrow strength
      read (5,*) djas    !jastrow healing distance
      if (djas.gt.0.5_r8*el) djas=0.5_r8*el
      read (5,*) quench  !quencher of potential
      write (6,'(''toy-model potential, 1/(1+x**2)'')')
! wavefunction parameters:
      write (6,'(''Jastrow parameters:'')')
      write (6,'(''Jastrow strength ='',t40,f10.5)') cjas
      write (6,'(''healing distance ='',t40,f10.5)') djas
      write (6,'(''interaction quencher ='',t40,f10.5)') quench
   endif
   call bcast(cjas)
   call bcast(djas)
   call bcast(quench)
   v0=quench
   allocate(spin(npart),mass(npart))
! assign spin 1 to up particles and -1 to down ones
   spin(1:nbin(1))=1
   spin(nbin(1)+1:nbin(1)+nbin(2))=-1
   mass(1:nbin(1))=massin(1)
   mass(nbin(1)+1:npart)=massin(2)
   if (el.ne.0.0_r8) range=djas
   dr=range/ntab
   scale=1.0_r8/dr
   allocate(vtab(0:ntab))
   do i=0,ntab
      r=i*dr
      vtab(i)=vpot(r)
   enddo
   allocate(utab(0:ntab),dutab(0:ntab),d2utab(0:ntab))
! use an effective interaction for Jastrow
   if (cjas.ne.0.0_r8) then
      if (ndim.ne.3) then
         if (myrank().eq.0) write(6,'(''numerical jastrow working only in 3D'')')
         stop
      endif
      v0=quench
      call orbital(utab,dutab,d2utab,hbar,ntab,dr,vpot,0)
      v0=1.0_r8
   endif
   if (cjas.ne.0.0_r8) then
 dutab=dutab/utab(ntab)
 d2utab=d2utab/utab(ntab)
 utab=utab/utab(ntab)
      dutab=dutab/utab
      d2utab=d2utab/utab-dutab**2
      utab=log(utab)
!     call polymatch3(utab,dutab,d2utab,ntab)
!     call polymatch(utab,dutab,d2utab,ntab)
   endif
   utab=cjas*utab
   dutab=cjas*dutab
   d2utab=cjas*d2utab
   if (myrank().eq.0) then
      do i=0,ntab
         r=i*dr
         if (myrank().eq.0) write(78,'(f10.5,5e15.7)') &
            i*dr,utab(i),dutab(i),d2utab(i),vpot(r)
      enddo
   endif
   v0=1.0_r8
   end subroutine setphi

   subroutine polymatch(f,df,d2f,notab)
! match f with u=a+b*r+c*r**2+d*r**3+e*r**4
! and solve to have f(rbar)=u(rbar) and the same for f' and f''
! at some distance h we have u'(h)=u''(h)=0
   use matrixmod
   integer(kind=i4) :: notab
   real(kind=r8) :: f(0:notab),df(0:notab),d2f(0:notab)
   real(kind=r8) :: r
   real(kind=r8) :: rbar,r1,r2,r3,r4,mat(3,3),vec(3),dum
   real(kind=r8) :: shift,dr,r0
   real(kind=r8) :: c1,c2,c3,c4
   integer(kind=i4) :: i,is,index
   r0=range
   r1=2.0_r8*r0
   rbar=0.75_r8*r0
   r4=rbar**4
   r3=rbar**3
   r2=rbar**2
   r=rbar
   mat(1,1)=1.0_r8
   mat(1,2)=3.0_r8*r0**2*r-3.0_r8*r0*r2+r3
   mat(1,3)=r1**3*r-3.0_r8*r1*r0*r2+r4
   mat(2,1)=0.0_r8
   mat(2,2)=3.0_r8*r0**2-3.0_r8*r1*r+3.0_r8*r2
   mat(2,3)=r1**3-3.0_r8*r1**2*r+4.0_r8*r3
   mat(3,1)=0.0_r8
   mat(3,2)=-3.0_r8*r1+6.0_r8*r
   mat(3,3)=-3.0_r8*r1**2+12.0_r8*r2
   call rmatinv(mat,dum,is,3)
   dr=scale*rbar
   index=dr
   dr=dr-index
   c1=-dr*(dr-1.0_r8)*(dr-2.0_r8)/6.0_r8
   c2=(dr+1.0_r8)*(dr-1.0_r8)*(dr-2.0_r8)/2.0_r8
   c3=-(dr+1.0_r8)*dr*(dr-2.0_r8)/2.0_r8
   c4=(dr+1.0_r8)*dr*(dr-1.0_r8)/6.0_r8
   vec(1)=c1*f(index-1)+c2*f(index) &
         +c3*f(index+1)+c4*f(index+2)
   vec(2)=c1*df(index-1)+c2*df(index) &
         +c3*df(index+1)+c4*df(index+2)
   vec(3)=c1*d2f(index-1)+c2*d2f(index) &
         +c3*d2f(index+1)+c4*d2f(index+2)
   vec=matmul(mat,vec)
   do i=1,notab
      r=i*r0/notab
      if (r.gt.rbar) then
         f(i)=vec(1)+(3.0_r8*vec(2)*r0**2+vec(3)*r1**3)*r    &
             +(-3.0_r8*vec(2)*r0-3.0_r8*vec(3)*r1*r0)*r**2  &
             +vec(2)*r**3+vec(3)*r**4
         df(i)=3.0_r8*vec(2)*r0**2+vec(3)*r1**3          &
              +(-3.0_r8*vec(2)*r1-3.0_r8*vec(3)*r1**2)*r  &
              +3.0_r8*vec(2)*r**2+4.0_r8*vec(3)*r**3
         d2f(i)=-3.0_r8*vec(2)*r1-3.0_r8*vec(3)*r1**2  &
               +6.0_r8*vec(2)*r+12.0_r8*vec(3)*r**2
     endif
   enddo
   shift=f(notab)
   f=f-shift
   end subroutine

   subroutine polymatch2(f,df,d2f,notab)
! match f with u=a+b*r+c*r**2+d*r**3+e*r**4+f*r**5
! and solve to have f(rbar)=u(rbar) and the same for f' and f''
! at some distance h we have u'(h)=u''(h)=u'''(h)=0
   use matrixmod
   integer(kind=i4) :: notab
   real(kind=r8) :: f(0:notab),df(0:notab),d2f(0:notab)
   real(kind=r8) :: r
   real(kind=r8) :: rbar,h,h2,h3,h4,mat(3,3),vec(3),dum
   real(kind=r8) :: shift,dr
   real(kind=r8) :: c1,c2,c3,c4
   integer(kind=i4) :: i,is,index
   rbar=0.75_r8*range
   h=range
   h4=h**4
   h3=h**3
   h2=h**2
   r=rbar
   mat(1,1)=1.0_r8
   mat(1,2)=-4.0_r8*h3*r+6.0_r8*h2*r**2-4.0_r8*h*r**3+r**4
   mat(1,3)=-15.0_r8*h4*r+20.0_r8*h3*r**2-10.0_r8*h2*r**3+r**5
   mat(2,1)=0.0_r8
   mat(2,2)=-4.0_r8*h3+12.0_r8*h2*r-12.0_r8*h*r**2+4.0_r8*r**3
   mat(2,3)=-15.0_r8*h4+40.0_r8*h3*r-30.0_r8*h2*r**2+5.0_r8*r**4
   mat(3,1)=0.0_r8
   mat(3,2)=12.0_r8*h2-24.0_r8*h*r+12.0_r8*r**2
   mat(3,3)=40.0_r8*h3-60.0_r8*h2*r+20.0_r8*r**3
   call rmatinv(mat,dum,is,3)
   dr=scale*rbar
   index=dr
   dr=dr-index
   c1=-dr*(dr-1.0_r8)*(dr-2.0_r8)/6.0_r8
   c2=(dr+1.0_r8)*(dr-1.0_r8)*(dr-2.0_r8)/2.0_r8
   c3=-(dr+1.0_r8)*dr*(dr-2.0_r8)/2.0_r8
   c4=(dr+1.0_r8)*dr*(dr-1.0_r8)/6.0_r8
   vec(1)=c1*f(index-1)+c2*f(index) &
         +c3*f(index+1)+c4*f(index+2)
   vec(2)=c1*df(index-1)+c2*df(index) &
         +c3*df(index+1)+c4*df(index+2)
   vec(3)=c1*d2f(index-1)+c2*d2f(index) &
         +c3*d2f(index+1)+c4*d2f(index+2)
   vec=matmul(mat,vec)
   do i=1,notab
      r=i*h/notab
      if (r.gt.rbar) then
         f(i)=vec(1)+(-4.0_r8*vec(2)*h3-15.0_r8*vec(3)*h4)*r  &
             +(6.0_r8*vec(2)*h2+20.0_r8*vec(3)*h3)*r**2       &
             +(-4.0_r8*vec(2)*h-10.0_r8*vec(3)*h2)*r**3         &
             +vec(2)*r**4+vec(3)*r**5
         df(i)=-4.0_r8*vec(2)*h3-15.0_r8*vec(3)*h4      &
              +(12.0_r8*vec(2)*h2+40.0_r8*vec(3)*h3)*r  &
              +(-12.0_r8*vec(2)*h-30.0_r8*vec(3)*h2)*r**2  &
              +4.0_r8*vec(2)*r**3+5.0_r8*vec(3)*r**4
         d2f(i)=12.0_r8*vec(2)*h2+40.0_r8*vec(3)*h3    &
               +(-24.0_r8*vec(2)*h-60.0_r8*vec(3)*h2)*r  &
               +12.0_r8*vec(2)*r**2+20.0_r8*vec(3)*r**3
     endif
   enddo
   shift=f(notab)
   f=f-shift
   end subroutine

   subroutine polymatch3(f,df,d2f,notab)
! match f with u=a+b*r+c*r**2+d*r**3+e*r**4+f*r**5
! and solve to have f(rbar)=u(rbar) and the same for f', f'' and f'''
! at some distance h we have u'(h)=u''(h)=0
   use matrixmod
   integer(kind=i4) :: notab
   real(kind=r8) :: f(0:notab),df(0:notab),d2f(0:notab)
   real(kind=r8) :: r
   real(kind=r8) :: rbar,h,h2,h3,h4,mat(4,4),vec(4),dum
   real(kind=r8) :: shift,dr
   real(kind=r8) :: c1,c2,c3,c4
   integer(kind=i4) :: i,is,index
   rbar=0.4_r8*range
   h=range
   h4=h**4
   h3=h**3
   h2=h**2
   r=rbar
   mat(1,1)=1.0_r8
   mat(1,2)=3.0_r8*h2*r-3.0_r8*h*r**2+r**3
   mat(1,3)=8.0_r8*h3*r-6.0_r8*h2*r**2+r**4
   mat(1,4)=15.0_r8*h4*r-10.0_r8*h3*r**2+r**5
   mat(2,1)=0.0_r8
   mat(2,2)=3.0_r8*h2-6.0_r8*h*r+3.0_r8*r**2
   mat(2,3)=8.0_r8*h3-12.0_r8*h2*r+4.0_r8*r**3
   mat(2,4)=15.0_r8*h4-20.0_r8*h3*r+5.0_r8*r**4
   mat(3,1)=0.0_r8
   mat(3,2)=-6.0_r8*h+6.0_r8*r
   mat(3,3)=-12.0_r8*h2+12.0_r8*r**2
   mat(3,4)=-20.0_r8*h3+20.0_r8*r**3
   mat(4,1)=0.0_r8
   mat(4,2)=6.0_r8
   mat(4,3)=24.0_r8*r
   mat(4,4)=60.0_r8*r**2
   call rmatinv(mat,dum,is,4)
   dr=scale*rbar
   index=dr
   dr=dr-index
   c1=-dr*(dr-1.0_r8)*(dr-2.0_r8)/6.0_r8
   c2=(dr+1.0_r8)*(dr-1.0_r8)*(dr-2.0_r8)/2.0_r8
   c3=-(dr+1.0_r8)*dr*(dr-2.0_r8)/2.0_r8
   c4=(dr+1.0_r8)*dr*(dr-1.0_r8)/6.0_r8
   vec(1)=c1*f(index-1)+c2*f(index) &
         +c3*f(index+1)+c4*f(index+2)
   vec(2)=c1*df(index-1)+c2*df(index) &
         +c3*df(index+1)+c4*df(index+2)
   vec(3)=c1*d2f(index-1)+c2*d2f(index) &
         +c3*d2f(index+1)+c4*d2f(index+2)
   vec(4)=(d2f(index+1)-d2f(index-1))/(2.0_r8*h/notab)
   vec=matmul(mat,vec)
   do i=1,notab
      r=i*h/notab
      if (r.gt.rbar) then
         f(i)=vec(1)+(3.0_r8*vec(2)*h2+8.0_r8*vec(3)*h3+15.0_r8*vec(4)*h4)*r  &
             +(-3.0_r8*vec(2)*h-6.0_r8*vec(3)*h2-10.0_r8*vec(4)*h3)*r**2      &
             +vec(2)*r**3+vec(3)*r**4+vec(4)*r**5
         df(i)=3.0_r8*vec(2)*h2+8.0_r8*vec(3)*h3+15.0_r8*vec(4)*h4      &
              +(-6.0_r8*vec(2)*h-12.0_r8*vec(3)*h2-20.0_r8*vec(4)*h3)*r  &
              +3.0_r8*vec(2)*r**2+4.0_r8*vec(3)*r**3+5.0_r8*vec(4)*r**4
         d2f(i)=-6.0_r8*vec(2)*h-12.0_r8*vec(3)*h2-20.0_r8*vec(4)*h3  &
               +6.0_r8*vec(2)*r+12.0_r8*vec(3)*r**2+20.0_r8*vec(4)*r**3
     endif
   enddo
   shift=f(notab)
   f=f-shift
   end subroutine


   function vpot(r)
   real(kind=r8) :: vpot,r
   vpot=v0/(1.0_r8+r**2)  !-v0/(1.0_r8+0.25_r8*el**2)
   return
   end function vpot

   subroutine getphi(x,u,is,du,d2u,v)
   real(kind=r8), dimension(ndim,npart) :: x,du
   real(kind=r8) :: v,u,d2u,uj,duj,d2uj
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
            uj=c1*utab(index-1)+c2*utab(index) &
              +c3*utab(index+1)+c4*utab(index+2)
            duj=c1*dutab(index-1)+c2*dutab(index) &
               +c3*dutab(index+1)+c4*dutab(index+2)
            d2uj=c1*d2utab(index-1)+c2*d2utab(index) &
                +c3*d2utab(index+1)+c4*d2utab(index+2)
            v=v+vpot(r)
            u=u+uj
            du(:,i)=du(:,i)+dx(:)*duj/(r*sqrt(mass(i)))
            du(:,j)=du(:,j)-dx(:)*duj/(r*sqrt(mass(j)))
            mu=mass(i)*mass(j)/(mass(i)+mass(j))
            d2u=d2u+(d2uj+(ndim-1.0_r8)*duj/r)/mu
         endif
      enddo
   enddo
   is=1
   d2u=d2u+sum(du*du)
   end subroutine getphi
end module potential
