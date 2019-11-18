module ewald
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, save :: ntab,nk
   real(kind=r8), private, save :: scale,dr,range
   real(kind=r8), private, save, allocatable :: phir(:),dphir(:),d2phir(:),phik(:)
   real(kind=r8), private, save, pointer :: ak(:,:)
contains
   subroutine setewald(nsh,ntabin,el2,phi,dphi,d2phi)
   use kshell
   real(kind=r8) :: pi,r,el2,phi(0:ntabin),dphi(0:ntabin),d2phi(0:ntabin)
   integer(kind=i4) :: nsh,ntabin,i
   real(kind=r8), pointer :: dummy2(:)
   real(kind=r8) :: sig,sigmax,sigmin,sigtmp
   integer(kind=i4) :: dummy1(nsh)
   real(kind=r8), parameter :: convergence=1.0e-8_r8
   real(kind=r8) :: cval,ak2,expo,f,df,d2f
   ntab=ntabin
   range=el2
   dr=range/ntab
   scale=1.0_r8/dr
   pi=4.0_r8*atan(1.0_r8)
   if (.not.allocated(phir)) then
      allocate(phir(0:ntab),dphir(0:ntab),d2phir(0:ntab))
      call setupk(2.0_r8*el2,nsh,dummy1,ak,dummy2,nk,3)
      allocate(phik(nk))
   endif
   phir=phi
   dphir=dphi
   d2phir=d2phi
   phik=0.0_r8
   if (dphir(ntab).gt.0.0_r8) then
      sigmin=-5000000.0_r8
      sigmax=0.0_r8
   else
      sigmin=0.0_r8
      sigmax=5000000.0_r8
   endif
   sig=0.5_r8*(sigmax+sigmin)
   sigtmp=sig
   do while (.true.)
      df=exp(-0.5_r8*el2**2/sig**2)*sqrt(2.0_r8)/(sqrt(pi)*sig*el2) &
         -erf(0.5_r8*el2*sqrt(2.0_r8)/sig)/el2**2
      if (abs(dphir(ntab)-df).lt.convergence) exit
      if (dphir(ntab).gt.df) then
         sigmin=sig
      else
         sigmax=sig
      endif
      sig=0.5_r8*(sigmax+sigmin)
      if (sig.eq.sigtmp) then
         write(6,'(''No convergence in Ewald setup'')')
         stop
         exit
      endif
      sigtmp=sig
   enddo
   cval=erf(el2/(sqrt(2.0_r8)*sig))/el2
   phir(0)=cval+phir(0)-sqrt(2.0_r8/pi)/sig
   dphir(0)=dphir(0)
   d2phir(0)=d2phir(0)+sqrt(2.0_r8/pi)/(3.0_r8*sig**3)
!  write (70,('(f10.5,3f18.10)')) 0.0,phir(0),dphir(0),d2phir(0)
   do i=1,ntab
      r=i*dr
      f=erf(r/(sqrt(2.0_r8)*sig))/r
      phir(i)=cval+phir(i)-f
      expo=exp(-0.5_r8*r**2/sig**2)*sqrt(2.0_r8)
      df=expo/(sqrt(pi)*sig*r)-f/r
      dphir(i)=dphir(i)-df
      d2f=-expo/(sqrt(pi)*sig**3)-2.0_r8*df/r
      d2phir(i)=d2phir(i)-d2f
!     write (70,('(f10.5,3f18.10)')) r,phir(i),dphir(i),d2phir(i)
   enddo
   do i=2,nk
      ak2=sum(ak(:,i)**2)
      phik(i)=4.0_r8*pi/(2.0_r8*el2)**3*exp(-0.5_r8*sig**2*ak2)/ak2
!     write(71,*) i,phik(i)
   enddo
   end subroutine setewald

   subroutine phiewald(x,phi,dphi,d2phi)
   real(kind=r8) :: x(3),phi,dphi(3),d2phi
   real(kind=r8) :: r,akr,dph,d2ph
   integer(kind=i4) :: idx,ik
   real(kind=r8) :: c1,c2,c3,c4
! Real space potential
   r=sqrt(sum(x**2))
   if (r.ge.range) then
      phi=phir(ntab)
      dph=0.0_r8
      d2ph=0.0_r8
   else
      dr=scale*r
      idx=dr
      idx=max(1,min(idx,ntab-2))
      dr=dr-idx
      c1=-dr*(dr-1.0_r8)*(dr-2.0_r8)/6.0_r8
      c2=(dr+1.0_r8)*(dr-1.0_r8)*(dr-2.0_r8)/2.0_r8
      c3=-(dr+1.0_r8)*dr*(dr-2.0_r8)/2.0_r8
      c4=(dr+1.0_r8)*dr*(dr-1.0_r8)/6.0_r8
      phi=c1*phir(idx-1)+c2*phir(idx)+c3*phir(idx+1)+c4*phir(idx+2)
      dph=c1*dphir(idx-1)+c2*dphir(idx)+c3*dphir(idx+1)+c4*dphir(idx+2)
      d2ph=c1*d2phir(idx-1)+c2*d2phir(idx)+c3*d2phir(idx+1)+c4*d2phir(idx+2)
   endif
   dphi=dph*x/r
   d2phi=d2ph+2.0_r8*dph/r
! K space
   do ik=1,nk
      akr=sum(ak(:,ik)*x(:))
      phi=phi+2.0_r8*phik(ik)*cos(akr)
      dphi=dphi-2.0_r8*ak(:,ik)*phik(ik)*sin(akr)
      d2phi=d2phi-2.0_r8*sum(ak(:,ik)**2)*phik(ik)*cos(akr)
   enddo
   end subroutine phiewald
end module ewald
