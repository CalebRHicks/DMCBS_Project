module ylmmod
implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)

   complex(kind=r8), private, parameter :: cone =(1.0_r8,0.0_r8)
   complex(kind=r8), private, parameter :: czero=(0.0_r8,0.0_r8)
   complex(kind=r8), private, parameter :: ci=(0.0_r8,1.0_r8)

contains
   function plmcal(x,lmax)
!
! routine to calculate the associated legendre functions P_lm
! from the recursion
! relation from l=0 to lmax. Note that plm should have dimension
! (lmax+1)*(lmax+2)/2 and l,m is indexed as (l*(l+1))/2+m.
! The conventions of J.D. Jackson, Classical Electrodynamics are used.
!
   integer(kind=i4) :: lmax,l,j,m
   real(kind=r8) :: x,s,fac
   real(kind=r8) :: plmcal(0:((lmax+1)*(lmax+2))/2-1)
!
! calculate P_ll and P_l,l-1
!
   s=sqrt(1.0_r8-x*x)
   plmcal(0)=1.0_r8
   fac=1.0_r8
   j=-1
   do l=1,lmax
      j=j+2
      fac=-j*s*fac
      plmcal((l*(l+1))/2+l)=fac
      plmcal((l*(l+1))/2+l-1)=j*x*plmcal(((l-1)*l)/2+l-1)
   enddo
!
! Use recurence to get others
!
   do m=0,lmax-2
      j=2*m+1
      do l=m+2,lmax
         j=j+2
         plmcal((l*(l+1))/2+m)=(x*j*plmcal((l*(l-1))/2+m) &
            -(l+m-1)*plmcal(((l-2)*(l-1))/2+m))/(l-m)
      enddo
   enddo
   return
   end function plmcal

   subroutine ylmcal(x,lmax,ylm,dylm)
!
! calculate r^l Y_lm using notation of J.D. Jackson, Classical
! Electrodynamics. These are indexed starting with 0 as lm=l*(l+1)+m
!
   real(kind=r8) :: x(0:2)
   integer(kind=i4) :: lmax
   complex(kind=r8) :: ylm(0:(lmax+1)*(lmax+1)-1)
   complex(kind=r8) :: dylm(0:2,0:(lmax+1)*(lmax+1)-1),ephi,ephim,dp,dm,d0
   real(kind=r8) :: plm(0:((lmax+1)*(lmax+2))/2-1)
   real(kind=r8) :: pi,pi4,x2,r,rpow,cphi,sphi,c
   real(kind=r8) :: fac,fac2,rxi,cp,cm,c0,xl
   integer(kind=i4) :: l,m,lm,lmm,isgn
   pi=4.0_r8*atan(1.0_r8)
   pi4=4.0_r8*pi
   x2=x(0)*x(0)+x(1)*x(1)
   r=sqrt(x2+x(2)*x(2))
   if (x2.eq.0.0_r8) then
      cphi=0.0_r8
      sphi=1.0_r8
      if (x(2).lt.0.0_r8) then
         c=-1.0_r8
      else
         c=1.0_r8
         endif
   else
      c=x(2)/r
      rxi=1.0_r8/sqrt(x2)
      cphi=x(0)*rxi
      sphi=x(1)*rxi
   endif
   ephi=cmplx(cphi,sphi,r8)
   plm=plmcal(c,lmax)
   do l=1,lmax
      ylm(l*(l+1):(l+1)*(l+1)-1)=plm((l*(l+1))/2:((l+1)*(l+2))/2-1)
   enddo
   ylm(0)=cmplx(1.0_r8/sqrt(pi4),0.0_r8,r8)
   rpow=1.0_r8
   do l=1,lmax
      rpow=r*rpow
      fac=rpow*sqrt((2*l+1)/pi4)
      ylm(l*(l+1))=fac*ylm(l*(l+1))
      ephim=cone
      do m=1,l
         fac2=(l+m)*(l+1-m)
         fac=fac/sqrt(fac2)
         ephim=ephi*ephim
         lm=l*(l+1)+m
         ylm(lm)=fac*real(ylm(lm),r8)*ephim
      enddo
   enddo
!
! now calculate gradient -- See Edmonds Angular Momentum in Quantum
! Mechanics for a derivation of the gradient formula
!
   dylm(:,0)=czero
   do l=1,lmax
      xl=(2*l+1)
      fac=sqrt(xl/(2*(2*l-1)))
      do m=0,l
         cp=fac*sqrt(0.5_r8*(l-m-1)*(l-m))
         cm=fac*sqrt(0.5_r8*(l+m-1)*(l+m))
         c0=fac*sqrt(2.0_r8*(l-m)*(l+m))
         if (m.gt.l-2) then
            dp=czero
         else
            dp=cp*ylm(l*(l-1)+m+1)
            endif
         if (m.eq.0) then
            if (l.gt.1) then
                dm=-cm*conjg(ylm(l*(l-1)+1))
             else
                dm=czero
             endif
         else
            dm=cm*ylm(l*(l-1)+m-1)
         endif
         if (m.gt.l-1) then
            d0=czero
         else
            d0=c0*ylm(l*(l-1)+m)
            endif
         lm=l*(l+1)+m
         dylm(0,lm)=dp-dm
         dylm(1,lm)=-ci*(dp+dm)
         dylm(2,lm)=d0
      enddo
   enddo
!
! The original code was for m >= 0. Fill in negative m here.
!
   do l=1,lmax
      isgn=1
      do m=1,l
         lm=l*(l+1)+m
         lmm=l*(l+1)-m
         isgn=-isgn
         ylm(lmm)=isgn*conjg(ylm(lm))
         dylm(:,lmm)=isgn*conjg(dylm(:,lm))
      enddo
   enddo
   end subroutine ylmcal
end module ylmmod
