module phijoe
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   real(kind=r8), private, parameter :: pi=4.0_r8*atan(1.0_r8)
   real(kind=r8), private, save :: el,kbes,reff,cutoff
   integer(kind=i4), private, save :: npart,nsh,ibox
   real(kind=r8), private, save, allocatable :: ck(:)
contains
   subroutine setupphijoe(elin,npartin)
   use mympi
   integer(kind=i4) :: npartin
   real(kind=r8) :: elin
   real(kind=r8) :: kf,xi,c,amu,dk,xkx,xky,xkz,q
   integer(kind=i4) :: k,ikx,iky,ikz
   real(kind=r8) :: daw1,daw2
   el=elin
   npart=npartin
   if (myrank().eq.0) then
      read(5,*) xi
      read(5,*) amu
      read(5,*) nsh
      read(5,*) ibox
      read(5,*) c 
      cutoff=c/el
      write(6,'(''xi parameter in coshfun ='',t40,f10.5)') xi
      write(6,'(''amu in coshfun ='',t40,f10.5)') amu
      write(6,'(''number of shell in the long-range ='',t40,i10)') nsh
      write(6,'(''number of boxes in the short-range ='',t40,i10)') ibox
      write(6,'(''cutoff of short-range part of coshfun ='',t40,f10.5)') cutoff
   endif
   call bcast(xi)
   call bcast(amu)
   call bcast(nsh)
   call bcast(ibox)
   call bcast(cutoff)
   kf=(3.0_r8*pi**2*npart)**(1.0_r8/3.0_r8)/el
   reff=2.0_r8/amu
   kbes=kf*sqrt(3.0_r8*xi/20.0_r8)
   if (myrank().eq.0) write(6,'(''k**2 ='',t40,f10.5)') kbes**2
   allocate(ck((2*nsh+1)**3))
   dk=(2.0_r8*pi)/el
   k=1
   do ikx=-nsh,nsh
      xkx=ikx*dk
      do iky=-nsh,nsh
         xky=iky*dk
         do ikz=-nsh,nsh
            xkz=ikz*dk
            q=max(1.d-10,sqrt(xkx**2+xky**2+xkz**2))
            call dawsonf((kbes-q)/(2.0_r8*cutoff),daw1)
            call dawsonf((kbes+q)/(2.0_r8*cutoff),daw2)
            ck(k)=(2.0_r8*pi*(-(1.0_r8/(exp((kbes-q)**2/(4.0_r8*cutoff**2)) &
              *(kbes-q)))+1.0_r8/(exp((kbes+q)**2/(4.0_r8*cutoff**2))*(kbes+q))))/q &
              +(2.0_r8*kbes*sqrt(pi)*reff*(daw1/(kbes-q) &
              -daw2/(kbes + q)))/q
            k=k+1
         enddo
      enddo
   enddo
   end subroutine setupphijoe

   subroutine getphilr(x,phi,dphi,d2phi)   ! calculate the long-range part of phi
   real(kind=r8) :: x(:),phi,dphi(:),d2phi
   real(kind=r8) :: c,s
   real(kind=r8) :: a,dk,xkx,xkxx,xky,xkyy,xkz,xkzz,r
   integer(kind=i4) :: ikx,iky,ikz,k
   a=1.0_r8/el**3
   r=sqrt(sum(x**2))
   dk=(2.0_r8*pi)/el
   phi=0.0_r8
   dphi=0.0_r8
   d2phi=0.0_r8
   k=1
   do ikx=-nsh,nsh
      xkx=ikx*dk
      xkxx=xkx*x(1)
      do iky=-nsh,nsh
         xky=iky*dk
         xkyy=xky*x(2)
         do ikz=-nsh,nsh
            xkz=ikz*dk
            xkzz=xkz*x(3)
            c=cos(xkxx+xkyy+xkzz)
            s=sin(xkxx+xkyy+xkzz)
            phi=phi+ck(k)*a*c
            dphi(1)=dphi(1)-ck(k)*a*s*xkx
            dphi(2)=dphi(2)-ck(k)*a*s*xky
            dphi(3)=dphi(3)-ck(k)*a*s*xkz
            d2phi=d2phi-ck(k)*a*c*(xkx**2+xky**2+xkz**2)
            k=k+1
         enddo
      enddo
   enddo
   return
   end subroutine getphilr

   subroutine getphisr(x,phi,dphi,d2phi)
   real(kind=r8) :: x(:)
   real(kind=r8) :: phi,dphi(:),d2phi,ph,dph,d2ph
   real(kind=r8) :: ph1,dph1,d2ph1,ph2,dph2,d2ph2,ex,dex,d2ex
   real(kind=r8) :: r,s,c,th,x0(3)
   integer(kind=i4) :: inx,iny,inz
   phi=0.0_r8
   dphi=0.0_r8
   d2phi=0.0_r8
   do inx=-ibox,ibox
      do iny=-ibox,ibox
         do inz=-ibox,ibox
            x0(1)=x(1)+inx*el
            x0(2)=x(2)+iny*el
            x0(3)=x(3)+inz*el
            r=sqrt(sum(x0**2))
            s=sin(kbes*r)
            c=cos(kbes*r)
            th=tanh(2.0_r8*r/reff)
            ph1=kbes*reff/2.0_r8*s+c
            dph1=kbes*(kbes*reff/2.0_r8*c-s)
            d2ph1=-kbes**2*ph1
            ph2=c*(th-1.0_r8)
            dph2=-kbes*s*(th-1.0_r8)+2.0_r8/reff*c*(1.0_r8-th**2)
            d2ph2=-kbes**2*c*(th-1.0_r8)-kbes*4.0_r8/reff*s*(1.0_r8-th**2)-8.0_r8/reff**2*c*th*(1.0_r8-th**2)
            ex=erfc(cutoff*r)
            dex=-2.0_r8*cutoff*exp(-(cutoff*r)**2)/sqrt(pi)
            d2ex=-2.0_r8*cutoff**2*r*dex
            ph=ex*ph1+ph2
            dph=dph1*ex+ph1*dex+dph2
            d2ph=d2ph1*ex+2.0_r8*dph1*dex+ph1*d2ex+d2ph2
            d2ph=d2ph/r-2.0_r8*dph/r**2+2.0_r8*ph/r**3
            dph=dph/r-ph/r**2
            ph=ph/r
            phi=phi+ph
            dphi(:)=dphi(:)+dph*x0(:)/r
            d2phi=d2phi+d2ph+2.0_r8/r*dph
         enddo
      enddo
   enddo
   return
   end subroutine getphisr

   subroutine getphijoe(x,phi,dphi,d2phi)
   real (kind=r8) :: x(:),phi,dphi(:),d2phi
   real(kind=r8) :: philr,dphilr(3),d2philr
   real(kind=r8) :: phisr,dphisr(3),d2phisr
   call getphilr (x,philr,dphilr,d2philr)
   call getphisr (x,phisr,dphisr,d2phisr)
   phi=philr+phisr
   dphi=dphilr+dphisr
   d2phi=d2philr+d2phisr
   return
   end subroutine getphijoe


   subroutine dawsonf(xx,daw)
   integer :: i
   real(kind=r8) :: &
     &     DAW,FRAC,HALF,ONE,ONE225,P1,P2,P3,P4,Q1,Q2,Q3,Q4,SIX25, &
     &     SUMP,SUMQ,TWO5,W2,X,XX,Y,XLARGE,XMAX,XSMALL,ZERO
      DIMENSION P1(10),P2(10),P3(10),P4(10),Q1(10),Q2(9),Q3(9),Q4(9)
      DATA ZERO,HALF,ONE/0.0D0,0.5D0,1.0D0/, &
     &     SIX25,ONE225,TWO5/6.25D0,12.25D0,25.0D0/
      DATA XSMALL/1.05D-08/, XLARGE/9.49D+07/, XMAX/2.24D+307/
      DATA P1/-2.69020398788704782410D-12, 4.18572065374337710778D-10, &
     &        -1.34848304455939419963D-08, 9.28264872583444852976D-07, &
     &        -1.23877783329049120592D-05, 4.07205792429155826266D-04, &
     &        -2.84388121441008500446D-03, 4.70139022887204722217D-02, &
     &        -1.38868086253931995101D-01, 1.00000000000000000004D+00/
      DATA Q1/ 1.71257170854690554214D-10, 1.19266846372297253797D-08, &
     &         4.32287827678631772231D-07, 1.03867633767414421898D-05, &
     &         1.78910965284246249340D-04, 2.26061077235076703171D-03, &
     &         2.07422774641447644725D-02, 1.32212955897210128811D-01, &
     &         5.27798580412734677256D-01, 1.00000000000000000000D+00/
      DATA P2/-1.70953804700855494930D+00,-3.79258977271042880786D+01, &
     &         2.61935631268825992835D+01, 1.25808703738951251885D+01, &
     &        -2.27571829525075891337D+01, 4.56604250725163310122D+00, &
     &        -7.33080089896402870750D+00, 4.65842087940015295573D+01, &
     &        -1.73717177843672791149D+01, 5.00260183622027967838D-01/
      DATA Q2/ 1.82180093313514478378D+00, 1.10067081034515532891D+03, &
     &        -7.08465686676573000364D+00, 4.53642111102577727153D+02, &
     &         4.06209742218935689922D+01, 3.02890110610122663923D+02, &
     &         1.70641269745236227356D+02, 9.51190923960381458747D+02, &
     &         2.06522691539642105009D-01/
      DATA P3/-4.55169503255094815112D+00,-1.86647123338493852582D+01, &
     &        -7.36315669126830526754D+00,-6.68407240337696756838D+01, &
     &         4.84507265081491452130D+01, 2.69790586735467649969D+01, &
     &        -3.35044149820592449072D+01, 7.50964459838919612289D+00, &
     &        -1.48432341823343965307D+00, 4.99999810924858824981D-01/
      DATA Q3/ 4.47820908025971749852D+01, 9.98607198039452081913D+01, &
     &         1.40238373126149385228D+01, 3.48817758822286353588D+03, &
     &        -9.18871385293215873406D+00, 1.24018500009917163023D+03, &
     &        -6.88024952504512254535D+01,-2.31251575385145143070D+00, &
     &         2.50041492369922381761D-01/
      DATA P4/-8.11753647558432685797D+00,-3.84043882477454453430D+01, &
     &        -2.23787669028751886675D+01,-2.88301992467056105854D+01, &
     &        -5.99085540418222002197D+00,-1.13867365736066102577D+01, &
     &        -6.52828727526980741590D+00,-4.50002293000355585708D+00, &
     &        -2.50000000088955834952D+00, 5.00000000000000488400D-01/
      DATA Q4/ 2.69382300417238816428D+02, 5.04198958742465752861D+01, &
     &         6.11539671480115846173D+01, 2.08210246935564547889D+02, &
     &         1.97325365692316183531D+01,-1.22097010558934838708D+01, &
     &        -6.99732735041547247161D+00,-2.49999970104184464568D+00, &
     &         7.49999999999027092188D-01/
      X = XX
      IF (ABS(X) .GT. XLARGE) THEN
            IF (ABS(X) .LE. XMAX) THEN
                  DAW = HALF / X
               ELSE
                  DAW = ZERO
            END IF
         ELSE IF (ABS(X) .LT. XSMALL) THEN
            DAW = X
         ELSE
            Y = X * X
            IF (Y .LT. SIX25) THEN
                  SUMP = P1(1)
                  SUMQ = Q1(1)
                  DO 100 I = 2, 10
                     SUMP = SUMP * Y + P1(I)
                     SUMQ = SUMQ * Y + Q1(I)
  100             CONTINUE
                  DAW = X * SUMP / SUMQ
               ELSE IF (Y .LT. ONE225) THEN
                  FRAC = ZERO
                  DO 200 I = 1, 9
  200                FRAC = Q2(I) / (P2(I) + Y + FRAC)
                  DAW = (P2(10) + FRAC) / X
               ELSE IF (Y .LT. TWO5) THEN
                  FRAC = ZERO
                  DO 300 I = 1, 9
  300                FRAC = Q3(I) / (P3(I) + Y + FRAC)
                  DAW = (P3(10) + FRAC) / X
               ELSE
                  W2 = ONE / X / X
                  FRAC = ZERO
                  DO 400 I = 1, 9
  400                FRAC = Q4(I) / (P4(I) + Y + FRAC)
                  FRAC = P4(10) + FRAC
                  DAW = (HALF + HALF * W2 * FRAC) / X
            END IF
      endif
      return
      end subroutine



end module phijoe
