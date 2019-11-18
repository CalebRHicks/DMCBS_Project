module betafun
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, save :: ntabb,ndim
   real(kind=r8), private, save :: el,rangeb,scaleb
   real(kind=r8), private, save :: b,c,d,small,b0,d0
   real(kind=r8), private, save :: c0,c1,c2,c3,c4
   real(kind=r8), private, save, allocatable :: swbtab(:),swdbtab(:),swd2btab(:)
   logical :: nobeta
contains
   subroutine setbeta(ndin,ntab,elin)
   use mympi
   real(kind=r8) :: bfact,bbeta,elin,dbeta,r
   integer(kind=i4) :: ntab,ndin,i
   el=elin
   ndim=ndin
   if (myrank().eq.0) then
      read (5,*) bbeta ! b parameter of s-wave beta-function
      read (5,*) dbeta ! d parameter of s-wave beta-function
      read (5,*) bfact ! multiply s-wave beta by this, if -1 read from file
      write (6,'(''b parameter of up-down beta-function ='',t40,f10.5)') bbeta
      write (6,'(''d parameter of up-down beta-function ='',t40,f10.5)') dbeta
      write (6,'(''multiply up-down beta-function by this factor ='',t40,f10.5)') bfact
   endif
   call bcast(bbeta)
   call bcast(dbeta)
   call bcast(bfact)
   if (bfact.eq.0.0_r8) then
      nobeta=.true.
      return
   endif
   nobeta=.false.
   if (bfact.ge.0.0_r8) then
! use the usual beta-function
      ntabb=ntab
      allocate(swbtab(0:ntabb),swdbtab(0:ntabb),swd2btab(0:ntabb))
      rangeb=0.5_r8*el
      scaleb=ntabb/rangeb
      call initbeta(bbeta,dbeta)
   else
! read the function from a file
      rewind 80
      read (80,*) ntabb
      ntabb=ntabb-1
      allocate(swbtab(0:ntabb),swdbtab(0:ntabb),swd2btab(0:ntabb))
      swbtab=0.0_r8
      swdbtab=0.0_r8
      swd2btab=0.0_r8
      do i=0,ntabb
         read (80,('(4e25.15)')) r,swbtab(i),swdbtab(i),swd2btab(i)
! make table for the laplacian
         if (r.ne.0.0_r8) swd2btab(i)=swd2btab(i)+(ndim-1.0_r8)*swdbtab(i)/r
      enddo
      rangeb=r
      if (rangeb.gt.0.5_r8*el) then
         write (6,'(''rangebeta > el2, stopping... '',2e15.7)') rangeb,0.5_r8*el
         stop
      endif
      scaleb=real(ntabb)/rangeb
   endif
   swbtab=abs(bfact)*swbtab
   swdbtab=abs(bfact)*swdbtab
   swd2btab=abs(bfact)*swd2btab
!  do i=0,ntabb
!     if ((myrank().eq.0).and.(mod(i,10000).eq.0)) then
!        r=i/scaleb
!!        write (30,'(f10.5,3e15.7)') 2.0_r8*r/el,swbtab(i),swdbtab(i),swd2btab(i)
!        write (30,'(f10.5,3e15.7)') r,swbtab(i),swdbtab(i),swd2btab(i)
!     endif
!  enddo
   end subroutine setbeta

   subroutine initbeta(bbeta,dbeta)
   real(kind=r8) :: dum1,dum2,betan,b1,b1p,b1pp,b2,b2p,b2pp,b3,b3p,b3pp,dr,r
   integer(kind=i4) :: nj,i,j
   real(kind=r8) :: bbeta,dbeta
   b=bbeta
   d=dbeta
   c=0.5_r8*(2.0_r8+2.0_r8*d*b*el+(d*b*el)**2*exp(b*el*(1.0_r8+d)) &
    +2.0_r8*d*b**2*el**2*exp(b*el*(1.0_r8+d))+2.0_r8*b*el &
    -2.0_r8*exp(d*b*el)*b*el-2.0_r8*exp(d*b*el)) &
    /(el**2*(exp(d*b*el)-1.0_r8-d+d*exp(b*el*(1.0_r8+d)))*b**2)
   c0=d
! use an expansion at small distances
   small=1.0e-3_r8/b
   c1=-d+c*d-d**2/2.0_r8
   c2=d/2.0_r8+d**2/2.0_r8+d**3/6.0_r8+c*(-d-d**2/2.0_r8)
   c3=-d/6.0_r8-d**2/4.0_r8-d**3/6.0_r8-d**4/24.0_r8 &
      +c*(d/2.0_r8+d**2/2.0_r8+d**3/6.0_r8)
   c4=d/24.0_r8+d**2/12.0_r8+d**3/12.0_r8+d**4/24.0_r8+d**5/120.0_r8 &
      -c*(d/6.0_r8+d**2/4.0_r8+d**3/6.0_r8+d**4/24.0_r8)
   dr=0.5_r8*el/ntabb
   nj=0
   betan=10.0_r8
   do while (abs(betan).gt.1.0d-10)
      nj=nj+1
      call beta(el*nj,betan,dum1,dum2)
   enddo
   do i=0,ntabb
      r=i*dr
      call beta(r,b1,b1p,b1pp)
      swbtab(i)=b1
      swdbtab(i)=b1p
      swd2btab(i)=b1pp
      do j=1,nj
         call beta(j*el-r,b2,b2p,b2pp)
         call beta(j*el+r,b3,b3p,b3pp)
         swbtab(i)=swbtab(i)+b2+b3
         swdbtab(i)=swdbtab(i)-b2p+b3p
         swd2btab(i)=swd2btab(i)+b2pp+b3pp
      enddo
   enddo
   do i=0,ntabb
      r=i*dr
      swbtab(i)=swbtab(i)-swbtab(ntabb)
      if (i.ne.0) then
         swd2btab(i)=swd2btab(i)+(ndim-1.0_r8)*swdbtab(i)/r
      else
         swd2btab(i)=3.0_r8*swd2btab(i) ! assume quadratic at origin
      endif
   enddo
   b0=b
   d0=d
   end subroutine initbeta

   subroutine beta0(r,b1,b1p,b1pp)
   real(kind=r8) :: r,b1,b1p,b1pp,ex,ex10,x
   if (r.lt.small) then
      x=b*r
      b1=c0+x*(c1+x*(c2+x*(c3+x*c4)))
      b1p=c1+x*(2.0_r8*c2+x*(3.0_r8*c3+4.0_r8*c4*x))
      b1pp=2.0_r8*c2+x*(6.0_r8*c3+12.0_r8*c4*x)
   else
      x=b*r
      ex=exp(-x)
      ex10=exp(-d*x)
      b1=(1.0_r8+c*x)*(1.0_r8-ex10)*ex/x
      b1p=ex*(-1.0_r8-x*(1.0_r8+c*x) &
         +ex10*(1.0_r8+x*(1.0_r8+d)*(1.0_r8+c*x)))/x**2
      b1pp=ex*(2.0_r8+x*(2.0_r8+x*(1.0_r8+c*x)) &
         -ex10*(2.0_r8+x*(1.0_r8+d)*(2.0_r8+x*(1.0_r8+d)*(1.0_r8+x*c))))/x**3
   endif
   b1p=b*b1p
   b1pp=b**2*b1pp
   b1=b1/d
   b1p=b1p/d
   b1pp=b1pp/d
   end subroutine beta0

   subroutine beta(r,b1,b1p,b1pp)
   real(kind=r8) :: r,b1,b1p,b1pp
   real(kind=r8) :: bet1,bet1p,bet1pp,bet2,bet2p,bet2pp,bet0
   real(kind=r8) :: tiny=1.0e-14_r8
   if (r.gt.el*0.5_r8) then
      b1=0.0_r8
      b1p=0.0_r8
      b1pp=0.0_r8
   else
      call beta0(el*0.5_r8,bet0,b1p,b1pp)
      call beta0(r,bet1,bet1p,bet1pp)
      call beta0(el-r,bet2,bet2p,bet2pp)
      b1=bet1+bet2-2.0_r8*bet0
      if (r.lt.tiny) then
         b1p=0.0_r8
      else
         b1p=bet1p-bet2p
      endif
      b1pp=bet1pp+bet2pp
   endif
   end subroutine beta

   subroutine radialfun(r,f,df,d2f)
   real(kind=r8) :: r,f,df,d2f
   real(kind=r8) :: ddr,a1,a2,a3,a4
   integer(kind=i4) :: index
   if (nobeta.or.r.gt.rangeb) then
      f=0.0_r8
      df=0.0_r8
      d2f=0.0_r8
      return
   endif
   ddr=scaleb*r
   index=ddr
   index=max(1,min(index,ntabb-2))
   ddr=ddr-index
   a1=-ddr*(ddr-1.0_r8)*(ddr-2.0_r8)/6.0_r8
   a2=(ddr+1.0_r8)*(ddr-1.0_r8)*(ddr-2.0_r8)/2.0_r8
   a3=-(ddr+1.0_r8)*ddr*(ddr-2.0_r8)/2.0_r8
   a4=(ddr+1.0_r8)*ddr*(ddr-1.0_r8)/6.0_r8
   f=a1*swbtab(index-1)+a2*swbtab(index)+a3*swbtab(index+1)+a4*swbtab(index+2)
   df=a1*swdbtab(index-1)+a2*swdbtab(index)+a3*swdbtab(index+1)+a4*swdbtab(index+2)
   d2f=a1*swd2btab(index-1)+a2*swd2btab(index)+a3*swd2btab(index+1)+a4*swd2btab(index+2)
   end subroutine radialfun

   subroutine derbeta(r,dbet)
   real(kind=r8) :: r,dbet(:)
   real(kind=r8) :: t1,t2,t3,t4,t5,t6,t10,t12,t13,t14,t15,t16,t17,t18,t21,t22,t25,t29,t30,t31,t32,t33,t34,t36
   real(kind=r8) :: t37,t38,t39,t40,t41,t42,t43,t44,t46,t47,t48,t49,t50,t54,t55,t56,t57,t58,t59,t60,t62,t63,t65,t66
   real(kind=r8) :: t67,t71,t72,t75,t77,t78,t79,t82,t83,t84,t88,t90,t91,t93,t94,t96,t97,t101,t106,t107,t109
   real(kind=r8) :: t110,t112,t117,t118,t119,t121,t122,t129
   real(kind=r8) :: t130,t133,t140,t141,t144,t151,t156
   dbet=0.0_r8
   if (nobeta.or.r.gt.rangeb) return
   t1 = el ** 2
   t3 = d0 * b0
   t4 = t3 * el
   t5 = exp(t4)
   t6 = t5 * b0
   t12 = d0 ** 2
   t14 = b0 * el
   t15 = d0 + 1
   t17 = exp(t14 * t15)
   t18 = t1 * t17
   t21 = b0 ** 2
   t22 = t12 * t21
   t25 = t1 * el * t15 * t17
   t29 = d0 * t21
   t32 = d0 * el
   t33 = t32 * t5
   t36 = 2 * d0 * t1 * t6 + 2 * t5 * el - 2 * el - 2 * t12 * b0 * t18 &
   - t22 * t25 - 4 * t3 * t18 - 2 * t29 * t25 + 2 * t33 - 2 * t32
   t37 = 0.1D1 / t1
   t38 = t36 * t37
   t40 = d0 * t17 - 1 - d0 + t5
   t41 = 0.1D1 / t40
   t42 = 0.1D1 / b0
   t43 = t41 * t42
   t44 = t43 * r
   t54 = 2 * t6 * el - 2 - 2 * t14 - t22 * t18 - 2 * t29 * t18 + 2 * t5 - 2 * t4
   t55 = t54 * t37
   t56 = t40 ** 2
   t57 = 0.1D1 / t56
   t58 = t55 * t57
   t62 = t32 * t15 * t17 + t33
   t65 = 0.1D1 / t21
   t66 = t41 * t65
   t71 = exp(-t3 * r)
   t72 = 1 - t71
   t75 = exp(-b0 * r)
   t77 = 1 / d0
   t78 = t77 * t42
   t79 = 0.1D1 / r
   t84 = 1 - t55 * t44 / 2
   t88 = t84 * t72
   t93 = t77 * t65
   t96 = el - r
   t97 = t43 * t96
   t106 = exp(-t3 * t96)
   t107 = 1 - t106
   t110 = exp(-b0 * t96)
   t112 = 0.1D1 / t96
   t117 = 1 - t55 * t97 / 2
   t121 = t117 * t107
   t130 = 0.1D1 / el
   t133 = t54 * t130
   t140 = exp(-t4 / 2)
   t141 = 1 - t140
   t144 = exp(-t14 / 2)
   t151 = 1 - t133 * t43 / 4
   t156 = t151 * t141
   dbet(1) = (-t38 * t44 + t58 * t42 * r * t62 + t55 * t66 * r) * t72 *   &
     t75 * t78 * t79 / 2 + t84 * t71 * t75 * t42 - t88 * t75 * t77 * t42 &
     - t88 * t75 * t93 * t79 + (-t38 * t97 + t58 * t42 * t96 * t62 +  &
     t55 * t66 * t96) * t107 * t110 * t78 * t112 / 2 + t117 * t106 * t110 &
     * t42 - t121 * t96 * t110 * t77 * t42 * t112 - t121 * t110 * t93 * &
     t112 - (-t36 * t130 * t43 + t133 * t57 * t42 * t62 + t133 * t66) * &
     t141 * t144 * t78 * t130 - 2 * t151 * t140 * t144 * t42 + 2 *t156 * &
     t144 * t77 * t42 + 4 * t156 * t144 * t93 * t130
   t1 = b0 ** 2
   t2 = el ** 2
   t3 = t1 * t2
   t4 = d0 * b0
   t5 = t4 * el
   t6 = exp(t5)
   t10 = b0 * el
   t13 = exp(t10 * (d0 + 1))
   t14 = t2 * t13
   t16 = 2 * d0 * t1 * t14
   t17 = d0 ** 2
   t18 = t1 * b0
   t21 = t2 * el * t13
   t29 = t6 * b0 * el
   t30 = 2 * t29
   t31 = 2 * t10
   t32 = 2 * t3 * t6 - t16 - t17 * t18 * t21 - 2 * t3 * t13 - 2 * d0 &
   * t18 * t21 + t30 - t31
   t33 = 0.1D1 / t2
   t34 = t32 * t33
   t36 = d0 * t13 - 1 - d0 + t6
   t38 = 0.1D1 / b0
   t39 = 0.1D1 / t36 * t38
   t40 = t39 * r
   t46 = t30 - 2 - t31 - t1 * t17 * t14 - t16 + 2 * t6 - 2 * t5
   t47 = t46 * t33
   t48 = t36 ** 2
   t49 = 0.1D1 / t48
   t50 = t47 * t49
   t54 = t13 + t4 * el * t13 - 1 + t29
   t59 = exp(-t4 * r)
   t60 = 1 - t59
   t63 = exp(-b0 * r)
   t65 = 1 / d0
   t66 = t65 * t38
   t67 = 0.1D1 / r
   t72 = 1 - t47 * t40 / 2
   t79 = 1 / t17 * t38
   t82 = el - r
   t83 = t39 * t82
   t90 = exp(-t4 * t82)
   t91 = 1 - t90
   t94 = exp(-b0 * t82)
   t96 = 0.1D1 / t82
   t101 = 1 - t47 * t83 / 2
   t109 = 0.1D1 / el
   t112 = t46 * t109
   t118 = exp(-t5 / 2)
   t119 = 1 - t118
   t122 = exp(-t10 / 2)
   t129 = 1 - t112 * t39 / 4
   dbet(2) = (-t34 * t40 + t50 * t38 * r * t54) * t60 * t63 * t66 * t67 &
   / 2 + t72 * t59 * t63 * t65 - t72 * t60 * t63 * t79 * t67 + (-t34 &
   * t83 + t50 * t38 * t82 * t54) * t91 * t94 * t66 * t96 / 2 + t101 &
   * t90 * t94 * t65 - t101 * t91 * t94 * t79 * t96 - (-t32 * t109 * &
   t39 + t112 * t49 * t38 * t54) * t119 * t122 * t66 * t109 - 2 * t129 &
   * t118 * t122 * t65 + 4 * t129 * t119 * t122 * t79 * t109
   end subroutine derbeta

   subroutine getbetadata(bbeta,dbeta)
   real(kind=r8) :: bbeta,dbeta
   bbeta=b0
   dbeta=d0
   end subroutine getbetadata
end module betafun
