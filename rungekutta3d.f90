module rungekutta
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
contains
   subroutine orbital(psi,dpsi,d2psi,hb,ntab,dr,vpot,nodes)
   integer(kind=i4) :: nodes,i,ntab,ncross,iter
   real(kind=r8), dimension(0:ntab) :: psi,dpsi,d2psi
   real(kind=r8) :: r,ene,hb,der,emax,emin,dr
   real(kind=r8), parameter :: convergence=1.0e-8_r8
   integer(kind=i4), parameter :: itmax=100
   external :: vpot
   real(kind=r8) :: vpot
   emax=vpot(0.0_r8)
   emin=emax
   do i=1,ntab
      emax=max(emax,vpot(i*dr))
      emin=min(emin,vpot(i*dr))
   enddo 
   emax=0.001_r8
   call rk(psi,dpsi,emax,hb,dr,ntab,ncross,vpot)
   do while (ncross.le.nodes)
      emax=2.0_r8*emax
      call rk(psi,dpsi,emax,hb,dr,ntab,ncross,vpot)
   enddo
   ene=0.5_r8*(emin+emax)
   iter=0
   do while(.true.)
      call rk(psi,dpsi,ene,hb,dr,ntab,ncross,vpot)
      r=dr*ntab
      der=(dpsi(ntab)-psi(ntab)/r)/r
      der=der/psi(ntab)
      if (der.lt.0.0_r8) then
         emax=ene
      else
         emin=ene
      endif
      if (abs(der).lt.convergence) exit
      ene=0.5_r8*(emin+emax)
      iter=iter+1
      if (iter.gt.itmax) then
         write (6,'(''orbital: no convergence'')') 
         stop
      endif
   enddo
   psi(0)=dpsi(0)
   d2psi(0)=(vpot(0.0_r8)-ene)*dpsi(0)/(3.0_r8*hb)
   dpsi(0)=0.0_r8
   do i=1,ntab
      r=i*dr
      dpsi(i)=(dpsi(i)-psi(i)/r)/r
      psi(i)=psi(i)/r
      d2psi(i)=(vpot(r)-ene)*psi(i)/hb-2.0_r8*dpsi(i)/r
   enddo
   end subroutine orbital

   subroutine rk(psi,dpsi,ene,hb,dr,ntab,ncross,vpot)
   integer(kind=i4) :: ntab,ncross,i
   real(kind=r8) :: psi(0:ntab),dpsi(0:ntab),ene,hb,dr,r
   real(kind=r8) :: k1a,k2a,k3a,k4a,k1b,k2b,k3b,k4b
   external :: vpot
   real(kind=r8) :: vpot
   psi(0)=0.0_r8
   dpsi(0)=1.0_r8
   ncross=0
   do i=1,ntab
      r=(i-1)*dr
      k1a=-(ene-vpot(r))*psi(i-1)/hb
      k1b=dpsi(i-1)
      k2a=-(ene-vpot(r+0.5_r8*dr))*(psi(i-1)+0.5_r8*dr*k1b)/hb
      k2b=dpsi(i-1)+0.5_r8*dr*k1a
      k3a=-(ene-vpot(r+0.5_r8*dr))*(psi(i-1)+0.5_r8*dr*k2b)/hb
      k3b=dpsi(i-1)+0.5_r8*dr*k2a
      k4a=-(ene-vpot(r+dr))*(psi(i-1)+dr*k3b)/hb
      k4b=dpsi(i-1)+dr*k3a
      dpsi(i)=dpsi(i-1)+(k1a+2.0_r8*k2a+2.0_r8*k3a+k4a)*dr/6.0_r8
      psi(i)=psi(i-1)+(k1b+2.0_r8*k2b+2.0_r8*k3b+k4b)*dr/6.0_r8
      if (psi(i).ne.sign(psi(i),psi(i-1))) ncross=ncross+1
   enddo
   end subroutine rk
end module rungekutta
