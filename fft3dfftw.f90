!$Id: fft3dfftw.f90,v 1.1 2005/08/01 15:52:35 schmidt Exp $
module fft3d
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer, private, parameter :: i8=selected_int_kind(15)
   real(kind=r8), allocatable, save, private :: swork(:),cwork(:)
   real(kind=r8), private, save :: dr,dk,pi
   integer(kind=i4), private, save :: n
   integer(kind=i8), private, save :: plancos,plansin

contains
   subroutine fft3dinit(nin,rmax)
   include "fftw3.f"
   integer(kind=i4) :: nin
   real(kind=r8) :: rmax
   integer(kind=i4) :: iflags
!
! routine to setup fft
!
   iflags=FFTW_ESTIMATE
   iflags=FFTW_EXHAUSTIVE
   iflags=FFTW_MEASURE
   pi=4.0_r8*atan(1.0_r8)
   n=nin
   dr=rmax/n
   dk=pi/rmax
   if (allocated(cwork)) then
      deallocate(cwork,swork)
      call dfftw_destroy_plan(plancos)
      call dfftw_destroy_plan(plansin)
   endif
   allocate(cwork(0:n),swork(n-1))
   call dfftw_plan_r2r_1d(plancos,n+1,cwork,cwork,FFTW_REDFT00,iflags)
   call dfftw_plan_r2r_1d(plansin,n-1,swork,swork,FFTW_RODFT00,iflags)
   return
   end subroutine fft3dinit

   function rintegral(fr)
!
! Integrate fr(r) over the volume. This is the same as the k=0
! Fourier component.
!
   real (kind=r8) :: fr(0:n),r,rintegral
   integer(kind=i4) :: i
   rintegral=0.0_r8
   do i=1,n
      r=i*dr
      rintegral=rintegral+r**2*fr(i)
   enddo
   rintegral=4.0_r8*pi*dr*rintegral
   end function rintegral

   function xfmr2k0(fr)
!
! routine to do 3-d l=0 transform from r to k.
!
   real(kind=r8) :: fr(0:n),xfmr2k0(0:n)
   real(kind=r8) :: r,ak
   integer(kind=i4) :: i,l
   do i=1,n-1
      r=i*dr
      swork(i)=2.0_r8*pi*dr*r*fr(i)
   enddo
   xfmr2k0(0)=0.0_r8
   do i=1,n-1
      r=i*dr
      xfmr2k0(0)=xfmr2k0(0)+r*swork(i)
   enddo
   xfmr2k0(0)=2.0_r8*xfmr2k0(0)
   call dfftw_execute(plansin)
   do i=1,n-1
      ak=dk*i
      xfmr2k0(i)=swork(i)/ak
   enddo
   xfmr2k0(n)=0.0_r8
   return
   end  function xfmr2k0

   function xfmk2r0(fk)
!
! routine to do 3-d l=0 transform from k to r
!
   real(kind=r8) :: fk(0:n),xfmk2r0(0:n)
   real(kind=r8) :: ak,rfac
   integer(kind=i4) :: i
   do i=1,n-1
      ak=i*dk
      swork(i)=dk*ak*fk(i)
   enddo
   xfmk2r0(0)=0.0_r8
   do i=1,n-1
      ak=i*dk
      xfmk2r0(0)=xfmk2r0(0)+ak*swork(i)
   enddo
   xfmk2r0(0)=xfmk2r0(0)/(2.0_r8*pi**2)
   call dfftw_execute(plansin)
   do i=1,n-1
      rfac=dr*i*4.0_r8*pi**2
      xfmk2r0(i)=swork(i)/rfac
   enddo
   return
   end function xfmk2r0

   function xfmr2k1(fr)
!
! routine to do 3-d l=1 transform from r to k.
!
   real(kind=r8) :: fr(0:n),xfmr2k1(0:n)
   real(kind=r8) :: r,ak
   integer(kind=i4) :: i,l
   do i=0,n
      r=i*dr
      cwork(i)=2.0_r8*pi*dr*r*fr(i)
   enddo
   call dfftw_execute(plancos)
   xfmr2k1(0)=0.0_r8
   do i=1,n
      ak=dk*i
      xfmr2k1(i)=-cwork(i)/ak
   enddo
   do i=1,n-1
      r=i*dr
      swork(i)=2.0_r8*pi*dr*fr(i)
   enddo
   call dfftw_execute(plansin)
   do i=1,n-1
      ak=dk*i
      xfmr2k1(i)=xfmr2k1(i)+swork(i)/ak**2
   enddo
   end  function xfmr2k1

   function xfmk2r1(fk)
!
! routine to do 3-d l=1 transform from k to r
!
   real(kind=r8) :: fk(0:n),xfmk2r1(0:n)
   real(kind=r8) :: r,ak
   integer(kind=i4) :: i
   do i=0,n
      ak=i*dk
      cwork(i)=dk*ak*fk(i)/(4.0_r8*pi**2)
   enddo
   call dfftw_execute(plancos)
   xfmk2r1(0)=0.0_r8
   do i=1,n
      r=dr*i
      xfmk2r1(i)=-cwork(i)/r
   enddo
   do i=1,n-1
      ak=i*dk
      swork(i)=dk*fk(i)/(4.0_r8*pi**2)
   enddo
   call dfftw_execute(plansin)
   do i=1,n-1
      r=dr*i
      xfmk2r1(i)=xfmk2r1(i)+swork(i)/r**2
   enddo
   end function xfmk2r1

end module fft3d
