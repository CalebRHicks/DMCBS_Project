module vexternal
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9) 
   integer, private, save :: ndim,nbin(2),npart,idx
   real(kind=r8), private, save :: vext
contains
   subroutine setupvext(ndimin,nbinin,massin)
   integer(kind=i4) :: ndimin,nbinin(:)
   real(kind=r8) :: massin(:),dummy(2)
   ndim=ndimin
   nbin=nbinin
   npart=sum(nbin)
   vext=0.0_r8
   idx=0
   dummy=massin
   end subroutine setupvext

   subroutine addvext(x,wt)
   use wavefunction
   real(kind=r8) :: x(:,:),wt
   real(kind=r8) :: d1,d2(ndim,npart),d3,v,dummy
   call getonebody(x,d1,d2,d3,v)
   vext=vext+v
   idx=idx+1
   dummy=wt
   end subroutine addvext

   subroutine updatevext
   use mympi
   real(kind=r8) :: v
   integer(kind=i4) :: idxt
   call addall(vext,v)
   call addall(idx,idxt)
   if (myrank().eq.0) then
      v=v/idxt
      write (6,'(''vext ='',t30,1p,e14.5)') v
   endif
   vext=0.0_r8
   idx=0
   end subroutine
end module
