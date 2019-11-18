module radii
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9) 
   integer, private, save :: ndim,nbin(2),npart,idx
   real(kind=r8), private, save :: radius,radup,raddn
   real(kind=r8), private, save, allocatable :: mass(:)
contains
   subroutine setupradii(ndimin,nbinin,massin)
   integer(kind=i4) :: ndimin,nbinin(:)
   real(kind=r8) :: massin(:)
   ndim=ndimin
   nbin=nbinin
   npart=sum(nbin)
   allocate(mass(npart))
   mass(1:nbin(1))=massin(1)
   if (nbin(1).ne.npart) mass(nbin(1)+1:npart)=massin(2)
   radius=0.0_r8
   radup=0.0_r8
   raddn=0.0_r8
   idx=0
   end subroutine setupradii

   subroutine addradii(x,wt)
   real(kind=r8) :: x(:,:),wt,xcm(3)
   real(kind=r8) :: r2,massi
   integer(kind=i4) :: i
   massi=1.0_r8/sum(mass(:))
   xcm(1)=(sum(x(1,:)*mass(:)))*massi
   xcm(2)=(sum(x(2,:)*mass(:)))*massi
   xcm(3)=(sum(x(3,:)*mass(:)))*massi
   do i=1,nbin(1)
      r2=sum((x(:,i)-xcm(:))**2)
      radup=radup+r2*wt
      radius=radius+r2*wt
   enddo
   do i=nbin(1)+1,nbin(1)+nbin(2)
      r2=sum((x(:,i)-xcm(:))**2)
      raddn=raddn+r2*wt
      radius=radius+r2*wt
   enddo
   idx=idx+1
   end subroutine addradii

   subroutine updateradii
   use mympi
   real(kind=r8) :: r,rup,rdn
   integer(kind=i4) :: idxt
   call addall(radius,r)
   call addall(radup,rup)
   call addall(raddn,rdn)
   call addall(idx,idxt)
   if (myrank().eq.0) then
      r=r/npart/idxt
      rup=rup/nbin(1)/idxt
      rdn=rdn/nbin(2)/idxt
      write (6,'(''rms radius ='',t30,1p,e14.5)') r
      write (6,'(''rmsup radius ='',t30,1p,e14.5)') rup
      write (6,'(''rmsdn radius ='',t30,1p,e14.5)') rdn
   endif
   radius=0.0_r8
   radup=0.0_r8
   raddn=0.0_r8
   idx=0
   end subroutine
end module
