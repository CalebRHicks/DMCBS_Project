   program testloff
   implicit none
   integer, parameter :: i4=selected_int_kind(9)
   integer, parameter :: r8=selected_real_kind(15,9)
   real(kind=r8), pointer :: akup(:,:),akdn(:,:),ak2up(:),ak2dn(:)
   integer(kind=i4), allocatable :: numshup(:),numshdn(:)
   integer(kind=i4) :: norb,nkup,nkdn
   integer(kind=i4) :: i,j,k,nqfin
   real(kind=r8) :: pi,el,c
   real(kind=r8), allocatable :: q(:,:,:),qfin(:,:,:),qfin2(:,:,:)
   real(kind=r8) :: tmp(3,3)
   pi=4.0_r8*atan(1.0_r8)
   el=2.0_r8*pi
   norb=10
   allocate(numshup(norb),numshdn(norb))
   call setupk(el,norb,numshup,akup,ak2up,nkup,3)
   call setupk(el,norb,numshdn,akdn,ak2dn,nkdn,3)
   allocate(q(3,3,2*nkup*nkdn),qfin(3,3,2*nkup*nkdn))
   allocate(qfin2(3,3,2*nkup*nkdn))
   q=0
   qfin=0
   qfin2=0
   k=0
   do i=1,nkup
      do j=1,nkdn
         k=k+1
         q(1,:,k)=akup(:,i)+akdn(:,j)
         q(2,:,k)=akup(:,i)
         q(3,:,k)=akdn(:,j)
!write(6,'(6f10.0,''  a'')') q(2,:,k),q(3,:,k)
         if ((sum(akup(:,i)**2).ne.0).and.(sum(akdn(:,j)**2).ne.0)) then
            k=k+1
            q(1,:,k)=akup(:,i)-akdn(:,j)
            q(2,:,k)=akup(:,i)
            q(3,:,k)=-akdn(:,j)
!write(6,'(6f10.0,''  b'')') q(2,:,k),q(3,:,k)
!print*,sum(akup(:,i)**2),sum(akdn(:,i)**2)
         endif
      enddo
   enddo
   j=0
   do i=1,k
      if (q(1,2,i).ge.0.and.q(1,3,i).ge.0) then
!     if (q(1,1,i).ge.0.and.q(1,2,i).ge.0.and.q(1,3,i).ge.0) then
!        write(6,'(i5,3f10.0)') i,q(1,:,i)
         j=j+1
         qfin(:,:,j)=q(:,:,i)
      endif
   enddo
   nqfin=j
   do i=1,nqfin-1
      do j=i+1,nqfin
         if (sum(qfin(1,:,i)**2).gt.sum(qfin(1,:,j)**2)) then
            tmp=qfin(:,:,i)
            qfin(:,:,i)=qfin(:,:,j)
            qfin(:,:,j)=tmp
         endif
      enddo
   enddo
   j=0
   do i=1,nqfin
      if (sum((qfin(2,:,i)+qfin(3,:,i))**2).eq.0.0_r8) then
         c=0.01_r8
      else
         c=0.00001_r8
      endif
!     write(6,'(f10.5,6f10.0)') c,qfin(2,:,i),qfin(3,:,i)
      write(6,'(f10.5,9f10.0,f10.5)') c,qfin(1,:,i),qfin(2,:,i),qfin(3,:,i),sum(qfin(1,:,i)**2)
      if (sum(qfin(1,:,i)**2).eq.1) then
         j=j+1
         qfin2(:,:,j)=qfin(:,:,i)
      endif
   enddo
stop
   nqfin=j
   do i=1,nqfin-1
      do j=i+1,nqfin
         if (qfin2(2,1,i).gt.qfin2(2,1,j)) then
            tmp=qfin2(:,:,i)
            qfin2(:,:,i)=qfin2(:,:,j)
            qfin2(:,:,j)=tmp
         else if ((qfin2(2,1,i).eq.qfin2(2,1,j)).and.(qfin2(2,2,i).gt.qfin2(2,2,j))) then
            tmp=qfin2(:,:,i)
            qfin2(:,:,i)=qfin2(:,:,j)
            qfin2(:,:,j)=tmp
         endif
      enddo
   enddo
   do i=1,nqfin
      write(6,'(i5,9f10.0,f10.5)') i,qfin2(1,:,i),qfin2(2,:,i),qfin2(3,:,i),sum(qfin2(1,:,i)**2)
   enddo
       
         
         
contains
   subroutine setupk(el,norb,numsh,ak,ak2,nk,ndim)
   real(kind=r8) :: el
   integer(kind=i4) :: norb,numsh(:),nk,ndim
   real(kind=r8), pointer :: ak(:,:),ak2(:)
   if (ndim.eq.3) then
      call kstates3d(el,norb,numsh,ak,ak2,nk)
   elseif (ndim.eq.2) then
      call kstates2d(el,norb,numsh,ak,ak2,nk)
   elseif (ndim.eq.1) then
      write(6,'(''setupk not working yet for 1D case'')')
      stop
   endif
   end subroutine

   subroutine kstates3d(elbox,nsh0,numsh0,ak,ak2,nk)
   integer(kind=i4) :: i2max,imax,ix,iy,iz,ix2,iy2,iz2,i2,i2m
   integer(kind=i4) :: nk,i,j,tmp(3),nsh,k
   integer(kind=i4), allocatable :: iktmp(:,:),numsh(:)
   integer(kind=i4) :: nsh0,numsh0(:)
   real(kind=r8) :: elbox,pi,piel
   real(kind=r8), pointer :: ak(:,:),ak2(:)
   i2max=2*nsh0
   imax=nint(sqrt(real(i2max)))
   nk=0
   do ix=0,imax
      ix2=ix*ix
      do iy=0,imax
         iy2=iy*iy
         do iz=0,imax
            iz2=iz*iz
            i2=ix2+iy2+iz2
            if (i2.le.i2max) then
               nk=nk+1
               if (ix.eq.0) then
                  if (iy.ne.0.and.iz.ne.0) nk=nk+1
               endif
               if (iy.eq.0) then
                  if (ix.ne.0.and.iz.ne.0) nk=nk+1
               endif
               if (iz.eq.0) then
                  if (ix.ne.0.and.iy.ne.0) nk=nk+1
               endif
               if (ix.ne.0.and.iy.ne.0.and.iz.ne.0) nk=nk+3
            endif
         enddo
      enddo
   enddo
   allocate(iktmp(nk,3))
   iktmp=0
   nk=0 
   do ix=0,imax
      ix2=ix*ix
      do iy=0,imax
         iy2=iy*iy
         do iz=0,imax
            iz2=iz*iz
            i2=ix2+iy2+iz2
            if (i2.le.i2max) then
               nk=nk+1
               iktmp(nk,1)=ix
               iktmp(nk,2)=iy
               iktmp(nk,3)=iz
               if (ix.eq.0) then
                  if (iy.ne.0.and.iz.ne.0) then
                     nk=nk+1
                     iktmp(nk,1)=ix
                     iktmp(nk,2)=iy
                     iktmp(nk,3)=-iz
                  endif
               endif
               if (iy.eq.0) then
                  if (ix.ne.0.and.iz.ne.0) then
                     nk=nk+1
                     iktmp(nk,1)=ix
                     iktmp(nk,2)=iy
                     iktmp(nk,3)=-iz
                  endif
               endif
               if (iz.eq.0) then
                  if (ix.ne.0.and.iy.ne.0) then
                     nk=nk+1
                     iktmp(nk,1)=ix
                     iktmp(nk,2)=-iy
                     iktmp(nk,3)=iz
                  endif
               endif
               if (ix.ne.0.and.iy.ne.0.and.iz.ne.0) then
                  nk=nk+1
                  iktmp(nk,1)=ix
                  iktmp(nk,2)=iy
                  iktmp(nk,3)=-iz
                  nk=nk+1
                  iktmp(nk,1)=ix
                  iktmp(nk,2)=-iy
                  iktmp(nk,3)=iz
                  nk=nk+1
                  iktmp(nk,1)=ix
                  iktmp(nk,2)=-iy
                  iktmp(nk,3)=-iz
               endif
             endif
         enddo
      enddo
   enddo
   do i=1,nk-1
      do j=i+1,nk
         if (sum(iktmp(i,:)**2).ge.sum(iktmp(j,:)**2)) then
            tmp=iktmp(i,:)
            iktmp(i,:)=iktmp(j,:)
            iktmp(j,:)=tmp
         endif
      enddo
   enddo
   allocate(numsh(nk)) ! this is too large, but it works
   numsh=0
   i2m=0
   nsh=1
   do i=1,nk
      i2=sum(iktmp(i,:)**2)
      if (i2.gt.i2m) then
         nsh=nsh+1
         i2m=i2
      endif
      numsh(nsh)=numsh(nsh)+1
   enddo        
   do i=1,nsh
      numsh(i)=2*numsh(i)
   enddo
   pi=4.0_r8*atan(1.0_r8)
   piel=2.0_r8*pi/elbox
   nk=1
   if (nsh0.ge.nsh) then
      write (6,'(''Too many states, stopping...'')')
      stop
   endif
   do i=2,nsh0
      nk=nk+numsh(i)/2
   enddo
   allocate(ak(3,nk),ak2(nk))
   k=0
   do i=1,nsh0
      numsh0(i)=numsh(i)/2
      do j=1,max(1,numsh(i)/2)
         k=k+1
         ak(1,k)=piel*iktmp(k,1)
         ak(2,k)=piel*iktmp(k,2)
         ak(3,k)=piel*iktmp(k,3)
         ak2(k)=sum(ak(:,k)**2)
      enddo
   enddo
   end subroutine kstates3d

   subroutine kstates2d(elbox,nsh0,numsh0,ak,ak2,nk)
   integer(kind=i4) :: i2max,imax,ix,iy,ix2,iy2,i2,i2m
   integer(kind=i4) :: nk,i,j,tmp(2),nsh,k
   integer(kind=i4), allocatable :: iktmp(:,:),numsh(:)
   integer(kind=i4) :: nsh0,numsh0(:)
   real(kind=r8) :: elbox,pi,piel
   real(kind=r8), pointer :: ak(:,:),ak2(:)
   i2max=3*nsh0
   imax=nint(sqrt(real(i2max)))
   nk=0
   do ix=0,imax
      ix2=ix*ix
      do iy=0,imax
         iy2=iy*iy
         i2=ix2+iy2
         if (i2.le.i2max) then
            nk=nk+1
            if (ix.eq.0) then
               if (iy.ne.0) nk=nk+1
            endif
            if (iy.eq.0) then
               if (ix.ne.0) nk=nk+1
            endif
            if (ix.ne.0.and.iy.ne.0) nk=nk+2
         endif
      enddo
   enddo
   allocate(iktmp(nk,2))
   iktmp=0
   nk=0 
   do ix=0,imax
      ix2=ix*ix
      do iy=0,imax
         iy2=iy*iy
         i2=ix2+iy2
         if (i2.le.i2max) then
            nk=nk+1
            iktmp(nk,1)=ix
            iktmp(nk,2)=iy
            if (ix.ne.0.and.iy.ne.0) then
               nk=nk+1
               iktmp(nk,1)=ix
               iktmp(nk,2)=-iy
             endif
         endif
      enddo
   enddo
   do i=1,nk-1
      do j=i+1,nk
         if (sum(iktmp(i,:)**2).ge.sum(iktmp(j,:)**2)) then
            tmp=iktmp(i,:)
            iktmp(i,:)=iktmp(j,:)
            iktmp(j,:)=tmp
         endif
      enddo
   enddo
   allocate(numsh(nk)) ! this is too large, but it works
   numsh=0
   i2m=0
   nsh=1
   do i=1,nk
      i2=sum(iktmp(i,:)**2)
      if (i2.gt.i2m) then
         nsh=nsh+1
         i2m=i2
      endif
      numsh(nsh)=numsh(nsh)+1
   enddo        
   do i=1,nsh
      numsh(i)=2*numsh(i)
   enddo
   pi=4.0_r8*atan(1.0_r8)
   piel=2.0_r8*pi/elbox
   nk=1
   if (nsh0.ge.nsh) then
      write (6,'(''Too many states, stopping...'')')
      stop
   endif
   do i=2,nsh0
      nk=nk+numsh(i)/2
   enddo
   allocate(ak(2,nk),ak2(nk))
   k=0
   do i=1,nsh0
      numsh0(i)=numsh(i)/2
      do j=1,max(1,numsh(i)/2)
         k=k+1
         ak(1,k)=piel*iktmp(k,1)
         ak(2,k)=piel*iktmp(k,2)
         ak2(k)=sum(ak(:,k)**2)
      enddo
   enddo
   end subroutine kstates2d
   end program
