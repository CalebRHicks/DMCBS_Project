   program makepairs
   use random
   use random2
   implicit none
   integer, parameter :: i4=selected_int_kind(9)
   integer, parameter :: i8=selected_int_kind(15)
   integer, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4) :: nwalk,npart,nbin(2),i,j
   integer(kind=i8) :: irn,idx
   real(kind=r8) :: rho,el,dummy,wt,ran(1),dx(3)
   real(kind=r8), allocatable :: x(:,:),xup(:,:),xdn(:,:)
   open(unit=5,file='pairs.dat',status='unknown')
   read(5,*) npart
   read(5,*) rho
   read(5,*) nbin(1)
   read(5,*) nbin(2)
   if (nbin(1).lt.nbin(2)) then
      write(6,'(''This works only for nup.gt.ndown'')')
      stop
   endif
   read(5,*) nwalk
   allocate(x(3,npart),xup(3,nbin(1)),xdn(3,nbin(2)))
   el=(npart/rho)**(1.0_r8/3.0_r8)
   wt=1.0_r8
   rewind 9
   write(9,'(i10)') nwalk
   do i=1,nwalk
      xup=el*reshape(randn(size(xup)),shape(xup))
      xup=xup-el*nint(xup/el)
      xdn(:,1:nbin(2))=xup(:,1:nbin(2))+0.0001_r8*el*reshape(randn(3*nbin(2)),shape(xdn))
      x(:,1:nbin(1))=xup
      x(:,nbin(1)+1:npart)=xdn
      call ran2(dummy,irn)
      write (9,'(6e15.7)') x
      write (9,'(1e15.7)') wt
      write (9,'(i20)') irn
   enddo
   close(9)
   end program
