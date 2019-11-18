module stack
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: i8=selected_int_kind(15)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   type :: walker
      real(kind=r8) :: psil,d2psi,v,weight,psigl,d2psig,psil0,psigl0,vext
      real(kind=r8), pointer :: x(:,:),dpsi(:,:),dpsig(:,:)
      integer(kind=i4) :: is,isg
      integer(kind=i8) :: irn
      real(kind=r8), pointer :: x0(:,:)
      real(kind=i4) :: wt0
   end type
   type (walker), private, allocatable, save :: s(:,:)
   integer(kind=i4), private, allocatable, save :: ist(:) 
   integer(kind=i4), private, save :: nstack,nwalk
   integer(kind=i4), private, save :: npart,ndim

   type (walker), public, save :: w1,w2

   interface assignment (=) ! overload equal operator for deep copy
      module procedure copywalker
   end interface

contains
   subroutine copywalker(wl,wr)
!
! overloaded = operator for deep copy of walker
!
   type (walker), intent(out) :: wl
   type (walker), intent(in) :: wr
   wl%psil=wr%psil
   wl%d2psi=wr%d2psi
   wl%v=wr%v
   wl%weight=wr%weight
   wl%psigl=wr%psigl
   wl%d2psig=wr%d2psig
   wl%psil0=wr%psil0
   wl%psigl0=wr%psigl0
   wl%vext=wr%vext
   wl%x(:,:)=wr%x(:,:)
   wl%dpsi(:,:)=wr%dpsi(:,:)
   wl%dpsig(:,:)=wr%dpsig(:,:)
   wl%is=wr%is
   wl%isg=wr%isg
   wl%irn=wr%irn
   wl%x0(:,:)=wr%x0(:,:)
   wl%wt0=wr%wt0
   end subroutine copywalker

   subroutine create(nst,nw,npartin,ndimin,euclidean)
   integer(kind=i4) :: nst,nw,npartin,ndimin
   integer(kind=i4) :: i,j
   integer :: err
   logical :: euclidean
   nwalk=nw
   nstack=nst
   npart=npartin
   ndim=ndimin
   allocate(s(nwalk,nstack),ist(nstack),stat=err)
   do i=1,nstack
      ist(i)=0
      do j=1,nwalk
         allocate(s(j,i)%x(ndim,npart),s(j,i)%dpsi(ndim,npart),s(j,i)%dpsig(ndim,npart),stat=err)
         if (euclidean) then
            allocate(s(j,i)%x0(ndim,npart),stat=err)
         else
            allocate(s(j,i)%x0(1,1)) ! save memory and communication if not using euclidean stuff
         endif
      enddo
   enddo
   allocate(w1%x(ndim,npart),w1%dpsi(ndim,npart),w1%dpsig(ndim,npart),stat=err)
   allocate(w2%x(ndim,npart),w2%dpsi(ndim,npart),w2%dpsig(ndim,npart),stat=err)
   if (euclidean) then
      allocate(w1%x0(ndim,npart),w2%x0(ndim,npart))
   else
      allocate(w1%x0(1,1),w2%x0(1,1))
   endif
   return
   end subroutine create

   subroutine push(i,w)
!
! push a walker onto stack i
!
   integer(kind=i4) :: i
   type(walker) :: w
   if (ist(i).ge.nwalk) then
      write (6,'(1x,'' stack overflow'',2i10)') ist(i),nwalk
      stop
      endif
   ist(i)=ist(i)+1
   s(ist(i),i)=w
   return
   end subroutine push

   subroutine pop(i,w,empty)
!
! pop a walker off stack i. If stack is empty, variabel empty is true.
!
   integer(kind=i4) :: i
   type(walker) :: w
   logical :: empty
   empty=ist(i).eq.0
   if (.not.empty) then
      w=s(ist(i),i)
      ist(i)=ist(i)-1
      endif
   return
   end subroutine pop

   function numstack(i)
!
! return number of walkers left on stack i
!
   integer(kind=i4) :: i,numstack
   numstack=ist(i)
   return
   end function numstack
end module stack
