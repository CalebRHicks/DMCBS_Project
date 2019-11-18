module estimator
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   character(len=30), private, dimension(:), allocatable, save :: label
   real(kind=r8), private, dimension(:), allocatable, save :: &
      valtot,valblk,valnow,avbad,av2bad,wttot,wtblk
   real(kind=r8), private, dimension(:), allocatable, save :: valblk2,errnow
   integer(kind=i4), private, save :: nblock,nest
!
! label = label for estimator
! valtot = sum of the valblk values for blocks
! valblk = current sum of value for block
! valnow = current average for block
! avbad = sum of (valblk/wtblk)
! av2bad = sum of (valblk/wtblk)^2
! wttot = sum of wtblk
! wtblk = current sum of weights for block
! nblock = number of blocks averaged
! nest = number of estimators
!

contains
   subroutine setestnum(n)
!
! set up arrays for n estimators
!
   integer(kind=i4) :: n
   nest=n
   allocate(valtot(n),valblk(n),valnow(n),avbad(n),av2bad(n))
   allocate(valblk2(n),errnow(n))
   allocate(wttot(n),wtblk(n))
   allocate(label(n))
   return
   end subroutine setestnum

   subroutine zerest
   integer(kind=i4) :: i
!
! zero all estimators
!
   do i=1,nest
      valtot(i)=0.0_r8
      valblk(i)=0.0_r8
      valblk2(i)=0.0_r8
      valnow(i)=0.0_r8
      errnow(i)=0.0_r8
      avbad(i)=0.0_r8
      av2bad(i)=0.0_r8
      wtblk(i)=0.0_r8
      wttot(i)=0.0_r8
   enddo
   nblock=0
   return
   end subroutine zerest

   subroutine addest(i,l)
!
! set label l for estimator i
!
   integer(kind=i4) :: i
   character(len=*) :: l
   integer(kind=i4) :: ln,j
   ln=min(len(l),len(label(i)))
   label(i)(1:ln)=l(1:ln)
   do j=ln+1,len(label(i))
      label(i)(j:j)=" "
   enddo
   return
   end subroutine addest

   subroutine addval(i,val,wt)
!
! add another value, val, with weight, wt to estimator i
!
   integer(kind=i4) :: i
   real(kind=r8) :: val,wt
   valblk(i)=valblk(i)+val*wt
   valblk2(i)=valblk2(i)+val*val*wt
   wtblk(i)=wtblk(i)+wt
   return
   end subroutine addval

   subroutine update
   use mympi
!
! block all estimators
!
   integer(kind=i4) :: i
   real(kind=r8) :: valsum(nest),wtsum(nest),valsum2(nest)
   call addall(valblk,valsum)
   call addall(valblk2,valsum2)
   call addall(wtblk,wtsum)
   if (myrank().eq.0) then
      valblk=valsum
      valblk2=valsum2
      wtblk=wtsum
      do i=1,nest
         valnow(i)=valblk(i)/wtblk(i)
         valblk2(i)=valblk2(i)/wtblk(i)
         errnow(i)=sqrt(abs(valblk2(i)-valnow(i)**2)/wtblk(i))
         avbad(i)=avbad(i)+valnow(i)
         av2bad(i)=av2bad(i)+valnow(i)**2
         wttot(i)=wttot(i)+wtblk(i)
         valtot(i)=valtot(i)+valblk(i)
      enddo
   endif
   valblk=0.0_r8
   valblk2=0.0_r8
   wtblk=0.0_r8
   nblock=nblock+1
   return
   end subroutine update
 
   subroutine result(i,vnow,enow,val,err)
!
! return current result, value for current block in vnow, and
! current average and error in val and err
!
   integer(kind=i4) :: i
   real(kind=r8) :: vnow,val,err,enow
!  real(kind=r8) :: av,av2
   vnow=valnow(i)
   enow=errnow(i)
!  val=valtot(i)/wttot(i)
!  av=avbad(i)/nblock
!  av2=av2bad(i)/nblock
!  err=sqrt(abs(av2-av**2)/max(1,nblock-1))
! prints block values, not total averages
   val=vnow
   err=enow
   return
   end subroutine result

   function resstring(tau)
!
! return all results in a string for printing
!
   character(len=90), dimension(nest) :: resstring
   integer(kind=i4) :: i
   real(kind=r8) :: valn,errn,val,err
   real(kind=r8) :: tau
   do i=1,nest
      call result(i,valn,errn,val,err)
      write (resstring(i),'(a30,t30,1x,1p,e15.8,e19.10,'' +- '',e19.10)') &
         label(i),tau,val,err
   enddo
   return
   end function resstring

end module estimator
