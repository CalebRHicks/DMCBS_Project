module matrixmod
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   complex(kind=r8), private, parameter :: cone = (1.0_r8,0.0_r8)
   complex(kind=r8), private, parameter :: czero = (0.0_r8,0.0_r8)

interface matinv
   module procedure rmatinv,cmatinv
end interface

contains
      subroutine rmatinv(a,detl,is,n)
!
! calculate inverse, log(abs(det)), and sign.
!
      integer(kind=i4) :: n,is,ipiv(n),i,info
      real(kind=r8) :: a(n,n),detl,work(n*n)
      real(kind=r8), parameter :: minushuge=-1e30_r8
      call dgetrf(n,n,a,n,ipiv,info)
      if (info.ne.0) then
!        write (6,'(1x,''error in dgetrf'',i10)') info
!        stop
         detl=minushuge
         return
      endif
      is=1
      detl=0.0_r8
      do i=1,n
         detl=detl+log(abs(a(i,i)))
         is=is*sign(1.0_r8,a(i,i))
         if (ipiv(i).ne.i) is=-is
      enddo
      call dgetri(n,a,n,ipiv,work,n*n,info)
      if (info.ne.0) then
!        write (6,'(1x,''error in dgetri'',i10)') info
!        stop
         detl=minushuge
      endif
      return
      end subroutine rmatinv

   subroutine cmatinv(a,detlog,n)
   integer(kind=i4) :: n,ipiv(n),info,i
   complex (kind=r8) :: a(n,n),detlog,cwork(n,n)
!
! lapack routine for lu factorization
!
   call zgetrf(n,n,a,n,ipiv,info)
   if (info.ne.0) then
      write (6,'(1x,''error in zgetrf'',i10)') info
      stop
   endif
!
! calculate determinant
!
   detlog=czero
   do i=1,n
      if (ipiv(i).ne.i) then
         detlog=detlog+log(-a(i,i))
      else
         detlog=detlog+log(a(i,i))
      endif
   enddo
!
! lapack routine to calculate inverse from factorization
!
   call zgetri(n,a,n,ipiv,cwork,n*n,info)
   if (info.ne.0) then
      write (6,'(1x,''error in zgetri'',i10)') info
      stop
   endif
   return
   end subroutine cmatinv
end module matrixmod
