module pfaffianmod
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)

contains
   subroutine pfafnopivot(a,n,p)
!
! no pivoting, destroys a, pfaffian returned in p
!
   integer(kind=i4) :: n,i,j
   real(kind=r8) :: a(n,n),p,fac
   if (mod(n,2).ne.0) then
      p=0.0_r8
      return
   endif
   p=1.0_r8
   do i=1,n,2
      do j=i+2,n
         fac=-a(i,j)/a(i,i+1)
         a(i+1:n,j)=a(i+1:n,j)+fac*a(i+1:n,i+1)
         a(j,i+1:n)=a(j,i+1:n)+fac*a(i+1,i+1:n)
      enddo
      p=p*a(i,i+1)
   enddo
   return
   end subroutine pfafnopivot

   subroutine pfaf(a,n,p)
!
! stupid pivoting, destroys a, pfaffian returned in p
!
   integer(kind=i4) :: n,i,j,ipiv
   real(kind=r8) :: a(n,n),p,fac,temp
   real(kind=r8) :: amag
   if (mod(n,2).ne.0) then
      p=0.0_r8
      return
   endif
   p=1.0_r8
   do i=1,n,2
      amag=abs(a(i,i+1))
      ipiv=i+1
      do j=i+2,n
        if (abs(a(i,j)).gt.amag) then
           ipiv=j
           amag=abs(a(i,j))
        endif
      end do
      if (ipiv.ne.i+1) then
         p=-p
         do j=i,n
            temp=a(j,ipiv)
            a(j,ipiv)=a(j,i+1)
            a(j,i+1)=temp
         enddo
         do j=i,n
            temp=a(ipiv,j)
            a(ipiv,j)=a(i+1,j)
            a(i+1,j)=temp
         enddo
      endif
      do j=i+2,n
         fac=-a(i,j)/a(i,i+1)
         a(i+1:n,j)=a(i+1:n,j)+fac*a(i+1:n,i+1)
         a(j,i+1:n)=a(j,i+1:n)+fac*a(i+1,i+1:n)
      enddo
      p=p*a(i,i+1)
   enddo
   return
   end subroutine pfaf

   subroutine updateinv(a,row,i,n)
   real(kind=r8) :: a(n,n),row(n),prod(n)
   integer(kind=i4) :: i,j,n
!
! update row and column using skew-symmetry of matrix
! 
    prod=matmul(a,row)
    a(i,:)=-a(i,:)/prod(i)
    a(:,i)=-a(:,i)/prod(i)
    do j=1,n
       if (j.ne.i) then
          a(j,:)=a(j,:)+prod(j)*a(i,:)
          a(:,j)=a(:,j)+prod(j)*a(:,i)
       endif
    enddo
   return
   end subroutine updateinv

   subroutine updateinv0(a,row,i,n)
   real(kind=r8) :: a(n,n),row(n),prod(n)
   integer(kind=i4) :: i,j,n
!
! update column (note row is -column)
! 
    prod=matmul(a,row)
    a(i,:)=-a(i,:)/prod(i)
    do j=1,n
       if (j.ne.i) then
          a(j,:)=a(j,:)+prod(j)*a(i,:)
       endif
    enddo
!
! update row
!
    prod=matmul(transpose(a),row)
    a(:,i)=a(:,i)/prod(i)
    do j=1,n
       if (j.ne.i) then
          a(:,j)=a(:,j)-prod(j)*a(:,i)
       endif
    enddo
   return
   end subroutine updateinv0

end module pfaffianmod
