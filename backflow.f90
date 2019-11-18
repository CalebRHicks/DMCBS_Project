module backflow
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
contains
! this is a fake module used to compile any code
   subroutine setbf(ntabin,elin)
   use mympi
   integer(kind=i4) :: ntabin,dummy
   real(kind=r8) :: elin,dummy1
   dummy=ntabin
   dummy1=elin
   if (myrank().eq.0) then
      write(6,'(''No Backflow correlations'')')
   endif
   return
   end subroutine setbf

   subroutine getbf(r,bf,dbf,d2bf)
   real(kind=r8) :: r,bf,dbf,d2bf,dummy
   dummy=r
   bf=0.0_r8
   dbf=0.0_r8
   d2bf=0.0_r8
   return
   end subroutine getbf

   subroutine initbf(p)
   real(kind=r8) :: p(:),dummy
   dummy=p(1)
   end subroutine initbf

   subroutine getbfpar(p)
   real(kind=r8) :: p(:),dummy
   dummy=p(1)
   end subroutine getbfpar
end module backflow
