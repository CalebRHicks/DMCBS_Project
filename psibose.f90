module bose
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
contains
   subroutine setupbose(rho)
   use mympi
   real(kind=r8) :: rho
   if (myrank().eq.0) then
      read (5,*) rho
      write (6,'(''rho ='',t40,f10.5)') rho ! number density
   endif
   call bcast(rho)
   end subroutine setupbose
end module bose
