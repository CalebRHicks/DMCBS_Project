module backflow
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, save :: ntab
   real(kind=r8), private, save :: el,range,scale
   real(kind=r8), private, save, allocatable :: bftab(:),dbftab(:),d2bftab(:)
   real(kind=r8), private, save :: lbf,sbf,rbf,wbf
contains
   subroutine setbf(ntabin,elin)
   use mympi
   integer(kind=i4) :: ntabin,i
   real(kind=r8) :: elin,dr,r
   real(kind=r8) :: b0,db0,d2b0,b1,db1,d2b1
   ntab=ntabin
   el=elin
   if (myrank().eq.0) then
      read(5,*) lbf
      read(5,*) sbf
      read(5,*) rbf
      read(5,*) wbf
      write(6,'(''Backflow parameters:'')')
      write(6,'(''lambda_B ='',t40,f10.5)') lbf
      write(6,'(''s_B ='',t40,f10.5)') sbf
      write(6,'(''r_B ='',t40,f10.5)') rbf
      write(6,'(''w_B ='',t40,f10.5)') wbf
   endif
   call bcast(lbf)
   call bcast(sbf)
   call bcast(rbf)
   call bcast(wbf)
   allocate(bftab(0:ntab),dbftab(0:ntab),d2bftab(0:ntab))
   range=0.5_r8*el
   dr=range/ntab
   scale=1.0_r8/dr
   call back(0.5_r8*el,b0,db0,d2b0)
   do i=0,ntab
      r=i*dr
      call back(r,bftab(i),dbftab(i),d2bftab(i))
      call back(el-r,b1,db1,d2b1)
      bftab(i)=bftab(i)+b1-2.0_r8*b0
      dbftab(i)=dbftab(i)-db1
      d2bftab(i)=d2bftab(i)+d2b1
!     if (myrank().eq.0) write(79,'(f10.5,3e15.7)') r,bftab(i),dbftab(i),d2bftab(i)
   enddo
   end subroutine setbf

   subroutine back(r,bf,dbf,d2bf)
   real(kind=r8) :: r,bf,dbf,d2bf
   real(kind=r8) :: r32,r52,r72
   r32=r**(3.0_r8/2.0_r8)
   r52=r32*r
   r72=r52*r
   bf=lbf*(1.0_r8+sbf*r)/(rbf+wbf*r+r72)
   dbf=lbf*(2.0_r8*sbf*rbf-5.0_r8*sbf*r72-2.0_r8*wbf-7.0_r8*r52) &
      /(2.0_r8*(rbf+wbf*r+r72)**2)
   d2bf=lbf*(-8.0_r8*sbf*wbf*rbf-15.0_r8*sbf*wbf*r72-63.0_r8*sbf*r52*rbf &
       +35.0_r8*sbf*r**6+8.0_r8*wbf**2+21.0_r8*wbf*r52+63.0_r8*r**5-35.0_r8*r32*rbf) &
       /(4.0_r8*(rbf+wbf*r+r72)**3)
   end subroutine back
 
   subroutine getbf(r,bf,dbf,d2bf)
   real(kind=r8) :: r,bf,dbf,d2bf
   real(kind=r8) :: dr,c1,c2,c3,c4
   integer(kind=i4) :: index
   if (r.gt.range) then
      bf=0.0_r8
      dbf=0.0_r8
      d2bf=0.0_r8
      return
   endif
   dr=scale*r
   index=dr
   index=max(1,min(index,ntab-2))
   dr=dr-index
   c1=-dr*(dr-1.0_r8)*(dr-2.0_r8)/6.0_r8
   c2=(dr+1.0_r8)*(dr-1.0_r8)*(dr-2.0_r8)/2.0_r8
   c3=-(dr+1.0_r8)*dr*(dr-2.0_r8)/2.0_r8
   c4=(dr+1.0_r8)*dr*(dr-1.0_r8)/6.0_r8
   bf=c1*bftab(index-1)+c2*bftab(index)+c3*bftab(index+1)+c4*bftab(index+2)
   dbf=c1*dbftab(index-1)+c2*dbftab(index)+c3*dbftab(index+1)+c4*dbftab(index+2)
   d2bf=c1*d2bftab(index-1)+c2*d2bftab(index)+c3*d2bftab(index+1)+c4*d2bftab(index+2)
   end subroutine getbf
end module backflow
