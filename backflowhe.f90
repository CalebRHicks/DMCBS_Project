module backflow
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, save :: ntab
   real(kind=r8), private, save :: el,range,scale
   real(kind=r8), private, save :: lbf,sbf,rbf,wbf
   real(kind=r8) :: b0,db0,d2b0
contains
   subroutine setbf(ntabin,elin)
   use mympi
   integer(kind=i4) :: ntabin
   real(kind=r8) :: elin,pbf(4)
   ntab=ntabin
   el=elin
   if (myrank().eq.0) then
      read(5,*) lbf
      read(5,*) sbf
      read(5,*) rbf
      read(5,*) wbf
      write(6,'(''Helium-backflow parameters:'')')
      write(6,'(''lambda_B ='',t40,f10.5)') lbf
      write(6,'(''s_B ='',t40,f10.5)') sbf
      write(6,'(''r_B ='',t40,f10.5)') rbf
      write(6,'(''w_B ='',t40,f10.5)') wbf
   endif
   call bcast(lbf)
   call bcast(sbf)
   call bcast(rbf)
   call bcast(wbf)
   range=0.5_r8*el
   pbf(1)=lbf
   pbf(2)=sbf
   pbf(3)=rbf
   pbf(4)=wbf
   call initbf(pbf)
   end subroutine setbf

   subroutine initbf(pbf)
   real(kind=r8) :: pbf(:)
   lbf=pbf(1)
   sbf=pbf(2)
   rbf=pbf(3)
   wbf=pbf(4)
   call back(0.5_r8*el,b0,db0,d2b0)
   end subroutine initbf

   subroutine getbfpar(pbf)
   real(kind=r8) :: pbf(:)
   pbf(1)=lbf
   pbf(2)=sbf
   pbf(3)=rbf
   pbf(4)=wbf
   end subroutine getbfpar

   subroutine back(r,bf,dbf,d2bf)
   real(kind=r8) :: r,bf,dbf,d2bf
   real(kind=r8) :: ex,dr
   if (r.gt.rbf) then
      bf=0.0_r8
      dbf=0.0_r8
      d2bf=0.0_r8
      return
   endif
   dr=r-rbf
   ex=exp(-((r-sbf)/wbf)**2)
   bf=lbf*(dr/rbf)**3*ex
   dbf=lbf*(3.0_r8*dr**2*ex-2.0_r8*dr**3*(r-sbf)*ex/wbf**2)/rbf**3
   d2bf=lbf*(6.0_r8*dr*ex-12.0_r8*dr**2*(r-sbf)*ex/wbf**2-2.0_r8*dr**3*ex/wbf**2 &
             +4.0_r8*dr**3*(r-sbf)**2*ex/wbf**4)/rbf**3
   end subroutine back
 
   subroutine getbf(r,bf,dbf,d2bf)
   real(kind=r8) :: r,bf,dbf,d2bf
   real(kind=r8) :: b1,db1,d2b1
   if (r.gt.range) then
      bf=0.0_r8
      dbf=0.0_r8
      d2bf=0.0_r8
      return
   endif
   call back(r,bf,dbf,d2bf)
   call back(el-r,b1,db1,d2b1)
   bf=bf+b1-2.0_r8*b0
   dbf=dbf-db1
   d2bf=d2bf+d2b1
   end subroutine getbf
end module backflow
