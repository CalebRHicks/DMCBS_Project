module potential
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, save, allocatable :: nbin(:),spin(:)
   integer(kind=i4), private, save :: npart,ndim,ntabu,ntabv
   real(kind=r8), private, save :: el,eli
   real(kind=r8), private, save :: rangev,scalev,rangeu,scaleu
   real(kind=r8), private, save, allocatable :: utab(:),dutab(:),d2utab(:)
   real(kind=r8), private, save, allocatable :: vtab(:)
   real(kind=r8), private, save, allocatable :: mass(:)
contains
   subroutine setphi(ndin,ntabin,nbinin,elin,massin,hbin)
   use mympi
   real(kind=r8) :: elin,massin(2),hbin
   real(kind=r8) :: r,dru,drv,dummy,cjas
   integer(kind=i4) :: ntabin,nbinin(:),i,ndin
   dummy=ntabin
   dummy=hbin
   ndim=ndin
   allocate(nbin(2))
   nbin=nbinin
   el=elin
   if (el.ne.0.0_r8) then
      eli=1.0_r8/el
   else
      eli=0.0_r8
   endif
   npart=sum(nbin)
   if (myrank().eq.0) then
      read (5,*) cjas !jastrow strength
      write (6,'(''potential from unit 90'')')
      write (6,'(''Jastrow from unit 91'')')
      write (6,'(''Jastrow strength ='',t40,f10.5)') cjas
   endif
   call bcast(cjas)
   allocate(spin(npart),mass(npart))
! assign spin 1 to up particles and -1 to down ones
   spin(1:nbin(1))=1
   spin(nbin(1)+1:nbin(1)+nbin(2))=-1
   mass(1:nbin(1))=massin(1)
   mass(nbin(1)+1:npart)=massin(2)
! read potential table from unit 90
! first value is the number of points
   read (90,*) ntabv
   ntabv=ntabv-1
   allocate(vtab(0:ntabv))
   read (90,*) r,vtab(0)
   read (90,*) drv,vtab(1)
   drv=drv-r
   scalev=1.0_r8/drv
   do i=2,ntabv
      read (90,*) r,vtab(i)
   enddo
   rangev=r
! read u and derivatives from unit 91
   read (91,*) ntabu
   ntabu=ntabu-1
   allocate(utab(0:ntabu),dutab(0:ntabu),d2utab(0:ntabu))
   read (91,*) r,utab(0),dutab(0),d2utab(0)
   read (91,*) dru,utab(1),dutab(1),d2utab(1)
   dru=dru-r
   scaleu=1.0_r8/dru
   do i=2,ntabu
      read (91,*) r,utab(i),dutab(i),d2utab(i)
   enddo
   utab=-cjas*utab
   dutab=-cjas*dutab
   d2utab=-cjas*d2utab
   rangeu=r
   end subroutine setphi

   subroutine getphi(x,u,is,du,d2u,v)
   real(kind=r8), dimension(ndim,npart) :: x,du
   real(kind=r8) :: v,u,d2u,uj,duj,d2uj
   integer(kind=i4) :: is,i,j,index
   real(kind=r8), dimension(ndim) :: dx
   real(kind=r8) :: r,dr,c1,c2,c3,c4,mu
   real(kind=r8), parameter :: tiny=1.0e-14_r8
   u=0.0_r8
   du=0.0_r8
   d2u=0.0_r8
   v=0.0_r8
   do i=1,npart-1
      do j=i+1,npart
         dx(:)=x(:,i)-x(:,j)
         dx=dx-el*nint(dx*eli)
         r=max(tiny,sqrt(sum(dx**2)))
         if (r.lt.rangeu) then
            dr=scaleu*r
            index=dr
            index=max(1,min(index,ntabu-2))
            dr=dr-index
            c1=-dr*(dr-1.0_r8)*(dr-2.0_r8)/6.0_r8
            c2=(dr+1.0_r8)*(dr-1.0_r8)*(dr-2.0_r8)/2.0_r8
            c3=-(dr+1.0_r8)*dr*(dr-2.0_r8)/2.0_r8
            c4=(dr+1.0_r8)*dr*(dr-1.0_r8)/6.0_r8
            uj=c1*utab(index-1)+c2*utab(index) &
              +c3*utab(index+1)+c4*utab(index+2)
            duj=c1*dutab(index-1)+c2*dutab(index) &
               +c3*dutab(index+1)+c4*dutab(index+2)
            d2uj=c1*d2utab(index-1)+c2*d2utab(index) &
                +c3*d2utab(index+1)+c4*d2utab(index+2)
!index=nint(r*scaleu)
!uj=utab(index)
!duj=dutab(index)
!d2uj=d2utab(index)
            u=u+uj
            du(:,i)=du(:,i)+dx(:)*duj/(r*sqrt(mass(i)))
            du(:,j)=du(:,j)-dx(:)*duj/(r*sqrt(mass(j)))
            mu=mass(i)*mass(j)/(mass(i)+mass(j))
            d2u=d2u+(d2uj+(ndim-1.0_r8)*duj/r)/mu
         endif
!        if (r.lt.rangev) then
!           dr=scalev*r
!           index=dr
!           index=max(1,min(index,ntabv-2))
!           dr=dr-index
!           c1=-dr*(dr-1.0_r8)*(dr-2.0_r8)/6.0_r8
!           c2=(dr+1.0_r8)*(dr-1.0_r8)*(dr-2.0_r8)/2.0_r8
!           c3=-(dr+1.0_r8)*dr*(dr-2.0_r8)/2.0_r8
!           c4=(dr+1.0_r8)*dr*(dr-1.0_r8)/6.0_r8
!           v=v+c1*vtab(index-1)+c2*vtab(index) &
!             +c3*vtab(index+1)+c4*vtab(index+2)
!        endif
 if (r.le.2.21976_r8) v=v+2993.0_r8*exp(-6.1_r8*r**2)-338.0_r8*exp(-1.17_r8*r**2)
      enddo
   enddo
   is=1
   d2u=d2u+sum(du*du)
   end subroutine getphi

   subroutine getphione(dummy,dummy1,dummy3,dummy2)
   real(kind=r8) :: dummy(:,:),dummy1(:),dummy2
   integer(kind=i4) :: dummy3
   end subroutine getphione
end module potential
