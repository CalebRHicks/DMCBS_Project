module atom
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   integer(kind=i4), private, save :: npart,ndim,ntab
   integer(kind=i4), private, save :: nup,ndn
   integer(kind=i4), private, save, allocatable :: nbin(:)
   real(kind=r8), private, save, allocatable :: mass(:)
   real(kind=r8), private, save, allocatable :: range(:),scale(:)
   real(kind=r8), private, save, allocatable :: phi(:,:),dphi(:,:),d2phi(:,:)
contains
   subroutine setuppsiatom(rho,npartin,nbinin,massin,ndin)
   use mympi
   integer(kind=i4) :: npartin,nbinin(:),ndin,i,iorb,norb
   real(kind=r8) :: rho,massin(:),r
   rho=0.0_r8
   npart=npartin
   ndim=ndin
   allocate(nbin(2))
   nbin=nbinin
   allocate(mass(npart))
   mass(1:nbin(1))=massin(1)
   if (nbin(1).ne.npart) mass(nbin(1)+1:npart)=massin(2)
   read (80,*) ntab
   ntab=ntab-1
! the order of radial orbitals is:
! 1s, 1p, 2s, 1d, 2p, 1f, 3s, 2d, 1g, 3p, 2f
   norb=1 ! 1s
   if (nbin(2).gt.1) norb=2 ! 1p
   if (nbin(2).gt.4) norb=3 ! 2s
   if (nbin(2).gt.5) norb=4 ! 1d
   if (nbin(2).gt.10) norb=5 ! 2p
   if (nbin(2).gt.13) norb=6 ! 1f
   if (nbin(2).gt.20) norb=7 ! 3s
   if (nbin(2).gt.21) norb=8 ! 2d
   if (nbin(2).gt.26) norb=9 ! 1g
   if (nbin(2).gt.35) norb=10 ! 3p
   if (nbin(2).gt.38) norb=11 ! 2f
   allocate(phi(norb,0:ntab),dphi(norb,0:ntab),d2phi(norb,0:ntab))
   allocate(range(norb),scale(norb))
   phi=0.0_r8
   dphi=0.0_r8
   d2phi=0.0_r8
   if (myrank().eq.0) then
      do iorb=1,norb
         do i=0,ntab
            read (80+iorb-1,('(4e25.15)')) r,phi(iorb,i),dphi(iorb,i),d2phi(iorb,i)
         enddo
         range(iorb)=r
         scale(iorb)=real(ntab)/range(iorb)
      enddo
! here I use 1s-light equal to 1s-heavy, fix it
   endif
   call bcast(phi)
   call bcast(dphi)
   call bcast(d2phi)
   call bcast(range)
   call bcast(scale)
   end subroutine setuppsiatom

   subroutine psiatom(x,phil,is,dphi,d2phi)
   use matrixmod
   real(kind=r8) :: x(:,:),phil,dphi(:,:),d2phi
   real(kind=r8) :: xcm(3),massi
   real(kind=r8) :: smati(npart,npart),dph(npart,3,npart)
   real(kind=r8) :: d2ph(npart,npart),d2cor,c1,c2,c3,c4,d2c(npart)
   real(kind=r8) :: dphisum(3),d2t(npart),d2tot,small(npart,3,npart)
   integer(kind=i4) :: is,i,ic,i1,i2
   massi=1.0_r8/sum(mass(:))
   xcm(1)=(sum(x(1,:)*mass(:)))*massi
   xcm(2)=(sum(x(2,:)*mass(:)))*massi
   xcm(3)=(sum(x(3,:)*mass(:)))*massi
   x(1,:)=x(1,:)-xcm(1)
   x(2,:)=x(2,:)-xcm(2)
   x(3,:)=x(3,:)-xcm(3)
   smati=0.0_r8
   dph=0.0_r8
   d2ph=0.0_r8
   call onebod(x(:,1:nbin(1)),nbin(1),smati(1:nbin(1),1:nbin(1)) &
          ,dph(1:nbin(1),:,1:nbin(1)),d2ph(1:nbin(1),1:nbin(1)),1)
   call onebod(x(:,nbin(1)+1:nbin(1)+nbin(2)),nbin(2)            &
          ,smati(1+nbin(1):npart,1+nbin(1):npart),dph(1+nbin(1):npart,:,1+nbin(1):npart) &
          ,d2ph(1+nbin(1):npart,1+nbin(1):npart),2)
   call matinv(smati,phil,is,npart)
   do i=1,npart
      dphi(1,i)=sum(smati(:,i)*dph(i,1,:))
      dphi(2,i)=sum(smati(:,i)*dph(i,2,:))
      dphi(3,i)=sum(smati(:,i)*dph(i,3,:))
   enddo
   dphisum(1)=sum(dphi(1,:))
   dphisum(2)=sum(dphi(2,:))
   dphisum(3)=sum(dphi(3,:))
   d2t=0.0_r8
   do i=1,npart
      dphi(1,i)=dphi(1,i)-dphisum(1)*massi*mass(i)
      dphi(2,i)=dphi(2,i)-dphisum(2)*massi*mass(i)
      dphi(3,i)=dphi(3,i)-dphisum(3)*massi*mass(i)
      d2t(i)=sum(smati(:,i)*d2ph(i,:))
   enddo
   d2tot=sum(d2t)
   small=reshape(matmul(smati,reshape(dph,(/npart,3*npart/))),shape(small))
   d2c=0.0_r8
   do ic=1,3
      do i1=1,npart
         c1=small(i1,ic,i1)
         do i2=1,npart
            if (i2.ne.i1) then
               c2=small(i2,ic,i2)
               c3=small(i1,ic,i2)
               c4=small(i2,ic,i1)
               d2c(i1)=d2c(i1)+c1*c2-c3*c4
            endif
         enddo
      enddo
   enddo
   d2cor=0.5_r8*sum(d2c)
   do i=1,npart
      d2t(i)=(1.0_r8-2.0_r8*massi*mass(i))*d2t(i)-2.0_r8*massi*mass(i)*d2c(i) &
            +massi**2*mass(i)**2*(2.0_r8*d2cor+d2tot)
   enddo
   dphi(1,:)=dphi(1,:)/sqrt(mass(:))
   dphi(2,:)=dphi(2,:)/sqrt(mass(:))
   dphi(3,:)=dphi(3,:)/sqrt(mass(:))
   d2phi=sum(d2t/mass)
   end subroutine psiatom

   subroutine onebod(x,np,smati,dph,d2ph,spin)
   integer(kind=i4) :: np,i,j,l,spin
   real(kind=r8), dimension(:,:) :: x,smati,d2ph
   real(kind=r8), dimension(:,:,:) :: dph
   real(kind=r8) :: dx(3),r,ri,ril,phir,dphr,dphir(3),d2phir
   real(kind=r8) :: ylm,dylm(3),d2ylm
   real(kind=r8) :: f,df,d2f,x2,y2,z2,x4,y4,z4
   smati=0.0_r8
   dph=0.0_r8
   d2ph=0.0_r8
   do i=1,np
      dx=x(:,i)
      r=max(1e-8_r8,sqrt(sum(dx**2)))
      ri=1.0_r8/r
      do j=1,np
         if (j.eq.1) then
            l=0
         elseif (j.le.4) then
            l=1
         elseif (j.le.5) then
            l=0
         elseif (j.le.10) then
            l=2
         elseif (j.le.13) then
            l=1
         elseif (j.le.20) then
            l=3
         elseif (j.le.21) then
            l=0
         elseif (j.le.26) then
            l=2
         elseif (j.le.35) then
            l=4
         elseif (j.le.38) then
            l=1
         elseif (j.le.45) then
            l=3
         elseif (j.gt.45) then
            write (6,*) 'error in onebod'
            stop
         endif
         ril=ri**l
         call radialfun(r,f,df,d2f,j,spin)
         phir=f*ril
         dphr=df*ril-l*f*ril*ri
         d2phir=(d2f-2.0_r8*l*df*ri+ri**2*l*f*(1.0_r8+l))*ril+2.0_r8*dphr*ri
         dphir=dphr*dx*ri
         if (j.gt.10) then
            x2=dx(1)*dx(1)
            y2=dx(2)*dx(2)
            z2=dx(3)*dx(3)
            x4=x2*x2
            y4=y2*y2
            z4=z2*z2
         endif
         select case (j)
         case (1,5,21) ! 1s,2s,3s
            ylm=1.0_r8
            dylm=0.0_r8
            d2ylm=0.0_r8
         case (2,11,36) ! 1p,2p
            ylm=dx(1)
            dylm(1)=1.0_r8
            dylm(2)=0.0_r8
            dylm(3)=0.0_r8
            d2ylm=0.0_r8
         case (3,12,37)
            ylm=dx(2)
            dylm(1)=0.0_r8
            dylm(2)=1.0_r8
            dylm(3)=0.0_r8
            d2ylm=0.0_r8
         case (4,13,38)
            ylm=dx(3)
            dylm(1)=0.0_r8
            dylm(2)=0.0_r8
            dylm(3)=1.0_r8
            d2ylm=0.0_r8
         case (6,22) ! 1d,2d
            ylm=2.0_r8*dx(3)**2-dx(1)**2-dx(2)**2
            dylm(1)=-2.0_r8*dx(1)
            dylm(2)=-2.0_r8*dx(2)
            dylm(3)=4.0_r8*dx(3)
            d2ylm=0.0_r8
         case (7,23)
            ylm=dx(1)*dx(3)
            dylm(1)=dx(3)
            dylm(2)=0.0_r8
            dylm(3)=dx(1)
            d2ylm=0.0_r8
         case (8,24)
            ylm=dx(2)*dx(3)
            dylm(1)=0.0_r8
            dylm(2)=dx(3)
            dylm(3)=dx(2)
            d2ylm=0.0_r8
         case (9,25)
            ylm=dx(1)**2-dx(2)**2
            dylm(1)=2.0_r8*dx(1)
            dylm(2)=-2.0_r8*dx(2)
            dylm(3)=0.0_r8
            d2ylm=0.0_r8
         case (10,26)
            ylm=dx(1)*dx(2)
            dylm(1)=dx(2)
            dylm(2)=dx(1)
            dylm(3)=0.0_r8
            d2ylm=0.0_r8
         case (14,39) ! ! 1f,2f
            ylm=dx(3)*(2.0_r8*z2-3.0_r8*x2-3.0_r8*y2)
            dylm(1)=-6.0_r8*dx(1)*dx(3)
            dylm(2)=-6.0_r8*dx(2)*dx(3)
            dylm(3)=6.0_r8*z2-3.0_r8*x2-3.0_r8*y2
            d2ylm=0.0_r8
         case (15,40) 
            ylm=dx(2)*(3.0_r8*x2-y2)
            dylm(1)=6.0_r8*dx(1)*dx(2)
            dylm(2)=3.0_r8*x2-3.0_r8*y2
            dylm(3)=0.0_r8
            d2ylm=0.0_r8
         case (16,41) 
            ylm=dx(1)*(x2-3.0_r8*y2)
            dylm(1)=3.0_r8*x2-3.0_r8*y2
            dylm(2)=-6.0_r8*dx(1)*dx(2)
            dylm(3)=0.0_r8
            d2ylm=0.0_r8
         case (17,42) 
            ylm=dx(3)*(x2-y2)
            dylm(1)=2.0_r8*dx(1)*dx(3)
            dylm(2)=-2.0_r8*dx(2)*dx(3)
            dylm(3)=x2-y2
            d2ylm=0.0_r8
         case (18,43)
            ylm=dx(1)*dx(2)*dx(3)
            dylm(1)=dx(2)*dx(3)
            dylm(2)=dx(1)*dx(3)
            dylm(3)=dx(1)*dx(2)
            d2ylm=0.0_r8
         case (19,44)
            ylm=dx(2)*(4.0_r8*z2-x2-y2)
            dylm(1)=-2.0_r8*dx(1)*dx(2)
            dylm(2)=4.0_r8*z2-x2-3.0_r8*y2
            dylm(3)=8.0_r8*dx(2)*dx(3)
            d2ylm=0.0_r8
         case (20,45)
            ylm=dx(1)*(4.0_r8*z2-x2-y2)
            dylm(1)=4.0_r8*z2-3.0_r8*x2-y2
            dylm(2)=-2.0_r8*dx(1)*dx(2)
            dylm(3)=8.0_r8*dx(1)*dx(3)
            d2ylm=0.0_r8
         case (27) ! 1g
            ylm=8.0_r8*z4-24.0_r8*x2*z2-24.0_r8*y2*z2+3.0_r8*x4+6.0_r8*x2*y2+3.0_r8*y4
            dylm(1)=-60.0_r8*dx(1)*z2+12.0_r8*(x2+y2+z2)*dx(1)
            dylm(2)=-60.0_r8*dx(2)*z2+12.0_r8*(x2+y2+z2)*dx(2)
            dylm(3)=80.0_r8*z2*dx(3)-48.0_r8*dx(3)*(x2+y2+z2)
            d2ylm=0.0_r8
         case (28)
            ylm=dx(1)*dx(3)*(4.0_r8*z2-3.0_r8*x2-3.0_r8*y2)
            dylm(1)=4.0_r8*dx(3)*z2-9.0_r8*x2*dx(3)-3.0_r8*dx(3)*y2
            dylm(2)=-6.0_r8*dx(1)*dx(2)*dx(3)
            dylm(3)=12.0_r8*dx(1)*z2-3.0_r8*dx(1)*x2-3.0_r8*dx(1)*y2
            d2ylm=0.0_r8
         case (29)
            ylm=dx(2)*dx(3)*(4.0_r8*z2-3.0_r8*x2-3.0_r8*y2)
            dylm(1)=-6.0_r8*dx(1)*dx(2)*dx(3)
            dylm(2)=4.0_r8*dx(3)*z2-3.0_r8*x2*dx(3)-9.0_r8*dx(3)*y2
            dylm(3)=12.0_r8*dx(2)*z2-3.0_r8*dx(2)*x2-3.0_r8*dx(2)*y2
            d2ylm=0.0_r8
         case (30)
            ylm=(x2-y2)*(6.0_r8*z2-x2-y2)
            dylm(1)=12.0_r8*dx(1)*z2-4.0_r8*dx(1)*x2
            dylm(2)=-12.0_r8*dx(2)*z2+4.0_r8*dx(2)*y2
            dylm(3)=12.0_r8*(x2-y2)*dx(3)
         case (31)
            ylm=dx(1)*dx(2)*(6.0_r8*z2-x2-y2)
            dylm(1)=6.0_r8*dx(2)*z2-3.0_r8*dx(2)*x2-dx(2)*y2
            dylm(2)=6.0_r8*dx(1)*z2-3.0_r8*dx(1)*y2-dx(1)*x2
            dylm(3)=12.0_r8*dx(1)*dx(2)*dx(3)
         case (32)
            ylm=dx(1)*dx(3)*(x2-3.0_r8*y2)
            dylm(1)=3.0_r8*x2*dx(3)-3.0_r8*dx(3)*y2
            dylm(2)=-6.0_r8*dx(1)*dx(2)*dx(3)
            dylm(3)=dx(1)*(x2-3.0_r8*y2)
         case (33)
            ylm=dx(2)*dx(3)*(3.0_r8*x2-y2)
            dylm(1)=6.0_r8*dx(1)*dx(2)*dx(3)
            dylm(2)=3.0_r8*x2*dx(3)-3.0_r8*y2*dx(3)
            dylm(3)=dx(2)*(3.0_r8*x2-y2)
         case (34)
            ylm=x2*(x2-3.0_r8*y2)-y2*(3.0_r8*x2-y2)
            dylm(1)=4.0_r8*dx(1)*x2-12.0_r8*dx(1)*y2
            dylm(2)=-12.0_r8*dx(2)*x2+4.0_r8*dx(2)*y2
            dylm(3)=0.0_r8
         case (35)
            ylm=dx(1)*dx(2)*(x2-y2)
            dylm(1)=3.0_r8*x2*dx(2)-dx(2)*y2
            dylm(2)=dx(1)*x2-3.0_r8*dx(1)*y2
            dylm(3)=0.0_r8
         case default
            write (6,*) 'I can use only 45 particles by now for each spin!!!'
            stop
         end select
         smati(i,j)=phir*ylm
         dph(i,:,j)=dphir(:)*ylm+phir*dylm(:)
         d2ph(i,j)=d2phir*ylm+phir*d2ylm+2.0_r8*dot_product(dphir,dylm)
      enddo
   enddo
   return
   end subroutine onebod

   subroutine radialfun(r,f,df,d2f,j,spin)
   use mympi
   real(kind=r8) :: r,f,df,d2f
   real(kind=r8) :: ddr,a1,a2,a3,a4
   integer(kind=i4) :: index,j,itab,spin,tmp
   tmp=spin
   if (j.eq.1) then
      itab=1
   elseif (j.le.4) then
      itab=2
   elseif (j.le.5) then
      itab=3
   elseif (j.le.10) then
      itab=4
   elseif (j.le.13) then
      itab=5
   elseif (j.le.20) then
      itab=6
   elseif (j.le.21) then
      itab=7
   elseif (j.le.26) then
      itab=8
   elseif (j.le.35) then
      itab=9
   elseif (j.le.38) then
      itab=10
   elseif (j.le.45) then
      itab=11
   elseif (j.gt.45) then
      write (6,*) 'error in radialfun!'
      stop
   endif
   if (r.gt.range(itab)) then
      f=0.0_r8
      df=0.0_r8
      d2f=0.0_r8
      return
   endif
   ddr=scale(itab)*r
   index=ddr
   index=max(1,min(index,ntab-2))
   ddr=ddr-index
   a1=-ddr*(ddr-1.0_r8)*(ddr-2.0_r8)/6.0_r8
   a2=(ddr+1.0_r8)*(ddr-1.0_r8)*(ddr-2.0_r8)/2.0_r8
   a3=-(ddr+1.0_r8)*ddr*(ddr-2.0_r8)/2.0_r8
   a4=(ddr+1.0_r8)*ddr*(ddr-1.0_r8)/6.0_r8
   f=a1*phi(itab,index-1)+a2*phi(itab,index)+a3*phi(itab,index+1)+a4*phi(itab,index+2)
   df=a1*dphi(itab,index-1)+a2*dphi(itab,index)+a3*dphi(itab,index+1)+a4*dphi(itab,index+2)
   d2f=a1*d2phi(itab,index-1)+a2*d2phi(itab,index)+a3*d2phi(itab,index+1)+a4*d2phi(itab,index+2)
   end subroutine radialfun
end module atom
