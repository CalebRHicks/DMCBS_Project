module operators
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)
   logical, private, save :: igofr,isofq,irho,irad,ienk,ivext,iqeuc
contains
   subroutine setoperators(ndim,nbin,el,oper,mass,iqeucin,kunit,hbar,etrial)
   use gofr
   use sofq
   use density
   use radii
   use pdistribution
   use vexternal
   use mympi
   use euclidean
   integer(kind=i4) :: ndim,nbin(:),npart
   integer(kind=i4) :: ntabg,ntabq,ntabr,ntabenk,ntabqeuc
   real(kind=r8) :: el,rmax,mass(:),kunit,hbar,etrial
   logical :: oper,iqeucin
   npart=sum(nbin)
   igofr=.false.
   isofq=.false.
   ienk=.false.
   iqeuc=.false.
   irho=.false.
   irad=.false.
   ntabg=0
   ntabq=0
   ntabenk=0
   ntabqeuc=0
   ntabr=0
   rmax=el
   if (myrank().eq.0) then
      if (el.gt.0.0_r8) then
         read (5,*) ntabg
         if (ntabg.gt.0) igofr=.true.
         read (5,*) ntabq
         if (ntabq.gt.0) isofq=.true.
         read (5,*) ntabenk
         if (ntabenk.gt.0) ienk=.true.
         read (5,*) ntabqeuc
         if (ntabqeuc.gt.0) iqeuc=.true.
      else
         read (5,*) ntabr
         if (ntabr.gt.0) then
            irho=.true.
            read (5,*) rmax
         endif
         read (5,*) irad
         read (5,*) ivext
      endif
   endif
   call bcast(igofr)
   call bcast(ntabg)
   call bcast(isofq)
   call bcast(ntabq)
   call bcast(ienk)
   call bcast(ntabenk)
   call bcast(iqeuc)
   call bcast(ntabqeuc)
   call bcast(irho)
   call bcast(ntabr)
   call bcast(rmax)
   call bcast(irad)
   call bcast(ivext)
   if (igofr) call setupgofr(ndim,ntabg,nbin,el,kunit)
   if (isofq) call setupsofq(ndim,nbin,el,ntabq,kunit)
   if (ienk) call setupenk(ndim,ntabenk,el,kunit)
   if (iqeuc) call setupqeuc(ndim,nbin,el,ntabqeuc,kunit,hbar,etrial)
   if (irho) call setuprho(ntabr,nbin,rmax,mass)
   if (irad) call setupradii(ndim,nbin,mass)
   if (ivext) call setupvext(ndim,nbin,mass)
   oper=igofr.or.isofq.or.irho.or.irad.or.ienk.or.ivext.or.iqeuc
   call zerooperators
   iqeucin=iqeuc
   end subroutine setoperators

   subroutine zerooperators
   use gofr
   use sofq
   use density
   use pdistribution
   use euclidean
   if (igofr) call zerogofr
   if (isofq) call zerosofq
   if (irho) call zerorho
   if (ienk) call zeroenk
   if (iqeuc) call zeroqeuc
   end subroutine zerooperators
 
   subroutine addoperators(w)
   use stack
   use gofr
   use sofq
   use density
   use radii
   use vexternal
   use pdistribution
   use euclidean
   type (walker) w
   real(kind=r8) :: wt
   wt=w%weight*exp(w%psil-w%psigl)*w%is*w%wt0
   if (igofr) call addgofr(w%x,wt)
   if (isofq) call addsofq(w%x,wt)
   if (ienk) call addenk(w,wt)
   if (iqeuc) call addqeuc(w)
   if (irho) call addrho(w%x,wt)
   if (irad) call addradii(w%x,wt)
   if (ivext) call addvext(w%x,wt)
   end subroutine addoperators

   subroutine updateoperators
   use gofr
   use sofq
   use density
   use radii
   use pdistribution
   use vexternal
   use euclidean
   if (igofr) call updategofr
   if (isofq) call updatesofq
   if (ienk) call updateenk
   if (iqeuc) call updateqeuc
   if (irho) call updaterho
   if (irad) call updateradii
   if (ivext) call updatevext
   end subroutine updateoperators

   subroutine writeoperators(tau)
   use gofr
   use sofq
   use density
   use pdistribution
   use vexternal
   use euclidean
   real(kind=r8) :: tau
   if (igofr) call writegofr('gofr.out',tau)
   if (isofq) call writesofq('sofq.out',tau)
   if (ienk) call writeenk('nofk.out',tau)
   if (irho) call writerho('rho.out')
   end subroutine writeoperators
end module
