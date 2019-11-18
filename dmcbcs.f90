   program dmcbcs
   use stack
   use lattice
   use random
   use random2
   use wavefunction
   use estimator
   use step
   use mympi
   use potential
   use optimizer
   use gofr
   use sofq
   use operators
   use euclidean
   implicit none
   integer, parameter :: i4=selected_int_kind(9)
   integer, parameter :: i8=selected_int_kind(15)
   integer, parameter :: r8=selected_real_kind(15,9)
   real(kind=r8) :: stmax ! population limit
   integer(kind=i8) :: irn,irn0
   integer(kind=i4) :: npart,nwalk,idmc,neq,nav,nstep,i,j,k,it,ndim
   integer(kind=i4) :: istin,istout,ntab,n1,ito,nt,nbin(2)
   real(kind=r8) :: hbar,dt,etrial,popup,popdn
   real(kind=r8) :: el,pi,rn(1),dcut
   real(kind=r8) :: dummy,eunit,nfac
   integer(kind=i4), allocatable :: numfig(:)
   real :: time0,time1,timebl,timebl0,timeop
   real(kind=r8) :: t0,t0t,tbl,tblt,top,topt,tot,tott
   real(kind=r8) :: wtcut,mass(2)
   character(len=90), dimension(:), allocatable :: answer
   logical :: isite,empty,iout,operat,units,fn
   real(kind=r8), allocatable :: xtemp(:)
   real(kind=r8) :: mu,ene,rho,kunit
   real(kind=r8) :: tau,enet,ene2,ene2t,psi,psit,psi2,psi2t,psig,psigt,psig2,psig2t
   real(kind=r8) :: sgn,sgnt,wt,tau0
   real(kind=r8) :: dobal,wtsum,wttot
   integer(kind=i4) :: optiter,iters
   logical :: euclid,ivmc
   character(len=9) :: string
   integer(kind=i4) :: counter
   call init0 ! mpi initialization that needs to be done before reading
   if (myrank().eq.0) then
      read (5,*) irn     !random number seed for sites
      irn0=irn
      read (5,*) idmc    !-1=opt, 0=vmc,1=dmc local weigth
      read (5,*) isite   !if true take initial walkers from lattice sites
      read (5,*) nwalk   !default number of walkers
      read (5,*) neq     !equilibration blocks
      read (5,*) nav     !averaging blocks
      read (5,*) nstep   !steps per block
      read (5,*) dt      !time step
      read (5,*) etrial  !trial energy
      read (5,*) units   !use units for the energy
      read (5,*) dcut    !drift cutoff
      read (5,*) fn      !fixed-node
      read (5,*) popup   !kill walkers if more than nwalk*popup
      read (5,*) popdn   !create walkers if less than nwalk*popdn
      read (5,*) dobal   !rebalance if walker number fluctuate more than this (per cpu)
      read (5,*) wtcut   !cutoff for weight
      read (5,*) ntab    !number of points for tables
      read (5,*) hbar    !hbar^2/2 without the mass
      read (5,*) ndim    !dimensions
      read (5,*) npart   !total particle number
      read (5,*) nbin(1),mass(1) !spin up
      read (5,*) nbin(2),mass(2) !spin up
      if (sum(nbin).ne.npart) then
         write (6,*) 'wrong number of particles'
         stop
      endif
   endif
   call bcast(idmc)
   call bcast(isite)
   call bcast(nwalk)
   call bcast(neq)
   call bcast(nav)
   call bcast(nstep)
   call bcast(dt)
   call bcast(etrial)
   call bcast(dcut)
   call bcast(fn)
   call bcast(popup)
   call bcast(popdn)
   call bcast(dobal)
   call bcast(wtcut)
   call bcast(ntab)
   call bcast(hbar)
   call bcast(ndim)
   call bcast(npart)
   call bcast(nbin)
   call bcast(mass)
   mu=mass(1)*mass(2)/(mass(1)+mass(2))
   counter=10
   if (myrank().eq.0) then
      iout=.true.
!      if (idmc.eq.-1.or.idmc.eq.-2) iout=.false.
      select case (idmc)
         case (-2)
            write (6,'(''Test Monte Carlo Run and optimization'')')
         case (-1)
            write (6,'(''Variational Monte Carlo Run and optimization'')')
         case (0)
            write (6,'(''Variational Monte Carlo Run'')')
         case (1)
            write (6,'(''Diffusion Monte Carlo Run -- elocal weight'')')
         case (2)
            write (6,'(''Diffusion Monte Carlo Run -- Metropolis'')')
         case (10)
            write (6,'(''Diffusion Monte Carlo Run -- without drift'')')
         case default
            write (6,'(''Illegal idmc value -- main'')')
            call abort
      end select
      write (6,'(''!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!'')')
      write (6,'(''Simulation parameters:'')')
      write (6,'(''random number seed ='',t30,i20)') irn
      write (6,'(''start from sites ='',t40,l10)') isite
      write (6,'(''nwalk ='',t40,i10)') nwalk
      write (6,'(''equilibration blocks ='',t40,i10)') neq
      write (6,'(''averaging blocks ='',t40,i10)') nav
      write (6,'(''nsteps/blocks ='',t40,i10)') nstep
      write (6,'(''time step ='',t35,f15.10)') dt
      write (6,'(''trial energy ='',t40,f10.5)') etrial
      write (6,'(''drift cutoff ='',t40,f10.5)') dcut
      write (6,'(''fixed-node ='',t40,l10)') fn
      write (6,'(''drastic population control up ='',t40,f10.5)') popup
      write (6,'(''drastic population control down ='',t40,f10.5)') popdn
      write (6,'(''rebalance if walkers fluctuate more than '',t45,f6.2,'' %'')') dobal
      write (6,'(''weight cut off ='',t40,f10.5)') wtcut
      write (6,'(''table size ='',t40,i10)') ntab
      write (6,'(''!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!'')')
      write (6,'(''System parameters:'')')
      write (6,'(''hbar^2/2 (without the mass) ='',t40,f10.5)') hbar
      write (6,'(''number of dimensions ='',t40,i10)') ndim
      write (6,'(''npart ='',t40,i10)') npart
      write (6,'(''number of spin up ='',t40,i10)') nbin(1)
      write (6,'(''number of spin down ='',t40,i10)') nbin(2)
      write (6,'(''mass of up particles ='',t40,f10.5)') mass(1)
      write (6,'(''mass of down particles ='',t40,f10.5)') mass(2)
   endif
   call bcast(iout)
   call setpsi(ndim,nbin,el,mass,ntab)
   call setphi(ndim,ntab,nbin,el,mass,hbar/(2.0_r8*mu))
   if (myrank().eq.0) then
      write (6,'(''!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!'')')
      if (el.ne.0.0_r8.and.units) then
!generalized units for 2D and 3D (not 1D) -Galea
         rho=npart/el**(ndim)
         pi=4.0_r8*atan(1.0_r8)
         kunit=(ndim*pi**(ndim-1)*rho)**(1.0_r8/ndim)
         eunit=3.0_r8/(8.0_r8-ndim)*hbar/(2.0_r8*mu)*kunit**2
         write (6,'(''k_Fermi (unbalanced) ='',t40,f10.5)') kunit
         write (6,'(''E_Fermi (unbalanced) ='',t40,f10.5)') eunit
!        eunit=3.0_r8/5.0_r8*hbar/(2.0_r8*mu)*(6.0_r8*pi**2*rho)**(2.0_r8/3.0_r8)
!        write (6,'(''E_Fermi (polarized) ='',t40,f10.5)') eunit
      else
         eunit=1.0_r8
         kunit=1.0_r8
      endif
      write (6,'(''Unit of energy ='',t40,f10.5)') eunit
      nfac=1.0_r8
      etrial=etrial*eunit
      if (el.ne.0.0_r8) then
         etrial=etrial*npart
         nfac=1.0_r8/npart
      endif
   endif
   call bcast(eunit)
   call bcast(etrial)
   call bcast(nfac)
   optiter=1
   call setoperators(ndim,nbin,el,operat,mass,euclid,kunit,hbar,etrial)
   if (idmc.eq.-1.or.idmc.eq.-2) then
      call setupopt(npart,optiter)
      operat=.false.  ! do not compute operators during optimization
   endif
   call setstep(idmc,ndim,npart,nbin,hbar,dt,etrial,el,mass,wtcut,dcut,eunit,fn)
   select case (idmc)
      case (-1:0)
         call setestnum(6)
         if (myrank().eq.0) allocate(answer(6))
         call addest(1,'total energy')
         call addest(2,'kinetic energy')
         call addest(3,'potential energy')
         call addest(4,'external potential')
         call addest(5,'jackson-feenberg energy')
         call addest(6,'acceptance =')
      case (-2,1,10)
         call setestnum(14)
         if (myrank().eq.0) allocate(answer(14))
         call addest(1,'total energy')
         call addest(2,'kinetic energy')
         call addest(3,'potential energy')
         call addest(4,'external potential')
         call addest(5,'node crossing')
         call addest(6,'growth energy')
         call addest(7,'|psit|/psig =')
         call addest(8,'psit =')
         call addest(9,'psig =')
         call addest(10,'sign(psit) =')
         call addest(11,'wt =')
         call addest(12,'wtg =')
         call addest(13,'sign(wtg) =')
         call addest(14,'sign(wtg)-weighted =')
      case (2)
         call setestnum(13)
         if (myrank().eq.0) allocate(answer(13))
         call addest(1,'total energy')
         call addest(2,'kinetic energy')
         call addest(3,'potential energy')
         call addest(4,'node crossing')
         call addest(5,'acceptance')
         call addest(6,'|psit|/psig =')
         call addest(7,'psit =')
         call addest(8,'psig =')
         call addest(9,'sign(psit) =')
         call addest(10,'wt =')
         call addest(11,'wtg =')
         call addest(12,'sign(wtg) =')
         call addest(13,'sign(wtg)-weighted =')
   end select
   ivmc=idmc.eq.0.or.idmc.eq.-1
   call init1(npart,dobal,ivmc) ! additional mpi initialization here
   stmax=10.0_r8
!  if (nproc().le.8) then
      stmax=stmax/nproc()
!  else
!     stmax=0.5_r8
!  endif
   call create(2,int(nwalk*stmax),npart,ndim,euclid) !create stacks
   allocate(numfig(0:nproc()-1))
   allocate(xtemp(ndim))
   call cpu_time(time0)
   do iters=1,optiter
      if (optiter.gt.1.and.myrank().eq.0) write (6,'(''optimization iteration number '',i10)') iters
      if (optiter.gt.1) call updateparams
      if (myrank().eq.0) then
         if (isite) then   !get initial walkers
            irn=irn0
            call setrn(irn)    !set seed for sites
            write (6,'(''Initial walkers from random configs'')')
!            write (6,'(''Initial walkers from sites'')')
!            if (el.ne.0.0_r8) then
!               w1%x=(bestcubic(npart)-0.5_r8)*el &
!                  +.02_r8*el*reshape(randn(size(w1%x)),shape(w1%x))
!               w1%x=(w1%x-el*nint(w1%x/el))
!            else
               w1%x=reshape(el*(randn(size(w1%x))-0.5_r8),shape(w1%x))
!            endif
            n1=nwalk
            tau0=0.0_r8
            write (6,'(''Configurations at tau ='',t35,f15.10)') tau0
         else
            write (6,'(''Initial walkers from file'')')
            rewind 9
            read (9,'(i10,f15.8)') n1,tau0
            write (6,'(''Number of walkers ='',t40,i10)') n1
            write (6,'(''Configurations at tau ='',t35,f15.10)') tau0
         endif
      endif
      call bcast(n1)
      call bcast(tau0)
      tau=tau0
      do k=1,n1
         if (myrank().eq.0) then
            if (isite) then
               do i=1,npart
                  rn=randn(1)
                  j=npart*rn(1)+1
                  xtemp(:)=w1%x(:,i)
                  w1%x(:,i)=w1%x(:,j)
                  w1%x(:,j)=xtemp(:)
               enddo
               w1%weight=1.0_r8
               w1%irn=irn
               call ran2(dummy,w1%irn)
               irn=w1%irn
               if (euclid) then
                  w1%x0=w1%x
               else
                  w1%x0=0.0_r8
               endif
            elseif (iters.eq.1) then
               read (9,'(6e15.7)') w1%x
               read (9,'(e15.7,2e25.15)') w1%weight,w1%psil0,w1%psigl0
!              read (9,'(e15.7)') w1%weight
               if (euclid) then
                  w1%x0=w1%x
               else
                  w1%x0=0.0_r8
               endif
               read (9,'(i20)') w1%irn
            endif
            if (iters.eq.1) call push(1,w1)
         endif
         call barrier
         if (iters.eq.1) then
            ito=mod(k-1,nproc())
            call movewalkers(1,1,0,ito,k)
         endif
      enddo
      istin=1
      istout=2
      n1=numstack(istin)
      ene=0.0_r8
      ene2=0.0_r8
      psi=0.0_r8
      psi2=0.0_r8
      psig=0.0_r8
      psig2=0.0_r8
      sgn=0.0_r8
      call zerest
      call barrier
      if (iters.eq.1) then
         wtsum=0.0_r8
         do i=1,n1
            call pop(istin,w1,empty)
            if (empty) exit
            call hpsi(w1)
!write (6,*) ' calling chkder '
!call chkder(w1,2.e-4_r8,1.e-4_r8,mass,nbin)
!write (6,*) ' calling chkderparams'
!if (idmc.eq.-1) call chkderp(w1,2.e-5_r8)
!stop
            ene=ene+(-hbar*w1%d2psi+w1%v)/eunit/npart
            ene2=ene2+((-hbar*w1%d2psi+w1%v)/eunit/npart)**2
            psi=psi+exp(2.0_r8*w1%psil)
            psi2=psi2+exp(4.0_r8*w1%psil)
            psig=psig+exp(2.0_r8*w1%psigl)
            psig2=psig2+exp(4.0_r8*w1%psigl)
            sgn=sgn+1.0_r8*w1%is
            w1%wt0=w1%is
            w1%weight=exp(w1%psigl-w1%psil0)
 w1%weight=1.0_r8
            wtsum=wtsum+w1%weight
            if (idmc.eq.-1.or.idmc.eq.0) w1%wt0=1.0_r8*w1%is
            if (tau.eq.0.0_r8.and.idmc.ne.-1.and.idmc.ne.0) then
               wt=exp(-(-hbar*w1%d2psig+w1%v-etrial)*dt)
               call addval(1,(-hbar*w1%d2psi+w1%v)/eunit/npart,wt)
            endif
            if (euclid) then
               call addqeuc(w1)
            endif
            call push(istout,w1)
         enddo
         call barrier
         if (euclid.and.idmc.ge.1) then
            call updateqeuc
            call writeqeuc2('euclidean','euclidean_sp','euclidean_tau',0.0_r8)
            call writeqeuc('euclideanall','euclideanall_sp','euclideanall_tau',0.0_r8)
         endif
         istout=istin
         istin=3-istin
         call addall(ene,enet)
         call addall(ene2,ene2t)
         call addall(psi,psit)
         call addall(psi2,psi2t)
         call addall(psig,psigt)
         call addall(psig2,psig2t)
         call addall(sgn,sgnt)
         call addall(n1,nt)
         call addall(wtsum,wttot)
         if (myrank().eq.0) then
            ene=enet/nt
            ene2=ene2t/nt
            ene2=sqrt(abs(ene2-ene*ene)/nt)
            write(6,'(''Energy of initial configurations = '',e12.5,'' +- '',e12.5)') ene,ene2
            psi=psit/nt
            psi2=psi2t/nt
            psi2=sqrt(abs(psi2-psi*psi)/nt)
            psig=psigt/nt
            psig2=psig2t/nt
            psig2=sqrt(abs(psig2-psig*psig)/nt)
            sgn=sgnt/nt
            write(6,'(''Psi**2 of initial configurations  = '',e12.5,'' +- '',e12.5)') psi,psi2
            write(6,'(''Psig**2 of initial configurations = '',e12.5,'' +- '',e12.5)') psig,psig2
            write(6,'(''sum of sign of Phi                = '',e12.5)') sgn
         endif
         if (tau.eq.0.0_r8.and.idmc.ne.-1.and.idmc.ne.0) then
            call update   !collect block averages
            if (myrank().eq.0) then
               answer=resstring(tau)
               write (6,'(a90)') answer(1)
            endif
            do while (.true.)
               call pop(istin,w1,empty)
               if (empty) exit
               call addoperators(w1)
               call push(istout,w1)
            enddo
            istout=istin
            istin=3-istin
            call updateoperators
            call writeoperators(tau)
         endif
      endif
      if (idmc.ge.1) then
         n1=numstack(istin)
         call addall(n1,nt)
         call bcast(nt)
         call bcast(wttot)
         if (myrank().eq.0) then
            write(6,'(''average weight at the beginning = '',t40,f20.5)') wttot/real(nt)
            write(6,'(''NOT doing branching with factor = '',t40,f25.15)') real(nt)/wttot
            write(6,'(''nwalk at the beginning = '',t40,i10)') nt
         endif
!         call branch(istin,istout,nt/wttot)  ! I commented this out - Galea
         call barrier
         call checkpop(istin)
         n1=numstack(istin)
         call addall(n1,nt)
         call bcast(nt)
         if (myrank().eq.0) write(6,'(''nwalk after initial branching = '',t40,i10)') nt
      endif
      call zerest
      call cpu_time(time1)
      t0=time1-time0
      call addall(t0,t0t)
      if (myrank().eq.0) then
         write (6,'(''!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!'')')
         write (6,'(''Time for starting (sec.) = '',t35,f10.2)') t0t/nproc()
      endif
      call barrier
      it=0
      call zerotimebal
      do i=1,nav+neq
         if (i.eq.neq+1) then
            if (myrank().eq.0.and.iout) write (6,'(/,''equilibration done'')')
            call zerest
            call zerooperators
            it=0
         endif
         it=it+1
         call cpu_time(timebl0)
         do j=1,nstep
            if (j.eq.nstep) then
               call step1(istin,istout,.true.)
            else
               call step1(istin,istout,.false.)
            endif
            if (idmc.gt.0.or.idmc.eq.-2) call branch(istin,istout,1.0_r8)
            n1=numstack(istin)
            call addall(n1,nt)
            call bcast(nt)
            if (nt.gt.popup*nwalk.or.nt.lt.popdn*nwalk) call branch(istin,istout,1.0_r8*nwalk/nt)
            if (idmc.ne.0.and.idmc.ne.-1) tau=tau+dt
            call barrier
            call checkpop(istin)
         enddo
         call cpu_time(timebl)
         n1=numstack(istin)
         call addall(n1,nt)
         call bcast(nt)
         call update   !collect block averages
         if (iout.or.i.eq.nav+neq) then
            if (myrank().eq.0) then
               if (idmc.ne.-1.and.idmc.ne.-2) then
                  write (6,'(/,''iteration ='',t30,i14)') it
                  write (6,'(''walker number = '',t30,i14)') nt
               endif
               answer=resstring(tau)
               write (6,'(a90)') (answer(k),k=1,size(answer))
            endif
         endif
         if (operat.and.i.gt.neq) then
            do while (.true.)
               call pop(istin,w1,empty)
               if (empty) exit
               call addoperators(w1)
               call push(istout,w1)
            enddo
            istout=istin
            istin=3-istin
            call updateoperators
         endif
         if (idmc.ge.1.and.euclid) call writeqeuc2('euclidean','euclidean_sp','euclidean_tau',tau-tau0)
         if (idmc.ge.1.and.euclid) call writeqeuc('euclideanall','euclideanall_sp','euclideanall_tau',tau-tau0)
         if (operat.and.i.gt.neq) call writeoperators(tau)
         call cpu_time(timeop)
         tbl=timebl-timebl0
         call addall(tbl,tblt)
         top=timeop-timebl
         call addall(top,topt)
         tot=timeop-time1
         call addall(tot,tott)
         if (myrank().eq.0.and.iout) then
            write (6,'(''Time for block (sec) (min) = '',t35,2f11.3)') tblt/nproc(),tblt/nproc()/60.0_r8
            write (6,'(''Time for operators (sec) (min) = '',t35,2f11.3)') topt/nproc(),topt/nproc()/60.0_r8
            write (6,'(''Elapsed block time (sec) (min) = '',t35,2f11.3)') tott/nproc(),tott/nproc()/60.0_r8
         endif
!         if ((idmc.eq.-1.or.idmc.eq.-2).and.i.gt.neq) then
         if (idmc.eq.-1.or.idmc.eq.-2) then
            do while(.true.)
               call pop(istin,w1,empty)
               if (empty) exit
               call hpsi(w1)
               ene=(-hbar*w1%d2psi+w1%v)/eunit/npart
               if (i.gt.neq) call calder(w1,ene)
               call push(istout,w1)
            enddo
            istout=istin
            istin=3-istout
         endif
!periodically write configurations to file -Galea
!         if (mod(i,6000).eq.0) then
!            counter=counter+1
!            n1=numstack(istin)
!            call gather(n1,numfig)
!            call bcast(numfig)
!            nt=sum(numfig)
!            write(string,'(a7,i2)') 'config.',counter
!            if (myrank().eq.0) write(6,*) 'Writing configurations to file: ', string, i
!            open (8,file=string)
!            if (myrank().eq.0) write(8,'(i10,f15.8)') nt,tau
!            close(8)
!            call barrier
!            do k=0,nproc()-1
!               if (k.eq.myrank()) then
!                  open(unit=8,file=string,position='append')
!                  do j=1,numfig(k)
!                     call pop(istin,w1,empty)
!                     if (.not.empty) then
!                        write (8,'(6e15.7)') w1%x
!                        write (8,'(e15.7,2e25.15)') w1%weight,w1%psil,w1%psigl
!                        write (8,'(i20,2e15.7)') w1%irn
!                        call push(istout,w1)
!                     endif
!                  enddo
!                  close(8)
!               endif
!               call barrier
!            enddo
!            istout=istin
!            istin=3-istout
!         end if
!end write to files section -Galea
      enddo
      if (idmc.eq.-1.or.idmc.eq.-2) call updatepar(ene) ! parameters optimization
      call cpu_time(timeop)
      t0=timeop-time1
      call addall(t0,t0t)
      if (myrank().eq.0) write (6,'(/,''Time for steps ='',f10.3,'' seconds, '',f10.3,'' minutes'')') t0t/nproc(),t0t/60.0/nproc()
      n1=numstack(istin)
      call gather(n1,numfig)
      call bcast(numfig)
      nt=sum(numfig)
   enddo
   if (nav+neq.ne.0) then
      rewind 9
      if (myrank().eq.0) write(9,'(i10,f15.8)') nt,tau
      close(9)
      call barrier
      do k=0,nproc()-1
         if (k.eq.myrank()) then
            open(unit=9,position='append')
            do i=1,numfig(k)
               call pop(istin,w1,empty)
               if (.not.empty) then
                  write (9,'(6e15.7)') w1%x
                  write (9,'(e15.7,2e25.15)') w1%weight,w1%psil,w1%psigl
                  write (9,'(i20,2e15.7)') w1%irn
               endif
            enddo
            close(9)
         endif
         call barrier
      enddo
      call writewtisto
   endif
   call cpu_time(time1)
   t0=time1-time0
   call addall(t0,t0t)
   if (myrank().eq.0) then
      write (6,'(''Total job time ='',f10.3,'' seconds, '',f10.3,'' minutes, ''  &
         ,f10.3,''hours'')') t0t/nproc(),t0t/nproc()/60.0_r8,t0t/nproc()/3600.0_r8
      write (6,'(''Total cpu-time ='',f10.3,'' minutes, '',f10.3,'' hours'')') t0t/60.0_r8,t0t/3600.0_r8
      write (6,'(''Number of cpus ='',i10)') nproc()
   endif
   call barrier
   call printmpilog
   call done
   if (myrank().eq.0) write (6,*) 'Finished!'
   end program dmcbcs
