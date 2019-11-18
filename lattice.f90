module lattice
   implicit none
   integer, private, parameter :: i4=selected_int_kind(9)
   integer, private, parameter :: r8=selected_real_kind(15,9)

contains
    function simplecubic(nside)
!
! set up a simple cubic lattice with nside unit cells per side in
! the unit cube
!
    integer(kind=i4) :: nside
    real(kind=r8), dimension(3,nside**3) :: simplecubic
    real(kind=r8), dimension(3,1) :: basis = &
       reshape((/0.25_r8, 0.25_r8, 0.25_r8/),(/3,1/))
    simplecubic=basiscubic(nside,basis)
    return
    end function simplecubic

    function bodycenteredcubic(nside)
!
! set up a body-centered cubic lattice with nside unit cells per side in
! the unit cube
!
    integer(kind=i4) :: nside
    real(kind=r8), dimension(3,2*nside**3) :: bodycenteredcubic
    real(kind=r8), dimension(3,2) :: basis = &
       reshape((/0.25_r8, 0.25_r8, 0.25_r8, 0.75_r8, 0.75_r8, 0.75_r8/) &
       ,(/3,2/))
    bodycenteredcubic=basiscubic(nside,basis)
    end function bodycenteredcubic
   
     function facecenteredcubic(nside)
!
! set up a face-centered cubic lattice with nside unit cells per side in
! the unit cube
!
     integer(kind=i4) :: nside
     real(kind=r8), dimension(3,4*nside**3) :: facecenteredcubic
     real(kind=r8), dimension(3,4) :: basis = &
        reshape((/0.25_r8, 0.25_r8, 0.25_r8, 0.75_r8, 0.75_r8, 0.25_r8, &
        0.75_r8, 0.25_r8, 0.75_r8, 0.25_r8, 0.75_r8, 0.75_r8/),(/3,4/))
     facecenteredcubic=basiscubic(nside,basis)
     return
     end function facecenteredcubic

   function basiscubic(nside,basis)
!
! given a basis in the cell, set up a cubic lattice with nside unit cells
! per side in the unit cube
!
   integer(kind=i4) :: nside
   real(kind=r8), dimension(:,:) :: basis
   real(kind=r8), dimension(3,nside**3*size(basis)/3) :: basiscubic
   real(kind=r8) :: scale
   integer(kind=i4) :: i,i1,i2,i3,ib,nbasis
   nbasis=size(basis)/3
   scale=1.0_r8/nside
   i=0
   do i1=0,nside-1
      do i2=0,nside-1
         do i3=0,nside-1
            do ib=1,nbasis
               i=i+1
               basiscubic(1,i)=(i1+basis(1,ib))*scale
               basiscubic(2,i)=(i2+basis(2,ib))*scale
               basiscubic(3,i)=(i3+basis(3,ib))*scale
            enddo
         enddo
      enddo
   enddo
   end function basiscubic

   function bestcubic(npart)
!
! try simple, body-centered, and face-centered cubic and return npart
! positions in the unit cube that has the fewest vacancies for npart
! particles
!
   integer(kind=i4) :: npart
   real(kind=r8), dimension(3,npart) :: bestcubic
   integer(kind=i4) :: nsc,nbcc,nfcc,nsidesc,nsidebcc,nsidefcc
   real(kind=r8), dimension(:,:), allocatable :: x
   real(kind=r8) :: en
   en=npart
   nsidesc=nint(en**(1.0_r8/3.0_r8))
   if (nsidesc**3.lt.npart) nsidesc=nsidesc+1
   nsc=nsidesc**3
   nsidebcc=nint((en*0.5_r8)**(1.0_r8/3._r8))
   if (nsidebcc**3*2.lt.npart) nsidebcc=nsidebcc+1
   nbcc=nsidebcc**3*2
   nsidefcc=nint((en*0.25_r8)**(1.0_r8/3._r8))
   if (nsidefcc**3*4.lt.npart) nsidefcc=nsidefcc+1
   nfcc=nsidefcc**3*4
   if (nsc.le.nbcc.and.nsc.le.nfcc) then
      allocate(x(3,nsc))
      x = simplecubic(nsidesc)
      else if (nbcc.le.nfcc) then
      allocate(x(3,nbcc))
      x = bodycenteredcubic(nsidebcc)
      else   
      allocate(x(3,nfcc))
      x = facecenteredcubic(nsidefcc)
      endif
   bestcubic(:,:)=x(:,1:npart:1)
   deallocate(x)
   return
   end function bestcubic

end module lattice
