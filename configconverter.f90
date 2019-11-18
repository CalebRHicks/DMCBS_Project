        program confs_converter
        implicit none
        integer, parameter :: i4=selected_int_kind(9)
        integer, parameter :: i8=selected_int_kind(15)
        integer, parameter :: r8=selected_real_kind(15,9)
        integer(kind=i4) :: n,i,npart,j,ic,npart1,npart2
        real(kind=r8) :: dum1
        real(kind=r8), allocatable :: x(:,:),sg(:,:),sy(:,:)
        integer(kind=i8) :: dum2
        npart1=10
        npart2=10
        npart=npart1+npart2
        allocate(x(3,npart),sg(3,npart),sy(3,npart))
        read (9,*) n
        do i=1,n
           read (9,*) x
           read (9,*) dum1
           read (9,*) dum2
           write (10,*) npart
           write (10,*) ''
           do j=1,npart1
              write (10,*) 'H',x(:,j)
           enddo
           do j=npart1+1,npart1+npart2
              write (10,*) 'D',x(:,j)
           enddo
        enddo
        end program
