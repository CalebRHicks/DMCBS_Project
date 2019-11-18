module debug

contains
   subroutine printwalk(w,text)
   use stack
   type (walker) :: w
   character(len=*) :: text
   open(unit=6,file='debug.log',position='append')
   write(6,'(/''debug: printwalker'')')
   write(6,'(''message: '',a)') text
   write(6,'(''psil = '',f25.15)') w%psil
   write(6,'(''psil0 = '',f25.15)') w%psil0
   write(6,'(''psigl = '',f25.15)') w%psigl
   write(6,'(''psigl0 = '',f25.15)') w%psigl0
   write(6,'(''d2psi = '',f25.15)') w%d2psi
   write(6,'(''d2psig = '',f25.15)') w%d2psig
   write(6,'(''v = '',f25.15)') w%v
   write(6,'(''vext = '',f25.15)') w%vext
   write(6,'(''weight = '',f25.15)') w%weight
   write(6,'(''wt0 = '',f25.15)') w%wt0
   write(6,'(''is = '',i10)') w%is
   write(6,'(''isg = '',i10)') w%isg
   write(6,'(''irn = '',i20)') w%irn
   write(6,'(''coordinates'')')
   write(6,'(3f25.15)') w%x
   write(6,'(''dpsi'')')
   write(6,'(3f25.15)') w%dpsi
   write(6,'(''dpsig'')')
   write(6,'(3f25.15)') w%dpsig
   write(6,'(''end printwalker''//)')
   close(6)
   end subroutine
end module debug
