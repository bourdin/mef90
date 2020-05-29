Program TestATClass
#include "../MEF90/mef90.inc"
#include "finclude/petscdef.h"
   Use m_MEF90
   Use m_MEF90_DefMechAT

   IMPLICIT NONE

   PetscInt,Parameter      :: N = 11
   integer                 :: i
   PetscReal, dimension(N) :: alpha = [( (i-1.0_Kr)/(N-1.0_Kr), i=1,N)]

   class(MEF90_DefMechAT_Type), allocatable :: AT

   AT = MEF90_DefMechAT1_Type()
   write(*,*) 'AT1:'
   write(*,*) '       cw: ',AT%cw
   write(*,*) '   aorder: ',AT%aorder
   write(*,*) '   worder: ',AT%worder
   write(*,'(13x,7(A12,2x))') 'alpha','a', 'Da', 'D2a', 'w', 'Dw', 'D2w'
   Do i = 1, N
      write(*,'(13x,7(ES12.5,2x))')alpha(i),AT%a(alpha(i)),AT%Da(alpha(i)),AT%D2a(alpha(i)),AT%w(alpha(i)),AT%Dw(alpha(i)),AT%D2w(alpha(i))
   End Do

   AT = MEF90_DefMechAT2_Type()
   write(*,*) 'AT2:'
   write(*,*) '       cw: ',AT%cw
   write(*,*) '   aorder: ',AT%aorder
   write(*,*) '   worder: ',AT%worder
   write(*,'(13x,7(A12,2x))') 'alpha','a', 'Da', 'D2a', 'w', 'Dw', 'D2w'
   Do i = 1, N
      write(*,'(13x,7(ES12.5,2x))')alpha(i),AT%a(alpha(i)),AT%Da(alpha(i)),AT%D2a(alpha(i)),AT%w(alpha(i)),AT%Dw(alpha(i)),AT%D2w(alpha(i))
   End Do

End Program TestATClass
