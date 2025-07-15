program TestATClass
#include "../MEF90/mef90.inc"
#include "petsc/finclude/petsc.h"
   use m_MEF90
   use m_MEF90_DefMechAT

   implicit none(type, external)

   PetscInt, parameter      :: N = 11
   integer                 :: i
   PetscReal, dimension(N) :: alpha = [((i - 1.0_kr) / (N - 1.0_kr), i=1, N)]

   class(MEF90DefMechAT_Type), allocatable :: AT

   AT = MEF90DefMechAT1_Type()
   write (*, *) trim(AT % type)
   write (*, *) '       cw: ', AT % cw
   write (*, *) '   aorder: ', AT % aorder
   write (*, *) '   worder: ', AT % worder
   write (*, '(13x,7(A12,2x))') 'alpha', 'a', 'Da', 'D2a', 'w', 'Dw', 'D2w'
   do i = 1, N
      write (*, '(13x,7(ES12.5,2x))') alpha(i), AT % a(alpha(i)), AT % Da(alpha(i)), AT % D2a(alpha(i)), AT % w(alpha(i)), AT % Dw(alpha(i)), AT % D2w(alpha(i))
   end do

   AT = MEF90DefMechAT2_Type()
   write (*, *) trim(AT % type)
   write (*, *) '       cw: ', AT % cw
   write (*, *) '   aorder: ', AT % aorder
   write (*, *) '   worder: ', AT % worder
   write (*, '(13x,7(A12,2x))') 'alpha', 'a', 'Da', 'D2a', 'w', 'Dw', 'D2w'
   do i = 1, N
      write (*, '(13x,7(ES12.5,2x))') alpha(i), AT % a(alpha(i)), AT % Da(alpha(i)), AT % D2a(alpha(i)), AT % w(alpha(i)), AT % Dw(alpha(i)), AT % D2w(alpha(i))
   end do

   AT = MEF90DefMechAT1exp_Type(0.1_kr)
   write (*, *) trim(AT % type)
   write (*, *) '       cw: ', AT % cw
   write (*, *) '   aorder: ', AT % aorder
   write (*, *) '   worder: ', AT % worder
   write (*, '(13x,7(A12,2x))') 'alpha', 'a', 'Da', 'D2a', 'w', 'Dw', 'D2w'
   do i = 1, N
      write (*, '(13x,7(ES12.5,2x))') alpha(i), AT % a(alpha(i)), AT % Da(alpha(i)), AT % D2a(alpha(i)), AT % w(alpha(i)), AT % Dw(alpha(i)), AT % D2w(alpha(i))
   end do

   AT = MEF90DefMechATKKL_Type()
   write (*, *) trim(AT % type)
   write (*, *) '       cw: ', AT % cw
   write (*, *) '   aorder: ', AT % aorder
   write (*, *) '   worder: ', AT % worder
   write (*, '(13x,7(A12,2x))') 'alpha', 'a', 'Da', 'D2a', 'w', 'Dw', 'D2w'
   do i = 1, N
      write (*, '(13x,7(ES12.5,2x))') alpha(i), AT % a(alpha(i)), AT % Da(alpha(i)), AT % D2a(alpha(i)), AT % w(alpha(i)), AT % Dw(alpha(i)), AT % D2w(alpha(i))
   end do

   AT = MEF90DefMechATLinSoft_Type(0.1_kr)
   write (*, *) trim(AT % type)
   write (*, *) '       cw: ', AT % cw
   write (*, *) '   aorder: ', AT % aorder
   write (*, *) '   worder: ', AT % worder
   write (*, '(13x,7(A12,2x))') 'alpha', 'a', 'Da', 'D2a', 'w', 'Dw', 'D2w'
   do i = 1, N
      write (*, '(13x,7(ES12.5,2x))') alpha(i), AT % a(alpha(i)), AT % Da(alpha(i)), AT % D2a(alpha(i)), AT % w(alpha(i)), AT % Dw(alpha(i)), AT % D2w(alpha(i))
   end do

end program TestATClass
