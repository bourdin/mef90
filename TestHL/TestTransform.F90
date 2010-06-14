Program TestTransform

#include "finclude/petscdef.h"

   Use petsc
   Use m_MEF90
   Use m_VarFrac_Struct
   Implicit NONE

   PetscReal                            :: lambda, mu, mu1, mu2
   PetscReal                            :: theta
   Type(Tens4OS2D)                      :: A, B
   PetscReal, Dimension(:,:), Pointer   :: R

   Allocate(R(2,2))
   Write(*, 100, advance = 'no') '(lambda, mu1, mu2): '
   Read(*,*)    lambda, mu1, mu2
   A = 0.0_Kr
   A%XXXX = lambda + 2.0_Kr * mu1
   A%XXXY = 0.0_Kr
   A%XXYY = lambda
   A%XYXY = mu2
   A%XYYY = 0.0_Kr
   A%YYYY = lambda + 2.0_Kr * mu1

   Write(*, 100, advance = 'no') 'Theta: '
   Read(*,*) theta
   theta = theta * PETSC_PI / 180.0_Kr
   R(1,1) = cos(theta) ; R(1,2) = -sin(theta)
   R(2,1) = sin(theta) ; R(2,2) = cos(theta)

   B = Tens4OS2DTransform(A, R)
   Write(*,200) 'A: ', A
   Write(*,200) 'B: ', B
   Call GenHL_Ortho2D_LambdaMu(lambda, mu1, mu2, theta, A)
   Write(*,200) 'A: ', A

   100 format(A)
   200 format(A, 6(ES12.5, '  '))
   300 format(A, 3(ES12.5, '  '))
   DeAllocate(R)
End Program TestTransform