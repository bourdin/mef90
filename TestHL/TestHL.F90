Program TestHL

#include "finclude/petscdef.h"

Use petsc
Use m_MEF90
Use m_VarFrac_Struct
Implicit NONE

PetscReal                            :: lambda_A, mu_A, E_A, nu_A
PetscReal                            :: lambda_B, mu_B, E_B, nu_B
PetscReal                            :: alpha, theta
PetscInt                             :: nlayers, i
Type(Tens4OS2D)                      :: A, B, Astar
Type(Vect2D)                         :: k
Type(MatS2D)                         :: xi

Write(*, 100, advance = 'no') '(E_A, nu_A): '
Read(*,*)    E_A, nu_A
lambda_A = E_A * nu_A / (1.0_Kr - nu_A**2) 
mu_A     = E_A / (1.0_Kr + nu_A) * .5_Kr
Call GenHL_Iso2D_EnuPlaneStress(E_A, nu_A, Astar) 

Write(*, 100, advance = 'no') '(E_B, nu_B): '
Read(*,*)    E_B, nu_B
lambda_B = E_B * nu_B / (1.0_Kr - nu_B**2) 
mu_B     = E_B / (1.0_Kr + nu_B) * .5_Kr
!Call GenHL_Iso2D_EnuPlaneStress(E_B, nu_B, B) 

Write(*, 100, advance = 'no') 'Number of layers: '
Read(*,*) nlayers

Do i = 1, nlayers
   Write(*, 100, advance = 'no') '   Lamination angle in degrees: '
   Read(*,*) alpha
   alpha = alpha/180_Kr * PETSC_PI
   k%X = cos(alpha); k%Y = sin(alpha)

   Write(*, 100, advance = 'no') '   Volume fraction: '
   Read(*,*) theta 

   Call GenHL_Laminate_LambdaMu(Astar, k, lambda_B, mu_B, theta, Astar)
   Write(*,200) 'Astar:      ', Astar
End Do




100 format(A)
200 format(A, 6(ES12.5, '  '))
300 format(A, 3(ES12.5, '  '))

End Program TestHL