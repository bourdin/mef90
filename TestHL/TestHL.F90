Program TestHL

#include "finclude/petscdef.h"

Use petsc
Use m_MEF90
Use m_VarFrac_Struct
Implicit NONE

PetscReal                :: lambda_A, mu_A, E_A, nu_A
PetscReal                :: lambda_B, mu_B, E_B, nu_B
PetscReal                :: alpha, theta
Type(Tens4OS2D)          :: A, B, Astar
Type(Vect2D)             :: k
Type(MatS2D)             :: xi

Write(*, 100, advance = 'no') '(E_A, nu_A): '
Read(*,*)    E_A, nu_A
lambda_A = E_A * nu_A / (1.0_Kr - nu_A**2) 
mu_A     = E_A / (1.0_Kr + nu_A) * .5_Kr
!Write(*,*) '(lambda_A, mu_A): ', lambda_A, mu_A
Call GenHL_Iso2D_EnuPlaneStress(E_A, nu_A, A) 
!Write(*,200) 'A: ', A

!!!$ !!! Make sure that the Hooke law inversion works
!!!$ B = invert(A)
!!!$ Write(*,200) 'A^{-1}: ', B
!!!$
!!!$ xi = 0.0_Kr ; xi%XX = 1.0_Kr
!!!$ Write(*,300) 'xi:       ', xi
!!!$ Write(*,300) 'A (B xi): ', A * (B*xi)
!!!$ Write(*,300) 'B (A xi): ', B * (A*xi)
!!!$ 
!!!$ xi = 0.0_Kr ; xi%XY = 0.5_Kr
!!!$ Write(*,300) 'xi:       ', xi
!!!$ Write(*,300) 'A (B xi): ', A * (B*xi)
!!!$ Write(*,300) 'B (A xi): ', B * (A*xi)
!!!$ 
!!!$ xi = 0.0_Kr ; xi%YY = 1.0_Kr
!!!$ Write(*,300) 'xi:       ', xi
!!!$ Write(*,300) 'A (B xi): ', A * (B*xi)
!!!$ Write(*,300) 'B (A xi): ', B * (A*xi)
!!!$ STOP
!!!$ !!! It does work
 
Write(*, 100, advance = 'no') '(E_B, nu_B): '
Read(*,*)    E_B, nu_B
lambda_B = E_B * nu_B / (1.0_Kr - nu_B**2) 
mu_B     = E_B / (1.0_Kr + nu_B) * .5_Kr
!Write(*,*) '(lambda_B, mu_B): ', lambda_B, mu_B
Call GenHL_Iso2D_EnuPlaneStress(E_B, nu_B, B) 
!Write(*,200) 'B: ', B

Write(*, 100, advance = 'no') 'Lamination angle in degrees: '
Read(*,*) alpha
alpha = alpha/180_Kr * PETSC_PI
k%X = cos(alpha); k%Y = sin(alpha)

Write(*, 100, advance = 'no') 'Volume fraction: '
Read(*,*) theta 



!!!$ xi = 0.0_Kr ; xi%XX = 1.0_Kr
!!!$ Write(*,300) 'xi:       ', xi
!!!$ Write(*,*) 'fB(k)xi:xi: ', StrainProjectionComponent2D_LambdaMu(xi, k, lambda_B, mu_B)
!!!$ 
!!!$ xi = 0.0_Kr ; xi%YY = 1.0_Kr
!!!$ Write(*,300) 'xi:       ', xi
!!!$ Write(*,*) 'fB(k)xi:xi: ', StrainProjectionComponent2D_LambdaMu(xi, k, lambda_B, mu_B)
!!!$ 
!!!$ xi = 0.0_Kr ; xi%XY = 0.5_Kr
!!!$ Write(*,300) 'xi:       ', xi
!!!$ Write(*,*) 'fB(k)xi:xi: ', StrainProjectionComponent2D_LambdaMu(xi, k, lambda_B, mu_B)
!!!$ 
Astar = StrainProjection2D_LambdaMu(k, lambda_B, mu_B)

write(*,200) 'Projection: ', Astar

Call GenHL_Laminate_LambdaMu(A, k, lambda_B, mu_B, theta, Astar)
Write(*,200) 'Astar:      ', Astar


!!! Make sure that the Hooke law inversion works
B = invert(Astar)
Write(*,200) 'A*^{-1}: ', B

xi = 0.0_Kr ; xi%XX = 1.0_Kr
Write(*,300) 'xi:       ', xi
Write(*,300) 'Astar (B xi): ', Astar * (B*xi)
Write(*,300) 'B (Astar xi): ', B * (Astar*xi)

xi = 0.0_Kr ; xi%XY = 0.5_Kr
Write(*,300) 'xi:       ', xi
Write(*,300) 'Astar (B xi): ', Astar * (B*xi)
Write(*,300) 'B (Astar xi): ', B * (Astar*xi)

xi = 0.0_Kr ; xi%YY = 1.0_Kr
Write(*,300) 'xi:       ', xi
Write(*,300) 'Astar (B xi): ', Astar * (B*xi)
Write(*,300) 'B (Astar xi): ', B * (Astar*xi)
STOP
!!! It does work

100 format(A)
200 format(A, 6(ES12.5, '  '))
300 format(A, 3(ES12.5, '  '))

End Program TestHL