Program TestHLOrtho

#include "finclude/petscdef.h"

   Use petsc
   Use m_MEF90
   Use m_VarFrac_Struct
   Implicit NONE

   PetscReal                            :: lambda, mu, E, nu, G, mu1, mu2
   PetscReal                            :: theta
   Type(Tens4OS2D)                      :: A
   PetscInt                             :: i, N
   
   N=11

   !!!$ Write(*, 100, advance = 'no') '(E, nu): '
   !!!$ Read(*,*)    E, nu
   !!!$ Call GenHL_Iso2D_EnuPlaneStress(E, nu, A) 
   !!!$ Write(*, 200) "Isotropic Hookes law:", A
   !!!$ lambda = E * nu / (1.0_Kr - nu**2) 
   !!!$ mu     = E / (1.0_Kr + nu) * .5_Kr
   !!!$ G      = mu
   !!!$ Write(*, 300) "lambda / mu / G: ", lambda, mu, G

   !!!$ !Write(*, 100, advance = 'no') '(lambda, mu1, mu2): '
   !!!$ !Read(*,*) lambda, mu1, mu2
   !!!$ mu1 = mu
   !!!$ mu2 = mu
   !!!$ Write(*, 300) "lambda / mu1 / mu2: ", lambda, mu1, mu2
   !!!$ 
   !!!$ Do i = 0, N
   !!!$    theta = PETSC_PI * 0.5_Kr * real(i) / real(N) 
   !!!$    Call GenHL_Ortho2D_LambdaMu(lambda, mu1, mu2, theta, A)
   !!!$    Write(*, 400) "Orthotropic Hookes law:", real(i) * 90 / real(N), A
   !!!$ End Do

   Write(*, 100, advance = 'no') '(E,G,mu): '
   Read(*,*) E,G,mu
   Write(*, 300) "E / G / mu: ", E, G, mu

   Do i = 0, N
      theta = PETSC_PI * 0.5_Kr * real(i) / real(N) 
      Call GenHL_Ortho2D_EGmu(E,G,mu, theta, A)
      Write(*, 400) "Orthotropic Hookes law:", real(i) * 90 / real(N), A
   End Do

   !!!$ !Write(*, 100, advance = 'no') '(E,G,mu): '
   !!!$ !Read(*,*) E,G,nu
   !!!$ Write(*, 300) "E / G / nu: ", E, G, nu
   !!!$ 
   !!!$ Do i = 0, N
   !!!$    theta = PETSC_PI * 0.5_Kr * real(i) / real(N) 
   !!!$    Call GenHL_Ortho2D_EGnu(E,G,nu, theta, A)
   !!!$    Write(*, 400) "Orthotropic Hookes law:", real(i) * 90 / real(N), A
   !!!$ End Do

100 format(A)
200 format(A, T30, 6(ES12.5, '  '))
300 format(A, 3(ES12.5, '  '))
400 format(A, F5.1, T30, 6(ES12.5, '  '))

End Program TestHLOrtho