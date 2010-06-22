Program TestHLOrtho

#include "finclude/petscdef.h"

   Use petsc
   Use m_MEF90
   Use m_VarFrac_Struct
   Implicit NONE

   PetscReal                            :: lambda, mu
   PetscReal                            :: E, nu
   PetscReal                            :: C11, C12, C44
   PetscReal                            :: B, C, Cp
   PetscReal                            :: theta
   Type(Tens4OS2D)                      :: A
   PetscInt                             :: i, N
   
   N=10

   !!!$ Write(*, 100, advance = 'no') '(E, nu): '
   !!!$ Read(*,*)    E, nu
   !!!$ Call GenHL_Iso2D_EnuPlaneStress(E, nu, A) 
   !!!$ Write(*, 200) "Isotropic Hookes law:", A
   !!!$ lambda = E * nu / (1.0_Kr - nu**2) 
   !!!$ mu     = E / (1.0_Kr + nu) * .5_Kr
   !!!$ C11 = 2.0_Kr * lambda *mu / (lambda + 2.0_Kr * mu) + 2.0_Kr * mu
   !!!$ C12 = 2.0_Kr * lambda *mu / (lambda + 2.0_Kr * mu)
   !!!$ C44 = mu
   !!!$ Write(*, 500) "lambda / mu : ", lambda, mu
   !!!$ Write(*, 200) "C11, C12, C44: ", C11, C12, C44
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

   Write(*, 100, advance = 'no') '(C11, C12, C44): '
   Read(*,*) C11, C12, C44
   Write(*, 300) "C11 / C12 / C44: ", C11, C12, C44
   B  = (C11 + 2.0_Kr * C12) / 3.0_Kr 
   C  = C44
   Cp = (C11 - C12) / 2.0_Kr
   Write(*, 300) "B / C / Cp: ", B, C, Cp

   Do i = 0, N
      theta = PETSC_PI * 0.5_Kr * real(i) / real(N) 
      Call GenHL_Cubic2DPlaneStress_Voigt(C11, C12, C44, theta, A)
      Write(*, 400) "Cubic Hookes law:", real(i) * 90 / real(N), A
   End Do

   Write(*, 100, advance = 'no') '(B, C, C''): '
   Read(*,*) B, C, Cp

   Do i = 0, N
      theta = PETSC_PI * 0.5_Kr * real(i) / real(N) 
      Call GenHL_Cubic2DPlaneStress_Zener(B, C, Cp, theta, A)
      Write(*, 400) "Cubic Hookes law:", real(i) * 90 / real(N), A
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
500 format(A, 2(ES12.5, '  '))

End Program TestHLOrtho