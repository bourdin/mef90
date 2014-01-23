Program HokeLaw
#include "finclude/petscdef.h"
   Use m_MEF90

   PetscReal               :: E, nu
   Type(Tens4OS2D)         :: HookeLaw2D
   Type(Tens4OS3D)         :: HookeLaw3D
   PetscBool               :: flg,mef90
   
   Call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
   Call MEF90Initialize(ierr)
   
   E = 1.0_Kr
   Call PetscOptionsGetReal(PETSC_NULL_CHARACTER,'-E',E,flg,ierr);CHKERRQ(ierr);
   nu = .3_Kr
   Call PetscOptionsGetReal(PETSC_NULL_CHARACTER,'-nu',nu,flg,ierr);CHKERRQ(ierr);
   mef90 = PETSC_FALSE
   Call PetscOptionsGetBool(PETSC_NULL_CHARACTER,'-mef90',mef90,flg,ierr);CHKERRQ(ierr);
   
   Write(*,100) E,nu
   Call MEF90HookeLawIsoENu3D(HookeLaw3D,E,nu)
   
   If (mef90) Then
      Write(*,*)"3D isotropic"
      Write(*,201) HookeLaw3D
      Call MEF90HookeLawIsoEnu2DPlaneStress(HookeLaw2D,E,nu)
      Write(*,*)"2D plane stress"
      Write(*,202) HookeLaw2D
      Call MEF90HookeLawIsoEnu2DPlaneStrain(HookeLaw2D,E,nu)
      Write(*,*)"2D plane strain"
      Write(*,202) HookeLaw2D
   Else
      Write(*,*)"3D isotropic"
      Write(*,101) HookeLaw3D
      Call MEF90HookeLawIsoEnu2DPlaneStress(HookeLaw2D,E,nu)
      Write(*,*)"2D plane stress"
      Write(*,102) HookeLaw2D
      Call MEF90HookeLawIsoEnu2DPlaneStrain(HookeLaw2D,E,nu)
      Write(*,*)"2D plane strain"
      Write(*,102) HookeLaw2D
   End If   
100 Format('Hooke''s laws for E=',ES9.2,' nu='ES9.2)
101 Format(20(ES12.5,","),ES12.5)
102 Format(5(ES12.5,","),ES12.5)
200 Format('Hooke''s laws for E=',F8.5,' nu='F8.5)
201 Format("[",20(F8.5,"_Kr,"),F8.5,"_Kr ]")
202 Format("[",5(F8.5,"_Kr,"),F8.5,"_Kr ]")
   Call MEF90Finalize(ierr)
   Call PetscFinalize()
End Program HokeLaw