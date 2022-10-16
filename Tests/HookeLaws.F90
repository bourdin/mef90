Program HookeLaw
#include "petsc/finclude/petsc.h"
   Use m_MEF90
   IMPLICIT NONE

   PetscReal                        :: E, nu
   Type(Tens4OS2D)                  :: HookeLaw2D
   Type(Tens4OS3D)                  :: HookeLaw3D
   PetscBool                        :: flg,mef90
   Character(len=1024)              :: IOBuffer
   PetscInt                         :: ierr

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
      Write(IOBuffer,201) HookeLaw3D
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
      Call MEF90HookeLawIsoEnu2DPlaneStress(HookeLaw2D,E,nu)
      Write(*,*)"2D plane stress"
      Write(IOBuffer,202) HookeLaw2D
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
      Call MEF90HookeLawIsoEnu2DPlaneStrain(HookeLaw2D,E,nu)
      Write(*,*)"2D plane strain"
      Write(IOBuffer,202) HookeLaw2D
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
   Else
      Write(*,*)"3D isotropic"
      Write(IOBuffer,101) HookeLaw3D
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
      Call MEF90HookeLawIsoEnu2DPlaneStress(HookeLaw2D,E,nu)
      Write(*,*)"2D plane stress"
      Write(IOBuffer,102) HookeLaw2D
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
      Call MEF90HookeLawIsoEnu2DPlaneStrain(HookeLaw2D,E,nu)
      Write(*,*)"2D plane strain"
      Write(IOBuffer,102) HookeLaw2D
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
   End If   
100 Format('Hooke''s laws for E=',ES9.2,' nu='ES9.2)
101 Format(20(ES12.5,","),ES12.5,"\n")
102 Format(5(ES12.5,","),ES12.5,"\n")
200 Format('Hooke''s laws for E=',F8.5,' nu='F8.5)
201 Format("[",6(F8.5,"_Kr,"),"   & ! XXXX,XXYY,XXZZ,XXYZ,XXXZ,XXXY\n",     &
           13(" "),5(F8.5,"_Kr,"),"   & !      YYYY,YYZZ,YYYZ,YYXZ,YYXY\n", &
           25(" "),4(F8.5,"_Kr,"),"   & !           ZZZZ ZZYZ,ZZXZ,ZZXY\n", &
           37(" "),3(F8.5,"_Kr,"),"   & !                YXYX,YZXZ,YZXY\n", &
           49(" "),2(F8.5,"_Kr,"),"   & !                     XZXZ,XZXY\n", &
           61(" "),F8.5,"_Kr ]    !                          XYXY\n")
202 Format("[",3(F8.5,"_Kr,"),"   & ! XXXX,XXYY,XXXY\n", &
           13(" "),2(F8.5,"_Kr,"),"   & !      YYYY,YYXY\n",   &
           25(" "),F8.5,"_Kr ]    !           XYXY\n")
   Call MEF90Finalize(ierr)
   Call PetscFinalize()
End Program HookeLaw