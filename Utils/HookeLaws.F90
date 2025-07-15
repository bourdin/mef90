program HookeLaw
#include "petsc/finclude/petsc.h"
   use m_MEF90
   use petsc
   implicit none(type, external)

   PetscReal                        :: E, nu
   type(Tens4OS2D)                  :: HookeLaw2D
   type(Tens4OS3D)                  :: HookeLaw3D
   PetscBool                        :: flg, mef90
   character(len=1024)              :: IOBuffer
   PetscErrorCode                   :: ierr

   PetscCallA(PetscInitialize(ierr))
   PetscCallA(MEF90Initialize(PETSC_COMM_WORLD, ierr))

   E = 1.0_kr
   PetscCallA(PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-E', E, flg, ierr))
   nu = .3_kr
   PetscCallA(PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-nu', nu, flg, ierr))
   mef90 = PETSC_FALSE
   PetscCallA(PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-mef90', mef90, flg, ierr))

   write (*, 100) E, nu
   PetscCallA(MEF90HookeLawIsoENu3D(HookeLaw3D, E, nu))

   if (mef90) then
      write (*, *) "3D isotropic"
      write (IOBuffer, 201) HookeLaw3D
      PetscCallA(PetscPrintf(PETSC_COMM_WORLD, IOBuffer, ierr))
      PetscCallA(MEF90HookeLawIsoEnu2DPlaneStress(HookeLaw2D, E, nu))
      write (*, *) "2D plane stress"
      write (IOBuffer, 202) HookeLaw2D
      PetscCallA(PetscPrintf(PETSC_COMM_WORLD, IOBuffer, ierr))
      PetscCallA(MEF90HookeLawIsoEnu2DPlaneStrain(HookeLaw2D, E, nu))
      write (*, *) "2D plane strain"
      write (IOBuffer, 202) HookeLaw2D
      PetscCallA(PetscPrintf(PETSC_COMM_WORLD, IOBuffer, ierr))
   else
      write (*, *) "3D isotropic"
      write (IOBuffer, 101) HookeLaw3D
      PetscCallA(PetscPrintf(PETSC_COMM_WORLD, IOBuffer, ierr))
      PetscCallA(MEF90HookeLawIsoEnu2DPlaneStress(HookeLaw2D, E, nu))
      write (*, *) "2D plane stress"
      write (IOBuffer, 102) HookeLaw2D
      PetscCallA(PetscPrintf(PETSC_COMM_WORLD, IOBuffer, ierr))
      PetscCallA(MEF90HookeLawIsoEnu2DPlaneStrain(HookeLaw2D, E, nu))
      write (*, *) "2D plane strain"
      write (IOBuffer, 102) HookeLaw2D
      PetscCallA(PetscPrintf(PETSC_COMM_WORLD, IOBuffer, ierr))
   end if
100 format('Hooke''s laws for E=', ES9.2, ' nu='ES9.2)
101 format(20(ES12.5, ","), ES12.5, "\n")
102 format(5(ES12.5, ","), ES12.5, "\n")
201 format("[", 6(F8.5, "_Kr,"), "   & ! XXXX,XXYY,XXZZ,XXYZ,XXXZ,XXXY\n", &
          13(" "), 5(F8.5, "_Kr,"), "   & !      YYYY,YYZZ,YYYZ,YYXZ,YYXY\n", &
          25(" "), 4(F8.5, "_Kr,"), "   & !           ZZZZ ZZYZ,ZZXZ,ZZXY\n", &
          37(" "), 3(F8.5, "_Kr,"), "   & !                YXYX,YZXZ,YZXY\n", &
          49(" "), 2(F8.5, "_Kr,"), "   & !                     XZXZ,XZXY\n", &
          61(" "), F8.5, "_Kr ]    !                          XYXY\n")
202 format("[", 3(F8.5, "_Kr,"), "   & ! XXXX,XXYY,XXXY\n", &
          13(" "), 2(F8.5, "_Kr,"), "   & !      YYYY,YYXY\n", &
          25(" "), F8.5, "_Kr ]    !           XYXY\n")
   PetscCallA(MEF90Finalize(ierr))
   PetscCallA(PetscFinalize(ierr))
end program HookeLaw
