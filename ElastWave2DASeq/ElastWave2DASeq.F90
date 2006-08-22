Program Elast2DA

  Use m_MEF90
  Use m_Rupt_Struct


  Use m_Elast2DA_Vars
  Use m_Elast2DA_Proc
  Implicit NONE

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscao.h"


  Integer               :: iE, iS, NbIter
  PetscScalar           :: ener
  PetscScalar           :: value

  PetscReal             :: One = 1.0_Kr
  PetscReal             :: Zero = 0.0_Kr

  Real(Kind = Kr), Dimension(:), Pointer       :: BC_Ptr
  Real(Kind = Kr), Dimension(:), Pointer       :: EXO_Ptr

  Mat                   :: Loc_Mat

  Call PetscGetTime(TotalTS, iErr)
  Call PetscGetTime(InitTS, iErr)

  Call Init()

  Call PetscGetTime(InitTF, iErr)

  Call PetscGetTime(AssembTS, iErr)

  Call Assemb_Mat_Elast2DA(MR, Geom, Params, Elem_db, Node_db)

  Call PetscGetTime(AssembTF, iErr)

  Call Init_KSP()

  Write(CharBuffer,*) 'Number of timesteps:                        ',      &
       & Size(Params%Load), '\n'c
  Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)

  Do TimeStep = 1, Size(Params%Load)
     Write(CharBuffer,*) '=== processing time step                 ',         &
          & TimeStep, '\n'c
     Call PetscPrintf(PETSC_COMM_WORLD, CharBuffer, iErr)

     !!! EXO -> BC_Ptr -> BC_Master
     Call VecGetArrayF90(BC, BC_Ptr, iErr)
     Call Read_EXO_Result_Nodes(Geom, 4, TimeStep, EXO_Ptr, 1)
     BC_Ptr = EXO_Ptr
     Call VecRestoreArrayF90(BC, BC_Ptr, iErr)
     DeAllocate(EXO_Ptr)

     Call VecGetArrayF90(F, BC_Ptr, iErr)
     Call Read_EXO_Result_Nodes(Geom, 7, TimeStep, EXO_Ptr, 1)
     BC_Ptr = EXO_Ptr
     Call VecRestoreArrayF90(F, BC_Ptr, iErr)
     DeAllocate(EXO_Ptr)

     !!! BC_Loc -> RHS
     Call PetscGetTime(RHSTS, iErr)
     Call Assemb_RHS_Elast2DA(RHS, Geom, Params, Elem_db, Node_db, BC, F)
     Call PetscGetTime(RHSTF, iErr)
     Write(CharBuffer,*) 'Total time in Assemb_RHS_Elast:          ',      &
          & RHSTF - RHSTS, '\n'c
     Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
     
     Call PetscGetTime(SolveTS, iErr)
     Call KSPSolve(KSP_MR, RHS, Sol, iErr)

     Call PetscGetTime(SolveTF, iErr)
     Write(CharBuffer,*) 'Total time in KSP_Solve:                 ',      &
          & SolveTF- SolveTS, '\n'c
     Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)

     Call KSPGetIterationNumber(KSP_MR, NbIter, iErr)
     Write(CharBuffer,*) 'Number of iterations:                    ',      &
          & NbIter, '\n'c
     Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)

!!! Saving the results
     Call Export()
  End Do

  Call PetscGetTime(TotalTF, iErr)
  Write(CharBuffer,*) 'Total time                                  ',      &
       & TotalTF - TotalTS, '\n'c
  Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)

  Call Finalize()

  100 Format(A)
End Program Elast2DA
