#if defined PB_2D
Program Poisson2D
#elif defined PB_3D
Program Poisson3D
#endif

  Use m_MEF90
  Use m_Poisson_Struct
#if defined PB_2D
  Use m_Poisson2D_Vars
  Use m_Poisson2D_Proc
#elif defined PB_3D
  Use m_Poisson3D_Vars
  Use m_Poisson3D_Proc
#endif


  Implicit NONE

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscis.h"

  Integer :: i
  Integer :: iBlk

  Call Init()

  Call PetscGetTime(AssembTS, iErr)
  Call AssembMR_VLV(MR, Geom, Params, Elem_db, Node_db)
  Call PetscGetTime(AssembTF, iErr)
  If (MEF90_MyRank == 0) Then
     Write(CharBuffer,*) 'Total time in AssembMR:                  ',         &
          & AssembTF - AssembTS, '\n'c
     Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
  End If

  Call PetscGetTime(RHSTS, iErr)
  Call AssembRHS_VLV(RHS_U, Geom, Params, Elem_db, Node_db)
  Call PetscGetTime(RHSTF, iErr)
  If (MEF90_MyRank == 0) Then
     Write(CharBuffer,*) 'Total time in AssembRHS_VLV:              ',         &
          & RHSTF- RHSTS, '\n'c
     Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
  End If



  Call PetscGetTime(SolveTS, iErr)
!  Call MatView(MR, PETSC_VIEWER_STDOUT_WORLD, iErr)
  Call KSPSolve(KSP_U, RHS_U, U_Dist, iErr)
  Call PetscGetTime(SolveTF, iErr)
  Call KSPGetIterationNumber(KSP_U, NbIter, iErr)
  If (MEF90_MyRank == 0) Then
     Write(CharBuffer,*) 'Total time in KSP_Solve:                 ',         &
          & SolveTF- SolveTS, '\n'c
     Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
     Write(CharBuffer,*) 'Number of iterations:                    ',         &
          & NbIter, '\n'c
     Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
  End If

  Call SaveToEXO()
!  Call SaveToEnsight()
!  Call VecView(U_Dist, PETSC_VIEWER_STDOUT_WORLD, iErr)
  Call Finalize()



#if defined PB_2D
End Program Poisson2D
#elif defined PB_3D
End Program Poisson3D
#endif

