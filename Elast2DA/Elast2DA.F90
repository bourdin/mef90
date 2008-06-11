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


   Integer                                      :: iE, iS, NbIter
   PetscScalar                                  :: ener
   PetscScalar                                  :: value
                                                
   PetscReal                                    :: One = 1.0_Kr
   PetscReal                                    :: Zero = 0.0_Kr
   
   Real(Kind = Kr), Dimension(:), Pointer       :: IO_Ptr
   Real(Kind = Kr), Dimension(:), Pointer       :: EXO_Ptr
   
   Mat                                          :: Loc_Mat
   
   Call Init()
   Call PetscGetTime(TotalTS, iErr)
   

!!!    Timestep = 2
!!!    If (MEF90_MyRank == 0) Then
!!!       Call VecGetArrayF90(IO_N, IO_Ptr, iErr)
!!!       Call Read_EXO_Result_Nodes(Geom, 4, TimeStep, EXO_Ptr, 1)
!!!       IO_Ptr = EXO_Ptr
!!!       Call VecRestoreArrayF90(IO_N, IO_Ptr, iErr)
!!!       DeAllocate(EXO_Ptr)
!!!    End If
!!!    
!!!    Call VecView(IO_N, PETSC_VIEWER_STDOUT_WORLD, iErr)
!!!    Call VecScatterBegin(Layout%ToIO_N, IO_N, BC, INSERT_VALUES, SCATTER_REVERSE, iErr) 
!!!    Call VecScatterEnd  (Layout%ToIO_N, IO_N, BC, INSERT_VALUES, SCATTER_REVERSE, iErr) 
!!!    Call PetscPrintf(PETSC_COMM_WORLD, "BC ===========================\n"c, iErr)
!!!    Call VecView(BC, PETSC_VIEWER_STDOUT_WORLD, iErr)
!!!    Call PetscPrintf(PETSC_COMM_WORLD, "BC ===========================\n"c, iErr)
!!! 
!!!    Call VecGhostUpdateBegin(BC, INSERT_VALUES, SCATTER_FORWARD, iErr)
!!!    Call VecGhostUpdateEnd  (BC, INSERT_VALUES, SCATTER_FORWARD, iErr)
!!!    Call PetscPrintf(PETSC_COMM_WORLD, "BC_Loc ===========================\n"c, iErr)
!!!    Call VecView(BC_Local, PETSC_VIEWER_STDOUT_SELF, iErr)
!!!    Call PetscPrintf(PETSC_COMM_WORLD, "BC_Loc ===========================\n"c, iErr)
!!! 
!!! 
!!!    Call PetscFinalize(iErr)
!!!    STOP
   
   Call Show_EXO_Geom_Info(MyGeom, 200+MEF90_MyRank)
   Do iE = 1, Size(Elem_db)
      Write(200+MEF90_MyRank, *) iE, Elem_db(iE)%ID_DoF
   End Do
   
   Call PetscPrintf(PETSC_COMM_WORLD, "Assembling the matrix\n"c, ierr)
   Call Assemb_Mat_Elast2DA(MR, MyGeom, Params, Elem_db, Node_db)
!   Call MatView(MR, PETSC_VIEWER_STDOUT_WORLD, iErr)
   
!   Call PetscFinalize(iErr)
!   STOP
   Call Init_KSP()
   

   Do TimeStep = 1, Size(Params%Load)
      Write(CharBuffer,*) '=== processing time step                 ', TimeStep, '\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, CharBuffer, iErr)
      
      !!! Read BC
      Call Read_EXO_Result_Nodes(Geom, Layout, 4, TimeStep, BC)
      Call VecGhostUpdateBegin(BC, INSERT_VALUES, SCATTER_FORWARD, iErr)
      Call VecGhostUpdateEnd  (BC, INSERT_VALUES, SCATTER_FORWARD, iErr)

!      Call PetscPrintf(PETSC_COMM_WORLD, "======== BC ========\n"c, iErr)
!      Call VecView(BC, PETSC_VIEWER_STDOUT_WORLD, iErr)

      !!! Read F
      Call Read_EXO_Result_Nodes(Geom, Layout, 7, TimeStep, F)
      Call VecGhostUpdateBegin(F, INSERT_VALUES, SCATTER_FORWARD, iErr)
      Call VecGhostUpdateEnd  (F, INSERT_VALUES, SCATTER_FORWARD, iErr)
      
      !!! BC_Loc -> RHS
      Call Assemb_RHS_Elast2DA(RHS, MyGeom, Params, Elem_db, Node_db, BC_Local, F_Local)

      Call PetscGetTime(SolveTS, iErr)
      Call KSPSolve(KSP_MR, RHS, Sol, iErr)
!      Call MatView(MR, PETSC_VIEWER_STDOUT_WORLD, iErr)
!      Call VecView(RHS, PETSC_VIEWER_STDOUT_WORLD, iErr)
      
      Call PetscGetTime(SolveTF, iErr)
      Write(CharBuffer,*) 'Total time in KSP_Solve:                 ', SolveTF- SolveTS, '\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, CharBuffer, iErr)
      
      Call KSPGetIterationNumber(KSP_MR, NbIter, iErr)
      Write(CharBuffer,*) 'Number of iterations:                    ', NbIter, '\n'c
      Call PetscPrintf(PETSC_COMM_WORLD, CharBuffer, iErr)
      
      !!! Saving the results
      Call Export()
   End Do

   Call PetscGetTime(TotalTF, iErr)
   Write(CharBuffer,*) 'Total time                               ', TotalTF - TotalTS, '\n'c
   Call PetscPrintf(PETSC_COMM_WORLD, CharBuffer, iErr)

   Call MEF90_Finalize()
   STOP

   Call Finalize()
   
100 Format(A)
End Program Elast2DA
