Program Elast

  Use m_MEF90
  Use m_Rupt_Struct


#ifdef PB_2D
  Use m_Elast2D_Vars
  Use m_Elast2D_Proc
  Use m_Elast2D_Debug
#else
  Use m_Elast3D_Vars
  Use m_Elast3D_Proc
  Use m_Elast3D_Debug
#endif
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
  Call VecSet(SOL_Dist, 1.0_Kr, iErr)


!!$  Call GenVect(Geom, MySD, Node_db, Elem_db, SOL_Dist)
!!$  Call VecGhostUpdateBegin(SOL_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)
!!$  Call VecGhostUpdateEnd(SOL_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)
!!$  Call Calc_Ener(SOL_Loc, Geom, Params, Elem_db, Node_db, MySD,            &
!!$       & 1.0_Kr, Ener_Elast) 
!!$  TimeStep = 1
!!$  Call Export()
!!$
!!$     STOP


#ifdef MEF90_TIMING
  Call PetscGetTime(InitTF, iErr)
  If (MEF90_MyRank == 0) Then
     Write(CharBuffer, *) 'Time in Init():                            ',      &
          & InitTF - InitTS, '\n'c
     Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
  End If
#endif


  Call PetscGetTime(AssembTS, iErr)

  Call Assemb_Mat_Elast(MR, Geom, Params, Elem_db, Node_db)

!  If (MEF90_MyRank == 0) Then
!     Call MatGetLocalMat(MR, MAT_INITIAL_MATRIX, Loc_Mat, iErr)
!  End If
!!$  Print*, 'BC_Ptr: ', Associated(BC_Ptr), Size(BC_Ptr)
!!$  Allocate(BC_Ptr(Size(Node_db)))
!!$  Print*, 'BC_Ptr: ', Associated(BC_Ptr), Size(BC_Ptr)
!!$  BC_Ptr = 0.0_Kr
!!$  Call Assemb_RHS_Elast(RHS, Geom, Params, Elem_db, Node_db, BC_Dist)

!  Call MEF90_FINALIZE()
!  STOP
!  Call MatView(MR, PETSC_VIEWER_DRAW_WORLD, iErr)

#ifdef MEF90_TIMING
  Call PetscGetTime(AssembTF, iErr)
  If (MEF90_MyRank == 0) Then
     Write(CharBuffer,*) 'Total time in Assemb_Mat_Elast:        :    ',      &
          & AssembTF - AssembTS, '\n'c
     Call PetscPrintf(PETSC_COMM_WORLD, CharBuffer, iErr)
  End If
#endif


  Call Init_KSP()

  If (MEF90_MyRank ==0) Then
     Write(CharBuffer,*) 'Number of timesteps:                        ',      &
          & Size(Params%Load), '\n'c
     Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
  End If

  Do TimeStep = 1, Size(Params%Load)
     Write(CharBuffer,*) '=== processing time step                 ',         &
          & TimeStep, '\n'c
     Call PetscPrintf(PETSC_COMM_WORLD, CharBuffer, iErr)

     !!! EXO -> BC_Ptr -> BC_Master
     If (MEF90_MyRank == 0) Then
        Call VecGetArrayF90(BC_Master, BC_Ptr, iErr)
        Call Read_EXO_Result_Nodes(Geom, 2, TimeStep, EXO_Ptr, Geom%Num_Dim)
        BC_Ptr = EXO_Ptr
        Call VecRestoreArrayF90(BC_Master, BC_Ptr, iErr)
        DeAllocate(EXO_Ptr)
        Call VecGetArrayF90(F_Master, BC_Ptr, iErr)
        Call Read_EXO_Result_Nodes(Geom, 5, TimeStep, EXO_Ptr, Geom%Num_Dim)
        BC_Ptr = EXO_Ptr
        Call VecRestoreArrayF90(F_Master, BC_Ptr, iErr)
        DeAllocate(EXO_Ptr)
        Call VecGetArrayF90(Temp_Master, BC_Ptr, iErr)
        Call Read_EXO_Result_Nodes(Geom, 8, TimeStep, EXO_Ptr,1)
        BC_Ptr = EXO_Ptr
        Call VecRestoreArrayF90(Temp_Master, BC_Ptr, iErr)
        DeAllocate(EXO_Ptr)
     EndIf

     !!! XXX_Master -> XXX_Dist
     Call VecScatterBegin(BC_Master, BC_Dist, INSERT_VALUES, SCATTER_REVERSE, &
          & MySD_Vect%ToMaster, iErr)
     Call VecScatterEnd(BC_Master, BC_Dist, INSERT_VALUES, SCATTER_REVERSE,   &
          & MySD_Vect%ToMaster, iErr)

     Call VecScatterBegin(F_Master, F_Dist, INSERT_VALUES, SCATTER_REVERSE,   &
          & MySD_Vect%ToMaster, iErr)
     Call VecScatterEnd(F_Master, F_Dist, INSERT_VALUES, SCATTER_REVERSE,     &
          & MySD_Vect%ToMaster, iErr)

     Call VecScatterBegin(Temp_Master, Temp_Dist, INSERT_VALUES,              &
     	  & SCATTER_REVERSE, MySD_Scal%ToMaster, iErr)
     Call VecScatterEnd(Temp_Master, Temp_Dist, INSERT_VALUES,                &
     	  & SCATTER_REVERSE, MySD_Scal%ToMaster, iErr)

     !!! XXX_Dist -> XXX_Loc
     Call VecGhostUpdateBegin(BC_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)
     Call VecGhostUpdateEnd(BC_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)

     Call VecGhostUpdateBegin(F_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)
     Call VecGhostUpdateEnd(F_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)

     Call VecGhostUpdateBegin(Temp_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)
     Call VecGhostUpdateEnd(Temp_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)

     !!! BC_Loc -> RHS
     Call PetscGetTime(RHSTS, iErr)
     Call Assemb_RHS_Elast(RHS, Geom, Params, Elem_db, Node_db, Elem_Scal,    &
          & Node_Scal, BC_Loc, F_Loc, Temp_Loc)
     Call PetscGetTime(RHSTF, iErr)
     If (MEF90_MyRank == 0) Then
        Write(CharBuffer,*) 'Total time in Assemb_RHS_Elast:          ',      &
             & RHSTF - RHSTS, '\n'c
        Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
     End If
!     Call VecView(RHS, PETSC_VIEWER_DRAW_WORLD, iErr)

     Call PetscGetTime(SolveTS, iErr)
     Call KSPSolve(KSP_MR, RHS, Sol_Dist, iErr)

!     Call VecView(SOL_Dist, PETSC_VIEWER_DRAW_WORLD, iErr)
#ifdef MEF90_TIMING
     Call PetscGetTime(SolveTF, iErr)
     If (MEF90_MyRank == 0) Then
        Write(CharBuffer,*) 'Total time in KSP_Solve:                 ',      &
             & SolveTF- SolveTS, '\n'c
        Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
     End If
#endif

     If (MEF90_MyRank ==0) Then
        Call KSPGetIterationNumber(KSP_MR, NbIter, iErr)
        Write(CharBuffer,*) 'Number of iterations:                    ',      &
             & NbIter, '\n'c
        Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
     End If

     Call VecGhostUpdateBegin(SOL_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)
     Call VecGhostUpdateEnd(SOL_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)
     Call Calc_Ener(SOL_Loc, Geom, Params, Elem_db, Node_db, MySD_Vect,       &
          & F_Loc, Ener_Elast) 
     Write(CharBuffer, *) '=== Elastic Energy: ', Ener_Elast, '\n'c
     Call PETScPrintf(PETSC_COMM_WORLD, CharBuffer, iErr)

     Call PetscGetTime(ExportTS, iErr)
     Call Export()

#ifdef MEF90_TIMING
     Call PetscGetTime(ExportTF, iErr)
     If (MEF90_MyRank ==0) Then
        Write(CharBuffer,*) 'Total time in Export:                    ',      &
             & ExportTF - ExportTS, '\n'c
        Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
     End If
#endif
  End Do

  Call PetscGetTime(TotalTF, iErr)
  If (MEF90_MyRank ==0) Then
     Write(CharBuffer,*) 'Total time                                  ',      &
          & TotalTF - TotalTS, '\n'c
     Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)
  End If

  Call Finalize()

  100 Format(A)

End Program Elast
