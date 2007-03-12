Program Rupt2DA


  Use m_MEF90
  Use m_Rupt_Struct


  Use m_Rupt2DA_Vars	
  Use m_Rupt2DA_U
  Use m_Rupt2DA_V
  Use m_Rupt2DA_Ener
  Use m_Rupt2DA_Proc

  Implicit NONE

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscis.h"

  Logical                                               :: Is_BackTracking
  Integer                                               :: iE
  Real(Kind = Kr)                                       :: Err_Rel_Ener, Tmp_Ener 
  Integer                                               :: TotIter = 0

  Call PetscGetTime(TotalTS, iErr)
  Call PetscGetTime(InitTS, iErr)
  
  Call Init()
  
  Call PetscGetTime(InitTF, iErr)
  If (MEF90_MyRank == 0) Then
     Write(Log_Unit, *) 'Time in Init():                            ', InitTF - InitTS
  End If
  
  Is_BackTracking  = .False.
  
!!! OUTER LOOP: TIME STEPS
  While_TS: Do while (TimeStep <= Size(Params%Load))     
     TotIter = TotIter + 1
     If (MEF90_MyRank == 0) Then
        Write(Log_Unit, 400) TimeStep, Params%Load(TimeStep)
     End If
     
     Call Update_BC_U(TimeStep)
     
!!! UPDATE BC FOR IRREVERSIBILITY
     If (Params%Do_Irrev) Then
        Call Update_BC_V(Geom, Params, MySD_V, Node_db_V, V_Old)
     End If
     Call Assemb_RHS_V(RHS_V, Geom, Params, MySD_V, Elem_db_V, Node_db_V)

!!! INITIALIZE V
     Select Case (Params%Init_V)
     Case (Init_V_ONE)
        If (.NOT. Is_BackTracking) Then
          Call VecSet(V_Dist, 1.0_Kr, iErr)
          Call VecGhostUpdateBegin(V_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)
          Call VecGhostUpdateEnd(V_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)
	End If 	
     Case (Init_V_RND)
        If (.NOT. Is_BackTracking) Then
          Call Init_V_Cracks(Geom, Params, MySD_V, Elem_db_V,  Node_db_V, V_Dist)
	End If  
     Case (Init_V_SPH)
        If (.NOT. Is_BackTracking) Then
          Call Init_V_Spheres(Geom, Params, MySD_V, Elem_db_V,  Node_db_V, V_Dist)
	End If  
     Case (Init_V_PREV)        
	Continue
     Case Default
        If (MEF90_MyRank ==0) Then
           Write(Log_Unit, *) '[ERROR] Unknown Params/Init_V: ', Params%Init_V
        End If
        STOP
     End Select

     Call VecCopy(V_Dist, V_Old, iErr)

!!! INNER LOOP: ALTERNATE MINIMIZATIONS
     ErrV = 2.0_Kr * Params%TolRelax   
     Do_Iter: Do iIter = 1, Params%MaxIterRelax
        If (MEF90_MyRank == 0) Then
           Write(Log_Unit, 410) TimeStep, iIter
        End If

        Call Solve_U(TimeStep)
        Call Solve_V(TimeStep)

!!$!!! Caltech 2006-04 save all intermediate steps
!!$        If (MEF90_MyRank == 0) Then
!!$           Write(Log_Unit, *) 'EXPORTING TMP RESULT'
!!$        End If
!!$        Call Export(TimeStep+iIter+1)
!!$!!!
        Call VecCopy(V_Dist, V_Old, iErr)

        If ( (Mod(iIter, SaveInt) == 0) .OR. (ErrV <= Params%TolRelax) )Then 
           If (MEF90_MyRank ==0) Then
              Write(Log_Unit, *) '    * Saving time step ',  TimeStep
              Call PetscGetTime(TotalTF, iErr)
              Write(Log_Unit, 930) TotalTF - TotalTS
           End If
!           Call VecSet(U_Loc, 1.0_Kr, iErr)
           Call Comp_Bulk_Ener(Bulk_Ener(TimeStep), U_Loc, V_Loc, Geom, Params, MySD_U, MySD_V, Elem_db_U, Elem_db_V, Node_db_U, Node_db_V, Params%Load(TimeStep))
           
           Call Comp_Surf_Ener(Surf_Ener(TimeSTep), V_Loc, Geom, Params, MySD_V, Elem_db_V, Node_db_V)
           Tot_Ener(TimeStep) = Bulk_Ener(TimeStep) + Surf_Ener(TimeSTep)
           If (MEF90_MyRank ==0) Then
              Write(Log_Unit, 910) Bulk_Ener(TimeStep), Surf_Ener(TimeStep), Tot_Ener(TimeStep)
              If (ErrV <= Params%TolRelax) Then
                 Open(File = Ener_Str, Unit =  Ener_Unit, Position = 'append')
                 Write(Ener_Unit, 920) TimeStep, Params%Load(TimeStep), Bulk_Ener(TimeStep), Surf_Ener(TimeStep), Tot_Ener(TimeStep)
              End If
           End If

          Call Export (TimeStep)
!!! Remove this to stop saving -all- intermediate results
!!! This is NOT optimal
!           Call Export2(Size(Params%Load) + TotIter, Params%Load(TimeStep))


           Call PetscLogPrintSummary(PETSC_COMM_WORLD, "petsc_log_summary.log", iErr)

           Is_BackTracking = .FALSE.
           If (Params%Do_BackTrack) Then
              If (MEF90_MyRank == 0) Then  
                 Do iE = 1, TimeStep - 1
                    Tmp_Ener = Bulk_Ener(TimeStep) / Params%Load(TimeStep)**2 * Params%Load(iE)**2 + Surf_Ener(TimeStep)
     	              Err_Rel_Ener = (Tmp_Ener - Tot_Ener(iE)) / ABS(Tot_Ener(TimeStep))	     
                    If (Err_Rel_Ener < -Params%Tol_Ener) Then         
                       Is_BackTracking = .True.
                       EXIT
                    Else
                    End If
                 End Do
                 If (Is_BackTracking) Then
!!! This is a bit silly, but I have already saved the energy if I am not BT''ing
                    If (ErrV > Params%TolRelax) Then
                       Open(File = Ener_Str, Unit =  Ener_Unit, Position = 'append')
                       Write(Ener_Unit, 920) TimeStep, Params%Load(TimeStep), Bulk_Ener(TimeStep), Surf_Ener(TimeStep),            &
					   &     Tot_Ener(TimeStep)
                    End If
                    TimeStep = iE
                    Write(Log_Unit, *) '********** Going back to step ', TimeStep
                    Open(File = Ener_Str, Unit =  Ener_Unit, Position = 'append')
!!$                    Write(Ener_Unit, 920) TimeStep, Params%Load(TimeStep),          &
!!$                         & Bulk_Ener(TimeStep), Surf_Ener(TimeStep),                &
!!$                         & Tot_Ener(TimeStep)
                    Write(Ener_Unit, *) '  '
                    Close(Ener_Unit)
                 End If
              End If
           End If
           Call MPI_BCAST(Is_BackTracking , 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, iErr)
           Call MPI_BCAST(TimeStep, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, iErr)
           If ( (Is_Backtracking) .OR. (ErrV <= Params%TolRelax) ) Then
              EXIT
           End If        
        End If
     End Do Do_Iter

     If (.NOT. Is_BackTracking) Then
        TimeStep = TimeStep + 1
     End If
  End Do While_TS

  Call PetscLogPrintSummary(PETSC_COMM_WORLD, "petsc_log_summary.log", iErr)
  Call Finalize()
  STOP

400 Format('**** Time step: ', I5, ' Load = ', ES12.5)
410 Format('==== Iteration: ', I5, ' /', I5)
901 Format('     Surface:      ', T24, ES12.5)
902 Format('     Bulk:         ', T24, ES12.5)
910 Format('     Energies:     ', T24, 2(ES10.3, ' '), 'Total: ', ES10.3)
920 Format(I4, 4(ES13.5,'  '))
930 Format('     Cumulated time: ', T24, ES12.5)

End Program Rupt2DA


