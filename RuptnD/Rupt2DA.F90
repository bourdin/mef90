#if defined PB_2D
Program  Rupt2D
#elif defined PB_3D
Program Rupt3D
#else 
Program Rupt2DA
#endif

  Use m_MEF90
  Use m_Rupt_Struct

#if defined PB_2D
  Use m_Rupt2D_Vars
  Use m_Rupt2D_U
  Use m_Rupt2D_V
  Use m_Rupt2D_Ener
  Use m_Rupt2D_Proc	
#elif defined PB_3D
  Use m_Rupt3D_Vars
  Use m_Rupt3D_U
  Use m_Rupt3D_V
  Use m_Rupt3D_Ener
  Use m_Rupt3D_Proc
#else 
  Use m_Rupt2DA_Vars	
  Use m_Rupt2DA_U
  Use m_Rupt2DA_V
  Use m_Rupt2DA_Ener
  Use m_Rupt2DA_Proc
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

  Integer                                               :: VMinPos, VMaxPos
  Logical                                               :: Is_BackTracking
  Integer                                               :: NbIterKSP, iE
  Real(Kind = Kr)                                       :: Err_Rel_Ener, Tmp_Ener 
  PetscLogDouble                                        :: MyFlops, TotalFlops


  Call PetscGetTime(TotalTS, iErr)
  Call PetscGetTime(InitTS, iErr)
  
  Call Init()

  Call PetscGetTime(InitTF, iErr)
  If (MEF90_MyRank == 0) Then
     Write(Log_Unit, *) 'Time in Init():                            ',      &
          & InitTF - InitTS, '\n'c
  End If

  Is_BackTracking  = .False.
  
  While_TS: Do while (TimeStep <= Size(Params%Load))
     
     If (MEF90_MyRank == 0) Then
        Write(Log_Unit, 400) TimeStep, Params%Load(TimeStep)
     End If

     Call Update_BC_U(TimeStep)

     Call Update_F(TimeStep)

     Call Update_Temp(TimeStep)

     Call Assemb_RHS_U(RHS_U, BCU_loc, Geom, Params, MySD_U, Elem_db_U,      &
          & Node_db_U, MySD_V, Elem_db_V, Node_db_V, V_Loc, F_Loc, Temp_Loc)

     Select Case (Params%Init_U)
     Case (Init_U_ZERO)
        Call VecSet(0.0_Kr, U_Dist, iErr)
        Call VecGhostUpdateBegin(U_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)
        Call VecGhostUpdateEnd(U_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)
     Case (Init_U_PREV)
        Continue
     Case Default
        If (MEF90_MyRank ==0) Then
           Write(Log_Unit, *) '[ERROR] Unknown Params/Init_U: ', Params%Init_U
        End If
        STOP
     End Select
     
     If (Params%Do_Irrev) Then
        Call Update_BC_V(Geom, Params, MySD_V, Node_db_V, V_Dist, V_Old, TimeStep-1)
!        Call Apply_BC_V(Geom, Params, MySD_V, Node_db_V, V_Dist)
!        Call Export(max(TimeStep-1, 1)) 
     End If
     Call Assemb_RHS_V(RHS_V, Geom, Params, MySD_V, Elem_db_V, Node_db_V)
     
     Select Case (Params%Init_V)
     Case (Init_V_ONE)
        If (.NOT. Is_BackTracking) Then
          Call VecSet(1.0_Kr, V_Dist, iErr)
          Call VecGhostUpdateBegin(V_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)
          Call VecGhostUpdateEnd(V_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)
	End If 	
     Case (Init_V_RND)
        If (.NOT. Is_BackTracking) Then
          Call Init_V_Cracks(Geom, Params, MySD_V, Elem_db_V,  Node_db_V, V_Dist)
!           Call Export (TimeStep+Size(Params%Load)+1)
	End If  
     Case (Init_V_PREV)        
	Continue
     Case Default
        If (MEF90_MyRank ==0) Then
           Write(Log_Unit, *) '[ERROR] Unknown Params/Init_V: ', Params%Init_V
        End If
        STOP
     End Select
     Call Export_V(TimeStep)

     Call VecCopy(V_Dist, V_Old, iErr)

     Do_Iter: Do iIter = 1, Params%MaxIterRelax
        If (MEF90_MyRank == 0) Then
           Write(Log_Unit, 410) TimeStep, iIter
        End If

        Call PetscGetTime(InitTS, iErr)
        
        Call Solve_U()
        Call KSPGetConvergedReason(KSP_U, KSP_TestCVG, iErr)
        If (KSP_TestCVG <= 0) Then
           If (MEF90_MyRank == 0) Then
              Write(Log_Unit, *) '[ERROR] KSPConvergedReason returned ', KSP_TestCVG
           End If
        End If
        Call PetscGetTime(InitTF, iErr)

        Call KSPGetIterationNumber(KSP_U, NbIterKSP, iErr)

        If (MEF90_MyRank == 0) Then
           Write(Log_Unit, 500) NbIterKSP, InitTF - InitTS, KSP_TestCVG
        End If

        Call PetscGetTime(InitTS, iErr)
        Call Solve_V()
        Call KSPGetConvergedReason(KSP_V, KSP_TestCVG, iErr)
        If (KSP_TestCVG <= 0) Then
           If (MEF90_MyRank == 0) Then
              Write(Log_Unit, *) '[ERROR] KSPConvergedReason returned ', KSP_TestCVG
           End If
        End If
        Call PetscGetTime(InitTF, iErr)

        Call KSPGetIterationNumber(KSP_V, NbIterKSP, iErr)
        If (MEF90_MyRank == 0) Then
           Write(Log_Unit, 600) NbIterKSP, InitTF - InitTS, KSP_TestCVG
        End If
        

        Call VecMax(V_Dist, VMaxPos, VMax, iErr)
        Call VecMin(V_Dist, VMinPos, VMin, iErr)
        If (MEF90_MyRank ==0) Then
           Write(Log_Unit, 700) VMin, VMax
        End If
        
        Call VecCopy(V_Old, V_Change, iErr)
        Call VecAxPy(-1.0_Kr, V_Dist, V_Change, iErr)
        Call VecNorm(V_Change, NORM_INFINITY, ErrV, iErr)
        If (MEF90_MyRank ==0) Then
           Write(Log_Unit, 800) ErrV
        End If
      
        If (Mod(iIter, SaveInt) == 0) Then 
           If (MEF90_MyRank ==0) Then
              Write(Log_Unit, *) 'Saving intermediate result for time step ', &
                   & TimeStep
              Call PetscGetTime(TotalTF, iErr)
              Write(Log_Unit, 930) TotalTF - TotalTS
           End If
           Call Export (TimeStep)
           Call Comp_Bulk_Ener(Bulk_Ener(TimeStep), U_Loc, V_Loc, Geom,       &
                & Params, MySD_U, MySD_V, Elem_db_U, Elem_db_V, Node_db_U,    &
                & Node_db_V, F_Loc, Temp_Loc )
           
           Call Comp_Surf_Ener(Surf_Ener(TimeSTep), V_Loc, Geom, Params,      &
                & MySD_V, Elem_db_V, Node_db_V)
           Tot_Ener(TimeStep) = Bulk_Ener(TimeStep) + Surf_Ener(TimeSTep)
           If (MEF90_MyRank ==0) Then
              Write(Log_Unit, 910) Bulk_Ener(TimeStep), Surf_Ener(TimeStep),  &
                   & Tot_Ener(TimeStep)
           End If
        End If


        Call VecCopy(V_Dist, V_Old, iErr)
        If ((ErrV .LE. Params%TolRelax) .AND. (iIter > 1) ) Then
           If (MEF90_MyRank ==0) Then
              Write(Log_Unit, *) 'Saving final result for time step ', &
                   & TimeStep
              Call PetscGetTime(TotalTF, iErr)
              Write(Log_Unit, 930) TotalTF - TotalTS
           End If
           Call Export(TimeStep) 
           Call Comp_Bulk_Ener(Bulk_Ener(TimeStep), U_Loc, V_Loc, Geom,       &
                & Params, MySD_U, MySD_V, Elem_db_U, Elem_db_V, Node_db_U,    &
                & Node_db_V, F_Loc, Temp_Loc )

           Call Comp_Surf_Ener(Surf_Ener(TimeStep), V_Loc, Geom, Params,      &
                & MySD_V, Elem_db_V, Node_db_V)
           Tot_Ener = Bulk_Ener + Surf_Ener

           If (MEF90_MyRank == 0) Then
              Write(Log_Unit, 910) Bulk_Ener(TimeStep), Surf_Ener(TimeStep),   &
                   & Tot_Ener(TimeStep)
              Open(File = Ener_Str, Unit =  Ener_Unit, Position = 'append')
              Write(Ener_Unit, 920) TimeStep, Params%Load(TimeStep),          &
                   & Bulk_Ener(TimeStep), Surf_Ener(TimeStep),                &
                   & Tot_Ener(TimeStep)
              Close(Ener_Unit)
           End If
           EXIT 
        End If
     End Do Do_Iter

     Is_BackTracking = .FALSE.
     If ( .NOT. Params%Do_BackTrack) Then
        TimeStep = TimeStep + 1 		 
     Else
        If (MEF90_MyRank == 0) Then  
           Write(Log_Unit, *) 'Current energy: ', Tot_Ener(TimeStep)
           Do iE = 1, TimeStep - 1
!              Write(Log_Unit, *) 'Comparing with iteration ', iE
              Tmp_Ener = Bulk_Ener(TimeStep) / Params%Load(TimeStep)**2 *     &
                   & Params%Load(iE)**2 + Surf_Ener(TimeStep)
!              Write(Log_Unit, *) 'Rescaled energy: ', Tmp_Ener
!              Write(Log_Unit, *) 'Previous energy: ', Tot_Ener(iE)
              !	      Err_Rel_Ener = (Tot_Ener(TimeStep) - Tot_Ener(iE))/Tot_Ener(TimeStep)	     
	      Err_Rel_Ener = (Tmp_Ener - Tot_Ener(iE)) / ABS(Tot_Ener(TimeStep))	     
!              Write(Log_Unit, *) 'Change Ratio:    ', Err_Rel_Ener
              If (Err_Rel_Ener < -Params%Tol_Ener) Then          
                 Is_BackTracking = .True.
                 EXIT
              End If
           End Do
           If (Is_BackTracking) Then
              TimeStep = iE
              Write(Log_Unit, *) '********** Going back to step ', TimeStep
              Open(File = Ener_Str, Unit =  Ener_Unit, Position = 'append')
              Write(Ener_Unit, *)
              Close(Ener_Unit)
           Else
              TimeStep = TimeStep + 1
           End If
        End If
        Call MPI_BCAST(Is_BackTracking , 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, iErr)
        Call MPI_BCAST(TimeStep, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, iErr)
     End If
  End Do While_TS

  Call PetscGetTime(TotalTF, iErr)
  Call PetscGetFlops(MyFlops, iErr)
  Call MPI_Reduce(MyFlops, TotalFlops, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0,   &
       & MPI_COMM_WORLD, iErr)
  If (MEF90_MyRank == 0) Then
     Write(Log_Unit, *) 'Total Time:                            ',      &
          & TotalTF - TotalTS
     Write(Log_Unit, *) 'Total FLOPs:                           ',      &
          & TotalFlops
     Write(Log_Unit, *) 'FLOP.s^-1:                             ',      &
          & TotalFlops / (TotalTF - TotalTS)
  End If

  Call Finalize()
  STOP

400 Format('**** Time step: ', I5, ' Load = ', ES12.5)
410 Format('==== Iteration: ', I5, ' /', I5)
500 Format('     Solve_U: ', I5, ' iterations in ', F7.2, 's. Return value',  &
         & I2)
600 Format('     Solve_V: ', I5, ' iterations in ', F7.2, 's. Return value',  &
         & I2)
700 Format('     VMin / Max:   ', T24, 2(ES12.5, '  '))
800 Format('     Max change V: ', T24, ES12.5)
901 Format('     Surface:      ', T24, ES12.5)
902 Format('     Bulk:         ', T24, ES12.5)
910 Format('     Energies:     ', T24, 2(ES10.3, ' '), 'Total: ', ES10.3)
920 Format(I4, 4(ES13.5,'  '))
930 Format('     Cumulated time: ', T24, ES12.5)
#if defined PB_2D
End Program  Rupt2D
#elif defined PB_3D
End Program Rupt3D
#else 
End Program Rupt2DA
#endif

