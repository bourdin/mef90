Module m_PrepVarFrac
#include "finclude/petscdef.h"

   Use m_MEF90

   Implicit NONE

Contains
   Subroutine AskInt(val, msg, ArgUnit, IsBatch)
      PetscInt                                  :: Val
      Character(len=*)                          :: msg 
      PetscInt                                  :: argunit
      PetscBool                                 :: IsBatch

      Character(len=MEF90_MXSTRLEN)             :: prefix, IOBuffer   
      PetscInt                                  :: iErr   
      If (IsBatch) Then
         If (MEF90_MyRank == 0) Then
            Read(ArgUnit,*) Val
         End If
         Call MPI_BCast(Val, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, iErr)
      Else
         Write(IOBuffer, "(A, t60,':  ')") Trim(msg)
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         If (MEF90_MyRank == 0) Then
            Read(*,*) Val
            Write(ArgUnit, "(I4, t60, A)") val, Trim(msg)
         End If
         Call MPI_BCast(Val, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, iErr)
      End If
   End Subroutine AskInt   
   
   Subroutine AskReal(val, msg, ArgUnit, IsBatch)
      PetscReal                                 :: Val
      Character(len=*)                          :: msg 
      PetscInt                                  :: argunit
      PetscBool                                 :: IsBatch

      Character(len=MEF90_MXSTRLEN)             :: prefix, IOBuffer      
      PetscInt                                  :: iErr
      If (IsBatch) Then
         If (MEF90_MyRank == 0) Then
            Read(ArgUnit,*) Val
         End If
         Call MPI_BCast(Val, 1, MPIU_SCALAR, 0, PETSC_COMM_WORLD, iErr)
      Else
         Write(IOBuffer, "(A, t60,':  ')") Trim(msg)
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         If (MEF90_MyRank == 0) Then
            Read(*,*) Val
            Write(ArgUnit, "(ES12.5, t60, A)") val, Trim(msg)
         End If
         Call MPI_BCast(Val, 1, MPIU_SCALAR, 0, PETSC_COMM_WORLD, iErr)
      End If
   End Subroutine AskReal

	subroutine EXOEBProperty_AskWithBatchFilm(dEXO, dMeshTopology, BatchUnit, IsBatch)
	
		Type(EXO_Type)					:: dEXO
		Type(MeshTopology_Type)				:: dMeshTopology
		PetscInt					:: BatchUnit
		PetscBool					:: IsBatch
		PetscInt					:: iErr
		PetscInt					:: i, j, IntBuffer
		PetscInt					:: NumEB
		PetscInt					:: EXO_MyRank
		Character(len=MEF90_MXSTRLEN)			:: IOBuffer
		PetscReal					:: TmpEBProperty
			
		Do i = 1, dMeshTopology%Num_Elem_Blks_Global
			Write(IOBuffer, 100) i
			Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
			Do j = 1, dEXO%Num_EBProperties
				Write(IOBuffer, 200) i, Trim(dEXO%EBProperty(j)%Name)
				Call AskInt(dEXO%EBProperty(j)%Value(i), IOBuffer, BatchUnit, IsBatch)
			End Do
			If (.NOT. IsBatch) Then
				Write(BatchUnit, *)
			End If
		End Do
		
		100 Format('*** Element Block ', T24, I3, '\n')
		200 Format('EB', I4.4, ': ', A)
		
	end subroutine EXOEBProperty_AskWithBatchFilm

	subroutine EXONSProperty_AskWithBatchFilm(dEXO, dMeshTopology, BatchUnit, IsBatch)
	
		Type(EXO_Type)					:: dEXO
		Type(MeshTopology_Type)				:: dMeshTopology
		PetscInt					:: BatchUnit
		PetscBool					:: IsBatch
		PetscInt					:: iErr
		PetscInt					:: i, j, IntBuffer
		PetscInt					:: NumEB
		PetscInt					:: EXO_MyRank
		Character(len=MEF90_MXSTRLEN)			:: IOBuffer
		PetscReal					:: TmpEBProperty
		
		
		Do i = 1, dMeshTopology%Num_Node_Sets_Global
			Write(IOBuffer, 102) i
			Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
			Do j = 1, dEXO%Num_NSProperties
				Write(IOBuffer, 200) i, Trim(dEXO%NSProperty(j)%Name)
				Call AskInt(dEXO%NSProperty(j)%Value(i), IOBuffer, BatchUnit, IsBatch)
			End Do
			If (.NOT. IsBatch) Then
				Write(BatchUnit, *)
			End If
		End Do	
			
		102 Format('*** Node Set ', T24, I3, '\n')
		202 Format('NS', I4.4, ': ', A)
	
	end subroutine EXONSProperty_AskWithBatchFilm

	subroutine getmaterialprop(dMeshTopology, dMatProp2D, BatchUnit, IsBatch)
		Type(MeshTopology_Type)				:: dMeshTopology
		Type(MatProp2D_Type)				:: dMatProp2D
		PetscInt					:: BatchUnit
		PetscBool					:: IsBatch
		
		Character(len=MEF90_MXSTRLEN)			:: IOBuffer
		PetscReal					:: fractough, deltough, ksubst, E, nu
		PetscReal					:: thermalexpxx, thermalexpyy, thermalexpxy
		
		Write(IOBuffer, 300) 'Fracture toughness'			
		Call AskReal(fractough, IOBuffer, BatchUnit, IsBatch)
		Write(IOBuffer, 300) 'Delamination toughness'
		Call AskReal(deltough, IOBuffer, BatchUnit, IsBatch)
		Write(IOBuffer, 300) 'Substrate stiffness'
		Call AskReal(ksubst, IOBuffer, BatchUnit, IsBatch)
		Write(IOBuffer, 300) 'Young Modulus'
		Call AskReal(E, IOBuffer, BatchUnit, IsBatch)
		Write(IOBuffer, 300) 'Poisson Ratio'
		Call AskReal(nu, IOBuffer, BatchUnit, IsBatch)
		Write(IOBuffer, 300) 'Therm Exp xx'
		Call AskReal(thermalexpxx, IOBuffer, BatchUnit, IsBatch)	
		Write(IOBuffer, 300) 'Therm Exp yy'
		Call AskReal(thermalexpyy, IOBuffer, BatchUnit, IsBatch)	
		Write(IOBuffer, 300) 'Therm Exp xy'
		Call AskReal(thermalexpxy, IOBuffer, BatchUnit, IsBatch)	
	
		Select Case(MeshTopology%Num_Dim)
		Case(2)
			dMatProp2D%FracToughness = fractough
			dMatProp2D%DelamToughness = deltough
			dMatProp2D%Ksubst = ksubst
			dMatProp2D%Therm_Exp    = 0.0_Kr
			dMatProp2D%Therm_Exp%XX = thermalexpxx
			dMatProp2D%Therm_Exp%YY = thermalexpyy
			dMatProp2D%Therm_Exp%XY = thermalexpxy
			Call GenHL_Iso2D_EnuPlaneStress(E, nu, MatProp2D(iBlock)%Hookes_Law)
		End Select
		
		300 Format(A)
	 
	end subroutine getmaterialprop

   Subroutine EXOEBProperty_AskWithBatchGrains(dEXO, dMeshTopology, BatchUnit, IsBatch, NumGrains)
      Type(EXO_Type)                                 :: dEXO
      Type(MeshTopology_Type)                        :: dMeshTopology
      PetscInt                                       :: BatchUnit
      PetscBool                                      :: IsBatch
      PetscInt                                       :: NumGrains

      PetscInt                                       :: iErr
      PetscInt                                       :: i, j, IntBuffer

      PetscInt                                       :: NumEB
      PetscInt                                       :: EXO_MyRank
      Character(len=MEF90_MXSTRLEN)                  :: IOBuffer
      PetscReal                                      :: TmpEBProperty

      Write(IOBuffer, 110) NumGrains
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      Do j = 1, dEXO%Num_EBProperties
         Write(IOBuffer, 210) Trim(dEXO%EBProperty(j)%Name)
         Call AskReal(TmpEBProperty, IOBuffer, BatchUnit, IsBatch)
         Do i = 1, NumGrains
            dEXO%EBProperty(j)%Value(i) = TmpEBProperty
         End Do
      End Do
      If ((.NOT. IsBatch) .AND. (MEF90_MyRank == 0))Then
         Write(BatchUnit, *)
      End If
      Do i = NumGrains+1, dMeshTopology%Num_Elem_Blks_Global
         Write(IOBuffer, 100) i
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Do j = 1, dEXO%Num_EBProperties
            Write(IOBuffer, 200) i, Trim(dEXO%EBProperty(j)%Name)
            Call AskInt(dEXO%EBProperty(j)%Value(i), IOBuffer, BatchUnit, IsBatch)
         End Do
         If (.NOT. IsBatch) Then
            Write(BatchUnit, *)
         End If
      End Do
100 Format('*** Element Block ', T24, I3, '\n')
110 Format('*** Grains 1-', I4.4, '\n')
200 Format('EB', I4.4, ': ', A)
210 Format('Grains: ', A)
   End Subroutine EXOEBProperty_AskWithBatchGrains
   
   Subroutine EXONSProperty_AskWithBatch(dEXO, dMeshTopology, BatchUnit, IsBatch)
      Type(EXO_Type)                                 :: dEXO
      Type(MeshTopology_Type)                        :: dMeshTopology
      PetscInt                                       :: BatchUnit
      PetscBool                                      :: IsBatch

      PetscInt                                       :: iErr
      PetscInt                                       :: i, j, IntBuffer

      PetscInt                                       :: NumNS
      PetscInt                                       :: EXO_MyRank
      Character(len=MEF90_MXSTRLEN)                  :: IOBuffer

      Do i = 1, dMeshTopology%Num_Node_Sets_Global
         Write(IOBuffer, 102) i
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Do j = 1, dEXO%Num_NSProperties
            Write(IOBuffer, 202) i, Trim(dEXO%NSProperty(j)%Name)
            Call AskInt(dEXO%NSProperty(j)%Value(i), IOBuffer, BatchUnit, IsBatch)
         End Do
         If ((.NOT. IsBatch) .AND. (MEF90_MyRank == 0)) Then
            Write(BatchUnit, *)
         End If
      End Do
102 Format('*** Node Set      ', T24, I3, '\n')
202 Format('NS', I4.4, ': ', A)
   End Subroutine EXONSProperty_AskWithBatch

   Subroutine EXOEBProperty_AskWithBatch(dEXO, dMeshTopology, BatchUnit, IsBatch)
      Type(EXO_Type)                                 :: dEXO
      Type(MeshTopology_Type)                        :: dMeshTopology
      PetscInt                                       :: BatchUnit
      PetscBool                                      :: IsBatch

      PetscInt                                       :: iErr
      PetscInt                                       :: i, j, IntBuffer

      PetscInt                                       :: NumEB
      PetscInt                                       :: EXO_MyRank
      Character(len=MEF90_MXSTRLEN)                  :: IOBuffer
      PetscReal                                      :: TmpEBProperty

      Do i = 1, dMeshTopology%Num_Elem_Blks_Global
         Write(IOBuffer, 100) i
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Do j = 1, dEXO%Num_EBProperties
            Write(IOBuffer, 200) i, Trim(dEXO%EBProperty(j)%Name)
            Call AskInt(dEXO%EBProperty(j)%Value(i), IOBuffer, BatchUnit, IsBatch)
         End Do
         If ((.NOT. IsBatch) .AND. (MEF90_MyRank == 0)) Then
            Write(BatchUnit, *)
         End If
      End Do
100 Format('*** Element Block ', T24, I3, '\n')
200 Format('EB', I4.4, ': ', A)
   End Subroutine EXOEBProperty_AskWithBatch
   
   Subroutine EXONSProperty_AskWithBatchGrains(dEXO, dMeshTopology, BatchUnit, IsBatch, NumGrains)
      Type(EXO_Type)                                 :: dEXO
      Type(MeshTopology_Type)                        :: dMeshTopology
      PetscInt                                       :: BatchUnit
      PetscBool                                      :: IsBatch
      PetscInt                                       :: NumGrains

      PetscInt                                       :: iErr
      PetscInt                                       :: i, j, IntBuffer

      PetscInt                                       :: NumNS
      PetscInt                                       :: EXO_MyRank
      Character(len=MEF90_MXSTRLEN)                  :: IOBuffer

      Do i = 1, dMeshTopology%Num_Node_Sets_Global
         Write(IOBuffer, 102) i
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
         Do j = 1, dEXO%Num_NSProperties
            Write(IOBuffer, 202) i, Trim(dEXO%NSProperty(j)%Name)
            Call AskInt(dEXO%NSProperty(j)%Value(i), IOBuffer, BatchUnit, IsBatch)
         End Do
         If ((.NOT. IsBatch) .AND. (MEF90_MyRank == 0)) Then
            Write(BatchUnit, *)
         End If
      End Do
102 Format('*** Node Set      ', T24, I3, '\n')
202 Format('NS', I4.4, ': ', A)
   End Subroutine EXONSProperty_AskWithBatchGrains
End Module m_PrepVarFrac