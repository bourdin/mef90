Module m_PrepFilm
#include "finclude/petscdef.h"

Use m_MEF90
Use m_Film_Struct

Implicit NONE

Contains

Subroutine getmaterialprop(dMeshTopology, dMatProp2D, BatchUnit, IsBatch)
	Type(MeshTopology_Type)				:: dMeshTopology
	Type(MatProp2D_Type), Intent(OUT)		:: dMatProp2D
	PetscInt					:: BatchUnit
	PetscBool					:: IsBatch
	
	Character(len=MEF90_MXSTRLEN)			:: IOBuffer
	PetscReal					:: fractough, deltough, ksubst, E, nu
	PetscReal					:: thermalexpxx, thermalexpyy, thermalexpxy
	
	Write(IOBuffer, 300) 'Fracture toughness'			
	Call MEF90_AskReal(fractough, IOBuffer, BatchUnit, IsBatch)
	Write(IOBuffer, 300) 'Delamination toughness'
	Call MEF90_AskReal(deltough, IOBuffer, BatchUnit, IsBatch)
	Write(IOBuffer, 300) 'Substrate stiffness'
	Call MEF90_AskReal(ksubst, IOBuffer, BatchUnit, IsBatch)
	Write(IOBuffer, 300) 'Young Modulus'
	Call MEF90_AskReal(E, IOBuffer, BatchUnit, IsBatch)
	Write(IOBuffer, 300) 'Poisson Ratio'
	Call MEF90_AskReal(nu, IOBuffer, BatchUnit, IsBatch)
	Write(IOBuffer, 300) 'Therm Exp xx'
	Call MEF90_AskReal(thermalexpxx, IOBuffer, BatchUnit, IsBatch)	
	Write(IOBuffer, 300) 'Therm Exp yy'
	Call MEF90_AskReal(thermalexpyy, IOBuffer, BatchUnit, IsBatch)	
	Write(IOBuffer, 300) 'Therm Exp xy'
	Call MEF90_AskReal(thermalexpxy, IOBuffer, BatchUnit, IsBatch)	

	Select Case(dMeshTopology%Num_Dim)
	Case(2)
		dMatProp2D%FracToughness = fractough
		dMatProp2D%DelamToughness = deltough
		dMatProp2D%Ksubst = ksubst
		dMatProp2D%Therm_Exp    = 0.0_Kr
		dMatProp2D%Therm_Exp%XX = thermalexpxx
		dMatProp2D%Therm_Exp%YY = thermalexpyy
		dMatProp2D%Therm_Exp%XY = thermalexpxy
		Call GenHL_Iso2D_EnuPlaneStress(E, nu, dMatProp2D%Hookes_Law)
	End Select
	
	300 Format(A)
 
End Subroutine getmaterialprop
	
End Module m_PrepFilm
