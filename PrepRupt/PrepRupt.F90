Program PrepRupt

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscviewerdef.h"
#include "finclude/petscmeshdef.h"

   Use m_MEF90
   Use m_RuptStruct

   Use petsc
   Use petscvec
   Use petscmat
   Use petscksp
   Use petscmesh

   Implicit NONE   
   
   Type TestCase_Type
      PetscInt                                  :: Index
      Character(len=MEF90_MXSTRLEN)             :: Description
   End Type
   
   
   Character(len=MEF90_MXSTRLEN)                :: prefix, IOBuffer, filename
   Type(MeshTopology_Type)                      :: MeshTopology
   Type(EXO_Type)                               :: EXO, MyEXO
   Type(Mesh)                                   :: Tmp_Mesh
   PetscTruth                                   :: HasPrefix, verbose
   PetscInt                                     :: iErr, i, iCase

   PetscInt, Parameter                          :: NumTestCase=2
   Type(TestCase_Type), Dimension(NumTestCase)  :: TestCase

   Call MEF90_Initialize()
   
   Do i = 1, NumTestCase
      TestCase(i)%Index = i
   End Do
   TestCase(1)%Description = "MIL 2D/3D elasticity"
   TestCase(2)%Description = "MIL antiplane elasticity"
   
   
   Call PetscOptionsHasName(PETSC_NULL_CHARACTER, '-verbose', verbose, iErr)    
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-p',       prefix, HasPrefix, iErr); CHKERRQ(iErr)
   If (.NOT. HasPrefix) Then
      Call PetscPrintf(PETSC_COMM_WORLD, "No input file prefix given\n", iErr)
      Call MEF90_Finalize()
      STOP
   End If
   
   Call Write_EXO_Case(prefix, '%0.4d', MEF90_NumProcs)
   EXO%Comm = PETSC_COMM_WORLD
   EXO%filename = Trim(prefix)//'.gen'
   
   Call EXO_Check_Numbering(EXO, iErr)
   If (iErr /= 0) Then
      SETERRQ(PETSC_ERR_SUP, 'Unsupported numbering of the element blocks, side sets or node sets\n'c, iErr)
   End If

   !!! Reading and distributing sequential mesh
   Call MeshCreateExodus(PETSC_COMM_WORLD, EXO%filename, Tmp_mesh, ierr); CHKERRQ(iErr)
   Call MeshDistribute(Tmp_mesh, PETSC_NULL_CHARACTER, MeshTopology%mesh, ierr); CHKERRQ(iErr)
   Call MeshDestroy(Tmp_mesh, iErr); CHKERRQ(iErr)

   Call MeshTopologyReadEXO(MeshTopology, EXO)

   MyEXO%comm = PETSC_COMM_SELF
   MyEXO%exoid = EXO%exoid
   Write(MyEXO%filename, 200) trim(prefix), MEF90_MyRank
200 Format(A, '-', I4.4, '.gen')

   Call RuptEXOProperty_Init(MyEXO)   
   Call RuptEXOVariable_Init(MyEXO)   


   !!! Get Elem_Type first
!   Call Write_MeshTopologyGlobal(MeshTopology, MyEXO, PETSC_COMM_WORLD)
!   Call EXO_Variable_Write(MyEXO)
!   Call EXO_Property_Write(MyEXO)
!!!
!!!   Do i = 1, NumTestCase
!!!      Write(IOBuffer, "('[',I2.2,'] ',A)"), TestCase(i)%Index, TestCase(i)%Description//'\n'c
!!!      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
!!!   End Do
!!!   
!!!   Write(IOBuffer, "(A, t32)") 'Test Case'
!!!   Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
!!!   If (MEF90_MyRank == 0) Then
!!!      Read(*,*) ICase
!!!   End If
!!!   Call MPI_BCast(iCase, 1, MPI_INTEGER, 0, EXO%Comm, iErr)
!!!
!!!   Select Case(iCase)
!!!   
!!!   Case (1)
!!!
!!!   Case (2)
!!!
!!!   Case Default
!!!      SETERRQ(PETSC_ERR_SUP, 'Unknown test case\n'c, iErr)      
!!!   End Select
   
   Call MEF90_Finalize()
End Program PrepRupt