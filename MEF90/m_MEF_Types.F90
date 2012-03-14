Module m_MEF_Types
#include "finclude/petscdef.h"
   Use m_MEF_LinAlg
   Use petsc
   IMPLICIT NONE
   
   Private
   Public :: Field
   Public :: Flag
   Public :: Element1D
   Public :: Element2D,Element2D_Scal,Element2D_Elast 
   Public :: Element3D,Element3D_Scal,Element3D_Elast 

   Public :: CellSet_Type,MeshTopology_Type
   Public :: EXO_Type,EXO_Property_Type,EXO_Variable_Type
   
   Public :: EXOView
   Public :: MeshTopologyDestroy
      
   Type Field
      !!! A Field is a general container containing a section,a vector of sections for its components
      !!! a (global) vec,a local Vec,and the associated scatter.
      !!! All the Vec and SectionReal share the same memory space for their data storage
      Type(SectionReal)                            :: Sec
      Type(SectionReal),Dimension(:),Pointer       :: Component_Sec
      Type(Vec)                                    :: Vec
      Type(Vec)                                    :: LocalVec
      Type(VecScatter)                             :: Scatter
      PetscInt                                     :: num_components
      PetscInt,Dimension(:),Pointer                :: component_size
   End Type Field
   
   Type Flag
      !!! a Flag is to field what a SectionInt is to a SectionReal...
      !!! except that there is not VecInt,so the whole Vec and scatter part is ignored.
      Type(SectionInt)                             :: Sec
      Type(SectionInt),Dimension(:),Pointer        :: Component_Sec
      PetscInt                                     :: num_components
      PetscInt,Dimension(:),Pointer                :: component_size
   End Type Flag


   Type Element1D
      PetscReal,Dimension(:,:),Pointer             :: BF
      PetscReal,Dimension(:,:),Pointer             :: Der_BF
      PetscReal,Dimension(:),Pointer               :: Gauss_C
   End Type Element1D
 
   Type Element2D_Scal
      PetscReal,Dimension(:,:),Pointer             :: BF
      Type(Vect2D),Dimension(:,:),Pointer          :: Grad_BF
      PetscReal,Dimension(:),Pointer               :: Gauss_C
   End Type Element2D_Scal
 
   Type Element2D
      Type (Vect2D),Dimension(:,:),Pointer         :: BF
      Type (Mat2D),Dimension(:,:),Pointer          :: Der_BF
      PetscReal,Dimension(:),Pointer               :: Gauss_C
   End Type Element2D
 
   Type Element2D_Elast
      Type (Vect2D),Dimension(:,:),Pointer         :: BF
      Type (MatS2D),Dimension(:,:),Pointer         :: GradS_BF
      PetscReal,Dimension(:),Pointer               :: Gauss_C
   End Type Element2D_Elast
 
   Type Element3D
      Type (Vect3D),Dimension(:,:),Pointer         :: BF
      Type (Mat3D),Dimension(:,:),Pointer          :: Der_BF
      PetscReal,Dimension(:),Pointer               :: Gauss_C
   End Type Element3D
 
   Type Element3D_Scal
      PetscReal,Dimension(:,:),Pointer             :: BF
      Type (Vect3D),Dimension(:,:),Pointer         :: Grad_BF
      PetscReal,Dimension(:),Pointer               :: Gauss_C
   End Type Element3D_Scal
 
   Type Element3D_Elast
      Type (Vect3D),Dimension(:,:),Pointer         :: BF
      Type (MatS3D),Dimension(:,:),Pointer         :: GradS_BF
      PetscReal,Dimension(:),Pointer               :: Gauss_C
   End Type Element3D_Elast
   
   Type CellSet_Type
      PetscInt                                     :: ElemType
      PetscInt,Dimension(4)                        :: dofLocation
      !!! DoF location is Cells,faces,edges,vertices in 3D and
      !!!                  Cell,unused,edges,vertices in2D
      PetscInt                                     :: numDoF !! = sum(DoF_Location)
      PetscInt                                     :: coDimension
   End Type CellSet_Type
 
   Type MeshTopology_Type
      !!! One could debate if these should really be part of MeshTopology
      Type(IS)                                     :: cellSetGlobalIS,vertexSetGlobalIS
      Type(CellSet_Type),Dimension(:),Pointer      :: cellSet
      Type(DM)                                     :: mesh
   End Type MeshTopology_Type
   
!!! Hopefully, all of this will be moved into the PETSc exodus viewer   
   Type EXO_Type
      MPI_Comm                                     :: comm      
      Integer                                      :: exoid
      Character(len=MXLNLN)                        :: filename  
      Character(len=MXLNLN)                        :: title     
      ! QA DATAS
      Integer                                      :: num_QA    
      Character(len=MXSTLN),Dimension(:,:),Pointer :: QA_rec
      ! Properties
      PetscInt                                     :: Num_EBProperties
      Type(EXO_Property_Type),Dimension(:),Pointer :: EBProperty
      PetscInt                                     :: Num_NSProperties
      Type(EXO_Property_Type),Dimension(:),Pointer :: NSProperty
      ! Variables
      PetscInt                                     :: Num_GlobVariables    
      Type(EXO_Variable_Type),Dimension(:),Pointer :: GlobVariable
      PetscInt                                     :: Num_CellVariables    
      Type(EXO_Variable_Type),Dimension(:),Pointer :: CellVariable
      PetscInt                                     :: Num_VertVariables    
      Type(EXO_Variable_Type),Dimension(:),Pointer :: VertVariable
   End Type EXO_Type
   
   Type EXO_Property_Type
      !!! Derived type used to store the values of a property at ALL NS,EB
      !!! in an EXO file
      Character(MXSTLN)                            :: Name
      PetscInt,Dimension(:),Pointer                :: Value
   End Type EXO_Property_Type
   
   Type EXO_Variable_Type
      !!! Links variable name and offset in an exodus file
      Character(MXSTLN)                            :: Name
      PetscInt                                     :: Offset
      !!! the position of the variable in the exo file
   End Type EXO_Variable_Type   
   
Contains
#undef __FUNCT__
#define __FUNCT__ "EXOView"
   Subroutine EXOView(dEXO,viewer)
      Type(EXO_Type)                               :: dEXO
      Type(PetscViewer)                            :: viewer
      
      PetscInt                                     :: i,j,iErr
      Character(len=MEF90_MXSTRLEN)                :: CharBuffer
   
      If (dEXO%comm == PETSC_COMM_WORLD) Then
         Write(CharBuffer,100) 'PETSC_COMM_WORLD'
      ElseIf (dEXO%comm == PETSC_COMM_SELF) Then
         Write(CharBuffer,100) 'PETSC_COMM_SELF'
      Else  
         Write(CharBuffer,105) 'unknown',dEXO%comm
      End If
      Call PetscViewerASCIIPrintf(viewer,CharBuffer,iErr); CHKERRQ(iErr)
      
      Write(CharBuffer,101) dEXO%exoid
      Call PetscViewerASCIIPrintf(viewer,CharBuffer,iErr); CHKERRQ(iErr)
      Write(CharBuffer,102) dEXO%filename
      Call PetscViewerASCIIPrintf(viewer,CharBuffer,iErr); CHKERRQ(iErr)
!      Write(CharBuffer,103) dEXO%title//' '
!      Call PetscViewerASCIIPrintf(viewer,CharBuffer,iErr); CHKERRQ(iErr)
!      Write(CharBuffer,'(A)') '\n'

      Write(CharBuffer,106) dEXO%Num_EBProperties
      Call PetscViewerASCIIPrintf(viewer,CharBuffer,iErr); CHKERRQ(iErr)
      Do i = 1,dEXO%Num_EBProperties
         Write(CharBuffer,109) dEXO%EBProperty(i)%Name,Size(dEXO%EBProperty(i)%Value)
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,iErr); CHKERRQ(iErr)
         Do j = 1,Size(dEXO%EBProperty(i)%Value)
            Write(CharBuffer,201) j,dEXO%EBProperty(i)%Value(j)
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,iErr); CHKERRQ(iErr)
         End Do
      End Do

      Write(CharBuffer,108) dEXO%Num_NSProperties
      Call PetscViewerASCIIPrintf(viewer,CharBuffer,iErr); CHKERRQ(iErr)
      Do i = 1,dEXO%Num_NSProperties
         Write(CharBuffer,109) dEXO%NSProperty(i)%Name,Size(dEXO%NSProperty(i)%Value)
         Call PetscViewerASCIIPrintf(viewer,CharBuffer,iErr); CHKERRQ(iErr)
         Do j = 1,Size(dEXO%NSProperty(i)%Value)
            Write(CharBuffer,203) j,dEXO%NSProperty(i)%Value(j)
            Call PetscViewerASCIIPrintf(viewer,CharBuffer,iErr); CHKERRQ(iErr)
         End Do
      End Do
!      Do i = 1,dEXO%Num_QA
!         Write(CharBuffer,104) i,dEXO%QA_rec(i,:)
!         Call PetscViewerASCIIPrintf(viewer,CharBuffer,iErr); CHKERRQ(iErr)
!      End Do
      Write(CharBuffer,110) dEXO%Num_GlobVariables
      Do i = 1,dEXO%Num_GlobVariables
         Write(CharBuffer,210) dEXO%GlobVariable(i)%Name,dEXO%GlobVariable(i)%offset
      End Do

      Write(CharBuffer,111) dEXO%Num_CellVariables
      Do i = 1,dEXO%Num_CellVariables
         Write(CharBuffer,210) dEXO%CellVariable(i)%Name,dEXO%CellVariable(i)%offset
      End Do

      Write(CharBuffer,112) dEXO%Num_VertVariables
      Do i = 1,dEXO%Num_VertVariables
         Write(CharBuffer,210) dEXO%VertVariable(i)%Name,dEXO%VertVariable(i)%offset
      End Do
    
      
 100 Format('Communicator:       ',A,      '\n')
 101 Format('exo ID:             ',I3,     '\n')
 102 Format('filename:           ',A,      '\n')
 105 Format('Communicator:       ',A,I3,  '\n')
 106 Format('Number of EB Properties: ',I3,'\n')
 108 Format('Number of NS Properties: ',I3,'\n')
 109 Format('   ',A,'num values:',I3,    '\n')
 110 Format('Global Variables: ',I3,       '\n')
 111 Format('Cell Variables:   ',I3,       '\n')
 112 Format('Vertex Variables: ',I3,       '\n')
 201 Format('   Element Block ',I3,' value ',I3,'\n')
 203 Format('   Node Set      ',I3,' value ',I3,'\n')
 210 Format(A,I3,'\n')
   End Subroutine EXOView


#undef __FUNCT__
#define __FUNCT__ "MeshTopologyDestroy"
   Subroutine MeshTopologyDestroy(dMeshTopology)
      Type(MeshTopology_Type)                      :: dMeshTopology
      PetscErrorCode                               :: ierr
      
      Call ISDestroy(dMeshTopology%cellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call ISDestroy(dMeshTopology%vertexSetGlobalIS,ierr);CHKERRQ(ierr)
      DeAllocate(dMeshTopology%cellSet)
   End Subroutine MeshTopologyDestroy

End Module m_MEF_Types
