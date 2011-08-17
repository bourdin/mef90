Module m_MEF_Types
#include "finclude/petscdef.h"

   Use m_MEF_LinAlg
   Use petsc
   
   IMPLICIT NONE
   Private
   
   Public :: Field
   Public :: Flag
   Public :: Element1D
   Public :: Element2D, Element2D_Scal, Element2D_Elast 
   Public :: Element3D, Element3D_Scal, Element3D_Elast 
   Public :: BoundaryElement2D, BoundaryElement2D_Scal, BoundaryElement2D_Elast 
   Public :: BoundaryElement3D, BoundaryElement3D_Scal, BoundaryElement3D_Elast 

   Public :: Elem_Blk_Type, Side_Set_Type, Node_Set_Type, MeshTopology_Type
   Public :: EXO_Type, EXO_Property_Type, EXO_Variable_Type
   
   Public :: EXOView
   Public :: MeshTopologyDestroy, MeshTopologyView
      
   Type Field
      Type(SectionReal)                              :: Sec
      Type(SectionReal), Dimension(:), Pointer       :: Component_Sec
      Type(Vec)                                      :: Vec
      Type(Vec)                                      :: LocalVec
      Type(VecScatter)                               :: Scatter
      PetscInt                                       :: num_components
      PetscInt, Dimension(:), Pointer                :: component_size
   End Type Field
   
   Type Flag
      Type(SectionInt)                               :: Sec
      Type(SectionInt), Dimension(:), Pointer        :: Component_Sec
      PetscInt                                       :: num_components
      PetscInt, Dimension(:), Pointer                :: component_size
   End Type Flag


   Type Element1D
      PetscReal, Dimension(:,:), Pointer             :: BF
      PetscReal, Dimension(:,:), Pointer             :: Der_BF
      PetscReal, Dimension(:), Pointer               :: Gauss_C
   End Type Element1D
 
   Type Element2D_Scal
      PetscReal, Dimension(:,:), Pointer             :: BF
      Type(Vect2D), Dimension(:,:), Pointer          :: Grad_BF
      PetscReal, Dimension(:), Pointer               :: Gauss_C
   End Type Element2D_Scal
 
   Type Element2D
      Type (Vect2D), Dimension(:,:), Pointer         :: BF
      Type (Mat2D), Dimension(:,:), Pointer          :: Der_BF
      PetscReal, Dimension(:), Pointer               :: Gauss_C
   End Type Element2D
 
   Type Element2D_Elast
      Type (Vect2D), Dimension(:,:), Pointer         :: BF
      Type (MatS2D), Dimension(:,:), Pointer         :: GradS_BF
      PetscReal, Dimension(:), Pointer               :: Gauss_C
   End Type Element2D_Elast
 
   Type Element3D
      Type (Vect3D), Dimension(:,:), Pointer         :: BF
      Type (Mat3D), Dimension(:,:), Pointer          :: Der_BF
      PetscReal, Dimension(:), Pointer               :: Gauss_C
   End Type Element3D
 
   Type Element3D_Scal
      PetscReal, Dimension(:,:), Pointer             :: BF
      Type (Vect3D), Dimension(:,:), Pointer         :: Grad_BF
      PetscReal, Dimension(:), Pointer               :: Gauss_C
   End Type Element3D_Scal
 
   Type Element3D_Elast
      Type (Vect3D), Dimension(:,:), Pointer         :: BF
      Type (MatS3D), Dimension(:,:), Pointer         :: GradS_BF
      PetscReal, Dimension(:), Pointer               :: Gauss_C
   End Type Element3D_Elast
   
   Type BoundaryElement2D_Scal
      Type(Vect2D)                                   :: NormalVector
      PetscReal, Dimension(:,:), Pointer             :: BF
      Type(Vect2D), Dimension(:,:), Pointer          :: Grad_BF
      PetscReal, Dimension(:), Pointer               :: Gauss_C
   End Type BoundaryElement2D_Scal
    
   Type BoundaryElement2D
      Type(Vect2D)                                   :: NormalVector
      Type(Vect2D), Dimension(:,:), Pointer          :: BF
      Type(Mat2D), Dimension(:,:), Pointer           :: Der_BF
      PetscReal, Dimension(:), Pointer               :: Gauss_C
   End Type BoundaryElement2D
    
   Type BoundaryElement2D_Elast
      Type(Vect2D)                                   :: NormalVector
      Type(Vect2D), Dimension(:,:), Pointer          :: BF
      Type(MatS2D), Dimension(:,:), Pointer          :: GradS_BF
      PetscReal, Dimension(:), Pointer               :: Gauss_C
   End Type BoundaryElement2D_Elast
    
   Type BoundaryElement3D_Scal
      Type(Vect3D)                                   :: NormalVector
      PetscReal, Dimension(:,:), Pointer             :: BF
      Type(Vect3D), Dimension(:,:), Pointer          :: Grad_BF
      PetscReal, Dimension(:), Pointer               :: Gauss_C
   End Type BoundaryElement3D_Scal
    
   Type BoundaryElement3D
      Type(Vect3D)                                   :: NormalVector
      Type(Vect3D), Dimension(:,:), Pointer          :: BF
      Type(Mat3D), Dimension(:,:), Pointer           :: Der_BF
      PetscReal, Dimension(:), Pointer               :: Gauss_C
   End Type BoundaryElement3D
    
   Type BoundaryElement3D_Elast
      Type(Vect3D)                                   :: NormalVector
      Type(Vect3D), Dimension(:,:), Pointer          :: BF
      Type(MatS3D), Dimension(:,:), Pointer          :: GradS_BF
      PetscReal, Dimension(:), Pointer               :: Gauss_C
   End Type BoundaryElement3D_Elast
    
   Type Elem_Blk_Type
      PetscInt                                       :: ID
      PetscInt                                       :: Elem_Type
      PetscInt, Dimension(4)                         :: DoF_Location
      !!! DoF location is Cells, faces, edges, vertices in 3D and
      !!!                  Cell, unused, edges, vertices in2D
      PetscInt                                       :: Num_DoF !! = sum(DoF_Location)
      PetscInt                                       :: Num_Elems
      PetscInt, Dimension(:), Pointer                :: Elem_ID
      PetscInt                                       :: Num_Face      
      PetscInt                                       :: Num_Edge
      PetscInt                                       :: Num_Vert
   End Type Elem_Blk_Type
 
   Type Side_Set_Type
      PetscInt                                       :: ID
      PetscInt                                       :: Elem_Type
      PetscInt, Dimension(3)                         :: DoF_Location
      !!! DoF location is Faces, edges, vertices in 3D and
      !!!                 Edges, vertices in2D
      PetscInt                                       :: Num_DoF !! = sum(DoF_Location)
      PetscInt                                       :: Num_Elems
      PetscInt, Dimension(:), Pointer                :: Elem_ID
      PetscInt                                       :: Num_Edge
      PetscInt                                       :: Num_Vert
   End Type Side_Set_Type
 
   Type Node_Set_Type
      PetscInt                                       :: ID
      PetscInt                                       :: Num_Nodes
      PetscInt, Dimension(:), Pointer                :: Node_ID
   End Type Node_Set_Type
 
   Type MeshTopology_Type
!      Sequence
      ! Global datas
      PetscInt                                       :: num_dim
      PetscInt                                       :: num_verts
      PetscInt                                       :: num_elems
      ! Element Blocks datas
      PetscInt                                       :: num_elem_blks_global
      PetscInt                                       :: num_elem_blks
      Type(Elem_Blk_Type), Dimension(:), Pointer     :: elem_blk
      ! Node sets datas
      PetscInt                                       :: num_node_sets_global
      PetscInt                                       :: num_node_sets 
      Type(Node_Set_Type), Dimension(:), Pointer     :: node_set
      ! Side Sets DATAS
      PetscInt                                       :: num_side_sets_global
      PetscInt                                       :: num_side_sets
      Type(DM)                                       :: mesh
   End Type MeshTopology_Type
   
   Type EXO_Type
      MPI_Comm                                       :: comm      
      Integer                                        :: exoid
      Character(len=MXLNLN)                          :: filename  
      Character(len=MXLNLN)                          :: title     
      ! QA DATAS
      Integer                                        :: num_QA    
      Character(len=MXSTLN), Dimension(:,:), Pointer :: QA_rec
      ! Properties
      PetscInt                                       :: Num_EBProperties
      Type(EXO_Property_Type), Dimension(:), Pointer :: EBProperty
      PetscInt                                       :: Num_SSProperties
      Type(EXO_Property_Type), Dimension(:), Pointer :: SSProperty
      PetscInt                                       :: Num_NSProperties
      Type(EXO_Property_Type), Dimension(:), Pointer :: NSProperty
      ! Variables
      PetscInt                                       :: Num_GlobVariables    
      Type(EXO_Variable_Type), Dimension(:), Pointer :: GlobVariable
      PetscInt                                       :: Num_CellVariables    
      Type(EXO_Variable_Type), Dimension(:), Pointer :: CellVariable
      PetscInt                                       :: Num_VertVariables    
      Type(EXO_Variable_Type), Dimension(:), Pointer :: VertVariable
   End Type EXO_Type
   
   Type EXO_Property_Type
      Character(MXSTLN)                              :: Name
      PetscInt, Dimension(:), Pointer                :: Value
   End Type EXO_Property_Type
   
   Type EXO_Variable_Type
      Character(MXSTLN)                              :: Name
      PetscInt                                       :: Offset
      !!! the position of the variable in the exo file
   End Type EXO_Variable_Type   
   
Contains
   Subroutine EXOView(dEXO, viewer)
      Type(EXO_Type)                 :: dEXO
      Type(PetscViewer)              :: viewer
      
      PetscInt                       :: i, j, iErr
      Character(len=MEF90_MXSTRLEN)  :: CharBuffer
   
      If (dEXO%comm == PETSC_COMM_WORLD) Then
         Write(CharBuffer, 100) 'PETSC_COMM_WORLD'
      ElseIf (dEXO%comm == PETSC_COMM_SELF) Then
         Write(CharBuffer, 100) 'PETSC_COMM_SELF'
      Else  
         Write(CharBuffer, 105) 'unknown', dEXO%comm
      End If
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      
      Write(CharBuffer, 101) dEXO%exoid
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Write(CharBuffer, 102) dEXO%filename
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
!      Write(CharBuffer, 103) dEXO%title//' '
!      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
!      Write(CharBuffer, '(A)') '\n'

      Write(CharBuffer, 106) dEXO%Num_EBProperties
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Do i = 1, dEXO%Num_EBProperties
         Write(CharBuffer, 109) dEXO%EBProperty(i)%Name, Size(dEXO%EBProperty(i)%Value)
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do j = 1, Size(dEXO%EBProperty(i)%Value)
            Write(CharBuffer, 201) j, dEXO%EBProperty(i)%Value(j)
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
      End Do

      Write(CharBuffer, 107) dEXO%Num_SSProperties
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Do i = 1, dEXO%Num_SSProperties
         Write(CharBuffer, 109) dEXO%SSProperty(i)%Name, Size(dEXO%SSProperty(i)%Value)
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do j = 1, Size(dEXO%SSProperty(i)%Value)
            Write(CharBuffer, 202) j, dEXO%SSProperty(i)%Value(j)
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
      End Do

      Write(CharBuffer, 108) dEXO%Num_NSProperties
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Do i = 1, dEXO%Num_NSProperties
         Write(CharBuffer, 109) dEXO%NSProperty(i)%Name, Size(dEXO%NSProperty(i)%Value)
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do j = 1, Size(dEXO%NSProperty(i)%Value)
            Write(CharBuffer, 203) j, dEXO%NSProperty(i)%Value(j)
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
      End Do
!      Do i = 1, dEXO%Num_QA
!         Write(CharBuffer, 104) i, dEXO%QA_rec(i,:)
!         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
!      End Do
      Write(CharBuffer, 110) dEXO%Num_GlobVariables
      Do i = 1, dEXO%Num_GlobVariables
         Write(CharBuffer, 210) dEXO%GlobVariable(i)%Name, dEXO%GlobVariable(i)%offset
      End Do

      Write(CharBuffer, 111) dEXO%Num_CellVariables
      Do i = 1, dEXO%Num_CellVariables
         Write(CharBuffer, 210) dEXO%CellVariable(i)%Name, dEXO%CellVariable(i)%offset
      End Do

      Write(CharBuffer, 112) dEXO%Num_VertVariables
      Do i = 1, dEXO%Num_VertVariables
         Write(CharBuffer, 210) dEXO%VertVariable(i)%Name, dEXO%VertVariable(i)%offset
      End Do
    
      
 100 Format('Communicator:       ', A,       '\n')
 101 Format('exo ID:             ', I3,      '\n')
 102 Format('filename:           ', A,       '\n')
 105 Format('Communicator:       ', A, I3,   '\n')
 106 Format('Number of EB Properties: ', I3, '\n')
 107 Format('Number of SS Properties: ', I3, '\n')
 108 Format('Number of NS Properties: ', I3, '\n')
 109 Format('   ', A, 'num values:', I3,     '\n')
 110 Format('Global Variables: ', I3,        '\n')
 111 Format('Cell Variables:   ', I3,        '\n')
 112 Format('Vertex Variables: ', I3,        '\n')
 201 Format('   Element Block ', I3, ' value ', I3, '\n')
 202 Format('   Side Set      ', I3, ' value ', I3, '\n')
 203 Format('   Node Set      ', I3, ' value ', I3, '\n')
 210 Format(A, I3, '\n')
   End Subroutine EXOView


   Subroutine MeshTopologyDestroy(dMeshTopology)
      Type(MeshTopology_Type)         :: dMeshTopology
      PetscInt                        :: iSet, iBlk
      
      If (Size(dMeshTopology%Node_set) > 0) Then
         Do iSet = 1, Size(dMeshTopology%Node_set)
            If (dMeshTopology%Node_set(iSet)%num_nodes>0) Then
               Deallocate(dMeshTopology%Node_Set(iSet)%Node_ID)
            End If
         End Do
         Deallocate (dMeshTopology%Node_Set)
      End If
      If (Size(dMeshTopology%Elem_Blk) > 0) Then
         Do iBlk = 1, Size(dMeshTopology%Elem_Blk)
            If (dMeshTopology%Elem_Blk(iBlk)%Num_Elems > 0) Then
               Deallocate(dMeshTopology%Elem_blk(iBlk)%Elem_ID)
            End If
         End Do
         Deallocate(dMeshTopology%Elem_blk)
      End If
   End Subroutine MeshTopologyDestroy

   Subroutine MeshTopologyView(dMeshTopology, viewer)
      Type(MeshTopology_Type), Intent(IN)           :: dMeshTopology
      Type(PetscViewer)                             :: viewer
      
      Integer                                        :: iErr
      Character(len=MEF90_MXSTRLEN)                  :: CharBuffer
      Integer                                        :: i, j

      Write(CharBuffer, 103) dMeshTopology%num_dim
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Write(CharBuffer, 104) dMeshTopology%num_Verts
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Write(CharBuffer, 105) dMeshTopology%num_elems
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Write(CharBuffer, 106) dMeshTopology%Num_Elem_blks
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Write(CharBuffer, 107) dMeshTopology%Num_Node_Sets
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Write(CharBuffer, 108) dMeshTopology%Num_Side_Sets
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
   
      Write(CharBuffer, 600) '\n'
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Write(CharBuffer, 200)
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Write(CharBuffer, 601) dMeshTopology%num_elem_blks_global
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Write(CharBuffer, 201) dMeshTopology%num_elem_blks
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Do i = 1, dMeshTopology%Num_Elem_blks
         Write(CharBuffer, 203) dMeshTopology%Elem_Blk(i)%ID, dMeshTopology%Elem_blk(i)%Num_Elems
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Write(CharBuffer, 204) dMeshTopology%Elem_Blk(i)%ID, dMeshTopology%Elem_blk(i)%Elem_Type
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Write(CharBuffer, 205) dMeshTopology%Elem_Blk(i)%ID, dMeshTopology%Elem_blk(i)%DoF_Location
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Write(CharBuffer, 207) dMeshTopology%Elem_Blk(i)%ID
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do j = 1, dMeshTopology%Elem_blk(i)%Num_Elems
            Write(CharBuffer, 208) dMeshTopology%Elem_blk(i)%Elem_ID(j)
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, 600) '\n'
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      End Do
      
      
      Write(CharBuffer, 600) '\n'
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Write(CharBuffer, 300)
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Write(CharBuffer, 701) dMeshTopology%num_node_sets_global
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Write(CharBuffer, 301) dMeshTopology%num_node_sets
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Do i = 1, dMeshTopology%num_node_sets
         Write(CharBuffer, 302) dMeshTopology%Node_Set(i)%ID, dMeshTopology%Node_Set(i)%Num_Nodes
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Write(CharBuffer, 303) dMeshTopology%Node_Set(i)%ID
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         Do j = 1, dMeshTopology%Node_Set(i)%Num_Nodes
            Write(CharBuffer, 500) dMeshTopology%Node_Set(i)%Node_ID(j)
            Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
         End Do
         Write(CharBuffer, 600) '\n'
         Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      End Do
      
      
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Write(CharBuffer, 400)
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      Write(CharBuffer, 401) dMeshTopology%num_side_sets
      Call PetscViewerASCIIPrintf(viewer, CharBuffer, iErr); CHKERRQ(iErr)
      
103 Format('    Number of dimensions ============ ', I6, '\n')
104 Format('    Number of vertices ============== ', I6, '\n')
105 Format('    Number of elements ============== ', I6, '\n')
106 Format('    Number of elements blocks ======= ', I6, '\n')
107 Format('    Number of node sets ============= ', I6, '\n')
108 Format('    Number of side sets ============= ', I6, '\n')

200 Format('*** ELEMENT BLOCKS ***', '\n')
201 Format('    Number of blocks ================ ', I4, '\n')
203 Format('    Block ', I3, ' Number of elements ==== ', I4, '\n')
204 Format('    Block ', I3, ' Element type ========== ', I4, '\n')
205 Format('    Block ', I3, ' DoF location ========== ', 4(I4, ' '), '\n')
207 Format('    Block ', I3, ' IDs: ')
208 Format(' ', I4)

300 Format('*** NODE SETS ***', '\n')
301 Format('    Number of sets ================== ', I4, '\n')
302 Format('    Set ', I3, ' Number of nodes ========= ', I4, '\n')
303 Format('    Set ', I3, ' IDs: ')
!304 Format('    Set ', I3, ' Number of dist. factors = ', I4)

400 Format('*** SIDE SETS ***', '\n')
401 Format('    Number of side sets ============= ', I4, '\n')
    
500 Format(I4) 
600 Format(A)
601 Format('    Number of blocks global ========= ', I4, '\n')
701 Format('    Number of sets global =========== ', I4, '\n')
   End Subroutine MeshTopologyView
   
End Module m_MEF_Types
